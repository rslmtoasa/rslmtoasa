!------------------------------------------------------------------------------
! RS-LMTO-ASA
!
! SUBMODULE: Self reciprocal helpers
!
!> Reciprocal-space SCF output and spinor moment helpers for `self`.
!------------------------------------------------------------------------------

submodule(self_mod) self_reciprocal

   implicit none

contains

   !=========================================================================
   !       FINALIZE THE ACCEPTED RECIPROCAL k-SPACE SCF STATE
   !=========================================================================
   !> The SCF DOS pass diagonalizes the Hamiltonian before mixing the newly
   !> computed moments into the accepted potential.  At convergence the
   !> difference is normally small, but the production handoff must not rely
   !> on that approximation.  Rebuild and diagonalize once, on the already
   !> accepted final potential, then solve EF from that same eigensystem.
   !> This routine deliberately does not map a new density back into SCF or
   !> perform another mixer update: it only refreshes the state being handed
   !> to TDDFT.
   module subroutine finalize_kspace_scf_state(this)
      class(self), intent(inout) :: this
      integer :: ia

      if (.not. this%use_kspace) return
      if (.not. this%converged) then
         call g_logger%fatal('self%finalize_kspace_scf_state: a non-converged k-space state cannot be handed to TDDFT.', &
                             __FILE__, __LINE__)
      end if
      if (.not. allocated(this%reciprocal_scf_cache)) then
         call g_logger%fatal('self%finalize_kspace_scf_state: reciprocal SCF cache is missing.', __FILE__, __LINE__)
      end if

      select case (this%control%calctype)
      case ('B')
         do ia = 1, this%lattice%nrec
            call this%symbolic_atom(ia)%build_pot()
         end do
         if (this%control%has_soc()) call this%hamiltonian%build_lsham()
         call this%hamiltonian%build_bulkham()
      case ('S', 'I', 'L')
         do ia = 1, this%lattice%ntype
            call this%symbolic_atom(ia)%build_pot()
         end do
         if (this%control%has_soc()) call this%hamiltonian%build_lsham()
         call this%hamiltonian%build_bulkham()
         if (this%control%calctype == 'I') call this%hamiltonian%build_locham()
      case default
         call g_logger%fatal('self%finalize_kspace_scf_state: unsupported calculation type.', __FILE__, __LINE__)
      end select

      call this%reciprocal_scf_cache%build_kspace_hamiltonian()
      call this%reciprocal_scf_cache%diagonalize_hamiltonian()
      this%reciprocal_scf_cache%auto_find_fermi = .true.
      this%bands%eband = this%reciprocal_scf_cache%calculate_canonical_band_energy(.true.)
      this%en%fermi = this%reciprocal_scf_cache%fermi_level
      call this%reciprocal_scf_cache%require_replicated_k_workset('self%finalize_kspace_scf_state')
   end subroutine finalize_kspace_scf_state

   !=========================================================================
   !       SERIALIZE THE ACCEPTED k-SPACE SCF STATE FOR HANDOFF AUDIT
   !=========================================================================
   !> The state file intentionally records the complete mesh, eigenvalues,
   !> explicit occupations, and occupation-weighted one-particle density
   !> matrix.  The latter is a gauge-invariant eigenvector/subspace identity
   !> diagnostic, so phase choices or rotations inside an occupied degenerate
   !> block cannot create a false mismatch.
   module subroutine write_kspace_scf_state_artifact(this, filename)
      class(self), intent(in) :: this
      character(len=*), intent(in) :: filename
      integer :: unit, ik, ib, i, j, nbasis, nbands, nk
      real(rp) :: kT, occupation, argument, weight_sum, k_fingerprint(5), moment
      real(rp), allocatable :: occupations(:)
      complex(rp), allocatable :: density_matrix(:, :)
      character(len=32) :: state_source

      if (.not. allocated(this%reciprocal_scf_cache)) then
         call g_logger%fatal('self%write_kspace_scf_state_artifact: reciprocal SCF cache is missing.', __FILE__, __LINE__)
      end if
      call this%reciprocal_scf_cache%require_replicated_k_workset('self%write_kspace_scf_state_artifact')
      if (.not. allocated(this%reciprocal_scf_cache%eigenvalues) .or. &
          .not. allocated(this%reciprocal_scf_cache%eigenvectors) .or. &
          .not. allocated(this%reciprocal_scf_cache%k_workset%points) .or. &
          .not. allocated(this%reciprocal_scf_cache%k_workset%weights)) then
         call g_logger%fatal('self%write_kspace_scf_state_artifact: accepted reciprocal eigensystem is incomplete.', &
                             __FILE__, __LINE__)
      end if
      if (rank /= 0) return

      nbands = size(this%reciprocal_scf_cache%eigenvalues, 1)
      nk = size(this%reciprocal_scf_cache%eigenvalues, 2)
      nbasis = size(this%reciprocal_scf_cache%eigenvectors, 1)
      if (size(this%reciprocal_scf_cache%k_workset%points, 2) /= nk .or. &
          size(this%reciprocal_scf_cache%k_workset%weights) /= nk .or. &
          size(this%reciprocal_scf_cache%eigenvectors, 2) /= nbands .or. &
          size(this%reciprocal_scf_cache%eigenvectors, 3) /= nk) then
         call g_logger%fatal('self%write_kspace_scf_state_artifact: accepted state array shapes disagree.', __FILE__, __LINE__)
      end if

      weight_sum = sum(this%reciprocal_scf_cache%k_workset%weights)
      k_fingerprint = 0.0_rp
      do ik = 1, nk
         k_fingerprint(1) = k_fingerprint(1) + this%reciprocal_scf_cache%k_workset%weights(ik)
         k_fingerprint(2:4) = k_fingerprint(2:4) + this%reciprocal_scf_cache%k_workset%weights(ik) * &
                               this%reciprocal_scf_cache%k_workset%points(:, ik)
         k_fingerprint(5) = k_fingerprint(5) + real(ik, rp) * this%reciprocal_scf_cache%k_workset%weights(ik)
      end do
      moment = 0.0_rp
      do ik = 1, this%lattice%nrec
         moment = moment + this%symbolic_atom(this%lattice%nbulk + ik)%potential%mtot
      end do
      kT = max(this%reciprocal_scf_cache%temperature*6.3336814e-6_rp, 1.0e-10_rp)
      allocate(occupations(nbands), density_matrix(nbasis, nbasis))

      state_source = 'accepted_kspace_scf'
      open(newunit=unit, file=trim(filename), status='replace', action='write')
      write(unit, '(a)') '# reciprocal state-consistency artifact'
      write(unit, '(a,a)') '# state_role = ', trim(state_source)
      write(unit, '(a,l1)') '# scf_converged = ', this%converged
      write(unit, '(a,i0)') '# scf_iterations = ', this%converged_iteration
      write(unit, '(a,3(i0,1x))') '# requested_k_mesh = ', this%reciprocal_scf_cache%nk_mesh
      write(unit, '(a,3(i0,1x))') '# actual_k_mesh = ', this%reciprocal_scf_cache%nk_mesh
      write(unit, '(a,i0)') '# actual_k_count = ', nk
      write(unit, '(a,i0)') '# nbasis = ', nbasis
      write(unit, '(a,i0)') '# nbands = ', nbands
      write(unit, '(a,es24.16)') '# k_weight_sum = ', weight_sum
      write(unit, '(a,5(es24.16,1x))') '# k_fingerprint_checksums = ', k_fingerprint
      write(unit, '(a,es24.16)') '# fermi_level_Ry = ', this%reciprocal_scf_cache%fermi_level
      write(unit, '(a,es24.16)') '# target_electron_count = ', this%reciprocal_scf_cache%total_electrons
      write(unit, '(a,es24.16)') '# accepted_electron_count = ', this%reciprocal_scf_cache%canonical_electron_count
      write(unit, '(a,es24.16)') '# accepted_electron_count_error = ', &
         this%reciprocal_scf_cache%canonical_electron_count - this%reciprocal_scf_cache%total_electrons
      write(unit, '(a,es24.16)') '# scf_moment_muB = ', moment
      write(unit, '(a,es24.16)') '# scf_residual = ', this%mix%delta
      write(unit, '(a,es24.16)') '# scf_physical_total_energy_Ry = ', this%physical_total_energy
      write(unit, '(a,es24.16)') '# accepted_potential_checksum = ', this%potential_checksum()
      write(unit, '(a,es24.16)') '# temperature_K = ', this%reciprocal_scf_cache%temperature
      write(unit, '(a,a)') '# reciprocal_mode = ', trim(this%reciprocal_scf_cache%reciprocal_mode)
      write(unit, '(a,a)') '# kspace_hamiltonian_order = ', trim(this%reciprocal_scf_cache%kspace_ham_order)
      write(unit, '(a)') '# occupation_semantics = reciprocal Fermi-Dirac electron-count solver; kT=max(T*kB,1e-10 Ry)'
      write(unit, '(a)') '# eigenvector_identity = occupation-weighted one-particle density-matrix projector residual'
      write(unit, '(a)') '# columns: tag k_index band_or_row column eigenvalue_or_occupation real imag'

      do ik = 1, nk
         write(unit, '(a,1x,i0,1x,4(es24.16,1x))') 'K', ik, this%reciprocal_scf_cache%k_workset%points(:, ik), &
            this%reciprocal_scf_cache%k_workset%weights(ik)
         do ib = 1, nbands
            argument = (this%reciprocal_scf_cache%eigenvalues(ib, ik) - this%reciprocal_scf_cache%fermi_level)/kT
            if (argument >= 50.0_rp) then
               occupation = 0.0_rp
            else if (argument <= -50.0_rp) then
               occupation = 1.0_rp
            else
               occupation = 1.0_rp/(exp(argument) + 1.0_rp)
            end if
            occupations(ib) = occupation
            write(unit, '(a,1x,2(i0,1x),2(es24.16,1x))') 'O', ik, ib, &
               this%reciprocal_scf_cache%eigenvalues(ib, ik), occupation
         end do
         density_matrix = cmplx(0.0_rp, 0.0_rp, rp)
         do i = 1, nbasis
            do j = 1, nbasis
               density_matrix(i, j) = sum(occupations(:) * this%reciprocal_scf_cache%eigenvectors(i, :, ik) * &
                  conjg(this%reciprocal_scf_cache%eigenvectors(j, :, ik)))
               write(unit, '(a,1x,3(i0,1x),2(es24.16,1x))') 'D', ik, i, j, real(density_matrix(i, j), rp), &
                  aimag(density_matrix(i, j))
            end do
         end do
      end do
      close(unit)
      deallocate(occupations, density_matrix)
   end subroutine write_kspace_scf_state_artifact

   !=========================================================================
   !  LR-02N SCALAR-RELATIVISTIC -> PAULI GROUND-STATE PROJECTION
   !=========================================================================
   module subroutine quantify_pauli_projection(this)
      use pauli_ground_state_projection_mod, only: pauli_ground_state_projection
      class(self), intent(inout) :: this

      type(reciprocal) :: reciprocal_obj
      integer :: ia
      logical :: use_shifted_kmesh

      if (this%control%nsp /= 1 .or. this%control%has_soc()) then
         call g_logger%fatal('[self.quantify_pauli_projection]: LR-02N requires control%nsp=1 (scalar-relativistic collinear mode).', &
                             __FILE__, __LINE__)
      end if
      if (this%hamiltonian%ccor_2c .or. this%hamiltonian%hubbard_u_general_check .or. &
          this%hamiltonian%hubbard_u_impurity_check .or. this%hamiltonian%hubbard_u_sc_check .or. &
          this%hamiltonian%hubbard_v_check) then
         call g_logger%fatal('[self.quantify_pauli_projection]: LR-02N requires a Hamiltonian without additive CCOR or Hubbard terms.', &
                             __FILE__, __LINE__)
      end if

      ! The k-space SCF branch already owns the accepted-potential spectrum.
      ! Reuse it so the diagnostic cannot silently change its mesh, weights, or
      ! eigensolver.  The ordinary real-space SCF route constructs the same
      ! reciprocal object from the final accepted potential below.
      if (allocated(this%reciprocal_scf_cache)) then
         if (allocated(this%reciprocal_scf_cache%eigenvalues) .and. &
             allocated(this%reciprocal_scf_cache%eigenvectors)) then
            this%reciprocal_scf_cache%fermi_level = this%en%fermi
            call pauli_ground_state_projection(this%reciprocal_scf_cache, this%symbolic_atom, &
                                               this%lattice%nbulk, 'lr02n_pauli_projection')
            return
         end if
      end if

      do ia = 1, this%lattice%nrec
         call this%symbolic_atom(ia)%build_pot()
      end do
      call this%hamiltonian%build_bulkham()

      reciprocal_obj = reciprocal(this%hamiltonian)
      reciprocal_obj%fermi_level = this%en%fermi
      reciprocal_obj%auto_find_fermi = .false.
      use_shifted_kmesh = sum(abs(reciprocal_obj%k_offset)) > 1.0e-12_rp
      if (reciprocal_obj%use_symmetry_reduction) then
         call reciprocal_obj%generate_reduced_kpoint_mesh(reciprocal_obj%nk_mesh, use_shifted_kmesh)
      else
         call reciprocal_obj%generate_mp_mesh()
      end if
      call reciprocal_obj%build_kspace_hamiltonian()
      call reciprocal_obj%diagonalize_hamiltonian()
      call pauli_ground_state_projection(reciprocal_obj, this%symbolic_atom, this%lattice%nbulk, &
                                         'lr02n_pauli_projection')
   end subroutine quantify_pauli_projection

   !=========================================================================
   !      WRITE k-SPACE SCF DOS OUTPUTS IN LEGACY-COMPATIBLE FILE FORMAT
   !=========================================================================
   module subroutine write_kspace_scf_dos_outputs(this, reciprocal_obj)
      class(self), intent(inout) :: this
      type(reciprocal), intent(in) :: reciprocal_obj
      integer :: i, isite, iorb
      real(rp), allocatable :: dos_up_tot(:), dos_dw_tot(:), dos_mx_tot(:), dos_my_tot(:), dos_mz_tot(:), dos_nmag_tot(:)

      if (rank /= 0) return
      if (.not. allocated(reciprocal_obj%total_dos)) return
      if (.not. allocated(reciprocal_obj%dos_energy_grid)) return

      allocate(dos_up_tot(reciprocal_obj%n_energy_points))
      allocate(dos_dw_tot(reciprocal_obj%n_energy_points))
      allocate(dos_mx_tot(reciprocal_obj%n_energy_points))
      allocate(dos_my_tot(reciprocal_obj%n_energy_points))
      allocate(dos_mz_tot(reciprocal_obj%n_energy_points))
      allocate(dos_nmag_tot(reciprocal_obj%n_energy_points))
      dos_up_tot = 0.0_rp
      dos_dw_tot = 0.0_rp
      dos_mx_tot = 0.0_rp
      dos_my_tot = 0.0_rp
      dos_mz_tot = 0.0_rp
      dos_nmag_tot = 0.0_rp

      if (allocated(reciprocal_obj%projected_dos)) then
         do i = 1, reciprocal_obj%n_energy_points
            do isite = 1, reciprocal_obj%n_sites
               do iorb = 1, reciprocal_obj%n_orb_types
                  dos_up_tot(i) = dos_up_tot(i) + reciprocal_obj%projected_dos(isite, iorb, 1, i)
                  dos_dw_tot(i) = dos_dw_tot(i) + reciprocal_obj%projected_dos(isite, iorb, 2, i)
               end do
            end do
         end do
      else
         dos_up_tot = 0.5_rp*reciprocal_obj%total_dos
         dos_dw_tot = 0.5_rp*reciprocal_obj%total_dos
      end if

      dos_mz_tot = dos_up_tot - dos_dw_tot
      dos_nmag_tot = 0.5_rp*(dos_up_tot + dos_dw_tot)

      if (allocated(reciprocal_obj%dos_mx_tot)) dos_mx_tot = reciprocal_obj%dos_mx_tot
      if (allocated(reciprocal_obj%dos_my_tot)) dos_my_tot = reciprocal_obj%dos_my_tot
      if (allocated(reciprocal_obj%dos_mz_tot)) dos_mz_tot = reciprocal_obj%dos_mz_tot

      open(unit=125, file='totaldos.out', status='replace', action='write')
      do i = 1, reciprocal_obj%n_energy_points
         write(125, '(2f16.5)') reciprocal_obj%dos_energy_grid(i) - reciprocal_obj%fermi_level, reciprocal_obj%total_dos(i)
      end do
      close(125)

      open(unit=126, file='magneticdos.out', status='replace', action='write')
      write(126, '(a)') '# energy_minus_fermi dos_up dos_dw dos_mx dos_my dos_mz dos_nmag'
      do i = 1, reciprocal_obj%n_energy_points
         write(126, '(7f16.5)') reciprocal_obj%dos_energy_grid(i) - reciprocal_obj%fermi_level, &
            dos_up_tot(i), dos_dw_tot(i), dos_mx_tot(i), dos_my_tot(i), dos_mz_tot(i), dos_nmag_tot(i)
      end do
      close(126)

      deallocate(dos_up_tot, dos_dw_tot, dos_mx_tot, dos_my_tot, dos_mz_tot, dos_nmag_tot)
   end subroutine write_kspace_scf_dos_outputs

   !=========================================================================
   !  SPINOR-RIGOROUS SITE MOMENTS FROM k-SPACE EIGENVECTORS (nsp=2/SOC path)
   !=========================================================================
   !> The returned components are in the eigenvector spin basis.  For
   !> `gbt_single_q` this is the primitive rotating frame: no cell-dependent
   !> spiral phase is present in a primitive k-space eigenvector sum.  The
   !> lab-frame moment belongs to output/comparison code and must be obtained
   !> through `gbt_structure_mod`’s authoritative frame transformation.
   module subroutine compute_kspace_spin_moments_spinor(this, reciprocal_obj, site_mom)
      class(self), intent(in) :: this
      type(reciprocal), intent(in) :: reciprocal_obj
      real(rp), intent(out) :: site_mom(3, this%lattice%nrec)
      integer :: ik, ik_global, ib, ia, io, iup, idn, site_offset
      real(rp) :: wk, occ, e, kT, farg
      complex(rp) :: u, d, ud
      real(rp) :: mxs, mys, mzs
      real(rp), parameter :: kB_Ry_per_K = 6.3336814e-6_rp

      site_mom(:, :) = 0.0_rp
      kT = reciprocal_obj%temperature*kB_Ry_per_K

      do ik = 1, size(reciprocal_obj%eigenvalues, 2)
         ik_global = ik
         if (allocated(reciprocal_obj%k_l2g_map) .and. ik <= size(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
         wk = reciprocal_obj%k_weights(ik_global)
         do ib = 1, size(reciprocal_obj%eigenvalues, 1)
            e = reciprocal_obj%eigenvalues(ib, ik)
            if (kT > 1.0e-12_rp) then
               farg = (e - reciprocal_obj%fermi_level)/kT
               if (farg > 50.0_rp) then
                  occ = 0.0_rp
               else if (farg < -50.0_rp) then
                  occ = 1.0_rp
               else
                  occ = 1.0_rp/(exp(farg) + 1.0_rp)
               end if
            else
               occ = merge(1.0_rp, 0.0_rp, e <= reciprocal_obj%fermi_level)
            end if
            if (occ <= 1.0e-14_rp) cycle

            do ia = 1, this%lattice%nrec
               site_offset = (ia - 1)*nb
               mxs = 0.0_rp
               mys = 0.0_rp
               mzs = 0.0_rp
               do io = 1, norb
                  iup = site_offset + io
                  idn = site_offset + norb + io
                  if (idn > size(reciprocal_obj%eigenvectors, 1)) cycle
                  u = reciprocal_obj%eigenvectors(iup, ib, ik)
                  d = reciprocal_obj%eigenvectors(idn, ib, ik)
                  ud = conjg(u)*d
                  mxs = mxs + 2.0_rp*real(ud, rp)
                  mys = mys + 2.0_rp*aimag(ud)
                  mzs = mzs + real(conjg(u)*u - conjg(d)*d, rp)
               end do
               site_mom(1, ia) = site_mom(1, ia) + wk*occ*mxs
               site_mom(2, ia) = site_mom(2, ia) + wk*occ*mys
               site_mom(3, ia) = site_mom(3, ia) + wk*occ*mzs
            end do
         end do
      end do
#ifdef USE_MPI
      if (reciprocal_obj%k_mesh_distributed_active) then
         call MPI_ALLREDUCE(MPI_IN_PLACE, site_mom, product(shape(site_mom)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      end if
#endif
   end subroutine compute_kspace_spin_moments_spinor

end submodule self_reciprocal
