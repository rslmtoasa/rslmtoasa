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
