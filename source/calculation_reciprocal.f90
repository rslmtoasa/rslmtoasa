submodule(calculation_mod) calculation_reciprocal

contains

   !> @brief Build the shared object stack for a post-processing driver.
   module subroutine prepare_post_processing_stack(this, use_paoflow, use_exchange_pairs, energy_mesh_before_hamiltonian, &
                                                   stochastic_moments, control_obj, lattice_obj, charge_obj, mix_obj, &
                                                   energy_obj, hamiltonian_obj, recursion_obj, dos_obj, green_obj, bands_obj, &
                                                   preprocessing_route)
      class(calculation), intent(in) :: this
      logical, intent(in) :: use_paoflow
      logical, intent(in) :: use_exchange_pairs
      logical, intent(in) :: energy_mesh_before_hamiltonian
      logical, intent(in) :: stochastic_moments
      type(control), target, intent(inout) :: control_obj
      type(lattice), target, intent(inout) :: lattice_obj
      type(charge), target, intent(inout) :: charge_obj
      type(mix), target, intent(inout) :: mix_obj
      type(energy), target, intent(inout) :: energy_obj
      type(hamiltonian), target, intent(inout) :: hamiltonian_obj
      type(recursion), target, intent(inout) :: recursion_obj
      type(dos), target, intent(inout) :: dos_obj
      type(green), target, intent(inout) :: green_obj
      type(bands), target, intent(inout) :: bands_obj
      character(len=*), intent(in), optional :: preprocessing_route
      character(len=sl) :: route
      integer :: i

      control_obj = control(this%fname)
      lattice_obj = lattice(control_obj)
      route = trim(control_obj%calctype)
      if (present(preprocessing_route)) then
         select case (trim(preprocessing_route))
         case ('bravais')
            route = 'B'
         case ('buildsurf')
            route = 'S'
         case ('newclubulk')
            route = 'I_bulk'
         case ('newclusurf')
            route = 'I_surface'
         case ('buildinterface')
            route = 'L'
         end select
      end if

      call g_timer%start('pre-processing')
      if (use_paoflow) then
         call lattice_obj%build_data()
         call lattice_obj%bravais()
         call lattice_obj%structb(.false.)
      else
         select case (trim(route))
         case ('B')
            call lattice_obj%build_data()
            call lattice_obj%bravais()
            call lattice_obj%structb(.true.)
         case ('S')
            call lattice_obj%build_data()
            call lattice_obj%bravais()
            call lattice_obj%build_surf_full()
            call lattice_obj%structb(.true.)
         case ('I', 'I_surface')
            call lattice_obj%build_data()
            call lattice_obj%bravais()
            call lattice_obj%build_surf_full()
            call lattice_obj%newclu()
            call lattice_obj%structb(.true.)
         case ('I_bulk')
            call lattice_obj%build_data()
            call lattice_obj%bravais()
            call lattice_obj%newclu()
            call lattice_obj%structb(.true.)
         case ('L')
            ! B7.5: build_interface_full instead of build_surf_full -- the
            ! two-sided counterpart. Mirrors pre_processing_buildinterface.
            call lattice_obj%build_data()
            call lattice_obj%bravais()
            call lattice_obj%build_interface_full()
            call lattice_obj%structb(.true.)
         end select
      end if

      call lattice_obj%atomlist()
      if (present(preprocessing_route)) then
         ! Spin dynamics consumes the active recursion atoms, matching the
         ! normal preprocessing ownership and the original SD partition.
         call get_mpi_variables(rank, lattice_obj%nrec)
      else if (use_exchange_pairs) then
         call get_mpi_variables(rank, lattice_obj%njij)
      else
         call get_mpi_variables(rank, lattice_obj%ntype)
      end if

      charge_obj = charge(lattice_obj)
      if (use_paoflow) then
         call charge_obj%bulkmat()
      else
         select case (control_obj%calctype)
         case ('B')
            call charge_obj%bulkmat()
         case ('S')
            call charge_obj%build_alelay
            call charge_obj%surfmat
         case ('I')
            call charge_obj%impmad()
            if (present(preprocessing_route)) call charge_obj%get_charge_transf
         case ('L')
            ! B7.5: same Madelung matrices as ’S’ (no interfacemat exists, by
            ! design), plus the region reference charges and the genuinely
            ! two-sided registry that overwrites surfmat’s one-sided one.
            ! Mirrors pre_processing_buildinterface.
            call charge_obj%get_charge_transf
            call charge_obj%build_alelay
            call charge_obj%surfmat
            call charge_obj%build_interface_registry()
         end select
      end if
      call g_timer%stop('pre-processing')

      mix_obj = mix(lattice_obj, charge_obj)
      energy_obj = energy(lattice_obj)
      if (energy_mesh_before_hamiltonian) call energy_obj%e_mesh()

      hamiltonian_obj = hamiltonian(charge_obj)
      if (use_paoflow) then
         select case (control_obj%calctype)
         case ('B')
            call hamiltonian_obj%build_from_paoflow_opt()
         case ('S')
            call g_logger%fatal('Surface calculation not implemented!', __FILE__, __LINE__)
         case ('I')
            call g_logger%fatal('Imputiry calculation not implemented!', __FILE__, __LINE__)
         case ('L')
            call g_logger%fatal('Layered/interface calculation not implemented!', __FILE__, __LINE__)
         end select
      else
         select case (control_obj%calctype)
         case ('B')
            do i = 1, lattice_obj%nrec
               call lattice_obj%symbolic_atoms(i)%build_pot()
            end do
            if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham
            call hamiltonian_obj%build_bulkham()
         case ('S')
            do i = 1, lattice_obj%ntype
               call lattice_obj%symbolic_atoms(i)%build_pot()
            end do
            if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham
            call hamiltonian_obj%build_bulkham()
         case ('I')
            do i = 1, lattice_obj%ntype
               call lattice_obj%symbolic_atoms(i)%build_pot()
            end do
            if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham
            call hamiltonian_obj%build_bulkham()
            call hamiltonian_obj%build_locham()
         case ('L')
            ! B7.5: identical to ’S’ -- see the same clause in self.f90’s
            ! run_recursion for why the loop runs to ntype and there is no
            ! build_locham.
            do i = 1, lattice_obj%ntype
               call lattice_obj%symbolic_atoms(i)%build_pot()
            end do
            if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham
            call hamiltonian_obj%build_bulkham()
         end select
      end if

      recursion_obj = recursion(hamiltonian_obj, energy_obj, sparse(hamiltonian_obj))
      if (stochastic_moments) call recursion_obj%compute_moments_stochastic()
      dos_obj = dos(recursion_obj, energy_obj)
      green_obj = green(dos_obj)
      bands_obj = bands(green_obj)
      if (.not. energy_mesh_before_hamiltonian) call energy_obj%e_mesh()
   end subroutine prepare_post_processing_stack

   !> @brief Run the intersite recursion pass that produces the moments for G_ij.
   module subroutine run_intersite_moments(control_obj, recursion_obj)
      type(control), intent(in) :: control_obj
      type(recursion), intent(inout) :: recursion_obj

      select case (control_obj%recur)
      case ('block')
         call recursion_obj%recur_b_ij()
      case ('chebyshev')
         call recursion_obj%chebyshev_recur_ij()
      end select
   end subroutine run_intersite_moments

   !> @brief k-space band-structure post-processing entry point.
   !> @details The routine builds a bulk-only pre-processing stack directly
   !>          (bravais, structb, atomlist, bulkmat, potential build,
   !>          Hamiltonian build). It does not support surface, impurity, or
   !>          interface calctypes. It then calls
   !>          reciprocal%calculate_band_structure along an automatic
   !>          high-symmetry path and writes band_structure.dat.
   !> @param[in] this Calculation object. fname selects the namelist input.
   module subroutine post_processing_band_structure(this)
      class(calculation), intent(in) :: this
      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(charge), target :: charge_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(reciprocal), target :: reciprocal_obj
      integer :: i

      control_obj = control(this%fname)
      lattice_obj = lattice(control_obj)
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      call lattice_obj%structb(.true.)
      call lattice_obj%atomlist()
      charge_obj = charge(lattice_obj)
      call charge_obj%bulkmat()
      hamiltonian_obj = hamiltonian(charge_obj)
      do i = 1, lattice_obj%nrec
         call lattice_obj%symbolic_atoms(i)%build_pot()
      end do
      if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham()
      call hamiltonian_obj%build_bulkham()
      reciprocal_obj = reciprocal(hamiltonian_obj)
      call reciprocal_obj%calculate_band_structure(hamiltonian_obj, 'auto', reciprocal_obj%nk_per_segment, 'band_structure.dat')
   end subroutine post_processing_band_structure

   !> @brief Bloch spectral function A(k,E) post-processing (milestone B3).
   !> @details Builds the same stack as post_processing_band_structure (the BSF is a
   !>          thin consumer of the B2 k-space Green’s function on the band path),
   !>          then delegates to reciprocal%calculate_bsf. Broadening is the
   !>          &reciprocal green_eta; the energy grid/range are n_energy_points /
   !>          dos_energy_range; the path density is nk_per_segment. With sigma=0 the
   !>          resolvent equals backend E; a non-zero Sigma provider (B8/B10) would
   !>          broaden A(k,E) through the same routine unchanged.
   module subroutine post_processing_bsf(this)
      class(calculation), intent(in) :: this
      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(charge), target :: charge_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(reciprocal), target :: reciprocal_obj
      integer :: i

      control_obj = control(this%fname)
      lattice_obj = lattice(control_obj)
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      call lattice_obj%structb(.true.)
      call lattice_obj%atomlist()
      charge_obj = charge(lattice_obj)
      call charge_obj%bulkmat()
      hamiltonian_obj = hamiltonian(charge_obj)
      do i = 1, lattice_obj%nrec
         call lattice_obj%symbolic_atoms(i)%build_pot()
      end do
      if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham()
      call hamiltonian_obj%build_bulkham()
      reciprocal_obj = reciprocal(hamiltonian_obj)
      call reciprocal_obj%calculate_bsf('bsf.dat')
   end subroutine post_processing_bsf

   !> @brief k-space density-of-states post-processing entry point.
   !> @details The routine builds the same bulk-only stack as
   !>          post_processing_band_structure. It then calls
   !>          reciprocal%calculate_density_of_states with the method, energy
   !>          range, temperature, and Fermi-search settings from the
   !>          &reciprocal namelist, and writes dos_kspace.dat.
   !> @param[in] this Calculation object. fname selects the namelist input.
   module subroutine post_processing_density_of_states(this)
      class(calculation), intent(in) :: this
      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(charge), target :: charge_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(reciprocal), target :: reciprocal_obj
      integer :: i

      control_obj = control(this%fname)
      lattice_obj = lattice(control_obj)
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      call lattice_obj%structb(.true.)
      call lattice_obj%atomlist()
      charge_obj = charge(lattice_obj)
      call charge_obj%bulkmat()
      hamiltonian_obj = hamiltonian(charge_obj)
      do i = 1, lattice_obj%nrec
         call lattice_obj%symbolic_atoms(i)%build_pot()
      end do
      if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham()
      call hamiltonian_obj%build_bulkham()
      reciprocal_obj = reciprocal(hamiltonian_obj)
      call reciprocal_obj%calculate_density_of_states(hamiltonian_obj, reciprocal_obj%n_energy_points, &
           reciprocal_obj%dos_energy_range, reciprocal_obj%dos_method, reciprocal_obj%gaussian_sigma, &
           reciprocal_obj%temperature, reciprocal_obj%fermi_level, reciprocal_obj%total_electrons, &
           reciprocal_obj%auto_find_fermi, 'dos_kspace.dat')
   end subroutine post_processing_density_of_states

   !> @brief Dense k-space eigenpair export for Fermi-surface analysis.
   !> @details Rebuilds the converged bulk post-processing Hamiltonian, then
   !>          samples it on a mesh independent of the SCF mesh. The writer
   !>          always emits the complete BZ so that each k-point has its own
   !>          eigenvectors for later site/orbital/spin colouring.
   !> @param[in] this Calculation object. fname selects the namelist input.
   module subroutine post_processing_fermi_surface(this)
      class(calculation), intent(in) :: this
      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(charge), target :: charge_obj
      type(energy), target :: energy_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(reciprocal), target :: reciprocal_obj
      integer :: i
      integer :: fs_mesh(3)

      control_obj = control(this%fname)
      lattice_obj = lattice(control_obj)
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      call lattice_obj%structb(.true.)
      call lattice_obj%atomlist()
      charge_obj = charge(lattice_obj)
      call charge_obj%bulkmat()
      energy_obj = energy(lattice_obj)
      hamiltonian_obj = hamiltonian(charge_obj)
      do i = 1, lattice_obj%nrec
         call lattice_obj%symbolic_atoms(i)%build_pot()
      end do
      if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham()
      call hamiltonian_obj%build_bulkham()

      reciprocal_obj = reciprocal(hamiltonian_obj)
      ! With auto_find_fermi disabled, use the Fermi level loaded from the
      ! converged &energy state, just as the normal post-processing routes do.
      reciprocal_obj%fermi_level = energy_obj%fermi
      fs_mesh = reciprocal_obj%nk_mesh
      if (reciprocal_obj%fs_nk1 > 0) fs_mesh(1) = reciprocal_obj%fs_nk1
      if (reciprocal_obj%fs_nk2 > 0) fs_mesh(2) = reciprocal_obj%fs_nk2
      if (reciprocal_obj%fs_nk3 > 0) fs_mesh(3) = reciprocal_obj%fs_nk3
      if (any(fs_mesh <= 0)) then
         call g_logger%fatal('post_processing_fermi_surface: FS mesh dimensions must be positive.', __FILE__, __LINE__)
      end if
      reciprocal_obj%nk_mesh = fs_mesh

      ! A reduced mesh is useful for scalar integrals, but cannot provide the
      ! per-k eigenvectors needed by a Fermi-surface payload. Keep this route
      ! deterministic and full-BZ even if the compatibility knob is set.
      if (reciprocal_obj%fs_use_symmetry_reduction) then
         call g_logger%warning('post_processing_fermi_surface: fs_use_symmetry_reduction is ignored; using the full BZ mesh.', &
                               __FILE__, __LINE__)
      end if
      reciprocal_obj%use_symmetry_reduction = .false.
      call reciprocal_obj%generate_mp_mesh()
      call reciprocal_obj%write_kspace_eigenpairs()
   end subroutine post_processing_fermi_surface

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief B2 validation driver (C1/C3 + C2): cross-check the k-space Lehmann
   !>        on-site Green’s function against the real-space recursion route.
   !> @details Runs BOTH routes on the same converged potential and compares the
   !>        on-site density of states rho(E) = -1/pi * Im Tr G_ii(E), built from
   !>        the SAME green%gij on-site block for each route (convention 5: the
   !>        full-resolvent nb x nb block, not the -i*pi*rho spectral form). This
   !>        is the route-agnosticism acceptance for backend E: the strict-Lehmann
   !>        k-space engine must reproduce the recursion on-site spectral function
   !>        (C1 representation) with the B1 bond/phase convention (C3, on-site
   !>        phase = 1). The recursion route fills green%gij at the on-site
   !>        self-pair (ijpair(:,1)==ijpair(:,2)); we snapshot its on-site block,
   !>        then reciprocal%fill_green overwrites green%gij from the k-space
   !>        eigenpairs and we recompute the same quantities. Report-only: writes
   !>        kspace_green_c1.dat with both traces and difference/spectral-weight
   !>        metrics; the acceptance tolerance is a maintainer gate. Regression is
   !>        untouched -- this is a new post_processing key, off by default.
   !>
   !>        C2 (spin frame): the charge DOS -1/pi Im Tr G_ii is invariant under
   !>        G -> R^dagger G R, so it cannot see the spin frame at all. We ALSO
   !>        report the z-projected spin DOS m_z(E) = -1/pi Im[Tr G_uu - Tr G_dd],
   !>        which IS frame-sensitive, as a noncollinear cross-check. The key
   !>        finding (verified on Mn3Sn, 120 deg NC): the INTERSITE recursion
   !>        `recur_b_ij` never rotates to local axes (only the on-site DOS
   !>        recursion does), so the RS gij is stored in the GLOBAL spin frame.
   !>        Backend E therefore fills the global-frame block directly and m_z from
   !>        both routes agrees to ~4e-4 (for an in-plane moment both are ~0). If a
   !>        LOCAL-frame comparison were wanted, BOTH routes’ G would have to be
   !>        rotated by the same rotmag_loc primitive -- rotating only one side is
   !>        wrong (it drove m_z diff to ~20). We still restore the global ee
   !>        before the k-space build (rotate_from_local_axis) as a safety net in
   !>        case any on-site recursion left hamiltonian%ee rotated.
   !---------------------------------------------------------------------------
   module subroutine post_processing_kspace_green(this)
      use sigma_provider_mod, only: sigma_zero
      class(calculation), intent(in) :: this
      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(energy), target :: energy_obj
      type(charge), target :: charge_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(recursion), target :: recursion_obj
      type(green), target :: green_obj
      type(dos), target :: dos_obj
      type(bands), target :: bands_obj
      type(mix), target :: mix_obj
      type(reciprocal), target :: reciprocal_obj
      type(sigma_zero) :: sigma
      real(rp), allocatable :: dos_rs(:), dos_le(:), mz_rs(:), mz_le(:)
      complex(rp), allocatable :: blk_rs(:, :, :), blk_le(:, :, :)
      complex(rp), allocatable :: gij_rs(:, :, :, :), gij_le(:, :, :, :), gij_dyson(:, :, :, :)
      real(rp), allocatable :: max_pair_diff(:)
      integer :: i, ie, ipair, ne, nblk, norb_h, npair, p0, newunit, gfunit
      real(rp) :: max_abs_diff, rms_diff, int_rs, int_le, de
      real(rp) :: max_mz_diff, max_blk_diff, max_dyson_diff
      real(rp) :: onsite_mom(3)

      call prepare_post_processing_stack(this, .false., .true., .true., .false., control_obj, lattice_obj, &
                                         charge_obj, mix_obj, energy_obj, hamiltonian_obj, recursion_obj, dos_obj, &
                                         green_obj, bands_obj)
      do i = 1, lattice_obj%ntype
         call lattice_obj%symbolic_atoms(i)%predls(lattice_obj%wav*ang2au)
      end do

      ! --- Real-space recursion route: fill green%gij (on-site = self pair). ---
      call run_intersite_moments(control_obj, recursion_obj)
      call green_obj%calculate_intersite_gf()

      ! Locate an on-site self-pair (ijpair(p,1) == ijpair(p,2)) to read G_ii from.
      ! Prefer one whose moment is OFF the global z axis so the C2 local-frame
      ! rotation is actually exercised (an on-site block has a single, unambiguous
      ! moment -- the intersite i/=j frame question is deliberately out of scope).
      ! Falls back to the first on-site pair (e.g. collinear bcc Fe, all along z).
      p0 = 0
      do i = 1, lattice_obj%njij
         if (lattice_obj%ijpair(i, 1) == lattice_obj%ijpair(i, 2)) then
            if (p0 == 0) p0 = i
            onsite_mom = lattice_obj%symbolic_atoms(i)%potential%mom(1:3)
            if (onsite_mom(1)**2 + onsite_mom(2)**2 > 1.0e-8_rp) then
               p0 = i
               exit
            end if
         end if
      end do
      if (p0 == 0) then
         call g_logger%fatal('[calculation.post_processing_kspace_green]: no on-site self-pair '// &
                             '(ijpair(:,1)==ijpair(:,2)) in the &lattice ijpair list; add e.g. ijpair(1,:)=1,1.', &
                             __FILE__, __LINE__)
      end if

      ne = green_obj%en%channels_ldos + 10
      nblk = size(green_obj%gij, 1)
      npair = size(green_obj%gij, 4)
      norb_h = nblk/2   ! spin-up orbitals 1..norb_h, spin-down norb_h+1..nblk
      allocate (dos_rs(ne), dos_le(ne), mz_rs(ne), mz_le(ne))
      allocate (blk_rs(nblk, nblk, ne), blk_le(nblk, nblk, ne))
      allocate (gij_rs(nblk, nblk, ne, npair), gij_le(nblk, nblk, ne, npair), &
                gij_dyson(nblk, nblk, ne, npair), max_pair_diff(npair))
      blk_rs = green_obj%gij(:, :, :, p0)
      gij_rs = green_obj%gij
      call onsite_dos_mz(blk_rs, norb_h, dos_rs, mz_rs)

      ! Restore the global-frame ee: the recursion route rotates hamiltonian%ee in
      ! place under local_axis and never calls rotate_from_local_axis, so H(k)
      ! must be rebuilt from the restored global representation (no-op collinear).
      onsite_mom = lattice_obj%symbolic_atoms(p0)%potential%mom(1:3)
      call hamiltonian_obj%rotate_from_local_axis(onsite_mom)

      ! --- k-space Lehmann route: overwrite green%gij from the eigenpairs. ---
      reciprocal_obj = reciprocal(hamiltonian_obj)
      call reciprocal_obj%generate_mp_mesh()   ! full unreduced BZ mesh (backend E requirement)
      call reciprocal_obj%fill_green(green_obj, sigma)
      blk_le = green_obj%gij(:, :, :, p0)
      gij_le = green_obj%gij
      call onsite_dos_mz(blk_le, norb_h, dos_le, mz_le)

      ! --- Backend D (Dyson, Sigma=0) cross-check on the real H(k). ---------------
      ! The permanent B2.4 invariant ″D with Sigma=0 == E″ over EVERY pair (the
      ! full gij array, so the intersite e^{ik.dR} phase is exercised, not just
      ! the on-site block). Both routes use S=I (backend E’s zheev is orthonormal,
      ! so backend D inverts z*I - H(k)); the difference is solver-tolerance ripple.
      gij_le = green_obj%gij
      reciprocal_obj%green_backend = 'dyson'
      call reciprocal_obj%fill_green(green_obj, sigma)
      gij_dyson = green_obj%gij
      max_dyson_diff = maxval(abs(green_obj%gij - gij_le))
      green_obj%gij = gij_le   ! restore backend-E result for the report below

      ! Direct G validation is pair-wise.  Keep the full complex block metric,
      ! and also write representative spin-diagonal matrix elements below so a
      ! campaign can compare G_ii(z) and G_ij(z), not only a derived DOS.
      max_pair_diff = 0.0_rp
      do ipair = 1, npair
         max_pair_diff(ipair) = maxval(abs(gij_le(:, :, :, ipair) - gij_rs(:, :, :, ipair)))
      end do

      ! --- Compare + report (report-only; tolerance is a maintainer gate). ---
      max_abs_diff = 0.0_rp; rms_diff = 0.0_rp; int_rs = 0.0_rp; int_le = 0.0_rp
      max_mz_diff = 0.0_rp; max_blk_diff = 0.0_rp
      do ie = 1, ne
         max_abs_diff = max(max_abs_diff, abs(dos_le(ie) - dos_rs(ie)))
         max_mz_diff = max(max_mz_diff, abs(mz_le(ie) - mz_rs(ie)))
         max_blk_diff = max(max_blk_diff, maxval(abs(blk_le(:, :, ie) - blk_rs(:, :, ie))))
         rms_diff = rms_diff + (dos_le(ie) - dos_rs(ie))**2
      end do
      rms_diff = sqrt(rms_diff/real(ne, rp))
      do ie = 2, ne
         de = green_obj%en%ene(ie) - green_obj%en%ene(ie - 1)
         int_rs = int_rs + 0.5_rp*de*(dos_rs(ie) + dos_rs(ie - 1))
         int_le = int_le + 0.5_rp*de*(dos_le(ie) + dos_le(ie - 1))
      end do

      if (rank == 0) then
         open (newunit=newunit, file='kspace_green_c1.dat', status='replace', action='write')
         write (newunit, '(A)') '# B2 C1/C3/C2 validation: on-site G_ii(E), recursion (RS) vs k-space Lehmann'
         write (newunit, '(A)') '# dos = -1/pi Im Tr G_ii ; mz = -1/pi Im (Tr G_uu - Tr G_dd)  (global z; gij stored global frame)'
         write (newunit, '(A,I0,A,I0,A,I0,A,I0,A)') '# k-mesh ', reciprocal_obj%nk_mesh(1), ' x ', &
            reciprocal_obj%nk_mesh(2), ' x ', reciprocal_obj%nk_mesh(3), '  (', reciprocal_obj%nk_total, ' points)'
         write (newunit, '(A,ES14.6,A)') '# green_eta = ', reciprocal_obj%green_eta, ' Ry'
         write (newunit, '(A,3F10.6)') '# on-site moment = ', onsite_mom
         write (newunit, '(A,ES14.6)') '# C1  max|dos_lehmann - dos_rs|      = ', max_abs_diff
         write (newunit, '(A,ES14.6)') '#     rms dos_lehmann - dos_rs       = ', rms_diff
         write (newunit, '(A,ES14.6)') '# C2  max|mz_lehmann - mz_rs|        = ', max_mz_diff
         write (newunit, '(A,ES14.6)') '# C2  max|G_ii^lehmann - G_ii^rs|    = ', max_blk_diff
         write (newunit, '(A,ES14.6)') '# B2.4 max|gij^dyson - gij^lehmann|  = ', max_dyson_diff
         do ipair = 1, npair
            write (newunit, '(A,I0,A,2(1x,I0),A,ES14.6)') '# direct G pair ', ipair, ' (', &
               lattice_obj%ijpair(ipair, 1), lattice_obj%ijpair(ipair, 2), ') max|G_RS-G_Lehmann| = ', max_pair_diff(ipair)
         end do
         write (newunit, '(A,ES14.6,A,ES14.6)') '# spectral weight   RS = ', int_rs, '   Lehmann = ', int_le
         write (newunit, '(A)') '#     E(Ry)          dos_rs        dos_lehmann       dos_diff          mz_rs         mz_lehmann'
         do ie = 1, ne
            write (newunit, '(6ES16.8)') green_obj%en%ene(ie), dos_rs(ie), dos_le(ie), &
               dos_le(ie) - dos_rs(ie), mz_rs(ie), mz_le(ie)
         end do
         close (newunit)
         call g_logger%info('[kspace_green] C1 on-site DOS cross-check: max|diff|='// &
                            trim(real2str(max_abs_diff))//' rms='//trim(real2str(rms_diff))// &
                            ' weight_RS='//trim(real2str(int_rs))//' weight_Lehmann='//trim(real2str(int_le)), &
                            __FILE__, __LINE__)
         call g_logger%info('[kspace_green] C2 spin-structure cross-check (both routes '// &
                            'global frame): max|mz_diff|='//trim(real2str(max_mz_diff))// &
                            ' max|block_diff|='//trim(real2str(max_blk_diff))// &
                            ' (block_diff is single-element ripple at coarse mesh)', &
                            __FILE__, __LINE__)
         call g_logger%info('[kspace_green] B2.4 backend-D (Dyson, Sigma=0) == backend-E '// &
                            '(Lehmann) cross-check: max|gij_dyson - gij_lehmann|='// &
                            trim(real2str(max_dyson_diff))//' (solver-tolerance invariant)', &
                            __FILE__, __LINE__)

         ! Machine-readable direct-G record.  The full-block maxima above are
         ! the acceptance metric; these selected spin-diagonal elements make
         ! onsite G_ii(z) and intersite G_ij(z) auditable at several energies.
         open (newunit=gfunit, file='kspace_green_gf.dat', status='replace', action='write')
         write (gfunit, '(A)') '# VAL-05 direct Green-function samples: RS, Lehmann, and Dyson(Sigma=0)'
         write (gfunit, '(A)') '# columns: pair i j E ReG11_RS ImG11_RS ReG11_Leh ImG11_Leh ReG11_Dyson ImG11_Dyson ReGdd_RS ImGdd_RS ReGdd_Leh ImGdd_Leh'
         do ipair = 1, npair
            do ie = 1, ne
               write (gfunit, '(2(I0,1x),11(1x,ES24.16))') lattice_obj%ijpair(ipair, 1), &
                  lattice_obj%ijpair(ipair, 2), green_obj%en%ene(ie), &
                  real(gij_rs(1, 1, ie, ipair)), aimag(gij_rs(1, 1, ie, ipair)), &
                  real(gij_le(1, 1, ie, ipair)), aimag(gij_le(1, 1, ie, ipair)), &
                  real(gij_dyson(1, 1, ie, ipair)), aimag(gij_dyson(1, 1, ie, ipair)), &
                  real(gij_rs(norb_h + 1, norb_h + 1, ie, ipair)), aimag(gij_rs(norb_h + 1, norb_h + 1, ie, ipair)), &
                  real(gij_le(norb_h + 1, norb_h + 1, ie, ipair)), aimag(gij_le(norb_h + 1, norb_h + 1, ie, ipair))
            end do
         end do
         close (gfunit)
      end if

      deallocate (dos_rs, dos_le, mz_rs, mz_le, blk_rs, blk_le, gij_rs, gij_le, gij_dyson, max_pair_diff)
   end subroutine post_processing_kspace_green

   !> @brief On-site charge DOS and z-projected spin DOS from a G_ii block.
   !> @details dos(E)=-1/pi Im Tr G ; mz(E)=-1/pi Im (Tr G_uu - Tr G_dd), with the
   !>          spin-up orbitals in rows 1..norb and spin-down in norb+1..2*norb.
   subroutine onsite_dos_mz(blk, norb, dos, mz)
      complex(rp), intent(in) :: blk(:, :, :)
      integer, intent(in) :: norb
      real(rp), intent(out) :: dos(:), mz(:)
      integer :: ie, ib, ne
      complex(rp) :: tr_up, tr_dn

      ne = size(blk, 3)
      do ie = 1, ne
         tr_up = (0.0_rp, 0.0_rp); tr_dn = (0.0_rp, 0.0_rp)
         do ib = 1, norb
            tr_up = tr_up + blk(ib, ib, ie)
            tr_dn = tr_dn + blk(ib + norb, ib + norb, ie)
         end do
         dos(ie) = -aimag(tr_up + tr_dn)/pi
         mz(ie) = -aimag(tr_up - tr_dn)/pi
      end do
   end subroutine onsite_dos_mz

end submodule calculation_reciprocal
