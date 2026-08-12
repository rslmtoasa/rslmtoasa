submodule (reciprocal_mod) reciprocal_lifecycle
   use magnetic_representation_mod, only: gbt_single_q
   implicit none

contains

   !> @brief Print a root-rank informational message with source location.
   !> @param[in] message Message text.
   !> @param[in] file_name Source file name for diagnostics.
   !> @param[in] line_no Source line number for diagnostics.
   module subroutine root_info(message, file_name, line_no)
      character(len=*), intent(in) :: message, file_name
      integer, intent(in) :: line_no

      if (rank == 0) call g_logger%info(message, file_name, line_no)
   end subroutine root_info

   !> @brief Construct a reciprocal-space object from an initialized Hamiltonian.
   !> @details Wires reciprocal state to the Hamiltonian, lattice, and control
   !>          objects, restores defaults, and reads the &reciprocal namelist.
   !> @param[in] hamiltonian_obj Hamiltonian object providing real-space blocks and lattice state.
   !> @return Initialized reciprocal-space object.
   module function constructor(hamiltonian_obj) result(obj)
      type(reciprocal) :: obj
      type(hamiltonian), target, intent(in) :: hamiltonian_obj

      obj%hamiltonian => hamiltonian_obj
      obj%lattice => hamiltonian_obj%lattice
      obj%control => hamiltonian_obj%lattice%control

      call obj%restore_to_default()
      call obj%build_from_file()  ! Read parameters from input.nml
      call obj%generate_reciprocal_vectors()
      call obj%set_basis_sizes()
      call obj%symmetry_analysis%initialize(obj%lattice)
      
      ! Auto-calculate total electrons from valence if not set in input.
      ! For bulk systems nbulk can be zero; nbulk_bulk is the correct valence span.
      if (obj%total_electrons <= 1.0e-3_rp) then
         if (obj%lattice%nbulk_bulk > 0) then
            obj%total_electrons = real(sum(obj%lattice%symbolic_atoms(1:obj%lattice%nbulk_bulk)%element%valence), rp)
         else if (obj%lattice%nrec > 0) then
            obj%total_electrons = real(sum(obj%lattice%symbolic_atoms(1:obj%lattice%nrec)%element%valence), rp)
            if (rank == 0) call g_logger%warning('reciprocal%constructor: nbulk_bulk<=0, using nrec span for total_electrons.', __FILE__, __LINE__)
         end if
         call root_info('reciprocal%constructor: Auto-calculated total_electrons = ' // trim(real2str(obj%total_electrons)) // ' from valence', __FILE__, __LINE__)
      end if
   end function constructor

   !> @brief Finalize a reciprocal-space object.
   !> @details Releases k-point, reciprocal-vector, Hamiltonian, overlap, band,
   !>          DOS, tetrahedron, symmetry-map, and projection arrays.
   !> @param[inout] this Reciprocal object being finalized.
   module subroutine destructor(this)
      type(reciprocal) :: this
      if (allocated(this%execution_backend)) then
         call this%execution_backend%release()
         deallocate(this%execution_backend)
      end if
#ifdef USE_SAFE_ALLOC
      if (allocated(this%k_points)) call g_safe_alloc%deallocate('reciprocal.k_points', this%k_points)
      if (allocated(this%k_weights)) call g_safe_alloc%deallocate('reciprocal.k_weights', this%k_weights)
      if (allocated(this%k_l2g_map)) deallocate(this%k_l2g_map)
      if (allocated(this%k_g2l_map)) deallocate(this%k_g2l_map)
      if (allocated(this%hk_bulk)) call g_safe_alloc%deallocate('reciprocal.hk_bulk', this%hk_bulk)
      if (allocated(this%hk_so)) call g_safe_alloc%deallocate('reciprocal.hk_so', this%hk_so)
      if (allocated(this%hk_total)) call g_safe_alloc%deallocate('reciprocal.hk_total', this%hk_total)
      if (allocated(this%sk_overlap)) call g_safe_alloc%deallocate('reciprocal.sk_overlap', this%sk_overlap)
      if (allocated(this%basis_size)) call g_safe_alloc%deallocate('reciprocal.basis_size', this%basis_size)
      if (allocated(this%k_path)) call g_safe_alloc%deallocate('reciprocal.k_path', this%k_path)
      if (allocated(this%k_labels)) call g_safe_alloc%deallocate('reciprocal.k_labels', this%k_labels)
      if (allocated(this%k_distances)) call g_safe_alloc%deallocate('reciprocal.k_distances', this%k_distances)
      if (allocated(this%eigenvalues)) call g_safe_alloc%deallocate('reciprocal.eigenvalues', this%eigenvalues)
      if (allocated(this%eigenvalues_path)) call g_safe_alloc%deallocate('reciprocal.eigenvalues_path', this%eigenvalues_path)
      if (allocated(this%eigenvectors)) call g_safe_alloc%deallocate('reciprocal.eigenvectors', this%eigenvectors)
      if (allocated(this%eigenvectors_path)) call g_safe_alloc%deallocate('reciprocal.eigenvectors_path', this%eigenvectors_path)
      ! DOS arrays
      if (allocated(this%dos_energy_grid)) call g_safe_alloc%deallocate('reciprocal.dos_energy_grid', this%dos_energy_grid)
      if (allocated(this%total_dos)) call g_safe_alloc%deallocate('reciprocal.total_dos', this%total_dos)
      if (allocated(this%total_nos)) call g_safe_alloc%deallocate('reciprocal.total_nos', this%total_nos)
      if (allocated(this%projected_dos)) call g_safe_alloc%deallocate('reciprocal.projected_dos', this%projected_dos)
      if (allocated(this%band_moments)) call g_safe_alloc%deallocate('reciprocal.band_moments', this%band_moments)
      if (allocated(this%dos_mx_tot)) call g_safe_alloc%deallocate('reciprocal.dos_mx_tot', this%dos_mx_tot)
      if (allocated(this%dos_my_tot)) call g_safe_alloc%deallocate('reciprocal.dos_my_tot', this%dos_my_tot)
      if (allocated(this%dos_mz_tot)) call g_safe_alloc%deallocate('reciprocal.dos_mz_tot', this%dos_mz_tot)
      if (allocated(this%projected_dos_moments)) call g_safe_alloc%deallocate('reciprocal.projected_dos_moments', this%projected_dos_moments)
      if (allocated(this%tetrahedra)) call g_safe_alloc%deallocate('reciprocal.tetrahedra', this%tetrahedra)
      if (allocated(this%tetrahedron_volumes)) call g_safe_alloc%deallocate('reciprocal.tetrahedron_volumes', this%tetrahedron_volumes)
      if (allocated(this%ham_vec_type)) call g_safe_alloc%deallocate('reciprocal.ham_vec_type', this%ham_vec_type)
      if (allocated(this%ham_vec_type_direct)) call g_safe_alloc%deallocate('reciprocal.ham_vec_type_direct', this%ham_vec_type_direct)
      if (allocated(this%mesh_cache_q)) deallocate(this%mesh_cache_q)
#else
      if (allocated(this%k_points)) deallocate (this%k_points)
      if (allocated(this%k_weights)) deallocate (this%k_weights)
      if (allocated(this%k_l2g_map)) deallocate(this%k_l2g_map)
      if (allocated(this%k_g2l_map)) deallocate(this%k_g2l_map)
      if (allocated(this%hk_bulk)) deallocate (this%hk_bulk)
      if (allocated(this%hk_so)) deallocate (this%hk_so)
      if (allocated(this%hk_total)) deallocate (this%hk_total)
      if (allocated(this%sk_overlap)) deallocate (this%sk_overlap)
      if (allocated(this%basis_size)) deallocate (this%basis_size)
      if (allocated(this%k_path)) deallocate (this%k_path)
      if (allocated(this%k_labels)) deallocate (this%k_labels)
      if (allocated(this%k_distances)) deallocate (this%k_distances)
      if (allocated(this%eigenvalues)) deallocate(this%eigenvalues)
      if (allocated(this%eigenvectors)) deallocate(this%eigenvectors)
      ! DOS arrays
      if (allocated(this%dos_energy_grid)) deallocate (this%dos_energy_grid)
      if (allocated(this%total_dos)) deallocate (this%total_dos)
      if (allocated(this%total_nos)) deallocate (this%total_nos)
      if (allocated(this%projected_dos)) deallocate (this%projected_dos)
      if (allocated(this%band_moments)) deallocate (this%band_moments)
      if (allocated(this%dos_mx_tot)) deallocate (this%dos_mx_tot)
      if (allocated(this%dos_my_tot)) deallocate (this%dos_my_tot)
      if (allocated(this%dos_mz_tot)) deallocate (this%dos_mz_tot)
      if (allocated(this%projected_dos_moments)) deallocate (this%projected_dos_moments)
      if (allocated(this%tetrahedra)) deallocate (this%tetrahedra)
      if (allocated(this%tetrahedron_volumes)) deallocate (this%tetrahedron_volumes)
      if (allocated(this%ham_vec_type)) deallocate (this%ham_vec_type)
      if (allocated(this%ham_vec_type_direct)) deallocate (this%ham_vec_type_direct)
      if (allocated(this%full_to_irred_k)) deallocate(this%full_to_irred_k)
      if (allocated(this%irred_to_full_k)) deallocate(this%irred_to_full_k)
      if (allocated(this%mesh_cache_q)) deallocate(this%mesh_cache_q)
#endif
   end subroutine destructor

   !> @brief Reset reciprocal-space settings and owned storage to defaults.
   !> @details Restores k-mesh, solver-mode, DOS, band-path, symmetry, and
   !>          diagnostic options to baseline values and clears allocatables.
   !> @param[inout] this Reciprocal object to reset.
   module subroutine restore_to_default(this)
      class(reciprocal), intent(inout) :: this

      ! Default k-point mesh settings
      this%nk_mesh = [8, 8, 8]  ! Default 8x8x8 mesh
      this%nk_total = 0
      this%nk_local = 0
      this%k_start = 1
      this%k_end = 0
      this%k_mesh_distributed_active = .false.
      call this%k_workset%restore_to_default()
      this%use_time_reversal = .true.
      this%strict_symmetry_checks = .false.
      this%dump_symmetry_kmap = .false.
      this%tetra_symmetry_mode = 'irreducible_native'
      this%k_offset = [0.0_rp, 0.0_rp, 0.0_rp]  ! No shift by default
      this%include_so = .false.
      this%max_orbs = nb
      this%cached_operator_generation = -1

      ! WP8: full BZ stays the default/oracle for finite-q GBT (WP0). Opt-in
      ! to 'little_group' or 'little_group_common' via &reciprocal.
      this%q_symmetry_policy = 'full_bz'
      this%mesh_cache_valid = .false.
      this%mesh_cache_dims = [0, 0, 0]
      this%mesh_cache_offset = [0.0_rp, 0.0_rp, 0.0_rp]
      this%mesh_cache_lattice = 0.0_rp
      this%mesh_cache_policy = ''
      if (allocated(this%mesh_cache_q)) deallocate(this%mesh_cache_q)

   ! By default suppress internal verbose prints (can be enabled by user)
   this%suppress_internal_logs = .true.

      ! Initialize reciprocal lattice to zero
      this%reciprocal_vectors = 0.0_rp
      this%reciprocal_volume = 0.0_rp

      ! Default DOS settings
      ! All energy values in Ry (consistent with Hamiltonian)
      this%n_energy_points = 1000
      this%dos_energy_range = [-1.0_rp, 1.0_rp]  ! Default energy range in Ry
      this%dos_method = 'tetrahedron'  ! Default to tetrahedron method
      this%gaussian_sigma = 0.01_rp  ! Default Gaussian smearing in Ry
      this%temperature = 300.0_rp  ! Default temperature in Kelvin
      this%fermi_level = 0.0_rp  ! Default Fermi level
      this%total_electrons = 0.0_rp  ! 0 = auto-calculate from valence in constructor
      this%canonical_electron_count = 0.0_rp
      this%canonical_band_energy = 0.0_rp
      this%canonical_weight_sum = 0.0_rp
      this%canonical_energy_valid = .false.
      this%auto_find_fermi = .true.  ! Auto-find Fermi level from eigenvalue occupations
      this%reciprocal_mode = 'ham_only'
      this%kspace_ham_order = 'auto'
      this%kanpur_diagnostics = .true.
      this%gamma_bounds_diagnostics = .false.
      this%hall_diag_experimental = .false.
      this%reciprocal_tile_size = 16
      call this%workspace%restore_to_default()
      this%n_sites = 0
      this%n_orb_types = 4  ! s, p, d, f
      this%n_spin_components = 2  ! RS-LMTO basis is always spin-polarized
      this%n_tetrahedra = 0

      ! Default k-path settings
      this%auto_kpath = .true.  ! Use automatic k-path generation by default
      this%nk_per_segment = 40  ! Default 40 k-points per segment
      this%override_space_group = 0  ! 0 = auto-detect
      this%custom_kpath_spec = ''  ! Empty = use automatic
      this%use_symmetry_reduction = .true.  ! Use symmetry reduction by default

      ! k-space Green's-function engine (milestone B2, reciprocal_green).
      ! green_eta is the retarded broadening for z = E + i*green_eta (namelist key
      ! in &reciprocal). Default set to 0.02 Ry by gate G-B2-2 (Anders, 2026-07-16):
      ! the documented working point where the Lehmann/Dyson exchange J k-converges
      ! to ~1% by a 16^3 mesh (docs/dev/reciprocal_green_convergence.md). Supersedes
      ! the 0.01 Ry placeholder from gate G-B2-1.
      this%green_eta = 0.02_rp        ! Ry; gate G-B2-2 working point (was 0.01, G-B2-1)
      this%green_backend = 'lehmann'  ! backend E (Sigma = 0) by default
   end subroutine restore_to_default

   !> @brief Invalidate every spectrum-derived cache while retaining the k mesh.
   !> @details Frozen-magnon and SCF probes change the Hamiltonian repeatedly.
   !>          Keeping k_points/k_weights is safe, but eigenpairs, DOS arrays,
   !>          projections, tetrahedra, and their canonical energy are not.
   module subroutine invalidate_spectral_cache(this)
      class(reciprocal), intent(inout) :: this

#ifdef USE_SAFE_ALLOC
      if (allocated(this%hk_bulk)) call g_safe_alloc%deallocate('reciprocal.hk_bulk', this%hk_bulk)
      if (allocated(this%hk_so)) call g_safe_alloc%deallocate('reciprocal.hk_so', this%hk_so)
      if (allocated(this%hk_total)) call g_safe_alloc%deallocate('reciprocal.hk_total', this%hk_total)
      if (allocated(this%sk_overlap)) call g_safe_alloc%deallocate('reciprocal.sk_overlap', this%sk_overlap)
      if (allocated(this%eigenvalues)) call g_safe_alloc%deallocate('reciprocal.eigenvalues', this%eigenvalues)
      if (allocated(this%eigenvalues_path)) call g_safe_alloc%deallocate('reciprocal.eigenvalues_path', this%eigenvalues_path)
      if (allocated(this%eigenvectors)) call g_safe_alloc%deallocate('reciprocal.eigenvectors', this%eigenvectors)
      if (allocated(this%eigenvectors_path)) call g_safe_alloc%deallocate('reciprocal.eigenvectors_path', this%eigenvectors_path)
      if (allocated(this%total_dos)) call g_safe_alloc%deallocate('reciprocal.total_dos', this%total_dos)
      if (allocated(this%total_nos)) call g_safe_alloc%deallocate('reciprocal.total_nos', this%total_nos)
      if (allocated(this%projected_dos)) call g_safe_alloc%deallocate('reciprocal.projected_dos', this%projected_dos)
      if (allocated(this%band_moments)) call g_safe_alloc%deallocate('reciprocal.band_moments', this%band_moments)
      if (allocated(this%dos_mx_tot)) call g_safe_alloc%deallocate('reciprocal.dos_mx_tot', this%dos_mx_tot)
      if (allocated(this%dos_my_tot)) call g_safe_alloc%deallocate('reciprocal.dos_my_tot', this%dos_my_tot)
      if (allocated(this%dos_mz_tot)) call g_safe_alloc%deallocate('reciprocal.dos_mz_tot', this%dos_mz_tot)
      if (allocated(this%projected_dos_moments)) call g_safe_alloc%deallocate('reciprocal.projected_dos_moments', this%projected_dos_moments)
      if (allocated(this%tetrahedra)) call g_safe_alloc%deallocate('reciprocal.tetrahedra', this%tetrahedra)
      if (allocated(this%tetrahedron_volumes)) call g_safe_alloc%deallocate('reciprocal.tetrahedron_volumes', this%tetrahedron_volumes)
#else
      if (allocated(this%hk_bulk)) deallocate(this%hk_bulk)
      if (allocated(this%hk_so)) deallocate(this%hk_so)
      if (allocated(this%hk_total)) deallocate(this%hk_total)
      if (allocated(this%sk_overlap)) deallocate(this%sk_overlap)
      if (allocated(this%eigenvalues)) deallocate(this%eigenvalues)
      if (allocated(this%eigenvalues_path)) deallocate(this%eigenvalues_path)
      if (allocated(this%eigenvectors)) deallocate(this%eigenvectors)
      if (allocated(this%eigenvectors_path)) deallocate(this%eigenvectors_path)
      if (allocated(this%total_dos)) deallocate(this%total_dos)
      if (allocated(this%total_nos)) deallocate(this%total_nos)
      if (allocated(this%projected_dos)) deallocate(this%projected_dos)
      if (allocated(this%band_moments)) deallocate(this%band_moments)
      if (allocated(this%dos_mx_tot)) deallocate(this%dos_mx_tot)
      if (allocated(this%dos_my_tot)) deallocate(this%dos_my_tot)
      if (allocated(this%dos_mz_tot)) deallocate(this%dos_mz_tot)
      if (allocated(this%projected_dos_moments)) deallocate(this%projected_dos_moments)
      if (allocated(this%tetrahedra)) deallocate(this%tetrahedra)
      if (allocated(this%tetrahedron_volumes)) deallocate(this%tetrahedron_volumes)
#endif

      this%canonical_electron_count = 0.0_rp
      this%canonical_band_energy = 0.0_rp
      this%canonical_weight_sum = 0.0_rp
      this%canonical_energy_valid = .false.
   end subroutine invalidate_spectral_cache

   !> @brief Invalidate H(k), eigensystem, DOS, and density projections after
   !>        any shared real-space operator rebuild.
   !> @details build_bulkham advances hamiltonian%operator_generation before
   !>          consuming q, cone/reference frames, and potential parameters.
   !>          A generation mismatch therefore invalidates every downstream
   !>          reciprocal object without relying on floating-point fingerprints.
   module subroutine invalidate_if_operator_changed(this, context_tag, changed)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in) :: context_tag
      logical, intent(out), optional :: changed
      logical :: mismatch

      mismatch = .false.
      if (associated(this%hamiltonian)) then
         mismatch = this%cached_operator_generation /= this%hamiltonian%operator_generation
      end if
      if (present(changed)) changed = mismatch
      if (.not. mismatch) return

      call this%invalidate_spectral_cache()
      call root_info(trim(context_tag)//': invalidated operator-dependent H(k), eigensystem, DOS, and density caches', &
                     __FILE__, __LINE__)
   end subroutine invalidate_if_operator_changed

   !> @brief Read the &reciprocal namelist and install reciprocal-space options.
   !> @details Parses k-mesh, Fourier mode, band/DOS controls, symmetry options,
   !>          tetrahedron settings, and diagnostics from this%control%fname.
   !> @param[inout] this Reciprocal object whose input-facing options are populated.
   !> @note This is an input boundary and may raise fatal diagnostics for invalid options.
   module subroutine build_from_file(this)
      class(reciprocal), intent(inout) :: this

      ! Reading process variables
      integer :: iostatus, funit

      ! Include namelist declarations for reciprocal module
      include 'include_codes/namelists/reciprocal.f90'
      include 'include_codes/namelists/kpath.f90'

      ! Save previous values (from defaults or previous read)
      nk1 = this%nk_mesh(1)
      nk2 = this%nk_mesh(2)
      nk3 = this%nk_mesh(3)
      k_offset_x = this%k_offset(1)
      k_offset_y = this%k_offset(2)
      k_offset_z = this%k_offset(3)
      use_symmetry_reduction = this%use_symmetry_reduction
      use_time_reversal = this%use_time_reversal
      strict_symmetry_checks = this%strict_symmetry_checks
      dump_symmetry_kmap = this%dump_symmetry_kmap
      tetra_symmetry_mode = this%tetra_symmetry_mode
      q_symmetry_policy = this%q_symmetry_policy
      use_shift = .false.  ! Derived from k_offset
      n_energy_points = this%n_energy_points
      dos_energy_min = this%dos_energy_range(1)
      dos_energy_max = this%dos_energy_range(2)
      gaussian_sigma = this%gaussian_sigma
      temperature = this%temperature
      dos_method = this%dos_method
      auto_find_fermi = this%auto_find_fermi
      suppress_internal_logs = this%suppress_internal_logs
      reciprocal_mode = this%reciprocal_mode
      kspace_ham_order = this%kspace_ham_order
      kanpur_diagnostics = this%kanpur_diagnostics
      gamma_bounds_diagnostics = this%gamma_bounds_diagnostics
      hall_diag_experimental = this%hall_diag_experimental
      green_eta = this%green_eta
      green_backend = this%green_backend

      ! K-path settings
      auto_kpath = this%auto_kpath
      nk_per_segment = this%nk_per_segment
      override_space_group = this%override_space_group
      custom_kpath_spec = this%custom_kpath_spec

      ! Read reciprocal namelist
      open (newunit=funit, file=this%control%fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//trim(this%control%fname)//' not found', __FILE__, __LINE__)
      end if

      read (funit, nml=reciprocal, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         ! Namelist not found or error - use defaults
         call root_info('reciprocal namelist not found in input file, using defaults', __FILE__, __LINE__)
      end if
      
      ! Read kpath namelist
      rewind(funit)
      read (funit, nml=kpath, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         ! Namelist not found or error - use defaults
         call root_info('kpath namelist not found in input file, using defaults', __FILE__, __LINE__)
      end if
      close (funit)

      ! Assign values back to type members
      this%nk_mesh = [nk1, nk2, nk3]
      this%k_offset = [k_offset_x, k_offset_y, k_offset_z]
      this%use_symmetry_reduction = use_symmetry_reduction
      this%use_time_reversal = use_time_reversal
      this%strict_symmetry_checks = strict_symmetry_checks
      this%dump_symmetry_kmap = dump_symmetry_kmap
      this%tetra_symmetry_mode = lower(trim(tetra_symmetry_mode))
      if (this%tetra_symmetry_mode /= 'full_expand_ref' .and. &
          this%tetra_symmetry_mode /= 'irreducible_native') then
         call g_logger%warning("reciprocal%build_from_file: tetra_symmetry_mode must be 'full_expand_ref' or 'irreducible_native'. Falling back to irreducible_native.", __FILE__, __LINE__)
         this%tetra_symmetry_mode = 'irreducible_native'
      end if
      this%q_symmetry_policy = lower(trim(q_symmetry_policy))
      if (this%q_symmetry_policy /= 'full_bz' .and. this%q_symmetry_policy /= 'little_group' .and. &
          this%q_symmetry_policy /= 'little_group_common') then
         call g_logger%warning("reciprocal%build_from_file: q_symmetry_policy must be 'full_bz', "// &
                               "'little_group', or 'little_group_common'. Falling back to full_bz.", &
                               __FILE__, __LINE__)
         this%q_symmetry_policy = 'full_bz'
      end if
      this%n_energy_points = n_energy_points
      this%dos_energy_range = [dos_energy_min, dos_energy_max]
      this%gaussian_sigma = gaussian_sigma
      this%temperature = temperature
      this%dos_method = dos_method
      this%auto_find_fermi = auto_find_fermi
      this%suppress_internal_logs = suppress_internal_logs
      this%reciprocal_mode = lower(trim(reciprocal_mode))
      this%kanpur_diagnostics = kanpur_diagnostics
      this%gamma_bounds_diagnostics = gamma_bounds_diagnostics
      this%hall_diag_experimental = hall_diag_experimental
      this%green_eta = green_eta
      this%green_backend = lower(trim(green_backend))
      if (this%green_backend /= 'lehmann' .and. this%green_backend /= 'dyson') then
         call g_logger%warning("reciprocal%build_from_file: green_backend must be " // &
            "'lehmann' or 'dyson'. Falling back to lehmann.", __FILE__, __LINE__)
         this%green_backend = 'lehmann'
      end if
      if (this%reciprocal_mode == 'generalized_overlap') then
         this%reciprocal_mode = 'generalized_overlap_proxy'
         call g_logger%warning("reciprocal_mode='generalized_overlap' is deprecated alias; using 'generalized_overlap_proxy'.", __FILE__, __LINE__)
      end if
      if (this%reciprocal_mode /= 'ham_only' .and. &
          this%reciprocal_mode /= 'generalized_overlap_proxy' .and. &
          this%reciprocal_mode /= 'generalized_overlap_kanpur') then
         call g_logger%warning("reciprocal_mode must be 'ham_only', 'generalized_overlap_proxy', or 'generalized_overlap_kanpur'. Falling back to ham_only.", __FILE__, __LINE__)
         this%reciprocal_mode = 'ham_only'
      end if
	      this%kspace_ham_order = lower(trim(kspace_ham_order))
	      if (this%kspace_ham_order == 'proper') then
	         this%kspace_ham_order = 'second'
	         call g_logger%warning("kspace_ham_order='proper' is deprecated; using 'second'.", __FILE__, __LINE__)
	      end if
	      if (this%kspace_ham_order /= 'first' .and. this%kspace_ham_order /= 'second' .and. &
	          this%kspace_ham_order /= 'auto') then
	         call g_logger%warning("kspace_ham_order must be 'auto', 'first', 'second', or deprecated alias 'proper'. Falling back to auto.", __FILE__, __LINE__)
	         this%kspace_ham_order = 'auto'
	      end if

      ! WP0: a finite-q spin spiral breaks the crystal symmetries used for an
      ! irreducible wedge (and, in general, time reversal).  Keep the chemical
      ! cell, but always integrate its *full* BZ.  The same enforcement is made
      ! immediately before an H(k) build because frozen-magnon sweeps mutate q.
      if (this%has_nonzero_q_gbt()) then
         this%use_symmetry_reduction = .false.
         this%use_time_reversal = .false.
         call root_info('reciprocal%build_from_file: nonzero-q GBT forces the full chemical BZ; '// &
                        'symmetry and time-reversal reduction disabled.', __FILE__, __LINE__)
      end if

      ! K-path settings
      this%auto_kpath = auto_kpath
      this%nk_per_segment = nk_per_segment
      this%override_space_group = override_space_group
      this%custom_kpath_spec = custom_kpath_spec

      ! Log what was read
      call root_info('reciprocal%build_from_file: Read k-mesh = ' // &
                     trim(int2str(nk1)) // ' x ' // trim(int2str(nk2)) // ' x ' // trim(int2str(nk3)), &
                     __FILE__, __LINE__)
      
      ! if (sum(abs(this%k_offset)) > 1.0e-8_rp) then
      !    call g_logger%info('reciprocal%build_from_file: k-offset = [' // &
      !                      trim(real2str(k_offset_x, '(F8.4)')) // ', ' // &
      !                      trim(real2str(k_offset_y, '(F8.4)')) // ', ' // &
      !                      trim(real2str(k_offset_z, '(F8.4)')) // ']', &
      !                      __FILE__, __LINE__)
      ! end if
      
      if (this%use_symmetry_reduction) then
         call root_info('reciprocal%build_from_file: Symmetry reduction enabled', __FILE__, __LINE__)
      end if
      if (this%use_time_reversal) then
         call root_info('reciprocal%build_from_file: use_time_reversal = true', __FILE__, __LINE__)
      else
         call root_info('reciprocal%build_from_file: use_time_reversal = false', __FILE__, __LINE__)
      end if
      if (this%strict_symmetry_checks) then
         call root_info('reciprocal%build_from_file: strict_symmetry_checks = true', __FILE__, __LINE__)
      else
         call root_info('reciprocal%build_from_file: strict_symmetry_checks = false', __FILE__, __LINE__)
      end if
      call root_info('reciprocal%build_from_file: tetra_symmetry_mode = ' // trim(this%tetra_symmetry_mode), __FILE__, __LINE__)
      call root_info('reciprocal%build_from_file: q_symmetry_policy = ' // trim(this%q_symmetry_policy), __FILE__, __LINE__)
      
      if (this%auto_kpath) then
         call root_info('reciprocal%build_from_file: Automatic k-path generation enabled', __FILE__, __LINE__)
      end if
      call root_info('reciprocal%build_from_file: reciprocal_mode = ' // trim(this%reciprocal_mode), __FILE__, __LINE__)
      call root_info('reciprocal%build_from_file: kspace_ham_order = ' // trim(this%kspace_ham_order), __FILE__, __LINE__)
   end subroutine build_from_file

   module logical function has_nonzero_q_gbt(this)
      class(reciprocal), intent(in) :: this
      real(rp), parameter :: q_tolerance = 1.0e-12_rp

      has_nonzero_q_gbt = .false.
      if (.not. associated(this%hamiltonian)) return
      has_nonzero_q_gbt = trim(this%hamiltonian%magnetic_representation) == gbt_single_q .and. &
                          maxval(abs(this%hamiltonian%q_ss)) > q_tolerance
   end function has_nonzero_q_gbt

   !> Despite the name (kept from WP0, where the reciprocal guards were scoped
   !> to finite q), every check below fires for *any* gbt_single_q run,
   !> including q=0: WP6a widened the overlap rejection and none of these
   !> terms has an audited GBT form at q=0 either. has_nonzero_q_gbt is the
   !> genuinely q-scoped predicate; use that one when q matters. Renaming this
   !> routine belongs to WP10's documentation pass.
   module subroutine validate_nonzero_q_gbt(this, context)
      class(reciprocal), intent(in) :: this
      character(len=*), intent(in) :: context
      real(rp), parameter :: term_tolerance = 1.0e-12_rp

      if (.not. associated(this%hamiltonian)) return
      if (trim(this%hamiltonian%magnetic_representation) /= gbt_single_q) return

      ! These operators have no audited finite-q bond-gauge transformation.
      ! Fail before reciprocal operator construction can produce a
      ! plausible-looking, but invalid, spectrum.
      if (this%hamiltonian%hubbard_v_check) then
         call g_logger%fatal(trim(context)//': GBT with intersite Hubbard-V is unsupported.', &
                             __FILE__, __LINE__)
      end if
      if (trim(this%reciprocal_mode) /= 'ham_only') then
         call g_logger%fatal(trim(context)//': GBT with reciprocal_mode='// &
                             trim(this%reciprocal_mode)//' is unsupported: the available generalized-overlap '// &
                             'representation is not a complete formal GBT metric; use ham_only.', &
                             __FILE__, __LINE__)
      end if
      if (allocated(this%hamiltonian%lsham)) then
         if (maxval(abs(this%hamiltonian%lsham)) > term_tolerance) then
            call g_logger%fatal(trim(context)//': GBT with SOC is unsupported.', &
                                __FILE__, __LINE__)
         end if
      end if
   end subroutine validate_nonzero_q_gbt

   !> @brief Enforce the nonzero-q GBT BZ policy before every H(k) build.
   !> @details Called from build_kspace_hamiltonian on every probe, so this is
   !>          the actual per-Hamiltonian-build enforcement point -- the
   !>          calculation.f90 sweep-level mesh calls only prepare a mesh this
   !>          routine must not silently discard. 'full_bz' (default)
   !>          reproduces the original WP0 behaviour exactly: force the full
   !>          chemical BZ every call. 'little_group' delegates to
   !>          ensure_kpoint_mesh, which is a no-op once the current q_ss's
   !>          mesh is already cached (set by the calculation.f90 per-q
   !>          rebuild) and only rebuilds on an actual q change.
   !>          'little_group_common' must never rebuild a mesh a sweep already
   !>          proved valid for its whole q-list (WP8): if one is cached under
   !>          that policy for the current mesh dims, it is left untouched;
   !>          only a genuinely missing mesh (e.g. this routine invoked outside
   !>          a sweep) falls back to the single current q_ss's little group.
   module subroutine force_full_bz_for_nonzero_q_gbt(this, context)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in) :: context

      if (.not. this%has_nonzero_q_gbt()) return

      this%use_symmetry_reduction = .false.
      this%use_time_reversal = .false.

      ! A band path is not a BZ integration mesh.  Do not replace it here;
      ! the finite-q mesh callers below are forced through generate_mp_mesh.
      if (allocated(this%k_path)) return

      select case (trim(this%q_symmetry_policy))
      case ('little_group')
         call this%ensure_kpoint_mesh(this%nk_mesh, sum(abs(this%k_offset)) > 1.0e-12_rp)
         return
      case ('little_group_common')
         if (this%mesh_cache_valid .and. trim(this%mesh_cache_policy) == 'little_group_common' .and. &
             allocated(this%k_points) .and. all(this%mesh_cache_dims == this%nk_mesh)) then
            return   ! sweep already built (and owns) the common-subgroup mesh
         end if
         call root_info(trim(context)//': little_group_common has no pre-built sweep mesh; '// &
                        'falling back to the current q_ss little group.', __FILE__, __LINE__)
         call this%generate_little_group_kpoint_mesh(this%nk_mesh, sum(abs(this%k_offset)) > 1.0e-12_rp)
         return
      end select

      call root_info(trim(context)//': nonzero-q GBT rebuilding the full chemical BZ mesh.', __FILE__, __LINE__)
      call this%generate_mp_mesh()
   end subroutine force_full_bz_for_nonzero_q_gbt

   !> @brief Set Monkhorst-Pack mesh dimensions.
   !> @param[inout] this Reciprocal object whose nk_mesh is updated.
   !> @param[in] nk1 Number of k-points along reciprocal axis 1.
   !> @param[in] nk2 Number of k-points along reciprocal axis 2.
   !> @param[in] nk3 Number of k-points along reciprocal axis 3.
   module subroutine set_kpoint_mesh(this, nk1, nk2, nk3)
      class(reciprocal), intent(inout) :: this
      integer, intent(in) :: nk1, nk2, nk3

      this%nk_mesh = [nk1, nk2, nk3]
      call g_logger%info('reciprocal%set_kpoint_mesh: Set k-point mesh to ' // &
         trim(int2str(nk1)) // 'x' // trim(int2str(nk2)) // 'x' // trim(int2str(nk3)), __FILE__, __LINE__)
   end subroutine set_kpoint_mesh

   !> @brief Generate reciprocal lattice vectors from the real-space lattice.
   !> @details Computes the reciprocal basis and reciprocal-cell volume from
   !>          lattice%a and lattice%alat for later k-point and Fourier work.
   !> @param[inout] this Reciprocal object receiving reciprocal_vectors and volume.
   module subroutine generate_reciprocal_vectors(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      real(rp), dimension(3, 3) :: real_vectors
      real(rp) :: det
      integer :: i

      ! Get real-space lattice vectors from lattice%a
      real_vectors = this%lattice%a

      ! Calculate determinant (volume of unit cell)
      det = real_vectors(1, 1) * (real_vectors(2, 2) * real_vectors(3, 3) - real_vectors(2, 3) * real_vectors(3, 2)) &
          - real_vectors(1, 2) * (real_vectors(2, 1) * real_vectors(3, 3) - real_vectors(2, 3) * real_vectors(3, 1)) &
          + real_vectors(1, 3) * (real_vectors(2, 1) * real_vectors(3, 2) - real_vectors(2, 2) * real_vectors(3, 1))

      ! Calculate reciprocal lattice vectors: b_i = 2π * (a_j × a_k) / V
      ! b1 = 2π * (a2 × a3) / V
      this%reciprocal_vectors(1, 1) = two_pi * (real_vectors(2, 2) * real_vectors(3, 3) - real_vectors(2, 3) * real_vectors(3, 2)) / det
      this%reciprocal_vectors(2, 1) = two_pi * (real_vectors(2, 3) * real_vectors(3, 1) - real_vectors(2, 1) * real_vectors(3, 3)) / det
      this%reciprocal_vectors(3, 1) = two_pi * (real_vectors(2, 1) * real_vectors(3, 2) - real_vectors(2, 2) * real_vectors(3, 1)) / det

      ! b2 = 2π * (a3 × a1) / V
      this%reciprocal_vectors(1, 2) = two_pi * (real_vectors(3, 2) * real_vectors(1, 3) - real_vectors(3, 3) * real_vectors(1, 2)) / det
      this%reciprocal_vectors(2, 2) = two_pi * (real_vectors(3, 3) * real_vectors(1, 1) - real_vectors(3, 1) * real_vectors(1, 3)) / det
      this%reciprocal_vectors(3, 2) = two_pi * (real_vectors(3, 1) * real_vectors(1, 2) - real_vectors(3, 2) * real_vectors(1, 1)) / det

      ! b3 = 2π * (a1 × a2) / V
      this%reciprocal_vectors(1, 3) = two_pi * (real_vectors(1, 2) * real_vectors(2, 3) - real_vectors(1, 3) * real_vectors(2, 2)) / det
      this%reciprocal_vectors(2, 3) = two_pi * (real_vectors(1, 3) * real_vectors(2, 1) - real_vectors(1, 1) * real_vectors(2, 3)) / det
      this%reciprocal_vectors(3, 3) = two_pi * (real_vectors(1, 1) * real_vectors(2, 2) - real_vectors(1, 2) * real_vectors(2, 1)) / det

      this%reciprocal_volume = (two_pi)**3 / abs(det)

      call root_info('reciprocal%generate_reciprocal_vectors: Reciprocal lattice vectors generated', __FILE__, __LINE__)
      
      ! Debug output - but use info level for now since debug is not enabled
      call root_info('reciprocal%generate_reciprocal_vectors: Real cell volume = ' // real2str(det), __FILE__, __LINE__)
      call root_info('reciprocal%generate_reciprocal_vectors: Reciprocal b1 = [' // &
         real2str(this%reciprocal_vectors(1, 1)) // ', ' // real2str(this%reciprocal_vectors(2, 1)) // ', ' // &
         real2str(this%reciprocal_vectors(3, 1)) // ']', __FILE__, __LINE__)
      call root_info('reciprocal%generate_reciprocal_vectors: Reciprocal b2 = [' // &
         real2str(this%reciprocal_vectors(1, 2)) // ', ' // real2str(this%reciprocal_vectors(2, 2)) // ', ' // &
         real2str(this%reciprocal_vectors(3, 2)) // ']', __FILE__, __LINE__)
      call root_info('reciprocal%generate_reciprocal_vectors: Reciprocal b3 = [' // &
         real2str(this%reciprocal_vectors(1, 3)) // ', ' // real2str(this%reciprocal_vectors(2, 3)) // ', ' // &
         real2str(this%reciprocal_vectors(3, 3)) // ']', __FILE__, __LINE__)
   end subroutine generate_reciprocal_vectors

   !> @brief Generate the configured Monkhorst-Pack k-point mesh.
   !> @details Builds full or symmetry-reduced fractional k-points, weights, and
   !>          optional MPI ownership metadata according to reciprocal settings.
   !> @param[inout] this Reciprocal object receiving k_points/k_weights state.
   module subroutine generate_mp_mesh(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      integer :: ik, ix, iy, iz, nk_irred
      real(rp) :: kx, ky, kz
      real(rp), parameter :: tol = 1.0e-8_rp
      integer :: shift(3), nfull
      real(rp), allocatable :: kpoints_full(:,:)

      this%nk_total = this%nk_mesh(1) * this%nk_mesh(2) * this%nk_mesh(3)

      ! Allocate k-point arrays
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('reciprocal.k_points', this%k_points, [3, this%nk_total])
      call g_safe_alloc%allocate('reciprocal.k_weights', this%k_weights, [this%nk_total])
#else
       if (allocated(this%k_points)) deallocate(this%k_points)
       if (allocated(this%k_weights)) deallocate(this%k_weights)
      allocate(this%k_points(3, this%nk_total))
      allocate(this%k_weights(this%nk_total))
#endif
      if (allocated(this%full_to_irred_k)) deallocate(this%full_to_irred_k)
      if (allocated(this%irred_to_full_k)) deallocate(this%irred_to_full_k)
      allocate(this%full_to_irred_k(this%nk_total))
      allocate(this%irred_to_full_k(this%nk_total))

      shift = [0, 0, 0]
      if (abs(this%k_offset(1)) > tol) shift(1) = 1
      if (abs(this%k_offset(2)) > tol) shift(2) = 1
      if (abs(this%k_offset(3)) > tol) shift(3) = 1

#ifdef USE_SPGLIB
      if (this%symmetry_analysis%spglib%is_available()) then
         nfull = this%symmetry_analysis%spglib%get_full_kpoint_mesh_with_points(this%nk_mesh, shift, kpoints_full)
         if (nfull == this%nk_total) then
            this%k_points = kpoints_full
            this%k_weights = 1.0_rp / real(this%nk_total, rp)
            do ik = 1, this%nk_total
               this%full_to_irred_k(ik) = ik
               this%irred_to_full_k(ik) = ik
            end do
            deallocate(kpoints_full)
            call setup_k_mesh_distribution(this, this%nk_total, .false.)
            call g_logger%info('reciprocal%generate_mp_mesh: Generated full mesh via spglib grid convention (' // &
                               trim(int2str(this%nk_total)) // ' k-points)', __FILE__, __LINE__)
            return
         end if
         if (allocated(kpoints_full)) deallocate(kpoints_full)
      end if
#endif

      ! Fallback: internal MP mesh formula
      ik = 0
      do iz = 1, this%nk_mesh(3)
         do iy = 1, this%nk_mesh(2)
            do ix = 1, this%nk_mesh(1)
               ik = ik + 1
               kx = (2.0_rp * ix - this%nk_mesh(1) - 1.0_rp) / (2.0_rp * this%nk_mesh(1)) + this%k_offset(1)
               ky = (2.0_rp * iy - this%nk_mesh(2) - 1.0_rp) / (2.0_rp * this%nk_mesh(2)) + this%k_offset(2)
               kz = (2.0_rp * iz - this%nk_mesh(3) - 1.0_rp) / (2.0_rp * this%nk_mesh(3)) + this%k_offset(3)
               this%k_points(1, ik) = kx
               this%k_points(2, ik) = ky
               this%k_points(3, ik) = kz
               this%k_weights(ik) = 1.0_rp / real(this%nk_total, rp)
               this%full_to_irred_k(ik) = ik
               this%irred_to_full_k(ik) = ik
            end do
         end do
      end do

      call setup_k_mesh_distribution(this, this%nk_total, .false.)
      call root_info('reciprocal%generate_mp_mesh: Generated Monkhorst-Pack mesh with ' // trim(int2str(this%nk_total)) // ' k-points', __FILE__, __LINE__)
   end subroutine generate_mp_mesh

   !> @brief Determine per-type basis sizes for reciprocal matrices.
   !> @details Maps each lattice atom type to sp/spd/spdf orbital counts used
   !>          when packing site-major k-space Hamiltonian and projection blocks.
   !> @param[inout] this Reciprocal object whose basis_size/max_orbs are updated.
   module subroutine set_basis_sizes(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      integer :: ntype
      character(len=10) :: basis_type

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('reciprocal.basis_size', this%basis_size, [this%lattice%ntype])
#else
    if (allocated(this%basis_size)) deallocate(this%basis_size)
      allocate(this%basis_size(this%lattice%ntype))
#endif

      ! Determine basis size for each atom type from active global basis.
      do ntype = 1, this%lattice%ntype
         this%basis_size(ntype) = norb
      end do

      this%max_orbs = nb

      call root_info('reciprocal%set_basis_sizes: Basis sizes set: max_orb_channels = ' // trim(int2str(this%max_orbs)), __FILE__, __LINE__)
   end subroutine set_basis_sizes

   !> @brief Build real-space neighbor vectors used by Fourier transforms.
   !> @details Converts lattice neighbor geometry into Cartesian and direct
   !>          per-type R-vector tables aligned with Hamiltonian neighbor slots.
   !> @param[inout] this Reciprocal object receiving ham_vec_type arrays.
   module subroutine build_neighbor_vectors(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      integer :: ntype, ia, nr, nn_max_loc, kk
      real(rp) :: r2
      real(rp), dimension(3, this%lattice%kk) :: cralat

      if (allocated(this%ham_vec_type) .and. allocated(this%ham_vec_type_direct)) then
         if (size(this%ham_vec_type, 1) == 3 .and. &
             size(this%ham_vec_type, 2) == this%lattice%nn_max .and. &
             size(this%ham_vec_type, 3) == this%lattice%ntype) then
            return
         end if
      end if

      call root_info('reciprocal%build_neighbor_vectors: Building neighbor vectors for each atom type', __FILE__, __LINE__)

      r2 = this%lattice%r2
      kk = this%lattice%kk
      cralat(1:3, 1:kk) = this%lattice%cr(1:3, 1:kk) * this%lattice%alat

      ! Allocate storage for all atom types
      ! Use nn_max as maximum neighbors
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('reciprocal.ham_vec_type', this%ham_vec_type, &
                                 [3, this%lattice%nn_max, this%lattice%ntype])
      call g_safe_alloc%allocate('reciprocal.ham_vec_type_direct', this%ham_vec_type_direct, &
                                 [3, this%lattice%nn_max, this%lattice%ntype])
#else
      if (allocated(this%ham_vec_type)) deallocate(this%ham_vec_type)
      allocate(this%ham_vec_type(3, this%lattice%nn_max, this%lattice%ntype))
      if (allocated(this%ham_vec_type_direct)) deallocate(this%ham_vec_type_direct)
      allocate(this%ham_vec_type_direct(3, this%lattice%nn_max, this%lattice%ntype))
#endif

      this%ham_vec_type = 0.0_rp
      this%ham_vec_type_direct = 0.0_rp

      ! Build neighbor vectors for each atom type
      do ntype = 1, this%lattice%ntype
         ia = this%lattice%atlist(ntype)
         nr = this%lattice%nn(ia, 1)
         nn_max_loc = nr

         ! Use clusba directly - no need for lattice%sbarvec storage here
         ! We're building type-specific ham_vec_type arrays instead
         call this%lattice%clusba(r2, cralat, ia, kk, kk, nn_max_loc, &
                                  this%ham_vec_type(:, 1:nr, ntype))

         ! Update nr with actual number of neighbors found
         nr = nn_max_loc

         ! Convert to fractional coordinates for ham_vec_type_direct
         ! Note: ham_vec_type from clusba is already in absolute Cartesian units (cralat = cr * alat)
         do nn_max_loc = 1, nr
            if (this%lattice%a_cart_inv_ready) then
               ! a_cart_inv expects input in units of alat, so divide by alat
               this%ham_vec_type_direct(:, nn_max_loc, ntype) = &
                  matmul(this%lattice%a_cart_inv, this%ham_vec_type(:, nn_max_loc, ntype) ) !/ this%lattice%alat
            else
               ! Fallback: use inverse of lattice vectors directly.
               this%ham_vec_type_direct(:, nn_max_loc, ntype) = &
                  matmul(inverse_3x3(this%lattice%a), this%ham_vec_type(:, nn_max_loc, ntype) / this%lattice%alat)
            end if
         end do

         call root_info('reciprocal%build_neighbor_vectors: Built ' // trim(int2str(nr)) // &
                        ' neighbor vectors for atom type ' // trim(int2str(ntype)), __FILE__, __LINE__)
         
      end do

      call root_info('reciprocal%build_neighbor_vectors: Completed neighbor vector build for all types', __FILE__, __LINE__)

   end subroutine build_neighbor_vectors

   !> @brief Configure MPI ownership for a k-point set.
   !> @param[inout] this Reciprocal object receiving local/global k maps.
   !> @param[in] nk_global Number of global k-points.
   !> @param[in] enable_distribution If true, distribute k-points over MPI ranks.
   module subroutine setup_k_mesh_distribution(this, nk_global, enable_distribution)
      class(reciprocal), intent(inout) :: this
      integer, intent(in) :: nk_global
      logical, intent(in) :: enable_distribution
      integer :: local_count
      integer :: ik
      logical :: have_k_points, have_k_path
      real(rp), allocatable :: unit_weights(:)

      ! Fortran does not require logical expressions to short-circuit.  Check
      ! allocation before querying a compatibility array's shape: a band path
      ! intentionally has no k_points array.
      have_k_points = .false.
      if (allocated(this%k_points)) have_k_points = size(this%k_points, 2) == nk_global
      have_k_path = .false.
      if (allocated(this%k_path)) have_k_path = size(this%k_path, 2) == nk_global

      if (have_k_points) then
         if (allocated(this%k_weights)) then
            if (size(this%k_weights) /= nk_global) then
               call g_logger%fatal('setup_k_mesh_distribution: k-point weight shape is invalid.', __FILE__, __LINE__)
            end if
            this%k_workset = make_kpoint_workset(this%k_points, this%k_weights, g_parallel_context, enable_distribution)
         else
         ! Hamiltonian-only mesh callers historically did not need integration
         ! weights.  Give their ownership object neutral unit weights without
         ! materializing a second compatibility array.
            allocate(unit_weights(nk_global)); unit_weights = 1.0_rp
            this%k_workset = make_kpoint_workset(this%k_points, unit_weights, g_parallel_context, enable_distribution)
            deallocate(unit_weights)
         end if
      else if (have_k_path) then
         ! A symmetry-generated band path is not a mesh and therefore has no
         ! k_points compatibility view.  It still needs an authoritative
         ! replicated workset for the shared reciprocal assembly path.
         allocate(unit_weights(nk_global)); unit_weights = 1.0_rp
         this%k_workset = make_kpoint_workset(this%k_path, unit_weights, g_parallel_context, .false.)
         deallocate(unit_weights)
      else
         call g_logger%fatal('setup_k_mesh_distribution: no complete k-point mesh or path is available.', __FILE__, __LINE__)
      end if

      if (allocated(this%k_l2g_map)) deallocate(this%k_l2g_map)
      if (allocated(this%k_g2l_map)) deallocate(this%k_g2l_map)

      if (enable_distribution) then
         call get_mpi_range(rank, nk_global, this%k_start, this%k_end, local_count, this%k_l2g_map, this%k_g2l_map, 'k')
         this%nk_local = local_count
         this%k_mesh_distributed_active = .true.
      else
         this%k_start = 1
         this%k_end = nk_global
         this%nk_local = nk_global
         this%k_mesh_distributed_active = .false.
         allocate(this%k_l2g_map(this%nk_local))
         allocate(this%k_g2l_map(nk_global))
         do ik = 1, nk_global
            this%k_l2g_map(ik) = ik
            this%k_g2l_map(ik) = ik
         end do
      end if
      ! Transitional mappings mirror the authoritative workset in one direction
      ! while legacy matrix/DOS storage is migrated.  The workset is validated
      ! at this public boundary and is never reconstructed from these maps.
      if (nk_global /= this%k_workset%nk_global .or. this%nk_local /= this%k_workset%nk_local .or. &
          this%k_start /= this%k_workset%global_start .or. this%k_end /= this%k_workset%global_end .or. &
          this%k_mesh_distributed_active .neqv. this%k_workset%distributed .or. &
          any(this%k_l2g_map /= this%k_workset%local_to_global) .or. any(this%k_g2l_map /= this%k_workset%global_to_local)) then
         call g_logger%fatal('setup_k_mesh_distribution: transitional k ownership view disagrees with workset.', __FILE__, __LINE__)
      end if
   end subroutine setup_k_mesh_distribution

   !> @brief Convert a local k-point index to a global k-point index.
   !> @param[in] this Reciprocal object containing distribution maps.
   !> @param[in] ik_local Local k-point index.
   !> @return Global k-point index.
   module integer function local_k_index_to_global(this, ik_local) result(ik_global)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: ik_local

      ik_global = this%k_workset%global_index(ik_local)
   end function local_k_index_to_global

   module subroutine require_replicated_k_workset(this, consumer)
      class(reciprocal), intent(in) :: this
      character(len=*), intent(in) :: consumer
      if (this%k_workset%distributed) then
         call g_logger%fatal(trim(consumer)//': requires a replicated kpoint_workset; distributed ownership is not supported here.', &
                             __FILE__, __LINE__)
      end if
   end subroutine require_replicated_k_workset

end submodule reciprocal_lifecycle
