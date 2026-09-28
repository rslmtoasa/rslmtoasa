!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Clean post-SCF production adapter for the rebuilt reciprocal TD-DFT
!> response services.
!>
!> This module owns input validation, accepted-state handoff, service request
!> assembly, and provenance output.  It deliberately owns no response
!> equation, radial augmentation formula, XC-kernel formula, Goldstone algebra,
!> or mode-extraction logic.
!------------------------------------------------------------------------------
submodule (linear_response_mod) linear_response_run

   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use mpi_mod, only: rank, numprocs
   use string_mod, only: lower, int2str
   use control_mod, only: control
   use lattice_mod, only: lattice
   use hamiltonian_mod, only: hamiltonian
   use energy_mod, only: energy
   use reciprocal_mod, only: reciprocal
   use recursion_mod, only: recursion
   use green_mod, only: green
   use symbolic_atom_mod, only: symbolic_atom
   use radial_ground_state_mod, only: radial_ground_state, RADIAL_PI, radial_simpson_weight
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use pauli_ground_state_projection_mod, only: compute_accepted_pauli_magnetization
   use logger_mod, only: g_logger
   implicit none
#ifdef VERSION
   character(len=*), parameter :: tddft_build_version = VERSION
#else
   character(len=*), parameter :: tddft_build_version = 'unavailable'
#endif

   ! The orchestration workers retain their validated internal routing shape;
   ! lr_run constructs it from the orthogonal &linear_response axes.  It is
   ! deliberately private so the retired TD-DFT configuration is not an API.
   character(len=*), parameter :: tddft_driver_backend_lehmann = 'lehmann'
   character(len=*), parameter :: tddft_driver_backend_native_rsgf = 'native_rsgf'
   character(len=*), parameter :: tddft_driver_backend_product_lehmann = 'product_lehmann'
   character(len=*), parameter :: tddft_driver_backend_projected_chi0 = 'projected_chi0'
   character(len=*), parameter :: tddft_driver_backend_compact_dyson = 'compact_dyson'
   character(len=*), parameter :: tddft_driver_backend_projected_mills = 'projected_mills'
   character(len=*), parameter :: tddft_driver_backend_projected_juelich = 'projected_juelich'
   character(len=*), parameter :: tddft_driver_route_direct_alsda = 'direct_alsda'
   character(len=*), parameter :: tddft_driver_route_goldstone_sumrule = 'goldstone_sumrule'

   type :: tddft_runtime_config
      logical :: enabled = .false.
      integer :: nq = 1
      integer :: nfrequency = 1
      character(len=32) :: channel = 'chi_plus'
      real(rp), allocatable :: q_list(:, :)
      real(rp), allocatable :: frequencies(:)
      real(rp), allocatable :: eta_values(:)
      real(rp) :: eta = 0.01_rp
      integer :: response_lmax = -1
      character(len=48) :: interaction_route = tddft_driver_route_direct_alsda
      logical :: goldstone_correction = .false.
      character(len=32) :: backend = tddft_driver_backend_lehmann
      character(len=8) :: projected_selector = 'spd'
      character(len=32) :: native_rsgf_provider = 'auto'
      integer :: gf_integration_points = 2001
      real(rp) :: gf_integration_eta = 0.0_rp
      real(rp) :: gf_energy_margin = 1.0_rp
      logical :: dyson_static_audit = .false.
      logical :: validate_interacting_covariance = .false.
      logical :: write_full_matrix = .true.
      character(len=256) :: output_file = 'tddft_response.dat'
   end type tddft_runtime_config

contains

   module function tddft_capability_is_supported(state, reason) result(ok)
      type(tddft_capability_state), intent(in) :: state
      character(len=*), intent(out) :: reason
      logical :: ok

      ok = .false.
      reason = 'unknown unsupported feature'
      if (.not. state%bulk_cell) then
         reason = 'non-bulk calculation geometry'
      else if (.not. state%collinear) then
         reason = 'noncollinear spin mode'
      else if (state%has_soc) then
         reason = 'spin-orbit coupling'
      else if (state%generalized_overlap .or. trim(state%reciprocal_mode) /= 'ham_only' .or. &
               .not. state%orthogonal) then
         reason = 'generalized overlap / non-orthogonal representation'
      else if (trim(state%hamiltonian_order) /= 'second') then
         reason = 'first-order Hamiltonian; validated TDDFT requires second order'
      else if (state%has_extra_operator) then
         reason = 'unsupported additive Hamiltonian operator'
      else if (state%basis_lmax /= 1 .and. state%basis_lmax /= 2) then
         reason = 'basis outside supported sp/spd (lmax=1 or 2)'
      else
         ok = .true.
         reason = 'supported collinear ham_only sp/spd baseline'
      end if
   end function tddft_capability_is_supported

   module subroutine require_tddft_capability(state)
      type(tddft_capability_state), intent(in) :: state
      character(len=256) :: reason
      if (.not. tddft_capability_is_supported(state, reason)) then
         error stop 'TDDFT capability gate: unsupported feature: '//trim(reason)
      end if
   end subroutine require_tddft_capability

   subroutine derive_runtime_config(config, runtime)
      type(linear_response_config), intent(in) :: config
      type(tddft_runtime_config), intent(out) :: runtime
      integer :: n_q, n_frequency, n_eta

      runtime%enabled = .true.
      n_q = size(config%q_list, 2)
      n_frequency = size(config%omega_grid)
      n_eta = size(config%eta_grid)
      runtime%nq = n_q
      runtime%nfrequency = n_frequency
      runtime%channel = trim(config%channel)
      runtime%eta = config%eta
      runtime%response_lmax = config%response_lmax
      runtime%projected_selector = trim(config%projection)
      runtime%native_rsgf_provider = trim(config%native_rsgf_provider)
      if (trim(config%realspace_solver) == 'block') runtime%native_rsgf_provider = 'block'
      if (trim(config%realspace_solver) == 'chebyshev') runtime%native_rsgf_provider = 'chebyshev'
      runtime%gf_integration_points = config%gf_integration_points
      runtime%gf_integration_eta = config%gf_integration_eta
      runtime%gf_energy_margin = config%gf_energy_margin
      runtime%write_full_matrix = config%write_full_matrix
      runtime%output_file = trim(config%output_file)
      runtime%dyson_static_audit = trim(config%diagnostics) == 'invariants'
      runtime%validate_interacting_covariance = trim(config%diagnostics) == 'invariants'
      runtime%goldstone_correction = .false.
      allocate(runtime%q_list(3, n_q), runtime%frequencies(n_frequency), runtime%eta_values(n_eta))
      runtime%q_list = config%q_list
      runtime%frequencies = config%omega_grid
      runtime%eta_values = config%eta_grid

      select case (trim(config%formulation))
      case ('tddft')
         select case (trim(config%representation)//':'//trim(config%bare_response)//':'//trim(config%interaction))
         case ('radial_points:lehmann:alsda')
            runtime%backend = tddft_driver_backend_lehmann
            runtime%interaction_route = tddft_driver_route_direct_alsda
         case ('radial_points:lehmann:lcmm')
            runtime%backend = tddft_driver_backend_lehmann
            runtime%interaction_route = tddft_driver_route_goldstone_sumrule
         case ('radial_points:realspace_gf:alsda')
            runtime%backend = tddft_driver_backend_native_rsgf
            runtime%interaction_route = tddft_driver_route_direct_alsda
         case ('radial_points:realspace_gf:lcmm')
            runtime%backend = tddft_driver_backend_native_rsgf
            runtime%interaction_route = tddft_driver_route_goldstone_sumrule
         case ('product_compact:lehmann:alsda')
            runtime%backend = tddft_driver_backend_compact_dyson
            runtime%interaction_route = tddft_driver_route_direct_alsda
         case ('product_compact:lehmann:none')
            runtime%backend = tddft_driver_backend_product_lehmann
            runtime%interaction_route = tddft_driver_route_direct_alsda
         end select
      case ('projected')
         runtime%interaction_route = tddft_driver_route_direct_alsda
         select case (trim(config%interaction))
         case ('none')
            runtime%backend = tddft_driver_backend_projected_chi0
         case ('stoner_fit')
            runtime%backend = tddft_driver_backend_projected_mills
         case ('lcmm')
            runtime%backend = tddft_driver_backend_projected_juelich
         end select
      end select
      call validate_tddft_runtime_config(runtime)
   end subroutine derive_runtime_config

   subroutine validate_tddft_runtime_config(config)
      type(tddft_runtime_config), intent(in) :: config
      integer :: eta_selected_index, eta_holdout_index

      if (.not. config%enabled) return
      if (trim(config%channel) /= 'chi_plus' .and. trim(config%channel) /= 'chi_minus') then
         error stop 'TDDFT input: unsupported response channel; use chi_plus or chi_minus'
      end if
      if (config%eta <= 0.0_rp .or. .not. allocated(config%eta_values) .or. any(config%eta_values <= 0.0_rp)) then
         error stop 'TDDFT input: every physical response eta must be positive'
      end if
      if (config%response_lmax < -1 .or. config%response_lmax > 4) then
         error stop 'TDDFT input: response_lmax is outside the validated sp/spd product baseline'
      end if
      if (trim(config%backend) == tddft_driver_backend_projected_mills .or. &
          trim(config%backend) == tddft_driver_backend_projected_juelich) then
         if (trim(config%projected_selector) /= 'd' .and. trim(config%projected_selector) /= 'spd' .and. &
             trim(config%projected_selector) /= 'both') error stop 'TDDFT input: projected selector is invalid'
         if (config%response_lmax >= 0 .and. config%response_lmax /= 4) then
            error stop 'TDDFT input: projected interaction requires response_lmax=4'
         end if
         if (trim(config%backend) == tddft_driver_backend_projected_juelich) then
            call select_projected_juelich_eta_indices(config%eta_values, eta_selected_index, eta_holdout_index)
            if (.not. any(sum(abs(config%q_list), dim=1) <= 1.0e-12_rp)) &
               error stop 'TDDFT input: projected_juelich requires Gamma'
         end if
      else if (trim(config%backend) /= tddft_driver_backend_projected_mills .and. &
               trim(config%backend) /= tddft_driver_backend_projected_juelich .and. size(config%eta_values) /= 1) then
         error stop 'TDDFT input: n_eta greater than one is restricted to projected validation backends'
      end if
      if (trim(config%backend) == tddft_driver_backend_native_rsgf) then
         if (config%gf_integration_points < 3 .or. mod(config%gf_integration_points,2) == 0) then
            error stop 'TDDFT input: realspace_gf requires an odd gf_integration_points value >= 3'
         end if
         if (config%gf_energy_margin <= 0.0_rp) error stop 'TDDFT input: realspace_gf requires gf_energy_margin positive'
         if (trim(config%native_rsgf_provider) /= 'auto' .and. trim(config%native_rsgf_provider) /= 'block' .and. &
             trim(config%native_rsgf_provider) /= 'block_recursion' .and. trim(config%native_rsgf_provider) /= 'chebyshev') then
            error stop 'TDDFT input: realspace_solver provider is invalid'
         end if
      end if
      if (trim(config%backend) == tddft_driver_backend_compact_dyson) then
         if (config%dyson_static_audit .and. find_gamma_q_index(config%q_list) == 0) then
            error stop 'TDDFT input: diagnostics=invariants requires Gamma for compact Dyson'
         end if
         if (config%validate_interacting_covariance) call validate_covariance_q_pair(config%q_list)
      end if
      if (trim(config%interaction_route) == tddft_driver_route_goldstone_sumrule .and. &
          .not. any(sum(abs(config%q_list), dim=1) <= 1.0e-12_rp)) then
         error stop 'TDDFT input: lcmm requires q=(0,0,0) in q_list for the static reference'
      end if
   end subroutine validate_tddft_runtime_config

   module subroutine lr_run(this, control_obj, lattice_obj, hamiltonian_obj, energy_obj, self_obj, reciprocal_obj, &
                            recursion_obj, green_obj, scf_converged)
      class(linear_response), intent(in) :: this
      type(control), intent(in) :: control_obj
      type(lattice), intent(inout) :: lattice_obj
      type(hamiltonian), intent(inout) :: hamiltonian_obj
      type(energy), intent(in) :: energy_obj
      type(self), intent(inout) :: self_obj
      type(reciprocal), intent(inout) :: reciprocal_obj
      type(recursion), target, intent(inout), optional :: recursion_obj
      type(green), target, intent(inout), optional :: green_obj
      logical, intent(in), optional :: scf_converged
      type(tddft_runtime_config) :: runtime
      logical :: converged

      converged = .true.
      if (present(scf_converged)) converged = scf_converged
      select case (trim(this%config%formulation))
      case ('rotation')
         call lr_run_rotation(this, control_obj, lattice_obj, hamiltonian_obj, energy_obj, self_obj, reciprocal_obj)
      case ('tddft', 'projected')
         if (.not. present(recursion_obj) .or. .not. present(green_obj)) then
            error stop "linear_response: tddft/projected formulations require recursion and green services"
         end if
         call derive_runtime_config(this%config, runtime)
         call run_tddft_production(runtime, control_obj, lattice_obj, hamiltonian_obj, energy_obj, reciprocal_obj, &
            recursion_obj, green_obj, converged, self_obj%use_kspace)
      case default
         error stop 'linear_response: unsupported formulation'
      end select
   end subroutine lr_run

   subroutine validate_tddft_production_capability(config, control_obj, lattice_obj, hamiltonian_obj, reciprocal_obj)
      type(tddft_runtime_config), intent(in) :: config
      type(control), intent(in) :: control_obj
      type(lattice), intent(in) :: lattice_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj
      type(reciprocal), intent(in) :: reciprocal_obj
      type(tddft_capability_state) :: state
      character(len=256) :: reason

      if (.not. config%enabled) return
      state%bulk_cell = control_obj%calctype == 'B'
      state%collinear = control_obj%is_collinear()
      ! The lsham allocation is ordinary Hamiltonian storage and may exist as
      ! a zero-filled cache even when SOC is disabled; nsp is the capability
      ! authority for this production baseline.
      state%has_soc = control_obj%has_soc()
      state%reciprocal_mode = trim(reciprocal_obj%reciprocal_mode)
      state%generalized_overlap = trim(state%reciprocal_mode) /= 'ham_only' .or. allocated(reciprocal_obj%sk_overlap)
      state%orthogonal = .not. state%generalized_overlap
      state%hamiltonian_order = trim(reciprocal_obj%kspace_ham_order)
      if (trim(state%hamiltonian_order) == 'auto') then
         if (hamiltonian_obj%hoh) then
            state%hamiltonian_order = 'second'
         else
            state%hamiltonian_order = 'first'
         end if
      end if
      state%basis_lmax = -1
      if (allocated(lattice_obj%symbolic_atoms) .and. size(lattice_obj%symbolic_atoms) > 0) then
         state%basis_lmax = lattice_obj%symbolic_atoms(1)%potential%lmax
      end if
      state%has_extra_operator = hamiltonian_obj%local_axis .or. hamiltonian_obj%orb_pol .or. &
         hamiltonian_obj%ccor_2c .or. hamiltonian_obj%hubbard_u_general_check .or. &
         hamiltonian_obj%hubbard_u_impurity_check .or. hamiltonian_obj%hubbard_u_sc_check .or. &
         hamiltonian_obj%hubbard_v_check .or. control_obj%constraints_enable
      if (.not. tddft_capability_is_supported(state, reason)) then
         error stop 'TDDFT capability gate: unsupported feature: '//trim(reason)
      end if
      if (config%response_lmax > 2*state%basis_lmax) then
         error stop 'TDDFT capability gate: response_lmax exceeds the accepted angular-product cutoff'
      end if
      if (trim(config%backend) == tddft_driver_backend_native_rsgf) then
         if (lattice_obj%njij < 1) then
            error stop 'TDDFT capability gate: native_rsgf requires a complete lattice%ijpair set'
         end if
         if (numprocs /= 1) then
            error stop 'TDDFT capability gate: native_rsgf registered baseline requires a serial replicated pair workset'
         end if
         if (trim(config%native_rsgf_provider) == 'block' .or. trim(config%native_rsgf_provider) == 'block_recursion') then
            if (trim(lower(control_obj%recur)) /= 'block') then
               error stop 'TDDFT capability gate: native block provider requires control%recur=block'
            end if
         else if (trim(config%native_rsgf_provider) == 'chebyshev') then
            if (trim(lower(control_obj%recur)) /= 'chebyshev') then
               error stop 'TDDFT capability gate: native Chebyshev provider requires control%recur=chebyshev'
            end if
         else if (trim(lower(control_obj%recur)) /= 'block' .and. trim(lower(control_obj%recur)) /= 'chebyshev') then
            error stop 'TDDFT capability gate: native_rsgf requires control%recur=block or chebyshev'
         end if
      end if
   end subroutine validate_tddft_production_capability

   !> Fail-closed validation for the direct k-space-SCF -> TDDFT handoff.
   !> No response service is allowed to repair or reinterpret this state.
   subroutine validate_accepted_kspace_scf_handoff(reciprocal_obj, energy_obj)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(energy), intent(in) :: energy_obj
      integer :: expected_nk
      real(rp) :: scale, weight_sum

      if (.not. reciprocal_obj%auto_find_fermi) then
         error stop 'TDDFT k-space-SCF handoff: accepted state has fixed-EF semantics'
      end if
      if (.not. reciprocal_obj%canonical_energy_valid) then
         error stop 'TDDFT k-space-SCF handoff: accepted canonical occupations are unavailable'
      end if
      if (.not. allocated(reciprocal_obj%k_workset%points) .or. &
          .not. allocated(reciprocal_obj%k_workset%weights) .or. &
          .not. allocated(reciprocal_obj%eigenvalues) .or. &
          .not. allocated(reciprocal_obj%eigenvectors)) then
         error stop 'TDDFT k-space-SCF handoff: accepted mesh/eigensystem is incomplete'
      end if
      call reciprocal_obj%require_replicated_k_workset('TDDFT k-space-SCF handoff')
      expected_nk = product(reciprocal_obj%nk_mesh)
      if (expected_nk < 1 .or. reciprocal_obj%k_workset%nk_global /= expected_nk .or. &
          reciprocal_obj%k_workset%nk_local /= expected_nk .or. .not. reciprocal_obj%k_workset%complete_bz) then
         error stop 'TDDFT k-space-SCF handoff: accepted state is not the complete requested Monkhorst-Pack mesh'
      end if
      if (size(reciprocal_obj%k_workset%points, 2) /= expected_nk .or. &
          size(reciprocal_obj%k_workset%weights) /= expected_nk .or. &
          size(reciprocal_obj%eigenvalues, 2) /= expected_nk .or. &
          size(reciprocal_obj%eigenvectors, 3) /= expected_nk) then
         error stop 'TDDFT k-space-SCF handoff: accepted state arrays do not cover the complete mesh'
      end if
      if (allocated(reciprocal_obj%k_points)) then
         if (size(reciprocal_obj%k_points, 2) /= expected_nk) then
            error stop 'TDDFT k-space-SCF handoff: compatibility k-points differ from the authoritative workset'
         end if
         if (maxval(abs(reciprocal_obj%k_points - reciprocal_obj%k_workset%points)) > 2.0e-14_rp) then
            error stop 'TDDFT k-space-SCF handoff: compatibility k-points differ from the authoritative workset'
         end if
      end if
      weight_sum = sum(reciprocal_obj%k_workset%weights)
      if (abs(weight_sum - reciprocal_obj%canonical_weight_sum) > 2.0e-12_rp*max(1.0_rp, weight_sum)) then
         error stop 'TDDFT k-space-SCF handoff: canonical weight sum differs from accepted mesh weights'
      end if
      scale = max(1.0_rp, abs(energy_obj%fermi), abs(reciprocal_obj%fermi_level))
      if (abs(energy_obj%fermi - reciprocal_obj%fermi_level) > 2.0e-12_rp*scale) then
         error stop 'TDDFT k-space-SCF handoff: TDDFT energy EF differs from accepted SCF EF'
      end if
      if (reciprocal_obj%total_electrons > 0.0_rp .and. &
          abs(reciprocal_obj%canonical_electron_count - reciprocal_obj%total_electrons) > &
          1.0e-9_rp*max(1.0_rp, reciprocal_obj%total_electrons)) then
         error stop 'TDDFT k-space-SCF handoff: accepted reciprocal state violates electron conservation'
      end if
   end subroutine validate_accepted_kspace_scf_handoff

   !> Write the exact left-state data consumed by the response routines.  The
   !> occupation-weighted density matrix is a gauge-invariant eigenvector
   !> identity diagnostic and is intentionally serialized for the harness.
   subroutine write_tddft_state_artifact(filename, reciprocal_obj, left_state, direct_handoff)
      character(len=*), intent(in) :: filename
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lr_electronic_state), intent(in) :: left_state
      logical, intent(in) :: direct_handoff
      integer :: unit, ik, ib, i, j, nbasis, nbands, nk
      real(rp) :: weight_sum, k_fingerprint(5)
      real(rp) :: eigenvalue_checksum, occupation_checksum, eigenvector_unitarity_residual
      complex(rp) :: eigenvector_checksum
      real(rp), allocatable :: occupations(:)
      complex(rp) :: density_element
      character(len=64) :: state_source

      if (rank /= 0) return
      nbands = left_state%nbands
      nbasis = left_state%nbasis
      nk = left_state%nk
      if (.not. left_state%occupations_are_explicit .or. .not. allocated(left_state%occupations)) then
         error stop 'TDDFT state artifact: left-state occupations are not explicit'
      end if
      if (direct_handoff) then
         state_source = 'accepted_kspace_scf_cache'
      else
         state_source = 'diagnostic_frozen_post_scf_rebuild'
      end if
      weight_sum = sum(left_state%k_weights)
      call state_identity_checksums(left_state, eigenvalue_checksum, eigenvector_checksum, occupation_checksum, &
         eigenvector_unitarity_residual)
      k_fingerprint = 0.0_rp
      do ik = 1, nk
         k_fingerprint(1) = k_fingerprint(1) + left_state%k_weights(ik)
         k_fingerprint(2:4) = k_fingerprint(2:4) + left_state%k_weights(ik)*left_state%k_points(:, ik)
         k_fingerprint(5) = k_fingerprint(5) + real(ik, rp)*left_state%k_weights(ik)
      end do
      allocate(occupations(nbands))
      open(newunit=unit, file=trim(filename), status='replace', action='write')
      write(unit, '(a)') '# reciprocal state-consistency artifact'
      write(unit, '(a)') '# state_role = tddft_left_state'
      write(unit, '(a,a)') '# state_source = ', trim(state_source)
      write(unit, '(a,l1)') '# direct_accepted_state_handoff = ', direct_handoff
      write(unit, '(a,3(i0,1x))') '# requested_k_mesh = ', reciprocal_obj%nk_mesh
      write(unit, '(a,3(i0,1x))') '# actual_k_mesh = ', reciprocal_obj%nk_mesh
      write(unit, '(a,i0)') '# actual_k_count = ', nk
      write(unit, '(a,i0)') '# nbasis = ', nbasis
      write(unit, '(a,i0)') '# nbands = ', nbands
      write(unit, '(a,es24.16)') '# k_weight_sum = ', weight_sum
      write(unit, '(a,5(es24.16,1x))') '# k_fingerprint_checksums = ', k_fingerprint
      write(unit, '(a,es24.16)') '# eigenvalue_checksum_Ry = ', eigenvalue_checksum
      write(unit, '(a,2(es24.16,1x))') '# eigenvector_checksum_real_imag = ', real(eigenvector_checksum, rp), &
         aimag(eigenvector_checksum)
      write(unit, '(a,es24.16)') '# occupation_checksum = ', occupation_checksum
      write(unit, '(a,es24.16)') '# eigenvector_unitarity_max_abs = ', eigenvector_unitarity_residual
      write(unit, '(a,es24.16)') '# fermi_level_Ry = ', left_state%fermi_level
      write(unit, '(a,es24.16)') '# target_electron_count = ', reciprocal_obj%total_electrons
      write(unit, '(a,es24.16)') '# accepted_electron_count = ', reciprocal_obj%canonical_electron_count
      write(unit, '(a,es24.16)') '# temperature_K = ', left_state%temperature
      write(unit, '(a,a)') '# reciprocal_mode = ', trim(left_state%reciprocal_mode)
      write(unit, '(a,a)') '# kspace_hamiltonian_order = ', trim(left_state%hamiltonian_order)
      write(unit, '(a)') '# occupation_semantics = immutable lr_snapshot occupations from the accepted reciprocal EF and temperature'
      write(unit, '(a)') '# eigenvector_identity = occupation-weighted one-particle density-matrix projector residual'
      write(unit, '(a)') '# columns: tag k_index band_or_row column eigenvalue_or_occupation real imag'
      do ik = 1, nk
         write(unit, '(a,1x,i0,1x,4(es24.16,1x))') 'K', ik, left_state%k_points(:, ik), left_state%k_weights(ik)
         do ib = 1, nbands
            occupations(ib) = left_state%occupations(ib, ik)
            write(unit, '(a,1x,2(i0,1x),2(es24.16,1x))') 'O', ik, ib, left_state%eigenvalues(ib, ik), occupations(ib)
         end do
         do i = 1, nbasis
            do j = 1, nbasis
               density_element = sum(occupations(:)*left_state%eigenvectors(i, :, ik)* &
                                     conjg(left_state%eigenvectors(j, :, ik)))
               write(unit, '(a,1x,3(i0,1x),2(es24.16,1x))') 'D', ik, i, j, real(density_element, rp), &
                  aimag(density_element)
            end do
         end do
      end do
      close(unit)
      deallocate(occupations)
   end subroutine write_tddft_state_artifact

   !> Deterministic, machine-readable fingerprints for the immutable response
   !> state.  These are diagnostics rather than cryptographic hashes: the
   !> indexed weights make accidental state replacement visible while the
   !> unitary residual checks the eigenvector normalization independently.
   subroutine state_identity_checksums(state, eigenvalue_checksum, eigenvector_checksum, occupation_checksum, &
                                       eigenvector_unitarity_residual)
      type(lr_electronic_state), intent(in) :: state
      real(rp), intent(out) :: eigenvalue_checksum, occupation_checksum, eigenvector_unitarity_residual
      complex(rp), intent(out) :: eigenvector_checksum
      integer :: ik, ib, i, jb
      complex(rp) :: overlap, expected

      eigenvalue_checksum = 0.0_rp
      occupation_checksum = 0.0_rp
      eigenvector_checksum = cmplx(0.0_rp, 0.0_rp, rp)
      eigenvector_unitarity_residual = 0.0_rp
      do ik = 1, state%nk
         do ib = 1, state%nbands
            eigenvalue_checksum = eigenvalue_checksum + (real(ib, rp) + 0.001_rp*real(ik, rp))* &
               state%eigenvalues(ib, ik)
            occupation_checksum = occupation_checksum + (real(ib, rp) + 0.001_rp*real(ik, rp))* &
               state%occupations(ib, ik)
            do i = 1, state%nbasis
               eigenvector_checksum = eigenvector_checksum + &
                  cmplx(real(i, rp) + 0.01_rp*real(ib, rp) + 0.0001_rp*real(ik, rp), 0.0_rp, rp)* &
                  state%eigenvectors(i, ib, ik)
            end do
         end do
         do ib = 1, state%nbands
            do jb = 1, state%nbands
               overlap = sum(conjg(state%eigenvectors(:, ib, ik))*state%eigenvectors(:, jb, ik))
               if (ib == jb) then
                  expected = cmplx(1.0_rp, 0.0_rp, rp)
               else
                  expected = cmplx(0.0_rp, 0.0_rp, rp)
               end if
               eigenvector_unitarity_residual = max(eigenvector_unitarity_residual, abs(overlap - expected))
            end do
         end do
      end do
   end subroutine state_identity_checksums

   !> Deterministic checksum/norm for the accepted radial/LMTO snapshot used
   !> to build the DRESP-01 site projection.  The arrays are included rather
   !> than just their scalar moments so a radial-state replacement is visible.
   subroutine radial_identity_checksums(states, radial_checksum, radial_l2_norm)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(out) :: radial_checksum, radial_l2_norm
      integer :: isite, ir
      real(rp) :: coefficient

      radial_checksum = 0.0_rp
      radial_l2_norm = 0.0_rp
      do isite = 1, size(states)
         coefficient = real(isite, rp)
         radial_checksum = radial_checksum + coefficient*(states(isite)%a + states(isite)%b + states(isite)%rmax + &
            sum(states(isite)%r) + sum(states(isite)%rho_weighted_up) + sum(states(isite)%rho_weighted_down) + &
            sum(states(isite)%pauli_large) + sum(states(isite)%pauli_large_dot) + sum(states(isite)%pauli_small) + &
            sum(states(isite)%pauli_sr_correction) + sum(states(isite)%pauli_enu) + &
            states(isite)%integrated_n_up + states(isite)%integrated_n_down + states(isite)%integrated_spin_number + &
            states(isite)%integrated_moment_muB)
         do ir = 1, size(states(isite)%r)
            radial_checksum = radial_checksum + (coefficient + 0.001_rp*real(ir, rp))* &
               (states(isite)%r(ir) + states(isite)%rho_weighted_up(ir) + states(isite)%rho_weighted_down(ir))
         end do
         radial_l2_norm = radial_l2_norm + sum(states(isite)%r**2) + &
            sum(states(isite)%rho_weighted_up**2) + sum(states(isite)%rho_weighted_down**2) + &
            sum(states(isite)%pauli_large**2) + sum(states(isite)%pauli_large_dot**2) + &
            sum(states(isite)%pauli_small**2) + sum(states(isite)%pauli_sr_correction**2) + &
            sum(states(isite)%pauli_enu**2)
      end do
      radial_l2_norm = sqrt(max(0.0_rp, radial_l2_norm))
   end subroutine radial_identity_checksums

   !> Compare the immutable left snapshot against its reciprocal source before
   !> any response contraction.  The density-matrix residual is invariant under
   !> eigenvector phase choices and rotations within equal-occupation subspaces.
   subroutine state_consistency_metrics(reciprocal_obj, left_state, mesh_max, weight_max, ef_diff, eigen_max, occ_max, &
                                        projector_max, projector_frobenius)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lr_electronic_state), intent(in) :: left_state
      real(rp), intent(out) :: mesh_max, weight_max, ef_diff, eigen_max, occ_max, projector_max, projector_frobenius
      integer :: ik, ib, i, j, nbasis, nbands, nk
      real(rp) :: expected_occupation
      real(rp), allocatable :: occupations(:)
      complex(rp) :: source_density, snapshot_density, difference

      nk = left_state%nk
      nbands = left_state%nbands
      nbasis = left_state%nbasis
      mesh_max = maxval(abs(reciprocal_obj%k_workset%points - left_state%k_points))
      weight_max = maxval(abs(reciprocal_obj%k_workset%weights - left_state%k_weights))
      ef_diff = abs(reciprocal_obj%fermi_level - left_state%fermi_level)
      eigen_max = maxval(abs(reciprocal_obj%eigenvalues - left_state%eigenvalues))
      occ_max = 0.0_rp
      projector_max = 0.0_rp
      projector_frobenius = 0.0_rp
      allocate(occupations(nbands))
      do ik = 1, nk
         do ib = 1, nbands
            expected_occupation = lr_fermi_dirac_occupation(reciprocal_obj%eigenvalues(ib, ik), &
               reciprocal_obj%fermi_level, reciprocal_obj%temperature)
            occupations(ib) = expected_occupation
            occ_max = max(occ_max, abs(expected_occupation - left_state%occupations(ib, ik)))
         end do
         do i = 1, nbasis
            do j = 1, nbasis
               source_density = sum(occupations(:)*reciprocal_obj%eigenvectors(i, :, ik)* &
                                    conjg(reciprocal_obj%eigenvectors(j, :, ik)))
               snapshot_density = sum(left_state%occupations(:, ik)*left_state%eigenvectors(i, :, ik)* &
                                       conjg(left_state%eigenvectors(j, :, ik)))
               difference = source_density - snapshot_density
               projector_max = max(projector_max, abs(difference))
               projector_frobenius = projector_frobenius + abs(difference)**2
            end do
         end do
      end do
      projector_frobenius = sqrt(projector_frobenius)
      deallocate(occupations)
   end subroutine state_consistency_metrics

   !> Run the response after SCF has accepted its state. The input objects are
   !> read-only except for the dedicated reciprocal eigenpair service cache.
   !> When accepted_kspace_scf is true, reciprocal_obj is the live cache owned
   !> by self and is consumed without a mesh generation, Hamiltonian build, or
   !> diagonalization in this driver.
   subroutine run_tddft_production(config, control_obj, lattice_obj, hamiltonian_obj, energy_obj, reciprocal_obj, &
                                   recursion_obj, green_obj, scf_converged, accepted_kspace_scf)
      type(tddft_runtime_config), intent(in) :: config
      type(control), intent(in) :: control_obj
      type(lattice), target, intent(in) :: lattice_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj
      type(energy), intent(in) :: energy_obj
      type(reciprocal), intent(inout) :: reciprocal_obj
      type(recursion), target, intent(inout) :: recursion_obj
      type(green), target, intent(inout) :: green_obj
      logical, intent(in) :: scf_converged
      logical, intent(in), optional :: accepted_kspace_scf

      type(radial_ground_state), pointer :: ground_states(:)
      type(lmto_radial_basis), allocatable, target :: radial_bases(:)
      type(response_space_layout), target :: response_space
      type(lr_electronic_state), target :: left_state
      type(lr_electronic_state), allocatable, target :: endpoints(:)
      type(tddft_production_result) :: result
      type(tddft_native_rsgf_provider), target :: native_provider
      real(rp), allocatable, target :: native_site_positions(:, :)
      real(rp), allocatable :: accepted_pauli_magnetization(:, :)
      real(rp) :: ignore_real
      integer :: first, last, nsite, response_lmax, isite, iq
      logical :: use_accepted_kspace_scf, need_complete_sr

      use_accepted_kspace_scf = .false.
      if (present(accepted_kspace_scf)) use_accepted_kspace_scf = accepted_kspace_scf
      need_complete_sr = trim(config%backend) == tddft_driver_backend_product_lehmann

      if (.not. config%enabled) return
      call validate_tddft_production_capability(config, control_obj, lattice_obj, hamiltonian_obj, reciprocal_obj)
      if (.not. scf_converged) error stop 'TDDFT production driver: an accepted converged ground state is required'
      if (lattice_obj%nrec < 1) error stop 'TDDFT production driver: no accepted response sites are available'

      first = lattice_obj%nbulk + 1
      last = lattice_obj%nbulk + lattice_obj%nrec
      if (.not. allocated(lattice_obj%symbolic_atoms) .or. last > size(lattice_obj%symbolic_atoms)) then
         error stop 'TDDFT production driver: accepted symbolic-atom state is incomplete'
      end if
      ground_states => lattice_obj%symbolic_atoms(first:last)%radial_ground_state
      nsite = size(ground_states)
      do isite = 1, nsite
         if (.not. ground_states(isite)%valid .or. .not. ground_states(isite)%accepted .or. &
             .not. ground_states(isite)%pauli_basis_valid) then
            error stop 'TDDFT production driver: accepted LR-01 radial snapshot is incomplete'
         end if
      end do
      ! The response angular space is the product space of two LMTO
      ! endpoints, so its complete cutoff is twice the radial orbital cutoff.
      response_lmax = 2*ground_states(1)%pauli_lmax
      if (config%response_lmax >= 0) response_lmax = config%response_lmax
      if (response_lmax < 0 .or. response_lmax > 2*ground_states(1)%pauli_lmax) then
         error stop 'TDDFT production driver: response angular cutoff is outside the accepted radial basis'
      end if
      allocate(radial_bases(nsite))
      do isite = 1, nsite
         if (size(ground_states(isite)%r) /= size(ground_states(1)%r) .or. &
             any(abs(ground_states(isite)%r - ground_states(1)%r) > 1.0e-12_rp)) then
            error stop 'TDDFT production driver: response sites do not share the accepted radial mesh'
         end if
         if (ground_states(isite)%pauli_lmax /= ground_states(1)%pauli_lmax) then
            error stop 'TDDFT production driver: response sites do not share the accepted Pauli cutoff'
         end if
         call radial_basis_from_snapshot(ground_states(isite), radial_bases(isite), need_complete_sr)
         call radial_bases(isite)%require_supported(trim(reciprocal_obj%reciprocal_mode), &
            effective_hamiltonian_order(reciprocal_obj, hamiltonian_obj), .true., .true., &
            control_obj%has_soc(), .false.)
      end do
      call response_space%initialize(nsite, response_lmax, ground_states(1)%r, ground_states(1)%a, ground_states(1)%b, 1)

      if (use_accepted_kspace_scf) then
         call validate_accepted_kspace_scf_handoff(reciprocal_obj, energy_obj)
      else
         ! Diagnostic-only legacy path: construct a reciprocal response state
         ! from the accepted real-space potential at the externally accepted
         ! EF.  This remains available for historical comparisons, but is not
         ! the production k-space-SCF handoff.
         reciprocal_obj%use_symmetry_reduction = .false.
         reciprocal_obj%use_time_reversal = .false.
         reciprocal_obj%dos_method = 'tetrahedron'
         reciprocal_obj%auto_find_fermi = .false.
         reciprocal_obj%fermi_level = energy_obj%fermi
         reciprocal_obj%kspace_ham_order = effective_hamiltonian_order(reciprocal_obj, hamiltonian_obj)
         call reciprocal_obj%generate_mp_mesh()
         call reciprocal_obj%build_kspace_hamiltonian()
         call reciprocal_obj%diagonalize_hamiltonian()
         reciprocal_obj%fermi_level = energy_obj%fermi
         ignore_real = reciprocal_obj%calculate_canonical_band_energy(.false.)
      end if
      call reciprocal_obj%require_replicated_k_workset('TDDFT production driver')
      call lr_snapshot_from_reciprocal(reciprocal_obj, left_state)
      call write_tddft_state_artifact(trim(config%output_file)//'.state', reciprocal_obj, left_state, use_accepted_kspace_scf)

      allocate(endpoints(size(config%q_list, 2)))
      do iq = 1, size(config%q_list, 2)
         call lr_q_endpoint_from_reciprocal(reciprocal_obj, left_state, config%q_list(:, iq), endpoints(iq))
      end do
      if (trim(config%backend) == tddft_driver_backend_projected_chi0) then
         call run_tddft_projected_chi0(config, response_space, radial_bases, ground_states, left_state, endpoints, &
            reciprocal_obj, lattice_obj, hamiltonian_obj)
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_projected_mills) then
         if (.not. use_accepted_kspace_scf) then
            error stop 'DRESP-04 projected_mills requires the accepted k-space SCF handoff'
         end if
         call run_tddft_projected_mills(config, response_space, radial_bases, ground_states, left_state, endpoints, &
            reciprocal_obj, lattice_obj)
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_projected_juelich) then
         if (.not. use_accepted_kspace_scf) then
            error stop 'DRESP-05 projected_juelich requires the accepted k-space SCF handoff'
         end if
         call run_tddft_projected_juelich(config, response_space, radial_bases, ground_states, left_state, endpoints, &
            reciprocal_obj, lattice_obj)
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_product_lehmann) then
         ! TDVK-02R2 is a bare-response validation seam only.  It stops at
         ! the naturally prepared reciprocal handoff and never enters KXC,
         ! Goldstone, Dyson, loss, or the dense point-space result container.
         call run_tddft_product_bare_smoke(config, response_space, radial_bases, ground_states, left_state, endpoints, &
            reciprocal_obj)
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_compact_dyson) then
         if (.not. use_accepted_kspace_scf) then
            error stop 'BLOCKED — ACCEPTED-STATE CONTINUITY'
         end if
         call run_tddft_compact_dyson(config, response_space, radial_bases, ground_states, left_state, endpoints, &
            reciprocal_obj, lattice_obj, control_obj)
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_native_rsgf) then
         call native_provider%initialize(config%native_rsgf_provider, green_obj, recursion_obj, hamiltonian_obj, &
            lattice_obj, reciprocal_obj, radial_bases)
         allocate(native_site_positions(3, nsite))
         native_site_positions = 0.0_rp
         if (trim(config%interaction_route) == tddft_driver_route_direct_alsda) then
            call compute_accepted_pauli_magnetization(reciprocal_obj, lattice_obj%symbolic_atoms, lattice_obj%nbulk, &
               accepted_pauli_magnetization)
         end if
         call evaluate_tddft_runtime_sweep(config, response_space, radial_bases, ground_states, left_state, endpoints, result, &
            native_provider, native_provider%pairs, native_site_positions, accepted_pauli_magnetization, &
            lr_kxc_magnetization_source_pauli_accepted)
      else
         if (trim(config%interaction_route) == tddft_driver_route_direct_alsda) then
            call compute_accepted_pauli_magnetization(reciprocal_obj, lattice_obj%symbolic_atoms, lattice_obj%nbulk, &
               accepted_pauli_magnetization)
         end if
         call evaluate_tddft_runtime_sweep(config, response_space, radial_bases, ground_states, left_state, endpoints, result, &
            accepted_pauli_magnetization=accepted_pauli_magnetization, &
            accepted_pauli_magnetization_source=lr_kxc_magnetization_source_pauli_accepted)
      end if
      call write_tddft_production_output(config, result, control_obj, lattice_obj, ground_states, response_space, reciprocal_obj, &
         use_accepted_kspace_scf)
      ! The production result has already been serialized. Release the dense
      ! response matrices before the accepted-state owner leaves scope; this
      ! is important for compact smoke runs with a full radial mesh.
      if (allocated(result%ks_susceptibility)) deallocate(result%ks_susceptibility)
      if (allocated(result%enhanced_susceptibility)) deallocate(result%enhanced_susceptibility)
      if (allocated(result%loss_matrix)) deallocate(result%loss_matrix)
      if (allocated(result%q_list)) deallocate(result%q_list)
      if (allocated(result%frequencies)) deallocate(result%frequencies)
   end subroutine run_tddft_production

   !> DRESP-04 material seam.  The accepted projected site chi0 is consumed
   !> directly by the Mills interaction and site-space Dyson wrapper.
   subroutine run_tddft_projected_mills(config, response_space, radial_bases, ground_states, left_state, endpoints, &
                                        reciprocal_obj, lattice_obj)
      type(tddft_runtime_config), intent(in) :: config
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(radial_ground_state), target, intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj

      type(projected_site_spin_contract), target :: contract_d, contract_spd
      type(lmto_product_response_basis), target :: product_plus, product_minus
      type(projected_chi0_request) :: bare_request, opposite_request, raw_request
      type(projected_chi0_result) :: bare_result, opposite_bare_result, raw_bare_result
      type(projected_mills_interaction_result) :: mills_result
      type(projected_dyson_request) :: dyson_request, opposite_dyson_request, raw_dyson_request
      type(projected_dyson_result) :: dyson_result, opposite_dyson_result, raw_dyson_result
      type(projected_site_spin_contract), pointer :: contract
      type(lmto_product_response_basis), pointer :: product
      real(rp), allocatable :: moment(:)
      real(rp) :: eta_run, covariance_error, raw_mode_residual
      integer :: projection_index, nprojection, ieta, iq, ifrequency, i, j, unit
      integer :: gamma_index, positive_index, negative_index
      character(len=8) :: selector
      logical :: gamma_found, covariance_found

      if (lattice_obj%nrec /= size(ground_states)) then
         error stop 'DRESP-04 material seam: lattice/site provenance mismatch'
      end if
      if (response_space%response_lmax /= 4) then
         error stop 'DRESP-04 material seam: complete response_lmax=4 is required'
      end if
      call contract_d%initialize(response_space, radial_bases, 'd')
      call contract_spd%initialize(response_space, radial_bases, 'spd')
      call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      call product_minus%initialize(response_space, radial_bases, lmto_product_channel_minus, .true.)
      allocate(moment(size(ground_states)))
      gamma_index = find_gamma_q_index(config%q_list)
      gamma_found = gamma_index > 0
      positive_index = 0
      do iq = 1, size(config%q_list, 2)
         if (sum(abs(config%q_list(:, iq))) > 2.0e-12_rp) then
            positive_index = iq
            exit
         end if
      end do
      negative_index = 0
      if (positive_index > 0) negative_index = find_matching_q(config%q_list, -config%q_list(:, positive_index), 2.0e-12_rp)
      covariance_found = positive_index > 0 .and. negative_index > 0
      if (trim(config%projected_selector) == 'both') then
         nprojection = 2
      else
         nprojection = 1
      end if

      open(newunit=unit, file=trim(config%output_file), status='replace', action='write')
      write(unit, '(a)') '# DRESP-04 projected Mills/Stoner site-space RPA'
      write(unit, '(a)') '# accepted_state = one converged reciprocal k-space SCF state; no ladder re-SCF'
      write(unit, '(a)') '# backend = projected_mills'
      write(unit, '(a)') '# chi0_backend = DRESP-02 Lehmann direct site matrix'
      write(unit, '(a)') '# interaction_route = Mills/Stoner mean-field splitting projection'
      write(unit, '(a)') '# dyson_equation = (I-chi0*U) chi=chi0; certified LAPACK solve; no explicit inverse in production'
      write(unit, '(a)') '# splitting_convention = H=H0 I+B_sigma sigma_z; H_up-H_down=2 B_sigma; Vz=sigma_z'
      write(unit, '(a)') '# loss_convention = L=-(chi-chi^dagger)/(2*i*pi), ordinary site matrix'
      write(unit, '(a)') '# dresp03tg_comparison = not asserted; no same-q DRESP-03TG artifact is consumed by this bounded run'
      write(unit, '(a,l1)') '# goldstone_correction = ', .false.
      write(unit, '(a,l1)') '# juelich_sumrule = ', .false.
      write(unit, '(a,l1)') '# alsda_kernel = ', .false.
      write(unit, '(a)') '# columns: selector eta_Ry q_index omega_Ry row col bare_Re bare_Im enhanced_Re enhanced_Im loss_Re loss_Im min_sv max_sv cond min_abs_eig dyson_residual loss_trace -pi_im_trace'

      do projection_index = 1, nprojection
         if (nprojection == 1 .and. trim(config%projected_selector) == 'd') then
            selector = 'd'
            contract => contract_d
            product => product_plus
         else if (nprojection == 1 .and. trim(config%projected_selector) == 'spd') then
            selector = 'spd'
            contract => contract_spd
            product => product_plus
         else if (projection_index == 1) then
            selector = 'd'
            contract => contract_d
            product => product_plus
         else
            selector = 'spd'
            contract => contract_spd
            product => product_plus
         end if
         if (trim(config%channel) == lr_channel_minus) product => product_minus
         call contract%moment_from_operator(left_state%eigenvalues, left_state%eigenvectors, left_state%k_weights, &
            left_state%fermi_level, left_state%temperature, ground_states, moment)
         call evaluate_projected_mills_from_reciprocal(contract, radial_bases, reciprocal_obj, &
            left_state%fermi_level, moment, mills_result)
         if (trim(mills_result%classification) == 'UNSUPPORTED') then
            error stop 'DRESP-04 material seam: projected Mills mapping is unsupported'
         end if
         write(unit, '(a,a)') '# selector = ', trim(selector)
         write(unit, '(a,*(es24.16,1x))') '# projected_moment = ', mills_result%projected_moment
         write(unit, '(a,*(es24.16,1x))') '# projected_splitting_B_sigma = ', mills_result%projected_splitting
         write(unit, '(a,*(es24.16,1x))') '# U_Mills = ', mills_result%interaction_U
         write(unit, '(a,a)') '# scalarization_classification = ', trim(mills_result%classification)
         write(unit, '(a,es24.16)') '# scalarization_residual = ', mills_result%scalarization_residual
         write(unit, '(a,es24.16)') '# locality_residual = ', mills_result%locality_residual
         write(unit, '(a,es24.16)') '# scalar_fit_condition_number = ', mills_result%fit_condition_number
         write(unit, '(a,a)') '# interaction_provenance = ', trim(mills_result%provenance)
         write(unit, '(a,l1)') '# raw_gamma_diagnostic = ', gamma_found
         write(unit, '(a,l1)') '# q_channel_covariance_requested = ', covariance_found

         do ieta = 1, size(config%eta_values)
            eta_run = config%eta_values(ieta)
            do iq = 1, size(config%q_list, 2)
               bare_request%q = config%q_list(:, iq)
               bare_request%frequencies = config%frequencies
               bare_request%eta = eta_run
               bare_request%channel = config%channel
               bare_request%contract => contract
               bare_request%product_basis => product
               bare_request%electronic_state => left_state
               bare_request%q_endpoint_state => endpoints(iq)
               bare_request%diagnostics = .false.
               call evaluate_projected_lehmann_chi0(bare_request, bare_result)

               dyson_request%selector = selector
               dyson_request%q = config%q_list(:, iq)
               dyson_request%frequencies = config%frequencies
               dyson_request%eta = eta_run
               dyson_request%channel = config%channel
               dyson_request%interaction_U = mills_result%interaction_U
               dyson_request%bare_chi = bare_result%susceptibility
               dyson_request%interaction_provenance = mills_result%provenance
               dyson_request%bare_provenance = bare_result%provenance
               call evaluate_projected_dyson(dyson_request, dyson_result)
               do ifrequency = 1, size(config%frequencies)
                  do j = 1, contract%nsite
                     do i = 1, contract%nsite
                        write(unit, '(a,1x,a,1x,es24.16,1x,i0,1x,es24.16,1x,2(i0,1x),13(es24.16,1x))') &
                           'ROW', trim(selector), eta_run, iq, config%frequencies(ifrequency), i, j, &
                           real(bare_result%susceptibility(i,j,ifrequency),rp), aimag(bare_result%susceptibility(i,j,ifrequency)), &
                           real(dyson_result%enhanced_chi(i,j,ifrequency),rp), aimag(dyson_result%enhanced_chi(i,j,ifrequency)), &
                           real(dyson_result%loss_matrix(i,j,ifrequency),rp), aimag(dyson_result%loss_matrix(i,j,ifrequency)), &
                           dyson_result%denominator_min_singular_value(ifrequency), dyson_result%denominator_max_singular_value(ifrequency), &
                           dyson_result%condition_number(ifrequency), dyson_result%minimum_magnitude_eigenvalue(ifrequency), &
                           dyson_result%dyson_residual(ifrequency), dyson_result%loss_trace(ifrequency), &
                           dyson_result%minus_im_trace_over_pi(ifrequency)
                     end do
                  end do
               end do
            end do

            if (gamma_found) then
               raw_request%q = config%q_list(:, gamma_index)
               raw_request%frequencies = [0.0_rp]
               raw_request%eta = eta_run
               raw_request%channel = config%channel
               raw_request%contract => contract
               raw_request%product_basis => product
               raw_request%electronic_state => left_state
               raw_request%q_endpoint_state => endpoints(gamma_index)
               raw_request%diagnostics = .false.
               call evaluate_projected_lehmann_chi0(raw_request, raw_bare_result)
               raw_dyson_request%selector = selector
               raw_dyson_request%q = 0.0_rp
               raw_dyson_request%frequencies = [0.0_rp]
               raw_dyson_request%eta = eta_run
               raw_dyson_request%channel = config%channel
               raw_dyson_request%interaction_U = mills_result%interaction_U
               raw_dyson_request%bare_chi = raw_bare_result%susceptibility
               raw_dyson_request%interaction_provenance = mills_result%provenance
               raw_dyson_request%bare_provenance = raw_bare_result%provenance
               call evaluate_projected_dyson(raw_dyson_request, raw_dyson_result)
               raw_mode_residual = maxval(abs(matmul(raw_dyson_result%denominator(:, :, 1), &
                  cmplx(moment, 0.0_rp, rp))))
               write(unit, '(a,1x,a,1x,es24.16,1x,4(es24.16,1x))') 'RAW_GAMMA', trim(selector), eta_run, &
                  raw_dyson_result%denominator_min_singular_value(1), raw_dyson_result%minimum_magnitude_eigenvalue(1), &
                  raw_mode_residual, raw_dyson_result%dyson_residual(1)
            end if
         end do

         if (covariance_found .and. trim(config%channel) == 'chi_plus') then
            opposite_request%q = config%q_list(:, negative_index)
            opposite_request%frequencies = -config%frequencies
            opposite_request%eta = config%eta
            opposite_request%channel = lr_channel_minus
            opposite_request%contract => contract
            opposite_request%product_basis => product_minus
            opposite_request%electronic_state => left_state
            opposite_request%q_endpoint_state => endpoints(negative_index)
            opposite_request%diagnostics = .false.
            bare_request%q = config%q_list(:, positive_index)
            bare_request%frequencies = config%frequencies
            bare_request%eta = config%eta
            bare_request%channel = lr_channel_plus
            bare_request%product_basis => product_plus
            bare_request%q_endpoint_state => endpoints(positive_index)
            call evaluate_projected_lehmann_chi0(bare_request, bare_result)
            call evaluate_projected_lehmann_chi0(opposite_request, opposite_bare_result)
            dyson_request%selector = selector
            dyson_request%q = config%q_list(:, positive_index)
            dyson_request%frequencies = config%frequencies
            dyson_request%eta = config%eta
            dyson_request%channel = lr_channel_plus
            dyson_request%interaction_U = mills_result%interaction_U
            dyson_request%bare_chi = bare_result%susceptibility
            call evaluate_projected_dyson(dyson_request, dyson_result)
            opposite_dyson_request%selector = selector
            opposite_dyson_request%q = config%q_list(:, negative_index)
            opposite_dyson_request%frequencies = -config%frequencies
            opposite_dyson_request%eta = config%eta
            opposite_dyson_request%channel = lr_channel_minus
            opposite_dyson_request%interaction_U = mills_result%interaction_U
            opposite_dyson_request%bare_chi = opposite_bare_result%susceptibility
            call evaluate_projected_dyson(opposite_dyson_request, opposite_dyson_result)
            covariance_error = maxval(abs(dyson_result%enhanced_chi - conjg(opposite_dyson_result%enhanced_chi)))
            write(unit, '(a,1x,a,1x,es24.16)') 'Q_CHANNEL_COVARIANCE', trim(selector), covariance_error
         end if
      end do
      close(unit)
      deallocate(moment)
   end subroutine run_tddft_projected_mills

   !> DRESP-05 material seam.  The static Gamma response constructs a frozen
   !> real site-diagonal Juelich interaction; the dynamic loop then evaluates
   !> bare, Mills, and Juelich routes on exactly the same accepted state.
   subroutine run_tddft_projected_juelich(config, response_space, radial_bases, ground_states, left_state, endpoints, &
                                          reciprocal_obj, lattice_obj)
      type(tddft_runtime_config), intent(in) :: config
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(radial_ground_state), target, intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj

      type(projected_site_spin_contract), target :: contract_d, contract_spd
      type(lmto_product_response_basis), target :: product_plus, product_minus
      type(projected_chi0_request) :: chi_request, opposite_request
      type(projected_chi0_result) :: chi_result, opposite_result
      type(projected_mills_interaction_result) :: mills_result
      type(projected_juelich_request) :: juelich_request
      type(projected_juelich_result), allocatable :: static_ladder(:)
      type(projected_juelich_result) :: juelich_result
      type(projected_dyson_request) :: dyson_request, bare_dyson_request, opposite_dyson_request
      type(projected_dyson_result) :: dyson_result, bare_dyson_result, opposite_dyson_result
      type(projected_site_spin_contract), pointer :: contract
      type(lmto_product_response_basis), pointer :: product
      real(rp), allocatable :: moment(:)
      complex(rp), allocatable :: selected_static_chi(:, :, :)
      real(rp) :: eta_run, ward_mills, ward_juelich, covariance_mills, covariance_juelich
      real(rp) :: moment_norm
      integer :: projection_index, nprojection, ieta, iq, ifrequency, i, j, unit, gamma_index
      integer :: positive_index, negative_index, selected_eta_index, previous_eta_index
      character(len=8) :: selector
      logical :: covariance_found, eta_stable

      if (lattice_obj%nrec /= size(ground_states)) then
         error stop 'DRESP-05 material seam: lattice/site provenance mismatch'
      end if
      if (response_space%response_lmax /= 4) then
         error stop 'DRESP-05 material seam: complete response_lmax=4 is required'
      end if
      gamma_index = find_gamma_q_index(config%q_list)
      if (gamma_index == 0) error stop 'DRESP-05 material seam: Gamma is required for the projected Ward construction'
      call contract_d%initialize(response_space, radial_bases, 'd')
      call contract_spd%initialize(response_space, radial_bases, 'spd')
      call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      call product_minus%initialize(response_space, radial_bases, lmto_product_channel_minus, .true.)
      allocate(moment(size(ground_states)), selected_static_chi(size(ground_states), size(ground_states), size(config%eta_values)))
      call select_projected_juelich_eta_indices(config%eta_values, selected_eta_index, previous_eta_index)
      positive_index = 0
      do iq = 1, size(config%q_list, 2)
         if (sum(abs(config%q_list(:, iq))) > 2.0e-12_rp) then
            positive_index = iq
            exit
         end if
      end do
      negative_index = 0
      if (positive_index > 0) negative_index = find_matching_q(config%q_list, -config%q_list(:, positive_index), 2.0e-12_rp)
      covariance_found = positive_index > 0 .and. negative_index > 0
      if (trim(config%projected_selector) == 'both') then
         nprojection = 2
      else
         nprojection = 1
      end if

      open(newunit=unit, file=trim(config%output_file), status='replace', action='write')
      write(unit, '(a)') '# DRESP-05 projected Juelich/LCMM site-space response; same-state Mills comparison'
      write(unit, '(a)') '# accepted_state = one converged reciprocal k-space SCF state; no ladder re-SCF'
      write(unit, '(a)') '# backend = projected_juelich'
      write(unit, '(a)') '# selector_primary = spd; d is a controlled projection diagnostic'
      write(unit, '(a,i0)') '# nk = ', left_state%nk
      write(unit, '(a,es24.16)') '# temperature_K = ', left_state%temperature
      write(unit, '(a,es24.16)') '# fermi_level_Ry = ', left_state%fermi_level
      write(unit, '(a)') '# site_Ward_identity = Gamma(i,j)=chi0(i,j;0)*M(j); Gamma*U=M'
      write(unit, '(a)') '# radial_4pi_policy = not copied; site quantities use DRESP-01/02 normalization'
      write(unit, '(a)') '# U_Juelich is constructed once from the static eta ladder and frozen in dynamics'
      write(unit, '(a)') '# no empirical Goldstone correction; DRESP-03TG is frozen and not fitted'
      write(unit, '(a)') '# columns: route selector eta_Ry q_index omega_Ry row col bare_Re bare_Im chi_Re chi_Im loss_Re loss_Im min_sv max_sv cond min_abs_eig dyson_residual loss_trace minus_im_trace_over_pi'

      do projection_index = 1, nprojection
         if (nprojection == 1 .and. trim(config%projected_selector) == 'd') then
            selector = 'd'
            contract => contract_d
            product => product_plus
         else if (nprojection == 1 .and. trim(config%projected_selector) == 'spd') then
            selector = 'spd'
            contract => contract_spd
            product => product_plus
         else if (projection_index == 1) then
            selector = 'd'
            contract => contract_d
            product => product_plus
         else
            selector = 'spd'
            contract => contract_spd
            product => product_plus
         end if
         call contract%moment_from_operator(left_state%eigenvalues, left_state%eigenvectors, left_state%k_weights, &
            left_state%fermi_level, left_state%temperature, ground_states, moment)
         call evaluate_projected_mills_from_reciprocal(contract, radial_bases, reciprocal_obj, &
            left_state%fermi_level, moment, mills_result)
         if (trim(mills_result%classification) == 'UNSUPPORTED') then
            error stop 'DRESP-05 material seam: projected Mills comparison is unsupported'
         end if
         allocate(static_ladder(size(config%eta_values)))
         do ieta = 1, size(config%eta_values)
            chi_request%q = config%q_list(:, gamma_index)
            chi_request%frequencies = [0.0_rp]
            chi_request%eta = config%eta_values(ieta)
            chi_request%channel = lr_channel_plus
            chi_request%contract => contract
            chi_request%product_basis => product_plus
            chi_request%electronic_state => left_state
            chi_request%q_endpoint_state => endpoints(gamma_index)
            chi_request%diagnostics = .false.
            call evaluate_projected_lehmann_chi0(chi_request, chi_result)
            selected_static_chi(:, :, ieta) = chi_result%susceptibility(:, :, 1)
            juelich_request%selector = selector
            juelich_request%projected_moment = moment
            juelich_request%static_chi0 = chi_result%susceptibility(:, :, 1)
            juelich_request%static_eta = config%eta_values(ieta)
            juelich_request%q = config%q_list(:, gamma_index)
            juelich_request%channel = lr_channel_plus
            juelich_request%state_provenance = 'accepted reciprocal state; DRESP-01 projected moment'
            juelich_request%chi0_provenance = chi_result%provenance
            call evaluate_projected_juelich_interaction(juelich_request, static_ladder(ieta))
         end do
         juelich_result = static_ladder(selected_eta_index)
         if (previous_eta_index > 0) then
            call assess_projected_juelich_eta_stability(juelich_result, static_ladder(previous_eta_index), eta_stable)
            call evaluate_projected_juelich_holdout(selected_static_chi(:, :, previous_eta_index), moment, &
               juelich_result%interaction_U_real, juelich_result%holdout_residual, juelich_result%holdout_relative_residual)
         else
            juelich_result%eta_stable = .false.
            juelich_result%eta_stability = 'not independently assessed: one static eta'
            juelich_result%holdout_residual = juelich_result%real_constrained_residual
            juelich_result%holdout_relative_residual = juelich_result%relative_real_constrained_residual
         end if
         if (trim(juelich_result%classification) == projected_juelich_rank_deficient .or. &
             trim(juelich_result%classification) == projected_juelich_unsupported) then
            error stop 'DRESP-05 material seam: projected Juelich static solve is unsupported'
         end if
         moment_norm = sqrt(sum(moment**2))
         ward_mills = sqrt(sum(abs(cmplx(moment, 0.0_rp, rp) - matmul(chi_result%susceptibility(:, :, 1), &
            cmplx(mills_result%interaction_U*moment, 0.0_rp, rp)))**2))/max(moment_norm, tiny(1.0_rp))
         ward_juelich = juelich_result%relative_real_constrained_residual
         write(unit, '(a,a)') '# selector = ', trim(selector)
         write(unit, '(a,*(es24.16,1x))') '# projected_moment = ', moment
         write(unit, '(a,*(es24.16,1x))') '# U_Mills = ', mills_result%interaction_U
         write(unit, '(a,es24.16)') '# Mills_scalarization_residual = ', mills_result%scalarization_residual
         write(unit, '(a,a)') '# Mills_classification = ', trim(mills_result%classification)
         write(unit, '(a,*(es24.16,1x))') '# U_Juelich_complex_Re = ', real(juelich_result%interaction_U_complex, rp)
         write(unit, '(a,*(es24.16,1x))') '# U_Juelich_complex_Im = ', aimag(juelich_result%interaction_U_complex)
         write(unit, '(a,*(es24.16,1x))') '# U_Juelich_real = ', juelich_result%interaction_U_real
         write(unit, '(a,a)') '# Juelich_classification = ', trim(juelich_result%classification)
         write(unit, '(a,i0)') '# Juelich_rank = ', juelich_result%rank
         write(unit, '(a,*(es24.16,1x))') '# Juelich_singular_values = ', juelich_result%singular_values
         write(unit, '(a,*(es24.16,1x))') '# Juelich_real_singular_values = ', juelich_result%real_singular_values
         write(unit, '(a,es24.16)') '# Juelich_condition_number = ', juelich_result%condition_number
         write(unit, '(a,es24.16)') '# Juelich_imaginary_U_ratio = ', juelich_result%imaginary_U_ratio
         write(unit, '(a,es24.16)') '# Juelich_complex_construction_residual = ', juelich_result%relative_residual
         write(unit, '(a,es24.16)') '# Juelich_construction_residual = ', juelich_result%relative_real_constrained_residual
         write(unit, '(a,es24.16)') '# Juelich_holdout_residual = ', juelich_result%holdout_relative_residual
         write(unit, '(a,l1)') '# Juelich_eta_stable = ', juelich_result%eta_stable
         write(unit, '(a,a)') '# Juelich_eta_stability = ', trim(juelich_result%eta_stability)
         write(unit, '(a,*(es24.16,1x))') '# relative_U_difference_Juelich_minus_Mills = ', &
            (juelich_result%interaction_U_real - mills_result%interaction_U)/ &
            sign(max(abs(mills_result%interaction_U), tiny(1.0_rp)), mills_result%interaction_U)
         write(unit, '(a,es24.16)') '# Mills_raw_Gamma_residual = ', ward_mills
         write(unit, '(a,es24.16)') '# Juelich_raw_Gamma_residual = ', ward_juelich
         write(unit, '(a,i0)') '# static_eta_count = ', size(config%eta_values)
         write(unit, '(a)') '# JUELICH_STATIC columns: selector eta U_real U_complex_Re U_complex_Im imaginary_ratio real_residual condition_number'
         do ieta = 1, size(config%eta_values)
            write(unit, '(a,1x,a,1x,es24.16,1x,*(es24.16,1x))') 'JUELICH_STATIC', trim(selector), &
               config%eta_values(ieta), static_ladder(ieta)%interaction_U_real, &
               real(static_ladder(ieta)%interaction_U_complex, rp), aimag(static_ladder(ieta)%interaction_U_complex), &
               static_ladder(ieta)%imaginary_U_ratio, static_ladder(ieta)%relative_real_constrained_residual, &
               static_ladder(ieta)%condition_number
         end do

         do ieta = 1, size(config%eta_values)
            eta_run = config%eta_values(ieta)
            do iq = 1, size(config%q_list, 2)
               chi_request%q = config%q_list(:, iq)
               chi_request%frequencies = config%frequencies
               chi_request%eta = eta_run
               chi_request%channel = config%channel
               chi_request%contract => contract
               if (trim(config%channel) == lr_channel_minus) then
                  chi_request%product_basis => product_minus
               else
                  chi_request%product_basis => product_plus
               end if
               chi_request%electronic_state => left_state
               chi_request%q_endpoint_state => endpoints(iq)
               chi_request%diagnostics = .false.
               call evaluate_projected_lehmann_chi0(chi_request, chi_result)

               bare_dyson_request%selector = selector
               bare_dyson_request%q = config%q_list(:, iq)
               bare_dyson_request%frequencies = config%frequencies
               bare_dyson_request%eta = eta_run
               bare_dyson_request%channel = config%channel
               allocate(bare_dyson_request%interaction_U(contract%nsite))
               bare_dyson_request%interaction_U = 0.0_rp
               bare_dyson_request%bare_chi = chi_result%susceptibility
               bare_dyson_request%interaction_provenance = 'bare U=0 comparison route'
               bare_dyson_request%bare_provenance = chi_result%provenance
               call evaluate_projected_dyson(bare_dyson_request, bare_dyson_result)
               call write_projected_juelich_rows(unit, 'bare', selector, eta_run, iq, contract%nsite, &
                  chi_result, bare_dyson_result)
               deallocate(bare_dyson_request%interaction_U)

               dyson_request%selector = selector
               dyson_request%q = config%q_list(:, iq)
               dyson_request%frequencies = config%frequencies
               dyson_request%eta = eta_run
               dyson_request%channel = config%channel
               dyson_request%interaction_U = mills_result%interaction_U
               dyson_request%bare_chi = chi_result%susceptibility
               dyson_request%interaction_provenance = 'same-state DRESP-04 Mills interaction'
               dyson_request%bare_provenance = chi_result%provenance
               call evaluate_projected_dyson(dyson_request, dyson_result)
               call write_projected_juelich_rows(unit, 'Mills', selector, eta_run, iq, contract%nsite, &
                  chi_result, dyson_result)

               dyson_request%interaction_U = juelich_result%interaction_U_real
               dyson_request%interaction_provenance = juelich_result%provenance
               call evaluate_projected_dyson(dyson_request, dyson_result)
               call write_projected_juelich_rows(unit, 'Juelich', selector, eta_run, iq, contract%nsite, &
                  chi_result, dyson_result)
            end do
         end do

         if (positive_index > 0) then
            call run_projected_frequency_refinement(unit, selector, config, positive_index, contract, product_plus, &
               left_state, endpoints(positive_index), mills_result%interaction_U, juelich_result%interaction_U_real)
         end if

         if (covariance_found .and. trim(config%channel) == lr_channel_plus) then
            chi_request%q = config%q_list(:, positive_index)
            chi_request%frequencies = config%frequencies
            chi_request%eta = config%eta
            chi_request%channel = lr_channel_plus
            chi_request%contract => contract
            chi_request%product_basis => product_plus
            chi_request%electronic_state => left_state
            chi_request%q_endpoint_state => endpoints(positive_index)
            call evaluate_projected_lehmann_chi0(chi_request, chi_result)
            opposite_request%q = config%q_list(:, negative_index)
            opposite_request%frequencies = -config%frequencies
            opposite_request%eta = config%eta
            opposite_request%channel = lr_channel_minus
            opposite_request%contract => contract
            opposite_request%product_basis => product_minus
            opposite_request%electronic_state => left_state
            opposite_request%q_endpoint_state => endpoints(negative_index)
            call evaluate_projected_lehmann_chi0(opposite_request, opposite_result)
            dyson_request%selector = selector
            dyson_request%q = config%q_list(:, positive_index)
            dyson_request%frequencies = config%frequencies
            dyson_request%eta = config%eta
            dyson_request%channel = lr_channel_plus
            dyson_request%bare_chi = chi_result%susceptibility
            dyson_request%interaction_U = mills_result%interaction_U
            call evaluate_projected_dyson(dyson_request, dyson_result)
            opposite_dyson_request%selector = selector
            opposite_dyson_request%q = config%q_list(:, negative_index)
            opposite_dyson_request%frequencies = -config%frequencies
            opposite_dyson_request%eta = config%eta
            opposite_dyson_request%channel = lr_channel_minus
            opposite_dyson_request%bare_chi = opposite_result%susceptibility
            opposite_dyson_request%interaction_U = mills_result%interaction_U
            call evaluate_projected_dyson(opposite_dyson_request, opposite_dyson_result)
            covariance_mills = maxval(abs(dyson_result%enhanced_chi - conjg(opposite_dyson_result%enhanced_chi)))
            dyson_request%interaction_U = juelich_result%interaction_U_real
            call evaluate_projected_dyson(dyson_request, dyson_result)
            opposite_dyson_request%interaction_U = juelich_result%interaction_U_real
            call evaluate_projected_dyson(opposite_dyson_request, opposite_dyson_result)
            covariance_juelich = maxval(abs(dyson_result%enhanced_chi - conjg(opposite_dyson_result%enhanced_chi)))
            write(unit, '(a,1x,a,1x,es24.16)') 'Q_CHANNEL_COVARIANCE', trim(selector)//'_Mills', covariance_mills
            write(unit, '(a,1x,a,1x,es24.16)') 'Q_CHANNEL_COVARIANCE', trim(selector)//'_Juelich', covariance_juelich
         end if
         deallocate(static_ladder)
      end do
      close(unit)
      deallocate(moment, selected_static_chi)
   end subroutine run_tddft_projected_juelich

   !> DRESP-06A same-state comparison seam.  Mills and Juelich remain site
   !> scalar comparison routes; the ALSDA route is solved in the complete
   !> weighted-orthonormal product space and is projected only for observables.

   !> Refine the ALSDA loss peak on a window selected from the ALSDA coarse
   !> response itself.  This is intentionally independent from the Mills and
   !> Juelich projected-site windows.



   !> One representative finite-q bare product-GF spot for the ALSDA
   !> comparison.  It is a closure audit of the frozen accepted state, not an
   !> interaction-specific GF construction.

   subroutine write_projected_juelich_rows(unit, route, selector, eta, q_index, nsite, chi0, dyson)
      integer, intent(in) :: unit, q_index, nsite
      character(len=*), intent(in) :: route, selector
      real(rp), intent(in) :: eta
      type(projected_chi0_result), intent(in) :: chi0
      type(projected_dyson_result), intent(in) :: dyson
      integer :: ifrequency, i, j

      do ifrequency = 1, size(chi0%frequencies)
         do j = 1, nsite
            do i = 1, nsite
               write(unit, '(a,1x,a,1x,a,1x,es24.16,1x,i0,1x,es24.16,1x,2(i0,1x),13(es24.16,1x))') &
                  'ROW', trim(route), trim(selector), eta, q_index, chi0%frequencies(ifrequency), i, j, &
                  real(chi0%susceptibility(i,j,ifrequency),rp), aimag(chi0%susceptibility(i,j,ifrequency)), &
                  real(dyson%enhanced_chi(i,j,ifrequency),rp), aimag(dyson%enhanced_chi(i,j,ifrequency)), &
                  real(dyson%loss_matrix(i,j,ifrequency),rp), aimag(dyson%loss_matrix(i,j,ifrequency)), &
                  dyson%denominator_min_singular_value(ifrequency), dyson%denominator_max_singular_value(ifrequency), &
                  dyson%condition_number(ifrequency), dyson%minimum_magnitude_eigenvalue(ifrequency), &
                  dyson%dyson_residual(ifrequency), dyson%loss_trace(ifrequency), dyson%minus_im_trace_over_pi(ifrequency)
            end do
         end do
      end do
   end subroutine write_projected_juelich_rows

   subroutine run_projected_frequency_refinement(unit, selector, config, q_index, contract, product_plus, &
                                                 left_state, endpoint, mills_U, juelich_U)
      integer, intent(in) :: unit, q_index
      character(len=*), intent(in) :: selector
      type(tddft_runtime_config), intent(in) :: config
      type(projected_site_spin_contract), target, intent(in) :: contract
      type(lmto_product_response_basis), target, intent(in) :: product_plus
      type(lr_electronic_state), target, intent(in) :: left_state, endpoint
      real(rp), intent(in) :: mills_U(:), juelich_U(:)
      type(projected_chi0_request) :: request
      type(projected_chi0_result) :: coarse_chi, refined_chi_mills, refined_chi_juelich
      type(projected_dyson_request) :: dyson_request
      type(projected_dyson_result) :: mills_coarse, juelich_coarse, mills_refined, juelich_refined
      real(rp), allocatable :: refined_frequencies_mills(:), refined_frequencies_juelich(:)
      real(rp) :: step, lower_mills, upper_mills, lower_juelich, upper_juelich, eta
      integer :: nref, i, coarse_peak_index, coarse_juelich_peak_index, refined_peak_index
      integer :: selected_eta_index, holdout_eta_index

      if (size(config%frequencies) < 2) return
      call select_projected_juelich_eta_indices(config%eta_values, selected_eta_index, holdout_eta_index)
      eta = config%eta_values(selected_eta_index)
      request%q = config%q_list(:, q_index)
      request%frequencies = config%frequencies
      request%eta = eta
      request%channel = lr_channel_plus
      request%contract => contract
      request%product_basis => product_plus
      request%electronic_state => left_state
      request%q_endpoint_state => endpoint
      request%diagnostics = .false.
      call evaluate_projected_lehmann_chi0(request, coarse_chi)
      dyson_request%selector = selector
      dyson_request%q = config%q_list(:, q_index)
      dyson_request%frequencies = config%frequencies
      dyson_request%eta = eta
      dyson_request%channel = lr_channel_plus
      dyson_request%bare_chi = coarse_chi%susceptibility
      dyson_request%interaction_U = mills_U
      call evaluate_projected_dyson(dyson_request, mills_coarse)
      dyson_request%interaction_U = juelich_U
      call evaluate_projected_dyson(dyson_request, juelich_coarse)
      coarse_peak_index = maxloc(mills_coarse%loss_trace, dim=1)
      coarse_juelich_peak_index = maxloc(juelich_coarse%loss_trace, dim=1)
      step = minval(abs(config%frequencies(2:) - config%frequencies(:size(config%frequencies)-1)))
      lower_mills = max(minval(config%frequencies), config%frequencies(coarse_peak_index) - step)
      upper_mills = min(maxval(config%frequencies), config%frequencies(coarse_peak_index) + step)
      lower_juelich = max(minval(config%frequencies), config%frequencies(coarse_juelich_peak_index) - step)
      upper_juelich = min(maxval(config%frequencies), config%frequencies(coarse_juelich_peak_index) + step)
      if (upper_mills <= lower_mills .or. upper_juelich <= lower_juelich) return
      nref = 41
      allocate(refined_frequencies_mills(nref), refined_frequencies_juelich(nref))
      do i = 1, nref
         refined_frequencies_mills(i) = lower_mills + real(i - 1, rp)*(upper_mills - lower_mills)/real(nref - 1, rp)
         refined_frequencies_juelich(i) = lower_juelich + real(i - 1, rp)*(upper_juelich - lower_juelich)/real(nref - 1, rp)
      end do
      ! Each route receives its own refinement grid.  The two coarse maxima
      ! need not coincide, so reusing the Mills window for Juelich would be a
      ! route-dependent frequency-selection error.
      request%frequencies = refined_frequencies_mills
      call evaluate_projected_lehmann_chi0(request, refined_chi_mills)
      dyson_request%frequencies = refined_frequencies_mills
      dyson_request%bare_chi = refined_chi_mills%susceptibility
      dyson_request%interaction_U = mills_U
      call evaluate_projected_dyson(dyson_request, mills_refined)
      request%frequencies = refined_frequencies_juelich
      call evaluate_projected_lehmann_chi0(request, refined_chi_juelich)
      dyson_request%frequencies = refined_frequencies_juelich
      dyson_request%bare_chi = refined_chi_juelich%susceptibility
      dyson_request%interaction_U = juelich_U
      call evaluate_projected_dyson(dyson_request, juelich_refined)
      refined_peak_index = maxloc(mills_refined%loss_trace, dim=1)
      write(unit, '(a,1x,a,1x,i0,1x,4(es24.16,1x))') 'FREQUENCY_REFINEMENT', trim(selector)//'_Mills', q_index, &
         config%frequencies(coarse_peak_index), refined_frequencies_mills(refined_peak_index), &
         minval(mills_refined%denominator_min_singular_value), mills_refined%loss_trace(refined_peak_index)
      refined_peak_index = maxloc(juelich_refined%loss_trace, dim=1)
      write(unit, '(a,1x,a,1x,i0,1x,4(es24.16,1x))') 'FREQUENCY_REFINEMENT', trim(selector)//'_Juelich', q_index, &
         config%frequencies(coarse_juelich_peak_index), refined_frequencies_juelich(refined_peak_index), &
         minval(juelich_refined%denominator_min_singular_value), juelich_refined%loss_trace(refined_peak_index)
      deallocate(refined_frequencies_mills, refined_frequencies_juelich)
   end subroutine run_projected_frequency_refinement

   !> DRESP-02 material seam.  Both projections consume the same accepted
   !> reciprocal state and stop at bare chi0; no interaction or Dyson object
   !> is constructed here.
   subroutine run_tddft_projected_chi0(config, response_space, radial_bases, ground_states, left_state, endpoints, &
                                      reciprocal_obj, lattice_obj, hamiltonian_obj)
      type(tddft_runtime_config), intent(in) :: config
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(radial_ground_state), target, intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(reciprocal), intent(inout) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj

      type(projected_site_spin_contract), target :: contract_d, contract_spd
      type(lmto_product_response_basis), target :: product_plus, product_minus
      type(projected_chi0_request) :: request, minus_request
      type(projected_chi0_result) :: lehmann_result, minus_result
      type(projected_site_spin_contract), pointer :: contract
      type(lmto_product_response_basis), pointer :: product
      real(rp), allocatable :: moment(:), accepted_moment(:)
      real(rp) :: difference, accepted_total, moment_residual
      real(rp) :: endpoint_eigenvalue_checksum, endpoint_occupation_checksum, endpoint_unitarity_residual
      complex(rp) :: endpoint_eigenvector_checksum
      real(rp) :: state_eigenvalue_checksum, state_occupation_checksum, state_unitarity_residual
      complex(rp) :: state_eigenvector_checksum
      real(rp) :: radial_checksum, radial_l2_norm
      integer :: projection_index, iq, ifrequency, i, j, orbital, unit, positive_index, negative_index
      logical :: covariance_saved
      character(len=8) :: projection

      ! lattice_obj is part of the material provenance boundary even though
      ! the site-space result itself needs only the accepted response sites.
      if (lattice_obj%nrec /= size(ground_states)) then
         error stop 'DRESP-02 material seam: lattice/site provenance mismatch'
      end if
      if (response_space%response_lmax /= 4) then
         error stop 'DRESP-02 material seam: complete response_lmax=4 is required'
      end if
      call contract_d%initialize(response_space, radial_bases, 'd')
      call contract_spd%initialize(response_space, radial_bases, 'spd')
      call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      call product_minus%initialize(response_space, radial_bases, lmto_product_channel_minus, .true.)
      allocate(moment(size(ground_states)))
      allocate(accepted_moment(size(ground_states)))
      if (.not. allocated(reciprocal_obj%band_moments)) then
         ! finalize_kspace_scf_state intentionally invalidates DOS-derived
         ! arrays. Recreate only this diagnostic projection on the same
         ! accepted eigensystem; calculate_density_of_states does not rebuild
         ! H(k) when those eigenpairs are resident.
         call reciprocal_obj%calculate_density_of_states(hamiltonian_obj, &
            n_energy_points=reciprocal_obj%n_energy_points, &
            energy_range=reciprocal_obj%dos_energy_range, &
            method=reciprocal_obj%dos_method, gaussian_sigma=reciprocal_obj%gaussian_sigma, &
            temperature=reciprocal_obj%temperature, fermi_level=left_state%fermi_level, &
            total_electrons=reciprocal_obj%total_electrons, auto_find_fermi=.false., &
            output_file=trim(config%output_file)//'.accepted_state_dos')
      end if
      if (size(reciprocal_obj%band_moments, 1) < size(ground_states) .or. &
          size(reciprocal_obj%band_moments, 2) < 3 .or. size(reciprocal_obj%band_moments, 3) < 2 .or. &
          size(reciprocal_obj%band_moments, 4) < 1) then
         error stop 'DRESP-02 material seam: accepted band moment dimensions are incomplete'
      end if
      accepted_total = 0.0_rp
      do i = 1, size(ground_states)
         ! The accepted reciprocal SCF state stores its integrated moment on
         ! the accepted potential. reported_moment is populated only by the
         ! later human-readable report, and the radial quadrature integral is
         ! retained as a separate cross-check rather than the gate source.
         accepted_total = accepted_total + lattice_obj%symbolic_atoms(lattice_obj%nbulk + i)%potential%mtot
      end do

      positive_index = 0
      do iq = 1, size(config%q_list, 2)
         if (sum(abs(config%q_list(:, iq))) > 2.0e-12_rp) then
            positive_index = iq
            exit
         end if
      end do
      negative_index = 0
      if (positive_index > 0) negative_index = find_matching_q(config%q_list, -config%q_list(:, positive_index), 2.0e-12_rp)
      covariance_saved = positive_index > 0 .and. negative_index > 0

      open(newunit=unit, file=trim(config%output_file), status='replace', action='write')
      write(unit, '(a)') '# DRESP-02 projected reciprocal bare susceptibility'
      write(unit, '(a)') '# verdict = algebraic projected-site Lehmann bare susceptibility'
      write(unit, '(a)') '# state_source = one accepted k-space SCF reciprocal state'
      write(unit, '(a,l1)') '# accepted_state_cache_reused = ', .true.
      write(unit, '(a,3(i0,1x))') '# accepted_k_mesh = ', reciprocal_obj%nk_mesh
      write(unit, '(a,i0)') '# accepted_k_count = ', left_state%nk
      write(unit, '(a,es24.16)') '# structure_alat = ', lattice_obj%alat
      write(unit, '(a,i0)') '# structure_ntype = ', lattice_obj%ntype
      write(unit, '(a,i0)') '# structure_nrec = ', lattice_obj%nrec
      write(unit, '(a,i0)') '# accepted_response_lmax = ', response_space%response_lmax
      write(unit, '(a,i0)') '# accepted_radial_sites = ', size(ground_states)
      write(unit, '(a)') '# accepted_hamiltonian = reciprocal ham_only eigensystem from the accepted k-space SCF cache'
      write(unit, '(a)') '# accepted_lmto_radial_state = accepted radial/Pauli snapshots; no radial rebuild in the GF ladders'
      write(unit, '(a)') '# spin_state = collinear, orthogonal, no SOC, no extra operator'
      write(unit, '(a)') '# material_state_gate = structure, SCF cache, Hamiltonian/eigenpairs, LMTO/radial, k mesh/weights, EF/occupations, and DRESP projection are frozen and shared'
      write(unit, '(a,es24.16)') '# EF_Ry = ', left_state%fermi_level
      write(unit, '(a,es24.16)') '# temperature_K = ', left_state%temperature
      write(unit, '(a,es24.16)') '# accepted_state_integrated_moment_muB = ', accepted_total
      call state_identity_checksums(left_state, state_eigenvalue_checksum, state_eigenvector_checksum, &
         state_occupation_checksum, state_unitarity_residual)
      call radial_identity_checksums(ground_states, radial_checksum, radial_l2_norm)
      write(unit, '(a,es24.16)') '# accepted_state_eigenvalue_checksum_Ry = ', state_eigenvalue_checksum
      write(unit, '(a,2(es24.16,1x))') '# accepted_state_eigenvector_checksum_real_imag = ', &
         real(state_eigenvector_checksum, rp), aimag(state_eigenvector_checksum)
      write(unit, '(a,es24.16)') '# accepted_state_occupation_checksum = ', state_occupation_checksum
      write(unit, '(a,es24.16)') '# accepted_state_eigenvector_unitarity_max_abs = ', state_unitarity_residual
      write(unit, '(a,es24.16)') '# accepted_radial_checksum = ', radial_checksum
      write(unit, '(a,es24.16)') '# accepted_radial_l2_norm = ', radial_l2_norm
      write(unit, '(a,a)') '# reciprocal_mode = ', trim(left_state%reciprocal_mode)
      write(unit, '(a,a)') '# hamiltonian_order = ', trim(left_state%hamiltonian_order)
      write(unit, '(a,es24.16)') '# eta_response_Ry = ', config%eta
      write(unit, '(a)') '# q_convention = exact folded reciprocal k+q endpoint; no extra DRESP site phase'
      write(unit, '(a)') '# q_endpoint_gate = exact folded k+q endpoint state; endpoint checksums below are compared before each response sample'
      write(unit, '(a)') '# scf_during_ladder = F'
      write(unit, '(a)') '# endpoint_identity columns: q_index qx qy qz eigenvalue_checksum_Ry eigenvector_checksum_real eigenvector_checksum_imag occupation_checksum unitarity_max_abs'
      do iq = 1, size(config%q_list, 2)
         call state_identity_checksums(endpoints(iq), endpoint_eigenvalue_checksum, endpoint_eigenvector_checksum, &
            endpoint_occupation_checksum, endpoint_unitarity_residual)
         write(unit, '(a,1x,i0,1x,3(es24.16,1x),5(es24.16,1x))') '# endpoint_identity', iq, config%q_list(:, iq), &
            endpoint_eigenvalue_checksum, real(endpoint_eigenvector_checksum, rp), aimag(endpoint_eigenvector_checksum), &
            endpoint_occupation_checksum, endpoint_unitarity_residual
      end do
      write(unit, '(a)') '# columns = projection q_index omega_Ry row col Lehmann_Re Lehmann_Im'

      do projection_index = 1, 2
         if (projection_index == 1) then
            projection = 'd'
            contract => contract_d
         else
            projection = 'spd'
            contract => contract_spd
         end if
         if (trim(config%channel) == 'chi_plus') then
            product => product_plus
         else
            product => product_minus
         end if
         call contract%moment_from_operator(left_state%eigenvalues, left_state%eigenvectors, left_state%k_weights, &
            left_state%fermi_level, left_state%temperature, ground_states, moment)
         accepted_moment = 0.0_rp
         do i = 1, size(ground_states)
            if (projection_index == 1) then
               accepted_moment(i) = reciprocal_obj%band_moments(i, 3, 1, 1) - &
                  reciprocal_obj%band_moments(i, 3, 2, 1)
            else
               do orbital = 1, 3
                  accepted_moment(i) = accepted_moment(i) + reciprocal_obj%band_moments(i, orbital, 1, 1) - &
                     reciprocal_obj%band_moments(i, orbital, 2, 1)
               end do
            end if
         end do
         if (.not. all(ieee_is_finite(accepted_moment))) then
            error stop 'DRESP-02 material seam: projected band moments are not finite'
         end if
         if (projection_index == 2) then
            moment_residual = abs(sum(accepted_moment) - accepted_total)
            if (moment_residual > 2.0e-5_rp) then
               error stop 'DRESP-02 material seam: accepted spd/total moment check failed'
            end if
         else
            moment_residual = 0.0_rp
         end if
         write(unit, '(a,a)') '# projection = ', trim(projection)
         write(unit, '(a,*(es24.16,1x))') '# projected_moment_muB = ', accepted_moment
         write(unit, '(a,*(es24.16,1x))') '# dresp_operator_moment_muB = ', moment
         write(unit, '(a,es24.16)') '# accepted_spd_total_residual_muB = ', moment_residual
         write(unit, '(a)') '# core_policy = valence-only; frozen core excluded'
         do iq = 1, size(config%q_list, 2)
            request%q = config%q_list(:, iq)
            request%frequencies = config%frequencies
            request%eta = config%eta
            request%channel = lr_channel_plus
            request%contract => contract
            request%product_basis => product
            request%electronic_state => left_state
            request%q_endpoint_state => endpoints(iq)
            request%channel = config%channel
            call evaluate_projected_lehmann_chi0(request, lehmann_result)
            do ifrequency = 1, size(config%frequencies)
               do j = 1, contract%nsite
                  do i = 1, contract%nsite
                     write(unit, '(a,1x,i0,1x,es24.16,1x,2(i0,1x),2(es24.16,1x))') trim(projection), iq, &
                        config%frequencies(ifrequency), i, j, real(lehmann_result%susceptibility(i, j, ifrequency), rp), &
                        aimag(lehmann_result%susceptibility(i, j, ifrequency))
                  end do
               end do
            end do
         end do

         if (covariance_saved) then
            minus_request%q = config%q_list(:, negative_index)
            minus_request%frequencies = -config%frequencies
            minus_request%eta = config%eta
            minus_request%channel = lr_channel_minus
            minus_request%contract => contract
            minus_request%product_basis => product_minus
            minus_request%electronic_state => left_state
            minus_request%q_endpoint_state => endpoints(negative_index)
            call evaluate_projected_lehmann_chi0(minus_request, minus_result)
            request%q = config%q_list(:, positive_index)
            request%frequencies = config%frequencies
            request%channel = config%channel
            request%product_basis => product
            request%q_endpoint_state => endpoints(positive_index)
            call evaluate_projected_lehmann_chi0(request, lehmann_result)
            difference = maxval(abs(lehmann_result%susceptibility - conjg(minus_result%susceptibility)))
            write(unit, '(a,a,es24.16)') '# q_minus_q_covariance_max_abs = ', trim(projection), difference
         end if

      end do
      close(unit)
      deallocate(moment)
      deallocate(accepted_moment)
   end subroutine run_tddft_projected_chi0

   !> Evaluate the compact direct-ALSDA Dyson response.
   !>
   !> This is the production worker for one accepted state.  The response
   !> services are evaluated q-by-q and serialized while their q-local result
   !> objects are live.  Validation audits are explicit opt-ins and add work
   !> around this worker; they do not define its production input contract.
   subroutine run_tddft_compact_dyson(config, response_space, radial_bases, ground_states, left_state, endpoints, &
                                      reciprocal_obj, lattice_obj, control_obj)
      type(tddft_runtime_config), intent(in) :: config
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(radial_ground_state), target, intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      type(control), intent(in) :: control_obj

      type(lmto_product_response_basis), target :: product_plus, product_minus
      type(lmto_product_response_basis), pointer :: product, opposite_product
      type(lr_product_ks_susceptibility_request) :: bare_request, opposite_bare_request
      type(lr_product_ks_susceptibility_result) :: bare_result, opposite_bare_result
      type(lr_alsda_kernel_request) :: kxc_request
      type(lr_alsda_kernel_result) :: kxc_result
      type(tddft_dyson_request) :: dyson_request, opposite_dyson_request
      type(tddft_dyson_result) :: dyson_result, opposite_dyson_result
      real(rp), allocatable :: magnetization(:, :), static_frequency(:)
      complex(rp), allocatable :: interaction(:, :), opposite_interaction(:, :)
      real(rp) :: static_eta(2), static_min_sv(2), static_max_sv(2), static_condition(2), static_min_eigen(2)
      real(rp) :: static_residual(2), static_residual_relative(2), static_residual_infinity(2)
      character(len=256) :: static_status(2)
      real(rp), allocatable :: covariance_residual(:), covariance_loss_difference(:), covariance_min_sv_difference(:), &
         covariance_condition_difference(:)
      logical :: rank_stable, finite_response, covariance_checked
      real(rp) :: state_mesh_max, state_weight_max, state_ef_diff, state_eigen_max, state_occ_max
      real(rp) :: state_projector_max, state_projector_frobenius, k_fingerprint(5), weight_sum, magnetic_moment
      real(rp) :: trace_loss, trace_loss_imag
      integer :: ndim, nfrequency, nq, iq, iw, i, j, unit, gamma_index, positive_q_index, negative_q_index
      character(len=256) :: state_file
      logical :: static_audit_enabled, covariance_enabled

      if (trim(config%interaction_route) /= tddft_driver_route_direct_alsda) then
         error stop 'compact_dyson requires direct ALSDA as the production interaction route'
      end if
      if (config%goldstone_correction) then
         error stop 'compact_dyson does not apply an implicit BES/GCR/Goldstone correction'
      end if
      call state_consistency_metrics(reciprocal_obj, left_state, state_mesh_max, state_weight_max, state_ef_diff, &
         state_eigen_max, state_occ_max, state_projector_max, state_projector_frobenius)
      if (max(state_mesh_max, state_weight_max, state_ef_diff, state_eigen_max, state_occ_max, state_projector_max, &
          state_projector_frobenius) > 5.0e-11_rp) then
         error stop 'compact_dyson: accepted-state continuity contract failed'
      end if

      static_audit_enabled = config%dyson_static_audit
      covariance_enabled = config%validate_interacting_covariance
      gamma_index = find_gamma_q_index(config%q_list)
      if (static_audit_enabled .and. gamma_index == 0) then
         error stop 'compact_dyson: dyson_static_audit requires Gamma; input should have been rejected during preflight'
      end if
      positive_q_index = 0
      negative_q_index = 0
      if (covariance_enabled) call find_covariance_q_indices(config%q_list, positive_q_index, negative_q_index)

      if (trim(config%channel) == 'chi_plus') then
         call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
         product => product_plus
      else
         call product_minus%initialize(response_space, radial_bases, lmto_product_channel_minus, .true.)
         product => product_minus
      end if
      rank_stable = .true.
      do i = 1, product%nsite
         do j = 0, product%response_lmax
            rank_stable = rank_stable .and. product%blocks(i, j)%rank_stable
         end do
      end do
      if (.not. rank_stable .or. product%product_dimension < 1) then
         error stop 'compact_dyson: retained compact product basis is empty or rank-unstable'
      end if
      if (response_space%response_lmax /= 2*ground_states(1)%pauli_lmax) then
         error stop 'compact_dyson: production requires the complete retained compact product basis'
      end if
      if (covariance_enabled) then
         if (trim(config%channel) == 'chi_plus') then
            call product_minus%initialize(response_space, radial_bases, lmto_product_channel_minus, .true.)
            opposite_product => product_minus
         else
            call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
            opposite_product => product_plus
         end if
         if (opposite_product%product_dimension /= product%product_dimension) then
            error stop 'compact_dyson: covariance requested but plus/minus product dimensions differ'
         end if
         do i = 1, opposite_product%nsite
            do j = 0, opposite_product%response_lmax
               if (.not. opposite_product%blocks(i, j)%rank_stable) then
                  error stop 'compact_dyson: covariance requested but opposite product basis is rank-unstable'
               end if
            end do
         end do
      end if
      ndim = product%product_dimension
      nfrequency = size(config%frequencies)
      nq = size(config%q_list, 2)

      allocate(magnetization(size(ground_states), response_space%npoint))
      call compute_accepted_pauli_magnetization(reciprocal_obj, lattice_obj%symbolic_atoms, lattice_obj%nbulk, magnetization)
      call prepare_direct_alsda_request(response_space, ground_states, magnetization, kxc_request)
      call evaluate_lr_alsda_kernel(kxc_request, kxc_result)
      allocate(interaction(ndim, ndim))
      call compact_project_local_operator(response_space, product, cmplx(kxc_result%pointwise_kernel, 0.0_rp, rp), interaction)
      static_eta = [0.01_rp, 0.005_rp]
      if (covariance_enabled) then
         allocate(opposite_interaction(ndim, ndim))
         call compact_project_local_operator(response_space, opposite_product, cmplx(kxc_result%pointwise_kernel, 0.0_rp, rp), &
            opposite_interaction)
         allocate(covariance_residual(nfrequency), covariance_loss_difference(nfrequency), &
            covariance_min_sv_difference(nfrequency), covariance_condition_difference(nfrequency))
         covariance_checked = .false.
      end if

      if (static_audit_enabled) allocate(static_frequency(1))
      magnetic_moment = 0.0_rp
      do i = 1, size(ground_states)
         magnetic_moment = magnetic_moment + ground_states(i)%integrated_moment_muB
      end do
      weight_sum = sum(left_state%k_weights)
      k_fingerprint = 0.0_rp
      do i = 1, left_state%nk
         k_fingerprint(1) = k_fingerprint(1) + left_state%k_weights(i)
         k_fingerprint(2:4) = k_fingerprint(2:4) + left_state%k_weights(i)*left_state%k_points(:, i)
         k_fingerprint(5) = k_fingerprint(5) + real(i, rp)*left_state%k_weights(i)
      end do
      state_file = trim(config%output_file)//'.state'
      if (rank == 0) then
         open(newunit=unit, file=trim(config%output_file), status='replace', action='write')
         write(unit, '(a)') '# TDDFT compact Dyson response'
         write(unit, '(a,a)') '# build_version = ', trim(tddft_build_version)
         write(unit, '(a)') '# backend = compact_dyson'
         write(unit, '(a,i0)') '# response_representation = orthonormal compact LMTO product representation; dimension=', ndim
         write(unit, '(a)') '# compact_dyson_equation = (I-chiKS*Kxc) chi = chiKS; LAPACK zgesv; no explicit inverse'
         write(unit, '(a)') '# loss_convention = L=-(chi-chi^dagger)/(2*i*pi), ordinary dagger in orthonormal compact space'
         write(unit, '(a)') '# state_source = accepted_kspace_scf_cache'
         write(unit, '(a)') '# accepted_state_cache_reused = T'
         write(unit, '(a)') '# reciprocal_rebuild_performed_for_tddft = F'
         write(unit, '(a,a)') '# structure_calctype = ', control_obj%calctype
         write(unit, '(a,es24.16)') '# structure_alat = ', lattice_obj%alat
         write(unit, '(a,i0)') '# structure_ntype = ', lattice_obj%ntype
         write(unit, '(a,i0)') '# structure_nrec = ', lattice_obj%nrec
         write(unit, '(a,a)') '# structure_symbols = ', trim(lattice_obj%symbolic_atoms(lattice_obj%nbulk + 1)%element%symbol)
         write(unit, '(a,3(i0,1x))') '# accepted_k_mesh = ', reciprocal_obj%nk_mesh
         write(unit, '(a,i0)') '# accepted_k_count = ', left_state%nk
         write(unit, '(a,i0)') '# product_dimension = ', ndim
         write(unit, '(a,i0)') '# response_angular_cutoff = ', response_space%response_lmax
         write(unit, '(a,a)') '# channel = ', trim(config%channel)
         write(unit, '(a,es24.16)') '# accepted_state_EF_Ry = ', left_state%fermi_level
         write(unit, '(a,es24.16)') '# scf_tdft_EF_difference_Ry = ', state_ef_diff
         write(unit, '(a,es24.16)') '# accepted_state_temperature_K = ', left_state%temperature
         write(unit, '(a,es24.16)') '# accepted_state_integrated_moment_muB = ', magnetic_moment
         write(unit, '(a,5(es24.16,1x))') '# k_fingerprint_checksums = ', k_fingerprint
         write(unit, '(a,es24.16)') '# k_weight_sum = ', weight_sum
         write(unit, '(a,es24.16)') '# state_consistency_mesh_max_abs = ', state_mesh_max
         write(unit, '(a,es24.16)') '# state_consistency_weight_max_abs = ', state_weight_max
         write(unit, '(a,es24.16)') '# state_consistency_eigenvalue_max_abs_Ry = ', state_eigen_max
         write(unit, '(a,es24.16)') '# state_consistency_occupation_max_abs = ', state_occ_max
         write(unit, '(a,es24.16)') '# state_consistency_projector_max_abs = ', state_projector_max
         write(unit, '(a,es24.16)') '# state_consistency_projector_frobenius = ', state_projector_frobenius
         write(unit, '(a,a)') '# state_artifact = ', trim(state_file)
         write(unit, '(a,a)') '# interaction_route = ', trim(config%interaction_route)
         write(unit, '(a)') '# interaction_provenance = KXC-01 direct ALSDA LR-03; compact projection U^H K_point U'
         write(unit, '(a,a)') '# magnetization_kind = ', trim(kxc_result%magnetization_kind)
         write(unit, '(a,a)') '# magnetization_source = ', trim(kxc_result%magnetization_source)
         write(unit, '(a)') '# BES_GCR_goldstone_correction = OFF'
         write(unit, '(a)') '# correction_status = no implicit Goldstone/BES/GCR correction'
         write(unit, '(a,es24.16)') '# physical_eta_Ry = ', config%eta
         write(unit, '(a,i0)') '# frequency_count = ', nfrequency
         do iw = 1, nfrequency
            write(unit, '(a,es24.16)') '# frequency_Ry = ', config%frequencies(iw)
         end do
         write(unit, '(a,i0)') '# q_count = ', nq
         write(unit, '(a)') '# q_convention = literal direct reciprocal coordinates'
         do iq = 1, nq
            write(unit, '(a,i0,3(1x,es24.16))') '# q_index = ', iq, config%q_list(:, iq)
         end do
         write(unit, '(a,l1)') '# dyson_static_audit = ', static_audit_enabled
         write(unit, '(a,l1)') '# validate_interacting_covariance = ', covariance_enabled
         if (static_audit_enabled) then
            write(unit, '(a)') '# static_denominator = validation diagnostic; fixed eta values are 0.01 and 0.005 Ry'
            write(unit, '(a)') '# static_denominator columns: eta_Ry min_singular_value max_singular_value condition_number min_magnitude_eigenvalue residual_F residual_relative residual_dInf'
         end if
         write(unit, '(a)') '# dynamic_metrics columns: q_index qx qy qz omega_Ry norm_chiKS norm_chi loss_trace_real loss_trace_imag min_sv max_sv condition residual_F residual_relative residual_dInf'
         if (config%write_full_matrix) then
            write(unit, '(a)') '# raw_matrix columns: q_index omega_Ry row column chiKS_real chiKS_imag chi_real chi_imag loss_real loss_imag'
         else
            write(unit, '(a)') '# output_storage_mode = streaming_summary'
            write(unit, '(a)') '# maximum_retained_q_batches = 1 (q-local bare, Dyson, and loss matrices)'
         end if
         if (config%write_full_matrix) then
            write(unit, '(a)') '# output_storage_mode = streaming_full_matrix'
            write(unit, '(a)') '# maximum_retained_q_batches = 1 (q-local bare, Dyson, and loss matrices)'
         end if
         if (covariance_enabled) then
            write(unit, '(a)') '# interacting_covariance = validation diagnostic; complete compact matrix transport'
            write(unit, '(a)') '# interacting_covariance columns: positive_q_index negative_q_index omega_Ry residual loss_difference min_sv_difference condition_difference'
         end if
      end if

      if (static_audit_enabled) then
         static_frequency(1) = 0.0_rp
         bare_request%q = config%q_list(:, gamma_index)
         bare_request%frequencies = static_frequency
         bare_request%channel = config%channel
         bare_request%product_basis => product
         bare_request%electronic_state => left_state
         bare_request%q_endpoint_state => endpoints(gamma_index)
         do i = 1, 2
            bare_request%eta = static_eta(i)
            call evaluate_lr_product_ks_susceptibility(bare_request, bare_result)
            dyson_request%response_space => response_space
            dyson_request%compact_orthonormal = .true.
            dyson_request%q = config%q_list(:, gamma_index)
            dyson_request%frequencies = static_frequency
            dyson_request%eta = static_eta(i)
            dyson_request%channel = config%channel
            dyson_request%ks_susceptibility = bare_result%susceptibility
            dyson_request%canonical_interaction = interaction
            dyson_request%interaction_route = tddft_driver_route_direct_alsda
            dyson_request%interaction_provenance = 'KXC-01 direct ALSDA LR-03'
            dyson_request%electronic_state_provenance = 'accepted k-space SCF state; SCF occupations and EF fixed'
            dyson_request%response_space_metadata = 'product_dimension='//trim(int2str(ndim))//'; orthonormal compact product space'
            call evaluate_tddft_dyson(dyson_request, dyson_result)
            static_min_sv(i) = dyson_result%denominator_min_singular_value(1)
            static_max_sv(i) = dyson_result%denominator_max_singular_value(1)
            static_condition(i) = dyson_result%denominator_condition_number(1)
            static_min_eigen(i) = dyson_result%denominator_min_magnitude_eigenvalue(1)
            static_residual(i) = dyson_result%dyson_residual_frobenius(1)
            static_residual_relative(i) = dyson_result%dyson_residual_relative(1)
            static_residual_infinity(i) = dyson_result%dyson_residual_infinity(1)
            static_status(i) = dyson_result%frequency_status(1)
            if (.not. dyson_result%solve_succeeded(1)) error stop 'compact_dyson static audit: Dyson solve failed'
            if (rank == 0) then
               write(unit, '(8(es24.16,1x))') static_eta(i), static_min_sv(i), static_max_sv(i), static_condition(i), &
                  static_min_eigen(i), static_residual(i), static_residual_relative(i), static_residual_infinity(i)
               write(unit, '(a,es24.16,1x,a)') '# static_solver_status eta=', static_eta(i), trim(static_status(i))
            end if
         end do
         if (rank == 0) write(unit, '(a)') '# static_rigid_vector_overlap = unavailable in the existing Dyson eigensolver API'
      end if

      do iq = 1, nq
         bare_request%q = config%q_list(:, iq)
         bare_request%frequencies = config%frequencies
         bare_request%eta = config%eta
         bare_request%channel = config%channel
         bare_request%product_basis => product
         bare_request%electronic_state => left_state
         bare_request%q_endpoint_state => endpoints(iq)
         call evaluate_lr_product_ks_susceptibility(bare_request, bare_result)
         finite_response = all(ieee_is_finite(real(bare_result%susceptibility, rp))) .and. &
            all(ieee_is_finite(aimag(bare_result%susceptibility)))
         if (.not. finite_response) error stop 'compact_dyson: bare compact response contains NaN or Inf'
         dyson_request%response_space => response_space
         dyson_request%compact_orthonormal = .true.
         dyson_request%q = config%q_list(:, iq)
         dyson_request%frequencies = config%frequencies
         dyson_request%eta = config%eta
         dyson_request%channel = config%channel
         dyson_request%ks_susceptibility = bare_result%susceptibility
         dyson_request%canonical_interaction = interaction
         dyson_request%interaction_route = tddft_driver_route_direct_alsda
         dyson_request%interaction_provenance = 'KXC-01 direct ALSDA LR-03'
         dyson_request%electronic_state_provenance = 'accepted k-space SCF state; SCF occupations and EF fixed'
         dyson_request%response_space_metadata = 'product_dimension='//trim(int2str(ndim))//'; orthonormal compact product space'
         call evaluate_tddft_dyson(dyson_request, dyson_result)
         if (any(.not. dyson_result%solve_succeeded)) error stop 'compact_dyson: Dyson denominator solve failed'
         finite_response = all(ieee_is_finite(real(dyson_result%enhanced_susceptibility, rp))) .and. &
            all(ieee_is_finite(aimag(dyson_result%enhanced_susceptibility))) .and. &
            all(ieee_is_finite(real(dyson_result%loss_matrix, rp))) .and. &
            all(ieee_is_finite(aimag(dyson_result%loss_matrix)))
         if (.not. finite_response) error stop 'compact_dyson: interacting compact response contains NaN or Inf'

         if (rank == 0) then
            do iw = 1, nfrequency
               trace_loss = 0.0_rp
               trace_loss_imag = 0.0_rp
               do i = 1, ndim
                  trace_loss = trace_loss + real(dyson_result%loss_matrix(i, i, iw), rp)
                  trace_loss_imag = trace_loss_imag + aimag(dyson_result%loss_matrix(i, i, iw))
               end do
               write(unit, '(i0,1x,4(es24.16,1x),10(es24.16,1x))') iq, config%q_list(:, iq), config%frequencies(iw), &
                  sqrt(sum(abs(bare_result%susceptibility(:, :, iw))**2)), &
                  sqrt(sum(abs(dyson_result%enhanced_susceptibility(:, :, iw))**2)), trace_loss, trace_loss_imag, &
                  dyson_result%denominator_min_singular_value(iw), dyson_result%denominator_max_singular_value(iw), &
                  dyson_result%denominator_condition_number(iw), dyson_result%dyson_residual_frobenius(iw), &
                  dyson_result%dyson_residual_relative(iw), dyson_result%dyson_residual_infinity(iw)
               write(unit, '(a,i0,1x,es24.16,1x,a)') '# dynamic_solver_status q_index=', iq, config%frequencies(iw), &
                  trim(dyson_result%frequency_status(iw))
               if (config%write_full_matrix) then
                  do j = 1, ndim
                     do i = 1, ndim
                        write(unit, '(i0,1x,es24.16,1x,2(i0,1x),6(es24.16,1x))') iq, config%frequencies(iw), i, j, &
                           real(bare_result%susceptibility(i, j, iw), rp), aimag(bare_result%susceptibility(i, j, iw)), &
                           real(dyson_result%enhanced_susceptibility(i, j, iw), rp), aimag(dyson_result%enhanced_susceptibility(i, j, iw)), &
                           real(dyson_result%loss_matrix(i, j, iw), rp), aimag(dyson_result%loss_matrix(i, j, iw))
                     end do
                  end do
               end if
            end do
         end if

         if (covariance_enabled .and. iq == positive_q_index) then
            opposite_bare_request%q = -config%q_list(:, positive_q_index)
            opposite_bare_request%frequencies = -config%frequencies
            opposite_bare_request%eta = config%eta
            if (trim(config%channel) == 'chi_plus') then
               opposite_bare_request%channel = 'chi_minus'
            else
               opposite_bare_request%channel = 'chi_plus'
            end if
            opposite_bare_request%product_basis => opposite_product
            opposite_bare_request%electronic_state => left_state
            opposite_bare_request%q_endpoint_state => endpoints(negative_q_index)
            call evaluate_lr_product_ks_susceptibility(opposite_bare_request, opposite_bare_result)
            opposite_dyson_request%response_space => response_space
            opposite_dyson_request%compact_orthonormal = .true.
            opposite_dyson_request%q = opposite_bare_request%q
            opposite_dyson_request%frequencies = opposite_bare_request%frequencies
            opposite_dyson_request%eta = config%eta
            opposite_dyson_request%channel = opposite_bare_request%channel
            opposite_dyson_request%ks_susceptibility = opposite_bare_result%susceptibility
            opposite_dyson_request%canonical_interaction = opposite_interaction
            opposite_dyson_request%interaction_route = tddft_driver_route_direct_alsda
            opposite_dyson_request%interaction_provenance = 'KXC-01 direct ALSDA LR-03'
            opposite_dyson_request%electronic_state_provenance = 'accepted k-space SCF state; SCF occupations and EF fixed'
            opposite_dyson_request%response_space_metadata = 'product_dimension='//trim(int2str(ndim))//'; orthonormal compact product space'
            call evaluate_tddft_dyson(opposite_dyson_request, opposite_dyson_result)
            if (any(.not. opposite_dyson_result%solve_succeeded)) then
               error stop 'compact_dyson covariance audit: opposite Dyson solve failed'
            end if
            covariance_checked = .true.
            do iw = 1, nfrequency
               if (trim(config%channel) == 'chi_plus') then
                  call compact_covariance_matrix_residual(product_plus, product_minus, dyson_result%enhanced_susceptibility(:, :, iw), &
                     opposite_dyson_result%enhanced_susceptibility(:, :, iw), covariance_residual(iw))
                  call compact_covariance_matrix_residual(product_plus, product_minus, dyson_result%loss_matrix(:, :, iw), &
                     -opposite_dyson_result%loss_matrix(:, :, iw), covariance_loss_difference(iw))
               else
                  call compact_covariance_matrix_residual(product_plus, product_minus, opposite_dyson_result%enhanced_susceptibility(:, :, iw), &
                     dyson_result%enhanced_susceptibility(:, :, iw), covariance_residual(iw))
                  call compact_covariance_matrix_residual(product_plus, product_minus, opposite_dyson_result%loss_matrix(:, :, iw), &
                     -dyson_result%loss_matrix(:, :, iw), covariance_loss_difference(iw))
               end if
               covariance_min_sv_difference(iw) = abs(dyson_result%denominator_min_singular_value(iw) - &
                  opposite_dyson_result%denominator_min_singular_value(iw))
               covariance_condition_difference(iw) = abs(dyson_result%denominator_condition_number(iw) - &
                  opposite_dyson_result%denominator_condition_number(iw))
               if (rank == 0) write(unit, '(2(i0,1x),5(es24.16,1x))') positive_q_index, negative_q_index, config%frequencies(iw), &
                  covariance_residual(iw), covariance_loss_difference(iw), covariance_min_sv_difference(iw), &
                  covariance_condition_difference(iw)
            end do
            if (.not. covariance_checked .or. maxval(covariance_residual) > 5.0e-8_rp .or. &
                maxval(covariance_loss_difference) > 5.0e-8_rp) then
               error stop 'compact_dyson covariance audit: interacting q/-q covariance failed'
            end if
         end if
      end do

      if (rank == 0) close(unit)
      write(*, '(a,i0,a,i0,a,i0)') 'TDDFT compact Dyson response: q_count=', nq, ' omega_count=', nfrequency, &
         ' product_dimension=', ndim
   end subroutine run_tddft_compact_dyson

   !> TDVK-06 static interaction diagnostics.  This backend consumes the
   !> accepted k-space SCF snapshot exactly once, evaluates the complete
   !> compact Gamma/omega=0 response, and stops after raw ALSDA and independent
   !> GSR diagnostics.  It never enters Dyson, loss, mode extraction, or any
   !> Goldstone repair/correction path.
   !> TDVK-02R2 validation-only handoff.  The accepted reciprocal snapshot is
   !> already naturally available here, so the live transition oracle and the
   !> compact bare response can be exercised without persisting eigenpairs or
   !> entering the full KXC/Dyson production lifecycle.
   subroutine run_tddft_product_bare_smoke(config, response_space, radial_bases, ground_states, left_state, endpoints, reciprocal_obj)
      type(tddft_runtime_config), intent(in) :: config
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lmto_product_response_basis), target :: product_plus, product_minus
      type(lr_product_ks_susceptibility_request) :: request
      type(lr_product_ks_susceptibility_result) :: result
      integer :: gamma_index, channel_kind
      real(rp) :: wall_start, wall_end, maximum_transition_error
      real(rp) :: accepted_moment
      logical :: finite_response
      integer :: isite

      gamma_index = find_gamma_q(config%q_list)
      if (gamma_index > size(endpoints)) error stop 'TDVK-02R2 product smoke: Gamma endpoint is unavailable'
      call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      call product_minus%initialize(response_space, radial_bases, lmto_product_channel_minus, .true.)
      if (product_plus%product_dimension /= product_plus%unpruned_dimension .or. &
          product_minus%product_dimension /= product_minus%unpruned_dimension) then
         error stop 'TDVK-02R2 product smoke: accepted Fe product basis is not full rank'
      end if

      maximum_transition_error = 0.0_rp
      do channel_kind = 1, 2
         if (channel_kind == 1) then
            call run_product_transition_oracle(response_space, radial_bases, product_plus, left_state, endpoints(gamma_index), &
               lmto_product_channel_plus, maximum_transition_error)
         else
            call run_product_transition_oracle(response_space, radial_bases, product_minus, left_state, endpoints(gamma_index), &
               lmto_product_channel_minus, maximum_transition_error)
         end if
      end do

      request%q = config%q_list(:, gamma_index)
      request%frequencies = config%frequencies
      request%eta = config%eta
      request%channel = config%channel
      request%electronic_state => left_state
      request%q_endpoint_state => endpoints(gamma_index)
      if (trim(config%channel) == 'chi_plus') then
         request%product_basis => product_plus
      else
         request%product_basis => product_minus
      end if
      call cpu_time(wall_start)
      call evaluate_lr_product_ks_susceptibility(request, result)
      call cpu_time(wall_end)
      finite_response = all(ieee_is_finite(real(result%susceptibility, rp))) .and. &
         all(ieee_is_finite(aimag(result%susceptibility)))
      if (.not. finite_response) error stop 'TDVK-02R2 product smoke: compact response contains NaN or Inf'

      accepted_moment = 0.0_rp
      do isite = 1, size(ground_states)
         accepted_moment = accepted_moment + ground_states(isite)%integrated_moment_muB
      end do
      call write_tddft_product_bare_smoke(config, result, reciprocal_obj, accepted_moment, maximum_transition_error, &
         wall_end - wall_start)
      write (*, '(a,i0,a,es12.4,a,es12.4)') 'TDVK-02R2 compact Fe smoke: Nprod=', result%product_dimension, &
         ' runtime_cpu_s=', wall_end - wall_start, ' max_transition_error=', maximum_transition_error
   end subroutine run_tddft_product_bare_smoke

   !> TDVK-04 accepted-Fe finite-q compact validation.  The Lehmann service is
   !> evaluated at every requested q, while the independent product-GF service
   !> is used only for two representative finite-q spots.  This routine owns
   !> diagnostics and serialization only; it does not introduce response
   !> physics or enter the point-space KXC/Dyson lifecycle.
   !> TDVK-05 compact Lehmann convergence campaign.  The caller supplies one
   !> accepted SCF state and its exact q endpoints; this routine varies only
   !> the physical response eta values and records complete compact-matrix
   !> diagnostics.  It deliberately stops before reciprocal GF, KXC,
   !> Goldstone, Dyson, loss, and mode interpretation.
   subroutine compact_covariance_residual(plus_product, minus_product, plus_result, minus_result, residual)
      type(lmto_product_response_basis), intent(in) :: plus_product, minus_product
      type(lr_product_ks_susceptibility_result), intent(in) :: plus_result, minus_result
      real(rp), intent(out) :: residual
      if (plus_product%product_dimension /= minus_product%product_dimension .or. &
          any(shape(plus_result%susceptibility) /= shape(minus_result%susceptibility))) then
         error stop 'TDVK-04 covariance: compact plus/minus dimensions differ'
      end if
      call compact_covariance_matrix_residual(plus_product, minus_product, plus_result%susceptibility(:, :, 1), &
         minus_result%susceptibility(:, :, 1), residual)
   end subroutine compact_covariance_residual

   !> Compare complete compact matrices under the established circular
   !> endpoint transport.  The first matrix is chi_plus(q,w), the second is
   !> chi_minus(-q,-w).  No direct element-wise conjugation or phase patch is
   !> substituted for this transport.
   subroutine compact_covariance_matrix_residual(plus_product, minus_product, plus_matrix, minus_matrix, residual)
      type(lmto_product_response_basis), intent(in) :: plus_product, minus_product
      complex(rp), intent(in) :: plus_matrix(:, :), minus_matrix(:, :)
      real(rp), intent(out) :: residual
      complex(rp), allocatable :: transport(:, :), mapped(:, :), block_transport(:, :)
      integer :: site, response_l, response_m, plus_first, minus_first, rank_plus, rank_minus
      real(rp) :: angular_sign

      if (plus_product%product_dimension /= minus_product%product_dimension .or. &
          any(shape(plus_matrix) /= [plus_product%product_dimension, plus_product%product_dimension]) .or. &
          any(shape(minus_matrix) /= [minus_product%product_dimension, minus_product%product_dimension])) then
         error stop 'compact covariance: compact plus/minus dimensions differ'
      end if
      allocate(transport(plus_product%product_dimension, minus_product%product_dimension))
      transport = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, plus_product%nsite
         do response_l = 0, plus_product%response_lmax
            rank_plus = plus_product%blocks(site, response_l)%rank
            rank_minus = minus_product%blocks(site, response_l)%rank
            allocate(block_transport(rank_plus, rank_minus))
            ! The weighted modes are the orthonormal compact radial bases.  The
            ! endpoint-adjoint relation transports minus(-q) coordinates from
            ! its -M block into the plus(q) M block without comparing SVD
            ! representatives directly.
            block_transport = matmul(conjg(transpose(plus_product%blocks(site, response_l)%weighted_modes)), &
               conjg(minus_product%blocks(site, response_l)%weighted_modes))
            do response_m = -response_l, response_l
               plus_first = plus_product%flat_index(site, response_l, response_m, 1)
               minus_first = minus_product%flat_index(site, response_l, -response_m, 1)
               angular_sign = merge(-1.0_rp, 1.0_rp, mod(abs(response_m), 2) == 1)
               transport(plus_first:plus_first + rank_plus - 1, minus_first:minus_first + rank_minus - 1) = &
                  angular_sign*block_transport
            end do
            deallocate(block_transport)
         end do
      end do
      allocate(mapped(plus_product%product_dimension, plus_product%product_dimension))
      mapped = matmul(transport, matmul(conjg(minus_matrix), conjg(transpose(transport))))
      residual = maxval(abs(plus_matrix - mapped))
      deallocate(mapped, transport)
   end subroutine compact_covariance_matrix_residual

   !> Project the established compact q/-q covariance transport into the
   !> plus-channel DRESP site observable. Direct conjugation of independently
   !> SVD-compressed plus/minus coordinates is not the covariance convention;
   !> circular angular and radial transport is applied first.

   subroutine finite_q_endpoint_metadata(left_state, endpoint, q, q_folded, endpoint_error, endpoint_first_k, endpoint_first, &
                                         endpoint_last, unique_count)
      type(lr_electronic_state), intent(in) :: left_state, endpoint
      real(rp), intent(in) :: q(3)
      real(rp), intent(out) :: q_folded(3), endpoint_error, endpoint_first_k(3), endpoint_first(3), endpoint_last(3)
      integer, intent(out) :: unique_count
      real(rp) :: expected(3)
      integer :: ik, iu

      if (left_state%nk /= endpoint%nk) error stop 'TDVK-04 metadata: endpoint k-point count differs'
      q_folded = fold_fractional_kpoint(q)
      endpoint_error = 0.0_rp
      endpoint_first_k = left_state%k_points(:, 1)
      endpoint_first = endpoint%k_points(:, 1)
      endpoint_last = endpoint%k_points(:, endpoint%nk)
      unique_count = 0
      do ik = 1, endpoint%nk
         expected = fold_fractional_kpoint(left_state%k_points(:, ik) + q)
         endpoint_error = max(endpoint_error, maxval(abs(expected - endpoint%k_points(:, ik))))
         if (ik == 1) endpoint_first = endpoint%k_points(:, ik)
         endpoint_last = endpoint%k_points(:, ik)
         iu = 1
         do while (iu < ik)
            if (all(endpoint%k_points(:, ik) == endpoint%k_points(:, iu))) exit
            iu = iu + 1
         end do
         if (iu == ik) unique_count = unique_count + 1
      end do
   end subroutine finite_q_endpoint_metadata

   integer function find_first_nonzero_q(q_list) result(index_nonzero)
      real(rp), intent(in) :: q_list(:, :)
      index_nonzero = find_first_nonzero_q_index(q_list)
      if (index_nonzero /= 0) return
      error stop 'TDVK-04 finite-q validation: no nonzero q was supplied'
   end function find_first_nonzero_q

   integer function find_first_nonzero_q_index(q_list) result(index_nonzero)
      real(rp), intent(in) :: q_list(:, :)
      integer :: iq

      index_nonzero = 0
      do iq = 1, size(q_list, 2)
         if (sum(abs(q_list(:, iq))) > 1.0e-12_rp) then
            index_nonzero = iq
            return
         end if
      end do
   end function find_first_nonzero_q_index

   subroutine validate_covariance_q_pair(q_list)
      real(rp), intent(in) :: q_list(:, :)
      integer :: positive_index, negative_index

      call find_covariance_q_indices(q_list, positive_index, negative_index)
   end subroutine validate_covariance_q_pair

   subroutine find_covariance_q_indices(q_list, positive_index, negative_index)
      real(rp), intent(in) :: q_list(:, :)
      integer, intent(out) :: positive_index, negative_index

      positive_index = find_first_nonzero_q_index(q_list)
      if (positive_index == 0) then
         error stop 'TDDFT input: validate_interacting_covariance requires a nonzero q and its exact -q partner'
      end if
      negative_index = find_matching_q(q_list, -q_list(:, positive_index), 2.0e-11_rp)
      if (negative_index == 0) then
         error stop 'TDDFT input: validate_interacting_covariance requires an exact +q/-q pair; rejected before SCF/response work'
      end if
   end subroutine find_covariance_q_indices

   integer function find_matching_q(q_list, target, tolerance) result(index_match)
      real(rp), intent(in) :: q_list(:, :), target(3), tolerance
      integer :: iq

      index_match = 0
      do iq = 1, size(q_list, 2)
         if (maxval(abs(q_list(:, iq) - target)) <= tolerance) then
            index_match = iq
            return
         end if
      end do
   end function find_matching_q

   integer function find_arbitrary_q(q_list, gamma_index, covariance_index, negative_index) result(index_arbitrary)
      real(rp), intent(in) :: q_list(:, :)
      integer, intent(in) :: gamma_index, covariance_index, negative_index
      integer :: iq

      index_arbitrary = 0
      do iq = 1, size(q_list, 2)
         if (iq /= gamma_index .and. iq /= covariance_index .and. iq /= negative_index .and. &
             sum(abs(q_list(:, iq))) > 1.0e-12_rp) then
            index_arbitrary = iq
            return
         end if
      end do
   end function find_arbitrary_q

   pure function fold_fractional_kpoint(k_point) result(folded)
      real(rp), intent(in) :: k_point(3)
      real(rp) :: folded(3)

      folded = k_point - floor(k_point + 0.5_rp)
   end function fold_fractional_kpoint

   !> TDVK-02R3 accepted-Fe compact reciprocal-GF smoke.  This is deliberately
   !> separate from the legacy point-grid reciprocal-GF backend and is not a
   !> production interaction/Dyson route.
   !> TDVK-03 accepted-state closure audit.  The SCF handoff above has already
   !> prepared one immutable reciprocal state and all exact folded endpoints.
   !> This routine deliberately evaluates the compact Lehmann reference once,
   !> then varies only GF integration controls while keeping those snapshots
   !> and the complete six-branch product basis fixed.
   subroutine run_product_transition_oracle(response_space, radial_bases, product, left_state, right_state, channel_kind, maximum_error)
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(lmto_product_response_basis), intent(in) :: product
      type(lr_electronic_state), intent(in) :: left_state, right_state
      integer, intent(in) :: channel_kind
      real(rp), intent(inout) :: maximum_error
      type(pauli_vertex_capabilities) :: capabilities
      type(pauli_endpoint_state) :: left_band, right_band
      complex(rp), allocatable :: point_transition(:), product_coordinates(:), reference_coordinates(:)
      complex(rp) :: operator_matrix(2, 2)
      integer :: selected_left(3), selected_right(3), ik, itransition
      real(rp) :: error

      call select_live_transitions(product, left_state, right_state, selected_left, selected_right, ik)
      allocate(point_transition(response_space%ndim), product_coordinates(product%product_dimension), &
         reference_coordinates(product%product_dimension))
      capabilities = pauli_vertex_capabilities()
      if (channel_kind == lmto_product_channel_plus) then
         operator_matrix = pauli_sigma_plus_matrix()
      else
         operator_matrix = pauli_sigma_minus_matrix()
      end if
      do itransition = 1, size(selected_left)
         call left_band%initialize(left_state%eigenvalues(selected_left(itransition), ik), &
            left_state%eigenvectors(:, selected_left(itransition), ik))
         call right_band%initialize(right_state%eigenvalues(selected_right(itransition), ik), &
            right_state%eigenvectors(:, selected_right(itransition), ik))
         call evaluate_pauli_transition_vertex(response_space, radial_bases, left_band, right_band, operator_matrix, &
            capabilities, point_transition)
         call product%transition_coordinates(left_band, right_band, product_coordinates)
         call project_point_transition(response_space, product, point_transition, reference_coordinates)
         error = sqrt(sum(abs(product_coordinates - reference_coordinates)**2))/ &
            max(sqrt(sum(abs(reference_coordinates)**2)), epsilon(1.0_rp))
         maximum_error = max(maximum_error, error)
         write (*, '(a,a,a,i0,a,es12.4)') 'TDVK-02R2 live transition channel=', &
            product_channel_name(channel_kind), ' index=', itransition, ' residual=', error
         if (error >= 1.0e-10_rp .or. .not. ieee_is_finite(error)) then
            error stop 'TDVK-02R2 live transition oracle failed'
         end if
      end do
      deallocate(point_transition, product_coordinates, reference_coordinates)
   end subroutine run_product_transition_oracle

   subroutine select_live_transitions(product, left_state, right_state, selected_left, selected_right, selected_k)
      type(lmto_product_response_basis), intent(in) :: product
      type(lr_electronic_state), intent(in) :: left_state, right_state
      integer, intent(out) :: selected_left(:), selected_right(:), selected_k
      type(pauli_endpoint_state) :: left_band, right_band
      complex(rp), allocatable :: coordinates(:)
      real(rp) :: score, norm, nearest_score, deepest_energy
      real(rp), parameter :: occupation_tolerance = 1.0e-8_rp, coordinate_tolerance = 1.0e-12_rp
      integer :: ik, ib, jb, nearest_left, nearest_right, deep_left, deep_right, first_left, first_right
      integer :: additional_left, additional_right

      if (size(selected_left) < 3 .or. size(selected_right) /= size(selected_left)) then
         error stop 'TDVK-02R2 transition selector: output shape mismatch'
      end if
      selected_k = 1
      do ik = 2, left_state%nk
         if (sum(left_state%k_points(:, ik)**2) < sum(left_state%k_points(:, selected_k)**2)) selected_k = ik
      end do
      allocate(coordinates(product%product_dimension))
      nearest_score = huge(1.0_rp)
      deepest_energy = huge(1.0_rp)
      nearest_left = 0
      nearest_right = 0
      deep_left = 0
      deep_right = 0
      first_left = 0
      first_right = 0
      do ib = 1, left_state%nbands
         do jb = 1, right_state%nbands
            if (left_state%occupations(ib, selected_k) <= 0.5_rp + occupation_tolerance .or. &
                right_state%occupations(jb, selected_k) >= 0.5_rp - occupation_tolerance .or. &
                left_state%occupations(ib, selected_k) - right_state%occupations(jb, selected_k) <= occupation_tolerance) cycle
            call left_band%initialize(left_state%eigenvalues(ib, selected_k), left_state%eigenvectors(:, ib, selected_k))
            call right_band%initialize(right_state%eigenvalues(jb, selected_k), right_state%eigenvectors(:, jb, selected_k))
            call product%transition_coordinates(left_band, right_band, coordinates)
            norm = sqrt(sum(abs(coordinates)**2))
            if (norm <= coordinate_tolerance) cycle
            if (first_left == 0) then
               first_left = ib
               first_right = jb
            end if
            score = abs(left_state%eigenvalues(ib, selected_k) - left_state%fermi_level) + &
               abs(right_state%eigenvalues(jb, selected_k) - right_state%fermi_level)
            if (score < nearest_score) then
               nearest_score = score
               nearest_left = ib
               nearest_right = jb
            end if
            if (left_state%eigenvalues(ib, selected_k) < left_state%fermi_level - 0.20_rp .and. &
                left_state%eigenvalues(ib, selected_k) < deepest_energy) then
               deepest_energy = left_state%eigenvalues(ib, selected_k)
               deep_left = ib
               deep_right = jb
            end if
         end do
      end do
      if (nearest_left == 0 .or. first_left == 0) error stop 'TDVK-02R2 transition selector: no nonzero Fe Gamma transition'
      if (deep_left == 0) then
         deep_left = first_left
         deep_right = first_right
      end if
      additional_left = 0
      additional_right = 0
      do ib = 1, left_state%nbands
         do jb = 1, right_state%nbands
            if (ib == nearest_left .and. jb == nearest_right) cycle
            if (ib == deep_left .and. jb == deep_right) cycle
            if (left_state%occupations(ib, selected_k) <= 0.5_rp + occupation_tolerance .or. &
                right_state%occupations(jb, selected_k) >= 0.5_rp - occupation_tolerance .or. &
                left_state%occupations(ib, selected_k) - right_state%occupations(jb, selected_k) <= occupation_tolerance) cycle
            call left_band%initialize(left_state%eigenvalues(ib, selected_k), left_state%eigenvectors(:, ib, selected_k))
            call right_band%initialize(right_state%eigenvalues(jb, selected_k), right_state%eigenvectors(:, jb, selected_k))
            call product%transition_coordinates(left_band, right_band, coordinates)
            if (sqrt(sum(abs(coordinates)**2)) > coordinate_tolerance) then
               additional_left = ib
               additional_right = jb
               exit
            end if
         end do
         if (additional_left /= 0) exit
      end do
      if (additional_left == 0) error stop 'TDVK-02R2 transition selector: fewer than three distinct transitions'
      selected_left = [nearest_left, deep_left, additional_left]
      selected_right = [nearest_right, deep_right, additional_right]
      deallocate(coordinates)
   end subroutine select_live_transitions

   subroutine project_point_transition(response_space, product, point_transition, product_coordinates)
      type(response_space_layout), intent(in) :: response_space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: point_transition(:)
      complex(rp), intent(out) :: product_coordinates(:)
      type(response_super_index) :: item
      integer :: product_flat, site, response_l, response_m, product_mode, ir, point_flat

      product_coordinates = cmplx(0.0_rp, 0.0_rp, rp)
      do product_flat = 1, product%product_dimension
         call product%unflatten_index(product_flat, site, response_l, response_m, product_mode)
         do ir = 1, response_space%npoint
            item = response_super_index(site, response_l, response_m, ir, 1)
            call response_flatten_superindex(item, response_space%nsite, response_space%response_lmax, &
               response_space%npoint, response_space%nchannel, point_flat)
            product_coordinates(product_flat) = product_coordinates(product_flat) + &
               conjg(product%blocks(site, response_l)%weighted_modes(ir, product_mode))* &
               sqrt(response_space%radial_weights(ir))*point_transition(point_flat)
         end do
      end do
   end subroutine project_point_transition

   subroutine write_tddft_product_bare_smoke(config, result, reciprocal_obj, accepted_moment, maximum_transition_error, runtime)
      type(tddft_runtime_config), intent(in) :: config
      type(lr_product_ks_susceptibility_result), intent(in) :: result
      type(reciprocal), intent(in) :: reciprocal_obj
      real(rp), intent(in) :: accepted_moment, maximum_transition_error, runtime
      integer :: unit, ifrequency, i
      real(rp) :: frobenius_norm, maximum_element, memory_bytes
      complex(rp) :: trace

      memory_bytes = 16.0_rp*real(result%product_dimension*result%product_dimension*size(result%frequencies), rp)
      open(newunit=unit, file=trim(config%output_file), status='replace', action='write')
      write(unit, '(a)') '# TDVK-02R2 compact Lehmann bare response; no KXC or Dyson invoked'
      write(unit, '(a,a)') '# build_version = ', trim(tddft_build_version)
      write(unit, '(a,i0)') '# product_dimension = ', result%product_dimension
      write(unit, '(a,es24.16)') '# response_matrix_memory_bytes = ', memory_bytes
      write(unit, '(a,es24.16)') '# response_matrix_memory_MiB = ', memory_bytes/(1024.0_rp**2)
      write(unit, '(a,es24.16)') '# runtime_cpu_seconds = ', runtime
      write(unit, '(a,es24.16)') '# accepted_moment_muB = ', accepted_moment
      write(unit, '(a,es24.16)') '# EF_Ry = ', reciprocal_obj%fermi_level
      write(unit, '(a,es24.16)') '# eta_Ry = ', config%eta
      write(unit, '(a,es24.16)') '# maximum_live_transition_residual = ', maximum_transition_error
      write(unit, '(a,l1)') '# finite_response = ', all(ieee_is_finite(real(result%susceptibility, rp))) .and. &
         all(ieee_is_finite(aimag(result%susceptibility)))
      write(unit, '(a)') '# columns: omega_Ry frobenius_norm max_abs_element trace_real trace_imag'
      do ifrequency = 1, size(result%frequencies)
         frobenius_norm = sqrt(sum(abs(result%susceptibility(:, :, ifrequency))**2))
         maximum_element = maxval(abs(result%susceptibility(:, :, ifrequency)))
         trace = cmplx(0.0_rp, 0.0_rp, rp)
         do i = 1, result%product_dimension
            trace = trace + result%susceptibility(i, i, ifrequency)
         end do
         write(unit, '(5(es24.16,1x))') result%frequencies(ifrequency), frobenius_norm, maximum_element, &
            real(trace, rp), aimag(trace)
      end do
      close(unit)
   end subroutine write_tddft_product_bare_smoke

   character(len=5) function product_channel_name(channel_kind) result(value)
      integer, intent(in) :: channel_kind

      if (channel_kind == lmto_product_channel_plus) then
         value = 'plus '
      else
         value = 'minus'
      end if
   end function product_channel_name

   !> Assemble the single production direct-ALSDA request.  The caller must
   !> have obtained the array from compute_accepted_pauli_magnetization; the
   !> provenance fields are deliberately assigned here so an SR density cannot
   !> be relabelled by an individual driver.
   subroutine prepare_direct_alsda_request(response_space, ground_states, pauli_magnetization, request)
      type(response_space_layout), target, intent(in) :: response_space
      type(radial_ground_state), target, intent(in) :: ground_states(:)
      real(rp), intent(in) :: pauli_magnetization(:, :)
      type(lr_alsda_kernel_request), intent(out) :: request
      real(rp) :: sr_difference, scale
      integer :: isite

      sr_difference = 0.0_rp
      scale = 0.0_rp
      do isite = 1, response_space%nsite
         sr_difference = max(sr_difference, maxval(abs(pauli_magnetization(isite, :) - &
            (ground_states(isite)%n_up - ground_states(isite)%n_down))))
         scale = max(scale, maxval(abs(pauli_magnetization(isite, :))))
      end do
      if (sr_difference <= 1.0e-12_rp*max(1.0_rp, scale)) then
         error stop 'TDDFT production driver: LR-01 SR n_up-n_down cannot masquerade as accepted Pauli magnetization'
      end if
      request%response_space => response_space
      request%ground_states => ground_states
      allocate(request%pauli_magnetization(size(pauli_magnetization, 1), size(pauli_magnetization, 2)))
      request%pauli_magnetization = pauli_magnetization
      request%magnetization_label = 'pauli_projected'
      request%magnetization_kind = lr_kxc_magnetization_kind_pauli_accepted
      request%magnetization_source = lr_kxc_magnetization_source_pauli_accepted
      request%production_contract = .true.
   end subroutine prepare_direct_alsda_request

   module subroutine evaluate_linear_response_sweep(config, response_space, radial_bases, ground_states, left_state, endpoints, result, &
                                                   native_provider, native_pairs, native_site_positions, accepted_pauli_magnetization, &
                                                   accepted_pauli_magnetization_source)
      type(linear_response_config), intent(in) :: config
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(radial_ground_state), target, intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(tddft_production_result), intent(out) :: result
      class(lr_rs_gf_provider), target, intent(inout), optional :: native_provider
      type(lr_rs_gf_pair), intent(in), optional :: native_pairs(:)
      real(rp), target, intent(in), optional :: native_site_positions(:, :)
      real(rp), allocatable, intent(in), optional :: accepted_pauli_magnetization(:, :)
      character(len=*), intent(in), optional :: accepted_pauli_magnetization_source
      type(tddft_runtime_config) :: runtime

      call derive_runtime_config(config, runtime)
      call evaluate_tddft_runtime_sweep(runtime, response_space, radial_bases, ground_states, left_state, endpoints, result, &
         native_provider, native_pairs, native_site_positions, accepted_pauli_magnetization, accepted_pauli_magnetization_source)
   end subroutine evaluate_linear_response_sweep

   !> Evaluate an already prepared request batch. This is the reproducibility
   !> seam: tests and validation can compare it directly with service calls,
   !> without reconstructing SCF state.
   subroutine evaluate_tddft_runtime_sweep(config, response_space, radial_bases, ground_states, left_state, endpoints, result, &
                                              native_provider, native_pairs, native_site_positions, accepted_pauli_magnetization, &
                                              accepted_pauli_magnetization_source)
      type(tddft_runtime_config), intent(in) :: config
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(radial_ground_state), target, intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(tddft_production_result), intent(out) :: result
      class(lr_rs_gf_provider), target, intent(inout), optional :: native_provider
      type(lr_rs_gf_pair), intent(in), optional :: native_pairs(:)
      real(rp), target, intent(in), optional :: native_site_positions(:, :)
      real(rp), allocatable, intent(in), optional :: accepted_pauli_magnetization(:, :)
      character(len=*), intent(in), optional :: accepted_pauli_magnetization_source

      type(lr_alsda_kernel_request) :: kxc_request
      type(lr_alsda_kernel_result) :: kxc_result
      type(lr_goldstone_sumrule_request) :: gsr_request
      type(lr_goldstone_sumrule_result) :: gsr_result
      type(lr_ks_susceptibility_result) :: bare_result, static_result
      type(tddft_dyson_request) :: dyson_request
      type(tddft_dyson_result) :: dyson_result
      complex(rp), allocatable :: interaction(:, :)
      real(rp), allocatable :: gsr_magnetization(:, :), static_frequency(:)
      integer :: iq, i0, ndim

      call validate_prepared_sweep(config, response_space, radial_bases, ground_states, left_state, endpoints, &
         native_provider, native_pairs)
      ndim = response_space%ndim
      result%nq = size(config%q_list, 2)
      result%nfrequency = size(config%frequencies)
      result%ndim = ndim
      result%response_lmax = response_space%response_lmax
      result%interaction_route = trim(config%interaction_route)
      result%backend = trim(config%backend)
      result%bare_response_provenance = ''
      result%q_list = config%q_list
      result%frequencies = config%frequencies
      allocate(result%ks_susceptibility(ndim, ndim, result%nfrequency, result%nq), &
               result%enhanced_susceptibility(ndim, ndim, result%nfrequency, result%nq), &
               result%loss_matrix(ndim, ndim, result%nfrequency, result%nq))
      result%ks_susceptibility = cmplx(0.0_rp, 0.0_rp, rp)
      result%enhanced_susceptibility = cmplx(0.0_rp, 0.0_rp, rp)
      result%loss_matrix = cmplx(0.0_rp, 0.0_rp, rp)

      select case (trim(config%interaction_route))
      case (tddft_driver_route_direct_alsda)
         if (.not. present(accepted_pauli_magnetization) .or. .not. allocated(accepted_pauli_magnetization)) then
            error stop 'TDDFT production driver: direct_alsda requires certified accepted Pauli magnetization'
         end if
         if (.not. present(accepted_pauli_magnetization_source) .or. &
             trim(accepted_pauli_magnetization_source) /= lr_kxc_magnetization_source_pauli_accepted) then
            error stop 'TDDFT production driver: direct_alsda Pauli magnetization provenance is not certified'
         end if
         call prepare_direct_alsda_request(response_space, ground_states, accepted_pauli_magnetization, kxc_request)
         call evaluate_lr_alsda_kernel(kxc_request, kxc_result)
         allocate(interaction(ndim, ndim))
         interaction = kxc_result%canonical_operator
         result%interaction_provenance = 'KXC-01 direct ALSDA LR-03'
         result%magnetization_kind = kxc_result%magnetization_kind
         result%magnetization_source = kxc_result%magnetization_source
         result%goldstone_correction_status = 'disabled'
      case (tddft_driver_route_goldstone_sumrule)
         allocate(gsr_magnetization(size(ground_states), response_space%npoint))
         do iq = 1, size(ground_states)
            gsr_magnetization(iq, :) = ground_states(iq)%n_up - ground_states(iq)%n_down
         end do
         i0 = find_gamma_q(config%q_list)
         allocate(static_frequency(1))
         static_frequency(1) = 0.0_rp
         call evaluate_bare_response(config, response_space, radial_bases, left_state, endpoints(i0), &
                                     config%q_list(:, i0), static_frequency, static_result, native_provider, native_pairs, &
                                     native_site_positions)
         gsr_request%response_space => response_space
         gsr_request%static_susceptibility = static_result%susceptibility(:, :, 1)
         gsr_request%magnetization = gsr_magnetization
         call evaluate_lr_goldstone_sumrule(gsr_request, gsr_result)
         if (gsr_result%blocked) error stop 'TDDFT production driver: Goldstone sum-rule interaction is blocked by its service contract'
         allocate(interaction(ndim, ndim))
         interaction = gsr_result%canonical_interaction
         result%interaction_provenance = 'GSR-01 independent Goldstone sum-rule interaction'
         result%magnetization_kind = 'SR_LR01'
         result%magnetization_source = 'accepted scalar-relativistic LR-01 radial density'
         result%goldstone_correction_status = 'not selected'
      case default
         error stop 'TDDFT production driver: interaction route was not validated'
      end select

      do iq = 1, size(config%q_list, 2)
         if (trim(config%backend) == tddft_driver_backend_native_rsgf) then
            call g_logger%info('TDRUN-02 native RSGF sweep: ndim='//int2str(ndim)//' radial_points='// &
               int2str(response_space%npoint)//' integration_points='//int2str(config%gf_integration_points), __FILE__, __LINE__)
         end if
         call evaluate_bare_response(config, response_space, radial_bases, left_state, endpoints(iq), &
                                     config%q_list(:, iq), config%frequencies, bare_result, native_provider, native_pairs, &
                                     native_site_positions)
         result%bare_response_provenance = bare_result%response_space_metadata
         dyson_request%response_space => response_space
         dyson_request%q = config%q_list(:, iq)
         dyson_request%frequencies = config%frequencies
         dyson_request%eta = config%eta
         dyson_request%channel = config%channel
         dyson_request%ks_susceptibility = bare_result%susceptibility
         dyson_request%canonical_interaction = interaction
         dyson_request%interaction_route = trim(config%interaction_route)
         dyson_request%interaction_provenance = result%interaction_provenance
         dyson_request%electronic_state_provenance = 'accepted reciprocal eigenpair snapshot; occupations and EF fixed at SCF handoff'
         dyson_request%response_space_metadata = bare_result%response_space_metadata
         call evaluate_tddft_dyson(dyson_request, dyson_result)
         result%ks_susceptibility(:, :, :, iq) = dyson_result%ks_susceptibility
         result%enhanced_susceptibility(:, :, :, iq) = dyson_result%enhanced_susceptibility
         result%loss_matrix(:, :, :, iq) = dyson_result%loss_matrix
         result%status = trim(dyson_result%status)
      end do
      result%initialized = .true.
   end subroutine evaluate_tddft_runtime_sweep

   subroutine validate_prepared_sweep(config, response_space, radial_bases, ground_states, left_state, endpoints, &
                                      native_provider, native_pairs)
      type(tddft_runtime_config), intent(in) :: config
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      class(lr_rs_gf_provider), intent(in), optional :: native_provider
      type(lr_rs_gf_pair), intent(in), optional :: native_pairs(:)
      integer :: iq

      call validate_tddft_runtime_config(config)
      if (.not. config%enabled) error stop 'TDDFT production driver: disabled configuration cannot run a sweep'
      if (size(ground_states) /= response_space%nsite .or. size(radial_bases) /= response_space%nsite) then
         error stop 'TDDFT production driver: response-site state shape mismatch'
      end if
      if (size(endpoints) /= size(config%q_list, 2)) error stop 'TDDFT production driver: q endpoint count mismatch'
      call left_state%validate('TDDFT production driver:left_state')
      do iq = 1, size(endpoints)
         call endpoints(iq)%validate('TDDFT production driver:q_endpoint_state')
      end do
      if (trim(config%backend) == tddft_driver_backend_native_rsgf) then
         if (.not. present(native_provider) .or. .not. present(native_pairs)) then
            error stop 'TDDFT production driver: native_rsgf requires a registered provider and complete pair set'
         end if
         if (size(native_pairs) < 1) error stop 'TDDFT production driver: native_rsgf pair set is empty'
         if (native_provider%nbasis_per_site < 1) then
            error stop 'TDDFT production driver: native_rsgf provider is not initialized'
         end if
      end if
   end subroutine validate_prepared_sweep

   subroutine evaluate_bare_response(config, response_space, radial_bases, left_state, endpoint, q, frequencies, result, &
                                     native_provider, native_pairs, native_site_positions)
      type(tddft_runtime_config), intent(in) :: config
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoint
      real(rp), intent(in) :: q(3), frequencies(:)
      type(lr_ks_susceptibility_result), intent(out) :: result
      class(lr_rs_gf_provider), target, intent(inout), optional :: native_provider
      type(lr_rs_gf_pair), intent(in), optional :: native_pairs(:)
      real(rp), target, intent(in), optional :: native_site_positions(:, :)
      type(lr_rs_gf_susceptibility_request) :: native_request

      if (trim(config%backend) == tddft_driver_backend_lehmann .or. trim(config%backend) == 'spectral') then
         call evaluate_lehmann_backend(response_space, radial_bases, left_state, endpoint, q, frequencies, config%eta, &
            config%channel, result)
      else if (trim(config%backend) == tddft_driver_backend_native_rsgf) then
         if (.not. present(native_provider) .or. .not. present(native_pairs)) then
            error stop 'TDDFT production driver: native_rsgf request lacks its registered provider/pair set'
         end if
         native_request%q = q
         native_request%frequencies = frequencies
         native_request%eta = config%eta
         native_request%channel = config%channel
         native_request%integration_points = config%gf_integration_points
         native_request%integration_eta = config%gf_integration_eta
         native_request%energy_margin = config%gf_energy_margin
         native_request%energy_min = minval(left_state%eigenvalues)
         native_request%energy_max = maxval(left_state%eigenvalues)
         native_request%fermi_level = left_state%fermi_level
         native_request%temperature = left_state%temperature
         native_request%response_space => response_space
         native_request%radial_bases => radial_bases
         native_request%provider => native_provider
         native_request%pairs = native_pairs
         if (present(native_site_positions)) native_request%site_positions => native_site_positions
         call evaluate_lr_rs_gf_susceptibility(native_request, result)
      else
         error stop 'TDDFT production driver: backend was not validated'
      end if
   end subroutine evaluate_bare_response

   subroutine evaluate_lehmann_backend(response_space, radial_bases, left_state, endpoint, q, frequencies, eta, channel, result)
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(lr_electronic_state), target, intent(in) :: left_state, endpoint
      real(rp), intent(in) :: q(3), frequencies(:), eta
      character(len=*), intent(in) :: channel
      type(lr_ks_susceptibility_result), intent(out) :: result
      type(lr_ks_susceptibility_request) :: request

      request%q = q
      request%frequencies = frequencies
      request%eta = eta
      request%channel = channel
      request%response_space => response_space
      request%radial_bases => radial_bases
      request%electronic_state => left_state
      request%q_endpoint_state => endpoint
      call evaluate_lr_ks_susceptibility(request, result)
   end subroutine evaluate_lehmann_backend

   subroutine radial_basis_from_snapshot(state, basis, complete_sr)
      type(radial_ground_state), intent(in) :: state
      type(lmto_radial_basis), intent(out) :: basis
      logical, intent(in), optional :: complete_sr
      logical :: use_complete_sr

      use_complete_sr = .false.
      if (present(complete_sr)) use_complete_sr = complete_sr

      if (use_complete_sr .and. (.not. allocated(state%pauli_small) .or. .not. state%sr_basis_valid .or. &
          .not. allocated(state%sr_large_dot) .or. &
          .not. allocated(state%sr_small_dot) .or. .not. allocated(state%sr_large_ddot) .or. &
          .not. allocated(state%sr_small_ddot) .or. .not. allocated(state%sr_gfac) .or. &
          .not. allocated(state%sr_tmc) .or. .not. allocated(state%sr_potential) .or. &
          .not. allocated(state%sr_enu))) then
         error stop 'TDDFT production driver: accepted radial snapshot lacks complete SR basis arrays'
      end if
      call basis%initialize(size(state%r), state%pauli_lmax, 2)
      basis%rofi = state%r
      basis%mesh_a = state%a
      basis%mesh_b = state%b
      basis%phi_large = state%pauli_large
      basis%phidot_large = state%pauli_large_dot
      basis%enu_radial = state%pauli_enu
      basis%enu_work = state%pauli_enu
      if (use_complete_sr) then
         basis%phidot_small = state%sr_small_dot
         basis%phiddot_large = state%sr_large_ddot
         basis%phiddot_small = state%sr_small_ddot
         basis%gfac = state%sr_gfac
         basis%tmc = state%sr_tmc
         basis%potential = state%sr_potential
         basis%enu_radial = state%sr_enu
         basis%enu_work = state%sr_enu
      end if
      basis%channel_present = .true.
      if (use_complete_sr) basis%phi_small = state%pauli_small
   end subroutine radial_basis_from_snapshot

   pure function effective_hamiltonian_order(reciprocal_obj, hamiltonian_obj) result(order)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj
      character(len=16) :: order

      order = trim(reciprocal_obj%kspace_ham_order)
      if (order == 'auto') then
         if (hamiltonian_obj%hoh) then
            order = 'second'
         else
            order = 'first'
         end if
      end if
   end function effective_hamiltonian_order

   integer function find_gamma_q(q_list) result(index_gamma)
      real(rp), intent(in) :: q_list(:, :)
      index_gamma = find_gamma_q_index(q_list)
      if (index_gamma /= 0) return
      error stop 'TDDFT production driver: gamma q was not supplied for the requested static route'
   end function find_gamma_q

   integer function find_gamma_q_index(q_list) result(index_gamma)
      real(rp), intent(in) :: q_list(:, :)
      integer :: iq

      index_gamma = 0
      do iq = 1, size(q_list, 2)
         if (sum(abs(q_list(:, iq))) <= 1.0e-12_rp) then
            index_gamma = iq
            return
         end if
      end do
   end function find_gamma_q_index

   subroutine write_tddft_production_output(config, result, control_obj, lattice_obj, ground_states, response_space, reciprocal_obj, &
                                             direct_handoff)
      type(tddft_runtime_config), intent(in) :: config
      type(tddft_production_result), intent(in) :: result
      type(control), intent(in) :: control_obj
      type(lattice), intent(in) :: lattice_obj
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(response_space_layout), intent(in) :: response_space
      type(reciprocal), intent(in) :: reciprocal_obj
      logical, intent(in) :: direct_handoff
      integer :: unit, iq, iw, i, j
      real(rp) :: magnetic_moment
      complex(rp) :: loss_trace, ks_trace, enhanced_trace

      if (rank /= 0) return
      magnetic_moment = 0.0_rp
      do i = 1, size(ground_states)
         magnetic_moment = magnetic_moment + ground_states(i)%integrated_moment_muB
      end do
      open(newunit=unit, file=trim(config%output_file), status='replace', action='write')
      write(unit, '(a)') '# TDDFT production response; matrices are LR-04 canonical service outputs'
      write(unit, '(a,a)') '# git_version = ', trim(tddft_build_version)
      write(unit, '(a,a)') '# structure_calctype = ', control_obj%calctype
      write(unit, '(a,es24.16)') '# structure_alat = ', lattice_obj%alat
      write(unit, '(a,i0)') '# structure_ntype = ', lattice_obj%ntype
      write(unit, '(a,i0)') '# structure_nrec = ', lattice_obj%nrec
      write(unit, '(a,a)') '# structure_symbols = ', trim(lattice_obj%symbolic_atoms(lattice_obj%nbulk + 1)%element%symbol)
      write(unit, '(a,a)') '# xc_identity = ', trim(ground_states(1)%xc_provenance%functional_name)
      write(unit, '(a,a)') '# xc_backend = ', trim(ground_states(1)%xc_provenance%backend_name)
      write(unit, '(a,i0)') '# xc_txc = ', ground_states(1)%xc_provenance%txc
      write(unit, '(a,es24.16)') '# magnetic_moment_muB = ', magnetic_moment
      write(unit, '(a,es24.16)') '# fermi_level_Ry = ', reciprocal_obj%fermi_level
      write(unit, '(a,es24.16)') '# temperature_K = ', reciprocal_obj%temperature
      write(unit, '(a)') '# response_capability = collinear,no_soc,ham_only,orthogonal,sp/spd,no_extra_operator'
      write(unit, '(a,a)') '# q = ', 'one row block per requested q; values are reduced coordinates'
      write(unit, '(a,a)') '# omega_Ry = ', 'one row per requested frequency'
      write(unit, '(a,es24.16)') '# eta_Ry = ', config%eta
      write(unit, '(a,i0)') '# response_angular_cutoff = ', result%response_lmax
      if (direct_handoff) then
         write(unit, '(a)') '# state_source = accepted_kspace_scf_cache'
         write(unit, '(a,l1)') '# direct_accepted_state_handoff = ', .true.
      else
         write(unit, '(a)') '# state_source = diagnostic_frozen_post_scf_rebuild'
         write(unit, '(a,l1)') '# direct_accepted_state_handoff = ', .false.
         write(unit, '(a)') '# route = frozen-potential diagnostic rebuild from accepted real-space SCF potential; input EF retained'
      end if
      write(unit, '(a,i0,2(es24.16,1x),a,es24.16,a,es24.16)') '# radial_mesh_identity = ', size(ground_states(1)%r), &
         ground_states(1)%a, ground_states(1)%b, 'rmax=', ground_states(1)%rmax, 'sum_r=', sum(ground_states(1)%r)
      write(unit, '(a,a)') '# interaction_route = ', trim(config%interaction_route)
      write(unit, '(a,a)') '# interaction_provenance = ', trim(result%interaction_provenance)
      write(unit, '(a,a)') '# magnetization_kind = ', trim(result%magnetization_kind)
      write(unit, '(a,a)') '# magnetization_source = ', trim(result%magnetization_source)
      write(unit, '(a,a)') '# goldstone_correction = ', trim(result%goldstone_correction_status)
      write(unit, '(a,a)') '# backend = ', trim(config%backend)
      if (trim(config%backend) == tddft_driver_backend_native_rsgf) then
         write(unit, '(a,a)') '# native_rsgf_provider = ', trim(config%native_rsgf_provider)
         write(unit, '(a,a)') '# native_rsgf_provenance = ', trim(result%bare_response_provenance)
         write(unit, '(a)') '# native_rsgf_route = coefficient-GF -> RSGF endpoint augmentation -> LR-04 canonical -> common KXC/Dyson'
      end if
      if (config%write_full_matrix) then
         write(unit, '(a)') '# output_matrix_mode = full'
         write(unit, '(a)') '# columns: q_index omega_Ry matrix_i matrix_j loss_real loss_imag chiKS_real chiKS_imag enhanced_real enhanced_imag'
      else
         write(unit, '(a)') '# output_matrix_mode = trace (complete matrices retained by the driver result)'
         write(unit, '(a)') '# columns: q_index omega_Ry loss_trace_real loss_trace_imag chiKS_trace_real chiKS_trace_imag enhanced_trace_real enhanced_trace_imag'
      end if
      do iq = 1, result%nq
         write(unit, '(a,i0,3(1x,es24.16))') '# q_index = ', iq, result%q_list(:, iq)
         do iw = 1, result%nfrequency
            if (config%write_full_matrix) then
               do j = 1, result%ndim
                  do i = 1, result%ndim
                     write(unit, '(i0,1x,es24.16,1x,2(i0,1x),8(es24.16,1x))') iq, result%frequencies(iw), i, j, &
                        real(result%loss_matrix(i, j, iw, iq), rp), aimag(result%loss_matrix(i, j, iw, iq)), &
                        real(result%ks_susceptibility(i, j, iw, iq), rp), aimag(result%ks_susceptibility(i, j, iw, iq)), &
                        real(result%enhanced_susceptibility(i, j, iw, iq), rp), aimag(result%enhanced_susceptibility(i, j, iw, iq))
                  end do
               end do
            else
               loss_trace = response_operator_trace(response_space, result%loss_matrix(:, :, iw, iq))
               ks_trace = response_operator_trace(response_space, result%ks_susceptibility(:, :, iw, iq))
               enhanced_trace = response_operator_trace(response_space, result%enhanced_susceptibility(:, :, iw, iq))
               write(unit, '(i0,1x,es24.16,1x,6(es24.16,1x))') iq, result%frequencies(iw), &
                  real(loss_trace, rp), aimag(loss_trace), real(ks_trace, rp), aimag(ks_trace), &
                  real(enhanced_trace, rp), aimag(enhanced_trace)
            end if
         end do
      end do
      close(unit)
   end subroutine write_tddft_production_output

   pure real(rp) function radial_volume_integral(state, values) result(value)
      type(radial_ground_state), intent(in) :: state
      real(rp), intent(in) :: values(:)
      integer :: ir

      if (size(values) /= size(state%r)) then
         value = 0.0_rp
         return
      end if
      value = 0.0_rp
      do ir = 1, size(values)
         value = value + radial_simpson_weight(ir, size(values))*state%a*(state%r(ir) + state%b)* &
            4.0_rp*RADIAL_PI*state%r(ir)**2*values(ir)
      end do
   end function radial_volume_integral

   pure real(rp) function radial_volume_norm(state, values) result(value)
      type(radial_ground_state), intent(in) :: state
      real(rp), intent(in) :: values(:)
      integer :: ir

      if (size(values) /= size(state%r)) then
         value = 0.0_rp
         return
      end if
      value = 0.0_rp
      do ir = 1, size(values)
         value = value + radial_simpson_weight(ir, size(values))*state%a*(state%r(ir) + state%b)* &
            4.0_rp*RADIAL_PI*state%r(ir)**2*abs(values(ir))**2
      end do
      value = sqrt(max(value, 0.0_rp))
   end function radial_volume_norm

   real(rp) function relative_site_difference(left, right) result(value)
      complex(rp), intent(in) :: left(:, :), right(:, :)
      real(rp) :: left_norm, right_norm

      left_norm = sqrt(sum(abs(left)**2))
      right_norm = sqrt(sum(abs(right)**2))
      value = sqrt(sum(abs(left - right)**2))/max(left_norm, right_norm, tiny(1.0_rp))
   end function relative_site_difference

end submodule linear_response_run
