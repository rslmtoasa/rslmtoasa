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
   use lr_ward_diagnostics_mod, only: dump_lr_ward_inputs
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
   character(len=*), parameter :: tddft_driver_backend_mills_1u = 'mills_1u'
   character(len=*), parameter :: tddft_driver_backend_projected_juelich = 'projected_juelich'
   character(len=*), parameter :: tddft_driver_route_direct_alsda = 'direct_alsda'

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
      character(len=32) :: backend = tddft_driver_backend_lehmann
      character(len=8) :: projected_selector = 'spd'
      character(len=32) :: native_rsgf_provider = 'auto'
      integer :: gf_integration_points = 2001
      real(rp) :: gf_integration_eta = 0.0_rp
      real(rp) :: gf_energy_margin = 1.0_rp
      logical :: dyson_static_audit = .false.
      logical :: validate_interacting_covariance = .false.
      logical :: write_full_matrix = .true.
      ! Product-space response routes consume the complete accepted
      ! scalar-relativistic radial snapshot, including second derivatives.
      ! This is a capability of the requested representation, not a
      ! property inferred from an individual worker/backend label.
      logical :: requires_second_order_product_basis = .false.
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

      call validate_linear_response_route(config)
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
      runtime%requires_second_order_product_basis = trim(config%formulation) == 'projected' .or. &
         trim(config%representation) == 'product_compact'
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
         case ('radial_points:realspace_gf:alsda')
            runtime%backend = tddft_driver_backend_native_rsgf
            runtime%interaction_route = tddft_driver_route_direct_alsda
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
         case ('mills_1u')
            runtime%backend = tddft_driver_backend_mills_1u
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
      if (trim(config%backend) == tddft_driver_backend_mills_1u) then
         if (trim(config%projected_selector) /= 'd') &
            error stop 'TDDFT input: mills_1u is defined only for the local d shell'
         if (find_gamma_q_index(config%q_list) == 0) &
            error stop 'TDDFT input: mills_1u requires Gamma for the raw static denominator report'
      else if (trim(config%backend) == tddft_driver_backend_projected_juelich) then
         if (trim(config%projected_selector) /= 'd' .or. trim(config%channel) /= 'chi_plus') &
            error stop 'Juelich production requires projection=d and channel=chi_plus'
         if (config%response_lmax >= 0 .and. config%response_lmax /= 4) then
            error stop 'TDDFT input: projected interaction requires response_lmax=4'
         end if
         call select_projected_juelich_eta_indices(config%eta_values, eta_selected_index, eta_holdout_index)
         if (.not. any(sum(abs(config%q_list), dim=1) <= 1.0e-12_rp)) &
            error stop 'TDDFT input: projected_juelich requires Gamma'
      else if (size(config%eta_values) /= 1) then
         error stop 'TDDFT input: n_eta greater than one is restricted to projected validation backends'
      end if
      if (trim(config%backend) == tddft_driver_backend_native_rsgf) then
         if (config%gf_integration_points < 3 .or. mod(config%gf_integration_points,2) == 0) then
            error stop 'TDDFT input: realspace_gf requires an odd gf_integration_points value >= 3'
         end if
         if (config%gf_energy_margin <= 0.0_rp) error stop 'TDDFT input: realspace_gf requires gf_energy_margin positive'
         if (trim(config%native_rsgf_provider) /= 'auto' .and. trim(config%native_rsgf_provider) /= 'block' .and. &
             trim(config%native_rsgf_provider) /= 'chebyshev') then
            error stop 'TDDFT input: realspace_solver provider is invalid'
         end if
      end if
      if (trim(config%backend) == tddft_driver_backend_compact_dyson) then
         if (config%dyson_static_audit .and. find_gamma_q_index(config%q_list) == 0) then
            error stop 'TDDFT input: diagnostics=invariants requires Gamma for compact Dyson'
         end if
         if (config%validate_interacting_covariance) call validate_covariance_q_pair(config%q_list)
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
      call validate_linear_response_route(this%config)
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
      if ((trim(config%backend) == tddft_driver_backend_mills_1u .or. &
           trim(config%backend) == tddft_driver_backend_projected_juelich) .and. state%basis_lmax /= 2) &
         error stop 'projected Juelich-d/Mills-1U require the accepted full spd electronic basis'
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
         if (trim(config%native_rsgf_provider) == 'block') then
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
         state_source = 'accepted_realspace_potential_adapter'
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
      write(unit, '(a,es24.16)') '# electron_count_difference = ', &
         reciprocal_obj%canonical_electron_count - reciprocal_obj%total_electrons
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
      logical :: use_accepted_kspace_scf, requires_second_order_product_basis

      use_accepted_kspace_scf = .false.
      if (present(accepted_kspace_scf)) use_accepted_kspace_scf = accepted_kspace_scf
      requires_second_order_product_basis = config%requires_second_order_product_basis

      if (.not. config%enabled) return
      call validate_tddft_production_capability(config, control_obj, lattice_obj, hamiltonian_obj, reciprocal_obj)
      if (.not. scf_converged) error stop 'TDDFT production driver: an accepted converged ground state is required'
      if (trim(config%backend) == tddft_driver_backend_mills_1u .and. .not. use_accepted_kspace_scf) then
         error stop 'MILLS-1U requires the accepted reciprocal state and will not rebuild or re-SCF it'
      end if
      if (lattice_obj%nrec < 1) error stop 'TDDFT production driver: no accepted response sites are available'

      if (trim(config%backend) == tddft_driver_backend_mills_1u) then
         call validate_accepted_kspace_scf_handoff(reciprocal_obj, energy_obj)
         call reciprocal_obj%require_replicated_k_workset('MILLS-1U production driver')
         call lr_snapshot_from_reciprocal(reciprocal_obj, left_state)
         call write_tddft_state_artifact(trim(config%output_file)//'.state', reciprocal_obj, left_state, .true.)
         allocate(endpoints(size(config%q_list, 2)))
         do iq = 1, size(config%q_list, 2)
            call lr_q_endpoint_from_reciprocal(reciprocal_obj, left_state, config%q_list(:, iq), endpoints(iq))
         end do
         call run_tddft_mills_1u(config, left_state, endpoints, reciprocal_obj, lattice_obj, scf_converged)
         return
      end if

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
         call radial_basis_from_snapshot(ground_states(isite), radial_bases(isite), &
            requires_second_order_product_basis, trim(config%backend))
         call radial_bases(isite)%require_supported(trim(reciprocal_obj%reciprocal_mode), &
            effective_hamiltonian_order(reciprocal_obj, hamiltonian_obj), .true., .true., &
            control_obj%has_soc(), .false.)
      end do
      call response_space%initialize(nsite, response_lmax, ground_states(1)%r, ground_states(1)%a, ground_states(1)%b, 1)

      if (use_accepted_kspace_scf) then
         call validate_accepted_kspace_scf_handoff(reciprocal_obj, energy_obj)
      else
         ! Accepted-real-space-potential adapter used by the existing spherical
         ! direct-ALSDA and bare-response fixtures. Preserve the accepted EF.
         ! This is shared state preparation, not another interaction method;
         ! compact ALSDA and Juelich-d require reciprocal SCF below.
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
            reciprocal_obj, lattice_obj, control_obj, hamiltonian_obj)
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

   !> Production Mills-1U response built from the accepted coefficient state.
   subroutine run_tddft_mills_1u(config, left_state, endpoints, reciprocal_obj, lattice_obj, scf_converged)
      type(tddft_runtime_config), intent(in) :: config
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(reciprocal), intent(inout) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      logical, intent(in) :: scf_converged

      type(mills_1u_bare_request) :: bare_request, opposite_request, raw_request
      type(mills_1u_bare_result) :: bare_result, opposite_bare_result, raw_bare_result
      type(mills_1u_model_reduction_result) :: reduction
      type(projected_dyson_result) :: dyson_result, opposite_dyson_result, raw_dyson_result
      real(rp), allocatable :: center_band(:, :, :), scalarity(:), delta_d(:), moment_d(:), physical_u(:)
      complex(rp), allocatable :: cx(:, :, :), spin_difference(:, :), first_spin_difference(:, :), &
         enim_spin_difference(:, :), hoh_spin_difference(:, :), first_h(:, :), &
         spin_difference_k(:, :, :), spin_difference_mean(:, :)
      complex(rp) :: value
      real(rp) :: accepted_total_moment, weight, residual2, mills2, dnorm2, local2, dsector2, nond2, offdiag2
      real(rp) :: hoh2, second_order2, kdependent2, covariance_error, action_norm, mode_norm, scalarity_limit
      integer :: nsite, norb, nbasis, nk, first_atom, last_atom, site, orbital, m, spin
      integer :: k, i, j, site_i, site_j, orbital_i, orbital_j, up_i, down_i, up_j, down_j
      integer :: atom, itype, iq, iq_negative, ifrequency, ieta, unit, gamma_index
      character(len=8) :: selector

      if (.not. scf_converged) error stop 'MILLS-1U: a converged accepted state is required'
      if (.not. left_state%initialized .or. left_state%has_soc .or. .not. left_state%collinear .or. &
          .not. left_state%orthogonal .or. trim(left_state%reciprocal_mode) /= 'ham_only' .or. &
          trim(left_state%hamiltonian_order) /= 'second') then
         error stop 'MILLS-1U: require accepted scalar-relativistic collinear second-order ham_only state'
      end if
      nsite = lattice_obj%nrec
      norb = 9
      nbasis = 2*nsite*norb
      nk = left_state%nk
      if (.not. allocated(reciprocal_obj%hk_bulk) .or. .not. allocated(reciprocal_obj%k_points)) then
         error stop 'MILLS-1U: the accepted reciprocal Hamiltonian and k mesh are required'
      end if
      if (nsite < 1 .or. left_state%nbasis /= nbasis .or. size(reciprocal_obj%hk_bulk, 1) /= nbasis .or. &
          size(reciprocal_obj%hk_bulk, 2) /= nbasis .or. size(reciprocal_obj%hk_bulk, 3) /= nk .or. &
          size(endpoints) /= size(config%q_list, 2) .or. size(reciprocal_obj%k_points, 2) /= nk) then
         error stop 'MILLS-1U: accepted reciprocal state must contain the complete spd k-space eigensystem'
      end if
      if (.not. associated(reciprocal_obj%hamiltonian)) then
         error stop 'MILLS-1U: accepted Hamiltonian owner is unavailable'
      end if
      if (.not. allocated(reciprocal_obj%hamiltonian%enim) .or. .not. allocated(lattice_obj%ib)) then
         error stop 'MILLS-1U: accepted RS-LMTO second-order Hamiltonian is incomplete'
      end if
      if (size(lattice_obj%ib) < nsite) then
         error stop 'MILLS-1U: accepted site-to-potential type map is incomplete'
      end if
      if (size(reciprocal_obj%hamiltonian%enim, 1) < 2*norb .or. &
          size(reciprocal_obj%hamiltonian%enim, 2) < 2*norb .or. &
          size(reciprocal_obj%hamiltonian%enim, 3) < 1) then
         error stop 'MILLS-1U: accepted onsite second-order Hamiltonian dimensions are inconsistent'
      end if
      first_atom = lattice_obj%nbulk + 1
      last_atom = first_atom + nsite - 1
      if (.not. allocated(lattice_obj%symbolic_atoms)) then
         error stop 'MILLS-1U: accepted transformed symbolic-atom potentials are unavailable'
      end if
      if (last_atom > size(lattice_obj%symbolic_atoms)) then
         error stop 'MILLS-1U: accepted transformed symbolic-atom potentials are unavailable'
      end if
      allocate(center_band(3, 2, nsite), cx(norb, 2, nsite), scalarity(nsite), delta_d(nsite), &
         moment_d(nsite), physical_u(nsite))
      do site = 1, nsite
         atom = first_atom + site - 1
         if (lattice_obj%symbolic_atoms(atom)%potential%lmax /= 2 .or. &
             .not. allocated(lattice_obj%symbolic_atoms(atom)%potential%center_band) .or. &
             .not. allocated(lattice_obj%symbolic_atoms(atom)%potential%cx)) then
            error stop 'MILLS-1U: accepted native potential must contain the full spherical spd basis'
         end if
         if (size(lattice_obj%symbolic_atoms(atom)%potential%center_band, 1) < 3 .or. &
             size(lattice_obj%symbolic_atoms(atom)%potential%center_band, 2) /= 2 .or. &
             size(lattice_obj%symbolic_atoms(atom)%potential%cx, 1) < 9 .or. &
             size(lattice_obj%symbolic_atoms(atom)%potential%cx, 2) /= 2) then
            error stop 'MILLS-1U: accepted native potential has inconsistent full-spd dimensions'
         end if
         center_band(:, :, site) = lattice_obj%symbolic_atoms(atom)%potential%center_band(1:3, 1:2)
         cx(:, :, site) = lattice_obj%symbolic_atoms(atom)%potential%cx(1:9, 1:2)
      end do
      call mills_1u_splitting_from_centers(center_band, cx, delta_d, scalarity)
      scalarity_limit = 100.0_rp*epsilon(1.0_rp)*max(1.0_rp, maxval(abs(center_band(3, :, :))))
      if (maxval(scalarity) > scalarity_limit) then
         error stop 'MILLS-1U: spherical ASA d center is not identical across all five m channels'
      end if
      call mills_1u_coefficient_moment(left_state, nsite, 2, moment_d)
      physical_u = delta_d/moment_d
      if (any(.not. ieee_is_finite(physical_u))) error stop 'MILLS-1U: Delta_d/m_d is non-finite'

      accepted_total_moment = 0.0_rp
      do site = 1, nsite
         atom = first_atom + site - 1
         accepted_total_moment = accepted_total_moment + lattice_obj%symbolic_atoms(atom)%potential%mtot
      end do

      ! Decompose the accepted full spd H_down-H_up at every k. The first-order
      ! Fourier transform reuses the accepted live real-space blocks; subtracting
      ! it from accepted H^(2) and removing the explicit e_nu term isolates the
      ! spin difference of -HOH without duplicating the RS-LMTO assembler.
      allocate(spin_difference(nsite*norb, nsite*norb), first_spin_difference(nsite*norb, nsite*norb), &
         enim_spin_difference(nsite*norb, nsite*norb), hoh_spin_difference(nsite*norb, nsite*norb), &
         first_h(nbasis, nbasis), spin_difference_k(nsite*norb, nsite*norb, nk), &
         spin_difference_mean(nsite*norb, nsite*norb))
      enim_spin_difference = cmplx(0.0_rp, 0.0_rp, rp)
      do site_i = 1, nsite
         itype = lattice_obj%ib(site_i)
         if (itype < 1 .or. itype > size(reciprocal_obj%hamiltonian%enim, 3)) then
            error stop 'MILLS-1U: accepted site-to-potential type map is invalid'
         end if
         do site_j = 1, nsite
            if (site_i /= site_j) cycle
            do orbital_i = 1, norb
               i = (site_i - 1)*norb + orbital_i
               do orbital_j = 1, norb
                  j = (site_j - 1)*norb + orbital_j
                  enim_spin_difference(i, j) = reciprocal_obj%hamiltonian%enim(norb + orbital_i, norb + orbital_j, itype) - &
                     reciprocal_obj%hamiltonian%enim(orbital_i, orbital_j, itype)
               end do
            end do
         end do
      end do
      residual2 = 0.0_rp
      mills2 = 0.0_rp
      dnorm2 = 0.0_rp
      local2 = 0.0_rp
      dsector2 = 0.0_rp
      nond2 = 0.0_rp
      offdiag2 = 0.0_rp
      hoh2 = 0.0_rp
      second_order2 = 0.0_rp
      spin_difference_mean = cmplx(0.0_rp, 0.0_rp, rp)
      do k = 1, nk
         spin_difference = cmplx(0.0_rp, 0.0_rp, rp)
         call reciprocal_obj%fourier_transform_hamiltonian(reciprocal_obj%k_points(:, k), first_h)
         first_spin_difference = cmplx(0.0_rp, 0.0_rp, rp)
         do site_i = 1, nsite
            do orbital_i = 1, norb
               i = (site_i - 1)*norb + orbital_i
               up_i = (site_i - 1)*2*norb + orbital_i
               down_i = up_i + norb
               do site_j = 1, nsite
                  do orbital_j = 1, norb
                     j = (site_j - 1)*norb + orbital_j
                     up_j = (site_j - 1)*2*norb + orbital_j
                     down_j = up_j + norb
                     spin_difference(i, j) = reciprocal_obj%hk_bulk(down_i, down_j, k) - &
                        reciprocal_obj%hk_bulk(up_i, up_j, k)
                     first_spin_difference(i, j) = first_h(down_i, down_j) - first_h(up_i, up_j)
                  end do
               end do
            end do
         end do
         spin_difference_k(:, :, k) = spin_difference
         spin_difference_mean = spin_difference_mean + left_state%k_weights(k)*spin_difference
         call evaluate_mills_1u_model_reduction(spin_difference, delta_d, nsite, norb, reduction)
         weight = left_state%k_weights(k)
         residual2 = residual2 + weight*reduction%residual_norm**2
         mills2 = mills2 + weight*reduction%mills_local_d_norm**2
         dnorm2 = dnorm2 + weight*reduction%spin_difference_norm**2
         local2 = local2 + weight*reduction%local_d_shell_residual_norm**2
         dsector2 = dsector2 + weight*reduction%remaining_d_sector_norm**2
         nond2 = nond2 + weight*reduction%non_d_sector_norm**2
         offdiag2 = offdiag2 + weight*reduction%site_offdiagonal_norm**2
         hoh_spin_difference = spin_difference - first_spin_difference - enim_spin_difference
         hoh2 = hoh2 + weight*sum(abs(hoh_spin_difference)**2)
         second_order2 = second_order2 + weight*sum(abs(spin_difference - first_spin_difference)**2)
      end do
      kdependent2 = 0.0_rp
      do k = 1, nk
         kdependent2 = kdependent2 + left_state%k_weights(k)* &
            sum(abs(spin_difference_k(:, :, k) - spin_difference_mean)**2)
      end do

      open(newunit=unit, file=trim(config%output_file), status='replace', action='write')
      write(unit, '(a)') '# MILLS-1U — CONTROLLED RS-LMTO MODEL PROJECTION'
      write(unit, '(a,a)') '# build_version = ', trim(tddft_build_version)
      write(unit, '(a)') '# state_provenance = accepted converged reciprocal SCF eigenpairs + same-run live transformed potential'
      write(unit, '(a,l1)') '# scf_converged = ', scf_converged
      write(unit, '(a)') '# hidden_rescf = false'
      write(unit, '(a,3(i0,1x))') '# k_mesh = ', reciprocal_obj%nk_mesh
      write(unit, '(a,i0)') '# nk = ', nk
      write(unit, '(a,es24.16)') '# k_weight_sum = ', sum(left_state%k_weights)
      write(unit, '(a,es24.16)') '# temperature_K = ', left_state%temperature
      write(unit, '(a,es24.16)') '# fermi_level_Ry = ', left_state%fermi_level
      write(unit, '(a,es24.16)') '# accepted_total_moment_muB = ', accepted_total_moment
      write(unit, '(a)') '# interaction = one scalar physical U_i=Delta_i^d/m_i^d acting only on local d orbitals'
      write(unit, '(a)') '# moment = sum_k w_k f_n sum_mu(d)(|c_up|^2-|c_down|^2); accepted orthogonal coefficients, muB'
      write(unit, '(a)') '# Dyson mapping = chi_code=2 chi_Kubo; solver K=-U/2 gives I-chi_code*K=I+U*chi_Kubo'
      write(unit, '(a)') '# no radial-overlap factors; full spd bands; d-only local S+ vertex'
      write(unit, '(a,es24.16)') '# d_center_spherical_identity_max_abs = ', maxval(scalarity)
      write(unit, '(a)') '# MILLS_SITE fields: site C_d_up_Ry C_d_down_Ry Delta_d_Ry m_d_muB physical_U_Ry d_center_identity_abs'
      do site = 1, nsite
         write(unit, '(a,1x,i0,1x,6(es24.16,1x))') 'MILLS_SITE', site, center_band(3, 1, site), &
            center_band(3, 2, site), delta_d(site), moment_d(site), physical_u(site), scalarity(site)
      end do
      write(unit, '(a,7(es24.16,1x))') 'RS_LMTO_MILLS_RESIDUAL', sqrt(dnorm2), sqrt(mills2), sqrt(residual2), &
         sqrt(residual2)/max(sqrt(dnorm2), tiny(1.0_rp)), sqrt(local2), sqrt(dsector2), sqrt(nond2)
      write(unit, '(a,3(es24.16,1x))') 'RS_LMTO_NONLOCAL_AND_HOH', sqrt(kdependent2), sqrt(offdiag2), sqrt(hoh2)
      write(unit, '(a,es24.16)') '# model_reduction_residual_over_retained_Mills_term = ', &
         sqrt(residual2)/max(sqrt(mills2), tiny(1.0_rp))
      write(unit, '(a)') '# residual fields: full_D_RMS Mills_d_RMS R_RMS R_over_D local_d_remainder_RMS d_sector_RMS non_d_RMS'
      write(unit, '(a)') '# nonlocal fields: k_dependent_D_RMS unit_cell_site_offdiagonal_D_RMS'
      write(unit, '(a)') '# second-order fields: HOH_spin_difference_RMS'
      write(unit, '(a,es24.16)') '# second_order_total_spin_difference_RMS = ', sqrt(second_order2)
      write(unit, '(a)') '# residual is D(k)-sum_i Delta_i^d P_i^d; no fit, threshold, or interaction retuning'
      write(unit, '(a)') '# columns: route eta_Ry q_index omega_Ry row col bare_Re bare_Im chi_Re chi_Im loss_Re loss_Im min_sv max_sv cond min_abs_eig dyson_residual loss_trace minus_im_trace_over_pi'
      write(unit, '(a)') '# RAW_GAMMA fields: eta_Ry min_sv min_abs_eig moment_action_norm moment_action_relative'
      do ieta = 1, size(config%eta_values)
         do iq = 1, size(config%q_list, 2)
            bare_request%q = config%q_list(:, iq)
            bare_request%frequencies = config%frequencies
            bare_request%eta = config%eta_values(ieta)
            bare_request%channel = config%channel
            bare_request%nsite = nsite
            bare_request%orbital_lmax = 2
            bare_request%electronic_state => left_state
            bare_request%q_endpoint_state => endpoints(iq)
            call evaluate_mills_1u_bare(bare_request, bare_result)
            call evaluate_mills_1u_dyson(bare_result, physical_u, dyson_result)
            do ifrequency = 1, size(config%frequencies)
               do j = 1, nsite
                  do i = 1, nsite
                     write(unit, '(a,1x,a,1x,es24.16,1x,i0,1x,es24.16,1x,2(i0,1x),13(es24.16,1x))') &
                        'ROW', 'd', config%eta_values(ieta), iq, config%frequencies(ifrequency), i, j, &
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

            gamma_index = find_gamma_q_index(config%q_list)
            if (gamma_index == iq) then
               raw_request = bare_request
               raw_request%q = config%q_list(:, iq)
               raw_request%frequencies = [0.0_rp]
               raw_request%channel = 'chi_plus'
               call evaluate_mills_1u_bare(raw_request, raw_bare_result)
               call evaluate_mills_1u_dyson(raw_bare_result, physical_u, raw_dyson_result)
               action_norm = sqrt(sum(abs(matmul(raw_dyson_result%denominator(:, :, 1), &
                  cmplx(moment_d, 0.0_rp, rp)))**2))
               mode_norm = sqrt(sum(moment_d**2))
               write(unit, '(a,1x,es24.16,1x,4(es24.16,1x))') 'RAW_GAMMA', config%eta_values(ieta), &
                  raw_dyson_result%denominator_min_singular_value(1), &
                  raw_dyson_result%minimum_magnitude_eigenvalue(1), action_norm, action_norm/max(mode_norm,tiny(1.0_rp))
            end if
         end do
      end do

      if (trim(config%channel) == 'chi_plus') then
         do iq = 1, size(config%q_list, 2)
            if (sum(abs(config%q_list(:, iq))) <= 2.0e-12_rp) cycle
            iq_negative = find_matching_q(config%q_list, -config%q_list(:, iq), 2.0e-12_rp)
            if (iq_negative == 0 .or. iq_negative < iq) cycle
            bare_request%q = config%q_list(:, iq)
            bare_request%frequencies = config%frequencies
            bare_request%eta = minval(config%eta_values)
            bare_request%channel = 'chi_plus'
            bare_request%nsite = nsite
            bare_request%orbital_lmax = 2
            bare_request%electronic_state => left_state
            bare_request%q_endpoint_state => endpoints(iq)
            opposite_request = bare_request
            opposite_request%q = config%q_list(:, iq_negative)
            opposite_request%frequencies = -config%frequencies
            opposite_request%channel = 'chi_minus'
            opposite_request%q_endpoint_state => endpoints(iq_negative)
            call evaluate_mills_1u_bare(bare_request, bare_result)
            call evaluate_mills_1u_dyson(bare_result, physical_u, dyson_result)
            call evaluate_mills_1u_bare(opposite_request, opposite_bare_result)
            call evaluate_mills_1u_dyson(opposite_bare_result, physical_u, opposite_dyson_result)
            covariance_error = maxval(abs(dyson_result%enhanced_chi - conjg(opposite_dyson_result%enhanced_chi)))
            write(unit, '(a,1x,2(i0,1x),es24.16)') 'Q_MINUS_Q_COVARIANCE', iq, iq_negative, covariance_error
         end do
      end if
      close(unit)
   end subroutine run_tddft_mills_1u

   !> DRESP-05 material seam.  The static Gamma response constructs a frozen
   !> Literature-traceable Juelich-d route. The projector is frozen at EF,
   !> the interaction comes from the eta->0+ static response, and that one U
   !> is then used for the requested site-space Dyson response.
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

      type(juelich_d_projector) :: projector
      type(projected_juelich_request) :: interaction_request
      type(projected_juelich_result), allocatable :: eta_results(:)
      type(projected_juelich_result) :: static_result
      type(projected_dyson_request) :: dyson_request
      type(projected_dyson_result) :: dyson_result
      real(rp), allocatable :: moment(:)
      complex(rp), allocatable :: chi_eta(:, :, :), chi_static(:, :), chi_dynamic(:, :, :), denominator(:, :), action(:)
      real(rp) :: relative_estimate, fit_residual, imaginary_ratio, u_eta_relative, u_eta_imaginary
      real(rp) :: moment_norm, eta, value
      integer :: gamma_index, nsite, neta, ieta, iq, unit, i, j, iw, isite

      nsite = size(radial_bases)
      neta = size(config%eta_values)
      if (trim(config%projected_selector) /= 'd' .or. trim(config%channel) /= lr_channel_plus) then
         error stop 'Juelich-d requires projection=d and channel=chi_plus'
      end if
      if (nsite < 1 .or. lattice_obj%nrec /= nsite .or. response_space%nsite /= nsite .or. &
          size(endpoints) /= size(config%q_list, 2) .or. neta < 4 .or. size(ground_states) /= nsite) then
         error stop 'Juelich-d material seam: inconsistent accepted-state/site dimensions or eta ladder'
      end if
      if (any(.not. ieee_is_finite(config%eta_values)) .or. any(config%eta_values <= 0.0_rp)) then
         error stop 'Juelich-d static interaction requires finite positive eta values'
      end if
      do ieta = 1, neta - 1
         if (config%eta_values(ieta) <= config%eta_values(ieta + 1)) then
            error stop 'Juelich-d static interaction requires a strictly decreasing eta ladder'
         end if
      end do
      gamma_index = find_gamma_q_index(config%q_list)
      if (gamma_index == 0) error stop 'Juelich-d static interaction requires Gamma'
      if (size(config%frequencies) < 1) error stop 'Juelich-d requires a nonempty frequency grid'

      call projector%initialize(radial_bases, left_state%fermi_level)
      allocate(moment(nsite), chi_eta(nsite, nsite, neta), chi_static(nsite, nsite), eta_results(neta))
      call projector%moment_from_state(radial_bases, left_state, moment)
      do ieta = 1, neta
         eta = config%eta_values(ieta)
         call evaluate_juelich_d_lehmann_chi0(projector, radial_bases, left_state, endpoints(gamma_index), &
            config%q_list(:, gamma_index), [0.0_rp], eta, chi_eta(:, :, ieta:ieta))
         interaction_request%selector = 'd'
         interaction_request%projected_moment = moment
         interaction_request%static_chi0 = chi_eta(:, :, ieta)
         interaction_request%static_eta = eta
         interaction_request%q = config%q_list(:, gamma_index)
         interaction_request%channel = lr_channel_plus
         interaction_request%state_provenance = 'accepted reciprocal state; frozen normalized R_d(EF) projection'
         interaction_request%chi0_provenance = 'Lounis Eq. 2/22 Lehmann bubble at positive eta; diagnostic U(eta) only'
         call evaluate_projected_juelich_interaction(interaction_request, eta_results(ieta))
      end do

      call extrapolate_juelich_static_chi0(config%eta_values, chi_eta, chi_static, relative_estimate, &
         fit_residual, imaginary_ratio)
      interaction_request%static_chi0 = chi_static
      interaction_request%static_eta = 0.0_rp
      interaction_request%state_provenance = 'accepted reciprocal state; eta->0+ extrapolated frozen d projection'
      interaction_request%chi0_provenance = 'Re chi0 fitted linearly in eta^2 from final four retarded eta samples'
      call evaluate_projected_juelich_interaction(interaction_request, static_result)
      if (static_result%rank /= nsite .or. static_result%real_rank /= nsite .or. &
          static_result%condition_number > projected_juelich_condition_limit .or. &
          static_result%real_condition_number > projected_juelich_condition_limit .or. &
          static_result%relative_real_constrained_residual > projected_juelich_residual_tolerance) then
         error stop 'BLOCKED - JUELICH STATIC LIMIT NOT CLOSED: static Gamma U=M solve is unsupported'
      end if
      u_eta_relative = maxval(abs(real(eta_results(neta)%interaction_U_complex, rp) - static_result%interaction_U_real))/ &
         max(maxval(abs(static_result%interaction_U_real)), tiny(1.0_rp))
      u_eta_imaginary = maxval(abs(aimag(eta_results(neta)%interaction_U_complex)))/ &
         max(maxval(abs(real(eta_results(neta)%interaction_U_complex, rp))), tiny(1.0_rp))
      if (relative_estimate > 5.0e-3_rp .or. fit_residual > 1.0e-2_rp .or. imaginary_ratio > 5.0e-3_rp .or. &
          u_eta_relative > 5.0e-3_rp .or. u_eta_imaginary > 5.0e-3_rp) then
         error stop 'BLOCKED - JUELICH STATIC LIMIT NOT CLOSED: eta extrapolation or U(eta) has not converged'
      end if

      allocate(denominator(nsite, nsite), action(nsite))
      denominator = -chi_static*spread(static_result%interaction_U_real, 1, nsite)
      do i = 1, nsite
         denominator(i, i) = denominator(i, i) + cmplx(1.0_rp, 0.0_rp, rp)
      end do
      action = matmul(denominator, cmplx(moment, 0.0_rp, rp))
      moment_norm = sqrt(sum(moment**2))
      value = sqrt(sum(abs(action)**2))/max(moment_norm, tiny(1.0_rp))

      open(newunit=unit, file=trim(config%output_file), status='replace', action='write')
      write(unit, '(a)') '# method = Juelich-d (Lounis et al. PRB 83, 035109 (2011), Eqs. 2, 22, 40-47)'
      write(unit, '(a)') '# state_source = accepted_kspace_scf_cache; accepted state is reused without re-SCF'
      write(unit, '(a)') '# accepted_state_cache_reused = T'
      write(unit, '(a,i0)') '# nk = ', left_state%nk
      write(unit, '(a,3(i0,1x))') '# k_mesh = ', reciprocal_obj%nk_mesh
      write(unit, '(a,es24.16)') '# fermi_level_Ry = ', left_state%fermi_level
      write(unit, '(a,es24.16)') '# temperature_K = ', left_state%temperature
      write(unit, '(a)') '# spin convention = M_d = mu_B*(n_up_d - n_down_d); code moment unit is mu_B'
      write(unit, '(a)') '# frozen radial projector = phi_d + (EF-enu_work)*phidot_d, separately for up/down'
      write(unit, '(a)') '# projector measure = log-mesh Simpson dr with scalar-relativistic large/small metric at EF'
      write(unit, '(a)') '# angular convention = unit-normalized Y_lm; Eq. 47 site chi0: no extra 4*pi; spin-flip factor = 1'
      write(unit, '(a,*(es24.16,1x))') '# projector_norm_up = ', projector%normalization(:, 1)
      write(unit, '(a,*(es24.16,1x))') '# projector_norm_down = ', projector%normalization(:, 2)
      write(unit, '(a,*(es24.16,1x))') '# projected_d_moment_muB = ', moment
      write(unit, '(a,es24.16)') '# static_eta_fit_relative_estimate = ', relative_estimate
      write(unit, '(a,es24.16)') '# static_eta_fit_relative_residual = ', fit_residual
      write(unit, '(a,es24.16)') '# smallest_eta_imaginary_chi_ratio = ', imaginary_ratio
      write(unit, '(a,es24.16)') '# smallest_eta_U_relative_to_static = ', u_eta_relative
      write(unit, '(a,es24.16)') '# smallest_eta_U_imaginary_ratio = ', u_eta_imaginary
      write(unit, '(a,i0)') '# static_sumrule_rank = ', static_result%rank
      write(unit, '(a,es24.16)') '# static_sumrule_condition = ', static_result%condition_number
      write(unit, '(a,*(es24.16,1x))') '# U_Juelich_d_static_energy_per_muB = ', static_result%interaction_U_real
      write(unit, '(a,es24.16)') '# static_sumrule_relative_residual = ', static_result%relative_real_constrained_residual
      write(unit, '(a,es24.16)') '# Dyson_denominator_magnetic_mode_relative_action = ', value
      write(unit, '(a)') '# CHI0_STATIC columns: row col Re Im'
      do i = 1, nsite
         do j = 1, nsite
            write(unit, '(a,1x,i0,1x,i0,1x,2(es24.16,1x))') 'CHI0_STATIC', i, j, &
               real(chi_static(i, j), rp), aimag(chi_static(i, j))
         end do
      end do
      write(unit, '(a)') '# JUELICH_ETA columns: eta site U_Re U_Im rank condition complex_residual real_constraint_residual'
      do ieta = 1, neta
         do isite = 1, nsite
            write(unit, '(a,1x,es24.16,1x,i0,1x,2(es24.16,1x),i0,1x,3(es24.16,1x))') 'JUELICH_ETA', &
               config%eta_values(ieta), isite, real(eta_results(ieta)%interaction_U_complex(isite), rp), &
               aimag(eta_results(ieta)%interaction_U_complex(isite)), eta_results(ieta)%rank, &
               eta_results(ieta)%condition_number, eta_results(ieta)%relative_residual, &
               eta_results(ieta)%relative_real_constrained_residual
         end do
      end do
      write(unit, '(a)') '# DYNAMICS columns: eta q_index omega row col chi0_Re chi0_Im chi_Re chi_Im min_sv condition dyson_residual'
      allocate(chi_dynamic(nsite, nsite, size(config%frequencies)), dyson_request%interaction_U(nsite))
      dyson_request%selector = 'd'
      dyson_request%channel = lr_channel_plus
      dyson_request%interaction_U = static_result%interaction_U_real
      dyson_request%interaction_provenance = 'Eq. 47 interaction from extrapolated Juelich-d static chi0'
      do ieta = 1, neta
         do iq = 1, size(config%q_list, 2)
            call evaluate_juelich_d_lehmann_chi0(projector, radial_bases, left_state, endpoints(iq), &
               config%q_list(:, iq), config%frequencies, config%eta_values(ieta), chi_dynamic)
            dyson_request%q = config%q_list(:, iq)
            dyson_request%frequencies = config%frequencies
            dyson_request%eta = config%eta_values(ieta)
            dyson_request%bare_chi = chi_dynamic
            dyson_request%bare_provenance = 'Eq. 2 spin-flip Lehmann bubble with frozen normalized d projector'
            call evaluate_projected_dyson(dyson_request, dyson_result)
            do iw = 1, size(config%frequencies)
               do i = 1, nsite
                  do j = 1, nsite
                     write(unit, '(a,1x,5(i0,1x),7(es24.16,1x),i0)') 'JUELICH_DYNAMICS', ieta, iq, iw, i, j, &
                        real(chi_dynamic(i, j, iw), rp), aimag(chi_dynamic(i, j, iw)), &
                        real(dyson_result%enhanced_chi(i, j, iw), rp), aimag(dyson_result%enhanced_chi(i, j, iw)), &
                        dyson_result%denominator_min_singular_value(iw), dyson_result%condition_number(iw), &
                        dyson_result%dyson_residual_relative(iw), dyson_result%solve_info(iw)
                  end do
               end do
            end do
         end do
      end do
      close(unit)
      deallocate(moment, chi_eta, chi_static, chi_dynamic, denominator, action, eta_results, dyson_request%interaction_U)
   end subroutine run_tddft_projected_juelich


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
         if (projection_index == 1 .and. trim(config%projected_selector) == 'spd') cycle
         if (projection_index == 2 .and. trim(config%projected_selector) == 'd') cycle
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
                                      reciprocal_obj, lattice_obj, control_obj, hamiltonian_obj)
      type(tddft_runtime_config), intent(in) :: config
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(radial_ground_state), target, intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      type(control), intent(in) :: control_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj

      type(lmto_product_response_basis), target :: product_plus, product_minus
      type(lmto_product_response_basis), pointer :: product, opposite_product
      type(lr_product_ks_susceptibility_request) :: bare_request, opposite_bare_request
      type(lr_product_ks_susceptibility_result) :: bare_result, opposite_bare_result
      type(lr_alsda_kernel_request) :: kxc_request
      type(lr_alsda_kernel_result) :: kxc_result
      type(tddft_dyson_request) :: dyson_request, opposite_dyson_request
      type(tddft_dyson_result) :: dyson_result, opposite_dyson_result
      real(rp), allocatable :: pauli_response_magnetization(:, :), static_frequency(:)
      complex(rp), allocatable :: interaction(:, :), opposite_interaction(:, :)
      real(rp) :: static_eta(5), static_min_sv(5), static_max_sv(5), static_condition(5), static_min_eigen(5)
      real(rp) :: static_residual(5), static_residual_relative(5), static_residual_infinity(5)
      character(len=256) :: static_status(5)
      real(rp), allocatable :: covariance_residual(:), covariance_loss_difference(:), covariance_min_sv_difference(:), &
         covariance_condition_difference(:)
      logical :: rank_stable, finite_response, covariance_checked
      real(rp) :: state_mesh_max, state_weight_max, state_ef_diff, state_eigen_max, state_occ_max
      real(rp) :: state_projector_max, state_projector_frobenius, k_fingerprint(5), weight_sum, magnetic_moment
      real(rp), allocatable :: ground_state_magnetization(:, :), accepted_xc_field(:, :)
      complex(rp), allocatable :: ward_field(:), ward_target(:), ward_response(:), ward_residual(:)
      complex(rp) :: ward_overlap
      real(rp) :: ward_abs, ward_rel, sr_norm, difference_norm, integrated_sr, integrated_pauli
      real(rp) :: magnetic_weight, magnetic_difference, magnetic_norm
      integer :: ir, magnetic_last
      real(rp) :: trace_loss, trace_loss_imag
      integer :: ndim, nfrequency, nq, iq, iw, i, j, unit, gamma_index, positive_q_index, negative_q_index
      character(len=1024) :: ward_dump_path
      integer :: ward_dump_status
      character(len=256) :: state_file
      logical :: static_audit_enabled, covariance_enabled

      if (trim(config%interaction_route) /= tddft_driver_route_direct_alsda) then
         error stop 'compact_dyson requires direct ALSDA as the production interaction route'
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
      ! Explicit read-only investigation hook; no diagnostic feeds production.
      call get_environment_variable('RSLMTO_LR_WARD_DUMP',ward_dump_path,status=ward_dump_status)
      if (ward_dump_status == 0 .and. len_trim(ward_dump_path) > 0 .and. rank == 0) &
         call dump_lr_ward_inputs(trim(ward_dump_path),response_space,radial_bases,ground_states,left_state, &
            product,reciprocal_obj,lattice_obj,hamiltonian_obj)
      ndim = product%product_dimension
      nfrequency = size(config%frequencies)
      nq = size(config%q_list, 2)

      allocate(pauli_response_magnetization(size(ground_states), response_space%npoint))
      call compute_accepted_pauli_magnetization(reciprocal_obj, lattice_obj%symbolic_atoms, lattice_obj%nbulk, pauli_response_magnetization)
      call prepare_direct_alsda_request(response_space, ground_states, pauli_response_magnetization, kxc_request)
      call evaluate_lr_alsda_kernel(kxc_request, kxc_result)
      allocate(interaction(ndim, ndim))
      call compact_project_local_operator(response_space, product, cmplx(kxc_result%pointwise_kernel, 0.0_rp, rp), interaction)
      static_eta = [0.01_rp, 0.005_rp, 0.0025_rp, 0.00125_rp, 0.000625_rp]
      allocate(ground_state_magnetization(size(ground_states), response_space%npoint), &
         accepted_xc_field(size(ground_states), response_space%npoint))
      do i = 1, size(ground_states)
         ground_state_magnetization(i, :) = ground_states(i)%n_up - ground_states(i)%n_down
         accepted_xc_field(i, :) = 0.5_rp*(ground_states(i)%vxc_up - ground_states(i)%vxc_down)
      end do
      allocate(ward_field(ndim), ward_target(ndim), ward_response(ndim), ward_residual(ndim))
      ! Same spherical-harmonic/metric mapping for an independent SCF field
      ! and the response magnetization; Kxc is not used to construct the field.
      call compact_project_magnetization(response_space, product, accepted_xc_field, ward_field)
      call compact_project_magnetization(response_space, product, pauli_response_magnetization, ward_target)
      integrated_sr = 0.0_rp
      integrated_pauli = 0.0_rp
      sr_norm = 0.0_rp
      difference_norm = 0.0_rp
      magnetic_difference = 0.0_rp
      magnetic_norm = 0.0_rp
      do i = 1, size(ground_states)
         integrated_sr = integrated_sr + sum(response_space%radial_weights*ground_state_magnetization(i, :))
         integrated_pauli = integrated_pauli + sum(response_space%radial_weights*pauli_response_magnetization(i, :))
         sr_norm = sr_norm + sum(response_space%radial_weights*ground_state_magnetization(i, :)**2)
         difference_norm = difference_norm + sum(response_space%radial_weights* &
            (ground_state_magnetization(i, :) - pauli_response_magnetization(i, :))**2)
         magnetic_weight = 0.0_rp
         magnetic_last = response_space%npoint
         do ir = 1, response_space%npoint
            magnetic_weight = magnetic_weight + response_space%radial_weights(ir)*abs(ground_state_magnetization(i, ir))
            if (magnetic_weight >= 0.9_rp*sum(response_space%radial_weights*abs(ground_state_magnetization(i, :)))) then
               magnetic_last = ir
               exit
            end if
         end do
         magnetic_difference = magnetic_difference + sum(response_space%radial_weights(:magnetic_last)* &
            (ground_state_magnetization(i, :magnetic_last) - pauli_response_magnetization(i, :magnetic_last))**2)
         magnetic_norm = magnetic_norm + sum(response_space%radial_weights(:magnetic_last)* &
            ground_state_magnetization(i, :magnetic_last)**2)
      end do
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
         write(unit, '(a)') '# interaction_provenance = KXC-01 LR-METHOD-05R accepted SR direct ALSDA; compact projection U^H K_point U'
         write(unit, '(a,a)') '# ground_state_magnetization_kind = ', trim(kxc_result%magnetization_kind)
         write(unit, '(a,a)') '# ground_state_magnetization_source = ', trim(kxc_result%magnetization_source)
         write(unit, '(a,3(es24.16,1x))') '# integrated_M_SR_M_P_Delta_M = ', &
            4.0_rp*response_angular_pi*integrated_sr, 4.0_rp*response_angular_pi*integrated_pauli, &
            4.0_rp*response_angular_pi*(integrated_sr-integrated_pauli)
         write(unit, '(a,es24.16)') '# relative_L2_SR_minus_Pauli = ', sqrt(difference_norm/sr_norm)
         write(unit, '(a,es24.16)') '# magnetic_region_relative_L2 = ', sqrt(magnetic_difference/magnetic_norm)
         write(unit, '(a)') '# magnetic_region = first radius enclosing 90 percent of integrated abs(m_SR)'
         write(unit, '(a,2(es24.16,1x))') '# Kxc_SR_min_max = ', &
            minval(kxc_result%pointwise_kernel), maxval(kxc_result%pointwise_kernel)
         write(unit, '(a,a)') '# pauli_response_magnetization_source = ', &
            trim(kxc_result%pauli_response_magnetization_source)
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
            write(unit, '(a)') '# static_denominator = finite-eta diagnostic; eta ladder is not the exact static identity'
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
         do i = 1, size(static_eta)
            bare_request%eta = static_eta(i)
            call evaluate_lr_product_ks_susceptibility(bare_request, bare_result)
            ward_response = matmul(bare_result%susceptibility(:, :, 1), ward_field)
            ward_residual = ward_response - ward_target
            ward_abs = sqrt(sum(abs(ward_residual)**2))
            ward_rel = ward_abs/sqrt(sum(abs(ward_target)**2))
            ward_overlap = dot_product(ward_target, ward_response)/sum(abs(ward_target)**2)
            if (rank == 0) write(unit, '(a,6(es24.16,1x))') '# raw_Ward_eta_abs_rel_inf_overlap = ', &
               static_eta(i), ward_abs, ward_rel, maxval(abs(ward_residual)), real(ward_overlap, rp), aimag(ward_overlap)
            dyson_request%response_space => response_space
            dyson_request%compact_orthonormal = .true.
            dyson_request%q = config%q_list(:, gamma_index)
            dyson_request%frequencies = static_frequency
            dyson_request%eta = static_eta(i)
            dyson_request%channel = config%channel
            dyson_request%ks_susceptibility = bare_result%susceptibility
            dyson_request%canonical_interaction = interaction
            dyson_request%interaction_route = tddft_driver_route_direct_alsda
            dyson_request%interaction_provenance = 'KXC-01 LR-METHOD-05R accepted SR direct ALSDA'
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
         if (rank == 0) write(unit, '(a)') '# Ward action = raw chiKS times independently projected accepted SR bxc; no correction'
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
         dyson_request%interaction_provenance = 'KXC-01 LR-METHOD-05R accepted SR direct ALSDA'
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
            opposite_dyson_request%interaction_provenance = 'KXC-01 LR-METHOD-05R accepted SR direct ALSDA'
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


   pure function fold_fractional_kpoint(k_point) result(folded)
      real(rp), intent(in) :: k_point(3)
      real(rp) :: folded(3)

      folded = k_point - floor(k_point + 0.5_rp)
   end function fold_fractional_kpoint

   !> TDVK-02R3 accepted-Fe compact reciprocal-GF smoke.  This is deliberately
   !> separate from the point-grid reciprocal-GF service and is not a
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

   !> Assemble the SR functional contract and separate accepted Pauli response provenance.
   subroutine prepare_direct_alsda_request(response_space, ground_states, pauli_magnetization, request)
      type(response_space_layout), target, intent(in) :: response_space
      type(radial_ground_state), target, intent(in) :: ground_states(:)
      real(rp), intent(in), optional :: pauli_magnetization(:, :)
      type(lr_alsda_kernel_request), intent(out) :: request
      request%response_space => response_space
      request%ground_states => ground_states
      if (present(pauli_magnetization)) then
         allocate(request%pauli_magnetization(size(pauli_magnetization, 1), size(pauli_magnetization, 2)))
         request%pauli_magnetization = pauli_magnetization
         request%magnetization_label = 'pauli_projected'
         request%magnetization_kind = lr_kxc_magnetization_kind_pauli_accepted
         request%magnetization_source = lr_kxc_magnetization_source_pauli_accepted
      end if
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
      type(lr_ks_susceptibility_result) :: bare_result
      type(tddft_dyson_request) :: dyson_request
      type(tddft_dyson_result) :: dyson_result
      complex(rp), allocatable :: interaction(:, :)
      integer :: iq, ndim

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
         ! Pauli density is optional response metadata, never the kernel denominator.
         if (present(accepted_pauli_magnetization)) then
            if (allocated(accepted_pauli_magnetization)) then
               if (.not. present(accepted_pauli_magnetization_source)) &
                  error stop 'TDDFT response diagnostic: missing accepted Pauli provenance'
               if (trim(accepted_pauli_magnetization_source) /= lr_kxc_magnetization_source_pauli_accepted) &
                  error stop 'TDDFT response diagnostic: uncertified accepted Pauli provenance'
            end if
         end if
         call prepare_direct_alsda_request(response_space, ground_states, accepted_pauli_magnetization, kxc_request)
         call evaluate_lr_alsda_kernel(kxc_request, kxc_result)
         allocate(interaction(ndim, ndim))
         interaction = kxc_result%canonical_operator
         result%interaction_provenance = 'KXC-01 LR-METHOD-05R accepted SR direct ALSDA'
         result%magnetization_kind = kxc_result%magnetization_kind
         result%magnetization_source = kxc_result%magnetization_source
         result%goldstone_correction_status = 'disabled'
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

   subroutine radial_basis_from_snapshot(state, basis, complete_sr, requester)
      type(radial_ground_state), intent(in) :: state
      type(lmto_radial_basis), intent(out) :: basis
      logical, intent(in) :: complete_sr
      character(len=*), intent(in) :: requester
      logical :: use_complete_sr

      use_complete_sr = complete_sr

      if (use_complete_sr .and. (.not. state%pauli_basis_valid .or. .not. allocated(state%pauli_large) .or. &
          .not. allocated(state%pauli_large_dot) .or. .not. allocated(state%pauli_small) .or. &
          .not. allocated(state%pauli_enu) .or. .not. state%sr_basis_valid .or. &
          .not. allocated(state%sr_large_dot) .or. &
          .not. allocated(state%sr_small_dot) .or. .not. allocated(state%sr_large_ddot) .or. &
          .not. allocated(state%sr_small_ddot) .or. .not. allocated(state%sr_gfac) .or. &
          .not. allocated(state%sr_tmc) .or. .not. allocated(state%sr_potential) .or. &
          .not. allocated(state%sr_enu))) then
         error stop 'TDDFT production driver: '//trim(requester)// &
            ' requires the complete accepted SR radial basis, including phiddot'
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
         basis%phidot_large = state%sr_large_dot
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
      write(unit, '(a,3(es24.16,1x))') '# reciprocal_target_delta_electron_count = ', &
         reciprocal_obj%canonical_electron_count, reciprocal_obj%total_electrons, &
         reciprocal_obj%canonical_electron_count-reciprocal_obj%total_electrons
      write(unit, '(a)') '# response_capability = collinear,no_soc,ham_only,orthogonal,sp/spd,no_extra_operator'
      write(unit, '(a,a)') '# q = ', 'one row block per requested q; values are reduced coordinates'
      write(unit, '(a,a)') '# omega_Ry = ', 'one row per requested frequency'
      write(unit, '(a,es24.16)') '# eta_Ry = ', config%eta
      write(unit, '(a,i0)') '# response_angular_cutoff = ', result%response_lmax
      if (direct_handoff) then
         write(unit, '(a)') '# state_source = accepted_kspace_scf_cache'
         write(unit, '(a,l1)') '# direct_accepted_state_handoff = ', .true.
      else
         write(unit, '(a)') '# state_source = accepted_realspace_potential_adapter'
         write(unit, '(a,l1)') '# direct_accepted_state_handoff = ', .false.
         write(unit, '(a)') '# route = accepted real-space SCF potential adapter; accepted EF retained'
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
