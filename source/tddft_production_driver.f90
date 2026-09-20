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
module tddft_production_driver_mod

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
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex
   use lr_response_space_mod, only: response_space_layout, response_operator_trace
   use lr_ks_susceptibility_mod, only: lr_electronic_state, lr_ks_susceptibility_request, &
      lr_ks_susceptibility_result, lr_snapshot_from_reciprocal, lr_q_endpoint_from_reciprocal, &
      evaluate_lr_ks_susceptibility, lr_product_ks_susceptibility_request, &
      lr_product_ks_susceptibility_result, evaluate_lr_product_ks_susceptibility, lr_channel_plus, lr_channel_minus, &
      lr_fermi_dirac_occupation
      ! `lr_fermi_dirac_occupation` is used only to audit that the immutable
      ! TDDFT snapshot reproduces the reciprocal SCF occupation semantics.
   use lr_pauli_transition_vertex_mod, only: pauli_endpoint_state, pauli_vertex_capabilities, &
      pauli_sigma_plus_matrix, pauli_sigma_minus_matrix, evaluate_pauli_transition_vertex
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, lmto_product_channel_plus, &
      lmto_product_channel_minus
   use lr_compact_static_interaction_mod, only: lr_compact_gsr_result, lr_compact_representation, lr_compact_mapping_contract, &
      compact_project_magnetization, compact_project_local_operator, compact_apply_local_operator, &
      compact_reconstruct_point_vector, compact_weighted_projection_diagnostics, evaluate_compact_goldstone_sumrule, &
      lr_compact_gsr_action_tolerance
   use pauli_ground_state_projection_mod, only: compute_accepted_pauli_magnetization
   use lr_gf_susceptibility_mod, only: lr_gf_susceptibility_request, evaluate_lr_gf_susceptibility
   use lr_product_gf_susceptibility_mod, only: lr_product_gf_susceptibility_request, &
      lr_product_gf_susceptibility_result, evaluate_lr_product_gf_susceptibility
   use lr_projected_site_spin_mod, only: projected_site_spin_contract
   use lr_projected_reciprocal_chi0_mod, only: projected_chi0_request, projected_chi0_result, &
      evaluate_projected_lehmann_chi0, evaluate_projected_finite_width_chi0, evaluate_projected_gf_chi0
   use lr_projected_interacting_response_mod, only: projected_mills_interaction_result, &
      projected_dyson_request, projected_dyson_result, evaluate_projected_mills_from_reciprocal, &
      evaluate_projected_dyson
   use lr_projected_juelich_interaction_mod, only: projected_juelich_request, projected_juelich_result, &
      evaluate_projected_juelich_interaction, evaluate_projected_juelich_holdout, &
      assess_projected_juelich_eta_stability, projected_juelich_eta_limited, projected_juelich_rank_deficient, &
      projected_juelich_unsupported, select_projected_juelich_eta_indices
   use lr_dresp06a_bridge_mod, only: project_compact_response_to_sites, project_compact_loss_to_sites, &
      compact_loss_matrix, site_loss_matrix, compact_apply_local_operator_independent, dresp06a_relative_matrix_difference, &
      dresp06a_projection_tolerance, build_compact_covariance_transport, transport_compact_matrix
   use lr_dresp07_bridge_mod, only: run_dresp07_exact_ks_ward
   use lr_dresp08_bridge_mod, only: run_dresp08_native_second_order
   use lr_dresp09r_bridge_mod, only: run_dresp09r_pauli_projected_native
   use lr_dresp09s_bridge_mod, only: run_dresp09s_scalar_relativistic
   use lr_dresp09t_bridge_mod, only: run_dresp09t_moving_basis
   use lr_dresp09u_bridge_mod, only: run_dresp09u_representation_tangent
   use lr_dresp09v_bridge_mod, only: run_dresp09v_density_moment_tangent
   use lr_dresp09w_bridge_mod, only: run_dresp09w_radial_observable_provenance
   use lr_ward_mode_analysis_mod, only: lr_ward_mode_analysis_result, analyze_lr_ward_mode
   use lr_rs_gf_susceptibility_mod, only: lr_rs_gf_provider, lr_rs_gf_pair, lr_rs_gf_susceptibility_request, &
      evaluate_lr_rs_gf_susceptibility
   use tddft_native_rsgf_provider_mod, only: tddft_native_rsgf_provider
   use lr_alsda_kernel_mod, only: lr_alsda_kernel_request, lr_alsda_kernel_result, evaluate_lr_alsda_kernel, &
      lr_kxc_magnetization_kind_pauli_accepted, lr_kxc_magnetization_source_pauli_accepted
   use lr_goldstone_sumrule_mod, only: lr_goldstone_sumrule_request, lr_goldstone_sumrule_result, &
      evaluate_lr_goldstone_sumrule
   use tddft_dyson_mod, only: tddft_dyson_request, tddft_dyson_result, evaluate_tddft_dyson, &
      lr_dyson_route_direct_alsda, lr_dyson_route_goldstone_sumrule
   use logger_mod, only: g_logger
   implicit none
   private

#ifdef VERSION
   character(len=*), parameter :: tddft_build_version = VERSION
#else
   character(len=*), parameter :: tddft_build_version = 'unavailable'
#endif

   character(len=*), parameter, public :: tddft_driver_backend_lehmann = 'lehmann'
   character(len=*), parameter, public :: tddft_driver_backend_reciprocal_gf = 'reciprocal_gf'
   character(len=*), parameter, public :: tddft_driver_backend_native_rsgf = 'native_rsgf'
   character(len=*), parameter, public :: tddft_driver_backend_product_lehmann = 'product_lehmann'
   character(len=*), parameter, public :: tddft_driver_backend_product_gf = 'product_gf'
   character(len=*), parameter, public :: tddft_driver_backend_product_finite_q = 'product_finite_q'
   character(len=*), parameter, public :: tddft_driver_backend_product_convergence = 'product_convergence'
   character(len=*), parameter, public :: tddft_driver_backend_projected_chi0 = 'projected_chi0'
   character(len=*), parameter, public :: tddft_driver_backend_static_interactions = 'static_interactions'
   character(len=*), parameter, public :: tddft_driver_backend_compact_dyson = 'compact_dyson'
   character(len=*), parameter, public :: tddft_driver_backend_projected_mills = 'projected_mills'
   character(len=*), parameter, public :: tddft_driver_backend_projected_juelich = 'projected_juelich'
   character(len=*), parameter, public :: tddft_driver_backend_alsda_compare = 'alsda_compare'
   character(len=*), parameter, public :: tddft_driver_backend_exact_ks_ward = 'exact_ks_ward'
   character(len=*), parameter, public :: tddft_driver_backend_native_second_order = 'native_second_order'
   character(len=*), parameter, public :: tddft_driver_backend_pauli_projected_native = 'pauli_projected_native'
   character(len=*), parameter, public :: tddft_driver_backend_sr_l0_scalar_relativistic = 'sr_l0_scalar_relativistic'
   character(len=*), parameter, public :: tddft_driver_backend_moving_basis_tangent = 'moving_basis_tangent'
   character(len=*), parameter, public :: tddft_driver_backend_representation_tangent = 'representation_tangent'
   character(len=*), parameter, public :: tddft_driver_backend_density_moment_tangent = 'density_moment_tangent'
   character(len=*), parameter, public :: tddft_driver_backend_radial_observable_provenance = 'radial_observable_provenance'
   character(len=*), parameter, public :: tddft_driver_route_direct_alsda = lr_dyson_route_direct_alsda
   character(len=*), parameter, public :: tddft_driver_route_goldstone_sumrule = lr_dyson_route_goldstone_sumrule

   !> Parsed, minimal production input. Frequencies and q points are expanded
   !> once at the parser boundary so the production path never reparses input.
   type, public :: tddft_production_config
      logical :: present = .false.
      logical :: enabled = .false.
      integer :: nq = 1
      integer :: nfrequency = 1
      character(len=32) :: channel = 'chi_plus'
      real(rp), allocatable :: q_list(:, :) ! (3,nq)
      real(rp), allocatable :: frequencies(:)
      real(rp), allocatable :: eta_values(:) ! physical response eta values for validation campaigns
      real(rp) :: eta = 0.01_rp
      integer :: response_lmax = -1
      character(len=48) :: interaction_route = tddft_driver_route_direct_alsda
      logical :: goldstone_correction = .false.
      character(len=32) :: backend = tddft_driver_backend_lehmann
      character(len=8) :: projected_selector = 'spd'
      logical :: reciprocal_backend_crosscheck = .false.
      character(len=32) :: native_rsgf_provider = 'auto'
      integer :: gf_integration_points = 2001
      real(rp) :: gf_integration_eta = 0.0_rp
      real(rp) :: gf_energy_margin = 1.0_rp
      logical :: gf_closure_audit = .false.
      logical :: dyson_static_audit = .false.
      logical :: validate_interacting_covariance = .false.
      logical :: write_full_matrix = .true.
      character(len=256) :: output_file = 'tddft_response.dat'
   contains
      procedure :: restore_to_default => tddft_config_restore
   end type tddft_production_config

   !> Small capability contract used both by the real-object gate and focused
   !> rejection tests.  It makes the reason for a rejection observable without
   !> requiring a material fixture.
   type, public :: tddft_capability_state
      character(len=32) :: reciprocal_mode = 'ham_only'
      character(len=16) :: hamiltonian_order = 'second'
      integer :: basis_lmax = 2
      logical :: collinear = .true.
      logical :: has_soc = .false.
      logical :: orthogonal = .true.
      logical :: generalized_overlap = .false.
      logical :: has_extra_operator = .false.
      logical :: bulk_cell = .true.
   end type tddft_capability_state

   !> Driver result. The complete matrices are retained; no spectrum or mode
   !> reduction is invented at this orchestration layer.
   type, public :: tddft_production_result
      logical :: initialized = .false.
      integer :: nq = 0
      integer :: nfrequency = 0
      integer :: ndim = 0
      integer :: response_lmax = -1
      real(rp), allocatable :: q_list(:, :)
      real(rp), allocatable :: frequencies(:)
      complex(rp), allocatable :: ks_susceptibility(:, :, :, :)
      complex(rp), allocatable :: enhanced_susceptibility(:, :, :, :)
      complex(rp), allocatable :: loss_matrix(:, :, :, :)
      logical :: compact_orthonormal = .false.
      logical :: reciprocal_backend_crosscheck = .false.
      logical, allocatable :: reciprocal_crosscheck_valid(:, :) ! (frequency,q)
      real(rp), allocatable :: reciprocal_crosscheck_norm_lehmann(:, :) ! (frequency,q)
      real(rp), allocatable :: reciprocal_crosscheck_norm_gf(:, :) ! (frequency,q)
      real(rp), allocatable :: reciprocal_crosscheck_difference_frobenius(:, :) ! (frequency,q)
      real(rp), allocatable :: reciprocal_crosscheck_relative_frobenius(:, :) ! (frequency,q)
      real(rp), allocatable :: reciprocal_crosscheck_difference_infinity(:, :) ! (frequency,q)
      complex(rp), allocatable :: reciprocal_crosscheck_delta(:, :, :, :) ! (I,J,frequency,q)
      character(len=128) :: status = 'not evaluated'
      character(len=48) :: interaction_route = ''
      character(len=32) :: backend = ''
      character(len=64) :: goldstone_correction_status = 'not selected'
      character(len=256) :: interaction_provenance = ''
      character(len=256) :: bare_response_provenance = ''
      character(len=32) :: magnetization_kind = ''
      character(len=256) :: magnetization_source = ''
   end type tddft_production_result

   public :: load_tddft_config
   public :: tddft_capability_is_supported
   public :: require_tddft_capability
   public :: validate_tddft_production_capability
   public :: run_tddft_production
   public :: evaluate_tddft_production_sweep

contains

   subroutine tddft_config_restore(this)
      class(tddft_production_config), intent(out) :: this
      this%present = .false.
      this%enabled = .false.
      this%nq = 1
      this%nfrequency = 1
      this%channel = 'chi_plus'
      this%eta = 0.01_rp
      this%response_lmax = -1
      this%interaction_route = tddft_driver_route_direct_alsda
      this%goldstone_correction = .false.
      this%backend = tddft_driver_backend_lehmann
      this%projected_selector = 'spd'
      this%reciprocal_backend_crosscheck = .false.
      this%native_rsgf_provider = 'auto'
      this%gf_integration_points = 2001
      this%gf_integration_eta = 0.0_rp
      this%gf_energy_margin = 1.0_rp
      this%gf_closure_audit = .false.
      this%dyson_static_audit = .false.
      this%validate_interacting_covariance = .false.
      this%write_full_matrix = .true.
      this%output_file = 'tddft_response.dat'
      if (allocated(this%q_list)) deallocate(this%q_list)
      if (allocated(this%frequencies)) deallocate(this%frequencies)
      if (allocated(this%eta_values)) deallocate(this%eta_values)
      allocate(this%q_list(3, 1), this%frequencies(1), this%eta_values(1))
      this%q_list = 0.0_rp
      this%frequencies = 0.0_rp
      this%eta_values = this%eta
   end subroutine tddft_config_restore

   !> Read only the new minimal &tddft group. No legacy TD-DFT type is
   !> constructed and an absent group is the ordinary feature-off default.
   subroutine load_tddft_config(filename, config)
      character(len=*), intent(in) :: filename
      type(tddft_production_config), intent(out) :: config
      integer :: unit, ios, n_q_local = 1, n_omega_local = 1, n_eta_local = 1, iq, iw
      logical :: found
      character(len=512) :: line

      include 'include_codes/namelists/tddft.f90'

      call config%restore_to_default()

      enabled = .false.
      channel = 'chi_plus'
      projected_selector = 'spd'
      n_q = 1
      q_list = 0.0_rp
      n_omega = 1
      use_omega_grid = .false.
      omega_grid = 0.0_rp
      omega_min = 0.0_rp
      omega_max = 0.0_rp
      eta = 0.01_rp
      n_eta = 1
      eta_grid = 0.0_rp
      eta_grid(1) = eta
      response_lmax = -1
      interaction_route = tddft_driver_route_direct_alsda
      goldstone_correction = .false.
      backend = tddft_driver_backend_lehmann
      reciprocal_backend_crosscheck = .false.
      native_rsgf_provider = 'auto'
      gf_integration_points = 2001
      gf_integration_eta = 0.0_rp
      gf_energy_margin = 1.0_rp
      gf_closure_audit = .false.
      dyson_static_audit = .false.
      validate_interacting_covariance = .false.
      write_full_matrix = .true.
      output_file = 'tddft_response.dat'

      found = .false.
      open(newunit=unit, file=filename, action='read', status='old', iostat=ios)
      if (ios /= 0) return
      do
         read(unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         if (index(adjustl(lower(line)), '&tddft') == 1) then
            found = .true.
            exit
         end if
      end do
      close(unit)
      if (.not. found) return

      open(newunit=unit, file=filename, action='read', status='old', iostat=ios)
      if (ios /= 0) error stop 'TDDFT parser: input file cannot be reopened'
      read(unit, nml=tddft, iostat=ios)
      close(unit)
      if (ios /= 0) error stop 'TDDFT parser: invalid &tddft namelist'

      config%present = .true.
      config%enabled = enabled
      config%channel = trim(lower(channel))
      config%eta = eta
      config%response_lmax = response_lmax
      config%interaction_route = trim(lower(interaction_route))
      config%goldstone_correction = goldstone_correction
      config%backend = trim(lower(backend))
      config%projected_selector = trim(lower(projected_selector))
      config%reciprocal_backend_crosscheck = reciprocal_backend_crosscheck
      config%native_rsgf_provider = trim(lower(native_rsgf_provider))
      config%gf_integration_points = gf_integration_points
      config%gf_integration_eta = gf_integration_eta
      config%gf_energy_margin = gf_energy_margin
      config%gf_closure_audit = gf_closure_audit
      config%dyson_static_audit = dyson_static_audit
      config%validate_interacting_covariance = validate_interacting_covariance
      config%write_full_matrix = write_full_matrix
      config%output_file = trim(output_file)

      n_q_local = n_q
      n_omega_local = n_omega
      n_eta_local = n_eta
      if (n_q_local < 1 .or. n_q_local > tddft_max_q) error stop 'TDDFT parser: n_q is outside [1,tddft_max_q]'
      if (n_omega_local < 1 .or. n_omega_local > tddft_max_omega) then
         error stop 'TDDFT parser: n_omega is outside [1,tddft_max_omega]'
      end if
      if (n_eta_local < 1 .or. n_eta_local > 16) error stop 'TDDFT parser: n_eta is outside [1,16]'
      if (allocated(config%q_list)) deallocate(config%q_list)
      if (allocated(config%frequencies)) deallocate(config%frequencies)
      if (allocated(config%eta_values)) deallocate(config%eta_values)
      allocate(config%q_list(3, n_q_local), config%frequencies(n_omega_local), config%eta_values(n_eta_local))
      config%nq = n_q_local
      config%nfrequency = n_omega_local
      config%q_list = q_list(:, 1:n_q_local)
      if (use_omega_grid) then
         config%frequencies = omega_grid(1:n_omega_local)
      else if (n_omega_local == 1) then
         config%frequencies(1) = omega_min
      else
         do iw = 1, n_omega_local
            config%frequencies(iw) = omega_min + real(iw - 1, rp)*(omega_max - omega_min)/real(n_omega_local - 1, rp)
         end do
      end if
      if (n_eta_local == 1) then
         config%eta_values(1) = eta
      else
         config%eta_values = eta_grid(1:n_eta_local)
         config%eta = config%eta_values(1)
      end if
      do iq = 1, n_q_local
         if (any(config%q_list(:, iq) /= config%q_list(:, iq))) error stop 'TDDFT parser: q contains NaN'
      end do
      call validate_tddft_config(config)
   end subroutine load_tddft_config

   subroutine validate_tddft_config(config)
      type(tddft_production_config), intent(in) :: config
      character(len=128) :: route
      integer :: eta_selected_index, eta_holdout_index

      if (.not. config%enabled) return
      if (trim(config%channel) /= 'chi_plus' .and. trim(config%channel) /= 'chi_minus') then
         error stop 'TDDFT input: unsupported response channel; use chi_plus or chi_minus'
      end if
      if (config%eta <= 0.0_rp) error stop 'TDDFT input: eta must be positive'
      if (.not. allocated(config%eta_values) .or. size(config%eta_values) < 1 .or. &
          any(config%eta_values <= 0.0_rp)) then
         error stop 'TDDFT input: every physical response eta must be positive'
      end if
      if (config%response_lmax < -1) error stop 'TDDFT input: response_lmax must be -1 or a nonnegative cutoff'
      if (config%response_lmax > 4) error stop 'TDDFT input: response_lmax above 4 is outside the validated sp/spd product baseline'
      route = trim(config%interaction_route)
      if (route /= tddft_driver_route_direct_alsda .and. route /= tddft_driver_route_goldstone_sumrule) then
         error stop 'TDDFT input: unsupported interaction_route; use direct_alsda or goldstone_sumrule'
      end if
      if (config%goldstone_correction) then
         error stop 'TDDFT input: Goldstone correction requires separate validated TDVAL evidence; it is not silently applied'
      end if
      if (trim(config%backend) /= tddft_driver_backend_lehmann .and. trim(config%backend) /= 'spectral' .and. &
          trim(config%backend) /= tddft_driver_backend_reciprocal_gf .and. &
          trim(config%backend) /= tddft_driver_backend_native_rsgf .and. &
          trim(config%backend) /= tddft_driver_backend_product_lehmann .and. &
          trim(config%backend) /= tddft_driver_backend_product_gf .and. &
          trim(config%backend) /= tddft_driver_backend_product_finite_q .and. &
          trim(config%backend) /= tddft_driver_backend_product_convergence .and. &
          trim(config%backend) /= tddft_driver_backend_projected_chi0 .and. &
          trim(config%backend) /= tddft_driver_backend_static_interactions .and. &
          trim(config%backend) /= tddft_driver_backend_compact_dyson .and. &
          trim(config%backend) /= tddft_driver_backend_projected_mills .and. &
          trim(config%backend) /= tddft_driver_backend_projected_juelich .and. &
          trim(config%backend) /= tddft_driver_backend_alsda_compare .and. &
          trim(config%backend) /= tddft_driver_backend_exact_ks_ward .and. &
          trim(config%backend) /= tddft_driver_backend_native_second_order .and. &
          trim(config%backend) /= tddft_driver_backend_pauli_projected_native .and. &
          trim(config%backend) /= tddft_driver_backend_sr_l0_scalar_relativistic .and. &
          trim(config%backend) /= tddft_driver_backend_moving_basis_tangent .and. &
          trim(config%backend) /= tddft_driver_backend_representation_tangent .and. &
          trim(config%backend) /= tddft_driver_backend_density_moment_tangent .and. &
          trim(config%backend) /= tddft_driver_backend_radial_observable_provenance) then
         error stop 'TDDFT input: unsupported backend; use spectral/lehmann, reciprocal_gf, native_rsgf, product_lehmann, product_gf, product_finite_q, product_convergence, projected_chi0, static_interactions, compact_dyson, projected_mills, projected_juelich, alsda_compare, exact_ks_ward, native_second_order, pauli_projected_native, sr_l0_scalar_relativistic, moving_basis_tangent, representation_tangent, density_moment_tangent, or radial_observable_provenance'
      end if
      if (trim(config%backend) == tddft_driver_backend_projected_mills .or. &
          trim(config%backend) == tddft_driver_backend_projected_juelich .or. &
          trim(config%backend) == tddft_driver_backend_alsda_compare) then
         if (trim(config%projected_selector) /= 'd' .and. trim(config%projected_selector) /= 'spd' .and. &
             trim(config%projected_selector) /= 'both') then
            error stop 'TDDFT input: projected_selector must be d, spd, or both'
         end if
         if (config%response_lmax >= 0 .and. config%response_lmax /= 4) then
            error stop 'TDDFT input: projected Mills/Juelich requires the complete DRESP-01 response_lmax=4 contract'
         end if
         if (trim(config%backend) == tddft_driver_backend_projected_juelich .and. &
             .not. any(sum(abs(config%q_list), dim=1) <= 1.0e-12_rp)) then
            error stop 'TDDFT input: projected_juelich requires Gamma for its static Ward construction'
         end if
         if (trim(config%backend) == tddft_driver_backend_alsda_compare .and. &
             .not. any(sum(abs(config%q_list), dim=1) <= 1.0e-12_rp)) then
            error stop 'TDDFT input: alsda_compare requires Gamma for the Juelich same-state construction'
         end if
      end if
      if (trim(config%backend) /= tddft_driver_backend_product_convergence .and. &
          trim(config%backend) /= tddft_driver_backend_static_interactions .and. &
          trim(config%backend) /= tddft_driver_backend_projected_mills .and. &
          trim(config%backend) /= tddft_driver_backend_projected_juelich .and. &
          trim(config%backend) /= tddft_driver_backend_alsda_compare .and. &
          trim(config%backend) /= tddft_driver_backend_exact_ks_ward .and. &
          trim(config%backend) /= tddft_driver_backend_native_second_order .and. &
          trim(config%backend) /= tddft_driver_backend_pauli_projected_native .and. &
          trim(config%backend) /= tddft_driver_backend_sr_l0_scalar_relativistic .and. &
          trim(config%backend) /= tddft_driver_backend_moving_basis_tangent .and. &
          trim(config%backend) /= tddft_driver_backend_representation_tangent .and. &
          trim(config%backend) /= tddft_driver_backend_density_moment_tangent .and. &
          trim(config%backend) /= tddft_driver_backend_radial_observable_provenance .and. size(config%eta_values) /= 1) then
         error stop 'TDDFT input: n_eta greater than one is only supported by product_convergence, static_interactions, projected_mills, projected_juelich, alsda_compare, exact_ks_ward, native_second_order, pauli_projected_native, sr_l0_scalar_relativistic, moving_basis_tangent, representation_tangent, density_moment_tangent, or radial_observable_provenance'
      end if
      if (trim(config%backend) == tddft_driver_backend_exact_ks_ward) then
         if (size(config%q_list, 2) /= 1 .or. sum(abs(config%q_list(:, 1))) > 1.0e-12_rp) then
            error stop 'TDDFT input: exact_ks_ward requires one Gamma q point'
         end if
         if (size(config%frequencies) /= 1 .or. abs(config%frequencies(1)) > 1.0e-12_rp) then
            error stop 'TDDFT input: exact_ks_ward requires one static omega=0 point'
         end if
         if (config%response_lmax >= 0 .and. config%response_lmax /= 4) then
            error stop 'TDDFT input: exact_ks_ward requires the complete response_lmax=4 product space'
         end if
      end if
      if (trim(config%backend) == tddft_driver_backend_native_second_order) then
         if (size(config%q_list, 2) /= 1 .or. sum(abs(config%q_list(:, 1))) > 1.0e-12_rp) then
            error stop 'TDDFT input: native_second_order requires one Gamma q point'
         end if
         if (size(config%frequencies) /= 1 .or. abs(config%frequencies(1)) > 1.0e-12_rp) then
            error stop 'TDDFT input: native_second_order requires one static omega=0 point'
         end if
         if (config%response_lmax >= 0 .and. config%response_lmax /= 4) then
            error stop 'TDDFT input: native_second_order requires the complete response_lmax=4 product space'
         end if
      end if
      if (trim(config%backend) == tddft_driver_backend_pauli_projected_native) then
         if (size(config%q_list, 2) /= 1 .or. sum(abs(config%q_list(:, 1))) > 1.0e-12_rp) then
            error stop 'TDDFT input: pauli_projected_native requires one Gamma q point'
         end if
         if (size(config%frequencies) /= 1 .or. abs(config%frequencies(1)) > 1.0e-12_rp) then
            error stop 'TDDFT input: pauli_projected_native requires one static omega=0 point'
         end if
         if (config%response_lmax >= 0 .and. config%response_lmax /= 4) then
            error stop 'TDDFT input: pauli_projected_native requires the complete response_lmax=4 product space'
         end if
      end if
      if (trim(config%backend) == tddft_driver_backend_sr_l0_scalar_relativistic) then
         if (size(config%q_list, 2) /= 1 .or. sum(abs(config%q_list(:, 1))) > 1.0e-12_rp) then
            error stop 'TDDFT input: sr_l0_scalar_relativistic requires one Gamma q point'
         end if
         if (size(config%frequencies) /= 1 .or. abs(config%frequencies(1)) > 1.0e-12_rp) then
            error stop 'TDDFT input: sr_l0_scalar_relativistic requires one static omega=0 point'
         end if
         if (config%response_lmax /= 4) then
            error stop 'TDDFT input: sr_l0_scalar_relativistic requires response_lmax=4'
         end if
      end if
      if (trim(config%backend) == tddft_driver_backend_moving_basis_tangent) then
         if (size(config%q_list, 2) /= 1 .or. sum(abs(config%q_list(:, 1))) > 1.0e-12_rp) then
            error stop 'TDDFT input: moving_basis_tangent requires one Gamma q point'
         end if
         if (size(config%frequencies) /= 1 .or. abs(config%frequencies(1)) > 1.0e-12_rp) then
            error stop 'TDDFT input: moving_basis_tangent requires one static omega=0 point'
         end if
         if (config%response_lmax /= 4) then
            error stop 'TDDFT input: moving_basis_tangent requires response_lmax=4'
         end if
      end if
      if (trim(config%backend) == tddft_driver_backend_representation_tangent) then
         if (size(config%q_list, 2) /= 1 .or. sum(abs(config%q_list(:, 1))) > 1.0e-12_rp) then
            error stop 'TDDFT input: representation_tangent requires one Gamma q point'
         end if
         if (size(config%frequencies) /= 1 .or. abs(config%frequencies(1)) > 1.0e-12_rp) then
            error stop 'TDDFT input: representation_tangent requires one static omega=0 point'
         end if
         if (config%response_lmax /= 4) then
            error stop 'TDDFT input: representation_tangent requires response_lmax=4'
         end if
      end if
      if (trim(config%backend) == tddft_driver_backend_density_moment_tangent) then
         if (size(config%q_list, 2) /= 1 .or. sum(abs(config%q_list(:, 1))) > 1.0e-12_rp) then
            error stop 'TDDFT input: density_moment_tangent requires one Gamma q point'
         end if
         if (size(config%frequencies) /= 1 .or. abs(config%frequencies(1)) > 1.0e-12_rp) then
            error stop 'TDDFT input: density_moment_tangent requires one static omega=0 point'
         end if
         if (config%response_lmax /= 4) then
            error stop 'TDDFT input: density_moment_tangent requires response_lmax=4'
         end if
      end if
      if (trim(config%backend) == tddft_driver_backend_radial_observable_provenance) then
         if (size(config%q_list, 2) /= 1 .or. sum(abs(config%q_list(:, 1))) > 1.0e-12_rp) then
            error stop 'TDDFT input: radial_observable_provenance requires one Gamma q point'
         end if
         if (size(config%frequencies) /= 1 .or. abs(config%frequencies(1)) > 1.0e-12_rp) then
            error stop 'TDDFT input: radial_observable_provenance requires one static omega=0 point'
         end if
         if (config%response_lmax /= 4) then
            error stop 'TDDFT input: radial_observable_provenance requires response_lmax=4'
         end if
      end if
      if (trim(config%backend) == tddft_driver_backend_projected_juelich .or. &
          trim(config%backend) == tddft_driver_backend_alsda_compare) then
         call select_projected_juelich_eta_indices(config%eta_values, eta_selected_index, eta_holdout_index)
      end if
      if (trim(config%backend) == tddft_driver_backend_static_interactions) then
         if (size(config%q_list, 2) /= 1 .or. sum(abs(config%q_list(:, 1))) > 1.0e-12_rp) then
            error stop 'TDDFT input: static_interactions requires one Gamma q point'
         end if
         if (size(config%frequencies) /= 1 .or. abs(config%frequencies(1)) > 1.0e-12_rp) then
            error stop 'TDDFT input: static_interactions requires one static omega=0 point'
         end if
         if (size(config%eta_values) > 2) then
            error stop 'TDDFT input: static_interactions supports at most eta and optional eta=.005 Ry'
         end if
      end if
      if (trim(config%backend) == tddft_driver_backend_reciprocal_gf .or. &
          trim(config%backend) == tddft_driver_backend_native_rsgf .or. config%reciprocal_backend_crosscheck) then
         if (config%gf_integration_points < 3 .or. mod(config%gf_integration_points, 2) == 0) then
            error stop 'TDDFT input: reciprocal_gf or backend crosscheck requires an odd gf_integration_points value >= 3'
         end if
         if (config%gf_energy_margin <= 0.0_rp) error stop 'TDDFT input: reciprocal-GF crosscheck requires gf_energy_margin positive'
      end if
      if (trim(config%backend) == tddft_driver_backend_product_gf) then
         if (config%gf_integration_points < 3 .or. mod(config%gf_integration_points, 2) == 0) then
            error stop 'TDDFT input: product_gf requires an odd gf_integration_points value >= 3'
         end if
         if (config%gf_energy_margin <= 0.0_rp) error stop 'TDDFT input: product_gf requires gf_energy_margin positive'
      end if
      if (trim(config%backend) == tddft_driver_backend_projected_chi0) then
         if (config%response_lmax >= 0 .and. config%response_lmax /= 4) then
            error stop 'TDDFT input: projected_chi0 requires the complete response_lmax=4 product space'
         end if
         if (config%gf_integration_points < 3 .or. mod(config%gf_integration_points, 2) == 0) then
            error stop 'TDDFT input: projected_chi0 requires an odd gf_integration_points value >= 3'
         end if
         if (config%gf_integration_eta < 0.0_rp .or. config%gf_integration_eta >= config%eta) then
            error stop 'TDDFT input: projected_chi0 requires 0 <= gf_integration_eta < eta'
         end if
         if (config%gf_energy_margin <= 0.0_rp) then
            error stop 'TDDFT input: projected_chi0 requires gf_energy_margin positive'
         end if
      end if
      if (trim(config%backend) == tddft_driver_backend_product_finite_q) then
         if (size(config%q_list, 2) < 4) then
            error stop 'TDDFT input: product_finite_q requires Gamma, q, -q, and an arbitrary q'
         end if
         if (config%gf_integration_points < 3 .or. mod(config%gf_integration_points, 2) == 0) then
            error stop 'TDDFT input: product_finite_q requires an odd gf_integration_points value >= 3'
         end if
         if (config%gf_energy_margin <= 0.0_rp) error stop 'TDDFT input: product_finite_q requires gf_energy_margin positive'
      end if
      if (trim(config%backend) == tddft_driver_backend_product_convergence) then
         if (size(config%q_list, 2) < 1) error stop 'TDDFT input: product_convergence requires at least one q point'
         if (.not. any(sum(abs(config%q_list), dim=1) <= 1.0e-12_rp)) then
            error stop 'TDDFT input: product_convergence requires Gamma for the full-operator artifact'
         end if
         if (.not. any(abs(config%frequencies) <= 1.0e-12_rp)) then
            error stop 'TDDFT input: product_convergence requires omega=0 for the full-operator artifact'
         end if
         if (config%gf_closure_audit) then
            if (config%gf_integration_points < 3 .or. mod(config%gf_integration_points, 2) == 0) then
               error stop 'TDDFT input: product_convergence GF audit requires an odd gf_integration_points value >= 3'
            end if
            if (config%gf_integration_eta <= 0.0_rp .or. config%gf_integration_eta >= config%eta) then
               error stop 'TDDFT input: product_convergence GF audit requires 0 < gf_integration_eta < eta'
            end if
            if (config%gf_energy_margin <= 0.0_rp) then
               error stop 'TDDFT input: product_convergence GF audit requires gf_energy_margin positive'
            end if
         end if
      end if
      if (trim(config%backend) == tddft_driver_backend_native_rsgf) then
         if (trim(config%native_rsgf_provider) /= 'auto' .and. trim(config%native_rsgf_provider) /= 'block' .and. &
             trim(config%native_rsgf_provider) /= 'block_recursion' .and. trim(config%native_rsgf_provider) /= 'chebyshev') then
            error stop 'TDDFT input: native_rsgf_provider must be auto, block or chebyshev'
         end if
      end if
      if (trim(config%backend) == tddft_driver_backend_compact_dyson .or. &
          trim(config%backend) == tddft_driver_backend_alsda_compare) then
         if (route /= tddft_driver_route_direct_alsda) then
            error stop 'TDDFT input: compact_dyson/alsda_compare requires interaction_route=direct_alsda; rejected before SCF/response work'
         end if
         if (config%dyson_static_audit .and. find_gamma_q_index(config%q_list) == 0) then
            error stop 'TDDFT input: dyson_static_audit requires Gamma in q_list; rejected before SCF/response work'
         end if
         if (config%validate_interacting_covariance) then
            call validate_covariance_q_pair(config%q_list)
         end if
         if (config%gf_closure_audit) then
            if (find_gamma_q_index(config%q_list) == 0) then
               error stop 'TDDFT input: gf_closure_audit requires Gamma in q_list; rejected before SCF/response work'
            end if
            if (config%gf_integration_points < 3 .or. mod(config%gf_integration_points, 2) == 0) then
               error stop 'TDDFT input: compact GF audit requires an odd gf_integration_points value >= 3'
            end if
            if (config%gf_integration_eta <= 0.0_rp .or. config%gf_integration_eta >= config%eta) then
               error stop 'TDDFT input: compact GF audit requires 0 < gf_integration_eta < eta'
            end if
            if (config%gf_energy_margin <= 0.0_rp) then
               error stop 'TDDFT input: compact GF audit requires gf_energy_margin positive'
            end if
         end if
      else if (config%dyson_static_audit .or. config%validate_interacting_covariance) then
         error stop 'TDDFT input: compact Dyson validation switches require backend=compact_dyson'
      end if
      if (len_trim(config%output_file) == 0) error stop 'TDDFT input: output_file must not be empty'
      if (route == tddft_driver_route_goldstone_sumrule) then
         if (.not. any(sum(abs(config%q_list), dim=1) <= 1.0e-12_rp)) then
            error stop 'TDDFT input: goldstone_sumrule requires q=(0,0,0) in q_list for the static reference'
         end if
      end if
   end subroutine validate_tddft_config

   logical function tddft_capability_is_supported(state, reason) result(ok)
      type(tddft_capability_state), intent(in) :: state
      character(len=*), intent(out) :: reason

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

   subroutine require_tddft_capability(state)
      type(tddft_capability_state), intent(in) :: state
      character(len=256) :: reason
      if (.not. tddft_capability_is_supported(state, reason)) then
         error stop 'TDDFT capability gate: unsupported feature: '//trim(reason)
      end if
   end subroutine require_tddft_capability

   subroutine validate_tddft_production_capability(config, control_obj, lattice_obj, hamiltonian_obj, reciprocal_obj)
      type(tddft_production_config), intent(in) :: config
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
      if (size(reciprocal_obj%k_workset%points, 2) /= nk .or. &
          size(reciprocal_obj%k_workset%weights) /= nk .or. &
          any(shape(reciprocal_obj%eigenvalues) /= [nbands, nk]) .or. &
          any(shape(reciprocal_obj%eigenvectors) /= [nbasis, nbands, nk])) then
         error stop 'TDDFT state consistency: reciprocal and left-state shapes differ'
      end if
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
      type(tddft_production_config), intent(in) :: config
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
      need_complete_sr = trim(config%backend) == tddft_driver_backend_radial_observable_provenance

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
      if (trim(config%backend) == tddft_driver_backend_exact_ks_ward) then
         if (.not. use_accepted_kspace_scf) then
            error stop 'DRESP-07 exact_ks_ward requires the accepted k-space SCF handoff'
         end if
         call run_dresp07_exact_ks_ward(trim(config%output_file), response_space, radial_bases, ground_states, &
            reciprocal_obj, lattice_obj, config%eta_values)
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_radial_observable_provenance) then
         if (.not. use_accepted_kspace_scf) then
            error stop 'DRESP-09W radial_observable_provenance requires the accepted k-space SCF handoff'
         end if
         call run_dresp09w_radial_observable_provenance(trim(config%output_file), response_space, radial_bases, ground_states, &
            reciprocal_obj, lattice_obj, hamiltonian_obj)
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_moving_basis_tangent) then
         if (.not. use_accepted_kspace_scf) then
            error stop 'DRESP-09T moving_basis_tangent requires the accepted k-space SCF handoff'
         end if
         call run_dresp09t_moving_basis(trim(config%output_file), response_space, radial_bases, ground_states, &
            reciprocal_obj, lattice_obj, hamiltonian_obj)
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_representation_tangent) then
         if (.not. use_accepted_kspace_scf) then
            error stop 'DRESP-09U representation_tangent requires the accepted k-space SCF handoff'
         end if
         call run_dresp09u_representation_tangent(trim(config%output_file), reciprocal_obj, lattice_obj, hamiltonian_obj)
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_density_moment_tangent) then
         if (.not. use_accepted_kspace_scf) then
            error stop 'DRESP-09V density_moment_tangent requires the accepted k-space SCF handoff'
         end if
         call run_dresp09v_density_moment_tangent(trim(config%output_file), response_space, radial_bases, ground_states, &
            reciprocal_obj, lattice_obj, hamiltonian_obj)
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_native_second_order) then
         if (.not. use_accepted_kspace_scf) then
            error stop 'DRESP-08 native_second_order requires the accepted k-space SCF handoff'
         end if
         call run_dresp08_native_second_order(trim(config%output_file), response_space, radial_bases, ground_states, &
            reciprocal_obj, lattice_obj, hamiltonian_obj, config%eta_values)
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_pauli_projected_native) then
         if (.not. use_accepted_kspace_scf) then
            error stop 'DRESP-09R pauli_projected_native requires the accepted k-space SCF handoff'
         end if
         call run_dresp09r_pauli_projected_native(trim(config%output_file), response_space, radial_bases, ground_states, &
            reciprocal_obj, lattice_obj, hamiltonian_obj)
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_sr_l0_scalar_relativistic) then
         if (.not. use_accepted_kspace_scf) then
            error stop 'DRESP-09S sr_l0_scalar_relativistic requires the accepted k-space SCF handoff'
         end if
         call run_dresp09s_scalar_relativistic(trim(config%output_file), response_space, radial_bases, ground_states, &
            reciprocal_obj, lattice_obj, hamiltonian_obj)
         return
      end if
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
      if (trim(config%backend) == tddft_driver_backend_alsda_compare) then
         if (.not. use_accepted_kspace_scf) then
            error stop 'DRESP-06A alsda_compare requires the accepted k-space SCF handoff'
         end if
         call run_tddft_alsda_compare(config, response_space, radial_bases, ground_states, left_state, endpoints, &
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
      if (trim(config%backend) == tddft_driver_backend_static_interactions) then
         if (.not. use_accepted_kspace_scf) then
            error stop 'TDVK-06 static interactions: direct accepted k-space SCF handoff is required'
         end if
         call run_tddft_static_interactions(config, response_space, radial_bases, ground_states, left_state, endpoints, &
            reciprocal_obj, lattice_obj)
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_product_finite_q) then
         ! TDVK-04 is a compact bare-response validation seam.  It evaluates
         ! all requested q points, the prescribed covariance pair, and only
         ! representative finite-q GF spots.  It never enters KXC, Goldstone,
         ! Dyson, loss, or the dense point-space result container.
         call run_tddft_product_finite_q(config, response_space, radial_bases, ground_states, left_state, endpoints)
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_product_convergence) then
         ! TDVK-05 is a compact numerical-convergence seam.  It evaluates the
         ! Lehmann service over the prescribed q/omega/physical-eta grid and
         ! stops before GF, KXC, Goldstone, Dyson, loss, or mode fitting.
         call run_tddft_product_convergence(config, response_space, radial_bases, ground_states, left_state, endpoints, &
            reciprocal_obj, use_accepted_kspace_scf)
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_product_gf) then
         ! Product-GF remains a bare-response validation seam.  The optional
         ! closure audit reuses this already accepted reciprocal state for the
         ! Lehmann reference and every GF control sample.
         if (config%gf_closure_audit) then
            call run_tddft_product_gf_closure(config, response_space, radial_bases, left_state, endpoints)
         else
            call run_tddft_product_gf_smoke(config, response_space, radial_bases, left_state, endpoints)
         end if
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
         call evaluate_tddft_production_sweep(config, response_space, radial_bases, ground_states, left_state, endpoints, result, &
            native_provider, native_provider%pairs, native_site_positions, accepted_pauli_magnetization, &
            lr_kxc_magnetization_source_pauli_accepted)
      else
         if (trim(config%interaction_route) == tddft_driver_route_direct_alsda) then
            call compute_accepted_pauli_magnetization(reciprocal_obj, lattice_obj%symbolic_atoms, lattice_obj%nbulk, &
               accepted_pauli_magnetization)
         end if
         call evaluate_tddft_production_sweep(config, response_space, radial_bases, ground_states, left_state, endpoints, result, &
            accepted_pauli_magnetization=accepted_pauli_magnetization, &
            accepted_pauli_magnetization_source=lr_kxc_magnetization_source_pauli_accepted)
      end if
      call write_tddft_production_output(config, result, control_obj, lattice_obj, ground_states, response_space, reciprocal_obj)
      ! The production result has already been serialized. Release the dense
      ! response matrices before the accepted-state owner leaves scope; this
      ! is important for compact smoke runs with a full radial mesh.
      if (allocated(result%ks_susceptibility)) deallocate(result%ks_susceptibility)
      if (allocated(result%enhanced_susceptibility)) deallocate(result%enhanced_susceptibility)
      if (allocated(result%loss_matrix)) deallocate(result%loss_matrix)
      if (allocated(result%reciprocal_crosscheck_delta)) deallocate(result%reciprocal_crosscheck_delta)
      if (allocated(result%reciprocal_crosscheck_valid)) deallocate(result%reciprocal_crosscheck_valid)
      if (allocated(result%reciprocal_crosscheck_norm_lehmann)) deallocate(result%reciprocal_crosscheck_norm_lehmann)
      if (allocated(result%reciprocal_crosscheck_norm_gf)) deallocate(result%reciprocal_crosscheck_norm_gf)
      if (allocated(result%reciprocal_crosscheck_difference_frobenius)) deallocate(result%reciprocal_crosscheck_difference_frobenius)
      if (allocated(result%reciprocal_crosscheck_relative_frobenius)) deallocate(result%reciprocal_crosscheck_relative_frobenius)
      if (allocated(result%reciprocal_crosscheck_difference_infinity)) deallocate(result%reciprocal_crosscheck_difference_infinity)
      if (allocated(result%q_list)) deallocate(result%q_list)
      if (allocated(result%frequencies)) deallocate(result%frequencies)
   end subroutine run_tddft_production

   !> DRESP-04 material seam.  The accepted projected site chi0 is consumed
   !> directly by the Mills interaction and site-space Dyson wrapper.
   subroutine run_tddft_projected_mills(config, response_space, radial_bases, ground_states, left_state, endpoints, &
                                        reciprocal_obj, lattice_obj)
      type(tddft_production_config), intent(in) :: config
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
      type(tddft_production_config), intent(in) :: config
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
   subroutine run_tddft_alsda_compare(config, response_space, radial_bases, ground_states, left_state, endpoints, &
                                      reciprocal_obj, lattice_obj)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(radial_ground_state), target, intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj

      type(projected_site_spin_contract), target :: contract
      type(lmto_product_response_basis), target :: product_plus, product_minus
      type(projected_chi0_request) :: site_request, opposite_site_request
      type(projected_chi0_result) :: site_bare, opposite_site_bare
      type(projected_mills_interaction_result) :: mills_result
      type(projected_juelich_request) :: juelich_request
      type(projected_juelich_result), allocatable :: static_ladder(:)
      type(projected_juelich_result) :: juelich_result
      type(lr_product_ks_susceptibility_request) :: product_request, opposite_product_request
      type(lr_product_ks_susceptibility_result) :: product_bare, opposite_product_bare
      type(lr_alsda_kernel_request) :: kxc_request
      type(lr_alsda_kernel_result) :: kxc_result
      type(tddft_dyson_request) :: alsda_request, opposite_alsda_request
      type(tddft_dyson_result) :: alsda_result, opposite_alsda_result
      type(projected_dyson_request) :: site_dyson_request
      type(projected_dyson_result) :: mills_result_dynamic, juelich_result_dynamic
      type(lr_product_gf_susceptibility_request) :: gf_request
      type(lr_product_gf_susceptibility_result) :: gf_result
      type(lr_product_ks_susceptibility_request) :: gf_lehmann_request
      type(lr_product_ks_susceptibility_result) :: gf_lehmann_result
      type(lr_alsda_kernel_result) :: old_sr_kxc_result
      real(rp), allocatable :: moment(:), magnetization(:, :), sr_magnetization(:, :), bxc_pointwise(:, :)
      complex(rp), allocatable :: test_vector(:)
      complex(rp), allocatable :: selected_static_chi(:, :, :), interaction(:, :), opposite_interaction(:, :), old_interaction(:, :)
      complex(rp), allocatable :: mcompact(:), bcompact(:), action_a(:), action_b(:), point_m(:), reconstructed_m(:)
      complex(rp), allocatable :: covariance_site_plus(:, :, :), covariance_site_minus(:, :, :)
      complex(rp), allocatable :: covariance_loss_plus(:, :, :), covariance_loss_minus(:, :, :)
      complex(rp), allocatable :: covariance_site_transport(:, :, :), covariance_loss_transport(:, :, :)
      complex(rp), allocatable :: site_alsda(:, :, :), site_loss_alsda(:, :, :), site_loss_bare(:, :), site_loss_tmp(:, :)
      real(rp) :: magnetization_norm, magnetization_projection_residual, magnetization_projection_relative
      real(rp) :: action_residual, ward_residual, ward_relative, kxc_identity_abs, kxc_identity_rel
      real(rp) :: magnetization_sr_pauli_relative, sr_moment, pauli_moment, radial_difference_max
      real(rp) :: old_min_kxc, old_max_kxc, old_maxabs_kxc, kxc_pointwise_relative, compact_kxc_relative
      real(rp) :: field_closure, direct_field_relative, kxc_ward_relative, kxc_field_seam_relative
      real(rp) :: eta_run, site_closure, loss_closure, r_alsda_mills, r_alsda_juelich
      real(rp) :: loss_trace_bare, loss_trace_mills, loss_trace_juelich, loss_trace_alsda
      real(rp) :: min_kxc, max_kxc, maxabs_kxc, denominator_min
      complex(rp) :: keff
      real(rp) :: covariance_site_error, covariance_loss_error
      real(rp) :: covariance_direct_site_error
      real(rp) :: covariance_compact_error
      real(rp) :: covariance_one
      real(rp) :: gf_norm_lehmann, gf_norm_gf, gf_difference, gf_relative, gf_infinity
      integer :: gf_q_index, gf_frequency_index
      integer :: gamma_index, iq, ieta, iw, unit, selected_eta_index, holdout_eta_index, i, j
      integer :: positive_index, negative_index
      logical :: rank_stable, eta_stable, covariance_found

      if (lattice_obj%nrec /= size(ground_states)) error stop 'DRESP-06A: lattice/site provenance mismatch'
      if (response_space%response_lmax /= 4) error stop 'DRESP-06A: complete response_lmax=4 is required'
      gamma_index = find_gamma_q_index(config%q_list)
      if (gamma_index == 0) error stop 'DRESP-06A: Gamma is required for the Juelich static construction'
      call contract%initialize(response_space, radial_bases, 'spd')
      call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      call product_minus%initialize(response_space, radial_bases, lmto_product_channel_minus, .true.)
      rank_stable = .true.
      do i = 1, product_plus%nsite
         do j = 0, product_plus%response_lmax
            rank_stable = rank_stable .and. product_plus%blocks(i, j)%rank_stable .and. &
               product_minus%blocks(i, j)%rank_stable
         end do
      end do
      if (.not. rank_stable) error stop 'DRESP-06A: complete product basis is rank-unstable'
      if (product_plus%product_dimension /= product_minus%product_dimension) then
         error stop 'DRESP-06A: plus/minus product dimensions differ'
      end if
      if (size(ground_states) == 1 .and. product_plus%product_dimension /= 232) then
         error stop 'DRESP-06A: Fe spd product dimension is not 232'
      end if

      allocate(moment(size(ground_states)), magnetization(size(ground_states), response_space%npoint), &
         sr_magnetization(size(ground_states), response_space%npoint), bxc_pointwise(size(ground_states), response_space%npoint))
      do i = 1, size(ground_states)
         sr_magnetization(i, :) = ground_states(i)%n_up - ground_states(i)%n_down
      end do
      call compute_accepted_pauli_magnetization(reciprocal_obj, lattice_obj%symbolic_atoms, lattice_obj%nbulk, magnetization)
      call contract%moment_from_operator(left_state%eigenvalues, left_state%eigenvectors, left_state%k_weights, &
         left_state%fermi_level, left_state%temperature, ground_states, moment)

      call prepare_direct_alsda_request(response_space, ground_states, magnetization, kxc_request)
      call evaluate_lr_alsda_kernel(kxc_request, kxc_result)
      call evaluate_lr_alsda_kernel(response_space, ground_states, sr_magnetization, old_sr_kxc_result)
      do i = 1, size(ground_states)
         bxc_pointwise(i, :) = 0.5_rp*(ground_states(i)%vxc_up - ground_states(i)%vxc_down)
      end do
      allocate(interaction(product_plus%product_dimension, product_plus%product_dimension), &
         opposite_interaction(product_minus%product_dimension, product_minus%product_dimension), &
         old_interaction(product_plus%product_dimension, product_plus%product_dimension))
      call compact_project_local_operator(response_space, product_plus, cmplx(kxc_result%pointwise_kernel, 0.0_rp, rp), interaction)
      call compact_project_local_operator(response_space, product_minus, cmplx(kxc_result%pointwise_kernel, 0.0_rp, rp), opposite_interaction)
      call compact_project_local_operator(response_space, product_plus, cmplx(old_sr_kxc_result%pointwise_kernel, 0.0_rp, rp), old_interaction)

      ! Independent Kxc action tests: deterministic, deterministic pseudo-random,
      ! and the actual compact ground-state magnetization vector.
      allocate(test_vector(product_plus%product_dimension), action_a(product_plus%product_dimension), &
         action_b(product_plus%product_dimension), mcompact(product_plus%product_dimension), bcompact(product_plus%product_dimension), &
         point_m(response_space%ndim), reconstructed_m(response_space%ndim))
      do i = 1, size(test_vector)
         test_vector(i) = real(mod(37*i + 11, 101), rp)/101.0_rp
      end do
      call compact_apply_local_operator(response_space, product_plus, cmplx(kxc_result%pointwise_kernel, 0.0_rp, rp), &
         test_vector, action_a)
      call compact_apply_local_operator_independent(response_space, product_plus, kxc_result%pointwise_kernel, test_vector, action_b)
      action_residual = maxval(abs(action_a - action_b))
      do i = 1, size(test_vector)
         test_vector(i) = sin(real(13*i, rp))*cos(real(7*i + 3, rp))
      end do
      call compact_apply_local_operator(response_space, product_plus, cmplx(kxc_result%pointwise_kernel, 0.0_rp, rp), &
         test_vector, action_a)
      call compact_apply_local_operator_independent(response_space, product_plus, kxc_result%pointwise_kernel, test_vector, action_b)
      action_residual = max(action_residual, maxval(abs(action_a - action_b)))
      call compact_project_magnetization(response_space, product_plus, magnetization, mcompact, point_m)
      call compact_project_magnetization(response_space, product_plus, bxc_pointwise, bcompact)
      call compact_reconstruct_point_vector(response_space, product_plus, mcompact, reconstructed_m)
      call compact_weighted_projection_diagnostics(response_space, point_m, reconstructed_m, magnetization_norm, &
         magnetization_projection_residual, magnetization_projection_relative)
      call compact_apply_local_operator(response_space, product_plus, cmplx(kxc_result%pointwise_kernel, 0.0_rp, rp), &
         mcompact, action_a)
      call compact_apply_local_operator_independent(response_space, product_plus, kxc_result%pointwise_kernel, mcompact, action_b)
      action_residual = max(action_residual, maxval(abs(action_a - action_b)))

      min_kxc = minval(kxc_result%pointwise_kernel(:, 2:))
      max_kxc = maxval(kxc_result%pointwise_kernel(:, 2:))
      maxabs_kxc = maxval(abs(kxc_result%pointwise_kernel(:, 2:)))
      kxc_identity_abs = 0.0_rp
      do i = 2, response_space%npoint
         do j = 1, size(ground_states)
            kxc_identity_abs = max(kxc_identity_abs, abs(kxc_result%pointwise_kernel(j, i)*magnetization(j, i) - &
               0.5_rp*(ground_states(j)%vxc_up(i) - ground_states(j)%vxc_down(i))))
         end do
      end do
      kxc_identity_rel = kxc_identity_abs/max(maxabs_kxc*maxval(abs(magnetization(:, 2:))), tiny(1.0_rp))
      old_min_kxc = minval(old_sr_kxc_result%pointwise_kernel(:, 2:))
      old_max_kxc = maxval(old_sr_kxc_result%pointwise_kernel(:, 2:))
      old_maxabs_kxc = maxval(abs(old_sr_kxc_result%pointwise_kernel(:, 2:)))
      kxc_pointwise_relative = maxval(abs(kxc_result%pointwise_kernel(:, 2:) - old_sr_kxc_result%pointwise_kernel(:, 2:))) / &
         max(old_maxabs_kxc, maxabs_kxc, tiny(1.0_rp))
      compact_kxc_relative = sqrt(sum(abs(interaction - old_interaction)**2))/ &
         max(sqrt(sum(abs(old_interaction)**2)), tiny(1.0_rp))
      sr_moment = 0.0_rp
      pauli_moment = 0.0_rp
      magnetization_sr_pauli_relative = 0.0_rp
      radial_difference_max = 0.0_rp
      do i = 1, size(ground_states)
         sr_moment = sr_moment + radial_volume_integral(ground_states(i), sr_magnetization(i, :))
         pauli_moment = pauli_moment + radial_volume_integral(ground_states(i), magnetization(i, :))
         magnetization_sr_pauli_relative = magnetization_sr_pauli_relative + &
            radial_volume_norm(ground_states(i), magnetization(i, :) - sr_magnetization(i, :))**2
         radial_difference_max = max(radial_difference_max, maxval(abs(magnetization(i, :) - sr_magnetization(i, :))))
      end do
      magnetization_sr_pauli_relative = sqrt(magnetization_sr_pauli_relative)/ &
         max(sqrt(sum([(radial_volume_norm(ground_states(i), magnetization(i, :))**2, i=1,size(ground_states))])), tiny(1.0_rp))

      call select_projected_juelich_eta_indices(config%eta_values, selected_eta_index, holdout_eta_index)
      allocate(static_ladder(size(config%eta_values)), selected_static_chi(size(moment), size(moment), size(config%eta_values)))
      do ieta = 1, size(config%eta_values)
         site_request%q = config%q_list(:, gamma_index)
         site_request%frequencies = [0.0_rp]
         site_request%eta = config%eta_values(ieta)
         site_request%channel = lr_channel_plus
         site_request%contract => contract
         site_request%product_basis => product_plus
         site_request%electronic_state => left_state
         site_request%q_endpoint_state => endpoints(gamma_index)
         call evaluate_projected_lehmann_chi0(site_request, site_bare)
         selected_static_chi(:, :, ieta) = site_bare%susceptibility(:, :, 1)
         juelich_request%selector = 'spd'
         juelich_request%projected_moment = moment
         juelich_request%static_chi0 = site_bare%susceptibility(:, :, 1)
         juelich_request%static_eta = config%eta_values(ieta)
         juelich_request%q = config%q_list(:, gamma_index)
         juelich_request%channel = lr_channel_plus
         juelich_request%state_provenance = 'DRESP-06A accepted state; DRESP-01 site moment'
         juelich_request%chi0_provenance = site_bare%provenance
         call evaluate_projected_juelich_interaction(juelich_request, static_ladder(ieta))
      end do
      juelich_result = static_ladder(selected_eta_index)
      if (holdout_eta_index > 0) then
         call assess_projected_juelich_eta_stability(juelich_result, static_ladder(holdout_eta_index), eta_stable)
         call evaluate_projected_juelich_holdout(selected_static_chi(:, :, holdout_eta_index), moment, &
            juelich_result%interaction_U_real, juelich_result%holdout_residual, juelich_result%holdout_relative_residual)
      else
         eta_stable = .false.
      end if
      call evaluate_projected_mills_from_reciprocal(contract, radial_bases, reciprocal_obj, left_state%fermi_level, moment, mills_result)
      if (trim(mills_result%classification) == 'UNSUPPORTED') error stop 'DRESP-06A: Mills scalar comparison unsupported'

      ! Static raw ALSDA Ward diagnostic in the same compact product space.
      product_request%q = config%q_list(:, gamma_index)
      product_request%frequencies = [0.0_rp]
      product_request%eta = config%eta_values(selected_eta_index)
      product_request%channel = lr_channel_plus
      product_request%product_basis => product_plus
      product_request%electronic_state => left_state
      product_request%q_endpoint_state => endpoints(gamma_index)
      call evaluate_lr_product_ks_susceptibility(product_request, product_bare)
      ward_residual = sqrt(sum(abs(matmul(product_bare%susceptibility(:, :, 1), matmul(interaction, mcompact)) - mcompact)**2))
      ward_relative = ward_residual/max(sqrt(sum(abs(mcompact)**2)), tiny(1.0_rp))
      denominator_min = huge(1.0_rp)
      site_closure = 0.0_rp
      loss_closure = 0.0_rp
      r_alsda_mills = 0.0_rp
      r_alsda_juelich = 0.0_rp

      open(newunit=unit, file=trim(config%output_file), status='replace', action='write')
      write(unit, '(a)') '# DRESP-06A spatial ALSDA same-state comparison'
      write(unit, '(a)') '# accepted-state provenance = one accepted reciprocal k-space SCF state; no interaction-specific re-SCF'
      write(unit, '(a)') '# backend = alsda_compare'
      write(unit, '(a,a)') '# magnetization_kind = ', trim(kxc_result%magnetization_kind)
      write(unit, '(a,a)') '# magnetization_source = ', trim(kxc_result%magnetization_source)
      write(unit, '(a,3(i0,1x))') '# k_mesh = ', reciprocal_obj%nk_mesh
      write(unit, '(a,i0)') '# k_count = ', left_state%nk
      write(unit, '(a,es24.16)') '# temperature_K = ', left_state%temperature
      write(unit, '(a,es24.16)') '# fermi_level_Ry = ', left_state%fermi_level
      write(unit, '(a,a)') '# selector = spd'
      write(unit, '(a,i0)') '# product_dimension = ', product_plus%product_dimension
      write(unit, '(a,i0)') '# response_lmax = ', response_space%response_lmax
      write(unit, '(a,a)') '# XC functional/backend/TXC = '//trim(kxc_result%xc_provenance%functional_name)//'/'// &
         trim(kxc_result%xc_provenance%backend_name)//'/'//trim(kxc_result%xc_provenance%txch)
      write(unit, '(a)') '# ALSDA kernel contract = Kxc=B_xc_sigma,SR/magnetization_Pauli; mixed P<-SR; not exact relativistic derivative'
      write(unit, '(a,3(es24.16,1x))') '# Kxc_min_max_maxabs = ', min_kxc, max_kxc, maxabs_kxc
      write(unit, '(a,2(es24.16,1x))') '# Kxc_m_identity_abs_rel = ', kxc_identity_abs, kxc_identity_rel
      write(unit, '(a,3(es24.16,1x))') '# SR_vs_Pauli_integrated_moments = ', sr_moment, pauli_moment, sr_moment-pauli_moment
      write(unit, '(a,2(es24.16,1x))') '# SR_vs_Pauli_weighted_relative_difference_max_radial = ', &
         magnetization_sr_pauli_relative, radial_difference_max
      write(unit, '(a,3(es24.16,1x))') '# OldSR_Kxc_min_max_maxabs = ', old_min_kxc, old_max_kxc, old_maxabs_kxc
      write(unit, '(a,2(es24.16,1x))') '# OldSR_vs_CorrectedPauli_Kxc_pointwise_compact_relative = ', &
         kxc_pointwise_relative, compact_kxc_relative
      write(unit, '(a)') '# SUPERSEDED material kernel = LR-01 scalar-relativistic n_up-n_down denominator; retained only as audit diagnostic'
      write(unit, '(a)') '# SR_PAULI_RADIAL columns: site radial_index m_SR m_Pauli difference'
      do i = 1, size(ground_states)
         do j = 2, response_space%npoint
            write(unit, '(a,1x,2(i0,1x),3(es24.16,1x))') 'SR_PAULI_RADIAL', i, j, sr_magnetization(i, j), &
               magnetization(i, j), magnetization(i, j)-sr_magnetization(i, j)
         end do
      end do
      write(unit, '(a)') '# SR_PAULI_BLOCK columns: site integrated_m_SR integrated_m_Pauli difference'
      do i = 1, size(ground_states)
         write(unit, '(a,1x,i0,1x,3(es24.16,1x))') 'SR_PAULI_BLOCK', i, &
            radial_volume_integral(ground_states(i), sr_magnetization(i, :)), &
            radial_volume_integral(ground_states(i), magnetization(i, :)), &
            radial_volume_integral(ground_states(i), magnetization(i, :) - sr_magnetization(i, :))
      end do
      write(unit, '(a,3(es24.16,1x))') '# magnetization_projection_norm_residual_relative = ', magnetization_norm, &
         magnetization_projection_residual, magnetization_projection_relative
      write(unit, '(a,es24.16)') '# Kxc_operator_action_oracle_max = ', action_residual
      write(unit, '(a,3(es24.16,1x))') '# U_Mills = ', mills_result%interaction_U
      write(unit, '(a,es24.16)') '# Mills_scalarization_residual = ', mills_result%scalarization_residual
      write(unit, '(a,a)') '# Mills_classification = ', trim(mills_result%classification)
      write(unit, '(a,*(es24.16,1x))') '# U_Juelich_real = ', juelich_result%interaction_U_real
      write(unit, '(a,a)') '# Juelich_classification = ', trim(juelich_result%classification)
      write(unit, '(a,es24.16)') '# Juelich_static_eta_selected = ', config%eta_values(selected_eta_index)
      write(unit, '(a,l1)') '# Juelich_eta_stable = ', eta_stable
      write(unit, '(a,es24.16)') '# Juelich_construction_relative_residual = ', juelich_result%relative_real_constrained_residual
      write(unit, '(a,es24.16)') '# Juelich_holdout_relative_residual = ', juelich_result%holdout_relative_residual
      write(unit, '(a,2(es24.16,1x))') '# ALSDA_raw_Gamma_Ward_abs_rel = ', ward_residual, ward_relative
      write(unit, '(a)') '# ALSDA correction = OFF; no BES/GCR/eigenvalue shifting/kernel rescaling'
      write(unit, '(a)') '# goldstone_correction = false'
      write(unit, '(a)') '# BES = false'
      write(unit, '(a)') '# GCR = false'
      write(unit, '(a)') '# kernel_rescaling = false'
      write(unit, '(a)') '# eigenvalue_pinning = false'
      write(unit, '(a)') '# columns: route eta_Ry q_index omega_Ry row col chi_Re chi_Im loss_Re loss_Im -ImTr/pi min_sv max_sv condition min_abs_eig Dyson_residual'
      write(unit, '(a)') '# FREQUENCY_REFINEMENT columns: route q_index coarse_peak_Ry refined_peak_Ry min_sv refined_loss_trace'
      write(unit, '(a)') '# GF_SPOT columns: quadrature q_index frequency_index N_E energy_min energy_max h h_over_integration_eta norm_Lehmann norm_GF dF relative_dF dInf'
      write(unit, '(a)') '# GF_BASE_FINE columns: q_index frequency_index dF relative_dF dInf'
      write(unit, '(a)') '# Q_CHANNEL_COVARIANCE_STAGE columns: stage positive_q_index negative_q_index frequency_index residual min_sv_difference condition_difference'
      write(unit, '(a)') '# KXC_TRANSPORT columns: positive_q_index negative_q_index relative_frobenius'
      write(unit, '(a)') '# WARD_ETA columns: eta rK_abs rK_rel rB_abs rB_rel seam_abs seam_rel field_closure R_raw min_sv max_sv condition min_abs_eig R_G0_over_R_raw geometric_alignment subtraction_ratio'
      write(unit, '(a)') '# WARD_MODE_TRACK columns: candidate lambda_real lambda_imag lambda_abs right_overlap biorthogonal_weight m_overlap svd_overlap R_raw R_G0'
      write(unit, '(a)') '# WARD_MODE columns: eta mode lambda_real lambda_imag lambda_abs right_overlap biorthogonal_weight m_overlap nonadditive_modal_norm_diagnostic'
      write(unit, '(a)') '# WARD_SVD columns: eta min_sv max_sv condition smallest_right_overlap'
      write(unit, '(a)') '# WARD_DM_RECONSTRUCTION columns: eta eigen_reconstruction_relative dm_reconstruction_relative ward_reconstruction_relative'

      call run_alsda_ward_eta_ladder(unit, config, product_plus, left_state, endpoints(gamma_index), interaction, mcompact, bcompact)

      allocate(site_alsda(size(moment), size(moment), size(config%frequencies)), &
         site_loss_alsda(size(moment), size(moment), size(config%frequencies)), site_loss_bare(size(moment), size(moment)), &
         site_loss_tmp(size(moment), size(moment)))
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

      if (positive_index > 0) then
         call run_projected_frequency_refinement(unit, 'spd', config, positive_index, contract, product_plus, &
            left_state, endpoints(positive_index), mills_result%interaction_U, juelich_result%interaction_U_real)
         call run_alsda_frequency_refinement(unit, config, positive_index, product_plus, left_state, endpoints(positive_index), &
            interaction)
      end if
      if (config%gf_closure_audit) then
         call run_alsda_product_gf_spot(unit, config, product_plus, left_state, endpoints)
      end if

      do ieta = 1, size(config%eta_values)
         eta_run = config%eta_values(ieta)
         do iq = 1, size(config%q_list, 2)
            site_request%q = config%q_list(:, iq)
            site_request%frequencies = config%frequencies
            site_request%eta = eta_run
            site_request%channel = config%channel
            if (trim(config%channel) == lr_channel_minus) then
               site_request%product_basis => product_minus
            else
               site_request%product_basis => product_plus
            end if
            site_request%q_endpoint_state => endpoints(iq)
            call evaluate_projected_lehmann_chi0(site_request, site_bare)
            product_request%q = config%q_list(:, iq)
            product_request%frequencies = config%frequencies
            product_request%eta = eta_run
            product_request%channel = config%channel
            if (trim(config%channel) == lr_channel_minus) then
               product_request%product_basis => product_minus
            else
               product_request%product_basis => product_plus
            end if
            product_request%electronic_state => left_state
            product_request%q_endpoint_state => endpoints(iq)
            call evaluate_lr_product_ks_susceptibility(product_request, product_bare)
            call project_compact_response_to_sites(contract, product_request%product_basis, product_bare%susceptibility, site_alsda)
            site_closure = max(site_closure, maxval(abs(site_bare%susceptibility - site_alsda)))
            alsda_request%compact_orthonormal = .true.
            alsda_request%q = config%q_list(:, iq)
            alsda_request%frequencies = config%frequencies
            alsda_request%eta = eta_run
            alsda_request%channel = config%channel
            alsda_request%ks_susceptibility = product_bare%susceptibility
            if (trim(config%channel) == lr_channel_minus) then
               alsda_request%canonical_interaction = opposite_interaction
            else
               alsda_request%canonical_interaction = interaction
            end if
            alsda_request%interaction_route = lr_dyson_route_direct_alsda
            alsda_request%interaction_provenance = 'KXC-01 direct ALSDA LR-03; DRESP-06A full product space'
            alsda_request%electronic_state_provenance = 'accepted reciprocal state; SCF occupations and EF fixed'
            alsda_request%response_space_metadata = 'complete weighted-orthonormal product response; no point Dyson'
            call evaluate_tddft_dyson(alsda_request, alsda_result)
            if (any(.not. alsda_result%solve_succeeded)) error stop 'DRESP-06A: ALSDA Dyson solve failed'
            call project_compact_response_to_sites(contract, product_request%product_basis, alsda_result%enhanced_susceptibility, site_alsda)
            call project_compact_loss_to_sites(contract, product_request%product_basis, alsda_result%enhanced_susceptibility, site_loss_alsda)

            allocate(site_dyson_request%frequencies(size(config%frequencies)), site_dyson_request%interaction_U(size(moment)), &
               site_dyson_request%bare_chi(size(moment), size(moment), size(config%frequencies)))
            site_dyson_request%selector = 'spd'
            site_dyson_request%q = config%q_list(:, iq)
            site_dyson_request%frequencies = config%frequencies
            site_dyson_request%eta = eta_run
            site_dyson_request%channel = config%channel
            site_dyson_request%bare_chi = site_bare%susceptibility
            site_dyson_request%interaction_U = mills_result%interaction_U
            site_dyson_request%interaction_provenance = 'same-state DRESP-04 Mills interaction'
            call evaluate_projected_dyson(site_dyson_request, mills_result_dynamic)
            site_dyson_request%interaction_U = juelich_result%interaction_U_real
            site_dyson_request%interaction_provenance = juelich_result%provenance
            call evaluate_projected_dyson(site_dyson_request, juelich_result_dynamic)

            do iw = 1, size(config%frequencies)
               call site_loss_matrix(site_bare%susceptibility(:, :, iw), site_loss_bare)
               call site_loss_matrix(site_alsda(:, :, iw), site_loss_tmp)
               loss_closure = max(loss_closure, maxval(abs(site_loss_alsda(:, :, iw) - site_loss_tmp)))
               call site_loss_matrix(mills_result_dynamic%enhanced_chi(:, :, iw), site_loss_tmp)
               loss_trace_mills = real(sum([(mills_result_dynamic%loss_matrix(i, i, iw), i=1,size(moment))]), rp)
               loss_trace_juelich = real(sum([(juelich_result_dynamic%loss_matrix(i, i, iw), i=1,size(moment))]), rp)
               loss_trace_alsda = real(sum([(site_loss_alsda(i, i, iw), i=1,size(moment))]), rp)
               loss_trace_bare = real(sum([(site_loss_bare(i, i), i=1,size(moment))]), rp)
               r_alsda_mills = max(r_alsda_mills, relative_site_difference(site_alsda(:, :, iw), mills_result_dynamic%enhanced_chi(:, :, iw)))
               r_alsda_juelich = max(r_alsda_juelich, relative_site_difference(site_alsda(:, :, iw), juelich_result_dynamic%enhanced_chi(:, :, iw)))
               if (size(moment) == 1 .and. abs(site_bare%susceptibility(1, 1, iw)) > sqrt(tiny(1.0_rp)) .and. &
                   abs(site_alsda(1, 1, iw)) > sqrt(tiny(1.0_rp))) then
                  keff = 1.0_rp/site_bare%susceptibility(1, 1, iw) - 1.0_rp/site_alsda(1, 1, iw)
                  write(unit, '(a,1x,es24.16,1x,i0,1x,es24.16,1x,3(es24.16,1x))') 'KEFF_ALSDA', eta_run, iq, &
                     config%frequencies(iw), real(keff, rp), aimag(keff), abs(keff)
               end if
               do j = 1, size(moment)
                  do i = 1, size(moment)
                     write(unit, '(a,1x,a,1x,es24.16,1x,i0,1x,es24.16,1x,2(i0,1x),12(es24.16,1x))') 'ROW', 'BARE_SITE', eta_run, iq, &
                        config%frequencies(iw), i, j, real(site_bare%susceptibility(i,j,iw),rp), aimag(site_bare%susceptibility(i,j,iw)), &
                        real(site_loss_bare(i,j),rp), aimag(site_loss_bare(i,j)), loss_trace_bare
                     write(unit, '(a,1x,a,1x,es24.16,1x,i0,1x,es24.16,1x,2(i0,1x),12(es24.16,1x))') 'ROW', 'MILLS_SITE', eta_run, iq, &
                        config%frequencies(iw), i, j, real(mills_result_dynamic%enhanced_chi(i,j,iw),rp), aimag(mills_result_dynamic%enhanced_chi(i,j,iw)), &
                        real(mills_result_dynamic%loss_matrix(i,j,iw),rp), aimag(mills_result_dynamic%loss_matrix(i,j,iw)), loss_trace_mills, &
                        mills_result_dynamic%denominator_min_singular_value(iw), mills_result_dynamic%denominator_max_singular_value(iw), &
                        mills_result_dynamic%condition_number(iw), mills_result_dynamic%minimum_magnitude_eigenvalue(iw), &
                        mills_result_dynamic%dyson_residual(iw)
                     write(unit, '(a,1x,a,1x,es24.16,1x,i0,1x,es24.16,1x,2(i0,1x),12(es24.16,1x))') 'ROW', 'JUELICH_SITE', eta_run, iq, &
                        config%frequencies(iw), i, j, real(juelich_result_dynamic%enhanced_chi(i,j,iw),rp), aimag(juelich_result_dynamic%enhanced_chi(i,j,iw)), &
                        real(juelich_result_dynamic%loss_matrix(i,j,iw),rp), aimag(juelich_result_dynamic%loss_matrix(i,j,iw)), loss_trace_juelich, &
                        juelich_result_dynamic%denominator_min_singular_value(iw), juelich_result_dynamic%denominator_max_singular_value(iw), &
                        juelich_result_dynamic%condition_number(iw), juelich_result_dynamic%minimum_magnitude_eigenvalue(iw), &
                        juelich_result_dynamic%dyson_residual(iw)
                     write(unit, '(a,1x,a,1x,es24.16,1x,i0,1x,es24.16,1x,2(i0,1x),12(es24.16,1x))') 'ROW', 'ALSDA_SITE', eta_run, iq, &
                        config%frequencies(iw), i, j, real(site_alsda(i,j,iw),rp), aimag(site_alsda(i,j,iw)), &
                        real(site_loss_alsda(i,j,iw),rp), aimag(site_loss_alsda(i,j,iw)), loss_trace_alsda, &
                        alsda_result%denominator_min_singular_value(iw), alsda_result%denominator_max_singular_value(iw), &
                        alsda_result%denominator_condition_number(iw), alsda_result%denominator_min_magnitude_eigenvalue(iw), &
                        alsda_result%dyson_residual_frobenius(iw)
                  end do
               end do
            end do
            write(unit, '(a,1x,i0,1x,es24.16,1x,3(es24.16,1x))') 'BARE_PRODUCT_SITE_CLOSURE', iq, site_closure, &
               r_alsda_mills, r_alsda_juelich, maxval(alsda_result%denominator_min_singular_value)
            deallocate(site_dyson_request%frequencies, site_dyson_request%interaction_U, site_dyson_request%bare_chi)
         end do
      end do

      if (covariance_found .and. config%validate_interacting_covariance) then
         call run_alsda_covariance_audit(unit, config, contract, product_plus, product_minus, left_state, endpoints, &
            interaction, opposite_interaction, positive_index, negative_index)
      else if (config%validate_interacting_covariance) then
         error stop 'DRESP-06A: requested q/-q covariance pair is unavailable'
      end if

      if (action_residual > 5.0e-10_rp) error stop 'DRESP-06A: independent Kxc action oracle failed'
      if (site_closure > dresp06a_projection_tolerance .or. loss_closure > dresp06a_projection_tolerance) then
         error stop 'DRESP-06A: compact-to-site or loss projection closure failed'
      end if

      write(unit, '(a,2(es24.16,1x))') '# ALSDA_vs_Mills_and_Juelich_site_relative = ', r_alsda_mills, r_alsda_juelich
      write(unit, '(a,es24.16)') '# bare_compact_to_site_closure_max = ', site_closure
      write(unit, '(a,es24.16)') '# loss_projection_commutation = ', loss_closure
      if (.not. config%validate_interacting_covariance) write(unit, '(a)') '# q/channel covariance = not requested'
      close(unit)
      if (allocated(covariance_site_plus)) then
         deallocate(covariance_site_plus, covariance_site_minus, covariance_loss_plus, covariance_loss_minus, &
            covariance_site_transport, covariance_loss_transport)
      end if
      deallocate(moment, magnetization, sr_magnetization, bxc_pointwise, static_ladder, selected_static_chi, interaction, &
         opposite_interaction, old_interaction, test_vector, action_a, action_b, mcompact, bcompact, point_m, reconstructed_m, &
         site_alsda, site_loss_alsda, site_loss_bare, site_loss_tmp)
   end subroutine run_tddft_alsda_compare

   !> Refine the ALSDA loss peak on a window selected from the ALSDA coarse
   !> response itself.  This is intentionally independent from the Mills and
   !> Juelich projected-site windows.
   subroutine run_alsda_ward_eta_ladder(unit, config, product, left_state, endpoint, interaction, magnetization, direct_field)
      integer, intent(in) :: unit
      type(tddft_production_config), intent(in) :: config
      type(lmto_product_response_basis), target, intent(in) :: product
      type(lr_electronic_state), target, intent(in) :: left_state, endpoint
      complex(rp), intent(in) :: interaction(:, :)
      complex(rp), intent(in) :: magnetization(:)
      complex(rp), intent(in) :: direct_field(:)
      type(lr_product_ks_susceptibility_request) :: request
      type(lr_product_ks_susceptibility_result) :: bare
      type(lr_ward_mode_analysis_result) :: analysis
      complex(rp), allocatable :: denominator(:, :), denominator_zeroed(:, :), previous_vector(:)
      complex(rp), allocatable :: r_k(:), r_b(:), r_g0(:), ward_mode(:), k_field(:)
      real(rp) :: overlap, previous_norm, current_norm, r_k_abs, r_b_abs, seam_abs, closure
      real(rp) :: raw_relative, zeroed_relative, zeroed_ratio, geometric_alignment, subtraction_ratio
      real(rp) :: magnetization_norm, direct_field_norm, seam_relative
      complex(rp) :: lambda_g
      integer :: ieta, i, j, mode, candidate
      character(len=64) :: ward_classification
      logical :: have_previous

      if (any(shape(interaction) /= [product%product_dimension, product%product_dimension]) .or. &
          size(magnetization) /= product%product_dimension .or. size(direct_field) /= product%product_dimension) then
         error stop 'DRESP-06A Ward analysis: compact shape mismatch'
      end if
      allocate(denominator(product%product_dimension, product%product_dimension), &
         denominator_zeroed(product%product_dimension, product%product_dimension), previous_vector(product%product_dimension), &
         r_k(product%product_dimension), r_b(product%product_dimension), r_g0(product%product_dimension), &
         ward_mode(product%product_dimension), k_field(product%product_dimension))
      magnetization_norm = sqrt(sum(abs(magnetization)**2))
      direct_field_norm = sqrt(sum(abs(direct_field)**2))
      if (magnetization_norm <= tiny(1.0_rp) .or. direct_field_norm <= tiny(1.0_rp)) then
         error stop 'DRESP-06A Ward analysis: magnetization or direct field is zero'
      end if
      have_previous = .false.
      ward_classification = 'UNSET'
      do ieta = 1, size(config%eta_values)
         request%q = 0.0_rp
         request%frequencies = [0.0_rp]
         request%eta = config%eta_values(ieta)
         request%channel = lr_channel_plus
         request%product_basis => product
         request%electronic_state => left_state
         request%q_endpoint_state => endpoint
         call evaluate_lr_product_ks_susceptibility(request, bare)
         denominator = -matmul(bare%susceptibility(:, :, 1), interaction)
         do i = 1, product%product_dimension
            denominator(i, i) = denominator(i, i) + cmplx(1.0_rp, 0.0_rp, rp)
         end do
         call analyze_lr_ward_mode(denominator, magnetization, analysis)
         if (.not. analysis%succeeded) error stop 'DRESP-06A Ward analysis: non-Hermitian eigensystem failed'

         if (.not. have_previous) then
            candidate = analysis%nearest_zero_index
         else
            candidate = 1
            current_norm = sqrt(sum(abs(analysis%right_eigenvectors(:, 1))**2))
            overlap = abs(dot_product(analysis%right_eigenvectors(:, 1), previous_vector))/ &
               max(current_norm*previous_norm, tiny(1.0_rp))
            do mode = 2, product%product_dimension
               current_norm = sqrt(sum(abs(analysis%right_eigenvectors(:, mode))**2))
               if (abs(dot_product(analysis%right_eigenvectors(:, mode), previous_vector))/ &
                   max(current_norm*previous_norm, tiny(1.0_rp)) > overlap) then
                  candidate = mode
                  overlap = abs(dot_product(analysis%right_eigenvectors(:, mode), previous_vector))/ &
                     max(current_norm*previous_norm, tiny(1.0_rp))
               end if
            end do
         end if
         previous_norm = sqrt(sum(abs(analysis%right_eigenvectors(:, candidate))**2))
         previous_vector = analysis%right_eigenvectors(:, candidate)/max(previous_norm, tiny(1.0_rp))
         have_previous = .true.

         ! The current residual is evaluated in its physical sign convention;
         ! the eigensystem uses D*m = -(chi*K*m-m), so the norms agree.
         k_field = matmul(interaction, magnetization)
         r_k = matmul(bare%susceptibility(:, :, 1), k_field) - magnetization
         r_b = matmul(bare%susceptibility(:, :, 1), direct_field) - magnetization
         r_k_abs = sqrt(sum(abs(r_k)**2))
         r_b_abs = sqrt(sum(abs(r_b)**2))
         seam_abs = sqrt(sum(abs(r_k - r_b)**2))
         seam_relative = seam_abs/max(magnetization_norm, tiny(1.0_rp))
         closure = sqrt(sum(abs(k_field - direct_field)**2))/direct_field_norm

         ! Diagnostic-only direct one-mode removal.  Since the left vectors
         ! are stored as kets, v*l^dagger is v*conjg(l) in this array form.
         lambda_g = analysis%eigenvalues(candidate)
         ward_mode = lambda_g*analysis%coefficients(candidate)*analysis%right_eigenvectors(:, candidate)
         denominator_zeroed = denominator
         do i = 1, product%product_dimension
            do j = 1, product%product_dimension
               denominator_zeroed(i, j) = denominator(i, j) - lambda_g*analysis%right_eigenvectors(i, candidate)* &
                  conjg(analysis%left_eigenvectors(j, candidate))
            end do
         end do
         r_g0 = matmul(denominator_zeroed, magnetization)
         if (maxval(abs(r_g0 - (analysis%ward_vector - ward_mode))) > 1.0e-8_rp*max(1.0_rp, analysis%ward_residual)) then
            error stop 'DRESP-06A Ward analysis: one-mode-zeroed denominator identity failed'
         end if
         raw_relative = analysis%ward_residual/max(magnetization_norm, tiny(1.0_rp))
         zeroed_relative = sqrt(sum(abs(r_g0)**2))/max(magnetization_norm, tiny(1.0_rp))
         zeroed_ratio = zeroed_relative/max(raw_relative, tiny(1.0_rp))
         geometric_alignment = abs(dot_product(ward_mode, analysis%ward_vector))/ &
            max(sqrt(sum(abs(ward_mode)**2))*analysis%ward_residual, tiny(1.0_rp))
         subtraction_ratio = sqrt(sum(abs(analysis%ward_vector - ward_mode)**2))/ &
            max(analysis%ward_residual, tiny(1.0_rp))

         write(unit, '(a,1x,*(es24.16,1x))') 'WARD_ETA', config%eta_values(ieta), r_k_abs, &
            r_k_abs/max(magnetization_norm, tiny(1.0_rp)), r_b_abs, r_b_abs/max(magnetization_norm, tiny(1.0_rp)), &
            seam_abs, seam_relative, closure, &
            analysis%ward_residual/max(magnetization_norm, tiny(1.0_rp)), analysis%min_singular_value, analysis%max_singular_value, &
            analysis%condition_number, analysis%minimum_magnitude_eigenvalue, zeroed_ratio, geometric_alignment, subtraction_ratio
         write(unit, '(a,1x,i0,1x,*(es24.16,1x))') 'WARD_MODE_TRACK', candidate, real(lambda_g, rp), aimag(lambda_g), &
            abs(lambda_g), analysis%right_overlap(candidate), analysis%biorthogonal_weight(candidate), &
            analysis%magnetization_mode_fraction(candidate), analysis%singular_vector_overlap(minloc(analysis%singular_values, dim=1)), &
            raw_relative, zeroed_relative
         do mode = 1, min(6, product%product_dimension)
            i = analysis%magnitude_order(mode)
            write(unit, '(a,1x,i0,1x,*(es24.16,1x))') 'WARD_MODE', i, real(analysis%eigenvalues(i), rp), &
               aimag(analysis%eigenvalues(i)), abs(analysis%eigenvalues(i)), analysis%right_overlap(i), &
               analysis%biorthogonal_weight(i), analysis%magnetization_mode_fraction(i), analysis%ward_modal_norm_diagnostic(i)
         end do
         mode = minloc(analysis%singular_values, dim=1)
         write(unit, '(a,1x,*(es24.16,1x))') 'WARD_SVD', config%eta_values(ieta), analysis%min_singular_value, &
            analysis%max_singular_value, analysis%condition_number, analysis%singular_vector_overlap(mode)
         write(unit, '(a,1x,*(es24.16,1x))') 'WARD_DM_RECONSTRUCTION', config%eta_values(ieta), &
            analysis%eigen_reconstruction_residual, analysis%ward_reconstruction_residual, raw_relative, zeroed_relative
         if (r_k_abs/max(magnetization_norm, tiny(1.0_rp)) > 1.0e-1_rp .and. &
             r_b_abs/max(magnetization_norm, tiny(1.0_rp)) > 1.0e-1_rp) then
            ward_classification = 'RESPONSE_FIELD_INCONSISTENCY'
         else if (closure > 1.0e-1_rp) then
            ward_classification = 'PRODUCT_KERNEL_CLOSURE_INCONSISTENCY'
         else if (zeroed_ratio < 1.0e-1_rp) then
            ward_classification = 'HALLE_LIKE_SINGLE_MODE'
         else if (zeroed_ratio < 5.0e-1_rp) then
            ward_classification = 'SINGLE_MODE_PARTIAL'
         else
            ward_classification = 'MIXED'
         end if
      end do
      write(unit, '(a,1x,a)') 'WARD_FINAL_CLASSIFICATION', trim(ward_classification)
      write(unit, '(a)') 'WARD_BES_ELIGIBILITY diagnostic-only; production BES/Halle correction = NO'
      deallocate(denominator, denominator_zeroed, previous_vector, r_k, r_b, r_g0, ward_mode, k_field)
   end subroutine run_alsda_ward_eta_ladder

   subroutine run_alsda_covariance_audit(unit, config, contract, plus_product, minus_product, left_state, endpoints, &
                                         interaction, opposite_interaction, positive_index, negative_index)
      integer, intent(in) :: unit, positive_index, negative_index
      type(tddft_production_config), intent(in) :: config
      type(projected_site_spin_contract), target, intent(in) :: contract
      type(lmto_product_response_basis), target, intent(in) :: plus_product, minus_product
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      complex(rp), intent(in) :: interaction(:, :), opposite_interaction(:, :)
      type(lr_product_ks_susceptibility_request) :: plus_request, minus_request
      type(lr_product_ks_susceptibility_result) :: plus_bare, minus_bare
      type(tddft_dyson_request) :: plus_request_dyson, minus_request_dyson
      type(tddft_dyson_result) :: plus_dyson, minus_dyson
      complex(rp), allocatable :: transport(:, :), mapped(:, :), mapped_loss(:, :), site_plus(:, :), site_minus(:, :)
      complex(rp), allocatable :: site_transport(:, :), site_loss_plus(:, :), site_loss_transport(:, :)
      complex(rp), allocatable :: kxc_transport(:, :), denominator_transport(:, :), identity(:, :)
      real(rp) :: kxc_residual, residual, loss_residual, site_residual, site_loss_residual
      real(rp) :: direct_site_residual, denominator_residual, sv_difference, condition_difference
      integer :: i, iw
      logical :: pass

      allocate(transport(plus_product%product_dimension, minus_product%product_dimension), &
         mapped(plus_product%product_dimension, plus_product%product_dimension), &
         mapped_loss(plus_product%product_dimension, plus_product%product_dimension), &
         denominator_transport(plus_product%product_dimension, plus_product%product_dimension), &
         identity(plus_product%product_dimension, plus_product%product_dimension), &
         site_plus(contract%nsite, contract%nsite), site_minus(contract%nsite, contract%nsite), &
         site_transport(contract%nsite, contract%nsite), site_loss_plus(contract%nsite, contract%nsite), &
         site_loss_transport(contract%nsite, contract%nsite), kxc_transport(plus_product%product_dimension, plus_product%product_dimension))
      call build_compact_covariance_transport(plus_product, minus_product, transport)
      identity = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, plus_product%product_dimension
         identity(i, i) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
      call transport_compact_matrix(plus_product, minus_product, opposite_interaction, kxc_transport, .false.)
      kxc_residual = dresp06a_relative_matrix_difference(interaction, kxc_transport)
      write(unit, '(a,1x,2(i0,1x),es24.16)') 'KXC_TRANSPORT', positive_index, negative_index, kxc_residual
      pass = kxc_residual <= 5.0e-8_rp

      plus_request%q = config%q_list(:, positive_index)
      plus_request%frequencies = config%frequencies
      plus_request%eta = config%eta
      plus_request%channel = lr_channel_plus
      plus_request%product_basis => plus_product
      plus_request%electronic_state => left_state
      plus_request%q_endpoint_state => endpoints(positive_index)
      call evaluate_lr_product_ks_susceptibility(plus_request, plus_bare)
      minus_request%q = config%q_list(:, negative_index)
      minus_request%frequencies = -config%frequencies
      minus_request%eta = config%eta
      minus_request%channel = lr_channel_minus
      minus_request%product_basis => minus_product
      minus_request%electronic_state => left_state
      minus_request%q_endpoint_state => endpoints(negative_index)
      call evaluate_lr_product_ks_susceptibility(minus_request, minus_bare)

      plus_request_dyson%compact_orthonormal = .true.
      plus_request_dyson%q = plus_request%q
      plus_request_dyson%frequencies = plus_request%frequencies
      plus_request_dyson%eta = config%eta
      plus_request_dyson%channel = lr_channel_plus
      plus_request_dyson%ks_susceptibility = plus_bare%susceptibility
      plus_request_dyson%canonical_interaction = interaction
      plus_request_dyson%interaction_route = lr_dyson_route_direct_alsda
      plus_request_dyson%interaction_provenance = 'KXC-01 direct ALSDA LR-03; DRESP-06A covariance'
      plus_request_dyson%electronic_state_provenance = 'accepted reciprocal state; SCF occupations and EF fixed'
      plus_request_dyson%response_space_metadata = 'complete weighted-orthonormal product response; covariance audit'
      call evaluate_tddft_dyson(plus_request_dyson, plus_dyson)
      if (any(.not. plus_dyson%solve_succeeded)) error stop 'DRESP-06A covariance: plus Dyson solve failed'
      minus_request_dyson%compact_orthonormal = .true.
      minus_request_dyson%q = minus_request%q
      minus_request_dyson%frequencies = minus_request%frequencies
      minus_request_dyson%eta = config%eta
      minus_request_dyson%channel = lr_channel_minus
      minus_request_dyson%ks_susceptibility = minus_bare%susceptibility
      minus_request_dyson%canonical_interaction = opposite_interaction
      minus_request_dyson%interaction_route = lr_dyson_route_direct_alsda
      minus_request_dyson%interaction_provenance = 'KXC-01 direct ALSDA LR-03; DRESP-06A covariance'
      minus_request_dyson%electronic_state_provenance = 'accepted reciprocal state; SCF occupations and EF fixed'
      minus_request_dyson%response_space_metadata = 'complete weighted-orthonormal product response; covariance audit'
      call evaluate_tddft_dyson(minus_request_dyson, minus_dyson)
      if (any(.not. minus_dyson%solve_succeeded)) error stop 'DRESP-06A covariance: minus Dyson solve failed'

      do iw = 1, size(config%frequencies)
         call transport_compact_matrix(plus_product, minus_product, minus_bare%susceptibility(:, :, iw), mapped)
         residual = maxval(abs(plus_bare%susceptibility(:, :, iw) - mapped))
         write(unit, '(a,1x,a,1x,2(i0,1x),i0,1x,3(es24.16,1x))') 'Q_CHANNEL_COVARIANCE_STAGE', 'BARE', positive_index, &
            negative_index, iw, residual, 0.0_rp, 0.0_rp
         pass = pass .and. residual <= 5.0e-8_rp

         denominator_transport = matmul(transport, matmul(conjg(minus_dyson%denominator(:, :, iw)), &
            conjg(transpose(transport))))
         denominator_residual = maxval(abs(plus_dyson%denominator(:, :, iw) - denominator_transport))
         sv_difference = abs(plus_dyson%denominator_min_singular_value(iw) - minus_dyson%denominator_min_singular_value(iw))
         condition_difference = abs(plus_dyson%denominator_condition_number(iw) - minus_dyson%denominator_condition_number(iw))
         write(unit, '(a,1x,a,1x,2(i0,1x),i0,1x,3(es24.16,1x))') 'Q_CHANNEL_COVARIANCE_STAGE', 'DENOMINATOR', &
            positive_index, negative_index, iw, denominator_residual, sv_difference, condition_difference
         pass = pass .and. denominator_residual <= 5.0e-8_rp .and. sv_difference <= 5.0e-8_rp .and. &
            condition_difference <= 5.0e-6_rp*max(1.0_rp, plus_dyson%denominator_condition_number(iw))

         call transport_compact_matrix(plus_product, minus_product, minus_dyson%enhanced_susceptibility(:, :, iw), mapped)
         residual = maxval(abs(plus_dyson%enhanced_susceptibility(:, :, iw) - mapped))
         write(unit, '(a,1x,a,1x,2(i0,1x),i0,1x,3(es24.16,1x))') 'Q_CHANNEL_COVARIANCE_STAGE', 'INTERACTING_COMPACT', &
            positive_index, negative_index, iw, residual, 0.0_rp, 0.0_rp
         pass = pass .and. residual <= 5.0e-8_rp

         call transport_compact_matrix(plus_product, minus_product, -minus_dyson%loss_matrix(:, :, iw), mapped_loss)
         loss_residual = maxval(abs(plus_dyson%loss_matrix(:, :, iw) - mapped_loss))
         write(unit, '(a,1x,a,1x,2(i0,1x),i0,1x,3(es24.16,1x))') 'Q_CHANNEL_COVARIANCE_STAGE', 'COMPACT_LOSS', &
            positive_index, negative_index, iw, loss_residual, 0.0_rp, 0.0_rp
         pass = pass .and. loss_residual <= 5.0e-8_rp

         call project_compact_response_to_sites(contract, plus_product, plus_dyson%enhanced_susceptibility(:, :, iw), site_plus)
         call project_compact_response_to_sites(contract, minus_product, minus_dyson%enhanced_susceptibility(:, :, iw), site_minus)
         call project_compact_response_to_sites(contract, plus_product, mapped, site_transport)
         site_residual = maxval(abs(site_plus - site_transport))
         direct_site_residual = maxval(abs(site_plus - conjg(site_minus)))
         write(unit, '(a,1x,a,1x,2(i0,1x),i0,1x,3(es24.16,1x))') 'Q_CHANNEL_COVARIANCE_STAGE', 'PROJECTED_SITE', &
            positive_index, negative_index, iw, site_residual, direct_site_residual, 0.0_rp
         pass = pass .and. site_residual <= 5.0e-8_rp

         call project_compact_loss_to_sites(contract, plus_product, plus_dyson%enhanced_susceptibility(:, :, iw), site_loss_plus)
         call project_compact_response_to_sites(contract, plus_product, mapped_loss, site_loss_transport)
         site_loss_residual = maxval(abs(site_loss_plus - site_loss_transport))
         write(unit, '(a,1x,a,1x,2(i0,1x),i0,1x,3(es24.16,1x))') 'Q_CHANNEL_COVARIANCE_STAGE', 'PROJECTED_LOSS', &
            positive_index, negative_index, iw, site_loss_residual, direct_site_residual, 0.0_rp
         pass = pass .and. site_loss_residual <= 5.0e-8_rp
      end do
      if (pass) then
         write(unit, '(a)') '# q/channel covariance status = PASS'
      else
         write(unit, '(a)') '# q/channel covariance status = BLOCKED — covariance stage failed'
         error stop 'BLOCKED — DRESP-06A commensurate covariance stage failed'
      end if
      deallocate(transport, mapped, mapped_loss, denominator_transport, identity, site_plus, site_minus, site_transport, &
         site_loss_plus, site_loss_transport, kxc_transport)
   end subroutine run_alsda_covariance_audit

   subroutine run_alsda_frequency_refinement(unit, config, q_index, product, left_state, endpoint, interaction)
      integer, intent(in) :: unit, q_index
      type(tddft_production_config), intent(in) :: config
      type(lmto_product_response_basis), target, intent(in) :: product
      type(lr_electronic_state), target, intent(in) :: left_state, endpoint
      complex(rp), intent(in) :: interaction(:, :)
      type(lr_product_ks_susceptibility_request) :: request
      type(lr_product_ks_susceptibility_result) :: coarse_bare, refined_bare
      type(tddft_dyson_request) :: dyson_request
      type(tddft_dyson_result) :: coarse_result, refined_result
      real(rp), allocatable :: refined_frequencies(:), coarse_loss(:), refined_loss(:)
      real(rp) :: eta, step, lower, upper
      integer :: selected_eta_index, holdout_eta_index, coarse_peak, refined_peak, nref, iw, i

      if (size(config%frequencies) < 2) return
      call select_projected_juelich_eta_indices(config%eta_values, selected_eta_index, holdout_eta_index)
      eta = config%eta_values(selected_eta_index)
      request%q = config%q_list(:, q_index)
      request%frequencies = config%frequencies
      request%eta = eta
      request%channel = config%channel
      request%product_basis => product
      request%electronic_state => left_state
      request%q_endpoint_state => endpoint
      call evaluate_lr_product_ks_susceptibility(request, coarse_bare)
      dyson_request%compact_orthonormal = .true.
      dyson_request%q = request%q
      dyson_request%frequencies = request%frequencies
      dyson_request%eta = eta
      dyson_request%channel = request%channel
      dyson_request%ks_susceptibility = coarse_bare%susceptibility
      dyson_request%canonical_interaction = interaction
      dyson_request%interaction_route = lr_dyson_route_direct_alsda
      dyson_request%interaction_provenance = 'KXC-01 direct ALSDA; DRESP-06A route-specific refinement'
      dyson_request%electronic_state_provenance = 'accepted reciprocal state; SCF occupations and EF fixed'
      dyson_request%response_space_metadata = 'complete weighted-orthonormal product response; no point Dyson'
      call evaluate_tddft_dyson(dyson_request, coarse_result)
      allocate(coarse_loss(size(config%frequencies)))
      do iw = 1, size(config%frequencies)
         coarse_loss(iw) = 0.0_rp
         do i = 1, product%product_dimension
            coarse_loss(iw) = coarse_loss(iw) + real(coarse_result%loss_matrix(i, i, iw), rp)
         end do
      end do
      coarse_peak = maxloc(coarse_loss, dim=1)
      step = minval(abs(config%frequencies(2:) - config%frequencies(:size(config%frequencies)-1)))
      lower = max(minval(config%frequencies), config%frequencies(coarse_peak) - step)
      upper = min(maxval(config%frequencies), config%frequencies(coarse_peak) + step)
      if (upper <= lower) then
         deallocate(coarse_loss)
         return
      end if
      nref = 41
      allocate(refined_frequencies(nref), refined_loss(nref))
      do iw = 1, nref
         refined_frequencies(iw) = lower + real(iw - 1, rp)*(upper - lower)/real(nref - 1, rp)
      end do
      request%frequencies = refined_frequencies
      call evaluate_lr_product_ks_susceptibility(request, refined_bare)
      dyson_request%frequencies = refined_frequencies
      dyson_request%ks_susceptibility = refined_bare%susceptibility
      call evaluate_tddft_dyson(dyson_request, refined_result)
      do iw = 1, nref
         refined_loss(iw) = 0.0_rp
         do i = 1, product%product_dimension
            refined_loss(iw) = refined_loss(iw) + real(refined_result%loss_matrix(i, i, iw), rp)
         end do
      end do
      refined_peak = maxloc(refined_loss, dim=1)
      write(unit, '(a,1x,i0,1x,4(es24.16,1x))') 'FREQUENCY_REFINEMENT ALSDA', q_index, &
         config%frequencies(coarse_peak), refined_frequencies(refined_peak), &
         minval(refined_result%denominator_min_singular_value), refined_loss(refined_peak)
      deallocate(coarse_loss, refined_frequencies, refined_loss)
   end subroutine run_alsda_frequency_refinement

   !> One representative finite-q bare product-GF spot for the ALSDA
   !> comparison.  It is a closure audit of the frozen accepted state, not an
   !> interaction-specific GF construction.
   subroutine run_alsda_product_gf_spot(unit, config, product, left_state, endpoints)
      integer, intent(in) :: unit
      type(tddft_production_config), intent(in) :: config
      type(lmto_product_response_basis), target, intent(in) :: product
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(lr_product_gf_susceptibility_request) :: gf_request
      type(lr_product_gf_susceptibility_result) :: gf_result, gf_fine_result
      type(lr_product_ks_susceptibility_request) :: lehmann_request
      type(lr_product_ks_susceptibility_result) :: lehmann_result
      integer :: q_index, frequency_index, base_points, fine_points, base_intervals
      real(rp) :: norm_lehmann, norm_gf, difference, relative, infinity
      real(rp) :: fine_norm_gf, base_fine_difference, base_fine_relative, base_fine_infinity
      real(rp) :: energy_min, energy_max, span

      q_index = find_first_nonzero_q(config%q_list)
      if (q_index == 0) error stop 'DRESP-06A: GF audit requires a finite-q spot'
      frequency_index = 1
      do while (frequency_index <= size(config%frequencies))
         if (abs(config%frequencies(frequency_index)) > 1.0e-12_rp) exit
         frequency_index = frequency_index + 1
      end do
      if (frequency_index > size(config%frequencies)) frequency_index = 1

      if (config%gf_integration_eta <= 0.0_rp) error stop 'DRESP-06A: GF audit requires positive integration_eta'
      energy_min = min(minval(left_state%eigenvalues), minval(endpoints(q_index)%eigenvalues)) - config%gf_energy_margin
      energy_max = max(maxval(left_state%eigenvalues), maxval(endpoints(q_index)%eigenvalues)) + config%gf_energy_margin
      span = energy_max - energy_min
      base_intervals = max(2, ceiling(span/(0.4_rp*config%gf_integration_eta)))
      if (mod(base_intervals, 2) /= 0) base_intervals = base_intervals + 1
      base_points = base_intervals + 1
      fine_points = 2*base_intervals + 1

      gf_request%q = config%q_list(:, q_index)
      gf_request%frequencies = [config%frequencies(frequency_index)]
      gf_request%eta = config%eta
      gf_request%channel = config%channel
      gf_request%integration_points = base_points
      gf_request%integration_eta = config%gf_integration_eta
      gf_request%energy_margin = config%gf_energy_margin
      gf_request%contraction_backend = 'factorized'
      gf_request%product_basis => product
      gf_request%electronic_state => left_state
      gf_request%q_endpoint_state => endpoints(q_index)
      call evaluate_lr_product_gf_susceptibility(gf_request, gf_result)

      gf_request%integration_points = fine_points
      call evaluate_lr_product_gf_susceptibility(gf_request, gf_fine_result)

      lehmann_request%q = gf_request%q
      lehmann_request%frequencies = gf_request%frequencies
      lehmann_request%eta = config%eta
      lehmann_request%channel = config%channel
      lehmann_request%product_basis => product
      lehmann_request%electronic_state => left_state
      lehmann_request%q_endpoint_state => endpoints(q_index)
      call evaluate_lr_product_ks_susceptibility(lehmann_request, lehmann_result)
      norm_lehmann = sqrt(sum(abs(lehmann_result%susceptibility(:, :, 1))**2))
      norm_gf = sqrt(sum(abs(gf_result%susceptibility(:, :, 1))**2))
      difference = sqrt(sum(abs(lehmann_result%susceptibility(:, :, 1) - gf_result%susceptibility(:, :, 1))**2))
      relative = difference/max(norm_lehmann, norm_gf, tiny(1.0_rp))
      infinity = maxval(abs(lehmann_result%susceptibility(:, :, 1) - gf_result%susceptibility(:, :, 1)))
      if (.not. ieee_is_finite(relative)) error stop 'DRESP-06A: product-GF spot is not finite'
      write(unit, '(a,1x,a,1x,2(i0,1x),i0,1x,6(es24.16,1x),3(es24.16,1x))') 'GF_SPOT', 'base', q_index, &
         frequency_index, base_points, gf_result%energy_min, gf_result%energy_max, gf_result%energy_spacing, &
         gf_result%spacing_over_integration_eta, norm_lehmann, norm_gf, difference, relative, infinity

      fine_norm_gf = sqrt(sum(abs(gf_fine_result%susceptibility(:, :, 1))**2))
      base_fine_difference = sqrt(sum(abs(gf_result%susceptibility(:, :, 1) - gf_fine_result%susceptibility(:, :, 1))**2))
      base_fine_relative = base_fine_difference/max(norm_gf, fine_norm_gf, tiny(1.0_rp))
      base_fine_infinity = maxval(abs(gf_result%susceptibility(:, :, 1) - gf_fine_result%susceptibility(:, :, 1)))
      write(unit, '(a,1x,a,1x,2(i0,1x),i0,1x,6(es24.16,1x),3(es24.16,1x))') 'GF_SPOT', 'fine', q_index, &
         frequency_index, fine_points, gf_fine_result%energy_min, gf_fine_result%energy_max, gf_fine_result%energy_spacing, &
         gf_fine_result%spacing_over_integration_eta, norm_lehmann, fine_norm_gf, &
         sqrt(sum(abs(lehmann_result%susceptibility(:, :, 1) - gf_fine_result%susceptibility(:, :, 1))**2)), &
         sqrt(sum(abs(lehmann_result%susceptibility(:, :, 1) - gf_fine_result%susceptibility(:, :, 1))**2))/ &
            max(norm_lehmann, fine_norm_gf, tiny(1.0_rp)), &
         maxval(abs(lehmann_result%susceptibility(:, :, 1) - gf_fine_result%susceptibility(:, :, 1)))
      write(unit, '(a,1x,2(i0,1x),3(es24.16,1x))') 'GF_BASE_FINE', q_index, frequency_index, &
         base_fine_difference, base_fine_relative, base_fine_infinity
   end subroutine run_alsda_product_gf_spot

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
      type(tddft_production_config), intent(in) :: config
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
      type(tddft_production_config), intent(in) :: config
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
      type(projected_chi0_result) :: lehmann_result, gf_result, finite_width_result, minus_result, gamma_lehmann
      type(projected_site_spin_contract), pointer :: contract
      type(lmto_product_response_basis), pointer :: product
      real(rp), allocatable :: moment(:), accepted_moment(:)
      real(rp) :: norm_lehmann, norm_gf, norm_finite, difference, finite_difference, relative, finite_relative, accepted_total, moment_residual
      real(rp) :: integration_eta
      real(rp) :: endpoint_eigenvalue_checksum, endpoint_occupation_checksum, endpoint_unitarity_residual
      complex(rp) :: endpoint_eigenvector_checksum
      real(rp) :: state_eigenvalue_checksum, state_occupation_checksum, state_unitarity_residual
      complex(rp) :: state_eigenvector_checksum
      real(rp) :: radial_checksum, radial_l2_norm
      integer :: projection_index, iq, ifrequency, i, j, orbital, unit, gamma_index, positive_index, negative_index
      logical :: gamma_saved, covariance_saved
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

      integration_eta = config%gf_integration_eta
      if (integration_eta <= 0.0_rp) integration_eta = config%eta/40.0_rp
      gamma_index = find_gamma_q_index(config%q_list)
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
      write(unit, '(a)') '# verdict = algebraic projected-site backends; material GF closure is reported below'
      write(unit, '(a)') '# state_source = one accepted k-space SCF reciprocal state; shared by Lehmann and GF'
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
      write(unit, '(a,es24.16)') '# integration_eta_Ry = ', integration_eta
      write(unit, '(a)') '# eta_response_role = physical retarded response broadening'
      write(unit, '(a)') '# integration_eta_role = numerical one-electron real-axis regulator; converged away in Track B'
      write(unit, '(a)') '# intrinsic_linewidth = not determined by DRESP-02; neither eta is Landau damping'
      write(unit, '(a)') '# finite_width_oracle = direct DRESP-01 transition spectral Kubo integral; no eta_eff substitution'
      write(unit, '(a,i0)') '# gf_energy_points = ', config%gf_integration_points
      write(unit, '(a,es24.16)') '# gf_energy_margin_Ry = ', config%gf_energy_margin
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
      write(unit, '(a)') '# columns = projection q_index omega_Ry row col Lehmann_Re Lehmann_Im GF_Re GF_Im abs_diff rel_diff'
      write(unit, '(a)') '# finite_width_columns = projection q_index omega_Ry row col Lehmann_Re Lehmann_Im finite_width_Re finite_width_Im GF_Re GF_Im GF_minus_finite_Re GF_minus_finite_Im GF_minus_finite_abs GF_minus_finite_rel finite_width_minus_Lehmann_Re finite_width_minus_Lehmann_Im finite_width_minus_Lehmann_abs finite_width_minus_Lehmann_rel'

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
         gamma_saved = .false.

         do iq = 1, size(config%q_list, 2)
            request%q = config%q_list(:, iq)
            request%frequencies = config%frequencies
            request%eta = config%eta
            request%channel = lr_channel_plus
            request%integration_points = config%gf_integration_points
            request%integration_eta = integration_eta
            request%energy_margin = config%gf_energy_margin
            request%contract => contract
            request%product_basis => product
            request%electronic_state => left_state
            request%q_endpoint_state => endpoints(iq)
            request%channel = config%channel
            request%diagnostics = config%gf_closure_audit .and. iq == gamma_index
            call evaluate_projected_lehmann_chi0(request, lehmann_result)
            call evaluate_projected_gf_chi0(request, gf_result)
            request%diagnostics = .false.
            call evaluate_projected_finite_width_chi0(request, finite_width_result)
            write(unit, '(a,1x,i0,1x,2(es24.16,1x),i0,1x,es24.16,1x,a)') &
               '# gf_controls q=', iq, gf_result%energy_min, gf_result%energy_max, &
               gf_result%integration_points, gf_result%integration_eta, 'quadrature=Simpson'
            norm_lehmann = sqrt(sum(abs(lehmann_result%susceptibility)**2))
            norm_gf = sqrt(sum(abs(gf_result%susceptibility)**2))
            norm_finite = sqrt(sum(abs(finite_width_result%susceptibility)**2))
            difference = sqrt(sum(abs(lehmann_result%susceptibility - gf_result%susceptibility)**2))
            relative = difference/max(norm_lehmann, tiny(1.0_rp))
            finite_difference = sqrt(sum(abs(finite_width_result%susceptibility - lehmann_result%susceptibility)**2))
            finite_relative = finite_difference/max(norm_lehmann, tiny(1.0_rp))
            write (*, '(a,a,a,i0,a,es12.4,a,es12.4,a,es12.4)') 'DRESP-02 Fe ', trim(projection), &
               ' q=', iq, ' norm_Lehmann=', norm_lehmann, ' norm_GF=', norm_gf, ' dF=', difference
            do ifrequency = 1, size(config%frequencies)
               do j = 1, contract%nsite
                  do i = 1, contract%nsite
                     write(unit, '(a,1x,i0,1x,es24.16,1x,2(i0,1x),6(es24.16,1x))') trim(projection), iq, &
                        config%frequencies(ifrequency), i, j, real(lehmann_result%susceptibility(i, j, ifrequency), rp), &
                        aimag(lehmann_result%susceptibility(i, j, ifrequency)), real(gf_result%susceptibility(i, j, ifrequency), rp), &
                        aimag(gf_result%susceptibility(i, j, ifrequency)), abs(lehmann_result%susceptibility(i, j, ifrequency) - &
                        gf_result%susceptibility(i, j, ifrequency)), relative
                     write(unit, '(a,1x,a,1x,i0,1x,es24.16,1x,2(i0,1x),16(es24.16,1x))') '# finite_width', trim(projection), iq, &
                        config%frequencies(ifrequency), i, j, real(lehmann_result%susceptibility(i, j, ifrequency), rp), &
                        aimag(lehmann_result%susceptibility(i, j, ifrequency)), real(finite_width_result%susceptibility(i, j, ifrequency), rp), &
                        aimag(finite_width_result%susceptibility(i, j, ifrequency)), real(gf_result%susceptibility(i, j, ifrequency), rp), &
                        aimag(gf_result%susceptibility(i, j, ifrequency)), real(gf_result%susceptibility(i, j, ifrequency) - &
                        finite_width_result%susceptibility(i, j, ifrequency), rp), aimag(gf_result%susceptibility(i, j, ifrequency) - &
                        finite_width_result%susceptibility(i, j, ifrequency)), abs(gf_result%susceptibility(i, j, ifrequency) - &
                        finite_width_result%susceptibility(i, j, ifrequency)), abs(gf_result%susceptibility(i, j, ifrequency) - &
                        finite_width_result%susceptibility(i, j, ifrequency))/max(abs(finite_width_result%susceptibility(i, j, ifrequency)), tiny(1.0_rp)), &
                        real(finite_width_result%susceptibility(i, j, ifrequency) - &
                        lehmann_result%susceptibility(i, j, ifrequency), rp), aimag(finite_width_result%susceptibility(i, j, ifrequency) - &
                        lehmann_result%susceptibility(i, j, ifrequency)), abs(finite_width_result%susceptibility(i, j, ifrequency) - &
                        lehmann_result%susceptibility(i, j, ifrequency)), abs(finite_width_result%susceptibility(i, j, ifrequency) - &
                        lehmann_result%susceptibility(i, j, ifrequency))/max(abs(lehmann_result%susceptibility(i, j, ifrequency)), tiny(1.0_rp))
                  end do
               end do
            end do
            if (iq == gamma_index) then
               gamma_lehmann = lehmann_result
               gamma_saved = .true.
            end if
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

         if (config%gf_closure_audit .and. gamma_saved) then
            call write_projected_gf_closure_audit(unit, projection, contract, product, left_state, &
               endpoints(gamma_index), gamma_lehmann, config, integration_eta)
         end if
      end do
      close(unit)
      deallocate(moment)
      deallocate(accepted_moment)
   end subroutine run_tddft_projected_chi0

   !> DRESP-02R material closure campaign.  Every sample uses the same
   !> accepted left state and exact folded Gamma endpoint; only the real-axis
   !> quadrature mesh, integration broadening, or finite energy window changes.
   subroutine write_projected_gf_closure_audit(unit, projection, contract, product, left_state, endpoint, &
                                               gamma_lehmann, config, integration_eta)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: projection
      type(projected_site_spin_contract), intent(in), target :: contract
      type(lmto_product_response_basis), intent(in), target :: product
      type(lr_electronic_state), intent(in), target :: left_state, endpoint
      type(projected_chi0_result), intent(in) :: gamma_lehmann
      type(tddft_production_config), intent(in) :: config
      real(rp), intent(in) :: integration_eta

      type(projected_chi0_request) :: request
      type(projected_chi0_result) :: gf_result, finite_width_result
      integer, parameter :: n_campaign = 3, n_eta_campaign = 5
      integer :: mesh_points(n_campaign), eta_points(n_eta_campaign), window_points, i, ik
      real(rp) :: eta_values(n_eta_campaign), margins(n_campaign), width, target_ratio
      real(rp) :: norm_lehmann, norm_gf, delta_re, delta_im, delta_abs, relative
      real(rp) :: sample_eta, sample_margin
      complex(rp), allocatable :: delta(:, :, :)

      if (size(config%frequencies) < 1) error stop 'DRESP-02R: no frequency available for GF audit'
      width = maxval(endpoint%eigenvalues) - minval(endpoint%eigenvalues) + 2.0_rp*config%gf_energy_margin
      width = max(width, maxval(left_state%eigenvalues) - minval(left_state%eigenvalues) + &
         2.0_rp*config%gf_energy_margin)
      target_ratio = 0.40_rp
      mesh_points = [max(101, config%gf_integration_points/4), max(201, config%gf_integration_points/2), &
         max(401, config%gf_integration_points)]
      do i = 1, n_campaign
         if (mod(mesh_points(i), 2) == 0) mesh_points(i) = mesh_points(i) + 1
      end do
      eta_values = integration_eta*[4.0_rp, 2.0_rp, 1.0_rp, 0.5_rp, 0.25_rp]
      do i = 1, n_eta_campaign
         eta_points(i) = odd_at_least(ceiling(width/(target_ratio*eta_values(i))) + 1)
      end do
      margins = [0.30_rp, 0.60_rp, 1.00_rp]
      window_points = odd_at_least(max(3, config%gf_integration_points))

      write(unit, '(a)') '# DRESP-02R controlled GF closure audit; all samples reuse one frozen accepted state'
      write(unit, '(a)') '# gf_audit_samples columns: projection campaign sample eta_response eta_int margin Emin Emax NE h h_over_eta chiL chiGF delta_Re delta_Im delta_abs relative finite_width_norm GF_minus_finite_Re GF_minus_finite_Im GF_minus_finite_abs GF_minus_finite_relative finite_width_minus_Lehmann_Re finite_width_minus_Lehmann_Im finite_width_minus_Lehmann_abs finite_width_minus_Lehmann_relative wall_seconds'
      write(unit, '(a)') '# gf_profile columns: projection campaign sample implementation allocation_s vertex_s endpoint_transform_s resolvent_s accumulator_s diagnostic_s wall_s resolvent_calls accumulator_calls'

      ! Fixed integration eta: this isolates the real-axis mesh error.
      do i = 1, n_campaign
         sample_eta = integration_eta
         sample_margin = config%gf_energy_margin
         call evaluate_projected_gf_audit_sample(gf_result, finite_width_result, contract, product, left_state, endpoint, config, &
            sample_eta, sample_margin, mesh_points(i), .true.)
         call write_projected_gf_audit_row(unit, projection, 'fixed_eta_mesh', i, gamma_lehmann, gf_result, &
            finite_width_result, config%eta, sample_eta, sample_margin)
         if (i == n_campaign) then
            write(unit, '(a,1x,a,1x,a,1x,7(es24.16,1x))') '# gf_spectral_moments', trim(projection), 'fine_mesh', &
               gf_result%left_spectral_zeroth_residual, gf_result%right_spectral_zeroth_residual, &
               gf_result%left_spectral_first_residual, gf_result%right_spectral_first_residual, &
               gf_result%left_spectral_fermi_residual, gf_result%right_spectral_fermi_residual, &
               gf_result%spacing_over_integration_eta
            write(unit, '(a,1x,a,1x,a,1x,7(es24.16,1x))') '# gf_kubo_terms', trim(projection), 'fine_mesh', &
               sqrt(sum(abs(gf_result%kubo_term_one)**2)), real(sum(gf_result%kubo_term_one), rp), &
               aimag(sum(gf_result%kubo_term_one)), sqrt(sum(abs(gf_result%kubo_term_two)**2)), &
               real(sum(gf_result%kubo_term_two), rp), aimag(sum(gf_result%kubo_term_two)), &
               sqrt(sum(abs(gf_result%kubo_term_one + gf_result%kubo_term_two - gf_result%susceptibility)**2))
            write(unit, '(a)') '# gf_k columns: projection k_index chiL_k chiGF_k delta_Re delta_Im delta_abs relative chiGF_Re chiGF_Im'
            do ik = 1, left_state%nk
               call write_projected_gf_k_row(unit, projection, ik, gamma_lehmann, gf_result)
            end do
            write(unit, '(a)') '# dominant_lehmann columns: projection rank k_index left_band right_band left_E right_E left_f right_f transition_E matrix_element_weight score'
            do ik = 1, gamma_lehmann%n_dominant_transition_records
               write(unit, '(a,1x,a,1x,i0,1x,3(i0,1x),7(es24.16,1x))') '# dominant_lehmann', trim(projection), ik, &
                  gamma_lehmann%dominant_transitions(ik)%k_index, gamma_lehmann%dominant_transitions(ik)%left_band, &
                  gamma_lehmann%dominant_transitions(ik)%right_band, gamma_lehmann%dominant_transitions(ik)%left_energy, &
                  gamma_lehmann%dominant_transitions(ik)%right_energy, gamma_lehmann%dominant_transitions(ik)%left_occupation, &
                  gamma_lehmann%dominant_transitions(ik)%right_occupation, gamma_lehmann%dominant_transitions(ik)%transition_energy, &
                  gamma_lehmann%dominant_transitions(ik)%matrix_element_weight, gamma_lehmann%dominant_transitions(ik)%score
            end do
         end if
      end do

      ! Integration eta ladder with h/eta held near target_ratio.  This is
      ! independent of the fixed-eta mesh campaign above.
      do i = 1, n_eta_campaign
         if (eta_values(i) >= config%eta) cycle
         call evaluate_projected_gf_audit_sample(gf_result, finite_width_result, contract, product, left_state, endpoint, config, &
            eta_values(i), config%gf_energy_margin, eta_points(i), .false.)
         call write_projected_gf_audit_row(unit, projection, 'controlled_eta', i, gamma_lehmann, gf_result, &
            finite_width_result, config%eta, eta_values(i), config%gf_energy_margin)
      end do

      ! Window ladder: broadening and mesh count remain fixed while the
      ! finite spectral interval is changed explicitly.
      do i = 1, n_campaign
         call evaluate_projected_gf_audit_sample(gf_result, finite_width_result, contract, product, left_state, endpoint, config, &
            integration_eta, margins(i), window_points, .false.)
         call write_projected_gf_audit_row(unit, projection, 'energy_window', i, gamma_lehmann, gf_result, &
            finite_width_result, config%eta, integration_eta, margins(i))
      end do

      deallocate(gf_result%susceptibility)
      deallocate(finite_width_result%susceptibility)
   end subroutine write_projected_gf_closure_audit

   subroutine evaluate_projected_gf_audit_sample(result, finite_width_result, contract, product, left_state, endpoint, config, &
                                                 integration_eta, margin, integration_points, diagnostics)
      type(projected_chi0_result), intent(out) :: result
      type(projected_chi0_result), intent(out) :: finite_width_result
      type(projected_site_spin_contract), intent(in), target :: contract
      type(lmto_product_response_basis), intent(in), target :: product
      type(lr_electronic_state), intent(in), target :: left_state, endpoint
      type(tddft_production_config), intent(in) :: config
      real(rp), intent(in) :: integration_eta, margin
      integer, intent(in) :: integration_points
      logical, intent(in) :: diagnostics
      type(projected_chi0_request) :: request

      request%q = [0.0_rp, 0.0_rp, 0.0_rp]
      request%frequencies = config%frequencies
      request%eta = config%eta
      request%channel = config%channel
      request%integration_points = integration_points
      request%integration_eta = integration_eta
      request%energy_margin = margin
      request%diagnostics = diagnostics
      request%contract => contract
      request%product_basis => product
      request%electronic_state => left_state
      request%q_endpoint_state => endpoint
      call evaluate_projected_gf_chi0(request, result)
      request%diagnostics = .false.
      call evaluate_projected_finite_width_chi0(request, finite_width_result)
   end subroutine evaluate_projected_gf_audit_sample

   subroutine write_projected_gf_audit_row(unit, projection, campaign, sample, gamma_lehmann, gf_result, &
                                           finite_width_result, response_eta, integration_eta, margin)
      integer, intent(in) :: unit, sample
      character(len=*), intent(in) :: projection, campaign
      type(projected_chi0_result), intent(in) :: gamma_lehmann, gf_result, finite_width_result
      real(rp), intent(in) :: response_eta, integration_eta, margin
      complex(rp), allocatable :: delta(:, :, :)
      complex(rp), allocatable :: finite_delta(:, :, :), oracle_target_delta(:, :, :)
      real(rp) :: norm_lehmann, norm_gf, norm_finite, delta_re, delta_im, delta_abs, relative
      real(rp) :: finite_delta_re, finite_delta_im, finite_delta_abs, finite_relative
      real(rp) :: oracle_target_re, oracle_target_im, oracle_target_abs, oracle_target_relative

      allocate(delta, mold=gamma_lehmann%susceptibility)
      allocate(finite_delta, mold=gamma_lehmann%susceptibility)
      allocate(oracle_target_delta, mold=gamma_lehmann%susceptibility)
      delta = gamma_lehmann%susceptibility - gf_result%susceptibility
      finite_delta = gf_result%susceptibility - finite_width_result%susceptibility
      oracle_target_delta = finite_width_result%susceptibility - gamma_lehmann%susceptibility
      norm_lehmann = sqrt(sum(abs(gamma_lehmann%susceptibility)**2))
      norm_gf = sqrt(sum(abs(gf_result%susceptibility)**2))
      norm_finite = sqrt(sum(abs(finite_width_result%susceptibility)**2))
      delta_re = sqrt(sum(real(delta, rp)**2))
      delta_im = sqrt(sum(aimag(delta)**2))
      delta_abs = sqrt(sum(abs(delta)**2))
      relative = delta_abs/max(norm_lehmann, tiny(1.0_rp))
      finite_delta_re = sqrt(sum(real(finite_delta, rp)**2))
      finite_delta_im = sqrt(sum(aimag(finite_delta)**2))
      finite_delta_abs = sqrt(sum(abs(finite_delta)**2))
      finite_relative = finite_delta_abs/max(norm_finite, tiny(1.0_rp))
      oracle_target_re = sqrt(sum(real(oracle_target_delta, rp)**2))
      oracle_target_im = sqrt(sum(aimag(oracle_target_delta)**2))
      oracle_target_abs = sqrt(sum(abs(oracle_target_delta)**2))
      oracle_target_relative = oracle_target_abs/max(norm_lehmann, tiny(1.0_rp))
      write(unit, '(a,1x,a,1x,a,1x,i0,1x,5(es24.16,1x),i0,1x,18(es24.16,1x))') '# gf_audit', trim(projection), trim(campaign), sample, &
         response_eta, integration_eta, margin, gf_result%energy_min, gf_result%energy_max, gf_result%integration_points, &
         gf_result%energy_spacing, gf_result%spacing_over_integration_eta, norm_lehmann, norm_gf, delta_re, delta_im, &
         delta_abs, relative, norm_finite, finite_delta_re, finite_delta_im, finite_delta_abs, finite_relative, &
         oracle_target_re, oracle_target_im, oracle_target_abs, oracle_target_relative, gf_result%wall_time_seconds
      write(unit, '(a,1x,a,1x,a,1x,i0,1x,a,1x,7(es24.16,1x),2(i0,1x))') '# gf_profile', trim(projection), trim(campaign), sample, &
         trim(gf_result%implementation), gf_result%allocation_seconds, gf_result%vertex_seconds, &
         gf_result%endpoint_transform_seconds, gf_result%resolvent_seconds, gf_result%accumulator_seconds, &
         gf_result%diagnostic_seconds, gf_result%wall_time_seconds, gf_result%resolvent_calls, gf_result%accumulator_calls
      write(unit, '(a,1x,a,1x,a,1x,i0,1x,a,1x,4(es24.16,1x))') '# finite_width_profile', trim(projection), trim(campaign), sample, &
         trim(finite_width_result%implementation), finite_width_result%wall_time_seconds, finite_width_result%energy_spacing, &
         finite_width_result%spacing_over_integration_eta, real(finite_width_result%ntransitions_evaluated, rp)
      deallocate(delta)
      deallocate(finite_delta, oracle_target_delta)
   end subroutine write_projected_gf_audit_row

   subroutine write_projected_gf_k_row(unit, projection, k_index, gamma_lehmann, gf_result)
      integer, intent(in) :: unit, k_index
      character(len=*), intent(in) :: projection
      type(projected_chi0_result), intent(in) :: gamma_lehmann, gf_result
      complex(rp), allocatable :: delta(:, :, :)
      real(rp) :: norm_lehmann, norm_gf, delta_re, delta_im, delta_abs, relative

      allocate(delta, mold=gamma_lehmann%k_susceptibility(:, :, :, k_index))
      delta = gamma_lehmann%k_susceptibility(:, :, :, k_index) - gf_result%k_susceptibility(:, :, :, k_index)
      norm_lehmann = sqrt(sum(abs(gamma_lehmann%k_susceptibility(:, :, :, k_index))**2))
      norm_gf = sqrt(sum(abs(gf_result%k_susceptibility(:, :, :, k_index))**2))
      delta_re = sqrt(sum(real(delta, rp)**2))
      delta_im = sqrt(sum(aimag(delta)**2))
      delta_abs = sqrt(sum(abs(delta)**2))
      relative = delta_abs/max(norm_lehmann, tiny(1.0_rp))
      write(unit, '(a,1x,a,1x,i0,1x,8(es24.16,1x))') '# gf_k', trim(projection), k_index, norm_lehmann, norm_gf, &
         delta_re, delta_im, delta_abs, relative, real(gf_result%k_susceptibility(1, 1, 1, k_index), rp), &
         aimag(gf_result%k_susceptibility(1, 1, 1, k_index))
      deallocate(delta)
   end subroutine write_projected_gf_k_row

   pure integer function odd_at_least(value) result(odd_value)
      integer, intent(in) :: value
      odd_value = max(3, value)
      if (mod(odd_value, 2) == 0) odd_value = odd_value + 1
   end function odd_at_least

   !> Evaluate the compact direct-ALSDA Dyson response.
   !>
   !> This is the production worker for one accepted state.  The response
   !> services are evaluated q-by-q and serialized while their q-local result
   !> objects are live.  Validation audits are explicit opt-ins and add work
   !> around this worker; they do not define its production input contract.
   subroutine run_tddft_compact_dyson(config, response_space, radial_bases, ground_states, left_state, endpoints, &
                                      reciprocal_obj, lattice_obj, control_obj)
      type(tddft_production_config), intent(in) :: config
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
      type(lr_product_gf_susceptibility_request) :: gf_request
      type(lr_product_gf_susceptibility_result) :: gf_result
      type(lr_product_ks_susceptibility_request) :: gf_lehmann_request
      type(lr_product_ks_susceptibility_result) :: gf_lehmann_result
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
      real(rp) :: gf_frequency, gf_norm_lehmann, gf_norm_gf, gf_d_frobenius, gf_relative_frobenius, gf_d_infinity
      real(rp) :: gf_integration_eta, gf_spacing_ratio, gf_energy_min, gf_energy_max, gf_energy_spacing, gf_wall_seconds
      integer :: ndim, nfrequency, nq, iq, iw, i, j, unit, gamma_index, positive_q_index, negative_q_index
      integer :: gf_frequency_index, gf_q_index, gf_integration_points
      character(len=256) :: state_file
      logical :: static_audit_enabled, covariance_enabled, gf_audit_enabled

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
      gf_audit_enabled = config%gf_closure_audit
      gamma_index = find_gamma_q_index(config%q_list)
      if (static_audit_enabled .and. gamma_index == 0) then
         error stop 'compact_dyson: dyson_static_audit requires Gamma; input should have been rejected during preflight'
      end if
      if (gf_audit_enabled .and. gamma_index == 0) then
         error stop 'compact_dyson: gf_closure_audit requires Gamma; input should have been rejected during preflight'
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
         write(unit, '(a,l1)') '# gf_closure_audit = ', gf_audit_enabled
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
         if (gf_audit_enabled) then
            write(unit, '(a)') '# gf_spots columns: q_index qx qy qz omega_Ry integration_points integration_eta spacing_over_eta energy_min_Ry energy_max_Ry h_Ry norm_lehmann norm_gf dF rF dInf wall_seconds'
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

      if (gf_audit_enabled) then
         gf_q_index = gamma_index
         gf_frequency_index = 1
         do iw = 1, nfrequency
            if (abs(config%frequencies(iw)) > 1.0e-12_rp) then
               gf_frequency_index = iw
               exit
            end if
         end do
         gf_frequency = config%frequencies(gf_frequency_index)
         gf_request%q = config%q_list(:, gf_q_index)
         gf_request%frequencies = [gf_frequency]
         gf_request%eta = config%eta
         gf_request%channel = config%channel
         gf_request%integration_points = config%gf_integration_points
         gf_request%integration_eta = config%gf_integration_eta
         gf_request%energy_margin = config%gf_energy_margin
         gf_request%product_basis => product
         gf_request%electronic_state => left_state
         gf_request%q_endpoint_state => endpoints(gf_q_index)
         call evaluate_lr_product_gf_susceptibility(gf_request, gf_result)
         gf_lehmann_request%q = gf_request%q
         gf_lehmann_request%frequencies = gf_request%frequencies
         gf_lehmann_request%eta = config%eta
         gf_lehmann_request%channel = config%channel
         gf_lehmann_request%product_basis => product
         gf_lehmann_request%electronic_state => left_state
         gf_lehmann_request%q_endpoint_state => endpoints(gf_q_index)
         call evaluate_lr_product_ks_susceptibility(gf_lehmann_request, gf_lehmann_result)
         gf_norm_lehmann = sqrt(sum(abs(gf_lehmann_result%susceptibility(:, :, 1))**2))
         gf_norm_gf = sqrt(sum(abs(gf_result%susceptibility(:, :, 1))**2))
         gf_d_frobenius = sqrt(sum(abs(gf_lehmann_result%susceptibility(:, :, 1) - gf_result%susceptibility(:, :, 1))**2))
         gf_relative_frobenius = gf_d_frobenius/max(gf_norm_lehmann, gf_norm_gf, tiny(1.0_rp))
         gf_d_infinity = maxval(abs(gf_lehmann_result%susceptibility(:, :, 1) - gf_result%susceptibility(:, :, 1)))
         gf_integration_eta = gf_result%actual_integration_eta
         gf_spacing_ratio = gf_result%spacing_over_integration_eta
         gf_energy_min = gf_result%energy_min
         gf_energy_max = gf_result%energy_max
         gf_energy_spacing = gf_result%energy_spacing
         gf_wall_seconds = gf_result%wall_time_seconds
         gf_integration_points = gf_result%integration_points
         if (rank == 0) then
            write(unit, '(i0,1x,4(es24.16,1x),i0,1x,11(es24.16,1x))') gf_q_index, config%q_list(:, gf_q_index), &
               gf_frequency, gf_integration_points, gf_integration_eta, gf_spacing_ratio, gf_energy_min, gf_energy_max, &
               gf_energy_spacing, gf_norm_lehmann, gf_norm_gf, gf_d_frobenius, gf_relative_frobenius, gf_d_infinity, gf_wall_seconds
         end if
      end if

      if (rank == 0) close(unit)
      write(*, '(a,i0,a,i0,a,i0)') 'TDDFT compact Dyson response: q_count=', nq, ' omega_count=', nfrequency, &
         ' product_dimension=', ndim
   end subroutine run_tddft_compact_dyson

   !> TDVK-06 static interaction diagnostics.  This backend consumes the
   !> accepted k-space SCF snapshot exactly once, evaluates the complete
   !> compact Gamma/omega=0 response, and stops after raw ALSDA and independent
   !> GSR diagnostics.  It never enters Dyson, loss, mode extraction, or any
   !> Goldstone repair/correction path.
   subroutine run_tddft_static_interactions(config, response_space, radial_bases, ground_states, left_state, endpoints, &
                                            reciprocal_obj, lattice_obj)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj

      type(lmto_product_response_basis), target :: product_plus, product_minus
      type(lmto_product_response_basis), pointer :: product
      type(lr_product_ks_susceptibility_request) :: request
      type(lr_product_ks_susceptibility_result) :: response_result
      type(lr_alsda_kernel_result) :: kernel_result
      type(lr_compact_gsr_result) :: gsr_result
      real(rp), allocatable :: magnetization(:, :), magnetization_valence(:, :), magnetization_core(:, :)
      complex(rp), allocatable :: magnetization_compact(:), magnetization_point(:), magnetization_projected_point(:)
      complex(rp), allocatable :: valence_compact(:), valence_point(:), valence_projected_point(:)
      complex(rp), allocatable :: core_compact(:), core_point(:), core_projected_point(:)
      complex(rp), allocatable :: compact_kernel(:, :), field(:), response(:), residual(:)
      complex(rp), allocatable :: point_residual(:)
      type(lr_alsda_kernel_request) :: kxc_request
      real(rp) :: direct_norm, direct_relative, target_norm, state_mesh_max, state_weight_max, state_ef_diff
      real(rp) :: state_eigen_max, state_occ_max, state_projector_max, state_projector_frobenius
      real(rp) :: k_fingerprint(5), accepted_moment, weight_sum
      real(rp) :: kernel_min, kernel_max, kernel_max_abs, runtime_start, runtime_end
      real(rp) :: magnetization_weighted_norm, magnetization_projection_norm, magnetization_projection_relative
      real(rp) :: valence_weighted_norm, valence_projection_norm, valence_projection_relative
      real(rp) :: core_weighted_norm, core_projection_norm, core_projection_relative
      real(rp) :: gsr_runtime_start, gsr_runtime_end
      complex(rp) :: rigid_overlap, trace
      integer :: gamma_index, ieta, ifrequency, unit, site, response_l, response_m, ir, flat, ik
      type(response_super_index) :: item
      real(rp) :: radial_residual_norm
      logical :: rank_stable, finite_response
      character(len=512) :: state_file

      if (.not. reciprocal_obj%auto_find_fermi .or. .not. reciprocal_obj%canonical_energy_valid) then
         error stop 'TDVK-06 static interactions: accepted auto-EF k-space SCF state is required'
      end if
      if (size(endpoints) /= size(config%q_list, 2)) then
         error stop 'TDVK-06 static interactions: q endpoint count differs from q_list'
      end if
      if (size(endpoints) /= 1 .or. sum(abs(config%q_list(:, 1))) > 1.0e-12_rp) then
         error stop 'TDVK-06 static interactions: the sole endpoint must be Gamma'
      end if
      gamma_index = 1
      if (size(config%frequencies) /= 1 .or. abs(config%frequencies(1)) > 1.0e-12_rp) then
         error stop 'TDVK-06 static interactions: exactly omega=0 is required'
      end if

      call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      call product_minus%initialize(response_space, radial_bases, lmto_product_channel_minus, .true.)
      if (trim(config%channel) == 'chi_plus') then
         product => product_plus
      else
         product => product_minus
      end if
      rank_stable = .true.
      do site = 1, product%nsite
         do response_l = 0, product%response_lmax
            rank_stable = rank_stable .and. product%blocks(site, response_l)%rank_stable
         end do
      end do
      if (.not. rank_stable .or. product%product_dimension /= 232) then
         error stop 'TDVK-06 static interactions: complete compact spd product representation is not certified'
      end if

      call compute_accepted_pauli_magnetization(reciprocal_obj, lattice_obj%symbolic_atoms, lattice_obj%nbulk, &
         magnetization, magnetization_valence, magnetization_core)
      allocate(magnetization_compact(product%product_dimension), valence_compact(product%product_dimension), &
         core_compact(product%product_dimension), magnetization_projected_point(response_space%ndim), &
         valence_projected_point(response_space%ndim), core_projected_point(response_space%ndim))
      call compact_project_magnetization(response_space, product, magnetization, magnetization_compact, magnetization_point)
      call compact_reconstruct_point_vector(response_space, product, magnetization_compact, magnetization_projected_point)
      call compact_weighted_projection_diagnostics(response_space, magnetization_point, magnetization_projected_point, &
         magnetization_weighted_norm, magnetization_projection_norm, magnetization_projection_relative)
      call compact_project_magnetization(response_space, product, magnetization_valence, valence_compact, valence_point)
      call compact_reconstruct_point_vector(response_space, product, valence_compact, valence_projected_point)
      call compact_weighted_projection_diagnostics(response_space, valence_point, valence_projected_point, &
         valence_weighted_norm, valence_projection_norm, valence_projection_relative)
      call compact_project_magnetization(response_space, product, magnetization_core, core_compact, core_point)
      call compact_reconstruct_point_vector(response_space, product, core_compact, core_projected_point)
      call compact_weighted_projection_diagnostics(response_space, core_point, core_projected_point, &
         core_weighted_norm, core_projection_norm, core_projection_relative)
      target_norm = sqrt(sum(abs(magnetization_compact)**2))
      if (target_norm <= tiny(1.0_rp)) then
         error stop 'TDVK-06 static interactions: accepted Pauli magnetization has no retained compact component'
      end if
      call prepare_direct_alsda_request(response_space, ground_states, magnetization, kxc_request)
      call evaluate_lr_alsda_kernel(kxc_request, kernel_result)
      allocate(compact_kernel(product%product_dimension, product%product_dimension), field(product%product_dimension), &
         response(product%product_dimension), residual(product%product_dimension), point_residual(response_space%ndim))
      call compact_project_local_operator(response_space, product, cmplx(kernel_result%pointwise_kernel, 0.0_rp, rp), &
         compact_kernel)
      kernel_min = minval(kernel_result%pointwise_kernel)
      kernel_max = maxval(kernel_result%pointwise_kernel)
      kernel_max_abs = maxval(abs(kernel_result%pointwise_kernel))
      state_file = trim(config%output_file)//'.state'
      call state_consistency_metrics(reciprocal_obj, left_state, state_mesh_max, state_weight_max, state_ef_diff, &
         state_eigen_max, state_occ_max, state_projector_max, state_projector_frobenius)

      accepted_moment = 0.0_rp
      do site = 1, size(ground_states)
         accepted_moment = accepted_moment + ground_states(site)%integrated_moment_muB
      end do
      weight_sum = sum(left_state%k_weights)
      k_fingerprint = 0.0_rp
      do ik = 1, left_state%nk
         k_fingerprint(1) = k_fingerprint(1) + left_state%k_weights(ik)
         k_fingerprint(2:4) = k_fingerprint(2:4) + left_state%k_weights(ik)*left_state%k_points(:, ik)
         k_fingerprint(5) = k_fingerprint(5) + real(ik, rp)*left_state%k_weights(ik)
      end do

      if (rank == 0) then
         open(newunit=unit, file=trim(config%output_file), status='replace', action='write')
         write(unit, '(a)') '# TDVK-06 Fe static ALSDA and independent GSR diagnostics'
         write(unit, '(a,a)') '# build_version = ', trim(tddft_build_version)
         write(unit, '(a)') '# backend = static_interactions'
         write(unit, '(a)') '# scope = accepted self-consistent k-space SCF -> same-state TDDFT; static Gamma only'
         write(unit, '(a)') '# no Dyson, denominator diagnostics, spectrum, mode extraction, Goldstone repair, BES, GCR, rescale, shift, or sign tuning'
         write(unit, '(a,a)') '# state_source = accepted_kspace_scf_cache'
         write(unit, '(a)') '# reciprocal_rebuild_performed_for_tddft = F'
         write(unit, '(a)') '# reciprocal_hamiltonian_rebuild_performed_for_tddft = F'
         write(unit, '(a)') '# accepted_state_cache_reused = T'
         write(unit, '(a,3(i0,1x))') '# requested_k_mesh = ', reciprocal_obj%nk_mesh
         write(unit, '(a,3(i0,1x))') '# actual_k_mesh = ', reciprocal_obj%nk_mesh
         write(unit, '(a,i0)') '# actual_k_count = ', left_state%nk
         write(unit, '(a,es24.16)') '# accepted_state_EF_Ry = ', left_state%fermi_level
         write(unit, '(a,es24.16)') '# accepted_state_temperature_K = ', left_state%temperature
         write(unit, '(a,es24.16)') '# accepted_state_integrated_electron_count = ', reciprocal_obj%canonical_electron_count
         write(unit, '(a,es24.16)') '# accepted_state_target_electron_count = ', reciprocal_obj%total_electrons
         write(unit, '(a,es24.16)') '# accepted_state_integrated_moment_muB = ', accepted_moment
         write(unit, '(a,es24.16)') '# accepted_state_radial_residual_control = ', ground_states(1)%accepted_residual_control
         write(unit, '(a,5(es24.16,1x))') '# k_fingerprint_checksums = ', k_fingerprint
         write(unit, '(a,es24.16)') '# k_weight_sum = ', weight_sum
         write(unit, '(a,es24.16)') '# state_consistency_mesh_max_abs = ', state_mesh_max
         write(unit, '(a,es24.16)') '# state_consistency_weight_max_abs = ', state_weight_max
         write(unit, '(a,es24.16)') '# state_consistency_EF_max_abs_Ry = ', state_ef_diff
         write(unit, '(a,es24.16)') '# state_consistency_eigenvalue_max_abs_Ry = ', state_eigen_max
         write(unit, '(a,es24.16)') '# state_consistency_occupation_max_abs = ', state_occ_max
         write(unit, '(a,es24.16)') '# state_consistency_projector_max_abs = ', state_projector_max
         write(unit, '(a,es24.16)') '# state_consistency_projector_frobenius = ', state_projector_frobenius
         write(unit, '(a,a)') '# state_artifact = ', trim(state_file)
         write(unit, '(a,a)') '# response_backend = complete compact Lehmann on accepted immutable k+q endpoints'
         write(unit, '(a,a)') '# response_representation = ', trim(lr_compact_representation)
         write(unit, '(a,a)') '# compact_mapping_contract = ', trim(lr_compact_mapping_contract)
         write(unit, '(a)') '# compact_representation_certification = PASS: independent point/product projection oracle and strict rank-stable 232-mode runtime inventory'
         write(unit, '(a,i0)') '# product_unpruned_dimension = ', product%unpruned_dimension
         write(unit, '(a,i0)') '# product_dimension = ', product%product_dimension
         write(unit, '(a,i0)') '# response_lmax = ', response_space%response_lmax
         write(unit, '(a,a)') '# channel = ', trim(config%channel)
         write(unit, '(a,es24.16)') '# omega_Ry = ', config%frequencies(1)
         write(unit, '(a,a)') '# accepted_xc_functional = ', trim(kernel_result%xc_provenance%functional_name)
         write(unit, '(a,a)') '# accepted_xc_backend = ', trim(kernel_result%xc_provenance%backend_name)
         write(unit, '(a,a)') '# accepted_xc_mapping_quality = ', trim(kernel_result%xc_provenance%mapping_quality)
         write(unit, '(a,i0)') '# accepted_xc_txc = ', kernel_result%xc_provenance%txc
         write(unit, '(a,a)') '# pauli_magnetization_label = ', trim(kernel_result%magnetization_label)
         write(unit, '(a)') '# pauli_magnetization_provenance = occupations deterministically reconstructed from the accepted reciprocal state using the same EF, temperature and Fermi function; accepted eigenvectors + POTPAR large-component and frozen-core projection'
         write(unit, '(a,es24.16)') '# pauli_magnetization_compact_norm = ', target_norm
         write(unit, '(a,es24.16)') '# pauli_magnetization_weighted_norm_m00 = ', magnetization_weighted_norm
         write(unit, '(a,es24.16)') '# pauli_magnetization_projection_residual_norm_m00 = ', magnetization_projection_norm
         write(unit, '(a,es24.16)') '# pauli_magnetization_projection_relative_residual_m00 = ', magnetization_projection_relative
         write(unit, '(a,es24.16)') '# pauli_valence_weighted_norm_m00 = ', valence_weighted_norm
         write(unit, '(a,es24.16)') '# pauli_valence_projection_residual_norm_m00 = ', valence_projection_norm
         write(unit, '(a,es24.16)') '# pauli_valence_projection_relative_residual_m00 = ', valence_projection_relative
         write(unit, '(a,es24.16)') '# pauli_core_weighted_norm_m00 = ', core_weighted_norm
         write(unit, '(a,es24.16)') '# pauli_core_projection_residual_norm_m00 = ', core_projection_norm
         write(unit, '(a,es24.16)') '# pauli_core_projection_relative_residual_m00 = ', core_projection_relative
         write(unit, '(a,i0)') '# pauli_magnetization_point_sites = ', size(magnetization, 1)
         write(unit, '(a,i0)') '# pauli_magnetization_point_radial_points = ', size(magnetization, 2)
         write(unit, '(a,a)') '# direct_alsda_kernel_formula = Bxc_sigma/m_pauli with Bxc_sigma=(Vxc_up-Vxc_down)/2'
         write(unit, '(a,es24.16)') '# direct_alsda_kernel_min_Ry_bohr3 = ', kernel_min
         write(unit, '(a,es24.16)') '# direct_alsda_kernel_max_Ry_bohr3 = ', kernel_max
         write(unit, '(a,es24.16)') '# direct_alsda_kernel_max_abs_Ry_bohr3 = ', kernel_max_abs
         write(unit, '(a,a)') '# direct_alsda_operator_mapping = compact K=U^H K_point U; vector metric factors occur only in c=U^H sqrt(W)x'
         write(unit, '(a)') '# direct_alsda_columns = eta_Ry residual_norm residual_relative rigid_overlap_real rigid_overlap_imag field_norm response_norm finite status'
         write(unit, '(a)') '# direct_alsda_residual_definition = chiKS_compact(0,eta) * Kxc_compact * m00_compact - m00_compact'
         write(unit, '(a)') '# direct_alsda_residual_by_block_columns = eta_Ry site response_l residual_block_norm'
         write(unit, '(a)') '# direct_alsda_residual_by_radial_columns = eta_Ry site response_l radial_index residual_point_norm'
         write(unit, '(a,a)') '# compact_gsr_action_contract = ', &
            'f_j=Kc_j*c=P*K_j*R*c; Gamma_j=chiKS_compact*f_j; Kc(u)*c=sum_j u_j*f_j'
         write(unit, '(a,es24.16)') '# compact_gsr_action_consistency_tolerance = ', lr_compact_gsr_action_tolerance
         write(unit, '(a)') '# gsr_solve_policy = ZGELSS with RCOND=-1 (machine precision); no regularization, singular-value tuning, rescaling, constraint, or zero-mode enforcement'
         write(unit, '(a)') '# gsr_residual_difference_relative_scale = max(equation_rhs_norm, assembled_action_norm)'
         write(unit, '(a)') '# gsr_assembled_action_columns = sum_action_norm matrix_action_norm absolute_difference relative_difference max_component_difference'
         write(unit, '(a)') '# gsr_columns = eta_Ry compact_dimension equation_rows unknowns rank rank_deficient singular_min singular_max condition svd_rcond svd_cutoff coefficient_norm solve_residual_norm solve_relative reconstructed_residual_norm reconstructed_relative residual_difference_norm residual_difference_relative residual_difference_max max_residual_component rigid_overlap_real rigid_overlap_imag runtime_cpu_seconds blocked'
      end if

      do ieta = 1, size(config%eta_values)
         request%q = config%q_list(:, gamma_index)
         request%frequencies = config%frequencies
         request%eta = config%eta_values(ieta)
         request%channel = config%channel
         request%product_basis => product
         request%electronic_state => left_state
         request%q_endpoint_state => endpoints(gamma_index)
         call cpu_time(runtime_start)
         call evaluate_lr_product_ks_susceptibility(request, response_result)
         call cpu_time(runtime_end)
         finite_response = all(ieee_is_finite(real(response_result%susceptibility, rp))) .and. &
            all(ieee_is_finite(aimag(response_result%susceptibility)))
         if (.not. finite_response) error stop 'TDVK-06 static interactions: compact static chiKS contains NaN or Inf'
         field = matmul(compact_kernel, magnetization_compact)
         response = matmul(response_result%susceptibility(:, :, 1), field)
         residual = response - magnetization_compact
         call compact_reconstruct_point_vector(response_space, product, residual, point_residual)
         direct_norm = sqrt(sum(abs(residual)**2))
         direct_relative = direct_norm/target_norm
         rigid_overlap = sum(conjg(residual)*magnetization_compact)
         trace = cmplx(0.0_rp, 0.0_rp, rp)
         do flat = 1, product%product_dimension
            trace = trace + response_result%susceptibility(flat, flat, 1)
         end do
         call cpu_time(gsr_runtime_start)
         call evaluate_compact_goldstone_sumrule(response_space, product, response_result%susceptibility(:, :, 1), &
            magnetization, gsr_result)
         call cpu_time(gsr_runtime_end)
         if (rank == 0) then
            write(unit, '(7(es24.16,1x),l1,1x,a)') config%eta_values(ieta), direct_norm, &
               direct_relative, real(rigid_overlap, rp), aimag(rigid_overlap), sqrt(sum(abs(field)**2)), &
               sqrt(sum(abs(response)**2)), finite_response, 'EXECUTED'
            do site = 1, product%nsite
               do response_l = 0, product%response_lmax
                  residual_norm_by_block: block
                     complex(rp), allocatable :: block_values(:)
                     integer :: block_count, response_m, mode
                     block_count = (2*response_l + 1)*product%blocks(site, response_l)%rank
                     allocate(block_values(block_count))
                     block_values = cmplx(0.0_rp, 0.0_rp, rp)
                     do response_m = -response_l, response_l
                        do mode = 1, product%blocks(site, response_l)%rank
                           flat = product%flat_index(site, response_l, response_m, mode)
                           block_values((response_m + response_l)*product%blocks(site, response_l)%rank + mode) = residual(flat)
                        end do
                     end do
                     write(unit, '(es24.16,1x,2(i0,1x),es24.16)') config%eta_values(ieta), site, response_l, &
                        sqrt(sum(abs(block_values)**2))
                     deallocate(block_values)
                  end block residual_norm_by_block
               end do
            end do
            do site = 1, product%nsite
               do response_l = 0, product%response_lmax
                  do ir = 2, response_space%npoint
                     radial_residual_norm = 0.0_rp
                     do response_m = -response_l, response_l
                        item = response_super_index(site, response_l, response_m, ir, 1)
                        call response_flatten_superindex(item, response_space%nsite, response_space%response_lmax, &
                           response_space%npoint, response_space%nchannel, flat)
                        radial_residual_norm = radial_residual_norm + abs(point_residual(flat))**2
                     end do
                     write(unit, '(es24.16,1x,3(i0,1x),es24.16)') config%eta_values(ieta), site, response_l, ir, &
                        sqrt(radial_residual_norm)
                  end do
               end do
            end do
            write(unit, '(a,es24.16,1x,a,es24.16)') '# compact_chiKS_frobenius_eta_Ry = ', &
               sqrt(sum(abs(response_result%susceptibility(:, :, 1))**2)), '# runtime_cpu_seconds = ', runtime_end-runtime_start
            write(unit, '(a,es24.16,1x,es24.16)') '# compact_chiKS_trace_real_imag = ', real(trace, rp), aimag(trace)
            write(unit, '(es24.16,1x,4(i0,1x),l1,1x,17(es24.16,1x),l1)') config%eta_values(ieta), &
               product%product_dimension, gsr_result%equation_rows, gsr_result%unknowns, gsr_result%rank, &
               gsr_result%rank_deficient, minval(gsr_result%singular_values), maxval(gsr_result%singular_values), &
               gsr_result%condition_number, gsr_result%svd_rcond, gsr_result%svd_cutoff, gsr_result%coefficient_norm, &
               gsr_result%equation_residual_norm, gsr_result%equation_relative_residual, gsr_result%residual_norm, &
               gsr_result%relative_residual, gsr_result%residual_difference_norm, gsr_result%residual_difference_relative, &
               gsr_result%residual_difference_max_component, gsr_result%max_residual_component, &
               real(gsr_result%rigid_overlap, rp), aimag(gsr_result%rigid_overlap), gsr_runtime_end-gsr_runtime_start, &
               gsr_result%blocked
            write(unit, '(a,5(es24.16,1x))') '# gsr_assembled_action_norms = ', gsr_result%assembled_action_sum_norm, &
               gsr_result%assembled_action_matrix_norm, gsr_result%assembled_action_difference_norm, &
               gsr_result%assembled_action_difference_relative, gsr_result%assembled_action_difference_max_component
            write(unit, '(a,es24.16,1x,a,es24.16)') '# gsr_singular_values_min_max = ', minval(gsr_result%singular_values), &
               '#', maxval(gsr_result%singular_values)
            write(unit, '(a,a)') '# gsr_status = ', trim(gsr_result%status)
         end if
         if (allocated(response_result%susceptibility)) deallocate(response_result%susceptibility)
         if (allocated(response_result%frequencies)) deallocate(response_result%frequencies)
         if (allocated(gsr_result%singular_values)) deallocate(gsr_result%singular_values)
      end do
      if (rank == 0) then
         write(unit, '(a)') '# denominator_diagnostics = NOT PERFORMED: compact Dyson representation is not certified; leave for TDVK-07'
         write(unit, '(a)') '# TDVK-06 direct ALSDA static diagnostic = EXECUTED'
         write(unit, '(a)') '# TDVK-06 independent GSR raw solve = EXECUTED; status is reported without repair or correction'
         write(unit, '(a)') '# TDVK-06 compact representation certification = PASS'
         write(unit, '(a)') '# TDVK-06 GSR CLOSURE PASS CANDIDATE = action-consistency gate passed; raw GSR status remains reported above'
         write(unit, '(a)') '# TDVK-06 PASS CANDIDATE'
         close(unit)
         write(*, '(a)') 'TDVK-06 compact representation certification = PASS'
         write(*, '(a)') 'TDVK-06 ALSDA STATIC DIAGNOSTIC EXECUTED'
         write(*, '(a)') 'TDVK-06 GSR SOLVE EXECUTED status='//trim(gsr_result%status)
         write(*, '(a)') 'TDVK-06 GSR CLOSURE PASS CANDIDATE'
         write(*, '(a)') 'TDVK-06 PASS CANDIDATE'
      end if
      if (allocated(magnetization)) deallocate(magnetization)
      if (allocated(magnetization_valence)) deallocate(magnetization_valence)
      if (allocated(magnetization_core)) deallocate(magnetization_core)
      if (allocated(magnetization_compact)) deallocate(magnetization_compact)
      if (allocated(magnetization_point)) deallocate(magnetization_point)
      if (allocated(magnetization_projected_point)) deallocate(magnetization_projected_point)
      if (allocated(valence_compact)) deallocate(valence_compact)
      if (allocated(valence_point)) deallocate(valence_point)
      if (allocated(valence_projected_point)) deallocate(valence_projected_point)
      if (allocated(core_compact)) deallocate(core_compact)
      if (allocated(core_point)) deallocate(core_point)
      if (allocated(core_projected_point)) deallocate(core_projected_point)
      if (allocated(compact_kernel)) deallocate(compact_kernel)
      if (allocated(field)) deallocate(field)
      if (allocated(response)) deallocate(response)
      if (allocated(residual)) deallocate(residual)
      if (allocated(point_residual)) deallocate(point_residual)
   end subroutine run_tddft_static_interactions

   !> TDVK-02R2 validation-only handoff.  The accepted reciprocal snapshot is
   !> already naturally available here, so the live transition oracle and the
   !> compact bare response can be exercised without persisting eigenpairs or
   !> entering the full KXC/Dyson production lifecycle.
   subroutine run_tddft_product_bare_smoke(config, response_space, radial_bases, ground_states, left_state, endpoints, reciprocal_obj)
      type(tddft_production_config), intent(in) :: config
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
      if (product_plus%product_dimension /= 232 .or. product_minus%product_dimension /= 232) then
         error stop 'TDVK-02R2 product smoke: accepted Fe product dimension is not 232'
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
   subroutine run_tddft_product_finite_q(config, response_space, radial_bases, ground_states, left_state, endpoints)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(lmto_product_response_basis), target :: product_plus, product_minus
      type(lmto_product_response_basis), pointer :: product, opposite_product
      type(lr_product_ks_susceptibility_request) :: lehmann_request, covariance_request
      type(lr_product_ks_susceptibility_result) :: lehmann_result, covariance_result
      type(lr_product_gf_susceptibility_request) :: gf_request
      type(lr_product_gf_susceptibility_result) :: gf_result
      integer :: unit, iq, gamma_index, covariance_q_index, negative_q_index, arbitrary_q_index
      integer :: gf_q_index, gf_spot, n_gf_spots, ifrequency
      real(rp) :: accepted_moment, norm_lehmann, maximum_element, covariance_residual
      real(rp) :: difference_frobenius, relative_frobenius, difference_infinity
      real(rp) :: trace_real, trace_imag, q_folded(3), endpoint_first_k(3), endpoint_first(3), endpoint_last(3)
      real(rp) :: endpoint_error, q_tolerance
      integer :: endpoint_unique_count
      logical :: finite_response, covariance_checked
      character(len=32) :: channel, opposite_channel

      q_tolerance = 2.0e-11_rp
      gamma_index = find_gamma_q(config%q_list)
      covariance_q_index = find_first_nonzero_q(config%q_list)
      negative_q_index = find_matching_q(config%q_list, -config%q_list(:, covariance_q_index), q_tolerance)
      if (negative_q_index == 0) then
         error stop 'TDVK-04 finite-q validation: exact negative q is not present in q_list'
      end if
      arbitrary_q_index = find_arbitrary_q(config%q_list, gamma_index, covariance_q_index, negative_q_index)
      if (arbitrary_q_index == 0) then
         error stop 'TDVK-04 finite-q validation: arbitrary off-mesh q is not present in q_list'
      end if

      call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      call product_minus%initialize(response_space, radial_bases, lmto_product_channel_minus, .true.)
      if (product_plus%product_dimension /= 232 .or. product_minus%product_dimension /= 232) then
         error stop 'TDVK-04 finite-q validation: accepted Fe product dimension is not 232'
      end if
      if (trim(config%channel) == 'chi_plus') then
         product => product_plus
         opposite_product => product_minus
         channel = lr_channel_plus
         opposite_channel = lr_channel_minus
      else
         product => product_minus
         opposite_product => product_plus
         channel = lr_channel_minus
         opposite_channel = lr_channel_plus
      end if

      accepted_moment = 0.0_rp
      do iq = 1, size(ground_states)
         accepted_moment = accepted_moment + ground_states(iq)%integrated_moment_muB
      end do

      n_gf_spots = 2
      open(newunit=unit, file=trim(config%output_file), status='replace', action='write')
      write(unit, '(a)') '# TDVK-04 Fe finite-q compact bare-response validation'
      write(unit, '(a,a)') '# build_version = ', trim(tddft_build_version)
      write(unit, '(a)') '# backend = product_finite_q'
      write(unit, '(a)') '# no KXC, Goldstone, Dyson, loss, mode fitting, or point-space response allocation'
      write(unit, '(a,i0)') '# product_dimension = ', product%product_dimension
      write(unit, '(a,i0)') '# response_angular_cutoff = ', response_space%response_lmax
      write(unit, '(a,i0)') '# accepted_state_nbasis = ', left_state%nbasis
      write(unit, '(a,i0)') '# accepted_state_nbands = ', left_state%nbands
      write(unit, '(a,i0)') '# accepted_state_nk = ', left_state%nk
      write(unit, '(a,es24.16)') '# accepted_state_EF_Ry = ', left_state%fermi_level
      write(unit, '(a,es24.16)') '# accepted_state_temperature_K = ', left_state%temperature
      write(unit, '(a,es24.16)') '# accepted_state_moment_muB = ', accepted_moment
      write(unit, '(a,a)') '# accepted_state_reciprocal_mode = ', trim(left_state%reciprocal_mode)
      write(unit, '(a,a)') '# accepted_state_hamiltonian_order = ', trim(left_state%hamiltonian_order)
      write(unit, '(a,a)') '# accepted_state_provenance = ', &
         'same reciprocal eigenpairs, occupations, EF, temperature, radial mesh, product basis, channel, eta'
      write(unit, '(a,a)') '# q_convention = ', 'literal direct reciprocal coordinates; folded k endpoints use [-1/2,1/2)'
      write(unit, '(a,es24.16)') '# eta_Ry = ', config%eta
      write(unit, '(a,a)') '# q_metadata columns: q_index supplied_qx supplied_qy supplied_qz folded_qx folded_qy folded_qz endpoint_max_abs_error endpoint_first_kx endpoint_first_ky endpoint_first_kz endpoint_first_folded_kx endpoint_first_folded_ky endpoint_first_folded_kz endpoint_last_folded_kx endpoint_last_folded_ky endpoint_last_folded_kz endpoint_unique_count'
      write(unit, '(a,a)') '# lehmann_metrics columns: q_index supplied_qx supplied_qy supplied_qz frobenius_norm max_abs_element trace_real trace_imag transitions_evaluated finite'
      write(unit, '(a,a)') '# covariance columns: q_index negative_q_index omega_Ry residual channel negative_channel'
      write(unit, '(a,a)') '# covariance convention: chi(q,w) compared with C*conjg(chi(-q,-w))*C; C includes (L,M)->(L,-M), (-1)^M, and compact radial-basis transport'
      write(unit, '(a,a)') '# gf_metrics columns: q_index supplied_qx supplied_qy supplied_qz integration_points integration_eta h_over_integration_eta norm_lehmann norm_gf dF rF dInf wall_seconds finite'

      covariance_checked = .false.
      do iq = 1, size(config%q_list, 2)
         call finite_q_endpoint_metadata(left_state, endpoints(iq), config%q_list(:, iq), q_folded, endpoint_error, &
            endpoint_first_k, endpoint_first, endpoint_last, endpoint_unique_count)
         if (endpoint_error > q_tolerance) then
            error stop 'TDVK-04 finite-q validation: endpoint metadata failed exact folding check'
         end if

         lehmann_request%q = config%q_list(:, iq)
         lehmann_request%frequencies = config%frequencies
         lehmann_request%eta = config%eta
         lehmann_request%channel = channel
         lehmann_request%product_basis => product
         lehmann_request%electronic_state => left_state
         lehmann_request%q_endpoint_state => endpoints(iq)
         call evaluate_lr_product_ks_susceptibility(lehmann_request, lehmann_result)
         finite_response = all(ieee_is_finite(real(lehmann_result%susceptibility, rp))) .and. &
            all(ieee_is_finite(aimag(lehmann_result%susceptibility)))
         if (.not. finite_response) then
            error stop 'TDVK-04 finite-q validation: compact Lehmann response contains NaN or Inf'
         end if
         norm_lehmann = sqrt(sum(abs(lehmann_result%susceptibility(:, :, 1))**2))
         maximum_element = maxval(abs(lehmann_result%susceptibility(:, :, 1)))
         trace_real = 0.0_rp
         trace_imag = 0.0_rp
         do ifrequency = 1, product%product_dimension
            trace_real = trace_real + real(lehmann_result%susceptibility(ifrequency, ifrequency, 1), rp)
            trace_imag = trace_imag + aimag(lehmann_result%susceptibility(ifrequency, ifrequency, 1))
         end do
         write(unit, '(i0,1x,4(es24.16,1x),3(es24.16,1x),i0,1x,l1)') iq, config%q_list(:, iq), norm_lehmann, &
            maximum_element, trace_real, trace_imag, lehmann_result%ntransitions_evaluated, finite_response
         write(unit, '(i0,1x,3(es24.16,1x),3(es24.16,1x),es24.16,1x,9(es24.16,1x),i0)') iq, config%q_list(:, iq), &
            q_folded, endpoint_error, endpoint_first_k, endpoint_first, endpoint_last, endpoint_unique_count

         if (iq == covariance_q_index) then
            covariance_request%q = -config%q_list(:, iq)
            covariance_request%frequencies = -config%frequencies
            covariance_request%eta = config%eta
            covariance_request%channel = opposite_channel
            covariance_request%product_basis => opposite_product
            covariance_request%electronic_state => left_state
            covariance_request%q_endpoint_state => endpoints(negative_q_index)
            call evaluate_lr_product_ks_susceptibility(covariance_request, covariance_result)
            call compact_covariance_residual(product, opposite_product, lehmann_result, covariance_result, covariance_residual)
            if (.not. ieee_is_finite(covariance_residual)) then
               error stop 'TDVK-04 finite-q validation: covariance residual is not finite'
            end if
            write(unit, '(i0,1x,i0,1x,es24.16,1x,es24.16,1x,a,1x,a)') iq, negative_q_index, config%frequencies(1), &
               covariance_residual, trim(channel), trim(opposite_channel)
            write(*, '(a,3(es16.8,1x),a,3(es16.8,1x),a,es12.4)') 'TDVK-04 covariance q=', config%q_list(:, iq), &
               ' -q=', config%q_list(:, negative_q_index), ' max_abs_residual=', covariance_residual
            covariance_checked = .true.
         end if
      end do
      if (.not. covariance_checked) error stop 'TDVK-04 finite-q validation: covariance pair was not evaluated'

      do gf_spot = 1, n_gf_spots
         if (gf_spot == 1) then
            gf_q_index = covariance_q_index
         else
            gf_q_index = arbitrary_q_index
         end if
         gf_request%q = config%q_list(:, gf_q_index)
         gf_request%frequencies = config%frequencies
         gf_request%eta = config%eta
         gf_request%channel = channel
         gf_request%integration_points = config%gf_integration_points
         gf_request%integration_eta = config%gf_integration_eta
         gf_request%energy_margin = config%gf_energy_margin
         gf_request%contraction_backend = 'factorized'
         gf_request%product_basis => product
         gf_request%electronic_state => left_state
         gf_request%q_endpoint_state => endpoints(gf_q_index)
         call evaluate_lr_product_gf_susceptibility(gf_request, gf_result)
         finite_response = all(ieee_is_finite(real(gf_result%susceptibility, rp))) .and. &
            all(ieee_is_finite(aimag(gf_result%susceptibility)))
         if (.not. finite_response) error stop 'TDVK-04 finite-q validation: compact GF response contains NaN or Inf'
         lehmann_request%q = config%q_list(:, gf_q_index)
         lehmann_request%frequencies = config%frequencies
         lehmann_request%eta = config%eta
         lehmann_request%channel = channel
         lehmann_request%product_basis => product
         lehmann_request%electronic_state => left_state
         lehmann_request%q_endpoint_state => endpoints(gf_q_index)
         call evaluate_lr_product_ks_susceptibility(lehmann_request, lehmann_result)
         norm_lehmann = sqrt(sum(abs(lehmann_result%susceptibility(:, :, 1))**2))
         difference_frobenius = sqrt(sum(abs(lehmann_result%susceptibility(:, :, 1) - &
            gf_result%susceptibility(:, :, 1))**2))
         relative_frobenius = difference_frobenius/max(norm_lehmann, &
            sqrt(sum(abs(gf_result%susceptibility(:, :, 1))**2)), tiny(1.0_rp))
         difference_infinity = maxval(abs(lehmann_result%susceptibility(:, :, 1) - &
            gf_result%susceptibility(:, :, 1)))
         call finite_q_endpoint_metadata(left_state, endpoints(gf_q_index), config%q_list(:, gf_q_index), q_folded, &
            endpoint_error, endpoint_first_k, endpoint_first, endpoint_last, endpoint_unique_count)
         write(unit, '(i0,1x,3(es24.16,1x),i0,1x,8(es24.16,1x),l1)') gf_q_index, config%q_list(:, gf_q_index), &
            gf_result%integration_points, gf_result%actual_integration_eta, gf_result%spacing_over_integration_eta, &
            norm_lehmann, sqrt(sum(abs(gf_result%susceptibility(:, :, 1))**2)), difference_frobenius, &
            relative_frobenius, difference_infinity, gf_result%wall_time_seconds, finite_response
         write(*, '(a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') 'TDVK-04 GF q index=', gf_q_index, &
            ' dF=', difference_frobenius, ' rF=', relative_frobenius, ' dInf=', difference_infinity, &
            ' wall_s=', gf_result%wall_time_seconds
      end do
      close(unit)
      write(*, '(a,i0,a,i0,a,i0,a,i0)') 'TDVK-04 Fe finite-q compact validation: q_count=', size(config%q_list, 2), &
         ' product_dimension=', product%product_dimension, ' covariance_q_index=', covariance_q_index, &
         ' gf_spots=', n_gf_spots
   end subroutine run_tddft_product_finite_q

   !> TDVK-05 compact Lehmann convergence campaign.  The caller supplies one
   !> accepted SCF state and its exact q endpoints; this routine varies only
   !> the physical response eta values and records complete compact-matrix
   !> diagnostics.  It deliberately stops before reciprocal GF, KXC,
   !> Goldstone, Dyson, loss, and mode interpretation.
   subroutine run_tddft_product_convergence(config, response_space, radial_bases, ground_states, left_state, endpoints, &
                                            reciprocal_obj, accepted_kspace_scf)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      logical, intent(in), optional :: accepted_kspace_scf
      type(lmto_product_response_basis), target :: product_plus, product_minus
      type(lmto_product_response_basis), pointer :: product
      type(lr_product_ks_susceptibility_request) :: request
      type(lr_product_ks_susceptibility_result) :: result
      integer :: unit, kpoint_unit, basis_unit, matrix_unit, gf_unit
      integer :: iq, ieta, ifrequency, isite, response_l, ik, imode, ir, matrix_i, matrix_j
      integer :: gamma_index, static_index
      real(rp) :: accepted_moment, runtime_start, runtime_end, fixed_ef_electrons
      real(rp) :: diagnostic_ef, eigen_min, eigen_max, eigen_mean, weight_sum
      real(rp) :: frobenius_norm, maximum_element, trace_real, trace_imag
      real(rp) :: k_fingerprint(5)
      real(rp) :: state_mesh_max, state_weight_max, state_ef_diff, state_eigen_max, state_occ_max
      real(rp) :: state_projector_max, state_projector_frobenius
      real(rp) :: gf_norm_lehmann, gf_norm, gf_d_frobenius, gf_relative_frobenius, gf_d_infinity
      complex(rp) :: trace
      logical :: finite_response, rank_stable, matrix_written, direct_handoff
      character(len=512) :: kpoint_file, basis_file, matrix_file, state_file, gf_file
      type(lr_product_ks_susceptibility_request) :: gf_lehmann_request
      type(lr_product_ks_susceptibility_result) :: gf_lehmann_result
      type(lr_product_gf_susceptibility_request) :: gf_request
      type(lr_product_gf_susceptibility_result) :: gf_result

      if (size(endpoints) /= size(config%q_list, 2)) then
         error stop 'TDVK-05 convergence: q endpoint count differs from q_list'
      end if
      direct_handoff = .false.
      if (present(accepted_kspace_scf)) direct_handoff = accepted_kspace_scf
      call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      call product_minus%initialize(response_space, radial_bases, lmto_product_channel_minus, .true.)
      if (trim(config%channel) == 'chi_plus') then
         product => product_plus
      else
         product => product_minus
      end if
      if (product%product_dimension < 1) error stop 'TDVK-05 convergence: compact product space is empty'
      if (reciprocal_obj%k_workset%nk_global /= left_state%nk .or. &
          reciprocal_obj%k_workset%nk_local /= left_state%nk) then
         error stop 'TDVK-05 convergence: reciprocal workset and immutable state have different k counts'
      end if
      if (.not. reciprocal_obj%k_workset%complete_bz .or. reciprocal_obj%k_workset%distributed) then
         error stop 'TDVK-05 convergence: a complete replicated BZ workset is required'
      end if
      gamma_index = 0
      do iq = 1, size(config%q_list, 2)
         if (sum(abs(config%q_list(:, iq))) <= 1.0e-12_rp) gamma_index = iq
      end do
      static_index = 0
      do ifrequency = 1, size(config%frequencies)
         if (abs(config%frequencies(ifrequency)) <= 1.0e-12_rp) static_index = ifrequency
      end do
      if (gamma_index == 0 .or. static_index == 0) then
         error stop 'TDVK-05 convergence: Gamma/static indices were not found'
      end if

      accepted_moment = 0.0_rp
      do isite = 1, size(ground_states)
         accepted_moment = accepted_moment + ground_states(isite)%integrated_moment_muB
      end do
      rank_stable = .true.
      do response_l = 0, product%response_lmax
         do isite = 1, product%nsite
            rank_stable = rank_stable .and. product%blocks(isite, response_l)%rank_stable
         end do
      end do

      fixed_ef_electrons = reciprocal_obj%canonical_electron_count
      weight_sum = sum(left_state%k_weights)
      eigen_min = minval(left_state%eigenvalues)
      eigen_max = maxval(left_state%eigenvalues)
      eigen_mean = sum(left_state%eigenvalues)/real(size(left_state%eigenvalues), rp)
      diagnostic_ef = reciprocal_obj%fermi_level
      if (reciprocal_obj%total_electrons > 0.0_rp) then
         diagnostic_ef = reciprocal_obj%find_fermi_level_from_eigenvalues(reciprocal_obj%total_electrons)
      end if
      k_fingerprint = 0.0_rp
      do ik = 1, left_state%nk
         k_fingerprint(1) = k_fingerprint(1) + left_state%k_weights(ik)
         k_fingerprint(2:4) = k_fingerprint(2:4) + left_state%k_weights(ik)*left_state%k_points(:, ik)
         k_fingerprint(5) = k_fingerprint(5) + real(ik, rp)*left_state%k_weights(ik)
      end do
      kpoint_file = trim(config%output_file)//'.kpoints'
      basis_file = trim(config%output_file)//'.basis'
      matrix_file = trim(config%output_file)//'.matrix'
      state_file = trim(config%output_file)//'.state'
      gf_file = trim(config%output_file)//'.gf'
      call state_consistency_metrics(reciprocal_obj, left_state, state_mesh_max, state_weight_max, state_ef_diff, &
         state_eigen_max, state_occ_max, state_projector_max, state_projector_frobenius)

      open(newunit=unit, file=trim(config%output_file), status='replace', action='write')
      open(newunit=kpoint_unit, file=trim(kpoint_file), status='replace', action='write')
      open(newunit=basis_unit, file=trim(basis_file), status='replace', action='write')
      open(newunit=matrix_unit, file=trim(matrix_file), status='replace', action='write')
      write(unit, '(a)') '# TDVK-05 Fe compact bare-response numerical convergence'
      write(unit, '(a,a)') '# build_version = ', trim(tddft_build_version)
      write(unit, '(a)') '# backend = product_convergence'
      write(unit, '(a)') '# no reciprocal-GF ladder, KXC, Goldstone, Dyson, loss, mode fitting, or point-space response allocation'
      write(unit, '(a,i0)') '# accepted_state_nbasis = ', left_state%nbasis
      write(unit, '(a,i0)') '# accepted_state_nbands = ', left_state%nbands
      write(unit, '(a,i0)') '# accepted_state_nk = ', left_state%nk
      write(unit, '(a,es24.16)') '# accepted_state_EF_Ry = ', left_state%fermi_level
      write(unit, '(a,es24.16)') '# accepted_state_temperature_K = ', left_state%temperature
      write(unit, '(a,es24.16)') '# accepted_state_moment_muB = ', accepted_moment
      write(unit, '(a,es24.16)') '# accepted_state_radial_residual_control = ', ground_states(1)%accepted_residual_control
      write(unit, '(a,a)') '# accepted_state_reciprocal_mode = ', trim(left_state%reciprocal_mode)
      write(unit, '(a,a)') '# accepted_state_hamiltonian_order = ', trim(left_state%hamiltonian_order)
      if (direct_handoff) then
         write(unit, '(a)') '# accepted_state_source = accepted_kspace_scf_cache'
         write(unit, '(a)') '# accepted_state_provenance = self-consistent k-space SCF cache handed directly to TDDFT; no reciprocal rebuild'
      else
         write(unit, '(a)') '# accepted_state_source = diagnostic_frozen_post_scf_rebuild'
         write(unit, '(a)') '# accepted_state_provenance = accepted real-space SCF potential rebuilt for diagnostic reciprocal response only'
      end if
      write(unit, '(a,l1)') '# reciprocal_rebuild_performed_for_tddft = ', .not. direct_handoff
      write(unit, '(a,es24.16)') '# state_consistency_mesh_max_abs = ', state_mesh_max
      write(unit, '(a,es24.16)') '# state_consistency_weight_max_abs = ', state_weight_max
      write(unit, '(a,es24.16)') '# state_consistency_EF_max_abs_Ry = ', state_ef_diff
      write(unit, '(a,es24.16)') '# state_consistency_eigenvalue_max_abs_Ry = ', state_eigen_max
      write(unit, '(a,es24.16)') '# state_consistency_occupation_max_abs = ', state_occ_max
      write(unit, '(a,es24.16)') '# state_consistency_projector_max_abs = ', state_projector_max
      write(unit, '(a,es24.16)') '# state_consistency_projector_frobenius = ', state_projector_frobenius
      write(unit, '(a,3(i0,1x))') '# requested_k_mesh = ', reciprocal_obj%nk_mesh
      write(unit, '(a,3(i0,1x))') '# actual_generated_k_mesh = ', reciprocal_obj%nk_mesh
      write(unit, '(a,i0)') '# actual_k_count = ', reciprocal_obj%k_workset%nk_global
      write(unit, '(a,es24.16)') '# actual_k_weight_sum = ', weight_sum
      write(unit, '(a,3(es24.16,1x))') '# representative_k_first = ', left_state%k_points(:, 1)
      write(unit, '(a,3(es24.16,1x))') '# representative_k_middle = ', left_state%k_points(:, (left_state%nk + 1)/2)
      write(unit, '(a,3(es24.16,1x))') '# representative_k_last = ', left_state%k_points(:, left_state%nk)
      write(unit, '(a,5(es24.16,1x))') '# k_fingerprint_checksums = ', k_fingerprint
      write(unit, '(a,es24.16)') '# eigenvalue_min_Ry = ', eigen_min
      write(unit, '(a,es24.16)') '# eigenvalue_max_Ry = ', eigen_max
      write(unit, '(a,es24.16)') '# eigenvalue_mean_Ry = ', eigen_mean
      write(unit, '(a,es24.16)') '# target_electron_count = ', reciprocal_obj%total_electrons
      if (direct_handoff) then
         write(unit, '(a,es24.16)') '# accepted_state_integrated_electron_count = ', fixed_ef_electrons
         write(unit, '(a,es24.16)') '# accepted_state_integrated_electron_count_error = ', &
            fixed_ef_electrons - reciprocal_obj%total_electrons
         write(unit, '(a)') '# accepted_state_fermi_owner = reciprocal electron-number occupation solver'
      else
         write(unit, '(a,es24.16)') '# fixed_EF_electron_count = ', fixed_ef_electrons
         write(unit, '(a,es24.16)') '# fixed_EF_electron_count_error = ', fixed_ef_electrons - reciprocal_obj%total_electrons
         write(unit, '(a,es24.16)') '# diagnostic_mesh_EF_Ry = ', diagnostic_ef
         write(unit, '(a,es24.16)') '# diagnostic_mesh_EF_shift_Ry = ', diagnostic_ef - left_state%fermi_level
         write(unit, '(a)') '# accepted_state_fermi_owner = input EF retained for diagnostic frozen-potential response'
      end if
      write(unit, '(a,a)') '# kpoint_artifact = ', trim(kpoint_file)
      write(unit, '(a,a)') '# basis_artifact = ', trim(basis_file)
      write(unit, '(a,a)') '# matrix_artifact = ', trim(matrix_file)
      write(unit, '(a,a)') '# state_artifact = ', trim(state_file)
      if (config%gf_closure_audit) write(unit, '(a,a)') '# gf_spot_artifact = ', trim(gf_file)
      write(unit, '(a,a)') '# channel = ', trim(config%channel)
      write(unit, '(a,es24.16)') '# primary_eta_Ry = ', config%eta
      write(unit, '(a,i0)') '# response_lmax = ', response_space%response_lmax
      if (response_space%response_lmax == 4) then
         write(unit, '(a)') '# response_space_label = complete certified spd product span'
      else
         write(unit, '(a)') '# response_space_label = reduced approximate response cutoff'
      end if
      write(unit, '(a,i0)') '# product_unpruned_dimension = ', product%unpruned_dimension
      write(unit, '(a,i0)') '# product_dimension = ', product%product_dimension
      write(unit, '(a,l1)') '# product_rank_stable_all_L = ', rank_stable
      write(unit, '(a,i0)') '# product_basis_mode_count = ', product%product_dimension
      write(unit, '(a,i0,4(es24.16,1x))') '# radial_mesh_identity = ', size(ground_states(1)%r), ground_states(1)%a, &
         ground_states(1)%b, ground_states(1)%rmax, sum(ground_states(1)%r)
      write(unit, '(a,a)') '# radial_mesh_provenance = ', &
         'accepted direct LR-01 logarithmic mesh; no decimation or material-dependent radial truncation'
      write(unit, '(a,a)') '# xc_provenance = ', trim(ground_states(1)%xc_provenance%functional_name)//' / '// &
         trim(ground_states(1)%xc_provenance%backend_name)
      write(unit, '(a,i0)') '# eta_count = ', size(config%eta_values)
      write(unit, '(a,*(es24.16,1x))') '# eta_values_Ry = ', config%eta_values
      write(unit, '(a,a)') '# columns: eta_index eta_Ry q_index qx qy qz omega_Ry frobenius_norm max_abs_element trace_real trace_imag transitions_evaluated occupation_skips runtime_cpu_seconds finite'

      write(kpoint_unit, '(a,3(i0,1x))') '# requested_k_mesh = ', reciprocal_obj%nk_mesh
      write(kpoint_unit, '(a,i0)') '# actual_k_count = ', reciprocal_obj%k_workset%nk_global
      write(kpoint_unit, '(a,es24.16)') '# actual_k_weight_sum = ', weight_sum
      write(kpoint_unit, '(a,5(es24.16,1x))') '# k_fingerprint_checksums = ', k_fingerprint
      write(kpoint_unit, '(a)') '# columns: k_index kx ky kz weight'
      do ik = 1, left_state%nk
         write(kpoint_unit, '(i0,1x,4(es24.16,1x))') ik, left_state%k_points(:, ik), left_state%k_weights(ik)
      end do

      write(basis_unit, '(a,3(i0,1x))') '# requested_k_mesh = ', reciprocal_obj%nk_mesh
      write(basis_unit, '(a,i0)') '# basis_product_dimension = ', product%product_dimension
      write(basis_unit, '(a,i0)') '# basis_unpruned_dimension = ', product%unpruned_dimension
      write(basis_unit, '(a,l1)') '# basis_rank_stable_all_L = ', rank_stable
      write(basis_unit, '(a,i0)') '# basis_response_lmax = ', product%response_lmax
      write(basis_unit, '(a)') '# columns: site response_l radial_index mode real imag singular_value'
      do isite = 1, product%nsite
         do response_l = 0, product%response_lmax
            do imode = 1, product%blocks(isite, response_l)%rank
               do ir = 1, product%blocks(isite, response_l)%npoint
                  write(basis_unit, '(3(i0,1x),i0,1x,3(es24.16,1x))') isite, response_l, ir, imode, &
                     real(product%blocks(isite, response_l)%weighted_modes(ir, imode), rp), &
                     aimag(product%blocks(isite, response_l)%weighted_modes(ir, imode)), &
                     product%blocks(isite, response_l)%singular_values(imode)
               end do
            end do
         end do
      end do

      write(matrix_unit, '(a,3(i0,1x))') '# requested_k_mesh = ', reciprocal_obj%nk_mesh
      write(matrix_unit, '(a,i0)') '# basis_product_dimension = ', product%product_dimension
      write(matrix_unit, '(a,i0)') '# basis_response_lmax = ', product%response_lmax
      write(matrix_unit, '(a,l1)') '# basis_rank_stable_all_L = ', rank_stable
      write(matrix_unit, '(a)') '# matrix_scope = Gamma, omega=0, every requested physical eta'
      write(matrix_unit, '(a)') '# columns: eta_index eta_Ry omega_Ry row column real imag'
      matrix_written = .false.

      do ieta = 1, size(config%eta_values)
         do iq = 1, size(config%q_list, 2)
            request%q = config%q_list(:, iq)
            request%frequencies = config%frequencies
            request%eta = config%eta_values(ieta)
            request%channel = config%channel
            request%product_basis => product
            request%electronic_state => left_state
            request%q_endpoint_state => endpoints(iq)
            call cpu_time(runtime_start)
            call evaluate_lr_product_ks_susceptibility(request, result)
            call cpu_time(runtime_end)
            finite_response = all(ieee_is_finite(real(result%susceptibility, rp))) .and. &
               all(ieee_is_finite(aimag(result%susceptibility)))
            if (.not. finite_response) then
               error stop 'TDVK-05 convergence: compact Lehmann response contains NaN or Inf'
            end if
            do ifrequency = 1, size(result%frequencies)
               frobenius_norm = sqrt(sum(abs(result%susceptibility(:, :, ifrequency))**2))
               maximum_element = maxval(abs(result%susceptibility(:, :, ifrequency)))
               trace = cmplx(0.0_rp, 0.0_rp, rp)
               do isite = 1, result%product_dimension
                  trace = trace + result%susceptibility(isite, isite, ifrequency)
               end do
               trace_real = real(trace, rp)
               trace_imag = aimag(trace)
               write(unit, '(i0,1x,es24.16,1x,i0,1x,3(es24.16,1x),es24.16,1x,4(es24.16,1x),2(i0,1x),es24.16,1x,l1)') &
                  ieta, config%eta_values(ieta), iq, config%q_list(:, iq), result%frequencies(ifrequency), frobenius_norm, &
                  maximum_element, trace_real, trace_imag, result%ntransitions_evaluated, result%noccupation_skips, &
                  runtime_end - runtime_start, finite_response
               if (iq == gamma_index .and. ifrequency == static_index) then
                  matrix_written = .true.
                  do matrix_i = 1, product%product_dimension
                     do matrix_j = 1, product%product_dimension
                        write(matrix_unit, '(i0,1x,es24.16,1x,es24.16,1x,2(i0,1x),2(es24.16,1x))') &
                           ieta, config%eta_values(ieta), result%frequencies(ifrequency), matrix_i, matrix_j, &
                           real(result%susceptibility(matrix_i, matrix_j, ifrequency), rp), &
                           aimag(result%susceptibility(matrix_i, matrix_j, ifrequency))
                     end do
                  end do
               end if
            end do
         end do
      end do
      close(unit)
      close(kpoint_unit)
      close(basis_unit)
      if (.not. matrix_written) error stop 'TDVK-05 convergence: full Gamma/static matrix was not written'
      close(matrix_unit)
      if (config%gf_closure_audit) then
         ! One independent Gamma GF spot is evaluated after the Lehmann
         ! response, in this same process and from the same immutable state.
         gf_request%q = config%q_list(:, gamma_index)
         gf_request%frequencies = config%frequencies
         gf_request%eta = config%eta
         gf_request%channel = config%channel
         gf_request%integration_points = config%gf_integration_points
         gf_request%integration_eta = config%gf_integration_eta
         gf_request%energy_margin = config%gf_energy_margin
         gf_request%contraction_backend = 'factorized'
         gf_request%product_basis => product
         gf_request%electronic_state => left_state
         gf_request%q_endpoint_state => endpoints(gamma_index)
         call evaluate_lr_product_gf_susceptibility(gf_request, gf_result)
         finite_response = all(ieee_is_finite(real(gf_result%susceptibility, rp))) .and. &
            all(ieee_is_finite(aimag(gf_result%susceptibility)))
         if (.not. finite_response) error stop 'TDVK k-space handoff: Gamma GF spot contains NaN or Inf'

         gf_lehmann_request%q = config%q_list(:, gamma_index)
         gf_lehmann_request%frequencies = config%frequencies
         gf_lehmann_request%eta = config%eta
         gf_lehmann_request%channel = config%channel
         gf_lehmann_request%product_basis => product
         gf_lehmann_request%electronic_state => left_state
         gf_lehmann_request%q_endpoint_state => endpoints(gamma_index)
         call evaluate_lr_product_ks_susceptibility(gf_lehmann_request, gf_lehmann_result)
         gf_norm_lehmann = sqrt(sum(abs(gf_lehmann_result%susceptibility(:, :, static_index))**2))
         gf_norm = sqrt(sum(abs(gf_result%susceptibility(:, :, static_index))**2))
         gf_d_frobenius = sqrt(sum(abs(gf_lehmann_result%susceptibility(:, :, static_index) - &
            gf_result%susceptibility(:, :, static_index))**2))
         gf_relative_frobenius = gf_d_frobenius/max(gf_norm_lehmann, gf_norm, tiny(1.0_rp))
         gf_d_infinity = maxval(abs(gf_lehmann_result%susceptibility(:, :, static_index) - &
            gf_result%susceptibility(:, :, static_index)))
         open(newunit=gf_unit, file=trim(gf_file), status='replace', action='write')
         write(gf_unit, '(a)') '# one reciprocal-GF spot check on the accepted TDDFT state'
         write(gf_unit, '(a)') '# state_source = same accepted reciprocal state as the compact Lehmann response'
         write(gf_unit, '(a,3(es24.16,1x))') '# q = ', config%q_list(:, gamma_index)
         write(gf_unit, '(a,es24.16)') '# omega_Ry = ', config%frequencies(static_index)
         write(gf_unit, '(a,es24.16)') '# eta_Ry = ', config%eta
         write(gf_unit, '(a,i0)') '# integration_points = ', gf_result%integration_points
         write(gf_unit, '(a,es24.16)') '# integration_eta_Ry = ', gf_result%actual_integration_eta
         write(gf_unit, '(a,es24.16)') '# h_over_integration_eta = ', gf_result%spacing_over_integration_eta
         write(gf_unit, '(a)') '# columns: norm_lehmann norm_gf dF rF dInf finite'
         write(gf_unit, '(5(es24.16,1x),l1)') gf_norm_lehmann, gf_norm, gf_d_frobenius, gf_relative_frobenius, &
            gf_d_infinity, finite_response
         close(gf_unit)
      end if
      write(*, '(a,i0,a,i0,a,i0,a,i0)') 'TDVK-05 Fe compact convergence: q_count=', size(config%q_list, 2), &
         ' omega_count=', size(config%frequencies), ' eta_count=', size(config%eta_values), &
         ' product_dimension=', product%product_dimension
   end subroutine run_tddft_product_convergence

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
   subroutine project_compact_covariance_to_sites(contract, plus_product, minus_product, minus_matrix, site_matrix)
      type(projected_site_spin_contract), intent(in) :: contract
      type(lmto_product_response_basis), intent(in) :: plus_product, minus_product
      complex(rp), intent(in) :: minus_matrix(:, :)
      complex(rp), intent(out) :: site_matrix(:, :)
      complex(rp), allocatable :: transport(:, :), mapped(:, :), functionals(:, :), image(:), block_transport(:, :)
      integer :: site_index, response_l, response_m, plus_first, minus_first, rank_plus, rank_minus, i, j
      real(rp) :: angular_sign

      if (plus_product%product_dimension /= minus_product%product_dimension) then
         error stop 'DRESP-06A site covariance: compact dimensions differ'
      end if
      if (any(shape(minus_matrix) /= [minus_product%product_dimension, minus_product%product_dimension])) then
         error stop 'DRESP-06A site covariance: compact matrix shape mismatch'
      end if
      if (any(shape(site_matrix) /= [contract%nsite, contract%nsite])) then
         error stop 'DRESP-06A site covariance: site matrix shape mismatch'
      end if
      allocate(transport(plus_product%product_dimension, minus_product%product_dimension), &
         mapped(plus_product%product_dimension, plus_product%product_dimension), &
         functionals(plus_product%product_dimension, contract%nsite), image(plus_product%product_dimension))
      transport = cmplx(0.0_rp, 0.0_rp, rp)
      do site_index = 1, plus_product%nsite
         do response_l = 0, plus_product%response_lmax
            rank_plus = plus_product%blocks(site_index, response_l)%rank
            rank_minus = minus_product%blocks(site_index, response_l)%rank
            allocate(block_transport(rank_plus, rank_minus))
            block_transport = matmul(conjg(transpose(plus_product%blocks(site_index, response_l)%weighted_modes)), &
               conjg(minus_product%blocks(site_index, response_l)%weighted_modes))
            do response_m = -response_l, response_l
               plus_first = plus_product%flat_index(site_index, response_l, response_m, 1)
               minus_first = minus_product%flat_index(site_index, response_l, -response_m, 1)
               angular_sign = merge(-1.0_rp, 1.0_rp, mod(abs(response_m), 2) == 1)
               transport(plus_first:plus_first + rank_plus - 1, minus_first:minus_first + rank_minus - 1) = &
                  angular_sign*block_transport
            end do
            deallocate(block_transport)
         end do
      end do
      mapped = matmul(transport, matmul(conjg(minus_matrix), conjg(transpose(transport))))
      call contract%site_integration_functional(plus_product, functionals)
      do j = 1, contract%nsite
         image = matmul(mapped, functionals(:, j))
         do i = 1, contract%nsite
            site_matrix(i, j) = dot_product(functionals(:, i), image)
         end do
      end do
      deallocate(transport, mapped, functionals, image)
   end subroutine project_compact_covariance_to_sites

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
   subroutine run_tddft_product_gf_smoke(config, response_space, radial_bases, left_state, endpoints)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(lmto_product_response_basis), target :: product_plus, product_minus
      type(lr_product_gf_susceptibility_request) :: request
      type(lr_product_gf_susceptibility_result) :: result
      integer :: gamma_index, channel_kind, i
      real(rp) :: frobenius, maximum_element
      complex(rp) :: trace
      logical :: finite_response

      gamma_index = find_gamma_q(config%q_list)
      if (gamma_index > size(endpoints)) error stop 'TDVK-02R3 product GF smoke: Gamma endpoint is unavailable'
      call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      call product_minus%initialize(response_space, radial_bases, lmto_product_channel_minus, .true.)
      if (product_plus%product_dimension /= 232 .or. product_minus%product_dimension /= 232) then
         error stop 'TDVK-02R3 product GF smoke: accepted Fe product dimension is not 232'
      end if

      request%q = config%q_list(:, gamma_index)
      request%frequencies = [0.0_rp]
      request%eta = config%eta
      request%channel = config%channel
      request%integration_points = config%gf_integration_points
      request%integration_eta = config%gf_integration_eta
      request%energy_margin = config%gf_energy_margin
      request%electronic_state => left_state
      request%q_endpoint_state => endpoints(gamma_index)
      if (trim(config%channel) == 'chi_plus') then
         request%product_basis => product_plus
      else
         request%product_basis => product_minus
      end if
      call evaluate_lr_product_gf_susceptibility(request, result)
      finite_response = all(ieee_is_finite(real(result%susceptibility, rp))) .and. &
         all(ieee_is_finite(aimag(result%susceptibility)))
      if (.not. finite_response) error stop 'TDVK-02R3 product GF smoke: compact response contains NaN or Inf'

      frobenius = sqrt(sum(abs(result%susceptibility(:, :, 1))**2))
      maximum_element = maxval(abs(result%susceptibility(:, :, 1)))
      trace = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, result%product_dimension
         trace = trace + result%susceptibility(i, i, 1)
      end do
      write (*, '(a)') 'TDVK-02R3 compact reciprocal-GF Fe smoke: execution/performance evidence only — not quadrature convergence'
      write (*, '(a,l1,a,i0,a,i0,a,es12.4)') '  finite=', finite_response, ' Nprod=', result%product_dimension, &
         ' integration_points=', result%integration_points, ' integration_eta=', result%actual_integration_eta
      write (*, '(a,es12.4,a,es12.4,a,2(es12.4,1x))') '  Frobenius_norm=', frobenius, &
         ' max_element=', maximum_element, ' trace=', real(trace, rp), aimag(trace)
      write (*, '(a,es12.4,a,es12.4,a,es12.4,a,i0,a,i0,a,i0,a,i0)') '  wall_s=', result%wall_time_seconds, &
         ' cpu_s=', result%cpu_time_seconds, ' time_per_energy_kpoint_s=', result%time_per_energy_kpoint, &
         ' energy_points=', result%integration_points, &
         ' k_points=', result%nk, ' frequencies=', result%nfrequency, &
         ' component_vertex_bytes=', result%component_vertex_memory_bytes
      write (*, '(a,i0,a,i0,a,i0)') '  gf_matrix_bytes=', result%gf_matrix_memory_bytes, &
         ' susceptibility_bytes=', result%susceptibility_memory_bytes, ' product_dimension=', result%product_dimension
   end subroutine run_tddft_product_gf_smoke

   !> TDVK-03 accepted-state closure audit.  The SCF handoff above has already
   !> prepared one immutable reciprocal state and all exact folded endpoints.
   !> This routine deliberately evaluates the compact Lehmann reference once,
   !> then varies only GF integration controls while keeping those snapshots
   !> and the complete 232-coordinate product basis fixed.
   subroutine run_tddft_product_gf_closure(config, response_space, radial_bases, left_state, endpoints)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(lmto_product_response_basis), target :: product_plus, product_minus
      type(lr_product_ks_susceptibility_request) :: lehmann_request
      type(lr_product_ks_susceptibility_result) :: lehmann_result
      real(rp), parameter :: eta_scale(4) = [10.0_rp, 5.0_rp, 2.5_rp, 1.0_rp]
      real(rp), parameter :: margin_increment(3) = [0.0_rp, 0.4_rp, 1.4_rp]
      real(rp) :: integration_eta, base_integration_eta, margin
      integer :: gamma_index, i, integration_points
      logical :: finite_response

      gamma_index = find_gamma_q(config%q_list)
      if (gamma_index > size(endpoints)) error stop 'TDVK-03 closure: Gamma endpoint is unavailable'
      call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      call product_minus%initialize(response_space, radial_bases, lmto_product_channel_minus, .true.)
      if (product_plus%product_dimension /= 232 .or. product_minus%product_dimension /= 232) then
         error stop 'TDVK-03 closure: accepted Fe product dimension is not 232'
      end if

      lehmann_request%q = config%q_list(:, gamma_index)
      lehmann_request%frequencies = [0.0_rp]
      lehmann_request%eta = config%eta
      lehmann_request%channel = config%channel
      lehmann_request%electronic_state => left_state
      lehmann_request%q_endpoint_state => endpoints(gamma_index)
      if (trim(config%channel) == 'chi_plus') then
         lehmann_request%product_basis => product_plus
      else
         lehmann_request%product_basis => product_minus
      end if
      call evaluate_lr_product_ks_susceptibility(lehmann_request, lehmann_result)
      finite_response = all(ieee_is_finite(real(lehmann_result%susceptibility, rp))) .and. &
         all(ieee_is_finite(aimag(lehmann_result%susceptibility)))
      if (.not. finite_response) error stop 'TDVK-03 closure: Lehmann reference contains NaN or Inf'

      base_integration_eta = config%gf_integration_eta
      if (base_integration_eta <= 0.0_rp) base_integration_eta = config%eta/40.0_rp
      if (base_integration_eta <= 0.0_rp .or. base_integration_eta >= config%eta) then
         error stop 'TDVK-03 closure: invalid base integration_eta'
      end if
      write (*, '(a)') 'TDVK-03 Fe reciprocal-backend closure: one accepted state, complete compact product space'
      write (*, '(a,i0,a,i0,a,i0,a,es16.8,a,es16.8,a,a)') '  accepted_state nbasis=', left_state%nbasis, &
         ' nbands=', left_state%nbands, ' nk=', left_state%nk, ' EF_Ry=', left_state%fermi_level, &
         ' temperature_K=', left_state%temperature, ' moment provenance=accepted LR-01 radial snapshot'
      write (*, '(a,3(es16.8,1x),a,a,a,i0)') '  q=', config%q_list(:, gamma_index), ' channel=', trim(config%channel), &
         ' product_dimension=', product_plus%product_dimension

      ! Width ladder.  The point count is derived from the same accepted
      ! energy span for each width so h/integration_eta remains controlled.
      do i = 1, size(eta_scale)
         integration_eta = eta_scale(i)*base_integration_eta
         integration_points = resolved_simpson_points(left_state, endpoints(gamma_index), integration_eta, &
            config%gf_energy_margin)
         if (i == size(eta_scale)) integration_points = max(integration_points, config%gf_integration_points)
         call report_product_gf_closure_sample('eta_ladder', config, response_space, product_plus, product_minus, &
            left_state, endpoints(gamma_index), lehmann_result, integration_points, integration_eta, config%gf_energy_margin)
      end do

      ! Simpson-resolution ladder at the narrowest controlled width.  Both
      ! grids are resolved; their difference isolates remaining mesh error.
      integration_points = resolved_simpson_points(left_state, endpoints(gamma_index), base_integration_eta, &
         config%gf_energy_margin)
      integration_points = max(integration_points, config%gf_integration_points)
      call report_product_gf_closure_sample('simpson_base', config, response_space, product_plus, product_minus, &
         left_state, endpoints(gamma_index), lehmann_result, integration_points, base_integration_eta, config%gf_energy_margin)
      call report_product_gf_closure_sample('simpson_fine', config, response_space, product_plus, product_minus, &
         left_state, endpoints(gamma_index), lehmann_result, 2*(integration_points - 1) + 1, base_integration_eta, &
         config%gf_energy_margin)

      ! Energy-window ladder.  Choose a resolved grid independently at every
      ! margin so the window comparison is not contaminated by coarse Simpson
      ! spacing.
      do i = 1, size(margin_increment)
         margin = config%gf_energy_margin + margin_increment(i)
         integration_points = resolved_simpson_points(left_state, endpoints(gamma_index), base_integration_eta, margin)
         call report_product_gf_closure_sample('window_ladder', config, response_space, product_plus, product_minus, &
            left_state, endpoints(gamma_index), lehmann_result, integration_points, base_integration_eta, margin)
      end do
      write (*, '(a)') 'TDVK-03 Fe reciprocal-backend closure: audit complete; evidence returned without physical interpretation'
   end subroutine run_tddft_product_gf_closure

   integer function resolved_simpson_points(left_state, right_state, integration_eta, margin) result(points)
      type(lr_electronic_state), intent(in) :: left_state, right_state
      real(rp), intent(in) :: integration_eta, margin
      real(rp) :: span
      integer :: intervals

      span = max(maxval(left_state%eigenvalues), maxval(right_state%eigenvalues)) - &
         min(minval(left_state%eigenvalues), minval(right_state%eigenvalues)) + 2.0_rp*margin
      intervals = ceiling(span/(0.4_rp*integration_eta))
      if (intervals < 2) intervals = 2
      if (mod(intervals, 2) /= 0) intervals = intervals + 1
      points = intervals + 1
   end function resolved_simpson_points

   subroutine report_product_gf_closure_sample(tag, config, response_space, product_plus, product_minus, left_state, endpoint, &
                                               lehmann_result, integration_points, integration_eta, energy_margin)
      character(len=*), intent(in) :: tag
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), intent(in) :: response_space
      type(lmto_product_response_basis), target, intent(in) :: product_plus, product_minus
      type(lr_electronic_state), target, intent(in) :: left_state, endpoint
      type(lr_product_ks_susceptibility_result), intent(in) :: lehmann_result
      integer, intent(in) :: integration_points
      real(rp), intent(in) :: integration_eta, energy_margin
      type(lr_product_gf_susceptibility_request) :: gf_request
      type(lr_product_gf_susceptibility_result) :: gf_result
      real(rp) :: norm_lehmann, norm_gf, difference_frobenius, relative_frobenius, difference_infinity
      logical :: finite_response
      type(lmto_product_response_basis), pointer :: product

      if (trim(config%channel) == 'chi_plus') then
         product => product_plus
      else
         product => product_minus
      end if
      gf_request%q = config%q_list(:, find_gamma_q(config%q_list))
      gf_request%frequencies = [0.0_rp]
      gf_request%eta = config%eta
      gf_request%channel = config%channel
      gf_request%integration_points = integration_points
      gf_request%integration_eta = integration_eta
      gf_request%energy_margin = energy_margin
      gf_request%product_basis => product
      gf_request%electronic_state => left_state
      gf_request%q_endpoint_state => endpoint
      call evaluate_lr_product_gf_susceptibility(gf_request, gf_result)
      finite_response = all(ieee_is_finite(real(gf_result%susceptibility, rp))) .and. &
         all(ieee_is_finite(aimag(gf_result%susceptibility)))
      if (.not. finite_response) error stop 'TDVK-03 closure: GF response contains NaN or Inf'
      norm_lehmann = sqrt(sum(abs(lehmann_result%susceptibility(:, :, 1))**2))
      norm_gf = sqrt(sum(abs(gf_result%susceptibility(:, :, 1))**2))
      difference_frobenius = sqrt(sum(abs(lehmann_result%susceptibility(:, :, 1) - &
         gf_result%susceptibility(:, :, 1))**2))
      relative_frobenius = difference_frobenius/max(norm_lehmann, norm_gf, tiny(1.0_rp))
      difference_infinity = maxval(abs(lehmann_result%susceptibility(:, :, 1) - gf_result%susceptibility(:, :, 1)))
      write (*, '(a,a,a,i0,4(a,es16.8))') '  ', trim(tag), ' N=', integration_points, &
         ' eta=', config%eta, ' integration_eta=', integration_eta, ' margin=', energy_margin, &
         ' energy_min=', gf_result%energy_min
      write (*, '(a,8(a,es16.8))') '    ', ' energy_max=', gf_result%energy_max, ' h=', gf_result%energy_spacing, &
         ' h_over_eta=', gf_result%spacing_over_integration_eta, ' norm_lehmann=', norm_lehmann, &
         ' norm_gf=', norm_gf, ' dF=', difference_frobenius, ' rF=', relative_frobenius, &
         ' dInf=', difference_infinity
      write (*, '(a,a,es16.8)') '    ', ' wall_s=', gf_result%wall_time_seconds
   end subroutine report_product_gf_closure_sample

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

      if (size(point_transition) /= response_space%ndim .or. size(product_coordinates) /= product%product_dimension) then
         error stop 'TDVK-02R2 transition projection: vector shape mismatch'
      end if
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
      type(tddft_production_config), intent(in) :: config
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

      if (any(shape(pauli_magnetization) /= [response_space%nsite, response_space%npoint])) then
         error stop 'TDDFT production driver: accepted Pauli magnetization shape mismatch'
      end if
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

   !> Evaluate an already prepared request batch. This is the reproducibility
   !> seam: tests and validation can compare it directly with service calls,
   !> without reconstructing SCF state.
   subroutine evaluate_tddft_production_sweep(config, response_space, radial_bases, ground_states, left_state, endpoints, result, &
                                              native_provider, native_pairs, native_site_positions, accepted_pauli_magnetization, &
                                              accepted_pauli_magnetization_source)
      type(tddft_production_config), intent(in) :: config
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
      result%reciprocal_backend_crosscheck = config%reciprocal_backend_crosscheck
      result%bare_response_provenance = ''
      result%q_list = config%q_list
      result%frequencies = config%frequencies
      allocate(result%ks_susceptibility(ndim, ndim, result%nfrequency, result%nq), &
               result%enhanced_susceptibility(ndim, ndim, result%nfrequency, result%nq), &
               result%loss_matrix(ndim, ndim, result%nfrequency, result%nq))
      result%ks_susceptibility = cmplx(0.0_rp, 0.0_rp, rp)
      result%enhanced_susceptibility = cmplx(0.0_rp, 0.0_rp, rp)
      result%loss_matrix = cmplx(0.0_rp, 0.0_rp, rp)
      if (config%reciprocal_backend_crosscheck) then
         allocate(result%reciprocal_crosscheck_valid(result%nfrequency, result%nq), &
            result%reciprocal_crosscheck_norm_lehmann(result%nfrequency, result%nq), &
            result%reciprocal_crosscheck_norm_gf(result%nfrequency, result%nq), &
            result%reciprocal_crosscheck_difference_frobenius(result%nfrequency, result%nq), &
            result%reciprocal_crosscheck_relative_frobenius(result%nfrequency, result%nq), &
            result%reciprocal_crosscheck_difference_infinity(result%nfrequency, result%nq), &
            result%reciprocal_crosscheck_delta(ndim, ndim, result%nfrequency, result%nq))
         result%reciprocal_crosscheck_valid = .false.
         result%reciprocal_crosscheck_norm_lehmann = 0.0_rp
         result%reciprocal_crosscheck_norm_gf = 0.0_rp
         result%reciprocal_crosscheck_difference_frobenius = 0.0_rp
         result%reciprocal_crosscheck_relative_frobenius = 0.0_rp
         result%reciprocal_crosscheck_difference_infinity = 0.0_rp
         result%reciprocal_crosscheck_delta = cmplx(0.0_rp, 0.0_rp, rp)
      end if

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
         if (config%reciprocal_backend_crosscheck) then
            call evaluate_reciprocal_backend_crosscheck(config, response_space, radial_bases, left_state, endpoints(iq), &
               config%q_list(:, iq), config%frequencies, bare_result, result%reciprocal_crosscheck_valid(:, iq), &
               result%reciprocal_crosscheck_norm_lehmann(:, iq), result%reciprocal_crosscheck_norm_gf(:, iq), &
               result%reciprocal_crosscheck_difference_frobenius(:, iq), result%reciprocal_crosscheck_relative_frobenius(:, iq), &
               result%reciprocal_crosscheck_difference_infinity(:, iq), result%reciprocal_crosscheck_delta(:, :, :, iq))
         end if
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
   end subroutine evaluate_tddft_production_sweep

   subroutine validate_prepared_sweep(config, response_space, radial_bases, ground_states, left_state, endpoints, &
                                      native_provider, native_pairs)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      class(lr_rs_gf_provider), intent(in), optional :: native_provider
      type(lr_rs_gf_pair), intent(in), optional :: native_pairs(:)
      integer :: iq

      call validate_tddft_config(config)
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
      type(tddft_production_config), intent(in) :: config
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
      else if (trim(config%backend) == tddft_driver_backend_reciprocal_gf) then
         call evaluate_reciprocal_gf_backend(config, response_space, radial_bases, left_state, endpoint, q, frequencies, result)
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

   subroutine evaluate_reciprocal_gf_backend(config, response_space, radial_bases, left_state, endpoint, q, frequencies, result)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(lr_electronic_state), target, intent(in) :: left_state, endpoint
      real(rp), intent(in) :: q(3), frequencies(:)
      type(lr_ks_susceptibility_result), intent(out) :: result
      type(lr_gf_susceptibility_request) :: request

      request%q = q
      request%frequencies = frequencies
      request%eta = config%eta
      request%channel = config%channel
      request%integration_points = config%gf_integration_points
      request%integration_eta = config%gf_integration_eta
      request%energy_margin = config%gf_energy_margin
      request%response_space => response_space
      request%radial_bases => radial_bases
      request%electronic_state => left_state
      request%q_endpoint_state => endpoint
      call evaluate_lr_gf_susceptibility(request, result)
   end subroutine evaluate_reciprocal_gf_backend

   subroutine evaluate_reciprocal_backend_crosscheck(config, response_space, radial_bases, left_state, endpoint, q, frequencies, &
                                                     selected_result, valid, norm_lehmann, norm_gf, difference_frobenius, &
                                                     relative_frobenius, difference_infinity, delta)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(lr_electronic_state), target, intent(in) :: left_state, endpoint
      real(rp), intent(in) :: q(3), frequencies(:)
      type(lr_ks_susceptibility_result), intent(in) :: selected_result
      logical, intent(out) :: valid(:)
      real(rp), intent(out) :: norm_lehmann(:), norm_gf(:), difference_frobenius(:), relative_frobenius(:), difference_infinity(:)
      complex(rp), intent(out) :: delta(:, :, :)
      type(lr_ks_susceptibility_result) :: lehmann_result, gf_result
      integer :: ifrequency

      if (trim(config%backend) == tddft_driver_backend_lehmann .or. trim(config%backend) == 'spectral') then
         lehmann_result = selected_result
         call evaluate_reciprocal_gf_backend(config, response_space, radial_bases, left_state, endpoint, q, frequencies, gf_result)
      else if (trim(config%backend) == tddft_driver_backend_reciprocal_gf) then
         call evaluate_lehmann_backend(response_space, radial_bases, left_state, endpoint, q, frequencies, config%eta, &
            config%channel, lehmann_result)
         gf_result = selected_result
      else
         call evaluate_lehmann_backend(response_space, radial_bases, left_state, endpoint, q, frequencies, config%eta, &
            config%channel, lehmann_result)
         call evaluate_reciprocal_gf_backend(config, response_space, radial_bases, left_state, endpoint, q, frequencies, gf_result)
      end if

      if (size(valid) /= size(frequencies) .or. size(norm_lehmann) /= size(frequencies) .or. &
          size(norm_gf) /= size(frequencies) .or. size(difference_frobenius) /= size(frequencies) .or. &
          size(relative_frobenius) /= size(frequencies) .or. size(difference_infinity) /= size(frequencies) .or. &
          any(shape(delta) /= [response_space%ndim, response_space%ndim, size(frequencies)])) then
         error stop 'TDDFT production driver: crosscheck output shape mismatch'
      end if

      do ifrequency = 1, size(frequencies)
         delta(:, :, ifrequency) = lehmann_result%susceptibility(:, :, ifrequency) - gf_result%susceptibility(:, :, ifrequency)
         difference_frobenius(ifrequency) = sqrt(sum(abs(delta(:, :, ifrequency))**2))
         norm_lehmann(ifrequency) = sqrt(sum(abs(lehmann_result%susceptibility(:, :, ifrequency))**2))
         norm_gf(ifrequency) = sqrt(sum(abs(gf_result%susceptibility(:, :, ifrequency))**2))
         relative_frobenius(ifrequency) = difference_frobenius(ifrequency)/ &
            max(norm_lehmann(ifrequency), norm_gf(ifrequency), tiny(1.0_rp))
         difference_infinity(ifrequency) = maxval(abs(delta(:, :, ifrequency)))
      end do
      valid = .true.
   end subroutine evaluate_reciprocal_backend_crosscheck

   subroutine radial_basis_from_snapshot(state, basis, complete_sr)
      type(radial_ground_state), intent(in) :: state
      type(lmto_radial_basis), intent(out) :: basis
      logical, intent(in), optional :: complete_sr
      logical :: use_complete_sr

      use_complete_sr = .false.
      if (present(complete_sr)) use_complete_sr = complete_sr

      if (.not. allocated(state%r) .or. .not. allocated(state%pauli_large) .or. &
          .not. allocated(state%pauli_large_dot) .or. .not. allocated(state%pauli_enu)) then
         error stop 'TDDFT production driver: accepted radial snapshot lacks Pauli basis arrays'
      end if
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

   subroutine write_tddft_production_output(config, result, control_obj, lattice_obj, ground_states, response_space, reciprocal_obj)
      type(tddft_production_config), intent(in) :: config
      type(tddft_production_result), intent(in) :: result
      type(control), intent(in) :: control_obj
      type(lattice), intent(in) :: lattice_obj
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(response_space_layout), intent(in) :: response_space
      type(reciprocal), intent(in) :: reciprocal_obj
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
      write(unit, '(a,i0,2(es24.16,1x),a,es24.16,a,es24.16)') '# radial_mesh_identity = ', size(ground_states(1)%r), &
         ground_states(1)%a, ground_states(1)%b, 'rmax=', ground_states(1)%rmax, 'sum_r=', sum(ground_states(1)%r)
      write(unit, '(a,a)') '# interaction_route = ', trim(config%interaction_route)
      write(unit, '(a,a)') '# interaction_provenance = ', trim(result%interaction_provenance)
      write(unit, '(a,a)') '# magnetization_kind = ', trim(result%magnetization_kind)
      write(unit, '(a,a)') '# magnetization_source = ', trim(result%magnetization_source)
      write(unit, '(a,a)') '# goldstone_correction = ', trim(result%goldstone_correction_status)
      write(unit, '(a,a)') '# backend = ', trim(config%backend)
      if (config%reciprocal_backend_crosscheck) then
         write(unit, '(a)') '# reciprocal_backend_crosscheck = enabled (validation diagnostic; selected backend remains authoritative)'
         write(unit, '(a)') '# crosscheck_metrics columns: q_index omega_Ry valid norm_lehmann norm_reciprocal_gf difference_frobenius relative_frobenius difference_infinity'
         do iq = 1, result%nq
            do iw = 1, result%nfrequency
               write(unit, '(a,i0,1x,es24.16,1x,l1,5(1x,es24.16))') '# crosscheck_metrics ', iq, result%frequencies(iw), &
                  result%reciprocal_crosscheck_valid(iw, iq), result%reciprocal_crosscheck_norm_lehmann(iw, iq), &
                  result%reciprocal_crosscheck_norm_gf(iw, iq), result%reciprocal_crosscheck_difference_frobenius(iw, iq), &
                  result%reciprocal_crosscheck_relative_frobenius(iw, iq), result%reciprocal_crosscheck_difference_infinity(iw, iq)
            end do
         end do
         write(unit, '(a)') '# crosscheck_delta columns: q_index omega_Ry matrix_i matrix_j delta_real delta_imag'
         do iq = 1, result%nq
            do iw = 1, result%nfrequency
               do j = 1, result%ndim
                  do i = 1, result%ndim
                     write(unit, '(a,i0,1x,es24.16,1x,2(i0,1x),2(es24.16,1x))') '# crosscheck_delta ', iq, &
                        result%frequencies(iw), i, j, real(result%reciprocal_crosscheck_delta(i, j, iw, iq), rp), &
                        aimag(result%reciprocal_crosscheck_delta(i, j, iw, iq))
                  end do
               end do
            end do
         end do
      end if
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

      if (any(shape(left) /= shape(right))) error stop 'DRESP-06A site difference: shape mismatch'
      left_norm = sqrt(sum(abs(left)**2))
      right_norm = sqrt(sum(abs(right)**2))
      value = sqrt(sum(abs(left - right)**2))/max(left_norm, right_norm, tiny(1.0_rp))
   end function relative_site_difference

end module tddft_production_driver_mod
