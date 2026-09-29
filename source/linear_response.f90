!> Linear-response services and (from LR-REF-02b) the linear_response class.
!> Rotation dynamics: second-order local-rotation response of the accepted H2.
module linear_response_mod
   use precision_mod, only: rp
   use math_mod, only: pi
   use control_mod, only: control
   use energy_mod, only: energy_type => energy
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use reciprocal_mod, only: reciprocal
   use self_mod, only: self
   use string_mod, only: sl
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use radial_ground_state_mod, only: radial_ground_state, radial_xc_provenance
   use, intrinsic :: iso_fortran_env, only: int64
   use green_mod, only: green
   use recursion_mod, only: recursion
   use lmto_path_operator_mod, only: native_turek_contour_report
   implicit none
   private

   ! --- from lr_kl_hessian ---
   real(rp), parameter, public :: force_theorem_pi = pi
   type, public :: lmto_live_hamiltonian_fixture
      integer :: nsite = 0
      integer :: norb = 0
      integer :: nbond = 0
      logical :: hoh = .false.
      logical :: include_enu = .true.
      logical :: cartesian_to_spherical = .false.
      complex(rp), allocatable :: hhh(:, :, :)       ! orbital directed bonds
      integer, allocatable :: bond_source(:), bond_target(:)
      real(rp), allocatable :: bond_vector(:, :)     ! Cartesian lattice vector
      logical, allocatable :: onsite(:)
      real(rp), allocatable :: moments(:, :)         ! (3,nsite), unit moments
      ! Fractional basis positions.  A bond vector is the full target-source
      ! displacement R+tau_target-tau_source.  The endpoint phases below are
      ! relative to the source endpoint because the Bloch basis already carries
      ! the absolute basis position.
      real(rp), allocatable :: site_position(:, :)   ! (3,nsite), fractional
      complex(rp), allocatable :: wx0(:, :), wx1(:, :) ! (norb,nsite)
      complex(rp), allocatable :: c0(:, :), c1(:, :)   ! live onsite c channels
      complex(rp), allocatable :: obar0(:, :), obar1(:, :) ! live O channels
      complex(rp), allocatable :: enu0(:, :), enu1(:, :)   ! live e_nu channels
   contains
      procedure :: clear => lmto_fixture_clear
   end type lmto_live_hamiltonian_fixture

   ! --- from lr_kl_contour ---
   type, public :: finite_h_contour_options
      integer :: contour_points = 32
      character(len=16) :: contour_shape = 'ellipse'
      real(rp) :: contour_margin = 0.25_rp
      real(rp) :: contour_height_fraction = 0.35_rp
      logical :: account_fermi_poles = .true.
   end type finite_h_contour_options
   type, public :: finite_h_contour_report
      integer :: contour_points = 0
      integer :: fermi_poles = 0
      real(rp) :: solve_seconds = 0.0_rp
      real(rp) :: contour_seconds = 0.0_rp
      real(rp) :: pole_seconds = 0.0_rp
   end type finite_h_contour_report

   ! --- from lr_rotation_response ---
   type, public :: rotation_state
      type(lmto_live_hamiltonian_fixture), pointer :: fixture => null()
      type(reciprocal), pointer :: reciprocal_state => null()
      real(rp) :: q(3) = 0.0_rp
      integer :: nsite = 0, ncoord = 0, nmat = 0, nbands = 0, nk = 0
      real(rp) :: fermi = 0.0_rp, kT = 0.0_rp, weight_sum = 0.0_rp
      real(rp) :: magnetization = 0.0_rp, berry = 0.0_rp, berry_residual = 0.0_rp
      real(rp) :: electron_count = 0.0_rp, electron_residual = 0.0_rp
      logical :: prepared = .false.
      real(rp), allocatable :: values(:, :), endpoint_values(:, :), weights(:)
      real(rp), allocatable :: occupations(:, :), endpoint_occupations(:, :)
      complex(rp), allocatable :: torque_band(:, :, :, :), partner_band(:, :, :, :)
      complex(rp), allocatable :: contact(:, :)
   contains
      procedure :: clear => rotation_state_clear
   end type rotation_state
   type, public :: rotation_request
      type(rotation_state), pointer :: state => null()
      real(rp) :: omega = 0.0_rp
      real(rp) :: eta = 0.0_rp
      logical :: exact_static = .false.
      logical :: want_inverse = .false.
   end type rotation_request
   type, public :: rotation_result
      real(rp) :: q(3) = 0.0_rp
      real(rp) :: omega = 0.0_rp, eta = 0.0_rp
      real(rp) :: berry = 0.0_rp, magnetization = 0.0_rp, berry_residual = 0.0_rp
      real(rp) :: electron_count = 0.0_rp, electron_residual = 0.0_rp
      logical :: inverse_available = .false.
      complex(rp), allocatable :: bubble(:, :), contact(:, :), kernel(:, :), kernel_pm(:, :)
      complex(rp), allocatable :: inverse_kernel(:, :), inverse_kernel_pm(:, :)
   contains
      procedure :: clear => rotation_result_clear
   end type rotation_result

   integer, parameter, public :: linear_response_max_points = 2000

   type, public :: linear_response_config
      character(len=sl) :: fname = ''
      character(len=16) :: formulation = 'rotation'
      character(len=16) :: representation = 'radial_points'
      character(len=16) :: bare_response = 'lehmann'
      character(len=16) :: realspace_solver = 'auto'
      character(len=16) :: interaction = 'alsda'
      character(len=8) :: projection = 'spd'
      character(len=16) :: diagnostics = 'none'
      character(len=32) :: channel = 'chi_plus'
      character(len=16) :: q_coordinates = 'direct'
      character(len=sl) :: q_file = ''
      integer :: n_q = 1
      integer :: n_q_points = 0
      real(rp), allocatable :: q_list(:, :)
      integer :: n_omega = 1
      logical :: use_omega_grid = .false.
      real(rp), allocatable :: omega_grid(:)
      real(rp) :: omega_min = 0.0_rp
      real(rp) :: omega_max = 0.0_rp
      real(rp) :: eta = 0.01_rp
      integer :: n_eta = 1
      real(rp), allocatable :: eta_grid(:)
      integer :: response_lmax = -1
      logical :: write_full_matrix = .true.
      real(rp) :: rotation_axis(3) = [1.0_rp, 0.0_rp, 0.0_rp]
      character(len=24) :: finite_h_spectral_mode = 'metallic'
      character(len=16) :: finite_h_response_backend = 'spectral'
      integer :: contour_points = 32
      character(len=16) :: contour_shape = 'ellipse'
      real(rp) :: contour_margin = 0.25_rp
      real(rp) :: contour_height_fraction = 0.35_rp
      logical :: contour_account_fermi_poles = .true.
      logical :: native_turek = .false.
      logical :: native_crosscheck = .false.
      real(rp) :: native_green_eta = 1.0e-3_rp
      integer :: native_energy_points = 0
      integer :: native_contour_points = 64
      real(rp) :: native_contour_margin = 0.25_rp
      real(rp) :: native_contour_height_fraction = 0.35_rp
      logical :: native_contour_account_fermi_poles = .true.
      integer :: native_contour_target_fermi_poles = 0
      character(len=32) :: native_rsgf_provider = 'auto'
      integer :: gf_integration_points = 2001
      real(rp) :: gf_integration_eta = 0.0_rp
      real(rp) :: gf_energy_margin = 1.0_rp
      character(len=sl) :: output_file = 'rotation_dynamics.dat'
      ! Rotation campaign controls (energies in Ry; defaults preserve the
      ! LR-REF-02b certification campaign bit-for-bit):
      ! rotation_eta_ladder=[1e-4,2.5e-5,6.25e-6] Ry;
      ! rotation_probe_omega=1e-5 Ry (the canonical slope step);
      ! rotation_probe_eta=1e-9 Ry; rotation_slope_step is a deprecated
      ! input-compatibility field and is ignored; rotation_pole_window_floor=1e-3 Ry;
      ! rotation_pole_window_scale=2.5; rotation_pole_window_max=2e-2 Ry;
      ! rotation_pole_coarse_points=61; rotation_pole_fine_points=41;
      ! rotation_pole_refinement_half_width=2.0 coarse-grid steps.
      real(rp) :: rotation_eta_ladder(3) = [1.0e-4_rp, 2.5e-5_rp, 6.25e-6_rp]
      real(rp) :: rotation_probe_omega = 1.0e-5_rp
      real(rp) :: rotation_probe_eta = 1.0e-9_rp
      real(rp) :: rotation_slope_step = 1.0e-5_rp
      real(rp) :: rotation_pole_window_floor = 1.0e-3_rp
      real(rp) :: rotation_pole_window_scale = 2.5_rp
      real(rp) :: rotation_pole_window_max = 2.0e-2_rp
      integer :: rotation_pole_coarse_points = 61
      integer :: rotation_pole_fine_points = 41
      real(rp) :: rotation_pole_refinement_half_width = 2.0_rp
      ! Optional frequency grid: rotation_n_omega=0 disables the new file;
      ! otherwise rotation_omega_min=0 Ry, rotation_omega_max=2e-2 Ry,
      ! rotation_eta=1e-4 Ry, rotation_grid_file='rotation_response_grid.dat'.
      integer :: rotation_n_omega = 0
      real(rp) :: rotation_omega_min = 0.0_rp
      real(rp) :: rotation_omega_max = 2.0e-2_rp
      real(rp) :: rotation_eta = 1.0e-4_rp
      character(len=sl) :: rotation_grid_file = 'rotation_response_grid.dat'
   contains
      procedure :: restore_to_default => lr_config_restore_to_default
   end type linear_response_config

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
      character(len=128) :: status = 'not evaluated'
      character(len=48) :: interaction_route = ''
      character(len=32) :: backend = ''
      character(len=64) :: goldstone_correction_status = 'not selected'
      character(len=256) :: interaction_provenance = ''
      character(len=256) :: bare_response_provenance = ''
      character(len=32) :: magnetization_kind = ''
      character(len=256) :: magnetization_source = ''
   end type tddft_production_result

   public :: tddft_capability_is_supported, require_tddft_capability
   public :: evaluate_linear_response_sweep

   type, public :: linear_response
      type(linear_response_config) :: config
   contains
      procedure :: restore_to_default => lr_restore_to_default
      procedure :: load_config => lr_load_config
      procedure :: validate_capability => lr_validate_capability
      procedure :: run => lr_run
   end type linear_response


   ! --- from response_angular_basis_mod ---

   real(rp), parameter, public :: response_angular_pi = &
      3.1415926535897932384626433832795_rp

   ! --- from response_basis_mapping_mod ---

   !> Sign of the live reciprocal/Fourier convention exp(+i 2*pi*k.R).
   integer, parameter, public :: response_fourier_phase_sign = 1

   type, public :: response_super_index
      integer :: site = 0
      integer :: response_l = 0
      integer :: response_m = 0
      integer :: radial_point = 0
      integer :: channel = 0
   end type response_super_index

   ! --- from lr_response_space_mod ---

   type, public :: response_space_layout
      integer :: nsite = 0
      integer :: response_lmax = -1
      integer :: npoint = 0
      integer :: nchannel = 0
      integer :: ndim = 0
      integer :: active_dimension = 0
      real(rp) :: a = 0.0_rp
      real(rp) :: b = 0.0_rp
      real(rp), allocatable :: radius(:)
      real(rp), allocatable :: radial_weights(:)
      real(rp), allocatable :: metric_weights(:)
   contains
      procedure :: initialize => response_space_initialize
   end type response_space_layout

   ! --- from lr_lmto_endpoint_branches_mod ---

   integer, parameter, public :: lmto_product_branch_00 = 1
   integer, parameter, public :: lmto_product_branch_10 = 2
   integer, parameter, public :: lmto_product_branch_01 = 3
   integer, parameter, public :: lmto_product_branch_11 = 4
   integer, parameter, public :: lmto_product_branch_20 = 5
   integer, parameter, public :: lmto_product_branch_02 = 6
   integer, parameter, public :: lmto_product_nbranch = 6
   integer, parameter, public :: lmto_product_max_endpoint_power = 2
   integer, parameter, public :: lmto_product_max_gf_moment = 4

   ! --- from lr_pauli_transition_vertex_mod ---

   !> Capability tuple for the initial LR-05 production vertex.
   type, public :: pauli_vertex_capabilities
      character(len=32) :: reciprocal_mode = 'ham_only'
      character(len=16) :: hamiltonian_order = 'second'
      logical :: orthogonal = .true.
      logical :: collinear = .true.
      logical :: has_soc = .false.
      logical :: has_extra_operator = .false.
   end type pauli_vertex_capabilities

   !> One reciprocal endpoint eigenstate in production site-blocked ordering.
   !>
   !> The coefficient vector is packed as
   !> `(site, spin-up orbitals, spin-down orbitals)`, with the orbital order
   !> `(s),(p,-1:1),(d,-2:2),...` supplied by `basis_mod`/LMTO.
   type, public :: pauli_endpoint_state
      real(rp) :: energy = 0.0_rp
      complex(rp), allocatable :: coefficients(:)
   contains
      procedure :: initialize => pauli_endpoint_state_initialize
   end type pauli_endpoint_state

   ! --- from lr_lmto_product_response_basis_mod ---

   integer, parameter, public :: lmto_product_channel_plus = 1
   integer, parameter, public :: lmto_product_channel_minus = 2

   type, public :: lmto_product_candidate
      integer :: l = 0
      integer :: lp = 0
      integer :: p = 0
      integer :: q = 0
      integer :: branch = 0
   end type lmto_product_candidate

   type, public :: lmto_product_block
      integer :: site = 0
      integer :: response_l = 0
      integer :: channel = 0
      integer :: npoint = 0
      integer :: ncandidate = 0
      integer :: nsv = 0
      integer :: rank = 0
      integer :: rank_tau1 = 0
      integer :: rank_tau10 = 0
      integer :: rank_tau100 = 0
      real(rp) :: tau1 = 0.0_rp
      real(rp) :: tau10 = 0.0_rp
      real(rp) :: tau100 = 0.0_rp
      logical :: rank_stable = .false.
      type(lmto_product_candidate), allocatable :: candidates(:)
      real(rp), allocatable :: column_norms(:)
      real(rp), allocatable :: singular_values(:)
      !> Weighted orthonormal radial modes U, with columns retained by rank.
      complex(rp), allocatable :: weighted_modes(:, :)
      !> Retained rows of V^H from the direct SVD, kept for reconstruction
      !> audits and future operator projection.
      complex(rp), allocatable :: right_modes(:, :)
      !> Forward candidate-to-coordinate map F = Sigma V^H D.
      complex(rp), allocatable :: forward_transform(:, :)
   end type lmto_product_block

   type, public :: lmto_product_response_basis
      integer :: nsite = 0
      integer :: orbital_lmax = -1
      integer :: response_lmax = -1
      integer :: circular_channel = 0
      integer :: npoint = 0
      integer :: unpruned_dimension = 0
      integer :: product_dimension = 0
      character(len=16) :: product_radial_order = 'second'
      integer :: product_endpoint_branches = lmto_product_nbranch
      character(len=32) :: product_branch_labels = '00,10,01,11,20,02'
      integer :: maximum_gf_energy_moment = 4
      type(lmto_product_block), allocatable :: blocks(:, :)
   contains
      procedure :: initialize => lmto_product_response_basis_initialize
      procedure :: flat_index => lmto_product_flat_index
      procedure :: unflatten_index => lmto_product_unflatten_index
      procedure :: candidate_coefficients => lmto_product_candidate_coefficients
      procedure :: transition_coordinates => lmto_product_transition_coordinates
      procedure :: component_vertex_tensor => lmto_product_component_vertex_tensor
   end type lmto_product_response_basis

   ! --- from lr_gf_endpoint_augmentation_mod ---

   character(len=*), parameter, public :: lr_gf_endpoint_representation = &
      'Pauli/no-SOC spatial GF, factored LMTO augmentation'

   !> Capability tuple accepted by this endpoint adapter.
   type, public :: lr_gf_endpoint_capabilities
      character(len=32) :: reciprocal_mode = 'ham_only'
      character(len=16) :: hamiltonian_order = 'second'
      logical :: orthogonal = .true.
      logical :: collinear = .true.
      logical :: has_soc = .false.
      logical :: has_extra_operator = .false.
      ! These two explicit names make the fail-closed boundary readable to
      ! callers that distinguish Hubbard from other additive operators.
      logical :: has_hubbard = .false.
      logical :: has_additive_operator = .false.
   contains
      procedure :: supported => lr_gf_endpoint_capabilities_supported
      procedure :: require_supported => lr_gf_endpoint_require_supported
   end type lr_gf_endpoint_capabilities

   !> Callback for a production effective-Hamiltonian action.
   !>
   !> `seed` is a localized source-site block.  The implementation applies the
   !> same H_eff action used by the native solver and returns the destination
   !> site block.  For destination=source the caller subtracts E_nu_work to
   !> obtain h^gamma_aa; for an offsite block the returned H block is already
   !> h^gamma_ab.
   type, abstract, public :: lr_gf_effective_hamiltonian_provider
   contains
      procedure(lr_gf_apply_effective_action), deferred :: apply
   end type lr_gf_effective_hamiltonian_provider

   abstract interface
      subroutine lr_gf_apply_effective_action(this, source_site, destination_site, seed, destination_block)
         import :: lr_gf_effective_hamiltonian_provider, rp
         class(lr_gf_effective_hamiltonian_provider), intent(inout) :: this
         integer, intent(in) :: source_site, destination_site
         complex(rp), intent(in) :: seed(:, :)
         complex(rp), intent(out) :: destination_block(:, :)
      end subroutine lr_gf_apply_effective_action
   end interface

   !> Simple dense provider used by small finite fixtures and by callers that
   !> already hold the completed H_eff matrix.  Production callers can provide
   !> a wrapper around their native two-sweep/action routine instead.
   type, extends(lr_gf_effective_hamiltonian_provider), public :: lr_gf_dense_hamiltonian_provider
      complex(rp), allocatable :: h_eff(:, :)
      integer :: block_size = 0
   contains
      procedure :: initialize => lr_gf_dense_provider_initialize
      procedure :: apply => lr_gf_dense_provider_apply
   end type lr_gf_dense_hamiltonian_provider

   !> Factored endpoint-augmented Green-function block.
   !>
   !> The four coefficient-space members are the explicit radial branches in
   !>   Phi G Phi^dagger,
   !>   Phidot (hG) Phi^dagger,
   !>   Phi (Gh) Phidot^dagger,
   !>   Phidot (hGh) Phidot^dagger.
   !>
   !> `phi_left/right` and `phidot_left/right` contain large-component radial
   !> numerators repeated in orbital/spin ordering.  Angular harmonics are not
   !> duplicated in storage; `evaluate_point` inserts the certified
   !> response_angular_basis convention and divides U_l by r at the endpoint.
   type, public :: lr_gf_augmented_block
      integer :: left_site = 0
      integer :: right_site = 0
      integer :: norb = 0
      integer :: nspin = 0
      complex(rp) :: z = (0.0_rp, 0.0_rp)
      complex(rp), allocatable :: coefficient_gf(:, :)
      complex(rp), allocatable :: hgamma_ab(:, :)
      complex(rp), allocatable :: h_g(:, :)
      complex(rp), allocatable :: g_h(:, :)
      complex(rp), allocatable :: h_g_h(:, :)
      ! Shape: (radial point, orbital, spin), large component only.
      real(rp), allocatable :: phi_left(:, :, :), phidot_left(:, :, :)
      real(rp), allocatable :: phi_right(:, :, :), phidot_right(:, :, :)
      real(rp), allocatable :: radius_left(:), radius_right(:)
   contains
      procedure :: evaluate_point => lr_gf_augmented_block_evaluate_point
      procedure :: branch_at_point => lr_gf_augmented_block_branch_at_point
   end type lr_gf_augmented_block

   ! --- from lr_projected_site_spin_mod ---


   integer, parameter, public :: projected_selector_d = 1
   integer, parameter, public :: projected_selector_spd = 2

   integer, parameter, public :: projected_operator_plus = 1
   integer, parameter, public :: projected_operator_minus = 2
   integer, parameter, public :: projected_operator_z = 3

   character(len=1), parameter, public :: projected_selection_d = 'd'
   character(len=3), parameter, public :: projected_selection_spd = 'spd'
   character(len=40), parameter, public :: projected_core_policy = &
      'valence-only; frozen core excluded'

   !> Immutable-after-initialization DRESP-01 contract metadata and operations.
   type, public :: projected_site_spin_contract
      integer :: nsite = 0
      integer :: orbital_lmax = -1
      integer :: response_lmax = -1
      integer :: npoint = 0
      integer :: selector_kind = 0
      character(len=8) :: selector = ''
      character(len=40) :: core_policy = projected_core_policy
      real(rp) :: mesh_a = 0.0_rp
      real(rp) :: mesh_b = 0.0_rp
      real(rp) :: angular_scalar_integral = 0.0_rp
      real(rp), allocatable :: radius(:)
      real(rp), allocatable :: radial_weights(:)
      logical :: selected_l(0:2) = .false.
   contains
      procedure :: initialize => projected_site_spin_initialize
      procedure :: site_integration_functional => projected_site_spin_functional
      procedure :: selected_transition_coordinates => projected_selected_coordinates
      procedure :: transition_amplitudes => projected_transition_amplitudes
      procedure :: direct_operator_matrix => projected_direct_operator_matrix
      procedure :: direct_transition_amplitudes => projected_direct_transition_amplitudes
      procedure :: site_component_vertex_tensor => projected_site_component_vertex_tensor
      procedure :: moment_from_density => projected_moment_from_density
      procedure :: moment_from_operator => projected_moment_from_operator
      procedure :: core_spin_number => projected_core_spin_number
   end type projected_site_spin_contract

   ! --- from lr_ks_susceptibility_mod ---


   real(rp), parameter, public :: lr_kb_ry_per_kelvin = 6.3336814e-6_rp

   character(len=*), parameter, public :: lr_channel_plus = 'chi_plus'
   character(len=*), parameter, public :: lr_channel_minus = 'chi_minus'

   !> Immutable electronic-state snapshot consumed by LR-06.
   !>
   !> The initializer copies every array, including occupations and k weights.
   !> The susceptibility evaluator accepts this type with INTENT(IN), so it
   !> cannot recompute EF, occupations, or endpoint eigenvectors while summing.
   type, public :: lr_electronic_state
      integer :: nbands = 0
      integer :: nbasis = 0
      integer :: nk = 0
      real(rp), allocatable :: eigenvalues(:, :)       ! (band,k)
      complex(rp), allocatable :: eigenvectors(:, :, :) ! (basis,band,k)
      real(rp), allocatable :: k_points(:, :)          ! (3,k), folded points
      real(rp), allocatable :: k_weights(:)
      real(rp), allocatable :: occupations(:, :)       ! (band,k), explicit
      real(rp) :: fermi_level = 0.0_rp
      real(rp) :: temperature = 0.0_rp
      real(rp) :: energy_zero = 0.0_rp
      character(len=32) :: reciprocal_mode = 'ham_only'
      character(len=16) :: hamiltonian_order = 'second'
      logical :: orthogonal = .true.
      logical :: collinear = .true.
      logical :: has_soc = .false.
      logical :: has_extra_operator = .false.
      logical :: occupations_are_explicit = .false.
      logical :: initialized = .false.
   contains
      procedure :: initialize => lr_electronic_state_initialize
      procedure :: restore_to_default => lr_electronic_state_restore
      procedure :: validate => lr_electronic_state_validate
   end type lr_electronic_state

   !> Request metadata for one q/frequency/channel sweep.
   !>
   !> The pointer fields support the preferred request/result API.  An explicit
   !> overload is also provided for callers that keep these objects separately.
   type, public :: lr_ks_susceptibility_request
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = lr_channel_plus
      type(response_space_layout), pointer :: response_space => null()
      type(lmto_radial_basis), pointer :: radial_bases(:) => null()
      type(lr_electronic_state), pointer :: electronic_state => null()
      type(lr_electronic_state), pointer :: q_endpoint_state => null()
   end type lr_ks_susceptibility_request

   !> Result in the LR-04 canonical right-weighted representation.
   type, public :: lr_ks_susceptibility_result
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = ''
      character(len=128) :: response_representation = ''
      character(len=256) :: response_space_metadata = ''
      complex(rp), allocatable :: susceptibility(:, :, :) ! (I,J,frequency)
   end type lr_ks_susceptibility_result

   !> Request metadata for the compact weighted-orthonormal product response.
   !>
   !> This is deliberately a separate request type from LR-06.  A compact
   !> result must never be mistaken for an LR-04 point-space matrix whose
   !> dimension is `response_space%ndim`.
   type, public :: lr_product_ks_susceptibility_request
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = lr_channel_plus
      type(lmto_product_response_basis), pointer :: product_basis => null()
      type(lr_electronic_state), pointer :: electronic_state => null()
      type(lr_electronic_state), pointer :: q_endpoint_state => null()
      ! Optional DRESP-01 orbital selector.  An absent mask preserves the
      ! historical complete product-space oracle exactly.
      logical, allocatable :: selected_l(:)
   end type lr_product_ks_susceptibility_request

   !> Result in the weighted-orthonormal LMTO product representation.
   type, public :: lr_product_ks_susceptibility_result
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = ''
      character(len=128) :: response_representation = ''
      character(len=256) :: response_space_metadata = ''
      character(len=16) :: product_radial_order = 'second'
      integer :: product_endpoint_branches = lmto_product_nbranch
      character(len=32) :: product_branch_labels = '00,10,01,11,20,02'
      integer :: maximum_gf_energy_moment = lmto_product_max_gf_moment
      integer :: product_dimension = 0
      integer :: transition_dimension = 0
      integer :: ntransitions_evaluated = 0
      integer :: noccupation_skips = 0
      logical :: point_space_transition_allocated = .false.
      complex(rp), allocatable :: susceptibility(:, :, :) ! (product,product,frequency)
   end type lr_product_ks_susceptibility_result
   ! --- from lr_rs_gf_susceptibility_mod ---



   !> One directed real-space pair.  The coefficient-GF provider interprets
   !> gij as G_(left,right)(z), with the first index as the row/destination
   !> site and the second index as the column/source site.  gji is the
   !> separately evaluated reverse block G_(right,left)(z); it is not inferred
   !> by transposition.
   type, public :: lr_rs_gf_pair
      integer :: left_site = 0
      integer :: right_site = 0
      real(rp) :: translation(3) = 0.0_rp
   end type lr_rs_gf_pair

   !> Abstract native coefficient-GF provider contract.
   type, abstract, public :: lr_rs_gf_provider
      integer :: nbasis_per_site = 0
      character(len=32) :: provider_kind = ''
   contains
      procedure(lr_rs_get_pair), deferred :: get_pair
      procedure(lr_rs_describe), deferred :: describe
   end type lr_rs_gf_provider

   abstract interface
      subroutine lr_rs_get_pair(this, left_site, right_site, translation, z, gij, gji, hgamma_ij, hgamma_ji)
         import :: lr_rs_gf_provider, rp
         class(lr_rs_gf_provider), intent(inout) :: this
         integer, intent(in) :: left_site, right_site
         real(rp), intent(in) :: translation(3)
         complex(rp), intent(in) :: z
         complex(rp), intent(out) :: gij(:, :), gji(:, :), hgamma_ij(:, :), hgamma_ji(:, :)
      end subroutine lr_rs_get_pair

      function lr_rs_describe(this) result(description)
         import :: lr_rs_gf_provider
         class(lr_rs_gf_provider), intent(in) :: this
         character(len=256) :: description
      end function lr_rs_describe
   end interface

   !> Callback signature for a native block-recursion provider.  The callback
   !> is the seam to the production two-sweep/native block implementation.
   abstract interface
      subroutine lr_rs_native_callback(left_site, right_site, translation, z, gij, gji, hgamma_ij, hgamma_ji)
         import :: rp
         integer, intent(in) :: left_site, right_site
         real(rp), intent(in) :: translation(3)
         complex(rp), intent(in) :: z
         complex(rp), intent(out) :: gij(:, :), gji(:, :), hgamma_ij(:, :), hgamma_ji(:, :)
      end subroutine lr_rs_native_callback
   end interface

   !> Block-recursion native coefficient-GF provider.  It is kept distinct
   !> from the Chebyshev provider so recursion depth and terminator controls
   !> cannot be accidentally reported as polynomial controls.
   type, extends(lr_rs_gf_provider), public :: lr_rs_block_recursion_provider
      procedure(lr_rs_native_callback), pointer, nopass :: callback => null()
      integer :: recursion_depth = 0
      character(len=64) :: terminator = 'certified native block terminator'
   contains
      procedure :: initialize => lr_rs_block_recursion_initialize
      procedure :: get_pair => lr_rs_block_recursion_get_pair
      procedure :: describe => lr_rs_block_recursion_describe
   end type lr_rs_block_recursion_provider

   !> Chebyshev native coefficient-GF provider.  It has an independent
   !> polynomial-order control and the same directed-block contract.
   type, extends(lr_rs_gf_provider), public :: lr_rs_chebyshev_provider
      procedure(lr_rs_native_callback), pointer, nopass :: callback => null()
      integer :: polynomial_order = 0
      character(len=64) :: kernel = 'certified native Chebyshev kernel'
   contains
      procedure :: initialize => lr_rs_chebyshev_initialize
      procedure :: get_pair => lr_rs_chebyshev_get_pair
      procedure :: describe => lr_rs_chebyshev_describe
   end type lr_rs_chebyshev_provider

   !> Request for one real-space pair/Fourier/frequency sweep.
   type, public :: lr_rs_gf_susceptibility_request
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = lr_channel_plus
      integer :: integration_points = 2001
      real(rp) :: integration_eta = 0.0_rp
      real(rp) :: energy_margin = 1.0_rp
      real(rp) :: energy_min = -1.0_rp
      real(rp) :: energy_max = 1.0_rp
      real(rp) :: fermi_level = 0.0_rp
      real(rp) :: temperature = 0.0_rp
      type(response_space_layout), pointer :: response_space => null()
      type(lmto_radial_basis), pointer :: radial_bases(:) => null()
      class(lr_rs_gf_provider), pointer :: provider => null()
      type(lr_rs_gf_pair), allocatable :: pairs(:)
      ! Fractional/direct basis-site positions tau(:,site).  A null pointer
      ! means all endpoint positions are zero (the finite-cluster gauge).
      real(rp), pointer :: site_positions(:, :) => null()
   end type lr_rs_gf_susceptibility_request

   ! --- from tddft_native_rsgf_provider_mod ---


   type, extends(lr_rs_gf_provider), public :: tddft_native_rsgf_provider
      type(green), pointer :: green_obj => null()
      type(recursion), pointer :: recursion_obj => null()
      type(hamiltonian), pointer :: hamiltonian_obj => null()
      type(lattice), pointer :: lattice_obj => null()
      type(reciprocal), pointer :: reciprocal_obj => null()
      type(lmto_radial_basis), pointer :: radial_bases(:) => null()
      character(len=32) :: selection = ''
      integer :: recursion_depth = 0
      integer :: polynomial_order = 0
      type(lr_rs_gf_pair), allocatable :: pairs(:)
      integer, allocatable :: left_atoms(:), right_atoms(:)
      complex(rp), allocatable :: hgamma_ij(:, :, :), hgamma_ji(:, :, :)
      real(rp), allocatable :: a_inf(:, :, :, :), b_inf(:, :, :, :)
      complex(rp), allocatable :: psi(:, :, :), psi_out(:, :, :), identity(:, :)
   contains
      procedure :: initialize => tddft_native_rsgf_initialize
      procedure :: get_pair => tddft_native_rsgf_get_pair
      procedure :: describe => tddft_native_rsgf_describe
   end type tddft_native_rsgf_provider

   ! --- from lr_projected_reciprocal_chi0_mod ---



   character(len=*), parameter, public :: projected_backend_lehmann = 'lehmann'
   character(len=*), parameter, public :: projected_backend_gf = 'real-axis-gf'
   character(len=*), parameter, public :: projected_backend_finite_width = 'finite-width'


   !> Compact record for the largest direct Lehmann terms retained by the
   !> optional DRESP-02R audit.  The response equation is unchanged; this is
   !> an observable of the already certified transition sum.
   type, public :: projected_chi0_transition_record
      integer :: k_index = 0
      integer :: left_band = 0
      integer :: right_band = 0
      real(rp) :: left_energy = 0.0_rp
      real(rp) :: right_energy = 0.0_rp
      real(rp) :: left_occupation = 0.0_rp
      real(rp) :: right_occupation = 0.0_rp
      real(rp) :: transition_energy = 0.0_rp
      real(rp) :: matrix_element_weight = 0.0_rp
      real(rp) :: score = 0.0_rp
   end type projected_chi0_transition_record

   !> One projected reciprocal response request.  The DRESP-01 contract and
   !> product basis are separate references so the same certified operator can
   !> be compared with both compact product-space oracles.
   type, public :: projected_chi0_request
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = lr_channel_plus
      integer :: integration_points = 2001
      real(rp) :: integration_eta = 0.0_rp
      real(rp) :: energy_margin = 1.0_rp
      type(projected_site_spin_contract), pointer :: contract => null()
      type(lmto_product_response_basis), pointer :: product_basis => null()
      type(lr_electronic_state), pointer :: electronic_state => null()
      type(lr_electronic_state), pointer :: q_endpoint_state => null()
      ! DRESP-02R diagnostics are opt-in and do not alter the production
      ! response accumulation or its state handoff.
      logical :: diagnostics = .false.
   end type projected_chi0_request

   type, public :: projected_chi0_result
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      real(rp) :: integration_eta = 0.0_rp
      real(rp) :: energy_min = 0.0_rp
      real(rp) :: energy_max = 0.0_rp
      real(rp) :: energy_spacing = 0.0_rp
      real(rp) :: spacing_over_integration_eta = 0.0_rp
      character(len=32) :: channel = ''
      character(len=16) :: backend = ''
      character(len=8) :: selector = ''
      character(len=128) :: response_representation = ''
      character(len=512) :: provenance = ''
      integer :: nsite = 0
      integer :: product_dimension = 0
      integer :: integration_points = 0
      integer :: ntransitions_evaluated = 0
      integer :: noccupation_skips = 0
      logical :: point_response_allocated = .false.
      integer(int64) :: site_vertex_memory_bytes = 0_int64
      integer(int64) :: susceptibility_memory_bytes = 0_int64
      real(rp) :: wall_time_seconds = 0.0_rp
      real(rp) :: cpu_time_seconds = 0.0_rp
      ! Performance observables.  The reference backend reports the dense
      ! resolvent and scalar bubble counts; the optimized backend reports the
      ! eigenbasis transform and vectorized denominator work instead.
      character(len=16) :: implementation = ''
      real(rp) :: allocation_seconds = 0.0_rp
      real(rp) :: vertex_seconds = 0.0_rp
      real(rp) :: endpoint_transform_seconds = 0.0_rp
      real(rp) :: resolvent_seconds = 0.0_rp
      real(rp) :: accumulator_seconds = 0.0_rp
      real(rp) :: diagnostic_seconds = 0.0_rp
      integer(int64) :: resolvent_calls = 0_int64
      integer(int64) :: accumulator_calls = 0_int64
      integer(int64) :: endpoint_transform_calls = 0_int64
      complex(rp), allocatable :: susceptibility(:, :, :) ! (site,site,frequency)
      ! Optional DRESP-02R Kubo and k-resolved observables.
      complex(rp), allocatable :: kubo_term_one(:, :, :) ! (site,site,frequency)
      complex(rp), allocatable :: kubo_term_two(:, :, :) ! (site,site,frequency)
      complex(rp), allocatable :: k_susceptibility(:, :, :, :) ! (site,site,frequency,k)
      real(rp) :: left_spectral_zeroth_residual = 0.0_rp
      real(rp) :: right_spectral_zeroth_residual = 0.0_rp
      real(rp) :: left_spectral_first_residual = 0.0_rp
      real(rp) :: right_spectral_first_residual = 0.0_rp
      real(rp) :: left_spectral_fermi_residual = 0.0_rp
      real(rp) :: right_spectral_fermi_residual = 0.0_rp
      integer :: n_dominant_transition_records = 0
      type(projected_chi0_transition_record), allocatable :: dominant_transitions(:)
   end type projected_chi0_result

   ! --- from lr_alsda_kernel_mod (LR-REF-03c) ---

   character(len=*), parameter, public :: lr_kxc_magnetization_pauli = 'pauli_projected'
   character(len=*), parameter, public :: lr_kxc_magnetization_kind_pauli_accepted = 'PAULI_ACCEPTED'
   character(len=*), parameter, public :: lr_kxc_magnetization_source_pauli_accepted = &
      'accepted reciprocal eigensystem + POTPAR large component + frozen core'
   character(len=*), parameter, public :: lr_kxc_units = 'Ry bohr^3'
   character(len=*), parameter, public :: lr_kxc_representation = 'LR-04 canonical local operator'
   real(rp), parameter, public :: lr_kxc_default_low_m_relative = 1.0e-8_rp

   !> Request for one direct ALSDA kernel on a complete response layout.
   !>
   !> `pauli_magnetization(site,radial_point)` is the Pauli/no-SOC response
   !> magnetization selected by LR-03.  It is not reconstructed from a
   !> compressed potential field and it is not silently replaced by the SR
   !> density in the accepted radial snapshot.
   type, public :: lr_alsda_kernel_request
      type(response_space_layout), pointer :: response_space => null()
      type(radial_ground_state), pointer :: ground_states(:) => null()
      real(rp), allocatable :: pauli_magnetization(:, :)
      character(len=32) :: magnetization_label = lr_kxc_magnetization_pauli
      character(len=32) :: magnetization_kind = ''
      character(len=256) :: magnetization_source = ''
      logical :: production_contract = .false.
      character(len=512) :: requested_functional = ''
      character(len=32) :: requested_backend = ''
      integer :: requested_txc = -1
      ! This is a diagnostic threshold only.  It never changes the ratio.
      real(rp) :: low_m_diagnostic_relative = lr_kxc_default_low_m_relative
   end type lr_alsda_kernel_request

   !> Low-m diagnostic for the active positive-measure radial grid.
   type, public :: lr_alsda_low_m_diagnostic
      integer :: active_radial_points = 0
      integer :: low_m_points = 0
      integer :: exact_zero_active_points = 0
      integer :: null_measure_zero_points = 0
      real(rp) :: max_abs_magnetization = 0.0_rp
      real(rp) :: min_abs_magnetization = huge(1.0_rp)
      real(rp) :: min_relative_abs_magnetization = huge(1.0_rp)
      real(rp) :: diagnostic_relative_threshold = lr_kxc_default_low_m_relative
   end type lr_alsda_low_m_diagnostic

   !> Direct ALSDA result and provenance carried to later interaction code.
   type, public :: lr_alsda_kernel_result
      real(rp), allocatable :: pointwise_kernel(:, :) ! (site,radial point)
      complex(rp), allocatable :: canonical_operator(:, :) ! LR-04 B=K_raw*W
      type(lr_alsda_low_m_diagnostic) :: low_m
      type(radial_xc_provenance) :: xc_provenance
      character(len=32) :: magnetization_label = ''
      character(len=32) :: magnetization_kind = ''
      character(len=256) :: magnetization_source = ''
      character(len=64) :: units = lr_kxc_units
      character(len=128) :: response_representation = lr_kxc_representation
      logical :: origin_null_measure_extension = .false.
   end type lr_alsda_kernel_result

   ! --- from lr_goldstone_sumrule_mod (LR-REF-03c) ---

   character(len=*), parameter, public :: lr_gsr_magnetization_pauli = 'pauli_projected'
   character(len=*), parameter, public :: lr_gsr_route = 'goldstone_sumrule'
   character(len=*), parameter, public :: lr_gsr_units = 'Ry bohr^3'
   character(len=*), parameter, public :: lr_gsr_representation = &
      'LR-04 canonical interaction K_eff=4*pi*U_LCMM'
   real(rp), parameter, public :: lr_gsr_residual_tolerance = 1.0e-9_rp

   !> Request for the independent static sum-rule construction.
   type, public :: lr_goldstone_sumrule_request
      type(response_space_layout), pointer :: response_space => null()
      complex(rp), allocatable :: static_susceptibility(:, :) ! canonical LR-04 chiKS(0)
      real(rp), allocatable :: magnetization(:, :) ! physical m_z(site,r)
      character(len=32) :: magnetization_label = lr_gsr_magnetization_pauli
   end type lr_goldstone_sumrule_request

   !> Diagnostics and products of the sum-rule construction.
   type, public :: lr_goldstone_sumrule_result
      ! U_LCMM(site,r)=B_eff(site,r)/(4*pi*m_z(site,r)); origin is zero.
      complex(rp), allocatable :: u_lcmm(:, :)
      ! Pointwise LCMM energy splitting B_eff=4*pi*U_LCMM*m_z.
      complex(rp), allocatable :: effective_field(:, :)
      ! Canonical LR-04 local operator for B_eff,00=K_eff*m_00.
      complex(rp), allocatable :: canonical_interaction(:, :)
      ! Normalized response vector m_00=sqrt(4*pi)*m_z, origin excluded.
      complex(rp), allocatable :: magnetization_response(:)
      ! Full response-space B_eff,00 generated by K_eff acting on m_00.
      complex(rp), allocatable :: generated_field(:)
      complex(rp), allocatable :: generated_response(:)
      complex(rp), allocatable :: residual_vector(:)
      ! Gamma is the active spherical radial equation in the solve order.
      complex(rp), allocatable :: gamma(:, :)
      complex(rp), allocatable :: solution_active(:)
      real(rp), allocatable :: singular_values(:)
      integer :: rank = 0
      integer :: unknowns = 0
      real(rp) :: condition_number = huge(1.0_rp)
      real(rp) :: equation_residual_norm = huge(1.0_rp)
      real(rp) :: equation_relative_residual = huge(1.0_rp)
      real(rp) :: metric_residual_norm = huge(1.0_rp)
      real(rp) :: metric_relative_residual = huge(1.0_rp)
      complex(rp) :: rigid_overlap = cmplx(0.0_rp, 0.0_rp, rp)
      logical :: rank_deficient = .false.
      logical :: blocked = .true.
      character(len=32) :: route = lr_gsr_route
      character(len=32) :: magnetization_label = ''
      character(len=64) :: units = lr_gsr_units
      character(len=128) :: response_representation = lr_gsr_representation
      character(len=256) :: status = 'not evaluated'
   end type lr_goldstone_sumrule_result

   !> Result of the generic complex SVD solve Gamma*x=rhs.
   type, public :: lr_sumrule_linear_solve_result
      complex(rp), allocatable :: solution(:)
      real(rp), allocatable :: singular_values(:)
      integer :: rank = 0
      real(rp) :: condition_number = huge(1.0_rp)
      real(rp) :: residual_norm = huge(1.0_rp)
      real(rp) :: relative_residual = huge(1.0_rp)
      logical :: rank_deficient = .true.
      logical :: blocked = .true.
      character(len=256) :: status = 'not evaluated'
   end type lr_sumrule_linear_solve_result

   ! --- from lr_projected_interacting_response_mod (LR-REF-03c) ---

   character(len=*), parameter, public :: projected_mills_exact_scalar = 'EXACT_SCALAR'
   character(len=*), parameter, public :: projected_mills_projected_scalar = &
      'PROJECTED_SCALAR_APPROXIMATION'
   character(len=*), parameter, public :: projected_mills_unsupported = 'UNSUPPORTED'
   character(len=*), parameter, public :: projected_mills_convention = &
      'H=H0 I+B_sigma sigma_z; H_up-H_down=2 B_sigma; DRESP-01 Vz=sigma_z'
   character(len=*), parameter, public :: projected_mills_loss_convention = &
      'L=-(chi-chi^dagger)/(2*i*pi); ordinary site-space matrix'

   real(rp), parameter, public :: projected_mills_exact_tolerance = 1.0e-10_rp
   real(rp), parameter, public :: projected_mills_moment_floor = 1.0e-12_rp
   real(rp), parameter, public :: projected_mills_condition_limit = 1.0e12_rp

   !> Accepted-state samples for the local scalarization.  The matrices are
   !> site-major in the supplied coefficient subspace.  A sample is normally
   !> one accepted reciprocal H(k); sample_weights are relative weights and
   !> default to a uniform Frobenius metric.
   type, public :: projected_mills_interaction_request
      character(len=8) :: selector = ''
      character(len=80) :: splitting_convention = projected_mills_convention
      integer :: nsite = 0
      integer :: site_block_size = 0
      complex(rp), allocatable :: actual_pauli_field(:, :, :) ! (basis,basis,sample)
      complex(rp), allocatable :: site_vertices(:, :, :, :) ! (basis,basis,site,sample)
      real(rp), allocatable :: projected_moment(:)
      real(rp), allocatable :: sample_weights(:)
      character(len=512) :: provenance = ''
   end type projected_mills_interaction_request

   type, public :: projected_mills_interaction_result
      character(len=8) :: selector = ''
      character(len=80) :: splitting_convention = projected_mills_convention
      character(len=32) :: classification = projected_mills_unsupported
      character(len=512) :: provenance = ''
      integer :: nsite = 0
      integer :: nbasis = 0
      integer :: nsample = 0
      integer :: rank = 0
      real(rp) :: fit_condition_number = huge(1.0_rp)
      real(rp) :: scalarization_residual = huge(1.0_rp)
      real(rp) :: locality_residual = huge(1.0_rp)
      real(rp) :: splitting_norm = 0.0_rp
      real(rp) :: fit_residual_norm = 0.0_rp
      real(rp) :: coefficient_imaginary_residual = 0.0_rp
      real(rp), allocatable :: projected_moment(:)
      real(rp), allocatable :: projected_splitting(:) ! B_sigma coefficient Delta_i
      real(rp), allocatable :: interaction_U(:)
   end type projected_mills_interaction_result

   !> Site-space projected Dyson request.  bare_chi is directly the
   !> DRESP-02 site x site x frequency object; no product-space input exists.
   type, public :: projected_dyson_request
      character(len=8) :: selector = ''
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = ''
      real(rp), allocatable :: interaction_U(:)
      complex(rp), allocatable :: bare_chi(:, :, :) ! (site,site,frequency)
      character(len=512) :: interaction_provenance = ''
      character(len=512) :: bare_provenance = ''
   end type projected_dyson_request

   type, public :: projected_dyson_result
      character(len=8) :: selector = ''
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = ''
      character(len=512) :: interaction_provenance = ''
      character(len=512) :: bare_provenance = ''
      character(len=256) :: dyson_convention = 'D=I-chi0*U; solve D*chi=chi0 with certified LAPACK zgesv'
      character(len=256) :: loss_convention = projected_mills_loss_convention
      character(len=128) :: status = 'not evaluated'
      real(rp), allocatable :: interaction_U(:)
      complex(rp), allocatable :: bare_chi(:, :, :)
      complex(rp), allocatable :: interaction(:, :)
      complex(rp), allocatable :: denominator(:, :, :)
      complex(rp), allocatable :: enhanced_chi(:, :, :)
      complex(rp), allocatable :: loss_matrix(:, :, :)
      real(rp), allocatable :: denominator_min_singular_value(:)
      real(rp), allocatable :: denominator_max_singular_value(:)
      real(rp), allocatable :: condition_number(:)
      real(rp), allocatable :: minimum_magnitude_eigenvalue(:)
      real(rp), allocatable :: loss_trace(:)
      real(rp), allocatable :: minus_im_trace_over_pi(:)
      real(rp), allocatable :: dyson_residual(:)
      real(rp), allocatable :: dyson_residual_relative(:)
      real(rp), allocatable :: dyson_residual_infinity(:)
      integer, allocatable :: solve_info(:)
   end type projected_dyson_result

   ! --- from lr_projected_juelich_interaction_mod (LR-REF-03c) ---

   character(len=*), parameter, public :: projected_juelich_exact_local = 'EXACT_LOCAL_SUMRULE'
   character(len=*), parameter, public :: projected_juelich_projected_local = 'PROJECTED_LOCAL_SUMRULE'
   character(len=*), parameter, public :: projected_juelich_eta_limited = 'ETA_LIMITED_LOCAL_SUMRULE'
   character(len=*), parameter, public :: projected_juelich_rank_deficient = 'UNSUPPORTED_RANK_DEFICIENT'
   character(len=*), parameter, public :: projected_juelich_unsupported = 'UNSUPPORTED'

   real(rp), parameter, public :: projected_juelich_residual_tolerance = 2.0e-9_rp
   real(rp), parameter, public :: projected_juelich_condition_limit = 1.0e12_rp
   real(rp), parameter, public :: projected_juelich_eta_relative_tolerance = 5.0e-3_rp

   type, public :: projected_juelich_request
      character(len=8) :: selector = ''
      real(rp), allocatable :: projected_moment(:)
      complex(rp), allocatable :: static_chi0(:, :)
      real(rp) :: static_eta = 0.0_rp
      real(rp) :: q(3) = 0.0_rp
      character(len=32) :: channel = ''
      character(len=512) :: state_provenance = ''
      character(len=512) :: chi0_provenance = ''
   end type projected_juelich_request

   type, public :: projected_juelich_result
      character(len=8) :: selector = ''
      character(len=32) :: classification = projected_juelich_unsupported
      character(len=512) :: provenance = ''
      real(rp) :: static_eta = 0.0_rp
      real(rp) :: q(3) = 0.0_rp
      character(len=32) :: channel = ''
      integer :: nsite = 0
      integer :: rank = 0
      integer :: real_rank = 0
      real(rp) :: condition_number = huge(1.0_rp)
      real(rp) :: real_condition_number = huge(1.0_rp)
      real(rp) :: equation_residual = huge(1.0_rp)
      real(rp) :: relative_residual = huge(1.0_rp)
      real(rp) :: real_constrained_residual = huge(1.0_rp)
      real(rp) :: relative_real_constrained_residual = huge(1.0_rp)
      real(rp) :: imaginary_U_ratio = huge(1.0_rp)
      real(rp) :: holdout_residual = huge(1.0_rp)
      real(rp) :: holdout_relative_residual = huge(1.0_rp)
      logical :: eta_stable = .false.
      character(len=256) :: eta_stability = 'not evaluated'
      real(rp), allocatable :: projected_moment(:)
      real(rp), allocatable :: singular_values(:)
      real(rp), allocatable :: real_singular_values(:)
      complex(rp), allocatable :: gamma(:, :)
      complex(rp), allocatable :: interaction_matrix(:, :)
      complex(rp), allocatable :: interaction_U_complex(:)
      real(rp), allocatable :: interaction_U_real(:)
      real(rp), allocatable :: interaction_U(:)
   end type projected_juelich_result

   ! --- from lr_compact_static_interaction_mod (LR-REF-03c) ---

   real(rp), parameter, public :: lr_compact_gsr_residual_tolerance = 1.0e-9_rp
   real(rp), parameter, public :: lr_compact_gsr_action_tolerance = 1.0e-10_rp
   character(len=*), parameter, public :: lr_compact_representation = &
      'weighted-orthonormal LMTO product representation; U=sqrt(W)V'
   character(len=*), parameter, public :: lr_compact_mapping_contract = &
      'c=U^H sqrt(W) x, x=inv(sqrt(W)) U c, Kc=U^H K U'

   !> Independent compact LCMM/GSR result.  The unknown is a local scalar
   !> U_LCMM represented in the retained L=0 radial product span.  The solve
   !> equation is rows-by-unknowns, Gamma*u=target, and is deliberately
   !> allowed to return BLOCKED when the full compact response cannot satisfy
   !> the spherical local sum-rule ansatz.
   type, public :: lr_compact_gsr_result
      complex(rp), allocatable :: u_lcmm(:, :)       ! site, radial point
      complex(rp), allocatable :: effective_field(:, :) ! 4*pi*U*m_z
      complex(rp), allocatable :: canonical_interaction(:, :) ! compact local K_eff
      complex(rp), allocatable :: magnetization_response(:)
      complex(rp), allocatable :: generated_field(:)
      complex(rp), allocatable :: generated_response(:)
      complex(rp), allocatable :: residual_vector(:)
      complex(rp), allocatable :: gamma(:, :)
      complex(rp), allocatable :: solution(:)
      real(rp), allocatable :: singular_values(:)
      integer :: equation_rows = 0
      integer :: unknowns = 0
      integer :: rank = 0
      real(rp) :: condition_number = huge(1.0_rp)
      real(rp) :: svd_rcond = -1.0_rp
      real(rp) :: svd_cutoff = huge(1.0_rp)
      real(rp) :: coefficient_norm = huge(1.0_rp)
      real(rp) :: equation_residual_norm = huge(1.0_rp)
      real(rp) :: equation_relative_residual = huge(1.0_rp)
      real(rp) :: residual_norm = huge(1.0_rp)
      real(rp) :: relative_residual = huge(1.0_rp)
      real(rp) :: residual_difference_norm = huge(1.0_rp)
      real(rp) :: residual_difference_relative = huge(1.0_rp)
      real(rp) :: residual_difference_max_component = huge(1.0_rp)
      real(rp) :: max_residual_component = huge(1.0_rp)
      real(rp) :: assembled_action_sum_norm = huge(1.0_rp)
      real(rp) :: assembled_action_matrix_norm = huge(1.0_rp)
      real(rp) :: assembled_action_difference_norm = huge(1.0_rp)
      real(rp) :: assembled_action_difference_relative = huge(1.0_rp)
      real(rp) :: assembled_action_difference_max_component = huge(1.0_rp)
      real(rp) :: magnetization_weighted_norm = huge(1.0_rp)
      real(rp) :: magnetization_projection_residual_norm = huge(1.0_rp)
      real(rp) :: magnetization_projection_relative_residual = huge(1.0_rp)
      complex(rp) :: rigid_overlap = cmplx(0.0_rp, 0.0_rp, rp)
      logical :: rank_deficient = .true.
      logical :: blocked = .true.
      character(len=64) :: status = 'not evaluated'
   end type lr_compact_gsr_result

   ! --- from tddft_dyson_mod (LR-REF-03c) ---

   character(len=*), parameter, public :: lr_dyson_route_direct_alsda = 'direct_alsda'
   character(len=*), parameter, public :: lr_dyson_route_goldstone_sumrule = 'goldstone_sumrule'
   character(len=*), parameter, public :: lr_dyson_route_direct_alsda_goldstone_corrected = &
      'direct_alsda_goldstone_corrected'
   character(len=*), parameter, public :: lr_dyson_response_representation = &
      'LR-04 canonical right-weighted B=A*W'
   character(len=*), parameter, public :: lr_dyson_compact_response_representation = &
      'orthonormal compact LMTO product representation'
   character(len=*), parameter, public :: lr_dyson_loss_convention = &
      'L=-(chi-chi^dagger_W)/(2*i*pi); canonical metric-adjoint form'
   character(len=*), parameter, public :: lr_dyson_convention = &
      'solve (I-chi_KS*K) chi=chi_KS with LAPACK zgesv; no explicit inverse'

   real(rp), parameter, public :: lr_dyson_near_singular_condition = 1.0e8_rp
   real(rp), parameter, public :: lr_dyson_ill_conditioned_condition = 1.0e12_rp

   !> Request for one q/frequency/channel response batch.
   !>
   !> `canonical_interaction` must already be in the LR-04 representation.  In
   !> particular, this request has no raw-kernel escape hatch and no automatic
   !> route selection.  For the corrected route the supplied matrix must be
   !> the explicitly selected GCR-01 corrected Kxc.
   type, public :: tddft_dyson_request
      type(response_space_layout), pointer :: response_space => null()
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = ''
      ! When true, the matrices are already in the accepted orthonormal
      ! compact product representation.  The response_space pointer is then
      ! optional and no point-space metric is inserted into Dyson or loss.
      logical :: compact_orthonormal = .false.
      complex(rp), allocatable :: ks_susceptibility(:, :, :)
      complex(rp), allocatable :: canonical_interaction(:, :)
      character(len=48) :: interaction_route = lr_dyson_route_direct_alsda
      character(len=256) :: interaction_provenance = ''
      character(len=512) :: electronic_state_provenance = ''
      character(len=256) :: response_space_metadata = ''
   end type tddft_dyson_request

   !> Enhanced response, loss matrices, and per-frequency diagnostics.
   type, public :: tddft_dyson_result
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = ''
      logical :: compact_orthonormal = .false.
      character(len=48) :: interaction_route = ''
      character(len=256) :: interaction_provenance = ''
      character(len=512) :: electronic_state_provenance = ''
      character(len=256) :: response_space_metadata = ''
      character(len=128) :: response_representation = lr_dyson_response_representation
      character(len=256) :: dyson_convention = lr_dyson_convention
      character(len=256) :: loss_convention = lr_dyson_loss_convention
      character(len=256) :: status = 'not evaluated'
      complex(rp), allocatable :: ks_susceptibility(:, :, :)
      complex(rp), allocatable :: canonical_interaction(:, :)
      complex(rp), allocatable :: denominator(:, :, :)
      complex(rp), allocatable :: enhanced_susceptibility(:, :, :)
      complex(rp), allocatable :: loss_matrix(:, :, :)
      real(rp), allocatable :: denominator_min_singular_value(:)
      real(rp), allocatable :: denominator_max_singular_value(:)
      real(rp), allocatable :: denominator_condition_number(:)
      real(rp), allocatable :: denominator_min_magnitude_eigenvalue(:)
      integer, allocatable :: solve_info(:)
      logical, allocatable :: solve_succeeded(:)
      logical, allocatable :: near_singular_collective_pole(:)
      logical, allocatable :: numerically_singular(:)
      logical, allocatable :: ill_conditioned(:)
      real(rp), allocatable :: loss_metric_hermiticity_residual(:)
      real(rp), allocatable :: dyson_residual_frobenius(:)
      real(rp), allocatable :: dyson_residual_relative(:)
      real(rp), allocatable :: dyson_residual_infinity(:)
      character(len=256), allocatable :: frequency_status(:)
   end type tddft_dyson_result

   ! --- kernel/Dyson public procedures moved under linear_response_mod ---
   public :: evaluate_lr_alsda_kernel
   public :: evaluate_lr_alsda_static_residual
   public :: lr_alsda_provenance_matches
   public :: evaluate_lr_goldstone_sumrule
   public :: solve_lr_goldstone_equation
   public :: evaluate_projected_mills_interaction
   public :: evaluate_projected_mills_from_reciprocal
   public :: evaluate_projected_dyson
   public :: evaluate_projected_juelich_interaction
   public :: evaluate_projected_juelich_holdout
   public :: assess_projected_juelich_eta_stability
   public :: select_projected_juelich_eta_indices
   public :: compact_project_point_vector
   public :: compact_reconstruct_point_vector
   public :: compact_project_local_operator
   public :: compact_apply_local_operator
   public :: compact_project_magnetization
   public :: compact_weighted_projection_diagnostics
   public :: evaluate_compact_goldstone_sumrule
   public :: evaluate_tddft_dyson
   public :: solve_tddft_dyson_frequency
   public :: tddft_denominator_minimum_magnitude_eigenvalue
   public :: tddft_loss_matrix
   public :: loss_matrix_hermiticity_residual

   ! --- bare-response public procedures ---
   public :: lr_fermi_dirac_occupation
   public :: lr_snapshot_from_reciprocal
   public :: lr_q_endpoint_from_reciprocal
   public :: evaluate_lr_ks_susceptibility
   public :: evaluate_lr_product_ks_susceptibility
   public :: evaluate_lr_static_residual
   public :: evaluate_lr_rs_gf_susceptibility
   public :: evaluate_projected_lehmann_chi0
   public :: evaluate_projected_finite_width_chi0
   public :: evaluate_projected_gf_chi0
   public :: evaluate_projected_gf_chi0_optimized
   public :: evaluate_projected_gf_chi0_reference

   public :: response_lmax
   public :: response_product_space_complete
   public :: response_lm_index
   public :: response_lm_from_index
   public :: response_harmonic
   public :: response_gaunt
   public :: response_product_coefficients
   public :: response_superindex_size
   public :: response_flatten_superindex
   public :: response_unflatten_superindex
   public :: response_simpson_weight
   public :: response_log_mesh_jacobian
   public :: response_volume_measure
   public :: response_weighted_density
   public :: response_physical_density
   public :: response_log_mesh_integral
   public :: response_endpoint_phase
   public :: response_real_space_phase
   public :: response_site_gauge
   public :: response_apply_site_gauge
   public :: response_mapping_supported
   public :: response_build_radial_metric
   public :: response_metric_weights
   public :: response_vector_inner_product
   public :: response_vector_norm
   public :: response_apply_operator
   public :: response_compose_operators
   public :: response_identity_operator
   public :: response_local_operator
   public :: response_operator_adjoint
   public :: response_operator_trace
   public :: response_rigid_vector_norm
   public :: response_rigid_vector_overlap
   public :: response_raw_to_canonical
   public :: response_canonical_to_raw
   public :: lmto_product_branch_powers
   public :: lmto_product_branch_label
   public :: lmto_product_branch_valid
   public :: lmto_product_energy_power
   public :: lmto_product_second_order_radial_branch
   public :: lmto_product_apply_branch_action
   public :: pauli_charge_matrix
   public :: pauli_sigma_x_matrix
   public :: pauli_sigma_y_matrix
   public :: pauli_sigma_z_matrix
   public :: pauli_sigma_plus_matrix
   public :: pauli_sigma_minus_matrix
   public :: evaluate_pauli_transition_vertex
   public :: lmto_product_candidate_count
   public :: lmto_enumerate_product_candidates
   public :: augment_lr_gf_endpoint
   public :: augment_lr_gf_endpoint_pair
   public :: projected_selector_supported
   public :: lmto_fixture_init
   public :: lmto_fixture_from_hamiltonian
   public :: assemble_lmto_hamiltonian
   public :: assemble_lmto_finite_q_torque
   public :: assemble_lmto_finite_q_torques
   public :: assemble_lmto_finite_q_mixed_derivative
   public :: force_theorem_finite_q_hessian_from_eigenbasis
   public :: force_theorem_finite_q_hessian_from_eigenbasis_batch
   public :: force_theorem_finite_q_hessian_from_eigenbasis_metallic
   public :: force_theorem_finite_q_hessian_from_eigenbasis_metallic_batch
   public :: finite_temperature_occupation
   public :: fermi_divided_difference
   public :: lmto_fixture_adapter_residual
   public :: build_finite_temperature_contour
   public :: build_zero_temperature_occupied_contour
   public :: finite_temperature_complex_fermi
   public :: finite_temperature_regularized_fermi
   public :: force_theorem_finite_q_hessian_from_resolvent
   public :: force_theorem_finite_q_hessian_from_resolvent_batch
   public :: force_theorem_finite_q_hessian_from_zero_temperature_contour
   public :: prepare_rotation_response
   public :: evaluate_rotation_response
   public :: evaluate_rotation_response_oracle
   public :: reduce_static_rotation_kernel
   public :: rotation_axes
   public :: rotation_circular_unitary
   public :: compute_static_rotation_curvature
   public :: detect_commensurate_endpoint_map
   public :: validate_rotation_capability
   public :: validate_rotation_pole_workflow_capability




   interface evaluate_lr_alsda_kernel
      module procedure evaluate_lr_alsda_kernel_request
      module procedure evaluate_lr_alsda_kernel_explicit
   end interface evaluate_lr_alsda_kernel

   interface evaluate_lr_goldstone_sumrule
      module procedure evaluate_lr_goldstone_sumrule_request
      module procedure evaluate_lr_goldstone_sumrule_explicit
   end interface evaluate_lr_goldstone_sumrule

   interface compact_project_local_operator
      module procedure compact_project_local_operator_complex
   end interface compact_project_local_operator

   interface tddft_loss_matrix
      module procedure tddft_loss_matrix_ordinary
      module procedure tddft_loss_matrix_metric
   end interface tddft_loss_matrix

   interface loss_matrix_hermiticity_residual
      module procedure loss_matrix_hermiticity_residual_ordinary
      module procedure loss_matrix_hermiticity_residual_metric
   end interface loss_matrix_hermiticity_residual

   interface response_local_operator
      module procedure response_local_operator_site_radial
      module procedure response_local_operator_site_radial_channel
      module procedure response_local_operator_complex_site_radial
      module procedure response_local_operator_complex_site_radial_channel
   end interface response_local_operator

   interface evaluate_pauli_transition_vertex
      module procedure evaluate_pauli_transition_vertex_one
      module procedure evaluate_pauli_transition_vertex_channels
   end interface evaluate_pauli_transition_vertex

   interface evaluate_lr_ks_susceptibility
      module procedure evaluate_lr_ks_susceptibility_request
      module procedure evaluate_lr_ks_susceptibility_explicit
   end interface evaluate_lr_ks_susceptibility

   interface evaluate_lr_product_ks_susceptibility
      module procedure evaluate_lr_product_ks_susceptibility_request
      module procedure evaluate_lr_product_ks_susceptibility_explicit
   end interface evaluate_lr_product_ks_susceptibility
   interface
      module function tddft_capability_is_supported(state, reason) result(ok)
         type(tddft_capability_state), intent(in) :: state
         character(len=*), intent(out) :: reason
         logical :: ok
      end function tddft_capability_is_supported

      module subroutine require_tddft_capability(state)
         type(tddft_capability_state), intent(in) :: state
      end subroutine require_tddft_capability

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
      end subroutine evaluate_linear_response_sweep

      ! --- kernel/Dyson submodule procedures ---
      module subroutine evaluate_lr_alsda_kernel_request(request, result)
         type(lr_alsda_kernel_request), intent(in) :: request
         type(lr_alsda_kernel_result), intent(out) :: result
      end subroutine evaluate_lr_alsda_kernel_request

      module subroutine evaluate_lr_alsda_kernel_explicit(space, ground_states, pauli_magnetization, result, &
         requested_functional, requested_backend, requested_txc, magnetization_label, low_m_diagnostic_relative)
         type(response_space_layout), intent(in) :: space
         type(radial_ground_state), intent(in) :: ground_states(:)
         real(rp), intent(in) :: pauli_magnetization(:, :)
         type(lr_alsda_kernel_result), intent(out) :: result
         character(len=*), intent(in), optional :: requested_functional
         character(len=*), intent(in), optional :: requested_backend
         integer, intent(in), optional :: requested_txc
         character(len=*), intent(in), optional :: magnetization_label
         real(rp), intent(in), optional :: low_m_diagnostic_relative
      end subroutine evaluate_lr_alsda_kernel_explicit

      module subroutine evaluate_lr_alsda_static_residual(space, static_susceptibility, canonical_kernel, rigid_vector, &
         absolute_residual, relative_residual, rigid_overlap, residual_vector)
         type(response_space_layout), intent(in) :: space
         complex(rp), intent(in) :: static_susceptibility(:, :), canonical_kernel(:, :), rigid_vector(:)
         real(rp), intent(out) :: absolute_residual, relative_residual
         complex(rp), intent(out), optional :: rigid_overlap
         complex(rp), intent(out), optional :: residual_vector(:)
      end subroutine evaluate_lr_alsda_static_residual

      module logical function lr_alsda_provenance_matches(accepted, requested_functional, requested_backend, requested_txc) &
         result(matches)
         type(radial_xc_provenance), intent(in) :: accepted
         character(len=*), intent(in) :: requested_functional, requested_backend
         integer, intent(in) :: requested_txc
      end function lr_alsda_provenance_matches

      module subroutine evaluate_lr_goldstone_sumrule_request(request, result)
         type(lr_goldstone_sumrule_request), intent(in) :: request
         type(lr_goldstone_sumrule_result), intent(out) :: result
      end subroutine evaluate_lr_goldstone_sumrule_request

      module subroutine evaluate_lr_goldstone_sumrule_explicit(space, static_susceptibility, magnetization, result, &
         magnetization_label)
         type(response_space_layout), intent(in) :: space
         complex(rp), intent(in) :: static_susceptibility(:, :)
         real(rp), intent(in) :: magnetization(:, :)
         type(lr_goldstone_sumrule_result), intent(out) :: result
         character(len=*), intent(in), optional :: magnetization_label
      end subroutine evaluate_lr_goldstone_sumrule_explicit

      module subroutine solve_lr_goldstone_equation(gamma, rhs, result)
         complex(rp), intent(in) :: gamma(:, :), rhs(:)
         type(lr_sumrule_linear_solve_result), intent(out) :: result
      end subroutine solve_lr_goldstone_equation

      module subroutine evaluate_projected_mills_interaction(request, result)
         type(projected_mills_interaction_request), intent(in) :: request
         type(projected_mills_interaction_result), intent(out) :: result
      end subroutine evaluate_projected_mills_interaction

      module subroutine evaluate_projected_mills_from_reciprocal(contract, radial_bases, reciprocal_obj, energy, moment, result)
         type(projected_site_spin_contract), intent(in) :: contract
         type(lmto_radial_basis), intent(in) :: radial_bases(:)
         type(reciprocal), intent(in) :: reciprocal_obj
         real(rp), intent(in) :: energy, moment(:)
         type(projected_mills_interaction_result), intent(out) :: result
      end subroutine evaluate_projected_mills_from_reciprocal

      module subroutine evaluate_projected_dyson(request, result)
         type(projected_dyson_request), intent(in) :: request
         type(projected_dyson_result), intent(out) :: result
      end subroutine evaluate_projected_dyson

      module subroutine select_projected_juelich_eta_indices(eta_values, selected_index, holdout_index)
         real(rp), intent(in) :: eta_values(:)
         integer, intent(out) :: selected_index, holdout_index
      end subroutine select_projected_juelich_eta_indices

      module subroutine evaluate_projected_juelich_interaction(request, result)
         type(projected_juelich_request), intent(in) :: request
         type(projected_juelich_result), intent(out) :: result
      end subroutine evaluate_projected_juelich_interaction

      module subroutine evaluate_projected_juelich_holdout(static_chi0, projected_moment, interaction_U, residual, &
         relative_residual)
         complex(rp), intent(in) :: static_chi0(:, :)
         real(rp), intent(in) :: projected_moment(:), interaction_U(:)
         real(rp), intent(out) :: residual, relative_residual
      end subroutine evaluate_projected_juelich_holdout

      module subroutine assess_projected_juelich_eta_stability(result, previous_result, is_stable)
         type(projected_juelich_result), intent(inout) :: result
         type(projected_juelich_result), intent(in) :: previous_result
         logical, intent(out) :: is_stable
      end subroutine assess_projected_juelich_eta_stability

      module subroutine compact_project_point_vector(space, product, point_vector, compact_vector)
         type(response_space_layout), intent(in) :: space
         type(lmto_product_response_basis), intent(in) :: product
         complex(rp), intent(in) :: point_vector(:)
         complex(rp), intent(out) :: compact_vector(:)
      end subroutine compact_project_point_vector

      module subroutine compact_reconstruct_point_vector(space, product, compact_vector, point_vector)
         type(response_space_layout), intent(in) :: space
         type(lmto_product_response_basis), intent(in) :: product
         complex(rp), intent(in) :: compact_vector(:)
         complex(rp), intent(out) :: point_vector(:)
      end subroutine compact_reconstruct_point_vector

      module subroutine compact_project_local_operator_complex(space, product, values, compact_operator)
         type(response_space_layout), intent(in) :: space
         type(lmto_product_response_basis), intent(in) :: product
         complex(rp), intent(in) :: values(:, :)
         complex(rp), intent(out) :: compact_operator(:, :)
      end subroutine compact_project_local_operator_complex

      module subroutine compact_apply_local_operator(space, product, values, vector, result)
         type(response_space_layout), intent(in) :: space
         type(lmto_product_response_basis), intent(in) :: product
         complex(rp), intent(in) :: values(:, :), vector(:)
         complex(rp), intent(out) :: result(:)
      end subroutine compact_apply_local_operator

      module subroutine compact_project_magnetization(space, product, magnetization, compact_vector, point_vector)
         type(response_space_layout), intent(in) :: space
         type(lmto_product_response_basis), intent(in) :: product
         real(rp), intent(in) :: magnetization(:, :)
         complex(rp), intent(out) :: compact_vector(:)
         complex(rp), allocatable, intent(out), optional :: point_vector(:)
      end subroutine compact_project_magnetization

      module subroutine compact_weighted_projection_diagnostics(space, point_vector, projected_point_vector, weighted_norm, &
         residual_norm, relative_residual)
         type(response_space_layout), intent(in) :: space
         complex(rp), intent(in) :: point_vector(:), projected_point_vector(:)
         real(rp), intent(out) :: weighted_norm, residual_norm, relative_residual
      end subroutine compact_weighted_projection_diagnostics

      module subroutine evaluate_compact_goldstone_sumrule(space, product, static_susceptibility, magnetization, result)
         type(response_space_layout), intent(in) :: space
         type(lmto_product_response_basis), intent(in) :: product
         complex(rp), intent(in) :: static_susceptibility(:, :)
         real(rp), intent(in) :: magnetization(:, :)
         type(lr_compact_gsr_result), intent(out) :: result
      end subroutine evaluate_compact_goldstone_sumrule

      module subroutine evaluate_tddft_dyson(request, result)
         type(tddft_dyson_request), intent(in) :: request
         type(tddft_dyson_result), intent(out) :: result
      end subroutine evaluate_tddft_dyson

      module subroutine solve_tddft_dyson_frequency(chi_ks, kernel, chi, denominator, info, condition_number, &
         min_singular_value, max_singular_value)
         complex(rp), intent(in) :: chi_ks(:, :), kernel(:, :)
         complex(rp), intent(out) :: chi(:, :), denominator(:, :)
         integer, intent(out) :: info
         real(rp), intent(out), optional :: condition_number, min_singular_value, max_singular_value
      end subroutine solve_tddft_dyson_frequency

      module function tddft_loss_matrix_ordinary(chi) result(loss)
         complex(rp), intent(in) :: chi(:, :)
         complex(rp) :: loss(size(chi, 1), size(chi, 2))
      end function tddft_loss_matrix_ordinary

      module function tddft_loss_matrix_metric(space, chi) result(loss)
         type(response_space_layout), intent(in) :: space
         complex(rp), intent(in) :: chi(:, :)
         complex(rp) :: loss(size(chi, 1), size(chi, 2))
      end function tddft_loss_matrix_metric

      module pure real(rp) function loss_matrix_hermiticity_residual_ordinary(loss) result(residual)
         complex(rp), intent(in) :: loss(:, :)
      end function loss_matrix_hermiticity_residual_ordinary

      module real(rp) function loss_matrix_hermiticity_residual_metric(space, loss) result(residual)
         type(response_space_layout), intent(in) :: space
         complex(rp), intent(in) :: loss(:, :)
      end function loss_matrix_hermiticity_residual_metric

      module subroutine tddft_denominator_minimum_magnitude_eigenvalue(matrix, minimum_magnitude)
         complex(rp), intent(in) :: matrix(:, :)
         real(rp), intent(out) :: minimum_magnitude
      end subroutine tddft_denominator_minimum_magnitude_eigenvalue
      module subroutine lmto_fixture_clear(this)
         class(lmto_live_hamiltonian_fixture), intent(inout) :: this
      end subroutine lmto_fixture_clear

      module subroutine rotation_state_clear(this)
         class(rotation_state), intent(inout) :: this
      end subroutine rotation_state_clear

      module subroutine rotation_result_clear(this)
         class(rotation_result), intent(inout) :: this
      end subroutine rotation_result_clear

      ! --- generic specifics from lr_response_space_mod ---
      module subroutine response_local_operator_site_radial(space, values, operator)
      type(response_space_layout), intent(in) :: space
      real(rp), intent(in) :: values(:, :)
      complex(rp), intent(out) :: operator(:, :)
      end subroutine response_local_operator_site_radial

      module subroutine response_local_operator_site_radial_channel(space, values, operator)
      type(response_space_layout), intent(in) :: space
      real(rp), intent(in) :: values(:, :, :)
      complex(rp), intent(out) :: operator(:, :)
      end subroutine response_local_operator_site_radial_channel

      module subroutine response_local_operator_complex_site_radial(space, values, operator)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: values(:, :)
      complex(rp), intent(out) :: operator(:, :)
      end subroutine response_local_operator_complex_site_radial

      module subroutine response_local_operator_complex_site_radial_channel(space, values, operator)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: values(:, :, :)
      complex(rp), intent(out) :: operator(:, :)
      end subroutine response_local_operator_complex_site_radial_channel

      ! --- generic specifics from lr_pauli_transition_vertex_mod ---
      module subroutine evaluate_pauli_transition_vertex_one(space, radial_bases, left_state, right_state, operator_matrix, &
         capabilities, transition_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(in) :: operator_matrix(:, :)
      type(pauli_vertex_capabilities), intent(in) :: capabilities
      complex(rp), intent(out) :: transition_vector(:)
      end subroutine evaluate_pauli_transition_vertex_one

      module subroutine evaluate_pauli_transition_vertex_channels(space, radial_bases, left_state, right_state, operator_matrices, &
         capabilities, transition_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(in) :: operator_matrices(:, :, :)
      type(pauli_vertex_capabilities), intent(in) :: capabilities
      complex(rp), intent(out) :: transition_vector(:)
      end subroutine evaluate_pauli_transition_vertex_channels



      ! --- from lr_ks_susceptibility_mod ---
      module pure real(rp) function lr_fermi_dirac_occupation(eigenvalue, fermi_level, temperature) result(value)
         real(rp), intent(in) :: eigenvalue, fermi_level, temperature
      end function lr_fermi_dirac_occupation

      module subroutine lr_electronic_state_initialize(this, eigenvalues, eigenvectors, k_points, k_weights, occupations, &
         fermi_level, temperature, energy_zero, reciprocal_mode, hamiltonian_order, orthogonal, collinear, has_soc, &
         has_extra_operator)
         class(lr_electronic_state), intent(out) :: this
         real(rp), intent(in) :: eigenvalues(:, :), k_points(:, :), k_weights(:), occupations(:, :)
         complex(rp), intent(in) :: eigenvectors(:, :, :)
         real(rp), intent(in) :: fermi_level, temperature
         real(rp), intent(in), optional :: energy_zero
         character(len=*), intent(in), optional :: reciprocal_mode, hamiltonian_order
         logical, intent(in), optional :: orthogonal, collinear, has_soc, has_extra_operator
      end subroutine lr_electronic_state_initialize

      module subroutine lr_electronic_state_restore(this)
         class(lr_electronic_state), intent(inout) :: this
      end subroutine lr_electronic_state_restore

      module subroutine lr_electronic_state_validate(this, caller)
         class(lr_electronic_state), intent(in) :: this
         character(len=*), intent(in) :: caller
      end subroutine lr_electronic_state_validate

      module subroutine lr_snapshot_from_reciprocal(recip, state)
         type(reciprocal), intent(in) :: recip
         type(lr_electronic_state), intent(out) :: state
      end subroutine lr_snapshot_from_reciprocal

      module subroutine lr_q_endpoint_from_reciprocal(recip, left_state, q, endpoint_state)
         type(reciprocal), intent(inout) :: recip
         type(lr_electronic_state), intent(in) :: left_state
         real(rp), intent(in) :: q(3)
         type(lr_electronic_state), intent(out) :: endpoint_state
      end subroutine lr_q_endpoint_from_reciprocal

      module subroutine evaluate_lr_ks_susceptibility_request(request, result)
         type(lr_ks_susceptibility_request), intent(in) :: request
         type(lr_ks_susceptibility_result), intent(out) :: result
      end subroutine evaluate_lr_ks_susceptibility_request

      module subroutine evaluate_lr_ks_susceptibility_explicit(space, radial_bases, left_state, right_state, request, result)
         type(response_space_layout), intent(in) :: space
         type(lmto_radial_basis), intent(in) :: radial_bases(:)
         type(lr_electronic_state), intent(in) :: left_state, right_state
         type(lr_ks_susceptibility_request), intent(in) :: request
         type(lr_ks_susceptibility_result), intent(out) :: result
      end subroutine evaluate_lr_ks_susceptibility_explicit

      module subroutine evaluate_lr_product_ks_susceptibility_request(request, result)
         type(lr_product_ks_susceptibility_request), intent(in) :: request
         type(lr_product_ks_susceptibility_result), intent(out) :: result
      end subroutine evaluate_lr_product_ks_susceptibility_request

      module subroutine evaluate_lr_product_ks_susceptibility_explicit(product_basis, left_state, right_state, request, result)
         type(lmto_product_response_basis), intent(in) :: product_basis
         type(lr_electronic_state), intent(in) :: left_state, right_state
         type(lr_product_ks_susceptibility_request), intent(in) :: request
         type(lr_product_ks_susceptibility_result), intent(out) :: result
      end subroutine evaluate_lr_product_ks_susceptibility_explicit

      module subroutine evaluate_lr_static_residual(space, static_susceptibility, field, magnetization, absolute_residual, &
         relative_residual)
         type(response_space_layout), intent(in) :: space
         complex(rp), intent(in) :: static_susceptibility(:, :), field(:), magnetization(:)
         real(rp), intent(out) :: absolute_residual, relative_residual
      end subroutine evaluate_lr_static_residual

      ! --- from lr_rs_gf_susceptibility_mod ---
      module subroutine lr_rs_block_recursion_initialize(this, callback, nbasis_per_site, recursion_depth, terminator)
         class(lr_rs_block_recursion_provider), intent(out) :: this
         procedure(lr_rs_native_callback) :: callback
         integer, intent(in) :: nbasis_per_site, recursion_depth
         character(len=*), intent(in), optional :: terminator
      end subroutine lr_rs_block_recursion_initialize

      module subroutine lr_rs_block_recursion_get_pair(this, left_site, right_site, translation, z, gij, gji, hgamma_ij, &
         hgamma_ji)
         class(lr_rs_block_recursion_provider), intent(inout) :: this
         integer, intent(in) :: left_site, right_site
         real(rp), intent(in) :: translation(3)
         complex(rp), intent(in) :: z
         complex(rp), intent(out) :: gij(:, :), gji(:, :), hgamma_ij(:, :), hgamma_ji(:, :)
      end subroutine lr_rs_block_recursion_get_pair

      module function lr_rs_block_recursion_describe(this) result(description)
         class(lr_rs_block_recursion_provider), intent(in) :: this
         character(len=256) :: description
      end function lr_rs_block_recursion_describe

      module subroutine lr_rs_chebyshev_initialize(this, callback, nbasis_per_site, polynomial_order, kernel)
         class(lr_rs_chebyshev_provider), intent(out) :: this
         procedure(lr_rs_native_callback) :: callback
         integer, intent(in) :: nbasis_per_site, polynomial_order
         character(len=*), intent(in), optional :: kernel
      end subroutine lr_rs_chebyshev_initialize

      module subroutine lr_rs_chebyshev_get_pair(this, left_site, right_site, translation, z, gij, gji, hgamma_ij, &
         hgamma_ji)
         class(lr_rs_chebyshev_provider), intent(inout) :: this
         integer, intent(in) :: left_site, right_site
         real(rp), intent(in) :: translation(3)
         complex(rp), intent(in) :: z
         complex(rp), intent(out) :: gij(:, :), gji(:, :), hgamma_ij(:, :), hgamma_ji(:, :)
      end subroutine lr_rs_chebyshev_get_pair

      module function lr_rs_chebyshev_describe(this) result(description)
         class(lr_rs_chebyshev_provider), intent(in) :: this
         character(len=256) :: description
      end function lr_rs_chebyshev_describe

      module subroutine evaluate_lr_rs_gf_susceptibility(request, result)
         type(lr_rs_gf_susceptibility_request), intent(in) :: request
         type(lr_ks_susceptibility_result), intent(out) :: result
      end subroutine evaluate_lr_rs_gf_susceptibility

      ! --- from tddft_native_rsgf_provider_mod ---
      module subroutine tddft_native_rsgf_initialize(this, provider_name, green_obj, recursion_obj, hamiltonian_obj, &
         lattice_obj, reciprocal_obj, radial_bases)
         class(tddft_native_rsgf_provider), intent(out) :: this
         character(len=*), intent(in) :: provider_name
         type(green), target, intent(inout) :: green_obj
         type(recursion), target, intent(inout) :: recursion_obj
         type(hamiltonian), target, intent(in) :: hamiltonian_obj
         type(lattice), target, intent(in) :: lattice_obj
         type(reciprocal), target, intent(in) :: reciprocal_obj
         type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      end subroutine tddft_native_rsgf_initialize

      module subroutine tddft_native_rsgf_get_pair(this, left_site, right_site, translation, z, gij, gji, hgamma_ij, &
         hgamma_ji)
         class(tddft_native_rsgf_provider), intent(inout) :: this
         integer, intent(in) :: left_site, right_site
         real(rp), intent(in) :: translation(3)
         complex(rp), intent(in) :: z
         complex(rp), intent(out) :: gij(:, :), gji(:, :), hgamma_ij(:, :), hgamma_ji(:, :)
      end subroutine tddft_native_rsgf_get_pair

      module function tddft_native_rsgf_describe(this) result(description)
         class(tddft_native_rsgf_provider), intent(in) :: this
         character(len=256) :: description
      end function tddft_native_rsgf_describe

      ! --- from lr_projected_reciprocal_chi0_mod ---
      module subroutine evaluate_projected_lehmann_chi0(request, result)
         type(projected_chi0_request), intent(in) :: request
         type(projected_chi0_result), intent(out) :: result
      end subroutine evaluate_projected_lehmann_chi0

      module subroutine evaluate_projected_finite_width_chi0(request, result)
         type(projected_chi0_request), intent(in) :: request
         type(projected_chi0_result), intent(out) :: result
      end subroutine evaluate_projected_finite_width_chi0

      module subroutine evaluate_projected_gf_chi0(request, result)
         type(projected_chi0_request), intent(in) :: request
         type(projected_chi0_result), intent(out) :: result
      end subroutine evaluate_projected_gf_chi0

      module subroutine evaluate_projected_gf_chi0_optimized(request, result)
         type(projected_chi0_request), intent(in) :: request
         type(projected_chi0_result), intent(out) :: result
      end subroutine evaluate_projected_gf_chi0_optimized

      module subroutine evaluate_projected_gf_chi0_reference(request, result)
         type(projected_chi0_request), intent(in) :: request
         type(projected_chi0_result), intent(out) :: result
      end subroutine evaluate_projected_gf_chi0_reference

      ! --- from response_angular_basis_mod ---
      module pure integer function response_lmax(orbital_lmax) result(value)
      integer, intent(in) :: orbital_lmax
      end function response_lmax

      module pure logical function response_product_space_complete(orbital_lmax, response_lmax_in) result(ok)
      integer, intent(in) :: orbital_lmax, response_lmax_in
      end function response_product_space_complete

      module pure integer function response_lm_index(l, m) result(index)
      integer, intent(in) :: l, m
      end function response_lm_index

      module pure subroutine response_lm_from_index(index, l, m)
      integer, intent(in) :: index
      integer, intent(out) :: l, m
      end subroutine response_lm_from_index

      module pure function response_harmonic(l, m, theta, phi) result(value)
      integer, intent(in) :: l, m
      real(rp), intent(in) :: theta, phi
      complex(rp) :: value
      end function response_harmonic

      module pure real(rp) function response_gaunt(l1, m1, l2, m2, lout, mout) result(value)
      integer, intent(in) :: l1, m1, l2, m2, lout, mout
      end function response_gaunt

      module subroutine response_product_coefficients(l, m, lp, mp, coefficients)
      integer, intent(in) :: l, m, lp, mp
      real(rp), intent(out) :: coefficients(:)
      end subroutine response_product_coefficients

      ! --- from response_basis_mapping_mod ---
      module pure integer function response_superindex_size(nsite, response_lmax, npoint, nchannel) result(size_out)
      integer, intent(in) :: nsite, response_lmax, npoint, nchannel
      end function response_superindex_size

      module subroutine response_flatten_superindex(item, nsite, response_lmax, npoint, nchannel, flat)
      type(response_super_index), intent(in) :: item
      integer, intent(in) :: nsite, response_lmax, npoint, nchannel
      integer, intent(out) :: flat
      end subroutine response_flatten_superindex

      module subroutine response_unflatten_superindex(flat, nsite, response_lmax, npoint, nchannel, item)
      integer, intent(in) :: flat, nsite, response_lmax, npoint, nchannel
      type(response_super_index), intent(out) :: item
      end subroutine response_unflatten_superindex

      module pure real(rp) function response_simpson_weight(ir, npoint) result(weight)
      integer, intent(in) :: ir, npoint
      end function response_simpson_weight

      module pure real(rp) function response_log_mesh_jacobian(a, b, radius) result(value)
      real(rp), intent(in) :: a, b, radius
      end function response_log_mesh_jacobian

      module pure real(rp) function response_volume_measure(a, b, radius) result(value)
      real(rp), intent(in) :: a, b, radius
      end function response_volume_measure

      module pure real(rp) function response_weighted_density(radius, physical_density) result(value)
      real(rp), intent(in) :: radius, physical_density
      end function response_weighted_density

      module pure real(rp) function response_physical_density(radius, weighted_density) result(value)
      real(rp), intent(in) :: radius, weighted_density
      end function response_physical_density

      module function response_log_mesh_integral(values, radius, a, b) result(integral)
      real(rp), intent(in) :: values(:), radius(:), a, b
      real(rp) :: integral
      end function response_log_mesh_integral

      module pure complex(rp) function response_endpoint_phase(site_tau, reciprocal_vector) result(phase)
      real(rp), intent(in) :: site_tau(3), reciprocal_vector(3)
      end function response_endpoint_phase

      module pure complex(rp) function response_real_space_phase(reciprocal_vector, translation, tau_left, tau_right) result(phase)
      real(rp), intent(in) :: reciprocal_vector(3), translation(3), tau_left(3), tau_right(3)
      end function response_real_space_phase

      module subroutine response_site_gauge(site_tau, reciprocal_vector, gauge)
      real(rp), intent(in) :: site_tau(:, :), reciprocal_vector(3)
      complex(rp), intent(out) :: gauge(:)
      end subroutine response_site_gauge

      module subroutine response_apply_site_gauge(coefficients, site_tau, reciprocal_vector, transformed)
      complex(rp), intent(in) :: coefficients(:, :)
      real(rp), intent(in) :: site_tau(:, :), reciprocal_vector(3)
      complex(rp), intent(out) :: transformed(:, :)
      end subroutine response_apply_site_gauge

      module pure logical function response_mapping_supported(orbital_lmax, scalar_relativistic, collinear, no_soc, &
      no_extra_operator, full_scalar_relativistic_spin_angular) result(ok)
      integer, intent(in) :: orbital_lmax
      logical, intent(in) :: scalar_relativistic, collinear, no_soc, no_extra_operator
      logical, intent(in) :: full_scalar_relativistic_spin_angular
      end function response_mapping_supported

      ! --- from lr_response_space_mod ---
      module subroutine response_build_radial_metric(a, b, radius, weights)
      real(rp), intent(in) :: a, b, radius(:)
      real(rp), intent(out) :: weights(:)
      end subroutine response_build_radial_metric

      module subroutine response_metric_weights(space, weights)
      type(response_space_layout), intent(in) :: space
      real(rp), intent(out) :: weights(:)
      end subroutine response_metric_weights

      module function response_vector_inner_product(space, left, right) result(value)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: left(:), right(:)
      complex(rp) :: value
      end function response_vector_inner_product

      module function response_vector_norm(space, vector) result(value)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: vector(:)
      real(rp) :: value
      end function response_vector_norm

      module subroutine response_apply_operator(space, operator, field, result)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: operator(:, :), field(:)
      complex(rp), intent(out) :: result(:)
      end subroutine response_apply_operator

      module subroutine response_compose_operators(space, first, second, composed)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: first(:, :), second(:, :)
      complex(rp), intent(out) :: composed(:, :)
      end subroutine response_compose_operators

      module subroutine response_identity_operator(space, identity)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(out) :: identity(:, :)
      end subroutine response_identity_operator

      module subroutine response_operator_adjoint(space, operator, adjoint)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: operator(:, :)
      complex(rp), intent(out) :: adjoint(:, :)
      end subroutine response_operator_adjoint

      module function response_operator_trace(space, operator) result(value)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: operator(:, :)
      complex(rp) :: value
      end function response_operator_trace

      module function response_rigid_vector_norm(space, rigid_vector) result(value)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: rigid_vector(:)
      real(rp) :: value
      end function response_rigid_vector_norm

      module function response_rigid_vector_overlap(space, vector, rigid_vector) result(value)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: vector(:), rigid_vector(:)
      complex(rp) :: value
      end function response_rigid_vector_overlap

      module subroutine response_raw_to_canonical(space, raw, canonical)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: raw(:, :)
      complex(rp), intent(out) :: canonical(:, :)
      end subroutine response_raw_to_canonical

      module subroutine response_canonical_to_raw(space, canonical, raw)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: canonical(:, :)
      complex(rp), intent(out) :: raw(:, :)
      end subroutine response_canonical_to_raw

      module subroutine response_space_initialize(this, nsite, response_lmax, radius, a, b, nchannel)
      class(response_space_layout), intent(out) :: this
      integer, intent(in) :: nsite, response_lmax, nchannel
      real(rp), intent(in) :: radius(:), a, b
      end subroutine response_space_initialize

      ! --- from lr_lmto_endpoint_branches_mod ---
      module pure subroutine lmto_product_branch_powers(branch, p, q)
      integer, intent(in) :: branch
      integer, intent(out) :: p, q
      end subroutine lmto_product_branch_powers

      module pure character(len=2) function lmto_product_branch_label(branch) result(label)
      integer, intent(in) :: branch
      end function lmto_product_branch_label

      module pure logical function lmto_product_branch_valid(branch) result(valid)
      integer, intent(in) :: branch
      end function lmto_product_branch_valid

      module pure real(rp) function lmto_product_energy_power(energy, power) result(value)
      real(rp), intent(in) :: energy
      integer, intent(in) :: power
      end function lmto_product_energy_power

      module recursive subroutine lmto_product_second_order_radial_branch(radial, ir, l, lp, spin_left, spin_right, branch, value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, lp, spin_left, spin_right, branch
      real(rp), intent(out) :: value
      end subroutine lmto_product_second_order_radial_branch

      module subroutine lmto_product_apply_branch_action(hamiltonian, matrix, branch, action, stored_dual)
      complex(rp), intent(in) :: hamiltonian(:, :), matrix(:, :)
      integer, intent(in) :: branch
      complex(rp), intent(out) :: action(:, :)
      logical, intent(in), optional :: stored_dual
      end subroutine lmto_product_apply_branch_action

      ! --- from lr_pauli_transition_vertex_mod ---
      module pure function pauli_charge_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)
      end function pauli_charge_matrix

      module pure function pauli_sigma_x_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)
      end function pauli_sigma_x_matrix

      module pure function pauli_sigma_y_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)
      end function pauli_sigma_y_matrix

      module pure function pauli_sigma_z_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)
      end function pauli_sigma_z_matrix

      module pure function pauli_sigma_plus_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)
      end function pauli_sigma_plus_matrix

      module pure function pauli_sigma_minus_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)
      end function pauli_sigma_minus_matrix

      module subroutine pauli_endpoint_state_initialize(this, energy, coefficients)
      class(pauli_endpoint_state), intent(out) :: this
      real(rp), intent(in) :: energy
      complex(rp), intent(in) :: coefficients(:)
      end subroutine pauli_endpoint_state_initialize

      ! --- from lr_lmto_product_response_basis_mod ---
      module pure integer function lmto_product_candidate_count(lmax, response_l) result(count)
      integer, intent(in) :: lmax, response_l
      end function lmto_product_candidate_count

      module subroutine lmto_enumerate_product_candidates(lmax, response_l, candidates)
      integer, intent(in) :: lmax, response_l
      type(lmto_product_candidate), allocatable, intent(out) :: candidates(:)
      end subroutine lmto_enumerate_product_candidates

      module subroutine lmto_product_response_basis_initialize(this, space, radial_bases, circular_channel, strict_rank)
      class(lmto_product_response_basis), intent(out) :: this
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      integer, intent(in) :: circular_channel
      logical, intent(in), optional :: strict_rank
      end subroutine lmto_product_response_basis_initialize

      module integer function lmto_product_flat_index(this, site, response_l, response_m, product_mode) result(flat)
      class(lmto_product_response_basis), intent(in) :: this
      integer, intent(in) :: site, response_l, response_m, product_mode
      end function lmto_product_flat_index

      module subroutine lmto_product_unflatten_index(this, flat, site, response_l, response_m, product_mode)
      class(lmto_product_response_basis), intent(in) :: this
      integer, intent(in) :: flat
      integer, intent(out) :: site, response_l, response_m, product_mode
      end subroutine lmto_product_unflatten_index

      module subroutine lmto_product_candidate_coefficients(this, site, response_l, response_m, left_state, right_state, coefficients, &
      selected_l)
      class(lmto_product_response_basis), intent(in) :: this
      integer, intent(in) :: site, response_l, response_m
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(out) :: coefficients(:)
      logical, intent(in), optional :: selected_l(0:)
      end subroutine lmto_product_candidate_coefficients

      module subroutine lmto_product_transition_coordinates(this, left_state, right_state, coordinates, selected_l, response_l_filter)
      class(lmto_product_response_basis), intent(in) :: this
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(out) :: coordinates(:)
      logical, intent(in), optional :: selected_l(0:)
      integer, intent(in), optional :: response_l_filter
      end subroutine lmto_product_transition_coordinates

      module subroutine lmto_product_component_vertex_tensor(this, vertices, selected_l)
      class(lmto_product_response_basis), intent(in) :: this
      complex(rp), allocatable, intent(out) :: vertices(:, :, :, :)
      logical, intent(in), optional :: selected_l(0:)
      end subroutine lmto_product_component_vertex_tensor

      ! --- from lr_gf_endpoint_augmentation_mod ---
      module subroutine augment_lr_gf_endpoint(radial_left, radial_right, left_site, right_site, z, coefficient_gf, &
      augmented, capabilities, hgamma_block, hamiltonian_provider, &
      reverse_coefficient_gf, reverse_hgamma_block)
      type(lmto_radial_basis), intent(in) :: radial_left, radial_right
      integer, intent(in) :: left_site, right_site
      complex(rp), intent(in) :: z
      complex(rp), intent(in) :: coefficient_gf(:, :)
      type(lr_gf_augmented_block), intent(out) :: augmented
      type(lr_gf_endpoint_capabilities), intent(in), optional :: capabilities
      complex(rp), intent(in), optional :: hgamma_block(:, :)
      class(lr_gf_effective_hamiltonian_provider), intent(inout), optional :: hamiltonian_provider
      complex(rp), intent(in), optional :: reverse_coefficient_gf(:, :)
      complex(rp), intent(in), optional :: reverse_hgamma_block(:, :)
      end subroutine augment_lr_gf_endpoint

      module subroutine augment_lr_gf_endpoint_pair(radial_left, radial_right, left_site, right_site, z, coefficient_gf, &
      reverse_coefficient_gf, forward, reverse, capabilities, hgamma_block, &
      reverse_hgamma_block, hamiltonian_provider)
      type(lmto_radial_basis), intent(in) :: radial_left, radial_right
      integer, intent(in) :: left_site, right_site
      complex(rp), intent(in) :: z
      complex(rp), intent(in) :: coefficient_gf(:, :), reverse_coefficient_gf(:, :)
      type(lr_gf_augmented_block), intent(out) :: forward, reverse
      type(lr_gf_endpoint_capabilities), intent(in), optional :: capabilities
      complex(rp), intent(in), optional :: hgamma_block(:, :), reverse_hgamma_block(:, :)
      class(lr_gf_effective_hamiltonian_provider), intent(inout), optional :: hamiltonian_provider
      end subroutine augment_lr_gf_endpoint_pair

      module pure logical function lr_gf_endpoint_capabilities_supported(this) result(ok)
      class(lr_gf_endpoint_capabilities), intent(in) :: this
      end function lr_gf_endpoint_capabilities_supported

      module subroutine lr_gf_endpoint_require_supported(this, caller)
      class(lr_gf_endpoint_capabilities), intent(in) :: this
      character(len=*), intent(in), optional :: caller
      end subroutine lr_gf_endpoint_require_supported

      module subroutine lr_gf_dense_provider_initialize(this, h_eff, block_size)
      class(lr_gf_dense_hamiltonian_provider), intent(out) :: this
      complex(rp), intent(in) :: h_eff(:, :)
      integer, intent(in) :: block_size
      end subroutine lr_gf_dense_provider_initialize

      module subroutine lr_gf_dense_provider_apply(this, source_site, destination_site, seed, destination_block)
      class(lr_gf_dense_hamiltonian_provider), intent(inout) :: this
      integer, intent(in) :: source_site, destination_site
      complex(rp), intent(in) :: seed(:, :)
      complex(rp), intent(out) :: destination_block(:, :)
      end subroutine lr_gf_dense_provider_apply

      module subroutine lr_gf_augmented_block_evaluate_point(this, left_radial_index, left_theta, left_phi, &
      right_radial_index, right_theta, right_phi, point_gf)
      class(lr_gf_augmented_block), intent(in) :: this
      integer, intent(in) :: left_radial_index, right_radial_index
      real(rp), intent(in) :: left_theta, left_phi, right_theta, right_phi
      complex(rp), intent(out) :: point_gf(:, :)
      end subroutine lr_gf_augmented_block_evaluate_point

      module subroutine lr_gf_augmented_block_branch_at_point(this, left_radial_index, left_theta, left_phi, &
      right_radial_index, right_theta, right_phi, branches)
      class(lr_gf_augmented_block), intent(in) :: this
      integer, intent(in) :: left_radial_index, right_radial_index
      real(rp), intent(in) :: left_theta, left_phi, right_theta, right_phi
      complex(rp), intent(out) :: branches(:, :, :)
      complex(rp) :: phi_left(2, 2*this%norb), dot_left(2, 2*this%norb)
      complex(rp) :: phi_right(2, 2*this%norb), dot_right(2, 2*this%norb)
      end subroutine lr_gf_augmented_block_branch_at_point

      ! --- from lr_projected_site_spin_mod ---
      module pure logical function projected_selector_supported(selection) result(ok)
      character(len=*), intent(in) :: selection
      end function projected_selector_supported

      module subroutine projected_site_spin_initialize(this, space, radial_bases, selection)
      class(projected_site_spin_contract), intent(out) :: this
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      character(len=*), intent(in) :: selection
      end subroutine projected_site_spin_initialize

      module subroutine projected_site_spin_functional(this, product, functionals)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(out) :: functionals(:, :)
      end subroutine projected_site_spin_functional

      module subroutine projected_selected_coordinates(this, product, left_state, right_state, coordinates)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_product_response_basis), intent(in) :: product
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(out) :: coordinates(:)
      end subroutine projected_selected_coordinates

      module subroutine projected_transition_amplitudes(this, product, left_state, right_state, amplitudes)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_product_response_basis), intent(in) :: product
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(out) :: amplitudes(:)
      end subroutine projected_transition_amplitudes

      module subroutine projected_direct_operator_matrix(this, radial_bases, energy, operator_kind, matrix)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      real(rp), intent(in) :: energy
      integer, intent(in) :: operator_kind
      complex(rp), intent(out) :: matrix(:, :)
      end subroutine projected_direct_operator_matrix

      module subroutine projected_direct_transition_amplitudes(this, radial_bases, left_state, right_state, operator_kind, amplitudes)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      integer, intent(in) :: operator_kind
      complex(rp), intent(out) :: amplitudes(:)
      end subroutine projected_direct_transition_amplitudes

      module subroutine projected_site_component_vertex_tensor(this, product, vertices)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), allocatable, intent(out) :: vertices(:, :, :, :)
      end subroutine projected_site_component_vertex_tensor

      module subroutine projected_moment_from_density(this, eigenvalues, eigenvectors, k_weights, fermi_level, temperature, &
      ground_states, moment)
      class(projected_site_spin_contract), intent(in) :: this
      real(rp), intent(in) :: eigenvalues(:, :), k_weights(:), fermi_level, temperature
      complex(rp), intent(in) :: eigenvectors(:, :, :)
      type(radial_ground_state), intent(in) :: ground_states(:)
      real(rp), intent(out) :: moment(:)
      end subroutine projected_moment_from_density

      module subroutine projected_moment_from_operator(this, eigenvalues, eigenvectors, k_weights, fermi_level, temperature, &
      ground_states, moment)
      class(projected_site_spin_contract), intent(in) :: this
      real(rp), intent(in) :: eigenvalues(:, :), k_weights(:), fermi_level, temperature
      complex(rp), intent(in) :: eigenvectors(:, :, :)
      type(radial_ground_state), intent(in) :: ground_states(:)
      real(rp), intent(out) :: moment(:)
      end subroutine projected_moment_from_operator

      module subroutine projected_core_spin_number(this, ground_states, core_spin)
      class(projected_site_spin_contract), intent(in) :: this
      type(radial_ground_state), intent(in) :: ground_states(:)
      real(rp), intent(out) :: core_spin(:)
      end subroutine projected_core_spin_number

      ! --- from lr_kl_hessian ---
      module subroutine lmto_fixture_init(this, nsite, norb, nbond, hoh)
         class(lmto_live_hamiltonian_fixture), intent(out) :: this
         integer, intent(in) :: nsite, norb, nbond
         logical, intent(in), optional :: hoh
      end subroutine lmto_fixture_init

      module subroutine lmto_fixture_from_hamiltonian(source, fixture)
         type(hamiltonian), intent(in) :: source
         type(lmto_live_hamiltonian_fixture), intent(out) :: fixture
      end subroutine lmto_fixture_from_hamiltonian

      module subroutine assemble_lmto_hamiltonian(this, k_point, hamiltonian)
         type(lmto_live_hamiltonian_fixture), intent(in) :: this
         real(rp), intent(in) :: k_point(3)
         complex(rp), intent(out) :: hamiltonian(:, :)
      end subroutine assemble_lmto_hamiltonian

      module subroutine assemble_lmto_finite_q_torque(this, k_point, q_point, site, axis, torque)
         type(lmto_live_hamiltonian_fixture), intent(in) :: this
         real(rp), intent(in) :: k_point(3), q_point(3), axis(3)
         integer, intent(in) :: site
         complex(rp), intent(out) :: torque(:, :)
      end subroutine assemble_lmto_finite_q_torque

      module subroutine assemble_lmto_finite_q_torques(this, k_point, q_point, axes, torques)
         type(lmto_live_hamiltonian_fixture), intent(in) :: this
         real(rp), intent(in) :: k_point(3), q_point(3), axes(:, :)
         complex(rp), intent(out) :: torques(:, :, :)
      end subroutine assemble_lmto_finite_q_torques

      module subroutine assemble_lmto_finite_q_mixed_derivative(this, k_point, q_point, site_i, axis_i, site_j, axis_j, mixed)
         type(lmto_live_hamiltonian_fixture), intent(in) :: this
         real(rp), intent(in) :: k_point(3), q_point(3), axis_i(3), axis_j(3)
         integer, intent(in) :: site_i, site_j
         complex(rp), intent(out) :: mixed(:, :)
      end subroutine assemble_lmto_finite_q_mixed_derivative

      module subroutine force_theorem_finite_q_hessian_from_eigenbasis(eigenvalues, eigenvectors, endpoint_values, endpoint_vectors, &
         fermi, torques_q, torques_minus_q, mixed, hessian, torque_torque, mixed_contact, complete)
         real(rp), intent(in) :: eigenvalues(:), endpoint_values(:), fermi
         complex(rp), intent(in) :: eigenvectors(:, :), endpoint_vectors(:, :)
         complex(rp), intent(in) :: torques_q(:, :, :), torques_minus_q(:, :, :), mixed(:, :, :, :)
         complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      end subroutine force_theorem_finite_q_hessian_from_eigenbasis

      module subroutine force_theorem_finite_q_hessian_from_eigenbasis_batch(eigenvalues, eigenvectors, endpoint_values, endpoint_vectors, &
         fermi, weights, torques_q, torques_minus_q, mixed, hessian, torque_torque, mixed_contact, complete)
         real(rp), intent(in) :: eigenvalues(:, :), endpoint_values(:, :), fermi, weights(:)
         complex(rp), intent(in) :: eigenvectors(:, :, :), endpoint_vectors(:, :, :)
         complex(rp), intent(in) :: torques_q(:, :, :, :), torques_minus_q(:, :, :, :), mixed(:, :, :, :, :)
         complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      end subroutine force_theorem_finite_q_hessian_from_eigenbasis_batch

      module subroutine force_theorem_finite_q_hessian_from_eigenbasis_metallic(eigenvalues, eigenvectors, endpoint_values, endpoint_vectors, &
         fermi, kT, torques_q, torques_minus_q, mixed, hessian, torque_torque, mixed_contact, complete)
         real(rp), intent(in) :: eigenvalues(:), endpoint_values(:), fermi, kT
         complex(rp), intent(in) :: eigenvectors(:, :), endpoint_vectors(:, :)
         complex(rp), intent(in) :: torques_q(:, :, :), torques_minus_q(:, :, :), mixed(:, :, :, :)
         complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      end subroutine force_theorem_finite_q_hessian_from_eigenbasis_metallic

      module subroutine force_theorem_finite_q_hessian_from_eigenbasis_metallic_batch(eigenvalues, eigenvectors, endpoint_values, endpoint_vectors, &
         fermi, kT, weights, torques_q, torques_minus_q, mixed, hessian, torque_torque, mixed_contact, complete)
         real(rp), intent(in) :: eigenvalues(:, :), endpoint_values(:, :), fermi, kT, weights(:)
         complex(rp), intent(in) :: eigenvectors(:, :, :), endpoint_vectors(:, :, :)
         complex(rp), intent(in) :: torques_q(:, :, :, :), torques_minus_q(:, :, :, :), mixed(:, :, :, :, :)
         complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
      end subroutine force_theorem_finite_q_hessian_from_eigenbasis_metallic_batch

      module pure real(rp) function finite_temperature_occupation(eigenvalue, fermi, kT) result(occupation)
         real(rp), intent(in) :: eigenvalue, fermi, kT
      end function finite_temperature_occupation

      module pure real(rp) function fermi_divided_difference(e1, e2, fermi, kT) result(kernel)
         real(rp), intent(in) :: e1, e2, fermi, kT
      end function fermi_divided_difference

      module subroutine lmto_fixture_adapter_residual(source, fixture, k_points, max_error)
         type(hamiltonian), intent(in) :: source
         type(lmto_live_hamiltonian_fixture), intent(in) :: fixture
         real(rp), intent(in) :: k_points(:, :)
         real(rp), intent(out) :: max_error
      end subroutine lmto_fixture_adapter_residual

      ! --- from lr_kl_contour ---
      module subroutine build_finite_temperature_contour(h_source, h_endpoint, fermi, kT, options, nodes, weights, fermi_poles)
         complex(rp), intent(in) :: h_source(:, :), h_endpoint(:, :)
         real(rp), intent(in) :: fermi, kT
         type(finite_h_contour_options), intent(in) :: options
         complex(rp), allocatable, intent(out) :: nodes(:), weights(:), fermi_poles(:)
      end subroutine build_finite_temperature_contour

      module subroutine build_zero_temperature_occupied_contour(occupied_bounds, options, nodes, weights)
         real(rp), intent(in) :: occupied_bounds(2)
         type(finite_h_contour_options), intent(in) :: options
         complex(rp), allocatable, intent(out) :: nodes(:), weights(:)
      end subroutine build_zero_temperature_occupied_contour

      module pure complex(rp) function finite_temperature_complex_fermi(z, fermi, kT) result(value)
         complex(rp), intent(in) :: z
         real(rp), intent(in) :: fermi, kT
      end function finite_temperature_complex_fermi

      module pure complex(rp) function finite_temperature_regularized_fermi(z, fermi, kT, fermi_poles) result(value)
         complex(rp), intent(in) :: z, fermi_poles(:)
         real(rp), intent(in) :: fermi, kT
      end function finite_temperature_regularized_fermi

      module subroutine force_theorem_finite_q_hessian_from_resolvent(h_source, h_endpoint, fermi, kT, torques_q, torques_minus_q, mixed, &
         options, hessian, torque_torque, mixed_contact, complete, report)
         complex(rp), intent(in) :: h_source(:, :), h_endpoint(:, :)
         real(rp), intent(in) :: fermi, kT
         complex(rp), intent(in) :: torques_q(:, :, :), torques_minus_q(:, :, :), mixed(:, :, :, :)
         type(finite_h_contour_options), intent(in) :: options
         complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
         type(finite_h_contour_report), intent(out), optional :: report
      end subroutine force_theorem_finite_q_hessian_from_resolvent

      module subroutine force_theorem_finite_q_hessian_from_resolvent_batch(h_source, h_endpoint, fermi, kT, k_weights, torques_q, &
         torques_minus_q, mixed, options, hessian, torque_torque, mixed_contact, complete, report)
         complex(rp), intent(in) :: h_source(:, :, :), h_endpoint(:, :, :)
         real(rp), intent(in) :: fermi, kT, k_weights(:)
         complex(rp), intent(in) :: torques_q(:, :, :, :), torques_minus_q(:, :, :, :), mixed(:, :, :, :, :)
         type(finite_h_contour_options), intent(in) :: options
         complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
         type(finite_h_contour_report), intent(out), optional :: report
      end subroutine force_theorem_finite_q_hessian_from_resolvent_batch

      module subroutine force_theorem_finite_q_hessian_from_zero_temperature_contour(h_source, h_endpoint, occupied_bounds, torques_q, &
         torques_minus_q, mixed, options, hessian, torque_torque, mixed_contact, complete, report)
         complex(rp), intent(in) :: h_source(:, :), h_endpoint(:, :)
         real(rp), intent(in) :: occupied_bounds(2)
         complex(rp), intent(in) :: torques_q(:, :, :), torques_minus_q(:, :, :), mixed(:, :, :, :)
         type(finite_h_contour_options), intent(in) :: options
         complex(rp), intent(out) :: hessian(:, :), torque_torque(:, :), mixed_contact(:, :), complete(:, :)
         type(finite_h_contour_report), intent(out), optional :: report
      end subroutine force_theorem_finite_q_hessian_from_zero_temperature_contour

      ! --- from lr_rotation_response ---
      module subroutine prepare_rotation_response(fixture, recip, q, state)
         type(lmto_live_hamiltonian_fixture), target, intent(in) :: fixture
         type(reciprocal), target, intent(inout) :: recip
         real(rp), intent(in) :: q(3)
         type(rotation_state), intent(inout) :: state
      end subroutine prepare_rotation_response

      module subroutine evaluate_rotation_response(request,result)
         type(rotation_request), intent(in) :: request
         type(rotation_result), intent(inout) :: result
      end subroutine evaluate_rotation_response

      module subroutine evaluate_rotation_response_oracle(fixture,recip,q,omega,eta,bubble,contact,kernel)
         type(lmto_live_hamiltonian_fixture), intent(in) :: fixture
         type(reciprocal), intent(inout) :: recip
         real(rp), intent(in) :: q(3),omega,eta
         complex(rp), intent(out) :: bubble(:, :),contact(:, :),kernel(:, :)
      end subroutine evaluate_rotation_response_oracle

      module pure subroutine reduce_static_rotation_kernel(kernel_q,reduced)
         complex(rp), intent(in) :: kernel_q(:, :)
         complex(rp), intent(out) :: reduced(:, :)
      end subroutine reduce_static_rotation_kernel

      module pure subroutine rotation_axes(moments,axes)
         real(rp), intent(in) :: moments(:, :)
         real(rp), allocatable, intent(out) :: axes(:, :)
      end subroutine rotation_axes

      module pure subroutine rotation_circular_unitary(nsite,unitary)
         integer, intent(in) :: nsite
         complex(rp), intent(out) :: unitary(:, :)
      end subroutine rotation_circular_unitary

      module subroutine validate_rotation_capability(n_q_points, native_crosscheck, native_turek, control_obj, ham, self_obj, recip)
         integer, intent(in) :: n_q_points
         logical, intent(in) :: native_crosscheck, native_turek
         type(control), intent(in) :: control_obj
         type(hamiltonian), intent(in) :: ham
         type(self), intent(in) :: self_obj
         type(reciprocal), intent(in) :: recip
      end subroutine validate_rotation_capability

      module subroutine validate_rotation_pole_workflow_capability(fixture)
         type(lmto_live_hamiltonian_fixture), intent(in) :: fixture
      end subroutine validate_rotation_pole_workflow_capability

      module subroutine detect_commensurate_endpoint_map(k_points, k_weights, nk_mesh, q_point, endpoint_index, commensurate, residual)
         real(rp), intent(in) :: k_points(:, :), k_weights(:), q_point(3)
         integer, intent(in) :: nk_mesh(3)
         integer, intent(out) :: endpoint_index(:)
         logical, intent(out) :: commensurate
         real(rp), intent(out) :: residual
      end subroutine detect_commensurate_endpoint_map

      module subroutine compute_static_rotation_curvature(q_coordinates, q_list, rotation_axis, finite_h_spectral_mode, &
         finite_h_response_backend, contour_points, contour_shape, contour_margin, contour_height_fraction, &
         contour_account_fermi_poles, native_crosscheck, native_turek, native_contour_points, native_contour_margin, &
         native_contour_height_fraction, native_contour_account_fermi_poles, native_contour_target_fermi_poles, &
         lattice_obj, hamiltonian_obj, energy_obj, self_obj, reciprocal_obj, fixture, q_direct, q_cart, finite_total, &
         finite_tt, finite_contact, spectral_total, spectral_tt, spectral_contact, contour_total, contour_tt, contour_contact, &
         native_jq_ud, native_jq_du, native_jq_sym, native_delta_j, native_curvature, native_report, native_ready, &
         endpoint_mode, endpoint_reused, q_commensurate, endpoint_residual, endpoint_seconds, assembly_seconds, &
         contraction_seconds, hamiltonian_seconds, gf_seconds, solve_seconds, contour_seconds, total_response_seconds)
         character(len=*), intent(in) :: q_coordinates, finite_h_spectral_mode, finite_h_response_backend, contour_shape
         real(rp), intent(in) :: q_list(:, :), rotation_axis(3), contour_margin, contour_height_fraction
         logical, intent(in) :: contour_account_fermi_poles, native_crosscheck, native_turek
         integer, intent(in) :: contour_points, native_contour_points, native_contour_target_fermi_poles
         real(rp), intent(in) :: native_contour_margin, native_contour_height_fraction
         logical, intent(in) :: native_contour_account_fermi_poles
         type(lattice), intent(inout) :: lattice_obj
         type(hamiltonian), intent(inout) :: hamiltonian_obj
         type(energy_type), intent(in) :: energy_obj
         type(self), intent(inout) :: self_obj
         type(reciprocal), intent(inout) :: reciprocal_obj
         type(lmto_live_hamiltonian_fixture), intent(out) :: fixture
         real(rp), allocatable, intent(out) :: q_direct(:, :), q_cart(:, :)
         real(rp), allocatable, intent(out) :: finite_total(:, :, :), finite_tt(:, :, :), finite_contact(:, :, :)
         real(rp), allocatable, intent(out) :: spectral_total(:, :, :), spectral_tt(:, :, :), spectral_contact(:, :, :)
         real(rp), allocatable, intent(out) :: contour_total(:, :, :), contour_tt(:, :, :), contour_contact(:, :, :)
         real(rp), allocatable, intent(out) :: native_jq_ud(:), native_jq_du(:), native_jq_sym(:), native_delta_j(:), native_curvature(:)
         type(native_turek_contour_report), intent(out) :: native_report
         logical, intent(out) :: native_ready
         character(len=24), allocatable, intent(out) :: endpoint_mode(:)
         logical, allocatable, intent(out) :: endpoint_reused(:), q_commensurate(:)
         real(rp), allocatable, intent(out) :: endpoint_residual(:)
         real(rp), intent(out) :: endpoint_seconds, assembly_seconds, contraction_seconds, hamiltonian_seconds
         real(rp), intent(out) :: gf_seconds, solve_seconds, contour_seconds, total_response_seconds
      end subroutine compute_static_rotation_curvature

      module subroutine lr_config_restore_to_default(this)
         class(linear_response_config), intent(out) :: this
      end subroutine lr_config_restore_to_default

      module subroutine lr_restore_to_default(this)
         class(linear_response), intent(out) :: this
      end subroutine lr_restore_to_default

      module subroutine lr_load_config(this, filename, validate_request)
         class(linear_response), intent(inout) :: this
         character(len=*), intent(in) :: filename
         logical, intent(in), optional :: validate_request
      end subroutine lr_load_config

      module subroutine lr_validate_capability(this, control_obj, hamiltonian_obj, self_obj, reciprocal_obj)
         class(linear_response), intent(in) :: this
         type(control), intent(in) :: control_obj
         type(hamiltonian), intent(in) :: hamiltonian_obj
         type(self), intent(in) :: self_obj
         type(reciprocal), intent(in) :: reciprocal_obj
      end subroutine lr_validate_capability

      module subroutine lr_run(this, control_obj, lattice_obj, hamiltonian_obj, energy_obj, self_obj, reciprocal_obj, &
                               recursion_obj, green_obj, scf_converged)
         class(linear_response), intent(in) :: this
         type(control), intent(in) :: control_obj
         type(lattice), intent(inout) :: lattice_obj
         type(hamiltonian), intent(inout) :: hamiltonian_obj
         type(energy_type), intent(in) :: energy_obj
         type(self), intent(inout) :: self_obj
         type(reciprocal), intent(inout) :: reciprocal_obj
         type(recursion), target, intent(inout), optional :: recursion_obj
         type(green), target, intent(inout), optional :: green_obj
         logical, intent(in), optional :: scf_converged
      end subroutine lr_run

      module subroutine lr_run_rotation(this, control_obj, lattice_obj, hamiltonian_obj, energy_obj, self_obj, reciprocal_obj)
         class(linear_response), intent(in) :: this
         type(control), intent(in) :: control_obj
         type(lattice), intent(inout) :: lattice_obj
         type(hamiltonian), intent(inout) :: hamiltonian_obj
         type(energy_type), intent(in) :: energy_obj
         type(self), intent(inout) :: self_obj
         type(reciprocal), intent(inout) :: reciprocal_obj
      end subroutine lr_run_rotation
   end interface

end module linear_response_mod
