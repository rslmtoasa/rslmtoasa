!> Linear-response services and (from LR-REF-02b) the linear_response class.
!> Rotation dynamics: second-order local-rotation response of the accepted H2.
module linear_response_mod
   use precision_mod, only: rp
   use math_mod, only: pi
   use control_mod, only: control
   use energy_mod, only: energy
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use reciprocal_mod, only: reciprocal
   use self_mod, only: self
   use string_mod, only: sl
   use lr_lmto_turek_contour_mod, only: native_turek_contour_report
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
   end type rotation_result

   integer, parameter, public :: linear_response_max_points = 2000

   type, public :: linear_response_config
      character(len=sl) :: fname = ''
      character(len=16) :: formulation = 'rotation'
      character(len=16) :: q_coordinates = 'direct'
      character(len=sl) :: q_file = ''
      integer :: n_q_points = 0
      real(rp), allocatable :: q_list(:, :)
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
      character(len=sl) :: output_file = 'rotation_dynamics.dat'
   contains
      procedure :: restore_to_default => lr_config_restore_to_default
   end type linear_response_config

   type, public :: linear_response
      type(linear_response_config) :: config
   contains
      procedure :: restore_to_default => lr_restore_to_default
      procedure :: load_config => lr_load_config
      procedure :: validate_capability => lr_validate_capability
      procedure :: run => lr_run
   end type linear_response

   public :: lmto_fixture_init
   public :: lmto_fixture_from_hamiltonian
   public :: assemble_lmto_hamiltonian
   public :: assemble_lmto_torque
   public :: assemble_lmto_rotation_terms
   public :: assemble_lmto_mixed_derivative
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
   public :: force_theorem_integrand
   public :: force_theorem_hessian_from_green
   public :: force_theorem_hessian_from_eigenbasis
   public :: mixed_second_difference
   public :: grand_potential_from_eigenvalues
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


   interface
      module subroutine lmto_fixture_clear(this)
         class(lmto_live_hamiltonian_fixture), intent(inout) :: this
      end subroutine lmto_fixture_clear

      module subroutine rotation_state_clear(this)
         class(rotation_state), intent(inout) :: this
      end subroutine rotation_state_clear

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

      module subroutine assemble_lmto_torque(this, k_point, site, axis, torque)
         type(lmto_live_hamiltonian_fixture), intent(in) :: this
         real(rp), intent(in) :: k_point(3), axis(3)
         integer, intent(in) :: site
         complex(rp), intent(out) :: torque(:, :)
      end subroutine assemble_lmto_torque

      module subroutine assemble_lmto_rotation_terms(this, k_point, site, axis, b, q, enu, bi, qi, enui, torque)
         type(lmto_live_hamiltonian_fixture), intent(in) :: this
         real(rp), intent(in) :: k_point(3), axis(3)
         integer, intent(in) :: site
         complex(rp), intent(out) :: b(:, :), q(:, :), enu(:, :), bi(:, :), qi(:, :), enui(:, :), torque(:, :)
      end subroutine assemble_lmto_rotation_terms

      module subroutine assemble_lmto_mixed_derivative(this, k_point, site_i, axis_i, site_j, axis_j, mixed)
         type(lmto_live_hamiltonian_fixture), intent(in) :: this
         real(rp), intent(in) :: k_point(3), axis_i(3), axis_j(3)
         integer, intent(in) :: site_i, site_j
         complex(rp), intent(out) :: mixed(:, :)
      end subroutine assemble_lmto_mixed_derivative

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

      module pure subroutine force_theorem_integrand(torque_i, green, torque_j, mixed, torque_torque, mixed_contact, complete, include_contact)
         complex(rp), intent(in) :: torque_i(:, :), green(:, :), torque_j(:, :), mixed(:, :)
         real(rp), intent(out) :: torque_torque, mixed_contact, complete
         logical, intent(in), optional :: include_contact
      end subroutine force_theorem_integrand

      module subroutine force_theorem_hessian_from_green(greens, weights, torque_i, torque_j, mixed, &
                                                           torque_torque, mixed_contact, complete, include_contact)
         complex(rp), intent(in) :: greens(:, :, :), torque_i(:, :), torque_j(:, :), mixed(:, :)
         real(rp), intent(in) :: weights(:)
         real(rp), intent(out) :: torque_torque, mixed_contact, complete
         logical, intent(in), optional :: include_contact
      end subroutine force_theorem_hessian_from_green

      module subroutine force_theorem_hessian_from_eigenbasis(eigenvalues, eigenvectors, fermi, torque_i, torque_j, mixed, &
                                                               torque_torque, mixed_contact, complete)
         real(rp), intent(in) :: eigenvalues(:), fermi
         complex(rp), intent(in) :: eigenvectors(:, :), torque_i(:, :), torque_j(:, :), mixed(:, :)
         real(rp), intent(out) :: torque_torque, mixed_contact, complete
      end subroutine force_theorem_hessian_from_eigenbasis

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

      module pure function mixed_second_difference(omega_pp, omega_pm, omega_mp, omega_mm, delta_i, delta_j) result(hessian)
         real(rp), intent(in) :: omega_pp, omega_pm, omega_mp, omega_mm, delta_i, delta_j
         real(rp) :: hessian
      end function mixed_second_difference

      module pure function grand_potential_from_eigenvalues(eigenvalues, fermi) result(omega)
         real(rp), intent(in) :: eigenvalues(:), fermi
         real(rp) :: omega
      end function grand_potential_from_eigenvalues

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
         type(rotation_result), intent(out) :: result
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
         type(energy), intent(in) :: energy_obj
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

      module subroutine lr_run(this, control_obj, lattice_obj, hamiltonian_obj, energy_obj, self_obj, reciprocal_obj)
         class(linear_response), intent(in) :: this
         type(control), intent(in) :: control_obj
         type(lattice), intent(inout) :: lattice_obj
         type(hamiltonian), intent(inout) :: hamiltonian_obj
         type(energy), intent(in) :: energy_obj
         type(self), intent(inout) :: self_obj
         type(reciprocal), intent(inout) :: reciprocal_obj
      end subroutine lr_run
   end interface

end module linear_response_mod
