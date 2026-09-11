! TDRUN-01 production-adapter contract and tiny prepared-state integration test.
program test_tddft_production_driver
   use precision_mod, only: rp
   use tddft_production_driver_mod, only: tddft_production_config, tddft_capability_state, &
      tddft_capability_is_supported, require_tddft_capability, load_tddft_config, &
      evaluate_tddft_production_sweep, tddft_production_result
   use radial_ground_state_mod, only: radial_ground_state
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_response_space_mod, only: response_space_layout
   use lr_ks_susceptibility_mod, only: lr_electronic_state, lr_ks_susceptibility_request, &
      lr_ks_susceptibility_result, evaluate_lr_ks_susceptibility
   use lr_alsda_kernel_mod, only: lr_alsda_kernel_request, lr_alsda_kernel_result, evaluate_lr_alsda_kernel
   use tddft_dyson_mod, only: tddft_dyson_request, tddft_dyson_result, evaluate_tddft_dyson
   implicit none

   character(len=32) :: argument
   type(tddft_capability_state) :: capability
   type(tddft_production_config) :: config
   type(radial_ground_state), target :: ground(1)
   type(lmto_radial_basis), target :: radial(1)
   type(response_space_layout), target :: space
   type(lr_electronic_state), target :: state, endpoint
   type(lr_ks_susceptibility_request) :: ks_request
   type(lr_ks_susceptibility_result) :: ks_direct
   type(lr_alsda_kernel_request) :: kxc_request
   type(lr_alsda_kernel_result) :: kxc_direct
   type(tddft_dyson_request) :: dyson_request
   type(tddft_dyson_result) :: dyson_direct
   type(tddft_production_result) :: driver_result
   real(rp) :: mesh(5), evals(2, 1), weights(1), kpoints(3, 1), occupations(2, 1)
   complex(rp) :: vectors(8, 2, 1)
   real(rp), allocatable :: saved_vup(:), saved_vdn(:), saved_nup(:), saved_ndn(:), saved_m(:)
   logical :: ok
   character(len=128) :: reason

   call get_command_argument(1, argument)
   if (trim(argument) == 'soc') then
      capability%has_soc = .true.
      call require_tddft_capability(capability)
   else if (trim(argument) == 'overlap') then
      capability%generalized_overlap = .true.
      capability%orthogonal = .false.
      capability%reciprocal_mode = 'generalized_overlap_proxy'
      call require_tddft_capability(capability)
   end if

   call load_tddft_config('this-file-is-not-a-production-input', config)
   if (config%enabled .or. config%present) error stop 'feature-off TDDFT parser default changed'

   ok = tddft_capability_is_supported(capability, reason)
   if (.not. ok) error stop 'baseline capability unexpectedly rejected'

   call build_fixture(mesh, evals, vectors, weights, kpoints, occupations, ground, radial, space, state, endpoint)
   call config%restore_to_default()
   config%present = .true.
   config%enabled = .true.
   config%nq = 1
   config%nfrequency = 1
   config%q_list(:, 1) = 0.0_rp
   config%frequencies(1) = 0.17_rp
   config%eta = 0.03_rp
   config%response_lmax = 2
   config%interaction_route = 'direct_alsda'
   config%backend = 'lehmann'

   allocate(saved_vup(size(ground(1)%vxc_up)), saved_vdn(size(ground(1)%vxc_down)), &
            saved_nup(size(ground(1)%n_up)), saved_ndn(size(ground(1)%n_down)), saved_m(size(ground(1)%n_up)))
   saved_vup = ground(1)%vxc_up
   saved_vdn = ground(1)%vxc_down
   saved_nup = ground(1)%n_up
   saved_ndn = ground(1)%n_down
   saved_m = ground(1)%n_up - ground(1)%n_down

   ! Direct validated service calls for the same tiny accepted-state fixture.
   ks_request%q = config%q_list(:, 1)
   ks_request%frequencies = config%frequencies
   ks_request%eta = config%eta
   ks_request%channel = config%channel
   ks_request%response_space => space
   ks_request%radial_bases => radial
   ks_request%electronic_state => state
   ks_request%q_endpoint_state => endpoint
   call evaluate_lr_ks_susceptibility(ks_request, ks_direct)

   kxc_request%response_space => space
   kxc_request%ground_states => ground
   allocate(kxc_request%pauli_magnetization(1, space%npoint))
   kxc_request%pauli_magnetization(1, :) = saved_m
   call evaluate_lr_alsda_kernel(kxc_request, kxc_direct)

   dyson_request%response_space => space
   dyson_request%q = config%q_list(:, 1)
   dyson_request%frequencies = config%frequencies
   dyson_request%eta = config%eta
   dyson_request%channel = config%channel
   dyson_request%ks_susceptibility = ks_direct%susceptibility
   dyson_request%canonical_interaction = kxc_direct%canonical_operator
   dyson_request%interaction_route = 'direct_alsda'
   dyson_request%interaction_provenance = 'KXC-01 direct ALSDA LR-03'
   dyson_request%electronic_state_provenance = 'accepted reciprocal eigenpair snapshot; occupations and EF fixed at SCF handoff'
   dyson_request%response_space_metadata = ks_direct%response_space_metadata
   call evaluate_tddft_dyson(dyson_request, dyson_direct)

   call evaluate_tddft_production_sweep(config, space, radial, ground, state, [endpoint], driver_result)
   if (maxval(abs(driver_result%ks_susceptibility(:, :, :, 1) - dyson_direct%ks_susceptibility)) > 2.0e-12_rp) then
      error stop 'production driver and direct KS response disagree'
   end if
   if (maxval(abs(driver_result%enhanced_susceptibility(:, :, :, 1) - dyson_direct%enhanced_susceptibility)) > 2.0e-12_rp) then
      error stop 'production driver and direct Dyson response disagree'
   end if
   if (maxval(abs(driver_result%loss_matrix(:, :, :, 1) - dyson_direct%loss_matrix)) > 2.0e-12_rp) then
      error stop 'production driver and direct loss response disagree'
   end if
   if (maxval(abs(ground(1)%vxc_up - saved_vup)) > 0.0_rp .or. maxval(abs(ground(1)%vxc_down - saved_vdn)) > 0.0_rp .or. &
       maxval(abs(ground(1)%n_up - saved_nup)) > 0.0_rp .or. maxval(abs(ground(1)%n_down - saved_ndn)) > 0.0_rp) then
      error stop 'production response mutated accepted radial state'
   end if
   write(*, '(a)') 'UnitTddftProductionDriver: PASS (feature-off, direct reproducibility, no state mutation)'

contains

   subroutine build_fixture(mesh, evals, vecs, k_weights, points, occ, ground, radial, space, state, endpoint)
      real(rp), intent(out) :: mesh(:), evals(:, :), k_weights(:), points(:, :), occ(:, :)
      complex(rp), intent(out) :: vecs(:, :, :)
      type(radial_ground_state), intent(out) :: ground(:)
      type(lmto_radial_basis), intent(out) :: radial(:)
      type(response_space_layout), intent(out) :: space
      type(lr_electronic_state), intent(out) :: state, endpoint
      integer :: ir

      mesh(1) = 0.0_rp
      do ir = 2, size(mesh)
         mesh(ir) = 0.10_rp*(exp(0.03_rp*real(ir - 1, rp)) - 1.0_rp)
      end do
      call ground(1)%clear()
      ground(1)%a = 0.03_rp
      ground(1)%b = 0.10_rp
      ground(1)%rmax = mesh(size(mesh))
      allocate(ground(1)%r(size(mesh)), ground(1)%n_up(size(mesh)), ground(1)%n_down(size(mesh)), &
               ground(1)%vxc_up(size(mesh)), ground(1)%vxc_down(size(mesh)))
      ground(1)%r = mesh
      ground(1)%n_up = 0.8_rp
      ground(1)%n_down = 0.2_rp
      ground(1)%vxc_up = 0.20_rp
      ground(1)%vxc_down = -0.10_rp
      ground(1)%valid = .true.
      ground(1)%accepted = .true.
      ground(1)%xc_provenance%functional_name = 'fixture-xc'
      ground(1)%xc_provenance%backend_name = 'fixture-backend'
      ground(1)%xc_provenance%internal_energy_units = 'Ry'
      ground(1)%xc_provenance%txc = 1
      ground(1)%xc_provenance%spin_polarized = .true.
      call ground(1)%begin_pauli_basis(1)
      ground(1)%pauli_large(:, :, :) = 1.0_rp
      ground(1)%pauli_basis_valid = .true.

      call radial(1)%initialize(size(mesh), 1, 2)
      radial(1)%rofi = mesh
      radial(1)%mesh_a = 0.03_rp
      radial(1)%mesh_b = 0.10_rp
      radial(1)%channel_present = .true.
      radial(1)%phi_large(:, :, :) = 1.0_rp
      radial(1)%phidot_large = 0.0_rp
      radial(1)%enu_work = 0.0_rp
      call space%initialize(1, 2, mesh, 0.03_rp, 0.10_rp, 1)

      evals(:, 1) = [-0.20_rp, 0.20_rp]
      vecs = cmplx(0.0_rp, 0.0_rp, rp)
      vecs(1, 1, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      vecs(5, 2, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      k_weights = 1.0_rp
      points = 0.0_rp
      occ(:, 1) = [1.0_rp, 0.0_rp]
      call state%initialize(evals, vecs, points, k_weights, occ, 0.0_rp, 0.0_rp)
      call endpoint%initialize(evals, vecs, points, k_weights, occ, 0.0_rp, 0.0_rp)
   end subroutine build_fixture

end program test_tddft_production_driver
