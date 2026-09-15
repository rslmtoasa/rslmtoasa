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
   use lr_gf_susceptibility_mod, only: lr_gf_susceptibility_request, evaluate_lr_gf_susceptibility
   use lr_alsda_kernel_mod, only: lr_alsda_kernel_request, lr_alsda_kernel_result, evaluate_lr_alsda_kernel
   use lr_rs_gf_susceptibility_mod, only: lr_rs_dense_gf_provider, lr_rs_gf_pair, &
      lr_rs_gf_susceptibility_request, evaluate_lr_rs_gf_susceptibility
   use tddft_dyson_mod, only: tddft_dyson_request, tddft_dyson_result, evaluate_tddft_dyson
   implicit none

   character(len=32) :: argument
   type(tddft_capability_state) :: capability
   type(tddft_production_config) :: config
   type(tddft_production_config) :: parsed_config
   type(radial_ground_state), target :: ground(1)
   type(lmto_radial_basis), target :: radial(1)
   type(response_space_layout), target :: space
   type(lr_electronic_state), target :: state, endpoint
   type(lr_ks_susceptibility_request) :: ks_request
   type(lr_ks_susceptibility_result) :: ks_direct
   type(lr_gf_susceptibility_request) :: gf_request
   type(lr_ks_susceptibility_result) :: gf_direct
   type(lr_alsda_kernel_request) :: kxc_request
   type(lr_alsda_kernel_result) :: kxc_direct
   type(tddft_dyson_request) :: dyson_request
   type(tddft_dyson_result) :: dyson_direct
   type(tddft_production_result) :: driver_result
   type(tddft_production_result) :: crosscheck_result
   type(tddft_production_result) :: native_driver_result
   type(lr_rs_dense_gf_provider), target :: native_provider
   complex(rp) :: native_hamiltonian(8, 8), native_hgamma(8, 8)
   type(lr_rs_gf_susceptibility_request) :: native_request
   type(lr_ks_susceptibility_result) :: native_direct
   type(lr_rs_gf_pair) :: native_pairs(1)
   real(rp) :: mesh(5), evals(2, 1), weights(1), kpoints(3, 1), occupations(2, 1)
   complex(rp) :: vectors(8, 2, 1)
   real(rp), allocatable :: saved_vup(:), saved_vdn(:), saved_nup(:), saved_ndn(:), saved_m(:)
   complex(rp), allocatable :: expected_delta(:, :, :)
   real(rp) :: expected_norm_lehmann, expected_norm_gf, expected_difference_frobenius
   real(rp) :: expected_relative_frobenius, expected_difference_infinity
   logical :: ok
   integer :: i, parser_unit, parser_ios
   character(len=128) :: parser_fixture
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
   else if (trim(argument) == 'static-no-gamma') then
      parser_fixture = 'tdv_static_no_gamma.nml'
      open(newunit=parser_unit, file=parser_fixture, status='replace', action='write', iostat=parser_ios)
      if (parser_ios /= 0) error stop 'could not create static-audit negative fixture'
      write(parser_unit, '(a)') '&tddft'
      write(parser_unit, '(a)') ' enabled = .true., backend = ''compact_dyson'', dyson_static_audit = .true.'
      write(parser_unit, '(a)') ' n_q = 1, q_list = 0.03, 0.0, 0.0'
      write(parser_unit, '(a)') '/'
      close(parser_unit)
      call load_tddft_config(parser_fixture, parsed_config)
      error stop 'static-audit preflight unexpectedly accepted a q list without Gamma'
   else if (trim(argument) == 'covariance-no-pair') then
      parser_fixture = 'tdv_covariance_no_pair.nml'
      open(newunit=parser_unit, file=parser_fixture, status='replace', action='write', iostat=parser_ios)
      if (parser_ios /= 0) error stop 'could not create covariance negative fixture'
      write(parser_unit, '(a)') '&tddft'
      write(parser_unit, '(a)') ' enabled = .true., backend = ''compact_dyson'', validate_interacting_covariance = .true.'
      write(parser_unit, '(a)') ' n_q = 1, q_list = 0.03, 0.0, 0.0'
      write(parser_unit, '(a)') '/'
      close(parser_unit)
      call load_tddft_config(parser_fixture, parsed_config)
      error stop 'covariance preflight unexpectedly accepted a q list without -q'
   end if

   call load_tddft_config('this-file-is-not-a-production-input', config)
   if (config%enabled .or. config%present) error stop 'feature-off TDDFT parser default changed'

   parser_fixture = 'tdv_k01_parser_default.nml'
   open(newunit=parser_unit, file=parser_fixture, status='replace', action='write', iostat=parser_ios)
   if (parser_ios /= 0) error stop 'could not create TDVK-01 parser fixture'
   write(parser_unit, '(a)') '&tddft'
   write(parser_unit, '(a)') ' enabled = .true.'
   write(parser_unit, '(a)') '/'
   close(parser_unit)
   call load_tddft_config(parser_fixture, parsed_config)
   open(newunit=parser_unit, file=parser_fixture, status='old', iostat=parser_ios)
   if (parser_ios == 0) close(parser_unit, status='delete')
   if (.not. parsed_config%present .or. .not. parsed_config%enabled .or. parsed_config%reciprocal_backend_crosscheck .or. &
       parsed_config%gf_closure_audit .or. parsed_config%dyson_static_audit .or. &
       parsed_config%validate_interacting_covariance .or. .not. allocated(parsed_config%eta_values) .or. &
       size(parsed_config%eta_values) /= 1) then
      error stop 'optional TDDFT audit flags did not default to false'
   end if
   open(newunit=parser_unit, file=parser_fixture, status='replace', action='write', iostat=parser_ios)
   if (parser_ios /= 0) error stop 'could not create TDVK-03 parser fixture'
   write(parser_unit, '(a)') '&tddft'
   write(parser_unit, '(a)') ' enabled = .true., gf_closure_audit = .true.'
   write(parser_unit, '(a)') '/'
   close(parser_unit)
   call load_tddft_config(parser_fixture, parsed_config)
   open(newunit=parser_unit, file=parser_fixture, status='old', iostat=parser_ios)
   if (parser_ios == 0) close(parser_unit, status='delete')
   if (.not. parsed_config%gf_closure_audit) error stop 'TDVK-03 GF closure audit flag did not parse'

   open(newunit=parser_unit, file=parser_fixture, status='replace', action='write', iostat=parser_ios)
   if (parser_ios /= 0) error stop 'could not create compact-audit parser fixture'
   write(parser_unit, '(a)') '&tddft'
   write(parser_unit, '(a)') ' enabled = .true., backend = ''compact_dyson'', dyson_static_audit = .true., validate_interacting_covariance = .true.'
   write(parser_unit, '(a)') ' n_q = 3, q_list = 0.0, 0.0, 0.0, 0.03, 0.0, 0.0, -0.03, 0.0, 0.0'
   write(parser_unit, '(a)') '/'
   close(parser_unit)
   call load_tddft_config(parser_fixture, parsed_config)
   open(newunit=parser_unit, file=parser_fixture, status='old', iostat=parser_ios)
   if (parser_ios == 0) close(parser_unit, status='delete')
   if (.not. parsed_config%dyson_static_audit .or. .not. parsed_config%validate_interacting_covariance) then
      error stop 'compact Dyson validation switches did not parse'
   end if

   open(newunit=parser_unit, file=parser_fixture, status='replace', action='write', iostat=parser_ios)
   if (parser_ios /= 0) error stop 'could not create TDVK-05 eta-ladder parser fixture'
   write(parser_unit, '(a)') '&tddft'
   write(parser_unit, '(a)') ' enabled = .true., backend = ''product_convergence'''
   write(parser_unit, '(a)') ' n_q = 2, q_list = 0.0, 0.0, 0.0, 0.125, 0.0, 0.0'
   write(parser_unit, '(a)') ' n_omega = 2, use_omega_grid = .true., omega_grid = 0.0, 0.02'
   write(parser_unit, '(a)') ' n_eta = 3, eta_grid = 0.02, 0.01, 0.005'
   write(parser_unit, '(a)') '/'
   close(parser_unit)
   call load_tddft_config(parser_fixture, parsed_config)
   open(newunit=parser_unit, file=parser_fixture, status='old', iostat=parser_ios)
   if (parser_ios == 0) close(parser_unit, status='delete')
   if (.not. allocated(parsed_config%eta_values) .or. size(parsed_config%eta_values) /= 3 .or. &
       maxval(abs(parsed_config%eta_values - [0.02_rp, 0.01_rp, 0.005_rp])) > 1.0e-14_rp .or. &
       abs(parsed_config%eta - 0.02_rp) > 1.0e-14_rp) then
      error stop 'TDVK-05 physical eta ladder did not parse'
   end if

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
   config%reciprocal_backend_crosscheck = .false.
   config%gf_integration_points = 101
   config%gf_integration_eta = 0.005_rp
   config%gf_energy_margin = 1.0_rp

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

   gf_request%q = config%q_list(:, 1)
   gf_request%frequencies = config%frequencies
   gf_request%eta = config%eta
   gf_request%channel = config%channel
   gf_request%integration_points = config%gf_integration_points
   gf_request%integration_eta = config%gf_integration_eta
   gf_request%energy_margin = config%gf_energy_margin
   gf_request%response_space => space
   gf_request%radial_bases => radial
   gf_request%electronic_state => state
   gf_request%q_endpoint_state => endpoint
   call evaluate_lr_gf_susceptibility(gf_request, gf_direct)

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

   if (driver_result%reciprocal_backend_crosscheck) error stop 'crosscheck was enabled in the flag-off fixture'
   if (allocated(driver_result%reciprocal_crosscheck_valid)) then
      error stop 'flag-off production result unexpectedly contains crosscheck diagnostics'
   end if

   ! The opt-in diagnostic evaluates both reciprocal services but must leave
   ! the selected Lehmann production path, including Dyson and loss, alone.
   config%reciprocal_backend_crosscheck = .true.
   call evaluate_tddft_production_sweep(config, space, radial, ground, state, [endpoint], crosscheck_result)
   if (.not. crosscheck_result%reciprocal_backend_crosscheck .or. &
       .not. all(crosscheck_result%reciprocal_crosscheck_valid)) then
      error stop 'enabled reciprocal backend crosscheck did not produce valid diagnostics'
   end if
   if (any(crosscheck_result%reciprocal_crosscheck_norm_lehmann < 0.0_rp) .or. &
       any(crosscheck_result%reciprocal_crosscheck_norm_gf < 0.0_rp) .or. &
       any(crosscheck_result%reciprocal_crosscheck_difference_frobenius < 0.0_rp) .or. &
       any(crosscheck_result%reciprocal_crosscheck_relative_frobenius < 0.0_rp) .or. &
       any(crosscheck_result%reciprocal_crosscheck_difference_infinity < 0.0_rp)) then
      error stop 'reciprocal backend crosscheck diagnostics are negative'
   end if
   if (any(crosscheck_result%reciprocal_crosscheck_norm_lehmann /= crosscheck_result%reciprocal_crosscheck_norm_lehmann) .or. &
       any(crosscheck_result%reciprocal_crosscheck_norm_gf /= crosscheck_result%reciprocal_crosscheck_norm_gf) .or. &
       any(crosscheck_result%reciprocal_crosscheck_difference_frobenius /= &
          crosscheck_result%reciprocal_crosscheck_difference_frobenius) .or. &
       any(crosscheck_result%reciprocal_crosscheck_relative_frobenius /= &
          crosscheck_result%reciprocal_crosscheck_relative_frobenius) .or. &
       any(crosscheck_result%reciprocal_crosscheck_difference_infinity /= &
          crosscheck_result%reciprocal_crosscheck_difference_infinity)) then
      error stop 'reciprocal backend crosscheck diagnostics are not finite'
   end if
   if (maxval(abs(crosscheck_result%ks_susceptibility - driver_result%ks_susceptibility)) > 0.0_rp .or. &
       maxval(abs(crosscheck_result%enhanced_susceptibility - driver_result%enhanced_susceptibility)) > 0.0_rp .or. &
       maxval(abs(crosscheck_result%loss_matrix - driver_result%loss_matrix)) > 0.0_rp) then
      error stop 'enabling reciprocal backend crosscheck changed the selected production result'
   end if
   allocate(expected_delta(space%ndim, space%ndim, size(config%frequencies)))
   expected_delta = ks_direct%susceptibility - gf_direct%susceptibility
   if (maxval(abs(crosscheck_result%reciprocal_crosscheck_delta(:, :, :, 1) - expected_delta)) > 0.0_rp) then
      error stop 'crosscheck delta does not match independently evaluated reciprocal services'
   end if
   expected_norm_lehmann = sqrt(sum(abs(ks_direct%susceptibility(:, :, 1))**2))
   expected_norm_gf = sqrt(sum(abs(gf_direct%susceptibility(:, :, 1))**2))
   expected_difference_frobenius = sqrt(sum(abs(expected_delta(:, :, 1))**2))
   expected_relative_frobenius = expected_difference_frobenius/max(expected_norm_lehmann, expected_norm_gf, tiny(1.0_rp))
   expected_difference_infinity = maxval(abs(expected_delta(:, :, 1)))
   if (abs(crosscheck_result%reciprocal_crosscheck_norm_lehmann(1, 1) - expected_norm_lehmann) > 0.0_rp .or. &
       abs(crosscheck_result%reciprocal_crosscheck_norm_gf(1, 1) - expected_norm_gf) > 0.0_rp .or. &
       abs(crosscheck_result%reciprocal_crosscheck_difference_frobenius(1, 1) - expected_difference_frobenius) > 0.0_rp .or. &
       abs(crosscheck_result%reciprocal_crosscheck_relative_frobenius(1, 1) - expected_relative_frobenius) > 0.0_rp .or. &
       abs(crosscheck_result%reciprocal_crosscheck_difference_infinity(1, 1) - expected_difference_infinity) > 0.0_rp) then
      error stop 'crosscheck metrics do not match independent norm calculations'
   end if
   write(*, '(a,5(1x,es14.6))') 'TDVK-01 example crosscheck (normL normGF dF relF dInf) =', &
      crosscheck_result%reciprocal_crosscheck_norm_lehmann(1, 1), crosscheck_result%reciprocal_crosscheck_norm_gf(1, 1), &
      crosscheck_result%reciprocal_crosscheck_difference_frobenius(1, 1), &
      crosscheck_result%reciprocal_crosscheck_relative_frobenius(1, 1), &
      crosscheck_result%reciprocal_crosscheck_difference_infinity(1, 1)

   ! R3 reproducibility seam: the native production route must forward the
   ! same prepared request to the public native-RSGF service.  The fixture
   ! provider is deliberately independent of reciprocal GF code.
   config%backend = 'native_rsgf'
   config%native_rsgf_provider = 'auto'
   config%gf_integration_points = 3
   config%gf_integration_eta = 0.005_rp
   config%gf_energy_margin = 1.0_rp
   native_hamiltonian = cmplx(0.0_rp, 0.0_rp, rp)
   native_hgamma = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, 8
      native_hamiltonian(i, i) = cmplx(0.05_rp*real(i, rp), 0.0_rp, rp)
      native_hgamma(i, i) = cmplx(0.10_rp, 0.0_rp, rp)
   end do
   call native_provider%initialize(native_hamiltonian, native_hgamma, 8, .true.)
   native_pairs(1)%left_site = 1
   native_pairs(1)%right_site = 1
   native_pairs(1)%translation = 0.0_rp
   native_request%q = config%q_list(:, 1)
   native_request%frequencies = config%frequencies
   native_request%eta = config%eta
   native_request%channel = config%channel
   native_request%integration_points = config%gf_integration_points
   native_request%integration_eta = config%gf_integration_eta
   native_request%energy_margin = config%gf_energy_margin
   native_request%energy_min = minval(state%eigenvalues)
   native_request%energy_max = maxval(state%eigenvalues)
   native_request%fermi_level = state%fermi_level
   native_request%temperature = state%temperature
   native_request%response_space => space
   native_request%radial_bases => radial
   native_request%provider => native_provider
   native_request%pairs = native_pairs
   call evaluate_lr_rs_gf_susceptibility(native_request, native_direct)

   call evaluate_tddft_production_sweep(config, space, radial, ground, state, [endpoint], native_driver_result, &
      native_provider, native_pairs)
   if (maxval(abs(native_driver_result%ks_susceptibility(:, :, :, 1) - native_direct%susceptibility)) > 2.0e-12_rp) then
      error stop 'native production driver and direct RSGF service disagree'
   end if
   write(*, '(a)') 'UnitTddftProductionDriver: PASS (feature-off, direct reproducibility, native RSGF/service equivalence, no state mutation)'

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
