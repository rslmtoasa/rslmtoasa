! KXC-01 direct local radial ALSDA kernel oracles.
program test_lr_alsda_kernel
   use precision_mod, only: rp
   use radial_ground_state_mod, only: radial_ground_state, RADIAL_PI
   use lr_response_space_mod, only: response_space_layout, response_apply_operator, response_identity_operator
   use lr_alsda_kernel_mod, only: lr_alsda_kernel_request, lr_alsda_kernel_result, &
      lr_kxc_magnetization_pauli, evaluate_lr_alsda_kernel, evaluate_lr_alsda_static_residual
   implicit none

   integer, parameter :: nr = 9, nsite = 2, lmax = 0, nchannel = 2
   real(rp), parameter :: mesh_a = 0.03_rp, mesh_b = 0.10_rp
   real(rp), parameter :: tolerance = 2.0e-13_rp
   real(rp) :: radius(nr), rho(nr, 2), origin(2), vup(nr), vdn(nr), total_v(nr, 2)
   real(rp) :: pauli_m(nsite, nr), reversed_m(nsite, nr), expected(nsite, nr)
   real(rp) :: field_norm, abs_residual, rel_residual
   complex(rp), allocatable :: field(:), action(:), expected_action(:), identity(:, :), residual(:)
   complex(rp) :: overlap
   type(radial_ground_state), target :: states(nsite), reversed_states(nsite)
   type(response_space_layout), target :: space
   type(lr_alsda_kernel_request) :: request
   type(lr_alsda_kernel_result) :: result, reversed_result, low_m_result
   integer :: ir, isite, flat, i, j
   character(len=32) :: argument
   logical :: failed

   call get_command_argument(1, argument)
   if (trim(argument) == 'provenance') then
      call build_fixture(radius, rho, origin, vup, vdn, total_v, pauli_m)
      call initialize_states(states, radius, rho, origin, vup, vdn, total_v)
      call space%initialize(nsite, lmax, radius, mesh_a, mesh_b, nchannel)
      request%response_space => space
      request%ground_states => states
      allocate(request%pauli_magnetization(nsite, nr))
      request%pauli_magnetization = pauli_m
      request%requested_functional = 'a deliberately mismatched XC functional'
      call evaluate_lr_alsda_kernel(request, result)
      error stop 'provenance mismatch was accepted'
   end if
   if (trim(argument) == 'zero') then
      call build_fixture(radius, rho, origin, vup, vdn, total_v, pauli_m)
      call initialize_states(states, radius, rho, origin, vup, vdn, total_v)
      call space%initialize(nsite, lmax, radius, mesh_a, mesh_b, nchannel)
      request%response_space => space
      request%ground_states => states
      allocate(request%pauli_magnetization(nsite, nr))
      request%pauli_magnetization = pauli_m
      request%pauli_magnetization(1, nr - 1) = 0.0_rp
      call evaluate_lr_alsda_kernel(request, result)
      error stop 'active-grid zero magnetization was accepted'
   end if

   call build_fixture(radius, rho, origin, vup, vdn, total_v, pauli_m)
   call initialize_states(states, radius, rho, origin, vup, vdn, total_v)
   call space%initialize(nsite, lmax, radius, mesh_a, mesh_b, nchannel)

   request%response_space => space
   request%ground_states => states
   allocate(request%pauli_magnetization(nsite, nr))
   request%pauli_magnetization = pauli_m
   request%magnetization_label = lr_kxc_magnetization_pauli
   request%requested_functional = 'Barth-Hedin'
   request%requested_backend = 'legacy RS-LMTO'
   request%requested_txc = 1
   request%low_m_diagnostic_relative = 1.0e-6_rp
   call evaluate_lr_alsda_kernel(request, result)

   failed = .false.
   do isite = 1, nsite
      do ir = 1, nr
         if (pauli_m(isite, ir) /= 0.0_rp) then
            expected(isite, ir) = 0.5_rp*(vup(ir) - vdn(ir))/pauli_m(isite, ir)
         else
            expected(isite, ir) = 0.0_rp
         end if
      end do
   end do
   call check_real_matrix(result%pointwise_kernel, expected, tolerance, 'pointwise LR-01 reconstruction', failed)
   if (trim(result%units) /= 'Ry bohr^3') failed = .true.
   if (trim(result%response_representation) /= 'LR-04 canonical local operator') failed = .true.
   if (trim(result%magnetization_label) /= lr_kxc_magnetization_pauli) failed = .true.
   if (.not. result%origin_null_measure_extension) failed = .true.
   if (result%low_m%active_radial_points /= nsite*(nr - 1)) failed = .true.
   if (result%low_m%low_m_points /= 0) failed = .true.
   if (result%low_m%exact_zero_active_points /= 0) failed = .true.

   ! The canonical local operator is diagonal in site, L, M, radial point,
   ! and channel.  Its action is compared with an independent explicit loop.
   allocate(field(space%ndim), action(space%ndim), expected_action(space%ndim))
   do flat = 1, space%ndim
      field(flat) = cmplx(0.17_rp*real(flat, rp), -0.03_rp*real(flat + 1, rp), rp)
   end do
   call response_apply_operator(space, result%canonical_operator, field, action)
   expected_action = cmplx(0.0_rp, 0.0_rp, rp)
   do flat = 1, space%ndim
      call index_components(flat, isite, ir, i, j)
      expected_action(flat) = result%pointwise_kernel(isite, ir)*field(flat)
   end do
   call check_complex_vector(action, expected_action, tolerance, 'local operator action', failed)

   do i = 1, space%ndim
      do j = 1, space%ndim
         if (i /= j .and. abs(result%canonical_operator(i, j)) > tolerance) failed = .true.
      end do
   end do

   ! Deliberately explicit nested-loop quadrature oracle for the local action.
   ! The local operator is pointwise, so the source metric is already accounted
   ! for by the LR-04 delta representation and no extra W is inserted here.
   expected_action = cmplx(0.0_rp, 0.0_rp, rp)
   do flat = 1, space%ndim
      call index_components(flat, isite, ir, i, j)
      do j = 1, space%ndim
         if (flat == j) expected_action(flat) = expected_action(flat) + &
            result%pointwise_kernel(isite, ir)*field(j)
      end do
   end do
   call check_complex_vector(action, expected_action, tolerance, 'direct local quadrature action', failed)

   ! Global spin reversal swaps both the stored XC channels and the selected
   ! Pauli magnetization.  The mapped ratio must be invariant.
   reversed_m = -pauli_m
   call initialize_reversed_states(reversed_states, states)
   request%ground_states => reversed_states
   request%pauli_magnetization = reversed_m
   call evaluate_lr_alsda_kernel(request, reversed_result)
   call check_real_matrix(reversed_result%pointwise_kernel, result%pointwise_kernel, tolerance, &
                          'global spin-reversal kernel invariance', failed)
   request%ground_states => states
   request%pauli_magnetization = pauli_m

   ! A small but nonzero active magnetization is evaluated exactly, without a
   ! floor or clipping.  It is reported by the diagnostic and remains finite.
   request%pauli_magnetization(nsite, nr - 1) = 1.0e-12_rp
   call evaluate_lr_alsda_kernel(request, low_m_result)
   if (low_m_result%low_m%low_m_points /= 1) failed = .true.
   if (abs(low_m_result%pointwise_kernel(nsite, nr - 1) - &
       0.5_rp*(vup(nr - 1) - vdn(nr - 1))/1.0e-12_rp) > 1.0e-3_rp) failed = .true.
   request%pauli_magnetization = pauli_m

   ! Static Ward/Goldstone diagnostic: the supplied chi and direct kernel are
   ! applied unchanged; the result is compared with an independently supplied
   ! rigid vector and its metric overlap is reported.
   allocate(identity(space%ndim, space%ndim), residual(space%ndim))
   call response_identity_operator(space, identity)
   call evaluate_lr_alsda_static_residual(space, identity, result%canonical_operator, field, &
                                          abs_residual, rel_residual, overlap, residual)
   field_norm = sqrt(sum(abs(field)**2))
   if (abs_residual < 0.0_rp .or. rel_residual < 0.0_rp .or. abs(overlap) < 0.0_rp .or. &
       field_norm <= 0.0_rp .or. maxval(abs(residual)) <= 0.0_rp) then
      failed = .true.
   end if
   write (*, '(a,2(es16.8,1x))') 'KXC raw static residual (absolute, relative) = ', abs_residual, rel_residual
   write (*, '(a,2(es16.8,1x))') 'KXC rigid-vector overlap (real, imag) = ', real(overlap, rp), aimag(overlap)

   if (failed) then
      write (*, '(a)') 'UnitLrAlsdaKernel: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrAlsdaKernel: PASS (pointwise, metric action, angular diagonal, reversal, provenance, low-m)'

contains

   subroutine build_fixture(mesh, density, density_origin, up, down, total, magnetization)
      real(rp), intent(out) :: mesh(:), density(:, :), density_origin(:), up(:), down(:), total(:, :)
      real(rp), intent(out) :: magnetization(:, :)
      real(rp) :: nup, ndn
      integer :: i, site

      mesh(1) = 0.0_rp
      do i = 2, size(mesh)
         mesh(i) = mesh_b*(exp(mesh_a*real(i - 1, rp)) - 1.0_rp)
      end do
      do i = 1, size(mesh)
         nup = 0.9_rp*exp(-mesh(i))
         ndn = 0.5_rp*exp(-0.7_rp*mesh(i))
         density(i, 1) = 4.0_rp*RADIAL_PI*mesh(i)**2*nup
         density(i, 2) = 4.0_rp*RADIAL_PI*mesh(i)**2*ndn
         up(i) = 0.20_rp + 0.03_rp*exp(-mesh(i))
         down(i) = -0.10_rp + 0.01_rp*exp(-0.5_rp*mesh(i))
         total(i, 1) = 1.0_rp + up(i) + 0.07_rp
         total(i, 2) = 1.0_rp + down(i) - 0.07_rp
      end do
      density_origin = [0.9_rp, 0.5_rp]
      do site = 1, size(magnetization, 1)
         do i = 1, size(magnetization, 2)
            magnetization(site, i) = (0.8_rp + 0.2_rp*real(site, rp))*exp(-0.3_rp*mesh(i))
         end do
         magnetization(site, 1) = 0.0_rp
      end do
   end subroutine build_fixture

   subroutine initialize_states(target_states, mesh, density, density_origin, up, down, total)
      type(radial_ground_state), intent(out) :: target_states(:)
      real(rp), intent(in) :: mesh(:), density(:, :), density_origin(:), up(:), down(:), total(:, :)
      integer :: site

      do site = 1, size(target_states)
         call target_states(site)%capture(mesh_a, mesh_b, mesh, density, density_origin, up, down, total, 0.07_rp)
         call target_states(site)%mark_accepted(4, 1.0e-9_rp)
         target_states(site)%xc_provenance%backend_name = 'legacy RS-LMTO'
         target_states(site)%xc_provenance%txch = 'B-H'
         target_states(site)%xc_provenance%functional_name = 'Barth-Hedin'
         target_states(site)%xc_provenance%mapping_quality = 'REFERENCE_EQUIVALENT'
         target_states(site)%xc_provenance%txc = 1
         target_states(site)%xc_provenance%spin_polarized = .true.
      end do
   end subroutine initialize_states

   subroutine initialize_reversed_states(target_states, original_states)
      type(radial_ground_state), intent(out) :: target_states(:)
      type(radial_ground_state), intent(in) :: original_states(:)
      real(rp), allocatable :: density(:, :), total(:, :)
      real(rp) :: origin(2)
      integer :: site, ir

      allocate(density(size(original_states(1)%r), 2), total(size(original_states(1)%r), 2))
      do site = 1, size(target_states)
         density(:, 1) = original_states(site)%rho_weighted_down
         density(:, 2) = original_states(site)%rho_weighted_up
         origin = [original_states(site)%n_down(1), original_states(site)%n_up(1)]
         total(:, 1) = original_states(site)%vxc_down + 0.07_rp
         total(:, 2) = original_states(site)%vxc_up - 0.07_rp
         call target_states(site)%capture(original_states(site)%a, original_states(site)%b, original_states(site)%r, &
            density, origin, original_states(site)%vxc_down, original_states(site)%vxc_up, total, 0.07_rp)
         call target_states(site)%mark_accepted(original_states(site)%accepted_inner_iteration, &
            original_states(site)%accepted_residual_control)
         target_states(site)%xc_provenance = original_states(site)%xc_provenance
      end do
      deallocate(density, total)
   end subroutine initialize_reversed_states

   subroutine index_components(flat, site, radial_point, harmonic, channel)
      integer, intent(in) :: flat
      integer, intent(out) :: site, radial_point, harmonic, channel
      integer :: remaining

      remaining = flat - 1
      channel = mod(remaining, nchannel) + 1
      remaining = remaining/nchannel
      radial_point = mod(remaining, nr) + 1
      remaining = remaining/nr
      harmonic = mod(remaining, (lmax + 1)**2) + 1
      site = remaining/((lmax + 1)**2) + 1
   end subroutine index_components

   subroutine check_real_matrix(actual, expected_matrix, tolerance_in, label, test_failed)
      real(rp), intent(in) :: actual(:, :), expected_matrix(:, :), tolerance_in
      character(len=*), intent(in) :: label
      logical, intent(inout) :: test_failed
      real(rp) :: error

      error = maxval(abs(actual - expected_matrix))
      write (*, '(a,es12.4)') trim(label)//' max error = ', error
      if (error > tolerance_in) test_failed = .true.
   end subroutine check_real_matrix

   subroutine check_complex_vector(actual, expected_vector, tolerance_in, label, test_failed)
      complex(rp), intent(in) :: actual(:), expected_vector(:)
      real(rp), intent(in) :: tolerance_in
      character(len=*), intent(in) :: label
      logical, intent(inout) :: test_failed
      real(rp) :: error

      error = maxval(abs(actual - expected_vector))
      write (*, '(a,es12.4)') trim(label)//' max error = ', error
      if (error > tolerance_in) test_failed = .true.
   end subroutine check_complex_vector

end program test_lr_alsda_kernel
