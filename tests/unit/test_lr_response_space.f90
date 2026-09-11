!------------------------------------------------------------------------------
! LR-04 discrete response-space metric and operator-algebra oracles.
!
! The composition oracle below is deliberately written as nested quadrature
! loops.  It does not call the production matrix-composition helper.
!------------------------------------------------------------------------------
program test_lr_response_space
   use precision_mod, only: rp
   use response_basis_mapping_mod, only: response_super_index, response_unflatten_superindex
   use lr_response_space_mod, only: response_space_layout, response_build_radial_metric, &
      response_vector_inner_product, response_vector_norm, response_apply_operator, &
      response_compose_operators, response_identity_operator, response_local_operator, &
      response_operator_adjoint, response_operator_trace, response_rigid_vector_norm, &
      response_rigid_vector_overlap, response_raw_to_canonical, response_canonical_to_raw, &
      response_metric_weights
   implicit none

   integer, parameter :: nsite = 2, lmax = 1, nr = 7, nchannel = 2
   real(rp), parameter :: a = 0.035_rp, b = 0.08_rp
   real(rp), parameter :: tol = 5.0e-12_rp
   type(response_space_layout) :: space
   real(rp) :: radius(nr), radial_metric(nr), metric_full(nsite*(lmax + 1)**2*nr*nchannel)
   real(rp) :: local_values(nsite, nr), local_channel_values(nsite, nr, nchannel)
   complex(rp), allocatable :: x(:), y(:), z(:), rigid(:), tmp(:)
   complex(rp), allocatable :: identity(:, :), local(:, :), local_channel(:, :), raw_a(:, :), raw_b(:, :)
   complex(rp), allocatable :: canonical_a(:, :), canonical_b(:, :), composed(:, :)
   complex(rp), allocatable :: raw_composed(:, :), oracle(:)
   complex(rp) :: ip, expected_ip, lhs, rhs
   real(rp) :: expected_norm, expected_trace, max_error
   integer :: i, j, k, ir, site, channel, flat, n
   type(response_super_index) :: item_i, item_j
   logical :: failed

   failed = .false.
   radius(1) = 0.0_rp
   do ir = 2, nr
      radius(ir) = b*(exp(a*real(ir - 1, rp)) - 1.0_rp)
   end do
   call space%initialize(nsite, lmax, radius, a, b, nchannel)
   n = space%ndim
   allocate(x(n), y(n), z(n), rigid(n), tmp(n), oracle(n))
   allocate(identity(n, n), local(n, n), local_channel(n, n), raw_a(n, n), raw_b(n, n))
   allocate(canonical_a(n, n), canonical_b(n, n), composed(n, n), raw_composed(n, n))

   call response_build_radial_metric(a, b, radius, radial_metric)
   call response_metric_weights(space, metric_full)
   do i = 1, n
      call response_unflatten_superindex(i, nsite, lmax, nr, nchannel, item_i)
      if (metric_full(i) /= radial_metric(item_i%radial_point)) failed = .true.
   end do
   if (maxval(abs(radial_metric - space%radial_weights)) > 0.0_rp) then
      failed = .true.
      write (*, '(a)') 'Radial metric helper/layout mismatch'
   end if
   if (space%radial_weights(1) /= 0.0_rp .or. space%active_dimension /= nsite*(lmax + 1)**2*(nr - 1)*nchannel) then
      failed = .true.
      write (*, '(a)') 'Origin or active-dimension contract failed'
   end if

   do i = 1, n
      call response_unflatten_superindex(i, nsite, lmax, nr, nchannel, item_i)
      x(i) = cmplx(0.13_rp + 0.011_rp*real(item_i%site, rp) + 0.007_rp*real(item_i%radial_point, rp), &
         -0.02_rp*real(item_i%response_m + item_i%channel, rp), rp)
      y(i) = cmplx(-0.09_rp + 0.003_rp*real(i, rp), 0.04_rp - 0.001_rp*real(item_i%radial_point, rp), rp)
      rigid(i) = cmplx(0.2_rp + 0.01_rp*real(item_i%site, rp), 0.0_rp, rp)
   end do

   ! Analytic radial fields are contracted with the production volume metric,
   ! and independently with an explicit super-index loop.
   ip = response_vector_inner_product(space, x, y)
   expected_ip = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, n
      call response_unflatten_superindex(i, nsite, lmax, nr, nchannel, item_i)
      expected_ip = expected_ip + conjg(x(i))*radial_metric(item_i%radial_point)*y(i)
   end do
   call check_complex(ip, expected_ip, tol, 'analytic radial metric contraction', failed)
   expected_norm = sqrt(max(0.0_rp, real(sum(conjg(x)*space%metric_weights*x), rp)))
   call check_real(response_vector_norm(space, x), expected_norm, tol, 'metric vector norm', failed)

   call response_identity_operator(space, identity)
   call response_apply_operator(space, identity, x, z)
   call check_vector(z, x, tol, 'identity action', failed)
   call check_complex(response_operator_trace(space, identity), cmplx(real(space%active_dimension, rp), 0.0_rp, rp), &
      tol, 'positive-measure identity trace', failed)

   do site = 1, nsite
      do ir = 1, nr
         local_values(site, ir) = 0.17_rp*real(site, rp) - 0.03_rp*real(ir, rp) + 0.005_rp*radius(ir)
      end do
   end do
   call response_local_operator(space, local_values, local)
   call response_apply_operator(space, local, x, z)
   do i = 1, n
      call response_unflatten_superindex(i, nsite, lmax, nr, nchannel, item_i)
      tmp(i) = cmplx(local_values(item_i%site, item_i%radial_point), 0.0_rp, rp)*x(i)
   end do
   call check_vector(z, tmp, tol, 'pointwise local multiplication', failed)
   do i = 1, n
      call response_unflatten_superindex(i, nsite, lmax, nr, nchannel, item_i)
      do j = 1, n
         call response_unflatten_superindex(j, nsite, lmax, nr, nchannel, item_j)
         if (item_i%site /= item_j%site .and. abs(local(i, j)) > 0.0_rp) failed = .true.
      end do
   end do

   do site = 1, nsite
      do ir = 1, nr
         do channel = 1, nchannel
            local_channel_values(site, ir, channel) = local_values(site, ir) + 0.011_rp*real(channel, rp)
         end do
      end do
   end do
   call response_local_operator(space, local_channel_values, local_channel)
   call response_apply_operator(space, local_channel, x, z)
   do i = 1, n
      call response_unflatten_superindex(i, nsite, lmax, nr, nchannel, item_i)
      tmp(i) = cmplx(local_channel_values(item_i%site, item_i%radial_point, item_i%channel), 0.0_rp, rp)*x(i)
   end do
   call check_vector(z, tmp, tol, 'channel-resolved local multiplication', failed)

   ! Two site-separated nonlocal kernels are converted once and composed by
   ! the production helper.  The comparison uses independent nested loops.
   raw_a = cmplx(0.0_rp, 0.0_rp, rp)
   raw_b = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, n
      call response_unflatten_superindex(i, nsite, lmax, nr, nchannel, item_i)
      do j = 1, n
         call response_unflatten_superindex(j, nsite, lmax, nr, nchannel, item_j)
         if (item_i%site == item_j%site .and. item_j%radial_point > 1) then
            raw_a(i, j) = cmplx(0.001_rp*real(i + 2*j, rp), 0.0003_rp*real(j - i, rp), rp)
            raw_b(i, j) = cmplx(-0.0007_rp*real(2*i + j, rp), 0.0002_rp*real(i + j, rp), rp)
         end if
      end do
   end do
   call response_raw_to_canonical(space, raw_a, canonical_a)
   call response_raw_to_canonical(space, raw_b, canonical_b)
   call response_compose_operators(space, canonical_a, canonical_b, composed)
   call response_canonical_to_raw(space, canonical_a, raw_composed)
   call check_matrix(raw_composed, raw_a, tol, 'raw/canonical active-column round trip', failed)

   oracle = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, n
      do j = 1, n
         do k = 1, n
            oracle(i) = oracle(i) + raw_a(i, j)*space%metric_weights(j)*raw_b(j, k)* &
               space%metric_weights(k)*x(k)
         end do
      end do
   end do
   call response_apply_operator(space, composed, x, z)
   call check_vector(z, oracle, 10.0_rp*tol, 'nested quadrature composition oracle', failed)

   call response_operator_adjoint(space, canonical_a, composed)
   call response_apply_operator(space, canonical_a, y, z)
   lhs = response_vector_inner_product(space, x, z)
   call response_apply_operator(space, composed, x, tmp)
   rhs = response_vector_inner_product(space, tmp, y)
   call check_complex(lhs, rhs, 20.0_rp*tol, 'metric adjoint relation', failed)

   ! The mapping is site-major and the local operator cannot mix sites or
   ! response harmonics.  This checks multisite separation through the live
   ! LR-02R super-index rather than a parallel index implementation.
   do i = 1, n
      call response_unflatten_superindex(i, nsite, lmax, nr, nchannel, item_i)
      do j = 1, n
         call response_unflatten_superindex(j, nsite, lmax, nr, nchannel, item_j)
         if (item_i%site /= item_j%site .and. abs(canonical_a(i, j)) > 0.0_rp) failed = .true.
      end do
   end do

   expected_trace = 0.0_rp
   do flat = 1, n
      call response_unflatten_superindex(flat, nsite, lmax, nr, nchannel, item_i)
      if (item_i%radial_point > 1) expected_trace = expected_trace + local_values(item_i%site, item_i%radial_point)
   end do
   call check_complex(response_operator_trace(space, local), cmplx(expected_trace, 0.0_rp, rp), tol, &
      'local positive-measure trace', failed)

   call check_real(response_rigid_vector_norm(space, rigid), response_vector_norm(space, rigid), tol, &
      'rigid-vector metric norm', failed)
   call check_complex(response_rigid_vector_overlap(space, x, rigid), &
      response_vector_inner_product(space, x, rigid), tol, 'rigid-vector metric overlap', failed)

   max_error = maxval(abs(space%metric_weights(1:nchannel)))
   write (*, '(a,i0,a,es12.4)') 'LR-04 response dimension = ', n, ', origin metric maximum = ', max_error
   if (failed) then
      write (*, '(a)') 'UnitLrResponseSpace: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrResponseSpace: PASS (metric, origin, local, composition, adjoint, multisite, rigid-vector oracles)'

contains

   subroutine check_complex(actual, expected, tolerance, label, test_failed)
      complex(rp), intent(in) :: actual, expected
      real(rp), intent(in) :: tolerance
      character(len=*), intent(in) :: label
      logical, intent(inout) :: test_failed
      real(rp) :: error

      error = abs(actual - expected)
      write (*, '(a,es12.4)') trim(label)//' error = ', error
      if (error > tolerance) test_failed = .true.
   end subroutine check_complex

   subroutine check_real(actual, expected, tolerance, label, test_failed)
      real(rp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      logical, intent(inout) :: test_failed

      write (*, '(a,es12.4)') trim(label)//' error = ', abs(actual - expected)
      if (abs(actual - expected) > tolerance) test_failed = .true.
   end subroutine check_real

   subroutine check_vector(actual, expected, tolerance, label, test_failed)
      complex(rp), intent(in) :: actual(:), expected(:)
      real(rp), intent(in) :: tolerance
      character(len=*), intent(in) :: label
      logical, intent(inout) :: test_failed
      real(rp) :: error

      error = maxval(abs(actual - expected))
      write (*, '(a,es12.4)') trim(label)//' max error = ', error
      if (error > tolerance) test_failed = .true.
   end subroutine check_vector

   subroutine check_matrix(actual, expected, tolerance, label, test_failed)
      complex(rp), intent(in) :: actual(:, :), expected(:, :)
      real(rp), intent(in) :: tolerance
      character(len=*), intent(in) :: label
      logical, intent(inout) :: test_failed
      real(rp) :: error

      error = maxval(abs(actual - expected))
      write (*, '(a,es12.4)') trim(label)//' max error = ', error
      if (error > tolerance) test_failed = .true.
   end subroutine check_matrix

end program test_lr_response_space
