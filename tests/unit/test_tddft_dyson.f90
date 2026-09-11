!------------------------------------------------------------------------------
! TDDY-01 -- enhanced susceptibility and LR-04 loss matrix
!------------------------------------------------------------------------------
program test_tddft_dyson
   use precision_mod, only: rp
   use tddft_dyson_mod, only: tddft_dyson_request, tddft_dyson_result, &
      evaluate_tddft_dyson, solve_tddft_dyson_frequency, tddft_loss_matrix, &
      loss_matrix_hermiticity_residual, lr_dyson_route_direct_alsda, &
      lr_dyson_route_goldstone_sumrule, lr_dyson_route_direct_alsda_goldstone_corrected
   use lr_response_space_mod, only: response_space_layout
   implicit none

   integer, parameter :: nr = 5
   real(rp), parameter :: mesh_a = 0.04_rp, mesh_b = 0.08_rp
   real(rp), parameter :: tol = 2.0e-10_rp
   type(response_space_layout), target :: space
   real(rp) :: radius(nr)
   logical :: failed
   integer :: ir

   failed = .false.
   radius(1) = 0.0_rp
   do ir = 2, nr
      radius(ir) = mesh_b*(exp(mesh_a*real(ir - 1, rp)) - 1.0_rp)
   end do
   call space%initialize(1, 0, radius, mesh_a, mesh_b, 1)

   call test_loss_convention(failed)
   call test_scalar_dyson(failed)
   call test_small_noncommuting_dyson(failed)
   call test_metric_sensitive_dyson(failed)
   call test_goldstone_near_pole(failed)
   call test_covariance(failed)
   call test_interaction_separation(failed)

   if (failed) then
      write (*, '(a)') 'UnitTddftDyson: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitTddftDyson: PASS (Dyson, metric, loss, pole, covariance, route separation)'

contains

   subroutine test_loss_convention(test_failed)
      logical, intent(inout) :: test_failed
      complex(rp) :: chi(2, 2), loss(2, 2)
      complex(rp), allocatable :: canonical(:, :), metric_loss(:, :)
      integer :: i, j

      chi = cmplx(0.0_rp, 0.0_rp, rp)
      chi(1, 1) = cmplx(1.0_rp, -2.0_rp, rp)
      chi(2, 2) = cmplx(2.0_rp, -1.0_rp, rp)
      chi(1, 2) = cmplx(0.3_rp, 0.4_rp, rp)
      chi(2, 1) = cmplx(-0.2_rp, 0.1_rp, rp)
      loss = tddft_loss_matrix(chi)
      call check('ordinary loss is Hermitian', loss_matrix_hermiticity_residual(loss) < tol, test_failed)
      call check('ordinary positive-frequency diagonal sign', real(loss(1, 1), rp) > 0.0_rp .and. &
         real(loss(2, 2), rp) > 0.0_rp, test_failed)

      ! A canonical susceptibility is a raw rank-one response multiplied by
      ! the source-column metric.  The metric overload must recover the same
      ! physical loss without an extra W.
      allocate(canonical(space%ndim, space%ndim), metric_loss(space%ndim, space%ndim))
      canonical = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 2, space%ndim
         do j = 2, space%ndim
            canonical(i, j) = chi(1, 1)*cmplx(0.1_rp*real(i - 1, rp), 0.0_rp, rp)* &
               cmplx(0.07_rp*real(j - 1, rp)*space%metric_weights(j), 0.0_rp, rp)
         end do
      end do
      metric_loss = tddft_loss_matrix(space, canonical)
      call check('canonical loss is metric-Hermitian', loss_matrix_hermiticity_residual(space, metric_loss) < tol, &
         test_failed)
      call check('canonical loss has zero origin extension', maxval(abs(metric_loss(1, :))) < tol .and. &
         maxval(abs(metric_loss(:, 1))) < tol, test_failed)
      deallocate(canonical, metric_loss)
   end subroutine test_loss_convention

   subroutine test_scalar_dyson(test_failed)
      logical, intent(inout) :: test_failed
      complex(rp) :: chi_ks(1, 1), kernel(1, 1), chi(1, 1), denominator(1, 1), expected
      real(rp) :: condition_number, minimum_singular, maximum_singular
      integer :: info

      chi_ks(1, 1) = cmplx(2.0_rp, -0.5_rp, rp)
      kernel(1, 1) = cmplx(0.25_rp, 0.1_rp, rp)
      expected = chi_ks(1, 1)/(1.0_rp - chi_ks(1, 1)*kernel(1, 1))
      call solve_tddft_dyson_frequency(chi_ks, kernel, chi, denominator, info, condition_number, &
         minimum_singular, maximum_singular)
      call check('scalar Dyson solve succeeds', info == 0, test_failed)
      call check_complex('scalar Dyson closed form', chi(1, 1), expected, tol, test_failed)
      call check_complex('scalar Dyson denominator', denominator(1, 1), &
         1.0_rp - chi_ks(1, 1)*kernel(1, 1), tol, test_failed)
      call check('scalar condition number is reported', abs(condition_number - 1.0_rp) < tol .and. &
         minimum_singular > 0.0_rp .and. maximum_singular > 0.0_rp, test_failed)
   end subroutine test_scalar_dyson

   subroutine test_small_noncommuting_dyson(test_failed)
      logical, intent(inout) :: test_failed
      complex(rp) :: chi_ks(2, 2), kernel(2, 2), chi(2, 2), denominator(2, 2), expected(2, 2)
      complex(rp) :: d11, d12, d21, d22, determinant
      integer :: info

      chi_ks = reshape([cmplx(0.8_rp, 0.1_rp, rp), cmplx(0.2_rp, -0.3_rp, rp), &
         cmplx(-0.1_rp, 0.25_rp, rp), cmplx(0.45_rp, -0.05_rp, rp)], [2, 2])
      kernel = reshape([cmplx(0.3_rp, -0.1_rp, rp), cmplx(-0.2_rp, 0.15_rp, rp), &
         cmplx(0.4_rp, 0.05_rp, rp), cmplx(-0.1_rp, 0.2_rp, rp)], [2, 2])
      call solve_tddft_dyson_frequency(chi_ks, kernel, chi, denominator, info)
      call check('small noncommuting solve succeeds', info == 0, test_failed)
      call check('small fixture is noncommuting', maxval(abs(matmul(chi_ks, kernel) - matmul(kernel, chi_ks))) > 0.1_rp, &
         test_failed)

      d11 = denominator(1, 1); d12 = denominator(1, 2)
      d21 = denominator(2, 1); d22 = denominator(2, 2)
      determinant = d11*d22 - d12*d21
      expected(1, 1) = (d22*chi_ks(1, 1) - d12*chi_ks(2, 1))/determinant
      expected(2, 1) = (-d21*chi_ks(1, 1) + d11*chi_ks(2, 1))/determinant
      expected(1, 2) = (d22*chi_ks(1, 2) - d12*chi_ks(2, 2))/determinant
      expected(2, 2) = (-d21*chi_ks(1, 2) + d11*chi_ks(2, 2))/determinant
      call check_matrix('small matrix independent solve', chi, expected, tol, test_failed)
   end subroutine test_small_noncommuting_dyson

   subroutine test_metric_sensitive_dyson(test_failed)
      logical, intent(inout) :: test_failed
      type(tddft_dyson_request) :: request
      type(tddft_dyson_result) :: result
      complex(rp), allocatable :: raw_chi(:, :), raw_kernel(:, :), canonical_chi(:, :, :), canonical_kernel(:, :)
      complex(rp), allocatable :: expected_denominator(:, :), expected_chi(:, :), wrong_denominator(:, :), wrong_chi(:, :)
      complex(rp) :: unit_matrix(space%ndim, space%ndim)
      real(rp) :: frequencies(1), wrong_error
      integer :: i, j, k

      allocate(raw_chi(space%ndim, space%ndim), raw_kernel(space%ndim, space%ndim), &
         canonical_chi(space%ndim, space%ndim, 1), canonical_kernel(space%ndim, space%ndim), &
         expected_denominator(space%ndim, space%ndim), expected_chi(space%ndim, space%ndim), &
         wrong_denominator(space%ndim, space%ndim), wrong_chi(space%ndim, space%ndim))
      raw_chi = cmplx(0.0_rp, 0.0_rp, rp)
      raw_kernel = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 2, space%ndim
         do j = 2, space%ndim
            raw_chi(i, j) = cmplx(0.012_rp*real(i + 2*j, rp), 0.002_rp*real(j - i, rp), rp)
            raw_kernel(i, j) = cmplx(0.008_rp*real(2*i + j, rp), -0.001_rp*real(i + j, rp), rp)
         end do
      end do
      do j = 1, space%ndim
         canonical_chi(:, j, 1) = raw_chi(:, j)*space%metric_weights(j)
         canonical_kernel(:, j) = raw_kernel(:, j)*space%metric_weights(j)
      end do
      frequencies = [0.17_rp]
      request%response_space => space
      request%q = [0.13_rp, -0.07_rp, 0.02_rp]
      request%frequencies = frequencies
      request%eta = 0.004_rp
      request%channel = 'chi_plus'
      request%ks_susceptibility = canonical_chi
      request%canonical_interaction = canonical_kernel
      request%interaction_route = lr_dyson_route_direct_alsda
      request%interaction_provenance = 'KXC-01 direct ALSDA LR-03 synthetic fixture'
      request%electronic_state_provenance = 'LR-06 synthetic ham_only complete eigenpair snapshot'
      call evaluate_tddft_dyson(request, result)

      unit_matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, space%ndim
         unit_matrix(i, i) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
      expected_denominator = unit_matrix
      do i = 1, space%ndim
         do j = 1, space%ndim
            do k = 1, space%ndim
               expected_denominator(i, j) = expected_denominator(i, j) - raw_chi(i, k)* &
                  space%metric_weights(k)*raw_kernel(k, j)*space%metric_weights(j)
            end do
         end do
      end do
      expected_chi = canonical_chi(:, :, 1)
      call independent_solve(expected_denominator, expected_chi)
      call check_matrix('metric-sensitive Dyson matches explicit quadrature', result%enhanced_susceptibility(:, :, 1), &
         expected_chi, 10.0_rp*tol, test_failed)

      wrong_denominator = unit_matrix - matmul(raw_chi, raw_kernel)
      wrong_chi = raw_chi
      call independent_solve(wrong_denominator, wrong_chi)
      wrong_error = maxval(abs(result%enhanced_susceptibility(:, :, 1) - wrong_chi))
      call check('weight-free Dyson fixture differs from production result', wrong_error > 1.0e-5_rp, test_failed)
      call check('metric fixture carries LR-04 metadata', index(result%response_representation, 'LR-04') > 0 .and. &
         index(result%response_space_metadata, 'active_dimension=') > 0, test_failed)
      deallocate(raw_chi, raw_kernel, canonical_chi, canonical_kernel, expected_denominator, expected_chi, &
         wrong_denominator, wrong_chi)
   end subroutine test_metric_sensitive_dyson

   subroutine test_goldstone_near_pole(test_failed)
      logical, intent(inout) :: test_failed
      type(tddft_dyson_request) :: request
      type(tddft_dyson_result) :: result, singular_result
      complex(rp), allocatable :: chi(:, :, :), kernel(:, :)
      real(rp) :: eta
      integer :: i

      eta = 1.0e-10_rp
      allocate(chi(space%ndim, space%ndim, 1), kernel(space%ndim, space%ndim))
      chi = cmplx(0.0_rp, 0.0_rp, rp)
      kernel = cmplx(0.0_rp, 0.0_rp, rp)
      chi(2, 2, 1) = cmplx(1.0_rp, eta, rp)
      do i = 3, space%ndim
         chi(i, i, 1) = cmplx(0.5_rp/real(i - 1, rp), 0.0_rp, rp)
      end do
      do i = 2, space%ndim
         kernel(i, i) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
      request%response_space => space
      request%frequencies = [0.0_rp]
      request%eta = eta
      request%channel = 'chi_plus'
      request%ks_susceptibility = chi
      request%canonical_interaction = kernel
      request%interaction_route = lr_dyson_route_goldstone_sumrule
      request%interaction_provenance = 'GSR-01 exact static Goldstone fixture'
      request%electronic_state_provenance = 'LR-06 q=0 static complete eigenpair fixture'
      call evaluate_tddft_dyson(request, result)
      call check('Goldstone near-pole solve succeeds', result%solve_succeeded(1), test_failed)
      call check('Goldstone near-pole is reported, not regularized', result%near_singular_collective_pole(1) .and. &
         .not. result%numerically_singular(1) .and. result%denominator_condition_number(1) > 1.0e8_rp, test_failed)
      call check('Goldstone q and omega provenance retained', maxval(abs(result%q)) == 0.0_rp .and. &
         result%frequencies(1) == 0.0_rp .and. result%eta == eta, test_failed)

      ! Remove the controlled finite eta from the same static identity.  This
      ! is a numerical singularity and must be blocked rather than regularized
      ! or relabelled as a physical collective pole.
      chi(2, 2, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      request%ks_susceptibility = chi
      call evaluate_tddft_dyson(request, singular_result)
      call check('exact Goldstone denominator is numerically singular', singular_result%numerically_singular(1) .and. &
         .not. singular_result%solve_succeeded(1) .and. index(singular_result%status, 'BLOCKED') == 1, test_failed)
      deallocate(chi, kernel)
   end subroutine test_goldstone_near_pole

   subroutine test_covariance(test_failed)
      logical, intent(inout) :: test_failed
      type(tddft_dyson_request) :: plus_request, minus_request, reversed_request
      type(tddft_dyson_result) :: plus_result, minus_result, reversed_result
      complex(rp), allocatable :: chi(:, :, :), kernel(:, :)
      integer :: i, j

      allocate(chi(space%ndim, space%ndim, 2), kernel(space%ndim, space%ndim))
      chi = cmplx(0.0_rp, 0.0_rp, rp)
      kernel = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 2, space%ndim
         kernel(i, i) = cmplx(0.03_rp*real(i, rp), -0.01_rp*real(i - 1, rp), rp)
         do j = 2, space%ndim
            chi(i, j, 1) = cmplx(0.02_rp*real(i + j, rp), 0.003_rp*real(i - j, rp), rp)* &
               space%metric_weights(j)
            chi(i, j, 2) = conjg(chi(i, j, 1))
         end do
      end do

      plus_request%response_space => space
      plus_request%q = [0.11_rp, -0.03_rp, 0.0_rp]
      plus_request%frequencies = [0.08_rp, 0.16_rp]
      plus_request%eta = 0.003_rp
      plus_request%channel = 'chi_plus'
      plus_request%ks_susceptibility = chi
      plus_request%canonical_interaction = kernel
      plus_request%interaction_route = lr_dyson_route_direct_alsda
      plus_request%interaction_provenance = 'KXC-01 covariance fixture'
      plus_request%electronic_state_provenance = 'LR-06 covariance fixture'
      call evaluate_tddft_dyson(plus_request, plus_result)

      minus_request = plus_request
      minus_request%q = -plus_request%q
      minus_request%frequencies = -plus_request%frequencies
      minus_request%channel = 'chi_minus'
      minus_request%ks_susceptibility = conjg(chi)
      minus_request%canonical_interaction = conjg(kernel)
      call evaluate_tddft_dyson(minus_request, minus_result)
      call check_matrix('q/frequency/circular retarded covariance', plus_result%enhanced_susceptibility(:, :, 1), &
         conjg(minus_result%enhanced_susceptibility(:, :, 1)), 50.0_rp*tol, test_failed)
      call check_matrix('q/frequency/circular loss covariance', plus_result%loss_matrix(:, :, 2), &
         conjg(minus_result%loss_matrix(:, :, 2)), 50.0_rp*tol, test_failed)

      reversed_request = plus_request
      reversed_request%channel = 'chi_minus'
      reversed_request%interaction_provenance = 'KXC-01 spin-reversal fixture'
      call evaluate_tddft_dyson(reversed_request, reversed_result)
      call check_matrix('global spin-reversal circular covariance', plus_result%enhanced_susceptibility(:, :, 1), &
         reversed_result%enhanced_susceptibility(:, :, 1), 50.0_rp*tol, test_failed)
      deallocate(chi, kernel)
   end subroutine test_covariance

   subroutine test_interaction_separation(test_failed)
      logical, intent(inout) :: test_failed
      type(tddft_dyson_request) :: request
      type(tddft_dyson_result) :: direct_result, sumrule_result, corrected_result
      complex(rp), allocatable :: chi(:, :, :), direct_kernel(:, :), sumrule_kernel(:, :), corrected_kernel(:, :)

      allocate(chi(space%ndim, space%ndim, 1), direct_kernel(space%ndim, space%ndim), &
         sumrule_kernel(space%ndim, space%ndim), corrected_kernel(space%ndim, space%ndim))
      chi = cmplx(0.0_rp, 0.0_rp, rp)
      direct_kernel = cmplx(0.0_rp, 0.0_rp, rp)
      sumrule_kernel = cmplx(0.0_rp, 0.0_rp, rp)
      corrected_kernel = cmplx(0.0_rp, 0.0_rp, rp)
      chi(2, 2, 1) = cmplx(0.6_rp, -0.02_rp, rp)
      direct_kernel(2, 2) = cmplx(0.2_rp, 0.0_rp, rp)
      sumrule_kernel(2, 2) = cmplx(0.3_rp, 0.0_rp, rp)
      corrected_kernel(2, 2) = cmplx(0.4_rp, 0.0_rp, rp)

      request%response_space => space
      request%frequencies = [0.05_rp]
      request%eta = 0.002_rp
      request%channel = 'chi_plus'
      request%ks_susceptibility = chi
      request%electronic_state_provenance = 'LR-06 route-separation fixture'

      request%canonical_interaction = direct_kernel
      request%interaction_route = lr_dyson_route_direct_alsda
      request%interaction_provenance = 'KXC-01 direct ALSDA fixture'
      call evaluate_tddft_dyson(request, direct_result)
      request%canonical_interaction = sumrule_kernel
      request%interaction_route = lr_dyson_route_goldstone_sumrule
      request%interaction_provenance = 'GSR-01 sum-rule fixture'
      call evaluate_tddft_dyson(request, sumrule_result)
      request%canonical_interaction = corrected_kernel
      request%interaction_route = lr_dyson_route_direct_alsda_goldstone_corrected
      request%interaction_provenance = 'GCR-01 corrected ALSDA fixture'
      call evaluate_tddft_dyson(request, corrected_result)

      call check('direct interaction route is retained', direct_result%interaction_route == lr_dyson_route_direct_alsda, &
         test_failed)
      call check('sum-rule interaction route is retained', sumrule_result%interaction_route == &
         lr_dyson_route_goldstone_sumrule, test_failed)
      call check('corrected interaction route is retained', corrected_result%interaction_route == &
         lr_dyson_route_direct_alsda_goldstone_corrected, test_failed)
      call check('different interaction objects remain separate', &
         maxval(abs(direct_result%canonical_interaction-sumrule_result%canonical_interaction)) > 0.0_rp .and. &
         maxval(abs(sumrule_result%canonical_interaction-corrected_result%canonical_interaction)) > 0.0_rp, test_failed)
      deallocate(chi, direct_kernel, sumrule_kernel, corrected_kernel)
   end subroutine test_interaction_separation

   subroutine independent_solve(matrix, rhs)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp), intent(inout) :: rhs(:, :)
      complex(rp), allocatable :: a(:, :), row(:)
      complex(rp) :: factor, value
      integer :: n, nrhs, i, j, k, pivot

      n = size(matrix, 1)
      nrhs = size(rhs, 2)
      allocate(a(n, n), row(n))
      a = matrix
      do k = 1, n - 1
         pivot = k
         do i = k + 1, n
            if (abs(a(i, k)) > abs(a(pivot, k))) pivot = i
         end do
         if (pivot /= k) then
            row = a(k, :); a(k, :) = a(pivot, :); a(pivot, :) = row
            value = rhs(k, 1); rhs(k, 1) = rhs(pivot, 1); rhs(pivot, 1) = value
            if (nrhs > 1) then
               do j = 2, nrhs
                  value = rhs(k, j); rhs(k, j) = rhs(pivot, j); rhs(pivot, j) = value
               end do
            end if
         end if
         do i = k + 1, n
            factor = a(i, k)/a(k, k)
            a(i, k:n) = a(i, k:n) - factor*a(k, k:n)
            rhs(i, :) = rhs(i, :) - factor*rhs(k, :)
         end do
      end do
      do k = n, 1, -1
         rhs(k, :) = rhs(k, :)/a(k, k)
         do i = 1, k - 1
            rhs(i, :) = rhs(i, :) - a(i, k)*rhs(k, :)
         end do
      end do
      deallocate(a, row)
   end subroutine independent_solve

   subroutine check(label, condition, test_failed)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      logical, intent(inout) :: test_failed

      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', trim(label)
         test_failed = .true.
      end if
   end subroutine check

   subroutine check_complex(label, actual, expected, tolerance, test_failed)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual, expected
      real(rp), intent(in) :: tolerance
      logical, intent(inout) :: test_failed

      call check(label, abs(actual - expected) <= tolerance, test_failed)
   end subroutine check_complex

   subroutine check_matrix(label, actual, expected, tolerance, test_failed)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual(:, :), expected(:, :)
      real(rp), intent(in) :: tolerance
      logical, intent(inout) :: test_failed

      call check(label, maxval(abs(actual - expected)) <= tolerance, test_failed)
   end subroutine check_matrix

end program test_tddft_dyson
