!------------------------------------------------------------------------------
! DRESP-04 algebraic closure:
!   one-site and two-site self-consistent Stoner fixtures,
!   projected Mills recovery, raw Goldstone identity, independent Dyson
!   inverse, loss convention, and q/channel covariance.
!------------------------------------------------------------------------------
program test_lr_projected_interacting_response

   use precision_mod, only: rp
   use lr_projected_interacting_response_mod, only: projected_mills_interaction_request, &
      projected_mills_interaction_result, projected_dyson_request, projected_dyson_result, &
      evaluate_projected_mills_interaction, evaluate_projected_dyson, projected_mills_exact_scalar
   implicit none

   logical :: failed

   failed = .false.
   call one_site_fixture(failed)
   call two_site_fixture(failed)
   if (failed) then
      write (*, '(a)') 'UnitLrProjectedInteractingResponse: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrProjectedInteractingResponse: PASS (Mills, self-consistency, Goldstone, Dyson, loss, covariance)'

contains

   subroutine one_site_fixture(failed)
      logical, intent(inout) :: failed
      type(projected_mills_interaction_request) :: request
      type(projected_mills_interaction_result) :: result
      real(rp), parameter :: u_true = -0.6_rp, moment = 0.7_rp

      allocate(request%actual_pauli_field(1, 1, 1), request%site_vertices(1, 1, 1, 1), request%projected_moment(1))
      request%selector = 'd'
      request%nsite = 1
      request%site_block_size = 1
      request%actual_pauli_field = cmplx(u_true*moment, 0.0_rp, rp)
      request%site_vertices = cmplx(1.0_rp, 0.0_rp, rp)
      request%projected_moment = moment
      request%provenance = 'DRESP-04 one-site exact Stoner fixture'
      call evaluate_projected_mills_interaction(request, result)
      call check('one-site selector plumbing', trim(result%selector) == 'd', failed)
      call check('one-site exact classification', trim(result%classification) == projected_mills_exact_scalar, failed)
      call check('one-site U recovery', abs(result%interaction_U(1) - u_true) < 1.0e-12_rp, failed)
      call check('one-site splitting recovery', abs(result%projected_splitting(1) - u_true*moment) < 1.0e-12_rp, failed)
      call check('one-site scalar residual', result%scalarization_residual < 1.0e-13_rp, failed)
      write (*, '(a,3(es14.6,1x))') 'one_site U Delta residual = ', result%interaction_U(1), &
         result%projected_splitting(1), result%scalarization_residual
   end subroutine one_site_fixture

   subroutine two_site_fixture(failed)
      logical, intent(inout) :: failed
      type(projected_mills_interaction_request) :: mills_request
      type(projected_mills_interaction_result) :: mills_result, approximate_result
      type(projected_dyson_request) :: dyson_request, covariance_request
      type(projected_dyson_result) :: dyson_result, covariance_result
      real(rp), parameter :: u_true = -0.8_rp
      real(rp), parameter :: moment(2) = [0.5_rp, 0.5_rp]
      real(rp), parameter :: b_field = u_true*moment(1)
      real(rp) :: energies(4), occupations(4), q(3), chi_static(2, 2)
      complex(rp) :: eigenvectors(4, 4), plus_operator(4, 4, 2)
      complex(rp) :: chi_dynamic(2, 2), denominator(2, 2), inverse(2, 2), oracle(2, 2)
      complex(rp) :: t(2), residual_vector(2)
      real(rp) :: u_matrix(2, 2), goldstone_residual, inverse_error, covariance_error, loss_error
      integer :: n, m, i, j, info
      integer :: pivots(2), lwork
      complex(rp), allocatable :: work(:)
      complex(rp) :: work_query(1)
      logical :: ok
      external :: zgetrf, zgetri

      ! H0 has a genuine two-site hopping.  The prescribed mean field is
      ! H_MF=H0+B_sigma*sigma_z, B_sigma=U_true*M.  One occupied bonding
      ! up state gives M_i=1/2 on both sites.
      energies = [-0.75_rp, -0.45_rp, 0.05_rp, 0.35_rp]
      occupations = [1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]
      eigenvectors = cmplx(0.0_rp, 0.0_rp, rp)
      eigenvectors(1, 1) = cmplx(1.0_rp/sqrt(2.0_rp), 0.0_rp, rp)
      eigenvectors(2, 1) = eigenvectors(1, 1)
      eigenvectors(1, 2) = eigenvectors(1, 1)
      eigenvectors(2, 2) = -eigenvectors(1, 1)
      eigenvectors(3, 3) = eigenvectors(1, 1)
      eigenvectors(4, 3) = eigenvectors(1, 1)
      eigenvectors(3, 4) = eigenvectors(1, 1)
      eigenvectors(4, 4) = -eigenvectors(1, 1)

      plus_operator = cmplx(0.0_rp, 0.0_rp, rp)
      plus_operator(1, 3, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      plus_operator(2, 4, 2) = cmplx(1.0_rp, 0.0_rp, rp)

      ! Independent Lehmann construction of the DRESP-02 circular chi0.
      chi_static = 0.0_rp
      do n = 1, 4
         do m = 1, 4
            if (abs(occupations(n) - occupations(m)) <= tiny(1.0_rp)) cycle
            do i = 1, 2
               t(i) = dot_product(eigenvectors(:, n), matmul(plus_operator(:, :, i), eigenvectors(:, m)))
            end do
            do i = 1, 2
               do j = 1, 2
                  chi_static(i, j) = chi_static(i, j) + 2.0_rp*(occupations(n) - occupations(m))/ &
                     (energies(n) - energies(m))*real(t(i)*conjg(t(j)), rp)
               end do
            end do
         end do
      end do

      allocate(mills_request%actual_pauli_field(2, 2, 1), mills_request%site_vertices(2, 2, 2, 1), &
         mills_request%projected_moment(2))
      mills_request%selector = 'spd'
      mills_request%nsite = 2
      mills_request%site_block_size = 1
      mills_request%actual_pauli_field = cmplx(0.0_rp, 0.0_rp, rp)
      mills_request%actual_pauli_field(1, 1, 1) = cmplx(b_field, 0.0_rp, rp)
      mills_request%actual_pauli_field(2, 2, 1) = cmplx(b_field, 0.0_rp, rp)
      mills_request%site_vertices = cmplx(0.0_rp, 0.0_rp, rp)
      mills_request%site_vertices(1, 1, 1, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      mills_request%site_vertices(2, 2, 2, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      mills_request%projected_moment = moment
      mills_request%provenance = 'DRESP-04 two-site self-consistent Stoner fixture'
      call evaluate_projected_mills_interaction(mills_request, mills_result)
      call check('two-site selector plumbing', trim(mills_result%selector) == 'spd', failed)
      call check('two-site exact classification', trim(mills_result%classification) == projected_mills_exact_scalar, failed)
      call check('two-site U recovery', maxval(abs(mills_result%interaction_U - u_true)) < 1.0e-12_rp, failed)
      call check('two-site Delta recovery', maxval(abs(mills_result%projected_splitting - b_field)) < 1.0e-12_rp, failed)
      call check('two-site fit residual', mills_result%scalarization_residual < 1.0e-13_rp, failed)

      ! An off-site field is deliberately outside the local scalar Mills ansatz.
      ! The fit must expose that mismatch instead of silently repairing it.
      mills_request%actual_pauli_field(1, 2, 1) = cmplx(0.1_rp, 0.0_rp, rp)
      call evaluate_projected_mills_interaction(mills_request, approximate_result)
      call check('off-site field classified as projected approximation', &
         trim(approximate_result%classification) == 'PROJECTED_SCALAR_APPROXIMATION', failed)
      call check('off-site locality residual is exposed', approximate_result%locality_residual > 0.0_rp, failed)
      mills_request%actual_pauli_field(1, 2, 1) = cmplx(0.0_rp, 0.0_rp, rp)

      u_matrix = 0.0_rp
      u_matrix(1, 1) = u_true
      u_matrix(2, 2) = u_true
      residual_vector = cmplx(0.0_rp, 0.0_rp, rp)
      residual_vector = matmul(cmplx(identity_real(2), 0.0_rp, rp) - matmul(cmplx(chi_static, 0.0_rp, rp), &
         cmplx(u_matrix, 0.0_rp, rp)), cmplx(moment, 0.0_rp, rp))
      goldstone_residual = maxval(abs(residual_vector))
      call check('self-consistent uniform Goldstone mode', goldstone_residual < 1.0e-12_rp, failed)

      chi_dynamic = cmplx(chi_static, 0.0_rp, rp)
      do n = 1, 4
         do m = 1, 4
            if (abs(occupations(n) - occupations(m)) <= tiny(1.0_rp)) cycle
            do i = 1, 2
               t(i) = dot_product(eigenvectors(:, n), matmul(plus_operator(:, :, i), eigenvectors(:, m)))
            end do
            do i = 1, 2
               do j = 1, 2
                  chi_dynamic(i, j) = chi_dynamic(i, j) + 2.0_rp*(occupations(n) - occupations(m))*t(i)*conjg(t(j))* &
                     (1.0_rp/(cmplx(0.07_rp + energies(n) - energies(m), 0.02_rp, rp)) - &
                      1.0_rp/cmplx(energies(n) - energies(m), 0.0_rp, rp))
               end do
            end do
         end do
      end do

      allocate(dyson_request%frequencies(1), dyson_request%interaction_U(2), dyson_request%bare_chi(2, 2, 1))
      dyson_request%selector = 'spd'
      dyson_request%q = [0.137_rp, 0.0_rp, 0.0_rp]
      dyson_request%frequencies = [0.07_rp]
      dyson_request%eta = 0.02_rp
      dyson_request%channel = 'chi_plus'
      dyson_request%interaction_U = mills_result%interaction_U
      dyson_request%bare_chi(:, :, 1) = chi_dynamic
      dyson_request%interaction_provenance = 'Mills fixture U'
      dyson_request%bare_provenance = 'independent finite Hamiltonian Lehmann chi0'
      call evaluate_projected_dyson(dyson_request, dyson_result)
      call check('site Dyson solve succeeds away from pole', dyson_result%solve_info(1) == 0, failed)
      call check('Dyson residual', dyson_result%dyson_residual_relative(1) < 1.0e-12_rp, failed)
      call check('enhancement differs from bare', maxval(abs(dyson_result%enhanced_chi(:, :, 1) - chi_dynamic)) > 1.0e-5_rp, failed)

      denominator = cmplx(identity_real(2), 0.0_rp, rp) - matmul(chi_dynamic, cmplx(u_matrix, 0.0_rp, rp))
      inverse = denominator
      call zgetrf(2, 2, inverse, 2, pivots, info)
      call check('independent zgetrf succeeds', info == 0, failed)
      call zgetri(2, inverse, 2, pivots, work_query, -1, info)
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      call zgetri(2, inverse, 2, pivots, work, lwork, info)
      call check('independent zgetri succeeds', info == 0, failed)
      oracle = matmul(inverse, chi_dynamic)
      inverse_error = maxval(abs(oracle - dyson_result%enhanced_chi(:, :, 1)))
      call check('Dyson versus direct inverse oracle', inverse_error < 1.0e-12_rp, failed)

      loss_error = maxval(abs(dyson_result%loss_matrix(:, :, 1) - &
         conjg(transpose(dyson_result%loss_matrix(:, :, 1)))))
      call check('loss matrix Hermiticity convention', loss_error < 1.0e-12_rp, failed)
      call check('minus Im trace over pi naming/formula', &
         abs(dyson_result%minus_im_trace_over_pi(1) + aimag(dyson_result%enhanced_chi(1,1,1) + &
            dyson_result%enhanced_chi(2,2,1))/acos(-1.0_rp)) < 1.0e-12_rp, failed)

      allocate(covariance_request%frequencies(1), covariance_request%interaction_U(2), covariance_request%bare_chi(2, 2, 1))
      covariance_request%selector = 'spd'
      covariance_request%q = -dyson_request%q
      covariance_request%frequencies = -dyson_request%frequencies
      covariance_request%eta = dyson_request%eta
      covariance_request%channel = 'chi_minus'
      covariance_request%interaction_U = dyson_request%interaction_U
      covariance_request%bare_chi(:, :, 1) = conjg(dyson_request%bare_chi(:, :, 1))
      covariance_request%interaction_provenance = dyson_request%interaction_provenance
      covariance_request%bare_provenance = 'q covariance partner'
      call evaluate_projected_dyson(covariance_request, covariance_result)
      covariance_error = maxval(abs(covariance_result%enhanced_chi(:, :, 1) - &
         conjg(dyson_result%enhanced_chi(:, :, 1))))
      call check('non-self-inverse q/channel covariance', covariance_error < 1.0e-12_rp, failed)

      write (*, '(a,4(es14.6,1x))') 'two_site U Delta scalar_residual Goldstone = ', &
         mills_result%interaction_U(1), mills_result%projected_splitting(1), mills_result%scalarization_residual, &
         goldstone_residual
      write (*, '(a,4(es14.6,1x))') 'Dyson residual inverse loss covariance = ', &
         dyson_result%dyson_residual_relative(1), inverse_error, loss_error, covariance_error
      deallocate(work)
   end subroutine two_site_fixture

   function identity_real(n) result(matrix)
      integer, intent(in) :: n
      real(rp) :: matrix(n, n)
      integer :: i
      matrix = 0.0_rp
      do i = 1, n
         matrix(i, i) = 1.0_rp
      end do
   end function identity_real

   subroutine check(label, condition, failed)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      logical, intent(inout) :: failed
      if (.not. condition) then
         write (*, '(a)') 'FAIL: '//trim(label)
         failed = .true.
      end if
   end subroutine check

end program test_lr_projected_interacting_response
