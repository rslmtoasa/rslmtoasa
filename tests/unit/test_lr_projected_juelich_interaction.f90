!------------------------------------------------------------------------------
! DRESP-05 projected Juelich/LCMM algebra and negative controls.
!------------------------------------------------------------------------------
program test_lr_projected_juelich_interaction

   use precision_mod, only: rp
   use lr_projected_juelich_interaction_mod, only: projected_juelich_request, projected_juelich_result, &
      evaluate_projected_juelich_interaction, evaluate_projected_juelich_holdout, &
      projected_juelich_exact_local, projected_juelich_projected_local, projected_juelich_rank_deficient
   use lr_projected_interacting_response_mod, only: projected_mills_interaction_request, &
      projected_mills_interaction_result, projected_dyson_request, projected_dyson_result, &
      evaluate_projected_mills_interaction, evaluate_projected_dyson
   implicit none

   logical :: failed

   failed = .false.
   call one_site_oracle(failed)
   call two_site_oracles(failed)
   call complex_eta_oracle(failed)
   if (failed) then
      write (*, '(a)') 'UnitLrProjectedJuelichInteraction: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrProjectedJuelichInteraction: PASS (Ward, SVD, real constraint, eta, rank, holdout, LU)'

contains

   subroutine one_site_oracle(failed)
      logical, intent(inout) :: failed
      type(projected_juelich_request) :: request
      type(projected_juelich_result) :: result
      type(projected_mills_interaction_request) :: mills_request
      type(projected_mills_interaction_result) :: mills_result
      type(projected_dyson_request) :: dyson_request, covariance_request
      type(projected_dyson_result) :: dyson_result, covariance_result
      real(rp), parameter :: u_true = -0.8_rp, moment = 0.7_rp

      allocate(request%projected_moment(1), request%static_chi0(1, 1))
      request%selector = 'd'
      request%channel = 'chi_plus'
      request%static_eta = 0.01_rp
      request%projected_moment = moment
      request%static_chi0(1, 1) = cmplx(1.0_rp/u_true, 0.0_rp, rp)
      request%state_provenance = 'one-site self-consistent Stoner fixture'
      request%chi0_provenance = 'independent static finite-H response'
      call evaluate_projected_juelich_interaction(request, result)
      call check('one-site exact classification', trim(result%classification) == projected_juelich_exact_local, failed)
      call check('one-site known U', abs(result%interaction_U_real(1) - u_true) < 1.0e-12_rp, failed)
      call check('one-site complex U', abs(result%interaction_U_complex(1) - cmplx(u_true, 0.0_rp, rp)) < 1.0e-12_rp, failed)
      call check('one-site Ward residual', result%relative_residual < 1.0e-12_rp, failed)
      call check('one-site real residual', result%relative_real_constrained_residual < 1.0e-12_rp, failed)
      allocate(mills_request%actual_pauli_field(1, 1, 1), mills_request%site_vertices(1, 1, 1, 1), &
         mills_request%projected_moment(1))
      mills_request%selector = 'd'
      mills_request%nsite = 1
      mills_request%site_block_size = 1
      mills_request%actual_pauli_field = cmplx(u_true*moment, 0.0_rp, rp)
      mills_request%site_vertices = cmplx(1.0_rp, 0.0_rp, rp)
      mills_request%projected_moment = moment
      call evaluate_projected_mills_interaction(mills_request, mills_result)
      call check('one-site Mills equality', abs(result%interaction_U(1) - mills_result%interaction_U(1)) < 1.0e-12_rp, failed)
      allocate(dyson_request%frequencies(1), dyson_request%interaction_U(1), dyson_request%bare_chi(1, 1, 1))
      dyson_request%selector = 'd'
      dyson_request%q = [0.137_rp, 0.0_rp, 0.0_rp]
      dyson_request%frequencies = [0.04_rp]
      dyson_request%eta = 0.01_rp
      dyson_request%channel = 'chi_plus'
      dyson_request%interaction_U = result%interaction_U_real
      dyson_request%bare_chi(1, 1, 1) = cmplx(0.2_rp, 0.03_rp, rp)
      dyson_request%interaction_provenance = 'one-site Juelich U'
      dyson_request%bare_provenance = 'one-site Dyson reuse fixture'
      call evaluate_projected_dyson(dyson_request, dyson_result)
      call check('Juelich Dyson reuse', dyson_result%solve_info(1) == 0, failed)
      call check('Juelich Dyson residual', dyson_result%dyson_residual_relative(1) < 1.0e-12_rp, failed)
      covariance_request = dyson_request
      covariance_request%q = -dyson_request%q
      covariance_request%frequencies = -dyson_request%frequencies
      covariance_request%channel = 'chi_minus'
      covariance_request%bare_chi = conjg(dyson_request%bare_chi)
      call evaluate_projected_dyson(covariance_request, covariance_result)
      call check('Juelich q/channel covariance', maxval(abs(covariance_result%enhanced_chi - &
         conjg(dyson_result%enhanced_chi))) < 1.0e-12_rp, failed)
      write (*, '(a,3(es14.6,1x))') 'one_site U complex_ratio residual = ', result%interaction_U(1), &
         result%imaginary_U_ratio, result%relative_real_constrained_residual
   end subroutine one_site_oracle

   subroutine two_site_oracles(failed)
      logical, intent(inout) :: failed
      type(projected_juelich_request) :: request
      type(projected_juelich_result) :: result, equal_result, rank_result
      complex(rp) :: gamma(2, 2), equal_gamma(2, 2), inverse(2, 2), oracle(2), work_query(1)
      complex(rp), allocatable :: work(:), holdout_chi(:, :)
      real(rp) :: moment(2), unequal_u(2), equal_u(2), holdout, holdout_relative
      real(rp) :: inverse_error
      integer :: pivots(2), info, lwork
      external :: zgetrf, zgetri

      moment = [0.4_rp, 0.7_rp]
      unequal_u = [-0.6_rp, -0.9_rp]
      equal_u = [-0.75_rp, -0.75_rp]
      gamma = reshape([cmplx(1.3_rp, 0.0_rp, rp), cmplx(0.25_rp, 0.0_rp, rp), &
                       cmplx(-1.3111111111111111_rp, 0.0_rp, rp), cmplx(-0.9444444444444444_rp, 0.0_rp, rp)], [2, 2])

      call initialize_request(request, 'spd', moment, gamma, 0.01_rp)
      call evaluate_projected_juelich_interaction(request, result)
      call check('two-site unequal classification', trim(result%classification) == projected_juelich_exact_local, failed)
      call check('two-site unequal U recovery', maxval(abs(result%interaction_U_real - unequal_u)) < 1.0e-12_rp, failed)
      call check('two-site unequal real residual', result%relative_real_constrained_residual < 1.0e-12_rp, failed)

      ! Independent inverse/LU oracle for Gamma U=M.
      inverse = gamma
      call zgetrf(2, 2, inverse, 2, pivots, info)
      call check('two-site LU factorization', info == 0, failed)
      call zgetri(2, inverse, 2, pivots, work_query, -1, info)
      call check('two-site inverse workspace query', info == 0, failed)
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      call zgetri(2, inverse, 2, pivots, work, lwork, info)
      call check('two-site LU inverse', info == 0, failed)
      oracle = matmul(inverse, cmplx(moment, 0.0_rp, rp))
      inverse_error = maxval(abs(oracle - result%interaction_U_complex))
      call check('two-site Juelich versus LU', inverse_error < 1.0e-12_rp, failed)
      deallocate(work)

      ! Equal-U is a separate case; passing unequal U above prevents a uniform
      ! scalar from satisfying the test by construction.
      equal_gamma = reshape([cmplx(1.3_rp, 0.0_rp, rp), cmplx(0.4_rp, 0.0_rp, rp), &
                             cmplx(-1.8333333333333333_rp, 0.0_rp, rp), cmplx(-1.3333333333333333_rp, 0.0_rp, rp)], [2, 2])
      call initialize_request(request, 'spd', moment, equal_gamma, 0.01_rp)
      call evaluate_projected_juelich_interaction(request, equal_result)
      call check('two-site equal-U recovery', maxval(abs(equal_result%interaction_U_real - equal_u)) < 1.0e-12_rp, failed)

      ! Rank deficiency is visible, never hidden by a cutoff or a Goldstone
      ! adjustment.
      gamma(:, 2) = gamma(:, 1)
      call initialize_request(request, 'd', [1.0_rp, 1.0_rp], gamma, 0.01_rp)
      call evaluate_projected_juelich_interaction(request, rank_result)
      call check('rank-deficient status', trim(rank_result%classification) == projected_juelich_rank_deficient, failed)
      call check('rank-deficient rank', rank_result%rank < 2, failed)

      ! A frozen diagonal fit does not reproduce a deliberately changed
      ! holdout response; this is the negative control for a nonlocal/model
      ! mismatch that the square static solve alone cannot expose.
      allocate(holdout_chi(2, 2))
      holdout_chi = request%static_chi0
      holdout_chi(1, 2) = holdout_chi(1, 2) + cmplx(0.31_rp, 0.0_rp, rp)
      call evaluate_projected_juelich_holdout(holdout_chi, [1.0_rp, 1.0_rp], rank_result%interaction_U_real, holdout, holdout_relative)
      call check('holdout nonlocal/model residual', holdout_relative > 1.0e-3_rp, failed)
      call check('holdout is not construction residual', holdout > 1.0e-3_rp, failed)
      deallocate(holdout_chi)
      write (*, '(a,3(es14.6,1x))') 'two_site inverse_error holdout_relative holdout = ', inverse_error, holdout_relative, holdout
   end subroutine two_site_oracles

   subroutine complex_eta_oracle(failed)
      logical, intent(inout) :: failed
      type(projected_juelich_request) :: request
      type(projected_juelich_result) :: coarse, fine
      real(rp), parameter :: u_true = 0.9_rp, moment = 0.8_rp

      allocate(request%projected_moment(1), request%static_chi0(1, 1))
      request%selector = 'd'
      request%channel = 'chi_plus'
      request%projected_moment = moment
      request%state_provenance = 'finite-broadening Stoner fixture'
      request%chi0_provenance = 'known analytic eta ladder'
      request%static_eta = 0.04_rp
      request%static_chi0(1, 1) = cmplx(1.0_rp/u_true, 0.04_rp, rp)
      call evaluate_projected_juelich_interaction(request, coarse)
      request%static_eta = 0.005_rp
      request%static_chi0(1, 1) = cmplx(1.0_rp/u_true, 0.005_rp, rp)
      call evaluate_projected_juelich_interaction(request, fine)
      call check('finite eta complex U diagnostic', coarse%imaginary_U_ratio > 1.0e-3_rp, failed)
      call check('finite eta real constrained residual reported', coarse%relative_real_constrained_residual > 1.0e-5_rp, failed)
      call check('eta-to-zero real U convergence', abs(fine%interaction_U_real(1) - u_true) < &
         abs(coarse%interaction_U_real(1) - u_true), failed)
      call check('eta-to-zero imaginary convergence', fine%imaginary_U_ratio < coarse%imaginary_U_ratio, failed)
      write (*, '(a,4(es14.6,1x))') 'eta coarse/fine Ureal imag_ratio residual = ', coarse%interaction_U_real(1), &
         fine%interaction_U_real(1), coarse%imaginary_U_ratio, fine%imaginary_U_ratio
   end subroutine complex_eta_oracle

   subroutine initialize_request(request, selector, moment, gamma, eta)
      type(projected_juelich_request), intent(out) :: request
      character(len=*), intent(in) :: selector
      real(rp), intent(in) :: moment(:), eta
      complex(rp), intent(in) :: gamma(:, :)
      integer :: i

      allocate(request%projected_moment(size(moment)), request%static_chi0(size(moment), size(moment)))
      request%selector = selector
      request%channel = 'chi_plus'
      request%static_eta = eta
      request%projected_moment = moment
      do i = 1, size(moment)
         request%static_chi0(:, i) = gamma(:, i)/moment(i)
      end do
      request%state_provenance = 'two-site coupled finite-H fixture'
      request%chi0_provenance = 'independent eigenpair/static response construction'
   end subroutine initialize_request

   subroutine check(label, condition, failed)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      logical, intent(inout) :: failed
      if (.not. condition) then
         write (*, '(a)') 'FAIL: '//trim(label)
         failed = .true.
      end if
   end subroutine check

end program test_lr_projected_juelich_interaction
