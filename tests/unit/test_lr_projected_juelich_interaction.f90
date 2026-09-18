!------------------------------------------------------------------------------
! DRESP-05 projected Juelich/LCMM algebra and negative controls.
!------------------------------------------------------------------------------
program test_lr_projected_juelich_interaction

   use precision_mod, only: rp
   use lr_projected_juelich_interaction_mod, only: projected_juelich_request, projected_juelich_result, &
      evaluate_projected_juelich_interaction, evaluate_projected_juelich_holdout, &
      projected_juelich_exact_local, projected_juelich_projected_local, projected_juelich_rank_deficient, &
      select_projected_juelich_eta_indices
   use lr_projected_interacting_response_mod, only: projected_mills_interaction_request, &
      projected_mills_interaction_result, projected_dyson_request, projected_dyson_result, &
      evaluate_projected_mills_interaction, evaluate_projected_dyson
   implicit none

   logical :: failed

   failed = .false.
   call one_site_oracle(failed)
   call two_site_oracles(failed)
   call complex_eta_oracle(failed)
   call finite_h_known_u_oracle(failed)
   call eta_selection_oracle(failed)
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
      complex(rp) :: complex_holdout(1, 1), real_holdout(1, 1)
      real(rp) :: complex_residual, complex_relative, real_residual, real_relative
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
      complex_holdout(1, 1) = cmplx(1.0_rp/u_true, 0.04_rp, rp)
      real_holdout(1, 1) = cmplx(1.0_rp/u_true, 0.0_rp, rp)
      call evaluate_projected_juelich_holdout(complex_holdout, [moment], fine%interaction_U_real, &
         complex_residual, complex_relative)
      call evaluate_projected_juelich_holdout(real_holdout, [moment], fine%interaction_U_real, &
         real_residual, real_relative)
      call check('complex holdout uses imaginary response', complex_relative > real_relative + 1.0e-3_rp, failed)
      write (*, '(a,4(es14.6,1x))') 'eta coarse/fine Ureal imag_ratio residual = ', coarse%interaction_U_real(1), &
         fine%interaction_U_real(1), coarse%imaginary_U_ratio, fine%imaginary_U_ratio
   end subroutine complex_eta_oracle

   ! Build H_up/H_down with explicit hopping, obtain the moments and static
   ! circular response independently from their eigensystems, and only then
   ! give chi0 and M to the DRESP-05 solver.  No Gamma=M/U construction is
   ! used here.  The two field choices exercise equal and unequal local U.
   subroutine finite_h_known_u_oracle(failed)
      logical, intent(inout) :: failed
      type(projected_juelich_request) :: request
      type(projected_juelich_result) :: result
      real(rp) :: h(2, 2), bfield(2), eval_up(2), eval_down(2), vec_up(2, 2), vec_down(2, 2)
      real(rp) :: occupations_up(2), occupations_down(2), moment(2), known_u(2), fermi
      complex(rp) :: chi0(2, 2), transition(2), denominator
      real(rp), parameter :: eta = 1.0e-9_rp
      integer :: i, j, site

      h = reshape([0.0_rp, 0.23_rp, 0.23_rp, 0.0_rp], [2, 2])
      bfield = [0.31_rp, 0.31_rp]
      call diagonalize_collinear_pair(h, bfield, eval_up, vec_up, eval_down, vec_down)
      call fill_two_electron_occupations(eval_up, eval_down, occupations_up, occupations_down, fermi)
      call moments_from_eigenpairs(vec_up, vec_down, occupations_up, occupations_down, moment)
      known_u = bfield/moment
      call build_circular_chi0(eval_up, eval_down, vec_up, vec_down, occupations_up, occupations_down, eta, chi0)
      call initialize_finite_h_request(request, moment, chi0, eta)
      call evaluate_projected_juelich_interaction(request, result)
      call check('finite-H equal-U hopping is nonzero', abs(h(1, 2)) > 0.0_rp, failed)
      call check('finite-H equal-U moment is nonzero', minval(abs(moment)) > 1.0e-8_rp, failed)
      call check('finite-H equal-U recovery', maxval(abs(result%interaction_U_real - known_u)) < 2.0e-6_rp, failed)
      write (*, '(a,2(es14.6,1x))') 'finite-H equal-U known/recovered = ', known_u(1), result%interaction_U_real(1)

      h = reshape([0.0_rp, 0.19_rp, 0.19_rp, 0.13_rp], [2, 2])
      bfield = [0.23_rp, 0.41_rp]
      call diagonalize_collinear_pair(h, bfield, eval_up, vec_up, eval_down, vec_down)
      call fill_two_electron_occupations(eval_up, eval_down, occupations_up, occupations_down, fermi)
      call moments_from_eigenpairs(vec_up, vec_down, occupations_up, occupations_down, moment)
      known_u = bfield/moment
      call build_circular_chi0(eval_up, eval_down, vec_up, vec_down, occupations_up, occupations_down, eta, chi0)
      call initialize_finite_h_request(request, moment, chi0, eta)
      call evaluate_projected_juelich_interaction(request, result)
      call check('finite-H unequal-U moments are nonzero', minval(abs(moment)) > 1.0e-8_rp, failed)
      call check('finite-H unequal-U input differs', abs(known_u(1) - known_u(2)) > 1.0e-4_rp, failed)
      call check('finite-H unequal-U recovery', maxval(abs(result%interaction_U_real - known_u)) < 2.0e-6_rp, failed)
      write (*, '(a,4(es14.6,1x))') 'finite-H unequal-U known/recovered = ', known_u(1), known_u(2), &
         result%interaction_U_real(1), result%interaction_U_real(2)
   end subroutine finite_h_known_u_oracle

   subroutine initialize_finite_h_request(request, moment, chi0, eta)
      type(projected_juelich_request), intent(out) :: request
      real(rp), intent(in) :: moment(:), eta
      complex(rp), intent(in) :: chi0(:, :)
      allocate(request%projected_moment(size(moment)), request%static_chi0(size(chi0, 1), size(chi0, 2)))
      request%selector = 'spd'
      request%channel = 'chi_plus'
      request%static_eta = eta
      request%projected_moment = moment
      request%static_chi0 = chi0
      request%state_provenance = 'independent coupled finite-H two-site eigensystem'
      request%chi0_provenance = 'direct circular Lehmann sum from H_up/H_down eigenpairs'
   end subroutine initialize_finite_h_request

   subroutine diagonalize_collinear_pair(h, bfield, eval_up, vec_up, eval_down, vec_down)
      real(rp), intent(in) :: h(:, :), bfield(:)
      real(rp), intent(out) :: eval_up(:), eval_down(:), vec_up(:, :), vec_down(:, :)
      real(rp) :: hup(2, 2), hdn(2, 2)
      hup = h
      hdn = h
      hup(1, 1) = hup(1, 1) + bfield(1)
      hup(2, 2) = hup(2, 2) + bfield(2)
      hdn(1, 1) = hdn(1, 1) - bfield(1)
      hdn(2, 2) = hdn(2, 2) - bfield(2)
      call diagonalize_real_2x2(hup, eval_up, vec_up)
      call diagonalize_real_2x2(hdn, eval_down, vec_down)
   end subroutine diagonalize_collinear_pair

   subroutine diagonalize_real_2x2(matrix, values, vectors)
      real(rp), intent(in) :: matrix(2, 2)
      real(rp), intent(out) :: values(2), vectors(2, 2)
      real(rp) :: center, half_width, lambda, norm
      integer :: ib
      center = 0.5_rp*(matrix(1, 1) + matrix(2, 2))
      half_width = sqrt(0.25_rp*(matrix(1, 1) - matrix(2, 2))**2 + matrix(1, 2)**2)
      values = [center - half_width, center + half_width]
      do ib = 1, 2
         lambda = values(ib)
         vectors(:, ib) = [matrix(1, 2), lambda - matrix(1, 1)]
         norm = sqrt(sum(vectors(:, ib)**2))
         vectors(:, ib) = vectors(:, ib)/norm
      end do
   end subroutine diagonalize_real_2x2

   subroutine fill_two_electron_occupations(eval_up, eval_down, occ_up, occ_down, fermi)
      real(rp), intent(in) :: eval_up(:), eval_down(:)
      real(rp), intent(out) :: occ_up(:), occ_down(:), fermi
      real(rp) :: all_values(4), ordered(4)
      integer :: i, j, k
      all_values = [eval_up, eval_down]
      ordered = all_values
      do i = 1, 3
         do j = i + 1, 4
            if (ordered(j) < ordered(i)) then
               fermi = ordered(i); ordered(i) = ordered(j); ordered(j) = fermi
            end if
         end do
      end do
      fermi = 0.5_rp*(ordered(2) + ordered(3))
      occ_up = 0.0_rp
      occ_down = 0.0_rp
      do i = 1, 2
         if (eval_up(i) < fermi) occ_up(i) = 1.0_rp
         if (eval_down(i) < fermi) occ_down(i) = 1.0_rp
      end do
   end subroutine fill_two_electron_occupations

   subroutine moments_from_eigenpairs(vec_up, vec_down, occ_up, occ_down, moment)
      real(rp), intent(in) :: vec_up(:, :), vec_down(:, :), occ_up(:), occ_down(:)
      real(rp), intent(out) :: moment(:)
      integer :: site, ib
      moment = 0.0_rp
      do site = 1, 2
         do ib = 1, 2
            moment(site) = moment(site) + occ_up(ib)*vec_up(site, ib)**2 - occ_down(ib)*vec_down(site, ib)**2
         end do
      end do
   end subroutine moments_from_eigenpairs

   subroutine build_circular_chi0(eval_up, eval_down, vec_up, vec_down, occ_up, occ_down, eta, chi0)
      real(rp), intent(in) :: eval_up(:), eval_down(:), vec_up(:, :), vec_down(:, :), occ_up(:), occ_down(:), eta
      complex(rp), intent(out) :: chi0(:, :)
      complex(rp) :: factor
      real(rp) :: occupation_difference
      complex(rp) :: transition(2)
      integer :: i, j, site, other
      chi0 = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, 2
         do j = 1, 2
            occupation_difference = occ_up(i) - occ_down(j)
            if (occupation_difference == 0.0_rp) cycle
            factor = 2.0_rp*occupation_difference/cmplx(eval_up(i) - eval_down(j), eta, rp)
            do site = 1, 2
               transition(site) = cmplx(vec_up(site, i)*vec_down(site, j), 0.0_rp, rp)
            end do
            do site = 1, 2
               do other = 1, 2
                  chi0(site, other) = chi0(site, other) + factor*transition(site)*conjg(transition(other))
               end do
            end do
         end do
      end do
   end subroutine build_circular_chi0

   subroutine eta_selection_oracle(failed)
      logical, intent(inout) :: failed
      real(rp) :: scrambled(3)
      integer :: selected, holdout
      scrambled = [0.02_rp, 0.005_rp, 0.04_rp]
      call select_projected_juelich_eta_indices(scrambled, selected, holdout)
      call check('scrambled eta selects minimum by value', selected == 2, failed)
      call check('scrambled eta selects second-smallest holdout', holdout == 1, failed)
      write (*, '(a,2(i0,1x))') 'eta selection selected/holdout indices = ', selected, holdout
   end subroutine eta_selection_oracle

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
