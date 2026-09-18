!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief DRESP-05 projected Juelich/LCMM site-space sum-rule interaction.
!>
!> The site Ward identity is derived in the DRESP-01/02 normalization as
!>
!>   Gamma(i,j) = chi0(i,j;0) M(j),   Gamma U = M.
!>
!> This is intentionally independent of the radial GSR solver.  In particular,
!> no 4*pi, magnetic moment, or Goldstone shift is inserted here.
!------------------------------------------------------------------------------
module lr_projected_juelich_interaction_mod

   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   implicit none
   private

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

   public :: evaluate_projected_juelich_interaction
   public :: evaluate_projected_juelich_holdout
   public :: assess_projected_juelich_eta_stability
   public :: select_projected_juelich_eta_indices

   interface
      subroutine zgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info)
         import :: rp
         integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
         complex(rp), intent(inout) :: a(lda, *), b(ldb, *), work(*)
         real(rp), intent(out) :: s(*)
         real(rp), intent(in) :: rcond
         integer, intent(out) :: rank, info
         real(rp), intent(inout) :: rwork(*)
      end subroutine zgelss
      subroutine dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
         import :: rp
         integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
         real(rp), intent(inout) :: a(lda, *), b(ldb, *), work(*)
         real(rp), intent(out) :: s(*)
         real(rp), intent(in) :: rcond
         integer, intent(out) :: rank, info
      end subroutine dgelss
   end interface

contains

   !> Select the finest positive static eta and the distinct second-finest
   !> positive eta independently of the order supplied by the user.  The
   !> second value is the holdout/stability partner; an ambiguous duplicate is
   !> rejected instead of being treated as a stability comparison.
   subroutine select_projected_juelich_eta_indices(eta_values, selected_index, holdout_index)
      real(rp), intent(in) :: eta_values(:)
      integer, intent(out) :: selected_index, holdout_index
      integer :: i, j
      real(rp) :: second_eta, scale

      if (size(eta_values) < 1) error stop 'DRESP-05 eta selection: empty eta ladder'
      if (any(.not. ieee_is_finite(eta_values)) .or. any(eta_values <= 0.0_rp)) then
         error stop 'DRESP-05 eta selection: every eta must be finite and positive'
      end if
      do i = 1, size(eta_values) - 1
         do j = i + 1, size(eta_values)
            scale = max(1.0_rp, abs(eta_values(i)), abs(eta_values(j)))
            if (abs(eta_values(i) - eta_values(j)) <= 100.0_rp*epsilon(1.0_rp)*scale) then
               error stop 'DRESP-05 eta selection: duplicated eta makes stability comparison ambiguous'
            end if
         end do
      end do
      selected_index = 1
      do i = 2, size(eta_values)
         if (eta_values(i) < eta_values(selected_index)) selected_index = i
      end do
      holdout_index = 0
      second_eta = huge(1.0_rp)
      do i = 1, size(eta_values)
         if (i == selected_index) cycle
         if (eta_values(i) < second_eta) then
            second_eta = eta_values(i)
            holdout_index = i
         end if
      end do
   end subroutine select_projected_juelich_eta_indices

   subroutine evaluate_projected_juelich_interaction(request, result)
      type(projected_juelich_request), intent(in) :: request
      type(projected_juelich_result), intent(out) :: result

      complex(rp), allocatable :: gamma(:, :), rhs(:, :), a_work(:, :), b_work(:, :), work(:)
      real(rp), allocatable :: singular_values(:), rwork(:), real_a(:, :), real_b(:, :), real_singular_values(:), &
         real_work(:), real_work_query(:)
      complex(rp) :: complex_work_query(1)
      real(rp) :: complex_scale
      integer :: nsite, lwork, real_lwork, info, i, j
      integer :: rank, real_rank, mreal, ldb
      real(rp), parameter :: rcond = -1.0_rp

      call validate_request(request, nsite)
      result%selector = trim(request%selector)
      result%static_eta = request%static_eta
      result%q = request%q
      result%channel = trim(request%channel)
      result%provenance = trim(request%state_provenance)//' | chi0='//trim(request%chi0_provenance)
      result%nsite = nsite
      allocate(result%projected_moment(nsite), result%gamma(nsite, nsite), result%interaction_matrix(nsite, nsite), &
         result%interaction_U_complex(nsite), result%interaction_U_real(nsite), result%interaction_U(nsite), &
         result%singular_values(nsite), result%real_singular_values(nsite))
      result%projected_moment = request%projected_moment
      result%gamma = cmplx(0.0_rp, 0.0_rp, rp)
      do j = 1, nsite
         result%gamma(:, j) = request%static_chi0(:, j)*request%projected_moment(j)
      end do
      result%interaction_matrix = cmplx(0.0_rp, 0.0_rp, rp)
      result%interaction_U_complex = cmplx(0.0_rp, 0.0_rp, rp)
      result%interaction_U_real = 0.0_rp
      result%interaction_U = 0.0_rp

      allocate(a_work(nsite, nsite), b_work(nsite, 1), rwork(max(1, 5*nsite)))
      a_work = result%gamma
      b_work(:, 1) = cmplx(request%projected_moment, 0.0_rp, rp)
      call zgelss(nsite, nsite, 1, a_work, nsite, b_work, nsite, result%singular_values, rcond, rank, &
         complex_work_query, -1, rwork, info)
      if (info /= 0) error stop 'DRESP-05 Juelich: complex ZGELSS workspace query failed'
      lwork = max(1, nint(real(complex_work_query(1), rp)))
      allocate(work(lwork))
      a_work = result%gamma
      b_work(:, 1) = cmplx(request%projected_moment, 0.0_rp, rp)
      call zgelss(nsite, nsite, 1, a_work, nsite, b_work, nsite, result%singular_values, rcond, rank, &
         work, lwork, rwork, info)
      if (info /= 0) error stop 'DRESP-05 Juelich: complex ZGELSS failed'
      result%rank = rank
      result%interaction_U_complex = b_work(:, 1)
      complex_scale = sqrt(sum(abs(request%projected_moment)**2))
      result%equation_residual = sqrt(sum(abs(matmul(result%gamma, result%interaction_U_complex) - &
         cmplx(request%projected_moment, 0.0_rp, rp))**2))
      result%relative_residual = result%equation_residual/max(complex_scale, tiny(1.0_rp))
      result%condition_number = singular_condition(result%singular_values, result%rank)
      result%imaginary_U_ratio = maxval(abs(aimag(result%interaction_U_complex)))/ &
         max(maxval(abs(real(result%interaction_U_complex, rp))), epsilon(1.0_rp))

      ! Route B: the unknown is real.  Use DGELSS on the stacked real and
      ! imaginary equations, so taking REAL(U_complex) is never substituted for
      ! the physical constrained solve.
      mreal = 2*nsite
      ldb = mreal
      allocate(real_a(mreal, nsite), real_b(ldb, 1), real_singular_values(nsite), real_work_query(1))
      real_a(1:nsite, :) = real(result%gamma, rp)
      real_a(nsite+1:mreal, :) = aimag(result%gamma)
      real_b = 0.0_rp
      real_b(1:nsite, 1) = request%projected_moment
      call dgelss(mreal, nsite, 1, real_a, mreal, real_b, ldb, real_singular_values, rcond, real_rank, &
         real_work_query, -1, info)
      if (info /= 0) error stop 'DRESP-05 Juelich: real DGELSS workspace query failed'
      real_lwork = max(1, nint(real_work_query(1)))
      allocate(real_work(real_lwork))
      real_a(1:nsite, :) = real(result%gamma, rp)
      real_a(nsite+1:mreal, :) = aimag(result%gamma)
      real_b = 0.0_rp
      real_b(1:nsite, 1) = request%projected_moment
      call dgelss(mreal, nsite, 1, real_a, mreal, real_b, ldb, real_singular_values, rcond, real_rank, &
         real_work, real_lwork, info)
      if (info /= 0) error stop 'DRESP-05 Juelich: real DGELSS failed'
      result%real_rank = real_rank
      result%interaction_U_real = real_b(1:nsite, 1)
      result%real_constrained_residual = sqrt(sum((matmul(real(result%gamma, rp), result%interaction_U_real) - &
         request%projected_moment)**2) + sum((matmul(aimag(result%gamma), result%interaction_U_real))**2))
      result%relative_real_constrained_residual = result%real_constrained_residual/max(real_scale_from(request%projected_moment), &
         tiny(1.0_rp))
      result%real_condition_number = singular_condition(real_singular_values, result%real_rank)
      do i = 1, nsite
         result%interaction_matrix(i, i) = cmplx(result%interaction_U_real(i), 0.0_rp, rp)
      end do
      result%interaction_U = result%interaction_U_real

      if (rank < nsite .or. real_rank < nsite) then
         result%classification = projected_juelich_rank_deficient
      else if (result%condition_number > projected_juelich_condition_limit .or. &
               result%real_condition_number > projected_juelich_condition_limit) then
         result%classification = projected_juelich_unsupported
      else if (result%relative_real_constrained_residual > projected_juelich_residual_tolerance) then
         result%classification = projected_juelich_projected_local
      else
         result%classification = projected_juelich_exact_local
      end if

      deallocate(a_work, b_work, work, rwork, real_a, real_b, real_singular_values, real_work_query, real_work)
   end subroutine evaluate_projected_juelich_interaction

   subroutine evaluate_projected_juelich_holdout(static_chi0, projected_moment, interaction_U, residual, relative_residual)
      complex(rp), intent(in) :: static_chi0(:, :)
      real(rp), intent(in) :: projected_moment(:), interaction_U(:)
      real(rp), intent(out) :: residual, relative_residual
      complex(rp), allocatable :: vector(:)
      integer :: nsite

      nsite = size(projected_moment)
      if (any(shape(static_chi0) /= [nsite, nsite]) .or. size(interaction_U) /= nsite) then
         error stop 'DRESP-05 Juelich holdout: inconsistent dimensions'
      end if
      allocate(vector(nsite))
      ! Holdout means that U is frozen while chi0 is replaced.  The residual
      ! is therefore the Ward action (I-chi0*U)M, not the construction solve.
      vector = cmplx(projected_moment, 0.0_rp, rp) - matmul(static_chi0, &
         cmplx(interaction_U*projected_moment, 0.0_rp, rp))
      residual = sqrt(sum(abs(vector)**2))
      relative_residual = residual/max(sqrt(sum(projected_moment**2)), tiny(1.0_rp))
      deallocate(vector)
   end subroutine evaluate_projected_juelich_holdout

   subroutine assess_projected_juelich_eta_stability(result, previous_result, is_stable)
      type(projected_juelich_result), intent(inout) :: result
      type(projected_juelich_result), intent(in) :: previous_result
      logical, intent(out) :: is_stable
      real(rp) :: scale, difference

      if (.not. allocated(result%interaction_U_real) .or. .not. allocated(previous_result%interaction_U_real) .or. &
          size(result%interaction_U_real) /= size(previous_result%interaction_U_real)) then
         error stop 'DRESP-05 Juelich eta audit: incompatible results'
      end if
      difference = maxval(abs(result%interaction_U_real - previous_result%interaction_U_real))
      scale = max(maxval(abs(result%interaction_U_real)), maxval(abs(previous_result%interaction_U_real)), tiny(1.0_rp))
      is_stable = difference/scale <= projected_juelich_eta_relative_tolerance
      result%eta_stable = is_stable
      write(result%eta_stability, '(a,es12.4,a,es12.4)') 'relative_delta_U=', difference/scale, &
         '; threshold=', projected_juelich_eta_relative_tolerance
      if (.not. is_stable .and. trim(result%classification) /= projected_juelich_rank_deficient .and. &
          trim(result%classification) /= projected_juelich_unsupported) then
         result%classification = projected_juelich_eta_limited
      end if
   end subroutine assess_projected_juelich_eta_stability

   subroutine validate_request(request, nsite)
      type(projected_juelich_request), intent(in) :: request
      integer, intent(out) :: nsite

      nsite = size(request%projected_moment)
      if (nsite < 1 .or. any(shape(request%static_chi0) /= [nsite, nsite]) .or. request%static_eta <= 0.0_rp .or. &
          len_trim(request%selector) == 0 .or. len_trim(request%channel) == 0) then
         error stop 'DRESP-05 Juelich: invalid projected request'
      end if
      if (any(.not. ieee_is_finite(request%projected_moment)) .or. &
          any(.not. ieee_is_finite(real(request%static_chi0, rp))) .or. &
          any(.not. ieee_is_finite(aimag(request%static_chi0)))) then
         error stop 'DRESP-05 Juelich: non-finite projected input'
      end if
   end subroutine validate_request

   real(rp) function singular_condition(singular_values, rank) result(value)
      real(rp), intent(in) :: singular_values(:)
      integer, intent(in) :: rank
      real(rp) :: smallest

      if (rank < 1) then
         value = huge(1.0_rp)
         return
      end if
      smallest = minval(singular_values(1:rank))
      if (smallest <= 0.0_rp) then
         value = huge(1.0_rp)
      else
         value = maxval(singular_values)/smallest
      end if
   end function singular_condition

   real(rp) function real_scale_from(vector) result(value)
      real(rp), intent(in) :: vector(:)
      value = sqrt(sum(vector**2))
   end function real_scale_from

end module lr_projected_juelich_interaction_mod
