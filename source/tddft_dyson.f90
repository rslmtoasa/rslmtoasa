!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief LR-04 Dyson enhancement and anti-Hermitian loss matrix.
!>
!> The production representation is the LR-04 canonical right-weighted
!> operator B=A*W.  Consequently the Dyson denominator is formed by ordinary
!> multiplication of canonical operators,
!>
!>   D = I - chi_KS K,
!>
!> and the enhanced response is obtained from D*chi=chi_KS.  No radial metric
!> is inserted by this module: it is already present in every canonical
!> nonlocal source column and in every canonical interaction supplied by the
!> caller.  A local LR-04 interaction is canonical-native and is consumed in
!> exactly the same way.
!>
!> Loss is the LR-03 Hermitian spectral product.  In canonical LR-04 storage
!> the ordinary dagger is replaced by the LR-04 metric adjoint, which is the
!> canonical representation of the ordinary physical/raw dagger:
!>
!>   L = -(chi - chi^dagger_W)/(2 i pi).
!>
!> The null-measure origin is retained in matrices for pointwise Dyson action,
!> but the loss matrix uses the zero active-space extension there.  No mode
!> finder, pair-Xi route, empirical rescaling, or pole regularization belongs
!> in this module.
!------------------------------------------------------------------------------
module tddft_dyson_mod

   use precision_mod, only: rp
   use math_mod, only: pi
   use lr_response_space_mod, only: response_space_layout, &
      response_compose_operators, response_identity_operator, response_operator_adjoint
   implicit none
   private

   character(len=*), parameter, public :: lr_dyson_route_direct_alsda = 'direct_alsda'
   character(len=*), parameter, public :: lr_dyson_route_goldstone_sumrule = 'goldstone_sumrule'
   character(len=*), parameter, public :: lr_dyson_route_direct_alsda_goldstone_corrected = &
      'direct_alsda_goldstone_corrected'
   character(len=*), parameter, public :: lr_dyson_response_representation = &
      'LR-04 canonical right-weighted B=A*W'
   character(len=*), parameter, public :: lr_dyson_compact_response_representation = &
      'orthonormal compact LMTO product representation'
   character(len=*), parameter, public :: lr_dyson_loss_convention = &
      'L=-(chi-chi^dagger_W)/(2*i*pi); canonical metric-adjoint form'
   character(len=*), parameter, public :: lr_dyson_convention = &
      'solve (I-chi_KS*K) chi=chi_KS with LAPACK zgesv; no explicit inverse'

   real(rp), parameter, public :: lr_dyson_near_singular_condition = 1.0e8_rp
   real(rp), parameter, public :: lr_dyson_ill_conditioned_condition = 1.0e12_rp

   !> Request for one q/frequency/channel response batch.
   !>
   !> `canonical_interaction` must already be in the LR-04 representation.  In
   !> particular, this request has no raw-kernel escape hatch and no automatic
   !> route selection.  For the corrected route the supplied matrix must be
   !> the explicitly selected GCR-01 corrected Kxc.
   type, public :: tddft_dyson_request
      type(response_space_layout), pointer :: response_space => null()
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = ''
      ! When true, the matrices are already in the accepted orthonormal
      ! compact product representation.  The response_space pointer is then
      ! optional and no point-space metric is inserted into Dyson or loss.
      logical :: compact_orthonormal = .false.
      complex(rp), allocatable :: ks_susceptibility(:, :, :)
      complex(rp), allocatable :: canonical_interaction(:, :)
      character(len=48) :: interaction_route = lr_dyson_route_direct_alsda
      character(len=256) :: interaction_provenance = ''
      character(len=512) :: electronic_state_provenance = ''
      character(len=256) :: response_space_metadata = ''
   end type tddft_dyson_request

   !> Enhanced response, loss matrices, and per-frequency diagnostics.
   type, public :: tddft_dyson_result
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = ''
      logical :: compact_orthonormal = .false.
      character(len=48) :: interaction_route = ''
      character(len=256) :: interaction_provenance = ''
      character(len=512) :: electronic_state_provenance = ''
      character(len=256) :: response_space_metadata = ''
      character(len=128) :: response_representation = lr_dyson_response_representation
      character(len=256) :: dyson_convention = lr_dyson_convention
      character(len=256) :: loss_convention = lr_dyson_loss_convention
      character(len=256) :: status = 'not evaluated'
      complex(rp), allocatable :: ks_susceptibility(:, :, :)
      complex(rp), allocatable :: canonical_interaction(:, :)
      complex(rp), allocatable :: denominator(:, :, :)
      complex(rp), allocatable :: enhanced_susceptibility(:, :, :)
      complex(rp), allocatable :: loss_matrix(:, :, :)
      real(rp), allocatable :: denominator_min_singular_value(:)
      real(rp), allocatable :: denominator_max_singular_value(:)
      real(rp), allocatable :: denominator_condition_number(:)
      real(rp), allocatable :: denominator_min_magnitude_eigenvalue(:)
      integer, allocatable :: solve_info(:)
      logical, allocatable :: solve_succeeded(:)
      logical, allocatable :: near_singular_collective_pole(:)
      logical, allocatable :: numerically_singular(:)
      logical, allocatable :: ill_conditioned(:)
      real(rp), allocatable :: loss_metric_hermiticity_residual(:)
      real(rp), allocatable :: dyson_residual_frobenius(:)
      real(rp), allocatable :: dyson_residual_relative(:)
      real(rp), allocatable :: dyson_residual_infinity(:)
      character(len=256), allocatable :: frequency_status(:)
   end type tddft_dyson_result

   public :: evaluate_tddft_dyson
   public :: solve_tddft_dyson_frequency
   public :: tddft_loss_matrix
   public :: loss_matrix_hermiticity_residual

   interface tddft_loss_matrix
      module procedure tddft_loss_matrix_ordinary
      module procedure tddft_loss_matrix_metric
   end interface tddft_loss_matrix

   interface loss_matrix_hermiticity_residual
      module procedure loss_matrix_hermiticity_residual_ordinary
      module procedure loss_matrix_hermiticity_residual_metric
   end interface loss_matrix_hermiticity_residual

contains

   !> Evaluate the canonical LR-04 Dyson response for all requested frequencies.
   subroutine evaluate_tddft_dyson(request, result)
      type(tddft_dyson_request), intent(in) :: request
      type(tddft_dyson_result), intent(out) :: result

      type(response_space_layout), pointer :: space
      complex(rp), allocatable :: identity(:, :), product(:, :), denominator(:, :), enhanced(:, :), loss(:, :)
      complex(rp), allocatable :: residual(:, :)
      integer :: n, nw, iw
      logical :: compact

      compact = request%compact_orthonormal
      if (compact) then
         call validate_compact_request(request, n)
      else
         if (.not. associated(request%response_space)) then
            error stop 'evaluate_tddft_dyson: response-space reference is incomplete'
         end if
         space => request%response_space
         call validate_request(request, space)
         n = space%ndim
      end if

      nw = size(request%frequencies)
      allocate(result%frequencies(nw), result%ks_susceptibility(n, n, nw), &
         result%canonical_interaction(n, n), result%denominator(n, n, nw), &
         result%enhanced_susceptibility(n, n, nw), result%loss_matrix(n, n, nw), &
         result%denominator_min_singular_value(nw), result%denominator_max_singular_value(nw), &
         result%denominator_condition_number(nw), result%denominator_min_magnitude_eigenvalue(nw), &
         result%solve_info(nw), result%solve_succeeded(nw), &
         result%near_singular_collective_pole(nw), result%numerically_singular(nw), result%ill_conditioned(nw), &
         result%loss_metric_hermiticity_residual(nw), result%dyson_residual_frobenius(nw), &
         result%dyson_residual_relative(nw), result%dyson_residual_infinity(nw), result%frequency_status(nw))

      allocate(identity(n, n), product(n, n), denominator(n, n), enhanced(n, n), loss(n, n), residual(n, n))
      if (compact) then
         identity = cmplx(0.0_rp, 0.0_rp, rp)
         do iw = 1, n
            identity(iw, iw) = cmplx(1.0_rp, 0.0_rp, rp)
         end do
      else
         call response_identity_operator(space, identity)
      end if

      result%q = request%q
      result%frequencies = request%frequencies
      result%eta = request%eta
      result%channel = trim(request%channel)
      result%compact_orthonormal = compact
      if (compact) result%response_representation = lr_dyson_compact_response_representation
      result%interaction_route = trim(request%interaction_route)
      result%interaction_provenance = trim(request%interaction_provenance)
      result%electronic_state_provenance = trim(request%electronic_state_provenance)
      result%response_space_metadata = trim(request%response_space_metadata)
      if (len_trim(result%response_space_metadata) == 0 .and. compact) then
         write (result%response_space_metadata, '(a,i0,a)') 'dimension=', n, '; orthonormal compact product space'
      else if (len_trim(result%response_space_metadata) == 0) then
         write (result%response_space_metadata, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'ndim=', n, &
            ' active_dimension=', space%active_dimension, ' nsite=', space%nsite, &
            ' radial_points=', space%npoint, ' channels=', space%nchannel
      end if
      result%ks_susceptibility = request%ks_susceptibility
      result%canonical_interaction = request%canonical_interaction

      do iw = 1, nw
         ! Both operands are canonical LR-04 operators.  This ordinary
         ! multiplication is exactly the LR-04 composition rule; no W is
         ! inserted here.
         if (compact) then
            product = matmul(request%ks_susceptibility(:, :, iw), request%canonical_interaction)
         else
            call response_compose_operators(space, request%ks_susceptibility(:, :, iw), &
               request%canonical_interaction, product)
         end if
         denominator = identity - product
         result%denominator(:, :, iw) = denominator

         call solve_tddft_denominator(denominator, request%ks_susceptibility(:, :, iw), enhanced, result%solve_info(iw), &
            result%denominator_condition_number(iw), result%denominator_min_singular_value(iw), &
            result%denominator_max_singular_value(iw), result%numerically_singular(iw), &
            result%ill_conditioned(iw))
         call minimum_magnitude_eigenvalue(denominator, result%denominator_min_magnitude_eigenvalue(iw))
         result%enhanced_susceptibility(:, :, iw) = enhanced
         result%solve_succeeded(iw) = result%solve_info(iw) == 0 .and. .not. result%numerically_singular(iw)
         result%near_singular_collective_pole(iw) = result%solve_info(iw) == 0 .and. &
            result%denominator_condition_number(iw) >= lr_dyson_near_singular_condition .and. &
            .not. result%numerically_singular(iw)

         if (result%solve_info(iw) /= 0) then
            result%frequency_status(iw) = 'BLOCKED: numerically singular Dyson denominator; no regularization applied'
         else if (result%numerically_singular(iw)) then
            result%frequency_status(iw) = 'BLOCKED: singular/precision-limited Dyson denominator; no regularization applied'
         else if (result%near_singular_collective_pole(iw)) then
            result%frequency_status(iw) = 'PASS: solved near-singular collective-pole denominator; no regularization applied'
         else if (result%ill_conditioned(iw)) then
            result%frequency_status(iw) = 'WARNING: solved but Dyson denominator is numerically ill-conditioned'
         else
            result%frequency_status(iw) = 'PASS: Dyson solve'
         end if

         residual = matmul(denominator, enhanced) - request%ks_susceptibility(:, :, iw)
         result%dyson_residual_frobenius(iw) = sqrt(sum(abs(residual)**2))
         result%dyson_residual_relative(iw) = result%dyson_residual_frobenius(iw)/ &
            max(sqrt(sum(abs(request%ks_susceptibility(:, :, iw))**2)), tiny(1.0_rp))
         result%dyson_residual_infinity(iw) = maxval(abs(residual))

         if (result%solve_info(iw) == 0) then
            if (compact) then
               loss = tddft_loss_matrix(enhanced)
            else
               loss = tddft_loss_matrix(space, enhanced)
            end if
         else
            loss = cmplx(0.0_rp, 0.0_rp, rp)
         end if
         result%loss_matrix(:, :, iw) = loss
         if (compact) then
            result%loss_metric_hermiticity_residual(iw) = loss_matrix_hermiticity_residual(loss)
         else
            result%loss_metric_hermiticity_residual(iw) = loss_matrix_hermiticity_residual(space, loss)
         end if
      end do

      if (any(result%solve_info /= 0) .or. any(result%numerically_singular)) then
         result%status = 'BLOCKED: at least one Dyson frequency is numerically singular'
      else if (any(result%near_singular_collective_pole)) then
         result%status = 'PASS: enhanced susceptibility and loss; near-singular poles reported'
      else if (any(result%ill_conditioned)) then
         result%status = 'WARNING: enhanced susceptibility solved; ill-conditioned frequencies reported'
      else
         result%status = 'PASS: enhanced susceptibility and loss matrix'
      end if
   end subroutine evaluate_tddft_dyson

   !> Form D=I-chi_KS*K and solve D*chi=chi_KS for one frequency.
   !>
   !> This standalone primitive is representation agnostic only in the sense
   !> that both inputs must use the same already-canonical representation.  It
   !> does not reconstruct raw kernels and does not apply a response metric.
   subroutine solve_tddft_dyson_frequency(chi_ks, kernel, chi, denominator, info, condition_number, &
                                          min_singular_value, max_singular_value)
      complex(rp), intent(in) :: chi_ks(:, :), kernel(:, :)
      complex(rp), intent(out) :: chi(:, :), denominator(:, :)
      integer, intent(out) :: info
      real(rp), intent(out), optional :: condition_number, min_singular_value, max_singular_value
      integer :: n, i

      n = size(chi_ks, 1)
      if (n < 1 .or. size(chi_ks, 2) /= n .or. any(shape(kernel) /= [n, n]) .or. &
          any(shape(chi) /= [n, n]) .or. any(shape(denominator) /= [n, n])) then
         error stop 'solve_tddft_dyson_frequency: incompatible canonical response dimensions'
      end if
      denominator = -matmul(chi_ks, kernel)
      do i = 1, n
         denominator(i, i) = denominator(i, i) + cmplx(1.0_rp, 0.0_rp, rp)
      end do
      call solve_tddft_denominator(denominator, chi_ks, chi, info, condition_number, min_singular_value, &
         max_singular_value)
   end subroutine solve_tddft_dyson_frequency

   !> LR-03 loss product for an ordinary orthonormal response matrix.
   !>
   !> The no-space overload is retained for isolated algebraic fixtures.  A
   !> production LR-04 caller must use the overload with `response_space`.
   function tddft_loss_matrix_ordinary(chi) result(loss)
      complex(rp), intent(in) :: chi(:, :)
      complex(rp) :: loss(size(chi, 1), size(chi, 2))
      integer :: n, i, j

      n = size(chi, 1)
      if (size(chi, 2) /= n .or. n < 1) error stop 'tddft_loss_matrix: chi must be nonempty and square'
      do j = 1, n
         do i = 1, n
            loss(i, j) = -(chi(i, j) - conjg(chi(j, i)))/cmplx(0.0_rp, 2.0_rp*pi, rp)
         end do
      end do
   end function tddft_loss_matrix_ordinary

   !> LR-03 loss product in canonical LR-04 storage.
   !>
   !> `response_operator_adjoint` is the canonical metric-adjoint.  The
   !> origin row/column is set to the zero active-space extension because its
   !> volume measure is exactly zero and it has no physical loss weight.
   function tddft_loss_matrix_metric(space, chi) result(loss)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: chi(:, :)
      complex(rp) :: loss(size(chi, 1), size(chi, 2))
      complex(rp) :: adjoint(size(chi, 1), size(chi, 2))
      integer :: n, i, j

      n = space%ndim
      if (size(chi, 1) /= n .or. size(chi, 2) /= n) then
         error stop 'tddft_loss_matrix: canonical chi/response-space shape mismatch'
      end if
      call response_operator_adjoint(space, chi, adjoint)
      loss = -(chi - adjoint)/cmplx(0.0_rp, 2.0_rp*pi, rp)
      do i = 1, n
         if (space%metric_weights(i) == 0.0_rp) then
            do j = 1, n
               loss(i, j) = cmplx(0.0_rp, 0.0_rp, rp)
               loss(j, i) = cmplx(0.0_rp, 0.0_rp, rp)
            end do
         end if
      end do
   end function tddft_loss_matrix_metric

   pure real(rp) function loss_matrix_hermiticity_residual_ordinary(loss) result(residual)
      complex(rp), intent(in) :: loss(:, :)
      integer :: n

      n = size(loss, 1)
      if (size(loss, 2) /= n .or. n < 1) error stop 'loss_matrix_hermiticity_residual: loss must be square'
      residual = maxval(abs(loss - conjg(transpose(loss))))
   end function loss_matrix_hermiticity_residual_ordinary

   real(rp) function loss_matrix_hermiticity_residual_metric(space, loss) result(residual)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: loss(:, :)
      complex(rp) :: adjoint(size(loss, 1), size(loss, 2))

      if (size(loss, 1) /= space%ndim .or. size(loss, 2) /= space%ndim) then
         error stop 'loss_matrix_hermiticity_residual: canonical loss/response-space shape mismatch'
      end if
      call response_operator_adjoint(space, loss, adjoint)
      residual = maxval(abs(loss - adjoint))
   end function loss_matrix_hermiticity_residual_metric

   subroutine validate_request(request, space)
      type(tddft_dyson_request), intent(in) :: request
      type(response_space_layout), intent(in) :: space
      integer :: n, nw

      if (space%ndim < 1 .or. .not. allocated(space%metric_weights) .or. &
          size(space%metric_weights) /= space%ndim) then
         error stop 'evaluate_tddft_dyson: response space is not initialized'
      end if
      if (any(space%metric_weights < 0.0_rp)) error stop 'evaluate_tddft_dyson: response metric is invalid'
      if (.not. allocated(request%frequencies) .or. size(request%frequencies) < 1) then
         error stop 'evaluate_tddft_dyson: at least one frequency is required'
      end if
      if (request%eta <= 0.0_rp) error stop 'evaluate_tddft_dyson: retarded eta must be positive'
      if (len_trim(request%channel) == 0) error stop 'evaluate_tddft_dyson: circular channel provenance is required'
      select case (trim(request%channel))
      case ('plus', 'minus', 'chi_plus', 'chi_minus')
      case default
         error stop 'evaluate_tddft_dyson: channel must be plus/minus or chi_plus/chi_minus'
      end select
      select case (trim(request%interaction_route))
      case (lr_dyson_route_direct_alsda)
         if (index(trim(request%interaction_provenance), 'KXC-01') /= 1) then
            error stop 'evaluate_tddft_dyson: direct_alsda requires KXC-01 interaction provenance'
         end if
      case (lr_dyson_route_goldstone_sumrule)
         if (index(trim(request%interaction_provenance), 'GSR-01') /= 1) then
            error stop 'evaluate_tddft_dyson: goldstone_sumrule requires GSR-01 interaction provenance'
         end if
      case (lr_dyson_route_direct_alsda_goldstone_corrected)
         if (index(trim(request%interaction_provenance), 'GCR-01') /= 1) then
            error stop 'evaluate_tddft_dyson: corrected route requires GCR-01 interaction provenance'
         end if
      case default
         error stop 'evaluate_tddft_dyson: interaction route is not an explicit supported choice'
      end select
      if (len_trim(request%interaction_provenance) == 0) then
         error stop 'evaluate_tddft_dyson: interaction provenance is required'
      end if
      if (len_trim(request%electronic_state_provenance) == 0) then
         error stop 'evaluate_tddft_dyson: electronic-state provenance is required'
      end if
      if (.not. allocated(request%ks_susceptibility) .or. .not. allocated(request%canonical_interaction)) then
         error stop 'evaluate_tddft_dyson: canonical chiKS and interaction are required'
      end if
      n = space%ndim
      nw = size(request%frequencies)
      if (any(shape(request%ks_susceptibility) /= [n, n, nw])) then
         error stop 'evaluate_tddft_dyson: canonical chiKS/frequency shape mismatch'
      end if
      if (any(shape(request%canonical_interaction) /= [n, n])) then
         error stop 'evaluate_tddft_dyson: canonical interaction shape mismatch'
      end if
   end subroutine validate_request

   subroutine validate_compact_request(request, n)
      type(tddft_dyson_request), intent(in) :: request
      integer, intent(out) :: n

      if (.not. allocated(request%frequencies) .or. size(request%frequencies) < 1) then
         error stop 'evaluate_tddft_dyson: at least one frequency is required'
      end if
      if (request%eta <= 0.0_rp) error stop 'evaluate_tddft_dyson: retarded eta must be positive'
      if (len_trim(request%channel) == 0) error stop 'evaluate_tddft_dyson: circular channel provenance is required'
      select case (trim(request%channel))
      case ('plus', 'minus', 'chi_plus', 'chi_minus')
      case default
         error stop 'evaluate_tddft_dyson: channel must be plus/minus or chi_plus/chi_minus'
      end select
      select case (trim(request%interaction_route))
      case (lr_dyson_route_direct_alsda)
         if (index(trim(request%interaction_provenance), 'KXC-01') /= 1) then
            error stop 'evaluate_tddft_dyson: direct_alsda requires KXC-01 interaction provenance'
         end if
      case (lr_dyson_route_goldstone_sumrule)
         if (index(trim(request%interaction_provenance), 'GSR-01') /= 1) then
            error stop 'evaluate_tddft_dyson: goldstone_sumrule requires GSR-01 interaction provenance'
         end if
      case (lr_dyson_route_direct_alsda_goldstone_corrected)
         if (index(trim(request%interaction_provenance), 'GCR-01') /= 1) then
            error stop 'evaluate_tddft_dyson: corrected route requires GCR-01 interaction provenance'
         end if
      case default
         error stop 'evaluate_tddft_dyson: interaction route is not an explicit supported choice'
      end select
      if (len_trim(request%interaction_provenance) == 0) then
         error stop 'evaluate_tddft_dyson: interaction provenance is required'
      end if
      if (len_trim(request%electronic_state_provenance) == 0) then
         error stop 'evaluate_tddft_dyson: electronic-state provenance is required'
      end if
      if (.not. allocated(request%ks_susceptibility) .or. .not. allocated(request%canonical_interaction)) then
         error stop 'evaluate_tddft_dyson: compact chiKS and interaction are required'
      end if
      n = size(request%canonical_interaction, 1)
      if (n < 1 .or. size(request%canonical_interaction, 2) /= n) then
         error stop 'evaluate_tddft_dyson: compact interaction must be nonempty and square'
      end if
      if (any(shape(request%ks_susceptibility) /= [n, n, size(request%frequencies)])) then
         error stop 'evaluate_tddft_dyson: compact chiKS/frequency shape mismatch'
      end if
   end subroutine validate_compact_request

   !> Solve one already-formed Dyson denominator and report its singular values.
   subroutine solve_tddft_denominator(denominator, rhs, solution, info, condition_number, min_singular_value, &
                                      max_singular_value, numerically_singular, ill_conditioned)
      complex(rp), intent(in) :: denominator(:, :), rhs(:, :)
      complex(rp), intent(out) :: solution(:, :)
      integer, intent(out) :: info
      real(rp), intent(out), optional :: condition_number, min_singular_value, max_singular_value
      logical, intent(out), optional :: numerically_singular, ill_conditioned
      complex(rp), allocatable :: system(:, :), rhs_work(:, :)
      integer, allocatable :: pivots(:)
      real(rp) :: cond, smin, smax
      logical :: precision_limited, poorly_conditioned
      integer :: n
      external :: zgesv

      n = size(denominator, 1)
      if (n < 1 .or. size(denominator, 2) /= n .or. any(shape(rhs) /= [n, n]) .or. &
          any(shape(solution) /= [n, n])) then
         error stop 'solve_tddft_denominator: incompatible matrix dimensions'
      end if

      allocate(system(n, n), rhs_work(n, n), pivots(n))
      system = denominator
      rhs_work = rhs
      call singular_value_diagnostics(system, cond, smin, smax, info)
      precision_limited = info /= 0 .or. smax == 0.0_rp .or. smin <= &
         100.0_rp*epsilon(1.0_rp)*smax
      poorly_conditioned = cond >= lr_dyson_ill_conditioned_condition

      ! The SVD above operates on a private copy.  zgesv therefore receives the
      ! original denominator and all RHS columns, not a reconstructed inverse.
      system = denominator
      if (precision_limited) then
         info = 1
         solution = cmplx(0.0_rp, 0.0_rp, rp)
      else
         call zgesv(n, n, system, n, pivots, rhs_work, n, info)
         if (info == 0) then
            solution = rhs_work
         else
            solution = cmplx(0.0_rp, 0.0_rp, rp)
         end if
      end if

      if (present(condition_number)) condition_number = cond
      if (present(min_singular_value)) min_singular_value = smin
      if (present(max_singular_value)) max_singular_value = smax
      if (present(numerically_singular)) numerically_singular = precision_limited
      if (present(ill_conditioned)) ill_conditioned = poorly_conditioned
   end subroutine solve_tddft_denominator

   subroutine singular_value_diagnostics(matrix, condition_number, min_singular_value, max_singular_value, info)
      complex(rp), intent(in) :: matrix(:, :)
      real(rp), intent(out) :: condition_number, min_singular_value, max_singular_value
      integer, intent(out) :: info
      complex(rp), allocatable :: work_matrix(:,:), work(:)
      complex(rp) :: work_query(1), u_dummy(1, 1), vt_dummy(1, 1)
      real(rp), allocatable :: singular_values(:), rwork(:)
      integer :: n, lwork, svd_info
      external :: zgesvd

      n = size(matrix, 1)
      if (size(matrix, 2) /= n .or. n < 1) error stop 'singular_value_diagnostics: matrix must be nonempty and square'
      allocate(work_matrix(n, n), singular_values(n), rwork(5*n))
      work_matrix = matrix
      call zgesvd('N', 'N', n, n, work_matrix, n, singular_values, u_dummy, 1, vt_dummy, 1, work_query, -1, rwork, svd_info)
      if (svd_info /= 0) then
         info = svd_info
         condition_number = huge(1.0_rp)
         min_singular_value = 0.0_rp
         max_singular_value = 0.0_rp
         return
      end if
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      work_matrix = matrix
      call zgesvd('N', 'N', n, n, work_matrix, n, singular_values, u_dummy, 1, vt_dummy, 1, work, lwork, rwork, svd_info)
      info = svd_info
      if (info /= 0) then
         condition_number = huge(1.0_rp)
         min_singular_value = 0.0_rp
         max_singular_value = 0.0_rp
         return
      end if
      max_singular_value = maxval(singular_values)
      min_singular_value = minval(singular_values)
      if (min_singular_value > 0.0_rp) then
         condition_number = max_singular_value/min_singular_value
      else
         condition_number = huge(1.0_rp)
      end if
   end subroutine singular_value_diagnostics

   !> Return the smallest eigenvalue magnitude of a denominator for the raw
   !> static/pole diagnostic.  This is diagnostic only; the solve continues to
   !> use the unmodified denominator and LAPACK zgesv.
   subroutine minimum_magnitude_eigenvalue(matrix, minimum_magnitude)
      complex(rp), intent(in) :: matrix(:, :)
      real(rp), intent(out) :: minimum_magnitude
      complex(rp), allocatable :: work_matrix(:, :), eigenvalues(:), work(:)
      complex(rp) :: work_query(1), vl_dummy(1, 1), vr_dummy(1, 1)
      real(rp), allocatable :: rwork(:)
      integer :: n, lwork, info
      external :: zgeev

      n = size(matrix, 1)
      if (size(matrix, 2) /= n .or. n < 1) error stop 'minimum_magnitude_eigenvalue: matrix must be nonempty and square'
      allocate(work_matrix(n, n), eigenvalues(n), rwork(max(1, 2*n)))
      work_matrix = matrix
      call zgeev('N', 'N', n, work_matrix, n, eigenvalues, vl_dummy, 1, vr_dummy, 1, work_query, -1, rwork, info)
      if (info /= 0) then
         minimum_magnitude = huge(1.0_rp)
         deallocate(work_matrix, eigenvalues, rwork)
         return
      end if
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      work_matrix = matrix
      call zgeev('N', 'N', n, work_matrix, n, eigenvalues, vl_dummy, 1, vr_dummy, 1, work, lwork, rwork, info)
      if (info == 0) then
         minimum_magnitude = minval(abs(eigenvalues))
      else
         minimum_magnitude = huge(1.0_rp)
      end if
      deallocate(work_matrix, eigenvalues, work, rwork)
   end subroutine minimum_magnitude_eigenvalue

end module tddft_dyson_mod
