!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Compact LR-04 interaction/vector projection and independent GSR route.
!>
!> The retained product modes are the Euclidean modes of W**(1/2) times a
!> point-space product field.  If U denotes those stored weighted modes, the
!> compact coordinates and the reconstructed point field are
!>
!>   c = U^H W**(1/2) x,       x = W**(-1/2) U c.
!>
!> Consequently a local pointwise operator K maps as U^H K U.  The radial
!> volume metric is present in the vector projection/reconstruction; it is not
!> inserted a second time into the local operator.  This module contains no
!> ALSDA call and no direct-ALSDA fallback.  Its sum-rule solve is an
!> independent compact local-interaction solve.
!------------------------------------------------------------------------------
module lr_compact_static_interaction_mod

   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use response_angular_basis_mod, only: response_angular_pi
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex
   use lr_response_space_mod, only: response_space_layout
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis
   implicit none
   private

   real(rp), parameter, public :: lr_compact_gsr_residual_tolerance = 1.0e-9_rp
   real(rp), parameter, public :: lr_compact_gsr_action_tolerance = 1.0e-10_rp
   character(len=*), parameter, public :: lr_compact_representation = &
      'weighted-orthonormal LMTO product representation; U=sqrt(W)V'
   character(len=*), parameter, public :: lr_compact_mapping_contract = &
      'c=U^H sqrt(W) x, x=inv(sqrt(W)) U c, Kc=U^H K U'

   !> Independent compact LCMM/GSR result.  The unknown is a local scalar
   !> U_LCMM represented in the retained L=0 radial product span.  The solve
   !> equation is rows-by-unknowns, Gamma*u=target, and is deliberately
   !> allowed to return BLOCKED when the full compact response cannot satisfy
   !> the spherical local sum-rule ansatz.
   type, public :: lr_compact_gsr_result
      complex(rp), allocatable :: u_lcmm(:, :)       ! site, radial point
      complex(rp), allocatable :: effective_field(:, :) ! 4*pi*U*m_z
      complex(rp), allocatable :: canonical_interaction(:, :) ! compact local K_eff
      complex(rp), allocatable :: magnetization_response(:)
      complex(rp), allocatable :: generated_field(:)
      complex(rp), allocatable :: generated_response(:)
      complex(rp), allocatable :: residual_vector(:)
      complex(rp), allocatable :: gamma(:, :)
      complex(rp), allocatable :: solution(:)
      real(rp), allocatable :: singular_values(:)
      integer :: equation_rows = 0
      integer :: unknowns = 0
      integer :: rank = 0
      real(rp) :: condition_number = huge(1.0_rp)
      real(rp) :: svd_rcond = -1.0_rp
      real(rp) :: svd_cutoff = huge(1.0_rp)
      real(rp) :: coefficient_norm = huge(1.0_rp)
      real(rp) :: equation_residual_norm = huge(1.0_rp)
      real(rp) :: equation_relative_residual = huge(1.0_rp)
      real(rp) :: residual_norm = huge(1.0_rp)
      real(rp) :: relative_residual = huge(1.0_rp)
      real(rp) :: residual_difference_norm = huge(1.0_rp)
      real(rp) :: residual_difference_relative = huge(1.0_rp)
      real(rp) :: residual_difference_max_component = huge(1.0_rp)
      real(rp) :: max_residual_component = huge(1.0_rp)
      real(rp) :: assembled_action_sum_norm = huge(1.0_rp)
      real(rp) :: assembled_action_matrix_norm = huge(1.0_rp)
      real(rp) :: assembled_action_difference_norm = huge(1.0_rp)
      real(rp) :: assembled_action_difference_relative = huge(1.0_rp)
      real(rp) :: assembled_action_difference_max_component = huge(1.0_rp)
      real(rp) :: magnetization_weighted_norm = huge(1.0_rp)
      real(rp) :: magnetization_projection_residual_norm = huge(1.0_rp)
      real(rp) :: magnetization_projection_relative_residual = huge(1.0_rp)
      complex(rp) :: rigid_overlap = cmplx(0.0_rp, 0.0_rp, rp)
      logical :: rank_deficient = .true.
      logical :: blocked = .true.
      character(len=64) :: status = 'not evaluated'
   end type lr_compact_gsr_result

   public :: compact_project_point_vector
   public :: compact_reconstruct_point_vector
   public :: compact_project_local_operator
   public :: compact_apply_local_operator
   public :: compact_project_magnetization
   public :: compact_weighted_projection_diagnostics
   public :: evaluate_compact_goldstone_sumrule

   interface compact_project_local_operator
      module procedure compact_project_local_operator_complex
   end interface compact_project_local_operator

   interface
      subroutine zgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info)
         import :: rp
         integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
         complex(rp), intent(inout) :: a(lda, *), b(ldb, *), work(*)
         real(rp), intent(out) :: s(*), rwork(*)
         real(rp), intent(in) :: rcond
         integer, intent(out) :: rank, info
      end subroutine zgelss
   end interface

contains

   !> Project a point-space vector with the actual weighted modes stored by the
   !> product basis.  The origin is an exact zero-measure row and contributes
   !> zero without division or regularization.
   subroutine compact_project_point_vector(space, product, point_vector, compact_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: point_vector(:)
      complex(rp), intent(out) :: compact_vector(:)
      integer :: site, response_l, response_m, mode, ir, point_flat, compact_flat
      real(rp) :: sqrt_weight
      type(response_super_index) :: item

      call validate_point_to_compact_shapes(space, product, point_vector, compact_vector)
      compact_vector = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, product%nsite
         do response_l = 0, product%response_lmax
            do response_m = -response_l, response_l
               do mode = 1, product%blocks(site, response_l)%rank
                  compact_flat = product%flat_index(site, response_l, response_m, mode)
                  do ir = 2, space%npoint
                     sqrt_weight = sqrt(space%radial_weights(ir))
                     item = response_super_index(site, response_l, response_m, ir, 1)
                     call response_flatten_superindex(item, space%nsite, space%response_lmax, &
                        space%npoint, space%nchannel, point_flat)
                     compact_vector(compact_flat) = compact_vector(compact_flat) + &
                        conjg(product%blocks(site, response_l)%weighted_modes(ir, mode))*sqrt_weight*point_vector(point_flat)
                  end do
               end do
            end do
         end do
      end do
   end subroutine compact_project_point_vector

   !> Reconstruct the retained product-space component of a compact vector.
   !> Null-measure origin values are set to the finite zero extension.
   subroutine compact_reconstruct_point_vector(space, product, compact_vector, point_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: compact_vector(:)
      complex(rp), intent(out) :: point_vector(:)
      integer :: site, response_l, response_m, mode, ir, point_flat, compact_flat
      real(rp) :: sqrt_weight
      type(response_super_index) :: item

      call validate_compact_to_point_shapes(space, product, compact_vector, point_vector)
      point_vector = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, product%nsite
         do response_l = 0, product%response_lmax
            do response_m = -response_l, response_l
               do mode = 1, product%blocks(site, response_l)%rank
                  compact_flat = product%flat_index(site, response_l, response_m, mode)
                  do ir = 2, space%npoint
                     sqrt_weight = sqrt(space%radial_weights(ir))
                     if (sqrt_weight <= 0.0_rp) cycle
                     item = response_super_index(site, response_l, response_m, ir, 1)
                     call response_flatten_superindex(item, space%nsite, space%response_lmax, &
                        space%npoint, space%nchannel, point_flat)
                     point_vector(point_flat) = point_vector(point_flat) + &
                        product%blocks(site, response_l)%weighted_modes(ir, mode)/sqrt_weight*compact_vector(compact_flat)
                  end do
               end do
            end do
         end do
      end do
   end subroutine compact_reconstruct_point_vector

   !> Project a spherical local pointwise operator.  Since the point basis is
   !> V=W**(-1/2)U, the W factors from the metric projection and reconstruction
   !> cancel exactly, leaving U^H K U.
   subroutine compact_project_local_operator_complex(space, product, values, compact_operator)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: values(:, :)
      complex(rp), intent(out) :: compact_operator(:, :)
      integer :: site, response_l, response_m, mode_out, mode_in, ir
      integer :: flat_out, flat_in

      call validate_operator_shapes(space, product, compact_operator)
      if (any(shape(values) /= [space%nsite, space%npoint])) then
         error stop 'compact_project_local_operator: site/radial scalar shape mismatch'
      end if
      compact_operator = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, product%nsite
         do response_l = 0, product%response_lmax
            do response_m = -response_l, response_l
               do mode_out = 1, product%blocks(site, response_l)%rank
                  flat_out = product%flat_index(site, response_l, response_m, mode_out)
                  do mode_in = 1, product%blocks(site, response_l)%rank
                     flat_in = product%flat_index(site, response_l, response_m, mode_in)
                     do ir = 2, space%npoint
                        compact_operator(flat_out, flat_in) = compact_operator(flat_out, flat_in) + &
                           conjg(product%blocks(site, response_l)%weighted_modes(ir, mode_out))*values(site, ir)* &
                           product%blocks(site, response_l)%weighted_modes(ir, mode_in)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine compact_project_local_operator_complex

   subroutine compact_apply_local_operator(space, product, values, vector, result)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: values(:, :), vector(:)
      complex(rp), intent(out) :: result(:)
      complex(rp), allocatable :: operator(:, :)

      call validate_common_layout(space, product)
      if (size(vector) /= product%product_dimension .or. size(result) /= product%product_dimension) then
         error stop 'compact local operator: vector shape mismatch'
      end if
      allocate(operator(product%product_dimension, product%product_dimension))
      call compact_project_local_operator(space, product, values, operator)
      result = matmul(operator, vector)
   end subroutine compact_apply_local_operator

   !> Build the normalized L=0 magnetization vector used by both static
   !> diagnostics.  The point-space coefficient is m_00=sqrt(4*pi)*m_z.
   subroutine compact_project_magnetization(space, product, magnetization, compact_vector, point_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      real(rp), intent(in) :: magnetization(:, :)
      complex(rp), intent(out) :: compact_vector(:)
      complex(rp), allocatable, intent(out), optional :: point_vector(:)
      complex(rp), allocatable :: point(:)
      integer :: site, ir, flat
      type(response_super_index) :: item

      if (any(shape(magnetization) /= [space%nsite, space%npoint])) then
         error stop 'compact_project_magnetization: magnetization shape mismatch'
      end if
      allocate(point(space%ndim))
      point = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, space%nsite
         do ir = 1, space%npoint
            item = response_super_index(site, 0, 0, ir, 1)
            call response_flatten_superindex(item, space%nsite, space%response_lmax, space%npoint, &
               space%nchannel, flat)
            point(flat) = cmplx(sqrt(4.0_rp*response_angular_pi)*magnetization(site, ir), 0.0_rp, rp)
         end do
      end do
      call compact_project_point_vector(space, product, point, compact_vector)
      if (present(point_vector)) then
         allocate(point_vector(space%ndim))
         point_vector = point
      end if
   end subroutine compact_project_magnetization

   !> Report the weighted point-space norm of a vector and its defect from a
   !> compact reconstruction.  This is the same diagonal radial metric used
   !> by compact_project_point_vector and compact_reconstruct_point_vector.
   subroutine compact_weighted_projection_diagnostics(space, point_vector, projected_point_vector, weighted_norm, &
                                                      residual_norm, relative_residual)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: point_vector(:), projected_point_vector(:)
      real(rp), intent(out) :: weighted_norm, residual_norm, relative_residual
      integer :: site, response_l, response_m, ir, channel, flat
      real(rp) :: weight
      type(response_super_index) :: item

      if (size(point_vector) /= space%ndim .or. size(projected_point_vector) /= space%ndim) then
         error stop 'compact projection diagnostics: point vector shape mismatch'
      end if
      weighted_norm = 0.0_rp
      residual_norm = 0.0_rp
      do site = 1, space%nsite
         do response_l = 0, space%response_lmax
            do response_m = -response_l, response_l
               do ir = 1, space%npoint
                  weight = space%radial_weights(ir)
                  if (weight <= 0.0_rp) cycle
                  do channel = 1, space%nchannel
                     item = response_super_index(site, response_l, response_m, ir, channel)
                     call response_flatten_superindex(item, space%nsite, space%response_lmax, &
                        space%npoint, space%nchannel, flat)
                     weighted_norm = weighted_norm + weight*abs(point_vector(flat))**2
                     residual_norm = residual_norm + weight*abs(point_vector(flat) - projected_point_vector(flat))**2
                  end do
               end do
            end do
         end do
      end do
      weighted_norm = sqrt(weighted_norm)
      residual_norm = sqrt(residual_norm)
      relative_residual = residual_norm/max(weighted_norm, tiny(1.0_rp))
   end subroutine compact_weighted_projection_diagnostics

   !> Solve the independent compact local LCMM/GSR equation.  The published
   !> U_LCMM radial function is represented in the retained L=0 product span;
   !> each unknown basis function is reconstructed to point space and mapped
   !> through the normalized m_00 source.  No ALSDA quantity enters Gamma.
   subroutine evaluate_compact_goldstone_sumrule(space, product, static_susceptibility, magnetization, result)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: static_susceptibility(:, :)
      real(rp), intent(in) :: magnetization(:, :)
      type(lr_compact_gsr_result), intent(out) :: result
      complex(rp), allocatable :: target(:), target_point(:), target_point_compact(:), field_basis(:, :), gamma(:, :)
      complex(rp), allocatable :: point_u(:), point_field(:), values(:, :), field(:), response(:), residual(:)
      complex(rp), allocatable :: assembled_field(:), solve_residual(:), residual_difference(:)
      complex(rp), allocatable :: a_work(:, :), b_work(:, :), work(:), work_query(:)
      real(rp), allocatable :: rwork(:)
      integer, allocatable :: unknown_site(:), unknown_mode(:)
      integer :: nrow, nunknown, l0_rank, site, mode, ir, unknown, info, lwork, i
      integer :: ldb, min_dimension, point_flat
      real(rp), parameter :: rcond = -1.0_rp
      real(rp) :: four_pi, largest_singular, smallest_positive, target_norm
      real(rp) :: sqrt_weight
      type(response_super_index) :: item

      call validate_gsr_inputs(space, product, static_susceptibility, magnetization)
      nrow = product%product_dimension
      l0_rank = 0
      do site = 1, product%nsite
         l0_rank = l0_rank + product%blocks(site, 0)%rank
      end do
      nunknown = l0_rank
      if (nunknown < 1) error stop 'evaluate_compact_goldstone_sumrule: no retained L=0 interaction modes'
      four_pi = 4.0_rp*response_angular_pi
      allocate(target(nrow), target_point(space%ndim), target_point_compact(space%ndim), field_basis(nrow, nunknown), gamma(nrow, nunknown), &
         unknown_site(nunknown), unknown_mode(nunknown), point_u(space%ndim), point_field(space%ndim), &
         values(space%nsite, space%npoint), field(nrow), response(nrow), residual(nrow), assembled_field(nrow), &
         solve_residual(nrow), residual_difference(nrow))
      call compact_project_magnetization(space, product, magnetization, target, target_point)
      call compact_reconstruct_point_vector(space, product, target, target_point_compact)
      call compact_weighted_projection_diagnostics(space, target_point, target_point_compact, result%magnetization_weighted_norm, &
         result%magnetization_projection_residual_norm, result%magnetization_projection_relative_residual)
      field_basis = cmplx(0.0_rp, 0.0_rp, rp)
      unknown = 0
      do site = 1, product%nsite
         do mode = 1, product%blocks(site, 0)%rank
            unknown = unknown + 1
            unknown_site(unknown) = site
            unknown_mode(unknown) = mode
            point_field = cmplx(0.0_rp, 0.0_rp, rp)
            do ir = 2, space%npoint
               sqrt_weight = sqrt(space%radial_weights(ir))
               if (sqrt_weight <= 0.0_rp) cycle
               item = response_super_index(site, 0, 0, ir, 1)
               call response_flatten_superindex(item, space%nsite, space%response_lmax, &
                  space%npoint, space%nchannel, point_flat)
               ! U_LCMM = weighted-mode/sqrt(W) times the compact unknown.
               ! The GSR column must act on R*c, not on the unprojected
               ! physical magnetization.  This makes it the same action as
               ! Kc*c used after the solve.
               point_field(point_flat) = four_pi*product%blocks(site, 0)%weighted_modes(ir, mode)/sqrt_weight * &
                  target_point_compact(point_flat)
            end do
            call compact_project_point_vector(space, product, point_field, field_basis(:, unknown))
         end do
      end do
      gamma = matmul(static_susceptibility, field_basis)

      ldb = max(1, max(nrow, nunknown))
      min_dimension = min(nrow, nunknown)
      allocate(a_work(nrow, nunknown), b_work(ldb, 1), result%solution(nunknown), &
         result%singular_values(min_dimension), rwork(max(1, 5*min_dimension)), work_query(1))
      a_work = gamma
      b_work = cmplx(0.0_rp, 0.0_rp, rp)
      b_work(1:nrow, 1) = target
      call zgelss(nrow, nunknown, 1, a_work, nrow, b_work, ldb, result%singular_values, rcond, result%rank, &
         work_query, -1, rwork, info)
      if (info /= 0) error stop 'evaluate_compact_goldstone_sumrule: ZGELSS workspace query failed'
      lwork = max(1, int(real(work_query(1), rp)))
      allocate(work(lwork))
      a_work = gamma
      b_work = cmplx(0.0_rp, 0.0_rp, rp)
      b_work(1:nrow, 1) = target
      call zgelss(nrow, nunknown, 1, a_work, nrow, b_work, ldb, result%singular_values, rcond, result%rank, &
         work, lwork, rwork, info)
      if (info /= 0) error stop 'evaluate_compact_goldstone_sumrule: ZGELSS failed'
      result%solution = b_work(1:nunknown, 1)
      result%coefficient_norm = sqrt(sum(abs(result%solution)**2))

      allocate(result%gamma(nrow, nunknown), result%magnetization_response(nrow), result%generated_field(nrow), &
         result%generated_response(nrow), result%residual_vector(nrow), result%u_lcmm(space%nsite, space%npoint), &
         result%effective_field(space%nsite, space%npoint), result%canonical_interaction(nrow, nrow))
      result%gamma = gamma
      result%magnetization_response = target
      result%u_lcmm = cmplx(0.0_rp, 0.0_rp, rp)
      result%effective_field = cmplx(0.0_rp, 0.0_rp, rp)
      point_u = cmplx(0.0_rp, 0.0_rp, rp)
      do unknown = 1, nunknown
         site = unknown_site(unknown)
         mode = unknown_mode(unknown)
         do ir = 2, space%npoint
            sqrt_weight = sqrt(space%radial_weights(ir))
            if (sqrt_weight <= 0.0_rp) cycle
            item = response_super_index(site, 0, 0, ir, 1)
            call response_flatten_superindex(item, space%nsite, space%response_lmax, &
               space%npoint, space%nchannel, point_flat)
            point_u(point_flat) = point_u(point_flat) + &
               product%blocks(site, 0)%weighted_modes(ir, mode)/sqrt_weight*result%solution(unknown)
         end do
      end do
      values = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, space%nsite
         do ir = 2, space%npoint
            item = response_super_index(site, 0, 0, ir, 1)
            call response_flatten_superindex(item, space%nsite, space%response_lmax, &
               space%npoint, space%nchannel, point_flat)
            result%u_lcmm(site, ir) = point_u(point_flat)
            result%effective_field(site, ir) = four_pi*result%u_lcmm(site, ir)*magnetization(site, ir)
            values(site, ir) = four_pi*result%u_lcmm(site, ir)
         end do
      end do
      call compact_project_local_operator(space, product, values, result%canonical_interaction)
      assembled_field = cmplx(0.0_rp, 0.0_rp, rp)
      do unknown = 1, nunknown
         assembled_field = assembled_field + result%solution(unknown)*field_basis(:, unknown)
      end do
      result%generated_field = matmul(result%canonical_interaction, target)
      result%assembled_action_sum_norm = sqrt(sum(abs(assembled_field)**2))
      result%assembled_action_matrix_norm = sqrt(sum(abs(result%generated_field)**2))
      result%assembled_action_difference_norm = sqrt(sum(abs(assembled_field - result%generated_field)**2))
      result%assembled_action_difference_relative = result%assembled_action_difference_norm / &
         max(result%assembled_action_matrix_norm, tiny(1.0_rp))
      result%assembled_action_difference_max_component = maxval(abs(assembled_field - result%generated_field))
      result%generated_response = matmul(static_susceptibility, result%generated_field)
      residual = result%generated_response - target
      result%residual_vector = residual
      solve_residual = matmul(gamma, result%solution) - target
      residual_difference = solve_residual - residual
      result%equation_residual_norm = sqrt(sum(abs(solve_residual)**2))
      target_norm = sqrt(sum(abs(target)**2))
      result%equation_relative_residual = result%equation_residual_norm/max(target_norm, tiny(1.0_rp))
      result%residual_norm = sqrt(sum(abs(residual)**2))
      result%relative_residual = result%residual_norm/max(target_norm, tiny(1.0_rp))
      result%residual_difference_norm = sqrt(sum(abs(residual_difference)**2))
      ! Scale the residual difference by the larger of the equation right-hand
      ! side and the assembled action.  When both residuals are at roundoff,
      ! scaling by either residual itself would manufacture a meaningless O(1)
      ! relative difference; for an ill-conditioned solve, the assembled-action
      ! scale also exposes roundoff amplified by the large coefficients.
      result%residual_difference_relative = result%residual_difference_norm / &
         max(max(target_norm, result%assembled_action_matrix_norm), tiny(1.0_rp))
      result%residual_difference_max_component = maxval(abs(residual_difference))
      result%max_residual_component = maxval(abs(residual))
      result%rigid_overlap = sum(conjg(residual)*target)
      largest_singular = maxval(result%singular_values)
      smallest_positive = huge(1.0_rp)
      do i = 1, size(result%singular_values)
         if (result%singular_values(i) > 0.0_rp) smallest_positive = min(smallest_positive, result%singular_values(i))
      end do
      if (result%rank == nunknown .and. smallest_positive < huge(1.0_rp)) then
         result%condition_number = largest_singular/smallest_positive
      end if
      result%svd_rcond = rcond
      result%svd_cutoff = epsilon(1.0_rp)*largest_singular
      result%equation_rows = nrow
      result%unknowns = nunknown
      result%rank_deficient = result%rank < nunknown
      result%blocked = result%rank_deficient .or. result%assembled_action_difference_relative > lr_compact_gsr_action_tolerance .or. &
         result%equation_relative_residual > lr_compact_gsr_residual_tolerance .or. &
         result%relative_residual > lr_compact_gsr_residual_tolerance
      if (result%rank_deficient) then
         result%status = 'BLOCKED: GSR rank/conditioning'
      else if (result%assembled_action_difference_relative > lr_compact_gsr_action_tolerance) then
         result%status = 'BLOCKED: compact GSR operator action'
      else if (result%equation_relative_residual > lr_compact_gsr_residual_tolerance) then
         result%status = 'BLOCKED: GSR rank/conditioning'
      else if (result%relative_residual > lr_compact_gsr_residual_tolerance) then
         result%status = 'BLOCKED: GSR rank/conditioning'
      else
         result%blocked = .false.
         result%status = 'PASS: independent compact LCMM sum-rule identity'
      end if
   end subroutine evaluate_compact_goldstone_sumrule

   subroutine validate_common_layout(space, product)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product

      if (space%nchannel /= 1 .or. space%nsite /= product%nsite .or. space%response_lmax /= product%response_lmax) then
         error stop 'compact projection: response/product layout mismatch'
      end if
      if (.not. allocated(space%radial_weights) .or. size(space%radial_weights) /= space%npoint) then
         error stop 'compact projection: LR-04 radial metric is unavailable'
      end if
      if (.not. allocated(product%blocks) .or. product%product_dimension < 1) then
         error stop 'compact projection: product representation is unavailable'
      end if
   end subroutine validate_common_layout

   subroutine validate_point_to_compact_shapes(space, product, point_vector, compact_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: point_vector(:)
      complex(rp), intent(out) :: compact_vector(:)

      call validate_common_layout(space, product)
      if (size(point_vector) /= space%ndim .or. size(compact_vector) /= product%product_dimension) then
         error stop 'compact projection: point-to-compact vector shape mismatch'
      end if
   end subroutine validate_point_to_compact_shapes

   subroutine validate_compact_to_point_shapes(space, product, compact_vector, point_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: compact_vector(:)
      complex(rp), intent(out) :: point_vector(:)

      call validate_common_layout(space, product)
      if (size(compact_vector) /= product%product_dimension .or. size(point_vector) /= space%ndim) then
         error stop 'compact projection: compact-to-point vector shape mismatch'
      end if
   end subroutine validate_compact_to_point_shapes

   subroutine validate_operator_shapes(space, product, operator)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: operator(:, :)

      if (space%nchannel /= 1 .or. space%nsite /= product%nsite .or. space%response_lmax /= product%response_lmax) then
         error stop 'compact operator projection: response/product layout mismatch'
      end if
      if (any(shape(operator) /= [product%product_dimension, product%product_dimension])) then
         error stop 'compact operator projection: compact operator shape mismatch'
      end if
      if (.not. allocated(space%radial_weights) .or. size(space%radial_weights) /= space%npoint) then
         error stop 'compact operator projection: LR-04 radial metric is unavailable'
      end if
   end subroutine validate_operator_shapes

   subroutine validate_gsr_inputs(space, product, susceptibility, magnetization)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: susceptibility(:, :)
      real(rp), intent(in) :: magnetization(:, :)

      call validate_common_layout(space, product)
      if (any(shape(susceptibility) /= [product%product_dimension, product%product_dimension])) then
         error stop 'evaluate_compact_goldstone_sumrule: compact susceptibility shape mismatch'
      end if
      if (any(shape(magnetization) /= [space%nsite, space%npoint])) then
         error stop 'evaluate_compact_goldstone_sumrule: magnetization shape mismatch'
      end if
      if (.not. all(ieee_is_finite(real(susceptibility, rp))) .or. .not. all(ieee_is_finite(aimag(susceptibility)))) then
         error stop 'evaluate_compact_goldstone_sumrule: compact susceptibility contains NaN or Inf'
      end if
   end subroutine validate_gsr_inputs

end module lr_compact_static_interaction_mod
