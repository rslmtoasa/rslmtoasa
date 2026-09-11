!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Discrete metric and operator algebra for the direct LR response space.
!>
!> The response coordinates are the LR-02R super-index
!> (site,L,M,radial-point,channel).  This module owns only the finite-dimensional
!> coordinate algebra.  It deliberately has no transition vertices,
!> susceptibility, XC, or Dyson physics.
!>
!> A canonical operator is the right-weighted matrix B=A W, where A is the
!> pointwise kernel and W is the diagonal physical radial quadrature metric.
!> Thus B acts on point values with ordinary matrix multiplication while vector
!> inner products retain W explicitly.  The origin has exactly zero volume
!> weight; no epsilon is used.
!------------------------------------------------------------------------------
module lr_response_space_mod

   use precision_mod, only: rp
   use response_basis_mapping_mod, only: response_super_index, response_superindex_size, &
      response_unflatten_superindex, response_simpson_weight, &
      response_volume_measure
   implicit none
   private

   type, public :: response_space_layout
      integer :: nsite = 0
      integer :: response_lmax = -1
      integer :: npoint = 0
      integer :: nchannel = 0
      integer :: ndim = 0
      integer :: active_dimension = 0
      real(rp) :: a = 0.0_rp
      real(rp) :: b = 0.0_rp
      real(rp), allocatable :: radius(:)
      real(rp), allocatable :: radial_weights(:)
      real(rp), allocatable :: metric_weights(:)
   contains
      procedure :: initialize => response_space_initialize
   end type response_space_layout

   public :: response_build_radial_metric
   public :: response_metric_weights
   public :: response_vector_inner_product
   public :: response_vector_norm
   public :: response_apply_operator
   public :: response_compose_operators
   public :: response_identity_operator
   public :: response_local_operator
   public :: response_operator_adjoint
   public :: response_operator_trace
   public :: response_rigid_vector_norm
   public :: response_rigid_vector_overlap
   public :: response_raw_to_canonical
   public :: response_canonical_to_raw

   interface response_local_operator
      module procedure response_local_operator_site_radial
      module procedure response_local_operator_site_radial_channel
      module procedure response_local_operator_complex_site_radial
      module procedure response_local_operator_complex_site_radial_channel
   end interface response_local_operator

contains

   !> Build the physical radial metric for a pointwise density on the LR-01
   !> logarithmic mesh.  The angular integral is already represented by
   !> normalized response harmonics and is therefore not included here.
   subroutine response_build_radial_metric(a, b, radius, weights)
      real(rp), intent(in) :: a, b, radius(:)
      real(rp), intent(out) :: weights(:)
      integer :: ir, npoint

      npoint = size(radius)
      if (size(weights) /= npoint .or. npoint < 3 .or. mod(npoint, 2) == 0) then
         error stop 'response_build_radial_metric: Simpson mesh must have odd size'
      end if
      if (a <= 0.0_rp .or. b <= 0.0_rp) then
         error stop 'response_build_radial_metric: logarithmic mesh parameters must be positive'
      end if
      if (radius(1) /= 0.0_rp) then
         error stop 'response_build_radial_metric: production response mesh must include r=0'
      end if
      do ir = 2, npoint
         if (radius(ir) <= radius(ir - 1)) then
            error stop 'response_build_radial_metric: radius must be strictly increasing'
         end if
      end do

      do ir = 1, npoint
         weights(ir) = response_simpson_weight(ir, npoint)*response_volume_measure(a, b, radius(ir))
      end do
   end subroutine response_build_radial_metric

   subroutine response_space_initialize(this, nsite, response_lmax, radius, a, b, nchannel)
      class(response_space_layout), intent(out) :: this
      integer, intent(in) :: nsite, response_lmax, nchannel
      real(rp), intent(in) :: radius(:), a, b
      integer :: flat, radial_point

      if (nsite < 1 .or. response_lmax < 0 .or. nchannel < 1) then
         error stop 'response_space_layout%initialize: invalid response dimensions'
      end if
      if (size(radius) < 3 .or. mod(size(radius), 2) == 0) then
         error stop 'response_space_layout%initialize: invalid Simpson mesh size'
      end if

      this%nsite = nsite
      this%response_lmax = response_lmax
      this%npoint = size(radius)
      this%nchannel = nchannel
      this%ndim = response_superindex_size(nsite, response_lmax, size(radius), nchannel)
      this%active_dimension = response_superindex_size(nsite, response_lmax, size(radius) - 1, nchannel)
      this%a = a
      this%b = b
      allocate(this%radius(size(radius)), this%radial_weights(size(radius)), this%metric_weights(this%ndim))
      this%radius = radius
      call response_build_radial_metric(a, b, radius, this%radial_weights)

      ! The LR-02R order is site, harmonic, radial point, channel.  Repeating
      ! the radial metric in this order makes every contraction use the same
      ! super-index without introducing a second indexing implementation.
      do flat = 1, this%ndim
         radial_point = mod((flat - 1)/nchannel, this%npoint) + 1
         this%metric_weights(flat) = this%radial_weights(radial_point)
      end do
   end subroutine response_space_initialize

   subroutine response_metric_weights(space, weights)
      type(response_space_layout), intent(in) :: space
      real(rp), intent(out) :: weights(:)

      call require_layout(space)
      if (size(weights) /= space%ndim) error stop 'response_metric_weights: vector shape mismatch'
      weights = space%metric_weights
   end subroutine response_metric_weights

   function response_vector_inner_product(space, left, right) result(value)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: left(:), right(:)
      complex(rp) :: value

      call require_layout(space)
      call require_vector(space, left, 'response_vector_inner_product: left vector shape mismatch')
      call require_vector(space, right, 'response_vector_inner_product: right vector shape mismatch')
      value = sum(conjg(left)*space%metric_weights*right)
   end function response_vector_inner_product

   function response_vector_norm(space, vector) result(value)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: vector(:)
      real(rp) :: value
      complex(rp) :: inner

      inner = response_vector_inner_product(space, vector, vector)
      if (real(inner, rp) < -100.0_rp*epsilon(1.0_rp)*max(1.0_rp, abs(inner))) then
         error stop 'response_vector_norm: negative metric norm'
      end if
      value = sqrt(max(0.0_rp, real(inner, rp)))
   end function response_vector_norm

   subroutine response_apply_operator(space, operator, field, result)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: operator(:, :), field(:)
      complex(rp), intent(out) :: result(:)

      call require_operator(space, operator, 'response_apply_operator: operator shape mismatch')
      call require_vector(space, field, 'response_apply_operator: field shape mismatch')
      call require_vector(space, result, 'response_apply_operator: result shape mismatch')
      result = matmul(operator, field)
   end subroutine response_apply_operator

   !> Composition follows directly from B=A W.  If C=A W and D=D_raw W,
   !> then the composite canonical matrix is C D; no additional local metric
   !> factor is allowed here.
   subroutine response_compose_operators(space, first, second, composed)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: first(:, :), second(:, :)
      complex(rp), intent(out) :: composed(:, :)

      call require_operator(space, first, 'response_compose_operators: first shape mismatch')
      call require_operator(space, second, 'response_compose_operators: second shape mismatch')
      call require_operator(space, composed, 'response_compose_operators: result shape mismatch')
      composed = matmul(first, second)
   end subroutine response_compose_operators

   subroutine response_identity_operator(space, identity)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(out) :: identity(:, :)
      integer :: flat

      call require_operator(space, identity, 'response_identity_operator: result shape mismatch')
      identity = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         identity(flat, flat) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
   end subroutine response_identity_operator

   !> Spherical local radial scalars are diagonal in site, L, M, radial point,
   !> and channel.  A rank-2 input is shared by all channels.
   subroutine response_local_operator_site_radial(space, values, operator)
      type(response_space_layout), intent(in) :: space
      real(rp), intent(in) :: values(:, :)
      complex(rp), intent(out) :: operator(:, :)
      type(response_super_index) :: item
      integer :: flat

      call require_operator(space, operator, 'response_local_operator: result shape mismatch')
      if (size(values, 1) /= space%nsite .or. size(values, 2) /= space%npoint) then
         error stop 'response_local_operator: site/radial scalar shape mismatch'
      end if
      operator = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, &
            space%npoint, space%nchannel, item)
         operator(flat, flat) = cmplx(values(item%site, item%radial_point), 0.0_rp, rp)
      end do
   end subroutine response_local_operator_site_radial

   subroutine response_local_operator_site_radial_channel(space, values, operator)
      type(response_space_layout), intent(in) :: space
      real(rp), intent(in) :: values(:, :, :)
      complex(rp), intent(out) :: operator(:, :)
      type(response_super_index) :: item
      integer :: flat

      call require_operator(space, operator, 'response_local_operator: result shape mismatch')
      if (size(values, 1) /= space%nsite .or. size(values, 2) /= space%npoint .or. &
          size(values, 3) /= space%nchannel) then
         error stop 'response_local_operator: site/radial/channel scalar shape mismatch'
      end if
      operator = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, &
            space%npoint, space%nchannel, item)
         operator(flat, flat) = cmplx(values(item%site, item%radial_point, item%channel), 0.0_rp, rp)
      end do
   end subroutine response_local_operator_site_radial_channel

   !> Complex-valued spherical local scalar.  This has exactly the same
   !> pointwise/canonical action as the real overload and is needed by static
   !> response routes evaluated with a finite retarded broadening.
   subroutine response_local_operator_complex_site_radial(space, values, operator)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: values(:, :)
      complex(rp), intent(out) :: operator(:, :)
      type(response_super_index) :: item
      integer :: flat

      call require_operator(space, operator, 'response_local_operator: result shape mismatch')
      if (size(values, 1) /= space%nsite .or. size(values, 2) /= space%npoint) then
         error stop 'response_local_operator: site/radial scalar shape mismatch'
      end if
      operator = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, &
            space%npoint, space%nchannel, item)
         operator(flat, flat) = values(item%site, item%radial_point)
      end do
   end subroutine response_local_operator_complex_site_radial

   subroutine response_local_operator_complex_site_radial_channel(space, values, operator)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: values(:, :, :)
      complex(rp), intent(out) :: operator(:, :)
      type(response_super_index) :: item
      integer :: flat

      call require_operator(space, operator, 'response_local_operator: result shape mismatch')
      if (size(values, 1) /= space%nsite .or. size(values, 2) /= space%npoint .or. &
          size(values, 3) /= space%nchannel) then
         error stop 'response_local_operator: site/radial/channel scalar shape mismatch'
      end if
      operator = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, &
            space%npoint, space%nchannel, item)
         operator(flat, flat) = values(item%site, item%radial_point, item%channel)
      end do
   end subroutine response_local_operator_complex_site_radial_channel

   !> Metric adjoint is defined by <x,B y>_W=<B^dagger x,y>_W.
   !> The origin entries form a null subspace of the quadrature metric.  The
   !> canonical extension is the ordinary conjugate transpose within that
   !> null subspace, with no active/null mixing.  A finite quadrature operator
   !> must not map a null input at the origin into an active output.
   subroutine response_operator_adjoint(space, operator, adjoint)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: operator(:, :)
      complex(rp), intent(out) :: adjoint(:, :)
      integer :: i, j

      call require_operator(space, operator, 'response_operator_adjoint: operator shape mismatch')
      call require_operator(space, adjoint, 'response_operator_adjoint: result shape mismatch')
      do j = 1, space%ndim
         if (space%metric_weights(j) == 0.0_rp) then
            do i = 1, space%ndim
               if (space%metric_weights(i) > 0.0_rp .and. operator(i, j) /= cmplx(0.0_rp, 0.0_rp, rp)) then
                  error stop 'response_operator_adjoint: null-origin input is not quadrature compatible'
               end if
            end do
         end if
      end do

      adjoint = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, space%ndim
         do j = 1, space%ndim
            if (space%metric_weights(i) > 0.0_rp .and. space%metric_weights(j) > 0.0_rp) then
               adjoint(i, j) = conjg(operator(j, i))*space%metric_weights(j)/space%metric_weights(i)
            elseif (space%metric_weights(i) == 0.0_rp .and. space%metric_weights(j) == 0.0_rp) then
               adjoint(i, j) = conjg(operator(j, i))
            end if
         end do
      end do
   end subroutine response_operator_adjoint

   !> Trace is the trace on the positive-measure response space.  Null origin
   !> point values are intentionally excluded.
   function response_operator_trace(space, operator) result(value)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: operator(:, :)
      complex(rp) :: value
      integer :: flat

      call require_operator(space, operator, 'response_operator_trace: operator shape mismatch')
      value = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         if (space%metric_weights(flat) > 0.0_rp) value = value + operator(flat, flat)
      end do
   end function response_operator_trace

   function response_rigid_vector_norm(space, rigid_vector) result(value)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: rigid_vector(:)
      real(rp) :: value

      value = response_vector_norm(space, rigid_vector)
   end function response_rigid_vector_norm

   function response_rigid_vector_overlap(space, vector, rigid_vector) result(value)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: vector(:), rigid_vector(:)
      complex(rp) :: value

      value = response_vector_inner_product(space, vector, rigid_vector)
   end function response_rigid_vector_overlap

   !> Convert a raw pointwise kernel A to canonical B=A W.  The origin column
   !> is exactly zero because its quadrature measure is exactly zero.
   subroutine response_raw_to_canonical(space, raw, canonical)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: raw(:, :)
      complex(rp), intent(out) :: canonical(:, :)
      integer :: column

      call require_operator(space, raw, 'response_raw_to_canonical: raw shape mismatch')
      call require_operator(space, canonical, 'response_raw_to_canonical: result shape mismatch')
      do column = 1, space%ndim
         canonical(:, column) = raw(:, column)*space%metric_weights(column)
      end do
   end subroutine response_raw_to_canonical

   !> Convert canonical B back to the finite raw kernel on its active columns.
   !> A nonzero origin column has no finite pointwise-kernel representation and
   !> is rejected rather than hidden with a regularizer.
   subroutine response_canonical_to_raw(space, canonical, raw)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: canonical(:, :)
      complex(rp), intent(out) :: raw(:, :)
      integer :: column

      call require_operator(space, canonical, 'response_canonical_to_raw: canonical shape mismatch')
      call require_operator(space, raw, 'response_canonical_to_raw: result shape mismatch')
      raw = cmplx(0.0_rp, 0.0_rp, rp)
      do column = 1, space%ndim
         if (space%metric_weights(column) > 0.0_rp) then
            raw(:, column) = canonical(:, column)/space%metric_weights(column)
         elseif (any(canonical(:, column) /= cmplx(0.0_rp, 0.0_rp, rp))) then
            error stop 'response_canonical_to_raw: nonzero null-origin column is not representable'
         end if
      end do
   end subroutine response_canonical_to_raw

   subroutine require_layout(space)
      type(response_space_layout), intent(in) :: space

      if (space%ndim < 1 .or. .not. allocated(space%metric_weights) .or. &
          size(space%metric_weights) /= space%ndim) then
         error stop 'LR response space: layout is not initialized'
      end if
   end subroutine require_layout

   subroutine require_vector(space, vector, message)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: vector(:)
      character(len=*), intent(in) :: message

      call require_layout(space)
      if (size(vector) /= space%ndim) error stop message
   end subroutine require_vector

   subroutine require_operator(space, operator, message)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: operator(:, :)
      character(len=*), intent(in) :: message

      call require_layout(space)
      if (any(shape(operator) /= [space%ndim, space%ndim])) error stop message
   end subroutine require_operator

end module lr_response_space_mod
