!------------------------------------------------------------------------------
! DRESP-09 transverse spatial-field insertion.
!
! This module owns the representation algebra which is already certified by
! LR-04/LR-05.  It deliberately does not manufacture a scalar-relativistic
! mixed-spin radial metric.  The direct radial/angular route is therefore
! fail-closed until the lower-component spin-angular convention is exposed by
! the accepted augmentation.
!------------------------------------------------------------------------------
module lr_dresp09_field_insertion_mod

   use precision_mod, only: rp
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex, &
      response_unflatten_superindex
   use lr_response_space_mod, only: response_space_layout, response_vector_inner_product
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis
   use lr_pauli_transition_vertex_mod, only: pauli_endpoint_state
   implicit none
   private

   character(len=*), parameter, public :: dresp09_blocked_mixed_spin_metric = &
      'BLOCKED - MIXED-SPIN RADIAL METRIC NOT CERTIFIED'

   type, public :: dresp09_direct_route_capability
      logical :: mixed_spin_radial_metric_certified = .false.
      logical :: lower_component_angular_convention_certified = .false.
   contains
      procedure :: certified => dresp09_direct_route_is_certified
   end type dresp09_direct_route_capability

   public :: dresp09_compact_field_from_raw
   public :: dresp09_raw_field_from_compact
   public :: dresp09_compact_pairing
   public :: dresp09_raw_pairing
   public :: dresp09_product_transition_tensor
   public :: dresp09_metric_adjoint_operator
   public :: dresp09_source_operator
   public :: dresp09_product_metric_adjoint_operator
   public :: dresp09_product_source_operator
   public :: dresp09_adjoint_density_coordinates
   public :: dresp09_trace_density_operator
   public :: dresp09_require_direct_route

contains

   !> Project a raw covariant field B onto the compact orthonormal product
   !> modes.  `weighted_modes` already contains W^(1/2)U, so exactly one
   !> square-root metric appears here.
   subroutine dresp09_compact_field_from_raw(space, product, raw_field, compact_field)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: raw_field(:)
      complex(rp), intent(out) :: compact_field(:)
      type(response_super_index) :: item
      integer :: flat, ir, raw_flat, mode
      real(rp) :: weight

      call validate_projection_inputs(space, product, raw_field, compact_field)
      compact_field = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, product%product_dimension
         call product%unflatten_index(flat, item%site, item%response_l, item%response_m, mode)
         item%channel = 1
         do ir = 1, space%npoint
            item%radial_point = ir
            call response_flatten_superindex(item, space%nsite, space%response_lmax, space%npoint, 1, raw_flat)
            weight = space%radial_weights(ir)
            compact_field(flat) = compact_field(flat) + &
               conjg(product%blocks(item%site, item%response_l)%weighted_modes(ir, mode))* &
               sqrt(max(0.0_rp, weight))*raw_field(raw_flat)
         end do
      end do
   end subroutine dresp09_compact_field_from_raw

   !> Reconstruct the positive-measure part of a raw field from compact modes.
   !> The origin is the null point of the LR-04 quadrature and is set to zero
   !> rather than regularized.
   subroutine dresp09_raw_field_from_compact(space, product, compact_field, raw_field)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: compact_field(:)
      complex(rp), intent(out) :: raw_field(:)
      type(response_super_index) :: item
      integer :: flat, ir, raw_flat, mode
      real(rp) :: weight

      call validate_projection_inputs(space, product, raw_field, compact_field)
      raw_field = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, product%product_dimension
         call product%unflatten_index(flat, item%site, item%response_l, item%response_m, mode)
         item%channel = 1
         do ir = 2, space%npoint
            item%radial_point = ir
            call response_flatten_superindex(item, space%nsite, space%response_lmax, space%npoint, 1, raw_flat)
            weight = space%radial_weights(ir)
            if (weight > 0.0_rp) raw_field(raw_flat) = raw_field(raw_flat) + &
               product%blocks(item%site, item%response_l)%weighted_modes(ir, mode)* &
               compact_field(flat)/sqrt(weight)
         end do
      end do
   end subroutine dresp09_raw_field_from_compact

   pure function dresp09_compact_pairing(field, density) result(value)
      complex(rp), intent(in) :: field(:), density(:)
      complex(rp) :: value

      if (size(field) /= size(density)) error stop 'DRESP-09 compact pairing: shape mismatch'
      value = sum(conjg(field)*density)
   end function dresp09_compact_pairing

   function dresp09_raw_pairing(space, field, density) result(value)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: field(:), density(:)
      complex(rp) :: value

      value = response_vector_inner_product(space, field, density)
   end function dresp09_raw_pairing

   !> Build the full compact transition tensor from the certified product
   !> component tensor.  The tensor convention is
   !> T(mu,n,m)=<n|O_mu|m>, with component order 1+p+2*q and endpoint powers
   !> E_left**p E_right**q.
   subroutine dresp09_product_transition_tensor(product, left_state, right_state, transition)
      type(lmto_product_response_basis), intent(in) :: product
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(out) :: transition(:, :, :)
      complex(rp), allocatable :: components(:, :, :, :)
      integer :: n, mu, p, q, component
      real(rp) :: left_power, right_power

      n = 2*(product%orbital_lmax + 1)**2*product%nsite
      if (any(shape(transition) /= [n, n, product%product_dimension])) then
         error stop 'DRESP-09 transition tensor: output shape mismatch'
      end if
      if (.not. allocated(left_state%coefficients) .or. .not. allocated(right_state%coefficients) .or. &
          size(left_state%coefficients) /= n .or. size(right_state%coefficients) /= n) then
         error stop 'DRESP-09 transition tensor: endpoint shape mismatch'
      end if

      call product%component_vertex_tensor(components)
      transition = cmplx(0.0_rp, 0.0_rp, rp)
      do p = 0, 1
         do q = 0, 1
            component = 1 + p + 2*q
            left_power = 1.0_rp
            right_power = 1.0_rp
            if (p == 1) left_power = left_state%energy
            if (q == 1) right_power = right_state%energy
            do mu = 1, product%product_dimension
               transition(:, :, mu) = transition(:, :, mu) + &
                  left_power*right_power*components(:, :, component, mu)
            end do
         end do
      end do
      deallocate(components)
   end subroutine dresp09_product_transition_tensor

   !> Metric-adjoint field insertion in compact orthonormal coordinates.
   !>
   !> `transition(:,:,mu)` is the measurement matrix O_mu.  The covariant
   !> compact field b represents the physical pairing b^H d, so the source
   !> operator is
   !>
   !>   V = sum_mu conjg(b_mu) O_mu^H.
   !>
   !> For a chi_plus measurement O_mu=sigma_plus, O_mu^H is the certified
   !> chi_plus source sigma_minus.  No response-space W is inserted here.
   subroutine dresp09_metric_adjoint_operator(transition, field, operator)
      complex(rp), intent(in) :: transition(:, :, :), field(:)
      complex(rp), intent(out) :: operator(:, :)
      integer :: n, mu

      n = size(transition, 1)
      if (size(transition, 2) /= n .or. size(transition, 3) /= size(field) .or. &
          any(shape(operator) /= [n, n])) then
         error stop 'DRESP-09 metric adjoint: shape mismatch'
      end if
      operator = cmplx(0.0_rp, 0.0_rp, rp)
      do mu = 1, size(field)
         operator = operator + conjg(field(mu))*transpose(conjg(transition(:, :, mu)))
      end do
   end subroutine dresp09_metric_adjoint_operator

   !> Physical source form used by the retarded chi convention.  It differs
   !> from the covariant metric-adjoint form only in whether the supplied
   !> source coefficient is conjugated.
   subroutine dresp09_source_operator(transition, field, operator)
      complex(rp), intent(in) :: transition(:, :, :), field(:)
      complex(rp), intent(out) :: operator(:, :)
      integer :: n, mu

      n = size(transition, 1)
      if (size(transition, 2) /= n .or. size(transition, 3) /= size(field) .or. &
          any(shape(operator) /= [n, n])) then
         error stop 'DRESP-09 source operator: shape mismatch'
      end if
      operator = cmplx(0.0_rp, 0.0_rp, rp)
      do mu = 1, size(field)
         operator = operator + field(mu)*transpose(conjg(transition(:, :, mu)))
      end do
   end subroutine dresp09_source_operator

   !> Product-space version of the adjoint.  The endpoint-energy powers in
   !> the certified linearized LMTO vertex are reversed by the adjoint:
   !> (p,q) -> H^q O^H H^p.
   subroutine dresp09_product_metric_adjoint_operator(component_vertices, hamiltonian, field, operator)
      complex(rp), intent(in) :: component_vertices(:, :, :, :), hamiltonian(:, :), field(:)
      complex(rp), intent(out) :: operator(:, :)

      call product_adjoint_impl(component_vertices, hamiltonian, field, operator, .true.)
   end subroutine dresp09_product_metric_adjoint_operator

   subroutine dresp09_product_source_operator(component_vertices, hamiltonian, field, operator)
      complex(rp), intent(in) :: component_vertices(:, :, :, :), hamiltonian(:, :), field(:)
      complex(rp), intent(out) :: operator(:, :)

      call product_adjoint_impl(component_vertices, hamiltonian, field, operator, .false.)
   end subroutine dresp09_product_source_operator

   !> Coordinates dual to the transition tensor for the metric-adjoint
   !> convention.  For arbitrary rho this is the exact finite-dimensional
   !> algebra used by the energy oracle; no Hermiticity assumption is made.
   subroutine dresp09_adjoint_density_coordinates(transition, density_matrix, coordinates)
      complex(rp), intent(in) :: transition(:, :, :), density_matrix(:, :)
      complex(rp), intent(out) :: coordinates(:)
      integer :: n, mu, i, j

      n = size(transition, 1)
      if (size(transition, 2) /= n .or. any(shape(density_matrix) /= [n, n]) .or. &
          size(coordinates) /= size(transition, 3)) then
         error stop 'DRESP-09 density dual: shape mismatch'
      end if
      coordinates = cmplx(0.0_rp, 0.0_rp, rp)
      do mu = 1, size(coordinates)
         do i = 1, n
            do j = 1, n
               coordinates(mu) = coordinates(mu) + conjg(transition(i, j, mu))*density_matrix(i, j)
            end do
         end do
      end do
   end subroutine dresp09_adjoint_density_coordinates

   function dresp09_trace_density_operator(density_matrix, operator) result(value)
      complex(rp), intent(in) :: density_matrix(:, :), operator(:, :)
      complex(rp) :: value
      integer :: n, i, j

      n = size(density_matrix, 1)
      if (size(density_matrix, 2) /= n .or. any(shape(operator) /= [n, n])) then
         error stop 'DRESP-09 trace pairing: shape mismatch'
      end if
      value = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, n
         do j = 1, n
            value = value + density_matrix(i, j)*operator(j, i)
         end do
      end do
   end function dresp09_trace_density_operator

   subroutine dresp09_require_direct_route(capability)
      type(dresp09_direct_route_capability), intent(in) :: capability

      if (.not. capability%certified()) then
         write (*, '(a)') trim(dresp09_blocked_mixed_spin_metric)
         error stop trim(dresp09_blocked_mixed_spin_metric)
      end if
   end subroutine dresp09_require_direct_route

   pure logical function dresp09_direct_route_is_certified(this) result(value)
      class(dresp09_direct_route_capability), intent(in) :: this

      value = this%mixed_spin_radial_metric_certified .and. &
         this%lower_component_angular_convention_certified
   end function dresp09_direct_route_is_certified

   subroutine validate_projection_inputs(space, product, raw_field, compact_field)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: raw_field(:), compact_field(:)

      if (space%nchannel /= 1 .or. size(raw_field) /= space%ndim .or. &
          size(compact_field) /= product%product_dimension .or. product%nsite /= space%nsite .or. &
          product%response_lmax /= space%response_lmax .or. product%npoint /= space%npoint) then
         error stop 'DRESP-09 compact projection: representation shape mismatch'
      end if
   end subroutine validate_projection_inputs

   subroutine product_adjoint_impl(component_vertices, hamiltonian, field, operator, conjugate_field)
      complex(rp), intent(in) :: component_vertices(:, :, :, :), hamiltonian(:, :), field(:)
      complex(rp), intent(out) :: operator(:, :)
      logical, intent(in) :: conjugate_field
      complex(rp), allocatable :: term(:, :)
      integer :: n, mu, p, q, component
      complex(rp) :: coefficient

      n = size(hamiltonian, 1)
      if (size(hamiltonian, 2) /= n .or. size(component_vertices, 1) /= n .or. &
          size(component_vertices, 2) /= n .or. size(component_vertices, 3) /= 4 .or. &
          size(component_vertices, 4) /= size(field) .or. &
          any(shape(operator) /= [n, n])) then
         error stop 'DRESP-09 product adjoint: shape mismatch'
      end if
      allocate(term(n, n))
      operator = cmplx(0.0_rp, 0.0_rp, rp)
      do mu = 1, size(field)
         do p = 0, 1
            do q = 0, 1
               component = 1 + p + 2*q
               term = transpose(conjg(component_vertices(:, :, component, mu)))
               if (q == 1) term = matmul(hamiltonian, term)
               if (p == 1) term = matmul(term, hamiltonian)
               if (conjugate_field) then
                  coefficient = conjg(field(mu))
               else
                  coefficient = field(mu)
               end if
               operator = operator + coefficient*term
            end do
         end do
      end do
      deallocate(term)
   end subroutine product_adjoint_impl

end module lr_dresp09_field_insertion_mod
