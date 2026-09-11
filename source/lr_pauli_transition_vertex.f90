!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Production Pauli/no-SOC transition-vector evaluator.
!>
!> This module owns the one-electron response vertex only.  It consumes two
!> endpoint eigenstates and the production LMTO large-component augmentation,
!> and emits pointwise coefficients in the LR-04 response space.  It does not
!> know about occupations, frequency denominators, susceptibility, XC, or
!> Dyson equations.
!------------------------------------------------------------------------------
module lr_pauli_transition_vertex_mod

   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use response_angular_basis_mod, only: response_angular_pi, response_gaunt
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex, &
      response_unflatten_superindex
   use lr_response_space_mod, only: response_space_layout
   implicit none
   private

   !> Capability tuple for the initial LR-05 production vertex.
   type, public :: pauli_vertex_capabilities
      character(len=32) :: reciprocal_mode = 'ham_only'
      character(len=16) :: hamiltonian_order = 'second'
      logical :: orthogonal = .true.
      logical :: collinear = .true.
      logical :: has_soc = .false.
      logical :: has_extra_operator = .false.
   end type pauli_vertex_capabilities

   !> One reciprocal endpoint eigenstate in production site-blocked ordering.
   !>
   !> The coefficient vector is packed as
   !> `(site, spin-up orbitals, spin-down orbitals)`, with the orbital order
   !> `(s),(p,-1:1),(d,-2:2),...` supplied by `basis_mod`/LMTO.
   type, public :: pauli_endpoint_state
      real(rp) :: energy = 0.0_rp
      complex(rp), allocatable :: coefficients(:)
   contains
      procedure :: initialize => pauli_endpoint_state_initialize
   end type pauli_endpoint_state

   public :: pauli_charge_matrix
   public :: pauli_sigma_x_matrix
   public :: pauli_sigma_y_matrix
   public :: pauli_sigma_z_matrix
   public :: pauli_sigma_plus_matrix
   public :: pauli_sigma_minus_matrix
   public :: evaluate_pauli_transition_vertex

   interface evaluate_pauli_transition_vertex
      module procedure evaluate_pauli_transition_vertex_one
      module procedure evaluate_pauli_transition_vertex_channels
   end interface evaluate_pauli_transition_vertex

contains

   subroutine pauli_endpoint_state_initialize(this, energy, coefficients)
      class(pauli_endpoint_state), intent(out) :: this
      real(rp), intent(in) :: energy
      complex(rp), intent(in) :: coefficients(:)

      this%energy = energy
      allocate(this%coefficients(size(coefficients)))
      this%coefficients = coefficients
   end subroutine pauli_endpoint_state_initialize

   pure function pauli_charge_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      matrix(1, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      matrix(2, 2) = cmplx(1.0_rp, 0.0_rp, rp)
   end function pauli_charge_matrix

   pure function pauli_sigma_x_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      matrix(1, 2) = cmplx(1.0_rp, 0.0_rp, rp)
      matrix(2, 1) = cmplx(1.0_rp, 0.0_rp, rp)
   end function pauli_sigma_x_matrix

   pure function pauli_sigma_y_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      matrix(1, 2) = cmplx(0.0_rp, -1.0_rp, rp)
      matrix(2, 1) = cmplx(0.0_rp, 1.0_rp, rp)
   end function pauli_sigma_y_matrix

   pure function pauli_sigma_z_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      matrix(1, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      matrix(2, 2) = cmplx(-1.0_rp, 0.0_rp, rp)
   end function pauli_sigma_z_matrix

   !> The typed circular convention is encoded by matrices, not an integer
   !> selector: sigma_plus=(sigma_x+i sigma_y)/2 and sigma_minus=(sigma_x-i sigma_y)/2.
   pure function pauli_sigma_plus_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      matrix(1, 2) = cmplx(1.0_rp, 0.0_rp, rp)
   end function pauli_sigma_plus_matrix

   pure function pauli_sigma_minus_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      matrix(2, 1) = cmplx(1.0_rp, 0.0_rp, rp)
   end function pauli_sigma_minus_matrix

   !> Evaluate one operator channel.  The response layout must have one
   !> channel; use the rank-3 overload for a vector containing several explicit
   !> Pauli operators.
   subroutine evaluate_pauli_transition_vertex_one(space, radial_bases, left_state, right_state, operator_matrix, &
                                                   capabilities, transition_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(in) :: operator_matrix(:, :)
      type(pauli_vertex_capabilities), intent(in) :: capabilities
      complex(rp), intent(out) :: transition_vector(:)
      complex(rp), allocatable :: operator_channels(:, :, :)

      if (size(operator_matrix, 1) /= 2 .or. size(operator_matrix, 2) /= 2) then
         error stop 'evaluate_pauli_transition_vertex: Pauli operator must be 2x2'
      end if
      if (space%nchannel /= 1) then
         error stop 'evaluate_pauli_transition_vertex: one operator requires nchannel=1'
      end if

      allocate(operator_channels(2, 2, 1))
      operator_channels(:, :, 1) = operator_matrix
      call evaluate_pauli_transition_vertex_channels(space, radial_bases, left_state, right_state, operator_channels, &
                                                      capabilities, transition_vector)
   end subroutine evaluate_pauli_transition_vertex_one

   !> Evaluate one or more explicit 2x2 Pauli-space operators.  The third
   !> operator dimension is the LR-04 response channel, so the complete
   !> `(site,L,M,radial,mu)` vector is emitted in one call.
   subroutine evaluate_pauli_transition_vertex_channels(space, radial_bases, left_state, right_state, operator_matrices, &
                                                        capabilities, transition_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(in) :: operator_matrices(:, :, :)
      type(pauli_vertex_capabilities), intent(in) :: capabilities
      complex(rp), intent(out) :: transition_vector(:)

      integer :: nsite, norb, npoint, nresponse, nchannel
      integer :: isite, ir, response_l, response_m, channel, iorb, jorb, ispin, jspin
      integer :: orbital_l, orbital_lp, left_offset, right_offset, flat
      type(response_super_index) :: item
      complex(rp), allocatable :: left_radial(:, :, :, :), right_radial(:, :, :, :)
      complex(rp) :: value, coefficient_product
      real(rp) :: radial_product

      call validate_vertex_inputs(space, radial_bases, left_state, right_state, operator_matrices, capabilities, transition_vector)

      nsite = space%nsite
      norb = (radial_bases(1)%lmax + 1)**2
      npoint = space%npoint
      nresponse = space%response_lmax
      nchannel = space%nchannel
      allocate(left_radial(nsite, npoint, norb, 2), right_radial(nsite, npoint, norb, 2))

      ! Reuse the production large-component augmentation directly.  The
      ! scalar-relativistic small component and GFAC are intentionally absent
      ! from this Pauli/no-SOC primitive.
      do isite = 1, nsite
         left_offset = (isite - 1)*2*norb
         right_offset = left_offset
         do ispin = 1, 2
            do iorb = 1, norb
               orbital_l = lmto_orbital_l(iorb)
               left_radial(isite, :, iorb, ispin) = cmplx(0.0_rp, 0.0_rp, rp)
               right_radial(isite, :, iorb, ispin) = cmplx(0.0_rp, 0.0_rp, rp)
               do ir = 1, npoint
                  left_radial(isite, ir, iorb, ispin) = cmplx( &
                     radial_bases(isite)%phi_large(ir, orbital_l + 1, ispin) + &
                     (left_state%energy - radial_bases(isite)%enu_work(orbital_l + 1, ispin))* &
                     radial_bases(isite)%phidot_large(ir, orbital_l + 1, ispin), 0.0_rp, rp)
                  right_radial(isite, ir, iorb, ispin) = cmplx( &
                     radial_bases(isite)%phi_large(ir, orbital_l + 1, ispin) + &
                     (right_state%energy - radial_bases(isite)%enu_work(orbital_l + 1, ispin))* &
                     radial_bases(isite)%phidot_large(ir, orbital_l + 1, ispin), 0.0_rp, rp)
               end do
            end do
         end do
      end do

      transition_vector = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, nsite
         left_offset = (isite - 1)*2*norb
         right_offset = left_offset
         do response_l = 0, nresponse
            do response_m = -response_l, response_l
               do ir = 1, npoint
                  do channel = 1, nchannel
                     value = cmplx(0.0_rp, 0.0_rp, rp)
                     do iorb = 1, norb
                        orbital_l = lmto_orbital_l(iorb)
                       do jorb = 1, norb
                           orbital_lp = lmto_orbital_l(jorb)
                           do ispin = 1, 2
                              do jspin = 1, 2
                                 radial_product = radial_product_at_point(left_radial(isite, :, iorb, ispin), &
                                    right_radial(isite, :, jorb, jspin), radial_bases(isite)%rofi, orbital_l, orbital_lp, ir)
                                 if (radial_product /= 0.0_rp) then
                                    coefficient_product = conjg(left_state%coefficients(left_offset + &
                                       (ispin - 1)*norb + iorb))*right_state%coefficients(right_offset + &
                                       (jspin - 1)*norb + jorb)
                                    value = value + coefficient_product*operator_matrices(ispin, jspin, channel)* &
                                       radial_product*response_gaunt(orbital_l, orbital_m(iorb, orbital_l), &
                                       orbital_lp, orbital_m(jorb, orbital_lp), response_l, response_m)
                                 end if
                              end do
                           end do
                        end do
                     end do
                     item = response_super_index(isite, response_l, response_m, ir, channel)
                     call response_flatten_superindex(item, space%nsite, space%response_lmax, space%npoint, &
                                                      space%nchannel, flat)
                     transition_vector(flat) = value
                  end do
               end do
            end do
         end do
      end do

      deallocate(left_radial, right_radial)
   end subroutine evaluate_pauli_transition_vertex_channels

   subroutine validate_vertex_inputs(space, radial_bases, left_state, right_state, operator_matrices, capabilities, transition_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(in) :: operator_matrices(:, :, :)
      type(pauli_vertex_capabilities), intent(in) :: capabilities
      complex(rp), intent(out) :: transition_vector(:)

      integer :: isite, lmax, norb
      real(rp), parameter :: mesh_tolerance = 2.0e-13_rp

      ! This call is deliberately inside the public evaluator validation
      ! path.  Unsupported representations cannot be enabled by a caller-side
      ! convention or by accidentally supplying compatible-sized arrays.
      if (size(radial_bases) /= space%nsite) then
         error stop 'evaluate_pauli_transition_vertex: one radial basis is required per response site'
      end if
      do isite = 1, size(radial_bases)
         call radial_bases(isite)%require_supported(capabilities%reciprocal_mode, capabilities%hamiltonian_order, &
            capabilities%orthogonal, capabilities%collinear, capabilities%has_soc, capabilities%has_extra_operator)
      end do

      if (size(operator_matrices, 1) /= 2 .or. size(operator_matrices, 2) /= 2 .or. &
          size(operator_matrices, 3) /= space%nchannel) then
         error stop 'evaluate_pauli_transition_vertex: operator/channel shape mismatch'
      end if
      if (size(transition_vector) /= space%ndim) then
         error stop 'evaluate_pauli_transition_vertex: response vector shape mismatch'
      end if

      lmax = radial_bases(1)%lmax
      norb = (lmax + 1)**2
      if (space%response_lmax < 2*lmax) then
         error stop 'evaluate_pauli_transition_vertex: response angular product space is incomplete'
      end if
      if (.not. allocated(space%radius) .or. size(space%radius) /= radial_bases(1)%npoint .or. &
          space%npoint /= radial_bases(1)%npoint) then
         error stop 'evaluate_pauli_transition_vertex: response/radial mesh shape mismatch'
      end if
      do isite = 1, size(radial_bases)
         if (radial_bases(isite)%lmax /= lmax .or. radial_bases(isite)%nspin /= 2 .or. &
             radial_bases(isite)%npoint /= space%npoint .or. .not. allocated(radial_bases(isite)%rofi) .or. &
             .not. allocated(radial_bases(isite)%channel_present) .or. &
             .not. all(radial_bases(isite)%channel_present)) then
            error stop 'evaluate_pauli_transition_vertex: incomplete or inconsistent Pauli radial basis'
         end if
         if (abs(radial_bases(isite)%mesh_a - space%a) > mesh_tolerance*max(1.0_rp, abs(space%a)) .or. &
             abs(radial_bases(isite)%mesh_b - space%b) > mesh_tolerance*max(1.0_rp, abs(space%b)) .or. &
             maxval(abs(radial_bases(isite)%rofi - space%radius)) > mesh_tolerance* &
             max(1.0_rp, maxval(abs(space%radius)))) then
            error stop 'evaluate_pauli_transition_vertex: radial mesh provenance mismatch'
         end if
      end do

      if (.not. allocated(left_state%coefficients) .or. .not. allocated(right_state%coefficients) .or. &
          size(left_state%coefficients) /= 2*norb*space%nsite .or. &
          size(right_state%coefficients) /= 2*norb*space%nsite) then
         error stop 'evaluate_pauli_transition_vertex: endpoint coefficient shape mismatch'
      end if
   end subroutine validate_vertex_inputs

   pure integer function orbital_m(iorb, l) result(m)
      integer, intent(in) :: iorb, l
      m = iorb - l*l - l - 1
   end function orbital_m

   function radial_product_at_point(left_values, right_values, radius, l, lp, ir) result(value)
      complex(rp), intent(in) :: left_values(:), right_values(:)
      real(rp), intent(in) :: radius(:)
      integer, intent(in) :: l, lp, ir
      real(rp) :: value
      real(rp) :: first, second, rfirst, rsecond

      if (ir /= 1) then
         if (radius(ir) <= tiny(1.0_rp)) error stop 'evaluate_pauli_transition_vertex: invalid positive radial point'
         value = real(left_values(ir)*right_values(ir), rp)/(radius(ir)**2)
         return
      end if

      ! Regular LMTO radial numerators obey U_l=O(r**(l+1)).  Therefore every
      ! origin limit with l plus l prime greater than zero is zero.  The only
      ! nonzero possible limit is s-s; extrapolate its regular U_0 U_0/r^2 value
      ! from the first two
      ! positive mesh points, matching the existing origin-density contract.
      if (l /= 0 .or. lp /= 0) then
         value = 0.0_rp
         return
      end if
      rfirst = radius(2)
      rsecond = radius(3)
      first = real(left_values(2)*right_values(2), rp)/(rfirst**2)
      second = real(left_values(3)*right_values(3), rp)/(rsecond**2)
      value = (first*rsecond**2 - second*rfirst**2)/(rsecond**2 - rfirst**2)
   end function radial_product_at_point

end module lr_pauli_transition_vertex_mod
