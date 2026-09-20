!------------------------------------------------------------------------------
! DRESP-09S scalar-relativistic L=0 augmentation contract.
!
! The radial arrays are the unchanged RSEQSR/PHDFSR outputs.  This module
! supplies only the recovered bilinear and the independent L=0 source/density
! contractions; it does not call the DRESP-08 or DRESP-09R oracle routes.
!------------------------------------------------------------------------------
module lr_dresp09s_scalar_relativistic_mod

   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use response_angular_basis_mod, only: response_gaunt
   use response_basis_mapping_mod, only: response_super_index, response_unflatten_superindex
   use lr_response_space_mod, only: response_space_layout
   implicit none
   private

   ! The angular average is over the two scalar-relativistic spin-orbit
   ! partners.  For a transverse Pauli operator,
   ! <sigma_r sigma_+ sigma_r>_Omega = -sigma_+/3.
   real(rp), parameter, public :: sr_l0_lower_transverse_factor = -1.0_rp/3.0_rp
   real(rp), parameter, public :: sr_l0_lower_spin_factor = -1.0_rp/3.0_rp

   public :: sr_endpoint_component
   public :: sr_second_order_branch_point
   public :: sr_physical_spin_point
   public :: sr_same_spin_numerator
   public :: sr_mixed_transverse_point
   public :: sr_l0_source_components
   public :: sr_l0_source
   public :: sr_l0_source_with_flags
   public :: sr_l0_density_from_matrix
   public :: sr_l0_density_from_matrix_hamiltonian
   public :: sr_l0_density_tangent_from_hamiltonian
   public :: sr_l0_density_from_endpoint_branches
   public :: sr_l0_density_from_second_order_endpoint_branches
   public :: sr_l0_observable_matrix

contains

   !> Linearized endpoint component.  Power 0 is phi-E_nu_work*phidot;
   !> power 1 is phidot.  Component 1 is the K-H upper numerator G_l and
   !> component 2 is the subsidiary/minor radial numerator Phi_l.
   pure real(rp) function sr_endpoint_component(radial, ir, l, ispin, power, component) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, ispin, power, component

      if (component == 1) then
         if (power == 0) then
            value = radial%phi_large(ir, l + 1, ispin) - radial%enu_work(l + 1, ispin)* &
               radial%phidot_large(ir, l + 1, ispin)
         else
            value = radial%phidot_large(ir, l + 1, ispin)
         end if
      else if (component == 2) then
         if (power == 0) then
            value = radial%phi_small(ir, l + 1, ispin) - radial%enu_work(l + 1, ispin)* &
               radial%phidot_small(ir, l + 1, ispin)
         else
            value = radial%phidot_small(ir, l + 1, ispin)
         end if
      else
         error stop 'DRESP-09S: endpoint component must be 1 or 2'
      end if
   end function sr_endpoint_component

   !> Full second-order radial endpoint coefficient in the absolute-energy
   !> polynomial.  Branches are ordered 00,10,01,11,20,02.  `lower_factor`
   !> selects the observable: +1 is the scalar-relativistic probability
   !> metric and -1/3 is the physical Pauli-spin operator.
   pure real(rp) function sr_second_order_branch_point(radial, ir, l, spin_left, spin_right, branch, lower_factor) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, spin_left, spin_right, branch
      real(rp), intent(in) :: lower_factor
      real(rp) :: enu_left, enu_right
      real(rp) :: raw00, raw10, raw01, raw11, raw20, raw02

      if (branch < 1 .or. branch > 6) error stop 'DRESP-09X: invalid second-order radial branch'
      if (ir == 1 .or. radial%rofi(ir) <= tiny(1.0_rp)) then
         value = 0.0_rp
         return
      end if
      enu_left = radial%enu_work(l + 1, spin_left)
      enu_right = radial%enu_work(l + 1, spin_right)
      raw00 = sr_raw_bilinear(radial, ir, l, spin_left, spin_right, 0, 0, lower_factor)
      raw10 = sr_raw_bilinear(radial, ir, l, spin_left, spin_right, 1, 0, lower_factor)
      raw01 = sr_raw_bilinear(radial, ir, l, spin_left, spin_right, 0, 1, lower_factor)
      raw11 = sr_raw_bilinear(radial, ir, l, spin_left, spin_right, 1, 1, lower_factor)
      raw20 = 0.5_rp*sr_raw_bilinear(radial, ir, l, spin_left, spin_right, 2, 0, lower_factor)
      raw02 = 0.5_rp*sr_raw_bilinear(radial, ir, l, spin_left, spin_right, 0, 2, lower_factor)
      select case (branch)
      case (1)
         value = raw00 - enu_left*raw10 - enu_right*raw01 + enu_left*enu_right*raw11 + &
            0.5_rp*enu_left**2*sr_raw_bilinear(radial, ir, l, spin_left, spin_right, 2, 0, lower_factor) + &
            0.5_rp*enu_right**2*sr_raw_bilinear(radial, ir, l, spin_left, spin_right, 0, 2, lower_factor)
      case (2)
         value = raw10 - enu_right*raw11 - enu_left*sr_raw_bilinear(radial, ir, l, spin_left, spin_right, 2, 0, lower_factor)
      case (3)
         value = raw01 - enu_left*raw11 - enu_right*sr_raw_bilinear(radial, ir, l, spin_left, spin_right, 0, 2, lower_factor)
      case (4)
         value = raw11
      case (5)
         value = raw20
      case (6)
         value = raw02
      end select
      value = value/radial%rofi(ir)**2
   end function sr_second_order_branch_point

   !> Same-spin physical longitudinal Pauli-spin radial observable.  The
   !> sign is supplied by the caller (up=+1, down=-1), while the lower
   !> component always carries the angular Pauli factor -1/3.
   pure real(rp) function sr_physical_spin_point(radial, ir, l, ispin, branch) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, ispin, branch
      value = sr_second_order_branch_point(radial, ir, l, ispin, ispin, branch, sr_l0_lower_spin_factor)
   end function sr_physical_spin_point

   pure real(rp) function sr_raw_bilinear(radial, ir, l, spin_left, spin_right, left_power, right_power, lower_factor) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, spin_left, spin_right, left_power, right_power
      real(rp), intent(in) :: lower_factor
      real(rp) :: left_large, right_large, left_small, right_small, r, angular

      left_large = sr_raw_component(radial, ir, l, spin_left, left_power, 1)
      right_large = sr_raw_component(radial, ir, l, spin_right, right_power, 1)
      left_small = sr_raw_component(radial, ir, l, spin_left, left_power, 2)
      right_small = sr_raw_component(radial, ir, l, spin_right, right_power, 2)
      r = radial%rofi(ir)
      angular = real(l*(l + 1), rp)*left_large*right_large/ &
         (radial%tmc(ir, l + 1, spin_left)*radial%tmc(ir, l + 1, spin_right)*r**2)
      value = left_large*right_large + lower_factor*(left_small*right_small + angular)
   end function sr_raw_bilinear

   pure real(rp) function sr_raw_component(radial, ir, l, ispin, power, component) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, ispin, power, component
      select case (power)
      case (0)
         if (component == 1) then
            value = radial%phi_large(ir, l + 1, ispin)
         else
            value = radial%phi_small(ir, l + 1, ispin)
         end if
      case (1)
         if (component == 1) then
            value = radial%phidot_large(ir, l + 1, ispin)
         else
            value = radial%phidot_small(ir, l + 1, ispin)
         end if
      case (2)
         if (component == 1) then
            value = radial%phiddot_large(ir, l + 1, ispin)
         else
            value = radial%phiddot_small(ir, l + 1, ispin)
         end if
      case default
         error stop 'DRESP-09X: raw radial endpoint power must be 0, 1, or 2'
      end select
   end function sr_raw_component

   !> Numerator of the same-spin scalar-relativistic metric.  For unlike
   !> endpoints the product TMC_left*TMC_right is used only by the transverse
   !> lower angular term below; same-channel use reproduces GFAC exactly.
   pure real(rp) function sr_same_spin_numerator(radial, ir, l, ispin, p, q) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, ispin, p, q
      real(rp) :: left_large, right_large, left_small, right_small

      left_large = sr_endpoint_component(radial, ir, l, ispin, p, 1)
      right_large = sr_endpoint_component(radial, ir, l, ispin, q, 1)
      left_small = sr_endpoint_component(radial, ir, l, ispin, p, 2)
      right_small = sr_endpoint_component(radial, ir, l, ispin, q, 2)
      value = radial%gfac(ir, l + 1, ispin)*left_large*right_large + left_small*right_small
   end function sr_same_spin_numerator

   !> Physical pointwise matrix element for an L=0 transverse Pauli field
   !> between scalar channels.  The upper angular matrix element is +1.
   !> The lower object is sigma_r times the K-H bracket; its angular average
   !> gives -1/3 for sigma_+ or sigma_-.  The result includes the physical
   !> 1/r^2 wavefunction factor, so multiplying by response_space%radial_weights
   !> performs the radial volume integral.
   pure real(rp) function sr_mixed_transverse_point(radial, ir, l, spin_left, p, spin_right, q) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, spin_left, p, spin_right, q
      value = sr_mixed_transverse_point_with_flags(radial, ir, l, spin_left, p, spin_right, q, .true., .true.)
   end function sr_mixed_transverse_point

   pure real(rp) function sr_mixed_transverse_point_with_flags(radial, ir, l, spin_left, p, spin_right, q, &
                                                               include_lower_radial, include_lower_angular) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, spin_left, p, spin_right, q
      logical, intent(in) :: include_lower_radial, include_lower_angular
      real(rp) :: left_large, right_large, left_small, right_small, r
      real(rp) :: lower_angular

      if (ir == 1 .or. radial%rofi(ir) <= tiny(1.0_rp)) then
         value = 0.0_rp
         return
      end if
      r = radial%rofi(ir)
      left_large = sr_endpoint_component(radial, ir, l, spin_left, p, 1)
      right_large = sr_endpoint_component(radial, ir, l, spin_right, q, 1)
      left_small = sr_endpoint_component(radial, ir, l, spin_left, p, 2)
      right_small = sr_endpoint_component(radial, ir, l, spin_right, q, 2)
      lower_angular = 0.0_rp
      if (include_lower_radial) lower_angular = left_small*right_small
      if (include_lower_angular) lower_angular = lower_angular + real(l*(l + 1), rp)*left_large*right_large/ &
         (radial%tmc(ir, l + 1, spin_left)*radial%tmc(ir, l + 1, spin_right)*r**2)
      value = (left_large*right_large + sr_l0_lower_transverse_factor*lower_angular)/r**2
   end function sr_mixed_transverse_point_with_flags

   !> Four branch source components for an L=0 field.  The matrix is already
   !> the physical adjoint source: for chi_plus it lives in down,up and for
   !> chi_minus in up,down.  The branch ordering is 1+p+2*q.
   subroutine sr_l0_source_components(space, radial_bases, source_field, circular_channel, components)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: source_field(:)
      integer, intent(in) :: circular_channel
      complex(rp), intent(out) :: components(:, :, :)

      call sr_l0_source_components_with_flags(space, radial_bases, source_field, circular_channel, .true., .true., components)
   end subroutine sr_l0_source_components

   subroutine sr_l0_source_components_with_flags(space, radial_bases, source_field, circular_channel, &
                                                  include_lower_radial, include_lower_angular, components)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: source_field(:)
      integer, intent(in) :: circular_channel
      logical, intent(in) :: include_lower_radial, include_lower_angular
      complex(rp), intent(out) :: components(:, :, :)
      type(response_super_index) :: item
      integer :: nsite, norb, n, flat, ir, site, l, m, channel, iorb, p, q, component
      integer :: spin_left, spin_right, row, col
      real(rp) :: gaunt, pair
      complex(rp) :: coefficient

      call validate_inputs(space, radial_bases, source_field, components, circular_channel)
      nsite = space%nsite
      norb = (radial_bases(1)%lmax + 1)**2
      n = 2*norb*nsite
      components = cmplx(0.0_rp, 0.0_rp, rp)
      if (circular_channel == 1) then
         spin_left = 1
         spin_right = 2
      else
         spin_left = 2
         spin_right = 1
      end if

      do flat = 1, size(source_field)
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, 1, item)
         if (item%response_l /= 0 .or. item%response_m /= 0) cycle
         site = item%site
         ir = item%radial_point
         channel = item%channel
         if (channel /= 1 .or. space%radial_weights(ir) == 0.0_rp) cycle
         coefficient = source_field(flat)
         if (coefficient == cmplx(0.0_rp, 0.0_rp, rp)) cycle
         do iorb = 1, norb
            l = lmto_orbital_l(iorb)
            m = iorb - l*l - l - 1
            gaunt = response_gaunt(l, m, l, m, 0, 0)
            do p = 0, 1
               do q = 0, 1
                  component = 1 + p + 2*q
                  pair = sr_mixed_transverse_point_with_flags(radial_bases(site), ir, l, spin_right, q, spin_left, p, &
                     include_lower_radial, include_lower_angular)
                  row = (site - 1)*2*norb + (spin_right - 1)*norb + iorb
                  col = (site - 1)*2*norb + (spin_left - 1)*norb + iorb
                  components(row, col, component) = components(row, col, component) + coefficient* &
                     space%radial_weights(ir)*gaunt*pair
               end do
            end do
         end do
      end do
   end subroutine sr_l0_source_components_with_flags

   subroutine sr_l0_source(space, radial_bases, source_field, circular_channel, hamiltonian, operator)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: source_field(:), hamiltonian(:, :)
      integer, intent(in) :: circular_channel
      complex(rp), intent(out) :: operator(:, :)
      complex(rp), allocatable :: components(:, :, :), term(:, :)
      integer :: n, p, q, component

      n = size(hamiltonian, 1)
      if (size(hamiltonian, 2) /= n .or. any(shape(operator) /= [n, n])) then
         error stop 'DRESP-09S: source Hamiltonian/operator shape mismatch'
      end if
      allocate(components(n, n, 4), term(n, n))
      call sr_l0_source_components(space, radial_bases, source_field, circular_channel, components)
      operator = cmplx(0.0_rp, 0.0_rp, rp)
      do p = 0, 1
         do q = 0, 1
            component = 1 + p + 2*q
            term = components(:, :, component)
            if (q == 1) term = matmul(hamiltonian, term)
            if (p == 1) term = matmul(term, hamiltonian)
            operator = operator + term
         end do
      end do
      deallocate(components, term)
   end subroutine sr_l0_source

   subroutine sr_l0_source_with_flags(space, radial_bases, source_field, circular_channel, hamiltonian, &
                                      include_lower_radial, include_lower_angular, operator)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: source_field(:), hamiltonian(:, :)
      integer, intent(in) :: circular_channel
      logical, intent(in) :: include_lower_radial, include_lower_angular
      complex(rp), intent(out) :: operator(:, :)
      complex(rp), allocatable :: components(:, :, :), term(:, :)
      integer :: n, p, q, component

      n = size(hamiltonian, 1)
      if (size(hamiltonian, 2) /= n .or. any(shape(operator) /= [n, n])) then
         error stop 'DRESP-09S: flagged source Hamiltonian/operator shape mismatch'
      end if
      allocate(components(n, n, 4), term(n, n))
      call sr_l0_source_components_with_flags(space, radial_bases, source_field, circular_channel, &
         include_lower_radial, include_lower_angular, components)
      operator = cmplx(0.0_rp, 0.0_rp, rp)
      do p = 0, 1
         do q = 0, 1
            component = 1 + p + 2*q
            term = components(:, :, component)
            if (q == 1) term = matmul(hamiltonian, term)
            if (p == 1) term = matmul(term, hamiltonian)
            operator = operator + term
         end do
      end do
      deallocate(components, term)
   end subroutine sr_l0_source_with_flags

   !> Build the pointwise L=0 transverse observable from a coefficient-space
   !> density matrix.  The returned value is the coefficient of Y00 in the
   !> observable, i.e. it is the exact dual of source_field above.
   subroutine sr_l0_density_from_matrix(space, radial_bases, density_matrix, circular_channel, density, &
                                        include_lower_radial, include_lower_angular)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: density_matrix(:, :)
      integer, intent(in) :: circular_channel
      complex(rp), intent(out) :: density(:, :)
      logical, intent(in), optional :: include_lower_radial, include_lower_angular
      integer :: nsite, norb, site, ir, iorb, l, m, spin_left, spin_right, row, col
      real(rp) :: gaunt
      logical :: use_lower_radial, use_lower_angular

      use_lower_radial = .true.
      use_lower_angular = .true.
      if (present(include_lower_radial)) use_lower_radial = include_lower_radial
      if (present(include_lower_angular)) use_lower_angular = include_lower_angular

      nsite = space%nsite
      norb = (radial_bases(1)%lmax + 1)**2
      if (size(radial_bases) /= nsite .or. any(shape(density_matrix) /= [2*norb*nsite, 2*norb*nsite]) .or. &
          any(shape(density) /= [nsite, space%npoint])) then
         error stop 'DRESP-09S: density matrix/output shape mismatch'
      end if
      if (circular_channel == 1) then
         spin_left = 1
         spin_right = 2
      else if (circular_channel == 2) then
         spin_left = 2
         spin_right = 1
      else
         error stop 'DRESP-09S: invalid circular channel'
      end if
      density = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, nsite
         do ir = 2, space%npoint
            do iorb = 1, norb
               l = lmto_orbital_l(iorb)
               m = iorb - l*l - l - 1
               gaunt = response_gaunt(l, m, l, m, 0, 0)
               row = (site - 1)*2*norb + (spin_left - 1)*norb + iorb
               col = (site - 1)*2*norb + (spin_right - 1)*norb + iorb
               density(site, ir) = density(site, ir) + gaunt*density_matrix(row, col)* &
                  sr_mixed_transverse_point_with_flags(radial_bases(site), ir, l, spin_left, 0, spin_right, 0, &
                     use_lower_radial, use_lower_angular)
            end do
         end do
      end do
   end subroutine sr_l0_density_from_matrix

   !> Density-side dual of the complete endpoint expansion.  H is the same
   !> accepted H^(2) coefficient matrix used by sr_l0_source; no radial
   !> finite difference is involved.  For branch (p,q), the effective density
   !> is H**p rho H**q, which is the cyclic trace dual of H**q O^H H**p.
   subroutine sr_l0_density_from_matrix_hamiltonian(space, radial_bases, density_matrix, circular_channel, hamiltonian, density, &
                                                    include_lower_radial, include_lower_angular)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: density_matrix(:, :), hamiltonian(:, :)
      integer, intent(in) :: circular_channel
      complex(rp), intent(out) :: density(:, :)
      logical, intent(in), optional :: include_lower_radial, include_lower_angular
      complex(rp), allocatable :: effective(:, :)
      integer :: nsite, norb, n, site, ir, iorb, l, m, p, q, row, col
      integer :: spin_left, spin_right
      real(rp) :: gaunt, pair
      logical :: use_lower_radial, use_lower_angular

      nsite = space%nsite
      norb = (radial_bases(1)%lmax + 1)**2
      n = 2*norb*nsite
      if (any(shape(density_matrix) /= [n, n]) .or. any(shape(hamiltonian) /= [n, n]) .or. &
          any(shape(density) /= [nsite, space%npoint])) error stop 'DRESP-09S: Hamiltonian density shape mismatch'
      use_lower_radial = .true.
      use_lower_angular = .true.
      if (present(include_lower_radial)) use_lower_radial = include_lower_radial
      if (present(include_lower_angular)) use_lower_angular = include_lower_angular
      if (circular_channel == 1) then
         spin_left = 1
         spin_right = 2
      else if (circular_channel == 2) then
         spin_left = 2
         spin_right = 1
      else
         error stop 'DRESP-09S: invalid Hamiltonian density circular channel'
      end if

      allocate(effective(n, n))
      density = cmplx(0.0_rp, 0.0_rp, rp)
      do p = 0, 1
         do q = 0, 1
            effective = density_matrix
            if (p == 1) effective = matmul(hamiltonian, effective)
            if (q == 1) effective = matmul(effective, hamiltonian)
            do site = 1, nsite
               do ir = 2, space%npoint
                  do iorb = 1, norb
                     l = lmto_orbital_l(iorb)
                     m = iorb - l*l - l - 1
                     gaunt = response_gaunt(l, m, l, m, 0, 0)
                     row = (site - 1)*2*norb + (spin_left - 1)*norb + iorb
                     col = (site - 1)*2*norb + (spin_right - 1)*norb + iorb
                     pair = sr_mixed_transverse_point_with_flags(radial_bases(site), ir, l, spin_left, p, spin_right, q, &
                        use_lower_radial, use_lower_angular)
                     density(site, ir) = density(site, ir) + gaunt*effective(row, col)*pair
                  end do
               end do
            end do
         end do
      end do
      deallocate(effective)
   end subroutine sr_l0_density_from_matrix_hamiltonian

   !> Complete density-side tangent of the endpoint object.  Unlike
   !> sr_l0_density_from_matrix_hamiltonian, this routine differentiates both
   !> endpoint Hamiltonian powers and the density matrix.  The branch algebra
   !> is deliberately supplied by the caller so the same service can be used
   !> with the k-summed production M0/M1/M2 object.
   subroutine sr_l0_density_tangent_from_hamiltonian(space, radial_bases, density_matrix, circular_channel, hamiltonian, &
                                                     delta_hamiltonian, delta_density_matrix, density, &
                                                     include_lower_radial, include_lower_angular)
      use lr_lmto_density_moment_tangent_mod, only: endpoint_tangent_branches
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: density_matrix(:, :), hamiltonian(:, :), delta_hamiltonian(:, :), delta_density_matrix(:, :)
      integer, intent(in) :: circular_channel
      complex(rp), intent(out) :: density(:, :)
      logical, intent(in), optional :: include_lower_radial, include_lower_angular
      complex(rp) :: branches(size(hamiltonian, 1), size(hamiltonian, 2), 4)

      call endpoint_tangent_branches(hamiltonian, density_matrix, delta_hamiltonian, delta_density_matrix, branches)
      call sr_l0_density_from_endpoint_branches(space, radial_bases, branches, circular_channel, density, &
         include_lower_radial, include_lower_angular)
   end subroutine sr_l0_density_tangent_from_hamiltonian

   !> Contract already-built complete endpoint branches against the certified
   !> scalar-relativistic radial bilinears.  Branch order is 00,10,01,11.
   !> Optional l bounds are diagnostic only and do not alter the production
   !> branch convention.
   subroutine sr_l0_density_from_endpoint_branches(space, radial_bases, endpoint_branches, circular_channel, density, &
                                                   include_lower_radial, include_lower_angular, l_first, l_last)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: endpoint_branches(:, :, :)
      integer, intent(in) :: circular_channel
      complex(rp), intent(out) :: density(:, :)
      logical, intent(in), optional :: include_lower_radial, include_lower_angular
      integer, intent(in), optional :: l_first, l_last
      integer :: nsite, norb, n, site, ir, iorb, l, m, p, q, branch, row, col, first_l, last_l
      integer :: spin_left, spin_right
      real(rp) :: gaunt, pair
      logical :: use_lower_radial, use_lower_angular

      nsite = space%nsite
      norb = (radial_bases(1)%lmax + 1)**2
      n = 2*norb*nsite
      if (size(radial_bases) /= nsite .or. any(shape(endpoint_branches) /= [n, n, 4]) .or. &
          any(shape(density) /= [nsite, space%npoint])) then
         error stop 'DRESP-09V endpoint density: shape mismatch'
      end if
      use_lower_radial = .true.
      use_lower_angular = .true.
      if (present(include_lower_radial)) use_lower_radial = include_lower_radial
      if (present(include_lower_angular)) use_lower_angular = include_lower_angular
      first_l = 0
      last_l = radial_bases(1)%lmax
      if (present(l_first)) first_l = max(0, l_first)
      if (present(l_last)) last_l = min(radial_bases(1)%lmax, l_last)
      if (first_l > last_l) then
         density = cmplx(0.0_rp, 0.0_rp, rp)
         return
      end if
      if (circular_channel == 1) then
         spin_left = 1
         spin_right = 2
      else if (circular_channel == 2) then
         spin_left = 2
         spin_right = 1
      else
         error stop 'DRESP-09V endpoint density: invalid circular channel'
      end if

      density = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, nsite
         do ir = 2, space%npoint
            do iorb = 1, norb
               l = lmto_orbital_l(iorb)
               if (l < first_l .or. l > last_l) cycle
               m = iorb - l*l - l - 1
               gaunt = response_gaunt(l, m, l, m, 0, 0)
               row = (site - 1)*2*norb + (spin_left - 1)*norb + iorb
               col = (site - 1)*2*norb + (spin_right - 1)*norb + iorb
               do p = 0, 1
                  do q = 0, 1
                     branch = 1 + p + 2*q
                     pair = sr_mixed_transverse_point_with_flags(radial_bases(site), ir, l, spin_left, p, spin_right, q, &
                        use_lower_radial, use_lower_angular)
                     density(site, ir) = density(site, ir) + gaunt*endpoint_branches(row, col, branch)*pair
                  end do
               end do
            end do
         end do
      end do
   end subroutine sr_l0_density_from_endpoint_branches

   !> Density-side contraction for the complete absolute-energy endpoint
   !> expansion.  This is intentionally separate from the historical
   !> four-branch DRESP-09S routine: branch 20/02 use phiddot through the
   !> explicit coefficient algebra in sr_second_order_branch_point.
   subroutine sr_l0_density_from_second_order_endpoint_branches(space, radial_bases, endpoint_branches, circular_channel, &
                                                               density, lower_factor, l_first, l_last)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: endpoint_branches(:, :, :)
      integer, intent(in) :: circular_channel
      complex(rp), intent(out) :: density(:, :)
      real(rp), intent(in), optional :: lower_factor
      integer, intent(in), optional :: l_first, l_last
      integer :: nsite, norb, n, site, ir, iorb, l, m, branch, row, col, first_l, last_l
      integer :: spin_left, spin_right
      real(rp) :: gaunt, pair, factor

      nsite = space%nsite
      norb = (radial_bases(1)%lmax + 1)**2
      n = 2*norb*nsite
      if (size(radial_bases) /= nsite .or. any(shape(endpoint_branches) /= [n, n, 6]) .or. &
          any(shape(density) /= [nsite, space%npoint])) then
         error stop 'DRESP-09X endpoint density: shape mismatch'
      end if
      factor = sr_l0_lower_spin_factor
      if (present(lower_factor)) factor = lower_factor
      first_l = 0
      last_l = radial_bases(1)%lmax
      if (present(l_first)) first_l = max(0, l_first)
      if (present(l_last)) last_l = min(radial_bases(1)%lmax, l_last)
      if (first_l > last_l) then
         density = cmplx(0.0_rp, 0.0_rp, rp)
         return
      end if
      if (circular_channel == 1) then
         spin_left = 1; spin_right = 2
      else if (circular_channel == 2) then
         spin_left = 2; spin_right = 1
      else
         error stop 'DRESP-09X endpoint density: invalid circular channel'
      end if

      density = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, nsite
         do ir = 2, space%npoint
            do iorb = 1, norb
               l = lmto_orbital_l(iorb)
               if (l < first_l .or. l > last_l) cycle
               m = iorb - l*l - l - 1
               gaunt = response_gaunt(l, m, l, m, 0, 0)
               row = (site - 1)*2*norb + (spin_left - 1)*norb + iorb
               col = (site - 1)*2*norb + (spin_right - 1)*norb + iorb
               do branch = 1, 6
                  pair = sr_second_order_branch_point(radial_bases(site), ir, l, spin_left, spin_right, branch, factor)
                  density(site, ir) = density(site, ir) + gaunt*endpoint_branches(row, col, branch)*pair
               end do
            end do
         end do
      end do
   end subroutine sr_l0_density_from_second_order_endpoint_branches

   !> Local observable matrix at one radius, without a field or quadrature.
   !> This is used by the radial field/density pairing oracle.
   subroutine sr_l0_observable_matrix(radial, ir, circular_channel, matrix)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, circular_channel
      complex(rp), intent(out) :: matrix(:, :)
      integer :: norb, iorb, l, spin_left, spin_right
      real(rp) :: gaunt

      norb = (radial%lmax + 1)**2
      if (any(shape(matrix) /= [2*norb, 2*norb])) error stop 'DRESP-09S: observable matrix shape mismatch'
      if (circular_channel == 1) then
         spin_left = 1
         spin_right = 2
      else if (circular_channel == 2) then
         spin_left = 2
         spin_right = 1
      else
         error stop 'DRESP-09S: invalid circular channel'
      end if
      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do iorb = 1, norb
         l = lmto_orbital_l(iorb)
         gaunt = response_gaunt(l, iorb - l*l - l - 1, l, iorb - l*l - l - 1, 0, 0)
         matrix((spin_left - 1)*norb + iorb, (spin_right - 1)*norb + iorb) = &
            cmplx(gaunt*sr_mixed_transverse_point(radial, ir, l, spin_left, 0, spin_right, 0), 0.0_rp, rp)
      end do
   end subroutine sr_l0_observable_matrix

   subroutine validate_inputs(space, radial_bases, source_field, components, circular_channel)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: source_field(:), components(:, :, :)
      integer, intent(in) :: circular_channel
      integer :: n, norb, isite

      norb = (radial_bases(1)%lmax + 1)**2
      n = 2*norb*space%nsite
      if (space%nchannel /= 1 .or. size(source_field) /= space%ndim .or. size(radial_bases) /= space%nsite .or. &
          any(shape(components) /= [n, n, 4]) .or. (circular_channel /= 1 .and. circular_channel /= 2)) then
         error stop 'DRESP-09S: source representation shape mismatch'
      end if
      if (space%response_lmax < 0 .or. space%npoint < 3) error stop 'DRESP-09S: incomplete L=0 response space'
      do isite = 1, size(radial_bases)
         if (radial_bases(isite)%lmax /= radial_bases(1)%lmax .or. radial_bases(isite)%nspin /= 2 .or. &
             radial_bases(isite)%npoint /= space%npoint .or. .not. allocated(radial_bases(isite)%tmc)) then
            error stop 'DRESP-09S: incomplete scalar-relativistic provenance'
         end if
      end do
   end subroutine validate_inputs

end module lr_dresp09s_scalar_relativistic_mod
