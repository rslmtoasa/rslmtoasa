!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Authoritative second-order LMTO endpoint-product branch contract.
!
!  The production product basis is a Taylor polynomial in the two absolute
!  endpoint energies.  The expansion is truncated at total second order in
!  the two endpoint deviations; it is therefore a six-branch polynomial, not
!  the four affine endpoint products used by the historical representation.
!
!  This module is deliberately independent of the frozen DRESP-09X/Y oracle
!  implementations.  It is the live production algebra used by the compact
!  basis and by the reciprocal-GF vertex path.
!------------------------------------------------------------------------------
module lr_lmto_endpoint_branches_mod

   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   implicit none
   private

   integer, parameter, public :: lmto_product_branch_00 = 1
   integer, parameter, public :: lmto_product_branch_10 = 2
   integer, parameter, public :: lmto_product_branch_01 = 3
   integer, parameter, public :: lmto_product_branch_11 = 4
   integer, parameter, public :: lmto_product_branch_20 = 5
   integer, parameter, public :: lmto_product_branch_02 = 6
   integer, parameter, public :: lmto_product_nbranch = 6
   integer, parameter, public :: lmto_product_max_endpoint_power = 2
   integer, parameter, public :: lmto_product_max_gf_moment = 4

   public :: lmto_product_branch_powers
   public :: lmto_product_branch_label
   public :: lmto_product_branch_valid
   public :: lmto_product_energy_power
   public :: lmto_product_second_order_radial_branch

contains

   pure logical function lmto_product_branch_valid(branch) result(valid)
      integer, intent(in) :: branch

      valid = branch >= lmto_product_branch_00 .and. branch <= lmto_product_branch_02
   end function lmto_product_branch_valid

   pure subroutine lmto_product_branch_powers(branch, p, q)
      integer, intent(in) :: branch
      integer, intent(out) :: p, q

      p = -1
      q = -1
      select case (branch)
      case (lmto_product_branch_00)
         p = 0; q = 0
      case (lmto_product_branch_10)
         p = 1; q = 0
      case (lmto_product_branch_01)
         p = 0; q = 1
      case (lmto_product_branch_11)
         p = 1; q = 1
      case (lmto_product_branch_20)
         p = 2; q = 0
      case (lmto_product_branch_02)
         p = 0; q = 2
      end select
   end subroutine lmto_product_branch_powers

   pure character(len=2) function lmto_product_branch_label(branch) result(label)
      integer, intent(in) :: branch

      select case (branch)
      case (lmto_product_branch_00); label = '00'
      case (lmto_product_branch_10); label = '10'
      case (lmto_product_branch_01); label = '01'
      case (lmto_product_branch_11); label = '11'
      case (lmto_product_branch_20); label = '20'
      case (lmto_product_branch_02); label = '02'
      case default; label = '??'
      end select
   end function lmto_product_branch_label

   !> Stable real energy power for the endpoint and GF contracts.
   pure real(rp) function lmto_product_energy_power(energy, power) result(value)
      real(rp), intent(in) :: energy
      integer, intent(in) :: power
      integer :: n

      if (power < 0 .or. power > lmto_product_max_gf_moment) then
         error stop 'lmto_product_energy_power: unsupported power'
      end if
      value = 1.0_rp
      do n = 1, power
         value = value*energy
      end do
   end function lmto_product_energy_power

   !> Evaluate one radial coefficient B_pq of the certified second-order
   !> absolute-energy product polynomial.  The expression is
   !>
   !>   R_L R_R |_2 = phi_L phi_R
   !>       + dL dotphi_L phi_R + phi_L dR dotphi_R
   !>       + dL dR dotphi_L dotphi_R
   !>       + 1/2 dL^2 phiddot_L phi_R
   !>       + 1/2 dR^2 phi_L phiddot_R,
   !>
   !> with dL=E_L-Enu_L and dR=E_R-Enu_R.  Expanding in absolute energies
   !> gives exactly the six branches 00/10/01/11/20/02.
   recursive subroutine lmto_product_second_order_radial_branch(radial, ir, l, lp, spin_left, spin_right, branch, value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, lp, spin_left, spin_right, branch
      real(rp), intent(out) :: value
      real(rp) :: phi_l, phi_r, dot_l, dot_r, ddot_l, ddot_r
      real(rp) :: enu_l, enu_r, positive_value

      if (.not. lmto_product_branch_valid(branch)) error stop &
         'lmto_product_second_order_radial_branch: invalid branch'
      if (ir < 1 .or. ir > radial%npoint .or. l < 0 .or. l > radial%lmax .or. lp < 0 .or. lp > radial%lmax .or. &
          spin_left < 1 .or. spin_left > radial%nspin .or. spin_right < 1 .or. spin_right > radial%nspin) then
         error stop 'lmto_product_second_order_radial_branch: index outside radial basis'
      end if

      if (ir == 1) then
         if (l /= 0 .or. lp /= 0) then
            value = 0.0_rp
            return
         end if
         call lmto_product_second_order_radial_branch(radial, 2, l, lp, spin_left, spin_right, branch, positive_value)
         call lmto_product_second_order_radial_branch(radial, 3, l, lp, spin_left, spin_right, branch, value)
         ! The origin value is the same r^2 extrapolation used by the
         ! historical response-space product evaluator.
         value = (positive_value*radial%rofi(3)**2 - value*radial%rofi(2)**2)/ &
            (radial%rofi(3)**2 - radial%rofi(2)**2)
         return
      end if
      if (radial%rofi(ir) <= tiny(1.0_rp)) error stop &
         'lmto_product_second_order_radial_branch: invalid positive radius'

      phi_l = radial%phi_large(ir, l + 1, spin_left)
      phi_r = radial%phi_large(ir, lp + 1, spin_right)
      dot_l = radial%phidot_large(ir, l + 1, spin_left)
      dot_r = radial%phidot_large(ir, lp + 1, spin_right)
      ddot_l = radial%phiddot_large(ir, l + 1, spin_left)
      ddot_r = radial%phiddot_large(ir, lp + 1, spin_right)
      enu_l = radial%enu_work(l + 1, spin_left)
      enu_r = radial%enu_work(lp + 1, spin_right)

      select case (branch)
      case (lmto_product_branch_00)
         value = phi_l*phi_r - enu_l*dot_l*phi_r - enu_r*phi_l*dot_r + enu_l*enu_r*dot_l*dot_r + &
            0.5_rp*enu_l**2*ddot_l*phi_r + 0.5_rp*enu_r**2*phi_l*ddot_r
      case (lmto_product_branch_10)
         value = dot_l*phi_r - enu_r*dot_l*dot_r - enu_l*ddot_l*phi_r
      case (lmto_product_branch_01)
         value = phi_l*dot_r - enu_l*dot_l*dot_r - enu_r*phi_l*ddot_r
      case (lmto_product_branch_11)
         value = dot_l*dot_r
      case (lmto_product_branch_20)
         value = 0.5_rp*ddot_l*phi_r
      case (lmto_product_branch_02)
         value = 0.5_rp*phi_l*ddot_r
      end select
      value = value/radial%rofi(ir)**2
   end subroutine lmto_product_second_order_radial_branch

end module lr_lmto_endpoint_branches_mod
