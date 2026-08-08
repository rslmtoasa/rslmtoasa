module gbt_structure_mod
   use precision_mod, only: rp
   implicit none
   private

   public :: gbt_bond_phase
   public :: gbt_reference_frame
   public :: gbt_endpoint_link
   public :: gbt_lift_orbital_block
   public :: gbt_contract_collinear

contains

   pure function gbt_bond_phase(q_cart_2pi_over_alat, d_cart, alat) result(alpha)
      real(rp), intent(in) :: q_cart_2pi_over_alat(3), d_cart(3), alat
      real(rp) :: alpha

      if (alat <= 0.0_rp) then
         alpha = huge(1.0_rp)
      else
         alpha = 2.0_rp*acos(-1.0_rp)*(q_cart_2pi_over_alat(1)*d_cart(1) + &
                  q_cart_2pi_over_alat(2)*d_cart(2) + q_cart_2pi_over_alat(3)*d_cart(3))/alat
      end if
   end function gbt_bond_phase

   pure subroutine gbt_reference_frame(theta, phi, frame)
      real(rp), intent(in) :: theta, phi
      complex(rp), intent(out) :: frame(2, 2)
      real(rp) :: c, s
      complex(rp) :: eiphi

      c = cos(0.5_rp*theta)
      s = sin(0.5_rp*theta)
      eiphi = cmplx(cos(phi), sin(phi), rp)
      frame(1, 1) = cmplx(c, 0.0_rp, rp)
      frame(1, 2) = -conjg(eiphi)*s
      frame(2, 1) = eiphi*s
      frame(2, 2) = cmplx(c, 0.0_rp, rp)
   end subroutine gbt_reference_frame

   pure subroutine gbt_endpoint_link(q_cart_2pi_over_alat, d_cart, alat, theta_a, phi_a, theta_b, phi_b, link)
      real(rp), intent(in) :: q_cart_2pi_over_alat(3), d_cart(3), alat
      real(rp), intent(in) :: theta_a, phi_a, theta_b, phi_b
      complex(rp), intent(out) :: link(2, 2)
      complex(rp) :: ua(2, 2), ub(2, 2), rz(2, 2), tmp(2, 2)
      real(rp) :: alpha

      call gbt_reference_frame(theta_a, phi_a, ua)
      call gbt_reference_frame(theta_b, phi_b, ub)
      alpha = gbt_bond_phase(q_cart_2pi_over_alat, d_cart, alat)
      rz = cmplx(0.0_rp, 0.0_rp, rp)
      rz(1, 1) = cmplx(cos(-0.5_rp*alpha), sin(-0.5_rp*alpha), rp)
      rz(2, 2) = cmplx(cos( 0.5_rp*alpha), sin( 0.5_rp*alpha), rp)
      tmp = matmul(rz, ub)
      link = matmul(transpose(conjg(ua)), tmp)
   end subroutine gbt_endpoint_link

   pure subroutine gbt_lift_orbital_block(s_orb, link, s_spinor)
      complex(rp), intent(in) :: s_orb(:, :), link(2, 2)
      complex(rp), intent(out) :: s_spinor(:, :)
      integer :: nrow, ncol, i, j, sa, sb

      nrow = size(s_orb, 1)
      ncol = size(s_orb, 2)
      s_spinor = cmplx(0.0_rp, 0.0_rp, rp)
      do sa = 1, 2
         do sb = 1, 2
            do i = 1, nrow
               do j = 1, ncol
                  s_spinor(i + (sa - 1)*nrow, j + (sb - 1)*ncol) = s_orb(i, j)*link(sa, sb)
               end do
            end do
         end do
      end do
   end subroutine gbt_lift_orbital_block

   pure subroutine gbt_contract_collinear(s_spinor, wx0_a, wx1_a, wx0_b, wx1_b, sign_a, sign_b, h_spinor)
      complex(rp), intent(in) :: s_spinor(:, :)
      complex(rp), intent(in) :: wx0_a(:), wx1_a(:), wx0_b(:), wx1_b(:)
      real(rp), intent(in) :: sign_a, sign_b
      complex(rp), intent(out) :: h_spinor(:, :)
      integer :: norb_a, norb_b, i, j, sa, sb
      complex(rp) :: wa, wb

      norb_a = size(wx0_a)
      norb_b = size(wx0_b)
      h_spinor = cmplx(0.0_rp, 0.0_rp, rp)
      do sa = 1, 2
         do sb = 1, 2
            do i = 1, norb_a
               wa = wx0_a(i) + merge(sign_a, -sign_a, sa == 1)*wx1_a(i)
               do j = 1, norb_b
                  wb = wx0_b(j) + merge(sign_b, -sign_b, sb == 1)*wx1_b(j)
                  h_spinor(i + (sa - 1)*norb_a, j + (sb - 1)*norb_b) = &
                     wa*s_spinor(i + (sa - 1)*norb_a, j + (sb - 1)*norb_b)*wb
               end do
            end do
         end do
      end do
   end subroutine gbt_contract_collinear

end module gbt_structure_mod
