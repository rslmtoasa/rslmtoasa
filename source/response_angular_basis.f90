!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Angular algebra for the LR-02 response-basis audit.
!>
!> This module owns only the complex spherical-harmonic convention, Gaunt
!> algebra, and product-space cutoff.  It does not construct response
!> vertices or a susceptibility.
!>
!> The live `hcpx` transformation uses the complex convention
!> Y(code)_lm = conjugate(Y(Condon--Shortley)_lm), in the order
!> (l,m)=(0,0),(1,-1),(1,0),(1,1),... .  Keeping this convention here avoids
!> importing a generic Gaunt phase into the response layer.
!------------------------------------------------------------------------------
module response_angular_basis_mod

   use precision_mod, only: rp
   implicit none
   private

   real(rp), parameter, public :: response_angular_pi = &
      3.1415926535897932384626433832795_rp

   public :: response_lmax
   public :: response_product_space_complete
   public :: response_lm_index
   public :: response_lm_from_index
   public :: response_harmonic
   public :: response_gaunt
   public :: response_product_coefficients

contains

   !> The complete product of an orbital basis truncated at lmax has this
   !> response angular cutoff unless symmetry removes channels.
   pure integer function response_lmax(orbital_lmax) result(value)
      integer, intent(in) :: orbital_lmax
      value = 2*orbital_lmax
   end function response_lmax

   pure logical function response_product_space_complete(orbital_lmax, response_lmax_in) result(ok)
      integer, intent(in) :: orbital_lmax, response_lmax_in
      ok = orbital_lmax >= 0 .and. response_lmax_in >= response_lmax(orbital_lmax)
   end function response_product_space_complete

   !> One-based index in the live (l,m) order used by hcpx for complex sp/spd.
   pure integer function response_lm_index(l, m) result(index)
      integer, intent(in) :: l, m
      if (l < 0 .or. abs(m) > l) then
         index = -1
      else
         index = l*l + l + m + 1
      end if
   end function response_lm_index

   pure subroutine response_lm_from_index(index, l, m)
      integer, intent(in) :: index
      integer, intent(out) :: l, m
      integer :: trial

      l = -1
      m = 0
      if (index < 1) return
      do trial = 0, 64
         if (trial*trial + 1 <= index .and. index <= (trial + 1)*(trial + 1)) then
            l = trial
            m = index - trial*trial - trial - 1
            return
         end if
      end do
   end subroutine response_lm_from_index

   !> Complex spherical harmonic in the convention implied by hcpx.
   pure function response_harmonic(l, m, theta, phi) result(value)
      integer, intent(in) :: l, m
      real(rp), intent(in) :: theta, phi
      complex(rp) :: value
      complex(rp) :: standard
      real(rp) :: x, p, norm
      integer :: ma

      value = cmplx(0.0_rp, 0.0_rp, rp)
      if (l < 0 .or. abs(m) > l) return

      x = cos(theta)
      ma = abs(m)
      p = associated_legendre_without_condon(l, ma, x)
      norm = sqrt(real(2*l + 1, rp)/(4.0_rp*response_angular_pi) * &
         factorial_real(l - ma)/factorial_real(l + ma))
      if (m >= 0) then
         standard = cmplx(norm*minus_one_power(ma)*p*cos(real(m, rp)*phi), &
            norm*minus_one_power(ma)*p*sin(real(m, rp)*phi), rp)
      else
         ! Condon--Shortley Y(l,-m)=(-1)^m conjugate(Y(l,m)).
         standard = minus_one_power(ma)*conjg(cmplx(norm*minus_one_power(ma)*p*cos(real(ma, rp)*phi), &
            norm*minus_one_power(ma)*p*sin(real(ma, rp)*phi), rp))
      end if

      ! hcpx columns are the conjugates of the usual Condon--Shortley
      ! harmonics; see the explicit p and d columns in source/math.f90.
      value = conjg(standard)
   end function response_harmonic

   !> Gaunt coefficient for the live hcpx harmonic convention:
   !>
   !>   integral Y(code)_LM^* Y(code)_lm^* Y(code)_l'm' dOmega.
   !>
   !> The expression is the Condon--Shortley 3j formula transformed by
   !> Y(code)=conjugate(Y(CS)).
   pure real(rp) function response_gaunt(l1, m1, l2, m2, lout, mout) result(value)
      integer, intent(in) :: l1, m1, l2, m2, lout, mout
      real(rp) :: prefactor

      value = 0.0_rp
      if (l1 < 0 .or. l2 < 0 .or. lout < 0) return
      if (abs(m1) > l1 .or. abs(m2) > l2 .or. abs(mout) > lout) return
      if (mout /= m2 - m1) return
      if (lout < abs(l1 - l2) .or. lout > l1 + l2) return
      if (mod(lout + l1 + l2, 2) /= 0) return

      prefactor = sqrt(real((2*lout + 1)*(2*l1 + 1)*(2*l2 + 1), rp) / &
         (4.0_rp*response_angular_pi))
      value = minus_one_power(m2)*prefactor* &
         wigner_3j(lout, l1, l2, 0, 0, 0)*wigner_3j(lout, l1, l2, mout, m1, -m2)
   end function response_gaunt

   !> Fill coefficients for Y(code)_lm^* Y(code)_l'm' in response order.
   subroutine response_product_coefficients(l, m, lp, mp, coefficients)
      integer, intent(in) :: l, m, lp, mp
      real(rp), intent(out) :: coefficients(:)
      integer :: lout, mout, expected_size

      expected_size = (l + lp + 1)**2
      if (l < 0 .or. lp < 0 .or. abs(m) > l .or. abs(mp) > lp .or. &
          size(coefficients) /= expected_size) then
         error stop 'response_product_coefficients: invalid product dimensions'
      end if
      coefficients = 0.0_rp
      do lout = 0, l + lp
         do mout = -lout, lout
            coefficients(response_lm_index(lout, mout)) = response_gaunt(l, m, lp, mp, lout, mout)
         end do
      end do
   end subroutine response_product_coefficients

   pure real(rp) function factorial_real(n) result(value)
      integer, intent(in) :: n
      integer :: i

      if (n < 0) then
         value = 0.0_rp
         return
      end if
      value = 1.0_rp
      do i = 2, n
         value = value*real(i, rp)
      end do
   end function factorial_real

   pure real(rp) function minus_one_power(n) result(value)
      integer, intent(in) :: n
      if (mod(abs(n), 2) == 0) then
         value = 1.0_rp
      else
         value = -1.0_rp
      end if
   end function minus_one_power

   pure real(rp) function associated_legendre_without_condon(l, m, x) result(value)
      integer, intent(in) :: l, m
      real(rp), intent(in) :: x
      real(rp) :: pmm, pmmp1, pll, somx2
      integer :: ell

      value = 0.0_rp
      if (m < 0 .or. l < m) return
      pmm = 1.0_rp
      if (m > 0) then
         somx2 = sqrt(max(0.0_rp, (1.0_rp - x)*(1.0_rp + x)))
         do ell = 1, m
            pmm = pmm*real(2*ell - 1, rp)*somx2
         end do
      end if
      if (l == m) then
         value = pmm
         return
      end if
      pmmp1 = x*real(2*m + 1, rp)*pmm
      if (l == m + 1) then
         value = pmmp1
         return
      end if
      do ell = m + 2, l
         pll = (real(2*ell - 1, rp)*x*pmmp1 - real(ell + m - 1, rp)*pmm)/real(ell - m, rp)
         pmm = pmmp1
         pmmp1 = pll
      end do
      value = pmmp1
   end function associated_legendre_without_condon

   pure real(rp) function wigner_3j(j1, j2, j3, m1, m2, m3) result(value)
      integer, intent(in) :: j1, j2, j3, m1, m2, m3
      integer :: z, zmin, zmax
      real(rp) :: delta, norm, denominator, sum_term

      value = 0.0_rp
      if (j1 < 0 .or. j2 < 0 .or. j3 < 0) return
      if (abs(m1) > j1 .or. abs(m2) > j2 .or. abs(m3) > j3) return
      if (m1 + m2 + m3 /= 0) return
      if (j3 < abs(j1 - j2) .or. j3 > j1 + j2) return

      delta = sqrt(factorial_real(j1 + j2 - j3)*factorial_real(j1 - j2 + j3)* &
         factorial_real(-j1 + j2 + j3)/factorial_real(j1 + j2 + j3 + 1))
      norm = sqrt(factorial_real(j1 + m1)*factorial_real(j1 - m1)* &
         factorial_real(j2 + m2)*factorial_real(j2 - m2)* &
         factorial_real(j3 + m3)*factorial_real(j3 - m3))
      zmin = max(0, max(-j3 + j2 - m1, -j3 + j1 + m2))
      zmax = min(j1 + j2 - j3, min(j1 - m1, j2 + m2))
      if (zmin > zmax) return

      sum_term = 0.0_rp
      do z = zmin, zmax
         denominator = factorial_real(z)*factorial_real(j1 + j2 - j3 - z)* &
            factorial_real(j1 - m1 - z)*factorial_real(j2 + m2 - z)* &
            factorial_real(j3 - j2 + m1 + z)*factorial_real(j3 - j1 - m2 + z)
         sum_term = sum_term + minus_one_power(z)/denominator
      end do
      value = minus_one_power(j1 - j2 - m3)*delta*norm*sum_term
   end function wigner_3j

end module response_angular_basis_mod
