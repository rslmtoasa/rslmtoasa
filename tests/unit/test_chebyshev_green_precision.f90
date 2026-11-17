!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
!> @brief Double-precision DOS reconstruction against a literal evaluation of the Jackson-damped sum.
!> @details Moments come from a dense Chebyshev recurrence on a 4x4 Hermitian fixture. The reference sums
!>          G(E) = sum_n g_n c_n (-i) exp(-i n theta) mu_n / sqrt(a^2 - (E-b)^2), theta = acos((E-b)/a),
!>          c_0 = 1, c_n = 2, with the Jackson coefficients g_n of Weisse et al., Rev. Mod. Phys. 78, 275 (2006), Eq. (71).
!>          cheb_green_dp must agree to rounding; the single-precision cheb_green_fast error is printed for scale.
!>          The optional damping argument is checked with the Lorentz kernel (lambda = 3.5,
!>          g_n = sinh(lambda (1 - n/N))/sinh(lambda)) and with no kernel (g_n = 1), written here in closed
!>          form, for cheb_green_dp, cheb_green_fast and the complex-energy cheb_green_complex. A damping
!>          argument that is ignored, inverted or indexed off by one gives errors of order one.
!------------------------------------------------------------------------------
program test_chebyshev_green_precision
   use precision_mod, only: rp
   use chebyshev_fast_mod, only: cheb_green_dp, cheb_green_fast, cheb_green_complex
   use math_mod, only: lorentz_kernel
   implicit none

   integer, parameter :: nmat = 4, n_mom = 64, nv = 25
   real(rp), parameter :: a = 1.2_rp, b = 0.1_rp, pi = 3.14159265358979323846_rp
   real(rp), parameter :: tolerance = 1.0e-12_rp
   complex(rp) :: h(nmat, nmat), moments(nmat, nmat, n_mom, 1), g_dp(nmat, nmat, nv, 1), g_sp(nmat, nmat, nv, 1)
   complex(rp) :: g_ref(nmat, nmat), g_cplx(nmat, nmat), z
   real(rp) :: ene(nv), err_dp, err_sp, scale, err_k, err_c, kern(n_mom), kern_math(n_mom), err_math
   real(rp), parameter :: lambda = 3.5_rp
   integer :: ie, ik
   character(len=8), parameter :: kernel_names(3) = ['jackson ', 'lorentz ', 'none    ']

   call build_fixture(h)
   call dense_chebyshev_moments(h, moments(:, :, :, 1))
   do ie = 1, nv
      ene(ie) = b - 0.95_rp*a + 1.9_rp*a*real(ie - 1, rp)/real(nv - 1, rp)
   end do
   call cheb_green_dp(moments, nmat, n_mom, 1, ene, nv, a, b, g_dp)
   call cheb_green_fast(moments, nmat, n_mom, 1, ene, nv, a, b, g_sp)

   err_dp = 0.0_rp
   err_sp = 0.0_rp
   scale = 0.0_rp
   do ie = 1, nv
      call literal_green(moments(:, :, :, 1), ene(ie), closed_kernel('jackson '), g_ref)
      err_dp = max(err_dp, maxval(abs(g_dp(:, :, ie, 1) - g_ref)))
      err_sp = max(err_sp, maxval(abs(g_sp(:, :, ie, 1) - g_ref)))
      scale = max(scale, maxval(abs(g_ref)))
   end do
   write (*, '(a,es12.4,a,es12.4)') 'max |G_dp - ref| / max|ref| = ', err_dp/scale, '   (single precision: ', err_sp/scale
   if (err_dp/scale > tolerance) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if

   ! math_mod's Lorentz kernel against the closed form
   call lorentz_kernel(n_mom, kern_math, lambda)
   err_math = maxval(abs(kern_math - closed_kernel('lorentz ')))
   write (*, '(a,es12.4)') 'max |lorentz_kernel - closed form| = ', err_math
   if (err_math > 1.0e-14_rp) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if

   ! Selectable damping: dp (rounding), fast (single precision) and complex energy (rounding)
   z = cmplx(0.15_rp, 0.07_rp, rp)
   do ik = 2, 3
      kern = closed_kernel(kernel_names(ik))
      call cheb_green_dp(moments, nmat, n_mom, 1, ene, nv, a, b, g_dp, kern)
      call cheb_green_fast(moments, nmat, n_mom, 1, ene, nv, a, b, g_sp, kern)
      err_dp = 0.0_rp
      err_sp = 0.0_rp
      scale = 0.0_rp
      do ie = 1, nv
         call literal_green(moments(:, :, :, 1), ene(ie), kern, g_ref)
         err_dp = max(err_dp, maxval(abs(g_dp(:, :, ie, 1) - g_ref)))
         err_sp = max(err_sp, maxval(abs(g_sp(:, :, ie, 1) - g_ref)))
         scale = max(scale, maxval(abs(g_ref)))
      end do
      call literal_green_complex(moments(:, :, :, 1), z, kern, g_ref)
      call cheb_green_complex(moments(:, :, :, 1), nmat, n_mom, z, a, b, g_cplx, kern)
      err_c = maxval(abs(g_cplx - g_ref))/maxval(abs(g_ref))
      err_k = err_dp/scale
      write (*, '(a,a,a,es10.3,a,es10.3,a,es10.3)') 'kernel ', trim(kernel_names(ik)), ': dp ', err_k, &
         '  fast ', err_sp/scale, '  complex ', err_c
      if (err_k > tolerance .or. err_sp/scale > 1.0e-4_rp .or. err_c > tolerance) then
         write (*, '(a)') 'RESULT: FAIL'
         error stop 1
      end if
   end do
   write (*, '(a)') 'RESULT: PASS'

contains

   !> Damping coefficients g_n, n = 0..n_mom-1, in closed form (Weisse et al., Eqs. (71) and (80)).
   function closed_kernel(name) result(gk)
      character(len=8), intent(in) :: name
      real(rp) :: gk(n_mom)
      real(rp) :: angle
      integer :: n

      do n = 0, n_mom - 1
         select case (trim(name))
         case ('jackson')
            angle = pi*real(n, rp)/real(n_mom + 1, rp)
            gk(n + 1) = (real(n_mom - n + 1, rp)*cos(angle) + sin(angle)/tan(pi/real(n_mom + 1, rp)))/real(n_mom + 1, rp)
         case ('lorentz')
            gk(n + 1) = sinh(lambda*(1.0_rp - real(n, rp)/real(n_mom, rp)))/sinh(lambda)
         case default
            gk(n + 1) = 1.0_rp
         end select
      end do
   end function closed_kernel

   subroutine literal_green(mu, energy, gk, g)
      complex(rp), intent(in) :: mu(nmat, nmat, n_mom)
      real(rp), intent(in) :: energy, gk(n_mom)
      complex(rp), intent(out) :: g(nmat, nmat)
      real(rp) :: theta, gn, cn
      integer :: n

      theta = acos((energy - b)/a)
      g = (0.0_rp, 0.0_rp)
      do n = 0, n_mom - 1
         gn = gk(n + 1)
         cn = merge(1.0_rp, 2.0_rp, n == 0)
         g = g + gn*cn*cmplx(0.0_rp, -1.0_rp, rp)*exp(cmplx(0.0_rp, -real(n, rp)*theta, rp))*mu(:, :, n + 1) &
                 /sqrt(a*a - (energy - b)**2)
      end do
   end subroutine literal_green

   subroutine literal_green_complex(mu, zz, gk, g)
      complex(rp), intent(in) :: mu(nmat, nmat, n_mom), zz
      real(rp), intent(in) :: gk(n_mom)
      complex(rp), intent(out) :: g(nmat, nmat)
      complex(rp) :: theta
      integer :: n

      theta = acos((zz - b)/a)
      g = (0.0_rp, 0.0_rp)
      do n = 0, n_mom - 1
         g = g + gk(n + 1)*merge(1.0_rp, 2.0_rp, n == 0)*cmplx(0.0_rp, -1.0_rp, rp)*exp(-cmplx(0.0_rp, 1.0_rp, rp)*real(n, rp)*theta) &
                 *mu(:, :, n + 1)/sqrt(a*a - (zz - b)**2)
      end do
   end subroutine literal_green_complex

   subroutine build_fixture(matrix)
      complex(rp), intent(out) :: matrix(:, :)

      matrix = (0.0_rp, 0.0_rp)
      matrix(1, 1) = (-0.70_rp, 0.0_rp)
      matrix(2, 2) = (-0.20_rp, 0.0_rp)
      matrix(3, 3) = (0.35_rp, 0.0_rp)
      matrix(4, 4) = (0.80_rp, 0.0_rp)
      matrix(1, 2) = (0.12_rp, 0.08_rp); matrix(2, 1) = conjg(matrix(1, 2))
      matrix(1, 3) = (0.03_rp, 0.0_rp); matrix(3, 1) = conjg(matrix(1, 3))
      matrix(2, 3) = (0.10_rp, 0.0_rp); matrix(3, 2) = conjg(matrix(2, 3))
      matrix(2, 4) = (0.0_rp, 0.04_rp); matrix(4, 2) = conjg(matrix(2, 4))
      matrix(3, 4) = (0.09_rp, -0.02_rp); matrix(4, 3) = conjg(matrix(3, 4))
   end subroutine build_fixture

   subroutine dense_chebyshev_moments(matrix, moments_out)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp), intent(out) :: moments_out(:, :, :)
      complex(rp) :: scaled(nmat, nmat), identity(nmat, nmat)
      integer :: i, order

      identity = (0.0_rp, 0.0_rp)
      do i = 1, nmat
         identity(i, i) = (1.0_rp, 0.0_rp)
      end do
      scaled = (matrix - b*identity)/a
      moments_out(:, :, 1) = identity
      moments_out(:, :, 2) = scaled
      do order = 3, n_mom
         moments_out(:, :, order) = 2.0_rp*matmul(scaled, moments_out(:, :, order - 1)) - moments_out(:, :, order - 2)
      end do
   end subroutine dense_chebyshev_moments

end program test_chebyshev_green_precision
