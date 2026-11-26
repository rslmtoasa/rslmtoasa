!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: cheb_sigma_kernel_mod
!
!> @brief Chebyshev expansion of G(z) = <i|(z - H - Sigma)^-1|i> with a non-Hermitian
!> site-local Sigma, using a complex centre b and a real half-width a.
!> @details Plain-array kernels: spectral window, root of the generating function,
!> log of the expansion coefficient, and the log-scaled coefficient-moment sum.
!------------------------------------------------------------------------------
module cheb_sigma_kernel_mod
   use precision_mod, only: rp
   implicit none
   private
   public :: cheb_sigma_window, cheb_sigma_w, cheb_sigma_log_coeff, cheb_sigma_sum

contains

   !> Ellipse parameter |u + sqrt(u^2-1)| with the branch giving rho >= 1.
   function ellipse_rho(u) result(rho)
      complex(rp), intent(in) :: u
      real(rp) :: rho
      complex(rp) :: s
      s = sqrt(u*u - 1.0_rp)
      rho = max(abs(u + s), abs(u - s))
   end function ellipse_rho

   !> Smallest and largest eigenvalue of a Hermitian matrix (zheev, no vectors).
   subroutine hermitian_range(m, lo, hi)
      complex(rp), intent(in) :: m(:, :)
      real(rp), intent(out) :: lo, hi
      complex(rp) :: work_query(1)
      complex(rp), allocatable :: work(:), mc(:, :)
      real(rp), allocatable :: ev(:), rwork(:)
      integer :: n, info
      n = size(m, 1)
      mc = m
      allocate (ev(n), rwork(max(1, 3*n - 2)))
      call zheev('N', 'U', n, mc, n, ev, work_query, -1, rwork, info)
      allocate (work(max(1, int(real(work_query(1), rp)))))
      call zheev('N', 'U', n, mc, n, ev, work, size(work), rwork, info)
      if (info /= 0) error stop 'cheb_sigma_kernel: zheev failed'
      lo = minval(ev)
      hi = maxval(ev)
   end subroutine hermitian_range

   !> Scale a, complex centre b and ellipse parameter rho_s of the rectangle that
   !> encloses the numerical range of H + Sigma(z); sig has shape (nb,nb,ntype).
   subroutine cheb_sigma_window(emin, emax, sig, a, b, rho_s)
      real(rp), intent(in) :: emin, emax
      complex(rp), intent(in) :: sig(:, :, :)
      real(rp), intent(out) :: a, rho_s
      complex(rp), intent(out) :: b
      real(rp) :: re_lo, re_hi, im_lo, im_hi, lo, hi
      complex(rp) :: corner(4)
      integer :: t, k
      re_lo = huge(1.0_rp); re_hi = -huge(1.0_rp)
      im_lo = huge(1.0_rp); im_hi = -huge(1.0_rp)
      do t = 1, size(sig, 3)
         call hermitian_range(0.5_rp*(sig(:, :, t) + conjg(transpose(sig(:, :, t)))), lo, hi)
         re_lo = min(re_lo, emin + lo); re_hi = max(re_hi, emax + hi)
         call hermitian_range((0.0_rp, -0.5_rp)*(sig(:, :, t) - conjg(transpose(sig(:, :, t)))), lo, hi)
         im_lo = min(im_lo, lo); im_hi = max(im_hi, hi)
      end do
      a = (re_hi - re_lo)/(2.0_rp - 0.3_rp)
      b = cmplx(0.5_rp*(re_lo + re_hi), 0.5_rp*(im_lo + im_hi), rp)
      corner = [cmplx(re_lo, im_lo, rp), cmplx(re_lo, im_hi, rp), cmplx(re_hi, im_lo, rp), cmplx(re_hi, im_hi, rp)]
      rho_s = 0.0_rp
      do k = 1, 4
         rho_s = max(rho_s, ellipse_rho((corner(k) - b)/a))
      end do
   end subroutine cheb_sigma_window

   !> Root w = zt - sqrt(zt^2-1) with |w| <= 1, zt = (z-b)/a, and rho_z = 1/|w|.
   subroutine cheb_sigma_w(z, a, b, w, rho_z)
      complex(rp), intent(in) :: z, b
      real(rp), intent(in) :: a
      complex(rp), intent(out) :: w
      real(rp), intent(out) :: rho_z
      complex(rp) :: zt, s, big
      zt = (z - b)/a
      s = sqrt(zt*zt - 1.0_rp)
      big = zt + s
      if (abs(zt - s) > abs(big)) big = zt - s
      w = 1.0_rp/big
      rho_z = abs(big)
   end subroutine cheb_sigma_w

   !> log of c_n = (2 - delta_n0) 2 w^(n+1)/(1 - w^2).
   function cheb_sigma_log_coeff(n, w) result(lc)
      integer, intent(in) :: n
      complex(rp), intent(in) :: w
      complex(rp) :: lc
      lc = log(2.0_rp) + (n + 1)*log(w) - log(1.0_rp - w*w)
      if (n > 0) lc = lc + log(2.0_rp)
   end function cheb_sigma_log_coeff

   !> g = (1/a) sum_k c_(k-1) mu(:,:,k) exp(logscale(k)); the exponents are combined
   !> before exponentiating. tail is the relative geometric tail bound
   !> max|last term| r/(1-r)/max|g| with r = rho_s/rho_z (huge if r >= 1).
   subroutine cheb_sigma_sum(mu, logscale, w, a, r, g, tail)
      complex(rp), intent(in) :: mu(:, :, :), w
      real(rp), intent(in) :: logscale(:), a, r
      complex(rp), intent(out) :: g(:, :)
      real(rp), intent(out) :: tail
      complex(rp) :: lc, term(size(g, 1), size(g, 2))
      real(rp) :: last
      integer :: k, nmom
      nmom = size(mu, 3)
      g = (0.0_rp, 0.0_rp)
      last = 0.0_rp
      do k = 1, nmom
         lc = cheb_sigma_log_coeff(k - 1, w)
         if (real(lc, rp) + logscale(k) < log(tiny(1.0_rp))) cycle
         term = exp(cmplx(real(lc, rp) + logscale(k), aimag(lc), rp))*mu(:, :, k)/a
         g = g + term
         if (k == nmom) last = maxval(abs(term))
      end do
      if (r < 0.0_rp .or. r >= 1.0_rp) then
         tail = huge(1.0_rp)
      else
         tail = last*r/(1.0_rp - r)/max(maxval(abs(g)), tiny(1.0_rp))
      end if
   end subroutine cheb_sigma_sum

end module cheb_sigma_kernel_mod
