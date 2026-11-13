!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
!> @brief Reconstruction of linear-prediction-extended moments without a damping kernel.
!> @details Damped-pole moments mu_n = sum_k w_k rho_k^n cos(n theta_k) have the closed-form density
!>          A(x) = sum_k w_k [P(t-t_k) + P(t+t_k)]/(2 pi sin t), x = cos t, P(f) = (1-r^2)/(1-2 r cos f+r^2).
!>          N = 302 moments are extended to 604 by cheb_lp_extend and reconstructed by cheb_green_dp with the
!>          Jackson kernel (default) and with unit damping (cheb_kernel = 'none'). The moments decay
!>          exponentially, the regime of Ganahl et al., PRB 90, 045144, where a kernel only costs resolution:
!>          the RMS error on |x| < 0.8 without a kernel must be below 0.25 times the Jackson error (a numpy
!>          toy of the same kind gave 0.0044 against 0.2055). A damping argument that is ignored gives ratio 1.
!------------------------------------------------------------------------------
program test_cheb_kernel
   use precision_mod, only: rp
   use cheb_lp_mod, only: cheb_lp_extend
   use chebyshev_fast_mod, only: cheb_green_dp
   implicit none

   integer, parameter :: n_known = 302, n_ext = 604, nk = 3, nx = 401
   real(rp), parameter :: pi = 3.14159265358979323846_rp
   real(rp), parameter :: rho(nk) = [0.970_rp, 0.985_rp, 0.990_rp]
   real(rp), parameter :: th(nk) = [0.90_rp, 1.40_rp, 2.10_rp]
   real(rp), parameter :: w(nk) = [1.0_rp, 0.8_rp, 1.2_rp]
   complex(rp) :: mu(1, 1, n_ext, 1), g_j(1, 1, nx, 1), g_n(1, 1, nx, 1)
   real(rp) :: x(nx), ref(nx), t, err_j, err_n, lam_max, fit_res, imag_max, uw
   integer :: i, k, n, n_rank, n_unstable

   do i = 1, nx
      x(i) = -0.8_rp + 1.6_rp*real(i - 1, rp)/real(nx - 1, rp)
      t = acos(x(i))
      ref(i) = 0.0_rp
      do k = 1, nk
         ref(i) = ref(i) + w(k)*(poisson(rho(k), t - th(k)) + poisson(rho(k), t + th(k)))/(2.0_rp*pi*sin(t))
      end do
   end do
   mu = (0.0_rp, 0.0_rp)
   do n = 0, n_known - 1
      do k = 1, nk
         mu(1, 1, n + 1, 1) = mu(1, 1, n + 1, 1) + w(k)*rho(k)**n*cos(n*th(k))
      end do
   end do
   call cheb_lp_extend(mu(:, :, :, 1), n_known, n_rank, n_unstable, lam_max, fit_res, imag_max, uw)
   call cheb_green_dp(mu, 1, n_ext, 1, x, nx, 1.0_rp, 0.0_rp, g_j)
   call cheb_green_dp(mu, 1, n_ext, 1, x, nx, 1.0_rp, 0.0_rp, g_n, spread(1.0_rp, 1, n_ext))
   err_j = sqrt(sum((-aimag(g_j(1, 1, :, 1))/pi - ref)**2)/nx)
   err_n = sqrt(sum((-aimag(g_n(1, 1, :, 1))/pi - ref)**2)/nx)
   write (*, '(a,es10.3,a,es10.3,a,f7.4)') 'RMS density error on |x|<0.8, LP x2: Jackson = ', err_j, &
      '  no kernel = ', err_n, '  ratio = ', err_n/err_j
   if (.not. (err_n < 0.25_rp*err_j)) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   pure function poisson(r, phi) result(p)
      real(rp), intent(in) :: r, phi
      real(rp) :: p

      p = (1.0_rp - r*r)/(1.0_rp - 2.0_rp*r*cos(phi) + r*r)
   end function poisson

end program test_cheb_kernel
