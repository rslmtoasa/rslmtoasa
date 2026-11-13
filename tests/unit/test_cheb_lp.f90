!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
!> @brief Linear prediction of Chebyshev moments against closed-form sequences.
!> @details U1: scalar damped-pole sequence x(n) = sum_k w_k rho_k^n cos(n theta_k).
!>          U2: 3x3 complex block mu(n) = sum_k W_k rho_k^n cos(n theta_k) with
!>          Hermitian positive-semidefinite W_k. The closed forms are written
!>          in this test; the extension is compared on n = N..N_ext-1 against
!>          max|x| of the known part. Tolerance 1e-11: the Stage A toy study
!>          measured 3.6e-14 (noise-free, same p and window, N/2 and N/2) in numpy.
!>          Hermiticity of the U2 extension holds by construction.
!>          U3: the extended moments go through the production reconstruction
!>          cheb_green_dp. The closed-form density is the Poisson-kernel sum
!>          A(x) = sum_k w_k [P(t-t_k) + P(t+t_k)]/(2 pi sin t), x = cos t,
!>          P(f) = (1-r^2)/(1-2 r cos f+r^2). The RMS error on |x|<0.8 must be
!>          below 0.6 times the Jackson@N error (worst Stage A ratio: 0.51).
!>          U4: one damped and one growing pole (rho = 1.003). The growing pair must be
!>          dropped and carry the share w2 rho2^(N-1)/(w1 rho1^(N-1) + w2 rho2^(N-1))
!>          of the mode amplitude at the last known moment.
!>          U5: as U2 but with NON-Hermitian W_k = u v^H (u /= v), so the full-channel path is taken;
!>          the closed form is compared to the same 1e-11 as U2.
!------------------------------------------------------------------------------
program test_cheb_lp
   use precision_mod, only: rp
   use cheb_lp_mod, only: cheb_lp_extend
   use chebyshev_fast_mod, only: cheb_green_dp
   implicit none

   integer, parameter :: n_known = 302, n_total = 4*302
   real(rp), parameter :: tol = 1.0e-11_rp
   real(rp), parameter :: pi = 3.14159265358979323846_rp
   logical :: failed

   failed = .false.
   call test_scalar(failed)
   call test_block(failed)
   call test_spectral(failed)
   call test_unstable(failed)
   call test_nonhermitian(failed)

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   !> U1: catches a wrong recurrence index/sign or companion-matrix layout.
   subroutine test_scalar(bad)
      logical, intent(inout) :: bad
      integer, parameter :: nk = 6
      real(rp), parameter :: rho(nk) = [0.975_rp, 0.982_rp, 0.988_rp, 0.991_rp, 0.994_rp, 0.997_rp]
      real(rp), parameter :: th(nk) = [0.35_rp, 0.80_rp, 1.25_rp, 1.70_rp, 2.20_rp, 2.75_rp]
      real(rp), parameter :: w(nk) = [1.0_rp, 0.7_rp, 1.3_rp, 0.5_rp, 0.9_rp, 0.6_rp]
      complex(rp) :: mu(1, 1, n_total)
      real(rp) :: ex(n_total), err
      integer :: n, k, n_rank, n_unstable
      real(rp) :: lam_max, fit_res, imag_max, uw

      do n = 0, n_total - 1
         ex(n + 1) = 0.0_rp
         do k = 1, nk
            ex(n + 1) = ex(n + 1) + w(k)*rho(k)**n*cos(n*th(k))
         end do
      end do
      mu = (0.0_rp, 0.0_rp)
      mu(1, 1, 1:n_known) = ex(1:n_known)
      call cheb_lp_extend(mu, n_known, n_rank, n_unstable, lam_max, fit_res, imag_max, uw)
      err = maxval(abs(real(mu(1, 1, n_known + 1:)) - ex(n_known + 1:)))/maxval(abs(ex(1:n_known)))
      write (*, '(a,es10.3,a,i0,a,i0,a,f8.5,a,es9.2,a,es9.2,a,es9.2)') 'U1 scalar: max rel err = ', err, &
         '  rank=', n_rank, '  unstable=', n_unstable, '  lam_max=', lam_max, '  fit_res=', fit_res, '  imag=', imag_max, '  unst_weight=', uw
      if (err > tol) then
         write (*, '(a)') 'FAIL: U1 extrapolation differs from the closed form'
         bad = .true.
      end if
   end subroutine test_scalar

   !> U2: catches channel or real/imag mixing and per-element coefficients.
   subroutine test_block(bad)
      logical, intent(inout) :: bad
      integer, parameter :: nk = 4, nb = 3
      real(rp), parameter :: rho(nk) = [0.970_rp, 0.985_rp, 0.992_rp, 0.996_rp]
      real(rp), parameter :: th(nk) = [0.45_rp, 1.05_rp, 1.90_rp, 2.60_rp]
      complex(rp) :: wk(nb, nb, nk), mu(nb, nb, n_total), ex(nb, nb, n_total), u(nb), v(nb)
      real(rp) :: err, herm
      integer :: n, k, l, n_rank, n_unstable
      real(rp) :: lam_max, fit_res, imag_max, uw

      do k = 1, nk
         do l = 1, nb
            u(l) = cmplx(cos(0.7_rp*l + 1.3_rp*k), sin(1.1_rp*l - 0.4_rp*k), rp)
            v(l) = cmplx(sin(0.9_rp*l + 0.6_rp*k), cos(0.5_rp*l - 1.7_rp*k), rp)
         end do
         ! Hermitian PSD: u u^H + 0.3 v v^H
         wk(:, :, k) = spread(u, 2, nb)*conjg(spread(u, 1, nb)) + 0.3_rp*spread(v, 2, nb)*conjg(spread(v, 1, nb))
      end do
      do n = 0, n_total - 1
         ex(:, :, n + 1) = (0.0_rp, 0.0_rp)
         do k = 1, nk
            ex(:, :, n + 1) = ex(:, :, n + 1) + wk(:, :, k)*rho(k)**n*cos(n*th(k))
         end do
      end do
      mu = (0.0_rp, 0.0_rp)
      mu(:, :, 1:n_known) = ex(:, :, 1:n_known)
      call cheb_lp_extend(mu, n_known, n_rank, n_unstable, lam_max, fit_res, imag_max, uw)
      err = maxval(abs(mu(:, :, n_known + 1:) - ex(:, :, n_known + 1:)))/maxval(abs(ex(:, :, 1:n_known)))
      herm = 0.0_rp
      do n = n_known + 1, n_total
         herm = max(herm, maxval(abs(mu(:, :, n) - conjg(transpose(mu(:, :, n))))))
      end do
      write (*, '(a,es10.3,a,i0,a,i0,a,f8.5,a,es9.2,a,es9.2,a,es9.2)') 'U2 3x3 block: max rel err = ', err, &
         '  rank=', n_rank, '  unstable=', n_unstable, '  lam_max=', lam_max, '  fit_res=', fit_res, '  imag=', imag_max, '  unst_weight=', uw
      write (*, '(a,es9.2,a)') 'U2 Hermiticity of extension: max|mu - mu^H| = ', herm, '  (by construction, not counted)'
      if (err > tol) then
         write (*, '(a)') 'FAIL: U2 extension differs from the closed form'
         bad = .true.
      end if
   end subroutine test_block

   !> U3: catches a kernel applied at the wrong n_mom and an extension that does not sharpen the density.
   subroutine test_spectral(bad)
      logical, intent(inout) :: bad
      integer, parameter :: nk = 3, nx = 401, n_ext = n_total
      real(rp), parameter :: rho(nk) = [0.970_rp, 0.985_rp, 0.990_rp]
      real(rp), parameter :: th(nk) = [0.90_rp, 1.40_rp, 2.10_rp]
      real(rp), parameter :: w(nk) = [1.0_rp, 0.8_rp, 1.2_rp]
      complex(rp) :: mu(1, 1, n_ext), mu_n(1, 1, n_known, 1), mu_e(1, 1, n_ext, 1)
      complex(rp) :: g_n(1, 1, nx, 1), g_e(1, 1, nx, 1)
      real(rp) :: x(nx), ref(nx), t, err_n, err_e, lam_max, fit_res, imag_max, uw
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
            mu(1, 1, n + 1) = mu(1, 1, n + 1) + w(k)*rho(k)**n*cos(n*th(k))
         end do
      end do
      mu_n(1, 1, :, 1) = mu(1, 1, 1:n_known)
      call cheb_lp_extend(mu, n_known, n_rank, n_unstable, lam_max, fit_res, imag_max, uw)
      mu_e(:, :, :, 1) = mu
      call cheb_green_dp(mu_n, 1, n_known, 1, x, nx, 1.0_rp, 0.0_rp, g_n)
      call cheb_green_dp(mu_e, 1, n_ext, 1, x, nx, 1.0_rp, 0.0_rp, g_e)
      err_n = sqrt(sum((-aimag(g_n(1, 1, :, 1))/pi - ref)**2)/nx)
      err_e = sqrt(sum((-aimag(g_e(1, 1, :, 1))/pi - ref)**2)/nx)
      write (*, '(a,es10.3,a,es10.3,a,f6.3)') 'U3 density RMS error on |x|<0.8: Jackson@302 = ', err_n, &
         '  LP->1208 = ', err_e, '  ratio = ', err_e/err_n
      if (.not. (err_e < 0.6_rp*err_n)) then
         write (*, '(a)') 'FAIL: U3 LP density is not below 0.6 x the Jackson@N error'
         bad = .true.
      end if
   end subroutine test_spectral

   !> U4: catches a wrong unstable-root mask or a weight taken after the roots are zeroed.
   subroutine test_unstable(bad)
      logical, intent(inout) :: bad
      real(rp), parameter :: rho(2) = [0.985_rp, 1.003_rp], th(2) = [0.9_rp, 2.0_rp], w(2) = [1.0_rp, 0.05_rp]
      complex(rp) :: mu(1, 1, n_total)
      real(rp) :: lam_max, fit_res, imag_max, uw, expected, a1, a2
      integer :: n, k, n_rank, n_unstable

      mu = (0.0_rp, 0.0_rp)
      do n = 0, n_known - 1
         do k = 1, 2
            mu(1, 1, n + 1) = mu(1, 1, n + 1) + w(k)*rho(k)**n*cos(n*th(k))
         end do
      end do
      call cheb_lp_extend(mu, n_known, n_rank, n_unstable, lam_max, fit_res, imag_max, uw)
      a1 = w(1)*rho(1)**(n_known - 1)
      a2 = w(2)*rho(2)**(n_known - 1)
      expected = a2/(a1 + a2)
      write (*, '(a,i0,a,i0,a,f8.5,a,es10.3,a,es10.3)') 'U4 growing pole: unstable=', n_unstable, '  rank=', n_rank, &
         '  lam_max=', lam_max, '  weight=', uw, '  closed form=', expected
      if (n_unstable /= 2 .or. abs(uw - expected) > 1.0e-7_rp*expected) then
         write (*, '(a)') 'FAIL: U4 dropped roots or their weight differ from the closed form'
         bad = .true.
      end if
   end subroutine test_unstable

   !> U5: catches a silent Hermitian assumption (symmetrised input or mirrored output) and a wrong full-channel fallback.
   subroutine test_nonhermitian(bad)
      logical, intent(inout) :: bad
      integer, parameter :: nk = 4, nb = 3
      real(rp), parameter :: rho(nk) = [0.970_rp, 0.985_rp, 0.992_rp, 0.996_rp]
      real(rp), parameter :: th(nk) = [0.45_rp, 1.05_rp, 1.90_rp, 2.60_rp]
      complex(rp) :: wk(nb, nb, nk), mu(nb, nb, n_total), ex(nb, nb, n_total), u(nb), v(nb)
      real(rp) :: err, herm
      integer :: n, k, l, n_rank, n_unstable
      real(rp) :: lam_max, fit_res, imag_max, uw

      do k = 1, nk
         do l = 1, nb
            u(l) = cmplx(cos(0.7_rp*l + 1.3_rp*k), sin(1.1_rp*l - 0.4_rp*k), rp)
            v(l) = cmplx(sin(0.9_rp*l + 0.6_rp*k), cos(0.5_rp*l - 1.7_rp*k), rp)
         end do
         wk(:, :, k) = spread(u, 2, nb)*conjg(spread(v, 1, nb))
      end do
      do n = 0, n_total - 1
         ex(:, :, n + 1) = (0.0_rp, 0.0_rp)
         do k = 1, nk
            ex(:, :, n + 1) = ex(:, :, n + 1) + wk(:, :, k)*rho(k)**n*cos(n*th(k))
         end do
      end do
      herm = 0.0_rp
      do n = 1, n_known
         herm = max(herm, maxval(abs(ex(:, :, n) - conjg(transpose(ex(:, :, n))))))
      end do
      mu = (0.0_rp, 0.0_rp)
      mu(:, :, 1:n_known) = ex(:, :, 1:n_known)
      call cheb_lp_extend(mu, n_known, n_rank, n_unstable, lam_max, fit_res, imag_max, uw)
      err = maxval(abs(mu(:, :, n_known + 1:) - ex(:, :, n_known + 1:)))/maxval(abs(ex(:, :, 1:n_known)))
      write (*, '(a,es10.3,a,i0,a,i0,a,f8.5,a,es9.2,a,es9.2,a,es9.2)') 'U5 non-Hermitian block: max rel err = ', err, &
         '  rank=', n_rank, '  unstable=', n_unstable, '  lam_max=', lam_max, '  fit_res=', fit_res, '  imag=', imag_max, '  unst_weight=', uw
      write (*, '(a,es9.2,a)') 'U5 input non-Hermiticity max|mu - mu^H| = ', herm, '  (must be well above 1e-8 max|mu|)'
      if (herm < 1.0e-3_rp*maxval(abs(ex(:, :, 1:n_known)))) then
         write (*, '(a)') 'FAIL: U5 input is not non-Hermitian'
         bad = .true.
      end if
      if (err > tol) then
         write (*, '(a)') 'FAIL: U5 extension differs from the closed form'
         bad = .true.
      end if
   end subroutine test_nonhermitian

   pure function poisson(r, phi) result(p)
      real(rp), intent(in) :: r, phi
      real(rp) :: p

      p = (1.0_rp - r*r)/(1.0_rp - 2.0_rp*r*cos(phi) + r*r)
   end function poisson

end program test_cheb_lp
