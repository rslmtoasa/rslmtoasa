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
!>          measured 1.6e-15 (noise-free, same p and window) in numpy.
!>          Hermiticity of the U2 extension holds by construction.
!------------------------------------------------------------------------------
program test_cheb_lp
   use precision_mod, only: rp
   use cheb_lp_mod, only: cheb_lp_extend
   implicit none

   integer, parameter :: n_known = 302, n_total = 4*302
   real(rp), parameter :: tol = 1.0e-11_rp
   real(rp), parameter :: pi = 3.14159265358979323846_rp
   logical :: failed

   failed = .false.
   call test_scalar(failed)
   call test_block(failed)

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
      real(rp) :: lam_max, fit_res, imag_max

      do n = 0, n_total - 1
         ex(n + 1) = 0.0_rp
         do k = 1, nk
            ex(n + 1) = ex(n + 1) + w(k)*rho(k)**n*cos(n*th(k))
         end do
      end do
      mu = (0.0_rp, 0.0_rp)
      mu(1, 1, 1:n_known) = ex(1:n_known)
      call cheb_lp_extend(mu, n_known, n_rank, n_unstable, lam_max, fit_res, imag_max)
      err = maxval(abs(real(mu(1, 1, n_known + 1:)) - ex(n_known + 1:)))/maxval(abs(ex(1:n_known)))
      write (*, '(a,es10.3,a,i0,a,i0,a,f8.5,a,es9.2,a,es9.2)') 'U1 scalar: max rel err = ', err, &
         '  rank=', n_rank, '  unstable=', n_unstable, '  lam_max=', lam_max, '  fit_res=', fit_res, '  imag=', imag_max
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
      real(rp) :: lam_max, fit_res, imag_max

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
      call cheb_lp_extend(mu, n_known, n_rank, n_unstable, lam_max, fit_res, imag_max)
      err = maxval(abs(mu(:, :, n_known + 1:) - ex(:, :, n_known + 1:)))/maxval(abs(ex(:, :, 1:n_known)))
      herm = 0.0_rp
      do n = n_known + 1, n_total
         herm = max(herm, maxval(abs(mu(:, :, n) - conjg(transpose(mu(:, :, n))))))
      end do
      write (*, '(a,es10.3,a,i0,a,i0,a,f8.5,a,es9.2,a,es9.2)') 'U2 3x3 block: max rel err = ', err, &
         '  rank=', n_rank, '  unstable=', n_unstable, '  lam_max=', lam_max, '  fit_res=', fit_res, '  imag=', imag_max
      write (*, '(a,es9.2,a)') 'U2 Hermiticity of extension: max|mu - mu^H| = ', herm, '  (by construction, not counted)'
      if (err > tol) then
         write (*, '(a)') 'FAIL: U2 extension differs from the closed form'
         bad = .true.
      end if
   end subroutine test_block

end program test_cheb_lp
