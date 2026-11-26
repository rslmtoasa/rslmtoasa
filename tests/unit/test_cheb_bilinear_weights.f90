!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
!> @brief Moment-space energy integral of f(E) E^p Tr[A G_X(E) B G_Y(E)] against direct quadrature.
!> @details G_X(E) = sum_m w_m k_m X_m (-i exp(-i m theta))/sqrt(a^2-(E-b)^2), theta = acos((E-b)/a), w_0 = 1,
!>          w_m = 2. X_m = P_i^H T_m(Htilde) P_j and Y_m = P_j^H T_m(Htilde) P_i come from an eigen-decomposition
!>          of a random Hermitian H (dimension 10, two 3-orbital blocks, rows 1-3 and 3-5, sharing row 3 so that X_0 and Y_0 are non-zero), with an arbitrary real kernel k_m.
!>          A(E) = alpha_i + beta_i E and B(E) = alpha_j + beta_j E are diagonal. The reference integrates
!>          the explicit series with composite Gauss-Legendre (400 panels x 16 nodes); the moment side is
!>          sum_m Tr[A X_m B Ytilde_m] built from cheb_bilinear_weights / cheb_bilinear_transform.
!>          Cases: step (kbt = 0) and Fermi (kbt = 0.004), with alpha-alpha (p = 0), alpha-beta (p = 1),
!>          beta-beta (p = 2) and the full A, B. A final check scales beta_i by 1.01 in the reference only
!>          and requires the relative difference to exceed 1e-4, so the comparison resolves a p-column error.
!------------------------------------------------------------------------------
program test_cheb_bilinear_weights
   use precision_mod, only: rp
   use math_mod, only: gauss_legendre, cheb_bilinear_weights, cheb_bilinear_transform
   implicit none

   integer, parameter :: nd = 10, nb = 3, n = 40, npanel = 400, nq = 16
   real(rp), parameter :: pi = 3.14159265358979323846_rp
   real(rp), parameter :: tol = 1.0e-10_rp, tol_sens = 1.0e-4_rp
   complex(rp), parameter :: ci = (0.0_rp, 1.0_rp)
   real(rp), parameter :: ai0(nb) = [0.30_rp, -0.55_rp, 0.80_rp], bi0(nb) = [0.20_rp, 0.45_rp, -0.35_rp]
   real(rp), parameter :: aj0(nb) = [-0.40_rp, 0.65_rp, 0.15_rp], bj0(nb) = [0.50_rp, -0.25_rp, 0.60_rp]
   real(rp), parameter :: zero(nb) = 0.0_rp
   complex(rp) :: hm(nd, nd), v(nd, nd), xm(nb, nb, n), ym(nb, nb, n), yt(nb, nb, n, 0:2)
   complex(rp) :: hh(0:2*n - 2, 0:2)
   real(rp) :: lam(nd), rwork(3*nd - 2), kern(n), a, b, emin_win, emax_win, e_lo, e_f
   real(rp) :: kbt, err, err_sens
   complex(rp) :: work(64), xref, xmom, xref_s
   integer :: i, j, m, k, info, icase, p
   logical :: ok
   character(len=12) :: lab(4)

   ! deterministic Hermitian H
   do j = 1, nd
      do i = 1, nd
         hm(i, j) = cmplx(cos(0.7_rp*i*j + 0.3_rp*i), sin(0.5_rp*i + 1.3_rp*j*j), rp)
      end do
   end do
   hm = 0.5_rp*(hm + conjg(transpose(hm)))
   v = hm
   call zheev('V', 'U', nd, v, nd, lam, work, 64, rwork, info)
   if (info /= 0) error stop 'zheev failed'
   emin_win = lam(1) - 0.1_rp
   emax_win = lam(nd) + 0.1_rp
   a = (emax_win - emin_win)/(2.0_rp - 0.3_rp)
   b = 0.5_rp*(emax_win + emin_win)
   e_lo = b - 0.8_rp*a
   e_f = b + 0.137_rp*a

   do m = 0, n - 1
      kern(m + 1) = 1.0_rp - 0.7_rp*m/n + 0.05_rp*sin(real(m, rp))
      call moment_block(1, nb, m, xm(:, :, m + 1))
      call moment_block(nb, 1, m, ym(:, :, m + 1))
   end do

   lab = [character(len=12) :: 'alpha-alpha', 'alpha-beta', 'beta-beta', 'full']
   ok = .true.
   write (*, '(a,2f9.4,a,f9.4)') 'a, b = ', a, b, '   e_f - b = ', e_f - b
   write (*, '(a)') 'kbt [Ry-like]  case          |X_ref|        rel. error |X_ref - X_mom|/|X_ref|'
   do k = 1, 2
      kbt = merge(0.0_rp, 0.004_rp, k == 1)
      call cheb_bilinear_weights(n, a, b, e_lo, e_f, kbt, hh)
      do p = 0, 2
         call cheb_bilinear_transform(hh(:, p), kern, ym, yt(:, :, :, p))
      end do
      do icase = 1, 4
         select case (icase)
         case (1); call compare(ai0, zero, aj0, zero, xref, xmom)
         case (2); call compare(ai0, zero, zero, bj0, xref, xmom)
         case (3); call compare(zero, bi0, zero, bj0, xref, xmom)
         case (4); call compare(ai0, bi0, aj0, bj0, xref, xmom)
         end select
         err = abs(xref - xmom)/abs(xref)
         write (*, '(f9.4,5x,a,es12.4,es16.6)') kbt, lab(icase), abs(xref), err
         if (.not. (err < tol)) ok = .false.
      end do
      ! sensitivity: beta_i scaled by 1.01 in the reference only
      call compare(ai0, bi0, aj0, bj0, xref, xmom, 1.01_rp, xref_s)
      err_sens = abs(xref_s - xmom)/abs(xref_s)
      write (*, '(f9.4,5x,a,es12.4,es16.6)') kbt, 'beta_i*1.01  ', abs(xref_s), err_sens
      if (.not. (err_sens > tol_sens)) ok = .false.
   end do
   if (.not. ok) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   !> Block P_r^H T_m(Htilde) P_c with P_r = basis rows r0..r0+nb-1 and P_c = rows c0..c0+nb-1.
   subroutine moment_block(r0, c0, mm, blk)
      integer, intent(in) :: r0, c0, mm
      complex(rp), intent(out) :: blk(nb, nb)
      integer :: kk, ii, jj
      blk = (0.0_rp, 0.0_rp)
      do kk = 1, nd
         do jj = 1, nb
            do ii = 1, nb
               blk(ii, jj) = blk(ii, jj) + conjg(v(r0 + ii - 1, kk))*cos(mm*acos((lam(kk) - b)/a))*v(c0 + jj - 1, kk)
            end do
         end do
      end do
   end subroutine moment_block

   !> Reference (explicit series + quadrature) and moment-side value for diagonal A, B.
   subroutine compare(ai, bi, aj, bj, xr, xo, scale_bi, xr_s)
      real(rp), intent(in) :: ai(nb), bi(nb), aj(nb), bj(nb)
      complex(rp), intent(out) :: xr, xo
      real(rp), intent(in), optional :: scale_bi
      complex(rp), intent(out), optional :: xr_s
      real(rp) :: upper
      complex(rp) :: acc
      integer :: mm
      upper = e_f
      if (kbt > 0.0_rp) upper = e_f + 40.0_rp*kbt
      if (upper >= b + a) error stop 'upper limit outside Chebyshev window'
      xr = (0.0_rp, 0.0_rp)
      call integrate(e_lo, e_f, ai, bi, aj, bj, acc)
      xr = xr + acc
      if (kbt > 0.0_rp) then
         call integrate(e_f, upper, ai, bi, aj, bj, acc)
         xr = xr + acc
      end if
      if (present(scale_bi)) then
         xr_s = (0.0_rp, 0.0_rp)
         call integrate(e_lo, e_f, ai, scale_bi*bi, aj, bj, acc)
         xr_s = xr_s + acc
         if (kbt > 0.0_rp) then
            call integrate(e_f, upper, ai, scale_bi*bi, aj, bj, acc)
            xr_s = xr_s + acc
         end if
      end if
      xo = (0.0_rp, 0.0_rp)
      do mm = 1, n
         xo = xo + tr(dg(ai), xm(:, :, mm), dg(aj), yt(:, :, mm, 0)) &
              + tr(dg(bi), xm(:, :, mm), dg(aj), yt(:, :, mm, 1)) &
              + tr(dg(ai), xm(:, :, mm), dg(bj), yt(:, :, mm, 1)) &
              + tr(dg(bi), xm(:, :, mm), dg(bj), yt(:, :, mm, 2))
      end do
   end subroutine compare

   !> Composite Gauss-Legendre of f(E) Tr[A(E) G_X(E) B(E) G_Y(E)] on [lo, hi].
   subroutine integrate(lo, hi, ai, bi, aj, bj, res)
      real(rp), intent(in) :: lo, hi, ai(nb), bi(nb), aj(nb), bj(nb)
      complex(rp), intent(out) :: res
      real(rp) :: xq(nq), wq(nq), p0, p1, e, fe
      integer :: ip, iq
      complex(rp) :: gx(nb, nb), gy(nb, nb)
      call gauss_legendre(nq, -1.0_rp, 1.0_rp, xq, wq)
      res = (0.0_rp, 0.0_rp)
      do ip = 1, npanel
         p0 = lo + (hi - lo)*(ip - 1)/npanel
         p1 = lo + (hi - lo)*ip/npanel
         do iq = 1, nq
            e = 0.5_rp*(p1 - p0)*xq(iq) + 0.5_rp*(p1 + p0)
            if (kbt > 0.0_rp) then
               fe = 1.0_rp/(exp((e - e_f)/kbt) + 1.0_rp)
            else
               fe = 1.0_rp
            end if
            call green_series(e, xm, gx)
            call green_series(e, ym, gy)
            res = res + 0.5_rp*(p1 - p0)*wq(iq)*fe &
                  *tr(dg(ai + bi*e), gx, dg(aj + bj*e), gy)
         end do
      end do
   end subroutine integrate

   subroutine green_series(e, mom, g)
      real(rp), intent(in) :: e
      complex(rp), intent(in) :: mom(nb, nb, n)
      complex(rp), intent(out) :: g(nb, nb)
      real(rp) :: th, w
      integer :: mm
      th = acos((e - b)/a)
      g = (0.0_rp, 0.0_rp)
      do mm = 0, n - 1
         w = merge(1.0_rp, 2.0_rp, mm == 0)
         g = g + w*kern(mm + 1)*mom(:, :, mm + 1)*(-ci)*exp(-ci*mm*th)
      end do
      g = g/sqrt(a**2 - (e - b)**2)
   end subroutine green_series

   pure function dg(d) result(mat)
      real(rp), intent(in) :: d(nb)
      complex(rp) :: mat(nb, nb)
      integer :: ii
      mat = (0.0_rp, 0.0_rp)
      do ii = 1, nb
         mat(ii, ii) = d(ii)
      end do
   end function dg

   pure function tr(p, q, r, s) result(t)
      complex(rp), intent(in) :: p(nb, nb), q(nb, nb), r(nb, nb), s(nb, nb)
      complex(rp) :: t, prod(nb, nb)
      integer :: ii
      prod = matmul(matmul(p, q), matmul(r, s))
      t = (0.0_rp, 0.0_rp)
      do ii = 1, nb
         t = t + prod(ii, ii)
      end do
   end function tr

end program test_cheb_bilinear_weights
