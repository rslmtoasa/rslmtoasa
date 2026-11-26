!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
!> @brief Point weights of G_X(E) and d2G_X/dE2 from cheb_green_point_weights against the explicit series.
!> @details G_X(E) = sum_m w_m k_m X_m (-i exp(-i m theta))/sqrt(a^2-(E-b)^2), theta = acos((E-b)/a), w_0 = 1,
!>          w_m = 2. X_m = P_i^H T_m(Htilde) P_j come from an eigen-decomposition of a deterministic Hermitian H
!>          (dimension 10, 3-orbital blocks, rows 1-3 against 3-5, and the transposed pair), n = 40, arbitrary real kernel.
!>          Reference G: the explicit series in this test (the g0 comparison holds by construction, it repeats
!>          the same formula). Reference d2G/dE2: central second differences of that series with Richardson
!>          extrapolation over h, h/2, h/4 (O(h^6)); the spread between the two last extrapolation levels is the
!>          reported reference error and must stay below tol_d2 / 10.
!------------------------------------------------------------------------------
program test_cheb_green_point
   use precision_mod, only: rp
   use math_mod, only: cheb_green_point_weights
   implicit none

   integer, parameter :: nd = 10, nb = 3, n = 40
   real(rp), parameter :: tol_g = 1.0e-12_rp, tol_d2 = 1.0e-7_rp
   complex(rp), parameter :: ci = (0.0_rp, 1.0_rp)
   complex(rp) :: hm(nd, nd), v(nd, nd), xm(nb, nb, n, 2)
   complex(rp) :: work(64), g0(n), g2(n)
   real(rp) :: lam(nd), rwork(3*nd - 2), kern(n), a, b, emin_win, emax_win
   real(rp) :: es(2), e, h, err_g, err_d, spread
   complex(rp) :: gref(nb, nb), gw(nb, nb), d2ref(nb, nb), d2w(nb, nb), r2(nb, nb, 2), r3(nb, nb)
   complex(rp) :: dd(nb, nb, 3)
   integer :: i, j, m, info, ie, is, l
   logical :: ok

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
   es = [b + 0.137_rp*a, b - 0.41_rp*a]

   do m = 0, n - 1
      kern(m + 1) = 1.0_rp - 0.7_rp*m/n + 0.05_rp*sin(real(m, rp))
      call moment_block(1, nb, m, xm(:, :, m + 1, 1))
      call moment_block(nb, 1, m, xm(:, :, m + 1, 2))
   end do

   ok = .true.
   h = 2.0e-3_rp*a
   write (*, '(a,2f9.4,a,es9.2)') 'a, b = ', a, b, '   h = ', h
   write (*, '(a)') 'set  E-b/a    rel.err G (by construction)   rel.err d2G   FD ref spread (rel.)'
   do ie = 1, 2
      e = es(ie)
      call cheb_green_point_weights(kern, a, b, e, g0, g2)
      do is = 1, 2
         gref = green_series(e, is)
         gw = (0.0_rp, 0.0_rp)
         d2w = (0.0_rp, 0.0_rp)
         do m = 1, n
            gw = gw + xm(:, :, m, is)*g0(m)
            d2w = d2w + xm(:, :, m, is)*g2(m)
         end do
         ! second differences at h, h/2, h/4 and Richardson levels
         do l = 1, 3
            dd(:, :, l) = (green_series(e + h/2**(l - 1), is) - 2.0_rp*gref + green_series(e - h/2**(l - 1), is)) &
                          /(h/2**(l - 1))**2
         end do
         r2(:, :, 1) = (4.0_rp*dd(:, :, 2) - dd(:, :, 1))/3.0_rp
         r2(:, :, 2) = (4.0_rp*dd(:, :, 3) - dd(:, :, 2))/3.0_rp
         r3 = (16.0_rp*r2(:, :, 2) - r2(:, :, 1))/15.0_rp
         d2ref = r3
         spread = sqrt(sum(abs(r2(:, :, 2) - r3)**2))/sqrt(sum(abs(d2ref)**2))
         err_g = sqrt(sum(abs(gw - gref)**2))/sqrt(sum(abs(gref)**2))
         err_d = sqrt(sum(abs(d2w - d2ref)**2))/sqrt(sum(abs(d2ref)**2))
         write (*, '(i3,f9.4,3es18.4)') is, (e - b)/a, err_g, err_d, spread
         if (.not. (err_g < tol_g)) ok = .false.
         if (.not. (err_d < tol_d2)) ok = .false.
         if (.not. (spread < tol_d2/10.0_rp)) ok = .false.
      end do
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

   function green_series(en, iset) result(g)
      real(rp), intent(in) :: en
      integer, intent(in) :: iset
      complex(rp) :: g(nb, nb)
      real(rp) :: th, w
      integer :: mm
      th = acos((en - b)/a)
      g = (0.0_rp, 0.0_rp)
      do mm = 0, n - 1
         w = merge(1.0_rp, 2.0_rp, mm == 0)
         g = g + w*kern(mm + 1)*xm(:, :, mm + 1, iset)*(-ci)*exp(-ci*mm*th)
      end do
      g = g/sqrt(a**2 - (en - b)**2)
   end function green_series

end program test_cheb_green_point
