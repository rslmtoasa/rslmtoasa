!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
!> @brief Energy nodes and weights of cheb_energy_nodes against direct quadrature of the explicit series.
!> @details G_X(E) = sum_m w_m k_m X_m (-i exp(-i m theta))/sqrt(a^2-(E-b)^2), theta = acos((E-b)/a), w_0 = 1,
!>          w_m = 2. X_m = P_i^H T_m(Htilde) P_j from an eigen-decomposition of a deterministic Hermitian H
!>          (dimension 10, 3-orbital blocks: on-site rows 1-3 and the overlapping pair rows 1-3 / 3-5), n = 40,
!>          arbitrary real kernel. Reference: int_{e_lo}^{inf} f(E) E^p G_X(E) dE by composite Gauss-Legendre
!>          (400 panels x 16 nodes, panel boundary at e_f) of the explicit series, f = step or Fermi function
!>          written here. Test: sum_q w(q) e(q)^p G_X(e(q)) with the nodes of cheb_energy_nodes and the same series.
!>          Also the occupation -(1/pi) Im Tr of the on-site p = 0 result.
!------------------------------------------------------------------------------
program test_cheb_energy_nodes
   use precision_mod, only: rp
   use math_mod, only: cheb_energy_nodes
   implicit none

   integer, parameter :: nd = 10, nb = 3, n = 40, npan = 400, ngl = 16
   real(rp), parameter :: tol = 1.0e-10_rp, pi = 3.14159265358979323846264338327950_rp
   complex(rp), parameter :: ci = (0.0_rp, 1.0_rp)
   complex(rp) :: hm(nd, nd), v(nd, nd), xm(nb, nb, n, 2)
   complex(rp) :: work(64)
   real(rp) :: lam(nd), rwork(3*nd - 2), kern(n), a, b, emin_win, emax_win
   real(rp) :: e_lo, e_f, kbt, xg(ngl), wg(ngl), kbts(2)
   real(rp), allocatable :: en(:), wn(:)
   complex(rp) :: iref(nb, nb), imom(nb, nb)
   real(rp) :: err, nref, nmom, err_n
   integer :: i, j, m, info, is, ik, p, q
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
   e_lo = b - 0.8_rp*a
   e_f = b + 0.137_rp*a
   kbts = [0.0_rp, 0.004_rp]

   do m = 0, n - 1
      kern(m + 1) = 1.0_rp - 0.7_rp*m/n + 0.05_rp*sin(real(m, rp))
      call moment_block(1, 1, m, xm(:, :, m + 1, 1))
      call moment_block(1, 3, m, xm(:, :, m + 1, 2))
   end do
   call gauss_legendre(xg, wg)

   ok = .true.
   write (*, '(a,2f9.4,a,2f9.4)') 'a, b = ', a, b, '   e_lo, e_f = ', e_lo, e_f
   write (*, '(a)') 'set  kbt[Ry]  p   rel.err (Frobenius, complex 3x3)'
   do ik = 1, 2
      kbt = kbts(ik)
      call cheb_energy_nodes(n, a, b, e_lo, e_f, kbt, en, wn)
      do is = 1, 2
         do p = 0, 2
            iref = ref_integral(is, p)
            imom = (0.0_rp, 0.0_rp)
            do q = 1, size(en)
               imom = imom + wn(q)*en(q)**p*green_series(en(q), is)
            end do
            err = sqrt(sum(abs(imom - iref)**2))/sqrt(sum(abs(iref)**2))
            write (*, '(i3,f9.4,i3,es14.4)') is, kbt, p, err
            if (.not. (err < tol)) ok = .false.
         end do
      end do
      ! occupation -(1/pi) Im Tr, on-site p = 0
      iref = ref_integral(1, 0)
      imom = (0.0_rp, 0.0_rp)
      do q = 1, size(en)
         imom = imom + wn(q)*green_series(en(q), 1)
      end do
      nref = 0.0_rp
      nmom = 0.0_rp
      do i = 1, nb
         nref = nref - aimag(iref(i, i))/pi
         nmom = nmom - aimag(imom(i, i))/pi
      end do
      err_n = abs(nmom - nref)/abs(nref)
      write (*, '(a,f8.4,a,2f16.10,es12.3)') 'occupation kbt=', kbt, '  ref, nodes, rel.err:', nref, nmom, err_n
      if (.not. (err_n < tol)) ok = .false.
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

   function fermi(en) result(f)
      real(rp), intent(in) :: en
      real(rp) :: f
      if (kbt > 0.0_rp) then
         f = 1.0_rp/(exp((en - e_f)/kbt) + 1.0_rp)
      else
         f = merge(1.0_rp, 0.0_rp, en <= e_f)
      end if
   end function fermi

   !> int_{e_lo}^{inf} f(E) E^p G(E) dE by composite Gauss-Legendre; upper limit e_f (step) or e_f + 40 kbt.
   function ref_integral(iset, pw) result(res)
      integer, intent(in) :: iset, pw
      complex(rp) :: res(nb, nb)
      real(rp) :: lo, hi, h, x, e1
      integer :: ip, iq, nseg, np1, ipan
      res = (0.0_rp, 0.0_rp)
      if (kbt > 0.0_rp) then
         nseg = 2
      else
         nseg = 1
      end if
      do ip = 1, nseg
         if (nseg == 1) then
            lo = e_lo; hi = e_f; np1 = npan
         else if (ip == 1) then
            lo = e_lo; hi = e_f; np1 = npan/2
         else
            lo = e_f; hi = e_f + 40.0_rp*kbt; np1 = npan/2
         end if
         h = (hi - lo)/np1
         do ipan = 1, np1
            do iq = 1, ngl
               x = lo + (ipan - 1)*h + 0.5_rp*h*(1.0_rp + xg(iq))
               e1 = x
               res = res + 0.5_rp*h*wg(iq)*fermi(e1)*e1**pw*green_series(e1, iset)
            end do
         end do
      end do
   end function ref_integral

   !> Gauss-Legendre nodes and weights on [-1, 1] by Newton iteration.
   subroutine gauss_legendre(x, w)
      real(rp), intent(out) :: x(ngl), w(ngl)
      real(rp) :: z, p1, p2, p3, pp
      integer :: ii, jj, it
      do ii = 1, ngl
         z = cos(pi*(ii - 0.25_rp)/(ngl + 0.5_rp))
         do it = 1, 100
            p1 = 1.0_rp
            p2 = 0.0_rp
            do jj = 1, ngl
               p3 = p2
               p2 = p1
               p1 = ((2*jj - 1)*z*p2 - (jj - 1)*p3)/jj
            end do
            pp = ngl*(z*p1 - p2)/(z*z - 1.0_rp)
            z = z - p1/pp
         end do
         x(ii) = z
         w(ii) = 2.0_rp/((1.0_rp - z*z)*pp*pp)
      end do
   end subroutine gauss_legendre

end program test_cheb_energy_nodes
