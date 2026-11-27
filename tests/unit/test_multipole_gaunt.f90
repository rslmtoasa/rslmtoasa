!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
! PROGRAM: test_multipole_gaunt
!
!> @brief Angular part of the general on-site multipoles against numerical quadrature.
!> @details For a random Hermitian density matrix n in the orbital basis of the code (b_lm = (-1)^m Y_lm
!>          Condon-Shortley, m = +l..-l within each shell), rho(r_hat) = sum_ab n_ab b_a conj(b_b) is integrated
!>          against the real harmonics R_LM by Gauss-Legendre x uniform-phi quadrature written here, and compared
!>          with real_gaunt_table + add_multipole_block (unit radial weights) for all L <= 2*lmax, lmax = 2 and 3.
!>          Bug caught: wrong in-shell m order, missing (-1)^m, a and b exchanged, wrong real-harmonic sign.
!------------------------------------------------------------------------------
program test_multipole_gaunt
   use electrostatics_multipole_mod, only: real_gaunt_table, add_multipole_block
   use precision_mod, only: rp
   implicit none

   real(rp), parameter :: pi = 3.14159265358979323846_rp
   real(rp), parameter :: tol = 1.0e-12_rp
   logical :: failed

   failed = .false.
   call check_lmax(2)
   call check_lmax(3)
   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   subroutine check_lmax(lmax)
      integer, intent(in) :: lmax
      integer, parameter :: nth = 16, nph = 40
      integer :: nsh, ia, ib, l, m, it, ip, ll, mm
      real(rp) :: x(nth), wx(nth), phi, h, q(0:2*lmax, -2*lmax:2*lmax), qq(0:2*lmax, -2*lmax:2*lmax), diff, scale
      real(rp), allocatable :: w(:, :, :), re(:, :), im(:, :)
      complex(rp), allocatable :: n(:, :), a(:, :, :, :), b(:)
      complex(rp) :: rho

      nsh = (lmax + 1)**2
      allocate (n(nsh, nsh), b(nsh), re(nsh, nsh), im(nsh, nsh), w(0:lmax, 0:lmax, 0:2*lmax))
      call random_number(re)
      call random_number(im)
      n = cmplx(re - 0.5_rp, im - 0.5_rp, rp)
      n = 0.5_rp*(n + conjg(transpose(n)))
      w = 1.0_rp

      qq = 0.0_rp
      call gauss_legendre(nth, x, wx)
      h = 2.0_rp*pi/nph
      do it = 1, nth
         do ip = 1, nph
            phi = (ip - 1)*h
            do l = 0, lmax
               do m = l, -l, -1
                  b(l*l + (l - m) + 1) = (-1)**m*ylm(l, m, x(it), phi)
               end do
            end do
            rho = (0.0_rp, 0.0_rp)
            do ib = 1, nsh
               do ia = 1, nsh
                  rho = rho + n(ia, ib)*b(ia)*conjg(b(ib))
               end do
            end do
            if (abs(aimag(rho)) > 1.0e-13_rp) call fail('rho is not real')
            do ll = 0, 2*lmax
               do mm = -ll, ll
                  qq(ll, mm) = qq(ll, mm) + wx(it)*h*real(rho, rp)*rlm(ll, mm, x(it), phi)
               end do
            end do
         end do
      end do

      call real_gaunt_table(lmax, a)
      q = 0.0_rp
      call add_multipole_block(lmax, a, n, w, q)
      scale = maxval(abs(qq))
      diff = maxval(abs(q - qq))
      write (*, '(a,i0,a,es10.3,a,es10.3)') 'lmax=', lmax, ' max|q_code - q_quad| = ', diff, ' max|q_quad| = ', scale
      if (diff > tol*scale) call fail('angular contraction differs from quadrature')
   end subroutine check_lmax

   !> Complex Y_lm(theta, phi), Condon-Shortley phase, x = cos(theta), any sign of m.
   function ylm(l, m, x, phi) result(y)
      integer, intent(in) :: l, m
      real(rp), intent(in) :: x, phi
      complex(rp) :: y
      integer :: ma, k
      real(rp) :: pmm, pmk, pmk1, s, nrm

      ma = abs(m)
      pmm = 1.0_rp
      s = sqrt((1.0_rp - x)*(1.0_rp + x))
      do k = 1, ma
         pmm = -pmm*(2*k - 1)*s
      end do
      pmk = pmm
      if (l > ma) then
         pmk1 = x*(2*ma + 1)*pmm
         pmk = pmk1
         do k = ma + 2, l
            pmk1 = (x*(2*k - 1)*pmk - (k + ma - 1)*pmm)/(k - ma)
            pmm = pmk
            pmk = pmk1
         end do
      end if
      nrm = sqrt((2*l + 1)/(4.0_rp*pi)*exp(log_gamma(real(l - ma + 1, rp)) - log_gamma(real(l + ma + 1, rp))))
      y = nrm*pmk*cmplx(cos(ma*phi), sin(ma*phi), rp)
      if (m < 0) y = (-1)**ma*conjg(y)
   end function ylm

   !> Real harmonic R_LM: Y_L0, sqrt2 (-1)^M Re Y_LM (M > 0), sqrt2 (-1)^M Im Y_L|M| (M < 0).
   function rlm(l, m, x, phi) result(r)
      integer, intent(in) :: l, m
      real(rp), intent(in) :: x, phi
      real(rp) :: r

      if (m == 0) then
         r = real(ylm(l, 0, x, phi), rp)
      else if (m > 0) then
         r = sqrt(2.0_rp)*(-1)**m*real(ylm(l, m, x, phi), rp)
      else
         r = sqrt(2.0_rp)*(-1)**m*aimag(ylm(l, -m, x, phi))
      end if
   end function rlm

   !> Gauss-Legendre nodes and weights on [-1, 1] by Newton iteration on P_n.
   subroutine gauss_legendre(nn, x, w)
      integer, intent(in) :: nn
      real(rp), intent(out) :: x(nn), w(nn)
      integer :: i, j, it
      real(rp) :: z, p1, p2, p3, pp

      do i = 1, nn
         z = cos(pi*(i - 0.25_rp)/(nn + 0.5_rp))
         do it = 1, 100
            p1 = 1.0_rp
            p2 = 0.0_rp
            do j = 1, nn
               p3 = p2
               p2 = p1
               p1 = ((2*j - 1)*z*p2 - (j - 1)*p3)/j
            end do
            pp = nn*(z*p1 - p2)/(z*z - 1.0_rp)
            z = z - p1/pp
         end do
         x(i) = z
         w(i) = 2.0_rp/((1.0_rp - z*z)*pp*pp)
      end do
   end subroutine gauss_legendre

   subroutine fail(msg)
      character(len=*), intent(in) :: msg

      write (*, '(a)') 'FAIL: '//msg
      failed = .true.
   end subroutine fail

end program test_multipole_gaunt
