program test_strux_lmax3_pair_blocks
   use strux_support, only: build_common_tables
   use strux_pair_expansion, only: mstrx2
   implicit none

! Bare (kappa = 0) pair blocks S0 at lmax = 3 from the tables that strux_lib
! builds, against the canonical closed form
!   S0_{LL'}(R) = (-1)^l' 8 pi (2(l+l')-1)!! / ((2l-1)!! (2l'-1)!!)
!                 sum_m'' G(L,L',(l+l',m'')) Y_(l+l',m'')(R) / |R|^(l+l'+1)
! with G the real-harmonic Gaunt coefficient.  Here G and Y come from
! Legendre recursion and Gauss quadrature written in this file.
! Bug caught: cy filled only to l = 2 x 2 (25 entries) leaves the factors for
! l'' = 5, 6 uninitialised, so every f-containing block is wrong.

   double precision, parameter :: pi = 3.141592653589793238d0
   integer, parameter :: nv = 4, nt = 24, np = 64
   double precision :: alat, plat(3, 3), pos(3, 1), cy(100), cg(1200), vec(3, nv), s(16, 16), sd(16, 16)
   double precision :: xg(nt), wg(nt), err, errmax, ref, smax
   integer :: indxcg(350), jcg(1200), ips(1), lmxb(1), ntab(2), nttab, iv, l, lp, m, mp, a, b
   integer, allocatable :: iax(:, :)

   alat = 5.4069d0
   plat = reshape([-.5d0, .5d0, .5d0, .5d0, -.5d0, .5d0, .5d0, .5d0, -.5d0], [3, 3])
   pos = 0d0; ips = 1; lmxb = 3
   call build_common_tables(1, 1, alat, plat, pos, 5.45d0, ips, lmxb, nttab, ntab, iax, cy, cg, indxcg, jcg)
   call gauss_legendre()
   vec(:, 1) = [0.5d0, 0.5d0, 0.5d0]; vec(:, 2) = [0.3d0, -0.7d0, 1.1d0]
   vec(:, 3) = [-1.2d0, 0.4d0, 0.2d0]; vec(:, 4) = [0.5d0, 1.5d0, -0.9d0]
   errmax = 0d0; smax = 0d0
   do iv = 1, nv
      call mstrx2(0d0, vec(:, iv), 16, 16, 16, cg, indxcg, jcg, cy, 11, s, sd)
      do l = 0, 3; do lp = 0, 3; do m = -l, l; do mp = -lp, lp
         a = l*l + l + 1 + m; b = lp*lp + lp + 1 + mp
         ref = closed_form(l, m, lp, mp, vec(:, iv))
         err = abs(s(a, b) - ref)/max(1d0, abs(ref))
         errmax = max(errmax, err); smax = max(smax, abs(ref))
      end do; end do; end do; end do
   end do
   write (*, '(a,es10.3,a,es10.3)') 'max relative |S0 - closed form| over 16x16 x 4 vectors: ', errmax, '   max|S0| = ', smax
   if (.not. (errmax < 1d-9)) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   double precision function closed_form(l, m, lp, mp, r)
      integer, intent(in) :: l, m, lp, mp
      double precision, intent(in) :: r(3)
      double precision :: rr, g
      integer :: mm, lt, it, ip
      rr = sqrt(sum(r**2)); lt = l + lp
      closed_form = 0d0
      do mm = -lt, lt
         g = 0d0
         do it = 1, nt; do ip = 1, np
            g = g + wg(it)*(2*pi/np)*yr(l, m, xg(it), 2*pi*(ip - 0.5d0)/np)*yr(lp, mp, xg(it), 2*pi*(ip - 0.5d0)/np) &
                *yr(lt, mm, xg(it), 2*pi*(ip - 0.5d0)/np)
         end do; end do
         closed_form = closed_form + g*yr(lt, mm, r(3)/rr, atan2(r(2), r(1)))
      end do
      closed_form = (-1)**lp*8*pi*dfac(2*lt - 1)/(dfac(2*l - 1)*dfac(2*lp - 1))*closed_form/rr**(lt + 1)
   end function closed_form

   ! Real harmonic without Condon-Shortley phase: m > 0 cos(m phi), m < 0 sin(|m| phi).
   double precision function yr(l, m, x, ph)
      integer, intent(in) :: l, m
      double precision, intent(in) :: x, ph
      integer :: am, i
      double precision :: pmm, pm1, pl, f1, f2
      am = abs(m); pmm = 1d0
      do i = 1, am
         pmm = pmm*(2*i - 1)*sqrt(max(0d0, 1d0 - x*x))
      end do
      pl = pmm
      if (l > am) then
         pm1 = x*(2*am + 1)*pmm; pl = pm1
         do i = am + 2, l
            pl = (x*(2*i - 1)*pm1 - (i + am - 1)*pmm)/(i - am); pmm = pm1; pm1 = pl
         end do
      end if
      f1 = 1d0; f2 = 1d0
      do i = 2, l - am; f1 = f1*i; end do
      do i = 2, l + am; f2 = f2*i; end do
      yr = sqrt((2*l + 1)/(4*pi)*f1/f2)*pl
      if (m > 0) yr = sqrt(2d0)*yr*cos(am*ph)
      if (m < 0) yr = sqrt(2d0)*yr*sin(am*ph)
   end function yr

   double precision function dfac(n)
      integer, intent(in) :: n
      integer :: k
      dfac = 1d0
      do k = n, 1, -2; dfac = dfac*k; end do
   end function dfac

   subroutine gauss_legendre()
      integer :: i, j, it
      double precision :: z, z1, p1, p2, p3, pp
      do i = 1, nt
         z = cos(pi*(i - 0.25d0)/(nt + 0.5d0))
         do it = 1, 100
            p1 = 1d0; p2 = 0d0
            do j = 1, nt
               p3 = p2; p2 = p1; p1 = ((2*j - 1)*z*p2 - (j - 1)*p3)/j
            end do
            pp = nt*(z*p1 - p2)/(z*z - 1d0); z1 = z; z = z1 - p1/pp
            if (abs(z - z1) < 1d-15) exit
         end do
         xg(i) = z; wg(i) = 2d0/((1d0 - z*z)*pp*pp)
      end do
   end subroutine gauss_legendre
end program test_strux_lmax3_pair_blocks
