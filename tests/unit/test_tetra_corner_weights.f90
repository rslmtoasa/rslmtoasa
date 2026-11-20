!> Oracles for tetra_corner_weights (per-corner linear tetrahedron DOS weights and the Bloechl curvature correction).
!> Expected values come from formulas of P. E. Bloechl, O. Jepsen, O. K. Andersen, PRB 49, 16223 (1994) coded in this file
!> (A1)-(A5), (B1)-(B19), (C1)-(C4), (22), and from independent quadratures, never from the production routine:
!>  (a) single tetrahedron: sum of corner weights = D_T (C); sum of B weights = n (A); int^E w_i dE = (B) weight; int w_i E^m
!>      dE = <lambda_i eps^m> by a collapsed-coordinate Gauss rule on the simplex (m = 0 gives 1/4 per corner);
!>  (b) int^E dw_i = (1/40) D_T(E) sum_j (e_j - e_i), Eq. (22), and sum_i dw_i = 0;
!>  (c) shifted free-electron band eps = |k-k0|^2 on cubic meshes of [-1,1]^3 (6 tetrahedra per cell): at the Fermi energy
!>      where the linear count equals the ball volume (4 pi/3) R^3, the corner-weighted band energy int eps theta
!>      (exact 4 pi R^5/5) and an f2 integral (closed form via the ball moment int u_x^2 = 4 pi R^5/15) are compared for
!>      linear and corrected weights. The correction leaves the count unchanged (sum of corner corrections is zero).
program test_tetra_corner_weights
   use precision_mod, only: rp
   use reciprocal_mod, only: tetra_corner_weights
   implicit none

   real(rp), parameter :: pi = 3.14159265358979323846_rp
   integer, parameter :: ngl = 12
   real(rp) :: gx(ngl), gw(ngl)
   real(rp) :: etet(4, 4)
   integer :: it, m, i, ie
   real(rp) :: worst_a, worst_b, worst_sum

   call gauss_legendre(ngl, gx, gw)
   etet(:, 1) = [0.00_rp, 0.31_rp, 0.52_rp, 1.10_rp]
   etet(:, 2) = [0.70_rp, -0.20_rp, 0.15_rp, 0.40_rp]    ! unsorted corners
   etet(:, 3) = [0.10_rp, 0.10_rp, 0.60_rp, 0.90_rp]     ! e1 = e2
   etet(:, 4) = [0.00_rp, 0.25_rp, 0.80_rp, 0.80_rp]     ! e3 = e4

   worst_a = 0.0_rp
   worst_b = 0.0_rp
   worst_sum = 0.0_rp
   do it = 1, 4
      call check_sum_dos(etet(:, it), worst_sum)
      do m = 0, 3
         call check_moments(etet(:, it), m, worst_a)
      end do
      do ie = 1, 7
         call check_corr_integral(etet(:, it), minval(etet(:, it)) + real(ie, rp)*(maxval(etet(:, it)) - minval(etet(:, it)))/7.0_rp, &
                                  worst_b)
      end do
   end do
   print '(a,3es11.3)', 'worst |sum w - D_T (C)| or |sum B - n (A)|, |int w E^m - <lambda E^m>|, |int w - (B)| or |int dw - (22)| :', worst_sum, worst_a, worst_b
   if (.not. worst_sum <= 1.0e-12_rp) error stop 'sum of corner weights differs from tetrahedron DOS'
   if (.not. worst_a <= 1.0e-12_rp) error stop 'corner-weight moments differ from simplex quadrature'
   if (.not. worst_b <= 1.0e-12_rp) error stop 'integrated weights differ from (B) or (22)'

   call check_mesh()
   print '(a)', 'test_tetra_corner_weights: passed'

contains

   !> Tetrahedron DOS (C1)-(C4) per unit volume, sorted corner energies es.
   function dos_closed(es, e) result(d)
      real(rp), intent(in) :: es(4), e
      real(rp) :: d, e21, e31, e41, e32, e42, e43
      e21 = es(2) - es(1); e31 = es(3) - es(1); e41 = es(4) - es(1)
      e32 = es(3) - es(2); e42 = es(4) - es(2); e43 = es(4) - es(3)
      d = 0.0_rp
      if (e <= es(1) .or. e >= es(4)) return
      if (e < es(2)) then
         d = 3.0_rp*(e - es(1))**2/(e21*e31*e41)
      else if (e < es(3)) then
         d = (3.0_rp*e21 + 6.0_rp*(e - es(2)) - 3.0_rp*(e31 + e42)*(e - es(2))**2/(e32*e42))/(e31*e41)
      else
         d = 3.0_rp*(es(4) - e)**2/(e41*e42*e43)
      end if
   end function dos_closed

   !> Number of states (A1)-(A5) per unit volume, sorted corner energies es.
   function nos_closed(es, e) result(n)
      real(rp), intent(in) :: es(4), e
      real(rp) :: n, e21, e31, e41, e32, e42, e43
      e21 = es(2) - es(1); e31 = es(3) - es(1); e41 = es(4) - es(1)
      e32 = es(3) - es(2); e42 = es(4) - es(2); e43 = es(4) - es(3)
      if (e < es(1)) then
         n = 0.0_rp
      else if (e < es(2)) then
         n = (e - es(1))**3/(e21*e31*e41)
      else if (e < es(3)) then
         n = (e21**2 + 3.0_rp*e21*(e - es(2)) + 3.0_rp*(e - es(2))**2 - (e31 + e42)/(e32*e42)*(e - es(2))**3)/(e31*e41)
      else if (e < es(4)) then
         n = 1.0_rp - (es(4) - e)**3/(e41*e42*e43)
      else
         n = 1.0_rp
      end if
   end function nos_closed

   !> Uncorrected integration weights (B1)-(B19) per unit volume, returned in the caller's corner order.
   function bweights(ecall, e) result(wc)
      real(rp), intent(in) :: ecall(4), e
      real(rp) :: wc(4), w(4), es(4), c, c1, c2, c3, e21, e31, e41, e32, e42, e43
      integer :: ord(4), k
      call sort_idx(ecall, es, ord)
      e21 = es(2) - es(1); e31 = es(3) - es(1); e41 = es(4) - es(1)
      e32 = es(3) - es(2); e42 = es(4) - es(2); e43 = es(4) - es(3)
      if (e < es(1)) then
         w = 0.0_rp
      else if (e < es(2)) then
         c = 0.25_rp*(e - es(1))**3/(e21*e31*e41)
         w(1) = c*(4.0_rp - (e - es(1))*(1.0_rp/e21 + 1.0_rp/e31 + 1.0_rp/e41))
         w(2) = c*(e - es(1))/e21
         w(3) = c*(e - es(1))/e31
         w(4) = c*(e - es(1))/e41
      else if (e < es(3)) then
         c1 = 0.25_rp*(e - es(1))**2/(e41*e31)
         c2 = 0.25_rp*(e - es(1))*(e - es(2))*(es(3) - e)/(e41*e32*e31)
         c3 = 0.25_rp*(e - es(2))**2*(es(4) - e)/(e42*e32*e41)
         w(1) = c1 + (c1 + c2)*(es(3) - e)/e31 + (c1 + c2 + c3)*(es(4) - e)/e41
         w(2) = c1 + c2 + c3 + (c2 + c3)*(es(3) - e)/e32 + c3*(es(4) - e)/e42
         w(3) = (c1 + c2)*(e - es(1))/e31 + (c2 + c3)*(e - es(2))/e32
         w(4) = (c1 + c2 + c3)*(e - es(1))/e41 + c3*(e - es(2))/e42
      else if (e < es(4)) then
         c = 0.25_rp*(es(4) - e)**3/(e41*e42*e43)
         w(1) = 0.25_rp - c*(es(4) - e)/e41
         w(2) = 0.25_rp - c*(es(4) - e)/e42
         w(3) = 0.25_rp - c*(es(4) - e)/e43
         w(4) = 0.25_rp - c*(4.0_rp - (1.0_rp/e41 + 1.0_rp/e42 + 1.0_rp/e43)*(es(4) - e))
      else
         w = 0.25_rp
      end if
      do k = 1, 4
         wc(ord(k)) = w(k)
      end do
   end function bweights

   !> Correction Eq. (22) per unit volume at energy e, caller corner order.
   function corr22(ecall, e) result(dw)
      real(rp), intent(in) :: ecall(4), e
      real(rp) :: dw(4), es(4)
      integer :: ord(4)
      call sort_idx(ecall, es, ord)
      dw = dos_closed(es, e)/40.0_rp*(sum(ecall) - 4.0_rp*ecall)
   end function corr22

   subroutine sort_idx(e, es, ord)
      real(rp), intent(in) :: e(4)
      real(rp), intent(out) :: es(4)
      integer, intent(out) :: ord(4)
      integer :: i, j, t
      ord = [1, 2, 3, 4]
      do i = 2, 4
         do j = i, 2, -1
            if (e(ord(j - 1)) > e(ord(j))) then
               t = ord(j); ord(j) = ord(j - 1); ord(j - 1) = t
            end if
         end do
      end do
      es = e(ord)
   end subroutine sort_idx

   subroutine sorted(e, es)
      real(rp), intent(in) :: e(4)
      real(rp), intent(out) :: es(4)
      integer :: ord(4)
      call sort_idx(e, es, ord)
   end subroutine sorted

   subroutine check_sum_dos(e, worst)
      real(rp), intent(in) :: e(4)
      real(rp), intent(inout) :: worst
      real(rp) :: es(4), w(4), dw(4), en
      integer :: k
      call sorted(e, es)
      do k = 1, 40
         en = es(1) + (es(4) - es(1))*(real(k, rp) - 0.37_rp)/40.0_rp
         call tetra_corner_weights(e, en, w, dw)
         worst = max(worst, abs(sum(w) - dos_closed(es, en)), abs(sum(dw)), abs(sum(bweights(e, en)) - nos_closed(es, en)))
      end do
   end subroutine check_sum_dos

   !> int_{e1}^{e4} w_i(E) E^m dE against <lambda_i eps^m> over the unit-volume tetrahedron.
   subroutine check_moments(e, mpow, worst)
      real(rp), intent(in) :: e(4)
      integer, intent(in) :: mpow
      real(rp), intent(inout) :: worst
      real(rp) :: es(4), w(4), dw(4), lhs(4), rhs(4), a, b, en, lam(4), eps_k
      integer :: i, k, seg, p, q, r
      real(rp) :: u, v, ww, jac, x0, x1, x2, x3
      call sorted(e, es)
      lhs = 0.0_rp
      do seg = 1, 3
         a = es(seg); b = es(seg + 1)
         if (b - a <= 0.0_rp) cycle
         do k = 1, ngl
            en = 0.5_rp*(a + b) + 0.5_rp*(b - a)*gx(k)
            call tetra_corner_weights(e, en, w, dw)
            lhs = lhs + 0.5_rp*(b - a)*gw(k)*w*en**mpow
         end do
      end do
      rhs = 0.0_rp
      do p = 1, ngl
         do q = 1, ngl
            do r = 1, ngl
               u = 0.5_rp*(gx(p) + 1.0_rp); v = 0.5_rp*(gx(q) + 1.0_rp); ww = 0.5_rp*(gx(r) + 1.0_rp)
               x0 = u; x1 = (1.0_rp - u)*v; x2 = (1.0_rp - u)*(1.0_rp - v)*ww; x3 = 1.0_rp - x0 - x1 - x2
               lam = [x0, x1, x2, x3]
               jac = (1.0_rp - u)**2*(1.0_rp - v)*0.125_rp*gw(p)*gw(q)*gw(r)*6.0_rp
               eps_k = sum(lam*e)
               rhs = rhs + jac*lam*eps_k**mpow
            end do
         end do
      end do
      if (mpow == 0 .and. maxval(abs(rhs - 0.25_rp)) > 1.0e-13_rp) error stop 'simplex quadrature self-check failed'
      do i = 1, 4
         worst = max(worst, abs(lhs(i) - rhs(i)))
      end do
   end subroutine check_moments

   !> int_{-inf}^{upto} of the kernel's w_lin and dw_corr against (B1)-(B19) and Eq. (22).
   subroutine check_corr_integral(e, upto, worst)
      real(rp), intent(in) :: e(4), upto
      real(rp), intent(inout) :: worst
      real(rp) :: es(4), w(4), dw(4), accw(4), accd(4), a, b, en
      integer :: k, seg
      call sorted(e, es)
      accw = 0.0_rp
      accd = 0.0_rp
      do seg = 1, 3
         a = es(seg); b = min(es(seg + 1), upto)
         if (b <= a) cycle
         do k = 1, ngl
            en = 0.5_rp*(a + b) + 0.5_rp*(b - a)*gx(k)
            call tetra_corner_weights(e, en, w, dw)
            accw = accw + 0.5_rp*(b - a)*gw(k)*w
            accd = accd + 0.5_rp*(b - a)*gw(k)*dw
         end do
      end do
      worst = max(worst, maxval(abs(accw - bweights(e, upto))), maxval(abs(accd - corr22(e, upto))))
   end subroutine check_corr_integral

   !> Corner-weighted integrals sums(1:4, 1:2) = [count, f1, f2, f3] (linear, corrected) at Fermi energy efe on an n^3 mesh.
   subroutine mesh_sums(n, efe, sums)
      integer, intent(in) :: n
      real(rp), intent(in) :: efe
      real(rp), intent(out) :: sums(4, 2)
      real(rp), parameter :: k0(3) = [0.13_rp, -0.07_rp, 0.21_rp]
      integer, parameter :: perm(3, 6) = reshape([1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1], [3, 6])
      real(rp) :: h, vt, kc(3, 4), ec(4), fc(4, 4), w(4), dw(4), a, b, en, wint(4), cint(4), es(4)
      real(rp) :: g6x(6), g6w(6)
      integer :: ix, iy, iz, ip, c, ax, seg, k, q
      integer :: off(3)

      call gauss_legendre(6, g6x, g6w)
      h = 2.0_rp/real(n, rp)
      vt = h**3/6.0_rp
      sums = 0.0_rp
      do iz = 0, n - 1
         do iy = 0, n - 1
            do ix = 0, n - 1
               do ip = 1, 6
                  off = 0
                  kc(:, 1) = -1.0_rp + h*real([ix, iy, iz], rp)
                  do c = 2, 4
                     if (c < 4) then
                        ax = perm(c - 1, ip)
                        off(ax) = 1
                     else
                        off = 1
                     end if
                     kc(:, c) = -1.0_rp + h*real([ix, iy, iz] + off, rp)
                  end do
                  do c = 1, 4
                     ec(c) = sum((kc(:, c) - k0)**2)
                     fc(c, 1) = 1.0_rp
                     fc(c, 2) = 1.0_rp + 0.5_rp*kc(1, c) - 0.8_rp*kc(2, c) + 0.3_rp*kc(3, c)
                     fc(c, 3) = kc(1, c)**2 + 0.5_rp*kc(2, c)*kc(3, c)
                     fc(c, 4) = ec(c)
                  end do
                  call sorted(ec, es)
                  wint = 0.0_rp; cint = 0.0_rp
                  do seg = 1, 3
                     a = es(seg); b = min(es(seg + 1), efe)
                     if (b <= a) cycle
                     do k = 1, 6
                        en = 0.5_rp*(a + b) + 0.5_rp*(b - a)*g6x(k)
                        call tetra_corner_weights(ec, en, w, dw)
                        wint = wint + 0.5_rp*(b - a)*g6w(k)*w
                        cint = cint + 0.5_rp*(b - a)*g6w(k)*dw
                     end do
                  end do
                  do q = 1, 4
                     sums(q, 1) = sums(q, 1) + vt*sum(wint*fc(:, q))
                     sums(q, 2) = sums(q, 2) + vt*sum((wint + cint)*fc(:, q))
                  end do
               end do
            end do
         end do
      end do
   end subroutine mesh_sums

   !> At the Fermi energy where the linear count equals the exact count the band energy (f = eps) and f2 are compared with
   !> the ball integrals: R^2 = E_F, int f = ..., f1 linear exact f(k0)*volume.
   subroutine check_mesh()
      real(rp), parameter :: k0(3) = [0.13_rp, -0.07_rp, 0.21_rp], ef = 0.45_rp
      integer, parameter :: ns(4) = [4, 8, 16, 32]
      real(rp) :: ball, ex(4), sums(4, 2), lo, hi, mid, prev
      integer :: in, it

      ball = 4.0_rp*pi/3.0_rp*ef**1.5_rp
      ex(1) = ball
      ex(2) = (1.0_rp + 0.5_rp*k0(1) - 0.8_rp*k0(2) + 0.3_rp*k0(3))*ball
      ex(3) = (k0(1)**2 + 0.5_rp*k0(2)*k0(3))*ball + 4.0_rp*pi*ef**2.5_rp/15.0_rp
      ex(4) = 4.0_rp*pi*ef**2.5_rp/5.0_rp
      prev = huge(1.0_rp)
      print '(a)', '   n  efermi_lin   err_lin(count,f1,f2,eps)  |  err_corr(count,f1,f2,eps)   (linear count = exact count)'
      do in = 1, size(ns)
         lo = 0.0_rp; hi = 1.0_rp
         do it = 1, 40
            mid = 0.5_rp*(lo + hi)
            call mesh_sums(ns(in), mid, sums)
            if (sums(1, 1) < ex(1)) then
               lo = mid
            else
               hi = mid
            end if
            if (hi - lo < 1.0e-9_rp) exit
         end do
         call mesh_sums(ns(in), 0.5_rp*(lo + hi), sums)
         print '(i4,f10.6,4es13.3,a,4es13.3)', ns(in), 0.5_rp*(lo + hi), sums(:, 1) - ex, ' |', sums(:, 2) - ex
         if (abs(sums(1, 2) - ex(1)) > 1.0e-7_rp) error stop 'Fermi-level bisection did not reach the exact count'
         if (abs(sums(4, 2) - ex(4)) >= abs(sums(4, 1) - ex(4))) error stop 'corrected band energy not closer than linear'
         if (abs(sums(3, 2) - ex(3)) >= abs(sums(3, 1) - ex(3))) error stop 'corrected f2 integral not closer than linear'
         if (ns(in) >= 8 .and. abs(sums(4, 2) - ex(4)) > 0.1_rp*abs(sums(4, 1) - ex(4))) &
            error stop 'corrected band-energy error not 10x below linear'
         if (abs(sums(4, 2) - ex(4)) >= prev) error stop 'corrected band-energy error does not decrease with mesh size'
         prev = abs(sums(4, 2) - ex(4))
      end do
   end subroutine check_mesh

   subroutine gauss_legendre(n, x, w)
      integer, intent(in) :: n
      real(rp), intent(out) :: x(n), w(n)
      integer :: i, it
      real(rp) :: z, z1, p1, p2, p3, pp
      integer :: j
      do i = 1, n
         z = cos(pi*(real(i, rp) - 0.25_rp)/(real(n, rp) + 0.5_rp))
         do it = 1, 100
            p1 = 1.0_rp; p2 = 0.0_rp
            do j = 1, n
               p3 = p2; p2 = p1
               p1 = ((2.0_rp*real(j, rp) - 1.0_rp)*z*p2 - (real(j, rp) - 1.0_rp)*p3)/real(j, rp)
            end do
            pp = real(n, rp)*(z*p1 - p2)/(z*z - 1.0_rp)
            z1 = z
            z = z1 - p1/pp
            if (abs(z - z1) < 1.0e-15_rp) exit
         end do
         x(i) = z
         w(i) = 2.0_rp/((1.0_rp - z*z)*pp*pp)
      end do
   end subroutine gauss_legendre

end program test_tetra_corner_weights
