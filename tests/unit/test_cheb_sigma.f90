!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
!> @brief Chebyshev expansion with a complex centre against closed-form resolvents.
!> @details Moments T_n(Kt) are computed here by the three-term recurrence; the
!>          references are the closed-form resolvents written in this test.
!>          T1: scalar, Sigma = 0, z = 0.05i; a = (emax-emin)/1.7 and Im b = 0.
!>              Catches: wrong branch of w (|w|>1), wrong 1/a factor.
!>          T2: Sigma = diag(-0.6i, 0), scalar x = 0.2-0.6i, z = 0.5i, rho_z > rho_s (z = 0.1i lies inside the enclosing ellipse).
!>              Catches: missing imaginary part of the centre b (series diverges).
!>          T3: 2x2 non-normal K = H + diag(s1, 0), z = 0.1+0.2i (z = 0.1+0.05i lies inside the enclosing ellipse), against the 2x2 inverse.
!>              Catches: transposed or conjugated blocks, wrong window from the Sigma blocks.
!>          T4: Sigma = diag(-50i, 0), x = 0.1-50i, renormalised moments plus logscale.
!>              Catches: missing log scaling (overflow of T_n), exponent combined after exp.
!>          T5: truncated T1 and T3: measured error <= 10 tail and tail < 1e3 error.
!>              Catches: tail bound that is wrong in sign, size or normalisation.
!>          T6: non-causal Sigma (Im +0.5), z = 0.05i inside the region: rho_z <= rho_s and
!>              tail = huge. Catches: divergence not flagged.
!------------------------------------------------------------------------------
program test_cheb_sigma
   use precision_mod, only: rp
   use cheb_sigma_kernel_mod, only: cheb_sigma_window, cheb_sigma_w, cheb_sigma_sum
   implicit none

   logical :: failed
   failed = .false.

   call test_t1()
   call test_t2()
   call test_t3()
   call test_t4()
   call test_t6()
   if (failed) error stop 1

contains

   subroutine report(name, value, tol, ok)
      character(len=*), intent(in) :: name
      real(rp), intent(in) :: value, tol
      logical, intent(in) :: ok
      if (ok) then
         print '(a,a,es12.4,a,es10.2)', 'PASS ', name, value, '  tol ', tol
      else
         print '(a,a,es12.4,a,es10.2)', 'FAIL ', name, value, '  tol ', tol
         failed = .true.
      end if
   end subroutine report

   !> Number of moments with r^n < eps, at least 8.
   function n_for(r, eps) result(n)
      real(rp), intent(in) :: r, eps
      integer :: n
      n = max(8, ceiling(log(eps)/log(r)))
   end function n_for

   !> Scalar moments T_n(xt), n = 0..nmom-1, by plain recurrence.
   subroutine scalar_moments(xt, nmom, mu)
      complex(rp), intent(in) :: xt
      integer, intent(in) :: nmom
      complex(rp), intent(out) :: mu(1, 1, nmom)
      integer :: k
      mu(1, 1, 1) = 1.0_rp
      mu(1, 1, 2) = xt
      do k = 3, nmom
         mu(1, 1, k) = 2.0_rp*xt*mu(1, 1, k - 1) - mu(1, 1, k - 2)
      end do
   end subroutine scalar_moments

   subroutine test_t1()
      real(rp), parameter :: emin = -1.0_rp, emax = 1.0_rp
      complex(rp) :: sig(1, 1, 1), b, w, z, g(1, 1), ref
      complex(rp), allocatable :: mu(:, :, :)
      real(rp) :: a, rho_s, rho_z, tail, err
      real(rp), allocatable :: ls(:)
      integer :: nmom, pass
      sig = (0.0_rp, 0.0_rp)
      z = cmplx(0.0_rp, 0.05_rp, rp)
      call cheb_sigma_window(emin, emax, sig, a, b, rho_s)
      call report('T1 window a-(emax-emin)/1.7', abs(a - (emax - emin)/1.7_rp), 1.0e-15_rp, &
                  abs(a - (emax - emin)/1.7_rp) <= 1.0e-15_rp)
      call report('T1 window aimag(b)', abs(aimag(b)), 0.0_rp, aimag(b) == 0.0_rp)
      call cheb_sigma_w(z, a, b, w, rho_z)
      ref = 1.0_rp/(z - 0.3_rp)
      do pass = 1, 2
         nmom = n_for(rho_s/rho_z, merge(1.0e-14_rp, 1.0e-6_rp, pass == 1))
         allocate (mu(1, 1, nmom), ls(nmom))
         ls = 0.0_rp
         call scalar_moments((cmplx(0.3_rp, 0.0_rp, rp) - b)/a, nmom, mu)
         call cheb_sigma_sum(mu, ls, w, a, rho_s/rho_z, g, tail)
         if (pass == 1) then
            err = abs(g(1, 1) - ref)/abs(ref)
            call report('T1 scalar Sigma=0', err, 1.0e-10_rp, err < 1.0e-10_rp)
         else
            err = abs(g(1, 1) - ref)/abs(g(1, 1))
            call report('T5 T1 err', err, 10.0_rp*tail, err <= 10.0_rp*tail)
            call report('T5 T1 tail', tail, 1.0e3_rp*err, tail < 1.0e3_rp*err)
         end if
         deallocate (mu, ls)
      end do
   end subroutine test_t1

   subroutine test_t2()
      complex(rp) :: sig(2, 2, 1), b, w, z, g(1, 1), ref, x
      complex(rp), allocatable :: mu(:, :, :)
      real(rp) :: a, rho_s, rho_z, tail, err
      real(rp), allocatable :: ls(:)
      integer :: nmom
      sig = (0.0_rp, 0.0_rp)
      sig(1, 1, 1) = cmplx(0.0_rp, -0.6_rp, rp)
      x = cmplx(0.2_rp, -0.6_rp, rp)
      z = cmplx(0.0_rp, 0.5_rp, rp)
      call cheb_sigma_window(-1.0_rp, 1.0_rp, sig, a, b, rho_s)
      call cheb_sigma_w(z, a, b, w, rho_z)
      call report('T2 rho_z - rho_s', rho_z - rho_s, 0.0_rp, rho_z > rho_s)
      nmom = n_for(rho_s/rho_z, 1.0e-14_rp)
      allocate (mu(1, 1, nmom), ls(nmom))
      ls = 0.0_rp
      call scalar_moments((x - b)/a, nmom, mu)
      call cheb_sigma_sum(mu, ls, w, a, rho_s/rho_z, g, tail)
      ref = 1.0_rp/(z - x)
      err = abs(g(1, 1) - ref)/abs(ref)
      call report('T2 damped scalar', err, 1.0e-10_rp, err < 1.0e-10_rp)
   end subroutine test_t2

   subroutine test_t3()
      real(rp), parameter :: e1 = -0.2_rp, e2 = 0.4_rp, t = 0.3_rp
      complex(rp) :: sig(2, 2, 1), b, w, z, g(2, 2), kt(2, 2), k(2, 2), ref(2, 2), det
      complex(rp), allocatable :: mu(:, :, :)
      real(rp) :: a, rho_s, rho_z, tail, err, emid, erad
      real(rp), allocatable :: ls(:)
      integer :: nmom, n, pass
      emid = 0.5_rp*(e1 + e2)
      erad = sqrt((0.5_rp*(e1 - e2))**2 + t*t)
      sig = (0.0_rp, 0.0_rp)
      sig(1, 1, 1) = cmplx(0.1_rp, -0.4_rp, rp)
      z = cmplx(0.1_rp, 0.2_rp, rp)
      call cheb_sigma_window(emid - erad, emid + erad, sig, a, b, rho_s)
      call cheb_sigma_w(z, a, b, w, rho_z)
      call report('T3 rho_z - rho_s', rho_z - rho_s, 0.0_rp, rho_z > rho_s)
      k = reshape([cmplx(e1, 0.0_rp, rp) + sig(1, 1, 1), cmplx(t, 0.0_rp, rp), cmplx(t, 0.0_rp, rp), &
                   cmplx(e2, 0.0_rp, rp)], [2, 2])
      kt = k
      kt(1, 1) = kt(1, 1) - b
      kt(2, 2) = kt(2, 2) - b
      kt = kt/a
      det = (z - k(1, 1))*(z - k(2, 2)) - k(1, 2)*k(2, 1)
      ref = reshape([z - k(2, 2), k(2, 1), k(1, 2), z - k(1, 1)], [2, 2])/det
      do pass = 1, 2
         nmom = merge(n_for(rho_s/rho_z, 1.0e-14_rp), 28, pass == 1)
         allocate (mu(2, 2, nmom), ls(nmom))
         ls = 0.0_rp
         mu(:, :, 1) = reshape([1.0_rp, 0.0_rp, 0.0_rp, 1.0_rp], [2, 2])
         mu(:, :, 2) = kt
         do n = 3, nmom
            mu(:, :, n) = 2.0_rp*matmul(kt, mu(:, :, n - 1)) - mu(:, :, n - 2)
         end do
         call cheb_sigma_sum(mu, ls, w, a, rho_s/rho_z, g, tail)
         if (pass == 1) then
            err = maxval(abs(g - ref))/maxval(abs(ref))
            call report('T3 non-normal 2x2', err, 1.0e-10_rp, err < 1.0e-10_rp)
         else
            err = maxval(abs(g - ref))/maxval(abs(g))
            call report('T5 T3 err', err, 10.0_rp*tail, err <= 10.0_rp*tail)
            call report('T5 T3 tail', tail, 1.0e3_rp*err, tail < 1.0e3_rp*err)
         end if
         deallocate (mu, ls)
      end do
   end subroutine test_t3

   subroutine test_t4()
      complex(rp) :: sig(2, 2, 1), b, w, z, g(1, 1), ref, x, xt, p, q, pn
      complex(rp), allocatable :: mu(:, :, :)
      real(rp), allocatable :: ls(:)
      real(rp) :: a, rho_s, rho_z, tail, err, lscale, m
      integer :: nmom, k
      sig = (0.0_rp, 0.0_rp)
      sig(1, 1, 1) = cmplx(0.0_rp, -50.0_rp, rp)
      x = cmplx(0.1_rp, -50.0_rp, rp)
      z = cmplx(0.0_rp, 0.1_rp, rp)
      call cheb_sigma_window(-1.0_rp, 1.0_rp, sig, a, b, rho_s)
      call cheb_sigma_w(z, a, b, w, rho_z)
      nmom = n_for(rho_s/rho_z, 1.0e-10_rp)
      allocate (mu(1, 1, nmom), ls(nmom))
      xt = (x - b)/a
      p = 1.0_rp
      q = xt
      lscale = 0.0_rp
      do k = 1, nmom
         mu(1, 1, k) = p
         ls(k) = lscale
         pn = 2.0_rp*xt*q - p
         p = q
         q = pn
         m = max(abs(p), abs(q))
         p = p/m
         q = q/m
         lscale = lscale + log(m)
      end do
      call cheb_sigma_sum(mu, ls, w, a, rho_s/rho_z, g, tail)
      ref = 1.0_rp/(z - x)
      err = abs(g(1, 1) - ref)/abs(ref)
      print '(a,i0,a,es10.3)', 'T4 nmom ', nmom, '  final logscale ', lscale
      call report('T4 strong damping', err, 1.0e-8_rp, err < 1.0e-8_rp)
   end subroutine test_t4

   subroutine test_t6()
      complex(rp) :: sig(2, 2, 1), b, w, z, g(2, 2), mu(2, 2, 5)
      real(rp) :: a, rho_s, rho_z, tail, ls(5)
      sig = (0.0_rp, 0.0_rp)
      sig(1, 1, 1) = cmplx(0.0_rp, 0.5_rp, rp)
      z = cmplx(0.0_rp, 0.05_rp, rp)
      call cheb_sigma_window(-1.0_rp, 1.0_rp, sig, a, b, rho_s)
      call cheb_sigma_w(z, a, b, w, rho_z)
      mu = (1.0_rp, 0.0_rp)
      ls = 0.0_rp
      call cheb_sigma_sum(mu, ls, w, a, rho_s/rho_z, g, tail)
      call report('T6 rho_z - rho_s', rho_z - rho_s, 0.0_rp, rho_z <= rho_s)
      call report('T6 tail flag', tail, huge(1.0_rp), tail == huge(1.0_rp))
   end subroutine test_t6

end program test_cheb_sigma
