!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
! PROGRAM: test_bsf_sumrule
!
!> @brief Bloch spectral function known-answer tests (B3.1), dependency-free.
!> @details Exercises the BSF core `bsf_spectral_trace` on the retarded resolvent
!>          from `dyson_kspace_inverse` (backend D, Sigma=0), with NO LMTO
!>          machinery -- the same 1-band tight-binding fixtures the B2 tests use.
!>
!>          Pins:
!>          1. Sum rule: for the 1-band chain eps(k) = -2t cos k, at each k
!>             ∫ A(k,E) dE = N_orb = 1 (the diagonal-resolvent spectral sum rule,
!>             a Lorentzian of unit weight), trapezoid on a wide real-axis grid.
!>          2. Partial-trace additivity + per-channel sum rule on a coupled 2x2
!>             H(k): A_total(E) = A_{orb1}(E) + A_{orb2}(E) pointwise, and each
!>             orbital-projected ∫ A dE = 1 (so ∫ A_total dE = N_orb = 2).
!>          3. Delta limit: for small eta the single-k A(E) peaks at the
!>             eigenvalue, |argmax_E A - eps(k)| < 2*eta.
!>
!>          Exits non-zero (error stop) on any failure so ctest registers a fail.
!------------------------------------------------------------------------------
program test_bsf_sumrule
   use precision_mod, only: rp
   use math_mod, only: pi, two_pi
   use bsf_kernel_mod, only: bsf_spectral_trace
   use dyson_kernel_mod, only: dyson_kspace_inverse
   implicit none

   real(rp), parameter :: t = 1.0_rp
   real(rp), parameter :: tol_sum = 5.0e-3_rp      ! eta-tail truncation dominated
   real(rp), parameter :: tol_add = 1.0e-12_rp     ! pointwise additivity, exact
   logical :: failed

   failed = .false.

   call test_chain_sum_rule()
   call test_partial_trace()
   call test_delta_limit()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   !> Test 1: ∫ A(k,E) dE = 1 for the 1-band chain at several k-points.
   subroutine test_chain_sum_rule()
      integer, parameter :: nk = 5
      real(rp), parameter :: eta = 0.10_rp
      real(rp), parameter :: emin = -40.0_rp, emax = 40.0_rp
      integer, parameter :: ne = 16001
      complex(rp) :: hk(1, 1), sig0(1, 1), gk(1, 1)
      real(rp) :: kf, eps, de, e, a, a_prev, weight, max_err
      integer :: ik, ie

      sig0 = (0.0_rp, 0.0_rp)
      de = (emax - emin)/real(ne - 1, rp)
      max_err = 0.0_rp
      do ik = 1, nk
         kf = real(ik - 1, rp)/real(nk, rp)
         eps = -2.0_rp*t*cos(two_pi*kf)
         hk(1, 1) = cmplx(eps, 0.0_rp, rp)
         weight = 0.0_rp
         a_prev = 0.0_rp
         do ie = 1, ne
            e = emin + de*real(ie - 1, rp)
            call dyson_kspace_inverse(hk, cmplx(e, eta, rp), sig0, gk)
            a = bsf_spectral_trace(gk, [1])
            if (ie > 1) weight = weight + 0.5_rp*de*(a + a_prev)
            a_prev = a
         end do
         max_err = max(max_err, abs(weight - 1.0_rp))
      end do
      write (*, '(a,es12.4)') 'Test 1 (chain sum rule ∫A dE = 1)   max_err = ', max_err
      if (max_err > tol_sum) failed = .true.
   end subroutine test_chain_sum_rule

   !> Test 2: partial-trace additivity + per-orbital sum rule on a coupled 2x2 H.
   subroutine test_partial_trace()
      real(rp), parameter :: eta = 0.10_rp
      real(rp), parameter :: emin = -40.0_rp, emax = 40.0_rp
      integer, parameter :: ne = 16001
      complex(rp) :: hk(2, 2), sig0(2, 2), gk(2, 2)
      real(rp) :: de, e, a_tot, a1, a2, w1, w2, a1_prev, a2_prev
      real(rp) :: max_add_err, max_sum_err
      integer :: ie

      hk = (0.0_rp, 0.0_rp)
      hk(1, 1) = cmplx(-0.7_rp, 0.0_rp, rp)
      hk(2, 2) = cmplx(1.3_rp, 0.0_rp, rp)
      hk(1, 2) = cmplx(0.4_rp, 0.0_rp, rp)
      hk(2, 1) = cmplx(0.4_rp, 0.0_rp, rp)
      sig0 = (0.0_rp, 0.0_rp)

      de = (emax - emin)/real(ne - 1, rp)
      max_add_err = 0.0_rp
      w1 = 0.0_rp; w2 = 0.0_rp
      a1_prev = 0.0_rp; a2_prev = 0.0_rp
      do ie = 1, ne
         e = emin + de*real(ie - 1, rp)
         call dyson_kspace_inverse(hk, cmplx(e, eta, rp), sig0, gk)
         a_tot = bsf_spectral_trace(gk, [1, 2])
         a1 = bsf_spectral_trace(gk, [1])
         a2 = bsf_spectral_trace(gk, [2])
         max_add_err = max(max_add_err, abs(a_tot - (a1 + a2)))
         if (ie > 1) then
            w1 = w1 + 0.5_rp*de*(a1 + a1_prev)
            w2 = w2 + 0.5_rp*de*(a2 + a2_prev)
         end if
         a1_prev = a1; a2_prev = a2
      end do
      max_sum_err = max(abs(w1 - 1.0_rp), abs(w2 - 1.0_rp))
      write (*, '(a,es12.4)') 'Test 2 (partial-trace additivity)   max_err = ', max_add_err
      write (*, '(a,es12.4)') 'Test 2 (per-orbital sum rule = 1)   max_err = ', max_sum_err
      if (max_add_err > tol_add) failed = .true.
      if (max_sum_err > tol_sum) failed = .true.
   end subroutine test_partial_trace

   !> Test 3: delta limit -- A(E) peaks at eps(k) to within 2*eta.
   subroutine test_delta_limit()
      real(rp), parameter :: eta = 0.02_rp
      real(rp), parameter :: emin = -4.0_rp, emax = 4.0_rp
      integer, parameter :: ne = 8001
      complex(rp) :: hk(1, 1), sig0(1, 1), gk(1, 1)
      real(rp) :: kf, eps, de, e, a, a_max, e_at_max
      integer :: ie

      sig0 = (0.0_rp, 0.0_rp)
      kf = 0.30_rp
      eps = -2.0_rp*t*cos(two_pi*kf)
      hk(1, 1) = cmplx(eps, 0.0_rp, rp)
      de = (emax - emin)/real(ne - 1, rp)
      a_max = -1.0_rp; e_at_max = emin
      do ie = 1, ne
         e = emin + de*real(ie - 1, rp)
         call dyson_kspace_inverse(hk, cmplx(e, eta, rp), sig0, gk)
         a = bsf_spectral_trace(gk, [1])
         if (a > a_max) then
            a_max = a; e_at_max = e
         end if
      end do
      write (*, '(a,es12.4,a,es12.4)') 'Test 3 (delta limit) |argmax A - eps| = ', &
         abs(e_at_max - eps), '   2*eta = ', 2.0_rp*eta
      if (abs(e_at_max - eps) > 2.0_rp*eta) failed = .true.
   end subroutine test_delta_limit

end program test_bsf_sumrule
