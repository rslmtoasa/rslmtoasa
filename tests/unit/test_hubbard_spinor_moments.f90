! Spinor local-moment kernel (correct_d_shell_moments_spinor) and the spinor R05 residual.
!
! (a) block-diagonal inputs reproduce the collinear kernel; (c) q2' >= 0 for states built from the exact
! eigenpairs of H0 + V (random Hermitian, with spin-off-diagonal V); R05: a sequence whose only change is
! in the spin-off-diagonal blocks of P and V is not declared converged. The finite-H comparison with the
! direct definition of c' and q2' is in projector_spinor_frozen_h.f90. Every control must exceed its tolerance.
program test_hubbard_spinor_moments
   use precision_mod, only: rp
   use math_mod, only: correct_d_shell_moments, correct_d_shell_moments_spinor
   use mix_mod, only: mix, fixed_u_residual_metrics
   implicit none

   integer(8) :: seed = 20261027_8
   logical :: failed = .false.
   external :: zheev

   call test_collinear_reduction()
   call test_q2_nonnegative()
   call test_r05_off_diagonal_blocks()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   function uniform() result(x)
      real(rp) :: x

      seed = mod(seed*1103515245_8 + 12345_8, 2147483648_8)
      x = real(seed, rp)/2147483648.0_rp
   end function uniform

   function rand_hermitian(n, scale) result(h)
      integer, intent(in) :: n
      real(rp), intent(in) :: scale
      complex(rp) :: h(n, n)
      integer :: i, j

      do i = 1, n
         h(i, i) = cmplx(scale*(uniform() - 0.5_rp), 0.0_rp, rp)
         do j = i + 1, n
            h(i, j) = scale*cmplx(uniform() - 0.5_rp, uniform() - 0.5_rp, rp)
            h(j, i) = conjg(h(i, j))
         end do
      end do
   end function rand_hermitian

   subroutine eigenpairs(h, e, v)
      complex(rp), intent(in) :: h(:, :)
      real(rp), intent(out) :: e(size(h, 1))
      complex(rp), intent(out) :: v(size(h, 1), size(h, 1))
      complex(rp) :: work(4*size(h, 1))
      real(rp) :: rwork(3*size(h, 1))
      integer :: info

      v = h
      call zheev('V', 'U', size(h, 1), v, size(h, 1), e, work, 4*size(h, 1), rwork, info)
      if (info /= 0) failed = .true.
   end subroutine eigenpairs

   !> P, M1 and the per-spin q0, c, q2 of the occupied eigenstates (the lowest nocc) of h.
   subroutine moments_of(h, nocc, p, m1, q0, c, q2)
      complex(rp), intent(in) :: h(10, 10)
      integer, intent(in) :: nocc
      complex(rp), intent(out) :: p(10, 10), m1(10, 10)
      real(rp), intent(out) :: q0(2), c(2), q2(2)
      complex(rp) :: v(10, 10)
      real(rp) :: e(10), w
      integer :: n, s, k

      call eigenpairs(h, e, v)
      p = (0.0_rp, 0.0_rp)
      m1 = (0.0_rp, 0.0_rp)
      do n = 1, nocc
         p = p + spread(v(:, n), 2, 10)*spread(conjg(v(:, n)), 1, 10)
         m1 = m1 + e(n)*spread(v(:, n), 2, 10)*spread(conjg(v(:, n)), 1, 10)
      end do
      do s = 1, 2
         q0(s) = 0.0_rp
         c(s) = 0.0_rp
         do k = (s - 1)*5 + 1, s*5
            q0(s) = q0(s) + real(p(k, k), rp)
            c(s) = c(s) + real(m1(k, k), rp)
         end do
         c(s) = c(s)/q0(s)
         q2(s) = 0.0_rp
         do n = 1, nocc
            w = sum(abs(v((s - 1)*5 + 1:s*5, n))**2)
            q2(s) = q2(s) + (e(n) - c(s))**2*w
         end do
      end do
   end subroutine moments_of

   subroutine check(name, measured, tol, control_check)
      character(len=*), intent(in) :: name
      real(rp), intent(in) :: measured, tol
      logical, intent(in), optional :: control_check
      logical :: ok, is_control

      is_control = .false.
      if (present(control_check)) is_control = control_check
      ok = merge(measured > tol, measured <= tol, is_control)
      if (.not. ok) failed = .true.
      if (is_control) then
         write (*, '(a,1x,a,2(2x,es10.3))') merge('PASS', 'FAIL', ok), trim(name)//' [control: must exceed]', measured, tol
      else
         write (*, '(a,1x,a,2(2x,es10.3))') merge('PASS', 'FAIL', ok), trim(name), measured, tol
      end if
   end subroutine check

   subroutine test_collinear_reduction()
      complex(rp) :: h(10, 10), p(10, 10), m1(10, 10), v(10, 10), h0(10, 10), vsw(10, 10)
      real(rp) :: q0(2), c(2), q2(2), cp(2), q2p(2), t(2), cp2(2), q2p2(2), t2(2)
      real(rp) :: cs(2), q2s(2), ts(2), worst, worst_swap, worst_off
      integer :: s, status, lo, hi

      ! block-diagonal H0 and V: the spinor states are those of the two spin blocks
      h0 = (0.0_rp, 0.0_rp)
      v = (0.0_rp, 0.0_rp)
      h0(1:5, 1:5) = rand_hermitian(5, 1.0_rp)
      h0(6:10, 6:10) = rand_hermitian(5, 1.0_rp)
      v(1:5, 1:5) = rand_hermitian(5, 0.4_rp)
      v(6:10, 6:10) = rand_hermitian(5, 0.4_rp)
      h = h0 + v
      call moments_of(h, 6, p, m1, q0, c, q2)
      call correct_d_shell_moments_spinor(p, m1, v, q0, c, q2, cp, q2p, t, status)
      do s = 1, 2
         lo = (s - 1)*5 + 1
         hi = s*5
         call correct_d_shell_moments(p(lo:hi, lo:hi), m1(lo:hi, lo:hi), v(lo:hi, lo:hi), q0(s), c(s), q2(s), &
            cs(s), q2s(s), ts(s), status)
      end do
      worst = max(maxval(abs(cp - cs)), maxval(abs(q2p - q2s)), maxval(abs(t - ts)))
      call check('(a) block-diagonal spinor kernel vs correct_d_shell_moments', worst, 1.0e-14_rp)

      ! control: V with its spin blocks swapped (a wrong pairing of V with P)
      vsw(1:5, 1:5) = v(6:10, 6:10)
      vsw(6:10, 6:10) = v(1:5, 1:5)
      vsw(1:5, 6:10) = (0.0_rp, 0.0_rp)
      vsw(6:10, 1:5) = (0.0_rp, 0.0_rp)
      call correct_d_shell_moments_spinor(p, m1, vsw, q0, c, q2, cp2, q2p2, t2, status)
      worst_swap = max(maxval(abs(cp2 - cs)), maxval(abs(q2p2 - q2s)))
      call check('(a) control, V spin blocks swapped', worst_swap, 1.0e-14_rp, .true.)

      ! control: dropping the spin-off-diagonal blocks of a spinor V changes the result
      v(1:5, 6:10) = 0.3_rp*rand_hermitian(5, 1.0_rp)
      v(6:10, 1:5) = conjg(transpose(v(1:5, 6:10)))
      h = h0 + v
      call moments_of(h, 6, p, m1, q0, c, q2)
      call correct_d_shell_moments_spinor(p, m1, v, q0, c, q2, cp, q2p, t, status)
      v(1:5, 6:10) = (0.0_rp, 0.0_rp)
      v(6:10, 1:5) = (0.0_rp, 0.0_rp)
      call correct_d_shell_moments_spinor(p, m1, v, q0, c, q2, cp2, q2p2, t2, status)
      worst_off = max(maxval(abs(cp - cp2)), maxval(abs(q2p - q2p2)))
      call check('(a) control, spin-off-diagonal V blocks dropped', worst_off, 1.0e-14_rp, .true.)
   end subroutine test_collinear_reduction

   subroutine test_q2_nonnegative()
      complex(rp) :: h0(10, 10), v(10, 10), h(10, 10), p(10, 10), m1(10, 10), vp(10, 10), vm1c(10, 10), vpv(10, 10)
      real(rp) :: q0(2), c(2), q2(2), cp(2), q2p(2), t(2), q2_flipped, min_q2p, min_flipped
      integer :: trial, s, k, status

      min_q2p = huge(1.0_rp)
      min_flipped = huge(1.0_rp)
      do trial = 1, 40
         h0 = rand_hermitian(10, 1.5_rp)
         v = rand_hermitian(10, 0.8_rp)
         h = h0 + v
         call moments_of(h, 3 + mod(trial, 5), p, m1, q0, c, q2)
         call correct_d_shell_moments_spinor(p, m1, v, q0, c, q2, cp, q2p, t, status)
         min_q2p = min(min_q2p, minval(q2p))
         ! control: the positive Tr(V Pi_s V P) term dropped
         vp = matmul(v, p)
         vpv = matmul(vp, v)
         do s = 1, 2
            vm1c = matmul(v, m1 - c(s)*p)
            q2_flipped = q2(s)
            do k = (s - 1)*5 + 1, s*5
               q2_flipped = q2_flipped - 2.0_rp*real(vm1c(k, k), rp)
            end do
            q2_flipped = q2_flipped - t(s)*t(s)/q0(s)
            min_flipped = min(min_flipped, q2_flipped)
         end do
      end do
      call check('(c) q2 prime >= 0 on 40 random physical states (minimum, negated)', max(0.0_rp, -min_q2p), 1.0e-12_rp)
      call check('(c) control, Tr(V Pi_s V P) term dropped goes negative (minimum, negated)', max(0.0_rp, -min_flipped), 1.0e-12_rp, .true.)
   end subroutine test_q2_nonnegative

   subroutine test_r05_off_diagonal_blocks()
      complex(rp) :: pin(1, 1, 10, 10), pout(1, 1, 10, 10), vin(1, 1, 10, 10), vout(1, 1, 10, 10)
      complex(rp) :: pin_c(1, 2, 5, 5), pout_c(1, 2, 5, 5), vin_c(1, 2, 5, 5), vout_c(1, 2, 5, 5)
      type(mix) :: mixer_s, mixer_c
      real(rp) :: rp_abs, rp_rel, rv_abs, rv_rel, ah_p, ah_v
      logical :: physical, spinor_converged, collinear_converged, ever_spinor, ever_collinear
      integer :: step, worst

      ! a physical P and a V; each step changes only the spin-off-diagonal blocks
      pin(1, 1, :, :) = (0.0_rp, 0.0_rp)
      vin(1, 1, :, :) = (0.0_rp, 0.0_rp)
      do step = 1, 10
         pin(1, 1, step, step) = 0.2_rp + 0.05_rp*real(step, rp)/10.0_rp
      end do
      vin(1, 1, :, :) = pin(1, 1, :, :)*0.5_rp
      mixer_s%fixed_u_stable_updates = 0
      mixer_c%fixed_u_stable_updates = 0
      ever_spinor = .false.
      ever_collinear = .false.
      do step = 1, 6
         pout = pin
         vout = vin
         pout(1, 1, 1, 6) = cmplx(0.01_rp*real(step, rp), 0.004_rp, rp)
         pout(1, 1, 6, 1) = conjg(pout(1, 1, 1, 6))
         vout(1, 1, 2, 7) = cmplx(0.02_rp, -0.01_rp*real(step, rp), rp)
         vout(1, 1, 7, 2) = conjg(vout(1, 1, 2, 7))
         call fixed_u_residual_metrics(pin, pout, vin, vout, rp_abs, rp_rel, rv_abs, rv_rel, physical, ah_p, ah_v, worst)
         call mixer_s%advance_fixed_u_convergence(rp_abs, rv_abs, physical, spinor_converged, 5.0e-7_rp, 1.0e-7_rp)
         ever_spinor = ever_spinor .or. spinor_converged
         ! the collinear layout of the same step sees only the two diagonal 5x5 blocks
         pin_c(1, 1, :, :) = pin(1, 1, 1:5, 1:5)
         pin_c(1, 2, :, :) = pin(1, 1, 6:10, 6:10)
         pout_c(1, 1, :, :) = pout(1, 1, 1:5, 1:5)
         pout_c(1, 2, :, :) = pout(1, 1, 6:10, 6:10)
         vin_c(1, 1, :, :) = vin(1, 1, 1:5, 1:5)
         vin_c(1, 2, :, :) = vin(1, 1, 6:10, 6:10)
         vout_c(1, 1, :, :) = vout(1, 1, 1:5, 1:5)
         vout_c(1, 2, :, :) = vout(1, 1, 6:10, 6:10)
         call fixed_u_residual_metrics(pin_c, pout_c, vin_c, vout_c, rp_abs, rp_rel, rv_abs, rv_rel, physical, ah_p, ah_v, worst)
         call mixer_c%advance_fixed_u_convergence(rp_abs, rv_abs, physical, collinear_converged, 5.0e-7_rp, 1.0e-7_rp)
         ever_collinear = ever_collinear .or. collinear_converged
      end do
      call check('R05 spinor: sequence changing only off-diagonal blocks is declared converged (0 = never)', &
         merge(1.0_rp, 0.0_rp, ever_spinor), 0.5_rp)
      call check('R05 control, diagonal blocks only: the same sequence is declared converged (1 = yes)', &
         merge(1.0_rp, 0.0_rp, ever_collinear), 0.5_rp, .true.)
   end subroutine test_r05_off_diagonal_blocks

end program test_hubbard_spinor_moments
