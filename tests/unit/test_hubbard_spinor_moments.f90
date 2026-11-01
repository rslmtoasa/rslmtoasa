! Spinor local-moment kernel (correct_d_shell_moments_spinor) and the spinor R05 residual.
!
! (a) block-diagonal inputs reproduce the collinear kernel; (c) q2' >= 0 for states built from the exact
! eigenpairs of H0 + V (random Hermitian, with spin-off-diagonal V); R05: a sequence whose only change is
! in the spin-off-diagonal blocks of P and V is not declared converged; (f) the spin-frame axis n: n = z
! against the S4 z-block loop (bitwise), n against the z kernel on spin-rotated inputs, global covariance,
! and a +n / -n site pair; (g) the d-shell orbital moment Tr(P L.n) against the closed form for a|m=2> + b|m=1>. The finite-H comparison with the
! direct definition of c' and q2' is in projector_spinor_frozen_h.f90. Every control must exceed its tolerance.
program test_hubbard_spinor_moments
   use precision_mod, only: rp
   use math_mod, only: correct_d_shell_moments, correct_d_shell_moments_spinor, init_math_operators, spinor_orbital_moment
   use basis_mod, only: basis_init
   use mix_mod, only: mix, fixed_u_residual_metrics
   implicit none

   integer(8) :: seed = 20261027_8
   logical :: failed = .false.
   real(rp), parameter :: zaxis(3) = [0.0_rp, 0.0_rp, 1.0_rp]
   external :: zheev

   call basis_init(2)
   call init_math_operators()
   call test_collinear_reduction()
   call test_q2_nonnegative()
   call test_r05_off_diagonal_blocks()
   call test_frame_rotation()
   call test_orbital_moment()

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
      call correct_d_shell_moments_spinor(p, m1, v, zaxis, q0, c, q2, cp, q2p, t, status)
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
      call correct_d_shell_moments_spinor(p, m1, vsw, zaxis, q0, c, q2, cp2, q2p2, t2, status)
      worst_swap = max(maxval(abs(cp2 - cs)), maxval(abs(q2p2 - q2s)))
      call check('(a) control, V spin blocks swapped', worst_swap, 1.0e-14_rp, .true.)

      ! control: dropping the spin-off-diagonal blocks of a spinor V changes the result
      v(1:5, 6:10) = 0.3_rp*rand_hermitian(5, 1.0_rp)
      v(6:10, 1:5) = conjg(transpose(v(1:5, 6:10)))
      h = h0 + v
      call moments_of(h, 6, p, m1, q0, c, q2)
      call correct_d_shell_moments_spinor(p, m1, v, zaxis, q0, c, q2, cp, q2p, t, status)
      v(1:5, 6:10) = (0.0_rp, 0.0_rp)
      v(6:10, 1:5) = (0.0_rp, 0.0_rp)
      call correct_d_shell_moments_spinor(p, m1, v, zaxis, q0, c, q2, cp2, q2p2, t2, status)
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
         call correct_d_shell_moments_spinor(p, m1, v, zaxis, q0, c, q2, cp, q2p, t, status)
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

   !> 2x2 spin rotation taking the unit vector n to z, from the eigenvectors of n.sigma
   !> (rows <n+| and <n-|, |n+> = (cos(t/2), e^{i f} sin(t/2))), as a 10x10 spin x orbital matrix.
   function rotation_to_z(n) result(u)
      real(rp), intent(in) :: n(3)
      complex(rp) :: u(10, 10), r2(2, 2)
      real(rp) :: theta, phi

      theta = acos(n(3))
      phi = atan2(n(2), n(1))
      r2(1, :) = [cmplx(cos(0.5_rp*theta), 0.0_rp, rp), exp(cmplx(0.0_rp, -phi, rp))*sin(0.5_rp*theta)]
      r2(2, :) = [-exp(cmplx(0.0_rp, phi, rp))*sin(0.5_rp*theta), cmplx(cos(0.5_rp*theta), 0.0_rp, rp)]
      u = spin_matrix(r2)
   end function rotation_to_z

   !> W = exp(-i alpha m.sigma/2) x 1_orb; it rotates spin vectors by +alpha about the unit vector m.
   function spin_rotation(m, alpha) result(u)
      real(rp), intent(in) :: m(3), alpha
      complex(rp) :: u(10, 10), r2(2, 2)
      complex(rp), parameter :: im = (0.0_rp, 1.0_rp)

      r2(1, 1) = cos(0.5_rp*alpha) - im*sin(0.5_rp*alpha)*m(3)
      r2(1, 2) = -im*sin(0.5_rp*alpha)*cmplx(m(1), -m(2), rp)
      r2(2, 1) = -im*sin(0.5_rp*alpha)*cmplx(m(1), m(2), rp)
      r2(2, 2) = cos(0.5_rp*alpha) + im*sin(0.5_rp*alpha)*m(3)
      u = spin_matrix(r2)
   end function spin_rotation

   function spin_matrix(r2) result(u)
      complex(rp), intent(in) :: r2(2, 2)
      complex(rp) :: u(10, 10)
      integer :: i

      u = (0.0_rp, 0.0_rp)
      do i = 1, 5
         u(i, i) = r2(1, 1)
         u(i, i + 5) = r2(1, 2)
         u(i + 5, i) = r2(2, 1)
         u(i + 5, i + 5) = r2(2, 2)
      end do
   end function spin_matrix

   !> Rodrigues: rotation of v by alpha about the unit vector m.
   function rotate_vector(m, alpha, v) result(w)
      real(rp), intent(in) :: m(3), alpha, v(3)
      real(rp) :: w(3), mxv(3)

      mxv = [m(2)*v(3) - m(3)*v(2), m(3)*v(1) - m(1)*v(3), m(1)*v(2) - m(2)*v(1)]
      w = v*cos(alpha) + mxv*sin(alpha) + m*dot_product(m, v)*(1.0_rp - cos(alpha))
   end function rotate_vector

   function random_unit() result(n)
      real(rp) :: n(3)

      n = [uniform() - 0.5_rp, uniform() - 0.5_rp, uniform() - 0.5_rp]
      n = n/norm2(n)
   end function random_unit

   !> The S4 z-frame kernel loop (global-z diagonal blocks), kept as the regression reference for n = z.
   subroutine s4_z_kernel(p, m1, vin, q0, c, q2, cprime, q2prime, t)
      complex(rp), intent(in) :: p(10, 10), m1(10, 10), vin(10, 10)
      real(rp), intent(in) :: q0(2), c(2), q2(2)
      real(rp), intent(out) :: cprime(2), q2prime(2), t(2)
      complex(rp) :: vp(10, 10), vpv(10, 10), vm1c(10, 10)
      integer :: s, k

      vp = matmul(vin, p)
      vpv = matmul(vp, vin)
      do s = 1, 2
         vm1c = matmul(vin, m1 - c(s)*p)
         t(s) = 0.0_rp
         q2prime(s) = q2(s)
         do k = (s - 1)*5 + 1, s*5
            t(s) = t(s) + real(vp(k, k), rp)
            q2prime(s) = q2prime(s) - 2.0_rp*real(vm1c(k, k), rp) + real(vpv(k, k), rp)
         end do
         cprime(s) = c(s) - t(s)/q0(s)
         q2prime(s) = q2prime(s) - t(s)*t(s)/q0(s)
      end do
   end subroutine s4_z_kernel

   subroutine test_frame_rotation()
      complex(rp) :: h(10, 10), v(10, 10), pz(10, 10), m1z(10, 10), vz(10, 10), u(10, 10), w(10, 10)
      complex(rp) :: p(10, 10), m1(10, 10), pw(10, 10), m1w(10, 10), vw(10, 10)
      real(rp) :: q0(2), c(2), q2(2), cp(2), q2p(2), t(2), cpr(2), q2pr(2), tr(2), cpw(2), q2pw(2), tw(2)
      real(rp) :: n(3), m(3), nw(3), worst_z, worst_rot, worst_cov, worst_ctrl, worst_afm, worst_afm_ctrl, alpha
      integer :: trial, status, st

      ! n = z: bitwise against the S4 loop, on general (spin-mixed) inputs
      worst_z = 0.0_rp
      do trial = 1, 10
         h = rand_hermitian(10, 1.5_rp)
         v = rand_hermitian(10, 0.6_rp)
         call moments_of(h + v, 3 + mod(trial, 5), p, m1, q0, c, q2)
         call correct_d_shell_moments_spinor(p, m1, v, zaxis, q0, c, q2, cp, q2p, t, status)
         call s4_z_kernel(p, m1, v, q0, c, q2, cpr, q2pr, tr)
         worst_z = max(worst_z, maxval(abs(cp - cpr)), maxval(abs(q2p - q2pr)), maxval(abs(t - tr)))
      end do
      call check('(f) n = z vs the S4 z-block loop (bitwise)', worst_z, 0.0_rp)

      worst_rot = 0.0_rp
      worst_cov = 0.0_rp
      worst_ctrl = 0.0_rp
      worst_afm = 0.0_rp
      worst_afm_ctrl = 0.0_rp
      do trial = 1, 20
         n = random_unit()
         u = rotation_to_z(n)
         ! a physical state in the n frame: q0, c, q2 along n are the z-block moments of the rotated problem
         h = rand_hermitian(10, 1.5_rp)
         v = rand_hermitian(10, 0.6_rp)
         call moments_of(h + v, 3 + mod(trial, 5), pz, m1z, q0, c, q2)
         vz = v
         p = matmul(conjg(transpose(u)), matmul(pz, u))
         m1 = matmul(conjg(transpose(u)), matmul(m1z, u))
         v = matmul(conjg(transpose(u)), matmul(vz, u))
         call correct_d_shell_moments_spinor(p, m1, v, n, q0, c, q2, cp, q2p, t, status)
         if (status /= 0) failed = .true.
         call correct_d_shell_moments_spinor(pz, m1z, vz, zaxis, q0, c, q2, cpr, q2pr, tr, status)
         worst_rot = max(worst_rot, maxval(abs(cp - cpr)), maxval(abs(q2p - q2pr)), maxval(abs(t - tr)))
         ! control: the z projectors on the tilted state
         call correct_d_shell_moments_spinor(p, m1, v, zaxis, q0, c, q2, cpw, q2pw, tw, status)
         worst_ctrl = max(worst_ctrl, maxval(abs(cpw - cpr)), maxval(abs(q2pw - q2pr)))
         ! global covariance: P, M1, V and n rotated together by W
         m = random_unit()
         alpha = 6.0_rp*uniform()
         w = spin_rotation(m, alpha)
         pw = matmul(w, matmul(p, conjg(transpose(w))))
         m1w = matmul(w, matmul(m1, conjg(transpose(w))))
         vw = matmul(w, matmul(v, conjg(transpose(w))))
         nw = rotate_vector(m, alpha, n)
         call correct_d_shell_moments_spinor(pw, m1w, vw, nw, q0, c, q2, cpw, q2pw, tw, status)
         worst_cov = max(worst_cov, maxval(abs(cpw - cp)), maxval(abs(q2pw - q2p)), maxval(abs(tw - t)))
         ! +n / -n site pair: the second site is the first rotated by pi about an axis normal to n
         m = [n(2), -n(1), 0.0_rp]
         m = m/norm2(m)
         w = spin_rotation(m, acos(-1.0_rp))
         pw = matmul(w, matmul(p, conjg(transpose(w))))
         m1w = matmul(w, matmul(m1, conjg(transpose(w))))
         vw = matmul(w, matmul(v, conjg(transpose(w))))
         call correct_d_shell_moments_spinor(pw, m1w, vw, -n, q0, c, q2, cpw, q2pw, tw, status)
         worst_afm = max(worst_afm, maxval(abs(cpw - cp)), maxval(abs(q2pw - q2p)), maxval(abs(tw - t)))
         ! control (pre-S6b behaviour for -z): the -n site evaluated with the +n projectors
         call correct_d_shell_moments_spinor(pw, m1w, vw, n, q0, c, q2, cpw, q2pw, tw, status)
         worst_afm_ctrl = max(worst_afm_ctrl, maxval(abs(cpw - cp)), maxval(abs(q2pw - q2p)))
      end do
      call check('(f) random n vs the z kernel on spin-rotated P, M1, V', worst_rot, 1.0e-13_rp)
      call check('(f) control, z projectors with a tilted n', worst_ctrl, 1.0e-13_rp, .true.)
      call check('(f) global covariance: P, M1, V and n rotated together', worst_cov, 1.0e-13_rp)
      call check('(f) +n / -n site pair: the -n site equals the +n site', worst_afm, 1.0e-13_rp)
      call check('(f) control, the -n site with the +n projectors', worst_afm_ctrl, 1.0e-13_rp, .true.)
   end subroutine test_frame_rotation

   !> m_L = Tr(P L.n) against closed forms. The d basis of the code (cart2sph) has L_z = 3 - k on index k and the
   !> neighbouring-state elements of L_x, L_y carry the opposite sign to the Condon-Shortley convention (basis phase
   !> (-1)^m; read off the L matrices, so the sign of <Lx>, <Ly> is not independent of them; the magnitudes are).
   !> For psi = a|k> + b|k+1> (spin up) with m = 3 - k and c = sqrt(l(l+1) - m(m-1))/2 the |m-1>,|m> element of L_x:
   !> <Lx> = -2 c Re(a* b), <Ly> = -2 c Im(a* b), <Lz> = m|a|^2 + (m-1)|b|^2; a spin-down |k=1> adds 2 n.z.
   subroutine test_orbital_moment()
      complex(rp) :: psi(10), chi(10), p(10, 10), a, b, v(10)
      real(rp) :: n(3), expected, worst, worst_t, worst_diag, lvec(3), c
      integer :: trial, k, m

      worst_diag = 0.0_rp
      do k = 1, 5
         v = (0.0_rp, 0.0_rp)
         v(k) = (1.0_rp, 0.0_rp)
         p = spread(v, 2, 10)*spread(conjg(v), 1, 10)
         worst_diag = max(worst_diag, abs(spinor_orbital_moment(p, zaxis) - real(3 - k, rp)), &
                          abs(spinor_orbital_moment(p, [1.0_rp, 0.0_rp, 0.0_rp])), abs(spinor_orbital_moment(p, [0.0_rp, 1.0_rp, 0.0_rp])))
      end do
      call check('(g) orbital moment of |k> along z = 3 - k, zero along x and y', worst_diag, 1.0e-13_rp)

      worst = 0.0_rp
      worst_t = 0.0_rp
      do trial = 1, 8
         k = 1 + mod(trial, 2)
         m = 3 - k
         c = 0.5_rp*sqrt(6.0_rp - real(m*(m - 1), rp))
         a = cmplx(0.6_rp, 0.0_rp, rp)
         b = 0.8_rp*exp(cmplx(0.0_rp, 0.4_rp + 0.9_rp*real(trial, rp), rp))
         psi = (0.0_rp, 0.0_rp)
         psi(k) = a
         psi(k + 1) = b
         chi = (0.0_rp, 0.0_rp)
         chi(6) = (1.0_rp, 0.0_rp)
         p = spread(psi, 2, 10)*spread(conjg(psi), 1, 10) + spread(chi, 2, 10)*spread(conjg(chi), 1, 10)
         n = random_unit()
         lvec = [-2.0_rp*c*real(conjg(a)*b, rp), -2.0_rp*c*aimag(conjg(a)*b), &
                 real(m, rp)*abs(a)**2 + real(m - 1, rp)*abs(b)**2]
         expected = dot_product(n, lvec) + 2.0_rp*n(3)
         worst = max(worst, abs(spinor_orbital_moment(p, n) - expected))
         ! control: the transposed P flips the sign of <Ly>
         worst_t = max(worst_t, abs(spinor_orbital_moment(transpose(p), n) - expected))
      end do
      call check('(g) orbital moment Tr(P L.n) vs the closed form, 8 states', worst, 1.0e-13_rp)
      call check('(g) control, transposed P', worst_t, 1.0e-13_rp, .true.)
   end subroutine test_orbital_moment

end program test_hubbard_spinor_moments
