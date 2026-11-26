!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
! Catches in spin_response_rs: a spin-down offset other than spin_off (lmax = 3
! here, spin_off = 16), a d block taken from the wrong rows, a missing transpose
! in the trace, swapped spin roles of the two columns, and a missing conjugate or
! wrong sign in chi0_from_traces; in d_moment_trace a wrong orbital range or spin
! sign; in type_sums a seed taken from a non-self pair or from a later self pair,
! pairs of another seed summed in, grouping by the wrong atom, or the imaginary
! part included; in juelich_vectors the spin blocks, orbitals or p = 1/ppar^2
! mixed up; in juelich_green and juelich_moments a wrong or missing correction
! term (H_ab, delta_ab, z dependence) or a wrong moment shift. The Juelich
! references are spectral sums over the eigenpairs of a small Hermitian H, a
! route independent of the resolvent identity and the Chebyshev recurrence the
! kernels use. In type_sums_q a wrong phase sign, pairs of another seed or a
! wrong pair-type grouping. In chi0_dynamic_traces a wrong sign of the energy shift, swapped spin roles, a missing
! conjugate or transpose, or a wrong d block, against the Lehmann sum of the eigenpairs of a two-site model. Plumbing only: no physics is checked here.
!------------------------------------------------------------------------------
program test_spin_response_rs
   use precision_mod, only: rp
   use math_mod, only: pi
   use basis_mod, only: basis_init, nb
   use spin_response_rs_mod, only: chi0_traces, chi0_from_traces, d_moment_trace, type_sums, juelich_vectors, &
                                   juelich_green, juelich_moments, type_sums_q, chi0_dynamic_traces
   implicit none
   complex(rp), allocatable :: gij(:, :, :), gji(:, :, :), dm(:)
   complex(rp) :: t(2, 2), c
   integer :: mu, ijpair(8, 2), seed(2), npair(2)
   real(rp) :: m(2), chi0s(2, 2)
   real(rp), allocatable :: av(:), bv(:)
   logical :: failed

   failed = .false.
   call basis_init(3)
   ! Every element outside the four d spin blocks, and inside them except the entries set below, would add 100 or more
   allocate (gij(nb, nb, 2), gji(nb, nb, 2))
   gij = (100.0_rp, 0.0_rp)
   gji = (100.0_rp, 0.0_rp)
   gij(5:9, 5:9, :) = 0.0_rp      ! up d
   gij(21:25, 21:25, :) = 0.0_rp  ! down d
   gji(5:9, 5:9, :) = 0.0_rp
   gji(21:25, 21:25, :) = 0.0_rp
   ! Column 1: Tr[G^upup_ji G^dndn_ij] = gji(5,6) gij(22,21) = 1 * 2; column 2: Tr[G^dndn_ji G^upup_ij] = gji(22,21) gij(5,6) = 3 * 4
   gji(5, 6, 1) = 1.0_rp
   gij(22, 21, 1) = 2.0_rp
   gji(22, 21, 1) = 3.0_rp
   gij(5, 6, 1) = 4.0_rp
   ! Second index: every set entry doubled, so both traces are four times larger
   gji(:, :, 2) = gji(:, :, 1)
   gij(:, :, 2) = gij(:, :, 1)
   gji(5, 6, 2) = 2.0_rp
   gij(22, 21, 2) = 4.0_rp
   gji(22, 21, 2) = 6.0_rp
   gij(5, 6, 2) = 8.0_rp

   t = chi0_traces(gij, gji)
   call require(all(abs(t(:, 1) - [2.0_rp, 8.0_rp]) < 1.0e-14_rp), 'column 1 = 2, 8')
   call require(all(abs(t(:, 2) - [12.0_rp, 48.0_rp]) < 1.0e-14_rp), 'column 2 = 12, 48')

   ! chi0 = -(x - conj(y))/(2 pi i): x = y = 1 + 2i gives -2/pi; x = 1, y = 0 gives i/(2 pi)
   c = chi0_from_traces((1.0_rp, 2.0_rp), (1.0_rp, 2.0_rp))
   call require(abs(c - cmplx(-2.0_rp/pi, 0.0_rp, rp)) < 1.0e-15_rp, 'x = y = 1 + 2i gives -2/pi')
   c = chi0_from_traces((1.0_rp, 0.0_rp), (0.0_rp, 0.0_rp))
   call require(abs(c - cmplx(0.0_rp, 0.5_rp/pi, rp)) < 1.0e-15_rp, 'x = 1, y = 0 gives i/(2 pi)')

   ! d_moment_trace: d up diagonal 1..5, d down diagonal 0.5 each, every other element 100: 15 - 2.5 = 12.5
   gij(:, :, 1) = (100.0_rp, 0.0_rp)
   do mu = 5, 9
      gij(mu, mu, 1) = real(mu - 4, rp)
      gij(mu + 16, mu + 16, 1) = 0.5_rp
   end do
   dm = d_moment_trace(gij(:, :, 1:1))
   call require(abs(dm(1) - 12.5_rp) < 1.0e-14_rp, 'd_moment_trace = 12.5')

   ! type_sums: atoms 1..6 of types 1, 2, 1, 2, 1, 2. Pair 1 (3,4) is not a self pair, pair 5 (3,3) a later self pair of
   ! type 1: neither may become a seed. Seeds 1 (m 2.0) and 2 (m -1.5). Type 1 rows: pairs 2, 3, 4, 8 (j types 1, 2, 1, 1),
   ! chi0s(1,1) = -5 - 0.5 - 0.125 = -5.625, chi0s(1,2) = -1; type 2 rows: pairs 6, 7, chi0s(2,2) = -4, chi0s(2,1) = -0.25.
   ijpair(:, 1) = [3, 1, 1, 1, 3, 2, 2, 1]
   ijpair(:, 2) = [4, 1, 2, 3, 3, 2, 1, 5]
   call type_sums(ijpair, [1, 2, 1, 2, 1, 2], 2, &
                  [(10.0_rp, 0.0_rp), (-5.0_rp, 0.0_rp), (-1.0_rp, 0.0_rp), (-0.5_rp, 0.0_rp), (-7.0_rp, 0.0_rp), &
                   (-4.0_rp, 0.0_rp), (-0.25_rp, 0.0_rp), (-0.125_rp, 3.0_rp)], &
                  [0.0_rp, 2.0_rp, 0.0_rp, 0.0_rp, 9.0_rp, -1.5_rp, 0.0_rp, 0.0_rp], seed, npair, m, chi0s)
   call require(all(seed == [1, 2]), 'type_sums seed = 1, 2')
   call require(all(npair == [4, 2]), 'type_sums npair = 4, 2')
   call require(all(abs(m - [2.0_rp, -1.5_rp]) < 1.0e-14_rp), 'type_sums m = 2.0, -1.5')
   call require(all(abs(chi0s - reshape([-5.625_rp, -0.25_rp, -1.0_rp, -4.0_rp], [2, 2])) < 1.0e-14_rp), &
                'type_sums chi0s = [-5.625 -1; -0.25 -4]')

   ! juelich_vectors, closed form: up E_nu = 0, ppar = 1, E_F = 1 gives alpha = beta = 1/sqrt(2); down E_nu = 0.5, ppar = 2
   ! gives x = 0.5, p = 0.25, alpha = 0.9375/sqrt(1.0625), beta = 0.125/sqrt(1.0625); s, p, f orbitals alpha = 1, beta = 0
   allocate (av(nb), bv(nb))
   call juelich_vectors([0.0_rp, 0.5_rp], [1.0_rp, 2.0_rp], 1.0_rp, av, bv)
   call require(all(abs(av(5:9) - 1.0_rp/sqrt(2.0_rp)) < 1.0e-15_rp) .and. all(abs(bv(5:9) - 1.0_rp/sqrt(2.0_rp)) < 1.0e-15_rp), &
                'juelich_vectors up d')
   call require(all(abs(av(21:25) - 0.9375_rp/sqrt(1.0625_rp)) < 1.0e-15_rp) .and. &
                all(abs(bv(21:25) - 0.125_rp/sqrt(1.0625_rp)) < 1.0e-15_rp), 'juelich_vectors down d')
   call require(all(av(1:4) == 1.0_rp) .and. all(av(10:20) == 1.0_rp) .and. all(av(26:32) == 1.0_rp) .and. &
                all(bv(1:4) == 0.0_rp) .and. all(bv(10:20) == 0.0_rp) .and. all(bv(26:32) == 0.0_rp), 'juelich_vectors other orbitals')

   call check_juelich_kernels()
   call check_type_sums_q()
   call check_chi0_dynamic()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains
   !> Two sites of two orbitals: H (4x4, Hermitian, spectrum inside (b_win - a_win, b_win + a_win)) from fixed numbers.
   !> Reference G~_ab(z) = sum_n F_a(e_n) F_b(e_n) v_an conj(v_bn)/(z - e_n) and mu~_m = sum_n F_a F_b T_m((e_n - b)/a) v v*.
   subroutine check_juelich_kernels()
      integer, parameter :: n = 4, nz = 3, nm = 8
      real(rp), parameter :: a_win = 2.0_rp, b_win = 0.3_rp
      complex(rp) :: h(n, n), v(n, n), work(4*n), g(n, n, nz), gt(n, n, nz), gref(n, n, nz)
      complex(rp) :: mu(n, n, nm), mt(n, n, nm - 2), mref(n, n, nm - 2)
      real(rp) :: e(n), rwork(3*n), z(nz), alpha(n), beta(n), fa, fb, x, tm
      integer :: info, k, l, q, mm

      h = reshape([(0.31_rp, 0.0_rp), (0.12_rp, -0.05_rp), (-0.07_rp, 0.02_rp), (0.04_rp, 0.09_rp), &
                   (0.12_rp, 0.05_rp), (-0.22_rp, 0.0_rp), (0.06_rp, -0.11_rp), (0.08_rp, 0.0_rp), &
                   (-0.07_rp, -0.02_rp), (0.06_rp, 0.11_rp), (0.45_rp, 0.0_rp), (-0.13_rp, 0.03_rp), &
                   (0.04_rp, -0.09_rp), (0.08_rp, 0.0_rp), (-0.13_rp, -0.03_rp), (0.02_rp, 0.0_rp)], [n, n])
      alpha = [1.04_rp, 0.97_rp, 1.10_rp, 0.92_rp]
      beta = [0.27_rp, -0.15_rp, 0.33_rp, 0.08_rp]
      z = [-0.3_rp, 0.05_rp, 0.4_rp]
      v = h
      call zheev('V', 'U', n, v, n, e, work, size(work), rwork, info)
      if (info /= 0) error stop 'zheev failed'
      ! exact G and spectral references
      g = 0.0_rp
      gref = 0.0_rp
      mu = 0.0_rp
      mref = 0.0_rp
      do q = 1, n
         x = (e(q) - b_win)/a_win
         do l = 1, n
            do k = 1, n
               fa = alpha(k) + beta(k)*e(q)
               fb = alpha(l) + beta(l)*e(q)
               g(k, l, :) = g(k, l, :) + v(k, q)*conjg(v(l, q))/(z - e(q))
               gref(k, l, :) = gref(k, l, :) + fa*fb*v(k, q)*conjg(v(l, q))/(z - e(q))
               do mm = 0, nm - 1
                  tm = cos(mm*acos(x))
                  mu(k, l, mm + 1) = mu(k, l, mm + 1) + tm*v(k, q)*conjg(v(l, q))
                  if (mm <= nm - 3) mref(k, l, mm + 1) = mref(k, l, mm + 1) + fa*fb*tm*v(k, q)*conjg(v(l, q))
               end do
            end do
         end do
      end do
      ! self pair: the whole 4x4 block, with delta_ab
      call juelich_green(g, h, z, alpha, beta, alpha, beta, .true., gt)
      call require(maxval(abs(gt - gref)) < 1.0e-13_rp, 'juelich_green self pair = spectral sum')
      ! pair between site 1 (orbitals 1, 2) and site 2 (orbitals 3, 4): no delta term
      call juelich_green(g(1:2, 3:4, :), h(1:2, 3:4), z, alpha(1:2), beta(1:2), alpha(3:4), beta(3:4), .false., gt(1:2, 1:2, :))
      call require(maxval(abs(gt(1:2, 1:2, :) - gref(1:2, 3:4, :))) < 1.0e-13_rp, 'juelich_green inter-site block = spectral sum')
      call juelich_moments(mu, a_win, b_win, alpha, beta, alpha, beta, mt)
      call require(maxval(abs(mt - mref)) < 1.0e-13_rp, 'juelich_moments = spectral sum, m = 0..5')
   end subroutine check_juelich_kernels

   !> One type, seed 1. Pairs (1,1) -10 at r = 0, (1,2) -1 at r = (1/2, 0, 0), (1,3) -2 at r = (0, 1/4, 0), and (3,2) 100 of
   !> another first atom. q = (1,0,0): phases -1, 1 give -10 + 1 - 2 = -11; q = (0,1,0): phases 1, exp(+i pi/2) = i give -11 - 2i.
   subroutine check_type_sums_q()
      complex(rp) :: cq(1, 1, 2)
      integer :: pairs(4, 2)

      pairs(:, 1) = [1, 1, 1, 3]
      pairs(:, 2) = [1, 2, 3, 2]
      cq = type_sums_q(pairs, [1, 1, 1], [1], [(-10.0_rp, 0.0_rp), (-1.0_rp, 0.0_rp), (-2.0_rp, 0.0_rp), (100.0_rp, 0.0_rp)], &
                       reshape([0.0_rp, 0.0_rp, 0.0_rp, 0.5_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.25_rp, 0.0_rp, 0.3_rp, 0.3_rp, 0.3_rp], [3, 4]), &
                       reshape([1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp], [3, 2]))
      call require(abs(cq(1, 1, 1) - (-11.0_rp, 0.0_rp)) < 1.0e-13_rp, 'type_sums_q q = (1,0,0) gives -11')
      call require(abs(cq(1, 1, 2) - (-11.0_rp, -2.0_rp)) < 1.0e-13_rp, 'type_sums_q q = (0,1,0) gives -11 - 2i')
   end subroutine check_type_sums_q

   !> Two sites of five d orbitals, different Hermitian H for spin up and down (10x10, fixed numbers). Exact G_ij, G_ji on a
   !> uniform mesh E_n = -1 + 0.05 (n - 1), n = 1..41, at z = E_n + i eta. Reference: the Lehmann form evaluated from the
   !> eigenpairs, with A_a(E) = (eta/pi)/((E - e_a)^2 + eta^2) and weights w_ab = sum_{mu nu} vu(j mu, a) conj(vu(i nu, a))
   !> vd(i nu, b) conj(vd(j mu, b)):
   !>   term 1 = sum_n f_n sum_ab w_ab A^up_a(E_n) /(E_(n+k) + i eta - e^dn_b),
   !>   term 2 = sum_n f_n sum_ab w_ab /(E_(n-k) - i eta - e^up_a) A^dn_b(E_n),
   !> f_n = 1 for E_n <= 0.2, terms with n +- k outside the mesh dropped. A wrong sign of the shift, a swapped spin role, a
   !> missing conjugate or transpose, or a wrong d block changes the sum for k = 3 and k = -2 (the reference values differ).
   subroutine check_chi0_dynamic()
      integer, parameter :: ns = 10, nn = 41
      real(rp), parameter :: de = 0.05_rp, eta = 0.03_rp, ef = 0.2_rp
      complex(rp) :: h(ns, ns, 2), v(ns, ns, 2), work(4*ns), g(nb, nb, nn), gt(nb, nb, nn), z, w, kern, ref, xs(2)
      complex(rp) :: tr(nn, 2), td(nn, 2)
      real(rp) :: e(ns, 2), rwork(3*ns), en(nn), a_up, a_dn
      integer :: s, p, q, info, n, a, b, mu, nu, ks(2), ik, k, mm

      do s = 1, 2
         do q = 1, ns
            do p = 1, q - 1
               h(p, q, s) = cmplx(0.3_rp*sin(0.7_rp*p + 1.3_rp*q + s), 0.2_rp*cos(0.4_rp*p - 0.9_rp*q + s), rp)
               h(q, p, s) = conjg(h(p, q, s))
            end do
            h(q, q, s) = cmplx(0.12_rp*q - 0.7_rp + 0.15_rp*(s - 1), 0.0_rp, rp)
         end do
         v(:, :, s) = h(:, :, s)
         call zheev('V', 'U', ns, v(:, :, s), ns, e(:, s), work, size(work), rwork, info)
         if (info /= 0) error stop 'zheev failed'
      end do
      en = [(-1.0_rp + de*(n - 1), n=1, nn)]
      ! G_ij = rows of site i (1:5), columns of site j (6:10); G_ji the transposed roles; every other element 100
      g = (100.0_rp, 0.0_rp)
      gt = (100.0_rp, 0.0_rp)
      do s = 1, 2
         do n = 1, nn
            z = cmplx(en(n), eta, rp)
            do mu = 1, 5
               do nu = 1, 5
                  g(4 + nu + 16*(s - 1), 4 + mu + 16*(s - 1), n) = &
                     sum(v(nu, :, s)*conjg(v(5 + mu, :, s))/(z - e(:, s)))
                  gt(4 + mu + 16*(s - 1), 4 + nu + 16*(s - 1), n) = &
                     sum(v(5 + mu, :, s)*conjg(v(nu, :, s))/(z - e(:, s)))
               end do
            end do
         end do
      end do

      ks = [3, -2]
      do ik = 1, 2
         k = ks(ik)
         ! kernel: gij = g, gji = gt
         td = chi0_dynamic_traces(g, gt, k)
         kern = (0.0_rp, 0.0_rp)
         ref = (0.0_rp, 0.0_rp)
         do n = 1, nn
            if (en(n) > ef) cycle
            kern = kern + td(n, 1) + td(n, 2)
            do a = 1, ns
               do b = 1, ns
                  w = (0.0_rp, 0.0_rp)
                  do mu = 1, 5
                     do nu = 1, 5
                        w = w + v(5 + mu, a, 1)*conjg(v(nu, a, 1))*v(nu, b, 2)*conjg(v(5 + mu, b, 2))
                     end do
                  end do
                  a_up = eta/pi/((en(n) - e(a, 1))**2 + eta**2)
                  a_dn = eta/pi/((en(n) - e(b, 2))**2 + eta**2)
                  mm = n + k
                  if (mm >= 1 .and. mm <= nn) ref = ref + w*a_up/cmplx(en(mm) - e(b, 2), eta, rp)
                  mm = n - k
                  if (mm >= 1 .and. mm <= nn) ref = ref + w/cmplx(en(mm) - e(a, 1), -eta, rp)*a_dn
               end do
            end do
         end do
         xs(ik) = ref
         call require(abs(ref) > 1.0e-3_rp, 'chi0_dynamic reference is not small')
         call require(abs(kern - ref) < 1.0e-11_rp*abs(ref), 'chi0_dynamic_traces = Lehmann sum')
      end do
      call require(abs(xs(1) - xs(2)) > 1.0e-3_rp*abs(xs(1)), 'chi0_dynamic reference depends on the sign of k')

      ! k = 0 against the static traces (by construction: the algebra of A = (G - G^A)/(-2 pi i) reduces to chi0_from_traces)
      td = chi0_dynamic_traces(g, gt, 0)
      tr = chi0_traces(g, gt)
      kern = sum(td(:, 1) + td(:, 2))
      ref = chi0_from_traces(sum(tr(:, 1)), sum(tr(:, 2)))
      call require(abs(kern - ref) < 1.0e-11_rp*abs(ref), 'chi0_dynamic_traces at k = 0 = static (by construction)')
   end subroutine check_chi0_dynamic

   subroutine require(ok, label)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label
      if (.not. ok) then
         failed = .true.
         write (*, '(2a)') 'FAIL: ', label
      end if
   end subroutine require
end program test_spin_response_rs
