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
! wrong pair-type grouping. In cluster_sites a non-self seed taken as a site, a wrong
! site order, a missing or unconjugated fill of a pair listed once, a listed pair
! overwritten by the conjugate of its partner, a non-zero entry for an absent pair or
! the wrong self pair for m; in cluster_frequencies a wrong U or sign of Omega.
! In j0_local_traces (3 sites x 2 orbitals, H^s = C^s + W^s S W^s with a Hermitian S, complex z,
! reference the direct sum over all sites of Tr[de_i G^up_ij de_j G^dn_ji] with G from a linear solve): r = W_dn/W_up
! inverted or applied to the wrong spin, the a/b (z-dependence) split of the vertex swapped, a wrong spin_off, a dropped
! orbital.
! In chi0_dynamic_traces a wrong sign of the energy shift, swapped spin roles, a missing
! conjugate or transpose, or a wrong d block, against the Lehmann sum of the eigenpairs of a two-site model. Plumbing only: no physics is checked here.
!------------------------------------------------------------------------------
program test_spin_response_rs
   use precision_mod, only: rp
   use math_mod, only: pi
   use basis_mod, only: basis_init, nb, spin_off
   use spin_response_rs_mod, only: chi0_traces, chi0_from_traces, d_moment_trace, type_sums, juelich_vectors, &
                                   juelich_green, juelich_moments, type_sums_q, cluster_sites, cluster_frequencies, &
                                   j0_local_traces, chi0_dynamic_traces, lorentz_smooth, dynamic_propagator
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
   call check_cluster()
   call check_j0_local()
   call check_chi0_dynamic()
   call check_lorentz_smooth()

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

   !> Pairs (5,5), (9,5), (2,2), (5,2), (7,7), (7,5), (2,5). Atom 9 only has a non-self pair: not a site. Sites 5, 2, 7 in
   !> list order, m = 2, -1.5, 1 (md of the self pairs; the md of the other pairs, 7, must not enter).
   !> (5,2) = (-0.5, 0.25) listed, (2,5) = (-0.6, -0.3) listed too: both used as given (not the conjugate of each other).
   !> (7,5) = (-0.1, 0.2) only: (5,7) = conj = (-0.1, -0.2). (2,7) and (7,2) absent: 0.
   subroutine check_cluster()
      integer, allocatable :: site(:)
      complex(rp), allocatable :: c(:, :)
      real(rp), allocatable :: mc(:)
      complex(rp) :: ref(3, 3)
      integer :: pairs(7, 2)
      real(rp) :: u(2), om(2), omf(2), imag_max

      pairs(:, 1) = [5, 9, 2, 5, 7, 7, 2]
      pairs(:, 2) = [5, 5, 2, 2, 7, 5, 5]
      call cluster_sites(pairs, [(-3.0_rp, 0.0_rp), (-1.0_rp, 0.5_rp), (-4.0_rp, 0.0_rp), (-0.5_rp, 0.25_rp), &
                                 (-2.0_rp, 0.0_rp), (-0.1_rp, 0.2_rp), (-0.6_rp, -0.3_rp)], &
                         [2.0_rp, 7.0_rp, -1.5_rp, 7.0_rp, 1.0_rp, 7.0_rp, 7.0_rp], site, c, mc)
      call require(size(site) == 3, 'cluster_sites: 3 sites')
      if (size(site) /= 3) return
      call require(all(site == [5, 2, 7]), 'cluster_sites: sites 5, 2, 7')
      call require(all(abs(mc - [2.0_rp, -1.5_rp, 1.0_rp]) < 1.0e-14_rp), 'cluster_sites: m = 2, -1.5, 1')
      ref = reshape([(-3.0_rp, 0.0_rp), (-0.6_rp, -0.3_rp), (-0.1_rp, 0.2_rp), &
                     (-0.5_rp, 0.25_rp), (-4.0_rp, 0.0_rp), (0.0_rp, 0.0_rp), &
                     (-0.1_rp, -0.2_rp), (0.0_rp, 0.0_rp), (-2.0_rp, 0.0_rp)], [3, 3])
      call require(maxval(abs(c - ref)) < 1.0e-14_rp, 'cluster_sites: chi0 matrix (listed, conjugate fill, zero)')

      ! Two identical sites, chi0 = [[a, b], [b, a]] = [[-3, -1], [-1, -3]], m = 2. Goldstone: (a + b) U m = -m, U = -1/(a + b) = 0.25.
      ! Eigenvectors (1, 1) and (1, -1) of chi0 with eigenvalues a + b, a - b; Omega = -m (1/lambda + U):
      ! symmetric 0, antisymmetric -m (1/(a - b) - 1/(a + b)) = -2 m b/(a^2 - b^2) = 4/8 = 0.5.
      ! Force theorem, m U (1 + U lambda): symmetric 0, antisymmetric m U (1 + U (a - b)) = 2 * 0.25 * (1 - 0.5) = 0.25.
      call cluster_frequencies(reshape([(-3.0_rp, 0.0_rp), (-1.0_rp, 0.0_rp), (-1.0_rp, 0.0_rp), (-3.0_rp, 0.0_rp)], [2, 2]), &
                               [2.0_rp, 2.0_rp], u, om, imag_max, omf)
      call require(all(abs(u - 0.25_rp) < 1.0e-13_rp), 'cluster_frequencies: U = 0.25')
      call require(all(abs(om - [0.0_rp, 0.5_rp]) < 1.0e-12_rp), 'cluster_frequencies: Omega = 0, 0.5')
      call require(all(abs(omf - [0.0_rp, 0.25_rp]) < 1.0e-12_rp), 'cluster_frequencies: force = 0, 0.25')
      call require(imag_max < 1.0e-12_rp, 'cluster_frequencies: max |Im| = 0')
   end subroutine check_cluster
   !> H^s = C^s + W^s S W^s on 3 sites of 2 orbitals (S Hermitian, spin independent), G^s = (z - H^s)^-1 by zgesv at two complex z.
   !> Reference: sum_j Tr[de_i G^up_ij de_j G^dn_ji] over all 6 orbitals j, de = W_up W_dn (P_up - P_dn), P = (z - C)/W^2.
   !> Kernel: de_k = a_k + b_k z with a = (C_dn W_up^2 - C_up W_dn^2)/(W_up W_dn), b = (W_dn^2 - W_up^2)/(W_up W_dn).
   subroutine check_j0_local()
      integer, parameter :: n = 6
      complex(rp) :: s(n, n), gs(n, n, 2), mm(n, n), de(n), z, tk(1, 2), gloc(nb, nb, 1), ref
      real(rp) :: cu(n), cd(n), wu(n), wd(n), a(2), b(2), r(2)
      integer :: k, l, iz, is, site, io, ipiv(n), info

      do l = 1, n
         do k = 1, n
            s(k, l) = 0.3_rp*cmplx(sin(1.3_rp*k + 0.7_rp*l), cos(0.9_rp*k - 1.1_rp*l), rp)
         end do
      end do
      s = 0.5_rp*(s + conjg(transpose(s)))
      cu = [0.10_rp, -0.20_rp, 0.15_rp, 0.30_rp, -0.05_rp, 0.25_rp]
      cd = [0.35_rp, 0.10_rp, -0.10_rp, 0.55_rp, 0.20_rp, 0.05_rp]
      wu = [0.90_rp, 1.10_rp, 1.00_rp, 0.80_rp, 1.20_rp, 0.95_rp]
      wd = [1.05_rp, 0.85_rp, 1.10_rp, 0.90_rp, 1.00_rp, 1.15_rp]
      do iz = 1, 2
         z = merge((0.3_rp, 0.2_rp), (-0.5_rp, 0.4_rp), iz == 1)
         do is = 1, 2
            do l = 1, n
               do k = 1, n
                  mm(k, l) = -merge(wu(k)*wu(l), wd(k)*wd(l), is == 1)*s(k, l)
               end do
               mm(l, l) = mm(l, l) + z - merge(cu(l), cd(l), is == 1)
            end do
            gs(:, :, is) = (0.0_rp, 0.0_rp)
            do l = 1, n
               gs(l, l, is) = (1.0_rp, 0.0_rp)
            end do
            call zgesv(n, n, mm, n, ipiv, gs(:, :, is), n, info)
            if (info /= 0) error stop 'zgesv failed'
         end do
         de = wu*wd*((z - cu)/wu**2 - (z - cd)/wd**2)
         do site = 1, 3
            k = 2*(site - 1)
            gloc = (0.0_rp, 0.0_rp)
            gloc(1:2, 1:2, 1) = gs(k + 1:k + 2, k + 1:k + 2, 1)
            gloc(spin_off + 1:spin_off + 2, spin_off + 1:spin_off + 2, 1) = gs(k + 1:k + 2, k + 1:k + 2, 2)
            a = (cd(k + 1:k + 2)*wu(k + 1:k + 2)**2 - cu(k + 1:k + 2)*wd(k + 1:k + 2)**2)/(wu(k + 1:k + 2)*wd(k + 1:k + 2))
            b = (wd(k + 1:k + 2)**2 - wu(k + 1:k + 2)**2)/(wu(k + 1:k + 2)*wd(k + 1:k + 2))
            r = wd(k + 1:k + 2)/wu(k + 1:k + 2)
            tk = j0_local_traces(gloc, a, b, r)
            ref = (0.0_rp, 0.0_rp)
            do io = 1, 2
               do l = 1, n
                  ref = ref + de(k + io)*gs(k + io, l, 1)*de(l)*gs(l, k + io, 2)
               end do
            end do
            call require(abs(tk(1, 1) + z*tk(1, 2) - ref) < 1.0e-12_rp*abs(ref), 'j0_local_traces = direct cluster sum')
         end do
      end do
   end subroutine check_j0_local
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

   !> lorentz_smooth: G(E) = 1/(E + i g) + 0.5/(E - 0.3 + i g) on a uniform mesh; the Poisson average with eta must give the
   !> closed form G(E + i (g + eta)) away from the mesh ends (window error about 3e-4). Catches a wrong sign of eta, a missing
   !> pi or dE factor, a window in the wrong units. eta = 0 returns the input; dynamic_propagator changes only the d spin blocks.
   subroutine check_lorentz_smooth()
      integer, parameter :: ne = 8001
      real(rp), parameter :: g = 0.02_rp, eta = 0.01_rp, de = 0.0005_rp
      complex(rp) :: gm(1, 1, ne), gs(1, 1, ne), ex(ne)
      complex(rp), allocatable :: gfull(:, :, :), gp(:, :, :)
      real(rp) :: e(ne), err
      integer :: n

      e = [(-2.0_rp + (n - 1)*de, n=1, ne)]
      gm(1, 1, :) = 1.0_rp/cmplx(e, g, rp) + 0.5_rp/cmplx(e - 0.3_rp, g, rp)
      ex = 1.0_rp/cmplx(e, g + eta, rp) + 0.5_rp/cmplx(e - 0.3_rp, g + eta, rp)
      gs = lorentz_smooth(gm, de, eta)
      err = maxval(abs(gs(1, 1, :) - ex), mask=abs(e) < 1.0_rp)/maxval(abs(ex), mask=abs(e) < 1.0_rp)
      call require(err < 1.0e-3_rp, 'lorentz_smooth = G(E + i (gamma + eta)) within 1e-3')
      gs = lorentz_smooth(gm, de, 0.0_rp)
      call require(all(gs == gm), 'lorentz_smooth with eta = 0 returns the input')
      allocate (gfull(nb, nb, 50))
      gfull = (1.0_rp, -2.0_rp)
      do n = 1, 50
         gfull(:, :, n) = gfull(:, :, n)*n
      end do
      gp = dynamic_propagator(gfull, de, 4.0_rp*de)
      call require(all(gp(1:4, :, :) == gfull(1:4, :, :)) .and. all(gp(10:20, :, :) == gfull(10:20, :, :)) .and. &
                   all(gp(26:, :, :) == gfull(26:, :, :)) .and. all(gp(5:9, 1:4, :) == gfull(5:9, 1:4, :)), &
                   'dynamic_propagator leaves everything but the d spin blocks unchanged')
      call require(any(gp(5:9, 5:9, :) /= gfull(5:9, 5:9, :)) .and. any(gp(21:25, 21:25, :) /= gfull(21:25, 21:25, :)), &
                   'dynamic_propagator smooths both d spin blocks')
   end subroutine check_lorentz_smooth

   subroutine require(ok, label)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label
      if (.not. ok) then
         failed = .true.
         write (*, '(2a)') 'FAIL: ', label
      end if
   end subroutine require
end program test_spin_response_rs
