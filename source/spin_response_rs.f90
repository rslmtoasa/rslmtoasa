!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: spin_response_rs
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> Real-space transverse spin response. Static bare susceptibility of the d channel per pair (i, j),
!> with the coefficient vertex sum_mu c+_{i mu up} c_{i mu dn}:
!>
!>   chi0_ij = -(X - conj(Y))/(2 pi i),  X = int f Tr_d[G^upup_ji G^dndn_ij] dE,  Y = int f Tr_d[G^dndn_ji G^upup_ij] dE
!>
!> with retarded G. Only the spin-diagonal blocks of the spinor G enter, also with SOC; with time reversal X = Y and
!> chi0_ij = -Im X/pi. Sign and units (1/Ry) are those of the k-space spin_response with method = 'mills',
!> channels = 'd'. The kernels take plain arrays; the exchange pair loop (calculate_pairwise) calls them.
!>
!> Goldstone diagnostics per atom type a, from the first self pair (i, i) with type a (the seed of a): its d moment
!> m_a = -(1/pi) Im int f Tr_d[G^upup_ii - G^dndn_ii] dE, and chi0(q = 0)_ab = sum of Re chi0_ij over the listed pairs
!> of that seed with type(j) = b (truncated to the pair list). As in spin_response: U from the Goldstone condition
!> sum_b chi0_ab U_b m_b = -m_a (by construction), U_Mills = (C_dn - C_up)_d/m and the Mills residual.
!>
!> Juelich vertex: the amplitudes of the d orbitals are multiplied by F(e) = alpha + beta e at the energy e of each state,
!> F = (1 + (e - E_nu)(E_F - E_nu) p)/sqrt(1 + (E_F - E_nu)^2 p), p = 1/ppar^2 (as d_amplitudes_juelich in spin_response).
!> In Green-function form G is replaced by G~_ab = <a|F_a(H) G F_b(H)|b>, which with HG = GH = zG - 1 is
!>   G~_ab = F_a(z) F_b(z) G_ab - beta_a beta_b H_ab - (alpha_a beta_b + beta_a alpha_b + beta_a beta_b z) delta_ab,
!> and in Chebyshev moments (H = a Ht + b, F = alpha' + beta' Ht) mu~_m = <a|F(Ht) T_m(Ht) F(Ht)|b>.
!>
!> Adiabatic spectrum: chi0_ab(q) = sum over the listed pairs (i, j) of seed(a) with type(j) = b of
!> chi0_ij exp(+2 pi i q.r_ij), r_ij = r_j - r_i in direct coordinates; with every pair included this is the chi0_ab(q)
!> of the k-space spin_response (its Bloch gauge carries the phase exp(2 pi i k.(r_j - r_i))). Then
!> Omega(q) = -N (chi0(q)^-1 + U) with N = diag(m) and the Goldstone U (adiabatic_frequencies of spin_response).
!>
!> Cluster of magnetic sites: S = the atoms i with a self pair (i, i), in pair-list order, and chi0_S(a, b) = chi0 of the
!> listed pair (S_a, S_b). The static chi0^{+-} matrix is Hermitian: from the definition above, with the cyclic trace,
!> X_ji = Tr[G^upup_ij G^dndn_ji] = Tr[G^dndn_ji G^upup_ij] = Y_ij and Y_ji = X_ij, so chi0_ji = -(Y_ij - conj(X_ij))/(2 pi i)
!> = conj(chi0_ij). A pair listed in one order only is completed by the conjugate; a pair listed in neither order is 0
!> (truncation to the pair list). m_a is the d moment of the self pair of S_a, U the Goldstone U of Re chi0_S and m, and
!> Omega = -N (chi0_S^-1 + U), |S| frequencies, as for the lattice sum with one 'type' per site.
!> Local J0 sum rule (do_exchange): with H^s = C^s + W^s S W^s, S spin independent, g^s = (P^s - S)^-1 = W^s G^s W^s and
!> g_up - g_dn = g_up (P_dn - P_up) g_dn, the vertex de = W_up W_dn (P_up - P_dn) gives
!>   sum_j G^up_ij de_j G^dn_ji = (W_up)^-1 W_dn G^dn_ii - G^up_ii W_up (W_dn)^-1,
!> so sum_j Tr[de_i G^up_ij de_j G^dn_ji] = sum_k de_k [r_k G^dn_kk - G^up_kk / r_k], r = W_dn/W_up, k over the orbitals of
!> the vertex. The down-up term of the collinear J integrand (Gn Gn - Gz Gz = (up dn + dn up)/2 per trace) gives the same
!> expression, so J0_i = sum_j J_ij is the integral of Im of it. Holds for nsp = 1, hoh = .false. and the full cluster sum.
!------------------------------------------------------------------------------
module spin_response_rs_mod
   use precision_mod, only: rp
   use math_mod, only: pi, i_unit
   use basis_mod, only: spin_off
   use spin_response_mod, only: u_juelich, u_mills, mills_residual, adiabatic_frequencies
   use logger_mod, only: g_logger
   use string_mod, only: int2str
   implicit none
   private

   public :: check_chi0_input
   public :: chi0_traces
   public :: chi0_from_traces
   public :: write_chi0_pairs
   public :: d_moment_trace
   public :: type_sums
   public :: write_chi0_types
   public :: juelich_vectors
   public :: juelich_green
   public :: juelich_moments
   public :: type_sums_q
   public :: write_chi0_adiabatic
   public :: cluster_sites
   public :: cluster_frequencies
   public :: write_chi0_cluster
   public :: j0_local_traces
   public :: write_j0_sum_rule

   !> Angular momentum of the response channel.
   integer, parameter :: l_d = 2
   !> |m| (mu_B) below which a type is treated as non-magnetic and the U columns are n/a.
   real(rp), parameter :: m_tol = 1.0e-6_rp

contains

   !> @brief Stop unless the run is collinear (nsp = 1, or nsp = 2 with every moment along +-z) and has a d channel: the
   !>        spin blocks of G are then the up and down states of the transverse response.
   !> @param[in] nsp   Spin-block count of the Hamiltonian.
   !> @param[in] lmax  Highest angular momentum of the basis.
   !> @param[in] mom   Moment vector of each type, shape (3, ntype).
   subroutine check_chi0_input(nsp, lmax, mom)
      integer, intent(in) :: nsp, lmax
      real(rp), intent(in) :: mom(:, :)

      integer :: it

      if (nsp /= 1 .and. nsp /= 2) &
         call g_logger%fatal('do_chi0 needs nsp = 1 or 2 (collinear, with or without SOC)', __FILE__, __LINE__)
      if (lmax < 2) call g_logger%fatal('do_chi0 needs lmax >= 2 (d channel)', __FILE__, __LINE__)
      if (nsp == 1) return
      do it = 1, size(mom, 2)
         if (norm2(mom(1:2, it)) > 1.0e-8_rp*norm2(mom(:, it))) &
            call g_logger%fatal('do_chi0 with nsp = 2 needs every moment along +-z (type '//int2str(it)//')', __FILE__, __LINE__)
      end do
   end subroutine check_chi0_input

   !> @brief Tr_d[G^upup_ji G^dndn_ij] (column 1) and Tr_d[G^dndn_ji G^upup_ij] (column 2) for each index n of the third
   !>        dimension (energy point or moment).
   !> @param[in] gij, gji  Spinor blocks G_ij, G_ji, shape (nb, nb, n), spin-down block from spin_off + 1.
   pure function chi0_traces(gij, gji) result(t)
      complex(rp), intent(in) :: gij(:, :, :), gji(:, :, :)
      complex(rp) :: t(size(gij, 3), 2)

      integer :: n, u0, u1, d0, d1

      u0 = l_d**2 + 1
      u1 = (l_d + 1)**2
      d0 = u0 + spin_off
      d1 = u1 + spin_off
      do n = 1, size(gij, 3)
         t(n, 1) = sum(transpose(gji(u0:u1, u0:u1, n))*gij(d0:d1, d0:d1, n))
         t(n, 2) = sum(transpose(gji(d0:d1, d0:d1, n))*gij(u0:u1, u0:u1, n))
      end do
   end function chi0_traces

   !> @brief chi0_ij = -(x - conj(y))/(2 pi i) from the energy integrals x, y of the two columns of chi0_traces.
   elemental function chi0_from_traces(x, y) result(chi0)
      complex(rp), intent(in) :: x, y
      complex(rp) :: chi0

      chi0 = -(x - conjg(y))/(2.0_rp*pi*i_unit)
   end function chi0_from_traces

   !> @brief Write one row per pair: atoms, types, r_j - r_i and |r_j - r_i| (cluster coordinates), Re and Im chi0_ij, and
   !>        Re and Im of the Juelich chi0_ij (n/a when chi0j is absent).
   !> @param[in] ijpair  Pair list, shape (npair, 2).
   !> @param[in] iz      Type of each cluster atom.
   !> @param[in] cr      Cluster coordinates, shape (3, natom).
   !> @param[in] chi0    chi0_ij (1/Ry), shape (npair); chi0j optional, the same with the Juelich vertex.
   subroutine write_chi0_pairs(fname, ijpair, iz, cr, chi0, chi0j)
      character(len=*), intent(in) :: fname
      integer, intent(in) :: ijpair(:, :), iz(:)
      real(rp), intent(in) :: cr(:, :)
      complex(rp), intent(in) :: chi0(:)
      complex(rp), intent(in), optional :: chi0j(:)

      integer :: funit, p, i, j
      character(len=48) :: jcols

      open (newunit=funit, file=fname, action='write', status='replace')
      write (funit, '(a)') '# i j type_i type_j r_ij(3, cluster coordinates) |r_ij| re_chi0(1/Ry) im_chi0(1/Ry) '// &
         're_chi0_juelich(1/Ry) im_chi0_juelich(1/Ry)'
      do p = 1, size(chi0)
         i = ijpair(p, 1)
         j = ijpair(p, 2)
         if (present(chi0j)) then
            write (jcols, '(2es24.14)') real(chi0j(p), rp), aimag(chi0j(p))
         else
            write (jcols, '(2a24)') 'n/a', 'n/a'
         end if
         write (funit, '(4i8,4f12.6,2es24.14,a)') i, j, iz(i), iz(j), cr(:, j) - cr(:, i), norm2(cr(:, j) - cr(:, i)), &
            real(chi0(p), rp), aimag(chi0(p)), jcols
      end do
      close (funit)
   end subroutine write_chi0_pairs

   !> @brief Tr_d[G^upup - G^dndn] of a self pair for each index n of the third dimension (energy point or moment).
   pure function d_moment_trace(g) result(t)
      complex(rp), intent(in) :: g(:, :, :)
      complex(rp) :: t(size(g, 3))

      integer :: n, mu

      do n = 1, size(g, 3)
         t(n) = (0.0_rp, 0.0_rp)
         do mu = l_d**2 + 1, (l_d + 1)**2
            t(n) = t(n) + g(mu, mu, n) - g(mu + spin_off, mu + spin_off, n)
         end do
      end do
   end function d_moment_trace

   !> @brief Juelich coefficients alpha(nb), beta(nb) of one site: F = alpha + beta e on the d orbitals of each spin block,
   !>        alpha = 1, beta = 0 on all other orbitals.
   !> @param[in] enu, ppar  E_nu (Ry) and ppar of the d channel per physical spin, shape (2); ef the Fermi level (Ry).
   pure subroutine juelich_vectors(enu, ppar, ef, alpha, beta)
      real(rp), intent(in) :: enu(2), ppar(2), ef
      real(rp), intent(out) :: alpha(:), beta(:)

      integer :: s, mu
      real(rp) :: x, p, d

      alpha = 1.0_rp
      beta = 0.0_rp
      do s = 1, 2
         x = ef - enu(s)
         p = 1.0_rp/ppar(s)**2
         d = sqrt(1.0_rp + x*x*p)
         do mu = (s - 1)*spin_off + l_d**2 + 1, (s - 1)*spin_off + (l_d + 1)**2
            alpha(mu) = (1.0_rp - enu(s)*x*p)/d
            beta(mu) = x*p/d
         end do
      end do
   end subroutine juelich_vectors

   !> @brief G~ = <a|F_a(H) G F_b(H)|b> on real energies z from G, the Hamiltonian block h = H_ab and the row (a) and
   !>        column (b) Juelich coefficients; the delta_ab term only for a self pair.
   !> @param[in] g  G_ab, shape (nb, nb, n); z the energies (Ry), shape (n).
   pure subroutine juelich_green(g, h, z, alpha_a, beta_a, alpha_b, beta_b, self_pair, gt)
      complex(rp), intent(in) :: g(:, :, :), h(:, :)
      real(rp), intent(in) :: z(:), alpha_a(:), beta_a(:), alpha_b(:), beta_b(:)
      logical, intent(in) :: self_pair
      complex(rp), intent(out) :: gt(:, :, :)

      integer :: n, k, l

      do n = 1, size(g, 3)
         do l = 1, size(g, 2)
            do k = 1, size(g, 1)
               gt(k, l, n) = (alpha_a(k) + beta_a(k)*z(n))*(alpha_b(l) + beta_b(l)*z(n))*g(k, l, n) - beta_a(k)*beta_b(l)*h(k, l)
            end do
            if (self_pair) gt(l, l, n) = gt(l, l, n) - (alpha_a(l)*beta_b(l) + beta_a(l)*alpha_b(l) + beta_a(l)*beta_b(l)*z(n))
         end do
      end do
   end subroutine juelich_green

   !> @brief Moments mu~_m = <a|F_a T_m(Ht) F_b|b>, m = 0..n-3, from the raw moments mu_m = <a|T_m(Ht)|b>, m = 0..n-1, with
   !>        H = a_win Ht + b_win: F = alpha' + beta' Ht, alpha' = alpha + beta b_win, beta' = beta a_win, and
   !>        Ht T_m = (T_m+1 + T_|m-1|)/2, Ht^2 T_m = (T_m+2 + 2 T_m + T_|m-2|)/4.
   pure subroutine juelich_moments(mu, a_win, b_win, alpha_a, beta_a, alpha_b, beta_b, mt)
      complex(rp), intent(in) :: mu(:, :, :)
      real(rp), intent(in) :: a_win, b_win, alpha_a(:), beta_a(:), alpha_b(:), beta_b(:)
      complex(rp), intent(out) :: mt(:, :, :)

      integer :: m, k, l
      real(rp) :: aa(size(alpha_a)), ba(size(alpha_a)), ab(size(alpha_b)), bb(size(alpha_b))

      aa = alpha_a + beta_a*b_win
      ba = beta_a*a_win
      ab = alpha_b + beta_b*b_win
      bb = beta_b*a_win
      do m = 0, size(mu, 3) - 3
         do l = 1, size(mu, 2)
            do k = 1, size(mu, 1)
               mt(k, l, m + 1) = aa(k)*ab(l)*mu(k, l, m + 1) &
                                 + 0.5_rp*(aa(k)*bb(l) + ba(k)*ab(l))*(mu(k, l, m + 2) + mu(k, l, abs(m - 1) + 1)) &
                                 + 0.25_rp*ba(k)*bb(l)*(mu(k, l, m + 3) + 2.0_rp*mu(k, l, m + 1) + mu(k, l, abs(m - 2) + 1))
            end do
         end do
      end do
   end subroutine juelich_moments

   !> @brief Type sums of the pair list. seed(a) is the first atom i of a self pair (i, i) of type a (0: none), m(a) its
   !>        d moment md, npair(a) the number of pairs with first atom seed(a) and chi0s(a, b) the sum of Re chi0 over
   !>        those pairs with iz(j) = b. Rows of types without a seed stay 0.
   !> @param[in] ijpair  Pair list, shape (npair, 2); iz the type of each cluster atom.
   !> @param[in] chi0    chi0_ij (1/Ry) and md the d moment of each pair (used for self pairs only), shape (npair).
   pure subroutine type_sums(ijpair, iz, ntype, chi0, md, seed, npair, m, chi0s)
      integer, intent(in) :: ijpair(:, :), iz(:), ntype
      complex(rp), intent(in) :: chi0(:)
      real(rp), intent(in) :: md(:)
      integer, intent(out) :: seed(ntype), npair(ntype)
      real(rp), intent(out) :: m(ntype), chi0s(ntype, ntype)

      integer :: p, a

      seed = 0
      npair = 0
      m = 0.0_rp
      chi0s = 0.0_rp
      do p = 1, size(chi0)
         a = iz(ijpair(p, 1))
         if (ijpair(p, 1) == ijpair(p, 2) .and. seed(a) == 0) then
            seed(a) = ijpair(p, 1)
            m(a) = md(p)
         end if
      end do
      do p = 1, size(chi0)
         a = iz(ijpair(p, 1))
         if (seed(a) == 0 .or. ijpair(p, 1) /= seed(a)) cycle
         npair(a) = npair(a) + 1
         chi0s(a, iz(ijpair(p, 2))) = chi0s(a, iz(ijpair(p, 2))) + real(chi0(p), rp)
      end do
   end subroutine type_sums

   !> @brief Write one row per type with a seed: type, seed atom, pair count, m_d and m_d Juelich, U_Goldstone (coefficient and
   !>        Juelich vertex), U_Mills, Mills residual, chi0(q = 0)(a, b) for every b, and the same with the Juelich vertex.
   !>        The U columns of a vertex are n/a unless every type has a seed and |m| >= m_tol; Juelich columns are n/a when
   !>        mj is absent.
   !> @param[in] c_up, c_dn  d band centres (Ry) of each type, shape (ntype).
   !> @param[in] mj, chi0sj  Optional: m and chi0s with the Juelich vertex.
   subroutine write_chi0_types(fname, seed, npair, m, chi0s, c_up, c_dn, mj, chi0sj)
      character(len=*), intent(in) :: fname
      integer, intent(in) :: seed(:), npair(:)
      real(rp), intent(in) :: m(:), chi0s(:, :), c_up(:), c_dn(:)
      real(rp), intent(in), optional :: mj(:), chi0sj(:, :)

      integer :: funit, a, b
      real(rp) :: ug(size(m)), uj(size(m)), um(size(m)), delta(size(m))
      character(len=24), parameter :: na = repeat(' ', 21)//'n/a'
      character(len=24) :: c_mj, c_ug, c_uj, c_um, c_delta
      character(len=24), allocatable :: c_sj(:)
      logical :: have_u, have_uj

      have_u = all(seed > 0) .and. all(abs(m) >= m_tol)
      have_uj = present(mj)
      if (have_uj) have_uj = all(seed > 0) .and. all(abs(mj) >= m_tol)
      if (have_u) then
         call u_juelich(chi0s, m, ug)
         um = u_mills(c_up, c_dn, m)
         delta = mills_residual(chi0s, um, m)
      end if
      if (have_uj) call u_juelich(chi0sj, mj, uj)
      allocate (c_sj(size(m)))
      open (newunit=funit, file=fname, action='write', status='replace')
      write (funit, '(a)') '# type seed_atom n_pairs m_d(mu_B) m_d_juelich(mu_B) u_goldstone(Ry, by construction) '// &
         'u_goldstone_juelich(Ry, by construction) u_mills(Ry) mills_residual chi0_q0(type, b)(1/Ry) for b = 1..ntype '// &
         'chi0_q0_juelich(type, b)(1/Ry) for b = 1..ntype'
      do a = 1, size(m)
         if (seed(a) == 0) cycle
         c_mj = na; c_ug = na; c_uj = na; c_um = na; c_delta = na; c_sj = na
         if (present(mj)) then
            write (c_mj, '(es24.14)') mj(a)
            do b = 1, size(m)
               write (c_sj(b), '(es24.14)') chi0sj(a, b)
            end do
         end if
         if (have_u) then
            write (c_ug, '(es24.14)') ug(a)
            write (c_um, '(es24.14)') um(a)
            write (c_delta, '(es24.14)') delta(a)
         end if
         if (have_uj) write (c_uj, '(es24.14)') uj(a)
         write (funit, '(3i8,es24.14,5a24,*(es24.14))', advance='no') a, seed(a), npair(a), m(a), c_mj, c_ug, c_uj, c_um, &
            c_delta, chi0s(a, :)
         write (funit, '(*(a24))') c_sj
      end do
      close (funit)
   end subroutine write_chi0_types

   !> @brief chi0_ab(q) = sum over the pairs p with ijpair(p, 1) = seed(a) and iz(ijpair(p, 2)) = b of chi0(p) exp(2 pi i q.r(p)).
   !> @param[in] rdir  r_j - r_i of each pair in direct coordinates, shape (3, npair).
   !> @param[in] q     q points in direct coordinates, shape (3, nq).
   !> @return    cq    shape (ntype, ntype, nq); rows of types without a seed stay 0.
   pure function type_sums_q(ijpair, iz, seed, chi0, rdir, q) result(cq)
      integer, intent(in) :: ijpair(:, :), iz(:), seed(:)
      complex(rp), intent(in) :: chi0(:)
      real(rp), intent(in) :: rdir(:, :), q(:, :)
      complex(rp) :: cq(size(seed), size(seed), size(q, 2))

      integer :: p, a, b, iq

      cq = (0.0_rp, 0.0_rp)
      do p = 1, size(chi0)
         a = iz(ijpair(p, 1))
         if (seed(a) == 0 .or. ijpair(p, 1) /= seed(a)) cycle
         b = iz(ijpair(p, 2))
         do iq = 1, size(q, 2)
            cq(a, b, iq) = cq(a, b, iq) + chi0(p)*exp(2.0_rp*pi*i_unit*dot_product(q(:, iq), rdir(:, p)))
         end do
      end do
   end function type_sums_q

   !> @brief Write one row per q: q (direct, Cartesian in 2 pi/a), |q| (1/Angstrom), then for each vertex (coefficient, then
   !>        Juelich) the adiabatic frequencies (Ry, ascending, one per type), the largest |Im| of their eigenvalues (Ry) and
   !>        the force-theorem frequencies N U (I + chi0 U) (Ry, ascending). U is the Goldstone U of chi0s and m (by
   !>        construction Omega(0) m = 0). A vertex whose U is not defined (a type without a seed, |m| < m_tol, or the
   !>        Juelich vertex absent) gives n/a.
   !> @param[in] cq, chi0s, m        chi0(q) (ntype, ntype, nq), chi0(q = 0) (ntype, ntype) and m (ntype), coefficient vertex.
   !> @param[in] cqj, chi0sj, mj     Optional: the same with the Juelich vertex.
   subroutine write_chi0_adiabatic(fname, q, qcart, qabs, seed, cq, chi0s, m, cqj, chi0sj, mj)
      character(len=*), intent(in) :: fname
      real(rp), intent(in) :: q(:, :), qcart(:, :), qabs(:)
      integer, intent(in) :: seed(:)
      complex(rp), intent(in) :: cq(:, :, :)
      real(rp), intent(in) :: chi0s(:, :), m(:)
      complex(rp), intent(in), optional :: cqj(:, :, :)
      real(rp), intent(in), optional :: chi0sj(:, :), mj(:)

      integer :: funit, iq, nt
      character(len=:), allocatable :: line
      character(len=24), parameter :: na = repeat(' ', 21)//'n/a'
      logical :: have_u, have_uj

      nt = size(m)
      have_u = all(seed > 0) .and. all(abs(m) >= m_tol)
      have_uj = present(mj)
      if (have_uj) have_uj = all(seed > 0) .and. all(abs(mj) >= m_tol)
      open (newunit=funit, file=fname, action='write', status='replace')
      write (funit, '(a)') '# q_direct(3) q_cartesian(3, units 2pi/a) |q|(1/Angstrom) omega_adiabatic(Ry, ascending, one per type) '// &
         'max_abs_imag(Ry) omega_force(Ry, ascending) ; the same three groups with the Juelich vertex'
      do iq = 1, size(q, 2)
         line = ''
         if (have_u) then
            line = line//frequencies(cq(:, :, iq), chi0s, m)
         else
            line = line//repeat(na, 2*nt + 1)
         end if
         if (have_uj) then
            line = line//frequencies(cqj(:, :, iq), chi0sj, mj)
         else
            line = line//repeat(na, 2*nt + 1)
         end if
         write (funit, '(7es16.8,a)') q(:, iq), qcart(:, iq), qabs(iq), line
      end do
      close (funit)

   contains

      !> Omega, max |Im| and the force-theorem frequencies at one q as text, with the Goldstone U of chi0_static and mm.
      function frequencies(chi0q, chi0_static, mm) result(text)
         complex(rp), intent(in) :: chi0q(:, :)
         real(rp), intent(in) :: chi0_static(:, :), mm(:)
         character(len=:), allocatable :: text
         real(rp) :: u(size(mm)), kmat(size(mm), size(mm)), omega(size(mm)), omega_force(size(mm)), imag_max
         character(len=24*(2*size(mm) + 1)) :: buf
         integer :: k

         call u_juelich(chi0_static, mm, u)
         kmat = 0.0_rp
         do k = 1, size(mm)
            kmat(k, k) = u(k)
         end do
         call adiabatic_frequencies(chi0q, kmat, mm, omega, imag_max, omega_force)
         write (buf, '(*(es24.14))') omega, imag_max, omega_force
         text = buf
      end function frequencies
   end subroutine write_chi0_adiabatic

   !> @brief Sites of the cluster and their chi0 matrix. site(a) is the a-th atom with a self pair (first occurrence, pair-list
   !>        order), m(a) the md of that self pair and chi0c(a, b) the chi0 of the listed pair (site(a), site(b)); if only
   !>        (site(b), site(a)) is listed its conjugate, if neither is listed 0. Pairs with an atom outside S are ignored.
   !> @param[in] ijpair  Pair list, shape (npair, 2); chi0 (1/Ry) and md the d moment of each pair (self pairs only), shape (npair).
   pure subroutine cluster_sites(ijpair, chi0, md, site, chi0c, m)
      integer, intent(in) :: ijpair(:, :)
      complex(rp), intent(in) :: chi0(:)
      real(rp), intent(in) :: md(:)
      integer, allocatable, intent(out) :: site(:)
      complex(rp), allocatable, intent(out) :: chi0c(:, :)
      real(rp), allocatable, intent(out) :: m(:)

      integer :: p, a, b, nsite
      logical, allocatable :: filled(:, :)

      nsite = 0
      allocate (site(size(chi0)))
      do p = 1, size(chi0)
         if (ijpair(p, 1) /= ijpair(p, 2) .or. any(site(1:nsite) == ijpair(p, 1))) cycle
         nsite = nsite + 1
         site(nsite) = ijpair(p, 1)
      end do
      site = site(1:nsite)
      allocate (chi0c(nsite, nsite), m(nsite), filled(nsite, nsite))
      chi0c = (0.0_rp, 0.0_rp)
      m = 0.0_rp
      filled = .false.
      do p = 1, size(chi0)
         a = findloc(site, ijpair(p, 1), 1)
         b = findloc(site, ijpair(p, 2), 1)
         if (a == 0 .or. b == 0 .or. filled(a, b)) cycle
         chi0c(a, b) = chi0(p)
         filled(a, b) = .true.
         if (a == b) m(a) = md(p)
      end do
      do a = 1, nsite
         do b = 1, nsite
            if (.not. filled(a, b) .and. filled(b, a)) chi0c(a, b) = conjg(chi0c(b, a))
         end do
      end do
   end subroutine cluster_sites

   !> @brief Goldstone U of Re chi0c and m, then Omega = -N (chi0c^-1 + U) (omega), its largest |Im| and the force-theorem
   !>        frequencies (omega_force), both ascending.
   subroutine cluster_frequencies(chi0c, m, u, omega, imag_max, omega_force)
      complex(rp), intent(in) :: chi0c(:, :)
      real(rp), intent(in) :: m(:)
      real(rp), intent(out) :: u(:), omega(:), imag_max, omega_force(:)

      real(rp) :: kmat(size(m), size(m))
      integer :: k

      call u_juelich(real(chi0c, rp), m, u)
      kmat = 0.0_rp
      do k = 1, size(m)
         kmat(k, k) = u(k)
      end do
      call adiabatic_frequencies(chi0c, kmat, m, omega, imag_max, omega_force)
   end subroutine cluster_frequencies

   !> @brief Write chi0_cluster.dat in three blocks, each after its own header line: one row per site (index, atom, type, m_d
   !>        and Goldstone U of the coefficient and of the Juelich vertex), one row per frequency index (adiabatic and force-theorem
   !>        frequencies of both vertices) and one row with the largest |Im| of each vertex. A vertex is n/a when there are no sites,
   !>        a site has |m| < m_tol, or (Juelich) mj is absent; m is printed whenever it exists.
   !> @param[in] site, iz  Atoms of the cluster (cluster_sites) and the type of each cluster atom.
   !> @param[in] m, chi0c  d moments and chi0 matrix of the coefficient vertex; mj, chi0cj optional, the same with the Juelich vertex.
   subroutine write_chi0_cluster(fname, site, iz, m, chi0c, mj, chi0cj)
      character(len=*), intent(in) :: fname
      integer, intent(in) :: site(:), iz(:)
      real(rp), intent(in) :: m(:)
      complex(rp), intent(in) :: chi0c(:, :)
      real(rp), intent(in), optional :: mj(:)
      complex(rp), intent(in), optional :: chi0cj(:, :)

      integer :: funit, a, n
      real(rp) :: u(size(m)), uj(size(m)), om(size(m)), omf(size(m)), omj(size(m)), omfj(size(m)), mja(size(m)), im, imj
      logical :: have, havej

      n = size(m)
      u = 0.0_rp; uj = 0.0_rp; om = 0.0_rp; omf = 0.0_rp; omj = 0.0_rp; omfj = 0.0_rp; mja = 0.0_rp; im = 0.0_rp; imj = 0.0_rp
      if (present(mj)) mja = mj
      have = n > 0
      if (have) have = all(abs(m) >= m_tol)
      havej = present(mj)
      if (havej) havej = n > 0
      if (havej) havej = all(abs(mj) >= m_tol)
      if (have) call cluster_frequencies(chi0c, m, u, om, im, omf)
      if (havej) call cluster_frequencies(chi0cj, mj, uj, omj, imj, omfj)
      open (newunit=funit, file=fname, action='write', status='replace')
      write (funit, '(a)') '# site atom type m_d(mu_B) u_goldstone(Ry) m_d_juelich(mu_B) u_goldstone_juelich(Ry)'
      do a = 1, n
         write (funit, '(3i8,a)') a, site(a), iz(site(a)), real_cell(m(a), .true.)//real_cell(u(a), have)// &
            real_cell(mja(a), present(mj))//real_cell(uj(a), havej)
      end do
      write (funit, '(a)') '# k omega_adiabatic(Ry) omega_force(Ry) omega_adiabatic_juelich(Ry) omega_force_juelich(Ry)'
      do a = 1, n
         write (funit, '(i8,a)') a, real_cell(om(a), have)//real_cell(omf(a), have)//real_cell(omj(a), havej)// &
            real_cell(omfj(a), havej)
      end do
      write (funit, '(a)') '# max_abs_imag(Ry) max_abs_imag_juelich(Ry)'
      write (funit, '(a)') real_cell(im, have)//real_cell(imj, havej)
      close (funit)
   end subroutine write_chi0_cluster

   !> One 24-character column: x in es24.14 when ok, otherwise n/a right-justified.
   function real_cell(x, ok) result(text)
      real(rp), intent(in) :: x
      logical, intent(in) :: ok
      character(len=24) :: text

      if (ok) then
         write (text, '(es24.14)') x
      else
         text = repeat(' ', 21)//'n/a'
      end if
   end function real_cell
   !> @brief sum_k de_k [r_k G^dn_kk - G^up_kk/r_k] of a self pair with de_k = a_k + b_k e: column 1 the a part, column 2
   !>        the b part, for each index n of the third dimension (energy point or moment). k runs over size(a) orbitals.
   !> @param[in] g     G_ii, shape (nb, nb, n), spin-down block from spin_off + 1.
   !> @param[in] a, b  Vertex de_k(e) = a_k + b_k e as d_matrix(e), shape (norb); r = W_dn/W_up per orbital.
   pure function j0_local_traces(g, a, b, r) result(t)
      complex(rp), intent(in) :: g(:, :, :)
      real(rp), intent(in) :: a(:), b(:), r(:)
      complex(rp) :: t(size(g, 3), 2)

      integer :: n, k
      complex(rp) :: x

      do n = 1, size(g, 3)
         t(n, :) = (0.0_rp, 0.0_rp)
         do k = 1, size(a)
            x = r(k)*g(k + spin_off, k + spin_off, n) - g(k, k, n)/r(k)
            t(n, 1) = t(n, 1) + a(k)*x
            t(n, 2) = t(n, 2) + b(k)*x
         end do
      end do
   end function j0_local_traces

   !> @brief One row per self pair (i, i): atom, type, number of pairs with first atom i, local J0 (mRy), sum of the listed
   !>        J_ij of first atom i (mRy) and their difference.
   !> @param[in] ijpair  Pair list, shape (npair, 2); iz the type of each atom.
   !> @param[in] j0      Local J0 per pair (mRy), read at the self pairs; jij the J_ij (mRy) per pair.
   subroutine write_j0_sum_rule(fname, ijpair, iz, j0, jij)
      character(len=*), intent(in) :: fname
      integer, intent(in) :: ijpair(:, :), iz(:)
      real(rp), intent(in) :: j0(:), jij(:)

      integer :: funit, p, i, np
      real(rp) :: s

      open (newunit=funit, file=fname, action='write', status='replace')
      write (funit, '(a)') '# atom type n_pairs j0_local(mRy) sum_J_ij(mRy) difference(mRy)'
      do p = 1, size(j0)
         i = ijpair(p, 1)
         if (ijpair(p, 2) /= i) cycle
         np = count(ijpair(1:size(j0), 1) == i)
         s = sum(jij, mask=ijpair(1:size(j0), 1) == i)
         write (funit, '(3i8,3es24.14)') i, iz(i), np, j0(p), s, j0(p) - s
      end do
      close (funit)
   end subroutine write_j0_sum_rule

end module spin_response_rs_mod
