submodule (spin_response_mod) spin_response_kernels
   implicit none

contains

   !> @brief Add the Lehmann sum of one k set to chi0 (the caller zeroes chi0 and may chunk k).
   !> @details Pairs (n, m) are listed once per call: band n needs a nonzero spin-1 amplitude at k and band m a nonzero
   !>          spin-2 amplitude at k+q (otherwise the vertex is exactly zero), and f_n - f_m must be nonzero.
   !>          eta = 0 is valid only at omega = 0, where a pair with |e_n - e_m| < 1e-10 Ry takes the Fermi divided
   !>          difference -f(1-f)/kt and stays in the list even when f_n = f_m. Nothing is pruned by a threshold.
   !>          Omega is cut into blocks of omega_block; each block is one matrix product over the pair list, so the
   !>          result does not depend on the number of threads.
   !> @param[in]    kt     Temperature (Ry).
   !> @param[in]    w      k weights, shape (nk).
   !> @param[in]    e_k    Energies at k, shape (nk, nb); f_k occupations, same shape.
   !> @param[in]    a_k    Amplitudes at k, shape (nk, nb, norb, nsc, 2), row kappa of the response space.
   !> @param[in]    e_kq   Energies at k+q, shape (nk, nb); f_kq, a_kq likewise.
   !> @param[in]    omega  Energies, shape (nw) (Ry).
   !> @param[in]    eta    Broadening (Ry).
   !> @param[inout] chi0   Bare response, shape (nsc, nsc, nw).
   module subroutine accumulate_chi0(kt, w, e_k, f_k, a_k, e_kq, f_kq, a_kq, omega, eta, chi0)
      real(rp), intent(in) :: kt, w(:), e_k(:, :), f_k(:, :), e_kq(:, :), f_kq(:, :), omega(:), eta
      complex(rp), intent(in) :: a_k(:, :, :, :, :), a_kq(:, :, :, :, :)
      complex(rp), intent(inout) :: chi0(:, :, :)

      integer :: ik, n, m, i, j, p, np, nsc, ib, iw0, nbl, iw
      real(rp) :: de, x
      logical :: up(size(w), size(e_k, 2)), down(size(w), size(e_kq, 2))
      real(rp), allocatable :: dfp(:), dep(:), fdp(:)
      complex(rp), allocatable :: pm(:, :), cm(:, :), blk(:, :)
      complex(rp) :: t(size(a_k, 4))

      nsc = size(a_k, 4)
      do n = 1, size(e_k, 2)
         up(:, n) = [(any(a_k(ik, n, :, :, 1) /= (0.0_rp, 0.0_rp)), ik=1, size(w))]
      end do
      do m = 1, size(e_kq, 2)
         down(:, m) = [(any(a_kq(ik, m, :, :, 2) /= (0.0_rp, 0.0_rp)), ik=1, size(w))]
      end do
      np = sum(count(up, 2)*count(down, 2))
      allocate (pm(np, nsc*nsc), dfp(np), dep(np), fdp(np))
      np = 0
      do ik = 1, size(w)
         do n = 1, size(e_k, 2)
            if (.not. up(ik, n)) cycle
            do m = 1, size(e_kq, 2)
               if (.not. down(ik, m)) cycle
               de = e_k(ik, n) - e_kq(ik, m)
               if (f_k(ik, n) == f_kq(ik, m) .and. .not. (eta == 0.0_rp .and. abs(de) < degenerate_tol)) cycle
               np = np + 1
               dfp(np) = f_k(ik, n) - f_kq(ik, m)
               dep(np) = de
               fdp(np) = merge(-f_k(ik, n)*(1.0_rp - f_k(ik, n))/kt, 0.0_rp, abs(de) < degenerate_tol)
               do i = 1, nsc
                  t(i) = sum(conjg(a_k(ik, n, :, i, 1))*a_kq(ik, m, :, i, 2))
               end do
               do j = 1, nsc
                  do i = 1, nsc
                     pm(np, (j - 1)*nsc + i) = w(ik)*t(i)*conjg(t(j))
                  end do
               end do
            end do
         end do
      end do
      if (np == 0) return

      !$omp parallel do default(shared) private(ib, iw0, nbl, iw, p, x, cm, blk) schedule(static) if (size(omega) > omega_block)
      do ib = 1, (size(omega) + omega_block - 1)/omega_block
         iw0 = (ib - 1)*omega_block + 1
         nbl = min(omega_block, size(omega) - iw0 + 1)
         allocate (cm(np, nbl), blk(nbl, nsc*nsc))
         do iw = 1, nbl
            do p = 1, np
               x = omega(iw0 + iw - 1) + dep(p)
               if (eta == 0.0_rp .and. omega(iw0 + iw - 1) == 0.0_rp .and. abs(dep(p)) < degenerate_tol) then
                  cm(p, iw) = cmplx(fdp(p), 0.0_rp, rp)
               else if (dfp(p) == 0.0_rp) then
                  cm(p, iw) = (0.0_rp, 0.0_rp)
               else
                  cm(p, iw) = cmplx(dfp(p)*x, -dfp(p)*eta, rp)/(x*x + eta*eta)
               end if
            end do
         end do
         call zgemm('T', 'N', nbl, nsc*nsc, np, (1.0_rp, 0.0_rp), cm, np, pm, size(pm, 1), (0.0_rp, 0.0_rp), blk, nbl)
         do iw = 1, nbl
            chi0(:, :, iw0 + iw - 1) = chi0(:, :, iw0 + iw - 1) + reshape(blk(iw, :), [nsc, nsc])
         end do
         deallocate (cm, blk)
      end do
      !$omp end parallel do
   end subroutine accumulate_chi0

   !> @brief Juelich U from the Goldstone condition sum_j chi0_ij(0,0) M_j U_j = -M_i.
   !> @param[in]  chi0s  Static bare response chi0(0,0), real, shape (nsite, nsite).
   !> @param[in]  m      Site moments, shape (nsite).
   !> @param[out] u      Interaction (Ry), shape (nsite).
   module subroutine u_juelich(chi0s, m, u)
      real(rp), intent(in) :: chi0s(:, :), m(:)
      real(rp), intent(out) :: u(:)

      real(rp) :: amat(size(m), size(m))
      integer :: j, info, ipiv(size(m))

      do j = 1, size(m)
         amat(:, j) = chi0s(:, j)*m(j)
      end do
      u = -m
      call dgesv(size(m), 1, amat, size(m), ipiv, u, size(m), info)
      if (info /= 0) error stop 'u_juelich: singular Goldstone matrix'
   end subroutine u_juelich

   !> @brief Mills U = (C_dn - C_up)/M per site, from the d band centres.
   elemental module function u_mills(c_up, c_dn, m) result(u)
      real(rp), intent(in) :: c_up, c_dn, m
      real(rp) :: u

      u = (c_dn - c_up)/m
   end function u_mills

   !> @brief ALSDA transverse kernel matrix of one site, K_ab = (1/4 pi) int dr P_a P_b K(r)/r^2 over channels a, b, spherical vertex (L = 0).
   !> @details K(r) = dV/m with m = n_up - n_dn (core included), P_a = sqrt(gfac_up gfac_dn) G_up G_dn + S_up S_dn. The sum
   !>          starts at the second mesh point (dV = 0 at r = 0 while m is not), with the Simpson weights of the full mesh.
   !>          Where |m| <= m_floor max|m| (m = 0 on a mesh point), K is the cubic interpolation from the four nearest
   !>          points in mesh index. The Dyson equation takes U = -K.
   !> @param[in]  r, wq              Mesh (bohr) and quadrature weights including dr/di, shape (nr).
   !> @param[in]  m, dv              n_up - n_dn (bohr^-3) and V_xc,up - V_xc,dn (Ry), shape (nr).
   !> @param[in]  g_up, g_dn         Large components times r of each channel, shape (nr, nch); s_* the small components.
   !> @param[in]  gfac_up, gfac_dn   Scalar-relativistic metric factors, shape (nr, nch).
   !> @param[out] kmat               K_ab (Ry), shape (nch, nch).
   module subroutine radial_kernel_matrix(r, wq, m, dv, g_up, g_dn, s_up, s_dn, gfac_up, gfac_dn, kmat)
      real(rp), intent(in) :: r(:), wq(:), m(:), dv(:), g_up(:, :), g_dn(:, :), s_up(:, :), s_dn(:, :), gfac_up(:, :), gfac_dn(:, :)
      real(rp), intent(out) :: kmat(:, :)

      integer :: i, j, nr
      logical :: guarded(size(r))
      real(rp) :: k(size(r)), p(size(r), size(g_up, 2))

      nr = size(r)
      p = sqrt(gfac_up*gfac_dn)*g_up*g_dn + s_up*s_dn
      guarded = abs(m) <= m_floor*maxval(abs(m))
      k = 0.0_rp
      do i = 2, nr
         if (.not. guarded(i)) k(i) = dv(i)/m(i)
      end do
      do i = 2, nr
         if (.not. guarded(i)) cycle
         if (i < 4 .or. i > nr - 2) error stop 'radial_kernel_matrix: m = 0 within two mesh points of an end'
         if (count(guarded(i - 2:i + 2)) > 1) error stop 'radial_kernel_matrix: m = 0 on neighbouring mesh points'
         k(i) = (-k(i - 2) + 4.0_rp*k(i - 1) + 4.0_rp*k(i + 1) - k(i + 2))/6.0_rp
      end do
      do i = 1, size(kmat, 1)
         do j = 1, size(kmat, 2)
            kmat(i, j) = sum(wq(2:)*p(2:, i)*p(2:, j)*k(2:)/r(2:)**2)/(4.0_rp*pi)
         end do
      end do
   end subroutine radial_kernel_matrix

   !> @brief Mills Goldstone residual per site, delta_i = 1 + sum_j chi0s_ij u_j m_j / m_i (1 + chi0(0,0) U for one site).
   !> @param[in] chi0s  Static bare response chi0(0,0), real, shape (nsite, nsite).
   !> @param[in] u      Interaction (Ry), shape (nsite).
   !> @param[in] m      Site moments, shape (nsite).
   pure module function mills_residual(chi0s, u, m) result(delta)
      real(rp), intent(in) :: chi0s(:, :), u(:), m(:)
      real(rp) :: delta(size(m))

      delta = 1.0_rp + matmul(chi0s, u*m)/m
   end function mills_residual

   !> @brief Juelich U (Goldstone solve with the Juelich chi0 and M), Mills U (centre difference over the Mills M),
   !>        and the Mills residual (Mills chi0 and M with the Mills U).
   !> @param[in]  chi0s_mills, chi0s_juelich  Static chi0(0,0) with each amplitude kind, shape (nsite, nsite).
   !> @param[in]  m_mills, m_juelich          Site moments of each kind, shape (nsite).
   !> @param[in]  c_up, c_dn                  d band centres of the majority and minority label (Ry), shape (nsite).
   module subroutine interaction_values(chi0s_mills, chi0s_juelich, m_mills, m_juelich, c_up, c_dn, uj, um, delta_mills)
      real(rp), intent(in) :: chi0s_mills(:, :), chi0s_juelich(:, :), m_mills(:), m_juelich(:), c_up(:), c_dn(:)
      real(rp), intent(out) :: uj(:), um(:), delta_mills(:)

      call u_juelich(chi0s_juelich, m_juelich, uj)
      um = u_mills(c_up, c_dn, m_mills)
      delta_mills = mills_residual(chi0s_mills, um, m_mills)
   end subroutine interaction_values

   !> @brief Trace of a square matrix.
   pure module function sum_diagonal(a) result(t)
      complex(rp), intent(in) :: a(:, :)
      complex(rp) :: t
      integer :: i

      t = (0.0_rp, 0.0_rp)
      do i = 1, size(a, 1)
         t = t + a(i, i)
      end do
   end function sum_diagonal

   !> @brief Eigenvalues of I + chi0 K at one omega (gauge invariant for several sites).
   function eig_one_plus(chi0, kmat) result(w)
      complex(rp), intent(in) :: chi0(:, :)
      real(rp), intent(in) :: kmat(:, :)
      complex(rp) :: w(size(kmat, 1))

      integer :: n, j, info
      complex(rp) :: a(size(kmat, 1), size(kmat, 1)), vl(1, 1), vr(1, 1), work(4*size(kmat, 1))
      real(rp) :: rwork(2*size(kmat, 1))

      n = size(kmat, 1)
      a = matmul(chi0, kmat)
      do j = 1, n
         a(j, j) = a(j, j) + 1.0_rp
      end do
      call zgeev('N', 'N', n, a, n, w, vl, 1, vr, 1, work, size(work), rwork, info)
      if (info /= 0) error stop 'eig_one_plus: zgeev failed'
   end function eig_one_plus

   !> @brief Smallest |eigenvalue| of I + chi0 K at one omega.
   module function min_abs_eig(chi0, kmat) result(m)
      complex(rp), intent(in) :: chi0(:, :)
      real(rp), intent(in) :: kmat(:, :)
      real(rp) :: m

      m = minval(abs(eig_one_plus(chi0, kmat)))
   end function min_abs_eig

   !> @brief Smallest real part of the eigenvalues of I + chi0 K at one omega (signed; the static delta_q at omega = 0, eta = 0).
   module function min_real_eig(chi0, kmat) result(m)
      complex(rp), intent(in) :: chi0(:, :)
      real(rp), intent(in) :: kmat(:, :)
      real(rp) :: m

      m = minval(real(eig_one_plus(chi0, kmat), rp))
   end function min_real_eig

   !> @brief Solve (I + chi0 K) chi = chi0 at each omega by LU, with K a matrix over kappa.
   !> @details Never call at eta = 0 for q = 0: the Goldstone condition makes I + chi0 K singular there.
   !> @param[in]  chi0  Bare response, shape (nsc, nsc, nw).
   !> @param[in]  kmat  Interaction kernel (Ry), shape (nsc, nsc).
   !> @param[out] chi   Enhanced response, shape (nsc, nsc, nw).
   module subroutine solve_dyson_matrix(chi0, kmat, chi)
      complex(rp), intent(in) :: chi0(:, :, :)
      real(rp), intent(in) :: kmat(:, :)
      complex(rp), intent(out) :: chi(:, :, :)

      integer :: ns, iw, j, info, ipiv(size(kmat, 1))
      complex(rp) :: amat(size(kmat, 1), size(kmat, 1))

      ns = size(kmat, 1)
      do iw = 1, size(chi0, 3)
         amat = matmul(chi0(:, :, iw), kmat)
         do j = 1, ns
            amat(j, j) = amat(j, j) + 1.0_rp
         end do
         chi(:, :, iw) = chi0(:, :, iw)
         call zgesv(ns, ns, amat, ns, ipiv, chi(:, :, iw), ns, info)
         if (info /= 0) error stop 'solve_dyson: singular I + chi0 K'
      end do
   end subroutine solve_dyson_matrix

   !> @brief solve_dyson_matrix with the diagonal kernel K = diag(u).
   module subroutine solve_dyson_diagonal(chi0, u, chi)
      complex(rp), intent(in) :: chi0(:, :, :)
      real(rp), intent(in) :: u(:)
      complex(rp), intent(out) :: chi(:, :, :)

      real(rp) :: kmat(size(u), size(u))
      integer :: i

      kmat = 0.0_rp
      do i = 1, size(u)
         kmat(i, i) = u(i)
      end do
      call solve_dyson_matrix(chi0, kmat, chi)
   end subroutine solve_dyson_diagonal

   !> @brief Adiabatic frequencies Omega = -N (chi0^-1 + K) from the static chi0(q, 0), N = diag(m) the signed moments.
   !> @details For chi(w) = (w - Omega)^-1 N the static value is -Omega^-1 N, and (I + chi0 K) chi = chi0 gives
   !>          chi^-1 = chi0^-1 + K. Omega is similar to a Hermitian matrix, so its eigenvalues are real: positive
   !>          where m > 0 and negative where m < 0. At q = 0 the Goldstone condition makes Omega m = 0.
   !> @param[in]  chi0s       Static chi0(q, 0), shape (nsc, nsc).
   !> @param[in]  kmat        Interaction kernel (Ry), shape (nsc, nsc).
   !> @param[in]  m           Signed moment of each row kappa, shape (nsc).
   !> @param[out] omega       Real parts of the eigenvalues of Omega (Ry), ascending, shape (nsc).
   !> @param[out] imag_max    Largest |imaginary part| of those eigenvalues (Ry).
   !> @param[out] omega_force Optional: real parts of the eigenvalues of the force-theorem matrix N K (I + chi0 K) (Ry), ascending,
   !>                         from the same chi0s; -chi0^-1 of Omega is replaced by K.
   module subroutine adiabatic_frequencies(chi0s, kmat, m, omega, imag_max, omega_force)
      complex(rp), intent(in) :: chi0s(:, :)
      real(rp), intent(in) :: kmat(:, :), m(:)
      real(rp), intent(out) :: omega(:), imag_max
      real(rp), intent(out), optional :: omega_force(:)

      integer :: n, i, info, ipiv(size(m))
      complex(rp) :: a(size(m), size(m)), b(size(m), size(m)), w(size(m)), vl(1, 1), vr(1, 1), work(4*size(m))
      real(rp) :: rwork(2*size(m))

      n = size(m)
      a = chi0s
      b = (0.0_rp, 0.0_rp)
      do i = 1, n
         b(i, i) = (1.0_rp, 0.0_rp)
      end do
      call zgesv(n, n, a, n, ipiv, b, n, info)
      if (info /= 0) error stop 'adiabatic_frequencies: singular static chi0'
      do i = 1, n
         a(i, :) = -m(i)*(b(i, :) + kmat(i, :))
      end do
      call zgeev('N', 'N', n, a, n, w, vl, 1, vr, 1, work, size(work), rwork, info)
      if (info /= 0) error stop 'adiabatic_frequencies: zgeev failed'
      omega = real(w, rp)
      call sort_ascending(omega)
      imag_max = maxval(abs(aimag(w)))
      if (present(omega_force)) then
         a = cmplx(kmat, 0.0_rp, rp)
         a = a + matmul(a, matmul(chi0s, a))
         do i = 1, n
            a(i, :) = m(i)*a(i, :)
         end do
         call zgeev('N', 'N', n, a, n, w, vl, 1, vr, 1, work, size(work), rwork, info)
         if (info /= 0) error stop 'adiabatic_frequencies: zgeev failed (force theorem)'
         omega_force = real(w, rp)
         call sort_ascending(omega_force)
      end if
   end subroutine adiabatic_frequencies

   !> @brief In-place ascending selection sort of a short list.
   pure subroutine sort_ascending(x)
      real(rp), intent(inout) :: x(:)

      integer :: i, j
      real(rp) :: t

      do i = 1, size(x) - 1
         j = i - 1 + minloc(x(i:), 1)
         t = x(i)
         x(i) = x(j)
         x(j) = t
      end do
   end subroutine sort_ascending

   !> @brief tr L(omega) with L = -(chi - chi^H)/(2 pi i), i.e. -Im tr chi / pi.
   pure module function spectral_trace(chi) result(trl)
      complex(rp), intent(in) :: chi(:, :, :)
      real(rp) :: trl(size(chi, 3))

      integer :: iw, i

      do iw = 1, size(chi, 3)
         trl(iw) = 0.0_rp
         do i = 1, size(chi, 1)
            trl(iw) = trl(iw) - aimag(chi(i, i, iw))/pi
         end do
      end do
   end function spectral_trace

   !> @brief One-site pole estimates: parabolic peak of tr L and first zero of Re(1 + chi0 U) with interpolated omega > 0.
   !> @details A maximum on either grid end gives has_peak = .false.; no sign change gives
   !>          has_crossing = .false. peak and crossing are undefined unless their flag is set.
   !> @param[in]  omega  Increasing grid, shape (nw).
   !> @param[in]  trl    tr L on the grid, shape (nw).
   !> @param[in]  chi0   Bare response on the grid, shape (nw).
   !> @param[in]  u      Interaction (Ry).
   pure module subroutine pole_estimates(omega, trl, chi0, u, peak, crossing, has_peak, has_crossing)
      real(rp), intent(in) :: omega(:), trl(:), u
      complex(rp), intent(in) :: chi0(:)
      real(rp), intent(out) :: peak, crossing
      logical, intent(out) :: has_peak, has_crossing

      integer :: iw
      real(rp) :: g0, g1, root

      crossing = 0.0_rp
      call peak_estimate(omega, trl, peak, has_peak)

      has_crossing = .false.
      do iw = 1, size(omega) - 1
         g0 = 1.0_rp + u*real(chi0(iw), rp)
         g1 = 1.0_rp + u*real(chi0(iw + 1), rp)
         if (g0*g1 < 0.0_rp) then
            root = omega(iw) - g0*(omega(iw + 1) - omega(iw))/(g1 - g0)
            if (root > 0.0_rp) then
               crossing = root
               has_crossing = .true.
               exit
            end if
         end if
      end do
   end subroutine pole_estimates

   !> @brief Parabolic peak of tr L; a maximum on either grid end gives has_peak = .false. and peak undefined (0).
   pure module subroutine peak_estimate(omega, trl, peak, has_peak)
      real(rp), intent(in) :: omega(:), trl(:)
      real(rp), intent(out) :: peak
      logical, intent(out) :: has_peak

      integer :: imax
      real(rp) :: h

      peak = 0.0_rp
      imax = maxloc(trl, 1)
      has_peak = imax > 1 .and. imax < size(trl)
      if (has_peak) then
         h = 0.5_rp*(omega(imax + 1) - omega(imax - 1))
         peak = omega(imax) + h*(trl(imax - 1) - trl(imax + 1)) &
                /(2.0_rp*(trl(imax - 1) - 2.0_rp*trl(imax) + trl(imax + 1)))
      end if
   end subroutine peak_estimate

end submodule spin_response_kernels
