submodule (exchange_mod) exchange_dynamics
   implicit none

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculation of the Gilbert damping parameters in real-space
   !> based on the LMTO solution, torque-correlation model following
   !> PRB 108, 014433 (2023), Eq. (15), in the T -> 0 limit:
   !>   alpha_ij^{mu nu} = alpha_C Tr[T_i^mu A_ij (T_j^nu)^dagger A_ji] at E_F,
   !> with A_ij = (G_ij - G_ji^dagger)/(2i) built from G - G^dagger (not Im G),
   !> alpha_C = g/(pi m_t) = 2/(pi m_spin) since g = 2(1 + m_orb/m_spin).
   !> E_F is the energy-grid point nearest to en%fermi.
   !> The following calculation works better with hoh.
   !> Implemented by Ivan Miranda on 14.11.2023
   !> Revised by Ivan Miranda on 29.04.2025
   !---------------------------------------------------------------------------
   function damping_ready(this) result(ready)
      class(exchange) :: this
      logical :: ready

      ready = this%control%nsp == 2
      if (.not. ready) then
         call g_logger%warning('calculate_gilbert_damping: needs nsp = 2 (collinear with l.s); skipped.', __FILE__, __LINE__)
         return
      end if
      call this%hamiltonian%torque_operator_collinear()
   end function damping_ready

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Pair loop of the exchange consumers (J/D/A, two-index J/D/A, Gilbert damping, moment of inertia, static
   !> transverse chi0 of spin_response_rs).
   !> @details Every enabled consumer runs its per-pair kernel on the Green functions of the pair, then its finalize
   !>          writes the outputs. With stream = .false. the Green functions of all local pairs are already filled
   !>          (k-space routes) and the kernels read slot = local pair index. With stream = .true. (real-space recursion
   !>          route) the Green-function arrays hold one slot: each pair is reconstructed into slot 1 before the kernels.
   !---------------------------------------------------------------------------
   module subroutine calculate_pairwise(this, stream, do_exchange, do_damping, do_inertia, do_chi0)
      use mpi_mod
      class(exchange) :: this
      logical, intent(in) :: stream, do_exchange, do_damping, do_inertia, do_chi0
      type(energy_nodes) :: nodes
      integer :: i, j, njij_glob, slot, ief, nenergy, nres, nb
      real(rp) :: step_size
      logical :: damping, inertia, moments
      complex(rp), allocatable :: h(:, :), gij0(:, :), gji0(:, :), gij2(:, :), gji2(:, :)
      real(rp), allocatable :: kernel(:)
      real(rp), allocatable :: xc(:, :), xcso(:, :), xcfo(:, :), xcparts(:, :)
      ! damping: dpair(1:9) Re tensor at E_F, dpair(10:18) Im tensor at E_F, dpair(19) prefactor
      real(rp), allocatable :: dpair(:, :), dtotal(:, :), dcontrib(:, :)
      real(rp), allocatable :: ipair(:, :), ipair_b(:, :, :), itotal(:, :), icontrib(:, :)
      complex(rp), allocatable :: chi0(:), lin(:)
      real(rp), allocatable :: md(:)

      damping = do_damping
      if (damping) damping = damping_ready(this)
      inertia = do_inertia
      moments = this%control%exchange_integration == 'moments'
      if (inertia) inertia = inertia_ready(this, step_size, moments)
      nenergy = size(this%en%ene)
      ief = fermi_index(this%en%ene, this%en%fermi)
      if (do_exchange) then
         ! row = [jij dij aij]; xcparts row = [jcd, jsd, jcc, jsc, dcc(3), dsc(3), isd(9), isc(9)]
         allocate (xc(13, this%lattice%njij), xcso(13, this%lattice%njij), xcfo(13, this%lattice%njij), &
                   xcparts(28, this%lattice%njij), source=0.0_rp)
         nodes = real_mesh_nodes(this)
      end if
      if (do_chi0) then
         call chi0_check(this)
         allocate (chi0(this%lattice%njij), source=(0.0_rp, 0.0_rp))
         allocate (md(this%lattice%njij), source=0.0_rp)
         if (.not. do_exchange) nodes = real_mesh_nodes(this)
      end if
      ! Moment space: the stream slot 1 holds the raw moments of the pair, slot 2 the transformed ones
      if (moments .and. do_chi0) then
         call this%green%chebyshev_exchange_weights(h, kernel, lin)
      else if (moments .and. do_exchange) then
         call this%green%chebyshev_exchange_weights(h, kernel)
      end if
      ! Energy-resolved arrays: empty in moment space, where the energy-resolved files are not written
      nres = merge(0, nenergy, moments)
      if (damping) allocate (dtotal(9, nres), dcontrib(9, nres), dpair(19, this%lattice%njij), source=0.0_rp)
      if (inertia) allocate (itotal(9, nres), icontrib(9, nres), ipair_b(4, nres, this%lattice%njij), &
                             ipair(18, this%lattice%njij), source=0.0_rp)

      if (moments .and. (damping .or. inertia)) then
         nb = size(this%green%gij, 1)
         allocate (gij0(nb, nb), gji0(nb, nb), gij2(nb, nb), gji2(nb, nb))
      end if
      if (stream) call this%green%prepare_intersite_gf()
      ! Each rank holds the pairs start_atom..end_atom (local index g2l_map)
      do njij_glob = start_atom, end_atom
         i = this%lattice%ijpair(njij_glob, 1) ! Atom number in the clust file, atom i
         j = this%lattice%ijpair(njij_glob, 2) ! Atom number in the clust file, atom j
         slot = g2l_map(njij_glob)
         if (moments) then
            slot = 1
            call this%green%fill_intersite_moments_pair(njij_glob, slot)
         else if (stream) then
            slot = 1
            call this%green%fill_intersite_gf_pair(njij_glob, slot)
            call this%green%fill_intersite_gf_twoindex_pair(slot)
         end if
         ! Before the exchange kernels: in moment space they overwrite slot 2, which chi0_pair fills with the p = 0 transform
         if (do_chi0) then
            if (moments) then
               call chi0_pair_moments(this, h, kernel, lin, i == j, chi0(njij_glob), md(njij_glob))
            else
               call chi0_pair_mesh(this, nodes, slot, i == j, chi0(njij_glob), md(njij_glob))
            end if
         end if
         if (do_exchange .and. moments) then
            call moment_exchange_pair(this, h, kernel, i, j, xc(:, njij_glob), xcso(:, njij_glob), xcfo(:, njij_glob), &
                                      xcparts(:, njij_glob))
         else if (do_exchange) then
            call exchange_pair(this, nodes, slot, i, j, xc(:, njij_glob))
            call exchange_twoindex_pair(this, nodes, slot, i, j, xcso(:, njij_glob), xcfo(:, njij_glob), xcparts(:, njij_glob))
         end if
         if (moments .and. (damping .or. inertia)) then
            call this%green%chebyshev_point_pair(njij_glob, gij0, gji0, gij2, gji2)
            if (damping) call damping_point_pair(this, i, j, gij0, gji0, dpair(:, njij_glob))
            if (inertia) call inertia_point_pair(this, i, j, gij0, gji0, gij2, gji2, ipair(:, njij_glob))
         else
            if (damping) then
               call damping_pair(this, slot, i, j, ief, dpair(:, njij_glob), dcontrib)
               dtotal = dtotal + dcontrib
            end if
            if (inertia) then
               call inertia_pair(this, slot, i, j, ief, step_size, ipair(:, njij_glob), ipair_b(:, :, njij_glob), icontrib)
               itotal = itotal + icontrib
            end if
         end if
      end do

      if (do_exchange) then
         call exchange_finalize(this, xc)
         call exchange_twoindex_finalize(this, xcso, xcfo, xcparts)
      end if
      if (do_chi0) call chi0_finalize(this, chi0, md)
      if (damping) call damping_finalize(this, dpair, dtotal, ief, moments)
      if (inertia) call inertia_finalize(this, ipair, ipair_b, itotal, ief, moments)
   end subroutine calculate_pairwise

   !> Stop unless the run is collinear (nsp = 1, or nsp = 2 with every moment along +-z) and has a d channel: the spin
   !> blocks of G are then the up and down states of the transverse response of spin_response_rs.
   subroutine chi0_check(this)
      class(exchange) :: this
      integer :: it
      real(rp) :: m(3)

      if (this%control%nsp /= 1 .and. this%control%nsp /= 2) &
         call g_logger%fatal('do_chi0 needs nsp = 1 or 2 (collinear, with or without SOC)', __FILE__, __LINE__)
      if (lmax_basis < 2) call g_logger%fatal('do_chi0 needs lmax >= 2 (d channel)', __FILE__, __LINE__)
      if (this%control%nsp == 1) return
      do it = 1, this%lattice%ntype
         m = this%symbolic_atom(it)%potential%mom
         if (norm2(m(1:2)) > 1.0e-8_rp*norm2(m)) &
            call g_logger%fatal('do_chi0 with nsp = 2 needs every moment along +-z (type '//int2str(it)//')', __FILE__, __LINE__)
      end do
   end subroutine chi0_check

   !> chi0_ij of the pair in Green-function slot `slot` on the real energy mesh (Simpson with the Fermi function), and for a
   !> self pair its d moment md = -(1/pi) Im int f Tr_d[G^upup - G^dndn] dE (0 otherwise).
   subroutine chi0_pair_mesh(this, nodes, slot, self_pair, chi0, md)
      use spin_response_rs_mod, only: chi0_traces, chi0_from_traces, d_moment_trace
      class(exchange) :: this
      type(energy_nodes), intent(in) :: nodes
      integer, intent(in) :: slot
      logical, intent(in) :: self_pair
      complex(rp), intent(out) :: chi0
      real(rp), intent(out) :: md
      complex(rp) :: t(size(nodes%z), 2)
      real(rp) :: re(2), im(2)
      integer :: k, n

      n = size(nodes%z)
      t = chi0_traces(this%green%gij(:, :, 1:n, slot), this%green%gji(:, :, 1:n, slot))
      do k = 1, 2
         call nodes%integrate(real(t(:, k), rp), re(k))
         call nodes%integrate(aimag(t(:, k)), im(k))
      end do
      chi0 = chi0_from_traces(cmplx(re(1), im(1), rp), cmplx(re(2), im(2), rp))
      md = 0.0_rp
      if (self_pair) then
         call nodes%integrate(aimag(d_moment_trace(this%green%gij(:, :, 1:n, slot))), md)
         md = -md/pi
      end if
   end subroutine chi0_pair_mesh

   !> chi0_ij from the raw moments in slot 1 and their p = 0 bilinear transform, written to slot 2:
   !> int f Tr[A(E) B(E)] dE = sum_m Tr[A_m B~_m] (see moment_exchange_pair). For a self pair md is its d moment from the
   !> linear weights lin (chebyshev_exchange_weights), 0 otherwise.
   subroutine chi0_pair_moments(this, h, kernel, lin, self_pair, chi0, md)
      use spin_response_rs_mod, only: chi0_traces, chi0_from_traces, d_moment_trace
      class(exchange) :: this
      complex(rp), intent(in) :: h(0:, 0:), lin(:)
      real(rp), intent(in) :: kernel(:)
      logical, intent(in) :: self_pair
      complex(rp), intent(out) :: chi0
      real(rp), intent(out) :: md
      complex(rp), allocatable :: t(:, :)
      integer :: n

      n = size(kernel)
      md = 0.0_rp
      if (self_pair) md = -aimag(sum(lin*d_moment_trace(this%green%gij(:, :, 1:n, 1))))/pi
      call this%green%transform_intersite_moments(1, 2, h(:, 0), kernel)
      t = chi0_traces(this%green%gij(:, :, 1:n, 1), this%green%gji(:, :, 1:n, 2))
      chi0 = chi0_from_traces(sum(t(:, 1)), sum(t(:, 2)))
   end subroutine chi0_pair_moments

   !> Sum chi0 and md over the ranks; on rank 0 write chi0_ij.dat and the type diagnostics chi0_types.dat.
   subroutine chi0_finalize(this, chi0, md)
      use mpi_mod
      use spin_response_rs_mod, only: write_chi0_pairs, type_sums, write_chi0_types
      class(exchange) :: this
      complex(rp), intent(inout) :: chi0(:)
      real(rp), intent(inout) :: md(:)
      integer :: nt, it
      integer, allocatable :: seed(:), npair(:)
      real(rp), allocatable :: m(:), chi0s(:, :), c_up(:), c_dn(:)

#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, chi0, size(chi0), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, md, size(md), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
      if (rank /= 0) return
      call write_chi0_pairs('chi0_ij.dat', this%lattice%ijpair(1:size(chi0), :), this%lattice%iz, this%lattice%cr, chi0)
      nt = this%lattice%ntype
      allocate (seed(nt), npair(nt), m(nt), chi0s(nt, nt), c_up(nt), c_dn(nt))
      call type_sums(this%lattice%ijpair(1:size(chi0), :), this%lattice%iz, nt, chi0, md, seed, npair, m, chi0s)
      do it = 1, nt
         c_up(it) = this%symbolic_atom(it)%potential%center_band(3, 1)
         c_dn(it) = this%symbolic_atom(it)%potential%center_band(3, 2)
      end do
      call write_chi0_types('chi0_types.dat', seed, npair, m, chi0s, c_up, c_dn)
   end subroutine chi0_finalize

   !> J, D, A of pair (i, j) from the raw moments in Green-function slot 1 (see fill_intersite_moments_pair). With
   !> Delta = alpha + beta E the integral is the sum over the E^p terms (alpha_i, alpha_j; beta_i, alpha_j and alpha_i, beta_j;
   !> beta_i, beta_j) of sum_m Tr[Di X_m Dj Y_m], X the raw moments (slot 1, i side) and Y the p-transformed ones (slot 2, j side).
   !> The sums are contracted as the 12 x 12 family-pair traces T(X, Y); every output is a combination of those.
   subroutine moment_exchange_pair(this, h, kernel, i, j, xc, xcso, xcfo, xcparts)
      class(exchange) :: this
      complex(rp), intent(in) :: h(0:, 0:)
      real(rp), intent(in) :: kernel(:)
      integer, intent(in) :: i, j
      real(rp), intent(out) :: xc(13), xcso(13), xcfo(13), xcparts(28)
      ! Family index: 1 Gnmag, 1 + k Gk (k = x, y, z), 5 g00, 6 g01, 6 + k gk0, 9 + k gk1
      integer, parameter :: c0 = 5, c1 = 6, s0 = 6, s1 = 9
      real(rp), parameter :: f = 1.0e3_rp/(4.0_rp*pi)
      complex(rp) :: alpha(9, 9, 2), beta(9, 9, 2), tt(12, 12)
      complex(rp), allocatable :: xf(:, :, :, :), yf(:, :, :, :), yt(:, :), z(:, :)
      real(rp) :: jcd, jcc, jsd, jsc, dsc, dcc, isd, isc
      integer :: p, ti, tj, nmom, m, q, k, l, idx

      nmom = size(kernel)
      allocate (xf(9, 9, nmom, 12), yf(9, 9, nmom, 12), yt(81*nmom, 12), z(12, 81*nmom))
      beta = (0.0_rp, 0.0_rp)
      if (this%use_onsite_delta) then
         alpha(:, :, 1) = onsite_exchange_delta(this, i)
         alpha(:, :, 2) = onsite_exchange_delta(this, j)
      else
         call this%symbolic_atom(this%lattice%iz(i))%d_matrix(alpha(:, :, 1), 0.0_rp)
         call this%symbolic_atom(this%lattice%iz(j))%d_matrix(alpha(:, :, 2), 0.0_rp)
         call this%symbolic_atom(this%lattice%iz(i))%d_matrix(beta(:, :, 1), 1.0_rp)
         call this%symbolic_atom(this%lattice%iz(j))%d_matrix(beta(:, :, 2), 1.0_rp)
         beta = beta - alpha
      end if
      call stack_families(this, 1, .false., xf)
      tt = (0.0_rp, 0.0_rp)
      do p = 0, merge(0, 2, this%use_onsite_delta)
         call this%green%transform_intersite_moments(1, 2, h(:, p), kernel)
         call stack_families(this, 2, .true., yf)
         ! Tr[A X B Y] = sum_kl (A X B)(k, l) Y(l, k): z is (A X_m B)(k, l), yt is Y_m(l, k), both flattened over (k, l, m)
         do q = 1, 12
            do m = 1, nmom
               yt(81*(m - 1) + 1:81*m, q) = reshape(transpose(yf(:, :, m, q)), [81])
            end do
         end do
         do tj = 0, 1
            ti = p - tj
            if (ti < 0 .or. ti > 1) cycle
            do q = 1, 12
               do m = 1, nmom
                  z(q, 81*(m - 1) + 1:81*m) = reshape(matmul(merge(beta(:, :, 1), alpha(:, :, 1), ti == 1), &
                                                             matmul(xf(:, :, m, q), merge(beta(:, :, 2), alpha(:, :, 2), tj == 1))), [81])
               end do
            end do
            tt = tt + matmul(z, yt)
         end do
      end do

      ! J, D, A
      xc(1) = f*aimag(tt(1, 1) - tt(2, 2) - tt(3, 3) - tt(4, 4))
      do k = 1, 3
         xc(1 + k) = f*real(tt(1, 1 + k) - tt(1 + k, 1))
         do l = 1, 3
            xc(4 + k + 3*(l - 1)) = f*0.5_rp*aimag(tt(1 + k, 1 + l) + tt(1 + l, 1 + k))
         end do
      end do
      ! Two-index J (charge/spin density and current parts), D, A
      jcd = aimag(tt(c0, c0)); jcc = aimag(tt(c1, c1))
      jsd = aimag(tt(s0 + 1, s0 + 1) + tt(s0 + 2, s0 + 2) + tt(s0 + 3, s0 + 3))
      jsc = aimag(tt(s1 + 1, s1 + 1) + tt(s1 + 2, s1 + 2) + tt(s1 + 3, s1 + 3))
      xcso(1) = f*(jcd - jsd + jcc - jsc)
      xcfo(1) = f*(jcd + jsd - jcc - jsc)
      xcparts(1:4) = f*[jcd, jsd, jcc, jsc]
      do k = 1, 3
         dsc = real(tt(c0, s1 + k)); dcc = real(tt(c1, s0 + k))
         xcso(1 + k) = f*2.0_rp*(dsc + dcc)
         xcfo(1 + k) = f*2.0_rp*(dsc - dcc)
         xcparts(4 + k) = f*2.0_rp*dcc
         xcparts(7 + k) = f*2.0_rp*dsc
         do l = 1, 3
            idx = k + 3*(l - 1)
            isd = aimag(tt(s0 + k, s0 + l)); isc = aimag(tt(s1 + k, s1 + l))
            xcso(4 + idx) = f*(isd + isc)
            xcfo(4 + idx) = f*(isc - isd)
            xcparts(10 + idx) = f*isd
            xcparts(19 + idx) = f*isc
         end do
      end do
   end subroutine moment_exchange_pair

   !> The 12 families of Green-function slot `slot` (i side and "ij", or with j_side the j side and "ji") for moments 1:size(f, 3).
   subroutine stack_families(this, slot, j_side, f)
      class(exchange) :: this
      integer, intent(in) :: slot
      logical, intent(in) :: j_side
      complex(rp), intent(out) :: f(:, :, :, :)
      integer :: n

      n = size(f, 3)
      associate (g => this%green)
         if (j_side) then
            f(:, :, :, 1) = g%gjnmag(1:9, 1:9, 1:n, slot); f(:, :, :, 2) = g%gjx(1:9, 1:9, 1:n, slot)
            f(:, :, :, 3) = g%gjy(1:9, 1:9, 1:n, slot); f(:, :, :, 4) = g%gjz(1:9, 1:9, 1:n, slot)
            f(:, :, :, 5) = g%g00ji(1:9, 1:9, 1:n, slot); f(:, :, :, 6) = g%g01ji(1:9, 1:9, 1:n, slot)
            f(:, :, :, 7) = g%gx0ji(1:9, 1:9, 1:n, slot); f(:, :, :, 8) = g%gy0ji(1:9, 1:9, 1:n, slot)
            f(:, :, :, 9) = g%gz0ji(1:9, 1:9, 1:n, slot); f(:, :, :, 10) = g%gx1ji(1:9, 1:9, 1:n, slot)
            f(:, :, :, 11) = g%gy1ji(1:9, 1:9, 1:n, slot); f(:, :, :, 12) = g%gz1ji(1:9, 1:9, 1:n, slot)
         else
            f(:, :, :, 1) = g%ginmag(1:9, 1:9, 1:n, slot); f(:, :, :, 2) = g%gix(1:9, 1:9, 1:n, slot)
            f(:, :, :, 3) = g%giy(1:9, 1:9, 1:n, slot); f(:, :, :, 4) = g%giz(1:9, 1:9, 1:n, slot)
            f(:, :, :, 5) = g%g00ij(1:9, 1:9, 1:n, slot); f(:, :, :, 6) = g%g01ij(1:9, 1:9, 1:n, slot)
            f(:, :, :, 7) = g%gx0ij(1:9, 1:9, 1:n, slot); f(:, :, :, 8) = g%gy0ij(1:9, 1:9, 1:n, slot)
            f(:, :, :, 9) = g%gz0ij(1:9, 1:9, 1:n, slot); f(:, :, :, 10) = g%gx1ij(1:9, 1:9, 1:n, slot)
            f(:, :, :, 11) = g%gy1ij(1:9, 1:9, 1:n, slot); f(:, :, :, 12) = g%gz1ij(1:9, 1:9, 1:n, slot)
         end if
      end associate
   end subroutine stack_families

   !> Torque-correlation damping tensor of pair (i, j) read from Green-function slot `slot`: pair(1:9) Re and pair(10:18) Im
   !> at energy index ief, pair(19) prefactor; contrib(:, nv) is the pair's energy-resolved contribution.
   subroutine damping_pair(this, slot, i, j, ief, pair, contrib)
      class(exchange) :: this
      integer, intent(in) :: slot, i, j, ief
      real(rp), intent(out) :: pair(19), contrib(:, :)
      complex(rp), allocatable :: Aij(:, :), Aji(:, :), tjd(:, :, :)
      complex(rp) :: c(3, 3)
      integer :: k, l, m, nv, nb_ij, nenergy

      nenergy = size(this%en%ene)
      nb_ij = size(this%green%gij, 1)
      allocate (Aij(nb_ij, nb_ij), Aji(nb_ij, nb_ij))
      call damping_setup(this, i, j, tjd, pair(19))
      do nv = 1, nenergy
         ! Anti-Hermitian parts of the intersite Green functions (times 2i)
         Aij = this%green%gij(:, :, nv, slot) - transpose(conjg(this%green%gji(:, :, nv, slot)))
         Aji = this%green%gji(:, :, nv, slot) - transpose(conjg(this%green%gij(:, :, nv, slot)))
         c = torque_trace(this%hamiltonian%tmat(:, :, :, this%lattice%iz(i)), Aij, tjd, Aji)
         m = 0
         do k = 1, 3
            do l = 1, 3
               m = m + 1
               contrib(m, nv) = pair(19)*real(c(k, l), rp)
               if (nv == ief) then
                  pair(m) = real(c(k, l), rp)
                  pair(9 + m) = aimag(c(k, l))
               end if
            end do
         end do
      end do
   end subroutine damping_pair

   !> Torque-correlation damping tensor of pair (i, j) at E_F from the point values of G_ij, G_ji (moment space).
   subroutine damping_point_pair(this, i, j, gij0, gji0, pair)
      class(exchange) :: this
      integer, intent(in) :: i, j
      complex(rp), intent(in) :: gij0(:, :), gji0(:, :)
      real(rp), intent(out) :: pair(19)
      complex(rp), allocatable :: tjd(:, :, :)
      complex(rp) :: c(3, 3)

      call damping_setup(this, i, j, tjd, pair(19))
      c = torque_trace(this%hamiltonian%tmat(:, :, :, this%lattice%iz(i)), gij0 - transpose(conjg(gji0)), tjd, &
                       gji0 - transpose(conjg(gij0)))
      pair(1:9) = real(reshape(transpose(c), [9]), rp)
      pair(10:18) = aimag(reshape(transpose(c), [9]))
   end subroutine damping_point_pair

   !> (T_j^mu)^dagger (nb, nb, 3) and the damping prefactor -(1/4) 2/(pi m_spin) of atom i.
   subroutine damping_setup(this, i, j, tjd, prefactor)
      class(exchange) :: this
      integer, intent(in) :: i, j
      complex(rp), allocatable, intent(out) :: tjd(:, :, :)
      real(rp), intent(out) :: prefactor
      integer :: k, l
      real(rp) :: spin_i

      allocate (tjd(size(this%green%gij, 1), size(this%green%gij, 1), 3))
      ! Spin magnetic moment of atom i
      spin_i = 0.0_rp
      do k = 0, this%symbolic_atom(this%lattice%iz(i))%potential%lmax
         spin_i = spin_i + this%symbolic_atom(this%lattice%iz(i))%potential%ql(1, k, 1) - &
                  this%symbolic_atom(this%lattice%iz(i))%potential%ql(1, k, 2)
      end do
      prefactor = (-0.25_rp)*(2.0_rp/(pi*spin_i))
      do l = 1, 3
         tjd(:, :, l) = transpose(conjg(this%hamiltonian%tmat(:, :, l, this%lattice%iz(j))))
      end do
   end subroutine damping_setup

   subroutine damping_finalize(this, pair, total_damping, ief, moments)
      use mpi_mod
      class(exchange) :: this
      real(rp), intent(inout) :: pair(:, :), total_damping(:, :)
      integer, intent(in) :: ief
      logical, intent(in) :: moments
      integer :: i, j, njij_glob, nv, nenergy, unit_e, unit_p
      real(rp) :: factor, distance_alat, alpha(3, 3), mhat(3)

      nenergy = size(this%en%ene)
#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, total_damping, size(total_damping), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, pair, size(pair), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      if (rank == 0) then
         if (.not. moments) open (newunit=unit_e, FILE='damping-energy.out', STATUS='replace', ACTION='write')
         open (newunit=unit_p, FILE='alldampings.out', STATUS='replace', ACTION='write')
         write (unit_p, '(A)') '#     i      j  alpha: xx xy xz yx yy yz zx zy zz transverse (dimensionless)'// &
            '  dist rij=ri-rj (alat)'
         if (moments) then
            write (*, *) 'Fermi energy used = ', this%en%fermi
         else
            write (*, *) 'ief = ', this%en%ene(ief), 'fermi = ', this%en%fermi
         end if
         do njij_glob = 1, this%lattice%njij
            i = this%lattice%ijpair(njij_glob, 1)
            j = this%lattice%ijpair(njij_glob, 2)
            distance_alat = norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            factor = pair(19, njij_glob)
            write (*, '(A,I6,A,I6,A)') 'Damping tensor between pair', i, ' and ', j, ' is'
            write (*, '(3ES16.8)') factor*pair(1:9, njij_glob)
            write (*, '(A,I6,A,I6,A)') 'Imaginary part of the tensor between pair', i, ' and ', j, ' is'
            write (*, '(3ES16.8)') factor*pair(10:18, njij_glob)
            write (*, *) '----------------------'
            ! Transverse damping 0.5*(Tr alpha - m.alpha.m), m the moment direction of atom i
            ! (0.5*(xx+yy) for m along z); alpha(l, k) holds component kl, the scalar is unchanged.
            alpha = reshape(factor*pair(1:9, njij_glob), [3, 3])
            mhat = this%symbolic_atom(this%lattice%iz(i))%potential%mom
            write (unit_p, '(2I7,10ES16.8,4F11.6)') i, j, factor*pair(1:9, njij_glob), &
               0.5_rp*(alpha(1, 1) + alpha(2, 2) + alpha(3, 3) - dot_product(mhat, matmul(alpha, mhat))), &
               distance_alat, this%lattice%cr(:, i) - this%lattice%cr(:, j)
         end do

         ! Energy-resolved total: sum over all listed pairs, each with its own prefactor
         if (.not. moments) then
            write (unit_e, '(A)') '# E-Ef (Ry)  alpha: xx xy xz yx yy yz zx zy zz (dimensionless)'
            do nv = 1, nenergy
               write (unit_e, '(F14.9,9ES16.8)') this%en%ene(nv) - this%en%fermi, total_damping(:, nv)
            end do
            close (unit_e)
         end if
         close (unit_p)
      end if

   end subroutine damping_finalize

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculation of the moment of inertia parameters in real-space
   !> based on the LMTO solution (see, e.g., Sci. Rep. 7, 931 (2017)).
   !> The derivation in real-space is new and experimental by Nov. 2023.
   !> The model is based on the so-called torque-correlation
   !> The following calculation works better with hoh.
   !> Implemented by Ivan Miranda on 22.11.2023
   !---------------------------------------------------------------------------
   function inertia_ready(this, step_size, moments) result(ready)
      class(exchange) :: this
      real(rp), intent(out) :: step_size
      logical, intent(in) :: moments
      logical :: ready
      integer :: nenergy

      ! This is the experimental torque-correlation kernel described in
      ! Sci. Rep. 7, 931 (2017).  It is deliberately reported as the raw
      ! production kernel: unlike Gilbert damping above, this routine has no
      ! independently established normalization or production relation.
      ready = .false.
      step_size = 0.0_rp
      if (this%control%nsp /= 2) then
         call g_logger%warning('calculate_moment_of_inertia: needs nsp = 2 (collinear with l.s); skipped.', __FILE__, __LINE__)
         return
      end if
      nenergy = size(this%en%ene)
      if (moments) then
         call this%hamiltonian%torque_operator_collinear()
         ready = .true.
         return
      end if
      if (nenergy < 3) then
         write (*, *) 'Moment of inertia skipped: at least three energy points are required.'
         return
      end if
      step_size = this%en%ene(2) - this%en%ene(1)
      if (step_size .eq. 0.0_rp) then
         write (*, *) 'Moment of inertia skipped: zero energy step.'
         return
      end if
      call this%hamiltonian%torque_operator_collinear()
      ready = .true.
   end function inertia_ready

   !> Raw inertia kernel of pair (i, j) read from Green-function slot `slot`: pair(1:9) Re and pair(10:18) Im at energy index ief,
   !> pair_b the Bij(1, 1) trace and its second derivative, contrib(:, nv) the energy-resolved real part.
   subroutine inertia_pair(this, slot, i, j, ief, step_size, pair, pair_b, contrib)
      class(exchange) :: this
      integer, intent(in) :: slot, i, j, ief
      real(rp), intent(in) :: step_size
      real(rp), intent(out) :: pair(18), pair_b(:, :), contrib(:, :)
      complex(rp), allocatable :: Aij(:, :, :), Aji(:, :, :), Bij(:, :, :), Bji(:, :, :)
      complex(rp), allocatable :: sBij(:, :, :), sBji(:, :, :)
      complex(rp), allocatable :: tjd(:, :, :)
      complex(rp) :: c(3, 3)
      integer :: k, l, m, nv, nenergy, nb_ij

      nenergy = size(this%en%ene)
      nb_ij = size(this%green%gij, 1)
      allocate (Aij(nb_ij, nb_ij, nenergy), Aji(nb_ij, nb_ij, nenergy), Bij(nb_ij, nb_ij, nenergy), Bji(nb_ij, nb_ij, nenergy))
      allocate (sBij(nb_ij, nb_ij, nenergy), sBji(nb_ij, nb_ij, nenergy))
      allocate (tjd(nb_ij, nb_ij, 3))

      do nv = 1, nenergy
         Aij(:, :, nv) = this%green%gij(:, :, nv, slot) - transpose(conjg(this%green%gji(:, :, nv, slot)))
         Aji(:, :, nv) = this%green%gji(:, :, nv, slot) - transpose(conjg(this%green%gij(:, :, nv, slot)))
         Bij(:, :, nv) = this%green%gij(:, :, nv, slot) + transpose(conjg(this%green%gji(:, :, nv, slot)))
         Bji(:, :, nv) = this%green%gji(:, :, nv, slot) + transpose(conjg(this%green%gij(:, :, nv, slot)))
      end do
      sBij = czero; sBji = czero
      do k = 1, nb_ij
         do l = 1, nb_ij
            sBij(k, l, 2:nenergy - 1) = cmplx(second_derivative(real(Bij(k, l, :)), step_size), &
                                              second_derivative(aimag(Bij(k, l, :)), step_size), kind=rp)
            sBji(k, l, 2:nenergy - 1) = cmplx(second_derivative(real(Bji(k, l, :)), step_size), &
                                              second_derivative(aimag(Bji(k, l, :)), step_size), kind=rp)
         end do
      end do

      pair_b(1, :) = real(Bij(1, 1, :)); pair_b(2, :) = real(sBij(1, 1, :))
      pair_b(3, :) = aimag(Bij(1, 1, :)); pair_b(4, :) = aimag(sBij(1, 1, :))
      do l = 1, 3
         tjd(:, :, l) = transpose(conjg(this%hamiltonian%tmat(:, :, l, this%lattice%iz(j))))
      end do
      do nv = 1, nenergy
         c = torque_trace(this%hamiltonian%tmat(:, :, :, this%lattice%iz(i)), Aij(:, :, nv), tjd, sBji(:, :, nv)) + &
             torque_trace(this%hamiltonian%tmat(:, :, :, this%lattice%iz(i)), sBij(:, :, nv), tjd, Aji(:, :, nv))
         m = 0
         do k = 1, 3
            do l = 1, 3
               m = m + 1
               contrib(m, nv) = real(c(k, l), rp)
               if (nv == ief) then
                  pair(m) = real(c(k, l), rp)
                  pair(9 + m) = aimag(c(k, l))
               end if
            end do
         end do
      end do
   end subroutine inertia_pair

   !> Raw inertia kernel of pair (i, j) at E_F from the point values of G and d2G/dE2 (moment space).
   subroutine inertia_point_pair(this, i, j, gij0, gji0, gij2, gji2, pair)
      class(exchange) :: this
      integer, intent(in) :: i, j
      complex(rp), intent(in) :: gij0(:, :), gji0(:, :), gij2(:, :), gji2(:, :)
      real(rp), intent(out) :: pair(18)
      complex(rp), allocatable :: tjd(:, :, :)
      complex(rp) :: c(3, 3), Aij(size(gij0, 1), size(gij0, 1)), Aji(size(gij0, 1), size(gij0, 1))
      real(rp) :: unused

      call damping_setup(this, i, j, tjd, unused)
      Aij = gij0 - transpose(conjg(gji0))
      Aji = gji0 - transpose(conjg(gij0))
      c = torque_trace(this%hamiltonian%tmat(:, :, :, this%lattice%iz(i)), Aij, tjd, gji2 + transpose(conjg(gij2))) + &
          torque_trace(this%hamiltonian%tmat(:, :, :, this%lattice%iz(i)), gij2 + transpose(conjg(gji2)), tjd, Aji)
      pair(1:9) = real(reshape(transpose(c), [9]), rp)
      pair(10:18) = aimag(reshape(transpose(c), [9]))
   end subroutine inertia_point_pair

   subroutine inertia_finalize(this, pair, pair_b, total_inertia, ief, moments)
      use mpi_mod
      class(exchange) :: this
      real(rp), intent(inout) :: pair(:, :), pair_b(:, :, :), total_inertia(:, :)
      integer, intent(in) :: ief
      logical, intent(in) :: moments
      integer :: i, j, njij_glob, nv, nenergy, u_re, u_im, u_e, u_p
      real(rp) :: distance_alat

      nenergy = size(this%en%ene)
#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, total_inertia, size(total_inertia), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, pair_b, size(pair_b), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, pair, size(pair), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      if (rank == 0) then
         open (newunit=u_p, FILE='allinertias.out', STATUS='replace', ACTION='write')
         if (.not. moments) then
            open (newunit=u_re, FILE='example-real.out', STATUS='replace', ACTION='write')
            open (newunit=u_im, FILE='example-imag.out', STATUS='replace', ACTION='write')
            open (newunit=u_e, FILE='inertia-energy.out', STATUS='replace', ACTION='write')
            write (u_re, '(A)') '# i j E Re(Bij_11) Re(d2Bij_11/dE2)'
            write (u_im, '(A)') '# i j E Im(Bij_11) Im(d2Bij_11/dE2)'
            write (u_e, '(A)') '# E-Ef total raw inertia kernel: 9 real components'
         end if
         write (u_p, '(A)') '# i j 9 real components 9 imaginary components distance'
         do njij_glob = 1, this%lattice%njij
            i = this%lattice%ijpair(njij_glob, 1)
            j = this%lattice%ijpair(njij_glob, 2)
            if (.not. moments) then
               do nv = 1, nenergy
                  write (u_re, '(2I7,3ES18.8)') i, j, this%en%ene(nv), pair_b(1, nv, njij_glob), pair_b(2, nv, njij_glob)
                  write (u_im, '(2I7,3ES18.8)') i, j, this%en%ene(nv), pair_b(3, nv, njij_glob), pair_b(4, nv, njij_glob)
               end do
            end if
            distance_alat = norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            write (u_p, '(2I7,18ES18.8,F12.6)') i, j, pair(:, njij_glob), distance_alat
            write (*, '(A,I6,A,I6,A,ES14.6)') 'Raw inertia kernel pair ', i, ' and ', j, ' at Ef = ', merge(this%en%fermi, this%en%ene(ief), moments)
            write (*, '(3ES14.6)') pair(1:3, njij_glob)
            write (*, '(3ES14.6)') pair(4:6, njij_glob)
            write (*, '(3ES14.6)') pair(7:9, njij_glob)
         end do
         if (.not. moments) then
            do nv = 1, nenergy
               write (u_e, '(10ES18.8)') this%en%ene(nv) - this%en%fermi, total_inertia(:, nv)
            end do
            close (u_re); close (u_im); close (u_e)
         end if
         close (u_p)
      end if
   end subroutine inertia_finalize

   !> Index of the energy-grid point closest to the Fermi level.
   pure integer function fermi_index(ene, fermi) result(ief)
      real(rp), intent(in) :: ene(:), fermi

      ief = minloc(abs(ene - fermi), dim=1)
   end function fermi_index

   !> c(k, l) = Tr[ti(:, :, k) x tjd(:, :, l) y] for k, l = 1..3, using
   !> Tr(P Q) = sum(P * transpose(Q)) so only six matrix products are formed.
   pure function torque_trace(ti, x, tjd, y) result(c)
      complex(rp), intent(in) :: ti(:, :, :), x(:, :), tjd(:, :, :), y(:, :)
      complex(rp) :: c(3, 3)
      complex(rp) :: p(size(x, 1), size(x, 2), 3), q(size(y, 1), size(y, 2), 3)
      integer :: k, l

      do k = 1, 3
         p(:, :, k) = matmul(ti(:, :, k), x)
         q(:, :, k) = transpose(matmul(tjd(:, :, k), y))
      end do
      do l = 1, 3
         do k = 1, 3
            c(k, l) = sum(p(:, :, k)*q(:, :, l))
         end do
      end do
   end function torque_trace

end submodule exchange_dynamics
