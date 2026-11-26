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
   !> Pair loop of the exchange consumers (J/D/A, two-index J/D/A, Gilbert damping, moment of inertia).
   !> @details Every enabled consumer runs its per-pair kernel on the Green functions of the pair, then its finalize
   !>          writes the outputs. With stream = .false. the Green functions of all local pairs are already filled
   !>          (k-space routes) and the kernels read slot = local pair index. With stream = .true. (real-space recursion
   !>          route) the Green-function arrays hold one slot: each pair is reconstructed into slot 1 before the kernels.
   !---------------------------------------------------------------------------
   module subroutine calculate_pairwise(this, stream, do_exchange, do_damping, do_inertia)
      use mpi_mod
      class(exchange) :: this
      logical, intent(in) :: stream, do_exchange, do_damping, do_inertia
      type(energy_nodes) :: nodes
      integer :: i, j, njij_glob, slot, ief, nenergy
      real(rp) :: step_size
      logical :: damping, inertia, moments
      complex(rp), allocatable :: h(:, :)
      real(rp), allocatable :: kernel(:)
      real(rp), allocatable :: xc(:, :), xcso(:, :), xcfo(:, :), xcparts(:, :)
      ! damping: dpair(1:9) Re tensor at E_F, dpair(10:18) Im tensor at E_F, dpair(19) prefactor
      real(rp), allocatable :: dpair(:, :), dtotal(:, :), dcontrib(:, :)
      real(rp), allocatable :: ipair(:, :), ipair_b(:, :, :), itotal(:, :), icontrib(:, :)

      damping = do_damping
      if (damping) damping = damping_ready(this)
      inertia = do_inertia
      if (inertia) inertia = inertia_ready(this, step_size)
      nenergy = size(this%en%ene)
      ief = fermi_index(this%en%ene, this%en%fermi)
      if (do_exchange) then
         ! row = [jij dij aij]; xcparts row = [jcd, jsd, jcc, jsc, dcc(3), dsc(3), isd(9), isc(9)]
         allocate (xc(13, this%lattice%njij), xcso(13, this%lattice%njij), xcfo(13, this%lattice%njij), &
                   xcparts(28, this%lattice%njij), source=0.0_rp)
         nodes = real_mesh_nodes(this)
      end if
      ! Moment space: one node per Chebyshev moment; the stream slot 1 holds the raw moments of the pair, slot 2 the transformed ones
      moments = this%control%exchange_integration == 'moments'
      if (moments .and. do_exchange) then
         call this%green%chebyshev_exchange_weights(h, kernel)
         nodes = moment_nodes(size(kernel))
      end if
      if (damping) allocate (dtotal(9, nenergy), dpair(19, this%lattice%njij), source=0.0_rp)
      if (damping) allocate (dcontrib(9, nenergy))
      if (inertia) allocate (itotal(9, nenergy), ipair_b(4, nenergy, this%lattice%njij), ipair(18, this%lattice%njij), source=0.0_rp)
      if (inertia) allocate (icontrib(9, nenergy))

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
         if (do_exchange .and. moments) then
            call moment_exchange_pair(this, nodes, h, kernel, i, j, xc(:, njij_glob), xcso(:, njij_glob), xcfo(:, njij_glob), &
                                      xcparts(:, njij_glob))
         else if (do_exchange) then
            call exchange_pair(this, nodes, slot, slot, i, j, xc(:, njij_glob))
            call exchange_twoindex_pair(this, nodes, slot, slot, i, j, xcso(:, njij_glob), xcfo(:, njij_glob), xcparts(:, njij_glob))
         end if
         if (damping) then
            call damping_pair(this, slot, i, j, ief, dpair(:, njij_glob), dcontrib)
            dtotal = dtotal + dcontrib
         end if
         if (inertia) then
            call inertia_pair(this, slot, i, j, ief, step_size, ipair(:, njij_glob), ipair_b(:, :, njij_glob), icontrib)
            itotal = itotal + icontrib
         end if
      end do

      if (do_exchange) then
         call exchange_finalize(this, xc)
         call exchange_twoindex_finalize(this, xcso, xcfo, xcparts)
      end if
      if (damping) call damping_finalize(this, dpair, dtotal, ief)
      if (inertia) call inertia_finalize(this, ipair, ipair_b, itotal, ief)
   end subroutine calculate_pairwise

   !> J, D, A of pair (i, j) from the raw moments in Green-function slot 1 (see fill_intersite_moments_pair). With
   !> Delta = alpha + beta E the integral is the sum over the E^p terms (alpha_i, alpha_j; beta_i, alpha_j and alpha_i, beta_j;
   !> beta_i, beta_j) of the kernels run on the raw moments (slot 1, i side) and the p-transformed moments (slot 2, j side).
   subroutine moment_exchange_pair(this, nodes, h, kernel, i, j, xc, xcso, xcfo, xcparts)
      class(exchange) :: this
      type(energy_nodes), intent(in) :: nodes
      complex(rp), intent(in) :: h(0:, 0:)
      real(rp), intent(in) :: kernel(:)
      integer, intent(in) :: i, j
      real(rp), intent(out) :: xc(13), xcso(13), xcfo(13), xcparts(28)
      real(rp) :: t(13), tso(13), tfo(13), tparts(28)
      complex(rp) :: alpha(9, 9, 2), beta(9, 9, 2), dfix(9, 9, 2)
      integer :: p, ti, tj

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
      xc = 0.0_rp; xcso = 0.0_rp; xcfo = 0.0_rp; xcparts = 0.0_rp
      do p = 0, merge(0, 2, this%use_onsite_delta)
         call this%green%transform_intersite_moments(1, 2, h(:, p), kernel)
         do tj = 0, 1
            ti = p - tj
            if (ti < 0 .or. ti > 1) cycle
            dfix(:, :, 1) = merge(beta(:, :, 1), alpha(:, :, 1), ti == 1)
            dfix(:, :, 2) = merge(beta(:, :, 2), alpha(:, :, 2), tj == 1)
            call exchange_pair(this, nodes, 1, 2, i, j, t, dfix)
            call exchange_twoindex_pair(this, nodes, 1, 2, i, j, tso, tfo, tparts, dfix)
            xc = xc + t; xcso = xcso + tso; xcfo = xcfo + tfo; xcparts = xcparts + tparts
         end do
      end do
   end subroutine moment_exchange_pair

   !> Torque-correlation damping tensor of pair (i, j) read from Green-function slot `slot`: pair(1:9) Re and pair(10:18) Im
   !> at energy index ief, pair(19) prefactor; contrib(:, nv) is the pair's energy-resolved contribution.
   subroutine damping_pair(this, slot, i, j, ief, pair, contrib)
      class(exchange) :: this
      integer, intent(in) :: slot, i, j, ief
      real(rp), intent(out) :: pair(19), contrib(:, :)
      complex(rp), allocatable :: Aij(:, :), Aji(:, :), tjd(:, :, :)
      complex(rp) :: c(3, 3)
      integer :: k, l, m, nv, nb_ij, nenergy
      real(rp) :: spin_i

      nenergy = size(this%en%ene)
      nb_ij = size(this%green%gij, 1)
      allocate (Aij(nb_ij, nb_ij), Aji(nb_ij, nb_ij), tjd(nb_ij, nb_ij, 3))
      ! Spin magnetic moment of atom i
      spin_i = 0.0_rp
      do k = 0, this%symbolic_atom(this%lattice%iz(i))%potential%lmax
         spin_i = spin_i + this%symbolic_atom(this%lattice%iz(i))%potential%ql(1, k, 1) - &
                  this%symbolic_atom(this%lattice%iz(i))%potential%ql(1, k, 2)
      end do
      pair(19) = (-0.25_rp)*(2.0_rp/(pi*spin_i))
      do l = 1, 3
         tjd(:, :, l) = transpose(conjg(this%hamiltonian%tmat(:, :, l, this%lattice%iz(j))))
      end do
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

   subroutine damping_finalize(this, pair, total_damping, ief)
      use mpi_mod
      class(exchange) :: this
      real(rp), intent(inout) :: pair(:, :), total_damping(:, :)
      integer, intent(in) :: ief
      integer :: i, j, njij_glob, nv, nenergy, unit_e, unit_p
      real(rp) :: factor, distance_alat, alpha(3, 3), mhat(3)

      nenergy = size(this%en%ene)
#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, total_damping, size(total_damping), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, pair, size(pair), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      if (rank == 0) then
         open (newunit=unit_e, FILE='damping-energy.out', STATUS='replace', ACTION='write')
         open (newunit=unit_p, FILE='alldampings.out', STATUS='replace', ACTION='write')
         write (unit_p, '(A)') '#     i      j  alpha: xx xy xz yx yy yz zx zy zz transverse (dimensionless)'// &
            '  dist rij=ri-rj (alat)'
         write (*, *) 'ief = ', this%en%ene(ief), 'fermi = ', this%en%fermi
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
         write (unit_e, '(A)') '# E-Ef (Ry)  alpha: xx xy xz yx yy yz zx zy zz (dimensionless)'
         do nv = 1, nenergy
            write (unit_e, '(F14.9,9ES16.8)') this%en%ene(nv) - this%en%fermi, total_damping(:, nv)
         end do
         close (unit_e)
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
   function inertia_ready(this, step_size) result(ready)
      class(exchange) :: this
      real(rp), intent(out) :: step_size
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

   subroutine inertia_finalize(this, pair, pair_b, total_inertia, ief)
      use mpi_mod
      class(exchange) :: this
      real(rp), intent(inout) :: pair(:, :), pair_b(:, :, :), total_inertia(:, :)
      integer, intent(in) :: ief
      integer :: i, j, njij_glob, nv, nenergy, u_re, u_im, u_e, u_p
      real(rp) :: distance_alat

      nenergy = size(this%en%ene)
#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, total_inertia, size(total_inertia), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, pair_b, size(pair_b), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, pair, size(pair), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      if (rank == 0) then
         open (newunit=u_re, FILE='example-real.out', STATUS='replace', ACTION='write')
         open (newunit=u_im, FILE='example-imag.out', STATUS='replace', ACTION='write')
         open (newunit=u_e, FILE='inertia-energy.out', STATUS='replace', ACTION='write')
         open (newunit=u_p, FILE='allinertias.out', STATUS='replace', ACTION='write')
         write (u_re, '(A)') '# i j E Re(Bij_11) Re(d2Bij_11/dE2)'
         write (u_im, '(A)') '# i j E Im(Bij_11) Im(d2Bij_11/dE2)'
         write (u_e, '(A)') '# E-Ef total raw inertia kernel: 9 real components'
         write (u_p, '(A)') '# i j 9 real components 9 imaginary components distance'
         do njij_glob = 1, this%lattice%njij
            i = this%lattice%ijpair(njij_glob, 1)
            j = this%lattice%ijpair(njij_glob, 2)
            do nv = 1, nenergy
               write (u_re, '(2I7,3ES18.8)') i, j, this%en%ene(nv), pair_b(1, nv, njij_glob), pair_b(2, nv, njij_glob)
               write (u_im, '(2I7,3ES18.8)') i, j, this%en%ene(nv), pair_b(3, nv, njij_glob), pair_b(4, nv, njij_glob)
            end do
            distance_alat = norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            write (u_p, '(2I7,18ES18.8,F12.6)') i, j, pair(:, njij_glob), distance_alat
            write (*, '(A,I6,A,I6,A,ES14.6)') 'Raw inertia kernel pair ', i, ' and ', j, ' at Ef = ', this%en%ene(ief)
            write (*, '(3ES14.6)') pair(1:3, njij_glob)
            write (*, '(3ES14.6)') pair(4:6, njij_glob)
            write (*, '(3ES14.6)') pair(7:9, njij_glob)
         end do
         do nv = 1, nenergy
            write (u_e, '(10ES18.8)') this%en%ene(nv) - this%en%fermi, total_inertia(:, nv)
         end do
         close (u_re); close (u_im); close (u_e); close (u_p)
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
