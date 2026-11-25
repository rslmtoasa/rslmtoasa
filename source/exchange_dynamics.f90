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
   module subroutine calculate_gilbert_damping(this)
      use mpi_mod
      !
      class(exchange) :: this
      !
      integer :: i, j, njij, njij_glob, ief, nenergy
      ! pair(1:9) Re tensor at E_F, pair(10:18) Im tensor at E_F, pair(19) prefactor
      real(rp), allocatable :: pair(:, :), total_damping(:, :), contrib(:, :)

      if (this%control%nsp /= 2) then
         call g_logger%warning('calculate_gilbert_damping: needs nsp = 2 (collinear with l.s); skipped.', __FILE__, __LINE__)
         return
      end if

      call this%hamiltonian%torque_operator_collinear()
      nenergy = size(this%en%ene)
      allocate (total_damping(9, nenergy), pair(19, this%lattice%njij), source=0.0_rp)
      allocate (contrib(9, nenergy))
      ief = fermi_index(this%en%ene, this%en%fermi)

      ! Each rank holds the Green functions of its pairs start_atom..end_atom (local index g2l_map)
      do njij_glob = start_atom, end_atom
         njij = g2l_map(njij_glob)
         i = this%lattice%ijpair(njij_glob, 1)
         j = this%lattice%ijpair(njij_glob, 2)
         call damping_pair(this, njij, i, j, ief, pair(:, njij_glob), contrib)
         total_damping = total_damping + contrib
      end do

      call damping_finalize(this, pair, total_damping, ief)
   end subroutine calculate_gilbert_damping

   !> Torque-correlation damping tensor of pair (i, j) at local index njij: pair(1:9) Re and pair(10:18) Im
   !> at energy index ief, pair(19) prefactor; contrib(:, nv) is the pair's energy-resolved contribution.
   subroutine damping_pair(this, njij, i, j, ief, pair, contrib)
      class(exchange) :: this
      integer, intent(in) :: njij, i, j, ief
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
            Aij = this%green%gij(:, :, nv, njij) - transpose(conjg(this%green%gji(:, :, nv, njij)))
            Aji = this%green%gji(:, :, nv, njij) - transpose(conjg(this%green%gij(:, :, nv, njij)))
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
   module subroutine calculate_moment_of_inertia(this)
      use mpi_mod
      !
      class(exchange) :: this
      !
      ! This is the experimental torque-correlation kernel described in
      ! Sci. Rep. 7, 931 (2017).  It is deliberately reported as the raw
      ! production kernel: unlike Gilbert damping above, this routine has no
      ! independently established normalization or production relation.
      ! pair(1:9) Re kernel at E_F, pair(10:18) Im kernel at E_F
      real(rp), allocatable :: total_inertia(:, :), pair_b(:, :, :), pair(:, :), contrib(:, :)
      integer :: i, j, njij, njij_glob, ief, nenergy
      real(rp) :: step_size

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
      allocate (total_inertia(9, nenergy), pair_b(4, nenergy, this%lattice%njij), pair(18, this%lattice%njij), source=0.0_rp)
      allocate (contrib(9, nenergy))
      ief = fermi_index(this%en%ene, this%en%fermi)

      ! Each rank holds the Green functions of its pairs start_atom..end_atom (local index g2l_map)
      do njij_glob = start_atom, end_atom
         njij = g2l_map(njij_glob)
         i = this%lattice%ijpair(njij_glob, 1)
         j = this%lattice%ijpair(njij_glob, 2)
         call inertia_pair(this, njij, i, j, ief, step_size, pair(:, njij_glob), pair_b(:, :, njij_glob), contrib)
         total_inertia = total_inertia + contrib
      end do

      call inertia_finalize(this, pair, pair_b, total_inertia, ief)
   end subroutine calculate_moment_of_inertia

   !> Raw inertia kernel of pair (i, j) at local index njij: pair(1:9) Re and pair(10:18) Im at energy index ief,
   !> pair_b the Bij(1, 1) trace and its second derivative, contrib(:, nv) the energy-resolved real part.
   subroutine inertia_pair(this, njij, i, j, ief, step_size, pair, pair_b, contrib)
      class(exchange) :: this
      integer, intent(in) :: njij, i, j, ief
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
            Aij(:, :, nv) = this%green%gij(:, :, nv, njij) - transpose(conjg(this%green%gji(:, :, nv, njij)))
            Aji(:, :, nv) = this%green%gji(:, :, nv, njij) - transpose(conjg(this%green%gij(:, :, nv, njij)))
            Bij(:, :, nv) = this%green%gij(:, :, nv, njij) + transpose(conjg(this%green%gji(:, :, nv, njij)))
            Bji(:, :, nv) = this%green%gji(:, :, nv, njij) + transpose(conjg(this%green%gij(:, :, nv, njij)))
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
