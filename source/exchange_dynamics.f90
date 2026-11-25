submodule (exchange_mod) exchange_dynamics
   implicit none

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculation of the Gilbert damping parameters in real-space
   !> based on the LMTO solution (see, e.g., PRM 2, 013801 (2018)).
   !> The model is the so-called torque-correlation
   !> The following calculation works better with hoh.
   !> Implemented by Ivan Miranda on 14.11.2023
   !> Revised by Ivan Miranda on 29.04.2025
   !---------------------------------------------------------------------------
   module subroutine calculate_gilbert_damping(this)
      use mpi_mod
      !
      class(exchange) :: this
      !
      complex(rp), allocatable :: Aij(:, :), Aji(:, :), temp1(:, :), temp2(:, :), temp3(:, :)
      integer :: i, j, k, l, m, njij, njij_glob, nv, ief, nb_ij, nenergy, unit_e, unit_p
      real(rp) :: spin_i, factor, distance_alat
      ! pair(1:9) Re tensor at E_F, pair(10:18) Im tensor at E_F, pair(19) prefactor
      real(rp), allocatable :: pair(:, :), total_damping(:, :)

      if (this%control%nsp /= 2) then
         call g_logger%warning('calculate_gilbert_damping: needs nsp = 2 (collinear with l.s); skipped.', __FILE__, __LINE__)
         return
      end if

      call this%hamiltonian%torque_operator_collinear()
      nenergy = size(this%en%ene)
      nb_ij = size(this%green%gij, 1)
      allocate (total_damping(9, nenergy), pair(19, this%lattice%njij), source=0.0_rp)
      allocate (Aij(nb_ij, nb_ij), Aji(nb_ij, nb_ij), temp1(nb_ij, nb_ij), temp2(nb_ij, nb_ij), temp3(nb_ij, nb_ij))
      ief = fermi_index(this%en%ene, this%en%fermi)

      ! Each rank holds the Green functions of its pairs start_atom..end_atom (local index g2l_map)
      do njij_glob = start_atom, end_atom
         njij = g2l_map(njij_glob)
         i = this%lattice%ijpair(njij_glob, 1)
         j = this%lattice%ijpair(njij_glob, 2)
         ! Spin magnetic moment of atom i
         spin_i = 0.0_rp
         do k = 0, this%symbolic_atom(this%lattice%iz(i))%potential%lmax
            spin_i = spin_i + this%symbolic_atom(this%lattice%iz(i))%potential%ql(1, k, 1) - &
                     this%symbolic_atom(this%lattice%iz(i))%potential%ql(1, k, 2)
         end do
         do nv = 1, nenergy
            ! Anti-Hermitian parts of the intersite Green functions
            Aij = this%green%gij(:, :, nv, njij) - transpose(conjg(this%green%gji(:, :, nv, njij)))
            Aji = this%green%gji(:, :, nv, njij) - transpose(conjg(this%green%gij(:, :, nv, njij)))
            m = 0
            do k = 1, 3
               do l = 1, 3
                  temp1 = matmul(this%hamiltonian%tmat(:, :, k, this%lattice%iz(i)), Aij)
                  temp2 = matmul(transpose(conjg(this%hamiltonian%tmat(:, :, l, this%lattice%iz(j)))), Aji)
                  temp3 = matmul(temp1, temp2)
                  m = m + 1
                  total_damping(m, nv) = total_damping(m, nv) + rtrace(temp3)
                  if (nv == ief) then
                     pair(m, njij_glob) = rtrace(temp3)
                     pair(9 + m, njij_glob) = imtrace(temp3)
                  end if
               end do
            end do
         end do
         pair(19, njij_glob) = (-0.25_rp)*(2.0_rp/(pi*spin_i))
      end do

#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, total_damping, size(total_damping), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, pair, size(pair), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      if (rank == 0) then
         open (newunit=unit_e, FILE='damping-energy.out', STATUS='replace', ACTION='write')
         open (newunit=unit_p, FILE='alldampings.out', STATUS='replace', ACTION='write')
         write (unit_p, '(A)') '#     i      j  alpha: xx xy xz yx yy yz zx zy zz 0.5*(xx+yy) (dimensionless)'// &
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
            write (unit_p, '(2I7,10ES16.8,4F11.6)') i, j, factor*pair(1:9, njij_glob), &
               0.5*factor*(pair(1, njij_glob) + pair(5, njij_glob)), &
               distance_alat, this%lattice%cr(:, i) - this%lattice%cr(:, j)
         end do

         ! Energy-resolved total, scaled by the factor of the last pair (see tests/KNOWN_ISSUES.md)
         write (unit_e, '(A)') '# E-Ef (Ry)  alpha: xx xy xz yx yy yz zx zy zz (dimensionless)'
         do nv = 1, nenergy
            write (unit_e, '(F14.9,9ES16.8)') this%en%ene(nv) - this%en%fermi, factor*total_damping(:, nv)
         end do
         close (unit_e)
         close (unit_p)
      end if

   end subroutine calculate_gilbert_damping

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
      complex(rp), allocatable :: Aij(:, :, :), Aji(:, :, :), Bij(:, :, :), Bji(:, :, :)
      complex(rp), allocatable :: sBij(:, :, :), sBji(:, :, :)
      complex(rp), allocatable :: temp1(:, :), temp2(:, :), temp3(:, :)
      complex(rp), allocatable :: temp4(:, :), temp5(:, :), temp6(:, :)
      ! pair(1:9) Re kernel at E_F, pair(10:18) Im kernel at E_F
      real(rp), allocatable :: total_inertia(:, :), pair_b(:, :, :), pair(:, :)
      integer :: i, j, k, l, m, njij, njij_glob, nv, ief, nenergy, nb_ij, u_re, u_im, u_e, u_p
      real(rp) :: distance_alat, step_size

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
      nb_ij = size(this%green%gij, 1)
      allocate (total_inertia(9, nenergy), pair_b(4, nenergy, this%lattice%njij), pair(18, this%lattice%njij), source=0.0_rp)
      allocate (Aij(nb_ij, nb_ij, nenergy), Aji(nb_ij, nb_ij, nenergy), Bij(nb_ij, nb_ij, nenergy), Bji(nb_ij, nb_ij, nenergy))
      allocate (sBij(nb_ij, nb_ij, nenergy), sBji(nb_ij, nb_ij, nenergy))
      allocate (temp1(nb_ij, nb_ij), temp2(nb_ij, nb_ij), temp3(nb_ij, nb_ij))
      allocate (temp4(nb_ij, nb_ij), temp5(nb_ij, nb_ij), temp6(nb_ij, nb_ij))
      ief = fermi_index(this%en%ene, this%en%fermi)

      ! Each rank holds the Green functions of its pairs start_atom..end_atom (local index g2l_map)
      do njij_glob = start_atom, end_atom
         njij = g2l_map(njij_glob)
         i = this%lattice%ijpair(njij_glob, 1)
         j = this%lattice%ijpair(njij_glob, 2)

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

         pair_b(1, :, njij_glob) = real(Bij(1, 1, :)); pair_b(2, :, njij_glob) = real(sBij(1, 1, :))
         pair_b(3, :, njij_glob) = aimag(Bij(1, 1, :)); pair_b(4, :, njij_glob) = aimag(sBij(1, 1, :))
         do nv = 1, nenergy
            m = 0
            do k = 1, 3
               do l = 1, 3
                  temp1 = matmul(this%hamiltonian%tmat(:, :, k, this%lattice%iz(i)), Aij(:, :, nv))
                  temp2 = matmul(transpose(conjg(this%hamiltonian%tmat(:, :, l, this%lattice%iz(j)))), sBji(:, :, nv))
                  temp3 = matmul(this%hamiltonian%tmat(:, :, k, this%lattice%iz(i)), sBij(:, :, nv))
                  temp4 = matmul(transpose(conjg(this%hamiltonian%tmat(:, :, l, this%lattice%iz(j)))), Aji(:, :, nv))
                  temp5 = matmul(temp1, temp2)
                  temp6 = matmul(temp3, temp4)
                  m = m + 1
                  total_inertia(m, nv) = total_inertia(m, nv) + rtrace(temp5 + temp6)
                  if (nv == ief) then
                     pair(m, njij_glob) = rtrace(temp5 + temp6)
                     pair(9 + m, njij_glob) = imtrace(temp5 + temp6)
                  end if
               end do
            end do
         end do
      end do

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
   end subroutine calculate_moment_of_inertia

   !> Index of the energy-grid point closest to the Fermi level.
   pure integer function fermi_index(ene, fermi) result(ief)
      real(rp), intent(in) :: ene(:), fermi

      ief = minloc(abs(ene - fermi), dim=1)
   end function fermi_index

end submodule exchange_dynamics
