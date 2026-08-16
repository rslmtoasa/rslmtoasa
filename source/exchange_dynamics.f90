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
      !
      class(exchange) :: this
      !
      complex(rp), allocatable, dimension(:, :, :) :: Aij, Aji ! Spectral functions (non-Hermitian)
      complex(rp), allocatable, dimension(:, :, :) :: temp1, temp2, temp3 ! Temporary matrices to store
      complex(rp), allocatable, dimension(:, :, :) :: tmati, tmatj ! Torque matrices
      !
      ! Local Variables
      integer :: i, j, k, l, m, njij, nv, lmaxi, lmaxj, ief
      real(rp) :: spin_i, orbital_i, gfac, diff, factor, distance_alat
      real(rp), allocatable, dimension(:, :) :: dtott, dtottim, total_damping
      real(rp), dimension(3, 3) :: damping_tensor, damping_tensor_im

      if (this%control%nsp .eq. 2) then ! check if spin-orbit (l.s) is enabled

         call this%hamiltonian%torque_operator_collinear() ! calculate the torque operators for all NTYPE
         allocate (total_damping(9, size(this%en%ene)))
         total_damping(:, :) = 0.0_rp

         ! Open the files for writing the results
         open (UNIT=103, FILE='damping-energy.out', STATUS='replace', ACTION='write')
         open (UNIT=104, FILE='alldampings.out', STATUS='replace', ACTION='write')

         write (104, *) '    #i     #j   #xx          #xy           #xz           '// &
            '#yx           #yy           #yz           #zx           #zy           '// &
            '#zz          #0.5*(xx + yy)     #Dist           #rij = (ri - rj)'

         do njij = 1, this%lattice%njij
            ! Obtain which are the atoms in the current pair
            i = this%lattice%ijpair(njij, 1) ! Atom number in the clust file, atom i
            j = this%lattice%ijpair(njij, 2) ! Atom number in the clust file, atom j
            ! Find the distance between the atoms
            distance_alat = norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            ! Find the lmax value for each atom
            lmaxi = this%symbolic_atom(this%lattice%iz(i))%potential%lmax ! l_max for atom i
            lmaxj = this%symbolic_atom(this%lattice%iz(j))%potential%lmax ! l_max for atom j
            ! Allocate all necessary arrays and quantities
            allocate (Aij(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
            allocate (Aji(2*(lmaxj + 1)**2, 2*(lmaxj + 1)**2, size(this%en%ene)))
            allocate (temp1(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
            allocate (temp2(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
            allocate (temp3(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
            allocate (tmati(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, 3))
            allocate (tmatj(2*(lmaxj + 1)**2, 2*(lmaxj + 1)**2, 3))
            allocate (dtott(9, size(this%en%ene)))
            allocate (dtottim(9, size(this%en%ene)))
            ! Make all necessary values equal to zero
            Aij(:, :, :) = czero; Aji(:, :, :) = czero
            temp1(:, :, :) = czero; temp2(:, :, :) = czero; temp3(:, :, :) = czero
            tmati(:, :, :) = czero; tmatj(:, :, :) = czero
            dtott(:, :) = 0.0_rp; dtottim(:, :) = 0.0_rp
            ! Obtain the torque matrices for the atoms i and j in the collinear mode
            tmati(:, :, :) = this%hamiltonian%tmat(:, :, :, this%lattice%iz(i))
            tmatj(:, :, :) = this%hamiltonian%tmat(:, :, :, this%lattice%iz(j))
            ! Calculate the spin magnetic moment of the i-th atom
            spin_i = 0.0_rp   ! per-pair accumulator; feeds 1/spin_i in `factor` below
            do k = 0, lmaxi ! azimuthal quantum number index
               spin_i = spin_i + this%symbolic_atom(this%lattice%iz(i))%potential%ql(1, k, 1) - &
                        this%symbolic_atom(this%lattice%iz(i))%potential%ql(1, k, 2)
            end do
            ! Now calculate the damping value
            do nv = 1, size(this%en%ene)
               ! Calculate the anti-Hermitian parts of the GF's Aij and Aji
               Aij(:, :, nv) = this%green%gij(:, :, nv, njij) - transpose(conjg(this%green%gji(:, :, nv, njij)))
               Aji(:, :, nv) = this%green%gji(:, :, nv, njij) - transpose(conjg(this%green%gij(:, :, nv, njij)))
               m = 0
               do k = 1, 3 ! First component index
                  do l = 1, 3 ! Second component index
                     temp1(:, :, nv) = matmul(tmati(:, :, k), Aij(:, :, nv))
                     temp2(:, :, nv) = matmul(transpose(conjg(tmatj(:, :, l))), Aji(:, :, nv))
                     temp3(:, :, nv) = matmul(temp1(:, :, nv), temp2(:, :, nv))
                     m = m + 1
                     dtott(m, nv) = rtrace(temp3(:, :, nv))
                     dtottim(m, nv) = imtrace(temp3(:, :, nv))
                  end do
               end do
            end do
            do m = 1, 9
               do nv = 1, size(this%en%ene)
                  total_damping(m, nv) = total_damping(m, nv) + dtott(m, nv)
               end do
            end do
            ! Now find the closest energy point of the Fermi level
            ief = 0; diff = 1000.0_rp
            do nv = 1, size(this%en%ene)
               if (abs(this%en%ene(nv) - this%en%fermi) .lt. diff) then
                  diff = abs(this%en%ene(nv) - this%en%fermi)
                  ief = nv
               end if
            end do
            ! Calculates the prefactor of the damping formula
            factor = (-0.25_rp)*(2.0_rp/(pi*spin_i))
            ! Finally calculates the damping term and write the results on screen
            write (*, '(A,I6,A,I6,A)') 'Damping tensor between pair', i, ' and ', j, ' is'
            write (*, '(3F14.9)') factor*dtott(:, ief)
            write (*, '(A,I6,A,I6,A)') 'Imaginary part of the tensor between pair', i, ' and ', j, ' is'
            write (*, '(3F14.9)') factor*dtottim(:, ief)
            write (*, *) '----------------------'
            write (104, '(2I7,11F14.9,3F10.6)') i, j, factor*dtott(:, ief), 0.5*factor*(dtott(1, ief) + dtott(5, ief)), &
               distance_alat, this%lattice%cr(:, i) - this%lattice%cr(:, j)
            write (*, *) 'ief = ', this%en%ene(ief), 'fermi = ', this%en%fermi

            ! Deallocate all
            deallocate (Aij)
            deallocate (Aji)
            deallocate (temp1)
            deallocate (temp2)
            deallocate (temp3)
            deallocate (tmati)
            deallocate (tmatj)
            deallocate (dtott)
            deallocate (dtottim)

         end do

         ! Write the results in output files
         write (103, *) '#Energy (E-Ef)         #xx         #xy           #xz           '// &
            '#yx           #yy           #yz           #zx           #zy           #zz'
         do nv = 1, size(this%en%ene)
            write (103, '(10F14.9)') this%en%ene(nv) - this%en%fermi, factor*total_damping(:, nv)
         end do

         deallocate (total_damping)

         ! Close the files
         close (103)
         close (104)

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
      !
      class(exchange) :: this
      !
      complex(rp), allocatable, dimension(:, :, :) :: Aij, Aji ! Spectral functions (non-Hermitian)
      complex(rp), allocatable, dimension(:, :, :) :: Bij, Bji ! Hermitian part of the Green´s functions
      complex(rp), allocatable, dimension(:, :, :) :: sBij, sBji ! Second derivative w.r.t. energy of Bij and Bji
      complex(rp), allocatable, dimension(:, :, :) :: temp1, temp2, temp3 ! Temporary matrices to store
      complex(rp), allocatable, dimension(:, :, :) :: temp4, temp5, temp6 ! Temporary matrices to store
      complex(rp), allocatable, dimension(:, :, :) :: tmati, tmatj ! Torque matrices
      real(rp), allocatable, dimension(:, :, :) :: sBij_real, sBji_real, sBij_imag, sBji_imag ! Imaginary and real parts of sBij and sBji
      !
      ! Local Variables
      integer :: i, j, k, l, m, njij, nv, lmaxi, lmaxj, ief
      real(rp) :: spin_i, orbital_i, gfac, diff, factor, distance_alat, step_size
      real(rp), allocatable, dimension(:, :) :: itott, itottim, total_inertia
      real(rp), dimension(3, 3) :: inertia_tensor, inertia_tensor_im

      if (this%control%nsp .eq. 2) then

         call this%hamiltonian%torque_operator_collinear() ! calculate the torque operators for all NTYPE

         open (UNIT=105, FILE='example-real.out', STATUS='replace', ACTION='write')
         open (UNIT=106, FILE='example-imag.out', STATUS='replace', ACTION='write')

         ! Calculate the step size of energy
         step_size = this%en%ene(2) - this%en%ene(1)
         write (*, *) 'step size=', step_size

         do njij = 1, this%lattice%njij
            ! Obtain which are the atoms in the current pair
            i = this%lattice%ijpair(njij, 1) ! Atom number in the clust file, atom i
            j = this%lattice%ijpair(njij, 2) ! Atom number in the clust file, atom j
            ! Find the distance between the atoms
            distance_alat = norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            ! Find the lmax value for each atom
            lmaxi = this%symbolic_atom(this%lattice%iz(i))%potential%lmax ! l_max for atom i
            lmaxj = this%symbolic_atom(this%lattice%iz(j))%potential%lmax ! l_max for atom j
            ! Allocate all necessary arrays and quantities
            allocate (Aij(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
            allocate (Aji(2*(lmaxj + 1)**2, 2*(lmaxj + 1)**2, size(this%en%ene)))
            allocate (Bij(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
            allocate (Bji(2*(lmaxj + 1)**2, 2*(lmaxj + 1)**2, size(this%en%ene)))
            allocate (sBij(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
            allocate (sBji(2*(lmaxj + 1)**2, 2*(lmaxj + 1)**2, size(this%en%ene)))
            allocate (sBij_real(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
            allocate (sBji_real(2*(lmaxj + 1)**2, 2*(lmaxj + 1)**2, size(this%en%ene)))
            allocate (sBij_imag(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
            allocate (sBji_imag(2*(lmaxj + 1)**2, 2*(lmaxj + 1)**2, size(this%en%ene)))
            allocate (temp1(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
            allocate (temp2(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
            allocate (temp3(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
            allocate (temp4(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
            allocate (temp5(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
            allocate (temp6(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
            allocate (tmati(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, 3))
            allocate (tmatj(2*(lmaxj + 1)**2, 2*(lmaxj + 1)**2, 3))
            allocate (itott(9, size(this%en%ene)))
            allocate (itottim(9, size(this%en%ene)))
            ! Make all necessary values equal to zero
            Aij(:, :, :) = czero; Aji(:, :, :) = czero
            Bij(:, :, :) = czero; Bji(:, :, :) = czero
            sBij(:, :, :) = czero; sBji(:, :, :) = czero
            sBij_real(:, :, :) = 0.0_rp; sBji_real(:, :, :) = 0.0_rp
            sBij_imag(:, :, :) = 0.0_rp; sBji_imag(:, :, :) = 0.0_rp
            temp1(:, :, :) = czero; temp2(:, :, :) = czero; temp3(:, :, :) = czero
            temp4(:, :, :) = czero; temp5(:, :, :) = czero; temp6(:, :, :) = czero
            tmati(:, :, :) = czero; tmatj(:, :, :) = czero
            itott(:, :) = 0.0_rp; itottim(:, :) = 0.0_rp
            ! Obtain the torque matrices for the atoms i and j in the collinear mode
            tmati(:, :, :) = this%hamiltonian%tmat(:, :, :, this%lattice%iz(i))
            tmatj(:, :, :) = this%hamiltonian%tmat(:, :, :, this%lattice%iz(j))
            ! Obtain the magnetic moment (spin and orbital) of the i-th atom
            ! TO-DO: include orbital moment
            spin_i = 0.0_rp; orbital_i = 0.0_rp
            do k = 0, lmaxi ! azimuthal quantum number index
               spin_i = spin_i + this%symbolic_atom(this%lattice%iz(i))%potential%ql(1, k, 1) - &
                        this%symbolic_atom(this%lattice%iz(i))%potential%ql(1, k, 2)
               !do l = -k, k ! magnetic quantum number index
               !  orbital_i = orbital_i + l*(this%symbolic_atom(this%lattice%iz(i))%potential%ql(1, k, 1) + &
               !                             this%symbolic_atom(this%lattice%iz(i))%potential%ql(1, k, 2))
               !end do
            end do
            do nv = 1, size(this%en%ene)
               ! Calculate the anti-Hermitian parts of the GF´s Aij and Aji
               Aij(:, :, nv) = this%green%gij(:, :, nv, njij) - transpose(conjg(this%green%gji(:, :, nv, njij)))
               Aji(:, :, nv) = this%green%gji(:, :, nv, njij) - transpose(conjg(this%green%gij(:, :, nv, njij)))
               ! Calculate the Hermitian parts of the GF´s Bij and Bji
               Bij(:, :, nv) = this%green%gij(:, :, nv, njij) + transpose(conjg(this%green%gji(:, :, nv, njij)))
               Bji(:, :, nv) = this%green%gji(:, :, nv, njij) + transpose(conjg(this%green%gij(:, :, nv, njij)))
            end do
            ! Calculate the second derivative with respect to energy of Bij and Bji
            do k = 1, size(Bij, 1) ! Size of first index (Bji is the same)
               do l = 1, size(Bij, 2) ! Size of second index (Bji is the same)
                  sBij_real(k, l, :) = second_derivative(real(Bij(k, l, :)), step_size)
                  sBji_real(k, l, :) = second_derivative(real(Bji(k, l, :)), step_size)
                  sBij_imag(k, l, :) = second_derivative(aimag(Bij(k, l, :)), step_size)
                  sBji_imag(k, l, :) = second_derivative(aimag(Bji(k, l, :)), step_size)
               end do
            end do
            ! Construct the matrices sBij and sBji
            do k = 1, size(Bij, 1) ! Size of first index (Bji is the same)
               do l = 1, size(Bij, 2) ! Size of second index (Bji is the same)
                  do nv = 1, size(this%en%ene)
                     sBij(k, l, nv) = cmplx(sBij_real(k, l, nv), sBij_imag(k, l, nv))
                     sBji(k, l, nv) = cmplx(sBji_real(k, l, nv), sBji_imag(k, l, nv))
                  end do
               end do
            end do
            do nv = 1, size(this%en%ene)
               write (105, *) this%en%ene(nv), real(Bij(1, 1, nv)), real(sBij(1, 1, nv))
               write (106, *) this%en%ene(nv), aimag(Bij(1, 1, nv)), aimag(sBij(1, 1, nv))
            end do
            ! Now calculate the moment of inertia value
            m = 0
            do k = 1, 3 ! First component index
               do l = 1, 3 ! Second component index
                  temp1(:, :, nv) = matmul(tmati(:, :, k), Aij(:, :, nv))
                  temp2(:, :, nv) = matmul(transpose(conjg(tmatj(:, :, l))), sBji(:, :, nv))
                  temp3(:, :, nv) = matmul(tmati(:, :, k), sBij(:, :, nv))
                  temp4(:, :, nv) = matmul(transpose(conjg(tmatj(:, :, l))), Aji(:, :, nv))
                  temp5(:, :, nv) = matmul(temp1(:, :, nv), temp2(:, :, nv))
                  temp6(:, :, nv) = matmul(temp3(:, :, nv), temp4(:, :, nv))
                  m = m + 1
                  itott(m, nv) = rtrace(temp5(:, :, nv) + temp6(:, :, nv))
                  itottim(m, nv) = imtrace(temp5(:, :, nv) + temp6(:, :, nv))
               end do
            end do
            ! Deallocate all arrays
            deallocate (Aij)
            deallocate (Aji)
            deallocate (Bij)
            deallocate (Bji)
            deallocate (sBij)
            deallocate (sBji)
            deallocate (sBij_real)
            deallocate (sBji_real)
            deallocate (sBij_imag)
            deallocate (sBji_imag)
            deallocate (temp1)
            deallocate (temp2)
            deallocate (temp3)
            deallocate (temp4)
            deallocate (temp5)
            deallocate (temp6)
            deallocate (tmati)
            deallocate (tmatj)
            deallocate (itott)
            deallocate (itottim)

         end do

         close (105)
         close (106)
      end if

   end subroutine calculate_moment_of_inertia


end submodule exchange_dynamics
