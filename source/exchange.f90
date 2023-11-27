!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Exchange
!
!> @author
!> Angela Klautau
!> Ramon Cardias
!> Lucas P. Campagna
!> S. Frota-Pessôa
!> Pascoal R. Peduto
!> Anders Bergman
!> S. B. Legoas
!> H. M. Petrilli
!> Ivan P. Miranda
!
! DESCRIPTION:
!> Module to handle calculation of exchange interactions
!------------------------------------------------------------------------------

module exchange_mod
   use bands_mod
   use energy_mod
   use green_mod
   use lattice_mod
   use control_mod
   use symbolic_atom_mod
   use density_of_states_mod
   use recursion_mod
   use hamiltonian_mod
   use precision_mod, only: rp
   use math_mod
   use logger_mod, only: g_logger
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   use string_mod
#ifdef USE_MPI
   use mpi
#endif
   implicit none

   private
   !> Module's main structure
   type, public :: exchange
      !> Green
      class(green), pointer :: green
      !> Lattice
      class(lattice), pointer :: lattice
      !> Symbolic atom
      class(symbolic_atom), dimension(:), pointer :: symbolic_atom
      !> Density of states
      class(dos), pointer :: dos
      !> Energy
      class(energy), pointer :: en
      !> Recursion
      class(recursion), pointer :: recursion
      !> Bands
      class(bands), pointer :: bands
      !> Hamiltonian
      class(hamiltonian), pointer :: hamiltonian
      !> Control
      class(control), pointer :: control

      !> General variables
      !> Heisenberg exchange
      real(rp) :: jij
      !> Heisenberg exchange tensor (obtained with auxiliary GF's)
      real(rp), dimension(9) :: jij_aux
      !> Total Heisenberg exchange J0 = sum_j J0j (obtained with auxiliary GF's)
      real(rp) :: jij00_aux
      !> Spin-lattice couplings tensor
      real(rp), dimension(9) :: jijk
      !> Dzyaloshinskii-Moriya interaction
      real(rp), dimension(3) :: dmi
      !> Anisotropy
      real(rp), dimension(3, 3) :: aij
   contains
      procedure :: calculate_jij
      procedure :: calculate_dij
      procedure :: calculate_aij
      procedure :: calculate_exchange
      procedure :: calculate_jij_auxgreen
      procedure :: calculate_jijk
      procedure :: calculate_exchange_gauss_legendre
      procedure :: calculate_exchange_rs2pao
      procedure :: calculate_gilbert_damping
      procedure :: restore_to_default
      final destructor
   end type

   interface exchange
      procedure :: constructor
   end interface exchange

contains

   function constructor(bands_obj) result(obj)
      type(exchange) :: obj
      type(bands), target, intent(in) :: bands_obj

      obj%bands => bands_obj
      obj%green => bands_obj%green
      obj%lattice => bands_obj%lattice
      obj%symbolic_atom => bands_obj%symbolic_atom
      obj%en => bands_obj%en
      obj%control => bands_obj%lattice%control
      obj%recursion => bands_obj%recursion
      obj%hamiltonian => bands_obj%recursion%hamiltonian

      call obj%restore_to_default()
   end function constructor

   subroutine destructor(this)
      type(exchange) :: this
   end subroutine destructor

   subroutine restore_to_default(this)
      class(exchange) :: this

      this%jij = 0.0D0
      this%jij_aux(:) = 0.0_RP
      this%jij00_aux = 0.0_RP
      this%jijk(:) = 0.0_RP
      this%dmi(:) = 0.0D0
      this%aij(:, :) = 0.0D0
   end subroutine restore_to_default

   subroutine calculate_jij(this)
      class(exchange) :: this
      complex(rp), dimension(9, 9) :: dmat1, dmat2, jmat
      real(rp), dimension(this%en%channels_ldos + 10) :: jtot
      integer :: nv, i, j, njij

      jtot = 0.0D0
      do njij = 1, this%lattice%njij
         do nv = 1, this%en%channels_ldos + 10
            i = this%lattice%ijpair(njij, 1) ! Atom number in the clust file, atom i
            j = this%lattice%ijpair(njij, 2) ! Atom number in the clust file, atom j
            call this%symbolic_atom(this%lattice%iz(i))%d_matrix(dmat1, this%en%ene(nv))
            call this%symbolic_atom(this%lattice%iz(j))%d_matrix(dmat2, this%en%ene(nv))
            !J
            jmat = dGdG_Jnc(this, dmat1, dmat2, njij, nv)
            jtot(nv) = imtrace9(jmat)
         end do
         call simpson_f(this%jij, this%en%ene, this%en%fermi, this%en%nv1, jtot, .true., .false., 0.0D0)
         write (*, *) 'Jij between pair', i, 'and ', j, 'is ', this%jij * 1.0D3 / 4.0D0 / pi
      end do
   end subroutine calculate_jij

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculation of the Jij following the auxiliary GF formalism
   !> in LMTO (see, e.g., Eq.3 of PRB 64, 174402 (2001)).
   !> The following calculation works better with hoh, as
   !> it supposes an orthogonal representation.
   !> Implemented by Ivan Miranda on 11.10.2023
   !---------------------------------------------------------------------------
   subroutine calculate_jij_auxgreen(this)
      !
      class(exchange) :: this
      !
      complex(rp), allocatable, dimension(:, :, :) :: pmatrix_i, pmatrix_j ! P matrices (potential functions) for atoms i and j
      complex(rp), allocatable, dimension(:, :, :) :: deltap_i, deltap_j ! P_up - P_dw (half of the size, as it does not consider the spin index)
      complex(rp), allocatable, dimension(:, :, :) :: aux_gij, aux_gji ! Auxiliary Green's functions
      complex(rp), allocatable, dimension(:, :, :, :) :: int_all ! Result for the components of the Jij tensor (when i != j)
      complex(rp), allocatable, dimension(:, :, :) :: int_j00 ! Result for J00 (when i = j)
      complex(rp), allocatable, dimension(:, :, :) :: temp1, temp2, temp3, temp4 ! Temporary matrices to store
      real(rp), allocatable, dimension(:, :) :: jtot_aux ! Jij total from the auxiliary GF's formalism (when i != j)
      real(rp), allocatable, dimension(:) :: jtot_00 ! J00 total from the auxiliary GF's formalism (when i = j)
      real(rp), dimension(4, 9) :: angles ! For the definition of xx, xy, xz, yx, yy, ... components
      integer :: nv, i, j, k, njij, lmaxi, lmaxj ! Internal (local) variables

      ! Definition of each of the angles
      do k = 1, 9
         if (k .eq. 1) then ! xx
            angles(1, k) = 0.5_RP * pi ! theta
            angles(2, k) = 0.5_RP * pi ! theta'
            angles(3, k) = 0.0_RP ! phi
            angles(4, k) = 0.0_RP ! phi'
         else if (k .eq. 2) then ! xy
            angles(1, k) = 0.5_RP * pi ! theta
            angles(2, k) = 0.5_RP * pi ! theta'
            angles(3, k) = 0.0_RP ! phi
            angles(4, k) = 0.5_RP * pi ! phi'
         else if (k .eq. 3) then ! xz
            angles(1, k) = 0.5_RP * pi ! theta
            angles(2, k) = 0.0_RP ! theta'
            angles(3, k) = 0.0_RP ! phi
            angles(4, k) = 0.0_RP ! phi'
         else if (k .eq. 4) then ! yx
            angles(1, k) = 0.5_RP * pi ! theta
            angles(2, k) = 0.5_RP * pi ! theta'
            angles(3, k) = 0.5_RP * pi ! phi
            angles(4, k) = 0.0_RP ! phi'
         else if (k .eq. 5) then ! yy
            angles(1, k) = 0.5_RP * pi ! theta
            angles(2, k) = 0.5_RP * pi ! theta'
            angles(3, k) = 0.5_RP * pi ! phi
            angles(4, k) = 0.5_RP * pi ! phi'
         else if (k .eq. 6) then ! yz
            angles(1, k) = 0.5_RP * pi ! theta
            angles(2, k) = 0.0_RP ! theta'
            angles(3, k) = 0.5_RP * pi ! phi
            angles(4, k) = 0.0_RP ! phi'
         else if (k .eq. 7) then ! zx
            angles(1, k) = 0.0_RP ! theta
            angles(2, k) = 0.5_RP * pi ! theta'
            angles(3, k) = 0.0_RP ! phi
            angles(4, k) = 0.0_RP ! phi'
         else if (k .eq. 8) then ! zy
            angles(1, k) = 0.0_RP ! theta
            angles(2, k) = 0.5_RP * pi ! theta'
            angles(3, k) = 0.0_RP ! phi
            angles(4, k) = 0.5_RP * pi ! phi'
         else ! zz
            angles(1, k) = 0.0_RP ! theta
            angles(2, k) = 0.0_RP ! theta'
            angles(3, k) = 0.0_RP ! phi
            angles(4, k) = 0.0_RP ! phi'
         end if
      end do

      ! Loop over all pairs informed by the user
      do njij = 1, this%lattice%njij
         ! Obtain which are the atoms in the current pair
         i = this%lattice%ijpair(njij, 1) ! Atom number in the clust file, atom i
         j = this%lattice%ijpair(njij, 2) ! Atom number in the clust file, atom j
         ! Now obtain the lmax value for each atom
         lmaxi = this%symbolic_atom(this%lattice%iz(i))%potential%lmax ! l_max for atom i
         lmaxj = this%symbolic_atom(this%lattice%iz(j))%potential%lmax ! l_max for atom j
         ! Allocate all necessary vectors and matrices
         allocate (pmatrix_i(2 * (lmaxi + 1)**2, 2 * (lmaxi + 1)**2, size(this%en%ene)))
         allocate (pmatrix_j(2 * (lmaxj + 1)**2, 2 * (lmaxj + 1)**2, size(this%en%ene)))
         allocate (aux_gij(2 * (lmaxi + 1)**2, 2 * (lmaxi + 1)**2, size(this%en%ene)))
         allocate (aux_gji(2 * (lmaxj + 1)**2, 2 * (lmaxj + 1)**2, size(this%en%ene)))
         allocate (deltap_i((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (deltap_j((lmaxj + 1)**2, (lmaxj + 1)**2, size(this%en%ene)))
         allocate (temp1((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (temp2((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (temp3((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (temp4((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (int_j00((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (int_all((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene), 9))
         allocate (jtot_00(size(this%en%ene)))
         allocate (jtot_aux(size(this%en%ene), 9))
         ! Set all matrices and arrays to zero first
         jtot_aux(:, :) = 0.0_RP
         jtot_00(:) = 0.0_RP
         pmatrix_i = czero; pmatrix_j = czero
         aux_gij = czero; aux_gji = czero
         deltap_i = czero; deltap_j = czero
         int_all = czero; int_j00 = czero
         ! Now calculate the potential functions P for atoms i and j
         call this%symbolic_atom(this%lattice%iz(i))%p_matrix(pmatrix_i, lmaxi, this%en%ene)
         call this%symbolic_atom(this%lattice%iz(j))%p_matrix(pmatrix_j, lmaxj, this%en%ene)
         ! Now obtain the auxiliar Green's function for atoms i and j (and vice-versa)
         call this%green%auxiliary_gij(this%green%gij(:, :, :, njij), aux_gij, i, j)
         call this%green%auxiliary_gij(this%green%gji(:, :, :, njij), aux_gji, j, i)
         ! Now transform the pmatrices in DeltaP = P_up - P_dw
         deltap_i(:, :, :) = pmatrix_i(1:(lmaxi + 1)**2, 1:(lmaxi + 1)**2, :) - pmatrix_i((lmaxi + 1)**2 + 1:2 * (lmaxi + 1)**2, (lmaxi + 1)**2 + 1:2 * (lmaxi + 1)**2, :)
         deltap_j(:, :, :) = pmatrix_j(1:(lmaxj + 1)**2, 1:(lmaxj + 1)**2, :) - pmatrix_j((lmaxj + 1)**2 + 1:2 * (lmaxj + 1)**2, (lmaxj + 1)**2 + 1:2 * (lmaxj + 1)**2, :)
         ! Now do the calculation of all the components of the jtot tensor
         temp1 = czero; temp2 = czero
         temp3 = czero; temp4 = czero
         do nv = 1, size(this%en%ene)
            ! Calculates Jij tensor
            if (i .ne. j) then
               temp1(:, :, nv) = matmul(deltap_i(:, :, nv), aux_gij(1:(lmaxi + 1)**2, 1:(lmaxi + 1)**2, nv)) ! (P_up(i) - P_dw(i))*gij_upup
               temp2(:, :, nv) = matmul(deltap_j(:, :, nv), aux_gji((lmaxj + 1)**2 + 1:2 * (lmaxj + 1)**2, (lmaxj + 1)**2 + 1:2 * (lmaxj + 1)**2, nv)) ! (P_up(j) - P_dw(j))*gji_dwdw
               temp3(:, :, nv) = matmul(deltap_i(:, :, nv), aux_gij((lmaxi + 1)**2 + 1:2 * (lmaxi + 1)**2, (lmaxi + 1)**2 + 1:2 * (lmaxi + 1)**2, nv)) ! (P_up(i) - P_dw(i))*gij_dwdw
               temp4(:, :, nv) = matmul(deltap_j(:, :, nv), aux_gji(1:(lmaxj + 1)**2, 1:(lmaxj + 1)**2, nv)) ! (P_up(j) - P_dw(j))*gji_upup
               do k = 1, 9 ! run over the tensor elements
                  int_all(:, :, nv, k) = scalar_multiply(matmul(temp1(:, :, nv), temp4(:, :, nv)), cone * cos(angles(1, k)) * cos(angles(2, k))) + &
                                         scalar_multiply(matmul(temp3(:, :, nv), temp4(:, :, nv)), cone * sin(angles(1, k)) * sin(angles(2, k)) * exp(i_unit * (angles(4, k) - angles(3, k)))) + &
                                         scalar_multiply(matmul(temp1(:, :, nv), temp2(:, :, nv)), cone * sin(angles(1, k)) * sin(angles(2, k)) * exp(i_unit * (angles(3, k) - angles(4, k)))) + &
                                         scalar_multiply(matmul(temp3(:, :, nv), temp2(:, :, nv)), cone * cos(angles(1, k)) * cos(angles(2, k)))
                  jtot_aux(nv, k) = imtrace(int_all(:, :, nv, k)) * (0.5_RP)
               end do
               ! Calculates J00 otherwise
            else ! i = j
               temp1(:, :, nv) = matmul(deltap_i(:, :, nv), aux_gij(1:(lmaxi + 1)**2, 1:(lmaxi + 1)**2, nv)) ! (P_up(i) - P_dw(i))*gii_upup
               temp2(:, :, nv) = matmul(deltap_j(:, :, nv), aux_gji((lmaxj + 1)**2 + 1:2 * (lmaxj + 1)**2, (lmaxj + 1)**2 + 1:2 * (lmaxj + 1)**2, nv)) ! (P_up(i) - P_dw(i))*gii_dwdw
               ! (P_up(i) - P_dw(i))*(gii_upup - gii_dwdw)
               temp3(:, :, nv) = matmul(deltap_i(:, :, nv), aux_gij(1:(lmaxi + 1)**2, 1:(lmaxi + 1)**2, nv) - aux_gji((lmaxj + 1)**2 + 1:2 * (lmaxj + 1)**2, (lmaxj + 1)**2 + 1:2 * (lmaxj + 1)**2, nv))
               int_j00(:, :, nv) = matmul(temp1(:, :, nv), temp2(:, :, nv)) + temp3(:, :, nv)
               jtot_00(nv) = imtrace(int_j00(:, :, nv)) * (-1.0_RP)
            end if
         end do
         ! Integrate jtot components with respect to energy
         if (i .ne. j) then
            do k = 1, 9
               call simpson_f(this%jij_aux(k), this%en%ene, this%en%fermi, this%en%nv1, jtot_aux(:, k), .true., .false., 0.0D0)
            end do
         else ! i = j
            call simpson_f(this%jij00_aux, this%en%ene, this%en%fermi, this%en%nv1, jtot_00, .true., .false., 0.0D0)
         end if
         ! Write the result
         if (i .ne. j) then
            write (*, *) 'Jij_aux tensor between pair', i, 'and ', j, 'is '
            write (*, '(3F14.9)') this%jij_aux * 1.0D3 / 4.0D0 / pi
            write (*, *) 'Dij_zz_aux between pair', i, 'and', j, 'is ', 0.5_RP * (this%jij_aux(2) - this%jij_aux(4)) * (1.0D3 / 4.0D0 / pi)
         else ! i = j
            write (*, *) 'J0_aux is', this%jij00_aux * 1.0D3 / 4.0D0 / pi
         end if
         ! Deallocate all vectors and matrices previously allocated
         deallocate (pmatrix_i)
         deallocate (pmatrix_j)
         deallocate (aux_gij)
         deallocate (aux_gji)
         deallocate (deltap_i)
         deallocate (deltap_j)
         deallocate (temp1)
         deallocate (temp2)
         deallocate (temp3)
         deallocate (temp4)
         deallocate (int_all)
         deallocate (int_j00)
         deallocate (jtot_aux)
         deallocate (jtot_00)
      end do

   end subroutine calculate_jij_auxgreen

   !> Calculation of the spin-lattice interaction
   subroutine calculate_jijk(this)
      !
      class(exchange) :: this
      !
      complex(rp), allocatable, dimension(:, :, :) :: aux_gij, aux_gji, aux_gik, aux_gki, aux_gjk, aux_gkj ! Auxiliary GF's (in the orthogonal representaton)
      complex(rp), allocatable, dimension(:, :, :) :: aux0_gij, aux0_gji, aux0_gik, aux0_gki, aux0_gjk, aux0_gkj ! Same as above, but in the canonical representation
      complex(rp), allocatable, dimension(:, :, :) :: pmatrix_i, pmatrix_j, pmatrix_k ! P matrices (potential function) for atoms i, j, k (orthogonal representation)
      complex(rp), allocatable, dimension(:, :, :) :: pmatrix0_i, pmatrix0_j, pmatrix0_k ! P matrix (potential function) of atom k in the canonical representation
      complex(rp), allocatable, dimension(:, :, :) :: deltap_i, deltap_j ! P_up - P_dw matrices (half of the size, as it does not consider the spin index)
      complex(rp), allocatable, dimension(:, :, :) :: umatrix_k ! Full displacement matrix for atom k
      complex(rp), allocatable, dimension(:, :, :) :: temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10 ! Temporary matrices to store
      complex(rp), allocatable, dimension(:, :, :, :) :: int_all ! Result for all components on the Jijk tensor
      real(rp), allocatable, dimension(:, :) :: jijk_tot ! The value of Jijk for each energy point
      !
      ! Internal (local) variables
      real(rp), allocatable, dimension(:, :) :: screening_ort_i, screening_ort_j, screening_ort_k ! Screening constants (orthogonal representation)
      real(rp), allocatable, dimension(:, :) :: screening_can ! Screening constants (canonical representation)
      real(rp), dimension(4, 9) :: angles ! For the definition of xx, xy, xz, yx, yy, ... components
      real(rp), dimension(3) :: disp_direction, uni_disp
      real(rp) :: ws_radius ! The WS radius
      integer :: i, j, k, p, s, lll, lmaxi, lmaxj, lmaxk, njijk, nv

      ! Definition of each of the angles
      do p = 1, 9
         if (p .eq. 1) then ! xx
            angles(1, p) = 0.5_RP * pi ! theta
            angles(2, p) = 0.5_RP * pi ! theta'
            angles(3, p) = 0.0_RP ! phi
            angles(4, p) = 0.0_RP ! phi'
         else if (p .eq. 2) then ! xy
            angles(1, p) = 0.5_RP * pi ! theta
            angles(2, p) = 0.5_RP * pi ! theta'
            angles(3, p) = 0.0_RP ! phi
            angles(4, p) = 0.5_RP * pi ! phi'
         else if (p .eq. 3) then ! xz
            angles(1, p) = 0.5_RP * pi ! theta
            angles(2, p) = 0.0_RP ! theta'
            angles(3, p) = 0.0_RP ! phi
            angles(4, p) = 0.0_RP ! phi'
         else if (p .eq. 4) then ! yx
            angles(1, p) = 0.5_RP * pi ! theta
            angles(2, p) = 0.5_RP * pi ! theta'
            angles(3, p) = 0.5_RP * pi ! phi
            angles(4, p) = 0.0_RP ! phi'
         else if (p .eq. 5) then ! yy
            angles(1, p) = 0.5_RP * pi ! theta
            angles(2, p) = 0.5_RP * pi ! theta'
            angles(3, p) = 0.5_RP * pi ! phi
            angles(4, p) = 0.5_RP * pi ! phi'
         else if (p .eq. 6) then ! yz
            angles(1, p) = 0.5_RP * pi ! theta
            angles(2, p) = 0.0_RP ! theta'
            angles(3, p) = 0.5_RP * pi ! phi
            angles(4, p) = 0.0_RP ! phi'
         else if (p .eq. 7) then ! zx
            angles(1, p) = 0.0_RP ! theta
            angles(2, p) = 0.5_RP * pi ! theta'
            angles(3, p) = 0.0_RP ! phi
            angles(4, p) = 0.0_RP ! phi'
         else if (p .eq. 8) then ! zy
            angles(1, p) = 0.0_RP ! theta
            angles(2, p) = 0.5_RP * pi ! theta'
            angles(3, p) = 0.0_RP ! phi
            angles(4, p) = 0.5_RP * pi ! phi'
         else ! zz
            angles(1, p) = 0.0_RP ! theta
            angles(2, p) = 0.0_RP ! theta'
            angles(3, p) = 0.0_RP ! phi
            angles(4, p) = 0.0_RP ! phi'
         end if
      end do

      !--------------------
      ! Initial definitions
      !--------------------

      ! The scale factor for the Laplace equation as the wav informed by the user.
      ! Default units of wav in Å, so the resulting spin-lattice interaction will be in mRy/Å
      ws_radius = this%lattice%wav

      !---------------------------
      ! End of initial definitions
      !---------------------------

      ! Loop over all trios informed by the user
      do njijk = 1, this%lattice%njijk
         ! Obtain which atoms are in the current trio
         i = int(this%lattice%ijktrio(njijk, 1)) ! Atom number in the clust file, atom i
         j = int(this%lattice%ijktrio(njijk, 2)) ! Atom number in the clust file, atom j
         k = int(this%lattice%ijktrio(njijk, 3)) ! Atom number in the clust file, atom k
         ! Now obtain the lmax values for each of the atoms in the trio
         lmaxi = this%symbolic_atom(this%lattice%iz(i))%potential%lmax ! l_max for atom i
         lmaxj = this%symbolic_atom(this%lattice%iz(j))%potential%lmax ! l_max for atom j
         lmaxk = this%symbolic_atom(this%lattice%iz(k))%potential%lmax ! l_max for atom k
         ! Now allocate the necessary vectors and matrices
         allocate (screening_ort_i(0:lmaxi, 2))
         allocate (screening_ort_j(0:lmaxj, 2))
         allocate (screening_ort_k(0:lmaxk, 2))
         allocate (screening_can(0:lmaxi, 2))
         allocate (aux_gij(2 * (lmaxi + 1)**2, 2 * (lmaxi + 1)**2, size(this%en%ene)))
         allocate (aux_gji(2 * (lmaxj + 1)**2, 2 * (lmaxj + 1)**2, size(this%en%ene)))
         allocate (aux_gik(2 * (lmaxi + 1)**2, 2 * (lmaxi + 1)**2, size(this%en%ene)))
         allocate (aux_gki(2 * (lmaxk + 1)**2, 2 * (lmaxk + 1)**2, size(this%en%ene)))
         allocate (aux_gjk(2 * (lmaxj + 1)**2, 2 * (lmaxj + 1)**2, size(this%en%ene)))
         allocate (aux_gkj(2 * (lmaxk + 1)**2, 2 * (lmaxk + 1)**2, size(this%en%ene)))
         allocate (aux0_gij(2 * (lmaxi + 1)**2, 2 * (lmaxi + 1)**2, size(this%en%ene)))
         allocate (aux0_gji(2 * (lmaxj + 1)**2, 2 * (lmaxj + 1)**2, size(this%en%ene)))
         allocate (aux0_gik(2 * (lmaxi + 1)**2, 2 * (lmaxi + 1)**2, size(this%en%ene)))
         allocate (aux0_gki(2 * (lmaxk + 1)**2, 2 * (lmaxk + 1)**2, size(this%en%ene)))
         allocate (aux0_gjk(2 * (lmaxj + 1)**2, 2 * (lmaxj + 1)**2, size(this%en%ene)))
         allocate (aux0_gkj(2 * (lmaxk + 1)**2, 2 * (lmaxk + 1)**2, size(this%en%ene)))
         allocate (pmatrix_i(2 * (lmaxi + 1)**2, 2 * (lmaxi + 1)**2, size(this%en%ene)))
         allocate (pmatrix_j(2 * (lmaxj + 1)**2, 2 * (lmaxj + 1)**2, size(this%en%ene)))
         allocate (pmatrix_k(2 * (lmaxk + 1)**2, 2 * (lmaxk + 1)**2, size(this%en%ene)))
         allocate (pmatrix0_i(2 * (lmaxi + 1)**2, 2 * (lmaxi + 1)**2, size(this%en%ene)))
         allocate (pmatrix0_j(2 * (lmaxj + 1)**2, 2 * (lmaxj + 1)**2, size(this%en%ene)))
         allocate (pmatrix0_k(2 * (lmaxk + 1)**2, 2 * (lmaxk + 1)**2, size(this%en%ene)))
         allocate (umatrix_k(2 * (lmaxk + 1)**2, 2 * (lmaxk + 1)**2, size(this%en%ene)))
         allocate (deltap_i((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (deltap_j((lmaxj + 1)**2, (lmaxj + 1)**2, size(this%en%ene)))
         allocate (temp1((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (temp2((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (temp3((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (temp4((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (temp5((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (temp6((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (temp7((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (temp8((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (temp9((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (temp10((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene)))
         allocate (int_all((lmaxi + 1)**2, (lmaxi + 1)**2, size(this%en%ene), 9))
         allocate (jijk_tot(size(this%en%ene), 9))
         ! Screening constants for the orthogonal representation are the qpar (or gamma, in Turek's notation)
         screening_ort_i(:, :) = this%symbolic_atom(this%lattice%iz(i))%potential%qpar(:, :)
         screening_ort_j(:, :) = this%symbolic_atom(this%lattice%iz(j))%potential%qpar(:, :)
         screening_ort_k(:, :) = this%symbolic_atom(this%lattice%iz(k))%potential%qpar(:, :)
         ! Screening constants for the canonical representation (always alpha_rl = 0)
         screening_can(:, :) = 0.0_RP
         ! Set all matrices and arrays to zero first
         jijk_tot(:, :) = 0.0_RP
         disp_direction(:) = 0.0_RP
         uni_disp(:) = 0.0_RP
         aux_gij = czero; aux_gji = czero; aux_gik = czero
         aux_gki = czero; aux_gjk = czero; aux_gkj = czero
         aux0_gij = czero; aux0_gji = czero; aux0_gik = czero
         aux0_gki = czero; aux0_gjk = czero; aux0_gkj = czero
         pmatrix_i = czero; pmatrix_j = czero; pmatrix_k = czero
         pmatrix0_i = czero; pmatrix0_j = czero; pmatrix0_k = czero
         deltap_i = czero; deltap_j = czero; umatrix_k = czero
         int_all = czero
         ! Now get the direction for the respective trio (in Cartesian coordinates)
         disp_direction(1) = (1.0_RP * (this%lattice%ijktrio(njijk, 4))) ! x component
         disp_direction(2) = (1.0_RP * (this%lattice%ijktrio(njijk, 5))) ! y component
         disp_direction(3) = (1.0_RP * (this%lattice%ijktrio(njijk, 6))) ! z component
         ! Make the direction unitary
         uni_disp(1) = disp_direction(1) / norm2(disp_direction)
         uni_disp(2) = disp_direction(2) / norm2(disp_direction)
         uni_disp(3) = disp_direction(3) / norm2(disp_direction)
         ! Now calculate the potential functions P for atoms i, j and k
         ! in the orthogonal representation
         call this%symbolic_atom(this%lattice%iz(i))%p_matrix(pmatrix_i, lmaxi, this%en%ene)
         call this%symbolic_atom(this%lattice%iz(j))%p_matrix(pmatrix_j, lmaxj, this%en%ene)
         call this%symbolic_atom(this%lattice%iz(k))%p_matrix(pmatrix_k, lmaxk, this%en%ene)
         ! And next transform those P potential functions to the canonical representation
         ! (screening constants = 0)
         call this%symbolic_atom(this%lattice%iz(i))%transform_pmatrix(pmatrix_i, pmatrix0_i, screening_ort_i, screening_can)
         call this%symbolic_atom(this%lattice%iz(j))%transform_pmatrix(pmatrix_j, pmatrix0_j, screening_ort_j, screening_can)
         call this%symbolic_atom(this%lattice%iz(k))%transform_pmatrix(pmatrix_k, pmatrix0_k, screening_ort_k, screening_can)
         ! Calculate the displacement matrix U for atom k, with the canonical representation
         call this%symbolic_atom(this%lattice%iz(k))%udisp_matrix(umatrix_k, pmatrix0_k, uni_disp, lmaxk, ws_radius, size(this%en%ene))
         ! Now obtain the auxiliary Green's function for atoms i, j and k (and vice-versa), in the orthogonal representation
         call this%green%auxiliary_gij(this%green%gij(:, :, :, 3 * (njijk - 1) + 1), aux_gij, i, j)
         call this%green%auxiliary_gij(this%green%gji(:, :, :, 3 * (njijk - 1) + 1), aux_gji, j, i)
         call this%green%auxiliary_gij(this%green%gij(:, :, :, 3 * (njijk - 1) + 2), aux_gik, i, k)
         call this%green%auxiliary_gij(this%green%gji(:, :, :, 3 * (njijk - 1) + 2), aux_gki, k, i)
         call this%green%auxiliary_gij(this%green%gij(:, :, :, 3 * (njijk - 1) + 3), aux_gjk, j, k)
         call this%green%auxiliary_gij(this%green%gji(:, :, :, 3 * (njijk - 1) + 3), aux_gkj, k, j)
         ! Transform those auxiliary GF's to the canonical representation
         call this%green%transform_auxiliary_gij(pmatrix_i, pmatrix0_i, pmatrix_j, pmatrix0_j, aux_gij, aux0_gij, screening_ort_i, screening_can, i, j)
         call this%green%transform_auxiliary_gij(pmatrix_j, pmatrix0_j, pmatrix_i, pmatrix0_i, aux_gji, aux0_gji, screening_ort_j, screening_can, j, i)
         call this%green%transform_auxiliary_gij(pmatrix_i, pmatrix0_i, pmatrix_k, pmatrix0_k, aux_gik, aux0_gik, screening_ort_i, screening_can, i, k)
         call this%green%transform_auxiliary_gij(pmatrix_k, pmatrix0_k, pmatrix_i, pmatrix0_i, aux_gki, aux0_gki, screening_ort_k, screening_can, k, i)
         call this%green%transform_auxiliary_gij(pmatrix_j, pmatrix0_j, pmatrix_k, pmatrix0_k, aux_gjk, aux0_gjk, screening_ort_j, screening_can, j, k)
         call this%green%transform_auxiliary_gij(pmatrix_k, pmatrix0_k, pmatrix_j, pmatrix0_j, aux_gkj, aux0_gkj, screening_ort_k, screening_can, k, j)
         ! Now transform the pmatrices in DeltaP = P_up - P_dw (i and j are the atoms whose spin are tilted a bit)
         ! Everything in the canonical representation
         deltap_i(:, :, :) = pmatrix0_i(1:(lmaxi + 1)**2, 1:(lmaxi + 1)**2, :) - pmatrix0_i((lmaxi + 1)**2 + 1:2 * (lmaxi + 1)**2, (lmaxi + 1)**2 + 1:2 * (lmaxi + 1)**2, :)
         deltap_j(:, :, :) = pmatrix0_j(1:(lmaxj + 1)**2, 1:(lmaxj + 1)**2, :) - pmatrix0_j((lmaxj + 1)**2 + 1:2 * (lmaxj + 1)**2, (lmaxj + 1)**2 + 1:2 * (lmaxj + 1)**2, :)
         ! Now do the calculation of jijk_tot for all energy points
         temp1 = czero; temp2 = czero; temp3 = czero; temp4 = czero; temp5 = czero
         temp6 = czero; temp7 = czero; temp8 = czero; temp9 = czero; temp10 = czero
         do nv = 1, size(this%en%ene)
            temp1(:, :, nv) = matmul(umatrix_k((lmaxk + 1)**2 + 1:2 * (lmaxk + 1)**2, (lmaxk + 1)**2 + 1:2 * (lmaxk + 1)**2, nv), aux0_gki((lmaxk + 1)**2 + 1:2 * (lmaxk + 1)**2, (lmaxk + 1)**2 + 1:2 * (lmaxk + 1)**2, nv)) ! Uk_dwdw*gki_dwdw
            temp2(:, :, nv) = matmul(umatrix_k(1:(lmaxk + 1)**2, 1:(lmaxk + 1)**2, nv), aux0_gki(1:(lmaxk + 1)**2, 1:(lmaxk + 1)**2, nv)) ! Uk_upup*gki_upup
            temp3(:, :, nv) = matmul(deltap_i(:, :, nv), aux0_gij(1:(lmaxi + 1)**2, 1:(lmaxi + 1)**2, nv)) ! (P_up(i) - P_dw(i))*gij_upup
            temp4(:, :, nv) = matmul(deltap_j(:, :, nv), aux0_gjk(1:(lmaxj + 1)**2, 1:(lmaxj + 1)**2, nv)) ! (P_up(j) - P_dw(j))*gjk_upup
            temp5(:, :, nv) = matmul(umatrix_k(1:(lmaxk + 1)**2, 1:(lmaxk + 1)**2, nv), aux0_gkj(1:(lmaxk + 1)**2, 1:(lmaxk + 1)**2, nv)) ! Uk_upup*gkj_upup
            temp6(:, :, nv) = matmul(umatrix_k((lmaxk + 1)**2 + 1:2 * (lmaxk + 1)**2, (lmaxk + 1)**2 + 1:2 * (lmaxk + 1)**2, nv), aux0_gkj((lmaxk + 1)**2 + 1:2 * (lmaxk + 1)**2, (lmaxk + 1)**2 + 1:2 * (lmaxk + 1)**2, nv)) ! Uk_dwdw*gkj_dwdw
            temp7(:, :, nv) = matmul(deltap_j(:, :, nv), aux0_gji(1:(lmaxj + 1)**2, 1:(lmaxj + 1)**2, nv)) ! (P_up(j) - P_dw(j))*gji_upup
            temp8(:, :, nv) = matmul(deltap_i(:, :, nv), aux0_gij((lmaxi + 1)**2 + 1:2 * (lmaxi + 1)**2, (lmaxi + 1)**2 + 1:2 * (lmaxi + 1)**2, nv)) ! (P_up(i) - P_dw(i))*gij_dwdw
            temp9(:, :, nv) = matmul(deltap_j(:, :, nv), aux0_gjk((lmaxj + 1)**2 + 1:2 * (lmaxj + 1)**2, (lmaxj + 1)**2 + 1:2 * (lmaxj + 1)**2, nv)) ! (P_up(j) - P_dw(j))*gjk_dwdw
            temp10(:, :, nv) = matmul(deltap_j(:, :, nv), aux0_gji((lmaxj + 1)**2 + 1:2 * (lmaxj + 1)**2, (lmaxj + 1)**2 + 1:2 * (lmaxj + 1)**2, nv)) ! (P_up(j) - P_dw(j))*gji_dwdw
            ! Calculates the Jijk tensor
            do p = 1, 9
               int_all(:, :, nv, p) = scalar_multiply(matmul(temp3(:, :, nv), matmul(temp4(:, :, nv), temp2(:, :, nv))), cone * cos(angles(1, p)) * cos(angles(2, p))) + &
                                      scalar_multiply(matmul(temp8(:, :, nv), matmul(temp4(:, :, nv), temp2(:, :, nv))), cone * sin(angles(1, p)) * sin(angles(2, p)) * exp(i_unit * (angles(4, p) - angles(3, p)))) + &
                                      scalar_multiply(matmul(temp3(:, :, nv), matmul(temp9(:, :, nv), temp1(:, :, nv))), cone * sin(angles(1, p)) * sin(angles(2, p)) * exp(i_unit * (angles(3, p) - angles(4, p)))) + &
                                      scalar_multiply(matmul(temp8(:, :, nv), matmul(temp9(:, :, nv), temp1(:, :, nv))), cone * cos(angles(1, p)) * cos(angles(2, p))) + &
                                      scalar_multiply(matmul(temp3(:, :, nv), matmul(temp5(:, :, nv), temp10(:, :, nv))), cone * sin(angles(1, p)) * sin(angles(2, p)) * exp(i_unit * (angles(3, p) - angles(4, p)))) + &
                                      scalar_multiply(matmul(temp8(:, :, nv), matmul(temp6(:, :, nv), temp10(:, :, nv))), cone * cos(angles(1, p)) * cos(angles(2, p))) + &
                                      scalar_multiply(matmul(temp3(:, :, nv), matmul(temp5(:, :, nv), temp7(:, :, nv))), cone * cos(angles(1, p)) * cos(angles(2, p))) + &
                                      scalar_multiply(matmul(temp8(:, :, nv), matmul(temp6(:, :, nv), temp7(:, :, nv))), cone * sin(angles(1, p)) * sin(angles(2, p)) * exp(i_unit * (angles(4, p) - angles(3, p))))
               jijk_tot(nv, p) = imtrace(int_all(:, :, nv, p)) * (0.5_RP)
            end do
         end do
         ! Integrate jijk_tot with respect to energy
         do p = 1, 9
            call simpson_f(this%jijk(p), this%en%ene, this%en%fermi, this%en%nv1, jijk_tot(:, p), .true., .false., 0.0D0)
         end do
         ! Write the result
         write (*, *) 'Jijk tensor between trio ', i, ',', j, ' and ', k, 'is (in meV/a.u.)'
         write (*, '(A, "(", F7.4, ", ", F7.4, ", ", F7.4, ")")') ' Displacement vector: ', uni_disp(1), uni_disp(2), uni_disp(3)
         write (*, '(3F14.9)') this%jijk * (1.0D3 / 8.0D0 / pi) * (13.605693122994D0 / 1.8897261246D0)
         ! Now deallocate all vectors and matrices allocated
         deallocate (screening_ort_i)
         deallocate (screening_ort_j)
         deallocate (screening_ort_k)
         deallocate (screening_can)
         deallocate (aux_gij)
         deallocate (aux_gji)
         deallocate (aux_gik)
         deallocate (aux_gki)
         deallocate (aux_gjk)
         deallocate (aux_gkj)
         deallocate (aux0_gij)
         deallocate (aux0_gji)
         deallocate (aux0_gik)
         deallocate (aux0_gki)
         deallocate (aux0_gjk)
         deallocate (aux0_gkj)
         deallocate (pmatrix_i)
         deallocate (pmatrix_j)
         deallocate (pmatrix_k)
         deallocate (pmatrix0_i)
         deallocate (pmatrix0_j)
         deallocate (pmatrix0_k)
         deallocate (umatrix_k)
         deallocate (deltap_i)
         deallocate (deltap_j)
         deallocate (temp1)
         deallocate (temp2)
         deallocate (temp3)
         deallocate (temp4)
         deallocate (temp5)
         deallocate (temp6)
         deallocate (temp7)
         deallocate (temp8)
         deallocate (temp9)
         deallocate (temp10)
         deallocate (int_all)
         deallocate (jijk_tot)
      end do

   end subroutine calculate_jijk

   subroutine calculate_gilbert_damping(this)
      !
      class(exchange) :: this
      !
      complex(rp), allocatable, dimension(:, :, :) :: Aij, Aji !
      complex(rp), allocatable, dimension(:, :, :) :: temp1, temp2, temp3 ! Temporary matrices to store
      complex(rp), allocatable, dimension(:, :, :) :: tmati, tmatj ! Torque matrices

      ! Local Variables
      integer :: i, j, k, l, m, njij, nv, lmaxi, lmaxj, ief
      real(rp) :: spin_i, orbital_i, gfac, diff
      real(rp), allocatable, dimension(:, :) :: dtott, dtottim, total_damping
      real(rp), dimension(3, 3) :: damping_tensor, damping_tensor_im

      if (this%control%nsp .eq. 2) then ! check if spin-orbit (l.s) is enabled

         call this%hamiltonian%torque_operator_collinear() ! calculate the torque operator for all NTYPE
         allocate (total_damping(9, size(this%en%ene)))
         total_damping(:, :) = 0.0_RP

         ! Open the files for writing the results
         open (UNIT=103, FILE='damping-energy.out', STATUS='replace', ACTION='write')

         do njij = 1, this%lattice%njij
            ! Obtain which are the atoms in the current pair
            i = this%lattice%ijpair(njij, 1) ! Atom number in the clust file, atom i
            j = this%lattice%ijpair(njij, 2) ! Atom number in the clust file, atom j
            ! Find the lmax value for each atom
            lmaxi = this%symbolic_atom(this%lattice%iz(i))%potential%lmax ! l_max for atom i
            lmaxj = this%symbolic_atom(this%lattice%iz(j))%potential%lmax ! l_max for atom j
            ! Allocate all necessary arrays and quantities
            allocate (Aij(2 * (lmaxi + 1)**2, 2 * (lmaxi + 1)**2, size(this%en%ene)))
            allocate (Aji(2 * (lmaxj + 1)**2, 2 * (lmaxj + 1)**2, size(this%en%ene)))
            allocate (temp1(2 * (lmaxi + 1)**2, 2 * (lmaxi + 1)**2, size(this%en%ene)))
            allocate (temp2(2 * (lmaxi + 1)**2, 2 * (lmaxi + 1)**2, size(this%en%ene)))
            allocate (temp3(2 * (lmaxi + 1)**2, 2 * (lmaxi + 1)**2, size(this%en%ene)))
            allocate (tmati(2 * (lmaxi + 1)**2, 2 * (lmaxi + 1)**2, 3))
            allocate (tmatj(2 * (lmaxj + 1)**2, 2 * (lmaxj + 1)**2, 3))
            allocate (dtott(9, size(this%en%ene)))
            allocate (dtottim(9, size(this%en%ene)))
            ! Make all necessary values equal to zero
            Aij(:, :, :) = czero; Aji(:, :, :) = czero
            temp1(:, :, :) = czero; temp2(:, :, :) = czero; temp3(:, :, :) = czero
            tmati(:, :, :) = czero; tmatj(:, :, :) = czero
            ! Obtain the torque matrices for the atoms i and j in the collinear mode
            tmati(:, :, :) = this%hamiltonian%tmat(:, :, :, this%lattice%iz(i))
            tmatj(:, :, :) = this%hamiltonian%tmat(:, :, :, this%lattice%iz(j))
            ! Obtain the magnetic moment (spin and orbital) of the i-th atom
            ! TO-DO: include orbital moment
            spin_i = 0.0_RP; orbital_i = 0.0_RP
            do k = 0, lmaxi ! azimuthal quantum number index
               spin_i = spin_i + this%symbolic_atom(this%lattice%iz(i))%potential%ql(1, k, 1) - &
                        this%symbolic_atom(this%lattice%iz(i))%potential%ql(1, k, 2)
               !do l = -k, k ! magnetic quantum number index
               !  orbital_i = orbital_i + l*(this%symbolic_atom(this%lattice%iz(i))%potential%ql(1, k, 1) + &
               !                             this%symbolic_atom(this%lattice%iz(i))%potential%ql(1, k, 2))
               !end do
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
            ief = 0; diff = 1000.0_RP
            do nv = 1, size(this%en%ene)
               if (abs(this%en%ene(nv) - this%en%fermi) .lt. diff) then
                  diff = abs(this%en%ene(nv) - this%en%fermi)
                  ief = nv
               end if
            end do
            !write(*, *) 'closest ief=', ief, 'energy=', this%en%ene(ief), 'fermi=', this%en%fermi
            ! Finally calculates the damping term and write the results on screen
            write (*, *) 'Damping tensor between pair', i, 'and ', j, 'is '
            write (*, '(3F14.9)') (-0.25_RP) * (2.0_RP / (pi * spin_i)) * dtott(:, ief)
            write (*, *) 'Imaginary part of the tensor between pair', i, 'and ', j, 'is '
            write (*, '(3F14.9)') (-0.25_RP) * (2.0_RP / (pi * spin_i)) * dtottim(:, ief)

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
         write (103, *) '#Energy (E-Ef)         #xx         #xy           #xz           #yx           #yy           #yz           #zx           #zy           #zz'
         do nv = 1, size(this%en%ene)
            write (103, '(10F14.9)') this%en%ene(nv) - this%en%fermi, (-0.25_RP) * (2.0_RP / (pi * spin_i)) * total_damping(:, nv)
         end do

         ! Deallocate all
         deallocate (total_damping)

         ! Close the files
         close (103)

      end if

   end subroutine calculate_gilbert_damping

   subroutine calculate_dij(this)
      class(exchange) :: this
   end subroutine

   function dGdG_Jnc(this, di, dj, njij, nv) result(Jmat)

      implicit none
      class(exchange) :: this
      complex(rp), dimension(9, 9), intent(in) :: di
      complex(rp), dimension(9, 9), intent(in) :: dj
      integer, intent(in) :: njij
      integer, intent(in) :: nv
      complex(rp), dimension(9, 9) :: Jmat

      complex(rp), dimension(9, 9) :: tmat1, tmat2, tmat3

      tmat1 = matmul(di, this%green%Ginmag(1:9, 1:9, nv, njij))
      tmat2 = matmul(dj, this%green%Gjnmag(1:9, 1:9, nv, njij))
      tmat3 = matmul(tmat1, tmat2)
      tmat1 = matmul(di, this%green%Gix(1:9, 1:9, nv, njij))
      tmat2 = matmul(dj, this%green%Gjx(1:9, 1:9, nv, njij))
      tmat3 = tmat3 - matmul(tmat1, tmat2)
      tmat1 = matmul(di, this%green%Giy(1:9, 1:9, nv, njij))
      tmat2 = matmul(dj, this%green%Gjy(1:9, 1:9, nv, njij))
      tmat3 = tmat3 - matmul(tmat1, tmat2)
      tmat1 = matmul(di, this%green%Giz(1:9, 1:9, nv, njij))
      tmat2 = matmul(dj, this%green%Gjz(1:9, 1:9, nv, njij))
      Jmat = tmat3 - matmul(tmat1, tmat2)

      return
   end function dGdG_Jnc

   function dGdG_Dnc(this, di, dj, njij, nv) result(Dmats)

      implicit none
      class(exchange) :: this
      complex(rp), dimension(9, 9), intent(in) :: di
      complex(rp), dimension(9, 9), intent(in) :: dj
      integer, intent(in) :: njij
      integer, intent(in) :: nv
      complex(rp), dimension(9, 9, 3) :: Dmats

      integer :: k
      complex(rp), dimension(9, 9) :: tmat1, tmat2, tmat3, tmat4

      do k = 1, 3
         tmat1 = matmul(di, this%green%Ginmag(1:9, 1:9, nv, njij))
         if (k == 1) tmat2 = matmul(dj, this%green%Gjx(1:9, 1:9, nv, njij))
         if (k == 2) tmat2 = matmul(dj, this%green%Gjy(1:9, 1:9, nv, njij))
         if (k == 3) tmat2 = matmul(dj, this%green%Gjz(1:9, 1:9, nv, njij))
         tmat3 = matmul(tmat1, tmat2)
         tmat1 = matmul(dj, this%green%Gjnmag(1:9, 1:9, nv, njij))
         if (k == 1) tmat2 = matmul(di, this%green%Gix(1:9, 1:9, nv, njij))
         if (k == 2) tmat2 = matmul(di, this%green%Giy(1:9, 1:9, nv, njij))
         if (k == 3) tmat2 = matmul(di, this%green%Giz(1:9, 1:9, nv, njij))
         tmat4 = matmul(tmat1, tmat2)
         Dmats(:, :, k) = tmat3 - tmat4
         ! Take the imaginary part of the trace
      end do

      return
   end function dGdG_Dnc

   function dGdG_Anc(this, di, dj, njij, nv) result(Amats)

      implicit none
      class(exchange) :: this
      complex(rp), dimension(9, 9), intent(in) :: di
      complex(rp), dimension(9, 9), intent(in) :: dj
      integer, intent(in) :: njij
      integer, intent(in) :: nv
      complex(rp), dimension(9, 9, 3, 3) :: Amats

      integer :: k, l
      complex(rp), dimension(9, 9) :: tmat1, tmat2, tmat3, tmat4

      do k = 1, 3
         do l = 1, 3
            if (k == 1) tmat1 = matmul(di, this%green%Gix(1:9, 1:9, nv, njij))
            if (k == 2) tmat1 = matmul(di, this%green%Giy(1:9, 1:9, nv, njij))
            if (k == 3) tmat1 = matmul(di, this%green%Giz(1:9, 1:9, nv, njij))
            if (l == 1) tmat2 = matmul(dj, this%green%Gjx(1:9, 1:9, nv, njij))
            if (l == 2) tmat2 = matmul(dj, this%green%Gjy(1:9, 1:9, nv, njij))
            if (l == 3) tmat2 = matmul(dj, this%green%Gjz(1:9, 1:9, nv, njij))
            tmat3 = matmul(tmat1, tmat2)
            if (k == 1) tmat1 = matmul(dj, this%green%Gjx(1:9, 1:9, nv, njij))
            if (k == 2) tmat1 = matmul(dj, this%green%Gjy(1:9, 1:9, nv, njij))
            if (k == 3) tmat1 = matmul(dj, this%green%Gjz(1:9, 1:9, nv, njij))
            if (l == 1) tmat2 = matmul(di, this%green%Gix(1:9, 1:9, nv, njij))
            if (l == 2) tmat2 = matmul(di, this%green%Giy(1:9, 1:9, nv, njij))
            if (l == 3) tmat2 = matmul(di, this%green%Giz(1:9, 1:9, nv, njij))
            tmat4 = matmul(tmat1, tmat2)
            Amats(:, :, k, l) = 0.5D0 * (tmat3 + tmat4)
         end do
      end do

      return
   end function dGdG_Anc

   subroutine calculate_aij(this)
      class(exchange) :: this
   end subroutine

   subroutine calculate_exchange(this)
      use mpi_mod
      class(exchange) :: this
      complex(rp), dimension(9, 9) :: dmat1, dmat2, tmat1, tmat2, tmat3, tmat4, smat
      complex(rp), dimension(this%en%channels_ldos + 10) :: jtot, jcomp
      complex(rp), dimension(this%en%channels_ldos + 10, 3) :: jjtot
      complex(rp), dimension(this%en%channels_ldos + 10, 3, 3) :: itot
      real(rp), dimension(this%en%channels_ldos + 10) :: y
      real(rp), dimension(3, 3) :: jtens
      integer :: nv, i, j, k, l, njij, funit, iostatus
      logical :: isopen

      integer :: njij_glob
      complex(rp), dimension(9, 9) :: jmat
      complex(rp), dimension(9, 9, 3) :: dmmats
      complex(rp), dimension(9, 9, 3, 3) :: amats

      real(rp), dimension(:, :), allocatable :: T_comm_xc

      ! Communications array for MPI
      ! row = [ jij dij aij]
      allocate (T_comm_xc(13, this%lattice%njij))
      T_comm_xc = 0.0_RP

      inquire (unit=20, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file jij.out: Unit 20 is already open', __FILE__, __LINE__)
      else
         open (unit=20, file='jij.out')
      end if
      inquire (unit=30, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file dij.out: Unit 30 is already open', __FILE__, __LINE__)
      else
         open (unit=30, file='dij.out')
      end if
      inquire (unit=40, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file aij.out: Unit 40 is already open', __FILE__, __LINE__)
      else
         open (unit=40, file='aij.out')
      end if
      inquire (unit=60, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file jtens.out: Unit 40 is already open', __FILE__, __LINE__)
      else
         open (unit=60, file='jtens.out')
      end if

      ! First calculate exchange interactions in parallel, without printing
      jtot = 0.0D0; jjtot = 0.0D0
      do njij_glob = start_atom, end_atom
         njij = g2l_map(njij_glob)
         i = this%lattice%ijpair(njij_glob, 1) ! Atom number in the clust file, atom i
         j = this%lattice%ijpair(njij_glob, 2) ! Atom number in the clust file, atom j

         do nv = 1, this%en%channels_ldos + 10
            call this%symbolic_atom(this%lattice%iz(i))%d_matrix(dmat1, this%en%ene(nv))
            call this%symbolic_atom(this%lattice%iz(j))%d_matrix(dmat2, this%en%ene(nv))
            !J
            jmat = dGdG_Jnc(this, dmat1, dmat2, njij, nv)
            jtot(nv) = imtrace9(jmat)
            jcomp(nv) = trace9(jmat)

            ! DMI
            dmmats = dGdG_Dnc(this, dmat1, dmat2, njij, nv)
            do k = 1, 3
               jjtot(nv, k) = rtrace9(dmmats(:, :, k))
            end do

            ! I Tensor
            amats = dGdG_Anc(this, dmat1, dmat2, njij, nv)
            do k = 1, 3
               do l = 1, 3
                  itot(nv, k, l) = imtrace9(amats(:, :, k, l))
               end do
            end do
         end do
         ! Jij integration
         this%jij = 0.0D0
         call simpson_f(this%jij, this%en%ene, this%en%fermi, this%en%nv1, real(jtot), .true., .false., 0.0D0)
         T_comm_xc(1, njij_glob) = this%jij * 1.0D3 / 4.0D0 / pi

         ! Dij integration
         this%dmi = 0.0D0
         do k = 1, 3
            y(:) = 0.0D0
            y(:) = real(jjtot(:, k))
            call simpson_f(this%dmi(k), this%en%ene, this%en%fermi, this%en%nv1, y(:), .true., .false., 0.0D0)
         end do
         T_comm_xc(2:4, njij_glob) = this%dmi * 1.0D3 / 4.0D0 / pi

         ! I tensor integration
         this%aij = 0.0D0
         do k = 1, 3
            do l = 1, 3
               y(:) = 0.0D0
               y(:) = real(itot(:, k, l))
               call simpson_f(this%aij(k, l), this%en%ene, this%en%fermi, this%en%nv1, y(:), .true., .false., 0.0D0)
            end do
         end do
         T_comm_xc(5:13, njij_glob) = reshape(this%aij * 1.0D3 / 4.0D0 / pi, [9])

      end do

#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, T_comm_xc, product(shape(T_comm_xc)), &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

! Print exchange interactions only using MPI rank 0
      if (rank == 0) then
         do njij_glob = 1, this%lattice%njij
            !do njij_glob=start_atom, end_atom
            njij = g2l_map(njij_glob)
            i = this%lattice%ijpair(njij_glob, 1) ! Atom number in the clust file, atom i
            j = this%lattice%ijpair(njij_glob, 2) ! Atom number in the clust file, atom j

            write (*, *) 'Atom', i, 'coordinates:', this%lattice%cr(:, i)
            write (*, *) 'Atom', j, 'coordinates:', this%lattice%cr(:, j)

            ! Jij
            this%jij = T_comm_xc(1, njij_glob)
            write (*, *) 'Jij between pair', i, 'and ', j, 'is ', this%jij
            ! Dij
            this%dmi = T_comm_xc(2:4, njij_glob)
            write (*, *) 'Dij between pair', i, 'and ', j, 'is ', this%dmi
            ! Aij
            this%aij = reshape(T_comm_xc(5:13, njij_glob), [3, 3])
            write (*, *) 'Iij between pair', i, 'and ', j, 'is'
            print '(3f12.6)', this%aij(1, :)
            print '(3f12.6)', this%aij(2, :)
            print '(3f12.6)', this%aij(3, :)

            ! Write the jij.out file
            write (20, '(2i8,2x,3f12.6,2x,1f12.6,1x,f12.6)') &
            & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), this%jij, norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            ! Write the dij.out file
            write (30, '(2i8,2x,3f12.6,2x,3f12.6,1x,f12.6)') &
            & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), this%dmi, norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            ! Write the aij.out file
            write (40, '(2i8,2x,3f12.6,2x,9f12.6,1x,f12.6)') &
            & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), this%aij, norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            ! Compute and write the jtens.out file
            ! J on the diagonal
            jtens = 0.0D0
            jtens(1, 1) = this%jij
            jtens(2, 2) = jtens(1, 1)
            jtens(3, 3) = jtens(1, 1)
            ! DM on skew form
            jtens(1, 2) = this%dmi(3)
            jtens(2, 1) = -jtens(1, 2)
            jtens(1, 3) = -this%dmi(2)
            jtens(3, 1) = -jtens(1, 3)
            jtens(2, 3) = this%dmi(1)
            jtens(3, 2) = -jtens(2, 3)
            ! A as full tensor
            jtens(:, :) = jtens(:, :) + this%aij(:, :)
            write (*, *) 'J tensor between pair', i, 'and ', j, 'is'
            print '(3f12.6)', jtens(1, :)
            print '(3f12.6)', jtens(2, :)
            print '(3f12.6)', jtens(3, :)

         end do
      end if

      close (20)
      close (30)
      close (40)
      close (50)

      deallocate (T_comm_xc)

#ifdef USE_MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

   end subroutine calculate_exchange

   subroutine calculate_exchange_rs2pao(this)
      use mpi_mod
      class(exchange) :: this
      complex(rp), dimension(9, 9) :: dmat1, dmat2, tmat1, tmat2, tmat3, tmat4, smat
      complex(rp), dimension(this%en%channels_ldos + 10) :: jtot
      complex(rp), dimension(this%en%channels_ldos + 10, 3) :: jjtot
      real(rp), dimension(this%en%channels_ldos + 10) :: y
      integer :: nv, i, j, k, njij, l
      logical :: isopen
      integer :: njij_glob

      complex(rp), dimension(this%en%channels_ldos + 10, 3, 3) :: itot
      complex(rp), dimension(9, 9) :: jmat
      complex(rp), dimension(9, 9, 3) :: dmmats
      complex(rp), dimension(9, 9, 3, 3) :: amats
      real(rp), dimension(:, :), allocatable :: T_comm_xc

      ! Communications array for MPI
      ! row = [ jij dij aij]
      allocate (T_comm_xc(13, this%lattice%njij))
      T_comm_xc = 0.0_RP

      inquire (unit=20, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file jij.out: Unit 20 is already open', __FILE__, __LINE__)
      else
         open (unit=20, file='jij.out')
      end if
      inquire (unit=30, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file dij.out: Unit 30 is already open', __FILE__, __LINE__)
      else
         open (unit=30, file='dij.out')
      end if

      do njij_glob = start_atom, end_atom
         njij = g2l_map(njij_glob)
         jtot = 0.0D0; jjtot = 0.0D0
         dmat1 = 0.0D0; dmat2 = 0.0D0; tmat1 = 0.0D0; tmat2 = 0.0D0; tmat3 = 0.0D0; smat = 0.0D0

         i = this%lattice%ijpair(njij_glob, 1) ! Atom number in the clust file, atom i
         j = this%lattice%ijpair(njij_glob, 2) ! Atom number in the clust file, atom j
         dmat1 = real(this%hamiltonian%ee(1:9, 1:9, 1, this%lattice%iz(i)) - this%hamiltonian%ee(10:18, 10:18, 1, this%lattice%iz(i)))
         dmat2 = real(this%hamiltonian%ee(1:9, 1:9, 1, this%lattice%iz(j)) - this%hamiltonian%ee(10:18, 10:18, 1, this%lattice%iz(j)))

         do nv = 1, this%en%channels_ldos + 10
            !J
            jmat = dGdG_Jnc(this, dmat1, dmat2, njij, nv)
            jtot(nv) = imtrace9(jmat)

            ! DMI
            dmmats = dGdG_Dnc(this, dmat1, dmat2, njij, nv)
            do k = 1, 3
               jjtot(nv, k) = rtrace9(dmmats(:, :, k))
            end do

            ! I Tensor
            amats = dGdG_Anc(this, dmat1, dmat2, njij, nv)
            do k = 1, 3
               do l = 1, 3
                  itot(nv, k, l) = imtrace9(amats(:, :, k, l))
               end do
            end do
         end do
         this%jij = 0.0D0
      !!! do nv=1,this%en%channels_ldos+10
      !!!   call simpson_f(this%jij,this%en%ene,this%en%ene(nv),this%en%nv1,real(jtot),.true.,.false.,0.0d0)
      !!!   write(902,*) this%en%ene(nv)-this%en%fermi, this%jij*1.0d3/4.0d0/pi
      !!! end do
         call simpson_f(this%jij, this%en%ene, this%en%fermi, this%en%nv1, real(jtot), .true., .false., 0.0D0)
         T_comm_xc(1, njij_glob) = this%jij * 1.0D3 / 4.0D0 / pi

         this%dmi = 0.0D0
         do k = 1, 3
            y(:) = 0.0D0
            y(:) = real(jjtot(:, k))
            call simpson_f(this%dmi(k), this%en%ene, this%en%fermi, this%en%nv1, y(:), .true., .false., 0.0D0)
         end do
         T_comm_xc(2:4, njij_glob) = this%dmi * 1.0D3 / 4.0D0 / pi

         ! I tensor integration
         this%aij = 0.0D0
         do k = 1, 3
            do l = 1, 3
               y(:) = 0.0D0
               y(:) = real(itot(:, k, l))
               call simpson_f(this%aij(k, l), this%en%ene, this%en%fermi, this%en%nv1, y(:), .true., .false., 0.0D0)
            end do
         end do
         T_comm_xc(5:13, njij_glob) = reshape(this%aij * 1.0D3 / 4.0D0 / pi, [9])

      end do

#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, T_comm_xc, product(shape(T_comm_xc)), &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      if (rank == 0) then
         do njij_glob = 1, this%lattice%njij
            !do njij_glob=start_atom, end_atom
            ! njij = g2l_map(njij_glob)

            i = this%lattice%ijpair(njij_glob, 1) ! Atom number in the clust file, atom i
            j = this%lattice%ijpair(njij_glob, 2) ! Atom number in the clust file, atom j

            write (*, *) 'Atom', i, 'coordinates:', this%lattice%cr(:, i)
            write (*, *) 'Atom', j, 'coordinates:', this%lattice%cr(:, j)
            ! Jij
            this%jij = T_comm_xc(1, njij_glob)
            write (*, *) 'Jij between pair', i, 'and ', j, 'is ', this%jij
            ! Dij
            this%dmi = T_comm_xc(2:4, njij_glob)
            write (*, *) 'Dij between pair', i, 'and ', j, 'is ', this%dmi
            ! Aij
            this%aij = reshape(T_comm_xc(5:13, njij_glob), [3, 3])
            write (*, *) 'Iij between pair', i, 'and ', j, 'is'
            print '(3f12.6)', this%aij(1, :)
            print '(3f12.6)', this%aij(2, :)
            print '(3f12.6)', this%aij(3, :)

            write (*, *) 'Fermi Energy', this%en%fermi
            !
            write (20, '(2i8,2x,3f12.6,2x,1f12.6,1x,f12.6)') &
            & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), this%jij, norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            write (30, '(2i8,2x,3f12.6,2x,3f12.6,1x,f12.6)') &
            & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), this%dmi, norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            write (40, '(2i8,2x,3f12.6,2x,9f12.6,1x,f12.6)') &
            & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), this%aij, norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))

         end do
      end if
      close (40)
      deallocate (T_comm_xc)
#ifdef USE_MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
   end subroutine calculate_exchange_rs2pao

   subroutine calculate_exchange_gauss_legendre(this)
      use mpi_mod
      class(exchange) :: this
      complex(rp), dimension(9, 9) :: dmat1, dmat2, tmat1, tmat2, tmat3, tmat4, smat
      complex(rp), dimension(64) :: jtot, jcomp
      complex(rp), dimension(this%en%channels_ldos + 10, 3) :: jjtot
      real(rp), dimension(this%en%channels_ldos + 10) :: y
      integer :: nv, i, j, k, njij, orbj, orbm, l
      real(rp), dimension(64) :: x, w
      logical :: isopen

      integer :: njij_glob
      complex(rp), dimension(64, 3, 3) :: itot
      complex(rp), dimension(9, 9) :: jmat
      complex(rp), dimension(9, 9, 3) :: dmmats
      complex(rp), dimension(9, 9, 3, 3) :: amats
      real(rp), dimension(:, :), allocatable :: T_comm_xc

      ! Communications array for MPI
      ! row = [ jij dij aij]
      allocate (T_comm_xc(13, this%lattice%njij))
      T_comm_xc = 0.0_RP

      inquire (unit=20, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file jij.outt: Unit 20 is already open', __FILE__, __LINE__)
      else
         open (unit=20, file='jij.out')
      end if
      inquire (unit=30, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file dij.outt: Unit 30 is already open', __FILE__, __LINE__)
      else
         open (unit=30, file='dij.out')
      end if
      inquire (unit=40, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file aij.out: Unit 40 is already open', __FILE__, __LINE__)
      else
         open (unit=40, file='aij.out')
      end if
      inquire (unit=60, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file jtens.out: Unit 40 is already open', __FILE__, __LINE__)
      else
         open (unit=60, file='jtens.out')
      end if

      ! Clone gij_eta to gij to enable use of generalized routines
      call this%green%gij_eta_to_gij()

      ! Find the Gauss Legendre roots and weights
      call gauss_legendre(64, 0.0_RP, 1.0_RP, x, w)

      !do njij=1,this%lattice%njij
      do njij_glob = start_atom, end_atom
         njij = g2l_map(njij_glob)
         jtot = 0.0D0; jjtot = 0.0D0
         dmat1 = 0.0D0; dmat2 = 0.0D0; tmat1 = 0.0D0; tmat2 = 0.0D0; tmat3 = 0.0D0; smat = 0.0D0

         i = this%lattice%ijpair(njij_glob, 1) ! Atom number in the clust file, atom i
         j = this%lattice%ijpair(njij_glob, 2) ! Atom number in the clust file, atom j

         dmat1 = real(this%hamiltonian%ee(1:9, 1:9, 1, this%lattice%iz(i)) - this%hamiltonian%ee(10:18, 10:18, 1, this%lattice%iz(i)))
         dmat2 = real(this%hamiltonian%ee(1:9, 1:9, 1, this%lattice%iz(j)) - this%hamiltonian%ee(10:18, 10:18, 1, this%lattice%iz(j)))

         jtot = 0.0D0

         do nv = 1, 64
            !J
            jmat = dGdG_Jnc(this, dmat1, dmat2, njij, nv)
            jmat = (jmat * w(nv)) / (x(nv) * x(nv))
            jtot(nv) = rtrace9(jmat)

            ! DMI
            dmmats = dGdG_Dnc(this, dmat1, dmat2, njij, nv)
            dmmats = (dmmats * w(nv)) / (x(nv) * x(nv))
            do k = 1, 3
               jjtot(nv, k) = imtrace9(dmmats(:, :, k))
            end do

            ! I Tensor
            amats = dGdG_Anc(this, dmat1, dmat2, njij, nv)
            amats = (amats * w(nv)) / (x(nv) * x(nv))
            do k = 1, 3
               do l = 1, 3
                  itot(nv, k, l) = rtrace9(amats(:, :, k, l))
               end do
            end do

         end do

         ! Jij summation
         this%jij = -sum(real(jtot(:)))
         T_comm_xc(1, njij_glob) = this%jij * 1.0D3 / 4.0D0 / pi

         ! DMI
         do k = 1, 3
            this%dmi(k) = sum(jjtot(:, k))
         end do
         T_comm_xc(2:4, njij_glob) = this%dmi * 1.0D3 / 4.0D0 / pi

         ! I tensor integration
         this%aij = 0.0D0
         do k = 1, 3
            do l = 1, 3
               this%aij(k, l) = -sum(itot(:, k, l))
            end do
         end do
         T_comm_xc(5:13, njij_glob) = reshape(this%aij * 1.0D3 / 4.0D0 / pi, [9])

      end do

#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, T_comm_xc, product(shape(T_comm_xc)), &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      if (rank == 0) then
         do njij_glob = 1, this%lattice%njij
            !do njij_glob=start_atom, end_atom
            !njij = g2l_map(njij_glob)
            jtot = 0.0D0; jjtot = 0.0D0
            dmat1 = 0.0D0; dmat2 = 0.0D0; tmat1 = 0.0D0; tmat2 = 0.0D0; tmat3 = 0.0D0; smat = 0.0D0

            i = this%lattice%ijpair(njij_glob, 1) ! Atom number in the clust file, atom i
            j = this%lattice%ijpair(njij_glob, 2) ! Atom number in the clust file, atom j

            write (*, *) 'Atom', i, 'coordinates:', this%lattice%cr(:, i), 'Atom Type', this%lattice%iz(i)
            write (*, *) 'Atom', j, 'coordinates:', this%lattice%cr(:, j), 'Atom Type', this%lattice%iz(j)
            write (*, *) 'Distance = ', norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))

            ! Jij
            this%jij = T_comm_xc(1, njij_glob)
            write (*, *) 'Jij between pair', i, 'and ', j, 'is ', this%jij

            ! Dij
            this%dmi = T_comm_xc(2:4, njij_glob)
            write (*, *) 'Dij between pair', i, 'and ', j, 'is ', this%dmi

            ! Aij
            this%aij = reshape(T_comm_xc(5:13, njij_glob), [3, 3])
            write (*, *) 'Iij between pair', i, 'and ', j, 'is'
            print '(3f12.6)', this%aij(1, :)
            print '(3f12.6)', this%aij(2, :)
            print '(3f12.6)', this%aij(3, :)
            !
            write (20, '(2i8,2x,3f12.6,2x,1f12.6,1x,f12.6)') &
            & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), (this%jij), norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            write (30, '(2i8,2x,3f12.6,2x,3f12.6,1x,f12.6)') &
            & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), (this%dmi), norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            write (40, '(2i8,2x,3f12.6,2x,9f12.6,1x,f12.6)') &
            & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), this%aij, norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))

         end do
      end if
      close (20)
      close (30)
      close (40)
      deallocate (T_comm_xc)
#ifdef USE_MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
   end subroutine calculate_exchange_gauss_legendre

end module exchange_mod
