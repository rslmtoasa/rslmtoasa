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
   use basis_mod, only: lmax_basis
#ifdef USE_MPI
   use mpi
#endif
   use cfd
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   use string_mod
   implicit none

   private
   !> Linked from the exchange_dynamics submodule (calculate_pairwise)
   public :: exchange_pair, exchange_twoindex_pair, exchange_finalize, exchange_twoindex_finalize, real_mesh_nodes, onsite_exchange_delta, spin_blocks_hcpx
   !> Energy nodes at which G is evaluated and the way an integrand sampled on them is integrated.
   type, public :: energy_nodes
      !> Complex energies z (Ry) at which G is evaluated
      complex(rp), allocatable :: z(:)
      !> Contour (Gauss-Legendre) nodes; otherwise the real Simpson mesh
      logical :: contour = .false.
      !> Real mesh: the Fermi-function energy, temperature and point count passed to simpson_f
      real(rp) :: ef = 0.0_rp, temperature = 0.0_rp
      integer :: nv1 = 0
      !> Contour: Gauss-Legendre abscissae and weights, applied as (M*w)/(x*x) by the kernel
      real(rp), allocatable :: x(:), w(:)
   contains
      procedure :: integrate => nodes_integrate
   end type energy_nodes

   !> Module´s main structure
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

      logical :: use_onsite_delta = .false.

      !> General variables
      !> Heisenberg exchange
      real(rp) :: jij, jijcd, jijsd, jijcc, jijsc
      !> Heisenberg exchange tensor (obtained with auxiliary GF´s)
      real(rp), dimension(9) :: jij_aux
      !> Total Heisenberg exchange J0 = sum_j J0j (obtained with auxiliary GF´s)
      real(rp) :: jij00_aux
      !> Spin-lattice couplings tensor
      real(rp), dimension(9) :: jijk
      !> Dzyaloshinskii-Moriya interaction
      real(rp), dimension(3) :: dmi, dmisc, dmicc
      !> Anisotropy
      real(rp), dimension(3, 3) :: aij, aijsd, aijsc
   contains
      procedure :: calculate_jij
      procedure :: calculate_dij
      procedure :: calculate_aij
      procedure :: calculate_jij_auxgreen
      procedure :: calculate_jijk
      procedure :: calculate_exchange_gauss_legendre
      procedure :: calculate_exchange_rs2pao
      procedure :: calculate_pairwise
      procedure :: restore_to_default
      final destructor
   end type

   interface exchange
      procedure :: constructor
   end interface exchange

   interface
      module subroutine calculate_pairwise(this, stream, do_exchange, do_damping, do_inertia, do_chi0)
         class(exchange) :: this
         logical, intent(in) :: stream, do_exchange, do_damping, do_inertia, do_chi0
      end subroutine calculate_pairwise
   end interface

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
      ! The exchange kernels index the spin blocks as 1:9 / 10:18 (spd basis only).
      if (lmax_basis /= 2) call g_logger%fatal('exchange: implemented only for an spd basis (lmax_basis = 2), got lmax_basis = '// &
                                               int2str(lmax_basis), __FILE__, __LINE__)
      ! initialize constraining if requested for exchange calculations
      if (associated(obj%control)) then
         if (obj%control%constraints_enable) then
            call initialize_cfd(obj%lattice%nrec, 1, obj%control%constraints_i_cons, obj%control%constraints_code_prefac)
         end if
      end if
   end function constructor

   subroutine destructor(this)
      type(exchange) :: this
   end subroutine destructor

   subroutine restore_to_default(this)
      class(exchange) :: this

      this%jij = 0.0d0
      this%jijcd = 0.0d0
      this%jijsd = 0.0d0
      this%jijcc = 0.0d0
      this%jijsc = 0.0d0
      this%jij_aux(:) = 0.0_rp
      this%jij00_aux = 0.0_rp
      this%jijk(:) = 0.0_rp
      this%dmi(:) = 0.0d0
      this%dmisc(:) = 0.0d0
      this%dmicc(:) = 0.0d0
      this%aij(:, :) = 0.0d0
      this%aijsd(:, :) = 0.0d0
      this%aijsc(:, :) = 0.0d0
   end subroutine restore_to_default

   subroutine calculate_jij(this)
      class(exchange) :: this
      complex(rp), dimension(9, 9) :: dmat1, dmat2, jmat
      real(rp), dimension(this%en%channels_ldos + 10) :: jtot
      integer :: nv, i, j, njij

      jtot = 0.0d0
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
         call simpson_f(this%jij, this%en%ene, this%en%ene(this%en%enpt + 1), this%en%nv1, jtot, .true., .false., this%en%temperature)
         write (*, *) 'Jij between pair', i, 'and ', j, 'is ', this%jij*1.0d3/4.0d0/pi
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
      complex(rp), allocatable, dimension(:, :, :) :: aux_gij, aux_gji ! Auxiliary Green´s functions
      complex(rp), allocatable, dimension(:, :, :, :) :: int_all ! Result for the components of the Jij tensor (when i != j)
      complex(rp), allocatable, dimension(:, :, :) :: int_j00 ! Result for J00 (when i = j)
      complex(rp), allocatable, dimension(:, :, :) :: temp1, temp2, temp3, temp4 ! Temporary matrices to store
      real(rp), allocatable, dimension(:, :) :: jtot_aux ! Jij total from the auxiliary GF´s formalism (when i != j)
      real(rp), allocatable, dimension(:) :: jtot_00 ! J00 total from the auxiliary GF´s formalism (when i = j)
      real(rp), dimension(4, 9) :: angles ! For the definition of xx, xy, xz, yx, yy, ... components
      integer :: nv, i, j, k, njij, lmaxi, lmaxj ! Internal (local) variables

      ! Definition of each of the angles
      do k = 1, 9
         if (k .eq. 1) then ! xx
            angles(1, k) = 0.5_rp*pi ! theta
            angles(2, k) = 0.5_rp*pi ! theta´
            angles(3, k) = 0.0_rp ! phi
            angles(4, k) = 0.0_rp ! phi´
         else if (k .eq. 2) then ! xy
            angles(1, k) = 0.5_rp*pi ! theta
            angles(2, k) = 0.5_rp*pi ! theta´
            angles(3, k) = 0.0_rp ! phi
            angles(4, k) = 0.5_rp*pi ! phi´
         else if (k .eq. 3) then ! xz
            angles(1, k) = 0.5_rp*pi ! theta
            angles(2, k) = 0.0_rp ! theta´
            angles(3, k) = 0.0_rp ! phi
            angles(4, k) = 0.0_rp ! phi´
         else if (k .eq. 4) then ! yx
            angles(1, k) = 0.5_rp*pi ! theta
            angles(2, k) = 0.5_rp*pi ! theta´
            angles(3, k) = 0.5_rp*pi ! phi
            angles(4, k) = 0.0_rp ! phi´
         else if (k .eq. 5) then ! yy
            angles(1, k) = 0.5_rp*pi ! theta
            angles(2, k) = 0.5_rp*pi ! theta´
            angles(3, k) = 0.5_rp*pi ! phi
            angles(4, k) = 0.5_rp*pi ! phi´
         else if (k .eq. 6) then ! yz
            angles(1, k) = 0.5_rp*pi ! theta
            angles(2, k) = 0.0_rp ! theta´
            angles(3, k) = 0.5_rp*pi ! phi
            angles(4, k) = 0.0_rp ! phi´
         else if (k .eq. 7) then ! zx
            angles(1, k) = 0.0_rp ! theta
            angles(2, k) = 0.5_rp*pi ! theta´
            angles(3, k) = 0.0_rp ! phi
            angles(4, k) = 0.0_rp ! phi´
         else if (k .eq. 8) then ! zy
            angles(1, k) = 0.0_rp ! theta
            angles(2, k) = 0.5_rp*pi ! theta´
            angles(3, k) = 0.0_rp ! phi
            angles(4, k) = 0.5_rp*pi ! phi´
         else ! zz
            angles(1, k) = 0.0_rp ! theta
            angles(2, k) = 0.0_rp ! theta´
            angles(3, k) = 0.0_rp ! phi
            angles(4, k) = 0.0_rp ! phi´
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
         allocate (pmatrix_i(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
         allocate (pmatrix_j(2*(lmaxj + 1)**2, 2*(lmaxj + 1)**2, size(this%en%ene)))
         allocate (aux_gij(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
         allocate (aux_gji(2*(lmaxj + 1)**2, 2*(lmaxj + 1)**2, size(this%en%ene)))
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
         jtot_aux(:, :) = 0.0_rp
         jtot_00(:) = 0.0_rp
         pmatrix_i = czero; pmatrix_j = czero
         aux_gij = czero; aux_gji = czero
         deltap_i = czero; deltap_j = czero
         int_all = czero; int_j00 = czero
         ! Now calculate the potential functions P for atoms i and j
         call this%symbolic_atom(this%lattice%iz(i))%p_matrix(pmatrix_i, lmaxi, this%en%ene)
         call this%symbolic_atom(this%lattice%iz(j))%p_matrix(pmatrix_j, lmaxj, this%en%ene)
         ! Now obtain the auxiliar Green´s function for atoms i and j (and vice-versa)
         call this%green%auxiliary_gij(this%green%gij(:, :, :, njij), aux_gij, i, j)
         call this%green%auxiliary_gij(this%green%gji(:, :, :, njij), aux_gji, j, i)
         ! Now transform the pmatrices in DeltaP = P_up - P_dw
         deltap_i(:, :, :) = pmatrix_i(1:(lmaxi + 1)**2, 1:(lmaxi + 1)**2, :) - pmatrix_i((lmaxi + 1)**2 + 1:2*(lmaxi + 1)**2, (lmaxi + 1)**2 + 1:2*(lmaxi + 1)**2, :)
         deltap_j(:, :, :) = pmatrix_j(1:(lmaxj + 1)**2, 1:(lmaxj + 1)**2, :) - pmatrix_j((lmaxj + 1)**2 + 1:2*(lmaxj + 1)**2, (lmaxj + 1)**2 + 1:2*(lmaxj + 1)**2, :)
         ! Now do the calculation of all the components of the jtot tensor
         temp1 = czero; temp2 = czero
         temp3 = czero; temp4 = czero
         do nv = 1, size(this%en%ene)
            ! Calculates Jij tensor
            if (i .ne. j) then
               temp1(:, :, nv) = matmul(deltap_i(:, :, nv), aux_gij(1:(lmaxi + 1)**2, 1:(lmaxi + 1)**2, nv)) ! (P_up(i) - P_dw(i))*gij_upup
               temp2(:, :, nv) = matmul(deltap_j(:, :, nv), aux_gji((lmaxj + 1)**2 + 1:2*(lmaxj + 1)**2, (lmaxj + 1)**2 + 1:2*(lmaxj + 1)**2, nv)) ! (P_up(j) - P_dw(j))*gji_dwdw
               temp3(:, :, nv) = matmul(deltap_i(:, :, nv), aux_gij((lmaxi + 1)**2 + 1:2*(lmaxi + 1)**2, (lmaxi + 1)**2 + 1:2*(lmaxi + 1)**2, nv)) ! (P_up(i) - P_dw(i))*gij_dwdw
               temp4(:, :, nv) = matmul(deltap_j(:, :, nv), aux_gji(1:(lmaxj + 1)**2, 1:(lmaxj + 1)**2, nv)) ! (P_up(j) - P_dw(j))*gji_upup
               do k = 1, 9 ! run over the tensor elements
                  int_all(:, :, nv, k) = scalar_multiply(matmul(temp1(:, :, nv), temp4(:, :, nv)), cone*cos(angles(1, k))*cos(angles(2, k))) + &
                                         scalar_multiply(matmul(temp3(:, :, nv), temp4(:, :, nv)), cone*sin(angles(1, k))*sin(angles(2, k))*exp(i_unit*(angles(4, k) - angles(3, k)))) + &
                                         scalar_multiply(matmul(temp1(:, :, nv), temp2(:, :, nv)), cone*sin(angles(1, k))*sin(angles(2, k))*exp(i_unit*(angles(3, k) - angles(4, k)))) + &
                                         scalar_multiply(matmul(temp3(:, :, nv), temp2(:, :, nv)), cone*cos(angles(1, k))*cos(angles(2, k)))
                  jtot_aux(nv, k) = imtrace(int_all(:, :, nv, k))*(0.5_rp)
               end do
               ! Calculates J00 otherwise
            else ! i = j
               temp1(:, :, nv) = matmul(deltap_i(:, :, nv), aux_gij(1:(lmaxi + 1)**2, 1:(lmaxi + 1)**2, nv)) ! (P_up(i) - P_dw(i))*gii_upup
               temp2(:, :, nv) = matmul(deltap_j(:, :, nv), aux_gji((lmaxj + 1)**2 + 1:2*(lmaxj + 1)**2, (lmaxj + 1)**2 + 1:2*(lmaxj + 1)**2, nv)) ! (P_up(i) - P_dw(i))*gii_dwdw
               ! (P_up(i) - P_dw(i))*(gii_upup - gii_dwdw)
               temp3(:, :, nv) = matmul(deltap_i(:, :, nv), aux_gij(1:(lmaxi + 1)**2, 1:(lmaxi + 1)**2, nv) - aux_gji((lmaxj + 1)**2 + 1:2*(lmaxj + 1)**2, (lmaxj + 1)**2 + 1:2*(lmaxj + 1)**2, nv))
               int_j00(:, :, nv) = matmul(temp1(:, :, nv), temp2(:, :, nv)) + temp3(:, :, nv)
               jtot_00(nv) = imtrace(int_j00(:, :, nv))*(-1.0_rp)
            end if
         end do
         ! Integrate jtot components with respect to energy
         if (i .ne. j) then
            do k = 1, 9
               call simpson_f(this%jij_aux(k), this%en%ene, this%en%ene(this%en%enpt + 1), this%en%nv1, jtot_aux(:, k), .true., .false., this%en%temperature)
            end do
         else ! i = j
            call simpson_f(this%jij00_aux, this%en%ene, this%en%ene(this%en%enpt + 1), this%en%nv1, jtot_00, .true., .false., this%en%temperature)
         end if
         ! Write the result
         if (i .ne. j) then
            write (*, *) 'Jij_aux tensor between pair', i, 'and ', j, 'is '
            write (*, '(3F14.9)') this%jij_aux*1.0d3/4.0d0/pi
            write (*, *) 'Dij_zz_aux between pair', i, 'and', j, 'is ', 0.5_rp*(this%jij_aux(2) - this%jij_aux(4))*(1.0d3/4.0d0/pi)
         else ! i = j
            write (*, *) 'J0_aux is', this%jij00_aux*1.0d3/4.0d0/pi
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
      complex(rp), allocatable, dimension(:, :, :) :: aux_gij, aux_gji, aux_gik, aux_gki, aux_gjk, aux_gkj ! Auxiliary GF´s (in the orthogonal representaton)
      complex(rp), allocatable, dimension(:, :, :) :: aux0_gij, aux0_gji, aux0_gik, aux0_gki, aux0_gjk, aux0_gkj ! Same as above, but in the canonical representation
      complex(rp), allocatable, dimension(:, :, :) :: pmatrix_i, pmatrix_j, pmatrix_k ! P matrices (potential function) for atoms i, j, k (orthogonal representation)
      complex(rp), allocatable, dimension(:, :, :) :: pmatrix0_i, pmatrix0_j, pmatrix0_k ! P matrix (potential function) of atom k in the canonical representation
      complex(rp), allocatable, dimension(:, :, :) :: deltap_i, deltap_j ! P_up - P_dw matrices (half of the size, as it does not consider the spin index)
      complex(rp), allocatable, dimension(:, :, :) :: umatrix_k ! Full displacement matrix for atom k
      complex(rp), allocatable, dimension(:, :, :) :: pcubic_k ! P matrix of atom k in the cubic basis of disp_matrix
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
            angles(1, p) = 0.5_rp*pi ! theta
            angles(2, p) = 0.5_rp*pi ! theta´
            angles(3, p) = 0.0_rp ! phi
            angles(4, p) = 0.0_rp ! phi´
         else if (p .eq. 2) then ! xy
            angles(1, p) = 0.5_rp*pi ! theta
            angles(2, p) = 0.5_rp*pi ! theta´
            angles(3, p) = 0.0_rp ! phi
            angles(4, p) = 0.5_rp*pi ! phi´
         else if (p .eq. 3) then ! xz
            angles(1, p) = 0.5_rp*pi ! theta
            angles(2, p) = 0.0_rp ! theta´
            angles(3, p) = 0.0_rp ! phi
            angles(4, p) = 0.0_rp ! phi´
         else if (p .eq. 4) then ! yx
            angles(1, p) = 0.5_rp*pi ! theta
            angles(2, p) = 0.5_rp*pi ! theta´
            angles(3, p) = 0.5_rp*pi ! phi
            angles(4, p) = 0.0_rp ! phi´
         else if (p .eq. 5) then ! yy
            angles(1, p) = 0.5_rp*pi ! theta
            angles(2, p) = 0.5_rp*pi ! theta´
            angles(3, p) = 0.5_rp*pi ! phi
            angles(4, p) = 0.5_rp*pi ! phi´
         else if (p .eq. 6) then ! yz
            angles(1, p) = 0.5_rp*pi ! theta
            angles(2, p) = 0.0_rp ! theta´
            angles(3, p) = 0.5_rp*pi ! phi
            angles(4, p) = 0.0_rp ! phi´
         else if (p .eq. 7) then ! zx
            angles(1, p) = 0.0_rp ! theta
            angles(2, p) = 0.5_rp*pi ! theta´
            angles(3, p) = 0.0_rp ! phi
            angles(4, p) = 0.0_rp ! phi´
         else if (p .eq. 8) then ! zy
            angles(1, p) = 0.0_rp ! theta
            angles(2, p) = 0.5_rp*pi ! theta´
            angles(3, p) = 0.0_rp ! phi
            angles(4, p) = 0.5_rp*pi ! phi´
         else ! zz
            angles(1, p) = 0.0_rp ! theta
            angles(2, p) = 0.0_rp ! theta´
            angles(3, p) = 0.0_rp ! phi
            angles(4, p) = 0.0_rp ! phi´
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
         if (lmaxk > 2) call g_logger%fatal('exchange%calculate_jijk: disp_matrix f-shell order differs from the structure-constant f order; not supported', &
                                            __FILE__, __LINE__)
         ! Now allocate the necessary vectors and matrices
         allocate (screening_ort_i(0:lmaxi, 2))
         allocate (screening_ort_j(0:lmaxj, 2))
         allocate (screening_ort_k(0:lmaxk, 2))
         allocate (screening_can(0:lmaxi, 2))
         allocate (aux_gij(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
         allocate (aux_gji(2*(lmaxj + 1)**2, 2*(lmaxj + 1)**2, size(this%en%ene)))
         allocate (aux_gik(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
         allocate (aux_gki(2*(lmaxk + 1)**2, 2*(lmaxk + 1)**2, size(this%en%ene)))
         allocate (aux_gjk(2*(lmaxj + 1)**2, 2*(lmaxj + 1)**2, size(this%en%ene)))
         allocate (aux_gkj(2*(lmaxk + 1)**2, 2*(lmaxk + 1)**2, size(this%en%ene)))
         allocate (aux0_gij(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
         allocate (aux0_gji(2*(lmaxj + 1)**2, 2*(lmaxj + 1)**2, size(this%en%ene)))
         allocate (aux0_gik(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
         allocate (aux0_gki(2*(lmaxk + 1)**2, 2*(lmaxk + 1)**2, size(this%en%ene)))
         allocate (aux0_gjk(2*(lmaxj + 1)**2, 2*(lmaxj + 1)**2, size(this%en%ene)))
         allocate (aux0_gkj(2*(lmaxk + 1)**2, 2*(lmaxk + 1)**2, size(this%en%ene)))
         allocate (pmatrix_i(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
         allocate (pmatrix_j(2*(lmaxj + 1)**2, 2*(lmaxj + 1)**2, size(this%en%ene)))
         allocate (pmatrix_k(2*(lmaxk + 1)**2, 2*(lmaxk + 1)**2, size(this%en%ene)))
         allocate (pmatrix0_i(2*(lmaxi + 1)**2, 2*(lmaxi + 1)**2, size(this%en%ene)))
         allocate (pmatrix0_j(2*(lmaxj + 1)**2, 2*(lmaxj + 1)**2, size(this%en%ene)))
         allocate (pmatrix0_k(2*(lmaxk + 1)**2, 2*(lmaxk + 1)**2, size(this%en%ene)))
         allocate (umatrix_k(2*(lmaxk + 1)**2, 2*(lmaxk + 1)**2, size(this%en%ene)))
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
         ! Screening constants for the orthogonal representation are the qpar (or gamma, in Turek´s notation)
         screening_ort_i(:, :) = this%symbolic_atom(this%lattice%iz(i))%potential%qpar(:, :)
         screening_ort_j(:, :) = this%symbolic_atom(this%lattice%iz(j))%potential%qpar(:, :)
         screening_ort_k(:, :) = this%symbolic_atom(this%lattice%iz(k))%potential%qpar(:, :)
         ! Screening constants for the canonical representation (always alpha_rl = 0)
         screening_can(:, :) = 0.0_rp
         ! Set all matrices and arrays to zero first
         jijk_tot(:, :) = 0.0_rp
         disp_direction(:) = 0.0_rp
         uni_disp(:) = 0.0_rp
         aux_gij = czero; aux_gji = czero; aux_gik = czero
         aux_gki = czero; aux_gjk = czero; aux_gkj = czero
         aux0_gij = czero; aux0_gji = czero; aux0_gik = czero
         aux0_gki = czero; aux0_gjk = czero; aux0_gkj = czero
         pmatrix_i = czero; pmatrix_j = czero; pmatrix_k = czero
         pmatrix0_i = czero; pmatrix0_j = czero; pmatrix0_k = czero
         deltap_i = czero; deltap_j = czero; umatrix_k = czero
         int_all = czero
         ! Now get the direction for the respective trio (in Cartesian coordinates)
         disp_direction(1) = (1.0_rp*(this%lattice%ijktrio(njijk, 4))) ! x component
         disp_direction(2) = (1.0_rp*(this%lattice%ijktrio(njijk, 5))) ! y component
         disp_direction(3) = (1.0_rp*(this%lattice%ijktrio(njijk, 6))) ! z component
         ! Make the direction unitary
         uni_disp(1) = disp_direction(1)/norm2(disp_direction)
         uni_disp(2) = disp_direction(2)/norm2(disp_direction)
         uni_disp(3) = disp_direction(3)/norm2(disp_direction)
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
         ! disp_matrix is real and in cubic order (l <= 2) while P and the Green functions are spherical: build U in the cubic basis
         ! Calculate the displacement matrix U for atom k, with the canonical representation
         pcubic_k = pmatrix0_k
         call spin_blocks_hcpx(pcubic_k, 'sph2cart')
         call this%symbolic_atom(this%lattice%iz(k))%udisp_matrix(umatrix_k, pcubic_k, uni_disp, lmaxk, ws_radius, size(this%en%ene))
         call spin_blocks_hcpx(umatrix_k, 'cart2sph')
         ! Now obtain the auxiliary Green´s function for atoms i, j and k (and vice-versa), in the orthogonal representation
         call this%green%auxiliary_gij(this%green%gij(:, :, :, 3*(njijk - 1) + 1), aux_gij, i, j)
         call this%green%auxiliary_gij(this%green%gji(:, :, :, 3*(njijk - 1) + 1), aux_gji, j, i)
         call this%green%auxiliary_gij(this%green%gij(:, :, :, 3*(njijk - 1) + 2), aux_gik, i, k)
         call this%green%auxiliary_gij(this%green%gji(:, :, :, 3*(njijk - 1) + 2), aux_gki, k, i)
         call this%green%auxiliary_gij(this%green%gij(:, :, :, 3*(njijk - 1) + 3), aux_gjk, j, k)
         call this%green%auxiliary_gij(this%green%gji(:, :, :, 3*(njijk - 1) + 3), aux_gkj, k, j)
         ! Transform those auxiliary GF´s to the canonical representation
         call this%green%transform_auxiliary_gij(pmatrix_i, pmatrix0_i, pmatrix_j, pmatrix0_j, aux_gij, aux0_gij, screening_ort_i, screening_can, i, j)
         call this%green%transform_auxiliary_gij(pmatrix_j, pmatrix0_j, pmatrix_i, pmatrix0_i, aux_gji, aux0_gji, screening_ort_j, screening_can, j, i)
         call this%green%transform_auxiliary_gij(pmatrix_i, pmatrix0_i, pmatrix_k, pmatrix0_k, aux_gik, aux0_gik, screening_ort_i, screening_can, i, k)
         call this%green%transform_auxiliary_gij(pmatrix_k, pmatrix0_k, pmatrix_i, pmatrix0_i, aux_gki, aux0_gki, screening_ort_k, screening_can, k, i)
         call this%green%transform_auxiliary_gij(pmatrix_j, pmatrix0_j, pmatrix_k, pmatrix0_k, aux_gjk, aux0_gjk, screening_ort_j, screening_can, j, k)
         call this%green%transform_auxiliary_gij(pmatrix_k, pmatrix0_k, pmatrix_j, pmatrix0_j, aux_gkj, aux0_gkj, screening_ort_k, screening_can, k, j)
         ! Now transform the pmatrices in DeltaP = P_up - P_dw (i and j are the atoms whose spin are tilted a bit)
         ! Everything in the canonical representation
         deltap_i(:, :, :) = pmatrix0_i(1:(lmaxi + 1)**2, 1:(lmaxi + 1)**2, :) - pmatrix0_i((lmaxi + 1)**2 + 1:2*(lmaxi + 1)**2, (lmaxi + 1)**2 + 1:2*(lmaxi + 1)**2, :)
         deltap_j(:, :, :) = pmatrix0_j(1:(lmaxj + 1)**2, 1:(lmaxj + 1)**2, :) - pmatrix0_j((lmaxj + 1)**2 + 1:2*(lmaxj + 1)**2, (lmaxj + 1)**2 + 1:2*(lmaxj + 1)**2, :)
         ! Now do the calculation of jijk_tot for all energy points
         temp1 = czero; temp2 = czero; temp3 = czero; temp4 = czero; temp5 = czero
         temp6 = czero; temp7 = czero; temp8 = czero; temp9 = czero; temp10 = czero
         do nv = 1, size(this%en%ene)
            temp1(:, :, nv) = matmul(umatrix_k((lmaxk + 1)**2 + 1:2*(lmaxk + 1)**2, (lmaxk + 1)**2 + 1:2*(lmaxk + 1)**2, nv), aux0_gki((lmaxk + 1)**2 + 1:2*(lmaxk + 1)**2, (lmaxk + 1)**2 + 1:2*(lmaxk + 1)**2, nv)) ! Uk_dwdw*gki_dwdw
            temp2(:, :, nv) = matmul(umatrix_k(1:(lmaxk + 1)**2, 1:(lmaxk + 1)**2, nv), aux0_gki(1:(lmaxk + 1)**2, 1:(lmaxk + 1)**2, nv)) ! Uk_upup*gki_upup
            temp3(:, :, nv) = matmul(deltap_i(:, :, nv), aux0_gij(1:(lmaxi + 1)**2, 1:(lmaxi + 1)**2, nv)) ! (P_up(i) - P_dw(i))*gij_upup
            temp4(:, :, nv) = matmul(deltap_j(:, :, nv), aux0_gjk(1:(lmaxj + 1)**2, 1:(lmaxj + 1)**2, nv)) ! (P_up(j) - P_dw(j))*gjk_upup
            temp5(:, :, nv) = matmul(umatrix_k(1:(lmaxk + 1)**2, 1:(lmaxk + 1)**2, nv), aux0_gkj(1:(lmaxk + 1)**2, 1:(lmaxk + 1)**2, nv)) ! Uk_upup*gkj_upup
            temp6(:, :, nv) = matmul(umatrix_k((lmaxk + 1)**2 + 1:2*(lmaxk + 1)**2, (lmaxk + 1)**2 + 1:2*(lmaxk + 1)**2, nv), aux0_gkj((lmaxk + 1)**2 + 1:2*(lmaxk + 1)**2, (lmaxk + 1)**2 + 1:2*(lmaxk + 1)**2, nv)) ! Uk_dwdw*gkj_dwdw
            temp7(:, :, nv) = matmul(deltap_j(:, :, nv), aux0_gji(1:(lmaxj + 1)**2, 1:(lmaxj + 1)**2, nv)) ! (P_up(j) - P_dw(j))*gji_upup
            temp8(:, :, nv) = matmul(deltap_i(:, :, nv), aux0_gij((lmaxi + 1)**2 + 1:2*(lmaxi + 1)**2, (lmaxi + 1)**2 + 1:2*(lmaxi + 1)**2, nv)) ! (P_up(i) - P_dw(i))*gij_dwdw
            temp9(:, :, nv) = matmul(deltap_j(:, :, nv), aux0_gjk((lmaxj + 1)**2 + 1:2*(lmaxj + 1)**2, (lmaxj + 1)**2 + 1:2*(lmaxj + 1)**2, nv)) ! (P_up(j) - P_dw(j))*gjk_dwdw
            temp10(:, :, nv) = matmul(deltap_j(:, :, nv), aux0_gji((lmaxj + 1)**2 + 1:2*(lmaxj + 1)**2, (lmaxj + 1)**2 + 1:2*(lmaxj + 1)**2, nv)) ! (P_up(j) - P_dw(j))*gji_dwdw
            ! Calculates the Jijk tensor
            do p = 1, 9
               int_all(:, :, nv, p) = scalar_multiply(matmul(temp3(:, :, nv), matmul(temp4(:, :, nv), temp2(:, :, nv))), cone*cos(angles(1, p))*cos(angles(2, p))) + &
                                      scalar_multiply(matmul(temp8(:, :, nv), matmul(temp4(:, :, nv), temp2(:, :, nv))), cone*sin(angles(1, p))*sin(angles(2, p))*exp(i_unit*(angles(4, p) - angles(3, p)))) + &
                                      scalar_multiply(matmul(temp3(:, :, nv), matmul(temp9(:, :, nv), temp1(:, :, nv))), cone*sin(angles(1, p))*sin(angles(2, p))*exp(i_unit*(angles(3, p) - angles(4, p)))) + &
                                      scalar_multiply(matmul(temp8(:, :, nv), matmul(temp9(:, :, nv), temp1(:, :, nv))), cone*cos(angles(1, p))*cos(angles(2, p))) + &
                                      scalar_multiply(matmul(temp3(:, :, nv), matmul(temp5(:, :, nv), temp10(:, :, nv))), cone*sin(angles(1, p))*sin(angles(2, p))*exp(i_unit*(angles(3, p) - angles(4, p)))) + &
                                      scalar_multiply(matmul(temp8(:, :, nv), matmul(temp6(:, :, nv), temp10(:, :, nv))), cone*cos(angles(1, p))*cos(angles(2, p))) + &
                                      scalar_multiply(matmul(temp3(:, :, nv), matmul(temp5(:, :, nv), temp7(:, :, nv))), cone*cos(angles(1, p))*cos(angles(2, p))) + &
                                      scalar_multiply(matmul(temp8(:, :, nv), matmul(temp6(:, :, nv), temp7(:, :, nv))), cone*sin(angles(1, p))*sin(angles(2, p))*exp(i_unit*(angles(4, p) - angles(3, p))))
               jijk_tot(nv, p) = imtrace(int_all(:, :, nv, p))*(0.5_rp)
            end do
         end do
         ! Integrate jijk_tot with respect to energy
         do p = 1, 9
            call simpson_f(this%jijk(p), this%en%ene, this%en%ene(this%en%enpt + 1), this%en%nv1, jijk_tot(:, p), .true., .false., this%en%temperature)
         end do
         ! Write the result
         write (*, *) 'Jijk tensor between trio ', i, ',', j, ' and ', k, 'is (in meV/a.u.)'
         write (*, '(A, "(", F7.4, ", ", F7.4, ", ", F7.4, ")")') ' Displacement vector: ', uni_disp(1), uni_disp(2), uni_disp(3)
         write (*, '(3F14.9)') this%jijk*(1.0d3/8.0d0/pi)*(13.605693122994d0/1.8897261246d0)
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
         deallocate (pcubic_k)
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

   !> Applies hcpx to the four orbital blocks (spin up/down) of every matrix in mat(2n, 2n, :).
   subroutine spin_blocks_hcpx(mat, transformation)
      complex(rp), intent(inout) :: mat(:, :, :)
      character(len=*), intent(in) :: transformation
      complex(rp) :: blk(size(mat, 1)/2, size(mat, 1)/2)
      integer :: n, ie, a, b

      n = size(mat, 1)/2
      do ie = 1, size(mat, 3)
         do a = 0, 1
            do b = 0, 1
               blk = mat(a*n + 1:(a + 1)*n, b*n + 1:(b + 1)*n, ie)
               call hcpx(blk, transformation)
               mat(a*n + 1:(a + 1)*n, b*n + 1:(b + 1)*n, ie) = blk
            end do
         end do
      end do
   end subroutine spin_blocks_hcpx

   subroutine calculate_dij(this)
      class(exchange) :: this
   end subroutine

   function dGdG(this, di, dj, gij, gji) result(Jmat)

      implicit none
      class(exchange) :: this
      complex(rp), dimension(9, 9), intent(in) :: di
      complex(rp), dimension(9, 9), intent(in) :: dj
      complex(rp), dimension(9, 9), intent(in) :: gij, gji
      complex(rp), dimension(9, 9) :: Jmat, tmat1, tmat2

      tmat1 = matmul(di, gij)
      tmat2 = matmul(dj, gji)
      Jmat = matmul(tmat1, tmat2)

      return
   end function dGdG

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
            Amats(:, :, k, l) = 0.5d0*(tmat3 + tmat4)
         end do
      end do

      return
   end function dGdG_Anc

   subroutine calculate_aij(this)
      class(exchange) :: this
   end subroutine

   !> Use the real part of the zero-image onsite spin difference, as in exchange_p2rs.
   function onsite_exchange_delta(this, site) result(delta)
      class(exchange), intent(in) :: this
      integer, intent(in) :: site
      complex(rp) :: delta(9,9), onsite(18,18)
      integer :: itype

      itype = this%lattice%iz(site)
      if (site <= this%lattice%nmax) then
         onsite = this%hamiltonian%hall(:,:,1,site)
      else
         onsite = this%hamiltonian%ee(:,:,1,itype)
      end if
      onsite = onsite + this%hamiltonian%lsham(:,:,itype)
      delta = real(onsite(:9,:9) - onsite(10:,10:), rp)
   end function onsite_exchange_delta

   !> Spin-orbit (so) and first-order (fo) J, D, A and their parts of pair (i, j) read from Green-function slot `slot`, in meV.
   subroutine exchange_twoindex_pair(this, nodes, slot, i, j, xcso, xcfo, xcparts)
      class(exchange) :: this
      type(energy_nodes), intent(in) :: nodes
      integer, intent(in) :: slot, i, j
      real(rp), intent(out) :: xcso(13), xcfo(13), xcparts(28)
      complex(rp), dimension(9, 9) :: dmat1, dmat2, tmat1, tmat2, tmat3, tmat4
      complex(rp), dimension(size(nodes%z)) :: jtotso, jtotfo, jsd, jcd, jsc, jcc
      complex(rp), dimension(size(nodes%z), 3) :: jjtotso, jjtotfo, dsc, dcc
      complex(rp), dimension(size(nodes%z), 3, 3) :: itotso, itotfo, isd, isc
      real(rp), dimension(size(nodes%z)) :: y
      integer :: nv, k, l

      jtotso = 0.0d0; jtotfo = 0.0d0; jjtotso = 0.0d0; jjtotfo = 0.0d0; itotso = 0.0d0; itotfo = 0.0d0
      jcd = 0.0d0; jsd = 0.0d0; jcc = 0.0d0; jsc = 0.0d0; dsc = 0.0d0; dcc = 0.0d0; isd = 0.0d0; isc = 0.0d0
      if (this%use_onsite_delta) then
         dmat1 = onsite_exchange_delta(this, i)
         dmat2 = onsite_exchange_delta(this, j)
      end if
      do nv = 1, this%en%channels_ldos + 10
         if (.not. this%use_onsite_delta) then
            call this%symbolic_atom(this%lattice%iz(i))%d_matrix(dmat1, real(nodes%z(nv)))
            call this%symbolic_atom(this%lattice%iz(j))%d_matrix(dmat2, real(nodes%z(nv)))
         end if
         !J
         ! Charge density part of Jij
         tmat1 = 0.0d0
         tmat1 = dGdG(this, dmat1, dmat2, this%green%g00ij(1:9, 1:9, nv, slot), this%green%g00ji(1:9, 1:9, nv, slot))
         jcd(nv) = imtrace9(tmat1)
         ! Charge current part of jij
         tmat1 = 0.0d0
         tmat1 = dGdG(this, dmat1, dmat2, this%green%g01ij(1:9, 1:9, nv, slot), this%green%g01ji(1:9, 1:9, nv, slot))
         jcc(nv) = imtrace9(tmat1)
         ! Spin density part of jij
         tmat1 = 0.0d0
         tmat1 = dGdG(this, dmat1, dmat2, this%green%gx0ij(1:9, 1:9, nv, slot), this%green%gx0ji(1:9, 1:9, nv, slot))
         tmat2 = 0.0d0
         tmat2 = dGdG(this, dmat1, dmat2, this%green%gy0ij(1:9, 1:9, nv, slot), this%green%gy0ji(1:9, 1:9, nv, slot))
         tmat3 = 0.0d0
         tmat3 = dGdG(this, dmat1, dmat2, this%green%gz0ij(1:9, 1:9, nv, slot), this%green%gz0ji(1:9, 1:9, nv, slot))
         tmat4 = 0.0d0
         tmat4 = tmat1 + tmat2 + tmat3
         jsd(nv) = imtrace9(tmat4)
         ! Spin current part of jij
         tmat1 = 0.0d0
         tmat1 = dGdG(this, dmat1, dmat2, this%green%gx1ij(1:9, 1:9, nv, slot), this%green%gx1ji(1:9, 1:9, nv, slot))
         tmat2 = 0.0d0
         tmat2 = dGdG(this, dmat1, dmat2, this%green%gy1ij(1:9, 1:9, nv, slot), this%green%gy1ji(1:9, 1:9, nv, slot))
         tmat3 = 0.0d0
         tmat3 = dGdG(this, dmat1, dmat2, this%green%gz1ij(1:9, 1:9, nv, slot), this%green%gz1ji(1:9, 1:9, nv, slot))
         tmat4 = 0.0d0
         tmat4 = tmat1 + tmat2 + tmat3
         jsc(nv) = imtrace9(tmat4)
         ! Combining and calculating the first order and second order Jij(E)
         ! Second order J = CD - SD + CC - SC
         jtotso(nv) = jcd(nv) - jsd(nv) + jcc(nv) - jsc(nv)
         ! First order J = CD + SD - CC - SC
         jtotfo(nv) = jcd(nv) + jsd(nv) - jcc(nv) - jsc(nv)
         ! DMI
         do k = 1, 3
            select case (k)
            case (1) ! Dx
               tmat1 = 0.0d0; tmat2 = 0.0d0
               tmat1 = dGdG(this, dmat1, dmat2, this%green%g00ij(1:9, 1:9, nv, slot), this%green%gx1ji(1:9, 1:9, nv, slot))
               tmat2 = dGdG(this, dmat1, dmat2, this%green%g01ij(1:9, 1:9, nv, slot), this%green%gx0ji(1:9, 1:9, nv, slot))
            case (2) ! Dy
               tmat1 = 0.0d0; tmat2 = 0.0d0
               tmat1 = dGdG(this, dmat1, dmat2, this%green%g00ij(1:9, 1:9, nv, slot), this%green%gy1ji(1:9, 1:9, nv, slot))
               tmat2 = dGdG(this, dmat1, dmat2, this%green%g01ij(1:9, 1:9, nv, slot), this%green%gy0ji(1:9, 1:9, nv, slot))
            case (3) ! Dz
               tmat1 = 0.0d0; tmat2 = 0.0d0
               tmat1 = dGdG(this, dmat1, dmat2, this%green%g00ij(1:9, 1:9, nv, slot), this%green%gz1ji(1:9, 1:9, nv, slot))
               tmat2 = dGdG(this, dmat1, dmat2, this%green%g01ij(1:9, 1:9, nv, slot), this%green%gz0ji(1:9, 1:9, nv, slot))
            end select
            ! Spin current part of the DMI
            dsc(nv, k) = rtrace9(tmat1)
            ! Charge current part of the DMI
            dcc(nv, k) = rtrace9(tmat2)
            ! Combining and calculating the first and second order Dij(E)
            ! Second order D = SC + CC
            jjtotso(nv, k) = 2*(dsc(nv, k) + dcc(nv, k))
            ! First order D = SC - CC
            jjtotfo(nv, k) = 2*(dsc(nv, k) - dcc(nv, K))
         end do
         ! I Tensor
         do k = 1, 9
            select case (k)
            case (1) ! Axx
               tmat1 = 0.0d0; tmat2 = 0.0d0
               tmat1 = dGdG(this, dmat1, dmat2, this%green%gx0ij(1:9, 1:9, nv, slot), this%green%gx0ji(1:9, 1:9, nv, slot))
               tmat2 = dGdG(this, dmat1, dmat2, this%green%gx1ij(1:9, 1:9, nv, slot), this%green%gx1ji(1:9, 1:9, nv, slot))
               isd(nv, 1, 1) = imtrace9(tmat1)
               isc(nv, 1, 1) = imtrace9(tmat2)
            case (2) ! Axy
               tmat1 = 0.0d0; tmat2 = 0.0d0
               tmat1 = dGdG(this, dmat1, dmat2, this%green%gx0ij(1:9, 1:9, nv, slot), this%green%gy0ji(1:9, 1:9, nv, slot))
               tmat2 = dGdG(this, dmat1, dmat2, this%green%gx1ij(1:9, 1:9, nv, slot), this%green%gy1ji(1:9, 1:9, nv, slot))
               isd(nv, 1, 2) = imtrace9(tmat1)
               isc(nv, 1, 2) = imtrace9(tmat2)
            case (3) ! Axz
               tmat1 = 0.0d0; tmat2 = 0.0d0
               tmat1 = dGdG(this, dmat1, dmat2, this%green%gx0ij(1:9, 1:9, nv, slot), this%green%gz0ji(1:9, 1:9, nv, slot))
               tmat2 = dGdG(this, dmat1, dmat2, this%green%gx1ij(1:9, 1:9, nv, slot), this%green%gz1ji(1:9, 1:9, nv, slot))
               isd(nv, 1, 3) = imtrace9(tmat1)
               isc(nv, 1, 3) = imtrace9(tmat2)
            case (4) ! Ayx
               tmat1 = 0.0d0; tmat2 = 0.0d0
               tmat1 = dGdG(this, dmat1, dmat2, this%green%gy0ij(1:9, 1:9, nv, slot), this%green%gx0ji(1:9, 1:9, nv, slot))
               tmat2 = dGdG(this, dmat1, dmat2, this%green%gy1ij(1:9, 1:9, nv, slot), this%green%gx1ji(1:9, 1:9, nv, slot))
               isd(nv, 2, 1) = imtrace9(tmat1)
               isc(nv, 2, 1) = imtrace9(tmat2)
            case (5) ! Ayy
               tmat1 = 0.0d0; tmat2 = 0.0d0
               tmat1 = dGdG(this, dmat1, dmat2, this%green%gy0ij(1:9, 1:9, nv, slot), this%green%gy0ji(1:9, 1:9, nv, slot))
               tmat2 = dGdG(this, dmat1, dmat2, this%green%gy1ij(1:9, 1:9, nv, slot), this%green%gy1ji(1:9, 1:9, nv, slot))
               isd(nv, 2, 2) = imtrace9(tmat1)
               isc(nv, 2, 2) = imtrace9(tmat2)
            case (6) ! Ayz
               tmat1 = 0.0d0; tmat2 = 0.0d0
               tmat1 = dGdG(this, dmat1, dmat2, this%green%gy0ij(1:9, 1:9, nv, slot), this%green%gz0ji(1:9, 1:9, nv, slot))
               tmat2 = dGdG(this, dmat1, dmat2, this%green%gy1ij(1:9, 1:9, nv, slot), this%green%gz1ji(1:9, 1:9, nv, slot))
               isd(nv, 2, 3) = imtrace9(tmat1)
               isc(nv, 2, 3) = imtrace9(tmat2)
            case (7) ! Azx
               tmat1 = 0.0d0; tmat2 = 0.0d0
               tmat1 = dGdG(this, dmat1, dmat2, this%green%gz0ij(1:9, 1:9, nv, slot), this%green%gx0ji(1:9, 1:9, nv, slot))
               tmat2 = dGdG(this, dmat1, dmat2, this%green%gz1ij(1:9, 1:9, nv, slot), this%green%gx1ji(1:9, 1:9, nv, slot))
               isd(nv, 3, 1) = imtrace9(tmat1)
               isc(nv, 3, 1) = imtrace9(tmat2)
            case (8) ! Azy
               tmat1 = 0.0d0; tmat2 = 0.0d0
               tmat1 = dGdG(this, dmat1, dmat2, this%green%gz0ij(1:9, 1:9, nv, slot), this%green%gy0ji(1:9, 1:9, nv, slot))
               tmat2 = dGdG(this, dmat1, dmat2, this%green%gz1ij(1:9, 1:9, nv, slot), this%green%gy1ji(1:9, 1:9, nv, slot))
               isd(nv, 3, 2) = imtrace9(tmat1)
               isc(nv, 3, 2) = imtrace9(tmat2)
            case (9) ! Azz
               tmat1 = 0.0d0; tmat2 = 0.0d0
               tmat1 = dGdG(this, dmat1, dmat2, this%green%gz0ij(1:9, 1:9, nv, slot), this%green%gz0ji(1:9, 1:9, nv, slot))
               tmat2 = dGdG(this, dmat1, dmat2, this%green%gz1ij(1:9, 1:9, nv, slot), this%green%gz1ji(1:9, 1:9, nv, slot))
               isd(nv, 3, 3) = imtrace9(tmat1)
               isc(nv, 3, 3) = imtrace9(tmat2)
            end select
            ! Combining and calculating the first and second order A(E)
            ! Second order A = SD + SC
            itotso(nv, :, :) = isd(nv, :, :) + isc(nv, :, :)
            ! First order A = - SD + SC
            itotfo(nv, :, :) = - isd(nv, :, :) + isc(nv, :, :)
         end do
      end do
      ! Jij integration
      this%jij = 0.0d0
      call nodes%integrate(real(jtotso), this%jij)
      xcso(1) = this%jij*1.0d3/4.0d0/pi
      this%jij = 0.0d0
      call nodes%integrate(real(jtotfo), this%jij)
      xcfo(1) = this%jij*1.0d3/4.0d0/pi
      this%jijcd = 0.0d0 ; this%jijsd = 0.0d0 ; this%jijcc = 0.0d0 ; this%jijsc = 0.0d0 
      call nodes%integrate(real(jcd), this%jijcd)
      xcparts(1) = this%jijcd*1.0d3/4.0d0/pi
      call nodes%integrate(real(jsd), this%jijsd)
      xcparts(2) = this%jijsd*1.0d3/4.0d0/pi
      call nodes%integrate(real(jcc), this%jijcc)
      xcparts(3) = this%jijcc*1.0d3/4.0d0/pi
      call nodes%integrate(real(jsc), this%jijsc)
      xcparts(4) = this%jijsc*1.0d3/4.0d0/pi

      ! Dij integration
      this%dmi = 0.0d0
      do k = 1, 3
         y(:) = 0.0d0
         y(:) = real(jjtotso(:, k))
         call nodes%integrate(y(:), this%dmi(k))
      end do
      xcso(2:4) = this%dmi*1.0d3/4.0d0/pi
      this%dmi = 0.0d0
      do k = 1, 3
         y(:) = 0.0d0
         y(:) = real(jjtotfo(:, k))
         call nodes%integrate(y(:), this%dmi(k))
      end do
      xcfo(2:4) = this%dmi*1.0d3/4.0d0/pi
      this%dmicc = 0.0d0 ; this%dmisc = 0.0d0
      do k = 1, 3
         y(:) = 0.0d0
         y(:) = real(dsc(:, k))
         call nodes%integrate(y(:), this%dmisc(k))
         y(:) = 0.0d0
         y(:) = real(dcc(:, k))
         call nodes%integrate(y(:), this%dmicc(k))
      end do
      xcparts(5:7) = this%dmicc*2.0d3/4.0d0/pi
      xcparts(8:10) = this%dmisc*2.0d3/4.0d0/pi

      ! I tensor integration
      this%aij = 0.0d0
      do k = 1, 3
         do l = 1, 3
            y(:) = 0.0d0
            y(:) = real(itotso(:, k, l))
            call nodes%integrate(y(:), this%aij(k, l))
         end do
      end do
      xcso(5:13) = reshape(this%aij*1.0d3/4.0d0/pi, [9])
      this%aij = 0.0d0
      do k = 1, 3
         do l = 1, 3
            y(:) = 0.0d0
            y(:) = real(itotfo(:, k, l))
            call nodes%integrate(y(:), this%aij(k, l))
         end do
      end do
      xcfo(5:13) = reshape(this%aij*1.0d3/4.0d0/pi, [9])
      this%aijsd = 0.0d0 ; this%aijsc = 0.0d0  
      do k = 1, 3
         do l = 1, 3
            y(:) = 0.0d0
            y(:) = real(isd(:, k, l))
            call nodes%integrate(y(:), this%aijsd(k, l))
            y(:) = 0.0d0
            y(:) = real(isc(:, k, l))
            call nodes%integrate(y(:), this%aijsc(k, l))
         end do
      end do
      xcparts(11:19) = reshape(this%aijsd*1.0d3/4.0d0/pi, [9])
      xcparts(20:28) = reshape(this%aijsc*1.0d3/4.0d0/pi, [9])
   end subroutine exchange_twoindex_pair

   subroutine exchange_twoindex_finalize(this, T_comm_xcso, T_comm_xcfo, T_comm_xcparts)
      use mpi_mod
      class(exchange) :: this
      real(rp), intent(inout) :: T_comm_xcso(:, :), T_comm_xcfo(:, :), T_comm_xcparts(:, :)
      integer :: i, j, njij, njij_glob
      logical :: isopen

      inquire (unit=20, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file jijso.out: Unit 20 is already open', __FILE__, __LINE__)
      else
         open (unit=20, file='jijso.out')
      end if

      inquire (unit=25, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file jijfo.out: Unit 25 is already open', __FILE__, __LINE__)
      else
         open (unit=25, file='jijfo.out')
      end if

      inquire (unit=30, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file dijso.out: Unit 30 is already open', __FILE__, __LINE__)
      else
         open (unit=30, file='dijso.out')
      end if

      inquire (unit=35, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file dijfo.out: Unit 30 is already open', __FILE__, __LINE__)
      else
         open (unit=35, file='dijfo.out')
      end if

      inquire (unit=40, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file aijso.out: Unit 40 is already open', __FILE__, __LINE__)
      else
         open (unit=40, file='aijso.out')
      end if

      inquire (unit=45, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file aijfo.out: Unit 45 is already open', __FILE__, __LINE__)
      else
         open (unit=45, file='aijfo.out')
      end if

      inquire (unit=60, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file jtensso.out: Unit 60 is already open', __FILE__, __LINE__)
      else
         open (unit=60, file='jtensfo.out')
      end if

      inquire (unit=65, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file jtensfo.out: Unit 65 is already open', __FILE__, __LINE__)
      else
         open (unit=65, file='jtensso.out')
      end if

      inquire (unit=70, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file jijparts.out: Unit 70 is already open', __FILE__, __LINE__)
      else
         open (unit=70, file='jijparts.out')
      end if

      inquire (unit=75, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file dijparts.out: Unit 75 is already open', __FILE__, __LINE__)
      else
         open (unit=75, file='dijparts.out')
      end if

      inquire (unit=80, opened=isopen)
      if (isopen) then
         call g_logger%fatal('exchange%calculate_exchange, file aijparts.out: Unit 80 is already open', __FILE__, __LINE__)
      else
         open (unit=80, file='aijparts.out')
      end if


#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, T_comm_xcso, product(shape(T_comm_xcso)), &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, T_comm_xcfo, product(shape(T_comm_xcfo)), &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, T_comm_xcparts, product(shape(T_comm_xcparts)), &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

! Print exchange interactions only using MPI rank 0
      if (rank == 0) then
         do njij_glob = 1, this%lattice%njij
            !do njij_glob=start_atom, end_atom
            njij = g2l_map(njij_glob)
            i = this%lattice%ijpair(njij_glob, 1) ! Atom number in the clust file, atom i
            j = this%lattice%ijpair(njij_glob, 2) ! Atom number in the clust file, atom j

            !write(*,*) ´Atom´, i, ´coordinates:´, this%lattice%cr(:,i)
            !write(*,*) ´Atom´, j, ´coordinates:´, this%lattice%cr(:,j)

            ! Jij second order
            this%jij = T_comm_xcso(1, njij_glob)
            write (20, '(2i8,2x,3e20.11,2x,1es16.6,1x,f12.6)') &
      & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), this%jij, norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            ! Jij first order
            this%jij = T_comm_xcfo(1, njij_glob)
            write (25, '(2i8,2x,3e20.11,2x,1es16.6,1x,f12.6)') &
      & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), this%jij, norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            !  Jij parts
            this%jijcd = T_comm_xcparts(1, njij_glob) ; this%jijsd = T_comm_xcparts(2, njij_glob)    
            this%jijcc = T_comm_xcparts(3, njij_glob) ; this%jijsc = T_comm_xcparts(4, njij_glob) 
            write (70, '(2i8,2x,3e20.11,2x,4es16.6,1x,f12.6)') &
      & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), this%jijcd, this%jijsd, this%jijcc, this%jijsc, &
      & norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            ! Dij second order
            this%dmi = T_comm_xcso(2:4, njij_glob)
            write (30, '(2i8,2x,3e20.11,2x,3es16.6,1x,f12.6)') &
      & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), this%dmi, norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            ! Dij first order
            this%dmi = T_comm_xcfo(2:4, njij_glob)
            write (35, '(2i8,2x,3e20.11,2x,3es16.6,1x,f12.6)') &
      & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), this%dmi, norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            ! Dij parts
            this%dmicc = T_comm_xcparts(5:7, njij_glob) ; this%dmisc = T_comm_xcparts(8:10, njij_glob)
            write (75, '(2i8,2x,3e20.11,2x,6es16.6,1x,f12.6)') &
      & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), this%dmicc, this%dmisc, norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            ! Aij second order
            this%aij = reshape(T_comm_xcso(5:13, njij_glob), [3, 3])
            write (40, '(2i8,2x,3e20.11,2x,9es16.6,1x,f12.6)') &
      & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), this%aij, norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            ! Aij first order
            this%aij = reshape(T_comm_xcfo(5:13, njij_glob), [3, 3])
            write (45, '(2i8,2x,3e20.11,2x,9es16.6,1x,f12.6)') &
      & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), this%aij, norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
            ! Aij parts
            this%aijsd = reshape(T_comm_xcparts(11:19, njij_glob), [3, 3]); this%aijsc = reshape(T_comm_xcparts(20:28, njij_glob), [3, 3])
            write (80, '(2i8,2x,3e20.11,2x,18es16.6,1x,f12.6)') &
      & this%lattice%iz(i), this%lattice%iz(j), this%lattice%cr(:, j) - this%lattice%cr(:, i), this%aijsd, this%aijsc, norm2(this%lattice%cr(:, i) - this%lattice%cr(:, j))
         end do
      end if

      close (20)
      close (25)
      close (30)
      close (35)
      close (40)
      close (45)
      close (60)
      close (65)
      close (70)
      close (75) 
      close (80)

#ifdef USE_MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

   end subroutine exchange_twoindex_finalize

   !> Jij, Dij, Aij kernel of pair (i, j) read from Green-function slot `slot`, in meV, into xc(1:13).
   subroutine exchange_pair(this, nodes, slot, i, j, xc)
      class(exchange) :: this
      type(energy_nodes), intent(in) :: nodes
      integer, intent(in) :: slot, i, j
      real(rp), intent(out) :: xc(13)
      complex(rp), dimension(9, 9) :: dmat1, dmat2
      complex(rp), dimension(size(nodes%z)) :: jtot
      complex(rp), dimension(size(nodes%z), 3) :: jjtot
      complex(rp), dimension(size(nodes%z), 3, 3) :: itot
      real(rp), dimension(size(nodes%z)) :: y
      integer :: nv, k, l
      complex(rp), dimension(9, 9) :: jmat
      complex(rp), dimension(9, 9, 3) :: dmmats
      complex(rp), dimension(9, 9, 3, 3) :: amats

      jtot = 0.0d0; jjtot = 0.0d0
      if (this%use_onsite_delta) then
         dmat1 = onsite_exchange_delta(this, i)
         dmat2 = onsite_exchange_delta(this, j)
      end if
      do nv = 1, this%en%channels_ldos + 10
         if (.not. this%use_onsite_delta) then
            call this%symbolic_atom(this%lattice%iz(i))%d_matrix(dmat1, real(nodes%z(nv)))
            call this%symbolic_atom(this%lattice%iz(j))%d_matrix(dmat2, real(nodes%z(nv)))
         end if
         !J
         jmat = dGdG_Jnc(this, dmat1, dmat2, slot, nv)
         jtot(nv) = imtrace9(jmat)

         ! DMI
         dmmats = dGdG_Dnc(this, dmat1, dmat2, slot, nv)
         do k = 1, 3
            jjtot(nv, k) = rtrace9(dmmats(:, :, k))
         end do

         ! I Tensor
         amats = dGdG_Anc(this, dmat1, dmat2, slot, nv)
         do k = 1, 3
            do l = 1, 3
               itot(nv, k, l) = imtrace9(amats(:, :, k, l))
            end do
         end do
      end do
      ! Jij integration
      this%jij = 0.0d0
      call nodes%integrate(real(jtot), this%jij)
      xc(1) = this%jij*1.0d3/4.0d0/pi

      ! Dij integration
      this%dmi = 0.0d0
      do k = 1, 3
         y(:) = 0.0d0
         y(:) = real(jjtot(:, k))
         call nodes%integrate(y(:), this%dmi(k))
      end do
      xc(2:4) = this%dmi*1.0d3/4.0d0/pi

      ! I tensor integration
      this%aij = 0.0d0
      do k = 1, 3
         do l = 1, 3
            y(:) = 0.0d0
            y(:) = real(itot(:, k, l))
            call nodes%integrate(y(:), this%aij(k, l))
         end do
      end do
      xc(5:13) = reshape(this%aij*1.0d3/4.0d0/pi, [9])

   end subroutine exchange_pair

   subroutine exchange_finalize(this, T_comm_xc)
      use mpi_mod
      use Parameters
      class(exchange) :: this
      real(rp), intent(inout) :: T_comm_xc(:, :)
      real(rp), dimension(3, 3) :: jtens
      integer :: i, j, njij, njij_glob
      logical :: isopen

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
            write (99, *) 'null', (this%lattice%cr(:, j) + this%lattice%cr(:, i))/2, this%dmi/MAX(norm2(this%dmi), TINY(1.0_rp))
            ! Compute and write the jtens.out file
            ! J on the diagonal
            jtens = 0.0d0
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
         if (allocated(this%lattice%pair_rep)) call write_pair_stars(this%lattice, T_comm_xc)
      end if

      close (20)
      close (30)
      close (40)
      close (60)

#ifdef USE_MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

   end subroutine exchange_finalize

   !> Full star expansion of the pair_unique list, one row per member pair of every star: J, D' = Rt D, A' = Rt A Rt^T
   !> of the representative, Rt = det(R) R with R the Cartesian rotation that maps the representative onto the member.
   subroutine write_pair_stars(lat, xc)
      type(lattice), intent(in) :: lat
      real(rp), intent(in) :: xc(:, :)
      real(rp) :: r(3, 3), a(3, 3), d(3)
      integer :: f, i, j, p, unit

      open (newunit=unit, file='exchange_stars.out', action='write', status='replace')
      write (unit, '(a)') '# i j type_i type_j r_ij(3, cluster coordinates) |r_ij| J(mRy) D(3, mRy) A(9, mRy, as aij.out)'
      do f = 1, size(lat%pair_rep)
         i = lat%pair_full(f, 1)
         j = lat%pair_full(f, 2)
         p = lat%pair_rep(f)
         r = lat%pair_rot(:, :, f)
         ! Without spin-orbit coupling the spins are decoupled from the lattice: nothing to rotate
         if (lat%control%nsp == 1) r = reshape([1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp], [3, 3])
         if (r(1, 1)*(r(2, 2)*r(3, 3) - r(2, 3)*r(3, 2)) - r(1, 2)*(r(2, 1)*r(3, 3) - r(2, 3)*r(3, 1)) &
             + r(1, 3)*(r(2, 1)*r(3, 2) - r(2, 2)*r(3, 1)) < 0.0_rp) r = -r
         d = matmul(r, xc(2:4, p))
         a = matmul(r, matmul(reshape(xc(5:13, p), [3, 3]), transpose(r)))
         write (unit, '(4i8,2x,3f12.6,1x,f12.6,2x,f12.6,2x,3f12.6,2x,9f12.6)') i, j, lat%iz(i), lat%iz(j), &
            lat%cr(:, j) - lat%cr(:, i), norm2(lat%cr(:, j) - lat%cr(:, i)), xc(1, p), d, a
      end do
      close (unit)
   end subroutine write_pair_stars

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
      T_comm_xc = 0.0_rp

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
         jtot = 0.0d0; jjtot = 0.0d0
         dmat1 = 0.0d0; dmat2 = 0.0d0; tmat1 = 0.0d0; tmat2 = 0.0d0; tmat3 = 0.0d0; smat = 0.0d0

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
         this%jij = 0.0d0
      !!! do nv=1,this%en%channels_ldos+10
      !!!   call simpson_f(this%jij,this%en%ene,this%en%ene(nv),this%en%nv1,real(jtot),.true.,.false.,0.0d0)
      !!!   write(902,*) this%en%ene(nv)-this%en%fermi, this%jij*1.0d3/4.0d0/pi
      !!! end do
         call simpson_f(this%jij, this%en%ene, this%en%ene(this%en%enpt + 1), this%en%nv1, real(jtot), .true., .false., this%en%temperature)
         T_comm_xc(1, njij_glob) = this%jij*1.0d3/4.0d0/pi

         this%dmi = 0.0d0
         do k = 1, 3
            y(:) = 0.0d0
            y(:) = real(jjtot(:, k))
            call simpson_f(this%dmi(k), this%en%ene, this%en%ene(this%en%enpt + 1), this%en%nv1, y(:), .true., .false., this%en%temperature)
         end do
         T_comm_xc(2:4, njij_glob) = this%dmi*1.0d3/4.0d0/pi

         ! I tensor integration
         this%aij = 0.0d0
         do k = 1, 3
            do l = 1, 3
               y(:) = 0.0d0
               y(:) = real(itot(:, k, l))
               call simpson_f(this%aij(k, l), this%en%ene, this%en%ene(this%en%enpt + 1), this%en%nv1, y(:), .true., .false., this%en%temperature)
            end do
         end do
         T_comm_xc(5:13, njij_glob) = reshape(this%aij*1.0d3/4.0d0/pi, [9])

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
      type(energy_nodes) :: nodes
      integer :: i, j, njij, njij_glob
      real(rp), dimension(:, :), allocatable :: T_comm_xc
      logical :: isopen

      ! Communications array for MPI
      ! row = [ jij dij aij]
      allocate (T_comm_xc(13, this%lattice%njij))
      T_comm_xc = 0.0_rp

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
      nodes = contour_nodes(this)

      do njij_glob = start_atom, end_atom
         njij = g2l_map(njij_glob)
         i = this%lattice%ijpair(njij_glob, 1) ! Atom number in the clust file, atom i
         j = this%lattice%ijpair(njij_glob, 2) ! Atom number in the clust file, atom j
         call exchange_contour_pair(this, nodes, njij, i, j, T_comm_xc(:, njij_glob))
      end do

      call exchange_contour_finalize(this, T_comm_xc)
      deallocate (T_comm_xc)
   end subroutine calculate_exchange_gauss_legendre

   !> Gauss-Legendre contour kernel: J, D, A of pair (i, j) at local index njij, in meV, into xc(1:13).
   subroutine exchange_contour_pair(this, nodes, njij, i, j, xc)
      class(exchange) :: this
      type(energy_nodes), intent(in) :: nodes
      integer, intent(in) :: njij, i, j
      real(rp), intent(out) :: xc(13)
      complex(rp), dimension(9, 9) :: dmat1, dmat2
      complex(rp), dimension(64) :: jtot
      complex(rp), dimension(64, 3) :: jjtot
      complex(rp), dimension(64, 3, 3) :: itot
      integer :: nv, k, l
      complex(rp), dimension(9, 9) :: jmat
      complex(rp), dimension(9, 9, 3) :: dmmats
      complex(rp), dimension(9, 9, 3, 3) :: amats

      jtot = 0.0d0; jjtot = 0.0d0
      dmat1 = real(this%hamiltonian%ee(1:9, 1:9, 1, this%lattice%iz(i)) - this%hamiltonian%ee(10:18, 10:18, 1, this%lattice%iz(i)))
      dmat2 = real(this%hamiltonian%ee(1:9, 1:9, 1, this%lattice%iz(j)) - this%hamiltonian%ee(10:18, 10:18, 1, this%lattice%iz(j)))

      jtot = 0.0d0

      do nv = 1, 64
         !J
         jmat = dGdG_Jnc(this, dmat1, dmat2, njij, nv)
         jmat = (jmat*nodes%w(nv))/(nodes%x(nv)*nodes%x(nv))
         jtot(nv) = rtrace9(jmat)

         ! DMI
         dmmats = dGdG_Dnc(this, dmat1, dmat2, njij, nv)
         dmmats = (dmmats*nodes%w(nv))/(nodes%x(nv)*nodes%x(nv))
         do k = 1, 3
            jjtot(nv, k) = imtrace9(dmmats(:, :, k))
         end do

         ! I Tensor
         amats = dGdG_Anc(this, dmat1, dmat2, njij, nv)
         amats = (amats*nodes%w(nv))/(nodes%x(nv)*nodes%x(nv))
         do k = 1, 3
            do l = 1, 3
               itot(nv, k, l) = rtrace9(amats(:, :, k, l))
            end do
         end do

      end do

      ! Jij summation
      call nodes%integrate(real(jtot), this%jij)
      this%jij = -this%jij
      xc(1) = this%jij*1.0d3/4.0d0/pi

      ! DMI
      do k = 1, 3
         call nodes%integrate(real(jjtot(:, k)), this%dmi(k))
      end do
      xc(2:4) = this%dmi*1.0d3/4.0d0/pi

      ! I tensor integration
      this%aij = 0.0d0
      do k = 1, 3
         do l = 1, 3
            call nodes%integrate(real(itot(:, k, l)), this%aij(k, l))
            this%aij(k, l) = -this%aij(k, l)
         end do
      end do
      xc(5:13) = reshape(this%aij*1.0d3/4.0d0/pi, [9])

   end subroutine exchange_contour_pair

   subroutine exchange_contour_finalize(this, T_comm_xc)
      use mpi_mod
      class(exchange) :: this
      real(rp), intent(inout) :: T_comm_xc(:, :)
      integer :: i, j, njij_glob

#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, T_comm_xc, product(shape(T_comm_xc)), &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      if (rank == 0) then
         do njij_glob = 1, this%lattice%njij
            !do njij_glob=start_atom, end_atom
            !njij = g2l_map(njij_glob)
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
            write (99, *) 'null', (this%lattice%cr(:, j) + this%lattice%cr(:, i))/2, this%dmi/MAX(norm2(this%dmi), TINY(1.0_rp))
         end do
      end if
      close (20)
      close (30)
      close (40)
#ifdef USE_MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
   end subroutine exchange_contour_finalize

   !> Integrate the real integrand y sampled on the nodes: Simpson with the Fermi function
   !> on the real mesh, plain sum of the already weighted values on the contour.
   subroutine nodes_integrate(this, y, res)
      class(energy_nodes), intent(in) :: this
      real(rp), intent(in) :: y(:)
      real(rp), intent(out) :: res

      if (this%contour) then
         res = sum(y)
      else
         call simpson_f(res, real(this%z), this%ef, this%nv1, y, .true., .false., this%temperature)
      end if
   end subroutine nodes_integrate

   !> Real-energy mesh of en%ene (Simpson, Fermi function at ene(enpt + 1)).
   function real_mesh_nodes(this) result(nodes)
      class(exchange), intent(in) :: this
      type(energy_nodes) :: nodes

      nodes%z = cmplx(this%en%ene(1:this%en%channels_ldos + 10), 0.0_rp, rp)
      nodes%ef = this%en%ene(this%en%enpt + 1)
      nodes%nv1 = this%en%nv1
      nodes%temperature = this%en%temperature
   end function real_mesh_nodes

   !> 64 Gauss-Legendre nodes of the contour z = E(fermi_point) + i (1 - x)/x, as in green_block.
   function contour_nodes(this) result(nodes)
      class(exchange), intent(in) :: this
      type(energy_nodes) :: nodes
      integer :: i, fermi_point

      allocate (nodes%x(64), nodes%w(64))
      call gauss_legendre(64, 0.0_rp, 1.0_rp, nodes%x, nodes%w)
      fermi_point = 1
      do i = 1, this%en%channels_ldos + 10
         if ((this%en%ene(i) - this%en%fermi) .le. 0.000001d0) fermi_point = i
      end do
      nodes%z = cmplx(this%en%ene(fermi_point), (1.0_rp - nodes%x)/nodes%x, rp)
      nodes%contour = .true.
   end function contour_nodes

end module exchange_mod
