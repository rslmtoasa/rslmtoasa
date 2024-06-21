!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Charge
!
!> @author
!> Angela Klautau
!> Ramon Cardias
!> Lucas P. Campagna
!> Ivan P. Miranda
!> S. Frota-Pessôa
!> Pascoal R. Peduto
!> Anders Bergman
!> S. B. Legoas
!> H. M. Petrilli
!
! DESCRIPTION:
!> Module to handle system´s potential charge and energy
!------------------------------------------------------------------------------

module charge_mod
   use mpi_mod
   use lattice_mod
   use symbolic_atom_mod, only: symbolic_atom
   use string_mod, only: sl, fmt
   use math_mod, only: pi, sqrt_pi, cross_product, ang2au, distance, angle, pos, erodrigues, normalize
   use precision_mod, only: rp
   use namelist_generator_mod, only: namelist_generator
   use logger_mod, only: g_logger
   use string_mod
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   !> Module´s main structure
   type, public :: charge
      !> Lattice
      class(lattice), pointer :: lattice
      !> Symbolic atom
      class(symbolic_atom), dimension(:), pointer :: symbolic_atom

      !> Bulkmat variables
      real(rp), dimension(:), allocatable :: w
      !> Lattice distortion paramaters for bulkmat. Default is 1 for no distortion
      real(rp) :: gx, gy, gz, gt

      !> Surfmat variables
      !> Wigner Seitz radius per layer considered to receive charge. It has a
      !dimension of nbas=17 as default, but if the layer is big with different
      !elements etc, a higher number has to be considered
      real(rp), dimension(:), allocatable :: wssurf
      !> Monopole contribution to the eletrostatic potential.
      real(rp), dimension(:, :), allocatable :: dss
      !> Dipole contribution to the eletrostatic potential.
      real(rp), dimension(:, :), allocatable :: dsz
      !> High order contributions to the eletrostatic potential.
      real(rp), dimension(:, :), allocatable :: ds3z2, dsx2y2, dsxy, dzz, dz3z2, am, bm, pm
      !> Translational 2d vectors on real space
      real(rp), dimension(:), allocatable :: bsx, bsy, bsz
      !> Translational 2d vectors on reciprocal space
      real(rp), dimension(:), allocatable :: bkx, bky, bkz
      !> Position of atoms in the 3D cell
      real(rp), dimension(:), allocatable :: bssx, bssy, bssz
      !> Translational 2d vectors in reciprocal space
      real(rp), dimension(:), allocatable :: new_x, new_y, new_z
      !> New rotated coordinate system in the reciprocal space (z is perpendicular to plane)
      real(rp), dimension(:), allocatable :: new_kx, new_ky, new_kz
      !> Position of atoms in the 3D cell
      real(rp), dimension(:), allocatable :: qx3, qy3, qz3
      !> NQ basis vectors of the 2D cell
      real(rp), dimension(:), allocatable :: qx, qy, qz
      !> Auxiliar vectors real and reciprocal space
      real(rp), dimension(:), allocatable :: asx, asy, asz, akx, aky, akz, dr, dg
      !> Lattice parameter a, b and c
      real(rp) :: a, b, c
      !> Setup parameters
      real(rp) :: amax, bmax, alamda, rmax, gmax, ar2d
      !> Atomic Wigner-Seitz radius
      real(rp) :: sws
      !> Volume (a.b)
      real(rp) :: vol
      !> Number of atoms in the 3D cell
      integer :: nq3
      !> Numbers for surfmat
      integer :: nr0, numr, numg, numvr, numvg
      !> Charge neutrality
      real(rp), dimension(:), allocatable :: dq
      !> Charge transfer
      real(rp) :: cht, trq
      !> Charge transfer for impurity calculation
      real(rp), dimension(:), allocatable :: bulk_charge
      !> Mix pot
      real(rp) :: vmix

      !> Impmat variables
      !> Impurity eletrostatic potential
      real(rp), dimension(:, :), allocatable :: amad
      !> Wigner Seitz radius per atom shell considered to receive charge
      real(rp), dimension(:), allocatable :: wsimp
   contains
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: bulkpot
      procedure :: surfpot
      procedure :: bulkmat
      procedure :: surfmat
      procedure :: build_alelay
      procedure :: get_charge_transf
      procedure :: impmad
      procedure :: imppot
      procedure :: print_state
      procedure :: print_state_full
      procedure :: print_state_formatted
      final :: destructor
   end type charge

   interface charge
      procedure :: constructor
   end interface charge

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] lattice_obj Pointer to system´s lattice
   !> @return type(calculation)
   !---------------------------------------------------------------------------
   function constructor(lattice_obj) result(obj)
      type(charge) :: obj
      type(lattice), target, intent(in) :: lattice_obj

      obj%lattice => lattice_obj
      obj%symbolic_atom => lattice_obj%symbolic_atoms

      call obj%restore_to_default()
      call obj%build_from_file()
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   subroutine destructor(this)
      type(charge) :: this
#ifdef USE_SAFE_ALLOC
      if (allocated(this%wsimp)) call g_safe_alloc%deallocate('charge.wsimp', this%wsimp)
      if (allocated(this%wssurf)) call g_safe_alloc%deallocate('charge.wssurf', this%wssurf)
      if (allocated(this%dss)) call g_safe_alloc%deallocate('charge.dss', this%dss)
      if (allocated(this%dsz)) call g_safe_alloc%deallocate('charge.dsz', this%dsz)
      if (allocated(this%ds3z2)) call g_safe_alloc%deallocate('charge.ds3z2', this%ds3z2)
      if (allocated(this%dsx2y2)) call g_safe_alloc%deallocate('charge.dsx2y2', this%dsx2y2)
      if (allocated(this%dsxy)) call g_safe_alloc%deallocate('charge.dsxy', this%dsxy)
      if (allocated(this%dzz)) call g_safe_alloc%deallocate('charge.dzz', this%dzz)
      if (allocated(this%dz3z2)) call g_safe_alloc%deallocate('charge.dz3z2', this%dz3z2)
      if (allocated(this%am)) call g_safe_alloc%deallocate('charge.am', this%am)
      if (allocated(this%bm)) call g_safe_alloc%deallocate('charge.bm', this%bm)
      if (allocated(this%pm)) call g_safe_alloc%deallocate('charge.pm', this%pm)
      if (allocated(this%bsx)) call g_safe_alloc%deallocate('charge.bsx', this%bsx)
      if (allocated(this%bsy)) call g_safe_alloc%deallocate('charge.bsy', this%bsy)
      if (allocated(this%bsz)) call g_safe_alloc%deallocate('charge.bsz', this%bsz)
      if (allocated(this%bkx)) call g_safe_alloc%deallocate('charge.bkx', this%bkx)
      if (allocated(this%bky)) call g_safe_alloc%deallocate('charge.bky', this%bky)
      if (allocated(this%bkz)) call g_safe_alloc%deallocate('charge.bkz', this%bkz)
      if (allocated(this%bssx)) call g_safe_alloc%deallocate('charge.bssx', this%bssx)
      if (allocated(this%bssy)) call g_safe_alloc%deallocate('charge.bssy', this%bssy)
      if (allocated(this%bssz)) call g_safe_alloc%deallocate('charge.bssz', this%bssz)
      if (allocated(this%new_x)) call g_safe_alloc%deallocate('charge.new_x', this%new_x)
      if (allocated(this%new_y)) call g_safe_alloc%deallocate('charge.new_y', this%new_y)
      if (allocated(this%new_z)) call g_safe_alloc%deallocate('charge.new_z', this%new_z)
      if (allocated(this%new_kx)) call g_safe_alloc%deallocate('charge.new_kx', this%new_kx)
      if (allocated(this%new_ky)) call g_safe_alloc%deallocate('charge.new_ky', this%new_ky)
      if (allocated(this%new_kz)) call g_safe_alloc%deallocate('charge.new_kz', this%new_kz)
      if (allocated(this%qx3)) call g_safe_alloc%deallocate('charge.qx3', this%qx3)
      if (allocated(this%qy3)) call g_safe_alloc%deallocate('charge.qy3', this%qy3)
      if (allocated(this%qz3)) call g_safe_alloc%deallocate('charge.qz3', this%qz3)
      if (allocated(this%qx)) call g_safe_alloc%deallocate('charge.qx', this%qx)
      if (allocated(this%qy)) call g_safe_alloc%deallocate('charge.qy', this%qy)
      if (allocated(this%qz)) call g_safe_alloc%deallocate('charge.qz', this%qz)
      if (allocated(this%asx)) call g_safe_alloc%deallocate('charge.asx', this%asx)
      if (allocated(this%asy)) call g_safe_alloc%deallocate('charge.asy', this%asy)
      if (allocated(this%asz)) call g_safe_alloc%deallocate('charge.asz', this%asz)
      if (allocated(this%akx)) call g_safe_alloc%deallocate('charge.akx', this%akx)
      if (allocated(this%aky)) call g_safe_alloc%deallocate('charge.aky', this%aky)
      if (allocated(this%akz)) call g_safe_alloc%deallocate('charge.akz', this%akz)
      if (allocated(this%dr)) call g_safe_alloc%deallocate('charge.dr', this%dr)
      if (allocated(this%dg)) call g_safe_alloc%deallocate('charge.dg', this%dg)
      if (allocated(this%dq)) call g_safe_alloc%deallocate('charge.dq', this%dq)
#else
      if (allocated(this%wsimp)) deallocate (this%wsimp)
      if (allocated(this%wssurf)) deallocate (this%wssurf)
      if (allocated(this%dss)) deallocate (this%dss)
      if (allocated(this%dsz)) deallocate (this%dsz)
      if (allocated(this%ds3z2)) deallocate (this%ds3z2)
      if (allocated(this%dsx2y2)) deallocate (this%dsx2y2)
      if (allocated(this%dsxy)) deallocate (this%dsxy)
      if (allocated(this%dzz)) deallocate (this%dzz)
      if (allocated(this%dz3z2)) deallocate (this%dz3z2)
      if (allocated(this%am)) deallocate (this%am)
      if (allocated(this%bm)) deallocate (this%bm)
      if (allocated(this%pm)) deallocate (this%pm)
      if (allocated(this%bsx)) deallocate (this%bsx)
      if (allocated(this%bsy)) deallocate (this%bsy)
      if (allocated(this%bsz)) deallocate (this%bsz)
      if (allocated(this%bkx)) deallocate (this%bkx)
      if (allocated(this%bky)) deallocate (this%bky)
      if (allocated(this%bkz)) deallocate (this%bkz)
      if (allocated(this%bssx)) deallocate (this%bssx)
      if (allocated(this%bssy)) deallocate (this%bssy)
      if (allocated(this%bssz)) deallocate (this%bssz)
      if (allocated(this%new_x)) deallocate (this%new_x)
      if (allocated(this%new_y)) deallocate (this%new_y)
      if (allocated(this%new_z)) deallocate (this%new_z)
      if (allocated(this%new_kx)) deallocate (this%new_kx)
      if (allocated(this%new_ky)) deallocate (this%new_ky)
      if (allocated(this%new_kz)) deallocate (this%new_kz)
      if (allocated(this%qx3)) deallocate (this%qx3)
      if (allocated(this%qy3)) deallocate (this%qy3)
      if (allocated(this%qz3)) deallocate (this%qz3)
      if (allocated(this%qx)) deallocate (this%qx)
      if (allocated(this%qy)) deallocate (this%qy)
      if (allocated(this%qz)) deallocate (this%qz)
      if (allocated(this%asx)) deallocate (this%asx)
      if (allocated(this%asy)) deallocate (this%asy)
      if (allocated(this%asz)) deallocate (this%asz)
      if (allocated(this%akx)) deallocate (this%akx)
      if (allocated(this%aky)) deallocate (this%aky)
      if (allocated(this%akz)) deallocate (this%akz)
      if (allocated(this%dr)) deallocate (this%dr)
      if (allocated(this%dg)) deallocate (this%dg)
      if (allocated(this%dq)) deallocate (this%dq)
#endif
   end subroutine destructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Read parameters from input file
   !
   !> @param[in] fname Namelist file
   !---------------------------------------------------------------------------
   subroutine build_from_file(this, fname)
      class(charge), intent(inout) :: this
      character(len=*), intent(in), optional :: fname
      character(len=sl) :: fname_
      ! Local variables
      integer :: iostatus, funit

      include 'include_codes/namelists/charge.f90'

      if (present(fname)) then
         fname_ = fname
         this%lattice%control%fname = fname
      else
         fname_ = this%lattice%control%fname
      end if

      ! Save previous values
      gx = this%gx
      gy = this%gy
      gz = this%gz
      gt = this%gt
      vmix = this%vmix

      if (size(this%wssurf) .ne. this%lattice%nbas) then
#ifdef USE_SAFE_ALLOC
         call g_safe_alloc%deallocate('charge.wssurf', this%wssurf)
         call g_safe_alloc%allocate('charge.wssurf', this%wssurf, (/this%lattice%nbas/))
#else
         deallocate (this%wssurf)
         allocate (this%wssurf(this%lattice%nbas))
#endif
      end if

      call move_alloc(this%wssurf, wssurf)

      open (newunit=funit, file=fname_, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//fmt('A', fname_)//' not found', __FILE__, __LINE__)
      end if

      read (funit, nml=charge, iostat=iostatus)
      close (funit)

      if (size(wssurf) .ne. this%lattice%nbas) then
         deallocate (wssurf)
         allocate (wssurf(this%lattice%nbas))
      end if

      this%gx = gx
      this%gy = gy
      this%gz = gz
      this%gt = gt
      this%vmix = vmix
      call move_alloc(wssurf, this%wssurf)
   end subroutine build_from_file

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this)
      class(charge), intent(inout) :: this

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('charge.wssurf', this%wssurf, this%lattice%nbas)
      call g_safe_alloc%allocate('charge.dq', this%dq, this%lattice%nrec)
      call g_safe_alloc%allocate('charge.bulk_charge', this%bulk_charge, this%lattice%nbulk)
#else
      allocate (this%wssurf(this%lattice%nbas), this%dq(this%lattice%nrec))
      allocate (this%bulk_charge(this%lattice%nbulk))
#endif

      this%dq(:) = 0.0d0
      this%vmix = 1.0d0
      ! For surfmat
      this%wssurf(:) = this%lattice%wav*ang2au

      ! For bulkmat
      this%gx = 1.0d0
      this%gy = 1.0d0
      this%gz = 1.0d0
      this%gt = 1.0d0
   end subroutine

   subroutine bulkpot(this)
      class(charge), intent(inout) :: this
      ! Local variables
      integer :: I, IBAS, ICLAS, IM, JBAS, j, i_all, i_stat
      real(rp) :: DIF, VADD, VERR, VMADI
      real(rp), dimension(this%lattice%nrec) :: TDQ, RMAX
      real(rp), dimension(:, :), allocatable :: AMAD
      real(rp), dimension(this%lattice%nrec) :: VMAD0
      integer, parameter :: npot = 5000

      do i = 1, this%lattice%nrec
         rmax(i) = this%symbolic_atom(this%lattice%nbulk + i)%potential%ws_r
      end do
      !write (*, *) " NBAS=", NBAS, " NCLAS=", NCLAS
      !print *, ´ Sum of electrons from recursion: ´, sum(DQ)
      do im = 1, this%lattice%nrec
         VMAD0(IM) = this%symbolic_atom(this%lattice%nbulk + im)%potential%vmad
      end do
      if (this%lattice%NBAS > npot) then
         write (*, *) "# DE ATOMOS NA CELA FERE AS DIMENSOES ATRIBUIDAS"
         write (*, *) "             JOB ABORTADO                       "
         write (*, *) "NBAS=", this%lattice%NBAS, "DIFF NDIM=", NPOT, "IN MAD"
         !IFAIL = 1
         return
      end if
      DIF = 0.d0
      tdq(:) = 0.0d0
      tdq(:) = this%dq(:)
      dif = sum(tdq(:))
      !write (*, *) "DIF IN MAD PROGRAM = ", DIF
      if (ABS(DIF) > 5d-01) then
         !IFAIL = 1
         write (6, *) &
            "******NEUTRALIDADE DE CARGA NAO OBEDECIDA - 5E-04", &
            " PROGRAMA ABORTADO"
         return
      end if
      allocate (AMAD(this%lattice%nbas, this%lattice%nbas), stat=i_stat)
      open (15, file="mad.mat", form="unformatted", STATUS="OLD")
      do i = 1, this%lattice%nbas
         read (15) (AMAD(i, j), j=1, this%lattice%nbas)
      end do
      close (15)
      do IBAS = 1, this%lattice%NBAS
         VMADI = 0.d0
         do JBAS = 1, this%lattice%NBAS
            VMADI = VMADI + 2.d0*AMAD(JBAS, IBAS)*TDQ(this%lattice%iz(JBAS))
         end do
         this%symbolic_atom(this%lattice%nbulk + this%lattice%iz(IBAS))%potential%VMAD = VMADI
      end do
      !write (9, 10001)
      do ICLAS = 1, this%lattice%nrec
         VADD = 2.d0*TDQ(ICLAS)/RMAX(ICLAS)
         !write (9, 10000) ICLAS, TDQ(ICLAS), VMAD(ICLAS), VMAD(ICLAS)+VADD
         this%symbolic_atom(this%lattice%nbulk + iclas)%potential%VMAD = &
            this%symbolic_atom(this%lattice%nbulk + iclas)%potential%VMAD + VADD
         !   write (42, *) VMAD(ICLAS), VMAD0(ICLAS)
         this%symbolic_atom(this%lattice%nbulk + iclas)%potential%VMAD = &
            this%symbolic_atom(this%lattice%nbulk + iclas)%potential%VMAD*this%VMIX + VMAD0(ICLAS)*(1.0 - this%VMIX)
         VERR = this%symbolic_atom(this%lattice%nbulk + iclas)%potential%VMAD - VMAD0(ICLAS)
         !write (*, *) this%symbolic_atom(this%lattice%nbulk+iclas)%potential%VMAD, VERR
         !   write (42, *)
      end do
   end subroutine bulkpot

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Gets the charge transfer from the bulk calculation
   !---------------------------------------------------------------------------
   subroutine get_charge_transf(this)
      class(charge), intent(inout) :: this
      ! Local variables
      integer :: ia ! Atomic types

      do ia = 1, this%lattice%nbulk
         this%bulk_charge(ia) = sum(this%symbolic_atom(ia)%potential%ql(1, :, 1)) + sum(this%symbolic_atom(ia)%potential%ql(1, :, 2)) - this%symbolic_atom(ia)%element%valence ! occupation - valence
      end do
   end subroutine get_charge_transf

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the local madelung potential for the impurity mode
   !---------------------------------------------------------------------------
   subroutine imppot(this)
      implicit none
      class(charge), intent(inout) :: this
      ! Local variables
      integer :: I, IBAS, ICLAS, IEX, II, INEQ, J, JBAS, JJ, KK, LB, LC, N, NBULK, NRX, NSUM
      real(rp) :: DIF, SS, VERR
      integer, dimension(this%lattice%nbas + 1) :: IBU
      real(rp), dimension(this%lattice%nbas + 1) :: BUC, BUV, CTB, TDQ, VEB
      real(rp), dimension(this%lattice%nbas + 1) :: TOTQ
      real(rp), dimension(this%lattice%nrec + 1) :: VMAD0

      ! Initialization
      ineq = this%lattice%nbas + 1
      iex = this%lattice%nrec + 1
      j = 0
      do ii = 1, this%lattice%nbas
         j = j + 1
      end do
      nrx = j

      ! Initializing charge variables
      do j = iex, ineq
         tdq(j) = 0.0
         ctb(j) = 0.0
      end do
      dif = 0.d0

      ! Calculating charge transfer
      do iclas = 1, this%lattice%nrec
         tdq(iclas) = 0.d0
         tdq(iclas) = this%dq(iclas)
         tdq(iclas) = tdq(iclas) - this%bulk_charge(this%lattice%chargetrf_type(iclas))
         dif = dif + tdq(iclas)
      end do
      nsum = 0
      do j = iex, this%lattice%nbas
         nsum = nsum + 1 ! each atom correspond to only one atom type (replacing the NTIPO of the old code)
      end do
      if (rank == 0) call g_logger%info('nsum is '//fmt('i5', nsum)//' dif is '//fmt('f10.6', dif), __FILE__, __LINE__)
      do j = iex, this%lattice%nbas
         tdq(j) = -dif/nsum
      end do

      if (rank == 0) call g_logger%info('tdq. external atom= '//fmt('f10.6', tdq(iex))//' '//fmt('f10.6', dif), __FILE__, __LINE__)

      if (abs(tdq(iex)) > 0.5) then
         if (rank == 0) call g_logger%warning('too much charge in the external atom! Careful', __FILE__, __LINE__)
      end if

      jj = 0
      do i = 1, this%lattice%nbas
         jj = jj + 1
         totq(jj) = tdq(i)
      end do

      do jbas = 1, this%lattice%nrec
         ss = 0.0d0
         do n = 1, nrx
            ss = ss + totq(n)*this%amad(jbas, n)
         end do
         this%symbolic_atom(this%lattice%nbulk + jbas)%potential%vmad = ss
         this%symbolic_atom(this%lattice%nbulk + jbas)%potential%vmad = ss + this%symbolic_atom(this%lattice%chargetrf_type(jbas))%potential%vmad
      end do

      do iclas = 1, this%lattice%nrec
         vmad0(iclas) = this%symbolic_atom(this%lattice%nbulk + iclas)%potential%vmad
         this%symbolic_atom(this%lattice%nbulk + iclas)%potential%vmad = &
         & this%symbolic_atom(this%lattice%nbulk + iclas)%potential%vmad*this%vmix + vmad0(iclas)*(1 - this%vmix)
         verr = this%symbolic_atom(this%lattice%nbulk + iclas)%potential%vmad - vmad0(iclas)
         if (rank == 0) call g_logger%info('Class '//fmt('i4', iclas)//' Chg. Transfer= '//fmt('f10.6', tdq(iclas))//' VMAD= '&
         &//fmt('f10.6', this%symbolic_atom(this%lattice%nbulk + iclas)%potential%vmad), __FILE__, __LINE__)
      end do
   end subroutine imppot

   subroutine surfpot(this)
      class(charge), intent(inout) :: this
      ! Local variables
      integer :: I, IB, IBAS, ICLAS, IEX, INEQ, IQ, J, JQ, NQ, NRLX, i_all, i_stat
      real(rp) :: DIF, SUM1, SUMM, SUMN, VBULK, VM1, VMAD1, VMARD, VMN, wsm, twooverwsm, wsms
      real(rp), dimension(this%lattice%nrec + this%lattice%nbas) :: TDQ, VM, qst
      real(rp), dimension(:, :), allocatable :: DSS
      integer, parameter :: npot = 5000
      integer :: init, atomrec, k

      !this%vmix = 0.6d0
      wsm = this%lattice%wav*ang2au
      wsms = this%sws*this%lattice%alat*ang2au
      init = 6
      ineq = this%lattice%nbas - init
      nrlx = this%lattice%nbas - init
      iex = this%lattice%nlay + 1
      twooverwsm = 2.0d0/wsm

      do j = iex, ineq
         tdq(j) = 0.0d0
      end do

      qst(:) = 0.0d0
      dif = 0.0d0
      atomrec = 0.0d0
      do iclas = 1, this%lattice%nlay
         tdq(iclas) = 0.0d0
         do k = 1, this%lattice%natoms_layer(iclas)
            atomrec = atomrec + 1
            tdq(iclas) = tdq(iclas) + this%dq(atomrec)
         end do
         qst(iclas) = tdq(iclas)/this%lattice%natoms_layer(iclas)
         dif = tdq(iclas) + dif
      end do
      ! Transfer excess charge to next layer
      tdq(iex) = -dif
      !qst(iex) = -dif/2
      if (rank == 0) call g_logger%info('Excess charge to the next layer is '//fmt('f10.6', tdq(iex)), __FILE__, __LINE__)
      if (abs(tdq(iex)) > 0.5d0) then
         if (rank == 0) call g_logger%warning('Too much charge in the external layer!', __FILE__, __LINE__)
      end if
      ! Calculate the potential for a given charge distribution
      do ibas = 1, this%lattice%nlay
         summ = 0.0d0
         do j = 1, nrlx
            iq = init + ibas
            jq = init + j
            summ = summ + this%dss(iq, jq)*tdq(j) !- twooverwsm*tdq(j) + twooverwsm*qst(j)
         end do
         vm(ibas) = (summ/wsms) !+ twooverwsm*qst(t)
      end do
      ! Print charge and potentials
      sum1 = 0.0d0
      sumn = 0.0d0
      do j = 1, nrlx
         jq = init + j
         sum1 = sum1 + this%dss(1, jq)*tdq(j)! - twooverwsm*tdq(j) + twooverwsm*qst(j)
         sumn = sumn + this%dss(this%lattice%nbas, jq)*tdq(j) !- twooverwsm*tdq(j) + twooverwsm*qst(j)
      end do
      if (rank == 0) call g_logger%info('sum1= '//real2str(sum1)//' sumn= '//real2str(sumn), __FILE__, __LINE__)
      vm1 = (sum1/wsms)
      vmn = (sumn/wsms) !vm(this%lattice%nlay) + twooverwsm*tdq(this%lattice%nlay)!(sumn/wsms)
      vbulk = vmn !+ twooverwsm*qst(this%lattice%nlay)
      vmad1 = vm1 - vbulk
      if (rank == 0) call g_logger%info('vm1= '//real2str(vm1)//' vmn= '//real2str(vmn), __FILE__, __LINE__)

      atomrec = 0
      do ib = 1, this%lattice%nlay
         do k = 1, this%lattice%natoms_layer(ib)
            vmard = 0.0d0
            atomrec = atomrec + 1
            vmard = vm(ib) - vbulk !+ twooverwsm*this%dq(atomrec) !- twooverwsm*tdq(ib)
            this%symbolic_atom(this%lattice%nbulk + atomrec)%potential%vmad = &
            & vmard*this%vmix + this%symbolic_atom(this%lattice%nbulk + atomrec)%potential%vmad*(1 - this%vmix)
         end do
      end do
      do iclas = 1, this%lattice%nrec
         if (rank == 0) call g_logger%info('Class '//fmt('i4', iclas)//' Chg. Transfer= '//fmt('f10.6', this%dq(iclas))//' VMAD= '&
                             &//fmt('f10.6', this%symbolic_atom(this%lattice%nbulk + iclas)%potential%vmad), __FILE__, __LINE__)
      end do
   end subroutine surfpot

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the matrix elements for the electrostatic potential from
   !> information in ´clust´ and ´self´.
   !---------------------------------------------------------------------------
   subroutine bulkmat(this)
      class(charge), intent(inout) :: this
      real(rp) :: alatbulkmat, tol, alat0, vol, awald0, awald
      real(rp), dimension(3, this%lattice%ndim) :: tau
      real(rp), dimension(3, 3) :: rb, qb
      integer :: nsize, nbmx, nkrmx, nkdmx, j7rlat, j7dlat, j7work, j7amad, lmxst, nkr, nkd
      logical :: isopen
      integer :: i, j
      inquire (unit=10, opened=isopen)
      if (isopen) then
         call g_logger%fatal('charge%bulkmat, file ves.out: Unit 10 is already open', __FILE__, __LINE__)
         error stop
      else
         open (unit=10, file='ves.out')
      end if

      inquire (unit=11, opened=isopen)
      if (isopen) then
         call g_logger%fatal('charge%bulkmat, file mad.mat: Unit 11 is already open', __FILE__, __LINE__)
      else
         open (unit=11, file='mad.mat', form='unformatted')
      end if

      ! Defining some parameters. Need to find a way to set them automatically in
      ! the future
      nsize = 3000000
      nbmx = 5500
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('charge.w', this%w, nsize)
#else
      allocate (this%w(nsize))
#endif
      call wkinit(this)

      nkrmx = 3000
      nkdmx = 3000
      awald0 = 3
      tol = 1.0d-06
      alat0 = 0.0

      alatbulkmat = this%lattice%alat/0.52917721d0
      call RDISTN(this%lattice%crd, tau, this%lattice%ntot, this%gx, this%gy, this%gz, this%gt)
      call DEFRR(this, J7RLAT, 3*NKRMX)
      call DEFRR(this, J7DLAT, 3*NKDMX)
      call DEFRR(this, J7WORK, MAX0(NKDMX, NKRMX))
      LMXST = 5
      call LATTC(AWALD0, TOL, alatbulkmat, ALAT0, this%lattice%a, this%GX, this%GY, this%GZ, this%GT, RB, QB, LMXST,     &
   &  VOL, AWALD, this%W(J7DLAT), NKD, this%W(J7RLAT), NKR, NKDMX, NKRMX, this%W(J7WORK))
      call RLSE(this, J7WORK)
      call DEFRR(this, J7AMAD, this%lattice%ntot*this%lattice%ntot)
      call MADMAT(this%lattice%ntot, TAU, AWALD, alatbulkmat, VOL,                              &
   &   this%W(J7RLAT), NKR, this%W(J7DLAT), NKD, this%W(J7AMAD))
      close (10)
      close (11)
   end subroutine bulkmat

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the matrix elements for the electrostatic potential from
   !> information in ´alelay.dat´ and ´sizelay´
   !---------------------------------------------------------------------------
   subroutine surfmat(this)
      class(charge), intent(inout) :: this
      real(rp), dimension(3, 3) :: bs, bk
      real(rp), dimension(3) :: temp_bsx, temp_bsy, temp_bsz
      integer :: i, j

      ! locally defining translational vectors
      bs(1, 1) = this%bsx(1); bs(2, 1) = this%bsy(1); bs(3, 1) = this%bsz(1)
      bs(1, 2) = this%bsx(2); bs(2, 2) = this%bsy(2); bs(3, 2) = this%bsz(2)
      bs(1, 3) = this%bsx(3); bs(2, 3) = this%bsy(3); bs(3, 3) = this%bsz(3)

      ! find primitive translational vectors
      temp_bsx(:) = bs(1, :)
      temp_bsy(:) = bs(2, :)
      temp_bsz(:) = bs(3, :)
      bk(1, :) = cross_product(temp_bsy, temp_bsz)
      bk(2, :) = cross_product(temp_bsz, temp_bsx)
      bk(3, :) = cross_product(temp_bsx, temp_bsy)

      this%vol = abs(this%bsx(1)*bk(1, 1) + this%bsy(1)*bk(2, 1) + this%bsz(1)*bk(3, 1))

      bk(:, :) = (bk(:, :)/this%vol)*2*pi

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('charge.bkx', this%bkx, 3)
      call g_safe_alloc%allocate('charge.bky', this%bky, 3)
      call g_safe_alloc%allocate('charge.bkz', this%bkz, 3)
#else
      allocate (this%bkx(3), this%bky(3), this%bkz(3))
#endif
      this%bkx(1) = bk(1, 1); this%bky(1) = bk(2, 1); this%bkz(1) = bk(3, 1)
      this%bkx(2) = bk(1, 2); this%bky(2) = bk(2, 2); this%bkz(2) = bk(3, 2)
      this%bkx(3) = bk(1, 3); this%bky(3) = bk(2, 3); this%bkz(3) = bk(3, 3)

      ! Atomic Wigner-Seitz radius in units of the lattice spacing
      ! this%lattice%alat
      this%sws = (3.0d0*this%vol/4.0d0/pi/this%nq3)**(1.d0/3.d0)

      ! set up a few parameters
      this%rmax = this%amax/this%alamda
      this%gmax = 2.0d0*this%alamda*this%bmax

      call set2d(this)
      call latt2d(this)
      call madl2r(this)
      call madl2d(this)

      ! Include the on-site term to the dss matrix
      do i = 1, this%lattice%nbas
         this%dss(i, i) = this%dss(i, i) + 2.0d0*(this%sws*this%lattice%alat*ang2au/this%wssurf(i))
      end do

      do i = 1, this%lattice%nbas
         write (120, '(50f10.4)') (this%dss(i, j), j=1, this%lattice%nbas)
      end do
      close (120)
   end subroutine surfmat

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Biulds the primitive vectors which define the surface direction.
   !---------------------------------------------------------------------------
   subroutine build_alelay(this)
      class(charge), intent(inout) :: this
      real(rp), dimension(3) :: v1, v2, x, y, z
      real(rp), dimension(3) :: new_bsx, new_bsy, new_bsz
      real(rp) :: d1, dmin0, dmin2, dmin3
      real(rp) :: amin, amin2
      real(rp) :: zst, zn
      real(rp) :: phi
      real(rp), parameter :: diff = 1.d-4
      real(rp), parameter :: minidiff = 1.d-8
      real(rp), parameter :: minpi = 4*atan(1.0d0) - diff
      integer, dimension(1000) :: bsyat, bsyat2 ! arrays of nearest neighbors
      integer :: i, j, k, at, idx, idx2

      this%nq3 = this%lattice%nbulk_bulk
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('charge.bsx', this%bsx, 3)
      call g_safe_alloc%allocate('charge.bsy', this%bsy, 3)
      call g_safe_alloc%allocate('charge.bsz', this%bsz, 3)
      call g_safe_alloc%allocate('charge.bssx', this%bssx, 3)
      call g_safe_alloc%allocate('charge.bssy', this%bssy, 3)
      call g_safe_alloc%allocate('charge.bssz', this%bssz, 3)
      call g_safe_alloc%allocate('charge.new_x', this%new_x, 3)
      call g_safe_alloc%allocate('charge.new_y', this%new_y, 3)
      call g_safe_alloc%allocate('charge.new_z', this%new_z, 3)
      call g_safe_alloc%allocate('charge.new_kx', this%new_kx, 3)
      call g_safe_alloc%allocate('charge.new_ky', this%new_ky, 3)
      call g_safe_alloc%allocate('charge.new_kz', this%new_kz, 3)
      call g_safe_alloc%allocate('charge.qx3', this%qx3, this%nq3)
      call g_safe_alloc%allocate('charge.qy3', this%qy3, this%nq3)
      call g_safe_alloc%allocate('charge.qz3', this%qz3, this%nq3)
#else
      allocate (this%bsx(3), this%bsy(3), this%bsz(3))
      allocate (this%bssx(3), this%bssy(3), this%bssz(3))
      allocate (this%new_x(3), this%new_y(3), this%new_z(3))
      allocate (this%new_kx(3), this%new_ky(3), this%new_kz(3))
      allocate (this%qx3(this%nq3), this%qy3(this%nq3), this%qz3(this%nq3))
#endif

!     this%qx3 = 0.0d0; this%qy3 = 0.0d0; this%qz3 = 0.0d0

      do i = 1, this%lattice%nbulk_bulk
         this%qx3(i) = this%lattice%crd(1, i) - this%lattice%crd(1, 1)
         this%qy3(i) = this%lattice%crd(2, i) - this%lattice%crd(2, 1)
         this%qz3(i) = this%lattice%crd(3, i) - this%lattice%crd(3, 1)
      end do

      this%a = 1.0d0
      this%b = 1.0d0
      this%c = 1.0d0

      this%amax = 4.0d0
      this%bmax = 4.0d0
      this%alamda = 4.0d0

!   ...... Finds bsx, bsy and bsz for all surfaces ..........

      ! Initial values
      d1 = 1000.d0; dmin0 = d1; dmin2 = d1; dmin3 = d1
      amin = minpi; amin2 = minpi

      ! Finds the central atom
      at = 0
      do i = 1, this%lattice%kk
         do j = 1, 3
            v1(j) = this%lattice%cr(j, i)
            v2(j) = 0.0d0
         end do
         if (distance(v1, v2) .lt. d1) then
            d1 = distance(v1, v2); at = i
         end if
      end do

      ! Finds bsx(:) and bsz(:)
      do i = 1, this%lattice%kk
      if (this%lattice%num(i) == this%lattice%num(at)) then
         do j = 1, 3
            v1(j) = this%lattice%cr(j, i)
            v2(j) = this%lattice%cr(j, at)
         end do
         zn = this%lattice%dx*v1(1) + this%lattice%dy*v1(2) + this%lattice%dz*v1(3)
         zst = this%lattice%dx*v2(1) + this%lattice%dy*v2(2) + this%lattice%dz*v2(3)
         if ((distance(v1, v2) .le. dmin0) .and. (zn .eq. zst) .and. (distance(v1, v2) .gt. diff)) then
            do j = 1, 3
               this%bsx(j) = this%lattice%cr(j, i) - this%lattice%cr(j, at)
            end do
            dmin0 = distance(v1, v2)
         end if
         if ((distance(v1, v2) .gt. dmin0) .and. (zn .eq. zst) .and. (distance(v1, v2) .lt. dmin3) .and. &
   &        (angle(pos(v1, v2), this%bsx) .gt. diff) .and. (angle(pos(v1, v2), this%bsx) .lt. minpi)) then
            dmin3 = distance(v1, v2)
         end if
!       if ((distance(v1, v2).le.dmin2).and.(zn.le.(zst-this%lattice%zstep+minidiff)).and. &
! &         (zn.ge.(zst-this%lattice%zstep-minidiff)).and.(distance(v1, v2).gt.diff)) then
         if ((distance(v1, v2) .le. dmin2) .and. (zn .ne. zst) .and. (distance(v1, v2) .gt. diff)) then
            do j = 1, 3
               this%bsz(j) = this%lattice%cr(j, i) - this%lattice%cr(j, at)
            end do
            dmin2 = distance(v1, v2)
         end if
      end if
      end do
      ! Finds bsy(:)

      idx = 1; idx2 = 1
      bsyat(:) = -1; bsyat2(:) = -1

      ! First, try to find if there are nearest neighbors with the same
      ! minimum distance (as bsx), which is ´dmin0´ in the code

      do i = 1, this%lattice%kk
      if (this%lattice%num(i) == this%lattice%num(at)) then
         do j = 1, 3
            v1(j) = this%lattice%cr(j, i)
            v2(j) = this%lattice%cr(j, at)
         end do
         zn = this%lattice%dx*v1(1) + this%lattice%dy*v1(2) + this%lattice%dz*v1(3)
         zst = this%lattice%dx*v2(1) + this%lattice%dy*v2(2) + this%lattice%dz*v2(3)
         if ((zn .eq. zst) .and. (distance(v1, v2) .eq. dmin0) .and. &
   &     (abs(angle(pos(v1, v2), this%bsx)) .gt. diff) .and. (abs(angle(pos(v1, v2), this%bsx)) .lt. minpi)) then
            bsyat(idx) = i
            idx = idx + 1
         else if ((zn .eq. zst) .and. (distance(v1, v2) .eq. dmin3) .and. &
   &     (abs(angle(pos(v1, v2), this%bsx)) .gt. diff) .and. (abs(angle(pos(v1, v2), this%bsx)) .lt. minpi)) then
            bsyat2(idx) = i
            idx2 = idx2 + 1
         end if
      end if
      end do

      ! If atoms with those conditions are found, then find the one with
      ! the lowest angle

      if (bsyat(1) .gt. 0) then
         do j = 1, idx
            if (bsyat(j) .gt. 0) then
               do k = 1, 3
                  v1(k) = this%lattice%cr(k, bsyat(j))
                  v2(k) = this%lattice%cr(k, at)
               end do
               if (abs(angle(pos(v1, v2), this%bsx)) .lt. amin) then
                  do k = 1, 3
                     this%bsy(k) = this%lattice%cr(k, bsyat(j)) - this%lattice%cr(k, at)
                  end do
                  amin = abs(angle(pos(v1, v2), this%bsx))
               end if
            end if
         end do

         ! If not, then the atoms should be in the second minor distance,
         ! which is ´dmin3´ in the code

      else if ((bsyat2(1) .gt. 0) .and. (bsyat(1) .lt. 0)) then
         do j = 1, idx2
            if (bsyat2(j) .gt. 0) then
               do k = 1, 3
                  v1(k) = this%lattice%cr(k, bsyat2(j))
                  v2(k) = this%lattice%cr(k, at)
               end do
               if (abs(angle(pos(v1, v2), this%bsx)) .lt. amin2) then
                  do k = 1, 3
                     this%bsy(k) = this%lattice%cr(k, bsyat2(j)) - this%lattice%cr(k, at)
                  end do
                  amin2 = abs(angle(pos(v1, v2), this%bsx))
               end if
            end if
         end do
      end if

      ! Rotates bsx and bsy to the xy plane if perp axis is diff from [001]

      if ((this%lattice%dx .ne. 0) .or. (this%lattice%dy .ne. 0)) then

         x = [1.0d0, 0.0d0, 0.0d0]
         y = [0.0d0, 1.0d0, 0.0d0]
         z = [0.0d0, 0.0d0, 1.0d0]

         this%new_z(1) = this%lattice%dx; this%new_z(2) = this%lattice%dy; this%new_z(3) = this%lattice%dz
         phi = angle(z, this%new_z)
         this%new_x = normalize(erodrigues(x, cross_product(z, this%new_z), phi))
         this%new_y = normalize(erodrigues(y, cross_product(z, this%new_z), phi))
         this%new_z = normalize(this%new_z)

         new_bsx(1) = dot_product(this%bsx, this%new_x)
         new_bsx(2) = dot_product(this%bsx, this%new_y)
         new_bsx(3) = dot_product(this%bsx, this%new_z)

         this%bsx(1) = dot_product(new_bsx, x)
         this%bsx(2) = dot_product(new_bsx, y)
         this%bsx(3) = dot_product(new_bsx, z)

         new_bsy(1) = dot_product(this%bsy, this%new_x)
         new_bsy(2) = dot_product(this%bsy, this%new_y)
         new_bsy(3) = dot_product(this%bsy, this%new_z)

         this%bsy(1) = dot_product(new_bsy, x)
         this%bsy(2) = dot_product(new_bsy, y)
         this%bsy(3) = dot_product(new_bsy, z)

         new_bsz(1) = dot_product(this%bsz, this%new_x)
         new_bsz(2) = dot_product(this%bsz, this%new_y)
         new_bsz(3) = dot_product(this%bsz, this%new_z)

         this%bsz(1) = dot_product(new_bsz, x)
         this%bsz(2) = dot_product(new_bsz, y)
         this%bsz(3) = dot_product(new_bsz, z)

      end if

      ! Change to old definitions

      this%bssx(:) = this%bsx(:); this%bssy(:) = this%bsy(:); this%bssz(:) = this%bsz(:)
      this%bsx(1) = this%bssx(1); this%bsx(2) = this%bssy(1); this%bsx(3) = this%bssz(1)
      this%bsy(1) = this%bssx(2); this%bsy(2) = this%bssy(2); this%bsy(3) = this%bssz(2)
      this%bsz(1) = this%bssx(3); this%bsz(2) = this%bssy(3); this%bsz(3) = this%bssz(3)
   end subroutine build_alelay

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the matrix elements for the electrostatic potential from
   !> information in ´clust´, ´size´ and ´self´.
   !---------------------------------------------------------------------------
   subroutine impmad(this)
      implicit none
      class(charge) :: this
      !> Local variables
      real(rp), dimension(3) :: ddum
      real(rp) :: r2, dd
      integer :: i, j, jj, ii, l, kk
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('charge.wsimp', this%wsimp, this%lattice%nbas)
      call g_safe_alloc%allocate('charge.amad', this%amad, (/this%lattice%nbas, this%lattice%nbas/))
#else
      allocate (this%wsimp(this%lattice%nbas), this%amad(this%lattice%nbas, this%lattice%nbas))
#endif
      !Set as if all the atoms have the same WS radius. Can be improved later
      this%wsimp(:) = this%lattice%wav*ang2au
      r2 = 0
      j = 0
!---determining number of atoms
!    do ii = 1, self_obj%nbas
!      DO JJ = 1, TIPO(II)
!        J=J+1
!      END DO
!    END DO
!    NRX=J
!    if (self_obj%nbas>this%lattice%nrec) then
!      write(*, *) "Too many shells in ´self´"
!      STOP
!    end if
!    if (NRX>NRLX) then
!      write(*, *) "To many atoms in shells, increase NRLX"
!      STOP
!    end if
      !
!----1st site
      j = 1
      ii = 1
      this%amad(ii, j) = (2_rp/this%wsimp(ii))
      do i = 1, this%lattice%nbas
         if (i .ne. j) then
            r2 = 0.0d0
            do l = 1, 3
               ddum(l) = (this%lattice%cr(l, j) - this%lattice%cr(l, i))*this%lattice%alat
               r2 = ddum(l)*ddum(l) + r2
            end do
            r2 = r2*(ang2au**2)
            dd = 1_rp/(dsqrt(r2))
            this%amad(ii, i) = 2_rp*dd
         end if
      end do
!---end 1st site
      do ii = 1, this%lattice%nbas - 1
!        DO JJ = 1, TIPO(II)
         j = j + 1   !looping to find first atom in 'shell'
!        END DO
         this%amad(ii + 1, j) = (2_rp/this%wsimp(ii + 1))
         do i = 1, this%lattice%nbas
            if (i .ne. j) then
               r2 = 0.0d0
               do l = 1, 3
                  ddum(l) = (this%lattice%cr(l, j) - this%lattice%cr(l, i))*this%lattice%alat
                  r2 = ddum(l)*ddum(l) + r2
               end do
               r2 = r2*(ang2au**2)
!         print ´(a, 2i4, f10.3, 3f10.6, 2x, 3f10.6)´, ´II, I´, J, I, R2,
!    .        this%lattice%crd(:, J)/ALAT, this%lattice%crd(:, I)/ALAT
               dd = 1_rp/(dsqrt(r2))
               this%amad(ii + 1, i) = 2_rp*dd
            end if
         end do
      end do
!    do kk=1, this%lattice%nbas
!      write(7, 10) (this%amad(i, kk), I=1, this%lattice%nbas)
!    end do
!  10 FORMAT(2X, 3E16.6)
   end subroutine impmad

! Subroutines from surfmat library
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the Madelung potential matrix for a finite number
   !> of 2D layers of point charges and dipoles. This is used to
   !> set up the Madelung potential in the surface program INF.
   !>
   !> Matrix    (l´, l)
   !>
   !> DSS       (s, s)
   !> DSZ       (s, z)
   !> DS3Z2     (s, 3z2-1)
   !> DSX2Y2    (s, x2-y2)
   !> DSXY      (s, xy)
   !> DZZ       (z, z)
   !> DZ3Z2     (z, 3zz-1)
   !>
   !> PM    : Plate-condenser potential matrix
   !---------------------------------------------------------------------------
   subroutine madl2d(this)

      ! Input
      class(charge), intent(inout) :: this
      ! Local variables
      integer :: I, IQ, JQ, JQP, NQM, nbas
      real(rp) :: ALPHA, AQPPZ, BETA, BETAM, BETAP, BMDL, BMG, DGI, &
                  DMDL, DMG, DRI, DRI2, DRI3, DRI4, DZ, ERFCA, ERFCM, &
                  ERFCP, EXPA, EXPAA, EXPEK, EXPP, EXPZ, FACBET, &
                  FACDIF, FACDK, FACDR, FACERF, FACG1, FACG2, &
                  FACGAU, FACGUA, FACP, FACQ1, FACQ2, FACQ3, FACQ4, &
                  FACQUA, FACQUR, FACZZ
      real(rp) :: G0MDL, GAM52, GLH, GLK, GLM, GM0, PHASE, Q0MDL, &
                  QIMDL, QM0G, QMIG, QMRG, QPPX, QPPY, QPPZ, QPX, QPY, &
                  QPZ, QRMDL, SUM0G, SUM0R, SUM1G, SUM1R, SUM20G, &
                  SUM20R, SUM2IG, SUM2IR, SUM2RG, SUM2RR, SUM30G, &
                  SUM30R, SUMM, SUMQ0, SUMQI, SUMQR, TWOLAM, TWOS, X, &
                  XI, XR
      real(rp) :: Y, Z, ZU, exf, exh, expm

      ! Process described in H. L. Skriver and N. M. Rosengaard Phys. Rev. B 43, 9538

      nbas = this%lattice%nbas
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('charge.ds3z2', this%ds3z2, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.dsx2y2', this%dsx2y2, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.dsxy', this%dsxy, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.dss', this%dss, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.dsz', this%dsz, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.dzz', this%dzz, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.dz3z2', this%dz3z2, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.am', this%am, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.bm', this%bm, (/nbas, nbas/))
      call g_safe_alloc%allocate('charge.pm', this%pm, (/nbas, nbas/))
#else
      allocate (this%ds3z2(nbas, nbas), this%dsx2y2(nbas, nbas), this%dsxy(nbas, nbas), this%dss(nbas, nbas), this%dsz(nbas, nbas) &
                , this%dzz(nbas, nbas), this%dz3z2(nbas, nbas), this%am(nbas, nbas), this%bm(nbas, nbas), this%pm(nbas, nbas))
#endif

      TWOS = 2.*this%SWS
      TWOLAM = 2.*this%ALAMDA
      FACBET = PI/this%AR2D/TWOLAM
      FACGAU = -2.*SQRT_PI/this%AR2D/this%ALAMDA
      FACQUA = 2.*this%SWS*this%SWS*SQRT(5.d0)/3.
      FACGUA = 6.*this%SWS*this%SWS*this%SWS/SQRT(15.d0)
      FACQUR = FACQUA*SQRT(3.d0/2.d0)
      FACERF = 2.*PI/this%AR2D
      FACDIF = 8.*PI*this%SWS/this%AR2D
      FACP = this%SWS*SQRT(3.d0)
      FACDK = FACP*PI/this%AR2D
      FACDR = FACP*2./SQRT_PI
      FACQ1 = 4.*SQRT_PI*this%ALAMDA/this%AR2D
      FACQ2 = 1.5*PI/this%AR2D
      FACQ3 = 2.d0/SQRT_PI
      FACQ4 = 0.5*PI/this%AR2D
      FACG1 = 8.*SQRT_PI*this%ALAMDA*this%ALAMDA/this%AR2D
      FACG2 = 2.5*PI/this%AR2D
      FACZZ = -6.d0/SQRT(5.d0)

!    write(100, ´(18f10.4)´) twos, twolam, facbet, facgau, facqua, facqur, facerf, facdif, &
!               facp, facdk, facdr, facq1, facq2, facq3, facq4, facg1, facg2, faczz

      !
      !     Layer-diagonal terms (R=R´)
      !
      SUMM = 0.d0
      SUMQ0 = -FACQ1
      SUMQR = 0.d0
      SUMQI = 0.d0
      ! Summation in the reciprocal space
      do I = 2, this%NUMVG
         DGI = this%DG(I)
         BETA = DGI/TWOLAM
         ERFCA = ERFC(BETA)
         SUMM = SUMM + 2.*ERFCA/BETA
         SUMQ0 = SUMQ0 - FACQ1*EXP(-BETA*BETA) + 2.*FACQ2*DGI*ERFCA
         SUMQR = SUMQR - 2.*FACQ4*(this%AKX(I)*this%AKX(I) - this%AKY(I)*this%AKY(I))*ERFCA/DGI
         SUMQI = SUMQI - 4.*FACQ4*this%AKX(I)*this%AKY(I)*ERFCA/DGI
      end do
      BMDL = FACBET*SUMM
      Q0MDL = FACQUA*SUMQ0
      QRMDL = FACQUR*SUMQR
      QIMDL = FACQUR*SUMQI
      !
      SUMM = 0.d0
      SUMQ0 = 0.d0
      SUMQR = 0.d0
      SUMQI = 0.d0
      ! Summation in the translational vector of the 2D lattice
      do I = 2, this%NR0
         DRI = this%DR(I)
         DRI2 = DRI*DRI
         DRI3 = DRI2*DRI
         ALPHA = this%ALAMDA*DRI
         ERFCA = ERFC(ALPHA)
         EXPAA = EXP(-ALPHA*ALPHA)
         GAM52 = 1.5*(ERFCA/FACQ3 + ALPHA*EXPAA) + EXPAA*ALPHA**3
         GAM52 = GAM52/DRI3
         SUMM = SUMM + ERFCA/DRI
         SUMQ0 = SUMQ0 - GAM52
         SUMQR = SUMQR + GAM52*(this%ASX(I)*this%ASX(I) - this%ASY(I)*this%ASY(I))/DRI2
         SUMQI = SUMQI + GAM52*2.*this%ASX(I)*this%ASY(I)/DRI2
      end do
      BMDL = BMDL + SUMM - TWOLAM/SQRT_PI
      Q0MDL = Q0MDL + FACQUA*FACQ3*SUMQ0
      QRMDL = QRMDL + FACQUR*FACQ3*SUMQR
      QIMDL = QIMDL + FACQUR*FACQ3*SUMQI
      !
      do IQ = 1, this%lattice%nbas
         this%AM(IQ, IQ) = FACGAU
         this%BM(IQ, IQ) = BMDL
         this%DSZ(IQ, IQ) = 0.d0
         this%DS3Z2(IQ, IQ) = Q0MDL
         this%DSX2Y2(IQ, IQ) = QRMDL
         this%DSXY(IQ, IQ) = QIMDL
         this%DZ3Z2(IQ, IQ) = 0.d0
      end do
      !
      !     Off-diagonal terms
      !
      NQM = this%lattice%nbas - 1
      do JQ = 1, NQM
         JQP = JQ + 1
         QPX = this%QX(JQ)
         QPY = this%QY(JQ)
         QPZ = this%QZ(JQ)
         do IQ = JQP, this%lattice%nbas
            QPPX = this%QX(IQ) - QPX
            QPPY = this%QY(IQ) - QPY
            QPPZ = this%QZ(IQ) - QPZ
            !
            !     Contribution to potential from k-parallel equal to zero
            !
            DZ = this%ALAMDA*QPPZ
            ERFCP = ERFC(DZ)
            ERFCM = 2.-ERFCP
            if (DZ > 12.) then
               EXPZ = 0.d0
            else
               EXPZ = EXP(-DZ*DZ)
            end if
            this%AM(IQ, JQ) = FACGAU*EXPZ - QPPZ*FACERF*ERFCM
            this%AM(JQ, IQ) = FACGAU*EXPZ + QPPZ*FACERF*ERFCP
            SUM0G = 0.d0
            SUM1G = ERFCM - ERFCP
            SUM20G = -FACQ1*EXPZ
            SUM2RG = 0.d0
            SUM2IG = 0.d0
            SUM30G = FACG1*DZ*EXPZ
            !
            !     Contributions from k-parallel greater than zero
            !
            do I = 2, this%NUMVG
               DGI = this%DG(I)
               PHASE = COS(this%AKX(I)*QPPX + this%AKY(I)*QPPY)
               BETA = DGI/TWOLAM
               AQPPZ = this%ALAMDA*QPPZ
               BETAP = BETA + AQPPZ
               BETAM = BETA - AQPPZ
               ERFCP = ERFC(BETAP)
               ERFCM = ERFC(BETAM)
               ! alterado por Peduto para evitar alguns casos de "overflow" que podem
               ! ocorrer quando, principalmente, o sistema se tratar de um VAX.
               ! a variavel EXPP estoura, geralmente, quando ercp e erfcm sao nulos.
               ! Mas
               ! nestes casos nao eh necessario calcular EXPP. Assim evita-se que o
               ! progra
               ! ma seja abortado.
               if (erfcp == 0) then
                  expm = exp(-dgi*qppz)
                  exf = expm*erfcm
                  exh = -exh
               else
                  EXPP = EXP(DGI*QPPZ)
                  EXPM = 1./EXPP
                  EXF = EXPP*ERFCP + EXPM*ERFCM
                  EXH = EXPP*ERFCP - EXPM*ERFCM
               end if
               EXPEK = EXP(-AQPPZ*AQPPZ)*EXP(-BETA*BETA)
               BMG = PHASE*EXF/BETA
               DMG = PHASE*EXH
               QM0G = PHASE*(FACQ2*DGI*EXF - FACQ1*EXPEK)
               QMRG = PHASE*(this%AKX(I)*this%AKX(I) - this%AKY(I)*this%AKY(I))*EXF/DGI
               QMIG = PHASE*this%AKX(I)*this%AKY(I)*EXF/DGI
               GM0 = PHASE*(FACG2*DGI*DGI*EXH + FACG1*DZ*EXPEK)
               SUM0G = SUM0G + BMG
               SUM1G = SUM1G - DMG
               SUM20G = SUM20G + QM0G
               SUM2RG = SUM2RG + QMRG
               SUM2IG = SUM2IG + QMIG
               SUM30G = SUM30G + GM0
            end do
            BMDL = FACBET*SUM0G
            DMDL = FACDK*SUM1G
            Q0MDL = FACQUA*SUM20G
            QRMDL = FACQUR*FACQ4*SUM2RG
            QIMDL = FACQUR*FACQ4*2.*SUM2IG
            G0MDL = FACGUA*SUM30G
            !
            SUM0R = 0.d0
            SUM1R = 0.d0
            SUM20R = 0.d0
            SUM2RR = 0.d0
            SUM2IR = 0.d0
            SUM30R = 0.d0
            do I = 1, this%NUMVR
               X = this%ASX(I) + QPPX
               Y = this%ASY(I) + QPPY
               Z = QPPZ
               DRI2 = X*X + Y*Y + Z*Z
               DRI = SQRT(DRI2)
               DRI3 = DRI*DRI2
               DRI4 = DRI*DRI3
               ZU = -Z/DRI
               XR = (X*X - Y*Y)/DRI2
               XI = 2.*X*Y/DRI2
               if (DRI < this%RMAX) then
                  ALPHA = this%ALAMDA*DRI
                  ERFCA = ERFC(ALPHA)
                  SUM0R = SUM0R + ERFCA/DRI
                  EXPA = EXP(-ALPHA*ALPHA)
                  GLH = 0.5*SQRT_PI*ERFCA + ALPHA*EXPA
                  GLK = 1.5*GLH + EXPA*ALPHA**3
                  GLM = 2.5*GLK + EXPA*ALPHA**5
                  SUM1R = SUM1R - GLH*ZU/DRI2
                  SUM20R = SUM20R + GLK*(3.*ZU*ZU - 1.)/DRI3
                  SUM2RR = SUM2RR + GLK*XR/DRI3
                  SUM2IR = SUM2IR + GLK*XI/DRI3
                  SUM30R = SUM30R + GLM*ZU*(5.*ZU*ZU - 3.)/DRI4
               end if
            end do
            BMDL = BMDL + SUM0R
            DMDL = DMDL + FACDR*SUM1R
            Q0MDL = Q0MDL + FACQUA*FACQ3*SUM20R
            QRMDL = QRMDL + FACQUR*FACQ3*SUM2RR
            QIMDL = QIMDL + FACQUR*FACQ3*SUM2IR
            G0MDL = G0MDL + FACGUA*2.*FACQ3*SUM30R
            this%BM(IQ, JQ) = BMDL
            this%BM(JQ, IQ) = BMDL
            this%DSZ(IQ, JQ) = DMDL
            this%DSZ(JQ, IQ) = -DMDL
            this%DS3Z2(IQ, JQ) = Q0MDL
            this%DS3Z2(JQ, IQ) = Q0MDL
            this%DSX2Y2(IQ, JQ) = QRMDL
            this%DSX2Y2(JQ, IQ) = QRMDL
            this%DSXY(IQ, JQ) = QIMDL
            this%DSXY(JQ, IQ) = QIMDL
            this%DZZ(IQ, JQ) = FACZZ*Q0MDL
            this%DZZ(JQ, IQ) = FACZZ*Q0MDL
            this%DZ3Z2(IQ, JQ) = G0MDL
            this%DZ3Z2(JQ, IQ) = -G0MDL
            !
         end do
      end do
      !
      !     Potential in units of e**2/(2S)
      !
      do JQ = 1, this%lattice%nbas
         do IQ = 1, this%lattice%nbas
            this%AM(IQ, JQ) = TWOS*this%AM(IQ, JQ)
            this%BM(IQ, JQ) = TWOS*this%BM(IQ, JQ)
            this%DSS(IQ, JQ) = this%AM(IQ, JQ) + this%BM(IQ, JQ)
            this%DSZ(IQ, JQ) = TWOS*this%DSZ(IQ, JQ)
            this%DS3Z2(IQ, JQ) = TWOS*this%DS3Z2(IQ, JQ)
            this%DSX2Y2(IQ, JQ) = TWOS*this%DSX2Y2(IQ, JQ)
            this%DSXY(IQ, JQ) = TWOS*this%DSXY(IQ, JQ)
            this%DZZ(IQ, JQ) = FACZZ*this%DS3Z2(IQ, JQ)
            this%DZ3Z2(IQ, JQ) = TWOS*this%DZ3Z2(IQ, JQ)
         end do
      end do
      !
      !     Plate-condenser matrix
      !
      do JQ = 1, this%lattice%nbas
         do IQ = 1, this%lattice%nbas
            if (IQ > JQ) then
               this%PM(IQ, JQ) = FACDIF*(this%QZ(JQ) - this%QZ(IQ))
            else
               this%PM(IQ, JQ) = 0.d0
            end if
         end do
      end do
   end subroutine madl2d

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TODO
   !---------------------------------------------------------------------------
   subroutine madl2r(this)
      class(charge), intent(inout) :: this
      ! Local variables
      real(rp) :: alpha, arg, beta, dgi, dri, dz, epsz, erfca, erfcm, &
                  erfcp, erfcz, erfr, expp, facbet, facgau, phase, r, &
                  rx, ry, sumg, sumr, twolam, twos, x, y, z
      real(rp), dimension(51) :: vg, vr, vz, vz0, zz
      integer :: i, iz, nz

      TWOS = 2.*this%SWS
      TWOLAM = 2.*this%ALAMDA
      FACBET = TWOS*PI/this%AR2D
      FACGAU = 1./SQRT_PI/this%ALAMDA
      NZ = 21
      DZ = TWOS/(NZ - 1)
      X = 0.
      Y = 0.
      do IZ = 1, NZ
         Z = -this%SWS + (IZ - 1)*DZ
         ZZ(IZ) = Z
         !
         !     Diagonal terms
         !
         EPSZ = this%ALAMDA*Z
         ARG = -EPSZ*EPSZ
         ERFCZ = ERFC(-EPSZ)
         SUMG = -2.*(Z*ERFCZ + FACGAU*EXP(ARG))
         VZ0(IZ) = FACBET*SUMG
         do I = 2, this%NUMVG
            DGI = this%DG(I)
            ARG = this%AKX(I)*X + this%AKY(I)*Y
            PHASE = COS(ARG)
            BETA = DGI/TWOLAM
            ERFCP = ERFC(BETA + EPSZ)
            ERFCM = ERFC(BETA - EPSZ)
            ! alteracao feita por Peduto para evitar alguns casos de "overflow" que
            ! podem
            ! ocorrer principalmente quando o sistema se tratar de um VAX
            EXPP = EXP(DGI*Z)
            SUMG = SUMG + PHASE*(EXPP*ERFCP + ERFCM/EXPP)/DGI
         end do
         VG(IZ) = FACBET*SUMG
         R = SQRT(X*X + Y*Y + Z*Z)
         ARG = this%ALAMDA*R
         if (ABS(Z) < 1.d-6) then
            ERFR = TWOLAM/SQRT_PI
         else
            ERFR = (1.-ERFC(ARG))/R
         end if
         SUMR = -ERFR
         do I = 2, this%NR0
            RX = this%ASX(I) - X
            RY = this%ASY(I) - Y
            DRI = SQRT(RX*RX + RY*RY + Z*Z)
            ALPHA = this%ALAMDA*DRI
            ERFCA = ERFC(ALPHA)
            SUMR = SUMR + ERFCA/DRI
         end do
         VR(IZ) = TWOS*SUMR
         VZ(IZ) = FACBET*SUMG + TWOS*SUMR
      end do
   end subroutine madl2r

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TODO
   !---------------------------------------------------------------------------
   subroutine latt2d(this)
      class(charge), intent(inout) :: this
      ! Local variables
      integer, parameter :: mr = 50000
      real(rp) :: a, amin, b, da, db, ddm, dkm, dq, dx, g1, ga, gx, gy, &
                  pqx, pqy, pqz, r1, ra, sx, sy, x, y, z
      real(rp), dimension(3) :: dd, dk
      real(rp), dimension(mr) :: csx, csy, d
      integer :: i, iq, jq, k, l, m, n, n1, ng, nr, nsh, nshl, numg, numgh, numr, numrh

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('charge.asx', this%asx, mr)
      call g_safe_alloc%allocate('charge.asy', this%asy, mr)
      call g_safe_alloc%allocate('charge.asz', this%asz, mr)
      call g_safe_alloc%allocate('charge.akx', this%akx, mr)
      call g_safe_alloc%allocate('charge.aky', this%aky, mr)
      call g_safe_alloc%allocate('charge.akz', this%akz, mr)
      call g_safe_alloc%allocate('charge.dr', this%dr, mr)
      call g_safe_alloc%allocate('charge.dg', this%dg, mr)
#else
      allocate (this%asx(mr), this%asy(mr), this%asz(mr), this%akx(mr), this%aky(mr), this%akz(mr), this%dr(mr), this%dg(mr))
#endif

      ! Calculate the radius RA of the sphere holding all the vectors  used in the
      ! lattice summations. R1 is the longest basis vector
      r1 = 1.e-06
      do IQ = 1, this%lattice%nbas
         PQX = this%QX(IQ)
         PQY = this%QY(IQ)
         PQZ = this%QZ(IQ)
         do JQ = IQ, this%lattice%nbas
            X = PQX - this%QX(JQ)
            Y = PQY - this%QY(JQ)
            Z = PQZ - this%QZ(JQ)
            DQ = SQRT(X*X + Y*Y + Z*Z)
            if (DQ >= R1) then
               R1 = DQ
            end if
         end do
      end do
      R1 = R1*1.001
      RA = this%RMAX + R1
      G1 = 0.d0
      GA = this%GMAX + G1

      do I = 1, 3
         DD(I) = SQRT(this%BSX(I)**2 + this%BSY(I)**2 + this%BSZ(I)**2)
         DK(I) = SQRT(this%BKX(I)**2 + this%BKY(I)**2 + this%BKZ(I)**2)
      end do
      DDM = MAX(DD(1), DD(2), DD(3))
      DKM = MAX(DK(1), DK(2), DK(3))
      DDM = 2*PI/DDM
      DKM = 2*PI/DKM
      NUMR = 2*(INT(RA/DKM) + 1) + 1
      NUMG = 2*(INT(GA/DDM) + 1) + 1
      NUMRH = NUMR/2 + 1
      NUMGH = NUMG/2 + 1
      ! Real space
      NR = 0
      this%NR0 = 0
      do L = 1, NUMR
         A = L - NUMRH
         do M = 1, NUMR
            B = M - NUMRH
            SX = A*this%BSX(1) + B*this%BSX(2)
            SY = A*this%BSY(1) + B*this%BSY(2)
            DX = SQRT(SX*SX + SY*SY)
            if (DX <= RA) then
               if (DX <= this%RMAX) then
                  this%NR0 = this%NR0 + 1
               end if
               NR = NR + 1
               if (NR > MR) then
                  write (6, 10005) NR, DX, RA, this%RMAX
                  stop
               else
                  D(NR) = DX
                  CSX(NR) = SX
                  CSY(NR) = SY
               end if
            end if
         end do
      end do
      !
      !     Sort vectors in order of increasing length
      !
      DA = 1.d-06
      NSH = 0
      NSHL = -1
      do K = 1, NR
         AMIN = 1000.
         do N = 1, NR
            if (D(N) < AMIN) then
               AMIN = D(N)
               N1 = N
            end if
         end do
         NSHL = NSHL + 1
         this%ASX(K) = CSX(N1)
         this%ASY(K) = CSY(N1)
         this%ASZ(K) = 0.d0
         DB = D(N1)
         this%DR(K) = DB
         if (DB > DA + 1.d-06) then
            NSH = NSH + 1
            NSHL = 0
            DA = DB
         end if
         D(N1) = 1000.
      end do
      NSH = NSH + 1
      NSHL = NSHL + 1
      this%NUMVR = NR
      !
      !     Reciprocal space
      !
      NG = 0
      do L = 1, NUMG
         A = L - NUMGH
         do M = 1, NUMG
            B = M - NUMGH
            GX = A*this%BKX(1) + B*this%BKX(2)
            GY = A*this%BKY(1) + B*this%BKY(2)
            DX = SQRT(GX*GX + GY*GY)
            if (DX <= GA) then
               NG = NG + 1
               if (NG > MR) then
                  write (6, 10006) NG, DX, GA, this%GMAX
                  stop
               else
                  D(NG) = DX
                  CSX(NG) = GX
                  CSY(NG) = GY
               end if
            end if
         end do
      end do
      !
      !     Sort vectors in order of increasing length
      !
      DA = 1.e-06
      NSH = 0
      NSHL = -1
      do K = 1, NG
         AMIN = 1000.
         do N = 1, NG
            if (D(N) < AMIN) then
               AMIN = D(N)
               N1 = N
            end if
         end do
         NSHL = NSHL + 1
         this%AKX(K) = CSX(N1)
         this%AKY(K) = CSY(N1)
         this%AKZ(K) = 0.d0
         DB = D(N1)
         this%DG(K) = DB
         if (DB > DA*1.000001) then
            NSH = NSH + 1
            NSHL = 0
            DA = DB
         end if
         D(N1) = 1000.
      end do
      NSH = NSH + 1
      NSHL = NSHL + 1
      this%NUMVG = NG

10005 format( &
         /, " LATT2D:** NR =", i5, " exceeds MR. Decrease RMAX by", &
         " changing ALAMDA or AMAX.", /, 11x, "Last vector: (ASX, ASY, ", "ASZ) =" &
         , 3f10.4)
10006 format( &
         /, " LATT2D:** NG =", i5, " exceeds MK. Decrease GMAX by", &
         " changing ALAMDA or BMAX.", /, 11x, "Last vector: (AKX, AKY, ", "AKZ) =" &
         , 3f10.4)
   end subroutine latt2d

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TODO
   !---------------------------------------------------------------------------
   subroutine set2d(this)
      class(charge) :: this
      ! Local variables
      real(rp) :: bx, by, bz, dperp
      real(rp), dimension(this%lattice%nbas*this%nq3) :: pz, uqx, uqy, uqz
      integer, dimension(this%lattice%nbas*this%nq3) :: index
      integer :: i, ib, inx, iq, isrf, n, nlam, nlama, nlamb, nq

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('charge.qx', this%qx, this%lattice%nbas)
      call g_safe_alloc%allocate('charge.qy', this%qy, this%lattice%nbas)
      call g_safe_alloc%allocate('charge.qz', this%qz, this%lattice%nbas)
#else
      allocate (this%qx(this%lattice%nbas), this%qy(this%lattice%nbas), this%qz(this%lattice%nbas))
#endif

      NLAM = this%lattice%nbas
      NLAMB = NLAM/2
      if (2*NLAMB == NLAM) then
         NLAMA = NLAMB - 1
      else
         NLAMA = NLAMB
      end if
      this%AR2D = this%BSX(1)*this%BSY(2) - this%BSY(1)*this%BSX(2)
      this%AR2D = sqrt(this%ar2d**2)
      DPERP = this%VOL/this%AR2D

      !
      !     Establish the interface layers
      !
      N = 0
      do IB = -NLAMA, NLAMB
         BX = IB*this%BSX(3)
         BY = IB*this%BSY(3)
         BZ = IB*this%BSZ(3)
         do IQ = 1, this%NQ3
            N = N + 1
            UQX(N) = BX + this%QX3(IQ)
            UQY(N) = BY + this%QY3(IQ)
            UQZ(N) = BZ + this%QZ3(IQ)
            PZ(N) = UQZ(N)
         end do
      end do
      call QSORT(PZ, INDEX, N)
      ISRF = 0
      do I = 1, N
         INX = INDEX(I)
         if (ABS(PZ(I)) < 1.d-6 .and. ISRF == 0) then
            ISRF = I
         end if
      end do
      IQ = 0
      do I = ISRF - NLAMA, ISRF + NLAMB
         INX = INDEX(I)
         IQ = IQ + 1
         this%QX(IQ) = UQX(INX)
         this%QY(IQ) = UQY(INX)
         this%QZ(IQ) = UQZ(INX)
      end do
   end subroutine set2d

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TODO
   !
   !> @param[inout] m
   !> @param[inout] inx
   !> @param[inout] nl
   !---------------------------------------------------------------------------
   subroutine qsort(m, inx, nl)
      ! sorts the vector m in increasing order
      implicit none
      ! inputs and outputs
      integer, intent(in) :: NL
      real(rp), dimension(NL), intent(inout) :: M
      integer, dimension(NL), intent(inout) :: INX
      ! Local variables
      integer :: I, IND, INIC, J, K
      real(rp) :: FIM, Z

      IND = 1
      INIC = 2
      FIM = NL
      do I = 1, NL
         INX(I) = I
      end do
      do while ((IND == 1) .and. (INIC <= FIM))
         IND = 0
         do J = INT(FIM), INIC, -1
            if (M(J) < M(J - 1)) then
               Z = M(J)
               K = INX(J)
               M(J) = M(J - 1)
               INX(J) = INX(J - 1)
               M(J - 1) = Z
               INX(J - 1) = K
            end if
         end do
         FIM = FIM - 1
         do J = INIC, INT(FIM)
            if (M(J + 1) < M(J)) then
               Z = M(J + 1)
               K = INX(J + 1)
               M(J + 1) = M(J)
               INX(J + 1) = INX(J)
               M(J) = Z
               INX(J) = K
               IND = 1
            end if
         end do
         INIC = INIC + 1
      end do
   end subroutine qsort

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Subroutines from bulkmat library
   !
   !> @param[in] A
   !> @param[in] B
   !> @param[in] C
   !> @return double
   !---------------------------------------------------------------------------
   DOUBLE PRECISION function tripl(A, B, C)
      real(rp), dimension(3) :: A(3), B(3), C(3)
      TRIPL = A(1)*B(2)*C(3) + A(2)*B(3)*C(1) + A(3)*B(1)*C(2) &
     &     - A(3)*B(2)*C(1) - A(2)*B(1)*C(3) - A(1)*B(3)*C(2)
      RETURN
   end function tripl

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TODO
   !
   !> @param[in] X
   !> @return double
   !---------------------------------------------------------------------------
   DOUBLE PRECISION function DERFC(X)
      real(rp) :: t, x

      T = 1.D0/(1.D0 + 0.3275911D0*DABS(X))
      DERFC = DEXP(-X*X)*T*(0.254829592D0 + T*(-.284496736D0 + T* &
   &   (1.421413741D0 + T*(-1.453152027D0 + T*1.061405429D0))))
      IF (X .LT. 0.D0) DERFC = 2.D0 - DERFC
      RETURN
   end function derfc

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TODO
   !
   !> @param[inout] NBAS
   !> @param[inout] TAU
   !> @param[inout] A
   !> @param[inout] ALAT
   !> @param[inout] VOL
   !> @param[inout] RLAT
   !> @param[inout] NKR
   !> @param[inout] DLAT
   !> @param[inout] NKD
   !> @param[inout] AMAD
   !---------------------------------------------------------------------------
   subroutine madmat(NBAS, TAU, A, ALAT, VOL, RLAT, NKR, DLAT, NKD, AMAD)
! ....MAKES MADELUNG MATRIX
      ! Input and output
      integer, parameter :: nbmx = 550
      integer, intent(inout) :: nkr, nkd, nbas
      real(rp), dimension(3, nbmx), intent(inout) :: tau, dlat, rlat
      real(rp), dimension(nbas, nbas), intent(inout) :: amad
      real(rp), intent(inout) :: a, alat, vol
      ! Local variables
      real(rp), dimension(3) :: dtau
      real(rp) :: dl
      integer :: ibas, jbas, m, i, j
      DO IBAS = 1, NBAS
         DO JBAS = 1, NBAS
         DO M = 1, 3
            DTAU(M) = TAU(M, IBAS) - TAU(M, JBAS)
         END DO
         CALL SHORTN(DTAU, DTAU, DLAT, NKD)
         CALL STRX00(DTAU, A, ALAT, VOL, RLAT, NKR, DLAT, NKD, DL)
         WRITE (10, 995) IBAS, JBAS, DL
         AMAD(JBAS, IBAS) = DL
         END DO
      END DO
      DO I = 1, NBAS
         WRITE (11) (AMAD(I, J), J=1, NBAS)
      END DO
      CLOSE (11)
      RETURN
995   FORMAT(' IBAS, JBAS= ', 2I5, '  AMAD=', F12.6)
   end subroutine madmat

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Sets up the real and reciprocal space lattice vectors
   !
   !> @param[inout] AS
   !> @param[inout] TOL
   !> @param[inout] ALAT
   !> @param[inout] ALAT0
   !> @param[inout] RB0
   !> @param[inout] G1
   !> @param[inout] G2
   !> @param[inout] G3
   !> @param[inout] GT
   !> @param[inout] RB
   !> @param[inout] QB
   !> @param[inout] LMAX
   !> @param[inout] VOL
   !> @param[inout] AWALD
   !> @param[inout] DLAT
   !> @param[inout] NKD
   !> @param[inout] RLAT
   !> @param[inout] NKR
   !> @param[inout] NKDMX
   !> @param[inout] NKRMX
   !> @param[inout] WORK
   !> @return type(calculation)
   !---------------------------------------------------------------------------
   subroutine lattc(AS, TOL, ALAT, ALAT0, RB0, G1, G2, G3, GT, RB, QB, &
                    LMAX, VOL, AWALD, DLAT, NKD, RLAT, NKR, NKDMX, NKRMX, WORK)
      ! Input and output
      integer, intent(inout) :: lmax, nkd, nkr, nkdmx, nkrmx
      real(rp), dimension(3, nkrmx), intent(inout) :: rlat, dlat
      real(rp), dimension(3, 3), intent(inout) :: rb0, rb, qb
      real(rp), dimension(1), intent(inout) :: work
      real(rp) :: tol, alat, alat0, g1, g2, g3, gt, vol, as, awald
      ! Local variables
      real(rp), dimension(3, 3) :: qb0
      real(rp) :: tpiba, vol0, dstx, rdist0, qdist0, radd, qadd, alat1, a0, tol1, q0, r0
      integer :: nkdest, nkrest, iv, k, m
      TPIBA = 8.D0*DATAN(1.D0)/ALAT

      qb0(1, :) = cross_product(rb0(2, :), rb0(3, :))
      qb0(2, :) = cross_product(rb0(3, :), rb0(1, :))
      qb0(3, :) = cross_product(rb0(1, :), rb0(2, :))

      VOL0 = DABS(TRIPL(RB0(:, 1), RB0(:, 2), RB0(:, 3)))
      DO M = 1, 3
         DO K = 1, 3
            QB0(M, K) = QB0(M, K)*(1.D0/VOL0)
         END DO
      END DO
      DO K = 1, 3
         CALL RDIST(RB0(:, K), RB(:, K), G1, G2, G3, GT)
         CALL QDIST(QB0(:, K), QB(:, K), G1, G2, G3, GT)
      END DO
      VOL = TRIPL(RB(:, 1), RB(:, 2), RB(:, 3))
      VOL = DABS(VOL)*(ALAT**3)
      WRITE (10, 351)
351   FORMAT(/15X, 'RB0', 30X, 'QB0')
      WRITE (10, 350) ((RB0(M, K), M=1, 3), (QB0(M, K), M=1, 3), K=1, 3)
350   FORMAT(3F10.5, 5X, 3F10.5)
      DSTX = DMAX1(DABS(G1 - 1D0), DABS(G2 - 1D0), DABS(G3 - 1D0), DABS(GT - 1D0))
      IF (DSTX .GT. 1D-5) THEN
         WRITE (10, 451) G1, G2, G3, GT
451      FORMAT(' DISTORTED WITH ', 4F12.7, ':')
         WRITE (10, 350) ((RB(M, K), M=1, 3), (QB(M, K), M=1, 3), K=1, 3)
         CALL STRAIN(G1, G2, G3, GT)
      END IF
      WRITE (10, 998) VOL, VOL0*(ALAT**3)
998   FORMAT(' CELL VOLUME=', F12.6, '     VOL0=', F12.6)
! ------ SET UP REAL AND RECIP VECTORS ----
      RDIST0 = VOL0**(1.D0/3.D0)
      QDIST0 = 1.D0/RDIST0
      RADD = .7*RDIST0
      QADD = .7*QDIST0
      A0 = AS/RDIST0
      AWALD = A0/ALAT
      ALAT1 = ALAT0
      IF (ALAT1 .LE. 0.5D0) ALAT1 = ALAT
      IF (DABS(ALAT1/ALAT - 1.D0) .GT. 0.04D0) WRITE (10, 560)
560   FORMAT(/' *** WARNING: ALAT AND ALAT0 DEVIATE BY MORE THAN 4 %'/)
      TOL1 = TOL*ALAT1**(LMAX + 1)
      CALL LCTOFF(A0, VOL0, LMAX, TOL1, R0, Q0)
      NKDEST = 4.18879*(R0 + RADD)**3/VOL0 + .5
      NKREST = 4.18879*(Q0 + QADD)**3*VOL0 + .5
      !write(*, *)radd, r0
      WRITE (10, 340) AS, TOL, LMAX, AWALD, VOL0, ALAT1, NKDEST, NKREST
340   FORMAT(/' LATTC:  AS=', F6.3, '   TOL=', 1P, E8.2, '   LMAX=', I1,      &
                                                           &  '   AWALD=', 0P, F7.4, '   V0=', F8.5/' ALAT1=', F9.5,               &
                                                                   &  '   ESTIMATES:   NKD', I6, '   NKD', I6)
      CALL LGEN(RB0, R0 + RADD, NKD, NKDMX, DLAT, WORK)
      WRITE (10, 342) R0, R0*ALAT, RADD, NKD
342   FORMAT('  R0=', F9.4, '   RC=', F9.4, '   RADD=', F9.4, '   NKD=', I7)
      CALL LGEN(QB0, Q0 + QADD, NKR, NKRMX, RLAT, WORK)
      WRITE (10, 341) Q0, Q0*TPIBA, QADD, NKR
      DO IV = 1, NKD
         CALL RDIST(DLAT(1, IV), DLAT(1, IV), G1, G2, G3, GT)
      END DO
      DO IV = 1, NKR
         CALL QDIST(RLAT(1, IV), RLAT(1, IV), G1, G2, G3, GT)
      END DO
      RETURN
341   FORMAT('  Q0=', F9.4, '   QC=', F9.4, '   QADD=', F9.4, '   NKD=', I7)
   end subroutine lattc

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Widget to make structure constant dl for L=0, E=0, K=0.
   !
   !> @param[inout] TAU
   !> @param[inout] A
   !> @param[inout] ALAT
   !> @param[inout] VOL
   !> @param[inout] RLAT
   !> @param[inout] NKR
   !> @param[inout] DLAT
   !> @param[inout] NKD
   !> @param[inout] DL
   !---------------------------------------------------------------------------
   subroutine strx00(TAU, A, ALAT, VOL, RLAT, NKR, DLAT, NKD, DL)
      ! Input and output
      integer, intent(inout) :: nkr, nkd
      real(rp), dimension(3), intent(inout) :: tau
      real(rp), intent(inout) :: a, alat, vol, dl
      real(rp), dimension(3, nkr), intent(inout) :: rlat
      real(rp), dimension(3, nkd), intent(inout) :: dlat
      ! Local variables
      real(rp) :: tpi, tpiba, gamma, r2, scalp, r1
      integer :: ir, ir1

      TPI = 2.D0*PI
      GAMMA = 0.25D0/(A*A)
      TPIBA = TPI/ALAT
      DL = -GAMMA
      DO IR = 2, NKR
         R2 = TPIBA*TPIBA*(RLAT(1, IR)**2 + RLAT(2, IR)**2 + RLAT(3, IR)**2)
         SCALP = TPI*(RLAT(1, IR)*TAU(1) + RLAT(2, IR)*TAU(2) + RLAT(3, IR)*TAU(3))
         DL = DL + DCOS(SCALP)*DEXP(-GAMMA*R2)/R2
      END DO
      DL = DL*4.D0*PI/VOL
      IR1 = 2
      IF (TAU(1)**2 + TAU(2)**2 + TAU(3)**2 .GT. 1D-6) IR1 = 1
      DO IR = IR1, NKD
         R1 = ALAT*DSQRT((TAU(1) - DLAT(1, IR))**2 + (TAU(2) - DLAT(2, IR))**2      &
         &   + (TAU(3) - DLAT(3, IR))**2)
         DL = DL + DERFC(A*R1)/R1
      END DO
      IF (IR1 .EQ. 2) DL = DL - 2.D0*A/DSQRT(PI)
      RETURN
   end subroutine strx00

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Returns P1 = shortest vector such that P1-P is a lattice vector.
   !> a slightly skewed norm is used to make result unique.
   !> the first vector in the list must be the zero vector.
   !
   !> @param[inout] P
   !> @param[inout] P1
   !> @param[inout] DLAT
   !> @param[inout] NKD
   !---------------------------------------------------------------------------
   subroutine shortn(P, P1, DLAT, NKD)
      ! Input and output
      integer, intent(inout) :: nkd
      real(rp), dimension(3), intent(inout) :: p, p1
      real(rp), dimension(3, nkd), intent(inout) :: dlat
      ! Local variables
      real(rp) :: x, y, z, p2, critk0, dd, crit, anrm2
      integer :: irep, k0, k
      ANRM2(X, Y, Z) = X*X*1.00001D0 + Y*Y*1.00002D0 + Z*Z*1.00003D0            &
     &  - X*0.000004D0 - Y*0.000003D0 - Z*0.000002D0
      P1(1) = P(1)
      P1(2) = P(2)
      P1(3) = P(3)
      do irep = 1, 20
         p2 = anrm2(p1(1), p1(2), p1(3))
         k0 = 1
         critk0 = 1.d20
         do k = 1, nkd
            dd = dlat(1, k)**2 + dlat(2, k)**2 + dlat(3, k)**2
            if (dd > p2*4.d0) exit
            crit = anrm2(p1(1) + dlat(1, k), p1(2) + dlat(2, k), p1(3) + dlat(3, k))
            if (crit < critk0) then
               k0 = k
               critk0 = crit
            end if
         end do
         if (k0 == 1) exit
         p1(1) = p1(1) + dlat(1, k0)
         p1(2) = p1(2) + dlat(2, k0)
         p1(3) = p1(3) + dlat(3, k0)
      end do
      WRITE (10, *) '*** SHORTN: SHORTEST VECTOR NOT FOUND'
      RETURN
   end subroutine shortn

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Makes limits R0, Q0 for sums in real and recip space for a lattice
   !> with lattice constant 1.
   !
   !> @param[inout] A0
   !> @param[inout] V0
   !> @param[inout] LMAX
   !> @param[inout] TOL
   !> @param[inout] R0
   !> @param[inout] Q0
   !---------------------------------------------------------------------------
   subroutine lctoff(A0, V0, LMAX, TOL, R0, Q0)
      ! Input and output
      real(rp), intent(inout) :: a0, v0, tol, r0, q0
      integer, intent(inout) :: lmax
      ! Local variables
      real(rp), dimension(0:10) :: f, g
      real(rp) :: pi, q1, q2, gq0, gq1, r1, r2, try
      integer :: i
      PI = 4.D0*DATAN(1.D0)
      Q1 = 0.001D0
      IF (LMAX .GT. 2) Q1 = DSQRT(.5D0*(LMAX - 2))*A0/PI
      GQ1 = (2D0*PI*Q1)**(LMAX - 2)*DEXP(-(PI*Q1/A0)**2)*4D0*PI/V0
      IF (TOL .GT. GQ1) WRITE (10, *) '**** LCTOFF: TOL GT GQ1'
      Q2 = 50.D0
      Q0 = 5.D0
      DO I = 1, 25
         GQ0 = (2D0*PI*Q0)**(LMAX - 2)*DEXP(-(PI*Q0/A0)**2)*4D0*PI/V0
         IF (GQ0 .GT. TOL) THEN
            Q1 = Q0
         ELSE
            Q2 = Q0
         END IF
         Q0 = .5D0*(Q1 + Q2)
      END DO
! ---------------------------------------
      R1 = 0.1D0
      R2 = 50.D0
      R0 = 5.D0
      DO I = 1, 25
         CALL DLMTOR(R0, A0, LMAX, F, G)
         IF (F(LMAX) > TOL) R1 = R0
         IF (F(LMAX) <= TOL) R2 = R0
         R0 = 0.5D0*(R1 + R2)
      END DO
      TRY = (2D0*PI*Q0)**(LMAX - 2)*DEXP(-(PI*Q0/A0)**2)*4D0*PI/V0
      WRITE (10, 957) Q0, TRY, R0, F(LMAX)
957   FORMAT(' LCUT: Q0=', F12.6, '   TRY=', F12.6, '   R0=', F12.6,         &
                                                                                                                  &  '   F=', F12.6)
   end subroutine lctoff

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Radial part of damped lmtos F and FBAR, L=0 TO LMAX
   !
   !> @param[inout] R
   !> @param[inout] A
   !> @param[inout] LMAX
   !> @param[inout] F
   !> @param[inout] FBAR
   !---------------------------------------------------------------------------
   subroutine dlmtor(R, A, LMAX, F, FBAR)
      ! Input and output
      real(rp), intent(inout) :: r, a
      integer, intent(inout) :: lmax
      real(rp), dimension(0:lmax), intent(inout) :: f, fbar
      ! Local variables
      real(rp) :: obsrpi, obsqpi, z, emz2, erfc0, erfc1, erfc2, g, flm2, ta2r
      integer :: l

      OBSRPI = 0.564189835D0
      OBSQPI = 1.D0/DSQRT(PI)
      Z = A*R
      EMZ2 = DEXP(-Z*Z)
      ERFC0 = DERFC(Z)
      ERFC1 = -Z*ERFC0 + OBSRPI*EMZ2
      ERFC2 = -0.5D0*Z*ERFC1 + 0.25D0*ERFC0
      F(0) = ERFC0/R
      FBAR(0) = -ERFC2/(A*A*R)
      TA2R = 2.D0*A*A*R
      G = 2.D0*A*EMZ2*OBSRPI/R
      FLM2 = OBSRPI*EMZ2/Z - ERFC0
      DO L = 1, LMAX
         F(L) = ((L + L - 1)/R)*F(L - 1) + G
         FBAR(L) = ((L + L - 1)/R)*FBAR(L - 1) - FLM2
         FLM2 = F(L - 1)
         G = G*TA2R
      END DO
      RETURN
   end subroutine dlmtor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !>  Printout of strains E1...E6
   !
   !> @param[inout] GX
   !> @param[inout] GY
   !> @param[inout] GZ
   !> @param[inout] GT
   !---------------------------------------------------------------------------
   subroutine strain(GX, GY, GZ, GT)
      ! Input and output
      real(rp), intent(inout) :: gx, gy, gz, gt
      ! Local variables
      real(rp), dimension(3, 3) :: e, eps
      integer :: ixyz, m
      DATA E/1D0, 0D0, 0D0, 0D0, 1D0, 0D0, 0D0, 0D0, 1D0/

      DO IXYZ = 1, 3
         CALL RDIST(E(1, IXYZ), EPS(1, IXYZ), GX, GY, GZ, GT)
         DO M = 1, 3
            EPS(M, IXYZ) = EPS(M, IXYZ) - E(M, IXYZ)
         END DO
      END DO
      WRITE (10, 230) EPS(1, 1), EPS(2, 2), EPS(3, 3), EPS(2, 3) + EPS(3, 2),        &
         &  EPS(1, 3) + EPS(3, 1), EPS(1, 2) + EPS(2, 1)
230   FORMAT(/' STRAINS E1..E6:', 6F10.6)
      RETURN
   end subroutine strain

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !>  Generates lattice vectors.
   !
   !> @param[inout] BAS
   !> @param[in] BMAX
   !> @param[inout] NV
   !> @param[inout] NVMAX
   !> @param[inout] VECS
   !> @param[inout] WORK
   !> @return type(calculation)
   !---------------------------------------------------------------------------
   subroutine lgen(BAS, BMAX, NV, NVMAX, VECS, WORK)
      ! Input and output
      integer, intent(inout) :: nv, nvmax
      real(rp), dimension(3, 3), intent(inout) :: bas
      real(rp), dimension(3, nvmax), intent(inout) :: vecs
      real(rp), dimension(nvmax), intent(inout) :: work
      real(rp), intent(in) :: bmax
      ! Local variables
      real(rp), dimension(3) :: v
      real(rp) :: bmax2, v2, vsm, alow, xx, ddd
      integer :: i, j, k, m, iv, ilow, jv, imax, jmax, kmax
      !write(*, *)nvmax
      CALL LATLIM(BAS, BMAX, IMAX, JMAX, KMAX)
      BMAX2 = BMAX*BMAX
      !write(*, *)nvmax
      NV = 0
      DO I = -IMAX, IMAX
         DO J = -JMAX, JMAX
            DO K = -KMAX, KMAX
               DO M = 1, 3
                  V(M) = I*BAS(M, 1) + J*BAS(M, 2) + K*BAS(M, 3)
               END DO
               V2 = V(1)*V(1) + V(2)*V(2) + V(3)*V(3)
               IF (V2 .GT. BMAX2) CYCLE
               NV = NV + 1
               IF (NV .GT. NVMAX) THEN
                  WRITE (10, 633) NVMAX, I, IMAX
                  STOP '*** ERROR IN LGEN ***'
               END IF
               DO M = 1, 3
                  VECS(M, NV) = V(M)
               END DO
               VSM = DABS(V(1)) + DABS(V(2)) + DABS(V(3))
               WORK(NV) = V2 + VSM/1000.
            END DO
         END DO
      END DO
! --- SORT BY LENGTH -----------
      DO IV = 1, NV
         ILOW = IV
         ALOW = WORK(IV)
         DO JV = IV, NV
            IF (WORK(JV) < ALOW) THEN
               ALOW = WORK(JV)
               ILOW = JV
            END IF
         END DO
         IF (ILOW == IV) CYCLE
         DO M = 1, 3
            XX = VECS(M, IV)
            VECS(M, IV) = VECS(M, ILOW)
            VECS(M, ILOW) = XX
         END DO
         WORK(ILOW) = WORK(IV)
         XX = WORK(ILOW)
         !|      WRITE(6, 300) IV, (VECS(M, IV), M=1, 3), XX
         !|300   FORMAT(I6, 3X, 3F9.4, F12.4)
      END DO
! ---- PRINT A WARNING IF A BASIS VEC IS NOT IN VECTOR LIST ----
      DO 40 K = 1, 3
      DO 41 IV = 1, NV
         DDD = (BAS(1, K) - VECS(1, IV))**2 + (BAS(2, K) - VECS(2, IV))**2             &
                                                      &   + (BAS(3, K) - VECS(3, IV))**2
!     write(*, *)vecs(1, IV), vecs(2, IV), vecs(3, IV)
         IF (DDD .LT. 1.D-8) GOTO 42
41    END DO
      WRITE (10, 650) K
650   FORMAT(//' **** WARNING FROM LGEN **** BASIS VECTOR', I3,          &
                                                                         & '  NOT IN LIST'/' **** SUBROUTINE SHORTN MIGHT NOT WORK')
42    CONTINUE
40    END DO
      RETURN
633   FORMAT(/' NV=', I6, '  EXCEEDED,   I=', I3, '  IMAX=', I3)
   end subroutine lgen

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Limits in x y z direction, also initialize matrix
   !
   !> @param[in] BAS
   !> @param[in] VMAX
   !> @param[out] I1
   !> @param[out] I2
   !> @param[out] I3
   !---------------------------------------------------------------------------
   subroutine latlim(BAS, VMAX, I1, I2, I3)
      ! Input
      real(rp), dimension(3, 3), intent(in) :: bas
      real(rp), intent(in) :: vmax
      ! Output
      integer, intent(out) :: i1, i2, i3
      ! Local variables
      real(rp), dimension(3, 3) :: a
      real(rp) :: det, A11, A22, A33, A12, A13, A23
      integer :: i, j
      DO I = 1, 3
         DO J = I, 3
            A(I, J) = BAS(1, I)*BAS(1, J) + BAS(2, I)*BAS(2, J) + BAS(3, I)*BAS(3, J)
         END DO
      END DO
      A11 = A(1, 1)
      A12 = A(1, 2)
      A13 = A(1, 3)
      A22 = A(2, 2)
      A23 = A(2, 3)
      A33 = A(3, 3)
      DET = A11*A22*A33 + A12*A23*A13                                       &
      &   + A12*A23*A13 - A13*A22*A13                                       &
                                      &   - A23*A23*A11 - A12*A12*A33
      I1 = VMAX*DSQRT((A22*A33 - A23**2)/DET)
      I2 = VMAX*DSQRT((A11*A33 - A13**2)/DET)
      I3 = VMAX*DSQRT((A11*A22 - A12**2)/DET)
      RETURN
   end subroutine latlim
!# QDIST FORTRAN *

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Distorts a vector in recip space
   !
   !> @param[in] Q
   !> @param[inout] QOUT
   !> @param[in] GAMX
   !> @param[in] GAMY
   !> @param[in] GAMZ
   !> @param[in] GAMT
   !> @return type(calculation)
   !---------------------------------------------------------------------------
   subroutine qdist(Q, QOUT, GAMX, GAMY, GAMZ, GAMT)
      ! Input
      real(rp), intent(in) :: gamx, gamy, gamz, gamt
      real(rp), dimension(3), intent(in) :: q
      ! Output
      real(rp), dimension(3), intent(inout) :: qout
      ! Local variables
      real(rp) :: add

      ADD = -(Q(1) + Q(2) + Q(3))*(GAMT - 1.D0)/GAMT/3.D0
      QOUT(1) = (Q(1) + ADD)/GAMX
      QOUT(2) = (Q(2) + ADD)/GAMY
      QOUT(3) = (Q(3) + ADD)/GAMZ
      RETURN
   end subroutine qdist
!# RDIST FORTRAN *

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Distorts a vector in real space
   !
   !> @param[in] R
   !> @param[inout] ROUT
   !> @param[in] GAMX
   !> @param[in] GAMY
   !> @param[in] GAMZ
   !> @param[in] GAMT
   !> @return type(calculation)
   !---------------------------------------------------------------------------
   subroutine rdist(R, ROUT, GAMX, GAMY, GAMZ, GAMT)
      ! Input
      real(rp), intent(in) :: gamx, gamy, gamz, gamt
      real(rp), dimension(3), intent(in) :: r
      ! Output
      real(rp), dimension(3), intent(inout) :: rout
      ! Local variables
      real(rp) :: add, rx, ry, rz

      RX = R(1)*GAMX
      RY = R(2)*GAMY
      RZ = R(3)*GAMZ
      ADD = (RX + RY + RZ)*(GAMT - 1.D0)/3.D0
      ROUT(1) = RX + ADD
      ROUT(2) = RY + ADD
      ROUT(3) = RZ + ADD
      RETURN
   end subroutine rdist

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Distort NA real-space vectors in array A1 into array A2
   !
   !> @param[in] A1
   !> @param[out] A2
   !> @param[in] NA
   !> @param[in] GX
   !> @param[in] GY
   !> @param[in] GZ
   !> @param[in] GT
   !---------------------------------------------------------------------------
   subroutine rdistn(A1, A2, NA, GX, GY, GZ, GT)
      ! Input
      integer, intent(in) :: na
      real(rp), intent(in) :: gx, gy, gz, gt
      real(rp), dimension(3, na), intent(in) :: a1
      ! Output
      real(rp), dimension(3, na), intent(out) :: a2
      ! Local variables
      integer :: ia

      DO IA = 1, NA
         CALL RDIST(A1(:, IA), A2(:, IA), GX, GY, GZ, GT)
      END DO

      RETURN
   end subroutine rdistn

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TODO
   !---------------------------------------------------------------------------
   subroutine wkinit(this)
      class(charge), intent(inout) :: this
      character(len=1), dimension(60) :: string, str
      integer :: j70, ndefmx, j7max, ndfdmx, dfdsum, ipr, j7free, itot, limit, nsize, &
                 iprint, length, leng, j7name, imod, ndefd, lrest, idfd, iumx, iunw, &
                 j7x, j7y, i, j7ploc, nstr
      save
! ----- DEFINE STORAGE SIZE ------
!  START OF FIRST ARRAY AND MAX NUMBER TO BE DEFINED:
      nsize = 3000000
      J70 = 5
      NDEFMX = 100
      LIMIT = NSIZE
      J7MAX = 0
      NDFDMX = 0
      DFDSUM = 0.
      IPR = 0
      J7FREE = 5
      this%W(J7FREE - 1) = J7FREE
      ITOT = LIMIT*0.001 + 0.5
      WRITE (10, 399) ITOT
399   FORMAT(/' WKINIT:   size=', I8, ' K')
      RETURN
! ------ SUB TO SWITCH OUTPUT OF ARRAY DEF INFO ----
      ENTRY WKPRNT(IPRINT)
      IPR = IPRINT
      RETURN
! ------ SUBROUTINES TO DEFINE ARRAYS OF VARIOUS TYPES -----
      ENTRY DEFI(J7NAME, LENG)
      LENGTH = LENG
      GOTO 10
      ENTRY DEFR(J7NAME, LENG)
      LENGTH = LENG
      GOTO 10
      ENTRY DEFC(J7NAME, LENG)
      LENGTH = LENG*2
      GOTO 10
      ENTRY DEFRR(this, J7NAME, LENG)
      LENGTH = LENG*2
      GOTO 10
      ENTRY DEFCC(J7NAME, LENG)
      LENGTH = LENG*4
10    IF (LENGTH .LT. 0) STOP '**** LENGTH OF ARRAY NEGATIVE'
      IF (LENGTH .EQ. 0) LENGTH = 1
      IMOD = 1
      GOTO 83
84    J7NAME = J7FREE
      J7FREE = J7FREE + LENGTH + 1
      J7FREE = 4*((J7FREE + 2)/4) + 1
      J7MAX = MAX0(J7MAX, J7FREE)
      this%W(J7NAME - 1) = J7FREE
      NDEFD = NDEFD + 1
      NDFDMX = MAX0(NDFDMX, NDEFD)
      DFDSUM = DFDSUM + LENGTH*0.001
      IF (IPR .GT. 0) WRITE (10, 100) NDEFD, LENG, LENGTH, J7NAME, J7FREE - 1
100   FORMAT(' define array', I4, ':   els=', I8, '   length=', I8, ', ',      &
                                                                                                                 &   I8, '  to', I8)
      IF (J7FREE .LE. LIMIT) RETURN
      WRITE (10, 101) J7FREE
101   FORMAT(' **** WORKSPACE OVERFLOW,  NEED AT LEAST', I8)
      STOP
! ------- RETURN NUMBER OF WORDS LEFT ------
      ENTRY DEFASK(LREST)
      LREST = LIMIT - J7FREE - 2
      IF (IPR .GT. 0) WRITE (10, *) 'SPACE LEFT=', LREST, '  SINGLE WORDS'
      RETURN
! ------ RELEASE ---------------
      ENTRY RLSE(this, J7NAME)
      IF (J7NAME .GT. LIMIT) STOP '**** RESET POINTER GT LIMIT'
      IF (J7NAME .LT. 3) STOP '**** RESET POINTER LT 3'
      J7FREE = J7NAME
      IF (IPR .GT. 0) WRITE (10, *) 'J7FREE reset to', J7FREE
      RETURN
! ----- OUTPUT WORKSPACE INFO -----
      ENTRY WKINFO
      IMOD = 2
      GOTO 83
81    ITOT = LIMIT*0.001 + 0.5
      IDFD = DFDSUM + 0.5
      IUMX = (J7MAX - 1)*0.001 + 0.5
      IUNW = (J7FREE - 1)*0.001 + 0.5
      WRITE (10, 601) ITOT, IDFD, IUMX, IUNW, NDFDMX, NDEFD
601   FORMAT(/'  total workspace size=', I6, ' K'                         &
                                                              &  /'  total space defined =', I6, ' K'                              &
                                                            &  /'  workspace used:   max', I6, ' K   now', I6, ' K'                &
                                                                                    &  /'  arrays defined:   max', I8, '   now', I8)
      IF (J7FREE .EQ. J70) RETURN
      WRITE (10, 602)
602   FORMAT(/'  array', 6X, 'begin', 10X, 'end', 7X, 'length')
      J7X = J70
      DO I = 1, NDEFMX
         J7Y = this%W(J7X - 1)
         IF (J7X .EQ. J7FREE) RETURN
         WRITE (10, 540) I, J7X, J7Y - 1, J7Y - J7X
540      FORMAT(I5, 3I13)
         IF (J7Y .LT. J70 .OR. J7Y .GT. LIMIT) WRITE (10, *) '   . . . . . '
         IF (J7Y .LT. J70 .OR. J7Y .GT. LIMIT) RETURN
         J7X = J7Y
      END DO
      RETURN
      ENTRY WKCHK(STRING)
      IMOD = 0
      DO I = 1, 60
         STR(I) = STRING(I)
         NSTR = I
         IF (STRING(I) == 'J7') THEN
            EXIT
         END IF
      END DO
89    STR(NSTR) = '>'
      WRITE (10, 891) (STR(I), I=1, NSTR)
891   FORMAT('     WKCHK   <', 60A1)
83    NDEFD = 0
      J7X = J70
      J7PLOC = -999
      DO I = 1, NDEFMX
         IF (J7X < J70 .OR. J7X > LIMIT) THEN
            WRITE (10, 888) NDEFD, J7X, J7PLOC
888         FORMAT(' *** LINK DESTROYED AT START OF ARRAY', I3,             &
                                                                                                   &     ',  PTR=', I8, ' AT', I8)
            STOP
         END IF
         IF (J7X == J7FREE) GOTO 86
         NDEFD = NDEFD + 1
         J7PLOC = J7X - 1
         J7X = this%W(J7PLOC)
86    END DO
      IF (IMOD .EQ. 1) GOTO 84
      IF (IMOD .EQ. 2) GOTO 81
      WRITE (10, 360) NDEFD, J7FREE - 1
360   FORMAT('     LINKS OK   NDEFD=', I3, '   SPACE USED=', I7)
      RETURN
   end subroutine wkinit

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print class members values in namelist format
   !>
   !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
   !> @param[in] unit File unit used to write namelist
   !> @param[in] file File name used to write namelist
   !---------------------------------------------------------------------------
   subroutine print_state_full(this, unit, file)
      implicit none

      class(charge), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file
      integer :: newunit

      include 'include_codes/namelists/charge.f90'

      ! scalar

      gx = this%gx
      gy = this%gy
      gz = this%gz
      gt = this%gt
      a = this%a
      b = this%b
      c = this%c
      amax = this%amax
      bmax = this%bmax
      alamda = this%alamda
      rmax = this%rmax
      gmax = this%gmax
      ar2d = this%ar2d
      sws = this%sws
      vol = this%vol
      nq3 = this%nq3
      nr0 = this%nr0
      numr = this%numr
      numg = this%numg
      numvr = this%numvr
      numvg = this%numvg

      ! 1d allocatable

      if (allocated(this%w)) then
         allocate (w, mold=this%w)
         w = this%w
      else
         allocate (w(0))
      end if
      if (allocated(this%wssurf)) then
         allocate (wssurf, mold=this%wssurf)
         wssurf = this%wssurf
      else
         allocate (wssurf(0))
      end if
      if (allocated(this%bsx)) then
         allocate (bsx, mold=this%bsx)
         bsx = this%bsx
      else
         allocate (bsx(0))
      end if
      if (allocated(this%bsy)) then
         allocate (bsy, mold=this%bsy)
         bsy = this%bsy
      else
         allocate (bsy(0))
      end if
      if (allocated(this%bsz)) then
         allocate (bsz, mold=this%bsz)
         bsz = this%bsz
      else
         allocate (bsz(0))
      end if
      if (allocated(this%bkx)) then
         allocate (bkx, mold=this%bkx)
         bkx = this%bkx
      else
         allocate (bkx(0))
      end if
      if (allocated(this%bky)) then
         allocate (bky, mold=this%bky)
         bky = this%bky
      else
         allocate (bky(0))
      end if
      if (allocated(this%bkz)) then
         allocate (bkz, mold=this%bkz)
         bkz = this%bkz
      else
         allocate (bkz(0))
      end if
      if (allocated(this%qx3)) then
         allocate (qx3, mold=this%qx3)
         qx3 = this%qx3
      else
         allocate (qx3(0))
      end if
      if (allocated(this%qy3)) then
         allocate (qy3, mold=this%qy3)
         qy3 = this%qy3
      else
         allocate (qy3(0))
      end if
      if (allocated(this%qz3)) then
         allocate (qz3, mold=this%qz3)
         qz3 = this%qz3
      else
         allocate (qz3(0))
      end if
      if (allocated(this%qx)) then
         allocate (qx, mold=this%qx)
         qx = this%qx
      else
         allocate (qx(0))
      end if
      if (allocated(this%qy)) then
         allocate (qy, mold=this%qy)
         qy = this%qy
      else
         allocate (qy(0))
      end if
      if (allocated(this%qz)) then
         allocate (qz, mold=this%qz)
         qz = this%qz
      else
         allocate (qz(0))
      end if
      if (allocated(this%asx)) then
         allocate (asx, mold=this%asx)
         asx = this%asx
      else
         allocate (asx(0))
      end if
      if (allocated(this%asy)) then
         allocate (asy, mold=this%asy)
         asy = this%asy
      else
         allocate (asy(0))
      end if
      if (allocated(this%asz)) then
         allocate (asz, mold=this%asz)
         asz = this%asz
      else
         allocate (asz(0))
      end if
      if (allocated(this%akx)) then
         allocate (akx, mold=this%akx)
         akx = this%akx
      else
         allocate (akx(0))
      end if
      if (allocated(this%aky)) then
         allocate (aky, mold=this%aky)
         aky = this%aky
      else
         allocate (aky(0))
      end if
      if (allocated(this%akz)) then
         allocate (akz, mold=this%akz)
         akz = this%akz
      else
         allocate (akz(0))
      end if
      if (allocated(this%dr)) then
         allocate (dr, mold=this%dr)
         dr = this%dr
      else
         allocate (dr(0))
      end if
      if (allocated(this%dg)) then
         allocate (dg, mold=this%dg)
         dg = this%dg
      else
         allocate (dg(0))
      end if
      if (allocated(this%wsimp)) then
         allocate (wsimp, mold=this%wsimp)
         wsimp = this%wsimp
      else
         allocate (wsimp(0))
      end if

      if (allocated(this%dss)) then
         allocate (dss, mold=this%dss)
         dss = this%dss
      else
         allocate (dss(0, 0))
      end if
      if (allocated(this%dsz)) then
         allocate (dsz, mold=this%dsz)
         dsz = this%dsz
      else
         allocate (dsz(0, 0))
      end if
      if (allocated(this%ds3z2)) then
         allocate (ds3z2, mold=this%ds3z2)
         ds3z2 = this%ds3z2
      else
         allocate (ds3z2(0, 0))
      end if
      if (allocated(this%dsx2y2)) then
         allocate (dsx2y2, mold=this%dsx2y2)
         dsx2y2 = this%dsx2y2
      else
         allocate (dsx2y2(0, 0))
      end if
      if (allocated(this%dsxy)) then
         allocate (dsxy, mold=this%dsxy)
         dsxy = this%dsxy
      else
         allocate (dsxy(0, 0))
      end if
      if (allocated(this%dzz)) then
         allocate (dzz, mold=this%dzz)
         dzz = this%dzz
      else
         allocate (dzz(0, 0))
      end if
      if (allocated(this%dz3z2)) then
         allocate (dz3z2, mold=this%dz3z2)
         dz3z2 = this%dz3z2
      else
         allocate (dz3z2(0, 0))
      end if
      if (allocated(this%am)) then
         allocate (am, mold=this%am)
         am = this%am
      else
         allocate (am(0, 0))
      end if
      if (allocated(this%bm)) then
         allocate (bm, mold=this%bm)
         bm = this%bm
      else
         allocate (bm(0, 0))
      end if
      if (allocated(this%pm)) then
         allocate (pm, mold=this%pm)
         pm = this%pm
      else
         allocate (pm(0, 0))
      end if
      if (allocated(this%amad)) then
         allocate (amad, mold=this%amad)
         amad = this%amad
      else
         allocate (amad(0, 0))
      end if

      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         write (unit, nml=charge)
      else if (present(file)) then
         open (unit=newunit, file=file)
         write (newunit, nml=charge)
         close (newunit)
      else
         write (*, nml=charge)
      end if

   end subroutine print_state_full

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print class members values in namelist format
   !>
   !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
   !> @param[in] unit File unit used to write namelist
   !> @param[in] file File name used to write namelist
   !---------------------------------------------------------------------------
   subroutine print_state(this, unit, file)
      implicit none

      class(charge), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file
      integer :: newunit

      include 'include_codes/namelists/charge.f90'

      ! scalar

      gx = this%gx
      gy = this%gy
      gz = this%gz
      gt = this%gt

      ! 1d allocatable

      if (allocated(this%wssurf)) then
         allocate (wssurf, mold=this%wssurf)
         wssurf = this%wssurf
      else
         allocate (wssurf(0))
      end if

      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         write (unit, nml=charge)
      else if (present(file)) then
         open (unit=newunit, file=file)
         write (newunit, nml=charge)
         close (newunit)
      else
         write (*, nml=charge)
      end if

   end subroutine print_state

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print class members values in namelist format
   !>
   !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
   !> @param[in] unit File unit used to write namelist
   !> @param[in] file File name used to write namelist
   !---------------------------------------------------------------------------
   subroutine print_state_formatted(this, unit, file)
      implicit none

      class(charge), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file
      integer :: newunit

      type(namelist_generator) :: nml

      include 'include_codes/namelists/charge.f90'

      ! scalar

      nml = namelist_generator('charge')

      call nml%add('gx', this%gx)
      call nml%add('gy', this%gy)
      call nml%add('gz', this%gz)
      call nml%add('gt', this%gt)
      call nml%add('a', this%a)
      call nml%add('b', this%b)
      call nml%add('c', this%c)
      call nml%add('amax', this%amax)
      call nml%add('bmax', this%bmax)
      call nml%add('alamda', this%alamda)
      call nml%add('rmax', this%rmax)
      call nml%add('gmax', this%gmax)
      call nml%add('ar2d', this%ar2d)
      call nml%add('sws', this%sws)
      call nml%add('vol', this%vol)
      call nml%add('nq3', this%nq3)
      call nml%add('nr0', this%nr0)
      call nml%add('numr', this%numr)
      call nml%add('numg', this%numg)
      call nml%add('numvr', this%numvr)
      call nml%add('numvg', this%numvg)

      ! 1d allocatable
      ! TODO: implement test inside namelist_generator
      if (allocated(this%w)) call nml%add('w', this%w)
      if (allocated(this%wssurf)) call nml%add('wssurf', this%wssurf)
      if (allocated(this%bsx)) call nml%add('bsx', this%bsx)
      if (allocated(this%bsy)) call nml%add('bsy', this%bsy)
      if (allocated(this%bsz)) call nml%add('bsz', this%bsz)
      if (allocated(this%bkx)) call nml%add('bkx', this%bkx)
      if (allocated(this%bky)) call nml%add('bky', this%bky)
      if (allocated(this%bkz)) call nml%add('bkz', this%bkz)
      if (allocated(this%qx3)) call nml%add('qx3', this%qx3)
      if (allocated(this%qy3)) call nml%add('qy3', this%qy3)
      if (allocated(this%qz3)) call nml%add('qz3', this%qz3)
      if (allocated(this%qx)) call nml%add('qx', this%qx)
      if (allocated(this%qy)) call nml%add('qy', this%qy)
      if (allocated(this%qz)) call nml%add('qz', this%qz)
      if (allocated(this%asx)) call nml%add('asx', this%asx)
      if (allocated(this%asy)) call nml%add('asy', this%asy)
      if (allocated(this%asz)) call nml%add('asz', this%asz)
      if (allocated(this%akx)) call nml%add('akx', this%akx)
      if (allocated(this%aky)) call nml%add('aky', this%aky)
      if (allocated(this%akz)) call nml%add('akz', this%akz)
      if (allocated(this%dr)) call nml%add('dr', this%dr)
      if (allocated(this%dg)) call nml%add('dg', this%dg)
      if (allocated(this%wsimp)) call nml%add('wsimp', this%wsimp)

      ! 2d allocatable
      ! TODO: implement test inside namelist_generator
      if (allocated(this%dss)) call nml%add('dss', this%dss)
      if (allocated(this%dsz)) call nml%add('dsz', this%dsz)
      if (allocated(this%ds3z2)) call nml%add('ds3z2', this%ds3z2)
      if (allocated(this%dsx2y2)) call nml%add('dsx2y2', this%dsx2y2)
      if (allocated(this%dsxy)) call nml%add('dsxy', this%dsxy)
      if (allocated(this%dzz)) call nml%add('dzz', this%dzz)
      if (allocated(this%dz3z2)) call nml%add('dz3z2', this%dz3z2)
      if (allocated(this%am)) call nml%add('am', this%am)
      if (allocated(this%bm)) call nml%add('bm', this%bm)
      if (allocated(this%pm)) call nml%add('pm', this%pm)
      if (allocated(this%amad)) call nml%add('amad', this%amad)

      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         call nml%generate_namelist(unit=unit)
      else if (present(file)) then
         call nml%generate_namelist(file=file)
      else
         call nml%generate_namelist()
      end if
   end subroutine print_state_formatted

end module charge_mod
