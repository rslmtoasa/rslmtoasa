!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Hamiltonian
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
!> Module to handle procedures related to the Hamiltonian
!------------------------------------------------------------------------------

module hamiltonian_mod

   use control_mod
   use symbolic_atom_mod
   use element_mod
   use potential_mod
   use lattice_mod
   use charge_mod
   use precision_mod, only: rp
   use math_mod
   use string_mod
   use logger_mod, only: g_logger
   use timer_mod, only: g_timer
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   !> Module's main procedure
   type, public :: hamiltonian
      !> Charge
      class(charge), pointer :: charge
      !> Lattice
      class(lattice), pointer :: lattice
      !> Control
      class(control), pointer :: control

      !> Spin-orbit coupling Hamiltonian
      complex(rp), dimension(:, :, :), allocatable :: lsham
      !> Torque operator T=[o, Hso]
      complex(rp), dimension(:, :, :, :), allocatable :: tmat
      !> Bulk Hamiltonian
      complex(rp), dimension(:, :, :, :), allocatable :: ee, eeo
      !> Local Hamiltonian
      complex(rp), dimension(:, :, :, :), allocatable :: hall, hallo
      !> Hamiltonian built in chbar_nc (description to be improved)
      complex(rp), dimension(:, :, :, :), allocatable :: hmag
      !> Hamiltonian built in ham0m_nc (description to be improved
      complex(rp), dimension(:, :, :), allocatable :: hhmag
      !> Overlap Hamiltonian
      complex(rp), dimension(:, :, :), allocatable :: obarm
      !> Gravity center Hamiltonian
      complex(rp), dimension(:, :, :), allocatable :: enim
      !> Logical variable to include hoh term
      logical :: hoh
      !> Rotate Hamiltonian to local spin axis
      logical :: local_axis
      !> Add orbital polarization to Hamiltonian
      logical :: orb_pol
      !> Bulk Hamiltonian (backup for rotation)
      complex(rp), dimension(:, :, :, :), allocatable :: ee_glob, eeo_glob
      !> Local Hamiltonian (backup for rotation)
      complex(rp), dimension(:, :, :, :), allocatable :: hall_glob, hallo_glob
    !!> Spin-orbit coupling Hamiltonian (backup for rotation)
      !complex(rp), dimension(:, :, :), allocatable :: lsham
      !> Gravity center Hamiltonian (backup for rotation)
      complex(rp), dimension(:, :, :), allocatable :: enim_glob
   contains
      procedure :: build_lsham
      procedure :: build_bulkham
      procedure :: build_locham
      procedure :: build_obarm
      procedure :: build_enim
      procedure :: build_from_paoflow
      procedure :: build_from_paoflow_opt
      procedure :: torque_operator_collinear
      procedure :: rs2pao
      procedure :: chbar_nc
      procedure :: ham0m_nc
      procedure :: hmfind
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: rotate_to_local_axis
      procedure :: rotate_from_local_axis
   end type hamiltonian

   interface hamiltonian
      procedure :: constructor
   end interface

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] charge_obj Variable that holds charge_mod properties
   !> @return type(hamiltonian)
   !---------------------------------------------------------------------------
   function constructor(charge_obj) result(obj)
      type(hamiltonian) :: obj
      type(charge), target, intent(in) :: charge_obj

      obj%charge => charge_obj
      obj%lattice => charge_obj%lattice
      obj%control => charge_obj%lattice%control

      call obj%restore_to_default()
      call obj%build_from_file()
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   subroutine destructor(this)
      type(hamiltonian) :: this
#ifdef USE_SAFE_ALLOC
      if (allocated(this%lsham)) call g_safe_alloc%deallocate('hamiltonian.lsham', this%lsham)
      if (allocated(this%tmat)) call g_safe_alloc%deallocate('hamiltonian.tmat', this%tmat)
      if (allocated(this%ee)) call g_safe_alloc%deallocate('hamiltonian.ee', this%ee)
      if (allocated(this%hmag)) call g_safe_alloc%deallocate('hamiltonian.hmag', this%hmag)
      if (allocated(this%hhmag)) call g_safe_alloc%deallocate('hamiltonian.hhmag', this%hhmag)
      if (allocated(this%hall)) call g_safe_alloc%deallocate('hamiltonian.hall', this%hall)
      if (allocated(this%eeo)) call g_safe_alloc%deallocate('hamiltonian.eeo', this%eeo)
      if (allocated(this%hallo)) call g_safe_alloc%deallocate('hamiltonian.hallo', this%hallo)
      if (allocated(this%obarm)) call g_safe_alloc%deallocate('hamiltonian.obarm', this%obarm)
      if (allocated(this%enim)) call g_safe_alloc%deallocate('hamiltonian.enim', this%enim)
      if (allocated(this%ee_glob)) call g_safe_alloc%deallocate('hamiltonian.ee_glob', this%ee_glob)
      if (allocated(this%eeo_glob)) call g_safe_alloc%deallocate('hamiltonian.eeo_glob', this%eeo_glob)
      if (allocated(this%enim_glob)) call g_safe_alloc%deallocate('hamiltonian.enim_glob', this%enim_glob)
#else
      if (allocated(this%lsham)) deallocate (this%lsham)
      if (allocated(this%tmat)) deallocate (this%tmat)
      if (allocated(this%ee)) deallocate (this%ee)
      if (allocated(this%eeo)) deallocate (this%eeo)
      if (allocated(this%hmag)) deallocate (this%hmag)
      if (allocated(this%hhmag)) deallocate (this%hhmag)
      if (allocated(this%hall)) deallocate (this%hall)
      if (allocated(this%hallo)) deallocate (this%hallo)
      if (allocated(this%obarm)) deallocate (this%obarm)
      if (allocated(this%enim)) deallocate (this%enim)
      if (allocated(this%ee_glob)) deallocate (this%ee_glob)
      if (allocated(this%eeo_glob)) deallocate (this%eeo_glob)
      if (allocated(this%enim_glob)) deallocate (this%enim_glob)
#endif
   end subroutine destructor

   ! Member functions
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Read parameters from input file
   !---------------------------------------------------------------------------
   subroutine build_from_file(this)
      class(hamiltonian), intent(inout) :: this

      ! variables associated with the reading processes
      integer :: iostatus, funit

      include 'include_codes/namelists/hamiltonian.f90'

      hoh = this%hoh
      local_axis = this%local_axis
      orb_pol = this%orb_pol

      ! Reading
      open (newunit=funit, file=this%control%fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//trim(this%control%fname)//' not found', __FILE__, __LINE__)
      end if

      read (funit, nml=hamiltonian, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%error('Error while reading namelist', __FILE__, __LINE__)
         call g_logger%error('iostatus = '//int2str(iostatus), __FILE__, __LINE__)
      end if
      close (funit)

      this%hoh = hoh
      this%local_axis = local_axis
      this%orb_pol = orb_pol

   end subroutine build_from_file

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this)
      class(hamiltonian) :: this

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('hamiltonian.lsham', this%lsham, (/18, 18, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.tmat', this%tmat, (/18, 18, 3, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hhmag', this%hhmag, (/9, 9, 4/))
      call g_safe_alloc%allocate('hamiltonian.hmag', this%hmag, (/9, 9, this%charge%lattice%kk, 4/))
      call g_safe_alloc%allocate('hamiltonian.ee', this%ee, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hall', this%hall, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.hall_glob', this%hall_glob, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.ee_glob', this%ee_glob, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      if (hoh) then
         call g_safe_alloc%allocate('hamiltonian.eeo', this%eeo, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
         call g_safe_alloc%allocate('hamiltonian.hallo', this%hallo, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
         call g_safe_alloc%allocate('hamiltonian.obarm', this%obarm, (/18, 18, this%charge%lattice%ntype/))
         call g_safe_alloc%allocate('hamiltonian.enim', this%enim, (/18, 18, this%charge%lattice%ntype/))
      end if
      if (local_axis) then
         call g_safe_alloc%allocate('hamiltonian.hall_glob', this%hall_glob, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
         call g_safe_alloc%allocate('hamiltonian.ee_glob', this%ee_glob, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
         if (hoh) then
            call g_safe_alloc%allocate('hamiltonian.ee0_glob', this%eeo_glob, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
            call g_safe_alloc%allocate('hamiltonian.hallo_glob', this%hallo_glob, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
            call g_safe_alloc%allocate('hamiltonian.enim_glob', this%enim_glob, (/18, 18, this%charge%lattice%ntype/))
         end if
      end if
#else
      allocate (this%lsham(18, 18, this%charge%lattice%ntype))
      allocate (this%tmat(18, 18, 3, this%charge%lattice%ntype))
      allocate (this%hhmag(9, 9, 4), this%hmag(9, 9, this%charge%lattice%kk, 4))
      allocate (this%ee(18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype))
      allocate (this%hall(18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax))
      if (this%hoh) then
         allocate (this%eeo(18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype))
         allocate (this%hallo(18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax))
         allocate (this%obarm(18, 18, this%charge%lattice%ntype))
         allocate (this%enim(18, 18, this%charge%lattice%ntype))
      end if
      if (this%local_axis) then
         allocate (this%ee_glob(18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype))
         allocate (this%hall_glob(18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax))
         if (this%hoh) then
            allocate (this%eeo_glob(18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype))
            allocate (this%hallo_glob(18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax))
            allocate (this%enim_glob(18, 18, this%charge%lattice%ntype))
         end if
      end if
#endif

      this%lsham(:, :, :) = 0.0D0
      this%tmat(:, :, :, :) = 0.0D0
      this%hhmag(:, :, :) = 0.0D0
      this%hmag(:, :, :, :) = 0.0D0
      this%hall(:, :, :, :) = 0.0D0
      this%ee(:, :, :, :) = 0.0D0
      if (this%hoh) then
         this%hallo(:, :, :, :) = 0.0D0
         this%eeo(:, :, :, :) = 0.0D0
         this%obarm(:, :, :) = 0.0D0
         this%enim(:, :, :) = 0.0D0
      end if
      if (this%local_axis) then
         this%hall_glob(:, :, :, :) = 0.0D0
         this%ee_glob(:, :, :, :) = 0.0D0
         if (this%hoh) then
            this%hallo_glob(:, :, :, :) = 0.0D0
            this%eeo_glob(:, :, :, :) = 0.0D0
            this%enim_glob(:, :, :) = 0.0D0
         end if
      end if
      this%hoh = .false.
      this%local_axis = .false.
      this%orb_pol = .false.
   end subroutine restore_to_default

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Build the spin-orbit coupling hamiltonian according to Wu/Freeman PRB 54, 61 (1996)
   !---------------------------------------------------------------------------
   subroutine build_lsham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k
      complex(rp) :: prefac, sg
      real(rp) :: soc_p, soc_d
      complex(rp), dimension(2) :: rac
      complex(rp), dimension(9, 9) :: Lx, Ly, Lz
      real(rp) :: lz_loc
      !  Getting the angular momentum operators from the math_mod that are in cartesian coordinates
      Lx(:, :) = L_x(:, :)
      Ly(:, :) = L_y(:, :)
      Lz(:, :) = L_z(:, :)

      ! Transforming them into the spherical harmonics coordinates
      call hcpx(Lx, 'cart2sph')
      call hcpx(Ly, 'cart2sph')
      call hcpx(Lz, 'cart2sph')

      ! Writing the L.S hamiltonian
      this%lsham(:, :, :) = cmplx(0.0D0, 0.0D0)
      do k = 1, this%charge%lattice%ntype
         sg = cmplx(0.5D0, 0.0D0)
         soc_p = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_p(1) * this%charge%lattice%symbolic_atoms(k)%potential%xi_p(2))
         soc_d = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1) * this%charge%lattice%symbolic_atoms(k)%potential%xi_d(2))
         ! Check if orbital polarization is enabled
         if (this%orb_pol) then
            rac = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1) * this%charge%lattice%symbolic_atoms(k)%potential%rac)
            lz_loc = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1) * this%charge%lattice%symbolic_atoms(k)%potential%lmom(3))
         else
            rac = 0.0_RP
            lz_loc = 0.0_RP
         end if

         prefac = 0.0_RP
         do i = 1, 9
            do j = 1, 9
               if (i >= 2 .and. i <= 4 .and. j >= 2 .and. j <= 4) prefac = sg * soc_p
               if (i >= 5 .and. i <= 9 .and. j >= 5 .and. j <= 9) prefac = sg * soc_d
               this%lsham(j, i, k) = this%lsham(j, i, k) + prefac * Lz(j, i) + Lz(j, i) * rac(1) * lz_loc ! H11
               this%lsham(j, i + 9, k) = this%lsham(j, i + 9, k) + prefac * (Lx(j, i) - i_unit * Ly(j, i)) ! H12
               this%lsham(j + 9, i, k) = this%lsham(j + 9, i, k) + prefac * (Lx(j, i) + i_unit * Ly(j, i)) ! H21
               this%lsham(j + 9, i + 9, k) = this%lsham(j + 9, i + 9, k) - prefac * Lz(j, i) - Lz(j, i) * rac(2) * lz_loc ! H22
               !write(50,*) 'ntype=', k
               !write(51,*) 'ntype=', k
               !write(50,'(18f10.6)') real(this%lsham(:,:,k))
               !write(51,'(18f10.6)') aimag(this%lsham(:,:,k))
            end do
         end do
      end do
   end subroutine build_lsham

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Build the torque operator for the collinear case (SQA in z)
   !> based on the definition of PRM 2, 013801 (2018).
   !> Implemented by Ivan Miranda on 07/11/2023.
   !---------------------------------------------------------------------------
   subroutine torque_operator_collinear(this)
      !
      class(hamiltonian), intent(inout) :: this
      !
      ! Local variables
      integer :: i, j, k
      complex(rp) :: prefac, sg, soc_p, soc_d
      complex(rp), dimension(9, 9) :: Lx, Ly, Lz
      !  Getting the angular momentum operators from the math_mod that are in cartesian coordinates
      Lx(:, :) = L_x(:, :)
      Ly(:, :) = L_y(:, :)
      Lz(:, :) = L_z(:, :)

      ! Transforming them into the spherical harmonics coordinates
      call hcpx(Lx, 'cart2sph')
      call hcpx(Ly, 'cart2sph')
      call hcpx(Lz, 'cart2sph')

      ! Now write the torque operator matrix
      this%tmat(:, :, :, :) = cmplx(0.0_RP, 0.0_RP)
      do k = 1, this%charge%lattice%ntype
         sg = cmplx(0.5_RP, 0.0_RP)
         soc_p = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_p(1) * this%charge%lattice%symbolic_atoms(k)%potential%xi_p(2))
         soc_d = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1) * this%charge%lattice%symbolic_atoms(k)%potential%xi_d(2))
         prefac = 0.0_RP
         do i = 1, 9
            do j = 1, 9
               if (i >= 2 .and. i <= 4 .and. j >= 2 .and. j <= 4) prefac = sg * soc_p
               if (i >= 5 .and. i <= 9 .and. j >= 5 .and. j <= 9) prefac = sg * soc_d
               ! build Tx
               this%tmat(j, i, 1, k) = this%tmat(j, i, 1, k) + prefac * i_unit * Ly(j, i) * 2.0_RP ! Tx_11
               this%tmat(j, i + 9, 1, k) = this%tmat(j, i + 9, 1, k) - prefac * Lz(j, i) * 2.0_RP * cone ! Tx_12
               this%tmat(j + 9, i, 1, k) = this%tmat(j + 9, i, 1, k) + prefac * Lz(j, i) * 2.0_RP * cone ! Tx_21
               this%tmat(j + 9, i + 9, 1, k) = this%tmat(j + 9, i + 9, 1, k) - prefac * i_unit * Ly(j, i) * 2.0_RP ! Tx_22
               ! build Ty
               this%tmat(j, i, 2, k) = this%tmat(j, i, 2, k) - prefac * i_unit * Lx(j, i) * 2.0_RP ! Ty_11
               this%tmat(j, i + 9, 2, k) = this%tmat(j, i + 9, 2, k) + prefac * i_unit * Lz(j, i) * 2.0_RP ! Ty_12
               this%tmat(j + 9, i, 2, k) = this%tmat(j + 9, i, 2, k) + prefac * i_unit * Lz(j, i) * 2.0_RP ! Ty_21
               this%tmat(j + 9, i + 9, 2, k) = this%tmat(j + 9, i + 9, 2, k) + prefac * i_unit * Lx(j, i) * 2.0_RP ! Ty_22
               ! build Tz
               this%tmat(j, i + 9, 3, k) = this%tmat(j, i + 9, 3, k) + prefac * (Lx(j, i) - i_unit * Ly(j, i)) * 2.0_RP * cone ! Tz_12
               this%tmat(j + 9, i, 3, k) = this%tmat(j + 9, i, 3, k) + prefac * (Lx(j, i) + i_unit * Ly(j, i)) * 2.0_RP * cmone ! Tz_21
            end do
         end do
      end do

   end subroutine torque_operator_collinear

   subroutine build_obarm(this)
      implicit none
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      complex(rp), dimension(9, 9) :: obm0, obm1
      complex(rp), dimension(3) :: mom
      integer :: ntype ! Atom type index
      integer :: l, m ! Orbital index

      this%obarm = 0.D00

      do ntype = 1, this%lattice%ntype
         obm0 = cmplx(0.0D0); obm1 = cmplx(0.0D0)
         do m = 1, 9
            obm0(m, m) = this%lattice%symbolic_atoms(ntype)%potential%obx0(m)
            obm1(m, m) = this%lattice%symbolic_atoms(ntype)%potential%obx1(m)
         end do
         mom(:) = cmplx(this%lattice%symbolic_atoms(ntype)%potential%mom(:), 0.0D0)
         do m = 1, 9
            do l = 1, 9
               this%obarm(m, l, ntype) = obm0(m, l) + obm1(m, l) * mom(3)
               this%obarm(m + 9, l + 9, ntype) = obm0(m, l) - obm1(m, l) * mom(3)
               this%obarm(l, m + 9, ntype) = obm1(m, l) * mom(1) - i_unit * obm1(m, l) * mom(2)
               this%obarm(l + 9, m, ntype) = obm1(m, l) * mom(1) + i_unit * obm1(m, l) * mom(2)
            end do
         end do
         call hcpx(this%obarm(1:9, 1:9, ntype), 'cart2sph')
         call hcpx(this%obarm(10:18, 10:18, ntype), 'cart2sph')
         call hcpx(this%obarm(1:9, 10:18, ntype), 'cart2sph')
         call hcpx(this%obarm(10:18, 1:9, ntype), 'cart2sph')
      end do
   end subroutine build_obarm

   subroutine build_enim(this)
      implicit none
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      complex(rp), dimension(9, 9) :: em0, em1
      complex(rp), dimension(9) :: ex0, ex1
      complex(rp), dimension(3) :: mom
      complex(rp) :: eu, ed
      integer :: ntype ! Atom type index
      integer :: l, m ! Orbital index

      this%enim = 0.0D0

      do ntype = 1, this%lattice%ntype
         em0 = cmplx(0.0D0); em1 = cmplx(0.0D0)
         do m = 1, 9
            eu = this%lattice%symbolic_atoms(ntype)%potential%cx(m, 1) - this%lattice%symbolic_atoms(ntype)%potential%cex(m, 1)
            ed = this%lattice%symbolic_atoms(ntype)%potential%cx(m, 2) - this%lattice%symbolic_atoms(ntype)%potential%cex(m, 2)
            ex0(m) = 0.5 * (eu + ed)
            ex1(m) = 0.5 * (eu - ed)
            em0(m, m) = ex0(m)
            em1(m, m) = ex1(m)
         end do
         mom(:) = cmplx(this%lattice%symbolic_atoms(ntype)%potential%mom(:), 0.0D0)
         do m = 1, 9
            do l = 1, 9
               this%enim(m, l, ntype) = em0(m, l) + em1(m, l) * mom(3)
               this%enim(m + 9, l + 9, ntype) = em0(m, l) - em1(m, l) * mom(3)
               this%enim(l, m + 9, ntype) = em1(m, l) * mom(1) - i_unit * em1(m, l) * mom(2)
               this%enim(l + 9, m, ntype) = em1(m, l) * mom(1) + i_unit * em1(m, l) * mom(2)
            end do
         end do
         call hcpx(this%enim(1:9, 1:9, ntype), 'cart2sph')
         call hcpx(this%enim(10:18, 10:18, ntype), 'cart2sph')
         call hcpx(this%enim(1:9, 10:18, ntype), 'cart2sph')
         call hcpx(this%enim(10:18, 1:9, ntype), 'cart2sph')

         if (this%local_axis) then
            this%enim_glob = this%enim
         end if
      end do
   end subroutine build_enim

   subroutine build_bulkham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, m, ino, ja, ji, nr, ia
      integer :: ntype

      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype) ! Atom number in clust
         ino = this%charge%lattice%num(ia) ! Atom bravais type of ia
         nr = this%charge%lattice%nn(ia, 1) ! Number of neighbours considered
         !write(123, *)'bulkham'
         call this%chbar_nc(ia, nr, ino, ntype)
         do m = 1, nr
            do i = 1, 9
               do j = 1, 9
                  this%ee(j, i, m, ntype) = this%hmag(j, i, m, 4) + this%hmag(j, i, m, 3)        ! H0+Hz
                  this%ee(j + 9, i + 9, m, ntype) = this%hmag(j, i, m, 4) - this%hmag(j, i, m, 3)        ! H0-Hz
                  this%ee(j, i + 9, m, ntype) = this%hmag(j, i, m, 1) - i_unit * this%hmag(j, i, m, 2) ! Hx-iHy
                  this%ee(j + 9, i, m, ntype) = this%hmag(j, i, m, 1) + i_unit * this%hmag(j, i, m, 2) ! Hx+iHy
               end do ! end of orbital j loop
            end do ! end of orbital i loop
            !write(128, *) 'm=', m, 'ntype= ', ntype
            !write(128, '(18f10.6)') real(this%ee(:, :, m, ntype))
         end do ! end of neighbour number
         if (this%hoh) then
            call this%build_obarm()
            call this%build_enim()
            do m = 1, nr
               ji = 0
               if (m > 1) then
                  ja = this%charge%lattice%nn(ia, m)
                  if (ja .ne. 0) then
                     ji = this%charge%lattice%iz(ja)
                  end if
               else
                  ji = this%charge%lattice%iz(ia)
               end if
               ! Check if neighbour 'm' exists for atom 'ntype', otherwise fill HoH Hamiltonian with zeros.
               if (ji > 0) then
                  call zgemm('n', 'n', 18, 18, 18, cone, this%ee(:, :, m, ntype), 18, this%obarm(:, :, ji), 18, czero, this%eeo(:, :, m, ntype), 18)
               else
                  this%eeo(:, :, m, ntype) = 0.0D0
               end if
               !write(*,*) 'm=', m
               !write(*,'(18f10.6)') real(this%eeo(:,:,m,ntype))
               !write(*,*) 'ee', m
               !write(*,'(18f10.6)') real(this%ee(:,:,m,ntype))
            end do
         end if
      end do ! end of atom type number
      if (this%local_axis) then
         this%ee_glob = this%ee
         if (this%hoh) this%eeo_glob = this%eeo
      end if
   end subroutine build_bulkham

   subroutine build_locham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: ino, nr, nlim, m, i, j, ja, ji

      do nlim = 1, this%charge%lattice%nmax
         nr = this%charge%lattice%nn(nlim, 1) ! Number of neighbours considered
         ino = this%charge%lattice%num(nlim)
         call this%chbar_nc(nlim, nr, ino, nlim)
         do m = 1, nr
            do i = 1, 9
               do j = 1, 9
                  this%hall(j, i, m, nlim) = this%hmag(j, i, m, 4) + this%hmag(j, i, m, 3) ! H0+Hz
                  this%hall(j + 9, i + 9, m, nlim) = this%hmag(j, i, m, 4) - this%hmag(j, i, m, 3) ! H0-Hz
                  this%hall(j, i + 9, m, nlim) = this%hmag(j, i, m, 1) - i_unit * this%hmag(j, i, m, 2) ! Hx-iHy
                  this%hall(j + 9, i, m, nlim) = this%hmag(j, i, m, 1) + i_unit * this%hmag(j, i, m, 2) ! Hx+iHy
               end do
            end do
         end do
         if (this%hoh) then
            call this%build_obarm()
            call this%build_enim()
            do m = 1, nr
               ji = 0
               if (m > 1) then
                  ja = this%charge%lattice%nn(nlim, m)
                  if (ja .ne. 0) then
                     ji = this%charge%lattice%iz(ja)
                  end if
               else
                  ji = this%charge%lattice%iz(nlim)
               end if
               ! Check if neighbour 'm' exists for atom 'nlim', otherwise fill HoH Hamiltonian with zeros.
               if (ji > 0) then
                  call zgemm('n', 'n', 18, 18, 18, cone, this%hall(1, 1, m, nlim), 18, this%obarm(1, 1, ji), 18, czero, this%hallo(1, 1, m, nlim), 18)
               else
                  this%hallo(:, :, m, nlim) = 0.0D0
               end if
            end do
         end if
      end do
      if (this%local_axis) then
         this%hall_glob = this%hall
         if (this%hoh) this%hallo_glob = this%hallo
      end if

   end subroutine build_locham

   subroutine rs2pao(this)
      implicit none
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      real(rp), dimension(3) :: rij, rijtest
      integer :: i, j, k, idxi, idxj, idxk, ino, nr, ia, ipao, jpao
      integer :: jj, max_orbital, n_atoms
      integer :: ntype, iostatus
      real(rp), dimension(3) :: idx

      n_atoms = this%charge%lattice%ntype
      max_orbital = 9

      open (unit=92, file='rs2paoham.dat', action='write', iostat=iostatus, status='replace')
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype) ! Atom number in clust
         ino = this%charge%lattice%num(ia) ! Atom bravais type of ia
         nr = this%charge%lattice%nn(ia, 1) ! Number of neighbours considered
         do k = 1, nr
            jj = this%charge%lattice%nn(ia, k)
            !write(123, *)'ia, ii', ia, m, this%charge%lattice%nn(ia, m)
            if (k == 1) then
               jj = ia
            end if
            if (jj /= 0) then
               rij(:) = this%charge%lattice%cr(:, ia) - this%charge%lattice%cr(:, jj)

               rijtest(:) = 0.0D0
               do idxi = -5, 5
                  do idxj = -5, 5
                     do idxk = -5, 5
                        rijtest(:) = this%charge%lattice%cr(:, ia) - (this%charge%lattice%cr(:, this%charge%lattice%iz(jj)) &
                                                                      + idxi * this%charge%lattice%a(:, 1) + idxj * this%charge%lattice%a(:, 2) + idxk * this%charge%lattice%a(:, 3))
                        if (norm2(rij(:) - rijtest(:)) < 1.0D-3) then
                           idx(:) = [idxi, idxj, idxk]
                        end if
                     end do
                  end do
               end do
               !a_inv(:,:) = inverse_3x3(this%charge%lattice%a)
               !idx(:) = matmul(a_inv(:,:),b(:))

               !write(*,*) b(:), 'b'
               !write(*,*) matmul(this%charge%lattice%a(:,:),idx(:))

               !call hcpx(this%ee(1:9,1:9,k,ntype), 'sph2cart')
               !call hcpx(this%ee(1:9,10:18,k,ntype), 'sph2cart')
               !call hcpx(this%ee(10:18,1:9,k,ntype), 'sph2cart')
               !call hcpx(this%ee(10:18,10:18,k,ntype), 'sph2cart')
               if (k == 1) this%ee(:, :, k, ntype) = this%ee(:, :, k, ntype) + this%lsham(:, :, ntype)
               do i = 1, 18
                  do j = 1, 18
                     ipao = 0; jpao = 0
                     call site2orb(i, ia, ipao, n_atoms, max_orbital)
                     call site2orb(j, this%charge%lattice%iz(jj), jpao, n_atoms, max_orbital)
                     write (92, '(3I4,2I7,2F22.14)') int(idx(:)), ipao, jpao, real(this%ee(i, j, k, ntype)) * ry2ev, aimag(this%ee(i, j, k, ntype)) * ry2ev
                  end do
               end do
            end if
         end do
      end do
      close (92)
   end subroutine rs2pao

   subroutine build_from_paoflow_opt(this)
      implicit none
      type hamData
         integer :: idxi, idxj, idxk
         integer :: orbl, orbm
         real :: dumre, dumcmplx
      end type hamData
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, ino, nr, ia, iia, jja
      integer :: jj
      integer :: ntype, iostat1, iostatus, n_atoms, max_orbital, numLines
      real(rp), dimension(3) :: vet, vetpao

      integer, dimension(maxval(this%charge%lattice%nn(:, 1)) + 1, 3) :: idx
      type(hamData) :: ham
      type(hamData), allocatable :: hamArray(:)

      n_atoms = this%charge%lattice%ntype
      max_orbital = 9
      numLines = countLines('paoham.dat')
      allocate (hamArray(numLines))

      ! Reading Hamiltonian data once
      open (unit=92, file='paoham.dat', action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file containing the paoflow Hamiltonian not found', __FILE__, __LINE__)
      end if
      do i = 1, numLines
         read (92, *, iostat=iostat1) hamArray(i)%idxi, hamArray(i)%idxj, hamArray(i)%idxk, &
            hamArray(i)%orbl, hamArray(i)%orbm, hamArray(i)%dumre, hamArray(i)%dumcmplx
      end do
      close (92)
      write (*, *) 'PAOFLOW Hamiltonian has been read'
      call g_timer%start('Hamiltonian allocation')
      !$omp parallel do private(ntype, ia, ino, nr, jj, k, l, ham, vet, vetpao, i, j, iia, jja, idx) shared(this, numLines, hamArray, n_atoms, max_orbital)
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype)
         ino = this%charge%lattice%num(ia)
         nr = this%charge%lattice%nn(ia, 1)
         write (*, *) ia, nr, n_atoms, max_orbital
         do k = 1, nr
            jj = this%charge%lattice%nn(ia, k)
            if (k == 1) then
               jj = ia
            end if
            if (jj /= 0) then
               vet(:) = (this%charge%lattice%cr(:, ia) - this%charge%lattice%cr(:, jj))!*this%charge%lattice%alat
               do l = 1, numLines
                  ham = hamArray(l)
                  call orb2site(ham%orbl, i, iia, n_atoms, max_orbital)
                  call orb2site(ham%orbm, j, jja, n_atoms, max_orbital)
                  if (iia == ia) then
                     vetpao(:) = this%charge%lattice%cr(:, iia) - (this%charge%lattice%cr(:, jja) &
                                                                   + ham%idxi * this%charge%lattice%a(:, 1) + ham%idxj * this%charge%lattice%a(:, 2) + ham%idxk * this%charge%lattice%a(:, 3))!*this%charge%lattice%alat
                     if (norm2(vet(:) - vetpao(:)) < 1.0D-3) then
                        idx(k, :) = [ham%idxi, ham%idxj, ham%idxk]
                        this%ee(i, j, k, ntype) = cmplx(ham%dumre, ham%dumcmplx) / 13.605703976
                     end if
                  end if
               end do
            end if
            !write(128,*)'m=',k, 'Atom=', jj, 'Coordinates=', this%charge%lattice%cr(:, jj), 'Ntype=',ntype, 'Index=', idx(k,:)
            !write(129,*)'m=',k, 'Atom=', jj, 'Coordinates=', this%charge%lattice%cr(:, jj), 'Ntype=',ntype, 'Index=', idx(k,:)
            !write(128,'(18f10.6)') real(this%EE(1:18,1:18,k,ntype))*13.605703976
            !write(129,'(18f10.6)') aimag(this%EE(1:18,1:18,k,ntype))*13.605703976
            !write(128,*) sum(real(this%ee(:,:,k,ntype)))
            !write(129,*) sum(real(this%ee(:,:,k,ntype)))
         end do
      end do
      !$omp end parallel do
      !write(*,*) 'before rotation', k, ntype
      !write(*,'(18f10.6)') abs(this%ee(:,:,1,1))
      !call rotmag(this%ee,63*62)
      !write(*,*) 'before rotation'
      !write(*,'(18f10.6)') abs(this%ee(:,:,1,1))
      close (92)
      deallocate (hamArray)
      call g_timer%stop('Hamiltonian allocation')
   end subroutine build_from_paoflow_opt

   subroutine build_from_paoflow(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, ino, nr, ia, iia, jja
      integer :: jj, orbl, orbm, idxi, idxj, idxk
      integer :: ntype, iostat1, iostatus, n_atoms, max_orbital
      real(rp), dimension(3) :: vet, vetpao
      integer, dimension(maxval(this%charge%lattice%nn(:, 1)) + 1, 3) :: idx
      real(rp) :: dumre, dumcmplx

      n_atoms = this%charge%lattice%ntype
      max_orbital = 9
      open (unit=90, file='paoup.dat', action='read', iostat=iostatus, status='old')
      open (unit=91, file='paodw.dat', action='read', iostat=iostatus, status='old')
      open (unit=92, file='paoham.dat', action='read', iostat=iostatus, status='old')

      if (iostatus /= 0) then
         call g_logger%fatal('file containing the paoflow Hamiltonian not found', __FILE__, __LINE__)
      end if
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype) ! Atom number in clust
         ino = this%charge%lattice%num(ia) ! Atom bravais type of ia
         nr = this%charge%lattice%nn(ia, 1) ! Number of neighbours considered
         write (*, *) ia, nr, n_atoms, max_orbital
         do k = 1, nr
            jj = this%charge%lattice%nn(ia, k)
            !write(123, *)'ia, ii', ia, m, this%charge%lattice%nn(ia, m)
            if (k == 1) then
               jj = ia
            end if
            if (jj /= 0) then
               !cri_cart(:) = this%charge%lattice%cr(:, ia)
               !crj_cart(:) = this%charge%lattice%cr(:, jj)

               !cri_dir(:) = cartesian_to_direct(this%charge%lattice%a,cri_cart)
               !crj_dir(:) = cartesian_to_direct(this%charge%lattice%a,crj_cart)

               !vet(:) = cri_dir(:) - crj_dir(:)
               vet(:) = (this%charge%lattice%cr(:, ia) - this%charge%lattice%cr(:, jj))!*this%charge%lattice%alat
               do
                  read (92, *, iostat=iostat1) idxi, idxj, idxk, orbl, orbm, dumre, dumcmplx
                  if (iostat1 /= 0) then
                     exit
                  end if
                  if (orbl <= n_atoms * max_orbital) then
                     i = modulo(orbl - 1, max_orbital) + 1
                     iia = int((orbl - 1) / max_orbital) + 1
                  else
                     i = modulo(orbl - 1, max_orbital) + 10
                     iia = int((orbl - 1 - n_atoms * max_orbital) / max_orbital) + 1
                  end if

                  if (orbm <= n_atoms * max_orbital) then
                     j = modulo(orbm - 1, max_orbital) + 1
                     jja = int((orbm - 1) / max_orbital) + 1
                  else
                     j = modulo(orbm - 1, max_orbital) + 10
                     jja = int((orbm - 1 - n_atoms * max_orbital) / max_orbital) + 1
                  end if

                  !cri_cart(:) = this%charge%lattice%cr(:, iia)
                  !crj_cart(:) = this%charge%lattice%cr(:, jja)

                  !cri_dir(:) = cartesian_to_direct(this%charge%lattice%a,cri_cart)
                  !crj_dir(:) = cartesian_to_direct(this%charge%lattice%a,crj_cart)
                  if (iia == ia) then
                     !vetpao(:) = cri_dir(:) - (crj_dir(:) + [idxi,idxj,idxk])

                     vetpao(:) = this%charge%lattice%cr(:, iia) - (this%charge%lattice%cr(:, jja) &
                                                                   + idxi * this%charge%lattice%a(:, 1) + idxj * this%charge%lattice%a(:, 2) + idxk * this%charge%lattice%a(:, 3))!*this%charge%lattice%alat
                     if (norm2(vet(:) - vetpao(:)) < 1.0D-3) then
                        idx(k, :) = [idxi, idxj, idxk]
                        this%ee(i, j, k, ntype) = cmplx(dumre, dumcmplx) / 13.605703976
                     end if
                  end if
               end do
!          do
!            read(91,*,iostat=iostat2) idxi, idxj, idxk, orbl, orbm, dumre, dumcmplx
!            if (iostat2 /= 0) then
!              exit
!            end if
!            iia = orb_to_site(orbl,9)
!            jja = orb_to_site(orbm,9)
!
!            cri_cart(:) = this%charge%lattice%cr(:, iia)
!            crj_cart(:) = this%charge%lattice%cr(:, jja)
!
!            cri_dir(:) = cartesian_to_direct(this%charge%lattice%a,cri_cart)
!            crj_dir(:) = cartesian_to_direct(this%charge%lattice%a,crj_cart)
!            if(iia==ia)then
!              vetpao(:) = cri_dir(:) - (crj_dir(:) + [idxi,idxj,idxk])
!
!!              vetpao(:) = this%charge%lattice%cr(:,iia) - this%charge%lattice%cr(:,jja) &
!!                         + (idxi*this%charge%lattice%a(:,1) + idxj*this%charge%lattice%a(:,2) + idxk*this%charge%lattice%a(:,3))!*this%charge%lattice%alat
!              if(norm2(vet(:)-vetpao(:))<1.0d-3)then
!                idxdw(k,:) = [idxi,idxj,idxk]
!                this%ee(orbl+9-(iia-1)*9, orbm+9-(jja-1)*9,k,ntype) = cmplx(dumre,dumcmplx)/13.605703976
!              end if
!            end if
!          end do
            end if
            write (128, *) 'm=', k, 'Atom=', jj, 'Coordinates=', this%charge%lattice%cr(:, jj), 'Ntype=', ntype, 'Index=', idx(k, :)
            write (128, '(18f10.6)') real(this%EE(1:18, 1:18, k, ntype))!*13.605703976
            write (128, *) sum(real(this%ee(:, :, k, ntype)))
            !write(129,*)'m=',k, 'Atom=', jj, 'Coordinates=', this%charge%lattice%cr(:, jj), 'Ntype=',ntype, 'Index=', idx(k,:)
            !write(129,'(9f10.6)') real(this%EE(10:18,1:9,k,ntype))*13.605703976
            !write(129,*) sum(real(this%ee(:,:,k,ntype)))
            rewind (90)
            rewind (91)
            rewind (92)
         end do
      end do
   end subroutine build_from_paoflow

   subroutine ham0m_nc(this, it, jt, vet, hhh)
      class(hamiltonian), intent(inout) :: this
      ! Input
      integer, intent(in) :: it, jt ! Type of atom i and j
      real(rp), dimension(3), intent(in) :: vet
      real(rp), dimension(9, 9), intent(in) :: hhh
      ! Local Variables
      integer :: i, j, ilm, jlm, m
      complex(rp), dimension(3) :: cross
      complex(rp), dimension(9, 9) :: hhhc
      complex(rp), dimension(this%charge%lattice%ntype, 3) :: momc
      complex(rp) :: dot
      real(rp) :: vv

      this%hhmag(:, :, :) = 0.0D0

      vv = norm2(vet)

      ! Real to complex
      dot = cmplx(dot_product(this%charge%lattice%symbolic_atoms(it)%potential%mom, this%charge%lattice%symbolic_atoms(jt)%potential%mom), kind=kind(0.0D0))
      do i = 1, this%charge%lattice%ntype
         do j = 1, 3
            momc(i, j) = cmplx(this%charge%lattice%symbolic_atoms(i)%potential%mom(j), kind=kind(0.0D0))
         end do
      end do
      cross = cmplx(cross_product(this%charge%lattice%symbolic_atoms(it)%potential%mom, this%charge%lattice%symbolic_atoms(jt)%potential%mom), kind=kind(0.0D0))
      hhhc(:, :) = cmplx(hhh(:, :), kind=kind(0.0D0))

      do ilm = 1, 9
         do jlm = 1, 9
            this%hhmag(ilm, jlm, 4) = &
               this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm) * hhhc(ilm, jlm) * this%charge%lattice%symbolic_atoms(jt)%potential%wx0(jlm) + &
               this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm) * hhhc(ilm, jlm) * this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm) * dot
         end do
      end do

!    do ilm=1, 9
!      write(123, '(9f10.6)') (real(this%hhmag(ilm, jlm, 4)), jlm=1, 9)
!    end do

      if (vv <= 0.01D0) then
         do ilm = 1, 9
            if (this%hoh) then
               this%hhmag(ilm, ilm, 4) = this%hhmag(ilm, ilm, 4) + this%charge%lattice%symbolic_atoms(it)%potential%cex0(ilm)
            else
               this%hhmag(ilm, ilm, 4) = this%hhmag(ilm, ilm, 4) + this%charge%lattice%symbolic_atoms(it)%potential%cx0(ilm)
            end if
         end do
      end if

      do m = 1, 3
         do jlm = 1, 9
            do ilm = 1, 9
               this%hhmag(ilm, jlm, m) = &
                  (this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm) * hhhc(ilm, jlm) * this%charge%lattice%symbolic_atoms(jt)%potential%wx0(jlm)) * momc(it, m) + &
                  (this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm) * hhhc(ilm, jlm) * this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm)) * momc(jt, m) + &
                  i_unit * this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm) * hhhc(ilm, jlm) * this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm) * cross(m)
            end do
         end do
      end do

      if (vv > 0.01D0) return
      do m = 1, 3
         do ilm = 1, 9
            if (this%hoh) then
               this%hhmag(ilm, ilm, m) = this%hhmag(ilm, ilm, m) + this%charge%lattice%symbolic_atoms(it)%potential%cex1(ilm) * momc(it, m)
            else
               this%hhmag(ilm, ilm, m) = this%hhmag(ilm, ilm, m) + this%charge%lattice%symbolic_atoms(it)%potential%cx1(ilm) * momc(it, m)
            end if
         end do
      end do

      !do m=1, 3
      !  write(123, *)'m=', m
      !  do ilm=1, 9
      !    write(123, '(9f10.6)') (real(this%hhmag(ilm, jlm, m)), jlm=1, 9)
      !  end do
      !end do
   end subroutine ham0m_nc

   subroutine chbar_nc(this, ia, nr, ino, ntype)
      class(hamiltonian), intent(inout) :: this
      ! Input
      integer, intent(in) :: ia ! Atom number in clust
      integer, intent(in) :: nr ! Number of neighbours considered
      integer, intent(in) :: ino ! Atom bravais type of ia
      integer, intent(in) :: ntype ! Atom type
      ! Local variables
      real(rp) :: r2
      real(rp), dimension(3, size(this%charge%lattice%cr(1, :))) :: cralat ! Clust position times the lattice constant
      real(rp), dimension(3) :: vet
      real(rp), dimension(9, 9) :: hhh
      integer :: m, it, jt, jj, dummy
      integer :: ni, mdir
      integer :: kk ! Clust size number

      this%hmag(:, :, :, :) = 0.0D0

      r2 = this%charge%lattice%r2
      cralat(:, :) = this%charge%lattice%cr(:, :) * this%charge%lattice%alat
      kk = this%charge%lattice%kk

      call this%charge%lattice%clusba(r2, cralat, ia, kk, kk, dummy)

      !do m=1, nr
      !  print '(9f10.6)', real(this%charge%lattice%sbar(:, :, m, ino))
      !end do
      it = this%charge%lattice%iz(ia)
      do m = 1, nr
         jj = this%charge%lattice%nn(ia, m)
         !write(123, *)'ia, ii', ia, m, this%charge%lattice%nn(ia, m)
         if (m == 1) then
            jj = ia
         end if
         if (jj /= 0) then
            jt = this%charge%lattice%iz(jj)
            vet(:) = (this%charge%lattice%cr(:, jj) - this%charge%lattice%cr(:, ia)) * this%charge%lattice%alat
            !write(123, '(3f10.6)') vet(:)
            !write(123, '(3f10.6)') this%charge%lattice%sbarvec(:, m)
            !write(123, '(a, 3i4, 3f10.6)') 'nn ', IA, m, JJ, VET(:)
            call this%hmfind(vet, nr, hhh, m, ia, m, ni, ntype)
            if (ni == 0) then
               this%charge%lattice%nn(ia, m) = 0
            end if
            call this%ham0m_nc(it, jt, vet, hhh)
            do mdir = 1, 4
               call hcpx(this%hhmag(:, :, mdir), 'cart2sph')
               this%hmag(:, :, m, mdir) = this%hhmag(:, :, mdir)
            end do
         end if
      end do
      !do m=1, nr
      !  write(123, *)'m=', m
      !  do mdir=1, 4
      !    write(123, *)'mdir=', mdir
      !    do i=1, 9
      !      write(123, '(9f10.4)')(real(this%hmag(i, j, m, mdir)), j=1, 9)
      !    end do
      !  end do
      !end do
   end subroutine chbar_nc

   subroutine hmfind(this, vet, nr, hhh, m, ia, jn, ni, ntype)
      class(hamiltonian), intent(inout) :: this
      ! Input
      integer, intent(in) :: ntype ! Atom type
      integer, intent(in) :: m ! Number of the given neighbour
      integer, intent(in) :: ia ! Atom number in clust
      integer, intent(in) :: jn ! ?
      integer, intent(in) :: nr ! Number of neighbours
      real(rp), dimension(3), intent(in) :: vet
      ! Output
      integer, intent(out) :: ni
      real(rp), dimension(9, 9), intent(inout) :: hhh
      ! Local variables
      real(rp) :: a1, a2, a3, aaa, eps
      integer :: i, ilm, jlm

      eps = 0.0001D0
      ni = 1
      a1 = 0.0D0
      a2 = 0.0D0
      a3 = 0.0D0
      aaa = 0.0D0
      do i = 1, nr
         !write(123, '(a, i4, 3f10.4)')'i', i, this%charge%lattice%sbarvec(:, i)
         a1 = (vet(1) - this%charge%lattice%sbarvec(1, i))
         a2 = (vet(2) - this%charge%lattice%sbarvec(2, i))
         a3 = (vet(3) - this%charge%lattice%sbarvec(3, i))
         aaa = a1**2 + a2**2 + a3**2
         if (aaa < eps) goto 1000
      end do
      write (*, '(1x, a, i4, a, i4, a, 3f10.6)') ' Error in hamiltonian%hmfind: Neighbour vector not found for atom', ia, &
         ' neighbour', jn, 'vector', vet(:)

      ni = 0
1000  continue
      do ilm = 1, 9
         do jlm = 1, 9
            hhh(ilm, jlm) = real(this%charge%lattice%sbar(jlm, ilm, m, this%charge%lattice%num(ia)))
         end do
      end do

      !do ilm = 1, 9
      !    write(123, '(9f10.6)')(hhh(ilm, jlm), jlm=1, 9)
      !end do
   end subroutine hmfind

   subroutine orb2site(orb, i_out, ia_out, n_atoms, max_orbital)
      integer, intent(in) :: orb, n_atoms, max_orbital
      integer, intent(out) :: i_out, ia_out

      if (orb <= n_atoms * max_orbital) then
         i_out = modulo(orb - 1, max_orbital) + 1
         ia_out = int((orb - 1) / max_orbital) + 1
      else
         i_out = modulo(orb - 1, max_orbital) + 10
         ia_out = int((orb - 1 - n_atoms * max_orbital) / max_orbital) + 1
      end if
   end subroutine orb2site

   subroutine site2orb(i_in, ia_in, orb_out, n_atoms, max_orbital)
      integer, intent(in) :: i_in, ia_in, n_atoms, max_orbital
      integer, intent(out) :: orb_out

      if (i_in <= max_orbital) then
         orb_out = (ia_in - 1) * max_orbital + i_in
      else
         orb_out = (ia_in - 1) * max_orbital + i_in + (n_atoms - 1) * max_orbital
      end if
   end subroutine site2orb

   ! Rotate Hamiltonian to local axis
   subroutine rotate_to_local_axis(this, m_loc)
      use math_mod, only: rotmag_loc
      class(hamiltonian), intent(inout) :: this
      real(rp), dimension(3), intent(in) :: m_loc

      ! Local variables
      integer :: sdim
      ! Rotate Hamiltonian to local axis if wanted
      if (this%local_axis) then
         sdim = product(shape(this%hall)) / 18 / 18
         call rotmag_loc(this%hall, this%hall_glob, sdim, m_loc)
         sdim = product(shape(this%ee)) / 18 / 18
         call rotmag_loc(this%ee, this%ee_glob, sdim, m_loc)
         if (this%hoh) then
            sdim = product(shape(this%eeo)) / 18 / 18
            call rotmag_loc(this%eeo, this%eeo_glob, sdim, m_loc)
            sdim = product(shape(this%hallo)) / 18 / 18
            call rotmag_loc(this%hallo, this%hallo_glob, sdim, m_loc)
            sdim = product(shape(this%enim)) / 18 / 18
            call rotmag_loc(this%enim, this%enim_glob, sdim, m_loc)
         end if
      end if
   end subroutine rotate_to_local_axis

   ! Rotate Hamiltonian back from local axis
   subroutine rotate_from_local_axis(this, m_loc)
      use math_mod, only: rotmag_loc
      class(hamiltonian), intent(inout) :: this
      real(rp), dimension(3), intent(in) :: m_loc

      ! Local variables

      ! Rotate Hamiltonian to local axis if wanted
      if (this%local_axis) then
         this%hall = this%hall_glob
         this%ee = this%ee_glob
         if (this%hoh) then
            this%eeo = this%eeo_glob
            this%hallo = this%hallo_glob
            this%enim = this%enim_glob
         end if
      end if
   end subroutine rotate_from_local_axis
end module hamiltonian_mod
