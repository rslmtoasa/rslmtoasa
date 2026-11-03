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
   public :: select_spin_operator, select_orbital_operator, select_quadrupole_operator, select_octupole_operator

   !> Module´s main procedure
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
      complex(rp), dimension(:, :, :, :), allocatable :: ee, eeo, eeoee
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
      !> Kubo polarization components and velocity directions
      character(len=10) :: pol_alpha, pol_beta
      real(rp), dimension(3) :: v_alpha, v_beta
      real(rp), dimension(:), allocatable :: velocity_scale
      !> Sparse Real Space Hamiltonian
      complex(rp), dimension(:, :), allocatable :: h_sparse
   contains
      procedure :: build_lsham
      procedure :: build_bulkham
      procedure :: build_locham
      procedure :: build_obarm
      procedure :: build_enim
      procedure :: build_from_paoflow
      procedure :: build_from_paoflow_opt
      procedure :: block_to_sparse
      procedure :: torque_operator_collinear
      procedure :: rs2pao
      procedure :: chbar_nc
      procedure :: ham0m_nc
      procedure :: hmfind
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: rotate_to_local_axis
      procedure :: rotate_from_local_axis
      final     :: destructor
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
      if (allocated(this%eeoee)) call g_safe_alloc%deallocate('hamiltonian.eeoee', this%eeoee)
      if (allocated(this%hallo)) call g_safe_alloc%deallocate('hamiltonian.hallo', this%hallo)
      if (allocated(this%obarm)) call g_safe_alloc%deallocate('hamiltonian.obarm', this%obarm)
      if (allocated(this%enim)) call g_safe_alloc%deallocate('hamiltonian.enim', this%enim)
      if (allocated(this%ee_glob)) call g_safe_alloc%deallocate('hamiltonian.ee_glob', this%ee_glob)
      if (allocated(this%eeo_glob)) call g_safe_alloc%deallocate('hamiltonian.eeo_glob', this%eeo_glob)
      if (allocated(this%enim_glob)) call g_safe_alloc%deallocate('hamiltonian.enim_glob', this%enim_glob)
      if (allocated(this%h_sparse)) call g_safe_alloc%deallocate('hamiltonian.h_sparse', this%h_sparse)
      if (allocated(this%velocity_scale)) call g_safe_alloc%deallocate('hamiltonian.velocity_scale', this%velocity_scale)
#else
      if (allocated(this%lsham)) deallocate (this%lsham)
      if (allocated(this%tmat)) deallocate (this%tmat)
      if (allocated(this%ee)) deallocate (this%ee)
      if (allocated(this%eeo)) deallocate (this%eeo)
      if (allocated(this%eeoee)) deallocate (this%eeoee)
      if (allocated(this%hmag)) deallocate (this%hmag)
      if (allocated(this%hhmag)) deallocate (this%hhmag)
      if (allocated(this%hall)) deallocate (this%hall)
      if (allocated(this%hallo)) deallocate (this%hallo)
      if (allocated(this%obarm)) deallocate (this%obarm)
      if (allocated(this%enim)) deallocate (this%enim)
      if (allocated(this%ee_glob)) deallocate (this%ee_glob)
      if (allocated(this%eeo_glob)) deallocate (this%eeo_glob)
      if (allocated(this%enim_glob)) deallocate (this%enim_glob)
      if (allocated(this%h_sparse)) deallocate(this%h_sparse)
      if (allocated(this%velocity_scale)) deallocate(this%velocity_scale)
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
      integer :: iostatus, funit, i

      include 'include_codes/namelists/hamiltonian.f90'

      hoh = this%hoh
      local_axis = this%local_axis
      orb_pol = this%orb_pol
      v_alpha(:) = this%v_alpha(:)
      v_beta(:) = this%v_beta(:)
      pol_alpha = this%pol_alpha
      pol_beta  = this%pol_beta
      call move_alloc(this%velocity_scale, velocity_scale)

      ! Reading
      open (newunit=funit, file=this%control%fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//trim(this%control%fname)//' not found', __FILE__, __LINE__)
      end if

      read (funit, nml=hamiltonian, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%fatal('Invalid &hamiltonian namelist; iostatus = '//int2str(iostatus), __FILE__, __LINE__)
      end if
      close (funit)

      this%hoh = hoh
      this%local_axis = local_axis
      this%orb_pol = orb_pol
      this%v_alpha(:) = v_alpha(:)
      this%v_beta(:) = v_beta(:)
      this%pol_alpha = pol_alpha
      this%pol_beta  = pol_beta
      call move_alloc(velocity_scale, this%velocity_scale)

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
      !if (hoh) then
      call g_safe_alloc%allocate('hamiltonian.eeo', this%eeo, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.eeoee', this%eeoee, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hallo', this%hallo, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.obarm', this%obarm, (/18, 18, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.enim', this%enim, (/18, 18, this%charge%lattice%ntype/))
      !end if
      !if (local_axis)  then
      call g_safe_alloc%allocate('hamiltonian.hall_glob', this%hall_glob, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.ee_glob', this%ee_glob, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      !if (hoh) then
      call g_safe_alloc%allocate('hamiltonian.ee0_glob', this%eeo_glob, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hallo_glob', this%hallo_glob, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.enim_glob', this%enim_glob, (/18, 18, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.velocity_scale', this%velocity_scale, (/this%charge%lattice%ntype/))
      !end if
      !end if
#else
      allocate (this%lsham(18, 18, this%charge%lattice%ntype))
      allocate (this%tmat(18, 18, 3, this%charge%lattice%ntype))
      allocate (this%hhmag(9, 9, 4), this%hmag(9, 9, this%charge%lattice%kk, 4))
      allocate (this%ee(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%hall(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      !if (this%hoh) then
      allocate (this%eeo(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%eeoee(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%hallo(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      allocate (this%obarm(18, 18, this%charge%lattice%ntype))
      allocate (this%enim(18, 18, this%charge%lattice%ntype))
      !end if
      !if (this%local_axis) then
      allocate (this%ee_glob(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%hall_glob(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      !if (this%hoh) then
      allocate (this%eeo_glob(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%hallo_glob(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      allocate (this%enim_glob(18, 18, this%charge%lattice%ntype))
      ! Velocity operators
      allocate (this%velocity_scale(this%charge%lattice%ntype))
      !end if
      !end if
#endif

      this%lsham(:, :, :) = 0.0d0
      this%tmat(:, :, :, :) = 0.0d0
      this%hhmag(:, :, :) = 0.0d0
      this%hmag(:, :, :, :) = 0.0d0
      this%hall(:, :, :, :) = 0.0d0
      this%ee(:, :, :, :) = 0.0d0
      !if (this%hoh) then
      this%hallo(:, :, :, :) = 0.0d0
      this%eeo(:, :, :, :) = 0.0d0
      this%eeoee(:, :, :, :) = 0.0d0
      this%obarm(:, :, :) = 0.0d0
      this%enim(:, :, :) = 0.0d0
      !end if
      !if (this%local_axis) then
      this%hall_glob(:, :, :, :) = 0.0d0
      this%ee_glob(:, :, :, :) = 0.0d0
      !  if (this%hoh) then
      this%hallo_glob(:, :, :, :) = 0.0d0
      this%eeo_glob(:, :, :, :) = 0.0d0
      this%enim_glob(:, :, :) = 0.0d0
      !  end if
      !end if
      this%velocity_scale(:) = 1.0d0
      this%hoh = .false.
      this%local_axis = .false.
      this%orb_pol = .false.
      this%v_alpha(:) = [1, 0, 0]
      this%v_beta(:) = [1, 0, 0] 
      this%pol_alpha = 'z'
      this%pol_beta  = 'z'
   end subroutine restore_to_default

   subroutine block_to_sparse(this)
       !*************************************************************************
       !> @brief Constructs the sparse Hamiltonian matrix for a real-space cluster.
       !>
       !> This subroutine loops through all cluster atoms and their neighbors to
       !> assemble the sparse Hamiltonian matrix in a block-wise manner. The matrix
       !> size is determined by the number of atoms in the cluster and their spd basis.
       !>
       !> @param[in] this       Hamiltonian type derived
       !> @param[out] H_sparse  Sparse matrix structure to store the Hamiltonian.
       !*************************************************************************
       class(hamiltonian), intent(inout) :: this
   
       ! Local variables
       integer :: kk, m, nr, i, j, i_start, j_start
       integer :: neighbor  ! Neighbor atom index
       real(rp), dimension(3) :: rij  ! Displacement vector (unused for now but available)
  
       if (allocated(this%h_sparse)) deallocate(this%h_sparse) 
       allocate(this%h_sparse(this%lattice%kk*18,this%lattice%kk*18))

       ! Loop over cluster atoms
       do kk = 1, this%lattice%kk  ! Loop over all atoms in the cluster
           ! Number of neighbors for the kk-th atom
           nr = this%charge%lattice%nn(kk, 1)
           ! Loop over neighbors (including onsite, m = 1)
           do m = 1, nr
               ! Compute block indices in the global matrix
               i_start = 18 * (kk - 1) + 1
               if (m == 1) then
                   j_start = i_start  ! Onsite term: diagonal block
               else
                   neighbor = this%charge%lattice%nn(kk, m)  ! Neighbor atom index
                   j_start = 18 * (neighbor - 1) + 1
               end if
               ! Add the block to the sparse matrix
               if (neighbor .ne. 0) then
                  this%h_sparse(i_start:i_start+17, j_start:j_start+17) = 0.0d0
                  do i = 1, 18
                      do j = 1, 18
                         this%h_sparse(i_start + i - 1, j_start + j - 1) = this%h_sparse(i_start + i - 1, j_start + j - 1) &
                                                                           + this%ee(i, j, m, 1)
                      end do
                  end do
               end if
           end do
       end do
       ! Placeholder: Incorporate atom type and neighbor type logic in future if required.
   end subroutine block_to_sparse


   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Select spin operator S_x, S_y, or S_z from a polarization label.
   !---------------------------------------------------------------------------
   subroutine select_spin_operator(pol, S_op)
      character(len=*), intent(in) :: pol
      complex(rp), dimension(18, 18), intent(out) :: S_op

      S_op(:, :) = (0.0_rp, 0.0_rp)

      select case(trim(pol))
      case('x')
         S_op(:, :) = S_x(:, :)
      case('y')
         S_op(:, :) = S_y(:, :)
      case('z')
         S_op(:, :) = S_z(:, :)
      case default
         call g_logger%fatal('Invalid spin polarization: '//trim(pol), __FILE__, __LINE__)
      end select

   end subroutine select_spin_operator


   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Select orbital operator L_x, L_y, or L_z from a polarization label.
   !> The math_mod matrices are transformed from Cartesian to spherical basis.
   !---------------------------------------------------------------------------
   subroutine select_orbital_operator(pol, L_op)
      character(len=*), intent(in) :: pol
      complex(rp), dimension(18, 18), intent(out) :: L_op

      complex(rp), dimension(9, 9) :: mLx, mLy, mLz

      mLx(:, :) = L_x(:, :)
      mLy(:, :) = L_y(:, :)
      mLz(:, :) = L_z(:, :)

      call hcpx(mLx, 'cart2sph')
      call hcpx(mLy, 'cart2sph')
      call hcpx(mLz, 'cart2sph')

      L_op(:, :) = (0.0_rp, 0.0_rp)

      select case(trim(pol))
      case('x')
         L_op(1:9, 1:9)     = mLx(:, :)
         L_op(10:18,10:18)  = mLx(:, :)
      case('y')
         L_op(1:9, 1:9)     = mLy(:, :)
         L_op(10:18,10:18)  = mLy(:, :)
      case('z')
         L_op(1:9, 1:9)     = mLz(:, :)
         L_op(10:18,10:18)  = mLz(:, :)
      case default
         call g_logger%fatal('Invalid orbital polarization: '//trim(pol), __FILE__, __LINE__)
      end select

   end subroutine select_orbital_operator

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Select orbital quadrupole / OAP operator from a component label.
   !>
   !> Definitions:
   !>   Qxx = Lx^2 - L^2/3
   !>   Qyy = Ly^2 - L^2/3
   !>   Qzz = Lz^2 - L^2/3
   !>   Qxy = 1/2 {Lx,Ly}
   !>   Qyz = 1/2 {Ly,Lz}
   !>   Qzx = 1/2 {Lz,Lx}
   !>   Qx2y2  = (Qxx - Qyy)/2
   !>   Q3z2r2 = (2 Qzz - Qxx - Qyy)/2
   !>
   !> L^2 = Lx^2 + Ly^2 + Lz^2. Matrices use hbar = 1; the
   !> physical quadrupoles are angular-momentum bilinears divided by hbar.
   !> The diagonal Cartesian components now form a traceless tensor.
   !> Both accumulation and current builders use this selector.
   !>
   !> The math_mod L matrices are transformed from Cartesian to spherical basis
   !> before forming products, exactly as done for orbital accumulation/current.
   !---------------------------------------------------------------------------
   subroutine select_quadrupole_operator(pol, Q_op)
      character(len=*), intent(in) :: pol
      complex(rp), dimension(18, 18), intent(out) :: Q_op

      complex(rp), dimension(9, 9) :: mLx, mLy, mLz
      complex(rp), dimension(9, 9) :: Q9
      complex(rp), dimension(9, 9) :: Lx2, Ly2, Lz2, L2

      mLx(:, :) = L_x(:, :)
      mLy(:, :) = L_y(:, :)
      mLz(:, :) = L_z(:, :)

      call hcpx(mLx, 'cart2sph')
      call hcpx(mLy, 'cart2sph')
      call hcpx(mLz, 'cart2sph')

      Lx2 = matmul(mLx, mLx)
      Ly2 = matmul(mLy, mLy)
      Lz2 = matmul(mLz, mLz)
      L2 = Lx2 + Ly2 + Lz2

      Q9(:, :) = (0.0_rp, 0.0_rp)

      select case(trim(pol))

      case('Qxx', 'qxx', 'xx')
         Q9 = Lx2 - L2 / 3.0_rp

      case('Qyy', 'qyy', 'yy')
         Q9 = Ly2 - L2 / 3.0_rp

      case('Qzz', 'qzz', 'zz')
         Q9 = Lz2 - L2 / 3.0_rp

      case('Qxy', 'qxy', 'xy')
         Q9 = 0.5_rp * (matmul(mLx, mLy) + matmul(mLy, mLx))

      case('Qyz', 'qyz', 'yz')
         Q9 = 0.5_rp * (matmul(mLy, mLz) + matmul(mLz, mLy))

      case('Qzx', 'qzx', 'zx', 'Qxz', 'qxz', 'xz')
         Q9 = 0.5_rp * (matmul(mLz, mLx) + matmul(mLx, mLz))

      case('Qx2y2', 'qx2y2', 'x2y2')
         Q9 = 0.5_rp * (Lx2 - Ly2)

      case('Q3z2r2', 'q3z2r2', '3z2r2')
         Q9 = 0.5_rp * (2.0_rp * Lz2 - Lx2 - Ly2)

      case default
         call g_logger%fatal('Invalid quadrupole component: '//trim(pol), __FILE__, __LINE__)

      end select

      ! Clean tiny numerical anti-Hermitian noise.
      Q9 = 0.5_rp * (Q9 + transpose(conjg(Q9)))

      Q_op(:, :) = (0.0_rp, 0.0_rp)

      ! spin-up block
      Q_op(1:9, 1:9) = Q9(:, :)

      ! spin-down block
      Q_op(10:18, 10:18) = Q9(:, :)

      write(*,*) 'Selected Q component = ', trim(pol), ' norm = ', sqrt(real(sum(conjg(Q_op)*Q_op), rp))
      write(*,*) 'Hermiticity check Q = ', maxval(abs(Q_op - transpose(conjg(Q_op))))

   end subroutine select_quadrupole_operator

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Select a local orbital octupole, in hbar = 1 units.
   !> Oabb (a /= b) = {La,Lb^2}/2 - {La,L^2}/10 - La/10.
   !> These seven pure rank-three components include fully symmetrized Oxyz.
   !> Oa5b2r2 is an alias of Oabb, with the SAME normalization (no factor 5).
   !> Explicit Ma3b2r2 labels retain the old mixed dipole/octupole observable.
   !> Old Oa3b2r2 labels are rejected to prevent silent reinterpretation.
   !> Physical operators in angular-momentum units are these matrices times hbar.
   !---------------------------------------------------------------------------
   subroutine select_octupole_operator(pol, O_op)
      character(len=*), intent(in) :: pol
      complex(rp), dimension(18, 18), intent(out) :: O_op

      complex(rp), dimension(9, 9) :: mLx, mLy, mLz
      complex(rp), dimension(9, 9) :: Lx2, Ly2, Lz2, L2, O9

      mLx(:, :) = L_x(:, :)
      mLy(:, :) = L_y(:, :)
      mLz(:, :) = L_z(:, :)

      call hcpx(mLx, 'cart2sph')
      call hcpx(mLy, 'cart2sph')
      call hcpx(mLz, 'cart2sph')

      Lx2 = matmul(mLx, mLx)
      Ly2 = matmul(mLy, mLy)
      Lz2 = matmul(mLz, mLz)
      L2  = Lx2 + Ly2 + Lz2

      O9(:, :) = (0.0_rp, 0.0_rp)

      select case(trim(pol))

      case('Oxyz', 'oxyz', 'xyz')
         O9 = ( &
              matmul(matmul(mLx, mLy), mLz) &
            + matmul(matmul(mLx, mLz), mLy) &
            + matmul(matmul(mLy, mLx), mLz) &
            + matmul(matmul(mLy, mLz), mLx) &
            + matmul(matmul(mLz, mLx), mLy) &
            + matmul(matmul(mLz, mLy), mLx) ) / 6.0_rp

      case('Oxyy', 'oxyy', 'xyy', 'Ox5y2r2', 'ox5y2r2')
         O9 = 0.5_rp * (matmul(mLx, Ly2) + matmul(Ly2, mLx)) &
            - (matmul(mLx, L2) + matmul(L2, mLx)) / 10.0_rp - mLx / 10.0_rp

      case('Mx3y2r2', 'mx3y2r2')
         O9 = 0.5_rp * (matmul(mLx, Ly2) + matmul(Ly2, mLx)) &
            - (matmul(mLx, L2) + matmul(L2, mLx)) / 6.0_rp

      case('Oyzz', 'oyzz', 'yzz', 'Oy5z2r2', 'oy5z2r2')
         O9 = 0.5_rp * (matmul(mLy, Lz2) + matmul(Lz2, mLy)) &
            - (matmul(mLy, L2) + matmul(L2, mLy)) / 10.0_rp - mLy / 10.0_rp

      case('My3z2r2', 'my3z2r2')
         O9 = 0.5_rp * (matmul(mLy, Lz2) + matmul(Lz2, mLy)) &
            - (matmul(mLy, L2) + matmul(L2, mLy)) / 6.0_rp

      case('Ozxx', 'ozxx', 'zxx', 'Oz5x2r2', 'oz5x2r2')
         O9 = 0.5_rp * (matmul(mLz, Lx2) + matmul(Lx2, mLz)) &
            - (matmul(mLz, L2) + matmul(L2, mLz)) / 10.0_rp - mLz / 10.0_rp

      case('Mz3x2r2', 'mz3x2r2')
         O9 = 0.5_rp * (matmul(mLz, Lx2) + matmul(Lx2, mLz)) &
            - (matmul(mLz, L2) + matmul(L2, mLz)) / 6.0_rp

      case('Oxzz', 'oxzz', 'xzz', 'Ox5z2r2', 'ox5z2r2')
         O9 = 0.5_rp * (matmul(mLx, Lz2) + matmul(Lz2, mLx)) &
            - (matmul(mLx, L2) + matmul(L2, mLx)) / 10.0_rp - mLx / 10.0_rp

      case('Mx3z2r2', 'mx3z2r2')
         O9 = 0.5_rp * (matmul(mLx, Lz2) + matmul(Lz2, mLx)) &
            - (matmul(mLx, L2) + matmul(L2, mLx)) / 6.0_rp

      case('Oyxx', 'oyxx', 'yxx', 'Oy5x2r2', 'oy5x2r2')
         O9 = 0.5_rp * (matmul(mLy, Lx2) + matmul(Lx2, mLy)) &
            - (matmul(mLy, L2) + matmul(L2, mLy)) / 10.0_rp - mLy / 10.0_rp

      case('My3x2r2', 'my3x2r2')
         O9 = 0.5_rp * (matmul(mLy, Lx2) + matmul(Lx2, mLy)) &
            - (matmul(mLy, L2) + matmul(L2, mLy)) / 6.0_rp

      case('Ozyy', 'ozyy', 'zyy', 'Oz5y2r2', 'oz5y2r2')
         O9 = 0.5_rp * (matmul(mLz, Ly2) + matmul(Ly2, mLz)) &
            - (matmul(mLz, L2) + matmul(L2, mLz)) / 10.0_rp - mLz / 10.0_rp

      case('Mz3y2r2', 'mz3y2r2')
         O9 = 0.5_rp * (matmul(mLz, Ly2) + matmul(Ly2, mLz)) &
            - (matmul(mLz, L2) + matmul(L2, mLz)) / 6.0_rp

      case('Ox3y2r2', 'ox3y2r2')
         call g_logger%fatal('Ambiguous legacy component: use Oxyy for pure rank three or Mx3y2r2 for the old mixed operator', &
                             __FILE__, __LINE__)

      case('Oy3z2r2', 'oy3z2r2')
         call g_logger%fatal('Ambiguous legacy component: use Oyzz for pure rank three or My3z2r2 for the old mixed operator', &
                             __FILE__, __LINE__)

      case('Oz3x2r2', 'oz3x2r2')
         call g_logger%fatal('Ambiguous legacy component: use Ozxx for pure rank three or Mz3x2r2 for the old mixed operator', &
                             __FILE__, __LINE__)

      case('Ox3z2r2', 'ox3z2r2')
         call g_logger%fatal('Ambiguous legacy component: use Oxzz for pure rank three or Mx3z2r2 for the old mixed operator', &
                             __FILE__, __LINE__)

      case('Oy3x2r2', 'oy3x2r2')
         call g_logger%fatal('Ambiguous legacy component: use Oyxx for pure rank three or My3x2r2 for the old mixed operator', &
                             __FILE__, __LINE__)

      case('Oz3y2r2', 'oz3y2r2')
         call g_logger%fatal('Ambiguous legacy component: use Ozyy for pure rank three or Mz3y2r2 for the old mixed operator', &
                             __FILE__, __LINE__)

      case default
         call g_logger%fatal('Invalid octupole component: '//trim(pol), __FILE__, __LINE__)

      end select

      O9 = 0.5_rp * (O9 + transpose(conjg(O9)))

      O_op(:, :) = (0.0_rp, 0.0_rp)
      O_op(1:9, 1:9) = O9(:, :)
      O_op(10:18, 10:18) = O9(:, :)

   end subroutine select_octupole_operator

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Build the spin-orbit coupling Hamiltonian according to
   !> Wu/Freeman, Phys. Rev. B 54, 61 (1996).
   !>
   !> The SOC parameters xi_p and xi_d are combined using a signed
   !> geometric mean, allowing negative SOC values to be specified
   !> in the potential files.
   !---------------------------------------------------------------------------
   subroutine build_lsham(this)
   
      class(hamiltonian), intent(inout) :: this
   
      ! Local variables
      integer :: i, j, k
   
      complex(rp) :: prefac, sg
   
      real(rp) :: soc_p, soc_d
      real(rp) :: xi_p1, xi_p2
      real(rp) :: xi_d1, xi_d2
   
      complex(rp), dimension(2) :: rac
      complex(rp), dimension(9, 9) :: Lx, Ly, Lz
   
      real(rp) :: lz_loc
   
      !-----------------------------------------------------------------------
      ! Getting the angular momentum operators from math_mod.
      ! They are initially given in Cartesian coordinates.
      !-----------------------------------------------------------------------
      Lx(:, :) = L_x(:, :)
      Ly(:, :) = L_y(:, :)
      Lz(:, :) = L_z(:, :)
   
      !-----------------------------------------------------------------------
      ! Transform them into the spherical-harmonics basis.
      !-----------------------------------------------------------------------
      call hcpx(Lx, 'cart2sph')
      call hcpx(Ly, 'cart2sph')
      call hcpx(Lz, 'cart2sph')
   
      !-----------------------------------------------------------------------
      ! Initialize the L.S Hamiltonian.
      !-----------------------------------------------------------------------
      this%lsham(:, :, :) = cmplx(0.0_rp, 0.0_rp, kind=rp)
   
      !-----------------------------------------------------------------------
      ! Loop over atomic types.
      !-----------------------------------------------------------------------
      do k = 1, this%charge%lattice%ntype
   
         sg = cmplx(0.5_rp, 0.0_rp, kind=rp)
   
         !--------------------------------------------------------------------
         ! Read SOC parameters from the potential.
         !--------------------------------------------------------------------
         xi_p1 = this%charge%lattice%symbolic_atoms(k)%potential%xi_p(1)
         xi_p2 = this%charge%lattice%symbolic_atoms(k)%potential%xi_p(2)
   
         xi_d1 = this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1)
         xi_d2 = this%charge%lattice%symbolic_atoms(k)%potential%xi_d(2)
   
         !--------------------------------------------------------------------
         ! The two SOC parameters belonging to the same orbital channel
         ! must have the same sign.
         !--------------------------------------------------------------------
         if (xi_p1 * xi_p2 < 0.0_rp) then
            error stop &
               'build_lsham: xi_p(1) and xi_p(2) must have the same sign'
         end if
   
         if (xi_d1 * xi_d2 < 0.0_rp) then
            error stop &
               'build_lsham: xi_d(1) and xi_d(2) must have the same sign'
         end if
   
         !--------------------------------------------------------------------
         ! Signed geometric mean.
         !
         ! (+a,+b) --> +sqrt(a*b)
         ! (-a,-b) --> -sqrt(a*b)
         !
         ! This preserves the sign supplied in the potential file.
         !--------------------------------------------------------------------
         soc_p = sign(sqrt(abs(xi_p1 * xi_p2)), xi_p1)
         soc_d = sign(sqrt(abs(xi_d1 * xi_d2)), xi_d1)
   
         !--------------------------------------------------------------------
         ! WARNING:
         !
         ! The orbital-polarization part below is intentionally left
         ! unchanged from the original implementation.
         !
         ! If xi_d(1) is negative, expressions such as
         !
         !   sqrt(xi_d(1) * rac)
         !   sqrt(xi_d(1) * lmom(3))
         !
         ! may have negative real arguments and may therefore produce NaNs.
         !
         ! The negative-SOC modification above is therefore safe by itself,
         ! but calculations with orb_pol = .true. require separate analysis
         ! of the orbital-polarization terms.
         !--------------------------------------------------------------------
         if (this%orb_pol .and. xi_d1 < 0.0_rp) then
            write(*, '(a)') &
               'WARNING(build_lsham): negative xi_d with orb_pol=.true.;'
            write(*, '(a)') &
               'rac/lz_loc square roots are unchanged and may produce NaNs.'
         end if
   
         !--------------------------------------------------------------------
         ! Check if orbital polarization is enabled.
         !
         ! NOTE: This block is unchanged from the original implementation.
         !--------------------------------------------------------------------
         if (this%orb_pol) then
   
            rac = sqrt( &
               this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1) * &
               this%charge%lattice%symbolic_atoms(k)%potential%rac )
   
            lz_loc = sqrt( &
               this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1) * &
               this%charge%lattice%symbolic_atoms(k)%potential%lmom(3) )
   
         else
   
            rac = 0.0_rp
            lz_loc = 0.0_rp
   
         end if
   
         !--------------------------------------------------------------------
         ! Construct the L.S Hamiltonian.
         !--------------------------------------------------------------------
         do i = 1, 9
            do j = 1, 9
   
               ! Reset prefactor for every matrix element.
               prefac = cmplx(0.0_rp, 0.0_rp, kind=rp)
   
               ! p-orbital block
               if (i >= 2 .and. i <= 4 .and. &
                   j >= 2 .and. j <= 4) then
   
                  prefac = sg * soc_p
   
               end if
   
               ! d-orbital block
               if (i >= 5 .and. i <= 9 .and. &
                   j >= 5 .and. j <= 9) then
   
                  prefac = sg * soc_d
   
               end if
   
               ! H11
               this%lsham(j, i, k) = &
                  this%lsham(j, i, k) + &
                  prefac * Lz(j, i) + &
                  Lz(j, i) * rac(1) * lz_loc
   
               ! H12
               this%lsham(j, i + 9, k) = &
                  this%lsham(j, i + 9, k) + &
                  prefac * (Lx(j, i) - i_unit * Ly(j, i))
   
               ! H21
               this%lsham(j + 9, i, k) = &
                  this%lsham(j + 9, i, k) + &
                  prefac * (Lx(j, i) + i_unit * Ly(j, i))
   
               ! H22
               this%lsham(j + 9, i + 9, k) = &
                  this%lsham(j + 9, i + 9, k) - &
                  prefac * Lz(j, i) - &
                  Lz(j, i) * rac(2) * lz_loc
   
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
      this%tmat(:, :, :, :) = cmplx(0.0_rp, 0.0_rp)
      do k = 1, this%charge%lattice%ntype
         sg = cmplx(0.5_rp, 0.0_rp)
         soc_p = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_p(1)*this%charge%lattice%symbolic_atoms(k)%potential%xi_p(2))
         soc_d = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1)*this%charge%lattice%symbolic_atoms(k)%potential%xi_d(2))
         prefac = 0.0_rp
         do i = 1, 9
            do j = 1, 9
               if (i >= 2 .and. i <= 4 .and. j >= 2 .and. j <= 4) prefac = sg*soc_p
               if (i >= 5 .and. i <= 9 .and. j >= 5 .and. j <= 9) prefac = sg*soc_d
               ! build Tx
               this%tmat(j, i, 1, k) = this%tmat(j, i, 1, k) + prefac*i_unit*Ly(j, i)*2.0_rp ! Tx_11
               this%tmat(j, i + 9, 1, k) = this%tmat(j, i + 9, 1, k) - prefac*Lz(j, i)*2.0_rp*cone ! Tx_12
               this%tmat(j + 9, i, 1, k) = this%tmat(j + 9, i, 1, k) + prefac*Lz(j, i)*2.0_rp*cone ! Tx_21
               this%tmat(j + 9, i + 9, 1, k) = this%tmat(j + 9, i + 9, 1, k) - prefac*i_unit*Ly(j, i)*2.0_rp ! Tx_22
               ! build Ty
               this%tmat(j, i, 2, k) = this%tmat(j, i, 2, k) - prefac*i_unit*Lx(j, i)*2.0_rp ! Ty_11
               this%tmat(j, i + 9, 2, k) = this%tmat(j, i + 9, 2, k) + prefac*i_unit*Lz(j, i)*2.0_rp ! Ty_12
               this%tmat(j + 9, i, 2, k) = this%tmat(j + 9, i, 2, k) + prefac*i_unit*Lz(j, i)*2.0_rp ! Ty_21
               this%tmat(j + 9, i + 9, 2, k) = this%tmat(j + 9, i + 9, 2, k) + prefac*i_unit*Lx(j, i)*2.0_rp ! Ty_22
               ! build Tz
               this%tmat(j, i + 9, 3, k) = this%tmat(j, i + 9, 3, k) + prefac*(Lx(j, i) - i_unit*Ly(j, i))*2.0_rp*cone ! Tz_12
               this%tmat(j + 9, i, 3, k) = this%tmat(j + 9, i, 3, k) + prefac*(Lx(j, i) + i_unit*Ly(j, i))*2.0_rp*cmone ! Tz_21
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

      this%obarm = 0.d00

      do ntype = 1, this%lattice%ntype
         obm0 = cmplx(0.0d0); obm1 = cmplx(0.0d0)
         do m = 1, 9
            obm0(m, m) = this%lattice%symbolic_atoms(ntype)%potential%obx0(m)
            obm1(m, m) = this%lattice%symbolic_atoms(ntype)%potential%obx1(m)
         end do
         mom(:) = cmplx(this%lattice%symbolic_atoms(ntype)%potential%mom(:), 0.0d0)
         do m = 1, 9
            do l = 1, 9
               this%obarm(m, l, ntype) = obm0(m, l) + obm1(m, l)*mom(3)
               this%obarm(m + 9, l + 9, ntype) = obm0(m, l) - obm1(m, l)*mom(3)
               this%obarm(l, m + 9, ntype) = obm1(m, l)*mom(1) - i_unit*obm1(m, l)*mom(2)
               this%obarm(l + 9, m, ntype) = obm1(m, l)*mom(1) + i_unit*obm1(m, l)*mom(2)
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

      this%enim = 0.0d0

      do ntype = 1, this%lattice%ntype
         em0 = cmplx(0.0d0); em1 = cmplx(0.0d0)
         do m = 1, 9
            eu = this%lattice%symbolic_atoms(ntype)%potential%cx(m, 1) - this%lattice%symbolic_atoms(ntype)%potential%cex(m, 1)
            ed = this%lattice%symbolic_atoms(ntype)%potential%cx(m, 2) - this%lattice%symbolic_atoms(ntype)%potential%cex(m, 2)
            ex0(m) = 0.5*(eu + ed)
            ex1(m) = 0.5*(eu - ed)
            em0(m, m) = ex0(m)
            em1(m, m) = ex1(m)
         end do
         mom(:) = cmplx(this%lattice%symbolic_atoms(ntype)%potential%mom(:), 0.0d0)
         do m = 1, 9
            do l = 1, 9
               this%enim(m, l, ntype) = em0(m, l) + em1(m, l)*mom(3)
               this%enim(m + 9, l + 9, ntype) = em0(m, l) - em1(m, l)*mom(3)
               this%enim(l, m + 9, ntype) = em1(m, l)*mom(1) - i_unit*em1(m, l)*mom(2)
               this%enim(l + 9, m, ntype) = em1(m, l)*mom(1) + i_unit*em1(m, l)*mom(2)
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
      integer :: i, j, k, l, m, n, itype, ino, ja, jo, ji, nr, ia
      integer :: ntype

      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype) ! Atom number in clust
         ino = this%charge%lattice%num(ia) ! Atom bravais type of ia
         nr = this%charge%lattice%nn(ia, 1) ! Number of neighbours considered
         !write(123, *)´bulkham´
         call this%chbar_nc(ia, nr, ino, ntype)
         do m = 1, nr
            do i = 1, 9
               do j = 1, 9
                  this%ee(j, i, m, ntype) = this%hmag(j, i, m, 4) + this%hmag(j, i, m, 3)        ! H0+Hz
                  this%ee(j + 9, i + 9, m, ntype) = this%hmag(j, i, m, 4) - this%hmag(j, i, m, 3)        ! H0-Hz
                  this%ee(j, i + 9, m, ntype) = this%hmag(j, i, m, 1) - i_unit*this%hmag(j, i, m, 2) ! Hx-iHy
                  this%ee(j + 9, i, m, ntype) = this%hmag(j, i, m, 1) + i_unit*this%hmag(j, i, m, 2) ! Hx+iHy
               end do ! end of orbital j loop
            end do ! end of orbital i loop
            write(131,*) 'm=', m
            write(131,'(18f10.6)') real(this%ee(:,:,m,ntype))
            write(132,*) 'm=', m
            write(132,'(18f10.6)') aimag(this%ee(:,:,m,ntype))
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
               ! Check if neighbour ´m´ exists for atom ´ntype´, otherwise fill HoH Hamiltonian with zeros.
               if (ji > 0) then
                  call zgemm('n', 'n', 18, 18, 18, cone, this%ee(:, :, m, ntype), 18, this%obarm(:, :, ji), 18, czero, this%eeo(:, :, m, ntype), 18)
                  call zgemm('n', 'c', 18, 18, 18, cone, this%eeo(:, :, m, ntype), 18, this%ee(:, :, m, ntype), 18, czero, this%eeoee(:, :, m, ntype), 18)
               else
                  this%eeo(:, :, m, ntype) = 0.0d0
               end if
               !write(*,*) ´m=´, m
               !write(*,´(18f10.6)´) real(this%eeo(:,:,m,ntype))
               !write(*,*) ´ee´, m
               !write(*,´(18f10.6)´) real(this%ee(:,:,m,ntype))
            end do
         end if
      end do ! end of atom type number
      if (this%local_axis) then
         this%ee_glob = this%ee
         if (this%hoh) this%eeo_glob = this%eeo
      end if
      close (128)
   end subroutine build_bulkham

   subroutine build_locham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: it, ino, nr, nlim, m, i, j, ja, ji

      call g_timer%start('build local hamiltonian')
    !!$omp parallel do private(nlim, nr, ino, m, i, j, ji, ja, this)
      do nlim = 1, this%charge%lattice%nmax
         nr = this%charge%lattice%nn(nlim, 1) ! Number of neighbours considered
         ino = this%charge%lattice%num(nlim)
         call this%chbar_nc(nlim, nr, ino, nlim)
         do m = 1, nr
            do i = 1, 9
               do j = 1, 9
                  this%hall(j, i, m, nlim) = this%hmag(j, i, m, 4) + this%hmag(j, i, m, 3) ! H0+Hz
                  this%hall(j + 9, i + 9, m, nlim) = this%hmag(j, i, m, 4) - this%hmag(j, i, m, 3) ! H0-Hz
                  this%hall(j, i + 9, m, nlim) = this%hmag(j, i, m, 1) - i_unit*this%hmag(j, i, m, 2) ! Hx-iHy
                  this%hall(j + 9, i, m, nlim) = this%hmag(j, i, m, 1) + i_unit*this%hmag(j, i, m, 2) ! Hx+iHy
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
               ! Check if neighbour ´m´ exists for atom ´nlim´, otherwise fill HoH Hamiltonian with zeros.
               if (ji > 0) then
                  call zgemm('n', 'n', 18, 18, 18, cone, this%hall(1, 1, m, nlim), 18, this%obarm(1, 1, ji), 18, czero, this%hallo(1, 1, m, nlim), 18)
               else
                  this%hallo(:, :, m, nlim) = 0.0d0
               end if
            end do
         end if
      end do
    !!$omp end parallel do
      if (this%local_axis) then
         this%hall_glob = this%hall
         if (this%hoh) this%hallo_glob = this%hallo
      end if
      call g_timer%stop('build local hamiltonian')
   end subroutine build_locham

subroutine rs2pao(this)
   implicit none
   class(hamiltonian), intent(inout) :: this

   integer :: n_atoms, max_orbital
   integer :: ntype, ia, nr_i, jlab, ja, jtype, nr_j, klab, ka, ktype
   integer :: qlab_jk
   integer :: i, j, p, ipao, jpao
   integer :: max_pairs, npairs
   integer :: iostatus

   integer, allocatable :: pair_type(:)
   integer, allocatable :: pair_idx(:, :)
   complex(rp), allocatable :: pair_h(:, :, :)

   real(rp), dimension(3) :: rjk
   complex(rp), dimension(18,18) :: hij
   complex(rp), dimension(18,18) :: hjk
   complex(rp), dimension(18,18) :: himom
   complex(rp), dimension(18,18) :: pathblk
   complex(rp), dimension(18,18) :: dum

   n_atoms     = this%charge%lattice%ntype
   max_orbital = 9

   if (this%hoh) then
      call this%build_obarm()
      call this%build_enim()
   end if

   open(unit=92, file='rs2paoham.dat', action='write', status='replace', iostat=iostatus)
   if (iostatus /= 0) then
      write(*,*) 'rs2pao: could not open rs2paoham.dat'
      return
   end if

   do ntype = 1, this%charge%lattice%ntype

      ia   = this%charge%lattice%atlist(ntype)
      nr_i = this%charge%lattice%nn(ia, 1)

      ! Safe upper bound:
      ! original bare pairs + all double-loop i->j->k path destinations
      max_pairs = max(1, nr_i)
      if (this%hoh) then
         do jlab = 1, nr_i
            call get_atom_from_ee_label(ia, jlab, ja)
            if (ja /= 0) then
               max_pairs = max_pairs + max(1, this%charge%lattice%nn(ja, 1))
            end if
         end do
      end if

      allocate(pair_type(max_pairs))
      allocate(pair_idx(3, max_pairs))
      allocate(pair_h(18,18,max_pairs))

      pair_type(:)    = 0
      pair_idx(:, :)  = 0
      pair_h(:, :, :) = czero
      npairs          = 0

      ! ============================================================
      ! 1) Bare h-part: accumulate h_ik from the original one-hop set
      ! ============================================================
      do jlab = 1, nr_i
         call get_atom_from_ee_label(ia, jlab, ja)
         if (ja == 0) cycle

         hij(:, :) = this%ee(:, :, jlab, ntype)

         call accumulate_block_by_actual_pair(ia, ja, hij, &
                                              pair_type, pair_idx, pair_h, npairs)
      end do

      ! ============================================================
      ! 2) HoH part: for every path i -> j -> k, add
      !              - h_ij * obar_j * h_jk
      ! ============================================================
      if (this%hoh) then
         do jlab = 1, nr_i

            call get_atom_from_ee_label(ia, jlab, ja)
            if (ja == 0) cycle

            jtype = this%charge%lattice%iz(ja)
            if (jtype <= 0) cycle

            hij(:, :) = this%ee(:, :, jlab, ntype)

            ! himom = h_ij * obar_j
            call zgemm('n', 'n', 18, 18, 18, cone, hij, 18, &
                       this%obarm(:, :, jtype), 18, czero, himom, 18)

            nr_j = this%charge%lattice%nn(ja, 1)

            do klab = 1, nr_j

               call get_atom_from_ee_label(ja, klab, ka)
               if (ka == 0) cycle

               ktype = this%charge%lattice%iz(ka)
               if (ktype <= 0) cycle

               ! Actual vector of the second hop j -> k
               rjk(:) = this%charge%lattice%cr(:, ja) - this%charge%lattice%cr(:, ka)

               ! Find which representative ee(:, :, label, jtype)
               ! corresponds to the actual vector r_jk
               call find_ee_label_from_vector(jtype, rjk, qlab_jk)
               if (qlab_jk <= 0) cycle

               hjk(:, :) = this%ee(:, :, qlab_jk, jtype)

               ! pathblk = - (h_ij * obar_j) * h_jk
               call zgemm('n', 'n', 18, 18, 18, -cone, himom, 18, &
                          hjk, 18, czero, pathblk, 18)

               ! Accumulate into the FINAL effective pair (i,k),
               ! whose key is built from the ACTUAL rik = r_i - r_k
               call accumulate_block_by_actual_pair(ia, ka, pathblk, &
                                                    pair_type, pair_idx, pair_h, npairs)
            end do
         end do
      end if

      ! ============================================================
      ! 3) Add onsite-only one-site terms at the very end
      ! ============================================================
      call accumulate_block_by_actual_pair(ia, ia, this%lsham(:, :, ntype), &
                                           pair_type, pair_idx, pair_h, npairs)

      if (this%hoh) then
         call accumulate_block_by_actual_pair(ia, ia, this%enim(:, :, ntype), &
                                              pair_type, pair_idx, pair_h, npairs)
      end if

      ! ============================================================
      ! 4) Convert and print all final effective blocks
      ! ============================================================
      do p = 1, npairs
         dum(:, :) = pair_h(:, :, p)

         call hcpx(dum(1:9,   1:9   ), 'sph2cart')
         call hcpx(dum(1:9,   10:18 ), 'sph2cart')
         call hcpx(dum(10:18, 1:9   ), 'sph2cart')
         call hcpx(dum(10:18, 10:18 ), 'sph2cart')

         do i = 1, 18
            do j = 1, 18
               ipao = 0
               jpao = 0
               call site2orb(i, ia,           ipao, n_atoms, max_orbital)
               call site2orb(j, pair_type(p), jpao, n_atoms, max_orbital)
               write(92, '(3I4,2I7,2F22.14)') pair_idx(:, p), ipao, jpao, &
                    real(dum(i, j)) * ry2ev, aimag(dum(i, j)) * ry2ev
            end do
         end do
      end do

      deallocate(pair_type)
      deallocate(pair_idx)
      deallocate(pair_h)

   end do

   close(92)

contains

   subroutine get_atom_from_ee_label(center_atom, ee_label, atom_out)
      implicit none
      integer, intent(in)  :: center_atom, ee_label
      integer, intent(out) :: atom_out

      ! Same convention as your original rs2pao:
      !   ee(:,:,1, type_of_center) = onsite block of THIS ACTUAL center atom
      !   ee(:,:,label>1, ...)      = neighbor block associated with nn(center_atom, label)
      if (ee_label == 1) then
         atom_out = center_atom
      else
         atom_out = this%charge%lattice%nn(center_atom, ee_label)
      end if
   end subroutine get_atom_from_ee_label


   subroutine build_idx_from_actual_pair(icenter, kactual, ktype_out, idx_out)
      implicit none
      integer, intent(in)  :: icenter, kactual
      integer, intent(out) :: ktype_out
      integer, intent(out) :: idx_out(3)

      real(rp), parameter :: tol = 1.0d-3
      real(rp), dimension(3) :: rik, rijtest
      integer :: repk_atom
      integer :: idxi, idxj, idxk

      ! IMPORTANT:
      ! rik is built from ACTUAL atom numbers in the cluster
      ! ktype is only used afterwards to select the representative atom
      ktype_out = this%charge%lattice%iz(kactual)
      idx_out(:) = 0

      repk_atom = this%charge%lattice%atlist(ktype_out)

      rik(:) = this%charge%lattice%cr(:, icenter) - this%charge%lattice%cr(:, kactual)

      do idxi = -10, 10
         do idxj = -10, 10
            do idxk = -10, 10
               rijtest(:) = this%charge%lattice%cr(:, icenter) - &
                            ( this%charge%lattice%cr(:, repk_atom) &
                            + idxi * this%charge%lattice%a(:, 1) &
                            + idxj * this%charge%lattice%a(:, 2) &
                            + idxk * this%charge%lattice%a(:, 3) )

               if (norm2(rik(:) - rijtest(:)) < tol) then
                  idx_out(:) = (/ idxi, idxj, idxk /)
                  return
               end if
            end do
         end do
      end do

      write(*,'(a,2i8)') 'rs2pao warning: idx not found for actual pair ', icenter, kactual
   end subroutine build_idx_from_actual_pair


   subroutine find_ee_label_from_vector(center_type, r_target, label_out)
      implicit none
      integer, intent(in)  :: center_type
      real(rp), intent(in) :: r_target(3)
      integer, intent(out) :: label_out

      real(rp), parameter :: tol = 1.0d-3
      integer :: center_rep, nrrep, qlab, qatom
      real(rp), dimension(3) :: rrep

      label_out  = 0
      center_rep = this%charge%lattice%atlist(center_type)
      nrrep      = this%charge%lattice%nn(center_rep, 1)

      do qlab = 1, nrrep
         call get_atom_from_ee_label(center_rep, qlab, qatom)
         if (qatom == 0) cycle

         ! Vector represented by ee(:, :, qlab, center_type)
         rrep(:) = this%charge%lattice%cr(:, center_rep) - this%charge%lattice%cr(:, qatom)

         if (norm2(rrep(:) - r_target(:)) < tol) then
            label_out = qlab
            return
         end if
      end do

      write(*,'(a,i8,3f12.6)') 'rs2pao warning: label not found for center_type / r_target = ', &
                               center_type, r_target(1), r_target(2), r_target(3)
   end subroutine find_ee_label_from_vector


   subroutine accumulate_block_by_actual_pair(icenter, kactual, block_in, &
                                              pair_type_arr, pair_idx_arr, pair_h_arr, npairs_inout)
      implicit none
      integer, intent(in) :: icenter, kactual
      complex(rp), intent(in) :: block_in(18,18)
      integer, intent(inout) :: pair_type_arr(:)
      integer, intent(inout) :: pair_idx_arr(:, :)
      complex(rp), intent(inout) :: pair_h_arr(:, :, :)
      integer, intent(inout) :: npairs_inout

      integer :: ktype_loc, idx_loc(3)
      integer :: q, nmax

      call build_idx_from_actual_pair(icenter, kactual, ktype_loc, idx_loc)

      do q = 1, npairs_inout
         if (pair_type_arr(q) == ktype_loc) then
            if (all(pair_idx_arr(:, q) == idx_loc(:))) then
               pair_h_arr(:, :, q) = pair_h_arr(:, :, q) + block_in(:, :)
               return
            end if
         end if
      end do

      nmax = size(pair_type_arr)
      if (npairs_inout >= nmax) then
         write(*,*) 'rs2pao error: pair buffer too small.'
         stop
      end if

      npairs_inout = npairs_inout + 1
      pair_type_arr(npairs_inout)   = ktype_loc
      pair_idx_arr(:, npairs_inout) = idx_loc(:)
      pair_h_arr(:, :, npairs_inout)= block_in(:, :)
   end subroutine accumulate_block_by_actual_pair

end subroutine rs2pao

!   subroutine rs2pao(this)
!      implicit none
!      class(hamiltonian), intent(inout) :: this
!      ! Local variables
!      real(rp), dimension(3) :: rij, rijtest
!      integer :: i, j, k, l, idxi, idxj, idxk, itype, ino, ja, jo, ji, nr, ia, iia, jja, ipao, jpao
!      integer :: jj, jt, max_orbital, n_atoms
!      integer :: ntype, iostat1, iostat2, iostatus
!      real(rp), dimension(3) :: vet, vetpao, idx
!      real(rp), dimension(3, 3) :: a_inv
!      complex(rp), dimension(18, 18) :: dum
!      n_atoms = this%charge%lattice%ntype
!      max_orbital = 9 
!
!      open (unit=92, file='rs2paoham.dat', action='write', iostat=iostatus, status='replace')
!      do ntype = 1, this%charge%lattice%ntype
!         ia = this%charge%lattice%atlist(ntype) ! Atom number in clust
!         ino = this%charge%lattice%num(ia) ! Atom bravais type of ia
!         nr = this%charge%lattice%nn(ia, 1) ! Number of neighbours considered
!         do k = 1, nr
!            jj = this%charge%lattice%nn(ia, k)
!            !write(123, *)´ia, ii´, ia, m, this%charge%lattice%nn(ia, m)
!            if (k == 1) then
!               jj = ia
!            end if
!            if (jj /= 0) then
!               rij(:) = this%charge%lattice%cr(:, ia) - this%charge%lattice%cr(:, jj)
!
!               rijtest(:) = 0.0d0
!               do idxi = -5, 5
!                  do idxj = -5, 5
!                     do idxk = -5, 5
!                        rijtest(:) = this%charge%lattice%cr(:, ia) - (this%charge%lattice%cr(:, this%charge%lattice%iz(jj)) &
!                                                                      + idxi*this%charge%lattice%a(:, 1) + idxj*this%charge%lattice%a(:, 2) + idxk*this%charge%lattice%a(:, 3))
!                        if (norm2(rij(:) - rijtest(:)) < 1.0d-3) then
!                           idx(:) = [idxi, idxj, idxk]
!                        end if
!                     end do
!                  end do
!               end do
!
!               if (k == 1) this%ee(:, :, k, ntype) = this%ee(:, :, k, ntype) + this%lsham(:, :, ntype) !+ this%enim(:,:,ntype)
!
!               call hcpx(this%ee(1:9,1:9,k,ntype), 'sph2cart')
!               call hcpx(this%ee(1:9,10:18,k,ntype), 'sph2cart')
!               call hcpx(this%ee(10:18,1:9,k,ntype), 'sph2cart')
!               call hcpx(this%ee(10:18,10:18,k,ntype), 'sph2cart')
!               do i = 1, 18
!                  do j = 1, 18
!                     ipao = 0; jpao = 0
!                     call site2orb(i, ia, ipao, n_atoms, max_orbital)
!                     call site2orb(j, this%charge%lattice%iz(jj), jpao, n_atoms, max_orbital)
!                     write (92, '(3I4,2I7,2F22.14)') int(idx(:)), ipao, jpao, real(this%ee(i, j, k, ntype))*ry2ev, aimag(this%ee(i, j, k, ntype))*ry2ev
!                  end do
!               end do
!            end if
!         end do
!      end do
!      close (92)
!   end subroutine rs2pao

   subroutine build_from_paoflow_opt(this)
      implicit none
      type hamData
         integer :: idxi, idxj, idxk
         integer :: orbl, orbm
         real :: dumre, dumcmplx
      end type hamData
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, itype, ino, ja, jo, ji, nr, ia, iia, jja
      integer :: jj, jt, orbl, orbm, idxi, idxj, idxk
      integer :: ntype, iostat1, iostat2, iostatus, n_atoms, max_orbital, numLines
      real(rp), dimension(3) :: vet, vetpao
      real(rp) :: dumre, dumcmplx
      integer, dimension(maxval(this%charge%lattice%nn(:, 1)) + 1, 3) :: idxup, idxdw, idx
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
      idx = 0
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
                                                                   + ham%idxi*this%charge%lattice%a(:, 1) + ham%idxj*this%charge%lattice%a(:, 2) + ham%idxk*this%charge%lattice%a(:, 3))!*this%charge%lattice%alat
                     if (norm2(vet(:) - vetpao(:)) < 1.0d-3) then
                        idx(k, :) = [ham%idxi, ham%idxj, ham%idxk]
                        this%ee(i, j, k, ntype) = cmplx(ham%dumre, ham%dumcmplx)/13.605703976
                        ! if(ntype==1.or.ntype==2)then
                        !   if(i==j.and.k==1.and.i>=5.and.i<=9)then
                        !     this%ee(i, j, k, ntype) = cmplx(ham%dumre-1.113 , ham%dumcmplx)/13.605703976
                        !   else if(i==j.and.k==1.and.i>=14.and.i<=18)then
                        !     this%ee(i, j, k, ntype) = cmplx(ham%dumre+1.113 , ham%dumcmplx)/13.605703976
                        !   end if
                        ! end if
                     end if
                  end if
               end do
            end if
            write (128, *) 'm=', k, 'Atom=', jj, 'Coordinates=', this%charge%lattice%cr(:, jj), 'Ntype=', ntype, 'Index=', idx(k, :)
            write (129, *) 'm=', k, 'Atom=', jj, 'Coordinates=', this%charge%lattice%cr(:, jj), 'Ntype=', ntype, 'Index=', idx(k, :)
            write (128, '(18f10.6)') real(this%EE(1:18, 1:18, k, ntype))*13.605703976
            write (129, '(18f10.6)') aimag(this%EE(1:18, 1:18, k, ntype))*13.605703976
            write (128, *) sum(real(this%ee(:, :, k, ntype)))
            write (129, *) sum(real(this%ee(:, :, k, ntype)))
         end do
      end do
      !$omp end parallel do
      deallocate (hamArray)
      call g_timer%stop('Hamiltonian allocation')
   end subroutine build_from_paoflow_opt

   subroutine build_from_paoflow(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, itype, ino, ja, jo, ji, nr, ia, iia, jja
      integer :: jj, jt, orbl, orbm, idxi, idxj, idxk
      integer :: ntype, iostat1, iostat2, iostatus, n_atoms, max_orbital
      real(rp), dimension(3) :: vet, vetpao, cri_dir, crj_dir, cri_cart, crj_cart
      integer, dimension(maxval(this%charge%lattice%nn(:, 1)) + 1, 3) :: idxup, idxdw, idx
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
            !write(123, *)´ia, ii´, ia, m, this%charge%lattice%nn(ia, m)
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
                  if (orbl <= n_atoms*max_orbital) then
                     i = modulo(orbl - 1, max_orbital) + 1
                     iia = int((orbl - 1)/max_orbital) + 1
                  else
                     i = modulo(orbl - 1, max_orbital) + 10
                     iia = int((orbl - 1 - n_atoms*max_orbital)/max_orbital) + 1
                  end if

                  if (orbm <= n_atoms*max_orbital) then
                     j = modulo(orbm - 1, max_orbital) + 1
                     jja = int((orbm - 1)/max_orbital) + 1
                  else
                     j = modulo(orbm - 1, max_orbital) + 10
                     jja = int((orbm - 1 - n_atoms*max_orbital)/max_orbital) + 1
                  end if

                  !cri_cart(:) = this%charge%lattice%cr(:, iia)
                  !crj_cart(:) = this%charge%lattice%cr(:, jja)

                  !cri_dir(:) = cartesian_to_direct(this%charge%lattice%a,cri_cart)
                  !crj_dir(:) = cartesian_to_direct(this%charge%lattice%a,crj_cart)
                  if (iia == ia) then
                     !vetpao(:) = cri_dir(:) - (crj_dir(:) + [idxi,idxj,idxk])

                     vetpao(:) = this%charge%lattice%cr(:, iia) - (this%charge%lattice%cr(:, jja) &
                                                                   + idxi*this%charge%lattice%a(:, 1) + idxj*this%charge%lattice%a(:, 2) + idxk*this%charge%lattice%a(:, 3))!*this%charge%lattice%alat
                     if (norm2(vet(:) - vetpao(:)) < 1.0d-3) then
                        idx(k, :) = [idxi, idxj, idxk]
                        this%ee(i, j, k, ntype) = cmplx(dumre, dumcmplx)/13.605703976
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
            !write(129,*)´m=´,k, ´Atom=´, jj, ´Coordinates=´, this%charge%lattice%cr(:, jj), ´Ntype=´,ntype, ´Index=´, idx(k,:)
            !write(129,´(9f10.6)´) real(this%EE(10:18,1:9,k,ntype))*13.605703976
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

      this%hhmag(:, :, :) = 0.0d0

      vv = norm2(vet)

      ! Real to complex
      dot = cmplx(dot_product(this%charge%lattice%symbolic_atoms(it)%potential%mom, this%charge%lattice%symbolic_atoms(jt)%potential%mom), kind=kind(0.0d0))
      do i = 1, this%charge%lattice%ntype
         do j = 1, 3
            momc(i, j) = cmplx(this%charge%lattice%symbolic_atoms(i)%potential%mom(j), kind=kind(0.0d0))
         end do
      end do
      cross = cmplx(cross_product(this%charge%lattice%symbolic_atoms(it)%potential%mom, this%charge%lattice%symbolic_atoms(jt)%potential%mom), kind=kind(0.0d0))
      hhhc(:, :) = cmplx(hhh(:, :), kind=kind(0.0d0))

      do ilm = 1, 9
         do jlm = 1, 9
            this%hhmag(ilm, jlm, 4) = &
               this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx0(jlm) + &
               this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm)*dot
         end do
      end do

!    do ilm=1, 9
!      write(123, ´(9f10.6)´) (real(this%hhmag(ilm, jlm, 4)), jlm=1, 9)
!    end do

      if (vv <= 0.01d0) then
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
                  (this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx0(jlm))*momc(it, m) + &
                  (this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm))*momc(jt, m) + &
                  i_unit*this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm)*cross(m)
            end do
         end do
      end do

      if (vv > 0.01d0) return
      do m = 1, 3
         do ilm = 1, 9
            if (this%hoh) then
               this%hhmag(ilm, ilm, m) = this%hhmag(ilm, ilm, m) + this%charge%lattice%symbolic_atoms(it)%potential%cex1(ilm)*momc(it, m)
            else
               this%hhmag(ilm, ilm, m) = this%hhmag(ilm, ilm, m) + this%charge%lattice%symbolic_atoms(it)%potential%cx1(ilm)*momc(it, m)
            end if
         end do
      end do

      !do m=1, 3
      !  write(123, *)´m=´, m
      !  do ilm=1, 9
      !    write(123, ´(9f10.6)´) (real(this%hhmag(ilm, jlm, m)), jlm=1, 9)
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
      integer :: i, j, k, l, m, n, it, jt, jj, dummy
      integer :: ni, mdir
      integer :: kk ! Clust size number

      this%hmag(:, :, :, :) = 0.0d0

      r2 = this%charge%lattice%r2
      cralat(:, :) = this%charge%lattice%cr(:, :)*this%charge%lattice%alat
      kk = this%charge%lattice%kk

      call this%charge%lattice%clusba(r2, cralat, ia, kk, kk, dummy)

      !do m=1, nr
      !  print ´(9f10.6)´, real(this%charge%lattice%sbar(:, :, m, ino))
      !end do
      it = this%charge%lattice%iz(ia)
      do m = 1, nr
         jj = this%charge%lattice%nn(ia, m)
         !write(123, *)´ia, ii´, ia, m, this%charge%lattice%nn(ia, m)
         if (m == 1) then
            jj = ia
         end if
         if (jj /= 0) then
            jt = this%charge%lattice%iz(jj)
            if (this%lattice%pbc) then
               call this%lattice%f_wrap_coord_diff(this%lattice%kk,this%lattice%cr*this%lattice%alat,ia,jj,vet)
            else
               vet(:) = (this%charge%lattice%cr(:, jj) - this%charge%lattice%cr(:, ia))*this%charge%lattice%alat
            end if
            !write(123, ´(3f10.6)´) vet(:)
            !write(123, ´(3f10.6)´) this%charge%lattice%sbarvec(:, m)
            !write(123, ´(a, 3i4, 3f10.6)´) ´nn ´, IA, m, JJ, VET(:)
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
      !  write(123, *)´m=´, m
      !  do mdir=1, 4
      !    write(123, *)´mdir=´, mdir
      !    do i=1, 9
      !      write(123, ´(9f10.4)´)(real(this%hmag(i, j, m, mdir)), j=1, 9)
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

      eps = 0.0001d0
      ni = 1
      a1 = 0.0d0
      a2 = 0.0d0
      a3 = 0.0d0
      aaa = 0.0d0
      do i = 1, nr
         !write(123, ´(a, i4, 3f10.4)´)´i´, i, this%charge%lattice%sbarvec(:, i)
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
      !    write(123, ´(9f10.6)´)(hhh(ilm, jlm), jlm=1, 9)
      !end do
   end subroutine hmfind

   subroutine orb2site(orb, i_out, ia_out, n_atoms, max_orbital)
      integer, intent(in) :: orb, n_atoms, max_orbital
      integer, intent(out) :: i_out, ia_out

      if (orb <= n_atoms*max_orbital) then
         i_out = modulo(orb - 1, max_orbital) + 1
         ia_out = int((orb - 1)/max_orbital) + 1
      else
         i_out = modulo(orb - 1, max_orbital) + max_orbital + 1
         ia_out = int((orb - 1 - n_atoms*max_orbital)/max_orbital) + 1
      end if
   end subroutine orb2site

   subroutine site2orb(i_in, ia_in, orb_out, n_atoms, max_orbital)
      integer, intent(in) :: i_in, ia_in, n_atoms, max_orbital
      integer, intent(out) :: orb_out

      if (i_in <= max_orbital) then
         orb_out = (ia_in - 1)*max_orbital + i_in
      else
         orb_out = (ia_in - 1)*max_orbital + i_in + (n_atoms - 1)*max_orbital
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
         sdim = product(shape(this%hall))/18/18
         call rotmag_loc(this%hall, this%hall_glob, sdim, m_loc)
         sdim = product(shape(this%ee))/18/18
         call rotmag_loc(this%ee, this%ee_glob, sdim, m_loc)
         if (this%hoh) then
            sdim = product(shape(this%eeo))/18/18
            call rotmag_loc(this%eeo, this%eeo_glob, sdim, m_loc)
            sdim = product(shape(this%hallo))/18/18
            call rotmag_loc(this%hallo, this%hallo_glob, sdim, m_loc)
            sdim = product(shape(this%enim))/18/18
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
      integer :: sdim
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
