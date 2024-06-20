!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Mix
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
!> lorem ipsum
!------------------------------------------------------------------------------

module mix_mod
   use control_mod
   use lattice_mod
   use charge_mod
   use symbolic_atom_mod, only: symbolic_atom
   use precision_mod, only: rp
   use mpi_mod
   use string_mod
   use logger_mod, only: g_logger
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   !> Module's main structure
   type, public :: mix
      ! Lattice
      class(lattice), pointer :: lattice
      ! Charge
      class(charge), pointer :: charge
      ! Symbolic atom
      class(symbolic_atom), dimension(:), pointer :: symbolic_atom
      !> Description
      !>
      !> Description
      integer :: var
      !> Variables to save parameters for mixing
      real(rp), dimension(:, :), allocatable :: qia, qia_new, qia_old, qiaprev
      !> Variables to save magnetic moments for mixing
      real(rp), dimension(:, :), allocatable :: mag_old, mag_new, mag_mix
      !> Mixing parameter
      real(rp) :: beta
      !> Magnetic mixing parameter
      real(rp), dimension(:), allocatable :: magbeta
      !> Difference between interactions
      real(rp) :: delta
      !> Type of mixing. Can be linear or broyden
      character(len=7) :: mixtype
      !> Variables for Broyden mixing
      real(rp), dimension(:), allocatable :: v_broy, u_broy, fo_broy, muo_broy
      real(rp) :: fsqo
      integer :: itr, nmix
      !> Variable to managed induced magnetic moments
      logical, dimension(:), allocatable :: is_induced
   contains
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: print_state
      procedure :: save_to
      procedure :: mixpq
      procedure :: mix_magnetic_moments
      final :: destructor
   end type mix

   interface mix
      procedure :: constructor
   end interface mix

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] lattice_obj pointer
   !> @return type(mix)
   !---------------------------------------------------------------------------
   function constructor(lattice_obj, charge_obj) result(obj)
      type(mix) :: obj
      type(lattice), target, intent(in) :: lattice_obj
      type(charge), target, intent(in) :: charge_obj

      obj%lattice => lattice_obj
      obj%charge => charge_obj
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
      type(mix) :: this

#ifdef USE_SAFE_ALLOC
      if (allocated(this%qia)) call g_safe_alloc%deallocate('mix.qia', this%qia)
      if (allocated(this%qia_new)) call g_safe_alloc%deallocate('mix.qia_new', this%qia_new)
      if (allocated(this%qia_old)) call g_safe_alloc%deallocate('mix.qia_old', this%qia_old)
      if (allocated(this%qiaprev)) call g_safe_alloc%deallocate('mix.qiaprev', this%qiaprev)
      if (allocated(this%v_broy)) call g_safe_alloc%deallocate('mix.v_broy', this%v_broy)
      if (allocated(this%u_broy)) call g_safe_alloc%deallocate('mix.u_broy', this%u_broy)
      if (allocated(this%fo_broy)) call g_safe_alloc%deallocate('mix.fo_broy', this%fo_broy)
      if (allocated(this%muo_broy)) call g_safe_alloc%deallocate('mix.muo_broy', this%muo_broy)
      if (allocated(this%magbeta)) call g_safe_alloc%deallocate('mix.magbeta', this%magbeta)
      if (allocated(this%mag_new)) call g_safe_alloc%deallocate('mix.mag_new', this%mag_new)
      if (allocated(this%mag_old)) call g_safe_alloc%deallocate('mix.mag_old', this%mag_old)
      if (allocated(this%mag_mix)) call g_safe_alloc%deallocate('mix.mag_mix', this%mag_mix)
      if (allocated(this%is_induced)) call g_safe_alloc%deallocate('mix.is_induced', this%is_induced)
#else
      if (allocated(this%qia)) deallocate (this%qia)
      if (allocated(this%qia_new)) deallocate (this%qia_new)
      if (allocated(this%qia_old)) deallocate (this%qia_old)
      if (allocated(this%qiaprev)) deallocate (this%qiaprev)
      if (allocated(this%v_broy)) deallocate (this%v_broy)
      if (allocated(this%u_broy)) deallocate (this%u_broy)
      if (allocated(this%fo_broy)) deallocate (this%fo_broy)
      if (allocated(this%muo_broy)) deallocate (this%muo_broy)
      if (allocated(this%magbeta)) deallocate (this%magbeta)
      if (allocated(this%mag_new)) deallocate (this%mag_new)
      if (allocated(this%mag_old)) deallocate (this%mag_old)
      if (allocated(this%mag_mix)) deallocate (this%mag_mix)
      if (allocated(this%is_induced)) deallocate (this%is_induced)
#endif
   end subroutine destructor

   ! Member functions
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Read parameters from input file
   !
   !> @param[in] fname Namelist file
   !---------------------------------------------------------------------------
   subroutine build_from_file(this)
      class(mix), intent(inout) :: this
      character(len=sl) :: fname
      ! variables associated with the reading processes
      integer :: iostatus, funit

      include 'include_codes/namelists/mix.f90'

      fname = this%lattice%control%fname
      var = this%var
      mixtype = this%mixtype
      beta = this%beta
      call move_alloc(this%magbeta, magbeta)

      open (newunit=funit, file=fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//fmt('A', fname)//' not found', __FILE__, __LINE__)
      end if

      read (funit, nml=mix, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%error('Error while reading namelist', __FILE__, __LINE__)
         call g_logger%error('iostatus = '//fmt('I0', iostatus), __FILE__, __LINE__)
      end if
      close (funit)

      this%var = 0
      this%mixtype = mixtype
      this%beta = beta
      call move_alloc(magbeta, this%magbeta)
   end subroutine build_from_file

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this)
      implicit none
      class(mix), intent(inout) :: this

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('mix.qia', this%qia, (/this%lattice%nrec, 18/))
      call g_safe_alloc%allocate('mix.qia_new', this%qia_new, (/this%lattice%nrec, 18/))
      call g_safe_alloc%allocate('mix.qia_old', this%qia_old, (/this%lattice%nrec, 18/))
      call g_safe_alloc%allocate('mix.qiaprev', this%qiaprev, (/this%lattice%nrec, 18/))
      call g_safe_alloc%allocate('mix.v_broy', this%v_broy, (/this%lattice%nrec*18/))
      call g_safe_alloc%allocate('mix.u_broy', this%u_broy, (/this%lattice%nrec*18/))
      call g_safe_alloc%allocate('mix.fo_broy', this%fo_broy, (/this%lattice%nrec*18/))
      call g_safe_alloc%allocate('mix.muo_broy', this%muo_broy, (/this%lattice%nrec*18/))
      call g_safe_alloc%allocate('mix.magbeta', this%magbeta, (/this%lattice%nrec/))
      call g_safe_alloc%allocate('mix.mag_old', this%mag_old, (/this%lattice%nrec, 3/))
      call g_safe_alloc%allocate('mix.mag_new', this%mag_new, (/this%lattice%nrec, 3/))
      call g_safe_alloc%allocate('mix.mag_mix', this%mag_mix, (/this%lattice%nrec, 3/))
      call g_safe_alloc%allocate('mix.is_induced', this%is_induced, (/this%lattice%nrec/))
#else
      allocate (this%qia(this%lattice%nrec, 18), this%qia_new(this%lattice%nrec, 18), this%qia_old(this%lattice%nrec, 18))
      allocate (this%qiaprev(this%lattice%nrec, 18))
      allocate (this%v_broy(this%lattice%nrec*18), this%u_broy(this%lattice%nrec*18), this%fo_broy(this%lattice%nrec*18), this%muo_broy(this%lattice%nrec*18))
      allocate (this%magbeta(this%lattice%nrec), this%mag_old(this%lattice%nrec, 3), this%mag_new(this%lattice%nrec, 3), this%mag_mix(this%lattice%nrec, 3))
      allocate (this%is_induced(this%lattice%nrec))
#endif

      this%v_broy(:) = 0.0d0
      this%u_broy(:) = 0.0d0
      this%fo_broy(:) = 0.0d0
      this%muo_broy(:) = 0.0d0
      this%fsqo = 1.0d0
      this%itr = 0
      this%var = 0
      this%beta = 0.1d0
      this%magbeta(:) = 1.0d0
      this%nmix = 2
      this%mixtype = 'linear'
      this%is_induced = .false.
   end subroutine restore_to_default

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
      class(mix), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file
      integer :: newunit

      include 'include_codes/namelists/mix.f90'

      var = this%var

      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         write (unit, nml=mix)
      else if (present(file)) then
         open (unit=newunit, file=file)
         write (newunit, nml=mix)
         close (newunit)
      else
         write (*, nml=mix)
      end if
      close (newunit)
   end subroutine print_state

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Save potential parameters and band moments.
   !>
   !> Save potential parameters and band moments. Save the parameters read from the &par namelist into qia variables. To be used to mix later.
   !> @param[in] dummy variable of type mix
   !> @return type(mix)
   !---------------------------------------------------------------------------
   subroutine save_to(this, whereto)
      class(mix), intent(inout) :: this
      character(len=*), intent(in) :: whereto
      integer :: i, it

      select case (whereto)

      case ('old')
         do it = 1, this%lattice%nrec
            do I = 1, 3
               this%qia_old(IT, I) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(1, I - 1, 1)
               this%qia_old(IT, I + 6) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(3, I - 1, 1)
               this%qia_old(IT, I + 12) = this%symbolic_atom(this%lattice%nbulk + it)%potential%pl(I - 1, 1) !- (0.5d0 + INT(PL(I, 1)))
               !QI_OLD(IT, I+12) = ENU(I, 1)
               this%qia_old(IT, I + 3) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(1, I - 1, 2)
               this%qia_old(IT, I + 9) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(3, I - 1, 2)
               this%qia_old(IT, I + 15) = this%symbolic_atom(this%lattice%nbulk + it)%potential%pl(I - 1, 2) !- (0.5d0 + INT(PL(I, 2)))
               !QI_OLD(IT, I+15) = ENU(I, 2)
            end do
         end do
      case ('new')
         do it = 1, this%lattice%nrec
            do I = 1, 3
               this%qia_new(it, i) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(1, i - 1, 1)
               this%qia_new(it, i + 6) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(3, i - 1, 1)
               this%qia_new(it, i + 12) = this%symbolic_atom(this%lattice%nbulk + it)%potential%pl(i - 1, 1) !- (0.5d0 + INT(PL(I, 1)))
               !QI_OLD(IT, I+12) = ENU(I, 1)
               this%qia_new(it, i + 3) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(1, i - 1, 2)
               this%qia_new(it, i + 9) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(3, i - 1, 2)
               this%qia_new(it, i + 15) = this%symbolic_atom(this%lattice%nbulk + it)%potential%pl(i - 1, 2) !- (0.5d0 + INT(PL(I, 2)))
               !QI_OLD(IT, I+15) = ENU(I, 2)
            end do
         end do

      case ('prev')
         do it = 1, this%lattice%nrec
            do I = 1, 3
               this%qiaprev(it, i) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(1, i - 1, 1)
               this%qiaprev(it, i + 6) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(3, i - 1, 1)
               this%qiaprev(it, i + 12) = this%symbolic_atom(this%lattice%nbulk + it)%potential%pl(i - 1, 1) !- (0.5d0 + INT(PL(I, 1)))
               !QI_OLD(IT, I+12) = ENU(I, 1)
               this%qiaprev(it, i + 3) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(1, i - 1, 2)
               this%qiaprev(it, i + 9) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(3, i - 1, 2)
               this%qiaprev(it, i + 15) = this%symbolic_atom(this%lattice%nbulk + it)%potential%pl(i - 1, 2) !- (0.5d0 + INT(PL(I, 2)))
               !QI_OLD(IT, I+15) = ENU(I, 2)
            end do
         end do

      case ('current')
         do it = 1, this%lattice%nrec
            do I = 1, 3
               this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(1, i - 1, 1) = this%qia(it, i)
               this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(3, i - 1, 1) = this%qia(it, i + 6)
               this%symbolic_atom(this%lattice%nbulk + it)%potential%pl(i - 1, 1) = this%qia(it, i + 12)
               this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(1, i - 1, 2) = this%qia(it, i + 3)
               this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(3, i - 1, 2) = this%qia(it, i + 9)
               this%symbolic_atom(this%lattice%nbulk + it)%potential%pl(i - 1, 2) = this%qia(it, i + 15)
            end do
         end do
      end select
   end subroutine save_to

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Mix the magnetic moments.
   !>
   !> Mix the magnetic moments.
   !> @param[in] magnetic moments read from the potential file
   !> @return type(symbolic_atom): new magnetic moments after mixing.
   !---------------------------------------------------------------------------
   subroutine mix_magnetic_moments(this, mag_old, mag_new, mag_mix, mtot)
      implicit none
      class(mix), intent(inout) :: this
      ! Input variables
      real(rp), dimension(:, :), intent(in) :: mag_old, mag_new
      real(rp), dimension(:), intent(in) :: mtot
      ! Output variable
      real(rp), dimension(:, :), intent(out) :: mag_mix

      ! Local variable
      integer :: ia ! Self-consistent atom index

      ! Calculate the mixed magnetic moments
      do ia = 1, this%lattice%nrec
         if (mtot(ia + this%lattice%nbulk) < 0.5d0) then
            this%is_induced(ia) = .true.
            this%magbeta(ia) = 0.0d0
            if (rank == 0) call g_logger%info('Spin moment at atom '//fmt('i4', (ia))//' is being considered induced', __FILE__, __LINE__)
         end if
         mag_mix(ia, :) = (1.d0 - this%magbeta(ia))*mag_old(ia, :) + this%magbeta(ia)*mag_new(ia, :)
      end do
   end subroutine mix_magnetic_moments

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Mix the parameters QL and PL. Can be linear, andersen or broyden.
   !>
   !> Mix the parameters QL and PL. Can be linear, andersen or broyden.
   !> @param[in] dummy variable of type mix
   !> @param[in] old and new potential variables qia_old and qia_new, respectively.
   !> @return type(mix): qia, the mix potential variable.
   !---------------------------------------------------------------------------
   subroutine mixpq(this, qia_old, qia_new)
      use mpi_mod
      class(mix), intent(inout) :: this
      real(rp), dimension(this%lattice%nrec, 18), intent(in) :: qia_old, qia_new
      real(rp), dimension(this%lattice%nrec, 18) :: qi_to, qi_tn ! Local variables for broyden mixing
      real(rp) :: delta_atom
      integer :: ia ! Atom index
      logical :: reset
      real(rp) :: Bnorm

      select case (trim(this%mixtype))
      case ('linear')
         this%qia(:, :) = (1.0d0 - this%beta)*qia_old(:, :) + this%beta*qia_new(:, :)
      case ('broyden')
         ! call g_logger%fatal('Broyden mixing not implemented yet!', __FILE__, __LINE__)
         reset = .false.
         qi_to(:, :) = this%qia_old(:, :)
         qi_tn(:, :) = this%qia_new(:, :)
         call broydn(this%beta, this%beta, reset, qi_to, qi_tn, Bnorm, this%lattice%nrec*18, this%itr, this%fsqo, this%u_broy, this%v_broy, this%muo_broy, this%fo_broy, this%nmix)
         this%qia(:, :) = qi_to(:, :)
         bnorm = bnorm**0.5d0
      end select
      this%charge%dq(:) = 0.0d0
      do ia = 1, this%lattice%nrec
         this%charge%trq = 0.0d0
         this%charge%cht = 0.0d0
         this%charge%cht = sum(this%qia_new(ia, 1:6)) - this%symbolic_atom(this%lattice%nbulk + ia)%element%valence
         this%charge%dq(ia) = sum(this%qia(ia, 1:6)) - this%symbolic_atom(this%lattice%nbulk + ia)%element%valence
         this%charge%trq = this%charge%dq(ia)
         if (rank == 0) call g_logger%info('Valence at atom'//fmt('i4', (ia))//' '//'is'//' '//fmt('f10.6', (this%symbolic_atom(this%lattice%nbulk + ia)%element%valence)), __FILE__, __LINE__)
         if (rank == 0) call g_logger%info('Charge transfer at atom '//fmt('i4', (ia))//': '// &
                                           fmt('f10.6', this%charge%cht)//' '//fmt('f10.6', this%charge%trq)//' '//fmt('f10.6', this%charge%cht - this%charge%trq) &
                                           , __FILE__, __LINE__)
         !if(rank==0) call g_logger%info('Charge transfer at atom '//fmt('i4',(ia)), __FILE__, __LINE__)
         !if(rank==0) call g_logger%info(fmt('f10.6',this%charge%cht)//' '//fmt('f10.6',this%charge%trq)//' '//fmt('f10.6',this%charge%cht-this%charge%trq), __FILE__, __LINE__)
      end do

      do ia = 1, this%lattice%nrec
         delta_atom = sqrt(sum((this%qia_old(ia, 1:12) - this%qia_new(ia, 1:12))**2))/6.0d0
         if (rank == 0) call g_logger%info('Moment diff for atom:'//fmt('i4', (ia))//' '//'is'//' '//fmt('f10.6', delta_atom), __FILE__, __LINE__)
      end do
      this%delta = sqrt(sum((this%qia_old(:, 1:12) - this%qia_new(:, 1:12))**2))/6.0d0/this%lattice%nrec
   end subroutine mixpq
   subroutine broydn(pmix, amix, reset, mu, f, fsq, imu, itr, fsqo, u, v, muo, fo, nmix)
      ! mix potentials with broydens jacobian updating method as
      ! implemented by g.p.srivastava, j. phys a 17, l317(1984),
      ! changed 021015 by Lars Bergqvist in order to avoid any writing on files
      ! eqs. 12 through 16, and yesim n. darici 8/85.
      ! input  : mu =  input potential vector ( destroyed )
      !          f  = output potential vector ( destroyed )
      ! output : mu =  new input vector
      ! ierr :  0 normal processing
      !         1 input reset
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      logical, intent(inout) :: reset
      integer, intent(in) :: imu
      real(rp), intent(in) :: pmix
      real(rp), intent(in) :: amix
      real(rp), intent(out) :: fsq
      real(rp), dimension(imu), intent(inout) :: f, mu, u, v, muo, fo
      integer, intent(inout) :: itr
      real(rp), intent(inout) :: fsqo
      integer, intent(in) :: nmix
      !
      !.. Local Scalars ..
      logical :: wjac
      integer :: i, ierr, it, n, itrn, diverg, myid
      real(rp) :: df2, t
      !
      !.. Local Arrays ..
      real(rp), dimension(imu) :: &
         df, dmu, w1, w2
      !
      !.. External Functions ..
      real(rp), external :: ddot
      !
      ! ... Executable Statements ...
      !
      myid = 0
      !
      wjac = .true.
      diverg = 0
      ierr = 0
      !  n = 0
      n = imu
      !
      !  do i = 1, nvs
      !     n = n + nves(i)
      !  end do
      do i = 1, n
         f(i) = f(i) - mu(i)
      end do
      fsq = ddot(n, f, 1, f, 1)/n
      !print *, 'BROYDN:', fsq, fsqo, fsq>fsqo
      !print *, 'BROYDN:', itr, nmix, reset
      if (.not. reset) then
         if (itr == 0 .or. fsq > fsqo) then
            ! if (itr==0 .or. itr>nmix .or. fsq>fsqo) then
            ! if (itr==0 .or. itr>nmix) then
            ! if (itr==0 .or. itr>4 .or. fsq>fsqo) then
            ! if (itr==0 .or. itr>10 .or. fsq>fsqo) then
            reset = .true.
            if (fsq > fsqo) then
               diverg = 1
            end if
         end if
      end if
      if (reset) then
         itr = 0
         ierr = 1
      end if
      !print *, 'Bmix mid', itr, reset, diverg
      ! do i = 1, nvs
      !     pmix(i) = amix
      !     rmix(i) = pmix(i)
      ! end do
      if (itr /= 0) then
         ! read mu and f from previous iteration
         dmu = muo
         df = fo
      end if
    !!! if (myid==0) then
    !!!    write (6, "(/ a / a / 1x, i6, 1pe12.4, 2i12 /)") &
    !!!         " Broydn:", "    itr         fsq      length        ierr", itr, fsq, &
    !!!         n, ierr
    !!!    write (0, "(/ a / a / 1x, i6, 1pe12.4, 2i12 /)") &
    !!!         " Broydn:", "    itr         fsq      length        ierr", itr, fsq, &
    !!!         n, ierr
    !!! end if
      !write (876, "(/ a / a / 1x, i6, 1pe12.4, 2i12 /)") &
      !      " Broydn:", "    itr         fsq      length        ierr", itr, fsq, &
      !      n, ierr
      if (wjac) then
         !
         itrn = itr + 1
         muo = mu
         fo = f
      end if
      ! set mix vector
      ! k = 0
      ! do i = 1, nvs
      !     do j = 1, nves(i)
      !       k = k + 1
      !       mixv(k) = rmix(i)
      !     end do
      ! end do
      !
      !print *, 'BROYDEN', itr, itrn, n, nmix
      if (itr == 0) then
         do i = 1, n
            ! mu(i) = mu(i) + mixv(i)*f(i)
            mu(i) = mu(i) + pmix*f(i)
            ! mu(i) = mu(i) + 0.001d0*f(i)
         end do
      elseif (itr == 1) then
         do i = 1, n
            ! u(i) = mu(i) - dmu(i) + mixv(i)*(f(i)-df(i))
            u(i) = mu(i) - dmu(i) + amix*(f(i) - df(i))
            v(i) = f(i) - df(i)
         end do
         !
         df2 = ddot(n, v, 1, v, 1)
         do i = 1, n
            v(i) = v(i)/df2
         end do
         !
         t = ddot(n, v, 1, f, 1)
         do i = 1, n
            ! mu(i) = mu(i) + mixv(i)*f(i) - u(i)*t
            mu(i) = mu(i) + amix*f(i) - u(i)*t
         end do
      else
         ! itr > 1
         do i = 1, n
            dmu(i) = mu(i) - dmu(i)
            df(i) = f(i) - df(i)
            w1(i) = 0
            w2(i) = 0
         end do
         !
         do it = 1, itr - 1
            !
            t = ddot(n, v, 1, f, 1)
            do i = 1, n
               w1(i) = w1(i) + u(i)*t
            end do
            !
            t = ddot(n, v, 1, df, 1)
            do i = 1, n
               w2(i) = w2(i) + u(i)*t
            end do
            !
         end do
         !
         do i = 1, n
            ! u(i) = dmu(i) + mixv(i)*df(i) - w2(i)
            u(i) = dmu(i) + amix*df(i) - w2(i)
            v(i) = df(i)
         end do
         !
         df2 = ddot(n, v, 1, v, 1)
         do i = 1, n
            v(i) = v(i)/df2
         end do
         !
         t = ddot(n, v, 1, f, 1)
         do i = 1, n
            w1(i) = w1(i) + u(i)*t
         end do
         !
         do i = 1, n
            ! mu(i) = mu(i) + mixv(i)*f(i) - w1(i)
            !write(877, '(i4, 3f12.6)') i, mu(i), amix*f(i), w1(i)
            mu(i) = mu(i) + amix*f(i) - w1(i)
         end do
      end if
      itr = itrn
      fsqo = fsq
      if (itr > nmix) itr = 1
      !print *, 'Bmix done', itr, fsq
   end subroutine broydn
end module mix_mod

