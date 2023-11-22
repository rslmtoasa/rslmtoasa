!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Spin Dynamics
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
!> Module to handle the spin dynamics related processes
!------------------------------------------------------------------------------


module spin_dynamics_mod

  use mpi_mod
  use Depondt
  use RandomNumbers
  use control_mod
  use self_mod
  use energy_mod
  use lattice_mod
  use charge_mod
  use symbolic_atom_mod
  use hamiltonian_mod
  use recursion_mod
  use green_mod
  use density_of_states_mod
  use bands_mod
  use exchange_mod
  use mix_mod
  use math_mod
  use precision_mod
  use string_mod
  use self_mod
  use timer_mod, only: g_timer
  use logger_mod, only: g_logger
  implicit none

  private

  type, public :: spin_dynamics
    !> Control
    class(control), pointer :: control
    !> Lattice
    class(lattice), pointer :: lattice
    !> Self
    class(self), pointer :: self
    !> Symbolic atom
    class(symbolic_atom), dimension(:), pointer :: symbolic_atom
    !> Mix
    class(mix), pointer :: mix
    !> Bands
    class(bands), pointer :: bands

    !> Integrator
    character(len=5) :: integrator

    !> Number of time points
    integer :: nt

    !> SD step number
    integer :: asd_step

    !> Initial time in fs
    real(rp) :: t_i
   
    !> Final time in fs
    real(rp) :: t_f

    !> Time step in fs
    real(rp) :: dt  

    !> Gilbert damping
    real(rp) :: alpha
  
    !> SD temperature
    real(rp) :: sd_temp

    !> Stochastic magnetic field
    real(rp), dimension(:,:), allocatable :: b_stochastic

    !> Effective magnetic field
    real(rp), dimension(:,:), allocatable :: beff, b2eff
 
    !> Torque
    real(rp), dimension(:,:), allocatable :: btorque, she_btorque

    !> Magnetic moment variables
    real(rp), dimension(:,:), allocatable :: emom, emom2, emomM

    !> Manetic moment norm
    real(rp), dimension(:), allocatable :: mmom
  contains
    ! Destructror
    final :: destructor
    ! Procedures
    procedure :: build_from_file
    procedure :: restore_to_default
    procedure :: asd_pred
    procedure :: asd_corr
    procedure :: asd_pred_euler
    procedure :: asd_corr_euler
    procedure :: sd_run
  end type spin_dynamics

  interface spin_dynamics
    procedure :: constructor
  end interface spin_dynamics

contains
  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Constructor
  !
  !> @param[in] self_obj Pointer to system's self object
  !> @return type(spin_dynamics)
  !---------------------------------------------------------------------------
  function constructor(self_obj) result(obj)
    type(spin_dynamics) :: obj
    class(self), target, intent(in) :: self_obj

    obj%self          => self_obj
    obj%control       => self_obj%control
    obj%lattice       => self_obj%lattice
    obj%symbolic_atom => self_obj%symbolic_atom
    obj%mix           => self_obj%mix
    obj%bands         => self_obj%bands
 
    call obj%restore_to_default()
    call obj%build_from_file()
  end function 

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine destructor(this)
    type(spin_dynamics) :: this
#ifdef USE_SAFE_ALLOC
    if(allocated(this%b_stochastic)) call g_safe_alloc%deallocate('spin_dynamics.b_stochastic',this%b_stochastic)
    if(allocated(this%beff))         call g_safe_alloc%deallocate('spin_dynamics.beff',this%beff)
    if(allocated(this%b2eff))        call g_safe_alloc%deallocate('spin_dynamics.b2eff',this%b2eff)
    if(allocated(this%btorque))      call g_safe_alloc%deallocate('spin_dynamics.btorque',this%btorque)
    if(allocated(this%she_btorque))  call g_safe_alloc%deallocate('spin_dynamics.she_btorque',this%she_btorque)
    if(allocated(this%emom))         call g_safe_alloc%deallocate('spin_dynamics.emom',this%emom)
    if(allocated(this%emom2))        call g_safe_alloc%deallocate('spin_dynamics.emom2',this%emom2)
    if(allocated(this%emomM))        call g_safe_alloc%deallocate('spin_dynamics.emomM',this%emomM)
    if(allocated(this%mmom))         call g_safe_alloc%deallocate('spin_dynamics.mmom',this%mmom)
#else
    if(allocated(this%b_stochastic)) deallocate(this%b_stochastic)
    if(allocated(this%beff))         deallocate(this%beff)
    if(allocated(this%b2eff))        deallocate(this%b2eff)
    if(allocated(this%btorque))      deallocate(this%btorque)
    if(allocated(this%she_btorque))  deallocate(this%she_btorque)
    if(allocated(this%emom))         deallocate(this%emom)
    if(allocated(this%emom2))        deallocate(this%emom2)
    if(allocated(this%emomM))        deallocate(this%emomM)
    if(allocated(this%mmom))         deallocate(this%mmom)
#endif
  end subroutine destructor

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Read parameters from input file
  !---------------------------------------------------------------------------
  subroutine build_from_file(this)
    class(spin_dynamics), intent(inout) :: this

    ! variables associated with the reading processes
    integer :: iostatus, funit, i

    include 'include_codes/namelists/spin_dynamics.f90'

    integrator = this%integrator 
    nt         = this%nt
    t_i        = this%t_i
    t_f        = this%t_f
    dt         = this%dt
    alpha      = this%alpha
    asd_step   = this%asd_step
    sd_temp    = this%sd_temp

    call move_alloc(this%b_stochastic,b_stochastic)

    ! Reading
    open(newunit=funit, file=this%control%fname, action='read', iostat=iostatus, status='old')
    if(iostatus /= 0) then
      call g_logger%fatal('file '//trim(this%control%fname)//' not found', __FILE__, __LINE__)
    endif

    read(funit, nml=sd, iostat=iostatus)
    if(iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
      call g_logger%error('Error while reading namelist', __FILE__, __LINE__)
      call g_logger%error('iostatus = '//int2str(iostatus), __FILE__, __LINE__)
    endif
    close(funit)

    this%integrator = integrator 
    this%nt         = nt
    this%t_i        = t_i
    this%t_f        = t_f
    this%dt         = dt
    this%alpha      = alpha
    this%asd_step   = asd_step
    this%sd_temp    = sd_temp

    call move_alloc(b_stochastic,this%b_stochastic)
  end subroutine build_from_file

  subroutine restore_to_default(this)
    class(spin_dynamics), intent(inout) :: this

#ifdef USE_SAFE_ALLOC
    call g_safe_alloc%allocate('spin_dynamics.b_stochastic', this%b_stochastic, (/3,this%lattice%nrec/))
    call g_safe_alloc%allocate('spin_dynamics.beff'        , this%beff        , (/3,this%lattice%nrec/))
    call g_safe_alloc%allocate('spin_dynamics.b2eff'       , this%b2eff       , (/3,this%lattice%nrec/))
    call g_safe_alloc%allocate('spin_dynamics.btorque'     , this%btorque     , (/3,this%lattice%nrec/))
    call g_safe_alloc%allocate('spin_dynamics.she_btorque' , this%she_btorque , (/3,this%lattice%nrec/))
    call g_safe_alloc%allocate('spin_dynamics.emom'        , this%emom        , (/3,this%lattice%nrec/))
    call g_safe_alloc%allocate('spin_dynamics.emom2'       , this%emom2       , (/3,this%lattice%nrec/))
    call g_safe_alloc%allocate('spin_dynamics.emomM'       , this%emomM       , (/3,this%lattice%nrec/))
    call g_safe_alloc%allocate('spin_dynamics.mmom'        , this%emmom       ,   (/this%lattice%nrec/))
#else
    allocate(this%b_stochastic(3, this%lattice%nrec))
    allocate(this%beff(        3, this%lattice%nrec))
    allocate(this%b2eff(       3, this%lattice%nrec))
    allocate(this%btorque(     3, this%lattice%nrec))
    allocate(this%she_btorque( 3, this%lattice%nrec))
    allocate(this%emom(        3, this%lattice%nrec))
    allocate(this%emom2(       3, this%lattice%nrec))
    allocate(this%emomM(       3, this%lattice%nrec))
    allocate(this%mmom(           this%lattice%nrec))
#endif

    this%integrator   = 'none'
    this%t_i          = 0.0d0
    this%t_f          = 0.0d0
    this%dt           = 1.0d-16
    this%alpha        = 0.0d0
    this%b_stochastic = 0.0d0
    this%beff         = 0.0d0
    this%b2eff        = 0.0d0
    this%btorque      = 0.0d0
    this%she_btorque  = 0.0d0
    this%emom         = 0.0d0
    this%emom2        = 0.0d0
    this%emomM        = 0.0d0
    this%mmom         = 0.0d0
    this%sd_temp      = 0.0d0
    this%asd_step     = 0
  end subroutine restore_to_default

  ! Member subroutines and functions
  !---------------------------------------------------------------------------
  ! DESCRIPTION: 
  !> @brief
  !> Performs the prediction step of the spin dynamics calculation
  !>      
  !> Performs the prediction step of the spin dynamics calculation
  !---------------------------------------------------------------------------
  subroutine asd_pred(this,mom_in,field_in,temp_in,dt_in,na_in)
    implicit none
    class(spin_dynamics), intent(inout)      :: this
    real(rp), dimension(3,na_in), intent(in) :: mom_in
    real(rp), dimension(3,na_in), intent(in) :: field_in
    real(rp), intent(in)                     :: dt_in
    real(rp), intent(in)                     :: temp_in
    integer, intent(in)                      :: na_in
    ! Local variables
    integer                                  :: i, j, natom, nplusbulk
    real(rp)                                 :: mnorm, dummy 
    real(rp), dimension(3)                   :: thermal_field

    dummy = 1.0d0
    this%emom2 = 0.0d0
    thermal_field = 0.0d0
    this%b2eff = 0.0d0
    
    do i=1, na_in
      mnorm = norm2(mom_in(:,i))
      this%mmom(i) = mnorm
      this%emom(:,i) = mom_in(:,i)/mnorm              
      this%emomM(:,i) = mom_in(:,i)
      this%beff(:,i) = - field_in(:,i) ! Check units
    end do

    call   depondt_evolve_first(na_in,1, [this%alpha], this%beff, this%b2eff, [0.0d0,0.0d0,0.0d0], this%emom, this%emom2, this%emomM, this%mmom, &
                                this%dt, [this%sd_temp],dummy,'N',thermal_field,'N',[0.0d0,0.0d0,0.0d0])
    do i=1,na_in
      nplusbulk = i + this%lattice%nbulk
      !if(.not.this%mix%is_induced(i)) then
        mnorm = norm2(mom_in(:,i))
        this%symbolic_atom(nplusbulk)%potential%mom0(:) = this%emom(:,i)*mnorm
      !else
        !this%emom(:,i)= this%symbolic_atom(nplusbulk)%potential%mom0(:)/mnorm
      !end if
    end do
  end subroutine asd_pred

  !---------------------------------------------------------------------------
  ! DESCRIPTION: 
  !> @brief
  !> Performs the correction step of the spin dynamics calculation
  !>      
  !> Performs the correction step of the spin dynamics calculation
  !---------------------------------------------------------------------------
  subroutine asd_corr(this,mom_in,field_in,temp_in,dt_in,na_in)
    implicit none
    class(spin_dynamics), intent(inout)      :: this
    real(rp), dimension(3,na_in), intent(in) :: mom_in
    real(rp), dimension(3,na_in), intent(in) :: field_in
    real(rp), intent(in)                     :: dt_in
    real(rp), intent(in)                     :: temp_in
    integer, intent(in)                      :: na_in
    ! Local variables
    integer                                  :: i, j, natom, nplusbulk
    real(rp)                                 :: mnorm

    do i=1, na_in
      mnorm = norm2(mom_in(:,i))
      this%mmom(i) = mnorm
      this%emom(:,i) = mom_in(:,i)/mnorm
      this%emomM(:,i) = mom_in(:,i)
      this%beff(:,i) = - field_in(:,i) ! Check units
    end do
    call   depondt_evolve_second(na_in,1, [this%alpha], this%beff, this%b2eff, [0.0d0,0.0d0,0.0d0], this%emom, this%emom2, this%dt, 'N','N',[0.0d0,0.0d0,0.0d0])
    do i=1,na_in
      nplusbulk = i + this%lattice%nbulk
      mnorm = norm2(mom_in(:,i))
      this%mmom(i) = mnorm
      !if(.not.this%mix%is_induced(i)) then
        this%symbolic_atom(nplusbulk)%potential%mom0(:) = this%emom2(:,i)*mnorm
      !else
      !  this%emom2(:,i)= this%symbolic_atom(nplusbulk)%potential%mom0(:)/mnorm
      !end if
    end do
  end subroutine asd_corr

  subroutine asd_pred_euler(this,mom_in, field_in, temp_in, dt_in, na_in)
      implicit none
      class(spin_dynamics), intent(inout)      :: this
      real(rp), dimension(3,na_in), intent(in) :: mom_in
      real(rp), dimension(3,na_in), intent(in) :: field_in
      real(rp), intent(in)                     :: dt_in
      real(rp), intent(in)                     :: temp_in
      integer, intent(in)                      :: na_in
      !real(rp), dimension(3,na_in), intent(out) :: mom_pred
      integer :: i, nplusbulk
      real(rp) :: alpha
      real(rp), dimension(3) :: torque_term1, torque_term2, dm
  
      alpha = this%alpha 

      this%b2eff(:,:) = 0.0d0
       
      do i = 1, na_in
          nplusbulk = i + this%lattice%nbulk
          this%symbolic_atom(nplusbulk)%potential%mom0(:) = 0.0d0
          torque_term1 = -gama * cross_product(mom_in(:, i), field_in(:, i))
          torque_term2 = -alpha * gama * cross_product(mom_in(:, i), cross_product(mom_in(:, i), field_in(:, i)))
          dm = torque_term1 + torque_term2
          this%b2eff(:,i) = dm
          this%symbolic_atom(nplusbulk)%potential%mom0(:) = mom_in(:, i) + dt_in * dm
          this%emom(:,i) = this%symbolic_atom(nplusbulk)%potential%mom0(:)/norm2(this%symbolic_atom(nplusbulk)%potential%mom0(:))
      end do
  end subroutine asd_pred_euler

  subroutine asd_corr_euler(this,mom_pred, field_in, temp_in, dt_in, na_in)
      implicit none
      class(spin_dynamics), intent(inout)      :: this
      real(rp), dimension(3,na_in), intent(in) :: mom_pred
      real(rp), dimension(3,na_in), intent(in) :: field_in
      real(rp), intent(in)                     :: dt_in
      real(rp), intent(in)                     :: temp_in
      integer, intent(in)                      :: na_in
      !real(rp), dimension(3,na_in), intent(out) :: mom_corr
  
      integer :: i, nplusbulk
      real(rp) :: alpha
      real(rp), dimension(3) :: torque_term1, torque_term2, dm
  
      this%alpha = alpha
   
      do i = 1, na_in
          nplusbulk = i + this%lattice%nbulk
          this%symbolic_atom(nplusbulk)%potential%mom0(:) = 0.0d0
          ! Recalculating the torque terms using the predicted moments
          torque_term1 = -gama * cross_product(mom_pred(:, i), field_in(:, i))
          torque_term2 = -alpha * gama * cross_product(mom_pred(:, i), cross_product(mom_pred(:, i), field_in(:, i)))
          dm = 0.5d0*(this%b2eff(:,i) + (torque_term1 + torque_term2))
          ! Updating the magnetic moments based on the corrected torque
          this%symbolic_atom(nplusbulk)%potential%mom0(:) = mom_pred(:, i) + dt_in * dm
          this%emom2(:,i) = this%symbolic_atom(nplusbulk)%potential%mom0(:)/norm2(this%symbolic_atom(nplusbulk)%potential%mom0(:))
      end do
  end subroutine asd_corr_euler

  subroutine sd_run(this)
    implicit none
    class(spin_dynamics), intent(inout)      :: this
    ! Local variables
    integer                                  :: n_step
    real(rp), dimension(3,this%lattice%nrec) :: mom_in, mom_prev
    real(rp), dimension(3,this%lattice%nrec) :: field_in
    real(rp)                                 :: dt_in, timestep
    real(rp)                                 :: temp_in
    integer                                  :: na_in, na

    dt_in   = this%dt
    na_in   = this%lattice%nrec
    temp_in = this%sd_temp

    call allocate_depondtfields(1,1,1)
    call setup_rng_hb(1, 'Y','N')
    call allocate_randomwork(1,1,1,'N')
    call mt_ran_init_c(1) ! Monte-Carlo

    call this%self%run()
    do na=1,na_in
       mom_prev(:,na) = this%symbolic_atom(na+this%lattice%nbulk)%potential%mom0(:)
    end do

    timestep = 0.0d0

    do n_step=1, this%asd_step

      timestep = timestep + this%dt

      if(rank==0)call g_logger%info('Spin dynamics at step '//int2str(n_step), __FILE__, __LINE__)

      ! Calculate the effective field for the predictor step
      call this%self%run()
      call this%bands%calculate_magnetic_torques()
      do na=1,na_in
        !mom_in(:,na) = this%symbolic_atom(na+this%lattice%nbulk)%potential%mom0(:)
        mom_in(:,na) = mom_prev(:,na)
      end do
      field_in(:,:) = 0.0d0
      field_in(:,:) = -this%bands%mag_for(:,:)
      ! Do the predictor step
      call this%asd_pred_euler(mom_in,field_in,temp_in,dt_in,na_in)
      do na=1,na_in
        this%symbolic_atom(na+this%lattice%nbulk)%potential%mom(:) = this%emom(:,na)
        mom_prev(:,na) = this%symbolic_atom(na+this%lattice%nbulk)%potential%mom0(:)
      end do


      ! Update the effective field 
      !call this%self%run()
      !call this%bands%calculate_magnetic_torques()
      !do na=1,na_in
        !mom_in(:,na) = this%symbolic_atom(na+this%lattice%nbulk)%potential%mom0(:)
      !  mom_in(:,na) = mom_prev(:,na)
      !end do
      !field_in(:,:) = 0.0d0
      !field_in(:,:) = this%bands%mag_for(:,:)

      ! Do the corrector step
      !write(*,*) 'corr', this%emom
      !call this%asd_corr_euler(mom_in,field_in,temp_in,dt_in,na_in)
      !write(*,*) 'corr', this%emom2
      !do na=1,na_in
      !  this%symbolic_atom(na+this%lattice%nbulk)%potential%mom(:) = this%emom2(:,na)
      !  mom_prev(:,na) = this%symbolic_atom(na+this%lattice%nbulk)%potential%mom0(:)
      !end do
      write(32,*) timestep*1.d+15, this%emom 
    end do
    call allocate_depondtfields(1,1,-1)

  end subroutine sd_run
end module spin_dynamics_mod
