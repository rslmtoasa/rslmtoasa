!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Self
!
!> @author
!> Angela Klautau
!> Ramon Cardias
!> Lucas f. Campagna
!> S. Frota-Pessôa
!> Pascoal R. Peduto
!> Anders Bergman
!> S. B. Legoas
!> H. M. Petrilli
!> Ivan P. Miranda
!
! DESCRIPTION:
!> Module to handle self-consistent process
!------------------------------------------------------------------------------

module self_mod

   use mpi_mod
   use symbolic_atom_mod, only: symbolic_atom
   use logger_mod, only: g_logger
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   use string_mod, only: real2str, int2str, fmt
   use control_mod
   use lattice_mod
   use charge_mod
   use xc_mod
   use recursion_mod
   use density_of_states_mod
   use green_mod
   use bands_mod
   use energy_mod
   use hamiltonian_mod
   use mix_mod
   use math_mod
   use precision_mod, only: rp
   use timer_mod, only: g_timer
   use namelist_generator_mod, only: namelist_generator
#ifdef USE_MPI
   use mpi
#endif
   implicit none

   private

   !> Module´s main structure
   type, public :: self
      !> Lattice
      class(lattice), pointer :: lattice
      !> Charge
      class(charge), pointer :: charge
      !> Control
      class(control), pointer :: control
      !> Symbolic atom
      class(symbolic_atom), dimension(:), pointer :: symbolic_atom
      !> Recursion
      class(recursion), pointer :: recursion
      !> Density of states
      class(dos), pointer :: dos
      !> Bands
      class(bands), pointer :: bands
      !> Energy
      class(energy), pointer :: en
      !> Mix
      class(mix), pointer :: mix
      !> Hamiltonian
      class(hamiltonian), pointer :: hamiltonian
      !> Green
      class(green), pointer :: green

      !TODO: check description
      !> If true treats all atoms as inequivalents. Default: true.
      !>
      !> If true treats all atoms as inequivalents.
      !>
      !> Default: true for all kind of calculations. If false, the user may provide the final lines of old format self file in a different file.
      logical :: all_inequivalent

      ! Calculation type dependents
      !TODO: lack description
      !> Default: 8
      !>
      !> Default: 8
      integer :: init

      !> Use same value of @ref ws for all atoms. Default: true.
      !>
      !> Use same value of @ref ws for all atoms. If true, @ref ws´s size is one slot of memory, else it is @ref lattice.nrec
      !>
      !> Default: true.
      logical :: ws_all

      !> Wigner Seitz Radius
      !>
      !> Wigner Seitz Radius. If @ref ws_all is true, @ref ws´s size is one slot of memory, else it is @ref lattice.nrec
      real(rp), dimension(:), allocatable :: ws

      ! Mixing parameters

      !> Use same value of @ref mix for all atoms. Default: true.
      !>
      !> Use same value of @ref mix for all atoms. If true, @ref mix´s size is one slot of memory, else it is @ref lattice.nrec
      !>
      !> Default: true.
      logical :: mix_all

      ! Magnetic mixing parameters

      !> If true enable spin-(up/down) mixing. Default: false.
      !>
      !> If true enable spin-(up/down) mixing, else disable.
      !>
      !> Default: false.
      logical :: magnetic_mixing

      !> Use same value of @ref mixmag for all atoms. Default: true.
      !>
      !> Use same value of @ref mixmag for all atoms. If true, @ref mixmag´s size is one slot of memory, else it is @ref lattice.nrec
      !>
      !> Default: true.
      logical :: mixmag_all

      !> Spin-(up/down) occupation in self-consistent calculation. Default: 0.05.
      !>
      !> Spin-(up/down) occupation in self-consistent calculation.
      !>
      !> Default: 0.05.
      real(rp), dimension(:), allocatable :: mixmag

      ! Convergence parameters

      !> Number of steps in calculation. Default: 1.
      !>
      !> Number of steps in self-consistent calculation.
      !>
      !> Default: 1.
      integer :: nstep

      !> Convergency precision threshold. Default: 0.5d-10.
      !>
      !> Specify the precision required in charge density to stop the self-consistent process.
      !>
      !> Default: 0.5d-10.
      real(rp) :: conv_thr

      ! Constrains variables

      !> Freezes magnetic moment direction during nonlinear calculation. Default: false.
      !>
      !> Freezes magnetic moment direction during nonlinear calculation. Do not set true in collinear calculation.
      !>
      !> Default: false.
      logical :: freeze

      !> If true rigid band calculations will performed for each atom. Default: false.
      !>
      !> If true rigid band calculations will performed for each atom. Using this option you should specify the variable @ref rb.
      !>
      !> Default: false.
      logical :: rigid_band

      !> Specify the number of rigid band calculation used for each atom.
      !>
      !> Specify the number of rigid band calculation used for each atom. The number of values provided should be the same as @ref lattice.nrec.
      integer, dimension(:), allocatable :: rb

      !> If true, performs a calculation with moment orbital polarization. Default: false.
      !>
      !> If true, performs a calculation with moment orbital polarization.
      !>
      !> Default: false.
      logical :: orbital_polarization

      !> Sum of eigen
      real(rp) :: totsumev
      !> Band energy
      real(rp) :: ebsum
      !> Total energy
      real(rp) :: esumn, esum
      !> Total energy w.o Madelung
      real(rp) :: esumc

      !> TODO: from commons.commons_caro
      integer, dimension(:), allocatable :: ifc
      integer :: ifcc

      ! TODO
      real(rp), dimension(:, :), allocatable :: vtn, vzt
      real(rp), dimension(:, :, :), allocatable :: fun2

      ! TODO
      real(rp), dimension(:), allocatable :: bxc

      ! TODO max ws
      ! User in lmtst: MAKE POTENTIAL PARAMETERS TO SAME DNU VALUES
      real(rp) :: ws_max
      !> Logical variable to check if the calculation is converged.
      logical :: converged

   contains
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: print_state_formatted
      procedure :: print_state
      procedure :: print_state_full
      procedure :: run
      procedure :: report
      procedure :: lmtst
      procedure :: is_converged
      procedure, private :: atomsc
      procedure, private :: ftype
      procedure, private :: newrho
      procedure, private :: rhocor
      procedure, private :: vxc0sp
      procedure, private :: racsi
      procedure, private :: potpar
      final :: destructor
   end type self

   interface self
      procedure :: constructor
   end interface self

contains
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] fname Namelist file
   !> @param[in] lattice_obj Pointer to system´s lattice
   !> @return type(self)
   !---------------------------------------------------------------------------
   function constructor(bands_obj, mix_obj) result(obj)
      type(self) :: obj
      type(bands), target, intent(in) :: bands_obj
      type(mix), target, intent(in) :: mix_obj

      obj%lattice => bands_obj%lattice
      obj%charge => bands_obj%recursion%hamiltonian%charge
      obj%control => bands_obj%lattice%control
      obj%symbolic_atom => bands_obj%symbolic_atom
      obj%recursion => bands_obj%recursion
      obj%dos => bands_obj%dos
      obj%bands => bands_obj
      obj%en => bands_obj%en
      obj%mix => mix_obj
      obj%hamiltonian => bands_obj%recursion%hamiltonian
      obj%green => bands_obj%green

      call obj%restore_to_default()
      call obj%build_from_file()
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   subroutine destructor(this)
      type(self) :: this
#ifdef USE_SAFE_ALLOC
      if (allocated(this%ws)) call g_safe_alloc%deallocate('self.ws', this%ws)
      if (allocated(this%mixmag)) call g_safe_alloc%deallocate('self.mixmag', this%mixmag)
      if (allocated(this%rb)) call g_safe_alloc%deallocate('self.rb', this%rb)
      if (allocated(this%bxc)) call g_safe_alloc%deallocate('self.bxc', this%bxc)
      if (allocated(this%vtn)) call g_safe_alloc%deallocate('self.vtn', this%vtn)
      if (allocated(this%vzt)) call g_safe_alloc%deallocate('self.vzt', this%vzt)
      if (allocated(this%fun2)) call g_safe_alloc%deallocate('self.fun2', this%fun2)
#else
      if (allocated(this%ws)) deallocate (this%ws)
      if (allocated(this%mixmag)) deallocate (this%mixmag)
      if (allocated(this%rb)) deallocate (this%rb)
#endif
   end subroutine destructor

   ! Member functions

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Read parameters from input file
   !---------------------------------------------------------------------------
   subroutine build_from_file(this)
      class(self), intent(inout) :: this

      ! variables associated with the reading processes
      integer :: iostatus, funit, i

      include 'include_codes/namelists/self.f90'

      ! Save previous values
      all_inequivalent = this%all_inequivalent
      ws_all = this%ws_all
      mix_all = this%mix_all
      magnetic_mixing = this%magnetic_mixing
      mixmag_all = this%mixmag_all
      freeze = this%freeze
      rigid_band = this%rigid_band
      orbital_polarization = this%orbital_polarization
      init = this%init
      nstep = this%nstep
      conv_thr = this%conv_thr

      call move_alloc(this%ws, ws)
      call move_alloc(this%mixmag, mixmag)
      call move_alloc(this%rb, rb)

      ! Check if size is right
      if (size(mixmag) .ne. this%lattice%nrec) then
         call g_logger%error('resizing array "mixmag"', __FILE__, __LINE__)
         call g_logger%error('from '//int2str(size(mixmag))//' to '//int2str(this%lattice%nrec), __FILE__, __LINE__)
         call g_logger%error('If this variable is not in the namelist, it may cause errors', __FILE__, __LINE__)
         deallocate (mixmag)
         allocate (mixmag(this%lattice%nrec))
      end if
      if (size(ws) .ne. this%lattice%nrec) then
         call g_logger%error('resizing array "ws"', __FILE__, __LINE__)
         call g_logger%error('from '//int2str(size(ws))//' to '//int2str(this%lattice%nrec), __FILE__, __LINE__)
         call g_logger%error('If this variable is not in the namelist, it may cause errors', __FILE__, __LINE__)
         deallocate (ws)
         allocate (ws(this%lattice%nrec))
      end if
      if (size(rb) .ne. this%lattice%nrec) then
         call g_logger%error('resizing array "rb"', __FILE__, __LINE__)
         call g_logger%error('from '//int2str(size(rb))//' to '//int2str(this%lattice%nrec), __FILE__, __LINE__)
         call g_logger%error('If this variable is not in the namelist, it may cause errors', __FILE__, __LINE__)
         deallocate (rb)
         allocate (rb(this%lattice%nrec))
      end if

      ! Reading
      open (newunit=funit, file=this%lattice%control%fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//trim(this%lattice%control%fname)//' not found', __FILE__, __LINE__)
      end if

      read (funit, nml=self, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%error('Error while reading namelist', __FILE__, __LINE__)
         call g_logger%error('iostatus = '//int2str(iostatus), __FILE__, __LINE__)
      end if
      close (funit)

      ! Setting user values

      ! Control variables
      this%all_inequivalent = all_inequivalent

      ! Wigner Seitz Radius
      this%ws_all = ws_all
      if (ws_all) then
#ifdef USE_SAFE_ALLOC
         call g_safe_alloc%allocate('self.ws', this%ws, 1)
#else
         allocate (this%ws(1))
#endif
         this%ws(1) = ws(1)
      else
         call move_alloc(ws, this%ws)
      end if

      ! Mixing parameters
      this%mix_all = mix_all

      ! Magnetic mixing parameters
      this%magnetic_mixing = magnetic_mixing
      this%mixmag_all = mixmag_all
      if (magnetic_mixing) then
         if (mixmag_all) then
#ifdef USE_SAFE_ALLOC
            call g_safe_alloc%allocate('self.mixmag', this%mixmag, 1)
#else
            allocate (this%mixmag(1))
#endif
            this%mixmag(1) = mixmag(1)
         else
            call move_alloc(mixmag, this%mixmag)
         end if
      else
#ifdef USE_SAFE_ALLOC
         call g_safe_alloc%allocate('self.mixmag', this%mixmag, 0)
#else
         allocate (this%mixmag(0))
#endif
      end if

      ! Convergence parameters
      this%conv_thr = conv_thr
      this%nstep = nstep

      ! Constrains variables
      ! static magnetic momentum (non-linar calculation)
      this%freeze = freeze
      this%rigid_band = rigid_band

      ! number of loops to use rigid bands per atom
      call move_alloc(rb, this%rb)

      ! other variables
      this%orbital_polarization = orbital_polarization
      this%init = init
   end subroutine build_from_file

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this, full)
      class(self), intent(inout):: this
      logical, intent(in), optional :: full
      integer :: lmax, nsp

      ! Control variables
      ! if false force to read the original self file
      this%all_inequivalent = .true.

      ! Wigner Seitz Radius
      this%ws_all = .true.
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('self.ws', this%ws, this%lattice%nrec)
#else
      allocate (this%ws(this%lattice%nrec))
#endif
      this%ws(:) = this%lattice%wav*ang2au

      ! Mixing parameters
      this%mix_all = .true.

      ! Magnetic mixing parameters
      this%magnetic_mixing = .false.
      this%mixmag_all = .true.
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('self.mixmag', this%mixmag, this%lattice%nrec)
#else
      allocate (this%mixmag(this%lattice%nrec))
#endif
      this%mixmag(:) = 0.05

      ! Convergence parameters
      this%conv_thr = 0.5d-8
      this%nstep = 1

      ! Constrains variables
      ! static magnetic momentum (non-linar calculation)
      this%freeze = .false.
      this%rigid_band = .false.
      ! number of loops to use rigid bands per atom
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('self.rb', this%rb, this%lattice%nrec)
#else
      allocate (this%rb(this%lattice%nrec))
#endif
      this%rb(:) = 2

      this%orbital_polarization = .false.

      this%ws_max = 9.99d0

      if (associated(this%lattice)) then
         if (present(full)) then
            if (full) then
               call this%lattice%restore_to_default(full)
            end if
         end if
      end if

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('self.bxc', this%bxc, this%lattice%nrec)
      call g_safe_alloc%allocate('self.vtn', this%vtn, (/8001, 2/))
      call g_safe_alloc%allocate('self.vzt', this%vzt, (/8001, 2/))
      call g_safe_alloc%allocate('self.fun2', this%fun2, (/8001, 3, 2/))
#else
      allocate (this%bxc(this%lattice%nrec))
      allocate (this%vtn(8001, 2), this%vzt(8001, 2), this%fun2(8001, 3, 2))
#endif

   end subroutine restore_to_default

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print class members values formatted
   !>
   !> Print class members values formatted
   !---------------------------------------------------------------------------
   subroutine print_state_formatted(this)
      class(self), intent(in) :: this

      print *, 'Printing SELF object'
      print *, ''
      print *, '[Control Variables]'
      print *, 'all_inequivalent ', this%all_inequivalent
      print *, ''
      print *, '[Wigner Seitz Radius]'
      print *, 'ws_all ', this%ws_all
      print *, 'ws     ', this%ws
      print *, ''
      print *, '[Mixing Parameters]'
      print *, 'mix_all ', this%mix_all
      print *, ''
      print *, '[Magnetic Mixing Parameters]'
      print *, 'magnetic_mixing ', this%magnetic_mixing
      print *, 'mixmag_all      ', this%mixmag_all
      print *, 'mixmag          ', this%mixmag
      print *, ''
      print *, '[Convergence Parameters]'
      print *, 'conv_thr ', this%conv_thr
      print *, 'nstep    ', this%nstep
      print *, ''
      print *, '[Constrain Variables]'
      print *, 'freeze     ', this%freeze
      print *, 'rigid_band ', this%rigid_band
      print *, 'rb         ', this%rb
      print *, ''
      print *, '[Other Variables]'
      print *, 'orbital_polarization ', this%orbital_polarization
      print *, 'init                 ', this%init

   end subroutine print_state_formatted

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print class members values in namelist format
   !>
   !> Print class members values in namelist format
   !---------------------------------------------------------------------------
   subroutine print_state_full(this)
      class(self), intent(in) :: this

      include 'include_codes/namelists/self.f90'

      ! scalar

      conv_thr = this%conv_thr
      ws_all = this%ws_all
      rigid_band = this%rigid_band
      orbital_polarization = this%orbital_polarization
      mixmag_all = this%mixmag_all
      mix_all = this%mix_all
      magnetic_mixing = this%magnetic_mixing
      freeze = this%freeze
      all_inequivalent = this%all_inequivalent
      nstep = this%nstep
      init = this%init

      ! 1d allocatable

      if (allocated(this%ws)) then
         allocate (ws, mold=this%ws)
         ws = this%ws
      else
         allocate (ws(0))
      end if
      if (allocated(this%mixmag)) then
         allocate (mixmag, mold=this%mixmag)
         mixmag = this%mixmag
      else
         allocate (mixmag(0))
      end if
      if (allocated(this%rb)) then
         allocate (rb, mold=this%rb)
         rb = this%rb
      else
         allocate (rb(0))
      end if

      write (*, nml=self)

   end subroutine print_state_full

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print input class members values in namelist format
   !>
   !> Print input class members values in namelist format
   !---------------------------------------------------------------------------
   subroutine print_state(this)
      class(self), intent(in) :: this

      include 'include_codes/namelists/self.f90'

      ! scalar

      conv_thr = this%conv_thr
      ws_all = this%ws_all
      rigid_band = this%rigid_band
      orbital_polarization = this%orbital_polarization
      mixmag_all = this%mixmag_all
      mix_all = this%mix_all
      magnetic_mixing = this%magnetic_mixing
      freeze = this%freeze
      all_inequivalent = this%all_inequivalent
      nstep = this%nstep
      init = this%init
      ! 1d allocatable

      if (allocated(this%ws)) then
         allocate (ws, mold=this%ws)
         ws = this%ws
      else
         allocate (ws(0))
      end if
      if (allocated(this%mixmag)) then
         allocate (mixmag, mold=this%mixmag)
         mixmag = this%mixmag
      else
         allocate (mixmag(0))
      end if
      if (allocated(this%rb)) then
         allocate (rb, mold=this%rb)
         rb = this%rb
      else
         allocate (rb(0))
      end if

      write (*, nml=self)

   end subroutine print_state


   !> Row 425 in green.f90
   ! subroutine intersite_gf(this)
      ! use mpi_mod
      ! implicit none
      ! class(green), intent(inout) :: this
   ! end subroutine intersite_gf



   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Performs the self-consistent calculation.
   !>
   !> Performs the self-consistent calculation.
   !---------------------------------------------------------------------------
   subroutine run(this)
      class(self), intent(inout) :: this
      integer :: i, ia, niter, ie, l1, l2
      real(rp), dimension(6) :: QSL

      real(rp), dimension(:), allocatable :: pot_arr
      integer :: na_glob, pot_size
      real(rp), dimension(:, :), allocatable :: T_comm

      !===========================================================================
      !                              BEGIN SCF LOOP
      !===========================================================================
      niter = 0
      do i = 1, this%nstep
         !=========================================================================
         !                        PERFORM THE RECURSION
         !=========================================================================
         if (rank == 0) call g_logger%info('Perform recursion at step '//int2str(i), __FILE__, __LINE__)
         call g_timer%start('recursion')
         select case (this%control%calctype)
         case ('B')

            do ia = 1, this%lattice%nrec
               call this%symbolic_atom(ia)%build_pot() ! Build the potential matrix
            end do
            if (this%hamiltonian%hubbardU_check .and. this%hamiltonian%hubbardJ_check .and. i .gt. 1) then
               !> Initiate the LDA+U+J method by building the Hubbard U+J potential matrix
               call this%bands%spdf_Hubbard() ! Improved code
               ! call this%bands%build_hubbard_pot() ! Previous code (without +V)

               !> Initiate the +V intersite Coulomb correction to LDA+U+J
               if (this%recursion%hamiltonian%hubbardV_check) then
                  call this%recursion%recur_b_ij()
                  call this%green%calculate_intersite_gf()
                  call this%bands%Hubbard_V()
                  ! do ie = 1, this%en%channels_ldos + 10
                  !    print *, ''
                  !    print *, 'ie : ', ie
                     ! print *, 'gij : ', this%green%gij(:,:,ie,1)
                     ! print *, 'gji^T : ', transpose(this%green%gji(:,:,ie,1))
                  !    print *, ''
                  !    print *, '---------------------------------------------------------------------------'
                  ! end do
                  ! print *, ''
                  ! print *, 'gij : ', this%green%gij(1:2, 1:2 ,1000,1)
                  ! print *, ''
                  ! print *, 'gji^T : ', transpose(this%green%gji(1:2, 1:2 ,1000,1))
                  ! print *, ''
                  ! print *, 'gji : ', this%green%gji(1:2, 1:2 ,1000,1)
      
               end if
            end if
            if (this%control%nsp == 2 .or. this%control%nsp == 4) call this%hamiltonian%build_lsham ! Calculate the spin-orbit coupling Hamiltonian
            call this%hamiltonian%build_bulkham() ! Build the bulk Hamiltonian
         case ('S')
            do ia = 1, this%lattice%ntype
               call this%symbolic_atom(ia)%build_pot() ! Build the potential matrix
            end do
            if (this%control%nsp == 2 .or. this%control%nsp == 4) call this%hamiltonian%build_lsham ! Calculate the spin-orbit coupling Hamiltonian
            call this%hamiltonian%build_bulkham() ! Build the bulk Hamiltonian for the surface
         case ('I')
            do ia = 1, this%lattice%ntype
               call this%symbolic_atom(ia)%build_pot() ! Build the potential matrix
            end do
            if (this%control%nsp == 2 .or. this%control%nsp == 4) call this%hamiltonian%build_lsham ! Calculate the spin-orbit coupling Hamiltonian
            call this%hamiltonian%build_bulkham() ! Build the bulk Hamiltonian
            call this%hamiltonian%build_locham() ! Build the local Hamiltonian
         end select
         select case (this%control%recur)
         case ('lanczos')
            call this%recursion%recur()
         case ('chebyshev')
            call this%recursion%chebyshev_recur()
         case ('block')
            call this%recursion%recur_b()
         end select
         call g_timer%stop('recursion')
         !=========================================================================
         !               SAVE THE TOTAL ENERGY FROM PREVIOUS ITERATION
         !=========================================================================
         this%esumn = sum(this%symbolic_atom(:)%potential%etot)
         this%esum = this%esumn
         !=========================================================================
         !            SAVE THE PARAMETERS QL AND PL TO BE MIXED LATER
         !=========================================================================
         call this%mix%save_to('old') ! Save to qia_old to mix with qia_new.
         !=========================================================================
         !                      SAVE THE MAGNETIC MOMENTS
         !=========================================================================
         do ia = 1, this%lattice%nrec
            this%mix%mag_old(ia, :) = this%symbolic_atom(this%lattice%nbulk + ia)%potential%mom(:)
         end do
         !=========================================================================
         !                  CALCULATE THE DENSITY OF STATES
         !=========================================================================
         if (rank == 0) call g_logger%info('Calculating the density of states and the new moment bands', __FILE__, __LINE__)
         call g_timer%start('calculation-of-DOS')
         call this%en%e_mesh() ! Solve the energy mesh
         select case (this%control%recur)
         case ('lanczos')
            call this%green%sgreen() ! Calculate the density of states using the continued fraction
         case ('chebyshev')
            call this%green%chebyshev_green()
         case ('block')
            ! HERE
            call this%recursion%zsqr()
            call this%green%block_green()
         end select
         call this%bands%calculate_fermi() ! Calculate the fermi energy
         !=========================================================================
         !  MIX THE MAGNETIC MOMENTS BEFORE CALCULATING THE NEW BAND MOMENTS QL
         !=========================================================================
         call this%bands%calculate_magnetic_moments() ! Calculate the magnetic moments
         do ia = 1, this%lattice%nrec
            this%mix%mag_new(ia, :) = this%symbolic_atom(this%lattice%nbulk + ia)%potential%mom(:)
         end do
         call this%mix%mix_magnetic_moments(this%mix%mag_old, this%mix%mag_new, this%mix%mag_mix, this%symbolic_atom(:)%potential%mtot) ! Mix magnetic moments
         do ia = 1, this%lattice%nrec
            this%symbolic_atom(this%lattice%nbulk + ia)%potential%mom(:) = this%mix%mag_mix(ia, :)
         end do
         !=========================================================================
         !                  CALCULATE THE NEW BAND MOMENTS QL
         !=========================================================================
         call this%bands%calculate_moments() ! Integrate the DOS and calculate the band momends QL
         do ia = 1, this%lattice%nrec
            this%mix%mag_new(ia, :) = this%symbolic_atom(this%lattice%nbulk + ia)%potential%mom(:)
         end do
         call this%mix%save_to('new') ! Save the calculated PL and QL into the mix%qia_new
         call g_timer%stop('calculation-of-DOS')
         !=========================================================================
         !                   MIX OLD AND NEW CALCULATED PL AND QL
         !=========================================================================
         !call g_timer%start(´mixing´)
         if (rank == 0) call g_logger%info('Mixtype is '//trim(this%mix%mixtype), __FILE__, __LINE__)
         call this%mix%mixpq(this%mix%qia_old, this%mix%qia_new) ! Mix mix%qia_new with mix%qia_old and save into mix%qia
         !call g_timer%stop(´mixing´)
         !=========================================================================
         !         CALCULATE THE MADELUNG POTENTIAL (BULK ONLY IMPLEMENTED)
         !=========================================================================
         !call g_timer%start(´madelung-potential´)
         select case (this%control%calctype)
         case ('B')
            call this%charge%bulkpot()
         case ('S')
            call this%charge%surfpot()
         case ('I')
            call this%charge%imppot()
         end select
         !call g_timer%stop(´madelung-potential´)
         !=========================================================================
         !                        SAVE MIXED PARAMETERS
         !=========================================================================
         call this%mix%save_to('current') ! Save mixed parameters mix%qia into potential%pl and potential%ql
         !=========================================================================
         !                       MAKE SFC ATOMIC SPHERE
         !=========================================================================
         call g_timer%start('atomic-sfc')
         !do ia=1, this%lattice%nrec
         do ia = start_atom, end_atom
            qsl = this%lmtst(this%symbolic_atom(this%lattice%nbulk + ia)) ! Makes the atomic sphere self-consistent and caltulate the orthogonal pottential parameters
            call g_logger%info('Atomic SFC done for atom '//this%symbolic_atom(this%lattice%nbulk + ia)%element%symbol, __FILE__, __LINE__)
         end do

         !=========================================================================
         !      TRANSFER ATOMIC POTENTIAL DATA ACROSS MPI RANKS
         !=========================================================================
#ifdef USE_MPI
         pot_size = this%symbolic_atom(start_atom)%potential%sizeof_potential_full()
         allocate (T_comm(pot_size, this%lattice%nrec))
         T_comm = 0.0_rp
         do na_glob = start_atom, end_atom
            call this%symbolic_atom(this%lattice%nbulk + na_glob)%potential%flatten_potential_full(T_comm(:, na_glob))
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE, T_comm, product(shape(T_comm)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         do na_glob = 1, this%lattice%nrec
            call this%symbolic_atom(this%lattice%nbulk + na_glob)%potential%expand_potential_full(T_comm(:, na_glob))
         end do
         deallocate (T_comm)
#endif

         !=========================================================================
         !      TRANSFORM POTENTIAL BASIS FROM ORTHOGONAL TO TIGHT-BINDING
         !=========================================================================
         do ia = 1, this%lattice%nrec
            !do ia=1, start_atom, end_atom
            if (rank == 0) call g_logger%info('From orthogonal to TB basis for atom '//this%symbolic_atom(this%lattice%nbulk + ia)%element%symbol, __FILE__, __LINE__)
            call this%symbolic_atom(this%lattice%nbulk + ia)%predls(this%lattice%wav*ang2au) ! Transforms the potential from orthogonal to tight-binding basis
         end do
         call g_timer%stop('atomic-sfc')
         !=========================================================================
         !                TEST IF THE CALCULATION IS CONVERGED
         !=========================================================================
         this%converged = this%is_converged(this%mix%delta)
         if (this%converged) then
            if (rank == 0) call g_logger%info('Converged!'//fmt('f12.10', this%mix%delta), __FILE__, __LINE__)

            if (this%hamiltonian%hubbardU_check .and. this%hamiltonian%hubbardJ_check .and. i == 1) then
               call g_logger%info('LDA+U+J module initiates at iteration 2. Continuing calculation...', __FILE__, __LINE__)
               niter = niter + 1
            else if (this%hamiltonian%hubbardU_sc_check) then
               print *, 'Calculates U_eff'
               call this%bands%calc_hubbard_U() ! Calculates Hubbard U
               if ( this%bands%hubbard_u_converged ) then
                  print *, 'The calculated Hubbard U has converged. Exit calculation.'
                  exit
               end if
               niter = niter + 1
            else
               exit
            end if
         else
            if (rank == 0) call g_logger%info('Not converged! Diff= '//fmt('f12.10', this%mix%delta), __FILE__, __LINE__)
            niter = niter + 1
         end if
      end do
   end subroutine run

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Gathers measurables post-scf
   !>
   !> Gathers measurables post-scf
   !---------------------------------------------------------------------------
   subroutine report(this)
      class(self), intent(inout) :: this
      integer :: newunit, iostatus
      integer :: ia, ia_loc
      real(rp), dimension(this%lattice%nrec, 3) :: magmom
      real(rp), dimension(3, this%lattice%nrec) :: mag_for
      ! Open report.out file
      open (newunit=newunit, file='report.out', action='write', iostat=iostatus, status='replace')

      ! Calculate outputs that are not calculated during the SFC run
      call this%bands%calculate_magnetic_torques()

      ! Transfer values across ranks
      magmom = 0.0d0
      mag_for = 0.0d0
      do ia = start_atom, end_atom
         ia_loc = g2l_map(ia)
         magmom(ia, :) = this%symbolic_atom(this%lattice%nbulk + ia)%potential%mom0(:)
         mag_for(:, ia) = this%bands%mag_for(:, ia_loc)
      end do

#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, magmom, product(shape(magmom)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, mag_for, product(shape(mag_for)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      if (rank == 0) then
         call g_logger%info('Calculation finished. Report printed in report.out', __FILE__, __LINE__)
         !===========================================================================
         !                      Total Energy
         !===========================================================================
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(A)') '|                       Total Energy                                      |'
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(a,f20.10)') 'Total energy of system: ', sum(this%symbolic_atom(:)%potential%etot)
         !===========================================================================
         !                       Band Energy
         !===========================================================================
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(A)') '|                       Band Energy                                       |'
         write (newunit, '(A)') '==========================================================================='
         call this%bands%calculate_band_energy()
         write (newunit, '(a,f16.10)') 'Band energy of system: ', this%bands%eband
         !===========================================================================
         !                       Magnetization
         !===========================================================================
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(A)') '|                     Magnetization                                       |'
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(a,f16.10)') 'Total magnetization: ', sum(this%symbolic_atom(this%lattice%nbulk + 1:this%lattice%ntype)%potential%mtot)
         do ia = 1, this%lattice%nrec
            write (newunit, '(a,i4,a,f10.6)') 'Spin moment of atom', ia, ':', norm2(magmom(ia, :))
            write (newunit, '(a,i4,a,3f10.6)') 'Spin moment projections of atom', ia, ':', magmom(ia, :)
            write (newunit, '(a,i4,a,3f16.6)') 'Magnetic force on atom', ia, ':', mag_for(:, ia)
         end do
         !===========================================================================
         !                           Charge
         !===========================================================================
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(A)') '|                       Charge Transfer                                   |'
         write (newunit, '(A)') '==========================================================================='
         do ia = 1, this%lattice%nrec
            write (newunit, '(a,i4,a,f10.6)') 'Charge at atom', ia, ':', sum(this%mix%qia(ia, 1:6))
            write (newunit, '(a,i4,a,f10.6)') 'Charge transfer at atom', ia, ':', this%charge%dq(ia)
         end do
         !===========================================================================
         !                           Fermi Energy
         !===========================================================================
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(A)') '|                       Fermi Energy                                      |'
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(a,f10.6)') 'Fermi energy: ', this%en%fermi
      end if

      if (this%control%hyperfine) then
         do ia = 1, this%lattice%nrec
            call this%symbolic_atom(ia)%potential%print_hyperfine(ia)
         end do
      end if

   end subroutine report

   function is_converged(this, delta_en) result(l)
      class(self), intent(inout) :: this
      real(rp), intent(in) :: delta_en
      logical :: l

      l = delta_en < this%conv_thr
   end function is_converged
   !==============================================================================
   !----- FROM HERE WE HAVE SUBROUTINES AND FUCTIONS RAW COPIED FROM OLD CODE ----
   !==============================================================================
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> LMTO sefconsitency in the atomic part
   !>
   !> LMTO sefconsitency in the atomic part
   !> by Michael Methfessel
   !---------------------------------------------------------------------------
   function lmtst(this, atom) result(QSL)
      implicit none
      class(self), intent(inout) :: this
      class(symbolic_atom), intent(inout) :: atom
      real(rp), dimension(6) :: QSL

      real(rp), dimension(:, :), allocatable :: v
      real(rp), dimension(:), allocatable :: rofi

      ! Iterative variables
      integer :: LMAX, NSP
      integer :: ISP, LP1, I, M, L

      if (rank == 0) call g_logger%info(atom%element%symbol, __FILE__, __LINE__)

      call this%atomsc(atom, v, rofi, "RHO")

      this%VZT(1, 1) = this%VZT(2, 1)
      this%VZT(1, 2) = this%VZT(2, 2)
      LMAX = atom%potential%lmax
      NSP = 2

      call this%RACSI(atom, ROFI, QSL)

      ! ---- MAKE POTENTIAL PARAMETERS TO SAME DNU VALUES ------
      if (atom%potential%ws_r > this%ws_max) then
         do ISP = 1, NSP
            do LP1 = 0, LMAX
               atom%potential%C(LP1, ISP) = 0.d0
               atom%potential%SRDEL(LP1, ISP) = 0.d0
               atom%potential%QPAR(LP1, ISP) = 0.d0
               atom%potential%PPAR(LP1, ISP) = 0.d0
               atom%potential%ENU(LP1, ISP) = 0.d0
               atom%potential%VL(LP1, ISP) = 0.d0
            end do
         end do
      else
         do I = 1, NSP
            do L = 0, LMAX
               atom%potential%PNU(L, I) = atom%potential%PL(L, I)
            end do
         end do
         call this%POTPAR(atom, V, ROFI)
         do I = 1, NSP
            do L = 0, LMAX
               ! write (660, 10002) -L, ENU(L, I), VL(L, I), C(L, I), SRDEL(L, I), QPAR(L, I), PPAR(L, I)
               atom%potential%QPAR(L, I) = 1.0d0/atom%potential%QPAR(L, I)
            end do
         end do
      end if
   end function lmtst

   subroutine atomsc(this, atom, v, rofi, job)
      class(self), intent(inout) :: this
      class(symbolic_atom), intent(inout) :: atom
      type(xc) :: xc_obj
      real(rp), dimension(:, :), allocatable, intent(out) :: v
      real(rp), dimension(:), allocatable, intent(out) :: rofi
      character(LEN=3), intent(in) :: JOB

      integer, parameter :: NCMX = 50
      integer, parameter :: NVMX = 20

      real(rp), dimension(2) :: RVH, RHO0, REPS, RMU, SEV, SEC
      real(rp) :: B_fsm, B, deg, DFCORE, AMGM, EA, RPB, DL, DRHO, TOL, TOLRSQ, BETA, VHRMAX, VSUM, TL, BETA1, VNUCL, SUM, RHO0T, WGT, DRDI, RHOMU, RHOVH, ZVNUCL, OB4PI
      integer :: ISP, ncore, nval, l, nsp, lmax, konf, IFCORE, LCORE, KONFIG, IPR, NR, IR, ITER, NITER, IPR1, II
      logical :: LAST

      real(rp), dimension(NCMX) :: EC
      real(rp), dimension(NVMX) :: EV
      real(rp), dimension(2) :: VRMAX
      real(rp), dimension(:, :), allocatable :: rho_in, rho
      real(rp), dimension(2) :: qval

      ipr = 0
      nsp = 2
      lmax = atom%potential%lmax
      rho_in = atom%rho0(nsp)

      allocate (v, mold=rho_in)
      allocate (rho, mold=rho_in)
      allocate (rofi(size(rho_in(:, 1))))

      nr = size(rofi)
      B_fsm = merge(-atom%mag_cfield(3), real(0.0, rp), this%lattice%control%do_comom)
      xc_obj = xc(this%lattice%control)
      b = atom%B()

      call g_logger%info('XC functional: '//trim(xc_obj%TXCH)//' '//int2str(xc_obj%txc)//real2str(B), __FILE__, __LINE__)

      lmax = atom%potential%lmax
      niter = 80
      qval(:) = 0
      atom%qc = 0
      ncore = 0
      nval = 0
      vhrmax = 0.d0
      tol = 1.d-6
      tolrsq = 1.d-8
      beta = 0.3d0

      do l = 0, lmax
         deg = (2*(2*l + 1))/nsp
         do isp = 1, nsp
            konfig = int(atom%potential%pl(l, isp))
            do konf = l + 1, konfig - 1
               ncore = ncore + 1
               ec(ncore) = -5.d0
               atom%qc = atom%qc + deg
            end do
            nval = nval + 1
            ev(nval) = -0.5d0
            qval(isp) = qval(isp) + atom%potential%ql(1, l, isp)
         end do
      end do
      !----MODIFICATIONS TO INCLUDE FRACTIONARY F OCCUPATION IN THE CORE------
      IFCORE = atom%element%f_core
      DFCORE = REAL(IFCORE)
      call g_logger%info('F core check:'//int2str(ifcore), __FILE__, __LINE__)
      !write (9, *) ´F-core check:´, IFCORE, DFCORE
      if (IFCORE /= 0) then
         LCORE = 3
         DEG = (2*(2*LCORE + 1))/NSP
         do ISP = 1, NSP
            !------DEFINES DEGENERACY OF UP AND DW F-CORE--------------
            !------MAKES DEG=2*DEG se NAO MAGNETICO!!!-----------------
            if (NSP == 1) then
               DEG = DFCORE
            elseif (IFCORE <= 7) then
               if (ISP == 1) then
                  DEG = DFCORE
               end if
               if (ISP == 2) then
                  DEG = 0.0d0
               end if
            else
               if (ISP == 1) then
                  DEG = 7.0d0
               end if
               if (ISP == 2) then
                  DEG = DFCORE - 7.0d0
               end if
            end if
            !-----END OF DEG DEFINITION-------------------------------
            KONFIG = 5
            do KONF = LCORE + 1, KONFIG - 1
               NCORE = NCORE + 1
               EC(NCORE) = -5.d0
               atom%QC = atom%QC + DEG
            end do
         end do
      end if
      !-----------------------END OF MODIFICATIONS----------------------------
      atom%qv = atom%element%atomic_number - atom%qc + .001
      if (NCORE > NCMX) stop "*** CHANGE NCMX IN ATOMSC"
      if (NVAL > NVMX) stop "*** CHANGE NVMX IN ATOMSC"
      atom%dq = atom%qc + qval(1) + qval(2) - atom%element%atomic_number
      amgm = 0.d0
      if (nsp == 2) then
         amgm = qval(1) - qval(2)
      end if
      ! Print statement, if ipr different than 0. Obsolete, therefore commented.
!    if (IPR >= 1) then
!       write (6, 10001) atom%element%atomic_number, atom%potential%ws_r, atom%a, NR, JOB, NCORE, atom%QC, QVAL, atom%DQ, AMGM
!       do ISP = 1, NSP
!          write (6, *) " "
!          if (NSP == 2) then
!             write (6, *) "SPIN", ISP
!          end if
!          write (6, 10002)
!          do l = 1, LMAX
!             DL = TAN(PI*(.5d0-atom%potential%PL(l, ISP)))
!             write (6, 10003) l, atom%potential%PL(l, ISP), DL, (atom%potential%QL(II, l, ISP), II = 1, 3)
!          end do
!       end do
!    end if
      ea = exp(atom%a)
      rpb = b
      do IR = 1, NR
         rofi(IR) = rpb - b
         rpb = rpb*ea
      end do
      if (IPR >= 1) then
         write (6, 10005)
      end if
      if (JOB == "POT") then
         call this%newrho(atom, atom%element%atomic_number, lmax, atom%a, b, nr, rofi, v, rho_in, atom%potential%pl, atom%potential%ql, SEC, SEV, EC, EV, TOLRSQ, NSP, 0)
      elseif (JOB /= "RHO") then
         stop "*** ATOMSC EXPECTS JOB=RHO OR JOB=POT"
      end if
      ! -------- START SELF-CONSISTENCY LOOP ------
      DRHO = 100.d0
      LAST = .false.
      do ITER = 1, NITER
         TL = TOLRSQ
         if (ITER >= 2 .and. DRHO > 2d0) then
            TL = 1.d-3
         end if
         BETA1 = BETA
         if (MOD(ITER, 3) == 2 .and. DRHO < 1.d0) then
            !BETA1 = 0.8d0
            BETA1 = 0.5d0
         end if
       !!    IF(MOD(ITER, 3).EQ.2.AND.DRHO.LT.0.1D0) BETA1=1.0D0
         IPR1 = 0
         if (LAST) then
            IPR1 = IPR
         end if
         call POISS0(atom%element%atomic_number, atom%a, b, rofi, rho_in, NR, VHRMAX, V, RVH, VSUM, NSP)
         VNUCL = V(1, 1)
         !call VXC0SP_old(atom%element%atomic_number, atom%a, B, rofi, rho_in, NR, V, RHO0, REPS, RMU, NSP)
         call this%VXC0SP(xc_obj, atom%element%atomic_number, atom%a, b, rofi, rho_in, NR, V, RHO0, REPS, RMU, NSP, B_fsm)
         call this%NEWRHO(atom, atom%element%atomic_number, lmax, atom%a, b, nr, rofi, v, rho, atom%potential%PL, atom%potential%QL, SEC, SEV, EC, EV, TL, NSP, IPR1)
         DRHO = 0.d0
         SUM = 0.d0
         RHO0T = 0.d0
         do ISP = 1, NSP
            RHO0T = RHO0T + RHO0(ISP)
            do IR = 1, NR
               WGT = 2*(MOD(IR + 1, 2) + 1)/3.d0
               if (IR == 1 .or. IR == NR) then
                  WGT = 1.d0/3.d0
               end if
               DRDI = atom%a*(rofi(IR) + B)
               DRHO = DRHO + WGT*ABS(RHO(IR, ISP) - rho_in(IR, ISP))
               rho_in(IR, ISP) = BETA1*RHO(IR, ISP) + (1.d0 - BETA1)*rho_in(IR, ISP)
               SUM = SUM + WGT*DRDI*RHO(IR, ISP)
            end do
         end do
         if (LAST) exit
         if (IPR >= 3 .or. (IPR >= 2 .and. (DRHO < TOL .or. ITER == 1 .or. ITER == NITER - 1))) then
            write (6, 10004) ITER, SUM, DRHO, VNUCL, RHO0T, VSUM, BETA1
         end if
         if (DRHO < TOL .or. ITER == NITER - 1) then
            LAST = .true.
         end if
      end do
      ! -----------------------------
      atom%potential%RHOEPS = 0.d0
      RHOMU = 0.d0
      atom%potential%SUMEV = 0.d0
      atom%potential%SUMEC = 0.d0
      RHOVH = 0.d0
      do ISP = 1, NSP
         if (NSP == 2 .and. IPR >= 1) then
            write (6, 10006) &
               ISP, V(NR, ISP) - 2*atom%element%atomic_number/atom%potential%ws_r, SEV(ISP), SEC(ISP), RVH(ISP), REPS(ISP), &
               RMU(ISP)
         end if
         atom%potential%RHOEPS = atom%potential%RHOEPS + REPS(ISP)
         RHOMU = RHOMU + RMU(ISP)
         atom%potential%SUMEV = atom%potential%SUMEV + SEV(ISP)
         atom%potential%SUMEC = atom%potential%SUMEC + SEC(ISP)
         RHOVH = RHOVH + RVH(ISP)
      end do
      ZVNUCL = -atom%element%atomic_number*VNUCL
      atom%potential%UTOT = .5d0*(RHOVH + ZVNUCL)
      atom%potential%EKIN = atom%potential%SUMEV + atom%potential%SUMEC - RHOVH - RHOMU
      atom%potential%ETOT = atom%potential%EKIN + atom%potential%UTOT + atom%potential%RHOEPS

      if (IPR >= 1) then
         write (6, 10007) &
            atom%potential%SUMEV, atom%potential%SUMEC, VNUCL, RHOVH, ZVNUCL, atom%potential%UTOT, &
            RHOMU, atom%potential%RHOEPS, atom%potential%EKIN, atom%potential%ETOT
      end if
      VRMAX(1) = -2.d0*atom%element%atomic_number/atom%potential%ws_r
      do ISP = 1, NSP
         VRMAX(1) = VRMAX(1) + V(NR, ISP)/NSP
      end do
      VRMAX(2) = 0.d0
      if (NSP == 2) then
         VRMAX(2) = V(NR, 1) - V(NR, 2)
      end if
      this%TOTSUMEV = this%TOTSUMEV + atom%potential%SUMEV
      OB4PI = 1.d0/(4.d0*PI)
      return

10000 format(i5, f12.5)
10001 format( &
         /" ATOMSC:  Z=", f4.1, "   RMAX=", f9.5, "   A=", f6.3, "   NR=", i4, "  JOB=" &
         , a3/" NCORE=", i3, "   QCORE=", f7.2, "   QVAL=", 2f10.6/ &
         " CHARGE TRANSFER=", f10.6, /" MAGNETIC MOMENT=", f10.6)
10002 format("   l", 8x, "PL", 11x, "DL", 11x, "Q0", 11x, "Q1", 11x, "Q2")
10003 format(i4, 5f13.6)
10004 format(i5, f12.6, 1p, d12.3, 0p, f14.4, d14.4, f14.4, f7.2)
10005 format( &
         /"  ITER     QINT", 9x, "DRHO", 10x, "VH0", 10x, "RHO0", 10x, "VSUM", 5x, "BETA" &
         )
10006 format( &
         " SPIN", i2, ":", /" VRM=", f12.5, "    SEV=", f12.5, "    SEC=", f12.5, / &
         " RVH=", f12.5, "    REP=", f12.5, "    RMU=", f12.5)
10007 format( &
         /" SUMEV=", f12.5, "    SUMEC =", f12.5, "    VNUCL=", f12.5/" RHOVH=" &
         , f12.5, "    ZVNUCL=", f12.5, "    UTOT =", f12.5/" RHOMU=", f12.5, &
         "    RHOEPS=", f12.5, "    EKIN =", f12.5/" ETOT =", f12.5)
   end subroutine atomsc

   ! ----------------------INDICE PARA F DO CORE ----------------------
   function FTYPE(this) result(IFCORE)
      class(self), intent(in) :: this
      integer :: IFCORE

      if (this%IFC(this%IFCC) == 1) then
         IFCORE = 14
      else if (this%IFC(this%IFCC) == 2) then
         IFCORE = 7
      else if (this%IFC(this%IFCC) == 0) then
         IFCORE = 0
      else if ((this%IFC(this%IFCC) .gt. 2) .and. (this%IFC(this%IFCC) .le. 14)) then
         IFCORE = this%IFC(this%IFCC)
      else if (this%IFC(this%IFCC) == 15) then
         IFCORE = 1
      else if (this%IFC(this%IFCC) == 16) then
         IFCORE = 2
      else
         IFCORE = 0
      end if
   end function FTYPE

   subroutine NEWRHO(this, atom, Z, LMAX, A, B, NR, rofi, V, RHO, PL, QL, SUMEC, SUMEV, EC, EV, TOL, NSP, IPR)
      !  MAKES SPHERICAL CHARGE DENSITY FOR A GIVEN POTENTIAL.
      !  RHO IS DETERMINED BY BOUNDARY CONDITIONS PL AND MOMENTS QL.
      !  FOR RMAX GE 10 SETS PHIDOT PHIDOTDOT TO ZERO
      !  SCALAR-RELATIVISTIC VERSION
      !  +++++++  PHDTSR WAS REPLACED BY PHDFSR   9.12.87
      !  +++++++  TOLVAL INTRODUCED FOR NUM DIFFERENTIATION  24.3.88
      !====MODIFICADO PARA F FRACIONARIO NO CAROCO = 7 ELE. =GD ===4/3/94===
      !use common_pri ! REMOVIDO POR LUCAS
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Parameters ..
      integer, parameter :: NRMX = 5001
      !integer, parameter :: NRMX = 501
      !
      !.. Formal Arguments ..
      class(self), intent(inout) :: this
      class(symbolic_atom), intent(inout) :: atom
      integer, intent(in) :: IPR, LMAX, NR, NSP
      real(rp), intent(in) :: A, B, TOL, Z
      real(rp), dimension(*), intent(inout) :: EC, EV
      real(rp), dimension(NR), intent(in) :: rofi
      real(rp), dimension(NSP), intent(inout) :: SUMEC, SUMEV
      real(rp), dimension(LMAX + 1, NSP), intent(in) :: PL
      real(rp), dimension(NR, NSP), intent(in) :: V
      real(rp), dimension(NR, NSP), intent(inout) :: RHO
      real(rp), dimension(3, LMAX + 1, NSP), intent(in) :: QL
      !
      !.. Local Scalars ..
      logical :: FREE
      integer :: IFCORE, IIR, IR, ISI, ISP, IVAL, JR, KONFIG, L, LP1, NN, NRE
      real(rp) :: C, DL, DPHI, DPHIP, E, EB1, EB2, EVAL, FLLP1, GFAC, P, &
                  PHI, PHIP, PI, Q0, Q1, Q2, R, RMAX, RO, ROCRIT, SLO, &
                  SUM, TMC, TOLVAL, VAL
      ! HF additions
      real(rp) :: SHUP, SHDW, SUT, ALFAM1, AL2, HYP, OBPI, RT, WGTA, DRDI
      real(rp), dimension(2) :: SH
      real(rp), dimension(NR) :: HG, DETH
      integer :: IJ
      !
      !.. Local Arrays ..
      integer, dimension(10) :: KONF
      real(rp), dimension(2*NRMX) :: G, GP, GPP
      !
      !.. External Calls ..
      ! external PHDFSR, RHOCOR, RSEQSR
      !
      !.. Intrinsic Functions ..
      intrinsic ATAN, SQRT, TAN
      !
      ! ... Executable Statements ...
      !
      !--------DEFINES POTENTIAL and TRANSFERS it via COMMON---------
      do ISI = 1, NSP
         do IIR = 2, NR
            this%VTN(IIR, ISI) = -2.0d0*Z/rofi(IIR)
            this%VZT(IIR, ISI) = V(IIR, ISI) + this%VTN(IIR, ISI)
         end do
      end do
      !---------------------------------------------------------------
      ROCRIT = 0.002d0
      PI = 4.d0*ATAN(1.d0)
      OBPI = 4.0D0*PI
      EB1 = -50.d0
      EB2 = 50.d0
      TOLVAL = 1.d-12
      C = 274.074d0
      RMAX = rofi(NR)
      FREE = (RMAX > 9.99d0)
      do L = 0, LMAX
         KONF(L + 1) = PL(L + 1, 1)
      end do
      !======MODIFICADO PARA INCLUIR F NO CAROCO================
      IFCORE = atom%element%f_core
      if (IFCORE /= 0) then
         KONF(LMAX + 2) = 5
      end if
      !========FIM  DAS MODIFICACOES NESTE TRECHO================
      do ISP = 1, NSP
         do IR = 1, NR
            RHO(IR, ISP) = 0.d0
         end do
      end do
      call this%RHOCOR(atom, Z, LMAX, KONF, A, B, NR, rofi, V, RHO, G, SUMEC, EC, TOL, NSP, IPR)
      ! ------ LOOP OVER VALENCE STATES -------
      IVAL = 0
      do ISP = 1, NSP
         SUMEV(ISP) = 0.d0
         do LP1 = 1, LMAX + 1
            L = LP1 - 1
            Q0 = QL(1, LP1, ISP)
            Q1 = QL(2, LP1, ISP)
            Q2 = QL(3, LP1, ISP)
            if (Q0 >= 1.d-5) then
               KONFIG = PL(LP1, ISP)
               DL = TAN(PI*(0.5d0 - PL(LP1, ISP)))
               NN = KONFIG - LP1
               IVAL = IVAL + 1
               EVAL = EV(IVAL)
               VAL = RMAX
               SLO = DL + 1.d0
               if (FREE) then
                  VAL = 1.d-30
               end if
               if (FREE) then
                  SLO = -VAL
               end if
               call RSEQSR(EB1, EB2, EVAL, TOLVAL, Z, L, NN, VAL, SLO, V(1, ISP), G, SUM, A, B, &
                           rofi, NR, NRE, 0)
               EV(IVAL) = EVAL
               SUMEV(ISP) = SUMEV(ISP) + EVAL*Q0 + Q1
               !write(*, *) ´atorb: EVAL, Q0, Q1´, EVAL, Q0, Q1
               RO = G(NR)**2
               if (.not. FREE .and. RO < ROCRIT) then
                  !         write (6, 10000) L, NN, NRE, E, RO
                  write (6, 10000) L, NN, NRE, EVAL, RO
               end if
               if (FREE .or. RO < ROCRIT) then
                  do IR = 1, NR
                     GP(IR) = 0.d0
                     GPP(IR) = 0.d0
                  end do
               else
                  VAL = VAL/SQRT(SUM)
                  SLO = SLO/SQRT(SUM)
                  call PHDFSR(Z, L, V(1, ISP), EVAL, A, B, rofi, NR, G, VAL, SLO, GP, GPP, PHI, DPHI, &
                              PHIP, DPHIP, P, TOLVAL, NN)
               end if
               FLLP1 = L*(L + 1)
               do IR = 2, NRE
                  JR = IR + NR
                  R = rofi(IR)
                  TMC = C - (V(IR, ISP) - 2.d0*Z/R - EVAL)/C
                  GFAC = 1.d0 + FLLP1/(TMC*R)**2
                  RHO(IR, ISP) = &
                     RHO(IR, ISP) + Q0*(GFAC*G(IR)**2 + G(JR)**2) + &
                     2.d0*Q1*(GFAC*G(IR)*GP(IR) + G(JR)*GP(JR)) + &
                     Q2*(GFAC*(GP(IR)**2 + G(IR)*GPP(IR)) + GP(JR)**2 + G(JR)*GPP(JR))
                  !-------- FUN2 GIVES the probability density and NOT the charge density---
                  !----of valence orbital LP1(s, p or d) with spin ISP at point r=IR----
                  this%FUN2(IR, LP1, ISP) = (GFAC*G(IR)**2 + G(JR)**2)
               end do
            end if
            !------------BEGIN-ORBITAL CONT. TO HYPERFINE FIELD------------
            if (this%control%hyperfine) then
               if (L == 0) then
                  SHUP = 0.0D0
                  SHDW = 0.0D0
                  SUT = 0.0D0
                  !ALFAM1 = 1.0D0 / 137.036D0
                  ALFAM1 = 2.0D0/C
                  AL2 = ALFAM1**2
                  RT = Z*AL2
                  do IJ = 2, NR
                     WGTA = 2.0D0*(mod(IJ + 1, 2) + 1)/3.0D0
                     if (IJ == 1 .or. IJ == NR) WGTA = 1.0D0/3.0D0
                     DRDI = A*(ROFI(IJ) + B)
                     HG(IJ) = Q0*(GFAC*G(IJ)**2 + 2.0D0*Q1* &
                                  (GFAC*G(IJ)*GP(IJ)) + Q2*(GFAC*(GP(IJ)**2 + G(IJ)*GPP(IJ))))
                     HG(IJ) = HG(IJ)/(OBPI*ROFI(IJ)**2)
                     DETH(IJ) = (ROFI(IJ) + (RT/2.0D0))**2
                     DETH(IJ) = (RT/2.0D0)/DETH(IJ)
                     SUT = SUT + WGTA*DRDI*DETH(IJ)
                     HG(IJ) = WGTA*DRDI*DETH(IJ)*HG(IJ)
                     if (ISP == 1) then
                        SHUP = SHUP + HG(IJ)
                        SH(1) = SHUP
                     else
                        SHDW = SHDW + HG(IJ)
                        SH(2) = SHDW
                     end if
                  end do
                  HYP = 52.42_rp*(SH(1) - SH(2))
                  if (ISP == 2) then
                     !write(8, *) SH(1), SH(2), HYP
                     !write(8, ´(a7, f10.6)´) "Hval:", HYP
                     atom%potential%hyper_field(2) = hyp
                  end if
               end if
            end if
            !--------------END-ORBITAL CONT. TO HYPERFINE FIELD------------
         end do
      end do
      return
      !
      ! ... Format Declarations ...
      !
10000 format(" ** PHP, PHPP SET TO ZERO,  L, NN, NRE, E, RHO=", 3i5, 2f8.4)
   end subroutine NEWRHO

   subroutine RHOCOR(this, atom, Z, LMAX, KONFIG, A, B, NR, rofi, V, RHO, G, SUMEC, EC, TOL, NSP, IPR)
      !  ADDS THE (SPHERICAL) CHARGE DENSITY OF THE CORE STATES
      !  SCALAR RELATIVISTIC VERSION
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      class(self), intent(inout) :: this
      class(symbolic_atom), intent(inout) :: atom
      integer, intent(in) :: IPR, LMAX, NR, NSP
      real(rp), intent(in) :: A, B, TOL, Z
      integer, dimension(*), intent(in) :: KONFIG
      real(rp), dimension(*), intent(inout) :: EC
      real(rp), dimension(NR), intent(in) :: rofi
      real(rp), dimension(NSP), intent(inout) :: SUMEC
      real(rp), dimension(NR, 2), intent(inout) :: G
      real(rp), dimension(NR, NSP), intent(in) :: V
      real(rp), dimension(NR, NSP), intent(inout) :: RHO
      !
      !.. Local Scalars ..
      integer :: ICORE, IFCORE, IR, ISP, KONF, L, LL, LP1, NODES, NRE
      real(rp) :: C, DEG, DFCORE, DLML, E1, E2, ECOR0, ECORE, FLLP1, &
                  GFAC, QCORE, R, RHORIM, RMAX, RORIM, SLO, SUM, TMC, &
                  VAL, YYY
      ! HF additions
      real(rp) :: SHUP, SHDW, SUT, ALFAM1, AL2, HYP, OBPI, RT, WGTA, DRDI, HCORE, PI
      real(rp), dimension(10, 2) :: SH
      real(rp), dimension(10) :: HFC
      real(rp), dimension(NR) :: HG, DETH
      integer :: ij, ii
      !
      !.. External Calls ..
      ! external RSEQSR
      !
      !.. Intrinsic Functions ..
      ! intrinsic REAL, SQRT
      !
      ! ... Executable Statements ...
      !
      !    .   PHI(0:5), PSI(0:5), SUMEC(NSP), EC(1)
      RMAX = rofi(NR)
      E1 = -2.5*Z*Z - 5.d0
      E2 = 20.d0
      C = 274.074d0
      PI = 4.d0*ATAN(1.d0)
      OBPI = 4.0D0*PI
      !C = 274.071979d00
      ! ----------------------------
      ICORE = 0
      do ISP = 1, NSP
         SUMEC(ISP) = 0.d0
         QCORE = 0.d0
         RHORIM = 0.d0
         if (IPR > 0) then
            write (6, 10000) ISP
         end if
         do LP1 = 1, LMAX + 1
            L = LP1 - 1
            DEG = (2*(2*L + 1))/NSP
            do KONF = LP1, KONFIG(LP1) - 1
               ICORE = ICORE + 1
               NODES = KONF - LP1
               ECOR0 = EC(ICORE)
               VAL = 1.d-30
               SLO = -VAL
               call RSEQSR(E1, E2, ECOR0, TOL, Z, L, NODES, VAL, SLO, V(1, ISP), G, SUM, A, B, rofi, NR, NRE, 0)
               ECORE = ECOR0
               YYY = ECOR0 - V(NR, ISP) + 2*Z/RMAX
               if (NRE == NR .and. YYY < 0.d0) then
                  DLML = -1.d0 - SQRT(-YYY)*RMAX
                  do LL = 1, L
                     DLML = -YYY*RMAX*RMAX/DLML - (2*LL + 1)
                  end do
                  SLO = VAL*(DLML + L + 1)/RMAX
                  call RSEQSR(E1, E2, ECORE, TOL, Z, L, NODES, VAL, SLO, V(1, ISP), G, SUM, A, B, &
                              rofi, NR, NRE, 0)
               end if
               EC(ICORE) = ECORE
               FLLP1 = L*(L + 1)
               do IR = 2, NRE
                  R = rofi(IR)
                  TMC = C - (V(IR, ISP) - 2.d0*Z/R - ECORE)/C
                  GFAC = 1.d0 + FLLP1/(TMC*R)**2
                  RHO(IR, ISP) = RHO(IR, ISP) + DEG*(GFAC*G(IR, 1)**2 + G(IR, 2)**2)
               end do
               RORIM = 0.d0
               if (NRE == NR) then
                  RORIM = DEG*(GFAC*G(NR, 1)**2 + G(NR, 2)**2)/RMAX**2
               end if
               RHORIM = RHORIM + RORIM
               QCORE = QCORE + DEG
               SUMEC(ISP) = SUMEC(ISP) + DEG*ECORE
               if (IPR > 0) then
                  write (6, 10001) KONF, L, NODES, ECOR0, ECORE, NRE, RORIM
               end if
               if (this%control%hyperfine) then
                  !------------BEGIN-CORE CONT. TO HYPERFINE FIELD------------
                  IF (L == 0) THEN
                     SHUP = 0.0D0
                     SHDW = 0.0D0
                     SUT = 0.0D0
                     ALFAM1 = 2.0/C
                     AL2 = ALFAM1**2
                     RT = Z*AL2
                     DO 22 IJ = 2, NRE
                        !      IF(IJ.EQ.NRE)WRITE(8,758) KONF,L,NODES,ECOR0,ECORE,NRE,RORIM
                        WGTA = 2*(MOD(IJ + 1, 2) + 1)/3.0D0
                        IF (IJ == 1 .OR. IJ == NRE) WGTA = 1.0D0/3.0D0
                        DRDI = A*(ROFI(IJ) + B)
                        HG(IJ) = (GFAC*G(IJ, 1)**2)
                        HG(IJ) = HG(IJ)/(OBPI*ROFI(IJ)**2)
                        DETH(IJ) = (ROFI(IJ) + (RT/2.0D0))**2
                        DETH(IJ) = (RT/2.0D0)/DETH(IJ)
                        SUT = SUT + WGTA*DRDI*DETH(IJ)
                        HG(IJ) = WGTA*DRDI*DETH(IJ)*HG(IJ)
                        IF (ISP == 1) THEN
                           !      SHUP=SHUP+WGTA*DRDI*DETH(IJ)*HG(IJ)
                           SHUP = SHUP + HG(IJ)
                           IF (IJ == NRE) SH(KONF, ISP) = SHUP
                        ELSE
                           !      SHDW=SHDW+WGTA*DRDI*DETH(IJ)*HG(IJ)
                           SHDW = SHDW + HG(IJ)
                           IF (IJ == NRE) SH(KONF, ISP) = SHDW
                        END IF
                        ! WRITE(8,777)ROFI(IJ),HG(IJ),SUT,DETH(IJ)
22                   END DO
                  END IF
                  !--------------END-CORE CONT. TO HYPERFINE FIELD------------
               end if

            end do
         end do

         if (this%control%hyperfine) then
            IF (ISP == 2) THEN
               HCORE = 0.0D0
               DO II = 1, KONFIG(1) - 1
                  HFC(II) = SH(II, 1) - SH(II, 2)
                  WRITE (8, *) II, HFC(II)
                  HCORE = HCORE + HFC(II)
               END DO
               HCORE = 52.42_rp*HCORE
               atom%potential%hyper_field(1) = HCORE
               !WRITE(8,´(a7,3f12.6)´)"Hcore: ",HCORE,(SH(1,1)-SH(1,2)) &
               !,(SH(1,1)-SH(1,2))*914.7744
            END IF
         end if

         !====MODIFICADO PARA INCLUIR F FRACIONARIO NO CAROCO==========
         IFCORE = atom%element%f_core
         DFCORE = REAL(IFCORE)
         if (IFCORE /= 0) then
            LP1 = LMAX + 2
            L = LP1 - 1
            !     DEG=(2*(2*L+1))/NSP
            !=======DEFINE DEGENERESCENCIA DEG PARA SPINS UP E DW DO F-CORE=====
            !------------MAKES DEG=2*DEG se nao magnetico!!!----------------
            if (NSP == 1) then
               DEG = DFCORE
            elseif (IFCORE <= 7) then
               if (ISP == 1) DEG = DFCORE
               if (ISP == 2) DEG = 0.0d0
            else
               if (ISP == 1) DEG = 7.0d0
               if (ISP == 2) DEG = DFCORE - 7.0d0
            end if
            !==============FIM DA DEFINICAO DE DEG ===================
            do KONF = LP1, KONFIG(LP1) - 1
               ICORE = ICORE + 1
               NODES = KONF - LP1
               ECOR0 = EC(ICORE)
               VAL = 1.d-30
               SLO = -VAL
               call RSEQSR(E1, E2, ECOR0, TOL, Z, L, NODES, VAL, SLO, V(1, ISP), G, SUM, A, B, rofi, NR, NRE, 0)
               ECORE = ECOR0
               YYY = ECOR0 - V(NR, ISP) + 2*Z/RMAX
               if (NRE == NR .and. YYY < 0.d0) then
                  DLML = -1.d0 - SQRT(-YYY)*RMAX
                  do LL = 1, L
                     DLML = -YYY*RMAX*RMAX/DLML - (2*LL + 1)
                  end do
                  SLO = VAL*(DLML + L + 1)/RMAX
                  call RSEQSR(E1, E2, ECORE, TOL, Z, L, NODES, VAL, SLO, V(1, ISP), G, SUM, A, B, &
                              rofi, NR, NRE, 0)
               end if
               EC(ICORE) = ECORE
               FLLP1 = L*(L + 1)
               do IR = 2, NRE
                  R = rofi(IR)
                  TMC = C - (V(IR, ISP) - 2.d0*Z/R - ECORE)/C
                  GFAC = 1.d0 + FLLP1/(TMC*R)**2
                  RHO(IR, ISP) = RHO(IR, ISP) + DEG*(GFAC*G(IR, 1)**2 + G(IR, 2)**2)
               end do
               RORIM = 0.d0
               if (NRE == NR) then
                  RORIM = DEG*(GFAC*G(NR, 1)**2 + G(NR, 2)**2)/RMAX**2
               end if
               RHORIM = RHORIM + RORIM
               QCORE = QCORE + DEG
               SUMEC(ISP) = SUMEC(ISP) + DEG*ECORE
               if (IPR > 0) then
                  write (6, 10001) KONF, L, NODES, ECOR0, ECORE, NRE, RORIM
               end if
            end do
         end if
         !====== FIM DAS MODIFICACOES NESTE TRECHO ============================
         if (IPR > 0) then
            write (6, 10002) QCORE, RHORIM, SUMEC(ISP)
         end if
      end do
      return
      !
      ! ... Format Declarations ...
      !
10000 format( &
         /" ROCORE:  SPIN", i2/" KONF  L  NN       ECOR0", 9x, "ECORE", 9x, "NRE" &
         , 7x, "RHONR")
10001 format(3i4, 2f14.5, i10, 1pe14.2)
10002 format( &
         " CORE CHARGE=", f9.3, "    DENSITY AT RMAX=", 1pe11.3, "    SUMEC=" &
         , 0pf12.5)
   end subroutine RHOCOR

   subroutine RSEQSR(EB1, EB2, E, TOL, Z, L, NOD, VAL, SLO, V, G, Q, A, B, rofi, NR, NRE, IPR)
      !  SOLVES RADIAL SCALAR RELATIVISTIC EQ. TO GIVEN BCS AND NODES.
      !  VAL, SLO ARE BC´S FOR LARGE COMPONENT U(R) WITH PSI=(U/R)*YLM.
      !  OUTPUT WAVEFCT IS NORMALIZED TO 1.    (M.METHFESSEL 2/87)
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      integer, intent(in) :: IPR, L, NOD, NR
      integer, intent(inout) :: NRE
      real(rp), intent(in) :: A, B, EB1, EB2, SLO, TOL, VAL, Z
      real(rp), intent(inout) :: E, Q
      real(rp), dimension(*), intent(in) :: rofi, V
      real(rp), dimension(NR, 2), intent(inout) :: G
      !
      !.. Local Scalars ..
      integer :: K, K2, KC, NCTP, NCTP0, NIT, NITMAX, NOD1, NOD2, NODE, NSAVE
      real(rp) :: C, DE, E1, E2, FAC, FLLP1, R, RATIO, RE, RHOK, SLO1, &
                  SLO2, SLOP, TMCR, VAL1, VAL2, VALU, WGT, XMIN, XRIM
      !
      !.. External Calls ..
      ! external FCTP, FCTP0, RSQSR1, RSQSR2
      !
      !.. Intrinsic Functions ..
      intrinsic ABS, LOG, MAX, MIN, MOD, SQRT
      !
      ! ... Executable Statements ...
      !
      NITMAX = 400
      C = 274.074d0
      if (IPR >= 1) then
         write (6, 10000) Z, L, NOD, VAL, SLO, EB1, EB2
      end if
      NIT = 0
      E1 = EB1
      E2 = EB2
      call FCTP0(L, rofi, V, Z, NR, A, B, NCTP0, XRIM, XMIN, NSAVE)
      do
         ! ----- START ITERATIONS TO FIND ENERGY ------
         NIT = NIT + 1
         if (NIT > NITMAX) then
            write (6, 10001) NITMAX, L, NOD, NODE, E
         end if
         if (NIT > NITMAX) return
         if (E <= E1 .or. E >= E2) then
            E = .5d0*(E1 + E2)
         end if
         call FCTP(E, NCTP, NCTP0, XRIM, XMIN, NSAVE, L, rofi, V, Z, NR, A, B)
         RE = 15.d0*rofi(NCTP)
         NRE = LOG(RE/B + 1.d0)/A + 1.d0
         NRE = (NRE/2)*2 + 1
         NRE = MAX(35, MIN(NRE, NR))
         VALU = VAL
         SLOP = SLO
         if (NRE < NR) then
            VALU = 1.d-5
         end if
         if (NRE < NR) then
            SLOP = -1.d-5
         end if
         K2 = 30
         if (NOD == 0) then
            K2 = NRE/3
         end if
         if (VALU*SLOP > 0.d0 .and. NOD == 0) then
            K2 = NRE - 10
         end if
         call RSQSR2(E, L, Z, V, NRE, K2, VALU, SLOP, G, VAL2, SLO2, NOD2, KC, A, B, rofi, NR)
         call RSQSR1(E, L, Z, V, KC, G, VAL1, SLO1, NOD1, A, B, rofi, NR)
         NODE = NOD1 + NOD2
         if (NODE /= NOD) then
            if (IPR >= 2 .or. NIT >= NITMAX - 5) then
               write (6, 10002) NIT, L, NODE, NRE, KC, E1, E, E2
            end if
            if (NODE > NOD) then
               E2 = E
            end if
            if (NODE < NOD) then
               E1 = E
            end if
            E = .5d0*(E1 + E2)
         else
            ! ---- CORRECTION TO EIGENVALUE ----
            RATIO = VAL2/VAL1
            Q = 0.d0
            do K = 2, KC
               Q = Q + (rofi(K) + B)*G(K, 1)**2
            end do
            Q = Q*RATIO*RATIO
            do K = KC + 1, NRE
               Q = Q + (rofi(K) + B)*G(K, 1)**2
            end do
            Q = A*(Q - .5d0*(rofi(NRE) + B)*G(NRE, 1)**2)
            DE = -VAL2*(SLO2 - RATIO*SLO1)/Q
            if (IPR >= 2 .or. NIT >= NITMAX - 5) then
               write (6, 10002) NIT, L, NODE, NRE, KC, E1, E, E2, DE
            end if
            if (DE > 0d0) then
               E1 = E
            end if
            if (DE < 0d0) then
               E2 = E
            end if
            E = E + DE
            if (ABS(DE) <= TOL .or. NIT >= NITMAX) exit
         end if
      end do
      ! -----  NORMALIZE G -------
      FLLP1 = L*(L + 1)
      E = E - DE
      do K = 1, KC
         G(K, 1) = G(K, 1)*RATIO
         G(K, 2) = G(K, 2)*RATIO
      end do
      Q = 0.d0
      do K = 2, NRE
         R = rofi(K)
         WGT = (MOD(K + 1, 2) + 1)*(R + B)
         TMCR = (C - (V(K) - 2.d0*Z/R - E)/C)*R
         RHOK = G(K, 1)**2*(1.d0 + FLLP1/TMCR**2) + G(K, 2)**2
         Q = Q + WGT*RHOK
      end do
      Q = (Q - 0.5d0*WGT*RHOK)*A*2.d0/3.d0
      FAC = 1.d0/SQRT(Q)
      do K = 1, NRE
         G(K, 1) = G(K, 1)*FAC
         G(K, 2) = G(K, 2)*FAC
      end do
      do K = NRE + 1, NR
         G(K, 1) = 0.d0
         G(K, 2) = 0.d0
      end do
      if (IPR >= 1 .or. NIT >= NITMAX) then
         write (123, 10003) E, Q, NR, NRE, KC, DE, NIT
      end if
      return
      !
      ! ... Format Declarations ...
      !
10000 format( &
         " RSEQSR:  Z=", f5.1, "  L=", i1, "  NOD=", i1, "  BC=", 2f7.3, "  E1, E2=" &
         , 2f8.2)
10001 format( &
         " *** RSEQ  NIT GT", i4, " AND BAD NODES, L=", i2, "  NOD=", i3, "  NODE=" &
         , i3, "   E=", f14.5)
10002 format(3i3, 2i4, 3f16.8, 1p, d13.4)
10003 format( &
         " E=", f13.5, "  Q=", f10.5, "   NR/NRE/KC=", 3i4, "   DE=", 1p, d10.2, &
         "  NIT=", i2)
   end subroutine RSEQSR

   subroutine PHDFSR(Z, L, V, E, A, B, rofi, NR, G, VAL, SLO, GP, GPP, PHI, DPHI, PHIP, DPHIP, P, TOL, NN)
      !  MAKE PHIDOT, PHIDOTDOT TO GIVEN E (SCALAR RELATIVISTIC)
      !  BY NUMERICAL DIFFERENTIATION
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      integer, intent(in) :: L, NN, NR
      real(rp), intent(in) :: A, B, E, SLO, TOL, VAL, Z
      real(rp), intent(out) :: DPHI, DPHIP, P, PHI, PHIP
      real(rp), dimension(*), intent(in) :: G, rofi, V
      real(rp), dimension(*), intent(inout) :: GP, GPP
      !
      !.. Local Scalars ..
      integer :: IR, NRE
      real(rp) :: DDDE, DDL, DELE, DEN, E1, E2, EB1, EB2, GPIR, RMAX, &
                  SLO1, SLO2, SLP, SRDRDI, SUM1, SUM2, VAL1, VAL2, &
                  VLP, WP0, WP1, WP2, WPP0, WPP1, WPP2, X1, X2
      !
      !.. External Calls ..
      ! external GINTSR, RSEQSR
      !
      !.. Intrinsic Functions ..
      intrinsic SQRT
      !
      ! ... Executable Statements ...
      !
      RMAX = rofi(NR)
      EB1 = -50.d0
      EB2 = 15.d0
      DELE = 0.003
      DDDE = -RMAX/G(NR)**2
      DDL = DELE*DDDE
      SLO1 = SLO - DDL*VAL/RMAX
      SLO2 = SLO + DDL*VAL/RMAX
      E1 = E
      E2 = E
      call RSEQSR(EB1, EB2, E1, TOL, Z, L, NN, VAL, SLO1, V, GP, SUM1, A, B, rofi, NR, NRE, 0)
      VAL1 = VAL/SQRT(SUM1)
      SLO1 = SLO1/SQRT(SUM1)
      call RSEQSR(EB1, EB2, E2, TOL, Z, L, NN, VAL, SLO2, V, GPP, SUM2, A, B, rofi, NR, NRE, 0)
      VAL2 = VAL/SQRT(SUM2)
      SLO2 = SLO2/SQRT(SUM2)
      X1 = E1 - E
      X2 = E2 - E
      DEN = X1*X2*(X1 - X2)
      WP0 = (X2**2 - X1**2)/DEN
      WP1 = -X2**2/DEN
      WP2 = X1**2/DEN
      WPP0 = 2.d0*(X1 - X2)/DEN
      WPP1 = 2.d0*X2/DEN
      WPP2 = -2.d0*X1/DEN
      do IR = 1, NR*2
         GPIR = WP0*G(IR) + WP1*GP(IR) + WP2*GPP(IR)
         GPP(IR) = WPP0*G(IR) + WPP1*GP(IR) + WPP2*GPP(IR)
         GP(IR) = GPIR
      end do
      VLP = WP0*VAL + WP1*VAL1 + WP2*VAL2
      SLP = WP0*SLO + WP1*SLO1 + WP2*SLO2
      call GINTSR(GP, GP, A, B, NR, Z, E, L, V, rofi, P)
      SRDRDI = SQRT(A*(RMAX + B))
      PHI = VAL/RMAX
      DPHI = SLO/RMAX - VAL/RMAX/RMAX
      PHIP = VLP/RMAX
      DPHIP = (SLP - VLP/RMAX)/RMAX
   end subroutine PHDFSR

   subroutine GINTSR(G1, G2, A, B, NR, Z, E, L, V, rofi, SUM)
      !  EVALUATES SR-SP OF G1, G2
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      integer, intent(in) :: L, NR
      real(rp), intent(in) :: A, B, E, Z
      real(rp), intent(out) :: SUM
      real(rp), dimension(NR), intent(in) :: rofi, V
      real(rp), dimension(NR, 2), intent(in) :: G1, G2
      !
      !.. Local Scalars ..
      integer :: IR
      real(rp) :: C, FLLP1, GFAC, R, TMC
      !
      ! ... Executable Statements ...
      !
      FLLP1 = L*(L + 1)
      C = 274.074d0
      SUM = 0.d0
      do IR = 2, NR - 1, 2
         R = rofi(IR)
         TMC = C - (V(IR) - 2.d0*Z/R - E)/C
         GFAC = 1.d0 + FLLP1/(TMC*R)**2
         SUM = SUM + (R + B)*(G1(IR, 1)*G2(IR, 1)*GFAC + G1(IR, 2)*G2(IR, 2))
      end do
      SUM = SUM + SUM
      do IR = 3, NR - 2, 2
         R = rofi(IR)
         TMC = C - (V(IR) - 2.d0*Z/R - E)/C
         GFAC = 1.d0 + FLLP1/(TMC*R)**2
         SUM = SUM + (R + B)*(G1(IR, 1)*G2(IR, 1)*GFAC + G1(IR, 2)*G2(IR, 2))
      end do
      SUM = SUM + SUM
      R = rofi(NR)
      TMC = C - (V(NR) - 2.d0*Z/R - E)/C
      GFAC = 1.d0 + FLLP1/(TMC*R)**2
      SUM = SUM + (R + B)*(G1(NR, 1)*G2(NR, 1)*GFAC + G1(NR, 2)*G2(NR, 2))
      SUM = SUM*A/3.d0
   end subroutine GINTSR

   subroutine FCTP0(L, rofi, V, Z, NR, A, B, NCTP0, XRIM, XMIN, NSAVE)
      !  INITIALIZE THINGS FOR FCTP, WHICH FINDS CLASSICAL TURNING POINT
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      integer, intent(in) :: L, NR
      integer, intent(out) :: NCTP0, NSAVE
      real(rp), intent(in) :: Z
      real(rp), intent(out) :: XMIN, XRIM
      real(rp) :: A, B
      real(rp), dimension(NR), intent(in) :: rofi, V
      !
      !.. Local Scalars ..
      integer :: IR
      real(rp) :: FLLP1, R, X, XLAST
      !
      ! ... Executable Statements ...
      !
      FLLP1 = L*(L + 1)
      R = rofi(10)
      X = FLLP1/R/R - 2.d0*Z/R + V(10)
      IR = 10
      do
         IR = IR + 1
         XLAST = X
         R = rofi(IR)
         X = FLLP1/R/R - 2.d0*Z/R + V(IR)
         if (X > XLAST .or. IR >= NR) exit
      end do
      NCTP0 = IR - 1
      XMIN = XLAST
      R = rofi(NR)
      XRIM = FLLP1/R/R - 2.d0*Z/R + V(NR)
      if (XMIN >= XRIM - 3.d0) then
         NCTP0 = NR
      end if
      if (XMIN >= XRIM - 3.d0) then
         XMIN = XRIM
      end if
      NSAVE = (NCTP0 + NR)/2
   end subroutine FCTP0

   subroutine FCTP(E, NCTP, NCTP0, XRIM, XMIN, NSAVE, L, rofi, V, Z, NR, A, B)
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      integer, intent(in) :: L, NCTP0, NR
      integer, intent(inout) :: NSAVE
      integer, intent(out) :: NCTP
      real(rp), intent(in) :: A, B, E, XMIN, XRIM, Z
      real(rp), dimension(*), intent(in) :: rofi, V
      ! real(rp), dimension(NR), intent(in) :: rofi, V
      !
      !.. Local Scalars ..
      integer :: IREP, N1, N2, NLAST, NTRY
      real(rp) :: DFDR, DVDR, FLLP1, FNTRY, FOFR, R, RTRY, VME
      !
      !.. Intrinsic Functions ..
      intrinsic LOG, MAX
      !
      ! ... Executable Statements ...
      !
      FLLP1 = L*(L + 1)
      if (NCTP0 == NR) then
         NCTP = NR
         return
      end if
      if (E > XRIM) then
         NCTP = NR
         return
      end if
      if (E < XMIN) then
         NCTP = 2
         return
      end if
      N1 = NCTP0
      N2 = NR
      NCTP = NSAVE
      NLAST = -10
      do IREP = 1, 20
         if (NCTP > N2 .or. NCTP < N1) then
            NCTP = (N1 + N2 + 1)/2
         end if
         ! ++++++++ NCTP CAN BE BIGGER THAN NR !!!!!
         R = rofi(NCTP)
         VME = (V(NCTP) - E)
         DVDR = (V(NCTP + 1) - V(NCTP - 1))/(2.d0*A*(R + B))
         FOFR = FLLP1/R/R - 2.d0*Z/R + VME
         DFDR = -2.d0*FLLP1/R/R/R + 2.d0*Z/R/R + DVDR
         RTRY = R - FOFR/DFDR
         RTRY = MAX(RTRY, rofi(2))
         FNTRY = LOG(RTRY/B + 1.d0)/A + 1.d0
         NTRY = FNTRY + .5d0
         !|    WRITE(6, 810) IREP, N1, NCTP, N2, FNTRY, NTRY
         !|810 FORMAT(I6, ´   N1, NCTP, N2=´, 3I5, ´   NTRY=´, F8.3, I6)
         if (NLAST == NCTP) exit
         if (FOFR > 0.d0) then
            N2 = NCTP
         end if
         if (FOFR < 0.d0) then
            N1 = NCTP
         end if
         NLAST = NCTP
         NCTP = NTRY
      end do
      if (NCTP == NCTP0 + 1) then
         NCTP = 2
      end if
      NSAVE = NCTP
   end subroutine FCTP

   subroutine RSQSR1(E, L, Z, V, KR, G, VAL, SLO, NN, A, B, rofi, NR)
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      integer, intent(in) :: KR, L, NR
      integer, intent(out) :: NN
      real(rp), intent(in) :: A, B, E, Z
      real(rp), intent(out) :: SLO, VAL
      real(rp), dimension(*), intent(in) :: rofi, V
      real(rp), dimension(NR, 2), intent(inout) :: G
      !
      !.. Local Scalars ..
      integer :: K
      real(rp) :: AA, B1, B2, C, DET, DF1, DF2, DF3, DG1, DG2, DG3, DRDI, &
                  F0, FLLP1, G0, H83, PHI, R, R1, R2, R3, R83SQ, S, SF, U, &
                  X, Y, ZZ
      !
      !.. Local Arrays ..
      real(rp), dimension(2, 3) :: D
      !
      !.. Intrinsic Functions ..
      intrinsic SQRT
      !
      ! ... Executable Statements ...
      !
      NN = 0
      ZZ = Z + Z
      C = 274.074d0
      FLLP1 = L*(L + 1.d0)
      R83SQ = 64.d0/9.d0
      R1 = 1.d0/9.d0
      R2 = -5.d0*R1
      R3 = 19.d0*R1
      H83 = 8.d0/3.d0
      ! ------ APPROXIMATE G, F BY LEADING TERM NEAR ZERO ----
      if (Z < 0.9d0) then
         S = L + 1
         SF = L
         G0 = 1.d0
         F0 = L/C
      else
         AA = ZZ/C
         S = SQRT(FLLP1 + 1.d0 - AA*AA)
         SF = S
         G0 = 1.d0
         F0 = G0*(S - 1.d0)/AA
      end if
      G(1, 1) = 0.d0
      G(1, 2) = 0.d0
      do K = 2, 4
         R = rofi(K)
         DRDI = A*(R + B)
         G(K, 1) = (R**S)*G0
         G(K, 2) = (R**SF)*F0
         D(1, K - 1) = DRDI*G(K, 1)*S/R
         D(2, K - 1) = DRDI*G(K, 2)*SF/R
      end do
      ! ----- INTEGRATE OVER REST OF POINTS ------
      DG1 = D(1, 1)
      DG2 = D(1, 2)
      DG3 = D(1, 3)
      DF1 = D(2, 1)
      DF2 = D(2, 2)
      DF3 = D(2, 3)
      do K = 5, KR
         R = rofi(K)
         DRDI = A*(R + B)
         PHI = (E + ZZ/R - V(K))*DRDI/C
         U = DRDI*C + PHI
         X = -DRDI/R
         Y = -FLLP1*X*X/U + PHI
         DET = R83SQ - X*X + U*Y
         B1 = G(K - 1, 1)*H83 + R1*DG1 + R2*DG2 + R3*DG3
         B2 = G(K - 1, 2)*H83 + R1*DF1 + R2*DF2 + R3*DF3
         G(K, 1) = (B1*(H83 - X) + B2*U)/DET
         G(K, 2) = (B2*(H83 + X) - B1*Y)/DET
         if (G(K, 1)*G(K - 1, 1) < 0d0) then
            NN = NN + 1
         end if
         DG1 = DG2
         DG2 = DG3
         DG3 = U*G(K, 2) - X*G(K, 1)
         DF1 = DF2
         DF2 = DF3
         DF3 = X*G(K, 2) - Y*G(K, 1)
      end do
      VAL = G(KR, 1)
      SLO = DG3/(A*(rofi(KR) + B))
   end subroutine RSQSR1

   subroutine RSQSR2(E, L, Z, V, K1, K2, VAL1, SLO1, G, VAL, SLO, NN, KC, A, B, rofi, NR)
      !  INTEGRATE THE SCALAR RELATIVISTIC EQUATION INWARD FROM K1.
      !  CUTOFF KC IS CHOSEN AT FIRST MAXIMUM (BUT KC GE K2)
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      integer, intent(in) :: K1, K2, L, NR
      integer, intent(out) :: KC, NN
      real(rp), intent(in) :: A, B, E, SLO1, VAL1, Z
      real(rp), intent(out) :: SLO, VAL
      real(rp), dimension(*), intent(in) :: rofi, V
      real(rp), dimension(NR, 2), intent(inout) :: G
      !
      !.. Local Scalars ..
      integer :: I, K, KP1
      real(rp) :: AF1, AF2, AF3, AG1, AG2, AG3, B1, B2, C, DET, DF1, DF2, &
                  DF3, DG1, DG2, DG3, DR, EA, FF, FLLP1, GG, H83, PHI, Q, &
                  R, R1, R2, R3, R83SQ, RPB, U, VB, X, Y, ZZ
      !
      !.. Local Arrays ..
      real(rp), dimension(2, 3) :: D
      !
      !.. Intrinsic Functions ..
      intrinsic EXP, MOD, SQRT
      !
      ! ... Executable Statements ...
      !
      NN = 0
      ZZ = Z + Z
      C = 274.074d0
      FLLP1 = L*(L + 1.)
      R83SQ = 64.d0/9.d0
      R1 = 1.d0/9.d0
      R2 = -5.d0*R1
      R3 = 19.d0*R1
      H83 = -8.d0/3.d0
      EA = EXP(A)
      ! ---- FIRST POINT ------
      RPB = B*EXP(A*K1 - A)
      R = RPB - B
      DR = A*RPB
      PHI = (E + ZZ/R - V(K1))*DR/C
      U = DR*C + PHI
      X = -DR/R
      Y = -FLLP1*X*X/U + PHI
      G(K1, 1) = VAL1
      G(K1, 2) = (SLO1*DR + X*VAL1)/U
      Q = 1.d0/SQRT(EA)
      AG1 = SLO1*DR
      AF1 = X*G(K1, 2) - Y*G(K1, 1)
      K = K1
      DG3 = AG1
      if (K2 /= K1) then
         ! ----- RUNGE-KUTTA FOR NEXT THREE POINTS -----
         do I = 1, 3
            KP1 = K
            K = K - 1
            RPB = RPB*Q
            DR = RPB*A
            R = RPB - B
            GG = G(KP1, 1) - .5d0*AG1
            FF = G(KP1, 2) - .5d0*AF1
            VB = (3.d0*V(KP1) + 6.d0*V(K) - V(K - 1))*.125d0
            PHI = (E + ZZ/R - VB)*DR/C
            U = DR*C + PHI
            X = -DR/R
            Y = -FLLP1*X*X/U + PHI
            AG2 = U*FF - X*GG
            AF2 = X*FF - Y*GG
            GG = G(KP1, 1) - .5d0*AG2
            FF = G(KP1, 2) - .5d0*AF2
            AG3 = U*FF - X*GG
            AF3 = X*FF - Y*GG
            RPB = RPB*Q
            DR = A*RPB
            R = RPB - B
            PHI = (E + ZZ/R - V(K))*DR/C
            U = DR*C + PHI
            X = -DR/R
            Y = -FLLP1*X*X/U + PHI
            GG = G(KP1, 1) - AG3
            FF = G(KP1, 2) - AF3
            G(K, 1) = G(KP1, 1) - (AG1 + 2.d0*(AG2 + AG3) + U*FF - X*GG)/6.d0
            G(K, 2) = G(KP1, 2) - (AF1 + 2.d0*(AF2 + AF3) + X*FF - Y*GG)/6.d0
            if (G(K, 1)*G(KP1, 1) < 0d0) then
               NN = NN + 1
            end if
            AG1 = U*G(K, 2) - X*G(K, 1)
            AF1 = X*G(K, 2) - Y*G(K, 1)
            if (K == K2) goto 1000
            D(1, I) = AG1
            D(2, I) = AF1
         end do
         ! ------ ALL REMAINING POINTS -----
         Q = 1.d0/EA
         DG1 = D(1, 1)
         DG2 = D(1, 2)
         DG3 = D(1, 3)
         DF1 = D(2, 1)
         DF2 = D(2, 2)
         DF3 = D(2, 3)
         do
            KP1 = K
            K = K - 1
            RPB = RPB*Q
            DR = A*RPB
            R = RPB - B
            PHI = (E + ZZ/R - V(K))*DR/C
            U = DR*C + PHI
            X = -DR/R
            Y = -FLLP1*X*X/U + PHI
            DET = R83SQ - X*X + U*Y
            B1 = G(KP1, 1)*H83 + R1*DG1 + R2*DG2 + R3*DG3
            B2 = G(KP1, 2)*H83 + R1*DF1 + R2*DF2 + R3*DF3
            G(K, 1) = (B1*(H83 - X) + B2*U)/DET
            G(K, 2) = (B2*(H83 + X) - B1*Y)/DET
            if (G(K, 1)*G(KP1, 1) < 0.d0) then
               NN = NN + 1
            end if
            DG1 = DG2
            DF1 = DF2
            DG2 = DG3
            DF2 = DF3
            DG3 = U*G(K, 2) - X*G(K, 1)
            DF3 = X*G(K, 2) - Y*G(K, 1)
            if (MOD(K, 2) /= 0) then
               if (K <= K2 .or. G(K, 1)*DG3 >= 0.d0) exit
            end if
         end do
      end if
      ! -----------------------
1000  KC = K
      VAL = G(KC, 1)
      SLO = DG3/(A*(rofi(KC) + B))
   end subroutine RSQSR2

   subroutine POISS0(Z, A, B, rofi, RHO, NR, VHRMAX, V, RHOVH, VSUM, NSP)
      !  HARTREE POTENTIAL FOR SPHERICAL RHO. VHRMAX=VALUE AT RMAX.
      !  RETURNS VSUM = INTEGRAL OVER THAT POT WHICH IS ZERO AT RMAX.
      !  RHO = ´SPHERICAL CHDEN´ = 4PI*R*R*RHOTRUE. V=VTRUE.
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      integer, intent(in) :: NR, NSP
      real(rp), intent(in) :: A, B, VHRMAX, Z
      real(rp), intent(out) :: VSUM
      real(rp), dimension(NR), intent(inout) :: rofi
      real(rp), dimension(NSP), intent(inout) :: RHOVH
      real(rp), dimension(NR, NSP), intent(in) :: RHO
      real(rp), dimension(NR, NSP), intent(inout) :: V
      !
      !.. Local Scalars ..
      integer :: IR, ISP
      real(rp) :: A2B4, BB, CC, DD, DF, DRDI, EA, F, F2, F3, F4, G, PI, R, &
                  R2, R3, R4, RMAX, RO, RPB, SRDRDI, VHAT0, VNOW, WGT, &
                  X23, X34, Y2, Y3, Y4
      !
      !.. Intrinsic Functions ..
      intrinsic ATAN, EXP, MOD, SQRT
      !
      ! ... Executable Statements ...
      !
      PI = 4.d0*ATAN(1.d0)
      EA = EXP(A)
      RPB = B
      do IR = 1, NR
         rofi(IR) = RPB - B
         RPB = RPB*EA
      end do
      RMAX = rofi(NR)
      ! ----- APPROXIMATE RHO/R**2 BY CC+BB*R+DD*R*R NEAR ZERO ---
      R2 = rofi(2)
      R3 = rofi(3)
      R4 = rofi(4)
      F2 = 0.d0
      F3 = 0.d0
      F4 = 0.d0
      do ISP = 1, NSP
         F2 = F2 + RHO(2, ISP)/R2**2
         F3 = F3 + RHO(3, ISP)/R3**2
         F4 = F4 + RHO(4, ISP)/R4**2
      end do
      X23 = (R3*R3*F2 - R2*R2*F3)/(R3 - R2)
      X34 = (R4*R4*F3 - R3*R3*F4)/(R4 - R3)
      CC = (R2*X34 - R4*X23)/(R3*(R2 - R4))
      BB = ((R2 + R3)*X34 - (R3 + R4)*X23)/(R3*R3*(R4 - R2))
      DD = (F2 - BB*R2 - CC)/R2**2
      ! ------ NUMEROV FOR INHOM SOLUTION --------
      A2B4 = A*A/4.d0
      V(1, 1) = 1.d0
      DF = 0.d0
      do IR = 2, 3
         R = rofi(IR)
         DRDI = A*(R + B)
         SRDRDI = SQRT(DRDI)
         V(IR, 1) = V(1, 1) - R*R*(CC/3.d0 + R*BB/6.d0 + R*R*DD/10.d0)
         G = V(IR, 1)*R/SRDRDI
         F = G*(1.d0 - A2B4/12.d0)
         if (IR == 2) then
            Y2 = -2.d0*F2*R2*DRDI*SRDRDI
         end if
         if (IR == 3) then
            Y3 = -2.d0*F3*R3*DRDI*SRDRDI
         end if
         DF = F - DF
      end do
      IR = 3
      do
         IR = IR + 1
         R = rofi(IR)
         DRDI = A*(R + B)
         SRDRDI = SQRT(DRDI)
         RO = 0.d0
         do ISP = 1, NSP
            RO = RO + RHO(IR, ISP)
         end do
         Y4 = -2.d0*DRDI*SRDRDI*RO/R
         DF = DF + G*A2B4 + (Y4 + 10.d0*Y3 + Y2)/12.d0
         F = F + DF
         G = F/(1.d0 - A2B4/12.d0)
         V(IR, 1) = G*SRDRDI/R
         Y2 = Y3
         Y3 = Y4
         if (IR >= NR) exit
      end do
      ! ------ ADD CONSTANT TO GET V(NR)=VHRMAX -------
      VNOW = V(NR, 1) - 2.d0*Z/RMAX
      do IR = 1, NR
         V(IR, 1) = V(IR, 1) + (VHRMAX - VNOW)
      end do
      do ISP = 1, NSP
         RHOVH(ISP) = 0.d0
      end do
      VSUM = 0.d0
      VHAT0 = 0.d0
      do IR = 2, NR
         R = rofi(IR)
         DRDI = A*(R + B)
         WGT = 2*(MOD(IR + 1, 2) + 1)/3.d0
         if (IR == NR) then
            WGT = 1.d0/3.d0
         end if
         RO = 0.d0
         do ISP = 1, NSP
            RHOVH(ISP) = RHOVH(ISP) + WGT*DRDI*RHO(IR, ISP)*(V(IR, 1) - 2.d0*Z/R)
            RO = RO + RHO(IR, ISP)
         end do
         VHAT0 = VHAT0 + WGT*DRDI*RO*(1.d0/R - 1.d0/RMAX)
         VSUM = VSUM + WGT*DRDI*R*R*(V(IR, 1) - VHRMAX)
      end do
      VSUM = 4.d0*PI*(VSUM - Z*RMAX*RMAX)
      VHAT0 = 2.d0*VHAT0 + 2.d0*Z/RMAX + VHRMAX
      V(1, 1) = VHAT0
      if (NSP /= 1) then
         do IR = 1, NR
            V(IR, 2) = V(IR, 1)
         end do
      end if
   end subroutine POISS0

   subroutine VXC0SP(this, xc_obj, Z, A, B, rofi, RHO, NR, V, RHO0, RHOEPS, RHOMU, NSP, B_fsm)
      !  ADDS XC PART TO SPHERICAL POTENTIAL, MAKES INTEGRALS RHOMU AND RHOEP
      !
      ! use xcdata
      ! use common_asd, only : asd_atom, bxc
      ! use control_mod, only : control_obj
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      class(self), intent(inout) :: this
      type(xc), intent(in) :: xc_obj
      integer, intent(in) :: NR, NSP
      real(rp), intent(in) :: A, B
      real(rp) :: Z
      real(rp), dimension(2), intent(inout) :: RHO0
      real(rp), dimension(NR), intent(in) :: rofi
      real(rp), dimension(NSP), intent(inout) :: RHOEPS, RHOMU
      real(rp), dimension(NR, NSP), intent(in) :: RHO
      real(rp), dimension(NR, NSP), intent(inout) :: V
      real(rp), intent(in)  :: B_fsm
      !
      !.. Local Scalars ..
      integer :: IR, ISP, IXC
      real(rp) :: DRDI, EXC, EXC1, EXC2, OB4PI, PI, RHO2, RHO3, RHOT1, &
                  RHOT2, RHOTRU, VXC, VXC1, VXC2, WGT, &
                  RHO1, R, RCE
      real(rp), dimension(NR, NSP) :: RHOP, RHOPP, tRHO
      real(rp), dimension(2) :: RHOD, RHODD
      real(rp) :: Bxc_up, Bxc_dw, Bxc_tot
      !.. External Calls ..
      ! external EVXC
      !
      !.. Intrinsic Functions ..
      intrinsic ATAN, MOD
      !
      ! ... Executable Statements ...
      !
      PI = 4.d0*ATAN(1.d0)
      OB4PI = 1.d0/(4.d0*PI)
      IXC = xc_obj%txc
      !
      ! Constraining field related hacks below
      Bxc_up = 0.0d0
      Bxc_dw = 0.0d0
      !
      ! Extrapolate density to core point
      do ISP = 1, NSP
         RHOEPS(ISP) = 0.d0
         RHOMU(ISP) = 0.d0
         RHO2 = RHO(2, ISP)/rofi(2)**2
         RHO3 = RHO(3, ISP)/rofi(3)**2
         RHO0(ISP) = OB4PI*(RHO2*rofi(3) - RHO3*rofi(2))/(rofi(3) - rofi(2))
      end do
      ! Calculate gradients in case of GGA
      do ISP = 1, NSP
         tRHO(1, ISP) = RHO0(ISP)
         do IR = 2, NR
            tRHO(IR, 1) = RHO(IR, 1)*OB4PI/rofi(IR)**2
            tRHO(IR, ISP) = RHO(IR, ISP)*OB4PI/rofi(IR)**2
         end do
         if (IXC == 5 .or. IXC >= 8) then
            !subroutine radgra(a, b, nr, rofi, f, gradf)
            !call DIFFN(tRHO(1, ISP), RHOP(1, ISP), RHOPP(1, ISP), NR, A)
            call radgra(a, b, NR, rofi, tRHO(1, ISP), RHOP(1, ISP))
            call radgra(a, b, NR, rofi, RHOP(1, ISP), RHOPP(1, ISP))
            !Comment line below or not?
            !call DIFFN(tRHO(2, ISP), RHOP(2, ISP), RHOPP(2, ISP), NR, A)
            !call radgra(a, b, NR, rofi, tRHO(2, ISP), RHOP(2, ISP))
            !call radgra(a, b, NR, rofi, RHOP(2, ISP), RHOPP(2, ISP))
         else
            RHOP = 0.0d0
            RHOPP = 0.0d0
         end if
      end do
      if (NSP == 1) then
         RHO1 = 0.5*tRHO(1, 1)
         RHO2 = RHO1
         R = (rofi(3) - rofi(2))
         RCE = R*R
         if (IXC == 5 .or. IXC >= 8) then
            ! RHODD(2) = RHODD(1)
            RHOD(1) = 0.5*RHOP(1, 1)/R
            RHODD(1) = 0.5*(RHOPP(1, 1) - RHOP(1, 1))/RCE
            RHOD(2) = RHOD(1)
            RHODD(2) = RHODD(1)
         end if
         ! CALL EVXC(RHO0(1), 0.5D0*RHO0(1), EXC, VXC)
         call xc_obj%XCPOT(RHO2, RHO1, tRHO(1, 1), RHOD, RHODD, R, VXC2, VXC1, EXC)
         V(1, 1) = V(1, 1) + VXC1
         do IR = 2, NR
            RHO1 = 0.5*tRHO(1, 1)
            RHO2 = RHO1
            R = (rofi(3) - rofi(2))
            RCE = R*R
            if (IXC == 5 .or. IXC >= 8) then
               RHOD(1) = 0.5*RHOP(IR, 1)/R
               RHODD(1) = 0.5*(RHOPP(IR, 1) - RHOP(IR, 1))/RCE
               RHOD(2) = RHOD(1)
               RHODD(2) = RHODD(1)
            end if
            call xc_obj%XCPOT(RHO2, RHO1, RHO(1, 1), RHOD, RHODD, R, VXC2, VXC1, EXC1)
            !RHOTRU = RHO(IR, 1) * OB4PI / rofi(IR)**2
        !!       Call EVXC(RHOTRU, 0.5D0*RHOTRU, EXC, VXC)
            !call EVXC(RHOTRU, 0.5*RHOTRU, EXC, VXC)
            V(IR, 1) = V(IR, 1) + VXC1
            WGT = 2*(MOD(IR + 1, 2) + 1)/3.d0
            if (IR == 1 .or. IR == NR) then
               WGT = 1.d0/3.d0
            end if
            DRDI = A*(rofi(IR) + B)
            RHOEPS(1) = RHOEPS(1) + WGT*DRDI*RHO(IR, 1)*EXC1
            RHOMU(1) = RHOMU(1) + WGT*DRDI*RHO(IR, 1)*VXC1
         end do
      else
         ! call EVXC(RHO0(1)+RHO0(2), RHO0(1), EXC1, VXC1)
         ! call EVXC(RHO0(1)+RHO0(2), RHO0(2), EXC2, VXC2)
         RHO1 = tRHO(1, 1)
         RHO2 = tRHO(1, 2)
         R = (rofi(3) - rofi(2))
         RCE = R*R
         if (IXC == 5 .or. IXC >= 8) then
            ! RHOD(2) = RHOP(IR, 1) / R
            ! RHOD(1) = RHOP(IR, 2) / R
            ! RHODD(2) = (RHOPP(IR, 1)-RHOP(IR, 1)) / RCE
            ! RHODD(1) = (RHOPP(IR, 2)-RHOP(IR, 2)) / RCE
            RHOD(2) = RHOP(1, 1)/R
            RHOD(1) = RHOP(1, 2)/R
            RHODD(2) = (RHOPP(1, 1) - RHOP(1, 1))/RCE
            RHODD(1) = (RHOPP(1, 2) - RHOP(1, 2))/RCE
         end if
         call xc_obj%XCPOT(RHO2, RHO1, RHO(1, 1), RHOD, RHODD, R, VXC2, VXC1, EXC1)
         !V(1, 1) = V(1, 1) + VXC1
         !V(1, 2) = V(1, 2) + VXC2
         V(1, 1) = V(1, 1) + VXC1 + B_fsm
         V(1, 2) = V(1, 2) + VXC2 - B_fsm
         do IR = 2, NR
            R = rofi(IR)
            RCE = R*R
            RHO3 = tRHO(IR, 1) + tRHO(IR, 2)
            RHO1 = tRHO(IR, 1)
            RHO2 = tRHO(IR, 2)
            if (IXC == 5 .or. IXC >= 8) then
               RHOD(2) = RHOP(IR, 1)/R
               RHOD(1) = RHOP(IR, 2)/R
               RHODD(2) = (RHOPP(IR, 1) - RHOP(IR, 1))/RCE
               RHODD(1) = (RHOPP(IR, 2) - RHOP(IR, 2))/RCE
            end if
            call xc_obj%XCPOT(RHO2, RHO1, RHO3, RHOD, RHODD, R, VXC2, VXC1, EXC1)
            !
            !RHOT1 = RHO(IR, 1) * (OB4PI/rofi(IR)**2)
            !RHOT2 = RHO(IR, 2) * (OB4PI/rofi(IR)**2)
            !call EVXC(RHOT1+RHOT2, RHOT1, EXC1, VXC1)
            !call EVXC(RHOT1+RHOT2, RHOT2, EXC2, VXC2)
            !V(IR, 1) = V(IR, 1) + VXC1
            !V(IR, 2) = V(IR, 2) + VXC2
            ! Insertion of constraining field for FSM
            V(IR, 1) = V(IR, 1) + VXC1 + B_fsm
            V(IR, 2) = V(IR, 2) + VXC2 - B_fsm
            WGT = 2*(MOD(IR + 1, 2) + 1)/3.d0
            if (IR == 1 .or. IR == NR) then
               WGT = 1.d0/3.d0
            end if
            DRDI = A*(rofi(IR) + B)
            RHOEPS(1) = RHOEPS(1) + WGT*DRDI*RHO(IR, 1)*EXC1
            RHOMU(1) = RHOMU(1) + WGT*DRDI*RHO(IR, 1)*(VXC1 + B_fsm)
            !RHOMU(1) = RHOMU(1) + WGT*DRDI*RHO(IR, 1)*VXC1
            RHOEPS(2) = RHOEPS(2) + WGT*DRDI*RHO(IR, 2)*EXC1
            RHOMU(2) = RHOMU(2) + WGT*DRDI*RHO(IR, 2)*(VXC2 - B_fsm)
            !RHOMU(2) = RHOMU(2) + WGT*DRDI*RHO(IR, 2)*VXC2
            !Bxc_up=Bxc_up + WGT*DRDI*VXC1
            !Bxc_dw=Bxc_dw + WGT*DRDI*VXC2
            Bxc_up = Bxc_up + WGT*DRDI*VXC1
            Bxc_dw = Bxc_dw + WGT*DRDI*VXC2
         end do
         !write(*, *)wgt, drdi, exc1, vxc1, vxc2
         Bxc_tot = Bxc_up - Bxc_dw
         ! RHOEPS=RHOEPS/4.0d0
         ! AB comment
      !!! if(this%lattice%control%do_asd) this%bxc(this%lattice%control%asd_atom)=Bxc_tot
         !print *, ´B_XC´, Bxc_tot, rhomu(1)-rhomu(2)
         !print *, ´b_XC´, 235e3*Bxc_tot, 235e3*(rhomu(1)-rhomu(2))
      end if
   end subroutine VXC0SP

   subroutine radgra(a, b, nr, rofi, f, gradf)
      !- radial gradient
      !-----------------------------------------------------------------------
      !i Inputs:
      !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
      !i   b     :                 -//-
      !i   nr    :number of mesh points
      !i   rofi  :radial mesh points
      !i   f     :given function defined in a mesh rofi(i)
      !o Outputs:
      !o  gradf  :derivative of the function f defined in a mesh rofi(i)
      !-----------------------------------------------------------------------
      implicit none
      ! Passed variables:
      integer, intent(in) ::  nr
      real(rp), intent(in) :: a, b
      real(rp), dimension(nr), intent(in) :: rofi
      real(rp), dimension(nr), intent(in) :: f
      real(rp), dimension(nr), intent(out) :: gradf

      ! Local variables:
      integer :: nm2, i
      !
      nm2 = nr - 2
      ! for the first and second point are used the Forward edifferences
      ! (Handbook, 25.3.9 with 25.1.1)
      gradf(1) = ((6.d0*f(2) + 20.d0/3.d0*f(4) + 1.2d0*f(6)) &
                  - (2.45d0*f(1) + 7.5d0*f(3) + 3.75d0*f(5) + 1.d0/6.d0*f(7)))/a
      !
      gradf(2) = ((6.d0*f(3) + 20.d0/3.d0*f(5) + 1.2d0*f(7)) &
                  - (2.45d0*f(2) + 7.5d0*f(4) + 3.75d0*f(6) + 1.d0/6.d0*f(8)))/a
      !
      ! --- Five points´ formula  (25.3.6)
      do i = 3, nm2
         gradf(i) = ((f(i - 2) + 8.d0*f(i + 1)) - (8.d0*f(i - 1) + f(i + 2)))/12.d0/a
      end do
      ! --- Five points´ formula  (25.3.6)
      gradf(nr - 1) = (-1.d0/12.d0*f(nr - 4) + 0.5d0*f(nr - 3) - 1.5d0*f(nr - 2) &
                       + 5.d0/6.d0*f(nr - 1) + 0.25d0*f(nr))/a
      gradf(nr) = (0.25d0*f(nr - 4) - 4.d0/3.d0*f(nr - 3) + 3.d0*f(nr - 2) &
                   - 4.d0*f(nr - 1) + 25.d0/12.d0*f(nr))/a
      !
      ! --- Three points´ formula  (25.3.4)
      !     gradf(nr-1)=(f(nr)-f(nr-2))/2.d0/a
      !     gradf(nr)=(f(nr-2)/2.d0-2.d0*f(nr-1)+1.5d0*f(nr))/a
      do i = 1, nr
         gradf(i) = gradf(i)/(rofi(i) + b)
      end do
      return
   end subroutine radgra

   subroutine RACSI(this, atom, ROFI, QSL)
      ! use common_pri
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      class(self), intent(inout) :: this
      type(symbolic_atom), intent(inout) :: atom
      real(rp), dimension(6), intent(inout) :: QSL
      real(rp), dimension(:), intent(in) :: ROFI
      !real(rp), dimension(500), intent(in) :: ROFI
      !
      !.. Local Scalars ..
      integer :: IDW, II, INUM, INUM1, IR, IR1, ISI, IUP, JM, JP, NSP
      real(rp) :: C2, DRDI, FAK2, FAK4, S12, SUM, SUM1, SUM2, WGT
      !
      !.. Local Arrays ..
      real(rp), dimension(:, :), allocatable :: DVDR, DVM, DVP
      real(rp) :: B
      integer :: NR
      !
      !.. Intrinsic Functions ..
      intrinsic MOD
      !
      ! ... Executable Statements ...
      !
      ! open (15, FILE = "dracsi")
      NR = size(ROFI)
      B = atom%B()
      NSP = 2
      C2 = 274.074d0**2
      allocate (DVDR(NR, 2), DVM(NR, 2), DVP(NR, 2))
      !     interpolate the derivatives
      do ISI = 1, 2
         do II = 3, NR - 1
            JP = II + 1
            JM = II - 1
            DVP(II, ISI) = ((this%VZT(JP, ISI) - this%VZT(II, ISI))/(ROFI(JP) - ROFI(II)))
            DVM(II, ISI) = ((this%VZT(JM, ISI) - this%VZT(II, ISI))/(ROFI(JM) - ROFI(II)))
            DVDR(II, ISI) = (DVP(II, ISI) + DVM(II, ISI))/2.0d0
         end do
         DVDR(2, ISI) = DVDR(3, ISI)
         DVDR(NR, ISI) = DVDR(NR - 1, ISI)
      end do
      !-- calculates QSRD(IT, 1) and QSRD(IT, 2), prefactors qssi for L.S (SPIN UP)
      !-- calculates QSRD(IT, 4) and QSRD(IT, 5), prefactors qssi for L.S (SPIN DW)
      !---QSRD(IT, 1or4) for p-band and QSRD(IT, 2or5) for d-band-------------------
      do INUM = 2, 3
         do ISI = 1, NSP
            SUM = 0.0d0
            do IR = 2, NR
               WGT = 2*(MOD(IR + 1, 2) + 1)/3.0d0
               if (IR == 1 .or. IR == NR) then
                  WGT = 1.0d0/3.0d0
               end if
               DRDI = atom%A*(ROFI(IR) + B)
               SUM = SUM + &
                     WGT*DRDI*this%FUN2(IR, INUM, ISI)*2.0d0*DVDR(IR, ISI)/(ROFI(IR)*C2)
               IUP = INUM - 1
               IDW = INUM + 2
            end do
            if (ISI == 1) then
               QSL(IUP) = SUM
            else
               QSL(IDW) = SUM
            end if
         end do
      end do
      !----CALCULATES RACAH COEF QSR(IT, 3)-(spin up) and QSRD(IT, 4)-(spin-down)
      do ISI = 1, NSP
         do INUM = 2, 4, 2
            INUM1 = INUM + 1
            SUM = 0.0d0
            do IR = 2, NR
               SUM1 = 0.0d0
               !---------first part of integral
               do IR1 = 2, IR
                  WGT = 2*(MOD(IR1 + 1, 2) + 1)/3.0d0
                  if (IR1 == 1 .or. IR1 == IR) then
                     WGT = 1.0d0/3.0d0
                  end if
                  DRDI = atom%A*(ROFI(IR1) + B)
                  SUM1 = SUM1 + WGT*DRDI*this%FUN2(IR1, 3, ISI)*(ROFI(IR1)**INUM)/(ROFI(IR)**INUM1)
               end do
               !-----------second part of integral
               SUM2 = 0.0d0
               do IR1 = IR, NR
                  WGT = 2*(MOD(IR1 + 1, 2) + 1)/3.0d0
                  if (IR1 == IR .or. IR1 == NR) then
                     WGT = 1.0d0/3.0d0
                  end if
                  DRDI = atom%A*(ROFI(IR1) + B)
                  SUM2 = SUM2 + WGT*DRDI*this%FUN2(IR1, 3, ISI)*(ROFI(IR)**INUM)/(ROFI(IR1)**INUM1)
               end do
               !-----------------------------------------
               S12 = SUM1 + SUM2
               WGT = 2*(MOD(IR + 1, 2) + 1)/3.0d0
               if (IR == 1 .or. IR == NR) then
                  WGT = 1.0d0/3.0d0
               end if
               DRDI = atom%A*(ROFI(IR) + B)
               SUM = SUM + WGT*DRDI*S12*this%FUN2(IR, 3, ISI)
            end do
            if (INUM == 2) then
               FAK2 = SUM/49.d0
               FAK4 = 0.0d0
            else
               FAK4 = SUM/441.d0
            end if
            if (ISI == 1) then
               QSL(3) = 2.d0*(FAK2 - 5*FAK4)
            else
               QSL(6) = 2.d0*(FAK2 - 5*FAK4)
            end if
         end do
         atom%potential%xi_p(:) = [qsl(1), qsl(4)]
         atom%potential%xi_d(:) = [qsl(2), qsl(5)]
         atom%potential%rac(:) = [qsl(3), qsl(6)]
         ! WRITE(17, 56)FAK2, FAK4
         ! WRITE(17, *)RCH
      end do
      !write (15, 10000) QSL(1), QSL(2), QSL(3)
      !write (15, 10000) QSL(4), QSL(5), QSL(6)
      return
      !
      ! ... Format Declarations ...
      !
10000 format(3f13.5)
      !     DO IR=2, 322, 5
      !     WRITE(1, 55)IR, DIF(IR), DIF(IR+1), DIF(IR+2), DIF(IR+3), DIF(IR+4)
      !     END DO
10001 format(i5, 5e13.5)
10002 format(i5, 1pe15.6, 1pe15.6, 1pe15.6, 1pe15.6, 1pe15.6)
10003 format( &
         8x, "  ROFI  ", "       V(R)        ", "     FUN2(S)     ", &
         "    FUN2(P)    ", "    FUN2(D)    ")
10004 format(i5, 2f15.6)
   end subroutine RACSI

   subroutine POTPAR(this, atom, V, ROFI)
      !  MAKES POTENTIAL PARAMETERS FOR POTENTIAL V. BOUNDARY CONDITION
      !  FOR PHI IS PNU=.5-ATAN(DNU)/PI+(PRINC.QUANT.NUMBER).
      !  C=CENTER OF BAND=ENU+OMEGA-, QPAR=1/Q, PPAR=1/SQRT(P),
      !  SRDEL=SQRT(DELTA) BUT WITH PROPER SIGN (THAT OF PHI-).
      !  THIS IS THE SCALAR-RELATIVISTIC VERSION
      !  ++++++ PHDTSR REPLACED BY PHDFSR  9.12.87
      !  ++++++ TOL SET SMALLER BECAUSE OF NUM DIFFERENTIATION 24.3.88
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      class(self) :: this
      type(symbolic_atom) :: atom
      real(rp), dimension(:), intent(in) :: ROFI
      real(rp), dimension(:, :), intent(in) :: V
      !
      !.. Local Scalars ..
      integer :: I, L, LP1, NN, NRE, LMAX
      real(rp) :: B, DLPHI, DLPHIP, DPHI, DPHIP, E, EB1, EB2, OMEGAM, OMEGAP, P, PHI, PHIP, PHMINS, PHPLUS, PI, Q, SLO, SUMM, TOL, VAL, RMAX
      !
      !.. Local Arrays ..
      integer, dimension(0:9) :: KONFIG
      real(rp), dimension(0:10) :: DNU
      !
      !.. External Calls ..
      ! external PHDFSR, RSEQSR
      !
      !.. Intrinsic Functions ..
      intrinsic ATAN, EXP, SQRT, TAN
      !
      real(rp), dimension(:, :), allocatable :: G, GP, GPP
      integer :: NSP, NR
      ! ... Executable Statements ...
      !
      TOL = 1.d-12
      EB1 = -10.d0
      EB2 = 10.d0
      PI = 4.d0*ATAN(1.d0)
      B = atom%B()
      nsp = 2
      rmax = atom%potential%ws_r
      nr = size(ROFI)
      lmax = atom%potential%lmax
      allocate (G(NR, NSP), GP(NR, NSP), GPP(NR, NSP))
      do I = 1, nsp
         do L = 0, lmax
            KONFIG(L) = int(atom%potential%PNU(L, I))
            DNU(L) = TAN(PI*(0.5d0 - atom%potential%PNU(L, I)))
         end do
         !   write (6, 10000) (PNU(L, I), L = 0, LMAX)
         !   write (6, 10001) (DNU(L), L = 0, LMAX)
         !   write (6, 10003)
         do L = 0, LMAX
            LP1 = L + 1
            NN = KONFIG(L) - L - 1
            E = -0.5d0
            VAL = RMAX
            SLO = DNU(L) + 1.d0
            call RSEQSR(EB1, EB2, E, TOL, atom%element%atomic_number, L, NN, VAL, SLO, V(:, I), G, SUMM, atom%A, B, ROFI, NR, NRE, 0)
            VAL = VAL/SQRT(SUMM)
            SLO = SLO/SQRT(SUMM)
            call PHDFSR(atom%element%atomic_number, L, V(:, I), E, atom%A, B, ROFI, NR, G, VAL, SLO, GP, GPP, PHI, DPHI, PHIP, DPHIP, P, TOL, NN)
        !!write(876, ´(f18.8)´) G(1:NR*2)
            !print ´(a, 6f12.6)´, ´PHI´, PHI, DPHI, PHIP, DPHIP, val
            !print ´(1x, a, 2i5, 3f12.6)´, ´ENU old, new´, L, I, ENU(L, I), E, RMAX
            atom%potential%ENU(L, I) = E
            DLPHI = RMAX*DPHI/PHI
            DLPHIP = RMAX*DPHIP/PHIP
            OMEGAM = -(PHI/PHIP)*(-L - 1 - DLPHI)/(-L - 1 - DLPHIP)
            OMEGAP = -(PHI/PHIP)*(L - DLPHI)/(L - DLPHIP)
            PHPLUS = PHI + OMEGAP*PHIP
            PHMINS = PHI + OMEGAM*PHIP
            atom%potential%C(L, I) = E + OMEGAM
            ! C = C for canonical scaling
            ! omegam= -wk/wkdot
            atom%potential%VL(L, I) = E + OMEGAP
            ! Vl =V for canonical scaling
            atom%potential%SRDEL(L, I) = PHMINS*SQRT(0.5d0*RMAX)
            !if(L==0) SRDEL(L, I)=-1.0d0*SRDEL(L, I)
            !print ´(2x, a, 8f12.6)´, ´potpars:´, E, omegam, rmax*phi*phi, phmins/phplus, 1.0d0/abs(phip), E+OMEGAM, E+omegap
            !print *, ´rmax´, rmax
            !print *, ´1/phmi´, 1.0d0/phmins
            !print *, ´srdel´, SRDEL(L, I)
        !! phmins= phi + omegam*phip = phi - (phi/phip)*(-L-1-DLPHI)/(-L-1-DLPHIP)
        !! phmins= phi + omegam*phip = phi - (phi/phip)*(-L-1-DLPHI)/-L-1-DLPHIP)
        !! wkdot= ((dphidt-l)*gi(l) +    (l+l+1)*gi(l+1))*phidt
            Q = PHMINS/(2*(2*L + 1)*PHPLUS)
            !Q=gamma
            atom%potential%QPAR(L, I) = 1.d0/Q
            atom%potential%PPAR(L, I) = 1.d0/SQRT(P)
            !      write (660, 10002) L, E, VL(L, I), C(L, I), SRDEL(L, I), QPAR(L, I), PPAR(L, I)
         end do
      end do
      return
      !
      ! ... Format Declarations ...
      !
      ! 10000 format (/" POTPAR:  INPUT LOG DERIVATIVES ARE"/" PNU=", 6f12.6)
      ! 10001 format (" DNU=", 6f12.6)
      ! 10002 format (i3, 4f12.5, 2f12.4)
      ! 10003 format (/"  L", 7x, "ENU", 9x, "V", 11x, "C", 10x, "SRDEL", 8x, "1/Q", 9x, "1/SRP")
   end subroutine POTPAR

end module self_mod
