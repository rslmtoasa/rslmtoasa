!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Calculation
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
!> Module to handle pre-processing, processing and post-processing calculations
!------------------------------------------------------------------------------

module calculation_mod

   use mpi_mod
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
   use spin_dynamics_mod
   use mix_mod
   use math_mod
   use precision_mod, only: rp
   use string_mod, only: sl, fmt, real2str
   use timer_mod, only: g_timer
   use logger_mod, only: g_logger
   implicit none

   private

   type, public :: calculation
      !> Pre-processing. Options are:
      !> 'none' (default)
      !> 'bravais' : Builds the bulk clust
      !> 'buildsurf' : Builds the surface clust
      !> 'newclubulk' : Builds the imputiry clust from the bluk clust
      !> 'newclusurf' : Builds the impurity clust from the surface clust
      character(len=sl) :: pre_processing

      !> Processing. Options are
      !> 'none' (default)
      character(len=sl) :: processing

      !> Post-processing. Options are
      !> 'none' (default)
      character(len=sl) :: post_processing

      !> Controller for preprocessing verbosity.
      !>
      !> Controller for preprocessing verbosity. If true, call the subroutines:
      !> @ref control.print_state, @ref lattice.print_state,
      !> @ref self.print_state and @ref charge.print_state after pre-processing calls.
      logical :: verbose

      !> name list input file
      character(len=sl) :: fname
   contains
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure, private :: pre_processing_bravais
      procedure, private :: pre_processing_buildsurf
      procedure, private :: pre_processing_newclubulk
      procedure, private :: pre_processing_newclusurf
      procedure, private :: processing_sd
      procedure, private :: post_processing_paoflow2rs
      procedure, private :: post_processing_exchange
      procedure, private :: post_processing_exchange_p2rs
      procedure :: process
      final :: destructor
   end type calculation

   interface calculation
      procedure :: constructor
   end interface calculation

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] fname Namelist file
   !> @return type(calculation)
   !---------------------------------------------------------------------------
   function constructor(fname) result(obj)
      type(calculation) :: obj
      character(len=*), intent(in) :: fname
      call obj%restore_to_default()
      call obj%build_from_file(fname)
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   subroutine destructor(this)
      type(calculation) :: this
   end subroutine destructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Read parameters from input file
   !
   !> @param[in] fname Namelist file with '&calculation' table
   !---------------------------------------------------------------------------
   subroutine build_from_file(this, fname)
      class(calculation), intent(inout) :: this
      character(len=*), intent(in) :: fname
      ! variables associated with the reading processes
      integer :: iostatus, funit

      include 'include_codes/namelists/calculation.f90'

      verbose = this%verbose
      pre_processing = this%pre_processing
      processing = this%processing
      post_processing = this%post_processing

      open (newunit=funit, file=fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//fmt('A', fname)//' not found', __FILE__, __LINE__)
      end if

      read (funit, nml=calculation, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%error('Error while reading namelist', __FILE__, __LINE__)
         call g_logger%error('iostatus = '//fmt('I0', iostatus), __FILE__, __LINE__)
      end if

      ! Pre-processing
      call check_pre_processing(trim(pre_processing))
      ! Processing
      call check_processing(trim(processing))
      ! Post-processing
      call check_post_processing(trim(post_processing))

      this%verbose = verbose
      this%fname = fname
      this%pre_processing = pre_processing
      this%processing = processing
      this%post_processing = post_processing

      close (funit)
   end subroutine build_from_file

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Handle for general process
   !---------------------------------------------------------------------------
   subroutine process(this)
      class(calculation), intent(in) :: this

      ! Pre-processing
      select case (this%pre_processing)
      case ('bravais')
         call this%pre_processing_bravais()
      case ('buildsurf')
         call this%pre_processing_buildsurf()
      case ('newclubulk')
         call this%pre_processing_newclubulk()
      case ('newclusurf')
         call this%pre_processing_newclusurf()
      end select

      ! Processing
      select case (this%processing)
      case ('sd')
         call this%processing_sd()
      end select

      ! Post-processing
      select case (this%post_processing)
      case ('paoflow2rs')
         call this%post_processing_paoflow2rs()
      case ('exchange')
         call this%post_processing_exchange()
      case ('exchange_p2rs')
         call this%post_processing_exchange_p2rs()
      end select
   end subroutine

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Process for the spin dynamics calculation
   !---------------------------------------------------------------------------
   subroutine processing_sd(this)
      class(calculation), intent(in) :: this

      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(energy), target :: energy_obj
      type(self), target :: self_obj
      type(charge), target :: charge_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(recursion), target :: recursion_obj
      type(green), target :: green_obj
      type(dos), target :: dos_obj
      type(bands), target :: bands_obj
      type(mix), target :: mix_obj
      type(spin_dynamics), target :: sd_obj

      ! Constructing control object
      control_obj = control(this%fname)

      ! Constructing lattice object
      lattice_obj = lattice(control_obj)

      ! Running the pre-calculation
      call g_timer%start('pre-processing')
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      call lattice_obj%newclu()
      call lattice_obj%structb()

      ! Creating the symbolic_atom object
      call lattice_obj%atomlist()

      ! Initializing MPI lookup tables and info.
      call get_mpi_variables(rank, lattice_obj%nrec)

      ! Constructing the charge object
      charge_obj = charge(lattice_obj)
      call charge_obj%impmad()
      call charge_obj%get_charge_transf
      call g_timer%stop('pre-processing')

      ! Constructing mixing object
      mix_obj = mix(lattice_obj, charge_obj)

      ! Creating the energy object
      energy_obj = energy(lattice_obj)

      ! Creating hamiltonian object
      hamiltonian_obj = hamiltonian(charge_obj)

      ! Creating recursion object
      recursion_obj = recursion(hamiltonian_obj, energy_obj)

      ! Creating density of states object
      dos_obj = dos(recursion_obj, energy_obj)

      ! Creating Green function object
      green_obj = green(dos_obj)

      ! Creating bands object
      bands_obj = bands(green_obj)

      ! Constructing mixing object
      mix_obj = mix(lattice_obj, charge_obj)

      ! Creating the energy object
      energy_obj = energy(lattice_obj)

      ! Creating hamiltonian object
      hamiltonian_obj = hamiltonian(charge_obj)

      ! Creating recursion object
      recursion_obj = recursion(hamiltonian_obj, energy_obj)

      ! Creating density of states object
      dos_obj = dos(recursion_obj, energy_obj)

      ! Creating Green function object
      green_obj = green(dos_obj)

      ! Creating bands object
      bands_obj = bands(green_obj)

      ! Creating the self object
      self_obj = self(bands_obj, mix_obj)

      ! Creating spin dynamics object
      sd_obj = spin_dynamics(self_obj)

      ! Spin dynamics routine
      call sd_obj%sd_run

      if (this%verbose) then
         call save_state(lattice_obj%symbolic_atoms)
      end if

   end subroutine processing_sd
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Pre-process for new clust bulk calculation
   !---------------------------------------------------------------------------
   subroutine pre_processing_newclubulk(this)
      class(calculation), intent(in) :: this

      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(energy), target :: energy_obj
      type(self), target :: self_obj
      type(charge), target :: charge_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(recursion), target :: recursion_obj
      type(green), target :: green_obj
      type(dos), target :: dos_obj
      type(bands), target :: bands_obj
      type(mix), target :: mix_obj

      ! Constructing control object
      control_obj = control(this%fname)

      ! Constructing lattice object
      lattice_obj = lattice(control_obj)

      ! Running the pre-calculation
      call g_timer%start('pre-processing')
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      call lattice_obj%newclu()
      call lattice_obj%structb()

      ! Creating the symbolic_atom object
      call lattice_obj%atomlist()

      ! Initializing MPI lookup tables and info.
      call get_mpi_variables(rank, lattice_obj%nrec)

      ! Constructing the charge object
      charge_obj = charge(lattice_obj)
      call charge_obj%impmad()
      call charge_obj%get_charge_transf
      call g_timer%stop('pre-processing')

      ! Constructing mixing object
      mix_obj = mix(lattice_obj, charge_obj)

      ! Creating the energy object
      energy_obj = energy(lattice_obj)

      ! Creating hamiltonian object
      hamiltonian_obj = hamiltonian(charge_obj)

      ! Creating recursion object
      recursion_obj = recursion(hamiltonian_obj, energy_obj)

      ! Creating density of states object
      dos_obj = dos(recursion_obj, energy_obj)

      ! Creating Green function object
      green_obj = green(dos_obj)

      ! Creating bands object
      bands_obj = bands(green_obj)

      ! Creating the self object
      self_obj = self(bands_obj, mix_obj)
      call g_timer%start('self-consistency')
      call self_obj%run()
      call g_timer%stop('self-consistency')

      if (this%verbose) then
         !call control_obj%print_state_formatted()
         !call lattice_obj%print_state_formatted()
         !call self_obj%print_state()
         !call charge_obj%print_state()
         call save_state(lattice_obj%symbolic_atoms)
      end if
   end subroutine pre_processing_newclubulk

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Pre-process for new clust surface calculation
   !---------------------------------------------------------------------------
   subroutine pre_processing_newclusurf(this)
      class(calculation), intent(in) :: this

      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(energy), target :: energy_obj
      type(self), target :: self_obj
      type(charge), target :: charge_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(recursion), target :: recursion_obj
      type(green), target :: green_obj
      type(dos), target :: dos_obj
      type(bands), target :: bands_obj
      type(mix), target :: mix_obj

      ! Constructing control object
      control_obj = control(this%fname)

      ! Constructing lattice object
      lattice_obj = lattice(control_obj)

      ! Running the pre-calculation
      call g_timer%start('pre-processing')
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      call lattice_obj%build_surf_full()
      call lattice_obj%newclu()
      call lattice_obj%structb()

      ! Creating the symbolic_atom object
      call lattice_obj%atomlist()

      ! Initializing MPI lookup tables and info.
      call get_mpi_variables(rank, lattice_obj%nrec)

      ! Constructing the charge object
      charge_obj = charge(lattice_obj)
      call charge_obj%impmad()
      !write(*,*) lattice_obj%nbulk
      call charge_obj%get_charge_transf
      call g_timer%stop('pre-processing')

      ! Constructing mixing object
      mix_obj = mix(lattice_obj, charge_obj)

      ! Creating the energy object
      energy_obj = energy(lattice_obj)

      ! Creating hamiltonian object
      hamiltonian_obj = hamiltonian(charge_obj)

      ! Creating recursion object
      recursion_obj = recursion(hamiltonian_obj, energy_obj)

      ! Creating density of states object
      dos_obj = dos(recursion_obj, energy_obj)

      ! Creating Green function object
      green_obj = green(dos_obj)

      ! Creating bands object
      bands_obj = bands(green_obj)

      ! Creating the self object
      self_obj = self(bands_obj, mix_obj)
      call g_timer%start('self-consistency')
      call self_obj%run()
      call g_timer%stop('self-consistency')

      if (this%verbose) then
         !call control_obj%print_state_formatted()
         !call lattice_obj%print_state_formatted()
         !call self_obj%print_state()
         !call charge_obj%print_state()
         call save_state(lattice_obj%symbolic_atoms)
      end if
   end subroutine pre_processing_newclusurf

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Pre-process for build surface calculation
   !---------------------------------------------------------------------------
   subroutine pre_processing_buildsurf(this)
      class(calculation), intent(in) :: this

      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(self), target :: self_obj
      type(energy), target :: energy_obj
      type(charge), target :: charge_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(recursion), target :: recursion_obj
      type(green), target :: green_obj
      type(dos), target :: dos_obj
      type(bands), target :: bands_obj
      type(mix), target :: mix_obj


      ! Constructing control object
      control_obj = control(this%fname)
      ! Constructing lattice object
      lattice_obj = lattice(control_obj)

      ! Running the pre-calculation
      call g_timer%start('pre-processing')
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      !call lattice_obj%build_clusup()
      call lattice_obj%build_surf_full()
      call lattice_obj%structb()

      ! Creating the symbolic_atom object
      call lattice_obj%atomlist()

      ! Initializing MPI lookup tables and info.
      call get_mpi_variables(rank, lattice_obj%nrec)

      ! Constructing the charge object
      charge_obj = charge(lattice_obj)
      call charge_obj%build_alelay
      call charge_obj%surfmat
      call g_timer%stop('pre-processing')
      ! Constructing mixing object
      mix_obj = mix(lattice_obj, charge_obj)

      ! Creating the energy object
      energy_obj = energy(lattice_obj)

      ! Creating hamiltonian object
      hamiltonian_obj = hamiltonian(charge_obj)

      ! Creating recursion object
      recursion_obj = recursion(hamiltonian_obj, energy_obj)

      ! Creating density of states object
      dos_obj = dos(recursion_obj, energy_obj)

      ! Creating Green function object
      green_obj = green(dos_obj)

      ! Creating bands object
      bands_obj = bands(green_obj)

      ! Creating the self object
      self_obj = self(bands_obj, mix_obj)
      call g_timer%start('self-consistency')
      call self_obj%run()
      call g_timer%stop('self-consistency')

      if (this%verbose) then
         !call control_obj%print_state_formatted()
         !call lattice_obj%print_state_formatted()
         !call self_obj%print_state()
         !call charge_obj%print_state()
         call save_state(lattice_obj%symbolic_atoms)
      end if
   end subroutine pre_processing_buildsurf

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Pre-process for bravais calculation
   !---------------------------------------------------------------------------
   subroutine pre_processing_bravais(this)
      class(calculation), intent(in) :: this

      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(energy), target :: energy_obj
      type(self), target :: self_obj
      type(charge), target :: charge_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(recursion), target :: recursion_obj
      type(green), target :: green_obj
      type(dos), target :: dos_obj
      type(bands), target :: bands_obj
      type(mix), target :: mix_obj


      ! Constructing control object
      control_obj = control(this%fname)

      ! Constructing lattice object
      lattice_obj = lattice(control_obj)

      ! Running the pre-calculation
      call g_timer%start('pre-processing')
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      call lattice_obj%structb()

      ! Creating the symbolic_atom object
      call lattice_obj%atomlist()

      ! Initializing MPI lookup tables and info.
      call get_mpi_variables(rank, lattice_obj%nrec)

      ! Constructing the charge object
      charge_obj = charge(lattice_obj)
      call charge_obj%bulkmat()
      call g_timer%stop('pre-processing')

      ! Constructing mixing object
      mix_obj = mix(lattice_obj, charge_obj)

      ! Creating the energy object
      energy_obj = energy(lattice_obj)

      ! Creating hamiltonian object
      hamiltonian_obj = hamiltonian(charge_obj)

      ! Creating recursion object
      recursion_obj = recursion(hamiltonian_obj, energy_obj)

      ! Creating density of states object
      dos_obj = dos(recursion_obj, energy_obj)

      ! Creating Green function object
      green_obj = green(dos_obj)

      ! Creating bands object
      bands_obj = bands(green_obj)

      ! Creating the self object
      self_obj = self(bands_obj, mix_obj)
      call g_timer%start('self-consistency')
      call self_obj%run()
      call g_timer%stop('self-consistency')

      call self_obj%report()

      if (this%verbose) then
         !call control_obj%print_state_formatted()
         !call lattice_obj%print_state_formatted()
         !call self_obj%print_state()
         !call charge_obj%print_state()
         call save_state(lattice_obj%symbolic_atoms)
      end if

      call hamiltonian_obj%rs2pao()
      !call bands_obj%calculate_moments_gauss_legendre()
   end subroutine pre_processing_bravais

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Post-process for to calculate things with Paoflow Hamiltonian
   !---------------------------------------------------------------------------
   subroutine post_processing_paoflow2rs(this)
      class(calculation), intent(in) :: this

      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(energy), target :: energy_obj
      type(self), target :: self_obj
      type(charge), target :: charge_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(recursion), target :: recursion_obj
      type(green), target :: green_obj
      type(dos), target :: dos_obj
      type(bands), target :: bands_obj
      type(mix), target :: mix_obj


      ! Constructing control object
      control_obj = control(this%fname)

      ! Constructing lattice object
      lattice_obj = lattice(control_obj)

      ! Running the pre-calculation
      call g_timer%start('pre-processing')
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      call lattice_obj%structb(.false.)

      ! Creating the symbolic_atom object
      call lattice_obj%atomlist()

      ! Initializing MPI lookup tables and info.
      call get_mpi_variables(rank, lattice_obj%nrec)

      ! Constructing the charge object
      charge_obj = charge(lattice_obj)
      call charge_obj%bulkmat()
      call g_timer%stop('pre-processing')

      ! Constructing mixing object
      mix_obj = mix(lattice_obj, charge_obj)

      ! Creating the energy object
      energy_obj = energy(lattice_obj)
      call energy_obj%e_mesh()

      ! Creating hamiltonian object
      hamiltonian_obj = hamiltonian(charge_obj)
      call hamiltonian_obj%build_from_paoflow_opt()

      ! Creating recursion object
      recursion_obj = recursion(hamiltonian_obj, energy_obj)

      select case (control_obj%recur)
      case ('lanczos')
         call recursion_obj%recur()
      case ('chebyshev')
         call recursion_obj%chebyshev_recur()
      case ('block')
         call recursion_obj%recur_b()
      end select
      ! Creating density of states object
      dos_obj = dos(recursion_obj, energy_obj)

      ! Creating Green function object
      green_obj = green(dos_obj)

      ! Creating bands object
      bands_obj = bands(green_obj)

      ! Creating the self object
      self_obj = self(bands_obj, mix_obj)

      select case (control_obj%recur)
      case ('lanczos')
         call green_obj%sgreen()
      case ('chebyshev')
         call green_obj%chebyshev_green()
      case ('block')
         call recursion_obj%zsqr()
         call green_obj%block_green()
      end select
      call bands_obj%calculate_fermi()
      !call energy_obj%e_mesh()
      !call bands_obj%calculate_moments()
      !call bands_obj%calculate_moments_chebgauss()
      !energy_obj%fermi = 0.0001d0
      !call energy_obj%e_mesh()
      call bands_obj%calculate_moments_gauss_legendre()
   end subroutine post_processing_paoflow2rs

   subroutine post_processing_exchange_p2rs(this)
      class(calculation), intent(in) :: this

      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(energy), target :: energy_obj
      type(self), target :: self_obj
      type(charge), target :: charge_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(recursion), target :: recursion_obj
      type(green), target :: green_obj
      type(dos), target :: dos_obj
      type(bands), target :: bands_obj
      type(mix), target :: mix_obj
      type(exchange), target :: exchange_obj



      ! Constructing control object
      control_obj = control(this%fname)

      ! Constructing lattice object
      lattice_obj = lattice(control_obj)

      ! Running the pre-calculation
      call g_timer%start('pre-processing')
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      call lattice_obj%structb(.false.)

      ! Creating the symbolic_atom object
      call lattice_obj%atomlist()

      ! Initializing MPI lookup tables and info.
      call get_mpi_variables(rank, lattice_obj%njij)

      ! Constructing the charge object
      charge_obj = charge(lattice_obj)
      call charge_obj%bulkmat()
      call g_timer%stop('pre-processing')

      ! Constructing mixing object
      mix_obj = mix(lattice_obj, charge_obj)

      ! Creating the energy object
      energy_obj = energy(lattice_obj)

      ! Creating hamiltonian object
      hamiltonian_obj = hamiltonian(charge_obj)
      select case (control_obj%calctype)
      case ('B')
         call hamiltonian_obj%build_from_paoflow_opt()
      case ('S')
         call g_logger%fatal('Surface calculation not implemented!', __FILE__, __LINE__)
      case ('I')
         call g_logger%fatal('Imputiry calculation not implemented!', __FILE__, __LINE__)
      end select

      ! Creating recursion object
      recursion_obj = recursion(hamiltonian_obj, energy_obj)

      ! Creating density of states object
      dos_obj = dos(recursion_obj, energy_obj)

      ! Creating Green function object
      green_obj = green(dos_obj)

      ! Creating bands object
      bands_obj = bands(green_obj)
      ! Calculate the Fermi energy
      call energy_obj%e_mesh()

      ! Creating the self object
      self_obj = self(bands_obj, mix_obj)

      ! Creating the exchange object
      exchange_obj = exchange(bands_obj)

      ! Calculating the recursion coefficients
      select case (control_obj%recur)
      case ('block')
         call recursion_obj%recur_b_ij()
      case ('chebyshev')
         call recursion_obj%chebyshev_recur_ij()
      end select

      ! Calculating the intersite GFs
      call green_obj%calculate_intersite_gf_eta()
      ! Calculate the heisenberg exhange Jij
      call exchange_obj%calculate_exchange_gauss_legendre()
   end subroutine post_processing_exchange_p2rs

   subroutine post_processing_exchange(this)
      class(calculation), intent(in) :: this

      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(energy), target :: energy_obj
      type(self), target :: self_obj
      type(charge), target :: charge_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(recursion), target :: recursion_obj
      type(green), target :: green_obj
      type(dos), target :: dos_obj
      type(bands), target :: bands_obj
      type(mix), target :: mix_obj
      type(exchange), target :: exchange_obj

      integer :: i

      ! Constructing control object
      control_obj = control(this%fname)

      ! Constructing lattice object
      lattice_obj = lattice(control_obj)

      ! Running the pre-calculation
      call g_timer%start('pre-processing')
      select case (control_obj%calctype)
      case ('B')
         call lattice_obj%build_data()
         call lattice_obj%bravais()
         call lattice_obj%structb()
      case ('S')
         call lattice_obj%build_data()
         call lattice_obj%bravais()
         call lattice_obj%build_surf_full()
         call lattice_obj%structb()
      case ('I')
         call g_logger%fatal('Imputiry calculation not implemented!', __FILE__, __LINE__)
      end select
      ! Creating the symbolic_atom object
      call lattice_obj%atomlist()

      ! Initializing MPI lookup tables and info.
      call get_mpi_variables(rank, lattice_obj%njij)

      ! Constructing the charge object
      charge_obj = charge(lattice_obj)

      select case (control_obj%calctype)
      case ('B')
         call charge_obj%bulkmat()
      case ('S')
         call charge_obj%build_alelay
         call charge_obj%surfmat
      case ('I')
         call g_logger%fatal('Imputiry calculation not implemented!', __FILE__, __LINE__)
      end select
      call g_timer%stop('pre-processing')

      ! Constructing mixing object
      mix_obj = mix(lattice_obj, charge_obj)

      ! Creating the energy object
      energy_obj = energy(lattice_obj)
      call energy_obj%e_mesh()

      ! Creating hamiltonian object
      hamiltonian_obj = hamiltonian(charge_obj)
      select case (control_obj%calctype)
      case ('B')
         do i = 1, lattice_obj%nrec
            call lattice_obj%symbolic_atoms(i)%build_pot() ! Build the potential matrix
         end do
         if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham ! Calculate the spin-orbit coupling Hamiltonian
         call hamiltonian_obj%build_bulkham() ! Build the bulk Hamiltonian
      case ('S')
         do i = 1, lattice_obj%ntype
            call lattice_obj%symbolic_atoms(i)%build_pot() ! Build the potential matrix
         end do
         if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham ! Calculate the spin-orbit coupling Hamiltonian
         call hamiltonian_obj%build_bulkham() ! Build the bulk Hamiltonian for the surface
      case ('I')
         call g_logger%fatal('Imputiry calculation not implemented!', __FILE__, __LINE__)
      end select

      ! Creating recursion object
      recursion_obj = recursion(hamiltonian_obj, energy_obj)

      ! Creating density of states object
      dos_obj = dos(recursion_obj, energy_obj)

      ! Creating Green function object
      green_obj = green(dos_obj)

      ! Creating bands object
      bands_obj = bands(green_obj)

      ! Creating the exchange object
      exchange_obj = exchange(bands_obj)

      ! Calculating the orthogonal parameters
      do i = 1, lattice_obj%ntype
         call lattice_obj%symbolic_atoms(i)%predls(lattice_obj%wav * ang2au)
      end do

      ! Calculating the recursion coefficients
      select case (control_obj%recur)
      case ('block')
         call recursion_obj%recur_b_ij()
      case ('chebyshev')
         call recursion_obj%chebyshev_recur_ij()
      end select

      ! Calculating the intersite GFs
      call green_obj%calculate_intersite_gf()
      ! Calculate the heisenberg exhange Jij
      if ((lattice_obj%njij .ne. 0) .and. (lattice_obj%njijk .eq. 0)) then
         call exchange_obj%calculate_exchange()
         call exchange_obj%calculate_gilbert_damping()
         ! call exchange_obj%calculate_jij_auxgreen()
      else if ((lattice_obj%njijk .ne. 0)) then
         call exchange_obj%calculate_jijk()
      end if
   end subroutine
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default ('none') value
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this)
      class(calculation) :: this
      this%pre_processing = 'none'
      this%processing = 'none'
      this%post_processing = 'none'
   end subroutine restore_to_default

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Check availability for post-processing
   !
   !> @param post_processing Type of processing. Allowed values: 'none'
   !---------------------------------------------------------------------------
   subroutine check_post_processing(post_processing)
      character(len=*), intent(in) :: post_processing
      if (post_processing /= 'none' &
          .and. post_processing /= 'paoflow2rs' &
          .and. post_processing /= 'exchange' &
          .and. post_processing /= 'exchange_p2rs') then
         call g_logger%fatal('[calculation.check_post_processing]: '// &
                             'calculation%post_processing must be one of: ''none'', ''paoflow2rs'', ''exchange'', ''exchange_p2rs''', __FILE__, __LINE__)
      end if
   end subroutine check_post_processing

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Check availability for pre-processing
   !
   !> @param[in] pre_processing Type of pre-processing. Allowed values:
   !> 'bravais', 'buildsurf', 'newclubulk', 'newclusurf', 'none'
   !---------------------------------------------------------------------------
   subroutine check_pre_processing(pre_processing)
      character(len=*), intent(in) :: pre_processing
      if (pre_processing /= 'none' &
          .and. pre_processing /= 'bravais' &
          .and. pre_processing /= 'buildsurf' &
          .and. pre_processing /= 'newclubulk' &
          .and. pre_processing /= 'newclusurf') then
         call g_logger%fatal("[calculation.check_pre_processing]:"// &
                             "calculation%pre_processing must be one of: ''none'', ''bravais'', "// &
                             "''buildsurf'', ''newclusurf'', ''newcluimp''", __FILE__, __LINE__)
      end if
   end subroutine check_pre_processing

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Check availability for processing
   !
   !> @param[in] processing Type of processing. Allowed values: 'none'
   !---------------------------------------------------------------------------
   subroutine check_processing(processing)
      character(len=*), intent(in) :: processing
      if (processing /= 'none' &
          .and. processing /= 'sd') then
         call g_logger%fatal("[calculation.check_processing]: "// &
                             "calculation%processing must be one of: 'none', 'sd' ", __FILE__, __LINE__)
      end if
   end subroutine check_processing
end module calculation_mod
