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
   use conductivity_mod
   use reciprocal_mod
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
      !> ´none´ (default)
      !> ´bravais´ : Builds the bulk clust
      !> ´buildsurf´ : Builds the surface clust
      !> ´newclubulk´ : Builds the imputiry clust from the bluk clust
      !> ´newclusurf´ : Builds the impurity clust from the surface clust
      character(len=sl) :: pre_processing

      !> Processing. Options are
      !> ´none´ (default)
      character(len=sl) :: processing

      !> Post-processing. Options are
      !> ´none´ (default)
      character(len=sl) :: post_processing

      !> Controller for preprocessing verbosity.
      !>
      !> Controller for preprocessing verbosity. If true, call the subroutines:
      !> @ref control.print_state, @ref lattice.print_state,
      !> @ref self.print_state and @ref charge.print_state after pre-processing calls.
      logical :: verbose

      !> name list input file
      character(len=sl) :: fname

      !> K-space and DOS calculation parameters (using &reciprocal namelist for consistency)
      ! K-point mesh settings
      integer :: nk1, nk2, nk3
      real(rp) :: k_offset_x, k_offset_y, k_offset_z
      logical :: use_symmetry_reduction
      logical :: use_time_reversal
      logical :: use_shift
      ! DOS settings
      character(len=sl) :: dos_method
      real(rp) :: gaussian_sigma
      integer :: n_energy_points
      real(rp) :: dos_energy_min
      real(rp) :: dos_energy_max
      real(rp) :: temperature
      real(rp) :: total_electrons
      logical :: auto_find_fermi
      logical :: suppress_internal_logs
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
      procedure, private :: post_processing_conductivity_p2rs
      procedure, private :: post_processing_conductivity
      procedure, private :: post_processing_orbital_modern
      procedure, private :: post_processing_band_structure
      procedure, private :: post_processing_density_of_states
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
   !> @param[in] fname Namelist file with ´&calculation´ table
   !---------------------------------------------------------------------------
   subroutine build_from_file(this, fname)
      class(calculation), intent(inout) :: this
      character(len=*), intent(in) :: fname
      ! variables associated with the reading processes
      integer :: iostatus, funit

      include 'include_codes/namelists/calculation.f90'
      include 'include_codes/namelists/reciprocal.f90'

      verbose = this%verbose
      pre_processing = this%pre_processing
      processing = this%processing
      post_processing = this%post_processing
      ! K-space parameters
      nk1 = this%nk1
      nk2 = this%nk2
      nk3 = this%nk3
      k_offset_x = this%k_offset_x
      k_offset_y = this%k_offset_y
      k_offset_z = this%k_offset_z
      use_symmetry_reduction = this%use_symmetry_reduction
      use_time_reversal = this%use_time_reversal
      use_shift = this%use_shift
      ! DOS parameters
      dos_method = this%dos_method
      gaussian_sigma = this%gaussian_sigma
      n_energy_points = this%n_energy_points
      dos_energy_min = this%dos_energy_min
      dos_energy_max = this%dos_energy_max
      temperature = this%temperature
      total_electrons = this%total_electrons
      auto_find_fermi = this%auto_find_fermi
      suppress_internal_logs = this%suppress_internal_logs

      open (newunit=funit, file=fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//fmt('A', fname)//' not found', __FILE__, __LINE__)
      end if

      read (funit, nml=calculation, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%error('Error while reading namelist', __FILE__, __LINE__)
         call g_logger%error('iostatus = '//fmt('I0', iostatus), __FILE__, __LINE__)
      end if

      ! Read reciprocal space namelist (used by k-space SCF and post-processing)
      read (funit, nml=reciprocal, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%info('Reciprocal namelist not found or error reading it, using defaults', __FILE__, __LINE__)
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
      ! K-space parameters
      this%nk1 = nk1
      this%nk2 = nk2
      this%nk3 = nk3
      this%k_offset_x = k_offset_x
      this%k_offset_y = k_offset_y
      this%k_offset_z = k_offset_z
      this%use_symmetry_reduction = use_symmetry_reduction
      this%use_time_reversal = use_time_reversal
      this%use_shift = use_shift
      ! DOS parameters
      this%dos_method = dos_method
      this%gaussian_sigma = gaussian_sigma
      this%n_energy_points = n_energy_points
      this%dos_energy_min = dos_energy_min
      this%dos_energy_max = dos_energy_max
      this%temperature = temperature
      this%total_electrons = total_electrons
      this%auto_find_fermi = auto_find_fermi
      this%suppress_internal_logs = suppress_internal_logs

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
      case ('conductivity_p2rs')
         call this%post_processing_conductivity_p2rs()
      case ('conductivity')
         call this%post_processing_conductivity()
      case ('orbital_modern')
         call this%post_processing_orbital_modern()
      case ('band_structure')
         call this%post_processing_band_structure()
      case ('density_of_states')
         call this%post_processing_density_of_states()
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
      call lattice_obj%build_surf_full()
      call lattice_obj%newclu()
      call lattice_obj%structb(.true.)

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

      call save_state(lattice_obj%symbolic_atoms)

      call self_obj%report()
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
      call lattice_obj%structb(.true.)

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

      call save_state(lattice_obj%symbolic_atoms)

      call self_obj%report()
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
      call lattice_obj%structb(.true.)

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

      call save_state(lattice_obj%symbolic_atoms)

      call self_obj%report()
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
      call lattice_obj%build_surf_full()
      call lattice_obj%structb(.true.)

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

      call save_state(lattice_obj%symbolic_atoms)

      call self_obj%report()
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
      integer :: i

      ! Constructing control object
      control_obj = control(this%fname)

      ! Constructing lattice object
      lattice_obj = lattice(control_obj)

      ! Running the pre-calculation
      call g_timer%start('pre-processing')
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      call lattice_obj%structb(.true.)

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

      call save_state(lattice_obj%symbolic_atoms)

      call hamiltonian_obj%rs2txt()
      !call hamiltonian_obj%rs2pao()
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
      integer :: i

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

      call self_obj%report()

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
      real(rp), dimension(6) :: QSL
      integer :: i, ia

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
      real(rp), dimension(6) :: QSL
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
         call lattice_obj%structb(.true.)
      case ('S')
         call lattice_obj%build_data()
         call lattice_obj%bravais()
         call lattice_obj%build_surf_full()
         call lattice_obj%structb(.true.)
      case ('I')
         call lattice_obj%build_data()
         call lattice_obj%bravais()
         call lattice_obj%build_surf_full()
         call lattice_obj%newclu()
         call lattice_obj%structb(.true.)
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
         call charge_obj%impmad()
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
         do i = 1, lattice_obj%ntype
            call lattice_obj%symbolic_atoms(i)%build_pot() ! Build the potential matrix
         end do
         if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham ! Calculate the spin-orbit coupling Hamiltonian
         call hamiltonian_obj%build_bulkham() ! Build the bulk Hamiltonian
         call hamiltonian_obj%build_locham() ! Build the local Hamiltonian
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
         call lattice_obj%symbolic_atoms(i)%predls(lattice_obj%wav*ang2au)
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
      call green_obj%calculate_intersite_gf_twoindex()
      ! Calculate the heisenberg exhange Jij
      if ((lattice_obj%njij .ne. 0) .and. (lattice_obj%njijk .eq. 0)) then
         call exchange_obj%calculate_exchange()
         call exchange_obj%calculate_exchange_twoindex()
         !call exchange_obj%calculate_gilbert_damping()
         !call exchange_obj%calculate_moment_of_inertia()
         ! call exchange_obj%calculate_jij_auxgreen()
         !else if ((lattice_obj%njijk.ne.0)) then
         !  call exchange_obj%calculate_jijk()
      end if
   end subroutine


   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Post-process for band structure calculation using reciprocal space module
   !---------------------------------------------------------------------------
   subroutine post_processing_band_structure(this)
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
      type(reciprocal), target :: reciprocal_obj
      integer :: i

      call g_logger%info('post_processing_band_structure: Starting band structure calculation', __FILE__, __LINE__)

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
         call lattice_obj%structb(.true.)
      case ('S')
         call lattice_obj%build_data()
         call lattice_obj%bravais()
         call lattice_obj%build_surf_full()
         call lattice_obj%structb(.true.)
      case ('I')
         call lattice_obj%build_data()
         call lattice_obj%bravais()
         call lattice_obj%build_surf_full()
         call lattice_obj%newclu()
         call lattice_obj%structb(.true.)
      end select

      ! Creating the symbolic_atom object
      call lattice_obj%atomlist()

      ! Initializing MPI lookup tables and info.
      call get_mpi_variables(rank, lattice_obj%nrec)

      ! Constructing the charge object
      charge_obj = charge(lattice_obj)

      select case (control_obj%calctype)
      case ('B')
         call charge_obj%bulkmat()
      case ('S')
         call charge_obj%build_alelay
         call charge_obj%surfmat
      case ('I')
         call charge_obj%impmad()
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
         do i = 1, lattice_obj%ntype
            call lattice_obj%symbolic_atoms(i)%build_pot() ! Build the potential matrix
         end do
         if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham ! Calculate the spin-orbit coupling Hamiltonian
         call hamiltonian_obj%build_bulkham() ! Build the bulk Hamiltonian
         call hamiltonian_obj%build_locham() ! Build the local Hamiltonian
      end select

      ! Check if Hamiltonian was built successfully
      if (.not. allocated(hamiltonian_obj%ee)) then
         call g_logger%error('post_processing_band_structure: Bulk Hamiltonian not allocated after build', __FILE__, __LINE__)
         return
      end if

      call g_logger%info('post_processing_band_structure: Hamiltonian dimensions: ' // &
         fmt('I0', size(hamiltonian_obj%ee, 1)) // 'x' // &
         fmt('I0', size(hamiltonian_obj%ee, 2)) // 'x' // &
         fmt('I0', size(hamiltonian_obj%ee, 3)) // 'x' // &
         fmt('I0', size(hamiltonian_obj%ee, 4)), __FILE__, __LINE__)

      ! Creating reciprocal space object
      reciprocal_obj = reciprocal(hamiltonian_obj)

      call g_logger%info('post_processing_band_structure: Setting up reciprocal space', __FILE__, __LINE__)

      ! Debug lattice information (sbarvec is populated on-demand, not a requirement here)
      if (allocated(lattice_obj%sbarvec)) then
         call g_logger%info('post_processing_band_structure: Lattice sbarvec already populated with ' // &
            fmt('I0', size(lattice_obj%sbarvec, 2)) // ' neighbors', __FILE__, __LINE__)
      else
         call g_logger%info('post_processing_band_structure: Lattice sbarvec not yet allocated (will be populated if needed)', __FILE__, __LINE__)
      end if

      ! Set k-point mesh parameters (use &reciprocal namelist values)
      reciprocal_obj%nk_mesh = [this%nk1, this%nk2, this%nk3]
      reciprocal_obj%k_offset = [this%k_offset_x, this%k_offset_y, this%k_offset_z]
      reciprocal_obj%use_time_reversal = this%use_time_reversal
      reciprocal_obj%use_symmetry_reduction = this%use_symmetry_reduction

      call g_logger%info('post_processing_band_structure: K-mesh = ' // &
         fmt('I0', this%nk1) // 'x' // fmt('I0', this%nk2) // 'x' // fmt('I0', this%nk3), __FILE__, __LINE__)

      ! Generate k-point mesh for general calculations
      call reciprocal_obj%generate_mp_mesh()

      ! Build k-space Hamiltonian
      call reciprocal_obj%build_kspace_hamiltonian()

      call g_logger%info('post_processing_band_structure: Calculating band structure along high-symmetry path', __FILE__, __LINE__)

      ! Automatically determine crystal structure
      ! In future, this could be read from input file
      block
         character(len=10) :: detected_crystal_type
         detected_crystal_type = reciprocal_obj%symmetry_analysis%determine_crystal_structure()
         call g_logger%info('post_processing_band_structure: Detected crystal type: ' // trim(detected_crystal_type), __FILE__, __LINE__)
         call reciprocal_obj%calculate_band_structure(hamiltonian_obj, detected_crystal_type, 50, 'band_structure.dat')
      end block

      call g_logger%info('post_processing_band_structure: Band structure calculation completed', __FILE__, __LINE__)
      call g_logger%info('post_processing_band_structure: Results written to band_structure.dat and band_structure_kpath.dat', __FILE__, __LINE__)

      ! Clean up
      call reciprocal_obj%restore_to_default()

   end subroutine post_processing_band_structure

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Post-process for density of states calculation using reciprocal space module
   !---------------------------------------------------------------------------
   subroutine post_processing_density_of_states(this)
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
      type(reciprocal), target :: reciprocal_obj
      integer :: i

      call g_logger%info('post_processing_density_of_states: Starting DOS calculation', __FILE__, __LINE__)

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
         call lattice_obj%structb(.true.)
      case ('S')
         call lattice_obj%build_data()
         call lattice_obj%bravais()
         call lattice_obj%build_surf_full()
         call lattice_obj%structb(.true.)
      case ('I')
         call lattice_obj%build_data()
         call lattice_obj%bravais()
         call lattice_obj%build_surf_full()
         call lattice_obj%newclu()
         call lattice_obj%structb(.true.)
      end select

      ! Creating the symbolic_atom object
      call lattice_obj%atomlist()

      ! Initializing MPI lookup tables and info.
      call get_mpi_variables(rank, lattice_obj%nrec)

      ! Constructing the charge object
      charge_obj = charge(lattice_obj)

      select case (control_obj%calctype)
      case ('B')
         call charge_obj%bulkmat()
      case ('S')
         call charge_obj%build_alelay
         call charge_obj%surfmat
      case ('I')
         call charge_obj%impmad()
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
         do i = 1, lattice_obj%ntype
            call lattice_obj%symbolic_atoms(i)%build_pot() ! Build the potential matrix
         end do
         if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham ! Calculate the spin-orbit coupling Hamiltonian
         call hamiltonian_obj%build_bulkham() ! Build the bulk Hamiltonian
         call hamiltonian_obj%build_locham() ! Build the local Hamiltonian
      end select

      ! Check if Hamiltonian was built successfully
      if (.not. allocated(hamiltonian_obj%ee)) then
         call g_logger%error('post_processing_density_of_states: Bulk Hamiltonian not allocated after build', __FILE__, __LINE__)
         return
      end if

      call g_logger%info('post_processing_density_of_states: Hamiltonian dimensions: ' // &
         fmt('I0', size(hamiltonian_obj%ee, 1)) // 'x' // &
         fmt('I0', size(hamiltonian_obj%ee, 2)) // 'x' // &
         fmt('I0', size(hamiltonian_obj%ee, 3)) // 'x' // &
         fmt('I0', size(hamiltonian_obj%ee, 4)), __FILE__, __LINE__)

      ! Creating reciprocal space object
      reciprocal_obj = reciprocal(hamiltonian_obj)

      call g_logger%info('post_processing_density_of_states: Setting up reciprocal space', __FILE__, __LINE__)

      ! Debug lattice information
      ! Debug lattice information (sbarvec is populated on-demand, not a requirement here)
      if (allocated(lattice_obj%sbarvec)) then
         call g_logger%info('post_processing_density_of_states: Lattice sbarvec already populated with ' // &
            fmt('I0', size(lattice_obj%sbarvec, 2)) // ' neighbors', __FILE__, __LINE__)
      else
         call g_logger%info('post_processing_density_of_states: Lattice sbarvec not yet allocated (will be populated if needed)', __FILE__, __LINE__)
      end if

      ! Set k-point mesh parameters from &reciprocal namelist
      reciprocal_obj%nk_mesh = [this%nk1, this%nk2, this%nk3]
      reciprocal_obj%k_offset = [this%k_offset_x, this%k_offset_y, this%k_offset_z]
      reciprocal_obj%use_time_reversal = this%use_time_reversal
      reciprocal_obj%use_symmetry_reduction = this%use_symmetry_reduction
      
      ! Set DOS parameters from &reciprocal namelist
      reciprocal_obj%dos_method = this%dos_method
      reciprocal_obj%gaussian_sigma = this%gaussian_sigma
      reciprocal_obj%n_energy_points = this%n_energy_points
      reciprocal_obj%dos_energy_range = [this%dos_energy_min, this%dos_energy_max]
      reciprocal_obj%temperature = this%temperature
      
      ! Enable automatic Fermi level finding by default for post-processing
      ! If user specified total_electrons in namelist, use it; otherwise auto-calculate from valence
      if (this%total_electrons > 0.0_rp) then
         reciprocal_obj%total_electrons = this%total_electrons
         call g_logger%info('post_processing_density_of_states: Using user-specified total_electrons = ' // &
                           fmt('F10.5', this%total_electrons), __FILE__, __LINE__)
      else
         ! Auto-calculate from valence electrons
         reciprocal_obj%total_electrons = real(sum(lattice_obj%symbolic_atoms(1:lattice_obj%nbulk_bulk)%element%valence), rp)
         call g_logger%info('post_processing_density_of_states: Auto-calculated total_electrons = ' // &
                           fmt('F10.5', reciprocal_obj%total_electrons) // ' from valence', __FILE__, __LINE__)
      end if
      
      ! Enable Fermi level auto-finding unless explicitly disabled
      reciprocal_obj%auto_find_fermi = .true.
      if (.not. this%auto_find_fermi) then
         reciprocal_obj%auto_find_fermi = .false.
         call g_logger%info('post_processing_density_of_states: Fermi level auto-finding disabled by user', __FILE__, __LINE__)
      end if
      
      reciprocal_obj%suppress_internal_logs = this%suppress_internal_logs

      call g_logger%info('post_processing_density_of_states: K-mesh = ' // &
         fmt('I0', this%nk1) // 'x' // fmt('I0', this%nk2) // 'x' // fmt('I0', this%nk3), __FILE__, __LINE__)
      call g_logger%info('post_processing_density_of_states: DOS method = ' // trim(this%dos_method), __FILE__, __LINE__)

      ! Generate k-point mesh for DOS calculations
      call reciprocal_obj%generate_mp_mesh()

      ! Build k-space Hamiltonian
      call reciprocal_obj%build_kspace_hamiltonian()

      call g_logger%info('post_processing_density_of_states: Calculating density of states', __FILE__, __LINE__)

      ! Calculate DOS - parameters already set on reciprocal_obj above
      call reciprocal_obj%calculate_density_of_states(hamiltonian_obj, &
                                                     fermi_level=energy_obj%fermi, &
                                                     output_file='density_of_states.dat')

      call g_logger%info('post_processing_density_of_states: DOS calculation completed', __FILE__, __LINE__)
      call g_logger%info('post_processing_density_of_states: Results written to density_of_states.dat', __FILE__, __LINE__)

      ! Clean up
      call reciprocal_obj%restore_to_default()

   end subroutine post_processing_density_of_states


   subroutine post_processing_conductivity(this)
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
      type(conductivity), target :: conductivity_obj
      real(rp), dimension(6) :: QSL
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
         call lattice_obj%structb(.true.)
      case ('S')
         call lattice_obj%build_data()
         call lattice_obj%bravais()
         call lattice_obj%build_surf_full()
         call lattice_obj%structb(.true.)
      case ('I')
         call lattice_obj%build_data()
         call lattice_obj%bravais()
         call lattice_obj%build_surf_full()
         call lattice_obj%newclu()
         call lattice_obj%structb(.true.)
      end select
      ! Creating the symbolic_atom object
      call lattice_obj%atomlist()

      ! Initializing MPI lookup tables and info.
      call get_mpi_variables(rank, lattice_obj%ntype)

      ! Constructing the charge object
      charge_obj = charge(lattice_obj)

      select case (control_obj%calctype)
      case ('B')
         call charge_obj%bulkmat()
      case ('S')
         call charge_obj%build_alelay
         call charge_obj%surfmat
      case ('I')
         call charge_obj%impmad()
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
         do i = 1, lattice_obj%ntype
            call lattice_obj%symbolic_atoms(i)%build_pot() ! Build the potential matrix
         end do
         if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham ! Calculate the spin-orbit coupling Hamiltonian
         call hamiltonian_obj%build_bulkham() ! Build the bulk Hamiltonian
         call hamiltonian_obj%build_locham() ! Build the local Hamiltonian
      end select

      ! Creating recursion object
      recursion_obj = recursion(hamiltonian_obj, energy_obj)
      
      call recursion_obj%compute_moments_stochastic()

      ! Creating density of states object
      dos_obj = dos(recursion_obj, energy_obj)

      ! Creating Green function object
      green_obj = green(dos_obj)

      ! Creating bands object
      bands_obj = bands(green_obj)
   
      ! Creating the self object
      self_obj = self(bands_obj, mix_obj)

      ! Creating the conductivity object
      conductivity_obj = conductivity(self_obj)

      call conductivity_obj%calculate_gamma_nm()
      call conductivity_obj%calculate_conductivity_tensor()
   end subroutine post_processing_conductivity


   subroutine post_processing_conductivity_p2rs(this)
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
      type(conductivity), target :: conductivity_obj
      real(rp), dimension(6) :: QSL
      integer :: i

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
      call get_mpi_variables(rank, lattice_obj%ntype)

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
      
      call recursion_obj%compute_moments_stochastic()

      ! Creating density of states object
      dos_obj = dos(recursion_obj, energy_obj)

      ! Creating Green function object
      green_obj = green(dos_obj)

      ! Creating bands object
      bands_obj = bands(green_obj)
   
      ! Creating the self object
      self_obj = self(bands_obj, mix_obj)

      ! Creating the conductivity object
      conductivity_obj = conductivity(self_obj)

      call conductivity_obj%calculate_gamma_nm()
      call conductivity_obj%calculate_conductivity_tensor()
   end subroutine post_processing_conductivity_p2rs


   subroutine post_processing_orbital_modern(this)
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
      type(conductivity), target :: conductivity_obj
      real(rp), dimension(6) :: QSL
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
         call lattice_obj%structb(.true.)
      case ('S')
         call lattice_obj%build_data()
         call lattice_obj%bravais()
         call lattice_obj%build_surf_full()
         call lattice_obj%structb(.true.)
      case ('I')
         call lattice_obj%build_data()
         call lattice_obj%bravais()
         call lattice_obj%build_surf_full()
         call lattice_obj%newclu()
         call lattice_obj%structb(.true.)
      end select
      ! Creating the symbolic_atom object
      call lattice_obj%atomlist()

      ! Initializing MPI lookup tables and info.
      call get_mpi_variables(rank, lattice_obj%ntype)

      ! Constructing the charge object
      charge_obj = charge(lattice_obj)

      select case (control_obj%calctype)
      case ('B')
         call charge_obj%bulkmat()
      case ('S')
         call charge_obj%build_alelay
         call charge_obj%surfmat
      case ('I')
         call charge_obj%impmad()
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
         do i = 1, lattice_obj%ntype
            call lattice_obj%symbolic_atoms(i)%build_pot() ! Build the potential matrix
         end do
         if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham ! Calculate the spin-orbit coupling Hamiltonian
         call hamiltonian_obj%build_bulkham() ! Build the bulk Hamiltonian
         call hamiltonian_obj%build_locham() ! Build the local Hamiltonian
      end select

      ! Creating recursion object
      recursion_obj = recursion(hamiltonian_obj, energy_obj)
      
      call recursion_obj%chebyshev_orbital_mod()

      ! Creating density of states object
      dos_obj = dos(recursion_obj, energy_obj)

      ! Creating Green function object
      green_obj = green(dos_obj)

      ! Creating bands object
      bands_obj = bands(green_obj)
   
      ! Creating the self object
      self_obj = self(bands_obj, mix_obj)
   end subroutine post_processing_orbital_modern



   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default (´none´) value
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this)
      class(calculation) :: this
      this%pre_processing = 'none'
      this%processing = 'none'
      this%post_processing = 'none'
      ! K-space defaults (matching reciprocal module defaults)
      this%nk1 = 8
      this%nk2 = 8
      this%nk3 = 8
      this%k_offset_x = 0.0_rp
      this%k_offset_y = 0.0_rp
      this%k_offset_z = 0.0_rp
      this%use_symmetry_reduction = .false.
      this%use_time_reversal = .true.
      this%use_shift = .false.
      ! DOS defaults (matching reciprocal module defaults)
      this%dos_method = 'gaussian'
      this%gaussian_sigma = 0.1_rp
      this%n_energy_points = 1000
      this%dos_energy_min = -10.0_rp
      this%dos_energy_max = 10.0_rp
      this%temperature = 300.0_rp
      this%total_electrons = 0.0_rp  ! 0 = auto-detect from valence
      this%auto_find_fermi = .true.
      this%suppress_internal_logs = .true.
   end subroutine restore_to_default

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Check availability for post-processing
   !
   !> @param post_processing Type of processing. Allowed values: ´none´
   !---------------------------------------------------------------------------
   subroutine check_post_processing(post_processing)
      character(len=*), intent(in) :: post_processing
      if (post_processing /= 'none' &
          .and. post_processing /= 'paoflow2rs' &
          .and. post_processing /= 'exchange' &
          .and. post_processing /= 'exchange_p2rs' &
          .and. post_processing /= 'conductivity' &
          .and. post_processing /= 'conductivity_p2rs' &
          .and. post_processing /= 'orbital_modern' &
          .and. post_processing /= 'band_structure' &
          .and. post_processing /= 'density_of_states') then 
         call g_logger%fatal('[calculation.check_post_processing]: '// &
                             "calculation%post_processing must be one of: ''none'', ''paoflow2rs'', ''exchange'', ''exchange_p2rs''," // &
                             " 'conductivity', 'conductivity_p2rs', 'orbital_modern', 'band_structure', 'density_of_states'", __FILE__, __LINE__)
      end if
   end subroutine check_post_processing

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Check availability for pre-processing
   !
   !> @param[in] pre_processing Type of pre-processing. Allowed values:
   !> ´bravais´, ´buildsurf´, ´newclubulk´, ´newclusurf´, ´none´
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
   !> @param[in] processing Type of processing. Allowed values: ´none´
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
