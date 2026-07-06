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
   use sparse_mod, only: sparse
   use recursion_mod
   use green_mod
   use density_of_states_mod
   use bands_mod
   use exchange_mod
   use spin_dynamics_mod
   use conductivity_mod
   use reciprocal_mod
   use mix_mod
   use frozen_magnon_mod
   use math_mod
   use precision_mod, only: rp
   use string_mod, only: sl, fmt, real2str
   use timer_mod, only: g_timer
   use logger_mod, only: g_logger
   use basis_mod, only: basis_init
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
      procedure, private :: post_processing_frozen_magnon
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
      case ('frozen_magnon')
         call this%post_processing_frozen_magnon()
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

      ! Initialize basis dimension parameters from lmax
      call basis_init(lattice_obj%symbolic_atoms(1)%potential%lmax)

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
      recursion_obj = recursion(hamiltonian_obj, energy_obj, sparse(hamiltonian_obj))

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
      recursion_obj = recursion(hamiltonian_obj, energy_obj, sparse(hamiltonian_obj))

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

      ! Initialize basis dimension parameters from lmax
      call basis_init(lattice_obj%symbolic_atoms(1)%potential%lmax)

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
      recursion_obj = recursion(hamiltonian_obj, energy_obj, sparse(hamiltonian_obj))

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

      ! Initialize basis dimension parameters from lmax
      call basis_init(lattice_obj%symbolic_atoms(1)%potential%lmax)

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
      recursion_obj = recursion(hamiltonian_obj, energy_obj, sparse(hamiltonian_obj))

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

      ! Initialize basis dimension parameters from lmax
      call basis_init(lattice_obj%symbolic_atoms(1)%potential%lmax)

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
      recursion_obj = recursion(hamiltonian_obj, energy_obj, sparse(hamiltonian_obj))

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
      recursion_obj = recursion(hamiltonian_obj, energy_obj, sparse(hamiltonian_obj))

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

      select case (trim(hamiltonian_obj%export))
      case ('rs2pao')
         call hamiltonian_obj%rs2pao()
      case ('python')
         call hamiltonian_obj%export_rs_tb_all()
      end select
      call bands_obj%calculate_orbital_quadrupoles()
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
      recursion_obj = recursion(hamiltonian_obj, energy_obj, sparse(hamiltonian_obj))

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
         call green_obj%chebyshev_dos_dispatch()  ! Dispatches to GPU/C++/legacy based on control flags
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

   subroutine prepare_post_processing_stack(this, use_paoflow, use_exchange_pairs, energy_mesh_before_hamiltonian, &
                                             stochastic_moments, control_obj, lattice_obj, charge_obj, mix_obj, &
                                             energy_obj, hamiltonian_obj, recursion_obj, dos_obj, green_obj, bands_obj)
      class(calculation), intent(in) :: this
      logical, intent(in) :: use_paoflow
      logical, intent(in) :: use_exchange_pairs
      logical, intent(in) :: energy_mesh_before_hamiltonian
      logical, intent(in) :: stochastic_moments
      type(control), target, intent(inout) :: control_obj
      type(lattice), target, intent(inout) :: lattice_obj
      type(charge), target, intent(inout) :: charge_obj
      type(mix), target, intent(inout) :: mix_obj
      type(energy), target, intent(inout) :: energy_obj
      type(hamiltonian), target, intent(inout) :: hamiltonian_obj
      type(recursion), target, intent(inout) :: recursion_obj
      type(dos), target, intent(inout) :: dos_obj
      type(green), target, intent(inout) :: green_obj
      type(bands), target, intent(inout) :: bands_obj
      integer :: i

      control_obj = control(this%fname)
      lattice_obj = lattice(control_obj)

      call g_timer%start('pre-processing')
      if (use_paoflow) then
         call lattice_obj%build_data()
         call lattice_obj%bravais()
         call lattice_obj%structb(.false.)
      else
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
      end if

      call lattice_obj%atomlist()
      if (use_exchange_pairs) then
         call get_mpi_variables(rank, lattice_obj%njij)
      else
         call get_mpi_variables(rank, lattice_obj%ntype)
      end if

      charge_obj = charge(lattice_obj)
      if (use_paoflow) then
         call charge_obj%bulkmat()
      else
         select case (control_obj%calctype)
         case ('B')
            call charge_obj%bulkmat()
         case ('S')
            call charge_obj%build_alelay
            call charge_obj%surfmat
         case ('I')
            call charge_obj%impmad()
         end select
      end if
      call g_timer%stop('pre-processing')

      mix_obj = mix(lattice_obj, charge_obj)
      energy_obj = energy(lattice_obj)
      if (energy_mesh_before_hamiltonian) call energy_obj%e_mesh()

      hamiltonian_obj = hamiltonian(charge_obj)
      if (use_paoflow) then
         select case (control_obj%calctype)
         case ('B')
            call hamiltonian_obj%build_from_paoflow_opt()
         case ('S')
            call g_logger%fatal('Surface calculation not implemented!', __FILE__, __LINE__)
         case ('I')
            call g_logger%fatal('Imputiry calculation not implemented!', __FILE__, __LINE__)
         end select
      else
         select case (control_obj%calctype)
         case ('B')
            do i = 1, lattice_obj%nrec
               call lattice_obj%symbolic_atoms(i)%build_pot()
            end do
            if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham
            call hamiltonian_obj%build_bulkham()
         case ('S')
            do i = 1, lattice_obj%ntype
               call lattice_obj%symbolic_atoms(i)%build_pot()
            end do
            if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham
            call hamiltonian_obj%build_bulkham()
         case ('I')
            do i = 1, lattice_obj%ntype
               call lattice_obj%symbolic_atoms(i)%build_pot()
            end do
            if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham
            call hamiltonian_obj%build_bulkham()
            call hamiltonian_obj%build_locham()
         end select
      end if

      recursion_obj = recursion(hamiltonian_obj, energy_obj, sparse(hamiltonian_obj))
      if (stochastic_moments) call recursion_obj%compute_moments_stochastic()
      dos_obj = dos(recursion_obj, energy_obj)
      green_obj = green(dos_obj)
      bands_obj = bands(green_obj)
      if (.not. energy_mesh_before_hamiltonian) call energy_obj%e_mesh()
   end subroutine prepare_post_processing_stack

   subroutine run_intersite_moments(control_obj, recursion_obj)
      type(control), intent(in) :: control_obj
      type(recursion), intent(inout) :: recursion_obj

      select case (control_obj%recur)
      case ('block')
         call recursion_obj%recur_b_ij()
      case ('chebyshev')
         call recursion_obj%chebyshev_recur_ij()
      end select
   end subroutine run_intersite_moments

   subroutine finish_conductivity_moments(green_obj, bands_obj)
      type(green), intent(inout) :: green_obj
      type(bands), intent(inout) :: bands_obj

      call green_obj%chebyshev_dos_dispatch()
      call bands_obj%calculate_fermi()
   end subroutine finish_conductivity_moments

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

      call prepare_post_processing_stack(this, .true., .true., .false., .false., control_obj, lattice_obj, &
                                         charge_obj, mix_obj, energy_obj, hamiltonian_obj, recursion_obj, dos_obj, &
                                         green_obj, bands_obj)
      self_obj = self(bands_obj, mix_obj)
      exchange_obj = exchange(bands_obj)

      call run_intersite_moments(control_obj, recursion_obj)
      call green_obj%calculate_intersite_gf_eta()
      call exchange_obj%calculate_exchange_gauss_legendre()
   end subroutine post_processing_exchange_p2rs

   subroutine post_processing_exchange(this)
      class(calculation), intent(in) :: this

      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(energy), target :: energy_obj
      type(charge), target :: charge_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(recursion), target :: recursion_obj
      type(green), target :: green_obj
      type(dos), target :: dos_obj
      type(bands), target :: bands_obj
      type(mix), target :: mix_obj
      type(exchange), target :: exchange_obj
      integer :: i

      call prepare_post_processing_stack(this, .false., .true., .true., .false., control_obj, lattice_obj, &
                                         charge_obj, mix_obj, energy_obj, hamiltonian_obj, recursion_obj, dos_obj, &
                                         green_obj, bands_obj)
      exchange_obj = exchange(bands_obj)

      do i = 1, lattice_obj%ntype
         call lattice_obj%symbolic_atoms(i)%predls(lattice_obj%wav*ang2au)
      end do

      call run_intersite_moments(control_obj, recursion_obj)
      call green_obj%calculate_intersite_gf()
      call green_obj%calculate_intersite_gf_twoindex()
      if ((lattice_obj%njij .ne. 0) .and. (lattice_obj%njijk .eq. 0)) then
         call exchange_obj%calculate_exchange()
         call exchange_obj%calculate_exchange_twoindex()
      end if
   end subroutine post_processing_exchange


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
      type(conductivity), target :: conductivity_obj

      call prepare_post_processing_stack(this, .false., .false., .true., .true., control_obj, lattice_obj, &
                                         charge_obj, mix_obj, energy_obj, hamiltonian_obj, recursion_obj, dos_obj, &
                                         green_obj, bands_obj)
      self_obj = self(bands_obj, mix_obj)
      call finish_conductivity_moments(green_obj, bands_obj)
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
      type(conductivity), target :: conductivity_obj

      call prepare_post_processing_stack(this, .true., .false., .true., .true., control_obj, lattice_obj, &
                                         charge_obj, mix_obj, energy_obj, hamiltonian_obj, recursion_obj, dos_obj, &
                                         green_obj, bands_obj)
      self_obj = self(bands_obj, mix_obj)
      call finish_conductivity_moments(green_obj, bands_obj)
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
      recursion_obj = recursion(hamiltonian_obj, energy_obj, sparse(hamiltonian_obj))
      
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

   subroutine post_processing_band_structure(this)
      class(calculation), intent(in) :: this
      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(charge), target :: charge_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(reciprocal), target :: reciprocal_obj
      integer :: i

      control_obj = control(this%fname)
      lattice_obj = lattice(control_obj)
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      call lattice_obj%structb(.true.)
      call lattice_obj%atomlist()
      charge_obj = charge(lattice_obj)
      call charge_obj%bulkmat()
      hamiltonian_obj = hamiltonian(charge_obj)
      do i = 1, lattice_obj%nrec
         call lattice_obj%symbolic_atoms(i)%build_pot()
      end do
      if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham()
      call hamiltonian_obj%build_bulkham()
      reciprocal_obj = reciprocal(hamiltonian_obj)
      call reciprocal_obj%calculate_band_structure(hamiltonian_obj, 'auto', reciprocal_obj%nk_per_segment, 'band_structure.dat')
   end subroutine post_processing_band_structure

   subroutine post_processing_density_of_states(this)
      class(calculation), intent(in) :: this
      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(charge), target :: charge_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(reciprocal), target :: reciprocal_obj
      integer :: i

      control_obj = control(this%fname)
      lattice_obj = lattice(control_obj)
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      call lattice_obj%structb(.true.)
      call lattice_obj%atomlist()
      charge_obj = charge(lattice_obj)
      call charge_obj%bulkmat()
      hamiltonian_obj = hamiltonian(charge_obj)
      do i = 1, lattice_obj%nrec
         call lattice_obj%symbolic_atoms(i)%build_pot()
      end do
      if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham()
      call hamiltonian_obj%build_bulkham()
      reciprocal_obj = reciprocal(hamiltonian_obj)
      call reciprocal_obj%calculate_density_of_states(hamiltonian_obj, reciprocal_obj%n_energy_points, &
           reciprocal_obj%dos_energy_range, reciprocal_obj%dos_method, reciprocal_obj%gaussian_sigma, &
           reciprocal_obj%temperature, reciprocal_obj%fermi_level, reciprocal_obj%total_electrons, &
           reciprocal_obj%auto_find_fermi, 'dos_kspace.dat')
   end subroutine post_processing_density_of_states

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Frozen-magnon post-processing: sweeps &hamiltonian's q_ss over a
   !> user-supplied list of points (see &frozen_magnon namelist) and reports
   !> the adiabatic magnon dispersion.
   !> @details For each q in the sweep, the Hamiltonian is rebuilt with that
   !> q_ss and the potential is either fully re-converged (mode='scf') or, for
   !> mode='mft' (the magnetic force theorem default), converged once at the
   !> reference point q_ss_list(:,1) and then re-evaluated with a single
   !> band-energy pass (self%nstep=1) at every other q, reusing the reference
   !> potential unchanged. The dispersion omega(q) = 4*(E(q)-E(q_ref)) /
   !> (M_tot*sin^2(theta_ss)) uses the band energy for mode='mft' and the
   !> total energy for mode='scf'. This is the single-acoustic-branch
   !> generalization of the Halilov frozen-magnon formula to multiple
   !> sublattices: M_tot is the reference point's total moment summed over
   !> sublattices, correct because
   !> every sublattice shares the same global (q_ss, theta_ss) cant (the
   !> uniform-cant ansatz already built into ham0m_nc/hamiltonian_ccor.f90).
   !> Independent per-sublattice magnon branches (e.g. optic modes in an
   !> antiferromagnetic cell) are out of scope; only the acoustic branch is
   !> captured here.
   !---------------------------------------------------------------------------
   subroutine post_processing_frozen_magnon(this)
      class(calculation), intent(in) :: this
      type(control), target :: control_obj
      type(lattice), target :: lattice_obj
      type(energy), target :: energy_obj
      type(charge), target :: charge_obj
      type(hamiltonian), target :: hamiltonian_obj
      type(recursion), target :: recursion_obj
      type(green), target :: green_obj
      type(dos), target :: dos_obj
      type(bands), target :: bands_obj
      type(mix), target :: mix_obj
      type(self), target :: self_obj
      type(frozen_magnon) :: fm_obj
      real(rp), allocatable :: etot_q(:), eband_q(:), mtot_q(:, :), omega_q(:), q_ss_cart(:, :)
      real(rp) :: sin2theta
      real(rp), dimension(3, 3) :: direct_to_cart
      character(len=200) :: fmt_str
      integer :: iq, i, newunit

      fm_obj = frozen_magnon(this%fname)

      call prepare_post_processing_stack(this, .false., .false., .true., .false., control_obj, lattice_obj, &
                                         charge_obj, mix_obj, energy_obj, hamiltonian_obj, recursion_obj, dos_obj, &
                                         green_obj, bands_obj)

      allocate (q_ss_cart(3, fm_obj%n_q))
      if (fm_obj%q_coordinates == 'direct') then
         direct_to_cart = transpose(inverse_3x3(lattice_obj%a))
         do iq = 1, fm_obj%n_q
            q_ss_cart(:, iq) = matmul(direct_to_cart, fm_obj%q_ss_list(:, iq))
         end do
      else
         q_ss_cart(:, :) = fm_obj%q_ss_list(:, :)
      end if

      allocate (etot_q(fm_obj%n_q), eband_q(fm_obj%n_q), mtot_q(lattice_obj%nrec, fm_obj%n_q), omega_q(fm_obj%n_q))

      do iq = 1, fm_obj%n_q
         hamiltonian_obj%q_ss(:) = q_ss_cart(:, iq)

         self_obj = self(bands_obj, mix_obj)
         if (fm_obj%mode == 'mft' .and. iq > 1) self_obj%nstep = 1
         call self_obj%run()

         etot_q(iq) = sum(lattice_obj%symbolic_atoms(:)%potential%etot)
         if (.not. self_obj%use_kspace) call bands_obj%calculate_band_energy()
         eband_q(iq) = bands_obj%eband
         do i = 1, lattice_obj%nrec
            mtot_q(i, iq) = lattice_obj%symbolic_atoms(lattice_obj%nbulk + i)%potential%mtot
         end do
      end do

      sin2theta = sin(hamiltonian_obj%theta_ss)**2
      if (sin2theta < 1.0e-8_rp) then
         call g_logger%fatal('[calculation.post_processing_frozen_magnon]: '// &
                             'theta_ss must be nonzero (a finite cone angle) to define omega(q)', __FILE__, __LINE__)
      end if
      do iq = 1, fm_obj%n_q
         if (fm_obj%mode == 'mft') then
            omega_q(iq) = 4.0_rp*(eband_q(iq) - eband_q(1))/(sum(mtot_q(:, 1))*sin2theta)
         else
            omega_q(iq) = 4.0_rp*(etot_q(iq) - etot_q(1))/(sum(mtot_q(:, 1))*sin2theta)
         end if
      end do

      open (newunit=newunit, file=trim(fm_obj%output_file), status='replace', action='write')
      write (newunit, '(A)') '# Frozen-magnon sweep (calculation%post_processing = "frozen_magnon")'
      write (newunit, '(A)') '# q_ss units: Cartesian 2*pi/alat (same convention as &hamiltonian q_ss); row 1 is the reference point'
      write (newunit, '(A,A)') '# q_file coordinates: ', trim(fm_obj%q_coordinates)
      write (newunit, '(A)') '# omega uses eband in mode="mft" and etot in mode="scf"'
      write (newunit, '(A,I0)') '# Number of sublattices (nrec): ', lattice_obj%nrec
      write (newunit, '(A)') '# Format: q1 q2 q3 etot eband mtot_1 .. mtot_nrec omega'
      write (fmt_str, '(A,I0,A)') '(3F12.6,1X,2(ES16.8E3,1X),', lattice_obj%nrec, '(ES16.8E3,1X),ES16.8E3)'
      do iq = 1, fm_obj%n_q
         write (newunit, fmt_str) q_ss_cart(:, iq), etot_q(iq), eband_q(iq), mtot_q(:, iq), omega_q(iq)
      end do
      close (newunit)
   end subroutine post_processing_frozen_magnon

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
          .and. post_processing /= 'density_of_states' &
          .and. post_processing /= 'frozen_magnon') then
         call g_logger%fatal('[calculation.check_post_processing]: '// &
                             "calculation%post_processing must be one of: ''none'', ''paoflow2rs'', ''exchange'', ''exchange_p2rs''," // &
                             " 'conductivity', 'conductivity_p2rs', 'orbital_modern', 'band_structure', 'density_of_states'," // &
                             " 'frozen_magnon'", __FILE__, __LINE__)
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
