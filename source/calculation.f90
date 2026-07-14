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
   use string_mod, only: sl, fmt, real2str, int2str
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
      procedure, private :: post_processing_kspace_green
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
      case ('kspace_green')
         call this%post_processing_kspace_green()
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
   !> @brief B2 validation driver (C1/C3 + C2): cross-check the k-space Lehmann
   !>        on-site Green's function against the real-space recursion route.
   !> @details Runs BOTH routes on the same converged potential and compares the
   !>        on-site density of states rho(E) = -1/pi * Im Tr G_ii(E), built from
   !>        the SAME green%gij on-site block for each route (convention 5: the
   !>        full-resolvent nb x nb block, not the -i*pi*rho spectral form). This
   !>        is the route-agnosticism acceptance for backend E: the strict-Lehmann
   !>        k-space engine must reproduce the recursion on-site spectral function
   !>        (C1 representation) with the B1 bond/phase convention (C3, on-site
   !>        phase = 1). The recursion route fills green%gij at the on-site
   !>        self-pair (ijpair(:,1)==ijpair(:,2)); we snapshot its on-site block,
   !>        then reciprocal%fill_green overwrites green%gij from the k-space
   !>        eigenpairs and we recompute the same quantities. Report-only: writes
   !>        kspace_green_c1.dat with both traces and difference/spectral-weight
   !>        metrics; the acceptance tolerance is a maintainer gate. Regression is
   !>        untouched -- this is a new post_processing key, off by default.
   !>
   !>        C2 (spin frame): the charge DOS -1/pi Im Tr G_ii is invariant under
   !>        G -> R^dagger G R, so it cannot see the spin frame at all. We ALSO
   !>        report the z-projected spin DOS m_z(E) = -1/pi Im[Tr G_uu - Tr G_dd],
   !>        which IS frame-sensitive, as a noncollinear cross-check. The key
   !>        finding (verified on Mn3Sn, 120 deg NC): the INTERSITE recursion
   !>        `recur_b_ij` never rotates to local axes (only the on-site DOS
   !>        recursion does), so the RS gij is stored in the GLOBAL spin frame.
   !>        Backend E therefore fills the global-frame block directly and m_z from
   !>        both routes agrees to ~4e-4 (for an in-plane moment both are ~0). If a
   !>        LOCAL-frame comparison were wanted, BOTH routes' G would have to be
   !>        rotated by the same rotmag_loc primitive -- rotating only one side is
   !>        wrong (it drove m_z diff to ~20). We still restore the global ee
   !>        before the k-space build (rotate_from_local_axis) as a safety net in
   !>        case any on-site recursion left hamiltonian%ee rotated.
   !---------------------------------------------------------------------------
   subroutine post_processing_kspace_green(this)
      use sigma_provider_mod, only: sigma_zero
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
      type(reciprocal), target :: reciprocal_obj
      type(sigma_zero) :: sigma
      real(rp), allocatable :: dos_rs(:), dos_le(:), mz_rs(:), mz_le(:)
      complex(rp), allocatable :: blk_rs(:, :, :), blk_le(:, :, :)
      integer :: i, ie, ne, nblk, norb_h, p0, newunit
      real(rp) :: max_abs_diff, rms_diff, int_rs, int_le, de
      real(rp) :: max_mz_diff, max_blk_diff
      real(rp) :: onsite_mom(3)

      call prepare_post_processing_stack(this, .false., .true., .true., .false., control_obj, lattice_obj, &
                                         charge_obj, mix_obj, energy_obj, hamiltonian_obj, recursion_obj, dos_obj, &
                                         green_obj, bands_obj)
      do i = 1, lattice_obj%ntype
         call lattice_obj%symbolic_atoms(i)%predls(lattice_obj%wav*ang2au)
      end do

      ! --- Real-space recursion route: fill green%gij (on-site = self pair). ---
      call run_intersite_moments(control_obj, recursion_obj)
      call green_obj%calculate_intersite_gf()

      ! Locate an on-site self-pair (ijpair(p,1) == ijpair(p,2)) to read G_ii from.
      ! Prefer one whose moment is OFF the global z axis so the C2 local-frame
      ! rotation is actually exercised (an on-site block has a single, unambiguous
      ! moment -- the intersite i/=j frame question is deliberately out of scope).
      ! Falls back to the first on-site pair (e.g. collinear bcc Fe, all along z).
      p0 = 0
      do i = 1, lattice_obj%njij
         if (lattice_obj%ijpair(i, 1) == lattice_obj%ijpair(i, 2)) then
            if (p0 == 0) p0 = i
            onsite_mom = lattice_obj%symbolic_atoms(i)%potential%mom(1:3)
            if (onsite_mom(1)**2 + onsite_mom(2)**2 > 1.0e-8_rp) then
               p0 = i
               exit
            end if
         end if
      end do
      if (p0 == 0) then
         call g_logger%fatal('[calculation.post_processing_kspace_green]: no on-site self-pair '// &
                             '(ijpair(:,1)==ijpair(:,2)) in the &lattice ijpair list; add e.g. ijpair(1,:)=1,1.', &
                             __FILE__, __LINE__)
      end if

      ne = green_obj%en%channels_ldos + 10
      nblk = size(green_obj%gij, 1)
      norb_h = nblk/2   ! spin-up orbitals 1..norb_h, spin-down norb_h+1..nblk
      allocate (dos_rs(ne), dos_le(ne), mz_rs(ne), mz_le(ne))
      allocate (blk_rs(nblk, nblk, ne), blk_le(nblk, nblk, ne))
      blk_rs = green_obj%gij(:, :, :, p0)
      call onsite_dos_mz(blk_rs, norb_h, dos_rs, mz_rs)

      ! Restore the global-frame ee: the recursion route rotates hamiltonian%ee in
      ! place under local_axis and never calls rotate_from_local_axis, so H(k)
      ! must be rebuilt from the restored global representation (no-op collinear).
      onsite_mom = lattice_obj%symbolic_atoms(p0)%potential%mom(1:3)
      call hamiltonian_obj%rotate_from_local_axis(onsite_mom)

      ! --- k-space Lehmann route: overwrite green%gij from the eigenpairs. ---
      reciprocal_obj = reciprocal(hamiltonian_obj)
      call reciprocal_obj%generate_mp_mesh()   ! full unreduced BZ mesh (backend E requirement)
      call reciprocal_obj%fill_green(green_obj, sigma)
      blk_le = green_obj%gij(:, :, :, p0)
      call onsite_dos_mz(blk_le, norb_h, dos_le, mz_le)

      ! --- Compare + report (report-only; tolerance is a maintainer gate). ---
      max_abs_diff = 0.0_rp; rms_diff = 0.0_rp; int_rs = 0.0_rp; int_le = 0.0_rp
      max_mz_diff = 0.0_rp; max_blk_diff = 0.0_rp
      do ie = 1, ne
         max_abs_diff = max(max_abs_diff, abs(dos_le(ie) - dos_rs(ie)))
         max_mz_diff = max(max_mz_diff, abs(mz_le(ie) - mz_rs(ie)))
         max_blk_diff = max(max_blk_diff, maxval(abs(blk_le(:, :, ie) - blk_rs(:, :, ie))))
         rms_diff = rms_diff + (dos_le(ie) - dos_rs(ie))**2
      end do
      rms_diff = sqrt(rms_diff/real(ne, rp))
      do ie = 2, ne
         de = green_obj%en%ene(ie) - green_obj%en%ene(ie - 1)
         int_rs = int_rs + 0.5_rp*de*(dos_rs(ie) + dos_rs(ie - 1))
         int_le = int_le + 0.5_rp*de*(dos_le(ie) + dos_le(ie - 1))
      end do

      if (rank == 0) then
         open (newunit=newunit, file='kspace_green_c1.dat', status='replace', action='write')
         write (newunit, '(A)') '# B2 C1/C3/C2 validation: on-site G_ii(E), recursion (RS) vs k-space Lehmann'
         write (newunit, '(A)') '# dos = -1/pi Im Tr G_ii ; mz = -1/pi Im (Tr G_uu - Tr G_dd)  (global z; gij stored global frame)'
         write (newunit, '(A,I0,A,I0,A,I0,A,I0,A)') '# k-mesh ', reciprocal_obj%nk_mesh(1), ' x ', &
            reciprocal_obj%nk_mesh(2), ' x ', reciprocal_obj%nk_mesh(3), '  (', reciprocal_obj%nk_total, ' points)'
         write (newunit, '(A,ES14.6,A)') '# green_eta = ', reciprocal_obj%green_eta, ' Ry'
         write (newunit, '(A,3F10.6)') '# on-site moment = ', onsite_mom
         write (newunit, '(A,ES14.6)') '# C1  max|dos_lehmann - dos_rs|      = ', max_abs_diff
         write (newunit, '(A,ES14.6)') '#     rms dos_lehmann - dos_rs       = ', rms_diff
         write (newunit, '(A,ES14.6)') '# C2  max|mz_lehmann - mz_rs|        = ', max_mz_diff
         write (newunit, '(A,ES14.6)') '# C2  max|G_ii^lehmann - G_ii^rs|    = ', max_blk_diff
         write (newunit, '(A,ES14.6,A,ES14.6)') '# spectral weight   RS = ', int_rs, '   Lehmann = ', int_le
         write (newunit, '(A)') '#     E(Ry)          dos_rs        dos_lehmann       dos_diff          mz_rs         mz_lehmann'
         do ie = 1, ne
            write (newunit, '(6ES16.8)') green_obj%en%ene(ie), dos_rs(ie), dos_le(ie), &
               dos_le(ie) - dos_rs(ie), mz_rs(ie), mz_le(ie)
         end do
         close (newunit)
         call g_logger%info('[kspace_green] C1 on-site DOS cross-check: max|diff|='// &
                            trim(real2str(max_abs_diff))//' rms='//trim(real2str(rms_diff))// &
                            ' weight_RS='//trim(real2str(int_rs))//' weight_Lehmann='//trim(real2str(int_le)), &
                            __FILE__, __LINE__)
         call g_logger%info('[kspace_green] C2 spin-structure cross-check (both routes '// &
                            'global frame): max|mz_diff|='//trim(real2str(max_mz_diff))// &
                            ' max|block_diff|='//trim(real2str(max_blk_diff))// &
                            ' (block_diff is single-element ripple at coarse mesh)', &
                            __FILE__, __LINE__)
      end if

      deallocate (dos_rs, dos_le, mz_rs, mz_le, blk_rs, blk_le)
   end subroutine post_processing_kspace_green

   !> @brief On-site charge DOS and z-projected spin DOS from a G_ii block.
   !> @details dos(E)=-1/pi Im Tr G ; mz(E)=-1/pi Im (Tr G_uu - Tr G_dd), with the
   !>          spin-up orbitals in rows 1..norb and spin-down in norb+1..2*norb.
   subroutine onsite_dos_mz(blk, norb, dos, mz)
      complex(rp), intent(in) :: blk(:, :, :)
      integer, intent(in) :: norb
      real(rp), intent(out) :: dos(:), mz(:)
      integer :: ie, ib, ne
      complex(rp) :: tr_up, tr_dn

      ne = size(blk, 3)
      do ie = 1, ne
         tr_up = (0.0_rp, 0.0_rp); tr_dn = (0.0_rp, 0.0_rp)
         do ib = 1, norb
            tr_up = tr_up + blk(ib, ib, ie)
            tr_dn = tr_dn + blk(ib + norb, ib + norb, ie)
         end do
         dos(ie) = -aimag(tr_up + tr_dn)/pi
         mz(ie) = -aimag(tr_up - tr_dn)/pi
      end do
   end subroutine onsite_dos_mz

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
      real(rp), allocatable :: q_ss_cart(:, :)
      real(rp), dimension(3, 3) :: direct_to_cart
      integer :: iq

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

      if (fm_obj%branch_mode == 'auto') then
         call post_processing_frozen_magnon_auto(fm_obj, q_ss_cart, lattice_obj, hamiltonian_obj, bands_obj, &
                                                mix_obj, self_obj, energy_obj)
      else
         call post_processing_frozen_magnon_acoustic(fm_obj, q_ss_cart, lattice_obj, hamiltonian_obj, bands_obj, &
                                                    mix_obj, self_obj, energy_obj)
      end if
   end subroutine post_processing_frozen_magnon

   subroutine post_processing_frozen_magnon_acoustic(fm_obj, q_ss_cart, lattice_obj, hamiltonian_obj, bands_obj, mix_obj, self_obj, energy_obj)
      type(frozen_magnon), intent(in) :: fm_obj
      real(rp), intent(in) :: q_ss_cart(:, :)
      type(lattice), target, intent(inout) :: lattice_obj
      type(hamiltonian), target, intent(inout) :: hamiltonian_obj
      type(bands), target, intent(inout) :: bands_obj
      type(mix), target, intent(inout) :: mix_obj
      type(self), target, intent(inout) :: self_obj
      type(energy), target, intent(inout) :: energy_obj
      type(reciprocal) :: recip_obj
      real(rp), allocatable :: etot_q(:), eband_q(:), mtot_q(:, :), omega_q(:)
      real(rp) :: sin2theta, etot_ref
      logical :: use_kspace
      character(len=200) :: fmt_str
      integer :: iq, i, newunit

      allocate (etot_q(fm_obj%n_q), eband_q(fm_obj%n_q), mtot_q(lattice_obj%nrec, fm_obj%n_q), omega_q(fm_obj%n_q))

      ! Reference point (row 1): converge the flat-spiral cone potential once. Its
      ! magnitudes/moments define the normalization and, in mft mode, are held FIXED
      ! for every q (magnetic force theorem).
      self_obj = self(bands_obj, mix_obj)
      use_kspace = self_obj%use_kspace
      ! k-space runs use the GBT twist (reciprocal module), so ham0m_nc must NOT
      ! rotate the moments by q_ss; real-space runs do the rotation (explicit spiral).
      hamiltonian_obj%gbt_kspace = use_kspace
      hamiltonian_obj%q_ss(:) = q_ss_cart(:, 1)
      call self_obj%run()
      etot_ref = sum(lattice_obj%symbolic_atoms(:)%potential%etot)
      do i = 1, lattice_obj%nrec
         mtot_q(i, :) = lattice_obj%symbolic_atoms(lattice_obj%nbulk + i)%potential%mtot
      end do

      if (fm_obj%mode == 'mft') then
         ! Force theorem: rebuild ONLY the Hamiltonian per q and take the band energy
         ! at the fixed reference potential (frozen_magnon_probe_energy). Running a
         ! single self-consistency step here instead would mix the potential once per
         ! q, letting the reference random-walk across the sweep -- so E(q) would fail
         ! to return to its value at symmetry-equivalent q (e.g. the two Gamma points).
         if (use_kspace) then
            recip_obj = reciprocal(hamiltonian_obj)
            if (.not. allocated(recip_obj%k_points)) then
               if (recip_obj%use_symmetry_reduction) then
                  call recip_obj%generate_reduced_kpoint_mesh(recip_obj%nk_mesh, &
                                                              sum(abs(recip_obj%k_offset)) > 1.0e-12_rp)
               else
                  call recip_obj%generate_mp_mesh()
               end if
            end if
         end if
         do iq = 1, fm_obj%n_q
            hamiltonian_obj%q_ss(:) = q_ss_cart(:, iq)
            eband_q(iq) = frozen_magnon_probe_energy(bands_obj, recip_obj, energy_obj, use_kspace)
            etot_q(iq) = etot_ref   ! potential frozen; total energy not re-evaluated
         end do
      else
         ! scf: fully self-consistent spiral at each q (each q independently relaxed).
         do iq = 1, fm_obj%n_q
            hamiltonian_obj%q_ss(:) = q_ss_cart(:, iq)
            self_obj = self(bands_obj, mix_obj)
            call self_obj%run()
            etot_q(iq) = sum(lattice_obj%symbolic_atoms(:)%potential%etot)
            if (.not. self_obj%use_kspace) call bands_obj%calculate_band_energy()
            eband_q(iq) = bands_obj%eband
            do i = 1, lattice_obj%nrec
               mtot_q(i, iq) = lattice_obj%symbolic_atoms(lattice_obj%nbulk + i)%potential%mtot
            end do
         end do
      end if

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
   end subroutine post_processing_frozen_magnon_acoustic

   !> @brief Multi-sublattice adiabatic magnon branches via the direct GBT frozen-magnon
   !>        method (Essenberger et al., PRB 84, 174425 (2011), Eqs. 24-27; Sandratskii,
   !>        Carva & Silkin, PRB 111, 184436 (2025)).
   !> @details For each q the magnon dynamical matrix is the second derivative of the
   !>        frozen-magnon energy surface with respect to the sublattice cone angles,
   !>        Re[J~_uv^q] = (1/M_u M_v) d^2 E_q/dtheta_u dtheta_v |_0 (Eq. 26), and the
   !>        magnon energies are the eigenvalues of the real symmetric matrix
   !>        sqrt(M_u M_v) Re[J~_uv^q] (Eq. 21). The energy surface is evaluated with the
   !>        magnetic force theorem: the band energy at the FIXED converged reference
   !>        potential, with the moment directions rotated to the spiral/tilt config
   !>        (imposed through the GBT Hamiltonian's q_ss + per-sublattice theta/phi). The
   !>        q-dependence enters through the GBT bond phase already carried by the
   !>        Hamiltonian, so the matrix is real symmetric (no separate azimuthal-phase
   !>        probe) and the Goldstone/acoustic mode is exact at Gamma by construction
   !>        (the Weiss field is intrinsic to the diagonal of the energy surface).
   !>        Works for the k-space (primary) and real-space recursion band-energy paths.
   !>        Assumptions: near-harmonic tilt angle (theta_probe ~ 5-30 deg); Hermitian/real
   !>        treatment appropriate for collinear FM/ferri references without SOC (the
   !>        antisymmetric DMI part is neglected).
   subroutine post_processing_frozen_magnon_auto(fm_obj, q_ss_cart, lattice_obj, hamiltonian_obj, bands_obj, mix_obj, self_obj, energy_obj)
      type(frozen_magnon), intent(in) :: fm_obj
      real(rp), intent(in) :: q_ss_cart(:, :)
      type(lattice), target, intent(inout) :: lattice_obj
      type(hamiltonian), target, intent(inout) :: hamiltonian_obj
      type(bands), target, intent(inout) :: bands_obj
      type(mix), target, intent(inout) :: mix_obj
      type(self), target, intent(inout) :: self_obj
      type(energy), target, intent(inout) :: energy_obj
      type(reciprocal) :: recip_obj
      complex(rp), allocatable :: omega_mat(:, :), eigvec(:, :)
      real(rp), allocatable :: eval(:), mtot_ref(:), ref_theta(:), ref_phi(:), single_energy(:)
      real(rp), allocatable :: active_moment(:)
      integer, allocatable :: active_rec(:), active_type(:)
      real(rp) :: etot_ref, energy_ref, theta_probe, dmz_fac, e_pair, off
      logical :: use_kspace
      integer :: iq, iact, jact, nactive, irec, branch, newunit_branch, newunit_modes

      ! theta_probe: SMALL cone angle for the harmonic finite-difference probes that
      ! build the magnon dynamical matrix. Must stay in the near-harmonic regime
      ! (Sandratskii PRB 111, 184436, Fig. 5: theta ~ 5-30 deg, 20 typical). It is
      ! deliberately decoupled from the flat-spiral theta_ss (typically 90 deg).
      theta_probe = fm_obj%theta_probe
      if (theta_probe <= 0.0_rp) theta_probe = 20.0_rp*pi/180.0_rp
      ! Magnon energy is normalized to omega = Delta E / Delta m_z (Eq. 1 of PRB 111,
      ! 184436). For a rigid tilt of sublattice mu by theta_probe, Delta m_z =
      ! M_mu (1 - cos theta_probe); dmz_fac is the per-unit-moment factor.
      dmz_fac = 1.0_rp - cos(theta_probe)
      if (dmz_fac < 1.0e-10_rp) then
         call g_logger%fatal('[calculation.post_processing_frozen_magnon_auto]: theta_probe must define a finite cone angle', &
                             __FILE__, __LINE__)
      end if

      ! Reference sublattice axes from the collinear ground state (theta=0 up, pi down).
      allocate (ref_theta(lattice_obj%ntype), ref_phi(lattice_obj%ntype), mtot_ref(lattice_obj%nrec))
      call set_reference_sublattice_angles(lattice_obj, ref_theta, ref_phi)

      ! --- Converge the collinear reference potential; it is held FIXED afterwards. ---
      if (allocated(hamiltonian_obj%theta_ss_sublattice)) deallocate (hamiltonian_obj%theta_ss_sublattice)
      if (allocated(hamiltonian_obj%phi_ss_sublattice)) deallocate (hamiltonian_obj%phi_ss_sublattice)
      hamiltonian_obj%theta_ss = 0.0_rp
      hamiltonian_obj%q_ss(:) = 0.0_rp
      self_obj = self(bands_obj, mix_obj)
      hamiltonian_obj%gbt_kspace = self_obj%use_kspace
      call self_obj%run()
      use_kspace = self_obj%use_kspace
      etot_ref = sum(lattice_obj%symbolic_atoms(:)%potential%etot)

      ! Reference moments (fixed) and the set of active magnetic sublattices.
      do irec = 1, lattice_obj%nrec
         mtot_ref(irec) = lattice_obj%symbolic_atoms(lattice_obj%nbulk + irec)%potential%mtot
      end do
      nactive = count(abs(mtot_ref(:)) >= fm_obj%active_moment_threshold)
      if (nactive <= 0) then
         call g_logger%fatal('[calculation.post_processing_frozen_magnon_auto]: no active magnetic sublattices above threshold', &
                             __FILE__, __LINE__)
      end if
      allocate (active_rec(nactive), active_type(nactive), active_moment(nactive))
      iact = 0
      do irec = 1, lattice_obj%nrec
         if (abs(mtot_ref(irec)) >= fm_obj%active_moment_threshold) then
            iact = iact + 1
            active_rec(iact) = irec
            active_type(iact) = lattice_obj%nbulk + irec
            active_moment(iact) = abs(mtot_ref(irec))
         end if
      end do

      ! --- Force-theorem evaluator: build a k-space mesh once if in k-space mode. ---
      if (use_kspace) then
         recip_obj = reciprocal(hamiltonian_obj)
         if (.not. allocated(recip_obj%k_points)) then
            if (recip_obj%use_symmetry_reduction) then
               call recip_obj%generate_reduced_kpoint_mesh(recip_obj%nk_mesh, &
                                                           sum(abs(recip_obj%k_offset)) > 1.0e-12_rp)
            else
               call recip_obj%generate_mp_mesh()
            end if
         end if
      end if

      ! Reference band energy (collinear) via the force theorem at the fixed potential;
      ! q-independent (no transverse components), so evaluated once and reused.
      if (allocated(hamiltonian_obj%theta_ss_sublattice)) deallocate (hamiltonian_obj%theta_ss_sublattice)
      if (allocated(hamiltonian_obj%phi_ss_sublattice)) deallocate (hamiltonian_obj%phi_ss_sublattice)
      hamiltonian_obj%theta_ss = 0.0_rp
      hamiltonian_obj%q_ss(:) = 0.0_rp
      energy_ref = frozen_magnon_probe_energy(bands_obj, recip_obj, energy_obj, use_kspace)

      open (newunit=newunit_branch, file='frozen_magnon_branches.dat', status='replace', action='write')
      open (newunit=newunit_modes, file='frozen_magnon_modes.dat', status='replace', action='write')
      call write_frozen_magnon_auto_headers(newunit_branch, newunit_modes, fm_obj, q_ss_cart, theta_probe, &
                                            active_rec, active_type, active_moment, etot_ref, energy_ref)

      allocate (omega_mat(nactive, nactive), eigvec(nactive, nactive), eval(nactive), single_energy(nactive))
      do iq = 1, fm_obj%n_q
         hamiltonian_obj%q_ss(:) = q_ss_cart(:, iq)
         omega_mat(:, :) = cmplx(0.0_rp, 0.0_rp, kind=rp)

         ! Diagonal: single-sublattice tilt. d^2E/dtheta_i^2 / M_i, in omega = dE/dm_z units.
         do iact = 1, nactive
            call set_probe_angles(hamiltonian_obj, ref_theta, ref_phi, active_type, theta_probe, [iact])
            single_energy(iact) = frozen_magnon_probe_energy(bands_obj, recip_obj, energy_obj, use_kspace)
            omega_mat(iact, iact) = cmplx((single_energy(iact) - energy_ref)/(active_moment(iact)*dmz_fac), &
                                          0.0_rp, kind=rp)
         end do

         ! Off-diagonal: two-sublattice tilt at the NATURAL GBT spiral phase (the q-phase
         ! is carried by the Hamiltonian's bond rotation, so no azimuthal-phase probe is
         ! needed). Cross second derivative d^2E/dtheta_i dtheta_j, real symmetric:
         !   omega_ij = (E_ij - E_i - E_j + E_ref) / [2 (1-cos theta) sqrt(M_i M_j)].
         do iact = 1, nactive - 1
            do jact = iact + 1, nactive
               call set_probe_angles(hamiltonian_obj, ref_theta, ref_phi, active_type, theta_probe, [iact, jact])
               e_pair = frozen_magnon_probe_energy(bands_obj, recip_obj, energy_obj, use_kspace)
               off = (e_pair - single_energy(iact) - single_energy(jact) + energy_ref) / &
                     (2.0_rp*dmz_fac*sqrt(active_moment(iact)*active_moment(jact)))
               omega_mat(iact, jact) = cmplx(off, 0.0_rp, kind=rp)
               omega_mat(jact, iact) = cmplx(off, 0.0_rp, kind=rp)
            end do
         end do

         ! Magnon energies = eigenvalues of the real symmetric matrix (Goldstone exact at
         ! Gamma by construction; no ad-hoc sum-rule correction).
         eigvec(:, :) = omega_mat(:, :)
         call diagonalize_frozen_magnon_matrix(eigvec, eval)
         do branch = 1, nactive
            write (newunit_branch, '(3F12.6,1X,I6,1X,ES16.8E3,1X,I6)') q_ss_cart(:, iq), branch, eval(branch), nactive
            do iact = 1, nactive
               write (newunit_modes, '(3F12.6,1X,I6,1X,I6,1X,I6,1X,ES16.8E3,1X,ES16.8E3,1X,ES16.8E3)') &
                  q_ss_cart(:, iq), branch, iact, active_rec(iact), abs(eigvec(iact, branch)), &
                  atan2(aimag(eigvec(iact, branch)), real(eigvec(iact, branch))), active_moment(iact)
            end do
         end do
      end do
      close (newunit_branch)
      close (newunit_modes)

      if (allocated(hamiltonian_obj%theta_ss_sublattice)) deallocate (hamiltonian_obj%theta_ss_sublattice)
      if (allocated(hamiltonian_obj%phi_ss_sublattice)) deallocate (hamiltonian_obj%phi_ss_sublattice)
   end subroutine post_processing_frozen_magnon_auto

   subroutine set_reference_sublattice_angles(lattice_obj, ref_theta, ref_phi)
      type(lattice), intent(in) :: lattice_obj
      real(rp), intent(out) :: ref_theta(:), ref_phi(:)
      integer :: itype

      ref_phi(:) = 0.0_rp
      do itype = 1, size(ref_theta)
         if (lattice_obj%symbolic_atoms(itype)%potential%mom(3) < 0.0_rp) then
            ref_theta(itype) = pi
         else
            ref_theta(itype) = 0.0_rp
         end if
      end do
   end subroutine set_reference_sublattice_angles

   !> @brief Impose small cone tilts on the selected active sublattices.
   !> @details All sublattices start on their reference axis (theta=0 up / pi down);
   !>        the ones listed in active_probe are tilted by theta_probe toward the xy-plane
   !>        at azimuth ref_phi (the natural relative spiral phase between sublattices is
   !>        supplied by the GBT bond rotation via q_ss, not here).
   subroutine set_probe_angles(hamiltonian_obj, ref_theta, ref_phi, active_type, theta_probe, active_probe)
      type(hamiltonian), intent(inout) :: hamiltonian_obj
      real(rp), intent(in) :: ref_theta(:), ref_phi(:), theta_probe
      integer, intent(in) :: active_type(:), active_probe(:)
      integer :: iprobe, itype

      if (allocated(hamiltonian_obj%theta_ss_sublattice)) deallocate (hamiltonian_obj%theta_ss_sublattice)
      if (allocated(hamiltonian_obj%phi_ss_sublattice)) deallocate (hamiltonian_obj%phi_ss_sublattice)
      allocate (hamiltonian_obj%theta_ss_sublattice(size(ref_theta)))
      allocate (hamiltonian_obj%phi_ss_sublattice(size(ref_phi)))
      hamiltonian_obj%theta_ss_sublattice(:) = ref_theta(:)
      hamiltonian_obj%phi_ss_sublattice(:) = ref_phi(:)

      do iprobe = 1, size(active_probe)
         itype = active_type(active_probe(iprobe))
         if (ref_theta(itype) > 0.5_rp*pi) then
            hamiltonian_obj%theta_ss_sublattice(itype) = pi - theta_probe
            hamiltonian_obj%phi_ss_sublattice(itype) = ref_phi(itype) + pi
         else
            hamiltonian_obj%theta_ss_sublattice(itype) = theta_probe
            hamiltonian_obj%phi_ss_sublattice(itype) = ref_phi(itype)
         end if
      end do
   end subroutine set_probe_angles

   !> @brief Magnetic force theorem: band energy at the FIXED reference potential for the
   !>        current spiral/tilt config (set on the hamiltonian via q_ss + theta/phi_ss_sublattice).
   !> @details Rebuilds only the Hamiltonian (rotated moment directions, fixed potential
   !>        parameters) and evaluates the single-particle band energy through the k-space
   !>        (build/diagonalize/DOS -> band energy from moments) or real-space recursion
   !>        (recursion -> Green -> fermi -> band energy) path. No charge/potential update
   !>        is done, so the reference potential is held fixed across every probe -- this is
   !>        the Liechtenstein/MFT energy surface whose second derivatives give the magnon
   !>        matrix (Essenberger PRB 84, 174425, Eq. 26).
   function frozen_magnon_probe_energy(bands_obj, recip_obj, energy_obj, use_kspace) result(e_band)
      type(bands), target, intent(inout) :: bands_obj
      type(reciprocal), intent(inout) :: recip_obj
      type(energy), intent(in) :: energy_obj
      logical, intent(in) :: use_kspace
      real(rp) :: e_band

      associate (ham => bands_obj%recursion%hamiltonian, rec => bands_obj%recursion, &
                 grn => bands_obj%green, ctl => bands_obj%lattice%control)
         if (ctl%nsp == 2 .or. ctl%nsp == 4) call ham%build_lsham()
         call ham%build_bulkham()
         if (use_kspace) then
            call recip_obj%build_kspace_hamiltonian()
            call recip_obj%diagonalize_hamiltonian()
            call recip_obj%calculate_density_of_states(ham, n_energy_points=energy_obj%channels_ldos + 10, &
                 energy_range=[energy_obj%energy_min, energy_obj%energy_max], fermi_level=energy_obj%fermi, &
                 auto_find_fermi=.true.)
            e_band = recip_obj%calculate_band_energy_from_moments()
         else
            select case (ctl%recur)
            case ('block')
               call rec%recur_b()
               call rec%zsqr()
               call grn%block_green()
            case ('chebyshev')
               call rec%chebyshev_recur()
               call grn%chebyshev_dos_dispatch()
            case ('lanczos')
               call rec%recur()
               call grn%sgreen()
            end select
            call bands_obj%calculate_fermi()
            call bands_obj%calculate_band_energy()
            e_band = bands_obj%eband
         end if
      end associate
   end function frozen_magnon_probe_energy

   subroutine diagonalize_frozen_magnon_matrix(mat, eval)
      complex(rp), intent(inout) :: mat(:, :)
      real(rp), intent(out) :: eval(:)
      complex(rp), allocatable :: work(:)
      real(rp), allocatable :: rwork(:)
      integer :: n, lwork, info
      external :: zheev

      n = size(mat, 1)
      if (size(mat, 2) /= n .or. size(eval) /= n) then
         call g_logger%fatal('[calculation.diagonalize_frozen_magnon_matrix]: inconsistent matrix dimensions', &
                             __FILE__, __LINE__)
      end if
      lwork = max(1, 2*n - 1)
      allocate (work(lwork), rwork(max(1, 3*n - 2)))
      call zheev('V', 'U', n, mat, n, eval, work, lwork, rwork, info)
      if (info /= 0) then
         call g_logger%fatal('[calculation.diagonalize_frozen_magnon_matrix]: zheev failed with info='// &
                             int2str(info), __FILE__, __LINE__)
      end if
      deallocate (work, rwork)
   end subroutine diagonalize_frozen_magnon_matrix

   subroutine write_frozen_magnon_auto_headers(branch_unit, modes_unit, fm_obj, q_ss_cart, theta_probe, active_rec, &
                                               active_type, active_moment, etot_ref, eband_ref)
      integer, intent(in) :: branch_unit, modes_unit
      type(frozen_magnon), intent(in) :: fm_obj
      real(rp), intent(in) :: q_ss_cart(:, :), theta_probe, active_moment(:), etot_ref, eband_ref
      integer, intent(in) :: active_rec(:), active_type(:)
      integer :: i

      write (branch_unit, '(A)') '# Frozen-magnon auto branches (calculation%post_processing = "frozen_magnon")'
      write (branch_unit, '(A)') '# branch_mode: auto'
      write (branch_unit, '(A,A)') '# q_file coordinates: ', trim(fm_obj%q_coordinates)
      write (branch_unit, '(A)') '# q_ss units: Cartesian 2*pi/alat'
      write (branch_unit, '(A,ES16.8E3)') '# theta_probe(rad): ', theta_probe
      write (branch_unit, '(A,ES16.8E3)') '# active_moment_threshold: ', fm_obj%active_moment_threshold
      write (branch_unit, '(A,ES16.8E3)') '# reference_etot: ', etot_ref
      write (branch_unit, '(A,ES16.8E3)') '# reference_eband: ', eband_ref
      write (branch_unit, '(A,I0)') '# active_sublattices: ', size(active_rec)
      do i = 1, size(active_rec)
         write (branch_unit, '(A,I0,A,I0,A,I0,A,ES16.8E3)') '# active ', i, ' rec=', active_rec(i), &
            ' type=', active_type(i), ' moment=', active_moment(i)
      end do
      write (branch_unit, '(A)') '# Format: q1 q2 q3 branch omega nactive'

      write (modes_unit, '(A)') '# Frozen-magnon auto branch eigenvectors'
      write (modes_unit, '(A)') '# branch_mode: auto'
      write (modes_unit, '(A,A)') '# q_file coordinates: ', trim(fm_obj%q_coordinates)
      write (modes_unit, '(A)') '# q_ss units: Cartesian 2*pi/alat'
      write (modes_unit, '(A,ES16.8E3)') '# theta_probe(rad): ', theta_probe
      write (modes_unit, '(A,ES16.8E3)') '# active_moment_threshold: ', fm_obj%active_moment_threshold
      do i = 1, size(active_rec)
         write (modes_unit, '(A,I0,A,I0,A,I0,A,ES16.8E3)') '# active ', i, ' rec=', active_rec(i), &
            ' type=', active_type(i), ' moment=', active_moment(i)
      end do
      write (modes_unit, '(A)') '# Format: q1 q2 q3 branch active_index sublattice_index amplitude phase moment'
   end subroutine write_frozen_magnon_auto_headers

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
          .and. post_processing /= 'kspace_green' &
          .and. post_processing /= 'frozen_magnon') then
         call g_logger%fatal('[calculation.check_post_processing]: '// &
                             "calculation%post_processing must be one of: ''none'', ''paoflow2rs'', ''exchange'', ''exchange_p2rs''," // &
                             " 'conductivity', 'conductivity_p2rs', 'orbital_modern', 'band_structure', 'density_of_states'," // &
                             " 'kspace_green', 'frozen_magnon'", __FILE__, __LINE__)
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
