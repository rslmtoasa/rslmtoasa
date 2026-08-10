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
   use vacuum_lead_mod, only: vacuum_lead, refresh_vacuum_region
   use math_mod
   use precision_mod, only: rp
   use string_mod, only: sl, fmt, real2str, int2str
   use timer_mod, only: g_timer
   use logger_mod, only: g_logger
   use basis_mod, only: basis_init, norb
   use magnetic_representation_mod, only: gbt_single_q
   use tddft_config_mod, only: tddft_config
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, build_chi_ks_from_eigenpairs, &
      build_static_chi_ks_from_eigenpairs, write_chi_ks_text
   use tddft_xi_mod, only: tddft_direct_xi_result, build_direct_xi_from_k_dependent_eigenpairs, &
      build_static_direct_xi_from_k_dependent_eigenpairs
   use tddft_chi0_green_mod, only: green_chi0_options, eigenpair_green_function_provider, &
      build_chi_ks_from_green_functions, build_four_component_chi_ks_from_green_functions
   use response_components_mod, only: RESPONSE_PLUS, RESPONSE_MINUS, RESPONSE_MZ
   use response_vertices_mod, only: response_channel
   use tddft_four_component_mod, only: build_four_component_chi_ks, build_four_component_kernel, &
      evaluate_four_component_zero_modes, tddft_four_component_zero_mode_diagnostics
   use tddft_dyson_mod, only: tddft_dyson_options, tddft_dyson_result, enhance_tddft_susceptibility, &
      enhance_tddft_susceptibility_from_xi, write_tddft_dyson_text
   use tddft_goldstone_mod, only: tddft_goldstone_options, tddft_goldstone_result, evaluate_goldstone, &
      tddft_goldstone_diagnostics, tddft_goldstone_column_correction, evaluate_raw_xi_diagnostics, &
      build_goldstone_column_correction, rescale_pair_potential_columns, spectral_weights_are_nonnegative, &
      write_goldstone_diagnostics_text, append_goldstone_column_correction_text
   use tddft_modes_mod, only: tddft_mode_options, tddft_mode_result, analyze_tddft_modes, write_tddft_modes_text
   use tddft_longitudinal_mod, only: tddft_longitudinal_options, tddft_longitudinal_static_result, &
      tddft_longitudinal_result, read_longitudinal_static_fields, build_longitudinal_static_response, &
      build_longitudinal_kernel, calibrate_longitudinal_response, write_longitudinal_report
#ifdef USE_MPI
   use mpi
#endif
   implicit none

   private

   !---------------------------------------------------------------------------
   ! B7.6 -- state backing the vacuum-refresh hook.
   !
   ! `charge%on_alignment_updated` is argument-free (see the abstract interface
   ! in charge.f90 for why), so the data the refresh needs is held here and the
   ! hook procedure closes over it. Module scope rather than procedure scope
   ! because the hook outlives `pre_processing_buildinterface`'s stack frame:
   ! it fires from inside `self%run()`, per SCF iteration.
   !
   ! Only ever populated on the `buildinterface` path with region_b_kind =
   ! 'vacuum'; a run that is not A | vacuum never installs the hook and never
   ! touches these.
   !---------------------------------------------------------------------------
   type(vacuum_lead), save :: vacuum_state
   type(self), pointer :: vacuum_solver => null()
   type(charge), pointer :: vacuum_charge => null()
   class(symbolic_atom), dimension(:), pointer :: vacuum_atoms => null()
   integer, save :: vacuum_nbulk_a = 0
   integer, save :: vacuum_nbulk = 0
   type(energy), pointer :: vacuum_energy => null()

   type, public :: calculation
      !> Pre-processing. Options are:
      !> ´none´ (default)
      !> ´bravais´ : Builds the bulk clust
      !> ´buildsurf´ : Builds the surface clust
      !> ´newclubulk´ : Builds the imputiry clust from the bluk clust
      !> ´newclusurf´ : Builds the impurity clust from the surface clust
      !> ´buildinterface´ : Builds the two-sided (region A | active | region B)
      !>   layered/interface clust (calctype='L', B7.5)
      character(len=sl) :: pre_processing

      !> Processing. Options are
      !> ´none´ (default)
      character(len=sl) :: processing

      !> Post-processing. Options are
      !> ´none´ (default)
      character(len=sl) :: post_processing

      !> Green's-function route for the intersite G_ij consumers (B2.5 dispatch).
      !> Selects HOW green%gij (+ the torque/eta families) is filled before an
      !> intersite consumer (exchange) runs -- the same canonical arrays, either
      !> route, so consumers are untouched by construction. Options:
      !> ´recursion´ (default) : real-space recursion (bit-identical legacy path)
      !> ´lehmann´             : k-space eigenpair / strict-Lehmann engine (Sigma=0)
      !> ´dyson´               : k-space direct Dyson inversion (Sigma via provider)
      !> The k-space routes read their mesh/eta from the &reciprocal namelist and
      !> override its green_backend with this key.
      character(len=sl) :: gf_route

      !> Gilbert-damping switch for the exchange post-processing (B5.3).
      !> When .true., `post_processing_exchange` also evaluates the on-site
      !> Kamberský torque-correlation Gilbert damping (`calculate_gilbert_damping`)
      !> from the same `green%gij` the J_ij consumer reads -- so it runs through
      !> whichever `gf_route` filled those arrays (recursion / lehmann / dyson),
      !> route-agnostic by construction. Requires SOC in the potential (xi_p/xi_d)
      !> and, for the k-space routes, kspace_ham_order='second'. Default .false.
      !> (bit-identical legacy: the damping routine is not invoked).
      logical :: do_damping

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
      procedure, private :: pre_processing_buildinterface
      procedure, private :: processing_sd
      procedure, private :: post_processing_paoflow2rs
      procedure, private :: post_processing_exchange
      procedure, private :: post_processing_exchange_p2rs
      procedure, private :: post_processing_conductivity_p2rs
      procedure, private :: post_processing_conductivity
      procedure, private :: post_processing_orbital_modern
      procedure, private :: post_processing_band_structure
      procedure, private :: post_processing_bsf
      procedure, private :: post_processing_density_of_states
      procedure, private :: post_processing_kspace_green
      procedure, private :: post_processing_frozen_magnon
      procedure, private :: post_processing_susceptibility
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
      gf_route = this%gf_route
      do_damping = this%do_damping

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
      ! Green's-function route (B2.5 dispatch)
      call check_gf_route(trim(gf_route))

      this%verbose = verbose
      this%fname = fname
      this%pre_processing = pre_processing
      this%processing = processing
      this%post_processing = post_processing
      this%gf_route = gf_route
      this%do_damping = do_damping

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
      case ('buildinterface')
         call this%pre_processing_buildinterface()
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
      case ('bsf')
         call this%post_processing_bsf()
      case ('density_of_states')
         call this%post_processing_density_of_states()
      case ('kspace_green')
         call this%post_processing_kspace_green()
      case ('frozen_magnon')
         call this%post_processing_frozen_magnon()
      case ('susceptibility')
         call this%post_processing_susceptibility()
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
   !> Pre-process for the layered/interface calculation (B7.5): region A |
   !> active zone | region B, calctype='L'. Mirrors pre_processing_buildsurf;
   !> buildsurf itself is untouched and remains the permanent regression
   !> oracle for the one-sided (vacuum|active|bulk) case (B7 §4 B7.5).
   !>
   !> Workflow chain: bulk-A -> bulk-B (-> vacuum generator) -> interface,
   !> mirroring bulk -> surf -> imp. Region A's and region B's converged
   !> parameter sets (potential%vmad included, CONTRACT_FROZEN_REGION.md) are
   !> loaded together through the existing &atoms database=/label(:)
   !> mechanism -- both sets of _out.nml files in one working directory,
   !> exactly as newclusurf/newclubulk already combine host + impurity
   !> labels. No parameter-set version stamp is added anywhere (G-B7-2,
   !> CONTRACT_FROZEN_REGION.md §4.1) -- old *_out.nml files keep working.
   !---------------------------------------------------------------------------
   subroutine pre_processing_buildinterface(this)
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
      call lattice_obj%build_interface_full()
      call lattice_obj%structb(.true.)

      ! Creating the symbolic_atom object
      call lattice_obj%atomlist()

      ! Initialize basis dimension parameters from lmax
      call basis_init(lattice_obj%symbolic_atoms(1)%potential%lmax)

      ! Initializing MPI lookup tables and info.
      call get_mpi_variables(rank, lattice_obj%nrec)

      ! Constructing the charge object
      charge_obj = charge(lattice_obj)
      ! Region reference charges (B7 §1.4, §2.4): bulk_charge per type, from
      ! the loaded symbolic atoms' occupation vs valence -- imppot's
      ! definition, generalized to two regions instead of one host.
      call charge_obj%get_charge_transf
      call charge_obj%build_alelay
      call charge_obj%surfmat
      ! surfmat (reused unchanged, B7 §2.10) ends by building a ONE-SIDED
      ! (vacuum|active|bulk) registry via build_region_registry. Overwrite it
      ! with the genuinely two-sided registry (B7 §1.2, §4 B7.1/B7.5) before
      ! anything else reads charge_obj%regions.
      call charge_obj%build_interface_registry()
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

      ! B7.6, A | vacuum: region B's frozen parameters are GENERATED, not read
      ! from an &atoms label. Install the refresh hook and run it once now, so
      ! the very first recursion sees real empty-lattice parameters rather than
      ! the defaults the unlabelled type was constructed with. The initial call
      ! runs at region_shift = 0, i.e. the vacuum level is the anchor's zero;
      ! every later call, driven from interfacepot, uses the solved level. That
      ! is the "generate once" behaviour arriving as iteration 0 of the
      ! self-consistent scheme rather than as a separate code path.
      if (trim(lattice_obj%region_b_kind) == 'vacuum') then
         vacuum_solver => self_obj
         vacuum_charge => charge_obj
         vacuum_atoms => lattice_obj%symbolic_atoms
         vacuum_energy => energy_obj
         vacuum_nbulk_a = lattice_obj%nbulk_bulk
         vacuum_nbulk = lattice_obj%nbulk
         vacuum_state = vacuum_lead(lattice_obj%symbolic_atoms(1)%potential%lmax)

         charge_obj%on_alignment_updated => vacuum_refresh_hook
         call vacuum_refresh_hook()
         if (rank == 0) then
            call g_logger%info('buildinterface: region B is VACUUM. Its frozen parameters '// &
                               'are generated per run by vacuum_lead and regenerated each '// &
                               'iteration at the solved vacuum level (B7 §1.6).', &
                               __FILE__, __LINE__)
            call vacuum_state%dump()
         end if
      end if

      call g_timer%start('self-consistency')
      call self_obj%run()
      call g_timer%stop('self-consistency')

      call save_state(lattice_obj%symbolic_atoms)

      call self_obj%report()
   end subroutine pre_processing_buildinterface

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> B7.6: the vacuum-refresh hook. Installed on
   !> `charge%on_alignment_updated` for A | vacuum runs and fired by
   !> `interfacepot` once per SCF iteration, after the alignment solver has
   !> updated the vacuum level.
   !>
   !> Argument-free by the hook's interface contract; the data it needs is the
   !> module state above, populated at installation time.
   !---------------------------------------------------------------------------
   subroutine vacuum_refresh_hook()
      if (.not. associated(vacuum_solver)) return
      if (.not. associated(vacuum_charge)) return
      if (.not. associated(vacuum_atoms)) return

      call refresh_vacuum_region(vacuum_state, vacuum_solver, vacuum_charge, &
                                 vacuum_atoms, vacuum_nbulk_a, vacuum_nbulk, &
                                 vacuum_energy%fermi)
   end subroutine vacuum_refresh_hook

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

   !> @brief Build the shared object stack for a post-processing driver.
   !> @details The routine builds the pre-processing stage for calctype B, S, I,
   !>          or L. If use_paoflow is true, it builds a PAOFLOW-imported
   !>          Hamiltonian instead. It also builds the potential and the
   !>          Hamiltonian. It does not run the recursion step or the
   !>          Green's-function step; the caller runs those after this routine
   !>          returns. Most post_processing_* drivers call this routine first
   !>          (post_processing_orbital_modern, post_processing_band_structure,
   !>          post_processing_bsf, and post_processing_density_of_states build
   !>          their own stack instead).
   !> @param[in] this Calculation object. Only fname is read.
   !> @param[in] use_paoflow True builds a PAOFLOW Hamiltonian. False builds the
   !>            normal RS-LMTO-ASA Hamiltonian for the current calctype.
   !> @param[in] use_exchange_pairs True sizes MPI ranks over lattice%njij (atom
   !>            pairs). False sizes them over lattice%ntype (atom types).
   !> @param[in] energy_mesh_before_hamiltonian True builds the energy mesh
   !>            before the Hamiltonian. False builds it after.
   !> @param[in] stochastic_moments True runs compute_moments_stochastic once
   !>            the recursion object exists.
   !> @param[out] control_obj Control object, built from this%fname.
   !> @param[out] lattice_obj Lattice object, built for the current calctype.
   !> @param[out] charge_obj Charge object, with the calctype-specific Madelung
   !>             matrices built.
   !> @param[out] mix_obj Mixing object for the caller's self-consistency step.
   !> @param[out] energy_obj Energy object; the energy mesh may already be built.
   !> @param[out] hamiltonian_obj Hamiltonian object, with the potential and the
   !>             bulk (and, for calctype I, local) blocks built.
   !> @param[out] recursion_obj Recursion object, ready for a recur/recur_b/
   !>             chebyshev_recur call.
   !> @param[out] dos_obj Density-of-states object built on recursion_obj.
   !> @param[out] green_obj Green object built on dos_obj.
   !> @param[out] bands_obj Bands object built on green_obj.
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
         case ('L')
            ! B7.5: build_interface_full instead of build_surf_full -- the
            ! two-sided counterpart. Mirrors pre_processing_buildinterface.
            call lattice_obj%build_data()
            call lattice_obj%bravais()
            call lattice_obj%build_interface_full()
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
         case ('L')
            ! B7.5: same Madelung matrices as 'S' (no interfacemat exists, by
            ! design), plus the region reference charges and the genuinely
            ! two-sided registry that overwrites surfmat's one-sided one.
            ! Mirrors pre_processing_buildinterface.
            call charge_obj%get_charge_transf
            call charge_obj%build_alelay
            call charge_obj%surfmat
            call charge_obj%build_interface_registry()
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
         case ('L')
            call g_logger%fatal('Layered/interface calculation not implemented!', __FILE__, __LINE__)
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
         case ('L')
            ! B7.5: identical to 'S' -- see the same clause in self.f90's
            ! run_recursion for why the loop runs to ntype and there is no
            ! build_locham.
            do i = 1, lattice_obj%ntype
               call lattice_obj%symbolic_atoms(i)%build_pot()
            end do
            if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham
            call hamiltonian_obj%build_bulkham()
         end select
      end if

      recursion_obj = recursion(hamiltonian_obj, energy_obj, sparse(hamiltonian_obj))
      if (stochastic_moments) call recursion_obj%compute_moments_stochastic()
      dos_obj = dos(recursion_obj, energy_obj)
      green_obj = green(dos_obj)
      bands_obj = bands(green_obj)
      if (.not. energy_mesh_before_hamiltonian) call energy_obj%e_mesh()
   end subroutine prepare_post_processing_stack

   !> @brief Run the intersite recursion pass that produces the moments for G_ij.
   !> @details The routine reads control_obj%recur. For 'block' it calls
   !>          recur_b_ij. For 'chebyshev' it calls chebyshev_recur_ij. The
   !>          lanczos route has no intersite form, so the routine does nothing
   !>          for it. Every exchange and kspace_green driver calls this routine
   !>          before it fills green%gij through the real-space recursion route.
   !> @param[in] control_obj Control object. Only recur is read.
   !> @param[inout] recursion_obj Recursion object. The routine fills its
   !>               pair-local moment arrays (a_b/b2_b or mu_n).
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

   !> @brief Finish the Chebyshev DOS step and resolve the Fermi level.
   !> @details The routine calls green%chebyshev_dos_dispatch to turn the
   !>          Chebyshev transport moments into a Green's function. It then
   !>          calls bands%calculate_fermi. Both conductivity drivers
   !>          (post_processing_conductivity and post_processing_conductivity_p2rs)
   !>          call this routine right before they build the conductivity tensor.
   !> @param[inout] green_obj Green object. The Chebyshev DOS dispatch fills its
   !>               Green's-function arrays.
   !> @param[inout] bands_obj Bands object. calculate_fermi sets its Fermi level.
   subroutine finish_conductivity_moments(green_obj, bands_obj)
      type(green), intent(inout) :: green_obj
      type(bands), intent(inout) :: bands_obj

      call green_obj%chebyshev_dos_dispatch()
      call bands_obj%calculate_fermi()
   end subroutine finish_conductivity_moments

   !> @brief Exchange (J_ij) post-processing for a PAOFLOW-imported Hamiltonian.
   !> @details The routine builds the post-processing stack with a PAOFLOW
   !>          Hamiltonian, runs the intersite recursion pass, fills the
   !>          eta-broadened intersite Green's function, and evaluates J_ij
   !>          with the Gauss-Legendre exchange integrator. It always uses the
   !>          real-space recursion route; it has no gf_route dispatch.
   !> @param[in] this Calculation object. fname selects the namelist input.
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

   !> @brief Exchange (J_ij) post-processing, with an optional Gilbert-damping
   !>        evaluation (B2.5/B5.3 dispatch).
   !> @details The routine builds the post-processing stack, then fills
   !>          green%gij (and the torque families) through the route named in
   !>          this%gf_route. Route 'recursion' runs the real-space intersite
   !>          recursion pass. Routes 'lehmann' and 'dyson' fill the same
   !>          arrays from the k-space engine (reciprocal%fill_green). Every
   !>          route fills the same canonical arrays, so
   !>          calculate_intersite_gf_twoindex and the exchange evaluation run
   !>          unchanged after this point. When this%do_damping is true, the
   !>          routine also evaluates the on-site Kambersky torque-correlation
   !>          Gilbert damping from the same green%gij.
   !> @param[in] this Calculation object. Reads fname, gf_route, and do_damping.
   !> @note The damping term needs SOC in the potential (xi_p/xi_d). The
   !>       k-space routes also need kspace_ham_order='second'.
   subroutine post_processing_exchange(this)
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
      type(exchange), target :: exchange_obj
      type(reciprocal), target :: reciprocal_obj
      type(sigma_zero) :: sigma
      integer :: i

      call prepare_post_processing_stack(this, .false., .true., .true., .false., control_obj, lattice_obj, &
                                         charge_obj, mix_obj, energy_obj, hamiltonian_obj, recursion_obj, dos_obj, &
                                         green_obj, bands_obj)
      exchange_obj = exchange(bands_obj)

      do i = 1, lattice_obj%ntype
         call lattice_obj%symbolic_atoms(i)%predls(lattice_obj%wav*ang2au)
      end do

      ! Fill green%gij (+ the torque families) via the selected route (B2.5
      ! dispatch): the real-space recursion route (default, bit-identical legacy
      ! path) or the k-space engine (reciprocal%fill_green, backend E/D). Both
      ! populate the SAME canonical arrays the consumers below read, so the
      ! downstream code is untouched by construction. The k-space routes read
      ! their mesh/eta from the &reciprocal namelist and override its
      ! green_backend with gf_route; they fill the global-frame blocks directly
      ! (C2 resolution: the intersite recursion never rotates to local axes), so
      ! calculate_intersite_gf_twoindex and the exchange consumers run unchanged.
      ! reciprocal_obj is kept at subroutine scope (not a nested helper) so its
      ! spglib-owning finalizer runs exactly as in post_processing_kspace_green.
      ! The zero-consumer-change exchange/damping acceptance is task B2.6.
      select case (trim(this%gf_route))
      case ('recursion')
         call run_intersite_moments(control_obj, recursion_obj)
         call green_obj%calculate_intersite_gf()
      case ('lehmann', 'dyson')
         ! B2.6 (done): the exchange consumer runs UNCHANGED on the k-space-filled
         ! arrays -- it reads the same gij/gji + intersite torque families
         ! (Ginmag/Gj{x,y,z}) the recursion route fills, so calculate_exchange
         ! produces a physical J_ij (correct tensor structure: isotropic J, zero
         ! DMI on collinear bcc Fe). The intersite normalization is CORRECT -- the
         ! former "1/sqrt2" worry is resolved: the fixed-broadening J difference vs
         ! the recursion route is a broadening / metallic-Fermi-surface k-convergence
         ! artifact (shell- and eta-dependent, not a global factor), NOT a
         ! normalization error. The kernel is pinned at machine precision by
         ! UnitGammaSupercell (intersite block == direct resolvent, <1e-12). See
         ! docs/dev/reciprocal_green_convergence.md for the J vs N_k / eta study
         ! (gate G-B2-2).
         call g_logger%info('[calculation.post_processing_exchange]: gf_route='// &
                            trim(this%gf_route)//' -- filling green%gij from the k-space '// &
                            'engine (reciprocal%fill_green); exchange runs unchanged on the '// &
                            'k-space-filled arrays. See docs/dev/reciprocal_green_convergence.md.', &
                            __FILE__, __LINE__)
         reciprocal_obj = reciprocal(hamiltonian_obj)
         reciprocal_obj%green_backend = trim(this%gf_route)
         call reciprocal_obj%generate_mp_mesh()   ! full unreduced BZ (backend E requirement)
         call reciprocal_obj%fill_green(green_obj, sigma)
      end select
      call green_obj%calculate_intersite_gf_twoindex()
      if ((lattice_obj%njij .ne. 0) .and. (lattice_obj%njijk .eq. 0)) then
         call exchange_obj%calculate_exchange()
         call exchange_obj%calculate_exchange_twoindex()
      end if
      ! B5.3: on-site Kamberský torque-correlation Gilbert damping, off by
      ! default. Consumes the same green%gij just filled by gf_route, so it is
      ! route-agnostic (recursion / lehmann / dyson) by construction. Needs SOC
      ! in the potential (xi_p/xi_d); the k-space routes additionally need
      ! kspace_ham_order='second'.
      if (this%do_damping) then
         call exchange_obj%calculate_gilbert_damping()
      end if
   end subroutine post_processing_exchange


   !> @brief Conductivity-tensor post-processing (B5.1 dispatch).
   !> @details The routine picks the source of the Chebyshev transport moments
   !>          (mu_nm_stochastic) from this%gf_route. Route 'recursion' runs
   !>          the stochastic moment generator inside prepare_post_processing_stack.
   !>          Routes 'lehmann' and 'dyson' fill the same array exactly from the
   !>          k-space eigenpairs (reciprocal%fill_moments), on the same
   !>          Chebyshev scaling window the recursion route uses
   !>          (resolve_chebyshev_window). Either way, calculate_gamma_nm and
   !>          calculate_conductivity_tensor then run unchanged on
   !>          mu_nm_stochastic.
   !> @param[in] this Calculation object. Reads fname and gf_route.
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
      type(reciprocal), target :: reciprocal_obj
      logical :: rec_moments
      real(rp) :: emin_win, emax_win, a_scale, b_scale

      ! Route the Chebyshev transport moments (mu_nm_stochastic) through the
      ! selected producer (B5.1 dispatch). 'recursion' (default) fills them
      ! stochastically inside the prepared stack -- the bit-identical legacy
      ! path. 'lehmann'/'dyson' fill the SAME array EXACTLY from the k-space
      ! eigenpairs (reciprocal%fill_moments), so calculate_conductivity_tensor
      ! runs unchanged. exact-vs-recursion on the same crystal is the direct KPM
      ! error bound (see docs/dev/route_agnostic_estimators.md).
      rec_moments = (trim(this%gf_route) == 'recursion')
      call prepare_post_processing_stack(this, .false., .false., .true., rec_moments, control_obj, lattice_obj, &
                                         charge_obj, mix_obj, energy_obj, hamiltonian_obj, recursion_obj, dos_obj, &
                                         green_obj, bands_obj)
      if (.not. rec_moments) then
         ! Exact k-space moments on the SAME Chebyshev window the recursion route
         ! and gamma_nm use (resolve_chebyshev_window -> a, b). The reciprocal
         ! object is kept at subroutine scope (spglib-owning finalizer) like the
         ! exchange/BSF drivers. Sigma=0: 'lehmann' and 'dyson' coincide here.
         call g_logger%info('[calculation.post_processing_conductivity]: gf_route='// &
                            trim(this%gf_route)//' -- filling mu_nm_stochastic from the '// &
                            'exact k-space moment generator (reciprocal%fill_moments); '// &
                            'conductivity runs unchanged on the k-space-filled moments.', &
                            __FILE__, __LINE__)
         call recursion_obj%resolve_chebyshev_window(emin_win, emax_win)
         a_scale = (emax_win - emin_win)/(2.0_rp - 0.3_rp)
         b_scale = (emax_win + emin_win)/2.0_rp
         reciprocal_obj = reciprocal(hamiltonian_obj)
         reciprocal_obj%green_backend = trim(this%gf_route)
         call reciprocal_obj%generate_mp_mesh()   ! full unreduced BZ
         call reciprocal_obj%fill_moments(recursion_obj%mu_nm_stochastic, a_scale, b_scale)
      end if
      self_obj = self(bands_obj, mix_obj)
      call finish_conductivity_moments(green_obj, bands_obj)
      conductivity_obj = conductivity(self_obj)
      call conductivity_obj%calculate_gamma_nm()
      call conductivity_obj%calculate_conductivity_tensor()
   end subroutine post_processing_conductivity


   !> @brief Conductivity-tensor post-processing for a PAOFLOW-imported
   !>        Hamiltonian.
   !> @details The routine builds the post-processing stack with a PAOFLOW
   !>          Hamiltonian and always runs the real-space stochastic moment
   !>          generator. There is no gf_route dispatch here: a PAOFLOW
   !>          Hamiltonian has no k-space reciprocal build.
   !> @param[in] this Calculation object. fname selects the namelist input.
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


   !> @brief Orbital-moment post-processing entry point.
   !> @details The routine builds its own control/lattice/charge/hamiltonian/
   !>          recursion/green/bands/mix stack for the current calctype (B, S,
   !>          I, or L). It predates prepare_post_processing_stack and is not
   !>          routed through it. It then calls
   !>          recursion%chebyshev_orbital_mod, which is always the Chebyshev
   !>          route and loops over every cluster atom rather than only the
   !>          recursion atoms.
   !> @param[in] this Calculation object. fname selects the namelist input.
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
      case ('L')
         ! B7.5: two-sided counterpart of 'S'.
         call lattice_obj%build_data()
         call lattice_obj%bravais()
         call lattice_obj%build_interface_full()
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
      case ('L')
         ! B7.5: 'S' matrices plus region reference charges and the two-sided
         ! registry overwriting surfmat's one-sided one.
         call charge_obj%get_charge_transf
         call charge_obj%build_alelay
         call charge_obj%surfmat
         call charge_obj%build_interface_registry()
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
      case ('L')
         ! B7.5: identical to 'S'; no build_locham (nmax = 0 for 'L').
         do i = 1, lattice_obj%ntype
            call lattice_obj%symbolic_atoms(i)%build_pot() ! Build the potential matrix
         end do
         if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham ! Calculate the spin-orbit coupling Hamiltonian
         call hamiltonian_obj%build_bulkham() ! Build the bulk Hamiltonian for the interface
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

   !> @brief k-space band-structure post-processing entry point.
   !> @details The routine builds a bulk-only pre-processing stack directly
   !>          (bravais, structb, atomlist, bulkmat, potential build,
   !>          Hamiltonian build). It does not support surface, impurity, or
   !>          interface calctypes. It then calls
   !>          reciprocal%calculate_band_structure along an automatic
   !>          high-symmetry path and writes band_structure.dat.
   !> @param[in] this Calculation object. fname selects the namelist input.
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

   !> @brief Bloch spectral function A(k,E) post-processing (milestone B3).
   !> @details Builds the same stack as post_processing_band_structure (the BSF is a
   !>          thin consumer of the B2 k-space Green's function on the band path),
   !>          then delegates to reciprocal%calculate_bsf. Broadening is the
   !>          &reciprocal green_eta; the energy grid/range are n_energy_points /
   !>          dos_energy_range; the path density is nk_per_segment. With sigma=0 the
   !>          resolvent equals backend E; a non-zero Sigma provider (B8/B10) would
   !>          broaden A(k,E) through the same routine unchanged.
   subroutine post_processing_bsf(this)
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
      call reciprocal_obj%calculate_bsf('bsf.dat')
   end subroutine post_processing_bsf

   !> @brief k-space density-of-states post-processing entry point.
   !> @details The routine builds the same bulk-only stack as
   !>          post_processing_band_structure. It then calls
   !>          reciprocal%calculate_density_of_states with the method, energy
   !>          range, temperature, and Fermi-search settings from the
   !>          &reciprocal namelist, and writes dos_kspace.dat.
   !> @param[in] this Calculation object. fname selects the namelist input.
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
      complex(rp), allocatable :: gij_le(:, :, :, :)
      integer :: i, ie, ne, nblk, norb_h, p0, newunit
      real(rp) :: max_abs_diff, rms_diff, int_rs, int_le, de
      real(rp) :: max_mz_diff, max_blk_diff, max_dyson_diff
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

      ! --- Backend D (Dyson, Sigma=0) cross-check on the real H(k). ---------------
      ! The permanent B2.4 invariant "D with Sigma=0 == E" over EVERY pair (the
      ! full gij array, so the intersite e^{ik.dR} phase is exercised, not just
      ! the on-site block). Both routes use S=I (backend E's zheev is orthonormal,
      ! so backend D inverts z*I - H(k)); the difference is solver-tolerance ripple.
      allocate (gij_le(nblk, nblk, ne, size(green_obj%gij, 4)))
      gij_le = green_obj%gij
      reciprocal_obj%green_backend = 'dyson'
      call reciprocal_obj%fill_green(green_obj, sigma)
      max_dyson_diff = maxval(abs(green_obj%gij - gij_le))
      green_obj%gij = gij_le   ! restore backend-E result for the report below
      deallocate (gij_le)

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
         write (newunit, '(A,ES14.6)') '# B2.4 max|gij^dyson - gij^lehmann|  = ', max_dyson_diff
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
         call g_logger%info('[kspace_green] B2.4 backend-D (Dyson, Sigma=0) == backend-E '// &
                            '(Lehmann) cross-check: max|gij_dyson - gij_lehmann|='// &
                            trim(real2str(max_dyson_diff))//' (solver-tolerance invariant)', &
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

   !---------------------------------------------------------------------------
   !> @brief Finite-q transverse KS susceptibility production route.
   !> @details This is intentionally MPI-over-q: each rank owns complete,
   !> caller-local k/k+q eigenpair batches for its assigned q points.  That
   !> keeps the existing reciprocal k mesh unmodified and makes every emitted
   !> q file independent, deterministic, and restart-friendly.  Frequency
   !> batching remains inside build_chi_ks_from_eigenpairs.
   !>
   !> The transverse kernel is refreshed through self%VXC0SP and passed from
   !> self%xc_response_provider into the Goldstone/Dyson core.  No kernel is
   !> inferred from hamiltonian%hxc/cx1 or another assembled Hamiltonian block.
   !---------------------------------------------------------------------------
   subroutine post_processing_susceptibility(this)
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
      type(self) :: self_obj
      type(reciprocal) :: reciprocal_obj
      type(tddft_config) :: config
      type(tddft_chi0_options) :: chi0_options
      type(green_chi0_options) :: green_options
      type(eigenpair_green_function_provider), target :: green_source
      type(tddft_chi0_result) :: chi0_result, chi0_static
      type(tddft_dyson_options) :: dyson_options
      type(tddft_dyson_result) :: dyson_result, dyson_pair_result, dyson_pair_corrected_result
      type(tddft_direct_xi_result) :: pair_xi_result, pair_xi_static, pair_xi_corrected_result
      type(tddft_goldstone_options) :: goldstone_options
      type(tddft_goldstone_result) :: goldstone_result
      type(tddft_goldstone_diagnostics) :: pair_goldstone
      type(tddft_goldstone_column_correction) :: pair_correction
      type(tddft_mode_options) :: mode_options
      type(tddft_mode_result) :: mode_result
      type(tddft_longitudinal_options) :: longitudinal_options
      type(tddft_longitudinal_static_result) :: longitudinal_static
      type(tddft_longitudinal_result) :: longitudinal_result
      type(tddft_four_component_zero_mode_diagnostics) :: full_zero_mode_diagnostics
      type(response_channel), allocatable :: left_channels(:), right_channels(:)
      real(rp), allocatable :: omega(:), omega_static(:), eigenvalues_k(:, :), eigenvalues_kq(:, :), kq_points(:, :)
      complex(rp), allocatable :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :), kernel(:, :), all_xi(:, :, :, :)
      complex(rp), allocatable :: pair_operators(:, :, :, :), pair_operators_static(:, :, :, :), &
         pair_operators_corrected(:, :, :, :), all_xi_pair(:, :, :, :)
      real(rp), allocatable :: all_trace_loss(:, :), all_trace_loss_pair(:, :), coulomb_site(:, :), magnetization(:, :), site_moments(:, :)
      real(rp), allocatable :: signed_moments(:)
      real(rp), allocatable :: m0(:), static_fields(:), static_moments(:, :)
      real(rp) :: response_eta, t_profile_start, t_profile_stop, kq_eigensolve_cpu_seconds
      real(rp) :: response_electron_count, response_band_energy, electron_count_tolerance
      real(rp) :: bare_gamma_peak, legacy_gamma_peak, pair_gamma_peak, pair_corrected_gamma_peak
      integer, allocatable :: site_orbital_counts(:), static_sources(:)
      integer :: iq, iq_start, iq_end, nq_per_rank, nq, nw, unit, ios, isite, nresponse
      logical :: has_soc, has_external_field, need_dyson, is_longitudinal, is_full_response, is_gamma, has_gamma
      logical :: pair_backend, legacy_backend, corrected_spectral_weight_ok
      character(len=sl) :: filename, chi0_filename, legacy_filename, pair_filename

      config = tddft_config(this%fname)
      if (.not. config%enabled) then
         if (rank == 0) call g_logger%info('TDDFT susceptibility route disabled by &tddft enabled=.false.', __FILE__, __LINE__)
         return
      end if
      call prepare_post_processing_stack(this, .false., .false., .true., .false., control_obj, lattice_obj, &
         charge_obj, mix_obj, energy_obj, hamiltonian_obj, recursion_obj, dos_obj, green_obj, bands_obj)
      if (control_obj%calctype /= 'B') then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: eigenpair TDDFT currently requires calctype=''B''.', &
            __FILE__, __LINE__)
      end if

      ! The response driver may be invoked after a normal SCF object has gone
      ! out of scope.  Refresh the provider from the same VXC0SP route, then
      ! rebuild the normal-state Hamiltonian from that ground-state potential.
      call get_mpi_variables(rank, lattice_obj%nrec)
      self_obj = self(bands_obj, mix_obj)
      call self_obj%refresh_xc_response_kernel()
      if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham()
      call hamiltonian_obj%build_bulkham()
      reciprocal_obj = reciprocal(hamiltonian_obj)
      ! Response uses a complete reciprocal mesh.  A reduced mesh cannot in
      ! general be paired with k+q without response-specific symmetry weights.
      reciprocal_obj%use_symmetry_reduction = .false.
      reciprocal_obj%use_time_reversal = .false.
      call reciprocal_obj%generate_mp_mesh()
      ! TDDFT inherits the canonical reciprocal occupation contract unless a
      ! value was explicitly placed in &tddft.  This closes the old hidden
      ! default (EF=0, T=300 K) which could silently differ from the SCF state.
      config%ground_state_fermi_level = reciprocal_obj%fermi_level
      config%ground_state_electronic_temperature = reciprocal_obj%temperature
      config%ground_state_electron_count = reciprocal_obj%total_electrons
      if (.not. config%fermi_level_overridden) config%fermi_level = config%ground_state_fermi_level
      if (.not. config%electronic_temperature_overridden) then
         config%electronic_temperature = config%ground_state_electronic_temperature
      end if
      if (config%electronic_temperature < 0.0_rp) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: response temperature is unresolved.', __FILE__, __LINE__)
      end if
      reciprocal_obj%fermi_level = config%fermi_level
      reciprocal_obj%temperature = config%electronic_temperature
      is_longitudinal = config%channel == 'longitudinal'
      is_full_response = config%channel == 'full'
      pair_backend = config%xi_backend == 'pair_potential' .or. config%xi_backend == 'compare'
      legacy_backend = config%xi_backend == 'legacy_site_scalar' .or. config%xi_backend == 'compare'
      if (pair_backend) then
         if (is_longitudinal .or. is_full_response .or. config%chi0_backend /= 'eigenpairs') then
            call g_logger%fatal('[calculation.post_processing_susceptibility]: xi_backend=pair_potential/compare is currently '// &
               'restricted to transverse chi0_backend=eigenpairs.', __FILE__, __LINE__)
         end if
         if (trim(reciprocal_obj%reciprocal_mode) /= 'ham_only') then
            call g_logger%fatal('[calculation.post_processing_susceptibility]: pair-potential Xi requires reciprocal_mode=ham_only; '// &
               'generalized-overlap response is explicitly unsupported.', __FILE__, __LINE__)
         end if
         if (.not. (config%output_xi .or. config%output_chi)) then
            call g_logger%fatal('[calculation.post_processing_susceptibility]: pair-potential shadow workflows require '// &
               'output_xi=.true. or output_chi=.true. so raw Xi is auditable.', __FILE__, __LINE__)
         end if
      end if
      if (config%goldstone_mode == 'correct' .and. .not. pair_backend) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: goldstone_mode=correct requires '// &
            'xi_backend=pair_potential or compare; legacy site-scalar correction was removed.', __FILE__, __LINE__)
      end if

      nq = size(config%q_points, 2)
      call get_mpi_range(rank, nq, iq_start, iq_end, nq_per_rank, region_tag='tddft-q')
      nw = config%nomega
      allocate(omega(nw))
      if (nw == 1) then
         omega(1) = config%omega_min
      else
         do iq = 1, nw
            omega(iq) = config%omega_min + real(iq-1, rp)*(config%omega_max-config%omega_min)/real(nw-1, rp)
         end do
      end if

      allocate(site_orbital_counts(lattice_obj%nrec), left_channels(lattice_obj%nrec), right_channels(lattice_obj%nrec))
      site_orbital_counts = norb
      do iq = 1, lattice_obj%nrec
         if (is_longitudinal) then
            ! sigma_z is diagonal in the collinear spin basis: the generic
            ! vertex engine therefore retains same-spin electron-hole pairs.
            left_channels(iq) = response_channel(iq, RESPONSE_MZ)
            right_channels(iq) = response_channel(iq, RESPONSE_MZ)
         else
            left_channels(iq) = response_channel(iq, RESPONSE_PLUS)
            right_channels(iq) = response_channel(iq, RESPONSE_MINUS)
         end if
      end do
      if (is_full_response) then
         nresponse = 4*lattice_obj%nrec
      else
         nresponse = lattice_obj%nrec
      end if
      chi0_options%eta = config%eta
      chi0_options%fermi_level = config%fermi_level
      chi0_options%electronic_temperature = config%electronic_temperature
      chi0_options%band_first = config%band_first
      chi0_options%band_last = config%band_last
      chi0_options%occupation_prune_tolerance = config%occupation_tolerance
      chi0_options%k_mesh_shape = reciprocal_obj%nk_mesh
      green_options%eta = config%eta
      green_options%green_eta = config%green_eta
      green_options%fermi_level = config%fermi_level
      green_options%electronic_temperature = config%electronic_temperature
      green_options%energy_min = config%green_energy_min
      green_options%energy_max = config%green_energy_max
      green_options%energy_points = config%green_energy_points
      green_options%k_mesh_shape = reciprocal_obj%nk_mesh
      has_soc = .false.
      do isite = 1, lattice_obj%ntype
         has_soc = has_soc .or. any(abs(lattice_obj%symbolic_atoms(isite)%potential%xi_p) > tiny(1.0_rp)) .or. &
            any(abs(lattice_obj%symbolic_atoms(isite)%potential%xi_d) > tiny(1.0_rp))
      end do
      has_external_field = control_obj%do_comom .or. control_obj%constraints_enable
      if (is_longitudinal .and. has_soc) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: TDDFT-08 longitudinal response is restricted to collinear no-SOC calculations.', &
            __FILE__, __LINE__)
      end if

      ! k eigenpairs are independent of q and are therefore reused on each q
      ! worker.  k+q eigenpairs remain caller-owned and exact at off-mesh q.
      call reciprocal_obj%calculate_eigenpairs_at_kpoints(reciprocal_obj%k_points, eigenvalues_k, eigenvectors_k)

      ! `*_out.nml` restart files retain the potential and its direction but
      ! not the scalar site moment `mtot`.  The XC radial projection recorded
      ! by refresh_xc_response_kernel must be normalized by the same occupied
      ! P_site sigma population used by the response vertices.  Reconstruct it
      ! from the complete, unreduced response mesh rather than relying on that
      ! non-serialized legacy cache.
      reciprocal_obj%eigenvalues = eigenvalues_k
      reciprocal_obj%eigenvectors = eigenvectors_k
      call reciprocal_obj%evaluate_eigenvalue_occupations(config%fermi_level, response_electron_count, response_band_energy)
      config%response_electron_count = response_electron_count
      electron_count_tolerance = 1.0e-8_rp*max(1.0_rp, config%ground_state_electron_count)
      if (abs(response_electron_count-config%ground_state_electron_count) > electron_count_tolerance) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: response electron count does not match the '// &
            'converged reciprocal ground state within 1e-8 relative tolerance.', __FILE__, __LINE__)
      end if
      allocate(site_moments(3, lattice_obj%nrec), signed_moments(lattice_obj%nrec))
      call self_obj%compute_kspace_spin_moments_spinor(reciprocal_obj, site_moments)
      do isite = 1, lattice_obj%nrec
         ! A transverse Goldstone vector is signed.  Replacing this by its
         ! magnitude breaks reversed and multi-sublattice reference states.
         call self_obj%xc_response_provider%set_site_spin_population(isite, site_moments(3, isite))
         call self_obj%xc_response_provider%set_site_signed_spin_population(isite, site_moments(3, isite))
         signed_moments(isite) = site_moments(3, isite)
         if (sqrt(sum(site_moments(:, isite)**2)) > tiny(1.0_rp)) then
            call self_obj%xc_response_provider%set_site_magnetization_direction(isite, site_moments(:, isite))
         end if
      end do
      deallocate(site_moments)

      ! The Gamma/static KS batch is shared by both channels.  Transverse uses
      ! it for Goldstone diagnostics; longitudinal instead calibrates its own
      ! kernel from independently converged symmetric +/- Bz calculations.
      allocate(omega_static(1))
      omega_static = 0.0_rp
      if (.not. is_longitudinal .and. .not. is_full_response) then
         if (config%chi0_backend /= 'eigenpairs') then
            call g_logger%fatal('[calculation.post_processing_susceptibility]: real static Ward diagnostics require '// &
               'chi0_backend=eigenpairs; the eigenpair-resolvent backend has no static-limit solver.', __FILE__, __LINE__)
         end if
         call build_static_chi_ks_from_eigenpairs(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, &
            site_orbital_counts, left_channels, right_channels, chi0_options, chi0_static)
      else if (config%chi0_backend == 'green') then
         call green_source%initialize(eigenvalues_k, eigenvectors_k, eigenvalues_k, eigenvectors_k)
         if (is_full_response) then
            call build_four_component_chi_ks_from_green_functions(green_source, reciprocal_obj%k_weights, &
               site_orbital_counts, omega_static, green_options, chi0_static)
         else
            call build_chi_ks_from_green_functions(green_source, reciprocal_obj%k_weights, site_orbital_counts, &
               left_channels, right_channels, omega_static, green_options, chi0_static)
         end if
      else if (is_full_response) then
         call build_four_component_chi_ks(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_k, &
            eigenvectors_k, site_orbital_counts, omega_static, chi0_options, chi0_static)
      else
         call build_chi_ks_from_eigenpairs(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_k, &
            eigenvectors_k, site_orbital_counts, left_channels, right_channels, omega_static, chi0_options, chi0_static)
      end if
      if (pair_backend) then
         call build_pair_potential_operators(reciprocal_obj, reciprocal_obj%k_points, signed_moments, &
            [0.0_rp, 0.0_rp, 0.0_rp], pair_operators_static)
         call build_static_direct_xi_from_k_dependent_eigenpairs(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, &
            site_orbital_counts, left_channels, pair_operators_static, chi0_options, pair_xi_static)
         call evaluate_raw_xi_diagnostics(pair_xi_static%xi(:, :, 1), cmplx(signed_moments, 0.0_rp, rp), pair_goldstone)
         if (config%goldstone_mode == 'correct') then
            if (has_soc .or. has_external_field) then
               pair_correction%requested = .true.
               pair_correction%rejected = .true.
               pair_correction%decision = 'Goldstone correction is unavailable with SOC or an external symmetry-breaking field'
            else
               call build_goldstone_column_correction(pair_xi_static%xi(:, :, 1), cmplx(signed_moments, 0.0_rp, rp), &
                  pair_correction)
            end if
         end if
         deallocate(pair_operators_static)
      end if
      allocate(kernel(nresponse, nresponse))
      kernel = cmplx(0.0_rp, 0.0_rp, rp)
      if (is_longitudinal) then
         longitudinal_options%pair_tolerance = config%longitudinal_pair_tolerance
         longitudinal_options%linearity_tolerance = config%longitudinal_linearity_tolerance
         longitudinal_options%static_agreement_tolerance = config%longitudinal_static_agreement_tolerance
         longitudinal_options%fit_omega_min = config%longitudinal_fit_omega_min
         longitudinal_options%fit_omega_max = config%longitudinal_fit_omega_max
         call read_longitudinal_static_fields(trim(config%longitudinal_static_file), lattice_obj%nrec, static_sources, &
            static_fields, static_moments)
         call build_longitudinal_static_response(static_sources, static_fields, static_moments, longitudinal_options, &
            longitudinal_static)
         allocate(m0(lattice_obj%nrec))
         do isite = 1, lattice_obj%nrec
            m0(isite) = self_obj%xc_response_provider%site(isite)%spin_population
         end do
         call build_longitudinal_kernel(chi0_static%chi(:, :, 1), cmplx(longitudinal_static%chi, 0.0_rp, rp), kernel)
      else if (is_full_response) then
         ! The common XC provider is the sole source of ALSDA derivatives;
         ! the charge response additionally uses the projected Madelung
         ! charge kernel.  Refuse a dimension mismatch instead of treating a
         ! real-space matrix as a different response projection.
         if (.not. allocated(charge_obj%amad) .or. size(charge_obj%amad, 1) /= lattice_obj%nrec .or. &
            size(charge_obj%amad, 2) /= lattice_obj%nrec) then
            call g_logger%fatal('[calculation.post_processing_susceptibility]: full TDDFT requires an nrec-by-nrec projected Coulomb kernel.', &
               __FILE__, __LINE__)
         end if
         allocate(coulomb_site(lattice_obj%nrec, lattice_obj%nrec), magnetization(3, lattice_obj%nrec))
         coulomb_site = charge_obj%amad
         do isite = 1, lattice_obj%nrec
            magnetization(:, isite) = self_obj%xc_response_provider%site(isite)%spin_population* &
               self_obj%xc_response_provider%site(isite)%magnetization_direction
         end do
         call build_four_component_kernel(self_obj%xc_response_provider, coulomb_site, kernel)
         call evaluate_four_component_zero_modes(chi0_static%chi(:, :, 1), kernel, magnetization, has_soc, &
            has_external_field, full_zero_mode_diagnostics)
         if (rank == 0 .and. full_zero_mode_diagnostics%applicable) then
            call g_logger%info('TDDFT full response rigid-rotation zero-mode residual = '// &
               real2str(maxval(full_zero_mode_diagnostics%residual)), __FILE__, __LINE__)
         end if
      else
         goldstone_options%goldstone_mode = config%goldstone_mode
         goldstone_options%has_soc = has_soc
         goldstone_options%has_external_field = has_external_field
         call evaluate_goldstone(chi0_static%chi(:, :, 1), self_obj%xc_response_provider, goldstone_options, goldstone_result)
         ! The finite-eta static sum-rule inverse remains a diagnostic only.
         ! It is generally complex; using it as a frequency-independent
         ! adiabatic kernel forces an artificial Gamma singularity and can
         ! produce negative spectral weight.  Production Dyson spectra use
         ! the real, XC-derived kernel below.
         if (legacy_backend) then
            do isite = 1, lattice_obj%nrec
               kernel(isite, isite) = goldstone_result%k_perp(isite)
            end do
         end if
         if (rank == 0) then
            write(filename, '(a,"_goldstone.dat")') trim(config%output_prefix)
            call write_goldstone_diagnostics_text(trim(filename), goldstone_result)
            if (pair_backend) call append_pair_goldstone_diagnostics(trim(filename), goldstone_result%raw, pair_goldstone)
            if (config%goldstone_mode == 'correct') call append_goldstone_column_correction_text(trim(filename), pair_correction)
            call append_tddft_metadata(trim(filename), config, 0, reciprocal_obj%nk_mesh, [0.0_rp, 0.0_rp, 0.0_rp], &
               rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'goldstone_compare')
         end if
         if (config%goldstone_mode == 'correct' .and. .not. pair_correction%applied) then
            call g_logger%fatal('[calculation.post_processing_susceptibility]: requested Goldstone correction rejected: '// &
               trim(pair_correction%decision), __FILE__, __LINE__)
         end if
      end if

      has_gamma = any(maxval(abs(config%q_points), dim=1) <= 1.0e-12_rp)
      if (is_longitudinal .and. .not. has_gamma) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: longitudinal response requires q=Gamma for static acceptance.', &
            __FILE__, __LINE__)
      end if
      ! Selecting a pair-potential backend is itself a request to construct
      ! the raw Xi shadow data at every requested (q,omega), even if the user
      ! has disabled the optional text products.
      need_dyson = config%output_xi .or. config%output_chi .or. config%output_modes .or. is_longitudinal .or. pair_backend
      dyson_options%diagonalize_xi = config%output_modes
      dyson_options%diagonalize_loss = config%output_modes
      if (config%output_modes) then
         allocate(all_xi(nresponse, nresponse, nw, nq), all_trace_loss(nw, nq))
         all_xi = cmplx(0.0_rp, 0.0_rp, rp)
         all_trace_loss = 0.0_rp
         if (pair_backend .and. legacy_backend) then
            allocate(all_xi_pair(nresponse, nresponse, nw, nq), all_trace_loss_pair(nw, nq))
            all_xi_pair = cmplx(0.0_rp, 0.0_rp, rp)
            all_trace_loss_pair = 0.0_rp
         end if
      end if
#ifdef USE_MPI
      ! Root writes the static Goldstone record above; q owners may append the
      ! independently observed dynamic Gamma peaks below.
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
      do iq = iq_start, iq_end
         is_gamma = maxval(abs(config%q_points(:, iq))) <= 1.0e-12_rp
         bare_gamma_peak = -1.0_rp; legacy_gamma_peak = -1.0_rp; pair_gamma_peak = -1.0_rp
         allocate(kq_points(3, reciprocal_obj%nk_total))
         kq_points = reciprocal_obj%k_points + spread(config%q_points(:, iq), dim=2, ncopies=reciprocal_obj%nk_total)
         call cpu_time(t_profile_start)
         call reciprocal_obj%calculate_eigenpairs_at_kpoints(kq_points, eigenvalues_kq, eigenvectors_kq)
         call cpu_time(t_profile_stop)
         kq_eigensolve_cpu_seconds = t_profile_stop-t_profile_start
         if (config%chi0_backend == 'green') then
            call green_source%initialize(eigenvalues_k, eigenvectors_k, eigenvalues_kq, eigenvectors_kq)
            if (is_full_response) then
               call build_four_component_chi_ks_from_green_functions(green_source, reciprocal_obj%k_weights, &
                  site_orbital_counts, omega, green_options, chi0_result)
            else
               call build_chi_ks_from_green_functions(green_source, reciprocal_obj%k_weights, site_orbital_counts, &
                  left_channels, right_channels, omega, green_options, chi0_result)
            end if
         else if (is_full_response) then
            call build_four_component_chi_ks(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_kq, &
               eigenvectors_kq, site_orbital_counts, omega, chi0_options, chi0_result)
         else
            call build_chi_ks_from_eigenpairs(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_kq, &
               eigenvectors_kq, site_orbital_counts, left_channels, right_channels, omega, chi0_options, chi0_result)
         end if
         chi0_result%metadata%arbitrary_kq_cpu_seconds = kq_eigensolve_cpu_seconds
         response_eta = chi0_result%metadata%eta
         if (is_gamma) bare_gamma_peak = observed_loss_peak(omega, chi0_result%trace_spectrum)
         if (config%output_chi0 .or. config%output_stoner) then
            write(filename, '(a,"_q",i6.6,"_chi0.dat")') trim(config%output_prefix), iq
            call write_chi_ks_text(trim(filename), omega, chi0_result)
            call append_tddft_metadata(trim(filename), config, iq, reciprocal_obj%nk_mesh, config%q_points(:, iq), rank, &
               has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'shared_chi_ks')
         end if
         if (need_dyson) then
            if (legacy_backend) then
               call enhance_tddft_susceptibility(chi0_result%chi, kernel, response_eta, dyson_options, dyson_result)
               if (is_gamma) legacy_gamma_peak = observed_loss_peak(omega, dyson_result%trace_spectral_weight)
               if (config%output_xi .or. config%output_chi) then
                  if (pair_backend) then
                     write(filename, '(a,"_q",i6.6,"_legacy_dyson.dat")') trim(config%output_prefix), iq
                  else
                     write(filename, '(a,"_q",i6.6,"_dyson.dat")') trim(config%output_prefix), iq
                  end if
                  call write_tddft_dyson_text(trim(filename), omega, dyson_result)
                  call append_tddft_metadata(trim(filename), config, iq, reciprocal_obj%nk_mesh, config%q_points(:, iq), &
                     rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'legacy_site_scalar_raw')
               end if
            end if
            if (pair_backend) then
               call build_pair_potential_operators(reciprocal_obj, reciprocal_obj%k_points, signed_moments, &
                  config%q_points(:, iq), pair_operators)
               call build_direct_xi_from_k_dependent_eigenpairs(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, &
                  eigenvalues_kq, eigenvectors_kq, site_orbital_counts, left_channels, pair_operators, omega, chi0_options, &
                  pair_xi_result)
               call enhance_tddft_susceptibility_from_xi(chi0_result%chi, pair_xi_result%xi, response_eta, dyson_options, &
                  dyson_pair_result)
               if (is_gamma) pair_gamma_peak = observed_loss_peak(omega, dyson_pair_result%trace_spectral_weight)
               if (config%output_xi .or. config%output_chi) then
                  write(filename, '(a,"_q",i6.6,"_pair_dyson.dat")') trim(config%output_prefix), iq
                  call write_tddft_dyson_text(trim(filename), omega, dyson_pair_result)
                  call append_tddft_metadata(trim(filename), config, iq, reciprocal_obj%nk_mesh, config%q_points(:, iq), &
                     rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'pair_potential_raw')
               end if
               if (config%goldstone_mode == 'correct') then
                  allocate(pair_operators_corrected, source=pair_operators)
                  call rescale_pair_potential_columns(pair_operators_corrected, pair_correction%scales)
                  call build_direct_xi_from_k_dependent_eigenpairs(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, &
                     eigenvalues_kq, eigenvectors_kq, site_orbital_counts, left_channels, pair_operators_corrected, omega, &
                     chi0_options, pair_xi_corrected_result)
                  call enhance_tddft_susceptibility_from_xi(chi0_result%chi, pair_xi_corrected_result%xi, response_eta, &
                     dyson_options, dyson_pair_corrected_result)
                  if (is_gamma) then
                     pair_corrected_gamma_peak = observed_loss_peak(omega, dyson_pair_corrected_result%trace_spectral_weight)
                  end if
                  corrected_spectral_weight_ok = spectral_weights_are_nonnegative(reshape( &
                     dyson_pair_corrected_result%site_spectral_weight, [size(dyson_pair_corrected_result%site_spectral_weight)]))
                  if (config%output_xi .or. config%output_chi) then
                     write(filename, '(a,"_q",i6.6,"_pair_corrected_dyson.dat")') trim(config%output_prefix), iq
                     call write_tddft_dyson_text(trim(filename), omega, dyson_pair_corrected_result)
                     call append_tddft_metadata(trim(filename), config, iq, reciprocal_obj%nk_mesh, config%q_points(:, iq), &
                        rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'pair_potential_corrected')
                     call append_goldstone_column_correction_text(trim(filename), pair_correction)
                     call append_corrected_spectral_weight_diagnostic(trim(filename), corrected_spectral_weight_ok, &
                        minval(dyson_pair_corrected_result%site_spectral_weight))
                  end if
                  deallocate(pair_operators_corrected)
                  if (.not. corrected_spectral_weight_ok) then
                     call g_logger%fatal('[calculation.post_processing_susceptibility]: corrected pair-potential spectrum has '// &
                        'negative spectral weight beyond tolerance.', __FILE__, __LINE__)
                  end if
               end if
               deallocate(pair_operators)
            end if
            if (is_longitudinal .and. is_gamma) then
               call calibrate_longitudinal_response(m0, chi0_static%chi(:, :, 1), cmplx(longitudinal_static%chi, 0.0_rp, rp), &
                  omega, dyson_result%chi, longitudinal_options, longitudinal_result)
               write(filename, '(a,"_q",i6.6,"_longitudinal.dat")') trim(config%output_prefix), iq
               call write_longitudinal_report(trim(filename), omega, response_eta, longitudinal_static, longitudinal_result)
               call append_tddft_metadata(trim(filename), config, iq, reciprocal_obj%nk_mesh, config%q_points(:, iq), rank, &
                  has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'longitudinal_static')
            end if
            if (config%output_modes) then
               if (legacy_backend) then
                  all_xi(:, :, :, iq) = dyson_result%xi
                  all_trace_loss(:, iq) = dyson_result%trace_spectral_weight
               else
                  all_xi(:, :, :, iq) = dyson_pair_result%xi
                  all_trace_loss(:, iq) = dyson_pair_result%trace_spectral_weight
               end if
               if (pair_backend .and. legacy_backend) then
                  all_xi_pair(:, :, :, iq) = dyson_pair_result%xi
                  all_trace_loss_pair(:, iq) = dyson_pair_result%trace_spectral_weight
               end if
            end if
         end if
         if (is_gamma .and. .not. is_longitudinal .and. .not. is_full_response) then
            write(filename, '(a,"_goldstone.dat")') trim(config%output_prefix)
            call append_dynamic_gamma_peaks(trim(filename), bare_gamma_peak, legacy_gamma_peak, pair_gamma_peak, &
               pair_corrected_gamma_peak)
         end if
         deallocate(kq_points, eigenvalues_kq, eigenvectors_kq)
      end do

      if (config%output_modes) then
#ifdef USE_MPI
         call MPI_ALLREDUCE(MPI_IN_PLACE, all_xi, size(all_xi), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, all_trace_loss, size(all_trace_loss), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         if (pair_backend .and. legacy_backend) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, all_xi_pair, size(all_xi_pair), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, all_trace_loss_pair, size(all_trace_loss_pair), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         end if
#endif
         if (rank == 0) then
            call analyze_tddft_modes(omega, all_xi, all_trace_loss, response_eta, mode_options, mode_result)
            if (pair_backend .and. legacy_backend) then
               write(filename, '(a,"_legacy_modes.dat")') trim(config%output_prefix)
            else if (pair_backend) then
               write(filename, '(a,"_pair_modes.dat")') trim(config%output_prefix)
            else
               write(filename, '(a,"_modes.dat")') trim(config%output_prefix)
            end if
            call write_tddft_modes_text(trim(filename), omega, response_eta, mode_result)
            if (legacy_backend) then
               call append_tddft_metadata(trim(filename), config, 0, reciprocal_obj%nk_mesh, [0.0_rp, 0.0_rp, 0.0_rp], &
                  rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'legacy_site_scalar_raw')
            else
               call append_tddft_metadata(trim(filename), config, 0, reciprocal_obj%nk_mesh, [0.0_rp, 0.0_rp, 0.0_rp], &
                  rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'pair_potential_raw')
            end if
            if (pair_backend .and. legacy_backend) then
               call analyze_tddft_modes(omega, all_xi_pair, all_trace_loss_pair, response_eta, mode_options, mode_result)
               write(filename, '(a,"_pair_modes.dat")') trim(config%output_prefix)
               call write_tddft_modes_text(trim(filename), omega, response_eta, mode_result)
               call append_tddft_metadata(trim(filename), config, 0, reciprocal_obj%nk_mesh, [0.0_rp, 0.0_rp, 0.0_rp], &
                  rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'pair_potential_raw')
            end if
         end if
      end if

      ! A deterministic manifest is written once; individual q files are
      ! never concurrently opened by more than one rank.
      if (rank == 0) then
         write(filename, '(a,"_manifest.dat")') trim(config%output_prefix)
         open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
         if (ios /= 0) then
            call g_logger%fatal('[calculation.post_processing_susceptibility]: cannot write manifest', __FILE__, __LINE__)
         end if
         if (is_longitudinal) then
            write(unit, '(a)') '# TDDFT longitudinal chi_KS manifest; static calibration is reported at Gamma'
         else
            write(unit, '(a)') '# TDDFT transverse chi_KS manifest; one self-describing output file per q point'
         end if
         if (pair_backend .and. legacy_backend .and. config%goldstone_mode == 'correct') then
            write(unit, '(a)') '# q_index q1 q2 q3 chi0_file legacy_raw_dyson_file pair_raw_dyson_file pair_corrected_dyson_file'
         else if (pair_backend .and. legacy_backend) then
            write(unit, '(a)') '# q_index q1 q2 q3 chi0_file legacy_raw_dyson_file pair_raw_dyson_file'
         else if (pair_backend .and. config%goldstone_mode == 'correct') then
            write(unit, '(a)') '# q_index q1 q2 q3 chi0_file pair_raw_dyson_file pair_corrected_dyson_file'
         else if (pair_backend) then
            write(unit, '(a)') '# q_index q1 q2 q3 chi0_file pair_raw_dyson_file'
         else
            write(unit, '(a)') '# q_index q1 q2 q3 chi0_file legacy_raw_dyson_file'
         end if
         do iq = 1, nq
            write(chi0_filename, '(a,"_q",i6.6,"_chi0.dat")') trim(config%output_prefix), iq
            if (pair_backend .and. legacy_backend .and. config%goldstone_mode == 'correct') then
               write(legacy_filename, '(a,"_q",i6.6,"_legacy_dyson.dat")') trim(config%output_prefix), iq
               write(pair_filename, '(a,"_q",i6.6,"_pair_dyson.dat")') trim(config%output_prefix), iq
               write(filename, '(a,"_q",i6.6,"_pair_corrected_dyson.dat")') trim(config%output_prefix), iq
               write(unit, '(i0,3(1x,es24.16),4(1x,a))') iq, config%q_points(:, iq), trim(chi0_filename), &
                  trim(legacy_filename), trim(pair_filename), trim(filename)
            else if (pair_backend .and. legacy_backend) then
               write(legacy_filename, '(a,"_q",i6.6,"_legacy_dyson.dat")') trim(config%output_prefix), iq
               write(pair_filename, '(a,"_q",i6.6,"_pair_dyson.dat")') trim(config%output_prefix), iq
               write(unit, '(i0,3(1x,es24.16),3(1x,a))') iq, config%q_points(:, iq), trim(chi0_filename), &
                  trim(legacy_filename), trim(pair_filename)
            else if (pair_backend .and. config%goldstone_mode == 'correct') then
               write(pair_filename, '(a,"_q",i6.6,"_pair_dyson.dat")') trim(config%output_prefix), iq
               write(filename, '(a,"_q",i6.6,"_pair_corrected_dyson.dat")') trim(config%output_prefix), iq
               write(unit, '(i0,3(1x,es24.16),3(1x,a))') iq, config%q_points(:, iq), trim(chi0_filename), &
                  trim(pair_filename), trim(filename)
            else if (pair_backend) then
               write(pair_filename, '(a,"_q",i6.6,"_pair_dyson.dat")') trim(config%output_prefix), iq
               write(unit, '(i0,3(1x,es24.16),2(1x,a))') iq, config%q_points(:, iq), trim(chi0_filename), trim(pair_filename)
            else
               write(legacy_filename, '(a,"_q",i6.6,"_dyson.dat")') trim(config%output_prefix), iq
               write(unit, '(i0,3(1x,es24.16),2(1x,a))') iq, config%q_points(:, iq), trim(chi0_filename), trim(legacy_filename)
            end if
         end do
         close(unit)
      end if
   end subroutine post_processing_susceptibility

   !> Report observed dynamic Gamma maxima separately from the real static
   !> Ward operator.  These are grid-resolved loss maxima, not a correction or
   !> a fitted frequency shift, and retain both legacy and pair raw routes.
   real(rp) function observed_loss_peak(omega, loss) result(peak)
      real(rp), intent(in) :: omega(:), loss(:)
      integer :: index(1)

      if (size(omega) /= size(loss) .or. size(omega) < 1) then
         error stop 'observed_loss_peak: incompatible omega/loss arrays'
      end if
      index = maxloc(loss)
      peak = omega(index(1))
   end function observed_loss_peak

   subroutine append_dynamic_gamma_peaks(filename, bare_peak, legacy_peak, pair_peak, pair_corrected_peak)
      character(len=*), intent(in) :: filename
      real(rp), intent(in) :: bare_peak, legacy_peak, pair_peak, pair_corrected_peak
      integer :: unit, ios

      open(newunit=unit, file=filename, status='old', position='append', action='write', iostat=ios)
      if (ios /= 0) call g_logger%fatal('[calculation.append_dynamic_gamma_peaks]: cannot append Gamma peaks', __FILE__, __LINE__)
      write(unit, '(a,es24.16)') '# dynamic_bare_gamma_loss_peak_Ry = ', bare_peak
      if (legacy_peak >= 0.0_rp) write(unit, '(a,es24.16)') '# dynamic_legacy_raw_gamma_loss_peak_Ry = ', legacy_peak
      if (pair_peak >= 0.0_rp) write(unit, '(a,es24.16)') '# dynamic_pair_raw_gamma_loss_peak_Ry = ', pair_peak
      if (pair_corrected_peak >= 0.0_rp) then
         write(unit, '(a,es24.16)') '# dynamic_pair_corrected_gamma_loss_peak_Ry = ', pair_corrected_peak
      end if
      write(unit, '(a)') '# dynamic Gamma peaks are observed loss-grid maxima; raw and corrected records remain distinct'
      close(unit)
   end subroutine append_dynamic_gamma_peaks

   subroutine append_corrected_spectral_weight_diagnostic(filename, acceptable, minimum_weight)
      character(len=*), intent(in) :: filename
      logical, intent(in) :: acceptable
      real(rp), intent(in) :: minimum_weight
      integer :: unit, ios

      open(newunit=unit, file=filename, status='old', position='append', action='write', iostat=ios)
      if (ios /= 0) call g_logger%fatal('[calculation.append_corrected_spectral_weight_diagnostic]: cannot append diagnostic', &
         __FILE__, __LINE__)
      write(unit, '(a,l1)') '# corrected_pair_spectral_weight_nonnegative = ', acceptable
      write(unit, '(a,es24.16)') '# corrected_pair_minimum_site_spectral_weight = ', minimum_weight
      close(unit)
   end subroutine append_corrected_spectral_weight_diagnostic

   !> Build the Q^- operators in the same ham_only coefficient representation
   !> as the two endpoint eigensystems.  Q is intentionally retained per k:
   !> replacing it by a site scalar or an average here would reintroduce the
   !> projection/multiplication ordering defect repaired by WR-03.
   subroutine build_pair_potential_operators(reciprocal_obj, k_points, signed_moments, q_point, operators)
      type(reciprocal), intent(inout) :: reciprocal_obj
      real(rp), intent(in) :: k_points(:, :), signed_moments(:), q_point(3)
      complex(rp), allocatable, intent(out) :: operators(:, :, :, :)
      complex(rp), allocatable :: qminus(:, :), qplus(:, :)
      integer :: nmat, ik, isite
      logical :: supported
      character(len=160) :: reason

      if (size(k_points, 1) /= 3 .or. size(signed_moments) /= reciprocal_obj%lattice%nrec) then
         call g_logger%fatal('[calculation.build_pair_potential_operators]: incompatible k-point or signed-moment shape.', &
            __FILE__, __LINE__)
      end if
      nmat = 2*norb*reciprocal_obj%lattice%nrec
      allocate(operators(nmat, nmat, reciprocal_obj%lattice%nrec, size(k_points, 2)), qminus(nmat, nmat), qplus(nmat, nmat))
      operators = cmplx(0.0_rp, 0.0_rp, rp)
      do ik = 1, size(k_points, 2)
         do isite = 1, reciprocal_obj%lattice%nrec
            call reciprocal_obj%build_lmto_pair_potential_at_kpoint(isite, k_points(:, ik), signed_moments(isite), &
               qminus, qplus, supported, reason, q_point)
            if (.not. supported) then
               call g_logger%fatal('[calculation.build_pair_potential_operators]: pair-potential construction rejected: '// &
                  trim(reason), __FILE__, __LINE__)
            end if
            operators(:, :, isite, ik) = qminus
         end do
      end do
      deallocate(qminus, qplus)
   end subroutine build_pair_potential_operators

   subroutine append_pair_goldstone_diagnostics(filename, legacy_diagnostics, diagnostics)
      character(len=*), intent(in) :: filename
      type(tddft_goldstone_diagnostics), intent(in) :: legacy_diagnostics
      type(tddft_goldstone_diagnostics), intent(in) :: diagnostics
      integer :: unit, ios

      open(newunit=unit, file=filename, status='old', position='append', action='write', iostat=ios)
      if (ios /= 0) call g_logger%fatal('[calculation.append_pair_goldstone_diagnostics]: cannot append pair diagnostics', &
         __FILE__, __LINE__)
      write(unit, '(a)') '# pair_potential_raw_goldstone_begin'
      write(unit, '(a,es24.16)') '# legacy_site_scalar_raw_residual = ', legacy_diagnostics%residual
      write(unit, '(a,es24.16)') '# legacy_site_scalar_raw_magnetization_overlap = ', legacy_diagnostics%magnetization_overlap
      write(unit, '(a,l1)') '# pair_potential_raw_available = ', diagnostics%available
      write(unit, '(a,es24.16)') '# pair_potential_raw_residual = ', diagnostics%residual
      write(unit, '(a,2(1x,es24.16))') '# pair_potential_raw_closest_eigenvalue = ', real(diagnostics%closest_eigenvalue, rp), &
         aimag(diagnostics%closest_eigenvalue)
      write(unit, '(a,es24.16)') '# pair_potential_raw_magnetization_overlap = ', diagnostics%magnetization_overlap
      write(unit, '(a,es24.16)') '# pair_potential_raw_left_magnetization_overlap = ', diagnostics%left_magnetization_overlap
      write(unit, '(a,es24.16)') '# pair_potential_raw_biorthogonal_magnetization_overlap = ', &
         diagnostics%biorthogonal_magnetization_overlap
      write(unit, '(a,es24.16)') '# pair_potential_raw_imaginary_norm = ', diagnostics%imaginary_norm
      write(unit, '(a)') '# pair_potential_static_solver = real q=0 omega=0 Fermi divided difference; dynamic eta excluded'
      write(unit, '(a)') '# pair_potential_provenance = analytic transverse rotation of ordinary LMTO ham_only operator'
      write(unit, '(a)') '# pair_potential_representation = k-resolved reciprocal ham_only coefficient basis'
      write(unit, '(a)') '# signed_moment_source = reconstructed occupied P_site sigma_z population'
      write(unit, '(a)') '# pair_potential_raw_goldstone_end'
      close(unit)
   end subroutine append_pair_goldstone_diagnostics

   subroutine append_tddft_metadata(filename, config, iq, k_mesh, q_point, mpi_rank, has_soc, has_external_field, &
      reciprocal_mode, xi_backend_label)
      character(len=*), intent(in) :: filename
      type(tddft_config), intent(in) :: config
      integer, intent(in) :: iq, k_mesh(3), mpi_rank
      real(rp), intent(in) :: q_point(3)
      logical, intent(in) :: has_soc, has_external_field
      character(len=*), intent(in) :: reciprocal_mode
      character(len=*), intent(in) :: xi_backend_label
      integer :: unit, ios

      open(newunit=unit, file=filename, status='old', position='append', action='write', iostat=ios)
      if (ios /= 0) call g_logger%fatal('[calculation.append_tddft_metadata]: cannot append response metadata', __FILE__, __LINE__)
      write(unit, '(a)') '# production_metadata_begin'
      write(unit, '(a,a)') '# channel = ', trim(config%channel)
      write(unit, '(a,a)') '# chi0_backend = ', trim(config%chi0_backend)
      write(unit, '(a,a)') '# xi_backend_requested = ', trim(config%xi_backend)
      write(unit, '(a,a)') '# xi_backend_output = ', trim(xi_backend_label)
      if (index(xi_backend_label, 'pair_potential') > 0) then
         write(unit, '(a)') '# pair_potential_provenance = analytic transverse rotation of ordinary LMTO ham_only operator'
         write(unit, '(a)') '# pair_potential_representation = k-resolved reciprocal ham_only coefficient basis'
         write(unit, '(a)') '# signed_moment_source = reconstructed occupied P_site sigma_z population'
         write(unit, '(a,a)') '# reciprocal_mode = ', trim(reciprocal_mode)
      else
         write(unit, '(a)') '# pair_potential_provenance = not used by this output'
         write(unit, '(a)') '# signed_moment_source = reconstructed occupied P_site sigma_z population'
         write(unit, '(a,a)') '# reciprocal_mode = ', trim(reciprocal_mode)
      end if
      write(unit, '(a,a)') '# response_projection = ', trim(config%response_projection)
      write(unit, '(a,a)') '# q_mode = ', trim(config%q_mode)
      write(unit, '(a,a)') '# q_coordinates = ', trim(config%q_coordinates)
      write(unit, '(a,i0)') '# q_index = ', iq
      write(unit, '(a,3(1x,es24.16))') '# q_direct =', q_point
      write(unit, '(a,3(1x,i0))') '# k_mesh =', k_mesh
      write(unit, '(a,es24.16)') '# omega_min_Ry = ', config%omega_min
      write(unit, '(a,es24.16)') '# omega_max_Ry = ', config%omega_max
      write(unit, '(a,i0)') '# nomega = ', config%nomega
      write(unit, '(a,es24.16)') '# eta_Ry = ', config%eta
      write(unit, '(a,es24.16)') '# ground_state_fermi_level_Ry = ', config%ground_state_fermi_level
      write(unit, '(a,es24.16)') '# electronic_temperature_K = ', config%electronic_temperature
      write(unit, '(a,es24.16)') '# ground_state_electronic_temperature_K = ', config%ground_state_electronic_temperature
      write(unit, '(a,es24.16)') '# fermi_level_Ry = ', config%fermi_level
      write(unit, '(a,l1)') '# response_fermi_level_overridden = ', config%fermi_level_overridden
      write(unit, '(a,l1)') '# response_electronic_temperature_overridden = ', config%electronic_temperature_overridden
      write(unit, '(a,2(1x,es24.16))') '# ground_state_response_electron_count = ', &
         config%ground_state_electron_count, config%response_electron_count
      write(unit, '(a)') '# static_ward_solver = real q=0 omega=0 Fermi divided difference; dynamic eta excluded'
      write(unit, '(a,2(1x,i0))') '# band_window_first_last = ', config%band_first, config%band_last
      write(unit, '(a,es24.16)') '# occupation_prune_tolerance = ', config%occupation_tolerance
      write(unit, '(a,es24.16)') '# green_eta_Ry = ', config%green_eta
      write(unit, '(a,2(1x,es24.16))') '# green_energy_window_Ry = ', config%green_energy_min, config%green_energy_max
      write(unit, '(a,i0)') '# green_energy_points = ', config%green_energy_points
      write(unit, '(a,a)') '# goldstone_mode = ', trim(config%goldstone_mode)
      write(unit, '(a,l1)') '# goldstone_mode_migrated_from_sum_rule = ', config%goldstone_mode_migrated_from_sum_rule
      write(unit, '(a,l1)') '# goldstone_correction_requested = ', config%goldstone_mode == 'correct'
      write(unit, '(a,l1)') '# goldstone_correction_applied = ', index(xi_backend_label, 'corrected') > 0
      write(unit, '(a,l1)') '# spin_orbit_present = ', has_soc
      write(unit, '(a,l1)') '# external_symmetry_breaking_field_present = ', has_external_field
      write(unit, '(a,5(1x,l1))') '# output_chi0 output_xi output_chi output_modes output_stoner =', &
         config%output_chi0, config%output_xi, config%output_chi, config%output_modes, config%output_stoner
      write(unit, '(a,i0)') '# mpi_q_owner_rank = ', mpi_rank
      write(unit, '(a)') '# production_metadata_end'
      close(unit)
   end subroutine append_tddft_metadata

   !> @brief Single acoustic-branch frozen-magnon sweep (fm_obj%branch_mode
   !>        /= 'auto' path of post_processing_frozen_magnon).
   !> @details For mode='mft' the routine converges the reference potential
   !>          once, at q_ss_list(:,1). It then reuses that fixed potential for
   !>          one force-theorem band-energy probe (frozen_magnon_probe_energy)
   !>          at every other q. A k-space run picks its q-mesh policy
   !>          (little_group, little_group_common, or full_bz) from
   !>          reciprocal%q_symmetry_policy. For mode='scf' the routine fully
   !>          re-converges the potential at every q instead. Either way it
   !>          builds the Halilov-style dispersion
   !>          omega(q) = 4*(E(q)-E(q_ref)) / (M_tot*sin^2(theta_ss)), with
   !>          M_tot summed over sublattices at the reference point, and writes
   !>          fm_obj%output_file.
   !> @param[in] fm_obj Frozen-magnon namelist settings (q list, mode, output
   !>            file).
   !> @param[in] q_ss_cart Cartesian spin-spiral q vectors for the sweep, one
   !>            column per q, in 2*pi/alat units. Row 1 is the reference point.
   !> @param[inout] lattice_obj Lattice object from prepare_post_processing_stack.
   !> @param[inout] hamiltonian_obj Hamiltonian object; q_ss is set to each
   !>               sweep point in turn.
   !> @param[inout] bands_obj Bands object used by the force-theorem probe.
   !> @param[inout] mix_obj Mixing object for the self-consistency runs.
   !> @param[inout] self_obj Self object; (re)constructed here and run once per
   !>               q in mode='scf'.
   !> @param[inout] energy_obj Energy object used by the force-theorem probe.
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
      ! GBT is a magnetic representation shared by both solvers.
      hamiltonian_obj%magnetic_representation = gbt_single_q
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
            ! WP8: the mesh must never be built once from row 1's q_ss and
            ! reused blindly for every other q in the sweep -- each policy
            ! below either shares one mesh proven valid for the whole sweep,
            ! or rebuilds per q inside the loop.
            select case (trim(recip_obj%q_symmetry_policy))
            case ('little_group_common')
               ! One mesh, reduced by the little group common to every q in
               ! the sweep (not just row 1's), valid for every probe below.
               call recip_obj%ensure_kpoint_mesh(recip_obj%nk_mesh, &
                                                 sum(abs(recip_obj%k_offset)) > 1.0e-12_rp, &
                                                 q_list_cart=q_ss_cart)
            case ('little_group')
               ! Built per q inside the loop below.
            case default   ! 'full_bz': WP0 oracle, unchanged from pre-WP8 behaviour.
               if (.not. allocated(recip_obj%k_points)) then
                  ! Every probe in one force-theorem difference must use the same
                  ! integration contract.  If the sweep contains finite q, make
                  ! the q=0 reference use the full chemical BZ as well.
                  if (any(abs(q_ss_cart) > 1.0e-12_rp)) then
                     recip_obj%use_symmetry_reduction = .false.
                     recip_obj%use_time_reversal = .false.
                     call recip_obj%generate_mp_mesh()
                  else if (recip_obj%use_symmetry_reduction) then
                     call recip_obj%generate_reduced_kpoint_mesh(recip_obj%nk_mesh, &
                                                                 sum(abs(recip_obj%k_offset)) > 1.0e-12_rp)
                  else
                     call recip_obj%generate_mp_mesh()
                  end if
               end if
            end select
         end if
         do iq = 1, fm_obj%n_q
            hamiltonian_obj%q_ss(:) = q_ss_cart(:, iq)
            if (use_kspace .and. trim(recip_obj%q_symmetry_policy) == 'little_group') then
               ! Rebuild (or, via the cache key, reuse if unchanged) for this
               ! specific q -- never reuse another q's little-group mesh.
               call recip_obj%ensure_kpoint_mesh(recip_obj%nk_mesh, sum(abs(recip_obj%k_offset)) > 1.0e-12_rp)
            end if
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
      hamiltonian_obj%magnetic_representation = gbt_single_q
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
      ! WP8: the mesh must never be built once (from whatever q_ss happens to be
      ! set here) and reused blindly for every other q in the sweep below --
      ! each policy either shares one mesh proven valid for the whole sweep, or
      ! rebuilds per q inside the loop.
      if (use_kspace) then
         recip_obj = reciprocal(hamiltonian_obj)
         select case (trim(recip_obj%q_symmetry_policy))
         case ('little_group_common')
            ! One mesh, reduced by the little group common to every q in the
            ! sweep (not just row 1's), valid for every probe below.
            call recip_obj%ensure_kpoint_mesh(recip_obj%nk_mesh, &
                                              sum(abs(recip_obj%k_offset)) > 1.0e-12_rp, &
                                              q_list_cart=q_ss_cart)
         case ('little_group')
            ! q_ss is 0 here (reset below); build the ordinary q=0 mesh now so
            ! the reference-energy probe has one. Per-nonzero-q rebuilds
            ! happen inside the sweep loop below via the same cache key.
            call recip_obj%ensure_kpoint_mesh(recip_obj%nk_mesh, sum(abs(recip_obj%k_offset)) > 1.0e-12_rp)
         case default   ! 'full_bz': WP0 oracle, unchanged from pre-WP8 behaviour.
            if (.not. allocated(recip_obj%k_points)) then
               if (any(abs(q_ss_cart) > 1.0e-12_rp)) then
                  recip_obj%use_symmetry_reduction = .false.
                  recip_obj%use_time_reversal = .false.
                  call recip_obj%generate_mp_mesh()
               else if (recip_obj%use_symmetry_reduction) then
                  call recip_obj%generate_reduced_kpoint_mesh(recip_obj%nk_mesh, &
                                                              sum(abs(recip_obj%k_offset)) > 1.0e-12_rp)
               else
                  call recip_obj%generate_mp_mesh()
               end if
            end if
         end select
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
         if (use_kspace .and. trim(recip_obj%q_symmetry_policy) == 'little_group') then
            ! Rebuild (or, via the cache key, reuse if unchanged) for this
            ! specific q -- never reuse another q's little-group mesh (WP8).
            call recip_obj%ensure_kpoint_mesh(recip_obj%nk_mesh, sum(abs(recip_obj%k_offset)) > 1.0e-12_rp)
         end if
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

   !> @brief Read the collinear reference axis for every sublattice type.
   !> @details The routine sets ref_theta to 0 (moment up) or pi (moment down)
   !>          from the sign of the converged z-moment. It sets ref_phi to 0
   !>          for every type. post_processing_frozen_magnon_auto uses these
   !>          reference angles as the zero point for its small cone-angle tilt
   !>          probes.
   !> @param[in] lattice_obj Lattice object. Reads symbolic_atoms(:)%potential%mom(3).
   !> @param[out] ref_theta Reference polar angle per sublattice type, radians.
   !> @param[out] ref_phi Reference azimuthal angle per sublattice type, radians.
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
   !>        (build/diagonalize -> normalized eigenvalue occupations and EBAND) or real-space recursion
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
            ! MFT must not depend on DOS window/grid/projections.  This solves
            ! EF from the target electron count and uses exactly those Fermi
            ! occupations in sum_k,n w_k f_nk epsilon_nk.
            e_band = recip_obj%calculate_canonical_band_energy(find_fermi=.true.)
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

   !> @brief Diagonalize the real symmetric magnon dynamical matrix.
   !> @details The routine calls LAPACK zheev on mat. It returns the magnon
   !>          branch energies in eval, ascending order, and the mode
   !>          eigenvectors in mat, in place. It raises a fatal error if mat is
   !>          not square, if eval does not match its dimension, or if zheev
   !>          reports a nonzero info.
   !> @param[inout] mat Magnon dynamical matrix on input; eigenvectors (mode
   !>               amplitudes and phases) on output.
   !> @param[out] eval Magnon branch energies, ascending order.
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

   !> @brief Write the metadata header lines for the auto branch-mode output files.
   !> @details The routine writes the sweep configuration (q-file coordinates,
   !>          theta_probe, active-moment threshold, reference energies) and
   !>          the active-sublattice list to frozen_magnon_branches.dat and
   !>          frozen_magnon_modes.dat. post_processing_frozen_magnon_auto
   !>          writes the per-q data rows after these headers.
   !> @param[in] branch_unit Open file unit for frozen_magnon_branches.dat.
   !> @param[in] modes_unit Open file unit for frozen_magnon_modes.dat.
   !> @param[in] fm_obj Frozen-magnon namelist settings.
   !> @param[in] q_ss_cart Cartesian spin-spiral q vectors for the sweep (only
   !>            its coordinate convention is recorded here).
   !> @param[in] theta_probe Cone-angle used for the finite-difference probes,
   !>            radians.
   !> @param[in] active_rec Sublattice index of each active sublattice.
   !> @param[in] active_type Atom-type index of each active sublattice.
   !> @param[in] active_moment Reference moment magnitude of each active
   !>            sublattice.
   !> @param[in] etot_ref Reference total energy at the collinear ground state.
   !> @param[in] eband_ref Reference band energy at the collinear ground state.
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
      this%gf_route = 'recursion'
      this%do_damping = .false.
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
          .and. post_processing /= 'bsf' &
          .and. post_processing /= 'density_of_states' &
          .and. post_processing /= 'kspace_green' &
          .and. post_processing /= 'frozen_magnon' &
          .and. post_processing /= 'susceptibility') then
         call g_logger%fatal('[calculation.check_post_processing]: '// &
                             "calculation%post_processing must be one of: ''none'', ''paoflow2rs'', ''exchange'', ''exchange_p2rs''," // &
                             " 'conductivity', 'conductivity_p2rs', 'orbital_modern', 'band_structure', 'bsf', 'density_of_states'," // &
                             " 'kspace_green', 'frozen_magnon', 'susceptibility'", __FILE__, __LINE__)
      end if
   end subroutine check_post_processing

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Check availability for the intersite Green's-function route (B2.5)
   !
   !> @param gf_route Route selector. Allowed: ´recursion´, ´lehmann´, ´dyson´
   !---------------------------------------------------------------------------
   subroutine check_gf_route(gf_route)
      character(len=*), intent(in) :: gf_route
      if (gf_route /= 'recursion' &
          .and. gf_route /= 'lehmann' &
          .and. gf_route /= 'dyson') then
         call g_logger%fatal('[calculation.check_gf_route]: '// &
                             "calculation%gf_route must be one of: 'recursion', 'lehmann', 'dyson'", &
                             __FILE__, __LINE__)
      end if
   end subroutine check_gf_route

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Check availability for pre-processing
   !
   !> @param[in] pre_processing Type of pre-processing. Allowed values:
   !> ´bravais´, ´buildsurf´, ´newclubulk´, ´newclusurf´, ´buildinterface´, ´none´
   !---------------------------------------------------------------------------
   subroutine check_pre_processing(pre_processing)
      character(len=*), intent(in) :: pre_processing
      if (pre_processing /= 'none' &
          .and. pre_processing /= 'bravais' &
          .and. pre_processing /= 'buildsurf' &
          .and. pre_processing /= 'newclubulk' &
          .and. pre_processing /= 'newclusurf' &
          .and. pre_processing /= 'buildinterface') then
         call g_logger%fatal("[calculation.check_pre_processing]:"// &
                             "calculation%pre_processing must be one of: ''none'', ''bravais'', "// &
                             "''buildsurf'', ''newclusurf'', ''newcluimp'', ''buildinterface''", __FILE__, __LINE__)
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
