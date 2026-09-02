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
   use kpoint_workset_mod, only: kpoint_workset
   use mix_mod
   use frozen_magnon_mod
   use vacuum_lead_mod, only: vacuum_lead, refresh_vacuum_region
   use math_mod
   use precision_mod, only: rp
   use string_mod, only: sl, fmt, real2str, int2str
   use timer_mod, only: g_timer
   use kpm_profile_mod, only: g_kpm_profile
   use logger_mod, only: g_logger
   use basis_mod, only: basis_init, norb
   use magnetic_representation_mod, only: periodic_nc, gbt_single_q
   use tddft_config_mod, only: tddft_config
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, tddft_chi0_batch_result, build_chi_ks_from_eigenpairs, &
      build_static_chi_ks_from_eigenpairs, write_chi_ks_text
   use tddft_xi_mod, only: tddft_direct_xi_result, build_direct_xi_from_operator_source, &
      build_static_direct_xi_from_operator_source
   use tddft_transition_engine_mod, only: pair_operator_tile_source
   use tddft_chi0_green_mod, only: green_chi0_options, eigenpair_green_function_provider, &
      build_chi_ks_from_green_functions, build_static_chi_ks_from_green_functions, &
      build_static_four_component_chi_ks_from_green_functions, build_four_component_chi_ks_from_green_functions
   use tddft_backend_mod, only: tddft_chi0_backend, tddft_eigenpair_backend, tddft_kspace_lehmann_backend, &
      tddft_realspace_gf_backend, canonical_tddft_backend_name, make_tddft_chi0_backend
   use tddft_chi0_realspace_mod, only: tddft_realspace_chi0_options, tddft_native_realspace_gf_provider
   use response_components_mod, only: RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel
   use tddft_four_component_mod, only: build_four_component_chi_ks, build_four_component_kernel, &
      evaluate_four_component_zero_modes, tddft_four_component_zero_mode_diagnostics
   use tddft_dyson_mod, only: tddft_dyson_options, tddft_dyson_result, enhance_tddft_susceptibility, &
      enhance_tddft_susceptibility_from_xi, write_tddft_dyson_text
   use tddft_goldstone_mod, only: tddft_goldstone_options, tddft_goldstone_result, evaluate_goldstone, &
      tddft_goldstone_diagnostics, tddft_goldstone_column_correction, evaluate_raw_xi_diagnostics, &
      build_goldstone_column_correction, spectral_weights_are_nonnegative, &
      spectral_weight_correction_is_acceptable, &
      write_goldstone_diagnostics_text, append_goldstone_column_correction_text
   use tddft_modes_mod, only: tddft_mode_options, tddft_mode_result, analyze_tddft_modes, write_tddft_modes_text
   use tddft_performance_mod, only: tddft_mpi_plan, make_tddft_mpi_plan
   use tddft_longitudinal_mod, only: build_charge_longitudinal_channels, build_charge_longitudinal_kernel, &
      append_longitudinal_response_metadata
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
   ! because the hook outlives `pre_processing_buildinterface`’s stack frame:
   ! it fires from inside `self%run()`, per SCF iteration.
   !
   ! Only ever populated on the `buildinterface` path with region_b_kind =
   ! ’vacuum’; a run that is not A | vacuum never installs the hook and never
   ! touches these.
   !---------------------------------------------------------------------------
   type(vacuum_lead), save :: vacuum_state
   type(self), pointer :: vacuum_solver => null()
   type(charge), pointer :: vacuum_charge => null()
   class(symbolic_atom), dimension(:), pointer :: vacuum_atoms => null()
   integer, save :: vacuum_nbulk_a = 0
   integer, save :: vacuum_nbulk = 0
   type(energy), pointer :: vacuum_energy => null()

   !> Streamed LMTO Q^- producer for TD-DFT pair-potential Xi.  It owns its
   !> small site-data vectors and two reusable construction matrices; the
   !> transition provider owns the one-k `(nmat,nmat,nright)` tile it fills.
   type, extends(pair_operator_tile_source) :: lmto_pair_operator_tile_source
      type(reciprocal), pointer :: reciprocal_obj => null()
      real(rp), allocatable :: signed_moments(:), column_scales(:)
      real(rp) :: q_point(3) = 0.0_rp
      complex(rp), allocatable :: qminus(:, :), qplus(:, :)
      logical :: use_qplus = .false.
      integer :: fetch_count = 0
   contains
      procedure :: initialize => initialize_lmto_pair_operator_source
      procedure :: channel_dimension => lmto_pair_operator_channel_dimension
      procedure :: fill_operator_tile => fill_lmto_pair_operator_tile
      procedure :: clear => clear_lmto_pair_operator_source
      final :: lmto_pair_operator_source_destructor
   end type lmto_pair_operator_tile_source

   type, public :: calculation
      !> Pre-processing. Options are:
      !> ´none´ (default)
      !> ´bravais´ : Builds the bulk clust
      !> ´buildsurf´ : Builds the surface clust
      !> ´newclubulk´ : Builds the imputiry clust from the bluk clust
      !> ´newclusurf´ : Builds the impurity clust from the surface clust
      !> ´buildinterface´ : Builds the two-sided (region A | active | region B)
      !>   layered/interface clust (calctype=’L’, B7.5)
      character(len=sl) :: pre_processing

      !> Processing. Options are
      !> ´none´ (default)
      character(len=sl) :: processing

      !> Post-processing. Options are
      !> ´none´ (default)
      character(len=sl) :: post_processing

      !> Green’s-function route for the intersite G_ij consumers (B2.5 dispatch).
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
      !> and, for the k-space routes, kspace_ham_order=’second’. Default .false.
      !> (bit-identical legacy: the damping routine is not invoked).
      logical :: do_damping

      !> Experimental magnetic moment-of-inertia switch for exchange post-processing.
      !> When true, evaluate the finite-difference torque-correlation kernel
      !> (`calculate_moment_of_inertia`) on the same canonical Green functions.
      !> This is diagnostic output only: no independent production relation or
      !> physical normalization is implied. Default .false.
      logical :: do_inertia

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
      procedure, private :: post_processing_fermi_surface
      procedure, private :: post_processing_kspace_green
      procedure, private :: post_processing_frozen_magnon
      procedure, private :: post_processing_susceptibility
      procedure :: process
      final :: destructor
   end type calculation

   interface calculation
      procedure :: constructor
   end interface calculation

   interface
      module subroutine pre_processing_bravais(this)
         class(calculation), intent(in) :: this
      end subroutine pre_processing_bravais

      module subroutine pre_processing_buildsurf(this)
         class(calculation), intent(in) :: this
      end subroutine pre_processing_buildsurf

      module subroutine pre_processing_newclubulk(this)
         class(calculation), intent(in) :: this
      end subroutine pre_processing_newclubulk

      module subroutine pre_processing_newclusurf(this)
         class(calculation), intent(in) :: this
      end subroutine pre_processing_newclusurf
   end interface

   interface
      module subroutine post_processing_band_structure(this)
         class(calculation), intent(in) :: this
      end subroutine post_processing_band_structure

      module subroutine post_processing_bsf(this)
         class(calculation), intent(in) :: this
      end subroutine post_processing_bsf

      module subroutine post_processing_density_of_states(this)
         class(calculation), intent(in) :: this
      end subroutine post_processing_density_of_states

      module subroutine post_processing_fermi_surface(this)
         class(calculation), intent(in) :: this
      end subroutine post_processing_fermi_surface

      module subroutine post_processing_kspace_green(this)
         class(calculation), intent(in) :: this
      end subroutine post_processing_kspace_green

      module subroutine prepare_post_processing_stack(this, use_paoflow, use_exchange_pairs, &
                                                       energy_mesh_before_hamiltonian, stochastic_moments, &
                                                       control_obj, lattice_obj, charge_obj, mix_obj, energy_obj, &
                                                       hamiltonian_obj, recursion_obj, dos_obj, green_obj, bands_obj, &
                                                       preprocessing_route, native_realspace_pairs, native_pair_rmax)
         class(calculation), intent(in) :: this
         logical, intent(in) :: use_paoflow, use_exchange_pairs
         logical, intent(in) :: energy_mesh_before_hamiltonian, stochastic_moments
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
         character(len=*), intent(in), optional :: preprocessing_route
         logical, intent(in), optional :: native_realspace_pairs
         real(rp), intent(in), optional :: native_pair_rmax
      end subroutine prepare_post_processing_stack

      module subroutine run_intersite_moments(control_obj, recursion_obj)
         type(control), intent(in) :: control_obj
         type(recursion), intent(inout) :: recursion_obj
      end subroutine run_intersite_moments
   end interface

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
      pre_processing = this%pre_processing(:len(pre_processing))
      processing = this%processing(:len(processing))
      post_processing = this%post_processing(:len(post_processing))
      gf_route = this%gf_route(:len(gf_route))
      do_damping = this%do_damping
      do_inertia = this%do_inertia

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
      ! Green’s-function route (B2.5 dispatch)
      call check_gf_route(trim(gf_route))

      this%verbose = verbose
      this%fname = fname
      this%pre_processing = pre_processing
      this%processing = processing
      this%post_processing = post_processing
      this%gf_route = gf_route
      this%do_damping = do_damping
      this%do_inertia = do_inertia

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
      case ('fermi_surface')
         call this%post_processing_fermi_surface()
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

      ! The selected preprocessing route has already converged and written the
      ! atom state.  Rebuild only the consumer stack from that route; do not
      ! manufacture a surface/impurity geometry here.
      call prepare_post_processing_stack(this, .false., .false., .true., .false., control_obj, lattice_obj, charge_obj, mix_obj, &
                                         energy_obj, hamiltonian_obj, recursion_obj, dos_obj, green_obj, bands_obj, &
                                         this%pre_processing)

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
   !> Pre-process for the layered/interface calculation (B7.5): region A |
   !> active zone | region B, calctype=’L’. Mirrors pre_processing_buildsurf;
   !> buildsurf itself is untouched and remains the permanent regression
   !> oracle for the one-sided (vacuum|active|bulk) case (B7 §4 B7.5).
   !>
   !> Workflow chain: bulk-A -> bulk-B (-> vacuum generator) -> interface,
   !> mirroring bulk -> surf -> imp. Region A’s and region B’s converged
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
      ! the loaded symbolic atoms’ occupation vs valence -- imppot’s
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

      ! B7.6, A | vacuum: region B’s frozen parameters are GENERATED, not read
      ! from an &atoms label. Install the refresh hook and run it once now, so
      ! the very first recursion sees real empty-lattice parameters rather than
      ! the defaults the unlabelled type was constructed with. The initial call
      ! runs at region_shift = 0, i.e. the vacuum level is the anchor’s zero;
      ! every later call, driven from interfacepot, uses the solved level. That
      ! is the ″generate once″ behaviour arriving as iteration 0 of the
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
   !> Argument-free by the hook’s interface contract; the data it needs is the
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


   !> @brief Finish the Chebyshev DOS step and resolve the Fermi level.
   !> @details The routine calls green%chebyshev_dos_dispatch to turn the
   !>          Chebyshev transport moments into a Green’s function. It then
   !>          calls bands%calculate_fermi. Both conductivity drivers
   !>          (post_processing_conductivity and post_processing_conductivity_p2rs)
   !>          call this routine right before they build the conductivity tensor.
   !> @param[inout] green_obj Green object. The Chebyshev DOS dispatch fills its
   !>               Green’s-function arrays.
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
   !>          eta-broadened intersite Green’s function, and evaluates J_ij
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
   !>          this%gf_route. Route ’recursion’ runs the real-space intersite
   !>          recursion pass. Routes ’lehmann’ and ’dyson’ fill the same
   !>          arrays from the k-space engine (reciprocal%fill_green). Every
   !>          route fills the same canonical arrays, so
   !>          calculate_intersite_gf_twoindex and the exchange evaluation run
   !>          unchanged after this point. When this%do_damping is true, the
   !>          routine also evaluates the on-site Kambersky torque-correlation
   !>          Gilbert damping from the same green%gij.
   !> @param[in] this Calculation object. Reads fname, gf_route, and do_damping.
   !> @note The damping term needs SOC in the potential (xi_p/xi_d). The
   !>       k-space routes also need kspace_ham_order=’second’.
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
         ! former ″1/sqrt2″ worry is resolved: the fixed-broadening J difference vs
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
      ! kspace_ham_order=’second’.
      if (this%do_damping) then
         call exchange_obj%calculate_gilbert_damping()
      end if
      if (this%do_inertia) then
         call exchange_obj%calculate_moment_of_inertia()
      end if
   end subroutine post_processing_exchange


   !> @brief Conductivity-tensor post-processing (B5.1 dispatch).
   !> @details The routine picks the source of the Chebyshev transport moments
   !>          (mu_nm_stochastic) from this%gf_route. Route ’recursion’ runs
   !>          the stochastic moment generator inside prepare_post_processing_stack.
   !>          Routes ’lehmann’ and ’dyson’ fill the same array exactly from the
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

      call g_kpm_profile%reset()
      call g_kpm_profile%start('T_transport_total')

      ! Route the Chebyshev transport moments (mu_nm_stochastic) through the
      ! selected producer (B5.1 dispatch). ’recursion’ (default) fills them
      ! stochastically inside the prepared stack -- the bit-identical legacy
      ! path. ’lehmann’/’dyson’ fill the SAME array EXACTLY from the k-space
      ! eigenpairs (reciprocal%fill_moments), so calculate_conductivity_tensor
      ! runs unchanged. exact-vs-recursion on the same crystal is the direct KPM
      ! error bound (see docs/dev/route_agnostic_estimators.md).
      rec_moments = (trim(this%gf_route) == 'recursion')
      call g_kpm_profile%start('P_stack_setup')
      call prepare_post_processing_stack(this, .false., .false., .true., .false., control_obj, lattice_obj, &
                                         charge_obj, mix_obj, energy_obj, hamiltonian_obj, recursion_obj, dos_obj, &
                                         green_obj, bands_obj)
      call g_kpm_profile%stop('P_stack_setup')
      if (rec_moments) call recursion_obj%compute_moments_stochastic()
      if (.not. rec_moments) then
         ! Exact k-space moments on the SAME Chebyshev window the recursion route
         ! and gamma_nm use (resolve_chebyshev_window -> a, b). The reciprocal
         ! object is kept at subroutine scope (spglib-owning finalizer) like the
         ! exchange/BSF drivers. Sigma=0: ’lehmann’ and ’dyson’ coincide here.
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
      call g_kpm_profile%start('P_moment_finalize')
      call finish_conductivity_moments(green_obj, bands_obj)
      call g_kpm_profile%stop('P_moment_finalize')
      conductivity_obj = conductivity(self_obj)
      call conductivity_obj%calculate_gamma_nm()
      call conductivity_obj%calculate_conductivity_tensor()
      call g_kpm_profile%stop('T_transport_total')
      call g_kpm_profile%emit()
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

      call g_kpm_profile%reset()
      call g_kpm_profile%start('T_transport_total')

      call g_kpm_profile%start('P_stack_setup')
      call prepare_post_processing_stack(this, .true., .false., .true., .false., control_obj, lattice_obj, &
                                         charge_obj, mix_obj, energy_obj, hamiltonian_obj, recursion_obj, dos_obj, &
                                         green_obj, bands_obj)
      call g_kpm_profile%stop('P_stack_setup')
      call recursion_obj%compute_moments_stochastic()
      self_obj = self(bands_obj, mix_obj)
      call g_kpm_profile%start('P_moment_finalize')
      call finish_conductivity_moments(green_obj, bands_obj)
      call g_kpm_profile%stop('P_moment_finalize')
      conductivity_obj = conductivity(self_obj)
      call conductivity_obj%calculate_gamma_nm()
      call conductivity_obj%calculate_conductivity_tensor()
      call g_kpm_profile%stop('T_transport_total')
      call g_kpm_profile%emit()
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
         ! B7.5: two-sided counterpart of ’S’.
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
         ! B7.5: ’S’ matrices plus region reference charges and the two-sided
         ! registry overwriting surfmat’s one-sided one.
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
         ! B7.5: identical to ’S’; no build_locham (nmax = 0 for ’L’).
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

   ! DESCRIPTION:
   !> @brief
   !> Frozen-magnon post-processing: sweeps &hamiltonian’s q_ss over a
   !> user-supplied list of points (see &frozen_magnon namelist) and reports
   !> the adiabatic magnon dispersion.
   !> @details For each q in the sweep, the Hamiltonian is rebuilt with that
   !> q_ss and the potential is either fully re-converged (mode=’scf’) or, for
   !> mode=’mft’ (the magnetic force theorem default), converged once at the
   !> reference point q_ss_list(:,1) and then re-evaluated with a single
   !> band-energy pass (self%nstep=1) at every other q, reusing the reference
   !> potential unchanged. The dispersion omega(q) = 4*(E(q)-E(q_ref)) /
   !> (M_tot*sin^2(theta_ss)) uses the band energy for mode=’mft’ and the
   !> total energy for mode=’scf’. This is the single-acoustic-branch
   !> generalization of the Halilov frozen-magnon formula to multiple
   !> sublattices: M_tot is the reference point’s total moment summed over
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
         if (fm_obj%mode == 'mft_constrained') then
            call g_logger%fatal('[calculation.post_processing_frozen_magnon]: mode=''mft_constrained'' is currently supported for branch_mode=''acoustic'' only.', &
                               __FILE__, __LINE__)
         end if
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
      type(reciprocal), target :: reciprocal_obj
      type(tddft_config) :: config
      type(tddft_chi0_options) :: chi0_options, chi0_options_reverse
      type(green_chi0_options) :: green_options, green_options_reverse
      type(eigenpair_green_function_provider), target :: green_source
      type(tddft_realspace_chi0_options) :: realspace_options, realspace_options_reverse
      type(tddft_native_realspace_gf_provider), target :: realspace_source, realspace_source_reverse
      type(tddft_chi0_result) :: chi0_result, chi0_result_reverse, chi0_static, chi0_static_reverse
      type(tddft_chi0_batch_result) :: realspace_batch, realspace_batch_reverse
      type(tddft_dyson_options) :: dyson_options, dyson_options_reverse
      type(tddft_dyson_result) :: dyson_result, dyson_pair_result, dyson_pair_corrected_result
      type(tddft_dyson_result) :: dyson_result_reverse, dyson_pair_result_reverse, dyson_pair_corrected_result_reverse
      type(tddft_direct_xi_result) :: pair_xi_result, pair_xi_static, pair_xi_corrected_result
      type(tddft_direct_xi_result) :: pair_xi_result_reverse, pair_xi_static_reverse, pair_xi_corrected_result_reverse
      type(tddft_goldstone_options) :: goldstone_options
      type(tddft_goldstone_result) :: goldstone_result, goldstone_result_reverse
      type(tddft_goldstone_diagnostics) :: pair_goldstone, pair_goldstone_reverse
      type(tddft_goldstone_column_correction) :: pair_correction
      type(tddft_mode_options) :: mode_options
      type(tddft_mode_result) :: mode_result
      type(tddft_four_component_zero_mode_diagnostics) :: full_zero_mode_diagnostics
      class(tddft_chi0_backend), allocatable :: chi0_backend
      type(tddft_mpi_plan) :: tddft_plan
      type(lmto_pair_operator_tile_source), target :: pair_operator_source
      type(response_channel), allocatable :: left_channels(:), right_channels(:), left_channels_reverse(:), right_channels_reverse(:)
      real(rp), allocatable :: omega(:), omega_static(:), eigenvalues_k(:, :), eigenvalues_kq(:, :)
      type(kpoint_workset) :: kq_workset
      complex(rp), allocatable :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :), kernel(:, :), all_xi(:, :, :, :), all_loss(:, :, :, :)
      complex(rp), allocatable :: all_xi_pair(:, :, :, :), all_loss_pair(:, :, :, :)
      real(rp), allocatable :: all_trace_loss(:, :), all_trace_loss_pair(:, :), all_trace_loss_reverse(:, :), all_trace_loss_pair_reverse(:, :)
      real(rp), allocatable :: coulomb_site(:, :), magnetization(:, :), site_moments(:, :)
      real(rp), allocatable :: signed_moments(:)
      complex(rp), allocatable :: all_xi_reverse(:, :, :, :), all_loss_reverse(:, :, :, :)
      complex(rp), allocatable :: all_xi_pair_reverse(:, :, :, :), all_loss_pair_reverse(:, :, :, :)
      real(rp) :: response_eta, t_profile_start, t_profile_stop, kq_eigensolve_cpu_seconds
      real(rp) :: response_electron_count, response_band_energy, electron_count_tolerance
      real(rp) :: bare_gamma_peak, legacy_gamma_peak, pair_gamma_peak, pair_corrected_gamma_peak
      real(rp) :: bare_gamma_peak_reverse, legacy_gamma_peak_reverse, pair_gamma_peak_reverse
      real(rp) :: pair_corrected_gamma_peak_reverse
      real(rp) :: raw_pair_minimum_spectral_weight, corrected_pair_minimum_spectral_weight
      integer, allocatable :: site_orbital_counts(:)
      integer :: iq, iq_start, iq_end, nq_per_rank, nq, nw, unit, ios, isite, nresponse
      integer :: corrected_minimum_location(2)
      logical :: has_soc, has_external_field, need_dyson, is_longitudinal, is_full_response, is_gamma, has_gamma
      logical :: pair_backend, legacy_backend, raw_pair_spectral_weight_ok, corrected_spectral_weight_ok, &
         correction_spectral_weight_ok, full_response_supported, circular_reverse, primary_minus_plus
      character(len=sl) :: filename, chi0_filename, legacy_filename, pair_filename
      character(len=256) :: full_response_capability_reason
      character(len=32) :: canonical_chi0_backend
      character(len=384) :: electron_count_message
      character(len=640) :: spectral_weight_message

      config = tddft_config(this%fname)
      if (.not. config%enabled) then
         if (rank == 0) call g_logger%info('TDDFT susceptibility route disabled by &tddft enabled=.false.', __FILE__, __LINE__)
         return
      end if
      is_longitudinal = config%channel == 'longitudinal'
      is_full_response = config%channel == 'full'
      canonical_chi0_backend = canonical_tddft_backend_name(config%chi0_backend)
      if (len_trim(canonical_chi0_backend) == 0) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: unknown chi0 backend '''// &
            trim(config%chi0_backend)//'''.', __FILE__, __LINE__)
      end if
      call prepare_post_processing_stack(this, .false., .false., .true., .false., control_obj, lattice_obj, &
         charge_obj, mix_obj, energy_obj, hamiltonian_obj, recursion_obj, dos_obj, green_obj, bands_obj, &
         native_realspace_pairs=(canonical_chi0_backend == 'realspace_gf'), &
         native_pair_rmax=config%realspace_rmax)
      if (control_obj%calctype /= 'B') then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: eigenpair TDDFT currently requires calctype=''B''.', &
            __FILE__, __LINE__)
      end if

      ! TDDFT-03 establishes one deliberately narrow response boundary.  The
      ! transverse eigenpair formula below is derived for the collinear,
      ! scalar-relativistic, orthogonal Hamiltonian only.  Reject every
      ! branch whose response operator or metric derivative has not been
      ! derived, before building eigenpairs or accepting a kernel.
      if (control_obj%is_noncollinear()) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: noncollinear TDDFT is unsupported: '// &
            'nsp=3/4 requires the full spinor response, torque terms, and a noncollinear kernel.', __FILE__, __LINE__)
      end if
      if (control_obj%has_soc()) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: relativistic TDDFT is unsupported: '// &
            'nsp=2/4 requires SOC response, anisotropy, and torque terms.', __FILE__, __LINE__)
      end if
      if (control_obj%nsp /= 1) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: eigenpair TDDFT baseline requires nsp=1 '// &
            '(collinear, SOC-free); relativistic and noncollinear response is not silently approximated.', __FILE__, __LINE__)
      end if
      if (control_obj%do_comom .or. control_obj%constraints_enable) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: eigenpair TDDFT baseline rejects common-moment or '// &
            'constrained external fields until their response derivatives are implemented.', __FILE__, __LINE__)
      end if
      if (trim(hamiltonian_obj%magnetic_representation) /= periodic_nc) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: eigenpair TDDFT baseline requires magnetic_representation=periodic_nc; '// &
            'GBT and explicit-texture response states are unsupported.', __FILE__, __LINE__)
      end if
      if (hamiltonian_obj%local_axis) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: local-axis Hamiltonian response is unsupported; '// &
            'the response basis must remain in the global spin frame.', __FILE__, __LINE__)
      end if
      if (hamiltonian_obj%hoh) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: eigenpair TDDFT baseline rejects HOH/second-order Hamiltonians.', &
            __FILE__, __LINE__)
      end if
      if (hamiltonian_obj%orb_pol) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: eigenpair TDDFT baseline rejects orbital-polarization Hamiltonians '// &
            'until the response operator is differentiated.', __FILE__, __LINE__)
      end if
      if (hamiltonian_obj%ccor_2c .or. abs(hamiltonian_obj%ccor_elin) > tiny(1.0_rp)) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: eigenpair TDDFT baseline rejects CCOR-modified Hamiltonians '// &
            'until the response derivative is implemented.', __FILE__, __LINE__)
      end if
      if (hamiltonian_obj%hubbard_u_general_check .or. hamiltonian_obj%hubbard_u_impurity_check .or. &
          hamiltonian_obj%hubbard_u_sc_check .or. hamiltonian_obj%hubbard_v_check) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: eigenpair TDDFT baseline rejects Hubbard-corrected Hamiltonians '// &
            'until the Hubbard response is implemented.', __FILE__, __LINE__)
      end if

      has_soc = .false.
      do isite = 1, lattice_obj%ntype
         has_soc = has_soc .or. any(abs(lattice_obj%symbolic_atoms(isite)%potential%xi_p) > tiny(1.0_rp)) .or. &
            any(abs(lattice_obj%symbolic_atoms(isite)%potential%xi_d) > tiny(1.0_rp))
      end do
      if (has_soc) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: eigenpair TDDFT baseline rejects nonzero SOC '// &
            'until the relativistic response and anisotropy terms are derived.', __FILE__, __LINE__)
      end if

      ! The response driver may be invoked after a normal SCF object has gone
      ! out of scope.  Refresh the provider from the same VXC0SP route, then
      ! rebuild the normal-state Hamiltonian from that ground-state potential.
      call get_mpi_variables(rank, lattice_obj%nrec)
      self_obj = self(bands_obj, mix_obj)
      self_obj%xc_response_derivatives_requested = is_longitudinal
      call self_obj%refresh_xc_response_kernel()
      if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham()
      call hamiltonian_obj%build_bulkham()
      reciprocal_obj = reciprocal(hamiltonian_obj)
      if (trim(reciprocal_obj%reciprocal_mode) /= 'ham_only') then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: eigenpair TDDFT baseline requires reciprocal_mode=ham_only; '// &
            'generalized-overlap response is unsupported.', __FILE__, __LINE__)
      end if
      if (trim(reciprocal_obj%kspace_ham_order) == 'second' .or. reciprocal_obj%include_so) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: eigenpair TDDFT baseline requires first-order, SOC-free '// &
            'Hamiltonian assembly.', __FILE__, __LINE__)
      end if
      ! Response uses a complete reciprocal mesh.  A reduced mesh cannot in
      ! general be paired with k+q without response-specific symmetry weights.
      reciprocal_obj%use_symmetry_reduction = .false.
      reciprocal_obj%use_time_reversal = .false.
      call reciprocal_obj%generate_mp_mesh()
      ! TDDFT inherits the canonical reciprocal occupation contract.  Its
      ! Fermi level is resolved only after complete response-mesh eigenpairs
      ! exist, so a coarser SCF-mesh value cannot change the response count.
      config%ground_state_fermi_level = reciprocal_obj%fermi_level
      config%ground_state_electronic_temperature = reciprocal_obj%temperature
      config%ground_state_electron_count = reciprocal_obj%total_electrons
      if (.not. config%electronic_temperature_overridden) then
         config%electronic_temperature = config%ground_state_electronic_temperature
      end if
      if (config%electronic_temperature < 0.0_rp) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: response temperature is unresolved.', __FILE__, __LINE__)
      end if
      reciprocal_obj%temperature = config%electronic_temperature
      circular_reverse = .not. is_longitudinal .and. .not. is_full_response .and. &
         trim(config%circular_channel) == 'both'
      primary_minus_plus = trim(config%circular_channel) == 'minus_plus'
      pair_backend = config%xi_backend == 'pair_potential' .or. config%xi_backend == 'compare'
      legacy_backend = config%xi_backend == 'legacy_site_scalar' .or. config%xi_backend == 'compare'
      if (canonical_chi0_backend == 'realspace_gf' .and. (is_longitudinal .or. is_full_response)) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: native real-space GF currently supports the '// &
            'collinear transverse response only; longitudinal/full response combinations remain explicitly unsupported.', __FILE__, __LINE__)
      end if
      if (canonical_chi0_backend == 'realspace_gf' .and. control_obj%recur /= 'block' .and. control_obj%recur /= 'chebyshev') then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: native real-space GF requires block or chebyshev '// &
            'intersite recursion storage; the lanczos route has no native G(R,z) source.', __FILE__, __LINE__)
      end if
      if (config%goldstone_policy /= 'diagnose' .and. (is_longitudinal .or. is_full_response .or. .not. legacy_backend)) then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: goldstone_policy sum_rule/projected currently requires '// &
            'the transverse legacy active response basis.', __FILE__, __LINE__)
      end if
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
      nw = config%nomega
      tddft_plan = make_tddft_mpi_plan(canonical_chi0_backend, config%gf_integration, nq, nw, &
         product(reciprocal_obj%nk_mesh), lattice_obj%njij, rank, numprocs, .true., config%green_energy_points)
      ! The planner owner range is used for q-labelled output.  The native
      ! R-GF compute range is intentionally wider: all q points are transformed
      ! after the local R blocks have been reduced, preserving q amortization.
      iq_start = tddft_plan%owner_q%first
      iq_end = tddft_plan%owner_q%last
      nq_per_rank = tddft_plan%owner_q%count
      if (rank == 0) call tddft_plan%emit(label='production')
      allocate(omega(nw))
      if (nw == 1) then
         omega(1) = config%omega_min
      else
         do iq = 1, nw
            omega(iq) = config%omega_min + real(iq-1, rp)*(config%omega_max-config%omega_min)/real(nw-1, rp)
         end do
      end if

      if (is_longitudinal) then
         nresponse = 2*lattice_obj%nrec
      else if (is_full_response) then
         nresponse = 4*lattice_obj%nrec
      else
         nresponse = lattice_obj%nrec
      end if
      allocate(site_orbital_counts(lattice_obj%nrec), left_channels(nresponse), right_channels(nresponse), &
         left_channels_reverse(lattice_obj%nrec), right_channels_reverse(lattice_obj%nrec))
      site_orbital_counts = norb
      if (is_longitudinal) then
         call build_charge_longitudinal_channels(lattice_obj%nrec, left_channels)
         right_channels = left_channels
      else
         do iq = 1, lattice_obj%nrec
            if (primary_minus_plus) then
               left_channels(iq) = response_channel(iq, RESPONSE_MINUS)
               right_channels(iq) = response_channel(iq, RESPONSE_PLUS)
               left_channels_reverse(iq) = response_channel(iq, RESPONSE_PLUS)
               right_channels_reverse(iq) = response_channel(iq, RESPONSE_MINUS)
            else
               left_channels(iq) = response_channel(iq, RESPONSE_PLUS)
               right_channels(iq) = response_channel(iq, RESPONSE_MINUS)
               left_channels_reverse(iq) = response_channel(iq, RESPONSE_MINUS)
               right_channels_reverse(iq) = response_channel(iq, RESPONSE_PLUS)
            end if
         end do
      end if
      chi0_options%eta = config%eta
      chi0_options%electronic_temperature = config%electronic_temperature
      chi0_options%band_first = config%band_first
      chi0_options%band_last = config%band_last
      chi0_options%occupation_prune_tolerance = config%occupation_tolerance
      chi0_options%k_mesh_shape = reciprocal_obj%nk_mesh
      chi0_options%response_projection = config%response_projection
      if (is_longitudinal) then
         chi0_options%circular_channel = 'charge_mz'
      else if (primary_minus_plus) then
         chi0_options%circular_channel = 'minus_plus'
      else
         chi0_options%circular_channel = 'plus_minus'
      end if
      chi0_options_reverse = chi0_options
      if (primary_minus_plus) then
         chi0_options_reverse%circular_channel = 'plus_minus'
      else
         chi0_options_reverse%circular_channel = 'minus_plus'
      end if
      green_options%eta = config%eta
      green_options%green_eta = config%green_eta
      green_options%electronic_temperature = config%electronic_temperature
      green_options%energy_min = config%green_energy_min
      green_options%energy_max = config%green_energy_max
      green_options%energy_points = config%green_energy_points
      green_options%energy_integration = config%gf_integration
      green_options%contour_points = config%contour_points
      green_options%contour_subdivisions = config%contour_subdivisions
      green_options%near_fermi_points = config%near_fermi_points
      green_options%contour_height = config%contour_height
      green_options%k_mesh_shape = reciprocal_obj%nk_mesh
      green_options%response_projection = config%response_projection
      green_options%circular_channel = chi0_options%circular_channel
      green_options_reverse = green_options
      green_options_reverse%circular_channel = chi0_options_reverse%circular_channel
      realspace_options%eta = config%eta
      realspace_options%green_eta = config%green_eta
      realspace_options%energy_integration = config%gf_integration
      realspace_options%electronic_temperature = config%electronic_temperature
      realspace_options%rmax = config%realspace_rmax
      realspace_options%tail_tolerance = config%realspace_tail_tolerance
      realspace_options%representation = config%realspace_representation
      realspace_options%fourier_axes = config%realspace_fourier_axes
      realspace_options%circular_channel = chi0_options%circular_channel
      realspace_options_reverse = realspace_options
      realspace_options_reverse%circular_channel = chi0_options_reverse%circular_channel
      if (canonical_chi0_backend == 'realspace_gf' .and. config%gf_integration /= 'direct') then
         call g_logger%fatal('[calculation.post_processing_susceptibility]: native real-space GF currently exposes sampled real-axis blocks only; '// &
            "gf_integration='mixed_contour' requires a complex-energy native source.", __FILE__, __LINE__)
      end if
      has_external_field = control_obj%do_comom .or. control_obj%constraints_enable
      if (canonical_chi0_backend == 'kspace_lehmann' .and. config%chi0_backend == 'green' .and. rank == 0) then
         call g_logger%warning('[calculation.post_processing_susceptibility]: chi0_backend=''green'' currently selects the '// &
            'K-space Lehmann backend backed by the eigenpair-resolvent reference, not a native RS Green-function provider; '// &
            'native RS response is selected with chi0_backend=''realspace_gf''.', &
            __FILE__, __LINE__)
      end if
      if (.not. is_full_response) call make_tddft_chi0_backend(config%chi0_backend, chi0_backend)

      if (canonical_chi0_backend == 'realspace_gf') then
         ! The common response setup above maps ownership to response sites so
         ! k-space spin moments and eigenpairs use the normal nrec layout.  The
         ! native source instead consumes every generated pair; restore the
         ! pair mapping immediately before intersite recursion.  This is a
         ! no-op in serial execution except for extending g2l_map from nrec to
         ! njij, but it prevents all non-self real-space pairs from remaining
         ! silently zero.
         call get_mpi_variables(rank, lattice_obj%njij)
         call run_intersite_moments(control_obj, recursion_obj)
         call green_obj%calculate_intersite_gf()
         call realspace_source%initialize_from_green(green_obj, lattice_obj, site_orbital_counts, left_channels, &
            right_channels, realspace_options)
         if (circular_reverse) then
            call realspace_source_reverse%initialize_from_green(green_obj, lattice_obj, site_orbital_counts, &
               left_channels_reverse, right_channels_reverse, realspace_options_reverse)
         end if
         select type (chi0_backend)
         type is (tddft_realspace_gf_backend)
            call chi0_backend%initialize(realspace_source)
         class default
            call g_logger%fatal('[calculation.post_processing_susceptibility]: native real-space provider attachment failed.', &
               __FILE__, __LINE__)
         end select
         call get_mpi_variables(rank, lattice_obj%nrec)
      end if

      ! k eigenpairs are independent of q and are therefore reused on each q
      ! worker.  k+q eigenpairs remain caller-owned and exact at off-mesh q.
      call reciprocal_obj%calculate_eigenpairs_at_kpoints(reciprocal_obj%k_workset%points, eigenvalues_k, eigenvectors_k)

      ! `*_out.nml` restart files retain the potential and its direction but
      ! not the scalar site moment `mtot`.  The XC radial projection recorded
      ! by refresh_xc_response_kernel must be normalized by the same occupied
      ! P_site sigma population used by the response vertices.  Reconstruct it
      ! from the complete, unreduced response mesh rather than relying on that
      ! non-serialized legacy cache.
      reciprocal_obj%eigenvalues = eigenvalues_k
      reciprocal_obj%eigenvectors = eigenvectors_k
      ! Resolve the chemical potential on the actual complete response mesh
      ! after its eigenpairs exist.  There is deliberately no &tddft EF input:
      ! the response remains at the reciprocal ground-state electron count.
      response_band_energy = reciprocal_obj%calculate_canonical_band_energy(find_fermi=.true., &
         electron_count=response_electron_count)
      config%fermi_level = reciprocal_obj%fermi_level
      config%response_electron_count = response_electron_count
      chi0_options%fermi_level = config%fermi_level
      green_options%fermi_level = config%fermi_level
      if (canonical_chi0_backend == 'realspace_gf') then
         realspace_source%options%fermi_level = config%fermi_level
         realspace_source%options%electronic_temperature = config%electronic_temperature
         if (circular_reverse) then
            realspace_source_reverse%options%fermi_level = config%fermi_level
            realspace_source_reverse%options%electronic_temperature = config%electronic_temperature
         end if
         select type (chi0_backend)
         type is (tddft_realspace_gf_backend)
            ! The backend owns a provider copy.  Reattach after the response
            ! mesh resolves EF so the native bubble uses the same occupations
            ! as the ground-state response contract.
            call chi0_backend%initialize(realspace_source)
         class default
            call g_logger%fatal('[calculation.post_processing_susceptibility]: native real-space backend reattachment failed.', &
               __FILE__, __LINE__)
         end select
      end if
      electron_count_tolerance = 1.0e-8_rp*max(1.0_rp, config%ground_state_electron_count)
      if (abs(response_electron_count-config%ground_state_electron_count) > electron_count_tolerance) then
         write(electron_count_message, '(a,es16.8,a,es16.8,a,es12.4,a,es16.8,a,es16.8,a)') &
            '[calculation.post_processing_susceptibility]: response electron count does not match target: target=', &
            config%ground_state_electron_count, ', recomputed=', response_electron_count, ', dN=', &
            response_electron_count-config%ground_state_electron_count, ', ground_EF=', config%ground_state_fermi_level, &
            ', response_EF=', config%fermi_level, ', response EF is derived from the response mesh'
         call g_logger%fatal(trim(electron_count_message)// &
            '. Check reciprocal total_electrons and whether the response band window can represent the target count.', __FILE__, __LINE__)
      end if
      allocate(site_moments(3, lattice_obj%nrec), signed_moments(lattice_obj%nrec))
      call self_obj%compute_kspace_spin_moments_spinor(reciprocal_obj, site_moments)
      do isite = 1, lattice_obj%nrec
         ! A transverse Goldstone vector is signed.  Replacing this by its
         ! magnitude breaks reversed and multi-sublattice reference states.
         call self_obj%xc_response_provider%set_site_spin_population(isite, abs(site_moments(3, isite)))
         call self_obj%xc_response_provider%set_site_signed_spin_population(isite, site_moments(3, isite))
         signed_moments(isite) = site_moments(3, isite)
         if (sqrt(sum(site_moments(:, isite)**2)) > tiny(1.0_rp)) then
            call self_obj%xc_response_provider%set_site_magnetization_direction(isite, site_moments(:, isite))
         end if
      end do
      deallocate(site_moments)
      if (is_full_response) then
         call self_obj%xc_response_provider%full_response_capability(full_response_supported, full_response_capability_reason)
         if (.not. full_response_supported) then
            call g_logger%fatal('[calculation.post_processing_susceptibility]: channel=''full'' is unavailable for XC route '''// &
               trim(self_obj%xc_response_provider%functional_label)//''': '//trim(full_response_capability_reason), __FILE__, __LINE__)
         end if
      end if

      ! The Gamma/static KS batch is shared by both channels.  Transverse uses
      ! it for Goldstone diagnostics; longitudinal uses it only as the exact
      ! zero-frequency bare response reference for the common Dyson path.
      allocate(omega_static(1))
      omega_static = 0.0_rp
      ! A native real-space source is dynamic-first in TDDFT-07.  Do not
      ! manufacture its static Ward input from reciprocal eigenpairs; the
      ! enhanced/Dyson route remains an explicit future provider requirement.
      need_dyson = config%output_xi .or. config%output_chi .or. config%output_modes .or. is_longitudinal .or. pair_backend
      if (canonical_chi0_backend == 'realspace_gf') then
         if (need_dyson) then
            call g_logger%fatal('[calculation.post_processing_susceptibility]: native real-space GF dynamic output currently '// &
               'does not provide an exact static kernel for Xi/Dyson enhancement; request bare chi0 output only.', __FILE__, __LINE__)
         end if
      else if (.not. is_longitudinal .and. .not. is_full_response) then
         if (canonical_chi0_backend /= 'eigenpairs' .and. canonical_chi0_backend /= 'kspace_lehmann') then
            call g_logger%fatal('[calculation.post_processing_susceptibility]: real static Ward diagnostics require '// &
               'an eigenpair or K-space Lehmann backend with an exact static-limit solver.', __FILE__, __LINE__)
         end if
         chi0_options%q_direct = 0.0_rp
         green_options%q_direct = 0.0_rp
         if (canonical_chi0_backend == 'kspace_lehmann') then
            call green_source%initialize(eigenvalues_k, eigenvectors_k, eigenvalues_k, eigenvectors_k)
            call build_static_chi_ks_from_green_functions(green_source, reciprocal_obj%k_weights, site_orbital_counts, &
               left_channels, right_channels, green_options, chi0_static)
         else
            call build_static_chi_ks_from_eigenpairs(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, &
               site_orbital_counts, left_channels, right_channels, chi0_options, chi0_static)
         end if
      else if (canonical_chi0_backend == 'kspace_lehmann') then
         call green_source%initialize(eigenvalues_k, eigenvectors_k, eigenvalues_k, eigenvectors_k)
         if (is_full_response) then
            call build_static_four_component_chi_ks_from_green_functions(green_source, reciprocal_obj%k_weights, &
               site_orbital_counts, green_options, chi0_static)
         else
            call build_static_chi_ks_from_green_functions(green_source, reciprocal_obj%k_weights, site_orbital_counts, &
               left_channels, right_channels, green_options, chi0_static)
         end if
      else if (is_full_response) then
         call build_four_component_chi_ks(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_k, &
            eigenvectors_k, site_orbital_counts, omega_static, chi0_options, chi0_static)
      else
         call build_chi_ks_from_eigenpairs(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_k, &
            eigenvectors_k, site_orbital_counts, left_channels, right_channels, omega_static, chi0_options, chi0_static)
      end if
      if (circular_reverse) then
         chi0_options_reverse%q_direct = 0.0_rp
         green_options_reverse%q_direct = 0.0_rp
         if (canonical_chi0_backend == 'kspace_lehmann') then
            call green_source%initialize(eigenvalues_k, eigenvectors_k, eigenvalues_k, eigenvectors_k)
            call build_static_chi_ks_from_green_functions(green_source, reciprocal_obj%k_weights, site_orbital_counts, &
               left_channels_reverse, right_channels_reverse, green_options_reverse, chi0_static_reverse)
         else
            call build_static_chi_ks_from_eigenpairs(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, &
               site_orbital_counts, left_channels_reverse, right_channels_reverse, chi0_options_reverse, chi0_static_reverse)
         end if
      end if
      if (pair_backend) then
         call pair_operator_source%initialize(reciprocal_obj, signed_moments, [0.0_rp, 0.0_rp, 0.0_rp], &
            use_qplus=primary_minus_plus)
         call build_static_direct_xi_from_operator_source(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, &
            site_orbital_counts, left_channels, pair_operator_source, chi0_options, pair_xi_static)
         call evaluate_raw_xi_diagnostics(pair_xi_static%xi(:, :, 1), cmplx(signed_moments, 0.0_rp, rp), pair_goldstone, &
            response_basis='site', kernel_provenance='direct LMTO ham_only transverse pair-potential Xi')
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
         if (circular_reverse) then
            call pair_operator_source%initialize(reciprocal_obj, signed_moments, [0.0_rp, 0.0_rp, 0.0_rp], &
               use_qplus=.not. primary_minus_plus)
            call build_static_direct_xi_from_operator_source(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, &
               site_orbital_counts, left_channels_reverse, pair_operator_source, chi0_options_reverse, pair_xi_static_reverse)
            call evaluate_raw_xi_diagnostics(pair_xi_static_reverse%xi(:, :, 1), cmplx(signed_moments, 0.0_rp, rp), &
               pair_goldstone_reverse, response_basis='site', kernel_provenance='direct LMTO reverse circular pair-potential Xi')
         end if
         call pair_operator_source%clear()
      end if
      allocate(kernel(nresponse, nresponse))
      kernel = cmplx(0.0_rp, 0.0_rp, rp)
      if (canonical_chi0_backend == 'realspace_gf') then
         ! Bare native chi0 output does not invent a static kernel or
         ! Goldstone diagnostic; those require the future native static-limit
         ! provider documented in TDDFT-07.
      else if (is_longitudinal) then
         ! The coupled longitudinal kernel is derived from the converged
         ! ground-state XC route and the same projected charge metric used by
         ! the charge response.  No independently supplied static-field file
         ! and no rigid-rotation/Goldstone constraint enter this branch.
         if (.not. allocated(charge_obj%amad) .or. size(charge_obj%amad, 1) /= lattice_obj%nrec .or. &
            size(charge_obj%amad, 2) /= lattice_obj%nrec) then
            call g_logger%fatal('[calculation.post_processing_susceptibility]: longitudinal TDDFT requires an nrec-by-nrec projected Coulomb kernel.', &
               __FILE__, __LINE__)
         end if
         allocate(coulomb_site(lattice_obj%nrec, lattice_obj%nrec))
         coulomb_site = charge_obj%amad
         call self_obj%xc_response_provider%longitudinal_response_capability(full_response_supported, full_response_capability_reason)
         if (.not. full_response_supported) then
            call g_logger%fatal('[calculation.post_processing_susceptibility]: channel=''longitudinal'' is unavailable for XC route '''// &
               trim(self_obj%xc_response_provider%functional_label)//''': '//trim(full_response_capability_reason), __FILE__, __LINE__)
         end if
         call build_charge_longitudinal_kernel(self_obj%xc_response_provider, coulomb_site, kernel)
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
         goldstone_options%goldstone_policy = config%goldstone_policy
         goldstone_options%has_soc = has_soc
         goldstone_options%has_external_field = has_external_field
         goldstone_options%circular_channel = chi0_options%circular_channel
         call evaluate_goldstone(chi0_static%chi(:, :, 1), self_obj%xc_response_provider, goldstone_options, goldstone_result)
         if (circular_reverse) then
            goldstone_options%circular_channel = chi0_options_reverse%circular_channel
            call evaluate_goldstone(chi0_static_reverse%chi(:, :, 1), self_obj%xc_response_provider, goldstone_options, &
               goldstone_result_reverse)
         end if
         ! A finite-eta inverse chi_KS remains diagnostic only.  It is generally
         ! complex; using it as a frequency-independent adiabatic kernel forces
         ! an artificial Gamma singularity.  An explicitly selected TDDFT-02
         ! sum_rule/projected policy instead operates on the real static Ward
         ! problem and records its raw and corrected residuals.
         if (legacy_backend) then
         if (config%goldstone_policy == 'sum_rule' .and. goldstone_result%lounis%applied) then
               kernel = goldstone_result%lounis%kernel
            else if (config%goldstone_policy == 'projected' .and. goldstone_result%projection%applied) then
               kernel = goldstone_result%kernel_corrected
            else
               do isite = 1, lattice_obj%nrec
                  kernel(isite, isite) = goldstone_result%k_perp(isite)
               end do
            end if
            if (config%goldstone_policy == 'sum_rule' .and. .not. goldstone_result%lounis%applied) then
               call g_logger%fatal('[calculation.post_processing_susceptibility]: requested Lounis Ward reconstruction rejected: '// &
                  trim(goldstone_result%lounis%decision), __FILE__, __LINE__)
            end if
            if (config%goldstone_policy == 'projected' .and. .not. goldstone_result%projection%applied) then
               call g_logger%fatal('[calculation.post_processing_susceptibility]: requested Halle Goldstone projection rejected: '// &
                  trim(goldstone_result%projection%decision), __FILE__, __LINE__)
            end if
         end if
         if (rank == 0) then
            write(filename, '(a,"_goldstone.dat")') trim(config%output_prefix)
            call write_goldstone_diagnostics_text(trim(filename), goldstone_result)
            if (pair_backend) call append_pair_goldstone_diagnostics(trim(filename), goldstone_result%raw, pair_goldstone)
            if (config%goldstone_mode == 'correct') call append_goldstone_column_correction_text(trim(filename), pair_correction)
            call append_tddft_metadata(trim(filename), config, 0, reciprocal_obj%nk_mesh, [0.0_rp, 0.0_rp, 0.0_rp], &
               rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'goldstone_compare')
            if (circular_reverse) then
               write(filename, '(a,"_minus_plus_goldstone.dat")') trim(config%output_prefix)
               call write_goldstone_diagnostics_text(trim(filename), goldstone_result_reverse)
               if (pair_backend) call append_pair_goldstone_diagnostics(trim(filename), goldstone_result_reverse%raw, pair_goldstone_reverse)
               call append_tddft_metadata(trim(filename), config, 0, reciprocal_obj%nk_mesh, [0.0_rp, 0.0_rp, 0.0_rp], &
                  rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'goldstone_reverse_compare')
            end if
         end if
         if (config%goldstone_mode == 'correct' .and. .not. pair_correction%applied) then
            call g_logger%fatal('[calculation.post_processing_susceptibility]: requested Goldstone correction rejected: '// &
               trim(pair_correction%decision), __FILE__, __LINE__)
         end if
      end if

      has_gamma = any(maxval(abs(config%q_points), dim=1) <= 1.0e-12_rp)
      ! Selecting a pair-potential backend is itself a request to construct
      ! the raw Xi shadow data at every requested (q,omega), even if the user
      ! has disabled the optional text products.
      dyson_options%diagonalize_xi = config%output_modes
      dyson_options%diagonalize_loss = config%output_modes
      if (is_longitudinal) then
         dyson_options%circular_channel = 'charge_mz'
      else
         dyson_options%circular_channel = chi0_options%circular_channel
      end if
      if (config%output_modes) then
         allocate(all_xi(nresponse, nresponse, nw, nq), all_loss(nresponse, nresponse, nw, nq), all_trace_loss(nw, nq))
         all_xi = cmplx(0.0_rp, 0.0_rp, rp)
         all_loss = cmplx(0.0_rp, 0.0_rp, rp)
         all_trace_loss = 0.0_rp
         if (circular_reverse) then
            allocate(all_xi_reverse(nresponse, nresponse, nw, nq), all_loss_reverse(nresponse, nresponse, nw, nq), &
               all_trace_loss_reverse(nw, nq))
            all_xi_reverse = cmplx(0.0_rp, 0.0_rp, rp)
            all_loss_reverse = cmplx(0.0_rp, 0.0_rp, rp)
            all_trace_loss_reverse = 0.0_rp
         end if
         if (pair_backend .and. legacy_backend) then
            allocate(all_xi_pair(nresponse, nresponse, nw, nq), all_loss_pair(nresponse, nresponse, nw, nq), all_trace_loss_pair(nw, nq))
            all_xi_pair = cmplx(0.0_rp, 0.0_rp, rp)
            all_loss_pair = cmplx(0.0_rp, 0.0_rp, rp)
            all_trace_loss_pair = 0.0_rp
            if (circular_reverse) then
               allocate(all_xi_pair_reverse(nresponse, nresponse, nw, nq), all_loss_pair_reverse(nresponse, nresponse, nw, nq), &
                  all_trace_loss_pair_reverse(nw, nq))
               all_xi_pair_reverse = cmplx(0.0_rp, 0.0_rp, rp)
               all_loss_pair_reverse = cmplx(0.0_rp, 0.0_rp, rp)
               all_trace_loss_pair_reverse = 0.0_rp
            end if
         end if
      end if
#ifdef USE_MPI
      ! Root writes the static Goldstone record above; q owners may append the
      ! independently observed dynamic Gamma peaks below.
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
      if (canonical_chi0_backend == 'realspace_gf') then
         ! chi0(R,omega) is q independent.  Build it once for the complete
         ! requested path, then consume only its cheap susceptibility FT in
         ! the output loop below.
         select type (chi0_backend)
         type is (tddft_realspace_gf_backend)
            call chi0_backend%evaluate_grid(config%q_points, omega, realspace_batch)
            call reduce_native_realspace_batch(realspace_batch)
            if (circular_reverse) then
               call chi0_backend%initialize(realspace_source_reverse)
               call chi0_backend%evaluate_grid(config%q_points, omega, realspace_batch_reverse)
               call reduce_native_realspace_batch(realspace_batch_reverse)
               call chi0_backend%initialize(realspace_source)
            end if
         class default
            call g_logger%fatal('[calculation.post_processing_susceptibility]: native real-space batch adapter is absent.', &
               __FILE__, __LINE__)
         end select
      end if
      do iq = iq_start, iq_end
         is_gamma = maxval(abs(config%q_points(:, iq))) <= 1.0e-12_rp
         bare_gamma_peak = -1.0_rp; legacy_gamma_peak = -1.0_rp; pair_gamma_peak = -1.0_rp
         pair_corrected_gamma_peak = -1.0_rp
         bare_gamma_peak_reverse = -1.0_rp; legacy_gamma_peak_reverse = -1.0_rp; pair_gamma_peak_reverse = -1.0_rp
         pair_corrected_gamma_peak_reverse = -1.0_rp
         if (canonical_chi0_backend == 'realspace_gf') then
            kq_eigensolve_cpu_seconds = 0.0_rp
            chi0_result = realspace_batch%q_response(iq)
         else if (is_gamma) then
            ! At exact Gamma, k+q is the same folded endpoint as k.  Reuse
            ! the canonical solve instead of diagonalizing the identical
            ! Hamiltonians a second time; this also preserves a common basis
            ! in degenerate subspaces for the explicit transition vertices.
            eigenvalues_kq = eigenvalues_k
            eigenvectors_kq = eigenvectors_k
            kq_eigensolve_cpu_seconds = 0.0_rp
         else
            kq_workset = reciprocal_obj%k_workset%shifted(config%q_points(:, iq))
            call cpu_time(t_profile_start)
            call reciprocal_obj%calculate_eigenpairs_at_kpoints(kq_workset%points, eigenvalues_kq, eigenvectors_kq)
            call cpu_time(t_profile_stop)
            kq_eigensolve_cpu_seconds = t_profile_stop-t_profile_start
         end if
         ! The workset folds k+q into the Hamiltonian BZ, while the requested
         ! direct-basis path coordinate is retained in output provenance.
         chi0_options%q_direct = config%q_points(:, iq)
         green_options%q_direct = config%q_points(:, iq)
         if (canonical_chi0_backend == 'realspace_gf') then
            ! The complete q batch was evaluated above.  No G(k), k+q, or
            ! per-q endpoint solve is used on the native path.
         else if (is_full_response) then
            if (canonical_chi0_backend == 'kspace_lehmann') then
               call green_source%initialize(eigenvalues_k, eigenvectors_k, eigenvalues_kq, eigenvectors_kq)
               call build_four_component_chi_ks_from_green_functions(green_source, reciprocal_obj%k_weights, &
                  site_orbital_counts, omega, green_options, chi0_result)
            else
               call build_four_component_chi_ks(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_kq, &
                  eigenvectors_kq, site_orbital_counts, omega, chi0_options, chi0_result)
            end if
         else
            select type (chi0_backend)
            type is (tddft_eigenpair_backend)
               call chi0_backend%initialize(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, &
                  reshape(eigenvalues_kq, [size(eigenvalues_kq, 1), size(eigenvalues_kq, 2), 1]), &
                  reshape(eigenvectors_kq, [size(eigenvectors_kq, 1), size(eigenvectors_kq, 2), &
                  size(eigenvectors_kq, 3), 1]), reshape(config%q_points(:, iq), [3, 1]), site_orbital_counts, &
                  left_channels, right_channels, chi0_options, green_options)
               call chi0_backend%evaluate_frequency_batch(config%q_points(:, iq), omega, chi0_result, q_index=1)
            type is (tddft_kspace_lehmann_backend)
               call chi0_backend%initialize(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, &
                  reshape(eigenvalues_kq, [size(eigenvalues_kq, 1), size(eigenvalues_kq, 2), 1]), &
                  reshape(eigenvectors_kq, [size(eigenvectors_kq, 1), size(eigenvectors_kq, 2), &
                  size(eigenvectors_kq, 3), 1]), reshape(config%q_points(:, iq), [3, 1]), site_orbital_counts, &
                  left_channels, right_channels, chi0_options, green_options)
               call chi0_backend%evaluate_frequency_batch(config%q_points(:, iq), omega, chi0_result, q_index=1)
            class default
               call g_logger%fatal('[calculation.post_processing_susceptibility]: selected chi0 backend has no '// &
                  'initialized periodic endpoint adapter.', __FILE__, __LINE__)
            end select
         end if
         ! Keep the public file metadata consistent across the reference and
         ! resolvent adapters.  The endpoint coordinates are the requested
         ! direct-basis q; the folded k+q points remain an implementation
         ! detail of the reciprocal workset.
         chi0_result%metadata%q_direct = config%q_points(:, iq)
         chi0_result%metadata%response_projection = config%response_projection
         chi0_result%metadata%arbitrary_kq_cpu_seconds = kq_eigensolve_cpu_seconds
         response_eta = chi0_result%metadata%eta
         if (is_gamma) bare_gamma_peak = observed_loss_peak(omega, chi0_result%trace_spectrum)
         if (config%output_chi0 .or. config%output_stoner) then
            write(filename, '(a,"_q",i6.6,"_chi0.dat")') trim(config%output_prefix), iq
            call write_chi_ks_text(trim(filename), omega, chi0_result)
            call append_tddft_metadata(trim(filename), config, iq, reciprocal_obj%nk_mesh, config%q_points(:, iq), rank, &
               has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'shared_chi_ks')
         end if
         if (circular_reverse) then
            chi0_options_reverse%q_direct = config%q_points(:, iq)
            green_options_reverse%q_direct = config%q_points(:, iq)
            if (canonical_chi0_backend == 'realspace_gf') then
               chi0_result_reverse = realspace_batch_reverse%q_response(iq)
            else if (canonical_chi0_backend == 'kspace_lehmann') then
               call green_source%initialize(eigenvalues_k, eigenvectors_k, eigenvalues_kq, eigenvectors_kq)
               call build_chi_ks_from_green_functions(green_source, reciprocal_obj%k_weights, site_orbital_counts, &
                  left_channels_reverse, right_channels_reverse, omega, green_options_reverse, chi0_result_reverse)
            else
               select type (chi0_backend)
               type is (tddft_eigenpair_backend)
                  call chi0_backend%initialize(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, &
                     reshape(eigenvalues_kq, [size(eigenvalues_kq, 1), size(eigenvalues_kq, 2), 1]), &
                     reshape(eigenvectors_kq, [size(eigenvectors_kq, 1), size(eigenvectors_kq, 2), &
                     size(eigenvectors_kq, 3), 1]), reshape(config%q_points(:, iq), [3, 1]), site_orbital_counts, &
                     left_channels_reverse, right_channels_reverse, chi0_options_reverse, green_options_reverse)
                  call chi0_backend%evaluate_frequency_batch(config%q_points(:, iq), omega, chi0_result_reverse, q_index=1)
               type is (tddft_kspace_lehmann_backend)
                  call chi0_backend%initialize(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, &
                     reshape(eigenvalues_kq, [size(eigenvalues_kq, 1), size(eigenvalues_kq, 2), 1]), &
                     reshape(eigenvectors_kq, [size(eigenvectors_kq, 1), size(eigenvectors_kq, 2), &
                     size(eigenvectors_kq, 3), 1]), reshape(config%q_points(:, iq), [3, 1]), site_orbital_counts, &
                     left_channels_reverse, right_channels_reverse, chi0_options_reverse, green_options_reverse)
                  call chi0_backend%evaluate_frequency_batch(config%q_points(:, iq), omega, chi0_result_reverse, q_index=1)
               class default
                  call g_logger%fatal('[calculation.post_processing_susceptibility]: reverse circular backend adapter is absent.', &
                     __FILE__, __LINE__)
               end select
            end if
            chi0_result_reverse%metadata%q_direct = config%q_points(:, iq)
            chi0_result_reverse%metadata%response_projection = config%response_projection
            chi0_result_reverse%metadata%circular_channel = 'minus_plus'
            chi0_result_reverse%metadata%arbitrary_kq_cpu_seconds = kq_eigensolve_cpu_seconds
            if (is_gamma) bare_gamma_peak_reverse = observed_loss_peak(omega, chi0_result_reverse%trace_spectrum)
            if (config%output_chi0 .or. config%output_stoner) then
               write(filename, '(a,"_q",i6.6,"_minus_plus_chi0.dat")') trim(config%output_prefix), iq
               call write_chi_ks_text(trim(filename), omega, chi0_result_reverse)
               call append_tddft_metadata(trim(filename), config, iq, reciprocal_obj%nk_mesh, config%q_points(:, iq), rank, &
                  has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'shared_reverse_chi_ks')
            end if
         end if
         if (need_dyson) then
            if (legacy_backend) then
               call enhance_tddft_susceptibility(chi0_result%chi, kernel, response_eta, dyson_options, dyson_result)
               if (is_gamma) legacy_gamma_peak = observed_loss_peak(omega, dyson_result%trace_spectral_weight)
               if (config%output_xi .or. config%output_chi) then
                  if (pair_backend) then
                     write(filename, '(a,"_q",i6.6,"_legacy_dyson.dat")') trim(config%output_prefix), iq
                  else if (is_longitudinal) then
                     write(filename, '(a,"_q",i6.6,"_longitudinal_dyson.dat")') trim(config%output_prefix), iq
                  else
                     write(filename, '(a,"_q",i6.6,"_dyson.dat")') trim(config%output_prefix), iq
                  end if
                  call write_tddft_dyson_text(trim(filename), omega, dyson_result)
                  call append_tddft_metadata(trim(filename), config, iq, reciprocal_obj%nk_mesh, config%q_points(:, iq), &
                     rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'legacy_site_scalar_raw')
                  if (is_longitudinal) call append_longitudinal_response_metadata(trim(filename), kernel, &
                     self_obj%xc_response_provider%functional_label)
               end if
            end if
            if (pair_backend) then
               call pair_operator_source%initialize(reciprocal_obj, signed_moments, config%q_points(:, iq), &
                  use_qplus=primary_minus_plus)
               call build_direct_xi_from_operator_source(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, &
                  eigenvalues_kq, eigenvectors_kq, site_orbital_counts, left_channels, pair_operator_source, omega, &
                  chi0_options, pair_xi_result)
               call enhance_tddft_susceptibility_from_xi(chi0_result%chi, pair_xi_result%xi, response_eta, dyson_options, &
                  dyson_pair_result)
               raw_pair_minimum_spectral_weight = minval(dyson_pair_result%site_spectral_weight)
               raw_pair_spectral_weight_ok = spectral_weights_are_nonnegative(reshape( &
                  dyson_pair_result%site_spectral_weight, [size(dyson_pair_result%site_spectral_weight)]))
               if (is_gamma) pair_gamma_peak = observed_loss_peak(omega, dyson_pair_result%trace_spectral_weight)
               if (config%output_xi .or. config%output_chi) then
                  write(filename, '(a,"_q",i6.6,"_pair_dyson.dat")') trim(config%output_prefix), iq
                  call write_tddft_dyson_text(trim(filename), omega, dyson_pair_result)
                  call append_tddft_metadata(trim(filename), config, iq, reciprocal_obj%nk_mesh, config%q_points(:, iq), &
                     rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'pair_potential_raw')
                  call append_pair_spectral_weight_diagnostic(trim(filename), 'raw_pair', raw_pair_spectral_weight_ok, &
                     raw_pair_minimum_spectral_weight)
               end if
               if (config%goldstone_mode == 'correct') then
                  call pair_operator_source%initialize(reciprocal_obj, signed_moments, config%q_points(:, iq), &
                     pair_correction%scales, use_qplus=primary_minus_plus)
                  call build_direct_xi_from_operator_source(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, &
                     eigenvalues_kq, eigenvectors_kq, site_orbital_counts, left_channels, pair_operator_source, omega, &
                     chi0_options, pair_xi_corrected_result)
                  call enhance_tddft_susceptibility_from_xi(chi0_result%chi, pair_xi_corrected_result%xi, response_eta, &
                     dyson_options, dyson_pair_corrected_result)
                  if (is_gamma) then
                     pair_corrected_gamma_peak = observed_loss_peak(omega, dyson_pair_corrected_result%trace_spectral_weight)
                  end if
                  corrected_pair_minimum_spectral_weight = minval(dyson_pair_corrected_result%site_spectral_weight)
                  corrected_spectral_weight_ok = spectral_weights_are_nonnegative(reshape( &
                     dyson_pair_corrected_result%site_spectral_weight, [size(dyson_pair_corrected_result%site_spectral_weight)]))
                  correction_spectral_weight_ok = spectral_weight_correction_is_acceptable(reshape( &
                     dyson_pair_result%site_spectral_weight, [size(dyson_pair_result%site_spectral_weight)]), reshape( &
                     dyson_pair_corrected_result%site_spectral_weight, [size(dyson_pair_corrected_result%site_spectral_weight)]))
                  if (config%output_xi .or. config%output_chi) then
                     write(filename, '(a,"_q",i6.6,"_pair_corrected_dyson.dat")') trim(config%output_prefix), iq
                     call write_tddft_dyson_text(trim(filename), omega, dyson_pair_corrected_result)
                     call append_tddft_metadata(trim(filename), config, iq, reciprocal_obj%nk_mesh, config%q_points(:, iq), &
                     rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'pair_potential_corrected')
                     call append_goldstone_column_correction_text(trim(filename), pair_correction)
                     call append_pair_spectral_weight_diagnostic(trim(filename), 'corrected_pair', corrected_spectral_weight_ok, &
                        corrected_pair_minimum_spectral_weight)
                     call append_pair_correction_spectral_weight_diagnostic(trim(filename), correction_spectral_weight_ok)
                  end if
                  if (.not. correction_spectral_weight_ok) then
                     corrected_minimum_location = minloc(dyson_pair_corrected_result%site_spectral_weight)
                     write(spectral_weight_message, '(a,3(1x,es12.4),a,es12.4,a,i0,a,es12.4,a,es12.4,a,es12.4)') &
                        '[calculation.post_processing_susceptibility]: pair-potential negative spectral weight at q=', &
                        config%q_points(:, iq), ', omega=', omega(corrected_minimum_location(2)), ', site=', &
                        corrected_minimum_location(1), ', raw_min=', raw_pair_minimum_spectral_weight, ', corrected_min=', &
                        corrected_pair_minimum_spectral_weight, ', correction_max_delta=', pair_correction%maximum_change
                     call g_logger%fatal(trim(spectral_weight_message)//'. The static column correction created or worsened '// &
                        'negative spectral weight relative to the raw pair response and is rejected.', __FILE__, __LINE__)
                  end if
               end if
               call pair_operator_source%clear()
            end if
            if (config%output_modes) then
               if (legacy_backend) then
                  all_xi(:, :, :, iq) = dyson_result%xi
                  all_loss(:, :, :, iq) = dyson_result%loss
                  all_trace_loss(:, iq) = dyson_result%trace_spectral_weight
               else
                  all_xi(:, :, :, iq) = dyson_pair_result%xi
                  all_loss(:, :, :, iq) = dyson_pair_result%loss
                  all_trace_loss(:, iq) = dyson_pair_result%trace_spectral_weight
               end if
               if (pair_backend .and. legacy_backend) then
                  all_xi_pair(:, :, :, iq) = dyson_pair_result%xi
                  all_loss_pair(:, :, :, iq) = dyson_pair_result%loss
                  all_trace_loss_pair(:, iq) = dyson_pair_result%trace_spectral_weight
               end if
            end if
            if (circular_reverse) then
               response_eta = chi0_result_reverse%metadata%eta
               dyson_options_reverse = dyson_options
               dyson_options_reverse%circular_channel = 'minus_plus'
               if (legacy_backend) then
                  call enhance_tddft_susceptibility(chi0_result_reverse%chi, kernel, response_eta, dyson_options_reverse, &
                     dyson_result_reverse)
                  if (is_gamma) legacy_gamma_peak_reverse = observed_loss_peak(omega, &
                     dyson_result_reverse%trace_spectral_weight)
                  if (config%output_xi .or. config%output_chi) then
                     write(filename, '(a,"_q",i6.6,"_minus_plus_legacy_dyson.dat")') trim(config%output_prefix), iq
                     call write_tddft_dyson_text(trim(filename), omega, dyson_result_reverse)
                     call append_tddft_metadata(trim(filename), config, iq, reciprocal_obj%nk_mesh, config%q_points(:, iq), &
                        rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'legacy_reverse_site_scalar_raw')
                  end if
               end if
               if (pair_backend) then
                  call pair_operator_source%initialize(reciprocal_obj, signed_moments, config%q_points(:, iq), &
                     use_qplus=.not. primary_minus_plus)
                  call build_direct_xi_from_operator_source(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, &
                     eigenvalues_kq, eigenvectors_kq, site_orbital_counts, left_channels_reverse, pair_operator_source, omega, &
                     chi0_options_reverse, pair_xi_result_reverse)
                  call enhance_tddft_susceptibility_from_xi(chi0_result_reverse%chi, pair_xi_result_reverse%xi, response_eta, &
                     dyson_options_reverse, dyson_pair_result_reverse)
                  if (is_gamma) pair_gamma_peak_reverse = observed_loss_peak(omega, &
                     dyson_pair_result_reverse%trace_spectral_weight)
                  if (config%output_xi .or. config%output_chi) then
                     write(filename, '(a,"_q",i6.6,"_minus_plus_pair_dyson.dat")') trim(config%output_prefix), iq
                     call write_tddft_dyson_text(trim(filename), omega, dyson_pair_result_reverse)
                     call append_tddft_metadata(trim(filename), config, iq, reciprocal_obj%nk_mesh, config%q_points(:, iq), &
                        rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'pair_reverse_potential_raw')
                  end if
                  if (config%goldstone_mode == 'correct') then
                     call pair_operator_source%initialize(reciprocal_obj, signed_moments, config%q_points(:, iq), &
                        pair_correction%scales, use_qplus=.not. primary_minus_plus)
                     call build_direct_xi_from_operator_source(reciprocal_obj%k_weights, eigenvalues_k, eigenvectors_k, &
                        eigenvalues_kq, eigenvectors_kq, site_orbital_counts, left_channels_reverse, pair_operator_source, omega, &
                        chi0_options_reverse, pair_xi_corrected_result_reverse)
                     call enhance_tddft_susceptibility_from_xi(chi0_result_reverse%chi, &
                        pair_xi_corrected_result_reverse%xi, response_eta, dyson_options_reverse, &
                        dyson_pair_corrected_result_reverse)
                     if (config%output_xi .or. config%output_chi) then
                        write(filename, '(a,"_q",i6.6,"_minus_plus_pair_corrected_dyson.dat")') trim(config%output_prefix), iq
                        call write_tddft_dyson_text(trim(filename), omega, dyson_pair_corrected_result_reverse)
                        call append_tddft_metadata(trim(filename), config, iq, reciprocal_obj%nk_mesh, config%q_points(:, iq), &
                           rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'pair_reverse_potential_corrected')
                        call append_goldstone_column_correction_text(trim(filename), pair_correction)
                     end if
                  end if
                  call pair_operator_source%clear()
               end if
               if (config%output_modes) then
                  if (legacy_backend) then
                     all_xi_reverse(:, :, :, iq) = dyson_result_reverse%xi
                     all_loss_reverse(:, :, :, iq) = dyson_result_reverse%loss
                     all_trace_loss_reverse(:, iq) = dyson_result_reverse%trace_spectral_weight
                  else
                     all_xi_reverse(:, :, :, iq) = dyson_pair_result_reverse%xi
                     all_loss_reverse(:, :, :, iq) = dyson_pair_result_reverse%loss
                     all_trace_loss_reverse(:, iq) = dyson_pair_result_reverse%trace_spectral_weight
                  end if
                  if (pair_backend .and. legacy_backend) then
                     all_xi_pair_reverse(:, :, :, iq) = dyson_pair_result_reverse%xi
                     all_loss_pair_reverse(:, :, :, iq) = dyson_pair_result_reverse%loss
                     all_trace_loss_pair_reverse(:, iq) = dyson_pair_result_reverse%trace_spectral_weight
                  end if
               end if
            end if
         end if
         if (is_gamma .and. .not. is_longitudinal .and. .not. is_full_response .and. &
             canonical_chi0_backend /= 'realspace_gf') then
            write(filename, '(a,"_goldstone.dat")') trim(config%output_prefix)
            call append_dynamic_gamma_peaks(trim(filename), bare_gamma_peak, legacy_gamma_peak, pair_gamma_peak, &
               pair_corrected_gamma_peak, bare_gamma_peak_reverse, legacy_gamma_peak_reverse, pair_gamma_peak_reverse, &
               pair_corrected_gamma_peak_reverse)
         end if
         if (allocated(eigenvalues_kq)) deallocate(eigenvalues_kq)
         if (allocated(eigenvectors_kq)) deallocate(eigenvectors_kq)
      end do

      if (config%output_modes) then
#ifdef USE_MPI
         call MPI_ALLREDUCE(MPI_IN_PLACE, all_xi, size(all_xi), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, all_loss, size(all_loss), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, all_trace_loss, size(all_trace_loss), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         if (circular_reverse) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, all_xi_reverse, size(all_xi_reverse), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, all_loss_reverse, size(all_loss_reverse), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, all_trace_loss_reverse, size(all_trace_loss_reverse), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         end if
         if (pair_backend .and. legacy_backend) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, all_xi_pair, size(all_xi_pair), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, all_loss_pair, size(all_loss_pair), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, all_trace_loss_pair, size(all_trace_loss_pair), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            if (circular_reverse) then
               call MPI_ALLREDUCE(MPI_IN_PLACE, all_xi_pair_reverse, size(all_xi_pair_reverse), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
               call MPI_ALLREDUCE(MPI_IN_PLACE, all_loss_pair_reverse, size(all_loss_pair_reverse), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
               call MPI_ALLREDUCE(MPI_IN_PLACE, all_trace_loss_pair_reverse, size(all_trace_loss_pair_reverse), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            end if
         end if
#endif
         if (rank == 0) then
            call analyze_tddft_modes(omega, all_xi, all_trace_loss, response_eta, mode_options, mode_result, all_loss)
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
               call analyze_tddft_modes(omega, all_xi_pair, all_trace_loss_pair, response_eta, mode_options, mode_result, all_loss_pair)
               write(filename, '(a,"_pair_modes.dat")') trim(config%output_prefix)
               call write_tddft_modes_text(trim(filename), omega, response_eta, mode_result)
               call append_tddft_metadata(trim(filename), config, 0, reciprocal_obj%nk_mesh, [0.0_rp, 0.0_rp, 0.0_rp], &
                  rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'pair_potential_raw')
            end if
            if (circular_reverse) then
               call analyze_tddft_modes(omega, all_xi_reverse, all_trace_loss_reverse, response_eta, mode_options, mode_result, &
                  all_loss_reverse)
               if (legacy_backend) then
                  write(filename, '(a,"_minus_plus_legacy_modes.dat")') trim(config%output_prefix)
               else if (pair_backend) then
                  write(filename, '(a,"_minus_plus_pair_modes.dat")') trim(config%output_prefix)
               else
                  write(filename, '(a,"_minus_plus_modes.dat")') trim(config%output_prefix)
               end if
               call write_tddft_modes_text(trim(filename), omega, response_eta, mode_result)
               call append_tddft_metadata(trim(filename), config, 0, reciprocal_obj%nk_mesh, [0.0_rp, 0.0_rp, 0.0_rp], &
                  rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'reverse_circular_raw')
               if (pair_backend .and. legacy_backend) then
                  call analyze_tddft_modes(omega, all_xi_pair_reverse, all_trace_loss_pair_reverse, response_eta, mode_options, &
                     mode_result, all_loss_pair_reverse)
                  write(filename, '(a,"_minus_plus_pair_modes.dat")') trim(config%output_prefix)
                  call write_tddft_modes_text(trim(filename), omega, response_eta, mode_result)
                  call append_tddft_metadata(trim(filename), config, 0, reciprocal_obj%nk_mesh, [0.0_rp, 0.0_rp, 0.0_rp], &
                     rank, has_soc, has_external_field, trim(reciprocal_obj%reciprocal_mode), 'reverse_pair_potential_raw')
               end if
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
            write(unit, '(a)') '# TDDFT coupled charge-longitudinal manifest; site-major (charge,m_z) response'
         else
            write(unit, '(a)') '# TDDFT transverse chi_KS manifest; one self-describing output file per q point'
         end if
         if (.not. is_longitudinal .and. .not. is_full_response) then
            write(unit, '(a,a)') '# primary_circular_channel = ', trim(chi0_options%circular_channel)
            if (circular_reverse) write(unit, '(a)') '# reverse_circular_channel = minus_plus; files carry the _minus_plus_ suffix'
         end if
         if (pair_backend .and. legacy_backend .and. config%goldstone_mode == 'correct') then
            write(unit, '(a)') '# q_index q1 q2 q3 chi0_file legacy_raw_dyson_file pair_raw_dyson_file pair_corrected_dyson_file'
         else if (pair_backend .and. legacy_backend) then
            write(unit, '(a)') '# q_index q1 q2 q3 chi0_file legacy_raw_dyson_file pair_raw_dyson_file'
         else if (pair_backend .and. config%goldstone_mode == 'correct') then
            write(unit, '(a)') '# q_index q1 q2 q3 chi0_file pair_raw_dyson_file pair_corrected_dyson_file'
         else if (pair_backend) then
            write(unit, '(a)') '# q_index q1 q2 q3 chi0_file pair_raw_dyson_file'
         else if (is_longitudinal) then
            write(unit, '(a)') '# q_index q1 q2 q3 chi0_file longitudinal_dyson_file'
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
            else if (is_longitudinal) then
               write(legacy_filename, '(a,"_q",i6.6,"_longitudinal_dyson.dat")') trim(config%output_prefix), iq
               write(unit, '(i0,3(1x,es24.16),2(1x,a))') iq, config%q_points(:, iq), trim(chi0_filename), trim(legacy_filename)
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

   subroutine append_dynamic_gamma_peaks(filename, bare_peak, legacy_peak, pair_peak, pair_corrected_peak, &
      bare_peak_reverse, legacy_peak_reverse, pair_peak_reverse, pair_corrected_peak_reverse)
      character(len=*), intent(in) :: filename
      real(rp), intent(in) :: bare_peak, legacy_peak, pair_peak, pair_corrected_peak
      real(rp), intent(in), optional :: bare_peak_reverse, legacy_peak_reverse, pair_peak_reverse, pair_corrected_peak_reverse
      integer :: unit, ios

      open(newunit=unit, file=filename, status='old', position='append', action='write', iostat=ios)
      if (ios /= 0) call g_logger%fatal('[calculation.append_dynamic_gamma_peaks]: cannot append Gamma peaks', __FILE__, __LINE__)
      write(unit, '(a,es24.16)') '# dynamic_bare_gamma_loss_peak_Ry = ', bare_peak
      if (legacy_peak >= 0.0_rp) write(unit, '(a,es24.16)') '# dynamic_legacy_raw_gamma_loss_peak_Ry = ', legacy_peak
      if (pair_peak >= 0.0_rp) write(unit, '(a,es24.16)') '# dynamic_pair_raw_gamma_loss_peak_Ry = ', pair_peak
      if (pair_corrected_peak >= 0.0_rp) then
         write(unit, '(a,es24.16)') '# dynamic_pair_corrected_gamma_loss_peak_Ry = ', pair_corrected_peak
      end if
      if (present(bare_peak_reverse)) then
         if (bare_peak_reverse >= 0.0_rp) write(unit, '(a,es24.16)') '# dynamic_minus_plus_bare_gamma_loss_peak_Ry = ', bare_peak_reverse
         if (legacy_peak_reverse >= 0.0_rp) write(unit, '(a,es24.16)') '# dynamic_minus_plus_legacy_raw_gamma_loss_peak_Ry = ', legacy_peak_reverse
         if (pair_peak_reverse >= 0.0_rp) write(unit, '(a,es24.16)') '# dynamic_minus_plus_pair_raw_gamma_loss_peak_Ry = ', pair_peak_reverse
         if (pair_corrected_peak_reverse >= 0.0_rp) write(unit, '(a,es24.16)') '# dynamic_minus_plus_pair_corrected_gamma_loss_peak_Ry = ', &
            pair_corrected_peak_reverse
      end if
      write(unit, '(a)') '# dynamic Gamma peaks are observed loss-grid maxima; raw and corrected records remain distinct'
      close(unit)
   end subroutine append_dynamic_gamma_peaks

   subroutine append_pair_spectral_weight_diagnostic(filename, label, acceptable, minimum_weight)
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: label
      logical, intent(in) :: acceptable
      real(rp), intent(in) :: minimum_weight
      integer :: unit, ios

      open(newunit=unit, file=filename, status='old', position='append', action='write', iostat=ios)
      if (ios /= 0) call g_logger%fatal('[calculation.append_pair_spectral_weight_diagnostic]: cannot append diagnostic', &
         __FILE__, __LINE__)
      write(unit, '(a,l1)') '# '//trim(label)//'_spectral_weight_nonnegative = ', acceptable
      write(unit, '(a,es24.16)') '# '//trim(label)//'_minimum_site_spectral_weight = ', minimum_weight
      close(unit)
   end subroutine append_pair_spectral_weight_diagnostic

   subroutine append_pair_correction_spectral_weight_diagnostic(filename, acceptable)
      character(len=*), intent(in) :: filename
      logical, intent(in) :: acceptable
      integer :: unit, ios

      open(newunit=unit, file=filename, status='old', position='append', action='write', iostat=ios)
      if (ios /= 0) call g_logger%fatal('[calculation.append_pair_correction_spectral_weight_diagnostic]: cannot append diagnostic', &
         __FILE__, __LINE__)
      write(unit, '(a,l1)') '# pair_correction_preserves_raw_spectral_weight = ', acceptable
      close(unit)
   end subroutine append_pair_correction_spectral_weight_diagnostic

   subroutine initialize_lmto_pair_operator_source(this, reciprocal_obj, signed_moments, q_point, column_scales, use_qplus)
      class(lmto_pair_operator_tile_source), intent(inout) :: this
      type(reciprocal), target, intent(inout) :: reciprocal_obj
      real(rp), intent(in) :: signed_moments(:)
      real(rp), intent(in) :: q_point(3)
      real(rp), intent(in), optional :: column_scales(:)
      logical, intent(in), optional :: use_qplus
      integer :: nmat, nright

      call this%clear()
      if (.not. associated(reciprocal_obj%lattice) .or. .not. allocated(reciprocal_obj%k_points)) then
         call g_logger%fatal('[calculation.lmto_pair_operator_source]: reciprocal mesh is unavailable.', __FILE__, __LINE__)
      end if
      nright = reciprocal_obj%lattice%nrec
      if (size(signed_moments) /= nright) then
         call g_logger%fatal('[calculation.lmto_pair_operator_source]: signed-moment shape is incompatible.', __FILE__, __LINE__)
      end if
      if (present(column_scales)) then
         if (size(column_scales) /= nright) then
            call g_logger%fatal('[calculation.lmto_pair_operator_source]: correction-scale shape is incompatible.', __FILE__, __LINE__)
         end if
         allocate(this%column_scales(nright))
         this%column_scales = column_scales
      end if
      nmat = 2*norb*nright
      this%reciprocal_obj => reciprocal_obj
      allocate(this%signed_moments(nright))
      this%signed_moments = signed_moments
      this%q_point = q_point
      if (present(use_qplus)) this%use_qplus = use_qplus
      allocate(this%qminus(nmat, nmat), this%qplus(nmat, nmat))
   end subroutine initialize_lmto_pair_operator_source

   integer function lmto_pair_operator_channel_dimension(this) result(nright)
      class(lmto_pair_operator_tile_source), intent(in) :: this
      if (.not. associated(this%reciprocal_obj)) error stop 'LMTO pair operator source is not configured'
      nright = this%reciprocal_obj%lattice%nrec
   end function lmto_pair_operator_channel_dimension

   subroutine fill_lmto_pair_operator_tile(this, ik, operator_tile)
      class(lmto_pair_operator_tile_source), intent(inout) :: this
      integer, intent(in) :: ik
      complex(rp), intent(out) :: operator_tile(:, :, :)
      integer :: isite, nmat, nright
      logical :: supported
      character(len=160) :: reason

      if (.not. associated(this%reciprocal_obj) .or. .not. allocated(this%signed_moments)) then
         error stop 'LMTO pair operator source is not configured'
      end if
      nright = this%reciprocal_obj%lattice%nrec
      nmat = 2*norb*nright
      if (ik < 1 .or. ik > size(this%reciprocal_obj%k_points, 2) .or. &
          any(shape(operator_tile) /= [nmat, nmat, nright])) then
         error stop 'LMTO pair operator source tile shape is incompatible'
      end if
      do isite = 1, nright
         call this%reciprocal_obj%build_lmto_pair_potential_at_kpoint(isite, this%reciprocal_obj%k_points(:, ik), &
            this%signed_moments(isite), this%qminus, this%qplus, supported, reason, this%q_point)
         if (.not. supported) then
            call g_logger%fatal('[calculation.lmto_pair_operator_source]: pair-potential construction rejected: '// &
               trim(reason), __FILE__, __LINE__)
         end if
         if (this%use_qplus) then
            operator_tile(:, :, isite) = this%qplus
         else
            operator_tile(:, :, isite) = this%qminus
         end if
         if (allocated(this%column_scales)) operator_tile(:, :, isite) = this%column_scales(isite)*operator_tile(:, :, isite)
      end do
      this%fetch_count = this%fetch_count + 1
   end subroutine fill_lmto_pair_operator_tile

   subroutine clear_lmto_pair_operator_source(this)
      class(lmto_pair_operator_tile_source), intent(inout) :: this
      if (allocated(this%qminus)) deallocate(this%qminus)
      if (allocated(this%qplus)) deallocate(this%qplus)
      if (allocated(this%signed_moments)) deallocate(this%signed_moments)
      if (allocated(this%column_scales)) deallocate(this%column_scales)
      nullify(this%reciprocal_obj)
      this%q_point = 0.0_rp
      this%use_qplus = .false.
      this%fetch_count = 0
   end subroutine clear_lmto_pair_operator_source

   subroutine lmto_pair_operator_source_destructor(this)
      type(lmto_pair_operator_tile_source), intent(inout) :: this
      call this%clear()
   end subroutine lmto_pair_operator_source_destructor

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
      write(unit, '(a,es24.16)') '# pair_potential_raw_ward_residual = ', diagnostics%ward%ward_residual
      write(unit, '(a,es24.16)') '# pair_potential_raw_dm_residual = ', diagnostics%ward%dm_residual
      write(unit, '(a,a)') '# pair_potential_raw_response_basis = ', trim(diagnostics%ward%response_basis)
      write(unit, '(a,a)') '# pair_potential_raw_kernel_provenance = ', trim(diagnostics%ward%kernel_provenance)
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
      if (trim(config%chi0_backend) == 'green') then
         write(unit, '(a)') '# chi0_backend_capability = eigenpair_resolvent reference; native RS Green-function provider unavailable'
      else
         write(unit, '(a)') '# chi0_backend_capability = eigenpair Lehmann response'
      end if
      if (trim(config%channel) == 'full') then
         write(unit, '(a)') '# full_response_capability = validated selected-XC ALSDA derivatives'
      else
         write(unit, '(a)') '# full_response_capability = not selected; production full response is capability-gated'
      end if
      if (trim(config%channel) == 'longitudinal') then
         write(unit, '(a)') '# longitudinal_capability = coupled charge-m_z response with ground-state ALSDA plus Hartree'
         write(unit, '(a)') '# longitudinal_response_basis = site-major (charge,m_z); Hartree is charge-charge only'
         write(unit, '(a)') '# longitudinal_goldstone_constraint = none'
         write(unit, '(a)') '# longitudinal_llb_status = susceptibility output only; no damping or LLB parameter is inferred'
      else
         write(unit, '(a)') '# longitudinal_capability = not selected; LLB readiness is not claimed'
      end if
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
      write(unit, '(a,es24.16)') '# response_fermi_level_Ry = ', config%fermi_level
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
      write(unit, '(a,a)') '# goldstone_policy = ', trim(config%goldstone_policy)
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
   !>        /= ’auto’ path of post_processing_frozen_magnon).
   !> @details For mode=’mft’ the routine converges the reference potential
   !>          once, at q_ss_list(:,1). It then reuses that fixed potential for
   !>          one force-theorem band-energy probe (frozen_magnon_probe_energy)
   !>          at every other q. A k-space run picks its q-mesh policy
   !>          (little_group, little_group_common, or full_bz) from
   !>          reciprocal%q_symmetry_policy. For mode=’scf’ the routine fully
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
   !>               q in mode=’scf’.
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
      real(rp), allocatable :: etot_q(:), eband_q(:), eband_gauge_q(:), mtot_q(:, :), omega_q(:)
      real(rp), allocatable :: fermi_q(:), electron_count_q(:), electron_error_q(:), target_electrons_q(:), weight_sum_q(:)
      integer, allocatable :: nk_total_q(:), nk_mesh_q(:, :)
      real(rp), allocatable :: constraint_metric_q(:), constraint_penalty_q(:), constraint_coupling_q(:)
      real(rp), allocatable :: potential_drift_q(:), potential_transient_drift_q(:), field_magnitude_q(:, :)
      integer, allocatable :: constraint_iteration_q(:)
      logical, allocatable :: constraint_converged_q(:)
      real(rp) :: sin2theta, etot_ref, theta_probe, ordinary_checksum
      real(rp) :: delta_raw, delta_gauge, delta_final, e_q0, e_ref0
      logical :: use_kspace, constraints_enabled, bare_mft, corrected_mft
      character(len=200) :: fmt_str
      integer :: iq, i, iter, newunit, diagnostic_unit

      bare_mft = fm_obj%mode == 'mft'
      corrected_mft = fm_obj%mode == 'mft_constrained'
      allocate (etot_q(fm_obj%n_q), eband_q(fm_obj%n_q), mtot_q(lattice_obj%nrec, fm_obj%n_q), &
                omega_q(fm_obj%n_q), constraint_metric_q(fm_obj%n_q), constraint_penalty_q(fm_obj%n_q), &
                constraint_coupling_q(fm_obj%n_q), potential_drift_q(fm_obj%n_q), &
                potential_transient_drift_q(fm_obj%n_q), &
                field_magnitude_q(lattice_obj%nrec, fm_obj%n_q), constraint_iteration_q(fm_obj%n_q), &
                constraint_converged_q(fm_obj%n_q), fermi_q(fm_obj%n_q), electron_count_q(fm_obj%n_q), &
                electron_error_q(fm_obj%n_q), target_electrons_q(fm_obj%n_q), weight_sum_q(fm_obj%n_q), nk_total_q(fm_obj%n_q), &
                nk_mesh_q(3, fm_obj%n_q))
      if (bare_mft) allocate (eband_gauge_q(fm_obj%n_q))
      constraint_metric_q = 0.0_rp
      constraint_penalty_q = 0.0_rp
      constraint_coupling_q = 0.0_rp
      potential_drift_q = 0.0_rp
      potential_transient_drift_q = 0.0_rp
      field_magnitude_q = 0.0_rp
      constraint_iteration_q = 0
      constraint_converged_q = .true.
      fermi_q = 0.0_rp
      electron_count_q = 0.0_rp
      electron_error_q = 0.0_rp
      target_electrons_q = 0.0_rp
      weight_sum_q = 0.0_rp
      nk_total_q = 0
      nk_mesh_q = 0

      ! Reference point (row 1): converge the flat-spiral cone potential once. Its
      ! magnitudes/moments define the normalization.  Corrected MFT first obtains
      ! this ordinary reference with the auxiliary constraint operator disabled.
      constraints_enabled = lattice_obj%control%constraints_enable
      if ((bare_mft .or. corrected_mft) .and. constraints_enabled) then
         ! Disable before constructing self_obj: its constructor initializes the
         ! constraint state, so delaying this switch would reject an ordinary
         ! reference whose input potential has no explicit constraint target.
         lattice_obj%control%constraints_enable = .false.
      end if
      self_obj = self(bands_obj, mix_obj)
      use_kspace = self_obj%use_kspace
      ! GBT is a magnetic representation shared by both solvers.
      hamiltonian_obj%magnetic_representation = gbt_single_q
      hamiltonian_obj%q_ss(:) = q_ss_cart(:, 1)
      if ((bare_mft .or. corrected_mft) .and. constraints_enabled) then
         ! Bare MFT and the ordinary reference of corrected MFT are explicitly
         ! constraint-free.  The corrected mode creates a fresh constrained self
         ! object after this run, retaining the converged ordinary potential.
         self_obj%control%constraints_enable = .false.
      else if (fm_obj%mode == 'scf' .and. constraints_enabled) then
         call set_constrained_spiral_targets(self_obj)
      end if
      call self_obj%run()
      etot_ref = sum(lattice_obj%symbolic_atoms(:)%potential%etot)
      ordinary_checksum = self_obj%potential_checksum()

      if (corrected_mft .and. constraints_enabled) then
         self_obj%control%constraints_enable = .true.
         self_obj = self(bands_obj, mix_obj)
         if (fm_obj%constraint_tolerance > 0.0_rp) then
            self_obj%constraint_tolerance = fm_obj%constraint_tolerance
         end if
      end if
      do i = 1, lattice_obj%nrec
         mtot_q(i, :) = lattice_obj%symbolic_atoms(lattice_obj%nbulk + i)%potential%mtot
      end do

      if (bare_mft .or. (corrected_mft .and. .not. constraints_enabled)) then
         ! Bare MFT and the corrected zero-field limit rebuild ONLY the Hamiltonian
         ! per q and take a band-energy probe at the fixed ordinary potential.
         if (use_kspace) then
            recip_obj = reciprocal(hamiltonian_obj)
            ! WP8: the mesh must never be built once from row 1’s q_ss and
            ! reused blindly for every other q in the sweep -- each policy
            ! below either shares one mesh proven valid for the whole sweep,
            ! or rebuilds per q inside the loop.
            select case (trim(recip_obj%q_symmetry_policy))
            case ('little_group_common')
               ! One mesh, reduced by the little group common to every q in
               ! the sweep (not just row 1’s), valid for every probe below.
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
               ! specific q -- never reuse another q’s little-group mesh.
               call recip_obj%ensure_kpoint_mesh(recip_obj%nk_mesh, sum(abs(recip_obj%k_offset)) > 1.0e-12_rp)
            end if

            ! At finite q the GBT phase is a pure gauge when the cone is
            ! collapsed onto the collinear reference axis.  In the continuum
            ! E(q,theta=0) = E(0,theta), but a finite k mesh does not preserve
            ! that identity when q is not a mesh translation.  Subtracting
            ! E(0,theta) therefore leaves a q-only integration artifact which
            ! is amplified by 1/sin^2(theta).  On the reciprocal route,
            ! evaluate the exact same-q, zero-cone reference on the same
            ! potential and same k-point set; this is a physical gauge
            ! identity, not an empirical theta correction.  The real-space
            ! path retains its established q=0 subtraction for compatibility.
            if (bare_mft .and. use_kspace) then
               theta_probe = hamiltonian_obj%theta_ss
               hamiltonian_obj%theta_ss = 0.0_rp
               eband_gauge_q(iq) = frozen_magnon_probe_energy(bands_obj, recip_obj, energy_obj, use_kspace)
               hamiltonian_obj%theta_ss = theta_probe
            end if
            eband_q(iq) = frozen_magnon_probe_energy(bands_obj, recip_obj, energy_obj, use_kspace)
            if (use_kspace) then
               fermi_q(iq) = recip_obj%fermi_level
               electron_count_q(iq) = recip_obj%canonical_electron_count
               target_electrons_q(iq) = recip_obj%total_electrons
               electron_error_q(iq) = recip_obj%canonical_electron_count - target_electrons_q(iq)
               weight_sum_q(iq) = recip_obj%canonical_weight_sum
               nk_total_q(iq) = recip_obj%nk_total
               nk_mesh_q(:, iq) = recip_obj%nk_mesh
            end if
            if (bare_mft .and. .not. use_kspace) eband_gauge_q(iq) = eband_q(1)
            etot_q(iq) = etot_ref   ! potential frozen; total energy not re-evaluated
            mtot_q(:, iq) = mtot_q(:, 1)
         end do
      else if (corrected_mft) then
         if (.not. constraints_enabled) then
            call g_logger%fatal('[calculation.post_processing_frozen_magnon]: corrected MFT requires constraints_enable = .true. or an explicit zero-field comparison path.', &
                               __FILE__, __LINE__)
         end if
         ! Corrected constrained MFT (Jacobsson convention): the ordinary
         ! potential is fixed, but the constraining field is converged separately
         ! at every q.  Only the final one-shot occupied-eigenvalue sum is used.
         if (use_kspace) then
            recip_obj = reciprocal(hamiltonian_obj)
            select case (trim(recip_obj%q_symmetry_policy))
            case ('little_group_common')
               call recip_obj%ensure_kpoint_mesh(recip_obj%nk_mesh, &
                                                 sum(abs(recip_obj%k_offset)) > 1.0e-12_rp, &
                                                 q_list_cart=q_ss_cart)
            case ('little_group')
               ! Built for the current q inside the loop.
            case default
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
         do iq = 1, fm_obj%n_q
            hamiltonian_obj%q_ss(:) = q_ss_cart(:, iq)
            if (use_kspace .and. trim(recip_obj%q_symmetry_policy) == 'little_group') then
               call recip_obj%ensure_kpoint_mesh(recip_obj%nk_mesh, sum(abs(recip_obj%k_offset)) > 1.0e-12_rp)
            end if
            call set_constrained_spiral_targets(self_obj)
            call self_obj%reset_constraint_for_fixed_potential(fm_obj%constraint_start_from_zero)
            do iter = 1, fm_obj%constraint_max_iterations
               call self_obj%run_fixed_potential_constraint_step(iter)
               if (self_obj%constraint_converged) exit
            end do
            if (.not. self_obj%constraint_converged) then
               call g_logger%fatal('[calculation.post_processing_frozen_magnon]: constraining field did not converge at q index '// &
                                  int2str(iq)//' within constraint_max_iterations.', __FILE__, __LINE__)
            end if
            eband_q(iq) = frozen_magnon_probe_energy(bands_obj, recip_obj, energy_obj, use_kspace)
            if (use_kspace) then
               fermi_q(iq) = recip_obj%fermi_level
               electron_count_q(iq) = recip_obj%canonical_electron_count
               target_electrons_q(iq) = recip_obj%total_electrons
               electron_error_q(iq) = recip_obj%canonical_electron_count - target_electrons_q(iq)
               weight_sum_q(iq) = recip_obj%canonical_weight_sum
               nk_total_q(iq) = recip_obj%nk_total
               nk_mesh_q(:, iq) = recip_obj%nk_mesh
            end if
            etot_q(iq) = etot_ref
            mtot_q(:, iq) = mtot_q(:, 1)
            constraint_metric_q(iq) = self_obj%constraint_metric
            constraint_penalty_q(iq) = self_obj%constraint_penalty_energy
            constraint_coupling_q(iq) = self_obj%constraint_field_coupling_energy
            potential_drift_q(iq) = self_obj%fixed_potential_max_drift
            potential_transient_drift_q(iq) = self_obj%fixed_potential_transient_drift
            constraint_iteration_q(iq) = self_obj%constraint_iteration
            constraint_converged_q(iq) = self_obj%constraint_converged
            do i = 1, lattice_obj%nrec
               field_magnitude_q(i, iq) = norm2(lattice_obj%symbolic_atoms(lattice_obj%nbulk + i)%mag_cfield)
            end do
         end do
      else
         ! scf: fully self-consistent spiral at each q (each q independently relaxed).
         do iq = 1, fm_obj%n_q
            hamiltonian_obj%q_ss(:) = q_ss_cart(:, iq)
            self_obj = self(bands_obj, mix_obj)
            if (self_obj%control%constraints_enable) call set_constrained_spiral_targets(self_obj)
            call self_obj%run()
            etot_q(iq) = sum(lattice_obj%symbolic_atoms(:)%potential%etot)
            if (.not. self_obj%use_kspace) call bands_obj%calculate_band_energy()
            eband_q(iq) = bands_obj%eband
            do i = 1, lattice_obj%nrec
               mtot_q(i, iq) = lattice_obj%symbolic_atoms(lattice_obj%nbulk + i)%potential%mtot
            end do
            if (self_obj%control%constraints_enable) then
               constraint_metric_q(iq) = self_obj%constraint_metric
               constraint_penalty_q(iq) = self_obj%constraint_penalty_energy
               constraint_coupling_q(iq) = self_obj%constraint_field_coupling_energy
               constraint_iteration_q(iq) = self_obj%constraint_iteration
               constraint_converged_q(iq) = self_obj%constraint_converged
               do i = 1, lattice_obj%nrec
                  field_magnitude_q(i, iq) = norm2(lattice_obj%symbolic_atoms(lattice_obj%nbulk + i)%mag_cfield)
               end do
            end if
         end do
      end if

      sin2theta = sin(hamiltonian_obj%theta_ss)**2
      if (sin2theta < 1.0e-8_rp) then
         call g_logger%fatal('[calculation.post_processing_frozen_magnon]: '// &
                             'theta_ss must be nonzero (a finite cone angle) to define omega(q)', __FILE__, __LINE__)
      end if
      if (abs(sum(mtot_q(:, 1))) < 1.0e-12_rp) then
         call g_logger%fatal('[calculation.post_processing_frozen_magnon]: total reference moment is zero; omega normalization is undefined.', &
                            __FILE__, __LINE__)
      end if
      do iq = 1, fm_obj%n_q
         if (bare_mft) then
            omega_q(iq) = 4.0_rp*(eband_q(iq) - eband_gauge_q(iq))/(sum(mtot_q(:, 1))*sin2theta)
         else if (corrected_mft) then
            omega_q(iq) = 4.0_rp*(eband_q(iq) - eband_q(1))/(sum(mtot_q(:, 1))*sin2theta)
         else
            omega_q(iq) = 4.0_rp*(etot_q(iq) - etot_q(1))/(sum(mtot_q(:, 1))*sin2theta)
         end if
      end do

      open (newunit=newunit, file=trim(fm_obj%output_file), status='replace', action='write')
      write (newunit, '(A)') '# Frozen-magnon sweep (calculation%post_processing = "frozen_magnon")'
      write (newunit, '(A)') '# q_ss units: Cartesian 2*pi/alat (same convention as &hamiltonian q_ss); row 1 is the reference point'
      write (newunit, '(A,A)') '# q_file coordinates: ', trim(fm_obj%q_coordinates)
      write (newunit, '(A)') '# omega uses eband in modes="mft"/"mft_constrained" and etot in mode="scf"'
      if (bare_mft .and. use_kspace) then
         write (newunit, '(A)') '# MFT k-space omega uses E(q,theta)-E(q,theta=0); E(q,theta=0) is the same-q GBT gauge reference'
      end if
      if (corrected_mft) then
         write (newunit, '(A)') '# corrected MFT: ordinary potential frozen; q-specific constraint field converged before one-shot band energy'
         write (newunit, '(A)') '# detailed field, residual, checksum, and energy-component records are in frozen_magnon_constrained_diagnostics.dat'
      end if
      write (newunit, '(A,I0)') '# Number of sublattices (nrec): ', lattice_obj%nrec
      write (newunit, '(A)') '# Format: q1 q2 q3 etot eband mtot_1 .. mtot_nrec omega'
      write (fmt_str, '(A,I0,A)') '(3F12.6,1X,2(ES16.8E3,1X),', lattice_obj%nrec, '(ES16.8E3,1X),ES16.8E3)'
      do iq = 1, fm_obj%n_q
         write (newunit, fmt_str) q_ss_cart(:, iq), etot_q(iq), eband_q(iq), mtot_q(:, iq), omega_q(iq)
      end do
      close (newunit)

      if (bare_mft .and. use_kspace) then
         ! Keep the two subtraction choices visible in a separate diagnostic
         ! file.  This makes the gauge-invariant production observable
         ! auditable without changing the long-standing frozen_magnon.dat
         ! column contract used by the quick fixtures.
         open (newunit=diagnostic_unit, file='frozen_magnon_diagnostics.dat', status='replace', action='write')
         write (diagnostic_unit, '(A)') '# Frozen-magnon MFT energy diagnostics'
         write (diagnostic_unit, '(A)') '# Format: q1 q2 q3 E(q,theta) E(0,theta) DeltaE_raw E(q,0) DeltaE_gauge sin2theta Mtot omega'
         do iq = 1, fm_obj%n_q
            delta_raw = eband_q(iq) - eband_q(1)
            delta_gauge = eband_q(iq) - eband_gauge_q(iq)
            write (diagnostic_unit, '(3F12.6,1X,8(ES20.12E3,1X))') q_ss_cart(:, iq), eband_q(iq), eband_q(1), &
               delta_raw, eband_gauge_q(iq), delta_gauge, sin2theta, sum(mtot_q(:, 1)), omega_q(iq)
         end do
         close (diagnostic_unit)
      end if

      if (corrected_mft) then
         open (newunit=diagnostic_unit, file='frozen_magnon_constrained_diagnostics.dat', status='replace', action='write')
         write (diagnostic_unit, '(A)') '# Frozen-magnon corrected constrained-MFT diagnostics'
         write (diagnostic_unit, '(A)') '# mode = mft_constrained'
         write (diagnostic_unit, '(A,3(1X,ES24.16))') '# reference_q =', q_ss_cart(:, 1)
         write (diagnostic_unit, '(A)') '# ordinary_potential_frozen = T'
         write (diagnostic_unit, '(A,ES24.16)') '# ordinary_potential_checksum_L1 = ', ordinary_checksum
         write (diagnostic_unit, '(A)') '# constraint_field_converged = per_q_column'
         write (diagnostic_unit, '(A,L1)') '# constraint_start_from_zero = ', fm_obj%constraint_start_from_zero
         write (diagnostic_unit, '(A)') '# gauge_reference_used = F'
         write (diagnostic_unit, '(A)') '# field_coupling_and_controller_penalty_are_diagnostics_only = T'
         write (diagnostic_unit, '(A)') '# energy_bookkeeping = occupied_eigenvalue_sum(fixed_ordinary_potential + converged_q_field)'
         write (diagnostic_unit, '(A)') '# cone_normalization = 4*DeltaE/(sum(M_reference)*sin(theta_ss)^2)'
         write (diagnostic_unit, '(A)') '# Format: q1 q2 q3 E_phys_reference E_band_raw E_band_reference DeltaE_raw DeltaE_final gauge_reference residual penalty field_coupling potential_checksum potential_max_drift potential_transient_drift controller_iteration field_converged sin2theta M_reference omega Bmag_1 .. Bmag_nrec'
         do iq = 1, fm_obj%n_q
            delta_raw = eband_q(iq) - eband_q(1)
            delta_final = delta_raw
            delta_gauge = 0.0_rp
            write (diagnostic_unit, *) q_ss_cart(:, iq), etot_q(iq), eband_q(iq), eband_q(1), delta_raw, &
               delta_final, delta_gauge, constraint_metric_q(iq), constraint_penalty_q(iq), constraint_coupling_q(iq), &
               ordinary_checksum, potential_drift_q(iq), potential_transient_drift_q(iq), constraint_iteration_q(iq), &
               merge(1, 0, constraint_converged_q(iq)), sin2theta, sum(mtot_q(:, 1)), omega_q(iq), &
               field_magnitude_q(:, iq)
         end do
         close (diagnostic_unit)
      end if

      if (use_kspace) then
         ! WP09 machine-readable contract.  Keep this file independent of the
         ! legacy frozen_magnon.dat and of the mode-specific constraint audit:
         ! every reciprocal acoustic run emits the raw q/theta energies plus
         ! the canonical occupation and mesh metadata needed to distinguish a
         ! physical cone response from finite-k quadrature noise.  For
         ! corrected MFT, the zero-cone columns are zero and gauge_available is
         ! false because the q-specific constraint field is not a pure-gauge
         ! reference at theta=0.
         open (newunit=diagnostic_unit, file='frozen_magnon_harmonic_diagnostics.dat', status='replace', action='write')
         write (diagnostic_unit, '(A)') '# GBT WP09 harmonic cone-angle and k-grid diagnostics'
         write (diagnostic_unit, '(A)') '# schema = gbt_wp09_harmonic_v1'
         write (diagnostic_unit, '(A,A)') '# mode = ', trim(fm_obj%mode)
         write (diagnostic_unit, '(A,A)') '# q_coordinates = ', trim(fm_obj%q_coordinates)
         write (diagnostic_unit, '(A)') '# q_units = Cartesian 2*pi/alat, matching hamiltonian%q_ss'
         write (diagnostic_unit, '(A,ES24.16)') '# alat_angstrom = ', lattice_obj%alat
         write (diagnostic_unit, '(A)') '# gauge_available = true only for bare reciprocal MFT same-q theta=0 probes'
         write (diagnostic_unit, '(A)') '# energy_bookkeeping = raw occupied-eigenvalue energies; no controller penalty or field coupling added'
         write (diagnostic_unit, '(A)') '# columns = q1 q2 q3 theta_deg E_q_theta E_q0 E_qref_theta E_qref0 DeltaE_raw DeltaE_gauge DeltaE_pure sin2theta Mtot omega fermi_level electron_count electron_error target_electrons weight_sum nk1 nk2 nk3 nk_total gauge_available'
         do iq = 1, fm_obj%n_q
            if (bare_mft .and. use_kspace) then
               e_q0 = eband_gauge_q(iq)
               e_ref0 = eband_gauge_q(1)
               delta_gauge = eband_q(iq) - eband_gauge_q(iq)
               delta_final = eband_gauge_q(iq) - eband_gauge_q(1)
            else
               e_q0 = 0.0_rp
               e_ref0 = 0.0_rp
               delta_gauge = 0.0_rp
               delta_final = 0.0_rp
            end if
            delta_raw = eband_q(iq) - eband_q(1)
            write (diagnostic_unit, *) q_ss_cart(:, iq), hamiltonian_obj%theta_ss*180.0_rp/acos(-1.0_rp), &
               eband_q(iq), e_q0, eband_q(1), e_ref0, delta_raw, delta_gauge, delta_final, &
               sin2theta, sum(mtot_q(:, 1)), omega_q(iq), fermi_q(iq), electron_count_q(iq), electron_error_q(iq), &
               target_electrons_q(iq), weight_sum_q(iq), nk_mesh_q(:, iq), nk_total_q(iq), &
               merge(1, 0, bare_mft .and. use_kspace)
         end do
         close (diagnostic_unit)
      end if
   end subroutine post_processing_frozen_magnon_acoustic

   !> @brief Set the controller targets to the collinear axis of the GBT frame.
   !> @details GBT keeps potential%mom as a collinear reference marker and puts
   !>          the physical cone in the bond rotation.  The constraint controller
   !>          therefore targets the local +/-z axis, not the lab-frame cone
   !>          vector.  This is also the axis consumed by the shared density
   !>          contract; using the cone vector here would ask a local onsite field
   !>          to rotate a moment that GBT has already rotated in the hopping.
   subroutine set_constrained_spiral_targets(self_obj)
      type(self), intent(inout) :: self_obj
      integer :: ia, itype
      real(rp) :: axis_sign

      if (.not. allocated(self_obj%constraint_reference)) then
         call g_logger%fatal('[calculation.set_constrained_spiral_targets]: constraint state is not initialized.', &
                            __FILE__, __LINE__)
      end if
      do ia = 1, self_obj%lattice%nrec
         itype = self_obj%lattice%nbulk + ia
         axis_sign = 1.0_rp
         if (self_obj%symbolic_atom(itype)%potential%mom(3) < 0.0_rp) axis_sign = -1.0_rp
         self_obj%constraint_reference(:, ia) = [0.0_rp, 0.0_rp, axis_sign]
         self_obj%symbolic_atom(itype)%constraint_target = self_obj%constraint_reference(:, ia)
      end do
   end subroutine set_constrained_spiral_targets

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
   !>        (imposed through the GBT Hamiltonian’s q_ss + per-sublattice theta/phi). The
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
      real(rp), allocatable :: eval(:), mtot_ref(:), ref_theta(:), ref_phi(:), single_energy(:), energy_gauge_q(:)
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
            ! sweep (not just row 1’s), valid for every probe below.
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

      allocate (omega_mat(nactive, nactive), eigvec(nactive, nactive), eval(nactive), single_energy(nactive), &
                energy_gauge_q(fm_obj%n_q))
      do iq = 1, fm_obj%n_q
         hamiltonian_obj%q_ss(:) = q_ss_cart(:, iq)
         if (use_kspace .and. trim(recip_obj%q_symmetry_policy) == 'little_group') then
            ! Rebuild (or, via the cache key, reuse if unchanged) for this
            ! specific q -- never reuse another q’s little-group mesh (WP8).
            call recip_obj%ensure_kpoint_mesh(recip_obj%nk_mesh, sum(abs(recip_obj%k_offset)) > 1.0e-12_rp)
         end if
         omega_mat(:, :) = cmplx(0.0_rp, 0.0_rp, kind=rp)

         ! The collinear theta=0 state is a q-dependent pure gauge in the
         ! continuum.  A finite k mesh can give it a small q-only integration
         ! offset, which otherwise contaminates every harmonic coefficient.
         ! Use the same-q collinear probe as the force-theorem reference on the
         ! reciprocal route.  This is an invariant-preserving subtraction, not
         ! a fitted diagonal shift; retain the historical RS reference energy.
         if (use_kspace) then
            call set_probe_angles(hamiltonian_obj, ref_theta, ref_phi, active_type, theta_probe, active_rec(1:0))
            energy_gauge_q(iq) = frozen_magnon_probe_energy(bands_obj, recip_obj, energy_obj, use_kspace)
         else
            energy_gauge_q(iq) = energy_ref
         end if

         ! Diagonal: single-sublattice tilt. d^2E/dtheta_i^2 / M_i, in omega = dE/dm_z units.
         do iact = 1, nactive
            call set_probe_angles(hamiltonian_obj, ref_theta, ref_phi, active_type, theta_probe, [iact])
            single_energy(iact) = frozen_magnon_probe_energy(bands_obj, recip_obj, energy_obj, use_kspace)
            omega_mat(iact, iact) = cmplx((single_energy(iact) - energy_gauge_q(iq))/(active_moment(iact)*dmz_fac), &
                                          0.0_rp, kind=rp)
         end do

         ! Off-diagonal: two-sublattice tilt at the NATURAL GBT spiral phase (the q-phase
         ! is carried by the Hamiltonian’s bond rotation, so no azimuthal-phase probe is
         ! needed). Cross second derivative d^2E/dtheta_i dtheta_j, real symmetric:
         !   omega_ij = (E_ij - E_i - E_j + E_ref) / [2 (1-cos theta) sqrt(M_i M_j)].
         do iact = 1, nactive - 1
            do jact = iact + 1, nactive
               call set_probe_angles(hamiltonian_obj, ref_theta, ref_phi, active_type, theta_probe, [iact, jact])
               e_pair = frozen_magnon_probe_energy(bands_obj, recip_obj, energy_obj, use_kspace)
               off = (e_pair - single_energy(iact) - single_energy(jact) + energy_gauge_q(iq)) / &
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
            ! MFT must not depend on DOS window/grid/projections.  Seed the
            ! reciprocal occupation evaluator from the converged SCF EF and
            ! let the reciprocal namelist decide whether EF is re-solved from
            ! the target electron count.  The previous unconditional
            ! find_fermi=.true. made auto_find_fermi=.false. ineffective on
            ! this force-theorem path and silently changed the occupation
            ! contract at every q.
            recip_obj%fermi_level = energy_obj%fermi
            e_band = recip_obj%calculate_canonical_band_energy()
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
      this%do_inertia = .false.
   end subroutine restore_to_default

   !> Reduce locally owned R-block contributions while retaining the complete
   !> q batch on every rank.  The native source is therefore parallel over its
   !> expensive real-space/energy work, but q Fourier phases are evaluated once
   !> per rank for the whole requested path and q-labelled files remain owned by
   !> the normal output range.
   subroutine reduce_native_realspace_batch(batch)
      type(tddft_chi0_batch_result), intent(inout) :: batch
#ifdef USE_MPI
      integer :: local_counts(2), global_counts(2)
      real(rp) :: local_tail, global_tail
      integer :: iq

      do iq = 1, size(batch%q_response)
         call MPI_ALLREDUCE(MPI_IN_PLACE, batch%q_response(iq)%chi, size(batch%q_response(iq)%chi), &
            MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%re_chi = real(batch%q_response(iq)%chi, rp)
         batch%q_response(iq)%im_chi = aimag(batch%q_response(iq)%chi)
         call MPI_ALLREDUCE(MPI_IN_PLACE, batch%q_response(iq)%site_diagonal_spectrum, &
            size(batch%q_response(iq)%site_diagonal_spectrum), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, batch%q_response(iq)%stoner_spectral_map, &
            size(batch%q_response(iq)%stoner_spectral_map), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, batch%q_response(iq)%trace_spectrum, &
            size(batch%q_response(iq)%trace_spectrum), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         local_counts = [batch%q_response(iq)%metadata%real_space_points, &
            batch%q_response(iq)%metadata%real_space_omitted_points]
         call MPI_ALLREDUCE(local_counts, global_counts, 2, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_points = global_counts(1)
         batch%q_response(iq)%metadata%real_space_omitted_points = global_counts(2)
         local_tail = batch%q_response(iq)%metadata%real_space_tail_norm
         call MPI_ALLREDUCE(local_tail, global_tail, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_tail_norm = global_tail
      end do
      if (size(batch%q_response) > 0) batch%metadata = batch%q_response(1)%metadata
#else
      ! Serial builds already hold the complete real-space source locally.
      continue
#endif
   end subroutine reduce_native_realspace_batch

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
          .and. post_processing /= 'fermi_surface' &
          .and. post_processing /= 'kspace_green' &
          .and. post_processing /= 'frozen_magnon' &
          .and. post_processing /= 'susceptibility') then
         call g_logger%fatal('[calculation.check_post_processing]: '// &
                             "calculation%post_processing must be one of: ''none'', ''paoflow2rs'', ''exchange'', ''exchange_p2rs''," // &
                             " 'conductivity', 'conductivity_p2rs', 'orbital_modern', 'band_structure', 'bsf', 'density_of_states'," // &
                             " 'fermi_surface', 'kspace_green', 'frozen_magnon', 'susceptibility'", __FILE__, __LINE__)
      end if
   end subroutine check_post_processing

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Check availability for the intersite Green’s-function route (B2.5)
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
