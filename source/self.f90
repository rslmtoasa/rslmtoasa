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
   use symbolic_atom_mod, only: symbolic_atom, save_state_scf
   use potential_mod, only: potential
   use logger_mod, only: g_logger
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   use string_mod, only: real2str, int2str, fmt, sl, log2str
   use control_mod
   use lattice_mod
   use charge_mod
   use xc_mod
   use xc_radial_mod, only: radgra
   use radial_ground_state_mod, only: radial_ground_state
   use recursion_mod
   use density_of_states_mod
   use green_mod
   use bands_mod
   use energy_mod
   use hamiltonian_mod
   use magnetic_representation_mod, only: gbt_single_q
   use spin_density_mod, only: sd_constrained_spiral, sd_relaxed_reference
   use mix_mod
   use reciprocal_mod
   use electrostatics_multipole_mod, only: compute_dipole_moments
   use math_mod
   use precision_mod, only: rp
   use timer_mod, only: g_timer
   use scf_benchmark_profile_mod, only: g_scf_benchmark_profile
   use scf_guard_mod, only: guard_q0_capacity
   use namelist_generator_mod, only: namelist_generator
      use cfd
      use Parameters
#ifdef USE_MPI
   use mpi
#endif
   use basis_mod, only: nb, norb, spin_off, lmax_basis
   implicit none

   private

   ! Historical ATOMSC control predicates are kept as small pure procedures
   ! so the iteration policy can be tested independently of the physical
   ! residual observables. Production defaults remain historical unless the
   ! explicitly named test-only ATOMSC controls are supplied.
   public :: atomsc_solver_tolerance, atomsc_mixing_beta, atomsc_should_stop
   public :: fixed_u_r05_targets, atomsc_inner_controls
   public :: vxc0sp_pbe_origin_derivatives
   public :: vxc0sp
   !> Deterministic audit hook used by LR-BASIS-00 tests.  It calls the live
   !> RSEQSR/PHDFSR machinery; it does not provide an alternative radial solver.
   public :: legacy_radial_fixture
   !> Test-only seam that invokes the unchanged production NEWRHO routine.
   public :: legacy_newrho_fixture
   !> Private accepted-state replay seam: run the real bands%calculate_pl and
   !> unchanged NEWRHO on fresh scratch objects supplied from a hash-bound
   !> capture.  This is an audit/test hook, not a production calculation path.
   public :: legacy_accepted_pl_newrho_fixture
   public :: legacy_accepted_m03_newrho_fixture
   public :: legacy_accepted_m03_components_fixture
   !> Public audit seam for the unchanged scalar-relativistic GINTSR metric.
   public :: legacy_gintsr
   ! The atomic/XC radial representation has two local spin eigenchannels.
   ! This is independent of control%nsp, which selects the global Hamiltonian
   ! mode (and is not a channel-count field).
   integer, parameter, public :: RADIAL_XC_SPIN_CHANNELS = 2
   integer, parameter :: fixed_u_m0_a11_nshift = 8
   integer, parameter :: fixed_u_m0_a2_nshift = 14
   real(rp), parameter :: fixed_u_m0_a2_v(fixed_u_m0_a2_nshift) = &
      [-2.0e-2_rp, -1.0e-2_rp, -5.0e-3_rp, -2.5e-3_rp, &
        2.5e-3_rp,  5.0e-3_rp,  1.0e-2_rp,  2.0e-2_rp, &
       -4.0e-2_rp, -8.0e-2_rp, -1.2e-1_rp, &
        4.0e-2_rp,  8.0e-2_rp,  1.2e-1_rp]

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

      !> Fixed reference directions used by the constraining-field update.
      real(rp), dimension(:, :), allocatable :: constraint_reference
      !> Physical DFT total energy of the current density, excluding the
      !> auxiliary constraining-field operator (Ry).
      real(rp) :: physical_total_energy
      !> Numerical penalty/controller merit value returned by the most recent
      !> constraint update (Ry).  This is not part of physical_total_energy.
      real(rp) :: constraint_penalty_energy
      !> Expectation value <V_constraint> = sum_i B_i^con dot m_i (Ry).  This
      !> is a field/eigenvalue bookkeeping diagnostic, not a DFT energy term.
      real(rp) :: constraint_field_coupling_energy
      !> Deprecated compatibility alias for constraint_penalty_energy.
      real(rp) :: constraint_energy
      !> RMS angular constraint metric from the most recent update (radians).
      real(rp) :: constraint_metric
      !> Declared convergence tolerance for the angular constraint metric.
      real(rp) :: constraint_tolerance
      !> Whether the most recent constraint update satisfies its tolerance.
      logical :: constraint_converged
      !> Controller-update number represented by the stored diagnostics.
      integer :: constraint_iteration
      !> True after the optional machine-readable file has received its header.
      logical :: constraint_diagnostics_header_written
      !> True after the optional WP12 SCF trace has received its header.
      logical :: gbt_scf_diagnostics_header_written
      !> User-controlled temporary symmetry-breaking field state.  This is
      !> deliberately separate from the persistent constraint machinery.
      logical :: magnetic_seed_active
      integer :: magnetic_seed_steps
      real(rp) :: magnetic_seed_field
      logical :: magnetic_seed_enable
      !> Measured moment entering the current SCF iteration.  `mix%mag_old`
      !> intentionally stores the frame marker, not this physical vector.
      real(rp), dimension(:, :), allocatable :: gbt_scf_in_moment
      !> L1 checksum of the ordinary potential used by a fixed-potential pass.
      real(rp) :: fixed_potential_checksum = 0.0_rp
      !> Maximum ordinary-potential drift after restoration in the most recent
      !> fixed-potential constraint step; this is the frozen-potential audit value.
      real(rp) :: fixed_potential_max_drift = 0.0_rp
      !> Diagnostic magnitude of the temporary DOS-side update before the fixed
      !> potential snapshot is restored.
      real(rp) :: fixed_potential_transient_drift = 0.0_rp

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

      ! Optional R05 fixed-U tolerances. Nonpositive values select backend defaults.
      real(rp) :: r05_p_tolerance
      real(rp) :: r05_v_tolerance

      ! Test-only atomsc settings; all three must be supplied together.
      real(rp) :: atomsc_test_tol
      real(rp) :: atomsc_test_tolrsq
      integer :: atomsc_test_inner_cap

      ! Constrains variables

      !> Freezes magnetic moment direction during nonlinear calculation. Default: false.
      !>
      !> Freezes magnetic moment direction during nonlinear calculation. Do not set true in collinear calculation.
      !>
      !> Default: false.
      logical :: freeze

      !> Freezes the strength of the spin-orbit coupling interaction to the assigned value in the potential file
      !>
      !> Freezes the strength of the spin-orbit coupling interaction to the assigned value in the potential file  
      !>
      !> Default: false.
      logical :: fix_soc


      !> Freezes the strength of the spin-orbit coupling interaction to the assigned value in the potential file
      !>
      !> Freezes the strength of the spin-orbit coupling interaction to the assigned value in the potential file  
      !>
      !> Default: false.
      real(rp) :: soc_scale
       

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
      !> Signed normalized radial partial-wave amplitude per (mesh, l+1, spin),
      !> exported from NEWRHO for the B6 multipole (l=1 dipole) radial matrix
      !> elements. Runs parallel to fun2 (which stores the |amplitude|^2 density).
      real(rp), dimension(:, :, :), allocatable :: phi_amp
      ! Exact core pieces accumulated by RHOCOR for the LR-02N valence/core
      ! bookkeeping.  The first array retains the accepted scalar-relativistic
      ! core density; the second is its large-component-only diagnostic.
      real(rp), dimension(:, :), allocatable :: core_rho_weighted, core_pauli_weighted
      ! Optional l-resolved copy of the same frozen-core accumulation for the
      ! DRESP-09W core provenance audit.  It is diagnostic-only and does not
      ! feed the production density.
      real(rp), dimension(:, :, :), allocatable :: core_rho_weighted_l, core_pauli_weighted_l

      ! TODO
      real(rp), dimension(:), allocatable :: bxc

      ! TODO max ws
      ! User in lmtst: MAKE POTENTIAL PARAMETERS TO SAME DNU VALUES
      real(rp) :: ws_max
      !> Logical variable to check if the calculation is converged.
      logical :: converged
      !> Last SCF iteration reached by the accepted state.
      integer :: converged_iteration

      !> Logical variable to control if initial
      ! potential parameters are calculated from moments
      logical :: cold
      !> Experimental opt-in k-space reciprocal diagnostics path for SCF loop.
      logical :: use_kspace
      !> Cached reciprocal helper for k-space SCF branch (reused across iterations).
      type(reciprocal), allocatable :: reciprocal_scf_cache
      !> Opt-in source-connected fixed-Liechtenstein trace. Empty unless the
      !> RSLMTO_FIXED_U_TRACE environment variable names an output CSV file.
      character(len=1024), private :: fixed_u_trace_filename = ''
      character(len=128), private :: fixed_u_trace_run_id = 'run'
      real(rp), private :: fixed_u_trace_lambda = 1.0_rp
      logical, private :: fixed_u_trace_enabled = .false.
      logical, private :: fixed_u_trace_capture_only = .false.
      logical, private :: fixed_u_trace_moments_only = .false.
      logical, private :: fixed_u_trace_ql_only = .false.
      logical, private :: fixed_u_trace_a2_scan_enabled = .false.
      integer, private :: fixed_u_trace_iteration = 0
      integer, private :: fixed_u_trace_selected_iteration = 0
      integer, private :: fixed_u_trace_current_site = -1
      logical, private :: fixed_u_trace_capture_radial_state = .false.
      ! M0-only immutable snapshots from one selected outer iteration. These
      ! are diagnostic inputs for the later frozen radial replay and are never
      ! copied back into the live SCF atom or Hamiltonian.
      logical, private :: fixed_u_m0_snapshot_valid = .false.
      logical, private :: fixed_u_m0_p_available = .false.
      integer, private :: fixed_u_m0_snapshot_iteration = -1
      real(rp), allocatable, private :: fixed_u_m0_q0(:, :), fixed_u_m0_center(:, :), fixed_u_m0_q2(:, :)
      real(rp), allocatable, private :: fixed_u_m0_local_center(:, :), fixed_u_m0_local_q2(:, :), fixed_u_m0_trvp(:, :)
      real(rp), allocatable, private :: fixed_u_m0_vout_center(:, :), fixed_u_m0_vout_q2(:, :)
      real(rp), allocatable, private :: fixed_u_m0_pt_center(:, :), fixed_u_m0_pt_q2(:, :)
      real(rp), allocatable, private :: fixed_u_m0_t_vout(:, :), fixed_u_m0_t_pt(:, :)
      real(rp), allocatable, private :: fixed_u_m0_ql(:, :, :, :), fixed_u_m0_pl_diag(:, :, :, :)
      complex(rp), allocatable, private :: fixed_u_m0_vin(:, :, :)
      ! Full complex fixed-U convergence state.  These gates are inactive for
      ! no-U, SC-U, +V, and non-Liechtenstein workflows.
      logical, private :: fixed_u_residual_active = .false.
      logical, private :: fixed_u_residual_converged = .false.
      logical, private :: fixed_u_residual_physical = .false.
      real(rp), private :: fixed_u_rp_abs = huge(1.0_rp)
      real(rp), private :: fixed_u_rp_rel = huge(1.0_rp)
      real(rp), private :: fixed_u_rv_abs = huge(1.0_rp)
      real(rp), private :: fixed_u_rv_rel = huge(1.0_rp)
      real(rp), private :: fixed_u_antiherm_p = huge(1.0_rp)
      real(rp), private :: fixed_u_antiherm_v = huge(1.0_rp)
      integer, private :: fixed_u_worst_site = 0
      integer, private :: fixed_u_worst_l = -1

   contains
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: print_state_formatted
      procedure :: print_state
      procedure :: print_state_full
      procedure :: run
      procedure :: run_fixed_potential_constraint_step
      procedure :: reset_constraint_for_fixed_potential
      procedure :: potential_checksum
      procedure :: finalize_kspace_scf_state
      procedure :: write_kspace_scf_state_artifact
      procedure :: report
      procedure :: quantify_pauli_projection
      procedure :: lmtst
      procedure :: is_converged
      procedure, private :: initialize_constraint_state
      procedure, private :: initialize_gbt_scf_diagnostics
      procedure, private :: set_magnetic_seed_state
      procedure, private :: apply_constraints
      procedure, private :: mix_magnetic_state
      procedure, private :: write_gbt_scf_diagnostics
      procedure, private :: run_scf
      procedure, private :: initialize_fixed_u_trace
      procedure, private :: reset_fixed_u_trace
      procedure, private :: write_fixed_u_trace
      procedure, private :: write_fixed_u_radial_trace
      procedure, private :: write_fixed_u_kspace_trace
      procedure, private :: write_fixed_u_rs_bonds
      procedure, private :: write_fixed_u_csv_row
      procedure, private :: write_fixed_u_total_energy
      procedure, private :: write_hubbard_moment_energy_diagnostics
      procedure, private :: validate_hubbard_mode_provenance
      procedure, private :: is_frozen_hubbard_host
      procedure, private :: all_frozen_hubbard_hosts_local
      procedure, private :: is_eligible_local_hubbard
      procedure, private :: is_eligible_local_moments
      procedure, private :: fixed_hubbard_energy_pout
      procedure, private :: write_fixed_u_radial_replay_input
      procedure, private :: write_fixed_u_radial_replay_output
      procedure, private :: verify_fixed_u_radial_replay
      procedure, private :: snapshot_fixed_u_m0_moments
      procedure, private :: replay_fixed_u_m0_moments
      procedure, private :: update_fixed_u_residuals
      procedure, private :: write_fixed_u_rs_green
      procedure, private :: atomsc
      procedure, private :: ftype
      procedure, private :: newrho
      procedure, private :: rhocor
      procedure, private :: vxc0sp
      procedure, private :: racsi
      !> @note B7.2: visibility widened from `private` to public. This is
      !> an ACCESS SPECIFIER CHANGE ONLY -- not one character of POTPAR’s
      !> body, nor of the radial solver it drives (RSEQSR/PHDFSR/RSQSR1/
      !> RSQSR2/FCTP/FCTP0/GINTSR), has been touched, so behaviour is
      !> bit-identical.
      !>
      !> Reason: B7 §1.6 requires the vacuum-lead parameter generator
      !> (source/vacuum_lead.f90) to produce empty-lattice parameters with
      !> *the code’s own radial solver driven at V(r) = const*, explicitly
      !> NOT by hand-coded spherical-Bessel expressions -- so that the
      !> representation, the gamma/o convention, the E_nu convention and
      !> the normalization come out automatically consistent with every
      !> other parameter set in the code. POTPAR is the only entry point
      !> that produces that set, and the entire solver beneath it is
      !> module-private, so the alternative was duplicating ~600 lines of
      !> physics-dense legacy code into a parallel copy that could silently
      !> drift out of convention. That is precisely the failure mode §1.6
      !> is written to prevent.
      !>
      !> POTPAR references no component of `self` whatsoever -- it reads
      !> only its `atom`, `V` and `ROFI` arguments -- so exposing it
      !> implies no SCF context and grants no access to `self`’s state.
      !> @endnote
      procedure :: potpar
      procedure, private :: write_kspace_scf_dos_outputs
      !> Occupied-state site moments used by the optional k-space SCF branch.
      procedure :: compute_kspace_spin_moments_spinor
      final :: destructor
   end type self

   interface
      module subroutine quantify_pauli_projection(this)
         class(self), intent(inout) :: this
      end subroutine quantify_pauli_projection

      module subroutine write_kspace_scf_dos_outputs(this, reciprocal_obj)
         class(self), intent(inout) :: this
         type(reciprocal), intent(in) :: reciprocal_obj
      end subroutine write_kspace_scf_dos_outputs

      module subroutine compute_kspace_spin_moments_spinor(this, reciprocal_obj, site_mom)
         class(self), intent(in) :: this
         type(reciprocal), intent(in) :: reciprocal_obj
         real(rp), intent(out) :: site_mom(3, this%lattice%nrec)
      end subroutine compute_kspace_spin_moments_spinor

      module subroutine finalize_kspace_scf_state(this)
         class(self), intent(inout) :: this
      end subroutine finalize_kspace_scf_state

      module subroutine write_kspace_scf_state_artifact(this, filename)
         class(self), intent(in) :: this
         character(len=*), intent(in) :: filename
      end subroutine write_kspace_scf_state_artifact

   end interface

   interface self
      procedure :: constructor
   end interface self

contains

   subroutine initialize_fixed_u_trace(this)
      class(self), intent(inout) :: this
      character(len=1024) :: trace_path
      character(len=128) :: trace_value
      character(len=128) :: root_run_id
      integer :: ios, unit, root_trace_enabled, root_selected_iteration, root_capture_only, root_a2_scan
      integer :: root_moments_only, root_ql_only
      real(rp) :: root_trace_lambda

      this%fixed_u_trace_enabled = .false.
      this%fixed_u_trace_filename = ''
      this%fixed_u_trace_run_id = 'run'
      this%fixed_u_trace_lambda = 1.0_rp
      this%fixed_u_trace_selected_iteration = 0
      this%fixed_u_trace_moments_only = .false.
      this%fixed_u_trace_ql_only = .false.

      ! Rank zero owns environment parsing and file creation.  Every rank
      ! reaches the broadcasts below, including the no-trace and file-open
      ! failure cases, so optional M1 paths cannot disagree across ranks.
      root_trace_enabled = 0
      root_selected_iteration = 0
      root_capture_only = 0
      root_a2_scan = 0
      root_moments_only = 0
      root_ql_only = 0
      root_trace_lambda = 1.0_rp
      root_run_id = 'run'
      if (rank == 0) then
         trace_path = ''
         trace_value = ''
         call get_environment_variable('RSLMTO_FIXED_U_TRACE', trace_path, status=ios)
         if (ios == 0 .and. len_trim(trace_path) > 0 .and. &
             this%control%nsp == 1 .and. this%control%lmax == 2 .and. &
             (.not. this%hamiltonian%hoh .or. this%hamiltonian%fixed_u_test_allow_hoh_w1_u0)) then
            ! Explicit trace opt-in also enables the M0 no-U reference.
            ! Ordinary no-U runs still take the untouched path when the
            ! environment variable is absent. Other Hubbard modes remain out.
            if (this%hamiltonian%hubbard_u_general_check) then
               if (trim(this%hamiltonian%hubbard_u_potential_form) == 'liechtenstein') root_trace_enabled = 1
            else if (.not. (this%hamiltonian%hubbard_u_sc_check .or. this%hamiltonian%hubbard_v_check) .and. &
                     (.not. this%hamiltonian%hubbard_u_impurity_check .or. &
                      (trim(this%control%calctype) == 'I' .and. this%control%nsp == 1 .and. &
                       this%control%lmax == 2 .and. &
                       trim(this%hamiltonian%hubbard_radial_moments) == 'spectral'))) then
               root_trace_enabled = 1
            end if
         end if

         if (root_trace_enabled == 1) then
            call get_environment_variable('RSLMTO_FIXED_U_RUN_ID', trace_value, status=ios)
            if (ios == 0 .and. len_trim(trace_value) > 0) root_run_id = trim(trace_value)
            call get_environment_variable('RSLMTO_FIXED_U_LAMBDA', trace_value, status=ios)
            if (ios == 0 .and. len_trim(trace_value) > 0) then
               read(trace_value, *, iostat=ios) root_trace_lambda
               if (ios /= 0) root_trace_enabled = 0
            end if
            ! Zero selects rolling-last mode: replace the trace file at each
            ! iteration start, leaving the last completed iteration on disk.
            call get_environment_variable('RSLMTO_FIXED_U_TRACE_ITERATION', trace_value, status=ios)
            if (ios == 0 .and. len_trim(trace_value) > 0) then
               read(trace_value, *, iostat=ios) root_selected_iteration
               if (ios /= 0 .or. root_selected_iteration < 0) root_trace_enabled = 0
            end if
            trace_value = ''
            call get_environment_variable('RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY', trace_value, status=ios)
            if (ios == 0 .and. (trim(trace_value) == '1' .or. trim(trace_value) == 'true')) root_capture_only = 1
            trace_value = ''
            call get_environment_variable('RSLMTO_FIXED_U_A2_SCAN', trace_value, status=ios)
            if (ios == 0 .and. (trim(trace_value) == '1' .or. trim(trace_value) == 'true')) root_a2_scan = 1
            trace_value = ''
            call get_environment_variable('RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY', trace_value, status=ios)
            if (ios == 0 .and. (trim(trace_value) == '1' .or. trim(trace_value) == 'true')) then
               root_moments_only = 1
            end if
            trace_value = ''
            call get_environment_variable('RSLMTO_FIXED_U_TRACE_QL_ONLY', trace_value, status=ios)
            if (ios == 0 .and. (trim(trace_value) == '1' .or. trim(trace_value) == 'true')) root_ql_only = 1
         end if

         if (root_trace_enabled == 1) then
            open(newunit=unit, file=trim(trace_path), status='replace', action='write', iostat=ios)
            if (ios == 0) then
               write(unit, '(a)', iostat=ios) &
                  'kind,run_id,lambda,iteration,stage,site,type,l,spin,k,band,row,col,order,weight,energy,scalar,re,im'
               close(unit, iostat=ios)
            end if
            if (ios == 0) then
               this%fixed_u_trace_filename = trim(trace_path)
            else
               root_trace_enabled = 0
               call g_logger%warning('Could not create RSLMTO_FIXED_U_TRACE output; trace disabled.', __FILE__, __LINE__)
            end if
         end if
      end if

#ifdef USE_MPI
      call MPI_BCAST(root_trace_enabled, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(root_selected_iteration, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(root_capture_only, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(root_a2_scan, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(root_moments_only, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(root_ql_only, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(root_trace_lambda, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(root_run_id, len(root_run_id), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
#endif
      this%fixed_u_trace_enabled = root_trace_enabled == 1
      this%fixed_u_trace_capture_only = root_capture_only == 1
      this%fixed_u_trace_moments_only = root_moments_only == 1
      this%fixed_u_trace_ql_only = root_ql_only == 1
      this%fixed_u_trace_a2_scan_enabled = root_a2_scan == 1
      this%fixed_u_trace_selected_iteration = root_selected_iteration
      this%fixed_u_trace_lambda = root_trace_lambda
      this%fixed_u_trace_run_id = root_run_id
   end subroutine initialize_fixed_u_trace

   subroutine reset_fixed_u_trace(this)
      class(self), intent(in) :: this
      integer :: unit, ios
      if (.not. this%fixed_u_trace_enabled .or. rank /= 0) return
      if (this%fixed_u_trace_ql_only) return
      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='replace', action='write', iostat=ios)
      if (ios /= 0) return
      write(unit, '(a)') 'kind,run_id,lambda,iteration,stage,site,type,l,spin,k,band,row,col,order,weight,energy,scalar,re,im'
      close(unit)
   end subroutine reset_fixed_u_trace

   subroutine write_fixed_u_csv_row(this, unit, kind, iteration, stage, site, type_id, ell, spin, k, band, &
                                    row, col, order, weight, energy_value, scalar, re, im)
      class(self), intent(in) :: this
      integer, intent(in) :: unit, iteration, site, type_id, ell, spin, k, band, row, col, order
      character(len=*), intent(in) :: kind, stage
      real(rp), intent(in) :: weight, energy_value, scalar, re, im
      ! QL-only is the safe history mode used when diagnostic replay helpers
      ! are excluded.  Keep its artifact strictly limited to the accepted
      ! P_OUT occupations, including writes made by independent trace hooks.
      if (this%fixed_u_trace_ql_only .and. (trim(kind) /= 'ql' .or. trim(stage) /= 'P_OUT')) return
      write(unit, '(a,",",a,",",es24.16,",",i0,",",a,",",9(i0,","),es24.16,",",es24.16,",",es24.16,",",es24.16,",",es24.16)') &
         trim(kind), trim(this%fixed_u_trace_run_id), this%fixed_u_trace_lambda, iteration, trim(stage), &
         site, type_id, ell, spin, k, band, row, col, order, weight, energy_value, scalar, re, im
   end subroutine write_fixed_u_csv_row

   subroutine write_fixed_u_trace(this, stage, m1_out, pl_diag_out)
      class(self), intent(inout) :: this
      character(len=*), intent(in) :: stage
      complex(rp), intent(in), optional :: m1_out(:, :, :, :, :)
      real(rp), intent(in), optional :: pl_diag_out(:, :, :, :)
      integer :: unit, ios, ia, l, spin, i, j, site_type, plusbulk
      complex(rp) :: z, m1_trace, test_operator(5, 5)
      real(rp) :: center, q0, q2, test_q0_trace, test_trace_wp

      if (.not. this%fixed_u_trace_enabled .or. rank /= 0) return
      if (this%fixed_u_trace_selected_iteration > 0 .and. &
          this%fixed_u_trace_iteration /= this%fixed_u_trace_selected_iteration) return
      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios /= 0) return

      if (this%fixed_u_trace_ql_only) then
         if (trim(stage) == 'P_OUT') then
            do ia = 1, this%lattice%nrec
               plusbulk = this%lattice%nbulk + ia
               do l = 0, this%symbolic_atom(plusbulk)%potential%lmax
                  do spin = 1, 2
                     call this%write_fixed_u_csv_row(unit, 'ql', this%fixed_u_trace_iteration, 'P_OUT', ia - 1, &
                        this%lattice%iz(ia), l, spin, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, &
                        this%symbolic_atom(plusbulk)%potential%ql(1, l, spin), 0.0_rp, 0.0_rp)
                  end do
               end do
            end do
         end if
         close(unit)
         return
      end if

      test_operator = cmplx(0.0_rp, 0.0_rp, rp)
      if (this%fixed_u_trace_moments_only .and. this%hamiltonian%hubbard_test_operator_check) &
         call this%hamiltonian%get_hubbard_test_operator(test_operator)

      call this%write_fixed_u_csv_row(unit, 'ldm_availability', this%fixed_u_trace_iteration, trim(stage), &
         -1, -1, -1, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, &
         merge(1.0_rp, 0.0_rp, this%hamiltonian%hubbard_u_general_check), 0.0_rp, 0.0_rp)
      call this%write_fixed_u_csv_row(unit, 'vin_availability', this%fixed_u_trace_iteration, trim(stage), &
         -1, -1, -1, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, &
         merge(1.0_rp, 0.0_rp, allocated(this%hamiltonian%hubbard_u_pot)), 0.0_rp, 0.0_rp)

      do ia = 1, this%lattice%nrec
         plusbulk = this%lattice%nbulk + ia
         do l = 0, this%symbolic_atom(plusbulk)%potential%lmax
            do spin = 1, 2
               center = this%symbolic_atom(plusbulk)%potential%gravity_center(l + 1, spin) + &
                        this%symbolic_atom(plusbulk)%potential%vmad
               if (present(pl_diag_out)) then
                  if (l + 1 <= size(pl_diag_out, 2) .and. spin <= size(pl_diag_out, 3)) then
                     do i = 1, min(6, size(pl_diag_out, 4))
                        call this%write_fixed_u_csv_row(unit, 'pl_map_input', this%fixed_u_trace_iteration, &
                           'CALCULATE_PL', ia - 1, this%lattice%iz(ia), l, spin, -1, -1, -1, -1, i - 1, &
                           0.0_rp, 0.0_rp, pl_diag_out(ia, l + 1, spin, i), 0.0_rp, 0.0_rp)
                     end do
                  end if
               end if
               if (this%hamiltonian%hubbard_u_general_check) then
                  do i = 1, 2*l + 1
                     do j = 1, 2*l + 1
                        z = this%symbolic_atom(plusbulk)%potential%ldm(l + 1, spin, i, j)
                        call this%write_fixed_u_csv_row(unit, 'ldm', this%fixed_u_trace_iteration, trim(stage), ia - 1, &
                           this%lattice%iz(ia), l, spin, -1, -1, i - 1, j - 1, 0, 0.0_rp, 0.0_rp, 0.0_rp, &
                           real(z, rp), aimag(z))
                     end do
                  end do
               end if
               if (this%fixed_u_trace_moments_only .and. this%hamiltonian%hubbard_test_operator_check .and. &
                   trim(stage) == 'P_OUT' .and. l == 2) then
                  test_q0_trace = 0.0_rp
                  test_trace_wp = 0.0_rp
                  do i = 1, 5
                     do j = 1, 5
                        z = this%symbolic_atom(plusbulk)%potential%ldm(l + 1, spin, i, j)
                        call this%write_fixed_u_csv_row(unit, 'test_density_matrix', &
                           this%fixed_u_trace_iteration, 'P_OUT', ia - 1, this%lattice%iz(ia), l, spin - 1, &
                           -1, -1, i - 1, j - 1, 0, 0.0_rp, 0.0_rp, 0.0_rp, real(z, rp), aimag(z))
                        if (i == j) test_q0_trace = test_q0_trace + real(z, rp)
                        if (this%hamiltonian%hubbard_test_spin == 0 .or. &
                            this%hamiltonian%hubbard_test_spin == spin) then
                           test_trace_wp = test_trace_wp + real(test_operator(i, j) * &
                              this%symbolic_atom(plusbulk)%potential%ldm(l + 1, spin, j, i), rp)
                        end if
                     end do
                  end do
                  call this%write_fixed_u_csv_row(unit, 'test_q0_trace', this%fixed_u_trace_iteration, &
                     'P_OUT', ia - 1, this%lattice%iz(ia), l, spin - 1, -1, -1, -1, -1, 0, &
                     0.0_rp, 0.0_rp, test_q0_trace, 0.0_rp, 0.0_rp)
                  call this%write_fixed_u_csv_row(unit, 'test_trace_wp', this%fixed_u_trace_iteration, &
                     'P_OUT', ia - 1, this%lattice%iz(ia), l, spin - 1, -1, -1, -1, -1, 0, &
                     0.0_rp, 0.0_rp, test_trace_wp, 0.0_rp, 0.0_rp)
               end if
               if (present(m1_out)) then
                  if (l < size(m1_out, 2)) then
                     m1_trace = (0.0_rp, 0.0_rp)
                     do i = 1, min(2*l + 1, size(m1_out, 4))
                        do j = 1, min(2*l + 1, size(m1_out, 5))
                           z = m1_out(ia, l + 1, spin, i, j)
                           if (i == j) m1_trace = m1_trace + z
                           call this%write_fixed_u_csv_row(unit, 'ldm_m1', this%fixed_u_trace_iteration, 'M1_OUT', ia - 1, &
                              this%lattice%iz(ia), l, spin, -1, -1, i - 1, j - 1, 0, 0.0_rp, 0.0_rp, 0.0_rp, &
                              real(z, rp), aimag(z))
                        end do
                     end do
                     call this%write_fixed_u_csv_row(unit, 'm1_trace_consistency', &
                        this%fixed_u_trace_iteration, 'M1_OUT', ia - 1, this%lattice%iz(ia), l, spin, &
                        -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, &
                        real(m1_trace, rp) - this%symbolic_atom(plusbulk)%potential%ql(1, l, spin)*center, &
                        aimag(m1_trace), 0.0_rp)
                  end if
               end if
            end do
            do spin = 1, 2
               center = this%symbolic_atom(plusbulk)%potential%gravity_center(l + 1, spin) + &
                        this%symbolic_atom(plusbulk)%potential%vmad
               q0 = this%symbolic_atom(plusbulk)%potential%ql(1, l, spin)
               q2 = this%symbolic_atom(plusbulk)%potential%ql(3, l, spin)
               call this%write_fixed_u_csv_row(unit, 'ql', this%fixed_u_trace_iteration, trim(stage), ia - 1, &
                  this%lattice%iz(ia), l, spin, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, q0, 0.0_rp, 0.0_rp)
               call this%write_fixed_u_csv_row(unit, 'radial_moment', this%fixed_u_trace_iteration, trim(stage), ia - 1, &
                  this%lattice%iz(ia), l, spin, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, q0, 0.0_rp, 0.0_rp)
               call this%write_fixed_u_csv_row(unit, 'radial_moment', this%fixed_u_trace_iteration, trim(stage), ia - 1, &
                  this%lattice%iz(ia), l, spin, -1, -1, -1, -1, 1, 0.0_rp, center, q0*center, 0.0_rp, 0.0_rp)
               call this%write_fixed_u_csv_row(unit, 'radial_moment', this%fixed_u_trace_iteration, trim(stage), ia - 1, &
                  this%lattice%iz(ia), l, spin, -1, -1, -1, -1, 2, 0.0_rp, center**2 + q2, &
                  q0*(center**2 + q2), 0.0_rp, 0.0_rp)
               call this%write_fixed_u_csv_row(unit, 'radial_parameter', this%fixed_u_trace_iteration, trim(stage), ia - 1, &
                  this%lattice%iz(ia), l, spin, -1, -1, -1, -1, -1, 0.0_rp, &
                  this%symbolic_atom(plusbulk)%potential%gravity_center(l + 1, spin), 0.0_rp, 0.0_rp, 0.0_rp)
               call this%write_fixed_u_csv_row(unit, 'radial_parameter', this%fixed_u_trace_iteration, trim(stage), ia - 1, &
                  this%lattice%iz(ia), l, spin, -1, -1, -1, -1, 3, 0.0_rp, 0.0_rp, q2, 0.0_rp, 0.0_rp)
            end do
         end do
      end do

      if (allocated(this%hamiltonian%hubbard_u_pot)) then
         do site_type = 1, size(this%hamiltonian%hubbard_u_pot, 3)
            do i = 1, size(this%hamiltonian%hubbard_u_pot, 1)
               do j = 1, size(this%hamiltonian%hubbard_u_pot, 2)
                  z = this%hamiltonian%hubbard_u_pot(i, j, site_type)
                  call this%write_fixed_u_csv_row(unit, 'v_hub', this%fixed_u_trace_iteration, trim(stage), &
                     -1, site_type, -1, 0, -1, -1, i - 1, j - 1, 0, 0.0_rp, 0.0_rp, 0.0_rp, &
                     real(z, rp), aimag(z))
               end do
            end do
         end do
      end if
      if (allocated(this%hamiltonian%ee) .and. .not. this%fixed_u_trace_moments_only) then
         do site_type = 1, size(this%hamiltonian%ee, 4)
            do i = 1, size(this%hamiltonian%ee, 1)
               do j = 1, size(this%hamiltonian%ee, 2)
                  z = this%hamiltonian%ee(i, j, 1, site_type)
                  call this%write_fixed_u_csv_row(unit, 'h_total', this%fixed_u_trace_iteration, trim(stage), &
                     -1, site_type, -1, 0, -1, -1, i - 1, j - 1, 0, 0.0_rp, 0.0_rp, 0.0_rp, &
                     real(z, rp), aimag(z))
               end do
            end do
         end do
      end if
      if (stage == 'Pin_H' .and. allocated(this%hamiltonian%ee) .and. .not. this%fixed_u_trace_moments_only) &
         call this%write_fixed_u_rs_bonds(unit)
      if (stage == 'P_OUT' .and. .not. this%use_kspace) call this%write_fixed_u_rs_green(unit)
      close(unit)
   end subroutine write_fixed_u_trace

   subroutine write_fixed_u_rs_bonds(this, unit)
      class(self), intent(in) :: this
      integer, intent(in) :: unit
      integer :: ntype, isite, ia, ja, jsite, iclu, ineigh, nr, i, j, coord
      complex(rp) :: z
      real(rp) :: r_direct(3)

      ! Export fractional direct vectors without invoking clusba: clusba
      ! updates lattice%sbarvec, so even a diagnostic trace must not call it.

      if (.not. allocated(this%lattice%irec) .or. .not. allocated(this%lattice%iz)) return
      if (size(this%lattice%irec) < this%lattice%nrec) return
      if (.not. allocated(this%lattice%atlist) .or. .not. allocated(this%lattice%nn)) return

      do ntype = 1, this%lattice%ntype
         if (ntype > size(this%lattice%atlist)) cycle
         ia = this%lattice%atlist(ntype)
         if (ia < 1 .or. ia > size(this%lattice%nn, 1)) cycle
         if (size(this%lattice%nn, 2) < 1) cycle
         nr = min(this%lattice%nn(ia, 1), size(this%lattice%nn, 2), size(this%hamiltonian%ee, 3))
         do ineigh = 1, nr
            do i = 1, size(this%hamiltonian%ee, 1)
               do j = 1, size(this%hamiltonian%ee, 2)
                  z = this%hamiltonian%ee(i, j, ineigh, ntype)
                  call this%write_fixed_u_csv_row(unit, 'rs_hopping', this%fixed_u_trace_iteration, 'RS_OPERATOR', &
                     -1, ntype - 1, -1, -1, -1, ineigh - 1, i - 1, j - 1, 0, 0.0_rp, 0.0_rp, 0.0_rp, &
                     real(z, rp), aimag(z))
               end do
            end do
         end do
      end do

      do isite = 1, this%lattice%nrec
         iclu = this%lattice%irec(isite)
         if (iclu < 1 .or. iclu > size(this%lattice%iz)) cycle
         ntype = this%lattice%iz(iclu)
         if (ntype < 1 .or. ntype > this%lattice%ntype) cycle
         if (ntype > size(this%lattice%atlist)) cycle
         ia = this%lattice%atlist(ntype)
         if (ia < 1 .or. ia > size(this%lattice%nn, 1)) cycle
         if (size(this%lattice%nn, 2) < 1) cycle
         nr = min(this%lattice%nn(ia, 1), size(this%lattice%nn, 2))
         do ineigh = 1, nr
            if (ineigh == 1) then
               jsite = isite
               r_direct = 0.0_rp
            else
               ja = this%lattice%nn(ia, ineigh)
               if (ja < 1 .or. ja > this%lattice%kk) cycle
               jsite = this%lattice%iz(ja)
               if (jsite < 1 .or. jsite > this%lattice%nrec) cycle
               if (this%lattice%a_cart_inv_ready) then
                  r_direct = matmul(this%lattice%a_cart_inv, &
                                    this%lattice%cr(:, ja) - this%lattice%cr(:, ia))
               else
                  r_direct = matmul(inverse_3x3(this%lattice%a), &
                                    this%lattice%cr(:, ja) - this%lattice%cr(:, ia))
               end if
            end if
            do coord = 1, 3
               call this%write_fixed_u_csv_row(unit, 'rs_bond_map', this%fixed_u_trace_iteration, 'RS_OPERATOR', &
                  isite - 1, ntype - 1, -1, -1, -1, ineigh - 1, -1, jsite - 1, coord - 1, 0.0_rp, &
                  r_direct(coord), 0.0_rp, 0.0_rp, 0.0_rp)
            end do
         end do
      end do
   end subroutine write_fixed_u_rs_bonds

   subroutine write_fixed_u_radial_trace(this, atom, isite, v, rofi, rho_in, rho_out)
      class(self), intent(in) :: this
      class(symbolic_atom), intent(in) :: atom
      integer, intent(in) :: isite
      real(rp), intent(in) :: v(:, :), rofi(:)
      real(rp), intent(in) :: rho_in(:, :), rho_out(:, :)
      integer :: unit, ios, ir, spin

      if (.not. this%fixed_u_trace_enabled .or. rank /= 0) return
      if (this%fixed_u_trace_selected_iteration > 0 .and. &
          this%fixed_u_trace_iteration /= this%fixed_u_trace_selected_iteration) return
      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios /= 0) return
      do ir = 1, size(rofi)
         do spin = 1, size(v, 2)
            call this%write_fixed_u_csv_row(unit, 'radial_v', this%fixed_u_trace_iteration, 'POST_ATOMSC', &
               isite, this%lattice%iz(isite), -1, spin - 1, -1, -1, ir - 1, spin - 1, 0, 0.0_rp, &
               rofi(ir), v(ir, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_density_in', this%fixed_u_trace_iteration, 'POST_ATOMSC', &
               isite, this%lattice%iz(isite), -1, spin - 1, -1, -1, ir - 1, spin - 1, 0, 0.0_rp, &
               rofi(ir), rho_in(ir, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_density_out', this%fixed_u_trace_iteration, 'POST_ATOMSC', &
               isite, this%lattice%iz(isite), -1, spin - 1, -1, -1, ir - 1, spin - 1, 0, 0.0_rp, &
               rofi(ir), rho_out(ir, spin), 0.0_rp, 0.0_rp)
         end do
      end do
      call energy_row('sum_ev', atom%potential%SUMEV)
      call energy_row('sum_ec', atom%potential%SUMEC)
      call energy_row('ekin', atom%potential%EKIN)
      call energy_row('utot', atom%potential%UTOT)
      call energy_row('rhoeps', atom%potential%RHOEPS)
      call energy_row('etot_site', atom%potential%ETOT)
      close(unit)
   contains
      subroutine energy_row(label, energy_value)
         character(len=*), intent(in) :: label
         real(rp), intent(in) :: energy_value
         call this%write_fixed_u_csv_row(unit, label, this%fixed_u_trace_iteration, 'POST_ATOMSC', &
            isite, this%lattice%iz(isite), -1, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, energy_value, 0.0_rp, 0.0_rp)
      end subroutine energy_row
   end subroutine write_fixed_u_radial_trace

   subroutine write_fixed_u_radial_replay_input(this, atom, rofi, v, rho_in, pl, ql, ec_seed, ncore, &
                                               tol, vhrmax, b_fsm, txc, hyperfine)
      class(self), intent(in) :: this
      class(symbolic_atom), intent(in) :: atom
      real(rp), intent(in) :: rofi(:), v(:, :), rho_in(:, :), pl(:, :), ql(:, :, :), ec_seed(:), tol
      integer, intent(in) :: ncore
      real(rp), intent(in) :: vhrmax, b_fsm
      integer, intent(in) :: txc
      logical, intent(in) :: hyperfine
      integer :: unit, ios, ir, spin, ell, core_index

      if (.not. this%fixed_u_trace_enabled .or. rank /= 0) return
      if (this%fixed_u_trace_selected_iteration > 0 .and. &
          this%fixed_u_trace_iteration /= this%fixed_u_trace_selected_iteration) return
      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios /= 0) return
      call meta('radial_A', atom%A)
      call meta('radial_B', atom%B())
      call meta('radial_Z', atom%element%atomic_number)
      call meta('radial_tolerance', tol)
      call meta('radial_VHRMAX', vhrmax)
      call meta('radial_B_fsm', b_fsm)
      call meta('radial_txc', real(txc, rp))
      call meta('radial_hyperfine', merge(1.0_rp, 0.0_rp, hyperfine))
      call meta('radial_NR', real(size(rofi), rp))
      call meta('radial_lmax', real(atom%potential%lmax, rp))
      call meta('radial_nsp', real(size(v, 2), rp))
      call meta('radial_f_core', real(atom%element%f_core, rp))
      call meta('radial_replay_vmad', atom%potential%vmad)
      do ir = 1, size(rofi)
         do spin = 1, size(v, 2)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_mesh', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, -1, spin - 1, -1, -1, ir - 1, spin - 1, 0, &
               0.0_rp, rofi(ir), 0.0_rp, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_v', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, -1, spin - 1, -1, -1, ir - 1, spin - 1, 0, &
               0.0_rp, rofi(ir), v(ir, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_rho_in', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, -1, spin - 1, -1, -1, ir - 1, spin - 1, 0, &
               0.0_rp, rofi(ir), rho_in(ir, spin), 0.0_rp, 0.0_rp)
         end do
      end do
      do ell = 0, atom%potential%lmax
         do spin = 1, size(pl, 2)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_pl', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, pl(ell + 1, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_pnu', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, atom%potential%pnu(ell, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_enu', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, atom%potential%enu(ell, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_gravity_center', this%fixed_u_trace_iteration, &
               'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, atom%potential%gravity_center(ell + 1, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_band_center', this%fixed_u_trace_iteration, &
               'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, atom%potential%gravity_center(ell + 1, spin) + atom%potential%vmad, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_c', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, atom%potential%c(ell, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_qpar', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, atom%potential%qpar(ell, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_srdel', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, atom%potential%srdel(ell, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_ql', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 0, &
               0.0_rp, 0.0_rp, ql(1, ell + 1, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_ql', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 1, &
               0.0_rp, 0.0_rp, ql(2, ell + 1, spin), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_ql', this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
               this%fixed_u_trace_current_site - 1, -1, ell, spin - 1, -1, -1, ell, spin - 1, 2, &
               0.0_rp, 0.0_rp, ql(3, ell + 1, spin), 0.0_rp, 0.0_rp)
         end do
      end do
      if (size(ec_seed) < ncore) error stop 'fixed-U replay capture has incomplete core-energy seeds'
      do core_index = 1, ncore
         call this%write_fixed_u_csv_row(unit, 'radial_replay_ec_seed', this%fixed_u_trace_iteration, &
            'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, -1, -1, -1, -1, core_index - 1, &
            core_index - 1, -1, 0, 0.0_rp, 0.0_rp, ec_seed(core_index), 0.0_rp, 0.0_rp)
      end do
      close(unit)
   contains
      subroutine meta(kind, value)
         character(len=*), intent(in) :: kind
         real(rp), intent(in) :: value
         call this%write_fixed_u_csv_row(unit, kind, this%fixed_u_trace_iteration, 'RADIAL_REPLAY', &
            this%fixed_u_trace_current_site - 1, -1, -1, -1, -1, -1, -1, -1, 0, &
            0.0_rp, 0.0_rp, value, 0.0_rp, 0.0_rp)
      end subroutine meta
   end subroutine write_fixed_u_radial_replay_input

   subroutine write_fixed_u_radial_replay_output(this, rho, sumec, sumev)
      class(self), intent(in) :: this
      real(rp), intent(in) :: rho(:, :), sumec(:), sumev(:)
      integer :: unit, ios, ir, spin

      if (.not. this%fixed_u_trace_enabled .or. rank /= 0) return
      if (this%fixed_u_trace_selected_iteration > 0 .and. &
          this%fixed_u_trace_iteration /= this%fixed_u_trace_selected_iteration) return
      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios /= 0) return
      do ir = 1, size(rho, 1)
         do spin = 1, size(rho, 2)
            call this%write_fixed_u_csv_row(unit, 'radial_replay_rho_out', this%fixed_u_trace_iteration, &
               'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, -1, -1, spin - 1, -1, -1, &
               ir - 1, spin - 1, 0, 0.0_rp, 0.0_rp, rho(ir, spin), 0.0_rp, 0.0_rp)
         end do
      end do
      do spin = 1, size(rho, 2)
         call this%write_fixed_u_csv_row(unit, 'radial_replay_sumec', this%fixed_u_trace_iteration, &
            'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, -1, -1, spin - 1, -1, -1, &
            -1, -1, 0, 0.0_rp, 0.0_rp, sumec(spin), 0.0_rp, 0.0_rp)
         call this%write_fixed_u_csv_row(unit, 'radial_replay_sumev', this%fixed_u_trace_iteration, &
            'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, -1, -1, spin - 1, -1, -1, &
            -1, -1, 0, 0.0_rp, 0.0_rp, sumev(spin), 0.0_rp, 0.0_rp)
      end do
      close(unit)
   end subroutine write_fixed_u_radial_replay_output

      subroutine verify_fixed_u_radial_replay(this, atom, z, lmax, a, b, nr, rofi_base, v_base, rofi, v, rho_in, &
                                           rho_actual, vhrmax, b_fsm, rvh_actual, reps_actual, rmu_actual, &
                                           vnucl_actual, pl, ql, sumec_actual, sumev_actual, &
                                           ec_seed, ev_seed, ec_actual, ev_actual, ncore, nval, tol)
      class(self), intent(inout) :: this
      class(symbolic_atom), intent(in) :: atom
      integer, intent(in) :: lmax, nr, ncore, nval
      real(rp), intent(in) :: z, a, b, rofi_base(nr), v_base(nr, 2), rofi(nr), v(nr, 2)
      real(rp), intent(in) :: rho_in(nr, 2), rho_actual(nr, 2), vhrmax, b_fsm
      real(rp), intent(in) :: rvh_actual(2), reps_actual(2), rmu_actual(2), vnucl_actual
      real(rp), intent(in) :: pl(lmax + 1, 2), ql(3, lmax + 1, 2), sumec_actual(2), sumev_actual(2)
      real(rp), intent(in) :: ec_seed(:), ev_seed(:), ec_actual(:), ev_actual(:), tol
      type(self) :: replay_state
      type(symbolic_atom) :: replay_atom
      real(rp), allocatable :: rho_replay(:, :), sumec_replay(:), sumev_replay(:), ec_replay(:), ev_replay(:)
      integer :: unit, ios
      real(rp) :: rho_error, sumec_error, sumev_error, ev_error, ec_error

      if (.not. this%fixed_u_trace_enabled .or. rank /= 0) return
      if (this%fixed_u_trace_selected_iteration > 0 .and. &
          this%fixed_u_trace_iteration /= this%fixed_u_trace_selected_iteration) return
      if (size(ec_seed) < ncore .or. size(ev_seed) < nval .or. size(ec_actual) < ncore .or. size(ev_actual) < nval) &
         error stop 'verify_fixed_u_radial_replay: inconsistent seed/output sizes'

      ! Repeat the unchanged production NEWRHO call against the exact accepted
      ! atom and input arrays.  The independent scratch object and atom copy
      ! keep this diagnostic replay from changing the live SCF state.
      replay_atom = atom
      replay_state%control => this%control
      allocate(replay_state%vtn(nr, 2), replay_state%vzt(nr, 2), &
               replay_state%fun2(nr, lmax + 1, 2), replay_state%phi_amp(nr, lmax + 1, 2))
      replay_state%vtn = 0.0_rp
      replay_state%vzt = 0.0_rp
      replay_state%fun2 = 0.0_rp
      replay_state%phi_amp = 0.0_rp
      allocate(rho_replay(nr, 2), sumec_replay(2), sumev_replay(2), &
               ec_replay(max(ncore, 1)), ev_replay(max(nval, 1)))
      rho_replay = 0.0_rp
      sumec_replay = 0.0_rp
      sumev_replay = 0.0_rp
      if (ncore > 0) ec_replay(1:ncore) = ec_seed(1:ncore)
      if (nval > 0) ev_replay(1:nval) = ev_seed(1:nval)
      call replay_state%newrho(replay_atom, z, lmax, a, b, nr, rofi, v, rho_replay, pl, ql, &
                               sumec_replay, sumev_replay, ec_replay, ev_replay, tol, 2, 0)
      rho_error = maxval(abs(rho_replay - rho_actual))
      sumec_error = maxval(abs(sumec_replay - sumec_actual))
      sumev_error = maxval(abs(sumev_replay - sumev_actual))
      ev_error = 0.0_rp
      ec_error = 0.0_rp
      if (nval > 0) ev_error = maxval(abs(ev_replay(1:nval) - ev_actual(1:nval)))
      if (ncore > 0) ec_error = maxval(abs(ec_replay(1:ncore) - ec_actual(1:ncore)))

      deallocate(replay_state%vtn, replay_state%vzt, replay_state%fun2, replay_state%phi_amp)

      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios == 0) then
         call result('rho_max_abs', rho_error, 0)
         call result('sumec_max_abs_ry', sumec_error, 1)
         call result('sumev_max_abs_ry', sumev_error, 2)
         call result('ec_max_abs_ry', ec_error, 3)
         call result('ev_max_abs_ry', ev_error, 4)
         close(unit)
      end if
      call this%replay_fixed_u_m0_moments(atom, z, lmax, a, b, nr, rofi, v, rho_in, rho_actual, &
         pl, ql, ec_seed, ev_seed, ncore, nval, tol, rofi_base, v_base, vhrmax, b_fsm, &
         rvh_actual, reps_actual, rmu_actual, vnucl_actual, sumec_actual, sumev_actual)
   contains
      subroutine result(label, value, result_index)
         character(len=*), intent(in) :: label
         real(rp), intent(in) :: value
         integer, intent(in) :: result_index
         call this%write_fixed_u_csv_row(unit, 'radial_replay_'//trim(label), this%fixed_u_trace_iteration, &
            'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, -1, -1, -1, -1, -1, -1, -1, &
            result_index, 0.0_rp, 0.0_rp, value, 0.0_rp, 0.0_rp)
      end subroutine result
   end subroutine verify_fixed_u_radial_replay

   ! Capture M0 replay operands from one immutable outer-SCF iteration.  This
   ! routine is diagnostic-only: it reads P, M1, the Hamiltonian V_in, ql and
   ! the exact calculate_pl call inputs, then stores only replay scalars and
   ! snapshots.  It never writes to a live potential or Hamiltonian field.
   subroutine snapshot_fixed_u_m0_moments(this, m1, pl_diag, vin)
      class(self), intent(inout) :: this
      complex(rp), intent(in) :: m1(:, :, :, :, :), vin(:, :, :)
      real(rp), intent(in) :: pl_diag(:, :, :, :)
      complex(rp), allocatable :: density_override(:, :, :, :, :), vout(:, :, :)
      integer :: ia, plusbulk, itype, iclu, spin, i, j, k, rowoff, norb_l, unit, ios
      complex(rp) :: p_ij, t_complex, vm1_complex, v2p_complex, v2_ij, p_trace, m1_trace
      real(rp), allocatable :: im_resid(:, :), p_trace_error(:, :), m1_trace_error(:, :)
      real(rp) :: q0, center, q2, t, local_q2, t_vout, vm1_vout, v2_vout
      real(rp) :: t_pt, vm1_pt, v2_pt
      logical :: p_available

      this%fixed_u_m0_snapshot_valid = .false.
      if (.not. this%fixed_u_trace_enabled .or. rank /= 0) return
      if (this%fixed_u_trace_selected_iteration > 0 .and. &
          this%fixed_u_trace_iteration /= this%fixed_u_trace_selected_iteration) return
      if (size(m1, 1) /= this%lattice%nrec .or. size(pl_diag, 1) /= this%lattice%nrec .or. &
          size(vin, 1) < nb .or. size(vin, 2) < nb) return
      if (lmax_basis < 2) return
      if (.not. allocated(this%lattice%irec) .or. .not. allocated(this%lattice%iz)) return
      if (size(this%lattice%irec) < this%lattice%nrec) return

      if (allocated(this%fixed_u_m0_q0)) deallocate(this%fixed_u_m0_q0)
      if (allocated(this%fixed_u_m0_center)) deallocate(this%fixed_u_m0_center)
      if (allocated(this%fixed_u_m0_q2)) deallocate(this%fixed_u_m0_q2)
      if (allocated(this%fixed_u_m0_local_center)) deallocate(this%fixed_u_m0_local_center)
      if (allocated(this%fixed_u_m0_local_q2)) deallocate(this%fixed_u_m0_local_q2)
      if (allocated(this%fixed_u_m0_trvp)) deallocate(this%fixed_u_m0_trvp)
      if (allocated(this%fixed_u_m0_vout_center)) deallocate(this%fixed_u_m0_vout_center)
      if (allocated(this%fixed_u_m0_vout_q2)) deallocate(this%fixed_u_m0_vout_q2)
      if (allocated(this%fixed_u_m0_pt_center)) deallocate(this%fixed_u_m0_pt_center)
      if (allocated(this%fixed_u_m0_pt_q2)) deallocate(this%fixed_u_m0_pt_q2)
      if (allocated(this%fixed_u_m0_t_vout)) deallocate(this%fixed_u_m0_t_vout)
      if (allocated(this%fixed_u_m0_t_pt)) deallocate(this%fixed_u_m0_t_pt)
      if (allocated(this%fixed_u_m0_ql)) deallocate(this%fixed_u_m0_ql)
      if (allocated(this%fixed_u_m0_pl_diag)) deallocate(this%fixed_u_m0_pl_diag)
      if (allocated(this%fixed_u_m0_vin)) deallocate(this%fixed_u_m0_vin)
      allocate(this%fixed_u_m0_q0(this%lattice%nrec, 2), this%fixed_u_m0_center(this%lattice%nrec, 2), &
               this%fixed_u_m0_q2(this%lattice%nrec, 2), this%fixed_u_m0_local_center(this%lattice%nrec, 2), &
               this%fixed_u_m0_local_q2(this%lattice%nrec, 2), this%fixed_u_m0_trvp(this%lattice%nrec, 2), &
               this%fixed_u_m0_vout_center(this%lattice%nrec, 2), this%fixed_u_m0_vout_q2(this%lattice%nrec, 2), &
               this%fixed_u_m0_pt_center(this%lattice%nrec, 2), this%fixed_u_m0_pt_q2(this%lattice%nrec, 2), &
               this%fixed_u_m0_t_vout(this%lattice%nrec, 2), this%fixed_u_m0_t_pt(this%lattice%nrec, 2), &
               this%fixed_u_m0_ql(this%lattice%nrec, 3, lmax_basis + 1, 2), &
               this%fixed_u_m0_pl_diag(this%lattice%nrec, lmax_basis + 1, 2, 6), &
               this%fixed_u_m0_vin(size(vin, 1), size(vin, 2), size(vin, 3)))
      this%fixed_u_m0_ql = 0.0_rp
      this%fixed_u_m0_pl_diag = pl_diag
      this%fixed_u_m0_vin = vin
      this%fixed_u_m0_vout_center = 0.0_rp
      this%fixed_u_m0_vout_q2 = 0.0_rp
      this%fixed_u_m0_pt_center = 0.0_rp
      this%fixed_u_m0_pt_q2 = 0.0_rp
      this%fixed_u_m0_t_vout = 0.0_rp
      this%fixed_u_m0_t_pt = 0.0_rp
      allocate(im_resid(this%lattice%nrec, 2))
      allocate(p_trace_error(this%lattice%nrec, 2), m1_trace_error(this%lattice%nrec, 2))
      im_resid = 0.0_rp
      p_trace_error = 0.0_rp
      m1_trace_error = 0.0_rp
      p_available = this%hamiltonian%hubbard_u_general_check
      if (p_available) then
         do ia = 1, this%lattice%nrec
            plusbulk = this%lattice%nbulk + ia
            if (plusbulk < 1 .or. plusbulk > size(this%symbolic_atom)) then
               p_available = .false.
               exit
            end if
            if (.not. allocated(this%symbolic_atom(plusbulk)%potential%ldm)) then
               p_available = .false.
               exit
            end if
            if (size(this%symbolic_atom(plusbulk)%potential%ldm, 1) < 3 .or. &
                size(this%symbolic_atom(plusbulk)%potential%ldm, 2) < 2 .or. &
                size(this%symbolic_atom(plusbulk)%potential%ldm, 3) < 5 .or. &
                size(this%symbolic_atom(plusbulk)%potential%ldm, 4) < 5) then
               p_available = .false.
               exit
            end if
         end do
      end if
      this%fixed_u_m0_p_available = p_available

      ! Evaluate V(Pout) without changing the live Hamiltonian.  This is a
      ! deliberately wrong-pairing control for the frozen radial replay.
      allocate(density_override(this%lattice%ntype, 4, 2, 7, 7), vout(nb, nb, this%lattice%ntype))
      density_override = (0.0_rp, 0.0_rp)
      vout = (0.0_rp, 0.0_rp)
      if (p_available) then
         do ia = 1, this%lattice%nrec
            plusbulk = this%lattice%nbulk + ia
            iclu = this%lattice%irec(ia)
            if (iclu < 1 .or. iclu > size(this%lattice%iz)) then
               deallocate(im_resid, p_trace_error, m1_trace_error)
               deallocate(density_override, vout)
               return
            end if
            itype = this%lattice%iz(iclu)
            if (itype < 1 .or. itype > this%lattice%ntype) then
               deallocate(im_resid, p_trace_error, m1_trace_error)
               deallocate(density_override, vout)
               return
            end if
            do spin = 1, 2
               do i = 1, min(5, size(this%symbolic_atom(plusbulk)%potential%ldm, 3))
                  do j = 1, min(5, size(this%symbolic_atom(plusbulk)%potential%ldm, 4))
                     density_override(itype, 3, spin, i, j) = &
                        this%symbolic_atom(plusbulk)%potential%ldm(3, spin, i, j)
                  end do
               end do
            end do
         end do
         call this%hamiltonian%calculate_hubbard_u_potential_general( &
            density_override=density_override, potential_out=vout)
      end if

      do ia = 1, this%lattice%nrec
         plusbulk = this%lattice%nbulk + ia
         iclu = this%lattice%irec(ia)
         if (iclu < 1 .or. iclu > size(this%lattice%iz)) then
            deallocate(im_resid, p_trace_error, m1_trace_error)
            deallocate(density_override, vout)
            return
         end if
         itype = this%lattice%iz(iclu)
         if (itype < 1 .or. itype > size(vin, 3)) return
         do spin = 1, 2
            do i = 1, 3
               this%fixed_u_m0_ql(ia, i, 1:lmax_basis + 1, spin) = &
                  this%symbolic_atom(plusbulk)%potential%ql(i, 0:lmax_basis, spin)
            end do
            q0 = this%symbolic_atom(plusbulk)%potential%ql(1, 2, spin)
            center = this%symbolic_atom(plusbulk)%potential%gravity_center(3, spin) + &
                     this%symbolic_atom(plusbulk)%potential%vmad
            q2 = this%symbolic_atom(plusbulk)%potential%ql(3, 2, spin)
            if (q0 <= tiny(1.0_rp)) then
               this%fixed_u_m0_snapshot_valid = .false.
               return
            end if
            this%fixed_u_m0_q0(ia, spin) = q0
            this%fixed_u_m0_center(ia, spin) = center
            this%fixed_u_m0_q2(ia, spin) = q2
            t_complex = (0.0_rp, 0.0_rp)
            vm1_complex = (0.0_rp, 0.0_rp)
            v2p_complex = (0.0_rp, 0.0_rp)
            p_trace = (0.0_rp, 0.0_rp)
            m1_trace = (0.0_rp, 0.0_rp)
            rowoff = 4 + spin_off*(spin - 1)
            norb_l = 5
            do i = 1, norb_l
               m1_trace = m1_trace + m1(ia, 3, spin, i, i)
            end do
            if (p_available) then
            do i = 1, norb_l
               do j = 1, norb_l
                  p_ij = this%symbolic_atom(plusbulk)%potential%ldm(3, spin, i, j)
                  if (i == j) p_trace = p_trace + p_ij
                  t_complex = t_complex + vin(rowoff + i, rowoff + j, itype)* &
                     this%symbolic_atom(plusbulk)%potential%ldm(3, spin, j, i)
                  vm1_complex = vm1_complex + vin(rowoff + i, rowoff + j, itype)* &
                     (m1(ia, 3, spin, j, i) - center*this%symbolic_atom(plusbulk)%potential%ldm(3, spin, j, i))
                  v2_ij = (0.0_rp, 0.0_rp)
                  do k = 1, norb_l
                     v2_ij = v2_ij + vin(rowoff + i, rowoff + k, itype)*vin(rowoff + k, rowoff + j, itype)
                  end do
                  v2p_complex = v2p_complex + v2_ij*this%symbolic_atom(plusbulk)%potential%ldm(3, spin, j, i)
               end do
            end do
            end if
            t = real(t_complex, rp)
            local_q2 = q2 - 2.0_rp*real(vm1_complex, rp) + real(v2p_complex, rp) - t*t/q0
            im_resid(ia, spin) = max(abs(aimag(t_complex)), abs(aimag(vm1_complex)), abs(aimag(v2p_complex)))
            p_trace_error(ia, spin) = real(p_trace, rp) - q0
            m1_trace_error(ia, spin) = real(m1_trace, rp) - q0*center
            this%fixed_u_m0_trvp(ia, spin) = t
            this%fixed_u_m0_local_center(ia, spin) = center - t/q0
            this%fixed_u_m0_local_q2(ia, spin) = local_q2

            t_vout = 0.0_rp
            vm1_vout = 0.0_rp
            v2_vout = 0.0_rp
            t_pt = 0.0_rp
            vm1_pt = 0.0_rp
            v2_pt = 0.0_rp
            if (p_available) then
            do i = 1, norb_l
               do j = 1, norb_l
                  t_vout = t_vout + real(vout(rowoff+i, rowoff+j, itype)* &
                     this%symbolic_atom(plusbulk)%potential%ldm(3, spin, j, i), rp)
                  vm1_vout = vm1_vout + real(vout(rowoff+i, rowoff+j, itype)* &
                     (m1(ia, 3, spin, j, i) - center*this%symbolic_atom(plusbulk)%potential%ldm(3, spin, j, i)), rp)
                  t_pt = t_pt + real(vin(rowoff+i, rowoff+j, itype)* &
                     this%symbolic_atom(plusbulk)%potential%ldm(3, spin, i, j), rp)
                  vm1_pt = vm1_pt + real(vin(rowoff+i, rowoff+j, itype)* &
                     (m1(ia, 3, spin, j, i) - center*this%symbolic_atom(plusbulk)%potential%ldm(3, spin, i, j)), rp)
                  v2_ij = (0.0_rp, 0.0_rp)
                  do k = 1, norb_l
                     v2_ij = v2_ij + vout(rowoff+i, rowoff+k, itype)*vout(rowoff+k, rowoff+j, itype)
                  end do
                  v2_vout = v2_vout + real(v2_ij*this%symbolic_atom(plusbulk)%potential%ldm(3, spin, j, i), rp)
                  v2_ij = (0.0_rp, 0.0_rp)
                  do k = 1, norb_l
                     v2_ij = v2_ij + vin(rowoff+i, rowoff+k, itype)*vin(rowoff+k, rowoff+j, itype)
                  end do
                  v2_pt = v2_pt + real(v2_ij*this%symbolic_atom(plusbulk)%potential%ldm(3, spin, i, j), rp)
               end do
            end do
            end if
            this%fixed_u_m0_vout_center(ia, spin) = center - t_vout/q0
            this%fixed_u_m0_vout_q2(ia, spin) = q2 - 2.0_rp*vm1_vout + v2_vout - t_vout*t_vout/q0
            this%fixed_u_m0_pt_center(ia, spin) = center - t_pt/q0
            this%fixed_u_m0_pt_q2(ia, spin) = q2 - 2.0_rp*vm1_pt + v2_pt - t_pt*t_pt/q0
            this%fixed_u_m0_t_vout(ia, spin) = t_vout
            this%fixed_u_m0_t_pt(ia, spin) = t_pt
         end do
      end do
      this%fixed_u_m0_snapshot_iteration = this%fixed_u_trace_iteration
      this%fixed_u_m0_snapshot_valid = .true.

      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios /= 0) return
      do ia = 1, this%lattice%nrec
         do spin = 1, 2
            itype = this%lattice%ib(ia)
            call this%write_fixed_u_csv_row(unit, 'm0_carrier_input', this%fixed_u_m0_snapshot_iteration, &
               'FROZEN_REPLAY', ia - 1, itype - 1, 2, spin, -1, -1, -1, -1, 0, 0.0_rp, &
               this%fixed_u_m0_center(ia, spin), this%fixed_u_m0_q0(ia, spin), &
               this%fixed_u_m0_q2(ia, spin), 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'm0_p_availability', this%fixed_u_m0_snapshot_iteration, &
               'FROZEN_REPLAY', ia - 1, itype - 1, 2, spin, -1, -1, -1, -1, 0, 0.0_rp, &
               merge(1.0_rp, 0.0_rp, p_available), 0.0_rp, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'm0_carrier_local', this%fixed_u_m0_snapshot_iteration, &
               'FROZEN_REPLAY', ia - 1, itype - 1, 2, spin, -1, -1, -1, -1, 1, &
               this%fixed_u_m0_trvp(ia, spin), this%fixed_u_m0_local_center(ia, spin), &
               this%fixed_u_m0_q0(ia, spin), this%fixed_u_m0_local_q2(ia, spin), 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'm0_carrier_control', this%fixed_u_m0_snapshot_iteration, &
               'VOUT_CONTROL', ia - 1, itype - 1, 2, spin, -1, -1, -1, -1, 0, 0.0_rp, &
               this%fixed_u_m0_t_vout(ia, spin), this%fixed_u_m0_vout_center(ia, spin), &
               this%fixed_u_m0_q0(ia, spin), this%fixed_u_m0_vout_q2(ia, spin))
            call this%write_fixed_u_csv_row(unit, 'm0_carrier_control', this%fixed_u_m0_snapshot_iteration, &
               'P_TRANSPOSE_CONTROL', ia - 1, itype - 1, 2, spin, -1, -1, -1, -1, 0, 0.0_rp, &
               this%fixed_u_m0_t_pt(ia, spin), this%fixed_u_m0_pt_center(ia, spin), &
               this%fixed_u_m0_q0(ia, spin), this%fixed_u_m0_pt_q2(ia, spin))
            call this%write_fixed_u_csv_row(unit, 'm0_carrier_identity', this%fixed_u_m0_snapshot_iteration, &
               'FROZEN_REPLAY', ia - 1, itype - 1, 2, spin, -1, -1, -1, -1, 2, &
               0.0_rp, 0.0_rp, im_resid(ia, spin), 0.0_rp, 0.0_rp)
            if (p_available) then
               call this%write_fixed_u_csv_row(unit, 'm0_trace_consistency', this%fixed_u_m0_snapshot_iteration, &
                  'FROZEN_REPLAY', ia - 1, itype - 1, 2, spin, -1, -1, -1, -1, 0, &
                  0.0_rp, 0.0_rp, p_trace_error(ia, spin), 0.0_rp, 0.0_rp)
            else
               call this%write_fixed_u_csv_row(unit, 'm0_p_trace_unavailable', &
                  this%fixed_u_m0_snapshot_iteration, 'FROZEN_REPLAY', ia - 1, itype - 1, 2, spin, &
                  -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp)
            end if
            call this%write_fixed_u_csv_row(unit, 'm0_trace_consistency', this%fixed_u_m0_snapshot_iteration, &
               'FROZEN_REPLAY', ia - 1, itype - 1, 2, spin, -1, -1, -1, -1, 1, &
               0.0_rp, 0.0_rp, m1_trace_error(ia, spin), 0.0_rp, 0.0_rp)
         end do
      end do
      close(unit)
      deallocate(im_resid, p_trace_error, m1_trace_error)
      deallocate(density_override, vout)
   end subroutine snapshot_fixed_u_m0_moments

   ! Replay the same accepted radial state twice with production PL inversion
   ! and NEWRHO: legacy spectral moments and the proposed local-potential
   ! moments.  All edits are made on independent atom/bands copies.  This is
   ! an opt-in M0 diagnostic and cannot feed any live SCF quantity.
   subroutine replay_fixed_u_m0_moments(this, atom, z, lmax, a, b, nr, rofi, v, rho_in, rho_actual, &
                                        pl_actual, ql_actual, ec_seed, ev_seed, ncore, nval, tol, &
                                        rofi_base, v_base, vhrmax, b_fsm, rvh_actual, reps_actual, rmu_actual, &
                                        vnucl_actual, sumec_actual, sumev_actual)
      class(self), intent(inout) :: this
      class(symbolic_atom), intent(in) :: atom
      integer, intent(in) :: lmax, nr, ncore, nval
      real(rp), intent(in) :: z, a, b, rofi(nr), v(nr, 2), rho_in(nr, 2), rho_actual(nr, 2)
      real(rp), intent(in) :: rofi_base(nr), v_base(nr, 2), vhrmax, b_fsm
      real(rp), intent(in) :: rvh_actual(2), reps_actual(2), rmu_actual(2), vnucl_actual
      real(rp), intent(in) :: sumec_actual(2), sumev_actual(2)
      real(rp), intent(in) :: pl_actual(lmax + 1, 2), ql_actual(3, lmax + 1, 2)
      real(rp), intent(in) :: ec_seed(:), ev_seed(:), tol
      type(bands) :: pl_bands
      type(symbolic_atom), allocatable, target :: pl_atoms(:)
      real(rp), allocatable :: ql_raw(:, :, :), ql_local(:, :, :), ql_vout(:, :, :), ql_pt(:, :, :)
      real(rp), allocatable :: rho_case(:, :), sumec_case(:), sumev_case(:)
      real(rp), allocatable :: ec_case(:), ev_case(:)
      real(rp), allocatable :: pl_raw(:, :), pl_local(:, :), pl_a11(:, :, :)
      real(rp), allocatable :: baseline_ev(:)
      integer :: ia, plusbulk, l, spin, unit, ios, row, iv, nshift
      real(rp) :: pl_reproduction_error, ql_vs_live_error, q2_min
      real(rp) :: baseline_sumev(2), baseline_sumec(2)
      real(rp) :: component_snapshot(10), actual_component(10), output_component(10)
      real(rp) :: actual_rhovh, actual_rhomu
      real(rp) :: actual_rhoeps, actual_zvnucl, actual_ekin, actual_utot, actual_etot
      logical :: have_raw_energy

      if (.not. this%fixed_u_m0_snapshot_valid) return
      if (this%fixed_u_m0_snapshot_iteration /= this%fixed_u_trace_iteration) return
      if (this%fixed_u_trace_current_site < 1 .or. this%fixed_u_trace_current_site > this%lattice%nrec) return
      if (.not. allocated(this%fixed_u_m0_pl_diag) .or. .not. allocated(this%fixed_u_m0_ql)) return
      if (size(ec_seed) < ncore .or. size(ev_seed) < nval) return

      allocate(pl_atoms(size(this%symbolic_atom)))
      pl_atoms = this%symbolic_atom
      pl_bands%symbolic_atom => pl_atoms
      pl_bands%lattice => this%lattice
      pl_bands%control => this%control
      call load_captured_pl_inputs(pl_atoms)
      call pl_bands%calculate_pl()
      allocate(pl_raw(lmax + 1, 2), pl_local(lmax + 1, 2))
      ia = this%fixed_u_trace_current_site
      plusbulk = this%lattice%nbulk + ia
      do spin = 1, 2
         do l = 0, lmax
            pl_raw(l + 1, spin) = pl_atoms(plusbulk)%potential%pl(l, spin)
         end do
      end do

      nshift = fixed_u_m0_a11_nshift
      if (this%fixed_u_trace_a2_scan_enabled) nshift = fixed_u_m0_a2_nshift
      allocate(pl_a11(lmax + 1, 2, nshift))
      do iv = 1, nshift
         pl_atoms = this%symbolic_atom
         call load_captured_pl_inputs(pl_atoms)
         do spin = 1, 2
            pl_atoms(plusbulk)%potential%gravity_center(3, spin) = &
               this%fixed_u_m0_pl_diag(ia, 3, spin, 1) - fixed_u_m0_a2_v(iv)
         end do
         call pl_bands%calculate_pl()
         do spin = 1, 2
            do l = 0, lmax
               pl_a11(l + 1, spin, iv) = pl_atoms(plusbulk)%potential%pl(l, spin)
            end do
         end do
      end do
      pl_reproduction_error = maxval(abs(pl_raw - this%fixed_u_m0_pl_diag(ia, 1:lmax + 1, 1:2, 6)))

      ! Reset the scratch atom to its independent copy of the accepted live
      ! state, then restore the original calculate_pl inputs.  Only d-shell
      ! gravity_center differs in this second call.
      pl_atoms = this%symbolic_atom
      call load_captured_pl_inputs(pl_atoms)
      do spin = 1, 2
         pl_atoms(plusbulk)%potential%gravity_center(3, spin) = &
            this%fixed_u_m0_local_center(ia, spin) - pl_atoms(plusbulk)%potential%vmad
      end do
      call pl_bands%calculate_pl()
      do spin = 1, 2
         do l = 0, lmax
            pl_local(l + 1, spin) = pl_atoms(plusbulk)%potential%pl(l, spin)
         end do
      end do

      allocate(ql_raw(3, lmax + 1, 2), ql_local(3, lmax + 1, 2), &
               ql_vout(3, lmax + 1, 2), ql_pt(3, lmax + 1, 2))
      allocate(baseline_ev(max(nval, 1)))
      baseline_ev = 0.0_rp
      have_raw_energy = .false.
      baseline_sumev = 0.0_rp
      baseline_sumec = 0.0_rp
      ql_raw = this%fixed_u_m0_ql(ia, :, 1:lmax + 1, :)
      ql_local = ql_raw
      ql_vout = ql_raw
      ql_pt = ql_raw
      do spin = 1, 2
         ql_local(3, 3, spin) = this%fixed_u_m0_local_q2(ia, spin)
         ql_vout(3, 3, spin) = this%fixed_u_m0_vout_q2(ia, spin)
         ql_pt(3, 3, spin) = this%fixed_u_m0_pt_q2(ia, spin)
      end do
      ql_vs_live_error = maxval(abs(ql_raw - ql_actual))

      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios == 0) then
         call this%write_fixed_u_csv_row(unit, 'm0_replay_match', this%fixed_u_trace_iteration, 'FROZEN_INPUT', &
            ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 0, 0.0_rp, &
            pl_reproduction_error, ql_vs_live_error, 0.0_rp, 0.0_rp)
         q2_min = minval(this%fixed_u_m0_local_q2(ia, :))
         call this%write_fixed_u_csv_row(unit, 'm0_replay_q2_min', this%fixed_u_trace_iteration, 'LOCAL_PL', &
            ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 0, 0.0_rp, q2_min, 0.0_rp, 0.0_rp, 0.0_rp)
         close(unit)
      end if

      ! Replay the exact PL/QL supplied to this live NEWRHO call before using
      ! the unmixed spectral q0/q2 values below.  This is the same-state
      ! control for the copied radial solver; RAW_PL remains a separate
      ! diagnostic because its moments are from the current spectral output.
      call run_replay_case('LIVE_STATE', pl_actual, ql_actual)
      call run_replay_case('RAW_PL', pl_raw, ql_raw)
      do iv = 1, nshift
         call run_replay_case(a11_case_name(iv), pl_a11(:, :, iv), ql_raw, fixed_u_m0_a2_v(iv), iv)
      end do
      if (minval(this%fixed_u_m0_local_q2(ia, :)) >= 0.0_rp) then
         call run_replay_case('LOCAL_PL', pl_local, ql_local)
      else
         open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
              action='write', iostat=ios)
         if (ios == 0) then
            call this%write_fixed_u_csv_row(unit, 'm0_replay_skipped', this%fixed_u_trace_iteration, 'LOCAL_PL', &
               ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 0, 0.0_rp, &
               minval(this%fixed_u_m0_local_q2(ia, :)), 0.0_rp, 0.0_rp, 0.0_rp)
            close(unit)
         end if
      end if
      if (this%fixed_u_m0_p_available .and. minval(this%fixed_u_m0_vout_q2(ia, :)) >= 0.0_rp) &
         call run_replay_case('VOUT_CONTROL', pl_raw, ql_vout)
      if (this%fixed_u_m0_p_available .and. minval(this%fixed_u_m0_pt_q2(ia, :)) >= 0.0_rp) &
         call run_replay_case('P_TRANSPOSE_CONTROL', pl_raw, ql_pt)
      nullify(pl_bands%symbolic_atom, pl_bands%lattice)
      deallocate(pl_atoms, pl_raw, pl_local, pl_a11, ql_raw, ql_local, ql_vout, ql_pt, baseline_ev)
   contains
      function a11_case_name(shift_index) result(name)
         integer, intent(in) :: shift_index
         character(len=24) :: name
         if (shift_index <= fixed_u_m0_a11_nshift) then
            write(name, '("A1_SHIFT_",i2.2)') shift_index
         else
            write(name, '("A2_SHIFT_",i2.2)') shift_index
         end if
      end function a11_case_name

      function a11_stage_name(case_name) result(name)
         character(len=*), intent(in) :: case_name
         character(len=24) :: name
         name = case_name
      end function a11_stage_name

      logical function is_a11_shift_case(case_name)
         character(len=*), intent(in) :: case_name
         is_a11_shift_case = index(trim(case_name), 'A1_SHIFT_') == 1 .or. &
                             index(trim(case_name), 'A2_SHIFT_') == 1
      end function is_a11_shift_case

      subroutine load_captured_pl_inputs(atoms)
         type(symbolic_atom), intent(inout) :: atoms(:)
         integer :: site, pb, shell, channel
         do site = 1, this%lattice%nrec
            pb = this%lattice%nbulk + site
            do channel = 1, 2
               do shell = 0, lmax
                  atoms(pb)%potential%gravity_center(shell + 1, channel) = &
                     this%fixed_u_m0_pl_diag(site, shell + 1, channel, 1)
                  atoms(pb)%potential%c(shell, channel) = this%fixed_u_m0_pl_diag(site, shell + 1, channel, 2)
                  atoms(pb)%potential%qpar(shell, channel) = this%fixed_u_m0_pl_diag(site, shell + 1, channel, 3)
                  atoms(pb)%potential%srdel(shell, channel) = this%fixed_u_m0_pl_diag(site, shell + 1, channel, 4)
                  atoms(pb)%potential%pl(shell, channel) = this%fixed_u_m0_pl_diag(site, shell + 1, channel, 5)
               end do
            end do
         end do
      end subroutine load_captured_pl_inputs

      subroutine run_replay_case(case_name, pl_case, ql_case, a11_v, a11_index)
         character(len=*), intent(in) :: case_name
         real(rp), intent(in) :: pl_case(lmax + 1, 2), ql_case(3, lmax + 1, 2)
         real(rp), intent(in), optional :: a11_v
         integer, intent(in), optional :: a11_index
         type(self) :: replay_state
         type(symbolic_atom) :: replay_atom
      type(xc) :: replay_xc
      integer :: out_unit, out_ios, ir, channel, shell, moment, state_index, shift_index
      real(rp) :: center_input, shift_value, ev_raw_spin(2), ev_case_spin(2)
      real(rp) :: ev_raw_total, ev_case_total
      real(rp) :: charge(2), drdi, wgt
      logical :: a11_case
         shift_value = 0.0_rp
         shift_index = 0
         if (present(a11_v)) shift_value = a11_v
         if (present(a11_index)) shift_index = a11_index
         a11_case = trim(case_name) == 'RAW_PL' .or. is_a11_shift_case(case_name)
         allocate(rho_case(nr, 2), sumec_case(2), sumev_case(2), &
                  ec_case(max(ncore, 1)), ev_case(max(nval, 1)))
         allocate(replay_state%vtn(nr, 2), replay_state%vzt(nr, 2), &
                  replay_state%fun2(nr, lmax + 1, 2), replay_state%phi_amp(nr, lmax + 1, 2))
         replay_state%control => this%control
         replay_state%magnetic_seed_active = this%magnetic_seed_active
         replay_state%magnetic_seed_field = this%magnetic_seed_field
         replay_state%vtn = 0.0_rp
         replay_state%vzt = 0.0_rp
         replay_state%fun2 = 0.0_rp
         replay_state%phi_amp = 0.0_rp
         rho_case = 0.0_rp
         sumec_case = 0.0_rp
         sumev_case = 0.0_rp
         if (ncore > 0) ec_case(1:ncore) = ec_seed(1:ncore)
         if (nval > 0) ev_case(1:nval) = ev_seed(1:nval)
         replay_atom = atom
         call replay_state%newrho(replay_atom, z, lmax, a, b, nr, rofi, v, rho_case, pl_case, ql_case, &
                                  sumec_case, sumev_case, ec_case, ev_case, tol, 2, 0)
         if (trim(case_name) == 'RAW_PL') then
            baseline_sumev = sumev_case
            baseline_sumec = sumec_case
            if (nval > 0) baseline_ev(1:nval) = ev_case(1:nval)
            have_raw_energy = .true.
         end if
         replay_xc = xc(this%lattice%control)
         charge = 0.0_rp
         do channel = 1, 2
            do ir = 1, nr
               wgt = 2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp
               if (ir == 1 .or. ir == nr) wgt = 1.0_rp/3.0_rp
               drdi = a*(rofi_base(ir) + b)
               charge(channel) = charge(channel) + wgt*drdi*rho_case(ir, channel)
            end do
         end do
         open(newunit=out_unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
              action='write', iostat=out_ios)
         if (out_ios == 0) then
            ! Recompute both component sets on copied arrays and write them
            ! only after the trace stream is open. ATOMSC evaluates the
            ! actual electrostatic/XC terms on rho_in before NEWRHO; the
            ! output-rho set is a separately labelled diagnostic functional.
            call evaluate_frozen_components(this, replay_state, replay_xc, atom, z, a, b, nr, rofi_base, &
               v_base, vhrmax, b_fsm, rho_in, sumec_case, sumev_case, out_unit, ia, &
               'ATOMSC_INPUT', component_snapshot)
            ! Use the VNUCL saved by the actual ATOMSC call-site, not an
            ! inferred boundary value, for the direct production comparison.
            actual_rhovh = sum(rvh_actual)
            actual_rhomu = sum(rmu_actual)
            actual_rhoeps = sum(reps_actual)
            actual_zvnucl = -z*vnucl_actual
            actual_ekin = sum(sumev_actual)+sum(sumec_actual)-actual_rhovh-actual_rhomu
            actual_utot = 0.5_rp*(actual_rhovh+actual_zvnucl)
            actual_etot = actual_ekin+actual_utot+actual_rhoeps
            actual_component = [sum(sumev_actual), sum(sumec_actual), actual_rhovh, actual_rhomu, &
                                vnucl_actual, actual_zvnucl, actual_ekin, actual_utot, actual_rhoeps, actual_etot]
            call evaluate_frozen_components(this, replay_state, replay_xc, atom, z, a, b, nr, rofi_base, &
               v_base, vhrmax, b_fsm, rho_case, sumec_case, sumev_case, out_unit, ia, &
               'OUTPUT_RHO', output_component)
            ! Emit the exact ATOMSC site-energy reconstruction residual. This
            ! is a same-state replay check, not an energy correction.
            do shell = 0, 9
               call this%write_fixed_u_csv_row(out_unit, 'm0_replay_component_delta', &
                  this%fixed_u_trace_iteration, trim(case_name), ia-1, this%lattice%ib(ia)-1, 2, -1, &
                  -1, -1, -1, -1, shell, 0.0_rp, 0.0_rp, component_snapshot(shell+1)-actual_component(shell+1), &
                  component_snapshot(shell+1), actual_component(shell+1))
            end do
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_control', this%fixed_u_trace_iteration, &
               trim(case_name), ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 0, &
               0.0_rp, a, z, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_control', this%fixed_u_trace_iteration, &
               trim(case_name), ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 1, &
               0.0_rp, b, tol, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_control', this%fixed_u_trace_iteration, &
               trim(case_name), ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 2, &
               0.0_rp, real(nr, rp), real(ncore, rp), real(nval, rp), 0.0_rp)
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_control', this%fixed_u_trace_iteration, &
               trim(case_name), ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 3, &
               0.0_rp, shift_value, this%fixed_u_m0_trvp(ia, 1), &
               this%fixed_u_m0_trvp(ia, 2), 0.0_rp)
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_control', this%fixed_u_trace_iteration, &
               trim(case_name), ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 4, &
               0.0_rp, vhrmax, 0.0_rp, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_control', this%fixed_u_trace_iteration, &
               trim(case_name), ia - 1, this%lattice%ib(ia) - 1, 2, -1, -1, -1, -1, -1, 5, &
               0.0_rp, b_fsm, 0.0_rp, 0.0_rp, 0.0_rp)
            ! Serialize the exact PL tuple passed to the production scratch
            ! calculate_pl call and its returned PL for every shell/spin.
            do channel = 1, 2
               do shell = 0, lmax
                  call this%write_fixed_u_csv_row(out_unit, 'm0_replay_radial_parameter', &
                     this%fixed_u_trace_iteration, trim(case_name), ia - 1, this%lattice%ib(ia) - 1, &
                     shell, channel, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, &
                     atom%potential%c(shell, channel), 0.0_rp, 0.0_rp)
                  call this%write_fixed_u_csv_row(out_unit, 'm0_replay_radial_parameter', &
                     this%fixed_u_trace_iteration, trim(case_name), ia - 1, this%lattice%ib(ia) - 1, &
                     shell, channel, -1, -1, -1, -1, 1, 0.0_rp, 0.0_rp, &
                     atom%potential%qpar(shell, channel), 0.0_rp, 0.0_rp)
                  call this%write_fixed_u_csv_row(out_unit, 'm0_replay_radial_parameter', &
                     this%fixed_u_trace_iteration, trim(case_name), ia - 1, this%lattice%ib(ia) - 1, &
                     shell, channel, -1, -1, -1, -1, 2, 0.0_rp, 0.0_rp, &
                     atom%potential%srdel(shell, channel), 0.0_rp, 0.0_rp)
                  center_input = this%fixed_u_m0_pl_diag(ia, shell + 1, channel, 1)
                  if (shell == 2) then
                     select case (trim(case_name))
                     case ('LOCAL_PL')
                        center_input = this%fixed_u_m0_local_center(ia, channel) - atom%potential%vmad
                     case ('VOUT_CONTROL')
                        center_input = this%fixed_u_m0_vout_center(ia, channel) - atom%potential%vmad
                     case ('P_TRANSPOSE_CONTROL')
                        center_input = this%fixed_u_m0_pt_center(ia, channel) - atom%potential%vmad
                     end select
                     if (is_a11_shift_case(case_name)) &
                        center_input = this%fixed_u_m0_pl_diag(ia, shell + 1, channel, 1) - shift_value
                  end if
                  do moment = 0, 4
                     select case (moment)
                     case (0)
                        center_input = this%fixed_u_m0_pl_diag(ia, shell + 1, channel, 1)
                        if (shell == 2) then
                           select case (trim(case_name))
                           case ('LOCAL_PL')
                              center_input = this%fixed_u_m0_local_center(ia, channel) - atom%potential%vmad
                           case ('VOUT_CONTROL')
                              center_input = this%fixed_u_m0_vout_center(ia, channel) - atom%potential%vmad
                           case ('P_TRANSPOSE_CONTROL')
                              center_input = this%fixed_u_m0_pt_center(ia, channel) - atom%potential%vmad
                           end select
                           if (is_a11_shift_case(case_name)) &
                              center_input = this%fixed_u_m0_pl_diag(ia, shell + 1, channel, 1) - shift_value
                        end if
                     case (1)
                        center_input = this%fixed_u_m0_pl_diag(ia, shell + 1, channel, 2)
                     case (2)
                        center_input = this%fixed_u_m0_pl_diag(ia, shell + 1, channel, 3)
                     case (3)
                        center_input = this%fixed_u_m0_pl_diag(ia, shell + 1, channel, 4)
                     case (4)
                        center_input = this%fixed_u_m0_pl_diag(ia, shell + 1, channel, 5)
                     end select
                     call this%write_fixed_u_csv_row(out_unit, 'm0_replay_pl_input', this%fixed_u_trace_iteration, &
                        trim(case_name), ia - 1, this%lattice%ib(ia) - 1, shell, channel, -1, -1, -1, -1, moment, &
                        0.0_rp, 0.0_rp, center_input, 0.0_rp, 0.0_rp)
                  end do
                  call this%write_fixed_u_csv_row(out_unit, 'm0_replay_pl_return', this%fixed_u_trace_iteration, &
                     trim(case_name), ia - 1, this%lattice%ib(ia) - 1, shell, channel, -1, -1, -1, -1, 0, &
                     0.0_rp, 0.0_rp, pl_case(shell + 1, channel), 0.0_rp, 0.0_rp)
                  do moment = 1, 3
                     call this%write_fixed_u_csv_row(out_unit, 'm0_replay_ql_input', this%fixed_u_trace_iteration, &
                        trim(case_name), ia - 1, this%lattice%ib(ia) - 1, shell, channel, -1, -1, -1, -1, moment - 1, &
                        0.0_rp, 0.0_rp, ql_case(moment, shell + 1, channel), 0.0_rp, 0.0_rp)
                  end do
               end do
            end do
            ! Each case records the same immutable radial solver inputs. The
            ! independent checker compares these rows exactly across cases.
            do ir = 1, nr
               call this%write_fixed_u_csv_row(out_unit, 'm0_replay_rofi', this%fixed_u_trace_iteration, &
                  trim(case_name), ia - 1, this%lattice%ib(ia) - 1, -1, -1, -1, -1, ir - 1, -1, 0, &
                  0.0_rp, 0.0_rp, rofi(ir), 0.0_rp, 0.0_rp)
               do channel = 1, 2
                  call this%write_fixed_u_csv_row(out_unit, 'm0_replay_v_input', this%fixed_u_trace_iteration, &
                     trim(case_name), ia - 1, this%lattice%ib(ia) - 1, -1, channel, -1, -1, ir - 1, -1, 0, &
                     0.0_rp, 0.0_rp, v(ir, channel), 0.0_rp, 0.0_rp)
                  call this%write_fixed_u_csv_row(out_unit, 'm0_replay_rho_input', this%fixed_u_trace_iteration, &
                     trim(case_name), ia - 1, this%lattice%ib(ia) - 1, -1, channel, -1, -1, ir - 1, -1, 0, &
                     0.0_rp, 0.0_rp, rho_in(ir, channel), 0.0_rp, 0.0_rp)
               end do
            end do
            do row = 1, ncore
               call this%write_fixed_u_csv_row(out_unit, 'm0_replay_ec_seed', this%fixed_u_trace_iteration, &
                  trim(case_name), ia - 1, this%lattice%ib(ia) - 1, -1, -1, -1, -1, row - 1, -1, 0, &
                  0.0_rp, 0.0_rp, ec_seed(row), 0.0_rp, 0.0_rp)
            end do
            do row = 1, nval
               call this%write_fixed_u_csv_row(out_unit, 'm0_replay_ev_seed', this%fixed_u_trace_iteration, &
                  trim(case_name), ia - 1, this%lattice%ib(ia) - 1, -1, -1, -1, -1, row - 1, -1, 0, &
                  0.0_rp, 0.0_rp, ev_seed(row), 0.0_rp, 0.0_rp)
            end do
            do channel = 1, 2
               call this%write_fixed_u_csv_row(out_unit, 'm0_replay_energy', this%fixed_u_trace_iteration, &
                  trim(case_name), ia - 1, this%lattice%ib(ia) - 1, 2, channel, -1, -1, -1, -1, 0, &
                  0.0_rp, sumec_case(channel), sumev_case(channel), 0.0_rp, 0.0_rp)
               if (have_raw_energy) then
                  call this%write_fixed_u_csv_row(out_unit, 'm0_replay_delta', this%fixed_u_trace_iteration, &
                     trim(case_name), ia - 1, this%lattice%ib(ia) - 1, 2, channel, -1, -1, -1, -1, 0, &
                     0.0_rp, sumec_case(channel)-baseline_sumec(channel), &
                     sumev_case(channel)-baseline_sumev(channel), 0.0_rp, 0.0_rp)
               end if
            end do
            if (a11_case) then
               ev_raw_spin = 0.0_rp
               ev_case_spin = 0.0_rp
               ev_raw_total = 0.0_rp
               ev_case_total = 0.0_rp
               state_index = 0
               do channel = 1, 2
                  do shell = 0, lmax
                     if (ql_case(1, shell + 1, channel) < 1.0e-5_rp) cycle
                     state_index = state_index + 1
                     if (state_index > nval) error stop 'A1 trace: NEWRHO EV ordering exceeds nval'
                     ev_raw_spin(channel) = ev_raw_spin(channel) + baseline_ev(state_index)
                     ev_case_spin(channel) = ev_case_spin(channel) + ev_case(state_index)
                     call this%write_fixed_u_csv_row(out_unit, 'm0_a11_ev', this%fixed_u_trace_iteration, &
                        a11_stage_name(case_name), ia-1, this%lattice%ib(ia)-1, shell, channel, -1, &
                        state_index-1, -1, -1, shift_index, shift_value, baseline_ev(state_index), &
                        ev_case(state_index), ev_case(state_index)-baseline_ev(state_index), 0.0_rp)
                  end do
                  call this%write_fixed_u_csv_row(out_unit, 'm0_a11_ev_sum', this%fixed_u_trace_iteration, &
                     a11_stage_name(case_name), ia-1, this%lattice%ib(ia)-1, 2, channel, -1, -1, -1, -1, &
                     shift_index, shift_value, ev_raw_spin(channel), ev_case_spin(channel), &
                     ev_case_spin(channel)-ev_raw_spin(channel), 0.0_rp)
                  ev_raw_total = ev_raw_total + ev_raw_spin(channel)
                  ev_case_total = ev_case_total + ev_case_spin(channel)
               end do
               call this%write_fixed_u_csv_row(out_unit, 'm0_a11_ev_sum', this%fixed_u_trace_iteration, &
                  a11_stage_name(case_name), ia-1, this%lattice%ib(ia)-1, 2, -1, -1, -1, -1, -1, &
                  shift_index, shift_value, ev_raw_total, ev_case_total, ev_case_total-ev_raw_total, 0.0_rp)
               do channel = 1, 2
                  call this%write_fixed_u_csv_row(out_unit, 'm0_a11_case', this%fixed_u_trace_iteration, &
                     a11_stage_name(case_name), ia-1, this%lattice%ib(ia)-1, 2, channel, -1, -1, -1, -1, &
                     shift_index, shift_value, baseline_sumev(channel), sumev_case(channel), &
                     this%fixed_u_m0_q0(ia, channel), this%fixed_u_m0_ql(ia, 3, 3, channel))
                  call this%write_fixed_u_csv_row(out_unit, 'm0_a11_center', this%fixed_u_trace_iteration, &
                     a11_stage_name(case_name), ia-1, this%lattice%ib(ia)-1, 2, channel, -1, -1, -1, -1, &
                     shift_index, shift_value, this%fixed_u_m0_pl_diag(ia, 3, channel, 1), &
                     this%fixed_u_m0_pl_diag(ia, 3, channel, 1)-shift_value, &
                     pl_raw(3, channel), pl_case(3, channel))
               end do
               call this%write_fixed_u_csv_row(out_unit, 'm0_a11_case', this%fixed_u_trace_iteration, &
                  a11_stage_name(case_name), ia-1, this%lattice%ib(ia)-1, 2, -1, -1, -1, -1, -1, &
                  shift_index, shift_value, sum(baseline_sumev), sum(sumev_case), &
                  sum(this%fixed_u_m0_q0(ia, :)), sum(this%fixed_u_m0_ql(ia, 3, 3, :)))
               call this%write_fixed_u_csv_row(out_unit, 'm0_a11_complete', this%fixed_u_trace_iteration, &
                  a11_stage_name(case_name), ia-1, this%lattice%ib(ia)-1, 2, -1, -1, -1, -1, -1, &
                  shift_index, shift_value, real(nval, rp), 1.0_rp, 0.0_rp, 0.0_rp)
            end if
            ! NEWRHO returns rho/SUMEC/SUMEV, but the full ATOMSC ETOT
            ! components are recomputed independently below from copied arrays.
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_radial_charge', &
               this%fixed_u_trace_iteration, trim(case_name), ia - 1, this%lattice%ib(ia)-1, &
               -1, 1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, charge(1), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_radial_charge', &
               this%fixed_u_trace_iteration, trim(case_name), ia - 1, this%lattice%ib(ia)-1, &
               -1, 2, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, charge(2), 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_radial_charge', &
               this%fixed_u_trace_iteration, trim(case_name), ia - 1, this%lattice%ib(ia)-1, &
               -1, 0, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, charge(1)+charge(2), &
               charge(1)-charge(2), 0.0_rp)
            do ir = 1, nr
               do channel = 1, 2
                  call this%write_fixed_u_csv_row(out_unit, 'm0_replay_rho', this%fixed_u_trace_iteration, &
                     trim(case_name), ia - 1, this%lattice%ib(ia) - 1, -1, channel, -1, -1, ir - 1, -1, 0, &
                     0.0_rp, 0.0_rp, rho_case(ir, channel), 0.0_rp, 0.0_rp)
               end do
            end do
            close(out_unit)
         end if
         deallocate(replay_state%vtn, replay_state%vzt, replay_state%fun2, replay_state%phi_amp)
         deallocate(rho_case, sumec_case, sumev_case, ec_case, ev_case)
      end subroutine run_replay_case

      subroutine evaluate_frozen_components(this, replay_state, replay_xc, atom, z, a, b, nr, rofi_base, &
                                            v_base, vhrmax, b_fsm, rho_sample, sumec_case, sumev_case, &
                                            out_unit, site, component_stage, components)
         class(self), intent(inout) :: this, replay_state
         type(xc), intent(in) :: replay_xc
         class(symbolic_atom), intent(in) :: atom
         real(rp), intent(in) :: z, a, b, rofi_base(nr), v_base(nr,2), vhrmax, b_fsm
         integer, intent(in) :: nr, out_unit, site
         character(len=*), intent(in) :: component_stage
         real(rp), intent(in) :: rho_sample(nr,2), sumec_case(2), sumev_case(2)
         real(rp), intent(out) :: components(10)
         integer :: isp, comp
         real(rp) :: rhoeps_local, rhomu_local, rhovh_local, utot_local, ekin_local, etot_local
         real(rp) :: vnucl_local, zvnucl_local, vsum_local, value_local
         real(rp), allocatable :: rho_work(:,:), v_work(:,:), rofi_work(:)
         real(rp) :: rvh_local(2), rho0_local(2), reps_local(2), rmu_local(2)

         allocate(rho_work(nr,2), v_work(nr,2), rofi_work(nr))
         rho_work = rho_sample
         v_work = v_base
         rofi_work = rofi_base
         rvh_local = 0.0_rp
         rho0_local = 0.0_rp
         reps_local = 0.0_rp
         rmu_local = 0.0_rp
         vsum_local = 0.0_rp
         call POISS0(z, a, b, rofi_work, rho_work, nr, vhrmax, v_work, rvh_local, vsum_local, 2)
         vnucl_local = v_work(1,1)
         zvnucl_local = -z*vnucl_local
         call replay_state%VXC0SP(replay_xc, z, a, b, rofi_work, rho_sample, nr, v_work, &
                                  rho0_local, reps_local, rmu_local, 2, b_fsm)
         rhovh_local = sum(rvh_local)
         rhoeps_local = sum(reps_local)
         rhomu_local = sum(rmu_local)
         utot_local = 0.5_rp*(rhovh_local + zvnucl_local)
         ekin_local = sum(sumev_case) + sum(sumec_case) - rhovh_local - rhomu_local
         etot_local = ekin_local + utot_local + rhoeps_local
         components = [sum(sumev_case), sum(sumec_case), rhovh_local, rhomu_local, vnucl_local, &
                       zvnucl_local, ekin_local, utot_local, rhoeps_local, etot_local]
         do comp = 0, 9
            value_local = components(comp+1)
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_energy_component', &
               this%fixed_u_trace_iteration, trim(component_stage), site-1, this%lattice%ib(site)-1, &
               2, -1, -1, -1, -1, -1, comp, 0.0_rp, 0.0_rp, value_local, 0.0_rp, 0.0_rp)
         end do
         do isp = 1, 2
            call this%write_fixed_u_csv_row(out_unit, 'm0_replay_energy_component', &
               this%fixed_u_trace_iteration, trim(component_stage), site-1, this%lattice%ib(site)-1, &
               2, isp, -1, -1, -1, -1, 10, 0.0_rp, 0.0_rp, rvh_local(isp), 0.0_rp, 0.0_rp)
         end do
         deallocate(rho_work, v_work, rofi_work)
      end subroutine evaluate_frozen_components
   end subroutine replay_fixed_u_m0_moments

   subroutine update_fixed_u_residuals(this)
      class(self), intent(inout) :: this
      complex(rp), allocatable :: density_override(:, :, :, :, :), vout(:, :, :)
      complex(rp), allocatable :: p_in(:, :, :, :), p_out(:, :, :, :), v_in(:, :, :, :), v_out(:, :, :, :)
      integer, allocatable :: site_for_type(:), rec_type(:)
      integer :: nsite, ntype, block, isite, iclu, itype, il, ispin, mmax, i, j, flat, orbital_i, orbital_j, worst_block
      integer :: worst_l, worst_site, unit, ios
      real(rp) :: rp_abs, rp_rel, rv_abs, rv_rel, antiherm_p, antiherm_v, u_value, j_value
      real(rp) :: rp_target, rv_target
      logical :: active, physical

      this%fixed_u_residual_active = .false.
      this%fixed_u_residual_converged = .false.
      this%fixed_u_residual_physical = .false.
      if (.not. this%hamiltonian%hubbard_u_general_check .or. &
          trim(this%hamiltonian%hubbard_u_potential_form) /= 'liechtenstein' .or. &
          this%hamiltonian%hubbard_u_sc_check .or. this%hamiltonian%hubbard_v_check) then
         this%mix%fixed_u_stable_updates = 0
         return
      end if

      nsite = this%lattice%nrec
      ntype = this%lattice%ntype
      if (nsite < 1 .or. ntype < 1) return
      if (.not. allocated(this%lattice%irec)) error stop 'active fixed-U residual: recursion-site map is not allocated'
      if (.not. allocated(this%lattice%iz)) error stop 'active fixed-U residual: cluster type map is not allocated'
      if (.not. allocated(this%lattice%symbolic_atoms)) &
         error stop 'active fixed-U residual: symbolic atom table is not allocated'
      if (size(this%lattice%irec) < nsite) error stop 'active fixed-U residual: recursion-site map is too short'
      if (size(this%lattice%symbolic_atoms) < ntype) &
         error stop 'active fixed-U residual: symbolic atom table is too short'
      allocate(site_for_type(ntype), rec_type(nsite))
      site_for_type = 0
      do isite = 1, nsite
         iclu = this%lattice%irec(isite)
         if (iclu < 1) error stop 'active fixed-U residual: recursion site has no cluster mapping'
         if (iclu > size(this%lattice%iz)) error stop 'active fixed-U residual: cluster mapping is out of range'
         itype = this%lattice%iz(iclu)
         if (itype < 1 .or. itype > ntype) error stop 'active fixed-U residual: mapped atom type is out of range'
         rec_type(isite) = itype
         if (site_for_type(itype) == 0) site_for_type(itype) = isite
      end do

      active = .false.
      do itype = 1, ntype
         do il = 1, min(4, size(this%mix%ldm_new, 2))
            u_value = 0.0_rp
            j_value = 0.0_rp
            if (allocated(this%lattice%symbolic_atoms(itype)%potential%hubbard_u)) then
               if (il <= size(this%lattice%symbolic_atoms(itype)%potential%hubbard_u)) &
                  u_value = this%lattice%symbolic_atoms(itype)%potential%hubbard_u(il)
            end if
            if (allocated(this%lattice%symbolic_atoms(itype)%potential%hubbard_j)) then
               if (il <= size(this%lattice%symbolic_atoms(itype)%potential%hubbard_j)) &
                  j_value = this%lattice%symbolic_atoms(itype)%potential%hubbard_j(il)
            end if
            if (abs(u_value) > 1.0e-10_rp .or. abs(j_value) > 1.0e-10_rp) active = .true.
         end do
      end do
      if (.not. active) then
         this%mix%fixed_u_stable_updates = 0
         return
      end if
      this%fixed_u_residual_active = .true.
      if (.not. allocated(this%hamiltonian%hubbard_u_pot)) &
         error stop 'active fixed-U residual: Hubbard potential is not allocated'

      allocate(density_override(ntype, 4, 2, 7, 7), vout(nb, nb, ntype))
      allocate(p_in(nsite*4, 2, 7, 7), p_out(nsite*4, 2, 7, 7), &
               v_in(nsite*4, 2, 7, 7), v_out(nsite*4, 2, 7, 7))
      density_override = 0.0_rp
      vout = 0.0_rp
      p_in = 0.0_rp
      p_out = 0.0_rp
      v_in = 0.0_rp
      v_out = 0.0_rp

      ! The production Hubbard operator is type-local. Use the same first
      ! site-of-type density convention for its read-only Pout evaluation;
      ! residuals below remain site-resolved using each site's own P matrices.
      do itype = 1, ntype
         isite = site_for_type(itype)
         if (isite == 0) cycle
         ! Only traverse angular-momentum channels represented by this build's
         ! allocated LDM.  The fixed four-shell potential storage can exceed
         ! the active basis (e.g. a d-only build has ldm_new(:,1:3,...)).
         do il = 1, min(4, size(this%mix%ldm_new, 2))
            mmax = 2*il - 1
            do ispin = 1, 2
               do i = 1, mmax
                  do j = 1, mmax
                     flat = (i - 1)*mmax + j
                     density_override(itype, il, ispin, i, j) = this%mix%ldm_new(isite, il, ispin, flat)
                  end do
               end do
            end do
         end do
      end do
      call this%hamiltonian%calculate_hubbard_u_potential_general( &
         density_override=density_override, potential_out=vout)

      do isite = 1, nsite
         itype = rec_type(isite)
         do il = 1, min(4, size(this%mix%ldm_new, 2))
            u_value = 0.0_rp
            j_value = 0.0_rp
            if (allocated(this%lattice%symbolic_atoms(itype)%potential%hubbard_u)) then
               if (il <= size(this%lattice%symbolic_atoms(itype)%potential%hubbard_u)) &
                  u_value = this%lattice%symbolic_atoms(itype)%potential%hubbard_u(il)
            end if
            if (allocated(this%lattice%symbolic_atoms(itype)%potential%hubbard_j)) then
               if (il <= size(this%lattice%symbolic_atoms(itype)%potential%hubbard_j)) &
                  j_value = this%lattice%symbolic_atoms(itype)%potential%hubbard_j(il)
            end if
            if (abs(u_value) <= 1.0e-10_rp .and. abs(j_value) <= 1.0e-10_rp) cycle
            block = (il - 1)*nsite + isite
            mmax = 2*il - 1
            do ispin = 1, 2
               do i = 1, mmax
                  do j = 1, mmax
                     flat = (i - 1)*mmax + j
                     p_in(block, ispin, i, j) = this%mix%ldm_old(isite, il, ispin, flat)
                     p_out(block, ispin, i, j) = this%mix%ldm_new(isite, il, ispin, flat)
                     orbital_i = (il - 1)**2 + i + (ispin - 1)*spin_off
                     orbital_j = (il - 1)**2 + j + (ispin - 1)*spin_off
                     v_in(block, ispin, i, j) = this%hamiltonian%hubbard_u_pot(orbital_i, orbital_j, itype)
                     v_out(block, ispin, i, j) = vout(orbital_i, orbital_j, itype)
                  end do
               end do
            end do
         end do
      end do

      call fixed_u_residual_metrics(p_in, p_out, v_in, v_out, rp_abs, rp_rel, rv_abs, rv_rel, &
                                    physical, antiherm_p, antiherm_v, worst_block)
      if (worst_block > 0) then
         worst_site = mod(worst_block - 1, nsite) + 1
         worst_l = (worst_block - 1)/nsite
      else
         worst_site = 0
         worst_l = -1
      end if
      this%fixed_u_rp_abs = rp_abs
      this%fixed_u_rp_rel = rp_rel
      this%fixed_u_rv_abs = rv_abs
      this%fixed_u_rv_rel = rv_rel
      this%fixed_u_antiherm_p = antiherm_p
      this%fixed_u_antiherm_v = antiherm_v
      this%fixed_u_worst_site = worst_site
      this%fixed_u_worst_l = worst_l
      this%fixed_u_residual_physical = physical
      call fixed_u_r05_targets(this%use_kspace, this%r05_p_tolerance, this%r05_v_tolerance, &
         this%hamiltonian%fixed_u_test_tighten_r05, rp_target, rv_target)
      call this%mix%advance_fixed_u_convergence(rp_abs, rv_abs, physical, this%fixed_u_residual_converged, &
                                                rp_target, rv_target)
      if (rank == 0) then
         call g_logger%info('FIXED_U_RESIDUAL rp_abs='//fmt('es12.4',rp_abs)//' rp_rel='//fmt('es12.4',rp_rel)// &
            ' rv_abs_Ry='//fmt('es12.4',rv_abs)//' rv_rel='//fmt('es12.4',rv_rel)// &
            ' target_rp='//fmt('es12.4',rp_target)//' target_rv_Ry='//fmt('es12.4',rv_target)// &
            ' required_streak='//int2str(fixed_u_required_stable_updates)// &
            ' antiherm_P='//fmt('es12.4',antiherm_p)//' antiherm_V='//fmt('es12.4',antiherm_v)// &
            ' physical='//trim(log2str(physical))//' streak='//int2str(this%mix%fixed_u_stable_updates)// &
            ' worst_site='//int2str(worst_site)//' l='//int2str(worst_l), __FILE__, __LINE__)
      end if
      if (this%fixed_u_trace_enabled .and. rank == 0) then
         open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
              action='write', iostat=ios)
         if (ios == 0) then
            call this%write_fixed_u_csv_row(unit, 'fixed_u_rp_abs', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, rp_abs, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'fixed_u_rp_rel', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, rp_rel, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'fixed_u_rv_abs_ry', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, rv_abs, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'fixed_u_rv_rel', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, rv_rel, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'fixed_u_rp_target', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, rp_target, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'fixed_u_rv_target_ry', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, rv_target, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'fixed_u_antiherm_p', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, antiherm_p, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'fixed_u_antiherm_v', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, antiherm_v, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(unit, 'fixed_u_physical', this%fixed_u_trace_iteration, 'SCF_CONVERGENCE', &
               worst_site - 1, -1, worst_l, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, merge(1.0_rp,0.0_rp,physical), &
               0.0_rp, 0.0_rp)
            close(unit)
         end if
      end if
      deallocate(site_for_type, rec_type, density_override, vout, p_in, p_out, v_in, v_out)
   end subroutine update_fixed_u_residuals

   subroutine write_fixed_u_total_energy(this)
      class(self), intent(in) :: this
      integer :: unit, ios
      real(rp) :: value

      if (.not. this%fixed_u_trace_enabled .or. rank /= 0) return
      if (this%fixed_u_trace_selected_iteration > 0 .and. &
          this%fixed_u_trace_iteration /= this%fixed_u_trace_selected_iteration) return
      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios /= 0) return
      value = sum(this%symbolic_atom(:)%potential%sumev)
      call total_row('sum_ev', value)
      value = sum(this%symbolic_atom(:)%potential%sumec)
      call total_row('sum_ec', value)
      value = sum(this%symbolic_atom(:)%potential%ekin)
      call total_row('ekin', value)
      value = sum(this%symbolic_atom(:)%potential%utot)
      call total_row('utot', value)
      value = sum(this%symbolic_atom(:)%potential%rhoeps)
      call total_row('rhoeps', value)
      value = sum(this%symbolic_atom(:)%potential%etot)
      call total_row('etot_total', value)
      close(unit)
   contains
      subroutine total_row(label, energy_value)
         character(len=*), intent(in) :: label
         real(rp), intent(in) :: energy_value
         call this%write_fixed_u_csv_row(unit, label, this%fixed_u_trace_iteration, 'POST_SCF_RADIAL', &
            -1, -1, -1, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, energy_value, 0.0_rp, 0.0_rp)
      end subroutine total_row
   end subroutine write_fixed_u_total_energy

   !> Report provisional fixed-U energy diagnostics without modifying the
   !> physical total-energy field.  P_out is the unmixed density from this
   !> iteration; V_in is the operator actually held by the Hamiltonian.
   subroutine write_hubbard_moment_energy_diagnostics(this)
      class(self), intent(inout) :: this
      complex(rp), allocatable :: density(:, :, :, :, :), potential(:, :, :, :, :)
      complex(rp) :: p_d(5, 5), v_d(5, 5)
      real(rp), allocatable :: e_int(:, :), e_dc(:, :), e_hub(:, :), tr_vp(:, :)
      real(rp) :: site_int, site_dc, site_hub, site_t, site_external, direct_t, etot_site
      real(rp) :: total_int, total_dc, total_hub, total_t, total_external, total_t_actual, value_a, value_b, value_c
      real(rp) :: entropy, ts, value_a_f, value_b_f, value_c_f
      real(rp) :: test_energy, test_f_local, test_f_spectral
      real(rp) :: production_hub, base_physical_energy
      real(rp), allocatable :: local_weights(:)
      complex(rp) :: test_operator(5, 5)
      integer :: ik_global
      logical :: entropy_available
      integer :: ia, plusbulk, itype, spin, i, j, flat, status, orbital_i, orbital_j
      character(len=64) :: candidate_a, candidate_b, candidate_c
      character(len=64) :: entropy_text, ts_text, candidate_a_f_text, candidate_b_f_text, candidate_c_f_text
      character(len=64) :: test_energy_text, test_f_local_text, test_f_spectral_text

      if (.not. (this%hamiltonian%hubbard_u_general_check .or. this%hamiltonian%hubbard_u_impurity_check .or. &
                 this%hamiltonian%hubbard_test_operator_check)) return
      entropy = 0.0_rp
      ts = 0.0_rp
      test_energy = 0.0_rp
      entropy_available = .false.
      if (this%use_kspace .and. allocated(this%reciprocal_scf_cache)) then
         if (allocated(this%reciprocal_scf_cache%eigenvalues) .and. &
             allocated(this%reciprocal_scf_cache%k_weights)) then
            allocate(local_weights(size(this%reciprocal_scf_cache%eigenvalues, 2)))
            do i = 1, size(local_weights)
               ik_global = i
               if (allocated(this%reciprocal_scf_cache%k_l2g_map)) then
                  if (i <= size(this%reciprocal_scf_cache%k_l2g_map)) &
                     ik_global = this%reciprocal_scf_cache%k_l2g_map(i)
               end if
               if (ik_global < 1 .or. ik_global > size(this%reciprocal_scf_cache%k_weights)) &
                  call g_logger%fatal('Mermin entropy: invalid local/global k-point map.', __FILE__, __LINE__)
               local_weights(i) = this%reciprocal_scf_cache%k_weights(ik_global)
            end do
            entropy = fermi_dirac_entropy(this%reciprocal_scf_cache%eigenvalues, local_weights, &
               this%reciprocal_scf_cache%fermi_level, &
               max(this%reciprocal_scf_cache%temperature*kB_Ry_per_K, 1.0e-10_rp))
#ifdef USE_MPI
            if (this%reciprocal_scf_cache%k_mesh_distributed_active) &
               call MPI_ALLREDUCE(MPI_IN_PLACE, entropy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
            ts = this%reciprocal_scf_cache%temperature*kB_Ry_per_K*entropy
            entropy_available = .true.
            deallocate(local_weights)
         end if
      end if
      if (rank /= 0) return
      production_hub = this%fixed_hubbard_energy_pout()
      base_physical_energy = this%physical_total_energy - production_hub
      test_operator = cmplx(0.0_rp, 0.0_rp, rp)
      if (this%hamiltonian%hubbard_test_operator_check) &
         call this%hamiltonian%get_hubbard_test_operator(test_operator)
      allocate(density(this%lattice%ntype, 4, 2, 7, 7), potential(this%lattice%ntype, 4, 2, 7, 7))
      allocate(e_int(this%lattice%ntype, 4), e_dc(this%lattice%ntype, 4), &
               e_hub(this%lattice%ntype, 4), tr_vp(this%lattice%ntype, 4))
      density = cmplx(0.0_rp, 0.0_rp, rp)
      potential = cmplx(0.0_rp, 0.0_rp, rp)
      total_int = 0.0_rp
      total_dc = 0.0_rp
      total_hub = 0.0_rp
      total_t = 0.0_rp
      total_external = 0.0_rp
      total_t_actual = 0.0_rp

      do ia = 1, this%lattice%nrec
         plusbulk = this%lattice%nbulk + ia
         ! Recursive indices address cluster sites, while the energy
         ! evaluator arrays are keyed by symbolic-atom type. In an impurity
         ! cluster ib(ia) is a geometric representative index, not that type.
         itype = plusbulk
         density = cmplx(0.0_rp, 0.0_rp, rp)
         potential = cmplx(0.0_rp, 0.0_rp, rp)
         p_d = cmplx(0.0_rp, 0.0_rp, rp)
         v_d = cmplx(0.0_rp, 0.0_rp, rp)
         direct_t = 0.0_rp
         site_external = 0.0_rp
         do spin = 1, 2
            do i = 1, 5
               orbital_i = 4 + i + spin_off*(spin - 1)
               do j = 1, 5
                  orbital_j = 4 + j + spin_off*(spin - 1)
                  flat = (i - 1)*5 + j
                  p_d(i, j) = this%mix%ldm_new(ia, 3, spin, flat)
                  v_d(i, j) = this%hamiltonian%hubbard_u_pot(orbital_i, orbital_j, itype)
                  density(itype, 3, spin, i, j) = p_d(i, j)
                  potential(itype, 3, spin, i, j) = v_d(i, j)
                  if (this%hamiltonian%hubbard_test_operator_check .and. &
                      (this%hamiltonian%hubbard_test_spin == 0 .or. this%hamiltonian%hubbard_test_spin == spin)) &
                     potential(itype, 3, spin, i, j) = v_d(i, j) - &
                        this%hamiltonian%hubbard_test_lambda*test_operator(i, j)
               end do
            end do
            do i = 1, 5
               do j = 1, 5
                  direct_t = direct_t + real(v_d(i, j)*p_d(j, i), rp)
                  if (this%hamiltonian%hubbard_test_operator_check .and. &
                      (this%hamiltonian%hubbard_test_spin == 0 .or. this%hamiltonian%hubbard_test_spin == spin)) then
                     site_external = site_external + this%hamiltonian%hubbard_test_lambda * &
                        real(test_operator(i, j)*p_d(j, i), rp)
                  end if
               end do
            end do
         end do
         call this%hamiltonian%evaluate_hubbard_energy(density, potential, e_int, e_dc, e_hub, tr_vp)
         site_int = sum(e_int(itype, :))
         site_dc = sum(e_dc(itype, :))
         site_hub = sum(e_hub(itype, :))
         site_t = sum(tr_vp(itype, :))
         etot_site = this%symbolic_atom(plusbulk)%potential%etot
         value_a = 0.0_rp
         value_b = 0.0_rp
         value_c = 0.0_rp
         candidate_a = 'UNAVAILABLE'
         candidate_b = 'UNAVAILABLE'
         candidate_c = 'UNAVAILABLE'
         if (trim(this%hamiltonian%hubbard_radial_moments) == 'spectral') then
            value_a = etot_site - site_t + site_hub
            value_c = etot_site + site_hub
            candidate_a = trim(fmt('e24.16', value_a))
            candidate_c = trim(fmt('e24.16', value_c))
         else
            value_b = etot_site + site_hub + site_external
            candidate_b = trim(fmt('e24.16', value_b))
         end if
         call g_logger%info('HUBBARD_MOMENT_ENERGY_SITE site='//fmt('i0', ia)//' type='//fmt('i0', itype)// &
            ' E_int='//fmt('e24.16', site_int)//' E_DC='//fmt('e24.16', site_dc)// &
            ' E_Hub='//fmt('e24.16', site_hub)//' t_in_out='//fmt('e24.16', site_t)// &
            ' t_U='//fmt('e24.16', site_t)//' e_external='//fmt('e24.16', site_external)// &
            ' t_total_actual_Vin='//fmt('e24.16', direct_t)// &
            ' direct_trace_check='//fmt('e24.16', direct_t-site_t-site_external)// &
            ' ETOT_selected='//fmt('e24.16', etot_site)//' candidate_a='//trim(candidate_a)// &
            ' candidate_b='//trim(candidate_b)//' candidate_c='//trim(candidate_c)// &
            ' candidate_F=UNAVAILABLE_GLOBAL_TS_NOT_SITE_DECOMPOSABLE'// &
            ' candidate_status=PROVISIONAL_ONLY', __FILE__, __LINE__)
         total_int = total_int + site_int
         total_dc = total_dc + site_dc
         total_hub = total_hub + site_hub
         total_t = total_t + site_t
         total_external = total_external + site_external
         total_t_actual = total_t_actual + direct_t
         test_energy = test_energy + site_external
      end do

      candidate_a = 'UNAVAILABLE'
      candidate_b = 'UNAVAILABLE'
      candidate_c = 'UNAVAILABLE'
      value_a_f = 0.0_rp
      value_b_f = 0.0_rp
      value_c_f = 0.0_rp
      if (trim(this%hamiltonian%hubbard_radial_moments) == 'spectral') then
         value_a = this%physical_total_energy - total_t + total_hub
         value_c = this%physical_total_energy + total_hub
         candidate_a = trim(fmt('e24.16', value_a))
         candidate_c = trim(fmt('e24.16', value_c))
         value_a_f = value_a - ts
         value_c_f = value_c - ts
      else
         if (this%is_eligible_local_hubbard() .and. &
             trim(this%hamiltonian%hubbard_radial_moments) == 'local') then
            value_b = this%physical_total_energy + total_external
         else
            value_b = this%physical_total_energy + total_hub + total_external
         end if
         candidate_b = trim(fmt('e24.16', value_b))
         value_b_f = value_b - ts
      end if
      entropy_text = 'UNAVAILABLE_KS_ONLY'
      ts_text = 'UNAVAILABLE_KS_ONLY'
      candidate_a_f_text = 'UNAVAILABLE'
      candidate_b_f_text = 'UNAVAILABLE'
      candidate_c_f_text = 'UNAVAILABLE'
      test_energy_text = 'UNAVAILABLE'
      test_f_local_text = 'UNAVAILABLE'
      test_f_spectral_text = 'UNAVAILABLE'
      if (entropy_available) then
         entropy_text = trim(fmt('e24.16', entropy))
         ts_text = trim(fmt('e24.16', ts))
         if (candidate_a /= 'UNAVAILABLE') candidate_a_f_text = trim(fmt('e24.16', value_a_f))
         if (candidate_b /= 'UNAVAILABLE') candidate_b_f_text = trim(fmt('e24.16', value_b_f))
         if (candidate_c /= 'UNAVAILABLE') candidate_c_f_text = trim(fmt('e24.16', value_c_f))
         if (this%hamiltonian%hubbard_test_operator_check) then
            test_energy_text = trim(fmt('e24.16', test_energy))
            test_f_local = this%physical_total_energy + test_energy - ts
            test_f_spectral = this%physical_total_energy - ts
            if (trim(this%hamiltonian%hubbard_radial_moments) == 'local') then
               test_f_local_text = trim(fmt('e24.16', test_f_local))
            else
               test_f_spectral_text = trim(fmt('e24.16', test_f_spectral))
            end if
         end if
      end if
      call g_logger%info('HUBBARD_MOMENT_ENERGY_TOTAL E_int='//fmt('e24.16', total_int)// &
         ' E_DC='//fmt('e24.16', total_dc)//' E_Hub='//fmt('e24.16', total_hub)// &
         ' E_Hub_P_OUT='//fmt('e24.16', production_hub)// &
         ' t_in_out='//fmt('e24.16', total_t)//' t_U='//fmt('e24.16', total_t)// &
         ' e_external='//fmt('e24.16', total_external)//' t_total_actual_Vin='//fmt('e24.16', total_t_actual)// &
         ' base_physical_total_energy='//fmt('e24.16', base_physical_energy)// &
         ' physical_total_energy='// &
         fmt('e24.16', this%physical_total_energy)//' candidate_a='//trim(candidate_a)// &
         ' candidate_b='//trim(candidate_b)//' candidate_c='//trim(candidate_c)// &
         ' entropy='//trim(entropy_text)//' TS='//trim(ts_text)// &
         ' candidate_a_F='//trim(candidate_a_f_text)//' candidate_b_F='//trim(candidate_b_f_text)// &
         ' candidate_c_F='//trim(candidate_c_f_text)// &
         ' test_lambda_trace_wp='//trim(test_energy_text)//' test_F_local='//trim(test_f_local_text)// &
         ' test_F_spectral='//trim(test_f_spectral_text)// &
         ' candidate_status=PROVISIONAL_ONLY; unavailable candidates require the other radial representation', &
         __FILE__, __LINE__)
      deallocate(density, potential, e_int, e_dc, e_hub, tr_vp)
   end subroutine write_hubbard_moment_energy_diagnostics

   logical function is_eligible_local_hubbard(this) result(eligible)
      class(self), intent(in) :: this

      eligible = this%hamiltonian%hubbard_u_general_check .and. &
         trim(this%hamiltonian%hubbard_u_potential_form) == 'liechtenstein' .and. &
         this%control%nsp == 1 .and. this%control%lmax == 2 .and. &
         this%control%calctype == 'B' .and. &
         this%lattice%nbulk == 0 .and. this%lattice%nrec == this%lattice%ntype .and. &
         .not. this%control%constraints_enable .and. &
         .not. this%control%has_soc() .and. .not. this%control%is_noncollinear() .and. &
         .not. this%hamiltonian%hoh .and. .not. this%hamiltonian%hubbard_u_sc_check .and. &
         .not. this%hamiltonian%hubbard_v_check .and. .not. this%hamiltonian%hubbard_u_impurity_check .and. &
         (.not. this%hamiltonian%hubbard_test_operator_check .or. &
          this%hamiltonian%hubbard_test_lambda == 0.0_rp)
   end function is_eligible_local_hubbard

   logical function is_eligible_local_moments(this, reason) result(eligible)
      class(self), intent(in) :: this
      character(len=*), intent(out), optional :: reason
      integer :: ia, ntype
      logical :: scope_eligible, active
      logical, allocatable :: host_active(:)
      character(len=16), allocatable :: host_mode(:)
      character(len=128) :: local_reason

      if (this%control%calctype == 'B') then
         eligible = this%is_eligible_local_hubbard()
         local_reason = 'existing fixed-U clean-bulk eligibility'
         if (.not. eligible) local_reason = 'bulk is outside existing fixed-U clean-bulk eligibility'
         if (present(reason)) reason = trim(local_reason)
         return
      end if

      scope_eligible = (this%hamiltonian%hubbard_u_general_check .or. &
         this%hamiltonian%hubbard_u_impurity_check) .and. &
         trim(this%hamiltonian%hubbard_u_potential_form) == 'liechtenstein' .and. &
         this%control%nsp == 1 .and. this%control%lmax == 2 .and. &
         (this%control%calctype == 'B' .or. this%control%calctype == 'I' .or. &
          this%control%calctype == 'S') .and. &
         .not. this%control%constraints_enable .and. &
         .not. this%control%has_soc() .and. .not. this%control%is_noncollinear() .and. &
         .not. this%hamiltonian%hoh .and. .not. this%hamiltonian%hubbard_u_sc_check .and. &
         .not. this%hamiltonian%hubbard_v_check .and. &
         .not. this%hamiltonian%hubbard_test_operator_check

      ntype = size(this%symbolic_atom)
      allocate(host_active(ntype), host_mode(ntype))
      host_active = .false.
      host_mode = 'spectral'
      do ia = 1, ntype
         if (.not. this%is_frozen_hubbard_host(ia)) cycle
         if (.not. allocated(this%symbolic_atom(ia)%potential%hubbard_u) .or. &
             .not. allocated(this%symbolic_atom(ia)%potential%hubbard_j)) cycle
         active = any(abs(this%symbolic_atom(ia)%potential%hubbard_u) > 1.0e-10_rp) .or. &
                  any(abs(this%symbolic_atom(ia)%potential%hubbard_j) > 1.0e-10_rp)
         if (.not. active) cycle
         host_active(ia) = .true.
         if (.not. this%symbolic_atom(ia)%potential%hubbard_file_has_uj) cycle
         if (len_trim(this%symbolic_atom(ia)%potential%hubbard_radial_moments) > 0) &
            host_mode(ia) = trim(this%symbolic_atom(ia)%potential%hubbard_radial_moments)
      end do
      call auto_impurity_local_host_eligible(scope_eligible, host_active, host_mode, eligible, local_reason)
      if (present(reason)) reason = trim(local_reason)
   end function is_eligible_local_moments

   logical function is_frozen_hubbard_host(this, atom_type) result(frozen)
      class(self), intent(in) :: this
      integer, intent(in) :: atom_type
      integer :: isite, cluster_site, mapped_type

      frozen = .false.
      if (this%control%calctype /= 'I' .and. this%control%calctype /= 'S') return
      if (.not. allocated(this%lattice%irec) .or. .not. allocated(this%lattice%ib) .or. &
          .not. allocated(this%lattice%iz)) then
         call g_logger%fatal('Cannot classify fixed-U provenance without the SCF host/recursion ownership maps.', &
            __FILE__, __LINE__)
         return
      end if
      if (this%lattice%nrec > size(this%lattice%irec) .or. this%lattice%nbulk > size(this%lattice%ib)) then
         call g_logger%fatal('Fixed-U SCF ownership map is shorter than lattice%nrec/nbulk.', __FILE__, __LINE__)
         return
      end if
      ! The recursion ownership map is the same map consumed by the SCF
      ! update: recursion index -> cluster site -> symbolic atom type.
      do isite = 1, min(this%lattice%nrec, size(this%lattice%irec))
         cluster_site = this%lattice%irec(isite)
         if (cluster_site < 1 .or. cluster_site > size(this%lattice%iz)) then
            call g_logger%fatal('Fixed-U recursion ownership map has an invalid cluster-site index.', __FILE__, __LINE__)
            return
         end if
         mapped_type = this%lattice%iz(cluster_site)
         if (mapped_type < 1 .or. mapped_type > size(this%symbolic_atom)) then
            call g_logger%fatal('Fixed-U recursion ownership map has an invalid symbolic atom type.', __FILE__, __LINE__)
            return
         end if
         if (mapped_type == atom_type) return
      end do
      ! A type is a frozen host only if the lattice host ownership table
      ! contains it. Unreferenced symbolic types are not inferred as hosts.
      do isite = 1, min(this%lattice%nbulk, size(this%lattice%ib))
         cluster_site = this%lattice%ib(isite)
         if (cluster_site < 1 .or. cluster_site > size(this%lattice%iz)) then
            call g_logger%fatal('Fixed-U frozen-host ownership map has an invalid cluster-site index.', __FILE__, __LINE__)
            return
         end if
         mapped_type = this%lattice%iz(cluster_site)
         if (mapped_type < 1 .or. mapped_type > size(this%symbolic_atom)) then
            call g_logger%fatal('Fixed-U frozen-host ownership map has an invalid symbolic atom type.', __FILE__, __LINE__)
            return
         end if
         if (mapped_type == atom_type) then
            frozen = .true.
            return
         end if
      end do
   end function is_frozen_hubbard_host

   logical function all_frozen_hubbard_hosts_local(this) result(all_local)
      class(self), intent(in) :: this
      integer :: ia
      logical :: active, found

      all_local = .false.
      found = .false.
      do ia = 1, size(this%symbolic_atom)
         if (.not. this%is_frozen_hubbard_host(ia)) cycle
         if (.not. allocated(this%symbolic_atom(ia)%potential%hubbard_u) .or. &
             .not. allocated(this%symbolic_atom(ia)%potential%hubbard_j)) cycle
         active = any(abs(this%symbolic_atom(ia)%potential%hubbard_u) > 1.0e-10_rp) .or. &
                  any(abs(this%symbolic_atom(ia)%potential%hubbard_j) > 1.0e-10_rp)
         if (.not. active) cycle
         found = .true.
         if (.not. this%symbolic_atom(ia)%potential%hubbard_file_has_uj .or. &
             trim(this%symbolic_atom(ia)%potential%hubbard_radial_moments) /= 'local') return
      end do
      all_local = found
   end function all_frozen_hubbard_hosts_local

   subroutine validate_hubbard_mode_provenance(this)
      class(self), intent(inout) :: this
      integer :: ia
      logical :: active
      logical :: frozen
      character(len=16) :: stored_mode, backend

      backend = 'rs'
      if (this%use_kspace) backend = 'ks'
      do ia = 1, size(this%symbolic_atom)
         if (.not. allocated(this%symbolic_atom(ia)%potential%hubbard_u) .or. &
             .not. allocated(this%symbolic_atom(ia)%potential%hubbard_j)) cycle
         active = any(abs(this%symbolic_atom(ia)%potential%hubbard_u) > 1.0e-10_rp) .or. &
                  any(abs(this%symbolic_atom(ia)%potential%hubbard_j) > 1.0e-10_rp)
         if (.not. active) cycle

         frozen = this%is_frozen_hubbard_host(ia)

         if (frozen) then
            stored_mode = trim(this%symbolic_atom(ia)%potential%hubbard_radial_moments)
            if (len_trim(stored_mode) == 0) then
               stored_mode = 'spectral'
               if (rank == 0) call g_logger%warning('Legacy U-carrying potential for symbolic atom '//fmt('i0', ia)// &
                  ' has no hubbard_radial_moments field; its mode is unknown, treating it as spectral; '// &
                  'restart the host once with the current version to record its mode.', __FILE__, __LINE__)
            end if
            if (this%symbolic_atom(ia)%potential%hubbard_file_has_uj .and. &
                trim(stored_mode) /= trim(this%hamiltonian%hubbard_radial_moments)) then
               call g_logger%fatal('Fixed-U potential moment-mode mismatch for symbolic atom '//fmt('i0', ia)// &
                  ': file mode='//trim(stored_mode)//', run mode='//trim(this%hamiltonian%hubbard_radial_moments)// &
                  '. Regenerate the host potential in the selected mode or select the matching mode.', __FILE__, __LINE__)
            end if
         end if

         if (.not. frozen) then
            stored_mode = trim(this%symbolic_atom(ia)%potential%hubbard_radial_moments)
            if (len_trim(stored_mode) == 0) then
               stored_mode = 'spectral'
               if (rank == 0) call g_logger%info('restart seed in '//trim(stored_mode)//', run uses '// &
                  trim(this%hamiltonian%hubbard_radial_moments)//'; seed will be rewritten', __FILE__, __LINE__)
            else if (trim(stored_mode) /= trim(this%hamiltonian%hubbard_radial_moments) .and. rank == 0) then
               call g_logger%info('restart seed in '//trim(stored_mode)//', run uses '// &
                  trim(this%hamiltonian%hubbard_radial_moments)//'; seed will be rewritten', __FILE__, __LINE__)
            end if
         end if

         if (.not. frozen) then
            this%symbolic_atom(ia)%potential%hubbard_radial_moments = &
               trim(this%hamiltonian%hubbard_radial_moments)
            this%symbolic_atom(ia)%potential%hubbard_backend = trim(backend)
         end if
      end do
   end subroutine validate_hubbard_mode_provenance

   real(rp) function fixed_hubbard_energy_pout(this) result(total_hub)
      class(self), intent(in) :: this
      complex(rp), allocatable :: density(:, :, :, :, :), potential(:, :, :, :, :)
      real(rp), allocatable :: e_int(:, :), e_dc(:, :), e_hub(:, :), tr_vp(:, :)
      complex(rp) :: p_d, v_d
      integer :: isite, itype, spin, i, j, orbital_i, orbital_j

      total_hub = 0.0_rp
      if (.not. this%is_eligible_local_hubbard()) return
      if (trim(this%hamiltonian%hubbard_radial_moments) /= 'local') return

      allocate(density(this%lattice%ntype, 4, 2, 7, 7), potential(this%lattice%ntype, 4, 2, 7, 7))
      allocate(e_int(this%lattice%ntype, 4), e_dc(this%lattice%ntype, 4), &
               e_hub(this%lattice%ntype, 4), tr_vp(this%lattice%ntype, 4))
      do isite = 1, this%lattice%nrec
         itype = this%lattice%ib(isite)
         density = cmplx(0.0_rp, 0.0_rp, rp)
         potential = cmplx(0.0_rp, 0.0_rp, rp)
         do spin = 1, 2
            do i = 1, 5
               orbital_i = 4 + i + spin_off*(spin - 1)
               do j = 1, 5
                  orbital_j = 4 + j + spin_off*(spin - 1)
                  p_d = this%symbolic_atom(this%lattice%nbulk + isite)%potential%ldm(3, spin, i, j)
                  v_d = this%hamiltonian%hubbard_u_pot(orbital_i, orbital_j, itype)
                  density(itype, 3, spin, i, j) = p_d
                  potential(itype, 3, spin, i, j) = v_d
               end do
            end do
         end do
         call this%hamiltonian%evaluate_hubbard_energy(density, potential, e_int, e_dc, e_hub, tr_vp)
         total_hub = total_hub + sum(e_hub(itype, :))
      end do
      deallocate(density, potential, e_int, e_dc, e_hub, tr_vp)
   end function fixed_hubbard_energy_pout

   subroutine write_fixed_u_kspace_trace(this, capture_spectrum)
      class(self), intent(inout) :: this
      logical, intent(in) :: capture_spectrum
      integer :: unit, ios, ik, ig, ib, i, j, n, nsite, nper_spin, base, m, spin, ell
      real(rp) :: kt, occ
      complex(rp) :: z

      if (.not. this%fixed_u_trace_enabled) return
      if (rank /= 0) return
      if (this%fixed_u_trace_ql_only) return
      if (this%fixed_u_trace_moments_only) return
      if (this%fixed_u_trace_selected_iteration > 0 .and. &
          this%fixed_u_trace_iteration /= this%fixed_u_trace_selected_iteration) then
         call g_logger%warning('Fixed-U k-space trace skipped: selected iteration mismatch.', __FILE__, __LINE__)
         return
      end if
      if (.not. allocated(this%reciprocal_scf_cache)) then
         call g_logger%warning('Fixed-U k-space trace skipped: reciprocal cache is not allocated.', __FILE__, __LINE__)
         return
      end if
      if (capture_spectrum .and. (.not. allocated(this%reciprocal_scf_cache%eigenvalues) .or. &
          .not. allocated(this%reciprocal_scf_cache%eigenvectors))) then
         call g_logger%warning('Fixed-U k-space spectrum trace skipped: eigenpairs are not allocated.', __FILE__, __LINE__)
         return
      end if
      open(newunit=unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
           action='write', iostat=ios)
      if (ios /= 0) return
      kt = max(this%reciprocal_scf_cache%temperature*kB_Ry_per_K, 1.0e-10_rp)
      if (allocated(this%reciprocal_scf_cache%hk_bulk)) then
         n = size(this%reciprocal_scf_cache%hk_bulk, 1)
      else if (allocated(this%reciprocal_scf_cache%eigenvectors)) then
         n = size(this%reciprocal_scf_cache%eigenvectors, 1)
      else
         call g_logger%warning('Fixed-U k-space trace skipped: neither Hamiltonian nor eigenvectors are allocated.', &
                               __FILE__, __LINE__)
         close(unit)
         return
      end if
      nsite = this%lattice%nrec
      nper_spin = n/(2*nsite)
      if (capture_spectrum) then
         call this%write_fixed_u_csv_row(unit, 'fermi_level', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
            -1, -1, -1, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, &
            this%reciprocal_scf_cache%fermi_level, 0.0_rp, 0.0_rp)
         call this%write_fixed_u_csv_row(unit, 'temperature_k', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
            -1, -1, -1, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, &
            this%reciprocal_scf_cache%temperature, 0.0_rp, 0.0_rp)
         call this%write_fixed_u_csv_row(unit, 'temperature_ry', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
            -1, -1, -1, -1, -1, -1, -1, -1, 0, 0.0_rp, 0.0_rp, kt, 0.0_rp, 0.0_rp)
         do ik = 1, size(this%reciprocal_scf_cache%eigenvalues, 2)
            ig = k_global_index(ik)
            do ib = 1, size(this%reciprocal_scf_cache%eigenvalues, 1)
               occ = fermi_dirac_occupation(this%reciprocal_scf_cache%eigenvalues(ib, ik), &
                                            this%reciprocal_scf_cache%fermi_level, kt)
               call this%write_fixed_u_csv_row(unit, 'eigenvalue', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
                  -1, -1, -1, -1, ig - 1, ib - 1, -1, -1, 0, 0.0_rp, &
                  this%reciprocal_scf_cache%eigenvalues(ib, ik), 0.0_rp, 0.0_rp, 0.0_rp)
               call this%write_fixed_u_csv_row(unit, 'occupation', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
                  -1, -1, -1, -1, ig - 1, ib - 1, -1, -1, 0, 0.0_rp, 0.0_rp, occ, 0.0_rp, 0.0_rp)
               do i = 1, size(this%reciprocal_scf_cache%eigenvectors, 1)
                  z = this%reciprocal_scf_cache%eigenvectors(i, ib, ik)
                  call this%write_fixed_u_csv_row(unit, 'eigenvector', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
                     -1, -1, -1, -1, ig - 1, ib - 1, i - 1, -1, 0, 0.0_rp, 0.0_rp, 0.0_rp, &
                     real(z, rp), aimag(z))
               end do
            end do
         end do
      else
      do ik = 1, size(this%reciprocal_scf_cache%hk_bulk, 3)
         ig = ik
         if (allocated(this%reciprocal_scf_cache%k_l2g_map)) then
            if (ik <= size(this%reciprocal_scf_cache%k_l2g_map)) ig = this%reciprocal_scf_cache%k_l2g_map(ik)
         end if
         call write_k_header(ik)
         do i = 1, n
            do j = 1, n
               z = this%reciprocal_scf_cache%hk_bulk(i, j, ik)
               call this%write_fixed_u_csv_row(unit, 'h_total', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
                  -1, -1, -1, -1, ig - 1, -1, i - 1, j - 1, 0, 0.0_rp, 0.0_rp, 0.0_rp, real(z, rp), aimag(z))
               if (allocated(this%reciprocal_scf_cache%sk_overlap)) then
                  z = this%reciprocal_scf_cache%sk_overlap(i, j, ik)
                  call this%write_fixed_u_csv_row(unit, 's_overlap', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
                     -1, -1, -1, -1, ig - 1, -1, i - 1, j - 1, 0, 0.0_rp, 0.0_rp, 0.0_rp, &
                     real(z, rp), aimag(z))
               end if
            end do
         end do
      end do
      end if
      ! Basis ordering is site-major, then spin-major, then l,m; for the
      ! supported scalar-collinear spd layout, d orbitals are local indices 5:9.
      if (capture_spectrum) then
         do i = 1, nsite
            do spin = 1, 2
               base = (i - 1)*(2*nper_spin) + (spin - 1)*nper_spin
               do ell = 0, int(sqrt(real(nper_spin, rp))) - 1
                  do m = 0, 2*ell
                     call this%write_fixed_u_csv_row(unit, 'basis_map', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
                        i - 1, this%lattice%iz(i), ell, spin - 1, -1, -1, base + ell*ell + m, m - ell, 0, &
                        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp)
                  end do
               end do
            end do
         end do
      end if
      close(unit)
   contains
      subroutine write_k_header(k_local)
         integer, intent(in) :: k_local
         integer :: kg
         kg = k_global_index(k_local)
         if (allocated(this%reciprocal_scf_cache%k_points) .and. allocated(this%reciprocal_scf_cache%k_weights)) then
            do i = 1, 3
               call krow('k_coordinate', kg, i - 1, -1, 0.0_rp, 0.0_rp, &
                         this%reciprocal_scf_cache%k_points(i, kg), 'fractional')
            end do
            call this%write_fixed_u_csv_row(unit, 'k_weight', this%fixed_u_trace_iteration, 'EIGENPAIRS', &
               -1, -1, -1, -1, kg - 1, -1, -1, -1, 0, &
               this%reciprocal_scf_cache%k_weights(kg), 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp)
         end if
      end subroutine write_k_header
      integer function k_global_index(k_local)
         integer, intent(in) :: k_local
         k_global_index = k_local
         if (allocated(this%reciprocal_scf_cache%k_l2g_map)) then
            if (k_local <= size(this%reciprocal_scf_cache%k_l2g_map)) &
               k_global_index = this%reciprocal_scf_cache%k_l2g_map(k_local)
         end if
      end function k_global_index
      subroutine krow(kind, k_index, row, col, re, im, value, name)
         integer, intent(in) :: k_index, row, col
         character(len=*), intent(in) :: kind, name
         real(rp), intent(in) :: re, im, value
         call this%write_fixed_u_csv_row(unit, trim(kind), this%fixed_u_trace_iteration, 'EIGENPAIRS', &
            -1, -1, -1, -1, k_index - 1, -1, row, col, 0, 0.0_rp, 0.0_rp, value, re, im)
      end subroutine krow
   end subroutine write_fixed_u_kspace_trace

   subroutine write_fixed_u_rs_green(this, unit)
      class(self), intent(in) :: this
      integer, intent(in) :: unit
      integer :: ia, iloc, ie, i, j, npts, nmat, orbital, ell, m, spin
      real(rp), allocatable :: quadrature_weight(:)
      real(rp) :: h, tail
      complex(rp) :: z

      if (this%fixed_u_trace_moments_only) return
      if (.not. allocated(this%green%g0) .or. .not. allocated(this%en%ene)) return
      npts = this%bands%nv1
      if (npts < 3 .or. size(this%en%ene) < npts + 2) return
      nmat = size(this%green%g0, 1)
      allocate(quadrature_weight(size(this%en%ene)))
      quadrature_weight = 0.0_rp
      h = this%en%edel
      do i = 2, npts - 1, 2
         quadrature_weight(i - 1) = quadrature_weight(i - 1) + h/3.0_rp
         quadrature_weight(i) = quadrature_weight(i) + 4.0_rp*h/3.0_rp
         quadrature_weight(i + 1) = quadrature_weight(i + 1) + h/3.0_rp
      end do
      if (this%bands%e1 /= this%en%fermi) then
         tail = (this%en%fermi - this%bands%e1)/6.0_rp
         quadrature_weight(npts) = quadrature_weight(npts) + tail
         quadrature_weight(npts + 1) = quadrature_weight(npts + 1) + 4.0_rp*tail
         quadrature_weight(npts + 2) = quadrature_weight(npts + 2) + tail
      end if

      do ia = start_atom, end_atom
         iloc = g2l_map(ia)
         if (iloc < 1 .or. iloc > size(this%green%g0, 4)) cycle
         do i = 0, nmat - 1
            orbital = mod(i, spin_off)
            spin = i/spin_off
            ell = int(sqrt(real(orbital, rp)))
            m = orbital - ell*ell - ell
            call this%write_fixed_u_csv_row(unit, 'basis_map', this%fixed_u_trace_iteration, 'RS_GREEN', ia - 1, &
               this%lattice%iz(ia), ell, spin, -1, -1, i, m, 0, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp)
         end do
         do ie = 1, npts + 2
            do i = 1, nmat
               do j = 1, nmat
                  z = this%green%g0(i, j, ie, iloc)
                  call this%write_fixed_u_csv_row(unit, 'green_function', this%fixed_u_trace_iteration, 'RS_GREEN', &
                     ia - 1, this%lattice%iz(ia), -1, -1, ie - 1, -1, i - 1, j - 1, 0, &
                     quadrature_weight(ie), this%en%ene(ie), 0.0_rp, real(z, rp), aimag(z))
               end do
            end do
            call this%write_fixed_u_csv_row(unit, 'fermi_step', this%fixed_u_trace_iteration, 'RS_GREEN', &
               ia - 1, this%lattice%iz(ia), -1, -1, ie - 1, -1, -1, -1, 0, quadrature_weight(ie), &
               this%en%ene(ie), merge(1.0_rp, 0.0_rp, this%en%ene(ie) <= this%en%fermi), 0.0_rp, 0.0_rp)
         end do
      end do
      deallocate(quadrature_weight)
   end subroutine write_fixed_u_rs_green

   !> @brief Generate one deterministic scalar-relativistic radial channel
   !>        through the unchanged production RSEQSR/PHDFSR routines.
   !>
   !> This is an audit/test seam only.  The caller supplies a radial mesh and
   !> potential; all numerical work remains in the legacy routines below.
   subroutine legacy_radial_fixture(z, l, a, b, rofi, potential, energy, g, gp, gpp, boundary_slope)
      real(rp), intent(in) :: z, a, b
      integer, intent(in) :: l
      real(rp), intent(in) :: rofi(:), potential(:)
      real(rp), intent(out) :: energy
      real(rp), allocatable, intent(out) :: g(:, :), gp(:, :), gpp(:, :)
      real(rp), intent(in), optional :: boundary_slope
      real(rp) :: val, slo, sum_norm, phi, dphi, phip, dphip, p
      integer :: nr, nre, nn

      nr = size(rofi)
      if (nr < 35 .or. size(potential) /= nr .or. mod(nr, 2) == 0) then
         error stop 'legacy_radial_fixture: use an odd radial mesh with at least 35 points'
      end if
      allocate(g(nr, 2), gp(nr, 2), gpp(nr, 2))
      g = 0.0_rp
      gp = 0.0_rp
      gpp = 0.0_rp
      energy = -0.5_rp
      val = rofi(nr)
      slo = 1.0_rp
      if (present(boundary_slope)) slo = boundary_slope
      nn = 0
      call RSEQSR(-10.0_rp, 10.0_rp, energy, 1.0e-12_rp, z, l, nn, val, slo, potential, g, sum_norm, &
                  a, b, rofi, nr, nre, 0)
      val = val/sqrt(sum_norm)
      slo = slo/sqrt(sum_norm)
      call PHDFSR(z, l, potential, energy, a, b, rofi, nr, g, val, slo, gp, gpp, phi, dphi, phip, dphip, p, &
                  1.0e-12_rp, nn)
   end subroutine legacy_radial_fixture

   !> Invoke the unchanged production NEWRHO routine for a deterministic
   !> valence-only fixture.  This is an audit/test seam; it does not duplicate
   !> or alter the legacy NEWRHO numerical machinery.
   subroutine legacy_newrho_fixture(lmax, z, a, b, rofi, potential, energies, ql, rho, &
                                    sumev_out, sumec_out, hubbard_u_value, hubbard_j_value, pl_override, tol_value)
      integer, intent(in) :: lmax
      real(rp), intent(in) :: z, a, b, rofi(:), potential(:, :), energies(:, :), ql(:, :, :)
      real(rp), intent(out) :: rho(:, :)
      real(rp), intent(out), optional :: sumev_out(:), sumec_out(:)
      real(rp), intent(in), optional :: hubbard_u_value, hubbard_j_value
      real(rp), intent(in), optional :: pl_override(:, :)
      real(rp), intent(in), optional :: tol_value
      type(self) :: audit_state
      type(control), target :: audit_control
      type(symbolic_atom) :: audit_atom
      real(rp), allocatable :: pl(:, :), sumec(:), sumev(:), ec(:), ev(:)
      integer :: nr, nval, l, ispin, ival

      nr = size(rofi)
      if (lmax < 0 .or. lmax > 2 .or. size(potential, 1) /= nr .or. size(potential, 2) /= 2 .or. &
          any(shape(rho) /= [nr, 2]) .or. any(shape(energies) /= [lmax + 1, 2]) .or. &
          any(shape(ql) /= [3, lmax + 1, 2])) then
         error stop 'legacy_newrho_fixture: inconsistent fixture dimensions'
      end if

      audit_control%hyperfine = .false.
      audit_state%control => audit_control
      allocate(audit_state%vtn(nr, 2), audit_state%vzt(nr, 2), audit_state%fun2(nr, lmax + 1, 2), &
               audit_state%phi_amp(nr, lmax + 1, 2))
      audit_state%vtn = 0.0_rp
      audit_state%vzt = 0.0_rp
      audit_state%fun2 = 0.0_rp
      audit_state%phi_amp = 0.0_rp

      call audit_atom%restore_to_default()
      audit_atom%element%atomic_number = z
      audit_atom%element%f_core = 0
      audit_atom%potential%lmax = lmax
      ! Optional U/J metadata supports a frozen-radial test of metadata
      ! independence. NEWRHO itself is unchanged and receives the same V/QL/PL;
      ! this fixture does not inject a Hubbard Hamiltonian or potential.
      if (present(hubbard_u_value) .or. present(hubbard_j_value)) then
         if (allocated(audit_atom%potential%hubbard_u)) deallocate(audit_atom%potential%hubbard_u)
         if (allocated(audit_atom%potential%hubbard_j)) deallocate(audit_atom%potential%hubbard_j)
         allocate(audit_atom%potential%hubbard_u(lmax + 1), audit_atom%potential%hubbard_j(lmax + 1))
         audit_atom%potential%hubbard_u = 0.0_rp
         audit_atom%potential%hubbard_j = 0.0_rp
         if (present(hubbard_u_value)) audit_atom%potential%hubbard_u(lmax + 1) = hubbard_u_value
         if (present(hubbard_j_value)) audit_atom%potential%hubbard_j(lmax + 1) = hubbard_j_value
      end if
      allocate(pl(lmax + 1, 2), sumec(2), sumev(2), ec(10), ev(2*(lmax + 1)))
      do ispin = 1, 2
         do l = 0, lmax
            pl(l + 1, ispin) = real(l + 1, rp) + 0.25_rp
         end do
      end do
      if (present(pl_override)) then
         if (any(shape(pl_override) /= [lmax + 1, 2])) then
            error stop 'legacy_newrho_fixture: PL override has inconsistent dimensions'
         end if
         pl = pl_override
      end if
      ev = 0.0_rp
      ival = 0
      do ispin = 1, 2
         do l = 0, lmax
            ival = ival + 1
            ev(ival) = energies(l + 1, ispin)
         end do
      end do
      if (present(tol_value)) then
         call audit_state%newrho(audit_atom, z, lmax, a, b, nr, rofi, potential, rho, pl, ql, sumec, sumev, &
                                 ec, ev, tol_value, 2, 0)
      else
         call audit_state%newrho(audit_atom, z, lmax, a, b, nr, rofi, potential, rho, pl, ql, sumec, sumev, &
                                 ec, ev, 1.0e-12_rp, 2, 0)
      end if
      nval = 2*(lmax + 1)
      if (size(ev) /= nval) error stop 'legacy_newrho_fixture: invalid valence fixture'
      if (present(sumev_out)) then
         if (size(sumev_out) /= 2) error stop 'legacy_newrho_fixture: SUMEV output must have size 2'
         sumev_out = sumev
      end if
      if (present(sumec_out)) then
         if (size(sumec_out) /= 2) error stop 'legacy_newrho_fixture: SUMEC output must have size 2'
         sumec_out = sumec
      end if
   end subroutine legacy_newrho_fixture

   !> Rehydrate a captured single-site radial state and replay signed d-center
   !> shifts using the actual bands%calculate_pl method and unchanged NEWRHO.
   !> Every shift starts from fresh copies of the exact captured operands.
   subroutine legacy_accepted_pl_newrho_fixture(lmax, z, a, b, tol, f_core, rofi, potential, ec_seed, ev_seed, &
                                                gravity_center, c_in, qpar_in, srdel_in, pl_seed, ql_in, shifts, &
                                                pl_out, rho_out, sumec_out, sumev_out, ec_out, ev_out)
      integer, intent(in) :: lmax, f_core
      real(rp), intent(in) :: z, a, b, tol
      real(rp), intent(in) :: rofi(:), potential(:, :), ec_seed(:), ev_seed(:)
      real(rp), intent(in) :: gravity_center(:, :), c_in(:, :), qpar_in(:, :), srdel_in(:, :), pl_seed(:, :)
      real(rp), intent(in) :: ql_in(:, :, :), shifts(:)
      real(rp), intent(out) :: pl_out(:, :, :), rho_out(:, :, :), sumec_out(:, :), sumev_out(:, :)
      real(rp), intent(out) :: ec_out(:, :), ev_out(:, :)
      type(self) :: replay_state
      type(control), target :: replay_control
      type(symbolic_atom), allocatable, target :: replay_atoms(:)
      type(symbolic_atom) :: atom_template
      type(lattice), target :: replay_lattice
      type(bands) :: pl_bands
      integer :: nr, ncore, nval, nshifts, iv, prior_start, prior_end

      nr = size(rofi)
      ncore = size(ec_seed)
      nval = size(ev_seed)
      nshifts = size(shifts)
      if (lmax < 0 .or. lmax > 2 .or. nr < 3 .or. size(potential, 1) /= nr .or. &
          size(potential, 2) /= 2 .or. size(gravity_center, 1) /= lmax + 1 .or. &
          size(gravity_center, 2) /= 2 .or. size(c_in, 1) /= lmax + 1 .or. size(c_in, 2) /= 2 .or. &
          size(qpar_in, 1) /= lmax + 1 .or. size(qpar_in, 2) /= 2 .or. &
          size(srdel_in, 1) /= lmax + 1 .or. size(srdel_in, 2) /= 2 .or. &
          size(pl_seed, 1) /= lmax + 1 .or. size(pl_seed, 2) /= 2 .or. &
          size(ql_in, 1) /= 3 .or. size(ql_in, 2) /= lmax + 1 .or. size(ql_in, 3) /= 2 .or. &
          ncore < 0 .or. nval < 1 .or. &
          any(shape(pl_out) /= [lmax + 1, 2, nshifts]) .or. &
          any(shape(rho_out) /= [nr, 2, nshifts]) .or. &
          any(shape(sumec_out) /= [2, nshifts]) .or. any(shape(sumev_out) /= [2, nshifts]) .or. &
          any(shape(ec_out) /= [ncore, nshifts]) .or. any(shape(ev_out) /= [nval, nshifts])) then
         error stop 'legacy_accepted_pl_newrho_fixture: inconsistent captured operand dimensions'
      end if
      if (lmax_basis /= lmax) error stop 'legacy_accepted_pl_newrho_fixture: build basis does not match capture'

      call atom_template%restore_to_default()
      atom_template%element%atomic_number = z
      atom_template%element%f_core = f_core
      atom_template%potential%lmax = lmax
      if (lbound(atom_template%potential%pl, 1) /= 0 .or. ubound(atom_template%potential%pl, 1) /= lmax .or. &
          lbound(atom_template%potential%c, 1) /= 0 .or. ubound(atom_template%potential%c, 1) /= lmax .or. &
          lbound(atom_template%potential%ql, 2) /= 0 .or. ubound(atom_template%potential%ql, 2) /= lmax) then
         error stop 'legacy_accepted_pl_newrho_fixture: scratch atom potential bounds mismatch'
      end if
      atom_template%potential%gravity_center = gravity_center
      atom_template%potential%c = c_in
      atom_template%potential%qpar = qpar_in
      atom_template%potential%srdel = srdel_in
      atom_template%potential%pl = pl_seed
      atom_template%potential%ql = ql_in

      allocate(replay_atoms(1))
      replay_atoms(1) = atom_template
      replay_lattice%nbulk = 0
      replay_lattice%nrec = 1
      pl_bands%symbolic_atom => replay_atoms
      pl_bands%lattice => replay_lattice
      replay_control%hyperfine = .false.
      replay_state%control => replay_control
      allocate(replay_state%vtn(nr, 2), replay_state%vzt(nr, 2), &
               replay_state%fun2(nr, lmax + 1, 2), replay_state%phi_amp(nr, lmax + 1, 2))
      replay_state%vtn = 0.0_rp
      replay_state%vzt = 0.0_rp
      replay_state%fun2 = 0.0_rp
      replay_state%phi_amp = 0.0_rp

      prior_start = start_atom
      prior_end = end_atom
      start_atom = 1
      end_atom = 1
      do iv = 1, nshifts
         ! Intrinsic assignment deep-copies allocatable potential arrays.
         replay_atoms(1) = atom_template
         replay_state%vtn = 0.0_rp
         replay_state%vzt = 0.0_rp
         replay_state%fun2 = 0.0_rp
         replay_state%phi_amp = 0.0_rp
         replay_atoms(1)%potential%gravity_center(lmax + 1, :) = &
            gravity_center(lmax + 1, :) - shifts(iv)
         call pl_bands%calculate_pl()
         pl_out(:, :, iv) = replay_atoms(1)%potential%pl
         rho_out(:, :, iv) = 0.0_rp
         sumec_out(:, iv) = 0.0_rp
         sumev_out(:, iv) = 0.0_rp
         ec_out(:, iv) = ec_seed
         ev_out(:, iv) = ev_seed
         call replay_state%newrho(replay_atoms(1), z, lmax, a, b, nr, rofi, potential, rho_out(:, :, iv), &
                                  pl_out(:, :, iv), ql_in, sumec_out(:, iv), sumev_out(:, iv), &
                                  ec_out(:, iv), ev_out(:, iv), tol, 2, 0)
      end do
      start_atom = prior_start
      end_atom = prior_end
      nullify(pl_bands%symbolic_atom, pl_bands%lattice)
      deallocate(replay_atoms, replay_state%vtn, replay_state%vzt, replay_state%fun2, replay_state%phi_amp)
   end subroutine legacy_accepted_pl_newrho_fixture

   !> Private frozen-state M0.3 replay. Each case supplies its own per-spin
   !> center shift and complete QL moments; production calculate_pl and NEWRHO
   !> operate on a fresh deep copy for every case. The optional d-density
   !> output is the l=2 valence summand observed inside NEWRHO itself.
   subroutine legacy_accepted_m03_newrho_fixture(lmax, z, a, b, tol, f_core, rofi, potential, ec_seed, ev_seed, &
                                                 gravity_center, c_in, qpar_in, srdel_in, pl_seed, pl_live, &
                                                 ql_cases, center_shift, pl_from_calc, pl_out, rho_out, d_rho_out, &
                                                 sumec_out, sumev_out, ec_out, ev_out)
      integer, intent(in) :: lmax, f_core
      real(rp), intent(in) :: z, a, b, tol
      real(rp), intent(in) :: rofi(:), potential(:, :), ec_seed(:), ev_seed(:)
      real(rp), intent(in) :: gravity_center(:, :), c_in(:, :), qpar_in(:, :), srdel_in(:, :), pl_seed(:, :)
      real(rp), intent(in) :: pl_live(:, :)
      real(rp), intent(in) :: ql_cases(:, :, :, :), center_shift(:, :)
      logical, intent(in) :: pl_from_calc(:)
      real(rp), intent(out) :: pl_out(:, :, :), rho_out(:, :, :), d_rho_out(:, :, :)
      real(rp), intent(out) :: sumec_out(:, :), sumev_out(:, :), ec_out(:, :), ev_out(:, :)
      type(self) :: replay_state
      type(control), target :: replay_control
      type(symbolic_atom), allocatable, target :: replay_atoms(:)
      type(symbolic_atom) :: atom_template
      type(lattice), target :: replay_lattice
      type(bands) :: pl_bands
      integer :: nr, ncore, nval, ncase, iv, prior_start, prior_end, spin

      nr = size(rofi)
      ncore = size(ec_seed)
      nval = size(ev_seed)
      ncase = size(center_shift, 2)
      if (lmax /= 2 .or. nr < 3 .or. size(potential, 1) /= nr .or. size(potential, 2) /= 2 .or. &
          any(shape(gravity_center) /= [lmax + 1, 2]) .or. any(shape(c_in) /= [lmax + 1, 2]) .or. &
          any(shape(qpar_in) /= [lmax + 1, 2]) .or. any(shape(srdel_in) /= [lmax + 1, 2]) .or. &
          any(shape(pl_seed) /= [lmax + 1, 2]) .or. any(shape(pl_live) /= [lmax + 1, 2]) .or. &
          any(shape(ql_cases) /= [3, lmax + 1, 2, ncase]) .or. any(shape(center_shift) /= [2, ncase]) .or. &
          size(pl_from_calc) /= ncase .or. &
          ncore < 0 .or. nval < 1 .or. any(shape(pl_out) /= [lmax + 1, 2, ncase]) .or. &
          any(shape(rho_out) /= [nr, 2, ncase]) .or. any(shape(d_rho_out) /= [nr, 2, ncase]) .or. &
          any(shape(sumec_out) /= [2, ncase]) .or. any(shape(sumev_out) /= [2, ncase]) .or. &
          any(shape(ec_out) /= [ncore, ncase]) .or. any(shape(ev_out) /= [nval, ncase])) then
         error stop 'legacy_accepted_m03_newrho_fixture: inconsistent frozen-state dimensions'
      end if
      if (lmax_basis /= lmax) error stop 'legacy_accepted_m03_newrho_fixture: build basis differs from capture'

      call atom_template%restore_to_default()
      atom_template%element%atomic_number = z
      atom_template%element%f_core = f_core
      atom_template%potential%lmax = lmax
      atom_template%potential%gravity_center = gravity_center
      atom_template%potential%c = c_in
      atom_template%potential%qpar = qpar_in
      atom_template%potential%srdel = srdel_in
      atom_template%potential%pl = pl_seed
      atom_template%potential%ql = ql_cases(:, :, :, 1)
      allocate(replay_atoms(1))
      replay_atoms(1) = atom_template
      replay_lattice%nbulk = 0
      replay_lattice%nrec = 1
      pl_bands%symbolic_atom => replay_atoms
      pl_bands%lattice => replay_lattice
      replay_control%hyperfine = .false.
      replay_state%control => replay_control
      allocate(replay_state%vtn(nr, 2), replay_state%vzt(nr, 2), &
               replay_state%fun2(nr, lmax + 1, 2), replay_state%phi_amp(nr, lmax + 1, 2))
      replay_state%vtn = 0.0_rp
      replay_state%vzt = 0.0_rp
      replay_state%fun2 = 0.0_rp
      replay_state%phi_amp = 0.0_rp

      prior_start = start_atom
      prior_end = end_atom
      start_atom = 1
      end_atom = 1
      do iv = 1, ncase
         replay_atoms(1) = atom_template
         replay_atoms(1)%potential%ql = ql_cases(:, :, :, iv)
         do spin = 1, 2
            replay_atoms(1)%potential%gravity_center(lmax + 1, spin) = &
               gravity_center(lmax + 1, spin) - center_shift(spin, iv)
         end do
         replay_state%vtn = 0.0_rp
         replay_state%vzt = 0.0_rp
         replay_state%fun2 = 0.0_rp
         replay_state%phi_amp = 0.0_rp
         if (pl_from_calc(iv)) then
            call pl_bands%calculate_pl()
            pl_out(:, :, iv) = replay_atoms(1)%potential%pl
         else
            pl_out(:, :, iv) = pl_live
         end if
         rho_out(:, :, iv) = 0.0_rp
         d_rho_out(:, :, iv) = 0.0_rp
         sumec_out(:, iv) = 0.0_rp
         sumev_out(:, iv) = 0.0_rp
         ec_out(:, iv) = ec_seed
         ev_out(:, iv) = ev_seed
         call replay_state%newrho(replay_atoms(1), z, lmax, a, b, nr, rofi, potential, rho_out(:, :, iv), &
                                  pl_out(:, :, iv), ql_cases(:, :, :, iv), sumec_out(:, iv), sumev_out(:, iv), &
                                  ec_out(:, iv), ev_out(:, iv), tol, 2, 0, d_rho_out(:, :, iv))
      end do
      start_atom = prior_start
      end_atom = prior_end
      nullify(pl_bands%symbolic_atom, pl_bands%lattice)
      deallocate(replay_atoms, replay_state%vtn, replay_state%vzt, replay_state%fun2, replay_state%phi_amp)
   end subroutine legacy_accepted_m03_newrho_fixture

   !> Re-evaluate the real POISS0/VXC0SP component set on frozen accepted
   !> input and replay output densities. The two density stages stay separate.
   subroutine legacy_accepted_m03_components_fixture(z, a, b, nr, rofi_base, v_base, vhrmax, b_fsm, txc, &
                                                      rho_in, rho_out, sumec, sumev, components_in, &
                                                      components_out, rvh_in, rvh_out)
      real(rp), intent(in) :: z, a, b, vhrmax, b_fsm
      integer, intent(in) :: nr, txc
      real(rp), intent(in) :: rofi_base(nr), v_base(nr, 2), rho_in(nr, 2)
      real(rp), intent(in) :: rho_out(:, :, :), sumec(:, :), sumev(:, :)
      real(rp), intent(out) :: components_in(:, :), components_out(:, :), rvh_in(:, :), rvh_out(:, :)
      type(self) :: replay_state
      type(control), target :: replay_control
      type(xc) :: replay_xc
      integer :: ncase, icase

      ncase = size(rho_out, 3)
      if (size(rho_out, 1) /= nr .or. size(rho_out, 2) /= 2 .or. &
          any(shape(sumec) /= [2, ncase]) .or. any(shape(sumev) /= [2, ncase]) .or. &
          any(shape(components_in) /= [10, ncase]) .or. any(shape(components_out) /= [10, ncase]) .or. &
          any(shape(rvh_in) /= [2, ncase]) .or. any(shape(rvh_out) /= [2, ncase])) then
         error stop 'legacy_accepted_m03_components_fixture: inconsistent frozen component dimensions'
      end if
      replay_control%txc = txc
      replay_control%hyperfine = .false.
      replay_xc = xc(replay_control)
      do icase = 1, ncase
         call component_set(rho_in, sumec(:, icase), sumev(:, icase), &
                            components_in(:, icase), rvh_in(:, icase))
         call component_set(rho_out(:, :, icase), sumec(:, icase), sumev(:, icase), &
                            components_out(:, icase), rvh_out(:, icase))
      end do
   contains
      subroutine component_set(rho_sample, sumec_sample, sumev_sample, components, rvh_spin)
         real(rp), intent(in) :: rho_sample(nr, 2), sumec_sample(2), sumev_sample(2)
         real(rp), intent(out) :: components(10), rvh_spin(2)
         real(rp) :: rho_work(nr, 2), v_work(nr, 2), rofi_work(nr), rvh_local(2), vsum_local
         real(rp) :: rho0_local(2), reps_local(2), rmu_local(2), vnucl_local, rhovh_local
         real(rp) :: zvnucl_local, rhoeps_local, rhomu_local, ekin_local, utot_local, etot_local

         rho_work = rho_sample
         v_work = v_base
         rofi_work = rofi_base
         rvh_local = 0.0_rp
         vsum_local = 0.0_rp
         rho0_local = 0.0_rp
         reps_local = 0.0_rp
         rmu_local = 0.0_rp
         call POISS0(z, a, b, rofi_work, rho_work, nr, vhrmax, v_work, rvh_local, vsum_local, 2)
         vnucl_local = v_work(1, 1)
         call replay_state%VXC0SP(replay_xc, z, a, b, rofi_work, rho_sample, nr, v_work, &
                                  rho0_local, reps_local, rmu_local, 2, b_fsm)
         rhovh_local = sum(rvh_local)
         rhoeps_local = sum(reps_local)
         rhomu_local = sum(rmu_local)
         zvnucl_local = -z*vnucl_local
         utot_local = 0.5_rp*(rhovh_local + zvnucl_local)
         ekin_local = sum(sumev_sample) + sum(sumec_sample) - rhovh_local - rhomu_local
         etot_local = ekin_local + utot_local + rhoeps_local
         components = [sum(sumev_sample), sum(sumec_sample), rhovh_local, rhomu_local, vnucl_local, &
                       zvnucl_local, ekin_local, utot_local, rhoeps_local, etot_local]
         rvh_spin = rvh_local
      end subroutine component_set
   end subroutine legacy_accepted_m03_components_fixture
   pure subroutine fixed_u_r05_targets(use_kspace, p_override, v_override, test_tighten, p_target, v_target)
      logical, intent(in) :: use_kspace, test_tighten
      real(rp), intent(in) :: p_override, v_override
      real(rp), intent(out) :: p_target, v_target

      ! KS keeps its frozen 1e-7 defaults; RS uses the approved density
      ! tolerance while retaining the same potential tolerance.
      p_target = fixed_u_density_tolerance
      v_target = fixed_u_potential_tolerance
      if (.not. use_kspace) p_target = 5.0e-7_rp
      if (p_override > 0.0_rp) p_target = p_override
      if (v_override > 0.0_rp) v_target = v_override
      if (test_tighten) then
         p_target = fixed_u_test_tight_density_tolerance
         v_target = fixed_u_test_tight_potential_tolerance
      end if
   end subroutine fixed_u_r05_targets

   pure subroutine atomsc_inner_controls(test_tol, test_tolrsq, test_cap, tol, tolrsq, niter)
      real(rp), intent(in) :: test_tol, test_tolrsq
      integer, intent(in) :: test_cap
      real(rp), intent(out) :: tol, tolrsq
      integer, intent(out) :: niter

      tol = 1.0e-6_rp
      tolrsq = 1.0e-8_rp
      niter = 80
      if (test_cap > 0) then
         tol = test_tol
         tolrsq = test_tolrsq
         niter = test_cap
      end if
   end subroutine atomsc_inner_controls

   pure function atomsc_solver_tolerance(iteration, drho_control, tolrsq) result(tolerance)
      integer, intent(in) :: iteration
      real(rp), intent(in) :: drho_control, tolrsq
      real(rp) :: tolerance

      tolerance = tolrsq
      if (iteration >= 2 .and. drho_control > 2.0_rp) tolerance = 1.0e-3_rp
   end function atomsc_solver_tolerance

   pure function atomsc_mixing_beta(iteration, drho_control, beta) result(beta1)
      integer, intent(in) :: iteration
      real(rp), intent(in) :: drho_control, beta
      real(rp) :: beta1

      beta1 = beta
      if (mod(iteration, 3) == 2 .and. drho_control < 1.0_rp) beta1 = 0.5_rp
   end function atomsc_mixing_beta

   pure function atomsc_should_stop(iteration, niter, drho_control, tolerance) result(should_stop)
      integer, intent(in) :: iteration, niter
      real(rp), intent(in) :: drho_control, tolerance
      logical :: should_stop

      should_stop = drho_control < tolerance .or. iteration == niter - 1
   end function atomsc_should_stop

   pure subroutine vxc0sp_pbe_origin_derivatives(n_radial_spin_channels, rhopp, rhod, rhodd)
      integer, intent(in) :: n_radial_spin_channels
      real(rp), dimension(2), intent(in) :: rhopp
      real(rp), dimension(2), intent(out) :: rhod, rhodd

      ! RHOPP contains raw physical second radial derivatives.  PBEGGA owns
      ! the regular spherical factor of three; only one local radial channel’s
      ! total-to-per-spin
      ! conversion belongs at this VXC0SP boundary.
      rhod = 0.0_rp
      if (n_radial_spin_channels == 1) then
         rhodd = 0.5_rp*rhopp(1)
      else
         ! XCPOT/PBEGGA receives the historical down/up ordering.
         rhodd(1) = rhopp(2)
         rhodd(2) = rhopp(1)
      end if
   end subroutine vxc0sp_pbe_origin_derivatives

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
      if (allocated(this%phi_amp)) call g_safe_alloc%deallocate('self.phi_amp', this%phi_amp)
      if (allocated(this%core_rho_weighted)) call g_safe_alloc%deallocate('self.core_rho_weighted', this%core_rho_weighted)
      if (allocated(this%core_pauli_weighted)) call g_safe_alloc%deallocate('self.core_pauli_weighted', this%core_pauli_weighted)
      if (allocated(this%core_rho_weighted_l)) call g_safe_alloc%deallocate('self.core_rho_weighted_l', this%core_rho_weighted_l)
      if (allocated(this%core_pauli_weighted_l)) call g_safe_alloc%deallocate('self.core_pauli_weighted_l', this%core_pauli_weighted_l)
#else
      if (allocated(this%ws)) deallocate (this%ws)
      if (allocated(this%mixmag)) deallocate (this%mixmag)
      if (allocated(this%rb)) deallocate (this%rb)
      if (allocated(this%core_rho_weighted)) deallocate (this%core_rho_weighted)
      if (allocated(this%core_pauli_weighted)) deallocate (this%core_pauli_weighted)
      if (allocated(this%core_rho_weighted_l)) deallocate (this%core_rho_weighted_l)
      if (allocated(this%core_pauli_weighted_l)) deallocate (this%core_pauli_weighted_l)
#endif
      if (allocated(this%constraint_reference)) deallocate(this%constraint_reference)
      if (allocated(this%gbt_scf_in_moment)) deallocate(this%gbt_scf_in_moment)
      if (allocated(this%reciprocal_scf_cache)) deallocate(this%reciprocal_scf_cache)
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
      integer :: iostatus, funit, i, k
      logical :: radial_mode_defaulted, atomsc_control_active
      character(len=128) :: radial_mode_reason
      character(len=8) :: radial_mode_source

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
      r05_p_tolerance = this%r05_p_tolerance
      r05_v_tolerance = this%r05_v_tolerance
      atomsc_test_tol = this%atomsc_test_tol
      atomsc_test_tolrsq = this%atomsc_test_tolrsq
      atomsc_test_inner_cap = this%atomsc_test_inner_cap
      fix_soc = this%fix_soc
      soc_scale = this%soc_scale
      cold = this%cold
      use_kspace = this%use_kspace
      magnetic_seed_enable = this%magnetic_seed_enable
      magnetic_seed_steps = this%magnetic_seed_steps
      magnetic_seed_field = this%magnetic_seed_field

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
      this%r05_p_tolerance = r05_p_tolerance
      this%r05_v_tolerance = r05_v_tolerance
      this%atomsc_test_tol = atomsc_test_tol
      this%atomsc_test_tolrsq = atomsc_test_tolrsq
      this%atomsc_test_inner_cap = atomsc_test_inner_cap
      if (this%r05_p_tolerance > 0.0_rp .and. this%r05_p_tolerance > 1.0_rp) &
         call g_logger%fatal('self%r05_p_tolerance must be in (0,1].', __FILE__, __LINE__)
      if (this%r05_v_tolerance > 0.0_rp .and. this%r05_v_tolerance > 1.0_rp) &
         call g_logger%fatal('self%r05_v_tolerance must be in (0,1] Ry.', __FILE__, __LINE__)
      if ((this%atomsc_test_tol > 0.0_rp .or. this%atomsc_test_tolrsq > 0.0_rp .or. &
           this%atomsc_test_inner_cap > 0) .and. &
          (this%atomsc_test_tol <= 0.0_rp .or. this%atomsc_test_tolrsq <= 0.0_rp .or. &
           this%atomsc_test_inner_cap < 2)) then
         call g_logger%fatal('atomsc_test_tol, atomsc_test_tolrsq, and atomsc_test_inner_cap must be set together.', &
                             __FILE__, __LINE__)
      end if
      if (this%atomsc_test_tol > 0.0_rp) then
         if (use_kspace) call g_logger%fatal('atomsc test overrides are restricted to RS calculations.', __FILE__, __LINE__)
         atomsc_control_active = .false.
         do k = 1, size(this%symbolic_atom)
            if (allocated(this%symbolic_atom(k)%potential%hubbard_u)) &
               atomsc_control_active = atomsc_control_active .or. &
                  any(abs(this%symbolic_atom(k)%potential%hubbard_u) > 1.0e-10_rp)
            if (allocated(this%symbolic_atom(k)%potential%hubbard_j)) &
               atomsc_control_active = atomsc_control_active .or. &
                  any(abs(this%symbolic_atom(k)%potential%hubbard_j) > 1.0e-10_rp)
         end do
         if (atomsc_control_active) &
            call g_logger%fatal('atomsc test overrides are restricted to U=J=0 controls.', __FILE__, __LINE__)
      end if

      ! Constrains variables
      ! static magnetic momentum (non-linar calculation)
      this%freeze = freeze
      this%rigid_band = rigid_band
      this%fix_soc = fix_soc
      this%soc_scale = soc_scale

      ! Scaling SOC if soc_scale is defined. SOC scale default = 1.0d0
      
      do i = 1, this%lattice%ntype
         this%symbolic_atom(i)%potential%xi_p(:) = this%soc_scale*this%symbolic_atom(i)%potential%xi_p(:) 
         this%symbolic_atom(i)%potential%xi_d(:) = this%soc_scale*this%symbolic_atom(i)%potential%xi_d(:)
      end do
      ! number of loops to use rigid bands per atom
      call move_alloc(rb, this%rb)

      ! other variables
      this%orbital_polarization = orbital_polarization
      this%init = init
      this%cold = cold
      this%use_kspace = use_kspace
      radial_mode_defaulted = trim(this%hamiltonian%hubbard_radial_moments) == 'auto'
      if (trim(this%hamiltonian%hubbard_radial_moments) == 'auto') then
         if (this%is_eligible_local_moments(radial_mode_reason)) then
            this%hamiltonian%hubbard_radial_moments = 'local'
         else
            this%hamiltonian%hubbard_radial_moments = 'spectral'
            if (rank == 0 .and. (this%control%calctype == 'I' .or. this%control%calctype == 'S')) &
               call g_logger%info('HUBBARD_RADIAL_MODE auto retained spectral: '//trim(radial_mode_reason), __FILE__, __LINE__)
         end if
      end if
      radial_mode_source = 'explicit'
      if (radial_mode_defaulted) radial_mode_source = 'default'
      if (rank == 0) call g_logger%info('HUBBARD_RADIAL_MODE mode='//trim(this%hamiltonian%hubbard_radial_moments)// &
         ' source='//trim(radial_mode_source)//' use_kspace='//trim(log2str(this%use_kspace))// &
         ' local_energy_adopted='//trim(log2str(this%is_eligible_local_hubbard() .and. &
            trim(this%hamiltonian%hubbard_radial_moments) == 'local')), __FILE__, __LINE__)
      call this%validate_hubbard_mode_provenance()
      this%magnetic_seed_enable = magnetic_seed_enable
      this%magnetic_seed_steps = magnetic_seed_steps
      this%magnetic_seed_field = magnetic_seed_field
      this%magnetic_seed_active = .false.
      if (this%magnetic_seed_enable) then
         if (associated(this%control) .and. this%control%constraints_enable) then
            call g_logger%fatal('magnetic_seed_enable cannot be combined with persistent magnetic constraints.', __FILE__, __LINE__)
         end if
         if (this%magnetic_seed_steps < 1) then
            call g_logger%fatal('magnetic_seed_steps must be positive when magnetic_seed_enable is true.', __FILE__, __LINE__)
         end if
         if (abs(this%magnetic_seed_field) <= tiny(1.0_rp)) then
            call g_logger%fatal('magnetic_seed_field must be nonzero when magnetic_seed_enable is true.', __FILE__, __LINE__)
         end if
         if (this%nstep > 0 .and. this%magnetic_seed_steps >= this%nstep) then
            call g_logger%fatal('magnetic_seed_steps must leave at least one fully unconstrained SCF iteration.', __FILE__, __LINE__)
         end if
      end if
      ! Initialize constraining facility if requested in control namelist
      if (associated(this%control)) then
         if (this%control%constraints_enable) then
            call initialize_cfd(this%lattice%nrec, 1, this%control%constraints_i_cons, this%control%constraints_code_prefac, &
                                this%control%constraints_diagnostics)
            call this%initialize_constraint_state()
         end if
         if (this%control%gbt_scf_diagnostics) call this%initialize_gbt_scf_diagnostics()
      end if
      if (associated(this%hamiltonian) .and. associated(this%control)) then
         if (trim(this%hamiltonian%magnetic_representation) == gbt_single_q .and. &
             trim(this%control%density_policy) == sd_relaxed_reference) then
            call g_logger%fatal('gbt_single_q with density_policy="relaxed_reference" is unsupported: '// &
                                'the GBT Hamiltonian requires a fixed rotating-frame reference axis.', &
                                __FILE__, __LINE__)
         end if
      end if
   end subroutine build_from_file

   !> Initialize the opt-in WP12 end-to-end SCF trace.
   !> The trace is one row per site and iteration, assembled after the radial
   !> update so old density, newly accumulated density, mixer state, next
   !> radial channels, Hamiltonian exchange, and separated energies are visible
   !> in a single record.
   subroutine initialize_gbt_scf_diagnostics(this)
      class(self), intent(inout) :: this
      integer :: unit, io_status

      this%gbt_scf_diagnostics_header_written = .false.
      if (rank /= 0) return
      open(newunit=unit, file='gbt_scf_diagnostics.dat', status='replace', action='write', iostat=io_status)
      if (io_status /= 0) then
         call g_logger%warning('Could not open gbt_scf_diagnostics.dat; WP12 trace disabled for this run.', &
                               __FILE__, __LINE__)
         return
      end if
      write(unit, '(a)') '# frame=GBT rotating primitive frame for gbt_single_q; ordinary active frame otherwise; energies=Ry'
      write(unit, '(a)') '# one row per site and SCF iteration; radial columns are summed l-channel ql(1,:,spin) values'
      write(unit, '(a)') '# iteration site policy representation in_frame_x in_frame_y in_frame_z in_moment_x in_moment_y in_moment_z incoming_magnitude in_ql_up in_ql_dn output_moment_x output_moment_y output_moment_z output_magnitude output_charge rho_uu rho_dd rho_ud_re rho_ud_im output_ql_up output_ql_dn projection_axis_x projection_axis_y projection_axis_z longitudinal_moment transverse_magnitude mixer_old_x mixer_old_y mixer_old_z mixer_new_x mixer_new_y mixer_new_z mixer_mix_x mixer_mix_y mixer_mix_z next_ql_up next_ql_dn next_potential_mom_x next_potential_mom_y next_potential_mom_z exchange_c0 exchange_w0 target_x target_y target_z field_x field_y field_z physical_total_energy constraint_penalty_energy constraint_field_coupling_energy constraint_metric'
      close(unit)
      this%gbt_scf_diagnostics_header_written = .true.
   end subroutine initialize_gbt_scf_diagnostics

   !> Set or release the explicitly requested temporary B_fsm seed.
   !> No constraint target, mixer state, or magnetic moment is modified here.
   subroutine set_magnetic_seed_state(this, iteration)
      class(self), intent(inout) :: this
      integer, intent(in) :: iteration
      logical :: active_now

      active_now = this%magnetic_seed_enable .and. iteration <= this%magnetic_seed_steps
      if (active_now .and. .not. this%magnetic_seed_active) then
         if (rank == 0) call g_logger%info('MAGNETIC_SEED temporary symmetry-breaking seed active: outer iteration '// &
                                           int2str(iteration)//', B_fsm='//real2str(this%magnetic_seed_field), __FILE__, __LINE__)
      else if (.not. active_now .and. this%magnetic_seed_active) then
         if (rank == 0) call g_logger%info('MAGNETIC_SEED released; subsequent SCF iterations are fully unconstrained.', &
                                           __FILE__, __LINE__)
      end if
      this%magnetic_seed_active = active_now
   end subroutine set_magnetic_seed_state

   !> Capture the fixed reference directions and the optional initial field.
   !> The state is stored on the same symbolic atoms consumed by the RS and
   !> reciprocal Hamiltonian builders, so a field is neither lost at the call
   !> boundary nor applied twice through separate potential paths.
   subroutine initialize_constraint_state(this)
      class(self), intent(inout) :: this
      integer :: ia, plusbulk

      integer :: unit, io_status

      if (allocated(this%constraint_reference)) deallocate(this%constraint_reference)
      allocate(this%constraint_reference(3, this%lattice%nrec))
      this%constraint_reference = 0.0_rp
      if (allocated(this%gbt_scf_in_moment)) deallocate(this%gbt_scf_in_moment)
      allocate(this%gbt_scf_in_moment(3, this%lattice%nrec))
      this%gbt_scf_in_moment = 0.0_rp
      do ia = 1, this%lattice%nrec
         plusbulk = this%lattice%nbulk + ia
         if (allocated(this%control%constraints_mom_ref) .and. &
             .not. (size(this%control%constraints_mom_ref, 1) == 1 .and. &
                    size(this%control%constraints_mom_ref, 2) == 1 .and. &
                    abs(this%control%constraints_mom_ref(1, 1)) <= tiny(1.0_rp))) then
            if (size(this%control%constraints_mom_ref, 1) /= 3 .or. &
                size(this%control%constraints_mom_ref, 2) /= this%lattice%nrec) then
               call g_logger%fatal('constraints_mom_ref must have shape (3,nrec); got ('// &
                                   fmt('i0', size(this%control%constraints_mom_ref, 1))//','// &
                                   fmt('i0', size(this%control%constraints_mom_ref, 2))//').', __FILE__, __LINE__)
            end if
            this%constraint_reference(:, ia) = this%control%constraints_mom_ref(:, ia)
         else
            this%constraint_reference(:, ia) = this%symbolic_atom(plusbulk)%potential%mom0(:)
         end if
         if (norm2(this%constraint_reference(:, ia)) <= tiny(1.0_rp)) then
            call g_logger%fatal('Constraint target for site '//fmt('i0', ia)//' has zero length.', __FILE__, __LINE__)
         end if
         this%symbolic_atom(plusbulk)%constraint_target = this%constraint_reference(:, ia)

         this%symbolic_atom(plusbulk)%mag_cfield = 0.0_rp
         if (allocated(this%control%constraints_bfield) .and. &
             .not. (size(this%control%constraints_bfield, 1) == 1 .and. &
                    size(this%control%constraints_bfield, 2) == 1 .and. &
                    abs(this%control%constraints_bfield(1, 1)) <= tiny(1.0_rp))) then
            if (size(this%control%constraints_bfield, 1) /= 3 .or. &
                size(this%control%constraints_bfield, 2) /= this%lattice%nrec) then
               call g_logger%fatal('constraints_bfield must have shape (3,nrec); got ('// &
                                   fmt('i0', size(this%control%constraints_bfield, 1))//','// &
                                   fmt('i0', size(this%control%constraints_bfield, 2))//').', __FILE__, __LINE__)
            end if
            this%symbolic_atom(plusbulk)%mag_cfield = this%control%constraints_bfield(:, ia)
         end if
      end do
      this%constraint_tolerance = this%control%constraints_tolerance
      this%constraint_metric = huge(1.0_rp)
      this%constraint_converged = .false.
      this%constraint_iteration = 0
      this%constraint_diagnostics_header_written = .false.
      this%gbt_scf_diagnostics_header_written = .false.
      if (this%control%constraints_diagnostics .and. rank == 0) then
         open(newunit=unit, file='constraint_diagnostics.dat', status='replace', action='write', iostat=io_status)
         if (io_status == 0) then
            write(unit, '(a)') '# frame=active_spin_frame; field_units=Ry; angle_units=rad'
            write(unit, '(a)') '# iteration site target_x target_y target_z actual_x actual_y actual_z moment_magnitude angular_error transverse_x transverse_y transverse_z Bx By Bz B_magnitude B_dot_target torque_x torque_y torque_z rms_angle'
            close(unit)
            this%constraint_diagnostics_header_written = .true.
         end if
      end if
   end subroutine initialize_constraint_state

   !> Write one complete WP12 SCF record per self-consistent site.
   !> `mix%qia_old` and `mix%qia_new` retain the radial channels entering and
   !> leaving the density step, while the current potential contains the
   !> channels that will enter the next Hamiltonian after `run_scf`.
   subroutine write_gbt_scf_diagnostics(this, iteration)
      class(self), intent(in) :: this
      integer, intent(in) :: iteration

      integer :: ia, plusbulk, lcount, unit, io_status
      real(rp) :: in_ql_up, in_ql_dn, out_ql_up, out_ql_dn
      real(rp) :: next_ql_up, next_ql_dn, in_moment(3), out_moment(3)
      real(rp) :: in_frame(3), target(3), axis(3), field(3), next_mom(3)
      real(rp) :: longitudinal, transverse, exchange_c0, exchange_w0, axis_norm
      real(rp) :: output_magnitude, incoming_magnitude, output_charge
      complex(rp) :: rho_uu, rho_dd, rho_ud

      if (.not. this%control%gbt_scf_diagnostics .or. rank /= 0) return
      if (.not. this%gbt_scf_diagnostics_header_written) return

      open(newunit=unit, file='gbt_scf_diagnostics.dat', status='unknown', action='write', &
           position='append', iostat=io_status)
      if (io_status /= 0) return

      do ia = 1, this%lattice%nrec
         plusbulk = this%lattice%nbulk + ia
         lcount = min(this%symbolic_atom(plusbulk)%potential%lmax, lmax_basis) + 1
         in_frame(:) = this%mix%mag_old(ia, :)
         if (sqrt(sum(in_frame(:)**2)) <= tiny(1.0_rp)) in_frame(:) = [0.0_rp, 0.0_rp, 1.0_rp]
         in_moment(:) = this%gbt_scf_in_moment(:, ia)
         incoming_magnitude = sqrt(sum(in_moment(:)**2))
         out_moment(:) = this%symbolic_atom(plusbulk)%potential%mom0(:)
         output_magnitude = this%symbolic_atom(plusbulk)%potential%mtot

         target(:) = this%symbolic_atom(plusbulk)%constraint_target(:)
         if (sqrt(sum(target(:)**2)) <= tiny(1.0_rp)) target(:) = in_frame(:)
         axis_norm = sqrt(sum(target(:)**2))
         if (axis_norm <= tiny(1.0_rp)) then
            axis(:) = [0.0_rp, 0.0_rp, 1.0_rp]
         else
            axis(:) = target(:)/axis_norm
         end if
         longitudinal = sum(out_moment(:)*axis(:))
         transverse = sqrt(max(sum(out_moment(:)**2) - longitudinal*longitudinal, 0.0_rp))
         field(:) = this%symbolic_atom(plusbulk)%mag_cfield(:)
         next_mom(:) = this%symbolic_atom(plusbulk)%potential%mom(:)

         in_ql_up = sum(this%mix%qia_old(ia, 1:lcount))
         in_ql_dn = sum(this%mix%qia_old(ia, lcount + 1:2*lcount))
         out_ql_up = sum(this%mix%qia_new(ia, 1:lcount))
         out_ql_dn = sum(this%mix%qia_new(ia, lcount + 1:2*lcount))
         output_charge = out_ql_up + out_ql_dn
         next_ql_up = sum(this%symbolic_atom(plusbulk)%potential%ql(1, 0:this%symbolic_atom(plusbulk)%potential%lmax, 1))
         next_ql_dn = sum(this%symbolic_atom(plusbulk)%potential%ql(1, 0:this%symbolic_atom(plusbulk)%potential%lmax, 2))

         exchange_c0 = 0.0_rp
         exchange_w0 = 0.0_rp
         if (allocated(this%symbolic_atom(plusbulk)%potential%cx1)) exchange_c0 = real(this%symbolic_atom(plusbulk)%potential%cx1(1), rp)
         if (allocated(this%symbolic_atom(plusbulk)%potential%wx1)) exchange_w0 = real(this%symbolic_atom(plusbulk)%potential%wx1(1), rp)

         rho_uu = (0.0_rp, 0.0_rp)
         rho_dd = (0.0_rp, 0.0_rp)
         rho_ud = (0.0_rp, 0.0_rp)
         if (this%use_kspace) then
            if (allocated(this%reciprocal_scf_cache)) then
               if (allocated(this%reciprocal_scf_cache%rf_density%rho)) then
                  rho_uu = sum(this%reciprocal_scf_cache%rf_density%rho(1, 1, :, ia, 1))
                  rho_dd = sum(this%reciprocal_scf_cache%rf_density%rho(2, 2, :, ia, 1))
                  rho_ud = sum(this%reciprocal_scf_cache%rf_density%rho(1, 2, :, ia, 1))
               end if
            end if
         else if (allocated(this%bands%rf_density%rho)) then
            rho_uu = sum(this%bands%rf_density%rho(1, 1, :, ia, 1))
            rho_dd = sum(this%bands%rf_density%rho(2, 2, :, ia, 1))
            rho_ud = sum(this%bands%rf_density%rho(1, 2, :, ia, 1))
         end if

         write(unit, '(2i8,1x,a,1x,a,51(1x,es20.12))') iteration, ia, &
            trim(this%control%density_policy), trim(this%hamiltonian%magnetic_representation), &
            in_frame(:), in_moment(:), incoming_magnitude, in_ql_up, in_ql_dn, out_moment(:), output_magnitude, &
            output_charge, real(rho_uu, rp), real(rho_dd, rp), real(rho_ud, rp), aimag(rho_ud), &
            out_ql_up, out_ql_dn, axis(:), longitudinal, transverse, this%mix%mag_old(ia, :), &
            this%mix%mag_new(ia, :), this%mix%mag_mix(ia, :), next_ql_up, next_ql_dn, next_mom(:), &
            exchange_c0, exchange_w0, target(:), field(:), this%physical_total_energy, &
            this%constraint_penalty_energy, this%constraint_field_coupling_energy, this%constraint_metric
      end do
      close(unit)
   end subroutine write_gbt_scf_diagnostics

   !> Apply the SCF magnetic-state contract at the mixer boundary.
   !> Under `constrained_spiral`, `mag_mix` is only a frame marker.  The
   !> longitudinal magnitude is carried by the mixed radial ql channels; the
   !> transverse output remains the measured residual and never rotates the
   !> reference used by the next Hamiltonian/density projection.  The legacy
   !> Cartesian mixer remains the explicit `relaxed_reference` behavior.
   subroutine mix_magnetic_state(this)
      class(self), intent(inout) :: this
      integer :: ia, plusbulk
      real(rp) :: axis(3), axis_norm

      if (trim(this%control%density_policy) == sd_constrained_spiral) then
         do ia = 1, this%lattice%nrec
            plusbulk = this%lattice%nbulk + ia
            axis(:) = this%mix%mag_old(ia, :)
            if (allocated(this%constraint_reference)) then
               if (sqrt(sum(this%constraint_reference(:, ia)**2)) > tiny(1.0_rp)) then
                  axis(:) = this%constraint_reference(:, ia)
               end if
            end if
            axis_norm = sqrt(sum(axis(:)**2))
            if (axis_norm <= tiny(1.0_rp)) axis(:) = [0.0_rp, 0.0_rp, 1.0_rp]
            if (axis_norm > tiny(1.0_rp)) axis(:) = axis(:)/axis_norm
            this%mix%mag_mix(ia, :) = axis(:)
            ! Keep this assignment explicit: `potential%mom` is the marker
            ! consumed by the Hamiltonian, not the transverse measured moment.
            this%symbolic_atom(plusbulk)%potential%mom(:) = axis(:)
         end do
      else if (trim(this%control%density_policy) == sd_relaxed_reference) then
         call this%mix%mix_magnetic_moments(this%mix%mag_old, this%mix%mag_new, &
                                            this%mix%mag_mix, this%symbolic_atom(:)%potential%mtot)
         do ia = 1, this%lattice%nrec
            plusbulk = this%lattice%nbulk + ia
            this%symbolic_atom(plusbulk)%potential%mom(:) = this%mix%mag_mix(ia, :)
         end do
      else
         call g_logger%fatal('Unknown density policy at SCF magnetic mixing boundary: '// &
                             trim(this%control%density_policy), __FILE__, __LINE__)
      end if
   end subroutine mix_magnetic_state

   !> Compute the next constraining field from the mixed moments and retain it
   !> on the symbolic atoms for the next Hamiltonian rebuild.
   subroutine apply_constraints(this)
      class(self), intent(inout) :: this
      real(rp), allocatable :: mom_in(:, :), mom_ref(:, :), bfield(:, :)
      real(rp) :: etcon
      type(constraint_diagnostics) :: diagnostics
      integer :: ia, plusbulk, unit, io_status
      real(rp) :: marker_sign

      if (.not. this%control%constraints_enable) then
         this%constraint_penalty_energy = 0.0_rp
         this%constraint_field_coupling_energy = 0.0_rp
         this%constraint_energy = 0.0_rp
         return
      end if

      allocate(mom_in(3, this%lattice%nrec), mom_ref(3, this%lattice%nrec), bfield(3, this%lattice%nrec))
      do ia = 1, this%lattice%nrec
         plusbulk = this%lattice%nbulk + ia
         ! The Hamiltonian keeps a collinear reference marker in GBT so the
         ! operator boundary cannot mistake a lab-frame transverse moment for
         ! the rotating-frame reference.  The actual measured vector remains
         ! the controller input through this separate state slot.
         if (trim(this%hamiltonian%magnetic_representation) == gbt_single_q .and. &
             norm2(this%symbolic_atom(plusbulk)%constraint_actual) > tiny(1.0_rp)) then
            mom_in(:, ia) = this%symbolic_atom(plusbulk)%constraint_actual(:)
         else
            ! Pass the measured moment, not its normalized direction.  The
            ! controller normalizes where required; the un-normalized vector
            ! is needed for the field-coupling diagnostic.
            mom_in(:, ia) = [this%symbolic_atom(plusbulk)%potential%mx, &
                             this%symbolic_atom(plusbulk)%potential%my, &
                             this%symbolic_atom(plusbulk)%potential%mz]
            if (norm2(mom_in(:, ia)) <= tiny(1.0_rp)) then
               mom_in(:, ia) = this%symbolic_atom(plusbulk)%potential%mom(:)
            end if
         end if
         mom_ref(:, ia) = this%constraint_reference(:, ia)
         bfield(:, ia) = this%symbolic_atom(plusbulk)%mag_cfield
      end do

      call constrain(mom_in, mom_ref, bfield, this%lattice%nrec, etcon, diagnostics_out=diagnostics)
      do ia = 1, this%lattice%nrec
         plusbulk = this%lattice%nbulk + ia
         this%symbolic_atom(plusbulk)%mag_cfield = bfield(:, ia)
         this%symbolic_atom(plusbulk)%constraint_actual = mom_in(:, ia)
         if (trim(this%hamiltonian%magnetic_representation) == gbt_single_q) then
            marker_sign = 1.0_rp
            ! In the fixed-axis policy the GBT marker follows the target
            ! reference, never the measured transverse/longitudinal output.
            if (trim(this%control%density_policy) == sd_constrained_spiral) then
               if (this%constraint_reference(3, ia) < 0.0_rp) marker_sign = -1.0_rp
            else if (mom_in(3, ia) < 0.0_rp) then
               marker_sign = -1.0_rp
            end if
            this%symbolic_atom(plusbulk)%potential%mom(:) = [0.0_rp, 0.0_rp, marker_sign]
         end if
      end do
      this%constraint_penalty_energy = diagnostics%controller_penalty_energy
      this%constraint_field_coupling_energy = diagnostics%field_coupling_energy
      ! Keep the old member as a source-compatible alias.  It is explicitly
      ! the controller penalty, never a contribution to physical_total_energy.
      this%constraint_energy = etcon
      this%constraint_metric = diagnostics%convergence_metric
      this%constraint_converged = this%constraint_metric <= this%constraint_tolerance
      this%constraint_iteration = diagnostics%iteration
      if (rank == 0) then
         call g_logger%info('CONSTRAINT penalty/controller energy (Ry) = '//fmt('f16.10', this%constraint_penalty_energy)// &
                            ' field coupling <Vcon> (Ry) = '//fmt('f16.10', this%constraint_field_coupling_energy)// &
                            ' RMS-angle = '//fmt('es12.4', this%constraint_metric), __FILE__, __LINE__)
      end if

      if (this%control%constraints_diagnostics .and. rank == 0) then
         open(newunit=unit, file='constraint_diagnostics.dat', status='unknown', action='write', &
              position='append', iostat=io_status)
         if (io_status == 0) then
            if (.not. this%constraint_diagnostics_header_written) then
               write(unit, '(a)') '# frame=active_spin_frame; field_units=Ry; angle_units=rad'
               write(unit, '(a)') '# iteration site target_x target_y target_z actual_x actual_y actual_z moment_magnitude angular_error transverse_x transverse_y transverse_z Bx By Bz B_magnitude B_dot_target torque_x torque_y torque_z rms_angle'
               this%constraint_diagnostics_header_written = .true.
            end if
            do ia = 1, this%lattice%nrec
               write(unit, '(2i8,20(1x,es20.12))') diagnostics%iteration, ia, &
                  diagnostics%target(:, ia), diagnostics%actual(:, ia), diagnostics%moment_magnitude(ia), &
                  diagnostics%angular_error(ia), diagnostics%transverse_residual(:, ia), diagnostics%bfield(:, ia), &
                  diagnostics%bfield_magnitude(ia), diagnostics%bfield_longitudinal(ia), diagnostics%torque(:, ia), &
                  diagnostics%convergence_metric
            end do
            close(unit)
         end if
      end if

      deallocate(mom_in, mom_ref, bfield)
   end subroutine apply_constraints

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this, full)
      class(self), intent(inout):: this
      logical, intent(in), optional :: full
      integer :: lmax, nfun_l

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
      this%r05_p_tolerance = -1.0_rp
      this%r05_v_tolerance = -1.0_rp
      this%atomsc_test_tol = 0.0_rp
      this%atomsc_test_tolrsq = 0.0_rp
      this%atomsc_test_inner_cap = 0

      ! Constrains variables
      ! static magnetic momentum (non-linar calculation)
      this%freeze = .false.
      this%rigid_band = .false.
      this%fix_soc = .false. 
      this%soc_scale = 1.0d0
      ! number of loops to use rigid bands per atom
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('self.rb', this%rb, this%lattice%nrec)
#else
      allocate (this%rb(this%lattice%nrec))
#endif
      this%rb(:) = 2

      this%orbital_polarization = .false.

      this%ws_max = 9.99d0

      this%cold = .false.
      this%converged = .false.
      this%converged_iteration = 0
      this%use_kspace = .false.
      this%magnetic_seed_enable = .false.
      this%magnetic_seed_steps = 0
      this%magnetic_seed_field = 0.0_rp
      this%magnetic_seed_active = .false.
      if (allocated(this%constraint_reference)) deallocate(this%constraint_reference)
      allocate(this%constraint_reference(3, this%lattice%nrec))
      this%constraint_reference = 0.0_rp
      if (allocated(this%gbt_scf_in_moment)) deallocate(this%gbt_scf_in_moment)
      allocate(this%gbt_scf_in_moment(3, this%lattice%nrec))
      this%gbt_scf_in_moment = 0.0_rp
      this%physical_total_energy = 0.0_rp
      this%constraint_penalty_energy = 0.0_rp
      this%constraint_field_coupling_energy = 0.0_rp
      this%constraint_energy = 0.0_rp
      this%constraint_metric = 0.0_rp
      this%constraint_tolerance = 1.0e-6_rp
      this%constraint_converged = .false.
      this%constraint_iteration = 0
      this%constraint_diagnostics_header_written = .false.
      this%gbt_scf_diagnostics_header_written = .false.

      if (associated(this%lattice)) then
         if (present(full)) then
            if (full) then
               call this%lattice%restore_to_default(full)
            end if
         end if
      end if

#ifdef USE_SAFE_ALLOC
      nfun_l = max(3, this%control%lmax + 1)
      call g_safe_alloc%allocate('self.bxc', this%bxc, this%lattice%nrec)
      call g_safe_alloc%allocate('self.vtn', this%vtn, (/8001, 2/))
      call g_safe_alloc%allocate('self.vzt', this%vzt, (/8001, 2/))
      call g_safe_alloc%allocate('self.fun2', this%fun2, (/8001, nfun_l, 2/))
      call g_safe_alloc%allocate('self.phi_amp', this%phi_amp, (/8001, nfun_l, 2/))
      call g_safe_alloc%allocate('self.core_rho_weighted', this%core_rho_weighted, (/8001, 2/))
      call g_safe_alloc%allocate('self.core_pauli_weighted', this%core_pauli_weighted, (/8001, 2/))
      call g_safe_alloc%allocate('self.core_rho_weighted_l', this%core_rho_weighted_l, (/8001, nfun_l, 2/))
      call g_safe_alloc%allocate('self.core_pauli_weighted_l', this%core_pauli_weighted_l, (/8001, nfun_l, 2/))
#else
      nfun_l = max(3, this%control%lmax + 1)
      allocate (this%bxc(this%lattice%nrec))
      allocate (this%vtn(8001, 2), this%vzt(8001, 2), this%fun2(8001, nfun_l, 2))
      allocate (this%phi_amp(8001, nfun_l, 2))
      allocate (this%core_rho_weighted(8001, 2), this%core_pauli_weighted(8001, 2))
      allocate (this%core_rho_weighted_l(8001, nfun_l, 2), this%core_pauli_weighted_l(8001, nfun_l, 2))
#endif
      this%fun2(:, :, :) = 0.0_rp
      this%phi_amp(:, :, :) = 0.0_rp
      this%core_rho_weighted(:, :) = 0.0_rp
      this%core_pauli_weighted(:, :) = 0.0_rp
      this%core_rho_weighted_l(:, :, :) = 0.0_rp
      this%core_pauli_weighted_l(:, :, :) = 0.0_rp

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
      print *, 'cold                 ', this%cold

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
      cold = this%cold
      magnetic_seed_enable = this%magnetic_seed_enable
      magnetic_seed_steps = this%magnetic_seed_steps
      magnetic_seed_field = this%magnetic_seed_field

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
      cold = this%cold
      magnetic_seed_enable = this%magnetic_seed_enable
      magnetic_seed_steps = this%magnetic_seed_steps
      magnetic_seed_field = this%magnetic_seed_field
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

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Performs the self-consistent calculation.
   !>
   !> Performs the self-consistent calculation.
   !---------------------------------------------------------------------------
   subroutine run(this)
      class(self), intent(inout) :: this
      integer :: i, ia, plusbulk, niter
      real(rp), dimension(6) :: QSL
      real(rp), dimension(:), allocatable :: pot_arr
         complex(rp), allocatable :: fixed_u_m1(:, :, :, :, :)
         real(rp), allocatable :: fixed_u_pl_diag(:, :, :, :)
      integer :: na_glob, pot_size
      real(rp), dimension(:, :), allocatable :: T_comm
      complex(rp), allocatable :: fixed_u_vin(:, :, :)
      real(rp) :: total_energy, magnetic_moment

      call this%initialize_fixed_u_trace()
   
      !===========================================================================
      !                              BEGIN SCF LOOP
      !===========================================================================
      if (.not. g_scf_benchmark_profile%enabled) then
         call g_scf_benchmark_profile%configure_from_environment()
      end if
      niter = 0
      do i = 1, this%nstep
         call g_scf_benchmark_profile%start_iteration(i)
         call this%set_magnetic_seed_state(i)
         if (this%cold .and. i==1) then
            if (this%use_kspace) then
               call g_scf_benchmark_profile%start_stage('P_hamiltonian_prepare')
               call run_scf(this)
               call g_scf_benchmark_profile%stop_stage('P_hamiltonian_prepare')
            else
               call g_scf_benchmark_profile%start_stage('P_rs_hamiltonian_prepare')
               call run_scf(this)
               call g_scf_benchmark_profile%stop_stage('P_rs_hamiltonian_prepare')
            end if
         end if
         !=========================================================================
         !                        PERFORM THE RECURSION
         !=========================================================================
         call run_recursion(this, i)
         this%fixed_u_trace_iteration = i
         this%fixed_u_m0_snapshot_valid = .false.
         if (this%fixed_u_trace_selected_iteration == 0) call this%reset_fixed_u_trace()
         call this%write_fixed_u_trace('Pin_H')
   
         !=========================================================================
         !               SAVE THE TOTAL ENERGY FROM PREVIOUS ITERATION
         !=========================================================================
         ! Constraint penalty is a controller diagnostic, not part of the
         ! ordinary DFT total-energy convergence signal.
         this%esumn = sum(this%symbolic_atom(:)%potential%etot)
         this%esum = this%esumn

         !=========================================================================
         !            SAVE THE PARAMETERS QL AND PL TO BE MIXED LATER
         !=========================================================================
         call this%mix%save_to('old') ! Save to qia_old to mix with qia_new.
         if (this%hamiltonian%hubbard_u_general_check .or. this%hamiltonian%hubbard_u_sc_check .or. this%hamiltonian%hubbard_v_check) then
            call this%mix%save_ldm_to('old') ! Save LDA+U density matrices before DOS step.
         end if

         !=========================================================================
         !                      SAVE THE MAGNETIC MOMENTS
         !=========================================================================
         do ia = 1, this%lattice%nrec
            plusbulk = this%lattice%nbulk + ia
            if (trim(this%control%density_policy) == sd_constrained_spiral) then
               if (allocated(this%constraint_reference) .and. &
                   sqrt(sum(this%constraint_reference(:, ia)**2)) > tiny(1.0_rp)) then
                  this%mix%mag_old(ia, :) = this%constraint_reference(:, ia) / &
                                            sqrt(sum(this%constraint_reference(:, ia)**2))
               else
                  this%mix%mag_old(ia, :) = this%symbolic_atom(plusbulk)%potential%mom(:)
               end if
            else
               if (trim(this%hamiltonian%magnetic_representation) == gbt_single_q .and. &
                   norm2(this%symbolic_atom(plusbulk)%constraint_actual) > tiny(1.0_rp)) then
                  this%mix%mag_old(ia, :) = this%symbolic_atom(plusbulk)%constraint_actual(:) / &
                                            norm2(this%symbolic_atom(plusbulk)%constraint_actual)
               else
                  this%mix%mag_old(ia, :) = this%symbolic_atom(plusbulk)%potential%mom(:)
               end if
            end if
            this%gbt_scf_in_moment(:, ia) = this%symbolic_atom(plusbulk)%potential%mom0(:)
         end do
   
         !=========================================================================
         !                  CALCULATE THE DENSITY OF STATES
         !=========================================================================
         if (this%fixed_u_trace_enabled .and. &
             (this%fixed_u_trace_selected_iteration == 0 .or. this%fixed_u_trace_selected_iteration == i)) then
            allocate(fixed_u_m1(this%lattice%nrec, lmax_basis + 1, 2, 2*lmax_basis + 1, 2*lmax_basis + 1))
            allocate(fixed_u_pl_diag(this%lattice%nrec, lmax_basis + 1, 2, 6))
            if (allocated(this%hamiltonian%hubbard_u_pot)) then
               allocate(fixed_u_vin, mold=this%hamiltonian%hubbard_u_pot)
               fixed_u_vin = this%hamiltonian%hubbard_u_pot
            else
               allocate(fixed_u_vin(nb, nb, this%lattice%ntype))
               fixed_u_vin = (0.0_rp, 0.0_rp)
            end if
            call run_dos(this, fixed_u_m1, fixed_u_pl_diag, i)
            call this%write_fixed_u_trace('P_OUT', fixed_u_m1, fixed_u_pl_diag)
            if (.not. this%fixed_u_trace_ql_only) &
               call this%snapshot_fixed_u_m0_moments(fixed_u_m1, fixed_u_pl_diag, fixed_u_vin)
            deallocate(fixed_u_vin)
            deallocate(fixed_u_m1, fixed_u_pl_diag)
         else
            call run_dos(this, iteration=i)
            call this%write_fixed_u_trace('P_OUT')
         end if
   
         !=========================================================================
         !                   MIX OLD AND NEW CALCULATED PL AND QL
         !=========================================================================
         if (rank == 0) call g_logger%info('Mixtype is '//trim(this%mix%mixtype), __FILE__, __LINE__)
         call g_scf_benchmark_profile%start_stage('P_mixing')
         call this%mix%mixpq(this%mix%qia_old, this%mix%qia_new, iteration=i, &
                             guard_collinear=this%control%nsp == 1) ! Mix qia_new with qia_old
         if (this%hamiltonian%hubbard_u_general_check .or. this%hamiltonian%hubbard_u_sc_check .or. this%hamiltonian%hubbard_v_check) then
            call this%mix%save_ldm_to('new') ! Save freshly computed LDA+U density matrices.
            call this%update_fixed_u_residuals()
            call this%mix%mix_ldm_linear() ! Optional linear mixing controlled by mix%ldm_beta.
         end if
         call g_scf_benchmark_profile%stop_stage('P_mixing')
   
         !=========================================================================
         !         CALCULATE THE MADELUNG POTENTIAL (BULK ONLY IMPLEMENTED)
         !=========================================================================
         call g_scf_benchmark_profile%start_stage('P_potential_update')
         select case (this%control%calctype)
         case ('B')
            call this%charge%bulkpot()
         case ('S')
            call this%charge%surfpot()
         case ('I')
            call this%charge%imppot()
         case ('L')
            call this%charge%interfacepot()
         end select
   
         !=========================================================================
         !                        SAVE MIXED PARAMETERS
         !=========================================================================
         call this%mix%save_to('current') ! Save mixed parameters into potential%pl and potential%ql
         if (this%hamiltonian%hubbard_u_general_check .or. this%hamiltonian%hubbard_u_sc_check .or. this%hamiltonian%hubbard_v_check) then
            call this%mix%save_ldm_to('current') ! Save mixed density matrices into potential%ldm.
         end if
         call this%write_fixed_u_trace('P_MIX')
   
         !=========================================================================
         !                       MAKE SFC ATOMIC SPHERE
         !=========================================================================
         call run_scf(this)
         call g_scf_benchmark_profile%stop_stage('P_potential_update')
   
         !=========================================================================
         !             UPDATE FILES AND INFORMATION IN THE DIRECTORY
         !=========================================================================
         call g_scf_benchmark_profile%start_stage('P_scf_io')
         if (rank == 0) call update_fermi_in_input(this%en%fermi, this%control%fname)
         if (rank == 0) call save_state_scf(this%lattice%symbolic_atoms(:))
         call g_scf_benchmark_profile%stop_stage('P_scf_io')

         this%physical_total_energy = sum(this%symbolic_atom(:)%potential%etot)
         this%physical_total_energy = this%physical_total_energy + this%fixed_hubbard_energy_pout()
         call this%write_fixed_u_total_energy()
         call this%write_hubbard_moment_energy_diagnostics()
         call this%write_gbt_scf_diagnostics(i)
         total_energy = this%physical_total_energy
         magnetic_moment = sum(this%symbolic_atom(this%lattice%nbulk + 1:this%lattice%nbulk + this%lattice%nrec)%potential%mtot)
         call g_scf_benchmark_profile%finish_iteration(this%mix%delta, total_energy, this%en%fermi, magnetic_moment)

         !=========================================================================
         !                TEST IF THE CALCULATION IS CONVERGED
         !=========================================================================
         this%converged_iteration = i
         this%converged = this%is_converged(this%mix%delta)
         if (this%converged) then
            if (rank == 0) call g_logger%info('Converged!'//fmt('f12.10', this%mix%delta), __FILE__, __LINE__)
            if (this%hamiltonian%hubbard_u_sc_check .and. .not. this%bands%hubbard_u_converged) then
               if (rank == 0) call g_logger%info('SCF converged but hubbard_u_sc not converged; continuing.', __FILE__, __LINE__)
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

   !> @brief Perform one constraint-controller update without an SCF potential update.
   !> @details This is the inner loop of the corrected frozen-magnon MFT mode.
   !>          The ordinary LMTO potential is copied before the Hamiltonian/DOS
   !>          pass and restored afterwards.  Thus the density may respond to
   !>          the current constraining field, while QL/PL, radial potentials,
   !>          energies, and all other ordinary potential state remain fixed.
   subroutine run_fixed_potential_constraint_step(this, iteration)
      class(self), intent(inout) :: this
      integer, intent(in) :: iteration
      type(potential), allocatable :: frozen_potential(:)
      real(rp), allocatable :: reference_flat(:), current_flat(:)
      integer :: ia, pot_size
      real(rp) :: max_drift

      if (.not. this%control%constraints_enable) then
         call g_logger%fatal('run_fixed_potential_constraint_step requires constraints_enable = .true.', &
                            __FILE__, __LINE__)
      end if
      if (this%lattice%ntype <= 0) then
         call g_logger%fatal('run_fixed_potential_constraint_step requires at least one atom type.', &
                            __FILE__, __LINE__)
      end if

      allocate(frozen_potential(this%lattice%ntype))
      do ia = 1, this%lattice%ntype
         frozen_potential(ia) = this%symbolic_atom(ia)%potential
      end do
      pot_size = frozen_potential(1)%sizeof_potential_full()
      allocate(reference_flat(pot_size), current_flat(pot_size))
      reference_flat = 0.0_rp
      do ia = 1, this%lattice%ntype
         call frozen_potential(ia)%flatten_potential_full(current_flat)
         reference_flat = reference_flat + abs(current_flat)
      end do
      this%fixed_potential_checksum = sum(reference_flat)

      ! Deliberately bypass self%run(): it would mix the density and call
      ! run_scf(), which is precisely the q-dependent ordinary-potential drift
      ! that corrected MFT must exclude.
      call run_recursion(this, iteration)
      call run_dos(this, iteration=iteration)

      max_drift = 0.0_rp
      do ia = 1, this%lattice%ntype
         call this%symbolic_atom(ia)%potential%flatten_potential_full(current_flat)
         call frozen_potential(ia)%flatten_potential_full(reference_flat)
         max_drift = max(max_drift, maxval(abs(current_flat - reference_flat)))
      end do
      this%fixed_potential_transient_drift = max_drift

      ! Restore after the diagnostic comparison.  The controller state and the
      ! measured constraint_actual vectors intentionally survive this restore.
      do ia = 1, this%lattice%ntype
         this%symbolic_atom(ia)%potential = frozen_potential(ia)
      end do
      max_drift = 0.0_rp
      do ia = 1, this%lattice%ntype
         call this%symbolic_atom(ia)%potential%flatten_potential_full(current_flat)
         call frozen_potential(ia)%flatten_potential_full(reference_flat)
         max_drift = max(max_drift, maxval(abs(current_flat - reference_flat)))
      end do
      this%fixed_potential_max_drift = max_drift
      this%physical_total_energy = sum(frozen_potential(:)%etot)
      this%esumn = this%physical_total_energy
      this%esum = this%physical_total_energy

      deallocate(reference_flat, current_flat, frozen_potential)
   end subroutine run_fixed_potential_constraint_step

   !> @brief Reset the controller state between independent fixed-potential q passes.
   !> @param[in] zero_field When true (the default), start this q point with no
   !>        carried constraining field or controller history.  False preserves
   !>        the field and controller history for an explicit continuation run.
   subroutine reset_constraint_for_fixed_potential(this, zero_field)
      class(self), intent(inout) :: this
      logical, intent(in), optional :: zero_field
      logical :: zero_field_
      integer :: ia, plusbulk

      if (.not. this%control%constraints_enable) return
      zero_field_ = .true.
      if (present(zero_field)) zero_field_ = zero_field
      if (zero_field_) then
         do ia = 1, this%lattice%nrec
            plusbulk = this%lattice%nbulk + ia
            this%symbolic_atom(plusbulk)%mag_cfield = 0.0_rp
            this%symbolic_atom(plusbulk)%constraint_actual = 0.0_rp
         end do
         call initialize_cfd(this%lattice%nrec, 0)
         call initialize_cfd(this%lattice%nrec, 1, this%control%constraints_i_cons, &
                             this%control%constraints_code_prefac, this%control%constraints_diagnostics)
      end if
      this%constraint_metric = huge(1.0_rp)
      this%constraint_converged = .false.
      this%constraint_iteration = 0
      this%constraint_penalty_energy = 0.0_rp
      this%constraint_field_coupling_energy = 0.0_rp
      this%constraint_energy = 0.0_rp
   end subroutine reset_constraint_for_fixed_potential

   !> @brief Return a deterministic L1 checksum of all ordinary potential fields.
   !> @details The full flattening intentionally excludes the separate
   !>          constraining field, which lives on symbolic_atom%mag_cfield.
   function potential_checksum(this) result(checksum)
      class(self), intent(in) :: this
      real(rp) :: checksum
      real(rp), allocatable :: flat(:)
      integer :: ia, pot_size

      checksum = 0.0_rp
      if (this%lattice%ntype <= 0) return
      pot_size = this%symbolic_atom(1)%potential%sizeof_potential_full()
      allocate(flat(pot_size))
      do ia = 1, this%lattice%ntype
         call this%symbolic_atom(ia)%potential%flatten_potential_full(flat)
         checksum = checksum + sum(abs(flat))
      end do
      deallocate(flat)
   end function potential_checksum
   
   !=========================================================================
   !                   RUN RECURSION STEP
   !=========================================================================
   subroutine run_recursion(this, iter)
      class(self), intent(inout) :: this
      integer, intent(in) :: iter
      integer :: ia
   
      if (rank == 0) call g_logger%info('Perform recursion at step '//int2str(iter), __FILE__, __LINE__)
      call g_timer%start('recursion')

      if (this%use_kspace) then
         call g_scf_benchmark_profile%start_stage('P_hamiltonian_prepare')
      else
         call g_scf_benchmark_profile%start_stage('P_rs_hamiltonian_prepare')
      end if

      if (this%control%calctype == 'I' .or. this%control%calctype == 'S') then
         call this%hamiltonian%begin_rd1_vin_capture(iter, enable_rs_capture=.not. this%use_kspace)
      end if
   
      select case (this%control%calctype)
      case ('B')
         do ia = 1, this%lattice%nrec
            call this%symbolic_atom(ia)%build_pot()
         end do
         if (this%control%has_soc()) call this%hamiltonian%build_lsham()
         call this%hamiltonian%build_bulkham()
      case ('S')
         do ia = 1, this%lattice%ntype
            call this%symbolic_atom(ia)%build_pot()
         end do
         if (this%control%has_soc()) call this%hamiltonian%build_lsham()
         call this%hamiltonian%build_bulkham()
      case ('I')
         do ia = 1, this%lattice%ntype
            call this%symbolic_atom(ia)%build_pot()
         end do
         if (this%control%has_soc()) call this%hamiltonian%build_lsham()
         call this%hamiltonian%build_bulkham()
         call this%hamiltonian%build_locham()
      case ('L')
         ! B7.5: identical to ’S’. The loop runs to ntype, not nrec, because the
         ! frozen reference types need built potentials too -- for ’L’ those are
         ! both regions’ types (ntype = 2*nbulk_bulk + nlay, lattice_cluster.f90).
         ! No build_locham: nmax = 0 for ’L’ (lattice_lifecycle.f90), so the
         ! solver consumes ee, not hall.
         do ia = 1, this%lattice%ntype
            call this%symbolic_atom(ia)%build_pot()
         end do
         if (this%control%has_soc()) call this%hamiltonian%build_lsham()
         call this%hamiltonian%build_bulkham()
      end select

      if (this%use_kspace) then
         call g_scf_benchmark_profile%stop_stage('P_hamiltonian_prepare')
      else
         call g_scf_benchmark_profile%stop_stage('P_rs_hamiltonian_prepare')
      end if
   
      if (this%use_kspace) then
         if (rank == 0) call g_logger%info('run_recursion: use_kspace=.true., skipping recursion solver stage (Hamiltonian build kept).', __FILE__, __LINE__)
      else
         call g_scf_benchmark_profile%start_stage('P_rs_solver_kernel')
         select case (this%control%recur)
         case ('lanczos')
            call this%recursion%recur()
         case ('chebyshev')
            call this%recursion%chebyshev_recur()
         case ('block')
            call this%recursion%recur_b()
         end select
         call g_scf_benchmark_profile%stop_stage('P_rs_solver_kernel')
      end if
   
      call g_timer%stop('recursion')
   end subroutine run_recursion
   
   !=========================================================================
   !                   RUN DENSITY OF STATES CALCULATION
   !=========================================================================
   subroutine run_dos(this, m1_out, pl_diag_out, iteration)
      class(self), intent(inout) :: this
      complex(rp), intent(out), optional :: m1_out(:, :, :, :, :)
      real(rp), intent(out), optional :: pl_diag_out(:, :, :, :)
      integer, intent(in), optional :: iteration
      integer :: ia, l, lmax_site, li, iorb, plusbulk
      integer :: iteration_value, spin_check, l_check
      real(rp) :: q_up, q_dn, mz, mtot
      real(rp) :: reciprocal_h2d, reciprocal_solver, reciprocal_d2h, reciprocal_total
      real(rp), allocatable :: kspace_spin_mom(:,:)
      complex(rp), allocatable :: local_m1(:, :, :, :, :)
      complex(rp) :: p_d(5, 5), m1_d(5, 5), v_d(5, 5)
      real(rp) :: q0_d, center_d, q2_d, center_local, q2_local, t_local
      integer :: moment_status, site_type
      logical :: use_shifted_kmesh

      iteration_value = 0
      if (present(iteration)) iteration_value = iteration
   
      if (rank == 0) call g_logger%info('Calculating the density of states and the new moment bands', __FILE__, __LINE__)
      call g_timer%start('calculation-of-DOS')
   
      call this%en%e_mesh()

      if (this%use_kspace) then
         if (rank == 0) call g_logger%info('run_dos: use_kspace=.true. reciprocal SCF branch enabled (optional path).', __FILE__, __LINE__)
         if (.not. allocated(this%reciprocal_scf_cache)) then
            allocate(this%reciprocal_scf_cache)
            this%reciprocal_scf_cache = reciprocal(this%hamiltonian)
         end if
         if (rank == 0) then
            call g_logger%info('k-space SCF mesh: ' // &
                               int2str(this%reciprocal_scf_cache%nk_mesh(1)) // ' x ' // &
                               int2str(this%reciprocal_scf_cache%nk_mesh(2)) // ' x ' // &
                               int2str(this%reciprocal_scf_cache%nk_mesh(3)), __FILE__, __LINE__)
         end if
         if (.not. allocated(this%reciprocal_scf_cache%k_points)) then
            use_shifted_kmesh = (sum(abs(this%reciprocal_scf_cache%k_offset)) > 1.0e-12_rp)
            if (this%reciprocal_scf_cache%use_symmetry_reduction) then
               call this%reciprocal_scf_cache%generate_reduced_kpoint_mesh( &
                  this%reciprocal_scf_cache%nk_mesh, use_shifted_kmesh)
            else
               call this%reciprocal_scf_cache%generate_mp_mesh()
            end if
         end if
         call g_scf_benchmark_profile%start_stage('P_hk_assembly')
         call this%reciprocal_scf_cache%build_kspace_hamiltonian()
         call this%write_fixed_u_kspace_trace(.false.)
         call g_scf_benchmark_profile%stop_stage('P_hk_assembly')
         call g_timer%stop('calculation-of-DOS')
         call g_timer%start('diagonalization')
         call g_scf_benchmark_profile%start_stage('P_eigensolver')
         call this%reciprocal_scf_cache%diagonalize_hamiltonian()
         call g_scf_benchmark_profile%stop_stage('P_eigensolver')
         call g_timer%stop('diagonalization')

         reciprocal_h2d = 0.0_rp
         reciprocal_solver = 0.0_rp
         reciprocal_d2h = 0.0_rp
         reciprocal_total = 0.0_rp
         if (allocated(this%reciprocal_scf_cache%execution_backend)) then
            select type (backend => this%reciprocal_scf_cache%execution_backend)
            type is (cuda_reciprocal_backend)
               reciprocal_h2d = backend%h2d_seconds
               reciprocal_solver = backend%gpu_solve_seconds
               reciprocal_d2h = backend%d2h_seconds
               reciprocal_total = backend%total_reciprocal_seconds
            class default
               continue
            end select
         end if
         call g_scf_benchmark_profile%set_reciprocal_details( &
            reciprocal_h2d, reciprocal_solver, reciprocal_d2h, reciprocal_total)
         call g_timer%start('calculation-of-DOS')
         call this%reciprocal_scf_cache%calculate_density_of_states( &
            this%hamiltonian, &
            n_energy_points=this%en%channels_ldos + 10, &
            energy_range=[this%en%energy_min, this%en%energy_max], &
            fermi_level=this%en%fermi, &
            auto_find_fermi=.true.)
         ! Capture the same eigenpairs after the canonical Fermi level has
         ! been updated for this spectrum; Hamiltonian capture remains above.
         call this%write_fixed_u_kspace_trace(.true.)
         if (this%control%nsp == 1) then
            do ia = 1, this%lattice%nrec
               lmax_site = this%symbolic_atom(this%lattice%nbulk + ia)%potential%lmax
               do l_check = 0, lmax_site
                  do spin_check = 1, 2
                     q0_d = 0.0_rp
                     if (l_check < this%reciprocal_scf_cache%n_orb_types) then
                        q0_d = this%reciprocal_scf_cache%band_moments(ia, l_check + 1, spin_check, 1)
                     end if
                     call guard_q0_capacity(this%symbolic_atom(this%lattice%nbulk + ia)%element%symbol, &
                        ia, l_check, spin_check, iteration_value, q0_d)
                  end do
               end do
            end do
         end if
         if (this%control%uses_spinor_representation()) then
            allocate(kspace_spin_mom(3, this%lattice%nrec))
            call g_scf_benchmark_profile%start_stage('P_charge_spin_accumulate')
            call this%compute_kspace_spin_moments_spinor(this%reciprocal_scf_cache, kspace_spin_mom)
            call g_scf_benchmark_profile%stop_stage('P_charge_spin_accumulate')
         end if

         ! The opt-in local-moment route needs same-eigensystem P and M1
         ! before the band center/second moment are mapped to radial state.
         ! Keep the historical spectral ordering exactly as before.
         if (trim(this%hamiltonian%hubbard_radial_moments) == 'local' .and. &
             (this%hamiltonian%hubbard_u_general_check .or. this%hamiltonian%hubbard_test_operator_check)) then
            allocate(local_m1(this%lattice%nrec, lmax_basis + 1, 2, 2*lmax_basis + 1, 2*lmax_basis + 1))
            call this%reciprocal_scf_cache%calculate_ldm_from_projected_dos(this%lattice, &
               fixed_liechtenstein=.true., m1_out=local_m1)
            if (present(m1_out)) then
               m1_out(:, 1:size(local_m1, 2), 1:size(local_m1, 3), &
                  1:size(local_m1, 4), 1:size(local_m1, 5)) = local_m1
            end if
         end if

         ! Map reciprocal projected-DOS moments onto SCF quantities expected by mix/save_to.
         call g_scf_benchmark_profile%start_stage('P_density_build')
         do ia = 1, this%lattice%nrec
            lmax_site = this%symbolic_atom(this%lattice%nbulk + ia)%potential%lmax
            do l = 0, lmax_site
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%ql(:, l, :) = 0.0_rp
            end do

            do l = 0, min(lmax_site, this%reciprocal_scf_cache%n_orb_types - 1)
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%ql(1, l, 1) = this%reciprocal_scf_cache%band_moments(ia, l + 1, 1, 1)
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%ql(1, l, 2) = this%reciprocal_scf_cache%band_moments(ia, l + 1, 2, 1)
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%gravity_center(l + 1, 1) = &
                  this%reciprocal_scf_cache%band_moments(ia, l + 1, 1, 2) - this%symbolic_atom(this%lattice%nbulk + ia)%potential%vmad
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%gravity_center(l + 1, 2) = &
                  this%reciprocal_scf_cache%band_moments(ia, l + 1, 2, 2) - this%symbolic_atom(this%lattice%nbulk + ia)%potential%vmad
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%ql(3, l, 1) = &
                  this%reciprocal_scf_cache%band_moments(ia, l + 1, 1, 3)**2 * this%reciprocal_scf_cache%band_moments(ia, l + 1, 1, 1)
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%ql(3, l, 2) = &
                  this%reciprocal_scf_cache%band_moments(ia, l + 1, 2, 3)**2 * this%reciprocal_scf_cache%band_moments(ia, l + 1, 2, 1)
               if (l == 2 .and. allocated(local_m1)) then
                  site_type = this%lattice%ib(ia)
                  do li = 1, 2
                     q0_d = this%symbolic_atom(this%lattice%nbulk + ia)%potential%ql(1, 2, li)
                     center_d = this%reciprocal_scf_cache%band_moments(ia, 3, li, 2)
                     q2_d = this%symbolic_atom(this%lattice%nbulk + ia)%potential%ql(3, 2, li)
                     p_d = this%symbolic_atom(this%lattice%nbulk + ia)%potential%ldm(3, li, 1:5, 1:5)
                     m1_d = local_m1(ia, 3, li, 1:5, 1:5)
                     call extract_hubbard_d_shell_potential(this%hamiltonian%hubbard_u_pot, &
                        site_type, li, spin_off, v_d, moment_status)
                     if (moment_status /= 0) call g_logger%fatal( &
                        'Unable to map local d-shell Hubbard potential in reciprocal path.', __FILE__, __LINE__)
                     if (maxval(abs(v_d)) == 0.0_rp) cycle
                     call correct_d_shell_moments(p_d, m1_d, v_d, q0_d, center_d, q2_d, &
                        center_local, q2_local, t_local, moment_status)
                     if (moment_status /= 0) call g_logger%fatal( &
                        'Invalid local d-shell Hubbard moment correction in reciprocal path; status='// &
                        fmt('i0', moment_status), __FILE__, __LINE__)
                     this%symbolic_atom(this%lattice%nbulk + ia)%potential%gravity_center(3, li) = &
                        center_local - this%symbolic_atom(this%lattice%nbulk + ia)%potential%vmad
                     this%symbolic_atom(this%lattice%nbulk + ia)%potential%ql(3, 2, li) = q2_local
                  end do
               end if
            end do

            if (this%control%uses_spinor_representation()) then
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%mx = kspace_spin_mom(1, ia)
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%my = kspace_spin_mom(2, ia)
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%mz = kspace_spin_mom(3, ia)
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%mom0(:) = kspace_spin_mom(:, ia)
               mtot = sqrt(sum(kspace_spin_mom(:, ia)**2)) + 1.0e-15_rp
            else
               ! Scalar-relativistic collinear proxy from spin occupations;
               ! this is a global-mode branch, not an XC channel count.
               q_up = sum(this%symbolic_atom(this%lattice%nbulk + ia)%potential%ql(1, 0:lmax_site, 1))
               q_dn = sum(this%symbolic_atom(this%lattice%nbulk + ia)%potential%ql(1, 0:lmax_site, 2))
               mz = q_up - q_dn
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%mx = 0.0_rp
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%my = 0.0_rp
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%mz = mz
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%mom0(:) = [0.0_rp, 0.0_rp, mz]
               mtot = abs(mz)
            end if
            this%symbolic_atom(this%lattice%nbulk + ia)%potential%mtot = mtot
            if (mtot > tiny(1.0_rp)) then
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%mom(:) = &
                  [this%symbolic_atom(this%lattice%nbulk + ia)%potential%mx, &
                   this%symbolic_atom(this%lattice%nbulk + ia)%potential%my, &
                   this%symbolic_atom(this%lattice%nbulk + ia)%potential%mz]/mtot
            else
               this%symbolic_atom(this%lattice%nbulk + ia)%potential%mom(:) = [0.0_rp, 0.0_rp, 1.0_rp]
            end if
            if (this%control%constraints_enable) then
               this%symbolic_atom(this%lattice%nbulk + ia)%constraint_actual(:) = &
                  [this%symbolic_atom(this%lattice%nbulk + ia)%potential%mx, &
                   this%symbolic_atom(this%lattice%nbulk + ia)%potential%my, &
                   this%symbolic_atom(this%lattice%nbulk + ia)%potential%mz]
            end if
            call this%symbolic_atom(this%lattice%nbulk + ia)%potential%copy_mom_to_scal()
            if (rank == 0) then
               call g_logger%info('Spin moment of atom'//fmt('i4', ia)//' is '// &
                                  fmt('f10.6', this%symbolic_atom(this%lattice%nbulk + ia)%potential%mtot), __FILE__, __LINE__)
               call g_logger%info('Spin moment projections of atom'//fmt('i4', ia)//' is '// &
                                  fmt('f10.6', this%symbolic_atom(this%lattice%nbulk + ia)%potential%mx)//' '// &
                                  fmt('f10.6', this%symbolic_atom(this%lattice%nbulk + ia)%potential%my)//' '// &
                                  fmt('f10.6', this%symbolic_atom(this%lattice%nbulk + ia)%potential%mz), __FILE__, __LINE__)
            end if
         end do
         if (allocated(kspace_spin_mom)) deallocate(kspace_spin_mom)

                        do ia = 1, this%lattice%nrec
                           this%mix%mag_new(ia, :) = this%symbolic_atom(this%lattice%nbulk + ia)%potential%mom(:)
                        end do

                        call this%mix_magnetic_state()

                        do ia = 1, this%lattice%nrec
                           this%symbolic_atom(this%lattice%nbulk + ia)%potential%mom(:) = this%mix%mag_mix(ia, :)
                        end do

         call this%apply_constraints()
         if (present(pl_diag_out)) then
            call this%bands%calculate_pl(pl_diag_out, iteration_value)
         else
            call this%bands%calculate_pl(iteration=iteration_value)
         end if
         this%en%fermi = this%reciprocal_scf_cache%fermi_level
         call g_scf_benchmark_profile%stop_stage('P_density_build')

         if ((this%hamiltonian%hubbard_u_general_check .or. this%hamiltonian%hubbard_u_sc_check .or. &
              this%hamiltonian%hubbard_v_check .or. this%hamiltonian%hubbard_test_operator_check) &
             .and. .not. allocated(local_m1)) then
            call this%reciprocal_scf_cache%calculate_ldm_from_projected_dos(this%lattice, &
               fixed_liechtenstein=(this%hamiltonian%hubbard_u_general_check .or. &
                  this%hamiltonian%hubbard_test_operator_check) .and. &
                  trim(this%hamiltonian%hubbard_u_potential_form) == 'liechtenstein', m1_out=m1_out)
            if (rank == 0) then
               do ia = 1, this%lattice%nrec
                  plusbulk = this%lattice%nbulk + ia
                  do li = 0, min(3, lmax_basis)
                     call g_logger%info('HUBBARD_BANDCMP atom='//fmt('i4', ia)//' l='//fmt('i2', li)// &
                                        ' ql_up='//fmt('f10.6', this%symbolic_atom(plusbulk)%potential%ql(1, li, 1))// &
                                        ' ql_dn='//fmt('f10.6', this%symbolic_atom(plusbulk)%potential%ql(1, li, 2))// &
                                        ' ldm_up='//fmt('f10.6', sum([(real(this%symbolic_atom(plusbulk)%potential%ldm(li + 1, 1, iorb, iorb), rp), iorb=1,2*li+1)]))// &
                                        ' ldm_dn='//fmt('f10.6', sum([(real(this%symbolic_atom(plusbulk)%potential%ldm(li + 1, 2, iorb, iorb), rp), iorb=1,2*li+1)])), __FILE__, __LINE__)
                  end do
               end do
            end if
         else if (present(m1_out)) then
            ! The no-U M0 reference still needs M1 from the canonical current
            ! eigenpairs, but must not populate/serialize stale LDM state.
            call this%reciprocal_scf_cache%calculate_ldm_from_projected_dos( &
               this%lattice, m1_out=m1_out, m1_only=.true.)
         end if
         call this%write_kspace_scf_dos_outputs(this%reciprocal_scf_cache)

         ! Canonical k-space EBAND is the projection-free occupied-eigenvalue
         ! sum computed with the same occupations that set EF/electron count.
         ! Projected m0*m1 remains a DOS diagnostic inside reciprocal_dos.
         this%bands%eband = this%reciprocal_scf_cache%canonical_band_energy

         call this%mix%save_to('new', iteration=iteration_value, guard_collinear=this%control%nsp == 1)
         call g_timer%stop('calculation-of-DOS')
         return
      end if
   
      select case (this%control%recur)
      case ('lanczos')
         call g_scf_benchmark_profile%start_stage('P_rs_green_function')
         call this%green%sgreen()
         call g_scf_benchmark_profile%stop_stage('P_rs_green_function')
      case ('chebyshev')
         call g_scf_benchmark_profile%start_stage('P_rs_spectral_reconstruct')
         call this%green%chebyshev_dos_dispatch()  ! Dispatches to GPU/C++/legacy based on control flags
         call g_scf_benchmark_profile%stop_stage('P_rs_spectral_reconstruct')
      case ('block')
         call g_scf_benchmark_profile%start_stage('P_rs_green_function')
         call this%recursion%zsqr()
         call this%green%block_green()
         call g_scf_benchmark_profile%stop_stage('P_rs_green_function')
      end select
   
      call g_scf_benchmark_profile%start_stage('P_rs_fermi')
      call this%bands%calculate_fermi() ! Calculate the Fermi energy
      if (this%hamiltonian%hubbard_u_sc_check) then
         call this%bands%calculate_hubbard_u_sc()
      end if
      call g_scf_benchmark_profile%stop_stage('P_rs_fermi')
      !=========================================================================
      !  MIX THE MAGNETIC MOMENTS BEFORE CALCULATING THE NEW BAND MOMENTS QL
      !=========================================================================
      call g_scf_benchmark_profile%start_stage('P_rs_charge_spin_accumulate')
      call this%bands%calculate_magnetic_moments()
      call g_scf_benchmark_profile%stop_stage('P_rs_charge_spin_accumulate')
   
      do ia = 1, this%lattice%nrec
         this%mix%mag_new(ia, :) = this%symbolic_atom(this%lattice%nbulk + ia)%potential%mom(:)
      end do
   
      call this%mix_magnetic_state()
   
      do ia = 1, this%lattice%nrec
         this%symbolic_atom(this%lattice%nbulk + ia)%potential%mom(:) = this%mix%mag_mix(ia, :)
      end do

      call this%apply_constraints()
   
      !=========================================================================
      !                  CALCULATE THE NEW BAND MOMENTS QL
      !=========================================================================
      call g_scf_benchmark_profile%start_stage('P_rs_density_build')
      if (present(m1_out) .or. present(pl_diag_out)) then
         call this%bands%calculate_moments(m1_out=m1_out, pl_diag_out=pl_diag_out, iteration=iteration_value)
      else
         call this%bands%calculate_moments(iteration=iteration_value)
      end if
      call g_scf_benchmark_profile%stop_stage('P_rs_density_build')
      call this%mix%save_to('new', iteration=iteration_value, guard_collinear=this%control%nsp == 1)

      !=========================================================================
      !        CALCULATE l=1 DIPOLE CHARGE MOMENTS (B6, surface only)
      !=========================================================================
      ! Uses the just-computed on-site density matrix (g0) and the radial
      ! partial waves exported by NEWRHO in the previous SCF iteration
      ! (standard SCF lag; identical at convergence). Result -> potential%q10,
      ! consumed by charge%surfpot below. Bulk-limit / disabled => q10 = 0.
      if (this%control%dipole_electrostatics) then
         call compute_dipole_moments(this%bands, this%phi_amp, this%control%dipole_mix)
      end if

      call g_timer%stop('calculation-of-DOS')
   end subroutine run_dos

   !=========================================================================
   !                   RUN SELF-CONSISTENT FIELD UPDATE
   !=========================================================================
   subroutine run_scf(this)
      class(self), intent(inout) :: this
      real(rp), dimension(6) :: QSL
      real(rp), dimension(:, :), allocatable :: T_comm
      integer :: ia, na_glob, pot_size
   
      call g_timer%start('atomic-scf')
   
      !=========================================================================
      !                       MAKE SFC ATOMIC SPHERE
      !=========================================================================
      do ia = start_atom, end_atom
         this%fixed_u_trace_current_site = ia
         qsl = this%lmtst(this%symbolic_atom(this%lattice%nbulk + ia), ia) ! Makes the atomic sphere self-consistent and caltulate the orthogonal pottential parameters
      end do
      if (rank == 0) call g_logger%info('Atomic SFC done for all atoms', __FILE__, __LINE__)

      !=========================================================================
      !      TRANSFER ATOMIC POTENTIAL DATA ACROSS MPI RANKS
      !=========================================================================
#ifdef USE_MPI
      ! Use a rank-independent reference atom for buffer sizing so idle MPI
      ! ranks (numprocs > nrec) never dereference an invalid start_atom.
      pot_size = this%symbolic_atom(this%lattice%nbulk + 1)%potential%sizeof_potential_full()
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
         if (rank == 0) call g_logger%info('From orthogonal to TB basis for atom '//this%symbolic_atom(this%lattice%nbulk + ia)%element%symbol, __FILE__, __LINE__)
         call this%symbolic_atom(this%lattice%nbulk + ia)%predls(this%lattice%wav*ang2au)
      end do
   
      call g_timer%stop('atomic-scf')
   end subroutine run_scf

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
      integer :: nb_slice, ncols
      real(rp) :: denom
      logical :: reciprocal_scf_active
      real(rp), dimension(this%lattice%nrec, 3) :: magmom, lmom
      real(rp), dimension(3, this%lattice%nrec) :: mag_for
      ! The reciprocal SCF branch has no real-space Green function.  Its
      ! legacy projected-DOS refresh is therefore both unnecessary and, under
      ! MPI, can dereference rank-local recursion storage that was never
      ! populated.  The reciprocal branch already supplies the band energy and
      ! site moments; retain the report’s zero-field convention instead.  The
      ! cache is the authoritative post-run marker too: it remains allocated
      ! through reporting even if a copied/stale self%use_kspace flag is not.
      reciprocal_scf_active = this%use_kspace .or. allocated(this%reciprocal_scf_cache)
      if (reciprocal_scf_active) then
         this%bands%mag_for(:, :) = 0.0_rp
      else
         call this%bands%calculate_magnetic_torques()
      end if

      ! Transfer values across ranks
      magmom = 0.0d0
      mag_for = 0.0d0
      lmom = 0.0d0
      do ia = start_atom, end_atom
         ia_loc = g2l_map(ia)
         ! `potential%mom0` is the integrated moment in the solver’s active
         ! spin basis.  For gbt_single_q that basis is the primitive rotating
         ! frame; lab-cell phases are reconstructed only by an explicit GBT
         ! frame transformation for output/comparison.
         magmom(ia, :) = this%symbolic_atom(this%lattice%nbulk + ia)%potential%mom0(:)
         lmom(ia, :) = this%symbolic_atom(this%lattice%nbulk + ia)%potential%lmom(:)
         mag_for(:, ia) = this%bands%mag_for(:, ia_loc)        
      end do

#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, magmom, product(shape(magmom)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, mag_for, product(shape(mag_for)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, lmom, product(shape(lmom)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      if (rank == 0) then
         ! These are shared files. Root owns their complete open/write/close
         ! lifetime; reductions above deliberately remain collective.
         open (newunit=newunit, file='report.out', action='write', iostat=iostatus, status='replace')
         open (unit=10, file='minfo.out', action='write', iostat=iostatus, status='replace')
         open (unit=20, file='linfo.out', action='write', iostat=iostatus, status='replace')
         call g_logger%info('Writing calculation report to report.out', __FILE__, __LINE__)
         !===========================================================================
         !                      Total Energy
         !===========================================================================
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(A)') '|                       Total Energy                                      |'
         write (newunit, '(A)') '==========================================================================='
!         write (newunit, ’(a,f20.10)’) ’Total energy of system: ’, sum(this%symbolic_atom(this%lattice%nbulk+1:this%lattice%ntype)%potential%etot)
         this%physical_total_energy = sum(this%symbolic_atom(:)%potential%etot)
         this%physical_total_energy = this%physical_total_energy + this%fixed_hubbard_energy_pout()
         write (newunit, '(a,f20.10)') 'Physical DFT total energy of system: ', this%physical_total_energy
         ! Keep the historical key for downstream parsers, but make the
         ! physical meaning explicit through the dedicated line above.
         write (newunit, '(a,f20.10)') 'Total energy of system: ', this%physical_total_energy
         write (newunit, '(a,f20.10)') 'Constraint penalty/controller energy (diagnostic only): ', this%constraint_penalty_energy
         write (newunit, '(a,f20.10)') 'Constraint-field coupling <V_con> (band diagnostic only): ', this%constraint_field_coupling_energy
         write (newunit, '(a,f20.10)') 'Constraint energy:      ', this%constraint_energy
         !===========================================================================
         !                       Band Energy
         !===========================================================================
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(A)') '|                       Band Energy                                       |'
         write (newunit, '(A)') '==========================================================================='
         if (associated(this%bands)) then
            if (.not. reciprocal_scf_active) call this%bands%calculate_band_energy()
            write (newunit, '(a,f16.10)') 'Band energy of system: ', this%bands%eband
         else
            write (newunit, '(a)') 'Band energy of system: unavailable (bands object not associated)'
         end if
         !===========================================================================
         !                       Magnetization
         !===========================================================================
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(A)') '|                     Spin moment                                         |'
         write (newunit, '(A)') '==========================================================================='
         if (trim(this%hamiltonian%magnetic_representation) == gbt_single_q) then
            write (newunit, '(A)') 'Spin-moment frame: rotating primitive frame (GBT; lab phases not applied)'
         else
            write (newunit, '(A)') 'Spin-moment frame: active calculation frame'
         end if
         write (newunit, '(a,3f16.10)') 'Total spin moment: ', sum(magmom(1:this%lattice%nrec, 1:3), dim=1)
         do ia = 1, this%lattice%nrec
            write (newunit, '(a,i4,a,f10.6)') 'Spin moment of atom', ia, ':', norm2(magmom(ia, :))
            write (newunit, '(a,i4,a,3f10.6)') 'Spin moment projections of atom', ia, ':', magmom(ia, :)
            write (newunit, '(a,i4,a,3f16.6)') 'Magnetic force on atom', ia, ':', mag_for(:, ia)
         end do
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(A)') '|                     Orbital moment                                      |'
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(a,3f16.10)') 'Total orbital moment: ', sum(lmom(1:this%lattice%nrec, 1:3), dim=1)
         do ia = 1, this%lattice%nrec
            write (newunit, '(a,i4,a,f10.6)') 'Orbital moment of atom', ia, ':', norm2(lmom(ia, :))
            write (newunit, '(a,i4,a,3f10.6)') 'Orbital moment projections of atom', ia, ':', lmom(ia, :)
         end do
         !===========================================================================
         !                           Charge
         !===========================================================================
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(A)') '|                       Charge Transfer                                   |'
         write (newunit, '(A)') '==========================================================================='
         do ia = 1, this%lattice%nrec
            nb_slice = min(this%symbolic_atom(this%lattice%nbulk + ia)%potential%lmax, lmax_basis) + 1
            write (newunit, '(a,i4,a,f10.6)') 'Occupation at atom', ia, ':', sum(this%mix%qia(ia, 1:2*nb_slice))
            write (newunit, '(a,i4,a,*(f10.6,1x))') 'Up orbital occupation at atom', ia, ':', this%mix%qia(ia, 1:nb_slice)
            write (newunit, '(a,i4,a,*(f10.6,1x))') 'Down orbital occupation at atom', ia, ':', this%mix%qia(ia, nb_slice+1:2*nb_slice)
            write (newunit, '(a,i4,a,f10.6)') 'Charge transfer at atom', ia, ':', this%charge%dq(ia)
         end do
         !===========================================================================
         !                           Fermi Energy
         !===========================================================================
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(A)') '|                       Fermi Energy                                      |'
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(a,f10.6)') 'Fermi energy: ', this%en%fermi
         
         do ia = 1, this%lattice%nrec
            write (10, '(a,i4,a,3f10.6)') 'Spin moment direction of atom', ia, ':', (magmom(ia, :))/MAX(norm2(magmom(ia, :)), TINY(1.0_rp))
            write (20, '(a,i4,a,3f10.6)') 'Orbital moment direction of atom', ia, ':', (lmom(ia, :))/MAX(norm2(lmom(ia, :)), TINY(1.0_rp))
         end do
         !===========================================================================
         !                           Log info     
         !===========================================================================
         write (newunit, '(A)') '==========================================================================='
         write (newunit, '(A)') '|                             Log info                                     |'
         write (newunit, '(A)') '==========================================================================='
         ncols = size(this%mix%qia, 2)
         nb_slice = min(6*(lmax_basis + 1), ncols)
         denom = max(1.0_rp, real(nb_slice)/2.0_rp)
         write (newunit, '(a,f10.6)') 'Total RMS Diff: ', this%mix%delta
         do ia = 1, this%lattice%nrec
            write (newunit, '(a,i4,a,f10.6)') 'RMS Diff of atom', ia, ':', sqrt(sum((this%mix%qia_old(ia, 1:nb_slice) - this%mix%qia_new(ia, 1:nb_slice))**2))/denom
         end do
         close(newunit)
         close(10)
         close(20)
         ! Persist the exact radial ground-state contract for audit and for
         ! consumers running after the atomic temporaries have expired.  The
         ! in-memory copy on symbolic_atom is the primary accessor; these
         ! files are a small, machine-readable evidence seam.
         do ia = 1, this%lattice%nrec
            call this%symbolic_atom(this%lattice%nbulk + ia)%radial_ground_state%set_reported_moment(magmom(ia, :))
            call this%symbolic_atom(this%lattice%nbulk + ia)%radial_ground_state%write_file( &
               'radial_ground_state_'//trim(this%symbolic_atom(this%lattice%nbulk + ia)%element%symbol)//'_'// &
               trim(fmt('i0', ia))//'.dat', ia)
         end do
      end if

      ! Print angle betweens magnetic and orbital moments
      call this%bands%calculate_angles(magmom, lmom)

      ! Print hyperfine structure
      if (this%control%hyperfine) then
         do ia = 1, this%lattice%nrec
            call this%symbolic_atom(ia)%potential%print_hyperfine(ia)
         end do
      end if

      if (rank == 0) call g_logger%info('Calculation finished. Report printed in report.out', __FILE__, __LINE__)
      
   end subroutine report

   function is_converged(this, delta_en) result(l)
      class(self), intent(inout) :: this
      real(rp), intent(in) :: delta_en
      logical :: l

      l = delta_en < this%conv_thr
      if (this%fixed_u_residual_active) l = l .and. this%fixed_u_residual_converged
      if (this%control%constraints_enable) then
         l = l .and. this%constraint_converged
         if (rank == 0 .and. .not. this%constraint_converged) then
            call g_logger%info('Constraint convergence pending: RMS-angle='// &
                               fmt('es12.4', this%constraint_metric)//' tolerance='// &
                               fmt('es12.4', this%constraint_tolerance), __FILE__, __LINE__)
         end if
      end if
   end function is_converged

   subroutine update_fermi_in_input(fermi_value, filename)
       use iso_fortran_env, only: stderr => error_unit
       real(rp), intent(in) :: fermi_value
       character(len=*), intent(in) :: filename
       character(len=256), allocatable :: lines(:)
       integer :: i, nlines, iunit, ounit, stat, eq_pos, com_pos
       logical :: in_energy, fermi_found
       character(len=256) :: temp_filename, key_part, rest, comment_part, new_rest
   
       ! Read the entire input file into memory
       open(newunit=iunit, file=filename, status='old', action='read', iostat=stat)
       if (stat /= 0) then
           write(stderr, '(a)') "Error: Failed to open input file."
           return
       end if
   
       ! Count lines
       nlines = 0
       do
           read(iunit, '(a)', iostat=stat) 
           if (stat /= 0) exit
           nlines = nlines + 1
       end do
       rewind(iunit)
   
       ! Read lines into buffer
       allocate(lines(nlines))
       do i = 1, nlines
           read(iunit, '(a)') lines(i)
       end do
       close(iunit)
   
       ! Find and update the fermi line in the &energy section
       in_energy = .false.
       fermi_found = .false.
       do i = 1, nlines
           lines(i) = adjustl(lines(i))
           if (lines(i) == '&energy') then
               in_energy = .true.
           else if (lines(i) == '/') then
               in_energy = .false.
           else if (in_energy .and. .not. fermi_found) then
               if (index(lines(i), 'fermi') == 1) then
                   ! Split line into key, value, and comment
                   eq_pos = index(lines(i), '=')
                   if (eq_pos == 0) cycle
                   key_part = lines(i)(1:eq_pos)
                   rest = adjustl(lines(i)(eq_pos+1:))
                   com_pos = index(rest, '!')
                   if (com_pos > 0) then
                       comment_part = rest(com_pos:)
                       rest = rest(:com_pos-1)
                   else
                       comment_part = ''
                   end if
   
                   ! Replace the value with the new Fermi energy
                   write(new_rest, '(f12.6)') fermi_value  ! Format to 6 decimal places
                   lines(i) = trim(key_part) // ' ' // trim(adjustl(new_rest)) // ' ' // trim(comment_part)
                   fermi_found = .true.
               end if
           end if
       end do
   
       ! Write the updated content to a temporary file
       temp_filename = trim(filename) // '.tmp'
       open(newunit=ounit, file=temp_filename, status='replace', action='write', iostat=stat)
       if (stat /= 0) then
           write(stderr, '(a)') "Error: Failed to create temporary file."
           return
       end if
   
       do i = 1, nlines
           write(ounit, '(a)') trim(lines(i))
       end do
       close(ounit)
   
       ! Replace the original file with the temporary file
       call rename(temp_filename, filename)
   
       deallocate(lines)
   end subroutine update_fermi_in_input
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
   function lmtst(this, atom, isite) result(QSL)
      implicit none
      class(self), intent(inout) :: this
      class(symbolic_atom), intent(inout) :: atom
      integer, intent(in), optional :: isite
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
      ! Local radial/XC channel count; independent of control%nsp.
      NSP = RADIAL_XC_SPIN_CHANNELS

      QSL(:) = 0.0_rp
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
      real(rp) :: B_fsm, B, deg, DFCORE, AMGM, EA, RPB, DL, DRHO_CONTROL, DRHO_INTEGRATED, TOL, TOLRSQ, BETA, VHRMAX, VSUM, TL, BETA1, VNUCL, SUM, RHO0T, WGT, DRDI, RHOMU, RHOVH, ZVNUCL, OB4PI
      integer :: ISP, ncore, nval, l, n_radial_spin_channels, lmax, konf, qn_default, IFCORE, LCORE, KONFIG, IPR, NR, IR, ITER, NITER, IPR1, II
      integer :: energy_trace_unit, energy_trace_ios
      logical :: LAST

      real(rp), dimension(NCMX) :: EC
      real(rp), dimension(NCMX) :: EC_SEED
      real(rp), dimension(NVMX) :: EV
      real(rp), dimension(NVMX) :: EV_SEED
      real(rp), dimension(2) :: VRMAX
      real(rp), dimension(:, :), allocatable :: rho_in, rho
      real(rp), dimension(:, :), allocatable :: v_trace_base
      real(rp), dimension(:), allocatable :: rofi_trace_base
      real(rp), dimension(2) :: qval

      ipr = 0
      ! This is a local radial/XC channel count.  It intentionally remains
      ! two for a scalar-relativistic collinear global mode (control%nsp=1).
      n_radial_spin_channels = RADIAL_XC_SPIN_CHANNELS
      lmax = atom%potential%lmax
      rho_in = atom%rho0(n_radial_spin_channels)

      allocate (v, mold=rho_in)
      allocate (rho, mold=rho_in)
      allocate (rofi(size(rho_in(:, 1))))
      if (allocated(this%fun2)) this%fun2(:, :, :) = 0.0_rp
      if (allocated(this%phi_amp)) this%phi_amp(:, :, :) = 0.0_rp

      nr = size(rofi)
      B_fsm = merge(-atom%mag_cfield(3), real(0.0, rp), this%lattice%control%do_comom)
      if (this%magnetic_seed_active) B_fsm = B_fsm + this%magnetic_seed_field
      xc_obj = xc(this%lattice%control)
      b = atom%B()

      call g_logger%info('XC functional: '//trim(xc_obj%TXCH)//' '//int2str(xc_obj%txc)//real2str(B), __FILE__, __LINE__)

      lmax = atom%potential%lmax
      qval(:) = 0
      atom%qc = 0
      ncore = 0
      nval = 0
      vhrmax = 0.d0
      tol = 1.d-6
      tolrsq = 1.d-8
      beta = 0.3d0
      call atomsc_inner_controls(this%atomsc_test_tol, this%atomsc_test_tolrsq, &
         this%atomsc_test_inner_cap, tol, tolrsq, niter)

      ! Backward-compatible guard: legacy spd potential files promoted to lmax=3
      ! may not provide valid f-channel PL. Seed invalid channels from element
      ! quantum numbers, otherwise enforce the minimum admissible n=l+1.
      do isp = 1, n_radial_spin_channels
         do l = 0, lmax
            if (atom%potential%pl(l, isp) < real(l + 1, rp)) then
               select case (l)
               case (0)
                  qn_default = atom%element%num_quant_s
               case (1)
                  qn_default = atom%element%num_quant_p
               case (2)
                  qn_default = atom%element%num_quant_d
               case (3)
                  qn_default = atom%element%num_quant_f
               case default
                  qn_default = l + 1
               end select
               if (qn_default < l + 1) qn_default = l + 1
               atom%potential%pl(l, isp) = real(qn_default, rp)
            end if
            ! For promoted high-l channels (e.g. f), enforce a minimal non-zero
            ! Q0 to avoid singular atomic initialization from legacy inputs.
            if (l >= 3) then
               if (abs(atom%potential%ql(1, l, isp)) <= tiny(1.0_rp)) then
                  atom%potential%ql(1, l, isp) = 1.0e-3_rp
               end if
            end if
         end do
      end do

      do l = 0, lmax
         deg = (2*(2*l + 1))/n_radial_spin_channels
         do isp = 1, n_radial_spin_channels
            if (atom%potential%pl(l, isp) /= atom%potential%pl(l, isp)) then
               call g_logger%fatal('Invalid PL before atomic core/valence setup: atom='//trim(atom%element%symbol)// &
                                   ' lmax='//int2str(lmax)//' l='//int2str(l)//' spin='//int2str(isp)// &
                                   ' PL=NaN', __FILE__, __LINE__)
            end if
            konfig = int(atom%potential%pl(l, isp))
            if (konfig < l + 1) then
               call g_logger%fatal('Invalid PL before atomic core/valence setup: atom='//trim(atom%element%symbol)// &
                                   ' lmax='//int2str(lmax)//' l='//int2str(l)//' spin='//int2str(isp)// &
                                   ' PL='//real2str(atom%potential%pl(l, isp))//' KONFIG='//int2str(konfig), __FILE__, __LINE__)
            end if
            do konf = l + 1, konfig - 1
               ncore = ncore + 1
               ec(ncore) = -5.d0
               atom%qc = atom%qc + deg
            end do
            nval = nval + 1
            ev(nval) = -0.5d0
            qval(isp) = qval(isp) + atom%potential%ql(1, l, isp)
         end do
         ! Debug output for orbital initialization
         !if (lmax >= 3) then
         !   write (6, ’(A, I1, A, F8.4, F8.4, A, I2, A, 2(F10.6,2X))’) &
         !      ″  L=″, l, ″  PL=″, atom%potential%pl(l, 1), atom%potential%pl(l, 2), &
         !      ″  KONFIG=″, int(atom%potential%pl(l, 1)), &
         !      ″  Q0(spin1,2)=″, atom%potential%ql(1, l, 1), atom%potential%ql(1, l, 2)
         !end if
      end do
      !----MODIFICATIONS TO INCLUDE FRACTIONARY F OCCUPATION IN THE CORE------
      IFCORE = atom%element%f_core
      DFCORE = REAL(IFCORE)
      call g_logger%info('F core check:'//int2str(ifcore), __FILE__, __LINE__)
      !write (9, *) ´F-core check:´, IFCORE, DFCORE
      ! Only apply f-core treatment when lmax < 3 (f is not in valence basis)
      if (IFCORE /= 0 .and. lmax < 3) then
         LCORE = 3
         DEG = (2*(2*LCORE + 1))/n_radial_spin_channels
         do ISP = 1, n_radial_spin_channels
            !------DEFINES DEGENERACY OF UP AND DW F-CORE--------------
            !------MAKES DEG=2*DEG se NAO MAGNETICO!!!-----------------
            if (n_radial_spin_channels == 1) then
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
      if (n_radial_spin_channels == 2) then
         amgm = qval(1) - qval(2)
      end if
      ! Print statement, if ipr different than 0. Obsolete, therefore commented.
!    if (IPR >= 1) then
!       write (6, 10001) atom%element%atomic_number, atom%potential%ws_r, atom%a, NR, JOB, NCORE, atom%QC, QVAL, atom%DQ, AMGM
!       do ISP = 1, NSP
!          write (6, *) ″ ″
!          if (NSP == 2) then
!             write (6, *) ″SPIN″, ISP
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
         call this%newrho(atom, atom%element%atomic_number, lmax, atom%a, b, nr, rofi, v, rho_in, atom%potential%pl, atom%potential%ql, SEC, SEV, EC, EV, TOLRSQ, n_radial_spin_channels, 0)
      elseif (JOB /= "RHO") then
         stop "*** ATOMSC EXPECTS JOB=RHO OR JOB=POT"
      end if
      ! -------- START SELF-CONSISTENCY LOOP ------
      ! DRHO_CONTROL is the historical mesh-coordinate L1 metric.  Keep it
      ! on the old nonlinear-control path; DRHO_INTEGRATED is an observable
      ! physical electron-number residual and must not change that path.
      DRHO_CONTROL = 100.d0
      LAST = .false.
      do ITER = 1, NITER
         TL = atomsc_solver_tolerance(ITER, DRHO_CONTROL, TOLRSQ)
         BETA1 = atomsc_mixing_beta(ITER, DRHO_CONTROL, BETA)
       !!    IF(MOD(ITER, 3).EQ.2.AND.DRHO_CONTROL.LT.0.1D0) BETA1=1.0D0
         IPR1 = 0
         if (LAST) then
            IPR1 = IPR
         end if
         if (LAST .and. this%fixed_u_trace_enabled .and. rank == 0 .and. &
             (this%fixed_u_trace_selected_iteration <= 0 .or. &
              this%fixed_u_trace_iteration == this%fixed_u_trace_selected_iteration)) then
            if (allocated(v_trace_base)) deallocate(v_trace_base)
            if (allocated(rofi_trace_base)) deallocate(rofi_trace_base)
            allocate(v_trace_base, mold=v)
            allocate(rofi_trace_base, mold=rofi)
            v_trace_base = v
            rofi_trace_base = rofi
         end if
         call POISS0(atom%element%atomic_number, atom%a, b, rofi, rho_in, NR, VHRMAX, V, RVH, VSUM, n_radial_spin_channels)
         VNUCL = V(1, 1)
         !call VXC0SP_old(atom%element%atomic_number, atom%a, B, rofi, rho_in, NR, V, RHO0, REPS, RMU, NSP)
         if (LAST) then
            ! LAST is entered only after the preceding density was mixed and
            ! the final accepted radial residual has been tested.  Capture
            ! this VXC evaluation, not an initial or pre-mixing state.
            call this%VXC0SP(xc_obj, atom%element%atomic_number, atom%a, b, rofi, rho_in, NR, V, RHO0, REPS, RMU, &
                             n_radial_spin_channels, B_fsm, atom%radial_ground_state)
         else
            call this%VXC0SP(xc_obj, atom%element%atomic_number, atom%a, b, rofi, rho_in, NR, V, RHO0, REPS, RMU, &
                             n_radial_spin_channels, B_fsm)
         end if
         this%fixed_u_trace_capture_radial_state = LAST .and. this%fixed_u_trace_enabled .and. &
            rank == 0 .and. .not. this%fixed_u_trace_moments_only .and. &
            .not. this%fixed_u_trace_ql_only .and. &
            (this%fixed_u_trace_selected_iteration <= 0 .or. &
                             this%fixed_u_trace_iteration == this%fixed_u_trace_selected_iteration)
         if (this%fixed_u_trace_capture_radial_state) then
            if (ncore > 0) EC_SEED(1:ncore) = EC(1:ncore)
            if (nval > 0) EV_SEED(1:nval) = EV(1:nval)
            call this%write_fixed_u_radial_replay_input(atom, rofi, v, rho_in, atom%potential%pl, &
                                                        atom%potential%ql, EC_SEED, ncore, TL, VHRMAX, B_fsm, &
                                                        xc_obj%txc, this%lattice%control%hyperfine)
         end if
         call this%NEWRHO(atom, atom%element%atomic_number, lmax, atom%a, b, nr, rofi, v, rho, atom%potential%PL, atom%potential%QL, SEC, SEV, EC, EV, TL, n_radial_spin_channels, IPR1)
         if (this%fixed_u_trace_capture_radial_state) then
            call this%write_fixed_u_radial_replay_output(rho, SEC, SEV)
            if (.not. this%fixed_u_trace_capture_only) then
               call this%verify_fixed_u_radial_replay(atom, atom%element%atomic_number, lmax, atom%a, b, nr, &
                  rofi_trace_base, v_trace_base, rofi, v, rho_in, rho, VHRMAX, B_fsm, RVH, REPS, RMU, VNUCL, &
                  atom%potential%pl, atom%potential%ql, SEC, SEV, EC_SEED, EV_SEED, EC, EV, ncore, nval, TL)
            end if
            this%fixed_u_trace_capture_radial_state = .false.
            if (allocated(v_trace_base)) deallocate(v_trace_base)
            if (allocated(rofi_trace_base)) deallocate(rofi_trace_base)
         end if
         if (LAST) then
            call atom%radial_ground_state%capture_core_density(this%core_rho_weighted(1:nr, :), &
                                                               this%core_pauli_weighted(1:nr, :), &
                                                               this%core_rho_weighted_l(1:nr, :, :), &
                                                               this%core_pauli_weighted_l(1:nr, :, :))
         end if
         DRHO_CONTROL = 0.d0
         DRHO_INTEGRATED = 0.d0
         SUM = 0.d0
         RHO0T = 0.d0
         do ISP = 1, n_radial_spin_channels
            RHO0T = RHO0T + RHO0(ISP)
            do IR = 1, NR
               WGT = 2*(MOD(IR + 1, 2) + 1)/3.d0
               if (IR == 1 .or. IR == NR) then
                  WGT = 1.d0/3.d0
               end if
               DRDI = atom%a*(rofi(IR) + B)
               DRHO_CONTROL = DRHO_CONTROL + WGT*ABS(RHO(IR, ISP) - rho_in(IR, ISP))
               DRHO_INTEGRATED = DRHO_INTEGRATED + WGT*DRDI*ABS(RHO(IR, ISP) - rho_in(IR, ISP))
               rho_in(IR, ISP) = BETA1*RHO(IR, ISP) + (1.d0 - BETA1)*rho_in(IR, ISP)
               SUM = SUM + WGT*DRDI*RHO(IR, ISP)
            end do
         end do
         if (LAST) then
            call atom%radial_ground_state%mark_accepted(ITER, DRHO_CONTROL)
         end if
         if (LAST) exit
         if (IPR >= 3 .or. (IPR >= 2 .and. (DRHO_CONTROL < TOL .or. ITER == 1 .or. ITER == NITER - 1))) then
            write (6, 10004) ITER, SUM, DRHO_CONTROL, VNUCL, RHO0T, VSUM, BETA1
         end if
         if (atomsc_should_stop(ITER, NITER, DRHO_CONTROL, TOL)) then
            LAST = .true.
         end if
      end do
      ! -----------------------------
      atom%potential%RHOEPS = 0.d0
      RHOMU = 0.d0
      atom%potential%SUMEV = 0.d0
      atom%potential%SUMEC = 0.d0
      RHOVH = 0.d0
      do ISP = 1, n_radial_spin_channels
         if (n_radial_spin_channels == 2 .and. IPR >= 1) then
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
      if (this%fixed_u_trace_enabled .and. rank == 0 .and. this%fixed_u_trace_current_site >= 1 .and. &
          this%fixed_u_trace_iteration > 0 .and. &
          (this%fixed_u_trace_selected_iteration <= 0 .or. &
           this%fixed_u_trace_iteration == this%fixed_u_trace_selected_iteration)) then
         open(newunit=energy_trace_unit, file=trim(this%fixed_u_trace_filename), status='old', &
              position='append', action='write', iostat=energy_trace_ios)
         if (energy_trace_ios == 0) then
            call this%write_fixed_u_csv_row(energy_trace_unit, 'energy_component_sumev', &
               this%fixed_u_trace_iteration, 'ATOMSC_ETOT', this%fixed_u_trace_current_site - 1, &
               this%lattice%ib(this%fixed_u_trace_current_site) - 1, -1, -1, -1, -1, -1, -1, &
               0, 0.0_rp, 0.0_rp, atom%potential%SUMEV, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(energy_trace_unit, 'energy_component_sumec', &
               this%fixed_u_trace_iteration, 'ATOMSC_ETOT', this%fixed_u_trace_current_site - 1, &
               this%lattice%ib(this%fixed_u_trace_current_site) - 1, -1, -1, -1, -1, -1, -1, &
               1, 0.0_rp, 0.0_rp, atom%potential%SUMEC, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(energy_trace_unit, 'energy_component_rhovh', &
               this%fixed_u_trace_iteration, 'ATOMSC_ETOT', this%fixed_u_trace_current_site - 1, &
               this%lattice%ib(this%fixed_u_trace_current_site) - 1, -1, -1, -1, -1, -1, -1, &
               2, 0.0_rp, 0.0_rp, RHOVH, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(energy_trace_unit, 'energy_component_rhomu', &
               this%fixed_u_trace_iteration, 'ATOMSC_ETOT', this%fixed_u_trace_current_site - 1, &
               this%lattice%ib(this%fixed_u_trace_current_site) - 1, -1, -1, -1, -1, -1, -1, &
               3, 0.0_rp, 0.0_rp, RHOMU, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(energy_trace_unit, 'energy_component_vnucl', &
               this%fixed_u_trace_iteration, 'ATOMSC_ETOT', this%fixed_u_trace_current_site - 1, &
               this%lattice%ib(this%fixed_u_trace_current_site) - 1, -1, -1, -1, -1, -1, -1, &
               4, 0.0_rp, 0.0_rp, VNUCL, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(energy_trace_unit, 'energy_component_zvnucl', &
               this%fixed_u_trace_iteration, 'ATOMSC_ETOT', this%fixed_u_trace_current_site - 1, &
               this%lattice%ib(this%fixed_u_trace_current_site) - 1, -1, -1, -1, -1, -1, -1, &
               5, 0.0_rp, 0.0_rp, ZVNUCL, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(energy_trace_unit, 'energy_component_ekin', &
               this%fixed_u_trace_iteration, 'ATOMSC_ETOT', this%fixed_u_trace_current_site - 1, &
               this%lattice%ib(this%fixed_u_trace_current_site) - 1, -1, -1, -1, -1, -1, -1, &
               6, 0.0_rp, 0.0_rp, atom%potential%EKIN, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(energy_trace_unit, 'energy_component_utot', &
               this%fixed_u_trace_iteration, 'ATOMSC_ETOT', this%fixed_u_trace_current_site - 1, &
               this%lattice%ib(this%fixed_u_trace_current_site) - 1, -1, -1, -1, -1, -1, -1, &
               7, 0.0_rp, 0.0_rp, atom%potential%UTOT, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(energy_trace_unit, 'energy_component_rhoeps', &
               this%fixed_u_trace_iteration, 'ATOMSC_ETOT', this%fixed_u_trace_current_site - 1, &
               this%lattice%ib(this%fixed_u_trace_current_site) - 1, -1, -1, -1, -1, -1, -1, &
               8, 0.0_rp, 0.0_rp, atom%potential%RHOEPS, 0.0_rp, 0.0_rp)
            call this%write_fixed_u_csv_row(energy_trace_unit, 'energy_component_etot', &
               this%fixed_u_trace_iteration, 'ATOMSC_ETOT', this%fixed_u_trace_current_site - 1, &
               this%lattice%ib(this%fixed_u_trace_current_site) - 1, -1, -1, -1, -1, -1, -1, &
               9, 0.0_rp, 0.0_rp, atom%potential%ETOT, 0.0_rp, 0.0_rp)
            close(energy_trace_unit)
         end if
      end if
      if (this%fixed_u_trace_current_site >= 1) then
         call this%write_fixed_u_radial_trace(atom, this%fixed_u_trace_current_site, v, rofi, rho_in, rho)
      end if

      if (IPR >= 1) then
         write (6, 10007) &
            atom%potential%SUMEV, atom%potential%SUMEC, VNUCL, RHOVH, ZVNUCL, atom%potential%UTOT, &
            RHOMU, atom%potential%RHOEPS, atom%potential%EKIN, atom%potential%ETOT
      end if
      VRMAX(1) = -2.d0*atom%element%atomic_number/atom%potential%ws_r
      do ISP = 1, n_radial_spin_channels
         VRMAX(1) = VRMAX(1) + V(NR, ISP)/n_radial_spin_channels
      end do
	      VRMAX(2) = 0.d0
	      if (n_radial_spin_channels == 2) then
	         VRMAX(2) = V(NR, 1) - V(NR, 2)
	      end if
	      atom%potential%vrmax(:) = VRMAX(:)
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
         /"  ITER     QINT", 9x, "DRHO_CONTROL", 10x, "VH0", 10x, "RHO0", 10x, "VSUM", 5x, "BETA" &
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

   subroutine NEWRHO(this, atom, Z, LMAX, A, B, NR, rofi, V, RHO, PL, QL, SUMEC, SUMEV, EC, EV, TOL, NSP, IPR, &
                    RHO_D_SHELL)
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
      ! NSP is the number of local radial density channels, not control%nsp’s
      ! global electronic-structure mode.
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
      real(rp), dimension(:, :), intent(out), optional :: RHO_D_SHELL
      !
      !.. Local Scalars ..
      logical :: FREE
      integer :: IFCORE, IIR, IR, ISI, ISP, IVAL, JR, KONFIG, L, LP1, NN, NRE
      real(rp) :: C, DL, DPHI, DPHIP, E, EB1, EB2, EVAL, FLLP1, GFAC, P, &
                  PHI, PHIP, PI, Q0, Q1, Q2, R, RMAX, RO, ROCRIT, SLO, &
                  SUM, TMC, TOLVAL, VAL, RHO_BEFORE_D
      ! HF additions
      real(rp) :: SHUP, SHDW, SUT, ALFAM1, AL2, HYP, OBPI, RT, WGTA, DRDI
      real(rp) :: wave_value
      real(rp), dimension(2) :: SH
      real(rp), dimension(NR) :: HG, DETH
      integer :: IJ
      integer :: trace_unit, trace_ios, radial_index, component
      logical :: trace_radial
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
      if (present(RHO_D_SHELL)) then
         if (size(RHO_D_SHELL, 1) /= NR .or. size(RHO_D_SHELL, 2) /= NSP) &
            error stop 'NEWRHO: d-shell observation shape mismatch'
         RHO_D_SHELL = 0.0_rp
      end if
      trace_radial = this%fixed_u_trace_capture_radial_state .and. this%fixed_u_trace_enabled .and. rank == 0
      trace_unit = -1
      if (trace_radial) then
         open(newunit=trace_unit, file=trim(this%fixed_u_trace_filename), status='old', position='append', &
              action='write', iostat=trace_ios)
         if (trace_ios /= 0) trace_radial = .false.
      end if
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
      ! Only apply f-core treatment when lmax < 3 (f is not in valence basis)
      if (IFCORE /= 0 .and. LMAX < 3) then
         KONF(LMAX + 2) = 5
      end if
      !========FIM  DAS MODIFICACOES NESTE TRECHO================
      do ISP = 1, NSP
         do IR = 1, NR
            RHO(IR, ISP) = 0.d0
         end do
      end do
      call this%RHOCOR(atom, Z, LMAX, KONF, A, B, NR, rofi, V, RHO, G, SUMEC, EC, TOL, NSP, IPR)
      if (trace_radial) then
         do ISP = 1, NSP
            do IR = 1, NR
               call this%write_fixed_u_csv_row(trace_unit, 'radial_replay_rho_core', &
                  this%fixed_u_trace_iteration, 'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, &
                  -1, -1, ISP - 1, -1, -1, IR - 1, ISP - 1, 0, 0.0_rp, rofi(IR), RHO(IR, ISP), 0.0_rp, 0.0_rp)
            end do
         end do
      end if
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
               if (PL(LP1, ISP) /= PL(LP1, ISP)) then
                  call g_logger%fatal('Invalid PL before NEWRHO RSEQSR: atom='//trim(atom%element%symbol)// &
                                      ' lmax='//int2str(LMAX)//' l='//int2str(L)//' spin='//int2str(ISP)// &
                                      ' PL=NaN Q0='//real2str(Q0), __FILE__, __LINE__)
               end if
               KONFIG = PL(LP1, ISP)
               if (KONFIG < LP1) then
                  call g_logger%fatal('Invalid PL before NEWRHO RSEQSR: atom='//trim(atom%element%symbol)// &
                                      ' lmax='//int2str(LMAX)//' l='//int2str(L)//' spin='//int2str(ISP)// &
                                      ' PL='//real2str(PL(LP1, ISP))//' KONFIG='//int2str(KONFIG)// &
                                      ' NN='//int2str(KONFIG - LP1)//' Q0='//real2str(Q0), __FILE__, __LINE__)
               end if
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
               if (trace_radial) then
                  call this%write_fixed_u_csv_row(trace_unit, 'radial_replay_eigen_seed', &
                     this%fixed_u_trace_iteration, 'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, &
                     -1, L, ISP - 1, -1, IVAL - 1, -1, -1, 0, 0.0_rp, 0.0_rp, EVAL, 0.0_rp, 0.0_rp)
               end if
               ! Debug output for orbital radial equation
               ! if (L >= 2) then
               !    write (6, ’(A, I1, A, I1, A, F8.4, A, I2, A, F10.6, A, F10.6)’) &
               !       ″NEWRHO: L=″, L, ″ ISP=″, ISP, ″ PL=″, PL(LP1, ISP), &
               !       ″ NN=″, NN, ″ EVAL=″, EVAL, ″ Q0=″, Q0
               ! end if
               call RSEQSR(EB1, EB2, EVAL, TOLVAL, Z, L, NN, VAL, SLO, V(1, ISP), G, SUM, A, B, &
                           rofi, NR, NRE, 0)
               EV(IVAL) = EVAL
               SUMEV(ISP) = SUMEV(ISP) + EVAL*Q0 + Q1
               if (trace_radial) then
                  call this%write_fixed_u_csv_row(trace_unit, 'radial_replay_nre', &
                     this%fixed_u_trace_iteration, 'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, &
                     -1, L, ISP - 1, -1, IVAL - 1, -1, -1, 0, 0.0_rp, 0.0_rp, real(NRE, rp), 0.0_rp, 0.0_rp)
               end if
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
               if (trace_radial) then
                  call this%write_fixed_u_csv_row(trace_unit, 'radial_replay_eigenvalue', &
                     this%fixed_u_trace_iteration, 'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, &
                     -1, L, ISP - 1, -1, IVAL - 1, -1, -1, 0, 0.0_rp, 0.0_rp, EVAL, 0.0_rp, 0.0_rp)
                  do radial_index = 1, NR
                     do component = 0, 5
                        select case (component)
                        case (0)
                           wave_value = G(radial_index)
                        case (1)
                           wave_value = GP(radial_index)
                        case (2)
                           wave_value = GPP(radial_index)
                        case (3)
                           wave_value = G(radial_index + NR)
                        case (4)
                           wave_value = GP(radial_index + NR)
                        case (5)
                           wave_value = GPP(radial_index + NR)
                        end select
                        call this%write_fixed_u_csv_row(trace_unit, 'radial_replay_wave', &
                           this%fixed_u_trace_iteration, 'RADIAL_REPLAY', this%fixed_u_trace_current_site - 1, &
                           -1, L, ISP - 1, -1, IVAL - 1, radial_index - 1, component, 0, &
                           0.0_rp, rofi(radial_index), wave_value, 0.0_rp, 0.0_rp)
                     end do
                  end do
               end if
               FLLP1 = L*(L + 1)
               do IR = 2, NRE
                  JR = IR + NR
                  R = rofi(IR)
                  TMC = C - (V(IR, ISP) - 2.d0*Z/R - EVAL)/C
                  GFAC = 1.d0 + FLLP1/(TMC*R)**2
                  if (present(RHO_D_SHELL) .and. L == 2) RHO_BEFORE_D = RHO(IR, ISP)
                  RHO(IR, ISP) = &
                     RHO(IR, ISP) + Q0*(GFAC*G(IR)**2 + G(JR)**2) + &
                     2.d0*Q1*(GFAC*G(IR)*GP(IR) + G(JR)*GP(JR)) + &
                     Q2*(GFAC*(GP(IR)**2 + G(IR)*GPP(IR)) + GP(JR)**2 + G(JR)*GPP(JR))
                  if (present(RHO_D_SHELL) .and. L == 2) &
                     RHO_D_SHELL(IR, ISP) = RHO(IR, ISP) - RHO_BEFORE_D
                  !-------- FUN2 GIVES the probability density and NOT the charge density---
                  !----of valence orbital LP1(s, p or d) with spin ISP at point r=IR----
                  this%FUN2(IR, LP1, ISP) = (GFAC*G(IR)**2 + G(JR)**2)
                  !-------- Signed large-component amplitude u_l(r) = r*phi_l(r) for the
                  !----B6 (l=1 dipole) radial matrix elements. Stored parallel to FUN2 so
                  !----phi_s, phi_p, phi_d coexist across LP1 for cross-l products. The
                  !----small component G(JR) is a relativistic correction not needed here.
                  this%phi_amp(IR, LP1, ISP) = G(IR)
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
                     !write(8, ´(a7, f10.6)´) ″Hval:″, HYP
                     atom%potential%hyper_field(2) = hyp
                  end if
               end if
            end if
            !--------------END-ORBITAL CONT. TO HYPERFINE FIELD------------
         end do
      end do
      if (trace_radial) close(trace_unit)
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
      ! NSP is the number of local radial density channels, not control%nsp’s
      ! global electronic-structure mode.
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
      ! Keep the core contribution separate for the LR-02N comparison.  The
      ! production scalar-relativistic density below still receives both
      ! components exactly as before.
      if (allocated(this%core_rho_weighted) .and. allocated(this%core_pauli_weighted)) then
         this%core_rho_weighted(1:NR, :) = 0.0_rp
         this%core_pauli_weighted(1:NR, :) = 0.0_rp
      end if
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
                  if (allocated(this%core_rho_weighted) .and. allocated(this%core_pauli_weighted)) then
                     this%core_rho_weighted(IR, ISP) = this%core_rho_weighted(IR, ISP) + &
                        DEG*(GFAC*G(IR, 1)**2 + G(IR, 2)**2)
                     this%core_pauli_weighted(IR, ISP) = this%core_pauli_weighted(IR, ISP) + &
                        DEG*G(IR, 1)**2
                     if (LP1 <= size(this%core_rho_weighted_l, 2)) then
                        this%core_rho_weighted_l(IR, LP1, ISP) = this%core_rho_weighted_l(IR, LP1, ISP) + &
                           DEG*(GFAC*G(IR, 1)**2 + G(IR, 2)**2)
                        this%core_pauli_weighted_l(IR, LP1, ISP) = this%core_pauli_weighted_l(IR, LP1, ISP) + &
                           DEG*G(IR, 1)**2
                     end if
                  end if
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
               !WRITE(8,´(a7,3f12.6)´)″Hcore: ″,HCORE,(SH(1,1)-SH(1,2)) &
               !,(SH(1,1)-SH(1,2))*914.7744
            END IF
         end if

         !====MODIFICADO PARA INCLUIR F FRACIONARIO NO CAROCO==========
         IFCORE = atom%element%f_core
         DFCORE = REAL(IFCORE)
         ! Only apply f-core treatment when lmax < 3 (f is not in valence basis)
         if (IFCORE /= 0 .and. LMAX < 3) then
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
                  if (allocated(this%core_rho_weighted) .and. allocated(this%core_pauli_weighted)) then
                     this%core_rho_weighted(IR, ISP) = this%core_rho_weighted(IR, ISP) + &
                        DEG*(GFAC*G(IR, 1)**2 + G(IR, 2)**2)
                     this%core_pauli_weighted(IR, ISP) = this%core_pauli_weighted(IR, ISP) + &
                        DEG*G(IR, 1)**2
                     if (LP1 <= size(this%core_rho_weighted_l, 2)) then
                        this%core_rho_weighted_l(IR, LP1, ISP) = this%core_rho_weighted_l(IR, LP1, ISP) + &
                           DEG*(GFAC*G(IR, 1)**2 + G(IR, 2)**2)
                        this%core_pauli_weighted_l(IR, LP1, ISP) = this%core_pauli_weighted_l(IR, LP1, ISP) + &
                           DEG*G(IR, 1)**2
                     end if
                  end if
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

   subroutine legacy_gintsr(G1, G2, A, B, NR, Z, E, L, V, rofi, SUM)
      !> Test-facing wrapper; the production implementation remains GINTSR.
      integer, intent(in) :: L, NR
      real(rp), intent(in) :: A, B, E, Z
      real(rp), intent(out) :: SUM
      real(rp), dimension(NR), intent(in) :: rofi, V
      real(rp), dimension(NR, 2), intent(in) :: G1, G2

      call GINTSR(G1, G2, A, B, NR, Z, E, L, V, rofi, SUM)
   end subroutine legacy_gintsr

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
      ! NSP is the local radial/XC channel count supplied by ATOMSC.
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

   subroutine VXC0SP(this, xc_obj, Z, A, B, rofi, RHO, NR, V, RHO0, RHOEPS, RHOMU, n_radial_spin_channels, B_fsm, snapshot)
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
      ! This argument is the number of local radial/XC spin channels, not
      ! control%nsp’s global collinear/noncollinear/SOC mode.
      integer, intent(in) :: NR, n_radial_spin_channels
      real(rp), intent(in) :: A, B
      real(rp) :: Z
      real(rp), dimension(2), intent(inout) :: RHO0
      real(rp), dimension(NR), intent(in) :: rofi
      real(rp), dimension(n_radial_spin_channels), intent(inout) :: RHOEPS, RHOMU
      real(rp), dimension(NR, n_radial_spin_channels), intent(in) :: RHO
      real(rp), dimension(NR, n_radial_spin_channels), intent(inout) :: V
      real(rp), intent(in)  :: B_fsm
      type(radial_ground_state), intent(inout), optional :: snapshot
      !
      !.. Local Scalars ..
      integer :: IR, ISP, IXC
      real(rp) :: DRDI, EXC, EXC1, EXC2, OB4PI, PI, RHO2, RHO3, RHOT1, &
                  RHOT2, RHOTRU, VXC, VXC1, VXC2, WGT, &
                  RHO1, R, RCE
      real(rp), dimension(NR, n_radial_spin_channels) :: RHOP, RHOPP, tRHO
      real(rp), dimension(2) :: RHOD, RHODD
      real(rp), allocatable :: vxc_up_radial(:), vxc_down_radial(:), exc_radial(:)
      logical :: use_libxc_gga, use_legacy_gga
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
      ! libxc_family is the aggregate route family: any active GGA component
      ! requires the complete radial evaluator, even when an LDA component is
      ! present in the same explicit combination.
      use_libxc_gga = xc_obj%use_libxc .and. xc_obj%libxc_has_gga
      ! Only these historical selectors own the legacy radial GGA machinery.
      ! In particular, a direct/predefined libXC LDA selector must not enter a
      ! legacy IXC>=8 branch merely because its TXC integer is large.
      use_legacy_gga = IXC == 5 .or. IXC == 8 .or. IXC == 9
      !
      ! Constraining field related hacks below
      !
      ! Extrapolate density to core point
      do ISP = 1, n_radial_spin_channels
         RHOEPS(ISP) = 0.d0
         RHOMU(ISP) = 0.d0
         RHO2 = RHO(2, ISP)/rofi(2)**2
         RHO3 = RHO(3, ISP)/rofi(3)**2
         RHO0(ISP) = OB4PI*(RHO2*rofi(3) - RHO3*rofi(2))/(rofi(3) - rofi(2))
      end do
      ! Calculate gradients in case of GGA
      do ISP = 1, n_radial_spin_channels
         tRHO(1, ISP) = RHO0(ISP)
         do IR = 2, NR
            tRHO(IR, 1) = RHO(IR, 1)*OB4PI/rofi(IR)**2
            tRHO(IR, ISP) = RHO(IR, ISP)*OB4PI/rofi(IR)**2
         end do
         if (use_libxc_gga .or. use_legacy_gga) then
            ! radgra returns d(tRHO)/dr; a second call returns d2(tRHO)/dr2.
            ! These are physical radial derivatives, not derivatives with
            ! respect to the logarithmic mesh coordinate.
            call radgra(a, b, NR, rofi, tRHO(1, ISP), RHOP(1, ISP))
            ! Enforce the regular spherical boundary condition before taking
            ! the second derivative, so RHOPP(1) is based on zero slope too.
            RHOP(1, ISP) = 0.d0
            call radgra(a, b, NR, rofi, RHOP(1, ISP), RHOPP(1, ISP))
         else
            RHOP = 0.0d0
            RHOPP = 0.0d0
         end if
      end do

      ! Keep the exact multiplicative XC arrays available until the caller
      ! decides whether this is the accepted final radial iteration.  They
      ! are otherwise automatic temporaries, which was the LR-01 lifetime
      ! blocker.
      allocate (vxc_up_radial(NR), vxc_down_radial(NR))
      if (use_libxc_gga) then
         allocate (exc_radial(NR))
         if (n_radial_spin_channels == 1) then
            call xc_obj%xcpot_libxc_gga_radial(A, B, rofi, 0.5d0*tRHO(:, 1), 0.5d0*tRHO(:, 1), &
                                               0.5d0*RHOP(:, 1), 0.5d0*RHOP(:, 1), &
                                               vxc_up_radial, vxc_down_radial, exc_radial)
         else
            ! The radial helper uses explicit libXC order (up, down); the
            ! historical VXC0SP outputs remain V1=up and V2=down.
            call xc_obj%xcpot_libxc_gga_radial(A, B, rofi, tRHO(:, 1), tRHO(:, 2), &
                                               RHOP(:, 1), RHOP(:, 2), &
                                               vxc_up_radial, vxc_down_radial, exc_radial)
         end if
      end if
      if (n_radial_spin_channels == 1) then
         RHO1 = 0.5d0*tRHO(1, 1)
         RHO2 = RHO1
         if (use_libxc_gga) then
            VXC1 = vxc_up_radial(1)
            VXC2 = vxc_down_radial(1)
            EXC = exc_radial(1)
         else if (use_legacy_gga) then
            ! The VXC0SP RHO1/RHO2 contract is down/up.  radgra already
            ! returned d n/dr and d2 n/dr2; do not transform them again.
            ! At r=0 use zero slope and the regular limit laplacian=3*d2n/dr2.
            call vxc0sp_pbe_origin_derivatives(n_radial_spin_channels, [RHOPP(1, 1), 0.0_rp], RHOD, RHODD)
            if (IXC == 8) then
               R = 0.d0
            else
               ! LAG enhancement is defined by its s->0 limit; use the
               ! first physical radius only as its non-singular coordinate.
               R = rofi(2)
            end if
            call xc_obj%XCPOT_hybrid(RHO2, RHO1, tRHO(1, 1), RHOD, RHODD, R, VXC2, VXC1, EXC)
         else
            ! Preserve the established pointwise LDA path byte-for-byte in
            ! its inputs; the radial derivative machinery is GGA-only.
            R = (rofi(3) - rofi(2))
            call xc_obj%XCPOT_hybrid(RHO2, RHO1, tRHO(1, 1), RHOD, RHODD, R, VXC2, VXC1, EXC)
         end if
         V(1, 1) = V(1, 1) + VXC1
         vxc_up_radial(1) = VXC1
         vxc_down_radial(1) = VXC2
         do IR = 2, NR
            if (use_libxc_gga) then
               VXC1 = vxc_up_radial(IR)
               VXC2 = vxc_down_radial(IR)
               EXC1 = exc_radial(IR)
            else if (use_legacy_gga) then
               RHO1 = 0.5d0*tRHO(IR, 1)
               RHO2 = RHO1
               R = rofi(IR)
               RHOD(1) = 0.5d0*RHOP(IR, 1)
               RHOD(2) = RHOD(1)
               RHODD(1) = 0.5d0*RHOPP(IR, 1)
               RHODD(2) = RHODD(1)
               call xc_obj%XCPOT_hybrid(RHO2, RHO1, tRHO(IR, 1), RHOD, RHODD, R, VXC2, VXC1, EXC1)
            else
               RHO1 = 0.5d0*tRHO(1, 1)
               RHO2 = RHO1
               R = (rofi(3) - rofi(2))
               call xc_obj%XCPOT_hybrid(RHO2, RHO1, RHO(1, 1), RHOD, RHODD, R, VXC2, VXC1, EXC1)
            end if
            V(IR, 1) = V(IR, 1) + VXC1
            vxc_up_radial(IR) = VXC1
            vxc_down_radial(IR) = VXC2
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
         if (use_libxc_gga) then
            VXC1 = vxc_up_radial(1)
            VXC2 = vxc_down_radial(1)
            EXC1 = exc_radial(1)
         else if (use_legacy_gga) then
            ! The radial derivative producer supplies raw dn/dr and d2n/dr2.
            ! PBEGGA applies laplacian(n)=3*d2n/dr2 at the origin.
            call vxc0sp_pbe_origin_derivatives(n_radial_spin_channels, [RHOPP(1, 1), RHOPP(1, 2)], RHOD, RHODD)
            if (IXC == 8) then
               R = 0.d0
            else
               R = rofi(2)
            end if
            call xc_obj%XCPOT_hybrid(RHO2, RHO1, RHO1 + RHO2, RHOD, RHODD, R, VXC2, VXC1, EXC1)
         else
            ! Preserve the established pointwise LDA path and its origin
            ! endpoint; no radial GGA operation is entered for LDA.
            R = (rofi(3) - rofi(2))
            call xc_obj%XCPOT_hybrid(RHO2, RHO1, RHO(1, 1), RHOD, RHODD, R, VXC2, VXC1, EXC1)
         end if
         !V(1, 1) = V(1, 1) + VXC1
         !V(1, 2) = V(1, 2) + VXC2
         V(1, 1) = V(1, 1) + VXC1 + B_fsm
         V(1, 2) = V(1, 2) + VXC2 - B_fsm
         vxc_up_radial(1) = VXC1
         vxc_down_radial(1) = VXC2
         do IR = 2, NR
            R = rofi(IR)
            RCE = R*R
            RHO3 = tRHO(IR, 1) + tRHO(IR, 2)
            RHO1 = tRHO(IR, 1)
            RHO2 = tRHO(IR, 2)
            if (use_libxc_gga) then
               VXC1 = vxc_up_radial(IR)
               VXC2 = vxc_down_radial(IR)
               EXC1 = exc_radial(IR)
            else if (use_legacy_gga) then
               ! RHOP/RHOPP are already dn/dr and d2n/dr2.  Passing a
               ! further 1/r or 1/r2 transform would differentiate in the
               ! logarithmic coordinate a second time.
               RHOD(2) = RHOP(IR, 1)
               RHOD(1) = RHOP(IR, 2)
               RHODD(2) = RHOPP(IR, 1)
               RHODD(1) = RHOPP(IR, 2)
               call xc_obj%XCPOT_hybrid(RHO2, RHO1, RHO3, RHOD, RHODD, R, VXC2, VXC1, EXC1)
            else
               call xc_obj%XCPOT_hybrid(RHO2, RHO1, RHO3, RHOD, RHODD, R, VXC2, VXC1, EXC1)
            end if
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
            vxc_up_radial(IR) = VXC1
            vxc_down_radial(IR) = VXC2
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
         end do
      end if
      if (present(snapshot)) then
         if (n_radial_spin_channels /= 2) then
            error stop 'VXC0SP snapshot requires two local radial spin channels'
         end if
         call snapshot%capture(A, B, rofi, RHO, RHO0, vxc_up_radial, vxc_down_radial, V, B_fsm)
         snapshot%xc_provenance%backend_name = xc_obj%backend_name
         snapshot%xc_provenance%txch = xc_obj%TXCH
         snapshot%xc_provenance%functional_name = xc_obj%functional_name
         snapshot%xc_provenance%mapping_quality = xc_obj%mapping_quality
         snapshot%xc_provenance%txc = xc_obj%txc
         snapshot%xc_provenance%libxc_family = xc_obj%libxc_family
         snapshot%xc_provenance%libxc_nspin = xc_obj%libxc_nspin
         snapshot%xc_provenance%use_libxc = xc_obj%use_libxc
         snapshot%xc_provenance%spin_polarized = (n_radial_spin_channels == 2)
         snapshot%xc_provenance%radial_gga = use_libxc_gga .or. use_legacy_gga
         if (allocated(xc_obj%libxc_func_id)) then
            snapshot%xc_provenance%n_components = size(xc_obj%libxc_func_id)
            allocate(snapshot%xc_provenance%component_ids(size(xc_obj%libxc_func_id)))
            snapshot%xc_provenance%component_ids = xc_obj%libxc_func_id
            if (allocated(xc_obj%libxc_component_family)) then
               allocate(snapshot%xc_provenance%component_families(size(xc_obj%libxc_component_family)))
               snapshot%xc_provenance%component_families = xc_obj%libxc_component_family
            end if
            if (allocated(xc_obj%libxc_component_kind)) then
               allocate(snapshot%xc_provenance%component_kinds(size(xc_obj%libxc_component_kind)))
               snapshot%xc_provenance%component_kinds = xc_obj%libxc_component_kind
            end if
         end if
      end if
   end subroutine VXC0SP

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
      integer :: IDW, II, INUM, INUM1, IR, IR1, ISI, IUP, JM, JP, LMAX, NSP
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
      ! open (15, FILE = ″dracsi″)
      NR = size(ROFI)
      B = atom%B()
      LMAX = atom%potential%lmax
      ! RACSI consumes two local spin eigenchannels even when control%nsp=1.
      NSP = RADIAL_XC_SPIN_CHANNELS
      QSL(:) = 0.0_rp
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
      do INUM = 2, min(3, LMAX + 1)
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
      if (LMAX >= 2) then
         do ISI = 1, NSP
            FAK2 = 0.0d0
            FAK4 = 0.0d0
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
               else
                  FAK4 = SUM/441.d0
               end if
            end do
            if (ISI == 1) then
               QSL(3) = 2.d0*(FAK2 - 5*FAK4)
            else
               QSL(6) = 2.d0*(FAK2 - 5*FAK4)
            end if
         end do
      end if
      if (.not. this%fix_soc) then
         atom%potential%xi_p(:) = [qsl(1), qsl(4)]
         if (LMAX >= 2) then
            atom%potential%xi_d(:) = [qsl(2), qsl(5)]
            atom%potential%rac(:) = [qsl(3), qsl(6)]
         else
            atom%potential%xi_d(:) = 0.0_rp
            atom%potential%rac(:) = 0.0_rp
         end if
      end if
      ! WRITE(17, 56)FAK2, FAK4
      ! WRITE(17, *)RCH
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
      integer :: I, L, LP1, NN, NRE, LMAX, IR
      real(rp) :: B, DLPHI, DLPHIP, DPHI, DPHIP, E, EB1, EB2, OMEGAM, OMEGAP, P, PHI, PHIP, PHMINS, PHPLUS, PI, Q, SLO, SUMM, TOL, TMC, VAL, RMAX
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
      real(rp), dimension(:), allocatable :: sr_correction, sr_gfac, sr_tmc
      integer :: NSP, NR
      ! ... Executable Statements ...
      !
      TOL = 1.d-12
      EB1 = -10.d0
      EB2 = 10.d0
      PI = 4.d0*ATAN(1.d0)
      B = atom%B()
      ! Potential parameters are built for the two local radial spin channels;
      ! this is not a query of the global control%nsp mode.
      nsp = RADIAL_XC_SPIN_CHANNELS
      rmax = atom%potential%ws_r
      nr = size(ROFI)
      lmax = atom%potential%lmax
      allocate (G(NR, NSP), GP(NR, NSP), GPP(NR, NSP))
      allocate (sr_correction(NR), sr_gfac(NR), sr_tmc(NR))
      if (atom%radial_ground_state%valid) call atom%radial_ground_state%begin_pauli_basis(lmax)
      do I = 1, nsp
         do L = 0, lmax
            if (atom%potential%PNU(L, I) /= atom%potential%PNU(L, I)) then
               call g_logger%fatal('Invalid PNU before POTPAR RSEQSR: atom='//trim(atom%element%symbol)// &
                                   ' lmax='//int2str(LMAX)//' l='//int2str(L)//' spin='//int2str(I)// &
                                   ' PNU=NaN', __FILE__, __LINE__)
            end if
            KONFIG(L) = int(atom%potential%PNU(L, I))
            if (KONFIG(L) < L + 1) then
               call g_logger%fatal('Invalid PNU before POTPAR RSEQSR: atom='//trim(atom%element%symbol)// &
                                   ' lmax='//int2str(LMAX)//' l='//int2str(L)//' spin='//int2str(I)// &
                                   ' PNU='//real2str(atom%potential%PNU(L, I))//' KONFIG='//int2str(KONFIG(L))// &
                                   ' NN='//int2str(KONFIG(L) - L - 1), __FILE__, __LINE__)
            end if
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
            sr_correction = 0.0_rp
            sr_gfac = 1.0_rp
            sr_tmc = 274.074_rp
            do IR = 2, NR
               if (abs(ROFI(IR)) > tiny(1.0_rp)) then
                  TMC = 274.074_rp - (V(IR, I) - 2.0_rp*atom%element%atomic_number/ROFI(IR) - E)/274.074_rp
                  sr_tmc(IR) = TMC
                  sr_gfac(IR) = 1.0_rp + real(L*(L + 1), rp)/(TMC*ROFI(IR))**2
                  sr_correction(IR) = real(L*(L + 1), rp)/(TMC*ROFI(IR))**2*G(IR, 1)**2 + G(IR, 2)**2
               end if
            end do
            if (atom%radial_ground_state%valid) then
               call atom%radial_ground_state%set_pauli_basis_channel(L, I, E, G(:, 1), GP(:, 1), G(:, 2), sr_correction)
               call atom%radial_ground_state%set_scalar_relativistic_basis_channel(L, I, E, GP(:, 1), GP(:, 2), &
                  GPP(:, 1), GPP(:, 2), sr_gfac, sr_tmc, V(:, I))
            end if
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
      deallocate (sr_correction, sr_gfac, sr_tmc)
      return
      !
      ! ... Format Declarations ...
      !
      ! 10000 format (/″ POTPAR:  INPUT LOG DERIVATIVES ARE″/″ PNU=″, 6f12.6)
      ! 10001 format (″ DNU=″, 6f12.6)
      ! 10002 format (i3, 4f12.5, 2f12.4)
      ! 10003 format (/″  L″, 7x, ″ENU″, 9x, ″V″, 11x, ″C″, 10x, ″SRDEL″, 8x, ″1/Q″, 9x, ″1/SRP″)
   end subroutine POTPAR

end module self_mod
