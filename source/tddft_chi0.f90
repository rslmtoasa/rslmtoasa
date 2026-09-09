!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Bare finite-temperature Kohn-Sham response from spinor eigenpairs.
!>
!> This is the explicit eigenpair backend for TDDFT.  Frequencies and eta are
!> energies in Rydberg and the denominator is retarded,
!> `omega + e_n - e_m + i*eta`.  The right transition factor is evaluated as
!> `<m|B|n>` rather than obtained by assuming that B is Hermitian.  That
!> distinction is essential for the circular `chi^{+-}` channel, where the
!> response operator and perturbing operator are adjoints of one another.
module tddft_chi0_mod
   use precision_mod, only: rp
   use math_mod, only: pi
   use response_vertices_mod, only: response_channel
   use tddft_occupation_mod, only: tddft_response_occupation_state, validate_response_occupation_fields
   use tddft_transition_engine_mod, only: tddft_transition_engine, site_channel_vertex_provider, &
      make_site_channel_vertex_provider
   implicit none

   private

   ! Keep this exactly aligned with reciprocal_occupations.f90.  TDDFT-06 can
   ! later obtain these values directly from the reciprocal/configuration
   ! objects; the array provider intentionally has no reciprocal-state side
   ! effects.
   real(rp), parameter, public :: tddft_kB_Ry_per_K = 6.3336814e-6_rp
   real(rp), parameter, public :: tddft_occupation_kT_floor = 1.0e-10_rp

   !> Inputs controlling the exact all-band transition path.  A zero
   !> occupation_prune_tolerance disables pruning; a positive value is an
   !> explicitly approximate performance option.  The shared transition engine
   !> has a selectable scalar reduction for deterministic equivalence tests and
   !> a batched BLAS path for production performance.
   type, public :: tddft_chi0_options
      real(rp) :: eta = 0.0_rp
      real(rp) :: fermi_level = huge(1.0_rp)
      real(rp) :: electronic_temperature = 0.0_rp
      integer :: band_first = 1
      integer :: band_last = 0
      real(rp) :: occupation_prune_tolerance = 0.0_rp
      integer :: k_mesh_shape(3) = 0
      logical :: use_batched_accumulation = .true.
      integer :: transition_batch_size = 128
      real(rp) :: q_direct(3) = 0.0_rp
      real(rp) :: q_cartesian(3) = 0.0_rp
      character(len=128) :: fourier_phase_convention = 'unresolved'
      logical :: kq_endpoint_folded = .false.
      integer :: kq_reciprocal_shift(3) = 0
      integer :: kq_folded_endpoint_count = 0
      character(len=32) :: response_projection = 'site'
      character(len=16) :: circular_channel = 'plus_minus'
      type(tddft_response_occupation_state) :: occupation_state
      character(len=48) :: fermi_source = 'unresolved'
      character(len=48) :: fermi_policy = 'unresolved'
   end type tddft_chi0_options

   !> Reproducibility metadata written with every chi_KS output.
   type, public :: tddft_chi0_metadata
      character(len=32) :: backend = 'eigenpairs'
      character(len=32) :: canonical_backend = 'eigenpairs'
      character(len=64) :: implementation = 'explicit eigenpair transitions'
      character(len=96) :: energy_integration = 'not applicable'
      character(len=16) :: energy_unit = 'Rydberg'
      character(len=32) :: susceptibility_unit = '1/Rydberg'
      character(len=96) :: frequency_convention = 'retarded: omega is energy; denominator omega+en-em+i*eta'
      character(len=96) :: spectral_convention = 'Stoner spectral weight = -Im chi_KS^{+-}/pi (positive excitation positive)'
      character(len=64) :: endpoint_provenance = 'separate k and k+q eigenpair arrays'
      character(len=32) :: response_projection = 'site'
      character(len=16) :: circular_channel = 'unspecified'
      character(len=80) :: eta_role = 'numerical broadening; not a physical linewidth'
      real(rp) :: eta = 0.0_rp
      real(rp) :: fermi_level = huge(1.0_rp)
      real(rp) :: electronic_temperature = 0.0_rp
      real(rp) :: electronic_kT = 0.0_rp
      real(rp) :: k_weight_sum = 0.0_rp
      integer :: k_mesh_shape(3) = 0
      integer :: nk = 0
      integer :: available_band_count = 0
      integer :: band_first = 0
      integer :: band_last = 0
      real(rp) :: occupation_prune_tolerance = 0.0_rp
      logical :: batched_accumulation = .false.
      integer :: transition_batch_size = 0
      real(rp) :: vertex_cpu_seconds = 0.0_rp
      real(rp) :: transition_preparation_cpu_seconds = 0.0_rp
      real(rp) :: denominator_cpu_seconds = 0.0_rp
      real(rp) :: accumulation_cpu_seconds = 0.0_rp
      real(rp) :: arbitrary_kq_cpu_seconds = 0.0_rp
      real(rp) :: green_energy_integration_cpu_seconds = 0.0_rp
      integer :: green_function_evaluations = 0
      real(rp) :: green_eta = 0.0_rp
      real(rp) :: integration_energy_min = 0.0_rp
      real(rp) :: integration_energy_max = 0.0_rp
      integer :: integration_energy_points = 0
      integer :: contour_points = 0
      integer :: contour_subdivisions = 0
      integer :: near_fermi_points = 0
      real(rp) :: contour_height = 0.0_rp
      real(rp) :: contour_max_imaginary_energy = 0.0_rp
      integer :: contour_gf_evaluations = 0
      logical :: contour_deformation = .false.
      real(rp) :: q_direct(3) = 0.0_rp
      real(rp) :: q_cartesian(3) = 0.0_rp
      character(len=128) :: fourier_phase_convention = 'unresolved'
      logical :: kq_endpoint_folded = .false.
      integer :: kq_reciprocal_shift(3) = 0
      integer :: kq_folded_endpoint_count = 0
      real(rp) :: omega_min = 0.0_rp
      real(rp) :: omega_max = 0.0_rp
      integer :: omega_points = 0
      logical :: static_limit = .false.
      logical :: eta_is_numerical = .true.
      character(len=48) :: convergence_status = 'not evaluated'
      logical :: converged = .false.
      logical :: goldstone_correction_applied = .false.
      character(len=96) :: goldstone_correction = 'none'
      integer :: q_batch_size = 1
      integer :: omega_batch_size = 0
      logical :: real_space_reuse = .false.
      integer :: real_space_points = 0
      real(rp) :: real_space_cutoff = 0.0_rp
      real(rp) :: real_space_requested_cutoff = 0.0_rp
      real(rp) :: real_space_source_radius = 0.0_rp
      character(len=16) :: real_space_representation = 'bulk'
      character(len=16) :: real_space_truncation_mode = 'full_tail'
      integer :: real_space_fourier_axes(3) = [1, 2, 3]
      integer :: real_space_omitted_points = 0
      integer :: real_space_shell_count = 0
      real(rp) :: real_space_selected_norm = 0.0_rp
      real(rp) :: real_space_tail_norm = 0.0_rp
      real(rp) :: real_space_total_norm = 0.0_rp
      real(rp) :: real_space_tail_ratio = 0.0_rp
      real(rp) :: real_space_last_shell_norm = 0.0_rp
      real(rp) :: real_space_pair_integration_cpu_seconds = 0.0_rp
      real(rp) :: real_space_fourier_cpu_seconds = 0.0_rp
      real(rp) :: real_space_green_cpu_seconds = 0.0_rp
      real(rp) :: real_space_total_cpu_seconds = 0.0_rp
      integer :: real_space_pair_response_integrations = 0
      real(rp) :: real_space_tail_tolerance = 0.0_rp
      logical :: real_space_tail_assessed = .false.
      logical :: real_space_source_covers_cutoff = .false.
      character(len=48) :: fermi_source = 'unresolved'
      character(len=48) :: fermi_policy = 'unresolved'
   end type tddft_chi0_metadata

   !> Response and directly consumable KS/Stoner spectral products.  The
   !> matrix dimensions are (left channel, right channel, frequency).  The
   !> site-diagonal array is (site channel, frequency), and is populated only
   !> for matching left/right site channels.
   type, public :: tddft_chi0_result
      complex(rp), allocatable :: chi(:, :, :)
      real(rp), allocatable :: re_chi(:, :, :)
      real(rp), allocatable :: im_chi(:, :, :)
      real(rp), allocatable :: site_diagonal_spectrum(:, :)
      real(rp), allocatable :: trace_spectrum(:)
      real(rp), allocatable :: stoner_spectral_map(:, :)
      type(tddft_chi0_metadata) :: metadata
   end type tddft_chi0_result

   !> Common request envelope consumed by every chi0 backend.  `q_points`
   !> uses the project-wide (3,nq) convention and may contain one point for a
   !> scalar request.  `q_indices` is optional backend provenance: a caller
   !> that has already built endpoint data can identify those points without
   !> requiring a backend to reverse-map floating-point coordinates.
   type, public :: tddft_chi0_request
      real(rp), allocatable :: q_points(:, :)
      integer, allocatable :: q_indices(:)
      real(rp), allocatable :: omega(:)
      character(len=24) :: batch_mode = 'q_and_omega'
      logical :: allow_real_space_reuse = .true.
   end type tddft_chi0_request

   !> Batch envelope around the compatibility-preserving per-q result.  The
   !> backend contract can therefore batch q points without changing existing
   !> Dyson and output consumers, which continue to receive one
   !> `tddft_chi0_result` at a time when appropriate.
   type, public :: tddft_chi0_batch_result
      type(tddft_chi0_result), allocatable :: q_response(:)
      real(rp), allocatable :: q_points(:, :)
      integer, allocatable :: q_indices(:)
      real(rp), allocatable :: omega(:)
      type(tddft_chi0_metadata) :: metadata
   end type tddft_chi0_batch_result

   public :: build_chi_ks_from_eigenpairs
   public :: build_static_chi_ks_from_eigenpairs
   public :: build_static_chi_ks_from_eigenpairs_at_q
   public :: tddft_fermi_occupation
   public :: tddft_static_divided_difference
   public :: write_chi_ks_text
   public :: install_tddft_occupation
   public :: validate_tddft_chi0_options

contains

   subroutine install_tddft_occupation(options, state)
      type(tddft_chi0_options), intent(inout) :: options
      type(tddft_response_occupation_state), intent(in) :: state

      call state%validate('install_tddft_occupation')
      options%occupation_state = state
      options%fermi_level = state%fermi_level
      options%electronic_temperature = state%electronic_temperature
      options%occupation_prune_tolerance = state%occupation_tolerance
      options%band_first = state%band_first
      options%band_last = state%band_last
      options%fermi_source = state%fermi_source
      options%fermi_policy = state%fermi_policy
   end subroutine install_tddft_occupation

   subroutine validate_tddft_chi0_options(options, context)
      type(tddft_chi0_options), intent(in) :: options
      character(len=*), intent(in) :: context

      call validate_response_occupation_fields(options%fermi_level, options%electronic_temperature, &
         options%occupation_prune_tolerance, options%band_first, options%band_last, context)
      if (options%occupation_state%fermi_is_resolved) then
         call options%occupation_state%validate(trim(context)//': occupation_state')
         if (options%fermi_level /= options%occupation_state%fermi_level .or. &
             options%electronic_temperature /= options%occupation_state%electronic_temperature .or. &
             options%occupation_prune_tolerance /= options%occupation_state%occupation_tolerance .or. &
             options%band_first /= options%occupation_state%band_first .or. options%band_last /= options%occupation_state%band_last .or. &
             trim(options%fermi_source) /= trim(options%occupation_state%fermi_source) .or. &
             trim(options%fermi_policy) /= trim(options%occupation_state%fermi_policy)) then
            error stop trim(context)//': response occupation state and backend fields diverge'
         end if
      end if
   end subroutine validate_tddft_chi0_options

   !> Stable Fermi occupation using the existing reciprocal-space temperature
   !> convention (temperature in K, energies in Ry).  It intentionally shares
   !> reciprocal_occupations.f90’s kT floor and exponential cutoffs.
   pure real(rp) function tddft_fermi_occupation(eigenvalue, fermi_level, temperature) result(occupation)
      real(rp), intent(in) :: eigenvalue, fermi_level, temperature
      real(rp) :: argument, kT

      kT = max(temperature*tddft_kB_Ry_per_K, tddft_occupation_kT_floor)
      argument = (eigenvalue - fermi_level)/kT
      if (argument >= 50.0_rp) then
         occupation = 0.0_rp
      else if (argument <= -50.0_rp) then
         occupation = 1.0_rp
      else
         occupation = 1.0_rp/(exp(argument) + 1.0_rp)
      end if
   end function tddft_fermi_occupation

   !> Build the generalized bare KS response for batches of energy-valued
   !> frequencies.  `eigenvalues_k` and `eigenvalues_kq` hold the n and m
   !> states, respectively; callers obtain the latter from TDDFT-01’s exact
   !> arbitrary-k service.  k weights may be normalized probabilities or raw
   !> multiplicities and are normalized by their explicit sum.
   subroutine build_chi_ks_from_eigenpairs(k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_kq, &
      eigenvectors_kq, site_orbital_counts, left_channels, right_channels, omega, options, result)
      real(rp), target, intent(in) :: k_weights(:)
      real(rp), intent(in) :: eigenvalues_k(:, :), eigenvalues_kq(:, :), omega(:)
      complex(rp), target, intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :)
      integer, target, intent(in) :: site_orbital_counts(:)
      type(response_channel), target, intent(in) :: left_channels(:), right_channels(:)
      type(tddft_chi0_options), intent(in) :: options
      type(tddft_chi0_result), intent(out) :: result

      integer :: nk, nbands, nspinor, nleft, nright, nw, batch_size
      integer :: band_first, band_last
      real(rp) :: weight_sum
      type(tddft_transition_engine) :: engine
      type(site_channel_vertex_provider) :: provider

      nk = size(k_weights)
      nbands = size(eigenvalues_k, 1)
      nleft = size(left_channels)
      nright = size(right_channels)
      nw = size(omega)
      nspinor = 2*sum(site_orbital_counts)
      call validate_chi_ks_inputs(nk, nbands, nspinor, nleft, nright, nw, k_weights, eigenvalues_k, &
         eigenvectors_k, eigenvalues_kq, eigenvectors_kq, options, require_eta=.true.)

      band_first = options%band_first
      band_last = options%band_last
      if (band_last == 0) band_last = nbands
      if (band_first < 1 .or. band_last < band_first .or. band_last > nbands) then
         error stop 'build_chi_ks_from_eigenpairs: invalid selected band window'
      end if
      weight_sum = sum(k_weights)

      allocate(result%chi(nleft, nright, nw), result%re_chi(nleft, nright, nw), &
         result%im_chi(nleft, nright, nw))
      result%chi = cmplx(0.0_rp, 0.0_rp, rp)
      batch_size = min(options%transition_batch_size, (band_last-band_first+1)**2)
      call make_site_channel_vertex_provider(provider, site_orbital_counts, left_channels, right_channels, &
         eigenvectors_k, eigenvectors_kq)
      call engine%accumulate_dynamic(k_weights, eigenvalues_k, eigenvalues_kq, omega, options%eta, options%fermi_level, &
         options%electronic_temperature, band_first, band_last, options%occupation_prune_tolerance, batch_size, &
         options%use_batched_accumulation, provider, result%chi, result%metadata%vertex_cpu_seconds, &
         result%metadata%transition_preparation_cpu_seconds, result%metadata%denominator_cpu_seconds, &
         result%metadata%accumulation_cpu_seconds)

      result%re_chi = real(result%chi, rp)
      result%im_chi = aimag(result%chi)
      call build_spectral_products(left_channels, right_channels, result)
      result%metadata%eta = options%eta
      result%metadata%fermi_level = options%fermi_level
      result%metadata%fermi_source = options%fermi_source
      result%metadata%fermi_policy = options%fermi_policy
      result%metadata%electronic_temperature = options%electronic_temperature
      result%metadata%electronic_kT = max(options%electronic_temperature*tddft_kB_Ry_per_K, &
         tddft_occupation_kT_floor)
      result%metadata%k_weight_sum = weight_sum
      result%metadata%k_mesh_shape = options%k_mesh_shape
      result%metadata%nk = nk
      result%metadata%available_band_count = nbands
      result%metadata%band_first = band_first
      result%metadata%band_last = band_last
      result%metadata%occupation_prune_tolerance = options%occupation_prune_tolerance
      result%metadata%batched_accumulation = options%use_batched_accumulation
      if (options%use_batched_accumulation) result%metadata%transition_batch_size = batch_size
      result%metadata%q_direct = options%q_direct
      result%metadata%q_cartesian = options%q_cartesian
      result%metadata%fourier_phase_convention = options%fourier_phase_convention
      result%metadata%kq_endpoint_folded = options%kq_endpoint_folded
      result%metadata%kq_reciprocal_shift = options%kq_reciprocal_shift
      result%metadata%kq_folded_endpoint_count = options%kq_folded_endpoint_count
      result%metadata%omega_min = minval(omega)
      result%metadata%omega_max = maxval(omega)
      result%metadata%omega_points = nw
      result%metadata%static_limit = .false.
      result%metadata%eta_is_numerical = .true.
      result%metadata%response_projection = options%response_projection
      result%metadata%circular_channel = trim(options%circular_channel)
      result%metadata%spectral_convention = 'Stoner spectral weight = -Im chi_KS^'//trim(options%circular_channel)//'/pi (positive excitation positive)'
      result%metadata%canonical_backend = 'eigenpairs'
      result%metadata%implementation = 'explicit eigenpair transitions'
      result%metadata%convergence_status = 'not assessed by backend'
      result%metadata%converged = .false.
      result%metadata%q_batch_size = 1
      result%metadata%omega_batch_size = nw
   end subroutine build_chi_ks_from_eigenpairs

   !> Real q=0, omega=0 Lehmann response used only for static Ward
   !> diagnostics.  It deliberately has no eta argument: the n=m and nearly
   !> degenerate limit is the derivative of the same finite-temperature Fermi
   !> function used by the dynamic response, (f_n-f_m)/(e_n-e_m) -> f’(e).
   subroutine build_static_chi_ks_from_eigenpairs(k_weights, eigenvalues, eigenvectors, site_orbital_counts, &
      left_channels, right_channels, options, result)
      real(rp), target, intent(in) :: k_weights(:), eigenvalues(:, :)
      complex(rp), target, intent(in) :: eigenvectors(:, :, :)
      integer, target, intent(in) :: site_orbital_counts(:)
      type(response_channel), target, intent(in) :: left_channels(:), right_channels(:)
      type(tddft_chi0_options), intent(in) :: options
      type(tddft_chi0_result), intent(out) :: result

      type(tddft_chi0_options) :: static_options

      static_options = options
      static_options%q_direct = 0.0_rp
      call build_static_chi_ks_from_eigenpairs_at_q(k_weights, eigenvalues, eigenvectors, eigenvalues, eigenvectors, &
         site_orbital_counts, left_channels, right_channels, static_options, result)
      result%metadata%endpoint_provenance = 'q=0: shared k endpoint eigenpair arrays'
   end subroutine build_static_chi_ks_from_eigenpairs

   !> Build the true static `chi_KS(q,0)` response from separate endpoint
   !> eigenpairs.  The static factor is the finite-temperature divided
   !> difference `(f_n(k)-f_m(k+q))/(e_n(k)-e_m(k+q))`; it is not obtained by
   !> evaluating the dynamic response at a large broadening.  This path is the
   !> static counterpart of the dynamic endpoint contract and is suitable for
   !> finite-q Ward/convergence diagnostics as well as q=0.
   subroutine build_static_chi_ks_from_eigenpairs_at_q(k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_kq, &
      eigenvectors_kq, site_orbital_counts, left_channels, right_channels, options, result)
      real(rp), target, intent(in) :: k_weights(:), eigenvalues_k(:, :), eigenvalues_kq(:, :)
      complex(rp), target, intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :)
      integer, target, intent(in) :: site_orbital_counts(:)
      type(response_channel), target, intent(in) :: left_channels(:), right_channels(:)
      type(tddft_chi0_options), intent(in) :: options
      type(tddft_chi0_result), intent(out) :: result

      integer :: nk, nbands, nspinor, nleft, nright, band_first, band_last, batch_size
      real(rp) :: weight_sum
      type(tddft_transition_engine) :: engine
      type(site_channel_vertex_provider) :: provider

      nk = size(k_weights)
      nbands = size(eigenvalues_k, 1)
      nleft = size(left_channels)
      nright = size(right_channels)
      nspinor = 2*sum(site_orbital_counts)
      call validate_chi_ks_inputs(nk, nbands, nspinor, nleft, nright, 1, k_weights, eigenvalues_k, eigenvectors_k, &
         eigenvalues_kq, eigenvectors_kq, options, require_eta=.false.)
      if (options%eta < 0.0_rp) error stop 'build_static_chi_ks_from_eigenpairs_at_q: eta must be non-negative'
      band_first = options%band_first
      band_last = options%band_last
      if (band_last == 0) band_last = nbands
      if (band_first < 1 .or. band_last < band_first .or. band_last > nbands) then
         error stop 'build_static_chi_ks_from_eigenpairs_at_q: invalid selected band window'
      end if
      weight_sum = sum(k_weights)
      allocate(result%chi(nleft, nright, 1), result%re_chi(nleft, nright, 1), result%im_chi(nleft, nright, 1))
      result%chi = cmplx(0.0_rp, 0.0_rp, rp)
      batch_size = min(options%transition_batch_size, (band_last-band_first+1)**2)
      call make_site_channel_vertex_provider(provider, site_orbital_counts, left_channels, right_channels, &
         eigenvectors_k, eigenvectors_kq)
      call engine%accumulate_static_shifted(k_weights, eigenvalues_k, eigenvalues_kq, options%fermi_level, &
         options%electronic_temperature, band_first, band_last, options%occupation_prune_tolerance, batch_size, provider, &
         result%chi, result%metadata%vertex_cpu_seconds, result%metadata%transition_preparation_cpu_seconds, &
         result%metadata%accumulation_cpu_seconds)
      result%re_chi = real(result%chi, rp)
      result%im_chi = aimag(result%chi)
      call build_spectral_products(left_channels, right_channels, result)
      result%metadata%backend = 'static_eigenpairs'
      result%metadata%canonical_backend = 'eigenpairs'
      result%metadata%implementation = 'explicit eigenpair static divided difference'
      result%metadata%energy_integration = 'not applicable'
      result%metadata%frequency_convention = 'static omega=0 divided difference; no dynamical eta'
      result%metadata%eta_role = 'not used by static divided-difference response'
      result%metadata%eta = 0.0_rp
      result%metadata%eta_is_numerical = .false.
      result%metadata%static_limit = .true.
      result%metadata%fermi_level = options%fermi_level
      result%metadata%fermi_source = options%fermi_source
      result%metadata%fermi_policy = options%fermi_policy
      result%metadata%electronic_temperature = options%electronic_temperature
      result%metadata%electronic_kT = max(options%electronic_temperature*tddft_kB_Ry_per_K, tddft_occupation_kT_floor)
      result%metadata%k_weight_sum = weight_sum
      result%metadata%k_mesh_shape = options%k_mesh_shape
      result%metadata%nk = nk
      result%metadata%available_band_count = nbands
      result%metadata%band_first = band_first
      result%metadata%band_last = band_last
      result%metadata%occupation_prune_tolerance = options%occupation_prune_tolerance
      result%metadata%batched_accumulation = .true.
      result%metadata%transition_batch_size = batch_size
      result%metadata%q_direct = options%q_direct
      result%metadata%q_cartesian = options%q_cartesian
      result%metadata%fourier_phase_convention = options%fourier_phase_convention
      result%metadata%kq_endpoint_folded = options%kq_endpoint_folded
      result%metadata%kq_reciprocal_shift = options%kq_reciprocal_shift
      result%metadata%kq_folded_endpoint_count = options%kq_folded_endpoint_count
      result%metadata%omega_min = 0.0_rp
      result%metadata%omega_max = 0.0_rp
      result%metadata%omega_points = 1
      result%metadata%response_projection = options%response_projection
      result%metadata%circular_channel = trim(options%circular_channel)
      result%metadata%spectral_convention = 'Stoner spectral weight = -Im chi_KS^'//trim(options%circular_channel)//'/pi (positive excitation positive)'
      result%metadata%convergence_status = 'not assessed by backend'
      result%metadata%converged = .false.
      result%metadata%q_batch_size = 1
      result%metadata%omega_batch_size = 1
   end subroutine build_static_chi_ks_from_eigenpairs_at_q

   pure real(rp) function tddft_static_divided_difference(energy_n, energy_m, fermi_level, temperature) result(value)
      real(rp), intent(in) :: energy_n, energy_m, fermi_level, temperature
      real(rp) :: delta, midpoint, scale, occupation, kT

      delta = energy_n-energy_m
      kT = max(temperature*tddft_kB_Ry_per_K, tddft_occupation_kT_floor)
      scale = max(1.0_rp, abs(energy_n), abs(energy_m), kT)
      if (abs(delta) > 32.0_rp*sqrt(epsilon(1.0_rp))*scale) then
         value = (tddft_fermi_occupation(energy_n, fermi_level, temperature) - &
                  tddft_fermi_occupation(energy_m, fermi_level, temperature))/delta
      else
         midpoint = 0.5_rp*(energy_n + energy_m)
         occupation = tddft_fermi_occupation(midpoint, fermi_level, temperature)
         if (abs((midpoint-fermi_level)/kT) >= 50.0_rp) then
            value = 0.0_rp
         else
            value = -occupation*(1.0_rp-occupation)/kT
         end if
      end if
   end function tddft_static_divided_difference

   !> Write a self-describing plain-text chi_KS result.  The `matrix` records
   !> contain Re and Im chi; `site_diagonal` and `trace` are the explicitly
   !> labelled KS/Stoner products using -Im(chi)/pi.
   subroutine write_chi_ks_text(filename, omega, result)
      character(len=*), intent(in) :: filename
      real(rp), intent(in) :: omega(:)
      type(tddft_chi0_result), intent(in) :: result
      integer :: unit, ios, iw, ileft, iright, idiag

      if (.not. allocated(result%chi) .or. size(omega) /= size(result%chi, 3)) then
         error stop 'write_chi_ks_text: omega/result shape mismatch'
      end if
      open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
      if (ios /= 0) error stop 'write_chi_ks_text: cannot open output file'
      write(unit, '(a)') '# quantity = bare chi_KS; no TDDFT enhancement; eta is numerical for dynamic response only'
      write(unit, '(a,a)') '# energy_unit = ', trim(result%metadata%energy_unit)
      write(unit, '(a,a)') '# susceptibility_unit = ', trim(result%metadata%susceptibility_unit)
      write(unit, '(a,a)') '# frequency_convention = ', trim(result%metadata%frequency_convention)
      write(unit, '(a,a)') '# spectral_convention = ', trim(result%metadata%spectral_convention)
      write(unit, '(a,a)') '# endpoint_provenance = ', trim(result%metadata%endpoint_provenance)
      write(unit, '(a,a)') '# response_projection = ', trim(result%metadata%response_projection)
      write(unit, '(a,a)') '# circular_channel = ', trim(result%metadata%circular_channel)
      write(unit, '(a,a)') '# eta_role = ', trim(result%metadata%eta_role)
      write(unit, '(a)') '# stoner_landau_map = site-diagonal KS loss; q resolution is supplied by the per-q output files'
      write(unit, '(a,3(1x,es24.16))') '# q_direct = ', result%metadata%q_direct
      write(unit, '(a,3(1x,es24.16))') '# q_cartesian = ', result%metadata%q_cartesian
      write(unit, '(a,a)') '# fourier_phase_convention = ', trim(result%metadata%fourier_phase_convention)
      write(unit, '(a,l1)') '# kq_endpoint_folded = ', result%metadata%kq_endpoint_folded
      write(unit, '(a,3(1x,i0))') '# kq_reciprocal_shift_first = ', result%metadata%kq_reciprocal_shift
      write(unit, '(a,i0)') '# kq_folded_endpoint_count = ', result%metadata%kq_folded_endpoint_count
      write(unit, '(a,l1)') '# static_limit = ', result%metadata%static_limit
      write(unit, '(a,l1)') '# eta_is_numerical = ', result%metadata%eta_is_numerical
      write(unit, '(a,a)') '# chi0_backend = ', trim(result%metadata%backend)
      write(unit, '(a,a)') '# chi0_backend_canonical = ', trim(result%metadata%canonical_backend)
      write(unit, '(a,a)') '# chi0_backend_implementation = ', trim(result%metadata%implementation)
      write(unit, '(a,a)') '# convergence_status = ', trim(result%metadata%convergence_status)
      write(unit, '(a,l1)') '# converged = ', result%metadata%converged
      write(unit, '(a,i0)') '# q_batch_size = ', result%metadata%q_batch_size
      write(unit, '(a,i0)') '# omega_batch_size = ', result%metadata%omega_batch_size
      write(unit, '(a,l1)') '# goldstone_correction_applied = ', result%metadata%goldstone_correction_applied
      write(unit, '(a,a)') '# goldstone_correction = ', trim(result%metadata%goldstone_correction)
      write(unit, '(a,l1)') '# real_space_reuse = ', result%metadata%real_space_reuse
      write(unit, '(a,i0)') '# real_space_points = ', result%metadata%real_space_points
      write(unit, '(a,es24.16)') '# real_space_cutoff = ', result%metadata%real_space_cutoff
      write(unit, '(a,es24.16)') '# real_space_requested_cutoff = ', result%metadata%real_space_requested_cutoff
      write(unit, '(a,es24.16)') '# real_space_source_radius = ', result%metadata%real_space_source_radius
      write(unit, '(a,a)') '# real_space_representation = ', trim(result%metadata%real_space_representation)
      write(unit, '(a,a)') '# real_space_truncation_mode = ', trim(result%metadata%real_space_truncation_mode)
      write(unit, '(a,3(1x,i0))') '# real_space_fourier_axes =', result%metadata%real_space_fourier_axes
      write(unit, '(a,i0)') '# real_space_omitted_points = ', result%metadata%real_space_omitted_points
      write(unit, '(a,i0)') '# real_space_shell_count = ', result%metadata%real_space_shell_count
      write(unit, '(a,es24.16)') '# real_space_selected_norm = ', result%metadata%real_space_selected_norm
      write(unit, '(a,es24.16)') '# real_space_tail_norm = ', result%metadata%real_space_tail_norm
      write(unit, '(a,es24.16)') '# real_space_total_norm = ', result%metadata%real_space_total_norm
      write(unit, '(a,es24.16)') '# real_space_tail_ratio = ', result%metadata%real_space_tail_ratio
      write(unit, '(a,es24.16)') '# real_space_last_shell_norm = ', result%metadata%real_space_last_shell_norm
      write(unit, '(a,es24.16)') '# real_space_tail_tolerance = ', result%metadata%real_space_tail_tolerance
      write(unit, '(a,l1)') '# real_space_tail_assessed = ', result%metadata%real_space_tail_assessed
      write(unit, '(a,l1)') '# real_space_source_covers_cutoff = ', result%metadata%real_space_source_covers_cutoff
      write(unit, '(a,i0)') '# profile_real_space_pair_response_integrations = ', &
         result%metadata%real_space_pair_response_integrations
      write(unit, '(a,es24.16)') '# profile_real_space_green_cpu_s = ', result%metadata%real_space_green_cpu_seconds
      write(unit, '(a,es24.16)') '# profile_real_space_pair_integration_cpu_s = ', &
         result%metadata%real_space_pair_integration_cpu_seconds
      write(unit, '(a,es24.16)') '# profile_real_space_fourier_cpu_s = ', result%metadata%real_space_fourier_cpu_seconds
      write(unit, '(a,es24.16)') '# profile_real_space_total_cpu_s = ', result%metadata%real_space_total_cpu_seconds
      write(unit, '(a,a)') '# energy_integration = ', trim(result%metadata%energy_integration)
      write(unit, '(a,es24.16)') '# eta_Ry = ', result%metadata%eta
      write(unit, '(a,es24.16)') '# fermi_level_Ry = ', result%metadata%fermi_level
      write(unit, '(a,a)') '# fermi_level_source = ', trim(result%metadata%fermi_source)
      write(unit, '(a,a)') '# fermi_level_policy = ', trim(result%metadata%fermi_policy)
      write(unit, '(a,es24.16)') '# electronic_temperature_K = ', result%metadata%electronic_temperature
      write(unit, '(a,es24.16)') '# electronic_kT_Ry = ', result%metadata%electronic_kT
      write(unit, '(a,3(1x,i0))') '# k_mesh_shape =', result%metadata%k_mesh_shape
      write(unit, '(a,i0)') '# nk = ', result%metadata%nk
      write(unit, '(a,2(1x,es24.16),1x,i0)') '# omega_grid_min_max_points = ', result%metadata%omega_min, &
         result%metadata%omega_max, result%metadata%omega_points
      write(unit, '(a,es24.16)') '# raw_k_weight_sum = ', result%metadata%k_weight_sum
      write(unit, '(a,i0)') '# available_band_count = ', result%metadata%available_band_count
      write(unit, '(a,2(1x,i0))') '# band_window_first_last = ', result%metadata%band_first, result%metadata%band_last
      write(unit, '(a,es24.16)') '# occupation_prune_tolerance = ', result%metadata%occupation_prune_tolerance
      write(unit, '(a,l1)') '# batched_accumulation = ', result%metadata%batched_accumulation
      write(unit, '(a,i0)') '# transition_batch_size = ', result%metadata%transition_batch_size
      write(unit, '(a,es24.16)') '# profile_transition_vertices_cpu_s = ', result%metadata%vertex_cpu_seconds
      write(unit, '(a,es24.16)') '# profile_transition_preparation_cpu_s = ', &
         result%metadata%transition_preparation_cpu_seconds
      write(unit, '(a,es24.16)') '# profile_frequency_denominators_cpu_s = ', result%metadata%denominator_cpu_seconds
      write(unit, '(a,es24.16)') '# profile_response_accumulation_cpu_s = ', result%metadata%accumulation_cpu_seconds
      write(unit, '(a,es24.16)') '# profile_arbitrary_kq_eigensolve_cpu_s = ', result%metadata%arbitrary_kq_cpu_seconds
      write(unit, '(a,es24.16)') '# profile_gf_energy_integration_cpu_s = ', &
         result%metadata%green_energy_integration_cpu_seconds
      write(unit, '(a,i0)') '# green_function_evaluations = ', result%metadata%green_function_evaluations
      write(unit, '(a,i0)') '# contour_points_per_segment = ', result%metadata%contour_points
      write(unit, '(a,i0)') '# contour_horizontal_subdivisions = ', result%metadata%contour_subdivisions
      write(unit, '(a,i0)') '# near_fermi_points = ', result%metadata%near_fermi_points
      write(unit, '(a,es24.16)') '# contour_height_Ry = ', result%metadata%contour_height
      write(unit, '(a,es24.16)') '# contour_max_imaginary_energy_Ry = ', result%metadata%contour_max_imaginary_energy
      write(unit, '(a,i0)') '# contour_gf_evaluations = ', result%metadata%contour_gf_evaluations
      write(unit, '(a,l1)') '# contour_deformation = ', result%metadata%contour_deformation
      if (result%metadata%integration_energy_points > 0) then
         write(unit, '(a,es24.16)') '# green_eta_Ry = ', result%metadata%green_eta
         write(unit, '(a,2(1x,es24.16))') '# integration_energy_window_Ry = ', &
            result%metadata%integration_energy_min, result%metadata%integration_energy_max
         write(unit, '(a,i0)') '# integration_energy_points = ', result%metadata%integration_energy_points
      end if
      write(unit, '(a)') '# record omega_Ry kind left_channel right_channel Re_chi_Ry^-1 Im_chi_Ry^-1 spectral_Ry^-1'
      do iw = 1, size(omega)
         do ileft = 1, size(result%chi, 1)
            do iright = 1, size(result%chi, 2)
               write(unit, '(es24.16,1x,a,2(1x,i0),3(1x,es24.16))') omega(iw), 'matrix', ileft, iright, &
                  result%re_chi(ileft, iright, iw), result%im_chi(ileft, iright, iw), &
                  -result%im_chi(ileft, iright, iw)/pi
            end do
         end do
         do idiag = 1, size(result%site_diagonal_spectrum, 1)
            write(unit, '(es24.16,1x,a,2(1x,i0),3(1x,es24.16))') omega(iw), 'site_diagonal', idiag, idiag, &
               0.0_rp, 0.0_rp, result%stoner_spectral_map(idiag, iw)
            write(unit, '(es24.16,1x,a,2(1x,i0),3(1x,es24.16))') omega(iw), 'stoner_landau', idiag, idiag, &
               0.0_rp, 0.0_rp, result%stoner_spectral_map(idiag, iw)
         end do
         write(unit, '(es24.16,1x,a,2(1x,i0),3(1x,es24.16))') omega(iw), 'trace', 0, 0, &
            0.0_rp, 0.0_rp, result%trace_spectrum(iw)
      end do
      close(unit)
   end subroutine write_chi_ks_text

   subroutine build_spectral_products(left_channels, right_channels, result)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_chi0_result), intent(inout) :: result
      integer :: nsite_diag, idiag, iw

      nsite_diag = min(size(left_channels), size(right_channels))
      allocate(result%site_diagonal_spectrum(nsite_diag, size(result%chi, 3)), &
         result%stoner_spectral_map(nsite_diag, size(result%chi, 3)), result%trace_spectrum(size(result%chi, 3)))
      result%site_diagonal_spectrum = 0.0_rp
      result%stoner_spectral_map = 0.0_rp
      do idiag = 1, nsite_diag
         if (left_channels(idiag)%site /= right_channels(idiag)%site) cycle
         do iw = 1, size(result%chi, 3)
            result%site_diagonal_spectrum(idiag, iw) = -aimag(result%chi(idiag, idiag, iw))/pi
            result%stoner_spectral_map(idiag, iw) = result%site_diagonal_spectrum(idiag, iw)
         end do
      end do
      do iw = 1, size(result%chi, 3)
         result%trace_spectrum(iw) = -sum(aimag([(result%chi(idiag, idiag, iw), idiag=1, nsite_diag)]))/pi
      end do
   end subroutine build_spectral_products

   subroutine validate_chi_ks_inputs(nk, nbands, nspinor, nleft, nright, nw, k_weights, eigenvalues_k, &
      eigenvectors_k, eigenvalues_kq, eigenvectors_kq, options, require_eta)
      integer, intent(in) :: nk, nbands, nspinor, nleft, nright, nw
      real(rp), intent(in) :: k_weights(:), eigenvalues_k(:, :), eigenvalues_kq(:, :)
      complex(rp), intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :)
      type(tddft_chi0_options), intent(in) :: options
      logical, intent(in), optional :: require_eta
      logical :: eta_required

      eta_required = .true.
      if (present(require_eta)) eta_required = require_eta
      call validate_tddft_chi0_options(options, 'build_chi_ks_from_eigenpairs')
      if (nk <= 0 .or. nbands <= 0 .or. nspinor <= 0 .or. nleft <= 0 .or. nright <= 0 .or. nw <= 0) then
         error stop 'build_chi_ks_from_eigenpairs: empty input is not valid'
      end if
      if (size(eigenvalues_k, 2) /= nk .or. any(shape(eigenvalues_kq) /= shape(eigenvalues_k))) then
         error stop 'build_chi_ks_from_eigenpairs: eigenvalue shapes are incompatible'
      end if
      if (size(eigenvectors_k, 1) /= nspinor .or. size(eigenvectors_k, 2) /= nbands .or. &
          size(eigenvectors_k, 3) /= nk .or. any(shape(eigenvectors_kq) /= shape(eigenvectors_k))) then
         error stop 'build_chi_ks_from_eigenpairs: eigenvector shapes are incompatible'
      end if
      if (any(k_weights < 0.0_rp) .or. sum(k_weights) <= tiny(1.0_rp)) then
         error stop 'build_chi_ks_from_eigenpairs: k weights must have a positive sum'
      end if
      if (options%eta < 0.0_rp .or. (eta_required .and. options%eta <= 0.0_rp)) then
         if (eta_required) then
            error stop 'build_chi_ks_from_eigenpairs: eta must be positive'
         else
            error stop 'build_chi_ks_from_eigenpairs: eta must be non-negative'
         end if
      end if
      if (options%electronic_temperature < 0.0_rp) then
         error stop 'build_chi_ks_from_eigenpairs: electronic temperature must be non-negative'
      end if
      if (options%occupation_prune_tolerance < 0.0_rp) then
         error stop 'build_chi_ks_from_eigenpairs: occupation pruning tolerance must be non-negative'
      end if
      if (options%transition_batch_size < 1) then
         error stop 'build_chi_ks_from_eigenpairs: transition_batch_size must be positive'
      end if
   end subroutine validate_chi_ks_inputs

end module tddft_chi0_mod
