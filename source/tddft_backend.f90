!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Common chi0 backend contract and response-grid data model.
!>
!> The contract is deliberately expressed in terms of response data, not in
!> terms of a particular one-electron storage object.  An implementation may
!> use eigenpairs, a K-space Lehmann Green function, or a native real-space
!> Green function provider.  In particular, no G(R,z)->G(k,z) conversion is
!> implied by this module.
module tddft_backend_mod
   use precision_mod, only: rp
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, tddft_chi0_request, &
      tddft_chi0_batch_result, tddft_chi0_metadata, build_chi_ks_from_eigenpairs, &
      build_static_chi_ks_from_eigenpairs_at_q
   use tddft_chi0_green_mod, only: green_chi0_options, eigenpair_green_function_provider, &
      build_chi_ks_from_green_functions, build_static_chi_ks_from_green_functions
   use response_vertices_mod, only: response_channel
   implicit none

   private

   character(len=*), parameter, public :: TDDFT_BACKEND_EIGENPAIRS = 'eigenpairs'
   character(len=*), parameter, public :: TDDFT_BACKEND_KSPACE_LEHMANN = 'kspace_lehmann'
   character(len=*), parameter, public :: TDDFT_BACKEND_REALSPACE_GF = 'realspace_gf'
   character(len=*), parameter, public :: TDDFT_BACKEND_MOCK = 'mock'

   !> Capabilities are data so the common driver can reject an unsupported
   !> request before allocating a large response tensor.
   type, public :: tddft_backend_capabilities
      logical :: supports_single_point = .true.
      logical :: supports_frequency_batch = .true.
      logical :: supports_q_batch = .true.
      logical :: supports_grid_batch = .true.
      logical :: supports_static_limit = .false.
      logical :: native_real_space = .false.
      logical :: reuses_real_space_response = .false.
      character(len=128) :: unsupported_reason = ''
   end type tddft_backend_capabilities

   !> A small selector record keeps input aliases out of backend code.
   type, public :: tddft_backend_selection
      character(len=32) :: requested_name = ''
      character(len=32) :: canonical_name = ''
      logical :: compatibility_alias = .false.
   end type tddft_backend_selection

   !> Abstract common chi0 provider.  `evaluate` receives a full q/omega grid;
   !> the convenience methods below expose the common scalar and batch forms
   !> without forcing callers to construct request envelopes themselves.
   type, abstract, public :: tddft_chi0_backend
   contains
      procedure(tddft_backend_evaluate), deferred :: evaluate
      procedure(tddft_backend_name), deferred :: backend_name
      procedure(tddft_backend_get_capabilities), deferred :: capabilities
      procedure :: evaluate_one => backend_evaluate_one
      procedure :: evaluate_frequency_batch => backend_evaluate_frequency_batch
      procedure :: evaluate_q_batch => backend_evaluate_q_batch
      procedure :: evaluate_grid => backend_evaluate_grid
      ! Static response is a separate operation.  Implementations must not
      ! obtain it by calling evaluate with omega=0 and finite eta.
      procedure :: evaluate_static_grid => backend_evaluate_static_grid
   end type tddft_chi0_backend

   !> Provider-side contract for native real-space response implementations.
   !> The provider returns chi0(R,omega) transformed to the requested q points
   !> (or a natural real-space/mixed representation) directly.  It never needs
   !> to expose a reciprocal-space Green function to the common layer.
   type, abstract, public :: tddft_realspace_chi0_provider
   contains
      procedure(tddft_realspace_evaluate), deferred :: evaluate_realspace
      procedure :: evaluate_static_realspace => unsupported_realspace_static
      procedure(tddft_realspace_describe), deferred :: describe
   end type tddft_realspace_chi0_provider

   !> Eigenpair adapter used by the reference production backend.  The q-end
   !> point arrays are retained as explicit data, so exact endpoint provenance
   !> is not lost when a q batch is evaluated.
   type, extends(tddft_chi0_backend), public :: tddft_eigenpair_backend
      real(rp), allocatable :: k_weights(:)
      real(rp), allocatable :: eigenvalues_k(:, :)
      real(rp), allocatable :: eigenvalues_kq(:, :, :)
      complex(rp), allocatable :: eigenvectors_k(:, :, :)
      complex(rp), allocatable :: eigenvectors_kq(:, :, :, :)
      real(rp), allocatable :: q_points(:, :)
      integer, allocatable :: site_orbital_counts(:)
      type(response_channel), allocatable :: left_channels(:), right_channels(:)
      type(tddft_chi0_options) :: options
      type(green_chi0_options) :: green_options
   contains
      procedure :: initialize => initialize_eigenpair_backend
      procedure :: evaluate => evaluate_eigenpair_backend
      procedure :: backend_name => eigenpair_backend_name
      procedure :: capabilities => eigenpair_backend_capabilities
   end type tddft_eigenpair_backend

   !> K-space Lehmann adapter.  It uses the same endpoint data as the
   !> eigenpair adapter, but exercises the one-particle GF bubble independently
   !> of the explicit transition denominator implementation.
   type, extends(tddft_eigenpair_backend), public :: tddft_kspace_lehmann_backend
   contains
      procedure :: evaluate => evaluate_kspace_lehmann_backend
      procedure :: backend_name => kspace_lehmann_backend_name
      procedure :: capabilities => kspace_lehmann_backend_capabilities
   end type tddft_kspace_lehmann_backend

   !> Native real-space adapter.  It can be selected before a provider is
   !> available, but evaluation then fails explicitly instead of silently
   !> converting through k space or using the wrong Green-function route.
   type, extends(tddft_chi0_backend), public :: tddft_realspace_gf_backend
      class(tddft_realspace_chi0_provider), allocatable :: provider
   contains
      procedure :: initialize => initialize_realspace_gf_backend
      procedure :: evaluate => evaluate_realspace_gf_backend
      procedure :: backend_name => realspace_gf_backend_name
      procedure :: capabilities => realspace_gf_backend_capabilities
   end type tddft_realspace_gf_backend

   !> Deterministic fixture backend.  It is intentionally small and has no
   !> electronic-structure assumptions; it verifies that Dyson and spectral
   !> consumers depend only on the common response result.
   type, extends(tddft_chi0_backend), public :: tddft_mock_chi0_backend
      integer :: response_dimension = 1
      complex(rp) :: scale = cmplx(1.0_rp, 0.0_rp, rp)
   contains
      procedure :: initialize => initialize_mock_backend
      procedure :: evaluate => evaluate_mock_backend
      procedure :: backend_name => mock_backend_name
      procedure :: capabilities => mock_backend_capabilities
   end type tddft_mock_chi0_backend

   abstract interface
      subroutine tddft_backend_evaluate(this, request, result)
         import :: tddft_chi0_backend, tddft_chi0_request, tddft_chi0_batch_result
         class(tddft_chi0_backend), intent(inout) :: this
         type(tddft_chi0_request), intent(in) :: request
         type(tddft_chi0_batch_result), intent(out) :: result
      end subroutine tddft_backend_evaluate

      function tddft_backend_name(this) result(name)
         import :: tddft_chi0_backend
         class(tddft_chi0_backend), intent(in) :: this
         character(len=32) :: name
      end function tddft_backend_name

      subroutine tddft_backend_get_capabilities(this, capabilities)
         import :: tddft_backend_capabilities, tddft_chi0_backend
         class(tddft_chi0_backend), intent(in) :: this
         type(tddft_backend_capabilities), intent(out) :: capabilities
      end subroutine tddft_backend_get_capabilities

      subroutine tddft_realspace_evaluate(this, request, result)
         import :: tddft_realspace_chi0_provider, tddft_chi0_request, tddft_chi0_batch_result
         class(tddft_realspace_chi0_provider), intent(inout) :: this
         type(tddft_chi0_request), intent(in) :: request
         type(tddft_chi0_batch_result), intent(out) :: result
      end subroutine tddft_realspace_evaluate

      subroutine tddft_realspace_describe(this, metadata)
         import :: tddft_backend_capabilities, tddft_realspace_chi0_provider
         class(tddft_realspace_chi0_provider), intent(in) :: this
         type(tddft_backend_capabilities), intent(out) :: metadata
      end subroutine tddft_realspace_describe
   end interface

   public :: canonical_tddft_backend_name
   public :: select_tddft_chi0_backend
   public :: make_tddft_chi0_backend

contains

   !> Normalize public input names.  `green` and `kspace_gf` are retained as
   !> compatibility aliases for the pre-interface branch; new input should use
   !> `kspace_lehmann`.  `realspace` and `rs_gf` are similarly explicit aliases.
   pure function canonical_tddft_backend_name(requested_name) result(canonical_name)
      character(len=*), intent(in) :: requested_name
      character(len=32) :: canonical_name
      character(len=64) :: name

      name = lower_ascii(trim(requested_name))
      select case (trim(name))
      case (TDDFT_BACKEND_EIGENPAIRS)
         canonical_name = TDDFT_BACKEND_EIGENPAIRS
      case (TDDFT_BACKEND_KSPACE_LEHMANN, 'kspace_gf', 'green')
         canonical_name = TDDFT_BACKEND_KSPACE_LEHMANN
      case (TDDFT_BACKEND_REALSPACE_GF, 'realspace', 'rs_gf', 'real_space_gf')
         canonical_name = TDDFT_BACKEND_REALSPACE_GF
      case (TDDFT_BACKEND_MOCK)
         canonical_name = TDDFT_BACKEND_MOCK
      case default
         canonical_name = ''
      end select
   end function canonical_tddft_backend_name

   pure function select_tddft_chi0_backend(requested_name) result(selection)
      character(len=*), intent(in) :: requested_name
      type(tddft_backend_selection) :: selection

      selection%requested_name = trim(requested_name)
      selection%canonical_name = canonical_tddft_backend_name(requested_name)
      selection%compatibility_alias = len_trim(selection%canonical_name) > 0 .and. &
         trim(lower_ascii(requested_name)) /= trim(selection%canonical_name)
   end function select_tddft_chi0_backend

   !> Allocate the selected backend.  Initialization is deliberately separate:
   !> the factory selects a contract, while a backend adapter owns its source
   !> data and can be initialized with the current calculation response mesh.
   subroutine make_tddft_chi0_backend(requested_name, backend)
      character(len=*), intent(in) :: requested_name
      class(tddft_chi0_backend), allocatable, intent(out) :: backend
      type(tddft_backend_selection) :: selection

      selection = select_tddft_chi0_backend(requested_name)
      if (len_trim(selection%canonical_name) == 0) then
         error stop 'make_tddft_chi0_backend: unknown chi0 backend'
      end if
      select case (trim(selection%canonical_name))
      case (TDDFT_BACKEND_EIGENPAIRS)
         allocate(tddft_eigenpair_backend :: backend)
      case (TDDFT_BACKEND_KSPACE_LEHMANN)
         allocate(tddft_kspace_lehmann_backend :: backend)
      case (TDDFT_BACKEND_REALSPACE_GF)
         allocate(tddft_realspace_gf_backend :: backend)
      case (TDDFT_BACKEND_MOCK)
         allocate(tddft_mock_chi0_backend :: backend)
      case default
         error stop 'make_tddft_chi0_backend: selector returned an invalid backend'
      end select
   end subroutine make_tddft_chi0_backend

   subroutine backend_evaluate_one(this, q_point, omega, result, q_index)
      class(tddft_chi0_backend), intent(inout) :: this
      real(rp), intent(in) :: q_point(3), omega
      type(tddft_chi0_result), intent(out) :: result
      integer, intent(in), optional :: q_index
      type(tddft_chi0_request) :: request
      type(tddft_chi0_batch_result) :: batch

      allocate(request%q_points(3, 1), request%omega(1))
      request%q_points(:, 1) = q_point
      request%omega(1) = omega
      if (present(q_index)) then
         allocate(request%q_indices(1))
         request%q_indices(1) = q_index
      end if
      request%batch_mode = 'single_point'
      call this%evaluate(request, batch)
      call require_batch_shape(batch, 1, 1, 'backend_evaluate_one')
      result = batch%q_response(1)
   end subroutine backend_evaluate_one

   subroutine backend_evaluate_frequency_batch(this, q_point, omega, result, q_index)
      class(tddft_chi0_backend), intent(inout) :: this
      real(rp), intent(in) :: q_point(3), omega(:)
      type(tddft_chi0_result), intent(out) :: result
      integer, intent(in), optional :: q_index
      type(tddft_chi0_request) :: request
      type(tddft_chi0_batch_result) :: batch

      allocate(request%q_points(3, 1), request%omega(size(omega)))
      request%q_points(:, 1) = q_point
      request%omega = omega
      if (present(q_index)) then
         allocate(request%q_indices(1))
         request%q_indices(1) = q_index
      end if
      request%batch_mode = 'frequency_batch'
      call this%evaluate(request, batch)
      call require_batch_shape(batch, 1, size(omega), 'backend_evaluate_frequency_batch')
      result = batch%q_response(1)
   end subroutine backend_evaluate_frequency_batch

   subroutine backend_evaluate_q_batch(this, q_points, omega, result, q_indices)
      class(tddft_chi0_backend), intent(inout) :: this
      real(rp), intent(in) :: q_points(:, :), omega
      type(tddft_chi0_batch_result), intent(out) :: result
      integer, intent(in), optional :: q_indices(:)
      type(tddft_chi0_request) :: request

      if (size(q_points, 1) /= 3 .or. size(q_points, 2) < 1) then
         error stop 'backend_evaluate_q_batch: q_points must have shape (3,nq)'
      end if
      allocate(request%q_points(3, size(q_points, 2)), request%omega(1))
      request%q_points = q_points
      request%omega(1) = omega
      if (present(q_indices)) then
         if (size(q_indices) /= size(q_points, 2)) error stop 'backend_evaluate_q_batch: q index count mismatch'
         allocate(request%q_indices(size(q_indices)))
         request%q_indices = q_indices
      end if
      request%batch_mode = 'q_batch'
      call this%evaluate(request, result)
      call require_batch_shape(result, size(q_points, 2), 1, 'backend_evaluate_q_batch')
   end subroutine backend_evaluate_q_batch

   subroutine backend_evaluate_grid(this, q_points, omega, result, q_indices)
      class(tddft_chi0_backend), intent(inout) :: this
      real(rp), intent(in) :: q_points(:, :), omega(:)
      type(tddft_chi0_batch_result), intent(out) :: result
      integer, intent(in), optional :: q_indices(:)
      type(tddft_chi0_request) :: request

      if (size(q_points, 1) /= 3 .or. size(q_points, 2) < 1 .or. size(omega) < 1) then
         error stop 'backend_evaluate_grid: invalid q/omega grid'
      end if
      allocate(request%q_points(3, size(q_points, 2)), request%omega(size(omega)))
      request%q_points = q_points
      request%omega = omega
      if (present(q_indices)) then
         if (size(q_indices) /= size(q_points, 2)) error stop 'backend_evaluate_grid: q index count mismatch'
         allocate(request%q_indices(size(q_indices)))
         request%q_indices = q_indices
      end if
      request%batch_mode = 'q_and_omega'
      call this%evaluate(request, result)
      call require_batch_shape(result, size(q_points, 2), size(omega), 'backend_evaluate_grid')
   end subroutine backend_evaluate_grid

   !> Evaluate the genuine zero-frequency response on a q batch.
   !>
   !> This entry point is intentionally separate from `evaluate`: a static
   !> implementation uses either the divided-difference Lehmann expression or
   !> the native real-space retarded/advanced contour identity.  In
   !> particular, omega=0 is never forwarded to a finite-eta dynamic route.
   subroutine backend_evaluate_static_grid(this, q_points, result, q_indices)
      class(tddft_chi0_backend), intent(inout) :: this
      real(rp), intent(in) :: q_points(:, :)
      type(tddft_chi0_batch_result), intent(out) :: result
      integer, intent(in), optional :: q_indices(:)
      type(tddft_chi0_request) :: request
      type(eigenpair_green_function_provider) :: source
      type(tddft_chi0_options) :: static_options
      type(green_chi0_options) :: green_options
      integer :: iq, q_index, nq

      if (size(q_points, 1) /= 3 .or. size(q_points, 2) < 1) then
         error stop 'backend_evaluate_static_grid: q_points must have shape (3,nq)'
      end if
      nq = size(q_points, 2)
      if (present(q_indices)) then
         if (size(q_indices) /= nq) error stop 'backend_evaluate_static_grid: q index count mismatch'
      end if
      allocate(request%q_points(3, nq), request%q_indices(nq), request%omega(1))
      request%q_points = q_points
      request%omega(1) = 0.0_rp
      do iq = 1, nq
         if (present(q_indices)) then
            request%q_indices(iq) = q_indices(iq)
         else
            request%q_indices(iq) = iq
         end if
      end do
      request%batch_mode = 'static_q_batch'

      ! The most-derived type must be selected first because the K-space
      ! Lehmann backend extends the eigenpair backend.
      select type (this)
      type is (tddft_kspace_lehmann_backend)
         if (.not. allocated(this%eigenvalues_k)) then
            error stop 'backend_evaluate_static_grid: K-space backend is not initialized'
         end if
         allocate(result%q_response(nq), result%q_points(3, nq), result%q_indices(nq), result%omega(1))
         result%q_points = q_points; result%q_indices = request%q_indices; result%omega = 0.0_rp
         green_options = this%green_options
         do iq = 1, nq
            q_index = request%q_indices(iq)
            if (q_index < 1 .or. q_index > size(this%q_points, 2)) then
               error stop 'backend_evaluate_static_grid: K-space q index is outside initialized endpoint data'
            end if
            green_options%q_direct = q_points(:, iq)
            green_options%fermi_level = this%options%fermi_level
            green_options%electronic_temperature = this%options%electronic_temperature
            green_options%k_mesh_shape = this%options%k_mesh_shape
            call source%initialize(this%eigenvalues_k, this%eigenvectors_k, this%eigenvalues_kq(:, :, q_index), &
               this%eigenvectors_kq(:, :, :, q_index))
            call build_static_chi_ks_from_green_functions(source, this%k_weights, this%site_orbital_counts, &
               this%left_channels, this%right_channels, green_options, result%q_response(iq))
            call annotate_static_result(result%q_response(iq), 'kspace_lehmann', q_points(:, iq), nq)
         end do
      type is (tddft_eigenpair_backend)
         if (.not. allocated(this%eigenvalues_k)) then
            error stop 'backend_evaluate_static_grid: eigenpair backend is not initialized'
         end if
         allocate(result%q_response(nq), result%q_points(3, nq), result%q_indices(nq), result%omega(1))
         result%q_points = q_points; result%q_indices = request%q_indices; result%omega = 0.0_rp
         do iq = 1, nq
            q_index = request%q_indices(iq)
            if (q_index < 1 .or. q_index > size(this%q_points, 2)) then
               error stop 'backend_evaluate_static_grid: eigenpair q index is outside initialized endpoint data'
            end if
            static_options = this%options
            static_options%eta = 0.0_rp
            static_options%q_direct = q_points(:, iq)
            call build_static_chi_ks_from_eigenpairs_at_q(this%k_weights, this%eigenvalues_k, this%eigenvectors_k, &
               this%eigenvalues_kq(:, :, q_index), this%eigenvectors_kq(:, :, :, q_index), this%site_orbital_counts, &
               this%left_channels, this%right_channels, static_options, result%q_response(iq))
            call annotate_static_result(result%q_response(iq), 'eigenpairs', q_points(:, iq), nq)
         end do
      type is (tddft_realspace_gf_backend)
         if (.not. allocated(this%provider)) then
            error stop 'backend_evaluate_static_grid: native real-space provider is absent'
         end if
         call this%provider%evaluate_static_realspace(request, result)
      class default
         error stop 'backend_evaluate_static_grid: selected backend has no exact static implementation'
      end select

      if (allocated(result%q_response) .and. size(result%q_response) > 0) then
         result%metadata = result%q_response(1)%metadata
         result%metadata%q_batch_size = nq
         result%metadata%omega_batch_size = 1
      end if
   end subroutine backend_evaluate_static_grid

   subroutine unsupported_realspace_static(this, request, result)
      class(tddft_realspace_chi0_provider), intent(inout) :: this
      type(tddft_chi0_request), intent(in) :: request
      type(tddft_chi0_batch_result), intent(out) :: result

      error stop 'native real-space provider: exact static chi0 is unavailable'
   end subroutine unsupported_realspace_static

   subroutine annotate_static_result(result, backend, q_point, nq)
      type(tddft_chi0_result), intent(inout) :: result
      character(len=*), intent(in) :: backend
      real(rp), intent(in) :: q_point(3)
      integer, intent(in) :: nq

      result%metadata%backend = trim(backend)
      result%metadata%canonical_backend = trim(backend)
      result%metadata%q_direct = q_point
      result%metadata%omega_min = 0.0_rp
      result%metadata%omega_max = 0.0_rp
      result%metadata%omega_points = 1
      result%metadata%eta = 0.0_rp
      result%metadata%green_eta = 0.0_rp
      result%metadata%eta_role = 'not used by exact static response'
      result%metadata%eta_is_numerical = .false.
      result%metadata%frequency_convention = 'static omega=0 divided difference; no dynamical eta'
      result%metadata%static_limit = .true.
      result%metadata%q_batch_size = nq
      result%metadata%omega_batch_size = 1
   end subroutine annotate_static_result

   subroutine initialize_eigenpair_backend(this, k_weights, eigenvalues_k, eigenvectors_k, eigenvalues_kq, &
      eigenvectors_kq, q_points, site_orbital_counts, left_channels, right_channels, options, green_options)
      class(tddft_eigenpair_backend), intent(inout) :: this
      real(rp), intent(in) :: k_weights(:), eigenvalues_k(:, :), eigenvalues_kq(:, :, :), q_points(:, :)
      complex(rp), intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :, :)
      integer, intent(in) :: site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_chi0_options), intent(in) :: options
      type(green_chi0_options), intent(in), optional :: green_options
      integer :: nbands, nk, nspinor, nq

      nbands = size(eigenvalues_k, 1); nk = size(eigenvalues_k, 2); nspinor = size(eigenvectors_k, 1)
      nq = size(q_points, 2)
      if (nk < 1 .or. nbands < 1 .or. nspinor < 1 .or. size(eigenvalues_kq, 1) /= nbands .or. &
          size(eigenvalues_kq, 2) /= nk .or. size(eigenvalues_kq, 3) /= nq .or. size(eigenvectors_k, 2) /= nbands .or. &
          size(eigenvectors_k, 3) /= nk .or. size(eigenvectors_kq, 1) /= nspinor .or. size(eigenvectors_kq, 2) /= nbands .or. &
          size(eigenvectors_kq, 3) /= nk .or. size(eigenvectors_kq, 4) /= nq .or. size(k_weights) /= nk .or. &
          size(q_points, 1) /= 3 .or. size(site_orbital_counts) < 1 .or. size(left_channels) < 1 .or. &
          size(left_channels) /= size(right_channels)) then
         error stop 'initialize_eigenpair_backend: incompatible endpoint or response arrays'
      end if
      if (any(k_weights < 0.0_rp) .or. sum(k_weights) <= tiny(1.0_rp)) then
         error stop 'initialize_eigenpair_backend: k weights must be non-negative with positive sum'
      end if
      this%k_weights = k_weights
      this%eigenvalues_k = eigenvalues_k
      this%eigenvalues_kq = eigenvalues_kq
      this%eigenvectors_k = eigenvectors_k
      this%eigenvectors_kq = eigenvectors_kq
      this%q_points = q_points
      this%site_orbital_counts = site_orbital_counts
      this%left_channels = left_channels
      this%right_channels = right_channels
      this%options = options
      if (present(green_options)) this%green_options = green_options
   end subroutine initialize_eigenpair_backend

   subroutine evaluate_eigenpair_backend(this, request, result)
      class(tddft_eigenpair_backend), intent(inout) :: this
      type(tddft_chi0_request), intent(in) :: request
      type(tddft_chi0_batch_result), intent(out) :: result
      type(tddft_chi0_options) :: options
      integer :: iq, q_index, nq

      call validate_request(request)
      if (.not. allocated(this%eigenvalues_k)) error stop 'eigenpair backend: backend is not initialized'
      nq = size(request%q_points, 2)
      allocate(result%q_response(nq), result%q_points(3, nq), result%q_indices(nq), result%omega(size(request%omega)))
      result%q_points = request%q_points
      result%omega = request%omega
      do iq = 1, nq
         q_index = request_q_index(request, iq)
         if (q_index < 1 .or. q_index > size(this%q_points, 2)) then
            error stop 'eigenpair backend: q index is outside initialized endpoint data'
         end if
         result%q_indices(iq) = q_index
         options = this%options
         options%q_direct = request%q_points(:, iq)
         call build_chi_ks_from_eigenpairs(this%k_weights, this%eigenvalues_k, this%eigenvectors_k, &
            this%eigenvalues_kq(:, :, q_index), this%eigenvectors_kq(:, :, :, q_index), this%site_orbital_counts, &
            this%left_channels, this%right_channels, request%omega, options, result%q_response(iq))
         call annotate_result(result%q_response(iq), TDDFT_BACKEND_EIGENPAIRS, 'explicit eigenpair transitions', &
            request%q_points(:, iq), nq, size(request%omega))
      end do
      call annotate_batch(result, TDDFT_BACKEND_EIGENPAIRS)
   end subroutine evaluate_eigenpair_backend

   function eigenpair_backend_name(this) result(name)
      class(tddft_eigenpair_backend), intent(in) :: this
      character(len=32) :: name
      name = TDDFT_BACKEND_EIGENPAIRS
   end function eigenpair_backend_name

   subroutine eigenpair_backend_capabilities(this, capabilities)
      class(tddft_eigenpair_backend), intent(in) :: this
      type(tddft_backend_capabilities), intent(out) :: capabilities
      capabilities = tddft_backend_capabilities()
      capabilities%supports_static_limit = .true.
   end subroutine eigenpair_backend_capabilities

   subroutine evaluate_kspace_lehmann_backend(this, request, result)
      class(tddft_kspace_lehmann_backend), intent(inout) :: this
      type(tddft_chi0_request), intent(in) :: request
      type(tddft_chi0_batch_result), intent(out) :: result
      type(eigenpair_green_function_provider) :: source
      type(green_chi0_options) :: options
      integer :: iq, q_index, nq

      call validate_request(request)
      if (.not. allocated(this%eigenvalues_k)) error stop 'K-space Lehmann backend: backend is not initialized'
      nq = size(request%q_points, 2)
      allocate(result%q_response(nq), result%q_points(3, nq), result%q_indices(nq), result%omega(size(request%omega)))
      result%q_points = request%q_points
      result%omega = request%omega
      options = this%green_options
      do iq = 1, nq
         q_index = request_q_index(request, iq)
         if (q_index < 1 .or. q_index > size(this%q_points, 2)) then
            error stop 'K-space Lehmann backend: q index is outside initialized endpoint data'
         end if
         result%q_indices(iq) = q_index
         options%q_direct = request%q_points(:, iq)
         call source%initialize(this%eigenvalues_k, this%eigenvectors_k, this%eigenvalues_kq(:, :, q_index), &
            this%eigenvectors_kq(:, :, :, q_index))
         options%fermi_level = this%options%fermi_level
         options%electronic_temperature = this%options%electronic_temperature
         options%k_mesh_shape = this%options%k_mesh_shape
         call build_chi_ks_from_green_functions(source, this%k_weights, this%site_orbital_counts, &
            this%left_channels, this%right_channels, request%omega, options, result%q_response(iq))
         call annotate_result(result%q_response(iq), TDDFT_BACKEND_KSPACE_LEHMANN, 'K-space Lehmann GF bubble', &
            request%q_points(:, iq), nq, size(request%omega))
         result%q_response(iq)%metadata%endpoint_provenance = 'K-space Lehmann resolvents at k and k+q endpoints'
      end do
      call annotate_batch(result, TDDFT_BACKEND_KSPACE_LEHMANN)
   end subroutine evaluate_kspace_lehmann_backend

   function kspace_lehmann_backend_name(this) result(name)
      class(tddft_kspace_lehmann_backend), intent(in) :: this
      character(len=32) :: name
      name = TDDFT_BACKEND_KSPACE_LEHMANN
   end function kspace_lehmann_backend_name

   subroutine kspace_lehmann_backend_capabilities(this, capabilities)
      class(tddft_kspace_lehmann_backend), intent(in) :: this
      type(tddft_backend_capabilities), intent(out) :: capabilities
      capabilities = tddft_backend_capabilities()
      capabilities%supports_static_limit = .true.
   end subroutine kspace_lehmann_backend_capabilities

   subroutine initialize_realspace_gf_backend(this, provider)
      class(tddft_realspace_gf_backend), intent(inout) :: this
      class(tddft_realspace_chi0_provider), intent(in) :: provider

      if (allocated(this%provider)) deallocate(this%provider)
      allocate(this%provider, source=provider)
   end subroutine initialize_realspace_gf_backend

   subroutine evaluate_realspace_gf_backend(this, request, result)
      class(tddft_realspace_gf_backend), intent(inout) :: this
      type(tddft_chi0_request), intent(in) :: request
      type(tddft_chi0_batch_result), intent(out) :: result
      integer :: iq

      call validate_request(request)
      if (.not. allocated(this%provider)) then
         error stop 'native real-space GF backend: no real-space chi0 provider is attached'
      end if
      call this%provider%evaluate_realspace(request, result)
      call require_batch_shape(result, size(request%q_points, 2), size(request%omega), 'native real-space provider')
      if (.not. allocated(result%q_indices)) then
         allocate(result%q_indices(size(request%q_points, 2)))
         do iq = 1, size(request%q_points, 2)
            result%q_indices(iq) = iq
         end do
      end if
      call annotate_batch(result, TDDFT_BACKEND_REALSPACE_GF)
      result%metadata%implementation = 'native real-space GF -> chi0(R,w) -> requested representation'
      result%metadata%real_space_reuse = request%allow_real_space_reuse
      do iq = 1, size(result%q_response)
         result%q_response(iq)%metadata%backend = TDDFT_BACKEND_REALSPACE_GF
         result%q_response(iq)%metadata%canonical_backend = TDDFT_BACKEND_REALSPACE_GF
         result%q_response(iq)%metadata%implementation = result%metadata%implementation
         result%q_response(iq)%metadata%real_space_reuse = request%allow_real_space_reuse
      end do
   end subroutine evaluate_realspace_gf_backend

   function realspace_gf_backend_name(this) result(name)
      class(tddft_realspace_gf_backend), intent(in) :: this
      character(len=32) :: name
      name = TDDFT_BACKEND_REALSPACE_GF
   end function realspace_gf_backend_name

   subroutine realspace_gf_backend_capabilities(this, capabilities)
      class(tddft_realspace_gf_backend), intent(in) :: this
      type(tddft_backend_capabilities), intent(out) :: capabilities
      type(tddft_backend_capabilities) :: provider_capabilities

      capabilities = tddft_backend_capabilities()
      capabilities%native_real_space = .true.
      capabilities%reuses_real_space_response = .true.
      if (.not. allocated(this%provider)) then
         capabilities%unsupported_reason = 'native real-space provider is not initialized'
         return
      end if
      call this%provider%describe(provider_capabilities)
      capabilities%supports_static_limit = provider_capabilities%supports_static_limit
      if (.not. capabilities%supports_static_limit) then
         capabilities%unsupported_reason = 'attached native real-space provider has no exact static operation'
      end if
   end subroutine realspace_gf_backend_capabilities

   subroutine initialize_mock_backend(this, response_dimension, scale)
      class(tddft_mock_chi0_backend), intent(inout) :: this
      integer, intent(in) :: response_dimension
      complex(rp), intent(in), optional :: scale

      if (response_dimension < 1) error stop 'mock chi0 backend: response dimension must be positive'
      this%response_dimension = response_dimension
      if (present(scale)) this%scale = scale
   end subroutine initialize_mock_backend

   subroutine evaluate_mock_backend(this, request, result)
      class(tddft_mock_chi0_backend), intent(inout) :: this
      type(tddft_chi0_request), intent(in) :: request
      type(tddft_chi0_batch_result), intent(out) :: result
      integer :: iq, iw, ich, nq, nw

      call validate_request(request)
      nq = size(request%q_points, 2); nw = size(request%omega)
      allocate(result%q_response(nq), result%q_points(3, nq), result%q_indices(nq), result%omega(nw))
      result%q_points = request%q_points; result%omega = request%omega
      do iq = 1, nq
         result%q_indices(iq) = request_q_index(request, iq)
         allocate(result%q_response(iq)%chi(this%response_dimension, this%response_dimension, nw), &
            result%q_response(iq)%re_chi(this%response_dimension, this%response_dimension, nw), &
            result%q_response(iq)%im_chi(this%response_dimension, this%response_dimension, nw), &
            result%q_response(iq)%trace_spectrum(nw), result%q_response(iq)%site_diagonal_spectrum(this%response_dimension, nw))
         result%q_response(iq)%chi = cmplx(0.0_rp, 0.0_rp, rp)
         do iw = 1, nw
            do ich = 1, this%response_dimension
               result%q_response(iq)%chi(ich, ich, iw) = this%scale*cmplx(real(iq, rp), request%omega(iw), rp)
               result%q_response(iq)%site_diagonal_spectrum(ich, iw) = &
                  -aimag(result%q_response(iq)%chi(ich, ich, iw))/real(max(1, this%response_dimension), rp)
            end do
            result%q_response(iq)%trace_spectrum(iw) = sum(result%q_response(iq)%site_diagonal_spectrum(:, iw))
         end do
         result%q_response(iq)%re_chi = real(result%q_response(iq)%chi, rp)
         result%q_response(iq)%im_chi = aimag(result%q_response(iq)%chi)
         call annotate_result(result%q_response(iq), TDDFT_BACKEND_MOCK, 'deterministic test fixture', &
            request%q_points(:, iq), nq, nw)
      end do
      call annotate_batch(result, TDDFT_BACKEND_MOCK)
      result%metadata%convergence_status = 'mock fixture'
      result%metadata%converged = .true.
   end subroutine evaluate_mock_backend

   function mock_backend_name(this) result(name)
      class(tddft_mock_chi0_backend), intent(in) :: this
      character(len=32) :: name
      name = TDDFT_BACKEND_MOCK
   end function mock_backend_name

   subroutine mock_backend_capabilities(this, capabilities)
      class(tddft_mock_chi0_backend), intent(in) :: this
      type(tddft_backend_capabilities), intent(out) :: capabilities
      capabilities = tddft_backend_capabilities()
      capabilities%unsupported_reason = 'mock backend has no exact static implementation'
   end subroutine mock_backend_capabilities

   subroutine validate_request(request)
      type(tddft_chi0_request), intent(in) :: request

      if (.not. allocated(request%q_points) .or. .not. allocated(request%omega)) then
         error stop 'chi0 backend request: q_points and omega must be allocated'
      end if
      if (size(request%q_points, 1) /= 3 .or. size(request%q_points, 2) < 1 .or. size(request%omega) < 1) then
         error stop 'chi0 backend request: invalid q_points or omega shape'
      end if
      if (allocated(request%q_indices)) then
         if (size(request%q_indices) /= size(request%q_points, 2)) then
            error stop 'chi0 backend request: q_indices count must equal q_points count'
         end if
      end if
   end subroutine validate_request

   integer function request_q_index(request, iq) result(q_index)
      type(tddft_chi0_request), intent(in) :: request
      integer, intent(in) :: iq

      if (allocated(request%q_indices)) then
         q_index = request%q_indices(iq)
      else
         q_index = iq
      end if
   end function request_q_index

   subroutine annotate_result(result, backend, implementation, q_point, nq, nw)
      type(tddft_chi0_result), intent(inout) :: result
      character(len=*), intent(in) :: backend, implementation
      real(rp), intent(in) :: q_point(3)
      integer, intent(in) :: nq, nw

      result%metadata%backend = backend
      result%metadata%canonical_backend = backend
      result%metadata%implementation = implementation
      result%metadata%q_direct = q_point
      result%metadata%q_batch_size = nq
      result%metadata%omega_batch_size = nw
      result%metadata%convergence_status = 'not assessed by backend'
      result%metadata%converged = .false.
   end subroutine annotate_result

   subroutine annotate_batch(result, backend)
      type(tddft_chi0_batch_result), intent(inout) :: result
      character(len=*), intent(in) :: backend
      integer :: nq, nw

      if (.not. allocated(result%q_response)) error stop 'annotate_batch: q response is absent'
      nq = size(result%q_response)
      nw = 0
      if (allocated(result%omega)) nw = size(result%omega)
      result%metadata%backend = backend
      result%metadata%canonical_backend = backend
      result%metadata%q_batch_size = nq
      result%metadata%omega_batch_size = nw
      if (allocated(result%q_points)) result%metadata%q_direct = result%q_points(:, 1)
      if (nq > 0) then
         result%metadata%eta = result%q_response(1)%metadata%eta
         result%metadata%energy_integration = result%q_response(1)%metadata%energy_integration
         result%metadata%convergence_status = result%q_response(1)%metadata%convergence_status
         result%metadata%converged = result%q_response(1)%metadata%converged
         result%metadata%real_space_reuse = result%q_response(1)%metadata%real_space_reuse
      end if
   end subroutine annotate_batch

   subroutine require_batch_shape(result, nq, nw, caller)
      type(tddft_chi0_batch_result), intent(in) :: result
      integer, intent(in) :: nq, nw
      character(len=*), intent(in) :: caller

      if (.not. allocated(result%q_response) .or. .not. allocated(result%omega)) then
         error stop trim(caller)//': backend returned an incompatible batch shape'
      end if
      if (size(result%q_response) /= nq .or. size(result%omega) /= nw) then
         error stop trim(caller)//': backend returned an incompatible batch shape'
      end if
      if (.not. allocated(result%q_points)) then
         error stop trim(caller)//': backend did not return q-point provenance'
      end if
      if (size(result%q_points, 1) /= 3 .or. size(result%q_points, 2) /= nq) then
         error stop trim(caller)//': backend returned incompatible q-point provenance'
      end if
      if (allocated(result%q_indices)) then
         if (size(result%q_indices) /= nq) then
            error stop trim(caller)//': backend returned incompatible q-index provenance'
         end if
      end if
      block
         integer :: iq
         do iq = 1, nq
            if (.not. allocated(result%q_response(iq)%chi)) then
               error stop trim(caller)//': backend returned an empty chi0 response'
            end if
            if (size(result%q_response(iq)%chi, 3) /= nw) then
               error stop trim(caller)//': backend returned an incompatible frequency response'
            end if
         end do
      end block
   end subroutine require_batch_shape

   pure function lower_ascii(input) result(output)
      character(len=*), intent(in) :: input
      character(len=len(input)) :: output
      integer :: i, code

      output = input
      do i = 1, len(input)
         code = iachar(output(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) output(i:i) = achar(code + 32)
      end do
   end function lower_ascii

end module tddft_backend_mod
