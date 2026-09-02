!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Native real-space Green-function bubble for TDDFT.
!>
!> This module owns the production-independent part of TDDFT-07:
!> `G_ab(R,z), G_ba(-R,z) -> chi0_ab(R,omega) -> chi0_ab(q,omega)`.
!> The real-space response is built once for a requested frequency batch and
!> all q points are then transformed from that same tensor.  The implementation
!> deliberately never constructs a reciprocal-space Green function.
module tddft_chi0_realspace_mod
   use precision_mod, only: rp
   use math_mod, only: pi, two_pi, i_unit, inverse_3x3, gauss_legendre
   use response_vertices_mod, only: response_channel, site_projected_operator
   use tddft_chi0_mod, only: tddft_chi0_result, tddft_chi0_request, tddft_chi0_batch_result, &
      tddft_fermi_occupation, tddft_kB_Ry_per_K, tddft_occupation_kT_floor
   use tddft_backend_mod, only: tddft_realspace_chi0_provider, tddft_backend_capabilities
   use green_mod, only: green
   use lattice_mod, only: lattice
   use mpi_mod, only: ierr, start_atom, end_atom, g2l_map, numprocs
#ifdef USE_MPI
   use mpi
#endif
   implicit none

   private

   type, public :: tddft_realspace_chi0_options
      real(rp) :: eta = 0.01_rp
      real(rp) :: green_eta = 0.0_rp
      character(len=24) :: energy_integration = 'direct'
      real(rp) :: fermi_level = 0.0_rp
      real(rp) :: electronic_temperature = 0.0_rp
      real(rp) :: rmax = huge(1.0_rp)
      real(rp) :: tail_tolerance = 1.0e-3_rp
      ! `full_tail` is the validation/reference mode.  `production` skips
      ! pair response work outside rmax and therefore never infers a tail
      ! estimate from the discarded pairs.
      character(len=16) :: truncation_mode = 'full_tail'
      character(len=16) :: representation = 'bulk'
      character(len=16) :: circular_channel = 'plus_minus'
      integer :: fourier_axes(3) = [1, 2, 3]
      integer :: contour_points = 64
      integer :: contour_subdivisions = 8
      integer :: near_fermi_points = 128
      real(rp) :: contour_height = 0.0_rp
      real(rp) :: metric(3, 3) = reshape([1.0_rp, 0.0_rp, 0.0_rp, &
         0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp], [3, 3])
   end type tddft_realspace_chi0_options

   type, public :: tddft_realspace_chi0_diagnostics
      integer :: input_points = 0
      integer :: selected_points = 0
      integer :: omitted_points = 0
      integer :: shell_count = 0
      real(rp) :: rmax = 0.0_rp
      real(rp) :: requested_rmax = 0.0_rp
      real(rp) :: source_radius = 0.0_rp
      real(rp) :: selected_norm = 0.0_rp
      real(rp) :: total_norm = 0.0_rp
      real(rp) :: omitted_tail_norm = 0.0_rp
      real(rp) :: tail_ratio = 0.0_rp
      real(rp) :: last_shell_norm = 0.0_rp
      real(rp) :: pair_integration_cpu_seconds = 0.0_rp
      real(rp) :: fourier_cpu_seconds = 0.0_rp
      real(rp) :: green_function_cpu_seconds = 0.0_rp
      real(rp) :: total_cpu_seconds = 0.0_rp
      integer :: pair_response_integrations = 0
      real(rp) :: integration_energy_min = 0.0_rp
      real(rp) :: integration_energy_max = 0.0_rp
      integer :: contour_points = 0
      integer :: contour_subdivisions = 0
      integer :: near_fermi_points = 0
      real(rp) :: contour_height = 0.0_rp
      real(rp) :: contour_max_imaginary_energy = 0.0_rp
      integer :: contour_gf_evaluations = 0
      logical :: contour_deformation = .false.
      logical :: tail_assessed = .false.
      logical :: converged = .false.
      logical :: source_covers_requested_cutoff = .false.
      character(len=48) :: status = 'not assessed'
   end type tddft_realspace_chi0_diagnostics

   type, public :: tddft_realspace_chi0_result
      complex(rp), allocatable :: chi_r(:, :, :, :)
      complex(rp), allocatable :: chi_q(:, :, :, :)
      real(rp), allocatable :: r_vectors(:, :)
      type(tddft_realspace_chi0_diagnostics) :: diagnostics
   end type tddft_realspace_chi0_result

   !> Provider-side contract for native G(R,z) values.  The source owns the
   !> RS recursion and returns all locally owned pairs for one complex-energy
   !> batch.  TDDFT consumes this interface without knowing how the recursion
   !> (block or Chebyshev) constructs the Green function.
   type, abstract, public :: tddft_realspace_complex_green_source
   contains
      procedure(tddft_realspace_complex_green_batch), deferred :: get_complex
   end type tddft_realspace_complex_green_source

   abstract interface
      subroutine tddft_realspace_complex_green_batch(this, z_grid, g_ab, g_ba)
         import :: tddft_realspace_complex_green_source, rp
         class(tddft_realspace_complex_green_source), intent(inout) :: this
         complex(rp), intent(in) :: z_grid(:)
         complex(rp), allocatable, intent(out) :: g_ab(:, :, :, :), g_ba(:, :, :, :)
      end subroutine tddft_realspace_complex_green_batch
   end interface

   !> Adapter from the mature RS Green object to the TDDFT source contract.
   !> It is intentionally only a thin bridge: recursion, terminators, and
   !> pair reconstruction remain implemented by green_mod.
   type, extends(tddft_realspace_complex_green_source), public :: tddft_native_green_source
      type(green), pointer :: green_obj => null()
   contains
      procedure :: get_complex => get_native_green_complex
   end type tddft_native_green_source

   !> Concrete provider used by tddft_realspace_gf_backend.  The generic array
   !> initializer is useful for tests and for future finite/impurity sources;
   !> initialize_from_green is the bridge from the native RS-LMTO storage.
   type, extends(tddft_realspace_chi0_provider), public :: tddft_native_realspace_gf_provider
      real(rp), allocatable :: energy_grid(:)
      complex(rp), allocatable :: g_ab(:, :, :, :), g_ba(:, :, :, :)
      real(rp), allocatable :: r_vectors(:, :)
      integer, allocatable :: pair_sites(:, :)
      integer, allocatable :: site_orbital_counts(:)
      type(response_channel), allocatable :: left_channels(:), right_channels(:)
      class(tddft_realspace_complex_green_source), allocatable :: complex_source
      type(tddft_realspace_chi0_options) :: options
      integer :: build_count = 0
      ! Time spent by the caller constructing the sampled native source.  It
      ! is set by the production driver; generic array providers leave it zero.
      real(rp) :: source_green_cpu_seconds = 0.0_rp
      logical :: initialized = .false.
   contains
      procedure :: initialize => initialize_native_realspace_provider
      procedure :: initialize_from_green => initialize_native_from_green
      procedure :: initialize_complex => initialize_native_with_complex_source
      procedure :: build_realspace => provider_build_realspace
      procedure :: evaluate_realspace => evaluate_native_realspace_provider
      procedure :: evaluate_static_realspace => evaluate_native_static_realspace_provider
      procedure :: describe => describe_native_realspace_provider
   end type tddft_native_realspace_gf_provider

   public :: build_chi0_from_realspace_gf
   public :: build_static_chi0_from_realspace_gf
   public :: fourier_transform_realspace_chi0
   public :: fourier_transform_realspace_green
   public :: check_realspace_pair_reversal
   public :: reduce_realspace_chi0_batch

contains

   subroutine initialize_native_realspace_provider(this, energy_grid, g_ab, g_ba, r_vectors, pair_sites, &
      site_orbital_counts, left_channels, right_channels, options)
      class(tddft_native_realspace_gf_provider), intent(inout) :: this
      real(rp), intent(in) :: energy_grid(:), r_vectors(:, :)
      complex(rp), intent(in) :: g_ab(:, :, :, :), g_ba(:, :, :, :)
      integer, intent(in) :: pair_sites(:, :), site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_realspace_chi0_options), intent(in), optional :: options

      call validate_realspace_source(energy_grid, g_ab, g_ba, r_vectors, pair_sites, site_orbital_counts, &
         left_channels, right_channels)
      this%energy_grid = energy_grid
      this%g_ab = g_ab
      this%g_ba = g_ba
      this%r_vectors = r_vectors
      this%pair_sites = pair_sites
      this%site_orbital_counts = site_orbital_counts
      this%left_channels = left_channels
      this%right_channels = right_channels
      if (allocated(this%complex_source)) deallocate(this%complex_source)
      this%options = tddft_realspace_chi0_options()
      if (present(options)) this%options = options
      this%build_count = 0
      this%source_green_cpu_seconds = 0.0_rp
      this%initialized = .true.
   end subroutine initialize_native_realspace_provider

   subroutine initialize_native_from_green(this, green_obj, lattice_obj, site_orbital_counts, left_channels, &
      right_channels, options)
      class(tddft_native_realspace_gf_provider), intent(inout) :: this
      type(green), target, intent(in) :: green_obj
      type(lattice), intent(in) :: lattice_obj
      integer, intent(in) :: site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_realspace_chi0_options), intent(in), optional :: options
      real(rp), allocatable :: rvec(:, :)
      integer, allocatable :: pairs(:, :), local_map(:)
      complex(rp), allocatable :: local_ab(:, :, :, :), local_ba(:, :, :, :)
      type(tddft_realspace_chi0_options) :: source_options
      type(tddft_native_green_source) :: native_source
      integer :: npair, ip, pglob, ilocal, iatom, jatom

      if (.not. allocated(green_obj%gij) .or. .not. allocated(green_obj%gji) .or. &
          .not. associated(green_obj%en)) error stop 'native real-space provider: Green source is incomplete'
      npair = size(green_obj%gij, 4)
      if (size(green_obj%gji, 4) /= npair .or. size(lattice_obj%ijpair, 1) < npair) then
         error stop 'native real-space provider: Green pair storage is incompatible with lattice pairs'
      end if
      allocate(rvec(3, npair), pairs(npair, 2), local_ab(size(green_obj%gij, 1), size(green_obj%gij, 2), &
         size(green_obj%gij, 3), npair), local_ba(size(green_obj%gji, 1), size(green_obj%gji, 2), &
         size(green_obj%gji, 3), npair))
      local_ab = green_obj%gij
      local_ba = green_obj%gji
      local_map = [(ip, ip=1, npair)]
      if (allocated(g2l_map) .and. end_atom >= start_atom .and. size(g2l_map) >= end_atom) then
         do pglob = start_atom, min(end_atom, size(lattice_obj%ijpair, 1))
            ilocal = g2l_map(pglob)
            if (ilocal >= 1 .and. ilocal <= npair) local_map(ilocal) = pglob
         end do
      end if
      do ip = 1, npair
         pglob = local_map(ip)
         iatom = lattice_obj%ijpair(pglob, 1)
         jatom = lattice_obj%ijpair(pglob, 2)
         if (iatom < 1 .or. iatom > size(lattice_obj%cr, 2) .or. jatom < 1 .or. jatom > size(lattice_obj%cr, 2)) then
            error stop 'native real-space provider: lattice pair atom is outside coordinates'
         end if
         ! `jatom` is a cluster-atom index, whereas response channels use the
         ! periodic unit-cell site index.  Keep the cluster coordinate for
         ! the Fourier phase, but fold only the channel label through iz.
         if (.not. allocated(lattice_obj%iz) .or. lattice_obj%iz(jatom) < 1 .or. &
             lattice_obj%iz(jatom) > size(site_orbital_counts)) then
            error stop 'native real-space provider: cluster atom has no valid response-site mapping'
         end if
         pairs(ip, :) = [iatom, lattice_obj%iz(jatom)]
         rvec(:, ip) = realspace_pair_vector(lattice_obj, iatom, jatom)
      end do
      source_options = tddft_realspace_chi0_options()
      if (present(options)) source_options = options
      source_options%metric = lattice_obj%alat*lattice_obj%a
      call this%initialize(green_obj%en%ene(1:size(green_obj%gij, 3)), local_ab, local_ba, rvec, pairs, &
         site_orbital_counts, left_channels, right_channels, source_options)
      native_source%green_obj => green_obj
      allocate(this%complex_source, source=native_source)
   end subroutine initialize_native_from_green

   subroutine initialize_native_with_complex_source(this, energy_grid, g_ab, g_ba, r_vectors, pair_sites, &
      site_orbital_counts, left_channels, right_channels, complex_source, options)
      class(tddft_native_realspace_gf_provider), intent(inout) :: this
      real(rp), intent(in) :: energy_grid(:), r_vectors(:, :)
      complex(rp), intent(in) :: g_ab(:, :, :, :), g_ba(:, :, :, :)
      integer, intent(in) :: pair_sites(:, :), site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      class(tddft_realspace_complex_green_source), intent(in) :: complex_source
      type(tddft_realspace_chi0_options), intent(in), optional :: options

      call this%initialize(energy_grid, g_ab, g_ba, r_vectors, pair_sites, site_orbital_counts, left_channels, &
         right_channels, options)
      allocate(this%complex_source, source=complex_source)
   end subroutine initialize_native_with_complex_source

   subroutine get_native_green_complex(this, z_grid, g_ab, g_ba)
      class(tddft_native_green_source), intent(inout) :: this
      complex(rp), intent(in) :: z_grid(:)
      complex(rp), allocatable, intent(out) :: g_ab(:, :, :, :), g_ba(:, :, :, :)

      if (.not. associated(this%green_obj)) error stop 'tddft_native_green_source: Green object is not associated'
      call this%green_obj%calculate_intersite_gf_complex(z_grid, g_ab, g_ba)
   end subroutine get_native_green_complex

   subroutine provider_build_realspace(this, omega, result)
      class(tddft_native_realspace_gf_provider), intent(inout) :: this
      real(rp), intent(in) :: omega(:)
      type(tddft_realspace_chi0_result), intent(out) :: result

      if (.not. this%initialized) error stop 'native real-space provider: provider is not initialized'
      if (trim(this%options%energy_integration) == 'mixed_contour') then
         if (.not. allocated(this%complex_source)) error stop &
            'native real-space provider: mixed_contour needs an arbitrary-complex-energy source'
         call build_chi0_from_realspace_gf_mixed_contour(this%complex_source, this%energy_grid, this%r_vectors, this%pair_sites, &
            this%site_orbital_counts, this%left_channels, this%right_channels, omega, this%options, result)
      else
         call build_chi0_from_realspace_gf(this%energy_grid, this%g_ab, this%g_ba, this%r_vectors, this%pair_sites, &
            this%site_orbital_counts, this%left_channels, this%right_channels, omega, this%options, result)
      end if
      this%build_count = this%build_count + 1
   end subroutine provider_build_realspace

   subroutine evaluate_native_realspace_provider(this, request, result)
      class(tddft_native_realspace_gf_provider), intent(inout) :: this
      type(tddft_chi0_request), intent(in) :: request
      type(tddft_chi0_batch_result), intent(out) :: result
      type(tddft_realspace_chi0_result) :: realspace_result
      integer :: iq, nq, nw

      if (.not. allocated(request%q_points) .or. .not. allocated(request%omega)) then
         error stop 'native real-space provider: request is incomplete'
      end if
      nq = size(request%q_points, 2); nw = size(request%omega)
      if (size(request%q_points, 1) /= 3 .or. nq < 1 .or. nw < 1) then
         error stop 'native real-space provider: invalid request dimensions'
      end if
      if (trim(this%options%energy_integration) == 'mixed_contour') then
         if (.not. allocated(this%complex_source)) error stop &
            'native real-space provider: mixed_contour needs an arbitrary-complex-energy source'
         call build_chi0_from_realspace_gf_mixed_contour(this%complex_source, this%energy_grid, this%r_vectors, this%pair_sites, &
            this%site_orbital_counts, this%left_channels, this%right_channels, request%omega, this%options, realspace_result, &
            q_points=request%q_points)
      else
         call build_chi0_from_realspace_gf(this%energy_grid, this%g_ab, this%g_ba, this%r_vectors, this%pair_sites, &
            this%site_orbital_counts, this%left_channels, this%right_channels, request%omega, this%options, realspace_result, &
            q_points=request%q_points)
      end if
      realspace_result%diagnostics%green_function_cpu_seconds = realspace_result%diagnostics%green_function_cpu_seconds + &
         this%source_green_cpu_seconds
      realspace_result%diagnostics%total_cpu_seconds = realspace_result%diagnostics%total_cpu_seconds + &
         this%source_green_cpu_seconds
      this%build_count = this%build_count + 1
      allocate(result%q_response(nq), result%q_points(3, nq), result%q_indices(nq), result%omega(nw))
      result%q_points = request%q_points
      result%omega = request%omega
      do iq = 1, nq
         result%q_indices(iq) = request_q_index(request, iq)
         call make_response_from_q(realspace_result%chi_q(:, :, :, iq), request%q_points(:, iq), request%omega, &
            this%left_channels, this%right_channels, this%options, realspace_result%diagnostics, result%q_response(iq))
      end do
      result%metadata = result%q_response(1)%metadata
      result%metadata%q_batch_size = nq
      result%metadata%omega_batch_size = nw
      result%metadata%q_direct = request%q_points(:, 1)
      result%metadata%real_space_reuse = request%allow_real_space_reuse
      do iq = 1, nq
         result%q_response(iq)%metadata%real_space_reuse = request%allow_real_space_reuse
         result%q_response(iq)%metadata%q_batch_size = nq
      end do
   end subroutine evaluate_native_realspace_provider

   subroutine evaluate_native_static_realspace_provider(this, request, result)
      class(tddft_native_realspace_gf_provider), intent(inout) :: this
      type(tddft_chi0_request), intent(in) :: request
      type(tddft_chi0_batch_result), intent(out) :: result
      type(tddft_realspace_chi0_result) :: realspace_result
      integer :: iq, nq

      if (.not. this%initialized) error stop 'native real-space provider: provider is not initialized'
      if (.not. allocated(request%q_points) .or. .not. allocated(request%omega) .or. size(request%omega) /= 1 .or. &
          abs(request%omega(1)) > 64.0_rp*epsilon(1.0_rp)) then
         error stop 'native real-space provider: static request must contain one zero frequency'
      end if
      nq = size(request%q_points, 2)
      if (size(request%q_points, 1) /= 3 .or. nq < 1) then
         error stop 'native real-space provider: invalid static q-point dimensions'
      end if
      call build_static_chi0_from_realspace_gf(this%energy_grid, this%g_ab, this%g_ba, this%r_vectors, this%pair_sites, &
         this%site_orbital_counts, this%left_channels, this%right_channels, this%options, realspace_result, &
         q_points=request%q_points)
      realspace_result%diagnostics%green_function_cpu_seconds = realspace_result%diagnostics%green_function_cpu_seconds + &
         this%source_green_cpu_seconds
      realspace_result%diagnostics%total_cpu_seconds = realspace_result%diagnostics%total_cpu_seconds + &
         this%source_green_cpu_seconds
      this%build_count = this%build_count + 1
      allocate(result%q_response(nq), result%q_points(3, nq), result%q_indices(nq), result%omega(1))
      result%q_points = request%q_points
      result%omega(1) = 0.0_rp
      do iq = 1, nq
         result%q_indices(iq) = request_q_index(request, iq)
         call make_response_from_q(realspace_result%chi_q(:, :, :, iq), request%q_points(:, iq), [0.0_rp], &
            this%left_channels, this%right_channels, this%options, realspace_result%diagnostics, result%q_response(iq))
         call mark_static_realspace_metadata(result%q_response(iq), request%q_points(:, iq), nq)
      end do
      result%metadata = result%q_response(1)%metadata
      result%metadata%q_batch_size = nq
      result%metadata%omega_batch_size = 1
   end subroutine evaluate_native_static_realspace_provider

   subroutine describe_native_realspace_provider(this, metadata)
      class(tddft_native_realspace_gf_provider), intent(in) :: this
      type(tddft_backend_capabilities), intent(out) :: metadata

      metadata = tddft_backend_capabilities()
      metadata%native_real_space = .true.
      metadata%reuses_real_space_response = .true.
      metadata%supports_static_limit = .true.
   end subroutine describe_native_realspace_provider

   subroutine build_chi0_from_realspace_gf(energy_grid, g_ab, g_ba, r_vectors, pair_sites, site_orbital_counts, &
      left_channels, right_channels, omega, options, result, q_points)
      real(rp), intent(in) :: energy_grid(:), r_vectors(:, :), omega(:)
      real(rp), intent(in), optional :: q_points(:, :)
      complex(rp), intent(in) :: g_ab(:, :, :, :), g_ba(:, :, :, :)
      integer, intent(in) :: pair_sites(:, :), site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_realspace_chi0_options), intent(in) :: options
      type(tddft_realspace_chi0_result), intent(out) :: result
      logical, allocatable :: keep(:)
      real(rp), allocatable :: radii(:)
      complex(rp), allocatable :: left_ops(:, :, :), right_ops(:, :, :), pair_response(:, :, :)
      real(rp), allocatable :: pair_norms(:)
      complex(rp), allocatable :: f0(:, :), f1(:, :), gr0(:, :), gr1(:, :), grm(:, :), rev0(:, :), &
         adv_reverse0(:, :), adv_forward0(:, :), advm(:, :), work1(:, :), work2(:, :), work3(:, :)
      integer :: ip, iw, ikeep, nkeep, nleft, nright, nw, nblock
      real(rp) :: effective_cutoff, source_radius, selected_norm, omitted_norm, total_norm
      real(rp) :: t_start, t_stop, total_start, total_stop
      logical :: full_tail_mode

      call validate_realspace_source(energy_grid, g_ab, g_ba, r_vectors, pair_sites, site_orbital_counts, &
         left_channels, right_channels)
      if (size(omega) < 1 .or. options%eta <= 0.0_rp .or. options%tail_tolerance < 0.0_rp) then
         error stop 'build_chi0_from_realspace_gf: invalid frequency or broadening options'
      end if
      if (trim(options%energy_integration) /= 'direct') then
         error stop 'build_chi0_from_realspace_gf: mixed_contour requires a complex-energy native real-space source; sampled real-axis source is direct-only'
      end if
      full_tail_mode = realspace_full_tail_mode(options%truncation_mode)
      nleft = size(left_channels); nright = size(right_channels); nw = size(omega)
      nblock = size(g_ab, 1)
      call cpu_time(total_start)
      allocate(radii(size(r_vectors, 2)), keep(size(r_vectors, 2)))
      do ip = 1, size(r_vectors, 2)
         radii(ip) = metric_norm(options%metric, r_vectors(:, ip))
         keep(ip) = radii(ip) <= options%rmax + 16.0_rp*epsilon(1.0_rp)*max(1.0_rp, radii(ip))
      end do
      nkeep = count(keep)
      if (nkeep < 1 .and. size(r_vectors, 2) > 0 .and. realspace_parallel_size() == 1) then
         error stop 'build_chi0_from_realspace_gf: Rmax removes every real-space pair'
      end if
      if (nkeep > 0) then
         effective_cutoff = maxval(pack(radii, keep))
      else
         ! A rank with no locally owned R block is a valid member of an
         ! MPI-distributed native source.  Its zero contribution is reduced
         ! with the other ranks before the common q Fourier batch is used.
         effective_cutoff = 0.0_rp
      end if
      if (size(radii) > 0) then
         source_radius = maxval(radii)
      else
         source_radius = 0.0_rp
      end if
      allocate(result%chi_r(nleft, nright, nw, nkeep), result%r_vectors(3, nkeep))
      result%chi_r = cmplx(0.0_rp, 0.0_rp, rp)
      ikeep = 0
      do ip = 1, size(keep)
         if (keep(ip)) then
            ikeep = ikeep + 1
            result%r_vectors(:, ikeep) = r_vectors(:, ip)
         end if
      end do
      allocate(left_ops(nblock, nblock, nleft), right_ops(nblock, nblock, nright))
      call make_local_operators(left_channels, site_orbital_counts, nblock, left_ops)
      call make_local_operators(right_channels, site_orbital_counts, nblock, right_ops)
      allocate(pair_response(nleft, nright, nw), pair_norms(size(keep)))
      allocate(f0(nleft, nright), f1(nleft, nright), gr0(nblock, nblock), gr1(nblock, nblock), &
         grm(nblock, nblock), rev0(nblock, nblock), adv_reverse0(nblock, nblock), &
         adv_forward0(nblock, nblock), advm(nblock, nblock), work1(nblock, nblock), &
         work2(nblock, nblock), work3(nblock, nblock))
      pair_response = cmplx(0.0_rp, 0.0_rp, rp)
      selected_norm = 0.0_rp; omitted_norm = 0.0_rp; total_norm = 0.0_rp; pair_norms = 0.0_rp
      ikeep = 0
      call cpu_time(t_start)
      do ip = 1, size(keep)
         if (.not. keep(ip) .and. .not. full_tail_mode) cycle
         pair_response = cmplx(0.0_rp, 0.0_rp, rp)
         do iw = 1, nw
            call integrate_pair_response(energy_grid, g_ab(:, :, :, ip), g_ba(:, :, :, ip), pair_sites(ip, :), &
               left_channels, right_channels, left_ops, right_ops, omega(iw), options, f0, f1, gr0, gr1, grm, rev0, &
               adv_reverse0, adv_forward0, advm, work1, work2, work3, pair_response(:, :, iw))
         end do
         pair_norms(ip) = sum(abs(pair_response))
         total_norm = total_norm + pair_norms(ip)
         if (keep(ip)) then
            ikeep = ikeep + 1
            result%chi_r(:, :, :, ikeep) = pair_response
            selected_norm = selected_norm + pair_norms(ip)
         else
            omitted_norm = omitted_norm + pair_norms(ip)
         end if
      end do
      call cpu_time(t_stop)
      result%diagnostics%pair_integration_cpu_seconds = t_stop-t_start
      call cpu_time(t_start)
      if (present(q_points)) then
         call fourier_transform_realspace_chi0(result%chi_r, result%r_vectors, q_points, options%representation, &
            options%fourier_axes, result%chi_q)
      else
         allocate(result%chi_q(nleft, nright, nw, 1))
         result%chi_q = cmplx(0.0_rp, 0.0_rp, rp)
      end if
      call cpu_time(t_stop)
      result%diagnostics%fourier_cpu_seconds = t_stop-t_start
      result%diagnostics%input_points = size(keep)
      result%diagnostics%integration_energy_min = energy_grid(1)
      result%diagnostics%integration_energy_max = energy_grid(size(energy_grid))
      result%diagnostics%selected_points = nkeep
      result%diagnostics%omitted_points = size(keep)-nkeep
      result%diagnostics%rmax = effective_cutoff
      result%diagnostics%requested_rmax = options%rmax
      result%diagnostics%source_radius = source_radius
      result%diagnostics%source_covers_requested_cutoff = source_radius + &
         16.0_rp*epsilon(1.0_rp)*max(1.0_rp, source_radius, abs(options%rmax)) >= options%rmax
      result%diagnostics%selected_norm = selected_norm
      result%diagnostics%omitted_tail_norm = omitted_norm
      result%diagnostics%total_norm = total_norm
      if (full_tail_mode .and. total_norm > tiny(1.0_rp)) result%diagnostics%tail_ratio = omitted_norm/total_norm
      result%diagnostics%tail_assessed = full_tail_mode .and. result%diagnostics%omitted_points > 0
      result%diagnostics%converged = result%diagnostics%tail_assessed .and. &
         result%diagnostics%tail_ratio <= options%tail_tolerance
      result%diagnostics%pair_response_integrations = nkeep*nw
      if (full_tail_mode) result%diagnostics%pair_response_integrations = size(keep)*nw
      if (.not. full_tail_mode .and. result%diagnostics%omitted_points > 0) then
         result%diagnostics%status = 'production truncation; omitted tail not assessed'
      else if (.not. result%diagnostics%source_covers_requested_cutoff) then
         result%diagnostics%status = 'source radius does not extend beyond requested R_max; tail not assessed'
      else if (.not. result%diagnostics%tail_assessed) then
         result%diagnostics%status = 'all supplied real-space pairs retained'
      else if (result%diagnostics%converged) then
         result%diagnostics%status = 'R-space tail below tolerance'
      else
         result%diagnostics%status = 'R-space tail exceeds tolerance'
      end if
      result%diagnostics%last_shell_norm = last_shell_norm(radii, keep, pair_norms)
      result%diagnostics%shell_count = count_realspace_shells(radii)
      call cpu_time(total_stop)
      result%diagnostics%total_cpu_seconds = total_stop-total_start
   end subroutine build_chi0_from_realspace_gf

   !> Native mixed-contour counterpart of build_chi0_from_realspace_gf.
   !>
   !> The two same-half-plane pieces are evaluated on upper/lower polylines
   !> and the mixed retarded/advanced piece remains on its finite real-energy
   !> interval.  Every contour node requests both directed pair Green
   !> functions from the source in one batch.  The resulting chi0(R,omega) is
   !> assembled before the one Fourier transform to all requested q points.
   subroutine build_chi0_from_realspace_gf_mixed_contour(source, energy_grid, r_vectors, pair_sites, site_orbital_counts, &
      left_channels, right_channels, omega, options, result, q_points)
      class(tddft_realspace_complex_green_source), intent(inout) :: source
      real(rp), intent(in) :: energy_grid(:), r_vectors(:, :), omega(:)
      real(rp), intent(in), optional :: q_points(:, :)
      integer, intent(in) :: pair_sites(:, :), site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_realspace_chi0_options), intent(in) :: options
      type(tddft_realspace_chi0_result), intent(out) :: result
      logical, allocatable :: keep(:), active(:)
      real(rp), allocatable :: radii(:), pair_norms(:)
      complex(rp), allocatable :: left_ops(:, :, :), right_ops(:, :, :), pair_response(:, :, :, :)
      complex(rp), allocatable :: rr_integral(:, :, :), ra_integral(:, :, :), aa_integral(:, :, :)
      integer :: ip, iw, iw_prev, ileft, iright, ikeep, nkeep, nleft, nright, nw, nblock
      real(rp) :: energy_min, energy_max, contour_height, same_min, same_max
      real(rp) :: thermal_window, effective_cutoff, source_radius, selected_norm, omitted_norm, total_norm
      real(rp) :: t_pair_start, t_pair_stop, t_fourier_start, t_fourier_stop, total_start, total_stop
      real(rp) :: green_function_cpu_seconds
      integer :: contour_gf_evaluations
      real(rp) :: max_contour_imaginary_energy
      logical :: repeated_frequency, full_tail_mode

      call validate_realspace_geometry(r_vectors, pair_sites, site_orbital_counts, left_channels, right_channels)
      if (size(omega) < 1 .or. options%eta <= 0.0_rp .or. options%tail_tolerance < 0.0_rp .or. &
          options%electronic_temperature < 0.0_rp) then
         error stop 'build_chi0_from_realspace_gf_mixed_contour: invalid options'
      end if
      if (size(energy_grid) < 2 .or. any(energy_grid(2:) <= energy_grid(:size(energy_grid)-1))) then
         error stop 'build_chi0_from_realspace_gf_mixed_contour: invalid real-axis reference grid'
      end if
      full_tail_mode = realspace_full_tail_mode(options%truncation_mode)
      nleft = size(left_channels); nright = size(right_channels); nw = size(omega)
      nblock = 2*site_orbital_counts(1)
      call cpu_time(total_start)
      ! The real-axis grid is retained by the native provider as the reference
      ! window.  The contour source itself is independent of that grid, but
      ! the finite-window identity must use exactly the same endpoints.
      energy_min = energy_grid(1)
      energy_max = energy_grid(size(energy_grid))
      if (energy_max <= energy_min) error stop 'build_chi0_from_realspace_gf_mixed_contour: invalid energy window'
      call mixed_native_contour_height(options, energy_max-energy_min, resolved_green_eta(options), contour_height)
      call mixed_native_same_interval(options, energy_min, energy_max, same_min, same_max)
      thermal_window = 50.0_rp*max(options%electronic_temperature*tddft_kB_Ry_per_K, tddft_occupation_kT_floor)

      allocate(radii(size(r_vectors, 2)), keep(size(r_vectors, 2)))
      do ip = 1, size(keep)
         radii(ip) = metric_norm(options%metric, r_vectors(:, ip))
         keep(ip) = radii(ip) <= options%rmax + 16.0_rp*epsilon(1.0_rp)*max(1.0_rp, radii(ip))
      end do
      allocate(active(size(keep)))
      active = keep
      if (full_tail_mode) active = .true.
      nkeep = count(keep)
      if (nkeep < 1 .and. size(r_vectors, 2) > 0 .and. realspace_parallel_size() == 1) then
         error stop 'build_chi0_from_realspace_gf_mixed_contour: Rmax removes every real-space pair'
      end if
      if (nkeep > 0) then
         effective_cutoff = maxval(pack(radii, keep))
      else
         ! A rank with no selected local pairs contributes a correctly shaped
         ! zero tensor; the production driver reduces it with other ranks.
         effective_cutoff = 0.0_rp
      end if
      if (size(radii) > 0) then
         source_radius = maxval(radii)
      else
         source_radius = 0.0_rp
      end if
      allocate(result%chi_r(nleft, nright, nw, nkeep), result%r_vectors(3, nkeep))
      result%chi_r = cmplx(0.0_rp, 0.0_rp, rp)
      ikeep = 0
      do ip = 1, size(keep)
         if (keep(ip)) then
            ikeep = ikeep+1
            result%r_vectors(:, ikeep) = r_vectors(:, ip)
         end if
      end do
      allocate(left_ops(nblock, nblock, nleft), right_ops(nblock, nblock, nright))
      call make_local_operators(left_channels, site_orbital_counts, nblock, left_ops)
      call make_local_operators(right_channels, site_orbital_counts, nblock, right_ops)
      allocate(pair_response(nleft, nright, nw, size(keep)), pair_norms(size(keep)), &
         rr_integral(nleft, nright, size(keep)), ra_integral(nleft, nright, size(keep)), &
         aa_integral(nleft, nright, size(keep)))
      pair_response = cmplx(0.0_rp, 0.0_rp, rp)
      contour_gf_evaluations = 0; max_contour_imaginary_energy = resolved_green_eta(options)
      green_function_cpu_seconds = 0.0_rp
      call cpu_time(t_pair_start)
      do iw = 1, nw
         repeated_frequency = .false.
         do iw_prev = 1, iw-1
            if (abs(omega(iw)-omega(iw_prev)) <= 32.0_rp*epsilon(1.0_rp)* &
               max(1.0_rp, abs(omega(iw)), abs(omega(iw_prev)))) then
               pair_response(:, :, iw, :) = pair_response(:, :, iw_prev, :)
               repeated_frequency = .true.
               exit
            end if
         end do
         if (repeated_frequency) cycle
         rr_integral = cmplx(0.0_rp, 0.0_rp, rp)
         ra_integral = cmplx(0.0_rp, 0.0_rp, rp)
         aa_integral = cmplx(0.0_rp, 0.0_rp, rp)
         if (same_max > same_min) then
            call integrate_native_same_contour_segment(source, pair_sites, active, left_channels, right_channels, left_ops, right_ops, &
               omega(iw), resolved_green_eta(options), options, 1, cmplx(same_min, 0.0_rp, rp), &
               cmplx(same_min, contour_height, rp), rr_integral, aa_integral, contour_gf_evaluations, &
               max_contour_imaginary_energy, green_function_cpu_seconds)
            call integrate_native_same_contour_segment(source, pair_sites, active, left_channels, right_channels, left_ops, right_ops, &
               omega(iw), resolved_green_eta(options), options, 1, cmplx(same_min, contour_height, rp), &
               cmplx(same_max, contour_height, rp), rr_integral, aa_integral, contour_gf_evaluations, &
               max_contour_imaginary_energy, green_function_cpu_seconds)
            call integrate_native_same_contour_segment(source, pair_sites, active, left_channels, right_channels, left_ops, right_ops, &
               omega(iw), resolved_green_eta(options), options, 1, cmplx(same_max, contour_height, rp), &
               cmplx(same_max, 0.0_rp, rp), rr_integral, aa_integral, contour_gf_evaluations, &
               max_contour_imaginary_energy, green_function_cpu_seconds)
            call integrate_native_same_contour_segment(source, pair_sites, active, left_channels, right_channels, left_ops, right_ops, &
               omega(iw), resolved_green_eta(options), options, -1, cmplx(same_min, 0.0_rp, rp), &
               cmplx(same_min, -contour_height, rp), rr_integral, aa_integral, contour_gf_evaluations, &
               max_contour_imaginary_energy, green_function_cpu_seconds)
            call integrate_native_same_contour_segment(source, pair_sites, active, left_channels, right_channels, left_ops, right_ops, &
               omega(iw), resolved_green_eta(options), options, -1, cmplx(same_min, -contour_height, rp), &
               cmplx(same_max, -contour_height, rp), rr_integral, aa_integral, contour_gf_evaluations, &
               max_contour_imaginary_energy, green_function_cpu_seconds)
            call integrate_native_same_contour_segment(source, pair_sites, active, left_channels, right_channels, left_ops, right_ops, &
               omega(iw), resolved_green_eta(options), options, -1, cmplx(same_max, -contour_height, rp), &
               cmplx(same_max, 0.0_rp, rp), rr_integral, aa_integral, contour_gf_evaluations, &
               max_contour_imaginary_energy, green_function_cpu_seconds)
         end if
         call integrate_native_cross_interval(source, pair_sites, active, left_channels, right_channels, left_ops, right_ops, omega(iw), &
            resolved_green_eta(options), options, energy_min, energy_max, thermal_window, ra_integral, contour_gf_evaluations, &
            max_contour_imaginary_energy, green_function_cpu_seconds)
         do ip = 1, size(keep)
            if (keep(ip)) pair_response(:, :, iw, ip) = -(rr_integral(:, :, ip)-ra_integral(:, :, ip)-aa_integral(:, :, ip))/ &
               (2.0_rp*pi*i_unit)
         end do
      end do
      call cpu_time(t_pair_stop)

      selected_norm = 0.0_rp; omitted_norm = 0.0_rp; total_norm = 0.0_rp; pair_norms = 0.0_rp
      ikeep = 0
      do ip = 1, size(keep)
         pair_norms(ip) = sum(abs(pair_response(:, :, :, ip)))
         total_norm = total_norm+pair_norms(ip)
         if (keep(ip)) then
            ikeep = ikeep+1
            result%chi_r(:, :, :, ikeep) = pair_response(:, :, :, ip)
            selected_norm = selected_norm+pair_norms(ip)
         else
            omitted_norm = omitted_norm+pair_norms(ip)
         end if
      end do
      call cpu_time(t_fourier_start)
      if (present(q_points)) then
         call fourier_transform_realspace_chi0(result%chi_r, result%r_vectors, q_points, options%representation, &
            options%fourier_axes, result%chi_q)
      else
         allocate(result%chi_q(nleft, nright, nw, 1)); result%chi_q = cmplx(0.0_rp, 0.0_rp, rp)
      end if
      call cpu_time(t_fourier_stop)
      result%diagnostics%input_points = size(keep)
      result%diagnostics%selected_points = nkeep
      result%diagnostics%omitted_points = size(keep)-nkeep
      result%diagnostics%shell_count = count_realspace_shells(radii)
      result%diagnostics%rmax = effective_cutoff
      result%diagnostics%requested_rmax = options%rmax
      result%diagnostics%source_radius = source_radius
      result%diagnostics%source_covers_requested_cutoff = source_radius + &
         16.0_rp*epsilon(1.0_rp)*max(1.0_rp, source_radius, abs(options%rmax)) >= options%rmax
      result%diagnostics%integration_energy_min = energy_min
      result%diagnostics%integration_energy_max = energy_max
      result%diagnostics%selected_norm = selected_norm
      result%diagnostics%omitted_tail_norm = omitted_norm
      result%diagnostics%total_norm = total_norm
      if (full_tail_mode .and. total_norm > tiny(1.0_rp)) result%diagnostics%tail_ratio = omitted_norm/total_norm
      result%diagnostics%tail_assessed = full_tail_mode .and. result%diagnostics%omitted_points > 0
      result%diagnostics%converged = result%diagnostics%tail_assessed .and. result%diagnostics%tail_ratio <= options%tail_tolerance
      result%diagnostics%pair_response_integrations = nkeep*count_unique_frequencies(omega)
      if (full_tail_mode) result%diagnostics%pair_response_integrations = size(keep)*count_unique_frequencies(omega)
      if (.not. full_tail_mode .and. result%diagnostics%omitted_points > 0) then
         result%diagnostics%status = 'mixed contour; production truncation; omitted tail not assessed'
      else if (.not. result%diagnostics%source_covers_requested_cutoff) then
         result%diagnostics%status = 'mixed contour; source radius does not extend beyond requested R_max; tail not assessed'
      else if (.not. result%diagnostics%tail_assessed) then
         result%diagnostics%status = 'mixed contour; all real-space pairs retained'
      else if (result%diagnostics%converged) then
         result%diagnostics%status = 'mixed contour; R-space tail below tolerance'
      else
         result%diagnostics%status = 'mixed contour; R-space tail exceeds tolerance'
      end if
      result%diagnostics%pair_integration_cpu_seconds = t_pair_stop-t_pair_start
      result%diagnostics%fourier_cpu_seconds = t_fourier_stop-t_fourier_start
      result%diagnostics%green_function_cpu_seconds = green_function_cpu_seconds
      result%diagnostics%contour_points = options%contour_points
      result%diagnostics%contour_subdivisions = options%contour_subdivisions
      result%diagnostics%near_fermi_points = options%near_fermi_points
      result%diagnostics%contour_height = contour_height
      result%diagnostics%contour_max_imaginary_energy = max_contour_imaginary_energy
      result%diagnostics%contour_gf_evaluations = contour_gf_evaluations
      result%diagnostics%contour_deformation = .true.
      call cpu_time(total_stop)
      result%diagnostics%total_cpu_seconds = total_stop-total_start
   end subroutine build_chi0_from_realspace_gf_mixed_contour

   subroutine integrate_native_same_contour_segment(source, pair_sites, active, left_channels, right_channels, left_ops, right_ops, &
      omega, green_eta, options, plane, z_start, z_end, rr_integral, aa_integral, gf_evaluations, max_imag, gf_seconds)
      class(tddft_realspace_complex_green_source), intent(inout) :: source
      integer, intent(in) :: pair_sites(:, :), plane
      logical, intent(in) :: active(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      complex(rp), intent(in) :: left_ops(:, :, :), right_ops(:, :, :)
      real(rp), intent(in) :: omega, green_eta
      type(tddft_realspace_chi0_options), intent(in) :: options
      complex(rp), intent(in) :: z_start, z_end
      complex(rp), intent(inout) :: rr_integral(:, :, :), aa_integral(:, :, :)
      integer, intent(inout) :: gf_evaluations
      real(rp), intent(inout) :: max_imag
      real(rp), intent(inout) :: gf_seconds
      real(rp) :: nodes(options%contour_points), weights(options%contour_points), parameter
      complex(rp) :: segment_start, segment_end, z, fermi, jacobian
      complex(rp), allocatable :: z_grid(:), g_ab_batch(:, :, :, :), g_ba_batch(:, :, :, :)
      complex(rp), allocatable :: work1(:, :), work2(:, :), work3(:, :)
      complex(rp) :: value
      real(rp) :: gf_start, gf_stop
      integer :: nsub, isub, inode, ip, il, ir

      nsub = 1
      if (abs(real(z_end-z_start, rp)) > abs(aimag(z_end-z_start))) nsub = options%contour_subdivisions
      allocate(work1(size(left_ops, 1), size(left_ops, 1)), work2(size(left_ops, 1), size(left_ops, 1)), &
         work3(size(left_ops, 1), size(left_ops, 1)))
      do isub = 1, nsub
         segment_start = z_start + (z_end-z_start)*real(isub-1, rp)/real(nsub, rp)
         segment_end = z_start + (z_end-z_start)*real(isub, rp)/real(nsub, rp)
         call gauss_legendre(options%contour_points, 0.0_rp, 1.0_rp, nodes, weights)
         allocate(z_grid(2*options%contour_points))
         do inode = 1, options%contour_points
            parameter = nodes(inode)
            z = segment_start + parameter*(segment_end-segment_start)
            if (plane > 0) then
               z_grid(2*inode-1) = z + cmplx(omega, green_eta, rp)
               z_grid(2*inode) = z + cmplx(0.0_rp, green_eta, rp)
            else
               z_grid(2*inode-1) = z - cmplx(0.0_rp, green_eta, rp)
               z_grid(2*inode) = z - cmplx(omega, green_eta, rp)
            end if
         end do
         call cpu_time(gf_start)
         call source%get_complex(z_grid, g_ab_batch, g_ba_batch)
         call cpu_time(gf_stop)
         gf_seconds = gf_seconds + gf_stop-gf_start
         if (size(g_ab_batch, 1) /= size(left_ops, 1) .or. size(g_ab_batch, 2) /= size(left_ops, 2) .or. &
             size(g_ba_batch, 1) /= size(left_ops, 1) .or. size(g_ba_batch, 2) /= size(left_ops, 2) .or. &
             size(g_ab_batch, 3) /= size(z_grid) .or. size(g_ba_batch, 3) /= size(z_grid) .or. &
             size(g_ab_batch, 4) /= size(pair_sites, 1) .or. size(g_ba_batch, 4) /= size(pair_sites, 1)) then
            error stop 'native mixed contour: complex source returned incompatible batch dimensions'
         end if
         gf_evaluations = gf_evaluations + 2*size(pair_sites, 1)*size(z_grid)
         max_imag = max(max_imag, maxval(abs(aimag(z_grid))))
         do inode = 1, options%contour_points
            z = segment_start + nodes(inode)*(segment_end-segment_start)
            if (abs(real(segment_end-segment_start, rp)) > abs(aimag(segment_end-segment_start))) then
               jacobian = weights(inode)*real(segment_end-segment_start, rp)
            else
               jacobian = weights(inode)*aimag(segment_end-segment_start)*i_unit
            end if
            if (options%electronic_temperature > 0.0_rp) then
               fermi = complex_fermi_occupation(z, options%fermi_level, options%electronic_temperature)
            else
               fermi = cmplx(1.0_rp, 0.0_rp, rp)
            end if
            do ip = 1, size(pair_sites, 1)
               if (.not. active(ip)) cycle
               do il = 1, size(left_channels)
                  if (left_channels(il)%site /= pair_sites(ip, 1)) cycle
                  do ir = 1, size(right_channels)
                     if (right_channels(ir)%site /= pair_sites(ip, 2)) cycle
                     value = trace_four_local(left_ops(:, :, il), g_ab_batch(:, :, 2*inode-1, ip), &
                        right_ops(:, :, ir), g_ba_batch(:, :, 2*inode, ip), &
                        work1, work2, work3)
                     if (plane > 0) then
                        rr_integral(il, ir, ip) = rr_integral(il, ir, ip) + jacobian*fermi*value
                     else
                        aa_integral(il, ir, ip) = aa_integral(il, ir, ip) + jacobian*fermi*value
                     end if
                  end do
               end do
            end do
         end do
         deallocate(z_grid, g_ab_batch, g_ba_batch)
      end do
   end subroutine integrate_native_same_contour_segment

   subroutine integrate_native_cross_interval(source, pair_sites, active, left_channels, right_channels, left_ops, right_ops, omega, &
      green_eta, options, energy_min, energy_max, thermal_window, ra_integral, gf_evaluations, max_imag, gf_seconds)
      class(tddft_realspace_complex_green_source), intent(inout) :: source
      integer, intent(in) :: pair_sites(:, :)
      logical, intent(in) :: active(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      complex(rp), intent(in) :: left_ops(:, :, :), right_ops(:, :, :)
      real(rp), intent(in) :: omega, green_eta, energy_min, energy_max, thermal_window
      type(tddft_realspace_chi0_options), intent(in) :: options
      complex(rp), intent(inout) :: ra_integral(:, :, :)
      integer, intent(inout) :: gf_evaluations
      real(rp), intent(inout) :: max_imag
      real(rp), intent(inout) :: gf_seconds
      real(rp) :: lower, upper

      if (abs(omega) <= tiny(1.0_rp)) return
      if (omega > 0.0_rp) then
         lower = max(energy_min, min(options%fermi_level, options%fermi_level-omega)-thermal_window)
         upper = min(energy_max-omega, max(options%fermi_level, options%fermi_level-omega)+thermal_window)
         call integrate_native_cross_subinterval(source, pair_sites, active, left_channels, right_channels, left_ops, right_ops, omega, &
            green_eta, options, energy_min-omega, energy_min, -1, ra_integral, gf_evaluations, max_imag, gf_seconds)
         call integrate_native_cross_subinterval(source, pair_sites, active, left_channels, right_channels, left_ops, right_ops, omega, &
            green_eta, options, lower, upper, 0, ra_integral, gf_evaluations, max_imag, gf_seconds)
         call integrate_native_cross_subinterval(source, pair_sites, active, left_channels, right_channels, left_ops, right_ops, omega, &
            green_eta, options, energy_max-omega, energy_max, 1, ra_integral, gf_evaluations, max_imag, gf_seconds)
      else
         lower = max(energy_min-omega, min(options%fermi_level, options%fermi_level-omega)-thermal_window)
         upper = min(energy_max, max(options%fermi_level, options%fermi_level-omega)+thermal_window)
         call integrate_native_cross_subinterval(source, pair_sites, active, left_channels, right_channels, left_ops, right_ops, omega, &
            green_eta, options, energy_min, energy_min-omega, 1, ra_integral, gf_evaluations, max_imag, gf_seconds)
         call integrate_native_cross_subinterval(source, pair_sites, active, left_channels, right_channels, left_ops, right_ops, omega, &
            green_eta, options, lower, upper, 0, ra_integral, gf_evaluations, max_imag, gf_seconds)
         call integrate_native_cross_subinterval(source, pair_sites, active, left_channels, right_channels, left_ops, right_ops, omega, &
            green_eta, options, energy_max, energy_max-omega, -1, ra_integral, gf_evaluations, max_imag, gf_seconds)
      end if
   end subroutine integrate_native_cross_interval

   subroutine integrate_native_cross_subinterval(source, pair_sites, active, left_channels, right_channels, left_ops, right_ops, omega, &
      green_eta, options, lower, upper, coefficient, ra_integral, gf_evaluations, max_imag, gf_seconds)
      class(tddft_realspace_complex_green_source), intent(inout) :: source
      integer, intent(in) :: pair_sites(:, :), coefficient
      logical, intent(in) :: active(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      complex(rp), intent(in) :: left_ops(:, :, :), right_ops(:, :, :)
      real(rp), intent(in) :: omega, green_eta, lower, upper
      type(tddft_realspace_chi0_options), intent(in) :: options
      complex(rp), intent(inout) :: ra_integral(:, :, :)
      integer, intent(inout) :: gf_evaluations
      real(rp), intent(inout) :: max_imag
      real(rp), intent(inout) :: gf_seconds
      real(rp) :: nodes(options%near_fermi_points), weights(options%near_fermi_points)
      real(rp) :: energy, occupation_difference
      complex(rp), allocatable :: z_grid(:), g_ab_batch(:, :, :, :), g_ba_batch(:, :, :, :)
      complex(rp), allocatable :: work1(:, :), work2(:, :), work3(:, :)
      complex(rp) :: value
      integer :: inode, ip, il, ir
      real(rp) :: gf_start, gf_stop

      if (upper <= lower) return
      call gauss_legendre(options%near_fermi_points, lower, upper, nodes, weights)
      allocate(work1(size(left_ops, 1), size(left_ops, 1)), work2(size(left_ops, 1), size(left_ops, 1)), &
         work3(size(left_ops, 1), size(left_ops, 1)))
      allocate(z_grid(2*options%near_fermi_points))
      do inode = 1, options%near_fermi_points
         energy = nodes(inode)
         z_grid(2*inode-1) = cmplx(energy+omega, green_eta, rp)
         z_grid(2*inode) = cmplx(energy, -green_eta, rp)
      end do
      call cpu_time(gf_start)
      call source%get_complex(z_grid, g_ab_batch, g_ba_batch)
      call cpu_time(gf_stop)
      gf_seconds = gf_seconds + gf_stop-gf_start
      if (size(g_ab_batch, 1) /= size(left_ops, 1) .or. size(g_ab_batch, 2) /= size(left_ops, 2) .or. &
          size(g_ba_batch, 1) /= size(left_ops, 1) .or. size(g_ba_batch, 2) /= size(left_ops, 2) .or. &
          size(g_ab_batch, 3) /= size(z_grid) .or. size(g_ba_batch, 3) /= size(z_grid) .or. &
          size(g_ab_batch, 4) /= size(pair_sites, 1) .or. size(g_ba_batch, 4) /= size(pair_sites, 1)) then
         error stop 'native mixed contour: complex source returned incompatible cross batch dimensions'
      end if
      gf_evaluations = gf_evaluations + 2*size(pair_sites, 1)*size(z_grid)
      max_imag = max(max_imag, maxval(abs(aimag(z_grid))))
      do inode = 1, options%near_fermi_points
         energy = nodes(inode)
         select case (coefficient)
         case (-1)
            occupation_difference = -tddft_fermi_occupation(energy+omega, options%fermi_level, options%electronic_temperature)
         case (0)
            occupation_difference = tddft_fermi_occupation(energy, options%fermi_level, options%electronic_temperature)- &
               tddft_fermi_occupation(energy+omega, options%fermi_level, options%electronic_temperature)
         case (1)
            occupation_difference = tddft_fermi_occupation(energy, options%fermi_level, options%electronic_temperature)
         case default
            error stop 'native mixed contour: unknown finite-window coefficient'
         end select
         if (abs(occupation_difference) <= tiny(1.0_rp)) cycle
         do ip = 1, size(pair_sites, 1)
            if (.not. active(ip)) cycle
            do il = 1, size(left_channels)
               if (left_channels(il)%site /= pair_sites(ip, 1)) cycle
               do ir = 1, size(right_channels)
                  if (right_channels(ir)%site /= pair_sites(ip, 2)) cycle
                  value = trace_four_local(left_ops(:, :, il), g_ab_batch(:, :, 2*inode-1, ip), &
                     right_ops(:, :, ir), g_ba_batch(:, :, 2*inode, ip), &
                     work1, work2, work3)
                  ra_integral(il, ir, ip) = ra_integral(il, ir, ip) + weights(inode)*occupation_difference*value
               end do
            end do
         end do
      end do
      deallocate(z_grid, g_ab_batch, g_ba_batch)
   end subroutine integrate_native_cross_subinterval

   !> Build the true static native-RGF response from the retarded real-space
   !> Green functions.  This is not the dynamic routine sampled at omega=0.
   !>
   !> For a directed pair (a,b,R), define
   !>
   !>   T^R(E) = Tr[A G^R_ab(R,E) B G^R_ba(-R,E)]
   !>   T^A(E) = Tr[A G^A_ab(R,E) B G^A_ba(-R,E)].
   !>
   !> The static contour identity used here is
   !>
   !>   chi^0_ab(R,0) = -1/(2*pi*i) integral dE f(E) [T^R(E)-T^A(E)].
   !>
   !> Inserting the spectral representation gives the controlled divided
   !> difference (f_n-f_m)/(e_n-e_m); the equal-energy limit is the Fermi derivative.
   !> Therefore no response eta, finite-eta dynamic surrogate, or empirical
   !> normalization is part of this definition. The native `green%gij/gji`
   !> source supplies the retarded real-axis representation used here.
   subroutine build_static_chi0_from_realspace_gf(energy_grid, g_ab, g_ba, r_vectors, pair_sites, site_orbital_counts, &
      left_channels, right_channels, options, result, q_points)
      real(rp), intent(in) :: energy_grid(:), r_vectors(:, :)
      real(rp), intent(in), optional :: q_points(:, :)
      complex(rp), intent(in) :: g_ab(:, :, :, :), g_ba(:, :, :, :)
      integer, intent(in) :: pair_sites(:, :), site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_realspace_chi0_options), intent(in) :: options
      type(tddft_realspace_chi0_result), intent(out) :: result
      logical, allocatable :: keep(:)
      real(rp), allocatable :: radii(:)
      complex(rp), allocatable :: left_ops(:, :, :), right_ops(:, :, :), pair_response(:, :)
      real(rp), allocatable :: pair_norms(:)
      complex(rp), allocatable :: f0(:, :), f1(:, :), gr0(:, :), rev0(:, :), adv_reverse0(:, :), &
         adv_forward0(:, :), work1(:, :), work2(:, :), work3(:, :)
      integer :: ip, ikeep, nkeep, nleft, nright, nblock
      real(rp) :: effective_cutoff, source_radius, selected_norm, omitted_norm, total_norm
      real(rp) :: t_start, t_stop, total_start, total_stop
      logical :: full_tail_mode

      call validate_realspace_source(energy_grid, g_ab, g_ba, r_vectors, pair_sites, site_orbital_counts, &
         left_channels, right_channels)
      if (options%tail_tolerance < 0.0_rp .or. options%electronic_temperature < 0.0_rp) then
         error stop 'build_static_chi0_from_realspace_gf: invalid static options'
      end if
      full_tail_mode = realspace_full_tail_mode(options%truncation_mode)
      ! The exact static route uses the sampled retarded/advanced source and is
      ! independent of the dynamic direct-versus-contour selection.
      nleft = size(left_channels); nright = size(right_channels); nblock = size(g_ab, 1)
      call cpu_time(total_start)
      allocate(radii(size(r_vectors, 2)), keep(size(r_vectors, 2)))
      do ip = 1, size(r_vectors, 2)
         radii(ip) = metric_norm(options%metric, r_vectors(:, ip))
         keep(ip) = radii(ip) <= options%rmax + 16.0_rp*epsilon(1.0_rp)*max(1.0_rp, radii(ip))
      end do
      nkeep = count(keep)
      if (nkeep < 1 .and. size(r_vectors, 2) > 0 .and. realspace_parallel_size() == 1) then
         error stop 'build_static_chi0_from_realspace_gf: Rmax removes every real-space pair'
      end if
      if (nkeep > 0) then
         effective_cutoff = maxval(pack(radii, keep))
      else
         effective_cutoff = 0.0_rp
      end if
      if (size(radii) > 0) then
         source_radius = maxval(radii)
      else
         source_radius = 0.0_rp
      end if
      allocate(result%chi_r(nleft, nright, 1, nkeep), result%r_vectors(3, nkeep))
      result%chi_r = cmplx(0.0_rp, 0.0_rp, rp)
      ikeep = 0
      do ip = 1, size(keep)
         if (keep(ip)) then
            ikeep = ikeep + 1
            result%r_vectors(:, ikeep) = r_vectors(:, ip)
         end if
      end do
      allocate(left_ops(nblock, nblock, nleft), right_ops(nblock, nblock, nright))
      call make_local_operators(left_channels, site_orbital_counts, nblock, left_ops)
      call make_local_operators(right_channels, site_orbital_counts, nblock, right_ops)
      allocate(pair_response(nleft, nright), pair_norms(size(keep)))
      allocate(f0(nleft, nright), f1(nleft, nright), gr0(nblock, nblock), rev0(nblock, nblock), &
         adv_reverse0(nblock, nblock), adv_forward0(nblock, nblock), work1(nblock, nblock), &
         work2(nblock, nblock), work3(nblock, nblock))
      selected_norm = 0.0_rp; omitted_norm = 0.0_rp; total_norm = 0.0_rp; pair_norms = 0.0_rp
      ikeep = 0
      call cpu_time(t_start)
      do ip = 1, size(keep)
         if (.not. keep(ip) .and. .not. full_tail_mode) cycle
         call integrate_static_pair_response(energy_grid, g_ab(:, :, :, ip), g_ba(:, :, :, ip), pair_sites(ip, :), &
            left_channels, right_channels, left_ops, right_ops, options, f0, f1, gr0, rev0, adv_reverse0, adv_forward0, &
            work1, work2, work3, pair_response)
         pair_norms(ip) = sum(abs(pair_response))
         total_norm = total_norm + pair_norms(ip)
         if (keep(ip)) then
            ikeep = ikeep + 1
            result%chi_r(:, :, 1, ikeep) = pair_response
            selected_norm = selected_norm + pair_norms(ip)
         else
            omitted_norm = omitted_norm + pair_norms(ip)
         end if
      end do
      call cpu_time(t_stop)
      result%diagnostics%pair_integration_cpu_seconds = t_stop-t_start
      call cpu_time(t_start)
      if (present(q_points)) then
         call fourier_transform_realspace_chi0(result%chi_r, result%r_vectors, q_points, options%representation, &
            options%fourier_axes, result%chi_q)
      else
         allocate(result%chi_q(nleft, nright, 1, 1)); result%chi_q = cmplx(0.0_rp, 0.0_rp, rp)
      end if
      call cpu_time(t_stop)
      result%diagnostics%fourier_cpu_seconds = t_stop-t_start
      result%diagnostics%input_points = size(keep)
      result%diagnostics%integration_energy_min = energy_grid(1)
      result%diagnostics%integration_energy_max = energy_grid(size(energy_grid))
      result%diagnostics%selected_points = nkeep
      result%diagnostics%omitted_points = size(keep)-nkeep
      result%diagnostics%rmax = effective_cutoff
      result%diagnostics%requested_rmax = options%rmax
      result%diagnostics%source_radius = source_radius
      result%diagnostics%source_covers_requested_cutoff = source_radius + &
         16.0_rp*epsilon(1.0_rp)*max(1.0_rp, source_radius, abs(options%rmax)) >= options%rmax
      result%diagnostics%selected_norm = selected_norm
      result%diagnostics%omitted_tail_norm = omitted_norm
      result%diagnostics%total_norm = total_norm
      if (full_tail_mode .and. total_norm > tiny(1.0_rp)) result%diagnostics%tail_ratio = omitted_norm/total_norm
      result%diagnostics%tail_assessed = full_tail_mode .and. result%diagnostics%omitted_points > 0
      result%diagnostics%converged = result%diagnostics%tail_assessed .and. &
         result%diagnostics%tail_ratio <= options%tail_tolerance
      result%diagnostics%pair_response_integrations = nkeep
      if (full_tail_mode) result%diagnostics%pair_response_integrations = size(keep)
      if (.not. full_tail_mode .and. result%diagnostics%omitted_points > 0) then
         result%diagnostics%status = 'production truncation; omitted static tail not assessed'
      else if (.not. result%diagnostics%source_covers_requested_cutoff) then
         result%diagnostics%status = 'source radius does not extend beyond requested R_max; tail not assessed'
      else if (.not. result%diagnostics%tail_assessed) then
         result%diagnostics%status = 'all supplied real-space pairs retained'
      else if (result%diagnostics%converged) then
         result%diagnostics%status = 'R-space static tail below tolerance'
      else
         result%diagnostics%status = 'R-space static tail exceeds tolerance'
      end if
      result%diagnostics%last_shell_norm = last_shell_norm(radii, keep, pair_norms)
      result%diagnostics%shell_count = count_realspace_shells(radii)
      call cpu_time(total_stop)
      result%diagnostics%total_cpu_seconds = total_stop-total_start
   end subroutine build_static_chi0_from_realspace_gf

   subroutine fourier_transform_realspace_chi0(chi_r, r_vectors, q_points, representation, fourier_axes, chi_q)
      complex(rp), intent(in) :: chi_r(:, :, :, :)
      real(rp), intent(in) :: r_vectors(:, :), q_points(:, :)
      character(len=*), intent(in) :: representation
      integer, intent(in) :: fourier_axes(:)
      complex(rp), allocatable, intent(out) :: chi_q(:, :, :, :)
      integer :: nleft, nright, nw, nr, nq, ir, iq, ia
      complex(rp) :: phase
      real(rp) :: dot_product_qr

      call validate_fourier_inputs(chi_r, r_vectors, q_points, representation, fourier_axes)
      nleft = size(chi_r, 1); nright = size(chi_r, 2); nw = size(chi_r, 3)
      nr = size(chi_r, 4); nq = size(q_points, 2)
      allocate(chi_q(nleft, nright, nw, nq)); chi_q = cmplx(0.0_rp, 0.0_rp, rp)
      do iq = 1, nq
         do ir = 1, nr
            dot_product_qr = 0.0_rp
            if (trim(representation) /= 'finite' .and. trim(representation) /= 'local') then
               do ia = 1, size(fourier_axes)
                  if (fourier_axes(ia) > 0) dot_product_qr = dot_product_qr + &
                     q_points(fourier_axes(ia), iq)*r_vectors(fourier_axes(ia), ir)
               end do
               phase = exp(-i_unit*two_pi*dot_product_qr)
            else
               phase = cmplx(1.0_rp, 0.0_rp, rp)
            end if
            chi_q(:, :, :, iq) = chi_q(:, :, :, iq) + phase*chi_r(:, :, :, ir)
         end do
      end do
   end subroutine fourier_transform_realspace_chi0

   !> This transform is a validation utility only.  Production TDDFT-07 uses
   !> the susceptibility transform above and does not route G(R,z) through it.
   subroutine fourier_transform_realspace_green(g_r, r_vectors, q_points, g_q)
      complex(rp), intent(in) :: g_r(:, :, :, :)
      real(rp), intent(in) :: r_vectors(:, :), q_points(:, :)
      complex(rp), allocatable, intent(out) :: g_q(:, :, :, :)
      integer :: ie, ir, iq
      complex(rp) :: phase
      real(rp) :: dot_product_qr

      if (size(r_vectors, 1) /= 3 .or. size(q_points, 1) /= 3 .or. size(g_r, 4) /= size(r_vectors, 2)) then
         error stop 'fourier_transform_realspace_green: incompatible real-space dimensions'
      end if
      allocate(g_q(size(g_r, 1), size(g_r, 2), size(g_r, 3), size(q_points, 2)))
      g_q = cmplx(0.0_rp, 0.0_rp, rp)
      do iq = 1, size(q_points, 2)
         do ir = 1, size(r_vectors, 2)
            dot_product_qr = dot_product(q_points(:, iq), r_vectors(:, ir))
            phase = exp(-i_unit*two_pi*dot_product_qr)
            do ie = 1, size(g_r, 3)
               g_q(:, :, ie, iq) = g_q(:, :, ie, iq) + phase*g_r(:, :, ie, ir)
            end do
         end do
      end do
   end subroutine fourier_transform_realspace_green

   !> @brief Replicate the complete native-RGF susceptibility on every rank.
   !> @details The native provider owns disjoint real-space pairs, so its
   !> `chi_q` is a rank-local partial sum after the local R-space integration
   !> and Fourier transform.  The Fourier phase is linear, hence reducing
   !> `chi_q` after the phase factors is algebraically equivalent to reducing
   !> `chi_r` first, without gathering variable-length R lists.  No division by
   !> the number of ranks is performed: each pair is owned by exactly one rank
   !> and contributes exactly once to the SUM.
   !>
   !> The derived spectral arrays are reduced alongside `chi` for consumers
   !> that use them directly.  `re_chi` and `im_chi` are regenerated from the
   !> reduced complex response so the complex susceptibility remains the sole
   !> normalization source of truth.
   subroutine reduce_realspace_chi0_batch(batch)
      type(tddft_chi0_batch_result), intent(inout) :: batch
#ifdef USE_MPI
      integer :: iq
      integer :: local_counts(2), global_counts(2)
      integer :: local_shell_count, global_shell_count
      integer :: local_pair_integrations, global_pair_integrations
      integer :: local_tail_assessed, global_tail_assessed
      real(rp) :: local_tail, global_tail, local_total, global_total
      real(rp) :: local_selected, global_selected, local_cutoff, global_cutoff
      real(rp) :: local_requested, global_requested, local_source, global_source
      real(rp) :: local_green_time, global_green_time, local_pair_time, global_pair_time
      real(rp) :: local_fourier_time, global_fourier_time, local_total_time, global_total_time
      real(rp) :: local_last_shell, global_last_shell
      integer :: local_source_covers, global_source_covers

      if (.not. allocated(batch%q_response)) return
      do iq = 1, size(batch%q_response)
         if (allocated(batch%q_response(iq)%chi) .and. size(batch%q_response(iq)%chi) > 0) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, batch%q_response(iq)%chi, size(batch%q_response(iq)%chi), &
               MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
            if (allocated(batch%q_response(iq)%re_chi)) batch%q_response(iq)%re_chi = real(batch%q_response(iq)%chi, rp)
            if (allocated(batch%q_response(iq)%im_chi)) batch%q_response(iq)%im_chi = aimag(batch%q_response(iq)%chi)
         end if
         if (allocated(batch%q_response(iq)%site_diagonal_spectrum) .and. &
             size(batch%q_response(iq)%site_diagonal_spectrum) > 0) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, batch%q_response(iq)%site_diagonal_spectrum, &
               size(batch%q_response(iq)%site_diagonal_spectrum), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         end if
         if (allocated(batch%q_response(iq)%stoner_spectral_map) .and. &
             size(batch%q_response(iq)%stoner_spectral_map) > 0) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, batch%q_response(iq)%stoner_spectral_map, &
               size(batch%q_response(iq)%stoner_spectral_map), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         end if
         if (allocated(batch%q_response(iq)%trace_spectrum) .and. size(batch%q_response(iq)%trace_spectrum) > 0) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, batch%q_response(iq)%trace_spectrum, &
               size(batch%q_response(iq)%trace_spectrum), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         end if

         ! Pair counts and tail norms are additive over disjoint local pair
         ! ownership.  The cutoff is a maximum, not a sum.  Shell counts and
         ! last-shell norms are retained as conservative rank-wide summaries;
         ! they are diagnostics only and do not enter chi0 normalization.
         local_counts = [batch%q_response(iq)%metadata%real_space_points, &
            batch%q_response(iq)%metadata%real_space_omitted_points]
         call MPI_ALLREDUCE(local_counts, global_counts, 2, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_points = global_counts(1)
         batch%q_response(iq)%metadata%real_space_omitted_points = global_counts(2)
         local_tail = batch%q_response(iq)%metadata%real_space_tail_norm
         call MPI_ALLREDUCE(local_tail, global_tail, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_tail_norm = global_tail
         local_selected = batch%q_response(iq)%metadata%real_space_selected_norm
         call MPI_ALLREDUCE(local_selected, global_selected, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_selected_norm = global_selected
         local_total = batch%q_response(iq)%metadata%real_space_total_norm
         call MPI_ALLREDUCE(local_total, global_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_total_norm = global_total
         if (global_total > tiny(1.0_rp)) then
            batch%q_response(iq)%metadata%real_space_tail_ratio = global_tail/global_total
         else
            batch%q_response(iq)%metadata%real_space_tail_ratio = 0.0_rp
         end if
         local_cutoff = batch%q_response(iq)%metadata%real_space_cutoff
         call MPI_ALLREDUCE(local_cutoff, global_cutoff, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_cutoff = global_cutoff
         local_requested = batch%q_response(iq)%metadata%real_space_requested_cutoff
         call MPI_ALLREDUCE(local_requested, global_requested, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_requested_cutoff = global_requested
         local_source = batch%q_response(iq)%metadata%real_space_source_radius
         call MPI_ALLREDUCE(local_source, global_source, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_source_radius = global_source
         local_source_covers = merge(1, 0, batch%q_response(iq)%metadata%real_space_source_covers_cutoff)
         call MPI_ALLREDUCE(local_source_covers, global_source_covers, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_source_covers_cutoff = global_source_covers /= 0
         local_pair_integrations = batch%q_response(iq)%metadata%real_space_pair_response_integrations
         call MPI_ALLREDUCE(local_pair_integrations, global_pair_integrations, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_pair_response_integrations = global_pair_integrations
         local_pair_time = batch%q_response(iq)%metadata%real_space_pair_integration_cpu_seconds
         call MPI_ALLREDUCE(local_pair_time, global_pair_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_pair_integration_cpu_seconds = global_pair_time
         local_green_time = batch%q_response(iq)%metadata%real_space_green_cpu_seconds
         call MPI_ALLREDUCE(local_green_time, global_green_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_green_cpu_seconds = global_green_time
         local_fourier_time = batch%q_response(iq)%metadata%real_space_fourier_cpu_seconds
         call MPI_ALLREDUCE(local_fourier_time, global_fourier_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_fourier_cpu_seconds = global_fourier_time
         local_total_time = batch%q_response(iq)%metadata%real_space_total_cpu_seconds
         call MPI_ALLREDUCE(local_total_time, global_total_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_total_cpu_seconds = global_total_time
         local_tail_assessed = merge(1, 0, batch%q_response(iq)%metadata%real_space_tail_assessed)
         call MPI_ALLREDUCE(local_tail_assessed, global_tail_assessed, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_tail_assessed = global_tail_assessed /= 0
         batch%q_response(iq)%metadata%converged = batch%q_response(iq)%metadata%real_space_tail_assessed .and. &
            batch%q_response(iq)%metadata%real_space_tail_ratio <= batch%q_response(iq)%metadata%real_space_tail_tolerance
         if (trim(batch%q_response(iq)%metadata%real_space_truncation_mode) == 'production' .or. &
             trim(batch%q_response(iq)%metadata%real_space_truncation_mode) == 'truncated') then
            if (global_counts(2) > 0) then
               batch%q_response(iq)%metadata%convergence_status = 'production truncation; omitted tail not assessed'
            else
               batch%q_response(iq)%metadata%convergence_status = &
                  'source radius does not extend beyond requested R_max; tail not assessed'
            end if
            batch%q_response(iq)%metadata%converged = .false.
         else if (.not. batch%q_response(iq)%metadata%real_space_tail_assessed) then
            batch%q_response(iq)%metadata%convergence_status = &
               'source radius does not extend beyond requested R_max; tail not assessed'
            batch%q_response(iq)%metadata%converged = .false.
         else if (batch%q_response(iq)%metadata%converged) then
            batch%q_response(iq)%metadata%convergence_status = 'R-space tail below tolerance'
         else
            batch%q_response(iq)%metadata%convergence_status = 'R-space tail exceeds tolerance'
         end if
         local_shell_count = batch%q_response(iq)%metadata%real_space_shell_count
         call MPI_ALLREDUCE(local_shell_count, global_shell_count, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_shell_count = global_shell_count
         local_last_shell = batch%q_response(iq)%metadata%real_space_last_shell_norm
         call MPI_ALLREDUCE(local_last_shell, global_last_shell, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_last_shell_norm = global_last_shell
      end do
      if (size(batch%q_response) > 0) batch%metadata = batch%q_response(1)%metadata
#else
      ! Serial builds already hold the complete real-space source locally.
      continue
#endif
   end subroutine reduce_realspace_chi0_batch

   subroutine check_realspace_pair_reversal(pair_sites, r_vectors, tolerance, all_pairs_have_reverse, max_residual)
      integer, intent(in) :: pair_sites(:, :)
      real(rp), intent(in) :: r_vectors(:, :), tolerance
      logical, intent(out) :: all_pairs_have_reverse
      real(rp), intent(out) :: max_residual
      integer :: ip, jp
      real(rp) :: residual

      if (size(pair_sites, 2) /= 2 .or. size(r_vectors, 1) /= 3 .or. size(pair_sites, 1) /= size(r_vectors, 2)) then
         error stop 'check_realspace_pair_reversal: incompatible pair dimensions'
      end if
      all_pairs_have_reverse = .true.; max_residual = 0.0_rp
      do ip = 1, size(pair_sites, 1)
         residual = huge(1.0_rp)
         do jp = 1, size(pair_sites, 1)
            if (pair_sites(jp, 1) == pair_sites(ip, 2) .and. pair_sites(jp, 2) == pair_sites(ip, 1)) then
               residual = min(residual, maxval(abs(r_vectors(:, jp)+r_vectors(:, ip))))
            end if
         end do
         max_residual = max(max_residual, residual)
         if (residual > tolerance) all_pairs_have_reverse = .false.
      end do
   end subroutine check_realspace_pair_reversal

   subroutine make_response_from_q(chi, q_point, omega, left_channels, right_channels, options, diagnostics, result)
      complex(rp), intent(in) :: chi(:, :, :)
      real(rp), intent(in) :: q_point(3), omega(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_realspace_chi0_options), intent(in) :: options
      type(tddft_realspace_chi0_diagnostics), intent(in) :: diagnostics
      type(tddft_chi0_result), intent(out) :: result
      integer :: nleft, nright, nw, ich, iw

      nleft = size(chi, 1); nright = size(chi, 2); nw = size(chi, 3)
      allocate(result%chi(nleft, nright, nw), result%re_chi(nleft, nright, nw), &
         result%im_chi(nleft, nright, nw), result%site_diagonal_spectrum(min(nleft, nright), nw), &
         result%trace_spectrum(nw), result%stoner_spectral_map(min(nleft, nright), nw))
      result%chi = chi; result%re_chi = real(chi, rp); result%im_chi = aimag(chi)
      result%site_diagonal_spectrum = 0.0_rp; result%stoner_spectral_map = 0.0_rp; result%trace_spectrum = 0.0_rp
      do iw = 1, nw
         do ich = 1, min(nleft, nright)
            if (left_channels(ich)%site /= right_channels(ich)%site) cycle
            result%site_diagonal_spectrum(ich, iw) = -aimag(chi(ich, ich, iw))/pi
            result%stoner_spectral_map(ich, iw) = result%site_diagonal_spectrum(ich, iw)
            result%trace_spectrum(iw) = result%trace_spectrum(iw) + result%site_diagonal_spectrum(ich, iw)
         end do
      end do
      result%metadata%backend = 'realspace_gf'
      result%metadata%canonical_backend = 'realspace_gf'
      result%metadata%implementation = 'native G(R,z) -> chi0(R,omega) -> susceptibility FT'
      result%metadata%energy_integration = 'direct near-real-axis trapezoid'
      if (diagnostics%contour_deformation) then
         result%metadata%implementation = 'native complex G(R,z) -> chi0(R,omega) -> susceptibility FT'
         result%metadata%energy_integration = 'mixed contour: RR/AA deformation plus near-Fermi RA'
      end if
      result%metadata%endpoint_provenance = 'native G_ab(R,z), G_ba(-R,z); no G(k) endpoint'
      result%metadata%response_projection = 'site'
      result%metadata%circular_channel = options%circular_channel
      result%metadata%spectral_convention = 'Stoner spectral weight = -Im chi_KS^'//trim(options%circular_channel)//'/pi (positive excitation positive)'
      result%metadata%eta = options%eta
      result%metadata%green_eta = resolved_green_eta(options)
      result%metadata%fermi_level = options%fermi_level
      result%metadata%electronic_temperature = options%electronic_temperature
      result%metadata%electronic_kT = max(options%electronic_temperature*tddft_kB_Ry_per_K, tddft_occupation_kT_floor)
      result%metadata%q_direct = q_point
      result%metadata%omega_min = minval(omega); result%metadata%omega_max = maxval(omega); result%metadata%omega_points = nw
      result%metadata%q_batch_size = 1; result%metadata%omega_batch_size = nw
      result%metadata%real_space_reuse = .true.
      result%metadata%real_space_points = diagnostics%selected_points
      result%metadata%real_space_cutoff = diagnostics%rmax
      result%metadata%real_space_requested_cutoff = diagnostics%requested_rmax
      result%metadata%real_space_source_radius = diagnostics%source_radius
      result%metadata%real_space_representation = options%representation
      result%metadata%real_space_truncation_mode = options%truncation_mode
      result%metadata%real_space_fourier_axes = options%fourier_axes
      result%metadata%real_space_omitted_points = diagnostics%omitted_points
      result%metadata%real_space_shell_count = diagnostics%shell_count
      result%metadata%real_space_selected_norm = diagnostics%selected_norm
      result%metadata%real_space_tail_norm = diagnostics%omitted_tail_norm
      result%metadata%real_space_total_norm = diagnostics%total_norm
      result%metadata%real_space_tail_ratio = diagnostics%tail_ratio
      result%metadata%real_space_last_shell_norm = diagnostics%last_shell_norm
      result%metadata%real_space_pair_integration_cpu_seconds = diagnostics%pair_integration_cpu_seconds
      result%metadata%real_space_fourier_cpu_seconds = diagnostics%fourier_cpu_seconds
      result%metadata%real_space_green_cpu_seconds = diagnostics%green_function_cpu_seconds
      result%metadata%real_space_total_cpu_seconds = diagnostics%total_cpu_seconds
      result%metadata%real_space_pair_response_integrations = diagnostics%pair_response_integrations
      result%metadata%real_space_tail_tolerance = options%tail_tolerance
      result%metadata%real_space_tail_assessed = diagnostics%tail_assessed
      result%metadata%real_space_source_covers_cutoff = diagnostics%source_covers_requested_cutoff
      result%metadata%convergence_status = diagnostics%status
      result%metadata%converged = diagnostics%converged
      result%metadata%integration_energy_min = diagnostics%integration_energy_min
      result%metadata%integration_energy_max = diagnostics%integration_energy_max
      result%metadata%integration_energy_points = 0
      result%metadata%green_function_evaluations = diagnostics%contour_gf_evaluations
      result%metadata%contour_points = diagnostics%contour_points
      result%metadata%contour_subdivisions = diagnostics%contour_subdivisions
      result%metadata%near_fermi_points = diagnostics%near_fermi_points
      result%metadata%contour_height = diagnostics%contour_height
      result%metadata%contour_max_imaginary_energy = diagnostics%contour_max_imaginary_energy
      result%metadata%contour_gf_evaluations = diagnostics%contour_gf_evaluations
      result%metadata%contour_deformation = diagnostics%contour_deformation
      if (diagnostics%contour_deformation) result%metadata%integration_energy_points = &
         2*(2+diagnostics%contour_subdivisions)*diagnostics%contour_points + 3*diagnostics%near_fermi_points
   end subroutine make_response_from_q

   subroutine mark_static_realspace_metadata(result, q_point, nq)
      type(tddft_chi0_result), intent(inout) :: result
      real(rp), intent(in) :: q_point(3)
      integer, intent(in) :: nq

      result%metadata%backend = 'realspace_gf'
      result%metadata%canonical_backend = 'realspace_gf'
      result%metadata%implementation = 'native G(R,E) retarded/advanced static contour identity'
      result%metadata%energy_integration = 'real-axis static retarded/advanced product'
      result%metadata%endpoint_provenance = 'native G_ab(R,E), G_ba(-R,E); exact static contour identity'
      result%metadata%frequency_convention = 'static omega=0 divided difference; no dynamical eta'
      result%metadata%q_direct = q_point
      result%metadata%eta = 0.0_rp
      result%metadata%green_eta = 0.0_rp
      result%metadata%eta_role = 'not used by exact static response'
      result%metadata%eta_is_numerical = .false.
      result%metadata%static_limit = .true.
      result%metadata%omega_min = 0.0_rp
      result%metadata%omega_max = 0.0_rp
      result%metadata%omega_points = 1
      result%metadata%q_batch_size = nq
      result%metadata%omega_batch_size = 1
      result%metadata%integration_energy_points = 0
   end subroutine mark_static_realspace_metadata

   subroutine integrate_static_pair_response(energy_grid, g_ab, g_ba, pair_sites, left_channels, right_channels, &
      left_ops, right_ops, options, f0, f1, gr0, rev0, adv_reverse0, adv_forward0, work1, work2, work3, response)
      real(rp), intent(in) :: energy_grid(:)
      complex(rp), intent(in) :: g_ab(:, :, :), g_ba(:, :, :), left_ops(:, :, :), right_ops(:, :, :)
      integer, intent(in) :: pair_sites(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_realspace_chi0_options), intent(in) :: options
      complex(rp), intent(inout) :: f0(:, :), f1(:, :), gr0(:, :), rev0(:, :), adv_reverse0(:, :), adv_forward0(:, :), &
         work1(:, :), work2(:, :), work3(:, :)
      complex(rp), intent(out) :: response(:, :)
      integer :: ie
      real(rp) :: h

      response = cmplx(0.0_rp, 0.0_rp, rp)
      do ie = 1, size(energy_grid)-1
         call static_response_integrand(energy_grid(ie), energy_grid, g_ab, g_ba, pair_sites, left_channels, right_channels, &
            left_ops, right_ops, options, gr0, rev0, adv_reverse0, adv_forward0, work1, work2, work3, f0)
         call static_response_integrand(energy_grid(ie+1), energy_grid, g_ab, g_ba, pair_sites, left_channels, right_channels, &
            left_ops, right_ops, options, gr0, rev0, adv_reverse0, adv_forward0, work1, work2, work3, f1)
         h = energy_grid(ie+1)-energy_grid(ie)
         response = response + 0.5_rp*h*(f0+f1)
      end do
   end subroutine integrate_static_pair_response

   subroutine static_response_integrand(energy, energy_grid, g_ab, g_ba, pair_sites, left_channels, right_channels, &
      left_ops, right_ops, options, gr0, rev0, adv_reverse0, adv_forward0, work1, work2, work3, response)
      real(rp), intent(in) :: energy, energy_grid(:)
      complex(rp), intent(in) :: g_ab(:, :, :), g_ba(:, :, :), left_ops(:, :, :), right_ops(:, :, :)
      integer, intent(in) :: pair_sites(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_realspace_chi0_options), intent(in) :: options
      complex(rp), intent(inout) :: gr0(:, :), rev0(:, :), adv_reverse0(:, :), adv_forward0(:, :), work1(:, :), work2(:, :), work3(:, :)
      complex(rp), intent(out) :: response(:, :)
      complex(rp) :: retarded_product, advanced_product
      integer :: il, ir
      real(rp) :: occupation

      call interpolate_block(energy_grid, g_ab, energy, gr0)
      call interpolate_block(energy_grid, g_ba, energy, rev0)
      adv_reverse0 = transpose(conjg(gr0))
      adv_forward0 = transpose(conjg(rev0))
      response = cmplx(0.0_rp, 0.0_rp, rp)
      occupation = tddft_fermi_occupation(energy, options%fermi_level, options%electronic_temperature)
      do il = 1, size(left_channels)
         if (left_channels(il)%site /= pair_sites(1)) cycle
         do ir = 1, size(right_channels)
            if (right_channels(ir)%site /= pair_sites(2)) cycle
            retarded_product = trace_four_local(left_ops(:, :, il), gr0, right_ops(:, :, ir), rev0, work1, work2, work3)
            advanced_product = trace_four_local(left_ops(:, :, il), adv_forward0, right_ops(:, :, ir), adv_reverse0, &
               work1, work2, work3)
            response(il, ir) = -occupation*(retarded_product-advanced_product)/(2.0_rp*pi*i_unit)
         end do
      end do
   end subroutine static_response_integrand

   subroutine integrate_pair_response(energy_grid, g_ab, g_ba, pair_sites, left_channels, right_channels, &
      left_ops, right_ops, omega, options, f0, f1, gr0, gr1, grm, rev0, adv_reverse0, adv_forward0, advm, work1, work2, &
      work3, response)
      real(rp), intent(in) :: energy_grid(:), omega
      complex(rp), intent(in) :: g_ab(:, :, :), g_ba(:, :, :), left_ops(:, :, :), right_ops(:, :, :)
      integer, intent(in) :: pair_sites(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_realspace_chi0_options), intent(in) :: options
      complex(rp), intent(inout) :: f0(:, :), f1(:, :), gr0(:, :), gr1(:, :), grm(:, :), rev0(:, :), &
         adv_reverse0(:, :), adv_forward0(:, :), advm(:, :), work1(:, :), work2(:, :), work3(:, :)
      complex(rp), intent(out) :: response(:, :)
      integer :: ie
      real(rp) :: lo, hi, h

      response = cmplx(0.0_rp, 0.0_rp, rp)
      lo = energy_grid(1) + abs(omega); hi = energy_grid(size(energy_grid)) - abs(omega)
      if (hi <= lo) error stop 'integrate_pair_response: Green energy grid does not cover omega'
      do ie = 1, size(energy_grid)-1
         if (energy_grid(ie) < lo .or. energy_grid(ie+1) > hi) cycle
         call response_integrand(energy_grid(ie), energy_grid, g_ab, g_ba, pair_sites, left_channels, right_channels, &
            left_ops, right_ops, omega, options, gr0, gr1, grm, rev0, adv_reverse0, adv_forward0, advm, work1, work2, work3, f0)
         call response_integrand(energy_grid(ie+1), energy_grid, g_ab, g_ba, pair_sites, left_channels, right_channels, &
            left_ops, right_ops, omega, options, gr0, gr1, grm, rev0, adv_reverse0, adv_forward0, advm, work1, work2, work3, f1)
         h = energy_grid(ie+1)-energy_grid(ie)
         response = response + 0.5_rp*h*(f0+f1)
      end do
   end subroutine integrate_pair_response

   subroutine response_integrand(energy, energy_grid, g_ab, g_ba, pair_sites, left_channels, right_channels, &
      left_ops, right_ops, omega, options, gr0, gr1, grm, rev0, adv_reverse0, adv_forward0, advm, work1, work2, work3, response)
      real(rp), intent(in) :: energy, energy_grid(:), omega
      complex(rp), intent(in) :: g_ab(:, :, :), g_ba(:, :, :), left_ops(:, :, :), right_ops(:, :, :)
      integer, intent(in) :: pair_sites(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_realspace_chi0_options), intent(in) :: options
      complex(rp), intent(inout) :: gr0(:, :), gr1(:, :), grm(:, :), rev0(:, :), adv_reverse0(:, :), adv_forward0(:, :), advm(:, :)
      complex(rp), intent(inout) :: work1(:, :), work2(:, :), work3(:, :)
      complex(rp), intent(out) :: response(:, :)
      complex(rp) :: bubble
      integer :: il, ir
      real(rp) :: occupation

      call interpolate_block(energy_grid, g_ab, energy, gr0)
      call interpolate_block(energy_grid, g_ba, energy, rev0)
      call interpolate_block(energy_grid, g_ab, energy+omega, gr1)
      call interpolate_block(energy_grid, g_ab, energy-omega, grm)
      adv_reverse0 = transpose(conjg(gr0))
      adv_forward0 = transpose(conjg(rev0))
      advm = transpose(conjg(grm))
      response = cmplx(0.0_rp, 0.0_rp, rp)
      occupation = tddft_fermi_occupation(energy, options%fermi_level, options%electronic_temperature)
      do il = 1, size(left_channels)
         if (left_channels(il)%site /= pair_sites(1)) cycle
         do ir = 1, size(right_channels)
            if (right_channels(ir)%site /= pair_sites(2)) cycle
            bubble = trace_four_local(left_ops(:, :, il), gr1, right_ops(:, :, ir), rev0-adv_reverse0, work1, work2, work3) + &
               trace_four_local(left_ops(:, :, il), gr0-adv_forward0, right_ops(:, :, ir), advm, work1, work2, work3)
            response(il, ir) = -occupation*bubble/(2.0_rp*pi*i_unit)
         end do
      end do
   end subroutine response_integrand

   subroutine interpolate_block(grid, values, x, block)
      real(rp), intent(in) :: grid(:), x
      complex(rp), intent(in) :: values(:, :, :)
      complex(rp), intent(out) :: block(:, :)
      integer :: ie
      real(rp) :: weight

      if (x < grid(1)-32.0_rp*epsilon(1.0_rp) .or. x > grid(size(grid))+32.0_rp*epsilon(1.0_rp)) then
         error stop 'interpolate_block: requested energy is outside native Green grid'
      end if
      if (x <= grid(1)) then
         block = values(:, :, 1); return
      else if (x >= grid(size(grid))) then
         block = values(:, :, size(grid, 1)); return
      end if
      ie = 1
      do while (ie < size(grid)-1 .and. grid(ie+1) < x)
         ie = ie + 1
      end do
      weight = (x-grid(ie))/(grid(ie+1)-grid(ie))
      block = (1.0_rp-weight)*values(:, :, ie) + weight*values(:, :, ie+1)
   end subroutine interpolate_block

   function trace_four_local(a, b, c, d, work1, work2, work3) result(value)
      complex(rp), intent(in) :: a(:, :), b(:, :), c(:, :), d(:, :)
      complex(rp), intent(out) :: work1(:, :), work2(:, :), work3(:, :)
      complex(rp) :: value
      integer :: ie

      work1 = matmul(a, b); work2 = matmul(work1, c); work3 = matmul(work2, d)
      value = cmplx(0.0_rp, 0.0_rp, rp)
      do ie = 1, size(work3, 1)
         value = value + work3(ie, ie)
      end do
   end function trace_four_local

   subroutine make_local_operators(channels, site_orbital_counts, nblock, local_ops)
      type(response_channel), intent(in) :: channels(:)
      integer, intent(in) :: site_orbital_counts(:), nblock
      complex(rp), intent(out) :: local_ops(:, :, :)
      complex(rp), allocatable :: full_op(:, :)
      integer :: ich, offset, site_nblock

      if (nblock /= 2*site_orbital_counts(1)) error stop &
         'make_local_operators: native provider requires equal orbital block sizes'
      do ich = 1, size(channels)
         full_op = site_projected_operator(channels(ich), site_orbital_counts)
         offset = site_block_offset(channels(ich)%site, site_orbital_counts)
         site_nblock = 2*site_orbital_counts(channels(ich)%site)
         if (site_nblock /= nblock) error stop 'make_local_operators: pair blocks have different sizes'
         local_ops(:, :, ich) = full_op(offset+1:offset+nblock, offset+1:offset+nblock)
      end do
   end subroutine make_local_operators

   integer function site_block_offset(site, counts) result(offset)
      integer, intent(in) :: site, counts(:)
      if (site < 1 .or. site > size(counts)) error stop 'site_block_offset: site is outside response layout'
      offset = 2*sum(counts(1:site-1))
   end function site_block_offset

   real(rp) function metric_norm(metric, vector) result(norm)
      real(rp), intent(in) :: metric(:, :), vector(:)
      real(rp) :: cartesian(3)
      cartesian = matmul(metric, vector)
      norm = sqrt(dot_product(cartesian, cartesian))
   end function metric_norm

   real(rp) function last_shell_norm(radii, keep, pair_norms) result(value)
      real(rp), intent(in) :: radii(:), pair_norms(:)
      logical, intent(in) :: keep(:)
      real(rp) :: outer_radius
      integer :: ip

      value = 0.0_rp
      if (count(keep) == 0) return
      outer_radius = maxval(pack(radii, keep))
      do ip = 1, size(keep)
         if (keep(ip) .and. radii(ip) >= outer_radius-16.0_rp*epsilon(1.0_rp)*max(1.0_rp, outer_radius)) then
            value = value + pair_norms(ip)
         end if
      end do
   end function last_shell_norm

   integer function count_realspace_shells(radii) result(count)
      real(rp), intent(in) :: radii(:)
      real(rp), allocatable :: unique_radii(:)
      integer :: ip, nu
      logical :: found

      allocate(unique_radii(size(radii))); nu = 0
      do ip = 1, size(radii)
         found = .false.
         if (nu > 0) found = any(abs(unique_radii(1:nu)-radii(ip)) <= 32.0_rp*epsilon(1.0_rp)*max(1.0_rp, radii(ip)))
         if (.not. found) then
            nu = nu + 1; unique_radii(nu) = radii(ip)
         end if
      end do
      count = nu
   end function count_realspace_shells

   subroutine validate_realspace_geometry(r_vectors, pair_sites, site_orbital_counts, left_channels, right_channels)
      real(rp), intent(in) :: r_vectors(:, :)
      integer, intent(in) :: pair_sites(:, :), site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      integer :: ich

      if (size(r_vectors, 1) /= 3 .or. size(pair_sites, 1) /= size(r_vectors, 2) .or. size(pair_sites, 2) /= 2 .or. &
          size(site_orbital_counts) < 1 .or. size(left_channels) < 1 .or. size(left_channels) /= size(right_channels)) then
         error stop 'native mixed contour: incompatible real-space geometry or response channels'
      end if
      if (any(site_orbital_counts < 1)) error stop 'native mixed contour: orbital counts must be positive'
      do ich = 1, size(pair_sites, 1)
         if (any(pair_sites(ich, :) < 1) .or. any(pair_sites(ich, :) > size(site_orbital_counts))) then
            error stop 'native mixed contour: pair site is outside response layout'
         end if
      end do
      do ich = 1, size(left_channels)
         if (left_channels(ich)%site < 1 .or. left_channels(ich)%site > size(site_orbital_counts) .or. &
             right_channels(ich)%site < 1 .or. right_channels(ich)%site > size(site_orbital_counts)) then
            error stop 'native mixed contour: response channel site is outside response layout'
         end if
      end do
   end subroutine validate_realspace_geometry

   subroutine mixed_native_same_interval(options, energy_min, energy_max, same_min, same_max)
      type(tddft_realspace_chi0_options), intent(in) :: options
      real(rp), intent(in) :: energy_min, energy_max
      real(rp), intent(out) :: same_min, same_max

      same_min = energy_min
      same_max = energy_max
      if (options%electronic_temperature <= 0.0_rp) same_max = min(energy_max, options%fermi_level)
   end subroutine mixed_native_same_interval

   subroutine mixed_native_contour_height(options, energy_width, green_eta, height)
      type(tddft_realspace_chi0_options), intent(in) :: options
      real(rp), intent(in) :: energy_width, green_eta
      real(rp), intent(out) :: height
      real(rp) :: kT

      kT = max(options%electronic_temperature*tddft_kB_Ry_per_K, tddft_occupation_kT_floor)
      if (options%contour_height > 0.0_rp) then
         height = options%contour_height
      else if (options%electronic_temperature > 0.0_rp) then
         height = min(0.25_rp*pi*kT, 0.25_rp*max(energy_width, green_eta))
      else
         height = 0.25_rp*max(energy_width, 4.0_rp*green_eta)
      end if
      if (height <= 0.0_rp) error stop 'native mixed contour: contour height must be positive'
      if (options%electronic_temperature > 0.0_rp .and. height >= 0.95_rp*pi*kT) then
         error stop 'native mixed contour: contour height crosses a Fermi-function pole; reduce contour_height'
      end if
   end subroutine mixed_native_contour_height

   pure complex(rp) function complex_fermi_occupation(energy, fermi_level, temperature) result(occupation)
      complex(rp), intent(in) :: energy
      real(rp), intent(in) :: fermi_level, temperature
      complex(rp) :: argument, exponential
      real(rp) :: kT

      kT = max(temperature*tddft_kB_Ry_per_K, tddft_occupation_kT_floor)
      argument = (energy-fermi_level)/kT
      if (real(argument, rp) >= 0.0_rp) then
         exponential = exp(-argument)
         occupation = exponential/(1.0_rp+exponential)
      else
         exponential = exp(argument)
         occupation = 1.0_rp/(1.0_rp+exponential)
      end if
   end function complex_fermi_occupation

   subroutine validate_realspace_source(energy_grid, g_ab, g_ba, r_vectors, pair_sites, site_orbital_counts, &
      left_channels, right_channels)
      real(rp), intent(in) :: energy_grid(:), r_vectors(:, :)
      complex(rp), intent(in) :: g_ab(:, :, :, :), g_ba(:, :, :, :)
      integer, intent(in) :: pair_sites(:, :), site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      integer :: ie

      if (size(energy_grid) < 2 .or. size(g_ab, 1) /= size(g_ab, 2) .or. size(g_ba, 1) /= size(g_ab, 1) .or. &
          size(g_ba, 2) /= size(g_ab, 2) .or. size(g_ab, 3) /= size(energy_grid) .or. &
          size(g_ba, 3) /= size(energy_grid) .or. size(g_ba, 4) /= size(g_ab, 4) .or. &
          size(r_vectors, 1) /= 3 .or. size(r_vectors, 2) /= size(g_ab, 4) .or. size(pair_sites, 1) /= size(g_ab, 4) .or. &
          size(pair_sites, 2) /= 2 .or. size(site_orbital_counts) < 1 .or. size(left_channels) < 1 .or. &
          size(left_channels) /= size(right_channels)) then
         error stop 'native real-space provider: incompatible source arrays'
      end if
      if (any(site_orbital_counts < 1)) error stop 'native real-space provider: orbital counts must be positive'
      do ie = 1, size(energy_grid)-1
         if (energy_grid(ie+1) <= energy_grid(ie)) error stop 'native real-space provider: energy grid is not increasing'
      end do
   end subroutine validate_realspace_source

   subroutine validate_fourier_inputs(chi_r, r_vectors, q_points, representation, fourier_axes)
      complex(rp), intent(in) :: chi_r(:, :, :, :)
      real(rp), intent(in) :: r_vectors(:, :), q_points(:, :)
      character(len=*), intent(in) :: representation
      integer, intent(in) :: fourier_axes(:)
      integer :: ia

      if (size(r_vectors, 1) /= 3 .or. size(q_points, 1) /= 3 .or. size(chi_r, 4) /= size(r_vectors, 2) .or. &
          (trim(representation) /= 'bulk' .and. trim(representation) /= 'film' .and. &
           trim(representation) /= 'finite' .and. trim(representation) /= 'local') .or. &
          any(fourier_axes < 0) .or. any(fourier_axes > 3)) then
         error stop 'fourier_transform_realspace_chi0: invalid representation or dimensions'
      end if
      if (trim(representation) == 'bulk' .and. (size(fourier_axes) /= 3 .or. any(fourier_axes == 0))) then
         error stop 'fourier_transform_realspace_chi0: bulk representation requires three Fourier axes'
      end if
      if (trim(representation) == 'film' .and. count(fourier_axes > 0) /= 2) then
         error stop 'fourier_transform_realspace_chi0: film representation requires two Fourier axes'
      end if
      if (trim(representation) == 'film' .and. any([(count(fourier_axes == ia), ia=1, 3)] > 1)) then
         error stop 'fourier_transform_realspace_chi0: film Fourier axes must be distinct'
      end if
   end subroutine validate_fourier_inputs

   real(rp) function resolved_green_eta(options) result(value)
      type(tddft_realspace_chi0_options), intent(in) :: options
      if (options%green_eta > 0.0_rp) then
         value = options%green_eta
      else
         value = 0.5_rp*options%eta
      end if
   end function resolved_green_eta

   logical function realspace_full_tail_mode(mode) result(full_tail)
      character(len=*), intent(in) :: mode

      select case (trim(mode))
      case ('full_tail', 'validation')
         full_tail = .true.
      case ('production', 'truncated')
         full_tail = .false.
      case default
         error stop "native real-space provider: truncation_mode must be 'full_tail'/'validation' or 'production'/'truncated'"
      end select
   end function realspace_full_tail_mode

   integer function count_unique_frequencies(omega) result(count)
      real(rp), intent(in) :: omega(:)
      integer :: iw, jw
      logical :: repeated

      count = 0
      do iw = 1, size(omega)
         repeated = .false.
         do jw = 1, iw-1
            if (abs(omega(iw)-omega(jw)) <= 32.0_rp*epsilon(1.0_rp)* &
                max(1.0_rp, abs(omega(iw)), abs(omega(jw)))) then
               repeated = .true.
               exit
            end if
         end do
         if (.not. repeated) count = count+1
      end do
   end function count_unique_frequencies

   integer function realspace_parallel_size() result(nparallel)
#ifdef USE_MPI
      logical :: initialized
      integer :: mpi_status

      call MPI_INITIALIZED(initialized, mpi_status)
      if (initialized) then
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nparallel, mpi_status)
         return
      end if
#endif
      nparallel = max(1, numprocs)
   end function realspace_parallel_size

   function realspace_pair_vector(lattice_obj, iatom, jatom) result(vector)
      type(lattice), intent(in) :: lattice_obj
      integer, intent(in) :: iatom, jatom
      real(rp) :: vector(3), dcart(3), lattice_inverse(3, 3)

      dcart = (lattice_obj%cr(:, iatom)-lattice_obj%cr(:, jatom))*lattice_obj%alat
      if (lattice_obj%a_cart_inv_ready) then
         vector = matmul(lattice_obj%a_cart_inv, dcart)
      else
         lattice_inverse = inverse_3x3(lattice_obj%a)
         vector = matmul(lattice_inverse, dcart/lattice_obj%alat)
      end if
   end function realspace_pair_vector

   integer function request_q_index(request, iq) result(index)
      type(tddft_chi0_request), intent(in) :: request
      integer, intent(in) :: iq
      if (allocated(request%q_indices)) then
         index = request%q_indices(iq)
      else
         index = iq
      end if
   end function request_q_index

end module tddft_chi0_realspace_mod
