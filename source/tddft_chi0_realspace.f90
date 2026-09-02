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
   use math_mod, only: pi, two_pi, i_unit, inverse_3x3
   use response_vertices_mod, only: response_channel, site_projected_operator
   use tddft_chi0_mod, only: tddft_chi0_result, tddft_chi0_request, tddft_chi0_batch_result, &
      tddft_fermi_occupation, tddft_kB_Ry_per_K, tddft_occupation_kT_floor
   use tddft_backend_mod, only: tddft_realspace_chi0_provider, tddft_backend_capabilities
   use green_mod, only: green
   use lattice_mod, only: lattice
   use mpi_mod, only: ierr, start_atom, end_atom, g2l_map
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
      character(len=16) :: representation = 'bulk'
      character(len=16) :: circular_channel = 'plus_minus'
      integer :: fourier_axes(3) = [1, 2, 3]
      real(rp) :: metric(3, 3) = reshape([1.0_rp, 0.0_rp, 0.0_rp, &
         0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp], [3, 3])
   end type tddft_realspace_chi0_options

   type, public :: tddft_realspace_chi0_diagnostics
      integer :: input_points = 0
      integer :: selected_points = 0
      integer :: omitted_points = 0
      integer :: shell_count = 0
      real(rp) :: rmax = 0.0_rp
      real(rp) :: selected_norm = 0.0_rp
      real(rp) :: omitted_tail_norm = 0.0_rp
      real(rp) :: tail_ratio = 0.0_rp
      real(rp) :: last_shell_norm = 0.0_rp
      real(rp) :: pair_integration_cpu_seconds = 0.0_rp
      real(rp) :: fourier_cpu_seconds = 0.0_rp
      logical :: tail_assessed = .false.
      logical :: converged = .false.
      character(len=48) :: status = 'not assessed'
   end type tddft_realspace_chi0_diagnostics

   type, public :: tddft_realspace_chi0_result
      complex(rp), allocatable :: chi_r(:, :, :, :)
      complex(rp), allocatable :: chi_q(:, :, :, :)
      real(rp), allocatable :: r_vectors(:, :)
      type(tddft_realspace_chi0_diagnostics) :: diagnostics
   end type tddft_realspace_chi0_result

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
      type(tddft_realspace_chi0_options) :: options
      integer :: build_count = 0
      logical :: initialized = .false.
   contains
      procedure :: initialize => initialize_native_realspace_provider
      procedure :: initialize_from_green => initialize_native_from_green
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
      this%options = tddft_realspace_chi0_options()
      if (present(options)) this%options = options
      this%build_count = 0
      this%initialized = .true.
   end subroutine initialize_native_realspace_provider

   subroutine initialize_native_from_green(this, green_obj, lattice_obj, site_orbital_counts, left_channels, &
      right_channels, options)
      class(tddft_native_realspace_gf_provider), intent(inout) :: this
      type(green), intent(in) :: green_obj
      type(lattice), intent(in) :: lattice_obj
      integer, intent(in) :: site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_realspace_chi0_options), intent(in), optional :: options
      real(rp), allocatable :: rvec(:, :)
      integer, allocatable :: pairs(:, :), local_map(:)
      complex(rp), allocatable :: local_ab(:, :, :, :), local_ba(:, :, :, :)
      type(tddft_realspace_chi0_options) :: source_options
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
   end subroutine initialize_native_from_green

   subroutine provider_build_realspace(this, omega, result)
      class(tddft_native_realspace_gf_provider), intent(inout) :: this
      real(rp), intent(in) :: omega(:)
      type(tddft_realspace_chi0_result), intent(out) :: result

      if (.not. this%initialized) error stop 'native real-space provider: provider is not initialized'
      call build_chi0_from_realspace_gf(this%energy_grid, this%g_ab, this%g_ba, this%r_vectors, this%pair_sites, &
         this%site_orbital_counts, this%left_channels, this%right_channels, omega, this%options, result)
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
      call build_chi0_from_realspace_gf(this%energy_grid, this%g_ab, this%g_ba, this%r_vectors, this%pair_sites, &
         this%site_orbital_counts, this%left_channels, this%right_channels, request%omega, this%options, realspace_result, &
         q_points=request%q_points)
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
      real(rp) :: effective_cutoff, selected_norm, omitted_norm, total_norm
      real(rp) :: t_start, t_stop

      call validate_realspace_source(energy_grid, g_ab, g_ba, r_vectors, pair_sites, site_orbital_counts, &
         left_channels, right_channels)
      if (size(omega) < 1 .or. options%eta <= 0.0_rp .or. options%tail_tolerance < 0.0_rp) then
         error stop 'build_chi0_from_realspace_gf: invalid frequency or broadening options'
      end if
      if (trim(options%energy_integration) /= 'direct') then
         error stop 'build_chi0_from_realspace_gf: mixed_contour requires a complex-energy native real-space source; sampled real-axis source is direct-only'
      end if
      nleft = size(left_channels); nright = size(right_channels); nw = size(omega)
      nblock = size(g_ab, 1)
      allocate(radii(size(r_vectors, 2)), keep(size(r_vectors, 2)))
      do ip = 1, size(r_vectors, 2)
         radii(ip) = metric_norm(options%metric, r_vectors(:, ip))
         keep(ip) = radii(ip) <= options%rmax + 16.0_rp*epsilon(1.0_rp)*max(1.0_rp, radii(ip))
      end do
      nkeep = count(keep)
      if (nkeep < 1 .and. size(r_vectors, 2) > 0) error stop 'build_chi0_from_realspace_gf: Rmax removes every real-space pair'
      if (nkeep > 0) then
         effective_cutoff = maxval(pack(radii, keep))
      else
         ! A rank with no locally owned R block is a valid member of an
         ! MPI-distributed native source.  Its zero contribution is reduced
         ! with the other ranks before the common q Fourier batch is used.
         effective_cutoff = 0.0_rp
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
      result%diagnostics%selected_points = nkeep
      result%diagnostics%omitted_points = size(keep)-nkeep
      result%diagnostics%rmax = effective_cutoff
      result%diagnostics%selected_norm = selected_norm
      result%diagnostics%omitted_tail_norm = omitted_norm
      if (total_norm > tiny(1.0_rp)) result%diagnostics%tail_ratio = omitted_norm/total_norm
      result%diagnostics%tail_assessed = result%diagnostics%omitted_points > 0
      result%diagnostics%converged = result%diagnostics%tail_assessed .and. &
         result%diagnostics%tail_ratio <= options%tail_tolerance
      if (.not. result%diagnostics%tail_assessed) then
         result%diagnostics%status = 'all supplied real-space pairs retained'
      else if (result%diagnostics%converged) then
         result%diagnostics%status = 'R-space tail below tolerance'
      else
         result%diagnostics%status = 'R-space tail exceeds tolerance'
      end if
      result%diagnostics%last_shell_norm = last_shell_norm(radii, keep, pair_norms)
      result%diagnostics%shell_count = count_realspace_shells(radii)
   end subroutine build_chi0_from_realspace_gf

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
      real(rp) :: effective_cutoff, selected_norm, omitted_norm, total_norm
      real(rp) :: t_start, t_stop

      call validate_realspace_source(energy_grid, g_ab, g_ba, r_vectors, pair_sites, site_orbital_counts, &
         left_channels, right_channels)
      if (options%tail_tolerance < 0.0_rp .or. options%electronic_temperature < 0.0_rp) then
         error stop 'build_static_chi0_from_realspace_gf: invalid static options'
      end if
      if (trim(options%energy_integration) /= 'direct') then
         error stop 'build_static_chi0_from_realspace_gf: native static source requires direct real-axis G(R,E) data'
      end if
      nleft = size(left_channels); nright = size(right_channels); nblock = size(g_ab, 1)
      allocate(radii(size(r_vectors, 2)), keep(size(r_vectors, 2)))
      do ip = 1, size(r_vectors, 2)
         radii(ip) = metric_norm(options%metric, r_vectors(:, ip))
         keep(ip) = radii(ip) <= options%rmax + 16.0_rp*epsilon(1.0_rp)*max(1.0_rp, radii(ip))
      end do
      nkeep = count(keep)
      if (nkeep < 1 .and. size(r_vectors, 2) > 0) then
         error stop 'build_static_chi0_from_realspace_gf: Rmax removes every real-space pair'
      end if
      if (nkeep > 0) then
         effective_cutoff = maxval(pack(radii, keep))
      else
         effective_cutoff = 0.0_rp
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
      result%diagnostics%selected_points = nkeep
      result%diagnostics%omitted_points = size(keep)-nkeep
      result%diagnostics%rmax = effective_cutoff
      result%diagnostics%selected_norm = selected_norm
      result%diagnostics%omitted_tail_norm = omitted_norm
      if (total_norm > tiny(1.0_rp)) result%diagnostics%tail_ratio = omitted_norm/total_norm
      result%diagnostics%tail_assessed = result%diagnostics%omitted_points > 0
      result%diagnostics%converged = result%diagnostics%tail_assessed .and. &
         result%diagnostics%tail_ratio <= options%tail_tolerance
      if (.not. result%diagnostics%tail_assessed) then
         result%diagnostics%status = 'all supplied real-space pairs retained'
      else if (result%diagnostics%converged) then
         result%diagnostics%status = 'R-space static tail below tolerance'
      else
         result%diagnostics%status = 'R-space static tail exceeds tolerance'
      end if
      result%diagnostics%last_shell_norm = last_shell_norm(radii, keep, pair_norms)
      result%diagnostics%shell_count = count_realspace_shells(radii)
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
      real(rp) :: local_tail, global_tail, local_cutoff, global_cutoff
      real(rp) :: local_last_shell, global_last_shell

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
         local_cutoff = batch%q_response(iq)%metadata%real_space_cutoff
         call MPI_ALLREDUCE(local_cutoff, global_cutoff, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
         batch%q_response(iq)%metadata%real_space_cutoff = global_cutoff
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
      result%metadata%real_space_representation = options%representation
      result%metadata%real_space_fourier_axes = options%fourier_axes
      result%metadata%real_space_omitted_points = diagnostics%omitted_points
      result%metadata%real_space_shell_count = diagnostics%shell_count
      result%metadata%real_space_tail_norm = diagnostics%omitted_tail_norm
      result%metadata%real_space_tail_ratio = diagnostics%tail_ratio
      result%metadata%real_space_last_shell_norm = diagnostics%last_shell_norm
      result%metadata%real_space_pair_integration_cpu_seconds = diagnostics%pair_integration_cpu_seconds
      result%metadata%real_space_fourier_cpu_seconds = diagnostics%fourier_cpu_seconds
      result%metadata%real_space_tail_tolerance = options%tail_tolerance
      result%metadata%real_space_tail_assessed = diagnostics%tail_assessed
      result%metadata%convergence_status = diagnostics%status
      result%metadata%converged = diagnostics%converged
      result%metadata%integration_energy_points = 0
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
