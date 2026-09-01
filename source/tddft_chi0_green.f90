!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Retarded Kohn-Sham response from one-particle Green functions.
!>
!> The real-axis reference implemented here evaluates the exact equilibrium
!> retarded bubble in the form
!>
!> chi_AB^R(w) = -1/(2*pi*i) int dE f(E) Tr[
!>   A G_q^R(E+w) B (G_k^R(E)-G_k^A(E))
!> + A (G_q^R(E)-G_q^A(E)) B G_k^A(E-w) ].
!>
!> Inserting G^R=sum_n |n><n|/(E-e_n+i*gamma) and taking gamma -> 0 gives
!> (f_n-f_m)<n|A|m><m|B|n>/(w+e_n-e_m+i*0+), exactly the convention of
!> tddft_chi0_mod.  With a common one-particle half-width gamma, the bubble
!> linewidth is 2*gamma; use green_eta=eta/2 to reproduce an eigenpair result
!> whose denominator broadening is eta.
!>
!> `green_chi0_provider` knows only the public `green_function_provider`
!> contract.  The eigenpair implementation below is a periodic reference
!> source; a recursion, surface, impurity, or Dyson source can implement the
!> same contract without exposing its storage layout to the response layer.
module tddft_chi0_green_mod
   use precision_mod, only: rp
   use math_mod, only: pi, i_unit, gauss_legendre
   use lehmann_kernel_mod, only: lehmann_kspace_resolvent
   use response_vertices_mod, only: response_channel, site_projected_operator
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, tddft_fermi_occupation, &
      tddft_occupation_kT_floor, build_static_chi_ks_from_eigenpairs_at_q
   implicit none

   private
   integer, parameter, public :: GREEN_BRANCH_K = 1
   integer, parameter, public :: GREEN_BRANCH_KQ = 2

   !> Minimum information a response bubble needs from a one-electron route.
   type, abstract, public :: green_function_provider
   contains
      procedure(green_retarded_matrix), deferred :: get_retarded
      !> Evaluate the same orthogonal Green function at an arbitrary complex
      !> energy.  Mixed contour integration requires this operation; keeping
      !> it separate from get_retarded prevents a contour implementation from
      !> mistaking a real-axis interpolation for analytic continuation.
      procedure :: get_complex => unsupported_green_complex
      procedure(green_spectral_bounds), deferred :: get_spectral_bounds
      procedure :: get_static_chi0 => unsupported_green_static
   end type green_function_provider

   abstract interface
      subroutine green_retarded_matrix(this, branch, ik, energy, eta, green_matrix)
         import :: green_function_provider, rp
         class(green_function_provider), intent(in) :: this
         integer, intent(in) :: branch, ik
         real(rp), intent(in) :: energy, eta
         complex(rp), intent(out) :: green_matrix(:, :)
      end subroutine green_retarded_matrix

      subroutine green_complex_matrix(this, branch, ik, energy, green_matrix)
         import :: green_function_provider, rp
         class(green_function_provider), intent(in) :: this
         integer, intent(in) :: branch, ik
         complex(rp), intent(in) :: energy
         complex(rp), intent(out) :: green_matrix(:, :)
      end subroutine green_complex_matrix

      subroutine green_spectral_bounds(this, energy_min, energy_max)
         import :: green_function_provider, rp
         class(green_function_provider), intent(in) :: this
         real(rp), intent(out) :: energy_min, energy_max
      end subroutine green_spectral_bounds

   end interface

   !> Controls for the validated real-axis path.  An unset energy window is
   !> inferred from the source spectrum plus thermal, frequency, and Lorentzian
   !> tails.  `green_eta=0` uses `eta/2`.
   type, public :: green_chi0_options
      real(rp) :: eta = 0.0_rp
      real(rp) :: green_eta = 0.0_rp
      real(rp) :: fermi_level = 0.0_rp
      real(rp) :: electronic_temperature = 0.0_rp
      real(rp) :: energy_min = huge(1.0_rp)
      real(rp) :: energy_max = -huge(1.0_rp)
      integer :: energy_points = 2001
      character(len=24) :: energy_integration = 'direct'
      integer :: contour_points = 64
      integer :: contour_subdivisions = 8
      integer :: near_fermi_points = 128
      real(rp) :: contour_height = 0.0_rp
      integer :: k_mesh_shape(3) = 0
      real(rp) :: q_direct(3) = 0.0_rp
      character(len=32) :: response_projection = 'site'
   end type green_chi0_options

   !> Canonical KS-response provider.  It is deliberately independent of
   !> reciprocal_green, green_mod, and sigma_provider storage types.
   type, public :: green_chi0_provider
      class(green_function_provider), allocatable :: one_particle
      type(green_chi0_options) :: options
   contains
      procedure :: initialize => initialize_green_chi0_provider
      procedure :: build => build_green_chi0
   end type green_chi0_provider

   !> Periodic reference one-particle provider.  It is also a useful adapter
   !> for validating a GF bubble against the eigenpair Lehmann backend.
   type, extends(green_function_provider), public :: eigenpair_green_function_provider
      real(rp), allocatable :: eigenvalues_k(:, :), eigenvalues_kq(:, :)
      complex(rp), allocatable :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :)
   contains
      procedure :: initialize => initialize_eigenpair_green_provider
      procedure :: get_retarded => eigenpair_retarded_green
      procedure :: get_complex => eigenpair_complex_green
      procedure :: get_spectral_bounds => eigenpair_spectral_bounds
      procedure :: get_static_chi0 => eigenpair_static_chi0
   end type eigenpair_green_function_provider

   public :: build_chi_ks_from_green_functions
   public :: build_static_chi_ks_from_green_functions
   public :: build_static_four_component_chi_ks_from_green_functions
   public :: build_four_component_chi_ks_from_green_functions

contains

   subroutine initialize_green_chi0_provider(this, one_particle, options)
      class(green_chi0_provider), intent(inout) :: this
      class(green_function_provider), intent(in) :: one_particle
      type(green_chi0_options), intent(in) :: options

      if (allocated(this%one_particle)) deallocate(this%one_particle)
      allocate(this%one_particle, source=one_particle)
      this%options = options
   end subroutine initialize_green_chi0_provider

   subroutine initialize_eigenpair_green_provider(this, eigenvalues_k, eigenvectors_k, eigenvalues_kq, eigenvectors_kq)
      class(eigenpair_green_function_provider), intent(inout) :: this
      real(rp), intent(in) :: eigenvalues_k(:, :), eigenvalues_kq(:, :)
      complex(rp), intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :)
      integer :: nmat, nk

      nmat = size(eigenvalues_k, 1)
      nk = size(eigenvalues_k, 2)
      if (nmat < 1 .or. any(shape(eigenvalues_kq) /= shape(eigenvalues_k)) .or. &
          any(shape(eigenvectors_k) /= [nmat, nmat, nk]) .or. any(shape(eigenvectors_kq) /= shape(eigenvectors_k))) then
         error stop 'initialize_eigenpair_green_provider: incompatible eigenpair arrays'
      end if
      this%eigenvalues_k = eigenvalues_k
      this%eigenvalues_kq = eigenvalues_kq
      this%eigenvectors_k = eigenvectors_k
      this%eigenvectors_kq = eigenvectors_kq
   end subroutine initialize_eigenpair_green_provider

   subroutine eigenpair_retarded_green(this, branch, ik, energy, eta, green_matrix)
      class(eigenpair_green_function_provider), intent(in) :: this
      integer, intent(in) :: branch, ik
      real(rp), intent(in) :: energy, eta
      complex(rp), intent(out) :: green_matrix(:, :)
      integer :: nmat
      complex(rp) :: z

      if (eta <= 0.0_rp) error stop 'eigenpair_retarded_green: eta must be positive'
      nmat = size(this%eigenvalues_k, 1)
      if (ik < 1 .or. ik > size(this%eigenvalues_k, 2) .or. any(shape(green_matrix) /= [nmat, nmat])) then
         error stop 'eigenpair_retarded_green: invalid Green-function request'
      end if
      z = cmplx(energy, eta, rp)
      call this%get_complex(branch, ik, z, green_matrix)
   end subroutine eigenpair_retarded_green

   subroutine eigenpair_complex_green(this, branch, ik, energy, green_matrix)
      class(eigenpair_green_function_provider), intent(in) :: this
      integer, intent(in) :: branch, ik
      complex(rp), intent(in) :: energy
      complex(rp), intent(out) :: green_matrix(:, :)
      integer :: nmat

      nmat = size(this%eigenvalues_k, 1)
      if (ik < 1 .or. ik > size(this%eigenvalues_k, 2) .or. any(shape(green_matrix) /= [nmat, nmat])) then
         error stop 'eigenpair_complex_green: invalid Green-function request'
      end if
      select case (branch)
      case (GREEN_BRANCH_K)
         call lehmann_kspace_resolvent(this%eigenvalues_k(:, ik), this%eigenvectors_k(:, :, ik), &
            energy, green_matrix)
      case (GREEN_BRANCH_KQ)
         call lehmann_kspace_resolvent(this%eigenvalues_kq(:, ik), this%eigenvectors_kq(:, :, ik), &
            energy, green_matrix)
      case default
         error stop 'eigenpair_complex_green: unknown Green-function branch'
      end select
   end subroutine eigenpair_complex_green

   subroutine eigenpair_spectral_bounds(this, energy_min, energy_max)
      class(eigenpair_green_function_provider), intent(in) :: this
      real(rp), intent(out) :: energy_min, energy_max

      if (.not. allocated(this%eigenvalues_k)) error stop 'eigenpair_spectral_bounds: provider is not initialized'
      energy_min = min(minval(this%eigenvalues_k), minval(this%eigenvalues_kq))
      energy_max = max(maxval(this%eigenvalues_k), maxval(this%eigenvalues_kq))
   end subroutine eigenpair_spectral_bounds

   !> The static limit is an actual zero-frequency divided difference, not a
   !> dynamic sample at finite eta.  The Lehmann GF provider owns the same
   !> spectral endpoint data as the dynamic GF source, so it can provide this
   !> exact static operation while retaining the common response vertices.
   subroutine eigenpair_static_chi0(this, k_weights, site_orbital_counts, left_channels, right_channels, options, result)
      class(eigenpair_green_function_provider), intent(in) :: this
      real(rp), intent(in) :: k_weights(:)
      integer, intent(in) :: site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_chi0_options), intent(in) :: options
      type(tddft_chi0_result), intent(out) :: result

      if (.not. allocated(this%eigenvalues_k)) then
         error stop 'eigenpair_static_chi0: provider is not initialized'
      end if
      call build_static_chi_ks_from_eigenpairs_at_q(k_weights, this%eigenvalues_k, this%eigenvectors_k, &
         this%eigenvalues_kq, this%eigenvectors_kq, site_orbital_counts, left_channels, right_channels, options, result)
   end subroutine eigenpair_static_chi0

   subroutine unsupported_green_static(this, k_weights, site_orbital_counts, left_channels, right_channels, options, result)
      class(green_function_provider), intent(in) :: this
      real(rp), intent(in) :: k_weights(:)
      integer, intent(in) :: site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_chi0_options), intent(in) :: options
      type(tddft_chi0_result), intent(out) :: result

      error stop 'green_function_provider: exact static chi0 is unavailable for this Green-function source'
   end subroutine unsupported_green_static

   subroutine unsupported_green_complex(this, branch, ik, energy, green_matrix)
      class(green_function_provider), intent(in) :: this
      integer, intent(in) :: branch, ik
      complex(rp), intent(in) :: energy
      complex(rp), intent(out) :: green_matrix(:, :)

      error stop 'green_function_provider: mixed contour requires a complex-energy Green-function source'
   end subroutine unsupported_green_complex

   !> Build chi_KS in the same `(left,right,omega)` basis and result type used
   !> by the eigenpair backend.
   subroutine build_chi_ks_from_green_functions(one_particle, k_weights, site_orbital_counts, left_channels, right_channels, &
      omega, options, result)
      class(green_function_provider), intent(in) :: one_particle
      real(rp), intent(in) :: k_weights(:), omega(:)
      integer, intent(in) :: site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(green_chi0_options), intent(in) :: options
      type(tddft_chi0_result), intent(out) :: result
      type(green_chi0_provider) :: provider

      call provider%initialize(one_particle, options)
      call provider%build(k_weights, site_orbital_counts, left_channels, right_channels, omega, result)
   end subroutine build_chi_ks_from_green_functions

   !> Build the exact static limit supported by the one-particle source.  A
   !> source without a spectral/static implementation fails explicitly rather
   !> than silently replacing the static response by a finite-eta dynamic one.
   subroutine build_static_chi_ks_from_green_functions(one_particle, k_weights, site_orbital_counts, left_channels, &
      right_channels, options, result)
      class(green_function_provider), intent(in) :: one_particle
      real(rp), intent(in) :: k_weights(:)
      integer, intent(in) :: site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(green_chi0_options), intent(in) :: options
      type(tddft_chi0_result), intent(out) :: result
      type(tddft_chi0_options) :: static_options

      static_options%eta = 0.0_rp
      static_options%fermi_level = options%fermi_level
      static_options%electronic_temperature = options%electronic_temperature
      static_options%k_mesh_shape = options%k_mesh_shape
      static_options%q_direct = options%q_direct
      static_options%response_projection = options%response_projection
      call one_particle%get_static_chi0(k_weights, site_orbital_counts, left_channels, right_channels, static_options, result)
      result%metadata%backend = 'green_static'
      result%metadata%canonical_backend = 'kspace_lehmann'
      result%metadata%implementation = 'K-space Lehmann GF exact static limit'
      result%metadata%energy_integration = 'exact provider static limit'
      result%metadata%endpoint_provenance = 'K and K+q Lehmann GF endpoint source; exact static limit'
      result%metadata%q_direct = options%q_direct
      result%metadata%response_projection = options%response_projection
      result%metadata%green_eta = 0.0_rp
      result%metadata%integration_energy_points = 0
      result%metadata%omega_min = 0.0_rp
      result%metadata%omega_max = 0.0_rp
      result%metadata%omega_points = 1
      result%metadata%q_batch_size = 1
      result%metadata%omega_batch_size = 1
   end subroutine build_static_chi_ks_from_green_functions

   subroutine build_green_chi0(this, k_weights, site_orbital_counts, left_channels, right_channels, omega, result)
      class(green_chi0_provider), intent(in) :: this
      real(rp), intent(in) :: k_weights(:), omega(:)
      integer, intent(in) :: site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_chi0_result), intent(out) :: result
      complex(rp), allocatable :: left_ops(:, :, :), right_ops(:, :, :)
      complex(rp), allocatable :: gr_k(:, :), ga_k(:, :), gr_q(:, :), ga_q(:, :), gr_shift(:, :), ga_shift(:, :)
      complex(rp), allocatable :: trace_work1(:, :), trace_work2(:, :), trace_work3(:, :)
      complex(rp) :: bubble_first, bubble_second
      real(rp) :: energy_min, energy_max, denergy, energy, weight_sum, green_eta, thermal_tail, omega_tail
      real(rp) :: quadrature_weight, prefactor, t_start, t_stop
      integer :: nmat, nk, nleft, nright, nw, ne, ik, ie, iw, ileft, iright

      if (.not. allocated(this%one_particle)) error stop 'green_chi0_provider%build: one-particle provider is absent'
      nk = size(k_weights); nleft = size(left_channels); nright = size(right_channels); nw = size(omega)
      nmat = 2*sum(site_orbital_counts); ne = this%options%energy_points
      if (nk < 1 .or. nleft < 1 .or. nright < 1 .or. nw < 1 .or. nmat < 1 .or. ne < 3) then
         error stop 'green_chi0_provider%build: empty response input or fewer than three energy points'
      end if
      if (any(k_weights < 0.0_rp) .or. sum(k_weights) <= tiny(1.0_rp) .or. this%options%eta <= 0.0_rp .or. &
          this%options%electronic_temperature < 0.0_rp) then
         error stop 'green_chi0_provider%build: invalid weights, broadening, or temperature'
      end if
      green_eta = this%options%green_eta
      if (green_eta == 0.0_rp) green_eta = 0.5_rp*this%options%eta
      if (green_eta <= 0.0_rp) error stop 'green_chi0_provider%build: green_eta must be positive'
      if (trim(this%options%energy_integration) /= 'direct' .and. &
          trim(this%options%energy_integration) /= 'mixed_contour') then
         error stop 'green_chi0_provider%build: energy_integration must be direct or mixed_contour'
      end if
      if (this%options%contour_points < 8 .or. this%options%contour_subdivisions < 1 .or. &
          this%options%near_fermi_points < 8 .or. &
          this%options%contour_height < 0.0_rp) then
         error stop 'green_chi0_provider%build: invalid mixed-contour quadrature settings'
      end if
      call resolve_green_energy_window(this%one_particle, this%options, omega, energy_min, energy_max, thermal_tail, omega_tail)
      if (energy_max <= energy_min) error stop 'green_chi0_provider%build: invalid integration window'
      if (trim(this%options%energy_integration) == 'mixed_contour') then
         call build_green_chi0_mixed_contour(this, k_weights, site_orbital_counts, left_channels, right_channels, omega, &
            energy_min, energy_max, green_eta, result)
         return
      end if
      denergy = (energy_max-energy_min)/real(ne-1, rp)
      weight_sum = sum(k_weights)

      allocate(left_ops(nmat, nmat, nleft), right_ops(nmat, nmat, nright))
      do ileft = 1, nleft
         left_ops(:, :, ileft) = site_projected_operator(left_channels(ileft), site_orbital_counts)
      end do
      do iright = 1, nright
         right_ops(:, :, iright) = site_projected_operator(right_channels(iright), site_orbital_counts)
      end do
      allocate(result%chi(nleft, nright, nw), result%re_chi(nleft, nright, nw), result%im_chi(nleft, nright, nw))
      result%chi = cmplx(0.0_rp, 0.0_rp, rp)
      allocate(gr_k(nmat, nmat), ga_k(nmat, nmat), gr_q(nmat, nmat), ga_q(nmat, nmat), &
         gr_shift(nmat, nmat), ga_shift(nmat, nmat), trace_work1(nmat, nmat), trace_work2(nmat, nmat), &
         trace_work3(nmat, nmat))

      call cpu_time(t_start)
      do ik = 1, nk
         prefactor = k_weights(ik)/weight_sum
         do ie = 1, ne
            energy = energy_min + real(ie-1, rp)*denergy
            quadrature_weight = denergy
            if (ie == 1 .or. ie == ne) quadrature_weight = 0.5_rp*quadrature_weight
            call this%one_particle%get_retarded(GREEN_BRANCH_K, ik, energy, green_eta, gr_k)
            ga_k = transpose(conjg(gr_k))
            call this%one_particle%get_retarded(GREEN_BRANCH_KQ, ik, energy, green_eta, gr_q)
            ga_q = transpose(conjg(gr_q))
            do iw = 1, nw
               call this%one_particle%get_retarded(GREEN_BRANCH_KQ, ik, energy+omega(iw), green_eta, gr_shift)
               call this%one_particle%get_retarded(GREEN_BRANCH_K, ik, energy-omega(iw), green_eta, ga_shift)
               ga_shift = transpose(conjg(ga_shift))
               do iright = 1, nright
                  do ileft = 1, nleft
                     bubble_first = trace_four(left_ops(:, :, ileft), gr_shift, right_ops(:, :, iright), gr_k-ga_k, trace_work1, &
                        trace_work2, trace_work3)
                     bubble_second = trace_four(left_ops(:, :, ileft), gr_q-ga_q, right_ops(:, :, iright), ga_shift, trace_work1, &
                        trace_work2, trace_work3)
                     result%chi(ileft, iright, iw) = result%chi(ileft, iright, iw) - prefactor*quadrature_weight* &
                        tddft_fermi_occupation(energy, this%options%fermi_level, this%options%electronic_temperature)* &
                        (bubble_first+bubble_second)/(2.0_rp*pi*i_unit)
                  end do
               end do
            end do
         end do
      end do
      call cpu_time(t_stop)
      result%re_chi = real(result%chi, rp)
      result%im_chi = aimag(result%chi)
      call build_spectral_products(left_channels, right_channels, result)
      result%metadata%backend = 'green'
      result%metadata%canonical_backend = 'kspace_lehmann'
      result%metadata%implementation = 'K-space Lehmann GF bubble'
      result%metadata%energy_integration = 'real-axis trapezoid'
      result%metadata%endpoint_provenance = 'K/K+q branches; exact folded endpoint supplied by caller'
      result%metadata%eta = 2.0_rp*green_eta
      result%metadata%green_eta = green_eta
      result%metadata%integration_energy_min = energy_min
      result%metadata%integration_energy_max = energy_max
      result%metadata%integration_energy_points = ne
      result%metadata%fermi_level = this%options%fermi_level
      result%metadata%electronic_temperature = this%options%electronic_temperature
      result%metadata%electronic_kT = max(this%options%electronic_temperature*6.3336814e-6_rp, tddft_occupation_kT_floor)
      result%metadata%k_weight_sum = weight_sum
      result%metadata%k_mesh_shape = this%options%k_mesh_shape
      result%metadata%nk = nk
      result%metadata%available_band_count = nmat
      result%metadata%band_first = 1
      result%metadata%band_last = nmat
      result%metadata%green_energy_integration_cpu_seconds = t_stop-t_start
      result%metadata%green_function_evaluations = nk*(2*ne+2*ne*nw)
      result%metadata%contour_points = 0
      result%metadata%contour_subdivisions = 0
      result%metadata%near_fermi_points = 0
      result%metadata%contour_height = 0.0_rp
      result%metadata%contour_max_imaginary_energy = green_eta
      result%metadata%contour_gf_evaluations = 0
      result%metadata%contour_deformation = .false.
      result%metadata%convergence_status = 'not assessed by backend'
      result%metadata%converged = .false.
      result%metadata%q_batch_size = 1
      result%metadata%omega_batch_size = nw
      result%metadata%q_direct = this%options%q_direct
      result%metadata%response_projection = this%options%response_projection
      result%metadata%omega_min = minval(omega)
      result%metadata%omega_max = maxval(omega)
      result%metadata%omega_points = nw
   end subroutine build_green_chi0

   subroutine build_four_component_chi_ks_from_green_functions(one_particle, k_weights, site_orbital_counts, omega, options, result)
      class(green_function_provider), intent(in) :: one_particle
      real(rp), intent(in) :: k_weights(:), omega(:)
      integer, intent(in) :: site_orbital_counts(:)
      type(green_chi0_options), intent(in) :: options
      type(tddft_chi0_result), intent(out) :: result
      type(response_channel), allocatable :: channels(:)
      integer :: isite, component

      allocate(channels(4*size(site_orbital_counts)))
      do isite = 1, size(site_orbital_counts)
         do component = 0, 3
            channels(4*(isite-1)+component+1) = response_channel(isite, component)
         end do
      end do
      call build_chi_ks_from_green_functions(one_particle, k_weights, site_orbital_counts, channels, channels, omega, options, result)
   end subroutine build_four_component_chi_ks_from_green_functions

   subroutine build_static_four_component_chi_ks_from_green_functions(one_particle, k_weights, site_orbital_counts, options, result)
      class(green_function_provider), intent(in) :: one_particle
      real(rp), intent(in) :: k_weights(:)
      integer, intent(in) :: site_orbital_counts(:)
      type(green_chi0_options), intent(in) :: options
      type(tddft_chi0_result), intent(out) :: result
      type(response_channel), allocatable :: channels(:)
      integer :: isite, component

      allocate(channels(4*size(site_orbital_counts)))
      do isite = 1, size(site_orbital_counts)
         do component = 0, 3
            channels(4*(isite-1)+component+1) = response_channel(isite, component)
         end do
      end do
      call build_static_chi_ks_from_green_functions(one_particle, k_weights, site_orbital_counts, channels, channels, &
         options, result)
   end subroutine build_static_four_component_chi_ks_from_green_functions

   !> Build the mixed Lounis/Mills response.  The two same-half-plane terms
   !> are deformed away from the real axis, while the retarded/advanced term
   !> is retained on the short interval where f(E)-f(E+omega) is nonzero.
   !> The algebra is
   !>
   !>   chi0 = -1/(2*pi*i) [ I_RR - I_RA - I_AA ],
   !>
   !> with I_RA = integral (f(E)-f(E+omega)) Tr[A Gq^R(E+w) B Gk^A(E)].
   !> No analytic continuation of chi0 is used.  `get_complex` is a required
   !> source operation here so a sampled real-axis GF cannot be mislabelled as
   !> a contour evaluator.
   subroutine build_green_chi0_mixed_contour(this, k_weights, site_orbital_counts, left_channels, right_channels, omega, &
      energy_min, energy_max, green_eta, result)
      class(green_chi0_provider), intent(in) :: this
      real(rp), intent(in) :: k_weights(:), omega(:), energy_min, energy_max, green_eta
      integer, intent(in) :: site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_chi0_result), intent(out) :: result
      complex(rp), allocatable :: left_ops(:, :, :), right_ops(:, :, :)
      complex(rp), allocatable :: rr_integral(:, :), ra_integral(:, :), aa_integral(:, :)
      complex(rp), allocatable :: gr_k(:, :), gr_q(:, :), gr_shift(:, :), ga_shift(:, :)
      complex(rp), allocatable :: trace_work1(:, :), trace_work2(:, :), trace_work3(:, :)
      real(rp) :: weight_sum, prefactor, contour_height, same_min, same_max, t_start, t_stop
      real(rp) :: max_contour_imaginary_energy
      integer :: nmat, nk, nleft, nright, nw, ik, iw, ileft, iright
      integer :: contour_gf_evaluations, near_gf_evaluations

      nk = size(k_weights); nleft = size(left_channels); nright = size(right_channels); nw = size(omega)
      nmat = 2*sum(site_orbital_counts); weight_sum = sum(k_weights)
      call mixed_contour_height(this%options, energy_max-energy_min, green_eta, contour_height)
      call mixed_contour_same_interval(this%options, energy_min, energy_max, same_min, same_max)

      allocate(left_ops(nmat, nmat, nleft), right_ops(nmat, nmat, nright))
      do ileft = 1, nleft
         left_ops(:, :, ileft) = site_projected_operator(left_channels(ileft), site_orbital_counts)
      end do
      do iright = 1, nright
         right_ops(:, :, iright) = site_projected_operator(right_channels(iright), site_orbital_counts)
      end do
      allocate(result%chi(nleft, nright, nw), result%re_chi(nleft, nright, nw), result%im_chi(nleft, nright, nw))
      result%chi = cmplx(0.0_rp, 0.0_rp, rp)
      allocate(rr_integral(nleft, nright), ra_integral(nleft, nright), aa_integral(nleft, nright), &
         gr_k(nmat, nmat), gr_q(nmat, nmat), gr_shift(nmat, nmat), ga_shift(nmat, nmat), &
         trace_work1(nmat, nmat), trace_work2(nmat, nmat), trace_work3(nmat, nmat))

      contour_gf_evaluations = 0; near_gf_evaluations = 0; max_contour_imaginary_energy = green_eta
      call cpu_time(t_start)
      do ik = 1, nk
         prefactor = k_weights(ik)/weight_sum
         do iw = 1, nw
            rr_integral = cmplx(0.0_rp, 0.0_rp, rp)
            aa_integral = cmplx(0.0_rp, 0.0_rp, rp)
            ra_integral = cmplx(0.0_rp, 0.0_rp, rp)
            if (same_max > same_min) then
               call integrate_same_contour_segment(this%one_particle, ik, omega(iw), green_eta, this%options, 1, &
                  cmplx(same_min, 0.0_rp, rp), cmplx(same_min, contour_height, rp), left_ops, right_ops, rr_integral, &
                  aa_integral, gr_k, gr_q, gr_shift, ga_shift, trace_work1, trace_work2, trace_work3, &
                  contour_gf_evaluations, max_contour_imaginary_energy)
               call integrate_same_contour_segment(this%one_particle, ik, omega(iw), green_eta, this%options, 1, &
                  cmplx(same_min, contour_height, rp), cmplx(same_max, contour_height, rp), left_ops, right_ops, rr_integral, &
                  aa_integral, gr_k, gr_q, gr_shift, ga_shift, trace_work1, trace_work2, trace_work3, &
                  contour_gf_evaluations, max_contour_imaginary_energy)
               call integrate_same_contour_segment(this%one_particle, ik, omega(iw), green_eta, this%options, 1, &
                  cmplx(same_max, contour_height, rp), cmplx(same_max, 0.0_rp, rp), left_ops, right_ops, rr_integral, &
                  aa_integral, gr_k, gr_q, gr_shift, ga_shift, trace_work1, trace_work2, trace_work3, &
                  contour_gf_evaluations, max_contour_imaginary_energy)
               call integrate_same_contour_segment(this%one_particle, ik, omega(iw), green_eta, this%options, -1, &
                  cmplx(same_min, 0.0_rp, rp), cmplx(same_min, -contour_height, rp), left_ops, right_ops, rr_integral, &
                  aa_integral, gr_k, gr_q, gr_shift, ga_shift, trace_work1, trace_work2, trace_work3, &
                  contour_gf_evaluations, max_contour_imaginary_energy)
               call integrate_same_contour_segment(this%one_particle, ik, omega(iw), green_eta, this%options, -1, &
                  cmplx(same_min, -contour_height, rp), cmplx(same_max, -contour_height, rp), left_ops, right_ops, rr_integral, &
                  aa_integral, gr_k, gr_q, gr_shift, ga_shift, trace_work1, trace_work2, trace_work3, &
                  contour_gf_evaluations, max_contour_imaginary_energy)
               call integrate_same_contour_segment(this%one_particle, ik, omega(iw), green_eta, this%options, -1, &
                  cmplx(same_max, -contour_height, rp), cmplx(same_max, 0.0_rp, rp), left_ops, right_ops, rr_integral, &
                  aa_integral, gr_k, gr_q, gr_shift, ga_shift, trace_work1, trace_work2, trace_work3, &
                  contour_gf_evaluations, max_contour_imaginary_energy)
            end if
            call integrate_cross_interval(this%one_particle, ik, omega(iw), green_eta, this%options, energy_min, energy_max, &
               left_ops, right_ops, ra_integral, gr_q, ga_shift, trace_work1, trace_work2, trace_work3, near_gf_evaluations)
            do iright = 1, nright
               do ileft = 1, nleft
                  result%chi(ileft, iright, iw) = result%chi(ileft, iright, iw) - prefactor* &
                     (rr_integral(ileft, iright)-ra_integral(ileft, iright)-aa_integral(ileft, iright))/ &
                     (2.0_rp*pi*i_unit)
               end do
            end do
         end do
      end do
      call cpu_time(t_stop)
      result%re_chi = real(result%chi, rp); result%im_chi = aimag(result%chi)
      call build_spectral_products(left_channels, right_channels, result)
      result%metadata%backend = 'green'
      result%metadata%canonical_backend = 'kspace_lehmann'
      result%metadata%implementation = 'K-space Lehmann GF mixed Lounis contour'
      result%metadata%energy_integration = 'mixed contour: RR/AA deformation plus near-Fermi RA'
      result%metadata%endpoint_provenance = 'K/K+q branches; complex-energy Lehmann source'
      result%metadata%eta = 2.0_rp*green_eta; result%metadata%green_eta = green_eta
      result%metadata%integration_energy_min = energy_min; result%metadata%integration_energy_max = energy_max
      result%metadata%integration_energy_points = 2*(2+this%options%contour_subdivisions)*this%options%contour_points + &
         3*this%options%near_fermi_points
      result%metadata%fermi_level = this%options%fermi_level
      result%metadata%electronic_temperature = this%options%electronic_temperature
      result%metadata%electronic_kT = max(this%options%electronic_temperature*6.3336814e-6_rp, tddft_occupation_kT_floor)
      result%metadata%k_weight_sum = weight_sum; result%metadata%k_mesh_shape = this%options%k_mesh_shape
      result%metadata%nk = nk; result%metadata%available_band_count = nmat
      result%metadata%band_first = 1; result%metadata%band_last = nmat
      result%metadata%green_energy_integration_cpu_seconds = t_stop-t_start
      result%metadata%contour_points = this%options%contour_points
      result%metadata%contour_subdivisions = this%options%contour_subdivisions
      result%metadata%near_fermi_points = this%options%near_fermi_points
      result%metadata%contour_height = contour_height
      result%metadata%contour_max_imaginary_energy = max_contour_imaginary_energy
      result%metadata%contour_gf_evaluations = contour_gf_evaluations + near_gf_evaluations
      result%metadata%green_function_evaluations = contour_gf_evaluations + near_gf_evaluations
      result%metadata%contour_deformation = .true.
      result%metadata%convergence_status = 'contour-node convergence not assessed by backend'
      result%metadata%converged = .false.
      result%metadata%q_batch_size = 1; result%metadata%omega_batch_size = nw
      result%metadata%q_direct = this%options%q_direct
      result%metadata%response_projection = this%options%response_projection
      result%metadata%omega_min = minval(omega); result%metadata%omega_max = maxval(omega); result%metadata%omega_points = nw
   end subroutine build_green_chi0_mixed_contour

   subroutine integrate_same_contour_segment(one_particle, ik, omega, green_eta, options, plane, z_start, z_end, left_ops, right_ops, &
      rr_integral, aa_integral, gr_k, gr_q, gr_shift, ga_shift, trace_work1, trace_work2, trace_work3, gf_evaluations, max_imag)
      class(green_function_provider), intent(in) :: one_particle
      integer, intent(in) :: ik
      real(rp), intent(in) :: omega, green_eta
      integer, intent(in) :: plane
      type(green_chi0_options), intent(in) :: options
      complex(rp), intent(in) :: z_start, z_end
      complex(rp), intent(in) :: left_ops(:, :, :), right_ops(:, :, :)
      complex(rp), intent(inout) :: rr_integral(:, :), aa_integral(:, :)
      complex(rp), intent(inout) :: gr_k(:, :), gr_q(:, :), gr_shift(:, :), ga_shift(:, :)
      complex(rp), intent(inout) :: trace_work1(:, :), trace_work2(:, :), trace_work3(:, :)
      integer, intent(inout) :: gf_evaluations
      real(rp), intent(inout) :: max_imag
      real(rp) :: nodes(options%contour_points), weights(options%contour_points), parameter
      complex(rp) :: segment_start, segment_end
      complex(rp) :: z, fermi, rr_value, jacobian
      integer :: inode, ileft, iright, nsub, isub

      nsub = 1
      if (abs(real(z_end-z_start, rp)) > abs(aimag(z_end-z_start))) nsub = options%contour_subdivisions
      do isub = 1, nsub
         segment_start = z_start + (z_end-z_start)*real(isub-1, rp)/real(nsub, rp)
         segment_end = z_start + (z_end-z_start)*real(isub, rp)/real(nsub, rp)
         call gauss_legendre(options%contour_points, 0.0_rp, 1.0_rp, nodes, weights)
         do inode = 1, options%contour_points
            parameter = nodes(inode)
            z = segment_start + parameter*(segment_end-segment_start)
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
            if (plane > 0) then
               call one_particle%get_complex(GREEN_BRANCH_KQ, ik, z+cmplx(omega, green_eta, rp), gr_shift)
               call one_particle%get_complex(GREEN_BRANCH_K, ik, z+cmplx(0.0_rp, green_eta, rp), gr_k)
               gf_evaluations = gf_evaluations + 2
               max_imag = max(max_imag, abs(aimag(z+cmplx(omega, green_eta, rp))), &
                  abs(aimag(z+cmplx(0.0_rp, green_eta, rp))))
               do iright = 1, size(right_ops, 3)
                  do ileft = 1, size(left_ops, 3)
                     rr_value = trace_four(left_ops(:, :, ileft), gr_shift, right_ops(:, :, iright), gr_k, trace_work1, &
                        trace_work2, trace_work3)
                     rr_integral(ileft, iright) = rr_integral(ileft, iright) + jacobian*fermi*rr_value
                  end do
               end do
            else
               call one_particle%get_complex(GREEN_BRANCH_KQ, ik, z-cmplx(0.0_rp, green_eta, rp), gr_q)
               call one_particle%get_complex(GREEN_BRANCH_K, ik, z-cmplx(omega, green_eta, rp), ga_shift)
               gf_evaluations = gf_evaluations + 2
               max_imag = max(max_imag, abs(aimag(z-cmplx(0.0_rp, green_eta, rp))), &
                  abs(aimag(z-cmplx(omega, green_eta, rp))))
               do iright = 1, size(right_ops, 3)
                  do ileft = 1, size(left_ops, 3)
                     rr_value = trace_four(left_ops(:, :, ileft), gr_q, right_ops(:, :, iright), ga_shift, trace_work1, &
                        trace_work2, trace_work3)
                     aa_integral(ileft, iright) = aa_integral(ileft, iright) + jacobian*fermi*rr_value
                  end do
               end do
            end if
         end do
      end do
   end subroutine integrate_same_contour_segment

   subroutine integrate_cross_interval(one_particle, ik, omega, green_eta, options, energy_min, energy_max, left_ops, right_ops, &
      ra_integral, gr_q, ga_shift, trace_work1, trace_work2, trace_work3, gf_evaluations)
      class(green_function_provider), intent(in) :: one_particle
      integer, intent(in) :: ik
      real(rp), intent(in) :: omega, green_eta, energy_min, energy_max
      type(green_chi0_options), intent(in) :: options
      complex(rp), intent(in) :: left_ops(:, :, :), right_ops(:, :, :)
      complex(rp), intent(inout) :: ra_integral(:, :), gr_q(:, :), ga_shift(:, :)
      complex(rp), intent(inout) :: trace_work1(:, :), trace_work2(:, :), trace_work3(:, :)
      integer, intent(inout) :: gf_evaluations
      real(rp) :: lower, upper, thermal_window

      if (abs(omega) <= tiny(1.0_rp)) return
      thermal_window = 0.0_rp
      if (options%electronic_temperature > 0.0_rp) then
         thermal_window = 50.0_rp*max(options%electronic_temperature*6.3336814e-6_rp, tddft_occupation_kT_floor)
      end if
      ! Keep the exact finite-window algebra.  After shifting the second
      ! mixed RA term by omega, its interval is [energy_min-omega,
      ! energy_max-omega], not [energy_min,energy_max].  The central interval
      ! is the near-Fermi Lounis piece; the two width-|omega| edge intervals
      ! are retained as inexpensive boundary corrections.  This avoids a
      ! hidden dependence on an arbitrarily wide artificial energy window.
      if (omega > 0.0_rp) then
         lower = max(energy_min, min(options%fermi_level, options%fermi_level-omega)-thermal_window)
         upper = min(energy_max-omega, max(options%fermi_level, options%fermi_level-omega)+thermal_window)
         call integrate_cross_subinterval(one_particle, ik, omega, green_eta, options, energy_min-omega, energy_min, -1, &
            left_ops, right_ops, ra_integral, gr_q, ga_shift, trace_work1, trace_work2, trace_work3, gf_evaluations)
         call integrate_cross_subinterval(one_particle, ik, omega, green_eta, options, lower, upper, 0, left_ops, right_ops, &
            ra_integral, gr_q, ga_shift, trace_work1, trace_work2, trace_work3, gf_evaluations)
         call integrate_cross_subinterval(one_particle, ik, omega, green_eta, options, energy_max-omega, energy_max, 1, &
            left_ops, right_ops, ra_integral, gr_q, ga_shift, trace_work1, trace_work2, trace_work3, gf_evaluations)
      else
         lower = max(energy_min-omega, min(options%fermi_level, options%fermi_level-omega)-thermal_window)
         upper = min(energy_max, max(options%fermi_level, options%fermi_level-omega)+thermal_window)
         call integrate_cross_subinterval(one_particle, ik, omega, green_eta, options, energy_min, energy_min-omega, 1, &
            left_ops, right_ops, ra_integral, gr_q, ga_shift, trace_work1, trace_work2, trace_work3, gf_evaluations)
         call integrate_cross_subinterval(one_particle, ik, omega, green_eta, options, lower, upper, 0, left_ops, right_ops, &
            ra_integral, gr_q, ga_shift, trace_work1, trace_work2, trace_work3, gf_evaluations)
         call integrate_cross_subinterval(one_particle, ik, omega, green_eta, options, energy_max, energy_max-omega, -1, &
            left_ops, right_ops, ra_integral, gr_q, ga_shift, trace_work1, trace_work2, trace_work3, gf_evaluations)
      end if
   end subroutine integrate_cross_interval

   subroutine integrate_cross_subinterval(one_particle, ik, omega, green_eta, options, lower, upper, coefficient, left_ops, &
      right_ops, ra_integral, gr_q, ga_shift, trace_work1, trace_work2, trace_work3, gf_evaluations)
      class(green_function_provider), intent(in) :: one_particle
      integer, intent(in) :: ik, coefficient
      real(rp), intent(in) :: omega, green_eta, lower, upper
      type(green_chi0_options), intent(in) :: options
      complex(rp), intent(in) :: left_ops(:, :, :), right_ops(:, :, :)
      complex(rp), intent(inout) :: ra_integral(:, :), gr_q(:, :), ga_shift(:, :)
      complex(rp), intent(inout) :: trace_work1(:, :), trace_work2(:, :), trace_work3(:, :)
      integer, intent(inout) :: gf_evaluations
      real(rp) :: nodes(options%near_fermi_points), weights(options%near_fermi_points)
      real(rp) :: energy, occupation_difference
      integer :: inode, ileft, iright
      complex(rp) :: value

      if (upper <= lower) return
      call gauss_legendre(options%near_fermi_points, lower, upper, nodes, weights)
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
            error stop 'integrate_cross_subinterval: unknown finite-window coefficient'
         end select
         if (abs(occupation_difference) <= tiny(1.0_rp)) cycle
         call one_particle%get_complex(GREEN_BRANCH_KQ, ik, cmplx(energy+omega, green_eta, rp), gr_q)
         call one_particle%get_complex(GREEN_BRANCH_K, ik, cmplx(energy, -green_eta, rp), ga_shift)
         gf_evaluations = gf_evaluations + 2
         do iright = 1, size(right_ops, 3)
            do ileft = 1, size(left_ops, 3)
               value = trace_four(left_ops(:, :, ileft), gr_q, right_ops(:, :, iright), ga_shift, trace_work1, trace_work2, &
                  trace_work3)
               ra_integral(ileft, iright) = ra_integral(ileft, iright) + weights(inode)*occupation_difference*value
            end do
         end do
      end do
   end subroutine integrate_cross_subinterval

   subroutine resolve_green_energy_window(one_particle, options, omega, energy_min, energy_max, thermal_tail, omega_tail)
      class(green_function_provider), intent(in) :: one_particle
      type(green_chi0_options), intent(in) :: options
      real(rp), intent(in) :: omega(:)
      real(rp), intent(out) :: energy_min, energy_max, thermal_tail, omega_tail
      real(rp) :: effective_green_eta

      call one_particle%get_spectral_bounds(energy_min, energy_max)
      effective_green_eta = options%green_eta
      if (effective_green_eta == 0.0_rp) effective_green_eta = 0.5_rp*options%eta
      thermal_tail = 50.0_rp*max(options%electronic_temperature*6.3336814e-6_rp, tddft_occupation_kT_floor)
      omega_tail = max(0.0_rp, maxval(abs(omega)))
      if (options%energy_min < huge(1.0_rp)/2.0_rp) energy_min = options%energy_min
      if (options%energy_max > -huge(1.0_rp)/2.0_rp) energy_max = options%energy_max
      if (options%energy_min >= huge(1.0_rp)/2.0_rp) energy_min = min(energy_min, options%fermi_level) - &
         max(20.0_rp*effective_green_eta, thermal_tail, omega_tail)
      if (options%energy_max <= -huge(1.0_rp)/2.0_rp) energy_max = max(energy_max, options%fermi_level) + &
         max(20.0_rp*effective_green_eta, thermal_tail, omega_tail)
   end subroutine resolve_green_energy_window

   subroutine mixed_contour_same_interval(options, energy_min, energy_max, same_min, same_max)
      type(green_chi0_options), intent(in) :: options
      real(rp), intent(in) :: energy_min, energy_max
      real(rp), intent(out) :: same_min, same_max

      same_min = energy_min
      same_max = energy_max
      if (options%electronic_temperature <= 0.0_rp) same_max = min(energy_max, options%fermi_level)
   end subroutine mixed_contour_same_interval

   subroutine mixed_contour_height(options, energy_width, green_eta, height)
      type(green_chi0_options), intent(in) :: options
      real(rp), intent(in) :: energy_width, green_eta
      real(rp), intent(out) :: height
      real(rp) :: kT

      kT = max(options%electronic_temperature*6.3336814e-6_rp, tddft_occupation_kT_floor)
      if (options%contour_height > 0.0_rp) then
         height = options%contour_height
      else if (options%electronic_temperature > 0.0_rp) then
         height = min(0.25_rp*pi*kT, 0.25_rp*max(energy_width, green_eta))
      else
         height = 0.25_rp*max(energy_width, 4.0_rp*green_eta)
      end if
      if (height <= 0.0_rp) error stop 'mixed contour: contour height must be positive'
      if (options%electronic_temperature > 0.0_rp .and. height >= 0.95_rp*pi*kT) then
         error stop 'mixed contour: contour height crosses a Fermi-function pole; reduce contour_height'
      end if
   end subroutine mixed_contour_height

   pure complex(rp) function complex_fermi_occupation(energy, fermi_level, temperature) result(occupation)
      complex(rp), intent(in) :: energy
      real(rp), intent(in) :: fermi_level, temperature
      complex(rp) :: argument, exponential
      real(rp) :: kT

      kT = max(temperature*6.3336814e-6_rp, tddft_occupation_kT_floor)
      argument = (energy-fermi_level)/kT
      if (real(argument, rp) >= 0.0_rp) then
         exponential = exp(-argument)
         occupation = exponential/(1.0_rp+exponential)
      else
         exponential = exp(argument)
         occupation = 1.0_rp/(1.0_rp+exponential)
      end if
   end function complex_fermi_occupation

   function trace_four(a, b, c, d, work1, work2, work3) result(value)
      complex(rp), intent(in) :: a(:, :), b(:, :), c(:, :), d(:, :)
      complex(rp), intent(out) :: work1(:, :), work2(:, :), work3(:, :)
      complex(rp) :: value
      integer :: i, n

      n = size(a, 1)
      if (size(a, 2) /= n .or. any(shape(b) /= [n, n]) .or. any(shape(c) /= [n, n]) .or. any(shape(d) /= [n, n]) .or. &
          any(shape(work1) /= [n, n]) .or. any(shape(work2) /= [n, n]) .or. any(shape(work3) /= [n, n])) then
         error stop 'trace_four: incompatible matrix dimensions'
      end if
      work1 = matmul(a, b)
      work2 = matmul(work1, c)
      work3 = matmul(work2, d)
      value = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, n
         value = value + work3(i, i)
      end do
   end function trace_four

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

end module tddft_chi0_green_mod
