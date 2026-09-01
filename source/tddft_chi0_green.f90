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
   use math_mod, only: pi, i_unit
   use tddft_conventions_mod, only: tddft_retarded_green_denominator
   use response_vertices_mod, only: response_channel, site_projected_operator
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, tddft_fermi_occupation, &
      tddft_occupation_kT_floor
   implicit none

   private
   integer, parameter, public :: GREEN_BRANCH_K = 1
   integer, parameter, public :: GREEN_BRANCH_KQ = 2

   !> Minimum information a response bubble needs from a one-electron route.
   type, abstract, public :: green_function_provider
   contains
      procedure(green_retarded_matrix), deferred :: get_retarded
      procedure(green_spectral_bounds), deferred :: get_spectral_bounds
   end type green_function_provider

   abstract interface
      subroutine green_retarded_matrix(this, branch, ik, energy, eta, green_matrix)
         import :: green_function_provider, rp
         class(green_function_provider), intent(in) :: this
         integer, intent(in) :: branch, ik
         real(rp), intent(in) :: energy, eta
         complex(rp), intent(out) :: green_matrix(:, :)
      end subroutine green_retarded_matrix

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
      integer :: k_mesh_shape(3) = 0
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
      procedure :: get_spectral_bounds => eigenpair_spectral_bounds
   end type eigenpair_green_function_provider

   public :: build_chi_ks_from_green_functions
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
      integer :: n, nmat, i, j

      if (eta <= 0.0_rp) error stop 'eigenpair_retarded_green: eta must be positive'
      nmat = size(this%eigenvalues_k, 1)
      if (ik < 1 .or. ik > size(this%eigenvalues_k, 2) .or. any(shape(green_matrix) /= [nmat, nmat])) then
         error stop 'eigenpair_retarded_green: invalid Green-function request'
      end if
      green_matrix = cmplx(0.0_rp, 0.0_rp, rp)
      select case (branch)
      case (GREEN_BRANCH_K)
         do n = 1, nmat
            do j = 1, nmat
               do i = 1, nmat
                  green_matrix(i, j) = green_matrix(i, j) + this%eigenvectors_k(i, n, ik)* &
                     conjg(this%eigenvectors_k(j, n, ik))/tddft_retarded_green_denominator(energy, &
                        this%eigenvalues_k(n, ik), eta)
               end do
            end do
         end do
      case (GREEN_BRANCH_KQ)
         do n = 1, nmat
            do j = 1, nmat
               do i = 1, nmat
                  green_matrix(i, j) = green_matrix(i, j) + this%eigenvectors_kq(i, n, ik)* &
                     conjg(this%eigenvectors_kq(j, n, ik))/tddft_retarded_green_denominator(energy, &
                        this%eigenvalues_kq(n, ik), eta)
               end do
            end do
         end do
      case default
         error stop 'eigenpair_retarded_green: unknown Green-function branch'
      end select
   end subroutine eigenpair_retarded_green

   subroutine eigenpair_spectral_bounds(this, energy_min, energy_max)
      class(eigenpair_green_function_provider), intent(in) :: this
      real(rp), intent(out) :: energy_min, energy_max

      if (.not. allocated(this%eigenvalues_k)) error stop 'eigenpair_spectral_bounds: provider is not initialized'
      energy_min = min(minval(this%eigenvalues_k), minval(this%eigenvalues_kq))
      energy_max = max(maxval(this%eigenvalues_k), maxval(this%eigenvalues_kq))
   end subroutine eigenpair_spectral_bounds

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

   subroutine build_green_chi0(this, k_weights, site_orbital_counts, left_channels, right_channels, omega, result)
      class(green_chi0_provider), intent(in) :: this
      real(rp), intent(in) :: k_weights(:), omega(:)
      integer, intent(in) :: site_orbital_counts(:)
      type(response_channel), intent(in) :: left_channels(:), right_channels(:)
      type(tddft_chi0_result), intent(out) :: result
      complex(rp), allocatable :: left_ops(:, :, :), right_ops(:, :, :)
      complex(rp), allocatable :: gr_k(:, :), ga_k(:, :), gr_q(:, :), ga_q(:, :), gr_shift(:, :), ga_shift(:, :)
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
      call this%one_particle%get_spectral_bounds(energy_min, energy_max)
      thermal_tail = 50.0_rp*max(this%options%electronic_temperature*6.3336814e-6_rp, tddft_occupation_kT_floor)
      omega_tail = max(0.0_rp, maxval(abs(omega)))
      if (this%options%energy_min < huge(1.0_rp)/2.0_rp) energy_min = this%options%energy_min
      if (this%options%energy_max > -huge(1.0_rp)/2.0_rp) energy_max = this%options%energy_max
      if (this%options%energy_min >= huge(1.0_rp)/2.0_rp) energy_min = min(energy_min, this%options%fermi_level) - &
         max(20.0_rp*green_eta, thermal_tail, omega_tail)
      if (this%options%energy_max <= -huge(1.0_rp)/2.0_rp) energy_max = max(energy_max, this%options%fermi_level) + &
         max(20.0_rp*green_eta, thermal_tail, omega_tail)
      if (energy_max <= energy_min) error stop 'green_chi0_provider%build: invalid integration window'
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
         gr_shift(nmat, nmat), ga_shift(nmat, nmat))

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
                     result%chi(ileft, iright, iw) = result%chi(ileft, iright, iw) - prefactor*quadrature_weight* &
                        tddft_fermi_occupation(energy, this%options%fermi_level, this%options%electronic_temperature)* &
                        (trace_four(left_ops(:, :, ileft), gr_shift, right_ops(:, :, iright), gr_k-ga_k) + &
                         trace_four(left_ops(:, :, ileft), gr_q-ga_q, right_ops(:, :, iright), ga_shift))/(2.0_rp*pi*i_unit)
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
      result%metadata%convergence_status = 'not assessed by backend'
      result%metadata%converged = .false.
      result%metadata%q_batch_size = 1
      result%metadata%omega_batch_size = nw
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

   function trace_four(a, b, c, d) result(value)
      complex(rp), intent(in) :: a(:, :), b(:, :), c(:, :), d(:, :)
      complex(rp) :: value
      complex(rp), allocatable :: work1(:, :), work2(:, :), work3(:, :)
      integer :: i, n

      n = size(a, 1)
      if (size(a, 2) /= n .or. any(shape(b) /= [n, n]) .or. any(shape(c) /= [n, n]) .or. any(shape(d) /= [n, n])) then
         error stop 'trace_four: incompatible matrix dimensions'
      end if
      allocate(work1(n, n), work2(n, n), work3(n, n))
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
