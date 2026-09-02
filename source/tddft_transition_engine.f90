!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Tiled Lehmann transition engine shared by chi_KS and direct Xi.
!>
!> The engine owns only transition selection and response-space accumulation.
!> Vertex physics remains in a provider.  A future device implementation can
!> replace a provider with one fed by RF-05 eigenpair/operator tiles, then call
!> `accumulate_dynamic` or `accumulate_static` with the same response result
!> matrix.  In particular, providers receive a complete `(n,m)` tile; there is
!> no polymorphic call in the `(k,n,m,omega)` inner loop.
module tddft_transition_engine_mod
   use precision_mod, only: rp
   use tddft_conventions_mod, only: tddft_retarded_denominator
   use response_vertices_mod, only: response_channel, response_transition_vectors, weighted_transition_vectors
   implicit none
   private

   real(rp), parameter, public :: tddft_kB_Ry_per_K = 6.3336814e-6_rp
   real(rp), parameter, public :: tddft_occupation_kT_floor = 1.0e-10_rp

   !> Bounded, reusable host tile.  Columns share deterministic n-major,
   !> m-minor order.  The arrays are deliberately response-space arrays, never
   !> a `(k,n,m,left,right)` tensor.
   type, public :: tddft_transition_workspace
      integer :: capacity = 0
      integer :: active_pairs = 0
      integer :: coefficient_dimension = 0
      integer :: storage_allocations = 0
      integer :: capacity_reuses = 0
      integer, allocatable :: band_n(:), band_m(:)
      real(rp), allocatable :: occupations(:), transition_energies(:)
      complex(rp), allocatable :: left_vertices(:, :), right_vertices(:, :)
      complex(rp), allocatable :: denominators(:), weighted_left(:, :)
      complex(rp), allocatable :: bra(:, :), ket(:, :)
   contains
      procedure :: ensure_capacity => transition_workspace_ensure_capacity
      procedure :: clear => transition_workspace_clear
      procedure :: restore_to_default => transition_workspace_restore_to_default
      final :: transition_workspace_destructor
   end type tddft_transition_workspace

   !> Deferred tile operation.  `band_n` and `band_m` have `npairs` entries
   !> and are all at one k point.  Implementations fill both ordered factors:
   !> left=<n|A|m>, right=<m|B|n>; right is not inferred by conjugation.
   type, abstract, public :: tddft_vertex_provider
   contains
      procedure(coefficient_dimension_interface), deferred :: coefficient_dimension
      procedure(begin_accumulation_interface), deferred :: begin_accumulation
      procedure(prepare_kpoint_interface), deferred :: prepare_kpoint
      procedure(fill_vertex_tile_interface), deferred :: fill_vertex_tile
   end type tddft_vertex_provider

   !> One local-k producer for pair-potential operators.  The producer fills
   !> caller-owned `(nmat,nmat,nright)` storage and never needs a complete
   !> k-resolved tensor.  It is deliberately separate from transition tiling:
   !> one source fetch serves every transition batch at that k point.
   type, abstract, public :: pair_operator_tile_source
   contains
      procedure(pair_operator_channel_dimension_interface), deferred :: channel_dimension
      procedure(fill_pair_operator_tile_interface), deferred :: fill_operator_tile
   end type pair_operator_tile_source

   !> Thin compatibility source for the pre-GC-05 constant operator API.
   type, extends(pair_operator_tile_source), public :: constant_pair_operator_tile_source
      complex(rp), pointer :: operators(:, :, :) => null()
   contains
      procedure :: channel_dimension => constant_pair_operator_channel_dimension
      procedure :: fill_operator_tile => fill_constant_pair_operator_tile
      procedure :: restore_to_default => restore_constant_pair_operator_source
      final :: constant_pair_operator_source_destructor
   end type constant_pair_operator_tile_source

   !> Thin compatibility source for the pre-GC-05 rank-four cached API.
   type, extends(pair_operator_tile_source), public :: cached_pair_operator_tile_source
      complex(rp), pointer :: operators(:, :, :, :) => null()
   contains
      procedure :: channel_dimension => cached_pair_operator_channel_dimension
      procedure :: fill_operator_tile => fill_cached_pair_operator_tile
      procedure :: restore_to_default => restore_cached_pair_operator_source
      final :: cached_pair_operator_source_destructor
   end type cached_pair_operator_tile_source

   !> Site/channel provider used by ordinary chi_KS and four-component paths.
   type, extends(tddft_vertex_provider), public :: site_channel_vertex_provider
      integer, pointer :: site_orbital_counts(:) => null()
      type(response_channel), pointer :: left_channels(:) => null(), right_channels(:) => null()
      complex(rp), pointer :: eigenvectors_k(:, :, :) => null(), eigenvectors_kq(:, :, :) => null()
   contains
      procedure :: coefficient_dimension => site_channel_coefficient_dimension
      procedure :: begin_accumulation => begin_site_channel_accumulation
      procedure :: prepare_kpoint => prepare_site_channel_kpoint
      procedure :: fill_vertex_tile => fill_site_channel_vertex_tile
      procedure :: restore_to_default => restore_site_channel_provider
      final :: site_channel_provider_destructor
   end type site_channel_vertex_provider

   !> Pair-potential provider.  It owns one reusable operator tile and borrows
   !> a source that fills it once per k point.  No rank-four k-resolved tensor
   !> is required by this path; legacy tensor inputs are source adapters.
   type, extends(tddft_vertex_provider), public :: pair_operator_vertex_provider
      integer, pointer :: site_orbital_counts(:) => null()
      type(response_channel), pointer :: left_channels(:) => null()
      complex(rp), pointer :: eigenvectors_k(:, :, :) => null(), eigenvectors_kq(:, :, :) => null()
      class(pair_operator_tile_source), pointer :: operator_source => null()
      type(constant_pair_operator_tile_source), pointer :: owned_constant_source => null()
      type(cached_pair_operator_tile_source), pointer :: owned_cached_source => null()
      complex(rp), allocatable :: operator_tile(:, :, :)
      integer :: prepared_kpoint = 0
      integer :: operator_tile_fetches = 0
   contains
      procedure :: coefficient_dimension => pair_operator_coefficient_dimension
      procedure :: begin_accumulation => begin_pair_operator_accumulation
      procedure :: prepare_kpoint => prepare_pair_operator_kpoint
      procedure :: fill_vertex_tile => fill_pair_operator_vertex_tile
      procedure :: restore_to_default => restore_pair_operator_provider
      final :: pair_operator_provider_destructor
   end type pair_operator_vertex_provider

   !> Concrete owner of a reusable workspace and the shared selection,
   !> denominator, static divided-difference, and GEMM accumulation kernels.
   type, public :: tddft_transition_engine
      type(tddft_transition_workspace) :: workspace
   contains
      procedure :: accumulate_dynamic => transition_engine_accumulate_dynamic
      procedure :: accumulate_static => transition_engine_accumulate_static
      procedure :: accumulate_static_shifted => transition_engine_accumulate_static_shifted
      procedure :: restore_to_default => transition_engine_restore_to_default
      final :: transition_engine_destructor
   end type tddft_transition_engine

   public :: make_site_channel_vertex_provider, make_pair_operator_vertex_provider
   public :: make_constant_pair_operator_tile_source, make_cached_pair_operator_tile_source
   public :: tddft_fermi_occupation, tddft_static_divided_difference

   abstract interface
      integer function coefficient_dimension_interface(this)
         import :: tddft_vertex_provider
         class(tddft_vertex_provider), intent(in) :: this
      end function coefficient_dimension_interface

      subroutine begin_accumulation_interface(this)
         import :: tddft_vertex_provider
         class(tddft_vertex_provider), intent(inout) :: this
      end subroutine begin_accumulation_interface

      subroutine prepare_kpoint_interface(this, ik)
         import :: tddft_vertex_provider
         class(tddft_vertex_provider), intent(inout) :: this
         integer, intent(in) :: ik
      end subroutine prepare_kpoint_interface

      subroutine fill_vertex_tile_interface(this, ik, band_n, band_m, npairs, bra, ket, left_vertices, right_vertices)
         import :: tddft_vertex_provider, rp
         class(tddft_vertex_provider), intent(inout) :: this
         integer, intent(in) :: ik, band_n(:), band_m(:), npairs
         complex(rp), intent(out) :: bra(:, :), ket(:, :), left_vertices(:, :), right_vertices(:, :)
      end subroutine fill_vertex_tile_interface

      integer function pair_operator_channel_dimension_interface(this)
         import :: pair_operator_tile_source
         class(pair_operator_tile_source), intent(in) :: this
      end function pair_operator_channel_dimension_interface

      subroutine fill_pair_operator_tile_interface(this, ik, operator_tile)
         import :: pair_operator_tile_source, rp
         class(pair_operator_tile_source), intent(inout) :: this
         integer, intent(in) :: ik
         complex(rp), intent(out) :: operator_tile(:, :, :)
      end subroutine fill_pair_operator_tile_interface
   end interface

   interface make_pair_operator_vertex_provider
      module procedure make_pair_operator_vertex_provider_k_resolved
      module procedure make_pair_operator_vertex_provider_constant
      module procedure make_pair_operator_vertex_provider_streamed
   end interface make_pair_operator_vertex_provider

contains

   pure real(rp) function tddft_fermi_occupation(eigenvalue, fermi_level, temperature) result(occupation)
      real(rp), intent(in) :: eigenvalue, fermi_level, temperature
      real(rp) :: argument, kT
      kT = max(temperature*tddft_kB_Ry_per_K, tddft_occupation_kT_floor)
      argument = (eigenvalue-fermi_level)/kT
      if (argument >= 50.0_rp) then
         occupation = 0.0_rp
      else if (argument <= -50.0_rp) then
         occupation = 1.0_rp
      else
         occupation = 1.0_rp/(exp(argument)+1.0_rp)
      end if
   end function tddft_fermi_occupation

   pure real(rp) function tddft_static_divided_difference(energy_n, energy_m, fermi_level, temperature) result(value)
      real(rp), intent(in) :: energy_n, energy_m, fermi_level, temperature
      real(rp) :: delta, midpoint, scale, occupation, kT
      delta = energy_n-energy_m
      kT = max(temperature*tddft_kB_Ry_per_K, tddft_occupation_kT_floor)
      scale = max(1.0_rp, abs(energy_n), abs(energy_m), kT)
      if (abs(delta) > 32.0_rp*sqrt(epsilon(1.0_rp))*scale) then
         value = (tddft_fermi_occupation(energy_n,fermi_level,temperature)- &
                  tddft_fermi_occupation(energy_m,fermi_level,temperature))/delta
      else
         midpoint = 0.5_rp*(energy_n+energy_m); occupation = tddft_fermi_occupation(midpoint,fermi_level,temperature)
         if (abs((midpoint-fermi_level)/kT) >= 50.0_rp) then
            value = 0.0_rp
         else
            value = -occupation*(1.0_rp-occupation)/kT
         end if
      end if
   end function tddft_static_divided_difference

   subroutine make_site_channel_vertex_provider(provider, site_orbital_counts, left_channels, right_channels, &
      eigenvectors_k, eigenvectors_kq)
      type(site_channel_vertex_provider), intent(out) :: provider
      integer, target, intent(in) :: site_orbital_counts(:)
      type(response_channel), target, intent(in) :: left_channels(:), right_channels(:)
      complex(rp), target, intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :)
      provider%site_orbital_counts => site_orbital_counts; provider%left_channels => left_channels
      provider%right_channels => right_channels; provider%eigenvectors_k => eigenvectors_k; provider%eigenvectors_kq => eigenvectors_kq
   end subroutine make_site_channel_vertex_provider

   subroutine make_pair_operator_vertex_provider_k_resolved(provider, site_orbital_counts, left_channels, eigenvectors_k, &
      eigenvectors_kq, operators_k, k_dependent)
      type(pair_operator_vertex_provider), intent(out) :: provider
      integer, target, intent(in) :: site_orbital_counts(:)
      type(response_channel), target, intent(in) :: left_channels(:)
      complex(rp), target, intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :), operators_k(:, :, :, :)
      logical, intent(in) :: k_dependent

      if (.not. k_dependent .or. size(operators_k, 4) /= size(eigenvectors_k, 3)) then
         error stop 'k-resolved pair provider requires one operator tile per k point'
      end if
      allocate(provider%owned_cached_source)
      call make_cached_pair_operator_tile_source(provider%owned_cached_source, operators_k)
      provider%operator_source => provider%owned_cached_source
      call configure_pair_operator_vertex_provider(provider, site_orbital_counts, left_channels, eigenvectors_k, eigenvectors_kq)
   end subroutine make_pair_operator_vertex_provider_k_resolved

   subroutine make_pair_operator_vertex_provider_constant(provider, site_orbital_counts, left_channels, eigenvectors_k, &
      eigenvectors_kq, operators)
      type(pair_operator_vertex_provider), intent(out) :: provider
      integer, target, intent(in) :: site_orbital_counts(:)
      type(response_channel), target, intent(in) :: left_channels(:)
      complex(rp), target, intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :), operators(:, :, :)
      allocate(provider%owned_constant_source)
      call make_constant_pair_operator_tile_source(provider%owned_constant_source, operators)
      provider%operator_source => provider%owned_constant_source
      call configure_pair_operator_vertex_provider(provider, site_orbital_counts, left_channels, eigenvectors_k, eigenvectors_kq)
   end subroutine make_pair_operator_vertex_provider_constant

   !> Construct a pair provider from an arbitrary one-k source.  The caller
   !> retains the source and its producer state for the duration of the
   !> accumulation; the provider owns only the reusable active operator tile.
   subroutine make_pair_operator_vertex_provider_streamed(provider, site_orbital_counts, left_channels, eigenvectors_k, &
      eigenvectors_kq, operator_source)
      type(pair_operator_vertex_provider), intent(out) :: provider
      integer, target, intent(in) :: site_orbital_counts(:)
      type(response_channel), target, intent(in) :: left_channels(:)
      complex(rp), target, intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :)
      class(pair_operator_tile_source), target, intent(inout) :: operator_source

      provider%operator_source => operator_source
      call configure_pair_operator_vertex_provider(provider, site_orbital_counts, left_channels, eigenvectors_k, eigenvectors_kq)
   end subroutine make_pair_operator_vertex_provider_streamed

   subroutine configure_pair_operator_vertex_provider(provider, site_orbital_counts, left_channels, eigenvectors_k, eigenvectors_kq)
      type(pair_operator_vertex_provider), intent(inout) :: provider
      integer, target, intent(in) :: site_orbital_counts(:)
      type(response_channel), target, intent(in) :: left_channels(:)
      complex(rp), target, intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :)
      integer :: nmat, nright

      if (.not. associated(provider%operator_source)) error stop 'pair operator source is not configured'
      nmat = size(eigenvectors_k, 1)
      nright = provider%operator_source%channel_dimension()
      if (nmat < 1 .or. nright < 1 .or. size(eigenvectors_kq, 1) /= nmat .or. &
          size(eigenvectors_kq, 3) /= size(eigenvectors_k, 3)) then
         error stop 'pair operator provider has incompatible source or eigenpair dimensions'
      end if
      provider%site_orbital_counts => site_orbital_counts
      provider%left_channels => left_channels
      provider%eigenvectors_k => eigenvectors_k
      provider%eigenvectors_kq => eigenvectors_kq
      allocate(provider%operator_tile(nmat, nmat, nright))
      provider%prepared_kpoint = 0
      provider%operator_tile_fetches = 0
   end subroutine configure_pair_operator_vertex_provider

   subroutine make_constant_pair_operator_tile_source(source, operators)
      type(constant_pair_operator_tile_source), intent(out) :: source
      complex(rp), target, intent(in) :: operators(:, :, :)
      source%operators => operators
   end subroutine make_constant_pair_operator_tile_source

   subroutine make_cached_pair_operator_tile_source(source, operators)
      type(cached_pair_operator_tile_source), intent(out) :: source
      complex(rp), target, intent(in) :: operators(:, :, :, :)
      source%operators => operators
   end subroutine make_cached_pair_operator_tile_source

   integer function site_channel_coefficient_dimension(this) result(ncoefficient)
      class(site_channel_vertex_provider), intent(in) :: this
      if (.not. associated(this%eigenvectors_k) .or. .not. associated(this%eigenvectors_kq)) then
         error stop 'site vertex provider is not configured'
      end if
      if (size(this%eigenvectors_k, 1) /= size(this%eigenvectors_kq, 1)) then
         error stop 'site vertex provider has incompatible coefficient dimensions'
      end if
      ncoefficient = size(this%eigenvectors_k, 1)
   end function site_channel_coefficient_dimension

   subroutine begin_site_channel_accumulation(this)
      class(site_channel_vertex_provider), intent(inout) :: this
   end subroutine begin_site_channel_accumulation

   subroutine prepare_site_channel_kpoint(this, ik)
      class(site_channel_vertex_provider), intent(inout) :: this
      integer, intent(in) :: ik
   end subroutine prepare_site_channel_kpoint

   integer function pair_operator_coefficient_dimension(this) result(ncoefficient)
      class(pair_operator_vertex_provider), intent(in) :: this
      if (.not. associated(this%eigenvectors_k) .or. .not. associated(this%eigenvectors_kq)) then
         error stop 'pair vertex provider is not configured'
      end if
      if (size(this%eigenvectors_k, 1) /= size(this%eigenvectors_kq, 1)) then
         error stop 'pair vertex provider has incompatible coefficient dimensions'
      end if
      ncoefficient = size(this%eigenvectors_k, 1)
   end function pair_operator_coefficient_dimension

   integer function constant_pair_operator_channel_dimension(this) result(nright)
      class(constant_pair_operator_tile_source), intent(in) :: this
      if (.not. associated(this%operators)) error stop 'constant pair operator source is not configured'
      nright = size(this%operators, 3)
   end function constant_pair_operator_channel_dimension

   integer function cached_pair_operator_channel_dimension(this) result(nright)
      class(cached_pair_operator_tile_source), intent(in) :: this
      if (.not. associated(this%operators)) error stop 'cached pair operator source is not configured'
      nright = size(this%operators, 3)
   end function cached_pair_operator_channel_dimension

   subroutine fill_constant_pair_operator_tile(this, ik, operator_tile)
      class(constant_pair_operator_tile_source), intent(inout) :: this
      integer, intent(in) :: ik
      complex(rp), intent(out) :: operator_tile(:, :, :)
      if (.not. associated(this%operators) .or. ik < 1 .or. any(shape(operator_tile) /= shape(this%operators))) then
         error stop 'constant pair operator source tile shape is incompatible'
      end if
      operator_tile = this%operators
   end subroutine fill_constant_pair_operator_tile

   subroutine fill_cached_pair_operator_tile(this, ik, operator_tile)
      class(cached_pair_operator_tile_source), intent(inout) :: this
      integer, intent(in) :: ik
      complex(rp), intent(out) :: operator_tile(:, :, :)
      if (.not. associated(this%operators) .or. ik < 1 .or. ik > size(this%operators, 4) .or. &
          any(shape(operator_tile) /= shape(this%operators(:,:,:,ik)))) then
         error stop 'cached pair operator source tile shape is incompatible'
      end if
      operator_tile = this%operators(:,:,:,ik)
   end subroutine fill_cached_pair_operator_tile

   subroutine fill_site_channel_vertex_tile(this, ik, band_n, band_m, npairs, bra, ket, left_vertices, right_vertices)
      class(site_channel_vertex_provider), intent(inout) :: this
      integer, intent(in) :: ik, band_n(:), band_m(:), npairs
      complex(rp), intent(out) :: bra(:, :), ket(:, :), left_vertices(:, :), right_vertices(:, :)
      integer :: ipair, nmat
      if (.not. associated(this%eigenvectors_k) .or. .not. associated(this%eigenvectors_kq)) error stop 'site vertex provider is not configured'
      nmat = size(this%eigenvectors_k,1)
      if (size(bra,1) /= nmat .or. size(ket,1) /= nmat .or. size(bra,2) < npairs .or. size(ket,2) < npairs) then
         error stop 'site vertex provider scratch shape is incompatible'
      end if
      do ipair=1,npairs
         bra(:,ipair)=this%eigenvectors_k(:,band_n(ipair),ik); ket(:,ipair)=this%eigenvectors_kq(:,band_m(ipair),ik)
      end do
      call response_transition_vectors(this%left_channels,this%site_orbital_counts,bra(:,1:npairs),ket(:,1:npairs), &
         left_vertices(:,1:npairs))
      do ipair=1,npairs
         bra(:,ipair)=this%eigenvectors_kq(:,band_m(ipair),ik); ket(:,ipair)=this%eigenvectors_k(:,band_n(ipair),ik)
      end do
      call response_transition_vectors(this%right_channels,this%site_orbital_counts,bra(:,1:npairs),ket(:,1:npairs), &
         right_vertices(:,1:npairs))
   end subroutine fill_site_channel_vertex_tile

   subroutine begin_pair_operator_accumulation(this)
      class(pair_operator_vertex_provider), intent(inout) :: this
      this%prepared_kpoint = 0
   end subroutine begin_pair_operator_accumulation

   subroutine prepare_pair_operator_kpoint(this, ik)
      class(pair_operator_vertex_provider), intent(inout) :: this
      integer, intent(in) :: ik
      integer :: nmat, nright

      if (.not. associated(this%operator_source) .or. .not. allocated(this%operator_tile)) then
         error stop 'pair vertex provider is not configured'
      end if
      if (ik < 1 .or. ik > size(this%eigenvectors_k, 3)) error stop 'pair vertex provider k index is incompatible'
      if (this%prepared_kpoint == ik) return
      nmat = size(this%eigenvectors_k, 1)
      nright = this%operator_source%channel_dimension()
      if (any(shape(this%operator_tile) /= [nmat, nmat, nright])) then
         error stop 'pair vertex provider operator tile shape is incompatible'
      end if
      call this%operator_source%fill_operator_tile(ik, this%operator_tile)
      this%prepared_kpoint = ik
      this%operator_tile_fetches = this%operator_tile_fetches + 1
   end subroutine prepare_pair_operator_kpoint

   subroutine fill_pair_operator_vertex_tile(this, ik, band_n, band_m, npairs, bra, ket, left_vertices, right_vertices)
      class(pair_operator_vertex_provider), intent(inout) :: this
      integer, intent(in) :: ik, band_n(:), band_m(:), npairs
      complex(rp), intent(out) :: bra(:, :), ket(:, :), left_vertices(:, :), right_vertices(:, :)
      integer :: ipair, nmat
      if (.not. associated(this%operator_source) .or. .not. allocated(this%operator_tile)) then
         error stop 'pair vertex provider is not configured'
      end if
      if (this%prepared_kpoint /= ik) error stop 'pair vertex provider tile was not prepared for this k point'
      nmat=size(this%eigenvectors_k,1)
      if (size(bra,1) /= nmat .or. size(ket,1) /= nmat .or. size(bra,2) < npairs .or. size(ket,2) < npairs) then
         error stop 'pair vertex provider scratch shape is incompatible'
      end if
      do ipair=1,npairs
         bra(:,ipair)=this%eigenvectors_k(:,band_n(ipair),ik); ket(:,ipair)=this%eigenvectors_kq(:,band_m(ipair),ik)
      end do
      call response_transition_vectors(this%left_channels,this%site_orbital_counts,bra(:,1:npairs),ket(:,1:npairs), &
         left_vertices(:,1:npairs))
      do ipair=1,npairs
         bra(:,ipair)=this%eigenvectors_kq(:,band_m(ipair),ik); ket(:,ipair)=this%eigenvectors_k(:,band_n(ipair),ik)
      end do
      call weighted_transition_vectors(this%operator_tile,bra(:,1:npairs),ket(:,1:npairs),right_vertices(:,1:npairs))
   end subroutine fill_pair_operator_vertex_tile

   subroutine transition_engine_accumulate_dynamic(this, k_weights, eigenvalues_k, eigenvalues_kq, omega, eta, fermi_level, &
      temperature, band_first, band_last, prune_tolerance, batch_size, use_batched, provider, response, vertex_seconds, &
      preparation_seconds, denominator_seconds, accumulation_seconds)
      class(tddft_transition_engine), intent(inout) :: this
      real(rp), intent(in) :: k_weights(:), eigenvalues_k(:, :), eigenvalues_kq(:, :), omega(:), eta, fermi_level, temperature, prune_tolerance
      integer, intent(in) :: band_first, band_last, batch_size
      logical, intent(in) :: use_batched
      class(tddft_vertex_provider), intent(inout) :: provider
      complex(rp), intent(inout) :: response(:, :, :)
      real(rp), intent(inout) :: vertex_seconds, preparation_seconds, denominator_seconds, accumulation_seconds
      integer :: ik,n,m,npairs,capacity,ncoefficient
      real(rp) :: weight_sum,prefactor,t0,t1
      if (eta <= 0.0_rp .or. batch_size < 1) error stop 'transition engine: invalid dynamic controls'
      weight_sum=sum(k_weights); capacity=merge(batch_size,1,use_batched)
      ncoefficient = provider%coefficient_dimension()
      call this%workspace%ensure_capacity(size(response,1),size(response,2),ncoefficient,capacity)
      call provider%begin_accumulation()
      do ik=1,size(k_weights)
         call provider%prepare_kpoint(ik)
         prefactor=k_weights(ik)/weight_sum; npairs=0
         ! Timer calls are deliberately outside the (n,m) loop.  The old
         ! per-transition clock sampling was measurable for small orbital
         ! tiles and obscured the actual preparation/contraction crossover.
         call cpu_time(t0)
         do n=band_first,band_last
            do m=band_first,band_last
               if (prune_tolerance > 0.0_rp) then
                  if (abs(tddft_fermi_occupation(eigenvalues_k(n,ik),fermi_level,temperature)- &
                     tddft_fermi_occupation(eigenvalues_kq(m,ik),fermi_level,temperature)) <= prune_tolerance) then
                     cycle
                  end if
               end if
               npairs=npairs+1; this%workspace%band_n(npairs)=n; this%workspace%band_m(npairs)=m
               this%workspace%occupations(npairs)=tddft_fermi_occupation(eigenvalues_k(n,ik),fermi_level,temperature)- &
                  tddft_fermi_occupation(eigenvalues_kq(m,ik),fermi_level,temperature)
               this%workspace%transition_energies(npairs)=eigenvalues_k(n,ik)-eigenvalues_kq(m,ik)
               if (npairs == capacity) then
                  call cpu_time(t1); preparation_seconds=preparation_seconds+t1-t0
                  call flush_dynamic(this,ik,npairs,omega,eta,prefactor,use_batched,provider,response,vertex_seconds,denominator_seconds,accumulation_seconds)
                  npairs=0; call cpu_time(t0)
               end if
            end do
         end do
         call cpu_time(t1); preparation_seconds=preparation_seconds+t1-t0
         if (npairs>0) call flush_dynamic(this,ik,npairs,omega,eta,prefactor,use_batched,provider,response,vertex_seconds,denominator_seconds,accumulation_seconds)
      end do
   end subroutine transition_engine_accumulate_dynamic

   subroutine transition_engine_accumulate_static(this, k_weights, eigenvalues, fermi_level, temperature, band_first, band_last, &
      prune_tolerance, batch_size, provider, response, vertex_seconds, preparation_seconds, accumulation_seconds)
      class(tddft_transition_engine), intent(inout) :: this
      real(rp), intent(in) :: k_weights(:), eigenvalues(:, :), fermi_level, temperature, prune_tolerance
      integer, intent(in) :: band_first, band_last, batch_size
      class(tddft_vertex_provider), intent(inout) :: provider
      complex(rp), intent(inout) :: response(:, :, :)
      real(rp), intent(inout) :: vertex_seconds, preparation_seconds, accumulation_seconds

      ! Preserve the original q=0 API while routing it through the general
      ! endpoint-aware implementation.  This is important for existing
      ! direct engine users and makes the static endpoint relationship
      ! explicit in one place.
      call this%accumulate_static_shifted(k_weights, eigenvalues, eigenvalues, fermi_level, temperature, band_first, band_last, &
         prune_tolerance, batch_size, provider, response, vertex_seconds, preparation_seconds, accumulation_seconds)
   end subroutine transition_engine_accumulate_static

   !> Accumulate the true static divided-difference response for separate
   !> `(n,k)` and `(m,k+q)` endpoint eigenvalue arrays.  The provider already
   !> owns the corresponding eigenvectors; keeping the energy endpoints
   !> separate here prevents a q=0-only static approximation from leaking into
   !> finite-q Ward and convergence tests.
   subroutine transition_engine_accumulate_static_shifted(this, k_weights, eigenvalues_k, eigenvalues_kq, fermi_level, temperature, &
      band_first, band_last, prune_tolerance, batch_size, provider, response, vertex_seconds, preparation_seconds, accumulation_seconds)
      class(tddft_transition_engine), intent(inout) :: this
      real(rp), intent(in) :: k_weights(:), eigenvalues_k(:, :), eigenvalues_kq(:, :), fermi_level, temperature, prune_tolerance
      integer, intent(in) :: band_first, band_last, batch_size
      class(tddft_vertex_provider), intent(inout) :: provider
      complex(rp), intent(inout) :: response(:, :, :)
      real(rp), intent(inout) :: vertex_seconds, preparation_seconds, accumulation_seconds
      integer :: ik,n,m,npairs,ncoefficient; real(rp)::weight_sum,prefactor,t0,t1,factor
      if (size(eigenvalues_k, 2) /= size(k_weights) .or. any(shape(eigenvalues_kq) /= shape(eigenvalues_k))) then
         error stop 'transition engine: incompatible static endpoint eigenvalues'
      end if
      if (band_first < 1 .or. band_last < band_first .or. band_last > size(eigenvalues_k, 1)) then
         error stop 'transition engine: invalid static band window'
      end if
      if (batch_size < 1) error stop 'transition engine: invalid static batch size'
      weight_sum=sum(k_weights); ncoefficient=provider%coefficient_dimension()
      call this%workspace%ensure_capacity(size(response,1),size(response,2),ncoefficient,batch_size)
      call provider%begin_accumulation()
      do ik=1,size(k_weights)
         call provider%prepare_kpoint(ik)
         prefactor=k_weights(ik)/weight_sum; npairs=0
         call cpu_time(t0)
         do n=band_first,band_last; do m=band_first,band_last
            factor=tddft_static_divided_difference(eigenvalues_k(n,ik),eigenvalues_kq(m,ik),fermi_level,temperature)
            if (prune_tolerance > 0.0_rp .and. abs(factor)<=prune_tolerance) then
               cycle
            end if
            npairs=npairs+1; this%workspace%band_n(npairs)=n; this%workspace%band_m(npairs)=m; this%workspace%occupations(npairs)=factor
            if (npairs==batch_size) then
               call cpu_time(t1); preparation_seconds=preparation_seconds+t1-t0
               call flush_static(this,ik,npairs,prefactor,provider,response,vertex_seconds,accumulation_seconds); npairs=0
               call cpu_time(t0)
            end if
         end do; end do
         call cpu_time(t1); preparation_seconds=preparation_seconds+t1-t0
         if(npairs>0) call flush_static(this,ik,npairs,prefactor,provider,response,vertex_seconds,accumulation_seconds)
      end do
   end subroutine transition_engine_accumulate_static_shifted

   subroutine flush_dynamic(this,ik,npairs,omega,eta,prefactor,use_batched,provider,response,vertex_seconds,denominator_seconds,accumulation_seconds)
      class(tddft_transition_engine),intent(inout)::this
      integer,intent(in)::ik,npairs; real(rp),intent(in)::omega(:),eta,prefactor
      logical,intent(in)::use_batched
      class(tddft_vertex_provider),intent(inout)::provider; complex(rp),intent(inout)::response(:,:,:)
      real(rp),intent(inout)::vertex_seconds,denominator_seconds,accumulation_seconds
      integer::iw,ipair,ileft,iright; real(rp)::t0,t1
      call cpu_time(t0); call provider%fill_vertex_tile(ik,this%workspace%band_n,this%workspace%band_m,npairs, &
         this%workspace%bra,this%workspace%ket,this%workspace%left_vertices,this%workspace%right_vertices)
      call cpu_time(t1); vertex_seconds=vertex_seconds+t1-t0
      do iw=1,size(omega)
         call cpu_time(t0)
         this%workspace%denominators(1:npairs)=tddft_retarded_denominator(omega(iw), &
            this%workspace%transition_energies(1:npairs), eta)
         do ipair=1,npairs
            this%workspace%weighted_left(:,ipair)=prefactor*this%workspace%occupations(ipair)*this%workspace%left_vertices(:,ipair)/this%workspace%denominators(ipair)
         end do
         call cpu_time(t1); denominator_seconds=denominator_seconds+t1-t0
         call cpu_time(t0)
         if (use_batched) then
            call zgemm('N','T',size(response,1),size(response,2),npairs,cmplx(1.0_rp,0.0_rp,rp), &
               this%workspace%weighted_left,size(this%workspace%weighted_left,1),this%workspace%right_vertices, &
               size(this%workspace%right_vertices,1),cmplx(1.0_rp,0.0_rp,rp),response(:,:,iw),size(response,1))
         else
            ! Independent scalar reduction oracle: deliberately no GEMM.
            do iright=1,size(response,2); do ileft=1,size(response,1)
               response(ileft,iright,iw)=response(ileft,iright,iw)+this%workspace%weighted_left(ileft,1)*this%workspace%right_vertices(iright,1)
            end do; end do
         end if
         call cpu_time(t1)
         accumulation_seconds=accumulation_seconds+t1-t0
      end do
   end subroutine flush_dynamic

   subroutine flush_static(this,ik,npairs,prefactor,provider,response,vertex_seconds,accumulation_seconds)
      class(tddft_transition_engine),intent(inout)::this; integer,intent(in)::ik,npairs; real(rp),intent(in)::prefactor
      class(tddft_vertex_provider),intent(inout)::provider; complex(rp),intent(inout)::response(:,:,:)
      real(rp),intent(inout)::vertex_seconds,accumulation_seconds; integer::ipair; real(rp)::t0,t1
      call cpu_time(t0); call provider%fill_vertex_tile(ik,this%workspace%band_n,this%workspace%band_m,npairs, &
         this%workspace%bra,this%workspace%ket,this%workspace%left_vertices,this%workspace%right_vertices)
      call cpu_time(t1); vertex_seconds=vertex_seconds+t1-t0
      do ipair=1,npairs
         this%workspace%weighted_left(:,ipair)=prefactor*this%workspace%occupations(ipair)*this%workspace%left_vertices(:,ipair)
      end do
      call cpu_time(t0); call zgemm('N','T',size(response,1),size(response,2),npairs,cmplx(1.0_rp,0.0_rp,rp), &
         this%workspace%weighted_left,size(this%workspace%weighted_left,1),this%workspace%right_vertices, &
         size(this%workspace%right_vertices,1),cmplx(1.0_rp,0.0_rp,rp),response(:,:,1),size(response,1)); call cpu_time(t1)
      accumulation_seconds=accumulation_seconds+t1-t0
   end subroutine flush_static

   subroutine transition_workspace_ensure_capacity(this,nleft,nright,ncoefficient,capacity)
      class(tddft_transition_workspace),intent(inout)::this; integer,intent(in)::nleft,nright,ncoefficient,capacity
      if (capacity<1 .or. nleft<1 .or. nright<1 .or. ncoefficient<1) error stop 'transition workspace: invalid shape'
      ! Fortran does not guarantee short-circuit evaluation of .and.; keep
      ! every allocation inquiry separate from the corresponding size inquiry.
      ! A partially initialized workspace is never reusable.
      if (this%capacity >= capacity) then
         if (allocated(this%band_n) .and. allocated(this%band_m) .and. allocated(this%occupations) .and. &
             allocated(this%transition_energies) .and. allocated(this%left_vertices) .and. &
             allocated(this%right_vertices) .and. allocated(this%denominators) .and. allocated(this%weighted_left) .and. &
             allocated(this%bra) .and. allocated(this%ket)) then
            if (size(this%band_n) >= capacity .and. size(this%band_m) >= capacity .and. &
                size(this%occupations) >= capacity .and. size(this%transition_energies) >= capacity .and. &
                size(this%denominators) >= capacity .and. size(this%left_vertices, 1) == nleft .and. &
                size(this%left_vertices, 2) >= capacity .and. size(this%right_vertices, 1) == nright .and. &
                size(this%right_vertices, 2) >= capacity .and. size(this%weighted_left, 1) == nleft .and. &
                size(this%weighted_left, 2) >= capacity .and. size(this%bra,1) == ncoefficient .and. &
                size(this%bra,2) >= capacity .and. size(this%ket,1) == ncoefficient .and. size(this%ket,2) >= capacity) then
               this%capacity_reuses=this%capacity_reuses+1
               return
            end if
         end if
      end if
      call this%clear(); this%capacity=capacity; this%coefficient_dimension=ncoefficient
      this%storage_allocations=this%storage_allocations+1
      allocate(this%band_n(capacity),this%band_m(capacity),this%occupations(capacity),this%transition_energies(capacity), &
         this%left_vertices(nleft,capacity),this%right_vertices(nright,capacity),this%denominators(capacity),this%weighted_left(nleft,capacity), &
         this%bra(ncoefficient,capacity),this%ket(ncoefficient,capacity))
   end subroutine transition_workspace_ensure_capacity
   subroutine transition_workspace_clear(this)
      class(tddft_transition_workspace),intent(inout)::this
      if(allocated(this%band_n))deallocate(this%band_n);if(allocated(this%band_m))deallocate(this%band_m);if(allocated(this%occupations))deallocate(this%occupations)
      if(allocated(this%transition_energies))deallocate(this%transition_energies);if(allocated(this%left_vertices))deallocate(this%left_vertices)
      if(allocated(this%right_vertices))deallocate(this%right_vertices);if(allocated(this%denominators))deallocate(this%denominators);if(allocated(this%weighted_left))deallocate(this%weighted_left)
      if(allocated(this%bra))deallocate(this%bra);if(allocated(this%ket))deallocate(this%ket)
      this%capacity=0;this%active_pairs=0;this%coefficient_dimension=0
   end subroutine transition_workspace_clear
   subroutine transition_workspace_restore_to_default(this);class(tddft_transition_workspace),intent(inout)::this;call this%clear();this%storage_allocations=0;this%capacity_reuses=0;end subroutine
   subroutine transition_workspace_destructor(this);type(tddft_transition_workspace)::this;call this%clear();end subroutine
   subroutine transition_engine_restore_to_default(this);class(tddft_transition_engine),intent(inout)::this;call this%workspace%restore_to_default();end subroutine
   subroutine transition_engine_destructor(this);type(tddft_transition_engine)::this;call this%workspace%clear();end subroutine
   subroutine restore_site_channel_provider(this);class(site_channel_vertex_provider),intent(inout)::this;nullify(this%site_orbital_counts,this%left_channels,this%right_channels,this%eigenvectors_k,this%eigenvectors_kq);end subroutine
   subroutine site_channel_provider_destructor(this);type(site_channel_vertex_provider)::this;call this%restore_to_default();end subroutine
   subroutine restore_pair_operator_provider(this)
      class(pair_operator_vertex_provider), intent(inout) :: this
      if (associated(this%owned_constant_source)) deallocate(this%owned_constant_source)
      if (associated(this%owned_cached_source)) deallocate(this%owned_cached_source)
      if (allocated(this%operator_tile)) deallocate(this%operator_tile)
      nullify(this%site_orbital_counts, this%left_channels, this%eigenvectors_k, this%eigenvectors_kq, this%operator_source)
      this%prepared_kpoint = 0
      this%operator_tile_fetches = 0
   end subroutine restore_pair_operator_provider
   subroutine pair_operator_provider_destructor(this);type(pair_operator_vertex_provider)::this;call this%restore_to_default();end subroutine
   subroutine restore_constant_pair_operator_source(this)
      class(constant_pair_operator_tile_source), intent(inout) :: this
      nullify(this%operators)
   end subroutine restore_constant_pair_operator_source
   subroutine constant_pair_operator_source_destructor(this)
      type(constant_pair_operator_tile_source) :: this
      call this%restore_to_default()
   end subroutine constant_pair_operator_source_destructor
   subroutine restore_cached_pair_operator_source(this)
      class(cached_pair_operator_tile_source), intent(inout) :: this
      nullify(this%operators)
   end subroutine restore_cached_pair_operator_source
   subroutine cached_pair_operator_source_destructor(this)
      type(cached_pair_operator_tile_source) :: this
      call this%restore_to_default()
   end subroutine cached_pair_operator_source_destructor
end module tddft_transition_engine_mod
