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
      integer, allocatable :: band_n(:), band_m(:)
      real(rp), allocatable :: occupations(:), transition_energies(:)
      complex(rp), allocatable :: left_vertices(:, :), right_vertices(:, :)
      complex(rp), allocatable :: denominators(:), weighted_left(:, :)
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
      procedure(fill_vertex_tile_interface), deferred :: fill_vertex_tile
   end type tddft_vertex_provider

   !> Site/channel provider used by ordinary chi_KS and four-component paths.
   type, extends(tddft_vertex_provider), public :: site_channel_vertex_provider
      integer, pointer :: site_orbital_counts(:) => null()
      type(response_channel), pointer :: left_channels(:) => null(), right_channels(:) => null()
      complex(rp), pointer :: eigenvectors_k(:, :, :) => null(), eigenvectors_kq(:, :, :) => null()
   contains
      procedure :: fill_vertex_tile => fill_site_channel_vertex_tile
      procedure :: restore_to_default => restore_site_channel_provider
      final :: site_channel_provider_destructor
   end type site_channel_vertex_provider

   !> Pair-potential provider. `operators_k` is a non-owning view of either a
   !> constant operator slice (extent one) or a streamed k slice.  Therefore
   !> the engine never creates a second full pair-operator tensor.
   type, extends(tddft_vertex_provider), public :: pair_operator_vertex_provider
      integer, pointer :: site_orbital_counts(:) => null()
      type(response_channel), pointer :: left_channels(:) => null()
      complex(rp), pointer :: eigenvectors_k(:, :, :) => null(), eigenvectors_kq(:, :, :) => null()
      complex(rp), pointer :: operators_k(:, :, :, :) => null()
      complex(rp), pointer :: constant_operators(:, :, :) => null()
      logical :: k_dependent = .false.
   contains
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
      procedure :: restore_to_default => transition_engine_restore_to_default
      final :: transition_engine_destructor
   end type tddft_transition_engine

   public :: make_site_channel_vertex_provider, make_pair_operator_vertex_provider
   public :: tddft_fermi_occupation, tddft_static_divided_difference

   abstract interface
      subroutine fill_vertex_tile_interface(this, ik, band_n, band_m, npairs, left_vertices, right_vertices)
         import :: tddft_vertex_provider, rp
         class(tddft_vertex_provider), intent(in) :: this
         integer, intent(in) :: ik, band_n(:), band_m(:), npairs
         complex(rp), intent(out) :: left_vertices(:, :), right_vertices(:, :)
      end subroutine fill_vertex_tile_interface
   end interface

   interface make_pair_operator_vertex_provider
      module procedure make_pair_operator_vertex_provider_k_resolved
      module procedure make_pair_operator_vertex_provider_constant
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
      provider%site_orbital_counts => site_orbital_counts; provider%left_channels => left_channels
      provider%eigenvectors_k => eigenvectors_k; provider%eigenvectors_kq => eigenvectors_kq; provider%operators_k => operators_k
      provider%k_dependent = k_dependent
   end subroutine make_pair_operator_vertex_provider_k_resolved

   subroutine make_pair_operator_vertex_provider_constant(provider, site_orbital_counts, left_channels, eigenvectors_k, &
      eigenvectors_kq, operators)
      type(pair_operator_vertex_provider), intent(out) :: provider
      integer, target, intent(in) :: site_orbital_counts(:)
      type(response_channel), target, intent(in) :: left_channels(:)
      complex(rp), target, intent(in) :: eigenvectors_k(:, :, :), eigenvectors_kq(:, :, :), operators(:, :, :)
      provider%site_orbital_counts => site_orbital_counts; provider%left_channels => left_channels
      provider%eigenvectors_k => eigenvectors_k; provider%eigenvectors_kq => eigenvectors_kq
      provider%constant_operators => operators; provider%k_dependent = .false.
   end subroutine make_pair_operator_vertex_provider_constant

   subroutine fill_site_channel_vertex_tile(this, ik, band_n, band_m, npairs, left_vertices, right_vertices)
      class(site_channel_vertex_provider), intent(in) :: this
      integer, intent(in) :: ik, band_n(:), band_m(:), npairs
      complex(rp), intent(out) :: left_vertices(:, :), right_vertices(:, :)
      complex(rp), allocatable :: bra(:, :), ket(:, :)
      integer :: ipair, nmat
      if (.not. associated(this%eigenvectors_k) .or. .not. associated(this%eigenvectors_kq)) error stop 'site vertex provider is not configured'
      nmat = size(this%eigenvectors_k,1); allocate(bra(nmat,npairs),ket(nmat,npairs))
      do ipair=1,npairs
         bra(:,ipair)=this%eigenvectors_k(:,band_n(ipair),ik); ket(:,ipair)=this%eigenvectors_kq(:,band_m(ipair),ik)
      end do
      call response_transition_vectors(this%left_channels,this%site_orbital_counts,bra,ket,left_vertices(:,1:npairs))
      do ipair=1,npairs
         bra(:,ipair)=this%eigenvectors_kq(:,band_m(ipair),ik); ket(:,ipair)=this%eigenvectors_k(:,band_n(ipair),ik)
      end do
      call response_transition_vectors(this%right_channels,this%site_orbital_counts,bra,ket,right_vertices(:,1:npairs))
   end subroutine fill_site_channel_vertex_tile

   subroutine fill_pair_operator_vertex_tile(this, ik, band_n, band_m, npairs, left_vertices, right_vertices)
      class(pair_operator_vertex_provider), intent(in) :: this
      integer, intent(in) :: ik, band_n(:), band_m(:), npairs
      complex(rp), intent(out) :: left_vertices(:, :), right_vertices(:, :)
      complex(rp), allocatable :: bra(:, :), ket(:, :)
      integer :: ipair, nmat, operator_k
      if (.not. associated(this%operators_k) .and. .not. associated(this%constant_operators)) then
         error stop 'pair vertex provider is not configured'
      end if
      operator_k=1; if (this%k_dependent) operator_k=ik
      nmat=size(this%eigenvectors_k,1); allocate(bra(nmat,npairs),ket(nmat,npairs))
      do ipair=1,npairs
         bra(:,ipair)=this%eigenvectors_k(:,band_n(ipair),ik); ket(:,ipair)=this%eigenvectors_kq(:,band_m(ipair),ik)
      end do
      call response_transition_vectors(this%left_channels,this%site_orbital_counts,bra,ket,left_vertices(:,1:npairs))
      do ipair=1,npairs
         bra(:,ipair)=this%eigenvectors_kq(:,band_m(ipair),ik); ket(:,ipair)=this%eigenvectors_k(:,band_n(ipair),ik)
      end do
      if (associated(this%constant_operators)) then
         call weighted_transition_vectors(this%constant_operators,bra,ket,right_vertices(:,1:npairs))
      else
         call weighted_transition_vectors(this%operators_k(:,:,:,operator_k),bra,ket,right_vertices(:,1:npairs))
      end if
   end subroutine fill_pair_operator_vertex_tile

   subroutine transition_engine_accumulate_dynamic(this, k_weights, eigenvalues_k, eigenvalues_kq, omega, eta, fermi_level, &
      temperature, band_first, band_last, prune_tolerance, batch_size, use_batched, provider, response, vertex_seconds, &
      preparation_seconds, denominator_seconds, accumulation_seconds)
      class(tddft_transition_engine), intent(inout) :: this
      real(rp), intent(in) :: k_weights(:), eigenvalues_k(:, :), eigenvalues_kq(:, :), omega(:), eta, fermi_level, temperature, prune_tolerance
      integer, intent(in) :: band_first, band_last, batch_size
      logical, intent(in) :: use_batched
      class(tddft_vertex_provider), intent(in) :: provider
      complex(rp), intent(inout) :: response(:, :, :)
      real(rp), intent(inout) :: vertex_seconds, preparation_seconds, denominator_seconds, accumulation_seconds
      integer :: ik,n,m,npairs,capacity,ipair,iw
      real(rp) :: weight_sum,prefactor,t0,t1
      if (eta <= 0.0_rp .or. batch_size < 1) error stop 'transition engine: invalid dynamic controls'
      weight_sum=sum(k_weights); capacity=merge(batch_size,1,use_batched)
      call this%workspace%ensure_capacity(size(response,1),size(response,2),capacity)
      do ik=1,size(k_weights)
         prefactor=k_weights(ik)/weight_sum; npairs=0
         do n=band_first,band_last
            do m=band_first,band_last
               call cpu_time(t0)
               if (prune_tolerance > 0.0_rp) then
                  if (abs(tddft_fermi_occupation(eigenvalues_k(n,ik),fermi_level,temperature)- &
                     tddft_fermi_occupation(eigenvalues_kq(m,ik),fermi_level,temperature)) <= prune_tolerance) then
                     call cpu_time(t1); preparation_seconds=preparation_seconds+t1-t0; cycle
                  end if
               end if
               npairs=npairs+1; this%workspace%band_n(npairs)=n; this%workspace%band_m(npairs)=m
               this%workspace%occupations(npairs)=tddft_fermi_occupation(eigenvalues_k(n,ik),fermi_level,temperature)- &
                  tddft_fermi_occupation(eigenvalues_kq(m,ik),fermi_level,temperature)
               this%workspace%transition_energies(npairs)=eigenvalues_k(n,ik)-eigenvalues_kq(m,ik)
               call cpu_time(t1); preparation_seconds=preparation_seconds+t1-t0
               if (npairs == capacity) then
                  call flush_dynamic(this,ik,npairs,omega,eta,prefactor,use_batched,provider,response,vertex_seconds,denominator_seconds,accumulation_seconds)
                  npairs=0
               end if
            end do
         end do
         if (npairs>0) call flush_dynamic(this,ik,npairs,omega,eta,prefactor,use_batched,provider,response,vertex_seconds,denominator_seconds,accumulation_seconds)
      end do
   end subroutine transition_engine_accumulate_dynamic

   subroutine transition_engine_accumulate_static(this, k_weights, eigenvalues, fermi_level, temperature, band_first, band_last, &
      prune_tolerance, batch_size, provider, response, vertex_seconds, preparation_seconds, accumulation_seconds)
      class(tddft_transition_engine), intent(inout) :: this
      real(rp), intent(in) :: k_weights(:), eigenvalues(:, :), fermi_level, temperature, prune_tolerance
      integer, intent(in) :: band_first, band_last, batch_size
      class(tddft_vertex_provider), intent(in) :: provider
      complex(rp), intent(inout) :: response(:, :, :)
      real(rp), intent(inout) :: vertex_seconds, preparation_seconds, accumulation_seconds
      integer :: ik,n,m,npairs; real(rp)::weight_sum,prefactor,t0,t1,factor
      weight_sum=sum(k_weights); call this%workspace%ensure_capacity(size(response,1),size(response,2),batch_size)
      do ik=1,size(k_weights)
         prefactor=k_weights(ik)/weight_sum; npairs=0
         do n=band_first,band_last; do m=band_first,band_last
            call cpu_time(t0); factor=tddft_static_divided_difference(eigenvalues(n,ik),eigenvalues(m,ik),fermi_level,temperature)
            if (prune_tolerance > 0.0_rp .and. abs(factor)<=prune_tolerance) then
               call cpu_time(t1); preparation_seconds=preparation_seconds+t1-t0; cycle
            end if
            npairs=npairs+1; this%workspace%band_n(npairs)=n; this%workspace%band_m(npairs)=m; this%workspace%occupations(npairs)=factor
            call cpu_time(t1); preparation_seconds=preparation_seconds+t1-t0
            if (npairs==batch_size) then
               call flush_static(this,ik,npairs,prefactor,provider,response,vertex_seconds,accumulation_seconds); npairs=0
            end if
         end do; end do
         if(npairs>0) call flush_static(this,ik,npairs,prefactor,provider,response,vertex_seconds,accumulation_seconds)
      end do
   end subroutine transition_engine_accumulate_static

   subroutine flush_dynamic(this,ik,npairs,omega,eta,prefactor,use_batched,provider,response,vertex_seconds,denominator_seconds,accumulation_seconds)
      class(tddft_transition_engine),intent(inout)::this
      integer,intent(in)::ik,npairs; real(rp),intent(in)::omega(:),eta,prefactor
      logical,intent(in)::use_batched
      class(tddft_vertex_provider),intent(in)::provider; complex(rp),intent(inout)::response(:,:,:)
      real(rp),intent(inout)::vertex_seconds,denominator_seconds,accumulation_seconds
      integer::iw,ipair,ileft,iright; real(rp)::t0,t1
      call cpu_time(t0); call provider%fill_vertex_tile(ik,this%workspace%band_n,this%workspace%band_m,npairs, &
         this%workspace%left_vertices,this%workspace%right_vertices); call cpu_time(t1); vertex_seconds=vertex_seconds+t1-t0
      do iw=1,size(omega)
         call cpu_time(t0); this%workspace%denominators(1:npairs)=cmplx(omega(iw)+this%workspace%transition_energies(1:npairs),eta,rp)
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
      class(tddft_vertex_provider),intent(in)::provider; complex(rp),intent(inout)::response(:,:,:)
      real(rp),intent(inout)::vertex_seconds,accumulation_seconds; integer::ipair; real(rp)::t0,t1
      call cpu_time(t0); call provider%fill_vertex_tile(ik,this%workspace%band_n,this%workspace%band_m,npairs, &
         this%workspace%left_vertices,this%workspace%right_vertices); call cpu_time(t1); vertex_seconds=vertex_seconds+t1-t0
      do ipair=1,npairs
         this%workspace%weighted_left(:,ipair)=prefactor*this%workspace%occupations(ipair)*this%workspace%left_vertices(:,ipair)
      end do
      call cpu_time(t0); call zgemm('N','T',size(response,1),size(response,2),npairs,cmplx(1.0_rp,0.0_rp,rp), &
         this%workspace%weighted_left,size(this%workspace%weighted_left,1),this%workspace%right_vertices, &
         size(this%workspace%right_vertices,1),cmplx(1.0_rp,0.0_rp,rp),response(:,:,1),size(response,1)); call cpu_time(t1)
      accumulation_seconds=accumulation_seconds+t1-t0
   end subroutine flush_static

   subroutine transition_workspace_ensure_capacity(this,nleft,nright,capacity)
      class(tddft_transition_workspace),intent(inout)::this; integer,intent(in)::nleft,nright,capacity
      if (capacity<1 .or. nleft<1 .or. nright<1) error stop 'transition workspace: invalid shape'
      ! Fortran does not guarantee short-circuit evaluation of .and.; keep
      ! every allocation inquiry separate from the corresponding size inquiry.
      ! A partially initialized workspace is never reusable.
      if (this%capacity >= capacity) then
         if (allocated(this%band_n) .and. allocated(this%band_m) .and. allocated(this%occupations) .and. &
             allocated(this%transition_energies) .and. allocated(this%left_vertices) .and. &
             allocated(this%right_vertices) .and. allocated(this%denominators) .and. allocated(this%weighted_left)) then
            if (size(this%band_n) >= capacity .and. size(this%band_m) >= capacity .and. &
                size(this%occupations) >= capacity .and. size(this%transition_energies) >= capacity .and. &
                size(this%denominators) >= capacity .and. size(this%left_vertices, 1) == nleft .and. &
                size(this%left_vertices, 2) >= capacity .and. size(this%right_vertices, 1) == nright .and. &
                size(this%right_vertices, 2) >= capacity .and. size(this%weighted_left, 1) == nleft .and. &
                size(this%weighted_left, 2) >= capacity) return
         end if
      end if
      call this%clear(); this%capacity=capacity
      allocate(this%band_n(capacity),this%band_m(capacity),this%occupations(capacity),this%transition_energies(capacity), &
         this%left_vertices(nleft,capacity),this%right_vertices(nright,capacity),this%denominators(capacity),this%weighted_left(nleft,capacity))
   end subroutine transition_workspace_ensure_capacity
   subroutine transition_workspace_clear(this)
      class(tddft_transition_workspace),intent(inout)::this
      if(allocated(this%band_n))deallocate(this%band_n);if(allocated(this%band_m))deallocate(this%band_m);if(allocated(this%occupations))deallocate(this%occupations)
      if(allocated(this%transition_energies))deallocate(this%transition_energies);if(allocated(this%left_vertices))deallocate(this%left_vertices)
      if(allocated(this%right_vertices))deallocate(this%right_vertices);if(allocated(this%denominators))deallocate(this%denominators);if(allocated(this%weighted_left))deallocate(this%weighted_left)
      this%capacity=0;this%active_pairs=0
   end subroutine transition_workspace_clear
   subroutine transition_workspace_restore_to_default(this);class(tddft_transition_workspace),intent(inout)::this;call this%clear();end subroutine
   subroutine transition_workspace_destructor(this);type(tddft_transition_workspace)::this;call this%clear();end subroutine
   subroutine transition_engine_restore_to_default(this);class(tddft_transition_engine),intent(inout)::this;call this%workspace%restore_to_default();end subroutine
   subroutine transition_engine_destructor(this);type(tddft_transition_engine)::this;call this%workspace%clear();end subroutine
   subroutine restore_site_channel_provider(this);class(site_channel_vertex_provider),intent(inout)::this;nullify(this%site_orbital_counts,this%left_channels,this%right_channels,this%eigenvectors_k,this%eigenvectors_kq);end subroutine
   subroutine site_channel_provider_destructor(this);type(site_channel_vertex_provider)::this;call this%restore_to_default();end subroutine
   subroutine restore_pair_operator_provider(this);class(pair_operator_vertex_provider),intent(inout)::this;nullify(this%site_orbital_counts,this%left_channels,this%eigenvectors_k,this%eigenvectors_kq,this%operators_k,this%constant_operators);this%k_dependent=.false.;end subroutine
   subroutine pair_operator_provider_destructor(this);type(pair_operator_vertex_provider)::this;call this%restore_to_default();end subroutine
end module tddft_transition_engine_mod
