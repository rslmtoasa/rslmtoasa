!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Native coefficient-GF to Pauli/no-SOC endpoint augmentation.
!>
!> This module is deliberately separate from green_mod.  `green%auxiliary_gij`
!> changes the LMTO Green-function representation by endpoint screening factors;
!> it does not construct the physical radial map
!>
!>     A = Phi + Phidot h^gamma.
!>
!> The adapter below implements that map for the certified orthogonal,
!> second-order/HOH, collinear, no-SOC, no-additive-operator sp/spd slice.
!------------------------------------------------------------------------------
module lr_gf_endpoint_augmentation_mod

   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use response_angular_basis_mod, only: response_harmonic
   implicit none
   private

   character(len=*), parameter, public :: lr_gf_endpoint_representation = &
      'Pauli/no-SOC spatial GF, factored LMTO augmentation'

   !> Capability tuple accepted by this endpoint adapter.
   type, public :: lr_gf_endpoint_capabilities
      character(len=32) :: reciprocal_mode = 'ham_only'
      character(len=16) :: hamiltonian_order = 'second'
      logical :: orthogonal = .true.
      logical :: collinear = .true.
      logical :: has_soc = .false.
      logical :: has_extra_operator = .false.
      ! These two explicit names make the fail-closed boundary readable to
      ! callers that distinguish Hubbard from other additive operators.
      logical :: has_hubbard = .false.
      logical :: has_additive_operator = .false.
   contains
      procedure :: supported => lr_gf_endpoint_capabilities_supported
      procedure :: require_supported => lr_gf_endpoint_require_supported
   end type lr_gf_endpoint_capabilities

   !> Callback for a production effective-Hamiltonian action.
   !>
   !> `seed` is a localized source-site block.  The implementation applies the
   !> same H_eff action used by the native solver and returns the destination
   !> site block.  For destination=source the caller subtracts E_nu_work to
   !> obtain h^gamma_aa; for an offsite block the returned H block is already
   !> h^gamma_ab.
   type, abstract, public :: lr_gf_effective_hamiltonian_provider
   contains
      procedure(lr_gf_apply_effective_action), deferred :: apply
   end type lr_gf_effective_hamiltonian_provider

   abstract interface
      subroutine lr_gf_apply_effective_action(this, source_site, destination_site, seed, destination_block)
         import :: lr_gf_effective_hamiltonian_provider, rp
         class(lr_gf_effective_hamiltonian_provider), intent(inout) :: this
         integer, intent(in) :: source_site, destination_site
         complex(rp), intent(in) :: seed(:, :)
         complex(rp), intent(out) :: destination_block(:, :)
      end subroutine lr_gf_apply_effective_action
   end interface

   !> Simple dense provider used by small finite fixtures and by callers that
   !> already hold the completed H_eff matrix.  Production callers can provide
   !> a wrapper around their native two-sweep/action routine instead.
   type, extends(lr_gf_effective_hamiltonian_provider), public :: lr_gf_dense_hamiltonian_provider
      complex(rp), allocatable :: h_eff(:, :)
      integer :: block_size = 0
   contains
      procedure :: initialize => lr_gf_dense_provider_initialize
      procedure :: apply => lr_gf_dense_provider_apply
   end type lr_gf_dense_hamiltonian_provider

   !> Factored endpoint-augmented Green-function block.
   !>
   !> The four coefficient-space members are the explicit radial branches in
   !>   Phi G Phi^dagger,
   !>   Phidot (hG) Phi^dagger,
   !>   Phi (Gh) Phidot^dagger,
   !>   Phidot (hGh) Phidot^dagger.
   !>
   !> `phi_left/right` and `phidot_left/right` contain large-component radial
   !> numerators repeated in orbital/spin ordering.  Angular harmonics are not
   !> duplicated in storage; `evaluate_point` inserts the certified
   !> response_angular_basis convention and divides U_l by r at the endpoint.
   type, public :: lr_gf_augmented_block
      integer :: left_site = 0
      integer :: right_site = 0
      integer :: norb = 0
      integer :: nspin = 0
      complex(rp) :: z = (0.0_rp, 0.0_rp)
      complex(rp), allocatable :: coefficient_gf(:, :)
      complex(rp), allocatable :: hgamma_ab(:, :)
      complex(rp), allocatable :: h_g(:, :)
      complex(rp), allocatable :: g_h(:, :)
      complex(rp), allocatable :: h_g_h(:, :)
      ! Shape: (radial point, orbital, spin), large component only.
      real(rp), allocatable :: phi_left(:, :, :), phidot_left(:, :, :)
      real(rp), allocatable :: phi_right(:, :, :), phidot_right(:, :, :)
      real(rp), allocatable :: radius_left(:), radius_right(:)
   contains
      procedure :: evaluate_point => lr_gf_augmented_block_evaluate_point
      procedure :: branch_at_point => lr_gf_augmented_block_branch_at_point
   end type lr_gf_augmented_block

   public :: augment_lr_gf_endpoint
   public :: augment_lr_gf_endpoint_pair

contains

   pure logical function lr_gf_endpoint_capabilities_supported(this) result(ok)
      class(lr_gf_endpoint_capabilities), intent(in) :: this

      ok = trim(this%reciprocal_mode) == 'ham_only' .and. trim(this%hamiltonian_order) == 'second' .and. &
           this%orthogonal .and. this%collinear .and. (.not. this%has_soc) .and. &
           (.not. this%has_extra_operator) .and. (.not. this%has_hubbard) .and. &
           (.not. this%has_additive_operator)
   end function lr_gf_endpoint_capabilities_supported

   subroutine lr_gf_endpoint_require_supported(this, caller)
      class(lr_gf_endpoint_capabilities), intent(in) :: this
      character(len=*), intent(in), optional :: caller
      character(len=96) :: message

      if (.not. this%supported()) then
         if (present(caller)) then
            message = trim(caller)//': native-GF endpoint augmentation capability is unsupported'
         else
            message = 'native-GF endpoint augmentation capability is unsupported'
         end if
         error stop trim(message)
      end if
   end subroutine lr_gf_endpoint_require_supported

   subroutine lr_gf_dense_provider_initialize(this, h_eff, block_size)
      class(lr_gf_dense_hamiltonian_provider), intent(out) :: this
      complex(rp), intent(in) :: h_eff(:, :)
      integer, intent(in) :: block_size

      if (block_size < 1 .or. size(h_eff, 1) /= size(h_eff, 2) .or. &
          mod(size(h_eff, 1), block_size) /= 0) then
         error stop 'lr_gf_dense_provider_initialize: inconsistent matrix/block dimensions'
      end if
      allocate(this%h_eff(size(h_eff, 1), size(h_eff, 2)))
      this%h_eff = h_eff
      this%block_size = block_size
   end subroutine lr_gf_dense_provider_initialize

   subroutine lr_gf_dense_provider_apply(this, source_site, destination_site, seed, destination_block)
      class(lr_gf_dense_hamiltonian_provider), intent(inout) :: this
      integer, intent(in) :: source_site, destination_site
      complex(rp), intent(in) :: seed(:, :)
      complex(rp), intent(out) :: destination_block(:, :)
      integer :: nsite, source_first, destination_first

      if (.not. allocated(this%h_eff) .or. this%block_size < 1) then
         error stop 'lr_gf_dense_provider_apply: provider is not initialized'
      end if
      nsite = size(this%h_eff, 1)/this%block_size
      if (source_site < 1 .or. source_site > nsite .or. destination_site < 1 .or. &
          destination_site > nsite .or. size(seed, 1) /= this%block_size .or. &
          size(destination_block, 1) /= this%block_size .or. size(destination_block, 2) /= size(seed, 2)) then
         error stop 'lr_gf_dense_provider_apply: invalid site or seed shape'
      end if

      source_first = (source_site - 1)*this%block_size + 1
      destination_first = (destination_site - 1)*this%block_size + 1
      destination_block = matmul(this%h_eff(destination_first:destination_first + this%block_size - 1, &
                                             source_first:source_first + this%block_size - 1), seed)
   end subroutine lr_gf_dense_provider_apply

   !> Promote one native coefficient-space G_ab(z) block to a factored
   !> Pauli/no-SOC endpoint representation.
   !>
   !> hgamma_block is optional only because a provider may obtain the block
   !> from the production H_eff action.  Exactly one of hgamma_block and
   !> hamiltonian_provider must be supplied.  `reverse_coefficient_gf` and
   !> `reverse_hgamma_block` are accepted for pair callers that need G_ba; the
   !> forward result itself uses only the explicitly supplied forward blocks.
   subroutine augment_lr_gf_endpoint(radial_left, radial_right, left_site, right_site, z, coefficient_gf, &
                                     augmented, capabilities, hgamma_block, hamiltonian_provider, &
                                     reverse_coefficient_gf, reverse_hgamma_block)
      type(lmto_radial_basis), intent(in) :: radial_left, radial_right
      integer, intent(in) :: left_site, right_site
      complex(rp), intent(in) :: z
      complex(rp), intent(in) :: coefficient_gf(:, :)
      type(lr_gf_augmented_block), intent(out) :: augmented
      type(lr_gf_endpoint_capabilities), intent(in), optional :: capabilities
      complex(rp), intent(in), optional :: hgamma_block(:, :)
      class(lr_gf_effective_hamiltonian_provider), intent(inout), optional :: hamiltonian_provider
      complex(rp), intent(in), optional :: reverse_coefficient_gf(:, :)
      complex(rp), intent(in), optional :: reverse_hgamma_block(:, :)

      type(lr_gf_endpoint_capabilities) :: caps
      complex(rp), allocatable :: hgamma(:, :)
      complex(rp), allocatable :: d_left(:, :), d_right(:, :), identity(:, :)
      complex(rp), allocatable :: seed(:, :), effective_block(:, :)
      integer :: nlocal, nspin, norb
      logical :: have_hgamma

      caps = lr_gf_endpoint_capabilities()
      if (present(capabilities)) caps = capabilities
      call caps%require_supported('augment_lr_gf_endpoint')

      call validate_radial_pair(radial_left, radial_right, coefficient_gf)
      if (left_site < 1 .or. right_site < 1) then
         error stop 'augment_lr_gf_endpoint: endpoint site indices must be positive'
      end if
      norb = (radial_left%lmax + 1)**2
      nspin = 2
      nlocal = nspin*norb
      have_hgamma = present(hgamma_block) .or. present(hamiltonian_provider)
      if (present(hgamma_block) .and. present(hamiltonian_provider)) then
         error stop 'augment_lr_gf_endpoint: supply hgamma_block or provider, not both'
      end if
      if (.not. have_hgamma) then
         error stop 'augment_lr_gf_endpoint: hgamma_ab requires an explicit block or effective-action provider'
      end if
      if (present(reverse_coefficient_gf)) then
         if (size(reverse_coefficient_gf, 1) /= nlocal .or. size(reverse_coefficient_gf, 2) /= nlocal) then
            error stop 'augment_lr_gf_endpoint: reverse GF block shape mismatch'
         end if
      end if
      if (present(reverse_hgamma_block)) then
         if (size(reverse_hgamma_block, 1) /= nlocal .or. size(reverse_hgamma_block, 2) /= nlocal) then
            error stop 'augment_lr_gf_endpoint: reverse hgamma block shape mismatch'
         end if
      end if

      allocate(hgamma(nlocal, nlocal), d_left(nlocal, nlocal), d_right(nlocal, nlocal), &
               identity(nlocal, nlocal), seed(nlocal, nlocal), effective_block(nlocal, nlocal))
      if (present(hgamma_block)) then
         if (any(shape(hgamma_block) /= [nlocal, nlocal])) then
            error stop 'augment_lr_gf_endpoint: hgamma block shape mismatch'
         end if
         hgamma = hgamma_block
      else
         identity = cmplx(0.0_rp, 0.0_rp, rp)
         do nlocal = 1, size(identity, 1)
            identity(nlocal, nlocal) = cmplx(1.0_rp, 0.0_rp, rp)
         end do
         seed = identity
         call hamiltonian_provider%apply(right_site, left_site, seed, effective_block)
         hgamma = effective_block
         if (left_site == right_site) hgamma = hgamma - diagonal_enu(radial_left)
      end if

      d_left = diagonal_difference(radial_left, z)
      d_right = diagonal_difference(radial_right, z)
      identity = cmplx(0.0_rp, 0.0_rp, rp)
      do nlocal = 1, size(identity, 1)
         identity(nlocal, nlocal) = cmplx(1.0_rp, 0.0_rp, rp)
      end do

      augmented%left_site = left_site
      augmented%right_site = right_site
      augmented%norb = norb
      augmented%nspin = nspin
      augmented%z = z
      allocate(augmented%coefficient_gf(nlocal, nlocal), augmented%hgamma_ab(nlocal, nlocal), &
               augmented%h_g(nlocal, nlocal), augmented%g_h(nlocal, nlocal), augmented%h_g_h(nlocal, nlocal))
      augmented%coefficient_gf = coefficient_gf
      augmented%hgamma_ab = hgamma
      ! Keep the contact terms visible.  They are not optional onsite
      ! corrections and they vanish algebraically only for offsite blocks.
      augmented%h_g = matmul(d_left, coefficient_gf)
      augmented%g_h = matmul(coefficient_gf, d_right)
      if (left_site == right_site) then
         augmented%h_g = augmented%h_g - identity
         augmented%g_h = augmented%g_h - identity
      end if
      augmented%h_g_h = matmul(d_left, matmul(coefficient_gf, d_right))
      if (left_site == right_site) augmented%h_g_h = augmented%h_g_h - d_left
      augmented%h_g_h = augmented%h_g_h - hgamma

      call copy_radial_maps(radial_left, radial_right, augmented)
      deallocate(hgamma, d_left, d_right, identity, seed, effective_block)
   end subroutine augment_lr_gf_endpoint

   !> Convenience pair entry point.  It makes the reverse G_ba/hgamma_ba
   !> requirement explicit while retaining the same four-branch representation
   !> for each direction.
   subroutine augment_lr_gf_endpoint_pair(radial_left, radial_right, left_site, right_site, z, coefficient_gf, &
                                          reverse_coefficient_gf, forward, reverse, capabilities, hgamma_block, &
                                          reverse_hgamma_block, hamiltonian_provider)
      type(lmto_radial_basis), intent(in) :: radial_left, radial_right
      integer, intent(in) :: left_site, right_site
      complex(rp), intent(in) :: z
      complex(rp), intent(in) :: coefficient_gf(:, :), reverse_coefficient_gf(:, :)
      type(lr_gf_augmented_block), intent(out) :: forward, reverse
      type(lr_gf_endpoint_capabilities), intent(in), optional :: capabilities
      complex(rp), intent(in), optional :: hgamma_block(:, :), reverse_hgamma_block(:, :)
      class(lr_gf_effective_hamiltonian_provider), intent(inout), optional :: hamiltonian_provider

      call augment_lr_gf_endpoint(radial_left, radial_right, left_site, right_site, z, coefficient_gf, forward, &
         capabilities, hgamma_block, hamiltonian_provider, reverse_coefficient_gf, reverse_hgamma_block)
      if (present(reverse_hgamma_block)) then
         call augment_lr_gf_endpoint(radial_right, radial_left, right_site, left_site, z, reverse_coefficient_gf, &
            reverse, capabilities, reverse_hgamma_block, reverse_coefficient_gf=coefficient_gf)
      else
         if (.not. present(hamiltonian_provider)) then
            error stop 'augment_lr_gf_endpoint_pair: reverse hgamma block is required when no provider is supplied'
         end if
         call augment_lr_gf_endpoint(radial_right, radial_left, right_site, left_site, z, reverse_coefficient_gf, &
            reverse, capabilities, hamiltonian_provider=hamiltonian_provider, reverse_coefficient_gf=coefficient_gf)
      end if
   end subroutine augment_lr_gf_endpoint_pair

   subroutine validate_radial_pair(radial_left, radial_right, coefficient_gf)
      type(lmto_radial_basis), intent(in) :: radial_left, radial_right
      complex(rp), intent(in) :: coefficient_gf(:, :)
      integer :: nlocal

      call radial_left%require_supported('ham_only', 'second', .true., .true., .false., .false.)
      call radial_right%require_supported('ham_only', 'second', .true., .true., .false., .false.)
      if (radial_left%nspin /= 2 .or. radial_right%nspin /= 2 .or. radial_left%lmax /= radial_right%lmax) then
         error stop 'augment_lr_gf_endpoint: endpoint radial bases are not a common Pauli sp/spd basis'
      end if
      if (.not. allocated(radial_left%channel_present) .or. .not. allocated(radial_right%channel_present)) then
         error stop 'augment_lr_gf_endpoint: endpoint radial channels are incomplete'
      end if
      if (.not. all(radial_left%channel_present) .or. .not. all(radial_right%channel_present)) then
         error stop 'augment_lr_gf_endpoint: endpoint radial channels are incomplete'
      end if
      nlocal = 2*(radial_left%lmax + 1)**2
      if (any(shape(coefficient_gf) /= [nlocal, nlocal])) then
         error stop 'augment_lr_gf_endpoint: coefficient GF block shape mismatch'
      end if
      if (.not. allocated(radial_left%rofi) .or. .not. allocated(radial_right%rofi)) then
         error stop 'augment_lr_gf_endpoint: endpoint radial meshes are not allocated'
      end if
      if (.not. allocated(radial_left%phi_large) .or. .not. allocated(radial_left%phidot_large) .or. &
          .not. allocated(radial_right%phi_large) .or. .not. allocated(radial_right%phidot_large) .or. &
          .not. allocated(radial_left%enu_work) .or. .not. allocated(radial_right%enu_work)) then
         error stop 'augment_lr_gf_endpoint: endpoint radial augmentation arrays are incomplete'
      end if
   end subroutine validate_radial_pair

   function diagonal_enu(radial) result(enu)
      type(lmto_radial_basis), intent(in) :: radial
      complex(rp) :: enu(2*(radial%lmax + 1)**2, 2*(radial%lmax + 1)**2)
      integer :: iorb, ispin, l, index

      enu = cmplx(0.0_rp, 0.0_rp, rp)
      do ispin = 1, 2
         do iorb = 1, (radial%lmax + 1)**2
            l = lmto_orbital_l(iorb)
            index = (ispin - 1)*(radial%lmax + 1)**2 + iorb
            enu(index, index) = cmplx(radial%enu_work(l + 1, ispin), 0.0_rp, rp)
         end do
      end do
   end function diagonal_enu

   function diagonal_difference(radial, z) result(difference)
      type(lmto_radial_basis), intent(in) :: radial
      complex(rp), intent(in) :: z
      complex(rp) :: difference(2*(radial%lmax + 1)**2, 2*(radial%lmax + 1)**2)
      integer :: i

      difference = -diagonal_enu(radial)
      do i = 1, size(difference, 1)
         difference(i, i) = z + difference(i, i)
      end do
   end function diagonal_difference

   subroutine copy_radial_maps(radial_left, radial_right, augmented)
      type(lmto_radial_basis), intent(in) :: radial_left, radial_right
      type(lr_gf_augmented_block), intent(inout) :: augmented
      integer :: iorb, ispin, l

      allocate(augmented%phi_left(radial_left%npoint, augmented%norb, 2), &
               augmented%phidot_left(radial_left%npoint, augmented%norb, 2), &
               augmented%phi_right(radial_right%npoint, augmented%norb, 2), &
               augmented%phidot_right(radial_right%npoint, augmented%norb, 2), &
               augmented%radius_left(radial_left%npoint), augmented%radius_right(radial_right%npoint))
      augmented%radius_left = radial_left%rofi
      augmented%radius_right = radial_right%rofi
      do ispin = 1, 2
         do iorb = 1, augmented%norb
            l = lmto_orbital_l(iorb)
            augmented%phi_left(:, iorb, ispin) = radial_left%phi_large(:, l + 1, ispin)
            augmented%phidot_left(:, iorb, ispin) = radial_left%phidot_large(:, l + 1, ispin)
            augmented%phi_right(:, iorb, ispin) = radial_right%phi_large(:, l + 1, ispin)
            augmented%phidot_right(:, iorb, ispin) = radial_right%phidot_large(:, l + 1, ispin)
         end do
      end do
   end subroutine copy_radial_maps

   !> Evaluate the full 2x2 Pauli spin GF at two radial/angular endpoints.
   subroutine lr_gf_augmented_block_evaluate_point(this, left_radial_index, left_theta, left_phi, &
                                                   right_radial_index, right_theta, right_phi, point_gf)
      class(lr_gf_augmented_block), intent(in) :: this
      integer, intent(in) :: left_radial_index, right_radial_index
      real(rp), intent(in) :: left_theta, left_phi, right_theta, right_phi
      complex(rp), intent(out) :: point_gf(:, :)
      complex(rp) :: branches(2, 2, 4)

      if (size(point_gf, 1) /= 2 .or. size(point_gf, 2) /= 2) then
         error stop 'lr_gf_augmented_block%evaluate_point: output must be 2x2'
      end if
      call this%branch_at_point(left_radial_index, left_theta, left_phi, right_radial_index, right_theta, right_phi, branches)
      point_gf = sum(branches, dim=3)
   end subroutine lr_gf_augmented_block_evaluate_point

   !> Return the four radial/angular branch matrices separately.  This is the
   !> diagnostic form used by the dense and negative-control oracles.
   subroutine lr_gf_augmented_block_branch_at_point(this, left_radial_index, left_theta, left_phi, &
                                                    right_radial_index, right_theta, right_phi, branches)
      class(lr_gf_augmented_block), intent(in) :: this
      integer, intent(in) :: left_radial_index, right_radial_index
      real(rp), intent(in) :: left_theta, left_phi, right_theta, right_phi
      complex(rp), intent(out) :: branches(:, :, :)
      complex(rp) :: phi_left(2, 2*this%norb), dot_left(2, 2*this%norb)
      complex(rp) :: phi_right(2, 2*this%norb), dot_right(2, 2*this%norb)
      integer :: iorb, ispin, l, index

      if (size(branches, 1) /= 2 .or. size(branches, 2) /= 2 .or. size(branches, 3) /= 4) then
         error stop 'lr_gf_augmented_block%branch_at_point: output must be 2x2x4'
      end if
      call validate_point_indices(this, left_radial_index, right_radial_index)
      phi_left = cmplx(0.0_rp, 0.0_rp, rp)
      dot_left = cmplx(0.0_rp, 0.0_rp, rp)
      phi_right = cmplx(0.0_rp, 0.0_rp, rp)
      dot_right = cmplx(0.0_rp, 0.0_rp, rp)
      do ispin = 1, 2
         do iorb = 1, this%norb
            l = lmto_orbital_l(iorb)
            index = (ispin - 1)*this%norb + iorb
            phi_left(ispin, index) = cmplx(radial_over_radius(this%phi_left(:, iorb, ispin), &
               this%radius_left, left_radial_index, l), 0.0_rp, rp)* &
               response_harmonic(l, orbital_m(iorb, l), left_theta, left_phi)
            dot_left(ispin, index) = cmplx(radial_over_radius(this%phidot_left(:, iorb, ispin), &
               this%radius_left, left_radial_index, l), 0.0_rp, rp)* &
               response_harmonic(l, orbital_m(iorb, l), left_theta, left_phi)
            phi_right(ispin, index) = cmplx(radial_over_radius(this%phi_right(:, iorb, ispin), &
               this%radius_right, right_radial_index, l), 0.0_rp, rp)* &
               response_harmonic(l, orbital_m(iorb, l), right_theta, right_phi)
            dot_right(ispin, index) = cmplx(radial_over_radius(this%phidot_right(:, iorb, ispin), &
               this%radius_right, right_radial_index, l), 0.0_rp, rp)* &
               response_harmonic(l, orbital_m(iorb, l), right_theta, right_phi)
         end do
      end do

      branches(:, :, 1) = matmul(phi_left, matmul(this%coefficient_gf, conjg(transpose(phi_right))))
      branches(:, :, 2) = matmul(dot_left, matmul(this%h_g, conjg(transpose(phi_right))))
      branches(:, :, 3) = matmul(phi_left, matmul(this%g_h, conjg(transpose(dot_right))))
      branches(:, :, 4) = matmul(dot_left, matmul(this%h_g_h, conjg(transpose(dot_right))))
   end subroutine lr_gf_augmented_block_branch_at_point

   subroutine validate_point_indices(this, left_radial_index, right_radial_index)
      class(lr_gf_augmented_block), intent(in) :: this
      integer, intent(in) :: left_radial_index, right_radial_index

      if (.not. allocated(this%phi_left) .or. .not. allocated(this%phi_right) .or. &
          left_radial_index < 1 .or. left_radial_index > size(this%radius_left) .or. &
          right_radial_index < 1 .or. right_radial_index > size(this%radius_right)) then
         error stop 'lr_gf_augmented_block: invalid or uninitialized radial endpoint'
      end if
   end subroutine validate_point_indices

   function radial_over_radius(values, radius, index, l) result(value)
      real(rp), intent(in) :: values(:), radius(:)
      integer, intent(in) :: index, l
      real(rp) :: value
      real(rp) :: first, second

      if (index /= 1) then
         if (radius(index) <= tiny(1.0_rp)) error stop 'lr_gf endpoint: nonpositive radial mesh point'
         value = values(index)/radius(index)
         return
      end if
      if (l /= 0) then
         value = 0.0_rp
      else
         if (size(radius) < 3 .or. radius(2) <= 0.0_rp .or. radius(3) <= radius(2)) then
            error stop 'lr_gf endpoint: cannot form the regular s-channel origin limit'
         end if
         first = values(2)/radius(2)
         second = values(3)/radius(3)
         value = (first*radius(3) - second*radius(2))/(radius(3) - radius(2))
      end if
   end function radial_over_radius

   pure integer function orbital_m(iorb, l) result(m)
      integer, intent(in) :: iorb, l
      m = iorb - l*l - l - 1
   end function orbital_m

end module lr_gf_endpoint_augmentation_mod
