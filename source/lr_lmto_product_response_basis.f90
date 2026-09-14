!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Production finite LMTO product-space representation.
!>
!> The representation is a numerical coordinate layer only.  It owns the
!> scaled radial product SVD and the forward transition map; it does not know
!> about occupations, frequency denominators, susceptibility accumulation,
!> XC kernels, or Dyson equations.
!------------------------------------------------------------------------------
module lr_lmto_product_response_basis_mod

   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use response_angular_basis_mod, only: response_gaunt
   use lr_response_space_mod, only: response_space_layout
   use lr_pauli_transition_vertex_mod, only: pauli_endpoint_state
   implicit none
   private

   integer, parameter, public :: lmto_product_channel_plus = 1
   integer, parameter, public :: lmto_product_channel_minus = 2

   type, public :: lmto_product_candidate
      integer :: l = 0
      integer :: lp = 0
      integer :: p = 0
      integer :: q = 0
   end type lmto_product_candidate

   type, public :: lmto_product_block
      integer :: site = 0
      integer :: response_l = 0
      integer :: channel = 0
      integer :: npoint = 0
      integer :: ncandidate = 0
      integer :: nsv = 0
      integer :: rank = 0
      integer :: rank_tau1 = 0
      integer :: rank_tau10 = 0
      integer :: rank_tau100 = 0
      real(rp) :: tau1 = 0.0_rp
      real(rp) :: tau10 = 0.0_rp
      real(rp) :: tau100 = 0.0_rp
      logical :: rank_stable = .false.
      type(lmto_product_candidate), allocatable :: candidates(:)
      real(rp), allocatable :: column_norms(:)
      real(rp), allocatable :: singular_values(:)
      !> Weighted orthonormal radial modes U, with columns retained by rank.
      complex(rp), allocatable :: weighted_modes(:, :)
      !> Retained rows of V^H from the direct SVD, kept for reconstruction
      !> audits and future operator projection.
      complex(rp), allocatable :: right_modes(:, :)
      !> Forward candidate-to-coordinate map F = Sigma V^H D.
      complex(rp), allocatable :: forward_transform(:, :)
   end type lmto_product_block

   type, public :: lmto_product_response_basis
      integer :: nsite = 0
      integer :: orbital_lmax = -1
      integer :: response_lmax = -1
      integer :: circular_channel = 0
      integer :: npoint = 0
      integer :: unpruned_dimension = 0
      integer :: product_dimension = 0
      type(lmto_product_block), allocatable :: blocks(:, :)
   contains
      procedure :: initialize => lmto_product_response_basis_initialize
      procedure :: flat_index => lmto_product_flat_index
      procedure :: unflatten_index => lmto_product_unflatten_index
      procedure :: candidate_coefficients => lmto_product_candidate_coefficients
      procedure :: transition_coordinates => lmto_product_transition_coordinates
   end type lmto_product_response_basis

   public :: lmto_product_candidate_count
   public :: lmto_enumerate_product_candidates

   interface
      subroutine zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
         import :: rp
         character(len=1), intent(in) :: jobu, jobvt
         integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
         complex(rp), intent(inout) :: a(lda, *)
         real(rp), intent(out) :: s(*)
         complex(rp), intent(out) :: u(ldu, *), vt(ldvt, *), work(*)
         real(rp), intent(inout) :: rwork(*)
         integer, intent(out) :: info
      end subroutine zgesvd
   end interface

contains

   pure integer function lmto_product_candidate_count(lmax, response_l) result(count)
      integer, intent(in) :: lmax, response_l
      integer :: l, lp

      count = 0
      do l = 0, lmax
         do lp = 0, lmax
            if (allowed_pair(l, lp, response_l)) count = count + 4
         end do
      end do
   end function lmto_product_candidate_count

   subroutine lmto_enumerate_product_candidates(lmax, response_l, candidates)
      integer, intent(in) :: lmax, response_l
      type(lmto_product_candidate), allocatable, intent(out) :: candidates(:)
      integer :: l, lp, p, q, k

      allocate(candidates(lmto_product_candidate_count(lmax, response_l)))
      k = 0
      do l = 0, lmax
         do lp = 0, lmax
            if (.not. allowed_pair(l, lp, response_l)) cycle
            do p = 0, 1
               do q = 0, 1
                  k = k + 1
                  candidates(k) = lmto_product_candidate(l, lp, p, q)
               end do
            end do
         end do
      end do
   end subroutine lmto_enumerate_product_candidates

   pure logical function allowed_pair(l, lp, response_l) result(is_allowed)
      integer, intent(in) :: l, lp, response_l

      is_allowed = abs(l - lp) <= response_l .and. response_l <= l + lp .and. &
         mod(l + lp + response_l, 2) == 0
   end function allowed_pair

   subroutine lmto_product_response_basis_initialize(this, space, radial_bases, circular_channel, strict_rank)
      class(lmto_product_response_basis), intent(out) :: this
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      integer, intent(in) :: circular_channel
      logical, intent(in), optional :: strict_rank
      logical :: strict
      integer :: isite, response_l

      strict = .true.
      if (present(strict_rank)) strict = strict_rank
      call validate_initialization_inputs(space, radial_bases, circular_channel)

      this%nsite = space%nsite
      this%orbital_lmax = radial_bases(1)%lmax
      this%response_lmax = space%response_lmax
      this%circular_channel = circular_channel
      this%npoint = space%npoint
      this%unpruned_dimension = 0
      this%product_dimension = 0
      allocate(this%blocks(this%nsite, 0:this%response_lmax))

      do isite = 1, this%nsite
         do response_l = 0, this%response_lmax
            this%unpruned_dimension = this%unpruned_dimension + &
               (2*response_l + 1)*lmto_product_candidate_count(this%orbital_lmax, response_l)
            call build_product_block(this%blocks(isite, response_l), space, radial_bases(isite), isite, &
               response_l, circular_channel, strict)
            this%product_dimension = this%product_dimension + &
               (2*response_l + 1)*this%blocks(isite, response_l)%rank
         end do
      end do
   end subroutine lmto_product_response_basis_initialize

   subroutine validate_initialization_inputs(space, radial_bases, circular_channel)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      integer, intent(in) :: circular_channel
      integer :: isite, lmax
      real(rp), parameter :: mesh_tolerance = 2.0e-13_rp

      if (space%nsite < 1 .or. space%response_lmax < 0 .or. space%nchannel /= 1 .or. &
          .not. allocated(space%radius) .or. .not. allocated(space%radial_weights)) then
         error stop 'lmto_product_response_basis: only one explicit response channel is supported'
      end if
      if (size(radial_bases) /= space%nsite) then
         error stop 'lmto_product_response_basis: one radial basis is required per site'
      end if
      if (circular_channel /= lmto_product_channel_plus .and. circular_channel /= lmto_product_channel_minus) then
         error stop 'lmto_product_response_basis: invalid circular channel'
      end if
      lmax = radial_bases(1)%lmax
      if (lmax < 0 .or. lmax > 2 .or. space%response_lmax /= 2*lmax) then
         error stop 'lmto_product_response_basis: unsupported or incomplete angular product space'
      end if
      if (space%npoint < 3 .or. size(space%radius) /= space%npoint .or. &
          size(space%radial_weights) /= space%npoint) then
         error stop 'lmto_product_response_basis: inconsistent radial layout'
      end if
      do isite = 1, space%nsite
         call radial_bases(isite)%require_supported('ham_only', 'second', .true., .true., .false., .false.)
         if (radial_bases(isite)%lmax /= lmax .or. radial_bases(isite)%nspin /= 2 .or. &
             radial_bases(isite)%npoint /= space%npoint .or. .not. allocated(radial_bases(isite)%rofi) .or. &
             .not. allocated(radial_bases(isite)%channel_present) .or. .not. all(radial_bases(isite)%channel_present)) then
            error stop 'lmto_product_response_basis: incomplete radial basis'
         end if
         if (abs(radial_bases(isite)%mesh_a - space%a) > mesh_tolerance*max(1.0_rp, abs(space%a)) .or. &
             abs(radial_bases(isite)%mesh_b - space%b) > mesh_tolerance*max(1.0_rp, abs(space%b)) .or. &
             maxval(abs(radial_bases(isite)%rofi - space%radius)) > mesh_tolerance* &
                max(1.0_rp, maxval(abs(space%radius)))) then
            error stop 'lmto_product_response_basis: radial mesh provenance mismatch'
         end if
      end do
      if (any(.not. ieee_is_finite(space%radial_weights)) .or. any(space%radial_weights < 0.0_rp)) then
         error stop 'lmto_product_response_basis: invalid LR-04 radial weights'
      end if
   end subroutine validate_initialization_inputs

   subroutine build_product_block(block, space, radial, site, response_l, circular_channel, strict)
      type(lmto_product_block), intent(out) :: block
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: site, response_l, circular_channel
      logical, intent(in) :: strict
      type(lmto_product_candidate), allocatable :: candidates(:)
      complex(rp), allocatable :: a(:, :), u(:, :), vt(:, :), work(:)
      complex(rp) :: work_query(1)
      real(rp), allocatable :: column_norms(:), singular_values(:), rwork(:)
      real(rp) :: tau1, tau10, tau100
      integer :: nr, ncandidate, nsv, lwork, info, ir, k, rank1, rank10, rank100

      nr = space%npoint
      call lmto_enumerate_product_candidates(radial%lmax, response_l, candidates)
      ncandidate = size(candidates)
      nsv = min(nr, ncandidate)
      allocate(a(nr, ncandidate), column_norms(ncandidate), singular_values(nsv), &
         u(nr, nsv), vt(nsv, ncandidate), rwork(max(1, 5*nsv)))
      a = cmplx(0.0_rp, 0.0_rp, rp)
      do k = 1, ncandidate
         do ir = 1, nr
            a(ir, k) = cmplx(radial_product(radial, ir, candidates(k)%l, candidates(k)%lp, &
               circular_channel, candidates(k)%p, candidates(k)%q), 0.0_rp, rp)
         end do
         column_norms(k) = sqrt(sum(space%radial_weights*real(a(:, k)*conjg(a(:, k)), rp)))
         if (.not. ieee_is_finite(column_norms(k)) .or. column_norms(k) <= tiny(1.0_rp)) then
            error stop 'lmto_product_response_basis: zero or non-finite candidate norm'
         end if
         a(:, k) = sqrt(space%radial_weights)*a(:, k)/column_norms(k)
      end do

      ! Direct SVD of W^(1/2) Btilde.  The candidate Gram matrix is not formed.
      call zgesvd('S', 'S', nr, ncandidate, a, nr, singular_values, u, nr, vt, nsv, work_query, -1, rwork, info)
      if (info /= 0) error stop 'lmto_product_response_basis: zgesvd workspace query failed'
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      call zgesvd('S', 'S', nr, ncandidate, a, nr, singular_values, u, nr, vt, nsv, work, lwork, rwork, info)
      deallocate(work)
      if (info /= 0) error stop 'lmto_product_response_basis: zgesvd failed'
      if (.not. all(ieee_is_finite(singular_values))) then
         error stop 'lmto_product_response_basis: non-finite singular value'
      end if

      tau1 = real(max(nr, ncandidate), rp)*epsilon(1.0_rp)*singular_values(1)
      tau10 = 10.0_rp*tau1
      tau100 = 100.0_rp*tau1
      rank1 = count(singular_values > tau1)
      rank10 = count(singular_values > tau10)
      rank100 = count(singular_values > tau100)
      if (rank1 < 1) error stop 'lmto_product_response_basis: no retained numerical product mode'
      if (strict .and. (rank1 /= rank10 .or. rank1 /= rank100)) then
         write (*, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'LMTO product rank sensitivity site=', site, ' L=', response_l, &
            ' tau1/tau10/tau100=', rank1, '/', rank10, '/', rank100
         error stop 'lmto_product_response_basis: rank sensitivity requires review'
      end if

      block%site = site
      block%response_l = response_l
      block%channel = circular_channel
      block%npoint = nr
      block%ncandidate = ncandidate
      block%nsv = nsv
      block%rank = rank1
      block%rank_tau1 = rank1
      block%rank_tau10 = rank10
      block%rank_tau100 = rank100
      block%tau1 = tau1
      block%tau10 = tau10
      block%tau100 = tau100
      block%rank_stable = rank1 == rank10 .and. rank1 == rank100
      allocate(block%candidates(ncandidate), block%column_norms(ncandidate), block%singular_values(nsv), &
         block%weighted_modes(nr, rank1), block%right_modes(rank1, ncandidate), block%forward_transform(rank1, ncandidate))
      block%candidates = candidates
      block%column_norms = column_norms
      block%singular_values = singular_values
      block%weighted_modes = u(:, 1:rank1)
      block%right_modes = vt(1:rank1, :)
      do k = 1, ncandidate
         block%forward_transform(:, k) = singular_values(1:rank1)*vt(1:rank1, k)*column_norms(k)
      end do
      deallocate(candidates, a, u, vt, column_norms, singular_values, rwork)
   end subroutine build_product_block

   integer function lmto_product_flat_index(this, site, response_l, response_m, product_mode) result(flat)
      class(lmto_product_response_basis), intent(in) :: this
      integer, intent(in) :: site, response_l, response_m, product_mode
      integer :: isite, l, rank

      call validate_product_index_inputs(this, site, response_l, response_m, product_mode)
      flat = 0
      do isite = 1, site - 1
         do l = 0, this%response_lmax
            flat = flat + (2*l + 1)*this%blocks(isite, l)%rank
         end do
      end do
      do l = 0, response_l - 1
         flat = flat + (2*l + 1)*this%blocks(site, l)%rank
      end do
      rank = this%blocks(site, response_l)%rank
      flat = flat + (response_m + response_l)*rank + product_mode
   end function lmto_product_flat_index

   subroutine lmto_product_unflatten_index(this, flat, site, response_l, response_m, product_mode)
      class(lmto_product_response_basis), intent(in) :: this
      integer, intent(in) :: flat
      integer, intent(out) :: site, response_l, response_m, product_mode
      integer :: remaining, block_size, isite, l

      if (flat < 1 .or. flat > this%product_dimension) then
         error stop 'lmto_product_response_basis: product flat index outside representation'
      end if
      remaining = flat
      do isite = 1, this%nsite
         do l = 0, this%response_lmax
            block_size = (2*l + 1)*this%blocks(isite, l)%rank
            if (remaining <= block_size) then
               site = isite
               response_l = l
               response_m = (remaining - 1)/this%blocks(isite, l)%rank - l
               product_mode = mod(remaining - 1, this%blocks(isite, l)%rank) + 1
               return
            end if
            remaining = remaining - block_size
         end do
      end do
      error stop 'lmto_product_response_basis: product unflatten failed'
   end subroutine lmto_product_unflatten_index

   subroutine validate_product_index_inputs(this, site, response_l, response_m, product_mode)
      class(lmto_product_response_basis), intent(in) :: this
      integer, intent(in) :: site, response_l, response_m, product_mode

      if (.not. allocated(this%blocks) .or. site < 1 .or. site > this%nsite .or. &
          response_l < 0 .or. response_l > this%response_lmax .or. abs(response_m) > response_l .or. &
          product_mode < 1 .or. product_mode > this%blocks(site, response_l)%rank) then
         error stop 'lmto_product_response_basis: product index outside representation'
      end if
   end subroutine validate_product_index_inputs

   subroutine lmto_product_candidate_coefficients(this, site, response_l, response_m, left_state, right_state, coefficients)
      class(lmto_product_response_basis), intent(in) :: this
      integer, intent(in) :: site, response_l, response_m
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(out) :: coefficients(:)
      integer :: norb, offset, iorb, jorb, spin_left, spin_right, k, ncandidate
      integer :: orbital_l, orbital_lp, orbital_m, orbital_mp
      real(rp) :: left_power, right_power

      if (.not. allocated(this%blocks)) error stop 'lmto_product_response_basis: representation is uninitialized'
      if (site < 1 .or. site > this%nsite .or. response_l < 0 .or. response_l > this%response_lmax .or. &
          abs(response_m) > response_l) error stop 'lmto_product_response_basis: invalid candidate block index'
      ncandidate = this%blocks(site, response_l)%ncandidate
      if (size(coefficients) /= ncandidate) then
         error stop 'lmto_product_response_basis: candidate coefficient shape mismatch'
      end if
      norb = (this%orbital_lmax + 1)**2
      if (.not. allocated(left_state%coefficients) .or. .not. allocated(right_state%coefficients) .or. &
          size(left_state%coefficients) /= 2*norb*this%nsite .or. size(right_state%coefficients) /= 2*norb*this%nsite) then
         error stop 'lmto_product_response_basis: endpoint coefficient shape mismatch'
      end if
      offset = (site - 1)*2*norb
      if (this%circular_channel == lmto_product_channel_plus) then
         spin_left = 1
         spin_right = 2
      else
         spin_left = 2
         spin_right = 1
      end if
      coefficients = cmplx(0.0_rp, 0.0_rp, rp)
      do k = 1, ncandidate
         left_power = 1.0_rp
         right_power = 1.0_rp
         if (this%blocks(site, response_l)%candidates(k)%p == 1) left_power = left_state%energy
         if (this%blocks(site, response_l)%candidates(k)%q == 1) right_power = right_state%energy
         do iorb = 1, norb
            orbital_l = lmto_orbital_l(iorb)
            if (orbital_l /= this%blocks(site, response_l)%candidates(k)%l) cycle
            orbital_m = iorb - orbital_l*orbital_l - orbital_l - 1
            do jorb = 1, norb
               orbital_lp = lmto_orbital_l(jorb)
               if (orbital_lp /= this%blocks(site, response_l)%candidates(k)%lp) cycle
               orbital_mp = jorb - orbital_lp*orbital_lp - orbital_lp - 1
               coefficients(k) = coefficients(k) + left_power*right_power* &
                  conjg(left_state%coefficients(offset + (spin_left - 1)*norb + iorb))* &
                  right_state%coefficients(offset + (spin_right - 1)*norb + jorb)* &
                  response_gaunt(orbital_l, orbital_m, orbital_lp, orbital_mp, response_l, response_m)
            end do
         end do
      end do
   end subroutine lmto_product_candidate_coefficients

   subroutine lmto_product_transition_coordinates(this, left_state, right_state, coordinates)
      class(lmto_product_response_basis), intent(in) :: this
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(out) :: coordinates(:)
      complex(rp), allocatable :: coefficients(:)
      integer :: site, response_l, response_m

      if (.not. allocated(this%blocks)) error stop 'lmto_product_response_basis: representation is uninitialized'
      if (size(coordinates) /= this%product_dimension) then
         error stop 'lmto_product_response_basis: product coordinate shape mismatch'
      end if
      coordinates = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, this%nsite
         do response_l = 0, this%response_lmax
            allocate(coefficients(this%blocks(site, response_l)%ncandidate))
            do response_m = -response_l, response_l
               call this%candidate_coefficients(site, response_l, response_m, left_state, right_state, coefficients)
               coordinates(this%flat_index(site, response_l, response_m, 1): &
                  this%flat_index(site, response_l, response_m, this%blocks(site, response_l)%rank)) = &
                  matmul(this%blocks(site, response_l)%forward_transform, coefficients)
            end do
            deallocate(coefficients)
         end do
      end do
   end subroutine lmto_product_transition_coordinates

   pure real(rp) function endpoint_component(radial, ir, l, spin, power) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, spin, power

      if (power == 0) then
         value = radial%phi_large(ir, l + 1, spin) - radial%enu_work(l + 1, spin)* &
            radial%phidot_large(ir, l + 1, spin)
      else
         value = radial%phidot_large(ir, l + 1, spin)
      end if
   end function endpoint_component

   real(rp) function radial_product(radial, ir, l, lp, circular_channel, p, q) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, lp, circular_channel, p, q
      integer :: spin_left, spin_right
      real(rp) :: first, second, rfirst, rsecond

      if (circular_channel == lmto_product_channel_plus) then
         spin_left = 1
         spin_right = 2
      else
         spin_left = 2
         spin_right = 1
      end if
      if (ir /= 1) then
         if (radial%rofi(ir) <= tiny(1.0_rp)) error stop 'lmto_product_response_basis: invalid radial point'
         value = endpoint_component(radial, ir, l, spin_left, p)*endpoint_component(radial, ir, lp, spin_right, q)/ &
            radial%rofi(ir)**2
         return
      end if
      if (l /= 0 .or. lp /= 0) then
         value = 0.0_rp
         return
      end if
      rfirst = radial%rofi(2)
      rsecond = radial%rofi(3)
      first = endpoint_component(radial, 2, l, spin_left, p)*endpoint_component(radial, 2, lp, spin_right, q)/rfirst**2
      second = endpoint_component(radial, 3, l, spin_left, p)*endpoint_component(radial, 3, lp, spin_right, q)/rsecond**2
      value = (first*rsecond**2 - second*rfirst**2)/(rsecond**2 - rfirst**2)
   end function radial_product

end module lr_lmto_product_response_basis_mod
