!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Site-resolved charge-spin response channels and transition vertices.
!>
!> The reciprocal Hamiltonian packs a unit-cell spinor as site-major blocks,
!> with each block ordered `(orbital up..., orbital down...)`.  A response
!> channel is consequently a site projector times one global Pauli component.
!> The optional `l_sector` and `orbital` selectors retain the information
!> required for future resolved response while the normal production use is a
!> site- and orbital-summed channel.
module response_vertices_mod
   use precision_mod, only: rp
   use math_mod, only: i_unit
   use response_components_mod, only: RESPONSE_CHARGE, RESPONSE_MX, RESPONSE_MY, RESPONSE_MZ, &
      RESPONSE_PLUS, RESPONSE_MINUS, is_response_component
   implicit none

   private

   !> Sentinel for a channel which sums every l sector or local orbital.
   integer, parameter, public :: RESPONSE_UNRESOLVED = -1

   !> A response channel P_site sigma^component.  `site` and `orbital` use
   !> one-based indices.  `orbital` is the local orbital index within a site,
   !> independent of spin.  A non-negative l_sector follows l=0,1,2,3.
   type, public :: response_channel
      integer :: site = 0
      integer :: component = RESPONSE_CHARGE
      integer :: l_sector = RESPONSE_UNRESOLVED
      integer :: orbital = RESPONSE_UNRESOLVED
   end type response_channel

   public :: site_projected_operator
   public :: response_transition_vertex
   public :: response_transition_vectors
   public :: build_site_charge_spin_channels

contains

   !> Construct the explicit dense P_site sigma operator for diagnostics and
   !> small reference calculations.  Production accumulation should use
   !> response_transition_vertex(s), which streams a response-space vector and
   !> never materializes a transition-pair response tensor.
   function site_projected_operator(channel, site_orbital_counts) result(op)
      type(response_channel), intent(in) :: channel
      integer, intent(in) :: site_orbital_counts(:)
      complex(rp), allocatable :: op(:, :)
      integer :: nmat, offset, nlocal, orbital

      call validate_channel(channel, site_orbital_counts)
      nmat = 2*sum(site_orbital_counts)
      allocate(op(nmat, nmat))
      op = cmplx(0.0_rp, 0.0_rp, rp)

      offset = site_block_offset(channel%site, site_orbital_counts)
      nlocal = site_orbital_counts(channel%site)
      do orbital = 1, nlocal
         if (.not. channel_selects_orbital(channel, orbital)) cycle
         call add_component_element(op, offset + orbital, offset + nlocal + orbital, channel%component)
      end do
   end function site_projected_operator

   !> Reference transition vertex v_A = <bra|P_site sigma^A|ket>.  The direct
   !> contraction accepts arbitrary complex spinors, including SOC and
   !> non-collinear eigenstates.
   function response_transition_vertex(channel, site_orbital_counts, bra_spinor, ket_spinor) result(vertex)
      type(response_channel), intent(in) :: channel
      integer, intent(in) :: site_orbital_counts(:)
      complex(rp), intent(in) :: bra_spinor(:), ket_spinor(:)
      complex(rp) :: vertex
      integer :: nmat, offset, nlocal, orbital, iup, idown

      call validate_channel(channel, site_orbital_counts)
      nmat = 2*sum(site_orbital_counts)
      if (size(bra_spinor) /= nmat .or. size(ket_spinor) /= nmat) then
         error stop 'response_transition_vertex: spinor size does not match site layout'
      end if

      vertex = cmplx(0.0_rp, 0.0_rp, rp)
      offset = site_block_offset(channel%site, site_orbital_counts)
      nlocal = site_orbital_counts(channel%site)
      do orbital = 1, nlocal
         if (.not. channel_selects_orbital(channel, orbital)) cycle
         iup = offset + orbital
         idown = offset + nlocal + orbital
         select case (channel%component)
         case (RESPONSE_CHARGE)
            vertex = vertex + conjg(bra_spinor(iup))*ket_spinor(iup) + &
               conjg(bra_spinor(idown))*ket_spinor(idown)
         case (RESPONSE_MX)
            vertex = vertex + conjg(bra_spinor(iup))*ket_spinor(idown) + &
               conjg(bra_spinor(idown))*ket_spinor(iup)
         case (RESPONSE_MY)
            vertex = vertex - i_unit*conjg(bra_spinor(iup))*ket_spinor(idown) + &
               i_unit*conjg(bra_spinor(idown))*ket_spinor(iup)
         case (RESPONSE_MZ)
            vertex = vertex + conjg(bra_spinor(iup))*ket_spinor(iup) - &
               conjg(bra_spinor(idown))*ket_spinor(idown)
         case (RESPONSE_PLUS)
            vertex = vertex + 2.0_rp*conjg(bra_spinor(iup))*ket_spinor(idown)
         case (RESPONSE_MINUS)
            vertex = vertex + 2.0_rp*conjg(bra_spinor(idown))*ket_spinor(iup)
         end select
      end do
   end function response_transition_vertex

   !> Batched response-space transition vectors.  Each input column is one
   !> `(n,m,k,q)` transition; `vertices(:,ipair)` is its v_A.  The caller can
   !> immediately accumulate w*v_A*conjg(v_B), so no `(k,n,m,A,B)` object is
   !> stored here.
   subroutine response_transition_vectors(channels, site_orbital_counts, bra_spinors, ket_spinors, vertices)
      type(response_channel), intent(in) :: channels(:)
      integer, intent(in) :: site_orbital_counts(:)
      complex(rp), intent(in) :: bra_spinors(:, :), ket_spinors(:, :)
      complex(rp), intent(out) :: vertices(:, :)
      integer :: ichannel, ipair, nmat

      nmat = 2*sum(site_orbital_counts)
      if (size(bra_spinors, 1) /= nmat .or. size(ket_spinors, 1) /= nmat .or. &
          size(bra_spinors, 2) /= size(ket_spinors, 2)) then
         error stop 'response_transition_vectors: incompatible spinor batch shape'
      end if
      if (size(vertices, 1) /= size(channels) .or. size(vertices, 2) /= size(bra_spinors, 2)) then
         error stop 'response_transition_vectors: vertices must have shape (nchannels,npairs)'
      end if

      do ipair = 1, size(bra_spinors, 2)
         do ichannel = 1, size(channels)
            vertices(ichannel, ipair) = response_transition_vertex(channels(ichannel), site_orbital_counts, &
               bra_spinors(:, ipair), ket_spinors(:, ipair))
         end do
      end do
   end subroutine response_transition_vectors

   !> Canonical active basis for a full charge-spin response: site-major with
   !> components (charge, mx, my, mz) at every site.  Circular labels remain
   !> views for the legacy transverse output and are intentionally excluded.
   subroutine build_site_charge_spin_channels(nsite, channels)
      integer, intent(in) :: nsite
      type(response_channel), allocatable, intent(out) :: channels(:)
      integer :: isite, component

      if (nsite < 1) error stop 'build_site_charge_spin_channels: nsite must be positive'
      allocate(channels(4*nsite))
      do isite = 1, nsite
         do component = RESPONSE_CHARGE, RESPONSE_MZ
            channels(4*(isite-1) + component + 1) = response_channel(isite, component)
         end do
      end do
   end subroutine build_site_charge_spin_channels

   subroutine add_component_element(op, iup, idown, component)
      complex(rp), intent(inout) :: op(:, :)
      integer, intent(in) :: iup, idown, component

      select case (component)
      case (RESPONSE_CHARGE)
         op(iup, iup) = cmplx(1.0_rp, 0.0_rp, rp)
         op(idown, idown) = cmplx(1.0_rp, 0.0_rp, rp)
      case (RESPONSE_MX)
         op(iup, idown) = cmplx(1.0_rp, 0.0_rp, rp)
         op(idown, iup) = cmplx(1.0_rp, 0.0_rp, rp)
      case (RESPONSE_MY)
         op(iup, idown) = -i_unit
         op(idown, iup) = i_unit
      case (RESPONSE_MZ)
         op(iup, iup) = cmplx(1.0_rp, 0.0_rp, rp)
         op(idown, idown) = cmplx(-1.0_rp, 0.0_rp, rp)
      case (RESPONSE_PLUS)
         op(iup, idown) = cmplx(2.0_rp, 0.0_rp, rp)
      case (RESPONSE_MINUS)
         op(idown, iup) = cmplx(2.0_rp, 0.0_rp, rp)
      end select
   end subroutine add_component_element

   pure logical function channel_selects_orbital(channel, orbital)
      type(response_channel), intent(in) :: channel
      integer, intent(in) :: orbital
      integer :: first_orbital, last_orbital

      channel_selects_orbital = .true.
      if (channel%l_sector /= RESPONSE_UNRESOLVED) then
         call l_sector_bounds(channel%l_sector, first_orbital, last_orbital)
         channel_selects_orbital = orbital >= first_orbital .and. orbital <= last_orbital
      end if
      if (channel%orbital /= RESPONSE_UNRESOLVED) then
         channel_selects_orbital = channel_selects_orbital .and. orbital == channel%orbital
      end if
   end function channel_selects_orbital

   pure integer function site_block_offset(site, site_orbital_counts)
      integer, intent(in) :: site, site_orbital_counts(:)
      integer :: isite

      site_block_offset = 0
      do isite = 1, site - 1
         site_block_offset = site_block_offset + 2*site_orbital_counts(isite)
      end do
   end function site_block_offset

   subroutine validate_channel(channel, site_orbital_counts)
      type(response_channel), intent(in) :: channel
      integer, intent(in) :: site_orbital_counts(:)
      integer :: first_orbital, last_orbital

      if (size(site_orbital_counts) == 0 .or. any(site_orbital_counts <= 0)) then
         error stop 'response vertices: site orbital counts must be positive'
      end if
      if (channel%site < 1 .or. channel%site > size(site_orbital_counts)) then
         error stop 'response vertices: channel site is outside the site layout'
      end if
      if (.not. is_response_component(channel%component)) then
         error stop 'response vertices: invalid response component'
      end if
      if (channel%l_sector /= RESPONSE_UNRESOLVED) then
         call l_sector_bounds(channel%l_sector, first_orbital, last_orbital)
         if (last_orbital > site_orbital_counts(channel%site)) then
            error stop 'response vertices: requested l sector is absent from this site basis'
         end if
      end if
      if (channel%orbital /= RESPONSE_UNRESOLVED) then
         if (channel%orbital < 1 .or. channel%orbital > site_orbital_counts(channel%site)) then
            error stop 'response vertices: requested local orbital is outside the site basis'
         end if
         if (channel%l_sector /= RESPONSE_UNRESOLVED) then
            call l_sector_bounds(channel%l_sector, first_orbital, last_orbital)
            if (channel%orbital < first_orbital .or. channel%orbital > last_orbital) then
               error stop 'response vertices: orbital is not in the requested l sector'
            end if
         end if
      end if
   end subroutine validate_channel

   pure subroutine l_sector_bounds(l_sector, first_orbital, last_orbital)
      integer, intent(in) :: l_sector
      integer, intent(out) :: first_orbital, last_orbital

      select case (l_sector)
      case (0)
         first_orbital = 1
         last_orbital = 1
      case (1)
         first_orbital = 2
         last_orbital = 4
      case (2)
         first_orbital = 5
         last_orbital = 9
      case (3)
         first_orbital = 10
         last_orbital = 16
      case default
         error stop 'response vertices: l sector must be 0 through 3'
      end select
   end subroutine l_sector_bounds

end module response_vertices_mod
