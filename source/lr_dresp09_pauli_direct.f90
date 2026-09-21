!------------------------------------------------------------------------------
! DRESP-09R independent Pauli projected field operator.
!
! This is an audit oracle.  It intentionally repeats the accepted Pauli
! large-component radial/angular contraction instead of calling the production
! transition evaluator or the compact metric-adjoint implementation.  The
! input field is the contravariant source field: for a compact covariant field
! b, its raw source representative is U conjg(b), not U b.
!------------------------------------------------------------------------------
module lr_dresp09_pauli_direct_mod

   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use response_angular_basis_mod, only: response_gaunt
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex, &
      response_unflatten_superindex
   use lr_response_space_mod, only: response_space_layout
   use lr_lmto_product_response_basis_mod, only: lmto_product_channel_plus, lmto_product_channel_minus
   use lr_lmto_endpoint_branches_mod, only: lmto_product_nbranch, lmto_product_second_order_radial_branch
   implicit none
   private

   public :: dresp09_pauli_direct_source
   public :: dresp09_pauli_direct_source_components

contains

   !> Build the six independent endpoint-power branches of the source
   !> operator.  The source matrix uses the reversed spin/orbital endpoints
   !> of the measurement tensor while retaining its authoritative branch map.
   subroutine dresp09_pauli_direct_source_components(space, radial_bases, source_field, circular_channel, components)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: source_field(:)
      integer, intent(in) :: circular_channel
      complex(rp), intent(out) :: components(:, :, :)

      integer :: norb, nbasis, nsite, npoint, lmax
      integer :: flat, ir, site, response_l, response_m, channel
      integer :: iorb, jorb, orbital_l, orbital_lp, orbital_m, orbital_mp
      integer :: spin_left, spin_right, branch, source_row, source_col
      integer :: left_offset, source_flat
      type(response_super_index) :: item
      real(rp) :: radial_pair, gaunt
      complex(rp) :: coefficient

      call validate_inputs(space, radial_bases, source_field, circular_channel, components)
      nsite = space%nsite
      npoint = space%npoint
      lmax = radial_bases(1)%lmax
      norb = (lmax + 1)**2
      nbasis = 2*norb*nsite
      components = cmplx(0.0_rp, 0.0_rp, rp)

      if (circular_channel == lmto_product_channel_plus) then
         spin_left = 1
         spin_right = 2
      else
         spin_left = 2
         spin_right = 1
      end if

      ! The source field is already the contravariant physical field.  Its
      ! coefficient is therefore not conjugated here.  The conjugation is in
      ! the Pauli/angular adjoint below, exactly as required by O^H.
      do flat = 1, size(source_field)
         call response_unflatten_superindex(flat, nsite, space%response_lmax, npoint, 1, item)
         site = item%site
         response_l = item%response_l
         response_m = item%response_m
         ir = item%radial_point
         channel = item%channel
         coefficient = source_field(flat)
         if (coefficient == cmplx(0.0_rp, 0.0_rp, rp)) cycle
         left_offset = (site - 1)*2*norb

         do iorb = 1, norb
            orbital_l = lmto_orbital_l(iorb)
            orbital_m = iorb - orbital_l*orbital_l - orbital_l - 1
            do jorb = 1, norb
               orbital_lp = lmto_orbital_l(jorb)
               orbital_mp = jorb - orbital_lp*orbital_lp - orbital_lp - 1
               gaunt = response_gaunt(orbital_l, orbital_m, orbital_lp, orbital_mp, response_l, response_m)
               if (gaunt == 0.0_rp) cycle

               source_row = left_offset + (spin_right - 1)*norb + jorb
               source_col = left_offset + (spin_left - 1)*norb + iorb
               do branch = 1, lmto_product_nbranch
                  ! The matrix indices are transposed below, but the branch
                  ! remains attached to the original measurement endpoints.
                  ! Keeping the original (l,spin)_left/right ordering is
                  ! essential for the asymmetric 20 and 02 second-order
                  ! branches.
                  call lmto_product_second_order_radial_branch(radial_bases(site), ir, orbital_l, orbital_lp, &
                     spin_left, spin_right, branch, radial_pair)
                  if (radial_pair == 0.0_rp) cycle
                  call response_flatten_superindex(response_super_index(site, response_l, response_m, ir, channel), &
                     nsite, space%response_lmax, npoint, 1, source_flat)
                  components(source_row, source_col, branch) = &
                     components(source_row, source_col, branch) + coefficient*space%metric_weights(source_flat)* &
                     cmplx(radial_pair*gaunt, 0.0_rp, rp)
               end do
            end do
         end do
      end do
   end subroutine dresp09_pauli_direct_source_components

   !> Contract the independent branch tensor with endpoint Hamiltonian powers:
   !> S = sum_(p,q) H**q S^(p,q) H**p.  This is the same endpoint orientation
   !> as the compact adjoint, but is implemented here without calling it.
   subroutine dresp09_pauli_direct_source(space, radial_bases, source_field, circular_channel, hamiltonian, operator)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: source_field(:), hamiltonian(:, :)
      integer, intent(in) :: circular_channel
      complex(rp), intent(out) :: operator(:, :)
      complex(rp), allocatable :: components(:, :, :), term(:, :)
      integer :: n, p, q, component

      n = size(hamiltonian, 1)
      if (size(hamiltonian, 2) /= n .or. any(shape(operator) /= [n, n])) then
         error stop 'DRESP-09R Pauli direct source: Hamiltonian/operator shape mismatch'
      end if
      allocate(components(n, n, lmto_product_nbranch), term(n, n))
      call dresp09_pauli_direct_source_components(space, radial_bases, source_field, circular_channel, components)
      operator = cmplx(0.0_rp, 0.0_rp, rp)
      do p = 0, 1
         do q = 0, 1
            component = 1 + p + 2*q
            term = components(:, :, component)
            if (q == 1) term = matmul(hamiltonian, term)
            if (p == 1) term = matmul(term, hamiltonian)
            operator = operator + term
         end do
      end do
      deallocate(components, term)
   end subroutine dresp09_pauli_direct_source

   subroutine validate_inputs(space, radial_bases, source_field, circular_channel, components)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: source_field(:), components(:, :, :)
      integer, intent(in) :: circular_channel
      integer :: isite, n, lmax
      real(rp), parameter :: mesh_tolerance = 2.0e-13_rp

      if (space%nchannel /= 1 .or. size(source_field) /= space%ndim .or. size(radial_bases) /= space%nsite .or. &
          size(components, 3) /= lmto_product_nbranch) then
         error stop 'DRESP-09R Pauli direct source: representation shape mismatch'
      end if
      if (circular_channel /= lmto_product_channel_plus .and. circular_channel /= lmto_product_channel_minus) then
         error stop 'DRESP-09R Pauli direct source: invalid circular channel'
      end if
      lmax = radial_bases(1)%lmax
      n = 2*(lmax + 1)**2*space%nsite
      if (any(shape(components) /= [n, n, lmto_product_nbranch])) then
         error stop 'DRESP-09R Pauli direct source: component shape mismatch'
      end if
      if (space%response_lmax < 2*lmax .or. space%npoint < 3) then
         error stop 'DRESP-09R Pauli direct source: incomplete response space'
      end if
      do isite = 1, size(radial_bases)
         if (radial_bases(isite)%lmax /= lmax .or. radial_bases(isite)%nspin /= 2 .or. &
             radial_bases(isite)%npoint /= space%npoint .or. .not. allocated(radial_bases(isite)%rofi) .or. &
             .not. allocated(radial_bases(isite)%phi_large) .or. .not. allocated(radial_bases(isite)%phidot_large) .or. &
             .not. allocated(radial_bases(isite)%enu_work)) then
            error stop 'DRESP-09R Pauli direct source: incomplete accepted large-component basis'
         end if
         if (abs(radial_bases(isite)%mesh_a - space%a) > mesh_tolerance*max(1.0_rp, abs(space%a)) .or. &
             abs(radial_bases(isite)%mesh_b - space%b) > mesh_tolerance*max(1.0_rp, abs(space%b)) .or. &
             maxval(abs(radial_bases(isite)%rofi - space%radius)) > mesh_tolerance* &
                max(1.0_rp, maxval(abs(space%radius)))) then
            error stop 'DRESP-09R Pauli direct source: radial mesh provenance mismatch'
         end if
      end do
   end subroutine validate_inputs

end module lr_dresp09_pauli_direct_mod
