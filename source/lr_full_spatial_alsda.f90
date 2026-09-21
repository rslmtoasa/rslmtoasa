!------------------------------------------------------------------------------
! DRESP-10 full-spatial scalar-relativistic transverse source insertion.
!
! The source is a pointwise Pauli field in the complete (site,L,M,r) response
! space.  The six endpoint branches and the four requested scalar-relativistic
! pieces are retained explicitly; the contracted operator is formed only after
! that decomposition.  In particular, this module does not call the historical
! L=0 source or the compact U^H K U projection.
!------------------------------------------------------------------------------
module lr_full_spatial_alsda_mod

   use precision_mod, only: rp
   use response_basis_mapping_mod, only: response_super_index, response_unflatten_superindex
   use response_angular_basis_mod, only: response_angular_pi
   use lr_response_space_mod, only: response_space_layout
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_lmto_product_response_basis_mod, only: lmto_product_nbranch, lmto_product_branch_powers
   use lr_sr_spatial_augmentation_tangent_mod, only: sr_spatial_branch_observable
   implicit none
   private

   integer, parameter, public :: lr_full_spatial_piece_upper = 1
   integer, parameter, public :: lr_full_spatial_piece_lower_small = 2
   integer, parameter, public :: lr_full_spatial_piece_lower_rank0 = 3
   integer, parameter, public :: lr_full_spatial_piece_lower_rank2 = 4
   integer, parameter, public :: lr_full_spatial_piece_total = 5
   integer, parameter, public :: lr_full_spatial_npiece = 5

   public :: lr_full_spatial_source_components
   public :: lr_full_spatial_contract_components
   public :: lr_full_spatial_source
   public :: lr_full_spatial_source_pieces

contains

   !> Build the six endpoint branches, split as
   !> upper/small/rank-0/rank-2/total.  The returned matrices are physical
   !> source matrices and include the response radial quadrature weight.
   subroutine lr_full_spatial_source_components(space, radial_bases, source_field, circular_channel, components)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: source_field(:)
      integer, intent(in) :: circular_channel
      complex(rp), intent(out) :: components(:, :, :, :)
      type(response_super_index) :: item
      complex(rp) :: upper_x(2,2), small_x(2,2), rank0_x(2,2), rank2_x(2,2), total_x(2,2)
      complex(rp) :: upper_y(2,2), small_y(2,2), rank0_y(2,2), rank2_y(2,2), total_y(2,2)
      complex(rp) :: source_matrix(2,2)
      integer :: norb, n, flat, site, ir, iorb, jorb, l, lp, m, mp, branch, piece
      integer :: row_up, row_down, col_up, col_down, row, col
      real(rp) :: radial_weight

      call validate_inputs(space, radial_bases, source_field, components, circular_channel)
      norb = (radial_bases(1)%lmax + 1)**2
      n = 2*norb*space%nsite
      components = cmplx(0.0_rp, 0.0_rp, rp)

      do flat = 1, size(source_field)
         if (source_field(flat) == cmplx(0.0_rp, 0.0_rp, rp)) cycle
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, &
            space%nchannel, item)
         if (item%channel /= 1) cycle
         ir = item%radial_point
         if (space%radial_weights(ir) <= 0.0_rp) cycle
         site = item%site
         radial_weight = space%radial_weights(ir)

         do iorb = 1, norb
            l = orbital_l_from_index(iorb)
            m = iorb - l*l - l - 1
            row_up = (site - 1)*2*norb + iorb
            row_down = row_up + norb
            do jorb = 1, norb
               lp = orbital_l_from_index(jorb)
               mp = jorb - lp*lp - lp - 1
               col_up = (site - 1)*2*norb + jorb
               col_down = col_up + norb

               do branch = 1, lmto_product_nbranch
                  call sr_spatial_branch_observable(radial_bases(site), ir, l, m, lp, mp, branch, &
                     item%response_l, item%response_m, 1, upper_x, small_x, rank0_x, rank2_x, total_x)
                  call sr_spatial_branch_observable(radial_bases(site), ir, l, m, lp, mp, branch, &
                     item%response_l, item%response_m, 2, upper_y, small_y, rank0_y, rank2_y, total_y)
                  call circular_source(upper_x, upper_y, circular_channel, source_matrix)
                  call add_spin_matrix(components, row_up, row_down, col_up, col_down, branch, &
                     lr_full_spatial_piece_upper, source_field(flat)*radial_weight, source_matrix)
                  call circular_source(small_x, small_y, circular_channel, source_matrix)
                  call add_spin_matrix(components, row_up, row_down, col_up, col_down, branch, &
                     lr_full_spatial_piece_lower_small, source_field(flat)*radial_weight, source_matrix)
                  call circular_source(rank0_x, rank0_y, circular_channel, source_matrix)
                  call add_spin_matrix(components, row_up, row_down, col_up, col_down, branch, &
                     lr_full_spatial_piece_lower_rank0, source_field(flat)*radial_weight, source_matrix)
                  call circular_source(rank2_x, rank2_y, circular_channel, source_matrix)
                  call add_spin_matrix(components, row_up, row_down, col_up, col_down, branch, &
                     lr_full_spatial_piece_lower_rank2, source_field(flat)*radial_weight, source_matrix)
                  call circular_source(total_x, total_y, circular_channel, source_matrix)
                  call add_spin_matrix(components, row_up, row_down, col_up, col_down, branch, &
                     lr_full_spatial_piece_total, source_field(flat)*radial_weight, source_matrix)
               end do
            end do
         end do
      end do
   end subroutine lr_full_spatial_source_components

   !> Contract all six branches and return the requested piece (total by
   !> default).  Branch ordering is 00/10/01/11/20/02.
   subroutine lr_full_spatial_source(space, radial_bases, source_field, circular_channel, hamiltonian, operator, &
                                      piece)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: source_field(:), hamiltonian(:, :)
      integer, intent(in) :: circular_channel
      complex(rp), intent(out) :: operator(:, :)
      integer, intent(in), optional :: piece
      complex(rp), allocatable :: components(:, :, :, :)
      integer :: selected_piece, branch, p, q

      selected_piece = lr_full_spatial_piece_total
      if (present(piece)) selected_piece = piece
      if (selected_piece < 1 .or. selected_piece > lr_full_spatial_npiece) then
         error stop 'lr_full_spatial_source: invalid scalar-relativistic piece'
      end if
      if (size(hamiltonian, 1) /= size(hamiltonian, 2) .or. any(shape(operator) /= shape(hamiltonian))) then
         error stop 'lr_full_spatial_source: Hamiltonian/operator shape mismatch'
      end if
      allocate(components(size(hamiltonian,1), size(hamiltonian,1), lmto_product_nbranch, lr_full_spatial_npiece))
      call lr_full_spatial_source_components(space, radial_bases, source_field, circular_channel, components)
      call contract_components(components, hamiltonian, selected_piece, operator)
      deallocate(components)
   end subroutine lr_full_spatial_source

   !> Contract a previously built, k-independent six-branch source map.
   !> This is the production path for a fixed point-space ALSDA field: all
   !> radial/angular work is performed once, while only the H^(p/q) branch
   !> algebra remains inside the k loop.
   subroutine lr_full_spatial_contract_components(components, hamiltonian, operators)
      complex(rp), intent(in) :: components(:, :, :, :), hamiltonian(:, :)
      complex(rp), intent(out) :: operators(:, :, :)
      integer :: piece

      if (size(hamiltonian,1) /= size(hamiltonian,2) .or. &
          any(shape(components) /= [size(hamiltonian,1), size(hamiltonian,1), lmto_product_nbranch, lr_full_spatial_npiece]) .or. &
          any(shape(operators) /= [size(hamiltonian,1), size(hamiltonian,1), lr_full_spatial_npiece])) then
         error stop 'lr_full_spatial_contract_components: shape mismatch'
      end if
      do piece = 1, lr_full_spatial_npiece
         call contract_components(components, hamiltonian, piece, operators(:,:,piece))
      end do
   end subroutine lr_full_spatial_contract_components

   !> Return one matrix per requested scalar-relativistic piece after branch
   !> contraction.  This is used for the DRESP-10 fixed/lower decomposition.
   subroutine lr_full_spatial_source_pieces(space, radial_bases, source_field, circular_channel, hamiltonian, operators)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: source_field(:), hamiltonian(:, :)
      integer, intent(in) :: circular_channel
      complex(rp), intent(out) :: operators(:, :, :)
      complex(rp), allocatable :: components(:, :, :, :)
      integer :: piece

      if (any(shape(operators) /= [size(hamiltonian,1), size(hamiltonian,2), lr_full_spatial_npiece])) then
         error stop 'lr_full_spatial_source_pieces: output shape mismatch'
      end if
      allocate(components(size(hamiltonian,1), size(hamiltonian,1), lmto_product_nbranch, lr_full_spatial_npiece))
      call lr_full_spatial_source_components(space, radial_bases, source_field, circular_channel, components)
      do piece = 1, lr_full_spatial_npiece
         call contract_components(components, hamiltonian, piece, operators(:,:,piece))
      end do
      deallocate(components)
   end subroutine lr_full_spatial_source_pieces

   subroutine contract_components(components, hamiltonian, piece, operator)
      complex(rp), intent(in) :: components(:, :, :, :), hamiltonian(:, :)
      integer, intent(in) :: piece
      complex(rp), intent(out) :: operator(:, :)
      complex(rp), allocatable :: term(:, :)
      integer :: branch, p, q

      allocate(term(size(hamiltonian,1), size(hamiltonian,1)))
      operator = cmplx(0.0_rp, 0.0_rp, rp)
      do branch = 1, lmto_product_nbranch
         call lmto_product_branch_powers(branch, p, q)
         term = components(:,:,branch,piece)
         if (q == 1) term = matmul(hamiltonian, term)
         if (p == 1) term = matmul(term, hamiltonian)
         operator = operator + term
      end do
      deallocate(term)
   end subroutine contract_components

   subroutine add_spin_matrix(components, row_up, row_down, col_up, col_down, branch, piece, coefficient, matrix)
      complex(rp), intent(inout) :: components(:, :, :, :)
      integer, intent(in) :: row_up, row_down, col_up, col_down, branch, piece
      complex(rp), intent(in) :: coefficient, matrix(2,2)

      components(row_up, col_up, branch, piece) = components(row_up, col_up, branch, piece) + coefficient*matrix(1,1)
      components(row_up, col_down, branch, piece) = components(row_up, col_down, branch, piece) + coefficient*matrix(1,2)
      components(row_down, col_up, branch, piece) = components(row_down, col_up, branch, piece) + coefficient*matrix(2,1)
      components(row_down, col_down, branch, piece) = components(row_down, col_down, branch, piece) + coefficient*matrix(2,2)
   end subroutine add_spin_matrix

   pure subroutine circular_source(x_matrix, y_matrix, circular_channel, source_matrix)
      complex(rp), intent(in) :: x_matrix(2,2), y_matrix(2,2)
      integer, intent(in) :: circular_channel
      complex(rp), intent(out) :: source_matrix(2,2)
      if (circular_channel == 1) then
         source_matrix = 0.5_rp*(x_matrix - cmplx(0.0_rp, 1.0_rp, rp)*y_matrix)
      else
         source_matrix = 0.5_rp*(x_matrix + cmplx(0.0_rp, 1.0_rp, rp)*y_matrix)
      end if
   end subroutine circular_source

   pure integer function orbital_l_from_index(iorb) result(l)
      integer, intent(in) :: iorb
      integer :: trial
      l = -1
      do trial = 0, 8
         if (trial*trial + 1 <= iorb .and. iorb <= (trial + 1)*(trial + 1)) then
            l = trial
            return
         end if
      end do
   end function orbital_l_from_index

   subroutine validate_inputs(space, radial_bases, source_field, components, circular_channel)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: source_field(:), components(:, :, :, :)
      integer, intent(in) :: circular_channel
      integer :: n

      if (space%nchannel /= 1 .or. space%response_lmax < 0) then
         error stop 'lr_full_spatial_source: a single Pauli source channel is required'
      end if
      if (size(radial_bases) /= space%nsite .or. size(source_field) /= space%ndim) then
         error stop 'lr_full_spatial_source: response/radial/source shape mismatch'
      end if
      if (circular_channel /= 1 .and. circular_channel /= 2) then
         error stop 'lr_full_spatial_source: invalid circular channel'
      end if
      n = 2*(radial_bases(1)%lmax + 1)**2*space%nsite
      if (any(shape(components) /= [n, n, lmto_product_nbranch, lr_full_spatial_npiece])) then
         error stop 'lr_full_spatial_source: component shape mismatch'
      end if
      if (space%response_lmax > 4) then
         error stop 'lr_full_spatial_source: DRESP-10 source map is certified only through L=4'
      end if
   end subroutine validate_inputs

end module lr_full_spatial_alsda_mod
