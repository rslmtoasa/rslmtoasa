!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief DRESP-06A response-space bridge and independent local-action oracle.
!>
!> DRESP-01 stores the site functional F as the compact vector whose Hermitian
!> contraction gives a site transition amplitude.  Consequently a compact
!> response is observed in the certified site basis as F^H chi F.  This module
!> owns only that representation bridge; it does not construct a susceptibility
!> or an interaction.
!------------------------------------------------------------------------------
module lr_dresp06a_bridge_mod

   use precision_mod, only: rp
   use response_basis_mapping_mod, only: response_super_index, response_unflatten_superindex, &
      response_flatten_superindex
   use lr_response_space_mod, only: response_space_layout
   use lr_projected_site_spin_mod, only: projected_site_spin_contract
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis
   use lr_compact_static_interaction_mod, only: compact_project_point_vector, compact_reconstruct_point_vector
   implicit none
   private

   real(rp), parameter, public :: dresp06a_projection_tolerance = 5.0e-10_rp

   public :: project_compact_response_to_sites
   public :: project_compact_loss_to_sites
   public :: compact_loss_matrix
   public :: site_loss_matrix
   public :: compact_apply_local_operator_independent
   public :: dresp06a_relative_matrix_difference

   interface project_compact_response_to_sites
      module procedure project_compact_response_matrix_to_sites
      module procedure project_compact_response_stack_to_sites
   end interface project_compact_response_to_sites

   interface project_compact_loss_to_sites
      module procedure project_compact_loss_matrix_to_sites
      module procedure project_compact_loss_stack_to_sites
   end interface project_compact_loss_to_sites

contains

   !> Project one compact response matrix with the exact DRESP-01 functional.
   subroutine project_compact_response_matrix_to_sites(contract, product, compact_response, site_response)
      type(projected_site_spin_contract), intent(in) :: contract
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: compact_response(:, :)
      complex(rp), intent(out) :: site_response(:, :)
      complex(rp), allocatable :: functionals(:, :), image(:)
      integer :: site, frequency_dummy

      if (any(shape(compact_response) /= [product%product_dimension, product%product_dimension])) then
         error stop 'DRESP-06A compact-to-site bridge: compact response shape mismatch'
      end if
      if (any(shape(site_response) /= [contract%nsite, contract%nsite])) then
         error stop 'DRESP-06A compact-to-site bridge: site response shape mismatch'
      end if
      allocate(functionals(product%product_dimension, contract%nsite), image(product%product_dimension))
      call contract%site_integration_functional(product, functionals)
      site_response = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, contract%nsite
         image = matmul(compact_response, functionals(:, site))
         do frequency_dummy = 1, contract%nsite
            site_response(frequency_dummy, site) = dot_product(functionals(:, frequency_dummy), image)
         end do
      end do
      deallocate(functionals, image)
   end subroutine project_compact_response_matrix_to_sites

   !> Project a compact response stack with the same F^H chi F operation.
   subroutine project_compact_response_stack_to_sites(contract, product, compact_response, site_response)
      type(projected_site_spin_contract), intent(in) :: contract
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: compact_response(:, :, :)
      complex(rp), intent(out) :: site_response(:, :, :)
      integer :: ifrequency

      if (size(compact_response, 1) /= product%product_dimension .or. &
          size(compact_response, 2) /= product%product_dimension) then
         error stop 'DRESP-06A compact-to-site bridge: compact response stack shape mismatch'
      end if
      if (any(shape(site_response) /= [contract%nsite, contract%nsite, size(compact_response, 3)])) then
         error stop 'DRESP-06A compact-to-site bridge: site response stack shape mismatch'
      end if
      do ifrequency = 1, size(compact_response, 3)
         call project_compact_response_matrix_to_sites(contract, product, compact_response(:, :, ifrequency), &
            site_response(:, :, ifrequency))
      end do
   end subroutine project_compact_response_stack_to_sites

   !> Construct compact loss and project it directly to the site basis.
   subroutine project_compact_loss_matrix_to_sites(contract, product, compact_response, site_loss)
      type(projected_site_spin_contract), intent(in) :: contract
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: compact_response(:, :)
      complex(rp), intent(out) :: site_loss(:, :)
      complex(rp), allocatable :: loss(:,:)

      allocate(loss(product%product_dimension, product%product_dimension))
      call compact_loss_matrix(compact_response, loss)
      call project_compact_response_matrix_to_sites(contract, product, loss, site_loss)
      deallocate(loss)
   end subroutine project_compact_loss_matrix_to_sites

   subroutine project_compact_loss_stack_to_sites(contract, product, compact_response, site_loss)
      type(projected_site_spin_contract), intent(in) :: contract
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: compact_response(:, :, :)
      complex(rp), intent(out) :: site_loss(:, :, :)
      integer :: ifrequency

      if (any(shape(site_loss) /= [contract%nsite, contract%nsite, size(compact_response, 3)])) then
         error stop 'DRESP-06A compact-to-site bridge: site loss stack shape mismatch'
      end if
      do ifrequency = 1, size(compact_response, 3)
         call project_compact_loss_matrix_to_sites(contract, product, compact_response(:, :, ifrequency), &
            site_loss(:, :, ifrequency))
      end do
   end subroutine project_compact_loss_stack_to_sites

   !> Ordinary Hermitian loss in the weighted-orthonormal compact basis.
   subroutine compact_loss_matrix(response, loss)
      complex(rp), intent(in) :: response(:, :)
      complex(rp), intent(out) :: loss(:, :)

      if (any(shape(loss) /= shape(response))) error stop 'DRESP-06A compact loss: shape mismatch'
      loss = -(response - conjg(transpose(response)))/cmplx(0.0_rp, 2.0_rp*acos(-1.0_rp), rp)
   end subroutine compact_loss_matrix

   subroutine site_loss_matrix(response, loss)
      complex(rp), intent(in) :: response(:, :)
      complex(rp), intent(out) :: loss(:, :)

      if (any(shape(loss) /= shape(response))) error stop 'DRESP-06A site loss: shape mismatch'
      loss = -(response - conjg(transpose(response)))/cmplx(0.0_rp, 2.0_rp*acos(-1.0_rp), rp)
   end subroutine site_loss_matrix

   !> Independent point/radial action oracle for a compact local operator.
   !>
   !> The production projection is intentionally not called here.  The vector
   !> is reconstructed, multiplied pointwise by the supplied local scalar, and
   !> projected back using only the LR-04 vector transforms.
   subroutine compact_apply_local_operator_independent(space, product, values, vector, result)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      real(rp), intent(in) :: values(:, :)
      complex(rp), intent(in) :: vector(:)
      complex(rp), intent(out) :: result(:)
      complex(rp), allocatable :: point(:), field(:)
      integer :: flat, site, response_l, response_m, radial_point, channel
      type(response_super_index) :: item

      if (any(shape(values) /= [space%nsite, space%npoint])) then
         error stop 'DRESP-06A independent local action: scalar shape mismatch'
      end if
      if (size(vector) /= product%product_dimension .or. size(result) /= product%product_dimension) then
         error stop 'DRESP-06A independent local action: compact vector shape mismatch'
      end if
      allocate(point(space%ndim), field(space%ndim))
      call compact_reconstruct_point_vector(space, product, vector, point)
      field = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, &
            space%nchannel, item)
         site = item%site
         radial_point = item%radial_point
         response_l = item%response_l
         response_m = item%response_m
         channel = item%channel
         if (response_l < 0 .or. response_m < -response_l .or. channel < 1) then
            error stop 'DRESP-06A independent local action: invalid response index'
         end if
         field(flat) = cmplx(values(site, radial_point), 0.0_rp, rp)*point(flat)
      end do
      call compact_project_point_vector(space, product, field, result)
      deallocate(point, field)
   end subroutine compact_apply_local_operator_independent

   real(rp) function dresp06a_relative_matrix_difference(left, right) result(value)
      complex(rp), intent(in) :: left(:, :), right(:, :)
      real(rp) :: left_norm, right_norm

      if (any(shape(left) /= shape(right))) error stop 'DRESP-06A matrix difference: shape mismatch'
      left_norm = sqrt(sum(abs(left)**2))
      right_norm = sqrt(sum(abs(right)**2))
      value = sqrt(sum(abs(left - right)**2))/max(left_norm, right_norm, tiny(1.0_rp))
   end function dresp06a_relative_matrix_difference

end module lr_dresp06a_bridge_mod
