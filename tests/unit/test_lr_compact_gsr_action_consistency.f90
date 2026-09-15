!------------------------------------------------------------------------------
! TDVK-06G compact GSR action-consistency regression.
!
! The point-space side of each action is assembled here by independent loops;
! it is deliberately not obtained from compact_project_local_operator or
! compact_reconstruct_point_vector.  The live GSR evaluator is also exercised
! with an identity susceptibility so its constructed Gamma columns are tested
! directly against the independent point-space actions.
!------------------------------------------------------------------------------
program test_lr_compact_gsr_action_consistency
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use logger_mod, only: g_logger
   use math_mod, only: init_math_operators
   use self_mod, only: legacy_radial_fixture
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex, &
      response_unflatten_superindex
   use response_angular_basis_mod, only: response_angular_pi
   use lr_response_space_mod, only: response_space_layout
   use lr_lmto_product_response_basis_mod, only: lmto_product_channel_plus, lmto_product_response_basis
   use lr_compact_static_interaction_mod, only: lr_compact_gsr_result, compact_project_magnetization, &
      compact_project_local_operator, compact_reconstruct_point_vector, compact_weighted_projection_diagnostics, &
      evaluate_compact_goldstone_sumrule
   implicit none

   integer, parameter :: nr = 51, orbital_lmax = 2, response_lmax = 2, nsite = 1, nchannel = 1
   real(rp), parameter :: mesh_a = 0.03_rp, mesh_b = 0.10_rp, nuclear_z = 1.0_rp
   real(rp), parameter :: tolerance = 2.0e-10_rp
   real(rp) :: radius(nr), magnetization(nsite, nr)
   type(lmto_radial_basis) :: radial
   type(response_space_layout) :: space
   type(lmto_product_response_basis) :: product
   type(lr_compact_gsr_result) :: gsr
   complex(rp), allocatable :: physical_compact(:), physical_point(:), physical_projected_point(:)
   complex(rp), allocatable :: mixed_compact(:), mixed_point(:)
   complex(rp), allocatable :: values(:, :), basis_values(:, :, :), basis_operator(:, :)
   complex(rp), allocatable :: compact_action(:, :), point_action(:, :)
   complex(rp), allocatable :: mixed_compact_action(:, :), mixed_point_action(:, :)
   complex(rp), allocatable :: coefficients(:), assembled_sum(:), assembled_operator(:, :), assembled_action(:)
   complex(rp), allocatable :: identity(:, :)
   real(rp) :: weighted_norm, projection_norm, projection_relative
   real(rp) :: single_action_error, mixed_action_error, gamma_error, assembled_error
   real(rp) :: assembled_sum_norm, assembled_operator_norm, assembled_relative
   real(rp) :: residual_difference_relative
   integer :: i, ir, unknown, nunknown, site, response_l, mode, flat, point_flat
   logical :: failed
   type(response_super_index) :: item

   call basis_init(orbital_lmax)
   call init_math_operators()
   call g_logger%init()
   call build_mesh(radius)
   call setup_radial_basis(radial, radius)
   call space%initialize(nsite, response_lmax, radius, mesh_a, mesh_b, nchannel)
   call product%initialize(space, [radial], lmto_product_channel_plus, .false.)

   nunknown = 0
   do site = 1, product%nsite
      nunknown = nunknown + product%blocks(site, 0)%rank
   end do
   if (nunknown < 3) error stop 'TDVK-06G fixture does not contain several retained L=0 GSR modes'

   allocate(physical_compact(product%product_dimension), physical_point(space%ndim), &
      physical_projected_point(space%ndim), mixed_compact(product%product_dimension), mixed_point(space%ndim), &
      values(nsite, nr), basis_values(nsite, nr, nunknown), basis_operator(product%product_dimension, product%product_dimension), &
      compact_action(product%product_dimension, nunknown), point_action(product%product_dimension, nunknown), &
      mixed_compact_action(product%product_dimension, nunknown), mixed_point_action(product%product_dimension, nunknown), &
      coefficients(nunknown), assembled_sum(product%product_dimension), &
      assembled_operator(product%product_dimension, product%product_dimension), assembled_action(product%product_dimension), &
      identity(product%product_dimension, product%product_dimension))

   do ir = 1, nr
      magnetization(1, ir) = 0.19_rp + 0.013_rp*real(ir, rp) + 0.00017_rp*real(ir*ir, rp)
   end do
   call compact_project_magnetization(space, product, magnetization, physical_compact, physical_point)
   call independent_reconstruct(space, product, physical_compact, physical_projected_point)
   call compact_weighted_projection_diagnostics(space, physical_point, physical_projected_point, weighted_norm, &
      projection_norm, projection_relative)
   failed = projection_norm <= 1.0e-12_rp

   do i = 1, product%product_dimension
      mixed_compact(i) = cmplx(0.071_rp + 0.0023_rp*real(i, rp), &
         -0.043_rp + 0.0017_rp*real(mod(i, 11), rp), rp)
   end do
   call independent_reconstruct(space, product, mixed_compact, mixed_point)

   basis_values = cmplx(0.0_rp, 0.0_rp, rp)
   unknown = 0
   do site = 1, product%nsite
      do mode = 1, product%blocks(site, 0)%rank
         unknown = unknown + 1
         do ir = 2, nr
            basis_values(site, ir, unknown) = 4.0_rp*response_angular_pi* &
               product%blocks(site, 0)%weighted_modes(ir, mode)/sqrt(space%radial_weights(ir))
         end do
      end do
   end do

   do unknown = 1, nunknown
      values = basis_values(:, :, unknown)
      call compact_project_local_operator(space, product, values, basis_operator)
      compact_action(:, unknown) = matmul(basis_operator, physical_compact)
      mixed_compact_action(:, unknown) = matmul(basis_operator, mixed_compact)
      call independent_operator_action(space, product, values, physical_projected_point, point_action(:, unknown))
      call independent_operator_action(space, product, values, mixed_point, mixed_point_action(:, unknown))
   end do
   single_action_error = maxval(abs(compact_action - point_action))
   mixed_action_error = maxval(abs(mixed_compact_action - mixed_point_action))
   if (single_action_error > tolerance .or. mixed_action_error > tolerance) failed = .true.

   coefficients = cmplx(0.23_rp, -0.17_rp, rp)
   do unknown = 1, nunknown
      coefficients(unknown) = coefficients(unknown) + cmplx(0.011_rp*real(unknown, rp), &
         0.007_rp*real(unknown*unknown, rp), rp)
   end do
   assembled_sum = cmplx(0.0_rp, 0.0_rp, rp)
   assembled_operator = cmplx(0.0_rp, 0.0_rp, rp)
   do unknown = 1, nunknown
      values = basis_values(:, :, unknown)
      call compact_project_local_operator(space, product, values, basis_operator)
      assembled_sum = assembled_sum + coefficients(unknown)*point_action(:, unknown)
      assembled_operator = assembled_operator + coefficients(unknown)*basis_operator
   end do
   assembled_action = matmul(assembled_operator, physical_compact)
   assembled_sum_norm = sqrt(sum(abs(assembled_sum)**2))
   assembled_operator_norm = sqrt(sum(abs(assembled_action)**2))
   assembled_error = maxval(abs(assembled_sum - assembled_action))
   assembled_relative = sqrt(sum(abs(assembled_sum - assembled_action)**2))/ &
      max(assembled_operator_norm, tiny(1.0_rp))
   if (assembled_error > tolerance) failed = .true.

   identity = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, product%product_dimension
      identity(i, i) = cmplx(1.0_rp, 0.0_rp, rp)
   end do
   call evaluate_compact_goldstone_sumrule(space, product, identity, magnetization, gsr)
   gamma_error = maxval(abs(gsr%gamma - point_action))
   residual_difference_relative = gsr%residual_difference_relative
   if (gamma_error > tolerance) failed = .true.
   if (residual_difference_relative > 1.0e-8_rp) failed = .true.

   write (*, '(a,es16.8)') 'magnetization_weighted_norm=', weighted_norm
   write (*, '(a,es16.8)') 'magnetization_projection_residual_norm=', projection_norm
   write (*, '(a,es16.8)') 'magnetization_projection_relative=', projection_relative
   write (*, '(a,es16.8)') 'single_gsr_action_max_error=', single_action_error
   write (*, '(a,es16.8)') 'mixed_compact_action_max_error=', mixed_action_error
   write (*, '(a,es16.8)') 'assembled_sum_norm=', assembled_sum_norm
   write (*, '(a,es16.8)') 'assembled_operator_norm=', assembled_operator_norm
   write (*, '(a,es16.8)') 'assembled_action_max_error=', assembled_error
   write (*, '(a,es16.8)') 'assembled_action_relative_error=', assembled_relative
   write (*, '(a,es16.8)') 'live_gsr_gamma_max_error=', gamma_error
   write (*, '(a,es16.8)') 'live_gsr_solve_residual_norm=', gsr%equation_residual_norm
   write (*, '(a,es16.8)') 'live_gsr_reconstructed_residual_norm=', gsr%residual_norm
   write (*, '(a,es16.8)') 'live_gsr_residual_difference_norm=', gsr%residual_difference_norm
   write (*, '(a,es16.8)') 'live_gsr_assembled_action_difference_relative=', gsr%assembled_action_difference_relative
   write (*, '(a,es16.8)') 'live_gsr_solve_reconstructed_relative=', residual_difference_relative
   if (failed) then
      write (*, '(a)') 'UnitLrCompactGsrActionConsistency: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrCompactGsrActionConsistency: PASS (single-action, assembled-action, live-GSR residual consistency)'

contains

   subroutine build_mesh(mesh)
      real(rp), intent(out) :: mesh(:)
      integer :: ir

      mesh(1) = 0.0_rp
      do ir = 2, size(mesh)
         mesh(ir) = mesh_b*(exp(mesh_a*real(ir - 1, rp)) - 1.0_rp)
      end do
   end subroutine build_mesh

   subroutine setup_radial_basis(basis, mesh)
      type(lmto_radial_basis), intent(out) :: basis
      real(rp), intent(in) :: mesh(:)
      real(rp) :: potential(size(mesh)), energy
      real(rp), allocatable :: g(:, :), gp(:, :), gpp(:, :)
      real(rp), allocatable :: gpack(:), gdotpack(:), gddotpack(:)
      integer :: ispin, l

      potential = 0.0_rp
      call basis%initialize(size(mesh), orbital_lmax, 2)
      do ispin = 1, 2
         do l = 0, orbital_lmax
            call legacy_radial_fixture(nuclear_z, l, mesh_a, mesh_b, mesh, potential, energy, g, gp, gpp, 2.0_rp)
            gpack = reshape(g, [2*size(mesh)])
            gdotpack = reshape(gp, [2*size(mesh)])
            gddotpack = reshape(gpp, [2*size(mesh)])
            call basis%capture_channel(l, ispin, energy, mesh, potential, mesh_a, mesh_b, nuclear_z, &
               gpack, gdotpack, gddotpack, energy)
         end do
      end do
   end subroutine setup_radial_basis

   subroutine independent_reconstruct(layout, product_basis, compact_vector, point_vector)
      type(response_space_layout), intent(in) :: layout
      type(lmto_product_response_basis), intent(in) :: product_basis
      complex(rp), intent(in) :: compact_vector(:)
      complex(rp), intent(out) :: point_vector(:)
      integer :: compact_flat, site_i, response_l_i, response_m_i, mode_i, radial_i, point_i
      real(rp) :: sqrt_weight
      type(response_super_index) :: local_item

      point_vector = cmplx(0.0_rp, 0.0_rp, rp)
      do compact_flat = 1, product_basis%product_dimension
         call product_basis%unflatten_index(compact_flat, site_i, response_l_i, response_m_i, mode_i)
         do radial_i = 2, layout%npoint
            sqrt_weight = sqrt(layout%radial_weights(radial_i))
            if (sqrt_weight <= 0.0_rp) cycle
            local_item = response_super_index(site_i, response_l_i, response_m_i, radial_i, 1)
            call response_flatten_superindex(local_item, layout%nsite, layout%response_lmax, layout%npoint, &
               layout%nchannel, point_i)
            point_vector(point_i) = point_vector(point_i) + &
               product_basis%blocks(site_i, response_l_i)%weighted_modes(radial_i, mode_i)/sqrt_weight*compact_vector(compact_flat)
         end do
      end do
   end subroutine independent_reconstruct

   subroutine independent_operator_action(layout, product_basis, local_values, point_vector, compact_result)
      type(response_space_layout), intent(in) :: layout
      type(lmto_product_response_basis), intent(in) :: product_basis
      complex(rp), intent(in) :: local_values(:, :), point_vector(:)
      complex(rp), intent(out) :: compact_result(:)
      integer :: flat_out, site_i, response_l_i, response_m_i, mode_i
      integer :: radial_i, point_i
      real(rp) :: sqrt_weight
      type(response_super_index) :: local_item

      compact_result = cmplx(0.0_rp, 0.0_rp, rp)
      do flat_out = 1, product_basis%product_dimension
         call product_basis%unflatten_index(flat_out, site_i, response_l_i, response_m_i, mode_i)
         do radial_i = 2, layout%npoint
            sqrt_weight = sqrt(layout%radial_weights(radial_i))
            if (sqrt_weight <= 0.0_rp) cycle
            local_item = response_super_index(site_i, response_l_i, response_m_i, radial_i, 1)
            call response_flatten_superindex(local_item, layout%nsite, layout%response_lmax, layout%npoint, &
               layout%nchannel, point_i)
            compact_result(flat_out) = compact_result(flat_out) + &
               conjg(product_basis%blocks(site_i, response_l_i)%weighted_modes(radial_i, mode_i))*sqrt_weight* &
               local_values(site_i, radial_i)*point_vector(point_i)
         end do
      end do
   end subroutine independent_operator_action

end program test_lr_compact_gsr_action_consistency
