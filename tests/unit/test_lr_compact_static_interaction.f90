!------------------------------------------------------------------------------
! TDVK-06 compact point/product representation and local-interaction oracle.
!
! This test deliberately evaluates the projection/reconstruction identities by
! explicit nested loops over the stored weighted SVD modes.  It does not call
! the adapter under test to construct the oracle.
!------------------------------------------------------------------------------
program test_lr_compact_static_interaction
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
   use lr_compact_static_interaction_mod, only: compact_project_point_vector, compact_reconstruct_point_vector, &
      compact_project_local_operator, compact_apply_local_operator, compact_project_magnetization, &
      lr_compact_mapping_contract
   implicit none

   integer, parameter :: nr = 51, orbital_lmax = 2, response_lmax = 2, nsite = 1, nchannel = 1
   real(rp), parameter :: mesh_a = 0.03_rp, mesh_b = 0.10_rp, nuclear_z = 1.0_rp
   real(rp), parameter :: tolerance = 2.0e-11_rp
   real(rp) :: radius(nr), local_values(nsite, nr), magnetization(nsite, nr)
   type(lmto_radial_basis) :: radial
   type(response_space_layout) :: space
   type(lmto_product_response_basis) :: product
   complex(rp), allocatable :: point(:), compact(:), reconstructed(:), point_oracle(:)
   complex(rp), allocatable :: compact_operator(:, :), operator_oracle(:, :), applied(:), applied_oracle(:)
   complex(rp), allocatable :: compact_magnetization(:), point_magnetization(:)
   real(rp) :: error, max_error
   integer :: i, j, site, response_l, response_m, mode, ir, flat, point_flat
   integer :: site_i, site_j, response_l_i, response_l_j, response_m_i, response_m_j, mode_i, mode_j
   type(response_super_index) :: item
   logical :: failed

   call basis_init(orbital_lmax)
   call init_math_operators()
   call g_logger%init()
   call build_mesh(radius)
   call setup_radial_basis(radial, radius)
   call space%initialize(nsite, response_lmax, radius, mesh_a, mesh_b, nchannel)
   call product%initialize(space, [radial], lmto_product_channel_plus, .false.)

   failed = .false.
   max_error = 0.0_rp
   if (product%product_dimension < 1) failed = .true.
   allocate(point(space%ndim), compact(product%product_dimension), reconstructed(space%ndim), point_oracle(space%ndim))
   allocate(compact_operator(product%product_dimension, product%product_dimension), &
      operator_oracle(product%product_dimension, product%product_dimension), applied(product%product_dimension), &
      applied_oracle(product%product_dimension))
   allocate(compact_magnetization(product%product_dimension), point_magnetization(space%ndim))

   do i = 1, space%ndim
      call response_unflatten_superindex(i, nsite, response_lmax, nr, nchannel, item)
      point(i) = cmplx(0.17_rp + 0.013_rp*real(item%radial_point, rp) + 0.007_rp*real(item%response_l, rp), &
         -0.021_rp*real(item%response_m, rp) + 0.003_rp*real(item%radial_point, rp), rp)
      magnetization(1, item%radial_point) = 0.25_rp + 0.01_rp*real(item%radial_point, rp)
   end do
   do site = 1, nsite
      do ir = 1, nr
         local_values(site, ir) = 0.09_rp + 0.017_rp*real(ir, rp) + 0.03_rp*real(site, rp)
      end do
   end do

   call compact_project_point_vector(space, product, point, compact)
   compact_magnetization = cmplx(0.0_rp, 0.0_rp, rp)
   do flat = 1, product%product_dimension
      call product%unflatten_index(flat, site, response_l, response_m, mode)
      do ir = 2, nr
         item = response_super_index(site, response_l, response_m, ir, 1)
         call response_flatten_superindex(item, nsite, response_lmax, nr, nchannel, point_flat)
         compact_magnetization(flat) = compact_magnetization(flat) + &
            conjg(product%blocks(site, response_l)%weighted_modes(ir, mode))*sqrt(space%radial_weights(ir))*point(point_flat)
      end do
   end do
   call check_vector(compact, compact_magnetization, 'point-to-compact projection', failed, max_error)

   call compact_reconstruct_point_vector(space, product, compact, reconstructed)
   point_oracle = cmplx(0.0_rp, 0.0_rp, rp)
   do flat = 1, product%product_dimension
      call product%unflatten_index(flat, site, response_l, response_m, mode)
      do ir = 2, nr
         item = response_super_index(site, response_l, response_m, ir, 1)
         call response_flatten_superindex(item, nsite, response_lmax, nr, nchannel, point_flat)
         point_oracle(point_flat) = point_oracle(point_flat) + &
            product%blocks(site, response_l)%weighted_modes(ir, mode)/sqrt(space%radial_weights(ir))*compact(flat)
      end do
   end do
   call check_vector(reconstructed, point_oracle, 'compact-to-point reconstruction', failed, max_error)
   if (abs(reconstructed(1)) > tolerance) failed = .true.

   call compact_project_local_operator(space, product, cmplx(local_values, 0.0_rp, rp), compact_operator)
   operator_oracle = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, product%product_dimension
      call product%unflatten_index(i, site_i, response_l_i, response_m_i, mode_i)
      do j = 1, product%product_dimension
         call product%unflatten_index(j, site_j, response_l_j, response_m_j, mode_j)
         if (site_i /= site_j .or. response_l_i /= response_l_j .or. response_m_i /= response_m_j) cycle
         do ir = 2, nr
            operator_oracle(i, j) = operator_oracle(i, j) + &
               conjg(product%blocks(site_i, response_l_i)%weighted_modes(ir, mode_i))*local_values(site_i, ir)* &
               product%blocks(site_i, response_l_i)%weighted_modes(ir, mode_j)
         end do
      end do
   end do
   call check_matrix(compact_operator, operator_oracle, 'local operator U^H K U', failed, max_error)
   call compact_apply_local_operator(space, product, cmplx(local_values, 0.0_rp, rp), compact, applied)
   applied_oracle = matmul(operator_oracle, compact)
   call check_vector(applied, applied_oracle, 'local operator action', failed, max_error)

   call compact_project_magnetization(space, product, magnetization, compact_magnetization, point_magnetization)
   do site = 1, nsite
      do ir = 1, nr
         item = response_super_index(site, 0, 0, ir, 1)
         call response_flatten_superindex(item, nsite, response_lmax, nr, nchannel, point_flat)
         if (abs(point_magnetization(point_flat) - cmplx(sqrt(4.0_rp*response_angular_pi)* &
            magnetization(site, ir), 0.0_rp, rp)) > tolerance) failed = .true.
      end do
   end do

   write (*, '(a)') 'mapping_contract='//trim(lr_compact_mapping_contract)
   write (*, '(a,i0)') 'compact_dimension=', product%product_dimension
   write (*, '(a,es12.4)') 'maximum_oracle_error=', max_error
   if (failed) then
      write (*, '(a)') 'UnitLrCompactStaticInteraction: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrCompactStaticInteraction: PASS (projection, reconstruction, local operator, action, m_00 oracle)'

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

   subroutine check_vector(actual, expected, label, test_failed, maximum)
      complex(rp), intent(in) :: actual(:), expected(:)
      character(len=*), intent(in) :: label
      logical, intent(inout) :: test_failed
      real(rp), intent(inout) :: maximum

      error = maxval(abs(actual - expected))
      maximum = max(maximum, error)
      write (*, '(a,es12.4)') trim(label)//' error=', error
      if (error > tolerance) test_failed = .true.
   end subroutine check_vector

   subroutine check_matrix(actual, expected, label, test_failed, maximum)
      complex(rp), intent(in) :: actual(:, :), expected(:, :)
      character(len=*), intent(in) :: label
      logical, intent(inout) :: test_failed
      real(rp), intent(inout) :: maximum

      error = maxval(abs(actual - expected))
      maximum = max(maximum, error)
      write (*, '(a,es12.4)') trim(label)//' error=', error
      if (error > tolerance) test_failed = .true.
   end subroutine check_matrix

end program test_lr_compact_static_interaction
