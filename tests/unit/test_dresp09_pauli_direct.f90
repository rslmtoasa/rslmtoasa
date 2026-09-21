! DRESP-09R independent Pauli projected direct-operator audit.
!
! The test constructs the accepted large-component radial representation, fills
! every L=0..4/M mode with a mixed complex compact field, reconstructs its
! contravariant raw source field independently, and compares all six endpoint
! branches against the compact adjoint.  It never calls the production
! Pauli transition evaluator.
program test_dresp09_pauli_direct
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use logger_mod, only: g_logger
   use math_mod, only: init_math_operators
   use self_mod, only: legacy_radial_fixture
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex
   use lr_response_space_mod, only: response_space_layout
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, lmto_product_channel_plus, &
      lmto_product_channel_minus
   use lr_dresp09_field_insertion_mod, only: dresp09_product_metric_adjoint_operator
   use lr_dresp09_pauli_direct_mod, only: dresp09_pauli_direct_source, dresp09_pauli_direct_source_components
   implicit none

   integer, parameter :: nr = 51, lmax = 2, response_lmax = 4
   real(rp), parameter :: mesh_a = 0.03_rp, mesh_b = 0.10_rp, nuclear_z = 1.0_rp
   real(rp), parameter :: tolerance = 3.0e-11_rp
   real(rp) :: radius(nr), branch_error(6), operator_error(2), coverage_error
   type(lmto_radial_basis) :: radial
   type(response_space_layout) :: space
   type(lmto_product_response_basis) :: product_plus, product_minus
   complex(rp), allocatable :: compact_plus(:), compact_minus(:), source_plus(:), source_minus(:)
   complex(rp), allocatable :: vertices(:, :, :, :), direct_components(:, :, :), expected(:, :, :)
   complex(rp), allocatable :: direct_operator(:, :), compact_operator(:, :), hamiltonian(:, :)
   logical :: failed

   call basis_init(2)
   call init_math_operators()
   call g_logger%init()
   call build_mesh(radius)
   call setup_radial_basis(radial, radius)
   call space%initialize(1, response_lmax, radius, mesh_a, mesh_b, 1)
   ! The independent oracle is an algebra audit; the production strict-rank
   ! policy remains enforced by the material handoff and is not duplicated by
   ! this synthetic fixture.
   call product_plus%initialize(space, [radial], lmto_product_channel_plus, .false.)
   call product_minus%initialize(space, [radial], lmto_product_channel_minus, .false.)

   allocate(hamiltonian(2*(lmax + 1)**2, 2*(lmax + 1)**2))
   call build_hamiltonian(hamiltonian)
   allocate(compact_plus(product_plus%product_dimension), compact_minus(product_minus%product_dimension), &
      source_plus(space%ndim), source_minus(space%ndim))
   call fill_compact(compact_plus)
   call fill_compact(compact_minus)
   call reconstruct_dual_source(space, product_plus, compact_plus, source_plus)
   call reconstruct_dual_source(space, product_minus, compact_minus, source_minus)

   allocate(direct_components(size(hamiltonian, 1), size(hamiltonian, 2), 6), &
      expected(size(hamiltonian, 1), size(hamiltonian, 2), 6), &
      direct_operator(size(hamiltonian, 1), size(hamiltonian, 2)), &
      compact_operator(size(hamiltonian, 1), size(hamiltonian, 2)))
   branch_error = 0.0_rp
   operator_error = 0.0_rp
   coverage_error = 0.0_rp

   call product_plus%component_vertex_tensor(vertices)
   call dresp09_pauli_direct_source_components(space, [radial], source_plus, lmto_product_channel_plus, direct_components)
   call compare_branches(vertices, compact_plus, direct_components, branch_error, coverage_error)
   call dresp09_pauli_direct_source(space, [radial], source_plus, lmto_product_channel_plus, hamiltonian, direct_operator)
   call dresp09_product_metric_adjoint_operator(vertices, hamiltonian, compact_plus, compact_operator)
   operator_error(1) = relative_error(direct_operator, compact_operator)
   deallocate(vertices)

   call product_minus%component_vertex_tensor(vertices)
   call dresp09_pauli_direct_source_components(space, [radial], source_minus, lmto_product_channel_minus, direct_components)
   call compare_branches(vertices, compact_minus, direct_components, branch_error, coverage_error)
   call dresp09_pauli_direct_source(space, [radial], source_minus, lmto_product_channel_minus, hamiltonian, direct_operator)
   call dresp09_product_metric_adjoint_operator(vertices, hamiltonian, compact_minus, compact_operator)
   operator_error(2) = relative_error(direct_operator, compact_operator)

   failed = maxval(branch_error) > tolerance .or. maxval(operator_error) > tolerance .or. coverage_error > tolerance
   if (failed) then
      write (*, '(a)') 'UnitDresp09PauliDirect: FAIL'
      write (*, '(a,6(es12.4,1x))') '  B00/B10/B01/B11/B20/B02=', branch_error
      write (*, '(a,2(es12.4,1x))') '  plus/minus_operator=', operator_error
      write (*, '(a,es12.4)') '  mixed_LM_coverage=', coverage_error
      error stop 1
   end if
   write (*, '(a,6(es12.4,1x))') '  B00/B10/B01/B11/B20/B02=', branch_error
   write (*, '(a,2(es12.4,1x))') '  plus/minus_operator=', operator_error
   write (*, '(a,es12.4)') '  mixed_LM_coverage=', coverage_error
   write (*, '(a)') 'UnitDresp09PauliDirect: PASS (independent Pauli projection; L=0..4; complex M; duality closure)'

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
      real(rp), allocatable :: g(:, :), gp(:, :), gpp(:, :), gpack(:), gdotpack(:), gddotpack(:)
      integer :: ispin, l

      potential = 0.0_rp
      call basis%initialize(size(mesh), lmax, 2)
      do ispin = 1, 2
         do l = 0, lmax
            call legacy_radial_fixture(nuclear_z, l, mesh_a, mesh_b, mesh, potential, energy, g, gp, gpp, 2.0_rp)
            gpack = reshape(g, [2*size(mesh)])
            gdotpack = reshape(gp, [2*size(mesh)])
            gddotpack = reshape(gpp, [2*size(mesh)])
            call basis%capture_channel(l, ispin, energy, mesh, potential, mesh_a, mesh_b, nuclear_z, &
               gpack, gdotpack, gddotpack, energy)
         end do
      end do
   end subroutine setup_radial_basis

   subroutine fill_compact(vector)
      complex(rp), intent(out) :: vector(:)
      integer :: i
      do i = 1, size(vector)
         vector(i) = cmplx(0.013_rp*real(i + 2, rp), -0.009_rp*real(3*i + 1, rp), rp)
      end do
   end subroutine fill_compact

   subroutine build_hamiltonian(matrix)
      complex(rp), intent(out) :: matrix(:, :)
      integer :: i, j
      do i = 1, size(matrix, 1)
         do j = 1, size(matrix, 2)
            matrix(i, j) = cmplx(0.007_rp*real(2*i + j, rp), -0.004_rp*real(i + 3*j, rp), rp)
            if (i == j) matrix(i, j) = matrix(i, j) + cmplx(0.21_rp*real(i, rp), 0.0_rp, rp)
         end do
      end do
   end subroutine build_hamiltonian

   subroutine reconstruct_dual_source(space, product, compact, source)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: compact(:)
      complex(rp), intent(out) :: source(:)
      integer :: flat, ir, raw_flat, site, response_l, response_m, mode
      type(response_super_index) :: item
      real(rp) :: weight

      source = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, product%product_dimension
         call product%unflatten_index(flat, site, response_l, response_m, mode)
         do ir = 2, space%npoint
            weight = space%radial_weights(ir)
            item%site = site
            item%response_l = response_l
            item%response_m = response_m
            item%radial_point = ir
            item%channel = 1
            call response_flatten_superindex(item, space%nsite, space%response_lmax, space%npoint, 1, raw_flat)
            source(raw_flat) = source(raw_flat) + product%blocks(site, response_l)%weighted_modes(ir, mode)* &
               conjg(compact(flat))/sqrt(weight)
         end do
      end do
   end subroutine reconstruct_dual_source

   subroutine compare_branches(vertices, compact, actual, errors, coverage)
      complex(rp), intent(in) :: vertices(:, :, :, :), compact(:), actual(:, :, :)
      real(rp), intent(inout) :: errors(:), coverage
      complex(rp) :: expected(size(actual, 1), size(actual, 2))
      integer :: mu, component

      do mu = 1, size(compact)
         if (abs(compact(mu)) == 0.0_rp) coverage = max(coverage, 1.0_rp)
      end do
      do component = 1, 6
         expected = cmplx(0.0_rp, 0.0_rp, rp)
         do mu = 1, size(compact)
            expected = expected + conjg(compact(mu))*transpose(conjg(vertices(:, :, component, mu)))
         end do
         errors(component) = max(errors(component), relative_error(actual(:, :, component), expected))
      end do
   end subroutine compare_branches

   pure real(rp) function relative_error(actual, expected) result(value)
      complex(rp), intent(in) :: actual(:, :), expected(:, :)
      value = sqrt(sum(abs(actual - expected)**2))/max(sqrt(sum(abs(expected)**2)), tiny(1.0_rp))
   end function relative_error

end program test_dresp09_pauli_direct
