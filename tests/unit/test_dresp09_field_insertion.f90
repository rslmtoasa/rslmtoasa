! DRESP-09 representation-level field insertion regression.
!
! The fixture covers every retained response harmonic L=0..4, complex raw and
! compact fields, the metric pairing, the compact adjoint, and the reversal of
! the endpoint-energy powers.  The scalar-relativistic direct route is tested
! only through its fail-closed capability boundary: its mixed-spin lower-
! component metric is not certified by the accepted production contract.
program test_dresp09_field_insertion
   use precision_mod, only: rp
   use response_angular_basis_mod, only: response_gaunt
   use lr_response_space_mod, only: response_space_layout
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis
   use lr_dresp09_field_insertion_mod, only: dresp09_direct_route_capability, &
      dresp09_compact_field_from_raw, dresp09_raw_field_from_compact, dresp09_compact_pairing, &
      dresp09_raw_pairing, dresp09_metric_adjoint_operator, dresp09_source_operator, &
      dresp09_adjoint_density_coordinates, dresp09_trace_density_operator, &
      dresp09_product_metric_adjoint_operator, dresp09_require_direct_route
   implicit none

   integer, parameter :: npoint = 5, response_lmax = 4, nproduct = (response_lmax + 1)**2, n = 4
   real(rp), parameter :: tolerance = 5.0e-13_rp
   type(response_space_layout) :: space
   type(lmto_product_response_basis) :: product
   type(dresp09_direct_route_capability) :: capability
   complex(rp), allocatable :: field(:), density(:), field_rebuilt(:), density_rebuilt(:)
   complex(rp), allocatable :: field_rebuilt_projection(:), density_rebuilt_projection(:)
   complex(rp), allocatable :: transition(:, :, :), operator(:, :), source(:, :), rho(:, :), coordinates(:)
   complex(rp), allocatable :: components(:, :, :, :), product_operator(:, :)
   real(rp) :: radius(npoint), max_projection_error, metric_error, trace_error, source_error, product_error
   real(rp) :: angular_error, expected_gaunt
   integer :: l, m, i, j, mu, p, q, component
   real(rp) :: norm
   character(len=32) :: command_mode

   command_mode = ''
   if (command_argument_count() > 0) call get_command_argument(1, command_mode)
   if (trim(command_mode) == 'blocked') call dresp09_require_direct_route(capability)

   radius = [0.0_rp, 0.07_rp, 0.16_rp, 0.29_rp, 0.45_rp]
   call space%initialize(1, response_lmax, radius, 0.03_rp, 0.10_rp, 1)
   call build_fixture_product(product, space)
   angular_error = 0.0_rp
   do l = 0, 2
      do m = -l, l
         expected_gaunt = 1.0_rp/sqrt(4.0_rp*acos(-1.0_rp))
         angular_error = max(angular_error, abs(response_gaunt(l, m, l, m, 0, 0) - expected_gaunt))
         if (m /= 0) angular_error = max(angular_error, abs(response_gaunt(l, m, l, -m, 0, 0)))
      end do
   end do
   angular_error = max(angular_error, abs(response_gaunt(0, 0, 0, 0, 1, 0)))
   angular_error = max(angular_error, abs(response_gaunt(0, 0, 1, 0, 0, 1)))

   allocate(field(nproduct), density(nproduct), field_rebuilt(space%ndim), density_rebuilt(space%ndim), &
      field_rebuilt_projection(nproduct), density_rebuilt_projection(nproduct))
   do mu = 1, nproduct
      field(mu) = cmplx(0.13_rp*real(mu, rp), -0.017_rp*real(2*mu + 1, rp), rp)
      density(mu) = cmplx(-0.021_rp*real(mu + 2, rp), 0.031_rp*real(mu, rp), rp)
   end do
   call dresp09_raw_field_from_compact(space, product, field, field_rebuilt)
   call dresp09_raw_field_from_compact(space, product, density, density_rebuilt)
   call dresp09_compact_field_from_raw(space, product, field_rebuilt, field_rebuilt_projection)
   call dresp09_compact_field_from_raw(space, product, density_rebuilt, density_rebuilt_projection)

   max_projection_error = maxval(abs(field - field_rebuilt_projection))
   max_projection_error = max(max_projection_error, maxval(abs(density - density_rebuilt_projection)))
   metric_error = abs(dresp09_compact_pairing(field_rebuilt_projection, density_rebuilt_projection) - &
      dresp09_raw_pairing(space, field_rebuilt, density_rebuilt))

   allocate(transition(n, n, nproduct), operator(n, n), source(n, n), rho(n, n), coordinates(nproduct))
   do mu = 1, nproduct
      do i = 1, n
         do j = 1, n
            transition(i, j, mu) = cmplx(0.003_rp*real(mu + 2*i + j, rp), &
               -0.002_rp*real(2*mu + i + 3*j, rp), rp)
         end do
      end do
   end do
   do i = 1, n
      do j = 1, n
         rho(i, j) = cmplx(0.017_rp*real(i + 2*j, rp), -0.011_rp*real(2*i + j, rp), rp)
      end do
   end do
   call dresp09_adjoint_density_coordinates(transition, rho, coordinates)
   call dresp09_metric_adjoint_operator(transition, field, operator)
   trace_error = abs(dresp09_trace_density_operator(rho, operator) - &
      dresp09_compact_pairing(field, coordinates))
   call dresp09_source_operator(transition, field, source)
   source_error = abs(dresp09_trace_density_operator(rho, source) - sum(field*coordinates))

   allocate(components(n, n, 4, nproduct), product_operator(n, n))
   do mu = 1, nproduct
      do component = 1, 4
         do i = 1, n
            do j = 1, n
               components(i, j, component, mu) = transition(i, j, mu)* &
                  cmplx(0.07_rp*real(component + i + j, rp), 0.0_rp, rp)
            end do
         end do
      end do
   end do
   do i = 1, n
      do j = 1, n
         if (i == j) then
            components(i, j, :, :) = components(i, j, :, :) + cmplx(0.23_rp*real(i, rp), 0.0_rp, rp)
         end if
      end do
   end do
   call dresp09_product_metric_adjoint_operator(components, diagonal_hamiltonian(), field, product_operator)
   product_error = product_adjoint_reference(components, diagonal_hamiltonian(), field, product_operator)

   if (max_projection_error > tolerance .or. metric_error > tolerance .or. trace_error > tolerance .or. &
       source_error > tolerance .or. product_error > tolerance .or. angular_error > tolerance .or. &
       capability%certified()) then
      write (*, '(a)') 'UnitDresp09FieldInsertion: FAIL'
      write (*, '(a,es12.4)') '  compact_projection_max=', max_projection_error
      write (*, '(a,es12.4)') '  compact_vs_raw_pairing=', metric_error
      write (*, '(a,es12.4)') '  metric_adjoint_trace=', trace_error
      write (*, '(a,es12.4)') '  source_trace=', source_error
      write (*, '(a,es12.4)') '  product_adjoint=', product_error
      write (*, '(a,es12.4)') '  angular_oracle=', angular_error
      error stop 1
   end if
   write (*, '(a,es12.4)') '  compact_projection_max=', max_projection_error
   write (*, '(a,es12.4)') '  compact_vs_raw_pairing=', metric_error
   write (*, '(a,es12.4)') '  metric_adjoint_trace=', trace_error
   write (*, '(a,es12.4)') '  source_trace=', source_error
   write (*, '(a,es12.4)') '  product_adjoint=', product_error
   write (*, '(a,es12.4)') '  angular_oracle=', angular_error
   write (*, '(a)') '  direct_route = BLOCKED - MIXED-SPIN RADIAL METRIC NOT CERTIFIED'
   write (*, '(a)') 'UnitDresp09FieldInsertion: PASS (compact algebra; direct SR route fail-closed)'

contains

   subroutine build_fixture_product(product, space)
      type(lmto_product_response_basis), intent(out) :: product
      type(response_space_layout), intent(in) :: space
      integer :: l, ir
      real(rp) :: local_norm

      product%nsite = 1
      product%orbital_lmax = 2
      product%response_lmax = response_lmax
      product%circular_channel = 1
      product%npoint = space%npoint
      product%product_dimension = nproduct
      allocate(product%blocks(1, 0:response_lmax))
      do l = 0, response_lmax
         product%blocks(1, l)%site = 1
         product%blocks(1, l)%response_l = l
         product%blocks(1, l)%channel = 1
         product%blocks(1, l)%npoint = space%npoint
         product%blocks(1, l)%rank = 1
         allocate(product%blocks(1, l)%weighted_modes(space%npoint, 1))
         do ir = 1, space%npoint
            if (ir == 1) then
               product%blocks(1, l)%weighted_modes(ir, 1) = cmplx(0.0_rp, 0.0_rp, rp)
            else
               product%blocks(1, l)%weighted_modes(ir, 1) = cmplx(real(ir + 2*l, rp), &
                  0.1_rp*real(ir - l, rp), rp)
            end if
         end do
         local_norm = sqrt(sum(abs(product%blocks(1, l)%weighted_modes(:, 1))**2))
         product%blocks(1, l)%weighted_modes(:, 1) = product%blocks(1, l)%weighted_modes(:, 1)/local_norm
      end do
   end subroutine build_fixture_product

   function diagonal_hamiltonian() result(matrix)
      complex(rp) :: matrix(n, n)
      integer :: index

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do index = 1, n
         matrix(index, index) = cmplx(-0.31_rp + 0.19_rp*real(index, rp), 0.0_rp, rp)
      end do
   end function diagonal_hamiltonian

   function product_adjoint_reference(vertices, hamiltonian, field, actual) result(error)
      complex(rp), intent(in) :: vertices(:, :, :, :), hamiltonian(:, :), field(:), actual(:, :)
      real(rp) :: error
      complex(rp) :: expected(n, n), term(n, n), coefficient
      integer :: mu, p, q, component

      expected = cmplx(0.0_rp, 0.0_rp, rp)
      do mu = 1, size(field)
         do p = 0, 1
            do q = 0, 1
               component = 1 + p + 2*q
               term = transpose(conjg(vertices(:, :, component, mu)))
               if (q == 1) term = matmul(hamiltonian, term)
               if (p == 1) term = matmul(term, hamiltonian)
               coefficient = conjg(field(mu))
               expected = expected + coefficient*term
            end do
         end do
      end do
      error = sqrt(sum(abs(expected - actual)**2))/max(sqrt(sum(abs(expected)**2)), tiny(1.0_rp))
   end function product_adjoint_reference

end program test_dresp09_field_insertion
