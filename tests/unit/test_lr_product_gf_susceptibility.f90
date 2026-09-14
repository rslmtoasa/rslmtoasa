!------------------------------------------------------------------------------
! TDVK-03R compact reciprocal-GF representation and contraction test.
!
! The point-grid LR-GF-02 result is projected independently with the LR-04
! weighted-orthonormal map.  The compact GF result must agree with that
! projection; the compact Lehmann comparison below is diagnostic only.
!------------------------------------------------------------------------------
program test_lr_product_gf_susceptibility
   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use response_angular_basis_mod, only: response_gaunt
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex
   use lr_response_space_mod, only: response_space_layout
   use lr_pauli_transition_vertex_mod, only: pauli_endpoint_state, pauli_sigma_plus_matrix, &
      pauli_sigma_minus_matrix, pauli_vertex_capabilities, evaluate_pauli_transition_vertex
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, lmto_product_channel_plus, &
      lmto_product_channel_minus
   use lr_ks_susceptibility_mod, only: lr_electronic_state, lr_ks_susceptibility_request, &
      lr_ks_susceptibility_result, lr_product_ks_susceptibility_request, lr_product_ks_susceptibility_result, &
      lr_channel_plus, lr_channel_minus, lr_fermi_dirac_occupation, evaluate_lr_product_ks_susceptibility
   use lr_gf_susceptibility_mod, only: lr_gf_susceptibility_request, evaluate_lr_gf_susceptibility
   use lr_product_gf_susceptibility_mod, only: lr_product_gf_susceptibility_request, &
      lr_product_gf_susceptibility_result, lr_product_gf_contraction_scalar, &
      lr_product_gf_contraction_optimized, lr_product_gf_contraction_factorized, &
      evaluate_lr_product_gf_susceptibility, build_lr_product_gf_transition_amplitudes
   implicit none

   integer, parameter :: nr = 7, orbital_lmax = 1, response_lmax = 2, nsite = 1
   integer, parameter :: norb = (orbital_lmax + 1)**2, nbasis = 2*norb*nsite
   integer, parameter :: nbands = nbasis, nk = 2, nfrequency = 2
   real(rp), parameter :: mesh_a = 0.03_rp, mesh_b = 0.10_rp
   real(rp), parameter :: eta = 0.04_rp, temperature = 0.03_rp, energy_margin = 0.60_rp
   real(rp), parameter :: q_gamma(3) = [0.0_rp, 0.0_rp, 0.0_rp]
   real(rp), parameter :: q_finite(3) = [0.23_rp, 0.0_rp, 0.0_rp]
   real(rp), parameter :: frequencies(nfrequency) = [0.0_rp, 0.17_rp]
   real(rp), parameter :: tolerance = 1.0e-10_rp

   real(rp) :: radius(nr), eigenvalues(nbands, nk), k_points(3, nk), endpoint_points(3, nk)
   real(rp) :: k_weights(nk), occupations(nbands, nk)
   complex(rp) :: eigenvectors(nbasis, nbands, nk)
   type(lmto_radial_basis), target :: radial(nsite)
   type(response_space_layout), target :: space
   type(lmto_product_response_basis), target :: product_plus, product_minus
   type(lr_electronic_state), target :: state, endpoint
   logical :: failed
   real(rp) :: maximum_component_residual, maximum_projection_residual, maximum_transition_residual
   character(len=32) :: mode

   call get_command_argument(1, mode)
   call basis_init(orbital_lmax)
   call build_mesh(radius)
   call setup_radial_basis(radial, radius)
   call space%initialize(nsite, response_lmax, radius, mesh_a, mesh_b, 1)
   call build_electronic_fixture(eigenvalues, eigenvectors, k_points, k_weights, occupations)
   call state%initialize(eigenvalues, eigenvectors, k_points, k_weights, occupations, 0.0_rp, temperature)
   endpoint_points = k_points
   call endpoint%initialize(eigenvalues, eigenvectors, endpoint_points, k_weights, occupations, 0.0_rp, temperature)
   call product_plus%initialize(space, radial, lmto_product_channel_plus, .false.)
   call product_minus%initialize(space, radial, lmto_product_channel_minus, .false.)

   failed = .false.
   maximum_component_residual = 0.0_rp
   maximum_projection_residual = 0.0_rp
   maximum_transition_residual = 0.0_rp
   call component_vertex_checks(product_plus, state, failed, maximum_component_residual, maximum_transition_residual)
   call component_vertex_checks(product_minus, state, failed, maximum_component_residual, maximum_transition_residual)

   if (trim(mode) == 'reject_integration_eta') then
      call expect_integration_eta_guard(product_plus, state, endpoint)
      error stop 'UnitLrProductGfSusceptibility: integration eta guard unexpectedly returned'
   else if (len_trim(mode) > 0) then
      error stop 'UnitLrProductGfSusceptibility: unknown command argument'
   end if

   call run_case(product_plus, lr_channel_plus, q_gamma, q_gamma, 21, state, endpoint, failed, &
      maximum_projection_residual, .true.)
   call run_case(product_plus, lr_channel_plus, q_gamma, q_gamma, 41, state, endpoint, failed, &
      maximum_projection_residual, .false.)
   call run_case(product_minus, lr_channel_minus, q_finite, q_finite, 21, state, endpoint, failed, &
      maximum_projection_residual, .true.)

   if (maxval(abs(state%eigenvalues - eigenvalues)) > 0.0_rp .or. &
       maxval(abs(state%eigenvectors - eigenvectors)) > 0.0_rp .or. &
       maxval(abs(state%occupations - occupations)) > 0.0_rp .or. &
       maxval(abs(endpoint%eigenvalues - eigenvalues)) > 0.0_rp .or. &
       maxval(abs(endpoint%eigenvectors - eigenvectors)) > 0.0_rp) then
      failed = .true.
      write (*, '(a)') 'immutable electronic-state snapshot check: FAIL'
   end if

   write (*, '(a,es12.4)') 'maximum compact component/R1 residual = ', maximum_component_residual
   write (*, '(a,es12.4)') 'maximum GF transition-factorization residual = ', maximum_transition_residual
   write (*, '(a,es12.4)') 'maximum point-GF projection residual = ', maximum_projection_residual
   write (*, '(a,i0,a,i0)') 'product dimensions plus/minus = ', product_plus%product_dimension, '/', &
      product_minus%product_dimension
   if (failed) then
      write (*, '(a)') 'UnitLrProductGfSusceptibility: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrProductGfSusceptibility: PASS (component oracle, point-GF projection, channels, q/omega, Simpson controls)'

contains

   subroutine build_mesh(mesh)
      real(rp), intent(out) :: mesh(:)
      integer :: ir

      mesh(1) = 0.0_rp
      do ir = 2, size(mesh)
         mesh(ir) = mesh_b*(exp(mesh_a*real(ir - 1, rp)) - 1.0_rp)
      end do
   end subroutine build_mesh

   subroutine setup_radial_basis(bases, mesh)
      type(lmto_radial_basis), intent(out) :: bases(:)
      real(rp), intent(in) :: mesh(:)
      integer :: isite, ispin, l, ir

      do isite = 1, size(bases)
         call bases(isite)%initialize(size(mesh), orbital_lmax, 2)
         bases(isite)%rofi = mesh
         bases(isite)%mesh_a = mesh_a
         bases(isite)%mesh_b = mesh_b
         bases(isite)%channel_present = .true.
         bases(isite)%enu_work = 0.0_rp
         bases(isite)%gfac = 1.0_rp
         do ispin = 1, 2
            do l = 0, orbital_lmax
               do ir = 1, size(mesh)
                  bases(isite)%phi_large(ir, l + 1, ispin) = &
                     (1.0_rp + 0.17_rp*real(l, rp) + 0.06_rp*real(ispin, rp))* &
                     (1.0_rp + 0.31_rp*mesh(ir))
                  bases(isite)%phidot_large(ir, l + 1, ispin) = &
                     (0.11_rp + 0.03_rp*real(l, rp) + 0.02_rp*real(ispin, rp))* &
                     (1.0_rp + 0.19_rp*mesh(ir)**2)
               end do
            end do
         end do
      end do
   end subroutine setup_radial_basis

   subroutine build_electronic_fixture(values, vectors, points, weights, occ)
      real(rp), intent(out) :: values(:, :), points(:, :), weights(:), occ(:, :)
      complex(rp), intent(out) :: vectors(:, :, :)
      integer :: ib, ik

      values = 0.0_rp
      do ik = 1, nk
         do ib = 1, nbands
            values(ib, ik) = -0.72_rp + 0.19_rp*real(ib - 1, rp) + 0.037_rp*real(ik - 1, rp)
         end do
      end do
      vectors = cmplx(0.0_rp, 0.0_rp, rp)
      do ik = 1, nk
         do ib = 1, nbands
            vectors(ib, ib, ik) = cmplx(1.0_rp, 0.0_rp, rp)
         end do
      end do
      points = 0.0_rp
      weights = [0.35_rp, 0.65_rp]
      do ik = 1, nk
         do ib = 1, nbands
            occ(ib, ik) = lr_fermi_dirac_occupation(values(ib, ik), 0.0_rp, temperature)
         end do
      end do
   end subroutine build_electronic_fixture

   subroutine component_vertex_checks(product, electronic_state, failed, maximum_residual, maximum_transition_residual)
      type(lmto_product_response_basis), intent(in) :: product
      type(lr_electronic_state), intent(in) :: electronic_state
      logical, intent(inout) :: failed
      real(rp), intent(inout) :: maximum_residual
      real(rp), intent(inout) :: maximum_transition_residual

      complex(rp), allocatable :: vertices(:, :, :, :), reference(:), r1_coordinates(:), &
         reference_coordinates(:), contraction(:), transitions(:, :, :)
      type(pauli_endpoint_state) :: left_state, right_state
      type(pauli_vertex_capabilities) :: capabilities
      complex(rp) :: operator_matrix(2, 2), value
      integer, parameter :: left_plus(3) = [1, 2, 3], right_plus(3) = [5, 6, 7]
      integer, parameter :: left_minus(3) = [5, 6, 7], right_minus(3) = [1, 2, 3]
      integer :: pair, left_index, right_index, alpha, p, q, component
      real(rp) :: residual, scale, left_power, right_power

      call product%component_vertex_tensor(vertices)
      if (size(vertices, 1) /= nbasis .or. size(vertices, 2) /= nbasis .or. &
          size(vertices, 3) /= 4 .or. size(vertices, 4) /= product%product_dimension) then
         failed = .true.
      end if
      capabilities = pauli_vertex_capabilities()
      if (product%circular_channel == lmto_product_channel_plus) then
         operator_matrix = pauli_sigma_plus_matrix()
      else
         operator_matrix = pauli_sigma_minus_matrix()
      end if
      allocate(reference(space%ndim), r1_coordinates(product%product_dimension), &
         reference_coordinates(product%product_dimension), contraction(product%product_dimension), &
         transitions(nbands, nbands, product%product_dimension))
      call build_lr_product_gf_transition_amplitudes(vertices, electronic_state, 1, electronic_state, 1, transitions)
      do pair = 1, 3
         if (product%circular_channel == lmto_product_channel_plus) then
            left_index = left_plus(pair)
            right_index = right_plus(pair)
         else
            left_index = left_minus(pair)
            right_index = right_minus(pair)
         end if
         call left_state%initialize(electronic_state%eigenvalues(left_index, 1), &
            electronic_state%eigenvectors(:, left_index, 1))
         call right_state%initialize(electronic_state%eigenvalues(right_index, 1), &
            electronic_state%eigenvectors(:, right_index, 1))
         call evaluate_pauli_transition_vertex(space, [radial], left_state, right_state, operator_matrix, &
            capabilities, reference)
         call product%transition_coordinates(left_state, right_state, r1_coordinates)
         call product_coordinates_from_point(product, reference, reference_coordinates)
         contraction = build_component_contraction(product, vertices, left_state, right_state)
         scale = max(sqrt(sum(abs(reference_coordinates)**2)), epsilon(1.0_rp))
         residual = max(sqrt(sum(abs(contraction - reference_coordinates)**2)), &
            sqrt(sum(abs(contraction - r1_coordinates)**2)))/scale
         maximum_residual = max(maximum_residual, residual)
         if (residual >= tolerance) failed = .true.
         residual = sqrt(sum(abs(transitions(left_index, right_index, :) - r1_coordinates)**2))/scale
         maximum_transition_residual = max(maximum_transition_residual, residual)
         if (residual >= tolerance) failed = .true.
      end do
      deallocate(vertices, reference, r1_coordinates, reference_coordinates, contraction, transitions)
   end subroutine component_vertex_checks

   subroutine run_case(product, channel, q, endpoint_q, integration_points, state, endpoint, failed, maximum_residual, &
                       compare_contractions)
      type(lmto_product_response_basis), target, intent(in) :: product
      character(len=*), intent(in) :: channel
      real(rp), intent(in) :: q(3), endpoint_q(3)
      integer, intent(in) :: integration_points
      type(lr_electronic_state), target, intent(in) :: state
      type(lr_electronic_state), target, intent(inout) :: endpoint
      logical, intent(inout) :: failed
      real(rp), intent(inout) :: maximum_residual

      type(lr_ks_susceptibility_request) :: point_request
      type(lr_gf_susceptibility_request) :: point_gf_request
      type(lr_product_gf_susceptibility_request) :: product_gf_request
      type(lr_product_ks_susceptibility_request) :: lehmann_request
      type(lr_ks_susceptibility_result) :: point_gf_result
      type(lr_product_gf_susceptibility_result) :: product_gf_result
      type(lr_product_gf_susceptibility_result) :: optimized_gf_result
      type(lr_product_ks_susceptibility_result) :: lehmann_result
      complex(rp), allocatable :: reference(:, :, :)
      real(rp) :: absolute, reference_norm, relative, diagnostic_relative, d_inf
      logical, intent(in) :: compare_contractions
      type(lr_product_gf_susceptibility_result) :: scalar_gf_result

      call endpoint%initialize(state%eigenvalues, state%eigenvectors, &
         spread(endpoint_q, 2, nk), state%k_weights, state%occupations, state%fermi_level, state%temperature)

      point_request%q = q
      point_request%frequencies = frequencies
      point_request%eta = eta
      point_request%channel = channel
      point_request%response_space => space
      point_request%radial_bases => radial
      point_request%electronic_state => state
      point_request%q_endpoint_state => endpoint

      point_gf_request%q = q
      point_gf_request%frequencies = frequencies
      point_gf_request%eta = eta
      point_gf_request%channel = channel
      point_gf_request%integration_points = integration_points
      point_gf_request%integration_eta = 0.0_rp
      point_gf_request%energy_margin = energy_margin
      point_gf_request%response_space => space
      point_gf_request%radial_bases => radial
      point_gf_request%electronic_state => state
      point_gf_request%q_endpoint_state => endpoint
      call evaluate_lr_gf_susceptibility(point_gf_request, point_gf_result)

      product_gf_request%q = q
      product_gf_request%frequencies = frequencies
      product_gf_request%eta = eta
      product_gf_request%channel = channel
      product_gf_request%integration_points = integration_points
      product_gf_request%integration_eta = 0.0_rp
      product_gf_request%energy_margin = energy_margin
      product_gf_request%contraction_backend = lr_product_gf_contraction_factorized
      product_gf_request%product_basis => product
      product_gf_request%electronic_state => state
      product_gf_request%q_endpoint_state => endpoint
      call evaluate_lr_product_gf_susceptibility(product_gf_request, product_gf_result)

      if (product_gf_result%actual_integration_eta /= eta/40.0_rp .or. &
          product_gf_result%point_response_allocated .or. product_gf_result%product_dimension /= product%product_dimension .or. &
          trim(product_gf_result%response_representation) /= 'weighted-orthonormal LMTO product representation') then
         failed = .true.
      end if
      if (trim(product_gf_result%contraction_backend) /= lr_product_gf_contraction_factorized .or. &
          product_gf_result%energy_spacing <= 0.0_rp .or. product_gf_result%spacing_over_integration_eta <= 0.0_rp) then
         failed = .true.
      end if
      if (compare_contractions) then
         product_gf_request%contraction_backend = lr_product_gf_contraction_optimized
         call evaluate_lr_product_gf_susceptibility(product_gf_request, optimized_gf_result)
         product_gf_request%contraction_backend = lr_product_gf_contraction_scalar
         call evaluate_lr_product_gf_susceptibility(product_gf_request, scalar_gf_result)
         call report_contraction_compare('scalar/optimized', scalar_gf_result, optimized_gf_result, failed)
         call report_contraction_compare('scalar/factorized', scalar_gf_result, product_gf_result, failed)
         call report_contraction_compare('optimized/factorized', optimized_gf_result, product_gf_result, failed)
         if (trim(scalar_gf_result%contraction_backend) /= lr_product_gf_contraction_scalar .or. &
             trim(optimized_gf_result%contraction_backend) /= lr_product_gf_contraction_optimized) then
            failed = .true.
         end if
      end if
      allocate(reference(product%product_dimension, product%product_dimension, nfrequency))
      call project_point_gf_to_product(product, point_gf_result%susceptibility, reference)
      call report_matrix_difference(product_gf_result%susceptibility, reference, absolute, reference_norm, relative)
      maximum_residual = max(maximum_residual, relative)
      write (*, '(a,a,a,i0,a,es12.4,a,es12.4,a,es12.4)') 'GF projection channel=', trim(channel), &
         ' integration_points=', integration_points, ' abs=', absolute, ' ref=', reference_norm, ' rel=', relative
      write (*, '(a,a,a,i0,a,i0,a,i0,a,i0)') 'GF compact memory channel=', trim(channel), &
         ' product_dim=', product_gf_result%product_dimension, ' component=', &
         product_gf_result%component_vertex_memory_bytes, ' gf=', product_gf_result%gf_matrix_memory_bytes, &
         ' susceptibility=', product_gf_result%susceptibility_memory_bytes
      if (relative >= tolerance) failed = .true.

      ! Independent compact Lehmann diagnostic.  It is reported but has no
      ! threshold and cannot affect the representation PASS decision.
      lehmann_request%q = q
      lehmann_request%frequencies = frequencies
      lehmann_request%eta = eta
      lehmann_request%channel = channel
      lehmann_request%product_basis => product
      lehmann_request%electronic_state => state
      lehmann_request%q_endpoint_state => endpoint
      call evaluate_lr_product_ks_susceptibility(lehmann_request, lehmann_result)
      call report_matrix_difference(lehmann_result%susceptibility, product_gf_result%susceptibility, &
         absolute, reference_norm, diagnostic_relative)
      d_inf = maxval(abs(lehmann_result%susceptibility - product_gf_result%susceptibility))
      write (*, '(a,a,a,es12.4,a,es12.4,a,es12.4)') 'Lehmann/GF diagnostic channel=', trim(channel), &
         ' dF=', absolute, ' rF=', diagnostic_relative, ' dInf=', d_inf
      deallocate(reference)
   end subroutine run_case

   subroutine report_contraction_compare(label, lhs, rhs, failed)
      character(len=*), intent(in) :: label
      type(lr_product_gf_susceptibility_result), intent(in) :: lhs, rhs
      logical, intent(inout) :: failed
      real(rp) :: lhs_norm, rhs_norm, difference, relative, difference_infinity, speedup

      lhs_norm = sqrt(sum(abs(lhs%susceptibility)**2))
      rhs_norm = sqrt(sum(abs(rhs%susceptibility)**2))
      difference = sqrt(sum(abs(lhs%susceptibility - rhs%susceptibility)**2))
      relative = difference/max(lhs_norm, rhs_norm, epsilon(1.0_rp))
      difference_infinity = maxval(abs(lhs%susceptibility - rhs%susceptibility))
      speedup = lhs%wall_time_seconds/max(rhs%wall_time_seconds, epsilon(1.0_rp))
      write (*, '(a,a,a,es16.8,a,es16.8,a,es16.8,a,es16.8,a,es16.8,a,es16.8,a,es16.8,a,es16.8)') &
         'GF three-way ', trim(label), ' norm_lhs=', lhs_norm, ' norm_rhs=', rhs_norm, &
         ' dF=', difference, ' rF=', relative, ' dInf=', difference_infinity, &
         ' wall_lhs=', lhs%wall_time_seconds, ' wall_rhs=', rhs%wall_time_seconds, ' lhs_over_rhs=', speedup
      if (relative >= tolerance .or. difference_infinity >= tolerance*max(1.0_rp, lhs_norm, rhs_norm)) failed = .true.
   end subroutine report_contraction_compare

   subroutine project_point_gf_to_product(product, canonical, projected)
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: canonical(:, :, :)
      complex(rp), intent(out) :: projected(:, :, :)
      complex(rp), allocatable :: b(:), y(:)
      integer :: frequency, alpha, beta, flat, site, response_l, response_m, product_mode, ir
      type(response_super_index) :: item
      real(rp) :: weight

      allocate(b(space%ndim), y(space%ndim))
      projected = cmplx(0.0_rp, 0.0_rp, rp)
      do frequency = 1, size(canonical, 3)
         do beta = 1, product%product_dimension
            b = cmplx(0.0_rp, 0.0_rp, rp)
            call product%unflatten_index(beta, site, response_l, response_m, product_mode)
            do ir = 1, space%npoint
               item = response_super_index(site, response_l, response_m, ir, 1)
               flat = item_to_flat(item)
               weight = space%metric_weights(flat)
               if (weight > 0.0_rp) b(flat) = product%blocks(site, response_l)%weighted_modes(ir, product_mode)/sqrt(weight)
            end do
            y = matmul(canonical(:, :, frequency), b)
            do alpha = 1, product%product_dimension
               projected(alpha, beta, frequency) = product_inner_point(product, alpha, y)
            end do
         end do
      end do
      deallocate(b, y)
   end subroutine project_point_gf_to_product

   function product_inner_point(product, alpha, vector) result(value)
      type(lmto_product_response_basis), intent(in) :: product
      integer, intent(in) :: alpha
      complex(rp), intent(in) :: vector(:)
      complex(rp) :: value
      integer :: site, response_l, response_m, product_mode, ir, flat
      type(response_super_index) :: item
      real(rp) :: weight

      call product%unflatten_index(alpha, site, response_l, response_m, product_mode)
      value = cmplx(0.0_rp, 0.0_rp, rp)
      do ir = 1, space%npoint
         item = response_super_index(site, response_l, response_m, ir, 1)
         flat = item_to_flat(item)
         weight = space%metric_weights(flat)
         if (weight > 0.0_rp) value = value + conjg(product%blocks(site, response_l)%weighted_modes(ir, product_mode))* &
            sqrt(weight)*vector(flat)
      end do
   end function product_inner_point

   integer function item_to_flat(item) result(flat)
      type(response_super_index), intent(in) :: item
      call response_flatten_superindex(item, space%nsite, space%response_lmax, space%npoint, space%nchannel, flat)
   end function item_to_flat

   subroutine report_matrix_difference(actual, reference, absolute, reference_norm, relative)
      complex(rp), intent(in) :: actual(:, :, :), reference(:, :, :)
      real(rp), intent(out) :: absolute, reference_norm, relative

      absolute = sqrt(sum(abs(actual - reference)**2))
      reference_norm = sqrt(sum(abs(reference)**2))
      relative = absolute/max(reference_norm, epsilon(1.0_rp))
   end subroutine report_matrix_difference

   subroutine expect_integration_eta_guard(product, state, endpoint)
      type(lmto_product_response_basis), target, intent(in) :: product
      type(lr_electronic_state), target, intent(in) :: state, endpoint
      type(lr_product_gf_susceptibility_request) :: request
      type(lr_product_gf_susceptibility_result) :: result

      request%frequencies = [0.0_rp]
      request%eta = eta
      request%integration_points = 3
      request%integration_eta = eta
      request%product_basis => product
      request%electronic_state => state
      request%q_endpoint_state => endpoint
      call evaluate_lr_product_gf_susceptibility(request, result)
   end subroutine expect_integration_eta_guard

   ! The following two helpers are retained as small compile-time/runtime
   ! guards for the component-ordering contract.  The production evaluator
   ! itself never calls the Lehmann route or constructs point response data.
   function build_component_contraction(product, vertices, left_state, right_state) result(coordinates)
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: vertices(:, :, :, :)
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp) :: coordinates(product%product_dimension)
      integer :: alpha, p, q, component
      real(rp) :: left_power, right_power

      coordinates = cmplx(0.0_rp, 0.0_rp, rp)
      do alpha = 1, product%product_dimension
         do p = 0, 1
            do q = 0, 1
               component = 1 + p + 2*q
               left_power = merge(left_state%energy, 1.0_rp, p == 1)
               right_power = merge(right_state%energy, 1.0_rp, q == 1)
               coordinates(alpha) = coordinates(alpha) + left_power*right_power*sum( &
                  conjg(left_state%coefficients)*matmul(vertices(:, :, component, alpha), right_state%coefficients))
            end do
         end do
      end do
   end function build_component_contraction

   subroutine product_coordinates_from_point(product, point_vector, coordinates)
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: point_vector(:)
      complex(rp), intent(out) :: coordinates(:)
      integer :: alpha

      do alpha = 1, product%product_dimension
         coordinates(alpha) = product_inner_point(product, alpha, point_vector)
      end do
   end subroutine product_coordinates_from_point

end program test_lr_product_gf_susceptibility
