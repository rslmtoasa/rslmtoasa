!------------------------------------------------------------------------------
! LR-06 spectral/Lehmann transverse KS susceptibility oracles.
!
! The reference accumulation below is deliberately independent of both the
! production susceptibility and production transition-vertex routines.  It
! uses a two-site s-only coefficient fixture, for which the Pauli transition
! vector is known analytically.
!------------------------------------------------------------------------------
program test_lr_ks_susceptibility
   use precision_mod, only: rp
   use response_angular_basis_mod, only: response_angular_pi, response_lm_index
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex, &
      response_unflatten_superindex
   use lr_response_space_mod, only: response_space_layout
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_ks_susceptibility_mod, only: lr_electronic_state, lr_ks_susceptibility_request, &
      lr_ks_susceptibility_result, lr_channel_plus, lr_channel_minus, &
      evaluate_lr_ks_susceptibility, evaluate_lr_static_residual
   implicit none

   integer, parameter :: nr = 9, nsite = 2, orbital_lmax = 1
   integer, parameter :: norb = (orbital_lmax + 1)**2, nbasis = 2*nsite*norb
   integer, parameter :: nbands = 4, nk = 2, nfrequency = 3
   real(rp), parameter :: mesh_a = 0.03_rp, mesh_b = 0.10_rp
   real(rp), parameter :: eta = 0.04_rp, tolerance = 1.0e-8_rp
   real(rp), parameter :: frequencies(nfrequency) = [-0.23_rp, 0.0_rp, 0.27_rp]
   real(rp), parameter :: eigenvalues(nbands, nk) = reshape([ &
      -0.40_rp, -0.40_rp, 0.30_rp, 0.30_rp, &
      -0.40_rp, -0.40_rp, 0.30_rp, 0.30_rp], [nbands, nk])
   real(rp), parameter :: occupations(nbands, nk) = reshape([ &
      1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, &
      1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp], [nbands, nk])
   real(rp), parameter :: weights(nk) = [1.0_rp, 3.0_rp]
   real(rp), parameter :: k_points(3, nk) = reshape([ &
      0.0_rp, 0.0_rp, 0.0_rp, 0.20_rp, 0.0_rp, 0.0_rp], [3, nk])

   real(rp) :: radius(nr), off_mesh_points(3, nk)
   complex(rp) :: eigenvectors(nbasis, nbands, nk), rotated_vectors(nbasis, nbands, nk)
   complex(rp) :: reversed_vectors(nbasis, nbands, nk)
   type(lmto_radial_basis), target :: radial(nsite)
   type(response_space_layout), target :: space
   type(lr_electronic_state), target :: state, endpoint, off_mesh_endpoint, rotated_state, rotated_endpoint, &
      reversed_state, reversed_endpoint, negative_endpoint
   type(lr_ks_susceptibility_request) :: request
   type(lr_ks_susceptibility_result) :: result, minus_result, reversed_result, &
      rotated_result, off_mesh_result
   complex(rp), allocatable :: reference(:, :, :), field(:), magnetization(:)
   real(rp) :: absolute_residual, relative_residual
   logical :: failed

   call build_mesh(radius)
   call setup_radial_basis(radial, radius)
   call space%initialize(nsite, 2*orbital_lmax, radius, mesh_a, mesh_b, 1)
   call build_coefficients(eigenvectors)
   call state%initialize(eigenvalues, eigenvectors, k_points, weights, occupations, 0.0_rp, 0.0_rp)
   call endpoint%initialize(eigenvalues, eigenvectors, k_points, weights, occupations, 0.0_rp, 0.0_rp)
   failed = .false.

   request%q = [0.0_rp, 0.0_rp, 0.0_rp]
   request%frequencies = frequencies
   request%eta = eta
   request%channel = lr_channel_plus
   request%response_space => space
   request%radial_bases => radial
   request%electronic_state => state
   request%q_endpoint_state => endpoint
   call evaluate_lr_ks_susceptibility(request, result)

   allocate(reference(space%ndim, space%ndim, nfrequency))
   call independent_explicit_loop(space, state, endpoint, request, reference)
   call compare_response(result%susceptibility, reference, tolerance, 'explicit-loop production oracle', failed)

   call analytic_two_level_check(space, state, result, failed)

   ! q=0 and an off-mesh q whose second endpoint crosses the half-open BZ edge.
   off_mesh_points(:, 1) = [0.37_rp, 0.0_rp, 0.0_rp]
   off_mesh_points(:, 2) = [-0.43_rp, 0.0_rp, 0.0_rp]
   call off_mesh_endpoint%initialize(eigenvalues, eigenvectors, off_mesh_points, weights, occupations, 0.0_rp, 0.0_rp)
   request%q = [0.37_rp, 0.0_rp, 0.0_rp]
   request%q_endpoint_state => off_mesh_endpoint
   call evaluate_lr_ks_susceptibility(request, off_mesh_result)
   if (maxval(abs(off_mesh_result%q - request%q)) > tolerance) failed = .true.
   write (*, '(a)') 'q=0 and off-mesh/folded k+q support: PASS'

   off_mesh_points(:, 1) = [-0.37_rp, 0.0_rp, 0.0_rp]
   off_mesh_points(:, 2) = [-0.17_rp, 0.0_rp, 0.0_rp]
   call negative_endpoint%initialize(eigenvalues, eigenvectors, off_mesh_points, weights, occupations, 0.0_rp, 0.0_rp)
   request%q = [0.37_rp, 0.0_rp, 0.0_rp]
   request%q_endpoint_state => off_mesh_endpoint
   request%frequencies = [frequencies(3)]
   request%channel = lr_channel_plus
   call evaluate_lr_ks_susceptibility(request, off_mesh_result)
   request%q = [-0.37_rp, 0.0_rp, 0.0_rp]
   request%q_endpoint_state => negative_endpoint
   request%frequencies = [-frequencies(3)]
   request%channel = lr_channel_minus
   call evaluate_lr_ks_susceptibility(request, minus_result)
   call compare_response(off_mesh_result%susceptibility, conjg(minus_result%susceptibility), tolerance, &
                         'q/-q retarded/advanced covariance', failed)

   ! Retarded/advanced covariance with circular-sector exchange from LR-03.
   request%q = [0.0_rp, 0.0_rp, 0.0_rp]
   request%q_endpoint_state => endpoint
   request%frequencies = [frequencies(3)]
   request%channel = lr_channel_plus
   call evaluate_lr_ks_susceptibility(request, result)
   request%frequencies = [-frequencies(3)]
   request%channel = lr_channel_minus
   call evaluate_lr_ks_susceptibility(request, minus_result)
   call compare_response(result%susceptibility, conjg(minus_result%susceptibility), tolerance, &
                         'retarded/advanced circular covariance', failed)

   ! Global spin reversal exchanges chi+ and chi-.  The radial fixture is spin
   ! symmetric, so no radial relabeling is hidden in this oracle.
   call spin_reverse_coefficients(eigenvectors, reversed_vectors)
   call reversed_state%initialize(eigenvalues, reversed_vectors, k_points, weights, occupations, 0.0_rp, 0.0_rp)
   call reversed_endpoint%initialize(eigenvalues, reversed_vectors, k_points, weights, occupations, 0.0_rp, 0.0_rp)
   request%frequencies = [frequencies(3)]
   request%channel = lr_channel_plus
   request%electronic_state => reversed_state
   request%q_endpoint_state => reversed_endpoint
   call evaluate_lr_ks_susceptibility(request, reversed_result)
   request%electronic_state => state
   request%q_endpoint_state => endpoint
   request%channel = lr_channel_minus
   call evaluate_lr_ks_susceptibility(request, minus_result)
   call compare_response(reversed_result%susceptibility, minus_result%susceptibility, tolerance, &
                         'global spin-reversal circular exchange', failed)

   ! Rotate both complete degenerate occupied/empty subspaces.  The summed
   ! susceptibility must not depend on the selected degenerate representatives.
   call rotate_degenerate_subspaces(eigenvectors, rotated_vectors)
   call rotated_state%initialize(eigenvalues, rotated_vectors, k_points, weights, occupations, 0.0_rp, 0.0_rp)
   call rotated_endpoint%initialize(eigenvalues, rotated_vectors, k_points, weights, occupations, 0.0_rp, 0.0_rp)
   request%frequencies = [frequencies(3)]
   request%channel = lr_channel_plus
   request%electronic_state => rotated_state
   request%q_endpoint_state => rotated_endpoint
   call evaluate_lr_ks_susceptibility(request, rotated_result)
   request%electronic_state => state
   request%q_endpoint_state => endpoint
   call evaluate_lr_ks_susceptibility(request, result)
   call compare_response(rotated_result%susceptibility, result%susceptibility, tolerance, &
                         'degenerate-subspace invariance', failed)

   ! The static identity diagnostic is intentionally raw: field and target are
   ! independent stored vectors and no correction is applied to chiKS.
   allocate(field(space%ndim), magnetization(space%ndim))
   call fill_diagnostic_vectors(field, magnetization)
   request%frequencies = [0.0_rp]
   request%channel = lr_channel_plus
   request%electronic_state => state
   request%q_endpoint_state => endpoint
   call evaluate_lr_ks_susceptibility(request, result)
   call evaluate_lr_static_residual(space, result%susceptibility(:, :, 1), field, magnetization, &
                                    absolute_residual, relative_residual)
   write (*, '(a,es16.8,a,es16.8)') 'Raw q=0 static residual (absolute, relative) = ', absolute_residual, ' ', &
      relative_residual
   write (*, '(a)') 'Raw static residual: diagnostic only; susceptibility unchanged'

   if (failed) then
      write (*, '(a)') 'UnitLrKsSusceptibility: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrKsSusceptibility: PASS (Lehmann, explicit loop, covariance, degeneracy)'

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
      integer :: isite, ispin

      do isite = 1, size(bases)
         call bases(isite)%initialize(size(mesh), orbital_lmax, 2)
         bases(isite)%rofi = mesh
         bases(isite)%mesh_a = mesh_a
         bases(isite)%mesh_b = mesh_b
         bases(isite)%channel_present = .true.
         bases(isite)%enu_work = 0.0_rp
         bases(isite)%phi_large = 0.0_rp
         bases(isite)%phidot_large = 0.0_rp
         bases(isite)%phiddot_large = 0.0_rp
         bases(isite)%gfac = 1.0_rp
         do ispin = 1, 2
            bases(isite)%phi_large(:, 1, ispin) = mesh
            bases(isite)%phi_large(:, 2, ispin) = mesh**2
         end do
      end do
   end subroutine setup_radial_basis

   subroutine build_coefficients(coefficients)
      complex(rp), intent(out) :: coefficients(:, :, :)
      integer :: ik

      coefficients = cmplx(0.0_rp, 0.0_rp, rp)
      do ik = 1, nk
         ! One occupied up and one empty down s state per site.  The two
         ! sites make the degenerate-subspace oracle nontrivial.
         coefficients(1, 1, ik) = cmplx(1.0_rp, 0.0_rp, rp)
         coefficients(5, 2, ik) = cmplx(1.0_rp, 0.0_rp, rp)
         coefficients(2, 3, ik) = cmplx(1.0_rp, 0.0_rp, rp)
         coefficients(6, 4, ik) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
   end subroutine build_coefficients

   subroutine independent_explicit_loop(space, left, right, request, canonical)
      type(response_space_layout), intent(in) :: space
      type(lr_electronic_state), intent(in) :: left, right
      type(lr_ks_susceptibility_request), intent(in) :: request
      complex(rp), intent(out) :: canonical(:, :, :)
      real(rp) :: metric(space%ndim), denominator_real, weight_sum, df
      complex(rp) :: transition(space%ndim), pair_factor
      integer :: ik, ib, jb, ifrequency, i, j

      call response_metric_for_test(space, metric)
      canonical = cmplx(0.0_rp, 0.0_rp, rp)
      weight_sum = sum(left%k_weights)
      do ik = 1, left%nk
         do ib = 1, left%nbands
            do jb = 1, right%nbands
               df = left%occupations(ib, ik) - right%occupations(jb, ik)
               if (df == 0.0_rp) cycle
               call independent_s_only_plus_transition(space, left, right, ik, ib, jb, transition)
               do ifrequency = 1, size(request%frequencies)
                  denominator_real = request%frequencies(ifrequency) + left%eigenvalues(ib, ik) - &
                                     right%eigenvalues(jb, ik)
                  pair_factor = cmplx(2.0_rp*left%k_weights(ik)/weight_sum, 0.0_rp, rp)*df/ &
                                cmplx(denominator_real, request%eta, rp)
                  do j = 1, space%ndim
                     do i = 1, space%ndim
                        canonical(i, j, ifrequency) = canonical(i, j, ifrequency) + &
                           pair_factor*transition(i)*conjg(transition(j))*metric(j)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine independent_explicit_loop

   subroutine independent_s_only_plus_transition(space, left, right, ik, ib, jb, transition)
      type(response_space_layout), intent(in) :: space
      type(lr_electronic_state), intent(in) :: left, right
      integer, intent(in) :: ik, ib, jb
      complex(rp), intent(out) :: transition(:)
      type(response_super_index) :: item
      integer :: flat, left_index, right_index
      real(rp) :: radial_factor

      transition = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, &
                                            space%nchannel, item)
         if (item%response_l /= 0 .or. item%response_m /= 0) cycle
         left_index = (item%site - 1)*2*norb + 1
         right_index = left_index + norb
         if (item%radial_point == 1) then
            radial_factor = s_only_origin_factor(space%radius)
         else
            radial_factor = 1.0_rp
         end if
         transition(flat) = conjg(left%eigenvectors(left_index, ib, ik))* &
            right%eigenvectors(right_index, jb, ik)*radial_factor/ &
            sqrt(4.0_rp*response_angular_pi)
      end do
   end subroutine independent_s_only_plus_transition

   real(rp) function s_only_origin_factor(mesh) result(value)
      real(rp), intent(in) :: mesh(:)
      real(rp) :: first, second

      first = real(mesh(2)*mesh(2), rp)/(mesh(2)**2)
      second = real(mesh(3)*mesh(3), rp)/(mesh(3)**2)
      value = (first*mesh(3)**2 - second*mesh(2)**2)/(mesh(3)**2 - mesh(2)**2)
   end function s_only_origin_factor

   subroutine response_metric_for_test(space, metric)
      type(response_space_layout), intent(in) :: space
      real(rp), intent(out) :: metric(:)

      if (size(metric) /= space%ndim) error stop 'response_metric_for_test: shape mismatch'
      metric = space%metric_weights
   end subroutine response_metric_for_test

   subroutine analytic_two_level_check(space, left, actual, failed)
      type(response_space_layout), intent(in) :: space
      type(lr_electronic_state), intent(in) :: left
      type(lr_ks_susceptibility_result), intent(in) :: actual
      logical, intent(inout) :: failed
      type(response_super_index) :: item
      integer :: target_flat
      complex(rp) :: expected
      real(rp) :: error_value, metric

      item = response_super_index(1, 0, 0, 2, 1)
      call response_flatten_superindex(item, space%nsite, space%response_lmax, space%npoint, &
                                       space%nchannel, target_flat)
      metric = space%metric_weights(target_flat)
      expected = cmplx(2.0_rp/sum(left%k_weights), 0.0_rp, rp)*sum(left%k_weights)/ &
                 cmplx(frequencies(1) - 0.70_rp, eta, rp)*(1.0_rp/(4.0_rp*response_angular_pi))*metric
      error_value = abs(actual%susceptibility(target_flat, target_flat, 1) - expected)
      write (*, '(a,es12.4)') 'Analytic two-level finite-model error = ', error_value
      if (error_value > tolerance) failed = .true.
   end subroutine analytic_two_level_check

   subroutine compare_response(actual, expected, tolerance_in, label, failed)
      complex(rp), intent(in) :: actual(:, :, :), expected(:, :, :)
      real(rp), intent(in) :: tolerance_in
      character(len=*), intent(in) :: label
      logical, intent(inout) :: failed
      real(rp) :: error_value

      error_value = maxval(abs(actual - expected))
      write (*, '(a,es12.4)') trim(label)//' error = ', error_value
      if (error_value > tolerance_in) failed = .true.
   end subroutine compare_response

   subroutine spin_reverse_coefficients(original, reversed)
      complex(rp), intent(in) :: original(:, :, :)
      complex(rp), intent(out) :: reversed(:, :, :)
      integer :: site, band, ik, offset

      reversed = cmplx(0.0_rp, 0.0_rp, rp)
      do ik = 1, size(original, 3)
         do band = 1, size(original, 2)
            do site = 1, nsite
               offset = (site - 1)*2*norb
               reversed(offset + 1:offset + norb, band, ik) = original(offset + norb + 1:offset + 2*norb, band, ik)
               reversed(offset + norb + 1:offset + 2*norb, band, ik) = -original(offset + 1:offset + norb, band, ik)
            end do
         end do
      end do
   end subroutine spin_reverse_coefficients

   subroutine rotate_degenerate_subspaces(original, rotated)
      complex(rp), intent(in) :: original(:, :, :)
      complex(rp), intent(out) :: rotated(:, :, :)
      real(rp), parameter :: theta = 0.37_rp
      real(rp) :: c, s
      integer :: ik

      c = cos(theta)
      s = sin(theta)
      rotated = original
      do ik = 1, size(original, 3)
         rotated(:, 1, ik) = c*original(:, 1, ik) + s*original(:, 2, ik)
         rotated(:, 2, ik) = -s*original(:, 1, ik) + c*original(:, 2, ik)
         rotated(:, 3, ik) = c*original(:, 3, ik) + s*original(:, 4, ik)
         rotated(:, 4, ik) = -s*original(:, 3, ik) + c*original(:, 4, ik)
      end do
   end subroutine rotate_degenerate_subspaces

   subroutine fill_diagnostic_vectors(field, magnetization)
      complex(rp), intent(out) :: field(:), magnetization(:)
      integer :: i

      do i = 1, size(field)
         field(i) = cmplx(0.11_rp + 0.002_rp*real(mod(i, 5), rp), &
                          -0.03_rp + 0.001_rp*real(mod(i, 3), rp), rp)
         magnetization(i) = cmplx(0.07_rp - 0.001_rp*real(mod(i, 4), rp), &
                                  0.02_rp + 0.0015_rp*real(mod(i, 6), rp), rp)
      end do
   end subroutine fill_diagnostic_vectors

end program test_lr_ks_susceptibility
