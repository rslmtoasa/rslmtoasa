!------------------------------------------------------------------------------
! LR-GF-02 reciprocal-GF susceptibility cross-check.
!
! The production-style fixture is intentionally small, but it retains the
! complete direct radial response space.  The reference is the LR-06 spectral
! route; the GF route is evaluated through its independent real-axis bubble.
!------------------------------------------------------------------------------
program test_lr_gf_susceptibility
   use precision_mod, only: rp
   use response_angular_basis_mod, only: response_angular_pi
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex
   use lr_response_space_mod, only: response_space_layout
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_ks_susceptibility_mod, only: lr_electronic_state, lr_ks_susceptibility_request, &
      lr_ks_susceptibility_result, lr_channel_plus, evaluate_lr_ks_susceptibility
   use lr_gf_susceptibility_mod, only: lr_gf_susceptibility_request, evaluate_lr_gf_susceptibility
   implicit none

   integer, parameter :: nr = 3, nsite = 1, orbital_lmax = 1
   integer, parameter :: norb = (orbital_lmax + 1)**2, nbasis = 2*nsite*norb, nbands = 4, nk = 1, nfrequency = 2
   real(rp), parameter :: mesh_a = 0.05_rp, mesh_b = 0.08_rp
   real(rp), parameter :: eta = 0.04_rp, tolerance = 5.0e-7_rp
   real(rp), parameter :: frequencies(nfrequency) = [0.17_rp, 0.0_rp]
   real(rp), parameter :: eigenvalues(nbands, nk) = reshape([ &
      -0.40_rp, -0.40_rp, 0.30_rp, 0.30_rp, &
      -0.40_rp, -0.40_rp, 0.30_rp, 0.30_rp], [nbands, nk])
   real(rp), parameter :: occupations(nbands, nk) = reshape([ &
      1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, &
      1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp], [nbands, nk])
   real(rp), parameter :: weights(nk) = [1.0_rp]
   real(rp), parameter :: k_points(3, nk) = reshape([0.0_rp, 0.0_rp, 0.0_rp], [3, nk])

   real(rp) :: radius(nr), difference_norm, reference_norm, relative_difference
   real(rp) :: q_endpoint_points(3, nk)
   complex(rp) :: eigenvectors(nbasis, nbands, nk)
   type(lmto_radial_basis), target :: radial(nsite)
   type(response_space_layout), target :: space
   type(lr_electronic_state), target :: state, endpoint
   type(lr_ks_susceptibility_request) :: spectral_request
   type(lr_gf_susceptibility_request) :: gf_request
   type(lr_ks_susceptibility_result) :: spectral_result, gf_coarse, gf_fine
   type(lr_ks_susceptibility_result) :: spectral_eta_result, gf_eta_result
   type(lr_ks_susceptibility_result) :: spectral_q_result, gf_q_result
   type(response_super_index) :: target_item
   integer :: target_flat
   real(rp) :: analytic_error, target_metric
   complex(rp) :: analytic_value
   logical :: failed

   call build_mesh(radius)
   call setup_radial_basis(radial, radius)
   call space%initialize(nsite, 2*orbital_lmax, radius, mesh_a, mesh_b, 1)
   call build_coefficients(eigenvectors)
   call state%initialize(eigenvalues, eigenvectors, k_points, weights, occupations, 0.0_rp, 0.0_rp)
   call endpoint%initialize(eigenvalues, eigenvectors, k_points, weights, occupations, 0.0_rp, 0.0_rp)
   failed = .false.

   spectral_request%q = [0.0_rp, 0.0_rp, 0.0_rp]
   spectral_request%frequencies = frequencies
   spectral_request%eta = eta
   spectral_request%channel = lr_channel_plus
   spectral_request%response_space => space
   spectral_request%radial_bases => radial
   spectral_request%electronic_state => state
   spectral_request%q_endpoint_state => endpoint
   call evaluate_lr_ks_susceptibility(spectral_request, spectral_result)

   gf_request%q = spectral_request%q
   gf_request%frequencies = frequencies
   gf_request%eta = eta
   gf_request%channel = lr_channel_plus
   gf_request%integration_points = 401
   gf_request%integration_eta = 0.010_rp
   gf_request%energy_margin = 1.0_rp
   gf_request%response_space => space
   gf_request%radial_bases => radial
   gf_request%electronic_state => state
   gf_request%q_endpoint_state => endpoint
   call evaluate_lr_gf_susceptibility(gf_request, gf_coarse)

   gf_request%integration_points = 2001
   gf_request%integration_eta = 0.001_rp
   call evaluate_lr_gf_susceptibility(gf_request, gf_fine)

   call report_difference(gf_coarse%susceptibility, spectral_result%susceptibility, &
      difference_norm, reference_norm, relative_difference)
   write (*, '(a,es12.4,a,es12.4)') 'GF coarse full-response difference (abs, rel) = ', &
      difference_norm, ' ', relative_difference
   if (difference_norm > tolerance) failed = .true.

   call report_difference(gf_fine%susceptibility, spectral_result%susceptibility, &
      difference_norm, reference_norm, relative_difference)
   write (*, '(a,es12.4,a,es12.4)') 'GF fine full-response difference (abs, rel) = ', &
      difference_norm, ' ', relative_difference
   if (difference_norm > tolerance) failed = .true.

   target_item = response_super_index(1, 0, 0, 2, 1)
   call response_flatten_superindex(target_item, space%nsite, space%response_lmax, space%npoint, &
      space%nchannel, target_flat)
   target_metric = space%metric_weights(target_flat)
   ! Nonzero phidot exercises the GF energy moments rather than only the
   ! energy-independent p=q=0 vertex component.
   analytic_value = cmplx(2.0_rp, 0.0_rp, rp)/cmplx(frequencies(1) - 0.70_rp, eta, rp)* &
      cmplx(((1.0_rp - 0.40_rp*0.20_rp)*(1.0_rp + 0.30_rp*0.20_rp))**2/ &
      (4.0_rp*response_angular_pi)*target_metric, 0.0_rp, rp)
   analytic_error = abs(gf_fine%susceptibility(target_flat, target_flat, 1) - analytic_value)
   write (*, '(a,es12.4)') 'Finite-spectrum analytic GF error = ', analytic_error
   if (analytic_error > tolerance) failed = .true.

   ! The second frequency is the independent q=0, omega=0 static comparison.
   write (*, '(a,es12.4)') 'Static q=0 full-response GF difference = ', &
      maxval(abs(gf_fine%susceptibility(:, :, 2) - spectral_result%susceptibility(:, :, 2)))
   if (maxval(abs(gf_fine%susceptibility(:, :, 2) - spectral_result%susceptibility(:, :, 2))) > tolerance) then
      failed = .true.
   end if

   ! Off-mesh reciprocal endpoint: the GF route must use the exact folded
   ! endpoint ordering rather than a nearest-mesh substitute.
   q_endpoint_points = k_points
   q_endpoint_points(1, 1) = 0.23_rp
   call endpoint%initialize(eigenvalues, eigenvectors, q_endpoint_points, weights, occupations, 0.0_rp, 0.0_rp)
   spectral_request%q = [0.23_rp, 0.0_rp, 0.0_rp]
   spectral_request%eta = eta
   call evaluate_lr_ks_susceptibility(spectral_request, spectral_q_result)
   gf_request%q = spectral_request%q
   gf_request%eta = eta
   gf_request%integration_points = 1001
   gf_request%integration_eta = 0.001_rp
   call evaluate_lr_gf_susceptibility(gf_request, gf_q_result)
   call report_difference(gf_q_result%susceptibility, spectral_q_result%susceptibility, &
      difference_norm, reference_norm, relative_difference)
   write (*, '(a,es12.4,a,es12.4)') 'GF off-mesh-q full-response difference (abs, rel) = ', &
      difference_norm, ' ', relative_difference
   if (difference_norm > tolerance) failed = .true.

   ! Repeat one response point with a different physical eta.  The GF route
   ! must track the requested retarded broadening, not only one fixed width.
   spectral_request%eta = 0.02_rp
   call evaluate_lr_ks_susceptibility(spectral_request, spectral_eta_result)
   gf_request%eta = 0.02_rp
   gf_request%integration_points = 1001
   gf_request%integration_eta = 0.001_rp
   call evaluate_lr_gf_susceptibility(gf_request, gf_eta_result)
   call report_difference(gf_eta_result%susceptibility, spectral_eta_result%susceptibility, &
      difference_norm, reference_norm, relative_difference)
   write (*, '(a,es12.4,a,es12.4)') 'GF alternate eta full-response difference (abs, rel) = ', &
      difference_norm, ' ', relative_difference
   if (difference_norm > tolerance) failed = .true.

   if (failed) then
      write (*, '(a)') 'UnitLrGfSusceptibility: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrGfSusceptibility: PASS (GF bubble, convergence, full response, static limit)'

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
            bases(isite)%phidot_large(:, 1, ispin) = 0.20_rp*mesh
         end do
      end do
   end subroutine setup_radial_basis

   subroutine build_coefficients(coefficients)
      complex(rp), intent(out) :: coefficients(:, :, :)
      integer :: ik

      coefficients = cmplx(0.0_rp, 0.0_rp, rp)
      do ik = 1, nk
         coefficients(1, 1, ik) = cmplx(1.0_rp, 0.0_rp, rp)
         coefficients(1 + norb, 3, ik) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
   end subroutine build_coefficients

   subroutine report_difference(actual, reference, absolute, reference_size, relative)
      complex(rp), intent(in) :: actual(:, :, :), reference(:, :, :)
      real(rp), intent(out) :: absolute, reference_size, relative

      if (any(shape(actual) /= shape(reference))) error stop 'report_difference: response shape mismatch'
      absolute = sqrt(sum(abs(actual - reference)**2))
      reference_size = sqrt(sum(abs(reference)**2))
      if (reference_size > tiny(1.0_rp)) then
         relative = absolute/reference_size
      else
         relative = absolute
      end if
   end subroutine report_difference

end program test_lr_gf_susceptibility
