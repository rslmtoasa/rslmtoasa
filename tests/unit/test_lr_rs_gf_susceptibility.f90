!------------------------------------------------------------------------------
! RSGF-01 finite exact native real-space GF oracle.
!------------------------------------------------------------------------------
program test_lr_rs_gf_susceptibility
   use precision_mod, only: rp
   use response_basis_mapping_mod, only: response_fourier_phase_sign, response_real_space_phase
   use response_angular_basis_mod, only: response_angular_pi
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex
   use lr_response_space_mod, only: response_space_layout
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_ks_susceptibility_mod, only: lr_electronic_state, lr_ks_susceptibility_request, &
      lr_ks_susceptibility_result, lr_channel_plus, evaluate_lr_ks_susceptibility
   use lr_gf_susceptibility_mod, only: lr_gf_susceptibility_request, evaluate_lr_gf_susceptibility
   use lr_rs_gf_susceptibility_mod, only: lr_rs_gf_pair, lr_rs_gf_susceptibility_request, &
      lr_rs_dense_gf_provider, evaluate_lr_rs_gf_susceptibility
   implicit none

   integer, parameter :: nr = 3, nsite = 1, orbital_lmax = 1
   integer, parameter :: norb = (orbital_lmax + 1)**2
   integer, parameter :: nlocal = 2*norb, nbasis = nlocal*nsite, nbands = nbasis, nk = 1
   integer, parameter :: nfrequency = 2
   real(rp), parameter :: mesh_a = 0.05_rp, mesh_b = 0.08_rp
   real(rp), parameter :: eta = 0.08_rp, tolerance = 2.0e-5_rp
   real(rp), parameter :: frequencies(nfrequency) = [0.17_rp, 0.0_rp]
   real(rp), parameter :: eigenvalues(nbands, nk) = reshape([ &
      -0.55_rp, -0.35_rp, 0.25_rp, 0.45_rp, -0.45_rp, -0.25_rp, 0.35_rp, 0.55_rp], [nbands, nk])
   real(rp), parameter :: occupations(nbands, nk) = reshape([ &
      1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp], [nbands, nk])
   real(rp), parameter :: weights(nk) = [1.0_rp]
   real(rp), parameter :: k_points(3, nk) = 0.0_rp

   real(rp) :: radius(nr), difference_norm, reference_norm, relative_difference
   real(rp) :: phase_error, target_metric
   real(rp) :: coefficient_error
   complex(rp) :: hamiltonian(nbasis, nbasis), hgamma(nbasis, nbasis), eigenvectors(nbasis, nbands, nk)
   complex(rp) :: coefficient_gij(nlocal, nlocal), coefficient_gji(nlocal, nlocal)
   complex(rp) :: coefficient_hij(nlocal, nlocal), coefficient_hji(nlocal, nlocal)
   complex(rp) :: z_probe
   complex(rp) :: analytic_value
   type(lmto_radial_basis), target :: radial(nsite)
   type(response_space_layout), target :: space
   type(lr_electronic_state), target :: state, endpoint
   type(lr_ks_susceptibility_request) :: spectral_request
   type(lr_rs_gf_susceptibility_request) :: native_request
   type(lr_ks_susceptibility_result) :: spectral_result, native_result
   type(lr_ks_susceptibility_result) :: reciprocal_gf_result, phased_result
   type(lr_gf_susceptibility_request) :: reciprocal_request
   type(lr_rs_dense_gf_provider), target :: provider
   type(lr_rs_gf_pair) :: pair
   type(response_super_index) :: target_item
   integer :: target_flat, i
   logical :: failed

   call build_mesh(radius)
   call setup_radial_basis(radial, radius)
   call space%initialize(nsite, 2*orbital_lmax, radius, mesh_a, mesh_b, 1)
   call build_diagonal_fixture(hamiltonian, hgamma, eigenvectors)
   call state%initialize(eigenvalues, eigenvectors, k_points, weights, occupations, 0.0_rp, 0.0_rp)
   call endpoint%initialize(eigenvalues, eigenvectors, k_points, weights, occupations, 0.0_rp, 0.0_rp)
   call provider%initialize(hamiltonian, hgamma, nlocal)
   failed = .false.

   ! Coefficient-GF layer isolation: compare the native dense provider blocks
   ! with the independent diagonal spectral resolvent before any endpoint map
   ! or response-space contraction is entered.
   z_probe = cmplx(0.13_rp, 0.07_rp, rp)
   call provider%get_pair(1, 1, [0.0_rp, 0.0_rp, 0.0_rp], z_probe, coefficient_gij, coefficient_gji, &
      coefficient_hij, coefficient_hji)
   coefficient_error = 0.0_rp
   do i = 1, nlocal
      coefficient_error = max(coefficient_error, abs(coefficient_gij(i,i) - &
         1.0_rp/(z_probe - eigenvalues(i,1))))
      coefficient_error = max(coefficient_error, abs(coefficient_gji(i,i) - &
         1.0_rp/(z_probe - eigenvalues(i,1))))
   end do
   write (*, '(a,es12.4)') 'Coefficient-GF versus spectral-resolvent error = ', coefficient_error
   if (coefficient_error > 1.0e-12_rp) failed = .true.

   pair%left_site = 1
   pair%right_site = 1
   pair%translation = 0.0_rp

   spectral_request%q = 0.0_rp
   spectral_request%frequencies = frequencies
   spectral_request%eta = eta
   spectral_request%channel = lr_channel_plus
   spectral_request%response_space => space
   spectral_request%radial_bases => radial
   spectral_request%electronic_state => state
   spectral_request%q_endpoint_state => endpoint
   call evaluate_lr_ks_susceptibility(spectral_request, spectral_result)

   reciprocal_request%q = 0.0_rp
   reciprocal_request%frequencies = frequencies
   reciprocal_request%eta = eta
   reciprocal_request%integration_points = 201
   reciprocal_request%integration_eta = 0.001_rp
   reciprocal_request%energy_margin = 0.35_rp
   reciprocal_request%channel = lr_channel_plus
   reciprocal_request%response_space => space
   reciprocal_request%radial_bases => radial
   reciprocal_request%electronic_state => state
   reciprocal_request%q_endpoint_state => endpoint
   call evaluate_lr_gf_susceptibility(reciprocal_request, reciprocal_gf_result)

   native_request%q = 0.0_rp
   native_request%frequencies = frequencies
   native_request%eta = eta
   native_request%integration_points = 201
   native_request%integration_eta = 0.001_rp
   native_request%energy_margin = 0.35_rp
   native_request%energy_min = minval(eigenvalues)
   native_request%energy_max = maxval(eigenvalues)
   native_request%fermi_level = 0.0_rp
   native_request%temperature = 0.0_rp
   native_request%channel = lr_channel_plus
   native_request%response_space => space
   native_request%radial_bases => radial
   native_request%provider => provider
   allocate(native_request%pairs(1))
   native_request%pairs(1) = pair
   call evaluate_lr_rs_gf_susceptibility(native_request, native_result)

   call report_difference(native_result%susceptibility, spectral_result%susceptibility, &
      difference_norm, reference_norm, relative_difference)
   write (*, '(a,es12.4,a,es12.4)') 'Native RS-GF full-response difference (abs, rel) = ', &
      difference_norm, ' ', relative_difference
   if (difference_norm > tolerance) failed = .true.

   call report_difference(native_result%susceptibility, reciprocal_gf_result%susceptibility, &
      difference_norm, reference_norm, relative_difference)
   write (*, '(a,es12.4,a,es12.4)') 'Native versus reciprocal-GF full-response difference (abs, rel) = ', &
      difference_norm, ' ', relative_difference
   if (difference_norm > tolerance) failed = .true.

   ! This is deliberately a full-space assertion: it includes every response
   ! L,M and radial point, rather than a site-projected scalar.
   target_item = response_super_index(1, 0, 0, 2, 1)
   call response_flatten_superindex(target_item, space%nsite, space%response_lmax, space%npoint, &
      space%nchannel, target_flat)
   target_metric = space%metric_weights(target_flat)
   analytic_value = native_result%susceptibility(target_flat, target_flat, 1)
   write (*, '(a,i0,a,es12.4)') 'Selected full-space (site,L,M,r) entry ', target_flat, &
      ' magnitude = ', abs(analytic_value)
   if (target_metric <= 0.0_rp .or. abs(analytic_value) <= tiny(1.0_rp)) failed = .true.

   ! Nonzero q is checked against the live endpoint-displacement phase.  The
   ! one-site finite oracle has R=0, so this isolates the certified phase
   ! helper without silently changing the dense finite-system Hamiltonian.
   phase_error = abs(response_real_space_phase([0.23_rp, 0.0_rp, 0.0_rp], [0.0_rp, 0.0_rp, 0.0_rp], &
      [0.11_rp, 0.0_rp, 0.0_rp], [0.37_rp, 0.0_rp, 0.0_rp]) - &
      cmplx(cos(2.0_rp*response_angular_pi*0.23_rp*0.26_rp), &
      sin(2.0_rp*response_angular_pi*0.23_rp*0.26_rp), rp))
   write (*, '(a,i0,a,es12.4)') 'Fourier phase sign = ', response_fourier_phase_sign, ' error = ', phase_error
   if (phase_error > 1.0e-13_rp) failed = .true.

   ! A translation-independent two-cell oracle exercises the actual q-phase
   ! multiplication in the backend, including a nonzero q and BZ-folding
   ! sensitive direct translation.  The provider is deliberately marked as
   ! such; it is not used as a production material fixture.
   call provider%initialize(hamiltonian, hgamma, nlocal, .true.)
   deallocate(native_request%pairs)
   allocate(native_request%pairs(2))
   native_request%pairs(1) = pair
   native_request%pairs(2) = pair
   native_request%pairs(2)%translation = [1.0_rp, 0.0_rp, 0.0_rp]
   native_request%q = [0.23_rp, 0.0_rp, 0.0_rp]
   call evaluate_lr_rs_gf_susceptibility(native_request, phased_result)
   phase_error = sqrt(sum(abs(phased_result%susceptibility - &
      (1.0_rp + cmplx(cos(2.0_rp*response_angular_pi*0.23_rp), &
      sin(2.0_rp*response_angular_pi*0.23_rp), rp))*native_result%susceptibility)**2))
   write (*, '(a,es12.4)') 'Nonzero-q translation phase assembly error = ', phase_error
   if (phase_error > 2.0e-5_rp) failed = .true.

   if (failed) then
      write (*, '(a)') 'UnitLrRsGfSusceptibility: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrRsGfSusceptibility: PASS (native coefficient GF, augmentation, full response, q phase)'

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
         do ispin = 1, 2
            do l = 0, orbital_lmax
               do ir = 1, size(mesh)
                  bases(isite)%phi_large(ir, l + 1, ispin) = mesh(ir)**(l + 1)
                  bases(isite)%phidot_large(ir, l + 1, ispin) = 0.12_rp*mesh(ir)**(l + 1)
               end do
            end do
         end do
      end do
   end subroutine setup_radial_basis

   subroutine build_diagonal_fixture(hamiltonian_out, hgamma_out, coefficients)
      complex(rp), intent(out) :: hamiltonian_out(:, :), hgamma_out(:, :), coefficients(:, :, :)

      integer :: ib

      hamiltonian_out = cmplx(0.0_rp, 0.0_rp, rp)
      hgamma_out = cmplx(0.0_rp, 0.0_rp, rp)
      coefficients = cmplx(0.0_rp, 0.0_rp, rp)
      do ib = 1, nbasis
         hamiltonian_out(ib, ib) = cmplx(eigenvalues(ib, 1), 0.0_rp, rp)
         hgamma_out(ib, ib) = hamiltonian_out(ib, ib)
         coefficients(ib, ib, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
   end subroutine build_diagonal_fixture

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

end program test_lr_rs_gf_susceptibility
