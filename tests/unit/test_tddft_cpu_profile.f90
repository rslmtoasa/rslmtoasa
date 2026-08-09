! TDDFT-11 CPU profiling fixture.  The dimensions use the 18-spinor
! (s,p,d) one-atom basis of the bcc Fe and fcc Ni examples; the deterministic
! synthetic spectra isolate response-kernel cost from an SCF calculation.
program test_tddft_cpu_profile
   use precision_mod, only: rp
   use response_components_mod, only: RESPONSE_CHARGE, RESPONSE_MZ, RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, build_chi_ks_from_eigenpairs
   use tddft_chi0_green_mod, only: green_chi0_options, eigenpair_green_function_provider, &
      build_chi_ks_from_green_functions
   use tddft_dyson_mod, only: tddft_dyson_options, tddft_dyson_result, enhance_tddft_susceptibility
   use tddft_modes_mod, only: tddft_mode_options, tddft_mode_result, analyze_tddft_modes
   implicit none

   call profile_fixture('bccFe', 16, 96)
   call profile_fixture('fccNi', 32, 192)
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine profile_fixture(label, nk, nw)
      character(len=*), intent(in) :: label
      integer, intent(in) :: nk, nw
      real(rp), allocatable :: weights(:), eigenvalues(:, :), omega(:), trace_loss(:, :)
      complex(rp), allocatable :: eigenvectors(:, :, :), kernel(:, :), xi(:, :, :, :)
      type(response_channel) :: left(2), right(2)
      type(tddft_chi0_options) :: chi_options
      type(green_chi0_options) :: green_options
      type(tddft_dyson_options) :: dyson_options
      type(tddft_mode_options) :: mode_options
      type(tddft_chi0_result) :: chi, green_chi
      type(eigenpair_green_function_provider) :: green_source
      type(tddft_dyson_result) :: dyson
      type(tddft_mode_result) :: modes
      integer :: ik, ib, ic
      real(rp) :: norm

      allocate(weights(nk), eigenvalues(18, nk), eigenvectors(18, 18, nk), omega(nw), kernel(2, 2))
      weights = 1.0_rp
      do ik = 1, nk
         do ib = 1, 18
            eigenvalues(ib, ik) = -0.42_rp + 0.048_rp*real(ib-1, rp) + 0.003_rp*real(ik-1, rp)
            do ic = 1, 18
               eigenvectors(ic, ib, ik) = cmplx(sin(real(3*ic+ib+ik, rp)), cos(real(ic-2*ib+ik, rp)), rp)
            end do
            norm = sqrt(sum(abs(eigenvectors(:, ib, ik))**2))
            eigenvectors(:, ib, ik) = eigenvectors(:, ib, ik)/norm
         end do
      end do
      do ik = 1, nw
         omega(ik) = 0.005_rp + 0.20_rp*real(ik-1, rp)/real(nw-1, rp)
      end do
      left = [response_channel(1, RESPONSE_PLUS), response_channel(1, RESPONSE_MZ)]
      right = [response_channel(1, RESPONSE_MINUS), response_channel(1, RESPONSE_CHARGE)]
      chi_options%eta = 0.006_rp
      chi_options%fermi_level = 0.0_rp
      chi_options%electronic_temperature = 300.0_rp
      chi_options%transition_batch_size = 128
      call build_chi_ks_from_eigenpairs(weights, eigenvalues, eigenvectors, eigenvalues, eigenvectors, [9], &
         left, right, omega, chi_options, chi)
      call green_source%initialize(eigenvalues, eigenvectors, eigenvalues, eigenvectors)
      green_options%eta = chi_options%eta
      green_options%fermi_level = chi_options%fermi_level
      green_options%electronic_temperature = chi_options%electronic_temperature
      green_options%energy_min = -0.50_rp
      green_options%energy_max = 0.50_rp
      green_options%energy_points = 13
      call build_chi_ks_from_green_functions(green_source, weights(:2), [9], left, right, omega(:8), green_options, green_chi)
      kernel = cmplx(0.0_rp, 0.0_rp, rp)
      kernel(1, 1) = cmplx(0.4_rp, 0.0_rp, rp)
      kernel(2, 2) = cmplx(0.2_rp, 0.0_rp, rp)
      call enhance_tddft_susceptibility(chi%chi, kernel, chi_options%eta, dyson_options, dyson)
      allocate(xi(2, 2, nw, 1), trace_loss(nw, 1))
      xi(:, :, :, 1) = dyson%xi; trace_loss(:, 1) = dyson%trace_spectral_weight
      call analyze_tddft_modes(omega, xi, trace_loss, chi_options%eta, mode_options, modes)
      write (*, '(a,1x,a,8(1x,es12.4))') 'PROFILE', trim(label), chi%metadata%vertex_cpu_seconds, &
         chi%metadata%denominator_cpu_seconds, chi%metadata%accumulation_cpu_seconds, &
         chi%metadata%arbitrary_kq_cpu_seconds, dyson%metadata%solve_cpu_seconds, &
         dyson%metadata%diagonalization_cpu_seconds, modes%analysis_cpu_seconds, &
         green_chi%metadata%green_energy_integration_cpu_seconds
   end subroutine profile_fixture

end program test_tddft_cpu_profile
