!------------------------------------------------------------------------------
! TDDFT-08 -- mixed Lounis contour versus the direct real-axis reference.
!------------------------------------------------------------------------------
program test_tddft_contour
   use precision_mod, only: rp
   use response_components_mod, only: RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel
   use tddft_chi0_mod, only: tddft_chi0_result
   use tddft_chi0_green_mod, only: green_chi0_options, eigenpair_green_function_provider, &
      build_chi_ks_from_green_functions
   implicit none

   logical :: failed

   failed = .false.
   call test_zero_temperature_mixed_contour()
   call test_finite_temperature_mixed_contour()
   if (failed) error stop 1
   write (*, '(a)') 'TDDFT mixed contour tests passed.'

contains

   subroutine test_zero_temperature_mixed_contour()
      type(eigenpair_green_function_provider) :: source
      type(green_chi0_options) :: direct_options, mixed_options, coarse_options
      type(tddft_chi0_result) :: direct, mixed, coarse
      type(response_channel) :: left(1), right(1)
      real(rp) :: eigenvalues(2, 1), omega(2), weights(1), t_start, t_stop, direct_seconds, mixed_seconds
      complex(rp) :: eigenvectors(2, 2, 1)

      call make_fixture(eigenvalues, eigenvectors, source, left, right, weights)
      omega = [0.20_rp, 0.40_rp]
      direct_options%eta = 0.04_rp; direct_options%green_eta = 0.02_rp
      direct_options%fermi_level = 0.0_rp; direct_options%energy_min = -2.0_rp; direct_options%energy_max = 2.0_rp
      direct_options%energy_points = 16001
      mixed_options = direct_options; mixed_options%energy_integration = 'mixed_contour'
      mixed_options%contour_points = 64; mixed_options%near_fermi_points = 256
      coarse_options = mixed_options; coarse_options%contour_points = 24

      call cpu_time(t_start)
      call build_chi_ks_from_green_functions(source, weights, [1], left, right, omega, direct_options, direct)
      call cpu_time(t_stop); direct_seconds = t_stop-t_start
      call cpu_time(t_start)
      call build_chi_ks_from_green_functions(source, weights, [1], left, right, omega, mixed_options, mixed)
      call cpu_time(t_stop); mixed_seconds = t_stop-t_start
      call build_chi_ks_from_green_functions(source, weights, [1], left, right, omega, coarse_options, coarse)

      call check_close('zero-temperature mixed contour', mixed%chi, direct%chi, 2.0e-5_rp)
      call check_close('coarse contour remains controlled', coarse%chi, direct%chi, 2.0e-3_rp)
      call check_true('mixed path is explicitly reported', mixed%metadata%contour_deformation .and. &
         index(trim(mixed%metadata%energy_integration), 'mixed contour') > 0)
      call check_true('contour moves GF arguments away from real poles', &
         mixed%metadata%contour_max_imaginary_energy > mixed%metadata%green_eta)
      call check_true('horizontal contour subdivision is reported', mixed%metadata%contour_subdivisions == 8)
      call check_true('mixed path uses fewer GF evaluations', &
         mixed%metadata%green_function_evaluations < direct%metadata%green_function_evaluations)
      call check_true('positive-frequency circular response is causal', mixed%im_chi(1, 1, 2) <= 1.0e-10_rp)
      call check_true('positive-frequency circular spectral weight is nonnegative', &
         mixed%site_diagonal_spectrum(1, 2) >= -1.0e-10_rp)
      write (*, '(a,4(1x,es16.8),2(1x,i0))') 'TDDFT08_T0 direct_s mixed_s direct_gf mixed_gf', direct_seconds, mixed_seconds, &
         real(direct%metadata%green_function_evaluations, rp), real(mixed%metadata%green_function_evaluations, rp), &
         direct%metadata%green_function_evaluations, mixed%metadata%green_function_evaluations
      write (*, '(a,3(1x,i0))') 'TDDFT08_T0 gf_calls direct mixed coarse', direct%metadata%green_function_evaluations, &
         mixed%metadata%green_function_evaluations, coarse%metadata%green_function_evaluations
      write (*, '(a,3(1x,es16.8))') 'TDDFT08_T0 maxerr_coarse fine contour_height max_imag', &
         maxval(abs(coarse%chi-direct%chi)), maxval(abs(mixed%chi-direct%chi)), mixed%metadata%contour_height, &
         mixed%metadata%contour_max_imaginary_energy
   end subroutine test_zero_temperature_mixed_contour

   subroutine test_finite_temperature_mixed_contour()
      type(eigenpair_green_function_provider) :: source
      type(green_chi0_options) :: direct_options, mixed_options
      type(tddft_chi0_result) :: direct, mixed
      type(response_channel) :: left(1), right(1)
      real(rp) :: eigenvalues(2, 1), omega(1), weights(1)
      complex(rp) :: eigenvectors(2, 2, 1)

      call make_fixture(eigenvalues, eigenvectors, source, left, right, weights)
      omega = [0.40_rp]
      direct_options%eta = 0.04_rp; direct_options%green_eta = 0.02_rp
      direct_options%fermi_level = 0.0_rp; direct_options%electronic_temperature = 2000.0_rp
      direct_options%energy_min = -2.0_rp; direct_options%energy_max = 2.0_rp; direct_options%energy_points = 20001
      mixed_options = direct_options; mixed_options%energy_integration = 'mixed_contour'
      mixed_options%contour_points = 64; mixed_options%near_fermi_points = 384
      call build_chi_ks_from_green_functions(source, weights, [1], left, right, omega, direct_options, direct)
      call build_chi_ks_from_green_functions(source, weights, [1], left, right, omega, mixed_options, mixed)
      call check_close('finite-temperature mixed contour', mixed%chi, direct%chi, 3.0e-4_rp)
      call check_true('finite-temperature contour remains below first Fermi pole', &
         mixed%metadata%contour_height < 0.95_rp*3.14159265358979323846_rp*direct_options%electronic_temperature*6.3336814e-6_rp)
      write (*, '(a,3(1x,es16.8),2(1x,i0))') 'TDDFT08_Tfinite err height max_imag direct_gf mixed_gf', &
         maxval(abs(mixed%chi-direct%chi))/max(1.0_rp, maxval(abs(direct%chi))), mixed%metadata%contour_height, &
         mixed%metadata%contour_max_imaginary_energy, direct%metadata%green_function_evaluations, &
         mixed%metadata%green_function_evaluations
   end subroutine test_finite_temperature_mixed_contour

   subroutine make_fixture(eigenvalues, eigenvectors, source, left, right, weights)
      real(rp), intent(out) :: eigenvalues(:, :), weights(:)
      complex(rp), intent(out) :: eigenvectors(:, :, :)
      type(eigenpair_green_function_provider), intent(out) :: source
      type(response_channel), intent(out) :: left(:), right(:)

      eigenvalues(:, 1) = [-0.20_rp, 0.20_rp]
      eigenvectors = cmplx(0.0_rp, 0.0_rp, rp)
      eigenvectors(1, 1, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      eigenvectors(2, 2, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      weights = 1.0_rp
      left(1) = response_channel(1, RESPONSE_PLUS)
      right(1) = response_channel(1, RESPONSE_MINUS)
      call source%initialize(eigenvalues, eigenvectors, eigenvalues, eigenvectors)
   end subroutine make_fixture

   subroutine check_close(label, actual, expected, tolerance)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual(:, :, :), expected(:, :, :)
      real(rp), intent(in) :: tolerance
      real(rp) :: error, scale

      scale = max(1.0_rp, maxval(abs(expected))); error = maxval(abs(actual-expected))/scale
      if (error > tolerance) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, error
         failed = .true.
      end if
   end subroutine check_close

   subroutine check_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine check_true

end program test_tddft_contour
