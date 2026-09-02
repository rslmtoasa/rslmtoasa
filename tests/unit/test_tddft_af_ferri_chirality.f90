!------------------------------------------------------------------------------
! TDDFT-12 -- signed two-sublattice AF/ferri circular-response regression.
!------------------------------------------------------------------------------
program test_tddft_af_ferri_chirality
   use precision_mod, only: rp
   use response_components_mod, only: RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, build_chi_ks_from_eigenpairs
   use tddft_chi0_green_mod, only: green_chi0_options, eigenpair_green_function_provider, &
      build_chi_ks_from_green_functions
   use tddft_goldstone_mod, only: tddft_goldstone_options, tddft_goldstone_result, evaluate_goldstone
   use xc_response_kernel_mod, only: xc_response_kernel_provider
   implicit none

   logical :: failed

   failed = .false.
   call test_ordered_circular_channels_and_backend_equivalence()
   call test_signed_multisublattice_goldstone_vector()
   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_ordered_circular_channels_and_backend_equivalence()
      real(rp) :: weights(1), eigenvalues(4, 1), omega(2)
      complex(rp) :: eigenvectors(4, 4, 1)
      type(response_channel) :: plus(2), minus(2)
      type(tddft_chi0_options) :: eigen_options
      type(green_chi0_options) :: green_options
      type(tddft_chi0_result) :: plus_result, minus_result, plus_green, minus_green
      type(eigenpair_green_function_provider) :: source

      ! Site 1 is spin-up polarized (occupied up, empty down); site 2 is
      ! spin-down polarized (occupied down, empty up).  The two positive
      ! spin-flip poles therefore belong to opposite circular channels.
      weights = 1.0_rp
      eigenvalues(:, 1) = [-0.20_rp, 0.20_rp, 0.30_rp, -0.30_rp]
      eigenvectors = cmplx(0.0_rp, 0.0_rp, rp)
      eigenvectors(1, 1, 1) = 1.0_rp
      eigenvectors(2, 2, 1) = 1.0_rp
      eigenvectors(3, 3, 1) = 1.0_rp
      eigenvectors(4, 4, 1) = 1.0_rp
      omega = [0.40_rp, 0.60_rp]
      plus = [response_channel(1, RESPONSE_PLUS), response_channel(2, RESPONSE_PLUS)]
      minus = [response_channel(1, RESPONSE_MINUS), response_channel(2, RESPONSE_MINUS)]

      eigen_options%eta = 0.01_rp
      eigen_options%fermi_level = 0.0_rp
      eigen_options%electronic_temperature = 0.0_rp
      eigen_options%k_mesh_shape = [1, 1, 1]
      eigen_options%circular_channel = 'plus_minus'
      call build_chi_ks_from_eigenpairs(weights, eigenvalues, eigenvectors, eigenvalues, eigenvectors, [1, 1], plus, minus, &
         omega, eigen_options, plus_result)
      eigen_options%circular_channel = 'minus_plus'
      call build_chi_ks_from_eigenpairs(weights, eigenvalues, eigenvectors, eigenvalues, eigenvectors, [1, 1], minus, plus, &
         omega, eigen_options, minus_result)

      call assert_true('plus-minus metadata is explicit', trim(plus_result%metadata%circular_channel) == 'plus_minus')
      call assert_true('minus-plus metadata is explicit', trim(minus_result%metadata%circular_channel) == 'minus_plus')
      call assert_true('AF site 1 has the positive plus-minus pole', &
         -aimag(plus_result%chi(1, 1, 1)) > 10.0_rp*(-aimag(plus_result%chi(2, 2, 1))))
      call assert_true('AF site 2 has the positive minus-plus pole', &
         -aimag(minus_result%chi(2, 2, 2)) > 10.0_rp*(-aimag(minus_result%chi(1, 1, 2))))
      call assert_true('opposite circular channels are not collapsed', &
         abs(plus_result%chi(1, 1, 1)-minus_result%chi(1, 1, 1)) > 1.0e-3_rp .and. &
         abs(plus_result%chi(2, 2, 2)-minus_result%chi(2, 2, 2)) > 1.0e-3_rp)

      ! The K-space GF bubble consumes the same ordered vertices and must
      ! reproduce both channels, not just the historical plus-minus route.
      green_options%eta = eigen_options%eta
      green_options%green_eta = 0.5_rp*eigen_options%eta
      green_options%fermi_level = eigen_options%fermi_level
      green_options%electronic_temperature = eigen_options%electronic_temperature
      green_options%energy_min = -1.0_rp
      green_options%energy_max = 1.0_rp
      green_options%energy_points = 4001
      green_options%k_mesh_shape = [1, 1, 1]
      green_options%circular_channel = 'plus_minus'
      call source%initialize(eigenvalues, eigenvectors, eigenvalues, eigenvectors)
      call build_chi_ks_from_green_functions(source, weights, [1, 1], plus, minus, omega, green_options, plus_green)
      green_options%circular_channel = 'minus_plus'
      call build_chi_ks_from_green_functions(source, weights, [1, 1], minus, plus, omega, green_options, minus_green)
      call assert_close_3d('K-space GF plus-minus backend equivalence', plus_green%chi, plus_result%chi, 0.25_rp)
      call assert_close_3d('K-space GF minus-plus backend equivalence', minus_green%chi, minus_result%chi, 0.25_rp)
      call assert_true('GF plus-minus metadata is explicit', trim(plus_green%metadata%circular_channel) == 'plus_minus')
      call assert_true('GF minus-plus metadata is explicit', trim(minus_green%metadata%circular_channel) == 'minus_plus')
   end subroutine test_ordered_circular_channels_and_backend_equivalence

   subroutine test_signed_multisublattice_goldstone_vector()
      type(xc_response_kernel_provider) :: provider
      type(tddft_goldstone_options) :: options
      type(tddft_goldstone_result) :: result
      complex(rp) :: chi(2, 2)

      call provider%initialize(2, 'TDDFT-12 unit AF/ferri fixture')
      ! Keep the old field magnitude positive while supplying a deliberately
      ! reversed signed sublattice population.  A hidden abs(moment) would
      ! fail the second Ward component below.
      call provider%set_site_spin_population(1, 2.0_rp)
      call provider%set_site_spin_population(2, 1.0_rp)
      call provider%set_site_signed_spin_population(1, 2.0_rp)
      call provider%set_site_signed_spin_population(2, -1.0_rp)
      provider%site(1)%k_perp_circular = 0.5_rp
      provider%site(2)%k_perp_circular = 0.25_rp
      provider%site(1)%has_k_perp_circular = .true.
      provider%site(2)%has_k_perp_circular = .true.
      chi = cmplx(0.0_rp, 0.0_rp, rp)
      chi(1, 1) = 2.0_rp
      chi(2, 2) = 4.0_rp
      options%goldstone_mode = 'diagnose'
      options%circular_channel = 'minus_plus'
      call evaluate_goldstone(chi, provider, options, result)
      call assert_true('signed AF/ferri Goldstone Ward record is available', result%raw%ward%available)
      call assert_close_vector('signed AF/ferri Goldstone vector', result%raw%ward%magnetization, &
         [cmplx(2.0_rp, 0.0_rp, rp), cmplx(-1.0_rp, 0.0_rp, rp)], 1.0e-12_rp)
      call assert_true('signed AF/ferri Goldstone Ward identity passes', result%raw%ward%ward_residual < 1.0e-12_rp)
      call assert_true('Goldstone report retains circular channel', trim(result%circular_channel) == 'minus_plus')
   end subroutine test_signed_multisublattice_goldstone_vector

   subroutine assert_close_3d(label, actual, expected, tolerance)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual(:, :, :), expected(:, :, :)
      real(rp), intent(in) :: tolerance

      if (maxval(abs(actual-expected)) > tolerance) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, maxval(abs(actual-expected))
         failed = .true.
      end if
   end subroutine assert_close_3d

   subroutine assert_close_vector(label, actual, expected, tolerance)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual(:), expected(:)
      real(rp), intent(in) :: tolerance

      if (maxval(abs(actual-expected)) > tolerance) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, maxval(abs(actual-expected))
         failed = .true.
      end if
   end subroutine assert_close_vector

   subroutine assert_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition

      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine assert_true

end program test_tddft_af_ferri_chirality
