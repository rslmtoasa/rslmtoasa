program test_tddft_realspace_gf
   use precision_mod, only: rp
   use math_mod, only: pi, i_unit
   use response_components_mod, only: RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel
   use tddft_chi0_mod, only: tddft_chi0_request, tddft_chi0_batch_result
   use tddft_chi0_mod, only: tddft_chi0_result
   use tddft_chi0_green_mod, only: green_chi0_options, eigenpair_green_function_provider, &
      build_chi_ks_from_green_functions
   use tddft_chi0_realspace_mod, only: tddft_realspace_chi0_options, tddft_realspace_chi0_result, &
      tddft_native_realspace_gf_provider, fourier_transform_realspace_chi0, fourier_transform_realspace_green, &
      check_realspace_pair_reversal
   implicit none

   call test_fourier_transforms()
   call test_native_provider_batch_and_cutoff()
   call test_periodic_kspace_agreement()
   call test_direction_reversal()
   write(*, '(a)') 'TDDFT native real-space Green-function tests passed.'

contains

   subroutine test_fourier_transforms()
      complex(rp) :: chi_r(1, 1, 1, 2), g_r(1, 1, 1, 2)
      complex(rp), allocatable :: chi_q(:, :, :, :), g_q(:, :, :, :)
      real(rp) :: r_vectors(3, 2), q_points(3, 2), expected

      chi_r = cmplx(0.0_rp, 0.0_rp, rp); chi_r(1, 1, 1, :) = [cmplx(2.0_rp, 0.0_rp, rp), cmplx(3.0_rp, 0.0_rp, rp)]
      g_r = chi_r
      r_vectors = 0.0_rp; r_vectors(1, 2) = 1.0_rp
      q_points = 0.0_rp; q_points(1, 2) = 0.25_rp
      call fourier_transform_realspace_chi0(chi_r, r_vectors, q_points, 'bulk', [1, 2, 3], chi_q)
      call assert_close(chi_q(1, 1, 1, 1), cmplx(5.0_rp, 0.0_rp, rp), 1.0e-12_rp, 'zero-q susceptibility FT')
      expected = 2.0_rp
      call assert_close(chi_q(1, 1, 1, 2), cmplx(expected, -3.0_rp, rp), 1.0e-12_rp, 'bulk susceptibility phase')
      call fourier_transform_realspace_green(g_r, r_vectors, q_points, g_q)
      call assert_close(g_q(1, 1, 1, 2), cmplx(expected, -3.0_rp, rp), 1.0e-12_rp, 'validation Green phase')
      call fourier_transform_realspace_chi0(chi_r, r_vectors, q_points, 'finite', [1, 2, 3], chi_q)
      call assert_close(chi_q(1, 1, 1, 2), cmplx(5.0_rp, 0.0_rp, rp), 1.0e-12_rp, 'finite direct response')
      r_vectors(3, 2) = 1.0_rp; q_points(3, 2) = 0.25_rp
      call fourier_transform_realspace_chi0(chi_r, r_vectors, q_points, 'film', [1, 2, 0], chi_q)
      call assert_close(chi_q(1, 1, 1, 2), cmplx(2.0_rp, -3.0_rp, rp), 1.0e-12_rp, 'film in-plane FT')
   end subroutine test_fourier_transforms

   subroutine test_native_provider_batch_and_cutoff()
      type(tddft_native_realspace_gf_provider) :: provider
      type(tddft_realspace_chi0_options) :: options
      type(tddft_chi0_request) :: request
      type(tddft_chi0_batch_result) :: batch_result
      type(tddft_chi0_batch_result) :: reversed_batch
      type(tddft_realspace_chi0_result) :: direct
      type(response_channel) :: left(1), right(1)
      real(rp) :: energy(401), r_vectors(3, 2), q_points(3, 3), omega(2)
      integer :: pair_sites(2, 2), counts(2), ie
      complex(rp) :: g_ab(2, 2, 401, 2), g_ba(2, 2, 401, 2)

      do ie = 1, size(energy)
         energy(ie) = -1.0_rp + 0.005_rp*real(ie-1, rp)
         g_ab(:, :, ie, :) = cmplx(0.0_rp, 0.0_rp, rp)
         g_ba(:, :, ie, :) = cmplx(0.0_rp, 0.0_rp, rp)
         g_ab(1, 1, ie, :) = 1.0_rp/cmplx(energy(ie)+0.2_rp, 0.05_rp, rp)
         g_ab(2, 2, ie, :) = 1.0_rp/cmplx(energy(ie)-0.2_rp, 0.05_rp, rp)
         g_ba(:, :, ie, :) = g_ab(:, :, ie, :)
      end do
      r_vectors = 0.0_rp; r_vectors(1, 2) = 1.0_rp
      pair_sites(1, :) = [1, 2]; pair_sites(2, :) = [1, 2]
      counts = [1, 1]
      left(1)%site = 1; left(1)%component = RESPONSE_PLUS
      right(1)%site = 2; right(1)%component = RESPONSE_MINUS
      omega = [0.1_rp, 0.2_rp]
      q_points = 0.0_rp; q_points(1, 2) = 0.25_rp; q_points(1, 3) = 0.5_rp
      options%eta = 0.1_rp; options%green_eta = 0.05_rp; options%fermi_level = 0.0_rp
      options%rmax = huge(1.0_rp); options%tail_tolerance = 1.0e-8_rp
      call provider%initialize(energy, g_ab, g_ba, r_vectors, pair_sites, counts, left, right, options)
      allocate(request%q_points(3, 3), request%omega(2)); request%q_points = q_points; request%omega = omega
      call provider%evaluate_realspace(request, batch_result)
      call assert(provider%build_count == 1, 'all q points reuse one chi0(R,omega) build')
      call assert(size(batch_result%q_response) == 3, 'native q batch shape')
      call assert(abs(batch_result%q_response(1)%chi(1, 1, 1)-batch_result%q_response(2)%chi(1, 1, 1)) > 1.0e-10_rp, &
         'susceptibility FT has q dependence')
      call assert(batch_result%metadata%real_space_points == 2, 'native pair count metadata')
      call assert(.not. batch_result%metadata%real_space_tail_assessed, 'untruncated source reports no tail estimate')
      g_ba = 1.7_rp*g_ab
      call provider%initialize(energy, g_ab, g_ba, r_vectors, pair_sites, counts, left, right, options)
      call provider%evaluate_realspace(request, reversed_batch)
      call assert(abs(batch_result%q_response(1)%chi(1, 1, 1)-reversed_batch%q_response(1)%chi(1, 1, 1)) > 1.0e-10_rp, &
         'native contraction uses the supplied reverse propagator')

      options%rmax = 0.5_rp; options%tail_tolerance = 0.0_rp
      call provider%initialize(energy, g_ab, g_ba, r_vectors, pair_sites, counts, left, right, options)
      call provider%build_realspace(omega, direct)
      call assert(direct%diagnostics%selected_points == 1, 'Rmax selects the near pair')
      call assert(direct%diagnostics%omitted_points == 1, 'Rmax exposes omitted tail points')
      call assert(direct%diagnostics%tail_assessed, 'Rmax enables tail assessment')
      call assert(direct%diagnostics%omitted_tail_norm > 0.0_rp, 'omitted real-space tail has a norm')
   end subroutine test_native_provider_batch_and_cutoff

   subroutine test_periodic_kspace_agreement()
      type(tddft_native_realspace_gf_provider) :: native_provider
      type(eigenpair_green_function_provider) :: kspace_provider
      type(tddft_realspace_chi0_options) :: native_options
      type(green_chi0_options) :: kspace_options
      type(tddft_chi0_request) :: request
      type(tddft_chi0_batch_result) :: native_batch
      type(tddft_chi0_result) :: kspace_result
      type(response_channel) :: left(1), right(1)
      real(rp) :: energy(401), r_vectors(3, 1), q_points(3, 1), omega(1), weights(1)
      real(rp) :: eigenvalues(2, 1)
      complex(rp) :: g_ab(2, 2, 401, 1), g_ba(2, 2, 401, 1), eigenvectors(2, 2, 1)
      integer :: counts(1), pair_sites(1, 2), ie

      do ie = 1, size(energy)
         energy(ie) = -1.0_rp + 0.005_rp*real(ie-1, rp)
         g_ab(:, :, ie, 1) = cmplx(0.0_rp, 0.0_rp, rp)
         g_ab(1, 1, ie, 1) = 1.0_rp/cmplx(energy(ie)+0.2_rp, 0.05_rp, rp)
         g_ab(2, 2, ie, 1) = 1.0_rp/cmplx(energy(ie)-0.2_rp, 0.05_rp, rp)
      end do
      g_ba = g_ab; r_vectors = 0.0_rp; pair_sites(1, :) = [1, 1]
      counts = [1]; left(1)%site = 1; left(1)%component = RESPONSE_PLUS
      right(1)%site = 1; right(1)%component = RESPONSE_MINUS
      omega = [0.1_rp]; q_points = 0.0_rp
      native_options%eta = 0.1_rp; native_options%green_eta = 0.05_rp
      call native_provider%initialize(energy, g_ab, g_ba, r_vectors, pair_sites, counts, left, right, native_options)
      allocate(request%q_points(3, 1), request%omega(1)); request%q_points = q_points; request%omega = omega
      call native_provider%evaluate_realspace(request, native_batch)

      eigenvalues(:, 1) = [-0.2_rp, 0.2_rp]; eigenvectors = cmplx(0.0_rp, 0.0_rp, rp)
      eigenvectors(1, 1, 1) = cmplx(1.0_rp, 0.0_rp, rp); eigenvectors(2, 2, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      weights = [1.0_rp]
      call kspace_provider%initialize(eigenvalues, eigenvectors, eigenvalues, eigenvectors)
      kspace_options%eta = 0.1_rp; kspace_options%green_eta = 0.05_rp; kspace_options%energy_min = -0.9_rp
      kspace_options%energy_max = 0.9_rp; kspace_options%energy_points = 361
      call build_chi_ks_from_green_functions(kspace_provider, weights, counts, left, right, omega, kspace_options, kspace_result)
      call assert_close(native_batch%q_response(1)%chi(1, 1, 1), kspace_result%chi(1, 1, 1), 2.0e-10_rp, &
         'periodic native real-space and K-space GF bubbles agree')
   end subroutine test_periodic_kspace_agreement

   subroutine test_direction_reversal()
      integer :: pair_sites(2, 2)
      real(rp) :: r_vectors(3, 2), residual
      logical :: complete

      pair_sites = reshape([1, 2, 2, 1], [2, 2])
      r_vectors = 0.0_rp; r_vectors(1, 1) = 1.0_rp; r_vectors(1, 2) = -1.0_rp
      call check_realspace_pair_reversal(pair_sites, r_vectors, 1.0e-12_rp, complete, residual)
      call assert(complete, 'real-space pair reversal map')
      call assert(residual < 1.0e-12_rp, 'real-space reversal vector residual')
   end subroutine test_direction_reversal

   subroutine assert(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) error stop 'FAIL: '//trim(message)
   end subroutine assert

   subroutine assert_close(value, expected, tolerance, message)
      complex(rp), intent(in) :: value, expected
      real(rp), intent(in) :: tolerance
      character(len=*), intent(in) :: message
      call assert(abs(value-expected) <= tolerance, message)
   end subroutine assert_close

end program test_tddft_realspace_gf
