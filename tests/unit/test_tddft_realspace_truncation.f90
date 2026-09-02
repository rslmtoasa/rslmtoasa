! TDDFT-R2-04 -- deterministic real-space tail/truncation and timing contract.
program test_tddft_realspace_truncation
   use precision_mod, only: rp
   use response_components_mod, only: RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel
   use tddft_chi0_mod, only: tddft_chi0_request, tddft_chi0_batch_result, tddft_chi0_result
   use tddft_chi0_realspace_mod, only: tddft_realspace_chi0_options, tddft_native_realspace_gf_provider
   implicit none

   call test_full_tail_and_production()
   call test_nested_source_radius()
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_full_tail_and_production()
      type(tddft_realspace_chi0_options) :: full_options, production_options
      type(tddft_chi0_result) :: full_result, production_result
      real(rp) :: q_points(3, 2), omega(3), full_seconds, production_seconds, speedup

      q_points = 0.0_rp
      q_points(1, 2) = 0.25_rp
      omega = [0.08_rp, 0.16_rp, 0.24_rp]
      full_options = truncation_options('full_tail')
      production_options = truncation_options('production')
      call run_source(192, full_options, q_points, omega, full_result)
      call run_source(192, production_options, q_points, omega, production_result)

      call assert(full_result%metadata%real_space_points == 9, 'full-tail retains the nine pairs inside R_max')
      call assert(full_result%metadata%real_space_omitted_points == 183, 'full-tail records every omitted source pair')
      call assert(full_result%metadata%real_space_tail_assessed, 'full-tail assesses the available response tail')
      call assert(full_result%metadata%real_space_source_covers_cutoff, 'full-tail source extends past R_max')
      call assert(trim(full_result%metadata%real_space_truncation_mode) == 'full_tail', 'full-tail mode is self-describing')
      call assert(production_result%metadata%real_space_points == 9, 'production retains the same selected pairs')
      call assert(production_result%metadata%real_space_omitted_points == 183, 'production reports discarded pairs')
      call assert(.not. production_result%metadata%real_space_tail_assessed, 'production does not infer an omitted tail')
      call assert(.not. production_result%metadata%converged, 'production never claims tail convergence')
      call assert(trim(production_result%metadata%real_space_truncation_mode) == 'production', &
         'production mode is self-describing')
      call assert(full_result%metadata%real_space_pair_response_integrations == 192*3, &
         'full-tail integrates every pair at every frequency')
      call assert(production_result%metadata%real_space_pair_response_integrations == 9*3, &
         'production skips pair integrations outside R_max')
      call assert(maxval(abs(full_result%chi-production_result%chi)) < 1.0e-12_rp, &
         'production and full-tail selected chi0 are identical')

      full_seconds = full_result%metadata%real_space_total_cpu_seconds
      production_seconds = production_result%metadata%real_space_total_cpu_seconds
      speedup = full_seconds/max(production_seconds, tiny(1.0_rp))
      write (*, '(a,4(1x,es16.8))') 'R2-04 timing full production speedup GF', full_seconds, production_seconds, speedup, &
         production_result%metadata%real_space_green_cpu_seconds
      write (*, '(a,4(1x,es16.8))') 'R2-04 phases full_pair production_pair full_FT production_FT', &
         full_result%metadata%real_space_pair_integration_cpu_seconds, &
         production_result%metadata%real_space_pair_integration_cpu_seconds, &
         full_result%metadata%real_space_fourier_cpu_seconds, production_result%metadata%real_space_fourier_cpu_seconds
   end subroutine test_full_tail_and_production

   subroutine test_nested_source_radius()
      type(tddft_realspace_chi0_options) :: options
      type(tddft_chi0_result) :: source_small, source_medium, source_large
      real(rp) :: q_points(3, 2), omega(2), medium_error, large_error

      q_points = 0.0_rp
      q_points(1, 2) = 0.25_rp
      omega = [0.08_rp, 0.16_rp]
      options = truncation_options('full_tail')
      call run_source(8, options, q_points, omega, source_small)
      call run_source(32, options, q_points, omega, source_medium)
      call run_source(192, options, q_points, omega, source_large)

      call assert(.not. source_small%metadata%real_space_source_covers_cutoff, &
         'short source is explicitly marked insufficient for R_max')
      call assert(.not. source_small%metadata%real_space_tail_assessed, 'short source cannot assess a tail')
      call assert(source_medium%metadata%real_space_source_radius > source_small%metadata%real_space_source_radius, &
         'nested source radius one is ordered')
      call assert(source_large%metadata%real_space_source_radius > source_medium%metadata%real_space_source_radius, &
         'nested source radius two is ordered')
      call assert(source_medium%metadata%real_space_tail_assessed, 'medium source provides an assessed tail')
      medium_error = maxval(abs(source_medium%chi-source_large%chi))/ &
         max(1.0e-14_rp, maxval(abs(source_large%chi)))
      large_error = maxval(abs(source_small%chi-source_large%chi))/ &
         max(1.0e-14_rp, maxval(abs(source_large%chi)))
      write (*, '(a,3(1x,es16.8))') 'R2-04 source radii small medium large', &
         source_small%metadata%real_space_source_radius, source_medium%metadata%real_space_source_radius, &
         source_large%metadata%real_space_source_radius
      write (*, '(a,2(1x,es16.8))') 'R2-04 chi0 source relative errors medium/short', medium_error, large_error
      call assert(medium_error < 5.0e-3_rp, 'nested source-radius chi0 convergence is demonstrated')
   end subroutine test_nested_source_radius

   function truncation_options(mode) result(options)
      character(len=*), intent(in) :: mode
      type(tddft_realspace_chi0_options) :: options

      options%eta = 0.10_rp
      options%green_eta = 0.05_rp
      options%fermi_level = 0.0_rp
      options%rmax = 1.0_rp
      options%tail_tolerance = 1.0e-3_rp
      options%truncation_mode = mode
      options%representation = 'bulk'
      options%metric = 0.0_rp
      options%metric(1, 1) = 1.0_rp
      options%metric(2, 2) = 1.0_rp
      options%metric(3, 3) = 1.0_rp
   end function truncation_options

   subroutine run_source(npair, options, q_points, omega, output)
      integer, intent(in) :: npair
      type(tddft_realspace_chi0_options), intent(in) :: options
      real(rp), intent(in) :: q_points(:, :), omega(:)
      type(tddft_chi0_result), intent(out) :: output
      type(tddft_native_realspace_gf_provider) :: provider
      type(tddft_chi0_request) :: request
      type(tddft_chi0_batch_result) :: batch
      type(response_channel) :: left(1), right(1)
      real(rp), allocatable :: energy(:), r_vectors(:, :)
      complex(rp), allocatable :: g_ab(:, :, :, :), g_ba(:, :, :, :)
      integer, allocatable :: pair_sites(:, :)
      integer :: counts(1), ie, ip, ne
      real(rp) :: radius, scale

      ne = 401
      allocate(energy(ne), r_vectors(3, npair), g_ab(2, 2, ne, npair), g_ba(2, 2, ne, npair), pair_sites(npair, 2))
      do ie = 1, ne
         energy(ie) = -1.0_rp + 2.0_rp*real(ie-1, rp)/real(ne-1, rp)
      end do
      r_vectors = 0.0_rp
      pair_sites = 1
      g_ab = cmplx(0.0_rp, 0.0_rp, rp)
      do ip = 1, npair
         radius = 0.125_rp*real(ip-1, rp)
         r_vectors(1, ip) = radius
         scale = exp(-1.2_rp*radius)
         do ie = 1, ne
            g_ab(1, 1, ie, ip) = scale/cmplx(energy(ie)+0.2_rp, 0.05_rp, rp)
            g_ab(2, 2, ie, ip) = scale/cmplx(energy(ie)-0.2_rp, 0.05_rp, rp)
         end do
      end do
      g_ba = 0.75_rp*g_ab
      counts = [1]
      left(1) = response_channel(1, RESPONSE_PLUS)
      right(1) = response_channel(1, RESPONSE_MINUS)
      call provider%initialize(energy, g_ab, g_ba, r_vectors, pair_sites, counts, left, right, options)
      allocate(request%q_points(3, size(q_points, 2)), request%omega(size(omega)))
      request%q_points = q_points
      request%omega = omega
      request%allow_real_space_reuse = .true.
      call provider%evaluate_realspace(request, batch)
      output = batch%q_response(1)
   end subroutine run_source

   subroutine assert(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) error stop 'FAIL: '//trim(message)
   end subroutine assert

end program test_tddft_realspace_truncation
