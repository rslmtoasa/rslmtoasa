!------------------------------------------------------------------------------
! TDDFT-09 -- three-backend equivalence campaign and invariant harness.
!------------------------------------------------------------------------------
! The fixture is deliberately built from a finite periodic spin-split model.
! Its K-space Lehmann source and native R-space source are generated from the
! same one-particle poles, but the three chi0 paths remain independent:
!
!   eigenpairs: explicit transition denominators;
!   K-GF:      energy-integrated Lehmann resolvents;
!   R-GF:      inverse discrete Fourier transformed G(R,z), followed by a
!              susceptibility transform (never a G(R)->G(k) runtime route).
!
! An optional first command-line argument names a JSON evidence file.  CTest
! runs the executable without an argument; developers can regenerate the
! committed evidence with:
!
!   build/bin/UnitTddftBackendEquivalence \
!      results/tddft_backend_equivalence.json
program test_tddft_backend_equivalence
   use precision_mod, only: rp
   use math_mod, only: pi
   use response_components_mod, only: RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, tddft_chi0_batch_result, &
      tddft_fermi_occupation, build_static_chi_ks_from_eigenpairs, build_static_chi_ks_from_eigenpairs_at_q
   use tddft_chi0_green_mod, only: green_chi0_options, eigenpair_green_function_provider, &
      build_chi_ks_from_green_functions, build_static_chi_ks_from_green_functions
   use tddft_chi0_realspace_mod, only: tddft_realspace_chi0_options, tddft_native_realspace_gf_provider
   use tddft_backend_mod, only: tddft_chi0_backend, tddft_eigenpair_backend, &
      tddft_kspace_lehmann_backend, tddft_realspace_gf_backend, make_tddft_chi0_backend, &
      tddft_backend_capabilities
   use tddft_ward_mod, only: tddft_ward_diagnostics, evaluate_static_ward_identity
   implicit none

   integer, parameter :: nq_main = 3, nw_main = 2
   integer, parameter :: nsum = 161
   real(rp), parameter :: main_eta = 0.02_rp
   real(rp), parameter :: pointwise_tolerance = 5.0e-2_rp
   real(rp), parameter :: kspace_realspace_tolerance = 1.0e-3_rp
   real(rp), parameter :: static_tolerance = 1.0e-11_rp
   real(rp), parameter :: native_static_tolerance = 8.0e-2_rp
   real(rp), parameter :: contour_tolerance = 5.0e-3_rp
   real(rp), parameter :: energy_window_min = -4.5_rp, energy_window_max = 4.5_rp
   real(rp), parameter :: pi2 = 2.0_rp*pi
   integer, parameter :: nk_levels(3) = [2, 4, 8]
   integer, parameter :: ne_levels(3) = [1001, 2001, 4001]
   integer, parameter :: contour_levels(3) = [16, 32, 64]
   integer, parameter :: eta_ne_levels(3) = [4001, 8001, 16001]
   real(rp), parameter :: rmax_levels(3) = [0.0_rp, 1.0_rp, 2.0_rp]
   real(rp), parameter :: eta_levels(3) = [0.04_rp, 0.02_rp, 0.01_rp]

   type :: campaign_result
      integer :: nk = 0, ne = 0, nq = 0, nw = 0
      real(rp) :: eta = 0.0_rp, rmax = 0.0_rp
      logical :: native_static_supported = .false.
      complex(rp), allocatable :: chi_eigen(:, :, :, :)
      complex(rp), allocatable :: chi_kspace(:, :, :, :)
      complex(rp), allocatable :: chi_realspace(:, :, :, :)
      real(rp), allocatable :: tail_ratio(:)
      integer :: native_builds = 0
   end type campaign_result

   real(rp) :: q_main(3, nq_main), omega_main(nw_main), q_zero(3, 1), q_static(3, 2), omega_zero(1)
   real(rp) :: q_two(3, 2), omega_two(nw_main), omega_sum(nsum)
   real(rp) :: k_errors(3), r_errors(3), energy_k_errors(3), energy_r_errors(3)
   real(rp) :: contour_errors(3), contour_relative_errors(3), eta_kr_errors(3), eta_ek_errors(3)
   real(rp) :: r_cutoff_errors(3), r_tail_ratios(3), sum_rule_values(3), sum_rule_residuals(3)
   real(rp) :: main_ek_max, main_er_max, main_kr_max, main_ek_eigen_max, main_er_eigen_max
   real(rp) :: static_error, static_k_error_q0, static_k_error_qfinite, static_r_error_q0, static_r_error_qfinite
   real(rp) :: static_ward_residuals(3), spectral_target, negative_sign_error, negative_factor_error
   complex(rp) :: bxc_static(2), magnetization(2)
   real(rp) :: zero_ward_residuals(3)
   logical :: negative_sign_detected, negative_factor_detected
   logical :: native_static_supported
   logical :: contour_passed, static_passed, sum_rule_passed
   character(len=512) :: evidence_path
   type(campaign_result) :: main_campaign, zero_campaign, sum_campaign, reference_campaign
   type(campaign_result) :: campaign
   type(tddft_chi0_result) :: static_eigen, static_kspace, static_realspace
   type(tddft_ward_diagnostics) :: static_eigen_ward, static_kspace_ward, static_realspace_ward
   type(tddft_ward_diagnostics) :: zero_eigen_ward, zero_kspace_ward, zero_realspace_ward
   logical :: failed
   integer :: iomega, ilevel, ios, command_length

   failed = .false.
   q_main = 0.0_rp
   q_main(1, :) = [0.0_rp, 0.25_rp, 0.5_rp]
   omega_main = [0.24_rp, 0.60_rp]
   q_zero = 0.0_rp
   q_static = 0.0_rp
   q_static(1, 2) = 0.25_rp
   q_two = 0.0_rp
   q_two(1, :) = [0.0_rp, 0.5_rp]
   omega_two = omega_main
   do iomega = 1, nsum
      omega_sum(iomega) = 0.80_rp*real(iomega-1, rp)/real(nsum-1, rp)
   end do

   ! Main pointwise campaign: all three backends see the same q/omega grid.
   call run_campaign(4, 8001, q_main, omega_main, main_eta, huge(1.0_rp), main_campaign)
   call compare_campaign(main_campaign, main_ek_max, main_er_max, main_kr_max, main_ek_eigen_max, main_er_eigen_max)

   ! Exact static is deliberately separate from finite-eta dynamic omega=0.
   ! Exercise all three independently: the R-GF source integrates its own
   ! retarded/advanced static contour identity, while the two K-space routes
   ! use their independent static divided-difference implementations.
   call make_static_fixture(static_eigen, static_kspace, q_zero, bxc_static, magnetization, static_error, &
      static_eigen_ward, static_kspace_ward)
   call make_native_static_fixture(static_realspace, q_static, bxc_static, magnetization, static_k_error_q0, &
      static_k_error_qfinite, static_r_error_q0, static_r_error_qfinite, static_realspace_ward)
   native_static_supported = main_campaign%native_static_supported
   static_passed = static_error <= static_tolerance .and. static_eigen_ward%ward_residual <= static_tolerance .and. &
      static_kspace_ward%ward_residual <= static_tolerance .and. static_k_error_q0 <= static_tolerance .and. &
      static_k_error_qfinite <= static_tolerance
   static_ward_residuals = [static_eigen_ward%ward_residual, static_kspace_ward%ward_residual, &
      static_realspace_ward%ward_residual]

   ! Zero-frequency dynamic diagnostics remain useful for the all-backend
   ! finite-eta Ward audit, but are reported separately from exact static.
   omega_zero = 0.0_rp
   call run_campaign(4, 2001, q_zero, omega_zero, main_eta, huge(1.0_rp), zero_campaign)
   call evaluate_zero_frequency_ward(zero_campaign, bxc_static, magnetization, zero_ward_residuals, &
      zero_eigen_ward, zero_kspace_ward, zero_realspace_ward)

   ! Positive-frequency integrated transverse spectral weight.  At T=0 the
   ! one-site fixture has target 4 in the selected m+/m- convention.
   call run_campaign(4, 2001, q_zero, omega_sum, main_eta, huge(1.0_rp), sum_campaign)
   spectral_target = 4.0_rp
   sum_rule_values(1) = integrate_positive_loss(sum_campaign%chi_eigen(:, :, :, 1), omega_sum)
   sum_rule_values(2) = integrate_positive_loss(sum_campaign%chi_kspace(:, :, :, 1), omega_sum)
   sum_rule_values(3) = integrate_positive_loss(sum_campaign%chi_realspace(:, :, :, 1), omega_sum)
   sum_rule_residuals = abs(sum_rule_values-spectral_target)/max(1.0_rp, abs(spectral_target))
   sum_rule_passed = maxval(sum_rule_residuals) <= 8.0e-2_rp

   ! Separate convergence axes.  Each row is a fresh run; no axis is hidden
   ! by reusing a converged output from another setting.
   do ilevel = 1, size(nk_levels)
      call run_campaign(nk_levels(ilevel), 1001, q_two, omega_two, main_eta, huge(1.0_rp), campaign)
      k_errors(ilevel) = max_backend_error(campaign%chi_kspace, campaign%chi_eigen)
      r_errors(ilevel) = max_backend_error(campaign%chi_realspace, campaign%chi_eigen)
   end do

   ! Native R cutoff convergence is compared with an independently generated
   ! full finite-source result.  R=0,1,2 are distinct shells for N=4.
   call run_campaign(4, 1001, q_two, omega_two, main_eta, huge(1.0_rp), reference_campaign)
   do ilevel = 1, size(rmax_levels)
      call run_campaign(4, 1001, q_two, omega_two, main_eta, rmax_levels(ilevel), campaign)
      r_cutoff_errors(ilevel) = max_backend_error(campaign%chi_realspace, reference_campaign%chi_realspace)
      r_tail_ratios(ilevel) = maxval(campaign%tail_ratio)
   end do

   ! Energy resolution is measured against a separately executed 4001-node
   ! reference for both GF realizations.
   call run_campaign(4, ne_levels(size(ne_levels)), q_zero, omega_two, main_eta, huge(1.0_rp), reference_campaign)
   do ilevel = 1, size(ne_levels)
      call run_campaign(4, ne_levels(ilevel), q_zero, omega_two, main_eta, huge(1.0_rp), campaign)
      energy_k_errors(ilevel) = max_backend_error(campaign%chi_kspace, reference_campaign%chi_kspace)
      energy_r_errors(ilevel) = max_backend_error(campaign%chi_realspace, reference_campaign%chi_realspace)
   end do

   ! gamma is tied to eta/2 in both GF sources.  The residual is measured at
   ! each ladder point instead of assuming broadening independence.
   do ilevel = 1, size(eta_levels)
      call run_campaign(4, eta_ne_levels(ilevel), q_two, omega_two, eta_levels(ilevel), huge(1.0_rp), campaign)
      eta_kr_errors(ilevel) = max_backend_error(campaign%chi_realspace, campaign%chi_kspace)
      eta_ek_errors(ilevel) = max_backend_error(campaign%chi_kspace, campaign%chi_eigen)
   end do

   call run_contour_convergence(contour_errors, contour_relative_errors)
   contour_passed = maxval(contour_relative_errors) <= contour_tolerance

   ! Test-only negative controls deliberately perturb a valid backend result.
   negative_sign_error = max_backend_error(-main_campaign%chi_kspace, main_campaign%chi_eigen)
   negative_factor_error = max_backend_error(2.0_rp*main_campaign%chi_kspace, main_campaign%chi_eigen)
   negative_sign_detected = negative_sign_error > pointwise_tolerance
   negative_factor_detected = negative_factor_error > pointwise_tolerance

   write (*, '(a,5(1x,es14.6))') 'TDDFT09 summary ek er kr eigk eigr', main_ek_max, main_er_max, main_kr_max, &
      main_ek_eigen_max, main_er_eigen_max
   write (*, '(a,3(1x,es14.6))') 'TDDFT09 eta kr', eta_kr_errors
   write (*, '(a,3(1x,es14.6))') 'TDDFT09 eta ek', eta_ek_errors
   write (*, '(a,3(1x,es14.6))') 'TDDFT09 contour rel', contour_relative_errors
   write (*, '(a,3(1x,es14.6))') 'TDDFT09 sum residual', sum_rule_residuals
   write (*, '(a,2(1x,es14.6))') 'TDDFT09 static K q0 qfinite', static_k_error_q0, static_k_error_qfinite
   write (*, '(a,2(1x,es14.6),a,3(1x,es14.6))') 'TDDFT09 static R q0 qfinite', static_r_error_q0, &
      static_r_error_qfinite, ' wards', static_ward_residuals

   call get_command_argument(1, evidence_path, length=command_length, status=ios)
   if (ios /= 0 .or. command_length == 0) evidence_path = ''
   if (len_trim(evidence_path) > 0) then
      call write_evidence(trim(evidence_path), main_campaign, omega_main, q_main, main_ek_max, main_er_max, main_kr_max, &
         main_ek_eigen_max, main_er_eigen_max, static_k_error_q0, static_k_error_qfinite, static_eigen_ward, &
         static_kspace_ward, static_realspace_ward, static_r_error_q0, static_r_error_qfinite, &
         native_static_supported, zero_ward_residuals, &
         sum_rule_values, sum_rule_residuals, spectral_target, &
         k_errors, r_errors, nk_levels, r_cutoff_errors, r_tail_ratios, rmax_levels, energy_k_errors, energy_r_errors, &
         ne_levels, contour_errors, contour_relative_errors, contour_levels, eta_kr_errors, eta_ek_errors, eta_levels, &
         eta_ne_levels, &
         negative_sign_error, negative_factor_error, negative_sign_detected, negative_factor_detected, &
         pointwise_tolerance, kspace_realspace_tolerance)
      write (*, '(a,1x,a)') 'TDDFT09 evidence written to', trim(evidence_path)
   end if

   call check_true('all main pointwise eigenpair/K-GF comparisons pass', main_ek_max <= pointwise_tolerance)
   call check_true('all main pointwise eigenpair/R-GF comparisons pass', main_er_max <= pointwise_tolerance)
   call check_true('all main K-GF/R-GF comparisons pass', main_kr_max <= kspace_realspace_tolerance)
   call check_true('main matrix eigenvalue comparisons pass', max(main_ek_eigen_max, main_er_eigen_max) <= pointwise_tolerance)
   call check_true('exact static eigenpair/K-GF comparison passes', static_passed)
   call check_true('all three backends advertise exact static support', native_static_supported)
   call check_true('native static metadata excludes dynamic eta', static_realspace%metadata%static_limit .and. &
      .not. static_realspace%metadata%eta_is_numerical .and. static_realspace%metadata%eta == 0.0_rp .and. &
      index(static_realspace%metadata%implementation, 'static contour identity') > 0)
   call check_true('exact static q=0 eigenpair/R-GF comparison passes', static_r_error_q0 <= native_static_tolerance)
   call check_true('exact static finite-q eigenpair/R-GF comparison passes', static_r_error_qfinite <= native_static_tolerance)
   call check_true('raw static Ward residuals agree across all three backends', &
      maxval(static_ward_residuals)-minval(static_ward_residuals) <= native_static_tolerance)
   call check_true('spectral sum rule is within finite-grid envelope', sum_rule_passed)
   call check_true('negative sign control is detected', negative_sign_detected)
   call check_true('negative factor-of-two control is detected', negative_factor_detected)
   call check_true('native real-space response is reused for q batches', main_campaign%native_builds == 1)
   call check_true('K-GF and native R-GF remain equivalent under eta/gamma ladder', &
      maxval(eta_kr_errors) <= kspace_realspace_tolerance)
   call check_true('contour-node convergence remains within justified envelope', contour_passed)

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine run_campaign(nk, ne, q_points, omega, eta_value, rmax, output)
      integer, intent(in) :: nk, ne
      real(rp), intent(in) :: q_points(:, :), omega(:), eta_value, rmax
      type(campaign_result), intent(out) :: output
      real(rp), allocatable :: weights(:), eval(:, :), evalq(:, :, :), energy(:), r_vectors(:, :)
      complex(rp), allocatable :: evec(:, :, :), evecq(:, :, :, :), g_ab(:, :, :, :), g_ba(:, :, :, :)
      integer, allocatable :: pair_sites(:, :)
      type(response_channel) :: left(2), right(2)
      type(tddft_chi0_options) :: eigen_options
      type(green_chi0_options) :: green_options
      type(tddft_realspace_chi0_options) :: realspace_options
      type(tddft_native_realspace_gf_provider) :: native_provider
      class(tddft_chi0_backend), allocatable :: eigen_backend, kspace_backend, realspace_backend
      type(tddft_chi0_batch_result) :: eigen_batch, kspace_batch, realspace_batch
      type(tddft_backend_capabilities) :: capabilities
      integer :: iq

      call make_fixture_arrays(nk, q_points, weights, eval, evalq, evec, evecq)
      call make_realspace_source(eval, ne, 0.5_rp*eta_value, energy, r_vectors, pair_sites, g_ab, g_ba)
      left = [response_channel(1, RESPONSE_PLUS), response_channel(1, RESPONSE_MINUS)]
      right = [response_channel(1, RESPONSE_MINUS), response_channel(1, RESPONSE_PLUS)]

      eigen_options%eta = eta_value
      eigen_options%fermi_level = 0.0_rp
      eigen_options%k_mesh_shape = [nk, 1, 1]
      eigen_options%response_projection = 'site'
      green_options%eta = eta_value
      green_options%green_eta = 0.5_rp*eta_value
      green_options%fermi_level = 0.0_rp
      green_options%energy_min = energy_window_min
      green_options%energy_max = energy_window_max
      green_options%energy_points = ne
      green_options%k_mesh_shape = [nk, 1, 1]
      realspace_options%eta = eta_value
      realspace_options%green_eta = 0.5_rp*eta_value
      realspace_options%fermi_level = 0.0_rp
      realspace_options%rmax = rmax
      realspace_options%tail_tolerance = 1.0e-2_rp
      realspace_options%representation = 'bulk'

      call native_provider%initialize(energy, g_ab, g_ba, r_vectors, pair_sites, [1], left, right, realspace_options)
      call make_tddft_chi0_backend('eigenpairs', eigen_backend)
      select type (eigen_backend)
      type is (tddft_eigenpair_backend)
         call eigen_backend%initialize(weights, eval, evec, evalq, evecq, q_points, [1], left, right, eigen_options)
      class default
         error stop 'TDDFT09: eigenpair factory returned wrong type'
      end select
      call make_tddft_chi0_backend('kspace_lehmann', kspace_backend)
      select type (kspace_backend)
      type is (tddft_kspace_lehmann_backend)
         call kspace_backend%initialize(weights, eval, evec, evalq, evecq, q_points, [1], left, right, eigen_options, green_options)
      class default
         error stop 'TDDFT09: K-GF factory returned wrong type'
      end select
      call make_tddft_chi0_backend('realspace_gf', realspace_backend)
      select type (realspace_backend)
      type is (tddft_realspace_gf_backend)
         call realspace_backend%initialize(native_provider)
         call realspace_backend%capabilities(capabilities)
      class default
         error stop 'TDDFT09: R-GF factory returned wrong type'
      end select

      call eigen_backend%evaluate_grid(q_points, omega, eigen_batch)
      call kspace_backend%evaluate_grid(q_points, omega, kspace_batch)
      call realspace_backend%evaluate_grid(q_points, omega, realspace_batch)
      output%nk = nk; output%ne = ne; output%nq = size(q_points, 2); output%nw = size(omega)
      output%eta = eta_value; output%rmax = rmax
      output%native_static_supported = capabilities%supports_static_limit
      allocate(output%chi_eigen(2, 2, output%nw, output%nq), output%chi_kspace(2, 2, output%nw, output%nq), &
         output%chi_realspace(2, 2, output%nw, output%nq), output%tail_ratio(output%nq))
      do iq = 1, output%nq
         output%chi_eigen(:, :, :, iq) = eigen_batch%q_response(iq)%chi
         output%chi_kspace(:, :, :, iq) = kspace_batch%q_response(iq)%chi
         output%chi_realspace(:, :, :, iq) = realspace_batch%q_response(iq)%chi
         output%tail_ratio(iq) = realspace_batch%q_response(iq)%metadata%real_space_tail_ratio
      end do
      select type (realspace_backend)
      type is (tddft_realspace_gf_backend)
         select type (provider => realspace_backend%provider)
         type is (tddft_native_realspace_gf_provider)
            output%native_builds = provider%build_count
         end select
      end select
   end subroutine run_campaign

   subroutine make_fixture_arrays(nk, q_points, weights, eval, evalq, evec, evecq)
      integer, intent(in) :: nk
      real(rp), intent(in) :: q_points(:, :)
      real(rp), allocatable, intent(out) :: weights(:), eval(:, :), evalq(:, :, :)
      complex(rp), allocatable, intent(out) :: evec(:, :, :), evecq(:, :, :, :)
      integer :: ik, iq
      real(rp) :: k, q

      allocate(weights(nk), eval(2, nk), evalq(2, nk, size(q_points, 2)), evec(2, 2, nk), &
         evecq(2, 2, nk, size(q_points, 2)))
      weights = 1.0_rp
      evec = cmplx(0.0_rp, 0.0_rp, rp); evecq = cmplx(0.0_rp, 0.0_rp, rp)
      do ik = 1, nk
         k = real(ik-1, rp)/real(nk, rp)
         eval(:, ik) = [dispersion(k, -1), dispersion(k, 1)]
         evec(1, 1, ik) = cmplx(1.0_rp, 0.0_rp, rp)
         evec(2, 2, ik) = cmplx(1.0_rp, 0.0_rp, rp)
         do iq = 1, size(q_points, 2)
            q = k + q_points(1, iq)
            evalq(:, ik, iq) = [dispersion(q, -1), dispersion(q, 1)]
            evecq(1, 1, ik, iq) = cmplx(1.0_rp, 0.0_rp, rp)
            evecq(2, 2, ik, iq) = cmplx(1.0_rp, 0.0_rp, rp)
         end do
      end do
   end subroutine make_fixture_arrays

   real(rp) function dispersion(k, spin) result(value)
      real(rp), intent(in) :: k
      integer, intent(in) :: spin
      value = real(spin, rp)*0.11_rp + 0.015_rp*cos(pi2*(k-floor(k)))
   end function dispersion

   subroutine make_realspace_source(eval, ne, gamma, energy, r_vectors, pair_sites, g_ab, g_ba)
      real(rp), intent(in) :: eval(:, :), gamma
      integer, intent(in) :: ne
      real(rp), allocatable, intent(out) :: energy(:), r_vectors(:, :)
      integer, allocatable, intent(out) :: pair_sites(:, :)
      complex(rp), allocatable, intent(out) :: g_ab(:, :, :, :), g_ba(:, :, :, :)
      integer :: nk, ir, ie, ik, nr
      real(rp) :: phase_argument, e
      complex(rp) :: phase

      nk = size(eval, 2); nr = nk
      allocate(energy(ne), r_vectors(3, nr), pair_sites(nr, 2), g_ab(2, 2, ne, nr), g_ba(2, 2, ne, nr))
      do ie = 1, ne
         energy(ie) = energy_window_min + (energy_window_max-energy_window_min)*real(ie-1, rp)/real(ne-1, rp)
      end do
      r_vectors = 0.0_rp
      do ir = 1, nr
         if (ir <= nk/2) then
            r_vectors(1, ir) = real(ir-1, rp)
         else
            r_vectors(1, ir) = real(ir-1-nk, rp)
         end if
      end do
      pair_sites = 1
      do ir = 1, nr
         do ie = 1, ne
            e = energy(ie)
            g_ab(:, :, ie, ir) = cmplx(0.0_rp, 0.0_rp, rp)
            g_ba(:, :, ie, ir) = cmplx(0.0_rp, 0.0_rp, rp)
            do ik = 1, nk
               phase_argument = pi2*real(ik-1, rp)*r_vectors(1, ir)/real(nk, rp)
               phase = exp(cmplx(0.0_rp, phase_argument, rp))
               g_ab(1, 1, ie, ir) = g_ab(1, 1, ie, ir) + phase/cmplx(e-eval(1, ik), gamma, rp)
               g_ab(2, 2, ie, ir) = g_ab(2, 2, ie, ir) + phase/cmplx(e-eval(2, ik), gamma, rp)
               phase = exp(cmplx(0.0_rp, -phase_argument, rp))
               g_ba(1, 1, ie, ir) = g_ba(1, 1, ie, ir) + phase/cmplx(e-eval(1, ik), gamma, rp)
               g_ba(2, 2, ie, ir) = g_ba(2, 2, ie, ir) + phase/cmplx(e-eval(2, ik), gamma, rp)
            end do
            g_ab(:, :, ie, ir) = g_ab(:, :, ie, ir)/real(nk, rp)
            g_ba(:, :, ie, ir) = g_ba(:, :, ie, ir)/real(nk, rp)
         end do
      end do
   end subroutine make_realspace_source

   subroutine make_static_fixture(static_eigen, static_kspace, q_zero, bxc, magnetization, error, eigen_ward, kspace_ward)
      type(tddft_chi0_result), intent(out) :: static_eigen, static_kspace
      real(rp), intent(in) :: q_zero(:, :)
      complex(rp), intent(out) :: bxc(:), magnetization(:)
      real(rp), intent(out) :: error
      type(tddft_ward_diagnostics), intent(out) :: eigen_ward, kspace_ward
      real(rp), allocatable :: weights(:), eval(:, :), evalq(:, :, :)
      complex(rp), allocatable :: evec(:, :, :), evecq(:, :, :, :)
      type(response_channel) :: left(2), right(2)
      type(tddft_chi0_options) :: options
      type(green_chi0_options) :: green_options
      type(eigenpair_green_function_provider) :: source

      call make_fixture_arrays(4, q_zero, weights, eval, evalq, evec, evecq)
      left = [response_channel(1, RESPONSE_PLUS), response_channel(1, RESPONSE_MINUS)]
      right = [response_channel(1, RESPONSE_MINUS), response_channel(1, RESPONSE_PLUS)]
      options%fermi_level = 0.0_rp; options%k_mesh_shape = [4, 1, 1]; options%q_direct = q_zero(:, 1)
      green_options%fermi_level = 0.0_rp; green_options%q_direct = q_zero(:, 1); green_options%k_mesh_shape = [4, 1, 1]
      call build_static_chi_ks_from_eigenpairs(weights, eval, evec, [1], left, right, options, static_eigen)
      call source%initialize(eval, evec, eval, evec)
      call build_static_chi_ks_from_green_functions(source, weights, [1], left, right, green_options, static_kspace)
      error = matrix_error(static_kspace%chi(:, :, 1), static_eigen%chi(:, :, 1))
      magnetization = [cmplx(1.0_rp, 0.0_rp, rp), cmplx(0.0_rp, 0.0_rp, rp)]
      bxc = cmplx(0.0_rp, 0.0_rp, rp); bxc(1) = 1.0_rp/static_eigen%chi(1, 1, 1)
      call evaluate_static_ward_identity(static_eigen%chi(:, :, 1), bxc, magnetization, eigen_ward, &
         response_basis='site circular fixture', bxc_provenance='static eigenpair fixture', &
         kernel_provenance='static fixture Ward field')
      call evaluate_static_ward_identity(static_kspace%chi(:, :, 1), bxc, magnetization, kspace_ward, &
         response_basis='site circular fixture', bxc_provenance='static eigenpair fixture', &
         kernel_provenance='static fixture Ward field')
   end subroutine make_static_fixture

   subroutine make_native_static_fixture(native_result, q_points, bxc, magnetization, k_error_q0, k_error_qfinite, &
      error_q0, error_qfinite, ward)
      type(tddft_chi0_result), intent(out) :: native_result
      real(rp), intent(in) :: q_points(:, :)
      complex(rp), intent(in) :: bxc(:), magnetization(:)
      real(rp), intent(out) :: k_error_q0, k_error_qfinite, error_q0, error_qfinite
      type(tddft_ward_diagnostics), intent(out) :: ward
      real(rp), allocatable :: weights(:), eval(:, :), evalq(:, :, :), energy(:), r_vectors(:, :)
      complex(rp), allocatable :: evec(:, :, :), evecq(:, :, :, :), g_ab(:, :, :, :), g_ba(:, :, :, :)
      integer, allocatable :: pair_sites(:, :)
      type(response_channel) :: left(2), right(2)
      type(tddft_chi0_options) :: static_options
      type(tddft_realspace_chi0_options) :: realspace_options
      type(tddft_native_realspace_gf_provider) :: native_provider
      class(tddft_chi0_backend), allocatable :: realspace_backend, kspace_backend
      type(tddft_chi0_batch_result) :: native_batch
      type(tddft_chi0_batch_result) :: kspace_batch
      type(green_chi0_options) :: green_options
      type(tddft_chi0_result) :: reference_q0, reference_qfinite
      integer :: iq

      if (size(q_points, 2) /= 2) error stop 'TDDFT09: native static fixture requires q=0 and finite-q'
      call make_fixture_arrays(4, q_points, weights, eval, evalq, evec, evecq)
      call make_realspace_source(eval, 32001, 5.0e-4_rp, energy, r_vectors, pair_sites, g_ab, g_ba)
      left = [response_channel(1, RESPONSE_PLUS), response_channel(1, RESPONSE_MINUS)]
      right = [response_channel(1, RESPONSE_MINUS), response_channel(1, RESPONSE_PLUS)]

      realspace_options%eta = main_eta
      realspace_options%green_eta = 0.0_rp
      realspace_options%fermi_level = 0.0_rp
      realspace_options%rmax = huge(1.0_rp)
      realspace_options%tail_tolerance = 1.0e-2_rp
      realspace_options%representation = 'bulk'
      call native_provider%initialize(energy, g_ab, g_ba, r_vectors, pair_sites, [1], left, right, realspace_options)
      static_options%eta = 0.0_rp
      static_options%fermi_level = 0.0_rp
      static_options%k_mesh_shape = [4, 1, 1]
      call make_tddft_chi0_backend('realspace_gf', realspace_backend)
      select type (realspace_backend)
      type is (tddft_realspace_gf_backend)
         call realspace_backend%initialize(native_provider)
         call realspace_backend%evaluate_static_grid(q_points, native_batch)
      class default
         error stop 'TDDFT09: R-GF static fixture factory returned wrong type'
      end select
      native_result = native_batch%q_response(1)

      green_options%eta = main_eta
      green_options%green_eta = 0.0_rp
      green_options%fermi_level = 0.0_rp
      green_options%k_mesh_shape = [4, 1, 1]
      call make_tddft_chi0_backend('kspace_lehmann', kspace_backend)
      select type (kspace_backend)
      type is (tddft_kspace_lehmann_backend)
         call kspace_backend%initialize(weights, eval, evec, evalq, evecq, q_points, [1], left, right, static_options, &
            green_options)
         call kspace_backend%evaluate_static_grid(q_points, kspace_batch)
      class default
         error stop 'TDDFT09: K-GF static fixture factory returned wrong type'
      end select

      do iq = 1, 2
         static_options%q_direct = q_points(:, iq)
         if (iq == 1) then
            call build_static_chi_ks_from_eigenpairs_at_q(weights, eval, evec, evalq(:, :, iq), evecq(:, :, :, iq), [1], &
               left, right, static_options, reference_q0)
         else
            call build_static_chi_ks_from_eigenpairs_at_q(weights, eval, evec, evalq(:, :, iq), evecq(:, :, :, iq), [1], &
               left, right, static_options, reference_qfinite)
         end if
      end do
      error_q0 = matrix_error(native_batch%q_response(1)%chi(:, :, 1), reference_q0%chi(:, :, 1))
      error_qfinite = matrix_error(native_batch%q_response(2)%chi(:, :, 1), reference_qfinite%chi(:, :, 1))
      k_error_q0 = matrix_error(kspace_batch%q_response(1)%chi(:, :, 1), reference_q0%chi(:, :, 1))
      k_error_qfinite = matrix_error(kspace_batch%q_response(2)%chi(:, :, 1), reference_qfinite%chi(:, :, 1))
      call evaluate_static_ward_identity(native_batch%q_response(1)%chi(:, :, 1), bxc, magnetization, ward, &
         response_basis='site circular fixture', bxc_provenance='exact static eigenpair fixture field', &
         kernel_provenance='unadjusted static fixture Ward field')
   end subroutine make_native_static_fixture

   subroutine evaluate_zero_frequency_ward(campaign, bxc, magnetization, residuals, eigen_ward, kspace_ward, realspace_ward)
      type(campaign_result), intent(in) :: campaign
      complex(rp), intent(in) :: bxc(:), magnetization(:)
      real(rp), intent(out) :: residuals(3)
      type(tddft_ward_diagnostics), intent(out) :: eigen_ward, kspace_ward, realspace_ward

      call evaluate_static_ward_identity(campaign%chi_eigen(:, :, 1, 1), bxc, magnetization, eigen_ward, &
         response_basis='site circular fixture', bxc_provenance='exact static fixture field', &
         kernel_provenance='finite-eta dynamic diagnostic')
      call evaluate_static_ward_identity(campaign%chi_kspace(:, :, 1, 1), bxc, magnetization, kspace_ward, &
         response_basis='site circular fixture', bxc_provenance='exact static fixture field', &
         kernel_provenance='finite-eta dynamic diagnostic')
      call evaluate_static_ward_identity(campaign%chi_realspace(:, :, 1, 1), bxc, magnetization, realspace_ward, &
         response_basis='site circular fixture', bxc_provenance='exact static fixture field', &
         kernel_provenance='finite-eta dynamic diagnostic')
      residuals = [eigen_ward%ward_residual, kspace_ward%ward_residual, realspace_ward%ward_residual]
   end subroutine evaluate_zero_frequency_ward

   subroutine run_contour_convergence(errors, relative_errors)
      real(rp), intent(out) :: errors(:), relative_errors(:)
      real(rp) :: q_zero_local(3, 1), omega_local(2)
      real(rp), allocatable :: weights(:), eval(:, :), evalq(:, :, :)
      complex(rp), allocatable :: evec(:, :, :), evecq(:, :, :, :)
      type(response_channel) :: left(2), right(2)
      type(eigenpair_green_function_provider) :: source
      type(green_chi0_options) :: direct_options, mixed_options
      type(tddft_chi0_result) :: direct, mixed
      integer :: ilevel

      q_zero_local = 0.0_rp; omega_local = [0.24_rp, 0.60_rp]
      call make_fixture_arrays(4, q_zero_local, weights, eval, evalq, evec, evecq)
      left = [response_channel(1, RESPONSE_PLUS), response_channel(1, RESPONSE_MINUS)]
      right = [response_channel(1, RESPONSE_MINUS), response_channel(1, RESPONSE_PLUS)]
      call source%initialize(eval, evec, eval, evec)
      direct_options%eta = main_eta; direct_options%green_eta = 0.5_rp*main_eta; direct_options%fermi_level = 0.0_rp
      direct_options%energy_min = energy_window_min; direct_options%energy_max = energy_window_max; direct_options%energy_points = 8001
      direct_options%q_direct = 0.0_rp; direct_options%k_mesh_shape = [4, 1, 1]
      call build_chi_ks_from_green_functions(source, weights, [1], left, right, omega_local, direct_options, direct)
      do ilevel = 1, size(contour_levels)
         mixed_options = direct_options; mixed_options%energy_integration = 'mixed_contour'
         mixed_options%contour_points = contour_levels(ilevel); mixed_options%contour_subdivisions = 8
         mixed_options%near_fermi_points = 128
         call build_chi_ks_from_green_functions(source, weights, [1], left, right, omega_local, mixed_options, mixed)
         errors(ilevel) = max_response_error(mixed%chi, direct%chi)
         relative_errors(ilevel) = errors(ilevel)/max(1.0_rp, maxval(abs(direct%chi)))
      end do
   end subroutine run_contour_convergence

   real(rp) function integrate_positive_loss(chi, omega) result(integral)
      complex(rp), intent(in) :: chi(:, :, :)
      real(rp), intent(in) :: omega(:)
      integer :: iw
      real(rp) :: h, previous, current

      integral = 0.0_rp
      do iw = 1, size(omega)-1
         h = omega(iw+1)-omega(iw)
         previous = -aimag(chi(1, 1, iw))/pi
         current = -aimag(chi(1, 1, iw+1))/pi
         integral = integral + 0.5_rp*h*(previous+current)
      end do
   end function integrate_positive_loss

   subroutine compare_campaign(campaign, max_ek, max_er, max_kr, max_ek_eigen, max_er_eigen)
      type(campaign_result), intent(in) :: campaign
      real(rp), intent(out) :: max_ek, max_er, max_kr, max_ek_eigen, max_er_eigen
      integer :: iq, iw

      max_ek = 0.0_rp; max_er = 0.0_rp; max_kr = 0.0_rp; max_ek_eigen = 0.0_rp; max_er_eigen = 0.0_rp
      do iq = 1, campaign%nq
         do iw = 1, campaign%nw
            max_ek = max(max_ek, matrix_error(campaign%chi_kspace(:, :, iw, iq), campaign%chi_eigen(:, :, iw, iq)))
            max_er = max(max_er, matrix_error(campaign%chi_realspace(:, :, iw, iq), campaign%chi_eigen(:, :, iw, iq)))
            max_kr = max(max_kr, matrix_error(campaign%chi_realspace(:, :, iw, iq), campaign%chi_kspace(:, :, iw, iq)))
            max_ek_eigen = max(max_ek_eigen, matrix_eigen_error(campaign%chi_kspace(:, :, iw, iq), &
               campaign%chi_eigen(:, :, iw, iq)))
            max_er_eigen = max(max_er_eigen, matrix_eigen_error(campaign%chi_realspace(:, :, iw, iq), &
               campaign%chi_eigen(:, :, iw, iq)))
         end do
      end do
   end subroutine compare_campaign

   real(rp) function max_backend_error(actual, reference) result(value)
      complex(rp), intent(in) :: actual(:, :, :, :), reference(:, :, :, :)
      integer :: iq, iw

      if (any(shape(actual) /= shape(reference))) error stop 'TDDFT09: comparison shape mismatch'
      value = 0.0_rp
      do iq = 1, size(actual, 4)
         do iw = 1, size(actual, 3)
            value = max(value, matrix_error(actual(:, :, iw, iq), reference(:, :, iw, iq)))
         end do
      end do
   end function max_backend_error

   real(rp) function max_response_error(actual, reference) result(value)
      complex(rp), intent(in) :: actual(:, :, :), reference(:, :, :)
      integer :: iw

      if (any(shape(actual) /= shape(reference))) error stop 'TDDFT09: response comparison shape mismatch'
      value = 0.0_rp
      do iw = 1, size(actual, 3)
         value = max(value, matrix_error(actual(:, :, iw), reference(:, :, iw)))
      end do
   end function max_response_error

   real(rp) function matrix_error(actual, reference) result(value)
      complex(rp), intent(in) :: actual(:, :), reference(:, :)
      value = sqrt(sum(abs(actual-reference)**2))/max(1.0_rp, sqrt(sum(abs(reference)**2)))
   end function matrix_error

   real(rp) function matrix_eigen_error(actual, reference) result(value)
      complex(rp), intent(in) :: actual(:, :), reference(:, :)
      complex(rp), allocatable :: actual_eigen(:), reference_eigen(:)
      real(rp) :: distance, scale
      integer :: i, j

      call matrix_eigenvalues(actual, actual_eigen)
      call matrix_eigenvalues(reference, reference_eigen)
      ! Use the same global matrix scale as the Frobenius comparison.  This
      ! checks every eigenmode without turning a small eigenvalue into an
      ! artificial large relative error.
      scale = max(1.0_rp, sqrt(sum(abs(reference)**2)))
      value = 0.0_rp
      do i = 1, size(actual_eigen)
         distance = huge(1.0_rp)
         do j = 1, size(reference_eigen)
            distance = min(distance, abs(actual_eigen(i)-reference_eigen(j)))
         end do
         value = max(value, distance/scale)
      end do
   end function matrix_eigen_error

   subroutine matrix_eigenvalues(matrix, eigenvalues)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp), allocatable, intent(out) :: eigenvalues(:)
      complex(rp), allocatable :: work_matrix(:, :), work(:)
      real(rp), allocatable :: rwork(:)
      complex(rp) :: work_query(1)
      integer :: n, info, lwork

      n = size(matrix, 1)
      if (n < 1 .or. size(matrix, 2) /= n) error stop 'TDDFT09: eigenvalue matrix is not square'
      allocate(work_matrix(n, n), eigenvalues(n), rwork(max(1, 2*n)))
      call zgeev('N', 'N', n, work_matrix, n, eigenvalues, work_matrix, 1, work_matrix, 1, work_query, -1, rwork, info)
      if (info /= 0) error stop 'TDDFT09: LAPACK eigenvalue workspace query failed'
      lwork = max(1, int(real(work_query(1), rp)))
      allocate(work(lwork))
      work_matrix = matrix
      call zgeev('N', 'N', n, work_matrix, n, eigenvalues, work_matrix, 1, work_matrix, 1, work, lwork, rwork, info)
      if (info /= 0) error stop 'TDDFT09: LAPACK eigenvalue solve failed'
   end subroutine matrix_eigenvalues

   subroutine check_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine check_true

   subroutine write_evidence(path, main, omega, q_points, max_ek, max_er, max_kr, max_ek_eigen, max_er_eigen, &
      static_k_error_q0, static_k_error_qfinite, static_eigen_ward, static_kspace_ward, static_realspace_ward, &
      static_r_error_q0, static_r_error_qfinite, &
      native_static_supported, zero_ward, sum_values, sum_residuals, &
      sum_target, k_errors, r_errors, nk_values, r_errors_cutoff, r_tail, r_values, energy_k, energy_r, ne_values, &
      contour_errors_local, contour_relative, contour_values, eta_kr, eta_ek, eta_values, eta_energy_points, negative_sign, negative_factor, &
      sign_detected, factor_detected, comparison_tol, kr_tol)
      character(len=*), intent(in) :: path
      type(campaign_result), intent(in) :: main
      real(rp), intent(in) :: omega(:), q_points(:, :), max_ek, max_er, max_kr, max_ek_eigen, max_er_eigen
      real(rp), intent(in) :: static_k_error_q0, static_k_error_qfinite
      type(tddft_ward_diagnostics), intent(in) :: static_eigen_ward, static_kspace_ward, static_realspace_ward
      real(rp), intent(in) :: static_r_error_q0, static_r_error_qfinite
      logical, intent(in) :: native_static_supported
      real(rp), intent(in) :: zero_ward(:), sum_values(:), sum_residuals(:), sum_target, k_errors(:), r_errors(:)
      integer, intent(in) :: nk_values(:)
      real(rp), intent(in) :: r_errors_cutoff(:), r_tail(:), r_values(:), energy_k(:), energy_r(:)
      integer, intent(in) :: ne_values(:)
      real(rp), intent(in) :: contour_errors_local(:), contour_relative(:)
      integer, intent(in) :: contour_values(:)
      real(rp), intent(in) :: eta_kr(:), eta_ek(:), eta_values(:), negative_sign, negative_factor
      integer, intent(in) :: eta_energy_points(:)
      logical, intent(in) :: sign_detected, factor_detected
      real(rp), intent(in) :: comparison_tol, kr_tol
      integer :: unit, ios, iq, iw, ilevel

      open(newunit=unit, file=trim(path), status='replace', action='write', iostat=ios)
      if (ios /= 0) error stop 'TDDFT09: cannot open evidence path'
      write(unit, '(a)') '{'
      write(unit, '(a)') '  "campaign": "TDDFT-09",'
      write(unit, '(a)') '  "fixture": {"kind": "periodic_one_site_spin_split", "response_basis": "site circular", "temperature_K": 0.0},'
      write(unit, '(a,i0,a,i0,a,i0,a,i0,a)') '  "settings": {"nk": ', main%nk, ', "energy_points": ', main%ne, &
         ', "q_points": ', main%nq, ', "omega_points": ', main%nw, ','
      write(unit, '(a,es24.16,a,es24.16,a)') '    "eta": ', main%eta, ', "green_eta": ', 0.5_rp*main%eta, ','
      write(unit, '(a,es24.16,a,es24.16,a)') '    "pointwise_relative_tolerance": ', comparison_tol, &
         ', "kspace_realspace_tolerance": ', kr_tol, ','
      write(unit, '(a,es24.16,a,es24.16,a)') '    "energy_window": [', energy_window_min, ', ', energy_window_max, &
         '], "energy_unit": "Rydberg", "eta_role": "numerical"},'
      write(unit, '(a)') '  "pointwise": ['
      do iq = 1, main%nq
         do iw = 1, main%nw
            call write_point(unit, iq, iw, q_points(1, iq), omega(iw), main%chi_eigen(:, :, iw, iq), &
               main%chi_kspace(:, :, iw, iq), main%chi_realspace(:, :, iw, iq), iw == main%nw .and. iq == main%nq)
         end do
      end do
      write(unit, '(a)') '  ],'
      write(unit, '(a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a)') &
         '  "summary": {"max_eigenpair_vs_kspace": ', max_ek, &
         ', "max_eigenpair_vs_realspace": ', max_er, ', "max_kspace_vs_realspace": ', max_kr, &
         ', "max_eigenvalue_vs_kspace": ', max_ek_eigen, ', "max_eigenvalue_vs_realspace": ', max_er_eigen, '},'
      write(unit, '(a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,a,a)') &
         '  "static": {"eigenpair_vs_kspace": ', static_k_error_q0, ', "eigenpair_vs_kspace_finite_q": ', &
         static_k_error_qfinite, ', "eigenpair_ward_residual": ', static_eigen_ward%ward_residual, &
         ', "kspace_ward_residual": ', static_kspace_ward%ward_residual, ', "eigenpair_vs_native_rgf_q0": ', &
         static_r_error_q0, ', "eigenpair_vs_native_rgf_finite_q": ', static_r_error_qfinite, &
         ', "native_ward_residual": ', static_realspace_ward%ward_residual, &
         ', "native_static_supported": ', bool_token(native_static_supported), '},'
      write(unit, '(a,i0,a,a,a)') '  "native_realspace": {"response_build_count": ', main%native_builds, &
         ', "q_batch_reuse": ', bool_token(main%native_builds == 1), '},'
      write(unit, '(a,es24.16,a,es24.16,a,es24.16,a)') '  "zero_frequency_dynamic_ward_residual": [', &
         zero_ward(1), ',', zero_ward(2), ',', zero_ward(3), '],'
      write(unit, '(a,es24.16,a)') '  "spectral_sum_rule": {"target": ', sum_target, ', "values": ['
      write(unit, '(a,es24.16,a,es24.16,a,es24.16,a)') '    ', sum_values(1), ',', sum_values(2), ',', sum_values(3), '],'
      write(unit, '(a,es24.16,a,es24.16,a,es24.16,a)') '    "relative_residuals": [', &
         sum_residuals(1), ',', sum_residuals(2), ',', sum_residuals(3), ']},'
      write(unit, '(a)') '  "convergence": {'
      write(unit, '(a)') '    "k_mesh": ['
      do ilevel = 1, size(nk_values)
         call write_convergence_row(unit, nk_values(ilevel), k_errors(ilevel), r_errors(ilevel), ilevel == size(nk_values))
      end do
      write(unit, '(a)') '    ],'
      write(unit, '(a)') '    "r_cutoff": ['
      do ilevel = 1, size(r_values)
         call write_r_row(unit, r_values(ilevel), r_errors_cutoff(ilevel), r_tail(ilevel), ilevel == size(r_values))
      end do
      write(unit, '(a)') '    ],'
      write(unit, '(a)') '    "energy_resolution": ['
      do ilevel = 1, size(ne_values)
         call write_energy_row(unit, ne_values(ilevel), energy_k(ilevel), energy_r(ilevel), ilevel == size(ne_values))
      end do
      write(unit, '(a)') '    ],'
      write(unit, '(a)') '    "contour": ['
      do ilevel = 1, size(contour_values)
         call write_contour_row(unit, contour_values(ilevel), contour_errors_local(ilevel), contour_relative(ilevel), &
            ilevel == size(contour_values))
      end do
      write(unit, '(a)') '    ],'
      write(unit, '(a)') '    "eta_gamma": ['
      do ilevel = 1, size(eta_values)
         call write_eta_row(unit, eta_values(ilevel), 0.5_rp*eta_values(ilevel), eta_energy_points(ilevel), &
            eta_kr(ilevel), eta_ek(ilevel), ilevel == size(eta_values))
      end do
      write(unit, '(a)') '    ]},'
      write(unit, '(a,es24.16,a,es24.16,a,a,a,a,a)') '  "negative_controls": {"sign_flip_error": ', negative_sign, &
         ', "factor_two_error": ', negative_factor, ', "sign_flip_detected": ', bool_token(sign_detected), &
         ', "factor_two_detected": ', bool_token(factor_detected), '}'
      write(unit, '(a)') '}'
      close(unit)
   end subroutine write_evidence

   subroutine write_point(unit, iq, iw, qx, omega, eigen, kspace, realspace, last)
      integer, intent(in) :: unit, iq, iw
      real(rp), intent(in) :: qx, omega
      complex(rp), intent(in) :: eigen(:, :), kspace(:, :), realspace(:, :)
      logical, intent(in) :: last
      character(len=:), allocatable :: comma

      if (last) then; comma = ''; else; comma = ','; end if
      write(unit, '(a,i0,a,i0,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a)') &
         '    {"q_index": ', iq, ', "omega_index": ', iw, ', "q_x": ', qx, ', "omega": ', omega, &
         ', "eigenpair_vs_kspace": ', matrix_error(kspace, eigen), ', "eigenpair_vs_realspace": ', &
         matrix_error(realspace, eigen), ', "kspace_vs_realspace": ', matrix_error(realspace, kspace), &
         ', "eigenvalue_vs_kspace": ', matrix_eigen_error(kspace, eigen), ', "eigenvalue_vs_realspace": ', &
         matrix_eigen_error(realspace, eigen), '}'//comma
   end subroutine write_point

   subroutine write_convergence_row(unit, n, eigen_k, eigen_r, last)
      integer, intent(in) :: unit, n
      real(rp), intent(in) :: eigen_k, eigen_r
      logical, intent(in) :: last
      character(len=:), allocatable :: comma
      if (last) then; comma = ''; else; comma = ','; end if
      write(unit, '(a,i0,a,es24.16,a,es24.16,a)') '      {"nk": ', n, ', "eigenpair_vs_kspace": ', eigen_k, &
         ', "eigenpair_vs_realspace": ', eigen_r, '}'//comma
   end subroutine write_convergence_row

   subroutine write_r_row(unit, rmax, error, tail, last)
      integer, intent(in) :: unit
      real(rp), intent(in) :: rmax, error, tail
      logical, intent(in) :: last
      character(len=:), allocatable :: comma
      if (last) then; comma = ''; else; comma = ','; end if
      write(unit, '(a,es24.16,a,es24.16,a,es24.16,a)') '      {"rmax": ', rmax, ', "relative_to_full": ', error, &
         ', "tail_ratio": ', tail, '}'//comma
   end subroutine write_r_row

   subroutine write_energy_row(unit, ne, error_k, error_r, last)
      integer, intent(in) :: unit, ne
      real(rp), intent(in) :: error_k, error_r
      logical, intent(in) :: last
      character(len=:), allocatable :: comma
      if (last) then; comma = ''; else; comma = ','; end if
      write(unit, '(a,i0,a,es24.16,a,es24.16,a)') '      {"energy_points": ', ne, ', "kspace_change": ', error_k, &
         ', "realspace_change": ', error_r, '}'//comma
   end subroutine write_energy_row

   subroutine write_contour_row(unit, points, error, relative, last)
      integer, intent(in) :: unit, points
      real(rp), intent(in) :: error, relative
      logical, intent(in) :: last
      character(len=:), allocatable :: comma
      if (last) then; comma = ''; else; comma = ','; end if
      write(unit, '(a,i0,a,es24.16,a,es24.16,a)') '      {"contour_points": ', points, ', "absolute_error": ', error, &
         ', "relative_error": ', relative, '}'//comma
   end subroutine write_contour_row

   subroutine write_eta_row(unit, eta_value, gamma_value, energy_points, error_kr, error_ek, last)
      integer, intent(in) :: unit
      integer, intent(in) :: energy_points
      real(rp), intent(in) :: eta_value, gamma_value, error_kr, error_ek
      logical, intent(in) :: last
      character(len=:), allocatable :: comma
      if (last) then; comma = ''; else; comma = ','; end if
      write(unit, '(a,es24.16,a,es24.16,a,i0,a,es24.16,a,es24.16,a,es24.16,a)') '      {"eta": ', eta_value, ', "gamma": ', gamma_value, &
         ', "energy_points": ', energy_points, ', "gamma_over_eta": ', gamma_value/eta_value, ', "kspace_vs_realspace": ', error_kr, &
         ', "eigenpair_vs_kspace": ', error_ek, '}'//comma
   end subroutine write_eta_row

   function bool_token(value) result(token)
      logical, intent(in) :: value
      character(len=5) :: token

      if (value) then
         token = 'true '
      else
         token = 'false'
      end if
   end function bool_token

end program test_tddft_backend_equivalence
