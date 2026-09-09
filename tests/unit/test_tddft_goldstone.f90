!------------------------------------------------------------------------------
! TDDFT-04 -- transverse ALSDA kernel and transparent Goldstone diagnostics
!------------------------------------------------------------------------------
program test_tddft_goldstone
   use precision_mod, only: rp
   use xc_response_kernel_mod, only: xc_response_kernel_provider, build_independent_circular_bxc_source
   use tddft_goldstone_mod, only: tddft_goldstone_options, tddft_goldstone_result, &
      tddft_goldstone_diagnostics, tddft_xi_compare_diagnostics, build_site_projected_k_perp, construct_transverse_xi, &
      evaluate_goldstone, compare_xi_goldstone_branches, &
      tddft_goldstone_column_correction, evaluate_raw_xi_diagnostics, build_goldstone_column_correction, &
      rescale_xi_columns, spectral_weights_are_nonnegative, spectral_weight_correction_is_acceptable, &
      write_goldstone_diagnostics_text
   implicit none

   real(rp), parameter :: tol = 2.0e-12_rp
   logical :: failed

   failed = .false.
   call test_xc_provider_kernel_and_raw_diagnostics()
   call test_off_mode_only_constructs_xi()
   call test_controlled_pair_column_correction()
   call test_correction_rejections_and_spectral_control()
   call test_symmetry_breaking_disables_sum_rule()
   call test_raw_residual_convergence_controls()
   call test_signed_magnetization_diagnostics()
   call test_explicit_goldstone_policies()
   call test_independent_bxc_source_and_ward()
   call test_anti_goldstone_is_phase_aware()
   call test_compare_branch_policy()
   call test_diagnostic_output_contract()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_xc_provider_kernel_and_raw_diagnostics()
      type(xc_response_kernel_provider) :: provider
      type(tddft_goldstone_options) :: options
      type(tddft_goldstone_result) :: result
      complex(rp), allocatable :: k_perp(:)
      complex(rp) :: chi(2, 2), xi(2, 2)

      call make_provider(provider, [2.0_rp, 1.0_rp], [3.0_rp, 4.0_rp])
      call build_site_projected_k_perp(provider, k_perp)
      call assert_real_vector('K_perp comes from xc_response_kernel provider', real(k_perp, rp), [3.0_rp, 4.0_rp])
      chi = cmplx(0.0_rp, 0.0_rp, rp)
      chi(1, 1) = cmplx(0.3_rp, 0.0_rp, rp)
      chi(2, 2) = cmplx(0.2_rp, 0.0_rp, rp)
      xi = construct_transverse_xi(chi, k_perp)
      call assert_complex('Xi(1,1) = chi_KS K_perp', xi(1, 1), cmplx(0.9_rp, 0.0_rp, rp))
      call assert_complex('Xi(2,2) = chi_KS K_perp', xi(2, 2), cmplx(0.8_rp, 0.0_rp, rp))

      ! diagnose is the default and records every raw eigenvalue, its closest
      ! unity mode, m overlap, raw residual, and supplied bare gap estimate.
      call evaluate_goldstone(chi, provider, options, result, bare_spectral_gap=0.10_rp)
      call assert_true('diagnose is the development default', result%raw%available)
      call assert_true('all raw Xi eigenvalues retained', size(result%raw%eigenvalues) == 2)
      call assert_complex('raw unity-nearest eigenvalue', result%raw%closest_eigenvalue, cmplx(0.9_rp, 0.0_rp, rp))
      call assert_real('unity-mode overlap with projected m', result%raw%magnetization_overlap, 2.0_rp/sqrt(5.0_rp))
      call assert_real('raw Goldstone residual', result%raw%residual, sqrt(0.08_rp/5.0_rp))
      call assert_true('bare spectral gap recorded', result%raw%has_bare_spectral_gap)
      call assert_real('bare spectral gap value', result%raw%bare_spectral_gap, 0.10_rp)
   end subroutine test_xc_provider_kernel_and_raw_diagnostics

   subroutine test_off_mode_only_constructs_xi()
      type(xc_response_kernel_provider) :: provider
      type(tddft_goldstone_options) :: options
      type(tddft_goldstone_result) :: result
      complex(rp) :: chi(1, 1)

      call make_provider(provider, [1.0_rp], [0.5_rp])
      chi(1, 1) = cmplx(0.5_rp, 0.0_rp, rp)
      options%goldstone_mode = 'off'
      call evaluate_goldstone(chi, provider, options, result)
      call assert_complex('off mode still constructs Xi', result%xi_raw(1, 1), cmplx(0.25_rp, 0.0_rp, rp))
      call assert_true('off mode does not calculate Goldstone diagnostics', .not. result%raw%available)
      call assert_true('off mode does not apply a sum rule', .not. result%sum_rule_applied)
   end subroutine test_off_mode_only_constructs_xi

   subroutine test_controlled_pair_column_correction()
      type(tddft_goldstone_column_correction) :: correction
      complex(rp) :: xi(2, 2), corrected(2, 2)
      complex(rp) :: moment(2)

      xi = cmplx(0.0_rp, 0.0_rp, rp)
      xi(1, 1) = 0.9_rp; xi(2, 2) = 0.95_rp
      moment = [cmplx(2.0_rp, 0.0_rp, rp), cmplx(1.0_rp, 0.0_rp, rp)]
      call build_goldstone_column_correction(xi, moment, correction)
      call assert_true('real static pair correction is accepted', correction%applied .and. .not. correction%rejected)
      call assert_real_vector('one-site and multisite scales use the real static Xi', correction%scales, &
         [1.0_rp/0.9_rp, 1.0_rp/0.95_rp])
      call assert_real('corrected pair Xi restores Ward residual', correction%corrected%residual, 0.0_rp)
      call rescale_xi_columns(xi, correction%scales, corrected)
      call assert_true('correct mode changes dynamic Xi columns', maxval(abs(corrected-xi)) > 1.0e-6_rp)
      call assert_real('raw pair Xi is retained', correction%raw%residual, sqrt((0.2_rp**2+0.05_rp**2)/5.0_rp))
   end subroutine test_controlled_pair_column_correction

   subroutine test_correction_rejections_and_spectral_control()
      type(tddft_goldstone_column_correction) :: correction
      complex(rp) :: xi(2, 2)
      complex(rp) :: moment(2)

      moment = [cmplx(1.0_rp, 0.0_rp, rp), cmplx(1.0_rp, 0.0_rp, rp)]
      xi = cmplx(0.0_rp, 0.0_rp, rp)
      xi(:, 1) = [cmplx(0.5_rp, 0.0_rp, rp), cmplx(1.0_rp, 0.0_rp, rp)]
      xi(:, 2) = [cmplx(0.5_rp, 0.0_rp, rp), cmplx(1.0_rp, 0.0_rp, rp)]
      call build_goldstone_column_correction(xi, moment, correction)
      call assert_true('rank-deficient static correction is rejected', correction%rejected .and. .not. correction%applied)

      xi = cmplx(0.0_rp, 0.0_rp, rp); xi(1, 1) = 1.0_rp; xi(2, 2) = 1.0e-9_rp
      call build_goldstone_column_correction(xi, moment, correction)
      call assert_true('ill-conditioned static correction is rejected', correction%rejected)

      xi = cmplx(0.0_rp, 0.0_rp, rp); xi(1, 1) = 1.0_rp; xi(2, 2) = 1.0_rp
      moment(2) = cmplx(1.0e-14_rp, 0.0_rp, rp)
      call build_goldstone_column_correction(xi, moment, correction)
      call assert_true('small response moment correction is rejected', correction%rejected)

      moment = [cmplx(1.0_rp, 0.0_rp, rp), cmplx(1.0_rp, 0.0_rp, rp)]
      xi = cmplx(0.0_rp, 0.0_rp, rp); xi(1, 1) = cmplx(0.95_rp, 1.0e-4_rp, rp); xi(2, 2) = 0.95_rp
      call build_goldstone_column_correction(xi, moment, correction)
      call assert_true('complex static correction is rejected', correction%rejected)
      call assert_true('negative spectral-weight control is detected', .not. spectral_weights_are_nonnegative([-1.0e-4_rp]))
      call assert_true('nonnegative spectral-weight control passes', spectral_weights_are_nonnegative([0.0_rp, 1.0e-8_rp]))
      call assert_true('unchanged finite-eta circular tail is accepted', &
         spectral_weight_correction_is_acceptable([-0.20_rp, 0.10_rp], [-0.20_rp, 0.10_rp]))
      call assert_true('correction-created negative weight is rejected', &
         .not. spectral_weight_correction_is_acceptable([0.10_rp, 0.10_rp], [0.10_rp, -0.01_rp]))
      call assert_true('correction-worsened negative weight is rejected', &
         .not. spectral_weight_correction_is_acceptable([-0.20_rp, 0.10_rp], [-0.21_rp, 0.10_rp]))
      call assert_true('roundoff-sized change near a large finite-eta pole is accepted', &
         spectral_weight_correction_is_acceptable([-100.0_rp], [-100.0_rp-5.0e-9_rp]))
      call assert_true('resolved worsening near a large finite-eta pole is rejected', &
         .not. spectral_weight_correction_is_acceptable([-100.0_rp], [-100.0_rp-2.0e-8_rp]))
   end subroutine test_correction_rejections_and_spectral_control

   subroutine test_symmetry_breaking_disables_sum_rule()
      type(xc_response_kernel_provider) :: provider
      type(tddft_goldstone_options) :: options
      type(tddft_goldstone_result) :: result
      complex(rp) :: chi(1, 1)

      call make_provider(provider, [1.0_rp], [0.5_rp])
      chi(1, 1) = cmplx(0.5_rp, 0.0_rp, rp)
      options%goldstone_mode = 'correct'
      options%has_soc = .true.
      call evaluate_goldstone(chi, provider, options, result)
      call assert_true('SOC disables zero-gap enforcement', result%sum_rule_disabled_by_symmetry_breaking)
      call assert_true('SOC does not apply a legacy correction', .not. result%sum_rule_applied)
      call assert_true('SOC still reports raw diagnostics', result%raw%available)

      options%has_soc = .false.
      options%has_external_field = .true.
      call evaluate_goldstone(chi, provider, options, result)
      call assert_true('external field disables zero-gap enforcement', result%sum_rule_disabled_by_symmetry_breaking)
   end subroutine test_symmetry_breaking_disables_sum_rule

   subroutine test_raw_residual_convergence_controls()
      ! These fixtures deliberately keep the physical static chi_KS fixed while
      ! changing response representation, a flat-band k mesh, and an explicit
      ! converged band window.  The raw residual must be invariant; a real
      ! calculation uses the same three controls as convergence axes.
      type(xc_response_kernel_provider) :: one_site, two_site
      type(tddft_goldstone_options) :: options
      type(tddft_goldstone_result) :: one_basis, two_basis, mesh4, mesh8, bands2, bands4
      complex(rp) :: chi1(1, 1), chi2(2, 2)

      call make_provider(one_site, [1.0_rp], [0.8_rp])
      chi1(1, 1) = cmplx(0.4_rp, 0.0_rp, rp)
      call evaluate_goldstone(chi1, one_site, options, one_basis)

      call make_provider(two_site, [1.0_rp, 1.0_rp], [0.8_rp, 0.8_rp])
      chi2 = cmplx(0.0_rp, 0.0_rp, rp)
      chi2(1, 1) = chi1(1, 1)
      chi2(2, 2) = chi1(1, 1)
      call evaluate_goldstone(chi2, two_site, options, two_basis)
      call assert_real('raw residual response-basis convergence', two_basis%raw%residual, one_basis%raw%residual)

      ! Flat-band 4 and 8 point meshes have the same exact BZ average.
      call evaluate_goldstone(chi1, one_site, options, mesh4)
      call evaluate_goldstone(chi1, one_site, options, mesh8)
      call assert_real('raw residual k-mesh convergence', mesh8%raw%residual, mesh4%raw%residual)

      ! The two-band and four-band fixtures here have already converged to the
      ! same active spin-flip response; the residual is intentionally raw.
      call evaluate_goldstone(chi1, one_site, options, bands2)
      call evaluate_goldstone(chi1, one_site, options, bands4)
      call assert_real('raw residual band-window convergence', bands4%raw%residual, bands2%raw%residual)
   end subroutine test_raw_residual_convergence_controls

   subroutine test_signed_magnetization_diagnostics()
      complex(rp) :: xi_one(1, 1), xi_two(2, 2)
      type(tddft_goldstone_diagnostics) :: reversed, ferro, antiferro

      xi_one(1, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      call evaluate_raw_xi_diagnostics(xi_one, [cmplx(-2.0_rp, 0.0_rp, rp)], reversed)
      call assert_real('reversed one-site signed moment is a zero mode', reversed%residual, 0.0_rp)
      call assert_real('reversed one-site right overlap is retained', reversed%magnetization_overlap, 1.0_rp)
      call assert_real('real static Xi has zero imaginary norm', reversed%imaginary_norm, 0.0_rp)

      xi_two = cmplx(0.0_rp, 0.0_rp, rp)
      xi_two(1, 1) = 1.0_rp; xi_two(2, 2) = 1.0_rp
      call evaluate_raw_xi_diagnostics(xi_two, [cmplx(1.0_rp, 0.0_rp, rp), cmplx(1.0_rp, 0.0_rp, rp)], ferro)
      call evaluate_raw_xi_diagnostics(xi_two, [cmplx(1.0_rp, 0.0_rp, rp), cmplx(-1.0_rp, 0.0_rp, rp)], antiferro)
      call assert_real('two-site ferro signed moment is a zero mode', ferro%residual, 0.0_rp)
      call assert_real('two-site antiferro signed moment is a zero mode', antiferro%residual, 0.0_rp)
   end subroutine test_signed_magnetization_diagnostics

   subroutine test_explicit_goldstone_policies()
      type(xc_response_kernel_provider) :: provider
      type(tddft_goldstone_options) :: options
      type(tddft_goldstone_result) :: result
      complex(rp) :: chi(2, 2)

      chi = cmplx(0.0_rp, 0.0_rp, rp)
      chi(1, 1) = 0.5_rp
      chi(2, 2) = 0.25_rp
      call make_provider(provider, [2.0_rp, 1.0_rp], [1.8_rp, 4.0_rp])
      options%goldstone_policy = 'sum_rule'
      call evaluate_goldstone(chi, provider, options, result)
      call assert_true('explicit Lounis policy is applied by Goldstone driver', result%lounis%applied)
      call assert_real('Goldstone driver Lounis kernel', real(result%k_perp_sum_rule(1), rp), 2.0_rp)
      call assert_real('Goldstone driver Lounis residual', result%lounis%corrected%ward_residual, 0.0_rp)

      chi = cmplx(0.0_rp, 0.0_rp, rp)
      chi(1, 1) = 1.0_rp
      chi(2, 2) = 1.0_rp
      call make_provider(provider, [1.0_rp, 0.0_rp], [0.9_rp, 0.2_rp])
      options%goldstone_policy = 'projected'
      call evaluate_goldstone(chi, provider, options, result)
      call assert_true('explicit Halle policy is applied by Goldstone driver', result%projection%applied)
      call assert_real('Goldstone driver projected kernel', real(result%kernel_corrected(1, 1), rp), 1.0_rp)
   end subroutine test_explicit_goldstone_policies

   subroutine test_independent_bxc_source_and_ward()
      type(xc_response_kernel_provider) :: provider
      type(tddft_goldstone_options) :: options
      type(tddft_goldstone_result) :: result
      complex(rp), allocatable :: source(:)
      complex(rp) :: chi(2, 2)
      real(rp) :: amplitudes(2), signed_magnetization(2)
      character(len=160) :: provenance

      call provider%initialize(2, 'unit independent VXC0SP')
      amplitudes = [2.0_rp, 1.0_rp]
      signed_magnetization = [2.0_rp, -1.0_rp]
      provider%site(1)%spin_population = amplitudes(1)
      provider%site(2)%spin_population = amplitudes(2)
      provider%site(1)%signed_spin_population = signed_magnetization(1)
      provider%site(2)%signed_spin_population = signed_magnetization(2)
      provider%site(:)%has_signed_spin_population = .true.
      provider%site(:)%has_radial_projection = .true.
      provider%site(1)%bxc_spin_moment = 8.0_rp
      provider%site(2)%bxc_spin_moment = 4.0_rp
      provider%site(1)%k_perp_circular = 1.0_rp
      provider%site(2)%k_perp_circular = 2.0_rp
      provider%site(:)%has_k_perp_circular = .true.
      call build_independent_circular_bxc_source(provider, amplitudes, signed_magnetization, source, provenance)
      call assert_real('independent +z circular source', real(source(1), rp), 2.0_rp)
      call assert_real('independent -z circular source carries sign', real(source(2), rp), -2.0_rp)
      call assert_true('independent source provenance is explicit', index(provenance, 'VXC0SP') > 0)

      chi = cmplx(0.0_rp, 0.0_rp, rp)
      chi(1, 1) = 1.0_rp
      chi(2, 2) = 0.5_rp
      options%require_independent_ward = .true.
      call evaluate_goldstone(chi, provider, options, result)
      call assert_true('production-style independent Ward record is available', result%raw%ward%independent_bxc)
      call assert_real('independent r_B is zero for the direct source', result%raw%ward%r_b, 0.0_rp)
      call assert_true('signed multisublattice Goldstone vector is retained', &
         real(result%raw%signed_magnetization(1), rp) > 0.0_rp .and. real(result%raw%signed_magnetization(2), rp) < 0.0_rp)
   end subroutine test_independent_bxc_source_and_ward

   subroutine test_anti_goldstone_is_phase_aware()
      type(tddft_goldstone_diagnostics) :: positive, negative
      complex(rp) :: xi_positive(1, 1), xi_negative(1, 1), m(1)

      xi_positive(1, 1) = 1.0_rp
      xi_negative(1, 1) = -1.0_rp
      m(1) = 2.0_rp
      call evaluate_raw_xi_diagnostics(xi_positive, m, positive)
      call evaluate_raw_xi_diagnostics(xi_negative, m, negative)
      call assert_real('positive Goldstone magnitude overlap', positive%normalized_magnitude_overlap, 1.0_rp)
      call assert_real('anti-Goldstone magnitude overlap alone is still one', negative%normalized_magnitude_overlap, 1.0_rp)
      call assert_complex('positive Goldstone phase-sensitive action', positive%phase_sensitive_overlap, cmplx(1.0_rp, 0.0_rp, rp))
      call assert_complex('anti-Goldstone phase-sensitive action is negative', negative%phase_sensitive_overlap, &
         cmplx(-1.0_rp, 0.0_rp, rp))
      call assert_complex('anti-Goldstone complex distance to +1', negative%closest_eigenvalue_delta, cmplx(-2.0_rp, 0.0_rp, rp))
      call assert_true('anti-Goldstone branch is not identified as +1', negative%closest_eigenvalue_distance > 1.0_rp)
   end subroutine test_anti_goldstone_is_phase_aware

   subroutine test_compare_branch_policy()
      type(tddft_xi_compare_diagnostics) :: comparison
      complex(rp) :: legacy(2, 2), pair(2, 2), anti(2, 2), m(2)

      legacy = cmplx(0.0_rp, 0.0_rp, rp)
      pair = legacy
      anti = legacy
      legacy(1, 1) = 0.98_rp; legacy(2, 2) = 0.97_rp
      pair(1, 1) = 1.0_rp; pair(2, 2) = 1.0_rp
      anti(1, 1) = -1.0_rp; anti(2, 2) = -1.0_rp
      m = [cmplx(1.0_rp, 0.0_rp, rp), cmplx(1.0_rp, 0.0_rp, rp)]
      call compare_xi_goldstone_branches(legacy, pair, m, comparison)
      call assert_true('compare policy accepts matching +1 branches', comparison%plus_one_branch_identified .and. &
         .not. comparison%hard_branch_failure)
      call compare_xi_goldstone_branches(legacy, anti, m, comparison)
      call assert_true('compare policy hard-fails opposite Goldstone branch', comparison%hard_branch_failure .and. &
         comparison%opposite_branch_detected)
      call assert_true('compare policy does not require magnitude equivalence by default', .not. comparison%strict_magnitude_gate)
   end subroutine test_compare_branch_policy

   subroutine test_diagnostic_output_contract()
      type(xc_response_kernel_provider) :: provider
      type(tddft_goldstone_options) :: options
      type(tddft_goldstone_result) :: result
      complex(rp) :: chi(1, 1)
      integer :: unit, ios
      logical :: found_identity, found_malformed_identity, found_rank, found_signed, found_provenance, found_r_b, found_alias
      character(len=512) :: line

      call provider%initialize(1, 'unit diagnostic output XC')
      call provider%set_site_spin_population(1, 2.0_rp)
      call provider%set_site_signed_spin_population(1, 2.0_rp)
      provider%site(1)%has_radial_projection = .true.
      provider%site(1)%bxc_spin_moment = 4.0_rp
      call provider%set_site_derivatives(1, k_perp_circular=0.5_rp)
      chi(1, 1) = 0.5_rp
      call evaluate_goldstone(chi, provider, options, result)
      call write_goldstone_diagnostics_text('unit_tddft_goldstone_diagnostics.dat', result)
      open(newunit=unit, file='unit_tddft_goldstone_diagnostics.dat', status='old', action='read', iostat=ios)
      found_identity = .false.; found_malformed_identity = .false.; found_rank = .false.; found_signed = .false.
      found_provenance = .false.; found_r_b = .false.; found_alias = .false.
      do while (ios == 0)
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         found_identity = found_identity .or. index(line, '# raw_identity_consistent = ') == 1
         found_malformed_identity = found_malformed_identity .or. index(line, 'raw_identity_consistentF') > 0
         found_rank = found_rank .or. index(line, '# response_space_rank = ') == 1
         found_signed = found_signed .or. index(line, 'raw_signed_magnetization') > 0
         found_provenance = found_provenance .or. index(line, 'VXC0SP') > 0
         found_r_b = found_r_b .or. index(line, '# raw_r_B = ') == 1
         found_alias = found_alias .or. index(line, 'ward_residual') > 0 .or. index(line, 'dm_residual') > 0 .or. &
            index(line, 'bxc_kernel_residual') > 0
      end do
      close(unit, status='delete')
      call assert_true('raw_identity_consistent output has a separated value', found_identity)
      call assert_true('raw_identity_consistent output is not concatenated with its value', .not. found_malformed_identity)
      call assert_true('response-space rank is present in the header', found_rank)
      call assert_true('signed site magnetization is written', found_signed)
      call assert_true('independent Bxc provenance is written', found_provenance)
      call assert_true('independent r_B is written', found_r_b)
      call assert_true('tautological residual aliases are absent from output', .not. found_alias)
   end subroutine test_diagnostic_output_contract

   subroutine make_provider(provider, moment, kernel)
      type(xc_response_kernel_provider), intent(out) :: provider
      real(rp), intent(in) :: moment(:), kernel(:)
      integer :: i

      call provider%initialize(size(moment), 'unit-test-XC')
      do i = 1, size(moment)
         call provider%set_site_derivatives(i, k_perp_circular=kernel(i))
         provider%site(i)%spin_population = moment(i)
      end do
   end subroutine make_provider

   subroutine assert_complex(label, actual, expected)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual, expected
      if (abs(actual - expected) > tol*max(1.0_rp, abs(expected))) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, abs(actual - expected)
         failed = .true.
      end if
   end subroutine assert_complex

   subroutine assert_real(label, actual, expected)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: actual, expected
      if (abs(actual - expected) > tol*max(1.0_rp, abs(expected))) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, abs(actual - expected)
         failed = .true.
      end if
   end subroutine assert_real

   subroutine assert_real_vector(label, actual, expected)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: actual(:), expected(:)
      call assert_real(label, maxval(abs(actual - expected)), 0.0_rp)
   end subroutine assert_real_vector

   subroutine assert_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine assert_true

end program test_tddft_goldstone
