!------------------------------------------------------------------------------
! TDDFT-04 -- transverse ALSDA kernel and transparent Goldstone diagnostics
!------------------------------------------------------------------------------
program test_tddft_goldstone
   use precision_mod, only: rp
   use xc_response_kernel_mod, only: xc_response_kernel_provider
   use tddft_goldstone_mod, only: tddft_goldstone_options, tddft_goldstone_result, &
      tddft_goldstone_diagnostics, build_site_projected_k_perp, construct_transverse_xi, evaluate_goldstone, &
      evaluate_raw_xi_diagnostics, write_goldstone_diagnostics_text
   implicit none

   real(rp), parameter :: tol = 2.0e-12_rp
   logical :: failed

   failed = .false.
   call test_xc_provider_kernel_and_raw_diagnostics()
   call test_off_mode_only_constructs_xi()
   call test_sum_rule_repair_preserves_raw_diagnostics()
   call test_symmetry_breaking_disables_sum_rule()
   call test_raw_residual_convergence_controls()
   call test_signed_magnetization_diagnostics()

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

   subroutine test_sum_rule_repair_preserves_raw_diagnostics()
      type(xc_response_kernel_provider) :: provider
      type(tddft_goldstone_options) :: options
      type(tddft_goldstone_result) :: result
      complex(rp) :: chi(2, 2)

      call make_provider(provider, [2.0_rp, 1.0_rp], [3.0_rp, 4.0_rp])
      chi = cmplx(0.0_rp, 0.0_rp, rp)
      chi(1, 1) = cmplx(0.3_rp, 0.0_rp, rp)
      chi(2, 2) = cmplx(0.2_rp, 0.0_rp, rp)
      options%goldstone_mode = 'sum_rule'
      call evaluate_goldstone(chi, provider, options, result)

      call assert_true('perturbed kernel has an exposed raw violation', result%raw%residual > 0.1_rp)
      call assert_true('sum-rule correction is marked as applied', result%sum_rule_applied)
      call assert_true('raw diagnostics survive correction', result%raw%available .and. size(result%raw%eigenvalues) == 2)
      call assert_real('raw residual was not overwritten', result%raw%residual, sqrt(0.08_rp/5.0_rp))
      call assert_real('sum-rule correction restores Goldstone condition', result%corrected%residual, 0.0_rp)
      call assert_real('corrected unity eigenvalue distance', result%corrected%closest_eigenvalue_distance, 0.0_rp)
      call assert_real_vector('documented static K_perp sum-rule correction', real(result%k_perp_sum_rule, rp), &
         [1.0_rp/0.3_rp, 5.0_rp])
      call write_goldstone_diagnostics_text('unit_goldstone_diagnostics.out', result)
      call assert_output_has_raw_and_corrected_records('unit_goldstone_diagnostics.out')
   end subroutine test_sum_rule_repair_preserves_raw_diagnostics

   subroutine test_symmetry_breaking_disables_sum_rule()
      type(xc_response_kernel_provider) :: provider
      type(tddft_goldstone_options) :: options
      type(tddft_goldstone_result) :: result
      complex(rp) :: chi(1, 1)

      call make_provider(provider, [1.0_rp], [0.5_rp])
      chi(1, 1) = cmplx(0.5_rp, 0.0_rp, rp)
      options%goldstone_mode = 'sum_rule'
      options%has_soc = .true.
      call evaluate_goldstone(chi, provider, options, result)
      call assert_true('SOC disables zero-gap enforcement', result%sum_rule_disabled_by_symmetry_breaking)
      call assert_true('SOC does not apply sum rule', .not. result%sum_rule_applied)
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

   subroutine make_provider(provider, moment, kernel)
      type(xc_response_kernel_provider), intent(out) :: provider
      real(rp), intent(in) :: moment(:), kernel(:)
      integer :: i

      call provider%initialize(size(moment), 'unit-test-XC')
      do i = 1, size(moment)
         call provider%set_site_derivatives(i, k_perp=kernel(i))
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

   subroutine assert_output_has_raw_and_corrected_records(filename)
      character(len=*), intent(in) :: filename
      character(len=256) :: line
      logical :: saw_raw_value, saw_raw_vector, saw_corrected_value
      integer :: unit, ios

      saw_raw_value = .false.
      saw_raw_vector = .false.
      saw_corrected_value = .false.
      open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call assert_true('Goldstone diagnostic output was written', .false.)
         return
      end if
      do
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, 'raw_eigenvalue') > 0) saw_raw_value = .true.
         if (index(line, 'raw_closest_eigenvector') > 0) saw_raw_vector = .true.
         if (index(line, 'sum_rule_corrected_eigenvalue') > 0) saw_corrected_value = .true.
      end do
      close(unit, status='delete')
      call assert_true('output preserves raw eigenvalues', saw_raw_value)
      call assert_true('output preserves raw unity eigenvector', saw_raw_vector)
      call assert_true('output includes corrected diagnostics separately', saw_corrected_value)
   end subroutine assert_output_has_raw_and_corrected_records

end program test_tddft_goldstone
