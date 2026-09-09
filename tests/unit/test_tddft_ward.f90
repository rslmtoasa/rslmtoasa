!------------------------------------------------------------------------------
! TDDFT-02 -- raw Ward identity and explicit repair operators
!------------------------------------------------------------------------------
program test_tddft_ward
   use precision_mod, only: rp
   use tddft_ward_mod, only: tddft_ward_diagnostics, tddft_lounis_repair, tddft_goldstone_projection, &
      tddft_static_dynamic_diagnostics, evaluate_static_ward_identity, evaluate_ward_from_xi, &
      evaluate_static_dynamic_consistency, reconstruct_lounis_kernel, project_goldstone_eigenvalue
   implicit none

   real(rp), parameter :: tol = 2.0e-11_rp
   logical :: failed

   failed = .false.
   call test_raw_ward_identity()
   call test_lounis_reconstruction()
   call test_halle_projection()
   call test_xi_only_provenance()
   call test_independent_residuals_fail_separately()
   call test_static_dynamic_eta_ladder()
   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_raw_ward_identity()
      type(tddft_ward_diagnostics) :: diagnostics
      complex(rp) :: chi(2, 2), kernel(2, 2), bxc(2), magnetization(2)

      chi = cmplx(0.0_rp, 0.0_rp, rp)
      chi(1, 1) = 0.5_rp
      chi(2, 2) = 0.25_rp
      kernel = cmplx(0.0_rp, 0.0_rp, rp)
      kernel(1, 1) = 2.0_rp
      kernel(2, 2) = 4.0_rp
      magnetization = [cmplx(2.0_rp, 0.0_rp, rp), cmplx(-1.0_rp, 0.0_rp, rp)]
      bxc = matmul(kernel, magnetization)
      call evaluate_static_ward_identity(chi, bxc, magnetization, diagnostics, kernel=kernel, response_basis='site', &
         bxc_provenance='unit ground-state XC field', kernel_provenance='unit physical kernel')
      call assert_true('raw Ward record is available', diagnostics%available)
      call assert_real('exact chi KS Bxc-m residual', diagnostics%ward_residual, 0.0_rp)
      call assert_real('exact Dm residual', diagnostics%dm_residual, 0.0_rp)
      call assert_true('Ward basis provenance is retained', trim(diagnostics%response_basis) == 'site')
      call assert_true('Bxc provenance is retained', index(diagnostics%bxc_provenance, 'ground-state') > 0)
   end subroutine test_raw_ward_identity

   subroutine test_lounis_reconstruction()
      type(tddft_lounis_repair) :: repair
      type(tddft_ward_diagnostics) :: raw
      complex(rp) :: chi(2, 2), physical(2), m(2)
      complex(rp), allocatable :: reconstructed(:, :)

      chi = cmplx(0.0_rp, 0.0_rp, rp)
      chi(1, 1) = 0.5_rp
      chi(2, 2) = 0.25_rp
      m = [cmplx(2.0_rp, 0.0_rp, rp), cmplx(-1.0_rp, 0.0_rp, rp)]
      physical = [cmplx(1.8_rp, 0.0_rp, rp), cmplx(4.0_rp, 0.0_rp, rp)]
      call evaluate_static_ward_identity(chi, physical*m, m, raw, response_basis='site', &
         bxc_provenance='unit physical field', kernel_provenance='unit physical kernel')
      call reconstruct_lounis_kernel(chi, m, reconstructed, repair, physical_kernel=physical, warning_threshold=0.01_rp)
      call assert_true('Lounis reconstruction is applied explicitly', repair%applied .and. .not. repair%rejected)
      call assert_real('Lounis first diagonal kernel', real(reconstructed(1, 1), rp), 2.0_rp)
      call assert_real('Lounis second diagonal kernel', real(reconstructed(2, 2), rp), 4.0_rp)
      call assert_real('Lounis corrected Ward residual', repair%corrected%ward_residual, 0.0_rp)
      call assert_true('Lounis material-change warning is recorded', repair%large_correction)
      call assert_true('Lounis raw residual is retained', repair%raw%available .and. repair%raw%ward_residual > 0.0_rp)
   end subroutine test_lounis_reconstruction

   subroutine test_halle_projection()
      type(tddft_goldstone_projection) :: repair
      complex(rp) :: xi(2, 2), m(2)
      complex(rp), allocatable :: projected(:, :)

      xi = cmplx(0.0_rp, 0.0_rp, rp)
      xi(1, 1) = 0.9_rp
      xi(2, 2) = 0.2_rp
      m = [cmplx(1.0_rp, 0.0_rp, rp), cmplx(0.0_rp, 0.0_rp, rp)]
      call project_goldstone_eigenvalue(xi, m, projected, repair)
      call assert_true('Halle projection is applied explicitly', repair%applied .and. .not. repair%rejected)
      call assert_real('projected Goldstone Ward residual', repair%corrected%ward_residual, 0.0_rp)
      call assert_real('projected Goldstone eigenvalue', real(projected(1, 1), rp), 1.0_rp)
      call assert_real('non-Goldstone eigenvalue is unchanged', real(projected(2, 2), rp), 0.2_rp)
      call assert_true('Halle raw residual is retained', repair%raw%ward_residual > 0.0_rp)
   end subroutine test_halle_projection

   subroutine test_xi_only_provenance()
      type(tddft_ward_diagnostics) :: diagnostics
      complex(rp) :: xi(1, 1), m(1)

      xi(1, 1) = 0.75_rp
      m(1) = 2.0_rp
      call evaluate_ward_from_xi(xi, m, diagnostics, response_basis='site-orbital', &
         kernel_provenance='direct LMTO pair potential')
      call assert_real('Xi-only Ward residual', diagnostics%ward_residual, 0.25_rp)
      call assert_real('Xi-only Dm residual', diagnostics%dm_residual, 0.25_rp)
      call assert_true('Xi-only operator provenance is retained', index(diagnostics%kernel_provenance, 'LMTO') > 0)
      call assert_real('Xi-only canonical derived residual', diagnostics%derived_identity_residual, 0.25_rp)
   end subroutine test_xi_only_provenance

   subroutine test_independent_residuals_fail_separately()
      type(tddft_ward_diagnostics) :: direct_good, direct_bad, xi_good, xi_bad
      complex(rp) :: chi(2, 2), m(2), bxc_good(2), bxc_bad(2), xi_good_matrix(2, 2), xi_bad_matrix(2, 2)

      chi = cmplx(0.0_rp, 0.0_rp, rp)
      chi(1, 1) = 0.5_rp
      chi(2, 2) = 0.25_rp
      m = [cmplx(2.0_rp, 0.0_rp, rp), cmplx(-1.0_rp, 0.0_rp, rp)]
      bxc_good = [cmplx(4.0_rp, 0.0_rp, rp), cmplx(-4.0_rp, 0.0_rp, rp)]
      bxc_bad = bxc_good
      bxc_bad(1) = bxc_bad(1) + cmplx(0.4_rp, 0.0_rp, rp)
      call evaluate_static_ward_identity(chi, bxc_good, m, direct_good, response_basis='site', &
         bxc_provenance='independent VXC0SP unit source', bxc_is_independent=.true.)
      call evaluate_static_ward_identity(chi, bxc_bad, m, direct_bad, response_basis='site', &
         bxc_provenance='perturbed independent VXC0SP unit source', bxc_is_independent=.true.)
      call assert_real('independent r_B passes before Bxc perturbation', direct_good%r_b, 0.0_rp)
      call assert_true('independent r_B fails under Bxc perturbation', direct_bad%r_b > 0.05_rp)
      call assert_true('independent provenance is not a kernel identity', direct_bad%independent_bxc .and. &
         .not. direct_bad%derived_identity)

      xi_good_matrix = cmplx(0.0_rp, 0.0_rp, rp)
      xi_good_matrix(1, 1) = 1.0_rp
      xi_good_matrix(2, 2) = 1.0_rp
      xi_bad_matrix = xi_good_matrix
      xi_bad_matrix(1, 1) = 0.8_rp
      call evaluate_ward_from_xi(xi_good_matrix, m, xi_good, response_basis='site', kernel_provenance='independent Xi unit source')
      call evaluate_ward_from_xi(xi_bad_matrix, m, xi_bad, response_basis='site', kernel_provenance='perturbed Xi unit source')
      call assert_real('r_Xi passes before Xi perturbation', xi_good%derived_identity_residual, 0.0_rp)
      call assert_true('r_Xi fails under Xi perturbation', xi_bad%derived_identity_residual > 0.05_rp)
      call assert_true('Bxc and Xi failures are separately observable', direct_bad%r_b > 0.05_rp .and. &
         xi_good%derived_identity_residual < tol)
   end subroutine test_independent_residuals_fail_separately

   subroutine test_static_dynamic_eta_ladder()
      type(tddft_static_dynamic_diagnostics) :: diagnostics
      complex(rp) :: chi_static(2, 2), chi_dynamic(2, 2, 3)
      real(rp) :: eta(3)

      chi_static = cmplx(0.0_rp, 0.0_rp, rp)
      chi_static(1, 1) = 2.0_rp
      chi_static(2, 2) = 1.0_rp
      chi_dynamic = spread(chi_static, 3, 3)
      chi_dynamic(1, 1, 1) = chi_dynamic(1, 1, 1) + cmplx(0.2_rp, 0.0_rp, rp)
      chi_dynamic(1, 1, 2) = chi_dynamic(1, 1, 2) + cmplx(0.02_rp, 0.0_rp, rp)
      chi_dynamic(1, 1, 3) = chi_dynamic(1, 1, 3) + cmplx(0.002_rp, 0.0_rp, rp)
      eta = [0.1_rp, 0.01_rp, 0.001_rp]
      call evaluate_static_dynamic_consistency(chi_static, chi_dynamic, eta, diagnostics, response_basis='site', &
         q_provenance='Gamma unit eta ladder')
      call assert_true('static/dynamic eta ladder is available', diagnostics%available .and. diagnostics%eta_ladder_valid)
      call assert_true('static/dynamic residual decreases along eta ladder', diagnostics%residual(1) > diagnostics%residual(2) .and. &
         diagnostics%residual(2) > diagnostics%residual(3))
      call assert_real('smallest eta residual is reported', diagnostics%smallest_eta_residual, diagnostics%residual(3))
      call assert_true('static/dynamic rank is reported', diagnostics%response_space_rank == 2)
   end subroutine test_static_dynamic_eta_ladder

   subroutine assert_real(label, actual, expected)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: actual, expected
      if (abs(actual-expected) > tol*max(1.0_rp, abs(expected))) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, abs(actual-expected)
         failed = .true.
      end if
   end subroutine assert_real

   subroutine assert_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine assert_true

end program test_tddft_ward
