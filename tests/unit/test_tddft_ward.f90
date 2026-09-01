!------------------------------------------------------------------------------
! TDDFT-02 -- raw Ward identity and explicit repair operators
!------------------------------------------------------------------------------
program test_tddft_ward
   use precision_mod, only: rp
   use tddft_ward_mod, only: tddft_ward_diagnostics, tddft_lounis_repair, tddft_goldstone_projection, &
      evaluate_static_ward_identity, evaluate_ward_from_xi, reconstruct_lounis_kernel, project_goldstone_eigenvalue
   implicit none

   real(rp), parameter :: tol = 2.0e-11_rp
   logical :: failed

   failed = .false.
   call test_raw_ward_identity()
   call test_lounis_reconstruction()
   call test_halle_projection()
   call test_xi_only_provenance()
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
   end subroutine test_xi_only_provenance

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
