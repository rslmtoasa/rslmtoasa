!------------------------------------------------------------------------------
! TDDFT-08 -- symmetric finite-field longitudinal calibration and relaxation.
!------------------------------------------------------------------------------
program test_tddft_longitudinal
   use precision_mod, only: rp
   use tddft_longitudinal_mod, only: tddft_longitudinal_options, tddft_longitudinal_static_result, &
      tddft_longitudinal_result, build_longitudinal_static_response, calibrate_longitudinal_response
   implicit none

   logical :: failed

   failed = .false.
   call test_symmetric_static_calibration_and_fit()
   if (failed) error stop 1
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_symmetric_static_calibration_and_fit()
      integer, parameter :: nsite = 2, nrecords = 8, nw = 5
      integer :: source(nrecords), ir, i, j, iw
      real(rp) :: field(nrecords), moments(nsite, nrecords), m0(nsite), omega(nw)
      real(rp) :: chi_static_real(nsite, nsite), tau(nsite)
      complex(rp) :: chi_ks(nsite, nsite), chi_static(nsite, nsite), chi_dynamic(nsite, nsite, nw)
      complex(rp) :: expected_u(nsite, nsite)
      type(tddft_longitudinal_options) :: options
      type(tddft_longitudinal_static_result) :: static_result
      type(tddft_longitudinal_result) :: result

      ! For each source site, two +/- amplitudes are supplied.  The exact
      ! linear slopes form a non-diagonal static site-site susceptibility.
      chi_static_real = reshape([4.0_rp, 0.3_rp, 0.3_rp, 3.0_rp], [nsite, nsite])
      source = [1, 1, 1, 1, 2, 2, 2, 2]
      field = [0.010_rp, -0.010_rp, 0.005_rp, -0.005_rp, 0.010_rp, -0.010_rp, 0.005_rp, -0.005_rp]
      do ir = 1, nrecords
         moments(:, ir) = matmul(chi_static_real, unit_source(source(ir), nsite))*field(ir)
      end do
      options%linearity_tolerance = 1.0e-12_rp
      options%static_agreement_tolerance = 1.0e-12_rp
      options%fit_omega_max = 0.040_rp
      call build_longitudinal_static_response(source, field, moments, options, static_result)
      call check_real_matrix('symmetric static chi', static_result%chi, chi_static_real)
      call check_true('linearity check passes', static_result%linearity_passed)

      m0 = [2.2_rp, 1.7_rp]
      chi_ks = reshape([2.0_rp, 0.1_rp, 0.1_rp, 1.5_rp], [nsite, nsite])
      chi_static = cmplx(chi_static_real, 0.0_rp, rp)
      expected_u = inverse_2x2(chi_ks) - inverse_2x2(chi_static)
      omega = [0.0_rp, 0.010_rp, 0.020_rp, 0.030_rp, 0.060_rp]
      tau = [12.0_rp, 7.0_rp]
      do iw = 1, nw
         do j = 1, nsite
            do i = 1, nsite
               chi_dynamic(i, j, iw) = chi_static(i, j)/cmplx(1.0_rp, omega(iw)*sqrt(tau(i)*tau(j)), rp)
            end do
         end do
      end do
      call calibrate_longitudinal_response(m0, chi_ks, chi_static, omega, chi_dynamic, options, result)
      call check_complex_matrix('inverse-susceptibility U', result%u_parallel, expected_u)
      call check_real_vector('relaxation times', result%t_parallel, tau)
      call check_real_vector('inverse relaxation rates', result%gamma_parallel, 1.0_rp/tau)
      call check_true('dynamic/static acceptance passes', result%static_agreement_passed)
      call check_true('fit excludes high frequency', result%fit_first == 2 .and. result%fit_last == 4)
      call check_true('zero residual synthetic relaxation', maxval(result%fit_residual) < 1.0e-12_rp)
   end subroutine test_symmetric_static_calibration_and_fit

   pure function unit_source(isource, n) result(vector)
      integer, intent(in) :: isource, n
      real(rp) :: vector(n)
      vector = 0.0_rp; vector(isource) = 1.0_rp
   end function unit_source

   pure function inverse_2x2(matrix) result(inverse)
      complex(rp), intent(in) :: matrix(2, 2)
      complex(rp) :: inverse(2, 2), determinant
      determinant = matrix(1, 1)*matrix(2, 2) - matrix(1, 2)*matrix(2, 1)
      inverse = reshape([matrix(2, 2), -matrix(2, 1), -matrix(1, 2), matrix(1, 1)], [2, 2])/determinant
   end function inverse_2x2

   subroutine check_real_matrix(label, actual, expected)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: actual(:, :), expected(:, :)
      if (maxval(abs(actual-expected)) > 1.0e-11_rp) then
         write (*, '(a,1x,a)') 'FAIL', label; failed = .true.
      end if
   end subroutine check_real_matrix

   subroutine check_complex_matrix(label, actual, expected)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual(:, :), expected(:, :)
      if (maxval(abs(actual-expected)) > 1.0e-11_rp) then
         write (*, '(a,1x,a)') 'FAIL', label; failed = .true.
      end if
   end subroutine check_complex_matrix

   subroutine check_real_vector(label, actual, expected)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: actual(:), expected(:)
      if (maxval(abs(actual-expected)) > 1.0e-11_rp) then
         write (*, '(a,1x,a)') 'FAIL', label; failed = .true.
      end if
   end subroutine check_real_vector

   subroutine check_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', label; failed = .true.
      end if
   end subroutine check_true

end program test_tddft_longitudinal
