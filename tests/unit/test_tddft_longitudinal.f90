!------------------------------------------------------------------------------
! TDDFT-13 -- coupled charge/longitudinal response, plus TDDFT-08 compatibility.
!------------------------------------------------------------------------------
program test_tddft_longitudinal
   use precision_mod, only: rp
   use response_components_mod, only: RESPONSE_CHARGE, RESPONSE_MZ
   use response_vertices_mod, only: response_channel
   use xc_response_kernel_mod, only: xc_response_kernel_provider, xc_response_radial_projection
   use tddft_longitudinal_mod, only: tddft_longitudinal_options, tddft_longitudinal_static_result, &
      tddft_longitudinal_result, build_longitudinal_static_response, calibrate_longitudinal_response, &
      longitudinal_index, build_charge_longitudinal_channels, build_charge_longitudinal_kernel, &
      build_charge_longitudinal_chi_ks
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result
   implicit none

   logical :: failed

   failed = .false.
   call test_coupled_basis_and_kernel()
   call test_coupled_eigenpair_adapter()
   call test_symmetric_static_calibration_and_fit()
   if (failed) error stop 1
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_coupled_basis_and_kernel()
      integer, parameter :: nsite = 2
      type(response_channel), allocatable :: channels(:)
      type(xc_response_kernel_provider) :: provider, projected_provider
      type(xc_response_radial_projection) :: projection
      complex(rp) :: kernel(2*nsite, 2*nsite)
      complex(rp) :: projected_kernel(2, 2)
      real(rp) :: coulomb(nsite, nsite)
      logical :: supported
      character(len=128) :: reason

      call build_charge_longitudinal_channels(nsite, channels)
      call check_true('charge index is site-major', longitudinal_index(1, RESPONSE_CHARGE) == 1 .and. &
         longitudinal_index(1, RESPONSE_MZ) == 2 .and. longitudinal_index(2, RESPONSE_CHARGE) == 3 .and. &
         longitudinal_index(2, RESPONSE_MZ) == 4)
      call check_true('coupled channel list has two entries per site', size(channels) == 2*nsite .and. &
         channels(1)%component == RESPONSE_CHARGE .and. channels(2)%component == RESPONSE_MZ .and. &
         channels(3)%component == RESPONSE_CHARGE .and. channels(4)%component == RESPONSE_MZ)

      call provider%initialize(nsite, 'synthetic-LDA')
      call provider%set_site_derivatives(1, 1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp)
      call provider%set_site_derivatives(2, 5.0_rp, 6.0_rp, 7.0_rp, 8.0_rp)
      coulomb = reshape([2.0_rp, 0.2_rp, 0.2_rp, 3.0_rp], [nsite, nsite])
      call build_charge_longitudinal_kernel(provider, coulomb, kernel)
      call check_complex_matrix('coupled local/Hartree kernel', kernel, expected_longitudinal_kernel())
      call provider%longitudinal_response_capability(supported, reason)
      call check_true('longitudinal derivative capability is explicit', supported .and. len_trim(reason) > 0)

      ! Production records keep raw radial moments until the response-site
      ! population is available; verify that normalization is not performed
      ! against an unrelated radial spin integral.
      call projected_provider%initialize(1, 'synthetic-LDA')
      projection%charge_population = 2.0_rp
      projection%spin_population = 1.0_rp
      projection%dvxc_dn_moment = 8.0_rp
      projection%dvxc_dm_moment = 6.0_rp
      projection%dbxc_dn_moment = 4.0_rp
      projection%dbxc_dm_moment = 2.0_rp
      projection%has_longitudinal_derivatives = .true.
      call projected_provider%record_radial_projection(1, projection)
      call projected_provider%set_site_spin_population(1, 1.0_rp)
      call build_charge_longitudinal_kernel(projected_provider, reshape([0.0_rp], [1, 1]), projected_kernel)
      call check_complex_matrix('radial longitudinal derivative normalization', projected_kernel, &
         cmplx(reshape([2.0_rp, 2.0_rp, 3.0_rp, 2.0_rp], [2, 2]), 0.0_rp, rp))
   end subroutine test_coupled_basis_and_kernel

   subroutine test_coupled_eigenpair_adapter()
      integer, parameter :: nsite = 1, nmat = 4, nk = 1, nw = 1
      real(rp) :: weights(nk), eigenvalues(nmat, nk), omega(nw)
      complex(rp) :: eigenvectors(nmat, nmat, nk)
      type(tddft_chi0_options) :: options
      type(tddft_chi0_result) :: result

      weights = 1.0_rp
      eigenvalues(:, 1) = [-1.0_rp, -0.5_rp, 0.5_rp, 1.0_rp]
      eigenvectors = cmplx(0.0_rp, 0.0_rp, rp)
      eigenvectors(:, :, 1) = cmplx(0.0_rp, 0.0_rp, rp)
      eigenvectors(1, 1, 1) = 1.0_rp
      eigenvectors(2, 2, 1) = 1.0_rp
      eigenvectors(3, 3, 1) = 1.0_rp
      eigenvectors(4, 4, 1) = 1.0_rp
      omega = [0.2_rp]
      options%eta = 0.01_rp
      options%fermi_level = 0.0_rp
      options%electronic_temperature = 0.0_rp
      call build_charge_longitudinal_chi_ks(weights, eigenvalues, eigenvectors, eigenvalues, eigenvectors, [2], omega, &
         options, result)
      call check_true('coupled eigenpair adapter returns (0,z) matrix', allocated(result%chi) .and. &
         all(shape(result%chi) == [2, 2, 1]))
   end subroutine test_coupled_eigenpair_adapter

   function expected_longitudinal_kernel() result(expected)
      complex(rp) :: expected(4, 4)

      expected = cmplx(0.0_rp, 0.0_rp, rp)
      expected(1, 1) = 3.0_rp
      expected(1, 2) = 2.0_rp
      expected(2, 1) = 3.0_rp
      expected(2, 2) = 4.0_rp
      expected(1, 3) = 0.2_rp
      expected(3, 1) = 0.2_rp
      expected(3, 3) = 8.0_rp
      expected(3, 4) = 6.0_rp
      expected(4, 3) = 7.0_rp
      expected(4, 4) = 8.0_rp
   end function expected_longitudinal_kernel

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
