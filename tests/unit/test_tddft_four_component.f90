!------------------------------------------------------------------------------
! TDDFT-09 -- full charge-spin response, local ALSDA frames, and SOC fixtures
!------------------------------------------------------------------------------
program test_tddft_four_component
   use precision_mod, only: rp
   use response_components_mod, only: RESPONSE_PLUS, RESPONSE_MINUS, RESPONSE_MZ
   use response_vertices_mod, only: response_channel
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, build_chi_ks_from_eigenpairs
   use tddft_four_component_mod, only: build_four_component_chi_ks, build_four_component_kernel, &
      evaluate_four_component_zero_modes, tddft_four_component_zero_mode_diagnostics, four_component_index
   use xc_response_kernel_mod, only: xc_response_kernel_provider
   implicit none

   real(rp), parameter :: tol = 2.0e-12_rp
   logical :: failed

   failed = .false.
   call test_collinear_reduction_and_soc_mixing_fixture()
   call test_full_response_capability_gate()
   call test_local_alsda_hartree_rotation_covariance()
   call test_rigid_rotation_diagnostics()
   if (failed) error stop 1
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_collinear_reduction_and_soc_mixing_fixture()
      real(rp) :: weights(1), eval(2, 1), omega(1), angle
      complex(rp) :: evec(2, 2, 1), evecq(2, 2, 1), transverse(1, 1, 1)
      type(tddft_chi0_options) :: options
      type(tddft_chi0_result) :: full, legacy, legacy_longitudinal
      type(response_channel) :: left(1), right(1)
      integer :: ix, iy, iz, ic

      weights = 1.0_rp; eval(:, 1) = [-0.10_rp, 0.10_rp]; omega = 0.07_rp
      options%eta = 0.003_rp
      evec = cmplx(0.0_rp, 0.0_rp, rp); evec(1, 1, 1) = 1.0_rp; evec(2, 2, 1) = 1.0_rp
      evecq = evec
      call build_four_component_chi_ks(weights, eval, evec, eval, evecq, [1], omega, options, full)
      left(1) = response_channel(1, RESPONSE_PLUS)
      right(1) = response_channel(1, RESPONSE_MINUS)
      call build_chi_ks_from_eigenpairs(weights, eval, evec, eval, evecq, [1], left, right, omega, options, legacy)
      ic = four_component_index(1, 0); ix = four_component_index(1, 1)
      iy = four_component_index(1, 2); iz = four_component_index(1, RESPONSE_MZ)
      transverse = full%chi(ix:ix, ix:ix, :) + full%chi(iy:iy, iy:iy, :) + &
         cmplx(0.0_rp, 1.0_rp, rp)*full%chi(iy:iy, ix:ix, :) - cmplx(0.0_rp, 1.0_rp, rp)*full%chi(ix:ix, iy:iy, :)
      call check_complex('collinear transverse is recovered from Cartesian full tensor', transverse(1, 1, 1), legacy%chi(1, 1, 1))
      left(1) = response_channel(1, RESPONSE_MZ)
      right(1) = response_channel(1, RESPONSE_MZ)
      call build_chi_ks_from_eigenpairs(weights, eval, evec, eval, evecq, [1], left, right, omega, options, legacy_longitudinal)
      call check_complex('collinear longitudinal is recovered from full tensor', full%chi(iz, iz, 1), &
         legacy_longitudinal%chi(1, 1, 1))
      call check_real('collinear charge-transverse block decouples', abs(full%chi(ic, ix, 1)) + abs(full%chi(ic, iy, 1)), 0.0_rp)

      ! Minimal SOC-like fixture: k and k+q spinors have different spin axes.
      ! The generalized vertices must retain charge-spin mixing rather than
      ! dropping it using a collinear selection rule.
      angle = 0.37_rp
      evecq = cmplx(0.0_rp, 0.0_rp, rp)
      evecq(1, 1, 1) = cos(angle); evecq(2, 1, 1) = sin(angle)
      evecq(1, 2, 1) = -sin(angle); evecq(2, 2, 1) = cos(angle)
      call build_four_component_chi_ks(weights, eval, evec, eval, evecq, [1], omega, options, full)
      call check_true('spin-mixed fixture retains charge-spin coupling', abs(full%chi(ic, ix, 1)) > tol)
   end subroutine test_collinear_reduction_and_soc_mixing_fixture

   subroutine test_full_response_capability_gate()
      type(xc_response_kernel_provider) :: provider
      logical :: supported
      character(len=256) :: reason

      call provider%initialize(1, 'unit-unvalidated-XC')
      call provider%set_site_magnetization_direction(1, [0.0_rp, 0.0_rp, 1.0_rp])
      call provider%set_site_derivatives(1, dvxc_dn=1.0_rp, dvxc_dm=0.2_rp, dbxc_dn=0.2_rp, dbxc_dm=1.4_rp, &
         k_perp_circular=0.7_rp)
      call provider%full_response_capability(supported, reason)
      call check_true('unvalidated full ALSDA derivatives are rejected', .not. supported .and. &
         index(reason, 'not been validated') > 0)
      call provider%set_site_derivatives(1, derivatives_validated=.true.)
      call provider%full_response_capability(supported, reason)
      call check_true('validated complete full ALSDA derivatives are accepted', supported)
   end subroutine test_full_response_capability_gate

   subroutine test_local_alsda_hartree_rotation_covariance()
      type(xc_response_kernel_provider) :: local_z, rotated
      complex(rp) :: kernel_z(4, 4), kernel_rotated(4, 4), expected(4, 4)
      real(rp) :: coulomb(1, 1), q(3, 3), transform(4, 4), theta

      coulomb(1, 1) = 0.6_rp
      call make_provider(local_z, [0.0_rp, 0.0_rp, 1.0_rp])
      call build_four_component_kernel(local_z, coulomb, kernel_z)
      call check_real('Hartree enters charge-charge only', real(kernel_z(1, 1), rp), 1.6_rp)
      call check_real('Cartesian transverse ALSDA x is twice circular scalar', real(kernel_z(2, 2), rp), 1.4_rp)
      call check_real('Cartesian transverse ALSDA y is twice circular scalar', real(kernel_z(3, 3), rp), 1.4_rp)
      call check_real('longitudinal ALSDA', real(kernel_z(4, 4), rp), 1.4_rp)

      theta = 0.31_rp
      q = reshape([cos(theta), 0.0_rp, -sin(theta), 0.0_rp, 1.0_rp, 0.0_rp, sin(theta), 0.0_rp, cos(theta)], [3, 3])
      call make_provider(rotated, matmul(q, [0.0_rp, 0.0_rp, 1.0_rp]))
      call build_four_component_kernel(rotated, coulomb, kernel_rotated)
      transform = 0.0_rp; transform(1, 1) = 1.0_rp; transform(2:4, 2:4) = q
      expected = matmul(cmplx(transform, 0.0_rp, rp), matmul(kernel_z, transpose(cmplx(transform, 0.0_rp, rp))))
      call check_complex_matrix('global spin rotation covariance of ALSDA kernel', kernel_rotated, expected)
   end subroutine test_local_alsda_hartree_rotation_covariance

   subroutine test_rigid_rotation_diagnostics()
      real(rp) :: magnetization(3, 1)
      complex(rp) :: chi(4, 4), kernel(4, 4)
      type(tddft_four_component_zero_mode_diagnostics) :: diagnostics

      magnetization(:, 1) = [0.0_rp, 0.0_rp, 2.0_rp]
      chi = cmplx(0.0_rp, 0.0_rp, rp); kernel = cmplx(0.0_rp, 0.0_rp, rp)
      chi = identity4(); kernel = identity4()
      call evaluate_four_component_zero_modes(chi, kernel, magnetization, .false., .false., diagnostics)
      call check_true('no-SOC rigid rotations are diagnosed', diagnostics%applicable .and. diagnostics%number_of_modes == 2)
      call check_real('no-SOC rigid-rotation residual', maxval(diagnostics%residual), 0.0_rp)
      call evaluate_four_component_zero_modes(chi, kernel, magnetization, .true., .false., diagnostics)
      call check_true('SOC preserves physical gap by disabling zero-mode enforcement', diagnostics%symmetry_broken .and. .not. diagnostics%applicable)
   end subroutine test_rigid_rotation_diagnostics

   subroutine make_provider(provider, direction)
      type(xc_response_kernel_provider), intent(out) :: provider
      real(rp), intent(in) :: direction(3)

      call provider%initialize(1, 'unit-full-ALSDA')
      call provider%set_site_magnetization_direction(1, direction)
      call provider%set_site_derivatives(1, dvxc_dn=1.0_rp, dvxc_dm=0.2_rp, dbxc_dn=0.2_rp, dbxc_dm=1.4_rp, &
         k_perp_circular=0.7_rp, derivatives_validated=.true.)
   end subroutine make_provider

   pure function identity4() result(matrix)
      complex(rp) :: matrix(4, 4)
      integer :: i
      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, 4
         matrix(i, i) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
   end function identity4

   subroutine check_complex(label, actual, expected)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual, expected
      call check_real(label, abs(actual-expected), 0.0_rp)
   end subroutine check_complex

   subroutine check_complex_matrix(label, actual, expected)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual(:, :), expected(:, :)
      call check_real(label, maxval(abs(actual-expected)), 0.0_rp)
   end subroutine check_complex_matrix

   subroutine check_real(label, actual, expected)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: actual, expected
      if (abs(actual-expected) > tol*max(1.0_rp, abs(expected))) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, abs(actual-expected)
         failed = .true.
      end if
   end subroutine check_real

   subroutine check_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine check_true

end program test_tddft_four_component
