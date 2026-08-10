!------------------------------------------------------------------------------
! TDDFT-05 -- Dyson enhancement, Xi diagnostics, branches, and damping policy
!------------------------------------------------------------------------------
program test_tddft_dyson_modes
   use precision_mod, only: rp
   use tddft_dyson_mod, only: tddft_dyson_options, tddft_dyson_result, enhance_tddft_susceptibility, &
      enhance_tddft_susceptibility_from_xi, solve_tddft_dyson_frequency, write_tddft_dyson_text
   use tddft_modes_mod, only: tddft_mode_options, tddft_mode_result, analyze_tddft_modes, &
      extrapolate_linewidth_zero_eta, TDDFT_MODE_COHERENT
   implicit none

   real(rp), parameter :: tol = 2.0e-10_rp
   logical :: failed

   failed = .false.
   call test_dyson_collective_pole_and_loss()
   call test_direct_xi_dyson_equivalence()
   call test_overlap_branch_tracking_and_fit_policy()
   call test_multi_eta_extrapolation()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_dyson_collective_pole_and_loss()
      integer, parameter :: nw = 81
      real(rp) :: omega(nw)
      complex(rp) :: chi_ks(2, 2, nw), kernel(2, 2), u(2), expected_chi(2, 2), xi(2, 2), chi(2, 2)
      type(tddft_dyson_options) :: options
      type(tddft_dyson_result) :: result
      integer :: ip, imode

      call make_collective_fixture(omega, chi_ks, kernel, u)
      options%diagonalize_loss = .true.
      options%diagonalize_xi = .true.
      call enhance_tddft_susceptibility(chi_ks, kernel, 0.001_rp, options, result)
      ip = 21
      call solve_tddft_dyson_frequency(chi_ks(:, :, ip), kernel, chi, xi, ip)
      call check_true('streaming Dyson solve succeeds', ip == 0)
      call check_complex_matrix('streaming and batch enhanced chi agree', chi, result%chi(:, :, 21))
      call check_complex_matrix('Xi retained by streaming solve', xi, result%xi(:, :, 21))
      call check_true('bare KS loss has no collective peak', maxval(result%chi_ks_site_spectral_weight(1, :)) / &
         minval(result%chi_ks_site_spectral_weight(1, :)) < 1.01_rp)
      call check_true('enhanced chi contains collective pole', result%trace_spectral_weight(21) > 10.0_rp*result%trace_spectral_weight(1))
      call check_true('Xi has unity eigenvalue at collective pole', minval(abs(result%xi_eigenvalues(:, 21)-cmplx(1.0_rp, 0.0_rp, rp))) < 0.006_rp)
      call check_true('loss matrix was diagonalized when requested', allocated(result%loss_eigenvalues))
      call check_true('loss eigenvalues are nonnegative in deterministic fixture', minval(result%loss_eigenvalues(:, 21)) > -tol)
      imode = minloc(abs(result%xi_eigenvalues(:, 21)-cmplx(1.0_rp, 0.0_rp, rp)), dim=1)
      call check_complex_matrix('known response eigenvector projector', projector_from_mode(result%xi_eigenvectors(:, imode, 21)), &
         projector_from_mode(u))
      call check_true('numerical eta is retained as metadata', abs(result%metadata%eta-0.001_rp) < tol)
      call write_tddft_dyson_text('unit_tddft_dyson.out', omega, result)
      call check_output_contains('unit_tddft_dyson.out', 'chi_KS')
      call check_output_contains('unit_tddft_dyson.out', 'Xi')
      call check_output_contains('unit_tddft_dyson.out', 'loss')
      call delete_output('unit_tddft_dyson.out')
   end subroutine test_dyson_collective_pole_and_loss

   subroutine test_direct_xi_dyson_equivalence()
      integer, parameter :: nw = 5
      complex(rp) :: chi_ks(2, 2, nw), kernel(2, 2), xi(2, 2, nw), u(2)
      real(rp) :: omega(nw)
      type(tddft_dyson_options) :: options
      type(tddft_dyson_result) :: from_kernel, from_xi
      integer :: iw

      call make_collective_fixture(omega, chi_ks, kernel, u)
      do iw = 1, nw
         xi(:, :, iw) = matmul(chi_ks(:, :, iw), kernel)
      end do
      call enhance_tddft_susceptibility(chi_ks, kernel, 0.001_rp, options, from_kernel)
      call enhance_tddft_susceptibility_from_xi(chi_ks, xi, 0.001_rp, options, from_xi)
      call check_complex_matrix('direct-Xi Dyson solve equals kernel route for uniform oracle', &
         reshape(from_xi%chi, [2, 2*nw]), reshape(from_kernel%chi, [2, 2*nw]))
      call check_complex_matrix('direct-Xi Dyson retains supplied Xi without reconstructing K', &
         reshape(from_xi%xi, [2, 2*nw]), reshape(xi, [2, 2*nw]))
   end subroutine test_direct_xi_dyson_equivalence

   subroutine test_overlap_branch_tracking_and_fit_policy()
      integer, parameter :: nw = 81, nq = 2
      real(rp) :: omega(nw), trace_loss(nw, nq)
      complex(rp) :: chi_ks(2, 2, nw), kernel(2, 2), u(2), xi(2, 2, nw, nq)
      type(tddft_dyson_options) :: dyson_options
      type(tddft_dyson_result) :: dyson
      type(tddft_mode_options) :: options
      type(tddft_mode_result) :: modes
      integer :: iq, iw

      call make_collective_fixture(omega, chi_ks, kernel, u)
      call enhance_tddft_susceptibility(chi_ks, kernel, 0.001_rp, dyson_options, dyson)
      do iq = 1, nq
         xi(:, :, :, iq) = dyson%xi
         trace_loss(:, iq) = dyson%trace_spectral_weight
      end do
      call analyze_tddft_modes(omega, xi, trace_loss, 0.001_rp, options, modes)
      call check_true('mode candidate follows Xi unity frequency', all(modes%candidate_frequency_index == 21))
      call check_true('q branch uses eigenvector overlap', modes%branch_overlap(2) > 0.999999_rp)
      call check_true('isolated enhanced mode has an accepted controlled fit', all(modes%fit%accepted))
      call check_true('fit reports explicit FWHM and HWHM', all(modes%fit%fwhm > 0.0_rp) .and. &
         maxval(abs([(modes%fit(iq)%fwhm-2.0_rp*modes%fit(iq)%hwhm, iq=1,nq)])) < tol)
      call check_true('accepted unity mode is classified collective', all(modes%classification == TDDFT_MODE_COHERENT))
      ! A deliberately distorted continuum-shaped peak must remain rejected.
      do iw = 1, nw
         trace_loss(iw, 1) = 1.0_rp + exp(-((omega(iw)-0.10_rp)/0.010_rp)**4)
      end do
      call analyze_tddft_modes(omega, xi(:, :, :, 1:1), trace_loss(:, 1:1), 0.001_rp, options, modes)
      call check_true('non-Lorentzian enhanced structure is not forced into a fit', .not. modes%fit(1)%accepted)
   end subroutine test_overlap_branch_tracking_and_fit_policy

   subroutine test_multi_eta_extrapolation()
      real(rp) :: eta(3), fwhm(3), intrinsic_fwhm, intrinsic_hwhm, slope, residual
      logical :: accepted

      eta = [0.001_rp, 0.002_rp, 0.003_rp]
      fwhm = 0.020_rp + 3.0_rp*eta
      call extrapolate_linewidth_zero_eta(eta, fwhm, intrinsic_fwhm, intrinsic_hwhm, slope, residual, accepted)
      call check_true('multi-eta extrapolation accepted', accepted)
      call check_real('zero-eta FWHM remains distinct from eta', intrinsic_fwhm, 0.020_rp)
      call check_real('zero-eta HWHM explicit', intrinsic_hwhm, 0.010_rp)
      call check_real('multi-eta slope', slope, 3.0_rp)
      call check_real('multi-eta residual', residual, 0.0_rp)
   end subroutine test_multi_eta_extrapolation

   subroutine make_collective_fixture(omega, chi_ks, kernel, u)
      real(rp), intent(out) :: omega(:)
      complex(rp), intent(out) :: chi_ks(:, :, :), kernel(:, :), u(:)
      complex(rp) :: scalar
      integer :: iw

      omega = [(0.050_rp + 0.0025_rp*real(iw-1, rp), iw=1,size(omega))]
      u = [cmplx(1.0_rp/sqrt(2.0_rp), 0.0_rp, rp), cmplx(0.0_rp, 1.0_rp/sqrt(2.0_rp), rp)]
      kernel = cmplx(0.0_rp, 0.0_rp, rp)
      kernel(1, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      kernel(2, 2) = cmplx(1.0_rp, 0.0_rp, rp)
      do iw = 1, size(omega)
         ! Bare loss is flat and small.  Xi=s*|u><u| reaches 1-i*0.005 at
         ! 0.10 Ry; Dyson therefore produces a pole absent in chi_KS.
         scalar = cmplx(1.0_rp + 0.4_rp*(omega(iw)-0.10_rp), -0.005_rp, rp)
         chi_ks(:, :, iw) = scalar*projector_from_mode(u)
      end do
   end subroutine make_collective_fixture

   function projector_from_mode(vector) result(projector)
      complex(rp), intent(in) :: vector(:)
      complex(rp) :: projector(size(vector), size(vector))
      integer :: i, j
      do j = 1, size(vector)
         do i = 1, size(vector)
            projector(i, j) = vector(i)*conjg(vector(j))
         end do
      end do
   end function projector_from_mode

   subroutine check_complex_matrix(label, actual, expected)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual(:, :), expected(:, :)
      call check_true(label, maxval(abs(actual-expected)) < tol)
   end subroutine check_complex_matrix

   subroutine check_real(label, actual, expected)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: actual, expected
      call check_true(label, abs(actual-expected) < tol*max(1.0_rp, abs(expected)))
   end subroutine check_real

   subroutine check_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine check_true

   subroutine check_output_contains(filename, pattern)
      character(len=*), intent(in) :: filename, pattern
      character(len=256) :: line
      integer :: unit, ios
      logical :: found

      found = .false.
      open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios == 0) then
         do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, pattern) > 0) found = .true.
         end do
         close(unit)
      end if
      call check_true('diagnostic text contains '//trim(pattern), found)
   end subroutine check_output_contains

   subroutine delete_output(filename)
      character(len=*), intent(in) :: filename
      integer :: unit, ios
      open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios == 0) close(unit, status='delete')
   end subroutine delete_output

end program test_tddft_dyson_modes
