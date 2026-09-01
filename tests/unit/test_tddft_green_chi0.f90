!------------------------------------------------------------------------------
! TDDFT-06 -- K-space Green-function chi_KS validation
!------------------------------------------------------------------------------
program test_tddft_green_chi0
   use precision_mod, only: rp
   use math_mod, only: pi
   use response_components_mod, only: RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, build_chi_ks_from_eigenpairs, &
      build_static_chi_ks_from_eigenpairs
   use tddft_chi0_green_mod, only: green_chi0_options, eigenpair_green_function_provider, &
      build_chi_ks_from_green_functions, build_static_chi_ks_from_green_functions
   use tddft_ward_mod, only: tddft_ward_diagnostics, evaluate_static_ward_identity
   implicit none

   logical :: failed

   failed = .false.
   call test_analytic_two_level_fixture()
   call test_periodic_green_bubble_matches_lehmann()
   call test_open_two_site_local_response()
   call test_static_limit_and_ward_residual()
   call test_circular_frequency_symmetry()
   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   ! A single spin-split level gives an analytic circular response.  With the
   ! frozen measurement convention sigma_+/- = sigma_x +/- i sigma_y, the
   ! matrix-element product is 2*2=4 and
   ! chi^{+-}(w) = 4/(w + e_up - e_down + i*eta).
   subroutine test_analytic_two_level_fixture()
      real(rp) :: weights(1) = 1.0_rp, omega(3) = [0.11_rp, 0.31_rp, 0.37_rp]
      real(rp) :: eval(2, 1), evalq(2, 1)
      complex(rp) :: evec(2, 2, 1), evecq(2, 2, 1), expected
      type(green_chi0_options) :: green_options
      type(tddft_chi0_result) :: result
      type(eigenpair_green_function_provider) :: source
      type(response_channel) :: left(1), right(1)
      integer :: iw

      eval(:, 1) = [-0.10_rp, 0.10_rp]
      evalq = eval
      evec = cmplx(0.0_rp, 0.0_rp, rp)
      evec(1, 1, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      evec(2, 2, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      evecq = evec
      left(1) = response_channel(1, RESPONSE_PLUS)
      right(1) = response_channel(1, RESPONSE_MINUS)
      green_options%eta = 0.003_rp
      green_options%fermi_level = 0.0_rp
      green_options%energy_min = -1.0_rp
      green_options%energy_max = 1.0_rp
      green_options%energy_points = 16001
      call source%initialize(eval, evec, evalq, evecq)
      call build_chi_ks_from_green_functions(source, weights, [1], left, right, omega, green_options, result)
      do iw = 1, size(omega)
         expected = 4.0_rp/cmplx(omega(iw)-0.20_rp, green_options%eta, rp)
         call check_close_scalar('analytic two-level GF bubble', result%chi(1, 1, iw), expected, 2.0e-2_rp)
         write (*, '(a,1x,i0,6(1x,es18.10))') 'POINTWISE chi0 analytic/GF', iw, omega(iw), &
            real(result%chi(1, 1, iw), rp), aimag(result%chi(1, 1, iw)), real(expected, rp), aimag(expected), &
            abs(result%chi(1, 1, iw)-expected)
      end do
   end subroutine test_analytic_two_level_fixture

   ! Several commensurate q points and frequencies are compared on a finite
   ! periodic spin-split chain.  The GF source only supplies resolvents, while
   ! the response path performs the real-axis convolution independently.
   subroutine test_periodic_green_bubble_matches_lehmann()
      integer, parameter :: nk = 4
      integer :: iq, iw
      integer :: q_shifts(3) = [0, 1, 2]
      real(rp) :: weights(nk), omega(4) = [0.025_rp, 0.075_rp, 0.110_rp, 0.170_rp]
      real(rp) :: eval(2, nk), evalq(2, nk)
      complex(rp) :: evec(2, 2, nk), evecq(2, 2, nk)
      type(tddft_chi0_options) :: eigen_options
      type(green_chi0_options) :: green_options
      type(tddft_chi0_result) :: eigen_result, green_result
      type(eigenpair_green_function_provider), target :: source
      type(response_channel) :: left(1), right(1)

      left(1) = response_channel(1, RESPONSE_PLUS)
      right(1) = response_channel(1, RESPONSE_MINUS)
      weights = [1.0_rp, 2.0_rp, 1.0_rp, 3.0_rp]
      eigen_options%eta = 0.001_rp
      eigen_options%fermi_level = 0.0_rp
      eigen_options%electronic_temperature = 2000.0_rp
      eigen_options%k_mesh_shape = [nk, 1, 1]
      green_options%eta = eigen_options%eta
      green_options%green_eta = 0.5_rp*eigen_options%eta
      green_options%fermi_level = eigen_options%fermi_level
      green_options%electronic_temperature = eigen_options%electronic_temperature
      green_options%energy_min = -0.35_rp
      green_options%energy_max = 0.35_rp
      green_options%energy_points = 16001
      green_options%k_mesh_shape = [nk, 1, 1]

      do iq = 1, size(q_shifts)
         call build_spin_split_fixture(q_shifts(iq), 0.030_rp, eval, evalq, evec, evecq)
         green_options%q_direct = [real(q_shifts(iq), rp)/real(nk, rp), 0.0_rp, 0.0_rp]
         call build_chi_ks_from_eigenpairs(weights, eval, evec, evalq, evecq, [1], left, right, omega, &
            eigen_options, eigen_result)
         call source%initialize(eval, evec, evalq, evecq)
         call build_chi_ks_from_green_functions(source, weights, [1], left, right, omega, green_options, green_result)
         call check_close_vector('periodic GF bubble vs eigenpair Lehmann', green_result%chi(1, 1, :), &
            eigen_result%chi(1, 1, :), 2.0e-2_rp)
         do iw = 1, size(omega)
            write (*, '(a,2(1x,i0),6(1x,es18.10))') 'POINTWISE chi0 GF/eigenpair', q_shifts(iq), iw, omega(iw), &
               real(green_result%chi(1, 1, iw), rp), aimag(green_result%chi(1, 1, iw)), &
               real(eigen_result%chi(1, 1, iw), rp), aimag(eigen_result%chi(1, 1, iw)), &
               abs(green_result%chi(1, 1, iw)-eigen_result%chi(1, 1, iw))
         end do
         call check_true('GF metadata identifies real-axis route', trim(green_result%metadata%backend) == 'green' .and. &
            green_result%metadata%integration_energy_points == green_options%energy_points)
         call check_true('GF metadata carries exact endpoint and q provenance', &
            index(trim(green_result%metadata%endpoint_provenance), 'K/K+q') > 0 .and. &
            maxval(abs(green_result%metadata%q_direct-green_options%q_direct)) < 1.0e-14_rp)
         call check_true('GF metadata carries frequency-grid provenance', green_result%metadata%omega_points == size(omega) .and. &
            abs(green_result%metadata%omega_min-minval(omega)) < 1.0e-14_rp .and. &
            abs(green_result%metadata%omega_max-maxval(omega)) < 1.0e-14_rp)
         call check_true('GF metadata carries integration controls', trim(green_result%metadata%energy_integration) == 'real-axis trapezoid' .and. &
            abs(green_result%metadata%green_eta-green_options%green_eta) < 1.0e-14_rp .and. &
            abs(green_result%metadata%integration_energy_min-green_options%energy_min) < 1.0e-14_rp .and. &
            abs(green_result%metadata%integration_energy_max-green_options%energy_max) < 1.0e-14_rp)
         call check_real('GF effective response eta', green_result%metadata%eta, eigen_options%eta, 64.0_rp*epsilon(1.0_rp))
      end do
   end subroutine test_periodic_green_bubble_matches_lehmann

   ! This is an open two-site (surface/impurity-like) local-response geometry:
   ! site 1 and site 2 have inequivalent on-site levels.  The same source API
   !> can be implemented by a real-space recursion/embedded solver, with no
   !> change to the bubble or TDDFT Dyson layers.
   subroutine test_open_two_site_local_response()
      real(rp) :: weights(1) = 1.0_rp, omega(2) = [0.06_rp, 0.12_rp]
      real(rp) :: eval(4, 1), evalq(4, 1)
      complex(rp) :: evec(4, 4, 1), evecq(4, 4, 1)
      type(tddft_chi0_options) :: eigen_options
      type(green_chi0_options) :: green_options
      type(tddft_chi0_result) :: eigen_result, green_result
      type(eigenpair_green_function_provider), target :: source
      type(response_channel) :: left(2), right(2)

      eval(:, 1) = [-0.09_rp, 0.07_rp, -0.03_rp, 0.13_rp]
      evalq = eval
      evec = cmplx(0.0_rp, 0.0_rp, rp)
      evecq = evec
      evec(1, 1, 1) = cmplx(1.0_rp, 0.0_rp, rp); evec(2, 2, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      evec(3, 3, 1) = cmplx(1.0_rp, 0.0_rp, rp); evec(4, 4, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      evecq = evec
      left(1) = response_channel(1, RESPONSE_PLUS); left(2) = response_channel(2, RESPONSE_PLUS)
      right(1) = response_channel(1, RESPONSE_MINUS); right(2) = response_channel(2, RESPONSE_MINUS)
      eigen_options%eta = 0.001_rp
      eigen_options%electronic_temperature = 2000.0_rp
      green_options%eta = eigen_options%eta
      green_options%green_eta = 0.0005_rp
      green_options%electronic_temperature = eigen_options%electronic_temperature
      green_options%energy_min = -0.35_rp; green_options%energy_max = 0.35_rp
      green_options%energy_points = 16001
      call build_chi_ks_from_eigenpairs(weights, eval, evec, evalq, evecq, [1, 1], left, right, omega, eigen_options, eigen_result)
      call source%initialize(eval, evec, evalq, evecq)
      call build_chi_ks_from_green_functions(source, weights, [1, 1], left, right, omega, green_options, green_result)
      call check_close_3d('open two-site local response', green_result%chi(:, :, :), eigen_result%chi(:, :, :), 2.0e-2_rp)
      call check_true('inequivalent local sites retain different response', &
         abs(green_result%chi(1, 1, 1)-green_result%chi(2, 2, 1)) > 1.0e-4_rp)
   end subroutine test_open_two_site_local_response

   subroutine test_static_limit_and_ward_residual()
      real(rp) :: weights(1) = 1.0_rp
      real(rp) :: eval(2, 1)
      complex(rp) :: evec(2, 2, 1), bxc(1)
      type(tddft_chi0_options) :: eigen_options
      type(green_chi0_options) :: green_options
      type(tddft_chi0_result) :: eigen_result, green_result
      type(tddft_ward_diagnostics) :: eigen_ward, green_ward
      type(eigenpair_green_function_provider) :: source
      type(response_channel) :: left(1), right(1)

      eval(:, 1) = [-0.10_rp, 0.10_rp]
      evec = cmplx(0.0_rp, 0.0_rp, rp)
      evec(1, 1, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      evec(2, 2, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      left(1) = response_channel(1, RESPONSE_PLUS)
      right(1) = response_channel(1, RESPONSE_MINUS)
      eigen_options%fermi_level = 0.0_rp
      eigen_options%electronic_temperature = 0.0_rp
      green_options%fermi_level = eigen_options%fermi_level
      green_options%electronic_temperature = eigen_options%electronic_temperature
      green_options%q_direct = 0.0_rp
      call source%initialize(eval, evec, eval, evec)
      call build_static_chi_ks_from_eigenpairs(weights, eval, evec, [1], left, right, eigen_options, eigen_result)
      call build_static_chi_ks_from_green_functions(source, weights, [1], left, right, green_options, green_result)
      call check_close_3d('static GF chi0 vs eigenpair chi0', green_result%chi, eigen_result%chi, 256.0_rp*epsilon(1.0_rp))
      call check_true('static GF metadata identifies an exact static limit', green_result%metadata%static_limit .and. &
         index(trim(green_result%metadata%energy_integration), 'exact provider') > 0 .and. &
         green_result%metadata%integration_energy_points == 0)

      ! The analytic static response is -20 for this fixture.  B_xc=-0.05
      ! consequently gives chi_KS B_xc = m=1 and a zero Ward residual.
      bxc(1) = cmplx(-0.05_rp, 0.0_rp, rp)
      call evaluate_static_ward_identity(eigen_result%chi(:, :, 1), bxc, [cmplx(1.0_rp, 0.0_rp, rp)], eigen_ward, &
         response_basis='site circular', bxc_provenance='analytic fixture', kernel_provenance='analytic fixture kernel')
      call evaluate_static_ward_identity(green_result%chi(:, :, 1), bxc, [cmplx(1.0_rp, 0.0_rp, rp)], green_ward, &
         response_basis='site circular', bxc_provenance='analytic fixture', kernel_provenance='analytic fixture kernel')
      call check_close_real('K-GF Ward residual vs eigenpair Ward residual', green_ward%ward_residual, &
         eigen_ward%ward_residual, 256.0_rp*epsilon(1.0_rp))
      call check_true('K-GF Ward residual is converged on analytic fixture', green_ward%ward_residual < 256.0_rp*epsilon(1.0_rp))
      write (*, '(a,5(1x,es18.10))') 'STATIC/WARD chi0_GF Re Im residual_GF residual_eigen abs_delta', &
         real(green_result%chi(1, 1, 1), rp), aimag(green_result%chi(1, 1, 1)), green_ward%ward_residual, &
         eigen_ward%ward_residual, abs(green_result%chi(1, 1, 1)-eigen_result%chi(1, 1, 1))
   end subroutine test_static_limit_and_ward_residual

   subroutine test_circular_frequency_symmetry()
      real(rp) :: weights(1) = 1.0_rp, omega_positive(2) = [0.13_rp, 0.31_rp], omega_negative(2)
      real(rp) :: eval(2, 1)
      complex(rp) :: evec(2, 2, 1)
      type(green_chi0_options) :: options
      type(tddft_chi0_result) :: positive, negative
      type(eigenpair_green_function_provider) :: source
      type(response_channel) :: left(2), right(2)
      integer :: iw

      eval(:, 1) = [-0.10_rp, 0.10_rp]
      evec = cmplx(0.0_rp, 0.0_rp, rp)
      evec(1, 1, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      evec(2, 2, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      left = [response_channel(1, RESPONSE_PLUS), response_channel(1, RESPONSE_MINUS)]
      right = [response_channel(1, RESPONSE_MINUS), response_channel(1, RESPONSE_PLUS)]
      options%eta = 0.003_rp
      options%energy_min = -1.0_rp; options%energy_max = 1.0_rp; options%energy_points = 12001
      omega_negative = -omega_positive
      call source%initialize(eval, evec, eval, evec)
      call build_chi_ks_from_green_functions(source, weights, [1], left, right, omega_positive, options, positive)
      call build_chi_ks_from_green_functions(source, weights, [1], left, right, omega_negative, options, negative)
      do iw = 1, size(omega_positive)
         call check_close_scalar('circular retarded frequency/channel symmetry', positive%chi(1, 2, iw), &
            conjg(negative%chi(2, 1, iw)), 3.0e-3_rp)
      end do
   end subroutine test_circular_frequency_symmetry

   subroutine build_spin_split_fixture(q_shift, hopping, eval, evalq, evec, evecq)
      integer, intent(in) :: q_shift
      real(rp), intent(in) :: hopping
      real(rp), intent(out) :: eval(:, :), evalq(:, :)
      complex(rp), intent(out) :: evec(:, :, :), evecq(:, :, :)
      integer :: ik, ikq
      real(rp) :: k, kq

      do ik = 1, size(eval, 2)
         ikq = 1 + modulo(ik - 1 + q_shift, size(eval, 2))
         k = real(ik-1, rp)/real(size(eval, 2), rp)
         kq = real(ikq-1, rp)/real(size(eval, 2), rp)
         eval(1, ik) = -0.050_rp-2.0_rp*hopping*cos(2.0_rp*pi*k)
         eval(2, ik) =  0.050_rp-2.0_rp*hopping*cos(2.0_rp*pi*k)
         evalq(1, ik) = -0.050_rp-2.0_rp*hopping*cos(2.0_rp*pi*kq)
         evalq(2, ik) =  0.050_rp-2.0_rp*hopping*cos(2.0_rp*pi*kq)
         evec(:, :, ik) = cmplx(0.0_rp, 0.0_rp, rp); evecq(:, :, ik) = evec(:, :, ik)
         evec(1, 1, ik) = cmplx(1.0_rp, 0.0_rp, rp); evec(2, 2, ik) = cmplx(1.0_rp, 0.0_rp, rp)
         evecq(1, 1, ik) = cmplx(1.0_rp, 0.0_rp, rp); evecq(2, 2, ik) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
   end subroutine build_spin_split_fixture

   subroutine check_close_vector(label, actual, expected, tolerance)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual(:), expected(:)
      real(rp), intent(in) :: tolerance
      real(rp) :: relative_error, scale

      scale = max(1.0_rp, maxval(abs(expected)))
      relative_error = maxval(abs(actual-expected))/scale
      if (relative_error > tolerance) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, relative_error
         failed = .true.
      end if
   end subroutine check_close_vector

   subroutine check_close_scalar(label, actual, expected, tolerance)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual, expected
      real(rp), intent(in) :: tolerance
      real(rp) :: relative_error, scale

      scale = max(1.0_rp, abs(expected))
      relative_error = abs(actual-expected)/scale
      if (relative_error > tolerance) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, relative_error
         failed = .true.
      end if
   end subroutine check_close_scalar

   subroutine check_close_real(label, actual, expected, tolerance)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: actual, expected, tolerance
      real(rp) :: relative_error, scale

      scale = max(1.0_rp, abs(expected))
      relative_error = abs(actual-expected)/scale
      if (relative_error > tolerance) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, relative_error
         failed = .true.
      end if
   end subroutine check_close_real

   subroutine check_close_3d(label, actual, expected, tolerance)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual(:, :, :), expected(:, :, :)
      real(rp), intent(in) :: tolerance
      real(rp) :: relative_error, scale

      scale = max(1.0_rp, maxval(abs(expected)))
      relative_error = maxval(abs(actual-expected))/scale
      if (relative_error > tolerance) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, relative_error
         failed = .true.
      end if
   end subroutine check_close_3d

   subroutine check_true(label, value)
      character(len=*), intent(in) :: label
      logical, intent(in) :: value
      if (.not. value) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine check_true

   subroutine check_real(label, actual, expected, tolerance)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: actual, expected, tolerance
      if (abs(actual-expected) > tolerance*max(1.0_rp, abs(expected))) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine check_real

end program test_tddft_green_chi0
