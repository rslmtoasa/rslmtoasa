!------------------------------------------------------------------------------
! TDDFT-10 -- real-axis Green-function chi_KS validation
!------------------------------------------------------------------------------
program test_tddft_green_chi0
   use precision_mod, only: rp
   use math_mod, only: pi
   use response_components_mod, only: RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, build_chi_ks_from_eigenpairs
   use tddft_chi0_green_mod, only: green_chi0_options, eigenpair_green_function_provider, &
      build_chi_ks_from_green_functions
   implicit none

   logical :: failed

   failed = .false.
   call test_periodic_green_bubble_matches_lehmann()
   call test_open_two_site_local_response()
   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   ! Several commensurate q points and frequencies are compared on a finite
   ! periodic spin-split chain.  The GF source only supplies resolvents, while
   ! the response path performs the real-axis convolution independently.
   subroutine test_periodic_green_bubble_matches_lehmann()
      integer, parameter :: nk = 4
      integer :: iq
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
         call build_chi_ks_from_eigenpairs(weights, eval, evec, evalq, evecq, [1], left, right, omega, &
            eigen_options, eigen_result)
         call source%initialize(eval, evec, evalq, evecq)
         call build_chi_ks_from_green_functions(source, weights, [1], left, right, omega, green_options, green_result)
         call check_close_vector('periodic GF bubble vs eigenpair Lehmann', green_result%chi(1, 1, :), &
            eigen_result%chi(1, 1, :), 2.0e-2_rp)
         call check_true('GF metadata identifies real-axis route', trim(green_result%metadata%backend) == 'green' .and. &
            green_result%metadata%integration_energy_points == green_options%energy_points)
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
