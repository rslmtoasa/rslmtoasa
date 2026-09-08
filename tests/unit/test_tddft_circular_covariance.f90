!------------------------------------------------------------------------------
! TDCOV-03 -- circular-response covariance and ordered-channel regression.
!------------------------------------------------------------------------------
program test_tddft_circular_covariance
   use precision_mod, only: rp
   use response_components_mod, only: RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel
   use tddft_circular_mod, only: TDDFT_CIRCULAR_BOTH, TDDFT_CIRCULAR_PLUS_MINUS, &
      TDDFT_CIRCULAR_MINUS_PLUS, circular_channel_code, circular_channel_name, circular_channel_components, &
      circular_channel_file_tag, opposite_circular_channel
   use tddft_occupation_mod, only: tddft_response_occupation_state
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, build_chi_ks_from_eigenpairs, &
      build_static_chi_ks_from_eigenpairs, install_tddft_occupation
   implicit none

   logical :: failed

   failed = .false.
   call test_channel_object_contract()
   call test_gamma_covariance_and_shared_occupation()
   call test_finite_q_frequency_covariance()
   call test_global_spin_reversal()
   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_channel_object_contract()
      integer :: left_component, right_component

      call assert_true('both channel code is explicit', circular_channel_code('both') == TDDFT_CIRCULAR_BOTH)
      call assert_true('plus-minus aliases resolve', circular_channel_code('+-') == TDDFT_CIRCULAR_PLUS_MINUS)
      call assert_true('minus-plus aliases resolve', circular_channel_code('-+') == TDDFT_CIRCULAR_MINUS_PLUS)
      call assert_true('opposite plus-minus channel', &
         opposite_circular_channel(TDDFT_CIRCULAR_PLUS_MINUS) == TDDFT_CIRCULAR_MINUS_PLUS)
      call assert_true('opposite minus-plus channel', &
         opposite_circular_channel(TDDFT_CIRCULAR_MINUS_PLUS) == TDDFT_CIRCULAR_PLUS_MINUS)

      call circular_channel_components(TDDFT_CIRCULAR_PLUS_MINUS, left_component, right_component)
      call assert_true('plus-minus retains ordered measurement/source components', &
         left_component == RESPONSE_PLUS .and. right_component == RESPONSE_MINUS)
      call circular_channel_components(TDDFT_CIRCULAR_MINUS_PLUS, left_component, right_component)
      call assert_true('minus-plus retains ordered measurement/source components', &
         left_component == RESPONSE_MINUS .and. right_component == RESPONSE_PLUS)
      call assert_true('historical plus-minus file tag is preserved', &
         trim(circular_channel_file_tag(TDDFT_CIRCULAR_PLUS_MINUS, historical_primary=.true.)) == '')
      call assert_true('explicit plus-minus file tag is available', &
         trim(circular_channel_file_tag(TDDFT_CIRCULAR_PLUS_MINUS)) == '_plus_minus')
      call assert_true('minus-plus file tag is explicit', &
         trim(circular_channel_file_tag(TDDFT_CIRCULAR_MINUS_PLUS)) == '_minus_plus')
      call assert_true('channel names are canonical', trim(circular_channel_name(TDDFT_CIRCULAR_PLUS_MINUS)) == 'plus_minus' .and. &
         trim(circular_channel_name(TDDFT_CIRCULAR_MINUS_PLUS)) == 'minus_plus')
   end subroutine test_channel_object_contract

   subroutine test_gamma_covariance_and_shared_occupation()
      real(rp) :: weights(1), eigenvalues(2, 1), omega(1)
      complex(rp) :: eigenvectors(2, 2, 1)
      type(response_channel) :: plus(1), minus(1)
      type(tddft_chi0_options) :: plus_options, minus_options
      type(tddft_chi0_result) :: plus_result, minus_result, plus_static, minus_static
      real(rp) :: dynamic_residual, static_residual

      weights = 1.0_rp
      eigenvalues(:, 1) = [-0.05_rp, 0.05_rp]
      eigenvectors = cmplx(0.0_rp, 0.0_rp, rp)
      eigenvectors(1, 1, 1) = 1.0_rp
      eigenvectors(2, 2, 1) = 1.0_rp
      omega = 0.0_rp
      plus = [response_channel(1, RESPONSE_PLUS)]
      minus = [response_channel(1, RESPONSE_MINUS)]

      call setup_options(plus_options, TDDFT_CIRCULAR_PLUS_MINUS, [1, 1, 1])
      minus_options = plus_options
      minus_options%circular_channel = 'minus_plus'
      call build_chi_ks_from_eigenpairs(weights, eigenvalues, eigenvectors, eigenvalues, eigenvectors, [1], plus, minus, &
         omega, plus_options, plus_result)
      call build_chi_ks_from_eigenpairs(weights, eigenvalues, eigenvectors, eigenvalues, eigenvectors, [1], minus, plus, &
         omega, minus_options, minus_result)
      dynamic_residual = maxval(abs(plus_result%chi(:, :, 1)-transpose(conjg(minus_result%chi(:, :, 1)))))
      write (*, '(a,1x,es12.4)') 'OBSERVED gamma_dynamic_covariance_residual', dynamic_residual
      call assert_true('Gamma dynamic covariance', dynamic_residual < 1.0e-12_rp)
      call assert_true('Gamma plus-minus metadata is actual channel', &
         trim(plus_result%metadata%circular_channel) == 'plus_minus')
      call assert_true('Gamma minus-plus metadata is actual channel', &
         trim(minus_result%metadata%circular_channel) == 'minus_plus')
      call assert_same_occupation('Gamma channels share dynamic occupation metadata', plus_result, minus_result)

      call build_static_chi_ks_from_eigenpairs(weights, eigenvalues, eigenvectors, [1], plus, minus, plus_options, plus_static)
      call build_static_chi_ks_from_eigenpairs(weights, eigenvalues, eigenvectors, [1], minus, plus, minus_options, minus_static)
      static_residual = maxval(abs(plus_static%chi(:, :, 1)-transpose(conjg(minus_static%chi(:, :, 1)))))
      write (*, '(a,1x,es12.4)') 'OBSERVED gamma_static_covariance_residual', static_residual
      call assert_true('Gamma static covariance', static_residual < 1.0e-12_rp)
      call assert_same_occupation('Gamma channels share static occupation metadata', plus_static, minus_static)
   end subroutine test_gamma_covariance_and_shared_occupation

   subroutine test_finite_q_frequency_covariance()
      integer, parameter :: nk = 4, nbands = 4, nsite = 2, nw = 3
      real(rp) :: weights(nk), eigenvalues(nbands, nk), eigenvalues_plus(nbands, nk), &
         eigenvalues_minus(nbands, nk), omega(nw), angle, cosine, sine
      complex(rp) :: eigenvectors(nbands, nbands, nk), eigenvectors_plus(nbands, nbands, nk), &
         eigenvectors_minus(nbands, nbands, nk)
      type(response_channel) :: plus(nsite), minus(nsite)
      type(tddft_chi0_options) :: plus_options, minus_options
      type(tddft_chi0_result) :: plus_result, minus_result
      real(rp) :: covariance_residual, off_diagonal
      integer :: ik, endpoint_plus, endpoint_minus, iw, mirror
      real(rp), parameter :: c = 0.8_rp, s = 0.6_rp

      weights = 1.0_rp
      omega = [-0.31_rp, 0.0_rp, 0.31_rp]
      eigenvectors = cmplx(0.0_rp, 0.0_rp, rp)
      do ik = 1, nk
         angle = 2.0_rp*acos(-1.0_rp)*real(ik-1, rp)/real(nk, rp)
         cosine = cos(angle)
         sine = sin(angle)
         eigenvalues(:, ik) = [-0.38_rp+0.02_rp*cosine, 0.18_rp+0.02_rp*cosine, &
            -0.18_rp+0.015_rp*sine, 0.39_rp+0.01_rp*sine]
         eigenvectors(1, 1, ik) = c
         eigenvectors(3, 1, ik) = s
         eigenvectors(1, 3, ik) = -s
         eigenvectors(3, 3, ik) = c
         eigenvectors(2, 2, ik) = c
         eigenvectors(4, 2, ik) = s
         eigenvectors(2, 4, ik) = -s
         eigenvectors(4, 4, ik) = c
      end do
      do ik = 1, nk
         endpoint_plus = mod(ik, nk)+1
         endpoint_minus = mod(ik-2+nk, nk)+1
         eigenvalues_plus(:, ik) = eigenvalues(:, endpoint_plus)
         eigenvalues_minus(:, ik) = eigenvalues(:, endpoint_minus)
         eigenvectors_plus(:, :, ik) = eigenvectors(:, :, endpoint_plus)
         eigenvectors_minus(:, :, ik) = eigenvectors(:, :, endpoint_minus)
      end do
      plus = [response_channel(1, RESPONSE_PLUS), response_channel(2, RESPONSE_PLUS)]
      minus = [response_channel(1, RESPONSE_MINUS), response_channel(2, RESPONSE_MINUS)]

      call setup_options(plus_options, TDDFT_CIRCULAR_PLUS_MINUS, [nk, 1, 1])
      plus_options%q_direct = [0.25_rp, 0.0_rp, 0.0_rp]
      minus_options = plus_options
      minus_options%circular_channel = 'minus_plus'
      minus_options%q_direct = [-0.25_rp, 0.0_rp, 0.0_rp]
      call build_chi_ks_from_eigenpairs(weights, eigenvalues, eigenvectors, eigenvalues_plus, eigenvectors_plus, [1, 1], &
         plus, minus, omega, plus_options, plus_result)
      call build_chi_ks_from_eigenpairs(weights, eigenvalues, eigenvectors, eigenvalues_minus, eigenvectors_minus, [1, 1], &
         minus, plus, omega, minus_options, minus_result)

      covariance_residual = 0.0_rp
      do iw = 1, nw
         mirror = nw-iw+1
         covariance_residual = max(covariance_residual, maxval(abs(plus_result%chi(:, :, iw)- &
            transpose(conjg(minus_result%chi(:, :, mirror))))))
      end do
      off_diagonal = maxval(abs(plus_result%chi(1, 2, :)))
      write (*, '(a,1x,es12.4)') 'OBSERVED finite_q_frequency_covariance_residual', covariance_residual
      write (*, '(a,1x,es12.4)') 'OBSERVED finite_q_rank2_offdiagonal', off_diagonal
      call assert_true('finite +q/-q and +/-omega covariance', covariance_residual < 1.0e-11_rp)
      call assert_true('finite-q fixture exercises rank greater than one', off_diagonal > 1.0e-6_rp)
      call assert_true('finite-q plus-minus metadata is actual channel', &
         trim(plus_result%metadata%circular_channel) == 'plus_minus')
      call assert_true('finite-q minus-plus metadata is actual channel', &
         trim(minus_result%metadata%circular_channel) == 'minus_plus')
      call assert_same_occupation('finite-q channels share occupation metadata', plus_result, minus_result)
   end subroutine test_finite_q_frequency_covariance

   subroutine test_global_spin_reversal()
      real(rp) :: weights(1), eigenvalues(2, 1), omega(3), covariance_residual
      complex(rp) :: eigenvectors(2, 2, 1), eigenvectors_reversed(2, 2, 1)
      type(response_channel) :: plus(1), minus(1)
      type(tddft_chi0_options) :: plus_options, minus_options
      type(tddft_chi0_result) :: plus_original, minus_original, plus_reversed, minus_reversed

      weights = 1.0_rp
      eigenvalues(:, 1) = [-0.20_rp, 0.20_rp]
      omega = [-0.17_rp, 0.0_rp, 0.17_rp]
      eigenvectors = cmplx(0.0_rp, 0.0_rp, rp)
      eigenvectors(1, 1, 1) = 1.0_rp
      eigenvectors(2, 2, 1) = 1.0_rp
      eigenvectors_reversed = cmplx(0.0_rp, 0.0_rp, rp)
      eigenvectors_reversed(2, 1, 1) = 1.0_rp
      eigenvectors_reversed(1, 2, 1) = 1.0_rp
      plus = [response_channel(1, RESPONSE_PLUS)]
      minus = [response_channel(1, RESPONSE_MINUS)]

      call setup_options(plus_options, TDDFT_CIRCULAR_PLUS_MINUS, [1, 1, 1])
      minus_options = plus_options
      minus_options%circular_channel = 'minus_plus'
      call build_chi_ks_from_eigenpairs(weights, eigenvalues, eigenvectors, eigenvalues, eigenvectors, [1], plus, minus, &
         omega, plus_options, plus_original)
      call build_chi_ks_from_eigenpairs(weights, eigenvalues, eigenvectors, eigenvalues, eigenvectors, [1], minus, plus, &
         omega, minus_options, minus_original)
      call build_chi_ks_from_eigenpairs(weights, eigenvalues, eigenvectors_reversed, eigenvalues, eigenvectors_reversed, [1], &
         plus, minus, omega, plus_options, plus_reversed)
      call build_chi_ks_from_eigenpairs(weights, eigenvalues, eigenvectors_reversed, eigenvalues, eigenvectors_reversed, [1], &
         minus, plus, omega, minus_options, minus_reversed)

      covariance_residual = max(maxval(abs(plus_original%chi-minus_reversed%chi)), &
         maxval(abs(minus_original%chi-plus_reversed%chi)))
      write (*, '(a,1x,es12.4)') 'OBSERVED global_spin_reversal_channel_swap_residual', covariance_residual
      call assert_true('global spin reversal swaps circular channels', covariance_residual < 1.0e-12_rp)
      call assert_same_occupation('spin-reversed channels share occupation metadata', plus_original, minus_reversed)
   end subroutine test_global_spin_reversal

   subroutine setup_options(options, channel, k_mesh_shape)
      type(tddft_chi0_options), intent(out) :: options
      integer, intent(in) :: channel, k_mesh_shape(3)
      type(tddft_response_occupation_state) :: state

      call state%set(0.0_rp, 0.0_rp, 0.0_rp, 1, 0, 'TDCOV-03 fixture', 'fixed Fermi fixture')
      call install_tddft_occupation(options, state)
      options%eta = 0.002_rp
      options%k_mesh_shape = k_mesh_shape
      options%response_projection = 'site'
      options%circular_channel = trim(circular_channel_name(channel))
      options%use_batched_accumulation = .true.
      options%transition_batch_size = 8
   end subroutine setup_options

   subroutine assert_same_occupation(label, actual, expected)
      character(len=*), intent(in) :: label
      type(tddft_chi0_result), intent(in) :: actual, expected
      logical :: same

      same = actual%metadata%fermi_level == expected%metadata%fermi_level .and. &
         actual%metadata%electronic_temperature == expected%metadata%electronic_temperature .and. &
         actual%metadata%electronic_kT == expected%metadata%electronic_kT .and. &
         actual%metadata%band_first == expected%metadata%band_first .and. &
         actual%metadata%band_last == expected%metadata%band_last .and. &
         actual%metadata%occupation_prune_tolerance == expected%metadata%occupation_prune_tolerance .and. &
         trim(actual%metadata%fermi_source) == trim(expected%metadata%fermi_source) .and. &
         trim(actual%metadata%fermi_policy) == trim(expected%metadata%fermi_policy)
      call assert_true(label, same)
   end subroutine assert_same_occupation

   subroutine assert_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition

      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine assert_true

end program test_tddft_circular_covariance
