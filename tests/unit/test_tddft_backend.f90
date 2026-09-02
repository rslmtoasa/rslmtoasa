!------------------------------------------------------------------------------
! TDDFT-04 -- common chi0 backend contract
!------------------------------------------------------------------------------
module tddft_backend_test_fixture_mod
   use precision_mod, only: rp
   use tddft_chi0_mod, only: tddft_chi0_request, tddft_chi0_batch_result
   use tddft_backend_mod, only: tddft_realspace_chi0_provider, tddft_backend_capabilities
   implicit none

   type, extends(tddft_realspace_chi0_provider) :: fixture_realspace_provider
   contains
      procedure :: evaluate_realspace => evaluate_fixture_realspace
      procedure :: describe => describe_fixture_realspace
   end type fixture_realspace_provider

contains

   subroutine evaluate_fixture_realspace(this, request, result)
      class(fixture_realspace_provider), intent(inout) :: this
      type(tddft_chi0_request), intent(in) :: request
      type(tddft_chi0_batch_result), intent(out) :: result
      integer :: iq, iw

      allocate(result%q_response(size(request%q_points, 2)), result%q_points(3, size(request%q_points, 2)), &
         result%q_indices(size(request%q_points, 2)), result%omega(size(request%omega)))
      result%q_points = request%q_points; result%omega = request%omega
      do iq = 1, size(request%q_points, 2)
         result%q_indices(iq) = iq
         allocate(result%q_response(iq)%chi(1, 1, size(request%omega)), &
            result%q_response(iq)%re_chi(1, 1, size(request%omega)), result%q_response(iq)%im_chi(1, 1, size(request%omega)))
         do iw = 1, size(request%omega)
            result%q_response(iq)%chi(1, 1, iw) = cmplx(real(iq, rp), -request%omega(iw), rp)
         end do
         result%q_response(iq)%re_chi = real(result%q_response(iq)%chi, rp)
         result%q_response(iq)%im_chi = aimag(result%q_response(iq)%chi)
         result%q_response(iq)%metadata%backend = 'fixture_realspace'
         result%q_response(iq)%metadata%canonical_backend = 'realspace_gf'
         result%q_response(iq)%metadata%implementation = 'fixture chi0(R,w) provider'
         result%q_response(iq)%metadata%real_space_reuse = .true.
      end do
   end subroutine evaluate_fixture_realspace

   subroutine describe_fixture_realspace(this, metadata)
      class(fixture_realspace_provider), intent(in) :: this
      type(tddft_backend_capabilities), intent(out) :: metadata

      metadata = tddft_backend_capabilities()
      metadata%native_real_space = .true.
      metadata%reuses_real_space_response = .true.
   end subroutine describe_fixture_realspace

end module tddft_backend_test_fixture_mod

program test_tddft_backend
   use precision_mod, only: rp
   use response_components_mod, only: RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, tddft_chi0_request, tddft_chi0_batch_result
   use tddft_chi0_green_mod, only: green_chi0_options
   use tddft_backend_test_fixture_mod, only: fixture_realspace_provider
   use tddft_backend_mod, only: tddft_chi0_backend, tddft_eigenpair_backend, tddft_kspace_lehmann_backend, &
      tddft_realspace_gf_backend, tddft_mock_chi0_backend, tddft_backend_capabilities, &
      tddft_backend_selection, canonical_tddft_backend_name, select_tddft_chi0_backend, make_tddft_chi0_backend
   use tddft_dyson_mod, only: tddft_dyson_options, tddft_dyson_result, enhance_tddft_susceptibility
   implicit none

   logical :: failed

   failed = .false.
   call test_selector_and_factory()
   call test_mock_batches_feed_common_dyson()
   call test_eigenpair_adapter_preserves_batch_shapes()
   call test_kspace_lehmann_adapter_uses_common_result()
   call test_native_realspace_provider_is_not_kspace_adapted()
   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_selector_and_factory()
      class(tddft_chi0_backend), allocatable :: backend
      type(tddft_backend_capabilities) :: capabilities
      type(tddft_backend_selection) :: selection

      call check_true('eigenpair name is canonical', trim(canonical_tddft_backend_name('eigenpairs')) == 'eigenpairs')
      call check_true('green remains a compatibility alias', trim(canonical_tddft_backend_name('green')) == 'kspace_lehmann')
      call check_true('native real-space name is explicit', trim(canonical_tddft_backend_name('rs_gf')) == 'realspace_gf')
      selection = select_tddft_chi0_backend('green')
      call check_true('selector records compatibility alias', selection%compatibility_alias)
      call check_true('selector retains requested spelling', trim(selection%requested_name) == 'green')

      call make_tddft_chi0_backend('eigenpairs', backend)
      call check_true('factory creates eigenpair contract', trim(backend%backend_name()) == 'eigenpairs')
      call backend%capabilities(capabilities)
      call check_true('eigenpair contract supports q batches', capabilities%supports_q_batch)
      deallocate(backend)

      call make_tddft_chi0_backend('kspace_lehmann', backend)
      call check_true('factory creates K-space Lehmann contract', trim(backend%backend_name()) == 'kspace_lehmann')
      deallocate(backend)

      call make_tddft_chi0_backend('realspace_gf', backend)
      call check_true('factory creates native real-space contract', trim(backend%backend_name()) == 'realspace_gf')
      call backend%capabilities(capabilities)
      call check_true('real-space capability advertises native reuse', capabilities%native_real_space .and. &
         capabilities%reuses_real_space_response)
      deallocate(backend)
   end subroutine test_selector_and_factory

   subroutine test_mock_batches_feed_common_dyson()
      type(tddft_mock_chi0_backend) :: backend
      type(tddft_chi0_batch_result) :: batch
      type(tddft_dyson_options) :: dyson_options
      type(tddft_dyson_result) :: dyson
      real(rp) :: q_points(3, 3), omega(4)
      complex(rp) :: kernel(2, 2)
      integer :: q_indices(3)

      q_points = reshape([0.0_rp, 0.0_rp, 0.0_rp, 0.1_rp, 0.0_rp, 0.0_rp, &
         0.2_rp, 0.0_rp, 0.0_rp], [3, 3])
      omega = [0.02_rp, 0.04_rp, 0.06_rp, 0.08_rp]
      q_indices = [7, 8, 9]
      kernel = cmplx(0.0_rp, 0.0_rp, rp)
      call backend%initialize(2)
      call backend%evaluate_grid(q_points, omega, batch, q_indices)
      call check_true('mock grid returns every q point', size(batch%q_response) == 3)
      call check_true('mock grid retains frequency batch', size(batch%q_response(1)%chi, 3) == 4)
      call check_true('mock provenance reports q and omega batch sizes', batch%metadata%q_batch_size == 3 .and. &
         batch%metadata%omega_batch_size == 4)
      call check_true('mock provenance reports completed fixture', batch%metadata%converged)

      ! The common Dyson layer receives only the ordinary chi0 result.  It has
      ! no branch for, or knowledge of, the backend that produced it.
      call enhance_tddft_susceptibility(batch%q_response(2)%chi, kernel, 0.001_rp, dyson_options, dyson)
      call check_true('mock response feeds backend-agnostic Dyson', maxval(abs(dyson%chi- &
         batch%q_response(2)%chi)) < 1.0e-12_rp)
   end subroutine test_mock_batches_feed_common_dyson

   subroutine test_eigenpair_adapter_preserves_batch_shapes()
      class(tddft_chi0_backend), allocatable :: backend
      type(tddft_chi0_options) :: options
      type(tddft_chi0_result) :: frequency_result
      type(tddft_chi0_batch_result) :: q_result
      type(response_channel) :: left(1), right(1)
      real(rp) :: weights(1), eval(2, 1), evalq(2, 1, 2), q_points(3, 2), omega(2)
      complex(rp) :: evec(2, 2, 1), evecq(2, 2, 1, 2)
      integer :: q_indices(2)

      weights = 1.0_rp
      eval(:, 1) = [-0.4_rp, 0.4_rp]
      evalq(:, 1, 1) = eval(:, 1); evalq(:, 1, 2) = eval(:, 1) + [0.03_rp, -0.02_rp]
      evec = cmplx(0.0_rp, 0.0_rp, rp)
      evec(1, 1, 1) = cmplx(1.0_rp, 0.0_rp, rp); evec(2, 2, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      evecq = 0.0_rp
      evecq(:, :, :, 1) = evec; evecq(:, :, :, 2) = evec
      q_points = reshape([0.0_rp, 0.0_rp, 0.0_rp, 0.25_rp, 0.0_rp, 0.0_rp], [3, 2])
      q_indices = [1, 2]
      omega = [0.05_rp, 0.10_rp]
      left(1) = response_channel(1, RESPONSE_PLUS)
      right(1) = response_channel(1, RESPONSE_MINUS)
      options%eta = 0.001_rp
      options%fermi_level = 0.0_rp
      options%k_mesh_shape = [1, 1, 1]

      call make_tddft_chi0_backend('eigenpairs', backend)
      select type (backend)
      type is (tddft_eigenpair_backend)
         call backend%initialize(weights, eval, evec, evalq, evecq, q_points, [1], left, right, options)
         call backend%evaluate_frequency_batch(q_points(:, 1), omega, frequency_result, q_index=1)
         call backend%evaluate_q_batch(q_points, omega(1), q_result, q_indices)
      class default
         call check_true('factory dynamic type is eigenpair adapter', .false.)
         return
      end select
      call check_true('eigenpair frequency batch retains both frequencies', size(frequency_result%chi, 3) == 2)
      call check_true('eigenpair q batch returns both endpoint responses', size(q_result%q_response) == 2)
      call check_true('eigenpair q batch returns endpoint indices', all(q_result%q_indices == q_indices))
      call check_true('eigenpair q batch metadata identifies canonical backend', &
         trim(q_result%metadata%canonical_backend) == 'eigenpairs')
      call check_true('eigenpair q batch preserves endpoint provenance', &
         index(trim(q_result%q_response(2)%metadata%endpoint_provenance), 'k and k+q') > 0)
      deallocate(backend)
   end subroutine test_eigenpair_adapter_preserves_batch_shapes

   subroutine test_kspace_lehmann_adapter_uses_common_result()
      class(tddft_chi0_backend), allocatable :: backend
      type(tddft_chi0_options) :: options
      type(green_chi0_options) :: green_options
      type(tddft_chi0_result) :: result
      type(response_channel) :: left(1), right(1)
      real(rp) :: weights(1), eval(2, 1), evalq(2, 1, 1), q_points(3, 1), omega(1)
      complex(rp) :: evec(2, 2, 1), evecq(2, 2, 1, 1)

      weights = 1.0_rp
      eval(:, 1) = [-0.4_rp, 0.4_rp]; evalq(:, 1, 1) = eval(:, 1)
      evec = cmplx(0.0_rp, 0.0_rp, rp)
      evec(1, 1, 1) = cmplx(1.0_rp, 0.0_rp, rp); evec(2, 2, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      evecq(:, :, :, 1) = evec
      q_points = 0.0_rp; omega = 0.05_rp
      left(1) = response_channel(1, RESPONSE_PLUS); right(1) = response_channel(1, RESPONSE_MINUS)
      options%eta = 0.001_rp; options%fermi_level = 0.0_rp
      green_options%eta = options%eta; green_options%green_eta = 0.0005_rp
      green_options%fermi_level = options%fermi_level
      green_options%energy_min = -0.6_rp; green_options%energy_max = 0.6_rp
      green_options%energy_points = 201

      call make_tddft_chi0_backend('kspace_lehmann', backend)
      select type (backend)
      type is (tddft_kspace_lehmann_backend)
         call backend%initialize(weights, eval, evec, evalq, evecq, q_points, [1], left, right, options, green_options)
         call backend%evaluate_frequency_batch(q_points(:, 1), omega, result, q_index=1)
      class default
         call check_true('factory dynamic type is K-space Lehmann adapter', .false.)
         return
      end select
      call check_true('K-space Lehmann adapter returns common chi0 result', allocated(result%chi))
      call check_true('K-space Lehmann adapter returns canonical provenance', &
         trim(result%metadata%canonical_backend) == 'kspace_lehmann')
      deallocate(backend)
   end subroutine test_kspace_lehmann_adapter_uses_common_result

   subroutine test_native_realspace_provider_is_not_kspace_adapted()
      type(fixture_realspace_provider) :: provider
      class(tddft_chi0_backend), allocatable :: backend
      type(tddft_chi0_batch_result) :: result
      type(tddft_backend_capabilities) :: capabilities
      real(rp) :: q_points(3, 2), omega(2)

      q_points = reshape([0.0_rp, 0.0_rp, 0.0_rp, 0.25_rp, 0.0_rp, 0.0_rp], [3, 2])
      omega = [0.05_rp, 0.10_rp]
      call make_tddft_chi0_backend('realspace_gf', backend)
      select type (backend)
      type is (tddft_realspace_gf_backend)
         call backend%capabilities(capabilities)
         call check_true('native real-space backend advertises real-space execution', capabilities%native_real_space)
         call check_true('uninitialized native backend does not advertise static support', &
            .not. capabilities%supports_static_limit)
         call backend%initialize(provider)
         call backend%capabilities(capabilities)
         call check_true('static support follows the attached provider capability', &
            .not. capabilities%supports_static_limit)
         call backend%evaluate_grid(q_points, omega, result)
      class default
         call check_true('factory dynamic type is native real-space adapter', .false.)
         return
      end select
      call check_true('native provider returns a q batch without k-space endpoints', size(result%q_response) == 2)
      call check_true('native provider provenance records reuse', result%metadata%real_space_reuse)
      call check_true('native provider provenance remains native real-space', &
         trim(result%metadata%canonical_backend) == 'realspace_gf')
      deallocate(backend)
   end subroutine test_native_realspace_provider_is_not_kspace_adapted

   subroutine check_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine check_true

end program test_tddft_backend
