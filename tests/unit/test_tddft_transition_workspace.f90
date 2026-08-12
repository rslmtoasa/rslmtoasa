!------------------------------------------------------------------------------
! GC-01 -- TD-DFT transition workspace lifecycle
!------------------------------------------------------------------------------
program test_tddft_transition_workspace
   use precision_mod, only: rp
   use response_components_mod, only: RESPONSE_MINUS, RESPONSE_PLUS
   use response_vertices_mod, only: response_channel
   use tddft_transition_engine_mod, only: tddft_transition_workspace, tddft_transition_engine, tddft_vertex_provider, &
      site_channel_vertex_provider, pair_operator_vertex_provider, make_site_channel_vertex_provider, &
      make_pair_operator_vertex_provider
   implicit none

   type(tddft_transition_workspace) :: workspace
   logical :: failed

   failed = .false.
   call check_workspace('default workspace is unallocated', workspace, 0, 0, 0, 0, .false.)

   call workspace%ensure_capacity(2, 3, 5, 4)
   call check_workspace('first use', workspace, 2, 3, 5, 4, .true.)

   workspace%band_n = 17
   call workspace%ensure_capacity(2, 3, 5, 4)
   call check_workspace('same-shape reuse', workspace, 2, 3, 5, 4, .true.)
   call check_true('same-shape reuse retains workspace storage', all(workspace%band_n == 17))

   call workspace%ensure_capacity(2, 3, 5, 6)
   call check_workspace('changed capacity', workspace, 2, 3, 5, 6, .true.)

   call workspace%ensure_capacity(4, 1, 2, 6)
   call check_workspace('changed response dimensions', workspace, 4, 1, 2, 6, .true.)

   deallocate(workspace%right_vertices)
   call workspace%ensure_capacity(4, 1, 2, 6)
   call check_workspace('partial workspace is rebuilt', workspace, 4, 1, 2, 6, .true.)

   call workspace%clear()
   call check_workspace('clear', workspace, 0, 0, 0, 0, .false.)

   call workspace%ensure_capacity(3, 2, 4, 5)
   call check_workspace('reuse after clear', workspace, 3, 2, 4, 5, .true.)

   call test_engine_reuse()

   if (failed) error stop 1
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_engine_reuse()
      integer, target :: site_orbital_counts(1)
      type(response_channel), target :: left_channels(1), right_channels(1)
      real(rp) :: k_weights(2), eigenvalues(2,2)
      complex(rp), target :: eigenvectors(2,2,2), operators(2,2,1)
      type(tddft_transition_engine) :: site_engine, pair_engine
      type(site_channel_vertex_provider) :: site_provider
      type(pair_operator_vertex_provider) :: pair_provider

      site_orbital_counts = [1]
      left_channels(1) = response_channel(1, RESPONSE_PLUS)
      right_channels(1) = response_channel(1, RESPONSE_MINUS)
      k_weights = [1.0_rp, 3.0_rp]
      eigenvalues(:,1) = [-0.10_rp, 0.10_rp]
      eigenvalues(:,2) = [-0.08_rp, 0.12_rp]
      eigenvectors = cmplx(0.0_rp, 0.0_rp, rp)
      eigenvectors(1,1,:) = cmplx(1.0_rp, 0.0_rp, rp)
      eigenvectors(2,2,:) = cmplx(1.0_rp, 0.0_rp, rp)
      operators = cmplx(0.0_rp, 0.0_rp, rp)
      operators(1,1,1) = cmplx(-0.03_rp, 0.01_rp, rp)
      operators(2,2,1) = cmplx(-0.02_rp, -0.01_rp, rp)

      call make_site_channel_vertex_provider(site_provider, site_orbital_counts, left_channels, right_channels, &
         eigenvectors, eigenvectors)
      call check_provider_reuse('site-channel chi_KS', site_engine, site_provider, k_weights, eigenvalues)

      call make_pair_operator_vertex_provider(pair_provider, site_orbital_counts, left_channels, eigenvectors, eigenvectors, operators)
      call check_provider_reuse('pair-operator Xi', pair_engine, pair_provider, k_weights, eigenvalues)
   end subroutine test_engine_reuse

   subroutine check_provider_reuse(label, engine, provider, k_weights, eigenvalues)
      character(len=*), intent(in) :: label
      class(tddft_transition_engine), intent(inout) :: engine
      class(tddft_vertex_provider), intent(in) :: provider
      real(rp), intent(in) :: k_weights(:), eigenvalues(:, :)
      complex(rp) :: dynamic(1,1,2), dynamic_repeat(1,1,2), static(1,1,1), static_repeat(1,1,1)
      real(rp) :: vertex_seconds, preparation_seconds, denominator_seconds, accumulation_seconds
      integer :: allocation_count

      dynamic = cmplx(0.0_rp, 0.0_rp, rp)
      vertex_seconds = 0.0_rp; preparation_seconds = 0.0_rp
      denominator_seconds = 0.0_rp; accumulation_seconds = 0.0_rp
      call engine%accumulate_dynamic(k_weights, eigenvalues, eigenvalues, [0.03_rp, 0.08_rp], 0.004_rp, 0.0_rp, 300.0_rp, &
         1, 2, 0.0_rp, 3, .true., provider, dynamic, vertex_seconds, preparation_seconds, denominator_seconds, accumulation_seconds)
      allocation_count = engine%workspace%storage_allocations
      call check_true(trim(label)//': coefficient scratch prepared', allocated(engine%workspace%bra) .and. &
         allocated(engine%workspace%ket) .and. all(shape(engine%workspace%bra) == [2, 3]))

      dynamic_repeat = cmplx(0.0_rp, 0.0_rp, rp)
      vertex_seconds = 0.0_rp; preparation_seconds = 0.0_rp
      denominator_seconds = 0.0_rp; accumulation_seconds = 0.0_rp
      call engine%accumulate_dynamic(k_weights, eigenvalues, eigenvalues, [0.03_rp, 0.08_rp], 0.004_rp, 0.0_rp, 300.0_rp, &
         1, 2, 0.0_rp, 3, .true., provider, dynamic_repeat, vertex_seconds, preparation_seconds, denominator_seconds, accumulation_seconds)
      call check_true(trim(label)//': repeated dynamic call reuses workspace', &
         engine%workspace%storage_allocations == allocation_count)
      call check_true(trim(label)//': repeated dynamic result is stable', &
         maxval(abs(dynamic_repeat-dynamic)) <= 1024.0_rp*epsilon(1.0_rp))

      static = cmplx(0.0_rp, 0.0_rp, rp)
      vertex_seconds = 0.0_rp; preparation_seconds = 0.0_rp; accumulation_seconds = 0.0_rp
      call engine%accumulate_static(k_weights, eigenvalues, 0.0_rp, 300.0_rp, 1, 2, 0.0_rp, 3, provider, static, &
         vertex_seconds, preparation_seconds, accumulation_seconds)
      static_repeat = cmplx(0.0_rp, 0.0_rp, rp)
      vertex_seconds = 0.0_rp; preparation_seconds = 0.0_rp; accumulation_seconds = 0.0_rp
      call engine%accumulate_static(k_weights, eigenvalues, 0.0_rp, 300.0_rp, 1, 2, 0.0_rp, 3, provider, static_repeat, &
         vertex_seconds, preparation_seconds, accumulation_seconds)
      call check_true(trim(label)//': repeated static call reuses workspace', &
         engine%workspace%storage_allocations == allocation_count .and. engine%workspace%capacity_reuses >= 3)
      call check_true(trim(label)//': repeated static result is stable', &
         maxval(abs(static_repeat-static)) <= 1024.0_rp*epsilon(1.0_rp))
   end subroutine check_provider_reuse

   subroutine check_workspace(label, value, nleft, nright, ncoefficient, capacity, expected_allocated)
      character(len=*), intent(in) :: label
      type(tddft_transition_workspace), intent(in) :: value
      integer, intent(in) :: nleft, nright, ncoefficient, capacity
      logical, intent(in) :: expected_allocated
      logical :: allocated_all

      allocated_all = allocated(value%band_n) .and. allocated(value%band_m) .and. allocated(value%occupations) .and. &
         allocated(value%transition_energies) .and. allocated(value%left_vertices) .and. allocated(value%right_vertices) .and. &
         allocated(value%denominators) .and. allocated(value%weighted_left) .and. allocated(value%bra) .and. allocated(value%ket)
      call check_true(trim(label)//': allocation state', allocated_all .eqv. expected_allocated)
      call check_true(trim(label)//': recorded capacity', value%capacity == capacity)
      call check_true(trim(label)//': recorded coefficient dimension', value%coefficient_dimension == ncoefficient)
      if (.not. expected_allocated) return
      call check_true(trim(label)//': band_n shape', size(value%band_n) == capacity)
      call check_true(trim(label)//': band_m shape', size(value%band_m) == capacity)
      call check_true(trim(label)//': occupations shape', size(value%occupations) == capacity)
      call check_true(trim(label)//': transition energy shape', size(value%transition_energies) == capacity)
      call check_true(trim(label)//': left vertex shape', all(shape(value%left_vertices) == [nleft, capacity]))
      call check_true(trim(label)//': right vertex shape', all(shape(value%right_vertices) == [nright, capacity]))
      call check_true(trim(label)//': denominator shape', size(value%denominators) == capacity)
      call check_true(trim(label)//': weighted left shape', all(shape(value%weighted_left) == [nleft, capacity]))
      call check_true(trim(label)//': bra shape', all(shape(value%bra) == [ncoefficient, capacity]))
      call check_true(trim(label)//': ket shape', all(shape(value%ket) == [ncoefficient, capacity]))
   end subroutine check_workspace

   subroutine check_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      if (.not. condition) then
         failed = .true.
         write (*, '(a)') 'FAIL: '//trim(label)
      end if
   end subroutine check_true

end program test_tddft_transition_workspace
