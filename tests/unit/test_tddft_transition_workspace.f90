!------------------------------------------------------------------------------
! GC-01 -- TD-DFT transition workspace lifecycle
!------------------------------------------------------------------------------
program test_tddft_transition_workspace
   use precision_mod, only: rp
   use tddft_transition_engine_mod, only: tddft_transition_workspace
   implicit none

   type(tddft_transition_workspace) :: workspace
   logical :: failed

   failed = .false.
   call check_workspace('default workspace is unallocated', workspace, 0, 0, 0, .false.)

   call workspace%ensure_capacity(2, 3, 4)
   call check_workspace('first use', workspace, 2, 3, 4, .true.)

   workspace%band_n = 17
   call workspace%ensure_capacity(2, 3, 4)
   call check_workspace('same-shape reuse', workspace, 2, 3, 4, .true.)
   call check_true('same-shape reuse retains workspace storage', all(workspace%band_n == 17))

   call workspace%ensure_capacity(2, 3, 6)
   call check_workspace('changed capacity', workspace, 2, 3, 6, .true.)

   call workspace%ensure_capacity(4, 1, 6)
   call check_workspace('changed response dimensions', workspace, 4, 1, 6, .true.)

   deallocate(workspace%right_vertices)
   call workspace%ensure_capacity(4, 1, 6)
   call check_workspace('partial workspace is rebuilt', workspace, 4, 1, 6, .true.)

   call workspace%clear()
   call check_workspace('clear', workspace, 0, 0, 0, .false.)

   call workspace%ensure_capacity(3, 2, 5)
   call check_workspace('reuse after clear', workspace, 3, 2, 5, .true.)

   if (failed) error stop 1
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine check_workspace(label, value, nleft, nright, capacity, expected_allocated)
      character(len=*), intent(in) :: label
      type(tddft_transition_workspace), intent(in) :: value
      integer, intent(in) :: nleft, nright, capacity
      logical, intent(in) :: expected_allocated
      logical :: allocated_all

      allocated_all = allocated(value%band_n) .and. allocated(value%band_m) .and. allocated(value%occupations) .and. &
         allocated(value%transition_energies) .and. allocated(value%left_vertices) .and. allocated(value%right_vertices) .and. &
         allocated(value%denominators) .and. allocated(value%weighted_left)
      call check_true(trim(label)//': allocation state', allocated_all .eqv. expected_allocated)
      call check_true(trim(label)//': recorded capacity', value%capacity == capacity)
      if (.not. expected_allocated) return
      call check_true(trim(label)//': band_n shape', size(value%band_n) == capacity)
      call check_true(trim(label)//': band_m shape', size(value%band_m) == capacity)
      call check_true(trim(label)//': occupations shape', size(value%occupations) == capacity)
      call check_true(trim(label)//': transition energy shape', size(value%transition_energies) == capacity)
      call check_true(trim(label)//': left vertex shape', all(shape(value%left_vertices) == [nleft, capacity]))
      call check_true(trim(label)//': right vertex shape', all(shape(value%right_vertices) == [nright, capacity]))
      call check_true(trim(label)//': denominator shape', size(value%denominators) == capacity)
      call check_true(trim(label)//': weighted left shape', all(shape(value%weighted_left) == [nleft, capacity]))
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
