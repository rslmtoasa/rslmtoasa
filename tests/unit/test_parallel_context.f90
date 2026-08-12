!------------------------------------------------------------------------------
! RF-02 -- deterministic parallel context, range splitting, and device mapping
!------------------------------------------------------------------------------
program test_parallel_context
   use mpi_mod, only: parallel_context, assignment(=), rank, numprocs, ierr, get_mpi_range, &
      split_range_for_worker, map_local_rank_to_device
#ifdef USE_MPI
   use mpi
#endif
   implicit none

   type(parallel_context) :: context
   integer :: first, last, count, device, global_count
   integer, allocatable :: l2g(:), g2l(:)
   logical :: valid, failed

#ifdef USE_MPI
   call MPI_INIT(ierr)
#endif
   context = parallel_context()
   failed = .false.

   call check_true('compatibility rank is synchronized', rank == context%rank, failed)
   call check_true('compatibility size is synchronized', numprocs == context%size, failed)
   call check_true('root predicate', context%is_root() .eqv. (context%rank == 0), failed)

   ! Explicit worker splitting is side-effect free and covers serial, no-work,
   ! and more-worker-than-item cases independently of this launch size.
   call split_range_for_worker(0, 1, 0, first, last, count)
   call check_range('serial zero items', first, last, count, 1, 0, 0, failed)
   call split_range_for_worker(0, 4, 2, first, last, count)
   call check_range('rank zero split', first, last, count, 1, 1, 1, failed)
   call split_range_for_worker(1, 4, 2, first, last, count)
   call check_range('rank one split', first, last, count, 2, 2, 1, failed)
   call split_range_for_worker(2, 4, 2, first, last, count)
   call check_range('oversubscribed rank split', first, last, count, 1, 0, 0, failed)
   call split_range_for_worker(0, 0, 2, first, last, count)
   call check_range('invalid worker count', first, last, count, 1, 0, 0, failed)

   call context%split_range(2, first, last, count)
   call check_true('context range count is nonnegative', count >= 0, failed)
   call check_true('context owns first local item', .not. context%owns_work(1, 2) .or. first == 1, failed)
   call get_mpi_range(rank, 0, first, last, count, l2g, g2l, 'RF-02 zero-work')
   call check_true('legacy zero-work range', first == 1 .and. last == 0 .and. count == 0, failed)
   call check_true('legacy zero-work maps', size(l2g) == 0 .and. size(g2l) == 0, failed)

   call map_local_rank_to_device(0, 4, device, valid)
   call check_true('four ranks/four devices rank zero', valid .and. device == 0, failed)
   call map_local_rank_to_device(3, 4, device, valid)
   call check_true('four ranks/four devices rank three', valid .and. device == 3, failed)
   call map_local_rank_to_device(3, 2, device, valid)
   call check_true('four ranks/two devices wraps', valid .and. device == 1, failed)
   call map_local_rank_to_device(1, 2, device, valid, override=0)
   call check_true('programmatic device override', valid .and. device == 0, failed)
   call map_local_rank_to_device(1, 2, device, valid, override=2)
   call check_true('invalid device override rejected', .not. valid .and. device == -1, failed)
   call map_local_rank_to_device(0, 0, device, valid)
   call check_true('zero device count rejected', .not. valid .and. device == -1, failed)

   call context%split_range(2, first, last, count)
#ifdef USE_MPI
   call MPI_ALLREDUCE(count, global_count, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   call check_true('collective range coverage', global_count == 2, failed)
#else
   global_count = count
   call check_true('serial range coverage', global_count == 2, failed)
#endif
   call context%barrier()

#ifdef USE_MPI
   call reduce_failure(failed)
#endif
   if (context%is_root()) then
      if (failed) then
         write (*, '(a)') 'RESULT: FAIL'
      else
         write (*, '(a)') 'RESULT: PASS'
      end if
   end if

   call context%restore_to_default()
#ifdef USE_MPI
   call MPI_FINALIZE(ierr)
#endif
   if (failed) error stop 1

contains

   subroutine check_range(label, actual_first, actual_last, actual_count, expected_first, expected_last, expected_count, test_failed)
      character(len=*), intent(in) :: label
      integer, intent(in) :: actual_first, actual_last, actual_count
      integer, intent(in) :: expected_first, expected_last, expected_count
      logical, intent(inout) :: test_failed

      if (actual_first /= expected_first .or. actual_last /= expected_last .or. actual_count /= expected_count) then
         if (rank == 0) write (*, '(a,1x,a)') 'FAIL range:', trim(label)
         test_failed = .true.
      end if
   end subroutine check_range

   subroutine check_true(label, condition, test_failed)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      logical, intent(inout) :: test_failed

      if (.not. condition) then
         if (rank == 0) write (*, '(a,1x,a)') 'FAIL:', trim(label)
         test_failed = .true.
      end if
   end subroutine check_true

#ifdef USE_MPI
   subroutine reduce_failure(test_failed)
      logical, intent(inout) :: test_failed
      integer :: local_failure, global_failure

      local_failure = merge(1, 0, test_failed)
      call MPI_ALLREDUCE(local_failure, global_failure, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      test_failed = global_failure /= 0
   end subroutine reduce_failure
#endif

end program test_parallel_context
