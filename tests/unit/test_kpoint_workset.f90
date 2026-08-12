! RF-03 ownership and reciprocal-coordinate contract tests.
program test_kpoint_workset
   use precision_mod, only: rp
   use mpi_mod, only: parallel_context
   use kpoint_workset_mod, only: kpoint_workset, make_kpoint_workset, make_replicated_kpoint_workset
   use logger_mod, only: g_logger
   implicit none

   real(rp), parameter :: points(3, 5) = reshape([ &
      -0.5_rp, 0.0_rp, 0.0_rp,  0.5_rp, 0.0_rp, 0.0_rp, &
      -1.25_rp, 0.2_rp, 0.0_rp,  0.25_rp, 0.2_rp, 0.0_rp, &
       1.25_rp, 0.2_rp, 0.0_rp], [3, 5])
   real(rp), parameter :: weights(5) = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp, 5.0_rp]
   type(parallel_context) :: context
   type(kpoint_workset) :: workset, shifted
   integer :: worker_size, worker_rank, ik, total_owned, isize
   integer, parameter :: worker_sizes(3) = [1, 2, 4]
   logical :: failed

   call g_logger%init()
   failed = .false.

   ! The simulated 1/2/4-rank matrix is independent of the test launch size.
   do isize = 1, size(worker_sizes)
      worker_size = worker_sizes(isize)
      total_owned = 0
      do worker_rank = 0, worker_size-1
         context%rank = worker_rank; context%size = worker_size
         workset = make_kpoint_workset(points, weights, context, .true.)
         call check('distributed count', workset%nk_local >= 0, failed)
         call check('distributed map shape', size(workset%global_to_local) == 5, failed)
         do ik = 1, workset%nk_local
            call check('local/global round trip', workset%local_index(workset%global_index(ik)) == ik, failed)
            call check('owned weight preserved', workset%weights(ik) == weights(workset%global_index(ik)), failed)
         end do
         total_owned = total_owned + workset%nk_local
      end do
      call check('distributed coverage', total_owned == 5, failed)
   end do

   context%rank = 3; context%size = 4
   workset = make_kpoint_workset(points(:, 1:2), weights(1:2), context, .true.)
   call check('zero-work rank count', workset%nk_local == 0, failed)
   call check('zero-work range', workset%global_start == 1 .and. workset%global_end == 0, failed)
   call check('zero-work allocated arrays', size(workset%points, 2) == 0 .and. size(workset%weights) == 0, failed)

   context%rank = 0; context%size = 4
   workset = make_replicated_kpoint_workset(points, weights, context)
   call check('replicated ownership', .not. workset%distributed .and. workset%nk_local == 5, failed)
   call check('replicated raw weight sum', abs(workset%weight_sum() - sum(weights)) < 1.0e-14_rp, failed)
   call workset%fold()
   call check('boundary folds to negative edge', workset%points(1, 2) == -0.5_rp, failed)
   call check('negative coordinate folds', workset%points(1, 3) == -0.25_rp, failed)
   call check('k plus G folds', workset%points(1, 5) == 0.25_rp, failed)
   call check('folded duplicates preserved in caller order', all(workset%points(:, 4) == workset%points(:, 5)), failed)
   shifted = workset%shifted([0.5_rp, 0.0_rp, 0.0_rp])
   call check('q shift leaves base untouched', workset%points(1, 1) == -0.5_rp, failed)
   call check('q shift follows fold convention', shifted%points(1, 1) == 0.0_rp, failed)
   call check('q shift preserves weights', all(shifted%weights == workset%weights), failed)

   if (failed) error stop 1
   write (*, '(a)') 'RESULT: PASS'
contains
   subroutine check(label, condition, failed)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      logical, intent(inout) :: failed
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAILED:', trim(label)
         failed = .true.
      end if
   end subroutine check
end program test_kpoint_workset
