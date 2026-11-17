!------------------------------------------------------------------------------
! ACC-08 -- exercise the existing RS CUDA context through MPI rank/device mapping
!------------------------------------------------------------------------------
program test_cuda_mpi_device_mapping
   use, intrinsic :: iso_c_binding, only: c_associated, c_int
   use logger_mod, only: g_logger
   use mpi_mod, only: parallel_context, g_parallel_context, ierr
   use rsrec_cuda_plugin_mod, only: rsrec_cuda_backend, get_gpu_context
#ifdef USE_MPI
   use mpi
#endif
   implicit none

   type(rsrec_cuda_backend), pointer :: gpu
   logical :: failed
   integer(c_int) :: device_count, status

   interface
      function cuda_device_count(count) bind(C, name='rsrec_cuda_device_count') result(status)
         import c_int
         integer(c_int), intent(out) :: count
         integer(c_int) :: status
      end function cuda_device_count
   end interface

#ifdef USE_MPI
   call MPI_INIT(ierr)
#endif
   call g_logger%init()
   g_parallel_context = parallel_context()
   status = cuda_device_count(device_count)
   if (status /= 0_c_int .or. device_count <= 0_c_int) then
      if (g_parallel_context%is_root()) write (*, '(a)') 'SKIP: no CUDA device is available'
#ifdef USE_MPI
      call MPI_FINALIZE(ierr)
#endif
      stop 77
   end if
   failed = .false.

   gpu => get_gpu_context(1, 1, 1, 1, 0)
   if (.not. c_associated(gpu%ctx)) then
      failed = .true.
      if (g_parallel_context%is_root()) write (*, '(a)') 'FAIL: CUDA context was not created'
   end if
   call gpu%destroy()

#ifdef USE_MPI
   call reduce_failure(failed)
#endif
   if (g_parallel_context%is_root()) then
      if (failed) then
         write (*, '(a)') 'RESULT: FAIL'
      else
         write (*, '(a)') 'RESULT: PASS'
      end if
   end if

#ifdef USE_MPI
   call MPI_FINALIZE(ierr)
#endif
   if (failed) error stop 1

contains

#ifdef USE_MPI
   subroutine reduce_failure(local_failed)
      logical, intent(inout) :: local_failed
      integer :: local_flag, global_flag

      local_flag = merge(1, 0, local_failed)
      call MPI_ALLREDUCE(local_flag, global_flag, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      local_failed = global_flag /= 0
   end subroutine reduce_failure
#endif

end program test_cuda_mpi_device_mapping
