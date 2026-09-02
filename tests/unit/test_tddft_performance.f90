! TDDFT-14 -- decomposition and profile contract.
program test_tddft_performance
   use, intrinsic :: iso_fortran_env, only: int64
   use precision_mod, only: rp
   use tddft_performance_mod, only: tddft_mpi_plan, tddft_performance_profile, tddft_work_range, &
      make_tddft_mpi_plan, split_tddft_work, tddft_checksum_complex
#ifdef USE_MPI
   use mpi
   use mpi_mod, only: ierr
#endif
   implicit none

   type(tddft_mpi_plan) :: plan
   type(tddft_performance_profile) :: profile
   type(tddft_work_range) :: work
   complex(rp) :: values(2)
   logical :: failed
   integer :: i, total
#ifdef USE_MPI
   integer :: mpi_rank
   integer :: mpi_size
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)
#endif

   failed = .false.
   call split_tddft_work(0, 4, 10, work)
   call check_true('balanced split gives first range', work%first == 1 .and. work%last == 3 .and. work%count == 3, failed)
   total = 0
   do i = 0, 3
      call split_tddft_work(i, 4, 10, work)
      total = total + work%count
   end do
   call check_true('balanced split covers every work item', total == 10, failed)
   call split_tddft_work(3, 4, 2, work)
   call check_true('oversubscribed split is empty and safe', work%count == 0 .and. work%first == 1 .and. work%last == 0, failed)

   plan = make_tddft_mpi_plan('realspace_gf', 'direct', 12, 64, 256, 48, 0, 4, .true., 1024)
   call check_true('R-GF computes the complete q batch', plan%q%first == 1 .and. plan%q%last == 12 .and. &
      plan%owner_q%count == 3, failed)
   call check_true('R-GF preserves q amortization', plan%preserves_realspace_reuse .and. plan%q_fourier_is_batched .and. &
      plan%requires_collective_reduction, failed)
   call check_true('R-GF partitions R blocks first', index(trim(plan%strategy), 'R-blocks') > 0 .and. plan%r%count == 12, failed)
#ifdef USE_MPI
   if (mpi_size > 1) then
      plan = make_tddft_mpi_plan('realspace_gf', 'direct', 12, 64, 256, 48, mpi_rank, mpi_size, .true., 1024)
      call check_true('MPI R-GF rank owns a valid output range', plan%owner_q%first >= 1 .and. &
         plan%owner_q%last <= 12 .and. plan%owner_q%count > 0, failed)
      call check_true('MPI R-GF rank keeps the shared q batch', plan%q%first == 1 .and. plan%q%last == 12 .and. &
         plan%preserves_realspace_reuse .and. plan%requires_collective_reduction, failed)
   end if
#endif

   plan = make_tddft_mpi_plan('eigenpairs', 'not_applicable', 12, 64, 256, 0, 2, 4, .true.)
   call check_true('reciprocal backend owns a q range', plan%q%first == 7 .and. plan%q%last == 9 .and. &
      trim(plan%strategy) == 'q-outer / omega-inner', failed)
   plan = make_tddft_mpi_plan('kspace_lehmann', 'direct', 1, 64, 256, 0, 2, 4, .true.)
   call check_true('frequency fallback keeps the q point complete', plan%q%count == 1 .and. plan%omega%count == 16 .and. &
      index(trim(plan%strategy), 'omega-outer') > 0, failed)

   values = [cmplx(1.0_rp, 2.0_rp, rp), cmplx(3.0_rp, 4.0_rp, rp)]
   call check_true('scientific checksum is deterministic', abs(tddft_checksum_complex(values)-23.0_rp) < 1.0e-12_rp, failed)
   call profile%configure(.true., 'eigenpairs', 'transitions')
   call profile%set_workload(128_int64, 0_int64, 1)
   call profile%set_memory(1000_int64, 2000_int64, 3000_int64, 6000_int64)
   call profile%set_checksum(tddft_checksum_complex(values))
   call profile%add_seconds('H_G_construction', 1.0e-3_rp)
   call profile%add_seconds('GF_contraction', 2.0e-3_rp)
   call profile%add_seconds('energy_integration', 3.0e-3_rp)
   call profile%add_seconds('R_q_FT', 4.0e-3_rp)
   call profile%add_seconds('Dyson', 5.0e-3_rp)
   call profile%add_seconds('spectral_analysis', 6.0e-3_rp)
   call profile%add_seconds('MPI_reduction', 7.0e-3_rp)
   call profile%emit(label='contract')

#ifdef USE_MPI
   call reduce_failure(failed)
   if (mpi_rank == 0) then
#endif
      if (failed) then
         write (*, '(a)') 'RESULT: FAIL'
      else
         write (*, '(a)') 'RESULT: PASS'
      end if
#ifdef USE_MPI
   end if
   call MPI_FINALIZE(ierr)
#endif
   if (failed) error stop 1

contains

   subroutine check_true(label, condition, test_failed)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      logical, intent(inout) :: test_failed
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL:', trim(label)
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

end program test_tddft_performance
