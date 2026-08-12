program main

#ifdef USE_MPI
   use mpi
#endif
   use mpi_mod
   use calculation_mod
   use os_mod
   use precision_mod, only: rp
   use logger_mod, only: g_logger
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc, safe_alloc
#endif
   use timer_mod, only: g_timer, timer

   implicit none
   type(calculation) :: calculation_obj
   type(argument_parser) :: args

   integer :: nomp = 1

#ifdef OpenMP_Fortran_FOUND
   ! External functions
   integer, external :: omp_get_num_threads
   integer, external :: omp_get_thread_num
#endif

#ifdef USE_MPI
   ! Initialize MPI
   call MPI_INIT(ierr)
#endif

   ! The context initializes deterministic serial state when MPI is absent and
   ! is the one synchronization point for legacy mpi_mod compatibility globals.
   g_parallel_context = parallel_context()
#ifdef USE_MPI
   if (g_parallel_context%is_root()) print *, 'Running with', g_parallel_context%size, 'MPI processes.'
#endif

#ifdef OpenMP_Fortran_FOUND
   !$omp parallel
   !!!!$omp master
   !!!nomp = omp_get_num_threads()
   !!!!$omp end master
   if (omp_get_thread_num() == 0) then
      nomp = omp_get_num_threads()
   end if
   !$omp end parallel
#endif

#ifdef USE_MPI
   if (g_parallel_context%is_root()) then
      print *, 'Each MPI process is using', nomp, 'OpenMP threads.'
      print *, ' OpenMP in play with ', nomp, ' cores.'
   end if
#endif

#ifdef USE_SAFE_ALLOC
   g_safe_alloc = safe_alloc()
#endif
   g_timer = timer()
   !call g_logger%debug(´Initializing with DEBUG=ON´, __FILE__, __LINE__)

   ! Input
   args = argument_parser()
   calculation_obj = calculation(args%input)

   ! MPI lookup table generation here
   !call get_mpi_variables(rank, <some-info> )

   ! Run
   call g_timer%start('Calculation')
   call calculation_obj%process
   call g_timer%stop('Calculation')

#ifdef USE_MPI
   if (rank == 0) then
#endif
      call g_timer%print_report()
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%print_report()
      write (*, '(1x,2(i0,a),i0)') g_safe_alloc%get_allocations_total(), ' allocations and ', g_safe_alloc%get_deallocations_total(), ' deallocations, remaining memory(B):', g_safe_alloc%get_remaining_memory()
#endif
#ifdef USE_MPI
   end if
#endif

#ifdef USE_MPI
   call g_parallel_context%restore_to_default()
   ! Finalize MPI
   call MPI_FINALIZE(ierr)
#endif
end program main
