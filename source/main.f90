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

#ifdef OpenMP_Fortran_FOUND
  ! External functions
  integer, external :: omp_get_num_threads
  integer :: nomp
#endif

#ifdef USE_MPI
  ! Initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
  if (rank == 0) then
    print *, 'Running with', numprocs, 'MPI processes.'
  end if
#endif

#ifdef OpenMP_Fortran_FOUND
  !$omp parallel
  !$omp master
  nomp=omp_get_num_threads()
  !$omp end master
  !$omp end parallel
#endif

#ifdef USE_MPI
  if (rank == 0) then
    print *, 'Each MPI process is using', nomp, 'OpenMP threads.'
    print *,' OpenMP in play with ' , nomp,' cores.'
  end if
#endif

#ifdef USE_SAFE_ALLOC
  g_safe_alloc = safe_alloc()
#endif
  g_timer = timer()
  !call g_logger%debug('Initializing with DEBUG=ON', __FILE__, __LINE__)

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
  write(*,'(1x,2(i0,a),i0)') g_safe_alloc%get_allocations_total(), ' allocations and ', g_safe_alloc%get_deallocations_total(), ' deallocations, remaining memory(B):', g_safe_alloc%get_remaining_memory()
#endif
#ifdef USE_MPI
  endif
#endif

#ifdef USE_MPI
  ! Finalize MPI
  call MPI_FINALIZE(ierr)
#endif
end program main
