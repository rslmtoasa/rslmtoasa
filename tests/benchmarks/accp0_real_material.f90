! ACC-P0 persistent real-material benchmark driver.
!
! This executable deliberately owns the benchmark lifecycle.  A single
! process reads a production input, builds the production lattice/charge/
! Hamiltonian/reciprocal stack, creates one execution backend, performs the
! requested warm-ups, resets interval counters, and then measures repetitions.
! It is opt-in and is not registered as a default CTest test.
program accp0_real_material
   use, intrinsic :: iso_c_binding, only: c_float, c_long_long
   use, intrinsic :: iso_fortran_env, only: int64
   use precision_mod, only: rp
   use basis_mod, only: basis_init, nb
   use math_mod, only: init_math_operators
   use mpi_mod, only: parallel_context, g_parallel_context
   use timer_mod, only: timer, g_timer
   use reciprocal_mod, only: reciprocal, reciprocal_execution_request, reciprocal_execution_result, &
      cuda_reciprocal_backend
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use control_mod, only: control
   use logger_mod, only: g_logger
   implicit none

   character(len=64) :: mode, backend, solver_strategy, fixture, workload
   character(len=128) :: source_name
   character(len=512) :: input_file, dump_file
   integer :: mesh(3), tile_size, eigenvectors, warmups, repetitions, supercell_l
   integer :: status

   mode = 'persistent'
   backend = 'lapack'
   solver_strategy = 'zheevd_serial'
   fixture = 'unknown'
   source_name = 'production'
   workload = 'crossover'
   input_file = 'input.nml'
   dump_file = ''
   mesh = [1, 1, 1]
   tile_size = 16
   eigenvectors = 0
   warmups = 2
   repetitions = 5
   supercell_l = 1

   call parse_arguments(mode, backend, solver_strategy, fixture, source_name, workload, input_file, dump_file, mesh, tile_size, &
      eigenvectors, warmups, repetitions, supercell_l)
   if (any(mesh < 1) .or. tile_size < 1 .or. eigenvectors < 0 .or. eigenvectors > 1 .or. warmups < 0 .or. repetitions < 1) then
      error stop 'ACCP0: invalid benchmark dimensions or repetition policy'
   end if
   if (trim(mode) /= 'persistent' .and. trim(mode) /= 'oracle' .and. trim(mode) /= 'preflight' .and. trim(mode) /= 'validate') then
      error stop 'ACCP0: mode must be persistent, oracle, preflight, or validate'
   end if
   if (trim(backend) /= 'lapack' .and. trim(backend) /= 'cuda') then
      error stop 'ACCP0: backend must be lapack or cuda'
   end if
   if (trim(backend) == 'cuda' .and. trim(solver_strategy) == 'lapack') then
      error stop 'ACCP1: CUDA backend requires zheevd_serial or zheevj_batched'
   end if
   call g_parallel_context%restore_to_default()
   g_parallel_context = parallel_context()
   g_timer = timer()
   call init_math_operators()
   call g_logger%init()

   if (trim(mode) == 'oracle') then
      call run_oracle(trim(input_file), trim(dump_file), trim(fixture), trim(backend), trim(solver_strategy), mesh, supercell_l)
   else if (trim(mode) == 'preflight') then
      call run_preflight(trim(input_file), trim(fixture), trim(backend), trim(solver_strategy), mesh, tile_size, eigenvectors == 1, supercell_l)
   else if (trim(mode) == 'validate') then
      call run_validation(trim(input_file), trim(fixture), trim(backend), trim(solver_strategy), mesh, tile_size, supercell_l)
   else
      call run_persistent(trim(input_file), trim(fixture), trim(source_name), trim(workload), trim(backend), trim(solver_strategy), mesh, &
         tile_size, eigenvectors == 1, warmups, repetitions, supercell_l)
   end if

contains

   subroutine parse_arguments(mode, backend, solver_strategy, fixture, source_name, workload, input_file, dump_file, mesh, tile_size, &
                              eigenvectors, warmups, repetitions, supercell_l)
      character(len=*), intent(inout) :: mode, backend, solver_strategy, fixture, source_name, workload, input_file, dump_file
      integer, intent(inout) :: mesh(3), tile_size, eigenvectors, warmups, repetitions, supercell_l
      integer :: i, narg
      character(len=128) :: key, value

      narg = command_argument_count()
      i = 1
      do while (i <= narg)
         call get_command_argument(i, key)
         if (i == narg) error stop 'ACCP0: option requires a value'
         call get_command_argument(i + 1, value)
         select case (trim(key))
         case ('--mode'); mode = trim(value)
         case ('--backend'); backend = trim(value)
         case ('--solver-strategy'); solver_strategy = trim(value)
         case ('--fixture'); fixture = trim(value)
         case ('--source'); source_name = trim(value)
         case ('--workload'); workload = trim(value)
         case ('--input'); input_file = trim(value)
         case ('--dump-eigenvalues'); dump_file = trim(value)
         case ('--mesh'); call parse_mesh(value, mesh)
         case ('--tile-size'); read(value, *) tile_size
         case ('--eigenvectors'); read(value, *) eigenvectors
         case ('--warmups'); read(value, *) warmups
         case ('--repetitions'); read(value, *) repetitions
         case ('--L'); read(value, *) supercell_l
         case default; error stop 'ACCP0: unknown option '//trim(key)
         end select
         i = i + 2
      end do
   end subroutine parse_arguments

   subroutine parse_mesh(text, mesh)
      character(len=*), intent(in) :: text
      integer, intent(out) :: mesh(3)
      character(len=128) :: normalized
      integer :: i

      normalized = text
      do i = 1, len_trim(normalized)
         if (normalized(i:i) == 'x' .or. normalized(i:i) == 'X' .or. normalized(i:i) == ',') normalized(i:i) = ' '
      end do
      read(normalized, *) mesh
   end subroutine parse_mesh

   subroutine setup_production(input_file, recip, ham, lat, chg, ctl)
      character(len=*), intent(in) :: input_file
      type(reciprocal), intent(out) :: recip
      type(hamiltonian), target, intent(out) :: ham
      type(lattice), target, intent(out) :: lat
      type(charge), target, intent(out) :: chg
      type(control), target, intent(out) :: ctl
      integer :: i

      ctl = control(input_file)
      lat = lattice(ctl)
      call lat%build_data()
      call lat%bravais()
      call lat%structb(.true.)
      call lat%atomlist()
      call basis_init(lat%symbolic_atoms(1)%potential%lmax)
      chg = charge(lat)
      call chg%bulkmat()
      ham = hamiltonian(chg)
      do i = 1, lat%nrec
         call lat%symbolic_atoms(i)%build_pot()
      end do
      if (ctl%nsp == 2 .or. ctl%nsp == 4) call ham%build_lsham()
      call ham%build_bulkham()
      recip = reciprocal(ham)
   end subroutine setup_production

   subroutine make_mesh(mesh, points)
      integer, intent(in) :: mesh(3)
      real(rp), allocatable, intent(out) :: points(:, :)
      integer :: i1, i2, i3, ik

      allocate(points(3, mesh(1) * mesh(2) * mesh(3)))
      ik = 0
      do i3 = 0, mesh(3) - 1
         do i2 = 0, mesh(2) - 1
            do i1 = 0, mesh(1) - 1
               ik = ik + 1
               points(:, ik) = [real(i1, rp) / real(mesh(1), rp), real(i2, rp) / real(mesh(2), rp), &
                  real(i3, rp) / real(mesh(3), rp)]
            end do
         end do
      end do
   end subroutine make_mesh

   subroutine assemble_all(recip, points, hamiltonians)
      type(reciprocal), intent(inout) :: recip
      real(rp), intent(in) :: points(:, :)
      complex(rp), intent(out) :: hamiltonians(:, :, :)
      integer :: ik

      do ik = 1, size(points, 2)
         call recip%build_hamiltonian_at_kpoint(points(:, ik), hamiltonians(:, :, ik))
      end do
   end subroutine assemble_all

   subroutine solve_tiles(recip, hamiltonians, tile_size, want_vectors, status, dump_unit)
      type(reciprocal), intent(inout) :: recip
      complex(rp), intent(in) :: hamiltonians(:, :, :)
      integer, intent(in) :: tile_size
      logical, intent(in) :: want_vectors
      integer, intent(out) :: status
      integer, intent(in), optional :: dump_unit
      type(reciprocal_execution_request) :: request
      type(reciprocal_execution_result) :: result
      integer :: first, last, length, ik, i

      status = 0
      do first = 1, size(hamiltonians, 3), tile_size
         last = min(size(hamiltonians, 3), first + tile_size - 1)
         length = last - first + 1
         request%assemble_hamiltonian = .false.
         request%assemble_overlap = .false.
         request%solve_eigensystem = .true.
         request%generalized = .false.
         request%request_eigenvectors = want_vectors
         request%operator_generation = 1
         if (allocated(request%input_hamiltonian)) deallocate(request%input_hamiltonian)
         allocate(request%input_hamiltonian(size(hamiltonians, 1), size(hamiltonians, 2), length), &
            source=hamiltonians(:, :, first:last))
         call recip%execution_backend%execute_batch(request, result)
         if (.not. result%eigenvalues_valid) status = 1
         if (result%local_point_count /= length) status = 1
         if (want_vectors .and. .not. result%eigenvectors_valid) status = 1
         if ((.not. want_vectors) .and. result%eigenvectors_valid) status = 1
         if (present(dump_unit)) then
            if (.not. result%eigenvalues_valid) status = 1
            do ik = 1, length
               do i = 1, size(result%eigenvalues, 1)
                  write(dump_unit, '(es26.17)') result%eigenvalues(i, ik)
               end do
            end do
         end if
         if (allocated(result%eigenvalues)) deallocate(result%eigenvalues)
         if (allocated(result%eigenvectors)) deallocate(result%eigenvectors)
      end do
   end subroutine solve_tiles

   subroutine run_oracle(input_file, dump_file, fixture, backend_name, solver_strategy, mesh, supercell_l)
      character(len=*), intent(in) :: input_file, dump_file, fixture, backend_name, solver_strategy
      integer, intent(in) :: mesh(3), supercell_l
      type(reciprocal) :: recip
      type(hamiltonian), target :: ham
      type(lattice), target :: lat
      type(charge), target :: chg
      type(control), target :: ctl
      real(rp), allocatable :: points(:, :)
      complex(rp), allocatable :: hamiltonians(:, :, :)
      integer :: dump_unit, status, unique_nk

      if (len_trim(dump_file) == 0) error stop 'ACCP0: oracle requires --dump-eigenvalues'
      call setup_production(input_file, recip, ham, lat, chg, ctl)
      call make_mesh(mesh, points)
      allocate(hamiltonians(nb * lat%nrec, nb * lat%nrec, size(points, 2)))
      call assemble_all(recip, points, hamiltonians)
      call recip%make_execution_backend(backend_name)
      call configure_solver_strategy(recip, backend_name, solver_strategy, status)
      if (status /= 0) error stop 'ACCP1: failed to configure oracle solver strategy'
      call recip%execution_backend%prepare_operator(ham%operator_generation)
      open (newunit=dump_unit, file=dump_file, status='replace', action='write')
      call solve_tiles(recip, hamiltonians, size(points, 2), .false., status, dump_unit)
      close(dump_unit)
      if (status /= 0) error stop 'ACCP0: oracle eigensolution failed'
      unique_nk = count_unique_points(recip, points)
      write (*, '(a)') 'ACCP0_ORACLE fixture='//trim(fixture)//' L='//trim(int_token(supercell_l))// &
         ' backend='//trim(backend_name)//' nmat='//trim(int_token(nb * lat%nrec))// &
         ' solver_strategy='//trim(solver_strategy)//' actual_unique_nk='//trim(int_token(unique_nk))
   end subroutine run_oracle

   subroutine run_persistent(input_file, fixture, source_name, workload, backend_name, solver_strategy, mesh, tile_size, want_vectors, &
                             warmups, repetitions, supercell_l)
      character(len=*), intent(in) :: input_file, fixture, source_name, workload, backend_name, solver_strategy
      integer, intent(in) :: mesh(3), tile_size, warmups, repetitions, supercell_l
      logical, intent(in) :: want_vectors
      type(reciprocal) :: recip
      type(hamiltonian), target :: ham
      type(lattice), target :: lat
      type(charge), target :: chg
      type(control), target :: ctl
      real(rp), allocatable :: points(:, :), steady_times(:), solve_times(:)
      complex(rp), allocatable :: hamiltonians(:, :, :)
      real(rp) :: start, finish, hk_start, hk_finish, solve_start, solve_finish
      real(rp) :: backend_init, prepare_seconds, first_solve, hk_cpu, total_median, solve_median
      real(rp) :: steady_min, steady_max, steady_spread, conversion, staging, h2d, solver, d2h, d2h_values, d2h_vectors
      real(rp) :: sync, widening, total_reciprocal, other_backend, post_seconds
      real(rp) :: memory_estimate, free_memory, total_memory
      integer :: status, i, nk, unique_nk, nmat, strategy_status
      logical :: supported
      character(len=128) :: unsupported_reason
      integer(c_long_long) :: free_bytes, total_bytes
      integer(c_long_long) :: h2d_bytes, d2h_values_bytes, d2h_vectors_bytes
      integer(c_long_long) :: cuda_malloc_count, cuda_free_count, workspace_query_count, workspace_reuse_count
      integer(c_long_long) :: event_create_count, event_destroy_count, pinned_alloc_count, pinned_free_count
      integer(c_long_long) :: cuda_malloc_before, cuda_free_before, workspace_query_before, workspace_reuse_before
      integer(c_long_long) :: event_create_before, event_destroy_before, pinned_alloc_before, pinned_free_before
      integer :: pinned_host_active
      character(len=2048) :: line

      call setup_production(input_file, recip, ham, lat, chg, ctl)
      write (*, '(a)') 'ACCP1_PRECISION cpu_solver=LAPACK_ZHEEV_ZHEGV_complex128 cuda_solver=cuDoubleComplex_eigenvalue_float64_'// &
         'cusolverDnZheevd_or_ZheevjBatched fp32_route_explicit_only=1 normal_physics_precision=complex128_real64'
      call make_mesh(mesh, points)
      nk = size(points, 2)
      unique_nk = count_unique_points(recip, points)
      nmat = nb * lat%nrec
      call wall_clock(start)
      call recip%make_execution_backend(backend_name)
      call wall_clock(finish)
      backend_init = finish - start
      call configure_solver_strategy(recip, backend_name, solver_strategy, strategy_status)
      if (strategy_status /= 0) error stop 'ACCP1: failed to configure CUDA solver strategy'

      call wall_clock(start)
      call recip%execution_backend%prepare_operator(ham%operator_generation)
      call wall_clock(finish)
      prepare_seconds = finish - start

      call backend_memory(recip, free_bytes, total_bytes)
      free_memory = real(free_bytes, rp) / (1024.0_rp * 1024.0_rp)
      total_memory = real(total_bytes, rp) / (1024.0_rp * 1024.0_rp)
      memory_estimate = estimate_memory_mib(nmat, tile_size, want_vectors, solver_strategy)
      if (trim(backend_name) == 'cuda' .and. total_memory > 0.0_rp .and. free_memory < 2.0_rp * memory_estimate) then
         error stop 'ACCP0: CUDA memory preflight rejected this workload'
      end if
      call backend_solver_supported(recip, nmat, tile_size, want_vectors, supported, unsupported_reason)
      if (.not. supported) then
         line = 'ACCP0_DIMENSIONS'
         call append_token(line, 'fixture='//trim(fixture))
         call append_token(line, 'source='//trim(source_name))
         call append_token(line, 'workload='//trim(workload))
         call append_token(line, 'backend='//trim(backend_name))
         call append_token(line, 'solver_strategy='//trim(solver_strategy))
         call append_token(line, 'L='//trim(int_token(supercell_l)))
         call append_token(line, 'natom='//trim(int_token(lat%nrec)))
         call append_token(line, 'nmat='//trim(int_token(nmat)))
         call append_token(line, 'nominal_mesh='//trim(mesh_token(mesh)))
         call append_token(line, 'actual_unique_nk='//trim(int_token(unique_nk)))
         call append_token(line, 'tile='//trim(int_token(tile_size)))
         call append_token(line, 'eigenvectors='//trim(int_token(merge(1, 0, want_vectors))))
         write (*, '(a)') trim(line)
         line = 'ACCP0_TIMING'
         call append_token(line, 'fixture='//trim(fixture))
         call append_token(line, 'backend='//trim(backend_name))
         call append_token(line, 'solver_strategy='//trim(solver_strategy))
         call append_token(line, 'support_status=unsupported')
         call append_token(line, 'unsupported_reason='//trim(underscored(unsupported_reason)))
         call append_real(line, 'cold_process_wall_s', 0.0_rp)
         call append_real(line, 'cuda_context_backend_init_s', backend_init)
         call append_real(line, 'operator_prepare_s', 0.0_rp)
         call append_real(line, 'first_solve_s', 0.0_rp)
         call append_real(line, 'steady_solve_median_s', 0.0_rp)
         call append_real(line, 'steady_solve_min_s', 0.0_rp)
         call append_real(line, 'steady_solve_spread_s', 0.0_rp)
         call append_token(line, 'metric_repetitions=0')
         call append_real(line, 'Hk_CPU_s', 0.0_rp)
         call append_real(line, 'H64_to_H32_s', 0.0_rp)
         call append_real(line, 'T_host_staging_s', 0.0_rp)
         call append_real(line, 'H2D_s', 0.0_rp)
         call append_real(line, 'solver_s', 0.0_rp)
         call append_real(line, 'D2H_s', 0.0_rp)
         call append_real(line, 'T_D2H_values_s', 0.0_rp)
         call append_real(line, 'T_D2H_vectors_s', 0.0_rp)
         call append_real(line, 'T_sync_s', 0.0_rp)
         call append_real(line, 'H32_to_H64_s', 0.0_rp)
         call append_real(line, 'T_other_backend_s', 0.0_rp)
         call append_real(line, 'T_total_s', 0.0_rp)
         call append_real(line, 'total_reciprocal_s', 0.0_rp)
         call append_token(line, 'H2D_bytes=0')
         call append_token(line, 'D2H_values_bytes=0')
         call append_token(line, 'D2H_vectors_bytes=0')
         call append_token(line, 'pinned_host_active=0')
         call append_token(line, 'cuda_malloc_count=0')
         call append_token(line, 'cuda_free_count=0')
         call append_token(line, 'workspace_query_count=0')
         call append_token(line, 'workspace_reuse_count=0')
         call append_token(line, 'event_create_count=0')
         call append_token(line, 'event_destroy_count=0')
         call append_token(line, 'pinned_alloc_count=0')
         call append_token(line, 'pinned_free_count=0')
         call append_real(line, 'post_s', 0.0_rp)
         call append_real(line, 'total_steady_s', 0.0_rp)
         call append_real(line, 'memory_estimate_mib', memory_estimate)
         call append_real(line, 'memory_free_before_mib', free_memory)
         call append_real(line, 'memory_total_mib', total_memory)
         write (*, '(a)') trim(line)
         write (*, '(a)') 'RESULT: UNSUPPORTED'
         return
      end if

      allocate(hamiltonians(nmat, nmat, nk), steady_times(repetitions), solve_times(repetitions))
      call wall_clock(start)
      call assemble_all(recip, points, hamiltonians)
      call wall_clock(finish)
      hk_cpu = finish - start

      call wall_clock(solve_start)
      call solve_tiles(recip, hamiltonians, tile_size, want_vectors, status)
      call wall_clock(solve_finish)
      first_solve = solve_finish - solve_start
      if (status /= 0) error stop 'ACCP0: first eigensolution failed'

      do i = 1, warmups
         call solve_tiles(recip, hamiltonians, tile_size, want_vectors, status)
         if (status /= 0) error stop 'ACCP0: warm-up eigensolution failed'
      end do
      call backend_resource_counters(recip, cuda_malloc_before, cuda_free_before, workspace_query_before, workspace_reuse_before, &
         event_create_before, event_destroy_before, pinned_alloc_before, pinned_free_before)
      call reset_backend_metrics(recip)

      do i = 1, repetitions
         call wall_clock(start)
         call wall_clock(hk_start)
         call assemble_all(recip, points, hamiltonians)
         call wall_clock(hk_finish)
         hk_cpu = hk_cpu + (hk_finish - hk_start)
         call wall_clock(solve_start)
         call solve_tiles(recip, hamiltonians, tile_size, want_vectors, status)
         call wall_clock(solve_finish)
         if (status /= 0) error stop 'ACCP0: measured eigensolution failed'
         solve_times(i) = solve_finish - solve_start
         call wall_clock(finish)
         steady_times(i) = finish - start
      end do
      hk_cpu = hk_cpu / real(repetitions + 1, rp)

      total_median = median(steady_times)
      solve_median = median(solve_times)
      steady_min = minval(steady_times)
      steady_max = maxval(steady_times)
      steady_spread = steady_max - steady_min
      call backend_timings(recip, conversion, staging, h2d, solver, d2h, d2h_values, d2h_vectors, sync, widening, &
         total_reciprocal, h2d_bytes, d2h_values_bytes, d2h_vectors_bytes, pinned_host_active, cuda_malloc_count, &
         cuda_free_count, workspace_query_count, workspace_reuse_count, event_create_count, event_destroy_count, &
         pinned_alloc_count, pinned_free_count)
      ! CUDA counters are reset immediately before this measured interval and
      ! accumulate over all measured repetitions.  Normalize them to the
      ! same per-repetition unit as the reported steady medians; retain the
      ! interval sample count in the machine-readable record.
      conversion = conversion / real(repetitions, rp)
      staging = staging / real(repetitions, rp)
      h2d = h2d / real(repetitions, rp)
      solver = solver / real(repetitions, rp)
      d2h = d2h / real(repetitions, rp)
      d2h_values = d2h_values / real(repetitions, rp)
      d2h_vectors = d2h_vectors / real(repetitions, rp)
      sync = sync / real(repetitions, rp)
      widening = widening / real(repetitions, rp)
      total_reciprocal = total_reciprocal / real(repetitions, rp)
      post_seconds = max(0.0_rp, total_median - hk_cpu - solve_median)
      other_backend = max(0.0_rp, solve_median - total_reciprocal)

      line = 'ACCP0_DIMENSIONS'
      call append_token(line, 'fixture='//trim(fixture))
      call append_token(line, 'source='//trim(source_name))
      call append_token(line, 'workload='//trim(workload))
      call append_token(line, 'backend='//trim(backend_name))
      call append_token(line, 'solver_strategy='//trim(solver_strategy))
      call append_token(line, 'L='//trim(int_token(supercell_l)))
      call append_token(line, 'natom='//trim(int_token(lat%nrec)))
      call append_token(line, 'nmat='//trim(int_token(nmat)))
      call append_token(line, 'nominal_mesh='//trim(mesh_token(mesh)))
      call append_token(line, 'actual_unique_nk='//trim(int_token(unique_nk)))
      call append_token(line, 'tile='//trim(int_token(tile_size)))
      call append_token(line, 'eigenvectors='//trim(int_token(merge(1, 0, want_vectors))))
      write (*, '(a)') trim(line)

      line = 'ACCP0_TIMING'
      call append_token(line, 'fixture='//trim(fixture))
      call append_token(line, 'backend='//trim(backend_name))
      call append_token(line, 'solver_strategy='//trim(solver_strategy))
      call append_token(line, 'support_status=supported')
      call append_real(line, 'cold_process_wall_s', 0.0_rp)
      call append_real(line, 'cuda_context_backend_init_s', backend_init)
      call append_real(line, 'operator_prepare_s', prepare_seconds)
      call append_real(line, 'first_solve_s', first_solve)
      call append_real(line, 'steady_solve_median_s', solve_median)
      call append_real(line, 'steady_solve_min_s', minval(solve_times))
      call append_real(line, 'steady_solve_spread_s', maxval(solve_times) - minval(solve_times))
      call append_token(line, 'metric_repetitions='//trim(int_token(repetitions)))
      call append_real(line, 'Hk_CPU_s', hk_cpu)
      call append_real(line, 'H64_to_H32_s', conversion)
      call append_real(line, 'T_Hk_CPU_s', hk_cpu)
      call append_real(line, 'T_host_staging_s', staging)
      call append_real(line, 'H2D_s', h2d)
      call append_real(line, 'solver_s', solver)
      call append_real(line, 'D2H_s', d2h)
      call append_real(line, 'T_D2H_values_s', d2h_values)
      call append_real(line, 'T_D2H_vectors_s', d2h_vectors)
      call append_real(line, 'T_sync_s', sync)
      call append_real(line, 'H32_to_H64_s', widening)
      call append_real(line, 'T_other_backend_s', other_backend)
      call append_real(line, 'T_total_s', total_median)
      call append_real(line, 'total_reciprocal_s', total_reciprocal)
      call append_token(line, 'H2D_bytes='//trim(int64_token(h2d_bytes / int(repetitions, c_long_long))))
      call append_token(line, 'D2H_values_bytes='//trim(int64_token(d2h_values_bytes / int(repetitions, c_long_long))))
      call append_token(line, 'D2H_vectors_bytes='//trim(int64_token(d2h_vectors_bytes / int(repetitions, c_long_long))))
      call append_token(line, 'pinned_host_active='//trim(int_token(pinned_host_active)))
      call append_token(line, 'cuda_malloc_count='//trim(int64_token(cuda_malloc_count)))
      call append_token(line, 'cuda_free_count='//trim(int64_token(cuda_free_count)))
      call append_token(line, 'workspace_query_count='//trim(int64_token(workspace_query_count)))
      call append_token(line, 'workspace_reuse_count='//trim(int64_token(workspace_reuse_count)))
      call append_token(line, 'event_create_count='//trim(int64_token(event_create_count)))
      call append_token(line, 'event_destroy_count='//trim(int64_token(event_destroy_count)))
      call append_token(line, 'pinned_alloc_count='//trim(int64_token(pinned_alloc_count)))
      call append_token(line, 'pinned_free_count='//trim(int64_token(pinned_free_count)))
      call append_token(line, 'cuda_malloc_count_before='//trim(int64_token(cuda_malloc_before)))
      call append_token(line, 'cuda_free_count_before='//trim(int64_token(cuda_free_before)))
      call append_token(line, 'workspace_query_count_before='//trim(int64_token(workspace_query_before)))
      call append_token(line, 'workspace_reuse_count_before='//trim(int64_token(workspace_reuse_before)))
      call append_token(line, 'event_create_count_before='//trim(int64_token(event_create_before)))
      call append_token(line, 'event_destroy_count_before='//trim(int64_token(event_destroy_before)))
      call append_token(line, 'pinned_alloc_count_before='//trim(int64_token(pinned_alloc_before)))
      call append_token(line, 'pinned_free_count_before='//trim(int64_token(pinned_free_before)))
      call append_real(line, 'post_s', post_seconds)
      call append_real(line, 'total_steady_s', total_median)
      call append_real(line, 'steady_total_min_s', steady_min)
      call append_real(line, 'steady_total_spread_s', steady_spread)
      call append_real(line, 'memory_estimate_mib', memory_estimate)
      call append_real(line, 'memory_free_before_mib', free_memory)
      call append_real(line, 'memory_total_mib', total_memory)
      write (*, '(a)') trim(line)
      write (*, '(a)') 'RESULT: PASS'
   end subroutine run_persistent

   ! ACC-P1 real-material correctness gate.  The CPU LAPACK result is kept as
   ! the oracle; CUDA is checked for eigenvalues, residuals, orthogonality,
   ! and projectors of degenerate eigenspaces.  Raw eigenvector phases are
   ! intentionally never compared.
   subroutine run_validation(input_file, fixture, backend_name, solver_strategy, mesh, tile_size, supercell_l)
      character(len=*), intent(in) :: input_file, fixture, backend_name, solver_strategy
      integer, intent(in) :: mesh(3), tile_size, supercell_l
      type(reciprocal) :: recip
      type(hamiltonian), target :: ham
      type(lattice), target :: lat
      type(charge), target :: chg
      type(control), target :: ctl
      real(rp), allocatable :: points(:, :), cpu_values(:, :), cuda_values(:, :)
      complex(rp), allocatable :: hamiltonians(:, :, :), cpu_vectors(:, :, :), cuda_vectors(:, :, :)
      integer :: nmat, nk, status, strategy_status
      logical :: supported
      character(len=128) :: reason
      real(rp) :: eigenvalue_error, residual_error, orthogonality_error, projector_error, hamiltonian_rounding_error
      character(len=2048) :: line

      if (trim(backend_name) /= 'cuda') error stop 'ACCP1: validation requires the CUDA backend'
      call setup_production(input_file, recip, ham, lat, chg, ctl)
      call make_mesh(mesh, points)
      nk = size(points, 2)
      nmat = nb * lat%nrec
      allocate(hamiltonians(nmat, nmat, nk))
      call assemble_all(recip, points, hamiltonians)

      call recip%make_execution_backend('lapack')
      call recip%execution_backend%prepare_operator(ham%operator_generation)
      call solve_tiles_capture(recip, hamiltonians, tile_size, cpu_values, cpu_vectors, status)
      if (status /= 0) error stop 'ACC-P1: CPU validation eigensolution failed'

      call recip%make_execution_backend('cuda')
      call configure_solver_strategy(recip, backend_name, solver_strategy, strategy_status)
      if (strategy_status /= 0) error stop 'ACC-P1: failed to configure validation solver strategy'
      call backend_solver_supported(recip, nmat, tile_size, .true., supported, reason)
      if (.not. supported) then
         line = 'ACCP1_VALIDATION'
         call append_token(line, 'fixture='//trim(fixture))
         call append_token(line, 'L='//trim(int_token(supercell_l)))
         call append_token(line, 'solver_strategy='//trim(solver_strategy))
         call append_token(line, 'nmat='//trim(int_token(nmat)))
         call append_token(line, 'nk='//trim(int_token(nk)))
         call append_token(line, 'status=UNSUPPORTED')
         call append_token(line, 'reason='//trim(underscored(reason)))
         write (*, '(a)') trim(line)
         return
      end if
      call recip%execution_backend%prepare_operator(ham%operator_generation)
      call solve_tiles_capture(recip, hamiltonians, tile_size, cuda_values, cuda_vectors, status)
      if (status /= 0) error stop 'ACC-P1: CUDA validation eigensolution failed'
      call compare_eigenpairs(hamiltonians, cpu_values, cpu_vectors, cuda_values, cuda_vectors, &
         eigenvalue_error, residual_error, orthogonality_error, projector_error, hamiltonian_rounding_error)

      line = 'ACCP1_VALIDATION'
      call append_token(line, 'fixture='//trim(fixture))
      call append_token(line, 'L='//trim(int_token(supercell_l)))
      call append_token(line, 'solver_strategy='//trim(solver_strategy))
      call append_token(line, 'nmat='//trim(int_token(nmat)))
      call append_token(line, 'nk='//trim(int_token(nk)))
      call append_real(line, 'eigenvalue_max_abs', eigenvalue_error)
      call append_real(line, 'residual_max', residual_error)
      call append_real(line, 'orthogonality_max', orthogonality_error)
      call append_real(line, 'degenerate_projector_max', projector_error)
      call append_real(line, 'H64_H32_relative_max', hamiltonian_rounding_error)
      ! Keep the established ACC-P1 numerical contract unchanged for the
      ! experimental FP32 route.  A failed FP32 result is data to classify,
      ! not a reason to loosen the production/reference tolerance.
      if (eigenvalue_error <= 1.0e-8_rp .and. residual_error <= 1.0e-8_rp .and. &
          orthogonality_error <= 1.0e-8_rp .and. projector_error <= 1.0e-7_rp) then
         call append_token(line, 'status=PASS')
      else
         call append_token(line, 'status=FAIL')
      end if
      write (*, '(a)') trim(line)
   end subroutine run_validation

   subroutine solve_tiles_capture(recip, hamiltonians, tile_size, eigenvalues, eigenvectors, status)
      type(reciprocal), intent(inout) :: recip
      complex(rp), intent(in) :: hamiltonians(:, :, :)
      integer, intent(in) :: tile_size
      real(rp), allocatable, intent(out) :: eigenvalues(:, :)
      complex(rp), allocatable, intent(out) :: eigenvectors(:, :, :)
      integer, intent(out) :: status
      type(reciprocal_execution_request) :: request
      type(reciprocal_execution_result) :: result
      integer :: nmat, nk, first, last, length

      nmat = size(hamiltonians, 1)
      nk = size(hamiltonians, 3)
      allocate(eigenvalues(nmat, nk), eigenvectors(nmat, nmat, nk))
      status = 0
      do first = 1, nk, tile_size
         last = min(nk, first + tile_size - 1)
         length = last - first + 1
         request%assemble_hamiltonian = .false.
         request%assemble_overlap = .false.
         request%solve_eigensystem = .true.
         request%generalized = .false.
         request%request_eigenvectors = .true.
         request%operator_generation = 1
         if (allocated(request%input_hamiltonian)) deallocate(request%input_hamiltonian)
         allocate(request%input_hamiltonian(nmat, nmat, length), source=hamiltonians(:, :, first:last))
         call recip%execution_backend%execute_batch(request, result)
         if (.not. result%eigenvalues_valid .or. .not. result%eigenvectors_valid .or. result%local_point_count /= length) then
            status = 1
            return
         end if
         eigenvalues(:, first:last) = result%eigenvalues
         eigenvectors(:, :, first:last) = result%eigenvectors
         if (allocated(result%eigenvalues)) deallocate(result%eigenvalues)
         if (allocated(result%eigenvectors)) deallocate(result%eigenvectors)
      end do
   end subroutine solve_tiles_capture

   subroutine compare_eigenpairs(hamiltonians, cpu_values, cpu_vectors, cuda_values, cuda_vectors, &
                                 eigenvalue_error, residual_error, orthogonality_error, projector_error, hamiltonian_rounding_error)
      complex(rp), intent(in) :: hamiltonians(:, :, :), cpu_vectors(:, :, :), cuda_vectors(:, :, :)
      real(rp), intent(in) :: cpu_values(:, :), cuda_values(:, :)
      real(rp), intent(out) :: eigenvalue_error, residual_error, orthogonality_error, projector_error, hamiltonian_rounding_error
      integer :: nmat, nk, ik, row, column, inner, group_start, group_end, group_size
      real(rp) :: hnorm, residual, group_tolerance, hdelta, hnorm2, real64, imag64
    complex(rp) :: left_value, projector_cpu, projector_cuda

      nmat = size(cpu_values, 1)
      nk = size(cpu_values, 2)
      eigenvalue_error = maxval(abs(cpu_values - cuda_values))
      residual_error = 0.0_rp
      orthogonality_error = 0.0_rp
      projector_error = 0.0_rp
      hamiltonian_rounding_error = 0.0_rp
      group_tolerance = 1.0e-8_rp
      do ik = 1, nk
         hnorm = max(1.0_rp, sqrt(sum(abs(hamiltonians(:, :, ik))**2)))
         hdelta = 0.0_rp
         hnorm2 = 0.0_rp
         do row = 1, nmat
            do inner = 1, nmat
               real64 = real(real(hamiltonians(row, inner, ik), kind=c_float), rp)
               imag64 = real(real(aimag(hamiltonians(row, inner, ik)), kind=c_float), rp)
               hdelta = hdelta + (real(hamiltonians(row, inner, ik)) - real64)**2 + &
                        (aimag(hamiltonians(row, inner, ik)) - imag64)**2
               hnorm2 = hnorm2 + abs(hamiltonians(row, inner, ik))**2
            end do
         end do
         hamiltonian_rounding_error = max(hamiltonian_rounding_error, sqrt(hdelta) / max(1.0_rp, sqrt(hnorm2)))
         do column = 1, nmat
            residual = 0.0_rp
            do row = 1, nmat
               left_value = (0.0_rp, 0.0_rp)
               do inner = 1, nmat
                  left_value = left_value + hamiltonians(row, inner, ik) * cuda_vectors(inner, column, ik)
               end do
               left_value = left_value - cuda_values(column, ik) * cuda_vectors(row, column, ik)
               residual = residual + abs(left_value)**2
            end do
            residual_error = max(residual_error, sqrt(residual) / hnorm)
         end do
         do row = 1, nmat
            do column = 1, nmat
               left_value = (0.0_rp, 0.0_rp)
               do inner = 1, nmat
                  left_value = left_value + conjg(cuda_vectors(inner, row, ik)) * cuda_vectors(inner, column, ik)
               end do
               if (row == column) left_value = left_value - (1.0_rp, 0.0_rp)
               orthogonality_error = max(orthogonality_error, abs(left_value))
            end do
         end do
         group_start = 1
         do while (group_start <= nmat)
            group_end = group_start
            do while (group_end < nmat .and. abs(cpu_values(group_end + 1, ik) - cpu_values(group_end, ik)) <= group_tolerance)
               group_end = group_end + 1
            end do
            group_size = group_end - group_start + 1
            if (group_size > 1) then
               do row = 1, nmat
                  do column = 1, nmat
                     projector_cpu = (0.0_rp, 0.0_rp)
                     projector_cuda = (0.0_rp, 0.0_rp)
                     do inner = group_start, group_end
                        projector_cpu = projector_cpu + cpu_vectors(row, inner, ik) * conjg(cpu_vectors(column, inner, ik))
                        projector_cuda = projector_cuda + cuda_vectors(row, inner, ik) * conjg(cuda_vectors(column, inner, ik))
                     end do
                     projector_error = max(projector_error, abs(projector_cpu - projector_cuda))
                  end do
               end do
            end if
            group_start = group_end + 1
         end do
      end do
   end subroutine compare_eigenpairs

   integer function count_unique_points(recip, points) result(unique_count)
      type(reciprocal), intent(in) :: recip
      real(rp), intent(in) :: points(:, :)
      real(rp), allocatable :: folded(:, :)
      real(rp) :: point(3)
      integer :: ik, iunique
      logical :: seen

      allocate(folded(3, size(points, 2)))
      unique_count = 0
      do ik = 1, size(points, 2)
         call recip%fold_kpoint(points(:, ik), point)
         seen = .false.
         do iunique = 1, unique_count
            if (all(abs(point - folded(:, iunique)) <= 1.0e-12_rp)) then
               seen = .true.
               exit
            end if
         end do
         if (.not. seen) then
            unique_count = unique_count + 1
            folded(:, unique_count) = point
         end if
      end do
      deallocate(folded)
   end function count_unique_points

   subroutine run_preflight(input_file, fixture, backend_name, solver_strategy, mesh, tile_size, want_vectors, supercell_l)
      character(len=*), intent(in) :: input_file, fixture, backend_name, solver_strategy
      integer, intent(in) :: mesh(3), tile_size, supercell_l
      logical, intent(in) :: want_vectors
      type(reciprocal) :: recip
      type(hamiltonian), target :: ham
      type(lattice), target :: lat
      type(charge), target :: chg
      type(control), target :: ctl
      integer(c_long_long) :: free_bytes, total_bytes
      real(rp) :: free_mib, total_mib
      integer :: nmat
      integer :: strategy_status
      logical :: supported
      character(len=128) :: reason, support_word

      call setup_production(input_file, recip, ham, lat, chg, ctl)
      nmat = nb * lat%nrec
      call recip%make_execution_backend(backend_name)
      call configure_solver_strategy(recip, backend_name, solver_strategy, strategy_status)
      if (strategy_status /= 0) error stop 'ACCP1: failed to configure CUDA solver strategy'
      call backend_solver_supported(recip, nmat, tile_size, want_vectors, supported, reason)
      if (supported) then
         support_word = 'supported'
      else
         support_word = 'unsupported'
      end if
      call backend_memory(recip, free_bytes, total_bytes)
      free_mib = real(free_bytes, rp) / (1024.0_rp * 1024.0_rp)
      total_mib = real(total_bytes, rp) / (1024.0_rp * 1024.0_rp)
      write (*, '(a)') 'ACCP0_PREFLIGHT fixture='//trim(fixture)//' backend='//trim(backend_name)// &
         ' solver_strategy='//trim(solver_strategy)//' support_status='//trim(support_word)// &
         ' L='//trim(int_token(supercell_l))//' nmat='//trim(int_token(nmat))// &
         ' tile='//trim(int_token(tile_size))//' reason='//trim(underscored(reason))//' memory_estimate_mib='// &
         trim(real_token(estimate_memory_mib(nmat, tile_size, want_vectors, solver_strategy)))// &
         ' free_mib='//trim(real_token(free_mib))//' total_mib='//trim(real_token(total_mib))
   end subroutine run_preflight

   subroutine configure_solver_strategy(recip, backend_name, solver_strategy, status)
      type(reciprocal), intent(inout) :: recip
      character(len=*), intent(in) :: backend_name, solver_strategy
      integer, intent(out) :: status

      status = 0
      if (trim(backend_name) /= 'cuda') return
      if (.not. allocated(recip%execution_backend)) then
         status = 1
         return
      end if
      select type (backend => recip%execution_backend)
      type is (cuda_reciprocal_backend)
         status = backend%set_solver_strategy(solver_strategy)
      class default
         status = 1
      end select
   end subroutine configure_solver_strategy

   subroutine backend_solver_supported(recip, nmat, batch_size, want_vectors, supported, reason)
      type(reciprocal), intent(in) :: recip
      integer, intent(in) :: nmat, batch_size
      logical, intent(in) :: want_vectors
      logical, intent(out) :: supported
      character(len=*), intent(out) :: reason

      supported = .true.
      reason = 'supported'
      if (.not. allocated(recip%execution_backend)) return
      select type (backend => recip%execution_backend)
      type is (cuda_reciprocal_backend)
         call backend%solver_strategy_supported(nmat, batch_size, want_vectors, supported, reason)
      class default
         continue
      end select
   end subroutine backend_solver_supported

   subroutine reset_backend_metrics(recip)
      type(reciprocal), intent(inout) :: recip
      if (.not. allocated(recip%execution_backend)) return
      select type (backend => recip%execution_backend)
      type is (cuda_reciprocal_backend)
         call backend%reset_timing_metrics()
      class default
         continue
      end select
   end subroutine reset_backend_metrics

   subroutine backend_resource_counters(recip, cuda_malloc_count, cuda_free_count, workspace_query_count, workspace_reuse_count, &
                                        event_create_count, event_destroy_count, pinned_alloc_count, pinned_free_count)
      type(reciprocal), intent(in) :: recip
      integer(c_long_long), intent(out) :: cuda_malloc_count, cuda_free_count, workspace_query_count, workspace_reuse_count
      integer(c_long_long), intent(out) :: event_create_count, event_destroy_count, pinned_alloc_count, pinned_free_count

      cuda_malloc_count = 0_c_long_long; cuda_free_count = 0_c_long_long
      workspace_query_count = 0_c_long_long; workspace_reuse_count = 0_c_long_long
      event_create_count = 0_c_long_long; event_destroy_count = 0_c_long_long
      pinned_alloc_count = 0_c_long_long; pinned_free_count = 0_c_long_long
      if (.not. allocated(recip%execution_backend)) return
      select type (backend => recip%execution_backend)
      type is (cuda_reciprocal_backend)
         cuda_malloc_count = backend%cuda_malloc_count
         cuda_free_count = backend%cuda_free_count
         workspace_query_count = backend%workspace_query_count
         workspace_reuse_count = backend%workspace_reuse_count
         event_create_count = backend%event_create_count
         event_destroy_count = backend%event_destroy_count
         pinned_alloc_count = backend%pinned_alloc_count
         pinned_free_count = backend%pinned_free_count
      class default
         continue
      end select
   end subroutine backend_resource_counters

   subroutine backend_timings(recip, conversion, staging, h2d, solver, d2h, d2h_values, d2h_vectors, sync, widening, &
                              total_reciprocal, h2d_bytes, d2h_values_bytes, d2h_vectors_bytes, pinned_host_active, &
                              cuda_malloc_count, cuda_free_count, workspace_query_count, workspace_reuse_count, &
                              event_create_count, event_destroy_count, pinned_alloc_count, pinned_free_count)
      type(reciprocal), intent(in) :: recip
      real(rp), intent(out) :: conversion, staging, h2d, solver, d2h, d2h_values, d2h_vectors, sync, widening, total_reciprocal
      integer(c_long_long), intent(out) :: h2d_bytes, d2h_values_bytes, d2h_vectors_bytes
      integer, intent(out) :: pinned_host_active
      integer(c_long_long), intent(out) :: cuda_malloc_count, cuda_free_count, workspace_query_count, workspace_reuse_count
      integer(c_long_long), intent(out) :: event_create_count, event_destroy_count, pinned_alloc_count, pinned_free_count
      integer :: execute_requests, combined_requests, assemble_only, input_solves

      conversion = 0.0_rp
      staging = 0.0_rp
      h2d = 0.0_rp
      solver = 0.0_rp
      d2h = 0.0_rp
      d2h_values = 0.0_rp
      d2h_vectors = 0.0_rp
      sync = 0.0_rp
      widening = 0.0_rp
      total_reciprocal = 0.0_rp
      h2d_bytes = 0_c_long_long
      d2h_values_bytes = 0_c_long_long
      d2h_vectors_bytes = 0_c_long_long
      pinned_host_active = 0
      cuda_malloc_count = 0_c_long_long; cuda_free_count = 0_c_long_long
      workspace_query_count = 0_c_long_long; workspace_reuse_count = 0_c_long_long
      event_create_count = 0_c_long_long; event_destroy_count = 0_c_long_long
      pinned_alloc_count = 0_c_long_long; pinned_free_count = 0_c_long_long
      execute_requests = 0; combined_requests = 0; assemble_only = 0; input_solves = 0
      if (.not. allocated(recip%execution_backend)) return
      call recip%execution_backend%execution_metrics(execute_requests, combined_requests, assemble_only, input_solves)
      select type (backend => recip%execution_backend)
      type is (cuda_reciprocal_backend)
         conversion = backend%host_conversion_seconds
         staging = backend%host_staging_seconds
         h2d = backend%h2d_seconds
         solver = backend%gpu_solve_seconds
         d2h = backend%d2h_seconds
         d2h_values = backend%d2h_values_seconds
         d2h_vectors = backend%d2h_vectors_seconds
         sync = backend%sync_seconds
         widening = backend%host_widen_seconds
         total_reciprocal = backend%total_reciprocal_seconds
         h2d_bytes = backend%h2d_bytes
         d2h_values_bytes = backend%d2h_values_bytes
         d2h_vectors_bytes = backend%d2h_vectors_bytes
         pinned_host_active = backend%pinned_host_active
         cuda_malloc_count = backend%cuda_malloc_count
         cuda_free_count = backend%cuda_free_count
         workspace_query_count = backend%workspace_query_count
         workspace_reuse_count = backend%workspace_reuse_count
         event_create_count = backend%event_create_count
         event_destroy_count = backend%event_destroy_count
         pinned_alloc_count = backend%pinned_alloc_count
         pinned_free_count = backend%pinned_free_count
      class default
         continue
      end select
   end subroutine backend_timings

   subroutine backend_memory(recip, free_bytes, total_bytes)
      type(reciprocal), intent(in) :: recip
      integer(c_long_long), intent(out) :: free_bytes, total_bytes

      free_bytes = 0_c_long_long
      total_bytes = 0_c_long_long
      if (.not. allocated(recip%execution_backend)) return
      select type (backend => recip%execution_backend)
      type is (cuda_reciprocal_backend)
         call backend%memory_info(free_bytes, total_bytes)
      class default
         continue
      end select
   end subroutine backend_memory

   function median(values) result(value)
      real(rp), intent(in) :: values(:)
      real(rp) :: value
      real(rp), allocatable :: sorted(:)
      integer :: n

      n = size(values)
      allocate(sorted(n), source=values)
      call sort_values(sorted)
      if (mod(n, 2) == 1) then
         value = sorted((n + 1) / 2)
      else
         value = 0.5_rp * (sorted(n / 2) + sorted(n / 2 + 1))
      end if
      deallocate(sorted)
   end function median

   function estimate_memory_mib(nmat, tile_size, want_vectors, solver_strategy) result(memory_mib)
      integer, intent(in) :: nmat, tile_size
      logical, intent(in) :: want_vectors
      character(len=*), intent(in), optional :: solver_strategy
      real(rp) :: memory_mib
      real(rp) :: complex_bytes, real_bytes

      complex_bytes = 16.0_rp
      real_bytes = 8.0_rp
      if (present(solver_strategy)) then
         if (index(trim(solver_strategy), 'fp32_') == 1) then
            complex_bytes = 8.0_rp
            real_bytes = 4.0_rp
         end if
      end if
      memory_mib = real(tile_size, rp) * real(nmat * nmat, rp) * complex_bytes
      if (want_vectors) memory_mib = memory_mib + real(tile_size, rp) * real(nmat * nmat, rp) * complex_bytes
      memory_mib = memory_mib + real(tile_size * nmat, rp) * real_bytes
      memory_mib = memory_mib / (1024.0_rp * 1024.0_rp)
   end function estimate_memory_mib

   subroutine sort_values(values)
      real(rp), intent(inout) :: values(:)
      integer :: i, j
      real(rp) :: key
      do i = 2, size(values)
         key = values(i)
         j = i - 1
         do while (j >= 1 .and. values(j) > key)
            values(j + 1) = values(j)
            j = j - 1
         end do
         values(j + 1) = key
      end do
   end subroutine sort_values

   subroutine wall_clock(seconds)
      real(rp), intent(out) :: seconds
      integer(int64) :: count, rate
      call system_clock(count, rate)
      seconds = real(count, rp) / real(rate, rp)
   end subroutine wall_clock

   function int_token(value) result(token)
      integer, intent(in) :: value
      character(len=64) :: token
      write(token, '(i0)') value
   end function int_token

   function int64_token(value) result(token)
      integer(c_long_long), intent(in) :: value
      character(len=64) :: token
      write(token, '(i0)') value
   end function int64_token

   function underscored(value) result(token)
      character(len=*), intent(in) :: value
      character(len=128) :: token
      integer :: i, n

      token = ' '
      n = min(len(token), len_trim(value))
      if (n > 0) token(:n) = value(:n)
      do i = 1, n
         if (token(i:i) == ' ') token(i:i) = '_'
      end do
   end function underscored

   function real_token(value) result(token)
      real(rp), intent(in) :: value
      character(len=64) :: token
      write(token, '(es24.16e3)') value
      token = adjustl(token)
   end function real_token

   function mesh_token(mesh) result(token)
      integer, intent(in) :: mesh(3)
      character(len=64) :: token
      write(token, '(i0,a,i0,a,i0)') mesh(1), 'x', mesh(2), 'x', mesh(3)
   end function mesh_token

   subroutine append_token(line, token)
      character(len=*), intent(inout) :: line
      character(len=*), intent(in) :: token
      line = trim(line)//' '//trim(token)
   end subroutine append_token

   subroutine append_real(line, key, value)
      character(len=*), intent(inout) :: line
      character(len=*), intent(in) :: key
      real(rp), intent(in) :: value
      call append_token(line, trim(key)//'='//trim(real_token(value)))
   end subroutine append_real

end program accp0_real_material
