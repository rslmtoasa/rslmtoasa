! ACC-P0 persistent real-material benchmark driver.
!
! This executable deliberately owns the benchmark lifecycle.  A single
! process reads a production input, builds the production lattice/charge/
! Hamiltonian/reciprocal stack, creates one execution backend, performs the
! requested warm-ups, resets interval counters, and then measures repetitions.
! It is opt-in and is not registered as a default CTest test.
program accp0_real_material
   use, intrinsic :: iso_c_binding, only: c_long_long
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

   character(len=64) :: mode, backend, fixture, workload
   character(len=128) :: source_name
   character(len=512) :: input_file, dump_file
   integer :: mesh(3), tile_size, eigenvectors, warmups, repetitions, supercell_l
   integer :: status

   mode = 'persistent'
   backend = 'lapack'
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

   call parse_arguments(mode, backend, fixture, source_name, workload, input_file, dump_file, mesh, tile_size, &
      eigenvectors, warmups, repetitions, supercell_l)
   if (any(mesh < 1) .or. tile_size < 1 .or. eigenvectors < 0 .or. eigenvectors > 1 .or. warmups < 0 .or. repetitions < 1) then
      error stop 'ACCP0: invalid benchmark dimensions or repetition policy'
   end if
   if (trim(mode) /= 'persistent' .and. trim(mode) /= 'oracle' .and. trim(mode) /= 'preflight') then
      error stop 'ACCP0: mode must be persistent, oracle, or preflight'
   end if
   if (trim(backend) /= 'lapack' .and. trim(backend) /= 'cuda') then
      error stop 'ACCP0: backend must be lapack or cuda'
   end if
   call g_parallel_context%restore_to_default()
   g_parallel_context = parallel_context()
   g_timer = timer()
   call init_math_operators()
   call g_logger%init()

   if (trim(mode) == 'oracle') then
      call run_oracle(trim(input_file), trim(dump_file), trim(fixture), trim(backend), mesh, supercell_l)
   else if (trim(mode) == 'preflight') then
      call run_preflight(trim(input_file), trim(fixture), trim(backend), mesh, tile_size, eigenvectors == 1, supercell_l)
   else
      call run_persistent(trim(input_file), trim(fixture), trim(source_name), trim(workload), trim(backend), mesh, &
         tile_size, eigenvectors == 1, warmups, repetitions, supercell_l)
   end if

contains

   subroutine parse_arguments(mode, backend, fixture, source_name, workload, input_file, dump_file, mesh, tile_size, &
                              eigenvectors, warmups, repetitions, supercell_l)
      character(len=*), intent(inout) :: mode, backend, fixture, source_name, workload, input_file, dump_file
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
      call recip%execution_backend%synchronize()
   end subroutine solve_tiles

   subroutine run_oracle(input_file, dump_file, fixture, backend_name, mesh, supercell_l)
      character(len=*), intent(in) :: input_file, dump_file, fixture, backend_name
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
      call recip%execution_backend%prepare_operator(ham%operator_generation)
      open (newunit=dump_unit, file=dump_file, status='replace', action='write')
      call solve_tiles(recip, hamiltonians, size(points, 2), .false., status, dump_unit)
      close(dump_unit)
      if (status /= 0) error stop 'ACCP0: oracle eigensolution failed'
      unique_nk = count_unique_points(recip, points)
      write (*, '(a)') 'ACCP0_ORACLE fixture='//trim(fixture)//' L='//trim(int_token(supercell_l))// &
         ' backend='//trim(backend_name)//' nmat='//trim(int_token(nb * lat%nrec))// &
         ' actual_unique_nk='//trim(int_token(unique_nk))
   end subroutine run_oracle

   subroutine run_persistent(input_file, fixture, source_name, workload, backend_name, mesh, tile_size, want_vectors, &
                             warmups, repetitions, supercell_l)
      character(len=*), intent(in) :: input_file, fixture, source_name, workload, backend_name
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
      real(rp) :: steady_min, steady_max, steady_spread, h2d, solver, d2h, post_seconds
      real(rp) :: memory_estimate, free_memory, total_memory
      integer :: status, i, nk, unique_nk, nmat
      integer(c_long_long) :: free_bytes, total_bytes
      character(len=2048) :: line

      call setup_production(input_file, recip, ham, lat, chg, ctl)
      call make_mesh(mesh, points)
      nk = size(points, 2)
      unique_nk = count_unique_points(recip, points)
      nmat = nb * lat%nrec
      call wall_clock(start)
      call recip%make_execution_backend(backend_name)
      call wall_clock(finish)
      backend_init = finish - start

      call wall_clock(start)
      call recip%execution_backend%prepare_operator(ham%operator_generation)
      call wall_clock(finish)
      prepare_seconds = finish - start

      call backend_memory(recip, free_bytes, total_bytes)
      free_memory = real(free_bytes, rp) / (1024.0_rp * 1024.0_rp)
      total_memory = real(total_bytes, rp) / (1024.0_rp * 1024.0_rp)
      memory_estimate = estimate_memory_mib(nmat, tile_size, want_vectors)
      if (trim(backend_name) == 'cuda' .and. total_memory > 0.0_rp .and. free_memory < 2.0_rp * memory_estimate) then
         error stop 'ACCP0: CUDA memory preflight rejected this workload'
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
      call backend_timings(recip, h2d, solver, d2h)
      ! CUDA counters are reset immediately before this measured interval and
      ! accumulate over all measured repetitions.  Normalize them to the
      ! same per-repetition unit as the reported steady medians; retain the
      ! interval sample count in the machine-readable record.
      h2d = h2d / real(repetitions, rp)
      solver = solver / real(repetitions, rp)
      d2h = d2h / real(repetitions, rp)
      post_seconds = max(0.0_rp, total_median - hk_cpu - solve_median)

      line = 'ACCP0_DIMENSIONS'
      call append_token(line, 'fixture='//trim(fixture))
      call append_token(line, 'source='//trim(source_name))
      call append_token(line, 'workload='//trim(workload))
      call append_token(line, 'backend='//trim(backend_name))
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
      call append_real(line, 'cold_process_wall_s', 0.0_rp)
      call append_real(line, 'cuda_context_backend_init_s', backend_init)
      call append_real(line, 'operator_prepare_s', prepare_seconds)
      call append_real(line, 'first_solve_s', first_solve)
      call append_real(line, 'steady_solve_median_s', solve_median)
      call append_real(line, 'steady_solve_min_s', minval(solve_times))
      call append_real(line, 'steady_solve_spread_s', maxval(solve_times) - minval(solve_times))
      call append_token(line, 'metric_repetitions='//trim(int_token(repetitions)))
      call append_real(line, 'Hk_CPU_s', hk_cpu)
      call append_real(line, 'H2D_s', h2d)
      call append_real(line, 'solver_s', solver)
      call append_real(line, 'D2H_s', d2h)
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

   subroutine run_preflight(input_file, fixture, backend_name, mesh, tile_size, want_vectors, supercell_l)
      character(len=*), intent(in) :: input_file, fixture, backend_name
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

      call setup_production(input_file, recip, ham, lat, chg, ctl)
      nmat = nb * lat%nrec
      call recip%make_execution_backend(backend_name)
      call backend_memory(recip, free_bytes, total_bytes)
      free_mib = real(free_bytes, rp) / (1024.0_rp * 1024.0_rp)
      total_mib = real(total_bytes, rp) / (1024.0_rp * 1024.0_rp)
      write (*, '(a)') 'ACCP0_PREFLIGHT fixture='//trim(fixture)//' backend='//trim(backend_name)// &
         ' L='//trim(int_token(supercell_l))//' nmat='//trim(int_token(nmat))// &
         ' tile='//trim(int_token(tile_size))//' memory_estimate_mib='//trim(real_token(estimate_memory_mib(nmat, tile_size, want_vectors)))// &
         ' free_mib='//trim(real_token(free_mib))//' total_mib='//trim(real_token(total_mib))
   end subroutine run_preflight

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

   subroutine backend_timings(recip, h2d, solver, d2h)
      type(reciprocal), intent(in) :: recip
      real(rp), intent(out) :: h2d, solver, d2h
      integer :: execute_requests, combined_requests, assemble_only, input_solves

      h2d = 0.0_rp
      solver = 0.0_rp
      d2h = 0.0_rp
      execute_requests = 0; combined_requests = 0; assemble_only = 0; input_solves = 0
      if (.not. allocated(recip%execution_backend)) return
      call recip%execution_backend%execution_metrics(execute_requests, combined_requests, assemble_only, input_solves)
      select type (backend => recip%execution_backend)
      type is (cuda_reciprocal_backend)
         h2d = backend%h2d_seconds
         solver = backend%gpu_solve_seconds
         d2h = backend%d2h_seconds
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

   function estimate_memory_mib(nmat, tile_size, want_vectors) result(memory_mib)
      integer, intent(in) :: nmat, tile_size
      logical, intent(in) :: want_vectors
      real(rp) :: memory_mib

      memory_mib = real(tile_size, rp) * real(nmat * nmat, rp) * 16.0_rp
      if (want_vectors) memory_mib = memory_mib + real(tile_size, rp) * real(nmat * nmat, rp) * 16.0_rp
      memory_mib = memory_mib + real(tile_size * nmat, rp) * 8.0_rp
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
