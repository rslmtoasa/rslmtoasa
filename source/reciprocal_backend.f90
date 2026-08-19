!------------------------------------------------------------------------------
! Reciprocal execution backend implementations.
! The module parent is private by default; only the typed contract exported by
! reciprocal_mod is visible to reciprocal SCF and response orchestration.
!------------------------------------------------------------------------------
submodule (reciprocal_mod) reciprocal_backend
   use lehmann_kernel_mod, only: lehmann_pair_block
   implicit none

contains

   module subroutine lapack_backend_contract_lehmann(this, request, result, status)
      class(lapack_reciprocal_backend), intent(inout) :: this
      type(reciprocal_lehmann_request), intent(in) :: request
      type(reciprocal_lehmann_result), intent(inout) :: result
      integer, intent(out) :: status
      integer :: nmat, nk, ne, npair, nblk, ipair, clock_start, clock_stop, clock_rate

      status = 0
      result%valid = .false.
      result%h2d_seconds = 0.0_rp
      result%contraction_seconds = 0.0_rp
      result%d2h_seconds = 0.0_rp
      if (allocated(result%blocks)) deallocate(result%blocks)

      if (.not. allocated(request%eigenvalues) .or. .not. allocated(request%eigenvectors) .or. &
          .not. allocated(request%k_points) .or. .not. allocated(request%z_contour) .or. &
          .not. allocated(request%dr) .or. .not. allocated(request%ioffset) .or. .not. allocated(request%joffset)) then
         status = 1
         return
      end if
      nmat = size(request%eigenvalues, 1)
      nk = size(request%eigenvalues, 2)
      ne = size(request%z_contour)
      npair = size(request%ioffset)
      nblk = request%nblk
      if (nmat < 1 .or. nk < 1 .or. ne < 1 .or. npair < 1 .or. nblk < 1 .or. &
          size(request%eigenvectors, 1) /= nmat .or. size(request%eigenvectors, 2) /= nmat .or. &
          size(request%eigenvectors, 3) /= nk .or. size(request%k_points, 1) /= 3 .or. &
          size(request%k_points, 2) /= nk .or. size(request%dr, 1) /= 3 .or. size(request%dr, 2) /= npair .or. &
          size(request%joffset) /= npair .or. nblk > nmat) then
         status = 1
         return
      end if

      allocate(result%blocks(nblk, nblk, ne, npair))
      call system_clock(count_rate=clock_rate)
      call system_clock(clock_start)
      do ipair = 1, npair
         call lehmann_pair_block(request%eigenvalues, request%eigenvectors, request%k_points, request%z_contour, &
                                 request%dr(:, ipair), request%ioffset(ipair), request%joffset(ipair), nblk, &
                                 result%blocks(:, :, :, ipair))
      end do
      call system_clock(clock_stop)
      result%contraction_seconds = real(clock_stop - clock_start, rp) / real(max(1, clock_rate), rp)
      result%valid = .true.
   end subroutine lapack_backend_contract_lehmann

   module subroutine cuda_backend_initialize(this, assembler)
      class(cuda_reciprocal_backend), intent(inout) :: this
      type(reciprocal_assembler), intent(in) :: assembler
#ifdef USE_CUDA_RECIPROCAL
      integer(c_int) :: device_count, status
      integer :: selected_device, override_device
      logical :: device_valid, override_configured, override_valid
      character(len=256) :: mapping_message
#endif

      call this%release()
      this%prepared_operator_generation = -1
      this%execute_batch_requests = 0
      this%input_hamiltonian_solve_requests = 0
      this%operator_prepare_requests = 0
      this%operator_prepare_reuses = 0
      this%h2d_seconds = 0.0_rp
      this%gpu_solve_seconds = 0.0_rp
      this%d2h_seconds = 0.0_rp
      this%timing_calls = 0

#ifdef USE_CUDA_RECIPROCAL
      device_count = 0_c_int
      status = rslmto_reciprocal_cuda_device_count(device_count)
      if (status /= 0_c_int .or. device_count <= 0_c_int) then
         call g_logger%error('reciprocal CUDA backend: no usable CUDA device was found.', __FILE__, __LINE__)
         return
      end if
      call get_cuda_device_override(override_device, override_configured, override_valid)
      if (override_configured .and. .not. override_valid) then
         call g_logger%error('reciprocal CUDA backend: RSLMTO_CUDA_DEVICE must be a non-negative integer.', __FILE__, __LINE__)
         return
      end if
      if (override_configured) then
         call g_parallel_context%device_index(int(device_count), selected_device, device_valid, override_device)
      else
         call g_parallel_context%device_index(int(device_count), selected_device, device_valid)
      end if
      if (.not. device_valid) then
         call g_logger%error('reciprocal CUDA backend: MPI local-rank/device mapping is invalid; '// &
                             'set RSLMTO_CUDA_DEVICE only for intentional sharing.', __FILE__, __LINE__)
         return
      end if
      write(mapping_message, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'CUDA_DEVICE_MAPPING world_rank=', &
         g_parallel_context%rank, ' local_rank=', g_parallel_context%local_rank, &
         ' local_size=', g_parallel_context%local_size, ' visible_devices=', int(device_count), &
         ' selected_device=', selected_device
      call g_logger%info(trim(mapping_message), __FILE__, __LINE__)
      this%context = rslmto_reciprocal_cuda_create(int(selected_device, c_int))
      if (.not. c_associated(this%context)) then
         call g_logger%error('reciprocal CUDA backend: context creation failed.', __FILE__, __LINE__)
         return
      end if
      this%device = selected_device
      this%initialized = .true.
#else
      call g_logger%error('reciprocal CUDA backend: executable was built without CUDA support.', __FILE__, __LINE__)
#endif
   end subroutine cuda_backend_initialize

   module subroutine cuda_backend_capabilities(this, capabilities)
      class(cuda_reciprocal_backend), intent(in) :: this
      type(reciprocal_execution_capabilities), intent(out) :: capabilities

      ! Do not advertise a usable numerical route after context creation
      ! failed (for example when a CUDA-enabled executable runs on a CPU-only
      ! node).  This lets callers reject the request at the typed seam instead
      ! of submitting an incomplete result and accidentally looking like a
      ! successful CPU fallback.
      ! ACC-09: first-order, second-order/HOH, and CCOR are all ordinary
      ! standard-Hermitian operators after the established CPU Fourier
      ! assembly.  CUDA owns only the dense eigensolution; leaving both
      ! assembly capabilities false deliberately selects that host handoff,
      ! rather than declaring those scientifically validated variants
      ! unsupported or duplicating their construction on the device.
      capabilities%standard_hermitian = this%initialized
      capabilities%generalized_hermitian = .false.
      capabilities%eigenvalues_only = this%initialized
      capabilities%eigenvectors = this%initialized
      capabilities%first_order_assembly = .false.
      capabilities%second_order_assembly = .false.
      capabilities%overlap = .false.
      capabilities%preferred_tile_size = 1
      capabilities%maximum_tile_size = huge(1)
      capabilities%residency = reciprocal_residency_device
   end subroutine cuda_backend_capabilities

   module subroutine cuda_backend_prepare_operator(this, operator_generation)
      class(cuda_reciprocal_backend), intent(inout) :: this
      integer, intent(in) :: operator_generation
#ifdef USE_CUDA_RECIPROCAL
      integer(c_int) :: status
#endif

      if (.not. this%initialized) then
         call g_logger%error('reciprocal CUDA backend: prepare requested before initialization.', __FILE__, __LINE__)
         return
      end if

#ifdef USE_CUDA_RECIPROCAL
      status = rslmto_reciprocal_cuda_prepare_operator(this%context, int(operator_generation, c_int))
      if (status < 0_c_int) then
         call g_logger%error('reciprocal CUDA backend: operator preparation failed.', __FILE__, __LINE__)
         return
      end if
      if (status == 0_c_int) then
         this%operator_prepare_reuses = this%operator_prepare_reuses + 1
      else
         this%operator_prepare_requests = this%operator_prepare_requests + 1
      end if
      this%prepared_operator_generation = operator_generation
#else
      call g_logger%error('reciprocal CUDA backend: operator preparation is unavailable in this build.', __FILE__, __LINE__)
#endif
   end subroutine cuda_backend_prepare_operator

   module subroutine cuda_backend_execute_batch(this, request, result)
      class(cuda_reciprocal_backend), intent(inout) :: this
      type(reciprocal_execution_request), intent(in) :: request
      type(reciprocal_execution_result), intent(inout) :: result
#ifdef USE_CUDA_RECIPROCAL
      integer(c_int) :: status
      integer :: nmat, nk
      complex(c_double_complex), target :: eigenvectors_dummy(1)
#endif

      call clear_execution_result(result)
      this%execute_batch_requests = this%execute_batch_requests + 1
#ifndef USE_CUDA_RECIPROCAL
      call g_logger%error('reciprocal CUDA backend: executable was built without CUDA support.', __FILE__, __LINE__)
#else
      if (.not. this%initialized) then
         call g_logger%error('reciprocal CUDA backend: execute requested before initialization.', __FILE__, __LINE__)
         return
      end if
      if (request%generalized .or. request%assemble_hamiltonian .or. request%assemble_overlap .or. &
          request%request_assembled_hamiltonian .or. request%request_assembled_overlap .or. &
          allocated(request%input_overlap) .or. .not. request%solve_eigensystem .or. &
          .not. allocated(request%input_hamiltonian)) then
         call g_logger%error('reciprocal CUDA backend: only standard-Hermitian eigensolution of host-assembled H(k) tiles is supported.', &
                             __FILE__, __LINE__)
         return
      end if
      nmat = size(request%input_hamiltonian, 1)
      nk = size(request%input_hamiltonian, 3)
      if (nmat < 1 .or. size(request%input_hamiltonian, 2) /= nmat .or. nk < 1) then
         call g_logger%error('reciprocal CUDA backend: input H(k) tile must be non-empty and square.', __FILE__, __LINE__)
         return
      end if
      call this%prepare_operator(request%operator_generation)
      this%input_hamiltonian_solve_requests = this%input_hamiltonian_solve_requests + 1
      allocate(result%eigenvalues(nmat, nk))
      if (request%request_eigenvectors) allocate(result%eigenvectors(nmat, nmat, nk))
      if (request%request_eigenvectors) then
         status = rslmto_reciprocal_cuda_solve_zheevd_batch(this%context, int(nmat, c_int), int(nk, c_int), &
            request%input_hamiltonian, result%eigenvalues, result%eigenvectors, 1_c_int)
      else
         status = rslmto_reciprocal_cuda_solve_zheevd_batch(this%context, int(nmat, c_int), int(nk, c_int), &
            request%input_hamiltonian, result%eigenvalues, eigenvectors_dummy, 0_c_int)
      end if
      if (status /= 0_c_int) then
         if (request%request_eigenvectors) then
            call g_logger%error('reciprocal CUDA backend: cuSOLVER eigensolution failed.', __FILE__, __LINE__)
         else
            call g_logger%error('reciprocal CUDA backend: cuSOLVER eigenvalue solve failed.', __FILE__, __LINE__)
         end if
         call clear_execution_result(result)
         return
      end if
      call rslmto_reciprocal_cuda_get_timings(this%context, this%h2d_seconds, this%gpu_solve_seconds, &
                                               this%d2h_seconds, this%timing_calls)
      result%local_point_count = nk
      result%operator_generation = request%operator_generation
      result%eigenvalues_valid = .true.
      result%eigenvectors_valid = request%request_eigenvectors
#endif
   end subroutine cuda_backend_execute_batch

   module subroutine cuda_backend_contract_lehmann(this, request, result, status)
      class(cuda_reciprocal_backend), intent(inout) :: this
      type(reciprocal_lehmann_request), intent(in) :: request
      type(reciprocal_lehmann_result), intent(inout) :: result
      integer, intent(out) :: status
#ifdef USE_CUDA_RECIPROCAL
      integer(c_int) :: c_status
      integer :: nmat, nk, ne, npair, nblk
#endif

      status = 1
      result%valid = .false.
      result%h2d_seconds = 0.0_rp
      result%contraction_seconds = 0.0_rp
      result%d2h_seconds = 0.0_rp
      if (allocated(result%blocks)) deallocate(result%blocks)
#ifndef USE_CUDA_RECIPROCAL
      call g_logger%error('reciprocal CUDA backend: Lehmann contraction is unavailable in this build.', __FILE__, __LINE__)
#else
      if (.not. this%initialized .or. .not. c_associated(this%context)) then
         call g_logger%error('reciprocal CUDA backend: Lehmann contraction requested before initialization.', __FILE__, __LINE__)
         return
      end if
      if (.not. allocated(request%eigenvalues) .or. .not. allocated(request%eigenvectors) .or. &
          .not. allocated(request%k_points) .or. .not. allocated(request%z_contour) .or. &
          .not. allocated(request%dr) .or. .not. allocated(request%ioffset) .or. .not. allocated(request%joffset)) then
         call g_logger%error('reciprocal CUDA backend: incomplete Lehmann request.', __FILE__, __LINE__)
         return
      end if
      nmat = size(request%eigenvalues, 1)
      nk = size(request%eigenvalues, 2)
      ne = size(request%z_contour)
      npair = size(request%ioffset)
      nblk = request%nblk
      if (nmat < 1 .or. nk < 1 .or. ne < 1 .or. npair < 1 .or. nblk < 1 .or. &
          size(request%eigenvectors, 1) /= nmat .or. size(request%eigenvectors, 2) /= nmat .or. &
          size(request%eigenvectors, 3) /= nk .or. size(request%k_points, 1) /= 3 .or. &
          size(request%k_points, 2) /= nk .or. size(request%dr, 1) /= 3 .or. size(request%dr, 2) /= npair .or. &
          size(request%joffset) /= npair .or. nblk > nmat) then
         call g_logger%error('reciprocal CUDA backend: invalid Lehmann request dimensions.', __FILE__, __LINE__)
         return
      end if

      allocate(result%blocks(nblk, nblk, ne, npair))
      c_status = rslmto_reciprocal_cuda_contract_lehmann(this%context, int(nmat, c_int), int(nk, c_int), &
         int(ne, c_int), int(npair, c_int), int(nblk, c_int), request%eigenvalues, request%eigenvectors, &
         request%k_points, request%z_contour, request%dr, request%ioffset, request%joffset, result%blocks, &
         result%h2d_seconds, result%contraction_seconds, result%d2h_seconds)
      status = int(c_status)
      if (status /= 0) then
         call g_logger%error('reciprocal CUDA backend: Lehmann contraction failed.', __FILE__, __LINE__)
         deallocate(result%blocks)
         return
      end if
      result%valid = .true.
#endif
   end subroutine cuda_backend_contract_lehmann

   module subroutine cuda_backend_execution_metrics(this, execute_requests, combined_requests, assemble_only, input_hamiltonian_solves)
      class(cuda_reciprocal_backend), intent(in) :: this
      integer, intent(out) :: execute_requests, combined_requests, assemble_only, input_hamiltonian_solves

      execute_requests = this%execute_batch_requests
      combined_requests = 0
      assemble_only = 0
      input_hamiltonian_solves = this%input_hamiltonian_solve_requests
   end subroutine cuda_backend_execution_metrics

   module subroutine cuda_backend_synchronize(this)
      class(cuda_reciprocal_backend), intent(inout) :: this
#ifdef USE_CUDA_RECIPROCAL
      integer(c_int) :: status
#endif

      if (.not. this%initialized) return
#ifdef USE_CUDA_RECIPROCAL
      status = rslmto_reciprocal_cuda_synchronize(this%context)
      if (status /= 0_c_int) then
         call g_logger%error('reciprocal CUDA backend: stream synchronization failed.', __FILE__, __LINE__)
      end if
#endif
   end subroutine cuda_backend_synchronize

   module subroutine cuda_backend_release(this)
      class(cuda_reciprocal_backend), intent(inout) :: this
#ifdef USE_CUDA_RECIPROCAL
      if (c_associated(this%context)) call rslmto_reciprocal_cuda_destroy(this%context)
#endif
      this%context = c_null_ptr
      this%device = -1
      this%prepared_operator_generation = -1
      this%h2d_seconds = 0.0_rp
      this%gpu_solve_seconds = 0.0_rp
      this%d2h_seconds = 0.0_rp
      this%timing_calls = 0
      this%initialized = .false.
   end subroutine cuda_backend_release

   module subroutine cuda_backend_destructor(this)
      type(cuda_reciprocal_backend), intent(inout) :: this
      call this%release()
   end subroutine cuda_backend_destructor

   module subroutine make_execution_backend(this, backend_name)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in), optional :: backend_name
      type(reciprocal_assembler) :: assembler
      character(len=16) :: selected_backend

      selected_backend = trim(this%reciprocal_backend)
      if (selected_backend == '') selected_backend = 'lapack'
      if (present(backend_name)) then
         selected_backend = trim(lower(backend_name))
         select case (selected_backend)
         case ('', 'lapack', 'cpu')
            selected_backend = 'lapack'
         case ('cuda')
         case default
            call g_logger%fatal('reciprocal backend factory: requested backend is unavailable.', __FILE__, __LINE__)
         end select
         ! An explicit programmatic selection is persistent for the following
         ! normal-mesh calls, which intentionally omit the optional argument.
         this%reciprocal_backend = selected_backend
      end if
      call this%make_reciprocal_assembler(assembler)
      if (allocated(this%execution_backend)) then
         select type (backend => this%execution_backend)
         type is (lapack_reciprocal_backend)
            if (selected_backend == 'cuda') then
               call backend%release()
               deallocate(this%execution_backend)
            else
               ! Geometry/order choices can change across SCF probes.  Refresh
               ! the non-owning assembler view without discarding its workspace.
               backend%assembler = assembler
               backend%prepared_operator_generation = -1
               return
            end if
         type is (cuda_reciprocal_backend)
            if (selected_backend == 'cuda') return
            call backend%release()
            deallocate(this%execution_backend)
         class default
            call g_logger%fatal('reciprocal backend factory: unknown configured backend.', __FILE__, __LINE__)
         end select
      end if

      select case (selected_backend)
      case ('lapack')
         allocate(lapack_reciprocal_backend :: this%execution_backend)
         call this%execution_backend%initialize(assembler)
      case ('cuda')
#ifdef USE_CUDA_RECIPROCAL
         allocate(cuda_reciprocal_backend :: this%execution_backend)
         call this%execution_backend%initialize(assembler)
#else
         call g_logger%fatal('reciprocal backend factory: CUDA backend requested but ENABLE_CUDA_RECIPROCAL=OFF.', &
                             __FILE__, __LINE__)
#endif
      end select
   end subroutine make_execution_backend

   module subroutine lapack_backend_initialize(this, assembler)
      class(lapack_reciprocal_backend), intent(inout) :: this
      type(reciprocal_assembler), intent(in) :: assembler

      call this%release()
      this%assembler = assembler
      this%execute_batch_requests = 0
      this%combined_assembly_solve_requests = 0
      this%assemble_only_requests = 0
      this%input_hamiltonian_solve_requests = 0
      this%host_assembly_seconds = 0.0_rp
      this%initialized = .true.
   end subroutine lapack_backend_initialize

   module subroutine lapack_backend_capabilities(this, capabilities)
      class(lapack_reciprocal_backend), intent(in) :: this
      type(reciprocal_execution_capabilities), intent(out) :: capabilities

      capabilities%standard_hermitian = .true.
      capabilities%generalized_hermitian = .true.
      capabilities%eigenvalues_only = .true.
      capabilities%eigenvectors = .true.
      capabilities%first_order_assembly = .true.
      capabilities%second_order_assembly = .true.
      capabilities%overlap = .true.
      capabilities%preferred_tile_size = max(1, this%workspace%tile_capacity)
      capabilities%maximum_tile_size = huge(1)
      capabilities%residency = reciprocal_residency_host
   end subroutine lapack_backend_capabilities

   module subroutine lapack_backend_prepare_operator(this, operator_generation)
      class(lapack_reciprocal_backend), intent(inout) :: this
      integer, intent(in) :: operator_generation

      if (.not. this%initialized) then
         call g_logger%fatal('reciprocal backend: prepare requested before initialization.', __FILE__, __LINE__)
      end if
      this%prepared_operator_generation = operator_generation
   end subroutine lapack_backend_prepare_operator

   module subroutine lapack_backend_execute_batch(this, request, result)
      use, intrinsic :: ieee_exceptions, only: ieee_divide_by_zero, &
                                               ieee_get_halting_mode, ieee_set_flag, ieee_set_halting_mode
      class(lapack_reciprocal_backend), intent(inout) :: this
      type(reciprocal_execution_request), intent(in) :: request
      type(reciprocal_execution_result), intent(inout) :: result

      integer :: nmat, nk, ik, nnmax, ntype, assembly_start, assembly_stop, assembly_rate
      logical :: halt_divide_by_zero, use_assembled_input
      character(len=1) :: jobz

      if (.not. this%initialized) then
         call g_logger%fatal('reciprocal backend: execute requested before initialization.', __FILE__, __LINE__)
      end if
      use_assembled_input = allocated(request%input_hamiltonian)
      if (request%assemble_hamiltonian .eqv. use_assembled_input) then
         call g_logger%fatal('reciprocal backend: request exactly one H(k) source (points or input matrices).', __FILE__, __LINE__)
      end if
      if (request%generalized .and. .not. request%assemble_overlap .and. .not. allocated(request%input_overlap)) then
         call g_logger%fatal('reciprocal backend: generalized solve requires an overlap source.', __FILE__, __LINE__)
      end if
      if (request%assemble_hamiltonian) then
         if (.not. allocated(request%k_points) .or. size(request%k_points,1) /= 3) then
            call g_logger%fatal('reciprocal backend: assembly request requires k_points(3,nk).', __FILE__, __LINE__)
         end if
         nk = size(request%k_points,2)
         if (.not. associated(this%assembler%lattice)) then
            call g_logger%fatal('reciprocal backend: assembler lattice state is unavailable.', __FILE__, __LINE__)
         end if
         nmat = nb*this%assembler%lattice%nrec
         nnmax = this%assembler%lattice%nn_max
         ntype = this%assembler%lattice%ntype
      else
         nmat = size(request%input_hamiltonian,1)
         nk = size(request%input_hamiltonian,3)
         if (size(request%input_hamiltonian,2) /= nmat) then
            call g_logger%fatal('reciprocal backend: input H(k) must be square.', __FILE__, __LINE__)
         end if
         nnmax = max(1, this%workspace%tile_capacity)
         ntype = 1
      end if
      if (nk < 1 .or. nmat < 1) then
         call g_logger%fatal('reciprocal backend: empty execution batch.', __FILE__, __LINE__)
      end if
      if (allocated(request%input_overlap)) then
         if (size(request%input_overlap,1) /= nmat .or. size(request%input_overlap,2) /= nmat .or. &
             size(request%input_overlap,3) /= nk) then
            call g_logger%fatal('reciprocal backend: input S(k) shape does not match H(k).', __FILE__, __LINE__)
         end if
      end if

      this%execute_batch_requests = this%execute_batch_requests + 1
      if (request%assemble_hamiltonian .and. request%solve_eigensystem) then
         this%combined_assembly_solve_requests = this%combined_assembly_solve_requests + 1
      else if (request%assemble_hamiltonian) then
         this%assemble_only_requests = this%assemble_only_requests + 1
      else if (request%solve_eigensystem) then
         this%input_hamiltonian_solve_requests = this%input_hamiltonian_solve_requests + 1
      end if

      call this%prepare_operator(request%operator_generation)
      call this%workspace%ensure_capacity(nmat, nk, request%generalized, request%operator_generation, nnmax, ntype)
      if (request%assemble_hamiltonian) then
         call system_clock(count_rate=assembly_rate)
         call system_clock(assembly_start)
         call this%assembler%assemble_batch(request%k_points, this%workspace)
         call system_clock(assembly_stop)
         this%host_assembly_seconds = this%host_assembly_seconds + &
                                      real(assembly_stop - assembly_start, rp) / real(assembly_rate, rp)
      else
         this%workspace%h(:,:,1:nk) = request%input_hamiltonian(:,:,1:nk)
         this%workspace%active_tile_length = nk
      end if
      if (request%assemble_overlap) then
         call this%assembler%assemble_overlap_batch(request%k_points, this%workspace)
      else if (allocated(request%input_overlap)) then
         this%workspace%s(:,:,1:nk) = request%input_overlap(:,:,1:nk)
      end if

      call clear_execution_result(result)
      result%local_point_count = nk
      result%operator_generation = request%operator_generation
      if (request%request_assembled_hamiltonian) then
         allocate(result%hamiltonian(nmat,nmat,nk))
         result%hamiltonian = this%workspace%h(:,:,1:nk)
         result%assembled_hamiltonian_valid = .true.
      end if
      if (request%request_assembled_overlap) then
         if (.not. request%assemble_overlap .and. .not. allocated(request%input_overlap)) then
            call g_logger%fatal('reciprocal backend: assembled overlap was requested without an overlap source.', __FILE__, __LINE__)
         end if
         allocate(result%overlap(nmat,nmat,nk))
         result%overlap = this%workspace%s(:,:,1:nk)
         result%assembled_overlap_valid = .true.
      end if
      if (.not. request%solve_eigensystem) return

      allocate(result%eigenvalues(nmat,nk))
      if (request%request_eigenvectors) allocate(result%eigenvectors(nmat,nmat,nk))
      jobz = merge('V', 'N', request%request_eigenvectors)
      do ik = 1, nk
         call assert_hermitian(this%workspace%h(:,:,ik), 'H(k)')
         if (request%generalized) call assert_overlap(this%workspace%s(:,:,ik), this%workspace%overlap_cholesky)
         this%workspace%eigenvector(:,:,ik) = this%workspace%h(:,:,ik)
         ! oneMKL's Hermitian eigensolvers can evaluate an intermediate
         ! divide-by-zero in their scaling path.  Keep application FPE
         ! diagnostics enabled, masking that external-library exception only.
         call ieee_get_halting_mode(ieee_divide_by_zero, halt_divide_by_zero)
         call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
         if (request%generalized) then
            call zhegv(1, jobz, 'U', nmat, this%workspace%eigenvector(:,:,ik), nmat, this%workspace%s(:,:,ik), nmat, &
                       result%eigenvalues(:,ik), this%workspace%lapack_work, this%workspace%lwork, &
                       this%workspace%lapack_rwork, this%workspace%info(ik))
         else
            call zheev(jobz, 'U', nmat, this%workspace%eigenvector(:,:,ik), nmat, result%eigenvalues(:,ik), &
                       this%workspace%lapack_work, this%workspace%lwork, this%workspace%lapack_rwork, this%workspace%info(ik))
         end if
         call ieee_set_flag(ieee_divide_by_zero, .false.)
         call ieee_set_halting_mode(ieee_divide_by_zero, halt_divide_by_zero)
         if (this%workspace%info(ik) /= 0) then
            call g_logger%fatal('reciprocal backend: LAPACK eigensolver failed.', __FILE__, __LINE__)
         end if
         if (request%request_eigenvectors) result%eigenvectors(:,:,ik) = this%workspace%eigenvector(:,:,ik)
      end do
      result%eigenvalues_valid = .true.
      result%eigenvectors_valid = request%request_eigenvectors
   end subroutine lapack_backend_execute_batch

   module subroutine lapack_backend_execution_metrics(this, execute_requests, combined_requests, assemble_only, input_hamiltonian_solves)
      class(lapack_reciprocal_backend), intent(in) :: this
      integer, intent(out) :: execute_requests, combined_requests, assemble_only, input_hamiltonian_solves

      execute_requests = this%execute_batch_requests
      combined_requests = this%combined_assembly_solve_requests
      assemble_only = this%assemble_only_requests
      input_hamiltonian_solves = this%input_hamiltonian_solve_requests
   end subroutine lapack_backend_execution_metrics

   module subroutine lapack_backend_synchronize(this)
      class(lapack_reciprocal_backend), intent(inout) :: this
      ! Host LAPACK is synchronous.  Kept as an explicit residency boundary.
   end subroutine lapack_backend_synchronize

   module subroutine lapack_backend_release(this)
      class(lapack_reciprocal_backend), intent(inout) :: this
      if (allocated(this%assembler%ham_vec_type_direct)) deallocate(this%assembler%ham_vec_type_direct)
      nullify(this%assembler%hamiltonian, this%assembler%lattice, this%assembler%control)
      call this%workspace%clear()
      this%prepared_operator_generation = -1
      this%initialized = .false.
   end subroutine lapack_backend_release

   module subroutine lapack_backend_destructor(this)
      type(lapack_reciprocal_backend), intent(inout) :: this
      call this%release()
   end subroutine lapack_backend_destructor

   subroutine clear_execution_result(result)
      type(reciprocal_execution_result), intent(inout) :: result
      if (allocated(result%hamiltonian)) deallocate(result%hamiltonian)
      if (allocated(result%overlap)) deallocate(result%overlap)
      if (allocated(result%eigenvalues)) deallocate(result%eigenvalues)
      if (allocated(result%eigenvectors)) deallocate(result%eigenvectors)
      result%local_point_count = 0
      result%operator_generation = -1
      result%assembled_hamiltonian_valid = .false.
      result%assembled_overlap_valid = .false.
      result%eigenvalues_valid = .false.
      result%eigenvectors_valid = .false.
   end subroutine clear_execution_result

   subroutine assert_hermitian(matrix, label)
      complex(rp), intent(in) :: matrix(:, :)
      character(len=*), intent(in) :: label
      real(rp) :: max_herm, scale
      max_herm = maxval(abs(matrix-transpose(conjg(matrix))))
      scale = max(1.0_rp,maxval(abs(matrix)))
      if (max_herm > 1.0e-10_rp*scale) then
         call g_logger%fatal('reciprocal backend: '//trim(label)//' is non-Hermitian before eigensolution.', __FILE__, __LINE__)
      end if
   end subroutine assert_hermitian

   subroutine assert_overlap(overlap, chol)
      complex(rp), intent(in) :: overlap(:, :)
      complex(rp), intent(inout) :: chol(:, :)
      integer :: info, n
      call assert_hermitian(overlap, 'S(k)')
      n = size(overlap,1)
      if (size(overlap, 2) /= n .or. any(shape(chol) /= [n, n])) then
         call g_logger%fatal('reciprocal backend: overlap validation scratch shape is invalid.', __FILE__, __LINE__)
      end if
      ! Validate a copy: ZPOTRF is destructive and the original S(k) is the
      ! input metric that the following ZHEGV call must receive unchanged.
      chol = overlap
      call zpotrf('U', n, chol, n, info)
      if (info /= 0) then
         call g_logger%fatal('reciprocal backend: S(k) is not positive definite for generalized eigensolution.', __FILE__, __LINE__)
      end if
   end subroutine assert_overlap

end submodule reciprocal_backend
