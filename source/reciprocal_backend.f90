!------------------------------------------------------------------------------
! Reciprocal execution backend implementations.
! The module parent is private by default; only the typed contract exported by
! reciprocal_mod is visible to reciprocal SCF and response orchestration.
!------------------------------------------------------------------------------
submodule (reciprocal_mod) reciprocal_backend
   implicit none

contains

   module subroutine make_execution_backend(this, backend_name)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in), optional :: backend_name
      type(reciprocal_assembler) :: assembler

      if (present(backend_name)) then
         select case (trim(lower(backend_name)))
         case ('', 'lapack', 'cpu')
         case default
            call g_logger%fatal('reciprocal backend factory: requested backend is unavailable.', __FILE__, __LINE__)
         end select
      end if
      call this%make_reciprocal_assembler(assembler)
      if (allocated(this%execution_backend)) then
         select type (backend => this%execution_backend)
         type is (lapack_reciprocal_backend)
            ! Geometry/order choices can change across SCF probes.  Refresh
            ! the non-owning assembler view without discarding its workspace.
            backend%assembler = assembler
            backend%prepared_operator_generation = -1
         class default
            call g_logger%fatal('reciprocal backend factory: unknown configured backend.', __FILE__, __LINE__)
         end select
         return
      end if
      allocate(lapack_reciprocal_backend :: this%execution_backend)
      call this%execution_backend%initialize(assembler)
   end subroutine make_execution_backend

   module subroutine lapack_backend_initialize(this, assembler)
      class(lapack_reciprocal_backend), intent(inout) :: this
      type(reciprocal_assembler), intent(in) :: assembler

      call this%release()
      this%assembler = assembler
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
      class(lapack_reciprocal_backend), intent(inout) :: this
      type(reciprocal_execution_request), intent(in) :: request
      type(reciprocal_execution_result), intent(inout) :: result

      integer :: nmat, nk, ik, nnmax, ntype
      logical :: use_assembled_input
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

      call this%prepare_operator(request%operator_generation)
      call this%workspace%ensure_capacity(nmat, nk, request%generalized, request%operator_generation, nnmax, ntype)
      if (request%assemble_hamiltonian) then
         call this%assembler%assemble_batch(request%k_points, this%workspace)
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
         if (request%generalized) call assert_overlap(this%workspace%s(:,:,ik))
         this%workspace%eigenvector(:,:,ik) = this%workspace%h(:,:,ik)
         if (request%generalized) then
            call zhegv(1, jobz, 'U', nmat, this%workspace%eigenvector(:,:,ik), nmat, this%workspace%s(:,:,ik), nmat, &
                       result%eigenvalues(:,ik), this%workspace%lapack_work, this%workspace%lwork, &
                       this%workspace%lapack_rwork, this%workspace%info(ik))
         else
            call zheev(jobz, 'U', nmat, this%workspace%eigenvector(:,:,ik), nmat, result%eigenvalues(:,ik), &
                       this%workspace%lapack_work, this%workspace%lwork, this%workspace%lapack_rwork, this%workspace%info(ik))
         end if
         if (this%workspace%info(ik) /= 0) then
            call g_logger%fatal('reciprocal backend: LAPACK eigensolver failed.', __FILE__, __LINE__)
         end if
         if (request%request_eigenvectors) result%eigenvectors(:,:,ik) = this%workspace%eigenvector(:,:,ik)
      end do
      result%eigenvalues_valid = .true.
      result%eigenvectors_valid = request%request_eigenvectors
   end subroutine lapack_backend_execute_batch

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

   subroutine assert_overlap(overlap)
      complex(rp), intent(in) :: overlap(:, :)
      complex(rp), allocatable :: chol(:, :)
      integer :: info, n
      call assert_hermitian(overlap, 'S(k)')
      n = size(overlap,1)
      allocate(chol(n,n), source=overlap)
      call zpotrf('U', n, chol, n, info)
      deallocate(chol)
      if (info /= 0) then
         call g_logger%fatal('reciprocal backend: S(k) is not positive definite for generalized eigensolution.', __FILE__, __LINE__)
      end if
   end subroutine assert_overlap

end submodule reciprocal_backend
