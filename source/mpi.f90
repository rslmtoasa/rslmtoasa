!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief MPI compatibility state and the process-local parallel context.
!> @details `rank` and `numprocs` remain compatibility views for legacy code.
!> `parallel_context` is their sole runtime source and is initialized by main
!> after MPI_Init (or deterministically serial when MPI is unavailable).
module mpi_mod
   use logger_mod, only: g_logger
#ifdef USE_MPI
   use mpi
#endif
   implicit none
   private

   integer, public :: ierr = 0
   integer, public :: rank = 0
   integer, public :: numprocs = 1
   integer, public :: remainder_atoms = 0
   integer, public :: atoms_per_process = 0
   integer, public :: start_atom = 1
   integer, public :: end_atom = 0
   logical, public :: oversubscribe_warned = .false.
   integer, dimension(:), allocatable, public :: l2g_map
   integer, dimension(:), allocatable, public :: g2l_map

   type, public :: parallel_context
      integer :: rank = 0
      integer :: size = 1
      integer :: local_rank = 0
      integer :: local_size = 1
      logical :: mpi_enabled = .false.
#ifdef USE_MPI
      ! Keep this independent of MPI constants: module/type initialization may
      ! occur before MPI_Init on some MPI implementations.
      integer :: local_comm = -1
#endif
   contains
      procedure :: is_root
      procedure :: owns_work
      procedure :: split_range
      procedure :: barrier
      procedure :: restore_to_default
      procedure :: device_index
      procedure, private :: initialize
      final :: destructor
   end type parallel_context

   type(parallel_context), public :: g_parallel_context

   interface parallel_context
      procedure :: constructor
   end interface parallel_context

   interface assignment(=)
      module procedure assign_parallel_context
   end interface assignment(=)

   public :: get_mpi_range, get_mpi_variables, get_mpi_mapping
   public :: split_range_for_worker, map_local_rank_to_device
   public :: assignment(=)

contains

   !> @brief Construct the runtime context after MPI_Init when MPI is enabled.
   function constructor() result(obj)
      type(parallel_context) :: obj

      call obj%restore_to_default()
      call obj%initialize()
   end function constructor

   !> @brief Copy a context without sharing ownership of its local communicator.
   !> @details Function-result assignment finalizes the temporary immediately
   !> afterwards.  Duplicate the communicator here so that the temporary and
   !> destination can independently release their owned handles.
   subroutine assign_parallel_context(lhs, rhs)
      type(parallel_context), intent(out) :: lhs
      type(parallel_context), intent(in) :: rhs
#ifdef USE_MPI
      logical :: mpi_is_initialized, mpi_is_finalized
      integer :: mpi_status
#endif

      lhs%rank = rhs%rank
      lhs%size = rhs%size
      lhs%local_rank = rhs%local_rank
      lhs%local_size = rhs%local_size
      lhs%mpi_enabled = rhs%mpi_enabled
#ifdef USE_MPI
      lhs%local_comm = -1
      if (rhs%local_comm == -1) return

      call MPI_Initialized(mpi_is_initialized, mpi_status)
      if (.not. mpi_is_initialized .or. mpi_status /= MPI_SUCCESS) return
      call MPI_Finalized(mpi_is_finalized, mpi_status)
      if (mpi_is_finalized .or. mpi_status /= MPI_SUCCESS) return

      call MPI_Comm_dup(rhs%local_comm, lhs%local_comm, mpi_status)
      if (mpi_status /= MPI_SUCCESS) then
         lhs%local_comm = -1
         call g_logger%warning('MPI_Comm_dup failed while copying parallel context; local communicator disabled.', __FILE__, __LINE__)
      end if
#endif
   end subroutine assign_parallel_context

   !> @brief Initialize from MPI_COMM_WORLD and its shared-memory subcommunicator.
   subroutine initialize(this)
      class(parallel_context), intent(inout) :: this
#ifdef USE_MPI
      logical :: mpi_is_initialized
      integer :: mpi_status

      call MPI_Initialized(mpi_is_initialized, mpi_status)
      if (.not. mpi_is_initialized .or. mpi_status /= MPI_SUCCESS) return

      call MPI_Comm_rank(MPI_COMM_WORLD, this%rank, mpi_status)
      if (mpi_status /= MPI_SUCCESS) return
      call MPI_Comm_size(MPI_COMM_WORLD, this%size, mpi_status)
      if (mpi_status /= MPI_SUCCESS) return

      this%mpi_enabled = .true.
      call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, this%rank, MPI_INFO_NULL, this%local_comm, mpi_status)
      if (mpi_status == MPI_SUCCESS) then
         call MPI_Comm_rank(this%local_comm, this%local_rank, mpi_status)
         call MPI_Comm_size(this%local_comm, this%local_size, mpi_status)
      else
         ! A usable world context is still preferable to aborting solely because
         ! node-local identity is unavailable on an older MPI implementation.
         this%local_rank = this%rank
         this%local_size = this%size
         this%local_comm = -1
         call g_logger%warning('MPI_Comm_split_type(MPI_COMM_TYPE_SHARED) failed; using world-local identity.', __FILE__, __LINE__)
      end if
#endif
      call synchronize_compatibility_globals(this)
   end subroutine initialize

   !> @brief Return whether this context owns root-only work/output.
   pure logical function is_root(this)
      class(parallel_context), intent(in) :: this
      is_root = this%rank == 0
   end function is_root

   !> @brief Return whether a one-based item belongs to this context's range.
   pure logical function owns_work(this, item_index, n_items)
      class(parallel_context), intent(in) :: this
      integer, intent(in) :: item_index, n_items
      integer :: first, last, count

      call this%split_range(n_items, first, last, count)
      owns_work = count > 0 .and. item_index >= first .and. item_index <= last
   end function owns_work

   !> @brief Split one-based work items over world ranks.
   pure subroutine split_range(this, n_items, start_idx, end_idx, items_per_rank)
      class(parallel_context), intent(in) :: this
      integer, intent(in) :: n_items
      integer, intent(out) :: start_idx, end_idx, items_per_rank

      call split_range_for_worker(this%rank, this%size, n_items, start_idx, end_idx, items_per_rank)
   end subroutine split_range

   !> @brief Synchronize all world ranks when MPI is active; serial is a no-op.
   subroutine barrier(this)
      class(parallel_context), intent(in) :: this
#ifdef USE_MPI
      integer :: mpi_status
      if (this%mpi_enabled) call MPI_Barrier(MPI_COMM_WORLD, mpi_status)
#endif
   end subroutine barrier

   !> @brief Reset this context and release only its owned local communicator.
   subroutine restore_to_default(this)
      class(parallel_context), intent(inout) :: this
#ifdef USE_MPI
      logical :: mpi_is_initialized, mpi_is_finalized
      integer :: mpi_status, communicator

      communicator = this%local_comm
      call MPI_Initialized(mpi_is_initialized, mpi_status)
      if (mpi_is_initialized .and. mpi_status == MPI_SUCCESS) then
         call MPI_Finalized(mpi_is_finalized, mpi_status)
         if (.not. mpi_is_finalized .and. communicator /= -1) then
            call MPI_Comm_free(communicator, mpi_status)
         end if
      end if
      this%local_comm = -1
#endif
      this%rank = 0
      this%size = 1
      this%local_rank = 0
      this%local_size = 1
      this%mpi_enabled = .false.
   end subroutine restore_to_default

   !> @brief Release context-owned communicators; MPI_Finalize remains program-owned.
   subroutine destructor(this)
      type(parallel_context), intent(inout) :: this
      call this%restore_to_default()
   end subroutine destructor

   !> @brief Map the node-local rank to a hardware-independent device index.
   subroutine device_index(this, available_device_count, index, valid, override)
      class(parallel_context), intent(in) :: this
      integer, intent(in) :: available_device_count
      integer, intent(out) :: index
      logical, intent(out) :: valid
      integer, intent(in), optional :: override

      call map_local_rank_to_device(this%local_rank, available_device_count, index, valid, override)
   end subroutine device_index

   !> @brief Deterministically split an item range for an explicit worker set.
   !> @details Invalid workers and non-positive item counts receive the empty
   !> range `[1,0]`; that keeps collective callers safe under oversubscription.
   pure subroutine split_range_for_worker(worker_rank, worker_size, n_items, start_idx, end_idx, items_per_worker)
      integer, intent(in) :: worker_rank, worker_size, n_items
      integer, intent(out) :: start_idx, end_idx, items_per_worker
      integer :: effective_items, remainder_items

      start_idx = 1
      end_idx = 0
      items_per_worker = 0
      if (worker_size <= 0 .or. worker_rank < 0 .or. worker_rank >= worker_size) return

      effective_items = max(0, n_items)
      items_per_worker = effective_items/worker_size
      remainder_items = mod(effective_items, worker_size)
      if (worker_rank < remainder_items) then
         items_per_worker = items_per_worker + 1
         start_idx = worker_rank*items_per_worker + 1
      else
         start_idx = worker_rank*items_per_worker + remainder_items + 1
      end if
      end_idx = start_idx + items_per_worker - 1
      if (items_per_worker == 0) then
         start_idx = 1
         end_idx = 0
      end if
   end subroutine split_range_for_worker

   !> @brief Select a device by local rank, optionally honoring a programmatic override.
   !> @details The default is round-robin `mod(local_rank, device_count)`. An
   !> invalid local rank, device count, or override is rejected with `valid=.false.`
   !> and `device_index=-1`; no accelerator runtime is queried here.
   pure subroutine map_local_rank_to_device(local_rank, available_device_count, device_index, valid, override)
      integer, intent(in) :: local_rank, available_device_count
      integer, intent(out) :: device_index
      logical, intent(out) :: valid
      integer, intent(in), optional :: override

      device_index = -1
      valid = .false.
      if (local_rank < 0 .or. available_device_count <= 0) return
      if (present(override)) then
         if (override < 0 .or. override >= available_device_count) return
         device_index = override
      else
         device_index = mod(local_rank, available_device_count)
      end if
      valid = .true.
   end subroutine map_local_rank_to_device

   !> @brief Synchronize compatibility globals from the one runtime context.
   subroutine synchronize_compatibility_globals(context)
      class(parallel_context), intent(in) :: context
      rank = context%rank
      numprocs = context%size
      ierr = 0
   end subroutine synchronize_compatibility_globals

   subroutine get_mpi_range(rank_in, n_items, start_idx, end_idx, items_per_process, l2g_map_out, g2l_map_out, region_tag)
      integer, intent(in) :: rank_in, n_items
      integer, intent(out) :: start_idx, end_idx, items_per_process
      integer, dimension(:), allocatable, intent(out), optional :: l2g_map_out
      integer, dimension(:), allocatable, intent(out), optional :: g2l_map_out
      character(len=*), intent(in), optional :: region_tag

      integer :: i_g, i_l
      character(len=160) :: warn_msg
      character(len=32) :: region_name

      region_name = 'work'
      if (present(region_tag)) region_name = trim(region_tag)
      if ((numprocs > n_items) .and. (rank_in == 0) .and. (.not. oversubscribe_warned)) then
         write (warn_msg, '(a,a,a,i0,a,i0,a)') 'MPI oversubscription in ', trim(region_name), &
            '-parallel region (ranks=', numprocs, ', items=', n_items, '). Some ranks will stay idle.'
         call g_logger%warning(trim(warn_msg), __FILE__, __LINE__)
         oversubscribe_warned = .true.
      end if

      call split_range_for_worker(rank_in, numprocs, n_items, start_idx, end_idx, items_per_process)

      if (present(l2g_map_out)) then
         if (allocated(l2g_map_out)) deallocate(l2g_map_out)
         allocate(l2g_map_out(items_per_process))
         l2g_map_out = 0
      end if
      if (present(g2l_map_out)) then
         if (allocated(g2l_map_out)) deallocate(g2l_map_out)
         allocate(g2l_map_out(max(0, n_items)))
         g2l_map_out = 0
      end if

      do i_g = start_idx, end_idx
         i_l = i_g - start_idx + 1
         if (present(l2g_map_out)) l2g_map_out(i_l) = i_g
         if (present(g2l_map_out)) g2l_map_out(i_g) = i_l
      end do
   end subroutine get_mpi_range

   subroutine get_mpi_variables(rank_in, number_of_atoms)
      integer, intent(in) :: number_of_atoms, rank_in
      call get_mpi_range(rank_in, number_of_atoms, start_atom, end_atom, atoms_per_process, region_tag='atom')
      remainder_atoms = max(0, number_of_atoms - numprocs*(number_of_atoms/max(1, numprocs)))
      call get_mpi_mapping(rank_in, number_of_atoms)
   end subroutine get_mpi_variables

   subroutine get_mpi_mapping(rank_in, number_of_atoms)
      integer, intent(in) :: number_of_atoms, rank_in
      call get_mpi_range(rank_in, number_of_atoms, start_atom, end_atom, atoms_per_process, l2g_map, g2l_map, 'atom')
   end subroutine get_mpi_mapping

end module mpi_mod
