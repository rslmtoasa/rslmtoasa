!------------------------------------------------------------------------------
!> Explicit ownership of fractional reciprocal k-points.
!> `points` and `weights` always contain the points owned by this execution
!> context.  A replicated workset therefore has nk_local == nk_global; a
!> distributed workset has one contiguous tile on each rank.  Empty ranks use
!> global_start=1, global_end=0 and allocated zero-length point/weight/maps.
module kpoint_workset_mod
   use precision_mod, only: rp
   use mpi_mod, only: parallel_context
   use logger_mod, only: g_logger
   implicit none
   private

   type, public :: kpoint_workset
      integer :: nk_global = 0
      integer :: nk_local = 0
      integer :: global_start = 1
      integer :: global_end = 0
      logical :: distributed = .false.
      real(rp), allocatable :: points(:, :)       ! (3,nk_local), fractional
      real(rp), allocatable :: weights(:)         ! (nk_local)
      integer, allocatable :: local_to_global(:)  ! (nk_local)
      integer, allocatable :: global_to_local(:)  ! (nk_global), 0 if remote
   contains
      procedure :: restore_to_default
      procedure :: validate
      procedure :: local_index
      procedure :: global_index
      procedure :: weight_sum
      procedure :: fold
      procedure :: shifted
      procedure :: select_tile
      final :: destructor
   end type kpoint_workset

   public :: make_kpoint_workset, make_replicated_kpoint_workset

contains

   !> Construct from a complete caller list.  The constructor copies precisely
   !> the selected owned tile, so there is no second mutable full-mesh copy.
   function make_kpoint_workset(points, weights, context, distributed) result(workset)
      real(rp), intent(in) :: points(:, :), weights(:)
      type(parallel_context), intent(in) :: context
      logical, intent(in) :: distributed
      type(kpoint_workset) :: workset

      if (size(points, 1) /= 3 .or. size(weights) /= size(points, 2)) then
         call g_logger%fatal('make_kpoint_workset: points must have shape (3,nk) and match weights.', __FILE__, __LINE__)
      end if
      workset%nk_global = size(weights)
      workset%distributed = distributed
      call workset%select_tile(points, weights, context, distributed)
      call workset%validate()
   end function make_kpoint_workset

   !> Arbitrary caller-owned points are replicated unless the caller explicitly
   !> requests otherwise through make_kpoint_workset.
   function make_replicated_kpoint_workset(points, weights, context) result(workset)
      real(rp), intent(in) :: points(:, :), weights(:)
      type(parallel_context), intent(in) :: context
      type(kpoint_workset) :: workset

      workset = make_kpoint_workset(points, weights, context, .false.)
   end function make_replicated_kpoint_workset

   subroutine restore_to_default(this)
      class(kpoint_workset), intent(inout) :: this
      if (allocated(this%points)) deallocate(this%points)
      if (allocated(this%weights)) deallocate(this%weights)
      if (allocated(this%local_to_global)) deallocate(this%local_to_global)
      if (allocated(this%global_to_local)) deallocate(this%global_to_local)
      this%nk_global = 0
      this%nk_local = 0
      this%global_start = 1
      this%global_end = 0
      this%distributed = .false.
   end subroutine restore_to_default

   subroutine destructor(this)
      type(kpoint_workset), intent(inout) :: this
      call this%restore_to_default()
   end subroutine destructor

   subroutine validate(this)
      class(kpoint_workset), intent(in) :: this
      integer :: ik

      if (this%nk_global < 0 .or. this%nk_local < 0) then
         call g_logger%fatal('kpoint_workset%validate: negative point count.', __FILE__, __LINE__)
      end if
      if (.not. allocated(this%points) .or. .not. allocated(this%weights) .or. &
          .not. allocated(this%local_to_global) .or. .not. allocated(this%global_to_local)) then
         call g_logger%fatal('kpoint_workset%validate: all ownership arrays must be allocated.', __FILE__, __LINE__)
      end if
      if (size(this%points, 1) /= 3 .or. size(this%points, 2) /= this%nk_local .or. &
          size(this%weights) /= this%nk_local .or. size(this%local_to_global) /= this%nk_local .or. &
          size(this%global_to_local) /= this%nk_global) then
         call g_logger%fatal('kpoint_workset%validate: inconsistent array shape.', __FILE__, __LINE__)
      end if
      if (this%nk_local == 0) then
         if (this%global_start /= 1 .or. this%global_end /= 0) then
            call g_logger%fatal('kpoint_workset%validate: empty worksets use global range [1,0].', __FILE__, __LINE__)
         end if
      else
         if (this%global_start < 1 .or. this%global_end > this%nk_global .or. &
             this%global_end-this%global_start+1 /= this%nk_local) then
            call g_logger%fatal('kpoint_workset%validate: invalid contiguous global range.', __FILE__, __LINE__)
         end if
      end if
      do ik = 1, this%nk_local
         if (this%local_to_global(ik) < 1 .or. this%local_to_global(ik) > this%nk_global .or. &
             this%global_to_local(this%local_to_global(ik)) /= ik) then
            call g_logger%fatal('kpoint_workset%validate: local/global maps are inconsistent.', __FILE__, __LINE__)
         end if
      end do
   end subroutine validate

   pure integer function local_index(this, global_index) result(index)
      class(kpoint_workset), intent(in) :: this
      integer, intent(in) :: global_index
      index = 0
      if (global_index >= 1 .and. global_index <= this%nk_global) index = this%global_to_local(global_index)
   end function local_index

   pure integer function global_index(this, local_index) result(index)
      class(kpoint_workset), intent(in) :: this
      integer, intent(in) :: local_index
      index = 0
      if (local_index >= 1 .and. local_index <= this%nk_local) index = this%local_to_global(local_index)
   end function global_index

   pure real(rp) function weight_sum(this) result(sum_weights)
      class(kpoint_workset), intent(in) :: this
      sum_weights = sum(this%weights)
   end function weight_sum

   !> Fold in the historic [-0.5,0.5) fractional reciprocal convention.
   subroutine fold(this)
      class(kpoint_workset), intent(inout) :: this
      if (this%nk_local > 0) this%points = this%points - floor(this%points + 0.5_rp)
   end subroutine fold

   !> Return k+q folded to the same convention without modifying this workset.
   function shifted(this, q_point) result(kq_workset)
      class(kpoint_workset), intent(in) :: this
      real(rp), intent(in) :: q_point(3)
      type(kpoint_workset) :: kq_workset

      kq_workset%nk_global = this%nk_global
      kq_workset%nk_local = this%nk_local
      kq_workset%global_start = this%global_start
      kq_workset%global_end = this%global_end
      kq_workset%distributed = this%distributed
      allocate(kq_workset%points(3, this%nk_local), kq_workset%weights(this%nk_local), &
               kq_workset%local_to_global(this%nk_local), kq_workset%global_to_local(this%nk_global))
      kq_workset%points = this%points + spread(q_point, dim=2, ncopies=this%nk_local)
      kq_workset%points = kq_workset%points - floor(kq_workset%points + 0.5_rp)
      kq_workset%weights = this%weights
      kq_workset%local_to_global = this%local_to_global
      kq_workset%global_to_local = this%global_to_local
   end function shifted

   !> Select a contiguous tile from a complete list.  This public boundary is
   !> intentionally the only point where global point storage is accepted.
   subroutine select_tile(this, all_points, all_weights, context, distributed)
      class(kpoint_workset), intent(inout) :: this
      real(rp), intent(in) :: all_points(:, :), all_weights(:)
      type(parallel_context), intent(in) :: context
      logical, intent(in) :: distributed
      integer :: ik

      if (size(all_points, 1) /= 3 .or. size(all_weights) /= size(all_points, 2)) then
         call g_logger%fatal('kpoint_workset%select_tile: invalid complete point list.', __FILE__, __LINE__)
      end if
      call this%restore_to_default()
      this%nk_global = size(all_weights)
      this%distributed = distributed
      if (distributed) then
         call context%split_range(this%nk_global, this%global_start, this%global_end, this%nk_local)
      else
         this%nk_local = this%nk_global
         if (this%nk_local > 0) then
            this%global_start = 1; this%global_end = this%nk_global
         end if
      end if
      allocate(this%points(3, this%nk_local), this%weights(this%nk_local), &
               this%local_to_global(this%nk_local), this%global_to_local(this%nk_global))
      this%global_to_local = 0
      do ik = 1, this%nk_local
         this%local_to_global(ik) = this%global_start + ik - 1
         this%global_to_local(this%local_to_global(ik)) = ik
         this%points(:, ik) = all_points(:, this%local_to_global(ik))
         this%weights(ik) = all_weights(this%local_to_global(ik))
      end do
   end subroutine select_tile
end module kpoint_workset_mod
