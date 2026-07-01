module rsk_cuda_plugin_mod

   use iso_c_binding
   use precision_mod, only: rp
   use logger_mod, only: g_logger
   implicit none

   private
   public :: rsk_cuda_backend
   public :: rsk_cuda_plugin_compiled

   type :: rsk_cuda_backend
      type(c_ptr) :: ctx = c_null_ptr
      integer(c_int) :: nb = 0
      integer(c_int) :: nsite = 0
      integer(c_int) :: ntype = 0
      integer(c_int) :: nnmax = 0
   contains
      procedure :: ensure_context
      procedure :: destroy
   end type rsk_cuda_backend

#ifdef USE_CUDA_PLUGIN
   interface
      function rsk_cuda_create(nb, nsite, ntype, nnmax, device) bind(C, name='rsk_cuda_create')
         import :: c_int, c_ptr
         integer(c_int), value :: nb, nsite, ntype, nnmax, device
         type(c_ptr) :: rsk_cuda_create
      end function rsk_cuda_create

      subroutine rsk_cuda_destroy(ctx) bind(C, name='rsk_cuda_destroy')
         import :: c_ptr
         type(c_ptr), value :: ctx
      end subroutine rsk_cuda_destroy

      function rsk_cuda_last_error() bind(C, name='rsk_cuda_last_error')
         import :: c_ptr
         type(c_ptr) :: rsk_cuda_last_error
      end function rsk_cuda_last_error
   end interface
#endif

contains

   logical function rsk_cuda_plugin_compiled()
#ifdef USE_CUDA_PLUGIN
      rsk_cuda_plugin_compiled = .true.
#else
      rsk_cuda_plugin_compiled = .false.
#endif
   end function rsk_cuda_plugin_compiled

   subroutine ensure_context(this, nb, nsite, ntype, nnmax)
      class(rsk_cuda_backend), intent(inout) :: this
      integer, intent(in) :: nb, nsite, ntype, nnmax
#ifdef USE_CUDA_PLUGIN
      if (c_associated(this%ctx)) then
         if (this%nb == nb .and. this%nsite == nsite .and. &
             this%ntype == ntype .and. this%nnmax == nnmax) return
         call this%destroy()
      end if

      this%ctx = rsk_cuda_create(int(nb, c_int), int(nsite, c_int), &
         int(ntype, c_int), int(nnmax, c_int), 0_c_int)
      if (.not. c_associated(this%ctx)) then
         call g_logger%fatal('Failed to create CUDA k-space plugin context.', &
            __FILE__, __LINE__)
      end if
      this%nb = int(nb, c_int)
      this%nsite = int(nsite, c_int)
      this%ntype = int(ntype, c_int)
      this%nnmax = int(nnmax, c_int)
#else
      call g_logger%fatal('gpu_plugin=.true. requested for k-space, but this executable '// &
         'was built without ENABLE_CUDA_PLUGIN.', __FILE__, __LINE__)
#endif
   end subroutine ensure_context

   subroutine destroy(this)
      class(rsk_cuda_backend), intent(inout) :: this
#ifdef USE_CUDA_PLUGIN
      if (c_associated(this%ctx)) call rsk_cuda_destroy(this%ctx)
#endif
      this%ctx = c_null_ptr
      this%nb = 0
      this%nsite = 0
      this%ntype = 0
      this%nnmax = 0
   end subroutine destroy

end module rsk_cuda_plugin_mod
