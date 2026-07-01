!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: hambuild_cuda_plugin_mod
!
!> @brief
!> iso_c_binding shim for the GPU Hamiltonian-assembly backend (hambuild).
!>
!> Mirrors rsk_cuda_plugin_mod / rsrec_cuda_plugin_mod. Owns an opaque C context
!> (device pointers live on the C side so assembled arrays stay resident and
!> feed the rsrec_cuda recursion backend without a host round-trip). Phase 0
!> only wires context lifecycle + backend selection; the numerical set_*/onsite/
!> bulk/local/ccor entry points are added phase by phase.
!>
!> The whole module compiles unconditionally. When ENABLE_CUDA_PLUGIN is OFF the
!> C symbols resolve against cuda/hambuild_stub.c (which returns a NULL context),
!> so requesting control%gpu_hambuild=.true. on a CPU-only binary fails cleanly.
!------------------------------------------------------------------------------
module hambuild_cuda_plugin_mod

   use iso_c_binding
   use precision_mod, only: rp
   use logger_mod, only: g_logger
   implicit none

   private
   public :: hambuild_cuda_backend
   public :: hambuild_cuda_plugin_compiled

   !> Backend selector values (must match enum in cuda/hambuild.h).
   integer, parameter, public :: HAMBUILD_BACKEND_CPU = 0
   integer, parameter, public :: HAMBUILD_BACKEND_GPU = 1

   type :: hambuild_cuda_backend
      type(c_ptr) :: ctx = c_null_ptr
      integer(c_int) :: nb = 0
      integer(c_int) :: norb = 0
      integer(c_int) :: nnmax = 0
      integer(c_int) :: ntype = 0
      integer(c_int) :: nmax = 0
   contains
      procedure :: ensure_context
      procedure :: set_backend
      procedure :: destroy
   end type hambuild_cuda_backend

#ifdef USE_CUDA_PLUGIN
   interface
      function hambuild_cuda_create(nb, norb, nnmax, ntype, nmax, device) &
         bind(C, name='hambuild_cuda_create')
         import :: c_int, c_ptr
         integer(c_int), value :: nb, norb, nnmax, ntype, nmax, device
         type(c_ptr) :: hambuild_cuda_create
      end function hambuild_cuda_create

      subroutine hambuild_cuda_destroy(ctx) bind(C, name='hambuild_cuda_destroy')
         import :: c_ptr
         type(c_ptr), value :: ctx
      end subroutine hambuild_cuda_destroy

      function hambuild_cuda_last_error() bind(C, name='hambuild_cuda_last_error')
         import :: c_ptr
         type(c_ptr) :: hambuild_cuda_last_error
      end function hambuild_cuda_last_error

      function hambuild_cuda_set_backend(ctx, backend) &
         bind(C, name='hambuild_cuda_set_backend')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx
         integer(c_int), value :: backend
         integer(c_int) :: hambuild_cuda_set_backend
      end function hambuild_cuda_set_backend
   end interface
#endif

contains

   !> Whether this binary was compiled with the GPU hambuild backend available.
   logical function hambuild_cuda_plugin_compiled()
#ifdef USE_CUDA_PLUGIN
      hambuild_cuda_plugin_compiled = .true.
#else
      hambuild_cuda_plugin_compiled = .false.
#endif
   end function hambuild_cuda_plugin_compiled

   !> Create (or reuse) the device context for the given dimensions.
   subroutine ensure_context(this, nb, norb, nnmax, ntype, nmax)
      class(hambuild_cuda_backend), intent(inout) :: this
      integer, intent(in) :: nb, norb, nnmax, ntype, nmax
#ifdef USE_CUDA_PLUGIN
      if (c_associated(this%ctx)) then
         if (this%nb == nb .and. this%norb == norb .and. &
             this%nnmax == nnmax .and. this%ntype == ntype .and. &
             this%nmax == nmax) return
         call this%destroy()
      end if

      this%ctx = hambuild_cuda_create(int(nb, c_int), int(norb, c_int), &
         int(nnmax, c_int), int(ntype, c_int), int(nmax, c_int), 0_c_int)
      if (.not. c_associated(this%ctx)) then
         call g_logger%fatal('Failed to create CUDA hambuild plugin context.', &
            __FILE__, __LINE__)
      end if
      this%nb = int(nb, c_int)
      this%norb = int(norb, c_int)
      this%nnmax = int(nnmax, c_int)
      this%ntype = int(ntype, c_int)
      this%nmax = int(nmax, c_int)
#else
      call g_logger%fatal('control%gpu_hambuild=.true. requested, but this '// &
         'executable was built without ENABLE_CUDA_PLUGIN.', __FILE__, __LINE__)
#endif
   end subroutine ensure_context

   !> Select the CPU-reference or GPU backend for assembly kernels.
   subroutine set_backend(this, backend)
      class(hambuild_cuda_backend), intent(inout) :: this
      integer, intent(in) :: backend
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr
      if (.not. c_associated(this%ctx)) then
         call g_logger%fatal('hambuild set_backend called before ensure_context.', &
            __FILE__, __LINE__)
      end if
      ierr = hambuild_cuda_set_backend(this%ctx, int(backend, c_int))
      if (ierr /= 0_c_int) then
         call g_logger%fatal('hambuild_cuda_set_backend failed.', __FILE__, __LINE__)
      end if
#else
      call g_logger%fatal('control%gpu_hambuild=.true. requested, but this '// &
         'executable was built without ENABLE_CUDA_PLUGIN.', __FILE__, __LINE__)
#endif
   end subroutine set_backend

   !> Release the device context.
   subroutine destroy(this)
      class(hambuild_cuda_backend), intent(inout) :: this
#ifdef USE_CUDA_PLUGIN
      if (c_associated(this%ctx)) call hambuild_cuda_destroy(this%ctx)
#endif
      this%ctx = c_null_ptr
      this%nb = 0
      this%norb = 0
      this%nnmax = 0
      this%ntype = 0
      this%nmax = 0
   end subroutine destroy

end module hambuild_cuda_plugin_mod
