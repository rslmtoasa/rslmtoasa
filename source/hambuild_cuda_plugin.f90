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
      procedure :: set_constants
      procedure :: set_potential_onsite
      procedure :: onsite
      procedure :: set_geometry
      procedure :: set_geometry_vet
      procedure :: build_geometry_maps
      procedure :: get_geometry_maps
      procedure :: set_potential_bulk
      procedure :: set_sbar
      procedure :: bulk
      procedure :: set_local_sites
      procedure :: build_local_geometry_maps
      procedure :: local
      procedure :: set_ccor
      procedure :: ccor_bulk
      procedure :: ccor_local
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

      function hambuild_cuda_set_constants(ctx, v, vc, lx, ly, lz) &
         bind(C, name='hambuild_cuda_set_constants')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx
         type(c_ptr), value :: v, vc, lx, ly, lz
         integer(c_int) :: hambuild_cuda_set_constants
      end function hambuild_cuda_set_constants

      function hambuild_cuda_set_potential_onsite(ctx, xi_p, xi_d, rac, obx0, &
         obx1, cx, cex, mom, lmom, orb_pol) &
         bind(C, name='hambuild_cuda_set_potential_onsite')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx
         type(c_ptr), value :: xi_p, xi_d, rac, obx0, obx1, cx, cex, mom, lmom
         integer(c_int), value :: orb_pol
         integer(c_int) :: hambuild_cuda_set_potential_onsite
      end function hambuild_cuda_set_potential_onsite

      function hambuild_cuda_onsite(ctx, lsham, obarm, enim) &
         bind(C, name='hambuild_cuda_onsite')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx
         type(c_ptr), value :: lsham, obarm, enim
         integer(c_int) :: hambuild_cuda_onsite
      end function hambuild_cuda_onsite

      function hambuild_cuda_set_geometry(ctx, cr, num, iz, nn, atlist, kk, &
         nn_max, ndi, alat, r2, pbc) bind(C, name='hambuild_cuda_set_geometry')
         import :: c_int, c_ptr, c_double
         type(c_ptr), value :: ctx
         type(c_ptr), value :: cr, num, iz, nn, atlist
         integer(c_int), value :: kk, nn_max, ndi, pbc
         real(c_double), value :: alat, r2
         integer(c_int) :: hambuild_cuda_set_geometry
      end function hambuild_cuda_set_geometry

      function hambuild_cuda_set_geometry_vet(ctx, vet) &
         bind(C, name='hambuild_cuda_set_geometry_vet')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx
         type(c_ptr), value :: vet
         integer(c_int) :: hambuild_cuda_set_geometry_vet
      end function hambuild_cuda_set_geometry_vet

      function hambuild_cuda_build_geometry_maps(ctx) &
         bind(C, name='hambuild_cuda_build_geometry_maps')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx
         integer(c_int) :: hambuild_cuda_build_geometry_maps
      end function hambuild_cuda_build_geometry_maps

      function hambuild_cuda_get_geometry_maps(ctx, valid, shell, ino, vet) &
         bind(C, name='hambuild_cuda_get_geometry_maps')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx
         type(c_ptr), value :: valid, shell, ino, vet
         integer(c_int) :: hambuild_cuda_get_geometry_maps
      end function hambuild_cuda_get_geometry_maps

      function hambuild_cuda_set_potential_bulk(ctx, wx0, wx1, cx0, cx1, cex0, &
         cex1, mom, q_ss, theta_ss) &
         bind(C, name='hambuild_cuda_set_potential_bulk')
         import :: c_int, c_ptr, c_double
         type(c_ptr), value :: ctx
         type(c_ptr), value :: wx0, wx1, cx0, cx1, cex0, cex1, mom, q_ss
         real(c_double), value :: theta_ss
         integer(c_int) :: hambuild_cuda_set_potential_bulk
      end function hambuild_cuda_set_potential_bulk

      function hambuild_cuda_set_sbar(ctx, sbar, nm_store, ntot) &
         bind(C, name='hambuild_cuda_set_sbar')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx
         type(c_ptr), value :: sbar
         integer(c_int), value :: nm_store, ntot
         integer(c_int) :: hambuild_cuda_set_sbar
      end function hambuild_cuda_set_sbar

      function hambuild_cuda_bulk(ctx, hoh, ee, hxc, eeo, eeoee) &
         bind(C, name='hambuild_cuda_bulk')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx
         integer(c_int), value :: hoh
         type(c_ptr), value :: ee, hxc, eeo, eeoee
         integer(c_int) :: hambuild_cuda_bulk
      end function hambuild_cuda_bulk

      function hambuild_cuda_set_local_sites(ctx, site_list, nmax) &
         bind(C, name='hambuild_cuda_set_local_sites')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx
         type(c_ptr), value :: site_list
         integer(c_int), value :: nmax
         integer(c_int) :: hambuild_cuda_set_local_sites
      end function hambuild_cuda_set_local_sites

      function hambuild_cuda_build_local_geometry_maps(ctx) &
         bind(C, name='hambuild_cuda_build_local_geometry_maps')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx
         integer(c_int) :: hambuild_cuda_build_local_geometry_maps
      end function hambuild_cuda_build_local_geometry_maps

      function hambuild_cuda_local(ctx, hoh, hall, hallo) &
         bind(C, name='hambuild_cuda_local')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx
         integer(c_int), value :: hoh
         type(c_ptr), value :: hall, hallo
         integer(c_int) :: hambuild_cuda_local
      end function hambuild_cuda_local

      function hambuild_cuda_set_ccor(ctx, ccd, lambda, sdot, avw) &
         bind(C, name='hambuild_cuda_set_ccor')
         import :: c_int, c_ptr, c_double
         type(c_ptr), value :: ctx
         type(c_ptr), value :: ccd, lambda, sdot
         real(c_double), value :: avw
         integer(c_int) :: hambuild_cuda_set_ccor
      end function hambuild_cuda_set_ccor

      function hambuild_cuda_ccor_bulk(ctx, eecc) &
         bind(C, name='hambuild_cuda_ccor_bulk')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx
         type(c_ptr), value :: eecc
         integer(c_int) :: hambuild_cuda_ccor_bulk
      end function hambuild_cuda_ccor_bulk

      function hambuild_cuda_ccor_local(ctx, hallcc) &
         bind(C, name='hambuild_cuda_ccor_local')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx
         type(c_ptr), value :: hallcc
         integer(c_int) :: hambuild_cuda_ccor_local
      end function hambuild_cuda_ccor_local
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

   !> Upload the basis-dependent constants (hcpx v/vc + cartesian L operators).
   !> All arrays are (norb x norb) complex, column-major, contiguous targets.
   subroutine set_constants(this, v, vc, lx, ly, lz)
      class(hambuild_cuda_backend), intent(inout) :: this
      complex(rp), dimension(:, :), intent(in), target, contiguous :: v, vc, lx, ly, lz
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr
      ierr = hambuild_cuda_set_constants(this%ctx, c_loc(v), c_loc(vc), &
         c_loc(lx), c_loc(ly), c_loc(lz))
      if (ierr /= 0_c_int) call fail('hambuild_cuda_set_constants')
#else
      call plugin_absent()
#endif
   end subroutine set_constants

   !> Upload the per-type on-site potential inputs. Arrays are contiguous
   !> targets, laid out per type as documented in cuda/hambuild.h. rac and lmom
   !> must already carry the CPU-precomputed polarization factors (see caller).
   subroutine set_potential_onsite(this, xi_p, xi_d, rac, obx0, obx1, cx, cex, &
      mom, lmom, orb_pol)
      class(hambuild_cuda_backend), intent(inout) :: this
      real(rp), dimension(:, :), intent(in), target, contiguous :: xi_p, xi_d, rac
      complex(rp), dimension(:, :), intent(in), target, contiguous :: obx0, obx1
      complex(rp), dimension(:, :, :), intent(in), target, contiguous :: cx, cex
      real(rp), dimension(:, :), intent(in), target, contiguous :: mom, lmom
      logical, intent(in) :: orb_pol
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr, op
      op = 0_c_int; if (orb_pol) op = 1_c_int
      ierr = hambuild_cuda_set_potential_onsite(this%ctx, c_loc(xi_p), &
         c_loc(xi_d), c_loc(rac), c_loc(obx0), c_loc(obx1), c_loc(cx), &
         c_loc(cex), c_loc(mom), c_loc(lmom), op)
      if (ierr /= 0_c_int) call fail('hambuild_cuda_set_potential_onsite')
#else
      call plugin_absent()
#endif
   end subroutine set_potential_onsite

   !> Build the on-site blocks on the GPU and copy them back. Pass unallocated
   !> optionals to skip a block; here all three are required outputs. Arrays are
   !> (nb x nb x ntype) complex, contiguous targets.
   subroutine onsite(this, lsham, obarm, enim)
      class(hambuild_cuda_backend), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(inout), target, contiguous :: lsham, obarm, enim
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr
      ierr = hambuild_cuda_onsite(this%ctx, c_loc(lsham), c_loc(obarm), c_loc(enim))
      if (ierr /= 0_c_int) call fail('hambuild_cuda_onsite')
#else
      call plugin_absent()
#endif
   end subroutine onsite

   !> Upload the static geometry (once per run; geometry is fixed across SCF).
   !> Integer arrays are the raw 1-based Fortran indices; the kernels convert to
   !> 0-based internally. cr is (3,kk), nn is (ndi,nn_max), etc.
   subroutine set_geometry(this, cr, num, iz, nn, atlist, nn_max, ndi, alat, r2, pbc)
      class(hambuild_cuda_backend), intent(inout) :: this
      real(rp), dimension(:, :), intent(in), target, contiguous :: cr
      integer(c_int), dimension(:), intent(in), target, contiguous :: num, iz, atlist
      integer(c_int), dimension(:, :), intent(in), target, contiguous :: nn
      integer, intent(in) :: nn_max, ndi
      real(rp), intent(in) :: alat, r2
      logical, intent(in) :: pbc
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr, pbc_i
      pbc_i = 0_c_int; if (pbc) pbc_i = 1_c_int
      ierr = hambuild_cuda_set_geometry(this%ctx, c_loc(cr), c_loc(num), &
         c_loc(iz), c_loc(nn), c_loc(atlist), int(size(cr, 2), c_int), &
         int(nn_max, c_int), int(ndi, c_int), real(alat, c_double), &
         real(r2, c_double), pbc_i)
      if (ierr /= 0_c_int) call fail('hambuild_cuda_set_geometry')
#else
      call plugin_absent()
#endif
   end subroutine set_geometry

   !> Upload precomputed displacement vectors (pbc / host-computed vet). vet is
   !> (3, nn_max, ntype) contiguous; when set, the map builder uses these.
   subroutine set_geometry_vet(this, vet)
      class(hambuild_cuda_backend), intent(inout) :: this
      real(rp), dimension(:, :, :), intent(in), target, contiguous :: vet
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr
      ierr = hambuild_cuda_set_geometry_vet(this%ctx, c_loc(vet))
      if (ierr /= 0_c_int) call fail('hambuild_cuda_set_geometry_vet')
#else
      call plugin_absent()
#endif
   end subroutine set_geometry_vet

   !> Build the resident neigh_map + vet tables from the uploaded geometry.
   subroutine build_geometry_maps(this)
      class(hambuild_cuda_backend), intent(inout) :: this
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr
      ierr = hambuild_cuda_build_geometry_maps(this%ctx)
      if (ierr /= 0_c_int) call fail('hambuild_cuda_build_geometry_maps')
#else
      call plugin_absent()
#endif
   end subroutine build_geometry_maps

   !> Copy the built maps back to the host (validation / inspection). Arrays are
   !> (nn_max, ntype) for valid/shell/ino and (3, nn_max, ntype) for vet.
   subroutine get_geometry_maps(this, valid, shell, ino, vet)
      class(hambuild_cuda_backend), intent(inout) :: this
      integer(c_int), dimension(:, :), intent(inout), target, contiguous :: valid, shell, ino
      real(rp), dimension(:, :, :), intent(inout), target, contiguous :: vet
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr
      ierr = hambuild_cuda_get_geometry_maps(this%ctx, c_loc(valid), &
         c_loc(shell), c_loc(ino), c_loc(vet))
      if (ierr /= 0_c_int) call fail('hambuild_cuda_get_geometry_maps')
#else
      call plugin_absent()
#endif
   end subroutine get_geometry_maps

   !> Upload per-type potential vectors for the neighbour Hamiltonian (ham0m_nc).
   !> wx0..cex1 are (norb, ntype) complex; mom is (3, ntype); q_ss is (3).
   subroutine set_potential_bulk(this, wx0, wx1, cx0, cx1, cex0, cex1, mom, q_ss, theta_ss)
      class(hambuild_cuda_backend), intent(inout) :: this
      complex(rp), dimension(:, :), intent(in), target, contiguous :: wx0, wx1, cx0, cx1, cex0, cex1
      real(rp), dimension(:, :), intent(in), target, contiguous :: mom
      real(rp), dimension(3), intent(in), target :: q_ss
      real(rp), intent(in) :: theta_ss
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr
      ierr = hambuild_cuda_set_potential_bulk(this%ctx, c_loc(wx0), c_loc(wx1), &
         c_loc(cx0), c_loc(cx1), c_loc(cex0), c_loc(cex1), c_loc(mom), &
         c_loc(q_ss), real(theta_ss, c_double))
      if (ierr /= 0_c_int) call fail('hambuild_cuda_set_potential_bulk')
#else
      call plugin_absent()
#endif
   end subroutine set_potential_bulk

   !> Upload the structure constants sbar (refreshed per SCF). sbar is
   !> (norb, norb, nm_store, ntot) complex; only real(sbar) is used on device.
   subroutine set_sbar(this, sbar)
      class(hambuild_cuda_backend), intent(inout) :: this
      complex(rp), dimension(:, :, :, :), intent(in), target, contiguous :: sbar
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr
      ierr = hambuild_cuda_set_sbar(this%ctx, c_loc(sbar), &
         int(size(sbar, 3), c_int), int(size(sbar, 4), c_int))
      if (ierr /= 0_c_int) call fail('hambuild_cuda_set_sbar')
#else
      call plugin_absent()
#endif
   end subroutine set_sbar

   !> Build the bulk neighbour Hamiltonian on-device and copy ee/hxc (and, when
   !> hoh, eeo/eeoee) back. Arrays are (nb, nb, nn_max, ntype) complex.
   subroutine bulk(this, hoh, ee, hxc, eeo, eeoee)
      class(hambuild_cuda_backend), intent(inout) :: this
      logical, intent(in) :: hoh
      complex(rp), dimension(:, :, :, :), intent(inout), target, contiguous :: ee, hxc, eeo, eeoee
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr, hoh_i
      hoh_i = 0_c_int; if (hoh) hoh_i = 1_c_int
      ierr = hambuild_cuda_bulk(this%ctx, hoh_i, c_loc(ee), c_loc(hxc), &
         c_loc(eeo), c_loc(eeoee))
      if (ierr /= 0_c_int) call fail('hambuild_cuda_bulk')
#else
      call plugin_absent()
#endif
   end subroutine bulk

   !> Upload the local (impurity interaction-zone) site list (1-based cluster
   !> indices, length nmax) and build the local geometry maps on-device.
   subroutine set_local_sites(this, site_list)
      class(hambuild_cuda_backend), intent(inout) :: this
      integer(c_int), dimension(:), intent(in), target, contiguous :: site_list
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr
      ierr = hambuild_cuda_set_local_sites(this%ctx, c_loc(site_list), &
         int(size(site_list), c_int))
      if (ierr /= 0_c_int) call fail('hambuild_cuda_set_local_sites')
#else
      call plugin_absent()
#endif
   end subroutine set_local_sites

   subroutine build_local_geometry_maps(this)
      class(hambuild_cuda_backend), intent(inout) :: this
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr
      ierr = hambuild_cuda_build_local_geometry_maps(this%ctx)
      if (ierr /= 0_c_int) call fail('hambuild_cuda_build_local_geometry_maps')
#else
      call plugin_absent()
#endif
   end subroutine build_local_geometry_maps

   !> Build the local neighbour Hamiltonian on-device and copy hall (and, when
   !> hoh, hallo) back. Arrays are (nb, nb, nn_max, nmax) complex.
   subroutine local(this, hoh, hall, hallo)
      class(hambuild_cuda_backend), intent(inout) :: this
      logical, intent(in) :: hoh
      complex(rp), dimension(:, :, :, :), intent(inout), target, contiguous :: hall, hallo
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr, hoh_i
      hoh_i = 0_c_int; if (hoh) hoh_i = 1_c_int
      ierr = hambuild_cuda_local(this%ctx, hoh_i, c_loc(hall), c_loc(hallo))
      if (ierr /= 0_c_int) call fail('hambuild_cuda_local')
#else
      call plugin_absent()
#endif
   end subroutine local

   !> Upload CCOR inputs: ccd (norb,0:2,ntype), lambda (ntype,ntype), sdot
   !> (norb,norb,nm_store,ntot; real part used), and the avw normalization scalar.
   subroutine set_ccor(this, ccd, lambda, sdot, avw)
      class(hambuild_cuda_backend), intent(inout) :: this
      real(rp), dimension(:, :, :), intent(in), target, contiguous :: ccd
      real(rp), dimension(:, :), intent(in), target, contiguous :: lambda
      complex(rp), dimension(:, :, :, :), intent(in), target, contiguous :: sdot
      real(rp), intent(in) :: avw
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr
      ierr = hambuild_cuda_set_ccor(this%ctx, c_loc(ccd), c_loc(lambda), &
         c_loc(sdot), real(avw, c_double))
      if (ierr /= 0_c_int) call fail('hambuild_cuda_set_ccor')
#else
      call plugin_absent()
#endif
   end subroutine set_ccor

   !> Build H_cc over the bulk sites into eecc (nb,nb,nn_max,ntype).
   subroutine ccor_bulk(this, eecc)
      class(hambuild_cuda_backend), intent(inout) :: this
      complex(rp), dimension(:, :, :, :), intent(inout), target, contiguous :: eecc
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr
      ierr = hambuild_cuda_ccor_bulk(this%ctx, c_loc(eecc))
      if (ierr /= 0_c_int) call fail('hambuild_cuda_ccor_bulk')
#else
      call plugin_absent()
#endif
   end subroutine ccor_bulk

   !> Build H_cc over the local sites into hallcc (nb,nb,nn_max,nmax).
   subroutine ccor_local(this, hallcc)
      class(hambuild_cuda_backend), intent(inout) :: this
      complex(rp), dimension(:, :, :, :), intent(inout), target, contiguous :: hallcc
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: ierr
      ierr = hambuild_cuda_ccor_local(this%ctx, c_loc(hallcc))
      if (ierr /= 0_c_int) call fail('hambuild_cuda_ccor_local')
#else
      call plugin_absent()
#endif
   end subroutine ccor_local

#ifdef USE_CUDA_PLUGIN
   !> Report a failed C entry point, appending the backend error string.
   subroutine fail(what)
      character(len=*), intent(in) :: what
      character(kind=c_char), pointer :: cstr(:)
      character(len=256) :: msg
      integer :: i
      type(c_ptr) :: cptr
      msg = ''
      cptr = hambuild_cuda_last_error()
      if (c_associated(cptr)) then
         call c_f_pointer(cptr, cstr, [256])
         do i = 1, 256
            if (cstr(i) == c_null_char) exit
            msg(i:i) = cstr(i)
         end do
      end if
      call g_logger%fatal(trim(what)//' failed: '//trim(msg), __FILE__, __LINE__)
   end subroutine fail
#endif

   subroutine plugin_absent()
      call g_logger%fatal('control%gpu_hambuild=.true. requested, but this '// &
         'executable was built without ENABLE_CUDA_PLUGIN.', __FILE__, __LINE__)
   end subroutine plugin_absent

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
