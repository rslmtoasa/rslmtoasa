module rsrec_cuda_plugin_mod

   use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, &
      c_double_complex, c_int, c_null_char, c_null_ptr, c_ptr, c_f_pointer, c_loc
   use precision_mod, only: rp
   use logger_mod, only: g_logger
   use sparse_mod, only: sparse
   implicit none

   private

   integer, parameter, public :: gpu_backend_csr = 0
   integer, parameter, public :: gpu_backend_bsr = 1
   integer, parameter, public :: gpu_backend_fft = 2
   integer, parameter, public :: gpu_backend_conv = 3

   public :: rsrec_cuda_backend
   public :: rsrec_cuda_plugin_compiled
   public :: decode_gpu_backend

   type :: rsrec_cuda_backend
      type(c_ptr) :: ctx = c_null_ptr
      integer(c_int) :: kk = 0
      integer(c_int) :: nb = 0
      integer(c_int) :: nnmax = 0
      integer(c_int) :: ntype = 0
      integer(c_int) :: nmax = -1
   contains
      procedure :: ensure_context
      procedure :: destroy
      procedure :: set_backend
      procedure :: set_periodic_lattice
      procedure :: set_hamiltonian
      procedure :: set_velocity
      procedure :: upload_bsr
      procedure :: chebyshev_moments
      procedure :: orbital_moments
      procedure :: block_lanczos
      procedure :: scalar_lanczos
      procedure :: stochastic_moments
   end type rsrec_cuda_backend

#ifdef USE_CUDA_PLUGIN
   interface
      function rsrec_cuda_create(kk, nb, nnmax, ntype, nmax, device) bind(C, name='rsrec_cuda_create')
         import :: c_int, c_ptr
         integer(c_int), value :: kk, nb, nnmax, ntype, nmax, device
         type(c_ptr) :: rsrec_cuda_create
      end function rsrec_cuda_create

      subroutine rsrec_cuda_destroy(ctx) bind(C, name='rsrec_cuda_destroy')
         import :: c_ptr
         type(c_ptr), value :: ctx
      end subroutine rsrec_cuda_destroy

      function rsrec_cuda_last_error() bind(C, name='rsrec_cuda_last_error')
         import :: c_ptr
         type(c_ptr) :: rsrec_cuda_last_error
      end function rsrec_cuda_last_error

      function rsrec_cuda_set_backend(ctx, backend) bind(C, name='rsrec_cuda_set_backend')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx
         integer(c_int), value :: backend
         integer(c_int) :: rsrec_cuda_set_backend
      end function rsrec_cuda_set_backend

      function rsrec_cuda_set_periodic_lattice(ctx, pbc, n1, n2, n3, a, crd, nbas) bind(C, name='rsrec_cuda_set_periodic_lattice')
         import :: c_double, c_int, c_ptr
         type(c_ptr), value :: ctx, a, crd
         integer(c_int), value :: pbc, n1, n2, n3, nbas
         integer(c_int) :: rsrec_cuda_set_periodic_lattice
      end function rsrec_cuda_set_periodic_lattice

      function rsrec_cuda_set_hamiltonian(ctx, ee, hall, lsham, nn, iz, eeo, hallo, enim) bind(C, name='rsrec_cuda_set_hamiltonian')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx, ee, hall, lsham, nn, iz
         type(c_ptr), value :: eeo, hallo, enim
         integer(c_int) :: rsrec_cuda_set_hamiltonian
      end function rsrec_cuda_set_hamiltonian

      function rsrec_cuda_set_hamiltonian_additive(ctx, ee_add, hall_add) bind(C, name='rsrec_cuda_set_hamiltonian_additive')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx, ee_add, hall_add
         integer(c_int) :: rsrec_cuda_set_hamiltonian_additive
      end function rsrec_cuda_set_hamiltonian_additive

      function rsrec_cuda_set_velocity(ctx, v_a, v_b, vo_a, vo_b) bind(C, name='rsrec_cuda_set_velocity')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx, v_a, v_b, vo_a, vo_b
         integer(c_int) :: rsrec_cuda_set_velocity
      end function rsrec_cuda_set_velocity

      function rsrec_cuda_orbital_moments(ctx, left, psiref, lld, a, b, mu) bind(C, name='rsrec_cuda_orbital_moments')
         import :: c_double, c_int, c_ptr
         type(c_ptr), value :: ctx, left, psiref, mu
         integer(c_int), value :: lld
         real(c_double), value :: a, b
         integer(c_int) :: rsrec_cuda_orbital_moments
      end function rsrec_cuda_orbital_moments

      function rsrec_cuda_chebyshev_moments(ctx, psi0, lld, a, b, mu_out) bind(C, name='rsrec_cuda_chebyshev_moments')
         import :: c_double, c_int, c_ptr
         type(c_ptr), value :: ctx, psi0, mu_out
         integer(c_int), value :: lld
         real(c_double), value :: a, b
         integer(c_int) :: rsrec_cuda_chebyshev_moments
      end function rsrec_cuda_chebyshev_moments

      function rsrec_cuda_block_lanczos(ctx, psi0, lld, a_b, b2_b, prec) &
         bind(C, name='rsrec_cuda_block_lanczos')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx, psi0, a_b, b2_b
         integer(c_int), value :: lld, prec
         integer(c_int) :: rsrec_cuda_block_lanczos
      end function rsrec_cuda_block_lanczos

      function rsrec_cuda_scalar_lanczos(ctx, site_j, lld, a_out, b2_out) &
         bind(C, name='rsrec_cuda_scalar_lanczos')
         import :: c_double, c_int, c_ptr
         type(c_ptr), value :: ctx, a_out, b2_out
         integer(c_int), value :: site_j, lld
         integer(c_int) :: rsrec_cuda_scalar_lanczos
      end function rsrec_cuda_scalar_lanczos

      function rsrec_cuda_stochastic_moments(ctx, psiref, lld, a, b, mu_nm) bind(C, name='rsrec_cuda_stochastic_moments')
         import :: c_double, c_int, c_ptr
         type(c_ptr), value :: ctx, psiref, mu_nm
         integer(c_int), value :: lld
         real(c_double), value :: a, b
         integer(c_int) :: rsrec_cuda_stochastic_moments
      end function rsrec_cuda_stochastic_moments
   end interface
#endif

contains

   logical function rsrec_cuda_plugin_compiled()
#ifdef USE_CUDA_PLUGIN
      rsrec_cuda_plugin_compiled = .true.
#else
      rsrec_cuda_plugin_compiled = .false.
#endif
   end function rsrec_cuda_plugin_compiled

   integer function decode_gpu_backend(name)
      character(len=*), intent(in) :: name

      select case (trim(name))
      case ('csr')
         decode_gpu_backend = gpu_backend_csr
      case ('bsr')
         decode_gpu_backend = gpu_backend_bsr
      case ('fft')
         decode_gpu_backend = gpu_backend_fft
      case ('conv')
         decode_gpu_backend = gpu_backend_conv
      case default
         decode_gpu_backend = -1
      end select
   end function decode_gpu_backend

   subroutine ensure_context(this, kk, nb, nnmax, ntype, nmax)
      class(rsrec_cuda_backend), intent(inout) :: this
      integer, intent(in) :: kk, nb, nnmax, ntype, nmax
#ifdef USE_CUDA_PLUGIN
      if (c_associated(this%ctx)) then
         if (this%kk == kk .and. this%nb == nb .and. this%nnmax == nnmax .and. &
             this%ntype == ntype .and. this%nmax == nmax) return
         call this%destroy()
      end if

      this%ctx = rsrec_cuda_create(int(kk, c_int), int(nb, c_int), int(nnmax, c_int), &
         int(ntype, c_int), int(nmax, c_int), 0_c_int)
      if (.not. c_associated(this%ctx)) then
         call g_logger%fatal('Failed to create CUDA plugin context: '// &
            trim(last_error_string()), __FILE__, __LINE__)
      end if
      this%kk = int(kk, c_int)
      this%nb = int(nb, c_int)
      this%nnmax = int(nnmax, c_int)
      this%ntype = int(ntype, c_int)
      this%nmax = int(nmax, c_int)
#else
      call g_logger%fatal('gpu_plugin=.true. requested, but this executable '// &
         'was built without ENABLE_CUDA_PLUGIN.', __FILE__, __LINE__)
#endif
   end subroutine ensure_context

   subroutine destroy(this)
      class(rsrec_cuda_backend), intent(inout) :: this
#ifdef USE_CUDA_PLUGIN
      if (c_associated(this%ctx)) call rsrec_cuda_destroy(this%ctx)
#endif
      this%ctx = c_null_ptr
      this%kk = 0
      this%nb = 0
      this%nnmax = 0
      this%ntype = 0
      this%nmax = -1
   end subroutine destroy

   subroutine set_backend(this, backend)
      class(rsrec_cuda_backend), intent(inout) :: this
      integer, intent(in) :: backend
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: status
      status = rsrec_cuda_set_backend(this%ctx, int(backend, c_int))
      call check_status(status, 'rsrec_cuda_set_backend')
#endif
   end subroutine set_backend

   subroutine set_periodic_lattice(this, pbc, n1, n2, n3, a, crd)
      class(rsrec_cuda_backend), intent(inout) :: this
      logical, intent(in) :: pbc
      integer, intent(in) :: n1, n2, n3
      real(rp), target, contiguous, intent(in) :: a(:, :)
      real(rp), target, contiguous, intent(in) :: crd(:, :)
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: status, pbc_i
      pbc_i = 0_c_int
      if (pbc) pbc_i = 1_c_int
      status = rsrec_cuda_set_periodic_lattice(this%ctx, pbc_i, int(n1, c_int), &
         int(n2, c_int), int(n3, c_int), c_loc(a), c_loc(crd), &
         int(size(crd, 2), c_int))
      call check_status(status, 'rsrec_cuda_set_periodic_lattice')
#endif
   end subroutine set_periodic_lattice

   subroutine set_hamiltonian(this, ee, hall, lsham, nn, iz, nmax, eeo, hallo, enim, &
                              ee_add, hall_add)
      class(rsrec_cuda_backend), intent(inout) :: this
      complex(rp), target, contiguous, intent(in) :: ee(:, :, :, :)
      complex(rp), target, contiguous, intent(in), optional :: hall(:, :, :, :)
      complex(rp), target, contiguous, intent(in), optional :: lsham(:, :, :)
      integer, target, contiguous, intent(in) :: nn(:, :)
      integer, target, contiguous, intent(in) :: iz(:)
      integer, intent(in) :: nmax
      ! Orthogonalisation (hoh) operands for the two-sweep apply. When eeo is
      ! present the backend builds H = (h - eeo*h + (enim+lsham) - b)/a.
      complex(rp), target, contiguous, intent(in), optional :: eeo(:, :, :, :)
      complex(rp), target, contiguous, intent(in), optional :: hallo(:, :, :, :)
      complex(rp), target, contiguous, intent(in), optional :: enim(:, :, :)
      complex(rp), target, contiguous, intent(in), optional :: ee_add(:, :, :, :)
      complex(rp), target, contiguous, intent(in), optional :: hall_add(:, :, :, :)
#ifdef USE_CUDA_PLUGIN
      type(c_ptr) :: hall_ptr, lsham_ptr, eeo_ptr, hallo_ptr, enim_ptr
      type(c_ptr) :: ee_add_ptr, hall_add_ptr
      integer(c_int) :: status

      call this%ensure_context(size(nn, 1), size(ee, 1), size(ee, 3), size(ee, 4), nmax)
      hall_ptr = c_null_ptr
      if (present(hall)) hall_ptr = c_loc(hall)
      lsham_ptr = c_null_ptr
      if (present(lsham)) lsham_ptr = c_loc(lsham)
      eeo_ptr = c_null_ptr
      if (present(eeo)) eeo_ptr = c_loc(eeo)
      hallo_ptr = c_null_ptr
      if (present(hallo)) hallo_ptr = c_loc(hallo)
      enim_ptr = c_null_ptr
      if (present(enim)) enim_ptr = c_loc(enim)
      status = rsrec_cuda_set_hamiltonian(this%ctx, c_loc(ee), hall_ptr, &
         lsham_ptr, c_loc(nn), c_loc(iz), eeo_ptr, hallo_ptr, enim_ptr)
      call check_status(status, 'rsrec_cuda_set_hamiltonian')
      if (present(ee_add)) then
         ee_add_ptr = c_loc(ee_add)
         hall_add_ptr = c_null_ptr
         if (present(hall_add)) hall_add_ptr = c_loc(hall_add)
         status = rsrec_cuda_set_hamiltonian_additive(this%ctx, ee_add_ptr, hall_add_ptr)
         call check_status(status, 'rsrec_cuda_set_hamiltonian_additive')
      end if
#else
      call this%ensure_context(size(nn, 1), size(ee, 1), size(ee, 3), size(ee, 4), nmax)
      if (present(eeo) .or. present(hallo) .or. present(enim) .or. &
          present(ee_add) .or. present(hall_add)) continue
#endif
   end subroutine set_hamiltonian

   subroutine set_velocity(this, v_a, v_b, vo_a, vo_b)
      class(rsrec_cuda_backend), intent(inout) :: this
      complex(rp), target, contiguous, intent(in) :: v_a(:, :, :, :)
      complex(rp), target, contiguous, intent(in) :: v_b(:, :, :, :)
      complex(rp), target, contiguous, intent(in), optional :: vo_a(:, :, :, :)
      complex(rp), target, contiguous, intent(in), optional :: vo_b(:, :, :, :)
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: status
      type(c_ptr) :: voa_ptr, vob_ptr
      voa_ptr = c_null_ptr
      vob_ptr = c_null_ptr
      if (present(vo_a)) voa_ptr = c_loc(vo_a)
      if (present(vo_b)) vob_ptr = c_loc(vo_b)
      status = rsrec_cuda_set_velocity(this%ctx, c_loc(v_a), c_loc(v_b), voa_ptr, vob_ptr)
      call check_status(status, 'rsrec_cuda_set_velocity')
#else
      if (present(vo_a) .or. present(vo_b)) continue
#endif
   end subroutine set_velocity

   subroutine orbital_moments(this, left, psiref, lld, a, b, mu)
      class(rsrec_cuda_backend), intent(inout) :: this
      complex(rp), target, contiguous, intent(in) :: left(:, :, :)
      complex(rp), target, contiguous, intent(in) :: psiref(:, :, :)
      integer, intent(in) :: lld
      real(rp), intent(in) :: a, b
      complex(rp), target, contiguous, intent(out) :: mu(:, :, :)
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: status
      status = rsrec_cuda_orbital_moments(this%ctx, c_loc(left), c_loc(psiref), &
         int(lld, c_int), real(a, c_double), real(b, c_double), c_loc(mu))
      call check_status(status, 'rsrec_cuda_orbital_moments')
#endif
   end subroutine orbital_moments

   subroutine upload_bsr(this, sparse_obj)
      class(rsrec_cuda_backend), intent(inout) :: this
      type(sparse), intent(inout) :: sparse_obj
      logical :: is_valid

      call sparse_obj%export_gpu_cusparse(is_valid)
      if (.not. is_valid) then
         call g_logger%fatal('Failed to export BSR data for CUDA backend.', __FILE__, __LINE__)
      end if
   end subroutine upload_bsr

   subroutine chebyshev_moments(this, psi0, lld, a, b, mu_out)
      class(rsrec_cuda_backend), intent(inout) :: this
      complex(rp), target, contiguous, intent(in) :: psi0(:, :, :)
      integer, intent(in) :: lld
      real(rp), intent(in) :: a, b
      complex(rp), target, contiguous, intent(out) :: mu_out(:, :, :)
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: status
      status = rsrec_cuda_chebyshev_moments(this%ctx, c_loc(psi0), int(lld, c_int), &
         real(a, c_double), real(b, c_double), c_loc(mu_out))
      call check_status(status, 'rsrec_cuda_chebyshev_moments')
#endif
   end subroutine chebyshev_moments

   subroutine block_lanczos(this, psi0, lld, a_b, b2_b, prec)
      class(rsrec_cuda_backend), intent(inout) :: this
      complex(rp), target, contiguous, intent(in) :: psi0(:, :, :)
      integer, intent(in) :: lld
      complex(rp), target, contiguous, intent(out) :: a_b(:, :, :)
      complex(rp), target, contiguous, intent(out) :: b2_b(:, :, :)
      !> Working precision for the hoh block engine: 0 = fp32 (mixed, fp64
      !> B-sqrt), 1 = fp64. Ignored by the ee-only (non-hoh) path.
      integer, intent(in), optional :: prec
      integer :: prec_
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: status
#endif
      prec_ = 0
      if (present(prec)) prec_ = prec
#ifdef USE_CUDA_PLUGIN
      status = rsrec_cuda_block_lanczos(this%ctx, c_loc(psi0), int(lld, c_int), &
         c_loc(a_b), c_loc(b2_b), int(prec_, c_int))
      call check_status(status, 'rsrec_cuda_block_lanczos')
#endif
   end subroutine block_lanczos

   subroutine scalar_lanczos(this, site_j, lld, a_out, b2_out)
      class(rsrec_cuda_backend), intent(inout) :: this
      integer, intent(in) :: site_j, lld
      real(rp), target, contiguous, intent(out) :: a_out(:, :)
      real(rp), target, contiguous, intent(out) :: b2_out(:, :)
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: status
      status = rsrec_cuda_scalar_lanczos(this%ctx, int(site_j, c_int), &
         int(lld, c_int), c_loc(a_out), c_loc(b2_out))
      call check_status(status, 'rsrec_cuda_scalar_lanczos')
#endif
   end subroutine scalar_lanczos

   subroutine stochastic_moments(this, psiref, lld, a, b, mu_nm)
      class(rsrec_cuda_backend), intent(inout) :: this
      complex(rp), target, contiguous, intent(in) :: psiref(:, :, :)
      integer, intent(in) :: lld
      real(rp), intent(in) :: a, b
      complex(rp), target, contiguous, intent(out) :: mu_nm(:, :, :, :)
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: status
      status = rsrec_cuda_stochastic_moments(this%ctx, c_loc(psiref), int(lld, c_int), &
         real(a, c_double), real(b, c_double), c_loc(mu_nm))
      call check_status(status, 'rsrec_cuda_stochastic_moments')
#endif
   end subroutine stochastic_moments

#ifdef USE_CUDA_PLUGIN
   subroutine check_status(status, where)
      integer(c_int), intent(in) :: status
      character(len=*), intent(in) :: where

      if (status /= 0_c_int) then
         call g_logger%fatal(trim(where)//' failed: '//trim(last_error_string()), __FILE__, __LINE__)
      end if
   end subroutine check_status

   function last_error_string() result(msg)
      character(len=512) :: msg
      type(c_ptr) :: err_ptr
      character(kind=c_char), pointer :: chars(:)
      integer :: i, n

      msg = 'unknown error'
      err_ptr = rsrec_cuda_last_error()
      if (.not. c_associated(err_ptr)) return
      call c_f_pointer(err_ptr, chars, [512])
      n = 0
      do i = 1, size(chars)
         if (chars(i) == c_null_char) exit
         n = n + 1
      end do
      if (n > 0) then
         msg = ''
         do i = 1, n
            msg(i:i) = chars(i)
         end do
      end if
   end function last_error_string
#endif

end module rsrec_cuda_plugin_mod
