module rsrec_cuda_plugin_mod

   use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, &
      c_double_complex, c_int, c_long_long, c_null_char, c_null_ptr, c_ptr, c_f_pointer, c_loc
   use precision_mod, only: rp
   use logger_mod, only: g_logger
   use mpi_mod, only: g_parallel_context, get_cuda_device_override
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
   public :: get_gpu_context

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
      procedure :: set_precision
      procedure :: chebyshev_dos
      procedure :: chebyshev_gf_eta
      procedure :: block_dos
      procedure :: block_gf_eta
      procedure :: chebyshev_moments
      procedure :: orbital_moments
      procedure :: block_lanczos
      procedure :: scalar_lanczos
      procedure :: stochastic_moments
      procedure :: stochastic_profile
   end type rsrec_cuda_backend

   type(rsrec_cuda_backend), target, save :: shared_gpu_context

#ifdef USE_CUDA_PLUGIN
   interface
      function rsrec_cuda_device_count(count) bind(C, name='rsrec_cuda_device_count')
         import :: c_int
         integer(c_int), intent(out) :: count
         integer(c_int) :: rsrec_cuda_device_count
      end function rsrec_cuda_device_count

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

      function rsrec_cuda_stochastic_profile(ctx, h2d_seconds, cheb_seconds, d2h_seconds, h2d_bytes, d2h_bytes) &
         bind(C, name='rsrec_cuda_stochastic_profile')
         import :: c_double, c_int, c_long_long, c_ptr
         type(c_ptr), value :: ctx
         real(c_double), intent(out) :: h2d_seconds, cheb_seconds, d2h_seconds
         integer(c_long_long), intent(out) :: h2d_bytes, d2h_bytes
         integer(c_int) :: rsrec_cuda_stochastic_profile
      end function rsrec_cuda_stochastic_profile

      function rsrec_cuda_set_precision(ctx, prec) bind(C, name='rsrec_cuda_set_precision')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx
         integer(c_int), value :: prec
         integer(c_int) :: rsrec_cuda_set_precision
      end function rsrec_cuda_set_precision

      function rsrec_cuda_chebyshev_dos(ctx, mu, n_mom, natoms, ene, nv, a, b, g0) &
         bind(C, name='rsrec_cuda_chebyshev_dos')
         import :: c_double, c_int, c_ptr
         type(c_ptr), value :: ctx, mu, ene, g0
         integer(c_int), value :: n_mom, natoms, nv
         real(c_double), value :: a, b
         integer(c_int) :: rsrec_cuda_chebyshev_dos
      end function rsrec_cuda_chebyshev_dos

      function rsrec_cuda_chebyshev_gf_eta(ctx, mu, n_mom, natoms, f, n_eta, g0) &
         bind(C, name='rsrec_cuda_chebyshev_gf_eta')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx, mu, f, g0
         integer(c_int), value :: n_mom, natoms, n_eta
         integer(c_int) :: rsrec_cuda_chebyshev_gf_eta
      end function rsrec_cuda_chebyshev_gf_eta

      function rsrec_cuda_block_dos(ctx, a_b, b2_b, a_inf, b_inf, ene, nv, &
                                    eta_re, eta_im, natoms, lld, sym, g0) &
         bind(C, name='rsrec_cuda_block_dos')
         import :: c_double, c_int, c_ptr
         type(c_ptr), value :: ctx, a_b, b2_b, a_inf, b_inf, ene, g0
         integer(c_int), value :: nv, natoms, lld, sym
         real(c_double), value :: eta_re, eta_im
         integer(c_int) :: rsrec_cuda_block_dos
      end function rsrec_cuda_block_dos

      function rsrec_cuda_block_gf_eta(ctx, a_b, b2_b, a_inf, b_inf, ef, &
                                       eta_re, eta_im, n_eta, natoms, lld, sym, g0) &
         bind(C, name='rsrec_cuda_block_gf_eta')
         import :: c_double, c_int, c_ptr
         type(c_ptr), value :: ctx, a_b, b2_b, a_inf, b_inf, eta_re, eta_im, g0
         integer(c_int), value :: n_eta, natoms, lld, sym
         real(c_double), value :: ef
         integer(c_int) :: rsrec_cuda_block_gf_eta
      end function rsrec_cuda_block_gf_eta
   end interface
#endif

contains

   function get_gpu_context(kk, nb, nnmax, ntype, nmax) result(ctx)
      integer, intent(in) :: kk, nb, nnmax, ntype, nmax
      type(rsrec_cuda_backend), pointer :: ctx

      call shared_gpu_context%ensure_context(kk, nb, nnmax, ntype, nmax)
      ctx => shared_gpu_context
   end function get_gpu_context

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
      integer(c_int) :: device_count, status
      integer :: selected_device, override_device
      logical :: device_valid, override_configured, override_valid
      character(len=256) :: mapping_message

      if (c_associated(this%ctx)) then
         if (this%kk == kk .and. this%nb == nb .and. this%nnmax == nnmax .and. &
             this%ntype == ntype .and. this%nmax == nmax) return
         call this%destroy()
      end if

      device_count = 0_c_int
      status = rsrec_cuda_device_count(device_count)
      if (status /= 0_c_int .or. device_count <= 0_c_int) then
         call g_logger%fatal('Failed to discover a usable CUDA device for the recursion plugin.', __FILE__, __LINE__)
      end if
      call get_cuda_device_override(override_device, override_configured, override_valid)
      if (override_configured .and. .not. override_valid) then
         call g_logger%fatal('RSLMTO_CUDA_DEVICE must be a non-negative integer.', __FILE__, __LINE__)
      end if
      if (override_configured) then
         call g_parallel_context%device_index(int(device_count), selected_device, device_valid, override_device)
      else
         call g_parallel_context%device_index(int(device_count), selected_device, device_valid)
      end if
      if (.not. device_valid) then
         call g_logger%fatal('CUDA recursion plugin refuses unconfigured MPI local-rank/device oversubscription.', &
            __FILE__, __LINE__)
      end if
      write(mapping_message, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'CUDA_DEVICE_MAPPING world_rank=', &
         g_parallel_context%rank, ' local_rank=', g_parallel_context%local_rank, &
         ' local_size=', g_parallel_context%local_size, ' visible_devices=', int(device_count), &
         ' selected_device=', selected_device
      call g_logger%info(trim(mapping_message), __FILE__, __LINE__)

      this%ctx = rsrec_cuda_create(int(kk, c_int), int(nb, c_int), int(nnmax, c_int), &
         int(ntype, c_int), int(nmax, c_int), int(selected_device, c_int))
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

   subroutine set_hamiltonian(this, ee, hall, lsham, nn, iz, nmax, eeo, hallo, enim)
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
#ifdef USE_CUDA_PLUGIN
      type(c_ptr) :: hall_ptr, lsham_ptr, eeo_ptr, hallo_ptr, enim_ptr
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
#else
      call this%ensure_context(size(nn, 1), size(ee, 1), size(ee, 3), size(ee, 4), nmax)
      if (present(eeo) .or. present(hallo) .or. present(enim)) continue
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

   subroutine set_precision(this, prec)
      class(rsrec_cuda_backend), intent(inout) :: this
      integer, intent(in) :: prec
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: status
      status = rsrec_cuda_set_precision(this%ctx, int(prec, c_int))
      call check_status(status, 'rsrec_cuda_set_precision')
#endif
   end subroutine set_precision

   subroutine chebyshev_dos(this, mu_n, ene, a, b, g0)
      class(rsrec_cuda_backend), intent(inout) :: this
      complex(rp), target, contiguous, intent(in) :: mu_n(:, :, :, :)
      real(rp), target, contiguous, intent(in) :: ene(:)
      real(rp), intent(in) :: a, b
      complex(rp), target, contiguous, intent(inout) :: g0(:, :, :, :)
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: status
      status = rsrec_cuda_chebyshev_dos(this%ctx, c_loc(mu_n), int(size(mu_n, 3), c_int), &
         int(size(mu_n, 4), c_int), c_loc(ene), int(size(ene), c_int), &
         real(a, c_double), real(b, c_double), c_loc(g0))
      call check_status(status, 'rsrec_cuda_chebyshev_dos')
#endif
   end subroutine chebyshev_dos

   subroutine chebyshev_gf_eta(this, mu_n, ef, etas, a, b, g0)
      class(rsrec_cuda_backend), intent(inout) :: this
      complex(rp), target, contiguous, intent(in) :: mu_n(:, :, :, :)
      real(rp), intent(in) :: ef, a, b
      complex(rp), intent(in) :: etas(:)
      complex(rp), target, contiguous, intent(inout) :: g0(:, :, :, :)
#ifdef USE_CUDA_PLUGIN
      complex(rp), allocatable, target :: f(:, :)
      complex(rp) :: ze, z, th, pref
      real(rp) :: thl, cot, gj, cc
      integer(c_int) :: n_mom, natoms, n_eta, status
      integer :: i, k

      n_mom = int(size(mu_n, 3), c_int)
      natoms = int(size(mu_n, 4), c_int)
      n_eta = int(size(etas), c_int)

      allocate (f(n_mom, n_eta))
      cot = cos(acos(-1.0_rp)/(n_mom + 1.0_rp))/sin(acos(-1.0_rp)/(n_mom + 1.0_rp))
      do k = 1, n_eta
         ze = cmplx(ef, 0.0_rp, kind=rp) + etas(k)
         z = (ze - b)/a
         th = acos(z)
         pref = 1.0_rp/sqrt(a*a - (ze - b)**2)
         do i = 1, n_mom
            thl = acos(-1.0_rp)*real(i - 1, rp)/(n_mom + 1.0_rp)
            gj = ((n_mom - (i - 1) + 1)*cos(thl) + sin(thl)*cot)/(n_mom + 1.0_rp)
            cc = merge(1.0_rp, 2.0_rp, i == 1)
            f(i, k) = gj*cc*pref*(-cmplx(0.0_rp, 1.0_rp, kind=rp))* &
                      exp(-cmplx(0.0_rp, 1.0_rp, kind=rp)*real(i - 1, rp)*th)
         end do
      end do

      status = rsrec_cuda_chebyshev_gf_eta(this%ctx, c_loc(mu_n), n_mom, natoms, &
         c_loc(f), n_eta, c_loc(g0))
      call check_status(status, 'rsrec_cuda_chebyshev_gf_eta')
      deallocate (f)
#endif
   end subroutine chebyshev_gf_eta

   subroutine block_dos(this, a_b, b2_b, a_inf, b_inf, ene, eta, sym, g0)
      class(rsrec_cuda_backend), intent(inout) :: this
      complex(rp), target, contiguous, intent(in) :: a_b(:, :, :, :), b2_b(:, :, :, :)
      real(rp), target, contiguous, intent(in) :: a_inf(:, :), b_inf(:, :)
      real(rp), target, contiguous, intent(in) :: ene(:)
      complex(rp), intent(in) :: eta
      logical, intent(in) :: sym
      complex(rp), target, contiguous, intent(inout) :: g0(:, :, :, :)
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: status
      status = rsrec_cuda_block_dos(this%ctx, c_loc(a_b), c_loc(b2_b), c_loc(a_inf), &
         c_loc(b_inf), c_loc(ene), int(size(ene), c_int), real(real(eta), c_double), &
         real(aimag(eta), c_double), int(size(a_b, 4), c_int), int(size(a_b, 3), c_int), &
         merge(1_c_int, 0_c_int, sym), c_loc(g0))
      call check_status(status, 'rsrec_cuda_block_dos')
#endif
   end subroutine block_dos

   subroutine block_gf_eta(this, a_b, b2_b, a_inf, b_inf, ef, etas, sym, g0)
      class(rsrec_cuda_backend), intent(inout) :: this
      complex(rp), target, contiguous, intent(in) :: a_b(:, :, :, :), b2_b(:, :, :, :)
      real(rp), target, contiguous, intent(in) :: a_inf(:, :), b_inf(:, :)
      real(rp), intent(in) :: ef
      complex(rp), intent(in) :: etas(:)
      logical, intent(in) :: sym
      complex(rp), target, contiguous, intent(inout) :: g0(:, :, :, :)
#ifdef USE_CUDA_PLUGIN
      real(c_double), allocatable, target :: etr(:), eti(:)
      integer(c_int) :: status
      integer :: k

      allocate (etr(size(etas)), eti(size(etas)))
      do k = 1, size(etas)
         etr(k) = real(real(etas(k)), c_double)
         eti(k) = real(aimag(etas(k)), c_double)
      end do
      status = rsrec_cuda_block_gf_eta(this%ctx, c_loc(a_b), c_loc(b2_b), c_loc(a_inf), &
         c_loc(b_inf), real(ef, c_double), c_loc(etr), c_loc(eti), int(size(etas), c_int), &
         int(size(a_b, 4), c_int), int(size(a_b, 3), c_int), merge(1_c_int, 0_c_int, sym), &
         c_loc(g0))
      call check_status(status, 'rsrec_cuda_block_gf_eta')
      deallocate (etr, eti)
#endif
   end subroutine block_gf_eta

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

   subroutine stochastic_profile(this, h2d_seconds, cheb_seconds, d2h_seconds, h2d_bytes, d2h_bytes)
      class(rsrec_cuda_backend), intent(inout) :: this
      real(rp), intent(out) :: h2d_seconds, cheb_seconds, d2h_seconds
      integer(c_long_long), intent(out) :: h2d_bytes, d2h_bytes
#ifdef USE_CUDA_PLUGIN
      integer(c_int) :: status
      status = rsrec_cuda_stochastic_profile(this%ctx, h2d_seconds, cheb_seconds, d2h_seconds, &
         h2d_bytes, d2h_bytes)
      call check_status(status, 'rsrec_cuda_stochastic_profile')
#else
      h2d_seconds = 0.0_rp
      cheb_seconds = 0.0_rp
      d2h_seconds = 0.0_rp
      h2d_bytes = 0
      d2h_bytes = 0
#endif
   end subroutine stochastic_profile

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
