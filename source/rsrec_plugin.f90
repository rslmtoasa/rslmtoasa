module rsrec_plugin_mod

   use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, &
      c_int, c_null_char, c_null_ptr, c_ptr, c_double_complex, c_f_pointer, c_loc
   use precision_mod, only: rp
   use logger_mod, only: g_logger
   implicit none

   private

   public :: rsrec_backend
   public :: rsrec_plugin_compiled

   type :: rsrec_backend
      type(c_ptr) :: ctx = c_null_ptr
      integer(c_int) :: kk = 0
      integer(c_int) :: nb = 0
      integer(c_int) :: nnmax = 0
      integer(c_int) :: ntype = 0
      integer(c_int) :: nmax = -1
   contains
      procedure :: ensure_context
      procedure :: destroy
      procedure :: set_hamiltonian
      procedure :: set_velocity
      procedure :: scalar_lanczos
      procedure :: block_lanczos
      procedure :: chebyshev_moments
      procedure :: chebyshev_dos
      procedure :: stochastic_moments
   end type rsrec_backend

#ifdef USE_CPP_PLUGIN
   interface
      function rsrec_create(kk, nb, nnmax, ntype, nmax, device) bind(C, name='rsrec_create')
         import :: c_int, c_ptr
         integer(c_int), value :: kk, nb, nnmax, ntype, nmax, device
         type(c_ptr) :: rsrec_create
      end function rsrec_create

      subroutine rsrec_destroy(ctx) bind(C, name='rsrec_destroy')
         import :: c_ptr
         type(c_ptr), value :: ctx
      end subroutine rsrec_destroy

      function rsrec_last_error() bind(C, name='rsrec_last_error')
         import :: c_ptr
         type(c_ptr) :: rsrec_last_error
      end function rsrec_last_error

      function rsrec_set_hamiltonian(ctx, ee, hall, lsham, nn, iz) bind(C, name='rsrec_set_hamiltonian')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx, ee, hall, lsham, nn, iz
         integer(c_int) :: rsrec_set_hamiltonian
      end function rsrec_set_hamiltonian

      function rsrec_set_velocity(ctx, v_a, v_b) bind(C, name='rsrec_set_velocity')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx, v_a, v_b
         integer(c_int) :: rsrec_set_velocity
      end function rsrec_set_velocity

      function rsrec_chebyshev_moments(ctx, psi0, lld, a, b, mu_out) bind(C, name='rsrec_chebyshev_moments')
         import :: c_double, c_int, c_ptr
         type(c_ptr), value :: ctx, psi0, mu_out
         integer(c_int), value :: lld
         real(c_double), value :: a, b
         integer(c_int) :: rsrec_chebyshev_moments
      end function rsrec_chebyshev_moments

      function rsrec_block_lanczos(ctx, psi0, lld, a_b, b2_b) bind(C, name='rsrec_block_lanczos')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx, psi0, a_b, b2_b
         integer(c_int), value :: lld
         integer(c_int) :: rsrec_block_lanczos
      end function rsrec_block_lanczos

      function rsrec_scalar_lanczos(ctx, site_j, lld, a_out, b2_out) bind(C, name='rsrec_scalar_lanczos')
         import :: c_int, c_ptr
         type(c_ptr), value :: ctx, a_out, b2_out
         integer(c_int), value :: site_j, lld
         integer(c_int) :: rsrec_scalar_lanczos
      end function rsrec_scalar_lanczos

      function rsrec_chebyshev_dos(ctx, mu, n_moments, natoms, ene, nv, a, b, g0) &
         bind(C, name='rsrec_chebyshev_dos')
         import :: c_double, c_int, c_ptr
         type(c_ptr), value :: ctx, mu, ene, g0
         integer(c_int), value :: n_moments, natoms, nv
         real(c_double), value :: a, b
         integer(c_int) :: rsrec_chebyshev_dos
      end function rsrec_chebyshev_dos

      function rsrec_stochastic_moments(ctx, psiref, lld, a, b, mu_nm) bind(C, name='rsrec_stochastic_moments')
         import :: c_double, c_int, c_ptr
         type(c_ptr), value :: ctx, psiref, mu_nm
         integer(c_int), value :: lld
         real(c_double), value :: a, b
         integer(c_int) :: rsrec_stochastic_moments
      end function rsrec_stochastic_moments
   end interface
#endif

contains

   logical function rsrec_plugin_compiled()
#ifdef USE_CPP_PLUGIN
      rsrec_plugin_compiled = .true.
#else
      rsrec_plugin_compiled = .false.
#endif
   end function rsrec_plugin_compiled

   subroutine ensure_context(this, kk, nb, nnmax, ntype, nmax)
      class(rsrec_backend), intent(inout) :: this
      integer, intent(in) :: kk, nb, nnmax, ntype, nmax
#ifdef USE_CPP_PLUGIN
      if (c_associated(this%ctx)) then
         if (this%kk == kk .and. this%nb == nb .and. this%nnmax == nnmax .and. this%ntype == ntype .and. this%nmax == nmax) return
         call this%destroy()
      end if

      this%ctx = rsrec_create(int(kk, c_int), int(nb, c_int), int(nnmax, c_int), &
         int(ntype, c_int), int(nmax, c_int), 0_c_int)
      if (.not. c_associated(this%ctx)) then
         call g_logger%fatal('Failed to create rsrec C++ plugin context: '// &
            trim(last_error_string()), __FILE__, __LINE__)
      end if

      this%kk = int(kk, c_int)
      this%nb = int(nb, c_int)
      this%nnmax = int(nnmax, c_int)
      this%ntype = int(ntype, c_int)
      this%nmax = int(nmax, c_int)
#else
      call g_logger%fatal('cpp_plugin=.true. requested, but this executable '// &
         'was built without ENABLE_CPP_PLUGIN.', __FILE__, __LINE__)
#endif
   end subroutine ensure_context

   subroutine destroy(this)
      class(rsrec_backend), intent(inout) :: this
#ifdef USE_CPP_PLUGIN
      if (c_associated(this%ctx)) call rsrec_destroy(this%ctx)
#endif
      this%ctx = c_null_ptr
      this%kk = 0
      this%nb = 0
      this%nnmax = 0
      this%ntype = 0
      this%nmax = -1
   end subroutine destroy

   subroutine set_hamiltonian(this, ee, hall, lsham, nn, iz, nmax)
      class(rsrec_backend), intent(inout) :: this
      complex(rp), target, contiguous, intent(in) :: ee(:, :, :, :)
      complex(rp), target, contiguous, intent(in), optional :: hall(:, :, :, :)
      complex(rp), target, contiguous, intent(in), optional :: lsham(:, :, :)
      integer, target, contiguous, intent(in) :: nn(:, :)
      integer, target, contiguous, intent(in) :: iz(:)
      integer, intent(in) :: nmax
#ifdef USE_CPP_PLUGIN
      type(c_ptr) :: hall_ptr, lsham_ptr
      integer(c_int) :: status

      call this%ensure_context(size(nn, 1), size(ee, 1), size(ee, 3), size(ee, 4), nmax)

      hall_ptr = c_null_ptr
      if (present(hall)) hall_ptr = c_loc(hall)
      lsham_ptr = c_null_ptr
      if (present(lsham)) lsham_ptr = c_loc(lsham)

      status = rsrec_set_hamiltonian(this%ctx, c_loc(ee), hall_ptr, lsham_ptr, &
         c_loc(nn), c_loc(iz))
      call check_status(status, 'rsrec_set_hamiltonian')
#else
      call this%ensure_context(size(nn, 1), size(ee, 1), size(ee, 3), size(ee, 4), nmax)
#endif
   end subroutine set_hamiltonian

   subroutine set_velocity(this, v_a, v_b)
      class(rsrec_backend), intent(inout) :: this
      complex(rp), target, contiguous, intent(in) :: v_a(:, :, :, :)
      complex(rp), target, contiguous, intent(in) :: v_b(:, :, :, :)
#ifdef USE_CPP_PLUGIN
      integer(c_int) :: status

      status = rsrec_set_velocity(this%ctx, c_loc(v_a), c_loc(v_b))
      call check_status(status, 'rsrec_set_velocity')
#endif
   end subroutine set_velocity

   subroutine scalar_lanczos(this, site_j, lld, a_out, b2_out)
      class(rsrec_backend), intent(inout) :: this
      integer, intent(in) :: site_j, lld
      real(rp), target, contiguous, intent(out) :: a_out(:, :)
      real(rp), target, contiguous, intent(out) :: b2_out(:, :)
#ifdef USE_CPP_PLUGIN
      integer(c_int) :: status

      status = rsrec_scalar_lanczos(this%ctx, int(site_j, c_int), int(lld, c_int), &
         c_loc(a_out), c_loc(b2_out))
      call check_status(status, 'rsrec_scalar_lanczos')
#endif
   end subroutine scalar_lanczos

   subroutine block_lanczos(this, psi0, lld, a_b, b2_b)
      class(rsrec_backend), intent(inout) :: this
      complex(rp), target, contiguous, intent(in) :: psi0(:, :, :)
      integer, intent(in) :: lld
      complex(rp), target, contiguous, intent(out) :: a_b(:, :, :)
      complex(rp), target, contiguous, intent(out) :: b2_b(:, :, :)
#ifdef USE_CPP_PLUGIN
      integer(c_int) :: status

      status = rsrec_block_lanczos(this%ctx, c_loc(psi0), int(lld, c_int), c_loc(a_b), c_loc(b2_b))
      call check_status(status, 'rsrec_block_lanczos')
#endif
   end subroutine block_lanczos

   subroutine chebyshev_moments(this, psi0, lld, a, b, mu_out)
      class(rsrec_backend), intent(inout) :: this
      complex(rp), target, contiguous, intent(in) :: psi0(:, :, :)
      integer, intent(in) :: lld
      real(rp), intent(in) :: a, b
      complex(rp), target, contiguous, intent(out) :: mu_out(:, :, :)
#ifdef USE_CPP_PLUGIN
      integer(c_int) :: status

      status = rsrec_chebyshev_moments(this%ctx, c_loc(psi0), int(lld, c_int), &
         real(a, c_double), real(b, c_double), c_loc(mu_out))
      call check_status(status, 'rsrec_chebyshev_moments')
#endif
   end subroutine chebyshev_moments

   subroutine chebyshev_dos(this, mu_n, ene, a, b, g0)
      class(rsrec_backend), intent(inout) :: this
      complex(rp), target, contiguous, intent(in) :: mu_n(:, :, :, :)
      real(rp), target, contiguous, intent(in) :: ene(:)
      real(rp), intent(in) :: a, b
      complex(rp), target, contiguous, intent(inout) :: g0(:, :, :, :)
#ifdef USE_CPP_PLUGIN
      integer(c_int) :: status, n_mom, natoms, nv

      n_mom = int(size(mu_n, 3), c_int)
      natoms = int(size(mu_n, 4), c_int)
      nv = int(size(ene), c_int)

      status = rsrec_chebyshev_dos(this%ctx, c_loc(mu_n), n_mom, natoms, &
         c_loc(ene), nv, real(a, c_double), real(b, c_double), c_loc(g0))
      call check_status(status, 'rsrec_chebyshev_dos')
#endif
   end subroutine chebyshev_dos

   subroutine stochastic_moments(this, psiref, lld, a, b, mu_nm)
      class(rsrec_backend), intent(inout) :: this
      complex(rp), target, contiguous, intent(in) :: psiref(:, :, :)
      integer, intent(in) :: lld
      real(rp), intent(in) :: a, b
      complex(rp), target, contiguous, intent(out) :: mu_nm(:, :, :, :)
#ifdef USE_CPP_PLUGIN
      integer(c_int) :: status

      status = rsrec_stochastic_moments(this%ctx, c_loc(psiref), int(lld, c_int), &
         real(a, c_double), real(b, c_double), c_loc(mu_nm))
      call check_status(status, 'rsrec_stochastic_moments')
#endif
   end subroutine stochastic_moments

#ifdef USE_CPP_PLUGIN
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
      err_ptr = rsrec_last_error()
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

end module rsrec_plugin_mod
