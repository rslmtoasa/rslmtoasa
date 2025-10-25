!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Spectrum Bounds
!
! DESCRIPTION:
!> Module to compute tight spectrum bounds for Tight-Binding Hamiltonians
!> Useful for KPM/Chebyshev polynomial expansion
!> Implements Gershgorin circle theorem and Sturm/exact diagonalization
!> Also provides an OO wrapper `bounds` type with constructor.
!> Configuration is typically done via the parent object's namelist
!> (e.g., hamiltonian reads bounds_algorithm and bounds_scaling from /hamiltonian/ namelist).
!------------------------------------------------------------------------------

module spectrum_bounds_mod

   use precision_mod, only: rp
   use logger_mod, only: g_logger
   implicit none

   real(rp), parameter :: default_scaling = 1.05_rp
   private

   public :: bounds, compute_spectrum_bounds, bounds_constructor

   !> Single type to store spectrum bounds and user options
   type, public :: bounds
      ! Numerical results (from Gershgorin/Sturm/exact)
      real(rp) :: e_min_sturm = 0.0_rp
      real(rp) :: e_max_sturm = 0.0_rp
      real(rp) :: e_min_gershgorin = 0.0_rp
      real(rp) :: e_max_gershgorin = 0.0_rp
      real(rp) :: e_min = 0.0_rp
      real(rp) :: e_max = 0.0_rp
      real(rp) :: sturm_width = 0.0_rp
      real(rp) :: gershgorin_width = 0.0_rp
      logical :: sturm_available = .false.
      logical :: use_sturm = .false.

      ! User-configurable options (set via parent object's namelist)
      character(len=16) :: algorithm = 'none'
      real(rp) :: scaling = default_scaling

   contains
      procedure :: restore_to_default
      ! read_from_namelist is DEPRECATED - use parent object's namelist instead
      ! (e.g., hamiltonian reads bounds_algorithm/bounds_scaling from /hamiltonian/)
      procedure :: read_from_namelist
   end type bounds

contains

   !---------------------------------------------------------------------------
   ! Compute spectrum bounds (Gershgorin + optional Sturm/Exact)
   subroutine compute_spectrum_bounds(H, res, method, verbose)
      complex(rp), dimension(:, :), intent(in) :: H
   type(bounds), intent(out) :: res
      character(len=*), intent(in), optional :: method
      logical, intent(in), optional :: verbose

      character(len=32) :: method_to_use
      logical :: verb, try_sturm
      integer :: n

      n = size(H, 1)
      if (size(H, 2) /= n) then
         call g_logger%error('spectrum_bounds: Hamiltonian must be square', __FILE__, __LINE__)
         return
      end if

      verb = .false.
      if (present(verbose)) verb = verbose

      method_to_use = 'both'
      if (present(method)) method_to_use = trim(adjustl(method))

   res%sturm_available = .false.
   res%use_sturm = .false.

      ! Always compute Gershgorin (fast, O(n^2))
   call compute_gershgorin_bounds(H, res%e_min_gershgorin, res%e_max_gershgorin)
   res%gershgorin_width = res%e_max_gershgorin - res%e_min_gershgorin

      ! Try Sturm if requested (diagonalization-based fallback)
      try_sturm = (method_to_use == 'sturm' .or. method_to_use == 'both')
      if (try_sturm) then
         call compute_sturm_bounds(H, res%e_min_sturm, res%e_max_sturm, res%sturm_available)
         if (res%sturm_available) then
            res%sturm_width = res%e_max_sturm - res%e_min_sturm
            res%use_sturm = .true.
         end if
      end if

      ! Select best bounds
      if (res%sturm_available) then
         res%e_min = res%e_min_sturm
         res%e_max = res%e_max_sturm
      else
         res%e_min = res%e_min_gershgorin
         res%e_max = res%e_max_gershgorin
      end if

      if (verb) then
         call print_bounds_info(res)
      end if

   end subroutine compute_spectrum_bounds


   !---------------------------------------------------------------------------
   ! Compute Gershgorin bounds
   subroutine compute_gershgorin_bounds(H, e_min, e_max)
      complex(rp), dimension(:, :), intent(in) :: H
      real(rp), intent(out) :: e_min, e_max

      integer :: i, j, n
      real(rp) :: diag, radius
      real(rp), parameter :: eps = 1.0d-15

      n = size(H, 1)
      e_max = -huge(1.0_rp)
      e_min = huge(1.0_rp)

      do i = 1, n
         diag = real(H(i, i))
         radius = 0.0_rp
         do j = 1, n
            if (i /= j) then
               radius = radius + abs(H(i, j))
            end if
         end do
         e_max = max(e_max, diag + radius)
         e_min = min(e_min, diag - radius)
      end do

      ! Small floating-point safety
      e_max = e_max * (1.0_rp + eps)
      e_min = e_min * (1.0_rp - eps)

   end subroutine compute_gershgorin_bounds


   !---------------------------------------------------------------------------
   ! Compute bounds via diagonalization (using LAPACK zheev)
   subroutine compute_sturm_bounds(H, e_min, e_max, success)
      complex(rp), dimension(:, :), intent(in) :: H
      real(rp), intent(out) :: e_min, e_max
      logical, intent(out) :: success

      integer :: n, info, lwork, alloc_stat
      complex(rp), allocatable :: work(:)
      real(rp), allocatable :: rwork(:)
      real(rp), allocatable :: evals(:)
      complex(rp), allocatable :: H_work(:,:)

      n = size(H, 1)
      success = .false.

      allocate(H_work(n, n), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call g_logger%warning('compute_sturm_bounds: failed to allocate H_work', __FILE__, __LINE__)
         return
      end if
      H_work = H

      allocate(evals(n), rwork(3*n - 2), stat=alloc_stat)
      if (alloc_stat /= 0) then
         if (allocated(H_work)) deallocate(H_work)
         call g_logger%warning('compute_sturm_bounds: failed to allocate eval arrays', __FILE__, __LINE__)
         return
      end if

      ! Query optimal complex workspace
      allocate(work(1), stat=alloc_stat)
      if (alloc_stat /= 0) then
         deallocate(H_work, evals, rwork)
         return
      end if
      call zheev('N', 'U', n, H_work, n, evals, work, -1, rwork, info)
      if (info /= 0) then
         if (allocated(work)) deallocate(work)
         if (allocated(rwork)) deallocate(rwork)
         if (allocated(evals)) deallocate(evals)
         if (allocated(H_work)) deallocate(H_work)
         return
      end if

      lwork = max(1, int(real(work(1))))
      deallocate(work)
      allocate(work(lwork))

      call zheev('N', 'U', n, H_work, n, evals, work, lwork, rwork, info)
      if (info == 0) then
         e_min = minval(evals)
         e_max = maxval(evals)
         success = .true.
      end if

      if (allocated(work)) deallocate(work)
      if (allocated(rwork)) deallocate(rwork)
      if (allocated(evals)) deallocate(evals)
      if (allocated(H_work)) deallocate(H_work)

   end subroutine compute_sturm_bounds


   !---------------------------------------------------------------------------
   ! Exact diagonalization wrapper (returns eigenvalues optionally)
   subroutine compute_exact_bounds(H, e_min, e_max, eigenvals)
      complex(rp), dimension(:, :), intent(in) :: H
      real(rp), intent(out) :: e_min, e_max
      real(rp), dimension(:), allocatable, intent(out), optional :: eigenvals

      integer :: n, info, lwork, alloc_stat
      complex(rp), allocatable :: work(:)
      real(rp), allocatable :: rwork(:)
      real(rp), allocatable :: evals(:)
      complex(rp), allocatable :: H_work(:,:)

      n = size(H, 1)
      allocate(evals(n), H_work(n,n), rwork(3*n - 2), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call g_logger%error('compute_exact_bounds: allocation failed', __FILE__, __LINE__)
         e_min = 0.0_rp
         e_max = 0.0_rp
         return
      end if

      H_work = H

      allocate(work(1), stat=alloc_stat)
      if (alloc_stat /= 0) then
         deallocate(H_work, evals, rwork)
         call g_logger%error('compute_exact_bounds: workspace alloc failed', __FILE__, __LINE__)
         e_min = 0.0_rp
         e_max = 0.0_rp
         return
      end if

      call zheev('N', 'U', n, H_work, n, evals, work, -1, rwork, info)
      if (info /= 0) then
         call g_logger%error('compute_exact_bounds: workspace query failed', __FILE__, __LINE__)
         e_min = 0.0_rp
         e_max = 0.0_rp
         if (allocated(work)) deallocate(work)
         deallocate(H_work, evals, rwork)
         return
      end if

      lwork = max(1, int(real(work(1))))
      deallocate(work)
      allocate(work(lwork))

      call zheev('N', 'U', n, H_work, n, evals, work, lwork, rwork, info)
      if (info == 0) then
         e_min = minval(evals)
         e_max = maxval(evals)
         if (present(eigenvals)) then
            allocate(eigenvals(n))
            eigenvals = evals
         end if
      else
         call g_logger%error('compute_exact_bounds: diagonalization failed', __FILE__, __LINE__)
         e_min = 0.0_rp
         e_max = 0.0_rp
      end if

      if (allocated(work)) deallocate(work)
      if (allocated(rwork)) deallocate(rwork)
      if (allocated(evals)) deallocate(evals)
      if (allocated(H_work)) deallocate(H_work)

   end subroutine compute_exact_bounds


   !---------------------------------------------------------------------------
   subroutine print_bounds_info(b)
      type(bounds), intent(in) :: b
      real(rp) :: ratio
      character(len=128) :: msg

      write(msg, '(A)') 'Spectrum Bounds Analysis:'
      call g_logger%info(msg, __FILE__, __LINE__)

      write(msg, '(A,F15.8,A,F15.8)') '  Gershgorin: [', b%e_min_gershgorin, ', ', b%e_max_gershgorin
      call g_logger%info(msg, __FILE__, __LINE__)
      write(msg, '(A,F15.8)') '    Width: ', b%gershgorin_width
      call g_logger%info(msg, __FILE__, __LINE__)

      if (b%sturm_available) then
         write(msg, '(A,F15.8,A,F15.8)') '  Sturm/Exact: [', b%e_min_sturm, ', ', b%e_max_sturm
         call g_logger%info(msg, __FILE__, __LINE__)
         write(msg, '(A,F15.8)') '    Width: ', b%sturm_width
         call g_logger%info(msg, __FILE__, __LINE__)

         ratio = b%sturm_width / b%gershgorin_width
         write(msg, '(A,F8.4)') '  Tightness ratio (Sturm/Gershgorin): ', ratio
         call g_logger%info(msg, __FILE__, __LINE__)
      else
         call g_logger%info('  Sturm sequence: NOT AVAILABLE', __FILE__, __LINE__)
      end if

      write(msg, '(A,F15.8,A,F15.8)') '  Selected: [', b%e_min, ', ', b%e_max
      call g_logger%info(msg, __FILE__, __LINE__)

   end subroutine print_bounds_info


   !---------------------------------------------------------------------------
   ! Bounds object constructor and helpers
   function bounds_constructor() result(obj)
      type(bounds), pointer :: obj
      allocate(obj)
      call obj%restore_to_default()
   end function bounds_constructor

   subroutine restore_to_default(this)
      class(bounds), intent(inout) :: this
      this%algorithm = 'none'
      this%scaling = default_scaling
      this%e_min = huge(1.0_rp)
      this%e_max = -huge(1.0_rp)
      this%e_min_gershgorin = huge(1.0_rp)
      this%e_max_gershgorin = -huge(1.0_rp)
      this%e_min_sturm = huge(1.0_rp)
      this%e_max_sturm = -huge(1.0_rp)
      this%sturm_available = .false.
      this%use_sturm = .false.
      this%sturm_width = 0.0_rp
      this%gershgorin_width = 0.0_rp
   end subroutine restore_to_default

   !> @brief Read bounds configuration from namelist (DEPRECATED)
   !> @details This method is deprecated. Configuration should be done via
   !>          the parent object's namelist (e.g., hamiltonian reads 
   !>          bounds_algorithm and bounds_scaling from /hamiltonian/ namelist).
   !>          This method is kept for backward compatibility but may be removed.
   subroutine read_from_namelist(this, fname)
      class(bounds), intent(inout) :: this
      character(len=*), intent(in) :: fname
      integer :: funit, iostat

   include 'include_codes/namelists/bounds.f90'

      call g_logger%warning('bounds%read_from_namelist is DEPRECATED. '// &
         'Use parent object namelist (e.g., bounds_algorithm in /hamiltonian/) instead.', __FILE__, __LINE__)

      open(newunit=funit, file=trim(fname), action='read', status='old', iostat=iostat)
      if (iostat /= 0) then
         call g_logger%warning('bounds%read_from_namelist: cannot open file '//trim(fname), __FILE__, __LINE__)
         return
      end if

      read(funit, nml=bounds, iostat=iostat)
      if (iostat /= 0 .and. .not. IS_IOSTAT_END(iostat)) then
         call g_logger%info('bounds%read_from_namelist: namelist /bounds/ not present or error', __FILE__, __LINE__)
      end if
      close(funit)

      ! Transfer namelist locals into object fields
      if (len_trim(bounds_algorithm) == 0) bounds_algorithm = 'none'
      if (bounds_scaling <= 0.0_rp) bounds_scaling = default_scaling
      this%algorithm = bounds_algorithm
      this%scaling = bounds_scaling

   end subroutine read_from_namelist

end module spectrum_bounds_mod
