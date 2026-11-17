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
!> Configuration is read by the parent object (hamiltonian reads
!> bounds_algorithm and bounds_scaling from the &hamiltonian namelist).
!------------------------------------------------------------------------------

module spectrum_bounds_mod

   use precision_mod, only: rp
   use logger_mod, only: g_logger
   use string_mod, only: lower
   implicit none

   real(rp), parameter :: default_scaling = 1.05_rp
   private

   public :: compute_spectrum_bounds
   public :: normalize_bounds_algorithm, select_bounds_interval, apply_bounds_scaling

   !> Single type to store spectrum bounds and user options
   type, public :: bounds
      ! Numerical results (from Gershgorin/Sturm/exact)
      real(rp) :: e_min_sturm = 0.0_rp
      real(rp) :: e_max_sturm = 0.0_rp
      real(rp) :: e_min_gershgorin = 0.0_rp
      real(rp) :: e_max_gershgorin = 0.0_rp
      real(rp) :: e_min = 0.0_rp
      real(rp) :: e_max = 0.0_rp
      logical :: sturm_available = .false.

      ! User-configurable options (set via parent object’s namelist)
      character(len=16) :: algorithm = 'none'
      real(rp) :: scaling = default_scaling

   contains
      procedure :: restore_to_default
   end type bounds

contains

   !---------------------------------------------------------------------------
   ! Normalize bounds algorithm keyword.
   subroutine normalize_bounds_algorithm(algo_in, algo_out)
      character(len=*), intent(in) :: algo_in
      character(len=*), intent(out) :: algo_out
      character(len=32) :: tmp

      tmp = adjustl(trim(algo_in))
      if (len_trim(tmp) == 0) tmp = 'none'
      tmp = lower(tmp)
      if (tmp == 'hybrid') tmp = 'both'
      if (tmp /= 'none' .and. tmp /= 'gershgorin' .and. tmp /= 'hgamma' .and. tmp /= 'sturm' .and. tmp /= 'both') then
         tmp = 'none'
      end if
      algo_out = tmp
   end subroutine normalize_bounds_algorithm

   !---------------------------------------------------------------------------
   ! Select active interval based on algorithm and available estimators.
   subroutine select_bounds_interval(algo, gmin, gmax, have_hgamma, hmin, hmax, e_min, e_max)
      character(len=*), intent(in) :: algo
      real(rp), intent(in) :: gmin, gmax, hmin, hmax
      logical, intent(in) :: have_hgamma
      real(rp), intent(out) :: e_min, e_max

      select case (trim(adjustl(algo)))
      case ('gershgorin')
         e_min = gmin
         e_max = gmax
      case ('hgamma', 'sturm')
         if (have_hgamma) then
            e_min = hmin
            e_max = hmax
         else
            e_min = gmin
            e_max = gmax
         end if
      case default
         if (have_hgamma) then
            e_min = min(gmin, hmin)
            e_max = max(gmax, hmax)
         else
            e_min = gmin
            e_max = gmax
         end if
      end select
   end subroutine select_bounds_interval

   !---------------------------------------------------------------------------
   ! Symmetric scaling around interval midpoint.
   subroutine apply_bounds_scaling(e_min, e_max, scaling)
      real(rp), intent(inout) :: e_min, e_max
      real(rp), intent(in) :: scaling
      real(rp) :: width, emid, scale_use

      width = e_max - e_min
      if (width <= 0.0_rp) return
      scale_use = scaling
      if (scale_use <= 0.0_rp) scale_use = default_scaling
      emid = 0.5_rp*(e_min + e_max)
      e_min = emid - 0.5_rp*scale_use*width
      e_max = emid + 0.5_rp*scale_use*width
   end subroutine apply_bounds_scaling

   !---------------------------------------------------------------------------
   ! Compute spectrum bounds (Gershgorin + optional Sturm/Exact)
   subroutine compute_spectrum_bounds(H, res, method)
      complex(rp), dimension(:, :), intent(in) :: H
   type(bounds), intent(out) :: res
      character(len=*), intent(in), optional :: method

      character(len=32) :: method_to_use
      logical :: try_sturm
      integer :: n

      n = size(H, 1)
      if (size(H, 2) /= n) then
         call g_logger%error('spectrum_bounds: Hamiltonian must be square', __FILE__, __LINE__)
         return
      end if

      method_to_use = 'both'
      if (present(method)) method_to_use = trim(adjustl(method))

   res%sturm_available = .false.

      ! Always compute Gershgorin (fast, O(n^2))
   call compute_gershgorin_bounds(H, res%e_min_gershgorin, res%e_max_gershgorin)

      ! Try Sturm if requested (diagonalization-based fallback)
      try_sturm = (method_to_use == 'sturm' .or. method_to_use == 'both')
      if (try_sturm) then
         call compute_sturm_bounds(H, res%e_min_sturm, res%e_max_sturm, res%sturm_available)
      end if

      ! Select best bounds
      if (res%sturm_available) then
         res%e_min = res%e_min_sturm
         res%e_max = res%e_max_sturm
      else
         res%e_min = res%e_min_gershgorin
         res%e_max = res%e_max_gershgorin
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
   end subroutine restore_to_default

end module spectrum_bounds_mod
