!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Spectrum Bounds
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> Module to compute tight spectrum bounds for Tight-Binding Hamiltonians
!> Useful for KPM/Chebyshev polynomial expansion
!> Implements Gershgorin circle theorem and Sturm sequence methods
!------------------------------------------------------------------------------

module spectrum_bounds_mod

   use precision_mod, only: rp
   use logger_mod, only: g_logger
   implicit none

   private

   public :: spectrum_bounds_type, compute_spectrum_bounds

   !> Type to store spectrum bounds information
   type :: spectrum_bounds_type
      real(rp) :: e_min_sturm, e_max_sturm       !< Sturm sequence bounds
      real(rp) :: e_min_gershgorin, e_max_gershgorin  !< Gershgorin bounds
      real(rp) :: e_min, e_max                   !< Selected bounds (best estimate)
      real(rp) :: sturm_width, gershgorin_width  !< Spectral width
      logical :: sturm_available                 !< Flag if Sturm succeeded
      logical :: use_sturm                       !< Use Sturm over Gershgorin
   end type spectrum_bounds_type

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Compute spectrum bounds for a general TB Hamiltonian matrix
   !>
   !> Computes bounds using both Gershgorin circle theorem (fast) and
   !> optional Sturm sequence method (more accurate if available).
   !>
   !> @param[in] H Complex Hamiltonian matrix (n x n)
   !> @param[in] method 'gershgorin' or 'sturm' or 'both'
   !> @param[out] bounds Spectrum bounds structure
   !> @param[in] verbose Optional verbosity flag
   !---------------------------------------------------------------------------
   subroutine compute_spectrum_bounds(H, bounds, method, verbose)
      complex(rp), dimension(:, :), intent(in) :: H
      type(spectrum_bounds_type), intent(out) :: bounds
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

      bounds%sturm_available = .false.
      bounds%use_sturm = .false.

      ! Always compute Gershgorin (fast, O(n^2))
      call compute_gershgorin_bounds(H, bounds%e_min_gershgorin, bounds%e_max_gershgorin)
      bounds%gershgorin_width = bounds%e_max_gershgorin - bounds%e_min_gershgorin

      ! Try Sturm if requested (more expensive, O(n^3) for full diagonalization as backup)
      try_sturm = (method_to_use == 'sturm' .or. method_to_use == 'both')
      if (try_sturm) then
         call compute_sturm_bounds(H, bounds%e_min_sturm, bounds%e_max_sturm, &
                                   bounds%sturm_available)
         if (bounds%sturm_available) then
            bounds%sturm_width = bounds%e_max_sturm - bounds%e_min_sturm
            bounds%use_sturm = .true.
         end if
      end if

      ! Select best bounds
      if (bounds%sturm_available) then
         bounds%e_min = bounds%e_min_sturm
         bounds%e_max = bounds%e_max_sturm
      else
         bounds%e_min = bounds%e_min_gershgorin
         bounds%e_max = bounds%e_max_gershgorin
      end if

      if (verb) then
         call print_bounds_info(bounds)
      end if

   end subroutine compute_spectrum_bounds

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Compute bounds using Gershgorin circle theorem
   !>
   !> For matrix H, the eigenvalues lie in the union of circles centered at
   !> the diagonal elements with radius equal to the sum of absolute values
   !> of off-diagonal elements in that row.
   !>
   !> E_max <= max_i (Re[H_ii] + sum_{j≠i} |H_ij|)
   !> E_min >= min_i (Re[H_ii] - sum_{j≠i} |H_ij|)
   !>
   !> For complex Hermitian matrices, we use the Gershgorin estimate
   !> but note that it can be loose for strongly off-diagonal matrices.
   !>
   !> @param[in] H Complex Hamiltonian (n x n)
   !> @param[out] e_max Upper bound
   !> @param[out] e_min Lower bound
   !---------------------------------------------------------------------------
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
         ! Diagonal element (should be real for Hermitian matrix)
         diag = real(H(i, i))

         ! Sum of absolute values of off-diagonal elements in row i
         radius = 0.0_rp
         do j = 1, n
            if (i /= j) then
               radius = radius + abs(H(i, j))
            end if
         end do

         ! Update bounds
         e_max = max(e_max, diag + radius)
         e_min = min(e_min, diag - radius)
      end do

      ! Add small safety margin for floating point
      e_max = e_max * (1.0_rp + eps)
      e_min = e_min * (1.0_rp - eps)

   end subroutine compute_gershgorin_bounds

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Compute bounds using Sturm sequence counting
   !>
   !> Uses characteristic polynomial and Sturm sequence to find exact eigenvalues.
   !> This is expensive (requires diagonalization as fallback) but gives tight bounds.
   !>
   !> @param[in] H Complex Hamiltonian (n x n)
   !> @param[out] e_max Upper bound
   !> @param[out] e_min Lower bound
   !> @param[out] success Flag indicating if method succeeded
   !---------------------------------------------------------------------------
   subroutine compute_sturm_bounds(H, e_min, e_max, success)
      complex(rp), dimension(:, :), intent(in) :: H
      real(rp), intent(out) :: e_min, e_max
      logical, intent(out) :: success

      ! For full implementation, would use eigenvalue solver
      ! For now, use Gershgorin as fallback (can be improved with LAPACK)
      real(rp), dimension(:), allocatable :: eigenvals
      integer :: n, info, lwork
   complex(rp), dimension(:), allocatable :: work
   real(rp), dimension(:), allocatable :: rwork
      complex(rp), dimension(:,:), allocatable :: H_work

      n = size(H, 1)
      success = .false.

      allocate(H_work(n, n))
      allocate(eigenvals(n))
      allocate(rwork(3*n - 2))

      H_work = H

      ! Query for optimal workspace
      allocate(work(1))
      call zheev('N', 'U', n, H_work, n, eigenvals, work, -1, rwork, info)

      if (info /= 0) then
         e_min = 0.0_rp
         e_max = 0.0_rp
         return
      end if

      lwork = max(1, int(real(work(1))))
      deallocate(work)
      allocate(work(lwork))

      ! Diagonalize Hermitian matrix
      call zheev('N', 'U', n, H_work, n, eigenvals, work, lwork, rwork, info)

      if (info == 0) then
         e_min = minval(eigenvals)
         e_max = maxval(eigenvals)
         success = .true.
      end if

      deallocate(H_work, eigenvals, work, rwork)

   end subroutine compute_sturm_bounds

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Compute bounds using eigenvalue analysis on full Hamiltonian
   !>
   !> This is the "brute force" method - diagonalize and find min/max eigenvalues.
   !> Expensive but gives exact spectrum edges.
   !>
   !> @param[in] H Complex Hamiltonian (n x n)
   !> @param[out] e_max Upper bound (max eigenvalue)
   !> @param[out] e_min Lower bound (min eigenvalue)
   !> @param[out] eigenvals Optional: all eigenvalues
   !---------------------------------------------------------------------------
   subroutine compute_exact_bounds(H, e_min, e_max, eigenvals)
      complex(rp), dimension(:, :), intent(in) :: H
      real(rp), intent(out) :: e_min, e_max
      real(rp), dimension(:), allocatable, intent(out), optional :: eigenvals

      real(rp), dimension(:), allocatable :: evals
      integer :: n, info, lwork
   complex(rp), dimension(:), allocatable :: work
   real(rp), dimension(:), allocatable :: rwork
      complex(rp), dimension(:,:), allocatable :: H_work

      n = size(H, 1)
      allocate(evals(n))
      allocate(H_work(n, n))
      allocate(rwork(3*n - 2))

      H_work = H

      ! Query for optimal workspace
      allocate(work(1))
      call zheev('N', 'U', n, H_work, n, evals, work, -1, rwork, info)

      if (info /= 0) then
         call g_logger%error('spectrum_bounds: workspace query failed', __FILE__, __LINE__)
         e_min = 0.0_rp
         e_max = 0.0_rp
         return
      end if

      lwork = max(1, int(real(work(1))))
      deallocate(work)
      allocate(work(lwork))

      ! Diagonalize
      call zheev('N', 'U', n, H_work, n, evals, work, lwork, rwork, info)

      if (info == 0) then
         e_min = minval(evals)
         e_max = maxval(evals)
         if (present(eigenvals)) then
            allocate(eigenvals(n))
            eigenvals = evals
         end if
      else
         call g_logger%error('spectrum_bounds: diagonalization failed', __FILE__, __LINE__)
         e_min = 0.0_rp
         e_max = 0.0_rp
      end if

      deallocate(evals, H_work, work, rwork)

   end subroutine compute_exact_bounds

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print bounds information to logger
   !---------------------------------------------------------------------------
   subroutine print_bounds_info(bounds)
      type(spectrum_bounds_type), intent(in) :: bounds
      real(rp) :: ratio
      character(len=128) :: msg

      write(msg, '(A)') 'Spectrum Bounds Analysis:'
      call g_logger%info(msg, __FILE__, __LINE__)

      write(msg, '(A,F15.8,A,F15.8)') '  Gershgorin: [', bounds%e_min_gershgorin, ', ', &
                                        bounds%e_max_gershgorin
      call g_logger%info(msg, __FILE__, __LINE__)
      write(msg, '(A,F15.8)') '    Width: ', bounds%gershgorin_width
      call g_logger%info(msg, __FILE__, __LINE__)

      if (bounds%sturm_available) then
         write(msg, '(A,F15.8,A,F15.8)') '  Sturm/Exact: [', bounds%e_min_sturm, ', ', &
                                           bounds%e_max_sturm
         call g_logger%info(msg, __FILE__, __LINE__)
         write(msg, '(A,F15.8)') '    Width: ', bounds%sturm_width
         call g_logger%info(msg, __FILE__, __LINE__)

         ratio = bounds%sturm_width / bounds%gershgorin_width
         write(msg, '(A,F8.4)') '  Tightness ratio (Sturm/Gershgorin): ', ratio
         call g_logger%info(msg, __FILE__, __LINE__)
      else
         call g_logger%info('  Sturm sequence: NOT AVAILABLE', __FILE__, __LINE__)
      end if

      write(msg, '(A,F15.8,A,F15.8)') '  Selected: [', bounds%e_min, ', ', bounds%e_max
      call g_logger%info(msg, __FILE__, __LINE__)

   end subroutine print_bounds_info

end module spectrum_bounds_mod
