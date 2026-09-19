!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Non-Hermitian Ward-denominator eigenmode diagnostics.
!>
!> This module is diagnostic only.  It forms no corrected kernel and never
!> changes a production response.  Right eigenvectors are obtained from ZGEEV;
!> the inverse of the right-eigenvector matrix supplies the biorthogonal left
!> bras used for modal coefficients.  An independent SVD is retained because
!> finite-eta Ward denominators need not be normal.
!------------------------------------------------------------------------------
module lr_ward_mode_analysis_mod

   use precision_mod, only: rp
   implicit none
   private

   type, public :: lr_ward_mode_analysis_result
      logical :: succeeded = .false.
      integer :: dimension = 0
      integer :: nearest_zero_index = 0
      real(rp) :: magnetization_norm = 0.0_rp
      real(rp) :: ward_residual = 0.0_rp
      real(rp) :: ward_relative = 0.0_rp
      real(rp) :: eigen_reconstruction_residual = huge(1.0_rp)
      real(rp) :: ward_reconstruction_residual = huge(1.0_rp)
      real(rp) :: min_singular_value = huge(1.0_rp)
      real(rp) :: max_singular_value = 0.0_rp
      real(rp) :: condition_number = huge(1.0_rp)
      real(rp) :: minimum_magnitude_eigenvalue = huge(1.0_rp)
      integer, allocatable :: magnitude_order(:)
      complex(rp), allocatable :: eigenvalues(:)
      complex(rp), allocatable :: right_eigenvectors(:, :)
      complex(rp), allocatable :: left_eigenvectors(:, :) ! biorthogonal left kets
      complex(rp), allocatable :: coefficients(:)
      complex(rp), allocatable :: ward_vector(:), ward_sum(:)
      real(rp), allocatable :: right_overlap(:)
      real(rp), allocatable :: biorthogonal_weight(:)
      real(rp), allocatable :: magnetization_mode_fraction(:)
      ! Individual modal Euclidean norm squared divided by the sum of modal
      ! norm squares.  This is intentionally not a fraction of ||D m||^2,
      ! because the non-normal modal vectors interfere under the sum.
      real(rp), allocatable :: ward_modal_norm_diagnostic(:)
      real(rp), allocatable :: singular_values(:)
      real(rp), allocatable :: singular_vector_overlap(:)
      complex(rp), allocatable :: smallest_singular_right_vector(:)
   end type lr_ward_mode_analysis_result

   public :: analyze_lr_ward_mode

   interface
      subroutine zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
         import :: rp
         character(len=1), intent(in) :: jobvl, jobvr
         integer, intent(in) :: n, lda, ldvl, ldvr, lwork
         complex(rp), intent(inout) :: a(lda, *)
         complex(rp), intent(out) :: w(*), vl(ldvl, *), vr(ldvr, *), work(*)
         real(rp), intent(out) :: rwork(*)
         integer, intent(out) :: info
      end subroutine zgeev

      subroutine zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
         import :: rp
         character(len=1), intent(in) :: jobu, jobvt
         integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
         complex(rp), intent(inout) :: a(lda, *)
         real(rp), intent(out) :: s(*)
         complex(rp), intent(out) :: u(ldu, *), vt(ldvt, *), work(*)
         real(rp), intent(inout) :: rwork(*)
         integer, intent(out) :: info
      end subroutine zgesvd

      subroutine zgetrf(m, n, a, lda, ipiv, info)
         import :: rp
         integer, intent(in) :: m, n, lda
         complex(rp), intent(inout) :: a(lda, *)
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
      end subroutine zgetrf

      subroutine zgetri(n, a, lda, ipiv, work, lwork, info)
         import :: rp
         integer, intent(in) :: n, lda, lwork
         complex(rp), intent(inout) :: a(lda, *), work(*)
         integer, intent(in) :: ipiv(*)
         integer, intent(out) :: info
      end subroutine zgetri
   end interface

contains

   subroutine analyze_lr_ward_mode(denominator, magnetization, result)
      complex(rp), intent(in) :: denominator(:, :)
      complex(rp), intent(in) :: magnetization(:)
      type(lr_ward_mode_analysis_result), intent(out) :: result

      complex(rp), allocatable :: matrix_work(:, :), left_raw(:, :), right(:, :), inverse_right(:, :)
      complex(rp), allocatable :: work(:), reconstruction(:, :), mode_vector(:)
      complex(rp), allocatable :: svd_work_matrix(:, :), vt(:, :), svd_work(:)
      complex(rp), allocatable :: ward_vector(:), ward_sum(:)
      complex(rp) :: query(1)
      real(rp), allocatable :: rwork(:), singular_values(:), svd_rwork(:)
      real(rp) :: mnorm, dnorm, mode_norm, ward_modal_norm_sum, scale
      integer, allocatable :: pivots(:)
      integer :: n, lwork, info, i, j, mode
      logical :: inverse_ok

      call initialize_result(result)
      n = size(denominator, 1)
      if (n < 1 .or. size(denominator, 2) /= n .or. size(magnetization) /= n) then
         result%succeeded = .false.
         return
      end if
      result%dimension = n
      allocate(result%eigenvalues(n), result%right_eigenvectors(n, n), result%left_eigenvectors(n, n), &
         result%coefficients(n), result%ward_vector(n), result%magnitude_order(n), &
         result%right_overlap(n), result%biorthogonal_weight(n), result%magnetization_mode_fraction(n), &
         result%ward_modal_norm_diagnostic(n), result%singular_values(n), result%singular_vector_overlap(n), &
         result%smallest_singular_right_vector(n), matrix_work(n, n), left_raw(n, n), right(n, n), &
         inverse_right(n, n), reconstruction(n, n), ward_vector(n), ward_sum(n))

      mnorm = sqrt(sum(abs(magnetization)**2))
      result%magnetization_norm = mnorm
      if (mnorm <= tiny(1.0_rp)) return

      matrix_work = denominator
      allocate(rwork(max(1, 2*n)))
      call zgeev('N', 'V', n, matrix_work, n, result%eigenvalues, left_raw, n, right, n, query, -1, rwork, info)
      if (info /= 0) return
      lwork = max(1, nint(real(query(1), rp)))
      allocate(work(lwork))
      matrix_work = denominator
      call zgeev('N', 'V', n, matrix_work, n, result%eigenvalues, left_raw, n, right, n, work, lwork, rwork, info)
      if (info /= 0) return
      result%right_eigenvectors = right

      call invert_square_matrix(right, inverse_right, inverse_ok)
      if (.not. inverse_ok) return
      ! If R is the right-eigenvector matrix, R^-1 supplies the rows of the
      ! biorthogonal left bras.  Store their corresponding left kets.
      result%left_eigenvectors = conjg(transpose(inverse_right))
      result%coefficients = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, n
         result%coefficients(i) = dot_product(result%left_eigenvectors(:, i), magnetization)
      end do

      reconstruction = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, n
         reconstruction(:, i) = result%right_eigenvectors(:, i)*result%coefficients(i)
      end do
      result%eigen_reconstruction_residual = sqrt(sum(abs(sum(reconstruction, dim=2) - magnetization)**2))/mnorm

      ward_vector = matmul(denominator, magnetization)
      result%ward_vector = ward_vector
      dnorm = sqrt(sum(abs(ward_vector)**2))
      result%ward_residual = dnorm
      result%ward_relative = dnorm/mnorm
      result%ward_reconstruction_residual = 0.0_rp
      ward_modal_norm_sum = 0.0_rp
      do i = 1, n
         mode_vector = result%right_eigenvectors(:, i)*result%coefficients(i)*result%eigenvalues(i)
         ward_modal_norm_sum = ward_modal_norm_sum + sum(abs(mode_vector)**2)
      end do
      if (ward_modal_norm_sum > tiny(1.0_rp)) then
         do i = 1, n
            mode_vector = result%right_eigenvectors(:, i)*result%coefficients(i)*result%eigenvalues(i)
            result%ward_modal_norm_diagnostic(i) = sum(abs(mode_vector)**2)/ward_modal_norm_sum
         end do
      end if
      ward_sum = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, n
         ward_sum = ward_sum + result%right_eigenvectors(:, i)*result%coefficients(i)*result%eigenvalues(i)
      end do
      result%ward_reconstruction_residual = sqrt(sum(abs(ward_sum - ward_vector)**2))/max(dnorm, tiny(1.0_rp))

      result%minimum_magnitude_eigenvalue = minval(abs(result%eigenvalues))
      result%nearest_zero_index = minloc(abs(result%eigenvalues), dim=1)
      do i = 1, n
         result%magnitude_order(i) = i
      end do
      call sort_mode_order(result%magnitude_order, result%eigenvalues)

      do i = 1, n
         mode_norm = sqrt(sum(abs(result%right_eigenvectors(:, i))**2))
         result%right_overlap(i) = abs(dot_product(result%right_eigenvectors(:, i), magnetization))/ &
            max(mode_norm*mnorm, tiny(1.0_rp))
         result%biorthogonal_weight(i) = abs(result%coefficients(i))*mode_norm/mnorm
         result%magnetization_mode_fraction(i) = sum(abs(result%right_eigenvectors(:, i)*result%coefficients(i))**2)/mnorm**2
      end do

      allocate(svd_work_matrix(n, n), vt(n, n), singular_values(n), svd_rwork(max(1, 5*n)), svd_work(1))
      svd_work_matrix = denominator
      call zgesvd('N', 'A', n, n, svd_work_matrix, n, singular_values, reconstruction, n, vt, n, svd_work, -1, svd_rwork, info)
      if (info /= 0) return
      lwork = max(1, nint(real(svd_work(1), rp)))
      deallocate(svd_work)
      allocate(svd_work(lwork))
      svd_work_matrix = denominator
      call zgesvd('N', 'A', n, n, svd_work_matrix, n, singular_values, reconstruction, n, vt, n, svd_work, lwork, svd_rwork, info)
      if (info /= 0) return
      result%singular_values = singular_values
      result%min_singular_value = minval(singular_values)
      result%max_singular_value = maxval(singular_values)
      result%condition_number = result%max_singular_value/max(result%min_singular_value, tiny(1.0_rp))
      mode = minloc(singular_values, dim=1)
      do i = 1, n
         result%smallest_singular_right_vector(i) = conjg(vt(mode, i))
      end do
      scale = sqrt(sum(abs(result%smallest_singular_right_vector)**2))
      result%singular_vector_overlap = 0.0_rp
      result%singular_vector_overlap(mode) = abs(dot_product(result%smallest_singular_right_vector, magnetization))/ &
         max(scale*mnorm, tiny(1.0_rp))

      result%succeeded = .true.
   end subroutine analyze_lr_ward_mode

   subroutine initialize_result(result)
      type(lr_ward_mode_analysis_result), intent(out) :: result
      result%succeeded = .false.
      result%dimension = 0
      result%nearest_zero_index = 0
      result%magnetization_norm = 0.0_rp
      result%ward_residual = 0.0_rp
      result%ward_relative = 0.0_rp
      result%eigen_reconstruction_residual = huge(1.0_rp)
      result%ward_reconstruction_residual = huge(1.0_rp)
      result%min_singular_value = huge(1.0_rp)
      result%max_singular_value = 0.0_rp
      result%condition_number = huge(1.0_rp)
      result%minimum_magnitude_eigenvalue = huge(1.0_rp)
   end subroutine initialize_result

   subroutine invert_square_matrix(input, inverse, success)
      complex(rp), intent(in) :: input(:, :)
      complex(rp), intent(out) :: inverse(:, :)
      logical, intent(out) :: success
      integer, allocatable :: pivots(:)
      complex(rp), allocatable :: work(:)
      integer :: n, info, lwork

      n = size(input, 1)
      inverse = input
      allocate(pivots(n), work(max(1, n*n)))
      call zgetrf(n, n, inverse, n, pivots, info)
      if (info /= 0) then
         success = .false.
         return
      end if
      lwork = max(1, n*n)
      call zgetri(n, inverse, n, pivots, work, lwork, info)
      success = info == 0
   end subroutine invert_square_matrix

   subroutine sort_mode_order(order, eigenvalues)
      integer, intent(inout) :: order(:)
      complex(rp), intent(in) :: eigenvalues(:)
      integer :: i, j, best, temporary
      do i = 1, size(order) - 1
         best = i
         do j = i + 1, size(order)
            if (abs(eigenvalues(order(j))) < abs(eigenvalues(order(best)))) best = j
         end do
         temporary = order(i)
         order(i) = order(best)
         order(best) = temporary
      end do
   end subroutine sort_mode_order

end module lr_ward_mode_analysis_mod
