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
      logical :: eigen_succeeded = .false.
      integer :: dimension = 0
      integer :: nearest_zero_index = 0
      real(rp) :: magnetization_norm = 0.0_rp
      real(rp) :: ward_residual = 0.0_rp
      real(rp) :: ward_relative = 0.0_rp
      real(rp) :: eigen_reconstruction_residual = huge(1.0_rp)
      real(rp) :: ward_reconstruction_residual = huge(1.0_rp)
      real(rp) :: min_singular_value = huge(1.0_rp)
      real(rp) :: second_singular_value = huge(1.0_rp)
      real(rp) :: max_singular_value = 0.0_rp
      real(rp) :: condition_number = huge(1.0_rp)
      real(rp) :: matrix_norm_2 = 0.0_rp
      real(rp) :: matrix_norm_frobenius = 0.0_rp
      real(rp) :: eigenvector_condition_estimate = huge(1.0_rp)
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
      complex(rp), allocatable :: singular_left_vectors(:, :), singular_right_vectors(:, :)
      complex(rp), allocatable :: smallest_singular_left_vector(:)
      complex(rp), allocatable :: smallest_singular_right_vector(:)
   end type lr_ward_mode_analysis_result

   public :: analyze_lr_ward_mode
   public :: lr_matrix_norms

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
      complex(rp), allocatable :: ward_vector(:), ward_sum(:)
      complex(rp) :: query(1)
      real(rp), allocatable :: rwork(:)
      real(rp) :: mnorm, dnorm, mode_norm, ward_modal_norm_sum
      integer, allocatable :: pivots(:)
      integer :: n, lwork, info, i
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
         result%singular_left_vectors(n, n), result%singular_right_vectors(n, n), &
         result%smallest_singular_left_vector(n), result%smallest_singular_right_vector(n), &
         matrix_work(n, n), left_raw(n, n), right(n, n), &
         inverse_right(n, n), reconstruction(n, n), ward_vector(n), ward_sum(n))
      result%eigenvalues = cmplx(0.0_rp, 0.0_rp, rp)
      result%right_eigenvectors = cmplx(0.0_rp, 0.0_rp, rp)
      result%left_eigenvectors = cmplx(0.0_rp, 0.0_rp, rp)
      result%coefficients = cmplx(0.0_rp, 0.0_rp, rp)
      result%ward_vector = cmplx(0.0_rp, 0.0_rp, rp)
      result%magnitude_order = 0
      result%right_overlap = 0.0_rp
      result%biorthogonal_weight = 0.0_rp
      result%magnetization_mode_fraction = 0.0_rp
      result%ward_modal_norm_diagnostic = 0.0_rp
      result%singular_values = 0.0_rp
      result%singular_vector_overlap = 0.0_rp
      result%singular_left_vectors = cmplx(0.0_rp, 0.0_rp, rp)
      result%singular_right_vectors = cmplx(0.0_rp, 0.0_rp, rp)
      result%smallest_singular_left_vector = cmplx(0.0_rp, 0.0_rp, rp)
      result%smallest_singular_right_vector = cmplx(0.0_rp, 0.0_rp, rp)

      mnorm = sqrt(sum(abs(magnetization)**2))
      result%magnetization_norm = mnorm
      if (mnorm <= tiny(1.0_rp)) return

      ! SVD is the primary diagnostic and must remain available even when the
      ! secondary non-Hermitian eigenvector problem is ill-conditioned.
      call analyze_svd(denominator, magnetization, result)
      if (.not. result%succeeded) return

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
      result%eigenvector_condition_estimate = sqrt(sum(abs(right)**2))*sqrt(sum(abs(inverse_right)**2))
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
      result%eigen_succeeded = .true.
   end subroutine analyze_lr_ward_mode

   subroutine analyze_svd(denominator, magnetization, result)
      complex(rp), intent(in) :: denominator(:, :), magnetization(:)
      type(lr_ward_mode_analysis_result), intent(inout) :: result
      complex(rp), allocatable :: work_matrix(:, :), u(:, :), vt(:, :), work(:)
      complex(rp) :: query(1)
      real(rp), allocatable :: singular_values(:), rwork(:)
      real(rp) :: scale
      integer :: n, lwork, info, i, mode

      n = size(denominator, 1)
      allocate(work_matrix(n,n), u(n,n), vt(n,n), singular_values(n), rwork(max(1,5*n)), work(1))
      work_matrix = denominator
      call zgesvd('A','A',n,n,work_matrix,n,singular_values,u,n,vt,n,work,-1,rwork,info)
      if (info /= 0) return
      lwork = max(1,nint(real(work(1),rp)))
      deallocate(work); allocate(work(lwork)); work_matrix = denominator
      call zgesvd('A','A',n,n,work_matrix,n,singular_values,u,n,vt,n,work,lwork,rwork,info)
      if (info /= 0) return
      result%singular_values = singular_values
      result%min_singular_value = minval(singular_values)
      result%max_singular_value = maxval(singular_values)
      result%matrix_norm_2 = result%max_singular_value
      result%matrix_norm_frobenius = sqrt(sum(singular_values**2))
      result%condition_number = result%max_singular_value/max(result%min_singular_value,tiny(1.0_rp))
      mode = minloc(singular_values,dim=1)
      result%second_singular_value = huge(1.0_rp)
      do i = 1,n
         if (i /= mode) result%second_singular_value = min(result%second_singular_value,singular_values(i))
         result%singular_left_vectors(:,i) = u(:,i)
         result%singular_right_vectors(:,i) = conjg(vt(i,:))
         scale = sqrt(sum(abs(result%singular_right_vectors(:,i))**2))
         result%singular_vector_overlap(i) = abs(dot_product(result%singular_right_vectors(:,i),magnetization))/ &
            max(scale*max(result%magnetization_norm,tiny(1.0_rp)),tiny(1.0_rp))
      end do
      result%smallest_singular_left_vector = result%singular_left_vectors(:,mode)
      result%smallest_singular_right_vector = result%singular_right_vectors(:,mode)
      result%succeeded = .true.
   end subroutine analyze_svd

   subroutine initialize_result(result)
      type(lr_ward_mode_analysis_result), intent(out) :: result
      result%succeeded = .false.
      result%eigen_succeeded = .false.
      result%dimension = 0
      result%nearest_zero_index = 0
      result%magnetization_norm = 0.0_rp
      result%ward_residual = 0.0_rp
      result%ward_relative = 0.0_rp
      result%eigen_reconstruction_residual = huge(1.0_rp)
      result%ward_reconstruction_residual = huge(1.0_rp)
      result%min_singular_value = huge(1.0_rp)
      result%max_singular_value = 0.0_rp
      result%second_singular_value = huge(1.0_rp)
      result%matrix_norm_2 = 0.0_rp
      result%matrix_norm_frobenius = 0.0_rp
      result%eigenvector_condition_estimate = huge(1.0_rp)
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

   subroutine lr_matrix_norms(matrix, norm_2, norm_frobenius)
      complex(rp), intent(in) :: matrix(:, :)
      real(rp), intent(out) :: norm_2, norm_frobenius
      complex(rp), allocatable :: work_matrix(:, :), u(:, :), vt(:, :), work(:)
      real(rp), allocatable :: singular_values(:), rwork(:)
      complex(rp) :: query(1)
      integer :: m, n, lwork, info

      m = size(matrix, 1); n = size(matrix, 2)
      if (m < 1 .or. n < 1) error stop 'lr_matrix_norms: empty matrix'
      allocate(work_matrix(m, n), u(m, min(m,n)), vt(min(m,n), n), singular_values(min(m,n)), &
         rwork(max(1, 5*min(m,n))), work(1))
      work_matrix = matrix
      call zgesvd('N', 'N', m, n, work_matrix, m, singular_values, u, m, vt, min(m,n), work, -1, rwork, info)
      if (info /= 0) error stop 'lr_matrix_norms: zgesvd workspace query failed'
      lwork = max(1, nint(real(work(1), rp)))
      deallocate(work); allocate(work(lwork)); work_matrix = matrix
      call zgesvd('N', 'N', m, n, work_matrix, m, singular_values, u, m, vt, min(m,n), work, lwork, rwork, info)
      if (info /= 0) error stop 'lr_matrix_norms: zgesvd failed'
      norm_2 = maxval(singular_values)
      norm_frobenius = sqrt(sum(singular_values**2))
   end subroutine lr_matrix_norms

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
