!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: block_tridiag_mod
!
!> First block column of (z - T)^-1 for a block-tridiagonal T, plain arrays.
!------------------------------------------------------------------------------
module block_tridiag_mod
   implicit none
   private
   public :: block_tridiag_first_column

contains

   !> @brief First block column X of (z I - T)^-1, T block tridiagonal with M diagonal blocks.
   !> @details T(m,m) = a(:,:,m), T(m+1,m) = b(:,:,m), T(m,m+1) = b(:,:,m)^H for m = 1..M-1
   !>          (b(:,:,M) is not used). Backward Schur complements D_M = z - a_M,
   !>          D_m = z - a_m - b_m^H D_{m+1}^-1 b_m; then X_1 = D_1^-1, X_{m+1} = D_{m+1}^-1 b_m X_m.
   !> @param[in] a Diagonal blocks (nb,nb,M).
   !> @param[in] b Sub-diagonal blocks (nb,nb,M); b(:,:,m) couples block m to block m+1.
   !> @param[in] z Complex energy.
   !> @param[out] x First block column (nb,nb,M).
   subroutine block_tridiag_first_column(a, b, z, x)
      complex(kind(0.0d0)), intent(in) :: a(:, :, :), b(:, :, :), z
      complex(kind(0.0d0)), intent(out) :: x(:, :, :)
      complex(kind(0.0d0)), allocatable :: d(:, :), w(:, :, :)
      integer, allocatable :: ipiv(:)
      integer :: nb, nblk, m, k, info

      nb = size(a, 1)
      nblk = size(a, 3)
      allocate (d(nb, nb), w(nb, nb, max(nblk - 1, 1)), ipiv(nb))
      d = -a(:, :, nblk)
      do k = 1, nb
         d(k, k) = d(k, k) + z
      end do
      do m = nblk - 1, 1, -1
         ! W_m = D_{m+1}^-1 b_m, then D_m = z - a_m - b_m^H W_m
         w(:, :, m) = b(:, :, m)
         call zgesv(nb, nb, d, nb, ipiv, w(:, :, m), nb, info)
         if (info /= 0) error stop 'block_tridiag_first_column: singular Schur complement'
         d = -a(:, :, m) - matmul(conjg(transpose(b(:, :, m))), w(:, :, m))
         do k = 1, nb
            d(k, k) = d(k, k) + z
         end do
      end do
      x(:, :, 1) = 0
      do k = 1, nb
         x(k, k, 1) = 1
      end do
      call zgesv(nb, nb, d, nb, ipiv, x(:, :, 1), nb, info)
      if (info /= 0) error stop 'block_tridiag_first_column: singular Schur complement'
      do m = 1, nblk - 1
         x(:, :, m + 1) = matmul(w(:, :, m), x(:, :, m))
      end do
   end subroutine block_tridiag_first_column

end module block_tridiag_mod
