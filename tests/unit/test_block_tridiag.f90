!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
! Catches: a wrong sub/super-diagonal convention (b vs b^H), a wrong Schur
! recursion order, or a misplaced block in X_m of block_tridiag_first_column.
! Reference: dense zgesv of (zI - T) E_1 assembled in the test (fixed seed,
! general complex b so that b and b^H differ).
!------------------------------------------------------------------------------
program test_block_tridiag
   use block_tridiag_mod, only: block_tridiag_first_column
   implicit none
   integer, parameter :: dp = kind(0.0d0), nb = 4
   integer, parameter :: nblks(3) = [1, 2, 5]
   complex(dp) :: zs(2)
   integer :: iz, in, nblk, m, n, k
   logical :: failed
   real(dp) :: err
   complex(dp), allocatable :: a(:, :, :), b(:, :, :), x(:, :, :), t(:, :), zt(:, :), rhs(:, :)
   integer, allocatable :: ipiv(:)
   integer :: info

   failed = .false.
   zs = [cmplx(0.3_dp, 0.7_dp, dp), cmplx(-1.1_dp, -0.4_dp, dp)]
   call srand(20260101)
   do in = 1, size(nblks)
      nblk = nblks(in)
      n = nb*nblk
      allocate (a(nb, nb, nblk), b(nb, nb, nblk), x(nb, nb, nblk), t(n, n), zt(n, n), rhs(n, nb), ipiv(n))
      call random_blocks(a, .true.)
      call random_blocks(b, .false.)
      t = 0
      do m = 1, nblk
         t((m - 1)*nb + 1:m*nb, (m - 1)*nb + 1:m*nb) = a(:, :, m)
         if (m < nblk) then
            t(m*nb + 1:(m + 1)*nb, (m - 1)*nb + 1:m*nb) = b(:, :, m)
            t((m - 1)*nb + 1:m*nb, m*nb + 1:(m + 1)*nb) = conjg(transpose(b(:, :, m)))
         end if
      end do
      do iz = 1, 2
         call block_tridiag_first_column(a, b, zs(iz), x)
         zt = -t
         rhs = 0
         do k = 1, n
            zt(k, k) = zt(k, k) + zs(iz)
         end do
         do k = 1, nb
            rhs(k, k) = 1
         end do
         call zgesv(n, nb, zt, n, ipiv, rhs, n, info)
         err = 0
         do m = 1, nblk
            err = max(err, maxval(abs(x(:, :, m) - rhs((m - 1)*nb + 1:m*nb, :))))
         end do
         err = err/maxval(abs(rhs))
         print '(a,i0,a,i0,a,es10.3)', 'nblk=', nblk, ' z#', iz, ' rel err=', err
         if (err > 1.0e-12_dp) failed = .true.
      end do
      deallocate (a, b, x, t, zt, rhs, ipiv)
   end do
   if (failed) error stop 'test_block_tridiag FAILED'

contains

   subroutine random_blocks(c, hermitian)
      complex(dp), intent(out) :: c(:, :, :)
      logical, intent(in) :: hermitian
      integer :: mm, i, j
      real(dp) :: re, im
      do mm = 1, size(c, 3)
         do j = 1, nb
            do i = 1, nb
               re = rand() - 0.5_dp
               im = rand() - 0.5_dp
               c(i, j, mm) = cmplx(re, im, dp)
            end do
         end do
         if (hermitian) c(:, :, mm) = 0.5_dp*(c(:, :, mm) + conjg(transpose(c(:, :, mm))))
      end do
   end subroutine random_blocks

end program test_block_tridiag
