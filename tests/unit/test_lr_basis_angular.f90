!------------------------------------------------------------------------------
! LR-BASIS-00 angular representation oracle
!------------------------------------------------------------------------------
program test_lr_basis_angular
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use math_mod, only: hcpx
   implicit none

   real(rp), parameter :: tol = 2.0e-13_rp
   complex(rp), allocatable :: original(:, :), transformed(:, :)
   integer :: n

   do n = 4, 9, 5
      call basis_init(merge(1, 2, n == 4))
      allocate(original(n, n), transformed(n, n))
      call make_hermitian(original)
      transformed = original
      call hcpx(transformed, 'cart2sph')
      call hcpx(transformed, 'sph2cart')
      write (*, '(a,i0,a,es16.8)') 'sp/spd angular round-trip n=', n, ' error = ', &
         maxval(abs(transformed - original))
      if (maxval(abs(transformed - original)) > tol) then
         error stop 'LR-BASIS angular round-trip failed'
      end if
      deallocate(original, transformed)
   end do

   call basis_init(3)
   n = 16
   allocate(original(n, n), transformed(n, n))
   call make_hermitian(original)
   transformed = original
   call hcpx(transformed, 'cart2sph')
   write (*, '(a,es16.8)') 'spdf guarded identity-transform error = ', maxval(abs(transformed - original))
   if (maxval(abs(transformed - original)) > tol) then
      error stop 'LR-BASIS spdf guard changed the matrix'
   end if
   write (*, '(a)') 'UnitLrBasisAngular: PASS (sp/spd supported; spdf explicitly guarded)'

contains

   subroutine make_hermitian(matrix)
      complex(rp), intent(out) :: matrix(:, :)
      integer :: i, j

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, size(matrix, 1)
         matrix(i, i) = cmplx(0.17_rp*real(i, rp), 0.0_rp, rp)
         do j = i + 1, size(matrix, 2)
            matrix(i, j) = cmplx(0.013_rp*real(2*i + j, rp), 0.007_rp*real(j - i, rp), rp)
            matrix(j, i) = conjg(matrix(i, j))
         end do
      end do
   end subroutine make_hermitian

end program test_lr_basis_angular
