!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
!> @brief Finite-system oracle for the native complex-energy Chebyshev GF.
!> @details Moments are generated independently by a dense matrix Chebyshev
!>          recurrence.  The production reconstruction is then compared with
!>          an independently computed dense inverse at several positive-eta
!>          complex energies.  This is a convergence test for the truncated
!>          Jackson reconstruction, not an assertion of exact finite-order
!>          equality.
!------------------------------------------------------------------------------
program test_chebyshev_gf_oracle
   use precision_mod, only: rp
   use chebyshev_fast_mod, only: cheb_green_complex
   use dyson_kernel_mod, only: dyson_kspace_inverse
   implicit none

   integer, parameter :: nmat = 4
   integer, parameter :: orders(4) = [16, 32, 64, 128]
   real(rp), parameter :: a = 1.2_rp, b = 0.0_rp
   real(rp), parameter :: tolerance = 0.15_rp
   complex(rp), parameter :: z_values(3) = [cmplx(0.10_rp, 0.20_rp, rp), &
      cmplx(-0.80_rp, 0.15_rp, rp), cmplx(1.30_rp, 0.20_rp, rp)]
   complex(rp) :: h(nmat, nmat), sigma(nmat, nmat), g_ref(nmat, nmat), g_cheb(nmat, nmat)
   complex(rp), allocatable :: moments(:, :, :)
   real(rp) :: errors(size(orders)), final_error
   integer :: i, n, iz
   logical :: failed

   failed = .false.
   call build_fixture(h)
   sigma = (0.0_rp, 0.0_rp)

   do i = 1, size(orders)
      n = orders(i)
      allocate(moments(nmat, nmat, n))
      call dense_chebyshev_moments(h, a, b, n, moments)
      errors(i) = 0.0_rp
      do iz = 1, size(z_values)
         call cheb_green_complex(moments, nmat, n, z_values(iz), a, b, g_cheb)
         call dyson_kspace_inverse(h, z_values(iz), sigma, g_ref)
         errors(i) = max(errors(i), maxval(abs(g_cheb-g_ref)))
      end do
      write (*, '(a,i0,a,es12.4)') 'Chebyshev order ', n, ' max dense-inverse error = ', errors(i)
      deallocate(moments)
   end do

   final_error = errors(size(errors))
   if (any(errors(2:) >= errors(:size(errors)-1))) then
      write (*, '(a)') 'FAIL: Chebyshev GF oracle did not converge monotonically on the fixture'
      failed = .true.
   end if
   if (final_error > tolerance) then
      write (*, '(a,es12.4)') 'FAIL: final Chebyshev GF error exceeds convergence envelope: ', final_error
      failed = .true.
   end if

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   subroutine build_fixture(matrix)
      complex(rp), intent(out) :: matrix(:, :)

      matrix = (0.0_rp, 0.0_rp)
      matrix(1, 1) = (-0.70_rp, 0.0_rp)
      matrix(2, 2) = (-0.20_rp, 0.0_rp)
      matrix(3, 3) = (0.35_rp, 0.0_rp)
      matrix(4, 4) = (0.80_rp, 0.0_rp)
      matrix(1, 2) = (0.12_rp, 0.08_rp); matrix(2, 1) = conjg(matrix(1, 2))
      matrix(1, 3) = (0.03_rp, 0.0_rp); matrix(3, 1) = conjg(matrix(1, 3))
      matrix(2, 3) = (0.10_rp, 0.0_rp); matrix(3, 2) = conjg(matrix(2, 3))
      matrix(2, 4) = (0.0_rp, 0.04_rp); matrix(4, 2) = conjg(matrix(2, 4))
      matrix(3, 4) = (0.09_rp, -0.02_rp); matrix(4, 3) = conjg(matrix(3, 4))
   end subroutine build_fixture

   subroutine dense_chebyshev_moments(matrix, scale, shift, n_mom, moments_out)
      complex(rp), intent(in) :: matrix(:, :)
      real(rp), intent(in) :: scale, shift
      integer, intent(in) :: n_mom
      complex(rp), intent(out) :: moments_out(:, :, :)
      complex(rp) :: scaled(nmat, nmat), identity(nmat, nmat)
      integer :: i, order

      identity = (0.0_rp, 0.0_rp)
      do i = 1, nmat
         identity(i, i) = (1.0_rp, 0.0_rp)
      end do
      scaled = (matrix-shift*identity)/scale
      moments_out(:, :, 1) = identity
      if (n_mom == 1) return
      moments_out(:, :, 2) = scaled
      do order = 3, n_mom
         moments_out(:, :, order) = 2.0_rp*matmul(scaled, moments_out(:, :, order-1)) - moments_out(:, :, order-2)
      end do
   end subroutine dense_chebyshev_moments

end program test_chebyshev_gf_oracle
