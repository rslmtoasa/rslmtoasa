! XCR-03: authoritative contract tests for the logarithmic radial derivative.
program test_radial_gga
   use precision_mod, only: rp
   use xc_radial_mod, only: radgra, radial_flux_divergence
   implicit none

   logical :: failed

   failed = .false.
   call check_derivative_convergence()
   call check_regular_origin()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine check_derivative_convergence()
      integer, parameter :: ncoarse = 80, nfine = 160
      real(rp) :: coarse_first, coarse_second, fine_first, fine_second

      call derivative_errors(ncoarse, coarse_first, coarse_second)
      call derivative_errors(nfine, fine_first, fine_second)
      write (*, '(a,4(es12.4,1x))') 'radgra max errors (coarse/fine, first/second): ', &
         coarse_first, fine_first, coarse_second, fine_second

      call require(fine_first < coarse_first, 'radgra first derivative converges under mesh refinement')
      call require(fine_second < coarse_second, 'radgra second derivative converges under mesh refinement')
      call require(fine_first < 2.e-6_rp, 'radgra first derivative analytic accuracy')
      call require(fine_second < 5.e-5_rp, 'radgra second derivative analytic accuracy')
   end subroutine check_derivative_convergence

   subroutine derivative_errors(nr, first_error, second_error)
      integer, intent(in) :: nr
      real(rp), intent(out) :: first_error, second_error
      real(rp), parameter :: b = 0.02_rp, rmax = 5.0_rp, alpha = 0.8_rp
      real(rp) :: a, r
      real(rp), allocatable :: rofi(:), f(:), fp(:), fpp(:), exact_fp(:), exact_fpp(:)
      integer :: i

      a = log(1.0_rp + rmax/b)/real(nr - 1, rp)
      allocate (rofi(nr), f(nr), fp(nr), fpp(nr), exact_fp(nr), exact_fpp(nr))
      do i = 1, nr
         rofi(i) = b*(exp(a*real(i - 1, rp)) - 1.0_rp)
         r = rofi(i)
         f(i) = exp(-alpha*r)
         exact_fp(i) = -alpha*f(i)
         exact_fpp(i) = alpha*alpha*f(i)
      end do

      call radgra(a, b, nr, rofi, f, fp)
      call radgra(a, b, nr, rofi, fp, fpp)
      first_error = maxval(abs(fp(3:nr-2) - exact_fp(3:nr-2)))
      second_error = maxval(abs(fpp(3:nr-2) - exact_fpp(3:nr-2)))

      ! Repeat the contract check for a function with a non-trivial radial
      ! prefactor, f(r)=r2*exp(-alpha*r).
      do i = 1, nr
         r = rofi(i)
         f(i) = r*r*exp(-alpha*r)
         exact_fp(i) = exp(-alpha*r)*(2.0_rp*r - alpha*r*r)
         exact_fpp(i) = exp(-alpha*r)*(2.0_rp - 4.0_rp*alpha*r + alpha*alpha*r*r)
      end do
      call radgra(a, b, nr, rofi, f, fp)
      call radgra(a, b, nr, rofi, fp, fpp)
      first_error = max(first_error, maxval(abs(fp(3:nr-2) - exact_fp(3:nr-2))))
      second_error = max(second_error, maxval(abs(fpp(3:nr-2) - exact_fpp(3:nr-2))))
      deallocate (rofi, f, fp, fpp, exact_fp, exact_fpp)
   end subroutine derivative_errors

   subroutine check_regular_origin()
      integer, parameter :: nr = 160
      real(rp), parameter :: a = 0.04_rp, b = 0.02_rp, beta = 0.7_rp
      real(rp) :: r
      real(rp), allocatable :: rofi(:), f(:), fp(:), flux(:), flux_div(:)
      integer :: i

      allocate (rofi(nr), f(nr), fp(nr), flux(nr), flux_div(nr))
      do i = 1, nr
         rofi(i) = b*(exp(a*real(i - 1, rp)) - 1.0_rp)
         r = rofi(i)
         f(i) = r*r*exp(-beta*r*r)
         flux(i) = 2.0_rp*r*exp(-beta*r*r)
      end do
      call radgra(a, b, nr, rofi, f, fp)
      call radial_flux_divergence(a, b, rofi, flux, flux_div)

      ! f'(0)=0.  For this regular radial field f'(r) rhat, the divergence
      ! has the analytic origin limit 3*f''(0)=6.
      call require(abs(fp(1)) < 2.e-8_rp, 'regular radial derivative is zero at the origin')
      call require(abs(flux_div(1) - 6.0_rp) < 2.e-3_rp, &
         'radial flux divergence uses the analytic origin limit')
      do i = 3, nr-2
         r = rofi(i)
         call require(abs(flux_div(i) - (6.0_rp*exp(-beta*r*r) - &
            4.0_rp*beta*r*r*exp(-beta*r*r))) < 6.e-3_rp, &
            'radial flux divergence agrees away from the origin')
      end do
      deallocate (rofi, f, fp, flux, flux_div)
   end subroutine check_regular_origin

   subroutine require(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description
      if (.not. condition) then
         write (*, '(a)') 'FAIL: '//trim(description)
         failed = .true.
      end if
   end subroutine require

end program test_radial_gga
