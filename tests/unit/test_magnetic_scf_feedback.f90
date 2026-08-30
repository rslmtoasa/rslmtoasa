! XCR-04: the ordinary collinear SCF radial feedback accumulator must retain
! signed and absolute spin/XC differences with the production quadrature.
program test_magnetic_scf_feedback
   use precision_mod, only: rp
   use xc_response_kernel_mod, only: xc_response_radial_projection
   implicit none

   type(xc_response_radial_projection) :: projection
   logical :: failed

   failed = .false.
   call projection%clear()
   call projection%accumulate(0.5_rp, 1.0_rp, 3.0_rp, -2.0_rp, 4.0_rp)
   call projection%accumulate(0.5_rp, 3.0_rp, 1.0_rp, 4.0_rp, -2.0_rp)

   call require_close(projection%charge_population, 4.0_rp, 'charge population')
   call require_close(projection%spin_population, 0.0_rp, 'signed spin population')
   call require_close(projection%spin_abs_population, 2.0_rp, 'absolute spin population')
   call require_close(projection%vxc_spin_difference, 0.0_rp, 'signed XC spin difference')
   call require_close(projection%vxc_spin_difference_abs, 6.0_rp, 'absolute XC spin difference')

   if (failed) error stop 1
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine require_close(actual, expected, description)
      real(rp), intent(in) :: actual, expected
      character(len=*), intent(in) :: description
      if (abs(actual - expected) > 1.0e-12_rp) then
         write (*, '(a,es16.8,a,es16.8)') 'FAIL: '//trim(description)//' actual=', actual, ' expected=', expected
         failed = .true.
      end if
   end subroutine require_close

end program test_magnetic_scf_feedback
