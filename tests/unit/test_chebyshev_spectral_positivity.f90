!------------------------------------------------------------------------------
! Unit checks for the local spectral-density projection used after real-axis
! Chebyshev reconstruction. A physical 2x2 spin spectral block is positive
! semidefinite: both diagonal DOS values are non-negative and the coherence is
! bounded by their geometric mean.
!------------------------------------------------------------------------------
program test_chebyshev_spectral_positivity
   use precision_mod, only: rp
   use math_mod, only: project_spectral_density_2x2
   implicit none

   real(rp) :: a_uu, a_dd, determinant
   complex(rp) :: a_ud
   logical :: corrected, failed

   failed = .false.

   ! A negative diagonal DOS must be removed, with coherence forced to zero
   ! when its remaining diagonal partner cannot support it.
   a_uu = -2.0e-7_rp
   a_dd = 0.4_rp
   a_ud = cmplx(0.03_rp, -0.01_rp, rp)
   call project_spectral_density_2x2(a_uu, a_ud, a_dd, corrected)
   call check_logical('negative diagonal is corrected', corrected, .true., failed)
   call check_real('negative diagonal becomes zero', a_uu, 0.0_rp, 0.0_rp, failed)
   call check_real('unsupported coherence becomes zero', abs(a_ud), 0.0_rp, 0.0_rp, failed)

   ! Retain diagonal DOS and reduce only an overlarge spin coherence to the
   ! positive-semidefinite boundary.
   a_uu = 0.09_rp
   a_dd = 0.16_rp
   a_ud = cmplx(0.20_rp, 0.0_rp, rp)
   call project_spectral_density_2x2(a_uu, a_ud, a_dd, corrected)
   determinant = a_uu*a_dd - abs(a_ud)**2
   call check_logical('overlarge coherence is corrected', corrected, .true., failed)
   call check_real('coherence is limited by sqrt(a_uu*a_dd)', abs(a_ud), 0.12_rp, 1.0e-14_rp, failed)
   call check_true('projected block is positive semidefinite', determinant >= -1.0e-14_rp, failed)

   ! A physical block must pass unchanged.
   a_uu = 0.09_rp
   a_dd = 0.16_rp
   a_ud = cmplx(0.06_rp, -0.04_rp, rp)
   call project_spectral_density_2x2(a_uu, a_ud, a_dd, corrected)
   call check_logical('physical block is unchanged', corrected, .false., failed)
   call check_real('physical a_uu is retained', a_uu, 0.09_rp, 0.0_rp, failed)
   call check_real('physical a_dd is retained', a_dd, 0.16_rp, 0.0_rp, failed)
   call check_real('physical coherence is retained', abs(a_ud), sqrt(0.0052_rp), 1.0e-14_rp, failed)

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine check_real(label, actual, expected, tolerance, failed)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: actual, expected, tolerance
      logical, intent(inout) :: failed

      if (abs(actual - expected) > tolerance) then
         write (*, '(a,1x,a,2(a,es14.6))') 'FAIL', trim(label), ' actual=', actual, ' expected=', expected
         failed = .true.
      end if
   end subroutine check_real

   subroutine check_logical(label, actual, expected, failed)
      character(len=*), intent(in) :: label
      logical, intent(in) :: actual, expected
      logical, intent(inout) :: failed

      if (actual .neqv. expected) then
         write (*, '(a,1x,a)') 'FAIL', trim(label)
         failed = .true.
      end if
   end subroutine check_logical

   subroutine check_true(label, condition, failed)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      logical, intent(inout) :: failed

      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', trim(label)
         failed = .true.
      end if
   end subroutine check_true

end program test_chebyshev_spectral_positivity
