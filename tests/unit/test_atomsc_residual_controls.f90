! XCF-05: the physical integrated residual is observable but cannot drive
! historical ATOMSC solver, mixing, or termination decisions.
program test_atomsc_residual_controls
   use precision_mod, only: rp
   use self_mod, only: atomsc_solver_tolerance, atomsc_mixing_beta, atomsc_should_stop
   implicit none

   real(rp), parameter :: tolrsq = 1.0e-8_rp, beta = 0.3_rp, tol = 1.0e-6_rp
   real(rp), parameter :: weights(4) = [1.0_rp/3.0_rp, 4.0_rp/3.0_rp, 2.0_rp/3.0_rp, 1.0_rp/3.0_rp]
   real(rp), parameter :: drdi(4) = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]
   real(rp), parameter :: density_delta(4) = [1.0_rp, -0.5_rp, 0.25_rp, -0.1_rp]
   real(rp) :: drho_control, drho_integrated, solver_lo, solver_hi, beta_lo, beta_hi
   logical :: stop_lo, stop_hi, failed

   failed = .false.

   ! These are the two production accumulators for one deterministic update.
   drho_control = sum(weights*abs(density_delta))
   drho_integrated = sum(weights*drdi*abs(density_delta))
   call require(drho_integrated > drho_control, 'fixture distinguishes integrated and control residuals')

   ! Every result below is held at the same historical control value while
   ! the physical observable is changed substantially.
   solver_lo = atomsc_solver_tolerance(2, drho_control, tolrsq)
   beta_lo = atomsc_mixing_beta(2, drho_control, beta)
   stop_lo = atomsc_should_stop(2, 80, drho_control, tol)
   drho_integrated = 100.0_rp*drho_control
   solver_hi = atomsc_solver_tolerance(2, drho_control, tolrsq)
   beta_hi = atomsc_mixing_beta(2, drho_control, beta)
   stop_hi = atomsc_should_stop(2, 80, drho_control, tol)
   call require_close(solver_lo, solver_hi, 'solver tolerance independent of integrated residual')
   call require_close(beta_lo, beta_hi, 'BETA1 independent of integrated residual')
   call require(stop_lo .eqv. stop_hi, 'termination independent of integrated residual')

   ! Pin all three strict historical thresholds: 2, 1, and 1e-6.
   call require_close(atomsc_solver_tolerance(2, 2.0_rp, tolrsq), tolrsq, 'solver threshold is strict at DRHO=2')
   call require_close(atomsc_solver_tolerance(2, 2.0_rp + 1.0e-12_rp, tolrsq), 1.0e-3_rp, &
                      'solver threshold switches above DRHO=2')
   call require_close(atomsc_mixing_beta(2, 1.0_rp, beta), beta, 'BETA1 threshold is strict at DRHO=1')
   call require_close(atomsc_mixing_beta(2, 1.0_rp - 1.0e-12_rp, beta), 0.5_rp, &
                      'BETA1 threshold switches below DRHO=1')
   call require(.not. atomsc_should_stop(2, 80, tol, tol), 'termination threshold is strict at DRHO=1e-6')
   call require(atomsc_should_stop(2, 80, tol - 1.0e-12_rp, tol), 'termination threshold stops below DRHO=1e-6')
   call require(atomsc_should_stop(79, 80, 100.0_rp, tol), 'historical final-iteration stop is retained')

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine require_close(actual, expected, description)
      real(rp), intent(in) :: actual, expected
      character(len=*), intent(in) :: description
      if (abs(actual - expected) > 1.0e-14_rp*max(1.0_rp, abs(expected))) then
         write (*, '(a,2(1x,es18.10))') 'FAIL: '//trim(description)//' actual expected', actual, expected
         failed = .true.
      end if
   end subroutine require_close

   subroutine require(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description
      if (.not. condition) then
         write (*, '(a)') 'FAIL: '//trim(description)
         failed = .true.
      end if
   end subroutine require

end program test_atomsc_residual_controls
