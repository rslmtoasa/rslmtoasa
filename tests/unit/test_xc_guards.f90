! XC fatal paths: every sub-case must stop the run (non-zero exit, no UNEXPECTED).
! Each sub-case runs in its own process (subcase_runner.f90); the test fails if any sub-case fails.
program test_xc_guards
   use subcase_runner_mod, only: subcase_case_name, run_all_subcases
   implicit none

   character(len=32), parameter :: cases(9) = [character(len=32) :: &
      'unknown_selector', 'polarized_density_txc6', 'polarized_density_txc7', 'polarization_rejected_txc6_case3', 'polarization_rejected_txc6_case4', 'polarization_rejected_txc6_case5', 'polarization_rejected_txc7_case3', 'polarization_rejected_txc7_case4', 'polarization_rejected_txc7_case5']
   logical, parameter :: expect_fatal(9) = [.true., .true., .true., .true., .true., .true., .true., .true., .true.]
   character(len=32), parameter :: workdirs(9) = [character(len=32) :: '', '', '', '', '', '', '', '', '']
   character(len=64) :: which

   which = subcase_case_name()
   if (trim(which) == 'all') then
      call run_all_subcases(cases, expect_fatal, workdirs)
   else
      call run_case(trim(which))
   end if

contains

   subroutine run_case(name)
      character(len=*), intent(in) :: name
      select case (name)
      case ('unknown_selector')
         call test_xc_unknown_selector_requires_fatal()
      case ('polarized_density_txc6')
         call test_xc_legacy_polarized_requires_fatal(6)
      case ('polarized_density_txc7')
         call test_xc_legacy_polarized_requires_fatal(7)
      case ('polarization_rejected_txc6_case3')
         call test_xc_legacy_polarization_tolerance(6, 3)
      case ('polarization_rejected_txc6_case4')
         call test_xc_legacy_polarization_tolerance(6, 4)
      case ('polarization_rejected_txc6_case5')
         call test_xc_legacy_polarization_tolerance(6, 5)
      case ('polarization_rejected_txc7_case3')
         call test_xc_legacy_polarization_tolerance(7, 3)
      case ('polarization_rejected_txc7_case4')
         call test_xc_legacy_polarization_tolerance(7, 4)
      case ('polarization_rejected_txc7_case5')
         call test_xc_legacy_polarization_tolerance(7, 5)
      case default
         write (*, '(a,a)') 'unknown sub-case ', name
         error stop 2
      end select
   end subroutine run_case

end program test_xc_guards
