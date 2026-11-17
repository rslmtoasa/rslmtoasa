! Spin-response checks on driver output (bcc Fe scratch, plus the AF scratch for s2af). Argument after the case: oracle directory.
! Each sub-case runs in its own process (subcase_runner.f90); the test fails if any sub-case fails.
program test_spin_response_driver_cases
   use subcase_runner_mod, only: subcase_case_name, run_all_subcases
   implicit none

   character(len=8), parameter :: cases(7) = [character(len=8) :: &
      'c1', 'c2', 'c3', 's2', 'channels', 'radial', 's2af']
   logical, parameter :: expect_fatal(7) = [.false., .false., .false., .false., .false., .false., .false.]
   character(len=5), parameter :: workdirs(7) = [character(len=5) :: '', '', '', '', '', '', '../af']
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
      case ('c1')
         call test_spin_response_c1()
      case ('c2')
         call test_spin_response_c2()
      case ('c3')
         call test_spin_response_c3()
      case ('s2')
         call test_spin_response_s2()
      case ('channels')
         call test_spin_response_channels()
      case ('radial')
         call test_spin_response_radial()
      case ('s2af')
         call test_spin_response_s2af()
      case default
         write (*, '(a,a)') 'unknown sub-case ', name
         error stop 2
      end select
   end subroutine run_case

end program test_spin_response_driver_cases
