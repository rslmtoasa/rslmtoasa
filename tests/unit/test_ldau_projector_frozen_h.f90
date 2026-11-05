! LDA+U projector frozen-H oracles: common and multisite.
! Each sub-case runs in its own process (subcase_runner.f90); the test fails if any sub-case fails.
program test_ldau_projector_frozen_h
   use subcase_runner_mod, only: subcase_case_name, run_all_subcases
   implicit none

   character(len=9), parameter :: cases(2) = [character(len=9) :: &
      'common', 'multisite']
   logical, parameter :: expect_fatal(2) = [.false., .false.]
   character(len=1), parameter :: workdirs(2) = [character(len=1) :: '', '']
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
      case ('common')
         call projector_common_frozen_h()
      case ('multisite')
         call projector_multisite_frozen_h()
      case default
         write (*, '(a,a)') 'unknown sub-case ', name
         error stop 2
      end select
   end subroutine run_case

end program test_ldau_projector_frozen_h
