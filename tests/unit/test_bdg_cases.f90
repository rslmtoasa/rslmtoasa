! BdG kernel oracles O1, O1b, O2, O3, O4 (no driver run needed). Argument after the case: oracle directory.
! Each sub-case runs in its own process (subcase_runner.f90); the test fails if any sub-case fails.
program test_bdg_cases
   use subcase_runner_mod, only: subcase_case_name, run_all_subcases
   use test_bdg_mod
   implicit none

   character(len=4), parameter :: cases(5) = [character(len=4) :: 'o1', 'o1b', 'o2', 'o3', 'o4']
   logical, parameter :: expect_fatal(5) = [.false., .false., .false., .false., .false.]
   character(len=1), parameter :: workdirs(5) = [character(len=1) :: '', '', '', '', '']
   character(len=64) :: which
   character(len=512) :: dir

   which = subcase_case_name()
   if (trim(which) == 'all') then
      call run_all_subcases(cases, expect_fatal, workdirs)
   else
      call get_command_argument(2, dir)
      if (len_trim(dir) == 0) error stop 'usage: UnitBdg <case> <oracle directory>'
      select case (trim(which))
      case ('o1')
         call test_bdg_o1(trim(dir))
      case ('o1b')
         call test_bdg_o1b(trim(dir))
      case ('o2')
         call test_bdg_o2(trim(dir))
      case ('o3')
         call test_bdg_o3(trim(dir))
      case ('o4')
         call test_bdg_o4(trim(dir))
      case default
         write (*, '(a,a)') 'unknown sub-case ', trim(which)
         error stop 2
      end select
   end if
end program test_bdg_cases
