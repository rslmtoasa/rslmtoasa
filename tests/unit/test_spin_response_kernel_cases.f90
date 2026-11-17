! Spin-response kernel oracles A1-A5 and R1 (no driver run needed). Argument after the case: oracle directory.
! Each sub-case runs in its own process (subcase_runner.f90); the test fails if any sub-case fails.
program test_spin_response_kernel_cases
   use subcase_runner_mod, only: subcase_case_name, run_all_subcases
   implicit none

   character(len=15), parameter :: cases(5) = [character(len=15) :: &
      'kernels', 'channels_kernel', 'a4', 'a5', 'r1']
   logical, parameter :: expect_fatal(5) = [.false., .false., .false., .false., .false.]
   character(len=1), parameter :: workdirs(5) = [character(len=1) :: '', '', '', '', '']
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
      case ('kernels')
         call test_spin_response_kernels()
      case ('channels_kernel')
         call test_spin_response_channels_kernel()
      case ('a4')
         call test_spin_response_a4()
      case ('a5')
         call test_spin_response_a5()
      case ('r1')
         call test_spin_response_r1()
      case default
         write (*, '(a,a)') 'unknown sub-case ', name
         error stop 2
      end select
   end subroutine run_case

end program test_spin_response_kernel_cases
