! LDA+U local-moment kernels: collinear oracle, spinor and spin-frame kernel, host-site swap. Arguments after the case: collinear fixture, site-swap fixture.
! Each sub-case runs in its own process (subcase_runner.f90); the test fails if any sub-case fails.
program test_ldau_local_moments
   use subcase_runner_mod, only: subcase_case_name, run_all_subcases
   implicit none

   character(len=19), parameter :: cases(3) = [character(len=19) :: &
      'collinear_kernel', 'spinor_frame_kernel', 'site_swap']
   logical, parameter :: expect_fatal(3) = [.false., .false., .false.]
   character(len=1), parameter :: workdirs(3) = [character(len=1) :: '', '', '']
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
      case ('collinear_kernel')
         call test_local_moments_kernel()
      case ('spinor_frame_kernel')
         call test_hubbard_spinor_moments()
      case ('site_swap')
         call test_local_moments_site_swap()
      case default
         write (*, '(a,a)') 'unknown sub-case ', name
         error stop 2
      end select
   end subroutine run_case

end program test_ldau_local_moments
