! libXC paths: baseline, GGA radial, production contract, composition, PBE comparison, LDA reconciliation, and the capability fatal paths.
! Each sub-case runs in its own process (subcase_runner.f90); the test fails if any sub-case fails.
program test_xc_libxc
   use subcase_runner_mod, only: subcase_case_name, run_all_subcases
   implicit none

   character(len=31), parameter :: cases(10) = [character(len=31) :: &
      'libxc_xc_baseline', 'libxc_gga_radial', 'libxc_production_contract', 'xc_composition_semantics', 'pbe_gga_comparison', 'xc_lda_reconciliation', 'unsupported_family', 'capability_kinetic', 'capability_wrong_dimensionality', 'capability_missing_vxc']
   logical, parameter :: expect_fatal(10) = [.false., .false., .false., .false., .false., .false., .true., .true., .true., .true.]
   character(len=31), parameter :: workdirs(10) = [character(len=31) :: '', '', '', '', '', '', '', '', '', '']
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
      case ('libxc_xc_baseline')
         call test_libxc_xc_baseline()
      case ('libxc_gga_radial')
         call test_libxc_gga_radial()
      case ('libxc_production_contract')
         call test_libxc_production_contract()
      case ('xc_composition_semantics')
         call test_xc_composition_semantics()
      case ('pbe_gga_comparison')
         call test_pbe_gga_comparison()
      case ('xc_lda_reconciliation')
         call test_xc_lda_reconciliation()
      case ('unsupported_family')
         call test_xc_unsupported_family_requires_fatal()
      case ('capability_kinetic')
         call test_xc_capability_requires_fatal(1050)
      case ('capability_wrong_dimensionality')
         call test_xc_capability_requires_fatal(1019)
      case ('capability_missing_vxc')
         call test_xc_capability_requires_fatal(1160)
      case default
         write (*, '(a,a)') 'unknown sub-case ', name
         error stop 2
      end select
   end subroutine run_case

end program test_xc_libxc
