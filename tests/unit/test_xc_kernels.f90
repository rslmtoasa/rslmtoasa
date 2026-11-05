! Built-in (non-libXC) XC kernels: selectors, legacy LDA/PBE kernels, radial GGA, origin and residual controls.
! Each sub-case runs in its own process (subcase_runner.f90); the test fails if any sub-case fails.
program test_xc_kernels
   use subcase_runner_mod, only: subcase_case_name, run_all_subcases
   implicit none

   character(len=28), parameter :: cases(8) = [character(len=28) :: &
      'selector_semantics', 'lda_legacy_kernel', 'legacy_unpolarized_density', 'radial_gga', 'legacy_pbe_gga', 'vxc0sp_gga_origin', 'atomsc_residual_controls', 'legacy_polarization_accepted']
   logical, parameter :: expect_fatal(8) = [.false., .false., .false., .false., .false., .false., .false., .false.]
   character(len=28), parameter :: workdirs(8) = [character(len=28) :: '', '', '', '', '', '', '', '']
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
      case ('selector_semantics')
         call test_xc_selector_semantics()
      case ('lda_legacy_kernel')
         call test_xc_lda_legacy_kernel()
      case ('legacy_unpolarized_density')
         call test_xc_legacy_unpolarized_density()
      case ('radial_gga')
         call test_radial_gga()
      case ('legacy_pbe_gga')
         call test_legacy_pbe_gga()
      case ('vxc0sp_gga_origin')
         call test_vxc0sp_gga_origin()
      case ('atomsc_residual_controls')
         call test_atomsc_residual_controls()
      case ('legacy_polarization_accepted')
         call test_xc_legacy_polarization_tolerance(0, -1)
      case default
         write (*, '(a,a)') 'unknown sub-case ', name
         error stop 2
      end select
   end subroutine run_case

end program test_xc_kernels
