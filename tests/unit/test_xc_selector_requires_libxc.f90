! This test is registered only in a build without libXC and must terminate
! during XC construction rather than falling through to legacy XCPOT.
program test_xc_selector_requires_libxc
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc
   implicit none

   type(control) :: ctl
   type(xc) :: functional

   call g_logger%init()
   call ctl%restore_to_default()
   ctl%nsp = 1
   ctl%txc = 1001
   functional = xc(ctl)

   write (*, '(a)') 'UNEXPECTED: TXC=1001 returned without libXC'
   error stop 1
end program test_xc_selector_requires_libxc
