! XCF-08 negative test driver.  The selector is supplied as a command-line
! argument so CTest can exercise several independent capability failures.
subroutine test_xc_capability_requires_fatal(selector)
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc
   implicit none

   type(control) :: ctl
   type(xc) :: functional
   integer, intent(in) :: selector

   call g_logger%init()

   call ctl%restore_to_default()
   ctl%nsp = 1
   ctl%txc = selector
   functional = xc(ctl)

   write (*, '(a,i0)') 'UNEXPECTED: incompatible TXC selector was accepted: ', selector
   error stop 1
end subroutine test_xc_capability_requires_fatal
