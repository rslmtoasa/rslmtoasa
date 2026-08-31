! XCF-08 negative test driver.  The selector is supplied as a command-line
! argument so CTest can exercise several independent capability failures.
program test_xc_capability_requires_fatal
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc
   implicit none

   type(control) :: ctl
   type(xc) :: functional
   character(len=32) :: selector_text
   integer :: selector, ios

   call g_logger%init()
   call get_command_argument(1, selector_text)
   read (selector_text, *, iostat=ios) selector
   if (ios /= 0) error stop 'missing or invalid TXC selector argument'

   call ctl%restore_to_default()
   ctl%nsp = 1
   ctl%txc = selector
   functional = xc(ctl)

   write (*, '(a,i0)') 'UNEXPECTED: incompatible TXC selector was accepted: ', selector
   error stop 1
end program test_xc_capability_requires_fatal
