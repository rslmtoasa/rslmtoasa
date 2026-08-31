! XCF-07: unknown legacy selectors must fail during XC construction for every
! spin mode instead of falling through to the Barth-Hedin implementation.
program test_xc_unknown_selector_requires_fatal
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc
   implicit none

   type(control) :: ctl
   type(xc) :: functional

   call g_logger%init()
   call ctl%restore_to_default()
   ctl%nsp = 1
   ctl%txc = 10
   functional = xc(ctl)

   write (*, '(a)') 'UNEXPECTED: unknown legacy selector was accepted'
   error stop 1
end program test_xc_unknown_selector_requires_fatal
