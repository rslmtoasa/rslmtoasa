! Negative XCF-06 test: native MGGA SCAN must fail during XC construction.
program test_xc_unsupported_family_requires_fatal
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc
   implicit none

   type(control) :: ctl
   type(xc) :: functional

   call g_logger%init()
   call ctl%restore_to_default()
   ctl%nsp = 2
   ctl%txc = 1263  ! TXC=1000+XC_MGGA_X_SCAN (native ID 263)
   functional = xc(ctl)

   write (*, '(a)') 'UNEXPECTED: unsupported MGGA functional was accepted'
   error stop 1
end program test_xc_unsupported_family_requires_fatal
