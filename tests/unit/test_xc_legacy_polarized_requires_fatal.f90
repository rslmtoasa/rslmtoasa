! XCF-08 negative test driver for the unpolarized-only legacy kernels.
subroutine test_xc_legacy_polarized_requires_fatal(selector)
   use precision_mod, only: rp
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc
   implicit none

   type(control) :: ctl
   type(xc) :: functional
   real(rp), parameter :: zero(2) = 0.0_rp
   real(rp) :: energy, v_down, v_up
   integer, intent(in) :: selector

   call g_logger%init()
   call ctl%restore_to_default()
   ctl%nsp = 1
   ctl%txc = selector
   functional = xc(ctl)
   call functional%XCPOT(0.10_rp, 0.30_rp, 0.40_rp, zero, zero, 1.0_rp, &
      v_down, v_up, energy)

   write (*, '(a,i0)') 'UNEXPECTED: polarized density accepted by TXC=', selector
   error stop 1
end subroutine test_xc_legacy_polarized_requires_fatal
