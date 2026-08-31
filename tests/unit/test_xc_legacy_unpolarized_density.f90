! XCF-08: TXC=6 and TXC=7 are unpolarized-only, but an equal two-channel
! density is valid even in the spin-capable global mode nsp=1.
program test_xc_legacy_unpolarized_density
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use precision_mod, only: rp
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc
   implicit none

   type(control) :: ctl
   type(xc) :: functional
   real(rp), parameter :: zero(2) = 0.0_rp
   real(rp) :: energy, v_down, v_up
   integer :: selector
   logical :: failed

   call g_logger%init()
   call ctl%restore_to_default()
   ctl%nsp = 1
   failed = .false.

   do selector = 6, 7
      ctl%txc = selector
      functional = xc(ctl)
      call functional%XCPOT(0.23_rp, 0.23_rp, 0.46_rp, zero, zero, 1.0_rp, &
         v_down, v_up, energy)
      call require(ieee_is_finite(energy) .and. ieee_is_finite(v_down) .and. ieee_is_finite(v_up), &
         'TXC='//int_string(selector)//' accepts an equal two-channel density')
      call require(v_down == v_up, 'TXC='//int_string(selector)//' retains its unpolarized potential')
   end do

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   function int_string(value) result(text)
      integer, intent(in) :: value
      character(len=16) :: text
      write (text, '(i0)') value
   end function int_string

   subroutine require(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description
      if (.not. condition) then
         write (*, '(a)') 'FAIL: '//trim(description)
         failed = .true.
      end if
   end subroutine require

end program test_xc_legacy_unpolarized_density
