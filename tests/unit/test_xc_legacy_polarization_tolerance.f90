! XCF-09: the TXC=6/7 density guard is relative to the supplied total
! density. Exercise both sides of the documented 1e-10 tolerance.
program test_xc_legacy_polarization_tolerance
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use precision_mod, only: rp
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc, LEGACY_UNPOLARIZED_REL_TOL
   implicit none

   integer :: selector, case_id, ios
   character(len=32) :: selector_text, case_text
   logical :: failed

   call g_logger%init()
   failed = .false.

   call get_command_argument(1, selector_text)
   call get_command_argument(2, case_text)
   if (len_trim(selector_text) == 0 .and. len_trim(case_text) == 0) then
      ! The default run covers all accepted density cases.
      do selector = 6, 7
         do case_id = 0, 2
            call check_accepted(selector, case_id)
         end do
         do case_id = 6, 8
            call check_accepted(selector, case_id)
         end do
      end do
   else
      read (selector_text, *, iostat=ios) selector
      if (ios /= 0) error stop 'missing or invalid TXC selector argument'
      read (case_text, *, iostat=ios) case_id
      if (ios /= 0) error stop 'missing or invalid tolerance case argument'
      if (case_id <= 2 .or. (case_id >= 6 .and. case_id <= 8)) then
         call check_accepted(selector, case_id)
      else
         call check_rejected(selector, case_id)
      end if
   end if

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine density_for_case(case_id, rho_down, rho_up, description)
      integer, intent(in) :: case_id
      real(rp), intent(out) :: rho_down, rho_up
      character(len=*), intent(out) :: description
      real(rp), parameter :: total = 0.40_rp
      real(rp) :: difference

      select case (case_id)
      case (0)
         difference = 0.0_rp
         description = 'exact equality'
      case (1)
         difference = 0.1_rp*LEGACY_UNPOLARIZED_REL_TOL*total
         description = '0.1 epsilon_spin'
      case (2)
         difference = 0.9_rp*LEGACY_UNPOLARIZED_REL_TOL*total
         description = '0.9 epsilon_spin'
      case (3)
         difference = 1.1_rp*LEGACY_UNPOLARIZED_REL_TOL*total
         description = '1.1 epsilon_spin'
      case (4)
         difference = 10.0_rp*LEGACY_UNPOLARIZED_REL_TOL*total
         description = '10 epsilon_spin'
      case (5)
         difference = 0.10_rp
         description = 'clearly polarized density'
      case (6)
         rho_down = 0.0_rp
         rho_up = 0.0_rp
         description = 'zero density'
         return
      case (7)
         rho_down = 1.0e-25_rp
         rho_up = 0.0_rp
         description = 'down-channel asymmetric tail'
         return
      case (8)
         rho_down = 0.0_rp
         rho_up = 1.0e-25_rp
         description = 'up-channel asymmetric tail'
         return
      case default
         error stop 'unknown tolerance case'
      end select

      ! Keep RHO equal to the channel sum while varying only the local spin
      ! asymmetry. XCPOT receives down in RHO1 and up in RHO2.
      rho_down = 0.5_rp*total + 0.5_rp*difference
      rho_up = 0.5_rp*total - 0.5_rp*difference
   end subroutine density_for_case

   subroutine check_accepted(selector, case_id)
      integer, intent(in) :: selector, case_id
      type(control) :: ctl
      type(xc) :: functional
      real(rp) :: rho_down, rho_up, energy, v_down, v_up
      character(len=64) :: description

      call density_for_case(case_id, rho_down, rho_up, description)
      call ctl%restore_to_default()
      ctl%nsp = 1
      ctl%txc = selector
      functional = xc(ctl)
      call functional%XCPOT(rho_down, rho_up, rho_down + rho_up, [0.0_rp, 0.0_rp], &
         [0.0_rp, 0.0_rp], 1.0_rp, v_down, v_up, energy)
      call require(ieee_is_finite(energy) .and. ieee_is_finite(v_down) .and. ieee_is_finite(v_up), &
         'TXC='//int_string(selector)//' accepts '//trim(description))
      if (case_id >= 6) then
         call require(v_down == 0.0_rp .and. v_up == 0.0_rp .and. energy == 0.0_rp, &
            'TXC='//int_string(selector)//' returns zero XC for '//trim(description))
      end if
   end subroutine check_accepted

   subroutine check_rejected(selector, case_id)
      integer, intent(in) :: selector, case_id
      type(control) :: ctl
      type(xc) :: functional
      real(rp) :: rho_down, rho_up, energy, v_down, v_up
      character(len=64) :: description

      call density_for_case(case_id, rho_down, rho_up, description)
      call ctl%restore_to_default()
      ctl%nsp = 1
      ctl%txc = selector
      functional = xc(ctl)
      call functional%XCPOT(rho_down, rho_up, rho_down + rho_up, [0.0_rp, 0.0_rp], &
         [0.0_rp, 0.0_rp], 1.0_rp, v_down, v_up, energy)
      write (*, '(a)') 'UNEXPECTED: TXC='//int_string(selector)//' accepted '//trim(description)
      error stop 1
   end subroutine check_rejected

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

end program test_xc_legacy_polarization_tolerance
