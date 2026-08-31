! XCF-08: control%nsp is a global electronic-structure mode, not an XC
! channel count.  The atomic/XC radial representation is tested separately
! through its explicit two-channel contract.
program test_control_spin_semantics
   use control_mod, only: control
   use logger_mod, only: g_logger
   use self_mod, only: RADIAL_XC_SPIN_CHANNELS
   implicit none

   type(control) :: ctl
   logical :: failed
   integer :: mode

   call g_logger%init()
   call ctl%restore_to_default()
   failed = .false.

   call require(RADIAL_XC_SPIN_CHANNELS == 2, 'radial/XC representation has two local channels')
   do mode = 1, 4
      ctl%nsp = mode
      call require(ctl%is_collinear() .eqv. (mode == 1 .or. mode == 2), &
         'nsp='//int_string(mode)//' collinear query')
      call require(ctl%is_noncollinear() .eqv. (mode == 3 .or. mode == 4), &
         'nsp='//int_string(mode)//' noncollinear query')
      call require(ctl%has_soc() .eqv. (mode == 2 .or. mode == 4), &
         'nsp='//int_string(mode)//' SOC query')
      call require(ctl%uses_spinor_representation() .eqv. (mode >= 2), &
         'nsp='//int_string(mode)//' spinor representation query')
      call require(ctl%is_spin_polarized_mode(), &
         'nsp='//int_string(mode)//' permits a spin-polarized physical density')
   end do

   ! This is the specific ambiguity guard: the ordinary scalar-relativistic
   ! collinear mode still has two local radial/XC spin-density channels.
   ctl%nsp = 1
   call require(ctl%is_collinear() .and. .not. ctl%has_soc() .and. &
      RADIAL_XC_SPIN_CHANNELS == 2, &
      'nsp=1 remains collinear while the radial/XC path has two channels')

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

end program test_control_spin_semantics
