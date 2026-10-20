program test_scf_b1_guard
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use precision_mod, only: rp
   use logger_mod, only: g_logger
   use scf_guard_mod, only: guard_q0_capacity, guard_nonfinite_scalar
   implicit none

   character(len=32) :: mode

   call g_logger%init()
   call get_command_argument(1, mode)

   select case (trim(mode))
   case ('capacity')
      ! A d-shell has five states per spin; test the approved upper-bound path.
      call guard_q0_capacity('Fe', 3, 2, 1, 17, 5.0_rp + 2.0e-6_rp)
      error stop 'over-capacity q0 was not rejected'
   case ('nonfinite')
      call guard_nonfinite_scalar('Co', 2, 1, 2, 'new q1', 9, &
         ieee_value(0.0_rp, ieee_quiet_nan))
      error stop 'non-finite q component was not rejected'
   case ('boundary')
      ! The stated allowance itself remains accepted; no extra lower bound.
      call guard_q0_capacity('Fe', 3, 2, 1, 17, 5.0_rp + 1.0e-6_rp)
      call guard_nonfinite_scalar('Co', 2, 1, 2, 'new q1', 9, 0.0_rp)
      write (*, '(a)') 'B1 guard boundary accepted'
   case default
      error stop 'expected test mode: capacity, nonfinite, or boundary'
   end select
end program test_scf_b1_guard
