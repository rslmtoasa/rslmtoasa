!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
! PROGRAM: test_basis_dimensions
!
!> @brief Check the supported normal orbital-basis dimension contract.
!> @details The dimensions are derived from lmax by basis_init.  Keep these
!>          checks explicit so an accidental return to spd-only sizing cannot
!>          silently pass through the unit suite.
!------------------------------------------------------------------------------
program test_basis_dimensions
   use basis_mod, only: basis_init, lmax_basis, norb, spin_off, nb, bdg_mode
   implicit none

   logical :: failed

   failed = .false.

   call check_basis(1, 'sp', 4, 4, 8)
   call check_basis(2, 'spd', 9, 9, 18)
   call check_basis(3, 'spdf', 16, 16, 32)
   call test_reinitialization()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   subroutine check_basis(lmax, label, expected_norb, expected_spin_off, expected_nb)
      integer, intent(in) :: lmax, expected_norb, expected_spin_off, expected_nb
      character(len=*), intent(in) :: label

      call basis_init(lmax)
      call require(lmax_basis == lmax, trim(label)//' lmax')
      call require(norb == expected_norb, trim(label)//' norb')
      call require(spin_off == expected_spin_off, trim(label)//' spin offset')
      call require(nb == expected_nb, trim(label)//' normal spinor block dimension')
      call require(.not. bdg_mode, trim(label)//' normal-mode flag')
   end subroutine check_basis

   subroutine test_reinitialization()
      ! Repeated calls must replace all protected module state, not retain the
      ! dimensions from the first supported basis that was selected.
      call check_basis(1, 'reinit sp', 4, 4, 8)
      call check_basis(2, 'reinit spd', 9, 9, 18)
      call check_basis(1, 'reinit sp again', 4, 4, 8)
   end subroutine test_reinitialization

   subroutine require(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description

      if (.not. condition) then
         write (*, '(a,a)') 'FAIL: ', trim(description)
         failed = .true.
      end if
   end subroutine require

end program test_basis_dimensions
