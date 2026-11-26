!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
! Catches in spin_response_rs: a spin-down offset other than spin_off (lmax = 3
! here, spin_off = 16), a d block taken from the wrong rows, a missing transpose
! in the trace, swapped spin roles of the two columns, and a missing conjugate or
! wrong sign in chi0_from_traces. Plumbing only: no physics is checked here.
!------------------------------------------------------------------------------
program test_spin_response_rs
   use precision_mod, only: rp
   use math_mod, only: pi
   use basis_mod, only: basis_init, nb
   use spin_response_rs_mod, only: chi0_traces, chi0_from_traces
   implicit none
   complex(rp), allocatable :: gij(:, :, :), gji(:, :, :)
   complex(rp) :: t(2, 2), c
   logical :: failed

   failed = .false.
   call basis_init(3)
   ! Every element outside the four d spin blocks, and inside them except the entries set below, would add 100 or more
   allocate (gij(nb, nb, 2), gji(nb, nb, 2))
   gij = (100.0_rp, 0.0_rp)
   gji = (100.0_rp, 0.0_rp)
   gij(5:9, 5:9, :) = 0.0_rp      ! up d
   gij(21:25, 21:25, :) = 0.0_rp  ! down d
   gji(5:9, 5:9, :) = 0.0_rp
   gji(21:25, 21:25, :) = 0.0_rp
   ! Column 1: Tr[G^upup_ji G^dndn_ij] = gji(5,6) gij(22,21) = 1 * 2; column 2: Tr[G^dndn_ji G^upup_ij] = gji(22,21) gij(5,6) = 3 * 4
   gji(5, 6, 1) = 1.0_rp
   gij(22, 21, 1) = 2.0_rp
   gji(22, 21, 1) = 3.0_rp
   gij(5, 6, 1) = 4.0_rp
   ! Second index: every set entry doubled, so both traces are four times larger
   gji(:, :, 2) = gji(:, :, 1)
   gij(:, :, 2) = gij(:, :, 1)
   gji(5, 6, 2) = 2.0_rp
   gij(22, 21, 2) = 4.0_rp
   gji(22, 21, 2) = 6.0_rp
   gij(5, 6, 2) = 8.0_rp

   t = chi0_traces(gij, gji)
   call require(all(abs(t(:, 1) - [2.0_rp, 8.0_rp]) < 1.0e-14_rp), 'column 1 = 2, 8')
   call require(all(abs(t(:, 2) - [12.0_rp, 48.0_rp]) < 1.0e-14_rp), 'column 2 = 12, 48')

   ! chi0 = -(x - conj(y))/(2 pi i): x = y = 1 + 2i gives -2/pi; x = 1, y = 0 gives i/(2 pi)
   c = chi0_from_traces((1.0_rp, 2.0_rp), (1.0_rp, 2.0_rp))
   call require(abs(c - cmplx(-2.0_rp/pi, 0.0_rp, rp)) < 1.0e-15_rp, 'x = y = 1 + 2i gives -2/pi')
   c = chi0_from_traces((1.0_rp, 0.0_rp), (0.0_rp, 0.0_rp))
   call require(abs(c - cmplx(0.0_rp, 0.5_rp/pi, rp)) < 1.0e-15_rp, 'x = 1, y = 0 gives i/(2 pi)')

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains
   subroutine require(ok, label)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label
      if (.not. ok) then
         failed = .true.
         write (*, '(2a)') 'FAIL: ', label
      end if
   end subroutine require
end program test_spin_response_rs
