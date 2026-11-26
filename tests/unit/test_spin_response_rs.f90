!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
! Catches in spin_response_rs: a spin-down offset other than spin_off (lmax = 3
! here, spin_off = 16), a d block taken from the wrong rows, a missing transpose
! in the trace, swapped spin roles of the two columns, and a missing conjugate or
! wrong sign in chi0_from_traces; in d_moment_trace a wrong orbital range or spin
! sign; in type_sums a seed taken from a non-self pair or from a later self pair,
! pairs of another seed summed in, grouping by the wrong atom, or the imaginary
! part included. Plumbing only: no physics is checked here.
!------------------------------------------------------------------------------
program test_spin_response_rs
   use precision_mod, only: rp
   use math_mod, only: pi
   use basis_mod, only: basis_init, nb
   use spin_response_rs_mod, only: chi0_traces, chi0_from_traces, d_moment_trace, type_sums
   implicit none
   complex(rp), allocatable :: gij(:, :, :), gji(:, :, :), dm(:)
   complex(rp) :: t(2, 2), c
   integer :: mu, ijpair(8, 2), seed(2), npair(2)
   real(rp) :: m(2), chi0s(2, 2)
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

   ! d_moment_trace: d up diagonal 1..5, d down diagonal 0.5 each, every other element 100: 15 - 2.5 = 12.5
   gij(:, :, 1) = (100.0_rp, 0.0_rp)
   do mu = 5, 9
      gij(mu, mu, 1) = real(mu - 4, rp)
      gij(mu + 16, mu + 16, 1) = 0.5_rp
   end do
   dm = d_moment_trace(gij(:, :, 1:1))
   call require(abs(dm(1) - 12.5_rp) < 1.0e-14_rp, 'd_moment_trace = 12.5')

   ! type_sums: atoms 1..6 of types 1, 2, 1, 2, 1, 2. Pair 1 (3,4) is not a self pair, pair 5 (3,3) a later self pair of
   ! type 1: neither may become a seed. Seeds 1 (m 2.0) and 2 (m -1.5). Type 1 rows: pairs 2, 3, 4, 8 (j types 1, 2, 1, 1),
   ! chi0s(1,1) = -5 - 0.5 - 0.125 = -5.625, chi0s(1,2) = -1; type 2 rows: pairs 6, 7, chi0s(2,2) = -4, chi0s(2,1) = -0.25.
   ijpair(:, 1) = [3, 1, 1, 1, 3, 2, 2, 1]
   ijpair(:, 2) = [4, 1, 2, 3, 3, 2, 1, 5]
   call type_sums(ijpair, [1, 2, 1, 2, 1, 2], 2, &
                  [(10.0_rp, 0.0_rp), (-5.0_rp, 0.0_rp), (-1.0_rp, 0.0_rp), (-0.5_rp, 0.0_rp), (-7.0_rp, 0.0_rp), &
                   (-4.0_rp, 0.0_rp), (-0.25_rp, 0.0_rp), (-0.125_rp, 3.0_rp)], &
                  [0.0_rp, 2.0_rp, 0.0_rp, 0.0_rp, 9.0_rp, -1.5_rp, 0.0_rp, 0.0_rp], seed, npair, m, chi0s)
   call require(all(seed == [1, 2]), 'type_sums seed = 1, 2')
   call require(all(npair == [4, 2]), 'type_sums npair = 4, 2')
   call require(all(abs(m - [2.0_rp, -1.5_rp]) < 1.0e-14_rp), 'type_sums m = 2.0, -1.5')
   call require(all(abs(chi0s - reshape([-5.625_rp, -0.25_rp, -1.0_rp, -4.0_rp], [2, 2])) < 1.0e-14_rp), &
                'type_sums chi0s = [-5.625 -1; -0.25 -4]')

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
