!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
! Catches: self pair not first (self), neighbour order changed (order), zero
! padding of nn kept as a pair (zeros), repeated neighbour kept (dup), wrong
! recursion-site range or irec indirection (range).
!------------------------------------------------------------------------------
program test_nonlocal_pairs
   use lattice_mod, only: lattice
   use nonlocal_density_mod, only: nonlocal_density
   implicit none
   type(lattice) :: lat
   type(nonlocal_density) :: nld
   logical :: failed

   failed = .false.
   ! Cluster atoms 1..5. Recursion sites 1,2 are cluster atoms 4 and 2.
   allocate (lat%irec(2), lat%nn(5, 6))
   lat%nn = 0
   lat%irec = [4, 2]
   lat%nn(4, :) = [5, 3, 0, 1, 3, 0]   ! four slots: 3, zero, 1, 3 again
   lat%nn(2, :) = [3, 5, 4, 0, 0, 0]

   call nld%neighbour_pairs(lat, 1, 2)
   call require(size(nld%pairs, 1) == 6, 'count')
   call require(all(nld%pairs(1, :) == [4, 4]), 'self')
   call require(all(nld%pairs(2:3, 2) == [3, 1]), 'order')
   call require(all(nld%pairs(1:3, 1) == 4), 'zeros')
   call require(all(nld%pairs(4, :) == [2, 2]) .and. all(nld%pairs(5:6, 2) == [5, 4]), 'range')

   call nld%neighbour_pairs(lat, 2, 2)
   call require(size(nld%pairs, 1) == 3 .and. all(nld%pairs(:, 1) == 2), 'dup')

   call nld%restore_to_default()
   call require(.not. allocated(nld%pairs), 'restore_to_default')

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
end program test_nonlocal_pairs
