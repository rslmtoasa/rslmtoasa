!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
! Catches: wrong CSR offsets / scatter indices in the seed-target grouping, or
! lost self, duplicate, reversed pairs and pair order.
!------------------------------------------------------------------------------
program test_intersite_request
   use intersite_request_mod, only: intersite_request
   implicit none
   type(intersite_request) :: req
   integer :: ijpair(9, 2)
   logical :: failed

   failed = .false.
   ! 1:(3,3) self  2:(3,5)  3:(1,2)  4:(3,5) dup of 2  5:(2,1) reversed of 3
   ! 6:(3,7) same seed as 2  7:(1,2) dup of 3  8:(5,5)  9:(4,4)
   ijpair(:, 1) = [3, 3, 1, 3, 2, 3, 1, 5, 4]
   ijpair(:, 2) = [3, 5, 2, 5, 1, 7, 2, 5, 4]

   ! Partial range 2..7 (start > 1, end < njij)
   call req%build(ijpair, 2, 7)
   ! seeds by first appearance: 3, 1, 2
   ! slots: (3,5) (1,2) (3,5) (2,1) (3,7) (1,2)
   call require(req%npair == 6, 'npair')
   call require(req%nseed == 3, 'nseed')
   call require(all(req%seed_site == [3, 1, 2]), 'seed_site')
   call require(all(req%target_ptr == [1, 3, 4, 5]), 'target_ptr')
   call require(all(req%target_site == [5, 7, 2, 1]), 'target_site')
   call require(all(req%pair_seed == [1, 2, 1, 3, 1, 2]), 'pair_seed')
   call require(all(req%pair_entry == [1, 3, 1, 4, 2, 3]), 'pair_entry')

   ! Full range keeps the self-pair as a target of its own seed
   call req%build(ijpair, 1, 9)
   call require(all(req%seed_site == [3, 1, 2, 5, 4]), 'full seed_site')
   call require(all(req%target_ptr == [1, 4, 5, 6, 7, 8]), 'full target_ptr')
   call require(all(req%target_site == [3, 5, 7, 2, 1, 5, 4]), 'full target_site')
   call require(all(req%pair_entry == [1, 2, 4, 2, 5, 3, 4, 6, 7]), 'full pair_entry')

   call req%restore_to_default()
   call require(req%nseed == 0 .and. .not. allocated(req%seed_site), 'restore_to_default')

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
end program test_intersite_request
