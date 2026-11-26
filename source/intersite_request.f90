!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: intersite_request_mod
!
!> @brief Group inter-site pair requests by seed site (CSR layout).
!> @details From ijpair(:,1:2) and the pair range owned by this rank, build the
!>          distinct seed sites (first atom i of each pair), the unique target
!>          sites j per seed, and a scatter map from each local pair slot to its
!>          (seed, target) entry.  Seeds and targets are kept in order of first
!>          appearance.  Reversed pairs (i,j) and (j,i) are separate seeds; no
!>          dedup across direction.  Local slot k corresponds to pair
!>          ij = start_pair + k - 1.  No MPI.
!------------------------------------------------------------------------------
module intersite_request_mod
   implicit none
   private

   type, public :: intersite_request
      integer :: nseed = 0                                !< number of distinct seed sites
      integer :: npair = 0                                !< number of local pair slots
      integer, allocatable :: seed_site(:)                !< (nseed) seed site id i
      integer, allocatable :: target_ptr(:)               !< (nseed+1) CSR offsets into target_site
      integer, allocatable :: target_site(:)              !< (ntarget) unique target site ids j
      integer, allocatable :: pair_seed(:)                !< (npair) seed index of local slot
      integer, allocatable :: pair_entry(:)               !< (npair) index into target_site of local slot
   contains
      procedure :: build
      procedure :: restore_to_default
      final :: destructor
   end type intersite_request

contains

   !> Build the grouping for pairs start_pair..end_pair of ijpair(:,1:2).
   subroutine build(this, ijpair, start_pair, end_pair)
      class(intersite_request), intent(inout) :: this
      integer, intent(in) :: ijpair(:, :)
      integer, intent(in) :: start_pair, end_pair
      integer :: k, ij, s, e, i, j, nt
      integer, allocatable :: tmp_site(:), count_s(:)

      call this%restore_to_default()
      this%npair = max(0, end_pair - start_pair + 1)
      allocate (this%pair_seed(this%npair), this%pair_entry(this%npair))
      allocate (tmp_site(this%npair), count_s(this%npair))
      count_s = 0

      ! Pass 1: seeds in order of first appearance, targets counted per seed.
      do k = 1, this%npair
         i = ijpair(start_pair + k - 1, 1)
         s = find_first(tmp_site, this%nseed, i)
         if (s == 0) then
            this%nseed = this%nseed + 1
            s = this%nseed
            tmp_site(s) = i
         end if
         this%pair_seed(k) = s
      end do
      allocate (this%seed_site(this%nseed), this%target_ptr(this%nseed + 1))
      this%seed_site = tmp_site(1:this%nseed)

      ! Pass 2: unique targets per seed, filled seed by seed (CSR).
      allocate (this%target_site(this%npair))
      nt = 0
      this%target_ptr(1) = 1
      do s = 1, this%nseed
         do k = 1, this%npair
            if (this%pair_seed(k) /= s) cycle
            ij = start_pair + k - 1
            j = ijpair(ij, 2)
            e = find_first(this%target_site(this%target_ptr(s):nt), nt - this%target_ptr(s) + 1, j)
            if (e == 0) then
               nt = nt + 1
               this%target_site(nt) = j
               e = nt - this%target_ptr(s) + 1
            end if
            this%pair_entry(k) = this%target_ptr(s) + e - 1
         end do
         this%target_ptr(s + 1) = nt + 1
      end do
      this%target_site = this%target_site(1:nt)
   end subroutine build

   !> Index of the first occurrence of val in list(1:n), or 0.
   pure function find_first(list, n, val) result(idx)
      integer, intent(in) :: list(:), n, val
      integer :: idx
      do idx = 1, n
         if (list(idx) == val) return
      end do
      idx = 0
   end function find_first

   subroutine restore_to_default(this)
      class(intersite_request), intent(inout) :: this
      this%nseed = 0
      this%npair = 0
      if (allocated(this%seed_site)) deallocate (this%seed_site)
      if (allocated(this%target_ptr)) deallocate (this%target_ptr)
      if (allocated(this%target_site)) deallocate (this%target_site)
      if (allocated(this%pair_seed)) deallocate (this%pair_seed)
      if (allocated(this%pair_entry)) deallocate (this%pair_entry)
   end subroutine restore_to_default

   subroutine destructor(this)
      type(intersite_request), intent(inout) :: this
      call this%restore_to_default()
   end subroutine destructor

end module intersite_request_mod
