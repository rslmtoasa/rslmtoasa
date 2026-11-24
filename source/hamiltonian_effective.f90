! Sparse row accumulation shared by periodic and finite-graph HOH export.
module hamiltonian_effective_mod
   use iso_fortran_env, only: int64
   use ieee_arithmetic, only: ieee_is_finite
   use precision_mod, only: rp
   implicit none
   private
   public :: effective_row, resolve_translation

   type :: effective_row
      integer :: count = 0
      integer, allocatable :: key(:, :), bucket(:)
      complex(rp), allocatable :: block(:, :, :)
   contains
      procedure :: reset
      procedure :: lookup
      procedure :: add
      procedure :: order
   end type effective_row

contains

   ! A key contains the column site followed by its three image indices.
   subroutine reset(this, n)
      class(effective_row), intent(inout) :: this
      integer, intent(in) :: n
      if (.not. allocated(this%bucket)) then
         allocate(this%key(4, 16), this%block(n, n, 16), this%bucket(32))
      end if
      this%count = 0
      this%bucket = 0
   end subroutine reset

   integer function hash_slot(key, capacity) result(slot)
      integer, intent(in) :: key(4), capacity
      integer(int64) :: value
      value = modulo(sum(int(key, int64)*[73856093_int64, 19349663_int64, 83492791_int64, 2654435761_int64]), &
                     int(capacity, int64))
      slot = int(value) + 1
   end function hash_slot

   integer function probe(this, key) result(slot)
      class(effective_row), intent(in) :: this
      integer, intent(in) :: key(4)
      integer :: record
      slot = hash_slot(key, size(this%bucket))
      do
         record = this%bucket(slot)
         if (record == 0) return
         if (all(this%key(:, record) == key)) return
         slot = modulo(slot, size(this%bucket)) + 1
      end do
   end function probe

   integer function lookup(this, key) result(record)
      class(effective_row), intent(in) :: this
      integer, intent(in) :: key(4)
      record = this%bucket(probe(this, key))
   end function lookup

   subroutine add(this, key, value)
      class(effective_row), intent(inout) :: this
      integer, intent(in) :: key(4)
      complex(rp), intent(in) :: value(:, :)
      integer :: slot, record, n, capacity
      integer, allocatable :: keys(:, :)
      complex(rp), allocatable :: blocks(:, :, :)

      slot = probe(this, key)
      record = this%bucket(slot)
      if (record /= 0) then
         this%block(:, :, record) = this%block(:, :, record) + value
         return
      end if
      if (this%count == size(this%key, 2)) then
         if (size(this%bucket) > huge(capacity)/2) error stop 'Effective Hamiltonian row exceeds index capacity'
         capacity = size(this%bucket)
         n = size(this%block, 1)
         allocate(keys(4, capacity), blocks(n, n, capacity))
         keys(:, :this%count) = this%key
         blocks(:, :, :this%count) = this%block
         call move_alloc(keys, this%key)
         call move_alloc(blocks, this%block)
         deallocate(this%bucket)
         allocate(this%bucket(2*capacity))
         this%bucket = 0
         do record = 1, this%count
            this%bucket(probe(this, this%key(:, record))) = record
         end do
         slot = probe(this, key)
      end if
      this%count = this%count + 1
      this%bucket(slot) = this%count
      this%key(:, this%count) = key
      this%block(:, :, this%count) = value
   end subroutine add

   function order(this) result(indices)
      class(effective_row), intent(in) :: this
      integer :: indices(this%count)
      integer :: i
      indices = [(i, i=1,this%count)]
      call sort(1, this%count)
   contains
      recursive subroutine sort(first, last)
         integer, intent(in) :: first, last
         integer :: left, right, pivot(4), temp
         if (first >= last) return
         pivot = this%key(:, indices(first + (last-first)/2))
         left = first
         right = last
         do
            do while (before(this%key(:, indices(left)), pivot))
               left = left + 1
            end do
            do while (before(pivot, this%key(:, indices(right))))
               right = right - 1
            end do
            if (left > right) exit
            temp = indices(left)
            indices(left) = indices(right)
            indices(right) = temp
            left = left + 1
            right = right - 1
            if (left > right) exit
         end do
         call sort(first, right)
         call sort(left, last)
      end subroutine sort
   end function order

   pure logical function before(left, right)
      integer, intent(in) :: left(4), right(4)
      integer :: i
      before = .false.
      do i = 1, 4
         if (left(i) == right(i)) cycle
         before = left(i) < right(i)
         return
      end do
   end function before

   ! The residual is measured in the same Cartesian units as cell and delta.
   subroutine resolve_translation(cell, inverse, delta, tolerance, translation, found)
      real(rp), intent(in) :: cell(3,3), inverse(3,3), delta(3), tolerance
      integer, intent(out) :: translation(3)
      logical, intent(out) :: found
      real(rp) :: fractional(3)
      translation = 0
      found = .false.
      fractional = matmul(inverse, delta)
      if (.not. all(ieee_is_finite(fractional))) return
      ! Reserve enough integer range to add two independently mapped hops.
      if (any(abs(fractional) >= real(huge(0)/2, rp))) return
      translation = nint(fractional)
      found = norm2(delta - matmul(cell, real(translation, rp))) <= tolerance
   end subroutine resolve_translation
end module hamiltonian_effective_mod
