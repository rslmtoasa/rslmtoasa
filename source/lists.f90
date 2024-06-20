module lists_mod
   use string_mod, only: sl
   implicit none

   private

   type list_characters
      character(len=sl), dimension(:), allocatable, private :: data
   contains
      procedure :: push => list_characters_push
      procedure :: pop => list_characters_pop
      ! procedure :: size => list_characters_size
      ! procedure :: sort => list_characters_sort
      ! procedure :: unique => list_characters_unique
      procedure :: copy => list_characters_copy
   end type list_characters

   type list_integers
      integer, dimension(:), allocatable, private :: data
   contains
      procedure :: push => list_integers_push
      ! procedure :: pop => list_integers_pop
      ! procedure :: size => list_integers_size
      ! procedure :: sort => list_integers_sort
      ! procedure :: unique => list_integers_unique
      procedure :: copy => list_integers_copy
   end type list_integers

   public :: list_characters
contains

! character list methods
   subroutine list_characters_push(this, new, index)
      class(list_characters), intent(inout) :: this
      character(len=*), intent(in) :: new
      character(len=sl), dimension(:), allocatable :: aux_data
      include 'list_abstract_methods/push.f90'
   end subroutine list_characters_push

   function list_characters_pop(this, index) result(val)
      class(list_characters), intent(inout) :: this
      character(len=sl) :: val
      include 'list_abstract_methods/pop.f90'
   end function list_characters_pop

   function list_characters_size(this) result(val)
      class(list_characters), intent(in) :: this
      integer :: val
      val = size(this%data)
   end function list_characters_size

   subroutine list_characters_sort(this)
      class(list_characters), intent(in) :: this
      character(len=sl), dimension(:), allocatable :: aux_data
      include 'list_abstract_methods/sort.f90'
   end subroutine list_characters_sort

   subroutine list_characters_unique(this)
      class(list_characters), intent(in) :: this
      character(len=sl), dimension(:), allocatable :: aux_data
      include 'list_abstract_methods/unique.f90'
   end subroutine list_characters_unique

   subroutine list_characters_copy(this, lst)
      class(list_characters), intent(inout) :: this
      character(len=sl), dimension(:), allocatable, intent(out) :: lst
      allocate (lst(size(this%data)))
      lst = this%data
   end subroutine list_characters_copy

! integer list methods
   subroutine list_integers_push(this, new, index)
      class(list_integers), intent(inout) :: this
      integer, intent(in) :: new
      integer, dimension(:), allocatable :: aux_data
      include 'list_abstract_methods/push.f90'
   end subroutine list_integers_push

   subroutine list_integers_copy(this, lst)
      class(list_integers), intent(inout) :: this
      integer, dimension(:), allocatable, intent(out) :: lst
      allocate (lst(size(this%data)))
      lst = this%data
   end subroutine list_integers_copy

! subroutine list_characters_size
!   class(list_characters), intent(in) :: this
! end subroutine list_characters_size

! subroutine list_characters_sort
!   class(list_characters), intent(inout) :: this

! end subroutine list_characters_sort

! subroutine list_characters_unique
!   class(list_characters), intent(inout) :: this
! end subroutine list_characters_unique

end module lists_mod
