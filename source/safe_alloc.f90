!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Safe_Alloc
!
!> @author
!> Angela Klautau
!> Ramon Cardias
!> Lucas P. Campagna
!> S. Frota-Pessôa
!> Pascoal R. Peduto
!> Anders Bergman
!> S. B. Legoas
!> H. M. Petrilli
!> Ivan P. Miranda
!
! DESCRIPTION:
!> Module to handle logging informations
!------------------------------------------------------------------------------

module safe_alloc_mod
   use string_mod, only: sl, int2str, fmt, lower
   use logger_mod, only: g_logger
   use string_mod
   use report_mod
   implicit none

   ! private

   type alloc_infos
      character(len=sl) :: label
      character(len=9) :: type
      integer :: size
   end type

   type :: safe_alloc
      type(alloc_infos), dimension(:), allocatable, private :: allocs_historic
   contains
      procedure, private :: push_new_info
      procedure, private :: allocate_int_0d, allocate_real4_0d, allocate_real8_0d, allocate_complex16_0d
      procedure, private :: allocate_int_1d, allocate_int_2d, allocate_int_3d, allocate_int_4d
      procedure, private :: allocate_real4_1d, allocate_real4_2d, allocate_real4_3d, allocate_real4_4d
      procedure, private :: allocate_real8_1d, allocate_real8_2d, allocate_real8_3d, allocate_real8_4d
      procedure, private :: allocate_complex16_1d, allocate_complex16_2d, allocate_complex16_3d, allocate_complex16_4d
      procedure, private :: deallocate_int_1d, deallocate_int_2d, deallocate_int_3d, deallocate_int_4d
      procedure, private :: deallocate_real4_1d, deallocate_real4_2d, deallocate_real4_3d, deallocate_real4_4d
      procedure, private :: deallocate_real8_1d, deallocate_real8_2d, deallocate_real8_3d, deallocate_real8_4d
      procedure, private :: deallocate_complex16_1d, deallocate_complex16_2d, deallocate_complex16_3d, deallocate_complex16_4d
      procedure, private :: move_alloc_int_1d, move_alloc_int_2d, move_alloc_int_3d, move_alloc_int_4d
      procedure, private :: move_alloc_real4_1d, move_alloc_real4_2d, move_alloc_real4_3d, move_alloc_real4_4d
      procedure, private :: move_alloc_real8_1d, move_alloc_real8_2d, move_alloc_real8_3d, move_alloc_real8_4d
      procedure, private :: move_alloc_complex16_1d, move_alloc_complex16_2d, move_alloc_complex16_3d, move_alloc_complex16_4d
      procedure, private :: report_allocate_int_1d, report_allocate_int_2d, report_allocate_int_3d, report_allocate_int_4d
      procedure, private :: report_allocate_real4_1d, report_allocate_real4_2d, report_allocate_real4_3d, report_allocate_real4_4d
      procedure, private :: report_allocate_real8_1d, report_allocate_real8_2d, report_allocate_real8_3d, report_allocate_real8_4d
      procedure, private :: report_allocate_complex16_1d, report_allocate_complex16_2d, report_allocate_complex16_3d, report_allocate_complex16_4d
      generic :: allocate => allocate_int_0d, allocate_real4_0d, allocate_real8_0d, allocate_complex16_0d, allocate_int_1d, allocate_int_2d, allocate_int_3d, allocate_int_4d, allocate_real4_1d, allocate_real4_2d, allocate_real4_3d, allocate_real4_4d, allocate_real8_1d, allocate_real8_2d, allocate_real8_3d, allocate_real8_4d, allocate_complex16_1d, allocate_complex16_2d, allocate_complex16_3d, allocate_complex16_4d
      generic :: deallocate => deallocate_int_1d, deallocate_int_2d, deallocate_int_3d, deallocate_int_4d, deallocate_real4_1d, deallocate_real4_2d, deallocate_real4_3d, deallocate_real4_4d, deallocate_real8_1d, deallocate_real8_2d, deallocate_real8_3d, deallocate_real8_4d, deallocate_complex16_1d, deallocate_complex16_2d, deallocate_complex16_3d, deallocate_complex16_4d
      generic :: move_alloc => move_alloc_int_1d, move_alloc_int_2d, move_alloc_int_3d, move_alloc_int_4d, move_alloc_real4_1d, move_alloc_real4_2d, move_alloc_real4_3d, move_alloc_real4_4d, move_alloc_real8_1d, move_alloc_real8_2d, move_alloc_real8_3d, move_alloc_real8_4d, move_alloc_complex16_1d, move_alloc_complex16_2d, move_alloc_complex16_3d, move_alloc_complex16_4d
      generic :: report_allocate => report_allocate_int_1d, report_allocate_int_2d, report_allocate_int_3d, report_allocate_int_4d, report_allocate_real4_1d, report_allocate_real4_2d, report_allocate_real4_3d, report_allocate_real4_4d, report_allocate_real8_1d, report_allocate_real8_2d, report_allocate_real8_3d, report_allocate_real8_4d, report_allocate_complex16_1d, report_allocate_complex16_2d, report_allocate_complex16_3d, report_allocate_complex16_4d
      procedure :: print_report => safe_alloc_report
      procedure :: get_allocations_total
      procedure :: get_allocations_memory
      procedure :: get_deallocations_total
      procedure :: get_deallocations_memory
      procedure :: get_remaining_total
      procedure :: get_remaining_memory
      final :: destructor
   end type

   interface safe_alloc
      procedure :: constructor
   end interface safe_alloc

   interface free
      procedure :: destructor
   end interface free

   interface get_type
      procedure :: get_type_int_1d, get_type_int_2d, get_type_int_3d, get_type_int_4d
      procedure :: get_type_real4_1d, get_type_real4_2d, get_type_real4_3d, get_type_real4_4d
      procedure :: get_type_real8_1d, get_type_real8_2d, get_type_real8_3d, get_type_real8_4d
      procedure :: get_type_complex16_1d, get_type_complex16_2d, get_type_complex16_3d, get_type_complex16_4d
   end interface get_type

   ! Global
   type(safe_alloc) :: g_safe_alloc

   ! Publics
   public :: g_safe_alloc, safe_alloc

contains

   function constructor() result(obj)
      type(safe_alloc) :: obj
      allocate (obj%allocs_historic(0))
   end function constructor

   subroutine destructor(this)
      type(safe_alloc) :: this
      character(len=sl) :: label
      integer :: i, j, allocated_size
      type(alloc_infos), dimension(:), allocatable :: new_allocs_historic
      do while (size(this%allocs_historic) > 1)
         label = this%allocs_historic(1)%label
         allocated_size = this%allocs_historic(1)%size
         ! list with new elements
         allocate (new_allocs_historic(size(this%allocs_historic)))
         j = 1
         do i = 2, size(this%allocs_historic)
            if (this%allocs_historic(i)%label == label) then
               allocated_size = allocated_size + this%allocs_historic(i)%size
            else
               new_allocs_historic(j) = this%allocs_historic(i)
               j = j + 1
            end if
         end do
         if (allocated_size > 0) then
            call g_logger%error('Possible memory leaky on list "'//trim(label)//'". Non-deallocated size: '//int2str(allocated_size), __FILE__, __LINE__)
         else if (allocated_size < 0) then
            call g_logger%error('More deallocations than allocations registered for list "'//trim(label)//'". Allocation-deallication balance: '//int2str(allocated_size), __FILE__, __LINE__)
         end if
         deallocate (this%allocs_historic)
         allocate (this%allocs_historic(j))
         this%allocs_historic = new_allocs_historic(:j)
      end do
   end subroutine destructor

   subroutine push_new_info(this, new_alloc)
      class(safe_alloc), intent(inout) :: this
      type(alloc_infos), intent(in) :: new_alloc
      type(alloc_infos), dimension(:), allocatable :: aux_allocs_historic
      integer :: new_size
      new_size = size(this%allocs_historic) + 1
      allocate (aux_allocs_historic(new_size - 1))
      aux_allocs_historic = this%allocs_historic
      deallocate (this%allocs_historic)
      allocate (this%allocs_historic(new_size))
      this%allocs_historic(:new_size - 1) = aux_allocs_historic
      this%allocs_historic(new_size) = new_alloc
      deallocate (aux_allocs_historic)
   end subroutine push_new_info

   function get_allocations_total(this) result(result)
      class(safe_alloc) :: this
      integer :: i, result, size_allocs_historic
      size_allocs_historic = size(this%allocs_historic)
      result = size_allocs_historic
      do i = 1, size_allocs_historic
         if (this%allocs_historic(i)%size > 0) result = result + 1
      end do
   end function get_allocations_total
   function get_allocations_memory(this, unit) result(result)
      class(safe_alloc), intent(in) :: this
      character, intent(in), optional :: unit
      character(len=2) :: unit_
      integer :: i, result, size_allocs_historic
      size_allocs_historic = size(this%allocs_historic)
      result = 0
      do i = 1, size_allocs_historic
         if (this%allocs_historic(i)%size > 0) result = result + this%allocs_historic(i)%size * get_size(this%allocs_historic(i)%type)
      end do
      !unit_ = merge(lower(unit), 'mb', present(unit))
      if (present(unit)) then
         unit_ = lower(unit)
      else
         unit_ = 'mb'
      end if
      select case (unit_)
      case ('kb')
         result = result / 1024
      case ('mb')
         result = result / 1024 / 1024
      case ('gb')
         result = result / 1024 / 1024 / 1024
      end select
   end function get_allocations_memory

   function get_deallocations_total(this) result(result)
      class(safe_alloc) :: this
      integer :: i, result, size_allocs_historic
      size_allocs_historic = size(this%allocs_historic)
      result = 0
      do i = 1, size_allocs_historic
         if (this%allocs_historic(i)%size < 0) result = result + 1
      end do
   end function get_deallocations_total
   function get_deallocations_memory(this, unit) result(result)
      class(safe_alloc), intent(in) :: this
      character(len=2), intent(in), optional :: unit
      character(len=2) :: unit_
      integer :: i, result, size_allocs_historic
      size_allocs_historic = size(this%allocs_historic)
      result = 0
      do i = 1, size_allocs_historic
         if (this%allocs_historic(i)%size < 0) result = result - this%allocs_historic(i)%size * get_size(this%allocs_historic(i)%type)
      end do
      !unit_ = merge(lower(unit), 'mb', present(unit))
      if (present(unit)) then
         unit_ = lower(unit)
      else
         unit_ = 'mb'
      end if
      select case (unit_)
      case ('kb')
         result = result / 1024
      case ('mb')
         result = result / 1024 / 1024
      case ('gb')
         result = result / 1024 / 1024 / 1024
      end select
   end function get_deallocations_memory

   function get_remaining_total(this) result(result)
      class(safe_alloc) :: this
      integer :: i, result, size_allocs_historic
      size_allocs_historic = size(this%allocs_historic)
      result = this%get_allocations_total() - this%get_deallocations_total()
   end function get_remaining_total
   function get_remaining_memory(this, unit) result(result)
      class(safe_alloc), intent(in) :: this
      character(len=2), intent(in), optional :: unit
      integer :: i, result, size_allocs_historic
      size_allocs_historic = size(this%allocs_historic)
      result = this%get_allocations_memory(unit) - this%get_deallocations_memory(unit)
   end function get_remaining_memory

   subroutine allocate_int_0d(this, label, list, list_size)
      integer, dimension(:), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/safe_alloc_allocate_0d.f90'
   end subroutine allocate_int_0d
   subroutine allocate_int_1d(this, label, list, list_size)
      integer, dimension(:), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/allocate_part1.f90'
      include 'include_codes/safe_alloc/allocate_part2_1d.f90'
      include 'include_codes/safe_alloc/allocate_part3.f90'
   end subroutine allocate_int_1d
   subroutine allocate_int_2d(this, label, list, list_size)
      integer, dimension(:, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/allocate_part1.f90'
      include 'include_codes/safe_alloc/allocate_part2_2d.f90'
      include 'include_codes/safe_alloc/allocate_part3.f90'
   end subroutine allocate_int_2d
   subroutine allocate_int_3d(this, label, list, list_size)
      integer, dimension(:, :, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/allocate_part1.f90'
      include 'include_codes/safe_alloc/allocate_part2_3d.f90'
      include 'include_codes/safe_alloc/allocate_part3.f90'
   end subroutine allocate_int_3d
   subroutine allocate_int_4d(this, label, list, list_size)
      integer, dimension(:, :, :, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/allocate_part1.f90'
      include 'include_codes/safe_alloc/allocate_part2_4d.f90'
      include 'include_codes/safe_alloc/allocate_part3.f90'
   end subroutine allocate_int_4d

   subroutine allocate_real4_0d(this, label, list, list_size)
      real(4), dimension(:), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/safe_alloc_allocate_0d.f90'
   end subroutine allocate_real4_0d
   subroutine allocate_real4_1d(this, label, list, list_size)
      real(4), dimension(:), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/allocate_part1.f90'
      include 'include_codes/safe_alloc/allocate_part2_1d.f90'
      include 'include_codes/safe_alloc/allocate_part3.f90'
   end subroutine allocate_real4_1d
   subroutine allocate_real4_2d(this, label, list, list_size)
      real(4), dimension(:, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/allocate_part1.f90'
      include 'include_codes/safe_alloc/allocate_part2_2d.f90'
      include 'include_codes/safe_alloc/allocate_part3.f90'
   end subroutine allocate_real4_2d
   subroutine allocate_real4_3d(this, label, list, list_size)
      real(4), dimension(:, :, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/allocate_part1.f90'
      include 'include_codes/safe_alloc/allocate_part2_3d.f90'
      include 'include_codes/safe_alloc/allocate_part3.f90'
   end subroutine allocate_real4_3d
   subroutine allocate_real4_4d(this, label, list, list_size)
      real(4), dimension(:, :, :, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/allocate_part1.f90'
      include 'include_codes/safe_alloc/allocate_part2_4d.f90'
      include 'include_codes/safe_alloc/allocate_part3.f90'
   end subroutine allocate_real4_4d

   subroutine allocate_real8_0d(this, label, list, list_size)
      real(8), dimension(:), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/safe_alloc_allocate_0d.f90'
   end subroutine allocate_real8_0d
   subroutine allocate_real8_1d(this, label, list, list_size)
      real(8), dimension(:), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/allocate_part1.f90'
      include 'include_codes/safe_alloc/allocate_part2_1d.f90'
      include 'include_codes/safe_alloc/allocate_part3.f90'
   end subroutine allocate_real8_1d
   subroutine allocate_real8_2d(this, label, list, list_size)
      real(8), dimension(:, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/allocate_part1.f90'
      include 'include_codes/safe_alloc/allocate_part2_2d.f90'
      include 'include_codes/safe_alloc/allocate_part3.f90'
   end subroutine allocate_real8_2d
   subroutine allocate_real8_3d(this, label, list, list_size)
      real(8), dimension(:, :, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/allocate_part1.f90'
      include 'include_codes/safe_alloc/allocate_part2_3d.f90'
      include 'include_codes/safe_alloc/allocate_part3.f90'
   end subroutine allocate_real8_3d
   subroutine allocate_real8_4d(this, label, list, list_size)
      real(8), dimension(:, :, :, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/allocate_part1.f90'
      include 'include_codes/safe_alloc/allocate_part2_4d.f90'
      include 'include_codes/safe_alloc/allocate_part3.f90'
   end subroutine allocate_real8_4d

   subroutine allocate_complex16_0d(this, label, list, list_size)
      complex(8), dimension(:), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/safe_alloc_allocate_0d.f90'
   end subroutine allocate_complex16_0d
   subroutine allocate_complex16_1d(this, label, list, list_size)
      complex(8), dimension(:), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/allocate_part1.f90'
      include 'include_codes/safe_alloc/allocate_part2_1d.f90'
      include 'include_codes/safe_alloc/allocate_part3.f90'
   end subroutine allocate_complex16_1d
   subroutine allocate_complex16_2d(this, label, list, list_size)
      complex(8), dimension(:, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/allocate_part1.f90'
      include 'include_codes/safe_alloc/allocate_part2_2d.f90'
      include 'include_codes/safe_alloc/allocate_part3.f90'
   end subroutine allocate_complex16_2d
   subroutine allocate_complex16_3d(this, label, list, list_size)
      complex(8), dimension(:, :, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/allocate_part1.f90'
      include 'include_codes/safe_alloc/allocate_part2_3d.f90'
      include 'include_codes/safe_alloc/allocate_part3.f90'
   end subroutine allocate_complex16_3d
   subroutine allocate_complex16_4d(this, label, list, list_size)
      complex(8), dimension(:, :, :, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/allocate_part1.f90'
      include 'include_codes/safe_alloc/allocate_part2_4d.f90'
      include 'include_codes/safe_alloc/allocate_part3.f90'
   end subroutine allocate_complex16_4d

   subroutine deallocate_int_1d(this, label, list)
      integer, dimension(:), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/deallocate.f90'
   end subroutine deallocate_int_1d
   subroutine deallocate_int_2d(this, label, list)
      integer, dimension(:, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/deallocate.f90'
   end subroutine deallocate_int_2d
   subroutine deallocate_int_3d(this, label, list)
      integer, dimension(:, :, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/deallocate.f90'
   end subroutine deallocate_int_3d
   subroutine deallocate_int_4d(this, label, list)
      integer, dimension(:, :, :, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/deallocate.f90'
   end subroutine deallocate_int_4d

   subroutine deallocate_real4_1d(this, label, list)
      real(4), dimension(:), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/deallocate.f90'
   end subroutine deallocate_real4_1d
   subroutine deallocate_real4_2d(this, label, list)
      real(4), dimension(:, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/deallocate.f90'
   end subroutine deallocate_real4_2d
   subroutine deallocate_real4_3d(this, label, list)
      real(4), dimension(:, :, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/deallocate.f90'
   end subroutine deallocate_real4_3d
   subroutine deallocate_real4_4d(this, label, list)
      real(4), dimension(:, :, :, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/deallocate.f90'
   end subroutine deallocate_real4_4d

   subroutine deallocate_real8_1d(this, label, list)
      real(8), dimension(:), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/deallocate.f90'
   end subroutine deallocate_real8_1d
   subroutine deallocate_real8_2d(this, label, list)
      real(8), dimension(:, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/deallocate.f90'
   end subroutine deallocate_real8_2d
   subroutine deallocate_real8_3d(this, label, list)
      real(8), dimension(:, :, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/deallocate.f90'
   end subroutine deallocate_real8_3d
   subroutine deallocate_real8_4d(this, label, list)
      real(8), dimension(:, :, :, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/deallocate.f90'
   end subroutine deallocate_real8_4d

   subroutine deallocate_complex16_1d(this, label, list)
      complex(8), dimension(:), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/deallocate.f90'
   end subroutine deallocate_complex16_1d
   subroutine deallocate_complex16_2d(this, label, list)
      complex(8), dimension(:, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/deallocate.f90'
   end subroutine deallocate_complex16_2d
   subroutine deallocate_complex16_3d(this, label, list)
      complex(8), dimension(:, :, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/deallocate.f90'
   end subroutine deallocate_complex16_3d
   subroutine deallocate_complex16_4d(this, label, list)
      complex(8), dimension(:, :, :, :), allocatable, intent(inout) :: list
      include 'include_codes/safe_alloc/deallocate.f90'
   end subroutine deallocate_complex16_4d

   subroutine move_alloc_int_1d(this, src, dst, label_src, label_dst)
      integer, dimension(:), allocatable, intent(inout) :: src, dst
      include 'include_codes/safe_alloc/move_alloc.f90'
   end subroutine move_alloc_int_1d
   subroutine move_alloc_int_2d(this, src, dst, label_src, label_dst)
      integer, dimension(:, :), allocatable, intent(inout) :: src, dst
      include 'include_codes/safe_alloc/move_alloc.f90'
   end subroutine move_alloc_int_2d
   subroutine move_alloc_int_3d(this, src, dst, label_src, label_dst)
      integer, dimension(:, :, :), allocatable, intent(inout) :: src, dst
      include 'include_codes/safe_alloc/move_alloc.f90'
   end subroutine move_alloc_int_3d
   subroutine move_alloc_int_4d(this, src, dst, label_src, label_dst)
      integer, dimension(:, :, :, :), allocatable, intent(inout) :: src, dst
      include 'include_codes/safe_alloc/move_alloc.f90'
   end subroutine move_alloc_int_4d

   subroutine move_alloc_real4_1d(this, src, dst, label_src, label_dst)
      real(4), dimension(:), allocatable, intent(inout) :: src, dst
      include 'include_codes/safe_alloc/move_alloc.f90'
   end subroutine move_alloc_real4_1d
   subroutine move_alloc_real4_2d(this, src, dst, label_src, label_dst)
      real(4), dimension(:, :), allocatable, intent(inout) :: src, dst
      include 'include_codes/safe_alloc/move_alloc.f90'
   end subroutine move_alloc_real4_2d
   subroutine move_alloc_real4_3d(this, src, dst, label_src, label_dst)
      real(4), dimension(:, :, :), allocatable, intent(inout) :: src, dst
      include 'include_codes/safe_alloc/move_alloc.f90'
   end subroutine move_alloc_real4_3d
   subroutine move_alloc_real4_4d(this, src, dst, label_src, label_dst)
      real(4), dimension(:, :, :, :), allocatable, intent(inout) :: src, dst
      include 'include_codes/safe_alloc/move_alloc.f90'
   end subroutine move_alloc_real4_4d

   subroutine move_alloc_real8_1d(this, src, dst, label_src, label_dst)
      real(8), dimension(:), allocatable, intent(inout) :: src, dst
      include 'include_codes/safe_alloc/move_alloc.f90'
   end subroutine move_alloc_real8_1d
   subroutine move_alloc_real8_2d(this, src, dst, label_src, label_dst)
      real(8), dimension(:, :), allocatable, intent(inout) :: src, dst
      include 'include_codes/safe_alloc/move_alloc.f90'
   end subroutine move_alloc_real8_2d
   subroutine move_alloc_real8_3d(this, src, dst, label_src, label_dst)
      real(8), dimension(:, :, :), allocatable, intent(inout) :: src, dst
      include 'include_codes/safe_alloc/move_alloc.f90'
   end subroutine move_alloc_real8_3d
   subroutine move_alloc_real8_4d(this, src, dst, label_src, label_dst)
      real(8), dimension(:, :, :, :), allocatable, intent(inout) :: src, dst
      include 'include_codes/safe_alloc/move_alloc.f90'
   end subroutine move_alloc_real8_4d

   subroutine move_alloc_complex16_1d(this, src, dst, label_src, label_dst)
      complex(8), dimension(:), allocatable, intent(inout) :: src, dst
      include 'include_codes/safe_alloc/move_alloc.f90'
   end subroutine move_alloc_complex16_1d
   subroutine move_alloc_complex16_2d(this, src, dst, label_src, label_dst)
      complex(8), dimension(:, :), allocatable, intent(inout) :: src, dst
      include 'include_codes/safe_alloc/move_alloc.f90'
   end subroutine move_alloc_complex16_2d
   subroutine move_alloc_complex16_3d(this, src, dst, label_src, label_dst)
      complex(8), dimension(:, :, :), allocatable, intent(inout) :: src, dst
      include 'include_codes/safe_alloc/move_alloc.f90'
   end subroutine move_alloc_complex16_3d
   subroutine move_alloc_complex16_4d(this, src, dst, label_src, label_dst)
      complex(8), dimension(:, :, :, :), allocatable, intent(inout) :: src, dst
      include 'include_codes/safe_alloc/move_alloc.f90'
   end subroutine move_alloc_complex16_4d

   subroutine safe_alloc_report(this)
      class(safe_alloc), intent(in) :: this
      type(report) :: rep
      integer :: i
      rep = report('total', integer_data=.true.)
      do i = 1, size(this%allocs_historic)
         if (this%allocs_historic(i)%size > 0) call rep%add_value(trim(rep%label)//'.'//this%allocs_historic(i)%label, this%allocs_historic(i)%size)
      end do
      call rep%print_report('MEMORY REPORT')
   end subroutine safe_alloc_report

   function get_size(type) result(size)
      integer :: size
      character(len=*) :: type
      select case (type)
      case ('integer')
         size = 4
      case ('real(4)')
         size = 4
      case ('real(8)')
         size = 8
      case ('complex(8)')
         size = 8
      end select
   end function get_size

   function get_type_int_1d(list) result(get_type)
      integer, dimension(:), allocatable :: list
      character(len=:), allocatable :: get_type
      get_type = 'integer'
   end function get_type_int_1d
   function get_type_int_2d(list) result(get_type)
      integer, dimension(:, :), allocatable :: list
      character(len=:), allocatable :: get_type
      get_type = 'integer'
   end function get_type_int_2d
   function get_type_int_3d(list) result(get_type)
      integer, dimension(:, :, :), allocatable :: list
      character(len=:), allocatable :: get_type
      get_type = 'integer'
   end function get_type_int_3d
   function get_type_int_4d(list) result(get_type)
      integer, dimension(:, :, :, :), allocatable :: list
      character(len=:), allocatable :: get_type
      get_type = 'integer'
   end function get_type_int_4d

   function get_type_real4_1d(list) result(get_type)
      real(4), dimension(:), allocatable :: list
      character(len=:), allocatable :: get_type
      get_type = 'real(4)'
   end function get_type_real4_1d
   function get_type_real4_2d(list) result(get_type)
      real(4), dimension(:, :), allocatable :: list
      character(len=:), allocatable :: get_type
      get_type = 'real(4)'
   end function get_type_real4_2d
   function get_type_real4_3d(list) result(get_type)
      real(4), dimension(:, :, :), allocatable :: list
      character(len=:), allocatable :: get_type
      get_type = 'real(4)'
   end function get_type_real4_3d
   function get_type_real4_4d(list) result(get_type)
      real(4), dimension(:, :, :, :), allocatable :: list
      character(len=:), allocatable :: get_type
      get_type = 'real(4)'
   end function get_type_real4_4d

   function get_type_real8_1d(list) result(get_type)
      real(8), dimension(:), allocatable :: list
      character(len=:), allocatable :: get_type
      get_type = 'real(8)'
   end function get_type_real8_1d
   function get_type_real8_2d(list) result(get_type)
      real(8), dimension(:, :), allocatable :: list
      character(len=:), allocatable :: get_type
      get_type = 'real(8)'
   end function get_type_real8_2d
   function get_type_real8_3d(list) result(get_type)
      real(8), dimension(:, :, :), allocatable :: list
      character(len=:), allocatable :: get_type
      get_type = 'real(8)'
   end function get_type_real8_3d
   function get_type_real8_4d(list) result(get_type)
      real(8), dimension(:, :, :, :), allocatable :: list
      character(len=:), allocatable :: get_type
      get_type = 'real(8)'
   end function get_type_real8_4d

   function get_type_complex16_1d(list) result(get_type)
      complex(8), dimension(:), allocatable :: list
      character(len=:), allocatable :: get_type
      get_type = 'complex(8)'
   end function get_type_complex16_1d
   function get_type_complex16_2d(list) result(get_type)
      complex(8), dimension(:, :), allocatable :: list
      character(len=:), allocatable :: get_type
      get_type = 'complex(8)'
   end function get_type_complex16_2d
   function get_type_complex16_3d(list) result(get_type)
      complex(8), dimension(:, :, :), allocatable :: list
      character(len=:), allocatable :: get_type
      get_type = 'complex(8)'
   end function get_type_complex16_3d
   function get_type_complex16_4d(list) result(get_type)
      complex(8), dimension(:, :, :, :), allocatable :: list
      character(len=:), allocatable :: get_type
      get_type = 'complex(8)'
   end function get_type_complex16_4d

   subroutine report_allocate_int_1d(this, label, list)
      integer, dimension(:), allocatable, intent(in) :: list
      include 'include_codes/safe_alloc/report_allocate.f90'
   end subroutine report_allocate_int_1d
   subroutine report_allocate_int_2d(this, label, list)
      integer, dimension(:, :), allocatable, intent(in) :: list
      include 'include_codes/safe_alloc/report_allocate.f90'
   end subroutine report_allocate_int_2d
   subroutine report_allocate_int_3d(this, label, list)
      integer, dimension(:, :, :), allocatable, intent(in) :: list
      include 'include_codes/safe_alloc/report_allocate.f90'
   end subroutine report_allocate_int_3d
   subroutine report_allocate_int_4d(this, label, list)
      integer, dimension(:, :, :, :), allocatable, intent(in) :: list
      include 'include_codes/safe_alloc/report_allocate.f90'
   end subroutine report_allocate_int_4d

   subroutine report_allocate_real4_1d(this, label, list)
      real(4), dimension(:), allocatable, intent(in) :: list
      include 'include_codes/safe_alloc/report_allocate.f90'
   end subroutine report_allocate_real4_1d
   subroutine report_allocate_real4_2d(this, label, list)
      real(4), dimension(:, :), allocatable, intent(in) :: list
      include 'include_codes/safe_alloc/report_allocate.f90'
   end subroutine report_allocate_real4_2d
   subroutine report_allocate_real4_3d(this, label, list)
      real(4), dimension(:, :, :), allocatable, intent(in) :: list
      include 'include_codes/safe_alloc/report_allocate.f90'
   end subroutine report_allocate_real4_3d
   subroutine report_allocate_real4_4d(this, label, list)
      real(4), dimension(:, :, :, :), allocatable, intent(in) :: list
      include 'include_codes/safe_alloc/report_allocate.f90'
   end subroutine report_allocate_real4_4d

   subroutine report_allocate_real8_1d(this, label, list)
      real(8), dimension(:), allocatable, intent(in) :: list
      include 'include_codes/safe_alloc/report_allocate.f90'
   end subroutine report_allocate_real8_1d
   subroutine report_allocate_real8_2d(this, label, list)
      real(8), dimension(:, :), allocatable, intent(in) :: list
      include 'include_codes/safe_alloc/report_allocate.f90'
   end subroutine report_allocate_real8_2d
   subroutine report_allocate_real8_3d(this, label, list)
      real(8), dimension(:, :, :), allocatable, intent(in) :: list
      include 'include_codes/safe_alloc/report_allocate.f90'
   end subroutine report_allocate_real8_3d
   subroutine report_allocate_real8_4d(this, label, list)
      real(8), dimension(:, :, :, :), allocatable, intent(in) :: list
      include 'include_codes/safe_alloc/report_allocate.f90'
   end subroutine report_allocate_real8_4d
   subroutine report_allocate_complex16_1d(this, label, list)
      complex(8), dimension(:), allocatable, intent(in) :: list
      include 'include_codes/safe_alloc/report_allocate.f90'
   end subroutine report_allocate_complex16_1d
   subroutine report_allocate_complex16_2d(this, label, list)
      complex(8), dimension(:, :), allocatable, intent(in) :: list
      include 'include_codes/safe_alloc/report_allocate.f90'
   end subroutine report_allocate_complex16_2d
   subroutine report_allocate_complex16_3d(this, label, list)
      complex(8), dimension(:, :, :), allocatable, intent(in) :: list
      include 'include_codes/safe_alloc/report_allocate.f90'
   end subroutine report_allocate_complex16_3d
   subroutine report_allocate_complex16_4d(this, label, list)
      complex(8), dimension(:, :, :, :), allocatable, intent(in) :: list
      include 'include_codes/safe_alloc/report_allocate.f90'
   end subroutine report_allocate_complex16_4d
end module safe_alloc_mod
