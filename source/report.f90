!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Report
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
!> Module to process and report information
!------------------------------------------------------------------------------

module report_mod
   use string_mod, only: sl, startswith, endswith, fcount_str, fjoin, split, str_contains, int2str, indent_str, upper, str_count, fmt
   use array_mod, only: push
   use precision_mod, only: rp
   implicit none

   private

   type report_data
      character(len=sl) :: label
      real(rp) :: value
   end type report_data

   type report
      character(len=sl) :: label
      real(rp), dimension(:), allocatable, private :: values
      integer, dimension(:), allocatable :: children_values_index
      type(report), dimension(:), allocatable :: children
      logical, private :: integer_data
   contains
      procedure, private :: get_labels_with_indented_lines => report_get_labels_with_indented_lines
      procedure, private :: report_add_value, report_add_label_value
      procedure, private :: report_add_value_r4, report_add_value_i
      procedure, private :: report_add_label_value_r4, report_add_label_value_i
      procedure, private :: push_child => report_push_child
      procedure, private :: report_get_values, report_get_values_with_label
      procedure, private :: get_child_values_private => report_get_child_values_private
      procedure, private :: report_get_child_values, report_get_child_values_label
      procedure, private :: report_get_num_children, report_get_num_children_label
      procedure, private :: report_get_size_values, report_get_size_values_label
      procedure, private :: set_values_on_leaf_alone => report_set_values_on_leaf_alone
      procedure, private :: report_sum_perc, report_sum, report_min, report_max, report_mean, report_total_ncalls, report_ncalls
      procedure, private :: report_has_child
      procedure, private :: get_values_private => report_get_values_private
      procedure :: is_leaf => report_is_leaf
      procedure :: assign => report_assign
      procedure :: print_report => report_print_report
      procedure :: add_child => report_add_child
      procedure :: get_labels => report_get_labels
      procedure :: get_child => report_get_child
      procedure :: get_valid_child => report_get_valid_child
      procedure :: get_parent => report_get_parent_label
      generic :: add_value => report_add_value, report_add_label_value, report_add_value_r4, report_add_value_i, report_add_label_value_r4, report_add_label_value_i
      generic :: get_values => report_get_values, report_get_values_with_label
      generic :: get_child_values => report_get_child_values, report_get_child_values_label
      generic :: get_num_children => report_get_num_children, report_get_num_children_label
      generic :: get_size_values => report_get_size_values, report_get_size_values_label
      generic :: sum_perc => report_sum_perc
      generic :: r_sum => report_sum
      generic :: r_min => report_min
      generic :: r_max => report_max
      generic :: mean => report_mean
      generic :: total_ncalls => report_total_ncalls
      generic :: ncalls => report_ncalls
      generic :: has_child => report_has_child
   end type report

   interface report
      procedure :: report_constructor
   end interface report

   public :: report

   character(len=sl), parameter :: default_fields = 'MEAN,MAX,MIN,TOTAL,TOTAL (%),NCALLS,TOTAL NCALLS'
   character(len=*), parameter :: report_category_label = 'CATEGORY'
   integer, parameter :: report_table_padding = 2
   integer, parameter :: report_field_margin = 3
   integer, parameter :: report_field_size = 7

contains
   recursive function report_constructor(label, integer_data) result(obj)
      character(len=*), intent(in) :: label
      character(len=sl), dimension(:), allocatable :: lst_labels
      type(report) :: obj
      logical, optional :: integer_data

      if (present(integer_data)) then
         obj%integer_data = integer_data
      else
         obj%integer_data = .false.
      end if
      !obj%integer_data = merge(integer_data, .False., present(integer_data))
      allocate (obj%values(0))
      allocate (obj%children_values_index(0))
      if (str_contains(label, '.')) then
         call split(label, '.', lst_labels)
         obj%label = lst_labels(1)
         allocate (obj%children(1))
         obj%children(1) = report(fjoin(lst_labels(2:), '.'))
      else
         allocate (obj%children(0))
         obj%label = label
      end if
   end function report_constructor

   !> Assigns obj to this
   recursive subroutine report_assign(this, obj)
      class(report), intent(inout) :: this
      type(report), intent(in) :: obj
      integer :: i
      this%label = obj%label
      this%integer_data = obj%integer_data
      if (allocated(obj%values)) then
         if (allocated(this%values)) deallocate (this%values)
         allocate (this%values(size(obj%values)))
         this%values = obj%values
      end if
      if (allocated(obj%children_values_index)) then
         if (allocated(this%children_values_index)) deallocate (this%children_values_index)
         allocate (this%children_values_index(size(obj%children_values_index)))
         this%children_values_index = obj%children_values_index
      end if
      if (allocated(obj%children)) then
         if (allocated(this%children)) deallocate (this%children)
         allocate (this%children(size(obj%children)))
         do i = 1, size(obj%children)
            call this%children(i)%assign(obj%children(i))
         end do
      end if
   end subroutine report_assign

   !> Append child argument variable to this report
   subroutine report_add_child(this, child)
      class(report), intent(inout) :: this
      type(report), intent(in) :: child
      call this%push_child(child)
   end subroutine report_add_child

   !> Private subroutine to push a new child into children list as it is
   subroutine report_push_child(this, child)
      class(report), intent(inout) :: this
      type(report), intent(in) :: child
      type(report), dimension(:), allocatable :: aux_children
      integer :: i, n_childrens
      n_childrens = size(this%children) + 1
      allocate (aux_children(n_childrens))
      do i = 1, n_childrens - 1
         call aux_children(i)%assign(this%children(i))
      end do
      call aux_children(n_childrens)%assign(child)
      deallocate (this%children)
      !TODO: try to use move_alloc
      allocate (this%children(n_childrens))
      do i = 1, n_childrens
         call this%children(i)%assign(aux_children(i))
      end do
   end subroutine report_push_child

   subroutine push_values(dst, src)
      real(rp), dimension(:), allocatable, intent(inout) :: dst
      real(rp), dimension(:), allocatable, intent(in) :: src
      real(rp), dimension(:), allocatable :: aux_values
      integer :: dst_size, src_size
      dst_size = size(dst)
      src_size = size(src)
      if (.not. allocated(dst)) then
         allocate (dst(src_size))
         dst = src
      else
         allocate (aux_values(dst_size + src_size))
         aux_values(:dst_size) = dst
         aux_values(dst_size + 1:) = src
         call move_alloc(aux_values, dst)
      end if
   end subroutine push_values

   recursive subroutine report_add_label_value(this, label, value)
      class(report), intent(inout) :: this
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: value
      type(report), pointer :: new_child, found_child
      type(report), target :: t_new_child
      character(len=sl):: new_label
      if (this%get_child(label, found_child)) then
         call found_child%add_value(value)
      else if (this%get_valid_child(label, found_child, new_label)) then
         call found_child%add_value(new_label, value)
      else
         new_label = label
         if (startswith(trim(label), trim(this%label)//'.')) new_label = label(len(trim(this%label)) + 2:)
         t_new_child = report(new_label)
         if (t_new_child%get_child(new_label, new_child)) call new_child%add_value(value)
         call this%add_child(t_new_child)
      end if
   end subroutine report_add_label_value

   subroutine report_add_value_r4(this, value)
      class(report), intent(inout) :: this
      real(4) :: value
      call this%add_value(real(value, rp))
   end subroutine report_add_value_r4

   subroutine report_add_value_i(this, value)
      class(report), intent(inout) :: this
      integer :: value
      call this%add_value(real(value, rp))
   end subroutine report_add_value_i

   subroutine report_add_label_value_r4(this, label, value)
      class(report), intent(inout) :: this
      character(len=*), intent(in) :: label
      real(4), intent(in) :: value
      call this%add_value(label, real(value, rp))
   end subroutine report_add_label_value_r4

   subroutine report_add_label_value_i(this, label, value)
      class(report), intent(inout) :: this
      character(len=*), intent(in) :: label
      integer, intent(in) :: value
      call this%add_value(label, real(value, rp))
   end subroutine report_add_label_value_i

   recursive function report_is_leaf(this, label) result(is_leaf)
      class(report), intent(in) :: this
      character(len=*), optional, intent(in) :: label
      type(report), pointer :: child
      logical :: is_leaf
      is_leaf = .false.
      if (present(label)) then
         if (this%get_child(label, child)) is_leaf = child%is_leaf()
      else
         is_leaf = .not. allocated(this%children) .or. size(this%children) == 0
      end if
   end function report_is_leaf

   subroutine report_set_values_on_leaf_alone(this)
      class(report), intent(inout) :: this
      type(report), pointer :: child
      type(report) :: child_aux
      character(len=sl), dimension(:), allocatable :: labels
      integer :: i, j, k
      call this%get_labels(labels)
      do i = size(labels), 1, -1
         if (this%get_child(labels(i), child)) then
            if (.not. child%is_leaf() .and. size(child%values) > 0) then
               child_aux = report(trim(child%label)//'*')
               j = 1
               do k = 1, child%children_values_index(j)
                  child%values(j) = child%values(j) - child%children(k)%r_sum()
               end do
               do j = 2, size(child%children_values_index)
                  do k = child%children_values_index(j - 1) + 1, child%children_values_index(j)
                     child%values(j) = child%values(j) - child%children(k)%r_sum()
                  end do
               end do
               call child%add_child(child%children(size(child%children)))
               do j = size(child%children), 2, -1
                  call child%children(j)%assign(child%children(j - 1))
               end do
               call move_alloc(child%values, child_aux%values)
               call child%children(1)%assign(child_aux)
               allocate (child%values(0))
            end if
         end if
      end do
   end subroutine report_set_values_on_leaf_alone

   subroutine report_add_value(this, value)
      class(report) :: this
      real(rp), intent(in) :: value
      call push(this%children_values_index, size(this%children))
      call push(this%values, value)
   end subroutine report_add_value

   recursive function report_get_child(this, label, child, asterisk) result(found)
      class(report), target, intent(in) :: this
      character(len=*), intent(in) :: label
      type(report), pointer, intent(out) :: child
      logical, optional, intent(in) :: asterisk
      character(len=sl), dimension(:), allocatable :: lst_labels
      logical :: asterisk_, found
      integer :: i
      found = .false.
      !asterisk_ = merge(asterisk, .False., present(asterisk))
      if (present(asterisk)) then
         asterisk_ = asterisk
      else
         asterisk_ = .false.
      end if
      if (this%label == label) then
         child => this
         found = .true.
         if (asterisk_) then
            call split(label, '.', lst_labels)
            call push(lst_labels, trim(lst_labels(size(lst_labels)))//'*')
            found = this%get_child(fjoin(lst_labels, '.'), child)
         end if
      else if (startswith(label, this%label) .and. allocated(this%children)) then
         call split(label, '.', lst_labels)
         if (size(lst_labels) > 1) then
            do i = 1, size(this%children)
               if (lst_labels(2) == this%children(i)%label) then
                  found = this%children(i)%get_child(fjoin(lst_labels(2:), '.'), child, asterisk_)
               end if
            end do
         end if
      end if
   end function report_get_child

   function report_get_valid_child(this, label, child, new_label) result(found)
      class(report), intent(in) :: this
      character(len=*), intent(in) :: label
      character(len=sl), intent(out), optional :: new_label
      type(report), intent(out), pointer :: child
      character(len=sl), dimension(:), allocatable :: lst_labels
      logical :: found
      integer :: i, init
      found = .false.
      call split(label, '.', lst_labels)
      init = merge(2, 1, lst_labels(1) == this%label)
      do i = size(lst_labels), init, -1
         if (this%get_child(fjoin(lst_labels(:i), '.'), child)) then
            found = .true.
            if (present(new_label)) then
               new_label = fjoin(lst_labels(i + 1:), '.')
            end if
            exit
         end if
      end do
   end function report_get_valid_child

   function report_get_num_children(this) result(num)
      class(report), intent(in) :: this
      integer :: num
      !num = merge(size(this%children), 0, allocated(this%children))
      if (allocated(this%children)) then
         num = size(this%children)
      else
         num = 0
      end if
   end function report_get_num_children

   function report_get_num_children_label(this, label) result(num)
      class(report), intent(in) :: this
      character(len=*), intent(in) :: label
      type(report), pointer :: child
      integer :: num
      num = 0
      if (this%get_child(label, child)) num = child%get_num_children()
   end function report_get_num_children_label

   function report_get_size_values(this) result(num)
      class(report), intent(in) :: this
      integer :: num
      !num = merge(size(this%values), 0, allocated(this%values))
      if (allocated(this%values)) then
         num = size(this%values)
      else
         num = 0
      end if
   end function report_get_size_values

   function report_get_size_values_label(this, label) result(num)
      class(report), intent(in) :: this
      character(len=*), intent(in) :: label
      type(report), pointer :: child
      integer :: num
      num = 0
      if (this%get_child(label, child)) num = child%get_size_values()
   end function report_get_size_values_label

   recursive subroutine report_get_labels(this, labels, root, append, depth, indent_level, indented_lines, levels)
      class(report), intent(in) :: this
      character(len=sl), dimension(:), allocatable, intent(inout) :: labels
      character(len=*), optional, intent(in) :: root
      integer, optional, intent(in) :: depth, indent_level
      logical, optional, intent(in) :: append, indented_lines
      integer, dimension(:), allocatable, optional, intent(out) :: levels
      integer, dimension(:), allocatable :: levels_
      integer :: depth_, indent_level_
      logical :: append_, indented_lines_
      character(len=sl) :: root_
      character(len=sl), dimension(:), allocatable :: lst_labels
      integer :: i, children_size
      !append_ = merge(append, .False., present(append))
      if (present(append)) then
         append_ = append
      else
         append_ = .false.
      end if
      !depth_ = merge(depth, -2, present(depth))
      if (present(depth)) then
         depth_ = depth
      else
         depth_ = -2
      end if
      !indent_level_ = merge(indent_level, -1, present(indent_level))
      if (present(indent_level)) then
         indent_level_ = indent_level
      else
         indent_level_ = -1
      end if
      !indented_lines_ = merge(indented_lines, .False., present(indented_lines))
      if (present(indented_lines)) then
         indented_lines_ = indented_lines
      else
         indented_lines_ = .false.
      end if
      if (indented_lines_) then
         call this%get_labels(labels, root=root, append=append, depth=depth, indent_level=0, levels=levels_)
         do i = 1, size(labels)
            labels(i) = this%get_labels_with_indented_lines(i, labels(i), levels_)
         end do
      else
         if (.not. append_ .and. allocated(labels)) deallocate (labels)
         if (.not. allocated(labels)) allocate (labels(0))
         if (depth_ == -1) return
         root_ = this%label
         if (present(root)) root_ = trim(root)//'.'//trim(this%label)
         call push_labels(labels, root_)
         if (allocated(this%children)) then
            children_size = size(this%children)
            do i = 1, children_size
               call this%children(i)%get_labels(labels, trim(root_), depth=depth_ - 1, append=.true.)
            end do
         end if
         if (present(levels) .or. indent_level_ > -1) then
            call this%get_labels(labels, root=root, append=append, depth=depth)
            if (allocated(levels_)) deallocate (levels_)
            allocate (levels_(size(labels)))
            do i = 1, size(levels_)
               levels_(i) = str_count(labels(i), '.')
            end do
         end if
         if (indent_level_ > -1) then
            do i = 1, size(labels)
               call split(labels(i), '.', lst_labels)
               if (levels_(i) > 0 .and. indent_level_ > 0) then
                  labels(i) = fmt(int2str(levels_(i) * indent_level_)//'X, A', lst_labels(size(lst_labels)))
               else
                  labels(i) = fmt('A', lst_labels(size(lst_labels)))
               end if
            end do
         end if
         if (present(levels)) call move_alloc(levels_, levels)
      end if
   end subroutine report_get_labels

   function report_get_labels_with_indented_lines(this, i, label, levels) result(result_label)
      class(report), intent(in) :: this
      integer, intent(in) :: i
      integer, dimension(:), allocatable :: levels
      character(len=*) :: label
      character(len=:), allocatable :: result_label
      logical :: has_next_item, is_last_item
      integer :: curr_level, lvl, k
      character(len=*), parameter :: blanck = '   '
      character(len=*), parameter :: no_item = '│  '
      character(len=*), parameter :: continue_item = '├─ '
      character(len=*), parameter :: final_item = '└─ '

      curr_level = levels(i)
      is_last_item = i == size(levels)
      result_label = ''
      do lvl = 1, curr_level
         if (lvl == curr_level) then
            has_next_item = .not. is_last_item
            if (has_next_item) then
               has_next_item = .false.
               do k = i + 1, size(levels)
                  if (levels(k) <= curr_level) then
                     has_next_item = (levels(k) == curr_level)
                     exit
                  end if
               end do
            end if
            if (has_next_item) then
               result_label = result_label//continue_item
            else
               result_label = result_label//final_item
            end if
         else
            has_next_item = .false.
            do k = i + 1, size(levels)
               if (levels(k) <= lvl) then
                  has_next_item = (levels(k) == lvl)
                  exit
               end if
            end do
            if (has_next_item) then
               result_label = result_label//no_item
            else
               result_label = result_label//blanck
            end if
         end if
      end do
      result_label = result_label//trim(label)
   end function report_get_labels_with_indented_lines

   subroutine push_labels(labels, label)
      character(len=sl), dimension(:), allocatable, intent(inout) :: labels
      character(len=*), intent(in) :: label
      character(len=sl), dimension(:), allocatable :: aux_labels
      integer :: labels_size
      labels_size = size(labels) + 1
      allocate (aux_labels(labels_size))
      aux_labels(:labels_size - 1) = labels
      aux_labels(labels_size) = label
      deallocate (labels)
      call move_alloc(aux_labels, labels)
   end subroutine push_labels

   !> Search for a child with 'label', if found returns .True. and fills 'values' array with child's values, otherwise returns .False.
   !> You should enter child levels in 'label', ex: 'calc1.subcalc2.subsubcalc3' for getting the values of 'subsubcalc3'
   subroutine report_get_values(this, values, append, root, children, depth)
      class(report), intent(in) :: this
      real(rp), dimension(:), allocatable, intent(inout) :: values
      logical, optional, intent(in) :: append, root, children
      integer, optional, intent(in) :: depth
      logical :: append_, root_, children_
      integer :: depth_

      if (present(depth)) then
         depth_ = depth
      else
         depth_ = -1
      end if
      !depth_ = merge(depth, -1, present(depth))
      if (depth_ == 0) return

      if (present(append)) then
         append_ = append
      else
         append_ = .false.
      end if
      !append_ = merge(append, .False., present(append))

      if (present(root)) then
         root_ = root
      else
         root_ = .true.
      end if
      !root_ = merge(root, .True., present(root))

      if (present(children)) then
         children_ = children
      else
         children_ = .false.
      end if
      !children_ = merge(children, .False., present(children))
      if (.not. append_ .or. .not. allocated(values)) then
         if (allocated(values)) deallocate (values)
         allocate (values(0))
      end if
      if (root_ .and. allocated(this%values)) call push_values(values, this%values)
      if (children_) call this%get_values_private(values, depth_)
   end subroutine report_get_values

   subroutine report_get_values_with_label(this, label, values, append, root, children, depth)
      class(report), intent(in) :: this
      character(len=*), intent(in) :: label
      real(rp), dimension(:), allocatable, intent(inout) :: values
      logical, optional, intent(in) :: append, root, children
      integer, optional, intent(in) :: depth
      type(report), pointer :: child
      if (this%get_child(label, child)) call child%get_values(values, append, root, children, depth)
   end subroutine report_get_values_with_label

   subroutine report_get_child_values(this, values, append)
      class(report), intent(in) :: this
      real(rp), dimension(:), allocatable, intent(inout) :: values
      logical, optional, intent(in) :: append
      if (.not. (present(append) .and. append)) then
         if (allocated(values)) deallocate (values)
         allocate (values(0))
      end if
      call this%get_child_values_private(values)
   end subroutine report_get_child_values

   subroutine report_get_child_values_label(this, label, values, append)
      class(report), intent(in) :: this
      character(len=*), intent(in) :: label
      real(rp), dimension(:), allocatable, intent(inout) :: values
      logical, optional, intent(in) :: append
      type(report), pointer :: child
      if (.not. (present(append) .and. append)) then
         if (allocated(values)) deallocate (values)
         allocate (values(0))
      end if
      if (this%get_child(label, child)) call child%get_child_values_private(values)
   end subroutine report_get_child_values_label

   recursive subroutine report_get_child_values_private(this, values)
      class(report), intent(in) :: this
      real(rp), dimension(:), allocatable, intent(inout) :: values
      integer :: i
      if (.not. allocated(values)) allocate (values(0))
      if (allocated(this%children)) then
         do i = 1, size(this%children)
            call push_values(values, this%children(i)%values)
            call this%children(i)%get_child_values_private(values)
         end do
      end if
   end subroutine report_get_child_values_private

   recursive subroutine report_get_values_private(this, values, depth)
      class(report), intent(in) :: this
      real(rp), dimension(:), allocatable, intent(inout) :: values
      integer, intent(in) :: depth
      integer :: i
      if (depth == 0) return
      if (allocated(this%children)) then
         do i = 1, size(this%children)
            call push_values(values, this%children(i)%values)
            call this%children(i)%get_values_private(values, depth - 1)
         end do
      end if
   end subroutine report_get_values_private

   function report_get_parent_label(this, label, parent) result(found)
      class(report), intent(in) :: this
      character(len=*), intent(in) :: label
      type(report), pointer, intent(out) :: parent
      character(len=sl) :: parent_label
      character(len=sl), dimension(:), allocatable :: lst_labels
      logical :: found
      call split(label, '.', lst_labels)
      parent_label = fjoin(lst_labels(:size(lst_labels) - 1), '.')
      found = this%get_child(parent_label, parent)
   end function report_get_parent_label

   subroutine get_short_labels(labels, indented_labels, levels)
      character(len=sl), dimension(:), allocatable, intent(in) :: labels
      character(len=sl), dimension(:), allocatable, intent(out) :: indented_labels
      integer, dimension(:), allocatable, intent(out) :: levels
      character(len=sl), dimension(:), allocatable :: lst_labels
      integer :: i
      allocate (indented_labels(size(labels)))
      allocate (levels(size(labels)))
      do i = 1, size(labels)
         call split(labels(i), '.', lst_labels)
         levels(i) = size(lst_labels)
         indented_labels(i) = lst_labels(size(lst_labels))
      end do
   end subroutine get_short_labels

   subroutine report_print_report(this, title, fields)
      class(report), intent(in) :: this
      character(len=*), intent(in) :: title
      character(len=*), intent(in), optional :: fields
      character(len=sl) :: fields_
      character(len=sl), dimension(:), allocatable :: lst_fields, labels, indented_labels
      type(report) :: this_copy
      integer :: i, j, indented_label_max_size, fields_max_size, report_width
      logical :: should_print_comment

      interface print_field
         procedure :: print_field_c
         procedure :: print_field_i
         procedure :: print_field_r4
         procedure :: print_field_r8
      end interface print_field

      interface get_field
         procedure :: get_field_c
         procedure :: get_field_i
         procedure :: get_field_r4
         procedure :: get_field_r8
      end interface get_field

      call this_copy%assign(this)
      call this_copy%set_values_on_leaf_alone()

      ! fields format
      if (present(fields)) then
         fields_ = fields
      else
         fields_ = default_fields
      end if
      call split(fields_, ',', lst_fields)
      fields_max_size = 0
      do i = 1, size(lst_fields)
         fields_max_size = max(fields_max_size, len(trim(lst_fields(i))))
      end do
      ! get labels
      call this_copy%get_labels(labels)
      call this_copy%get_labels(indented_labels, indented_lines=.true.)
      indented_label_max_size = 0
      do i = 1, size(indented_labels)
         indented_label_max_size = max(indented_label_max_size, label_len(indented_labels(i)))
      end do
      report_width = 2 * report_table_padding + indented_label_max_size + report_field_margin + size(lst_fields) * (report_field_margin + fields_max_size)
#ifdef DEBUG
      call print_full_info()
#endif
      call print_table_header()
      call print_table_content()
      call print_table_footer()

   contains
      subroutine print_full_info()
         real(rp), dimension(:), allocatable :: values
         integer :: label_max_size

         call print_dashed_line(report_width, title)
         label_max_size = 0
         do i = 1, size(labels)
            label_max_size = max(label_max_size, len(trim(labels(i))))
         end do
         do i = 1, size(labels)
            call this_copy%get_values(labels(i), values)
            write (*, '(A'//int2str(label_max_size)//', " ", '//int2str(label_max_size)//'F'//int2str(fields_max_size)//'.4)') adjustl(labels(i)), values
         end do
      end subroutine print_full_info

      subroutine print_dashed_line(width, title)
         integer, intent(in) :: width
         character(len=*), optional, intent(in) :: title
         integer :: with_, title_len
         if (present(title)) then
            title_len = len(trim(title))
            with_ = (width - 3 - title_len) / 2
            write (*, '("'//repeat('-', with_ - 1)//' ", A, " '//repeat('-', with_ - merge(1, 0, mod(title_len, 2) /= 0))//'")') trim(upper(title))
         else
            write (*, '("'//repeat('-', width - 3)//'")')
         end if
      end subroutine print_dashed_line

      subroutine print_table_header()
         call print_dashed_line(report_width, title)
         write (*, '('//int2str(report_table_padding)//'X, A, '//int2str(max(1, indented_label_max_size - len(report_category_label) + report_field_margin))//'X)', advance='no') trim(report_category_label)
         do i = 1, size(lst_fields)
            lst_fields(i) = upper(lst_fields(i))
            write (*, '(A'//int2str(fields_max_size)//', '//int2str(report_field_margin)//'X)', advance='no') trim(lst_fields(i))
         end do
         write (*, *)
      end subroutine print_table_header

      subroutine print_table_content()
         type(report), pointer :: child
         integer :: aux
         should_print_comment = .false.
         do i = 1, size(labels)
            if (endswith(labels(i), '*')) should_print_comment = .true.
            ! call this_copy%get_values(labels(i), values, children=.True.)
            write (*, '('//int2str(report_table_padding)//'X, A'//int2str(indented_label_max_size + label_extra_space(indented_labels(i)) + report_field_margin)//')', advance='no') indented_labels(i)
            do j = 1, size(lst_fields)
               select case (lst_fields(j))
               case ('MAX')
                  if (this_copy%ncalls(labels(i)) == 1 .or. endswith(labels(i), '*')) then
                     call print_field('-')
                  else
                     if (this_copy%integer_data) then
                        call print_field(int(this_copy%r_max(labels(i))))
                     else
                        call print_field(this_copy%r_max(labels(i)))
                     end if
                  end if
               case ('MIN')
                  if (this_copy%ncalls(labels(i)) == 1 .or. endswith(labels(i), '*')) then
                     call print_field('-')
                  else
                     if (this_copy%integer_data) then
                        call print_field(int(this_copy%r_min(labels(i))))
                     else
                        call print_field(this_copy%r_min(labels(i)))
                     end if
                  end if
               case ('MEAN')
                  if (this_copy%ncalls(labels(i)) == 1 .or. endswith(labels(i), '*')) then
                     call print_field('-')
                  else
                     call print_field(this_copy%mean(labels(i)))
                  end if
               case ('TOTAL')
                  call print_field(this_copy%r_sum(labels(i)))
               case ('TOTAL (%)')
                  call print_field(this_copy%sum_perc(labels(i), tottot=this_copy%r_sum(labels(1))), decimal=1)
               case ('NCALLS')
                  aux = this_copy%ncalls(labels(i))
                  if (endswith(labels(i), '*')) then
                     call print_field('-')
                  else
                     call print_field(aux)
                  end if
               case ('TOTAL NCALLS')
                  if (this_copy%get_child(labels(i), child, asterisk=.true.)) then
                     call print_field(this_copy%total_ncalls(labels(i)) - child%total_ncalls())
                  else
                     call print_field(this_copy%total_ncalls(labels(i)))
                  end if
               end select
            end do
            write (*, *) ! new line
         end do
      end subroutine print_table_content

      subroutine print_table_footer()
         call print_dashed_line(report_width)
         if (should_print_comment) then
            write (*, '('//int2str(report_table_padding)//'X, A)') '* exclusive values to this label'
            call print_dashed_line(report_width)
         end if
      end subroutine print_table_footer

      function label_extra_space(label) result(qtd)
         character(len=*) :: label
         integer :: qtd
         qtd = (len('├') - 1) * str_count(label, '├') + &
               (len('└') - 1) * str_count(label, '└') + &
               (len('─') - 1) * str_count(label, '─') + &
               (len('│') - 1) * str_count(label, '│')
      end function label_extra_space

      function label_len(label)
         character(len=*) :: label
         integer :: label_len
         label_len = len(trim(label)) - &
                     (len('├') - 1) * str_count(label, '├') - &
                     (len('└') - 1) * str_count(label, '└') - &
                     (len('─') - 1) * str_count(label, '─') - &
                     (len('│') - 1) * str_count(label, '│')
      end function label_len

      subroutine print_field_c(value)
         character(len=*), intent(in) :: value
         write (*, '(A'//int2str(fields_max_size)//', '//int2str(report_field_margin)//'X)', advance='no') value
      end subroutine print_field_c
      subroutine print_field_i(value)
         integer, intent(in) :: value
         write (*, '(I'//int2str(fields_max_size)//', '//int2str(report_field_margin)//'X)', advance='no') value
      end subroutine print_field_i
      subroutine print_field_r4(value, decimal)
         real(4), intent(in) :: value
         integer, intent(in), optional :: decimal
         integer :: decimal_
         if (present(decimal)) then
            decimal_ = decimal
         else
            decimal_ = 4
         end if
         !decimal_ = merge(decimal, 4, present(decimal))
         write (*, '(F'//int2str(fields_max_size)//'.'//int2str(decimal_)//', '//int2str(report_field_margin)//'X)', advance='no') value
      end subroutine print_field_r4
      subroutine print_field_r8(value, decimal)
         real(8), intent(in) :: value
         integer, intent(in), optional :: decimal
         integer :: decimal_
         if (present(decimal)) then
            decimal_ = decimal
         else
            decimal_ = 4
         end if
         !decimal_ = merge(decimal, 4, present(decimal))
         write (*, '(F'//int2str(fields_max_size)//'.'//int2str(decimal_)//', '//int2str(report_field_margin)//'X)', advance='no') value
      end subroutine print_field_r8

      function get_field_c(value) result(field)
         character(len=*), intent(in) :: value
         character(len=sl) :: field
         write (field, '(A'//int2str(fields_max_size)//', '//int2str(report_field_margin)//'X)') value
      end function get_field_c
      function get_field_i(value) result(field)
         integer, intent(in) :: value
         character(len=sl) :: field
         write (field, '(I'//int2str(fields_max_size)//', '//int2str(report_field_margin)//'X)') value
      end function get_field_i
      function get_field_r4(value) result(field)
         real(4), intent(in) :: value
         character(len=sl) :: field
         write (field, '(F'//int2str(fields_max_size)//'.4, '//int2str(report_field_margin)//'X)') value
      end function get_field_r4
      function get_field_r8(value) result(field)
         real(8), intent(in) :: value
         character(len=sl) :: field
         write (field, '(F'//int2str(fields_max_size)//'.4, '//int2str(report_field_margin)//'X)') value
      end function get_field_r8
   end subroutine report_print_report

   recursive function report_sum_perc(this, label, root, children, depth, tottot) result(val)
      class(report), intent(in) :: this
      character(len=*), optional, intent(in) :: label
      logical, optional, intent(in) :: root, children
      integer, optional, intent(in) :: depth
      real(rp), optional, intent(in) :: tottot  !< Optional value used to calculate total percentages
      real(rp) :: val, total_val

      type(report), pointer :: child, parent
      logical :: children_

      if (present(children)) then
         children_ = children
      else
         children_ = .true.
      end if
      !children_ = merge(children, .True., present(children))
      val = 100
      if (present(label)) then
         if (this%get_child(label, child) .and. (.not. children_ .or. this%get_parent(label, parent))) then
            val = 100 * child%r_sum(root=root, children=children_, depth=depth)
            if (children_) then
               if (present(depth)) then
                  total_val = parent%r_sum(depth=depth + 1)
               else
                  total_val = parent%r_sum()
               end if
            else
               total_val = child%r_sum(depth=depth)
            end if
            if (present(tottot)) then
               val = val / tottot
            else
               val = val / total_val
            end if
         end if
      else
         val = 100
      end if
   end function report_sum_perc

   recursive function report_sum(this, label, root, children, depth) result(val)
      class(report), intent(in) :: this
      character(len=*), optional, intent(in) :: label
      logical, optional, intent(in) :: root, children
      integer, optional, intent(in) :: depth
      real(rp) :: val
      real(rp), dimension(:), allocatable :: values
      type(report), pointer :: child
      logical :: children_

      if (present(children)) then
         children_ = children
      else
         children_ = .true.
      end if
      !children_ = merge(children, .True., present(children))
      if (present(label)) then
         if (this%get_child(label, child)) val = child%r_sum(root=root, children=children_, depth=depth)
      else
         call this%get_values(values, root=root, children=children_, depth=depth)
         val = sum(values)
      end if
   end function report_sum

   recursive function report_min(this, label, root, children, depth) result(val)
      class(report), intent(in) :: this
      character(len=*), optional, intent(in) :: label
      logical, optional, intent(in) :: root, children
      integer, optional, intent(in) :: depth
      real(rp) :: val
      real(rp), dimension(:), allocatable :: values
      type(report), pointer :: child
      logical :: children_

      if (present(children)) then
         children_ = children
      else
         children_ = .true.
      end if
      !children_ = merge(children, .True., present(children))
      if (present(label)) then
         if (this%get_child(label, child)) val = child%r_min(root=root, children=children_, depth=depth)
      else
         call this%get_values(values, root=root, children=children_, depth=depth)
         val = minval(values)
      end if
   end function report_min

   recursive function report_max(this, label, root, children, depth) result(val)
      class(report), intent(in) :: this
      character(len=*), optional, intent(in) :: label
      logical, optional, intent(in) :: root, children
      integer, optional, intent(in) :: depth
      real(rp) :: val
      real(rp), dimension(:), allocatable :: values
      type(report), pointer :: child
      logical :: children_

      if (present(children)) then
         children_ = children
      else
         children_ = .true.
      end if
      !children_ = merge(children, .True., present(children))
      if (present(label)) then
         if (this%get_child(label, child)) val = child%r_max(root=root, children=children_, depth=depth)
      else
         call this%get_values(values, root=root, children=children_, depth=depth)
         val = maxval(values)
      end if
   end function report_max

   recursive function report_mean(this, label, root, children, depth) result(val)
      class(report), intent(in) :: this
      character(len=*), optional, intent(in) :: label
      logical, optional, intent(in) :: root, children
      integer, optional, intent(in) :: depth
      real(rp) :: val
      real(rp), dimension(:), allocatable :: values
      type(report), pointer :: child
      logical :: children_

      if (present(children)) then
         children_ = children
      else
         children_ = .true.
      end if
      !children_ = merge(children, .True., present(children))
      if (present(label)) then
         if (this%get_child(label, child)) then
            call child%get_values(values, root=root, children=children_, depth=depth)
            val = sum(values) / size(values)
         else
            val = 100
         end if
      else
         call this%get_values(values, root=root, children=children_, depth=depth)
         val = sum(values) / size(values)
      end if
   end function report_mean

   recursive function report_total_ncalls(this, label, root, children, depth) result(val)
      class(report), intent(in) :: this
      character(len=*), optional, intent(in) :: label
      logical, optional, intent(in) :: root, children
      integer, optional, intent(in) :: depth
      integer :: val
      real(rp), dimension(:), allocatable :: values
      type(report), pointer :: child
      logical :: children_

      if (present(children)) then
         children_ = children
      else
         children_ = .true.
      end if
      !children_ = merge(children, .True., present(children))
      if (present(label)) then
         if (this%get_child(label, child)) val = child%total_ncalls(root=root, children=children_, depth=depth)
      else
         call this%get_values(values, root=root, children=children_, depth=depth)
         val = size(values)
      end if
   end function report_total_ncalls

   recursive function report_ncalls(this, label, root, children, depth) result(val)
      class(report), intent(in) :: this
      character(len=*), optional, intent(in) :: label
      logical, optional, intent(in) :: root, children
      integer, optional, intent(in) :: depth
      integer :: val

      type(report), pointer :: child
      logical :: children_
      type(report) :: this_copy

      if (present(children)) then
         children_ = children
      else
         children_ = .true.
      end if
      !children_ = merge(children, .True., present(children))
      if (present(label)) then
         if (this%get_child(label, child)) val = child%ncalls(root=root, children=children_, depth=depth)
      else
         call this_copy%assign(this)
         call this_copy%set_values_on_leaf_alone()

         if (this_copy%is_leaf(label)) then
            val = this_copy%total_ncalls(label, children=.false.)
         else
            if (this_copy%get_child(label, child, asterisk=.true.)) then
               val = child%total_ncalls()
            else
               val = 1
            end if
         end if
      end if
   end function report_ncalls

   function report_has_child(this, label) result(val)
      class(report), intent(in) :: this
      character(len=*), optional, intent(in) :: label
      type(report), pointer :: child
      logical :: val
      val = this%get_child(label, child)
   end function report_has_child
end module report_mod
