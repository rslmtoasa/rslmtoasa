class(namelist_generator), intent(inout) :: this
character(len=*), intent(in) :: name
integer, dimension(3) :: shape_value, ordered_index_shape_value, index_offset
integer :: i, j, rank
shape_value = shape(value)
rank = size(shape_value)
ordered_index_shape_value = fqsortloc(shape_value)
index_offset = 1
if (present(lower_bounds)) index_offset = lower_bounds
do i = 1, shape_value(ordered_index_shape_value(1))
   do j = 1, shape_value(ordered_index_shape_value(2))
      select case (ordered_index_shape_value(rank))
      case (1)
         if (ordered_index_shape_value(1) == 2) then
            call this%add(name//'(:, '//fmt('I0', i + index_offset(2) - 1)//', '//fmt('I0', j + index_offset(3) - 1)//')', value(:, i, j))
         else
            call this%add(name//'(:, '//fmt('I0', j + index_offset(2) - 1)//', '//fmt('I0', i + index_offset(3) - 1)//')', value(:, j, i))
         end if
      case (2)
         if (ordered_index_shape_value(1) == 1) then
            call this%add(name//'('//fmt('I0', i + index_offset(1) - 1)//', :, '//fmt('I0', j + index_offset(3) - 1)//')', value(i, :, j))
         else
            call this%add(name//'('//fmt('I0', j + index_offset(1) - 1)//', :, '//fmt('I0', i + index_offset(3) - 1)//')', value(j, :, i))
         end if
      case (3)
         if (ordered_index_shape_value(1) == 1) then
            call this%add(name//'('//fmt('I0', i + index_offset(1) - 1)//', '//fmt('I0', j + index_offset(2) - 1)//', :)', value(i, j, :))
         else
            call this%add(name//'('//fmt('I0', j + index_offset(1) - 1)//', '//fmt('I0', i + index_offset(2) - 1)//', :)', value(j, i, :))
         end if
      end select
   end do
end do
