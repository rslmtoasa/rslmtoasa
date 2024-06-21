class(namelist_generator), intent(inout) :: this
character(len=*), intent(in) :: name
integer, dimension(3) :: shape_value, ordered_index_shape_value
integer :: i, j, rank
shape_value = shape(value)
rank = size(shape_value)
ordered_index_shape_value = fqsortloc(shape_value)
do i = 1, shape_value(ordered_index_shape_value(1))
   do j = 1, shape_value(ordered_index_shape_value(2))
      select case (ordered_index_shape_value(rank))
      case (1)
         if (ordered_index_shape_value(1) == 2) then
            call this%add(name//'(:, '//fmt('I0', i)//', '//fmt('I0', j)//')', value(:, i, j))
         else
            call this%add(name//'(:, '//fmt('I0', j)//', '//fmt('I0', i)//')', value(:, j, i))
         end if
      case (2)
         if (ordered_index_shape_value(1) == 1) then
            call this%add(name//'('//fmt('I0', i)//', :, '//fmt('I0', j)//')', value(i, :, j))
         else
            call this%add(name//'('//fmt('I0', j)//', :, '//fmt('I0', i)//')', value(j, :, i))
         end if
      case (3)
         if (ordered_index_shape_value(1) == 1) then
            call this%add(name//'('//fmt('I0', i)//', '//fmt('I0', j)//', :)', value(i, j, :))
         else
            call this%add(name//'('//fmt('I0', j)//', '//fmt('I0', i)//', :)', value(j, i, :))
         end if
      end select
   end do
end do
