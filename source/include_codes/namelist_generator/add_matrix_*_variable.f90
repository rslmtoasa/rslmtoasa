class(namelist_generator), intent(inout) :: this
character(len=*), intent(in) :: name
integer, dimension(2) :: value_shape
integer :: i
value_shape = shape(value)
if (size(value) > 0 .and. minval(value_shape) < MAX_ROWS .and. maxval(value_shape) < MAX_COLS) then
   if (value_shape(1) <= value_shape(2)) then
      do i = 1, value_shape(1)
         call this%add(name//'('//fmt('I0', i)//', :)', value(i, :))
      end do
   else
      do i = 1, value_shape(2)
         call this%add(name//'(:, '//fmt('I0', i)//')', value(:, i))
      end do
   end if
end if
