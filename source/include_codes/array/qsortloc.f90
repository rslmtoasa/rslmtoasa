integer, dimension(:), allocatable :: arr_idx
integer :: i, j, arr_size
arr_size = size(arr)
allocate (arr_idx(arr_size))
do i = 1, arr_size
   arr_idx(i) = i
end do
do i = arr_size, 2, -1
   do j = 2, i
      if (arr(arr_idx(j - 1)) > arr(arr_idx(j))) then
         call swap(arr_idx, j - 1, j)
      end if
   end do
end do
