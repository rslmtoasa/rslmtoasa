integer :: i, j, arr_size
arr_size = size(arr)
do i=arr_size, 2, -1
  do j=2, i
    if (arr(j-1) > arr(j)) call swap(arr, j-1, j)
  enddo
enddo
