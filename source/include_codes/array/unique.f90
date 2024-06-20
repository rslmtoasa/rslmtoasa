integer :: i, j
integer :: arr_size, tmp_size
arr_size = size(arr)
allocate(tmp, mold=arr)
tmp(1) = arr(1)
tmp_size = 1
do i=2, arr_size
    if (.not.any(arr(i) == tmp(:i-1))) then
        tmp_size = tmp_size + 1
        tmp(tmp_size) = arr(i)
    endif
enddo
deallocate(arr)
allocate(arr(tmp_size))
arr = tmp(:tmp_size)