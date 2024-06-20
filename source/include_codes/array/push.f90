integer, intent(in), optional :: index
integer :: index_, ar_size
if(.not.allocated(arr)) allocate(arr(0))
ar_size = size(arr)
    if(present(index)) then
       index_ = index
    else
       index_ = ar_size+1
    end if
!index_ = merge(index, ar_size + 1, present(index))
call move_alloc(arr, tmp)
allocate(arr(ar_size + 1))
arr(:index_-1) = tmp(:index_-1)
arr(index_) = new
arr(index_+1:) = tmp(index_:)
