integer, intent(in), optional :: index
integer :: index_
index_ = merge(index, size(arr), present(index))
call move_alloc(arr, tmp)
allocate(arr(size(tmp) - 1))
arr(:index_ - 1) = tmp(:index_ - 1)
arr(index_:) = tmp(index_ + 1:)