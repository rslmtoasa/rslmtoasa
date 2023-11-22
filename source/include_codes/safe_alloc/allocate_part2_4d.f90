integer, dimension(4), intent(in) :: list_size
if (allocated(list)) then
    deallocate(list)
    call g_logger%warning('Reallocating list "' // trim(label) // '" of size ' // int2str(size(list)) // ' to new size ' // int2str(list_size(1)) // ', ' // int2str(list_size(2)) // ', ' // int2str(list_size(3)) // ', ' // int2str(list_size(4)), 'safe_alloc.f90', -1)
endif
allocate(list(list_size(1), list_size(2), list_size(3), list_size(4)))