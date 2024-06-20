integer, dimension(1), intent(in) :: list_size
if (allocated(list)) then
   deallocate (list)
   call g_logger%warning('Reallocating list "'//trim(label)//'" of size '//int2str(size(list))//' to new size '//int2str(list_size(1)), 'safe_alloc.f90', -1)
end if
allocate (list(list_size(1)))
