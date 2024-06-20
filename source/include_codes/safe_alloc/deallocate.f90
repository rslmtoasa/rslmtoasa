class(safe_alloc), intent(inout) :: this
character(len=*), intent(in) :: label
type(alloc_infos) :: new_alloc
if (allocated(list)) then
   new_alloc%label = label
   new_alloc%type = get_type(list)
   new_alloc%size = -product(shape(list))
   call this%push_new_info(new_alloc)
   deallocate (list)
else
   call g_logger%warning('Deallocating list "'//trim(label)//'" not allocated', '', 0)
end if
