class(safe_alloc), intent(inout) :: this
character(len=*), intent(in) :: label
type(alloc_infos) :: new_alloc
new_alloc%label = label
new_alloc%size = product(shape(list))
new_alloc%type = get_type(list)
call this%push_new_info(new_alloc)
