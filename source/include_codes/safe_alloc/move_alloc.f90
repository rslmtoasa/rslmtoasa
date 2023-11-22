character(len=*), intent(in), optional :: label_src, label_dst
class(safe_alloc), intent(inout) :: this
type(alloc_infos) :: new_alloc
call move_alloc(src, dst)
if(present(label_src)) then
    new_alloc%label = label_src
    new_alloc%type = get_type(src)
    new_alloc%size = -product(shape(src))
    call this%push_new_info(new_alloc)
endif
if(present(label_dst)) then
    new_alloc%label = label_dst
    new_alloc%type = get_type(dst)
    new_alloc%size = product(shape(dst))
    call this%push_new_info(new_alloc)
endif