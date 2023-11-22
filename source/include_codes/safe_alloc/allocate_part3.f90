new_alloc%label = label
new_alloc%size = product(list_size)
new_alloc%type = get_type(list)
call this%push_new_info(new_alloc)
