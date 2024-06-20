class(safe_alloc), intent(inout) :: this
character(len=*), intent(in) :: label
integer, intent(in) :: list_size
call this%allocate(label, list, (/list_size/))
