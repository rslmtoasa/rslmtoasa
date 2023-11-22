class(safe_alloc), intent(inout) :: this
character(len=*), intent(in) :: label
call this%allocate(label, list, shape(mold))