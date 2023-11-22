integer :: var
real(rp) :: beta
real(rp), dimension(:), allocatable :: magbeta
character(len=7) :: mixtype

namelist /mix/ var, beta, mixtype, magbeta
