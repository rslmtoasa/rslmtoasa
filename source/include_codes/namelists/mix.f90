integer :: var
real(rp) :: beta
real(rp), dimension(:), allocatable :: magbeta
character(len=7) :: mixtype
real(rp) :: ldm_beta

namelist /mix/ var, beta, mixtype, magbeta, ldm_beta
