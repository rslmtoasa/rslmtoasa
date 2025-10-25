! Namelist for spectrum bounds options
! Declares variables read by the bounds namelist
character(len=16) :: bounds_algorithm
real(rp) :: bounds_scaling

namelist /bounds/ bounds_algorithm, bounds_scaling
