real(rp) :: t_i, t_f, dt, alpha, sd_temp
real(rp), dimension(:,:), allocatable :: b_stochastic
integer :: nt, asd_step
character(len=5) :: integrator
namelist /sd/ t_i, t_f, dt, alpha, b_stochastic, nt, integrator, asd_step, sd_temp
