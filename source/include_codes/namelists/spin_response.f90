
character(len=16) :: method, kernel
character(len=8) :: channels
character(len=256) :: output_prefix, q_file
integer :: n_q, n_omega
real(rp) :: q_list(3, max_q), omega_min, omega_max, eta

! kernel = 'goldstone' (default, production: Juelich-d, Goldstone-fitted U); 'alsda' and 'alsda_scaled' are experimental (LDA only).
namelist /spin_response/ method, kernel, channels, output_prefix, n_q, q_list, q_file, omega_min, omega_max, n_omega, eta
