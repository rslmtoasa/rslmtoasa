
character(len=16) :: method
character(len=256) :: output_prefix, q_file
integer :: n_q, n_omega
real(rp) :: q_list(3, max_q), omega_min, omega_max, eta

namelist /spin_response/ method, output_prefix, n_q, q_list, q_file, omega_min, omega_max, n_omega, eta
