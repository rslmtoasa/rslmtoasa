logical :: hoh, local_axis, orb_pol, v_a, v_b
character(len=10) :: pol_alpha, pol_beta
real(rp), dimension(3) :: v_alpha, v_beta
real(rp), dimension(:), allocatable :: velocity_scale
namelist /hamiltonian/ hoh, local_axis, orb_pol, v_alpha, v_beta, pol_alpha, pol_beta, velocity_scale
