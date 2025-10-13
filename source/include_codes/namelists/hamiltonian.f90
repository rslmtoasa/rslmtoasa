logical :: hoh, local_axis, orb_pol, v_a, v_b
character(len=20) :: js_alpha, jl_alpha
real(rp), dimension(3) :: v_alpha, v_beta, q_ss
real(rp), dimension(:), allocatable :: velocity_scale
namelist /hamiltonian/ hoh, local_axis, orb_pol, v_alpha, v_beta, js_alpha, jl_alpha, velocity_scale, q_ss
