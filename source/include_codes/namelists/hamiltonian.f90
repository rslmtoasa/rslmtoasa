logical :: hoh, local_axis, orb_pol, v_a, v_b
character(len=10) :: js_alpha, jl_alpha
real(rp), dimension(3) :: v_alpha, v_beta
namelist /hamiltonian/ hoh, local_axis, orb_pol, v_alpha, v_beta, js_alpha, jl_alpha
