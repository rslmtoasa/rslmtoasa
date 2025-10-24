logical :: hoh, local_axis, orb_pol, v_a, v_b
character(len=20) :: js_alpha, jl_alpha
real(rp), dimension(3) :: v_alpha, v_beta, q_ss
real(rp), dimension(:), allocatable :: velocity_scale
real(rp), dimension(:,:), allocatable :: hubbard_u
character(len=10), dimension(:), allocatable :: uj_orb
real(rp), dimension(:,:), allocatable :: hubbard_j
real(rp), dimension(:,:,:,:), allocatable :: hubbard_v
real(rp), dimension(:,:), allocatable :: hubbard_u_impurity, hubbard_j_impurity
real(rp), dimension(:,:), allocatable :: hubbard_u_general, hubbard_j_general
integer, dimension(:,:), allocatable :: hubbard_u_sc
character(len=16) :: bounds_algorithm
real(rp) :: bounds_scaling

namelist /hamiltonian/ hoh, local_axis, orb_pol, v_alpha, v_beta, js_alpha, jl_alpha, velocity_scale, q_ss, uj_orb, hubbard_u, hubbard_j, hubbard_v, &
        hubbard_u_sc, hubbard_u_impurity, hubbard_j_impurity, hubbard_u_general, hubbard_j_general, bounds_algorithm, bounds_scaling
