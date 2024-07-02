logical :: hoh, local_axis, orb_pol
real(rp), dimension(:,:), allocatable :: hubbard_u
character(len=3), dimension(:), allocatable :: hubbard_orb
real(rp), dimension(:,:), allocatable :: hubbard_j
namelist /hamiltonian/ hoh, local_axis, orb_pol, hubbard_orb, hubbard_u, hubbard_j
