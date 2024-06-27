logical :: hoh, local_axis, orb_pol
real(rp), dimension(:,:), allocatable :: hubbard_u
character(len=3), dimension(:), allocatable :: hubbard_orb
namelist /hamiltonian/ hoh, local_axis, orb_pol, hubbard_u, hubbard_orb
