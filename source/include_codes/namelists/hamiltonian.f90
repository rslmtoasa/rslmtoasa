logical :: hoh, local_axis, orb_pol
real(rp), dimension(:,:), allocatable :: hubbard_u
character(len=10), dimension(:), allocatable :: uj_orb
real(rp), dimension(:,:), allocatable :: hubbard_j
real(rp), dimension(:,:,:,:), allocatable :: hubbard_v
namelist /hamiltonian/ hoh, local_axis, orb_pol, uj_orb, hubbard_u, hubbard_j, hubbard_v
