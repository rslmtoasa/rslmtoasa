logical :: hoh, local_axis, orb_pol
real(rp), dimension(:,:), allocatable :: hubbard_u
character(len=3), dimension(:), allocatable :: hubbard_orb
real(rp), dimension(:,:), allocatable :: hubbard_j
real(rp), dimension(:,:), allocatable :: F0
real(rp), dimension(:,:), allocatable :: F2
real(rp), dimension(:,:), allocatable :: F4
namelist /hamiltonian/ hoh, local_axis, orb_pol, hubbard_u, hubbard_j, hubbard_orb, F0, F2, F4
