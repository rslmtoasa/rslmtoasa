! Reciprocal space k-point mesh settings
integer :: nk1, nk2, nk3
real :: k_offset_x, k_offset_y, k_offset_z
logical :: use_symmetry_reduction, use_time_reversal, use_shift

! Density of states settings
integer :: n_energy_points
real :: dos_energy_min, dos_energy_max
real :: gaussian_sigma
real :: temperature
real :: total_electrons
character(len=20) :: dos_method
logical :: auto_find_fermi
logical :: suppress_internal_logs

namelist /reciprocal/ nk1, nk2, nk3, k_offset_x, k_offset_y, k_offset_z, &
   use_symmetry_reduction, use_time_reversal, use_shift, &
   n_energy_points, dos_energy_min, dos_energy_max, &
   gaussian_sigma, temperature, total_electrons, dos_method, &
   auto_find_fermi, suppress_internal_logs
