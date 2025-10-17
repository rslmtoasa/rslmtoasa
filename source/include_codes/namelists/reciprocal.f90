! Reciprocal space k-point mesh settings
integer :: nk1, nk2, nk3
real :: k_offset_x, k_offset_y, k_offset_z
logical :: use_symmetry_reduction, use_time_reversal, use_shift

! Density of states settings
! NOTE: All energy values in Rydberg (Ry) for consistency with Hamiltonian
!       dos_energy_min, dos_energy_max: Energy window in Ry (e.g., -1.0, 1.0)
!       gaussian_sigma: Gaussian smearing width in Ry (e.g., 0.01)
!       temperature: Temperature in Kelvin for Fermi-Dirac distribution
!       Example for BCC Fe: dos_energy_min = -0.8, dos_energy_max = 0.9 (Ry)
integer :: n_energy_points
real :: dos_energy_min, dos_energy_max    ! Energy range in Ry
real :: gaussian_sigma                     ! Gaussian smearing in Ry
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
