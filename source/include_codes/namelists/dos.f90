! ============================================================================
! DEPRECATED NAMELIST - DO NOT USE
! ============================================================================
! This &dos namelist is NO LONGER READ by the code!
! Use &reciprocal namelist instead for all k-space post-processing tasks.
!
! Migration: All parameters from &dos have been moved to &reciprocal.
! See NAMELIST_HOMOGENIZATION.md for migration guide.
! ============================================================================

character(len=30) :: dos_method
real :: gaussian_sigma
integer :: n_energy_points
real :: dos_energy_min
real :: dos_energy_max
real :: temperature
integer :: nk1, nk2, nk3
real :: total_electrons
logical :: auto_find_fermi

namelist /dos/ dos_method, gaussian_sigma, n_energy_points, dos_energy_min, dos_energy_max, temperature, nk1, nk2, nk3, total_electrons, auto_find_fermi