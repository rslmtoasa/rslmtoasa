! ============================================================================
! &reciprocal NAMELIST - K-space calculations and post-processing
! ============================================================================
! This namelist is used for:
!   - post_processing = 'band_structure'
!   - post_processing = 'density_of_states'
!   - K-space SCF (calctype = 'K')
!
! NOTE: Replaces deprecated &dos namelist (see NAMELIST_HOMOGENIZATION.md)
! ============================================================================

! K-point mesh settings
integer :: nk1, nk2, nk3                   ! K-mesh dimensions (e.g., 20, 20, 20)
real :: k_offset_x, k_offset_y, k_offset_z ! K-mesh offsets (0.0 = Gamma-centered)
logical :: use_symmetry_reduction          ! Use crystal symmetry (recommended: .true.)
logical :: use_time_reversal               ! Use time-reversal symmetry (recommended: .true.)
logical :: use_shift                       ! Additional Monkhorst-Pack shift (default: .false.)

! Density of states settings (for post_processing = 'density_of_states')
! **CRITICAL**: All energy values are in RYDBERGS (Ry), NOT eV!
!               Conversion: 1 Ry = 13.6057 eV,  1 eV = 0.07350 Ry
!
! Typical values for 3d transition metals (e.g., Fe, Co, Ni):
!   dos_energy_min = -0.8 Ry   (≈ -11 eV below Fermi level)
!   dos_energy_max =  0.9 Ry   (≈ +12 eV above Fermi level)
!   gaussian_sigma =  0.02 Ry  (≈ 0.27 eV smearing, for 'gaussian' method)
!
! DOS methods:
!   'tetrahedron' - Linear tetrahedron method (accurate, needs dense k-mesh)
!   'blochl'      - Blöchl improved tetrahedron (recommended, more stable)
!   'gaussian'    - Gaussian broadening (smoother, good for quick tests)
!
integer :: n_energy_points                 ! Number of DOS energy grid points (e.g., 2000)
real :: dos_energy_min, dos_energy_max     ! Energy window in Ry (NOT eV!)
real :: gaussian_sigma                     ! Gaussian smearing in Ry (only for 'gaussian')
real :: temperature                        ! Temperature in Kelvin (e.g., 300.0)
real :: total_electrons                    ! Total valence electrons (for auto_find_fermi)
character(len=20) :: dos_method            ! 'tetrahedron', 'blochl', or 'gaussian'
logical :: auto_find_fermi                 ! Find Fermi from charge conservation
logical :: suppress_internal_logs          ! Reduce verbosity (.true. recommended)

namelist /reciprocal/ nk1, nk2, nk3, k_offset_x, k_offset_y, k_offset_z, &
   use_symmetry_reduction, use_time_reversal, use_shift, &
   n_energy_points, dos_energy_min, dos_energy_max, &
   gaussian_sigma, temperature, total_electrons, dos_method, &
   auto_find_fermi, suppress_internal_logs
