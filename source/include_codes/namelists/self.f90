logical :: ws_all, rigid_band, orbital_polarization, mixmag_all, mix_all, magnetic_mixing, freeze, all_inequivalent, fix_soc
logical :: use_kspace !< Use k-space diagonalization instead of recursion for SCF
logical :: magnetic_scf_diagnostics !< Emit ordinary q=0 magnetic SCF feedback diagnostics.
logical :: magnetic_seed_enable !< Apply a temporary symmetry-breaking B_fsm seed.
logical :: magnetic_seed_active !< Runtime state; not intended for input.
integer :: nstep, init
integer :: magnetic_seed_steps !< Number of outer SCF iterations with the temporary seed.
integer, dimension(:), allocatable :: rb
real(rp) :: conv_thr, soc_scale
real(rp) :: magnetic_seed_field !< Temporary B_fsm seed in the existing Ry convention.
real(rp), dimension(:), allocatable :: ws, mixmag
logical :: cold !< Cold start: perform ASA to extract potential parameters.

namelist /self/ ws_all, all_inequivalent, &
   mix_all, magnetic_mixing, mixmag_all, freeze, orbital_polarization, &
   rigid_band, rb, nstep, init, soc_scale, use_kspace, &
   conv_thr, ws, mixmag, fix_soc , cold, magnetic_scf_diagnostics, &
   magnetic_seed_enable, magnetic_seed_steps, magnetic_seed_field
