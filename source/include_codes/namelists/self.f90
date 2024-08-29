logical :: ws_all, rigid_band, orbital_polarization, mixmag_all, mix_all, magnetic_mixing, freeze, all_inequivalent, fix_soc
integer :: nstep, init
integer, dimension(:), allocatable :: rb
real(rp) :: conv_thr, soc_scale
real(rp), dimension(:), allocatable :: ws, mixmag

namelist /self/ ws_all, all_inequivalent, &
   mix_all, magnetic_mixing, mixmag_all, freeze, orbital_polarization, &
   rigid_band, rb, nstep, init, soc_scale, &
   conv_thr, ws, mixmag, fix_soc
