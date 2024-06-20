logical :: ws_all, rigid_band, orbital_polarization, mixmag_all, mix_all, magnetic_mixing, freeze, all_inequivalent
integer :: nstep, init
integer, dimension(:), allocatable :: rb
real(rp) :: conv_thr
real(rp), dimension(:), allocatable :: ws, mixmag

namelist /self/ ws_all, all_inequivalent, &
   mix_all, magnetic_mixing, mixmag_all, freeze, orbital_polarization, &
   rigid_band, rb, nstep, init, &
   conv_thr, ws, mixmag
