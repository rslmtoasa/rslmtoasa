logical :: fix_fermi
integer :: channels_ldos
real(rp) :: fermi, energy_min, energy_max

namelist /energy/ fix_fermi, channels_ldos, fermi, energy_min, energy_max