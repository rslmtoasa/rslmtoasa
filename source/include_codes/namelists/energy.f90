logical :: fix_fermi
integer :: channels_ldos
real(rp) :: fermi, energy_min, energy_max
!> Electronic temperature (K) of the exchange energy integration; 0 = step function at E_F
real(rp) :: temperature

namelist /energy/ fix_fermi, channels_ldos, fermi, energy_min, energy_max, temperature
