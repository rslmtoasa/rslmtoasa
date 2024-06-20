character(len=10) :: symbol
integer :: f_core, num_quant_s, num_quant_p, num_quant_d
real(rp) :: atomic_number, core, valence

namelist /element/ symbol, atomic_number, core, valence, f_core, num_quant_s, num_quant_p, num_quant_d
