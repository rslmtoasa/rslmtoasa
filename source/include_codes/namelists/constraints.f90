logical :: constraints_enable
integer :: constraints_i_cons, constraints_code_prefac
real(rp), dimension(:, :), allocatable :: constraints_mom_ref
real(rp), dimension(:, :), allocatable :: constraints_bfield

namelist /constraints/ constraints_enable, constraints_i_cons, constraints_code_prefac, constraints_mom_ref, constraints_bfield
