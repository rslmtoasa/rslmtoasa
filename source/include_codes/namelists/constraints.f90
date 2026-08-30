logical :: constraints_enable
integer :: constraints_i_cons, constraints_code_prefac
logical :: constraints_diagnostics
real(rp) :: constraints_tolerance
real(rp), dimension(:, :), allocatable :: constraints_mom_ref
real(rp), dimension(:, :), allocatable :: constraints_bfield

namelist /constraints/ constraints_enable, constraints_i_cons, constraints_code_prefac, constraints_diagnostics, constraints_tolerance, &
   constraints_mom_ref, constraints_bfield
