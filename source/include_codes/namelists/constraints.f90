logical :: constraints_enable
integer :: constraints_i_cons, constraints_code_prefac
logical :: constraints_diagnostics
real(rp) :: constraints_tolerance
! Fixed-size buffers: a namelist read cannot allocate, so the columns actually set are copied out afterwards.
integer, parameter :: constraints_max_sites = 1024
real(rp), dimension(3, constraints_max_sites) :: constraints_mom_ref
real(rp), dimension(3, constraints_max_sites) :: constraints_bfield

namelist /constraints/ constraints_enable, constraints_i_cons, constraints_code_prefac, constraints_diagnostics, constraints_tolerance, &
   constraints_mom_ref, constraints_bfield
