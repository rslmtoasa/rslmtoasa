real(rp), allocatable :: lambda(:, :, :)
character(len=256) :: delta_file
integer :: max_iter, n_dos
real(rp) :: delta_seed, mix, tol, dos_window, dos_sigma

namelist /bdg/ lambda, delta_seed, delta_file, mix, tol, max_iter, dos_window, n_dos, dos_sigma
