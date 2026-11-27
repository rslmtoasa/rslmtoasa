character(len=30) :: pre_processing
character(len=30) :: processing
character(len=30) :: post_processing
character(len=30) :: gf_route, exchange_delta
logical :: verbose
logical :: do_damping
logical :: do_inertia
logical :: do_chi0
logical :: do_sigma
real(rp) :: sigma_tol

namelist /calculation/ verbose, pre_processing, processing, post_processing, gf_route, exchange_delta, do_damping, do_inertia, do_chi0, do_sigma, sigma_tol
