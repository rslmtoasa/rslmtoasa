! &tddft: finite-q, finite-frequency linear-response configuration.
! q coordinates are fractional reciprocal-lattice coordinates.  q_list is
! used by q_mode='list'; q_mode='path' interpolates q_start to q_end; and
! q_mode='mesh' builds a Gamma-centred uniform mesh.
integer, parameter :: tddft_max_q = 2000
logical :: enabled
character(len=16) :: channel, chi0_backend, response_projection, q_mode, q_coordinates
character(len=16) :: goldstone_mode
character(len=256) :: q_file, output_prefix
integer :: n_q_points, q_mesh1, q_mesh2, q_mesh3
integer :: nomega, band_first, band_last
real(rp) :: q_start(3), q_end(3), q_list(3, tddft_max_q)
real(rp) :: omega_min, omega_max, eta, electronic_temperature, fermi_level, occupation_tolerance
logical :: output_chi0, output_xi, output_chi, output_modes, output_stoner

namelist /tddft/ enabled, channel, chi0_backend, response_projection, q_mode, q_coordinates, goldstone_mode, &
   q_file, output_prefix, n_q_points, q_mesh1, q_mesh2, q_mesh3, q_start, q_end, q_list, &
   omega_min, omega_max, nomega, eta, electronic_temperature, fermi_level, band_first, band_last, &
   occupation_tolerance, output_chi0, output_xi, output_chi, output_modes, output_stoner
