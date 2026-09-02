! &tddft: finite-q, finite-frequency linear-response configuration.
! q coordinates are fractional reciprocal-lattice coordinates.  q_list is
! used by q_mode='list'; q_mode='path' interpolates q_start to q_end; and
! q_mode='mesh' builds a Gamma-centred uniform mesh.
integer, parameter :: tddft_max_q = 2000
logical :: enabled
character(len=16) :: channel, chi0_backend, response_projection, q_mode, q_coordinates, circular_channel
character(len=24) :: gf_integration
character(len=24) :: xi_backend
character(len=16) :: goldstone_mode, goldstone_policy
character(len=256) :: q_file, output_prefix, longitudinal_static_file
integer :: n_q_points, q_mesh1, q_mesh2, q_mesh3
integer :: nomega, band_first, band_last, green_energy_points, contour_points, contour_subdivisions, near_fermi_points
real(rp) :: q_start(3), q_end(3), q_list(3, tddft_max_q)
real(rp) :: omega_min, omega_max, eta, electronic_temperature, occupation_tolerance
real(rp) :: green_eta, green_energy_min, green_energy_max
real(rp) :: contour_height
real(rp) :: realspace_rmax, realspace_source_rmax, realspace_tail_tolerance
integer :: realspace_fourier_axes(3)
character(len=16) :: realspace_representation, realspace_truncation_mode
real(rp) :: longitudinal_pair_tolerance, longitudinal_linearity_tolerance, longitudinal_static_agreement_tolerance
real(rp) :: longitudinal_fit_omega_min, longitudinal_fit_omega_max
logical :: output_chi0, output_xi, output_chi, output_modes, output_stoner

namelist /tddft/ enabled, channel, chi0_backend, xi_backend, response_projection, circular_channel, q_mode, q_coordinates, goldstone_mode, goldstone_policy, &
   gf_integration, contour_points, contour_subdivisions, near_fermi_points, contour_height, &
   q_file, output_prefix, n_q_points, q_mesh1, q_mesh2, q_mesh3, q_start, q_end, q_list, &
   omega_min, omega_max, nomega, eta, electronic_temperature, band_first, band_last, &
   occupation_tolerance, green_eta, green_energy_min, green_energy_max, green_energy_points, &
   realspace_rmax, realspace_source_rmax, realspace_tail_tolerance, realspace_fourier_axes, realspace_representation, &
   realspace_truncation_mode, &
   longitudinal_static_file, longitudinal_pair_tolerance, longitudinal_linearity_tolerance, &
   longitudinal_static_agreement_tolerance, longitudinal_fit_omega_min, longitudinal_fit_omega_max, &
   output_chi0, output_xi, output_chi, output_modes, output_stoner
