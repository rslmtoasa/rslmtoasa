! ============================================================================
! Minimal production TD-DFT namelist.
! This is intentionally not the removed legacy TD-DFT object/parser.
! ============================================================================
integer, parameter :: tddft_max_q = 128
integer, parameter :: tddft_max_omega = 4096

logical :: enabled
character(len=32) :: channel
integer :: n_q
real(rp) :: q_list(3, tddft_max_q)
integer :: n_omega
logical :: use_omega_grid
real(rp) :: omega_grid(tddft_max_omega)
real(rp) :: omega_min, omega_max
real(rp) :: eta
integer :: n_eta
real(rp) :: eta_grid(16)
integer :: response_lmax
character(len=48) :: interaction_route
logical :: goldstone_correction
character(len=32) :: backend
logical :: reciprocal_backend_crosscheck
character(len=32) :: native_rsgf_provider
integer :: gf_integration_points
real(rp) :: gf_integration_eta
real(rp) :: gf_energy_margin
logical :: gf_closure_audit
logical :: dyson_static_audit
logical :: validate_interacting_covariance
logical :: write_full_matrix
character(len=256) :: output_file

namelist /tddft/ enabled, channel, n_q, q_list, n_omega, use_omega_grid, omega_grid, &
   omega_min, omega_max, eta, n_eta, eta_grid, response_lmax, interaction_route, goldstone_correction, &
   backend, reciprocal_backend_crosscheck, native_rsgf_provider, gf_integration_points, gf_integration_eta, gf_energy_margin, &
   gf_closure_audit, dyson_static_audit, validate_interacting_covariance, write_full_matrix, output_file
