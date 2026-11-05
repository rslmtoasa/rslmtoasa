logical :: hoh, local_axis, orb_pol, v_a, v_b
character(len=20) :: js_alpha, jl_alpha
real(rp), dimension(3) :: v_alpha, v_beta
! Generalized Bloch theorem spin-spiral input units in this namelist:
! - q_ss is given in Cartesian units of 2*pi/alat (q_ss=0.5 along a simple
!   cubic direction is the zone boundary at pi/alat); theta_ss (cone angle
!   from the z axis) is given in degrees.
! Internal Hamiltonian storage (hamiltonian%q_ss, hamiltonian%theta_ss) uses the
! same q_ss phase convention and radians, respectively -- see build_from_file.
real(rp), dimension(3) :: q_ss
real(rp) :: theta_ss
character(len=24) :: magnetic_representation
logical :: ccor_2c, ccor_debug, ccor_strict
real(rp) :: ccor_elin
character(len=16) :: ccor_vmt_mode
real(rp), dimension(:), allocatable :: velocity_scale
! Hubbard input units in this namelist:
! - hubbard_u / hubbard_j / hubbard_u_impurity / hubbard_j_impurity /
!   hubbard_u_general / hubbard_j_general / hubbard_v are read in eV.
! Internal Hamiltonian/potential storage uses Ry.
real(rp), dimension(:,:), allocatable :: hubbard_u
character(len=10), dimension(:), allocatable :: uj_orb
real(rp), dimension(:,:), allocatable :: hubbard_j
real(rp), dimension(:,:,:,:), allocatable :: hubbard_v
real(rp), dimension(:,:), allocatable :: hubbard_u_impurity, hubbard_j_impurity
real(rp), dimension(:,:), allocatable :: hubbard_u_general, hubbard_j_general
integer, dimension(:,:), allocatable :: hubbard_u_sc
character(len=16) :: bounds_algorithm
logical :: hubbard_experimental = .false. ! Opt in to experimental ACBN0 and +V paths.
character(len=16) :: hubbard_log_level = 'summary'
logical :: hubbard_spinor_dm_diagnostic = .false. ! Extract the spinor d-shell P and M1 at nsp>=2 (diagnostic).
logical :: hubbard_soc_experimental = .false. ! Experimental: spinor fixed U (nsp=2, bulk, spectral radial moments).
character(len=16) :: hubbard_u_potential_form
! Omitted input uses an internal auto sentinel, resolved after backend selection.
character(len=16) :: hubbard_radial_moments
!> Explicit A4 test-only onsite d-shell operator selector and scale in Ry.
character(len=16) :: hubbard_test_operator
real(rp) :: hubbard_test_lambda
!> Test-only W probe spin: 0 both, 1 up, 2 down.
integer :: hubbard_test_spin
logical :: fixed_u_test_tighten_r05
!> A5 test-only HOH control; valid only for explicit U=J=0, local W1 runs.
logical :: fixed_u_test_allow_hoh_w1_u0
real(rp) :: bounds_scaling
character(len=16) :: export

namelist /hamiltonian/ hoh, local_axis, orb_pol, v_alpha, v_beta, js_alpha, jl_alpha, velocity_scale, q_ss, theta_ss, magnetic_representation, ccor_2c, ccor_elin, &
        ccor_vmt_mode, ccor_debug, ccor_strict, uj_orb, hubbard_u, hubbard_j, hubbard_v, &
        hubbard_u_sc, hubbard_u_impurity, hubbard_j_impurity, hubbard_u_general, hubbard_j_general, bounds_algorithm, bounds_scaling, &
        hubbard_experimental, hubbard_log_level, hubbard_spinor_dm_diagnostic, hubbard_soc_experimental, hubbard_u_potential_form, hubbard_radial_moments, hubbard_test_operator, hubbard_test_lambda, hubbard_test_spin, &
        fixed_u_test_tighten_r05, fixed_u_test_allow_hoh_w1_u0, export
