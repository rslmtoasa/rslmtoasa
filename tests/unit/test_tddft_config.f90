!------------------------------------------------------------------------------
! TDDFT-06 -- focused input-boundary tests for &tddft configuration.
!------------------------------------------------------------------------------
program test_tddft_config
   use precision_mod, only: rp
   use tddft_config_mod, only: tddft_config
   implicit none

   logical :: failed

   failed = .false.
   call test_path_input()
   call test_mesh_input()
   call test_longitudinal_input()
   call test_full_input()
   call test_green_input()
   call test_xi_backend_input()
   call test_ground_state_provenance_defaults()
   if (failed) error stop 1
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_path_input()
      type(tddft_config) :: config
      integer :: unit

      open(newunit=unit, file='unit_tddft_config.nml', status='replace', action='write')
      write(unit, '(a)') '&tddft'
      write(unit, '(a)') " channel = 'transverse', chi0_backend = 'eigenpairs', response_projection = 'site'"
      write(unit, '(a)') " q_mode = 'path', q_coordinates = 'direct', n_q_points = 3"
      write(unit, '(a)') ' q_start = 0.0, 0.0, 0.0, q_end = 0.5, 0.0, 0.0'
      write(unit, '(a)') ' omega_min = 0.0, omega_max = 0.2, nomega = 5, eta = 0.01'
      write(unit, '(a)') ' electronic_temperature = 150.0, band_first = 2, band_last = 8'
      write(unit, '(a)') " goldstone_mode = 'diagnose', output_prefix = 'unit_response'"
      write(unit, '(a)') '/'
      close(unit)

      config = tddft_config('unit_tddft_config.nml')
      call assert_true('path produces three q points', size(config%q_points, 2) == 3)
      call assert_real('path midpoint is deterministic', config%q_points(1, 2), 0.25_rp)
      call assert_true('band window is read', config%band_first == 2 .and. config%band_last == 8)
      call assert_true('default bare output remains enabled', config%output_chi0 .and. config%output_stoner)
   end subroutine test_path_input

   subroutine test_mesh_input()
      type(tddft_config) :: config
      integer :: unit

      open(newunit=unit, file='unit_tddft_config.nml', status='replace', action='write')
      write(unit, '(a)') '&tddft'
      write(unit, '(a)') " q_mode = 'mesh', q_mesh1 = 2, q_mesh2 = 1, q_mesh3 = 2"
      write(unit, '(a)') ' output_xi = .false., output_chi = .false., output_modes = .false.'
      write(unit, '(a)') '/'
      close(unit)

      config = tddft_config('unit_tddft_config.nml')
      call assert_true('mesh produces product-sized q list', size(config%q_points, 2) == 4)
      call assert_real('mesh uses canonical negative edge', config%q_points(1, 1), -0.25_rp)
      call assert_real('mesh uses canonical positive edge', config%q_points(3, 4), 0.25_rp)
      open(newunit=unit, file='unit_tddft_config.nml', status='old')
      close(unit, status='delete')
   end subroutine test_mesh_input

   subroutine test_longitudinal_input()
      type(tddft_config) :: config
      integer :: unit

      open(newunit=unit, file='unit_tddft_config.nml', status='replace', action='write')
      write(unit, '(a)') '&tddft'
      write(unit, '(a)') " channel = 'longitudinal', longitudinal_static_file = 'static_fields.dat'"
      write(unit, '(a)') ' longitudinal_pair_tolerance = 1.0e-9, longitudinal_linearity_tolerance = 0.02'
      write(unit, '(a)') ' longitudinal_static_agreement_tolerance = 0.03'
      write(unit, '(a)') ' longitudinal_fit_omega_min = 0.001, longitudinal_fit_omega_max = 0.04'
      write(unit, '(a)') '/'
      close(unit)

      config = tddft_config('unit_tddft_config.nml')
      call assert_true('longitudinal channel is accepted', trim(config%channel) == 'longitudinal')
      call assert_true('longitudinal static driver path is retained', trim(config%longitudinal_static_file) == 'static_fields.dat')
      call assert_real('longitudinal fit upper bound is read', config%longitudinal_fit_omega_max, 0.04_rp)
      open(newunit=unit, file='unit_tddft_config.nml', status='old')
      close(unit, status='delete')
   end subroutine test_longitudinal_input

   subroutine test_full_input()
      type(tddft_config) :: config
      integer :: unit

      open(newunit=unit, file='unit_tddft_config.nml', status='replace', action='write')
      write(unit, '(a)') '&tddft'
      write(unit, '(a)') " channel = 'full', q_mode = 'list', n_q_points = 1"
      write(unit, '(a)') '/'
      close(unit)
      config = tddft_config('unit_tddft_config.nml')
      call assert_true('full charge-spin channel is accepted', trim(config%channel) == 'full')
      open(newunit=unit, file='unit_tddft_config.nml', status='old')
      close(unit, status='delete')
   end subroutine test_full_input

   subroutine test_green_input()
      type(tddft_config) :: config
      integer :: unit

      open(newunit=unit, file='unit_tddft_config.nml', status='replace', action='write')
      write(unit, '(a)') '&tddft'
      write(unit, '(a)') " chi0_backend = 'green', green_eta = 0.003, green_energy_min = -1.2"
      write(unit, '(a)') ' green_energy_max = 0.8, green_energy_points = 4001'
      write(unit, '(a)') '/'
      close(unit)
      config = tddft_config('unit_tddft_config.nml')
      call assert_true('Green chi0 backend is accepted', trim(config%chi0_backend) == 'green')
      call assert_real('Green eta is read', config%green_eta, 0.003_rp)
      call assert_true('Green energy mesh is read', config%green_energy_points == 4001)
      open(newunit=unit, file='unit_tddft_config.nml', status='old')
      close(unit, status='delete')
   end subroutine test_green_input

   subroutine test_xi_backend_input()
      type(tddft_config) :: config
      integer :: unit

      open(newunit=unit, file='unit_tddft_config.nml', status='replace', action='write')
      write(unit, '(a)') '&tddft'
      write(unit, '(a)') " xi_backend = 'compare', output_xi = .true., output_chi = .true."
      write(unit, '(a)') '/'
      close(unit)
      config = tddft_config('unit_tddft_config.nml')
      call assert_true('pair-potential shadow comparison backend is accepted', trim(config%xi_backend) == 'compare')
      open(newunit=unit, file='unit_tddft_config.nml', status='old')
      close(unit, status='delete')
   end subroutine test_xi_backend_input

   subroutine test_ground_state_provenance_defaults()
      type(tddft_config) :: inherited, overridden
      integer :: unit

      open(newunit=unit, file='unit_tddft_config.nml', status='replace', action='write')
      write(unit, '(a)') '&tddft'
      write(unit, '(a)') " output_prefix = 'unit_response'"
      write(unit, '(a)') '/'
      close(unit)
      inherited = tddft_config('unit_tddft_config.nml')
      call assert_true('omitted response EF inherits reciprocal provenance', .not. inherited%fermi_level_overridden)
      call assert_true('omitted response temperature inherits reciprocal provenance', .not. inherited%electronic_temperature_overridden)

      open(newunit=unit, file='unit_tddft_config.nml', status='replace', action='write')
      write(unit, '(a)') '&tddft'
      write(unit, '(a)') ' fermi_level = -0.071, electronic_temperature = 420.0'
      write(unit, '(a)') '/'
      close(unit)
      overridden = tddft_config('unit_tddft_config.nml')
      call assert_true('explicit response EF is auditable', overridden%fermi_level_overridden)
      call assert_true('explicit response temperature is auditable', overridden%electronic_temperature_overridden)
      call assert_real('explicit response EF is retained', overridden%fermi_level, -0.071_rp)
      open(newunit=unit, file='unit_tddft_config.nml', status='old')
      close(unit, status='delete')
   end subroutine test_ground_state_provenance_defaults

   subroutine assert_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine assert_true

   subroutine assert_real(label, actual, expected)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: actual, expected
      if (abs(actual-expected) > 1.0e-12_rp) then
         write (*, '(a,1x,a,2(1x,es12.4))') 'FAIL', label, actual, expected
         failed = .true.
      end if
   end subroutine assert_real

end program test_tddft_config
