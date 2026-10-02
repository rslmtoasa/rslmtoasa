program test_linear_response_config
   use linear_response_mod, only: linear_response, tddft_capability_state, require_tddft_capability
   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use radial_ground_state_mod, only: radial_ground_state
   use linear_response_mod, only: linear_response_config, response_space_layout, lr_electronic_state, &
      tddft_production_result, evaluate_linear_response_sweep
   implicit none

   character(len=64) :: argument
   call get_command_argument(1, argument)
   select case (trim(argument))
   case ('bad-row')
      call check_case('bad-row', 'tddft', 'radial_points', 'lehmann', 'auto', 'stoner_fit', 'spd', 'none', 1)
   case ('bad-bare')
      call check_case('bad-bare', 'tddft', 'radial_points', 'kspace_resolvent', 'auto', 'alsda', 'spd', 'none', 1)
   case ('bad-solver')
      call check_case('bad-solver', 'tddft', 'radial_points', 'lehmann', 'block', 'alsda', 'spd', 'none', 1)
   case ('bad-diagnostics')
      call check_case('bad-diagnostics', 'tddft', 'radial_points', 'lehmann', 'auto', 'alsda', 'spd', 'bad', 1)
   case ('bad-projection')
      call check_case('bad-projection', 'projected', 'radial_points', 'lehmann', 'auto', 'none', 'bad', 'none', 1)
   case ('bad-channel')
      call check_case('bad-channel', 'tddft', 'radial_points', 'lehmann', 'auto', 'alsda', 'spd', 'none', 1, 'bad_channel')
   case ('bad-eta')
      call check_case('bad-eta', 'tddft', 'radial_points', 'lehmann', 'auto', 'alsda', 'spd', 'none', 1, 'chi_plus', -1.0)
   case ('bad-invariants-gamma')
      call check_case('bad-invariants-gamma', 'tddft', 'product_compact', 'lehmann', 'auto', 'alsda', 'spd', 'invariants', 1, &
         q_variant='no_gamma')
   case ('bad-invariants-pair')
      call check_case('bad-invariants-pair', 'tddft', 'product_compact', 'lehmann', 'auto', 'alsda', 'spd', 'invariants', 1, &
         q_variant='no_pair')
   case ('capability-soc')
      block
         type(tddft_capability_state) :: capability
         capability%has_soc = .true.
         call require_tddft_capability(capability)
      end block
   case ('capability-overlap')
      block
         type(tddft_capability_state) :: capability
         capability%generalized_overlap = .true.
         capability%orthogonal = .false.
         call require_tddft_capability(capability)
      end block
   case ('radial-lcmm-prepared')
      block
         type(linear_response_config) :: config
         type(response_space_layout) :: space
         type(lmto_radial_basis) :: radial(1)
         type(radial_ground_state) :: ground(1)
         type(lr_electronic_state) :: left, endpoints(1)
         type(tddft_production_result) :: result
         config%formulation = 'tddft'
         config%representation = 'radial_points'
         config%interaction = 'lcmm'
         config%response_lmax = 0
         ! Reject before any uninitialized numerical object is consumed.
         call evaluate_linear_response_sweep(config, space, radial, ground, left, endpoints, result)
         error stop 'uncertified prepared request executed'
      end block
   case ('radial-lcmm-lehmann')
      call check_case('radial-lcmm-lehmann', 'tddft', 'radial_points', 'lehmann', 'auto', 'lcmm', 'spd', 'none', 1)
   case ('radial-lcmm-gf')
      call check_case('radial-lcmm-gf', 'tddft', 'radial_points', 'realspace_gf', 'block', 'lcmm', 'spd', 'none', 1)
   case ('juelich-spd')
      call check_case('juelich-spd', 'projected', 'radial_points', 'lehmann', 'auto', 'lcmm', 'spd', 'none', 4)
   case ('juelich-spdf')
      call check_case('juelich-spdf', 'projected', 'radial_points', 'lehmann', 'auto', 'lcmm', 'spdf', 'none', 4)
   case ('juelich-both')
      call check_case('juelich-both', 'projected', 'radial_points', 'lehmann', 'auto', 'lcmm', 'both', 'none', 4)
   case ('juelich-minus')
      call check_case('juelich-minus', 'projected', 'radial_points', 'lehmann', 'auto', 'lcmm', 'd', 'none', 4, 'chi_minus')
   case ('mills-spd')
      call check_case('mills-spd', 'projected', 'radial_points', 'lehmann', 'auto', 'mills_1u', 'spd', 'none', 4, 'chi_plus')
   case ('kl')
      call check_case('kl', 'kl', 'radial_points', 'lehmann', 'auto', 'alsda', 'spd', 'none', 4, 'chi_plus')
   case ('kl_dynamic')
      call check_case('kl_dynamic', 'kl_dynamic', 'radial_points', 'lehmann', 'auto', 'alsda', 'spd', 'none', 4, 'chi_plus')
   case ('halle')
      call check_case('halle', 'halle', 'radial_points', 'lehmann', 'auto', 'alsda', 'spd', 'none', 4, 'chi_plus')
   case ('longitudinal')
      call check_case('longitudinal', 'tddft', 'radial_points', 'lehmann', 'auto', 'alsda', 'spd', 'none', 4, 'chi_z')
   case ('charge')
      call check_case('charge', 'tddft', 'radial_points', 'lehmann', 'auto', 'alsda', 'spd', 'none', 4, 'charge')
   case ('native-alias')
      call check_case('native-alias', 'rotation', 'radial_points', 'lehmann', 'auto', 'alsda', 'spd', 'none', 1, &
         extra="native_crosscheck=.true.")
   case ('provider-alias')
      call check_case('provider-alias', 'rotation', 'radial_points', 'lehmann', 'auto', 'alsda', 'spd', 'none', 1, &
         extra="native_rsgf_provider='block_recursion'")
   case ('occupied')
      call check_case('occupied', 'rotation', 'radial_points', 'lehmann', 'auto', 'alsda', 'spd', 'none', 1, &
         extra="finite_h_spectral_mode='legacy_occupied'")
   case ('contour-only')
      call check_case('contour-only', 'rotation', 'radial_points', 'lehmann', 'auto', 'alsda', 'spd', 'none', 1, &
         extra="finite_h_response_backend='contour'")
   case ('radial-full')
      call check_case('radial-full', 'tddft', 'radial_points', 'lehmann', 'auto', 'alsda', 'spd', 'none', 1, &
         extra='response_lmax=2')
   case ('juelich-eta')
      call check_case('juelich-eta', 'projected', 'radial_points', 'lehmann', 'auto', 'lcmm', 'd', 'none', 1)
   case ('juelich-no-gamma')
      call check_case('juelich-no-gamma', 'projected', 'radial_points', 'lehmann', 'auto', 'lcmm', 'd', 'none', 4, &
         q_variant='no_gamma')
   case ('capability-noncollinear')
      block
         type(tddft_capability_state) :: capability
         capability%collinear = .false.
         call require_tddft_capability(capability)
      end block
   case ('mills_fit', 'projected_mills', 'projected_stoner', 'scalarized_mills')
      call check_case(trim(argument), 'projected', 'radial_points', 'lehmann', 'auto', trim(argument), 'd', 'none', 1)
   case ('row_rotation')
      call check_case('row_rotation', 'rotation', 'radial_points', 'lehmann', 'auto', 'alsda', 'spd', 'none', 1)
   case ('row_radial_lehmann_alsda')
      call check_case('row_radial_lehmann_alsda', 'tddft', 'radial_points', 'lehmann', 'auto', 'alsda', 'spd', 'none', 1)
   case ('row_radial_gf_alsda')
      call check_case('row_radial_gf_alsda', 'tddft', 'radial_points', 'realspace_gf', 'block', 'alsda', 'spd', 'none', 1)
   case ('row_product_alsda')
      call check_case('row_product_alsda', 'tddft', 'product_compact', 'lehmann', 'auto', 'alsda', 'spd', 'invariants', 1)
   case ('row_product_none')
      call check_case('row_product_none', 'tddft', 'product_compact', 'lehmann', 'auto', 'none', 'spd', 'none', 1)
   case ('row_projected_none')
      call check_case('row_projected_none', 'projected', 'radial_points', 'lehmann', 'auto', 'none', 'd', 'none', 1)
   case ('row_projected_mills_1u')
      call check_case('row_projected_mills_1u', 'projected', 'radial_points', 'lehmann', 'auto', 'mills_1u', 'd', 'none', 1)
   case ('row_projected_lcmm', '')
      call check_case('row_projected_lcmm', 'projected', 'radial_points', 'lehmann', 'auto', 'lcmm', 'd', 'none', 4)
   case default
      error stop 'unknown configuration test case'
   end select

contains

   subroutine check_case(name, formulation, representation, bare_response, solver, interaction, projection, diagnostics, n_eta, &
                         channel, eta, q_variant, extra)
      character(len=*), intent(in) :: name, formulation, representation, bare_response, solver, interaction, projection, diagnostics
      integer, intent(in) :: n_eta
      character(len=*), intent(in), optional :: channel
      real, intent(in), optional :: eta
      character(len=*), intent(in), optional :: q_variant, extra
      type(linear_response) :: response
      character(len=32) :: selected_channel
      real :: selected_eta

      selected_channel = 'chi_plus'
      if (present(channel)) selected_channel = channel
      selected_eta = 0.01
      if (present(eta)) selected_eta = eta
      call write_case(name, formulation, representation, bare_response, solver, interaction, projection, diagnostics, &
         selected_channel, selected_eta, n_eta, q_variant, extra)
      call response%load_config('lr_config_'//trim(name)//'.nml', .true.)
      if (trim(response%config%formulation) /= trim(formulation)) error stop 'formulation did not parse'
      if (trim(response%config%interaction) /= trim(interaction)) error stop 'interaction did not parse'
      if (size(response%config%q_list, 2) /= 3) error stop 'q list did not expand'
   end subroutine check_case

   subroutine load_case(name)
      character(len=*), intent(in) :: name
      type(linear_response) :: response
      call response%load_config('lr_config_'//trim(name)//'.nml', .true.)
   end subroutine load_case

   subroutine write_case(name, formulation, representation, bare_response, solver, interaction, projection, diagnostics, &
                         channel, eta, n_eta, q_variant, extra)
      character(len=*), intent(in) :: name, formulation, representation, bare_response, solver, interaction, projection, diagnostics
      character(len=*), intent(in), optional :: channel
      real, intent(in), optional :: eta
      integer, intent(in), optional :: n_eta
      character(len=*), intent(in), optional :: q_variant, extra
      integer :: unit
      character(len=32) :: selected_channel
      real :: selected_eta
      integer :: selected_n_eta

      selected_channel = 'chi_plus'
      if (present(channel)) selected_channel = channel
      selected_eta = 0.01
      if (present(eta)) selected_eta = eta
      selected_n_eta = 1
      if (present(n_eta)) selected_n_eta = n_eta
      open(newunit=unit, file='lr_config_'//trim(name)//'.nml', status='replace', action='write')
      write(unit, '(a)') '&linear_response'
      write(unit, '(a,a,a)') " formulation = '", trim(formulation), "'"
      write(unit, '(a,a,a)') " representation = '", trim(representation), "'"
      write(unit, '(a,a,a)') " bare_response = '", trim(bare_response), "'"
      write(unit, '(a,a,a)') " realspace_solver = '", trim(solver), "'"
      write(unit, '(a,a,a)') " interaction = '", trim(interaction), "'"
      write(unit, '(a,a,a)') " projection = '", trim(projection), "'"
      write(unit, '(a,a,a)') " diagnostics = '", trim(diagnostics), "'"
      write(unit, '(a,a,a)') " channel = '", trim(selected_channel), "'"
      if (present(q_variant)) then
         if (trim(q_variant) == 'no_gamma') then
            write(unit, '(a)') ' n_q = 3, q_list = 0.03, 0.0, 0.0, 0.06, 0.0, 0.0, 0.09, 0.0, 0.0'
         else if (trim(q_variant) == 'no_pair') then
            write(unit, '(a)') ' n_q = 3, q_list = 0.0, 0.0, 0.0, 0.03, 0.0, 0.0, 0.06, 0.0, 0.0'
         else
            write(unit, '(a)') ' n_q = 3, q_list = 0.0, 0.0, 0.0, 0.03, 0.0, 0.0, -0.03, 0.0, 0.0'
         end if
      else
         write(unit, '(a)') ' n_q = 3, q_list = 0.0, 0.0, 0.0, 0.03, 0.0, 0.0, -0.03, 0.0, 0.0'
      end if
      if (trim(formulation) == 'rotation') write(unit, '(a)') ' n_q_points=3'
      write(unit, '(a)') ' n_omega = 2, omega_min = 0.0, omega_max = 0.1'
      if (selected_n_eta > 1) then
         write(unit, '(a)') ' eta = 0.01, n_eta = 4, eta_grid = 0.04, 0.02, 0.01, 0.005'
      else
         write(unit, '(a)') ' eta = 0.01, n_eta = 1, eta_grid = 0.01'
      end if
      if (selected_eta /= 0.01) write(unit, '(a,es16.8)') ' eta = ', selected_eta
      write(unit, '(a)') " response_lmax = -1, write_full_matrix = .false., output_file = 'lr_config.dat'"
      if (trim(formulation) == 'tddft' .and. trim(representation) == 'radial_points') &
         write(unit, '(a)') ' response_lmax=0'
      write(unit, '(a)') ' gf_integration_points = 3, gf_energy_margin = 1.0'
      if (present(extra)) write(unit, '(a)') extra
      write(unit, '(a)') '/'
      close(unit)
   end subroutine write_case

end program test_linear_response_config
