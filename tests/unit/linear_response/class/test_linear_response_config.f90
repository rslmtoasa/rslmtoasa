program test_linear_response_config
   use linear_response_mod, only: linear_response, tddft_capability_state, require_tddft_capability
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
   case default
      call check_case('row_radial_lehmann_alsda', 'tddft', 'radial_points', 'lehmann', 'auto', 'alsda', 'spd', 'none', 1)
      call check_case('row_radial_lehmann_lcmm', 'tddft', 'radial_points', 'lehmann', 'auto', 'lcmm', 'spd', 'none', 1)
      call check_case('row_radial_gf_alsda', 'tddft', 'radial_points', 'realspace_gf', 'block', 'alsda', 'spd', 'none', 1)
      call check_case('row_radial_gf_lcmm', 'tddft', 'radial_points', 'realspace_gf', 'chebyshev', 'lcmm', 'spd', 'none', 1)
      call check_case('row_product_alsda', 'tddft', 'product_compact', 'lehmann', 'auto', 'alsda', 'spd', 'invariants', 1)
      call check_case('row_product_none', 'tddft', 'product_compact', 'lehmann', 'auto', 'none', 'spd', 'none', 1)
      call check_case('row_projected_none', 'projected', 'radial_points', 'lehmann', 'auto', 'none', 'd', 'none', 1)
      call check_case('row_projected_stoner', 'projected', 'radial_points', 'lehmann', 'auto', 'stoner_fit', 'both', 'none', 3)
      call check_case('row_projected_lcmm', 'projected', 'radial_points', 'lehmann', 'auto', 'lcmm', 'spd', 'invariants', 3)
   end select

contains

   subroutine check_case(name, formulation, representation, bare_response, solver, interaction, projection, diagnostics, n_eta, &
                         channel, eta, q_variant)
      character(len=*), intent(in) :: name, formulation, representation, bare_response, solver, interaction, projection, diagnostics
      integer, intent(in) :: n_eta
      character(len=*), intent(in), optional :: channel
      real, intent(in), optional :: eta
      character(len=*), intent(in), optional :: q_variant
      type(linear_response) :: response
      character(len=32) :: selected_channel
      real :: selected_eta

      selected_channel = 'chi_plus'
      if (present(channel)) selected_channel = channel
      selected_eta = 0.01
      if (present(eta)) selected_eta = eta
      call write_case(name, formulation, representation, bare_response, solver, interaction, projection, diagnostics, &
         selected_channel, selected_eta, n_eta, q_variant)
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
                         channel, eta, n_eta, q_variant)
      character(len=*), intent(in) :: name, formulation, representation, bare_response, solver, interaction, projection, diagnostics
      character(len=*), intent(in), optional :: channel
      real, intent(in), optional :: eta
      integer, intent(in), optional :: n_eta
      character(len=*), intent(in), optional :: q_variant
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
      write(unit, '(a)') ' n_omega = 2, omega_min = 0.0, omega_max = 0.1'
      if (selected_n_eta > 1) then
         write(unit, '(a)') ' eta = 0.01, n_eta = 3, eta_grid = 0.04, 0.02, 0.01'
      else
         write(unit, '(a)') ' eta = 0.01, n_eta = 1, eta_grid = 0.01'
      end if
      if (selected_eta /= 0.01) write(unit, '(a,es16.8)') ' eta = ', selected_eta
      write(unit, '(a)') " response_lmax = -1, write_full_matrix = .false., output_file = 'lr_config.dat'"
      write(unit, '(a)') ' gf_integration_points = 3, gf_energy_margin = 1.0'
      write(unit, '(a)') '/'
      close(unit)
   end subroutine write_case

end program test_linear_response_config
