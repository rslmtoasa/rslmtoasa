! TDCOV-01 -- shared response occupation/Fermi lifecycle.
program test_tddft_occupation_lifecycle
   use precision_mod, only: rp
   use response_components_mod, only: RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel
   use tddft_occupation_mod, only: tddft_response_occupation_state
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, install_tddft_occupation
   use tddft_chi0_green_mod, only: green_chi0_options, install_green_occupation
   use tddft_backend_mod, only: tddft_eigenpair_backend
   implicit none

   character(len=32) :: argument
   logical :: failed

   call get_command_argument(1, argument)
   if (trim(argument) == '--negative') then
      call unresolved_backend_must_fail()
      error stop 'negative unresolved-occupation test unexpectedly returned'
   end if

   failed = .false.
   call test_fixed_fermi_state()
   call test_auto_resolved_state()
   call test_circular_metadata_equality()
   if (failed) error stop 1
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_fixed_fermi_state()
      type(tddft_response_occupation_state) :: state
      type(tddft_chi0_options) :: primary, reverse
      type(green_chi0_options) :: green_primary, green_reverse
      type(tddft_eigenpair_backend) :: backend
      type(tddft_chi0_result) :: result

      call state%set(-0.123456_rp, 300.0_rp, 1.0e-14_rp, 1, 0, 'ground_state_inherited', 'fixed_ground_state')
      primary%eta = 0.001_rp
      primary%response_projection = 'site'
      primary%circular_channel = 'plus_minus'
      call install_tddft_occupation(primary, state)
      reverse = primary
      reverse%circular_channel = 'minus_plus'
      call install_green_occupation(green_primary, state)
      green_reverse = green_primary
      green_reverse%circular_channel = 'minus_plus'
      green_primary%eta = primary%eta
      green_primary%green_eta = 0.0005_rp
      green_primary%energy_min = -0.5_rp
      green_primary%energy_max = 0.5_rp
      green_primary%energy_points = 101
      green_reverse = green_primary
      green_reverse%circular_channel = 'minus_plus'

      call backend%initialize([1.0_rp], reshape([-0.2_rp, 0.2_rp], [2, 1]), identity(2, 1), &
         reshape([-0.2_rp, 0.2_rp], [2, 1, 1]), reshape(identity(2, 1), [2, 2, 1, 1]), &
         reshape([0.0_rp, 0.0_rp, 0.0_rp], [3, 1]), &
         [1], [response_channel(1, RESPONSE_PLUS)], [response_channel(1, RESPONSE_MINUS)], primary, green_primary)
      call backend%evaluate_one([0.0_rp, 0.0_rp, 0.0_rp], 0.05_rp, result)
      call check('fixed EF is preserved exactly in chi0 metadata', result%metadata%fermi_level == -0.123456_rp)
      call check('fixed EF source is inherited', trim(result%metadata%fermi_source) == 'ground_state_inherited')
      call check('fixed EF policy is reported', trim(result%metadata%fermi_policy) == 'fixed_ground_state')
      call check('primary/reverse chi occupation state matches', primary%occupation_state%same_as(reverse%occupation_state))
      call check('primary/reverse Green occupation state matches', green_primary%occupation_state%same_as(green_reverse%occupation_state))
      call check('fixed EF is not the unset sentinel', abs(primary%fermi_level) < 1.0_rp)
   end subroutine test_fixed_fermi_state

   subroutine test_auto_resolved_state()
      type(tddft_response_occupation_state) :: state
      type(tddft_chi0_options) :: primary, reverse
      type(green_chi0_options) :: green_primary, green_reverse

      ! The deterministic fixture stands for the one response-mesh result.  A
      ! single state object is copied only after that result is available.
      call state%set(0.0375_rp, 600.0_rp, 2.0e-12_rp, 2, 11, 'response_mesh_recomputed', 'auto_find_fermi')
      call install_tddft_occupation(primary, state)
      reverse = primary
      reverse%circular_channel = 'minus_plus'
      call install_green_occupation(green_primary, state)
      green_reverse = green_primary
      green_reverse%circular_channel = 'minus_plus'
      call check('auto-resolved state is marked resolved', state%fermi_is_resolved)
      call check('auto-resolved EF reaches primary/reverse chi', primary%fermi_level == reverse%fermi_level)
      call check('auto-resolved EF reaches primary/reverse Green', green_primary%fermi_level == green_reverse%fermi_level)
      call check('auto-resolved source is response mesh', trim(state%fermi_source) == 'response_mesh_recomputed')
   end subroutine test_auto_resolved_state

   subroutine test_circular_metadata_equality()
      type(tddft_response_occupation_state) :: state
      type(tddft_chi0_options) :: primary, reverse
      type(green_chi0_options) :: green_primary, green_reverse

      call state%set(-0.123456_rp, 300.0_rp, 1.0e-14_rp, 1, 0, 'ground_state_inherited', 'fixed_ground_state')
      call install_tddft_occupation(primary, state)
      reverse = primary
      reverse%circular_channel = 'minus_plus'
      call install_green_occupation(green_primary, state)
      green_reverse = green_primary
      green_reverse%circular_channel = 'minus_plus'
      call check('circular chi thermodynamic metadata is identical', primary%occupation_state%same_as(reverse%occupation_state))
      call check('circular Green thermodynamic metadata is identical', green_primary%occupation_state%same_as(green_reverse%occupation_state))
      call check('circular channel remains the only changed chi label', trim(primary%circular_channel) == 'plus_minus' .and. &
         trim(reverse%circular_channel) == 'minus_plus')
   end subroutine test_circular_metadata_equality

   subroutine unresolved_backend_must_fail()
      type(tddft_eigenpair_backend) :: backend
      type(tddft_chi0_options) :: options
      type(response_channel) :: plus(1), minus(1)

      plus(1) = response_channel(1, RESPONSE_PLUS)
      minus(1) = response_channel(1, RESPONSE_MINUS)
      call backend%initialize([1.0_rp], reshape([-0.2_rp, 0.2_rp], [2, 1]), identity(2, 1), &
         reshape([-0.2_rp, 0.2_rp], [2, 1, 1]), reshape(identity(2, 1), [2, 2, 1, 1]), &
         reshape([0.0_rp, 0.0_rp, 0.0_rp], [3, 1]), &
         [1], plus, minus, options)
   end subroutine unresolved_backend_must_fail

   function identity(n, nk) result(matrix)
      integer, intent(in) :: n, nk
      complex(rp) :: matrix(n, n, nk)
      integer :: i, ik

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do ik = 1, nk
         do i = 1, n
            matrix(i, i, ik) = cmplx(1.0_rp, 0.0_rp, rp)
         end do
      end do
   end function identity

   subroutine check(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition

      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine check

end program test_tddft_occupation_lifecycle
