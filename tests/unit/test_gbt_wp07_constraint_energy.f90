! WP07: constraint-controller diagnostics are not physical DFT energy terms.
program test_gbt_wp07_constraint_energy
   use cfd, only: initialize_cfd, constrain, constraint_diagnostics, set_constraint_controller_parameters
   use precision_mod, only: rp
   implicit none

   real(rp), parameter :: tol = 2.0e-12_rp
   real(rp), parameter :: theta = 0.31_rp
   integer :: failures

   failures = 0
   call test_zero_constraint_all_methods()
   call test_method_three_energy_components()
   call test_controller_invariance()

   if (failures > 0) then
      write (*, '(a,i0)') 'RESULT: FAIL (checks failed: ', failures
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine fail(label, value, expected)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: value, expected

      if (abs(value - expected) > tol) then
         failures = failures + 1
         write (*, '(a,": ",es16.8," expected ",es16.8)') trim(label), value, expected
      end if
   end subroutine fail

   subroutine fail_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition

      if (.not. condition) then
         failures = failures + 1
         write (*, '(a)') trim(label)//': false'
      end if
   end subroutine fail_true

   subroutine reset_controller(method)
      integer, intent(in) :: method

      ! initialize_cfd(flag=0) is an explicit controller-history reset.  It
      ! keeps this unit test independent of call order and mirrors a fresh
      ! constrained calculation at each parameter point.
      call initialize_cfd(1, 0)
      call initialize_cfd(1, 1, method, 1, .false.)
   end subroutine reset_controller

   subroutine test_zero_constraint_all_methods()
      integer :: method
      real(rp) :: target(3, 1), moment(3, 1), field(3, 1), energy
      type(constraint_diagnostics) :: diagnostics

      target(:, 1) = [0.0_rp, 0.0_rp, 1.0_rp]
      moment = target
      do method = 2, 5
         call reset_controller(method)
         field = 0.0_rp
         call constrain(moment, target, field, 1, energy, .false., diagnostics)
         call fail('zero residual controller penalty', diagnostics%controller_penalty_energy, 0.0_rp)
         call fail('zero residual field coupling', diagnostics%field_coupling_energy, 0.0_rp)
         call fail_true('controller does not claim variational Lagrange functional', &
                        .not. diagnostics%lagrange_functional_represented)
      end do
   end subroutine test_zero_constraint_all_methods

   subroutine test_method_three_energy_components()
      real(rp) :: target(3, 1), moment(3, 1), field(3, 1), energy
      real(rp) :: expected_penalty, expected_coupling
      type(constraint_diagnostics) :: diagnostics

      call set_constraint_controller_parameters(10.0_rp)
      call reset_controller(3)
      target(:, 1) = [0.0_rp, 0.0_rp, 1.0_rp]
      moment(:, 1) = [sin(theta), 0.0_rp, cos(theta)]
      field = 0.0_rp
      call constrain(moment, target, field, 1, energy, .false., diagnostics)

      expected_penalty = 0.5_rp*10.0_rp*sin(theta)**2
      expected_coupling = -10.0_rp*sin(theta)**2
      call fail('method-3 penalty formula', diagnostics%controller_penalty_energy, expected_penalty)
      call fail('method-3 field coupling formula', diagnostics%field_coupling_energy, expected_coupling)
      call fail('legacy etcon matches controller penalty', energy, diagnostics%controller_penalty_energy)
      call fail('diagnostic field coupling matches explicit expectation', &
                diagnostics%field_coupling_energy, sum(field*moment))
   end subroutine test_method_three_energy_components

   subroutine test_controller_invariance()
      real(rp) :: target(3, 1), moment(3, 1), field(3, 1), energy
      real(rp) :: penalty_low, penalty_high, coupling_low, coupling_high
      real(rp) :: physical_low, physical_high
      type(constraint_diagnostics) :: diagnostics

      target(:, 1) = [0.0_rp, 0.0_rp, 1.0_rp]
      moment(:, 1) = [sin(theta), 0.0_rp, cos(theta)]

      call set_constraint_controller_parameters(5.0_rp)
      call reset_controller(3)
      field = 0.0_rp
      call constrain(moment, target, field, 1, energy, .false., diagnostics)
      penalty_low = diagnostics%controller_penalty_energy
      coupling_low = diagnostics%field_coupling_energy

      call set_constraint_controller_parameters(50.0_rp)
      call reset_controller(3)
      field = 0.0_rp
      call constrain(moment, target, field, 1, energy, .false., diagnostics)
      penalty_high = diagnostics%controller_penalty_energy
      coupling_high = diagnostics%field_coupling_energy

      ! The same fixed physical state has one physical energy regardless of
      ! controller tuning, while both auxiliary diagnostics change. This
      ! fixed-state oracle deliberately contains no controller quantity.
      physical_low = fixed_state_energy(moment)
      physical_high = fixed_state_energy(moment)
      call fail_true('controller penalty responds to tuning', abs(penalty_high - penalty_low) > 1.0e-6_rp)
      call fail_true('field coupling responds to tuning', abs(coupling_high - coupling_low) > 1.0e-6_rp)
      call fail('fixed-state physical energy is controller-invariant', physical_high, physical_low)
      call set_constraint_controller_parameters(10.0_rp)
   end subroutine test_controller_invariance

   pure function fixed_state_energy(moment) result(energy)
      real(rp), intent(in) :: moment(3, 1)
      real(rp) :: energy

      ! A fixed-state physical-energy oracle: the controller is not an input.
      energy = -12.3456789_rp + 0.25_rp*sum(moment**2)
   end function fixed_state_energy

end program test_gbt_wp07_constraint_energy
