! WP06: nonzero constraining-field fixture and explicit commensurate GBT oracle.
program test_gbt_constraint_covariance
   use cfd, only: initialize_cfd, constrain, constraint_diagnostics
   use precision_mod, only: rp
   use gbt_structure_mod, only: gbt_frame_t, gbt_frame_for_cell, gbt_rotating_to_lab_vector
   implicit none

   real(rp), parameter :: tol = 2.0e-11_rp
   real(rp), parameter :: pi = acos(-1.0_rp)
   integer :: failures

   failures = 0
   call test_nonzero_field_fixture()
   call test_transverse_controller_variants()
   call test_commensurate_covariance()

   if (failures > 0) then
      write (*, '(a,i0)') 'RESULT: FAIL (checks failed: ', failures
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine check(value, label)
      real(rp), intent(in) :: value
      character(len=*), intent(in) :: label

      if (value > tol) then
         failures = failures + 1
         write (*, '(a,a,a,es14.6)') 'FAIL ', trim(label), ': ', value
      end if
   end subroutine check

   subroutine check_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         failures = failures + 1
         write (*, '(a,a)') 'FAIL ', trim(label)
      end if
   end subroutine check_true

   subroutine test_nonzero_field_fixture()
      type(constraint_diagnostics) :: diagnostics
      real(rp) :: target(3, 1), moment(3, 1), field(3, 1), energy
      real(rp) :: angle, first_metric, last_metric
      integer :: iteration

      call initialize_cfd(1, 1, 3, 1, .false.)
      target(:, 1) = [0.0_rp, 0.0_rp, 1.0_rp]

      ! A deterministic frozen-response fixture: the measured moment starts
      ! canted and approaches the fixed target across controller updates.
      first_metric = 0.0_rp
      last_metric = 0.0_rp
      do iteration = 1, 11
         angle = 0.36_rp*(0.25_rp**real(iteration - 1, rp))
         moment(:, 1) = [sin(angle), 0.0_rp, cos(angle)]
         field = 0.0_rp
         call constrain(moment, target, field, 1, energy, .false., diagnostics)
         if (iteration == 1) first_metric = diagnostics%convergence_metric
         last_metric = diagnostics%convergence_metric
         if (iteration == 1) call check_true(diagnostics%bfield_magnitude(1) > 1.0e-8_rp, &
                                              'canted fixture produces a nonzero field')
         call check(abs(diagnostics%moment_magnitude(1) - 1.0_rp), 'fixture moment magnitude')
         call check(abs(diagnostics%angular_error(1) - angle), 'fixture angular error')
         call check(abs(diagnostics%bfield_longitudinal(1)), 'fixture field is transverse')
         call check(abs(sqrt(sum(diagnostics%transverse_residual(:, 1)**2)) - sin(angle)), &
                     'fixture transverse residual')
         call check(abs(sqrt(sum(diagnostics%torque(:, 1)**2)) - sin(angle)), 'fixture torque')
      end do
      call check_true(first_metric > 0.1_rp, 'zero-field fixture has measurable orientation error')
      call check_true(last_metric < first_metric, 'constraint fixture reduces angular residual')
      call check_true(last_metric < 1.0e-6_rp, 'constraint fixture reaches the declared tolerance')
      call check_true(diagnostics%bfield_magnitude(1) > 1.0e-8_rp, 'converged fixture field remains nonzero')
   end subroutine test_nonzero_field_fixture

   subroutine test_transverse_controller_variants()
      type(constraint_diagnostics) :: diagnostics
      real(rp) :: target(3, 1), moment(3, 1), field(3, 1), energy

      target(:, 1) = [0.0_rp, 0.0_rp, 1.0_rp]
      moment(:, 1) = [sin(0.31_rp), 0.0_rp, cos(0.31_rp)]

      call initialize_cfd(1, 1, 4, 1, .false.)
      field(:, 1) = [0.0_rp, 0.0_rp, 0.2_rp]
      call constrain(moment, target, field, 1, energy, .false., diagnostics)
      call check(abs(diagnostics%bfield_longitudinal(1)), 'PID field is transverse')

      call initialize_cfd(1, 1, 5, 1, .false.)
      field(:, 1) = [0.0_rp, 0.0_rp, 0.2_rp]
      call constrain(moment, target, field, 1, energy, .false., diagnostics)
      call check(abs(diagnostics%bfield_longitudinal(1)), 'Ma-Dudarev field is transverse')
   end subroutine test_transverse_controller_variants

   subroutine test_commensurate_covariance()
      type(gbt_frame_t) :: frame
      type(constraint_diagnostics) :: diagnostics
      real(rp) :: q(3), r(3), local_target(3), local_moment(3), local_field(3)
      real(rp) :: target(3, 1), moment(3, 1), field(3, 1), expected(3), energy
      real(rp) :: phase, theta
      integer :: n

      call initialize_cfd(1, 1, 3, 1, .false.)
      q = [0.25_rp, 0.0_rp, 0.0_rp]
      local_target = [0.0_rp, 0.0_rp, 1.0_rp]
      ! Use the final, still-nonzero residual scale of the fixture above so
      ! the covariance comparison is a converged primitive-field snapshot.
      theta = 0.36_rp*(0.25_rp**10)
      local_moment = [sin(theta), 0.0_rp, cos(theta)]
      local_field = 0.0_rp
      target(:, 1) = local_target
      moment(:, 1) = local_moment
      field = 0.0_rp
      call constrain(moment, target, field, 1, energy, .false., diagnostics)
      local_field = field(:, 1)

      ! Four translations are one explicit period-four supercell.  The field
      ! is an onsite spin vector, so it must transform with the same R[U] as
      ! the target and measured moment, while its magnitude/angle remain fixed.
      do n = 0, 3
         r = [real(n, rp), 0.0_rp, 0.0_rp]
         phase = 2.0_rp*pi*q(1)*real(n, rp)
         call gbt_frame_for_cell(q, r, 1.0_rp, 0.27_rp, -0.18_rp, frame)
         call gbt_rotating_to_lab_vector(frame, local_target, target(:, 1))
         call gbt_rotating_to_lab_vector(frame, local_moment, moment(:, 1))
         field = 0.0_rp
         call constrain(moment, target, field, 1, energy, .false., diagnostics)
         expected = 0.0_rp
         call gbt_rotating_to_lab_vector(frame, local_field, expected)
         call check(maxval(abs(field(:, 1) - expected)), 'period-four field covariance')
         call check(abs(diagnostics%bfield_magnitude(1) - sqrt(sum(local_field**2))), &
                     'period-four field magnitude covariance')
         call check(abs(diagnostics%angular_error(1) - theta), 'period-four angle covariance')
         call check(abs(diagnostics%bfield_longitudinal(1)), 'period-four transversality')
         call check(abs(phase - 0.5_rp*pi*real(n, rp)), 'period-four phase')
      end do
   end subroutine test_commensurate_covariance

end program test_gbt_constraint_covariance
