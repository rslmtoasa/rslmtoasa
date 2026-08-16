program test_abspinlib_integrators
   !
   ! VAL-11: direct validation of the production abspinlib Depondt solver.
   ! The analytic constant-field solution is the oracle; no second numerical
   ! integrator is used here.
   !
   use Parameters, only: dblprec
   use Constants, only: gama
   use Depondt, only: allocate_depondtfields, depondt_evolve_first, depondt_evolve_second
   use RandomNumbers, only: allocate_randomwork, setup_rng_hb
   implicit none

   integer :: failures

   failures = 0
   call setup_solver()
   call test_constant_field_precession(failures)
   call test_constant_field_damping(failures)
   call test_two_spin_exchange(failures)
   call teardown_solver()

   if (failures == 0) then
      write (*, '(a)') 'VAL-11 PASS: abspinlib Depondt integrator validated against analytic LLG and exchange invariants'
   else
      write (*, '(a,i0,a)') 'VAL-11 FAIL: ', failures, ' check(s) failed'
      error stop 1
   end if

contains

   subroutine setup_solver()
      call setup_rng_hb(11, 'Y', 'N')
      call allocate_randomwork(2, 1, 1, 'N')
      call allocate_depondtfields(2, 1, 1)
   end subroutine setup_solver

   subroutine teardown_solver()
      call allocate_depondtfields(2, 1, -1)
      call allocate_randomwork(2, 1, -1, 'N')
   end subroutine teardown_solver

   subroutine test_constant_field_precession(nfail)
      integer, intent(inout) :: nfail
      integer, parameter :: nstep = 200
      real(dblprec), parameter :: dt = 1.0d-11
      real(dblprec), parameter :: field_magnitude = 1.0d-3
      real(dblprec) :: field(3), moment0(3), moment(3), expected(3)
      real(dblprec) :: norm_error, invariant_error, final_error, projection
      real(dblprec) :: max_norm_error, max_invariant_error

      field = [0.0d0, 0.0d0, field_magnitude]
      moment0 = [0.8d0, 0.0d0, 0.6d0]
      call evolve_constant_field(0.0d0, dt, nstep, field, moment0, moment)
      call exact_constant_field(0.0d0, dt*real(nstep, dblprec), field, moment0, expected)

      call exact_constant_field(0.0d0, dt*real(nstep, dblprec), field, moment0, expected)
      norm_error = abs(norm2(moment) - norm2(moment0))
      projection = dot_product(moment, field)
      invariant_error = abs(projection - dot_product(moment0, field))
      max_norm_error = norm_error
      max_invariant_error = invariant_error

      final_error = norm2(moment - expected)
      write (*, '(a,1x,3(es12.4,1x))') 'VAL-11 precession final moment', moment
      write (*, '(a,1x,es12.4,1x,a,1x,es12.4,1x,a,1x,es12.4)') &
         'VAL-11 constant-field precession: norm_err', max_norm_error, &
         'invariant_err', max_invariant_error, 'final_err', final_error

      if (max_norm_error > 2.0d-12) then
         call fail(nfail, 'constant-field precession does not preserve |m|')
      end if
      if (max_invariant_error > 2.0d-18) then
         call fail(nfail, 'zero-damping constant-field precession does not preserve m dot B')
      end if
      if (final_error > 2.0d-11) then
         call fail(nfail, 'constant-field precession disagrees with analytic frequency/direction')
      end if
      call evolve_constant_field(0.0d0, dt, 1, field, moment0, expected)
      if (expected(2) <= 0.0d0) then
         call fail(nfail, 'constant-field precession has the wrong direction')
      end if
   end subroutine test_constant_field_precession

   subroutine test_constant_field_damping(nfail)
      integer, intent(inout) :: nfail
      integer, parameter :: nstep = 2000
      real(dblprec), parameter :: dt = 1.0d-11
      real(dblprec), parameter :: alpha = 0.2d0
      real(dblprec), parameter :: field_magnitude = 1.0d-3
      real(dblprec) :: field(3), moment0(3), moment(3), updated(3), expected(3)
      real(dblprec) :: previous_projection, projection, min_increment
      real(dblprec) :: max_norm_error, final_error
      integer :: istep
      logical :: monotonic

      field = [0.0d0, 0.0d0, field_magnitude]
      moment0 = [0.8d0, 0.0d0, 0.6d0]
      moment = moment0
      previous_projection = dot_product(moment, field)
      min_increment = huge(1.0d0)
      max_norm_error = 0.0d0
      monotonic = .true.

      do istep = 1, nstep
         call evolve_constant_field(alpha, dt, 1, field, moment, updated)
         moment = updated
         projection = dot_product(moment, field)
         min_increment = min(min_increment, projection - previous_projection)
         if (projection < previous_projection - 2.0d-18) monotonic = .false.
         max_norm_error = max(max_norm_error, abs(norm2(moment) - norm2(moment0)))
         previous_projection = projection
      end do

      call exact_constant_field(alpha, dt*real(nstep, dblprec), field, moment0, expected)
      final_error = norm2(moment - expected)
      write (*, '(a,1x,3(es12.4,1x),a,1x,es12.4)') &
         'VAL-11 damping final moment', moment, 'norm', norm2(moment)
      write (*, '(a,1x,l1,1x,a,1x,es12.4,1x,a,1x,es12.4)') &
         'VAL-11 constant-field damping: monotonic', monotonic, &
         'min_increment', min_increment, 'final_err', final_error

      if (.not. monotonic) then
         call fail(nfail, 'damped constant-field motion is not monotonic toward the field')
      end if
      if (max_norm_error > 2.0d-12) then
         call fail(nfail, 'damped constant-field motion does not preserve |m|')
      end if
      if (final_error > 2.0d-6) then
         call fail(nfail, 'damped constant-field motion disagrees with analytic LLG')
      end if
      if (dot_product(moment, field) <= dot_product(moment0, field)) then
         call fail(nfail, 'damped constant-field motion does not approach the field')
      end if

      call check_timestep_convergence(nfail, field, moment0)
   end subroutine test_constant_field_damping

   subroutine check_timestep_convergence(nfail, field, moment0)
      integer, intent(inout) :: nfail
      real(dblprec), intent(in) :: field(3), moment0(3)
      real(dblprec), parameter :: alpha = 0.2d0
      real(dblprec), parameter :: final_time = 2.0d-8
      real(dblprec) :: coarse(3), medium(3), fine(3), expected(3)
      real(dblprec) :: err_coarse, err_medium, err_fine, order

      call evolve_constant_field(alpha, 2.0d-11, 1000, field, moment0, coarse)
      call evolve_constant_field(alpha, 1.0d-11, 2000, field, moment0, medium)
      call evolve_constant_field(alpha, 5.0d-12, 4000, field, moment0, fine)
      call exact_constant_field(alpha, final_time, field, moment0, expected)
      err_coarse = norm2(coarse - expected)
      err_medium = norm2(medium - expected)
      err_fine = norm2(fine - expected)
      order = log(err_medium/err_fine)/log(2.0d0)

      write (*, '(a,1x,3(es12.4,1x),a,1x,f8.4)') &
         'VAL-11 timestep convergence errors', err_coarse, err_medium, err_fine, 'order', order
      if (.not. (err_fine < err_medium .and. err_medium < err_coarse)) then
         call fail(nfail, 'timestep refinement does not reduce constant-field damping error')
      end if
      if (order < 1.5d0) then
         call fail(nfail, 'constant-field damping convergence is below second-order behavior')
      end if
   end subroutine check_timestep_convergence

   subroutine test_two_spin_exchange(nfail)
      integer, intent(inout) :: nfail
      integer, parameter :: nstep = 2000
      real(dblprec), parameter :: dt = 5.0d-12
      real(dblprec), parameter :: exchange = 1.0d-3
      real(dblprec) :: moments(3, 2), initial(3, 2), fields(3, 2)
      real(dblprec) :: energy0, energy, total0(3), total(3)
      real(dblprec) :: norm_error, energy_error, total_error
      integer :: istep

      initial(:, 1) = [1.0d0, 0.0d0, 0.0d0]
      initial(:, 2) = [0.0d0, 0.0d0, 1.0d0]
      moments = initial
      energy0 = -exchange*dot_product(initial(:, 1), initial(:, 2))
      total0 = initial(:, 1) + initial(:, 2)

      do istep = 1, nstep
         fields(:, 1) = exchange*moments(:, 2)
         fields(:, 2) = exchange*moments(:, 1)
         call evolve_exchange(0.0d0, dt, exchange, fields, moments)
      end do

      energy = -exchange*dot_product(moments(:, 1), moments(:, 2))
      total = moments(:, 1) + moments(:, 2)
      norm_error = max(abs(norm2(moments(:, 1)) - 1.0d0), abs(norm2(moments(:, 2)) - 1.0d0))
      energy_error = abs(energy - energy0)
      total_error = norm2(total - total0)
      write (*, '(a,1x,3(es12.4,1x))') &
         'VAL-11 two-spin exchange errors norm energy total', norm_error, energy_error, total_error

      if (norm_error > 2.0d-12) then
         call fail(nfail, 'two-spin exchange does not preserve spin norms')
      end if
      if (energy_error > 2.0d-12) then
         call fail(nfail, 'zero-damping two-spin exchange does not conserve energy')
      end if
      if (total_error > 2.0d-10) then
         call fail(nfail, 'two-spin exchange does not conserve total spin sufficiently')
      end if
   end subroutine test_two_spin_exchange

   subroutine evolve_constant_field(alpha, dt, nstep, field, initial, final)
      real(dblprec), intent(in) :: alpha, dt, field(3), initial(3)
      integer, intent(in) :: nstep
      real(dblprec), intent(out) :: final(3)
      real(dblprec) :: beff(3, 1, 1), b2eff(3, 1, 1), btorque(3, 1, 1)
      real(dblprec) :: emom(3, 1, 1), emom2(3, 1, 1), emomM(3, 1, 1), thermal(3, 1, 1)
      real(dblprec) :: mmom(1, 1), lambda(1), temp(1), temprescale
      integer :: istep

      ! The Depondt API takes the effective field directly.  The higher-level
      ! spin_dynamics wrapper applies its own field-sign convention before
      ! calling this API; VAL-11 intentionally validates this library seam.
      beff(:, 1, 1) = field
      btorque = 0.0d0
      thermal = 0.0d0
      lambda(1) = alpha
      temp(1) = 0.0d0
      temprescale = 1.0d0
      mmom(1, 1) = norm2(initial)
      emom(:, 1, 1) = initial/mmom(1, 1)
      emom2 = 0.0d0
      emomM = 0.0d0
      b2eff = 0.0d0

      do istep = 1, nstep
         call depondt_evolve_first(1, 1, lambda, beff, b2eff, btorque, emom, emom2, emomM, mmom, &
                                    dt, temp, temprescale, 'N', thermal, 'N', btorque)
         call depondt_evolve_second(1, 1, lambda, beff, b2eff, btorque, emom, emom2, dt, 'N', 'N', btorque)
         emom(:, 1, 1) = emom2(:, 1, 1)
      end do
      final = emom(:, 1, 1)*mmom(1, 1)
   end subroutine evolve_constant_field

   subroutine evolve_exchange(alpha, dt, exchange, fields, moments)
      real(dblprec), intent(in) :: alpha, dt, exchange, fields(3, 2)
      real(dblprec), intent(inout) :: moments(3, 2)
      real(dblprec) :: beff(3, 2, 1), b2eff(3, 2, 1), btorque(3, 2, 1)
      real(dblprec) :: emom(3, 2, 1), emom2(3, 2, 1), emomM(3, 2, 1), thermal(3, 2, 1)
      real(dblprec) :: mmom(2, 1), lambda(2), temp(2), temprescale
      integer :: i

      beff(:, :, 1) = fields
      btorque = 0.0d0
      thermal = 0.0d0
      lambda = alpha
      temp = 0.0d0
      temprescale = 1.0d0
      do i = 1, 2
         mmom(i, 1) = norm2(moments(:, i))
         emom(:, i, 1) = moments(:, i)/mmom(i, 1)
      end do
      emom2 = 0.0d0
      emomM = 0.0d0
      b2eff = 0.0d0

      call depondt_evolve_first(2, 1, lambda, beff, b2eff, btorque, emom, emom2, emomM, mmom, &
                                dt, temp, temprescale, 'N', thermal, 'N', btorque)
      ! Exchange fields are prescribed from the predictor configuration for
      ! the corrector, matching the production predictor/corrector semantics.
      beff(:, 1, 1) = exchange*emom(:, 2, 1)
      beff(:, 2, 1) = exchange*emom(:, 1, 1)
      call depondt_evolve_second(2, 1, lambda, beff, b2eff, btorque, emom, emom2, dt, 'N', 'N', btorque)
      do i = 1, 2
         moments(:, i) = emom2(:, i, 1)*mmom(i, 1)
      end do
   end subroutine evolve_exchange

   subroutine exact_constant_field(alpha, time, field, initial, exact)
      real(dblprec), intent(in) :: alpha, time, field(3), initial(3)
      real(dblprec), intent(out) :: exact(3)
      real(dblprec) :: bmag, u0, u, rate, angle
      real(dblprec) :: bhat(3), perpendicular(3), rotated(3)

      bmag = norm2(field)
      bhat = field/bmag
      u0 = dot_product(initial, bhat)
      rate = gama*alpha*bmag/(1.0d0 + alpha*alpha)
      u = tanh(atanh(u0) + rate*time)
      perpendicular = initial - u0*bhat
      angle = gama*bmag*time/(1.0d0 + alpha*alpha)
      rotated = cos(angle)*perpendicular + sin(angle)*cross_product(bhat, perpendicular)
      exact = u*bhat + rotated*sqrt((1.0d0 - u*u)/(1.0d0 - u0*u0))
   end subroutine exact_constant_field

   pure function norm2(vector) result(norm)
      real(dblprec), intent(in) :: vector(:)
      real(dblprec) :: norm
      norm = sqrt(dot_product(vector, vector))
   end function norm2

   pure function cross_product(a, b) result(c)
      real(dblprec), intent(in) :: a(3), b(3)
      real(dblprec) :: c(3)
      c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
   end function cross_product

   subroutine fail(nfail, message)
      integer, intent(inout) :: nfail
      character(len=*), intent(in) :: message
      nfail = nfail + 1
      write (*, '(a)') 'FAIL: '//trim(message)
   end subroutine fail

end program test_abspinlib_integrators
