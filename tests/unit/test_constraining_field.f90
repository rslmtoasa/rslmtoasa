program test_constraining_field
   use cfd, only: initialize_cfd, constrain
   use Parameters, only: dblprec
   implicit none

   real(dblprec), parameter :: tol = 2.0d-9
   real(dblprec) :: mom(3, 1), ref(3, 1), field(3, 1), energy
   real(dblprec) :: theta, h, eplus, eminus, derivative
   real(dblprec) :: rotation(3, 3), rotated(3), expected(3)
   integer :: failures

   failures = 0
   call initialize_cfd(1, 1, 3, 1)

   ! Aligned moments have no transverse constraint residual.
   mom(:, 1) = [0.0d0, 0.0d0, 1.0d0]
   ref(:, 1) = mom(:, 1)
   field = 0.0d0
   call constrain(mom, ref, field, 1, energy, .false.)
   call check('aligned field', maxval(abs(field)), 0.0d0, failures)
   call check('aligned energy', energy, 0.0d0, failures)

   ! A controlled canting in +x produces a Hamiltonian coefficient in -x.
   theta = 0.3d0
   mom(:, 1) = [sin(theta), 0.0d0, cos(theta)]
   ref(:, 1) = [0.0d0, 0.0d0, 1.0d0]
   field = 0.0d0
   call constrain(mom, ref, field, 1, energy, .false.)
   call check_true('canted field opposes deviation', field(1, 1) < 0.0d0, failures)
   call check('canted field longitudinal component', field(3, 1), 0.0d0, failures)
   call check('canted energy', energy, 5.0d0*sin(theta)**2, failures)

   ! A restart seed is retained only in the transverse subspace; a
   ! longitudinal component is not a constraint torque.
   field(:, 1) = [0.0d0, 0.0d0, 0.25d0]
   call constrain(mom, ref, field, 1, energy, .false.)
   call check('seed longitudinal component removed', field(3, 1), 0.0d0, failures)

   ! For E = lambda*|m_perp|^2/2, B = -dE/dm_perp.
   h = 1.0d-6
   mom(:, 1) = [sin(theta + h), 0.0d0, cos(theta + h)]
   field = 0.0d0
   call constrain(mom, ref, field, 1, eplus, .false.)
   mom(:, 1) = [sin(theta - h), 0.0d0, cos(theta - h)]
   field = 0.0d0
   call constrain(mom, ref, field, 1, eminus, .false.)
   derivative = (eplus - eminus)/(2.0d0*h)
   call check('finite-difference penalty derivative', derivative, 10.0d0*sin(theta)*cos(theta), failures)

   ! Rotate both input vectors in the SOC-free spin space; the field and energy
   ! must rotate covariantly and the scalar penalty must remain invariant.
   rotation = 0.0d0
   rotation(1, 1) = cos(0.47d0); rotation(1, 3) = sin(0.47d0)
   rotation(2, 2) = 1.0d0
   rotation(3, 1) = -sin(0.47d0); rotation(3, 3) = cos(0.47d0)
   mom(:, 1) = [sin(theta), 0.0d0, cos(theta)]
   ref(:, 1) = [0.0d0, 0.0d0, 1.0d0]
   field = 0.0d0
   call constrain(mom, ref, field, 1, energy, .false.)
   expected = field(:, 1)
   rotated = matmul(rotation, mom(:, 1))
   ref(:, 1) = matmul(rotation, ref(:, 1))
   mom(:, 1) = rotated
   field = 0.0d0
   call constrain(mom, ref, field, 1, eplus, .false.)
   call check('global spin rotation field', maxval(abs(field(:, 1) - matmul(rotation, expected))), 0.0d0, failures)
   call check('global spin rotation energy', eplus, energy, failures)

   if (failures > 0) error stop 1
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine check(label, value, target, count)
      character(len=*), intent(in) :: label
      real(dblprec), intent(in) :: value, target
      integer, intent(inout) :: count
      if (abs(value - target) > tol) then
         write (*, '(a,": ",es14.6," expected ",es14.6)') trim(label), value, target
         count = count + 1
      end if
   end subroutine check

   subroutine check_true(label, condition, count)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      integer, intent(inout) :: count
      if (.not. condition) then
         write (*, '(a)') trim(label)//': false'
         count = count + 1
      end if
   end subroutine check_true

end program test_constraining_field
