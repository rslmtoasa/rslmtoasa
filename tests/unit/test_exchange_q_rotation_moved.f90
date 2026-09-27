!------------------------------------------------------------------------------
! Negative input-boundary check for the retired exchange_q rotation flag.
!------------------------------------------------------------------------------
program test_exchange_q_rotation_moved
   use exchange_q_mod, only: exchange_q_config, load_exchange_q_config
   implicit none

   type(exchange_q_config) :: config
   integer :: unit

   open(newunit=unit, file='exchange_q_rotation_moved.nml', status='replace', action='write')
   write(unit, '(a)') '&exchange_q'
   write(unit, '(a)') 'rotation_dynamics = .true.'
   write(unit, '(a)') 'n_q_points = 2'
   write(unit, '(a)') 'q_list(:,1) = 0.0, 0.0, 0.0'
   write(unit, '(a)') 'q_list(:,2) = 0.0, 0.0, 0.25'
   write(unit, '(a)') '/'
   close(unit)

   call load_exchange_q_config('exchange_q_rotation_moved.nml', config, .true.)
   error stop 'UnitExchangeQRotationMoved: expected the legacy rotation flag to fail'
end program test_exchange_q_rotation_moved
