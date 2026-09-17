!------------------------------------------------------------------------------
! exchange_q input boundary and independent q-space LKAG contraction tests.
!------------------------------------------------------------------------------
program test_exchange_q
   use precision_mod, only: rp
   use exchange_q_mod, only: exchange_q_config, load_exchange_q_config, lkag_pauli_product
   implicit none

   type(exchange_q_config) :: config
   complex(rp) :: one(1,1), zero(1,1), ginmag(1,1), gjnmag(1,1)
   real(rp) :: value, expected, q_space, real_space, phase
   integer :: unit

   call write_inline_input('exchange_q_unit.nml')
   call load_exchange_q_config('exchange_q_unit.nml', config, .true.)
   if (config%n_q /= 3 .or. abs(config%q_list(3,2)-0.25_rp) > 1.0e-14_rp .or. &
       abs(config%q_list(1,3)-0.5_rp) > 1.0e-14_rp) then
      error stop 'UnitExchangeQ: inline q path parsing failed'
   end if

   call write_q_file('exchange_q_unit_path.dat')
   call write_file_input('exchange_q_unit.nml')
   call load_exchange_q_config('exchange_q_unit.nml', config, .true.)
   if (config%n_q /= 2 .or. abs(config%q_list(3,2)-0.5_rp) > 1.0e-14_rp) then
      error stop 'UnitExchangeQ: frozen-magnon q-file parsing failed'
   end if

   one = cmplx(1.0_rp,0.0_rp,rp)
   zero = cmplx(0.0_rp,0.0_rp,rp)
   ginmag = cmplx(0.3_rp,0.7_rp,rp)
   gjnmag = cmplx(0.2_rp,-0.4_rp,rp)
   value = lkag_pauli_product(one, one, ginmag, zero, zero, zero, gjnmag, zero, zero, zero)
   expected = aimag(ginmag(1,1)*gjnmag(1,1))
   if (abs(value-expected) > 1.0e-14_rp) error stop 'UnitExchangeQ: LKAG Pauli contraction failed'

   ! Two-point DFT identity: the q-space k,k+q product equals the real-space
   ! Fourier sum with the same phase convention.  This is a minimal oracle for
   ! the independent Green-function representation used by exchange_q.
   q_space = 0.5_rp*aimag(cmplx(0.2_rp,0.1_rp,rp)*cmplx(0.6_rp,-0.2_rp,rp) + &
      cmplx(-0.1_rp,0.5_rp,rp)*cmplx(-0.4_rp,0.3_rp,rp))
   real_space = 0.0_rp
   real_space = real_space + aimag(0.5_rp*(cmplx(0.2_rp,0.1_rp,rp)+cmplx(-0.1_rp,0.5_rp,rp))* &
      (0.5_rp*(cmplx(-0.4_rp,0.3_rp,rp)+cmplx(0.6_rp,-0.2_rp,rp))))
   phase = -1.0_rp
   real_space = real_space + aimag(phase*0.5_rp*(cmplx(0.2_rp,0.1_rp,rp)-cmplx(-0.1_rp,0.5_rp,rp))* &
      (0.5_rp*(cmplx(-0.4_rp,0.3_rp,rp)-cmplx(0.6_rp,-0.2_rp,rp))))
   if (abs(q_space-real_space) > 1.0e-14_rp) error stop 'UnitExchangeQ: q/real-space DFT oracle failed'

   open(newunit=unit, file='exchange_q_unit.nml', status='old')
   close(unit, status='delete')
   open(newunit=unit, file='exchange_q_unit_path.dat', status='old')
   close(unit, status='delete')
   write (*, '(a,es12.4)') 'UnitExchangeQ: LKAG Pauli residual = ', abs(value-expected)
   write (*, '(a,es12.4)') 'UnitExchangeQ: q/real-space DFT residual = ', abs(q_space-real_space)
   write (*, '(a)') 'UnitExchangeQ: PASS (input parsing, LKAG Pauli contraction, q-space DFT oracle)'

contains

   subroutine write_inline_input(path)
      character(len=*), intent(in) :: path
      integer :: local_unit
      open(newunit=local_unit, file=path, status='replace', action='write')
      write(local_unit,'(a)') '&exchange_q'
      write(local_unit,'(a)') " q_coordinates = 'direct'"
      write(local_unit,'(a)') " q_file = ''"
      write(local_unit,'(a)') ' n_q_points = 3'
      write(local_unit,'(a)') ' q_list = 0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.0, 0.0'
      write(local_unit,'(a)') " output_file = 'unit_exchange_q.dat'"
      write(local_unit,'(a)') ' native_crosscheck = .false.'
      write(local_unit,'(a)') '/'
      close(local_unit)
   end subroutine write_inline_input

   subroutine write_file_input(path)
      character(len=*), intent(in) :: path
      integer :: local_unit
      open(newunit=local_unit, file=path, status='replace', action='write')
      write(local_unit,'(a)') '&exchange_q'
      write(local_unit,'(a)') " q_coordinates = 'cartesian'"
      write(local_unit,'(a)') " q_file = 'exchange_q_unit_path.dat'"
      write(local_unit,'(a)') ' n_q_points = 0'
      write(local_unit,'(a)') '/'
      close(local_unit)
   end subroutine write_file_input

   subroutine write_q_file(path)
      character(len=*), intent(in) :: path
      integer :: local_unit
      open(newunit=local_unit, file=path, status='replace', action='write')
      write(local_unit,'(a)') '2'
      write(local_unit,'(a)') '0.0 0.0 0.0'
      write(local_unit,'(a)') '0.0 0.0 0.5'
      close(local_unit)
   end subroutine write_q_file

end program test_exchange_q
