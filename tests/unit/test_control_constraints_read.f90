! &constraints: constraints_mom_ref(:,i) = ... is read into control%constraints_mom_ref with shape (3, last set site).
! Expected values are the literals written into the input files below. Before the fix the indexed assignment was out of
! range for the (1,1) pre-allocation, the read failed and the array stayed unset.
program test_control_constraints_read
   use precision_mod, only: rp
   use control_mod, only: control
   use logger_mod, only: g_logger
   implicit none

   logical :: failed = .false.

   call g_logger%init()
   call read_case('three sites, third set along x', &
      'constraints_mom_ref(:,1) = 0.0, 0.0, 1.0'//new_line('a')// &
      'constraints_mom_ref(:,2) = 0.0, 0.0, -1.0'//new_line('a')// &
      'constraints_mom_ref(:,3) = 1.0, 0.0, 0.0', &
      3, reshape([0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, -1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp], [3, 3]))
   call read_case('one site, component form', &
      'constraints_mom_ref(1,1) = 0.0'//new_line('a')// &
      'constraints_mom_ref(2,1) = 0.6'//new_line('a')// &
      'constraints_mom_ref(3,1) = 0.8', &
      1, reshape([0.0_rp, 0.6_rp, 0.8_rp], [3, 1]))
   call read_case('not set: stays the unset (1,1) form', 'constraints_tolerance = 1.0e-6', 0, reshape([0.0_rp], [1, 1]))

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine read_case(name, body, ncols, expected)
      character(len=*), intent(in) :: name, body
      integer, intent(in) :: ncols
      real(rp), intent(in) :: expected(:, :)
      type(control) :: ctl
      integer :: unit
      real(rp) :: worst
      character(len=*), parameter :: fname = 'constraints_read_case.nml'

      open (newunit=unit, file=fname, status='replace', action='write')
      write (unit, '(a)') '&control'
      write (unit, '(a)') "calctype = 'B'"
      write (unit, '(a)') 'nsp = 2'
      write (unit, '(a)') 'lld = 20'
      write (unit, '(a)') "recur = 'block'"
      write (unit, '(a)') '/'
      write (unit, '(a)') '&constraints'
      write (unit, '(a)') 'constraints_enable = .true.'
      write (unit, '(a)') 'constraints_i_cons = 3'
      write (unit, '(a)') body
      write (unit, '(a)') '/'
      close (unit, status='keep')
      call ctl%restore_to_default()
      call ctl%build_from_file(fname)
      if (.not. allocated(ctl%constraints_mom_ref)) then
         write (*, '(a)') 'FAIL '//name//': constraints_mom_ref not allocated'
         failed = .true.
         return
      end if
      if (any(shape(ctl%constraints_mom_ref) /= shape(expected))) then
         write (*, '(a,2i6)') 'FAIL '//name//': shape', shape(ctl%constraints_mom_ref)
         failed = .true.
         return
      end if
      worst = maxval(abs(ctl%constraints_mom_ref - expected))
      write (*, '(a,1x,a,es10.3,a,i0)') merge('PASS', 'FAIL', worst == 0.0_rp), name//': max |read - written| ', worst, &
         '; columns ', size(ctl%constraints_mom_ref, 2)
      if (worst /= 0.0_rp) failed = .true.
      if (ncols == 0 .and. .not. ctl%constraints_enable) failed = .true.
   end subroutine read_case

end program test_control_constraints_read
