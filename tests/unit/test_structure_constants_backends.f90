program test_structure_constants_backends
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use control_mod, only: control
   use lattice_mod, only: lattice
   use timer_mod, only: g_timer, timer
   implicit none

! The two backends are called through lattice%structb; this test deliberately
! contains no structure-constant formula.  The manual screening values below
! are the legacy micha values before its internal fak=2 conversion.  Sdot is
! intentionally not requested: equivalence of the two Sdot conventions is a
! separate follow-up test.
!
! Observed on the first comparison run (Debug, serial build, bcc Fe fixture):
! sp  max(abs(delta Sbar))=1.4024e-8, relative=7.0387e-9;
! spd max(abs(delta Sbar))=1.2077e-12, relative=4.6697e-13.
! The tolerance is therefore 2e-8 absolute and 1e-8 relative: a narrow margin
! over the measured sp discrepancy, not a loose backend-normalizing bound.

   character(len=*), parameter :: database_directory = TEST_STRUCTURE_CONSTANT_DATABASE
   real(rp), parameter :: alpha_legacy(4) = [0.3485_rp, 0.05303_rp, 0.010714_rp, 0.00337_rp]
   real(rp), parameter :: abs_tol = 2.0e-8_rp
   real(rp), parameter :: rel_tol = 1.0e-8_rp

   type :: backend_state
      integer :: lmax = 0
      integer :: norb = 0
      integer :: kk = 0
      integer :: ntot = 0
      integer :: nn_max = 0
      integer :: sbar_dim = 0
      integer, allocatable :: nn(:, :)
      real(rp), allocatable :: sbarvec(:, :)
      complex(rp), allocatable :: sbar(:, :, :, :)
   end type backend_state

   logical :: failed

   g_timer = timer()
   failed = .false.
   call run_basis_case(1, 'sp')
   call run_basis_case(2, 'spd')

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine run_basis_case(lmax, basis_name)
      integer, intent(in) :: lmax
      character(len=*), intent(in) :: basis_name
      type(backend_state) :: legacy, strux

      call run_backend(lmax, basis_name, 'legacy', legacy)
      call run_backend(lmax, basis_name, 'strux_lib', strux)
      call compare_states(lmax, basis_name, legacy, strux)
   end subroutine run_basis_case

   subroutine run_backend(lmax, basis_name, backend_name, state)
      integer, intent(in) :: lmax
      character(len=*), intent(in) :: basis_name, backend_name
      type(backend_state), intent(out) :: state
      type(control), target :: control_obj
      type(lattice) :: lattice_obj
      character(len=256) :: input_name
      integer :: norb_case

      norb_case = (lmax + 1)**2
      write (input_name, '(a,a,a,a)') 'structure_constants_', trim(basis_name), '_', trim(backend_name)//'.nml'
      call write_fixture(trim(input_name), lmax, backend_name)

      ! control and lattice constructors are the production input boundary.
      ! Keep one lattice alive at a time: restore_to_default allocates the
      ! historical large cluster work arrays before reading ndim from input.
      control_obj = control(trim(input_name))
      call basis_init(lmax)
      lattice_obj = lattice(control_obj)
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      call lattice_obj%structb(.true.)

      state%lmax = lmax
      state%norb = norb_case
      state%kk = lattice_obj%kk
      state%ntot = lattice_obj%ntot
      state%nn_max = lattice_obj%nn_max
      state%sbar_dim = size(lattice_obj%sbar, 1)
      allocate (state%nn, source=lattice_obj%nn)
      allocate (state%sbarvec, source=lattice_obj%sbarvec)
      allocate (state%sbar, source=lattice_obj%sbar(1:norb_case, 1:norb_case, 1:lattice_obj%nn_max, &
                                                    1:lattice_obj%ntot))

      call close_structure_files()
      call delete_file(trim(input_name))
   end subroutine run_backend

   subroutine write_fixture(input_name, lmax, backend_name)
      character(len=*), intent(in) :: input_name, backend_name
      integer, intent(in) :: lmax
      integer :: unit

      open (newunit=unit, file=input_name, status='replace', action='write')
      write (unit, '(a)') '&lattice'
      write (unit, '(a)') '  ndim = 20000, rc = 12.0d0, ntype = 1,'
      write (unit, '(a)') '  alat = 2.86120d0, crystal_sym = ''bcc'', wav = 1.40880d0,'
      write (unit, '(a)') '  ct(1) = 3.0d0, r2 = 9.0d0, npe = 49,'
      write (unit, '(a,a,a)') '  strux_backend = ''', trim(backend_name), ''','
      write (unit, '(a)') '  screening = ''manual'', strux_want_sdot = .false., strux_solve_scale = 9.0d0,'
      write (unit, '(a,4(es14.6,1x))') '  screening_alpha = ', alpha_legacy
      write (unit, '(a)') '/'
      write (unit, '(a)') '&atoms'
      write (unit, '(a,a,a)') '  database = ''', trim(database_directory), ''','
      write (unit, '(a)') '  label(1) = ''Fe'''
      write (unit, '(a)') '/'
      write (unit, '(a)') '&control'
      write (unit, '(a,i0,a,i0)') '  calctype = ''B'', nsp = 1, lmax = ', lmax, ', npold = ', (lmax + 1)**2
      write (unit, '(a)') '  recur = ''block'''
      write (unit, '(a)') '/'
      close (unit)
   end subroutine write_fixture

   subroutine compare_states(lmax, basis_name, lhs, rhs)
      integer, intent(in) :: lmax
      character(len=*), intent(in) :: basis_name
      type(backend_state), intent(in) :: lhs, rhs
      real(rp) :: max_diff, relative_diff, reference_norm
      real(rp) :: onsite_diff, near_diff, vector_diff
      integer :: expected_norb

      expected_norb = (lmax + 1)**2
      call require(lhs%norb == expected_norb, trim(basis_name)//' legacy orbital dimension')
      call require(rhs%norb == expected_norb, trim(basis_name)//' strux orbital dimension')
      call require(lhs%kk == rhs%kk, trim(basis_name)//' cluster dimension')
      call require(lhs%ntot == rhs%ntot, trim(basis_name)//' representative-site dimension')
      call require(lhs%nn_max == rhs%nn_max, trim(basis_name)//' stored-neighbour dimension')
      call require(lhs%sbar_dim == rhs%sbar_dim, trim(basis_name)//' allocated Sbar dimension')

      if (all(shape(lhs%nn) == shape(rhs%nn))) then
         call require(all(lhs%nn == rhs%nn), trim(basis_name)//' neighbour mapping')
      else
         call require(.false., trim(basis_name)//' neighbour-map shape')
      end if
      if (all(shape(lhs%sbarvec) == shape(rhs%sbarvec))) then
         vector_diff = maxval(abs(lhs%sbarvec - rhs%sbarvec))
         call require(vector_diff == 0.0_rp, trim(basis_name)//' neighbour-vector mapping')
      else
         vector_diff = huge(1.0_rp)
         call require(.false., trim(basis_name)//' neighbour-vector shape')
      end if
      if (.not. all(shape(lhs%sbar) == shape(rhs%sbar))) then
         call require(.false., trim(basis_name)//' Sbar shape')
         return
      end if

      call require(all(ieee_is_finite(real(lhs%sbar))) .and. all(ieee_is_finite(aimag(lhs%sbar))), &
                    trim(basis_name)//' legacy Sbar finiteness')
      call require(all(ieee_is_finite(real(rhs%sbar))) .and. all(ieee_is_finite(aimag(rhs%sbar))), &
                    trim(basis_name)//' strux Sbar finiteness')
      call require(all(ieee_is_finite(lhs%sbarvec)) .and. all(ieee_is_finite(rhs%sbarvec)), &
                    trim(basis_name)//' neighbour-vector finiteness')

      max_diff = maxval(abs(lhs%sbar - rhs%sbar))
      reference_norm = sqrt(sum(abs(lhs%sbar)**2))
      relative_diff = sqrt(sum(abs(lhs%sbar - rhs%sbar)**2))/max(reference_norm, tiny(1.0_rp))
      onsite_diff = maxval(abs(lhs%sbar(:, :, 1, :) - rhs%sbar(:, :, 1, :)))
      near_diff = 0.0_rp
      if (lhs%nn_max >= 2) near_diff = maxval(abs(lhs%sbar(:, :, 2, :) - rhs%sbar(:, :, 2, :)))

      write (*, '(a,a,es12.4,a,es12.4,a,es12.4)') trim(basis_name), &
         ' Sbar max=', max_diff, ' relative=', relative_diff, ' vector=', vector_diff
      write (*, '(a,a,es12.4,a,es12.4)') trim(basis_name), &
         ' selected blocks: onsite=', onsite_diff, ' near=', near_diff
      call require(max_diff <= abs_tol .and. relative_diff <= rel_tol, trim(basis_name)//' Sbar agreement')
   end subroutine compare_states

   subroutine require(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         write (*, '(a,a)') 'FAIL: ', trim(label)
         failed = .true.
      end if
   end subroutine require

   subroutine close_structure_files()
      integer :: unit, ios

      do unit = 12, 17
         close (unit, iostat=ios)
      end do
   end subroutine close_structure_files

   subroutine delete_file(file_name)
      character(len=*), intent(in) :: file_name
      integer :: unit, ios

      open (newunit=unit, file=file_name, status='old', action='read', iostat=ios)
      if (ios == 0) close (unit, status='delete', iostat=ios)
   end subroutine delete_file

end program test_structure_constants_backends
