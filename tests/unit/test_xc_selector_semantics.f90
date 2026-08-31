! Focused tests for the TXC and native libXC-ID namespaces.
program test_xc_selector_semantics
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc
   implicit none

   logical :: failed

   failed = .false.
   call g_logger%init()

   call check_legacy(1, 'Barth-Hedin', 'REFERENCE_EQUIVALENT')
   call check_legacy(3, 'Barth-Hedin-Janak', 'APPROXIMATE_ANALOGUE')
   call check_legacy(4, 'Vosko-Wilk-Nusair', 'REFERENCE_EQUIVALENT')
   call check_legacy(5, 'Perdew-Burke-Enzerhof 96 LDA', 'REFERENCE_EQUIVALENT')
   call check_legacy(7, 'Perdew-Zunger', 'REFERENCE_EQUIVALENT_UNPOLARIZED')
   call check_legacy(8, 'Perdew-Burke-Enzerhof 96 GGA', 'REFERENCE_EQUIVALENT')
   call check_legacy(11, 'Barth-Hedin ASW variant', 'APPROXIMATE_ANALOGUE')

#ifdef HAVE_LIBXC
   call check_libxc(101, [1, 17], 'von Barth & Hedin', 'REFERENCE_EQUIVALENT')
   call check_libxc(104, [1, 9], 'Perdew & Zunger', 'REFERENCE_EQUIVALENT_UNPOLARIZED')
   call check_libxc(105, [1, 12], 'Perdew & Wang', 'REFERENCE_EQUIVALENT')
   call check_libxc(108, [101, 130], 'Perdew, Burke & Ernzerhof', 'REFERENCE_EQUIVALENT')
   ! A direct request remains exactly one native ID; it is never auto-paired.
   call check_libxc(1001, [1], 'Slater exchange', 'NO_EQUIVALENT')
#endif

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   endif
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine check_legacy(selector, expected_name, expected_quality)
      integer, intent(in) :: selector
      character(len=*), intent(in) :: expected_name, expected_quality
      type(control) :: ctl
      type(xc) :: functional
      integer, allocatable :: ids(:)

      call ctl%restore_to_default()
      ctl%nsp = 1
      ctl%txc = selector
      functional = xc(ctl)
      ids = functional%get_libxc_functional_ids()

      call require(.not. functional%use_libxc, 'TXC='//trim(int_string(selector))// &
         ' uses the legacy backend')
      call require(.not. functional%is_libxc_functional(), 'TXC='//trim(int_string(selector))// &
         ' is not a libXC selector')
      call require(.not. allocated(functional%libxc_func_id), 'TXC='//trim(int_string(selector))// &
         ' has no active native libXC IDs')
      call require(size(ids) == 0, 'TXC='//trim(int_string(selector))// &
         ' returns no active native libXC IDs')
      call require(trim(functional%backend_name) == 'legacy RS-LMTO', &
         'TXC='//trim(int_string(selector))//' backend provenance')
      call require(trim(functional%functional_name) == expected_name, &
         'TXC='//trim(int_string(selector))//' functional provenance')
      call require(trim(functional%mapping_quality) == expected_quality, &
         'TXC='//trim(int_string(selector))//' reference mapping quality')
   end subroutine check_legacy

#ifdef HAVE_LIBXC
   subroutine check_libxc(selector, expected_ids, name_fragment, expected_quality)
      integer, intent(in) :: selector, expected_ids(:)
      character(len=*), intent(in) :: name_fragment, expected_quality
      type(control) :: ctl
      type(xc) :: functional
      integer, allocatable :: ids(:)

      call ctl%restore_to_default()
      ctl%nsp = 2
      ctl%txc = selector
      functional = xc(ctl)
      ids = functional%get_libxc_functional_ids()

      call require(functional%use_libxc, 'TXC='//trim(int_string(selector))// &
         ' enables the libXC backend')
      call require(functional%is_libxc_functional(), 'TXC='//trim(int_string(selector))// &
         ' is a libXC selector')
      call require(allocated(functional%libxc_func_id), 'TXC='//trim(int_string(selector))// &
         ' has active native libXC IDs')
      call require(size(ids) == size(expected_ids), 'TXC='//trim(int_string(selector))// &
         ' native libXC ID count')
      if (size(ids) == size(expected_ids)) then
         call require(all(ids == expected_ids), 'TXC='//trim(int_string(selector))// &
            ' native libXC IDs')
      endif
      call require(trim(functional%backend_name) == 'libXC', &
         'TXC='//trim(int_string(selector))//' backend provenance')
      call require(index(functional%functional_name, trim(name_fragment)) > 0, &
         'TXC='//trim(int_string(selector))//' libXC name provenance')
      call require(trim(functional%mapping_quality) == expected_quality, &
         'TXC='//trim(int_string(selector))//' alias mapping quality')
   end subroutine check_libxc
#endif

   function int_string(value) result(text)
      integer, intent(in) :: value
      character(len=16) :: text
      write (text, '(i0)') value
   end function int_string

   subroutine require(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description
      if (.not. condition) then
         write (*, '(a)') 'FAIL: '//trim(description)
         failed = .true.
      endif
   end subroutine require

end program test_xc_selector_semantics
