! TEST-14: fixed-density baseline checks for the internal and libXC XC paths.
!
! This executable is registered only in a libXC-enabled build.  It calls the
! production legacy XCPOT routine and the production libXC wrapper directly;
! it does not contain a second implementation of either functional.  The two
! routes are checked independently against their own reference values.
program test_libxc_xc_baseline
   use precision_mod, only: rp
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc
   implicit none

   integer, parameter :: sample_count = 4
   ! These tolerances allow small compiler/libXC-version variation while
   ! catching changes in either production path and their fixed-density
   ! equivalence.
   real(rp), parameter :: abs_tolerance = 2.0e-7_rp
   real(rp), parameter :: rel_tolerance = 2.0e-6_rp
   real(rp), parameter :: radius = 1.25_rp
   real(rp), parameter :: rho_down(sample_count) = [0.31_rp, 0.40_rp, 0.20_rp, 0.55_rp]
   real(rp), parameter :: rho_up(sample_count) = [0.69_rp, 0.40_rp, 0.80_rp, 0.55_rp]
   real(rp), parameter :: rho_gradient(2) = [0.0_rp, 0.0_rp]
   real(rp), parameter :: rho_laplacian(2) = [0.0_rp, 0.0_rp]
   ! Generated from the production paths with libXC 5.2.3; internal and libXC
   ! values intentionally have separate baselines.
   real(rp), parameter :: internal_exc_reference(sample_count) = [ &
      -1.6606541999070092_rp, -1.5099683469349781_rp, &
      -1.7240724982621929_rp, -1.6687764735993575_rp]
   real(rp), parameter :: internal_v_down_reference(sample_count) = [ &
      -1.8821860386718519_rp, -1.9833919038375623_rp, &
      -1.6962292377842685_rp, -2.1936272288988179_rp]
   real(rp), parameter :: internal_v_up_reference(sample_count) = [ &
      -2.3208260903616056_rp, -1.9833919038375623_rp, &
      -2.4158688251756741_rp, -2.1936272288988179_rp]
   real(rp), parameter :: libxc_exc_reference(sample_count) = [ &
      -1.6606546889583951_rp, -1.5099688309059647_rp, &
      -1.7240729202852934_rp, -1.6687769941351385_rp]
   real(rp), parameter :: libxc_v_down_reference(sample_count) = [ &
      -1.8821868284272996_rp, -1.9833925006991866_rp, &
      -1.6962306928293545_rp, -2.1936278661632440_rp]
   real(rp), parameter :: libxc_v_up_reference(sample_count) = [ &
      -2.3208266062120977_rp, -1.9833925006991866_rp, &
      -2.4158691110126762_rp, -2.1936278661632440_rp]

   type(control) :: ctl
   type(xc) :: internal_functional, libxc_functional
   logical :: failed
   integer :: isample

   failed = .false.
   call g_logger%init()
   call ctl%restore_to_default()
   ctl%nsp = 2
   ctl%txc = 5
   internal_functional = xc(ctl)
   ctl%txc = 105
   libxc_functional = xc(ctl)

   call require(.not. internal_functional%use_libxc, 'txc=5 stays on the legacy production route')
   call require(libxc_functional%use_libxc, 'txc=105 enables the explicit libXC production route')
   call require(libxc_functional%libxc_nspin == 2, 'libXC route is initialized for two spin channels')
   call require(allocated(libxc_functional%libxc_func_id), 'PBE LDA libXC functional IDs are available')
   if (allocated(libxc_functional%libxc_func_id)) then
      call require(size(libxc_functional%libxc_func_id) == 2, 'PBE LDA route contains exchange and correlation')
      if (size(libxc_functional%libxc_func_id) == 2) then
         call require(all(libxc_functional%libxc_func_id == [1, 12]), &
            'txc=105 maps to native libXC LDA_X plus LDA_C_PW')
      end if
   end if

   do isample = 1, sample_count
      call compare_sample(isample)
   end do

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   subroutine compare_sample(isample)
      integer, intent(in) :: isample
      real(rp) :: rho_total
      real(rp) :: internal_exc, internal_v_down, internal_v_up
      real(rp) :: libxc_exc, libxc_v_down, libxc_v_up

      rho_total = rho_down(isample) + rho_up(isample)
      call internal_functional%XCPOT(rho_down(isample), rho_up(isample), rho_total, rho_gradient, rho_laplacian, &
         radius, internal_v_down, internal_v_up, internal_exc)
      call libxc_functional%xcpot_libxc_wrapper(rho_down(isample), rho_up(isample), rho_total, rho_gradient, rho_laplacian, &
         radius, libxc_v_down, libxc_v_up, libxc_exc)

      call require_close(internal_exc, internal_exc_reference(isample), &
         'internal XC energy density', isample)
      call require_close(internal_v_down, internal_v_down_reference(isample), &
         'internal spin-down XC potential', isample)
      call require_close(internal_v_up, internal_v_up_reference(isample), &
         'internal spin-up XC potential', isample)
      call require_close(internal_v_up - internal_v_down, &
         internal_v_up_reference(isample) - internal_v_down_reference(isample), &
         'internal XC spin splitting', isample)
      call require_close(libxc_exc, libxc_exc_reference(isample), &
         'libXC XC energy density', isample)
      call require_close(libxc_v_down, libxc_v_down_reference(isample), &
         'libXC spin-down XC potential', isample)
      call require_close(libxc_v_up, libxc_v_up_reference(isample), &
         'libXC spin-up XC potential', isample)
      call require_close(libxc_v_up - libxc_v_down, &
         libxc_v_up_reference(isample) - libxc_v_down_reference(isample), &
         'libXC XC spin splitting', isample)
      call require_close(libxc_exc, internal_exc, &
         'internal/libXC XC energy equivalence', isample)
      call require_close(libxc_v_down, internal_v_down, &
         'internal/libXC spin-down XC potential equivalence', isample)
      call require_close(libxc_v_up, internal_v_up, &
         'internal/libXC spin-up XC potential equivalence', isample)
   end subroutine compare_sample

   subroutine require_close(actual, expected, quantity, isample)
      real(rp), intent(in) :: actual, expected
      character(len=*), intent(in) :: quantity
      integer, intent(in) :: isample
      real(rp) :: scale

      scale = max(1.0_rp, abs(expected))
      if (abs(actual - expected) > abs_tolerance + rel_tolerance*scale) then
         write (*, '(a,i0,a,a,es16.8,a,es16.8)') 'FAIL: sample ', isample, ' ', trim(quantity)// &
            ' actual=', actual, ' reference=', expected
         failed = .true.
      end if
   end subroutine require_close

   subroutine require(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description

      if (.not. condition) then
         write (*, '(a,a)') 'FAIL: ', trim(description)
         failed = .true.
      end if
   end subroutine require

end program test_libxc_xc_baseline
