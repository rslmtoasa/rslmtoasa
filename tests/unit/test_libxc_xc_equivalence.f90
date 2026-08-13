! TEST-14: fixed-density equivalence checks for the internal and libXC XC paths.
!
! This executable is registered only in a libXC-enabled build.  It calls the
! production legacy XCPOT routine and the production libXC wrapper directly;
! it does not contain a second implementation of either functional.
program test_libxc_xc_equivalence
   use precision_mod, only: rp
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc
   implicit none

   integer, parameter :: sample_count = 4
   ! The legacy PBE path keeps decimal constants while libXC evaluates its
   ! reference implementation; the fixed-density differences are below 2e-6
   ! relative for the cases below.
   real(rp), parameter :: abs_tolerance = 2.0e-8_rp
   real(rp), parameter :: rel_tolerance = 2.0e-6_rp
   real(rp), parameter :: radius = 1.25_rp
   real(rp), parameter :: rho_down(sample_count) = [0.31_rp, 0.40_rp, 0.20_rp, 0.55_rp]
   real(rp), parameter :: rho_up(sample_count) = [0.69_rp, 0.40_rp, 0.80_rp, 0.55_rp]
   real(rp), parameter :: rho_gradient(2) = [0.0_rp, 0.0_rp]
   real(rp), parameter :: rho_laplacian(2) = [0.0_rp, 0.0_rp]

   type(control) :: ctl
   type(xc) :: functional
   logical :: failed
   integer :: isample

   failed = .false.
   call g_logger%init()
   call ctl%restore_to_default()
   ctl%nsp = 2
   ctl%txc = 5
   functional = xc(ctl)

   call require(functional%use_libxc, 'txc=5 enables the libXC production route')
   call require(functional%libxc_nspin == 2, 'libXC route is initialized for two spin channels')
   call require(allocated(functional%libxc_func_id), 'PBE LDA libXC functional IDs are available')
   if (allocated(functional%libxc_func_id)) then
      call require(size(functional%libxc_func_id) == 2, 'PBE LDA route contains exchange and correlation')
      if (size(functional%libxc_func_id) == 2) then
         call require(all(functional%libxc_func_id == [1, 12]), &
            'txc=5 maps to libXC LDA_X plus LDA_C_PW')
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
      call functional%XCPOT(rho_down(isample), rho_up(isample), rho_total, rho_gradient, rho_laplacian, &
         radius, internal_v_down, internal_v_up, internal_exc)
      call functional%xcpot_libxc_wrapper(rho_down(isample), rho_up(isample), rho_total, rho_gradient, rho_laplacian, &
         radius, libxc_v_down, libxc_v_up, libxc_exc)

      call require_close(libxc_exc, internal_exc, 'XC energy density', isample)
      call require_close(libxc_v_down, internal_v_down, 'spin-down XC potential', isample)
      call require_close(libxc_v_up, internal_v_up, 'spin-up XC potential', isample)
      call require_close(libxc_v_up - libxc_v_down, internal_v_up - internal_v_down, &
         'XC spin splitting', isample)
   end subroutine compare_sample

   subroutine require_close(actual, expected, quantity, isample)
      real(rp), intent(in) :: actual, expected
      character(len=*), intent(in) :: quantity
      integer, intent(in) :: isample
      real(rp) :: scale

      scale = max(1.0_rp, abs(expected))
      if (abs(actual - expected) > abs_tolerance + rel_tolerance*scale) then
         write (*, '(a,i0,a,a,es16.8,a,es16.8)') 'FAIL: sample ', isample, ' ', trim(quantity)// &
            ' libXC=', actual, ' internal=', expected
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

end program test_libxc_xc_equivalence
