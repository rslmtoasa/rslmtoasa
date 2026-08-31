! XCF-08: direct components and explicit compositions remain exact and are
! never completed by silently adding a partner.
program test_xc_composition_semantics
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use iso_c_binding, only: c_size_t
   use precision_mod, only: rp
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc, LIBXC_FAMILY_GGA
   use xc_f03_lib_m
   implicit none

   type(control) :: ctl
   type(xc) :: functional
   integer, allocatable :: ids(:)
   real(rp), parameter :: zero_gradient(2) = 0.0_rp
   real(rp), parameter :: zero_laplacian(2) = 0.0_rp
   real(rp), parameter :: a = 0.05_rp
   real(rp), parameter :: b = 0.02_rp
   integer, parameter :: nr = 40
   real(rp), allocatable :: rofi(:), rho_up(:), rho_down(:), zero(:)
   real(rp), allocatable :: v_up(:), v_down(:), exc(:)
   real(rp) :: rho(2), sigma(3), direct_exc(1), direct_vrho(2), direct_vsigma(3)
   integer :: ir
   logical :: failed

   call g_logger%init()
   call ctl%restore_to_default()
   ctl%nsp = 1
   failed = .false.

   ! Direct selection is one native component, even when it is exchange or
   ! correlation only.
   ctl%txc = 1001
   functional = xc(ctl)
   ids = functional%get_libxc_functional_ids()
   call require(size(ids) == 1 .and. ids(1) == 1, 'direct exchange selection is not auto-paired')
   call require(functional%libxc_n_exchange == 1 .and. functional%libxc_n_correlation == 0, &
      'direct exchange composition classification')
   call require(trim(functional%libxc_composition_status) == 'WARNING', &
      'direct exchange selection reports a warning')

   ctl%txc = 1012
   functional = xc(ctl)
   ids = functional%get_libxc_functional_ids()
   call require(size(ids) == 1 .and. ids(1) == 12, 'direct correlation selection is not auto-paired')
   call require(functional%libxc_n_exchange == 0 .and. functional%libxc_n_correlation == 1, &
      'direct correlation composition classification')

   ! A valid, ordinary user-defined composition gets explicit provenance even
   ! when it needs no scientific warning. It is not relabeled as a predefined
   ! bundle merely because its component list happens to match one.
   ctl%txc = 1001
   functional = xc(ctl)
   call functional%set_libxc_component_ids([101, 130])
   call require(functional%libxc_explicit_composition .and. .not. functional%is_predefined_bundle(), &
      'ordinary custom composition retains custom provenance')
   call require(trim(functional%libxc_composition_status) == 'OK', &
      'ordinary custom composition does not emit an unconventional warning')

   ! An explicitly supplied mixed composition is accepted, regardless of
   ! component ordering, and is routed by the complete list.
   ctl%txc = 1001
   functional = xc(ctl)
   call functional%set_libxc_component_ids([12, 117])
   ids = functional%get_libxc_functional_ids()
   call require(size(ids) == 2 .and. ids(1) == 12 .and. ids(2) == 117, &
      'explicit component order is retained')
   call require(functional%libxc_family == LIBXC_FAMILY_GGA .and. functional%libxc_requires_gradients, &
      'mixed LDA plus GGA composition uses the radial route')
   call require(functional%libxc_n_exchange == 1 .and. functional%libxc_n_correlation == 1, &
      'mixed composition counts components by kind')
   call require(functional%libxc_explicit_composition .and. trim(functional%mapping_quality) == 'NO_EQUIVALENT', &
      'explicit composition retains explicit provenance')
   call require(trim(functional%libxc_composition_status) == 'WARNING', &
      'unconventional mixed composition reports a warning')

   allocate(rofi(nr), rho_up(nr), rho_down(nr), zero(nr), v_up(nr), v_down(nr), exc(nr))
   do ir = 1, nr
      rofi(ir) = b*(exp(a*real(ir - 1, rp)) - 1.0_rp)
      rho_up(ir) = 0.37_rp
      rho_down(ir) = 0.11_rp
   end do
   zero = 0.0_rp
   call functional%xcpot_libxc_gga_radial(a, b, rofi, rho_up, rho_down, zero, zero, v_up, v_down, exc)
   ir = nr/2
   rho = [rho_up(ir), rho_down(ir)]
   sigma = 0.0_rp
   call direct_sum(ids, rho, sigma, direct_exc, direct_vrho, direct_vsigma)
   call require(ieee_is_finite(exc(ir)) .and. ieee_is_finite(v_up(ir)) .and. ieee_is_finite(v_down(ir)), &
      'explicit composition produces finite output')
   call require_close(exc(ir), 2.0_rp*direct_exc(1), 'explicit composition sums exact energy')
   call require_close(v_up(ir), 2.0_rp*direct_vrho(1), 'explicit composition sums exact up potential')
   call require_close(v_down(ir), 2.0_rp*direct_vrho(2), 'explicit composition sums exact down potential')

   deallocate(rofi, rho_up, rho_down, zero, v_up, v_down, exc)
   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine direct_sum(component_ids, rho, sigma, exc, vrho, vsigma)
      integer, intent(in) :: component_ids(:)
      real(rp), intent(in) :: rho(2), sigma(3)
      real(rp), intent(out) :: exc(1), vrho(2), vsigma(3)
      type(xc_f03_func_t) :: native
      real(rp) :: exc_tmp(1), vrho_tmp(2), vsigma_tmp(3)
      integer :: i, family

      exc = 0.0_rp
      vrho = 0.0_rp
      vsigma = 0.0_rp
      do i = 1, size(component_ids)
         call xc_f03_func_init(native, component_ids(i), 2)
         family = xc_f03_func_info_get_family(xc_f03_func_get_info(native))
         select case (family)
         case (1)
            call xc_f03_lda_exc_vxc(native, 1_c_size_t, rho, exc_tmp, vrho_tmp)
            exc = exc + exc_tmp
            vrho = vrho + vrho_tmp
         case (2)
            call xc_f03_gga_exc_vxc(native, 1_c_size_t, rho, sigma, exc_tmp, vrho_tmp, vsigma_tmp)
            exc = exc + exc_tmp
            vrho = vrho + vrho_tmp
            vsigma = vsigma + vsigma_tmp
         end select
         call xc_f03_func_end(native)
      end do
   end subroutine direct_sum

   subroutine require(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description
      if (.not. condition) then
         write (*, '(a)') 'FAIL: '//trim(description)
         failed = .true.
      end if
   end subroutine require

   subroutine require_close(actual, expected, description)
      real(rp), intent(in) :: actual, expected
      character(len=*), intent(in) :: description
      if (abs(actual - expected) > 2.0e-11_rp*max(1.0_rp, abs(expected))) then
         write (*, '(a,2(1x,es16.8))') 'FAIL: '//trim(description), actual, expected
         failed = .true.
      end if
   end subroutine require_close

end program test_xc_composition_semantics
