! XCF-06: compact production-contract regression for the libXC interface.
!
! The checks in this executable deliberately use fixed analytic densities and
! independently initialized libXC objects.  No material calculation is an
! acceptance oracle for this interface.
program test_libxc_production_contract
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use iso_c_binding, only: c_size_t
   use precision_mod, only: rp
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc, LIBXC_DENSITY_FLOOR, LIBXC_FAMILY_LDA, LIBXC_FAMILY_GGA
   use xc_radial_mod, only: radgra, radial_flux_divergence
   use xc_f03_lib_m
   implicit none

   real(rp), parameter :: pi = 3.1415926535897932384626433832795_rp
   real(rp), parameter :: zero_gradient(2) = 0.0_rp
   real(rp), parameter :: zero_laplacian(2) = 0.0_rp
   logical :: failed

   failed = .false.
   call g_logger%init()
   call check_aliases()
   call check_direct_native_namespace()
   call check_lda_grid_and_exchange_oracle()
   call check_spin_and_units()
   call check_gga_routes_and_components()
   call check_density_floor_and_zero_density()
   call check_object_lifecycle()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine check_aliases()
      call check_alias(101, [1, 17], LIBXC_FAMILY_LDA)
      call check_alias(102, [1, 24], LIBXC_FAMILY_LDA)
      call check_alias(103, [1], LIBXC_FAMILY_LDA)
      call check_alias(104, [1, 9], LIBXC_FAMILY_LDA)
      call check_alias(105, [1, 12], LIBXC_FAMILY_LDA)
      call check_alias(106, [1, 7], LIBXC_FAMILY_LDA)
      call check_alias(107, [1, 5], LIBXC_FAMILY_LDA)
      call check_alias(108, [101, 130], LIBXC_FAMILY_GGA)
      call check_alias(109, [117, 130], LIBXC_FAMILY_GGA)
   end subroutine check_aliases

   subroutine check_alias(selector, expected_ids, expected_family)
      integer, intent(in) :: selector, expected_ids(:), expected_family
      type(control) :: ctl
      type(xc) :: functional
      integer, allocatable :: ids(:)
      integer :: i, family, kind
      character(len=256) :: native_name

      call ctl%restore_to_default()
      ctl%nsp = 2
      ctl%txc = selector
      functional = xc(ctl)
      ids = functional%get_libxc_functional_ids()

      call require(functional%use_libxc, 'alias TXC='//int_string(selector)//' activates libXC')
      call require(size(ids) == size(expected_ids), 'alias TXC='//int_string(selector)//' native ID count')
      if (size(ids) == size(expected_ids)) then
         call require(all(ids == expected_ids), 'alias TXC='//int_string(selector)//' native IDs')
      end if
      call require(functional%libxc_family == expected_family, &
         'alias TXC='//int_string(selector)//' aggregate family route')
      call require(size(functional%libxc_component_family) == size(expected_ids), &
         'alias TXC='//int_string(selector)//' stores all component families')
      do i = 1, size(expected_ids)
         call direct_metadata(expected_ids(i), family, kind, native_name)
         call require(functional%libxc_component_family(i) == family, &
            'alias TXC='//int_string(selector)//' runtime family for component '//int_string(i))
         call require(index(functional%functional_name, trim(native_name)) > 0, &
            'alias TXC='//int_string(selector)//' runtime name provenance')
      end do
   end subroutine check_alias

   subroutine check_direct_native_namespace()
      call check_direct(1001, 1, LIBXC_FAMILY_LDA)
      call check_direct(1012, 12, LIBXC_FAMILY_LDA)
      call check_direct(1101, 101, LIBXC_FAMILY_GGA)
      call check_direct(1130, 130, LIBXC_FAMILY_GGA)
   end subroutine check_direct_native_namespace

   subroutine check_direct(selector, expected_id, expected_family)
      integer, intent(in) :: selector, expected_id, expected_family
      type(control) :: ctl
      type(xc) :: functional
      integer, allocatable :: ids(:)
      integer :: family, kind
      character(len=256) :: native_name

      call ctl%restore_to_default()
      ctl%nsp = 2
      ctl%txc = selector
      functional = xc(ctl)
      ids = functional%get_libxc_functional_ids()
      call direct_metadata(expected_id, family, kind, native_name)

      call require(size(ids) == 1, 'direct TXC='//int_string(selector)//' has one native ID')
      if (size(ids) == 1) call require(ids(1) == expected_id, &
         'direct TXC='//int_string(selector)//' selects exactly TXC-1000')
      call require(functional%libxc_family == expected_family, &
         'direct TXC='//int_string(selector)//' family route')
      call require(functional%libxc_component_family(1) == family, &
         'direct TXC='//int_string(selector)//' runtime family')
      call require(trim(functional%functional_name) == trim(native_name), &
         'direct TXC='//int_string(selector)//' runtime name')
   end subroutine check_direct

   subroutine check_lda_grid_and_exchange_oracle()
      integer, parameter :: lda_ids(6) = [1, 5, 7, 9, 12, 17]
      real(rp), parameter :: rho_up_grid(5) = [0.18_rp, 0.31_rp, 0.44_rp, 0.77_rp, 1.0e-16_rp]
      real(rp), parameter :: rho_down_grid(5) = [0.07_rp, 0.19_rp, 0.11_rp, 0.23_rp, 5.0e-17_rp]
      type(control) :: ctl
      type(xc) :: functional
      real(rp) :: rho(2), direct_exc(1), direct_vrho(2), direct_vup, direct_vdown
      real(rp) :: wrapper_exc, wrapper_vdown, wrapper_vup, total
      integer :: i, j
      logical :: finite_values

      do i = 1, size(lda_ids)
         call ctl%restore_to_default()
         ctl%nsp = 2
         ctl%txc = 1000 + lda_ids(i)
         functional = xc(ctl)
         do j = 1, size(rho_up_grid)
            total = rho_up_grid(j) + rho_down_grid(j)
            rho = [rho_up_grid(j), rho_down_grid(j)]
            call direct_lda(lda_ids(i), rho, direct_exc, direct_vrho)
            call functional%xcpot_libxc_wrapper(rho_down_grid(j), rho_up_grid(j), total, &
               zero_gradient, zero_laplacian, 1.0_rp, wrapper_vdown, wrapper_vup, wrapper_exc)
            finite_values = ieee_is_finite(wrapper_exc) .and. ieee_is_finite(wrapper_vup) .and. &
               ieee_is_finite(wrapper_vdown)
            call require(finite_values, 'native LDA grid output is finite')
            if (j < size(rho_up_grid)) then
               call require_close(wrapper_exc, 2.0_rp*direct_exc(1), 'LDA exc Hartree-to-Ry conversion')
               call require_close(wrapper_vup, 2.0_rp*direct_vrho(1), 'LDA v_up channel conversion')
               call require_close(wrapper_vdown, 2.0_rp*direct_vrho(2), 'LDA v_down channel conversion')
            end if
         end do
      end do

      ! Native LDA exchange is also checked against an independent analytic
      ! oracle, including an asymmetric spin density.
      call ctl%restore_to_default()
      ctl%nsp = 2
      ctl%txc = 1001
      functional = xc(ctl)
      call functional%xcpot_libxc_wrapper(0.11_rp, 0.37_rp, 0.48_rp, zero_gradient, zero_laplacian, &
         1.0_rp, wrapper_vdown, wrapper_vup, wrapper_exc)
      direct_vup = -2.0_rp*(6.0_rp/pi)**(1.0_rp/3.0_rp)*0.37_rp**(1.0_rp/3.0_rp)
      direct_vdown = -2.0_rp*(6.0_rp/pi)**(1.0_rp/3.0_rp)*0.11_rp**(1.0_rp/3.0_rp)
      direct_exc = [-1.5_rp*(6.0_rp/pi)**(1.0_rp/3.0_rp)* &
         (0.37_rp**(4.0_rp/3.0_rp) + 0.11_rp**(4.0_rp/3.0_rp))/0.48_rp]
      call require_close(wrapper_exc, direct_exc(1), 'analytic LDA exchange energy oracle')
      call require_close(wrapper_vup, direct_vup, 'analytic LDA exchange up oracle')
      call require_close(wrapper_vdown, direct_vdown, 'analytic LDA exchange down oracle')
   end subroutine check_lda_grid_and_exchange_oracle

   subroutine check_spin_and_units()
      type(control) :: ctl
      type(xc) :: functional
      real(rp) :: rho(2), exc_direct(1), vrho_direct(2)
      real(rp) :: exc_wrapper, vdown_wrapper, vup_wrapper, b_historical

      call ctl%restore_to_default()
      ctl%nsp = 2
      ctl%txc = 1001
      functional = xc(ctl)
      rho = [0.37_rp, 0.11_rp]
      call direct_lda(1, rho, exc_direct, vrho_direct)
      call functional%xcpot_libxc_wrapper(rho(2), rho(1), sum(rho), zero_gradient, zero_laplacian, 1.0_rp, &
         vdown_wrapper, vup_wrapper, exc_wrapper)

      call require_close(exc_wrapper, 2.0_rp*exc_direct(1), 'spin test exc is converted once')
      call require_close(vup_wrapper, 2.0_rp*vrho_direct(1), 'libXC rho(:,1)=up reaches historical V2/up')
      call require_close(vdown_wrapper, 2.0_rp*vrho_direct(2), 'libXC rho(:,2)=down reaches historical V1/down')
      b_historical = 0.5_rp*(vdown_wrapper - vup_wrapper)
      call require(b_historical > 0.0_rp, 'B_xc=(V_down-V_up)/2 has the expected asymmetric sign')
      call require_close(b_historical, 0.5_rp*(vdown_wrapper - vup_wrapper), &
         'historical B_xc conversion')
   end subroutine check_spin_and_units

   subroutine check_gga_routes_and_components()
      integer, parameter :: nr = 128
      real(rp), parameter :: a = log(1.0_rp + 5.0_rp/0.02_rp)/real(nr - 1, rp)
      real(rp), parameter :: b = 0.02_rp
      integer, parameter :: selectors(4) = [108, 109, 1101, 1130]
      real(rp), allocatable :: rofi(:), rho_up(:), rho_down(:), drho_up(:), drho_down(:)
      real(rp), allocatable :: rho_up_const(:), rho_down_const(:)
      real(rp), allocatable :: vup(:), vdown(:), exc(:), zero(:)
      real(rp), allocatable :: expected_vup(:), expected_vdown(:), expected_exc(:)
      real(rp), allocatable :: expected_flux_up(:), expected_flux_down(:)
      real(rp), allocatable :: expected_div_up(:), expected_div_down(:)
      type(control) :: ctl
      type(xc) :: functional, mixed
      real(rp) :: rho(2), sigma(3), direct_exc(1), direct_vrho(2), direct_vsigma(3)
      real(rp) :: rho_up_value, rho_down_value, grad_up_test, grad_down_test
      integer :: i, ir

      allocate(rofi(nr), rho_up(nr), rho_down(nr), drho_up(nr), drho_down(nr))
      allocate(rho_up_const(nr), rho_down_const(nr))
      allocate(vup(nr), vdown(nr), exc(nr), zero(nr))
      allocate(expected_vup(nr), expected_vdown(nr), expected_exc(nr))
      allocate(expected_flux_up(nr), expected_flux_down(nr), expected_div_up(nr), expected_div_down(nr))
      do ir = 1, nr
         rofi(ir) = b*(exp(a*real(ir - 1, rp)) - 1.0_rp)
         rho_up(ir) = 0.20_rp + 0.13_rp*exp(-0.30_rp*rofi(ir)**2)
         rho_down(ir) = 0.09_rp + 0.08_rp*exp(-0.45_rp*rofi(ir)**2)
      end do
      call radgra(a, b, nr, rofi, rho_up, drho_up)
      call radgra(a, b, nr, rofi, rho_down, drho_down)
      ir = nr/2
      rho = [rho_up(ir), rho_down(ir)]
      sigma = [drho_up(ir)**2, drho_up(ir)*drho_down(ir), drho_down(ir)**2]

      do i = 1, size(selectors)
         call ctl%restore_to_default()
         ctl%nsp = 2
         ctl%txc = selectors(i)
         functional = xc(ctl)
         call functional%xcpot_libxc_gga_radial(a, b, rofi, rho_up, rho_down, drho_up, drho_down, &
            vup, vdown, exc)
         call direct_sum(functional%libxc_func_id, rho, sigma, direct_exc, direct_vrho, direct_vsigma)
         call require(functional%libxc_family == LIBXC_FAMILY_GGA, &
            'GGA selector uses the radial route')
         call require_close(exc(ir), 2.0_rp*direct_exc(1), 'GGA energy density is converted once')

         ! With zero gradients the complete radial helper must reduce to the
         ! direct vrho result, including for exchange-only and correlation-only
         ! native component requests.
         rho_up_value = 0.23_rp
         rho_down_value = 0.17_rp
         rho_up_const = rho_up_value
         rho_down_const = rho_down_value
         zero = 0.0_rp
         call functional%xcpot_libxc_gga_radial(a, b, rofi, rho_up_const, rho_down_const, zero, zero, vup, vdown, exc)
         rho = [rho_up_value, rho_down_value]
         sigma = 0.0_rp
         call direct_sum(functional%libxc_func_id, rho, sigma, direct_exc, direct_vrho, direct_vsigma)
         call require_close(vup(ir), 2.0_rp*direct_vrho(1), 'GGA v_up route and unit conversion')
         call require_close(vdown(ir), 2.0_rp*direct_vrho(2), 'GGA v_down route and unit conversion')
         rho = [rho_up(ir), rho_down(ir)]
         sigma = [drho_up(ir)**2, drho_up(ir)*drho_down(ir), drho_down(ir)**2]
      end do

      ! Exercise an explicit mixed LDA+GGA active combination.  There is no
      ! implicit pairing in the selector layer; this is only a regression of
      ! the evaluator contract for an explicitly populated native-ID list.
      call ctl%restore_to_default()
      ctl%nsp = 2
      ctl%txc = 108
      mixed = xc(ctl)
      mixed%libxc_func_id = [1, 101, 130]

      ! Non-zero gradients exercise the mixed GGA vsigma flux and its radial
      ! divergence, while the LDA component contributes only pointwise
      ! exc/vrho terms.
      call mixed%xcpot_libxc_gga_radial(a, b, rofi, rho_up, rho_down, drho_up, drho_down, vup, vdown, exc)
      do ir = 1, nr
         rho = [rho_up(ir), rho_down(ir)]
         if (ir == 1) then
            grad_up_test = 0.0_rp
            grad_down_test = 0.0_rp
         else
            grad_up_test = drho_up(ir)
            grad_down_test = drho_down(ir)
         end if
         sigma = [grad_up_test**2, grad_up_test*grad_down_test, grad_down_test**2]
         call direct_sum(mixed%libxc_func_id, rho, sigma, direct_exc, direct_vrho, direct_vsigma)
         expected_exc(ir) = 2.0_rp*direct_exc(1)
         expected_flux_up(ir) = 2.0_rp*direct_vsigma(1)*grad_up_test + direct_vsigma(2)*grad_down_test
         expected_flux_down(ir) = 2.0_rp*direct_vsigma(3)*grad_down_test + direct_vsigma(2)*grad_up_test
         expected_vup(ir) = direct_vrho(1)
         expected_vdown(ir) = direct_vrho(2)
      end do
      call radial_flux_divergence(a, b, rofi, expected_flux_up, expected_div_up)
      call radial_flux_divergence(a, b, rofi, expected_flux_down, expected_div_down)
      expected_vup = 2.0_rp*(expected_vup - expected_div_up)
      expected_vdown = 2.0_rp*(expected_vdown - expected_div_down)
      call require(maxval(abs(exc - expected_exc)) < 2.0e-11_rp, &
         'mixed LDA+GGA energy includes all components on the radial mesh')
      call require(maxval(abs(vup - expected_vup)) < 2.0e-9_rp, &
         'mixed LDA+GGA up potential includes vrho and vsigma divergence')
      call require(maxval(abs(vdown - expected_vdown)) < 2.0e-9_rp, &
         'mixed LDA+GGA down potential includes vrho and vsigma divergence')

      ! With zero gradients the same active list reduces to the direct vrho
      ! values, explicitly checking the LDA contribution in both channels.
      rho_up_const = rho_up_value
      rho_down_const = rho_down_value
      zero = 0.0_rp
      call mixed%xcpot_libxc_gga_radial(a, b, rofi, rho_up_const, rho_down_const, zero, zero, vup, vdown, exc)
      ir = nr/2
      rho = [rho_up_value, rho_down_value]
      sigma = 0.0_rp
      call direct_sum(mixed%libxc_func_id, rho, sigma, direct_exc, direct_vrho, direct_vsigma)
      call require(mixed%libxc_family == LIBXC_FAMILY_GGA, 'mixed LDA+GGA route is aggregate GGA')
      call require_close(exc(ir), 2.0_rp*direct_exc(1), 'mixed LDA+GGA constant energy includes all components')
      call require_close(vup(ir), 2.0_rp*direct_vrho(1), 'mixed LDA+GGA up potential includes LDA vrho')
      call require_close(vdown(ir), 2.0_rp*direct_vrho(2), 'mixed LDA+GGA down potential includes LDA vrho')

      deallocate(rofi, rho_up, rho_down, drho_up, drho_down, rho_up_const, rho_down_const, &
         vup, vdown, exc, zero, expected_vup, expected_vdown, expected_exc, expected_flux_up, &
         expected_flux_down, expected_div_up, expected_div_down)
   end subroutine check_gga_routes_and_components

   subroutine check_density_floor_and_zero_density()
      type(control) :: ctl
      type(xc) :: functional
      real(rp) :: energy, energy_below, energy_at_floor, energy_above
      real(rp) :: vdown, vup, below, at_floor, above
      real(rp) :: vdown_below, vup_below, vdown_at_floor, vup_at_floor
      real(rp) :: vdown_above, vup_above

      call ctl%restore_to_default()
      ctl%nsp = 2
      ctl%txc = 1001
      functional = xc(ctl)
      call functional%xcpot_libxc_wrapper(0.0_rp, 0.0_rp, 0.0_rp, zero_gradient, zero_laplacian, 1.0_rp, &
         vdown, vup, energy)
      call require(energy == 0.0_rp .and. vdown == 0.0_rp .and. vup == 0.0_rp, &
         'total-zero-density behavior is exactly zero')

      call functional%xcpot_libxc_wrapper(0.0_rp, 0.4_rp, 0.4_rp, zero_gradient, zero_laplacian, 1.0_rp, &
         vdown, vup, energy)
      call require(ieee_is_finite(energy) .and. ieee_is_finite(vdown) .and. ieee_is_finite(vup), &
         'one-spin-channel-zero limit is finite')
      below = LIBXC_DENSITY_FLOOR*0.5_rp
      at_floor = LIBXC_DENSITY_FLOOR
      above = LIBXC_DENSITY_FLOOR*1.01_rp
      call functional%xcpot_libxc_wrapper(below, 0.4_rp, below + 0.4_rp, zero_gradient, zero_laplacian, 1.0_rp, &
         vdown_below, vup_below, energy_below)
      call functional%xcpot_libxc_wrapper(at_floor, 0.4_rp, at_floor + 0.4_rp, zero_gradient, zero_laplacian, 1.0_rp, &
         vdown_at_floor, vup_at_floor, energy_at_floor)
      call require(ieee_is_finite(energy_at_floor) .and. ieee_is_finite(vdown_at_floor) .and. &
         ieee_is_finite(vup_at_floor), 'density-floor energy and potential are finite')
      call functional%xcpot_libxc_wrapper(above, 0.4_rp, above + 0.4_rp, zero_gradient, zero_laplacian, 1.0_rp, &
         vdown_above, vup_above, energy_above)
      call require(ieee_is_finite(energy_above) .and. ieee_is_finite(vdown_above) .and. &
         ieee_is_finite(vup_above), 'above-floor energy and potential are finite')
      call require(abs(energy_above - energy_at_floor) < 1.0e-5_rp, &
         'density floor does not create an appreciable energy jump')
      call require(abs(energy_below - energy_at_floor) < 1.0e-5_rp, &
         'density floor does not create an appreciable below-threshold jump')
      call require(abs(vdown_below - vdown_at_floor) < 1.0e-5_rp .and. &
         abs(vup_below - vup_at_floor) < 1.0e-5_rp, &
         'density floor does not create an appreciable below-threshold potential jump')
      call require(abs(vdown_above - vdown_at_floor) < 1.0e-5_rp .and. &
         abs(vup_above - vup_at_floor) < 1.0e-5_rp, &
         'density floor does not create an appreciable threshold potential jump')
   end subroutine check_density_floor_and_zero_density

   subroutine check_object_lifecycle()
      type(control) :: ctl
      type(xc) :: functional
      integer :: i, selector

      do i = 1, 24
         selector = merge(1001, 108, mod(i, 2) == 1)
         call ctl%restore_to_default()
         ctl%nsp = 2
         ctl%txc = selector
         functional = xc(ctl)
         call require(functional%use_libxc, 'repeated libXC construction remains active')
         call functional%cleanup_libxc()
         call require(.not. functional%use_libxc, 'explicit libXC cleanup releases metadata')
      end do
   end subroutine check_object_lifecycle

   subroutine direct_metadata(id, family, kind, name)
      integer, intent(in) :: id
      integer, intent(out) :: family, kind
      character(len=*), intent(out) :: name
      type(xc_f03_func_t) :: xfunc
      type(xc_f03_func_info_t) :: info

      call xc_f03_func_init(xfunc, id, 2)
      info = xc_f03_func_get_info(xfunc)
      family = xc_f03_func_info_get_family(info)
      kind = xc_f03_func_info_get_kind(info)
      name = trim(xc_f03_func_info_get_name(info))
      call xc_f03_func_end(xfunc)
   end subroutine direct_metadata

   subroutine direct_lda(id, rho, exc, vrho)
      integer, intent(in) :: id
      real(rp), intent(in) :: rho(2)
      real(rp), intent(out) :: exc(1), vrho(2)
      type(xc_f03_func_t) :: xfunc

      call xc_f03_func_init(xfunc, id, 2)
      call xc_f03_lda_exc_vxc(xfunc, 1_c_size_t, rho, exc, vrho)
      call xc_f03_func_end(xfunc)
   end subroutine direct_lda

   subroutine direct_sum(ids, rho, sigma, exc, vrho, vsigma)
      integer, intent(in) :: ids(:)
      real(rp), intent(in) :: rho(2), sigma(3)
      real(rp), intent(out) :: exc(1), vrho(2), vsigma(3)
      type(xc_f03_func_t) :: xfunc
      type(xc_f03_func_info_t) :: info
      real(rp) :: exc_tmp(1), vrho_tmp(2), vsigma_tmp(3)
      integer :: id, family

      exc = 0.0_rp
      vrho = 0.0_rp
      vsigma = 0.0_rp
      do id = 1, size(ids)
         call xc_f03_func_init(xfunc, ids(id), 2)
         info = xc_f03_func_get_info(xfunc)
         family = xc_f03_func_info_get_family(info)
         if (family == LIBXC_FAMILY_LDA) then
            call xc_f03_lda_exc_vxc(xfunc, 1_c_size_t, rho, exc_tmp, vrho_tmp)
            vsigma_tmp = 0.0_rp
         else if (family == LIBXC_FAMILY_GGA) then
            call xc_f03_gga_exc_vxc(xfunc, 1_c_size_t, rho, sigma, exc_tmp, vrho_tmp, vsigma_tmp)
         else
            call require(.false., 'independent direct evaluator sees only LDA/GGA test components')
            exc_tmp = 0.0_rp
            vrho_tmp = 0.0_rp
            vsigma_tmp = 0.0_rp
         end if
         exc = exc + exc_tmp
         vrho = vrho + vrho_tmp
         vsigma = vsigma + vsigma_tmp
         call xc_f03_func_end(xfunc)
      end do
   end subroutine direct_sum

   subroutine require_close(actual, expected, description)
      real(rp), intent(in) :: actual, expected
      character(len=*), intent(in) :: description
      if (abs(actual - expected) > 3.0e-11_rp*max(1.0_rp, abs(expected))) then
         write (*, '(a,es16.8,a,es16.8)') 'FAIL: '//trim(description)//' actual=', actual, ' expected=', expected
         failed = .true.
      end if
   end subroutine require_close

   subroutine require(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description
      if (.not. condition) then
         write (*, '(a)') 'FAIL: '//trim(description)
         failed = .true.
      end if
   end subroutine require

   function int_string(value) result(text)
      integer, intent(in) :: value
      character(len=16) :: text
      write (text, '(i0)') value
   end function int_string

end program test_libxc_production_contract
