! XCR-03: libXC GGA channel, radial-flux, and functional-derivative tests.
program test_libxc_gga_radial
   use iso_c_binding, only: c_size_t
   use precision_mod, only: rp
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc
   use xc_radial_mod, only: radgra
   use xc_f03_lib_m
   implicit none

   integer, parameter :: nr = 160
   real(rp), parameter :: a = log(1.0_rp + 5.0_rp/0.02_rp)/real(nr - 1, rp)
   real(rp), parameter :: b = 0.02_rp, rmax = 5.0_rp
   real(rp), parameter :: pi = 3.1415926535897932384626433832795_rp
   real(rp), parameter :: rho_up0 = 0.45_rp, rho_down0 = 0.30_rp
   real(rp), parameter :: rho_up_floor = 0.05_rp, rho_down_floor = 0.04_rp

   type(control) :: ctl
   type(xc) :: functional
   real(rp), allocatable :: rofi(:), rho_up(:), rho_down(:), drho_up(:), drho_down(:)
   real(rp), allocatable :: v_up(:), v_down(:), exc(:), perturb_up(:), perturb_down(:)
   real(rp), allocatable :: rho_plus_up(:), rho_plus_down(:), rho_minus_up(:), rho_minus_down(:)
   real(rp), allocatable :: drho_plus_up(:), drho_plus_down(:), drho_minus_up(:), drho_minus_down(:)
   real(rp), allocatable :: vtmp_up(:), vtmp_down(:), exctmp(:)
   logical :: failed
   integer :: i

   failed = .false.
   call g_logger%init()
   call ctl%restore_to_default()
   ctl%nsp = 2
   ctl%txc = 108
   functional = xc(ctl)

   call require(functional%use_libxc, 'TXC=108 selects libXC')
   call require(functional%libxc_family == 2, 'TXC=108 is a GGA route')

   allocate (rofi(nr), rho_up(nr), rho_down(nr), drho_up(nr), drho_down(nr))
   allocate (v_up(nr), v_down(nr), exc(nr), perturb_up(nr), perturb_down(nr))
   allocate (rho_plus_up(nr), rho_plus_down(nr), rho_minus_up(nr), rho_minus_down(nr))
   allocate (drho_plus_up(nr), drho_plus_down(nr), drho_minus_up(nr), drho_minus_down(nr))
   allocate (vtmp_up(nr), vtmp_down(nr), exctmp(nr))

   do i = 1, nr
      rofi(i) = b*(exp(a*real(i - 1, rp)) - 1.0_rp)
      rho_up(i) = rho_up_floor + rho_up0*exp(-0.35_rp*rofi(i)*rofi(i))
      rho_down(i) = rho_down_floor + rho_down0*exp(-0.50_rp*rofi(i)*rofi(i))
      perturb_up(i) = 0.02_rp*rofi(i)**2*exp(-0.20_rp*rofi(i)**2)* &
         (1.0_rp - rofi(i)/rmax)**4
      perturb_down(i) = -0.015_rp*rofi(i)**2*exp(-0.15_rp*rofi(i)**2)* &
         (1.0_rp - rofi(i)/rmax)**4
   end do
   call radgra(a, b, nr, rofi, rho_up, drho_up)
   call radgra(a, b, nr, rofi, rho_down, drho_down)
   call functional%xcpot_libxc_gga_radial(a, b, rofi, rho_up, rho_down, drho_up, drho_down, &
                                           v_up, v_down, exc)

   call check_sigma_and_channel_mapping(functional, rofi, rho_up, rho_down, drho_up, drho_down, exc, v_up, v_down)
   call check_functional_derivative(functional, rofi, rho_up, rho_down, drho_up, drho_down, &
                                     perturb_up, perturb_down, v_up, v_down)

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine check_sigma_and_channel_mapping(functional, rofi, rho_up, rho_down, drho_up, drho_down, exc, v_up, v_down)
      type(xc), intent(in) :: functional
      real(rp), intent(in) :: rofi(:), rho_up(:), rho_down(:), drho_up(:), drho_down(:), exc(:), v_up(:), v_down(:)
      real(rp) :: rho(2), sigma(3), exc_direct, vup_direct, vdown_direct
      real(rp) :: exc_tmp(1), vrho_tmp(2), vsigma_tmp(3)
      type(xc_f03_func_t) :: xfunc
      integer :: ir

      ir = nr/2
      rho = [rho_up(ir), rho_down(ir)]
      sigma = [drho_up(ir)**2, drho_up(ir)*drho_down(ir), drho_down(ir)**2]
      call direct_pbe_gga(rho, sigma, exc_direct, vup_direct, vdown_direct)
      call require_close(exc(ir), 2.0_rp*exc_direct, 'libXC sigma ordering and GGA energy', 2.e-12_rp)

      ! With constant asymmetric spin densities and zero gradients, the radial
      ! flux vanishes and the helper must return libXC vrho in up/down order.
      call xcpot_constant_check(functional, rofi, rho(1), rho(2))
   end subroutine check_sigma_and_channel_mapping

   subroutine direct_pbe_gga(rho, sigma, exc_direct, vup_direct, vdown_direct)
      real(rp), intent(in) :: rho(2), sigma(3)
      real(rp), intent(out) :: exc_direct, vup_direct, vdown_direct
      real(rp) :: exc_tmp(1), vrho_tmp(2), vsigma_tmp(3)
      integer :: id

      exc_direct = 0.0_rp
      vup_direct = 0.0_rp
      vdown_direct = 0.0_rp
      do id = 101, 130, 29
         call direct_one_gga(id, rho, sigma, exc_tmp, vrho_tmp, vsigma_tmp)
         exc_direct = exc_direct + exc_tmp(1)
         vup_direct = vup_direct + vrho_tmp(1)
         vdown_direct = vdown_direct + vrho_tmp(2)
      end do
   end subroutine direct_pbe_gga

   subroutine direct_one_gga(id, rho, sigma, exc, vrho, vsigma)
      integer, intent(in) :: id
      real(rp), intent(in) :: rho(2), sigma(3)
      real(rp), intent(out) :: exc(1), vrho(2), vsigma(3)
      type(xc_f03_func_t) :: xfunc

      call xc_f03_func_init(xfunc, id, 2)
      call xc_f03_gga_exc_vxc(xfunc, 1_c_size_t, rho, sigma, exc, vrho, vsigma)
      call xc_f03_func_end(xfunc)
   end subroutine direct_one_gga

   subroutine xcpot_constant_check(functional, rofi, rho_up_value, rho_down_value)
      type(xc), intent(in) :: functional
      real(rp), intent(in) :: rofi(:), rho_up_value, rho_down_value
      real(rp), allocatable :: rho_up_constant(:), rho_down_constant(:), zero_gradient(:)
      real(rp), allocatable :: vup(:), vdown(:), e(:)
      real(rp) :: rho(2), sigma(3), expected_exc, expected_up, expected_down
      integer :: ir

      allocate (rho_up_constant(nr), rho_down_constant(nr), zero_gradient(nr), vup(nr), vdown(nr), e(nr))
      rho_up_constant = rho_up_value
      rho_down_constant = rho_down_value
      zero_gradient = 0.0_rp
      call functional%xcpot_libxc_gga_radial(a, b, rofi, rho_up_constant, rho_down_constant, &
                                              zero_gradient, zero_gradient, vup, vdown, e)
      rho = [rho_up_value, rho_down_value]
      sigma = 0.0_rp
      call direct_pbe_gga(rho, sigma, expected_exc, expected_up, expected_down)
      ir = nr/2
      call require_close(e(ir), 2.0_rp*expected_exc, 'constant-density GGA energy', 2.e-12_rp)
      call require_close(vup(ir), 2.0_rp*expected_up, 'legacy/libXC up-channel mapping', 2.e-11_rp)
      call require_close(vdown(ir), 2.0_rp*expected_down, 'legacy/libXC down-channel mapping', 2.e-11_rp)
      deallocate (rho_up_constant, rho_down_constant, zero_gradient, vup, vdown, e)
   end subroutine xcpot_constant_check

   subroutine check_functional_derivative(functional, rofi, rho_up, rho_down, drho_up, drho_down, &
                                           perturb_up, perturb_down, v_up, v_down)
      type(xc), intent(in) :: functional
      real(rp), intent(in) :: rofi(:), rho_up(:), rho_down(:), drho_up(:), drho_down(:)
      real(rp), intent(in) :: perturb_up(:), perturb_down(:), v_up(:), v_down(:)
      real(rp), parameter :: epsilon_values(4) = [1.e-2_rp, 3.e-3_rp, 1.e-3_rp, 3.e-4_rp]
      real(rp) :: rhs, derivative, error_value
      real(rp) :: energy_plus, energy_minus
      integer :: iepsilon

      rhs = radial_integral(rofi, v_up*perturb_up + v_down*perturb_down)
      write (*, '(a,es16.8)') 'functional derivative RHS: ', rhs
      do iepsilon = 1, size(epsilon_values)
         rho_plus_up = rho_up + epsilon_values(iepsilon)*perturb_up
         rho_plus_down = rho_down + epsilon_values(iepsilon)*perturb_down
         rho_minus_up = rho_up - epsilon_values(iepsilon)*perturb_up
         rho_minus_down = rho_down - epsilon_values(iepsilon)*perturb_down
         call radgra(a, b, nr, rofi, rho_plus_up, drho_plus_up)
         call radgra(a, b, nr, rofi, rho_plus_down, drho_plus_down)
         call radgra(a, b, nr, rofi, rho_minus_up, drho_minus_up)
         call radgra(a, b, nr, rofi, rho_minus_down, drho_minus_down)
         call functional%xcpot_libxc_gga_radial(a, b, rofi, rho_plus_up, rho_plus_down, &
                                                 drho_plus_up, drho_plus_down, vtmp_up, vtmp_down, exctmp)
         energy_plus = radial_integral(rofi, (rho_plus_up + rho_plus_down)*exctmp)
         call functional%xcpot_libxc_gga_radial(a, b, rofi, rho_minus_up, rho_minus_down, &
                                                 drho_minus_up, drho_minus_down, vtmp_up, vtmp_down, exctmp)
         energy_minus = radial_integral(rofi, (rho_minus_up + rho_minus_down)*exctmp)
         derivative = (energy_plus - energy_minus)/(2.0_rp*epsilon_values(iepsilon))
         error_value = abs(derivative - rhs)
         write (*, '(a,es10.3,a,es16.8,a,es10.3)') 'epsilon=', epsilon_values(iepsilon), &
            ' derivative=', derivative, ' error=', error_value
         if (iepsilon == size(epsilon_values)) then
            call require(error_value < 2.e-4_rp, 'libXC radial potential passes finite-difference functional derivative')
         end if
      end do
   end subroutine check_functional_derivative

   function radial_integral(rofi, values) result(integral)
      real(rp), intent(in) :: rofi(:), values(:)
      real(rp) :: integral, weight, drdi
      integer :: ir

      integral = 0.0_rp
      do ir = 1, size(rofi)
         weight = 2.0_rp*(mod(ir + 1, 2) + 1.0_rp)/3.0_rp
         if (ir == 1 .or. ir == size(rofi)) weight = 1.0_rp/3.0_rp
         drdi = a*(rofi(ir) + b)
         integral = integral + weight*drdi*4.0_rp*pi*rofi(ir)**2*values(ir)
      end do
   end function radial_integral

   subroutine require_close(actual, expected, description, tolerance)
      real(rp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: description
      if (abs(actual - expected) > tolerance*max(1.0_rp, abs(expected))) then
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

end program test_libxc_gga_radial
