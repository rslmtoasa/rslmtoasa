! XCR-03: fixed-density gradient audit for the legacy PBE route.
program test_legacy_pbe_gga
   use precision_mod, only: rp
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc
   use xc_radial_mod, only: radgra
   implicit none

   integer, parameter :: nr = 160
   real(rp), parameter :: a = log(1.0_rp + 5.0_rp/0.02_rp)/real(nr - 1, rp)
   real(rp), parameter :: b = 0.02_rp, rmax = 5.0_rp
   real(rp), parameter :: pi = 3.1415926535897932384626433832795_rp

   type(control) :: ctl
   type(xc) :: functional
   real(rp), allocatable :: rofi(:), rho_up(:), rho_down(:), drho_up(:), drho_down(:)
   real(rp), allocatable :: d2rho_up(:), d2rho_down(:), perturb_up(:), perturb_down(:)
   real(rp), allocatable :: rho_plus_up(:), rho_plus_down(:), rho_minus_up(:), rho_minus_down(:)
   real(rp), allocatable :: drho_plus_up(:), drho_plus_down(:), drho_minus_up(:), drho_minus_down(:)
   real(rp), allocatable :: d2rho_plus_up(:), d2rho_plus_down(:), d2rho_minus_up(:), d2rho_minus_down(:)
   real(rp), allocatable :: v_up(:), v_down(:), exc(:)
   real(rp) :: rhs, energy_plus, energy_minus, derivative, error_value
   real(rp), parameter :: epsilon_values(4) = [1.e-2_rp, 3.e-3_rp, 1.e-3_rp, 3.e-4_rp]
   logical :: failed
   integer :: i, iepsilon

   failed = .false.
   call g_logger%init()
   call ctl%restore_to_default()
   ctl%nsp = 2
   ctl%txc = 8
   functional = xc(ctl)

   allocate (rofi(nr), rho_up(nr), rho_down(nr), drho_up(nr), drho_down(nr))
   allocate (d2rho_up(nr), d2rho_down(nr), perturb_up(nr), perturb_down(nr))
   allocate (rho_plus_up(nr), rho_plus_down(nr), rho_minus_up(nr), rho_minus_down(nr))
   allocate (drho_plus_up(nr), drho_plus_down(nr), drho_minus_up(nr), drho_minus_down(nr))
   allocate (d2rho_plus_up(nr), d2rho_plus_down(nr), d2rho_minus_up(nr), d2rho_minus_down(nr))
   allocate (v_up(nr), v_down(nr), exc(nr))

   do i = 1, nr
      rofi(i) = b*(exp(a*real(i - 1, rp)) - 1.0_rp)
      rho_up(i) = 0.05_rp + 0.45_rp*exp(-0.35_rp*rofi(i)*rofi(i))
      rho_down(i) = 0.04_rp + 0.30_rp*exp(-0.50_rp*rofi(i)*rofi(i))
      perturb_up(i) = 0.02_rp*rofi(i)**2*exp(-0.20_rp*rofi(i)**2)* &
         (1.0_rp - rofi(i)/rmax)**4
      perturb_down(i) = -0.015_rp*rofi(i)**2*exp(-0.15_rp*rofi(i)**2)* &
         (1.0_rp - rofi(i)/rmax)**4
   end do
   call evaluate(functional, rho_up, rho_down, drho_up, drho_down, d2rho_up, d2rho_down, v_up, v_down, exc)
   rhs = radial_integral(rofi, v_up*perturb_up + v_down*perturb_down)
   write (*, '(a,es16.8)') 'legacy PBE functional derivative RHS: ', rhs

   do iepsilon = 1, size(epsilon_values)
      rho_plus_up = rho_up + epsilon_values(iepsilon)*perturb_up
      rho_plus_down = rho_down + epsilon_values(iepsilon)*perturb_down
      rho_minus_up = rho_up - epsilon_values(iepsilon)*perturb_up
      rho_minus_down = rho_down - epsilon_values(iepsilon)*perturb_down
      call evaluate(functional, rho_plus_up, rho_plus_down, drho_plus_up, drho_plus_down, &
                    d2rho_plus_up, d2rho_plus_down, v_up, v_down, exc)
      energy_plus = radial_integral(rofi, (rho_plus_up + rho_plus_down)*exc)
      call evaluate(functional, rho_minus_up, rho_minus_down, drho_minus_up, drho_minus_down, &
                    d2rho_minus_up, d2rho_minus_down, v_up, v_down, exc)
      energy_minus = radial_integral(rofi, (rho_minus_up + rho_minus_down)*exc)
      derivative = (energy_plus - energy_minus)/(2.0_rp*epsilon_values(iepsilon))
      error_value = abs(derivative - rhs)
      write (*, '(a,es10.3,a,es16.8,a,es10.3)') 'epsilon=', epsilon_values(iepsilon), &
         ' derivative=', derivative, ' error=', error_value
      if (iepsilon == size(epsilon_values)) then
         call require(error_value < 2.e-4_rp, 'legacy PBE gradient path passes fixed-density derivative test')
      end if
   end do

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine evaluate(functional, rho_up, rho_down, drho_up, drho_down, d2rho_up, d2rho_down, &
                       v_up, v_down, exc)
      type(xc), intent(in) :: functional
      real(rp), intent(in) :: rho_up(:), rho_down(:)
      real(rp), intent(out) :: drho_up(:), drho_down(:), d2rho_up(:), d2rho_down(:), v_up(:), v_down(:), exc(:)
      real(rp) :: n(2), nd(2), ndd(2)
      integer :: ir

      call radgra(a, b, nr, rofi, rho_up, drho_up)
      call radgra(a, b, nr, rofi, rho_down, drho_down)
      drho_up(1) = 0.0_rp
      drho_down(1) = 0.0_rp
      call radgra(a, b, nr, rofi, drho_up, d2rho_up)
      call radgra(a, b, nr, rofi, drho_down, d2rho_down)
      do ir = 1, nr
         n = [rho_up(ir), rho_down(ir)]
         nd = [drho_up(ir), drho_down(ir)]
         ndd = [d2rho_up(ir), d2rho_down(ir)]
         call functional%PBEGGA(1, 1, n, nd, ndd, rofi(ir), exc(ir), v_up(ir), v_down(ir))
      end do
   end subroutine evaluate

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

   subroutine require(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description
      if (.not. condition) then
         write (*, '(a)') 'FAIL: '//trim(description)
         failed = .true.
      end if
   end subroutine require

end program test_legacy_pbe_gga
