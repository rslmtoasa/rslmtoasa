! XCR-03: compare the legacy and libXC PBE radial potentials on one mesh.
program test_pbe_gga_comparison
   use precision_mod, only: rp
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc
   use xc_radial_mod, only: radgra
   implicit none

   integer, parameter :: nr = 160
   real(rp), parameter :: a = log(1.0_rp + 5.0_rp/0.02_rp)/real(nr - 1, rp)
   real(rp), parameter :: b = 0.02_rp, pi = 3.1415926535897932384626433832795_rp

   type(control) :: ctl
   type(xc) :: legacy, libxc
   real(rp), allocatable :: rofi(:), rho_up(:), rho_down(:), drho_up(:), drho_down(:)
   real(rp), allocatable :: d2rho_up(:), d2rho_down(:)
   real(rp), allocatable :: legacy_up(:), legacy_down(:), legacy_exc(:)
   real(rp), allocatable :: libxc_up(:), libxc_down(:), libxc_exc(:)
   real(rp) :: legacy_energy, libxc_energy, max_exc, max_scalar, max_magnetic
   real(rp) :: n(2), nd(2), ndd(2)
   logical :: failed
   integer :: i

   failed = .false.
   call g_logger%init()
   call ctl%restore_to_default()
   ctl%nsp = 2
   ctl%txc = 8
   legacy = xc(ctl)
   ctl%txc = 108
   libxc = xc(ctl)

   allocate (rofi(nr), rho_up(nr), rho_down(nr), drho_up(nr), drho_down(nr))
   allocate (d2rho_up(nr), d2rho_down(nr))
   allocate (legacy_up(nr), legacy_down(nr), legacy_exc(nr))
   allocate (libxc_up(nr), libxc_down(nr), libxc_exc(nr))

   do i = 1, nr
      rofi(i) = b*(exp(a*real(i - 1, rp)) - 1.0_rp)
      rho_up(i) = 0.05_rp + 0.45_rp*exp(-0.35_rp*rofi(i)*rofi(i))
      rho_down(i) = 0.04_rp + 0.30_rp*exp(-0.50_rp*rofi(i)*rofi(i))
   end do

   call radgra(a, b, nr, rofi, rho_up, drho_up)
   call radgra(a, b, nr, rofi, rho_down, drho_down)
   drho_up(1) = 0.0_rp
   drho_down(1) = 0.0_rp
   call radgra(a, b, nr, rofi, drho_up, d2rho_up)
   call radgra(a, b, nr, rofi, drho_down, d2rho_down)

   do i = 1, nr
      n = [rho_up(i), rho_down(i)]
      nd = [drho_up(i), drho_down(i)]
      ndd = [d2rho_up(i), d2rho_down(i)]
      call legacy%PBEGGA(1, 1, n, nd, ndd, rofi(i), legacy_exc(i), legacy_up(i), legacy_down(i))
   end do
   call libxc%xcpot_libxc_gga_radial(a, b, rofi, rho_up, rho_down, drho_up, drho_down, &
                                     libxc_up, libxc_down, libxc_exc)

   legacy_energy = radial_integral(rofi, (rho_up + rho_down)*legacy_exc)
   libxc_energy = radial_integral(rofi, (rho_up + rho_down)*libxc_exc)
   max_exc = maxval(abs(legacy_exc - libxc_exc))
   max_scalar = maxval(abs(0.5_rp*(legacy_up + legacy_down - libxc_up - libxc_down)))
   max_magnetic = maxval(abs(0.5_rp*(legacy_up - legacy_down - libxc_up + libxc_down)))

   write (*, '(a,es16.8)') 'TXC=8 energy:       ', legacy_energy
   write (*, '(a,es16.8)') 'TXC=108 energy:     ', libxc_energy
   write (*, '(a,es16.8)') 'max |exc(8)-exc(108)|:      ', max_exc
   write (*, '(a,es16.8)') 'max |scalar(8)-scalar(108)|: ', max_scalar
   write (*, '(a,es16.8)') 'max |magnetic(8)-magnetic(108)|: ', max_magnetic

   ! The two PBE implementations use the same functional but differ in their
   ! finite-difference and libXC internal details.  These bounds catch channel
   ! swaps, unit errors, and missing radial terms while allowing that roundoff.
   call require(abs(legacy_energy - libxc_energy) < 5.e-5_rp, 'TXC=8 and TXC=108 energies agree')
   call require(max_exc < 2.e-5_rp, 'TXC=8 and TXC=108 energy densities agree')
   call require(max_scalar < 5.e-4_rp, 'TXC=8 and TXC=108 scalar potentials agree')
   call require(max_magnetic < 5.e-4_rp, 'TXC=8 and TXC=108 magnetic potentials agree')

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

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

end program test_pbe_gga_comparison
