!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!
! PROGRAM: test_simpson_f_temperature
!
!> @brief Finite-temperature Fermi integral of simpson_f against the Sommerfeld
!>        closed form.
!> @details For g(E) = E the integral of g f(E; mu, kT) from a (f = 1 there) to
!>          infinity is (mu^2 - a^2)/2 + (pi^2/6) (kT)^2, with no further terms.
!>          The closed form is written here; kB is the constant simpson_f uses.
!------------------------------------------------------------------------------
program test_simpson_f_temperature
   use math_mod, only: simpson_f, kB_simpson_f
   use precision_mod, only: rp
   implicit none

   integer, parameter :: npts = 2001, n = npts + 9, imu = 1001
   real(rp), parameter :: h = 5.0e-4_rp, emin = -1.0_rp, tol = 1.0e-9_rp
   real(rp), parameter :: pi = 3.14159265358979323846_rp
   real(rp) :: ene(n), y(n), aint, expected, mu, kt, kt_wrong, wrong_expected
   integer :: i, it
   logical :: failed
   real(rp), parameter :: temps(2) = [300.0_rp, 600.0_rp]

   failed = .false.
   do i = 1, n
      ene(i) = emin + h*real(i - 1, rp)
   end do
   y = ene
   mu = ene(imu)

   do it = 1, size(temps)
      kt = kB_simpson_f*temps(it)
      call simpson_f(aint, ene, mu, npts, y, .true., .false., temps(it))
      expected = 0.5_rp*(mu**2 - ene(1)**2) + pi**2/6.0_rp*kt**2
      ! the same closed form with a 1 % wrong kB, to show the margin of the check
      kt_wrong = 1.01_rp*kt
      wrong_expected = 0.5_rp*(mu**2 - ene(1)**2) + pi**2/6.0_rp*kt_wrong**2
      write (*, '(a,f7.1,a,es22.14,a,es22.14,a,es10.2,a,es10.2)') 'T = ', temps(it), ' K  simpson_f ', aint, '  closed form ', expected, &
         '  |diff| ', abs(aint - expected), '  |diff| for kB + 1 % ', abs(aint - wrong_expected)
      if (abs(aint - expected) > tol) failed = .true.
      if (abs(aint - wrong_expected) <= tol) failed = .true.
   end do

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if
end program test_simpson_f_temperature
