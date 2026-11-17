!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> Radial finite-difference utilities used by the spherical XC paths.
module xc_radial_mod
   use precision_mod, only: rp
   implicit none
   private

   public :: radgra
   public :: radial_flux_divergence

contains

   !> Differentiate a function on the logarithmic atomic radial mesh.
   !>
   !> The logarithmic mesh coordinate is x=a*(i-1) and
   !>   r(x) = b*(exp(x)-1),  dr/dx = r+b.
   !> The finite-difference stencil below approximates df/dx.  Dividing by
   !> dr/dx therefore makes the returned quantity df/dr.  Calling radgra a
   !> second time differentiates that result and returns d2f/dr2.
   subroutine radgra(a, b, nr, rofi, f, gradf)
      integer, intent(in) :: nr
      real(rp), intent(in) :: a, b
      real(rp), dimension(nr), intent(in) :: rofi
      real(rp), dimension(nr), intent(in) :: f
      real(rp), dimension(nr), intent(out) :: gradf

      integer :: nm2, i

      if (nr < 8) error stop 'radgra requires at least 8 mesh points'
      if (abs(a) <= tiny(1.0_rp)) error stop 'radgra requires a non-zero mesh step'

      nm2 = nr - 2

      ! Forward differences for the first two points.
      gradf(1) = ((6.d0*f(2) + 20.d0/3.d0*f(4) + 1.2d0*f(6)) &
                  - (2.45d0*f(1) + 7.5d0*f(3) + 3.75d0*f(5) + 1.d0/6.d0*f(7)))/a
      gradf(2) = ((6.d0*f(3) + 20.d0/3.d0*f(5) + 1.2d0*f(7)) &
                  - (2.45d0*f(2) + 7.5d0*f(4) + 3.75d0*f(6) + 1.d0/6.d0*f(8)))/a

      ! Five-point centered differences in the interior.
      do i = 3, nm2
         gradf(i) = ((f(i - 2) + 8.d0*f(i + 1)) - &
                     (8.d0*f(i - 1) + f(i + 2)))/12.d0/a
      end do

      ! Five-point backward differences at the end of the table.
      gradf(nr - 1) = (-1.d0/12.d0*f(nr - 4) + 0.5d0*f(nr - 3) - &
                       1.5d0*f(nr - 2) + 5.d0/6.d0*f(nr - 1) + 0.25d0*f(nr))/a
      gradf(nr) = (0.25d0*f(nr - 4) - 4.d0/3.d0*f(nr - 3) + 3.d0*f(nr - 2) &
                   - 4.d0*f(nr - 1) + 25.d0/12.d0*f(nr))/a

      ! Convert the derivative with respect to x to the derivative with
      ! respect to the physical radius r.
      do i = 1, nr
         gradf(i) = gradf(i)/(rofi(i) + b)
      end do
   end subroutine radgra

   !> Compute div(F(r) rhat) from a radial flux F.
   !>
   !> For r>0 this evaluates d(r**2 F)/dr / r**2.  At the origin a regular
   !> spherical vector has F(r)=F'(0)r+O(r**3), so the analytic limit 3F'(0)
   !> is used.  No artificial origin radius is introduced.
   subroutine radial_flux_divergence(a, b, rofi, flux, divergence)
      real(rp), intent(in) :: a, b
      real(rp), dimension(:), intent(in) :: rofi, flux
      real(rp), dimension(:), intent(out) :: divergence

      integer :: i, nr
      real(rp), allocatable :: weighted_flux(:), dweighted_flux(:), dflux(:)

      nr = size(rofi)
      if (size(flux) /= nr .or. size(divergence) /= nr) then
         error stop 'radial_flux_divergence array sizes do not agree'
      end if

      allocate (weighted_flux(nr), dweighted_flux(nr), dflux(nr))
      weighted_flux = rofi*rofi*flux
      call radgra(a, b, nr, rofi, weighted_flux, dweighted_flux)
      call radgra(a, b, nr, rofi, flux, dflux)

      ! The first mesh point is the extrapolated r=0 point.  Use the regular
      ! radial limit rather than dividing by zero or inserting an epsilon.
      divergence(1) = 3.d0*dflux(1)
      do i = 2, nr
         divergence(i) = dweighted_flux(i)/(rofi(i)*rofi(i))
      end do
   end subroutine radial_flux_divergence

end module xc_radial_mod
