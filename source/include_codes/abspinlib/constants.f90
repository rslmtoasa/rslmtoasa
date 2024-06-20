!>Physical constants
!! @todo Change to double precision once regression tests are passed
module Constants
   use Parameters
   implicit none
   !.. Scalar parameters
   real(dblprec) :: gama = 1.76D11
   real(dblprec) :: k_bolt = 1.38D-23
   real(dblprec) :: k_bolt_ev = 8.61734D-5
   real(dblprec) :: mub = 9.274D-24
   real(dblprec) :: mu0 = 1.2566D-6
   real(dblprec) :: mry = 2.1799D-21
   real(dblprec) :: hbar = 6.62606D-34
   real(dblprec) :: hbar_mev = 6.58211899D-13
   real(dblprec) :: Joule_ev = 6.24150934D18
   real(dblprec) :: ry_ev = 13.605698066
   real(dblprec) :: ev = 1.602176565D-19
   real(dblprec) :: amu = 1.660539040D-27
   real(dblprec) :: angstrom = 1.0D-10
   real(dblprec), parameter :: pi = 3.141592653589793D0
end module Constants

