!
! Copyright (C) 2017
! Cyrille Barreteau <mailto:cyrille.barreteau@cea.fr>,
! Mathieu Cesar <mailto:mathieu.cesar@cea.fr>,
! Pascal Thibaudeau <mailto:pascal.thibaudeau@cea.fr>.
!
! This software is a computer program whose purpose is DyNaMol.
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and inRIA at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.
!
!  math.f90
!  DyNaMol
!
!
! !  Modified by Ivan Miranda and Ramon Cardias to include more functions.
module math_mod
#if defined(LAPACK95_FOUND)
!  use lapack95, only: getrf, getri
#endif
   use precision_mod, only: ip, rp, longint
   implicit none

   !> type for sorting routine
   type, public :: group
      integer :: order        ! original order of unsorted data
      real(rp) :: value       ! values to be sorted by
   end type group

   ! To do: precise the quantity (charge, energy, etc)
   real(rp), parameter :: epsilon = 1.0E-10_rp

   ! Irrational
   !> \f$ \sqrt{2} \f$
   real(rp), parameter :: sqrt_two = sqrt(2.0_rp)
   !> \f$ 1/\sqrt{2} \f$
   real(rp), parameter :: one_over_sqrt_two = 1.0_rp/sqrt_two
   !> \f$ \sqrt{3} \f$
   real(rp), parameter :: sqrt_three = sqrt(3.0_rp)
   !> \f$ 1/3 \f$
   real(rp), parameter :: one_third = 0.333333333333333333333_rp
   !> \f$ 2/3 \f$
   real(rp), parameter :: two_third = 0.666666666666666666666_rp

   ! Transcendental
   !> \f$ \pi \f$
   real(rp), parameter :: pi = 3.14159265358979323846_rp
   !> \f$ 2\times\pi \f$
   real(rp), parameter :: two_pi = 2.0_rp*pi
   !> \f$ 4\times\pi \f$
   real(rp), parameter :: four_pi = 4.0_rp*pi
   !> \f$ \sqrt{\pi} \f$
   real(rp), parameter :: sqrt_pi = sqrt(pi)
   !> \f$ 1/\sqrt{\pi} \f$
   real(rp), parameter :: one_over_sqrt_pi = 1.0_rp/sqrt_pi
   !> \f$ \sqrt{2\pi} \f$
   real(rp), parameter :: sqrt_two_pi = sqrt_two*sqrt_pi
   !> \f$ 1/\sqrt{2\pi} \f$
   real(rp), parameter :: one_over_sqrt_two_pi = 1.0_rp/sqrt_two_pi
   !> \f$ \pi/180 \f$
   real(rp), parameter :: deg2rad = pi/180.0_rp
   !> \f$ 180/\pi \f$
   real(rp), parameter :: rad2deg = 180.0_rp/pi
   !> \f$ angstrom to atomic units \f$
   real(rp), parameter :: ang2au = 1.8897259886_rp
   !> \f$ Rydberg to eV \f$
   real(rp), parameter :: ry2ev = 13.605703976_rp
   !> \f$ Rydberg to Tesla \f$
   real(rp), parameter :: ry2tesla = 2.35051754997e5_rp
   !> \f$ Gyromagnetic ratio in rad/(s·T) \f$
   real(rp), parameter :: gama = 1.76e11_rp

   ! Screening parameters
   !> Original screening (From Jepsen)
   real(rp), dimension(3), parameter :: qm_canonical = [.348485d0, .053030d0, .010714d0]
   !> Screening from LMTO47 (RWS = average WS for all)
   real(rp), dimension(3), parameter :: qm_lmto47 = [0.3500000d0, 0.0571667d0, 0.0168070d0]

   ! Complex
   !> Imaginary unit \f$ \mathrm{i} \f$
   complex(rp), parameter :: i_unit = cmplx(0.0_rp, 1.0_rp, kind=rp)
   complex(rp), parameter :: cone = (1.0d0, 0.0d0)
   complex(rp), parameter :: cmone = (-1.0d0, 0.0d0)
   complex(rp), parameter :: czero = (0.0d0, 0.0d0)
   !> Identity matrix \f$ \mathbf{I} \f$ (rank 2)
   complex(rp), dimension(2, 2), parameter :: sigma_0 &
                                              = reshape((/1.0_rp, 0.0_rp, 0.0_rp, 1.0_rp/), (/2, 2/))
   ! Pauli matrices
   !> Pauli matrix \f$ \sigma_x \f$
   complex(rp), dimension(2, 2), parameter :: sigma_x &
                                              = reshape((/0.0_rp, 1.0_rp, 1.0_rp, 0.0_rp/), (/2, 2/))
   !> Pauli matrix \f$ \sigma_y \f$
   complex(rp), dimension(2, 2), parameter :: sigma_y &
                                              = reshape((/0.0_rp, 1.0_rp, -1.0_rp, 0.0_rp/), (/2, 2/))*i_unit
   !> Pauli matrix \f$ \sigma_z \f$
   complex(rp), dimension(2, 2), parameter :: sigma_z &
                                              = reshape((/1.0_rp, 0.0_rp, 0.0_rp, -1.0_rp/), (/2, 2/))
   !> Angular momentum operators
   complex(rp), dimension(9, 9), parameter :: L_x &
                                              = reshape((/0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, -1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, -1.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, -1.0_rp, -sqrt_three, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, sqrt_three, 0.0_rp, 0.0_rp, 0.0_rp/), &
                                                        (/9, 9/))*i_unit
   complex(rp), dimension(9, 9), parameter :: L_y &
                                              = reshape((/0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, -1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, -1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, -1.0_rp, sqrt_three, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, -sqrt_three, 0.0_rp, 0.0_rp/), &
                                                        (/9, 9/))*i_unit
   complex(rp), dimension(9, 9), parameter :: L_z &
                                              = reshape((/0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, -1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 2.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, -1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, -2.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                                          0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp/), &
                                                        (/9, 9/))*i_unit

contains
   !> Create vector of all zeros
   function zeros1(n) result(v)
      integer, intent(in) :: n
      real(rp), dimension(n) :: v

      v = 0.0_rp
   end function zeros1

   !> Create matrix of all zeros
   function zeros2(n1, n2) result(m)
      integer, intent(in) :: n1, n2
      real(rp), dimension(n1, n2) :: m

      m = 0.0_rp
   end function zeros2

   !> Create vector of all ones
   function ones1(n) result(v)
      integer, intent(in) :: n
      real(rp), dimension(n) :: v

      v = 1.0_rp
   end function ones1

   !> Create matrix of all ones
   function ones2(n1, n2) result(m)
      integer, intent(in) :: n1, n2
      real(rp), dimension(n1, n2) :: m

      m = 1.0_rp
   end function ones2

   !> Create diagonal matrix
   function diag1(v) result(m)
      real(rp), dimension(:), intent(in) :: v
      real(rp), dimension(size(v), size(v)) :: m
      integer :: i, n

      m = 0.0_rp
      n = size(v)

      do i = 1, n
         m(i, i) = v(i)
      end do
   end function diag1

   !> Get diagonal element of matrix
   function diag2(m) result(v)
      real(rp), dimension(:, :), intent(in) :: m
      real(rp), dimension(min(size(m, 1), size(m, 2))) :: v
      integer :: i, n

      n = min(size(m, 1), size(m, 2))

      do i = 1, n
         v(i) = m(i, i)
      end do
   end function diag2

   !> Factorial function by Ivan Miranda (17.09.2023)
   !> This function takes factorial(n<0) = 0, as necessary for
   !> calculating the Wigner-3j symbol. If you want the conventional
   !> definition, just use the built-in function gamma.
   function fact(n) result(f)
      !
      implicit none
      !
      real(rp), intent(in) :: n
      real(rp) :: f
      !
      if (n .lt. 0) then
         f = 0.0_rp
      else
         f = gamma(n + 1.0_rp)
      end if

   end function fact

   !> Regular double factorial function by Ivan Miranda (18.09.2023)
   !> Can accept real arguments.
   function factorial2(n) result(f)
      !
      implicit none
      !
      real(rp), intent(in) :: n
      real(rp) :: f
      !
      ! Local variables
      real(rp) :: exponent
      !
      exponent = 0.25_rp*(1.0_rp - cos(pi*n))
      f = ((2.0_rp/pi)**(exponent))*(2.0_rp**(n/2.0_rp))*gamma((n/2.0_rp) + 1.0_rp)

   end function factorial2

   !> Calculates the Wigner-3j symbol for the integers (j1, j2, j3, m1, m2, m3)
   !> Implemented by Ivan Miranda (17.09.2023)
   !> Computed here using the Racah formula.
   function wigner3j(j1, j2, j3, m1, m2, m3) result(w)
      !
      implicit none
      !
      integer, intent(in) :: j1, j2, j3, m1, m2, m3
      real(rp) :: w ! final result
      !
      ! Local variables
      integer :: N, K, i
      integer, dimension(3) :: array1, array2
      real(rp) :: part1, part2, part3
      real(rp) :: temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8
      !

      array1(1) = 0; array1(2) = (j2 - j3 - m1); array1(3) = (j1 - j3 + m2)
      array2(1) = (j1 + j2 - j3); array2(2) = (j1 - m1); array2(3) = (j2 + m2)

      K = maxval(array1)
      N = minval(array2)

      part3 = 0.0_rp

      if ((m1 + m2 + m3) .ne. 0) then
         w = 0.0_rp
      else
         temp1 = fact(1.0_rp*real(j1 + j2 - j3))
         temp2 = fact(1.0_rp*real(j1 - j2 + j3))
         temp3 = fact(1.0_rp*real(j2 + j3 - j1))
         temp4 = fact(1.0_rp*real(j1 + j2 + j3 + 1))
         part1 = real((cmone)**(cone*real(j1 - j2 - m3)))*sqrt(temp1*temp2*temp3/temp4)
         temp5 = fact(1.0_rp*real(j1 - m1))*fact(1.0_rp*real(j1 + m1))*fact(1.0_rp*real(j2 - m2))
         temp6 = fact(1.0_rp*real(j2 + m2))*fact(1.0_rp*real(j3 - m3))*fact(1.0_rp*real(j3 + m3))
         part2 = sqrt(temp5*temp6)
         do i = K, N
            temp7 = fact(1.0_rp*real(i))*fact(1.0_rp*real(j1 + j2 - j3 - i))*fact(1.0_rp*real(j1 - m1 - i))
            temp8 = fact(1.0_rp*real(j2 + m2 - i))*fact(1.0_rp*real(j3 - j2 + m1 + i))*fact(1.0_rp*real(j3 - j1 - m2 + i))
            part3 = part3 + real((cmone)**(cone*real(i)))*(1.0_rp/(temp7*temp8))
         end do
         w = part1*part2*part3
      end if

   end function wigner3j

   !> Computes the Gaunt coefficients using the Cruzan-Racah expression
   !> More about: Didier Sébilleau 1998 J. Phys. A: Math. Gen. 31 7157
   !> Implemented by Ivan Miranda on 17.09.2023
   function gaunt(l1, l2, l3, m1, m2, m3) result(g)
      !
      implicit none
      !
      integer, intent(in) :: l1, l2, l3, m1, m2, m3
      real(rp) :: g
      !
      ! Local Variables
      real(rp) :: part1, part2, part3
      real(rp) :: temp1

      part1 = wigner3j(l1, l2, l3, 0, 0, 0)
      part2 = wigner3j(l1, l2, l3, m1, m2, m3)

      temp1 = (2.0_rp*real(l1) + 1.0_rp)*(2.0_rp*real(l2) + 1.0_rp)*(2.0_rp*real(l3) + 1.0_rp)
      part3 = sqrt(temp1/four_pi)

      g = part1*part2*part3

   end function gaunt

   !> Kronecker delta by Ivan Miranda
   function kron(i, j) result(delta)
      !
      implicit none
      !
      integer, intent(in) :: i, j
      real(rp) :: delta
      !
      if (i .eq. j) then
         delta = 1.0_rp
      else
         delta = 0.0_rp
      end if

   end function kron

   !> Computes the Gaunt coefficients, but for real spherical harmonics
   !> Uses Eqs. of Homeier, Journal of Molecular Structure (Theochem) 368 (1996) 31-37
   !> Implemented by Ivan Miranda on 17.09.2023
   !> Assumes the order u1 >= u2 >= u3
   function realgaunt_ord(l1, l2, l3, u1, u2, u3) result(g)
      !
      implicit none
      !
      integer, intent(in) :: l1, l2, l3, u1, u2, u3
      complex(rp) :: g
      !
      ! Local Variables
      complex(rp) :: factor
      !
      factor = (1.0_rp/sqrt(2.0_rp))*cone
      if ((u1 .gt. 0) .and. (u2 .gt. 0) .and. (u3 .gt. 0)) then
         if ((u2 + u3) .eq. u1) then
            g = factor*((cmone)**(cone*u1))*gaunt(l1, l2, l3, -u1, u2, u3)
         else if ((u1 + u3) .eq. u2) then
            g = factor*((cmone)**(cone*u2))*gaunt(l1, l2, l3, u1, -u2, u3)
         else if ((u1 + u2) .eq. u3) then
            g = factor*((cmone)**(cone*u3))*gaunt(l1, l2, l3, u1, u2, -u3)
         else
            g = czero
         end if
      else if ((u1 .gt. 0) .and. (u2 .lt. 0) .and. (u3 .lt. 0)) then
         if ((u1 + u3) .eq. u2) then
            g = factor*((cmone)**(cone*u3))*gaunt(l1, l2, l3, u1, -u2, u3)
         else if ((u1 + u2) .eq. u3) then
            g = factor*((cmone)**(cone*u2))*gaunt(l1, l2, l3, u1, u2, -u3)
         else if ((u1 + u2 + u3) .eq. 0) then
            g = factor*((cmone)**(cone*u1 + cone))*gaunt(l1, l2, l3, u1, u2, u3)
         else
            g = czero
         end if
      else if ((u3 .eq. 0) .and. (u1 .eq. u2) .and. (u1 .gt. 0)) then
         g = ((cmone)**(cone*u1))*gaunt(l1, l2, l3, u1, -u1, 0)*cone
      else if ((u1 .eq. 0) .and. (u2 .eq. u3) .and. (u3 .lt. 0)) then
         g = ((cmone)**(cone*u2))*gaunt(l1, l2, l3, 0, u2, -u2)*cone
      else if ((u1 .eq. 0) .and. (u2 .eq. 0) .and. (u3 .eq. 0)) then
         g = gaunt(l1, l2, l3, 0, 0, 0)*cone
      else
         g = czero
      end if

   end function realgaunt_ord

   !> Computes the Gaunt coefficients, but for real spherical harmonics
   !> Considering any order of u1, u2, u3 and the symmetry properties
   !> of the gaunt function
   function realgaunt(l1, l2, l3, u1, u2, u3) result(g)
      !
      implicit none
      !
      integer, intent(in) :: l1, l2, l3, u1, u2, u3
      complex(rp) :: g
      !
      if ((u1 .ge. u2) .and. (u2 .ge. u3)) then
         g = realgaunt_ord(l1, l2, l3, u1, u2, u3)
      else if ((u1 .ge. u3) .and. (u3 .ge. u2)) then
         g = realgaunt_ord(l1, l3, l2, u1, u3, u2)
      else if ((u2 .ge. u1) .and. (u1 .ge. u3)) then
         g = realgaunt_ord(l2, l1, l3, u2, u1, u3)
      else if ((u2 .ge. u3) .and. (u3 .ge. u1)) then
         g = realgaunt_ord(l2, l3, l1, u2, u3, u1)
      else if ((u3 .ge. u1) .and. (u1 .ge. u2)) then
         g = realgaunt_ord(l3, l1, l2, u3, u1, u2)
      else if ((u3 .ge. u2) .and. (u2 .ge. u1)) then
         g = realgaunt_ord(l3, l2, l1, u3, u2, u1)
      end if

   end function realgaunt

   !> The generalized Newton binomial factor
   !> Supports n as a real number, and it is well defined for k>=0
   !> Implemented by Ivan Miranda on 20.09.2023
   function genew(n, k) result(b)
      !
      implicit none
      !
      integer, intent(in) :: k
      real(rp), intent(in) :: n
      real(rp) :: b
      !
      ! Local Variables
      integer :: i
      real(rp) :: numerator, denominator
      !
      numerator = 1.0_rp
      denominator = 1.0_rp
      !
      do i = 0, (k - 1)
         numerator = numerator*(n - 1.0_rp*i)
         denominator = denominator*(1.0_rp*(k - i))
      end do

      b = numerator/denominator

   end function genew

   !> Calculate the associated Legendre polynomials Plm
   !> Implemented by Ivan Miranda on 20.09.2023
   !> Without the factor (-1)^m in the front
   function associated_lp(x, l, m) result(a)
      !
      implicit none
      !
      real(rp), intent(in) :: x
      integer, intent(in) :: l, m
      real(rp) :: a
      !
      ! Local Variables
      integer :: i, new_l, new_m
      real(rp) :: part1, part2, part3
      real(rp) :: temp1, temp2, temp3
      !
      if (abs(m) .gt. abs(l)) then
         a = 0.0_rp
      else
         part1 = 0.0_rp
         new_l = l
         new_m = m
         if (l .lt. 0) then; new_l = abs(l) - 1; end if
         if (m .lt. 0) then; new_m = abs(m); end if
         do i = new_m, new_l
            temp1 = gamma(1.0_rp*i + 1.0_rp)/gamma(1.0_rp*(i - new_m) + 1.0_rp)
            temp2 = x**(1.0_rp*(i - new_m))
            temp3 = genew(0.5_rp*(new_l + i - 1), new_l)*genew(1.0_rp*new_l, i)
            part1 = part1 + (temp1*temp2*temp3)
         end do
         part2 = (2.0_rp**(1.0_rp*new_l))
         part3 = (1.0_rp - x**(2.0_rp))**(0.5_rp*new_m)
         a = part1*part2*part3
         if (m .lt. 0) then
            a = a*(gamma(1.0_rp*(new_l - new_m) + 1.0_rp)/gamma(1.0_rp*(new_l + new_m) + 1.0_rp))
            a = a*real((cmone)**(cone*new_m))
         end if
      end if

   end function associated_lp

   !> Computes the complex spherical harmonics
   !> for a given unit vector in cartesian coordinates (v_x, v_y, v_z)
   !> (if the vector is not unitary, the function forces it)
   !> Implemented by Ivan Miranda on 20.09.2023
   function complex_spharm(v, l, m) result(s)
      !
      implicit none
      !
      integer, intent(in) :: l, m
      real(rp), dimension(3), intent(in) :: v ! Input in cartesian coordinates (vx, vy, vz)
      complex(rp) :: s
      !
      ! Local Variables
      real(rp), dimension(3) :: unit_v, unit_sph
      real(rp) :: temp1, temp2
      complex(rp) :: part1, part2, part3
      !
      if (abs(norm2(v)) .eq. 0) then
         unit_sph(:) = 0.0_rp
      else
         unit_v(1) = v(1)/norm2(v)
         unit_v(2) = v(2)/norm2(v)
         unit_v(3) = v(3)/norm2(v)
         unit_sph = cart2sph(unit_v)
      end if

      temp1 = (2.0_rp*l + 1.0_rp)/(four_pi)
      temp2 = gamma(1.0_rp*(l - m) + 1.0_rp)/gamma(1.0_rp*(l + m) + 1.0_rp)
      part1 = real((cmone)**(cone*m))*sqrt(temp1*temp2)*cone
      part2 = associated_lp(cos(unit_sph(2)), l, m)*cone
      part3 = zexp(i_unit*1.0_rp*m*unit_sph(3))

      s = part1*part2*part3*unit_sph(1)

   end function complex_spharm

   !> Computes the real spherical harmonics (orbitals)
   !> for a given unit vector in cartesian coordinates (v_x, v_y, v_z)
   !> (if the vector is not unitary, the function forces it)
   !> Implemented by Ivan Miranda on 21.09.2023
   function real_spharm(v, l, m) result(s)
      !
      implicit none
      !
      integer, intent(in) :: l, m
      real(rp), dimension(3), intent(in) :: v ! same as for the complex spherical harmonics
      real(rp) :: s
      !
      ! Local Variables
      complex(rp) :: temp1
      !
      if (m .lt. 0) then
         temp1 = complex_spharm(v, l, abs(m))
         s = sqrt(2.0_rp)*real((cmone)**(cone*m))*aimag(temp1)
      else if (m .eq. 0) then
         s = real(complex_spharm(v, l, 0))
      else
         temp1 = complex_spharm(v, l, m)
         s = sqrt(2.0_rp)*real((cmone)**(cone*m))*real(temp1)
      end if

   end function real_spharm

   !> Anticommutator of two operators A and B
   !> Accepts only real values in.
   !> Implemented by Ivan Miranda on 24.09.2023
   function anticommutator(A, B) result(res)
      !
      implicit none
      !
      real(rp), intent(in) :: A(:, :), B(:, :)
      real(rp) :: res(size(A, 1), size(A, 2))
      !
      res = matmul(A, B) + matmul(B, A)

   end function anticommutator

   !> Analog anticommutator of two operators A and B
   !> but now accepts complex values in.
   !> Implemented by Ivan Miranda on 24.09.2023
   function canticommutator(A, B) result(res)
      !
      implicit none
      !
      complex(rp), intent(in) :: A(:, :), B(:, :)
      complex(rp) :: res(size(A, 1), size(A, 2))
      !
      res = matmul(A, B) + matmul(B, A)

   end function canticommutator

   !> Complex identity matrix
   function ceye(n) result(m)
      integer, intent(in) :: n
      complex(rp), dimension(n, n) :: m
      integer :: i

      m = czero

      do i = 1, n
         m(i, i) = cone
      end do

   end function ceye

   !> Identity matrix
   function eye(n) result(m)
      integer, intent(in) :: n
      real(rp), dimension(n, n) :: m
      integer :: i

      m = 0.0_rp

      do i = 1, n
         m(i, i) = 1.0_rp
      end do
   end function eye

   !> Test if rank 4 matrix is Hermitian
   function is_hermitian_r4(m) result(l)
      complex(rp), dimension(:, :, :, :), intent(in) :: m
      logical :: l
      integer :: ia, io1, io2, is1, is2
      integer, dimension(2, 2), parameter :: iss2is = reshape((/1, 4, 3, 2/), (/2, 2/))

      l = .true.
      do ia = 1, size(m, 1)
         do io1 = 1, size(m, 2)
            do io2 = 1, size(m, 3)
               do is1 = 1, 2
                  do is2 = 1, 2
                     l = l .and. abs(m(ia, io1, io2, iss2is(is1, is2)) &
                                     - conjg(m(ia, io2, io1, iss2is(is2, is1)))) < epsilon
                     if (.not. l) return
                  end do !end is1
               end do !end is2
            end do !end io1
         end do !end io2
      end do !end ia
   end function is_hermitian_r4

   !> Vector cross product
   function cross_product(v1, v2) result(v3)
      real(rp), dimension(:), intent(in) :: v1, v2
      real(rp), dimension(3) :: v3

      v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
      v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
      v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
   end function cross_product

   !> Normalize a vector
   function normalize(v1) result(v2)
      real(rp), dimension(3), intent(in) :: v1
      real(rp), dimension(3) :: v2
      real(rp) :: norm
      norm = norm2(v1)
      v2(:) = v1(:)/norm
   end function normalize

   !> Euler-Rodrigues rotation
   function erodrigues(v1, axis, phi) result(v2)
      real(rp), dimension(3), intent(in) :: v1, axis
      real(rp), intent(in) :: phi
      real(rp), dimension(3) :: v2, norm
      real(rp), dimension(3, 3) :: matrix
      real(rp) :: a, b, c, d

      norm = normalize(axis)
      a = cos(phi/2)
      b = norm(1)*sin(phi/2)
      c = norm(2)*sin(phi/2)
      d = norm(3)*sin(phi/2)
      matrix(1, 1) = a**2 + b**2 - c**2 - d**2
      matrix(1, 2) = 2*(b*c - a*d)
      matrix(1, 3) = 2*(b*d + a*c)
      matrix(2, 1) = 2*(b*c + a*d)
      matrix(2, 2) = a**2 + c**2 - b**2 - d**2
      matrix(2, 3) = 2*(c*d - a*b)
      matrix(3, 1) = 2*(b*d - a*c)
      matrix(3, 2) = 2*(c*d + a*b)
      matrix(3, 3) = a**2 + d**2 - b**2 - c**2
      v2(:) = matmul(matrix, v1)

   end function erodrigues

   !> Matrix determinant
   function determinant(m) result(d)
      real(rp), dimension(3, 3), intent(in) :: m
      real(rp) :: d
      real(rp), dimension(3) :: v1, v2, v3

      v1(:) = m(1, :)
      v2(:) = m(2, :)
      v3(:) = m(3, :)

      d = dot_product(v1, cross_product(v2, v3))

   end function determinant

   !> Distance between two positions
   function distance(v1, v2) result(d)
      real(rp), dimension(3), intent(in) :: v1, v2
      real(rp), dimension(3) :: v3
      real(rp) :: d

      v3(1) = v1(1) - v2(1)
      v3(2) = v1(2) - v2(2)
      v3(3) = v1(3) - v2(3)

      d = sqrt(v3(1)**2 + v3(2)**2 + v3(3)**2)

   end function distance

   !> Angle between two vectors, in radians
   function angle(v1, v2) result(d)
      real(rp), dimension(3), intent(in) :: v1, v2
      real(rp) :: d

      d = acos(dot_product(v1, v2)/(norm2(v1)*norm2(v2)))

   end function angle

   !> Vector position between two atoms
   function pos(v1, v2) result(v3)
      real(rp), dimension(3), intent(in) :: v1, v2
      real(rp), dimension(3) :: v3

      v3(:) = v1(:) - v2(:)

   end function pos

   !> Matrix inverse
   function inverse(m1) result(m2)
      real(rp), dimension(:, :), intent(in) :: m1
      real(rp), dimension(size(m1, 1), size(m1, 2)) :: m2
      integer, dimension(min(size(m1, 1), size(m1, 2))) :: ipiv
#if !defined(LAPACK95_FOUND)
      integer :: m
      integer :: lda
      real(rp), dimension(max(1, size(m1, 1))) :: work
      integer :: lwork
      integer :: info
#endif
#if defined(LAPACK95_FOUND)
      m2 = m1
      call getrf(m2, ipiv)
      call getri(m2, ipiv)
#else
      m = size(m1, 1)
      lda = max(1, m)
      lwork = max(1, m)
      m2 = m1
      call dgetrf(m, m, m2, lda, ipiv, info)
      call dgetri(m, m2, lda, ipiv, work, lwork, info)
#endif
   end function inverse

   !> Calculates numerically the second derivative of a given function
   !> by the central difference method.
   !> Implemented by Ivan Miranda on 14.11.2023
   function second_derivative(x, h) result(second_deriv)
      !
      implicit none
      !
      real(rp), intent(in) :: x(:) ! Array of elements representing the function
      real(rp), intent(in) :: h ! Step size
      real(rp) :: second_deriv(size(x) - 2)
      !
      ! Local Variables
      integer :: i
      !
      do i = 2, size(x) - 1
         second_deriv(i - 1) = (x(i + 1) - 2.0_rp*x(i) + x(i - 1))/(h**2)
      end do

   end function second_derivative

   !> Matrix inverse of a 3x3 matrix
   function inverse_3x3(m1) result(m2)
      real(rp), dimension(3, 3), intent(in) :: m1
      real(rp), dimension(3, 3) :: m2
      real(rp) :: det

      det = determinant(m1)
      if (det >= epsilon) then

         m2(1, 1) = (m1(2, 2)*m1(3, 3) - m1(2, 3)*m1(3, 2))/det; 
         m2(1, 2) = -(m1(1, 2)*m1(3, 3) - m1(1, 3)*m1(3, 2))/det; 
         m2(1, 3) = (m1(1, 2)*m1(2, 3) - m1(1, 3)*m1(2, 2))/det; 
         m2(2, 1) = -(m1(2, 1)*m1(3, 3) - m1(2, 3)*m1(3, 1))/det; 
         m2(2, 2) = (m1(1, 1)*m1(3, 3) - m1(1, 3)*m1(3, 1))/det; 
         m2(2, 3) = -(m1(1, 1)*m1(2, 3) - m1(1, 3)*m1(2, 1))/det; 
         m2(3, 1) = (m1(2, 1)*m1(3, 2) - m1(2, 2)*m1(3, 1))/det; 
         m2(3, 2) = -(m1(1, 1)*m1(3, 2) - m1(1, 2)*m1(3, 1))/det; 
         m2(3, 3) = (m1(1, 1)*m1(2, 2) - m1(1, 2)*m1(2, 1))/det; 
      else
         write (*, *) 'math%inverse_3x3(): Warning - determinant is nearely zero'
         stop
      end if

   end function inverse_3x3

   !> Unique values c in an integer array a, index ia is such that c = a(ia) and
   !> index ic is such that a = c(ic)
   subroutine unique_int(a, c, ia, ic)
      integer, dimension(:), intent(in)  :: a
      integer, dimension(:), allocatable, intent(out) :: c, ia
      integer, dimension(size(a)), intent(out) :: ic
      integer :: sa, sc, ia1, ia2, ic1
      integer, dimension(size(a)) :: index_a
      logical, dimension(size(a)) :: mask

      sa = size(a)
      index_a = (/(ia1, ia1=1, sa)/)
      mask = .false.
      do ia1 = 1, sa - 1
         if (.not. mask(ia1)) then
            do ia2 = ia1 + 1, sa
               mask(ia2) = a(ia1) == a(ia2)
               if (mask(ia2)) index_a(ia2) = index_a(ia1)
            end do
         end if
      end do
      sc = sa - count(mask)
      allocate (c(sc), ia(sc))
      c = pack(a,.not. mask)
      ia = pack(index_a,.not. mask)
      do ia1 = 1, sa
         do ic1 = 1, sc
            if (index_a(ia1) == ia(ic1)) then
               ic(ia1) = ic1
               exit
            end if
         end do
      end do
   end subroutine unique_int

   ! Unique values c in a real array a down to precision epsilon, index ia is
   ! such that c = a(ia) andindex ic is such that a = c(ic)
   !subroutine unique_real(a, c, ia, ic)

   !> Returns a spherical vector \f$ v_{\mathrm{sph}} = (v_r, v_{\theta}, v_{\phi})
   !> \f$ from a cartesian vector \f$ v_{\mathrm{cart}} = (v_x, v_y, v_z) \f$ where
   !> \f$ 0 \leq v_{\theta} \leq \pi \f$ and
   !> \f$ -\pi \leq v_{\phi} < \pi \f$
   function cart2sph(v_cart) result(v_sph)
      real(rp), dimension(3), intent(in) :: v_cart
      real(rp), dimension(3) :: v_sph

      v_sph(1) = norm2(v_cart)

      if (v_sph(1) < epsilon) then
         v_sph(2) = 0.0_rp
      else
         v_sph(2) = acos(v_cart(3)/v_sph(1))
      end if

      if (sqrt(v_cart(1)*v_cart(1) + v_cart(2)*v_cart(2)) < epsilon) then
         v_sph(3) = 0.0_rp
      else
         if (abs(v_cart(2)) > epsilon) then
            if (v_cart(2) > 0.0_rp) then
               v_sph(3) = acos(v_cart(1)/sqrt(v_cart(1)*v_cart(1) + v_cart(2)*v_cart(2)))
            else
               v_sph(3) = -acos(v_cart(1)/sqrt(v_cart(1)*v_cart(1) + v_cart(2)*v_cart(2)))
            end if
         elseif (abs(v_cart(2)) < epsilon) then
            if (v_cart(1) > 0.0_rp) then
               v_sph(3) = 0.0_rp
            else
               v_sph(3) = -pi
            end if
         end if
      end if
      ! write(6, *) "v_cart=", v_cart(:)
      ! write(6, *) "v_sph=", v_sph(:)
   end function cart2sph

   !> Returns a cartesian vector \f$ v_{\mathrm{cart}} = (v_x, v_y, v_z) \f$ from a
   !> spherical vector \f$ v_{\mathrm{sph}} = (v_r, v_{\theta}, v_{\phi}) \f$
   function sph2cart(v_sph) result(v_cart)
      real(rp), dimension(3), intent(in) :: v_sph
      real(rp), dimension(3) :: v_cart

      v_cart(1) = v_sph(1)*sin(v_sph(2))*cos(v_sph(3))
      v_cart(2) = v_sph(1)*sin(v_sph(2))*sin(v_sph(3))
      v_cart(3) = v_sph(1)*cos(v_sph(2))
   end function sph2cart

   subroutine nm2rho(n, m_cart, rho)
      real(rp), intent(in) :: n
      real(rp), dimension(3), intent(in) :: m_cart
      complex(rp), dimension(2, 2), intent(out) :: rho

      rho(1, 1) = (n + m_cart(3))/2
      rho(2, 2) = (n - m_cart(3))/2
      rho(1, 2) = (m_cart(1) - i_unit*m_cart(2))/2
      rho(2, 1) = (m_cart(1) + i_unit*m_cart(2))/2
   end subroutine nm2rho

   subroutine rho2nm(rho, n, m_cart)
      complex(rp), dimension(2, 2), intent(in) :: rho
      real(rp), intent(out) :: n
      real(rp), dimension(3), intent(out) :: m_cart
      integer :: ispin, jspin

      do ispin = 1, 2
         do jspin = 1, 2
            if ((abs(rho(ispin, jspin) - conjg(rho(jspin, ispin)))) > epsilon) then
               write (*, *) 'math%rho2nm(): Warning - non Hermitian density'
            end if
         end do
      end do

      n = real(rho(1, 1) + rho(2, 2))
      m_cart(1) = real(rho(1, 2) + rho(2, 1))
      !m_cart(2) = i_unit*(rho(1, 2)-rho(2, 1))
      m_cart(2) = aimag(rho(2, 1) - rho(1, 2))
      m_cart(3) = real(rho(1, 1) - rho(2, 2))
   end subroutine rho2nm

   !> Compute Fermi(E, beta) occupation function
   function fermi_function(E, beta) result(f)
      real(rp), intent(in) :: E, beta
      real(rp) :: f
      real(rp), parameter :: maxarg = 100._rp
      real(rp) :: x

      x = beta*E

      if (x < -maxarg) then
         f = 1.0_rp
      elseif (x > maxarg) then
         f = 0.0_rp
      else
         f = 1.0_rp/(1.0_rp + exp(x))
      end if
   end function fermi_function

   function fermifun(e, ef, kbt)
      implicit none
      real(rp), intent(in) :: e, ef, kbt
      real(rp) :: fermifun

      fermifun = 1.0d0/(exp((e - ef)/kbt) + 1.0d0)
   end function fermifun

   function dfermifun(e, ef, kbt)
      implicit none
      real(rp), intent(in) :: e, ef, kbt
      real(rp) :: dfermifun
      if ((e - ef)/kbt > 100.0d0) then
         dfermifun = 0.0d0
      else
         dfermifun = -1.0d0/kbt*exp((e - ef)/kbt)/(exp((e - ef)/kbt) + 1.0d0)**2
      end if
   end function dfermifun

   !> Compute dFermi(E, beta)/dE
   function fermi_function_derivative(E, beta) result(f)
      real(rp), intent(in) :: E, beta
      real(rp) :: f
      real(rp), parameter :: maxarg = 100._rp
      real(rp) :: x

      x = beta*E

      if (x < -maxarg) then
         f = 0.0_rp
      elseif (x > maxarg) then
         f = 0.0_rp
      else
         f = -beta*exp(x)/(1.0_rp + exp(x))**2
      end if
   end function fermi_function_derivative

   !> fill an array arr and sort it to get the indx of sorted values
   subroutine indexx(n, arr, indx)
      ! INPUT
      integer, intent(in) :: n
      real(rp), dimension(0:n), intent(in) :: arr
      ! OUTPUT
      integer, dimension(0:n), intent(out) :: indx
      ! LOCAL
      integer, parameter :: m = 7, nstack = 50
      integer, dimension(:), allocatable :: istack
      integer :: i, indxt, ir, itemp, j, jstack, k, l
      real(rp) :: a
      if (.not. allocated(istack)) allocate (istack(nstack))
      indx = (/(j, j=0, n)/)
      jstack = 0
      l = 1
      ir = n
      do
         if (ir - l < m) then
            do j = l + 1, ir
               indxt = indx(j)
               a = arr(indxt)
               do i = j - 1, 1, -1
                  if (arr(indx(i)) <= a) exit
                  indx(i + 1) = indx(i)
               end do
               if (i == 1 .and. arr(indx(i)) > a) i = 0
               indx(i + 1) = indxt
            end do
            if (jstack == 0) return
            ir = istack(jstack)
            l = istack(jstack - 1)
            jstack = jstack - 2
         else
            k = (l + ir)/2
            itemp = indx(k)
            indx(k) = indx(l + 1)
            indx(l + 1) = itemp
            if (arr(indx(l + 1)) > arr(indx(ir))) then
               itemp = indx(l + 1)
               indx(l + 1) = indx(ir)
               indx(ir) = itemp
            end if
            if (arr(indx(l)) > arr(indx(ir))) then
               itemp = indx(l)
               indx(l) = indx(ir)
               indx(ir) = itemp
            end if
            if (arr(indx(l + 1)) > arr(indx(l))) then
               itemp = indx(l + 1)
               indx(l + 1) = indx(l)
               indx(l) = itemp
            end if
            i = l + 1
            j = ir
            indxt = indx(l)
            a = arr(indxt)
            do
               i = i + 1
               if (arr(indx(i)) < a) cycle
               do
                  j = j - 1
                  if (arr(indx(j)) <= a) exit
               end do
               if (j < i) exit
               itemp = indx(i)
               indx(i) = indx(j)
               indx(j) = itemp
            end do
            indx(l) = indx(j)
            indx(j) = indxt
            jstack = jstack + 2
            if (jstack > nstack) then
               write (*, *) "nstack too small in indexx"
               stop
            end if
            if (ir - i + 1 >= j - l) then
               istack(jstack) = ir
               istack(jstack - 1) = i
               ir = j - 1
            else
               istack(jstack) = j - 1
               istack(jstack - 1) = l
               l = i
            end if
         end if
      end do
      if (allocated(istack)) deallocate (istack)
   end subroutine indexx

   !> fill an array arr and sort it to get the indx of sorted values
   recursive subroutine QSort(na, a, indx)

      ! DUMMY ARGUMENTS
      integer, intent(in) :: nA
      type(group), dimension(nA), intent(in) :: A
      integer, dimension(nA), intent(out) :: indx

      ! LOCAL VARIABLES
      integer :: left, right
      real(rp) :: random
      real(rp) :: pivot
      type(group) :: temp
      type(group), dimension(nA) :: B
      integer :: i, marker

      forall (i=1:nA) B(i) = A(i)

      if (nA > 1) then
         call random_number(random)
         pivot = B(int(random*real(nA - 1)) + 1)%value   ! random pivot (not best performance, but avoids worst-case)
         left = 0
         right = nA + 1
         do while (left < right)
            right = right - 1
            do while (B(right)%value > pivot)
               right = right - 1
            end do
            left = left + 1
            do while (B(left)%value < pivot)
               left = left + 1
            end do
            if (left < right) then
               temp = B(left)
               B(left) = B(right)
               B(right) = temp
            end if
         end do

         if (left == right) then
            marker = left + 1
         else
            marker = left
         end if

         call QSort(marker - 1, B(:marker - 1), indx)
         call QSort(nA - marker + 1, B(marker:), indx)
      end if
      do i = 1, nA
         indx(i) = B(i)%order
      end do

   end subroutine QSort

   real(rp) function theta_function(x, smearing)
      !-----------------------------------------------------------------------
      !
      !     This function computes the approximate theta function for the
      !     given order n, at the point x.
      !
      ! --> 'fd': Fermi-Dirac smearing.
      !       1.0/(1.0+exp(-x))
      !
      ! --> 'g':  Gaussian smearing (n=0).
      !
      ! --> 'mp': Methfessel-Paxton smearing (n=1). See PRB 40, 3616 (1989).
      !
      ! --> 'mv': Marzari-Vanderbilt (cold) smearing. See PRL 82, 3296 (1999).
      !       1/2*erf(x-1/sqrt(2)) + 1/sqrt(2*pi)*exp(-(x-1/sqrt(2))**2) + 1/2
      !
      real(rp), intent(in) :: x
      character(len=*), intent(in) :: smearing
      ! output: the value of the function
      ! input: the argument of the function
      ! input: the order of the function
      !
      !    the local variables
      !
      real(rp) :: a, hp, arg, hd, xp
      ! the coefficient a_n
      ! the hermitean function
      ! the argument of the exponential
      ! the hermitean function
      ! auxiliary variable (cold smearing)
      integer :: i, ni, ngauss
      ! counter on the n indices
      ! counter on 2n
      real(rp), parameter :: maxarg = 100._rp
      ! maximum value for the argument of the exponential

      ! Fermi-Dirac smearing
      if (smearing == 'fd') then
         if (x < -maxarg) then
            theta_function = 0._rp
         elseif (x > maxarg) then
            theta_function = 1._rp
         else
            theta_function = 1.0_rp/(1.0_rp + exp(-x))
         end if
         return
      end if
      ! Cold smearing
      if (smearing == 'mv') then
         xp = x - one_over_sqrt_two
         arg = min(maxarg, xp**2)
         theta_function = 0.5_rp*erf_qe(xp) + one_over_sqrt_two_pi*exp(- &
                                                                       arg) + 0.5_rp
         return
      end if
      if (smearing == 'mp') ngauss = 1
      if (smearing == 'g') ngauss = 0

      ! Methfessel-Paxton
      theta_function = 0.5_rp*erfc_qe(-x)

      if (ngauss == 0) return
      hd = 0._rp
      arg = min(maxarg, x**2)
      hp = exp(-arg)
      ni = 0
      a = one_over_sqrt_pi
      do i = 1, ngauss
         hd = 2.0_rp*x*hp - 2.0_rp*dble(ni)*hd
         ni = ni + 1
         a = -a/(dble(i)*4.0_rp)
         theta_function = theta_function - a*hd
         hp = 2.0_rp*x*hd - 2.0_rp*dble(ni)*hp
         ni = ni + 1
      end do
      return
   end function theta_function

   !***********************************************************************|
   real(rp) function delta_function(x, smearing)
      !
      !     The derivative of wgauss: an approximation to the delta function
      !
      ! --> 'fd': derivative of the Fermi-Dirac smearing.
      !       0.5/(1.0+cosh(x))
      !
      ! --> 'g':  derivative of the Gaussian smearing (n=0).
      !
      ! --> 'mp': derivative of the Methfessel-Paxton smearing.
      !
      ! --> 'mv': derivative of the Marzari-Vanderbilt (cold) smearing.
      !       1/sqrt(pi)*exp(-(x-1/sqrt(2))**2)*(2-sqrt(2)*x)
      !
      real(rp), intent(in) :: x
      character(len=*), intent(in) :: smearing
      ! output: the value of the function
      ! input: the point where to compute the function

      ! input: the order of the smearing function
      !
      !    here the local variables
      !
      real(rp) :: a, arg, hp, hd
      ! the coefficients a_n
      ! the argument of the exponential
      ! the hermite function
      ! the hermite function

      integer :: i, ni, ngauss
      ! counter on n values
      ! counter on 2n values

      ! Fermi-Dirac smearing
      if (smearing == 'fd') then
         if (abs(x) <= 36.0_rp) then
            delta_function = 1.0_rp/(2.0_rp + exp(-x) + exp(+x))
            ! in order to avoid problems for large values of x in the e
         else
            delta_function = 0.0_rp
         end if
         return
      end if
      ! cold smearing  (Marzari-Vanderbilt)
      if (smearing == 'mv') then
         arg = min(200.0_rp, (x - one_over_sqrt_two)**2)
         delta_function = one_over_sqrt_pi*exp(-arg)*(2.0_rp - sqrt_two*x)
         return
      end if
      ! Methfessel-Paxton
      if (smearing == 'mp') ngauss = 1
      if (smearing == 'g') ngauss = 0

      arg = min(200.0_rp, x**2)
      delta_function = exp(-arg)*one_over_sqrt_pi
      if (ngauss == 0) return
      hd = 0.0_rp
      hp = exp(-arg)
      ni = 0
      a = one_over_sqrt_pi
      do i = 1, ngauss
         hd = 2.0_rp*x*hp - 2.0_rp*dble(ni)*hd
         ni = ni + 1
         a = -a/(dble(i)*4.0_rp)
         hp = 2.0_rp*x*hd - 2.0_rp*dble(ni)*hp
         ni = ni + 1
         delta_function = delta_function + a*hp
      end do
      return
   end function delta_function

   real(rp) function integrated_delta_function(x, smearing)
      !-----------------------------------------------------------------------
      !
      !    w1gauss(x, n) = \int_{-\infty}^x   y delta(y) dy
      !    where delta(x) is the current approximation for the delta function,
      !    as obtained from w0gauss(x, n)
      !
      ! --> 'fd.: Fermi-Dirac smearing (n=-99). In this case w1gauss corresponds
      !     to the negative of the electronic entropy.
      !
      ! --> 'mp': Methfessel-Paxton smearing (n>=0).
      !
      ! --> 'mv': Marzari-Vanderbilt (cold) smearing (n=-1).
      !       1/sqrt(2*pi)*(x-1/sqrt(2))*exp(-(x-1/sqrt(2))**2)
      !
      real(rp), intent(in) :: x
      character(len=*), intent(in) :: smearing
      ! output: the value of the function
      ! input: the point where to compute the function

      ! input: the order of the smearing function
      !
      !    here the local variables
      !

      real(rp) :: a, hp, arg, hpm1, hd, f, onemf, xp
      ! the coefficients a_n
      ! the hermite function
      ! the argument of the exponential
      ! the hermite function
      ! the hermite function
      ! Fermi-Dirac occupation number
      ! 1 - f
      ! auxiliary variable (cold smearing)

      integer :: i, ni, ngauss
      ! counter on n values
      ! counter on 2n values

      if (smearing == 'fd') then
         if (abs(x) <= 36.0_rp) then
            f = 1.0_rp/(1.0_rp + exp(-x))
            onemf = 1.0_rp - f
            integrated_delta_function = f*log(f) + onemf*log(onemf)
            ! in order to avoid problems for large values of x
         else
            ! neglect w1gauss when abs(w1gauss) < 1.0d-14
            integrated_delta_function = 0.0_rp
         end if
         return

      end if
      ! Cold smearing
      if (smearing == 'mv') then
         xp = x - one_over_sqrt_two
         arg = min(200._rp, xp**2)
         integrated_delta_function = one_over_sqrt_two_pi*xp*exp(-arg)
         return

      end if
      ! Methfessel-Paxton
      if (smearing == 'mp') ngauss = 1
      if (smearing == 'g') ngauss = 0

      arg = min(200.0_rp, x**2)
      integrated_delta_function = -0.5_rp*one_over_sqrt_pi*exp(-arg)
      if (ngauss == 0) return
      hd = 0.0_rp
      hp = exp(-arg)
      ni = 0
      a = one_over_sqrt_pi
      do i = 1, ngauss
         hd = 2.0_rp*x*hp - 2.0_rp*dble(ni)*hd
         ni = ni + 1
         hpm1 = hp
         hp = 2.0_rp*x*hd - 2.0_rp*dble(ni)*hp
         ni = ni + 1
         a = -a/(dble(i)*4.0_rp)
         integrated_delta_function = integrated_delta_function - a*(0.5_rp*hp + dble(ni)*hpm1)
      end do
      return
   end function integrated_delta_function

   !***********************************************************************|
   !
   ! Copyright (C) 2002-2009 Quantum ESPRESSO group
   ! This file is distributed under the terms of the
   ! GNU General Public License. See the file `License'
   ! in the root directory of the present distribution,
   ! or http://www.gnu.org/copyleft/gpl.txt .
   !
   !---------------------------------------------------------------------
   real(rp) function erf_qe(x)
      !---------------------------------------------------------------------
      !
      !     Error function - computed from the rational approximations of
      !     W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
      !
      !     for abs(x) le 0.47 erf is calculated directly
      !     for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
      !
      real(rp), intent(in) :: x
      real(rp) :: x2, p1(4), q1(4)
      data p1/2.426679552305318E2_rp, 2.197926161829415E1_rp, &
         6.996383488619136_rp, -3.560984370181538E-2_rp/
      data q1/2.150588758698612E2_rp, 9.116490540451490E1_rp, &
         1.508279763040779E1_rp, 1.000000000000000_rp/
      !
      if (abs(x) > 6.0_rp) then
         !
         !  erf(6)=1-10^(-17) cannot be distinguished from 1
         !
         erf_qe = sign(1.0_rp, x)
      else
         if (abs(x) <= 0.47_rp) then
            x2 = x**2
            erf_qe = x*(p1(1) + x2*(p1(2) + x2*(p1(3) + x2*p1(4)))) &
                     /(q1(1) + x2*(q1(2) + x2*(q1(3) + x2*q1(4))))
         else
            erf_qe = 1.0_rp - erfc_qe(x)
         end if
      end if
      !
      return
   end function erf_qe

   real(rp) function erfc_qe(x)
      !---------------------------------------------------------------------
      !
      !     Complementary error function
      !     erfc(x) = 1-erf(x)  - See comments in erf
      !
      real(rp), intent(in) :: x
      real(rp) :: ax, x2, xm2, p2(8), q2(8), p3(5), q3(5), pim1
      data p2/3.004592610201616E2_rp, 4.519189537118719E2_rp, &
         3.393208167343437E2_rp, 1.529892850469404E2_rp, &
         4.316222722205674E1_rp, 7.211758250883094_rp, &
         5.641955174789740E-1_rp, -1.368648573827167E-7_rp/
      data q2/3.004592609569833E2_rp, 7.909509253278980E2_rp, &
         9.313540948506096E2_rp, 6.389802644656312E2_rp, &
         2.775854447439876E2_rp, 7.700015293522947E1_rp, &
         1.278272731962942E1_rp, 1.000000000000000_rp/
      data p3/-2.996107077035422E-3_rp, -4.947309106232507E-2_rp, &
         -2.269565935396869E-1_rp, -2.786613086096478E-1_rp, &
         -2.231924597341847E-2_rp/
      data q3/1.062092305284679E-2_rp, 1.913089261078298E-1_rp, &
         1.051675107067932_rp, 1.987332018171353_rp, &
         1.000000000000000_rp/

      data pim1/0.56418958354775629_rp/
      !        ( pim1= sqrt(1/pi) )
      ax = abs(x)
      if (ax > 26.0_rp) then
         !
         !  erfc(26.0)=10^(-296); erfc( 9.0)=10^(-37);
         !
         erfc_qe = 0.0_rp
      elseif (ax > 4.0_rp) then
         x2 = x**2
         xm2 = (1.0_rp/ax)**2
         erfc_qe = (1.0_rp/ax)*exp(-x2)*(pim1 + xm2*(p3(1) &
                                                     + xm2*(p3(2) + xm2*(p3(3) + xm2*(p3(4) + xm2*p3(5) &
                                                                                      ))))/(q3(1) + xm2*(q3(2) + xm2*(q3(3) + xm2* &
                                                                                                                      (q3(4) + xm2*q3(5))))))
      elseif (ax > 0.47_rp) then
         x2 = x**2
         erfc_qe = exp(-x2)*(p2(1) + ax*(p2(2) + ax*(p2(3) &
                                                     + ax*(p2(4) + ax*(p2(5) + ax*(p2(6) + ax*(p2(7) &
                                                                                               + ax*p2(8))))))))/(q2(1) + ax*(q2(2) + ax* &
                                                                                                                              (q2(3) + ax*(q2(4) + ax*(q2(5) + ax*(q2(6) + ax* &
                                                                                                                                                                   (q2(7) + ax*q2(8))))))))
      else
         erfc_qe = 1.0_rp - erf_qe(ax)
      end if
      !
      ! erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)
      !
      if (x < 0.0_rp) erfc_qe = 2.0_rp - erfc_qe
      !
      return
   end function erfc_qe

   ! Transform the spd matrix from spherical harmonics/cartesian to cartesian/spherical harmonics
   subroutine hcpx(ham, transformation)
      ! Input
      character(len=*), intent(in) :: transformation
      complex(rp), dimension(9, 9), intent(inout) :: ham
      ! local variables
      integer :: i, ii, j, jj, k, l
      complex(rp) :: const, del, dum, pim, cone
      complex(rp), dimension(9, 9) :: hesf, hcart, v, vc, htmp
      cone = (1.0d0, 0.0d0)
      const = CMPLX(1.0d0/sqrt(2.0d0), KIND=Kind(.0d0))
      do ii = 1, 9
         do jj = 1, 9
            v(ii, jj) = (0.0d0, 0.0d0)
            vc(ii, jj) = (0.0d0, 0.0d0)
         end do
      end do
      !---base Y(lm) na ordem (00)(1-1)(10)(11)(2-2)(2-1)(20)(21)(22)---
      ! s
      v(1, 1) = cone
      vc(1, 1) = cone
      ! p
      v(2, 4) = -const
      vc(4, 2) = -const
      v(2, 2) = const
      vc(2, 2) = const
      v(3, 4) = i_unit*const
      vc(4, 3) = -i_unit*const
      v(3, 2) = i_unit*const
      vc(2, 3) = -i_unit*const
      v(4, 3) = cone
      vc(3, 4) = cone
      ! d
      v(5, 5) = i_unit*const
      v(5, 9) = -i_unit*const
      v(6, 6) = i_unit*const
      v(6, 8) = i_unit*const
      v(7, 6) = const
      v(7, 8) = -const
      v(8, 5) = const
      v(8, 9) = const
      v(9, 7) = cone
      vc(5, 5) = -i_unit*const
      vc(9, 5) = i_unit*const
      vc(6, 6) = -i_unit*const
      vc(8, 6) = -i_unit*const
      vc(6, 7) = const
      vc(8, 7) = -const
      vc(5, 8) = const
      vc(9, 8) = const
      vc(7, 9) = cone

      hesf = 0.0d0
      hcart = 0.0d0

      select case (transformation)

      case ('cart2sph')
         htmp = matmul(ham, v)
         hesf = matmul(vc, htmp)

         ham(:, :) = hesf(:, :)
      case ('sph2cart')
         htmp = matmul(ham, vc)
         hcart = matmul(v, htmp)

         ham(:, :) = hcart(:, :)
      end select
   end subroutine hcpx

   !> Integration of a given function Y by Simpson, modified to integrate
   !energy^(n)*f
   subroutine simpson_m(AINT, H, EF, NPTS, Y, EA, nexp, Ene)
      ! Input
      integer, intent(in) :: NPTS, nexp
      real(rp), intent(in) :: EA, EF, H
      real(rp), dimension(NPTS + 2), intent(in) :: Y
      real(rp), dimension(NPTS + 2), intent(in) :: Ene
      ! Output
      real(rp), intent(out) :: AINT
      ! Local variables
      integer :: I

      AINT = 0.d0
      do I = 2, NPTS - 1, 2
         AINT = AINT + Y(I - 1)*Ene(I - 1)**nexp + 4.d0*Y(I)*Ene(I)**nexp + Y(I + 1)*Ene(I + 1)**nexp
      end do
      AINT = H*AINT/3.d0
      if ((EA /= EF)) then
         AINT = AINT + (EF - EA)*(Y(NPTS)*Ene(NPTS)**nexp + 4.d0*Y(NPTS + 1)*Ene(NPTS + 1)**nexp + Y(NPTS + 2)*Ene(NPTS + 2)**nexp)/6.d0
      end if
   end subroutine simpson_m

   subroutine simpson_f(AINT, Ene, EF, NPTS, Y, fermi, dfermi, T)
      implicit none
      ! Input
      integer, intent(in) :: NPTS
      real(rp), intent(in) :: EF, T
      logical, intent(in) :: fermi, dfermi
      real(rp), intent(out) :: AINT
      real(rp), dimension(NPTS + 10), intent(in) :: Y, Ene
      ! Local variables
      integer :: I
      real(rp) :: H
      real(rp), parameter :: kB = 0.633362019d-5
      real(rp) :: kBT

      AINT = 0.d0
      H = Ene(2) - Ene(1)
      kBT = kB*T + 1.0d-15
      !
      if (fermi) then
         do I = 2, NPTS + 9, 2
            AINT = AINT + Y(I - 1)*fermifun(Ene(i - 1), Ef, kBT) + 4.d0*Y(I)*fermifun(Ene(i), Ef, kBT) + Y(I + 1)*fermifun(Ene(i + 1), Ef, kBT)
         end do
      else if (dfermi) then
         do I = 2, NPTS + 9, 2
            AINT = AINT + Y(I - 1)*dfermifun(Ene(i - 1), Ef, kBT) + 4.d0*Y(I)*dfermifun(Ene(i), Ef, kBT) + Y(I + 1)*dfermifun(Ene(i + 1), Ef, kBT)
         end do
      else
         do I = 2, NPTS - 1, 2
            AINT = AINT + Y(I - 1) + 4.d0*Y(I) + Y(I + 1)
         end do
      end if
      AINT = H*AINT/3.d0
   end subroutine simpson_f

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> The optimal Jackson (gaussian) kernel of order `N` for the Chebyshev
   !moments Eq. (71) in the KPM review
   ![https://doi.org/10.1103/RevModPhys.78.275]
   !---------------------------------------------------------------------------
   subroutine jackson_kernel(n, k)
      integer, intent(in) :: n
      real(rp), dimension(n), intent(inout) :: k
      ! Local variables
      integer ::  ll, bign
      real(rp) :: theta

      bign = size(k(:))

      do ll = 1, bign
         theta = pi*(real(ll) - 1)/(real(bign) + 1)
         k(ll) = (real(bign) - (real(ll) - 1) + 1)*cos(theta) + sin(theta)/tan(pi/(real(bign) + 1))
         k(ll) = k(ll)/(real(bign) + 1)
      end do
   end subroutine jackson_kernel

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> The Lorentz kernel, the best for the Green Function
   ! Eq.(79) of [https://doi.org/10.1103/RevModPhys.78.275]
   !---------------------------------------------------------------------------
   subroutine lorentz_kernel(n, k, lambda)
      integer, intent(in) :: n
      real(rp), dimension(2*n + 2), intent(inout) :: k
      real(rp), intent(in) :: lambda
      ! Local variables
      integer ::  ll, bign
      real(rp) :: theta

      bign = size(k(:))

      do ll = 1, bign
         theta = lambda*(1 - (real(ll) - 1)/real(bign))
         k(ll) = sinh(theta)/sinh(lambda)
      end do
   end subroutine lorentz_kernel

   subroutine t_polynomial(m, n, x, v)
      !*****************************************************************************80
      !
  !! T_POLYNOMIAL evaluates Chebyshev polynomials T(n, x).
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    01 March 2001
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Reference:
      !
      !    Milton Abramowitz, Irene Stegun,
      !    Handbook of Mathematical Functions,
      !    National Bureau of Standards, 1964,
      !    ISBN: 0-486-61272-4,
      !    LC: QA47.A34.
      !
      !  Parameters:
      !
      !    Input, integer M, the number of evaluation points.
      !
      !    Input, integer N, the highest polynomial to compute.
      !
      !    Input, real ( kind = rk ) X(1:M), the evaluation points.
      !
      !    Output, real ( kind = rk ) V(1:M, 0:N), the values of the polynomials.
      !*******************************************************************************

      implicit none

      integer, parameter :: rk = kind(1.0D+00)

      integer m
      integer n

      integer j
      real(kind=rk) v(1:m, 0:n)
      real(kind=rk) x(1:m)

      if (n < 0) then
         return
      end if

      v(1:m, 0) = 1.0D+00

      if (n < 1) then
         return
      end if

      v(1:m, 1) = x(1:m)

      do j = 2, n
         v(1:m, j) = 2.0D+00*x(1:m)*v(1:m, j - 1) - v(1:m, j - 2)
      end do

      return
   end subroutine t_polynomial

   subroutine chebyshev_gauss_quadrature(recursion_steps, x_i, w_i)
      ! Chebyshev-Gauss quadrature points and weights.
      ! Jie Shen, Tao Tang, Li-Lian Wang, Spectral methods (2011)
      integer, intent(inout) :: recursion_steps
      real(rp), dimension(recursion_steps), intent(inout) ::  x_i, w_i
      ! Local variables
      integer, dimension(recursion_steps) :: i
      integer :: j

      i(:) = 0
      do j = 2, recursion_steps
         i(j) = i(j - 1) + 1
      end do

      x_i(:) = -cos((2*i(:) + 1)*pi/(2*recursion_steps + 2))
      w_i(:) = (pi/(recursion_steps + 1))*sqrt(1 - x_i(:)**2)
   end subroutine chebyshev_gauss_quadrature

   subroutine gauss_legendre(n, a, b, x, w)
      implicit none
      integer, intent(in) :: n
      real(rp), intent(in) :: a, b
      real(rp), dimension(n), intent(out) :: x, w
      real(rp) :: z, z1, p1, p2, p3, pp
      real(rp), parameter :: eps = 1.0d-14
      integer :: m, i, j

      m = nint((n + 1)/2.0d0)

      do i = 1, m
         z = cos(pi*(i - 0.25d0)/(n + 0.5d0))
         do
            p1 = 1.0d0
            p2 = 0.0d0
            do j = 1, n
               p3 = p2
               p2 = p1
               p1 = ((2.0d0*j - 1.0d0)*z*p2 - (j - 1.0d0)*p3)/j
            end do
            pp = n*(z*p1 - p2)/(z*z - 1.0d0)
            z1 = z
            z = z1 - p1/pp
            if (abs(z - z1) .le. eps) exit
         end do
         x(i) = a + (b - a)*(z + 1)/2
         x(n + 1 - i) = a + (b - a)*(1 - z)/2
         w(i) = (b - a)*2.0d0/((1.0d0 - z*z)*pp*pp)/2
         w(n + 1 - i) = w(i)
      end do
   end subroutine gauss_legendre

   !> Calculates the imaginary part of the trace
   !> of a square matrix of any dimension.
   !> Implemented by Ivan Miranda on 04.10.2023
   function imtrace(mat) result(ires)
      !
      implicit none
      !
      complex(rp), intent(in) :: mat(:, :)
      real(rp) :: ires
      !
      ! Local Variables
      integer :: i

      ires = 0.0_rp
      do i = 1, size(mat, 1)
         ires = ires + aimag(mat(i, i))
      end do

   end function imtrace

   !> Calculated the real part of the trace
   !> of a square matrix of any dimension.
   !> Implemented by Ivan Miranda on 04.10.2023
   function rtrace(mat) result(rres)
      !
      implicit none
      !
      complex(rp), intent(in) :: mat(:, :)
      real(rp) :: rres
      !
      ! Local Variables
      integer :: i

      rres = 0.0_rp
      do i = 1, size(mat, 1)
         rres = rres + real(mat(i, i))
      end do

   end function rtrace

   !> Calculates the multiplication of a matrix
   !> by a scalar number (both complex).
   !> Implemented by Ivan Miranda on 16.10.2023
   function scalar_multiply(A, s) result(B)
      !
      implicit none
      !
      complex(rp), intent(in) :: A(:, :), s
      complex(rp), dimension(size(A, 1), size(A, 2)) :: B
      !
      ! Local Variables
      integer :: i, j

      do i = 1, size(A, 1)
         do j = 1, size(A, 2)
            B(i, j) = A(i, j)*s
         end do
      end do

   end function scalar_multiply

   function imtrace9(mat)
      implicit none
      complex(rp), dimension(9, 9), intent(in) :: mat
      real(rp) :: imtrace9

      imtrace9 = aimag(mat(1, 1) + mat(2, 2) + mat(3, 3) + mat(4, 4) + mat(5, 5) + mat(6, 6) + mat(7, 7) + mat(8, 8) + mat(9, 9))
   end function imtrace9

   function rtrace9(mat)
      implicit none
      complex(rp), dimension(9, 9), intent(in) :: mat
      real(rp) :: rtrace9

      rtrace9 = real(mat(1, 1) + mat(2, 2) + mat(3, 3) + mat(4, 4) + mat(5, 5) + mat(6, 6) + mat(7, 7) + mat(8, 8) + mat(9, 9))
   end function rtrace9

   function trace9(mat)
      implicit none
      complex(rp), dimension(9, 9), intent(in) :: mat
      complex(rp) :: trace9

      trace9 = (mat(1, 1) + mat(2, 2) + mat(3, 3) + mat(4, 4) + mat(5, 5) + mat(6, 6) + mat(7, 7) + mat(8, 8) + mat(9, 9))
   end function trace9

   function cartesian_to_direct(a, crd_car)
      implicit none
      real(rp), intent(in) :: a(3, 3), crd_car(3)
      real(rp) :: cartesian_to_direct(3)
      real(rp) :: a_inv(3, 3), work(3)
      integer :: pivot(3), info

      ! calculate inverse of a using LAPACK
      call DGETRF(3, 3, a, 3, pivot, info)
      if (info /= 0) stop "Error in DGETRF"
      a_inv = a
      call DGETRI(3, a_inv, 3, pivot, work, 3, info)
      if (info /= 0) stop "Error in DGETRI"

      ! transform from cartesian to direct coordinates
      cartesian_to_direct = matmul(a_inv, crd_car)
   end function cartesian_to_direct

   function direct_to_cartesian(a, crd_dir)
      implicit none
      real(rp), intent(in) :: a(3, 3), crd_dir(3)
      real(rp) :: direct_to_cartesian(3)

      ! transform from direct to cartesian coordinates
      direct_to_cartesian = matmul(a, crd_dir)
   end function direct_to_cartesian

   !
   ! Subroutine to rotate the matrix
   !
   subroutine rotmag(MAT, sdim)
      implicit none

      ! Formal Arguments
      integer, intent(in) :: sdim
      complex(rp), dimension(18, 18, sdim), intent(inout) :: MAT

      ! Local Scalars
      integer :: i, j, k, m, ndim
      real(rp) :: alpha, beta, gamma, pi

      ! Local Arrays
      complex(rp), dimension(18, 18) :: MATs, RMAT, tmpMAT

      ! Define constants and initialize variables
      pi = acos(-1.0_rp)
      ndim = 18
      RMAT = 0.0_rp
      alpha = -pi/2.0_rp  ! Rotation angle in radians
      beta = -pi/2.0_rp   ! Rotation angle in radians
      gamma = 0.0_rp        ! Rotation angle in radians

      ! Generate the rotation matrix
      call ROTMAT(RMAT, alpha, beta, gamma)

      ! Rotate each matrix
      do k = 1, sdim
         ! Copy MAT to MATs
         do m = 1, 18
            do i = 1, 18
               MATs(i, m) = MAT(i, m, k)
            end do
         end do

         ! Apply the rotation: h * Ux
         call zgemm("N", "C", ndim, ndim, ndim, cone, MATs, ndim, RMAT, ndim, czero, tmpMAT, ndim)

         ! Apply the rotation: Ux' * (h * Ux)
         call zgemm("N", "N", ndim, ndim, ndim, cone, RMAT, ndim, tmpMAT, ndim, czero, MATs, ndim)

         ! Copy MATs back to MAT
         do m = 1, 18
            do i = 1, 18
               MAT(i, m, k) = MATs(i, m)
            end do
         end do
      end do

      return
   end subroutine rotmag

!
   subroutine rotmag_loc(MATout, MATin, sdim, MOM) !
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Parameters ..
      complex(selected_real_kind(15)), parameter :: cone = (1.0d+0, 0.0d+0)
      complex(selected_real_kind(15)), parameter :: czero = (0.0d+0, 0.0d+0)
      !
      !.. Formal Arguments ..
      integer, intent(in) :: sdim
      complex(selected_real_kind(15)), &
         dimension(18, 18, sdim), intent(out) :: MATout
      complex(selected_real_kind(15)), &
         dimension(18, 18, sdim), intent(out) :: MATin
      real(selected_real_kind(15)), dimension(3), intent(in) :: MOM
      !
      !
      !.. Local Scalars ..
      integer :: I, J, K, M, ndim
      real(selected_real_kind(15)) :: ALFA, BETA, GAMA, n1, x1, y1, x2, z2, th, phi
      !
      !.. Local Arrays ..
      real(selected_real_kind(15)), dimension(3) :: v, sv
      complex(selected_real_kind(15)), dimension(18, 18) :: MATs, RMAT, tmpMAT
      !
    !!.. External Calls ..
      !external ROTMAT, zgemm
      !
      !.. Intrinsic Functions ..
      !intrinsic selected_real_kind
      !
      ! ... Executable Statements ...
      !
      ndim = 18
      RMAT = (0.0d0, 0.0d0)
      !
      v = MOM

      ! Option 1)  Spherical angles
      call car2sph(v, sv)
      alfa = sv(1)
      beta = sv(2)
      gama = 0.0d0

      call ROTMAT(RMAT, alfa, beta, gama)
      !  call ROTMAT_org(RMAT,alfa,beta,gama)

      do J = 1, sdim
         do M = 1, 18
            do K = 1, 18
               MATs(K, M) = MATin(K, M, J)
            end do
         end do
         !  h*Ux
         call zgemm("N", "N", ndim, ndim, ndim, cone, MATs, ndim, RMAT, ndim, czero, tmpMAT, &
                    ndim)
         !  Ux'*(h*Ux)
         call zgemm("C", "N", ndim, ndim, ndim, cone, RMAT, ndim, tmpMAT, ndim, czero, MATs, &
                    ndim)
         !
         do M = 1, 18
            do K = 1, 18
               MATout(K, M, J) = MATs(K, M)
            end do
         end do
      end do
      return
      !
      ! ... Format Declarations ...
      !
10000 format(3x, a8, f10.4, a8, f10.4, a8, f10.4)
   end subroutine rotmag_loc

   subroutine ROTMAT(MAT, A, B, G)
      implicit none

      ! Formal Arguments
      complex(rp), dimension(18, 18), intent(inout) :: MAT
      real(rp), intent(in) :: A, B, G

      ! Local Scalars
      integer :: J, M, Mprime, S
      real(rp) :: BETA
      complex(rp) :: ALFA, GAMA, IM

      ! Local Arrays
      complex(rp), dimension(2, 2) :: SM
      complex(rp), dimension(9, 9) :: MAT9

      ! Initialization
      SM = 0.0_rp
      MAT = 0.0_rp
      MAT9 = 0.0_rp
      IM = CMPLX(0.0_rp, 1.0_rp)
      ALFA = CMPLX(A, KIND=KIND(0.0_rp))
      BETA = B
      GAMA = CMPLX(G, KIND=KIND(0.0_rp))

      ! Generate the submatrix SM for the rotation
      SM(1, 1) = DSs(0.5_rp, 0.5_rp, 0.5_rp, BETA)*exp(-IM*(0.5_rp*ALFA + 0.5_rp*GAMA))
      SM(1, 2) = DSs(0.5_rp, 0.5_rp, -0.5_rp, BETA)*exp(-IM*(0.5_rp*ALFA - 0.5_rp*GAMA))
      SM(2, 1) = DSs(0.5_rp, -0.5_rp, 0.5_rp, BETA)*exp(-IM*(-0.5_rp*ALFA + 0.5_rp*GAMA))
      SM(2, 2) = DSs(0.5_rp, -0.5_rp, -0.5_rp, BETA)*exp(-IM*(-0.5_rp*ALFA - 0.5_rp*GAMA))

      ! Calculate MAT9 for J=0,1,2 and corresponding M, Mprime
      do J = 0, 2
         S = J*J + 1 + J
         do M = -J, J
            do Mprime = -J, J
               MAT9(S + M, S + Mprime) = DSs(J*1.0d0, M*1.0d0, Mprime*1.0d0, BETA)*exp(-IM*(M*ALFA + Mprime*GAMA))
            end do
         end do
      end do

      ! Populate the full matrix MAT using MAT9 and SM
      do M = 1, 9
         do Mprime = 1, 9
            MAT(Mprime, M) = MAT9(Mprime, M)*SM(1, 1)
            MAT(Mprime, M + 9) = MAT9(Mprime, M)*SM(1, 2)
            MAT(Mprime + 9, M) = MAT9(Mprime, M)*SM(2, 1)
            MAT(Mprime + 9, M + 9) = MAT9(Mprime, M)*SM(2, 2)
         end do
      end do
   end subroutine ROTMAT

   function DSs(J, M, Mp, beta)
      use precision_mod, only: rp
      implicit none
      real(rp), intent(in) :: J, M, Mp, beta
      real(rp) :: DSs, DSt, ds
      integer :: s, smax, smin

      smin = max(0, nint(-Mp - M))
      smax = min(nint(J - Mp), nint(J - M))
      DSs = 0.0_rp

      do s = smin, smax
         ds = s*1.0_rp
         DSt = binom(nint(J + M), nint(J - Mp - ds))*binom(nint(J - M), s)*(-1.0_rp)**nint(J - Mp - ds)
         DSs = DSs + DSt*cos(0.5_rp*beta)**(2*ds + Mp + M)*sin(0.5_rp*beta)**(2*J - 2*ds - Mp - M)
      end do

      DSs = DSs*sqrt(1.0d0*factint(nint(J + Mp))*factint(nint(J - Mp))/(factint(nint(J - M))*factint(nint(J + M))))
   end function DSs

   integer function binom(x, y)
      implicit none
      integer, intent(in) :: x, y

      if (y < 0 .or. y > x) then
         binom = 0
         return
      end if

      binom = factint(x)/(factint(y)*factint(x - y))
   end function binom

   integer function factint(n)
      implicit none
      integer, intent(in) :: n
      integer :: i

      if (n < 0) then
         factint = 0
         return
      end if

      factint = 1
      do i = 1, n
         factint = factint*i
      end do
   end function factint

   subroutine car2sph(C, S)
      ! transforms cartesian (x,y,z) to spherical (Theta,Phi,R) coordinates
      !
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !
      !.. Formal Arguments ..
      real(selected_real_kind(15)) :: X, Y, Z, THETA, PHI, D2, R2
      real(selected_real_kind(15)), dimension(3), intent(in) :: C
      real(selected_real_kind(15)), dimension(3), intent(out) :: S
      !
      !
      ! ... Executable Statements ...
      !
      !
      X = C(1)
      Y = C(2)
      Z = C(3)
      D2 = X*X + Y*Y
      R2 = X*X + Y*Y + Z*Z

      IF (D2 .EQ. 0D0) THEN
         THETA = 0D0
      ELSE
         THETA = ATAN2(Y, X)!-4.0d0*atan(1.0d0)
      END IF

      PHI = ACOS(Z/R2)

      S(1) = THETA
      S(2) = PHI
      S(3) = R2

      return
      !
   end subroutine car2sph

   subroutine sph2car(C, S)
      ! transforms spherical (Theta,Phi,R) to cartesian (x,y,z) coordinates
      !
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !
      !.. Formal Arguments ..
      real(selected_real_kind(15)), dimension(3), intent(out) :: C
      real(selected_real_kind(15)), dimension(3), intent(in) :: S
      !
      !
      ! ... Executable Statements ...
      !
      !
      C(1) = S(3)*cos(S(2))*cos(S(1))
      C(2) = S(3)*cos(S(2))*sin(S(1))
      C(3) = S(3)*sin(S(2))
      return
      !
   end subroutine sph2car

   subroutine updatrotmom_single(Mvect, MOM)
      !
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Parameters ..
      complex(selected_real_kind(15)), parameter :: cone = (1.0d+0, 0.0d+0)
      complex(selected_real_kind(15)), parameter :: czero = (0.0d+0, 0.0d+0)
      !
      !.. Formal Arguments ..
      real(selected_real_kind(15)), dimension(3), intent(inout) :: Mvect
      real(selected_real_kind(15)), dimension(3), intent(in) :: MOM
      !
      !
      !.. Local Scalars ..
      integer :: j
      real(selected_real_kind(15)) :: alfa, beta, theta
      !
      !.. Local Arrays ..
      real(selected_real_kind(15)), dimension(3) :: v, vout, sv, vrod, mz
      real(selected_real_kind(15)), dimension(3, 3) :: B, Rx, Ry, Rz, Rrod
      !complex(selected_real_kind(15)), dimension(3,3) :: Rz,Ry,R
      !
      !.. Intrinsic Functions ..
      intrinsic selected_real_kind
      !
      ! ... Executable Statements ...
      !
      !
      mz = (/0.0d0, 0.0d0, 1.0d0/)
      do j = 1, 3
         v(j) = MOM(j)
      end do
      ! Plan B) User Euler rotation
      !
      call car2sph(v, sv)
      alfa = 0.0d0*atan(1.0d0) + sv(2)
      beta = -0.0d0*atan(1.0d0) + sv(1)
      ! Rx
      Rx(1, 1) = 1.0d0
      Rx(2, 1) = 0.0d0
      Rx(3, 1) = 0.0d0
      Rx(1, 2) = 0.0d0
      Rx(2, 2) = cos(alfa)
      Rx(3, 2) = sin(alfa)
      Rx(1, 3) = 0.0d0
      Rx(2, 3) = -sin(alfa)
      Rx(3, 3) = cos(alfa)
      ! Ry
      Ry(1, 1) = cos(alfa)
      Ry(2, 1) = 0.0d0
      Ry(3, 1) = -sin(alfa)
      Ry(1, 2) = 0.0d0
      Ry(2, 2) = 1.0d0
      Ry(3, 2) = 0.0d0
      Ry(1, 3) = sin(alfa)
      Ry(2, 3) = 0.0d0
      Ry(3, 3) = cos(alfa)
    !! Rz
      Rz(1, 1) = cos(beta)
      Rz(2, 1) = sin(beta)
      Rz(3, 1) = 0.0d0
      Rz(1, 2) = -sin(beta)
      Rz(2, 2) = cos(beta)
      Rz(3, 2) = 0.0d0
      Rz(1, 3) = 0.0d0
      Rz(2, 3) = 0.0d0
      Rz(3, 3) = 1.0d0

      ! Rotate!
      v = mvect
      vout = matmul(Rz, matmul(Ry, v))

      do j = 1, 3
         mvect(j) = vout(j)
      end do
      return
      !
      ! ... Format Declarations ...
      !
10000 format(3x, a8, f10.4, a8, f10.4, a8, f10.4)
   end subroutine updatrotmom_single

   !----------------------------------------------------
   ! Set up the rotation matrix around a given axis
   ! based on Rodrigues' rotation formula
   !----------------------------------------------------
   subroutine setup_rotation_matrix(matrix, axis, theta)
      ! Declarations
      real(rp), intent(out) :: matrix(3, 3)
      real(rp), intent(in) :: axis(3), theta
      real(rp) :: ux, uy, uz, cosTheta, sinTheta

      ! Assign values
      ux = axis(1)
      uy = axis(2)
      uz = axis(3)
      cosTheta = cos(theta)
      sinTheta = sin(theta)

      ! Define the rotation matrix using Rodrigues' formula
      matrix(1, 1) = cosTheta + ux*ux*(1.0 - cosTheta)
      matrix(1, 2) = ux*uy*(1.0 - cosTheta) - uz*sinTheta
      matrix(1, 3) = ux*uz*(1.0 - cosTheta) + uy*sinTheta

      matrix(2, 1) = uy*ux*(1.0 - cosTheta) + uz*sinTheta
      matrix(2, 2) = cosTheta + uy*uy*(1.0 - cosTheta)
      matrix(2, 3) = uy*uz*(1.0 - cosTheta) - ux*sinTheta

      matrix(3, 1) = uz*ux*(1.0 - cosTheta) - uy*sinTheta
      matrix(3, 2) = uz*uy*(1.0 - cosTheta) + ux*sinTheta
      matrix(3, 3) = cosTheta + uz*uz*(1.0 - cosTheta)
   end subroutine

end module math_mod

