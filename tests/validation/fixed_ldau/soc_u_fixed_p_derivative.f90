! dE_Hub/dU and dE_Hub/dJ (Ry/Ry) of the spinor d-shell functional at a fixed P: central differences of
! hubbard_spinor_functional at U +/- eps and J +/- eps. E_Hub is linear in U and J at fixed P, so the
! result does not depend on eps. Input file: 100 lines "Re Im" of P(a,b), a major, a,b = 1..10 (a fast in b).
! Usage: soc_u_fixed_p_derivative P.txt U_Ry J_Ry eps_Ry
program soc_u_fixed_p_derivative
   use precision_mod, only: rp
   use hubbard_spinor_mod, only: hubbard_d_tensor, hubbard_spinor_functional
   implicit none
   complex(rp) :: p(10, 10), t(5, 5, 5, 5), v(10, 10)
   real(rp) :: u, j, eps, re, im, e_int, e_dc, ep(2), em(2), e0
   character(len=256) :: arg
   integer :: a, b, unit

   call get_command_argument(1, arg)
   open (newunit=unit, file=trim(arg), status='old', action='read')
   do a = 1, 10
      do b = 1, 10
         read (unit, *) re, im
         p(a, b) = cmplx(re, im, rp)
      end do
   end do
   close (unit)
   call get_command_argument(2, arg); read (arg, *) u
   call get_command_argument(3, arg); read (arg, *) j
   call get_command_argument(4, arg); read (arg, *) eps
   call energy(u, j, e0)
   call energy(u + eps, j, ep(1)); call energy(u - eps, j, em(1))
   call energy(u, j + eps, ep(2)); call energy(u, j - eps, em(2))
   write (*, '(a,es24.16)') 'E_Hub  ', e0
   write (*, '(a,es24.16)') 'dE_dU  ', (ep(1) - em(1))/(2.0_rp*eps)
   write (*, '(a,es24.16)') 'dE_dJ  ', (ep(2) - em(2))/(2.0_rp*eps)

contains

   subroutine energy(uu, jj, e)
      real(rp), intent(in) :: uu, jj
      real(rp), intent(out) :: e

      call hubbard_d_tensor(uu, jj, t)
      call hubbard_spinor_functional(t, uu, jj, p, e_int, e_dc, e, v)
   end subroutine energy

end program soc_u_fixed_p_derivative
