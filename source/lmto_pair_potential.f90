!------------------------------------------------------------------------------
! Pair-potential algebra shared by the reciprocal LMTO provider and its oracle.
! `signed_moment` is deliberately an explicit input: potential%mom is only a
! unit orientation and must never be used as the response normalization.
!------------------------------------------------------------------------------
module lmto_pair_potential_mod
   use precision_mod, only: rp
   use math_mod, only: i_unit
   implicit none
   private
   public :: lmto_circular_pair_potential
   public :: lmto_bloch_phase

contains

   subroutine lmto_circular_pair_potential(dh_dx, dh_dy, signed_moment, qminus, qplus, supported, reason)
      complex(rp), intent(in) :: dh_dx(:, :), dh_dy(:, :)
      real(rp), intent(in) :: signed_moment
      complex(rp), intent(out) :: qminus(:, :), qplus(:, :)
      logical, intent(out) :: supported
      character(len=*), intent(out), optional :: reason

      qminus = cmplx(0.0_rp, 0.0_rp, rp); qplus = cmplx(0.0_rp, 0.0_rp, rp)
      supported = .false.
      if (size(dh_dx, 1) /= size(dh_dx, 2) .or. any(shape(dh_dy) /= shape(dh_dx)) .or. &
          any(shape(qminus) /= shape(dh_dx)) .or. any(shape(qplus) /= shape(dh_dx))) then
         if (present(reason)) reason = 'incompatible Cartesian tangent matrix shapes'
         return
      end if
      if (abs(signed_moment) <= tiny(1.0_rp)) then
         if (present(reason)) reason = 'signed response moment is zero or unavailable'
         return
      end if
      qminus = (dh_dx - i_unit*dh_dy)/(2.0_rp*signed_moment)
      qplus = transpose(conjg(qminus))
      supported = .true.
      if (present(reason)) reason = 'ordinary ham_only LMTO pair potential'
   end subroutine lmto_circular_pair_potential

   pure function lmto_bloch_phase(k_direct, displacement_direct) result(phase)
      real(rp), intent(in) :: k_direct(3), displacement_direct(3)
      complex(rp) :: phase
      real(rp) :: argument
      argument = 2.0_rp*acos(-1.0_rp)*dot_product(k_direct, displacement_direct)
      phase = cmplx(cos(argument), sin(argument), rp)
   end function lmto_bloch_phase

end module lmto_pair_potential_mod
