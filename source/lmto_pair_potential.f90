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
   public :: lmto_circular_pair_potential_from_reverse
   public :: lmto_bloch_phase
   public :: lmto_endpoint_phases
   public :: lmto_transition_metadata
   public :: lmto_unfold_site_spinors

   type, public :: lmto_pair_transition_metadata
      real(rp) :: k(3) = 0.0_rp
      real(rp) :: kq_unfolded(3) = 0.0_rp
      real(rp) :: kq_folded(3) = 0.0_rp
      integer :: reciprocal_shift(3) = 0
      logical :: unfolded_gauge_required = .false.
   end type lmto_pair_transition_metadata

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

   !> Form the emission-side circular vertex from a separately assembled
   !> reverse transition.  `dh_*_reverse` has rows in the k coefficient space
   !> and columns in the K=k+q coefficient space.  It is deliberately a
   !> different input from the absorption-side tangent so Q+ is not created by
   !> assigning the adjoint of Q-.
   subroutine lmto_circular_pair_potential_from_reverse(dh_dx_reverse, dh_dy_reverse, signed_moment, qplus, supported, reason)
      complex(rp), intent(in) :: dh_dx_reverse(:, :), dh_dy_reverse(:, :)
      real(rp), intent(in) :: signed_moment
      complex(rp), intent(out) :: qplus(:, :)
      logical, intent(out) :: supported
      character(len=*), intent(out), optional :: reason

      qplus = cmplx(0.0_rp, 0.0_rp, rp)
      supported = .false.
      if (size(dh_dx_reverse, 1) /= size(dh_dx_reverse, 2) .or. &
          any(shape(dh_dy_reverse) /= shape(dh_dx_reverse)) .or. any(shape(qplus) /= shape(dh_dx_reverse))) then
         if (present(reason)) reason = 'incompatible reverse Cartesian tangent matrix shapes'
         return
      end if
      if (abs(signed_moment) <= tiny(1.0_rp)) then
         if (present(reason)) reason = 'signed response moment is zero or unavailable'
         return
      end if
      qplus = (dh_dx_reverse + i_unit*dh_dy_reverse)/(2.0_rp*signed_moment)
      supported = .true.
      if (present(reason)) reason = 'ordinary ham_only LMTO reverse pair potential'
   end subroutine lmto_circular_pair_potential_from_reverse

   pure function lmto_bloch_phase(k_direct, displacement_direct) result(phase)
      real(rp), intent(in) :: k_direct(3), displacement_direct(3)
      complex(rp) :: phase
      real(rp) :: argument
      argument = 2.0_rp*acos(-1.0_rp)*dot_product(k_direct, displacement_direct)
      phase = cmplx(cos(argument), sin(argument), rp)
   end function lmto_bloch_phase

   pure subroutine lmto_endpoint_phases(k_direct, q_direct, displacement_direct, left_phase, right_phase)
      real(rp), intent(in) :: k_direct(3), q_direct(3), displacement_direct(3)
      complex(rp), intent(out) :: left_phase, right_phase
      left_phase = lmto_bloch_phase(k_direct, displacement_direct)
      right_phase = lmto_bloch_phase(k_direct + q_direct, displacement_direct)
   end subroutine lmto_endpoint_phases

   pure function lmto_transition_metadata(k_direct, q_direct) result(metadata)
      real(rp), intent(in) :: k_direct(3), q_direct(3)
      type(lmto_pair_transition_metadata) :: metadata
      metadata%k = k_direct
      metadata%kq_unfolded = k_direct + q_direct
      metadata%kq_folded = metadata%kq_unfolded - floor(metadata%kq_unfolded + 0.5_rp)
      metadata%reciprocal_shift = nint(metadata%kq_unfolded - metadata%kq_folded)
      metadata%unfolded_gauge_required = any(metadata%reciprocal_shift /= 0)
   end function lmto_transition_metadata

   subroutine lmto_unfold_site_spinors(metadata, tau_direct, spinors_folded, spinors_unfolded)
      type(lmto_pair_transition_metadata), intent(in) :: metadata
      real(rp), intent(in) :: tau_direct(:, :)
      complex(rp), intent(in) :: spinors_folded(:, :)
      complex(rp), intent(out) :: spinors_unfolded(:, :)
      integer :: nsite, nblock, isite, ibeg, iend
      complex(rp) :: phase
      nsite = size(tau_direct,2)
      if (size(tau_direct,1) /= 3 .or. nsite < 1 .or. mod(size(spinors_folded,1),nsite) /= 0 .or. &
          any(shape(spinors_unfolded) /= shape(spinors_folded))) error stop 'lmto_unfold_site_spinors: incompatible site gauge shape'
      nblock = size(spinors_folded,1)/nsite
      spinors_unfolded = spinors_folded
      do isite=1,nsite
         ibeg=(isite-1)*nblock+1; iend=isite*nblock
         phase = exp(cmplx(0.0_rp, -2.0_rp*acos(-1.0_rp)*dot_product(real(metadata%reciprocal_shift,rp),tau_direct(:,isite)), rp))
         spinors_unfolded(ibeg:iend,:) = phase*spinors_folded(ibeg:iend,:)
      end do
   end subroutine lmto_unfold_site_spinors

end module lmto_pair_potential_mod
