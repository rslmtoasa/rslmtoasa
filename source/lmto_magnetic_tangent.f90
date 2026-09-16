!------------------------------------------------------------------------------
! Ordinary noncollinear LMTO bond-value algebra.
! This module is deliberately independent of hamiltonian_mod so the
! ground-state Hamiltonian assembly has a small, side-effect-free seam.
!------------------------------------------------------------------------------
module lmto_magnetic_tangent_mod
   use precision_mod, only: rp
   use math_mod, only: i_unit
   implicit none
   private

   public :: lmto_bond_value
   public :: lmto_bond_derivative
   public :: lmto_bond_mixed_derivative
   public :: lmto_hhmag_to_spinor
   public :: lmto_transform_potential

contains

   ! The hhmag convention is (:,:,1:3) = (Hx,Hy,Hz), (:,:,4) = H0.
   pure subroutine lmto_bond_value(hhh, wx0_i, wx1_i, wx0_j, wx1_j, c0_i, c1_i, &
                                   mom_i, mom_j, onsite, hhmag)
      complex(rp), intent(in) :: hhh(:, :), wx0_i(:), wx1_i(:), wx0_j(:), wx1_j(:)
      complex(rp), intent(in) :: c0_i(:), c1_i(:)
      real(rp), intent(in) :: mom_i(3), mom_j(3)
      logical, intent(in) :: onsite
      complex(rp), intent(out) :: hhmag(:, :, :)
      integer :: ilm, jlm, idir
      complex(rp) :: dot, cross(3), momc_i(3), momc_j(3)

      hhmag = cmplx(0.0_rp, 0.0_rp, rp)
      dot = cmplx(dot_product(mom_i, mom_j), 0.0_rp, rp)
      momc_i = cmplx(mom_i, 0.0_rp, rp)
      momc_j = cmplx(mom_j, 0.0_rp, rp)
      cross = cmplx(cross3(mom_i, mom_j), 0.0_rp, rp)

      do jlm = 1, size(hhh, 2)
         do ilm = 1, size(hhh, 1)
            hhmag(ilm, jlm, 4) = wx0_i(ilm)*hhh(ilm, jlm)*wx0_j(jlm) + &
                                   wx1_i(ilm)*hhh(ilm, jlm)*wx1_j(jlm)*dot
            do idir = 1, 3
               hhmag(ilm, jlm, idir) = &
                  (wx1_i(ilm)*hhh(ilm, jlm)*wx0_j(jlm))*momc_i(idir) + &
                  (wx0_i(ilm)*hhh(ilm, jlm)*wx1_j(jlm))*momc_j(idir) + &
                  i_unit*wx1_i(ilm)*hhh(ilm, jlm)*wx1_j(jlm)*cross(idir)
            end do
         end do
      end do
      if (.not. onsite) return
      do ilm = 1, size(hhh, 1)
         hhmag(ilm, ilm, 4) = hhmag(ilm, ilm, 4) + c0_i(ilm)
         do idir = 1, 3
            hhmag(ilm, ilm, idir) = hhmag(ilm, ilm, idir) + c1_i(ilm)*momc_i(idir)
         end do
      end do
   end subroutine lmto_bond_value

   !> First derivative of lmto_bond_value with respect to independent endpoint
   !> moment variations.  The two derivative vectors are deliberately passed
   !> independently: an onsite block uses the same local variation at both
   !> endpoints, while a directed bond normally varies only one endpoint.
   pure subroutine lmto_bond_derivative(hhh, wx0_i, wx1_i, wx0_j, wx1_j, c1_i, &
                                        mom_i, mom_j, dm_i, dm_j, onsite, dhmag)
      complex(rp), intent(in) :: hhh(:, :), wx0_i(:), wx1_i(:), wx0_j(:), wx1_j(:), c1_i(:)
      real(rp), intent(in) :: mom_i(3), mom_j(3), dm_i(3), dm_j(3)
      logical, intent(in) :: onsite
      complex(rp), intent(out) :: dhmag(:, :, :)
      real(rp) :: dcross(3), ddot
      integer :: ilm, jlm, idir

      dhmag = cmplx(0.0_rp, 0.0_rp, rp)
      ddot = dot_product(dm_i, mom_j) + dot_product(mom_i, dm_j)
      dcross = cross3(dm_i, mom_j) + cross3(mom_i, dm_j)
      do jlm = 1, size(hhh, 2)
         do ilm = 1, size(hhh, 1)
            dhmag(ilm, jlm, 4) = wx1_i(ilm)*hhh(ilm, jlm)*wx1_j(jlm)*ddot
            do idir = 1, 3
               dhmag(ilm, jlm, idir) = wx1_i(ilm)*hhh(ilm, jlm)*wx0_j(jlm)*dm_i(idir) + &
                  wx0_i(ilm)*hhh(ilm, jlm)*wx1_j(jlm)*dm_j(idir) + &
                  i_unit*wx1_i(ilm)*hhh(ilm, jlm)*wx1_j(jlm)*dcross(idir)
            end do
         end do
      end do
      if (.not. onsite) return
      do ilm = 1, min(size(hhh, 1), size(c1_i))
         do idir = 1, 3
            dhmag(ilm, ilm, idir) = dhmag(ilm, ilm, idir) + &
               c1_i(ilm)*cmplx(dm_i(idir), 0.0_rp, rp)
         end do
      end do
   end subroutine lmto_bond_derivative

   !> Mixed derivative for two independent local rotations.  This is the
   !> derivative of the complete dot/cross bond algebra, not a derivative of
   !> only its spin-diagonal block.  The onsite c1*m correction is linear in
   !> the local moment and therefore has no mixed derivative for i/=j.
   pure subroutine lmto_bond_mixed_derivative(hhh, wx0_i, wx1_i, wx0_j, wx1_j, &
                                               mom_i, mom_j, dm_i_a, dm_j_a, dm_i_b, dm_j_b, d2hmag)
      complex(rp), intent(in) :: hhh(:, :), wx0_i(:), wx1_i(:), wx0_j(:), wx1_j(:)
      real(rp), intent(in) :: mom_i(3), mom_j(3), dm_i_a(3), dm_j_a(3), dm_i_b(3), dm_j_b(3)
      complex(rp), intent(out) :: d2hmag(:, :, :)
      real(rp) :: d2cross(3), d2dot
      integer :: ilm, jlm, idir

      d2hmag = cmplx(0.0_rp, 0.0_rp, rp)
      d2dot = dot_product(dm_i_a, dm_j_b) + dot_product(dm_i_b, dm_j_a)
      d2cross = cross3(dm_i_a, dm_j_b) + cross3(dm_i_b, dm_j_a)
      do jlm = 1, size(hhh, 2)
         do ilm = 1, size(hhh, 1)
            d2hmag(ilm, jlm, 4) = wx1_i(ilm)*hhh(ilm, jlm)*wx1_j(jlm)*d2dot
            do idir = 1, 3
               d2hmag(ilm, jlm, idir) = i_unit*wx1_i(ilm)*hhh(ilm, jlm)*wx1_j(jlm)*d2cross(idir)
            end do
         end do
      end do
   end subroutine lmto_bond_mixed_derivative

   ! Convert the ordinary four-channel bond value into a spinor block for
   ! ground-state noncollinear and GBT Hamiltonian consumers.
   pure subroutine lmto_hhmag_to_spinor(hhmag, spinor)
      complex(rp), intent(in) :: hhmag(:, :, :)
      complex(rp), intent(out) :: spinor(:, :)
      integer :: i, j, nrow, ncol

      nrow = size(hhmag, 1)
      ncol = size(hhmag, 2)
      spinor = cmplx(0.0_rp, 0.0_rp, rp)
      do j = 1, ncol
         do i = 1, nrow
            spinor(i, j) = hhmag(i, j, 4) + hhmag(i, j, 3)
            spinor(i + nrow, j + ncol) = hhmag(i, j, 4) - hhmag(i, j, 3)
            spinor(i, j + ncol) = hhmag(i, j, 1) - i_unit*hhmag(i, j, 2)
            spinor(i + nrow, j) = hhmag(i, j, 1) + i_unit*hhmag(i, j, 2)
         end do
      end do
   end subroutine lmto_hhmag_to_spinor

   pure function cross3(a, b) result(c)
      real(rp), intent(in) :: a(3), b(3)
      real(rp) :: c(3)
      c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
   end function cross3

   !> Reproduce the live symbolic-atom predls transformation without creating
   !> a symbolic_atom object.  This is used by controlled DRESP fixtures so
   !> their wx/cx/obx channels originate from one LMTO potential data set.
   pure subroutine lmto_transform_potential(c, enu, srdel, qpar, vmad, ws_radius, wsm, qm, &
                                            center_band, width_band, shifted_band, obar, dele_out)
      real(rp), intent(in) :: c(:, :), enu(:, :), srdel(:, :), qpar(:, :), vmad, ws_radius, wsm, qm(:)
      real(rp), intent(out) :: center_band(:, :), width_band(:, :), shifted_band(:, :), obar(:, :), dele_out(:, :)
      real(rp) :: wow, x, y, qi
      integer :: i, s, lmax

      if (size(c, 1) /= size(enu, 1) .or. size(c, 1) /= size(srdel, 1) .or. size(c, 1) /= size(qpar, 1) .or. &
          size(c, 2) /= 2 .or. size(enu, 2) /= 2 .or. size(srdel, 2) /= 2 .or. size(qpar, 2) /= 2 .or. &
          size(center_band, 1) /= size(c, 1) .or. size(width_band, 1) /= size(c, 1) .or. &
          size(shifted_band, 1) /= size(c, 1) .or. size(obar, 1) /= size(c, 1) .or. &
          size(dele_out, 1) /= size(c, 1) .or. size(center_band, 2) /= 2 .or. size(width_band, 2) /= 2 .or. &
          size(shifted_band, 2) /= 2 .or. size(obar, 2) /= 2 .or. size(dele_out, 2) /= 2 .or. &
          size(qm) < size(c, 1) .or. ws_radius <= 0.0_rp .or. wsm <= 0.0_rp) then
         error stop 'lmto_transform_potential: inconsistent LMTO channel shapes'
      end if
      wow = wsm/ws_radius
      lmax = size(c, 1) - 1
      do s = 1, 2
         do i = 1, lmax + 1
            qi = qpar(i, s)*wow**(1 - 2*i)
            dele_out(i, s) = srdel(i, s)*wow**(0.5_rp - i)
            x = 1.0_rp - ((qi - qm(i))*(c(i, s) - enu(i, s))/dele_out(i, s)**2)
            y = (qi - qm(i))/(((c(i, s) - enu(i, s))*(qi - qm(i))) - dele_out(i, s)**2)
            center_band(i, s) = (c(i, s) - enu(i, s))*x + enu(i, s) + vmad
            shifted_band(i, s) = (c(i, s) - enu(i, s))*x
            width_band(i, s) = dele_out(i, s)*x
            obar(i, s) = y
         end do
      end do
   end subroutine lmto_transform_potential

end module lmto_magnetic_tangent_mod
