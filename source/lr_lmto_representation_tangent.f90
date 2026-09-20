!------------------------------------------------------------------------------
! DRESP-09U production LMTO representation map and analytic tangent.
!
! This module mirrors the scalar algebra in symbolic_atom%predls without
! changing predls itself.  The second half mirrors the structured spin lift
! used by build_pot, build_obarm, build_enim, and chbar_nc.  It is deliberately
! not a generic Lowdin/S^{-1/2} service: the live production map is a map of
! per-spin LMTO potential channels followed by a Pauli spin lift.
!------------------------------------------------------------------------------
module lr_lmto_representation_tangent_mod

   use precision_mod, only: rp
   use basis_mod, only: nb, norb, spin_off
   use math_mod, only: hcpx, i_unit
   use lmto_magnetic_tangent_mod, only: lmto_hhmag_to_spinor
   implicit none
   private

   public :: lmto_predls_transform
   public :: lmto_predls_transform_tangent
   public :: lmto_spin_parameter_tangent
   public :: lmto_representation_bond_tangent
   public :: lmto_representation_relative_residual

contains

   !> Exact scalar predls map, with the same 1-based local channel storage as
   !> symbolic_atom%predls after it copies c/enu/srdel/qpar(0:lmax,:).
   pure subroutine lmto_predls_transform(c, enu, srdel, qpar, vmad, ws_radius, wsm, qm, &
                                         center_band, shifted_band, width_band, obar, dele_out, qi_out)
      real(rp), intent(in) :: c(:, :), enu(:, :), srdel(:, :), qpar(:, :)
      real(rp), intent(in) :: vmad, ws_radius, wsm, qm(:)
      real(rp), intent(out) :: center_band(:, :), shifted_band(:, :), width_band(:, :), obar(:, :)
      real(rp), intent(out), optional :: dele_out(:, :), qi_out(:, :)
      real(rp) :: wow, dele, qi, x, y, delta
      integer :: i, spin, nchannel

      call check_predls_shapes(c, enu, srdel, qpar, qm, center_band, shifted_band, width_band, obar, &
         ws_radius, wsm)
      nchannel = size(c, 1)
      wow = wsm/ws_radius
      do spin = 1, 2
         do i = 1, nchannel
            dele = srdel(i, spin)*wow**(0.5_rp - real(i, rp))
            qi = qpar(i, spin)*wow**(1.0_rp - 2.0_rp*real(i, rp))
            delta = c(i, spin) - enu(i, spin)
            x = 1.0_rp - (qi - qm(i))*delta/(dele*dele)
            y = (qi - qm(i))/((delta*(qi - qm(i))) - dele*dele)
            center_band(i, spin) = delta*x + enu(i, spin) + vmad
            shifted_band(i, spin) = delta*x
            width_band(i, spin) = dele*x
            obar(i, spin) = y
            if (present(dele_out)) dele_out(i, spin) = dele
            if (present(qi_out)) qi_out(i, spin) = qi
         end do
      end do
   end subroutine lmto_predls_transform

   !> Analytic derivative of the exact live predls equations.
   !>
   !> D = C-E_nu, A = QI-QM, N = D*A-DELE^2,
   !> X = 1-A*D/DELE^2, Y=A/N.
   !> The screening constants and wow are fixed representation inputs here;
   !> their possible self-consistent response is intentionally outside the
   !> rigid local-spin-frame rotation gate.
   pure subroutine lmto_predls_transform_tangent(c, enu, srdel, qpar, dc, denu, dsrdel, dqpar, &
                                                 vmad, ws_radius, wsm, qm, center_band, shifted_band, width_band, obar, &
                                                 dcenter_band, dshifted_band, dwidth_band, dobar, dele_out, ddele_out, qi_out, dqi_out)
      real(rp), intent(in) :: c(:, :), enu(:, :), srdel(:, :), qpar(:, :)
      real(rp), intent(in) :: dc(:, :), denu(:, :), dsrdel(:, :), dqpar(:, :)
      real(rp), intent(in) :: vmad, ws_radius, wsm, qm(:)
      real(rp), intent(out) :: center_band(:, :), shifted_band(:, :), width_band(:, :), obar(:, :)
      real(rp), intent(out) :: dcenter_band(:, :), dshifted_band(:, :), dwidth_band(:, :), dobar(:, :)
      real(rp), intent(out), optional :: dele_out(:, :), ddele_out(:, :), qi_out(:, :), dqi_out(:, :)
      real(rp) :: wow, factor_dele, factor_qi
      real(rp) :: dele, ddele, qi, dqi, delta, ddelta, a, da, x, dx, denominator, ddenominator, y, dy
      integer :: i, spin, nchannel

      call check_predls_tangent_shapes(c, enu, srdel, qpar, dc, denu, dsrdel, dqpar, qm, center_band, shifted_band, &
         width_band, obar, dcenter_band, dshifted_band, dwidth_band, dobar, ws_radius, wsm)
      nchannel = size(c, 1)
      wow = wsm/ws_radius
      do spin = 1, 2
         do i = 1, nchannel
            factor_dele = wow**(0.5_rp - real(i, rp))
            factor_qi = wow**(1.0_rp - 2.0_rp*real(i, rp))
            dele = srdel(i, spin)*factor_dele
            ddele = dsrdel(i, spin)*factor_dele
            qi = qpar(i, spin)*factor_qi
            dqi = dqpar(i, spin)*factor_qi
            delta = c(i, spin) - enu(i, spin)
            ddelta = dc(i, spin) - denu(i, spin)
            a = qi - qm(i)
            da = dqi
            x = 1.0_rp - a*delta/(dele*dele)
            dx = -(da*delta + a*ddelta)/(dele*dele) + 2.0_rp*a*delta*ddele/(dele**3)
            denominator = delta*a - dele*dele
            ddenominator = ddelta*a + delta*da - 2.0_rp*dele*ddele
            y = a/denominator
            dy = da/denominator - a*ddenominator/(denominator*denominator)
            center_band(i, spin) = delta*x + enu(i, spin) + vmad
            shifted_band(i, spin) = delta*x
            width_band(i, spin) = dele*x
            obar(i, spin) = y
            dcenter_band(i, spin) = ddelta*x + delta*dx + denu(i, spin)
            dshifted_band(i, spin) = ddelta*x + delta*dx
            dwidth_band(i, spin) = ddele*x + dele*dx
            dobar(i, spin) = dy
            if (present(dele_out)) dele_out(i, spin) = dele
            if (present(ddele_out)) ddele_out(i, spin) = ddele
            if (present(qi_out)) qi_out(i, spin) = qi
            if (present(dqi_out)) dqi_out(i, spin) = dqi
         end do
      end do
   end subroutine lmto_predls_transform_tangent

   !> Lift two transformed spin channels into the exact production spinor and
   !> its tangent.  Channel derivatives are formed before this lift; for a
   !> rigid local-frame rotation they are zero and only dmoment is nonzero.
   subroutine lmto_spin_parameter_tangent(up, down, dup, ddown, moment, dmoment, matrix, delta_matrix)
      complex(rp), intent(in) :: up(:), down(:), dup(:), ddown(:)
      real(rp), intent(in) :: moment(3), dmoment(3)
      complex(rp), intent(out) :: matrix(:, :), delta_matrix(:, :)
      complex(rp) :: scalar, vector, dscalar, dvector
      integer :: i

      if (size(up) /= norb .or. size(down) /= norb .or. size(dup) /= norb .or. size(ddown) /= norb .or. &
          any(shape(matrix) /= [nb, nb]) .or. any(shape(delta_matrix) /= [nb, nb])) then
         error stop 'DRESP-09U spin lift: inconsistent dimensions'
      end if
      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      delta_matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, norb
         scalar = 0.5_rp*(up(i) + down(i))
         vector = 0.5_rp*(up(i) - down(i))
         dscalar = 0.5_rp*(dup(i) + ddown(i))
         dvector = 0.5_rp*(dup(i) - ddown(i))
         matrix(i, i) = scalar + vector*cmplx(moment(3), 0.0_rp, rp)
         matrix(i + spin_off, i + spin_off) = scalar - vector*cmplx(moment(3), 0.0_rp, rp)
         matrix(i, i + spin_off) = vector*cmplx(moment(1), 0.0_rp, rp) - i_unit*vector*cmplx(moment(2), 0.0_rp, rp)
         matrix(i + spin_off, i) = vector*cmplx(moment(1), 0.0_rp, rp) + i_unit*vector*cmplx(moment(2), 0.0_rp, rp)
         delta_matrix(i, i) = dscalar + dvector*cmplx(moment(3), 0.0_rp, rp) + vector*cmplx(dmoment(3), 0.0_rp, rp)
         delta_matrix(i + spin_off, i + spin_off) = dscalar - dvector*cmplx(moment(3), 0.0_rp, rp) - vector*cmplx(dmoment(3), 0.0_rp, rp)
         delta_matrix(i, i + spin_off) = dvector*cmplx(moment(1), 0.0_rp, rp) + vector*cmplx(dmoment(1), 0.0_rp, rp) - &
            i_unit*(dvector*cmplx(moment(2), 0.0_rp, rp) + vector*cmplx(dmoment(2), 0.0_rp, rp))
         delta_matrix(i + spin_off, i) = dvector*cmplx(moment(1), 0.0_rp, rp) + vector*cmplx(dmoment(1), 0.0_rp, rp) + &
            i_unit*(dvector*cmplx(moment(2), 0.0_rp, rp) + vector*cmplx(dmoment(2), 0.0_rp, rp))
      end do
      call hcpx(matrix(1:norb, 1:norb), 'cart2sph')
      call hcpx(matrix(norb+1:nb, norb+1:nb), 'cart2sph')
      call hcpx(matrix(1:norb, norb+1:nb), 'cart2sph')
      call hcpx(matrix(norb+1:nb, 1:norb), 'cart2sph')
      call hcpx(delta_matrix(1:norb, 1:norb), 'cart2sph')
      call hcpx(delta_matrix(norb+1:nb, norb+1:nb), 'cart2sph')
      call hcpx(delta_matrix(1:norb, norb+1:nb), 'cart2sph')
      call hcpx(delta_matrix(norb+1:nb, 1:norb), 'cart2sph')
   end subroutine lmto_spin_parameter_tangent

   !> Analytic derivative of the live lmto_bond_value algebra.  The endpoint
   !> factors are the transformed wx channels and c channels are cex for HOH
   !> or cx for the first-order path.  No native tangent is used here.
   subroutine lmto_representation_bond_tangent(hhh, wx0_i, wx1_i, wx0_j, wx1_j, c0_i, c1_i, &
                                                     dwx0_i, dwx1_i, dwx0_j, dwx1_j, dc0_i, dc1_i, &
                                                     mom_i, mom_j, dm_i, dm_j, onsite, delta_spinor)
      complex(rp), intent(in) :: hhh(:, :), wx0_i(:), wx1_i(:), wx0_j(:), wx1_j(:), c0_i(:), c1_i(:)
      complex(rp), intent(in) :: dwx0_i(:), dwx1_i(:), dwx0_j(:), dwx1_j(:), dc0_i(:), dc1_i(:)
      real(rp), intent(in) :: mom_i(3), mom_j(3), dm_i(3), dm_j(3)
      logical, intent(in) :: onsite
      complex(rp), intent(out) :: delta_spinor(:, :)
      complex(rp), allocatable :: dhmag(:, :, :)
      complex(rp) :: a, b, da, db, left, right, dleft, dright, crossc(3)
      real(rp) :: dcross(3)
      real(rp) :: dot, ddot
      integer :: i, j, idir

      if (size(wx0_i) /= size(hhh, 1) .or. size(wx1_i) /= size(hhh, 1) .or. size(wx0_j) /= size(hhh, 2) .or. &
          size(wx1_j) /= size(hhh, 2) .or. any(shape(delta_spinor) /= [2*size(hhh,1),2*size(hhh,2)])) then
         error stop 'DRESP-09U bond tangent: inconsistent dimensions'
      end if
      allocate(dhmag(size(hhh,1), size(hhh,2), 4))
      dhmag = cmplx(0.0_rp, 0.0_rp, rp)
      dot = dot_product(mom_i, mom_j)
      ddot = dot_product(dm_i, mom_j) + dot_product(mom_i, dm_j)
      crossc = cmplx([mom_i(2)*mom_j(3)-mom_i(3)*mom_j(2), mom_i(3)*mom_j(1)-mom_i(1)*mom_j(3), &
         mom_i(1)*mom_j(2)-mom_i(2)*mom_j(1)], 0.0_rp, rp)
      dcross = [dm_i(2)*mom_j(3)-dm_i(3)*mom_j(2)+mom_i(2)*dm_j(3)-mom_i(3)*dm_j(2), &
         dm_i(3)*mom_j(1)-dm_i(1)*mom_j(3)+mom_i(3)*dm_j(1)-mom_i(1)*dm_j(3), &
         dm_i(1)*mom_j(2)-dm_i(2)*mom_j(1)+mom_i(1)*dm_j(2)-mom_i(2)*dm_j(1)]
      do j = 1, size(hhh, 2)
         do i = 1, size(hhh, 1)
            a = wx0_i(i)*hhh(i,j)*wx0_j(j)
            b = wx1_i(i)*hhh(i,j)*wx1_j(j)
            da = dwx0_i(i)*hhh(i,j)*wx0_j(j) + wx0_i(i)*hhh(i,j)*dwx0_j(j)
            db = dwx1_i(i)*hhh(i,j)*wx1_j(j) + wx1_i(i)*hhh(i,j)*dwx1_j(j)
            dhmag(i,j,4) = da + db*dot + b*ddot
            left = wx1_i(i)*hhh(i,j)*wx0_j(j)
            right = wx0_i(i)*hhh(i,j)*wx1_j(j)
            dleft = dwx1_i(i)*hhh(i,j)*wx0_j(j) + wx1_i(i)*hhh(i,j)*dwx0_j(j)
            dright = dwx0_i(i)*hhh(i,j)*wx1_j(j) + wx0_i(i)*hhh(i,j)*dwx1_j(j)
            do idir = 1, 3
               dhmag(i,j,idir) = dleft*cmplx(mom_i(idir),0.0_rp,rp) + left*cmplx(dm_i(idir),0.0_rp,rp) + &
                  dright*cmplx(mom_j(idir),0.0_rp,rp) + right*cmplx(dm_j(idir),0.0_rp,rp) + &
                  i_unit*(db*crossc(idir) + b*cmplx(dcross(idir),0.0_rp,rp))
            end do
         end do
      end do
      if (onsite) then
         do i = 1, min(size(hhh,1),size(c0_i))
            dhmag(i,i,4) = dhmag(i,i,4) + dc0_i(i)
            do idir = 1, 3
               dhmag(i,i,idir) = dhmag(i,i,idir) + dc1_i(i)*cmplx(mom_i(idir),0.0_rp,rp) + &
                  c1_i(i)*cmplx(dm_i(idir),0.0_rp,rp)
            end do
         end do
      end if
      do idir = 1, 4
         call hcpx(dhmag(:,:,idir), 'cart2sph')
      end do
      call lmto_hhmag_to_spinor(dhmag, delta_spinor)
      deallocate(dhmag)
   end subroutine lmto_representation_bond_tangent

   pure real(rp) function lmto_representation_relative_residual(left, right) result(value)
      complex(rp), intent(in) :: left(:,:), right(:,:)
      value = sqrt(sum(abs(left-right)**2))/max(sqrt(sum(abs(right)**2)),tiny(1.0_rp))
   end function lmto_representation_relative_residual

   pure subroutine check_predls_shapes(c, enu, srdel, qpar, qm, center, shifted, width, obar, ws_radius, wsm)
      real(rp), intent(in) :: c(:,:), enu(:,:), srdel(:,:), qpar(:,:), qm(:), center(:,:), shifted(:,:), width(:,:), obar(:,:)
      real(rp), intent(in) :: ws_radius, wsm
      if (size(c,1) < 1 .or. size(c,2) /= 2 .or. any(shape(enu) /= shape(c)) .or. any(shape(srdel) /= shape(c)) .or. &
          any(shape(qpar) /= shape(c)) .or. size(qm) < size(c,1) .or. any(shape(center) /= shape(c)) .or. &
          any(shape(shifted) /= shape(c)) .or. any(shape(width) /= shape(c)) .or. any(shape(obar) /= shape(c)) .or. &
          ws_radius <= 0.0_rp .or. wsm <= 0.0_rp) error stop 'DRESP-09U predls: inconsistent channel shapes'
   end subroutine check_predls_shapes

   pure subroutine check_predls_tangent_shapes(c, enu, srdel, qpar, dc, denu, dsrdel, dqpar, qm, center, shifted, width, obar, &
                                               dcenter, dshifted, dwidth, dobar, ws_radius, wsm)
      real(rp), intent(in) :: c(:,:), enu(:,:), srdel(:,:), qpar(:,:), dc(:,:), denu(:,:), dsrdel(:,:), dqpar(:,:), qm(:)
      real(rp), intent(in) :: center(:,:), shifted(:,:), width(:,:), obar(:,:), dcenter(:,:), dshifted(:,:), dwidth(:,:), dobar(:,:)
      real(rp), intent(in) :: ws_radius, wsm
      call check_predls_shapes(c, enu, srdel, qpar, qm, center, shifted, width, obar, ws_radius, wsm)
      if (any(shape(dc) /= shape(c)) .or. any(shape(denu) /= shape(c)) .or. any(shape(dsrdel) /= shape(c)) .or. &
          any(shape(dqpar) /= shape(c)) .or. any(shape(dcenter) /= shape(c)) .or. any(shape(dshifted) /= shape(c)) .or. &
          any(shape(dwidth) /= shape(c)) .or. any(shape(dobar) /= shape(c))) error stop 'DRESP-09U predls tangent: inconsistent shapes'
   end subroutine check_predls_tangent_shapes

end module lr_lmto_representation_tangent_mod
