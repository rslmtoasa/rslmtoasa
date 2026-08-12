!------------------------------------------------------------------------------
! Endpoint-resolved first variations of the ordinary noncollinear LMTO bond.
! This module is deliberately independent of hamiltonian_mod so that the
! value-and-tangent algebra has a small, side-effect-free test seam.
!------------------------------------------------------------------------------
module lmto_magnetic_tangent_mod
   use precision_mod, only: rp
   use math_mod, only: i_unit
   implicit none
   private

   public :: lmto_bond_value
   public :: lmto_bond_tangent
   public :: lmto_hhmag_to_spinor
   public :: lmto_make_endpoint_record
   public :: lmto_ordinary_tangent_supported

   type, public :: lmto_endpoint_tangent_record
      ! Response identities are site identities, never chemical-type keys.
      integer :: source_site = 0
      integer :: neighbor_site = 0
      integer :: source_type = 0
      integer :: neighbor_type = 0
      integer :: directed_bond = 0
      integer :: operator_generation = 0
      real(rp) :: displacement(3) = 0.0_rp
      logical :: onsite_owned_by_source = .false.
      logical :: supported = .false.
      character(len=32) :: provenance = 'ham0m_nc_endpoint_tangent_v1'
   end type lmto_endpoint_tangent_record

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

   ! Complete directional derivative at fixed LMTO potential parameters.
   pure subroutine lmto_bond_tangent(hhh, wx0_i, wx1_i, wx0_j, wx1_j, c1_i, &
                                     mom_i, mom_j, delta_mom_i, delta_mom_j, onsite, delta_hhmag)
      complex(rp), intent(in) :: hhh(:, :), wx0_i(:), wx1_i(:), wx0_j(:), wx1_j(:), c1_i(:)
      real(rp), intent(in) :: mom_i(3), mom_j(3), delta_mom_i(3), delta_mom_j(3)
      logical, intent(in) :: onsite
      complex(rp), intent(out) :: delta_hhmag(:, :, :)
      integer :: ilm, jlm, idir
      complex(rp) :: delta_dot, delta_cross(3), delta_momc_i(3), delta_momc_j(3)

      delta_hhmag = cmplx(0.0_rp, 0.0_rp, rp)
      delta_dot = cmplx(dot_product(delta_mom_i, mom_j) + dot_product(mom_i, delta_mom_j), 0.0_rp, rp)
      delta_momc_i = cmplx(delta_mom_i, 0.0_rp, rp)
      delta_momc_j = cmplx(delta_mom_j, 0.0_rp, rp)
      delta_cross = cmplx(cross3(delta_mom_i, mom_j) + cross3(mom_i, delta_mom_j), 0.0_rp, rp)

      do jlm = 1, size(hhh, 2)
         do ilm = 1, size(hhh, 1)
            delta_hhmag(ilm, jlm, 4) = wx1_i(ilm)*hhh(ilm, jlm)*wx1_j(jlm)*delta_dot
            do idir = 1, 3
               delta_hhmag(ilm, jlm, idir) = &
                  (wx1_i(ilm)*hhh(ilm, jlm)*wx0_j(jlm))*delta_momc_i(idir) + &
                  (wx0_i(ilm)*hhh(ilm, jlm)*wx1_j(jlm))*delta_momc_j(idir) + &
                  i_unit*wx1_i(ilm)*hhh(ilm, jlm)*wx1_j(jlm)*delta_cross(idir)
            end do
         end do
      end do
      if (.not. onsite) return
      do ilm = 1, size(hhh, 1)
         do idir = 1, 3
            delta_hhmag(ilm, ilm, idir) = delta_hhmag(ilm, ilm, idir) + c1_i(ilm)*delta_momc_i(idir)
         end do
      end do
   end subroutine lmto_bond_tangent

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

   pure function lmto_make_endpoint_record(source_site, neighbor_site, source_type, neighbor_type, &
                                           directed_bond, operator_generation, displacement, supported) result(record)
      integer, intent(in) :: source_site, neighbor_site, source_type, neighbor_type
      integer, intent(in) :: directed_bond, operator_generation
      real(rp), intent(in) :: displacement(3)
      logical, intent(in) :: supported
      type(lmto_endpoint_tangent_record) :: record
      record%source_site = source_site; record%neighbor_site = neighbor_site
      record%source_type = source_type; record%neighbor_type = neighbor_type
      record%directed_bond = directed_bond; record%operator_generation = operator_generation
      record%displacement = displacement; record%onsite_owned_by_source = norm2(displacement) <= 0.01_rp
      record%supported = supported
   end function lmto_make_endpoint_record

   pure logical function lmto_ordinary_tangent_supported(is_gbt, has_hoh, has_ccor, has_hubbard, has_local_axis, &
                                                         has_soc, has_external_field)
      logical, intent(in) :: is_gbt, has_hoh, has_ccor, has_hubbard, has_local_axis, has_soc, has_external_field
      lmto_ordinary_tangent_supported = .not. (is_gbt .or. has_hoh .or. has_ccor .or. has_hubbard .or. has_local_axis .or. &
                                                has_soc .or. has_external_field)
   end function lmto_ordinary_tangent_supported

   pure function cross3(a, b) result(c)
      real(rp), intent(in) :: a(3), b(3)
      real(rp) :: c(3)
      c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
   end function cross3

end module lmto_magnetic_tangent_mod
