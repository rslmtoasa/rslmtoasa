!------------------------------------------------------------------------------
! Deterministic second-order radial data for algebraic LR unit tests.
!
! This fixture intentionally does not call the legacy radial ODE solver.  It
! supplies the complete accepted lmto_radial_basis contract so product-space
! and kernel tests exercise LR algebra without depending on RSEQSR/PHDFSR.
!------------------------------------------------------------------------------
module lr_test_radial_fixture_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, scalar_relativistic_c
   use linear_response_mod, only: response_space_layout, lmto_product_branch_00, lmto_product_branch_10, &
      lmto_product_branch_01, lmto_product_branch_11, lmto_product_branch_20, lmto_product_branch_02, &
      lmto_product_second_order_radial_branch
   implicit none
   private

   public :: build_lr_test_radial_basis
   public :: assert_lr_test_product_branch_norms

contains

   subroutine build_lr_test_radial_basis(basis, radius, lmax, mesh_a, mesh_b, nuclear_z)
      type(lmto_radial_basis), intent(out) :: basis
      real(rp), intent(in) :: radius(:), mesh_a, mesh_b, nuclear_z
      integer, intent(in) :: lmax
      integer :: ir, l, ispin
      real(rp) :: r, radial_scale, small_scale, energy, tmc
      real(rp) :: potential(size(radius), 2)

      if (size(radius) < 3 .or. mod(size(radius), 2) == 0 .or. lmax < 0 .or. lmax > 2) then
         error stop 'build_lr_test_radial_basis: invalid fixture dimensions'
      end if
      if (.not. all(ieee_is_finite(radius)) .or. any(radius(2:) <= 0.0_rp)) then
         error stop 'build_lr_test_radial_basis: invalid radial mesh'
      end if

      call basis%initialize(size(radius), lmax, 2)
      basis%rofi = radius
      basis%mesh_a = mesh_a
      basis%mesh_b = mesh_b
      basis%nuclear_z = nuclear_z
      basis%channel_present = .true.

      do ispin = 1, 2
         do ir = 1, size(radius)
            r = radius(ir)
            potential(ir, ispin) = 0.012_rp*real(ispin, rp) + 0.0017_rp*r + &
               0.00011_rp*real(ispin*ispin, rp)*r**2
         end do
         basis%potential(:, ispin) = potential(:, ispin)
         do l = 0, lmax
            energy = -0.31_rp + 0.043_rp*real(l, rp) + 0.017_rp*real(ispin - 1, rp)
            basis%enu_radial(l + 1, ispin) = energy
            basis%enu_work(l + 1, ispin) = energy
            radial_scale = 0.46_rp + 0.075_rp*real(l, rp) + 0.021_rp*real(ispin - 1, rp)
            small_scale = 0.016_rp + 0.003_rp*real(l, rp) + 0.001_rp*real(ispin, rp)
            do ir = 1, size(radius)
               r = radius(ir)
               basis%phi_large(ir, l + 1, ispin) = (0.71_rp + 0.037_rp*real(l, rp) + &
                  0.014_rp*real(ispin, rp))*r**l*(1.0_rp + 0.11_rp*r + 0.013_rp*r**2)* &
                  exp(-radial_scale*r)
               basis%phidot_large(ir, l + 1, ispin) = (0.19_rp + 0.023_rp*real(l, rp) + &
                  0.009_rp*real(ispin, rp))*r**l*(1.0_rp + 0.17_rp*r + 0.021_rp*r**2)* &
                  exp(-(0.63_rp*radial_scale + 0.025_rp)*r)
               basis%phiddot_large(ir, l + 1, ispin) = (0.061_rp + 0.011_rp*real(l, rp) + &
                  0.004_rp*real(ispin, rp))*r**l*(1.0_rp + 0.23_rp*r + 0.019_rp*r**2)* &
                  exp(-(0.39_rp*radial_scale + 0.041_rp)*r)

               basis%phi_small(ir, l + 1, ispin) = small_scale*r**(l + 1)* &
                  (1.0_rp + 0.07_rp*r)*exp(-(1.13_rp*radial_scale + 0.02_rp)*r)
               basis%phidot_small(ir, l + 1, ispin) = (0.42_rp*small_scale)*r**(l + 1)* &
                  (1.0_rp + 0.09_rp*r + 0.01_rp*r**2)*exp(-(0.81_rp*radial_scale + 0.03_rp)*r)
               basis%phiddot_small(ir, l + 1, ispin) = (0.17_rp*small_scale)*r**(l + 1)* &
                  (1.0_rp + 0.12_rp*r + 0.008_rp*r**2)*exp(-(0.57_rp*radial_scale + 0.05_rp)*r)

               if (ir == 1) then
                  basis%tmc(ir, l + 1, ispin) = scalar_relativistic_c
                  basis%gfac(ir, l + 1, ispin) = 1.0_rp
               else
                  tmc = scalar_relativistic_c - (potential(ir, ispin) - 2.0_rp*nuclear_z/r - energy) / &
                     scalar_relativistic_c
                  basis%tmc(ir, l + 1, ispin) = tmc
                  basis%gfac(ir, l + 1, ispin) = 1.0_rp + real(l*(l + 1), rp)/(tmc*r)**2
               end if
            end do
         end do
      end do

      if (.not. all(ieee_is_finite(basis%rofi)) .or. .not. all(ieee_is_finite(basis%potential)) .or. &
          .not. all(ieee_is_finite(basis%tmc)) .or. .not. all(ieee_is_finite(basis%gfac)) .or. &
          .not. all(ieee_is_finite(basis%phi_large)) .or. .not. all(ieee_is_finite(basis%phi_small)) .or. &
          .not. all(ieee_is_finite(basis%phidot_large)) .or. .not. all(ieee_is_finite(basis%phidot_small)) .or. &
          .not. all(ieee_is_finite(basis%phiddot_large)) .or. .not. all(ieee_is_finite(basis%phiddot_small))) then
         error stop 'build_lr_test_radial_basis: non-finite fixture data'
      end if
   end subroutine build_lr_test_radial_basis


   subroutine assert_lr_test_product_branch_norms(basis, space, label)
      type(lmto_radial_basis), intent(in) :: basis
      type(response_space_layout), intent(in) :: space
      character(len=*), intent(in), optional :: label
      integer, parameter :: branches(6) = [lmto_product_branch_00, lmto_product_branch_10, &
         lmto_product_branch_01, lmto_product_branch_11, lmto_product_branch_20, lmto_product_branch_02]
      real(rp) :: candidate(space%npoint), norm
      integer :: branch, ir, l, lp
      character(len=64) :: prefix

      if (space%npoint /= basis%npoint .or. size(space%radial_weights) /= basis%npoint) then
         error stop 'assert_lr_test_product_branch_norms: radial layout mismatch'
      end if
      l = min(2, basis%lmax)
      lp = l
      prefix = 'LR fixture'
      if (present(label)) prefix = label
      do branch = 1, size(branches)
         do ir = 1, basis%npoint
            call lmto_product_second_order_radial_branch(basis, ir, l, lp, 1, 2, branches(branch), candidate(ir))
         end do
         norm = sqrt(sum(space%radial_weights*candidate**2))
         write (*, '(a,1x,a,i0,a,es16.8)') trim(prefix), 'branch=', branches(branch), ' weighted_norm=', norm
         if (.not. ieee_is_finite(norm) .or. norm <= tiny(1.0_rp)) then
            error stop 'assert_lr_test_product_branch_norms: invalid second-order candidate norm'
         end if
      end do
   end subroutine assert_lr_test_product_branch_norms

end module lr_test_radial_fixture_mod
