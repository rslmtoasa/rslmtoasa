!------------------------------------------------------------------------------
! LR-BASIS-00 radial normalization oracle
!------------------------------------------------------------------------------
program test_lr_basis_radial
   use precision_mod, only: rp
   use self_mod, only: legacy_radial_fixture
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   implicit none

   integer, parameter :: nr = 101
   real(rp), parameter :: a = 0.03_rp, b = 0.10_rp, z = 1.0_rp
   real(rp), parameter :: norm_tol = 2.0e-10_rp
   real(rp), parameter :: derivative_tol = 5.0e-5_rp
   real(rp) :: rofi(nr), potential(nr), energy
   real(rp), allocatable :: g(:, :), gp(:, :), gpp(:, :)
   real(rp), allocatable :: gpack(:), gdotpack(:), gddotpack(:)
   type(lmto_radial_basis) :: radial
   real(rp) :: norm_phi, cross, norm_dot, second_identity
   integer :: ir

   rofi(1) = 0.0_rp
   do ir = 2, nr
      rofi(ir) = b*(exp(a*real(ir - 1, rp)) - 1.0_rp)
   end do
   potential = 0.0_rp

   ! This is a real radial solution produced by the production solver, not a
   ! synthetic array satisfying a normalization equation by construction.
   call legacy_radial_fixture(z, 0, a, b, rofi, potential, energy, g, gp, gpp)
   gpack = reshape(g, [2*nr])
   gdotpack = reshape(gp, [2*nr])
   gddotpack = reshape(gpp, [2*nr])
   call radial%initialize(nr, 0, 1)
   call radial%capture_channel(0, 1, energy, rofi, potential, a, b, z, gpack, gdotpack, gddotpack, energy)

   norm_phi = radial%normalization_integral(0, 1, 'phi')
   cross = radial%overlap_integral(0, 1, 'phidot')
   norm_dot = radial%normalization_integral(0, 1, 'phidot')
   second_identity = norm_dot + radial%overlap_integral(0, 1, 'phiddot')

   write (*, '(a,es16.8)') 'LR-BASIS radial energy = ', energy
   write (*, '(a,es16.8)') 'RSEQSR/GINTSR <phi|phi>-1 = ', norm_phi - 1.0_rp
   write (*, '(a,es16.8)') 'Re <phi|phidot> = ', cross
   write (*, '(a,es16.8)') '<phidot|phidot> + Re <phi|phiddot> = ', second_identity

   if (abs(norm_phi - 1.0_rp) > norm_tol) error stop 'LR-BASIS radial normalization failed'
   if (abs(cross) > derivative_tol) error stop 'LR-BASIS first normalization derivative failed'
   if (abs(second_identity) > derivative_tol) error stop 'LR-BASIS second normalization derivative failed'
   write (*, '(a)') 'UnitLrBasisRadial: PASS'
end program test_lr_basis_radial
