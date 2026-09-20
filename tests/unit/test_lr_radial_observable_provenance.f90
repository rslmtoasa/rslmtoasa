! DRESP-09W synthetic radial observable provenance fixture.
!
! The expected values below are formed directly from the synthetic radial
! functions and a single occupied state.  No production reconstruction is
! called to manufacture the reference values.
program test_lr_radial_observable_provenance

   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_radial_observable_provenance_mod, only: pauli_density_from_moments, sr_density_from_moments, &
      weighted_from_density, physical_from_weighted
   implicit none

   integer, parameter :: npoint = 7
   real(rp), parameter :: tol = 2.0e-12_rp
   type(lmto_radial_basis) :: radial(1)
   real(rp) :: moments(3, 1, 1, 2), pauli(1, npoint), sr(1, npoint)
   real(rp) :: weighted(1, npoint), recovered(1, npoint), direct_sr(npoint), direct_pauli(npoint)
   real(rp) :: expected_difference(npoint), delta, energy, coefficient
   integer :: ir
   logical :: failed

   failed = .false.
   call radial(1)%initialize(npoint, 0, 2)
   radial(1)%rofi = [0.0_rp, 0.10_rp, 0.25_rp, 0.45_rp, 0.70_rp, 1.00_rp, 1.35_rp]
   radial(1)%mesh_a = 0.1_rp
   radial(1)%mesh_b = 0.2_rp
   radial(1)%enu_work(:, :) = reshape([0.11_rp, 0.17_rp], [1, 2])
   radial(1)%enu_radial = radial(1)%enu_work
   radial(1)%channel_present = .true.

   do ir = 1, npoint
      radial(1)%phi_large(ir, 1, 1) = 0.30_rp + 0.07_rp*ir
      radial(1)%phi_small(ir, 1, 1) = 0.04_rp + 0.01_rp*ir
      radial(1)%phidot_large(ir, 1, 1) = -0.08_rp + 0.015_rp*ir
      radial(1)%phidot_small(ir, 1, 1) = 0.02_rp - 0.003_rp*ir
      radial(1)%phiddot_large(ir, 1, 1) = 0.11_rp - 0.004_rp*ir
      radial(1)%phiddot_small(ir, 1, 1) = -0.05_rp + 0.006_rp*ir
      radial(1)%gfac(ir, 1, 1) = 1.20_rp + 0.01_rp*ir
      radial(1)%tmc(ir, 1, 1) = 250.0_rp + ir
      ! Keep the down channel nonzero and distinct so spin-channel shape and
      ! sign handling are exercised, while its occupation remains zero.
      radial(1)%phi_large(ir, 1, 2) = 0.21_rp + 0.03_rp*ir
      radial(1)%phi_small(ir, 1, 2) = 0.03_rp + 0.005_rp*ir
      radial(1)%phidot_large(ir, 1, 2) = 0.06_rp - 0.007_rp*ir
      radial(1)%phidot_small(ir, 1, 2) = -0.01_rp + 0.002_rp*ir
      radial(1)%phiddot_large(ir, 1, 2) = -0.09_rp + 0.003_rp*ir
      radial(1)%phiddot_small(ir, 1, 2) = 0.02_rp - 0.001_rp*ir
      radial(1)%gfac(ir, 1, 2) = 1.35_rp + 0.008_rp*ir
      radial(1)%tmc(ir, 1, 2) = 260.0_rp + ir
   end do

   energy = 0.41_rp
   coefficient = 0.73_rp
   delta = energy - radial(1)%enu_work(1, 1)
   moments = 0.0_rp
   moments(:, 1, 1, 1) = coefficient*[1.0_rp, energy, energy**2]

   call pauli_density_from_moments(radial, moments, pauli)
   call sr_density_from_moments(radial, moments, sr)

   do ir = 1, npoint
      direct_pauli(ir) = coefficient*(radial(1)%phi_large(ir, 1, 1) + delta*radial(1)%phidot_large(ir, 1, 1))**2
      direct_sr(ir) = coefficient*((radial(1)%gfac(ir, 1, 1)*radial(1)%phi_large(ir, 1, 1)**2 + &
         radial(1)%phi_small(ir, 1, 1)**2) + 2.0_rp*delta*(radial(1)%gfac(ir, 1, 1)* &
         radial(1)%phi_large(ir, 1, 1)*radial(1)%phidot_large(ir, 1, 1) + radial(1)%phi_small(ir, 1, 1)* &
         radial(1)%phidot_small(ir, 1, 1)) + delta**2*(radial(1)%gfac(ir, 1, 1)* &
         (radial(1)%phidot_large(ir, 1, 1)**2 + radial(1)%phi_large(ir, 1, 1)*radial(1)%phiddot_large(ir, 1, 1)) + &
         radial(1)%phidot_small(ir, 1, 1)**2 + radial(1)%phi_small(ir, 1, 1)*radial(1)%phiddot_small(ir, 1, 1)))
      expected_difference(ir) = coefficient*delta**2*radial(1)%phi_large(ir, 1, 1)* &
         radial(1)%phiddot_large(ir, 1, 1)
   end do

   if (maxval(abs(sr(1, :) - direct_sr)) > tol) failed = .true.
   if (maxval(abs(pauli(1, :) - direct_pauli - expected_difference)) > tol) failed = .true.
   if (maxval(abs(pauli(1, :) - direct_pauli)) <= 1.0e-5_rp) failed = .true.

   call weighted_from_density(radial(1)%rofi, pauli, weighted)
   call physical_from_weighted(radial(1)%rofi, weighted, recovered)
   if (weighted(1, 1) /= 0.0_rp) failed = .true.
   if (maxval(abs(recovered(1, 2:npoint) - pauli(1, 2:npoint))) > tol) failed = .true.

   if (failed) then
      write (*, '(a)') 'UnitLrRadialObservableProvenance: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrRadialObservableProvenance: PASS (P/SR direct oracles, phi*phiddot, origin/measure)'

end program test_lr_radial_observable_provenance
