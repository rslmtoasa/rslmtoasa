! LR-01 snapshot algebra and independent logarithmic-mesh quadrature supplement
program test_lr_radial_ground_state
   use precision_mod, only: rp
   use radial_ground_state_mod, only: radial_ground_state, RADIAL_PI
   implicit none

   integer, parameter :: nr = 101
   real(rp), parameter :: a = 0.03_rp, b = 0.10_rp
   real(rp) :: r(nr), rho(nr, 2), nup(nr), ndn(nr), vup(nr), vdn(nr), total_v(nr, 2)
   real(rp) :: origin(2), independent_up, independent_down, bfield
   real(rp) :: max_decomp, max_field_invariance
   real(rp) :: w, drdi
   integer :: ir
   type(radial_ground_state) :: snapshot, field_snapshot

   r(1) = 0.0_rp
   do ir = 2, nr
      r(ir) = b*(exp(a*real(ir - 1, rp)) - 1.0_rp)
   end do
   do ir = 1, nr
      nup(ir) = 0.70_rp*exp(-r(ir))
      ndn(ir) = 0.40_rp*exp(-r(ir))
      rho(ir, 1) = 4.0_rp*RADIAL_PI*r(ir)**2*nup(ir)
      rho(ir, 2) = 4.0_rp*RADIAL_PI*r(ir)**2*ndn(ir)
      vup(ir) = -1.3_rp*nup(ir)
      vdn(ir) = -0.7_rp*ndn(ir)
      total_v(ir, 1) = 2.0_rp + vup(ir) + 0.125_rp
      total_v(ir, 2) = 2.0_rp + vdn(ir) - 0.125_rp
   end do
   origin = [nup(1), ndn(1)]
   bfield = 0.125_rp

   call snapshot%capture(a, b, r, rho, origin, vup, vdn, total_v, bfield)
   call snapshot%mark_accepted(7, 1.0e-8_rp)

   independent_up = 0.0_rp
   independent_down = 0.0_rp
   do ir = 1, nr
      w = 2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp
      if (ir == 1 .or. ir == nr) w = 1.0_rp/3.0_rp
      drdi = a*(r(ir) + b)
      independent_up = independent_up + w*drdi*(4.0_rp*RADIAL_PI*r(ir)**2*nup(ir))
      independent_down = independent_down + w*drdi*(4.0_rp*RADIAL_PI*r(ir)**2*ndn(ir))
   end do

   if (abs(snapshot%charge_integral(1) - independent_up) > 2.0e-13_rp) error stop 'charge up quadrature failed'
   if (abs(snapshot%charge_integral(2) - independent_down) > 2.0e-13_rp) error stop 'charge down quadrature failed'
   if (abs(snapshot%spin_number_integral() - (independent_up - independent_down)) > 2.0e-13_rp) &
      error stop 'spin quadrature failed'
   if (.not. snapshot%valid .or. .not. snapshot%accepted) error stop 'snapshot lifecycle failed'
   if (snapshot%accepted_inner_iteration /= 7) error stop 'accepted iteration was not retained'
   if (abs(snapshot%constraining_field_ry - bfield) > 2.0e-14_rp) error stop 'field metadata failed'

   max_decomp = maxval(abs(snapshot%vxc_up - (snapshot%vxc_scalar + snapshot%bxc_pauli)))
   max_decomp = max(max_decomp, maxval(abs(snapshot%vxc_down - (snapshot%vxc_scalar - snapshot%bxc_pauli))))
   if (max_decomp > 2.0e-14_rp) error stop 'Pauli decomposition failed'
   if (maxval(abs(snapshot%total_ks_spin_splitting - (snapshot%delta_vxc + 2.0_rp*bfield))) > 2.0e-14_rp) &
      error stop 'total/XC field separation failed'

   ! A second total potential with a different constraining field must leave
   ! the captured XC arrays unchanged while shifting only the total splitting.
   total_v(:, 1) = 2.0_rp + vup + 0.375_rp
   total_v(:, 2) = 2.0_rp + vdn - 0.375_rp
   call field_snapshot%capture(a, b, r, rho, origin, vup, vdn, total_v, 0.375_rp)
   max_field_invariance = maxval(abs(field_snapshot%vxc_up - snapshot%vxc_up))
   max_field_invariance = max(max_field_invariance, maxval(abs(field_snapshot%vxc_down - snapshot%vxc_down)))
   if (max_field_invariance > 2.0e-14_rp) error stop 'constraining field contaminated XC snapshot'
   if (maxval(abs(field_snapshot%total_ks_spin_splitting - snapshot%total_ks_spin_splitting - 0.5_rp)) > 2.0e-14_rp) &
      error stop 'constraining field shift failed'

   write (*, '(a,2(es16.8,1x))') 'LR-01 independent N(up/down) = ', independent_up, independent_down
   write (*, '(a,es16.8)') 'LR-01 independent spin number = ', independent_up - independent_down
   write (*, '(a,es16.8)') 'LR-01 Pauli decomposition max error = ', max_decomp
   write (*, '(a)') 'UnitLrRadialGroundState: PASS (algebraic supplement)'
end program test_lr_radial_ground_state
