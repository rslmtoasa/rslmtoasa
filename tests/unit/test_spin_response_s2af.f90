!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!> @brief spin_response S2' on a two-atom bcc-Fe cell with antiparallel moments (frozen, not self-consistent):
!>        the state (up, down) and the state (down, up) give the same gauge-invariant outputs.
!> @details Usage, from a scratch copy made by af_scratch.py: test_spin_response_s2af <oracle directory>.
!>          (down, up) is (up, down) with the spin labels exchanged and also with the sublattices exchanged, so
!>          when the two runs choose different spin_swapped the comparison is the S2 relabelling, and when they
!>          choose the same one it is the sublattice-exchange symmetry. The total moment is zero up to round-off,
!>          so spin_swapped is decided by round-off; both flags are printed.
!------------------------------------------------------------------------------
program test_spin_response_s2af
   use precision_mod, only: rp
   use control_mod, only: control
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use energy_mod, only: energy
   use hamiltonian_mod, only: hamiltonian
   use reciprocal_mod, only: reciprocal
   use spin_response_mod, only: spin_response, solve_dyson, adiabatic_frequencies
   use spin_response_fe_state
   implicit none

   integer, parameter :: nw = 41
   type(control), target :: control_x, control_y
   type(lattice), target :: lattice_x, lattice_y
   type(charge), target :: charge_x, charge_y
   type(energy), target :: energy_x, energy_y
   type(hamiltonian), target :: ham_x, ham_y
   type(reciprocal), target :: rec_x, rec_y
   type(spin_response) :: sr_x, sr_y
   character(len=512) :: dir
   real(rp) :: omega(nw), q(3, 1), wad_x(2), wad_y(2), imag_max
   complex(rp) :: chi0_x(2, 2, nw, 1), chi_x(2, 2, nw), chi0_y(2, 2, nw, 1), chi_y(2, 2, nw), st_x(2, 2, 1), st_y(2, 2, 1)
   complex(rp) :: tr0_x(nw), tr0_y(nw), tr_x(nw), tr_y(nw)
   integer :: iw

   call get_command_argument(1, dir)
   if (len_trim(dir) == 0) error stop 'usage: test_spin_response_s2af <oracle directory>'
   call read_tolerances(trim(dir))
   call require('s2_rel', s2_rel)

   call build_fe_state(.true., control_x, lattice_x, charge_x, energy_x, ham_x, 'input.nml', site=2)
   rec_x = reciprocal(ham_x)
   rec_x%fermi_level = energy_x%fermi
   sr_x = spin_response(rec_x)
   call sr_x%prepare()
   call sr_x%interaction()
   call build_fe_state(.true., control_y, lattice_y, charge_y, energy_y, ham_y, 'input.nml', site=1)
   rec_y = reciprocal(ham_y)
   rec_y%fermi_level = energy_y%fermi
   sr_y = spin_response(rec_y)
   call sr_y%prepare()
   call sr_y%interaction()

   q(:, 1) = [0.0_rp, 0.0_rp, 1.0_rp]/6.0_rp
   omega = [(0.03_rp*(iw - 1)/real(nw - 1, rp), iw=1, nw)]
   call sr_x%bare_response_list(q, omega, 3.0e-3_rp, .true., chi0_x, st_x)
   call sr_y%bare_response_list(q, omega, 3.0e-3_rp, .true., chi0_y, st_y)
   call solve_dyson(chi0_x(:, :, :, 1), sr_x%uj, chi_x)
   call solve_dyson(chi0_y(:, :, :, 1), sr_y%uj, chi_y)
   do iw = 1, nw
      tr0_x(iw) = chi0_x(1, 1, iw, 1) + chi0_x(2, 2, iw, 1)
      tr0_y(iw) = chi0_y(1, 1, iw, 1) + chi0_y(2, 2, iw, 1)
      tr_x(iw) = chi_x(1, 1, iw) + chi_x(2, 2, iw)
      tr_y(iw) = chi_y(1, 1, iw) + chi_y(2, 2, iw)
   end do

   write (*, '(a,2l2,a,2f9.5,a,2f9.5)') 'spin_swapped (up,down) and (down,up):', sr_x%spin_swapped, sr_y%spin_swapped, &
      ' ; moments ', sr_x%moment, ' ; ', sr_y%moment
   write (*, '(a,2es12.4)') 'static |chi0_11|, |chi0_12| (1/Ry):', abs(st_x(1, 1, 1)), abs(st_x(1, 2, 1))
   ! non-vacuity: antiparallel moments, and a cross-sublattice chi0 that a decoupled (local-frame) vertex would remove
   call check('S2p non-vacuity: moments of opposite sign, each above 0.5 mu_B', &
              sr_x%moment(1)*sr_x%moment(2) < 0.0_rp .and. minval(abs(sr_x%moment)) > 0.5_rp)
   call check('S2p non-vacuity: |chi0_12|/|chi0_11| at the static q above 1e-2', &
              abs(st_x(1, 2, 1)) > 1.0e-2_rp*abs(st_x(1, 1, 1)))
   ! sublattice symmetry inside one state: both sites have the same U (cross-sublattice Goldstone solve)
   call report('S2p U_Juelich site 1 vs site 2, relative', abs(sr_x%uj(1) - sr_x%uj(2))/abs(sr_x%uj(1)), s2_rel)
   call report('S2p U_Juelich (up,down) vs (down,up), site 1, relative', abs(sr_x%uj(1) - sr_y%uj(1))/abs(sr_x%uj(1)), s2_rel)
   call report('S2p tr chi0, max|a-b|/max|b|', maxval(abs(tr0_x - tr0_y))/maxval(abs(tr0_y)), s2_rel)
   call report('S2p tr chi, max|a-b|/max|b|', maxval(abs(tr_x - tr_y))/maxval(abs(tr_y)), s2_rel)
   call check('S2p non-vacuity: chi differs from chi0 by more than 1e-3 relative', &
              maxval(abs(tr_x - tr0_x))/maxval(abs(tr0_x)) > 1.0e-3_rp)
   call adiabatic_frequencies(st_x(:, :, 1), sr_x%kmat, sr_x%occ_juelich(:, 1) - sr_x%occ_juelich(:, 2), wad_x, imag_max)
   write (*, '(a,2es14.6,a,es10.2)') 'omega_adiabatic (Ry) (up,down):', wad_x, ' ; max |Im| ', imag_max
   call adiabatic_frequencies(st_y(:, :, 1), sr_y%kmat, sr_y%occ_juelich(:, 1) - sr_y%occ_juelich(:, 2), wad_y, imag_max)
   call report('S2p omega_adiabatic (up,down) vs (down,up), max|a-b|/max|b|', maxval(abs(wad_x - wad_y))/maxval(abs(wad_y)), s2_rel)
   call report('S2p omega_adiabatic pair: |w1 + w2|/max|w| (sublattice symmetry)', abs(wad_x(1) + wad_x(2))/maxval(abs(wad_x)), s2_rel)
   call check('S2p non-vacuity: one adiabatic frequency of each sign', wad_x(1) < 0.0_rp .and. wad_x(2) > 0.0_rp)
   call finish()

end program test_spin_response_s2af
