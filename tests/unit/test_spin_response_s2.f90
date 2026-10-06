!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!> @brief spin_response S2 on bcc Fe: the moment-reversed state gives the same U, tr chi0 and tr chi, and spin_swapped flips.
!> @details Usage, from a scratch copy of tests/spin_response/fe_bcc in which `rslmto.x input_driver.nml`
!>          has run (see fe_scratch.py): test_spin_response_s2 <oracle directory>. The reversed state swaps the spin
!>          index of every (l, spin) potential parameter (reverse_spin), so it is the normal state with the spin
!>          labels exchanged and needs no second SCF.
!------------------------------------------------------------------------------
program test_spin_response_s2
   use precision_mod, only: rp
   use control_mod, only: control
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use energy_mod, only: energy
   use hamiltonian_mod, only: hamiltonian
   use reciprocal_mod, only: reciprocal
   use spin_response_mod, only: spin_response, solve_dyson
   use spin_response_fe_state
   implicit none

   integer, parameter :: nw = 41
   type(control), target :: control_obj, control_rev
   type(lattice), target :: lattice_obj, lattice_rev
   type(charge), target :: charge_obj, charge_rev
   type(energy), target :: energy_obj, energy_rev
   type(hamiltonian), target :: ham, ham_rev
   type(reciprocal), target :: rec, rec_rev
   type(spin_response) :: sr, sr_rev
   character(len=512) :: dir
   real(rp) :: omega(nw), q(3)
   complex(rp) :: chi0(1, 1, nw), chi(1, 1, nw), chi0_rev(1, 1, nw), chi_rev(1, 1, nw)
   integer :: iw

   call get_command_argument(1, dir)
   if (len_trim(dir) == 0) error stop 'usage: test_spin_response_s2 <oracle directory>'
   call read_tolerances(trim(dir))
   call require('s2_rel', s2_rel)

   call build_fe_state(.false., control_obj, lattice_obj, charge_obj, energy_obj, ham, 'input_driver.nml')
   rec = reciprocal(ham)
   rec%fermi_level = energy_obj%fermi
   sr = spin_response(rec)
   call sr%prepare()
   call sr%interaction()
   call build_fe_state(.true., control_rev, lattice_rev, charge_rev, energy_rev, ham_rev, 'input_driver.nml')
   rec_rev = reciprocal(ham_rev)
   rec_rev%fermi_level = energy_rev%fermi
   sr_rev = spin_response(rec_rev)
   call sr_rev%prepare()
   call sr_rev%interaction()

   ! Mesh-commensurate q along Gamma-H (xi = 1/3), 0 to 0.02 Ry, eta = 5e-3 Ry, Juelich amplitudes.
   q = [-1.0_rp, 1.0_rp, 1.0_rp]/6.0_rp
   omega = [(0.02_rp*(iw - 1)/real(nw - 1, rp), iw=1, nw)]
   call sr%bare_response(q, omega, 5.0e-3_rp, .true., chi0)
   call solve_dyson(chi0, sr%uj, chi)
   call sr_rev%bare_response(q, omega, 5.0e-3_rp, .true., chi0_rev)
   call solve_dyson(chi0_rev, sr_rev%uj, chi_rev)

   ! S2: catches a reversed state used without the spin relabelling (negative M, spin blocks exchanged).
   call check('S2 spin_swapped: normal state F, reversed state T', (.not. sr%spin_swapped) .and. sr_rev%spin_swapped)
   call report('S2 U_Juelich reversed vs normal, relative', abs(sr_rev%uj(1) - sr%uj(1))/abs(sr%uj(1)), s2_rel)
   call report('S2 tr chi0 reversed vs normal, max|a-b|/max|b|', &
               maxval(abs(chi0_rev(1, 1, :) - chi0(1, 1, :)))/maxval(abs(chi0(1, 1, :))), s2_rel)
   call report('S2 tr chi reversed vs normal, max|a-b|/max|b|', &
               maxval(abs(chi_rev(1, 1, :) - chi(1, 1, :)))/maxval(abs(chi(1, 1, :))), s2_rel)
   write (*, '(a,f10.7,a,f10.7,a)') 'S2 U_Juelich normal ', sr%uj(1), ' Ry, reversed ', sr_rev%uj(1), ' Ry'
   call check('S2 non-vacuity: chi differs from chi0 by more than 1e-3 relative', &
              maxval(abs(chi(1, 1, :) - chi0(1, 1, :)))/maxval(abs(chi0(1, 1, :))) > 1.0e-3_rp)
   call finish()

end program test_spin_response_s2
