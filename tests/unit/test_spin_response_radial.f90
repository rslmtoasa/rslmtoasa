!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!> @brief ALSDA radial producer on the bcc-Fe state of the deck (spin_response step 3d).
!> @details Usage, from a scratch copy of tests/spin_response/fe_bcc: test_spin_response_radial <oracle directory>.
!>          The sphere is solved again by site_radial_data; its n_up - n_dn is integrated over the sphere and compared
!>          with sum_l (Q0_up - Q0_dn) of the deck potential (the moment of the SCF that wrote it), and its E_nu with
!>          the E_nu stored in the deck. Tolerances come from <dir>/tolerances.nml. Then the unscaled ALSDA kernel goes
!>          through the response (alsda_kernel, prepare, interaction): lambda must be > 0 and the q = 0 Goldstone
!>          Rayleigh eigenvalue of 1 + chi0 U along M must have magnitude < 1 (bounds fixed before the first measurement).
!------------------------------------------------------------------------------
program test_spin_response_radial
   use precision_mod, only: rp
   use math_mod, only: pi
   use control_mod, only: control
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use energy_mod, only: energy
   use hamiltonian_mod, only: hamiltonian
   use reciprocal_mod, only: reciprocal
   use spin_response_mod, only: spin_response
   use alsda_radial_mod, only: site_radial_data, alsda_kernel
   use spin_response_fe_state
   implicit none

   type(control), target :: control_obj
   type(lattice), target :: lattice_obj
   type(charge), target :: charge_obj
   type(energy), target :: energy_obj
   type(hamiltonian), target :: ham
   type(reciprocal), target :: rec
   type(spin_response) :: sr
   character(len=512) :: dir
   real(rp), allocatable :: r(:), wq(:), m(:), dv(:), enu(:, :), g_up(:, :), g_dn(:, :), s_up(:, :), s_dn(:, :), gfac_up(:, :), gfac_dn(:, :)
   real(rp) :: moment

   call get_command_argument(1, dir)
   if (len_trim(dir) == 0) error stop 'usage: test_spin_response_radial <oracle directory>'
   call read_tolerances(trim(dir))
   call require('radial_moment_rel', radial_moment_rel)
   call require('radial_enu_abs', radial_enu_abs)
   call build_fe_state(.false., control_obj, lattice_obj, charge_obj, energy_obj, ham)

   call site_radial_data(lattice_obj, 1, [0, 1, 2], 0.0_rp, .false., r, wq, m, dv, enu, g_up, g_dn, s_up, s_dn, gfac_up, gfac_dn)
   associate (pot => lattice_obj%symbolic_atoms(1)%potential)
      moment = sum(pot%ql(1, 0:2, 1) - pot%ql(1, 0:2, 2))
      ! catches: wrong weights or mesh; blind to the core (9.000 electrons per spin, net moment zero)
      call report('radial: sphere integral of n_up - n_dn vs sum_l (Q0_up - Q0_dn)', &
                  abs(sum(wq*4.0_rp*pi*r**2*m) - moment)/abs(moment), radial_moment_rel)
      ! catches: l or spin of the radial arrays mixed up
      call report('radial: E_nu of the re-run vs the deck, max over l and spin (Ry)', maxval(abs(enu - pot%enu(0:2, :))), radial_enu_abs)
   end associate

   rec = reciprocal(ham)
   rec%fermi_level = energy_obj%fermi
   sr = spin_response(rec)
   sr%kernel = 'alsda'
   call sr%prepare()
   call alsda_kernel(sr, lattice_obj)
   call sr%interaction()
   write (*, '(a,es12.4,a,es12.4)') 'ALSDA kernel on Fe: lambda ', sr%kernel_lambda, ', Goldstone eigenvalue ', sr%kernel_goldstone_eig
   ! catches: kernel sign wrong in U = -K (lambda < 0, eigenvalue 1 + 1/lambda)
   call check('ALSDA kernel on Fe: lambda > 0', sr%kernel_lambda > 0.0_rp)
   call check('ALSDA kernel on Fe: |Goldstone eigenvalue| < 1', abs(sr%kernel_goldstone_eig) < 1.0_rp)
   call finish()
end program test_spin_response_radial
