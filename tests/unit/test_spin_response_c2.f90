!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!> @brief spin_response C2 (d amplitudes and moments) on the accepted bcc-Fe state, plus the driver.
!> @details Usage, from a scratch copy of tests/spin_response/fe_bcc in which the executable has
!>          already run post_processing='spin_response': test_spin_response_c2 <oracle directory>.
!>          Tolerances come from <dir>/tolerances.nml.
!------------------------------------------------------------------------------
program test_spin_response_c2
   use precision_mod, only: rp
   use control_mod, only: control
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use energy_mod, only: energy
   use hamiltonian_mod, only: hamiltonian
   use reciprocal_mod, only: reciprocal
   use spin_response_mod, only: spin_response, d_amplitudes_coefficient
   use spin_response_fe_state
   implicit none

   type(control), target :: control_obj, control_rev
   type(lattice), target :: lattice_obj, lattice_rev
   type(charge), target :: charge_obj, charge_rev
   type(energy), target :: energy_obj, energy_rev
   type(hamiltonian), target :: ham, ham_rev
   type(reciprocal), target :: rec, rec_ref, rec_rev
   type(spin_response) :: sr, sr_rev
   character(len=512) :: dir
   real(rp), allocatable :: e(:, :), bm(:, :, :, :)
   complex(rp), allocatable :: v(:, :, :), a(:, :, :, :, :)
   real(rp) :: enu(2), p(2), x, q0, cen, w, n_juelich(2), dmax
   integer :: nk, first, last, s

   call get_command_argument(1, dir)
   if (len_trim(dir) == 0) error stop 'usage: test_spin_response_c2 <oracle directory>'
   call read_tolerances(trim(dir))
   call build_fe_state(.false., control_obj, lattice_obj, charge_obj, energy_obj, ham)

   rec = reciprocal(ham)
   rec%fermi_level = energy_obj%fermi
   sr = spin_response(rec)
   call sr%prepare()
   nk = rec%k_workset%nk_local

   ! C2.1: sum_n |v|^2 = 1 for every d orbital, site, spin and k.
   ! Catches the amplitude kernel reading outside an orbital block or one orbital twice (a pure
   ! shift of the orbital range is invisible here: the eigenvector matrix is unitary).
   dmax = 0.0_rp
   do first = 1, nk, 64
      last = min(nk, first + 63)
      call sr%eigenpairs_chunk([0.0_rp, 0.0_rp, 0.0_rp], first, last, e, v)
      call d_amplitudes_coefficient(v, 1, a)
      dmax = max(dmax, maxval(abs(sum(abs(a(:, :, :, 1, :))**2, dim=2) - 1.0_rp)))
   end do
   call report('C2.1 completeness, max |sum_n |v|^2 - 1|', dmax, eigen_contract_abs)

   ! Reference band moments from the pre-campaign path, on E_F from reciprocal's own Fermi solve on this mesh.
   call reference_band_moments(ham, rec_ref, bm)
   call report('C2.2 E_F used vs reciprocal Fermi solve (Ry)', abs(sr%fermi_used - rec_ref%fermi_level), eigen_contract_abs)

   ! C2.2: Mills d occupation per spin vs band_moments(site, d, spin, 1).
   ! Catches a shifted d range, wrong E_F or temperature, wrong k-weight normalization.
   do s = 1, 2
      call report('C2.2 Mills d occupation, spin '//achar(48 + s)//', |n - q0|/q0', &
                  abs(sr%d_occ_mills(1, s) - bm(1, 3, s, 1))/bm(1, 3, s, 1), moment_rel)
   end do

   ! C2.3: Juelich d occupation vs [m0 + 2 D p m1 + D^2 p^2 m2]/N, moments about E_nu from q0, <E>, w.
   ! Catches ppar used as p, E_nu of the other spin, and a wrong E_F in the Juelich factor.
   do s = 1, 2
      enu(s) = lattice_obj%symbolic_atoms(lattice_obj%nbulk + 1)%potential%enu(2, s)
      p(s) = 1.0_rp/lattice_obj%symbolic_atoms(lattice_obj%nbulk + 1)%potential%ppar(2, s)**2
      q0 = bm(1, 3, s, 1)
      cen = bm(1, 3, s, 2) - enu(s)
      w = bm(1, 3, s, 3)
      x = rec_ref%fermi_level - enu(s)
      n_juelich(s) = (q0 + 2.0_rp*x*p(s)*q0*cen + x*x*p(s)**2*q0*(w*w + cen*cen))/(1.0_rp + x*x*p(s))
      call report('C2.3 Juelich d occupation, spin '//achar(48 + s)//', |n - ref|/ref', &
                  abs(sr%d_occ_juelich(1, s) - n_juelich(s))/n_juelich(s), moment_rel)
   end do

   ! C2.4: reversed state (spin index of every potential parameter swapped): same |M|, spin_swapped flips.
   ! Catches a missing or wrongly applied spin relabelling.
   call build_fe_state(.true., control_rev, lattice_rev, charge_rev, energy_rev, ham_rev)
   rec_rev = reciprocal(ham_rev)
   rec_rev%fermi_level = energy_rev%fermi
   sr_rev = spin_response(rec_rev)
   call sr_rev%prepare()
   call check('C2.4 spin_swapped: normal state F, reversed state T', (.not. sr%spin_swapped) .and. sr_rev%spin_swapped)
   call report('C2.4 |M_juelich| reversed vs normal, relative', abs(sr_rev%moment(1) - sr%moment(1))/abs(sr%moment(1)), moment_rel)
   call report('C2.4 |M_mills| reversed vs normal, relative', &
               abs(sum(sr_rev%d_occ_mills(1, :)*[1.0_rp, -1.0_rp]) - sum(sr%d_occ_mills(1, :)*[1.0_rp, -1.0_rp])) &
               /abs(sum(sr%d_occ_mills(1, :)*[1.0_rp, -1.0_rp])), moment_rel)
   call report('C2.4 E_F used reversed vs normal (Ry)', abs(sr_rev%fermi_used - sr%fermi_used), eigen_contract_abs)

   ! Driver: the executable's <prefix>_state.dat vs the in-process values.
   call compare_driver_state()

   call finish()

contains

   !> spin_response_state.dat written by `rslmto.x` with post_processing='spin_response' in this directory.
   subroutine compare_driver_state()
      integer :: u, ios
      character(len=64) :: key
      character(len=256) :: line
      real(rp) :: ef, n_el, m_mills, m_juelich, row(4)
      logical :: seen_ef, seen_n, seen_m

      seen_ef = .false.
      seen_n = .false.
      seen_m = .false.
      open (newunit=u, file='spin_response_state.dat', status='old', action='read')
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (line(1:1) == '#') cycle
         read (line, *, iostat=ios) key
         select case (trim(key))
         case ('fermi_used')
            read (line(len('fermi_used') + 1:), *) ef
            seen_ef = .true.
         case ('electron_count')
            read (line(len('electron_count') + 1:), *) n_el
            seen_n = .true.
         case default
            if (.not. seen_m .and. key(1:1) >= '0' .and. key(1:1) <= '9') then
               read (line, *) s, m_mills, m_juelich, row
               seen_m = .true.
            end if
         end select
      end do
      close (u)
      call check('driver: state file has E_F, N and a site row', seen_ef .and. seen_n .and. seen_m)
      call report('driver: E_F used vs in-process (Ry)', abs(ef - sr%fermi_used), eigen_contract_abs)
      call report('driver: N vs in-process (electrons)', abs(n_el - sr%electron_count), electron_count_abs)
      call report('driver: M_mills vs in-process, relative', &
                  abs(m_mills - (sr%d_occ_mills(1, 1) - sr%d_occ_mills(1, 2)))/abs(m_mills), moment_rel)
      call report('driver: M_juelich vs in-process, relative', abs(m_juelich - sr%moment(1))/abs(m_juelich), moment_rel)
   end subroutine compare_driver_state

end program test_spin_response_c2
