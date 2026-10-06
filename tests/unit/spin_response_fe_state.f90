!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test helper
!> @brief Frozen bcc-Fe state for the spin_response C1/C2 tests, tolerances and reporting.
!> @details build_fe_state follows post_processing_spin_response and reads ./input.nml of the
!>          working directory (a scratch copy of tests/spin_response/fe_bcc). reverse_spin swaps the
!>          spin index of every (l, spin) potential parameter, giving the reversed state exactly.
!------------------------------------------------------------------------------
module spin_response_fe_state
   use precision_mod, only: rp
   use mpi_mod, only: g_parallel_context, parallel_context
   use timer_mod, only: g_timer, timer
   use control_mod, only: control
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use energy_mod, only: energy
   use hamiltonian_mod, only: hamiltonian
   use reciprocal_mod, only: reciprocal
   implicit none

   real(rp) :: electron_count_abs, eigen_contract_abs, moment_rel
   real(rp) :: s1_rel, s2_rel, e1_fe_rel, q0_pole_abs
   real(rp) :: region1_rel, region2_rel, baseline_repro_rel   ! read, used by the Stage 2 comparisons
   real(rp) :: a1_rel, a2_chi_rel, a2_identity_rel, static_eta0_rel, a3_rel, window_moment_rel, &
               dyson_u0_rel, pole_peak_rel, pole_crossing_rel
   namelist /spin_response_tolerances/ a1_rel, a2_chi_rel, a2_identity_rel, static_eta0_rel, a3_rel, &
      window_moment_rel, dyson_u0_rel, pole_peak_rel, pole_crossing_rel, &
      electron_count_abs, eigen_contract_abs, moment_rel, s1_rel, s2_rel, e1_fe_rel, q0_pole_abs, &
      region1_rel, region2_rel, baseline_repro_rel
   logical :: failed = .false.

contains

   subroutine read_tolerances(oracle_dir)
      character(*), intent(in) :: oracle_dir

      integer :: u

      g_parallel_context = parallel_context()
      g_timer = timer()
      a1_rel = -1.0_rp; a2_chi_rel = -1.0_rp; a2_identity_rel = -1.0_rp; static_eta0_rel = -1.0_rp
      a3_rel = -1.0_rp; window_moment_rel = -1.0_rp; dyson_u0_rel = -1.0_rp; pole_peak_rel = -1.0_rp
      pole_crossing_rel = -1.0_rp; electron_count_abs = -1.0_rp; eigen_contract_abs = -1.0_rp; moment_rel = -1.0_rp
      s1_rel = -1.0_rp; s2_rel = -1.0_rp; e1_fe_rel = -1.0_rp; q0_pole_abs = -1.0_rp
      open (newunit=u, file=trim(oracle_dir)//'/tolerances.nml', status='old', action='read')
      read (u, nml=spin_response_tolerances)
      close (u)
      call require('electron_count_abs', electron_count_abs)
      call require('eigen_contract_abs', eigen_contract_abs)
      call require('moment_rel', moment_rel)
      write (*, '(a)') 'test | measured | tolerance'
   end subroutine read_tolerances

   subroutine require(key, tol)
      character(*), intent(in) :: key
      real(rp), intent(in) :: tol

      if (.not. tol > 0.0_rp) then
         write (*, '(a)') 'tolerances.nml: key '//key//' missing or <= 0'
         error stop 1
      end if
   end subroutine require

   subroutine build_fe_state(reversed, control_obj, lattice_obj, charge_obj, energy_obj, hamiltonian_obj, fname)
      logical, intent(in) :: reversed
      character(*), intent(in), optional :: fname
      type(control), target, intent(inout) :: control_obj
      type(lattice), target, intent(inout) :: lattice_obj
      type(charge), target, intent(inout) :: charge_obj
      type(energy), target, intent(inout) :: energy_obj
      type(hamiltonian), target, intent(inout) :: hamiltonian_obj

      integer :: i

      if (present(fname)) then
         control_obj = control(fname)
      else
         control_obj = control('input.nml')
      end if
      lattice_obj = lattice(control_obj)
      call lattice_obj%build_data()
      call lattice_obj%bravais()
      call lattice_obj%structb(.true.)
      call lattice_obj%atomlist()
      if (reversed) call reverse_spin(lattice_obj)
      charge_obj = charge(lattice_obj)
      call charge_obj%bulkmat()
      energy_obj = energy(lattice_obj)
      hamiltonian_obj = hamiltonian(charge_obj)
      do i = 1, lattice_obj%nrec
         call lattice_obj%symbolic_atoms(i)%build_pot()
      end do
      if (control_obj%nsp == 2 .or. control_obj%nsp == 4) call hamiltonian_obj%build_lsham()
      call hamiltonian_obj%build_bulkham()
   end subroutine build_fe_state

   !> Pre-campaign mesh diagonalization on rec_ref, E_F from reciprocal's own solve, and the band moments of
   !> accumulate_spin_density_kspace (spin 1 along +z). band_moments(site, l+1, spin, 1:3) = q0, <E>, rms width.
   subroutine reference_band_moments(ham, rec_ref, bm)
      type(hamiltonian), target, intent(inout) :: ham
      type(reciprocal), intent(inout) :: rec_ref
      real(rp), allocatable, intent(out) :: bm(:, :, :, :)

      real(rp) :: axis(3, 1)

      rec_ref = reciprocal(ham)
      call rec_ref%generate_mp_mesh()
      call rec_ref%build_kspace_hamiltonian()
      call rec_ref%diagonalize_hamiltonian()
      rec_ref%fermi_level = rec_ref%find_fermi_level_from_eigenvalues(rec_ref%total_electrons)
      call rec_ref%accumulate_spin_density_kspace()
      call rec_ref%fill_band_moments_from_spin_density(trim(rec_ref%lattice%control%density_policy), &
                                                       reshape([0.0_rp, 0.0_rp, 1.0_rp], [3, 1]), axis)
      bm = rec_ref%band_moments
   end subroutine reference_band_moments

   subroutine reverse_spin(lattice_obj)
      type(lattice), intent(inout) :: lattice_obj

      integer :: i

      do i = 1, lattice_obj%nrec
         associate (pot => lattice_obj%symbolic_atoms(i)%potential)
            call swap_spin(pot%center_band)
            call swap_spin(pot%width_band)
            call swap_spin(pot%gravity_center)
            call swap_spin(pot%shifted_band)
            call swap_spin(pot%obar)
            call swap_spin(pot%pl)
            call swap_spin(pot%c)
            call swap_spin(pot%enu)
            call swap_spin(pot%ppar)
            call swap_spin(pot%qpar)
            call swap_spin(pot%srdel)
            call swap_spin(pot%vl)
            call swap_spin(pot%pnu)
            call swap_spin(pot%qi)
            call swap_spin(pot%dele)
            pot%ql(:, :, [1, 2]) = pot%ql(:, :, [2, 1])
         end associate
      end do
   end subroutine reverse_spin

   subroutine swap_spin(x)
      real(rp), allocatable, intent(inout) :: x(:, :)

      if (allocated(x)) x(:, [1, 2]) = x(:, [2, 1])
   end subroutine swap_spin

   subroutine report(label, measured, tol)
      character(*), intent(in) :: label
      real(rp), intent(in) :: measured, tol

      write (*, '(a,t74,es10.3,2x,es10.3,2x,a)') label, measured, tol, merge('ok  ', 'FAIL', measured <= tol)
      if (.not. measured <= tol) failed = .true.
   end subroutine report

   subroutine check(label, ok)
      character(*), intent(in) :: label
      logical, intent(in) :: ok

      write (*, '(a,t74,a)') label, merge('ok  ', 'FAIL', ok)
      if (.not. ok) failed = .true.
   end subroutine check

   subroutine finish()
      if (failed) then
         write (*, '(a)') 'RESULT: FAIL'
         error stop 1
      end if
      write (*, '(a)') 'RESULT: PASS'
   end subroutine finish

end module spin_response_fe_state
