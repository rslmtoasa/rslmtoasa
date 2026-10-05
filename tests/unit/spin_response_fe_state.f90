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
   implicit none

   real(rp) :: electron_count_abs, eigen_contract_abs, moment_rel
   real(rp) :: a1_rel, a2_chi_rel, a2_identity_rel, static_eta0_rel, a3_rel, window_moment_rel, &
               dyson_u0_rel, pole_peak_rel, pole_crossing_rel
   namelist /spin_response_tolerances/ a1_rel, a2_chi_rel, a2_identity_rel, static_eta0_rel, a3_rel, &
      window_moment_rel, dyson_u0_rel, pole_peak_rel, pole_crossing_rel, &
      electron_count_abs, eigen_contract_abs, moment_rel
   logical :: failed = .false.

contains

   subroutine read_tolerances(oracle_dir)
      character(*), intent(in) :: oracle_dir

      integer :: u

      g_parallel_context = parallel_context()
      g_timer = timer()
      open (newunit=u, file=trim(oracle_dir)//'/tolerances.nml', status='old', action='read')
      read (u, nml=spin_response_tolerances)
      close (u)
      write (*, '(a)') 'test | measured | tolerance'
   end subroutine read_tolerances

   subroutine build_fe_state(reversed, control_obj, lattice_obj, charge_obj, energy_obj, hamiltonian_obj)
      logical, intent(in) :: reversed
      type(control), target, intent(inout) :: control_obj
      type(lattice), target, intent(inout) :: lattice_obj
      type(charge), target, intent(inout) :: charge_obj
      type(energy), target, intent(inout) :: energy_obj
      type(hamiltonian), target, intent(inout) :: hamiltonian_obj

      integer :: i

      control_obj = control('input.nml')
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
