!------------------------------------------------------------------------------
! DRESP-03Q -- production bcc-Fe adapter gate.
!
! The state is built through the ordinary lattice -> charge -> Hamiltonian ->
! reciprocal stack.  The finite-q fixture is then extracted from that live
! Hamiltonian and compared with the production reciprocal assembler at several
! k points.  The source Hamiltonian is intent(in) at the adapter boundary.
!------------------------------------------------------------------------------
program test_dresp03q_production_adapter
   use basis_mod, only: basis_init
   use control_mod, only: control
   use charge_mod, only: charge
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use logger_mod, only: g_logger
   use lr_kl_hessian_mod, only: lmto_live_hamiltonian_fixture, lmto_fixture_from_hamiltonian, &
      lmto_fixture_adapter_residual, assemble_lmto_hamiltonian
   use math_mod, only: init_math_operators
   use precision_mod, only: rp
   use reciprocal_mod, only: reciprocal
   use timer_mod, only: g_timer, timer
   implicit none

   type(control) :: ctl
   type(lattice), target :: lat
   type(charge), target :: chg
   type(hamiltonian), target :: ham
   type(reciprocal) :: recip
   type(lmto_live_hamiltonian_fixture) :: fixture
   real(rp), parameter :: kpoints(3,5) = reshape([ &
      0.0_rp, 0.0_rp, 0.0_rp, &
      0.137_rp, -0.091_rp, 0.053_rp, &
      -0.211_rp, 0.083_rp, 0.119_rp, &
      0.317_rp, 0.019_rp, -0.071_rp, &
      -0.403_rp, 0.127_rp, 0.041_rp], [3,5])
   complex(rp), allocatable :: h_production(:, :), h_fixture(:, :)
   real(rp) :: max_error, adapter_error, point_error
   integer :: ik, i, nmat

   call init_math_operators()
   call g_logger%init()
   g_timer = timer()

   ctl = control('input.nml')
   lat = lattice(ctl)
   call lat%build_data()
   call lat%bravais()
   call lat%structb(.true.)
   call lat%atomlist()
   call basis_init(lat%symbolic_atoms(1)%potential%lmax)

   chg = charge(lat)
   ham = hamiltonian(chg)
   do i = 1, lat%ntype
      call lat%symbolic_atoms(i)%build_pot()
   end do
   if (ctl%nsp == 2 .or. ctl%nsp == 4) call ham%build_lsham()
   call ham%build_bulkham()

   if (ham%hoh .or. ham%ccor_2c .or. ham%hubbard_u_general_check .or. ham%hubbard_v_check) then
      error stop 'DRESP-03Q production adapter gate is not the clean scalar-relativistic capability'
   end if
   recip = reciprocal(ham)
   recip%reciprocal_mode = 'ham_only'
   call lmto_fixture_from_hamiltonian(ham, fixture)
   call lmto_fixture_adapter_residual(ham, fixture, kpoints, adapter_error)

   nmat = 2*fixture%norb*fixture%nsite
   allocate(h_production(nmat,nmat), h_fixture(nmat,nmat))
   max_error = 0.0_rp
   do ik = 1, size(kpoints,2)
      call recip%build_hamiltonian_at_kpoint(kpoints(:,ik), h_production)
      call assemble_lmto_hamiltonian(fixture, kpoints(:,ik), h_fixture)
      point_error = maxval(abs(h_production-h_fixture))
      max_error = max(max_error, point_error)
      write (*, '(a,i0,a,es12.4)') 'DRESP-03Q production H(k) adapter error [k=', ik, '] = ', point_error
   end do
   write (*, '(a,es12.4)') 'DRESP-03Q production real-space adapter residual = ', adapter_error
   write (*, '(a,es12.4)') 'DRESP-03Q production reciprocal max error       = ', max_error
   if (adapter_error > 3.0e-12_rp .or. max_error > 3.0e-12_rp) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'
   call fixture%clear()
   deallocate(h_production, h_fixture)
end program test_dresp03q_production_adapter
