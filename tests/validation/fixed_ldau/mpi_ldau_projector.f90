! MPI integration harness for the production reciprocal fixed-Liechtenstein
! occupation projector and its paired real/imaginary lattice LDM reduction.
program mpi_ldau_projector
   use mpi
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use logger_mod, only: g_logger
   use mpi_mod, only: parallel_context, rank, numprocs, ierr
   use kpoint_workset_mod, only: make_kpoint_workset
   use reciprocal_mod, only: reciprocal
   use lattice_mod, only: lattice
   implicit none
   integer, parameter :: nk = 2
   real(rp), parameter :: expected = -0.5_rp
   type(parallel_context) :: context
   type(reciprocal) :: rec
   type(lattice) :: lat
   real(rp) :: points(3,nk), weights(nk), error, restart_error
   integer :: local_nk
   logical :: failed

   call MPI_INIT(ierr)
   context = parallel_context()
   call g_logger%init()
   call basis_init(1)
   call rec%restore_to_default()
   nullify(rec%hamiltonian, rec%lattice, rec%control)

   lat%nbulk = 0
   lat%nrec = 1
   allocate(lat%symbolic_atoms(1))
   call lat%symbolic_atoms(1)%restore_to_default()
   lat%symbolic_atoms(1)%potential%lmax = 1
   lat%symbolic_atoms(1)%potential%ldm = cmplx(0.0_rp,0.0_rp,rp)

   points = 0.0_rp
   weights = 0.5_rp
   rec%k_mesh_distributed_active = numprocs > 1
   rec%nk_mesh = [nk,1,1]
   rec%nk_total = nk
   rec%use_symmetry_reduction = .false.
   rec%n_spin_components = 2
   rec%fermi_level = 0.0_rp
   rec%temperature = 0.0_rp
   rec%k_workset = make_kpoint_workset(points, weights, context, numprocs > 1, .true.)
   local_nk = rec%k_workset%nk_local
   allocate(rec%eigenvalues(1,local_nk), rec%eigenvectors(8,1,local_nk), rec%k_weights(nk))
   rec%k_weights = weights
   rec%eigenvalues = -1.0_rp
   rec%eigenvectors = cmplx(0.0_rp,0.0_rp,rp)
   rec%eigenvectors(2,1,:) = cmplx(1.0_rp/sqrt(2.0_rp),0.0_rp,rp)
   rec%eigenvectors(3,1,:) = cmplx(0.0_rp,1.0_rp/sqrt(2.0_rp),rp)

   call rec%calculate_ldm_from_projected_dos(lat, .true.)
   ! Exercise the same real+imag storage used by potential namelist restart.
   call lat%symbolic_atoms(1)%potential%flatten_ldm()
   lat%symbolic_atoms(1)%potential%ldm = cmplx(0.0_rp,0.0_rp,rp)
   call lat%symbolic_atoms(1)%potential%expand_ldm()
   restart_error = abs(aimag(lat%symbolic_atoms(1)%potential%ldm(2,1,1,2)) - expected)
   error = abs(aimag(lat%symbolic_atoms(1)%potential%ldm(2,1,1,2)) - expected)
   failed = error > 1.0e-12_rp .or. restart_error > 1.0e-12_rp .or. &
      abs(real(lat%symbolic_atoms(1)%potential%ldm(2,1,1,2),rp)) > 1.0e-12_rp .or. &
      abs(lat%symbolic_atoms(1)%potential%ldm(2,1,2,1) - conjg(lat%symbolic_atoms(1)%potential%ldm(2,1,1,2))) > 1.0e-12_rp
   if (rank == 0) then
      write(*,'(a,i0,a,i0,a,2es24.15,a,es12.4)') 'MPI_LDAU_PROJECTOR ranks=',numprocs, &
         ' nk_global=',nk,' P12=',real(lat%symbolic_atoms(1)%potential%ldm(2,1,1,2),rp), &
         aimag(lat%symbolic_atoms(1)%potential%ldm(2,1,1,2)),' error=',error
      write(*,'(a,es12.4)') 'restart_expand_error=',restart_error
      if (failed) write(*,'(a)') 'FAIL: production reciprocal LDM projector or complex reduction mismatch'
      if (.not. failed) write(*,'(a)') 'PASS: production reciprocal LDM projector and complex reduction'
   end if
   call MPI_FINALIZE(ierr)
   if (failed) error stop 1
end program mpi_ldau_projector
