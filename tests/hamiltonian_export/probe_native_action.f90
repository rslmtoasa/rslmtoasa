! Frozen-state action oracle for the six material decks. Run in the deck
! directory with 'native' or 'import'. The native action comes from the
! unchanged production kernels, with full-support masks and a=1, b=0.
program probe_native_action
   use precision_mod, only: rp
   use basis_mod, only: nb
   use control_mod, only: control
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use energy_mod, only: energy
   use hamiltonian_mod, only: hamiltonian
   use recursion_mod, only: recursion
   use mpi_mod, only: g_parallel_context, parallel_context
   use timer_mod, only: g_timer, timer
#ifdef USE_MPI
   use mpi
#endif
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc, safe_alloc
#endif
   implicit none
   type(control), target :: ctl
   type(lattice), target :: lat
   type(charge), target :: ch
   type(energy), target :: en
   type(hamiltonian), target :: ham
   type(recursion) :: rec
   complex(rp), allocatable :: psi(:,:,:), action(:,:,:), scratch(:,:,:,:)
   character(len=16) :: mode
   integer :: i, j, k, u, ierr
   real(rp) :: dir(3)

#ifdef USE_MPI
   call MPI_Init(ierr)
#endif
   g_parallel_context = parallel_context()
   g_timer = timer()
#ifdef USE_SAFE_ALLOC
   g_safe_alloc = safe_alloc()
#endif
   call get_command_argument(1,mode)
   if (mode /= 'native' .and. mode /= 'import') error stop 'Usage: HamiltonianActionProbe native|import'
   ctl = control('input.nml')
   lat = lattice(ctl)
   call lat%build_data()
   call lat%bravais()
   if (ctl%calctype == 'I') then
      if (lat%nlay > 0) call lat%build_surf_full()
      call lat%newclu()
   end if
   call lat%structb(mode == 'native')
   call lat%atomlist()
   ch = charge(lat)
   if (ctl%calctype == 'I') then
      call ch%impmad()
   else
      call ch%bulkmat()
   end if
   en = energy(lat)
   ham = hamiltonian(ch)
   if (mode == 'import') then
      call ham%import_effective(en%fermi)
   else
      do i = 1, lat%ntype
         call lat%symbolic_atoms(i)%build_pot()
      end do
      if (ctl%has_soc()) call ham%build_lsham()
      call ham%build_bulkham()
      if (ctl%calctype == 'I') call ham%build_locham()
      if (ham%hubbard_u_general_check) then
         print '(a,i0)', 'Total valence electrons: ', nint(sum(lat%symbolic_atoms%element%valence))
         print '(a,es23.15)', 'Maximum native |U block| in Ry: ', maxval(abs(ham%hubbard_u_pot))
      end if
   end if
   rec%lattice => lat
   rec%control => ctl
   rec%hamiltonian => ham
   allocate(rec%izero(0:lat%kk),rec%idum(0:lat%kk))
   allocate(psi(nb,3,lat%kk),action(nb,3,lat%kk),scratch(nb,3,lat%kk,3))
   do i = 1, lat%kk
      do j = 1, 3
         do k = 1, nb
            psi(k,j,i) = cmplx(sin(real(k+i*j,rp)),cos(real(i*k-j,rp)),rp)
         end do
      end do
   end do
   open(newunit=u,file=trim(mode)//'.action',status='replace',access='stream',form='unformatted')
   write(u) nb, lat%kk
   write(u) psi
   do i = 0, 3
      rec%izero = 1
      rec%izero(0) = 0
      if (i == 0) then
         call rec%kubo_ham_vec_matmul(psi,action,scratch)
      else
         dir = 0.0_rp
         dir(i) = 1.0_rp
         if (ham%hoh) then
            call rec%velo_hoh_vec_matmul(psi,action,scratch,dir,.false.)
         else
            call rec%velo_vec_matmul(psi,action,dir,.false.)
         end if
      end if
      write(u) action
   end do
   close(u)
#ifdef USE_MPI
   call MPI_Finalize(ierr)
#endif
end program probe_native_action
