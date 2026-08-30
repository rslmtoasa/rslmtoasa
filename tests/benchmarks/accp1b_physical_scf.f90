! SCF-B0C opt-in physical CPU/GPU SCF probe.
!
! This is deliberately a benchmark-only driver.  It constructs the same
! production SCF object graph as the bravais pre-processing path, installs an
! explicitly selected reciprocal backend on the cache before self%run(), and
! emits observables used by the FP32 precision study.  The normal rslmto
! executable and its namelist/default precision are not changed.
program accp1b_physical_scf
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: iso_c_binding, only: c_long_long
   use precision_mod, only: rp
   use basis_mod, only: basis_init, nb
   use math_mod, only: init_math_operators
   use mpi_mod, only: parallel_context, g_parallel_context, get_mpi_variables, rank
   use timer_mod, only: timer, g_timer
   use scf_benchmark_profile_mod, only: g_scf_benchmark_profile
   use logger_mod, only: g_logger
   use control_mod, only: control
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use energy_mod, only: energy
   use hamiltonian_mod, only: hamiltonian
   use sparse_mod, only: sparse
   use recursion_mod, only: recursion
   use density_of_states_mod, only: dos
   use green_mod, only: green
   use bands_mod, only: bands
   use mix_mod, only: mix
   use self_mod, only: self
   use reciprocal_mod, only: reciprocal, cuda_reciprocal_backend
   use rsrec_cuda_plugin_mod, only: rsrec_cuda_plugin_compiled
#ifdef USE_MPI
   use mpi
#endif
   implicit none

   character(len=512) :: input_file
   character(len=64) :: backend_name, solver_strategy, dos_method, benchmark_level
   character(len=64) :: scf_route, rs_solver, rs_backend
   real(rp) :: sigma, temperature
   integer :: status, nstep
   logical :: profile_requested
#ifdef USE_MPI
   logical :: mpi_is_initialized
   integer :: mpi_ierr
#endif

   input_file = 'input.nml'
   backend_name = 'lapack'
   solver_strategy = 'fp64_zheevd'
   dos_method = ''
   benchmark_level = 'scf_convergence'
   scf_route = 'reciprocal'
   rs_solver = ''
   rs_backend = ''
   sigma = -1.0_rp
   temperature = -1.0_rp
   nstep = 0
   profile_requested = .false.
   call parse_arguments(input_file, backend_name, solver_strategy, dos_method, benchmark_level, scf_route, rs_solver, rs_backend, &
                        sigma, temperature, nstep, profile_requested)
   if (trim(backend_name) /= 'lapack' .and. trim(backend_name) /= 'cuda') then
      error stop 'ACC-P1b: backend must be lapack or cuda'
   end if
   select case (trim(benchmark_level))
   case ('eigensolver', 'reciprocal_phase', 'rs_kernel', 'rs_electronic_structure', 'scf_iteration', 'scf_convergence')
   case default
      error stop 'ACC-P1b: invalid benchmark level'
   end select
   if (trim(scf_route) /= 'reciprocal' .and. trim(scf_route) /= 'real_space') then
      error stop 'SCF-B0C: scf_route must be reciprocal or real_space'
   end if
   if (trim(scf_route) == 'real_space') then
      if (len_trim(rs_solver) == 0) rs_solver = 'block'
      if (trim(rs_solver) /= 'block' .and. trim(rs_solver) /= 'chebyshev' .and. trim(rs_solver) /= 'lanczos') then
         error stop 'SCF-B0C: rs_solver must be block, chebyshev, or lanczos'
      end if
      if (len_trim(rs_backend) == 0) rs_backend = 'csr'
   end if

#ifdef USE_MPI
   call MPI_Initialized(mpi_is_initialized, mpi_ierr)
   if (.not. mpi_is_initialized) call MPI_Init(mpi_ierr)
#endif
   call g_parallel_context%restore_to_default()
   g_parallel_context = parallel_context()
   g_timer = timer()
   call init_math_operators()
   call g_logger%init()
   if (profile_requested) call g_scf_benchmark_profile%configure(.true.)

   call run_probe(trim(input_file), trim(backend_name), trim(solver_strategy), trim(dos_method), trim(benchmark_level), &
                  trim(scf_route), trim(rs_solver), trim(rs_backend), sigma, temperature, nstep, status)
   if (status /= 0) error stop 'ACC-P1b: physical SCF probe failed'
#ifdef USE_MPI
   call MPI_Finalize(mpi_ierr)
#endif

contains

   subroutine parse_arguments(input_file, backend_name, solver_strategy, dos_method, benchmark_level, scf_route, rs_solver, rs_backend, &
                              sigma, temperature, nstep, profile_requested)
      character(len=*), intent(inout) :: input_file, backend_name, solver_strategy, dos_method, benchmark_level
      character(len=*), intent(inout) :: scf_route, rs_solver, rs_backend
      real(rp), intent(inout) :: sigma, temperature
      integer, intent(inout) :: nstep
      logical, intent(inout) :: profile_requested
      integer :: i, narg
      character(len=128) :: key, value

      narg = command_argument_count()
      i = 1
      do while (i <= narg)
         call get_command_argument(i, key)
         if (trim(key) == '--profile') then
            profile_requested = .true.
            i = i + 1
            cycle
         end if
         if (i == narg) error stop 'ACC-P1b: option requires a value'
         call get_command_argument(i + 1, value)
         select case (trim(key))
         case ('--input')
            input_file = trim(value)
         case ('--backend')
            backend_name = trim(value)
         case ('--solver-strategy')
            solver_strategy = trim(value)
         case ('--dos-method')
            dos_method = trim(value)
         case ('--benchmark-level')
            benchmark_level = trim(value)
         case ('--scf-route')
            scf_route = trim(value)
         case ('--rs-solver')
            rs_solver = trim(value)
         case ('--rs-backend')
            rs_backend = trim(value)
         case ('--sigma')
            read(value, *) sigma
         case ('--temperature')
            read(value, *) temperature
         case ('--nstep')
            read(value, *) nstep
         case default
            error stop 'ACC-P1b: unknown option '//trim(key)
         end select
         i = i + 2
      end do
      if (sigma == 0.0_rp .or. temperature == 0.0_rp) error stop 'ACC-P1b: smearing values must be positive or omitted'
      if (nstep < 0) error stop 'ACC-P1b: nstep must be non-negative'
   end subroutine parse_arguments

   subroutine run_probe(input_file, backend_name, solver_strategy, dos_method, benchmark_level, scf_route, rs_solver, rs_backend, &
                        sigma, temperature, nstep, status)
      character(len=*), intent(in) :: input_file, backend_name, solver_strategy, dos_method, benchmark_level
      character(len=*), intent(in) :: scf_route, rs_solver, rs_backend
      real(rp), intent(in) :: sigma, temperature
      integer, intent(in) :: nstep
      integer, intent(out) :: status
      type(control), target :: ctl
      type(lattice), target :: lat
      type(charge), target :: chg
      type(energy), target :: en
      type(hamiltonian), target :: ham
      type(recursion), target :: rec
      type(dos), target :: dos_obj
      type(green), target :: green_obj
      type(bands), target :: bands_obj
      type(mix), target :: mix_obj
      type(self), target :: self_obj
      integer :: i, j, k, m, strategy_status
      real(rp) :: site_occupation, site_charge_transfer, near_ef_min, dos_states, value
      real(rp) :: final_fermi, total_charge, band_energy
      real(rp) :: near_ef_values(4), near_ef_distance(4)
      real(rp) :: backend_host_conversion, backend_host_staging, backend_host_widen, backend_sync
      integer(c_long_long) :: h2d_bytes, d2h_values_bytes, d2h_vectors_bytes
      integer(c_long_long) :: workspace_queries, workspace_reuses, gpu_free_bytes, gpu_total_bytes
      integer :: gpu_device, pinned_host_active
      character(len=32) :: converged_word
      logical :: rs_gpu_used, fallback_detected
      integer :: rs_nmat
      integer :: result_nmat, result_nk_unique, result_nk_mesh(3)

      status = 1
      backend_host_conversion = 0.0_rp
      backend_host_staging = 0.0_rp
      backend_host_widen = 0.0_rp
      backend_sync = 0.0_rp
      h2d_bytes = 0_c_long_long
      d2h_values_bytes = 0_c_long_long
      d2h_vectors_bytes = 0_c_long_long
      workspace_queries = 0_c_long_long
      workspace_reuses = 0_c_long_long
      gpu_free_bytes = 0_c_long_long
      gpu_total_bytes = 0_c_long_long
      gpu_device = -1
      pinned_host_active = 0
      rs_gpu_used = .false.
      fallback_detected = .false.
      rs_nmat = 0
      ctl = control(input_file)
      if (trim(scf_route) == 'real_space') then
         ctl%recur = trim(rs_solver)
         ctl%gpu_backend = trim(rs_backend)
         ctl%gpu_plugin = trim(backend_name) == 'cuda'
      end if
      lat = lattice(ctl)
      call lat%build_data()
      call lat%bravais()
      call lat%structb(.true.)
      call lat%atomlist()
      call basis_init(lat%symbolic_atoms(1)%potential%lmax)
      call get_mpi_variables(rank, lat%nrec)

      chg = charge(lat)
      call chg%bulkmat()
      en = energy(lat)
      ham = hamiltonian(chg)
      do i = 1, lat%nrec
         call lat%symbolic_atoms(i)%build_pot()
      end do
      if (ctl%nsp == 2 .or. ctl%nsp == 4) call ham%build_lsham()
      call ham%build_bulkham()
      mix_obj = mix(lat, chg)
      rec = recursion(ham, en, sparse(ham))
      dos_obj = dos(rec, en)
      green_obj = green(dos_obj)
      bands_obj = bands(green_obj)
      self_obj = self(bands_obj, mix_obj)
      self_obj%use_kspace = trim(scf_route) == 'reciprocal'
      if (nstep > 0) self_obj%nstep = nstep

      if (trim(scf_route) == 'reciprocal') then
         ! Build the cache early so the strategy is installed only for this
         ! explicit benchmark process.  The production self loop still owns
         ! all physics arrays and remains complex(rp)/real(rp).
         allocate(self_obj%reciprocal_scf_cache)
         self_obj%reciprocal_scf_cache = reciprocal(ham)
         call self_obj%reciprocal_scf_cache%make_execution_backend(backend_name)
         if (.not. allocated(self_obj%reciprocal_scf_cache%execution_backend)) then
            write (*, '(a)') 'SCF_B0C status=UNSUPPORTED reason=backend_initialization_failed'
            return
         end if
         if (trim(backend_name) == 'cuda') then
            select type (backend => self_obj%reciprocal_scf_cache%execution_backend)
            type is (cuda_reciprocal_backend)
               strategy_status = backend%set_solver_strategy(solver_strategy)
               if (strategy_status /= 0) then
                  write (*, '(a)') 'SCF_B0C status=UNSUPPORTED reason=solver_strategy_initialization_failed'
                  return
               end if
               gpu_device = backend%device
               call backend%memory_info(gpu_free_bytes, gpu_total_bytes)
            class default
               write (*, '(a)') 'SCF_B0C status=UNSUPPORTED reason=wrong_backend_type'
               return
            end select
         end if
         if (sigma > 0.0_rp) self_obj%reciprocal_scf_cache%gaussian_sigma = sigma
         if (temperature > 0.0_rp) self_obj%reciprocal_scf_cache%temperature = temperature
         if (len_trim(dos_method) > 0) self_obj%reciprocal_scf_cache%dos_method = trim(dos_method)
      else
         if (trim(backend_name) == 'cuda' .and. .not. rsrec_cuda_plugin_compiled()) then
            write (*, '(a)') 'SCF_B0C status=UNSUPPORTED reason=rs_cuda_plugin_not_compiled'
            return
         end if
      end if

      call self_obj%run()

      if (trim(scf_route) == 'reciprocal' .and. trim(backend_name) == 'cuda') then
         if (allocated(self_obj%reciprocal_scf_cache%execution_backend)) then
            select type (backend => self_obj%reciprocal_scf_cache%execution_backend)
            type is (cuda_reciprocal_backend)
               gpu_device = backend%device
               backend_host_conversion = backend%host_conversion_seconds
               backend_host_staging = backend%host_staging_seconds
               backend_host_widen = backend%host_widen_seconds
               backend_sync = backend%sync_seconds
               h2d_bytes = backend%h2d_bytes
               d2h_values_bytes = backend%d2h_values_bytes
               d2h_vectors_bytes = backend%d2h_vectors_bytes
               workspace_queries = backend%workspace_query_count
               workspace_reuses = backend%workspace_reuse_count
               pinned_host_active = backend%pinned_host_active
               call backend%memory_info(gpu_free_bytes, gpu_total_bytes)
            class default
               continue
            end select
         end if
      end if

      if (trim(scf_route) == 'real_space') then
         rs_gpu_used = associated(self_obj%recursion%gpu_backend)
         fallback_detected = trim(backend_name) == 'cuda' .and. .not. rs_gpu_used
         if (fallback_detected) then
            write (*, '(a)') 'SCF_B0C status=UNSUPPORTED reason=real_space_cuda_fallback'
            return
         end if
      end if

      site_occupation = sum(self_obj%symbolic_atom(lat%nbulk + 1)%potential%ql(1, :, :))
      site_charge_transfer = site_occupation - self_obj%symbolic_atom(lat%nbulk + 1)%element%valence
      final_fermi = self_obj%en%fermi
      total_charge = 0.0_rp
      do i = 1, lat%nrec
         total_charge = total_charge + sum(self_obj%symbolic_atom(lat%nbulk + i)%potential%ql(1, :, :))
      end do
      band_energy = self_obj%bands%eband
      near_ef_min = 0.0_rp
      near_ef_values = final_fermi
      near_ef_distance = huge(1.0_rp)
      dos_states = 0.0_rp
      if (trim(scf_route) == 'reciprocal') then
         final_fermi = self_obj%reciprocal_scf_cache%fermi_level
         total_charge = self_obj%reciprocal_scf_cache%canonical_electron_count
         band_energy = self_obj%reciprocal_scf_cache%canonical_band_energy
         near_ef_min = minval(abs(self_obj%reciprocal_scf_cache%eigenvalues - final_fermi))
         near_ef_values = final_fermi
         near_ef_distance = huge(1.0_rp)
         do k = 1, size(self_obj%reciprocal_scf_cache%eigenvalues, 2)
            do j = 1, size(self_obj%reciprocal_scf_cache%eigenvalues, 1)
               value = self_obj%reciprocal_scf_cache%eigenvalues(j, k)
               do m = 1, size(near_ef_values)
                  if (abs(value - final_fermi) < near_ef_distance(m)) then
                     if (m < size(near_ef_values)) then
                        near_ef_values(m+1:) = near_ef_values(m:size(near_ef_values)-1)
                        near_ef_distance(m+1:) = near_ef_distance(m:size(near_ef_distance)-1)
                     end if
                     near_ef_values(m) = value
                     near_ef_distance(m) = abs(value - final_fermi)
                     exit
                  end if
               end do
            end do
         end do
         if (allocated(self_obj%reciprocal_scf_cache%eigenvalues)) then
            dos_states = real(size(self_obj%reciprocal_scf_cache%eigenvalues, 1), rp)
            if (allocated(self_obj%reciprocal_scf_cache%k_weights)) then
               dos_states = dos_states * sum(self_obj%reciprocal_scf_cache%k_weights)
            end if
         end if
      else
         rs_nmat = nb*lat%nrec
         dos_states = total_charge
         ! The normal SCF loop does not evaluate the post-run band-energy
         ! diagnostic at this boundary.  Keep the field explicit rather than
         ! serializing an uninitialized bands%eband value.
         band_energy = 0.0_rp
      end if
      result_nmat = rs_nmat
      result_nk_unique = 0
      result_nk_mesh = 0
      if (trim(scf_route) == 'reciprocal') then
         result_nmat = size(self_obj%reciprocal_scf_cache%eigenvalues, 1)
         result_nk_unique = size(self_obj%reciprocal_scf_cache%eigenvalues, 2)
         result_nk_mesh = self_obj%reciprocal_scf_cache%nk_mesh
      end if
      converged_word = 'false'
      if (self_obj%converged) converged_word = 'true'
      write (*, '(*(a,1x))') &
         'SCF_B0C', &
         'status=PASS', &
         'scf_route='//trim(scf_route), &
         'backend='//trim(backend_name), &
         'solver_strategy='//trim(solver_strategy), &
         'converged='//trim(converged_word), &
         'electron_count='//real_token(total_charge), &
         'fermi_level='//real_token(final_fermi), &
         'band_energy='//real_token(band_energy), &
         'site_occupation='//real_token(site_occupation), &
         'site_charge_transfer='//real_token(site_charge_transfer), &
         'site_moment='//real_token(self_obj%symbolic_atom(lat%nbulk + 1)%potential%mtot), &
         'dos_state_count='//real_token(dos_states), &
         'near_ef_min_abs='//real_token(near_ef_min), &
         'near_ef_value_1='//real_token(near_ef_values(1)), &
         'near_ef_value_2='//real_token(near_ef_values(2)), &
         'near_ef_value_3='//real_token(near_ef_values(3)), &
         'near_ef_value_4='//real_token(near_ef_values(4)), &
         'scf_residual='//real_token(self_obj%mix%delta)
      write (*, '(A,1X,*(A,1X))') 'SCF_B0C_RESULT', &
         'benchmark_level='//trim(benchmark_level), &
         'scf_route='//trim(scf_route), &
         'backend='//trim(backend_name), &
         'solver_strategy='//trim(solver_strategy), &
         'rs_solver='//trim(merge(rs_solver, 'none', trim(scf_route) == 'real_space')), &
         'rs_backend='//trim(merge(rs_backend, 'none', trim(scf_route) == 'real_space')), &
         'natom='//int_token(lat%nrec), &
         'nmat='//int_token(result_nmat), &
         'nk_unique='//int_token(result_nk_unique), &
         'nk_mesh='//mesh_token(result_nk_mesh), &
         'final_total_energy='//real_token(self_obj%physical_total_energy), &
         'constraint_penalty_energy='//real_token(self_obj%constraint_penalty_energy), &
         'constraint_field_coupling_energy='//real_token(self_obj%constraint_field_coupling_energy), &
         'fermi_energy='//real_token(final_fermi), &
         'total_charge='//real_token(total_charge), &
         'electron_count='//real_token(total_charge), &
         'band_energy='//real_token(band_energy), &
         'site_charge_transfer='//real_token(site_charge_transfer), &
         'site_moment='//real_token(self_obj%symbolic_atom(lat%nbulk + 1)%potential%mtot), &
         'final_residual='//real_token(self_obj%mix%delta), &
         'converged='//trim(converged_word), &
         'starting_state_id=normal_initial', &
         'eigenvectors=required', &
         'recursion_depth='//int_token(merge(ctl%lld, 0, trim(scf_route) == 'real_space')), &
         'block_size='//int_token(merge(nb, 0, trim(scf_route) == 'real_space')), &
         'terminator='//int_token(merge(ctl%terminator, 0, trim(scf_route) == 'real_space')), &
         'chebyshev_order='//int_token(merge(ctl%lld, 0, trim(scf_route) == 'real_space' .and. trim(rs_solver) == 'chebyshev')), &
         'chebyshev_kernel='//trim(merge(ctl%cheb_backend, '                ', trim(scf_route) == 'real_space' .and. trim(rs_solver) == 'chebyshev')), &
         'spectral_bounds_policy='//trim(merge('resolved_runtime', 'none            ', trim(scf_route) == 'real_space' .and. trim(rs_solver) == 'chebyshev')), &
         'rs_gpu_used='//trim(merge('true ', 'false', rs_gpu_used)), &
         'fallback_detected='//trim(merge('true ', 'false', fallback_detected)), &
         'rs_kernel_correctness_status='//trim(merge('PASS', 'FAIL', .not. fallback_detected)), &
         'rs_kernel_invariant=finite_coefficients', &
         'rs_detail_timers_status='//trim(merge('not_exposed_by_backend', 'not_applicable        ', trim(scf_route) == 'real_space')), &
         'profile_enabled='//trim(merge('true ', 'false', g_scf_benchmark_profile%enabled)), &
         'gpu_device='//int_token(gpu_device), &
         'T_host_conversion='//real_token(backend_host_conversion), &
         'T_host_staging='//real_token(backend_host_staging), &
         'T_host_widen='//real_token(backend_host_widen), &
         'T_sync='//real_token(backend_sync), &
         'h2d_bytes='//int64_token(h2d_bytes), &
         'd2h_values_bytes='//int64_token(d2h_values_bytes), &
         'd2h_vectors_bytes='//int64_token(d2h_vectors_bytes), &
         'workspace_queries='//int64_token(workspace_queries), &
         'workspace_reuses='//int64_token(workspace_reuses), &
         'pinned_host_active='//int_token(pinned_host_active), &
         'gpu_free_bytes='//int64_token(gpu_free_bytes), &
         'gpu_total_bytes='//int64_token(gpu_total_bytes)
      status = 0
   end subroutine run_probe

   function real_token(value) result(token)
      real(rp), intent(in) :: value
      character(len=32) :: token
      write(token, '(es24.16)') value
      token = adjustl(token)
   end function real_token

   function int_token(value) result(token)
      integer, intent(in) :: value
      character(len=32) :: token
      write(token, '(i0)') value
      token = adjustl(token)
   end function int_token

   function int64_token(value) result(token)
      integer(c_long_long), intent(in) :: value
      character(len=32) :: token
      write(token, '(i0)') value
      token = adjustl(token)
   end function int64_token

   function mesh_token(value) result(token)
      integer, intent(in) :: value(3)
      character(len=64) :: token
      write(token, '(i0,a,i0,a,i0)') value(1), 'x', value(2), 'x', value(3)
      token = adjustl(token)
   end function mesh_token

end program accp1b_physical_scf
