! ACC-P1b opt-in physical reciprocal-SCF probe.
!
! This is deliberately a benchmark-only driver.  It constructs the same
! production SCF object graph as the bravais pre-processing path, installs an
! explicitly selected reciprocal backend on the cache before self%run(), and
! emits observables used by the FP32 precision study.  The normal rslmto
! executable and its namelist/default precision are not changed.
program accp1b_physical_scf
   use, intrinsic :: iso_fortran_env, only: int64
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use math_mod, only: init_math_operators
   use mpi_mod, only: parallel_context, g_parallel_context, get_mpi_variables, rank
   use timer_mod, only: timer, g_timer
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
   implicit none

   character(len=512) :: input_file
   character(len=64) :: backend_name, solver_strategy, dos_method
   real(rp) :: sigma, temperature
   integer :: status, nstep

   input_file = 'input.nml'
   backend_name = 'lapack'
   solver_strategy = 'fp64_zheevd'
   dos_method = ''
   sigma = -1.0_rp
   temperature = -1.0_rp
   nstep = 0
   call parse_arguments(input_file, backend_name, solver_strategy, dos_method, sigma, temperature, nstep)
   if (trim(backend_name) /= 'lapack' .and. trim(backend_name) /= 'cuda') then
      error stop 'ACC-P1b: backend must be lapack or cuda'
   end if

   call g_parallel_context%restore_to_default()
   g_parallel_context = parallel_context()
   g_timer = timer()
   call init_math_operators()
   call g_logger%init()

   call run_probe(trim(input_file), trim(backend_name), trim(solver_strategy), trim(dos_method), sigma, temperature, nstep, status)
   if (status /= 0) error stop 'ACC-P1b: physical SCF probe failed'

contains

   subroutine parse_arguments(input_file, backend_name, solver_strategy, dos_method, sigma, temperature, nstep)
      character(len=*), intent(inout) :: input_file, backend_name, solver_strategy, dos_method
      real(rp), intent(inout) :: sigma, temperature
      integer, intent(inout) :: nstep
      integer :: i, narg
      character(len=128) :: key, value

      narg = command_argument_count()
      i = 1
      do while (i <= narg)
         if (i == narg) error stop 'ACC-P1b: option requires a value'
         call get_command_argument(i, key)
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

   subroutine run_probe(input_file, backend_name, solver_strategy, dos_method, sigma, temperature, nstep, status)
      character(len=*), intent(in) :: input_file, backend_name, solver_strategy, dos_method
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
      real(rp) :: near_ef_values(4), near_ef_distance(4)
      character(len=32) :: converged_word

      status = 1
      ctl = control(input_file)
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
      ! The physical probe is explicitly a reciprocal-SCF workflow, even
      ! when the source fixture's default is a real-space example.
      self_obj%use_kspace = .true.
      if (nstep > 0) self_obj%nstep = nstep

      ! Build the cache early so the strategy is installed only for this
      ! explicit benchmark process.  The production self loop still owns all
      ! physics arrays and remains complex(rp)/real(rp).
      allocate(self_obj%reciprocal_scf_cache)
      self_obj%reciprocal_scf_cache = reciprocal(ham)
      call self_obj%reciprocal_scf_cache%make_execution_backend(backend_name)
      if (.not. allocated(self_obj%reciprocal_scf_cache%execution_backend)) then
         write (*, '(a)') 'ACCP1B_SCF status=UNSUPPORTED reason=backend_initialization_failed'
         return
      end if
      if (trim(backend_name) == 'cuda') then
         select type (backend => self_obj%reciprocal_scf_cache%execution_backend)
         type is (cuda_reciprocal_backend)
            strategy_status = backend%set_solver_strategy(solver_strategy)
            if (strategy_status /= 0) then
               write (*, '(a)') 'ACCP1B_SCF status=UNSUPPORTED reason=solver_strategy_initialization_failed'
               return
            end if
         class default
            write (*, '(a)') 'ACCP1B_SCF status=UNSUPPORTED reason=wrong_backend_type'
            return
         end select
      end if
      if (sigma > 0.0_rp) self_obj%reciprocal_scf_cache%gaussian_sigma = sigma
      if (temperature > 0.0_rp) self_obj%reciprocal_scf_cache%temperature = temperature
      if (len_trim(dos_method) > 0) self_obj%reciprocal_scf_cache%dos_method = trim(dos_method)

      call self_obj%run()

      site_occupation = sum(self_obj%symbolic_atom(lat%nbulk + 1)%potential%ql(1, :, :))
      site_charge_transfer = site_occupation - self_obj%symbolic_atom(lat%nbulk + 1)%element%valence
      near_ef_min = minval(abs(self_obj%reciprocal_scf_cache%eigenvalues - self_obj%reciprocal_scf_cache%fermi_level))
      near_ef_values = self_obj%reciprocal_scf_cache%fermi_level
      near_ef_distance = huge(1.0_rp)
      do k = 1, size(self_obj%reciprocal_scf_cache%eigenvalues, 2)
         do j = 1, size(self_obj%reciprocal_scf_cache%eigenvalues, 1)
            value = self_obj%reciprocal_scf_cache%eigenvalues(j, k)
            do m = 1, size(near_ef_values)
               if (abs(value - self_obj%reciprocal_scf_cache%fermi_level) < near_ef_distance(m)) then
                  if (m < size(near_ef_values)) then
                     near_ef_values(m+1:) = near_ef_values(m:size(near_ef_values)-1)
                     near_ef_distance(m+1:) = near_ef_distance(m:size(near_ef_distance)-1)
                  end if
                  near_ef_values(m) = value
                  near_ef_distance(m) = abs(value - self_obj%reciprocal_scf_cache%fermi_level)
                  exit
               end if
            end do
         end do
      end do
      dos_states = 0.0_rp
      if (allocated(self_obj%reciprocal_scf_cache%eigenvalues)) then
         dos_states = real(size(self_obj%reciprocal_scf_cache%eigenvalues, 1), rp)
         if (allocated(self_obj%reciprocal_scf_cache%k_weights)) then
            dos_states = dos_states * sum(self_obj%reciprocal_scf_cache%k_weights)
         end if
      end if
      converged_word = 'false'
      if (self_obj%converged) converged_word = 'true'
      write (*, '(18(a,1x))') &
         'ACCP1B_SCF', &
         'status=PASS', &
         'backend='//trim(backend_name), &
         'solver_strategy='//trim(solver_strategy), &
         'converged='//trim(converged_word), &
         'electron_count='//real_token(self_obj%reciprocal_scf_cache%canonical_electron_count), &
         'fermi_level='//real_token(self_obj%reciprocal_scf_cache%fermi_level), &
         'band_energy='//real_token(self_obj%reciprocal_scf_cache%canonical_band_energy), &
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
      status = 0
   end subroutine run_probe

   function real_token(value) result(token)
      real(rp), intent(in) :: value
      character(len=32) :: token
      write(token, '(es24.16)') value
      token = adjustl(token)
   end function real_token

end program accp1b_physical_scf
