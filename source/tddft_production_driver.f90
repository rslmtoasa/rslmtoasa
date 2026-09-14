!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Clean post-SCF production adapter for the rebuilt reciprocal TD-DFT
!> response services.
!>
!> This module owns input validation, accepted-state handoff, service request
!> assembly, and provenance output.  It deliberately owns no response
!> equation, radial augmentation formula, XC-kernel formula, Goldstone algebra,
!> or mode-extraction logic.
!------------------------------------------------------------------------------
module tddft_production_driver_mod

   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use mpi_mod, only: rank, numprocs
   use string_mod, only: lower, int2str
   use control_mod, only: control
   use lattice_mod, only: lattice
   use hamiltonian_mod, only: hamiltonian
   use energy_mod, only: energy
   use reciprocal_mod, only: reciprocal
   use recursion_mod, only: recursion
   use green_mod, only: green
   use symbolic_atom_mod, only: symbolic_atom
   use radial_ground_state_mod, only: radial_ground_state
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex
   use lr_response_space_mod, only: response_space_layout, response_operator_trace
   use lr_ks_susceptibility_mod, only: lr_electronic_state, lr_ks_susceptibility_request, &
      lr_ks_susceptibility_result, lr_snapshot_from_reciprocal, lr_q_endpoint_from_reciprocal, &
      evaluate_lr_ks_susceptibility, lr_product_ks_susceptibility_request, &
      lr_product_ks_susceptibility_result, evaluate_lr_product_ks_susceptibility
   use lr_pauli_transition_vertex_mod, only: pauli_endpoint_state, pauli_vertex_capabilities, &
      pauli_sigma_plus_matrix, pauli_sigma_minus_matrix, evaluate_pauli_transition_vertex
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, lmto_product_channel_plus, &
      lmto_product_channel_minus
   use lr_gf_susceptibility_mod, only: lr_gf_susceptibility_request, evaluate_lr_gf_susceptibility
   use lr_product_gf_susceptibility_mod, only: lr_product_gf_susceptibility_request, &
      lr_product_gf_susceptibility_result, evaluate_lr_product_gf_susceptibility
   use lr_rs_gf_susceptibility_mod, only: lr_rs_gf_provider, lr_rs_gf_pair, lr_rs_gf_susceptibility_request, &
      evaluate_lr_rs_gf_susceptibility
   use tddft_native_rsgf_provider_mod, only: tddft_native_rsgf_provider
   use lr_alsda_kernel_mod, only: lr_alsda_kernel_request, lr_alsda_kernel_result, evaluate_lr_alsda_kernel
   use lr_goldstone_sumrule_mod, only: lr_goldstone_sumrule_request, lr_goldstone_sumrule_result, &
      evaluate_lr_goldstone_sumrule
   use tddft_dyson_mod, only: tddft_dyson_request, tddft_dyson_result, evaluate_tddft_dyson, &
      lr_dyson_route_direct_alsda, lr_dyson_route_goldstone_sumrule
   use logger_mod, only: g_logger
   implicit none
   private

#ifdef VERSION
   character(len=*), parameter :: tddft_build_version = VERSION
#else
   character(len=*), parameter :: tddft_build_version = 'unavailable'
#endif

   character(len=*), parameter, public :: tddft_driver_backend_lehmann = 'lehmann'
   character(len=*), parameter, public :: tddft_driver_backend_reciprocal_gf = 'reciprocal_gf'
   character(len=*), parameter, public :: tddft_driver_backend_native_rsgf = 'native_rsgf'
   character(len=*), parameter, public :: tddft_driver_backend_product_lehmann = 'product_lehmann'
   character(len=*), parameter, public :: tddft_driver_backend_product_gf = 'product_gf'
   character(len=*), parameter, public :: tddft_driver_route_direct_alsda = lr_dyson_route_direct_alsda
   character(len=*), parameter, public :: tddft_driver_route_goldstone_sumrule = lr_dyson_route_goldstone_sumrule

   !> Parsed, minimal production input. Frequencies and q points are expanded
   !> once at the parser boundary so the production path never reparses input.
   type, public :: tddft_production_config
      logical :: present = .false.
      logical :: enabled = .false.
      integer :: nq = 1
      integer :: nfrequency = 1
      character(len=32) :: channel = 'chi_plus'
      real(rp), allocatable :: q_list(:, :) ! (3,nq)
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.01_rp
      integer :: response_lmax = -1
      character(len=48) :: interaction_route = tddft_driver_route_direct_alsda
      logical :: goldstone_correction = .false.
      character(len=32) :: backend = tddft_driver_backend_lehmann
      logical :: reciprocal_backend_crosscheck = .false.
      character(len=32) :: native_rsgf_provider = 'auto'
      integer :: gf_integration_points = 2001
      real(rp) :: gf_integration_eta = 0.0_rp
      real(rp) :: gf_energy_margin = 1.0_rp
      logical :: gf_closure_audit = .false.
      logical :: write_full_matrix = .true.
      character(len=256) :: output_file = 'tddft_response.dat'
   contains
      procedure :: restore_to_default => tddft_config_restore
   end type tddft_production_config

   !> Small capability contract used both by the real-object gate and focused
   !> rejection tests.  It makes the reason for a rejection observable without
   !> requiring a material fixture.
   type, public :: tddft_capability_state
      character(len=32) :: reciprocal_mode = 'ham_only'
      character(len=16) :: hamiltonian_order = 'second'
      integer :: basis_lmax = 2
      logical :: collinear = .true.
      logical :: has_soc = .false.
      logical :: orthogonal = .true.
      logical :: generalized_overlap = .false.
      logical :: has_extra_operator = .false.
      logical :: bulk_cell = .true.
   end type tddft_capability_state

   !> Driver result. The complete matrices are retained; no spectrum or mode
   !> reduction is invented at this orchestration layer.
   type, public :: tddft_production_result
      logical :: initialized = .false.
      integer :: nq = 0
      integer :: nfrequency = 0
      integer :: ndim = 0
      integer :: response_lmax = -1
      real(rp), allocatable :: q_list(:, :)
      real(rp), allocatable :: frequencies(:)
      complex(rp), allocatable :: ks_susceptibility(:, :, :, :)
      complex(rp), allocatable :: enhanced_susceptibility(:, :, :, :)
      complex(rp), allocatable :: loss_matrix(:, :, :, :)
      logical :: reciprocal_backend_crosscheck = .false.
      logical, allocatable :: reciprocal_crosscheck_valid(:, :) ! (frequency,q)
      real(rp), allocatable :: reciprocal_crosscheck_norm_lehmann(:, :) ! (frequency,q)
      real(rp), allocatable :: reciprocal_crosscheck_norm_gf(:, :) ! (frequency,q)
      real(rp), allocatable :: reciprocal_crosscheck_difference_frobenius(:, :) ! (frequency,q)
      real(rp), allocatable :: reciprocal_crosscheck_relative_frobenius(:, :) ! (frequency,q)
      real(rp), allocatable :: reciprocal_crosscheck_difference_infinity(:, :) ! (frequency,q)
      complex(rp), allocatable :: reciprocal_crosscheck_delta(:, :, :, :) ! (I,J,frequency,q)
      character(len=128) :: status = 'not evaluated'
      character(len=48) :: interaction_route = ''
      character(len=32) :: backend = ''
      character(len=64) :: goldstone_correction_status = 'not selected'
      character(len=256) :: interaction_provenance = ''
      character(len=256) :: bare_response_provenance = ''
   end type tddft_production_result

   public :: load_tddft_config
   public :: tddft_capability_is_supported
   public :: require_tddft_capability
   public :: validate_tddft_production_capability
   public :: run_tddft_production
   public :: evaluate_tddft_production_sweep

contains

   subroutine tddft_config_restore(this)
      class(tddft_production_config), intent(out) :: this
      this%present = .false.
      this%enabled = .false.
      this%nq = 1
      this%nfrequency = 1
      this%channel = 'chi_plus'
      this%eta = 0.01_rp
      this%response_lmax = -1
      this%interaction_route = tddft_driver_route_direct_alsda
      this%goldstone_correction = .false.
      this%backend = tddft_driver_backend_lehmann
      this%reciprocal_backend_crosscheck = .false.
      this%native_rsgf_provider = 'auto'
      this%gf_integration_points = 2001
      this%gf_integration_eta = 0.0_rp
      this%gf_energy_margin = 1.0_rp
      this%gf_closure_audit = .false.
      this%write_full_matrix = .true.
      this%output_file = 'tddft_response.dat'
      if (allocated(this%q_list)) deallocate(this%q_list)
      if (allocated(this%frequencies)) deallocate(this%frequencies)
      allocate(this%q_list(3, 1), this%frequencies(1))
      this%q_list = 0.0_rp
      this%frequencies = 0.0_rp
   end subroutine tddft_config_restore

   !> Read only the new minimal &tddft group. No legacy TD-DFT type is
   !> constructed and an absent group is the ordinary feature-off default.
   subroutine load_tddft_config(filename, config)
      character(len=*), intent(in) :: filename
      type(tddft_production_config), intent(out) :: config
      integer :: unit, ios, n_q_local = 1, n_omega_local = 1, iq, iw
      logical :: found
      character(len=512) :: line

      include 'include_codes/namelists/tddft.f90'

      call config%restore_to_default()

      enabled = .false.
      channel = 'chi_plus'
      n_q = 1
      q_list = 0.0_rp
      n_omega = 1
      use_omega_grid = .false.
      omega_grid = 0.0_rp
      omega_min = 0.0_rp
      omega_max = 0.0_rp
      eta = 0.01_rp
      response_lmax = -1
      interaction_route = tddft_driver_route_direct_alsda
      goldstone_correction = .false.
      backend = tddft_driver_backend_lehmann
      reciprocal_backend_crosscheck = .false.
      native_rsgf_provider = 'auto'
      gf_integration_points = 2001
      gf_integration_eta = 0.0_rp
      gf_energy_margin = 1.0_rp
      gf_closure_audit = .false.
      write_full_matrix = .true.
      output_file = 'tddft_response.dat'

      found = .false.
      open(newunit=unit, file=filename, action='read', status='old', iostat=ios)
      if (ios /= 0) return
      do
         read(unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         if (index(adjustl(lower(line)), '&tddft') == 1) then
            found = .true.
            exit
         end if
      end do
      close(unit)
      if (.not. found) return

      open(newunit=unit, file=filename, action='read', status='old', iostat=ios)
      if (ios /= 0) error stop 'TDDFT parser: input file cannot be reopened'
      read(unit, nml=tddft, iostat=ios)
      close(unit)
      if (ios /= 0) error stop 'TDDFT parser: invalid &tddft namelist'

      config%present = .true.
      config%enabled = enabled
      config%channel = trim(lower(channel))
      config%eta = eta
      config%response_lmax = response_lmax
      config%interaction_route = trim(lower(interaction_route))
      config%goldstone_correction = goldstone_correction
      config%backend = trim(lower(backend))
      config%reciprocal_backend_crosscheck = reciprocal_backend_crosscheck
      config%native_rsgf_provider = trim(lower(native_rsgf_provider))
      config%gf_integration_points = gf_integration_points
      config%gf_integration_eta = gf_integration_eta
      config%gf_energy_margin = gf_energy_margin
      config%gf_closure_audit = gf_closure_audit
      config%write_full_matrix = write_full_matrix
      config%output_file = trim(output_file)

      n_q_local = n_q
      n_omega_local = n_omega
      if (n_q_local < 1 .or. n_q_local > tddft_max_q) error stop 'TDDFT parser: n_q is outside [1,tddft_max_q]'
      if (n_omega_local < 1 .or. n_omega_local > tddft_max_omega) then
         error stop 'TDDFT parser: n_omega is outside [1,tddft_max_omega]'
      end if
      if (allocated(config%q_list)) deallocate(config%q_list)
      if (allocated(config%frequencies)) deallocate(config%frequencies)
      allocate(config%q_list(3, n_q_local), config%frequencies(n_omega_local))
      config%nq = n_q_local
      config%nfrequency = n_omega_local
      config%q_list = q_list(:, 1:n_q_local)
      if (use_omega_grid) then
         config%frequencies = omega_grid(1:n_omega_local)
      else if (n_omega_local == 1) then
         config%frequencies(1) = omega_min
      else
         do iw = 1, n_omega_local
            config%frequencies(iw) = omega_min + real(iw - 1, rp)*(omega_max - omega_min)/real(n_omega_local - 1, rp)
         end do
      end if
      do iq = 1, n_q_local
         if (any(config%q_list(:, iq) /= config%q_list(:, iq))) error stop 'TDDFT parser: q contains NaN'
      end do
      call validate_tddft_config(config)
   end subroutine load_tddft_config

   subroutine validate_tddft_config(config)
      type(tddft_production_config), intent(in) :: config
      character(len=128) :: route

      if (.not. config%enabled) return
      if (trim(config%channel) /= 'chi_plus' .and. trim(config%channel) /= 'chi_minus') then
         error stop 'TDDFT input: unsupported response channel; use chi_plus or chi_minus'
      end if
      if (config%eta <= 0.0_rp) error stop 'TDDFT input: eta must be positive'
      if (config%response_lmax < -1) error stop 'TDDFT input: response_lmax must be -1 or a nonnegative cutoff'
      if (config%response_lmax > 4) error stop 'TDDFT input: response_lmax above 4 is outside the validated sp/spd product baseline'
      route = trim(config%interaction_route)
      if (route /= tddft_driver_route_direct_alsda .and. route /= tddft_driver_route_goldstone_sumrule) then
         error stop 'TDDFT input: unsupported interaction_route; use direct_alsda or goldstone_sumrule'
      end if
      if (config%goldstone_correction) then
         error stop 'TDDFT input: Goldstone correction requires separate validated TDVAL evidence; it is not silently applied'
      end if
      if (trim(config%backend) /= tddft_driver_backend_lehmann .and. trim(config%backend) /= 'spectral' .and. &
          trim(config%backend) /= tddft_driver_backend_reciprocal_gf .and. &
          trim(config%backend) /= tddft_driver_backend_native_rsgf .and. &
          trim(config%backend) /= tddft_driver_backend_product_lehmann .and. &
          trim(config%backend) /= tddft_driver_backend_product_gf) then
         error stop 'TDDFT input: unsupported backend; use spectral/lehmann, reciprocal_gf, native_rsgf, product_lehmann or product_gf'
      end if
      if (trim(config%backend) == tddft_driver_backend_reciprocal_gf .or. &
          trim(config%backend) == tddft_driver_backend_native_rsgf .or. config%reciprocal_backend_crosscheck) then
         if (config%gf_integration_points < 3 .or. mod(config%gf_integration_points, 2) == 0) then
            error stop 'TDDFT input: reciprocal_gf or backend crosscheck requires an odd gf_integration_points value >= 3'
         end if
         if (config%gf_energy_margin <= 0.0_rp) error stop 'TDDFT input: reciprocal-GF crosscheck requires gf_energy_margin positive'
      end if
      if (trim(config%backend) == tddft_driver_backend_product_gf) then
         if (config%gf_integration_points < 3 .or. mod(config%gf_integration_points, 2) == 0) then
            error stop 'TDDFT input: product_gf requires an odd gf_integration_points value >= 3'
         end if
         if (config%gf_energy_margin <= 0.0_rp) error stop 'TDDFT input: product_gf requires gf_energy_margin positive'
      end if
      if (trim(config%backend) == tddft_driver_backend_native_rsgf) then
         if (trim(config%native_rsgf_provider) /= 'auto' .and. trim(config%native_rsgf_provider) /= 'block' .and. &
             trim(config%native_rsgf_provider) /= 'block_recursion' .and. trim(config%native_rsgf_provider) /= 'chebyshev') then
            error stop 'TDDFT input: native_rsgf_provider must be auto, block or chebyshev'
         end if
      end if
      if (len_trim(config%output_file) == 0) error stop 'TDDFT input: output_file must not be empty'
      if (route == tddft_driver_route_goldstone_sumrule) then
         if (.not. any(sum(abs(config%q_list), dim=1) <= 1.0e-12_rp)) then
            error stop 'TDDFT input: goldstone_sumrule requires q=(0,0,0) in q_list for the static reference'
         end if
      end if
   end subroutine validate_tddft_config

   logical function tddft_capability_is_supported(state, reason) result(ok)
      type(tddft_capability_state), intent(in) :: state
      character(len=*), intent(out) :: reason

      ok = .false.
      reason = 'unknown unsupported feature'
      if (.not. state%bulk_cell) then
         reason = 'non-bulk calculation geometry'
      else if (.not. state%collinear) then
         reason = 'noncollinear spin mode'
      else if (state%has_soc) then
         reason = 'spin-orbit coupling'
      else if (state%generalized_overlap .or. trim(state%reciprocal_mode) /= 'ham_only' .or. &
               .not. state%orthogonal) then
         reason = 'generalized overlap / non-orthogonal representation'
      else if (trim(state%hamiltonian_order) /= 'second') then
         reason = 'first-order Hamiltonian; validated TDDFT requires second order'
      else if (state%has_extra_operator) then
         reason = 'unsupported additive Hamiltonian operator'
      else if (state%basis_lmax /= 1 .and. state%basis_lmax /= 2) then
         reason = 'basis outside supported sp/spd (lmax=1 or 2)'
      else
         ok = .true.
         reason = 'supported collinear ham_only sp/spd baseline'
      end if
   end function tddft_capability_is_supported

   subroutine require_tddft_capability(state)
      type(tddft_capability_state), intent(in) :: state
      character(len=256) :: reason
      if (.not. tddft_capability_is_supported(state, reason)) then
         error stop 'TDDFT capability gate: unsupported feature: '//trim(reason)
      end if
   end subroutine require_tddft_capability

   subroutine validate_tddft_production_capability(config, control_obj, lattice_obj, hamiltonian_obj, reciprocal_obj)
      type(tddft_production_config), intent(in) :: config
      type(control), intent(in) :: control_obj
      type(lattice), intent(in) :: lattice_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj
      type(reciprocal), intent(in) :: reciprocal_obj
      type(tddft_capability_state) :: state
      character(len=256) :: reason

      if (.not. config%enabled) return
      state%bulk_cell = control_obj%calctype == 'B'
      state%collinear = control_obj%is_collinear()
      ! The lsham allocation is ordinary Hamiltonian storage and may exist as
      ! a zero-filled cache even when SOC is disabled; nsp is the capability
      ! authority for this production baseline.
      state%has_soc = control_obj%has_soc()
      state%reciprocal_mode = trim(reciprocal_obj%reciprocal_mode)
      state%generalized_overlap = trim(state%reciprocal_mode) /= 'ham_only' .or. allocated(reciprocal_obj%sk_overlap)
      state%orthogonal = .not. state%generalized_overlap
      state%hamiltonian_order = trim(reciprocal_obj%kspace_ham_order)
      if (trim(state%hamiltonian_order) == 'auto') then
         if (hamiltonian_obj%hoh) then
            state%hamiltonian_order = 'second'
         else
            state%hamiltonian_order = 'first'
         end if
      end if
      state%basis_lmax = -1
      if (allocated(lattice_obj%symbolic_atoms) .and. size(lattice_obj%symbolic_atoms) > 0) then
         state%basis_lmax = lattice_obj%symbolic_atoms(1)%potential%lmax
      end if
      state%has_extra_operator = hamiltonian_obj%local_axis .or. hamiltonian_obj%orb_pol .or. &
         hamiltonian_obj%ccor_2c .or. hamiltonian_obj%hubbard_u_general_check .or. &
         hamiltonian_obj%hubbard_u_impurity_check .or. hamiltonian_obj%hubbard_u_sc_check .or. &
         hamiltonian_obj%hubbard_v_check .or. control_obj%constraints_enable
      if (.not. tddft_capability_is_supported(state, reason)) then
         error stop 'TDDFT capability gate: unsupported feature: '//trim(reason)
      end if
      if (config%response_lmax > 2*state%basis_lmax) then
         error stop 'TDDFT capability gate: response_lmax exceeds the accepted angular-product cutoff'
      end if
      if (trim(config%backend) == tddft_driver_backend_native_rsgf) then
         if (lattice_obj%njij < 1) then
            error stop 'TDDFT capability gate: native_rsgf requires a complete lattice%ijpair set'
         end if
         if (numprocs /= 1) then
            error stop 'TDDFT capability gate: native_rsgf registered baseline requires a serial replicated pair workset'
         end if
         if (trim(config%native_rsgf_provider) == 'block' .or. trim(config%native_rsgf_provider) == 'block_recursion') then
            if (trim(lower(control_obj%recur)) /= 'block') then
               error stop 'TDDFT capability gate: native block provider requires control%recur=block'
            end if
         else if (trim(config%native_rsgf_provider) == 'chebyshev') then
            if (trim(lower(control_obj%recur)) /= 'chebyshev') then
               error stop 'TDDFT capability gate: native Chebyshev provider requires control%recur=chebyshev'
            end if
         else if (trim(lower(control_obj%recur)) /= 'block' .and. trim(lower(control_obj%recur)) /= 'chebyshev') then
            error stop 'TDDFT capability gate: native_rsgf requires control%recur=block or chebyshev'
         end if
      end if
   end subroutine validate_tddft_production_capability

   !> Run the response after SCF has accepted its state. The input objects are
   !> read-only except for the dedicated reciprocal eigenpair service cache.
   subroutine run_tddft_production(config, control_obj, lattice_obj, hamiltonian_obj, energy_obj, reciprocal_obj, &
                                   recursion_obj, green_obj, scf_converged)
      type(tddft_production_config), intent(in) :: config
      type(control), intent(in) :: control_obj
      type(lattice), target, intent(in) :: lattice_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj
      type(energy), intent(in) :: energy_obj
      type(reciprocal), intent(inout) :: reciprocal_obj
      type(recursion), target, intent(inout) :: recursion_obj
      type(green), target, intent(inout) :: green_obj
      logical, intent(in) :: scf_converged

      type(radial_ground_state), pointer :: ground_states(:)
      type(lmto_radial_basis), allocatable, target :: radial_bases(:)
      type(response_space_layout), target :: response_space
      type(lr_electronic_state), target :: left_state
      type(lr_electronic_state), allocatable, target :: endpoints(:)
      type(tddft_production_result) :: result
      type(tddft_native_rsgf_provider), target :: native_provider
      real(rp), allocatable, target :: native_site_positions(:, :)
      real(rp) :: ignore_real
      integer :: first, last, nsite, response_lmax, isite, iq

      if (.not. config%enabled) return
      call validate_tddft_production_capability(config, control_obj, lattice_obj, hamiltonian_obj, reciprocal_obj)
      if (.not. scf_converged) error stop 'TDDFT production driver: an accepted converged ground state is required'
      if (lattice_obj%nrec < 1) error stop 'TDDFT production driver: no accepted response sites are available'

      first = lattice_obj%nbulk + 1
      last = lattice_obj%nbulk + lattice_obj%nrec
      if (.not. allocated(lattice_obj%symbolic_atoms) .or. last > size(lattice_obj%symbolic_atoms)) then
         error stop 'TDDFT production driver: accepted symbolic-atom state is incomplete'
      end if
      ground_states => lattice_obj%symbolic_atoms(first:last)%radial_ground_state
      nsite = size(ground_states)
      do isite = 1, nsite
         if (.not. ground_states(isite)%valid .or. .not. ground_states(isite)%accepted .or. &
             .not. ground_states(isite)%pauli_basis_valid) then
            error stop 'TDDFT production driver: accepted LR-01 radial snapshot is incomplete'
         end if
      end do
      ! The response angular space is the product space of two LMTO
      ! endpoints, so its complete cutoff is twice the radial orbital cutoff.
      response_lmax = 2*ground_states(1)%pauli_lmax
      if (config%response_lmax >= 0) response_lmax = config%response_lmax
      if (response_lmax < 0 .or. response_lmax > 2*ground_states(1)%pauli_lmax) then
         error stop 'TDDFT production driver: response angular cutoff is outside the accepted radial basis'
      end if
      allocate(radial_bases(nsite))
      do isite = 1, nsite
         if (size(ground_states(isite)%r) /= size(ground_states(1)%r) .or. &
             any(abs(ground_states(isite)%r - ground_states(1)%r) > 1.0e-12_rp)) then
            error stop 'TDDFT production driver: response sites do not share the accepted radial mesh'
         end if
         if (ground_states(isite)%pauli_lmax /= ground_states(1)%pauli_lmax) then
            error stop 'TDDFT production driver: response sites do not share the accepted Pauli cutoff'
         end if
         call radial_basis_from_snapshot(ground_states(isite), radial_bases(isite))
         call radial_bases(isite)%require_supported(trim(reciprocal_obj%reciprocal_mode), &
            effective_hamiltonian_order(reciprocal_obj, hamiltonian_obj), .true., .true., &
            control_obj%has_soc(), .false.)
      end do
      call response_space%initialize(nsite, response_lmax, ground_states(1)%r, ground_states(1)%a, ground_states(1)%b, 1)

      ! The accepted EF is authoritative. This service rebuilds reciprocal
      ! eigenpairs from the accepted potential; it never solves a new SCF or
      ! searches for a different EF.
      reciprocal_obj%use_symmetry_reduction = .false.
      reciprocal_obj%use_time_reversal = .false.
      reciprocal_obj%dos_method = 'tetrahedron'
      reciprocal_obj%auto_find_fermi = .false.
      reciprocal_obj%fermi_level = energy_obj%fermi
      reciprocal_obj%kspace_ham_order = effective_hamiltonian_order(reciprocal_obj, hamiltonian_obj)
      call reciprocal_obj%generate_mp_mesh()
      call reciprocal_obj%build_kspace_hamiltonian()
      call reciprocal_obj%diagonalize_hamiltonian()
      reciprocal_obj%fermi_level = energy_obj%fermi
      ignore_real = reciprocal_obj%calculate_canonical_band_energy(.false.)
      call reciprocal_obj%require_replicated_k_workset('TDDFT production driver')
      call lr_snapshot_from_reciprocal(reciprocal_obj, left_state)

      allocate(endpoints(size(config%q_list, 2)))
      do iq = 1, size(config%q_list, 2)
         call lr_q_endpoint_from_reciprocal(reciprocal_obj, left_state, config%q_list(:, iq), endpoints(iq))
      end do
      if (trim(config%backend) == tddft_driver_backend_product_lehmann) then
         ! TDVK-02R2 is a bare-response validation seam only.  It stops at
         ! the naturally prepared reciprocal handoff and never enters KXC,
         ! Goldstone, Dyson, loss, or the dense point-space result container.
         call run_tddft_product_bare_smoke(config, response_space, radial_bases, ground_states, left_state, endpoints, &
            reciprocal_obj)
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_product_gf) then
         ! Product-GF remains a bare-response validation seam.  The optional
         ! closure audit reuses this already accepted reciprocal state for the
         ! Lehmann reference and every GF control sample.
         if (config%gf_closure_audit) then
            call run_tddft_product_gf_closure(config, response_space, radial_bases, left_state, endpoints)
         else
            call run_tddft_product_gf_smoke(config, response_space, radial_bases, left_state, endpoints)
         end if
         return
      end if
      if (trim(config%backend) == tddft_driver_backend_native_rsgf) then
         call native_provider%initialize(config%native_rsgf_provider, green_obj, recursion_obj, hamiltonian_obj, &
            lattice_obj, reciprocal_obj, radial_bases)
         allocate(native_site_positions(3, nsite))
         native_site_positions = 0.0_rp
         call evaluate_tddft_production_sweep(config, response_space, radial_bases, ground_states, left_state, endpoints, result, &
            native_provider, native_provider%pairs, native_site_positions)
      else
         call evaluate_tddft_production_sweep(config, response_space, radial_bases, ground_states, left_state, endpoints, result)
      end if
      call write_tddft_production_output(config, result, control_obj, lattice_obj, ground_states, response_space, reciprocal_obj)
      ! The production result has already been serialized. Release the dense
      ! response matrices before the accepted-state owner leaves scope; this
      ! is important for compact smoke runs with a full radial mesh.
      if (allocated(result%ks_susceptibility)) deallocate(result%ks_susceptibility)
      if (allocated(result%enhanced_susceptibility)) deallocate(result%enhanced_susceptibility)
      if (allocated(result%loss_matrix)) deallocate(result%loss_matrix)
      if (allocated(result%reciprocal_crosscheck_delta)) deallocate(result%reciprocal_crosscheck_delta)
      if (allocated(result%reciprocal_crosscheck_valid)) deallocate(result%reciprocal_crosscheck_valid)
      if (allocated(result%reciprocal_crosscheck_norm_lehmann)) deallocate(result%reciprocal_crosscheck_norm_lehmann)
      if (allocated(result%reciprocal_crosscheck_norm_gf)) deallocate(result%reciprocal_crosscheck_norm_gf)
      if (allocated(result%reciprocal_crosscheck_difference_frobenius)) deallocate(result%reciprocal_crosscheck_difference_frobenius)
      if (allocated(result%reciprocal_crosscheck_relative_frobenius)) deallocate(result%reciprocal_crosscheck_relative_frobenius)
      if (allocated(result%reciprocal_crosscheck_difference_infinity)) deallocate(result%reciprocal_crosscheck_difference_infinity)
      if (allocated(result%q_list)) deallocate(result%q_list)
      if (allocated(result%frequencies)) deallocate(result%frequencies)
   end subroutine run_tddft_production

   !> TDVK-02R2 validation-only handoff.  The accepted reciprocal snapshot is
   !> already naturally available here, so the live transition oracle and the
   !> compact bare response can be exercised without persisting eigenpairs or
   !> entering the full KXC/Dyson production lifecycle.
   subroutine run_tddft_product_bare_smoke(config, response_space, radial_bases, ground_states, left_state, endpoints, reciprocal_obj)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lmto_product_response_basis), target :: product_plus, product_minus
      type(lr_product_ks_susceptibility_request) :: request
      type(lr_product_ks_susceptibility_result) :: result
      integer :: gamma_index, channel_kind
      real(rp) :: wall_start, wall_end, maximum_transition_error
      real(rp) :: accepted_moment
      logical :: finite_response
      integer :: isite

      gamma_index = find_gamma_q(config%q_list)
      if (gamma_index > size(endpoints)) error stop 'TDVK-02R2 product smoke: Gamma endpoint is unavailable'
      call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      call product_minus%initialize(response_space, radial_bases, lmto_product_channel_minus, .true.)
      if (product_plus%product_dimension /= 232 .or. product_minus%product_dimension /= 232) then
         error stop 'TDVK-02R2 product smoke: accepted Fe product dimension is not 232'
      end if

      maximum_transition_error = 0.0_rp
      do channel_kind = 1, 2
         if (channel_kind == 1) then
            call run_product_transition_oracle(response_space, radial_bases, product_plus, left_state, endpoints(gamma_index), &
               lmto_product_channel_plus, maximum_transition_error)
         else
            call run_product_transition_oracle(response_space, radial_bases, product_minus, left_state, endpoints(gamma_index), &
               lmto_product_channel_minus, maximum_transition_error)
         end if
      end do

      request%q = config%q_list(:, gamma_index)
      request%frequencies = config%frequencies
      request%eta = config%eta
      request%channel = config%channel
      request%electronic_state => left_state
      request%q_endpoint_state => endpoints(gamma_index)
      if (trim(config%channel) == 'chi_plus') then
         request%product_basis => product_plus
      else
         request%product_basis => product_minus
      end if
      call cpu_time(wall_start)
      call evaluate_lr_product_ks_susceptibility(request, result)
      call cpu_time(wall_end)
      finite_response = all(ieee_is_finite(real(result%susceptibility, rp))) .and. &
         all(ieee_is_finite(aimag(result%susceptibility)))
      if (.not. finite_response) error stop 'TDVK-02R2 product smoke: compact response contains NaN or Inf'

      accepted_moment = 0.0_rp
      do isite = 1, size(ground_states)
         accepted_moment = accepted_moment + ground_states(isite)%integrated_moment_muB
      end do
      call write_tddft_product_bare_smoke(config, result, reciprocal_obj, accepted_moment, maximum_transition_error, &
         wall_end - wall_start)
      write (*, '(a,i0,a,es12.4,a,es12.4)') 'TDVK-02R2 compact Fe smoke: Nprod=', result%product_dimension, &
         ' runtime_cpu_s=', wall_end - wall_start, ' max_transition_error=', maximum_transition_error
   end subroutine run_tddft_product_bare_smoke

   !> TDVK-02R3 accepted-Fe compact reciprocal-GF smoke.  This is deliberately
   !> separate from the legacy point-grid reciprocal-GF backend and is not a
   !> production interaction/Dyson route.
   subroutine run_tddft_product_gf_smoke(config, response_space, radial_bases, left_state, endpoints)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(lmto_product_response_basis), target :: product_plus, product_minus
      type(lr_product_gf_susceptibility_request) :: request
      type(lr_product_gf_susceptibility_result) :: result
      integer :: gamma_index, channel_kind, i
      real(rp) :: frobenius, maximum_element
      complex(rp) :: trace
      logical :: finite_response

      gamma_index = find_gamma_q(config%q_list)
      if (gamma_index > size(endpoints)) error stop 'TDVK-02R3 product GF smoke: Gamma endpoint is unavailable'
      call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      call product_minus%initialize(response_space, radial_bases, lmto_product_channel_minus, .true.)
      if (product_plus%product_dimension /= 232 .or. product_minus%product_dimension /= 232) then
         error stop 'TDVK-02R3 product GF smoke: accepted Fe product dimension is not 232'
      end if

      request%q = config%q_list(:, gamma_index)
      request%frequencies = [0.0_rp]
      request%eta = config%eta
      request%channel = config%channel
      request%integration_points = config%gf_integration_points
      request%integration_eta = config%gf_integration_eta
      request%energy_margin = config%gf_energy_margin
      request%electronic_state => left_state
      request%q_endpoint_state => endpoints(gamma_index)
      if (trim(config%channel) == 'chi_plus') then
         request%product_basis => product_plus
      else
         request%product_basis => product_minus
      end if
      call evaluate_lr_product_gf_susceptibility(request, result)
      finite_response = all(ieee_is_finite(real(result%susceptibility, rp))) .and. &
         all(ieee_is_finite(aimag(result%susceptibility)))
      if (.not. finite_response) error stop 'TDVK-02R3 product GF smoke: compact response contains NaN or Inf'

      frobenius = sqrt(sum(abs(result%susceptibility(:, :, 1))**2))
      maximum_element = maxval(abs(result%susceptibility(:, :, 1)))
      trace = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, result%product_dimension
         trace = trace + result%susceptibility(i, i, 1)
      end do
      write (*, '(a)') 'TDVK-02R3 compact reciprocal-GF Fe smoke: execution/performance evidence only — not quadrature convergence'
      write (*, '(a,l1,a,i0,a,i0,a,es12.4)') '  finite=', finite_response, ' Nprod=', result%product_dimension, &
         ' integration_points=', result%integration_points, ' integration_eta=', result%actual_integration_eta
      write (*, '(a,es12.4,a,es12.4,a,2(es12.4,1x))') '  Frobenius_norm=', frobenius, &
         ' max_element=', maximum_element, ' trace=', real(trace, rp), aimag(trace)
      write (*, '(a,es12.4,a,es12.4,a,es12.4,a,i0,a,i0,a,i0,a,i0)') '  wall_s=', result%wall_time_seconds, &
         ' cpu_s=', result%cpu_time_seconds, ' time_per_energy_kpoint_s=', result%time_per_energy_kpoint, &
         ' energy_points=', result%integration_points, &
         ' k_points=', result%nk, ' frequencies=', result%nfrequency, &
         ' component_vertex_bytes=', result%component_vertex_memory_bytes
      write (*, '(a,i0,a,i0,a,i0)') '  gf_matrix_bytes=', result%gf_matrix_memory_bytes, &
         ' susceptibility_bytes=', result%susceptibility_memory_bytes, ' product_dimension=', result%product_dimension
   end subroutine run_tddft_product_gf_smoke

   !> TDVK-03 accepted-state closure audit.  The SCF handoff above has already
   !> prepared one immutable reciprocal state and all exact folded endpoints.
   !> This routine deliberately evaluates the compact Lehmann reference once,
   !> then varies only GF integration controls while keeping those snapshots
   !> and the complete 232-coordinate product basis fixed.
   subroutine run_tddft_product_gf_closure(config, response_space, radial_bases, left_state, endpoints)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(lmto_product_response_basis), target :: product_plus, product_minus
      type(lr_product_ks_susceptibility_request) :: lehmann_request
      type(lr_product_ks_susceptibility_result) :: lehmann_result
      real(rp), parameter :: eta_scale(4) = [10.0_rp, 5.0_rp, 2.5_rp, 1.0_rp]
      real(rp), parameter :: margin_increment(3) = [0.0_rp, 0.4_rp, 1.4_rp]
      real(rp) :: integration_eta, base_integration_eta, margin
      integer :: gamma_index, i, integration_points
      logical :: finite_response

      gamma_index = find_gamma_q(config%q_list)
      if (gamma_index > size(endpoints)) error stop 'TDVK-03 closure: Gamma endpoint is unavailable'
      call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      call product_minus%initialize(response_space, radial_bases, lmto_product_channel_minus, .true.)
      if (product_plus%product_dimension /= 232 .or. product_minus%product_dimension /= 232) then
         error stop 'TDVK-03 closure: accepted Fe product dimension is not 232'
      end if

      lehmann_request%q = config%q_list(:, gamma_index)
      lehmann_request%frequencies = [0.0_rp]
      lehmann_request%eta = config%eta
      lehmann_request%channel = config%channel
      lehmann_request%electronic_state => left_state
      lehmann_request%q_endpoint_state => endpoints(gamma_index)
      if (trim(config%channel) == 'chi_plus') then
         lehmann_request%product_basis => product_plus
      else
         lehmann_request%product_basis => product_minus
      end if
      call evaluate_lr_product_ks_susceptibility(lehmann_request, lehmann_result)
      finite_response = all(ieee_is_finite(real(lehmann_result%susceptibility, rp))) .and. &
         all(ieee_is_finite(aimag(lehmann_result%susceptibility)))
      if (.not. finite_response) error stop 'TDVK-03 closure: Lehmann reference contains NaN or Inf'

      base_integration_eta = config%gf_integration_eta
      if (base_integration_eta <= 0.0_rp) base_integration_eta = config%eta/40.0_rp
      if (base_integration_eta <= 0.0_rp .or. base_integration_eta >= config%eta) then
         error stop 'TDVK-03 closure: invalid base integration_eta'
      end if
      write (*, '(a)') 'TDVK-03 Fe reciprocal-backend closure: one accepted state, complete compact product space'
      write (*, '(a,i0,a,i0,a,i0,a,es16.8,a,es16.8,a,a)') '  accepted_state nbasis=', left_state%nbasis, &
         ' nbands=', left_state%nbands, ' nk=', left_state%nk, ' EF_Ry=', left_state%fermi_level, &
         ' temperature_K=', left_state%temperature, ' moment provenance=accepted LR-01 radial snapshot'
      write (*, '(a,3(es16.8,1x),a,a,a,i0)') '  q=', config%q_list(:, gamma_index), ' channel=', trim(config%channel), &
         ' product_dimension=', product_plus%product_dimension

      ! Width ladder.  The point count is derived from the same accepted
      ! energy span for each width so h/integration_eta remains controlled.
      do i = 1, size(eta_scale)
         integration_eta = eta_scale(i)*base_integration_eta
         integration_points = resolved_simpson_points(left_state, endpoints(gamma_index), integration_eta, &
            config%gf_energy_margin)
         if (i == size(eta_scale)) integration_points = max(integration_points, config%gf_integration_points)
         call report_product_gf_closure_sample('eta_ladder', config, response_space, product_plus, product_minus, &
            left_state, endpoints(gamma_index), lehmann_result, integration_points, integration_eta, config%gf_energy_margin)
      end do

      ! Simpson-resolution ladder at the narrowest controlled width.  Both
      ! grids are resolved; their difference isolates remaining mesh error.
      integration_points = resolved_simpson_points(left_state, endpoints(gamma_index), base_integration_eta, &
         config%gf_energy_margin)
      integration_points = max(integration_points, config%gf_integration_points)
      call report_product_gf_closure_sample('simpson_base', config, response_space, product_plus, product_minus, &
         left_state, endpoints(gamma_index), lehmann_result, integration_points, base_integration_eta, config%gf_energy_margin)
      call report_product_gf_closure_sample('simpson_fine', config, response_space, product_plus, product_minus, &
         left_state, endpoints(gamma_index), lehmann_result, 2*(integration_points - 1) + 1, base_integration_eta, &
         config%gf_energy_margin)

      ! Energy-window ladder.  Choose a resolved grid independently at every
      ! margin so the window comparison is not contaminated by coarse Simpson
      ! spacing.
      do i = 1, size(margin_increment)
         margin = config%gf_energy_margin + margin_increment(i)
         integration_points = resolved_simpson_points(left_state, endpoints(gamma_index), base_integration_eta, margin)
         call report_product_gf_closure_sample('window_ladder', config, response_space, product_plus, product_minus, &
            left_state, endpoints(gamma_index), lehmann_result, integration_points, base_integration_eta, margin)
      end do
      write (*, '(a)') 'TDVK-03 Fe reciprocal-backend closure: audit complete; evidence returned without physical interpretation'
   end subroutine run_tddft_product_gf_closure

   integer function resolved_simpson_points(left_state, right_state, integration_eta, margin) result(points)
      type(lr_electronic_state), intent(in) :: left_state, right_state
      real(rp), intent(in) :: integration_eta, margin
      real(rp) :: span
      integer :: intervals

      span = max(maxval(left_state%eigenvalues), maxval(right_state%eigenvalues)) - &
         min(minval(left_state%eigenvalues), minval(right_state%eigenvalues)) + 2.0_rp*margin
      intervals = ceiling(span/(0.4_rp*integration_eta))
      if (intervals < 2) intervals = 2
      if (mod(intervals, 2) /= 0) intervals = intervals + 1
      points = intervals + 1
   end function resolved_simpson_points

   subroutine report_product_gf_closure_sample(tag, config, response_space, product_plus, product_minus, left_state, endpoint, &
                                               lehmann_result, integration_points, integration_eta, energy_margin)
      character(len=*), intent(in) :: tag
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), intent(in) :: response_space
      type(lmto_product_response_basis), target, intent(in) :: product_plus, product_minus
      type(lr_electronic_state), target, intent(in) :: left_state, endpoint
      type(lr_product_ks_susceptibility_result), intent(in) :: lehmann_result
      integer, intent(in) :: integration_points
      real(rp), intent(in) :: integration_eta, energy_margin
      type(lr_product_gf_susceptibility_request) :: gf_request
      type(lr_product_gf_susceptibility_result) :: gf_result
      real(rp) :: norm_lehmann, norm_gf, difference_frobenius, relative_frobenius, difference_infinity
      logical :: finite_response
      type(lmto_product_response_basis), pointer :: product

      if (trim(config%channel) == 'chi_plus') then
         product => product_plus
      else
         product => product_minus
      end if
      gf_request%q = config%q_list(:, find_gamma_q(config%q_list))
      gf_request%frequencies = [0.0_rp]
      gf_request%eta = config%eta
      gf_request%channel = config%channel
      gf_request%integration_points = integration_points
      gf_request%integration_eta = integration_eta
      gf_request%energy_margin = energy_margin
      gf_request%product_basis => product
      gf_request%electronic_state => left_state
      gf_request%q_endpoint_state => endpoint
      call evaluate_lr_product_gf_susceptibility(gf_request, gf_result)
      finite_response = all(ieee_is_finite(real(gf_result%susceptibility, rp))) .and. &
         all(ieee_is_finite(aimag(gf_result%susceptibility)))
      if (.not. finite_response) error stop 'TDVK-03 closure: GF response contains NaN or Inf'
      norm_lehmann = sqrt(sum(abs(lehmann_result%susceptibility(:, :, 1))**2))
      norm_gf = sqrt(sum(abs(gf_result%susceptibility(:, :, 1))**2))
      difference_frobenius = sqrt(sum(abs(lehmann_result%susceptibility(:, :, 1) - &
         gf_result%susceptibility(:, :, 1))**2))
      relative_frobenius = difference_frobenius/max(norm_lehmann, norm_gf, tiny(1.0_rp))
      difference_infinity = maxval(abs(lehmann_result%susceptibility(:, :, 1) - gf_result%susceptibility(:, :, 1)))
      write (*, '(a,a,a,i0,4(a,es16.8))') '  ', trim(tag), ' N=', integration_points, &
         ' eta=', config%eta, ' integration_eta=', integration_eta, ' margin=', energy_margin, &
         ' energy_min=', gf_result%energy_min
      write (*, '(a,8(a,es16.8))') '    ', ' energy_max=', gf_result%energy_max, ' h=', gf_result%energy_spacing, &
         ' h_over_eta=', gf_result%spacing_over_integration_eta, ' norm_lehmann=', norm_lehmann, &
         ' norm_gf=', norm_gf, ' dF=', difference_frobenius, ' rF=', relative_frobenius, &
         ' dInf=', difference_infinity
      write (*, '(a,a,es16.8)') '    ', ' wall_s=', gf_result%wall_time_seconds
   end subroutine report_product_gf_closure_sample

   subroutine run_product_transition_oracle(response_space, radial_bases, product, left_state, right_state, channel_kind, maximum_error)
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(lmto_product_response_basis), intent(in) :: product
      type(lr_electronic_state), intent(in) :: left_state, right_state
      integer, intent(in) :: channel_kind
      real(rp), intent(inout) :: maximum_error
      type(pauli_vertex_capabilities) :: capabilities
      type(pauli_endpoint_state) :: left_band, right_band
      complex(rp), allocatable :: point_transition(:), product_coordinates(:), reference_coordinates(:)
      complex(rp) :: operator_matrix(2, 2)
      integer :: selected_left(3), selected_right(3), ik, itransition
      real(rp) :: error

      call select_live_transitions(product, left_state, right_state, selected_left, selected_right, ik)
      allocate(point_transition(response_space%ndim), product_coordinates(product%product_dimension), &
         reference_coordinates(product%product_dimension))
      capabilities = pauli_vertex_capabilities()
      if (channel_kind == lmto_product_channel_plus) then
         operator_matrix = pauli_sigma_plus_matrix()
      else
         operator_matrix = pauli_sigma_minus_matrix()
      end if
      do itransition = 1, size(selected_left)
         call left_band%initialize(left_state%eigenvalues(selected_left(itransition), ik), &
            left_state%eigenvectors(:, selected_left(itransition), ik))
         call right_band%initialize(right_state%eigenvalues(selected_right(itransition), ik), &
            right_state%eigenvectors(:, selected_right(itransition), ik))
         call evaluate_pauli_transition_vertex(response_space, radial_bases, left_band, right_band, operator_matrix, &
            capabilities, point_transition)
         call product%transition_coordinates(left_band, right_band, product_coordinates)
         call project_point_transition(response_space, product, point_transition, reference_coordinates)
         error = sqrt(sum(abs(product_coordinates - reference_coordinates)**2))/ &
            max(sqrt(sum(abs(reference_coordinates)**2)), epsilon(1.0_rp))
         maximum_error = max(maximum_error, error)
         write (*, '(a,a,a,i0,a,es12.4)') 'TDVK-02R2 live transition channel=', &
            product_channel_name(channel_kind), ' index=', itransition, ' residual=', error
         if (error >= 1.0e-10_rp .or. .not. ieee_is_finite(error)) then
            error stop 'TDVK-02R2 live transition oracle failed'
         end if
      end do
      deallocate(point_transition, product_coordinates, reference_coordinates)
   end subroutine run_product_transition_oracle

   subroutine select_live_transitions(product, left_state, right_state, selected_left, selected_right, selected_k)
      type(lmto_product_response_basis), intent(in) :: product
      type(lr_electronic_state), intent(in) :: left_state, right_state
      integer, intent(out) :: selected_left(:), selected_right(:), selected_k
      type(pauli_endpoint_state) :: left_band, right_band
      complex(rp), allocatable :: coordinates(:)
      real(rp) :: score, norm, nearest_score, deepest_energy
      real(rp), parameter :: occupation_tolerance = 1.0e-8_rp, coordinate_tolerance = 1.0e-12_rp
      integer :: ik, ib, jb, nearest_left, nearest_right, deep_left, deep_right, first_left, first_right
      integer :: additional_left, additional_right

      if (size(selected_left) < 3 .or. size(selected_right) /= size(selected_left)) then
         error stop 'TDVK-02R2 transition selector: output shape mismatch'
      end if
      selected_k = 1
      do ik = 2, left_state%nk
         if (sum(left_state%k_points(:, ik)**2) < sum(left_state%k_points(:, selected_k)**2)) selected_k = ik
      end do
      allocate(coordinates(product%product_dimension))
      nearest_score = huge(1.0_rp)
      deepest_energy = huge(1.0_rp)
      nearest_left = 0
      nearest_right = 0
      deep_left = 0
      deep_right = 0
      first_left = 0
      first_right = 0
      do ib = 1, left_state%nbands
         do jb = 1, right_state%nbands
            if (left_state%occupations(ib, selected_k) <= 0.5_rp + occupation_tolerance .or. &
                right_state%occupations(jb, selected_k) >= 0.5_rp - occupation_tolerance .or. &
                left_state%occupations(ib, selected_k) - right_state%occupations(jb, selected_k) <= occupation_tolerance) cycle
            call left_band%initialize(left_state%eigenvalues(ib, selected_k), left_state%eigenvectors(:, ib, selected_k))
            call right_band%initialize(right_state%eigenvalues(jb, selected_k), right_state%eigenvectors(:, jb, selected_k))
            call product%transition_coordinates(left_band, right_band, coordinates)
            norm = sqrt(sum(abs(coordinates)**2))
            if (norm <= coordinate_tolerance) cycle
            if (first_left == 0) then
               first_left = ib
               first_right = jb
            end if
            score = abs(left_state%eigenvalues(ib, selected_k) - left_state%fermi_level) + &
               abs(right_state%eigenvalues(jb, selected_k) - right_state%fermi_level)
            if (score < nearest_score) then
               nearest_score = score
               nearest_left = ib
               nearest_right = jb
            end if
            if (left_state%eigenvalues(ib, selected_k) < left_state%fermi_level - 0.20_rp .and. &
                left_state%eigenvalues(ib, selected_k) < deepest_energy) then
               deepest_energy = left_state%eigenvalues(ib, selected_k)
               deep_left = ib
               deep_right = jb
            end if
         end do
      end do
      if (nearest_left == 0 .or. first_left == 0) error stop 'TDVK-02R2 transition selector: no nonzero Fe Gamma transition'
      if (deep_left == 0) then
         deep_left = first_left
         deep_right = first_right
      end if
      additional_left = 0
      additional_right = 0
      do ib = 1, left_state%nbands
         do jb = 1, right_state%nbands
            if (ib == nearest_left .and. jb == nearest_right) cycle
            if (ib == deep_left .and. jb == deep_right) cycle
            if (left_state%occupations(ib, selected_k) <= 0.5_rp + occupation_tolerance .or. &
                right_state%occupations(jb, selected_k) >= 0.5_rp - occupation_tolerance .or. &
                left_state%occupations(ib, selected_k) - right_state%occupations(jb, selected_k) <= occupation_tolerance) cycle
            call left_band%initialize(left_state%eigenvalues(ib, selected_k), left_state%eigenvectors(:, ib, selected_k))
            call right_band%initialize(right_state%eigenvalues(jb, selected_k), right_state%eigenvectors(:, jb, selected_k))
            call product%transition_coordinates(left_band, right_band, coordinates)
            if (sqrt(sum(abs(coordinates)**2)) > coordinate_tolerance) then
               additional_left = ib
               additional_right = jb
               exit
            end if
         end do
         if (additional_left /= 0) exit
      end do
      if (additional_left == 0) error stop 'TDVK-02R2 transition selector: fewer than three distinct transitions'
      selected_left = [nearest_left, deep_left, additional_left]
      selected_right = [nearest_right, deep_right, additional_right]
      deallocate(coordinates)
   end subroutine select_live_transitions

   subroutine project_point_transition(response_space, product, point_transition, product_coordinates)
      type(response_space_layout), intent(in) :: response_space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: point_transition(:)
      complex(rp), intent(out) :: product_coordinates(:)
      type(response_super_index) :: item
      integer :: product_flat, site, response_l, response_m, product_mode, ir, point_flat

      if (size(point_transition) /= response_space%ndim .or. size(product_coordinates) /= product%product_dimension) then
         error stop 'TDVK-02R2 transition projection: vector shape mismatch'
      end if
      product_coordinates = cmplx(0.0_rp, 0.0_rp, rp)
      do product_flat = 1, product%product_dimension
         call product%unflatten_index(product_flat, site, response_l, response_m, product_mode)
         do ir = 1, response_space%npoint
            item = response_super_index(site, response_l, response_m, ir, 1)
            call response_flatten_superindex(item, response_space%nsite, response_space%response_lmax, &
               response_space%npoint, response_space%nchannel, point_flat)
            product_coordinates(product_flat) = product_coordinates(product_flat) + &
               conjg(product%blocks(site, response_l)%weighted_modes(ir, product_mode))* &
               sqrt(response_space%radial_weights(ir))*point_transition(point_flat)
         end do
      end do
   end subroutine project_point_transition

   subroutine write_tddft_product_bare_smoke(config, result, reciprocal_obj, accepted_moment, maximum_transition_error, runtime)
      type(tddft_production_config), intent(in) :: config
      type(lr_product_ks_susceptibility_result), intent(in) :: result
      type(reciprocal), intent(in) :: reciprocal_obj
      real(rp), intent(in) :: accepted_moment, maximum_transition_error, runtime
      integer :: unit, ifrequency, i
      real(rp) :: frobenius_norm, maximum_element, memory_bytes
      complex(rp) :: trace

      memory_bytes = 16.0_rp*real(result%product_dimension*result%product_dimension*size(result%frequencies), rp)
      open(newunit=unit, file=trim(config%output_file), status='replace', action='write')
      write(unit, '(a)') '# TDVK-02R2 compact Lehmann bare response; no KXC or Dyson invoked'
      write(unit, '(a,a)') '# build_version = ', trim(tddft_build_version)
      write(unit, '(a,i0)') '# product_dimension = ', result%product_dimension
      write(unit, '(a,es24.16)') '# response_matrix_memory_bytes = ', memory_bytes
      write(unit, '(a,es24.16)') '# response_matrix_memory_MiB = ', memory_bytes/(1024.0_rp**2)
      write(unit, '(a,es24.16)') '# runtime_cpu_seconds = ', runtime
      write(unit, '(a,es24.16)') '# accepted_moment_muB = ', accepted_moment
      write(unit, '(a,es24.16)') '# EF_Ry = ', reciprocal_obj%fermi_level
      write(unit, '(a,es24.16)') '# eta_Ry = ', config%eta
      write(unit, '(a,es24.16)') '# maximum_live_transition_residual = ', maximum_transition_error
      write(unit, '(a,l1)') '# finite_response = ', all(ieee_is_finite(real(result%susceptibility, rp))) .and. &
         all(ieee_is_finite(aimag(result%susceptibility)))
      write(unit, '(a)') '# columns: omega_Ry frobenius_norm max_abs_element trace_real trace_imag'
      do ifrequency = 1, size(result%frequencies)
         frobenius_norm = sqrt(sum(abs(result%susceptibility(:, :, ifrequency))**2))
         maximum_element = maxval(abs(result%susceptibility(:, :, ifrequency)))
         trace = cmplx(0.0_rp, 0.0_rp, rp)
         do i = 1, result%product_dimension
            trace = trace + result%susceptibility(i, i, ifrequency)
         end do
         write(unit, '(5(es24.16,1x))') result%frequencies(ifrequency), frobenius_norm, maximum_element, &
            real(trace, rp), aimag(trace)
      end do
      close(unit)
   end subroutine write_tddft_product_bare_smoke

   character(len=5) function product_channel_name(channel_kind) result(value)
      integer, intent(in) :: channel_kind

      if (channel_kind == lmto_product_channel_plus) then
         value = 'plus '
      else
         value = 'minus'
      end if
   end function product_channel_name

   !> Evaluate an already prepared request batch. This is the reproducibility
   !> seam: tests and validation can compare it directly with service calls,
   !> without reconstructing SCF state.
   subroutine evaluate_tddft_production_sweep(config, response_space, radial_bases, ground_states, left_state, endpoints, result, &
                                              native_provider, native_pairs, native_site_positions)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(radial_ground_state), target, intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      type(tddft_production_result), intent(out) :: result
      class(lr_rs_gf_provider), target, intent(inout), optional :: native_provider
      type(lr_rs_gf_pair), intent(in), optional :: native_pairs(:)
      real(rp), target, intent(in), optional :: native_site_positions(:, :)

      type(lr_alsda_kernel_request) :: kxc_request
      type(lr_alsda_kernel_result) :: kxc_result
      type(lr_goldstone_sumrule_request) :: gsr_request
      type(lr_goldstone_sumrule_result) :: gsr_result
      type(lr_ks_susceptibility_result) :: bare_result, static_result
      type(tddft_dyson_request) :: dyson_request
      type(tddft_dyson_result) :: dyson_result
      complex(rp), allocatable :: interaction(:, :)
      real(rp), allocatable :: magnetization(:, :), static_frequency(:)
      integer :: iq, i0, ndim

      call validate_prepared_sweep(config, response_space, radial_bases, ground_states, left_state, endpoints, &
         native_provider, native_pairs)
      ndim = response_space%ndim
      result%nq = size(config%q_list, 2)
      result%nfrequency = size(config%frequencies)
      result%ndim = ndim
      result%response_lmax = response_space%response_lmax
      result%interaction_route = trim(config%interaction_route)
      result%backend = trim(config%backend)
      result%reciprocal_backend_crosscheck = config%reciprocal_backend_crosscheck
      result%bare_response_provenance = ''
      result%q_list = config%q_list
      result%frequencies = config%frequencies
      allocate(result%ks_susceptibility(ndim, ndim, result%nfrequency, result%nq), &
               result%enhanced_susceptibility(ndim, ndim, result%nfrequency, result%nq), &
               result%loss_matrix(ndim, ndim, result%nfrequency, result%nq))
      result%ks_susceptibility = cmplx(0.0_rp, 0.0_rp, rp)
      result%enhanced_susceptibility = cmplx(0.0_rp, 0.0_rp, rp)
      result%loss_matrix = cmplx(0.0_rp, 0.0_rp, rp)
      if (config%reciprocal_backend_crosscheck) then
         allocate(result%reciprocal_crosscheck_valid(result%nfrequency, result%nq), &
            result%reciprocal_crosscheck_norm_lehmann(result%nfrequency, result%nq), &
            result%reciprocal_crosscheck_norm_gf(result%nfrequency, result%nq), &
            result%reciprocal_crosscheck_difference_frobenius(result%nfrequency, result%nq), &
            result%reciprocal_crosscheck_relative_frobenius(result%nfrequency, result%nq), &
            result%reciprocal_crosscheck_difference_infinity(result%nfrequency, result%nq), &
            result%reciprocal_crosscheck_delta(ndim, ndim, result%nfrequency, result%nq))
         result%reciprocal_crosscheck_valid = .false.
         result%reciprocal_crosscheck_norm_lehmann = 0.0_rp
         result%reciprocal_crosscheck_norm_gf = 0.0_rp
         result%reciprocal_crosscheck_difference_frobenius = 0.0_rp
         result%reciprocal_crosscheck_relative_frobenius = 0.0_rp
         result%reciprocal_crosscheck_difference_infinity = 0.0_rp
         result%reciprocal_crosscheck_delta = cmplx(0.0_rp, 0.0_rp, rp)
      end if

      allocate(magnetization(size(ground_states), response_space%npoint))
      do iq = 1, size(ground_states)
         magnetization(iq, :) = ground_states(iq)%n_up - ground_states(iq)%n_down
      end do

      kxc_request%response_space => response_space
      kxc_request%ground_states => ground_states
      allocate(kxc_request%pauli_magnetization(size(ground_states), response_space%npoint))
      kxc_request%pauli_magnetization = magnetization

      select case (trim(config%interaction_route))
      case (tddft_driver_route_direct_alsda)
         call evaluate_lr_alsda_kernel(kxc_request, kxc_result)
         allocate(interaction(ndim, ndim))
         interaction = kxc_result%canonical_operator
         result%interaction_provenance = 'KXC-01 direct ALSDA LR-03'
         result%goldstone_correction_status = 'disabled'
      case (tddft_driver_route_goldstone_sumrule)
         i0 = find_gamma_q(config%q_list)
         allocate(static_frequency(1))
         static_frequency(1) = 0.0_rp
         call evaluate_bare_response(config, response_space, radial_bases, left_state, endpoints(i0), &
                                     config%q_list(:, i0), static_frequency, static_result, native_provider, native_pairs, &
                                     native_site_positions)
         gsr_request%response_space => response_space
         gsr_request%static_susceptibility = static_result%susceptibility(:, :, 1)
         gsr_request%magnetization = magnetization
         call evaluate_lr_goldstone_sumrule(gsr_request, gsr_result)
         if (gsr_result%blocked) error stop 'TDDFT production driver: Goldstone sum-rule interaction is blocked by its service contract'
         allocate(interaction(ndim, ndim))
         interaction = gsr_result%canonical_interaction
         result%interaction_provenance = 'GSR-01 independent Goldstone sum-rule interaction'
         result%goldstone_correction_status = 'not selected'
      case default
         error stop 'TDDFT production driver: interaction route was not validated'
      end select

      do iq = 1, size(config%q_list, 2)
         if (trim(config%backend) == tddft_driver_backend_native_rsgf) then
            call g_logger%info('TDRUN-02 native RSGF sweep: ndim='//int2str(ndim)//' radial_points='// &
               int2str(response_space%npoint)//' integration_points='//int2str(config%gf_integration_points), __FILE__, __LINE__)
         end if
         call evaluate_bare_response(config, response_space, radial_bases, left_state, endpoints(iq), &
                                     config%q_list(:, iq), config%frequencies, bare_result, native_provider, native_pairs, &
                                     native_site_positions)
         if (config%reciprocal_backend_crosscheck) then
            call evaluate_reciprocal_backend_crosscheck(config, response_space, radial_bases, left_state, endpoints(iq), &
               config%q_list(:, iq), config%frequencies, bare_result, result%reciprocal_crosscheck_valid(:, iq), &
               result%reciprocal_crosscheck_norm_lehmann(:, iq), result%reciprocal_crosscheck_norm_gf(:, iq), &
               result%reciprocal_crosscheck_difference_frobenius(:, iq), result%reciprocal_crosscheck_relative_frobenius(:, iq), &
               result%reciprocal_crosscheck_difference_infinity(:, iq), result%reciprocal_crosscheck_delta(:, :, :, iq))
         end if
         result%bare_response_provenance = bare_result%response_space_metadata
         dyson_request%response_space => response_space
         dyson_request%q = config%q_list(:, iq)
         dyson_request%frequencies = config%frequencies
         dyson_request%eta = config%eta
         dyson_request%channel = config%channel
         dyson_request%ks_susceptibility = bare_result%susceptibility
         dyson_request%canonical_interaction = interaction
         dyson_request%interaction_route = trim(config%interaction_route)
         dyson_request%interaction_provenance = result%interaction_provenance
         dyson_request%electronic_state_provenance = 'accepted reciprocal eigenpair snapshot; occupations and EF fixed at SCF handoff'
         dyson_request%response_space_metadata = bare_result%response_space_metadata
         call evaluate_tddft_dyson(dyson_request, dyson_result)
         result%ks_susceptibility(:, :, :, iq) = dyson_result%ks_susceptibility
         result%enhanced_susceptibility(:, :, :, iq) = dyson_result%enhanced_susceptibility
         result%loss_matrix(:, :, :, iq) = dyson_result%loss_matrix
         result%status = trim(dyson_result%status)
      end do
      result%initialized = .true.
   end subroutine evaluate_tddft_production_sweep

   subroutine validate_prepared_sweep(config, response_space, radial_bases, ground_states, left_state, endpoints, &
                                      native_provider, native_pairs)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoints(:)
      class(lr_rs_gf_provider), intent(in), optional :: native_provider
      type(lr_rs_gf_pair), intent(in), optional :: native_pairs(:)
      integer :: iq

      call validate_tddft_config(config)
      if (.not. config%enabled) error stop 'TDDFT production driver: disabled configuration cannot run a sweep'
      if (size(ground_states) /= response_space%nsite .or. size(radial_bases) /= response_space%nsite) then
         error stop 'TDDFT production driver: response-site state shape mismatch'
      end if
      if (size(endpoints) /= size(config%q_list, 2)) error stop 'TDDFT production driver: q endpoint count mismatch'
      call left_state%validate('TDDFT production driver:left_state')
      do iq = 1, size(endpoints)
         call endpoints(iq)%validate('TDDFT production driver:q_endpoint_state')
      end do
      if (trim(config%backend) == tddft_driver_backend_native_rsgf) then
         if (.not. present(native_provider) .or. .not. present(native_pairs)) then
            error stop 'TDDFT production driver: native_rsgf requires a registered provider and complete pair set'
         end if
         if (size(native_pairs) < 1) error stop 'TDDFT production driver: native_rsgf pair set is empty'
         if (native_provider%nbasis_per_site < 1) then
            error stop 'TDDFT production driver: native_rsgf provider is not initialized'
         end if
      end if
   end subroutine validate_prepared_sweep

   subroutine evaluate_bare_response(config, response_space, radial_bases, left_state, endpoint, q, frequencies, result, &
                                     native_provider, native_pairs, native_site_positions)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(lr_electronic_state), target, intent(in) :: left_state
      type(lr_electronic_state), target, intent(in) :: endpoint
      real(rp), intent(in) :: q(3), frequencies(:)
      type(lr_ks_susceptibility_result), intent(out) :: result
      class(lr_rs_gf_provider), target, intent(inout), optional :: native_provider
      type(lr_rs_gf_pair), intent(in), optional :: native_pairs(:)
      real(rp), target, intent(in), optional :: native_site_positions(:, :)
      type(lr_rs_gf_susceptibility_request) :: native_request

      if (trim(config%backend) == tddft_driver_backend_lehmann .or. trim(config%backend) == 'spectral') then
         call evaluate_lehmann_backend(response_space, radial_bases, left_state, endpoint, q, frequencies, config%eta, &
            config%channel, result)
      else if (trim(config%backend) == tddft_driver_backend_reciprocal_gf) then
         call evaluate_reciprocal_gf_backend(config, response_space, radial_bases, left_state, endpoint, q, frequencies, result)
      else if (trim(config%backend) == tddft_driver_backend_native_rsgf) then
         if (.not. present(native_provider) .or. .not. present(native_pairs)) then
            error stop 'TDDFT production driver: native_rsgf request lacks its registered provider/pair set'
         end if
         native_request%q = q
         native_request%frequencies = frequencies
         native_request%eta = config%eta
         native_request%channel = config%channel
         native_request%integration_points = config%gf_integration_points
         native_request%integration_eta = config%gf_integration_eta
         native_request%energy_margin = config%gf_energy_margin
         native_request%energy_min = minval(left_state%eigenvalues)
         native_request%energy_max = maxval(left_state%eigenvalues)
         native_request%fermi_level = left_state%fermi_level
         native_request%temperature = left_state%temperature
         native_request%response_space => response_space
         native_request%radial_bases => radial_bases
         native_request%provider => native_provider
         native_request%pairs = native_pairs
         if (present(native_site_positions)) native_request%site_positions => native_site_positions
         call evaluate_lr_rs_gf_susceptibility(native_request, result)
      else
         error stop 'TDDFT production driver: backend was not validated'
      end if
   end subroutine evaluate_bare_response

   subroutine evaluate_lehmann_backend(response_space, radial_bases, left_state, endpoint, q, frequencies, eta, channel, result)
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(lr_electronic_state), target, intent(in) :: left_state, endpoint
      real(rp), intent(in) :: q(3), frequencies(:), eta
      character(len=*), intent(in) :: channel
      type(lr_ks_susceptibility_result), intent(out) :: result
      type(lr_ks_susceptibility_request) :: request

      request%q = q
      request%frequencies = frequencies
      request%eta = eta
      request%channel = channel
      request%response_space => response_space
      request%radial_bases => radial_bases
      request%electronic_state => left_state
      request%q_endpoint_state => endpoint
      call evaluate_lr_ks_susceptibility(request, result)
   end subroutine evaluate_lehmann_backend

   subroutine evaluate_reciprocal_gf_backend(config, response_space, radial_bases, left_state, endpoint, q, frequencies, result)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(lr_electronic_state), target, intent(in) :: left_state, endpoint
      real(rp), intent(in) :: q(3), frequencies(:)
      type(lr_ks_susceptibility_result), intent(out) :: result
      type(lr_gf_susceptibility_request) :: request

      request%q = q
      request%frequencies = frequencies
      request%eta = config%eta
      request%channel = config%channel
      request%integration_points = config%gf_integration_points
      request%integration_eta = config%gf_integration_eta
      request%energy_margin = config%gf_energy_margin
      request%response_space => response_space
      request%radial_bases => radial_bases
      request%electronic_state => left_state
      request%q_endpoint_state => endpoint
      call evaluate_lr_gf_susceptibility(request, result)
   end subroutine evaluate_reciprocal_gf_backend

   subroutine evaluate_reciprocal_backend_crosscheck(config, response_space, radial_bases, left_state, endpoint, q, frequencies, &
                                                     selected_result, valid, norm_lehmann, norm_gf, difference_frobenius, &
                                                     relative_frobenius, difference_infinity, delta)
      type(tddft_production_config), intent(in) :: config
      type(response_space_layout), target, intent(in) :: response_space
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)
      type(lr_electronic_state), target, intent(in) :: left_state, endpoint
      real(rp), intent(in) :: q(3), frequencies(:)
      type(lr_ks_susceptibility_result), intent(in) :: selected_result
      logical, intent(out) :: valid(:)
      real(rp), intent(out) :: norm_lehmann(:), norm_gf(:), difference_frobenius(:), relative_frobenius(:), difference_infinity(:)
      complex(rp), intent(out) :: delta(:, :, :)
      type(lr_ks_susceptibility_result) :: lehmann_result, gf_result
      integer :: ifrequency

      if (trim(config%backend) == tddft_driver_backend_lehmann .or. trim(config%backend) == 'spectral') then
         lehmann_result = selected_result
         call evaluate_reciprocal_gf_backend(config, response_space, radial_bases, left_state, endpoint, q, frequencies, gf_result)
      else if (trim(config%backend) == tddft_driver_backend_reciprocal_gf) then
         call evaluate_lehmann_backend(response_space, radial_bases, left_state, endpoint, q, frequencies, config%eta, &
            config%channel, lehmann_result)
         gf_result = selected_result
      else
         call evaluate_lehmann_backend(response_space, radial_bases, left_state, endpoint, q, frequencies, config%eta, &
            config%channel, lehmann_result)
         call evaluate_reciprocal_gf_backend(config, response_space, radial_bases, left_state, endpoint, q, frequencies, gf_result)
      end if

      if (size(valid) /= size(frequencies) .or. size(norm_lehmann) /= size(frequencies) .or. &
          size(norm_gf) /= size(frequencies) .or. size(difference_frobenius) /= size(frequencies) .or. &
          size(relative_frobenius) /= size(frequencies) .or. size(difference_infinity) /= size(frequencies) .or. &
          any(shape(delta) /= [response_space%ndim, response_space%ndim, size(frequencies)])) then
         error stop 'TDDFT production driver: crosscheck output shape mismatch'
      end if

      do ifrequency = 1, size(frequencies)
         delta(:, :, ifrequency) = lehmann_result%susceptibility(:, :, ifrequency) - gf_result%susceptibility(:, :, ifrequency)
         difference_frobenius(ifrequency) = sqrt(sum(abs(delta(:, :, ifrequency))**2))
         norm_lehmann(ifrequency) = sqrt(sum(abs(lehmann_result%susceptibility(:, :, ifrequency))**2))
         norm_gf(ifrequency) = sqrt(sum(abs(gf_result%susceptibility(:, :, ifrequency))**2))
         relative_frobenius(ifrequency) = difference_frobenius(ifrequency)/ &
            max(norm_lehmann(ifrequency), norm_gf(ifrequency), tiny(1.0_rp))
         difference_infinity(ifrequency) = maxval(abs(delta(:, :, ifrequency)))
      end do
      valid = .true.
   end subroutine evaluate_reciprocal_backend_crosscheck

   subroutine radial_basis_from_snapshot(state, basis)
      type(radial_ground_state), intent(in) :: state
      type(lmto_radial_basis), intent(out) :: basis

      if (.not. allocated(state%r) .or. .not. allocated(state%pauli_large) .or. &
          .not. allocated(state%pauli_large_dot) .or. .not. allocated(state%pauli_enu)) then
         error stop 'TDDFT production driver: accepted radial snapshot lacks Pauli basis arrays'
      end if
      call basis%initialize(size(state%r), state%pauli_lmax, 2)
      basis%rofi = state%r
      basis%mesh_a = state%a
      basis%mesh_b = state%b
      basis%phi_large = state%pauli_large
      basis%phidot_large = state%pauli_large_dot
      basis%enu_radial = state%pauli_enu
      basis%enu_work = state%pauli_enu
      basis%channel_present = .true.
   end subroutine radial_basis_from_snapshot

   pure function effective_hamiltonian_order(reciprocal_obj, hamiltonian_obj) result(order)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj
      character(len=16) :: order

      order = trim(reciprocal_obj%kspace_ham_order)
      if (order == 'auto') then
         if (hamiltonian_obj%hoh) then
            order = 'second'
         else
            order = 'first'
         end if
      end if
   end function effective_hamiltonian_order

   integer function find_gamma_q(q_list) result(index_gamma)
      real(rp), intent(in) :: q_list(:, :)
      integer :: iq
      index_gamma = 0
      do iq = 1, size(q_list, 2)
         if (sum(abs(q_list(:, iq))) <= 1.0e-12_rp) then
            index_gamma = iq
            return
         end if
      end do
      error stop 'TDDFT production driver: gamma q was not supplied for the requested static route'
   end function find_gamma_q

   subroutine write_tddft_production_output(config, result, control_obj, lattice_obj, ground_states, response_space, reciprocal_obj)
      type(tddft_production_config), intent(in) :: config
      type(tddft_production_result), intent(in) :: result
      type(control), intent(in) :: control_obj
      type(lattice), intent(in) :: lattice_obj
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(response_space_layout), intent(in) :: response_space
      type(reciprocal), intent(in) :: reciprocal_obj
      integer :: unit, iq, iw, i, j
      real(rp) :: magnetic_moment
      complex(rp) :: loss_trace, ks_trace, enhanced_trace

      if (rank /= 0) return
      magnetic_moment = 0.0_rp
      do i = 1, size(ground_states)
         magnetic_moment = magnetic_moment + ground_states(i)%integrated_moment_muB
      end do
      open(newunit=unit, file=trim(config%output_file), status='replace', action='write')
      write(unit, '(a)') '# TDDFT production response; matrices are LR-04 canonical service outputs'
      write(unit, '(a,a)') '# git_version = ', trim(tddft_build_version)
      write(unit, '(a,a)') '# structure_calctype = ', control_obj%calctype
      write(unit, '(a,es24.16)') '# structure_alat = ', lattice_obj%alat
      write(unit, '(a,i0)') '# structure_ntype = ', lattice_obj%ntype
      write(unit, '(a,i0)') '# structure_nrec = ', lattice_obj%nrec
      write(unit, '(a,a)') '# structure_symbols = ', trim(lattice_obj%symbolic_atoms(lattice_obj%nbulk + 1)%element%symbol)
      write(unit, '(a,a)') '# xc_identity = ', trim(ground_states(1)%xc_provenance%functional_name)
      write(unit, '(a,a)') '# xc_backend = ', trim(ground_states(1)%xc_provenance%backend_name)
      write(unit, '(a,i0)') '# xc_txc = ', ground_states(1)%xc_provenance%txc
      write(unit, '(a,es24.16)') '# magnetic_moment_muB = ', magnetic_moment
      write(unit, '(a,es24.16)') '# fermi_level_Ry = ', reciprocal_obj%fermi_level
      write(unit, '(a,es24.16)') '# temperature_K = ', reciprocal_obj%temperature
      write(unit, '(a)') '# response_capability = collinear,no_soc,ham_only,orthogonal,sp/spd,no_extra_operator'
      write(unit, '(a,a)') '# q = ', 'one row block per requested q; values are reduced coordinates'
      write(unit, '(a,a)') '# omega_Ry = ', 'one row per requested frequency'
      write(unit, '(a,es24.16)') '# eta_Ry = ', config%eta
      write(unit, '(a,i0)') '# response_angular_cutoff = ', result%response_lmax
      write(unit, '(a,i0,2(es24.16,1x),a,es24.16,a,es24.16)') '# radial_mesh_identity = ', size(ground_states(1)%r), &
         ground_states(1)%a, ground_states(1)%b, 'rmax=', ground_states(1)%rmax, 'sum_r=', sum(ground_states(1)%r)
      write(unit, '(a,a)') '# interaction_route = ', trim(config%interaction_route)
      write(unit, '(a,a)') '# interaction_provenance = ', trim(result%interaction_provenance)
      write(unit, '(a,a)') '# goldstone_correction = ', trim(result%goldstone_correction_status)
      write(unit, '(a,a)') '# backend = ', trim(config%backend)
      if (config%reciprocal_backend_crosscheck) then
         write(unit, '(a)') '# reciprocal_backend_crosscheck = enabled (validation diagnostic; selected backend remains authoritative)'
         write(unit, '(a)') '# crosscheck_metrics columns: q_index omega_Ry valid norm_lehmann norm_reciprocal_gf difference_frobenius relative_frobenius difference_infinity'
         do iq = 1, result%nq
            do iw = 1, result%nfrequency
               write(unit, '(a,i0,1x,es24.16,1x,l1,5(1x,es24.16))') '# crosscheck_metrics ', iq, result%frequencies(iw), &
                  result%reciprocal_crosscheck_valid(iw, iq), result%reciprocal_crosscheck_norm_lehmann(iw, iq), &
                  result%reciprocal_crosscheck_norm_gf(iw, iq), result%reciprocal_crosscheck_difference_frobenius(iw, iq), &
                  result%reciprocal_crosscheck_relative_frobenius(iw, iq), result%reciprocal_crosscheck_difference_infinity(iw, iq)
            end do
         end do
         write(unit, '(a)') '# crosscheck_delta columns: q_index omega_Ry matrix_i matrix_j delta_real delta_imag'
         do iq = 1, result%nq
            do iw = 1, result%nfrequency
               do j = 1, result%ndim
                  do i = 1, result%ndim
                     write(unit, '(a,i0,1x,es24.16,1x,2(i0,1x),2(es24.16,1x))') '# crosscheck_delta ', iq, &
                        result%frequencies(iw), i, j, real(result%reciprocal_crosscheck_delta(i, j, iw, iq), rp), &
                        aimag(result%reciprocal_crosscheck_delta(i, j, iw, iq))
                  end do
               end do
            end do
         end do
      end if
      if (trim(config%backend) == tddft_driver_backend_native_rsgf) then
         write(unit, '(a,a)') '# native_rsgf_provider = ', trim(config%native_rsgf_provider)
         write(unit, '(a,a)') '# native_rsgf_provenance = ', trim(result%bare_response_provenance)
         write(unit, '(a)') '# native_rsgf_route = coefficient-GF -> RSGF endpoint augmentation -> LR-04 canonical -> common KXC/Dyson'
      end if
      if (config%write_full_matrix) then
         write(unit, '(a)') '# output_matrix_mode = full'
         write(unit, '(a)') '# columns: q_index omega_Ry matrix_i matrix_j loss_real loss_imag chiKS_real chiKS_imag enhanced_real enhanced_imag'
      else
         write(unit, '(a)') '# output_matrix_mode = trace (complete matrices retained by the driver result)'
         write(unit, '(a)') '# columns: q_index omega_Ry loss_trace_real loss_trace_imag chiKS_trace_real chiKS_trace_imag enhanced_trace_real enhanced_trace_imag'
      end if
      do iq = 1, result%nq
         write(unit, '(a,i0,3(1x,es24.16))') '# q_index = ', iq, result%q_list(:, iq)
         do iw = 1, result%nfrequency
            if (config%write_full_matrix) then
               do j = 1, result%ndim
                  do i = 1, result%ndim
                     write(unit, '(i0,1x,es24.16,1x,2(i0,1x),8(es24.16,1x))') iq, result%frequencies(iw), i, j, &
                        real(result%loss_matrix(i, j, iw, iq), rp), aimag(result%loss_matrix(i, j, iw, iq)), &
                        real(result%ks_susceptibility(i, j, iw, iq), rp), aimag(result%ks_susceptibility(i, j, iw, iq)), &
                        real(result%enhanced_susceptibility(i, j, iw, iq), rp), aimag(result%enhanced_susceptibility(i, j, iw, iq))
                  end do
               end do
            else
               loss_trace = response_operator_trace(response_space, result%loss_matrix(:, :, iw, iq))
               ks_trace = response_operator_trace(response_space, result%ks_susceptibility(:, :, iw, iq))
               enhanced_trace = response_operator_trace(response_space, result%enhanced_susceptibility(:, :, iw, iq))
               write(unit, '(i0,1x,es24.16,1x,6(es24.16,1x))') iq, result%frequencies(iw), &
                  real(loss_trace, rp), aimag(loss_trace), real(ks_trace, rp), aimag(ks_trace), &
                  real(enhanced_trace, rp), aimag(enhanced_trace)
            end if
         end do
      end do
      close(unit)
   end subroutine write_tddft_production_output

end module tddft_production_driver_mod
