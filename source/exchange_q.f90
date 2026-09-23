!------------------------------------------------------------------------------
! Static finite-q exchange-curvature post-processing.
!
! This module deliberately owns the user-facing exchange_q input boundary and
! the production driver only.  The finite-H algebra remains in
! lr_kl_hessian_mod.  The native Turek contour evaluator is a separate
! path-operator route and does not consume finite-H resolvents or eigenvectors.
!------------------------------------------------------------------------------
module exchange_q_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use basis_mod, only: norb
   use control_mod, only: control
   use energy_mod, only: energy
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use lr_kl_hessian_mod, only: lmto_fixture_adapter_residual
   use lr_kl_contour_mod, only: finite_h_contour_options, finite_h_contour_report, &
      force_theorem_finite_q_hessian_from_resolvent_batch
   use lr_lmto_turek_contour_mod, only: native_turek_contour_options, native_turek_contour_report, &
      native_turek_static_reference
   use logger_mod, only: g_logger
   use lr_rotation_response_mod, only: rotation_state, rotation_request, rotation_result, &
      prepare_rotation_response, evaluate_rotation_response, reduce_static_rotation_kernel
   use lr_kl_hessian_mod, only: lmto_live_hamiltonian_fixture, lmto_fixture_from_hamiltonian, &
      assemble_lmto_hamiltonian, &
      assemble_lmto_finite_q_torques, assemble_lmto_finite_q_mixed_derivative, &
      force_theorem_finite_q_hessian_from_eigenbasis_batch, &
      force_theorem_finite_q_hessian_from_eigenbasis_metallic_batch
   use math_mod, only: ang2au, inverse_3x3, pi
   use precision_mod, only: rp
   use reciprocal_mod, only: reciprocal
   use self_mod, only: self
   use string_mod, only: sl, int2str, lower
   implicit none

   private

   integer, parameter, public :: exchange_q_max_points = 2000
   real(rp), parameter :: q_zero_tolerance = 1.0e-12_rp
   real(rp), parameter :: reciprocal_kb_ry_per_k = 6.3336814e-6_rp

   type, public :: exchange_q_config
      character(len=sl) :: fname = ''
      character(len=16) :: q_coordinates = 'direct'
      character(len=sl) :: q_file = ''
      integer :: n_q = 0
      real(rp), allocatable :: q_list(:, :)
      character(len=sl) :: output_file = 'exchange_q.dat'
      logical :: write_components = .true.
      logical :: rotation_dynamics = .false.
      character(len=sl) :: rotation_output_file = 'rotation_dynamics.dat'
      logical :: native_crosscheck = .false.
      logical :: native_turek = .false.
      character(len=24) :: finite_h_spectral_mode = 'metallic'
      character(len=16) :: finite_h_response_backend = 'spectral'
      real(rp) :: rotation_axis(3) = [1.0_rp, 0.0_rp, 0.0_rp]
      integer :: contour_points = 32
      character(len=16) :: contour_shape = 'ellipse'
      real(rp) :: contour_margin = 0.25_rp
      real(rp) :: contour_height_fraction = 0.35_rp
      logical :: contour_account_fermi_poles = .true.
      ! Deprecated aliases retained for old decks.  Native evaluation now
      ! uses the contour controls below; these no longer select a Simpson path.
      real(rp) :: native_green_eta = 1.0e-3_rp
      integer :: native_energy_points = 0
      integer :: native_contour_points = 64
      real(rp) :: native_contour_margin = 0.25_rp
      real(rp) :: native_contour_height_fraction = 0.35_rp
      logical :: native_contour_account_fermi_poles = .true.
      integer :: native_contour_target_fermi_poles = 0
   contains
      procedure :: clear => exchange_q_config_clear
   end type exchange_q_config

   public :: load_exchange_q_config
   public :: exchange_q_convert_points
   public :: detect_commensurate_endpoint_map
   public :: run_exchange_q
   public :: lkag_pauli_product

contains

   subroutine exchange_q_config_clear(this)
      class(exchange_q_config), intent(inout) :: this
      if (allocated(this%q_list)) deallocate(this%q_list)
      this%fname = ''
      this%q_coordinates = 'direct'
      this%q_file = ''
      this%n_q = 0
      this%output_file = 'exchange_q.dat'
      this%write_components = .true.
      this%rotation_dynamics = .false.
      this%rotation_output_file = 'rotation_dynamics.dat'
      this%native_crosscheck = .false.
      this%native_turek = .false.
      this%finite_h_spectral_mode = 'metallic'
      this%finite_h_response_backend = 'spectral'
      this%rotation_axis = [1.0_rp, 0.0_rp, 0.0_rp]
      this%contour_points = 32
      this%contour_shape = 'ellipse'
      this%contour_margin = 0.25_rp
      this%contour_height_fraction = 0.35_rp
      this%contour_account_fermi_poles = .true.
      this%native_green_eta = 1.0e-3_rp
      this%native_energy_points = 0
      this%native_contour_points = 64
      this%native_contour_margin = 0.25_rp
      this%native_contour_height_fraction = 0.35_rp
      this%native_contour_account_fermi_poles = .true.
      this%native_contour_target_fermi_poles = 0
   end subroutine exchange_q_config_clear

   subroutine load_exchange_q_config(filename, this, validate_request)
      character(len=*), intent(in) :: filename
      type(exchange_q_config), intent(inout) :: this
      logical, intent(in), optional :: validate_request
      logical :: validate
      integer :: iostatus, funit, n_q_points
      character(len=16) :: q_coordinates
      character(len=sl) :: q_file, output_file, rotation_output_file
      logical :: write_components, rotation_dynamics, native_crosscheck, native_turek, native_contour_account_fermi_poles
      character(len=24) :: finite_h_spectral_mode
      character(len=16) :: finite_h_response_backend, contour_shape
      real(rp) :: rotation_axis(3), native_green_eta
      real(rp) :: contour_margin, contour_height_fraction, native_contour_margin, native_contour_height_fraction
      integer :: native_energy_points, contour_points, native_contour_points, native_contour_target_fermi_poles
      logical :: contour_account_fermi_poles
      real(rp) :: q_list(3, exchange_q_max_points)

      namelist /exchange_q/ q_coordinates, q_file, n_q_points, q_list, output_file, &
         write_components, rotation_dynamics, rotation_output_file, native_crosscheck, native_turek, &
         finite_h_spectral_mode, finite_h_response_backend, rotation_axis, &
         contour_points, contour_shape, contour_margin, contour_height_fraction, contour_account_fermi_poles, &
         native_green_eta, native_energy_points, native_contour_points, native_contour_margin, &
         native_contour_height_fraction, native_contour_account_fermi_poles, native_contour_target_fermi_poles

      call this%clear()
      this%fname = filename
      q_coordinates = this%q_coordinates
      q_file = this%q_file
      n_q_points = 0
      q_list = 0.0_rp
      output_file = this%output_file
      write_components = this%write_components
      rotation_dynamics = this%rotation_dynamics
      rotation_output_file = this%rotation_output_file
      native_crosscheck = this%native_crosscheck
      native_turek = this%native_turek
      finite_h_spectral_mode = this%finite_h_spectral_mode
      finite_h_response_backend = this%finite_h_response_backend
      rotation_axis = this%rotation_axis
      contour_points = this%contour_points
      contour_shape = this%contour_shape
      contour_margin = this%contour_margin
      contour_height_fraction = this%contour_height_fraction
      contour_account_fermi_poles = this%contour_account_fermi_poles
      native_green_eta = this%native_green_eta
      native_energy_points = this%native_energy_points
      native_contour_points = this%native_contour_points
      native_contour_margin = this%native_contour_margin
      native_contour_height_fraction = this%native_contour_height_fraction
      native_contour_account_fermi_poles = this%native_contour_account_fermi_poles
      native_contour_target_fermi_poles = this%native_contour_target_fermi_poles

      open(newunit=funit, file=filename, action='read', status='old', iostat=iostatus)
      if (iostatus /= 0) then
         call g_logger%fatal('[exchange_q]: input file '//trim(filename)//' not found', __FILE__, __LINE__)
      end if
      read(funit, nml=exchange_q, iostat=iostatus)
      close(funit)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%fatal('[exchange_q]: error reading &exchange_q, iostatus='//int2str(iostatus), __FILE__, __LINE__)
      end if

      this%q_coordinates = lower(trim(q_coordinates))
      if (len_trim(this%q_coordinates) == 0) this%q_coordinates = 'direct'
      if (this%q_coordinates /= 'direct' .and. this%q_coordinates /= 'cartesian') then
         call g_logger%fatal("[exchange_q]: q_coordinates must be 'direct' or 'cartesian'", __FILE__, __LINE__)
      end if
      this%q_file = trim(q_file)
      this%output_file = trim(output_file)
      this%write_components = write_components
      this%rotation_dynamics = rotation_dynamics
      this%rotation_output_file = trim(rotation_output_file)
      this%native_crosscheck = native_crosscheck
      this%native_turek = native_turek
      this%finite_h_spectral_mode = lower(trim(finite_h_spectral_mode))
      this%finite_h_response_backend = lower(trim(finite_h_response_backend))
      this%rotation_axis = rotation_axis
      this%contour_points = contour_points
      this%contour_shape = lower(trim(contour_shape))
      this%contour_margin = contour_margin
      this%contour_height_fraction = contour_height_fraction
      this%contour_account_fermi_poles = contour_account_fermi_poles
      this%native_green_eta = native_green_eta
      this%native_energy_points = native_energy_points
      this%native_contour_points = native_contour_points
      this%native_contour_margin = native_contour_margin
      this%native_contour_height_fraction = native_contour_height_fraction
      this%native_contour_account_fermi_poles = native_contour_account_fermi_poles
      this%native_contour_target_fermi_poles = native_contour_target_fermi_poles
      validate = .false.
      if (present(validate_request)) validate = validate_request
      if (validate) then
         if (len_trim(this%output_file) == 0) then
            call g_logger%fatal('[exchange_q]: output_file must not be blank', __FILE__, __LINE__)
         end if
         if (this%rotation_dynamics .and. len_trim(this%rotation_output_file) == 0) then
            call g_logger%fatal('[exchange_q]: rotation_output_file must not be blank', __FILE__, __LINE__)
         end if
         if (this%native_green_eta <= 0.0_rp) then
            call g_logger%fatal('[exchange_q]: native_green_eta must be positive', __FILE__, __LINE__)
         end if
         if (this%native_energy_points < 0) then
            call g_logger%fatal('[exchange_q]: native_energy_points must be non-negative', __FILE__, __LINE__)
         end if
         if (this%native_contour_points < 8) then
            call g_logger%fatal('[exchange_q]: native_contour_points must be at least 8', __FILE__, __LINE__)
         end if
         if (this%native_contour_margin <= 0.0_rp .or. this%native_contour_height_fraction <= 0.0_rp) then
            call g_logger%fatal('[exchange_q]: native contour margin and height fraction must be positive', __FILE__, __LINE__)
         end if
         if (this%native_contour_target_fermi_poles < 0 .or. &
             (this%native_contour_target_fermi_poles > 0 .and. mod(this%native_contour_target_fermi_poles,2) /= 0)) then
            call g_logger%fatal('[exchange_q]: native_contour_target_fermi_poles must be zero or positive even', __FILE__, __LINE__)
         end if
         if (maxval(abs(this%rotation_axis)) <= tiny(1.0_rp)) then
            call g_logger%fatal('[exchange_q]: rotation_axis must be nonzero', __FILE__, __LINE__)
         end if
         if (this%finite_h_spectral_mode /= 'metallic' .and. this%finite_h_spectral_mode /= 'legacy_occupied') then
            call g_logger%fatal("[exchange_q]: finite_h_spectral_mode must be 'metallic' or 'legacy_occupied'", &
               __FILE__, __LINE__)
         end if
         if (this%finite_h_response_backend /= 'spectral' .and. this%finite_h_response_backend /= 'contour' .and. &
             this%finite_h_response_backend /= 'both') then
            call g_logger%fatal("[exchange_q]: finite_h_response_backend must be 'spectral', 'contour', or 'both'", &
               __FILE__, __LINE__)
         end if
         if (this%contour_points < 8) call g_logger%fatal('[exchange_q]: contour_points must be at least 8', __FILE__, __LINE__)
         if (this%contour_shape /= 'ellipse') call g_logger%fatal("[exchange_q]: contour_shape must be 'ellipse'", __FILE__, __LINE__)
         if (this%contour_margin <= 0.0_rp .or. this%contour_height_fraction <= 0.0_rp) then
            call g_logger%fatal('[exchange_q]: contour margin and height fraction must be positive', __FILE__, __LINE__)
         end if
      end if

      if (len_trim(this%q_file) > 0) then
         call read_exchange_q_file(this, this%q_file)
      else
         if (n_q_points < 1 .or. n_q_points > exchange_q_max_points) then
            if (validate) call g_logger%fatal('[exchange_q]: n_q_points must be in [1,'// &
               int2str(exchange_q_max_points)//']', __FILE__, __LINE__)
            return
         end if
         this%n_q = n_q_points
         allocate(this%q_list(3, this%n_q))
         this%q_list = q_list(:, 1:this%n_q)
      end if
      if (validate .and. this%n_q < 2) then
         call g_logger%fatal('[exchange_q]: provide at least q=0 and one finite-q point', __FILE__, __LINE__)
      end if
      if (validate .and. any(.not. ieee_is_finite(this%q_list))) then
         call g_logger%fatal('[exchange_q]: q path contains a non-finite coordinate', __FILE__, __LINE__)
      end if
   end subroutine load_exchange_q_config

   subroutine read_exchange_q_file(this, filename)
      type(exchange_q_config), intent(inout) :: this
      character(len=*), intent(in) :: filename
      integer :: funit, ios, iq, n_q

      open(newunit=funit, file=trim(filename), action='read', status='old', iostat=ios)
      if (ios /= 0) call g_logger%fatal('[exchange_q]: q_file '//trim(filename)//' not found', __FILE__, __LINE__)
      read(funit, *, iostat=ios) n_q
      if (ios /= 0 .or. n_q < 1 .or. n_q > exchange_q_max_points) then
         close(funit)
         call g_logger%fatal('[exchange_q]: malformed q_file count in '//trim(filename), __FILE__, __LINE__)
      end if
      allocate(this%q_list(3, n_q))
      this%n_q = n_q
      do iq = 1, n_q
         read(funit, *, iostat=ios) this%q_list(:, iq)
         if (ios /= 0) then
            close(funit)
            call g_logger%fatal('[exchange_q]: malformed q_file row '//int2str(iq)//' in '//trim(filename), &
               __FILE__, __LINE__)
         end if
      end do
      close(funit)
   end subroutine read_exchange_q_file

   ! Convert the input path to reciprocal-direct coordinates consumed by the
   ! finite-q API.  Cartesian input follows frozen_magnon: values are in units
   ! of 2*pi/alat.  The output Cartesian path is in those same units.
   subroutine exchange_q_convert_points(config, lattice_obj, q_direct, q_cart)
      type(exchange_q_config), intent(in) :: config
      type(lattice), intent(in) :: lattice_obj
      real(rp), allocatable, intent(out) :: q_direct(:, :), q_cart(:, :)
      real(rp) :: direct_to_cart(3,3)

      if (.not. allocated(config%q_list)) error stop 'exchange_q_convert_points: q path is not allocated'
      allocate(q_direct(3, config%n_q), q_cart(3, config%n_q))
      direct_to_cart = transpose(inverse_3x3(lattice_obj%a))
      if (config%q_coordinates == 'direct') then
         q_direct = config%q_list
         q_cart = matmul(direct_to_cart, q_direct)
      else
         q_cart = config%q_list
         q_direct = matmul(transpose(lattice_obj%a), q_cart)
      end if
   end subroutine exchange_q_convert_points

   !> Detect and map an exact q translation on a full uniform reciprocal mesh.
   !> The map is constructed from the stored mesh coordinates, not from an
   !> assumed flat ordering.  This makes it valid for both the internal MP
   !> generator and a full mesh returned by spglib.  No tolerance is used to
   !> round an incommensurate q: the routine returns commensurate=.false. in
   !> that case and the caller must use explicit endpoint diagonalization.
   subroutine detect_commensurate_endpoint_map(k_points, k_weights, nk_mesh, q_point, endpoint_index, commensurate, residual)
      real(rp), intent(in) :: k_points(:, :), k_weights(:), q_point(3)
      integer, intent(in) :: nk_mesh(3)
      integer, intent(out) :: endpoint_index(:)
      logical, intent(out) :: commensurate
      real(rp), intent(out) :: residual
      integer, allocatable :: grid_index(:, :, :)
      integer :: nk, ix, iy, iz, ik, target, qcell(3), cell(3), flat
      integer :: n1, n2, n3, nint_coord
      real(rp) :: origin(3), folded(3), delta, qscaled
      real(rp), parameter :: mesh_tolerance = 2.0e-10_rp

      nk = size(k_points,2); n1 = nk_mesh(1); n2 = nk_mesh(2); n3 = nk_mesh(3)
      commensurate = .false.; residual = huge(1.0_rp)
      if (size(k_points,1) /= 3 .or. size(k_weights) /= nk .or. size(endpoint_index) /= nk .or. &
          n1 < 1 .or. n2 < 1 .or. n3 < 1 .or. nk /= n1*n2*n3) return

      ! Gamma is always a safe identity map, including a symmetry-reduced
      ! mesh.  Nonzero reuse below deliberately requires uniform weights and
      ! a full mesh.
      if (sqrt(sum(q_point**2)) <= q_zero_tolerance) then
         do ik = 1, nk
            endpoint_index(ik) = ik
         end do
         commensurate = .true.; residual = 0.0_rp
         return
      end if
      if (nk > 1 .and. maxval(abs(k_weights-k_weights(1))) > mesh_tolerance*max(1.0_rp,abs(k_weights(1)))) return

      origin = fold_direct_point(k_points(:,1))
      allocate(grid_index(n1,n2,n3)); grid_index = 0
      do ik = 1, nk
         folded = fold_direct_point(k_points(:,ik))
         do ix = 1, 3
            delta = real(nk_mesh(ix),rp)*(folded(ix)-origin(ix))
            nint_coord = nint(delta)
            if (abs(delta-real(nint_coord,rp)) > mesh_tolerance) then
               deallocate(grid_index)
               return
            end if
            cell(ix) = modulo(nint_coord,nk_mesh(ix)) + 1
         end do
         flat = cell(1) + (cell(2)-1)*n1 + (cell(3)-1)*n1*n2
         ix = 1 + modulo(flat-1,n1)
         iy = 1 + modulo((flat-1)/n1,n2)
         iz = 1 + (flat-1)/(n1*n2)
         if (grid_index(ix,iy,iz) /= 0) then
            deallocate(grid_index)
            return
         end if
         grid_index(ix,iy,iz) = ik
      end do
      if (any(grid_index == 0)) then
         deallocate(grid_index)
         return
      end if

      residual = 0.0_rp
      do ix = 1, 3
         qscaled = q_point(ix)*real(nk_mesh(ix),rp)
         qcell(ix) = nint(qscaled)
         residual = max(residual,abs(qscaled-real(qcell(ix),rp))/real(nk_mesh(ix),rp))
      end do
      if (residual > mesh_tolerance) then
         deallocate(grid_index)
         return
      end if

      do ik = 1, nk
         folded = fold_direct_point(k_points(:,ik))
         do ix = 1, 3
            delta = real(nk_mesh(ix),rp)*(folded(ix)-origin(ix))
            cell(ix) = modulo(nint(delta)+qcell(ix),nk_mesh(ix)) + 1
         end do
         target = grid_index(cell(1),cell(2),cell(3))
         endpoint_index(ik) = target
         ! Check the actual folded endpoint, including the boundary crossing.
         folded = fold_direct_point(k_points(:,ik)+q_point)
         residual = max(residual,maxval(abs(folded-fold_direct_point(k_points(:,target)))))
         if (abs(k_weights(ik)-k_weights(target)) > mesh_tolerance*max(1.0_rp,abs(k_weights(ik)))) then
            deallocate(grid_index)
            return
         end if
      end do
      commensurate = residual <= mesh_tolerance
      deallocate(grid_index)
   end subroutine detect_commensurate_endpoint_map

   pure function fold_direct_point(point) result(folded)
      real(rp), intent(in) :: point(3)
      real(rp) :: folded(3)
      folded = point-floor(point+0.5_rp)
   end function fold_direct_point

   ! Exact scalar kernel copied algebraically from exchange%dGdG_Jnc:
   !   d_i G_inmag d_j G_jnmag
   ! - d_i G_ix    d_j G_jx
   ! - d_i G_iy    d_j G_jy
   ! - d_i G_iz    d_j G_jz.
   ! Keeping this as a small independent contraction makes the q-space
   ! evaluator auditable without importing or modifying exchange.f90.
   pure function lkag_pauli_product(di, dj, ginmag, gix, giy, giz, gjnmag, gjx, gjy, gjz) result(value)
      complex(rp), intent(in) :: di(:, :), dj(:, :), ginmag(:, :), gix(:, :), giy(:, :), giz(:, :)
      complex(rp), intent(in) :: gjnmag(:, :), gjx(:, :), gjy(:, :), gjz(:, :)
      real(rp) :: value
      complex(rp) :: product(size(di,1), size(di,2))
      integer :: i

      product = matmul(di, matmul(ginmag, matmul(dj, gjnmag))) - &
                matmul(di, matmul(gix, matmul(dj, gjx))) - &
                matmul(di, matmul(giy, matmul(dj, gjy))) - &
                matmul(di, matmul(giz, matmul(dj, gjz)))
      value = 0.0_rp
      do i = 1, min(size(product,1), size(product,2))
         value = value + aimag(product(i,i))
      end do
   end function lkag_pauli_product

   subroutine run_exchange_q(config, control_obj, lattice_obj, hamiltonian_obj, energy_obj, self_obj, reciprocal_obj)
      type(exchange_q_config), intent(in) :: config
      type(control), intent(in) :: control_obj
      type(lattice), intent(inout) :: lattice_obj
      type(hamiltonian), intent(inout) :: hamiltonian_obj
      type(energy), intent(inout) :: energy_obj
      type(self), intent(inout) :: self_obj
      type(reciprocal), intent(inout) :: reciprocal_obj
      type(lmto_live_hamiltonian_fixture) :: fixture
      real(rp), allocatable :: q_direct(:, :), q_cart(:, :)
      real(rp), allocatable :: evals(:, :), endpoint_evals(:, :), endpoint_eval_check(:, :), weights(:)
      complex(rp), allocatable :: evecs(:, :, :), endpoint_evecs(:, :, :), endpoint_evec_check(:, :, :)
      complex(rp), allocatable :: torques_q(:, :, :, :), torques_minus_q(:, :, :, :), mixed(:, :, :, :, :)
      complex(rp), allocatable :: torque_minus_check(:, :, :)
      complex(rp), allocatable :: hessian(:, :), torque_torque(:, :), contact(:, :), complete(:, :)
      real(rp), allocatable :: axes(:, :), finite_total(:, :, :), finite_tt(:, :, :), finite_contact(:, :, :)
      real(rp), allocatable :: spectral_total(:, :, :), spectral_tt(:, :, :), spectral_contact(:, :, :)
      real(rp), allocatable :: contour_total(:, :, :), contour_tt(:, :, :), contour_contact(:, :, :)
      complex(rp), allocatable :: h_source(:, :, :), h_endpoint(:, :, :)
      real(rp), allocatable :: native_jq_ud(:), native_jq_du(:), native_jq_sym(:), native_delta_j(:), native_curvature(:)
      integer, allocatable :: endpoint_index(:)
      logical, allocatable :: endpoint_reused(:), q_commensurate(:)
      real(rp), allocatable :: endpoint_residual(:)
      character(len=24), allocatable :: endpoint_mode(:)
      real(rp) :: adapter_before, adapter_after, axis_norm, finite_h_kT
      real(rp) :: native_gamma, native_q
      real(rp) :: endpoint_seconds, assembly_seconds, contraction_seconds, vertex_error, endpoint_identity_error
      real(rp) :: hamiltonian_seconds, gf_seconds, solve_seconds, contour_seconds, total_response_seconds
      type(finite_h_contour_options) :: contour_options
      type(finite_h_contour_report) :: contour_report
      type(native_turek_contour_options) :: native_contour_options
      type(native_turek_contour_report) :: native_contour_report
      integer :: nsite, nmat, nk, iq, ik, ia, ja, clock_start, clock_end, clock_rate, response_start
      logical :: native_ready, vertex_identity_checked, endpoint_identity_checked, do_spectral, do_contour

      call validate_exchange_q_capability(config, control_obj, hamiltonian_obj, self_obj, reciprocal_obj)
      call exchange_q_convert_points(config, lattice_obj, q_direct, q_cart)
      nsite = hamiltonian_obj%charge%lattice%nrec
      axis_norm = sqrt(sum(config%rotation_axis**2))
      allocate(axes(3,nsite))
      do ia = 1, nsite
         axes(:,ia) = config%rotation_axis/axis_norm
      end do

      ! Freeze the finite-H accepted state before representation-specific
      ! native-LKAG preparation.  predls() is allowed to mutate symbolic-atom
      ! P/S channels, but must not alter the production Hamiltonian consumed by
      ! finite-H or reciprocal endpoint solves.
      call lmto_fixture_from_hamiltonian(hamiltonian_obj, fixture)
      call lmto_fixture_adapter_residual(hamiltonian_obj, fixture, reciprocal_obj%k_points(:,1:1), adapter_before)
      do ik = 1, lattice_obj%ntype
         call lattice_obj%symbolic_atoms(ik)%predls(lattice_obj%wav*ang2au)
      end do
      call lmto_fixture_adapter_residual(hamiltonian_obj, fixture, reciprocal_obj%k_points(:,1:1), adapter_after)
      call g_logger%info('[exchange_q]: accepted-H adapter residual before/after='//trim(real_to_string(adapter_before))// &
         '/'//trim(real_to_string(adapter_after)), __FILE__, __LINE__)
      if (adapter_before > 3.0e-12_rp .or. adapter_after > 3.0e-12_rp .or. &
          abs(adapter_after-adapter_before) > 1.0e-13_rp) then
         call g_logger%fatal('[exchange_q]: predls() changed the accepted production-H adapter state', __FILE__, __LINE__)
      end if
      if (abs(self_obj%reciprocal_scf_cache%fermi_level-energy_obj%fermi) > 3.0e-12_rp) then
         call g_logger%fatal('[exchange_q]: accepted Fermi-level provenance mismatch', __FILE__, __LINE__)
      end if

      nsite = fixture%nsite
      nmat = 2*fixture%norb*nsite
      nk = size(reciprocal_obj%k_points,2)
      do_spectral = config%finite_h_response_backend == 'spectral' .or. config%finite_h_response_backend == 'both'
      do_contour = config%finite_h_response_backend == 'contour' .or. config%finite_h_response_backend == 'both'
      if (nsite < 1 .or. nk < 1 .or. (do_spectral .and. .not. allocated(reciprocal_obj%eigenvectors))) then
         call g_logger%fatal('[exchange_q]: accepted reciprocal eigensystem is incomplete', __FILE__, __LINE__)
      end if
      if (do_spectral .and. size(reciprocal_obj%eigenvalues,1) /= nmat) then
         call g_logger%fatal('[exchange_q]: accepted eigensystem/Hamiltonian dimensions disagree', __FILE__, __LINE__)
      end if
      allocate(weights(nk))
      weights = reciprocal_obj%k_weights
      if (do_spectral) then
         allocate(evals(nmat,nk), evecs(nmat,nmat,nk))
         evals = reciprocal_obj%eigenvalues
         evecs = reciprocal_obj%eigenvectors
      end if
      if (size(weights) /= nk .or. sum(weights) <= tiny(1.0_rp)) then
         call g_logger%fatal('[exchange_q]: accepted k-mesh weights are invalid', __FILE__, __LINE__)
      end if

      allocate(finite_total(nsite,nsite,config%n_q), finite_tt(nsite,nsite,config%n_q), &
         finite_contact(nsite,nsite,config%n_q))
      finite_total = 0.0_rp; finite_tt = 0.0_rp; finite_contact = 0.0_rp
      allocate(spectral_total(nsite,nsite,config%n_q), spectral_tt(nsite,nsite,config%n_q), spectral_contact(nsite,nsite,config%n_q), &
         contour_total(nsite,nsite,config%n_q), contour_tt(nsite,nsite,config%n_q), contour_contact(nsite,nsite,config%n_q))
      spectral_total = 0.0_rp; spectral_tt = 0.0_rp; spectral_contact = 0.0_rp
      contour_total = 0.0_rp; contour_tt = 0.0_rp; contour_contact = 0.0_rp
      allocate(hessian(nsite,nsite), torque_torque(nsite,nsite), contact(nsite,nsite), complete(nsite,nsite))
      allocate(endpoint_index(nk), endpoint_reused(config%n_q), q_commensurate(config%n_q), &
         endpoint_residual(config%n_q), endpoint_mode(config%n_q))
      endpoint_reused = .false.; q_commensurate = .false.; endpoint_residual = huge(1.0_rp)
      endpoint_mode = 'explicit_diagonalization'
      endpoint_seconds = 0.0_rp; assembly_seconds = 0.0_rp; contraction_seconds = 0.0_rp
      hamiltonian_seconds = 0.0_rp; gf_seconds = 0.0_rp; solve_seconds = 0.0_rp; contour_seconds = 0.0_rp; total_response_seconds = 0.0_rp
      call system_clock(count_rate=clock_rate)
      finite_h_kT = max(reciprocal_obj%temperature*reciprocal_kb_ry_per_k, 1.0e-10_rp)
      vertex_identity_checked = .false.
      endpoint_identity_checked = .false.
      native_ready = config%native_crosscheck .or. config%native_turek
      allocate(native_jq_ud(config%n_q), native_jq_du(config%n_q), native_jq_sym(config%n_q), &
         native_delta_j(config%n_q), native_curvature(config%n_q))
      native_jq_ud = 0.0_rp; native_jq_du = 0.0_rp; native_jq_sym = 0.0_rp
      native_delta_j = 0.0_rp; native_curvature = 0.0_rp
      contour_options%contour_points = config%contour_points
      contour_options%contour_shape = config%contour_shape
      contour_options%contour_margin = config%contour_margin
      contour_options%contour_height_fraction = config%contour_height_fraction
      contour_options%account_fermi_poles = config%contour_account_fermi_poles
      native_contour_options%contour_points = config%native_contour_points
      native_contour_options%contour_shape = 'ellipse'
      native_contour_options%contour_margin = config%native_contour_margin
      native_contour_options%contour_height_fraction = config%native_contour_height_fraction
      native_contour_options%account_fermi_poles = config%native_contour_account_fermi_poles
      native_contour_options%target_fermi_poles = config%native_contour_target_fermi_poles

      call system_clock(response_start)

      do iq = 1, config%n_q
         call detect_commensurate_endpoint_map(reciprocal_obj%k_points, reciprocal_obj%k_weights, reciprocal_obj%nk_mesh, &
            q_direct(:,iq), endpoint_index, q_commensurate(iq), endpoint_residual(iq))
         call system_clock(clock_start)
         if (do_spectral .and. q_commensurate(iq)) then
            allocate(endpoint_evals(nmat,nk), endpoint_evecs(nmat,nmat,nk))
            do ik = 1, nk
               endpoint_evals(:,ik) = evals(:,endpoint_index(ik))
               endpoint_evecs(:,:,ik) = evecs(:,:,endpoint_index(ik))
            end do
            endpoint_reused(iq) = .true.
            endpoint_mode(iq) = 'mesh_reuse'
            if (.not. endpoint_identity_checked .and. sqrt(sum(q_direct(:,iq)**2)) > q_zero_tolerance) then
               ! The integer-cell map is a performance optimization, so verify
               ! its spectral identity once against the exact endpoint solve.
               call solve_unfolded_endpoints(reciprocal_obj, reciprocal_obj%k_points + &
                  spread_q(q_direct(:,iq),nk), endpoint_eval_check, endpoint_evec_check)
               endpoint_identity_error = maxval(abs(endpoint_evals-endpoint_eval_check))
               call g_logger%info('[exchange_q]: commensurate endpoint eigenvalue residual='// &
                  trim(real_to_string(endpoint_identity_error)), __FILE__, __LINE__)
               if (endpoint_identity_error > 2.0e-10_rp) then
                  call g_logger%fatal('[exchange_q]: commensurate endpoint spectral identity check failed', __FILE__, __LINE__)
               end if
               deallocate(endpoint_eval_check, endpoint_evec_check)
               endpoint_identity_checked = .true.
            end if
         else if (do_spectral) then
            call solve_unfolded_endpoints(reciprocal_obj, reciprocal_obj%k_points + spread_q(q_direct(:,iq),nk), &
                                          endpoint_evals, endpoint_evecs)
         end if
         call system_clock(clock_end)
         endpoint_seconds = endpoint_seconds + elapsed_clock_seconds(clock_start,clock_end,clock_rate)
         allocate(torques_q(nmat,nmat,nsite,nk), torques_minus_q(nmat,nmat,nsite,nk), &
                  mixed(nmat,nmat,nsite,nsite,nk))
         if (do_contour) allocate(h_source(nmat,nmat,nk), h_endpoint(nmat,nmat,nk))
         call system_clock(clock_start)
         if (do_contour) then
            do ik = 1, nk
               call assemble_lmto_hamiltonian(fixture, reciprocal_obj%k_points(:,ik), h_source(:,:,ik))
               call assemble_lmto_hamiltonian(fixture, reciprocal_obj%k_points(:,ik)+q_direct(:,iq), h_endpoint(:,:,ik))
            end do
         end if
         call system_clock(clock_end)
         hamiltonian_seconds = hamiltonian_seconds + elapsed_clock_seconds(clock_start,clock_end,clock_rate)
         call system_clock(clock_start)
         do ik = 1, nk
            call assemble_lmto_finite_q_torques(fixture, reciprocal_obj%k_points(:,ik), q_direct(:,iq), axes, &
               torques_q(:,:,:,ik))
            do ia = 1, nsite
               torques_minus_q(:,:,ia,ik) = transpose(conjg(torques_q(:,:,ia,ik)))
            end do
            do ia = 1, nsite
               do ja = 1, nsite
                  call assemble_lmto_finite_q_mixed_derivative(fixture, reciprocal_obj%k_points(:,ik), q_direct(:,iq), &
                     ia, axes(:,ia), ja, axes(:,ja), mixed(:,:,ia,ja,ik))
               end do
            end do
         end do
         if (.not. vertex_identity_checked .and. sqrt(sum(q_direct(:,iq)**2)) > q_zero_tolerance) then
            allocate(torque_minus_check(nmat,nmat,nsite)); vertex_error = 0.0_rp
            call assemble_lmto_finite_q_torques(fixture, reciprocal_obj%k_points(:,1)+q_direct(:,iq), &
               -q_direct(:,iq), axes, torque_minus_check)
            do ia = 1, nsite
               vertex_error = max(vertex_error, maxval(abs(torque_minus_check(:,:,ia)-torques_minus_q(:,:,ia,1))))
            end do
            if (vertex_error > 5.0e-10_rp) call g_logger%fatal('[exchange_q]: finite-q vertex adjoint/gauge check failed', &
               __FILE__, __LINE__)
            vertex_identity_checked = .true.
            deallocate(torque_minus_check)
         end if
         call system_clock(clock_end)
         assembly_seconds = assembly_seconds + elapsed_clock_seconds(clock_start,clock_end,clock_rate)
         if (do_spectral) then
            call system_clock(clock_start)
            if (config%finite_h_spectral_mode == 'legacy_occupied') then
               call force_theorem_finite_q_hessian_from_eigenbasis_batch(evals, evecs, endpoint_evals, endpoint_evecs, &
                  reciprocal_obj%fermi_level, weights, torques_q, torques_minus_q, mixed, hessian, torque_torque, contact, complete)
            else
               call force_theorem_finite_q_hessian_from_eigenbasis_metallic_batch(evals, evecs, endpoint_evals, endpoint_evecs, &
                  reciprocal_obj%fermi_level, finite_h_kT, weights, torques_q, torques_minus_q, mixed, hessian, torque_torque, contact, complete)
            end if
            call system_clock(clock_end)
            contraction_seconds = contraction_seconds + elapsed_clock_seconds(clock_start,clock_end,clock_rate)
            spectral_tt(:,:,iq) = real(torque_torque,rp)
            spectral_contact(:,:,iq) = real(contact,rp)
            spectral_total(:,:,iq) = real(complete,rp)
         end if
         if (do_contour) then
            call system_clock(clock_start)
            call force_theorem_finite_q_hessian_from_resolvent_batch(h_source, h_endpoint, reciprocal_obj%fermi_level, finite_h_kT, &
               weights, torques_q, torques_minus_q, mixed, contour_options, hessian, torque_torque, contact, complete, contour_report)
            call system_clock(clock_end)
            gf_seconds = gf_seconds + elapsed_clock_seconds(clock_start,clock_end,clock_rate)
            solve_seconds = solve_seconds + contour_report%solve_seconds
            contour_seconds = contour_seconds + contour_report%contour_seconds
            contour_tt(:,:,iq) = real(torque_torque,rp)
            contour_contact(:,:,iq) = real(contact,rp)
            contour_total(:,:,iq) = real(complete,rp)
         end if
         if (config%finite_h_response_backend == 'contour') then
            finite_tt(:,:,iq) = contour_tt(:,:,iq); finite_contact(:,:,iq) = contour_contact(:,:,iq); finite_total(:,:,iq) = contour_total(:,:,iq)
         else
            finite_tt(:,:,iq) = spectral_tt(:,:,iq); finite_contact(:,:,iq) = spectral_contact(:,:,iq); finite_total(:,:,iq) = spectral_total(:,:,iq)
         end if
         if (allocated(h_source)) deallocate(h_source, h_endpoint)
         deallocate(torques_q, torques_minus_q, mixed)
         if (allocated(endpoint_evals)) deallocate(endpoint_evals, endpoint_evecs)
      end do
      call system_clock(clock_end)
      total_response_seconds = elapsed_clock_seconds(response_start,clock_end,clock_rate)

      if (native_ready) then
         call g_logger%info('[exchange_q]: native reference uses accepted reciprocal SCF k mesh, weights, EF, temperature, lattice and potential; order='// &
            trim(reciprocal_obj%kspace_ham_order),__FILE__,__LINE__)
         call native_exchange_q_driver(lattice_obj, reciprocal_obj, q_direct, native_contour_options, native_jq_ud, native_jq_du, &
            native_jq_sym, native_delta_j, native_curvature, native_contour_report)
         call g_logger%info('[exchange_q]: native Turek contour points/poles='//trim(int2str(native_contour_report%contour_points))//'/'// &
            trim(int2str(native_contour_report%fermi_poles))//' solve/contour/pole='// &
            trim(real_to_string(native_contour_report%solve_seconds))//'/'//trim(real_to_string(native_contour_report%contour_seconds))//'/'// &
            trim(real_to_string(native_contour_report%pole_seconds))//' s; native max ellipse='// &
            trim(real_to_string(native_contour_report%native_max_ellipse_value)), __FILE__, __LINE__)
      end if
      call g_logger%info('[exchange_q]: timing endpoint/hamiltonian+T/C/spectral='//trim(real_to_string(endpoint_seconds))//'/'// &
         trim(real_to_string(hamiltonian_seconds))//'/'//trim(real_to_string(assembly_seconds))//'/'//trim(real_to_string(contraction_seconds))//' s', __FILE__, __LINE__)
      if (config%rotation_dynamics) call run_native_rotation_dynamics_campaign(config,lattice_obj,reciprocal_obj,self_obj,fixture, &
         q_direct,q_cart,finite_total,native_delta_j,native_curvature,native_ready)
      call write_exchange_q_output(config, lattice_obj, reciprocal_obj, q_direct, q_cart, finite_tt, finite_contact, &
         finite_total, native_jq_ud, native_jq_du, native_jq_sym, native_delta_j, native_curvature, native_contour_report, native_ready, endpoint_mode, endpoint_reused, q_commensurate, endpoint_residual, &
         endpoint_seconds, assembly_seconds, contraction_seconds, hamiltonian_seconds, spectral_tt, spectral_contact, spectral_total, contour_tt, &
         contour_contact, contour_total, gf_seconds, solve_seconds, contour_seconds, total_response_seconds)
      if (native_ready) then
         call compute_native_gamma_and_report(native_delta_j, q_direct, native_gamma)
         native_q = maxval(abs(native_delta_j))
         call g_logger%info('[exchange_q]: native Turek J(q) completed; |DeltaJ|_max='//trim(real_to_string(native_q))// &
            ' Ry, gamma='//trim(real_to_string(native_gamma))//' Ry', __FILE__, __LINE__)
      end if
      call fixture%clear()
      if (allocated(evals)) deallocate(evals)
      if (allocated(evecs)) deallocate(evecs)
      deallocate(q_direct, q_cart, weights, finite_total, finite_tt, finite_contact, spectral_total, spectral_tt, spectral_contact, &
         contour_total, contour_tt, contour_contact, &
         axes, hessian, torque_torque, contact, complete, endpoint_index, endpoint_reused, q_commensurate, &
         endpoint_residual, endpoint_mode)
      deallocate(native_jq_ud, native_jq_du, native_jq_sym, native_delta_j, native_curvature)
   end subroutine run_exchange_q

   subroutine run_native_rotation_dynamics_campaign(config,lat,recip,self_obj,fixture,q_direct,q_cart,finite_h,native_delta,native_curv,native_ready)
      type(exchange_q_config), intent(in) :: config
      type(lattice), intent(in) :: lat
      type(reciprocal), target, intent(inout) :: recip
      type(self), intent(in) :: self_obj
      type(lmto_live_hamiltonian_fixture), target, intent(in) :: fixture
      real(rp), intent(in) :: q_direct(:, :),q_cart(:, :),finite_h(:, :, :),native_delta(:),native_curv(:)
      logical, intent(in) :: native_ready
      type(rotation_state), target :: state,minus_state
      type(rotation_request) :: request
      type(rotation_result) :: response,minus_response,plus_probe,minus_probe
      complex(rp), allocatable :: reduced(:, :)
      real(rp), parameter :: ry_to_mev=13605.693122994_rp
      real(rp), parameter :: eta_ladder(3)=[1.0e-4_rp,2.5e-5_rp,6.25e-6_rp]
      real(rp) :: bplus,bminus,slope_step,slope_eta,slope_rel,qmag,scf_moment,field_max,cov_q(3)
      real(rp) :: static_residual,q0_residual,circular_offdiag,covariance_residual,omega_cov,eta_cov
      real(rp) :: finite_mev,turek_mev,omega_max,pred_plus,pred_minus,expected,df,fit_d,fit_resid
      real(rp) :: pole_re,pole_min,loss_peak,loss_height,fwhm,fwhm_mev,resolution,re_k,im_k,abs_k,pole_slope
      real(rp) :: energy_eta(3),resid_eta(3),eta_move,window_min,window_max
      integer :: unit,iq,igamma,ik,site,channel,ieta,ipole,nfinite,clean_count
      logical :: static_ok,berry_ok,circular_ok,covariance_ok,causal_ok,pole_ok(3),resolved,fit_ok
      character(len=8) :: channel_name
      character(len=24) :: gate_status
      character(len=10) :: pole_status

      if (.not.native_ready) error stop 'rotation dynamics requires native_turek=true for the independent adiabatic reference'
      if (fixture%nsite/=1) error stop 'first Fe rotation pole driver currently reports the one-site primitive bcc state'
      if (recip%hamiltonian%ccor_2c) error stop 'first Fe rotation pole state requires CCOR off'
      if (abs(recip%temperature-300.0_rp)>1.0e-8_rp) error stop 'first Fe rotation pole state requires T=300 K'
      if (size(finite_h,3)/=size(q_direct,2) .or. size(native_delta)/=size(q_direct,2) .or. &
          size(native_curv)/=size(q_direct,2)) error stop 'rotation dynamics campaign q/reference shape mismatch'
      if (size(q_cart,2)/=size(q_direct,2)) error stop 'rotation dynamics campaign q Cartesian shape mismatch'
      if (len_trim(config%rotation_output_file)==0) error stop 'rotation dynamics output path is blank'
      open(newunit=unit,file=trim(config%rotation_output_file),status='replace',action='write')
      write(unit,'(a)') '# native second-order local-rotation dynamics; K^R = contact + unsymmetrized retarded bubble'
      write(unit,'(a)') '# energies are Ry internally; reported q_Ainv uses 2*pi/alat times the actual Cartesian reciprocal vector'
      write(unit,'(a)') '# q_fraction q_Ainv adiabatic_finiteH_meV adiabatic_Turek_meV pole_ReK_meV pole_minK_meV loss_peak_meV eta_Ry FWHM_meV max_minus_ImG_invRy circular_channel ReK_Ry ImK_Ry absK_Ry frequency_resolution_Ry status'
      write(*,'(a,3(i0,1x))') 'Rotation dynamics accepted k mesh = ',recip%nk_mesh
      write(*,'(a,es16.8)') 'Rotation dynamics EF (Ry) = ',recip%fermi_level
      write(*,'(a,es16.8)') 'Rotation dynamics temperature (K) = ',recip%temperature
      write(*,'(a,es16.8)') 'Rotation dynamics SCF constraining field max (Ry) = ',max_constraint_field(self_obj,lat%nrec,lat%nbulk)

      igamma=0
      do iq=1,size(q_direct,2)
         if (maxval(abs(q_direct(:,iq)))<q_zero_tolerance) igamma=iq
      end do
      if (igamma==0) error stop 'rotation dynamics q list has no Gamma point'
      allocate(reduced(2,2))
      call prepare_rotation_response(fixture,recip,q_direct(:,igamma),state)
      request%state=>state; request%omega=0.0_rp; request%eta=0.0_rp; request%exact_static=.true.; request%want_inverse=.false.
      call evaluate_rotation_response(request,response)
      call reduce_static_rotation_kernel(response%kernel,reduced)
      q0_residual=maxval(abs(reduced))
      request%exact_static=.false.; request%omega=1.0e-5_rp; request%eta=1.0e-9_rp
      call evaluate_rotation_response(request,plus_probe)
      request%omega=-1.0e-5_rp
      call evaluate_rotation_response(request,minus_probe)
      slope_step=1.0e-5_rp; slope_eta=1.0e-9_rp
      bplus=real((plus_probe%kernel_pm(1,1)-minus_probe%kernel_pm(1,1))/(2.0_rp*slope_step),rp)
      bminus=real((plus_probe%kernel_pm(2,2)-minus_probe%kernel_pm(2,2))/(2.0_rp*slope_step),rp)
      slope_rel=max(abs(abs(bplus)-abs(response%berry)),abs(abs(bminus)-abs(response%berry)),abs(bplus+bminus))/ &
         max(abs(response%berry),1.0e-12_rp)
      circular_offdiag=max(abs(reduced(1,1)-reduced(2,2)),abs(reduced(1,2)),abs(reduced(2,1)))
      scf_moment=0.0_rp
      do site=1,lat%nrec
         scf_moment=scf_moment+self_obj%symbolic_atom(lat%nbulk+site)%potential%mtot
      end do
      write(*,'(a,es16.8)') 'Rotation dynamics electron-number residual = ',state%electron_residual
      write(*,'(a,es16.8)') 'Rotation dynamics M_band = ',state%magnetization
      write(*,'(a,es16.8)') 'Rotation dynamics SCF magnetic moment = ',scf_moment
      write(*,'(a,es16.8)') 'Rotation dynamics Berry commutator = ',response%berry
      write(*,'(a,es12.4)') 'Rotation dynamics direct Berry vs M/2 residual = ',response%berry_residual
      write(*,'(a,2(es16.8,1x),a,es12.4)') 'Rotation dynamics circular slopes (+,-) = ',bplus,bminus, &
         ' relative residual=',slope_rel
      write(*,'(a,es12.4)') 'Rotation dynamics q0 Goldstone residual = ',q0_residual
      write(*,'(a,es12.4)') 'Rotation dynamics circular off-diagonal residual = ',max(abs(reduced(1,2)),abs(reduced(2,1)))
      berry_ok=slope_rel<=5.0e-2_rp .and. abs(abs(bplus)-abs(bminus))<=5.0e-2_rp*max(abs(response%berry),1.0e-12_rp)
      circular_ok=max(abs(reduced(1,2)),abs(reduced(2,1)))<=1.0e-7_rp
      q0_residual=maxval(abs(reduced))
      static_ok=q0_residual<=2.0e-7_rp
      if (abs(response%berry)<=1.0e-8_rp) berry_ok=.false.

      write(unit,'(a)') '# state_provenance = accepted second-order ham_only; SOC off; CCOR off; constraining field zero'
      write(unit,'(a,3(i0,1x))') '# k_mesh = ',recip%nk_mesh
      write(unit,'(a,es20.12)') '# EF_Ry = ',recip%fermi_level
      write(unit,'(a,es20.12)') '# electron_number_residual = ',state%electron_residual
      write(unit,'(a,es20.12)') '# M_band = ',state%magnetization
      write(unit,'(a,es20.12)') '# SCF_magnetic_moment = ',scf_moment
      write(unit,'(a,es20.12)') '# Berry_commutator = ',response%berry
      write(unit,'(a,es20.12)') '# direct_Berry_vs_Mband_over_2_residual = ',response%berry_residual
      write(unit,'(a,es20.12)') '# q0_slope_plus = ',bplus
      write(unit,'(a,es20.12)') '# q0_slope_minus = ',bminus
      write(unit,'(a,es20.12)') '# q0_Goldstone_residual = ',q0_residual
      write(unit,'(a,es20.12)') '# circular_offdiagonal_residual = ',max(abs(reduced(1,2)),abs(reduced(2,1)))

      pole_ok=.false.; energy_eta=-1.0_rp; resid_eta=huge(1.0_rp); nfinite=0; clean_count=0
      covariance_ok=.false.; covariance_residual=huge(1.0_rp); causal_ok=.false.
      do iq=1,size(q_direct,2)
         if (iq==igamma) cycle
         nfinite=nfinite+1
         call prepare_rotation_response(fixture,recip,q_direct(:,iq),state)
         request%state=>state; request%omega=0.0_rp; request%eta=0.0_rp; request%exact_static=.true.; request%want_inverse=.false.
         call evaluate_rotation_response(request,response)
         finite_mev=2.0_rp*real(finite_h(1,1,iq),rp)/max(abs(state%magnetization),1.0e-12_rp)*ry_to_mev
         turek_mev=4.0_rp*native_delta(iq)/max(abs(state%magnetization),1.0e-12_rp)*ry_to_mev
         qmag=2.0_rp*pi/lat%alat*sqrt(sum(q_cart(:,iq)**2))
         if (iq==2) then
            call prepare_rotation_response(fixture,recip,-q_direct(:,iq),minus_state)
            request%state=>minus_state; request%exact_static=.true.; request%omega=0.0_rp; request%eta=0.0_rp
            call evaluate_rotation_response(request,minus_response)
            call reduce_static_rotation_kernel(response%kernel,reduced)
            static_residual=maxval(abs(reduced(1:1,1:1)-finite_h(1:1,1:1,iq)))
            write(*,'(a,es12.4)') 'Rotation dynamics q,-q static reduction residual = ',static_residual
            static_ok=static_ok .and. static_residual<=max(1.0e-8_rp,1.0e-7_rp*maxval(abs(finite_h(:,:,iq))))
            ! The covariance identity is evaluated at a non-self-inverse q
            ! that translates the sampled mesh exactly; arbitrary q remains
            ! supported by the response API and the separate direct oracle.
            cov_q=[1.0_rp/real(recip%nk_mesh(1),rp),0.0_rp,0.0_rp]
            call state%clear(); call minus_state%clear()
            call prepare_rotation_response(fixture,recip,cov_q,state)
            call prepare_rotation_response(fixture,recip,-cov_q,minus_state)
            request%state=>state; request%exact_static=.false.; request%omega=2.0e-4_rp; request%eta=5.0e-5_rp
            call evaluate_rotation_response(request,plus_probe)
            request%state=>minus_state; request%omega=-2.0e-4_rp
            call evaluate_rotation_response(request,minus_probe)
            ! The live Bloch endpoint convention satisfies K_AB(q,w) =
            ! conjg(K_AB(-q,-w)) for real Cartesian coordinates.
            covariance_residual=maxval(abs(plus_probe%kernel-conjg(minus_probe%kernel)))/ &
               max(1.0_rp,maxval(abs(plus_probe%kernel)))
            write(*,'(a,3(es12.4,1x))') 'Rotation dynamics q/-q/-omega covariance q/residual = ',cov_q,covariance_residual
            covariance_ok=covariance_residual<=1.0e-7_rp
            call state%clear(); call minus_state%clear()
            call prepare_rotation_response(fixture,recip,q_direct(:,iq),state)
         else
            static_residual=abs(real(response%kernel(1,1),rp)-finite_h(1,1,iq))
            if (nfinite==1) static_ok=static_ok .and. static_residual<=max(1.0e-8_rp,1.0e-7_rp*abs(finite_h(1,1,iq)))
         end if
         if (iq==2) then
            pred_plus=-real(response%kernel_pm(1,1),rp)/bplus
            pred_minus=-real(response%kernel_pm(2,2),rp)/bminus
         else
            pred_plus=-real(response%kernel_pm(1,1),rp)/bplus
            pred_minus=-real(response%kernel_pm(2,2),rp)/bminus
         end if
         channel=0
         if (pred_plus>0.0_rp .and. pred_minus<=0.0_rp) channel=1
         if (pred_minus>0.0_rp .and. pred_plus<=0.0_rp) channel=2
         if (channel==0 .and. pred_plus>0.0_rp .and. pred_minus>0.0_rp) then
            if (abs(pred_plus-finite_mev/ry_to_mev)<=abs(pred_minus-finite_mev/ry_to_mev)) then
               channel=1
            else
               channel=2
            end if
         end if
         if (channel==1) then
            channel_name='+ '
            expected=pred_plus
         else if (channel==2) then
            channel_name='- '
            expected=pred_minus
         else
            channel_name='none'
            expected=0.0_rp
         end if
         omega_max=max(1.0e-3_rp,2.5_rp*max(abs(finite_mev),abs(turek_mev))/ry_to_mev)
         omega_max=min(omega_max,2.0e-2_rp)
         write(*,'(a,i0,a,3(es12.4,1x),a,es12.4,a,es12.4)') 'Rotation q index ',iq,' finite-H/Turek/meV=', &
            finite_mev,turek_mev,qmag,' A^-1 omega_window_Ry=',omega_max,' predicted_Ry=',expected
         if (channel>0) then
            if (channel==1) then
               ipole=1
            else
               ipole=2
            end if
            do ieta=1,size(eta_ladder)
               call scan_rotation_pole(state,ipole,omega_max,eta_ladder(ieta),pole_re,pole_min,loss_peak,loss_height, &
                  fwhm,resolution,re_k,im_k,abs_k,pole_slope,resolved)
               pole_ok(nfinite)=resolved .or. pole_ok(nfinite)
               causal_ok=causal_ok .or. (pole_re>0.0_rp .and. pole_min>0.0_rp .and. pole_slope*im_k>0.0_rp)
               if (ieta==size(eta_ladder)) then
                  energy_eta(nfinite)=pole_re*ry_to_mev
                  resid_eta(nfinite)=resolution*ry_to_mev
               end if
               if (resolved) then
                  pole_status='RESOLVED'
               else
                  pole_status='UNRESOLVED'
               end if
               fwhm_mev=-1.0_rp
               if (fwhm>=0.0_rp) fwhm_mev=fwhm*ry_to_mev
               write(unit,'(3(es14.6,1x),9(es14.6,1x),a,4(es14.6,1x),a)') q_direct(:,iq),qmag,finite_mev,turek_mev, &
                  pole_re*ry_to_mev,pole_min*ry_to_mev,loss_peak*ry_to_mev,eta_ladder(ieta),fwhm_mev,loss_height, &
                  trim(channel_name),re_k,im_k,abs_k,resolution,trim(pole_status)
               write(*,'(a,3(es14.6,1x),a,es12.4,a,a,a,l1)') '  pole ReK/minK/loss (meV) = ', &
                  pole_re*ry_to_mev,pole_min*ry_to_mev,loss_peak*ry_to_mev,' eta=',eta_ladder(ieta), &
                  ' channel=',trim(channel_name),' resolved=',resolved
            end do
            if (pole_ok(nfinite)) clean_count=clean_count+1
         else
            write(unit,'(3(es14.6,1x),9(es14.6,1x),a,4(es14.6,1x),a)') q_direct(:,iq),qmag,finite_mev,turek_mev, &
               -1.0_rp,-1.0_rp,-1.0_rp,eta_ladder(3),-1.0_rp,0.0_rp,'none',0.0_rp,0.0_rp,0.0_rp,0.0_rp,'NO_POSITIVE_CHANNEL'
         end if
         call state%clear()
         call minus_state%clear()
      end do

      ! covariance_ok is set only by the explicit q/-q/-omega comparison above.
      berry_ok=berry_ok .and. slope_rel<=5.0e-2_rp
      static_ok=static_ok .and. q0_residual<=2.0e-7_rp
      fit_ok=.false.; fit_d=0.0_rp; fit_resid=-1.0_rp; window_min=0.0_rp; window_max=0.0_rp
      if (nfinite>=3 .and. clean_count>=3) then
         fit_ok=.true.; expected=0.0_rp; window_min=huge(1.0_rp); window_max=0.0_rp
         do ipole=1,3
            if (energy_eta(ipole)<=0.0_rp) fit_ok=.false.
            if (ipole>1 .and. energy_eta(ipole)<=energy_eta(ipole-1)) fit_ok=.false.
            expected=expected+energy_eta(ipole)*(2.0_rp*pi/lat%alat*sqrt(sum(q_cart(:,ipole+1)**2)))**2
            fit_d=fit_d+(2.0_rp*pi/lat%alat*sqrt(sum(q_cart(:,ipole+1)**2)))**4
            window_min=min(window_min,2.0_rp*pi/lat%alat*sqrt(sum(q_cart(:,ipole+1)**2)))
            window_max=max(window_max,2.0_rp*pi/lat%alat*sqrt(sum(q_cart(:,ipole+1)**2)))
         end do
         if (fit_ok .and. fit_d>tiny(1.0_rp)) then
            fit_d=expected/fit_d
            fit_resid=0.0_rp
            do ipole=1,3
               qmag=2.0_rp*pi/lat%alat*sqrt(sum(q_cart(:,ipole+1)**2))
               fit_resid=fit_resid+(energy_eta(ipole)-fit_d*qmag*qmag)**2
            end do
            fit_resid=sqrt(fit_resid/3.0_rp)
         else
            fit_ok=.false.
         end if
      end if
      if (fit_ok .and. clean_count>=3 .and. static_ok .and. berry_ok .and. circular_ok .and. covariance_ok .and. causal_ok) then
         gate_status='PASS-A'
      else if (static_ok .and. berry_ok .and. circular_ok .and. covariance_ok .and. causal_ok) then
         gate_status='PASS-B'
      else
         gate_status='BLOCKED'
      end if
      write(*,'(a,l1)') 'Rotation dynamics causal positive-frequency pole = ',causal_ok
      write(*,'(a,a)') 'Rotation dynamics implementation = ',trim(gate_status)
      if (fit_ok) then
         write(*,'(a)') 'Fe material convergence = PRELIMINARY — MATERIAL CONVERGENCE OPEN'
      else
         write(*,'(a)') 'Fe material convergence = OPEN'
      end if
      if (fit_ok) write(*,'(a,es16.8,a,es16.8,a,2(es12.4,1x))') 'Preliminary D (meV A2) = ',fit_d, &
         ' fit residual (meV) = ',fit_resid,' q window (A^-1) = ',window_min,window_max
      write(*,'(a)') 'Goldstone correction = OFF; strict-ASA L0 dynamics = NOT RUN'
      close(unit)
      call state%clear(); call minus_state%clear()
      deallocate(reduced)
   end subroutine run_native_rotation_dynamics_campaign

   subroutine scan_rotation_pole(state,channel,omega_max,eta,pole_re,pole_min,loss_peak,loss_height,fwhm,resolution, &
      re_k,im_k,abs_k,pole_slope,resolved)
      type(rotation_state), target, intent(inout) :: state
      integer, intent(in) :: channel
      real(rp), intent(in) :: omega_max,eta
      real(rp), intent(out) :: pole_re,pole_min,loss_peak,loss_height,fwhm,resolution,re_k,im_k,abs_k,pole_slope
      logical, intent(out) :: resolved
      type(rotation_request) :: request
      type(rotation_result) :: response
      integer, parameter :: ncoarse=61,nfine=41
      real(rp) :: omega_c(ncoarse),kabs_c(ncoarse),kre_c(ncoarse),loss_c(ncoarse)
      real(rp) :: omega_f(nfine),kabs_f(nfine),kre_f(nfine),loss_f(nfine)
      real(rp) :: center,dw,left,right,half,peak,root_distance,x1,x2,t
      integer :: i,imin,ipeak,iroot
      logical :: have_root,have_left,have_right
      request%state=>state; request%eta=eta; request%exact_static=.false.; request%want_inverse=.true.
      do i=1,ncoarse
         omega_c(i)=omega_max*real(i-1,rp)/real(ncoarse-1,rp)
         request%omega=omega_c(i)
         call evaluate_rotation_response(request,response)
         kabs_c(i)=abs(response%kernel_pm(channel,channel))
         kre_c(i)=real(response%kernel_pm(channel,channel),rp)
         loss_c(i)=-aimag(response%inverse_kernel_pm(channel,channel))
      end do
      imin=minloc(kabs_c,dim=1); center=omega_c(imin); dw=omega_max/real(ncoarse-1,rp)
      left=max(0.0_rp,center-2.0_rp*dw); right=min(omega_max,center+2.0_rp*dw)
      if (right<=left) right=min(omega_max,left+dw)
      resolution=(right-left)/real(nfine-1,rp)
      do i=1,nfine
         omega_f(i)=left+resolution*real(i-1,rp)
         request%omega=omega_f(i)
         call evaluate_rotation_response(request,response)
         kabs_f(i)=abs(response%kernel_pm(channel,channel))
         kre_f(i)=real(response%kernel_pm(channel,channel),rp)
         loss_f(i)=-aimag(response%inverse_kernel_pm(channel,channel))
      end do
      imin=minloc(kabs_f,dim=1); pole_min=omega_f(imin); pole_re=-1.0_rp; have_root=.false.; iroot=0
      root_distance=huge(1.0_rp)
      do i=1,nfine-1
         if (kre_f(i)==0.0_rp) then
            t=0.0_rp
         else if (kre_f(i)*kre_f(i+1)<0.0_rp) then
            t=abs(kre_f(i))/(abs(kre_f(i))+abs(kre_f(i+1)))
         else
            cycle
         end if
         x1=omega_f(i)+t*(omega_f(i+1)-omega_f(i))
         if (abs(x1-center)<root_distance) then
            root_distance=abs(x1-center); pole_re=x1; iroot=i; have_root=.true.
         end if
      end do
      ipeak=maxloc(loss_f,dim=1); loss_peak=omega_f(ipeak); peak=loss_f(ipeak); loss_height=peak; fwhm=-1.0_rp
      half=0.5_rp*peak; have_left=.false.; have_right=.false.; x1=left; x2=right
      if (peak>0.0_rp) then
         do i=ipeak-1,1,-1
            if (loss_f(i)<=half .and. loss_f(i+1)>half) then
               t=(half-loss_f(i))/(loss_f(i+1)-loss_f(i))
               x1=omega_f(i)+t*(omega_f(i+1)-omega_f(i)); have_left=.true.; exit
            end if
         end do
         do i=ipeak,nfine-1
            if (loss_f(i)>half .and. loss_f(i+1)<=half) then
               t=(half-loss_f(i))/(loss_f(i+1)-loss_f(i))
               x2=omega_f(i)+t*(omega_f(i+1)-omega_f(i)); have_right=.true.; exit
            end if
         end do
         if (have_left .and. have_right) fwhm=x2-x1
      end if
      pole_slope=0.0_rp
      if (iroot>0) pole_slope=(kre_f(iroot+1)-kre_f(iroot))/(omega_f(iroot+1)-omega_f(iroot))
      ! Report the kernel components at the refined Re K=0 crossing.  The
      ! separate pole_min diagnostic remains the minimum of |K|.
      request%omega=merge(pole_re,pole_min,pole_re>0.0_rp)
      call evaluate_rotation_response(request,response)
      re_k=real(response%kernel_pm(channel,channel),rp)
      im_k=aimag(response%kernel_pm(channel,channel))
      abs_k=abs(response%kernel_pm(channel,channel))
      resolved=have_root .and. peak>0.0_rp .and. abs(pole_re-pole_min)<=max(2.0_rp*eta,2.0_rp*resolution)
      if (peak>0.0_rp) resolved=resolved .and. abs(loss_peak-pole_min)<=max(2.0_rp*eta,2.0_rp*resolution)
      if (iroot==0) resolved=.false.
   end subroutine scan_rotation_pole

   real(rp) function max_constraint_field(self_obj,nsite,nbulk) result(value)
      type(self), intent(in) :: self_obj
      integer, intent(in) :: nsite,nbulk
      integer :: site
      value=0.0_rp
      do site=1,nsite
         value=max(value,maxval(abs(self_obj%symbolic_atom(nbulk+site)%mag_cfield)))
      end do
   end function max_constraint_field

   subroutine validate_exchange_q_capability(config, control_obj, ham, self_obj, recip)
      type(exchange_q_config), intent(in) :: config
      type(control), intent(in) :: control_obj
      type(hamiltonian), intent(in) :: ham
      type(self), intent(in) :: self_obj
      type(reciprocal), intent(in) :: recip
      logical :: native_requested

      native_requested=config%native_crosscheck .or. config%native_turek
      if (config%n_q<2) call g_logger%fatal('[exchange_q]: q path is empty or has no finite-q point',__FILE__,__LINE__)
      if (trim(control_obj%calctype)/='B' .or. control_obj%nsp/=1 .or. control_obj%has_soc()) then
         call g_logger%fatal('[exchange_q]: capability gate requires bulk nsp=1 scalar-relativistic collinear state with SOC off', &
            __FILE__,__LINE__)
      end if
      if (.not.self_obj%use_kspace .or. .not.allocated(self_obj%reciprocal_scf_cache)) then
         call g_logger%fatal('[exchange_q]: requires self%use_kspace=.true. to consume the accepted reciprocal SCF state', &
            __FILE__,__LINE__)
      end if
      if (ham%ccor_2c .or. ham%hubbard_u_general_check .or. ham%hubbard_v_check .or. control_obj%constraints_enable) then
         call g_logger%fatal('[exchange_q]: unsupported CCOR/Hubbard/constraint combination',__FILE__,__LINE__)
      end if
      if (trim(recip%reciprocal_mode)/='ham_only') then
         call g_logger%fatal('[exchange_q]: capability gate requires reciprocal_mode=ham_only',__FILE__,__LINE__)
      end if
      if (trim(recip%kspace_ham_order)/='first' .and. trim(recip%kspace_ham_order)/='second') then
         call g_logger%fatal('[exchange_q]: reciprocal Hamiltonian order must be first or second',__FILE__,__LINE__)
      end if
      if (fixture_basis_size(ham)/=9) then
         call g_logger%fatal('[exchange_q]: capability gate requires the full spd production basis',__FILE__,__LINE__)
      end if
      if (native_requested) call validate_native_turek_capability(config,ham)
      call validate_finite_h_capability(ham,recip)
   end subroutine validate_exchange_q_capability

   ! Native P-S contour capability uses only lattice/potential and accepted
   ! reciprocal integration metadata; its gate is independent of HOH assembly.
   subroutine validate_native_turek_capability(config,ham)
      type(exchange_q_config), intent(in) :: config
      type(hamiltonian), intent(in) :: ham
      if ((config%native_crosscheck .or. config%native_turek) .and. ham%charge%lattice%nrec/=1) then
         call g_logger%fatal('[exchange_q]: native contour production is currently capability-gated to one sublattice', &
            __FILE__,__LINE__)
      end if
   end subroutine validate_native_turek_capability

   ! Finite-H torque/contact vertices require a live orthogonal H(k) adapter.
   ! The second-order branch is accepted only when reciprocal assembly inputs
   ! exist; the reciprocal assembler otherwise silently falls back to first order.
   subroutine validate_finite_h_capability(ham,recip)
      type(hamiltonian), intent(in) :: ham
      type(reciprocal), intent(in) :: recip
      if (trim(recip%kspace_ham_order)=='second') then
         if (.not.ham%hoh .or. .not.allocated(ham%eeo) .or. .not.allocated(ham%enim)) then
            call g_logger%fatal('[exchange_q]: SECOND_ORDER_STATE_NOT_ACTIVE in reciprocal finite-H path',__FILE__,__LINE__)
         end if
      end if
   end subroutine validate_finite_h_capability

   integer function fixture_basis_size(ham) result(size_orb)
      type(hamiltonian), intent(in) :: ham
      if (.not. associated(ham%charge)) then
         size_orb = 0
      else
         size_orb = size(ham%charge%lattice%sbar,1)
      end if
   end function fixture_basis_size

   function spread_q(q, count) result(points)
      real(rp), intent(in) :: q(3)
      integer, intent(in) :: count
      real(rp) :: points(3,count)
      integer :: i
      do i = 1, count
         points(:,i) = q
      end do
   end function spread_q

   subroutine solve_unfolded_endpoints(recip, points, eigenvalues, eigenvectors)
      type(reciprocal), intent(inout) :: recip
      real(rp), intent(in) :: points(:, :)
      real(rp), allocatable, intent(out) :: eigenvalues(:, :)
      complex(rp), allocatable, intent(out) :: eigenvectors(:, :, :)
      complex(rp), allocatable :: h(:, :), work(:)
      real(rp), allocatable :: rwork(:)
      integer :: nmat, nk, ik, lwork, info
      external :: zheev

      nk = size(points,2)
      nmat = size(recip%eigenvalues,1)
      allocate(eigenvalues(nmat,nk), eigenvectors(nmat,nmat,nk), h(nmat,nmat), rwork(max(1,3*nmat-2)))
      lwork = max(1, 2*nmat-1)
      allocate(work(lwork))
      do ik = 1, nk
         call recip%build_hamiltonian_at_kpoint(points(:,ik), h)
         call zheev('V','U',nmat,h,nmat,eigenvalues(:,ik),work,lwork,rwork,info)
         if (info /= 0) call g_logger%fatal('[exchange_q]: endpoint diagonalization failed, info='//int2str(info), &
            __FILE__, __LINE__)
         eigenvectors(:,:,ik) = h
      end do
      deallocate(h,work,rwork)
   end subroutine solve_unfolded_endpoints

   subroutine native_exchange_q_driver(lat, recip, q_direct, options, jq_ud_out, jq_du_out, jq_sym_out, delta_j, curvature, report)
      type(lattice), intent(inout) :: lat
      type(reciprocal), intent(in) :: recip
      real(rp), intent(in) :: q_direct(:, :)
      type(native_turek_contour_options), intent(in) :: options
      real(rp), intent(out) :: jq_ud_out(:), jq_du_out(:), jq_sym_out(:), delta_j(:), curvature(:)
      type(native_turek_contour_report), intent(out) :: report
      complex(rp), allocatable :: jq_ud(:, :, :), jq_du(:, :, :), jq_sym(:, :, :), native_delta(:, :, :), native_curvature(:, :, :)
      real(rp) :: kT
      integer :: nsite, nq

      nsite=lat%nrec; nq=size(q_direct,2)
      if (nsite/=1) error stop 'native_exchange_q_driver: production scalar reduction requires nsite=1'
      if (size(q_direct,1)/=3 .or. size(jq_ud_out)/=nq .or. size(jq_du_out)/=nq .or. &
          size(jq_sym_out)/=nq .or. size(delta_j)/=nq .or. size(curvature)/=nq) then
         error stop 'native_exchange_q_driver: shape mismatch'
      end if
      kT=max(recip%temperature*reciprocal_kb_ry_per_k,1.0e-10_rp)
      allocate(jq_ud(nsite,nsite,nq),jq_du(nsite,nsite,nq),jq_sym(nsite,nsite,nq), &
         native_delta(nsite,nsite,nq),native_curvature(nsite,nsite,nq))
      call native_turek_static_reference(lat,recip%k_points,recip%k_weights,q_direct,recip%fermi_level,kT,options, &
         jq_ud,jq_du,jq_sym,native_delta,native_curvature,report)
      jq_ud_out=real(jq_ud(1,1,:),rp); jq_du_out=real(jq_du(1,1,:),rp)
      jq_sym_out=real(jq_sym(1,1,:),rp); delta_j=real(native_delta(1,1,:),rp); curvature=real(native_curvature(1,1,:),rp)
      deallocate(jq_ud,jq_du,jq_sym,native_delta,native_curvature)
   end subroutine native_exchange_q_driver

   subroutine write_exchange_q_output(config, lattice_obj, recip, q_direct, q_cart, finite_tt, finite_contact, &
                                      finite_total, native_jq_ud, native_jq_du, native_jq_sym, native_delta_j, native_curvature, native_report, native_ready, endpoint_mode, endpoint_reused, &
                                      q_commensurate, endpoint_residual, endpoint_seconds, assembly_seconds, contraction_seconds, &
                                      hamiltonian_seconds, spectral_tt, spectral_contact, spectral_total, contour_tt, contour_contact, contour_total, &
                                      gf_seconds, solve_seconds, contour_seconds, total_response_seconds)
      type(exchange_q_config), intent(in) :: config
      type(lattice), intent(in) :: lattice_obj
      type(reciprocal), intent(in) :: recip
      real(rp), intent(in) :: q_direct(:, :), q_cart(:, :), finite_tt(:, :, :), finite_contact(:, :, :), finite_total(:, :, :)
      real(rp), intent(in) :: spectral_tt(:, :, :), spectral_contact(:, :, :), spectral_total(:, :, :)
      real(rp), intent(in) :: contour_tt(:, :, :), contour_contact(:, :, :), contour_total(:, :, :)
      real(rp), intent(in), allocatable :: native_jq_ud(:), native_jq_du(:), native_jq_sym(:), native_delta_j(:), native_curvature(:)
      type(native_turek_contour_report), intent(in) :: native_report
      logical, intent(in) :: native_ready
      character(len=*), intent(in) :: endpoint_mode(:)
      logical, intent(in) :: endpoint_reused(:), q_commensurate(:)
      real(rp), intent(in) :: endpoint_residual(:), endpoint_seconds, assembly_seconds, contraction_seconds, hamiltonian_seconds
      real(rp), intent(in) :: gf_seconds, solve_seconds, contour_seconds, total_response_seconds
      integer :: unit, iq, ia, ja, nsite
      real(rp) :: qmag, q2

      open(newunit=unit, file=trim(config%output_file), status='replace', action='write')
      write(unit,'(a)') '# exchange_q static force-theorem exchange curvature'
      write(unit,'(a)') '# observable = DeltaJ(q) = J(Gamma) - J(q)'
      if (native_ready) write(unit,'(a)') '# native_observables = Jud(q), Jdu(q), Jsym(q)=0.5*(Jud+Jdu), DeltaJ_sym=Jsym(Gamma)-Jsym(q), expected_curvature=2*DeltaJ_sym'
      write(unit,'(a)') '# k_mesh performs the electronic BZ integration; q_path is the independent magnetic perturbation path.'
      write(unit,'(a,3(i0,1x))') '# k_mesh = ', recip%nk_mesh
      write(unit,'(a,a)') '# q_coordinates_input = ', trim(config%q_coordinates)
      write(unit,'(a)') '# q_direct is in reciprocal-lattice coordinates; q_cart is in units of 2*pi/alat.'
      write(unit,'(a,a)') '# finite_h_spectral_mode = ', trim(config%finite_h_spectral_mode)
      write(unit,'(a,a)') '# finite_h_response_backend = ', trim(config%finite_h_response_backend)
      if (config%finite_h_response_backend == 'contour' .or. config%finite_h_response_backend == 'both') then
         write(unit,'(a)') '# contour_endpoint_provenance = direct live finite-H fixture assembly at every k and k+q; no endpoint eigenpair reuse'
      end if
      write(unit,'(a,es24.16)') '# electronic_temperature_K = ', recip%temperature
      write(unit,'(a,es24.16)') '# electronic_kT_Ry = ', max(recip%temperature*reciprocal_kb_ry_per_k,1.0e-10_rp)
      write(unit,'(a,es24.16)') '# fermi_level_Ry = ', recip%fermi_level
      write(unit,'(a,a)') '# hamiltonian_order = ', trim(recip%kspace_ham_order)
      if (native_ready) then
         write(unit,'(a)') '# native_state_provenance = accepted reciprocal SCF k points, weights, EF, temperature, lattice, and potential'
         write(unit,'(a,a)') '# native_state_hamiltonian_order = ',trim(recip%kspace_ham_order)
         write(unit,'(a)') '# native_capability = native P-S path operator; no finite-H resolvent or torque vertices consumed'
      end if
      write(unit,'(a)') '# endpoint_mode is mesh_reuse for exact full-mesh commensurate q and explicit_diagonalization otherwise.'
      write(unit,'(a)') '# endpoint_provenance columns: q_index mode reused commensurate coordinate_residual'
      do iq = 1, size(endpoint_mode)
         write(unit,'(a,1x,i0,1x,a,1x,l1,1x,l1,1x,es24.16)') '# endpoint_provenance', iq, trim(endpoint_mode(iq)), &
            endpoint_reused(iq), q_commensurate(iq), endpoint_residual(iq)
      end do
      write(unit,'(a,8(es24.16,1x))') '# timing_seconds endpoint T/C hamiltonian spectral GF_total GF_solve contour total_response = ', &
         endpoint_seconds, assembly_seconds, hamiltonian_seconds, contraction_seconds, gf_seconds, solve_seconds, contour_seconds, total_response_seconds
      write(unit,'(a,a,1x,i0,1x,a,1x,es16.8,1x,a,1x,es16.8)') '# contour_provenance shape=', trim(config%contour_shape), &
         config%contour_points, 'margin=', config%contour_margin, 'height_fraction=', config%contour_height_fraction
      write(unit,'(a,a)') '# fermi_pole_residues_accounted=', merge('true ','false',config%contour_account_fermi_poles)
      if (native_ready) then
         write(unit,'(a,a,1x,i0,1x,a,1x,es16.8,1x,a,1x,es16.8)') '# native_contour_provenance shape=', 'ellipse', &
            config%native_contour_points, 'margin=', config%native_contour_margin, &
            'height_fraction=', config%native_contour_height_fraction
         write(unit,'(a,a)') '# native_fermi_pole_residues_accounted=', &
            merge('true ','false',config%native_contour_account_fermi_poles)
         write(unit,'(a,i0)') '# native_target_fermi_poles = ', config%native_contour_target_fermi_poles
         write(unit,'(a,l1)') '# native_bounds_verified = ', native_report%native_bounds_verified
         write(unit,'(a,es24.16)') '# native_max_spectral_ellipse = ', native_report%native_max_ellipse_value
         write(unit,'(a,i0)') '# native_spectral_poles = ', native_report%native_spectral_poles
      end if
      write(unit,'(a)') '# units: exchange=Ry, q=1/A, stiffness diagnostic=Ry A^2'
      nsite = size(finite_total,1)
      if (nsite == 1) then
         if (config%finite_h_response_backend == 'both') then
            if (native_ready) then
               write(unit,'(a)') '# columns: q_index q1 q2 q3 qx_Ainv qy_Ainv qz_Ainv qmag_Ainv spectral_TT_Ry spectral_contact_Ry spectral_total_Ry contour_TT_Ry contour_contact_Ry contour_total_Ry abs_total_residual_Ry relative_total_residual native_Jud_Ry native_Jdu_Ry native_Jsym_Ry native_DeltaJ_sym_Ry native_DeltaJ_sym_over_q2_RyA2 native_expected_curvature_Ry finiteH_minus_native_curvature_Ry'
            else
               write(unit,'(a)') '# columns: q_index q1 q2 q3 qx_Ainv qy_Ainv qz_Ainv qmag_Ainv spectral_TT_Ry spectral_contact_Ry spectral_total_Ry contour_TT_Ry contour_contact_Ry contour_total_Ry abs_total_residual_Ry relative_total_residual'
            end if
         else if (config%write_components .and. native_ready) then
            write(unit,'(a)') '# columns: q_index q1 q2 q3 qx_Ainv qy_Ainv qz_Ainv qmag_Ainv finiteH_TT_Ry finiteH_contact_Ry finiteH_total_Ry finiteH_total_mRy finiteH_dJ_over_q2_RyA2 native_Jud_Ry native_Jdu_Ry native_Jsym_Ry native_DeltaJ_sym_Ry native_DeltaJ_sym_over_q2_RyA2 native_expected_curvature_Ry finiteH_minus_native_curvature_Ry'
         else if (config%write_components) then
            write(unit,'(a)') '# columns: q_index q1 q2 q3 qx_Ainv qy_Ainv qz_Ainv qmag_Ainv finiteH_TT_Ry finiteH_contact_Ry finiteH_total_Ry finiteH_total_mRy finiteH_dJ_over_q2_RyA2'
         else if (native_ready) then
            write(unit,'(a)') '# columns: q_index q1 q2 q3 qx_Ainv qy_Ainv qz_Ainv qmag_Ainv finiteH_total_Ry finiteH_total_mRy finiteH_dJ_over_q2_RyA2 native_Jud_Ry native_Jdu_Ry native_Jsym_Ry native_DeltaJ_sym_Ry native_DeltaJ_sym_over_q2_RyA2 native_expected_curvature_Ry finiteH_minus_native_curvature_Ry'
         else
            write(unit,'(a)') '# columns: q_index q1 q2 q3 qx_Ainv qy_Ainv qz_Ainv qmag_Ainv finiteH_total_Ry finiteH_total_mRy finiteH_dJ_over_q2_RyA2'
         end if
         do iq = 1, config%n_q
            qmag = 2.0_rp*pi/lattice_obj%alat*sqrt(sum(q_cart(:,iq)**2)); q2 = qmag*qmag
            if (q2 <= tiny(1.0_rp)) q2 = 0.0_rp
            write(unit,'(i0,1x,3(es24.16,1x),3(es24.16,1x),es24.16,1x)', advance='no') &
               iq, q_direct(:,iq), 2.0_rp*pi/lattice_obj%alat*q_cart(:,iq), qmag
            if (config%finite_h_response_backend == 'both') then
               write(unit,'(8(es24.16,1x))', advance='no') spectral_tt(1,1,iq), spectral_contact(1,1,iq), spectral_total(1,1,iq), &
                  contour_tt(1,1,iq), contour_contact(1,1,iq), contour_total(1,1,iq), &
                  abs(spectral_total(1,1,iq)-contour_total(1,1,iq)), &
                  safe_divide(abs(spectral_total(1,1,iq)-contour_total(1,1,iq)), max(abs(spectral_total(1,1,iq)),tiny(1.0_rp)))
            else if (config%write_components) then
               write(unit,'(3(es24.16,1x))', advance='no') finite_tt(1,1,iq), finite_contact(1,1,iq), finite_total(1,1,iq)
            end if
            if (config%finite_h_response_backend /= 'both') then
               if (.not. config%write_components) write(unit,'(es24.16,1x)', advance='no') finite_total(1,1,iq)
               write(unit,'(2(es24.16,1x))', advance='no') 1000.0_rp*finite_total(1,1,iq), safe_divide(finite_total(1,1,iq),q2)
            end if
            if (native_ready) write(unit,'(7(es24.16,1x))', advance='no') native_jq_ud(iq), native_jq_du(iq), native_jq_sym(iq), &
               native_delta_j(iq), safe_divide(native_delta_j(iq),q2), native_curvature(iq), &
               finite_total(1,1,iq)-native_curvature(iq)
            write(unit,*)
         end do
      else
         write(unit,'(a)') '# columns: q_index q1 q2 q3 qx_Ainv qy_Ainv qz_Ainv qmag_Ainv; following columns are row-major finite-H matrices'
         if (config%finite_h_response_backend == 'both') then
            write(unit,'(a,i0,a)') '# matrix columns: spectral TT[1,1..', nsite*nsite, '], spectral contact[...], spectral total[...], contour TT[...], contour contact[...], contour total[...]'
         else if (config%write_components) then
            write(unit,'(a,i0,a)') '# matrix columns: finiteH_TT[1,1..', nsite*nsite, '], finiteH_contact[...], finiteH_total[...]'
         else
            write(unit,'(a,i0,a)') '# matrix columns: finiteH_total[1,1..', nsite*nsite, ']'
         end if
         if (native_ready) write(unit,'(a)') '# optional trailing columns: native_Jud_Ry native_Jdu_Ry native_Jsym_Ry native_DeltaJ_sym_Ry native_expected_curvature_Ry'
         do iq = 1, config%n_q
            qmag = 2.0_rp*pi/lattice_obj%alat*sqrt(sum(q_cart(:,iq)**2))
            write(unit,'(i0,1x,3(es24.16,1x),3(es24.16,1x),es24.16,1x)', advance='no') iq, q_direct(:,iq), &
               2.0_rp*pi/lattice_obj%alat*q_cart(:,iq), qmag
            if (config%finite_h_response_backend == 'both') then
               do ia = 1, nsite
                  do ja = 1, nsite
                     write(unit,'(es24.16,1x)', advance='no') spectral_tt(ia,ja,iq)
                  end do
               end do
               do ia = 1, nsite
                  do ja = 1, nsite
                     write(unit,'(es24.16,1x)', advance='no') spectral_contact(ia,ja,iq)
                  end do
               end do
               do ia = 1, nsite
                  do ja = 1, nsite
                     write(unit,'(es24.16,1x)', advance='no') spectral_total(ia,ja,iq)
                  end do
               end do
               do ia = 1, nsite
                  do ja = 1, nsite
                     write(unit,'(es24.16,1x)', advance='no') contour_tt(ia,ja,iq)
                  end do
               end do
               do ia = 1, nsite
                  do ja = 1, nsite
                     write(unit,'(es24.16,1x)', advance='no') contour_contact(ia,ja,iq)
                  end do
               end do
               do ia = 1, nsite
                  do ja = 1, nsite
                     write(unit,'(es24.16,1x)', advance='no') contour_total(ia,ja,iq)
                  end do
               end do
            else if (config%write_components) then
               do ia = 1, nsite
                  do ja = 1, nsite
                     write(unit,'(es24.16,1x)', advance='no') finite_tt(ia,ja,iq)
                  end do
               end do
               do ia = 1, nsite
                  do ja = 1, nsite
                     write(unit,'(es24.16,1x)', advance='no') finite_contact(ia,ja,iq)
                  end do
               end do
            end if
            do ia = 1, nsite
               do ja = 1, nsite
                  write(unit,'(es24.16,1x)', advance='no') finite_total(ia,ja,iq)
               end do
            end do
            if (native_ready) write(unit,'(5(es24.16,1x))', advance='no') native_jq_ud(iq), native_jq_du(iq), native_jq_sym(iq), &
               native_delta_j(iq), native_curvature(iq)
            write(unit,*)
         end do
      end if
      close(unit)
   end subroutine write_exchange_q_output

   pure function safe_divide(value, denominator) result(result_value)
      real(rp), intent(in) :: value, denominator
      real(rp) :: result_value
      if (denominator <= tiny(1.0_rp)) then
         result_value = 0.0_rp
      else
         result_value = value/denominator
      end if
   end function safe_divide

   pure function elapsed_clock_seconds(start_count, end_count, rate) result(seconds)
      integer, intent(in) :: start_count, end_count, rate
      real(rp) :: seconds
      if (rate > 0) then
         seconds = real(end_count-start_count,rp)/real(rate,rp)
      else
         seconds = 0.0_rp
      end if
   end function elapsed_clock_seconds

   subroutine compute_native_gamma_and_report(native_delta_j, q_direct, gamma)
      real(rp), intent(in) :: native_delta_j(:), q_direct(:, :)
      real(rp), intent(out) :: gamma
      integer :: i
      gamma = 0.0_rp
      do i = 1, size(native_delta_j)
         if (sqrt(sum(q_direct(:,i)**2)) > q_zero_tolerance) gamma = max(gamma,abs(native_delta_j(i)))
      end do
   end subroutine compute_native_gamma_and_report

   pure function real_to_string(value) result(text)
      real(rp), intent(in) :: value
      character(len=48) :: text
      write(text,'(es24.16)') value
   end function real_to_string

end module exchange_q_mod
