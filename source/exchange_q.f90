!------------------------------------------------------------------------------
! Static finite-q exchange-curvature post-processing.
!
! This module deliberately owns the user-facing exchange_q input boundary and
! the production driver only.  The finite-H algebra remains in
! lr_kl_hessian_mod.  The native reciprocal LKAG evaluator below is a
! q-space reference derived from the unchanged dGdG_Jnc contract in
! source/exchange.f90; it is not a replacement for that module.
!------------------------------------------------------------------------------
module exchange_q_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use basis_mod, only: norb
   use control_mod, only: control
   use energy_mod, only: energy
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use lr_kl_hessian_mod, only: lmto_fixture_adapter_residual
   use logger_mod, only: g_logger
   use lr_kl_hessian_mod, only: lmto_live_hamiltonian_fixture, lmto_fixture_from_hamiltonian, &
      assemble_lmto_finite_q_torques, assemble_lmto_finite_q_mixed_derivative, &
      force_theorem_finite_q_hessian_from_eigenbasis_batch
   use math_mod, only: ang2au, inverse_3x3, pi, simpson_f
   use precision_mod, only: rp
   use reciprocal_mod, only: reciprocal
   use self_mod, only: self
   use string_mod, only: sl, int2str, lower
   implicit none

   private

   integer, parameter, public :: exchange_q_max_points = 2000
   real(rp), parameter :: q_zero_tolerance = 1.0e-12_rp

   type, public :: exchange_q_config
      character(len=sl) :: fname = ''
      character(len=16) :: q_coordinates = 'direct'
      character(len=sl) :: q_file = ''
      integer :: n_q = 0
      real(rp), allocatable :: q_list(:, :)
      character(len=sl) :: output_file = 'exchange_q.dat'
      logical :: write_components = .true.
      logical :: native_crosscheck = .false.
      real(rp) :: rotation_axis(3) = [1.0_rp, 0.0_rp, 0.0_rp]
      ! These controls are intentionally diagnostic-only.  They make the
      ! metallic LKAG reference reproducible without changing finite-H data.
      real(rp) :: native_green_eta = 1.0e-3_rp
      integer :: native_energy_points = 0
   contains
      procedure :: clear => exchange_q_config_clear
   end type exchange_q_config

   public :: load_exchange_q_config
   public :: exchange_q_convert_points
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
      this%native_crosscheck = .false.
      this%rotation_axis = [1.0_rp, 0.0_rp, 0.0_rp]
      this%native_green_eta = 1.0e-3_rp
      this%native_energy_points = 0
   end subroutine exchange_q_config_clear

   subroutine load_exchange_q_config(filename, this, validate_request)
      character(len=*), intent(in) :: filename
      type(exchange_q_config), intent(inout) :: this
      logical, intent(in), optional :: validate_request
      logical :: validate
      integer :: iostatus, funit, n_q_points
      character(len=16) :: q_coordinates
      character(len=sl) :: q_file, output_file
      logical :: write_components, native_crosscheck
      real(rp) :: rotation_axis(3), native_green_eta
      integer :: native_energy_points
      real(rp) :: q_list(3, exchange_q_max_points)

      namelist /exchange_q/ q_coordinates, q_file, n_q_points, q_list, output_file, &
         write_components, native_crosscheck, rotation_axis, native_green_eta, native_energy_points

      call this%clear()
      this%fname = filename
      q_coordinates = this%q_coordinates
      q_file = this%q_file
      n_q_points = 0
      q_list = 0.0_rp
      output_file = this%output_file
      write_components = this%write_components
      native_crosscheck = this%native_crosscheck
      rotation_axis = this%rotation_axis
      native_green_eta = this%native_green_eta
      native_energy_points = this%native_energy_points

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
      this%native_crosscheck = native_crosscheck
      this%rotation_axis = rotation_axis
      this%native_green_eta = native_green_eta
      this%native_energy_points = native_energy_points
      validate = .false.
      if (present(validate_request)) validate = validate_request
      if (validate) then
         if (len_trim(this%output_file) == 0) then
            call g_logger%fatal('[exchange_q]: output_file must not be blank', __FILE__, __LINE__)
         end if
         if (this%native_green_eta <= 0.0_rp) then
            call g_logger%fatal('[exchange_q]: native_green_eta must be positive', __FILE__, __LINE__)
         end if
         if (this%native_energy_points < 0) then
            call g_logger%fatal('[exchange_q]: native_energy_points must be non-negative', __FILE__, __LINE__)
         end if
         if (maxval(abs(this%rotation_axis)) <= tiny(1.0_rp)) then
            call g_logger%fatal('[exchange_q]: rotation_axis must be nonzero', __FILE__, __LINE__)
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
      real(rp), allocatable :: evals(:, :), endpoint_evals(:, :), weights(:)
      complex(rp), allocatable :: evecs(:, :, :), endpoint_evecs(:, :, :)
      complex(rp), allocatable :: torques_q(:, :, :, :), torques_minus_q(:, :, :, :), mixed(:, :, :, :, :)
      complex(rp), allocatable :: hessian(:, :), torque_torque(:, :), contact(:, :), complete(:, :)
      real(rp), allocatable :: axes(:, :), finite_total(:, :, :), finite_tt(:, :, :), finite_contact(:, :, :)
      real(rp), allocatable :: native_total(:)
      real(rp) :: adapter_before, adapter_after, axis_norm
      real(rp) :: native_gamma, native_q, native_eta
      integer :: nsite, nmat, nk, iq, ik, ia, ja
      logical :: native_ready

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
      if (adapter_before > 3.0e-12_rp .or. adapter_after > 3.0e-12_rp .or. &
          abs(adapter_after-adapter_before) > 1.0e-13_rp) then
         call g_logger%fatal('[exchange_q]: predls() changed the accepted production-H adapter state', __FILE__, __LINE__)
      end if
      call g_logger%info('[exchange_q]: accepted-H adapter residual before/after='//trim(real_to_string(adapter_before))// &
         '/'//trim(real_to_string(adapter_after)), __FILE__, __LINE__)
      if (abs(self_obj%reciprocal_scf_cache%fermi_level-energy_obj%fermi) > 3.0e-12_rp) then
         call g_logger%fatal('[exchange_q]: accepted Fermi-level provenance mismatch', __FILE__, __LINE__)
      end if

      nsite = fixture%nsite
      nmat = 2*fixture%norb*nsite
      nk = size(reciprocal_obj%eigenvalues,2)
      if (nsite < 1 .or. nk < 1 .or. .not. allocated(reciprocal_obj%eigenvectors)) then
         call g_logger%fatal('[exchange_q]: accepted reciprocal eigensystem is incomplete', __FILE__, __LINE__)
      end if
      if (size(reciprocal_obj%eigenvalues,1) /= nmat) then
         call g_logger%fatal('[exchange_q]: accepted eigensystem/Hamiltonian dimensions disagree', __FILE__, __LINE__)
      end if
      allocate(evals(nmat,nk), evecs(nmat,nmat,nk), weights(nk))
      evals = reciprocal_obj%eigenvalues
      evecs = reciprocal_obj%eigenvectors
      weights = reciprocal_obj%k_weights
      if (size(weights) /= nk .or. sum(weights) <= tiny(1.0_rp)) then
         call g_logger%fatal('[exchange_q]: accepted k-mesh weights are invalid', __FILE__, __LINE__)
      end if

      allocate(finite_total(nsite,nsite,config%n_q), finite_tt(nsite,nsite,config%n_q), &
         finite_contact(nsite,nsite,config%n_q))
      finite_total = 0.0_rp; finite_tt = 0.0_rp; finite_contact = 0.0_rp
      allocate(hessian(nsite,nsite), torque_torque(nsite,nsite), contact(nsite,nsite), complete(nsite,nsite))
      native_ready = config%native_crosscheck
      allocate(native_total(config%n_q)); native_total = 0.0_rp

      do iq = 1, config%n_q
         if (sqrt(sum(q_direct(:,iq)**2)) <= q_zero_tolerance .and. nsite == 1) then
            ! For the one-site production state the accepted scalar-relativistic
            ! collinear state has the exact global-rotation Goldstone reduction.
            ! Avoid a Fermi-edge spectral denominator at Gamma; the independent
            ! finite-difference and q=0 unit gates cover the API reduction.
            ! Multi-sublattice Gamma is evaluated below so its full Hessian
            ! matrix remains available (including optic/acoustic row sums).
            finite_total(:,:,iq) = 0.0_rp
            finite_tt(:,:,iq) = 0.0_rp
            finite_contact(:,:,iq) = 0.0_rp
         else
            call solve_unfolded_endpoints(reciprocal_obj, reciprocal_obj%k_points + spread_q(q_direct(:,iq),nk), &
                                          endpoint_evals, endpoint_evecs)
            allocate(torques_q(nmat,nmat,nsite,nk), torques_minus_q(nmat,nmat,nsite,nk), &
                     mixed(nmat,nmat,nsite,nsite,nk))
            do ik = 1, nk
               call assemble_lmto_finite_q_torques(fixture, reciprocal_obj%k_points(:,ik), q_direct(:,iq), axes, &
                  torques_q(:,:,:,ik))
               call assemble_lmto_finite_q_torques(fixture, reciprocal_obj%k_points(:,ik)+q_direct(:,iq), -q_direct(:,iq), &
                  axes, torques_minus_q(:,:,:,ik))
               do ia = 1, nsite
                  do ja = 1, nsite
                     call assemble_lmto_finite_q_mixed_derivative(fixture, reciprocal_obj%k_points(:,ik), q_direct(:,iq), &
                        ia, axes(:,ia), ja, axes(:,ja), mixed(:,:,ia,ja,ik))
                  end do
               end do
            end do
            call force_theorem_finite_q_hessian_from_eigenbasis_batch(evals, evecs, endpoint_evals, endpoint_evecs, &
               reciprocal_obj%fermi_level, weights, torques_q, torques_minus_q, mixed, hessian, torque_torque, contact, complete)
            finite_tt(:,:,iq) = real(torque_torque,rp)
            finite_contact(:,:,iq) = real(contact,rp)
            finite_total(:,:,iq) = real(complete,rp)
            deallocate(torques_q, torques_minus_q, mixed, endpoint_evals, endpoint_evecs)
         end if
      end do

      if (native_ready) then
         if (config%native_energy_points > 0) energy_obj%channels_ldos = config%native_energy_points
         call energy_obj%e_mesh()
         native_eta = config%native_green_eta
         call lkag_native_path(reciprocal_obj, lattice_obj, energy_obj, q_direct, native_eta, native_total)
      end if
      call write_exchange_q_output(config, lattice_obj, reciprocal_obj, q_direct, q_cart, finite_tt, finite_contact, &
         finite_total, native_total, native_ready)
      if (native_ready) then
         call compute_native_gamma_and_report(native_total, q_direct, native_gamma)
         native_q = maxval(abs(native_total))
         call g_logger%info('[exchange_q]: native LKAG-q completed; |dJ|_max='//trim(real_to_string(native_q))// &
            ' Ry, gamma='//trim(real_to_string(native_gamma))//' Ry', __FILE__, __LINE__)
      end if
      call fixture%clear()
      deallocate(q_direct, q_cart, evals, evecs, weights, finite_total, finite_tt, finite_contact, &
         axes, hessian, torque_torque, contact, complete)
      deallocate(native_total)
   end subroutine run_exchange_q

   subroutine validate_exchange_q_capability(config, control_obj, ham, self_obj, recip)
      type(exchange_q_config), intent(in) :: config
      type(control), intent(in) :: control_obj
      type(hamiltonian), intent(in) :: ham
      type(self), intent(in) :: self_obj
      type(reciprocal), intent(in) :: recip

      if (config%n_q < 2) call g_logger%fatal('[exchange_q]: q path is empty or has no finite-q point', __FILE__, __LINE__)
      if (trim(control_obj%calctype) /= 'B' .or. control_obj%nsp /= 1 .or. control_obj%has_soc()) then
         call g_logger%fatal('[exchange_q]: capability gate requires bulk nsp=1 scalar-relativistic collinear state with SOC off', &
            __FILE__, __LINE__)
      end if
      if (.not. self_obj%use_kspace .or. .not. allocated(self_obj%reciprocal_scf_cache)) then
         call g_logger%fatal('[exchange_q]: requires self%use_kspace=.true. to consume the accepted reciprocal SCF state', &
            __FILE__, __LINE__)
      end if
      if (ham%hoh .or. ham%ccor_2c .or. ham%hubbard_u_general_check .or. ham%hubbard_v_check .or. &
          control_obj%constraints_enable) then
         call g_logger%fatal('[exchange_q]: unsupported HOH/CCOR/Hubbard/constraint combination', __FILE__, __LINE__)
      end if
      if (trim(recip%reciprocal_mode) /= 'ham_only' .or. trim(recip%kspace_ham_order) /= 'first') then
         call g_logger%fatal('[exchange_q]: capability gate requires reciprocal_mode=ham_only and kspace_ham_order=first', &
            __FILE__, __LINE__)
      end if
      if (fixture_basis_size(ham) /= 9) then
         call g_logger%fatal('[exchange_q]: capability gate requires the full spd production basis', __FILE__, __LINE__)
      end if
   end subroutine validate_exchange_q_capability

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

   subroutine lkag_native_path(recip, lattice_obj, energy_obj, q_direct, eta, delta_j)
      type(reciprocal), intent(inout) :: recip
      type(lattice), intent(in) :: lattice_obj
      type(energy), intent(in) :: energy_obj
      real(rp), intent(in) :: q_direct(:, :), eta
      real(rp), intent(out) :: delta_j(:)
      real(rp), allocatable :: path_j(:), q_j(:, :)
      integer :: iq

      allocate(path_j(size(q_direct,2)), q_j(3,1))
      do iq = 1, size(q_direct,2)
         q_j(:,1) = q_direct(:,iq)
         call lkag_native_q(recip, lattice_obj, energy_obj, q_j, eta, path_j(iq))
      end do
      delta_j = (path_j(1)-path_j)
      deallocate(path_j,q_j)
   end subroutine lkag_native_path

   subroutine lkag_native_q(recip, lattice_obj, energy_obj, q_direct, eta, value)
      type(reciprocal), intent(inout) :: recip
      type(lattice), intent(in) :: lattice_obj
      type(energy), intent(in) :: energy_obj
      real(rp), intent(in) :: q_direct(:, :), eta
      real(rp), intent(out) :: value
      real(rp), allocatable :: endpoint_evals(:, :), integrand(:, :)
      complex(rp), allocatable :: endpoint_evecs(:, :, :), green_left(:, :), green_right(:, :)
      complex(rp), allocatable :: weighted_left(:, :), weighted_right(:, :)
      complex(rp), allocatable :: di(:, :, :)
      complex(rp), allocatable :: gni(:, :), gxi(:, :), gyi(:, :), gzi(:, :)
      complex(rp), allocatable :: gnj(:, :), gxj(:, :), gyj(:, :), gzj(:, :)
      complex(rp) :: z
      integer :: nsite, nmat, nk, ne, ie, ik, ia, ja
      real(rp) :: weighted, qnorm

      nsite = lattice_obj%nrec
      nmat = size(recip%eigenvalues,1)
      nk = size(recip%eigenvalues,2)
      ne = size(energy_obj%ene)
      qnorm = sqrt(sum(q_direct(:,1)**2))
      if (qnorm <= q_zero_tolerance) then
         endpoint_evals = recip%eigenvalues
         endpoint_evecs = recip%eigenvectors
      else
         call solve_unfolded_endpoints(recip, recip%k_points + spread_q(q_direct(:,1),nk), endpoint_evals, endpoint_evecs)
      end if
      allocate(integrand(ne,nsite*nsite), di(norb,norb,nsite), &
               green_left(nmat,nmat), green_right(nmat,nmat), weighted_left(nmat,nmat), weighted_right(nmat,nmat), &
               gni(norb,norb), gxi(norb,norb), &
               gyi(norb,norb), gzi(norb,norb), gnj(norb,norb), gxj(norb,norb), gyj(norb,norb), gzj(norb,norb))
      integrand = 0.0_rp
      do ie = 1, ne
         z = cmplx(energy_obj%ene(ie), eta, rp)
         do ia = 1, nsite
            call lattice_obj%symbolic_atoms(lattice_obj%nbulk+ia)%d_matrix(di(:,:,ia), energy_obj%ene(ie))
         end do
         do ik = 1, nk
            call spectral_green(recip%eigenvalues(:,ik), recip%eigenvectors(:,:,ik), z, green_left, weighted_left)
            call spectral_green(endpoint_evals(:,ik), endpoint_evecs(:,:,ik), z, green_right, weighted_right)
            do ia = 1, nsite
               do ja = 1, nsite
                  call pauli_site_block(green_left, ia, ja, gni, gxi, gyi, gzi)
                  call pauli_site_block(green_right, ja, ia, gnj, gxj, gyj, gzj)
                  weighted = lkag_pauli_product(di(:,:,ia), di(:,:,ja), gni, gxi, gyi, gzi, gnj, gxj, gyj, gzj)
                  integrand(ie, (ia-1)*nsite+ja) = integrand(ie, (ia-1)*nsite+ja) + &
                     recip%k_weights(ik)*weighted
               end do
            end do
         end do
      end do
      value = 0.0_rp
      do ia = 1, nsite*nsite
         call simpson_f(weighted, energy_obj%ene, energy_obj%fermi, energy_obj%nv1, integrand(:,ia), &
            .true., .false., 0.0_rp)
         value = value + weighted/real(sum(recip%k_weights),rp)/(4.0_rp*pi)
      end do
      deallocate(endpoint_evals,endpoint_evecs,integrand,di,green_left,green_right,weighted_left,weighted_right, &
         gni,gxi,gyi,gzi,gnj,gxj,gyj,gzj)
   end subroutine lkag_native_q

   subroutine spectral_green(eigenvalues, eigenvectors, z, green, weighted_vectors)
      real(rp), intent(in) :: eigenvalues(:)
      complex(rp), intent(in) :: eigenvectors(:, :), z
      complex(rp), intent(out) :: green(:, :), weighted_vectors(:, :)
      complex(rp) :: alpha, beta
      integer :: i, j, n
      external :: zgemm
      n = size(eigenvalues)
      do j = 1, n
         do i = 1, n
            weighted_vectors(i,j) = eigenvectors(i,j)/ (z-eigenvalues(j))
         end do
      end do
      alpha = cmplx(1.0_rp,0.0_rp,rp)
      beta = cmplx(0.0_rp,0.0_rp,rp)
      call zgemm('N','C',n,n,n,alpha,weighted_vectors,n,eigenvectors,n,beta,green,n)
   end subroutine spectral_green

   subroutine pauli_site_block(green, site_i, site_j, gn, gx, gy, gz)
      complex(rp), intent(in) :: green(:, :)
      integer, intent(in) :: site_i, site_j
      complex(rp), intent(out) :: gn(:, :), gx(:, :), gy(:, :), gz(:, :)
      integer :: i, j, n, no
      integer :: i0, j0, io, jo
      no = size(gn,1)
      i0 = (site_i-1)*2*no; j0 = (site_j-1)*2*no
      do j = 1, no
         do i = 1, no
            io = i0+i; jo = j0+j
            gn(i,j) = 0.5_rp*(green(io,jo)+green(io+no,jo+no))
            gz(i,j) = 0.5_rp*(green(io,jo)-green(io+no,jo+no))
            gy(i,j) = 0.5_rp*cmplx(0.0_rp,1.0_rp,rp)*(green(io,jo+no)-green(io+no,jo))
            gx(i,j) = 0.5_rp*(green(io,jo+no)+green(io+no,jo))
         end do
      end do
   end subroutine pauli_site_block

   subroutine write_exchange_q_output(config, lattice_obj, recip, q_direct, q_cart, finite_tt, finite_contact, &
                                      finite_total, native_total, native_ready)
      type(exchange_q_config), intent(in) :: config
      type(lattice), intent(in) :: lattice_obj
      type(reciprocal), intent(in) :: recip
      real(rp), intent(in) :: q_direct(:, :), q_cart(:, :), finite_tt(:, :, :), finite_contact(:, :, :), finite_total(:, :, :)
      real(rp), intent(in), allocatable :: native_total(:)
      logical, intent(in) :: native_ready
      integer :: unit, iq, ia, ja, nsite
      real(rp) :: qmag, q2

      open(newunit=unit, file=trim(config%output_file), status='replace', action='write')
      write(unit,'(a)') '# exchange_q static force-theorem exchange curvature'
      write(unit,'(a)') '# observable = DeltaJ(q) = J(Gamma) - J(q)'
      write(unit,'(a)') '# k_mesh performs the electronic BZ integration; q_path is the independent magnetic perturbation path.'
      write(unit,'(a,3(i0,1x))') '# k_mesh = ', recip%nk_mesh
      write(unit,'(a,a)') '# q_coordinates_input = ', trim(config%q_coordinates)
      write(unit,'(a)') '# q_direct is in reciprocal-lattice coordinates; q_cart is in units of 2*pi/alat.'
      write(unit,'(a)') '# q values are not rounded to the k mesh; endpoint diagonalization uses the exact k+q Hamiltonian.'
      write(unit,'(a)') '# units: exchange=Ry, q=1/A, stiffness diagnostic=Ry A^2'
      nsite = size(finite_total,1)
      if (nsite == 1) then
         if (config%write_components .and. native_ready) then
            write(unit,'(a)') '# columns: q_index q1 q2 q3 qx_Ainv qy_Ainv qz_Ainv qmag_Ainv finiteH_TT_Ry finiteH_contact_Ry finiteH_total_Ry finiteH_total_mRy finiteH_dJ_over_q2_RyA2 native_dJ_Ry native_dJ_mRy native_dJ_over_q2_RyA2 finiteH_minus_native_Ry'
         else if (config%write_components) then
            write(unit,'(a)') '# columns: q_index q1 q2 q3 qx_Ainv qy_Ainv qz_Ainv qmag_Ainv finiteH_TT_Ry finiteH_contact_Ry finiteH_total_Ry finiteH_total_mRy finiteH_dJ_over_q2_RyA2'
         else if (native_ready) then
            write(unit,'(a)') '# columns: q_index q1 q2 q3 qx_Ainv qy_Ainv qz_Ainv qmag_Ainv finiteH_total_Ry finiteH_total_mRy finiteH_dJ_over_q2_RyA2 native_dJ_Ry native_dJ_mRy native_dJ_over_q2_RyA2 finiteH_minus_native_Ry'
         else
            write(unit,'(a)') '# columns: q_index q1 q2 q3 qx_Ainv qy_Ainv qz_Ainv qmag_Ainv finiteH_total_Ry finiteH_total_mRy finiteH_dJ_over_q2_RyA2'
         end if
         do iq = 1, config%n_q
            qmag = 2.0_rp*pi/lattice_obj%alat*sqrt(sum(q_cart(:,iq)**2)); q2 = qmag*qmag
            if (q2 <= tiny(1.0_rp)) q2 = 0.0_rp
            write(unit,'(i0,1x,3(es24.16,1x),3(es24.16,1x),es24.16,1x)', advance='no') &
               iq, q_direct(:,iq), 2.0_rp*pi/lattice_obj%alat*q_cart(:,iq), qmag
            if (config%write_components) write(unit,'(3(es24.16,1x))', advance='no') finite_tt(1,1,iq), finite_contact(1,1,iq), finite_total(1,1,iq)
            if (.not. config%write_components) write(unit,'(es24.16,1x)', advance='no') finite_total(1,1,iq)
            write(unit,'(2(es24.16,1x))', advance='no') 1000.0_rp*finite_total(1,1,iq), safe_divide(finite_total(1,1,iq),q2)
            if (native_ready) write(unit,'(4(es24.16,1x))', advance='no') native_total(iq), 1000.0_rp*native_total(iq), &
               safe_divide(native_total(iq),q2), finite_total(1,1,iq)-native_total(iq)
            write(unit,*)
         end do
      else
         write(unit,'(a)') '# columns: q_index q1 q2 q3 qx_Ainv qy_Ainv qz_Ainv qmag_Ainv; following columns are row-major finite-H matrices'
         if (config%write_components) then
            write(unit,'(a,i0,a)') '# matrix columns: finiteH_TT[1,1..', nsite*nsite, '], finiteH_contact[...], finiteH_total[...]'
         else
            write(unit,'(a,i0,a)') '# matrix columns: finiteH_total[1,1..', nsite*nsite, ']'
         end if
         if (native_ready) write(unit,'(a)') '# optional trailing column: native_qspace_dJ_Ry (ordered site-pair sum)'
         do iq = 1, config%n_q
            qmag = 2.0_rp*pi/lattice_obj%alat*sqrt(sum(q_cart(:,iq)**2))
            write(unit,'(i0,1x,3(es24.16,1x),3(es24.16,1x),es24.16,1x)', advance='no') iq, q_direct(:,iq), &
               2.0_rp*pi/lattice_obj%alat*q_cart(:,iq), qmag
            if (config%write_components) then
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
            if (native_ready) write(unit,'(es24.16,1x)', advance='no') native_total(iq)
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

   subroutine compute_native_gamma_and_report(native_total, q_direct, gamma)
      real(rp), intent(in) :: native_total(:), q_direct(:, :)
      real(rp), intent(out) :: gamma
      integer :: i
      gamma = 0.0_rp
      do i = 1, size(native_total)
         if (sqrt(sum(q_direct(:,i)**2)) > q_zero_tolerance) gamma = max(gamma,abs(native_total(i)))
      end do
   end subroutine compute_native_gamma_and_report

   pure function real_to_string(value) result(text)
      real(rp), intent(in) :: value
      character(len=48) :: text
      write(text,'(es24.16)') value
   end function real_to_string

end module exchange_q_mod
