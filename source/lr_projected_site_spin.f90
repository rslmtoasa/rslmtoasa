!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief DRESP-01 projected site-spin operator and moment contract.
!>
!> This module is deliberately downstream of the accepted radial and product
!> representations.  It owns the projection selector, the site-integration
!> functional, and scalar spin-moment audits.  It does not construct a
!> susceptibility, exchange kernel, or Dyson problem.
!>
!> The selector is an orbital selector: `d` means l=2 and `spd` means
!> l=0,1,2.  It is never a response-harmonic selector.  Every Gaunt-allowed
!> response (L,M) product harmonic remains in the product basis; site
!> integration subsequently selects only L=0,M=0.
!------------------------------------------------------------------------------
module lr_projected_site_spin_mod

   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use radial_ground_state_mod, only: radial_ground_state, radial_simpson_weight
   use response_angular_basis_mod, only: response_angular_pi, response_harmonic
   use lr_response_space_mod, only: response_space_layout
   use lr_pauli_transition_vertex_mod, only: pauli_endpoint_state
   use lr_lmto_product_response_basis_mod, only: lmto_product_channel_plus, &
      lmto_product_channel_minus, lmto_product_response_basis
   implicit none
   private

   real(rp), parameter :: kB_Ry_per_K = 6.3336814e-6_rp
   real(rp), parameter :: occupation_kT_floor = 1.0e-10_rp
   real(rp), parameter :: mesh_tolerance = 2.0e-13_rp

   integer, parameter, public :: projected_selector_d = 1
   integer, parameter, public :: projected_selector_spd = 2

   integer, parameter, public :: projected_operator_plus = 1
   integer, parameter, public :: projected_operator_minus = 2
   integer, parameter, public :: projected_operator_z = 3

   character(len=1), parameter, public :: projected_selection_d = 'd'
   character(len=3), parameter, public :: projected_selection_spd = 'spd'
   character(len=40), parameter, public :: projected_core_policy = &
      'valence-only; frozen core excluded'

   !> Immutable-after-initialization DRESP-01 contract metadata and operations.
   type, public :: projected_site_spin_contract
      integer :: nsite = 0
      integer :: orbital_lmax = -1
      integer :: response_lmax = -1
      integer :: npoint = 0
      integer :: selector_kind = 0
      character(len=8) :: selector = ''
      character(len=40) :: core_policy = projected_core_policy
      real(rp) :: mesh_a = 0.0_rp
      real(rp) :: mesh_b = 0.0_rp
      real(rp) :: angular_scalar_integral = 0.0_rp
      real(rp), allocatable :: radius(:)
      real(rp), allocatable :: radial_weights(:)
      logical :: selected_l(0:2) = .false.
   contains
      procedure :: initialize => projected_site_spin_initialize
      procedure :: site_integration_functional => projected_site_spin_functional
      procedure :: selected_transition_coordinates => projected_selected_coordinates
      procedure :: transition_amplitudes => projected_transition_amplitudes
      procedure :: direct_operator_matrix => projected_direct_operator_matrix
      procedure :: direct_transition_amplitudes => projected_direct_transition_amplitudes
      procedure :: moment_from_density => projected_moment_from_density
      procedure :: moment_from_operator => projected_moment_from_operator
      procedure :: core_spin_number => projected_core_spin_number
   end type projected_site_spin_contract

   public :: projected_selector_supported

contains

   pure logical function projected_selector_supported(selection) result(ok)
      character(len=*), intent(in) :: selection

      select case (trim(selection))
      case (projected_selection_d, projected_selection_spd)
         ok = .true.
      case default
         ok = .false.
      end select
   end function projected_selector_supported

   subroutine projected_site_spin_initialize(this, space, radial_bases, selection)
      class(projected_site_spin_contract), intent(out) :: this
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      character(len=*), intent(in) :: selection

      if (.not. projected_selector_supported(selection)) then
         error stop 'DRESP-01: selector must be exactly d or spd; spdf is unsupported'
      end if
      if (space%nsite < 1 .or. space%nchannel /= 1 .or. space%response_lmax < 4 .or. &
          space%npoint < 3 .or. .not. allocated(space%radius) .or. &
          .not. allocated(space%radial_weights)) then
         error stop 'DRESP-01: response space is not a complete spd one-channel layout'
      end if
      if (size(radial_bases) /= space%nsite) then
         error stop 'DRESP-01: one complete radial basis is required per site'
      end if

      this%nsite = space%nsite
      this%orbital_lmax = radial_bases(1)%lmax
      this%response_lmax = space%response_lmax
      this%npoint = space%npoint
      this%mesh_a = space%a
      this%mesh_b = space%b
      this%selector = trim(selection)
      this%core_policy = projected_core_policy
      this%selected_l = .false.
      select case (trim(selection))
      case (projected_selection_d)
         this%selector_kind = projected_selector_d
         this%selected_l(2) = .true.
      case (projected_selection_spd)
         this%selector_kind = projected_selector_spd
         this%selected_l = .true.
      end select

      if (this%orbital_lmax /= 2 .or. this%response_lmax /= 2*this%orbital_lmax) then
         error stop 'DRESP-01: only the accepted lmax=2, response_lmax=4 spd contract is supported'
      end if
      allocate(this%radius(this%npoint), this%radial_weights(this%npoint))
      this%radius = space%radius
      this%radial_weights = space%radial_weights
      if (any(.not. ieee_is_finite(this%radial_weights)) .or. any(this%radial_weights < 0.0_rp)) then
         error stop 'DRESP-01: invalid response-space metric weights'
      end if
      this%angular_scalar_integral = 4.0_rp*response_angular_pi* &
         real(response_harmonic(0, 0, 0.0_rp, 0.0_rp), rp)
      if (.not. ieee_is_finite(this%angular_scalar_integral) .or. this%angular_scalar_integral <= 0.0_rp) then
         error stop 'DRESP-01: invalid Y00 full-sphere normalization'
      end if
      call validate_radial_bases(this, radial_bases)
   end subroutine projected_site_spin_initialize

   subroutine projected_site_spin_functional(this, product, functionals)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(out) :: functionals(:, :)

      integer :: flat, site, response_l, response_m, product_mode

      call validate_product(this, product, 'DRESP-01 site integration')
      if (any(shape(functionals) /= [product%product_dimension, this%nsite])) then
         error stop 'DRESP-01 site integration: functional shape mismatch'
      end if
      functionals = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, product%product_dimension
         call product%unflatten_index(flat, site, response_l, response_m, product_mode)
         if (response_l /= 0 .or. response_m /= 0) cycle
         ! The row integral is a^T z.  Store p=conjg(a), so p^H z is the
         ! Fortran DOT_PRODUCT(functional,z) contraction below.
         functionals(flat, site) = conjg(cmplx(this%angular_scalar_integral, 0.0_rp, rp)* &
            sum(sqrt(this%radial_weights)*product%blocks(site, response_l)%weighted_modes(:, product_mode)))
      end do
   end subroutine projected_site_spin_functional

   subroutine projected_selected_coordinates(this, product, left_state, right_state, coordinates)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_product_response_basis), intent(in) :: product
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(out) :: coordinates(:)

      complex(rp), allocatable :: coefficients(:)
      integer :: site, response_l, response_m, k, first, last

      call validate_product(this, product, 'DRESP-01 selected coordinates')
      if (size(coordinates) /= product%product_dimension) then
         error stop 'DRESP-01 selected coordinates: coordinate shape mismatch'
      end if
      coordinates = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, this%nsite
         do response_l = 0, product%response_lmax
            allocate(coefficients(product%blocks(site, response_l)%ncandidate))
            do response_m = -response_l, response_l
               call product%candidate_coefficients(site, response_l, response_m, left_state, right_state, coefficients)
               do k = 1, size(coefficients)
                  if (.not. this%selected_l(product%blocks(site, response_l)%candidates(k)%l) .or. &
                      .not. this%selected_l(product%blocks(site, response_l)%candidates(k)%lp)) then
                     coefficients(k) = cmplx(0.0_rp, 0.0_rp, rp)
                  end if
               end do
               first = product%flat_index(site, response_l, response_m, 1)
               last = product%flat_index(site, response_l, response_m, product%blocks(site, response_l)%rank)
               coordinates(first:last) = matmul(product%blocks(site, response_l)%forward_transform, coefficients)
            end do
            deallocate(coefficients)
         end do
      end do
   end subroutine projected_selected_coordinates

   subroutine projected_transition_amplitudes(this, product, left_state, right_state, amplitudes)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_product_response_basis), intent(in) :: product
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(out) :: amplitudes(:)

      complex(rp), allocatable :: coordinates(:), functionals(:, :)
      integer :: site

      call validate_product(this, product, 'DRESP-01 transition amplitudes')
      if (size(amplitudes) /= this%nsite) then
         error stop 'DRESP-01 transition amplitudes: site amplitude shape mismatch'
      end if
      allocate(coordinates(product%product_dimension), functionals(product%product_dimension, this%nsite))
      call this%selected_transition_coordinates(product, left_state, right_state, coordinates)
      call this%site_integration_functional(product, functionals)
      amplitudes = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, this%nsite
         amplitudes(site) = dot_product(functionals(:, site), coordinates)
      end do
   end subroutine projected_transition_amplitudes

   subroutine projected_direct_operator_matrix(this, radial_bases, energy, operator_kind, matrix)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      real(rp), intent(in) :: energy
      integer, intent(in) :: operator_kind
      complex(rp), intent(out) :: matrix(:, :)

      integer :: norb, site, iorb, l, up, down, offset
      real(rp) :: overlap

      call validate_radial_bases(this, radial_bases)
      call validate_operator_kind(operator_kind, 'DRESP-01 direct operator')
      norb = (this%orbital_lmax + 1)**2
      if (any(shape(matrix) /= [2*norb*this%nsite, 2*norb*this%nsite])) then
         error stop 'DRESP-01 direct operator: matrix shape mismatch'
      end if
      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, this%nsite
         offset = (site - 1)*2*norb
         do iorb = 1, norb
            l = lmto_orbital_l(iorb)
            if (.not. this%selected_l(l)) cycle
            up = offset + iorb
            down = offset + norb + iorb
            select case (operator_kind)
            case (projected_operator_plus)
               overlap = endpoint_radial_overlap(radial_bases(site), l, 1, energy, l, 2, energy)
               matrix(up, down) = cmplx(overlap, 0.0_rp, rp)
            case (projected_operator_minus)
               overlap = endpoint_radial_overlap(radial_bases(site), l, 2, energy, l, 1, energy)
               matrix(down, up) = cmplx(overlap, 0.0_rp, rp)
            case (projected_operator_z)
               matrix(up, up) = cmplx(endpoint_radial_overlap(radial_bases(site), l, 1, energy, l, 1, energy), 0.0_rp, rp)
               matrix(down, down) = cmplx(-endpoint_radial_overlap(radial_bases(site), l, 2, energy, l, 2, energy), &
                                          0.0_rp, rp)
            end select
         end do
      end do
   end subroutine projected_direct_operator_matrix

   subroutine projected_direct_transition_amplitudes(this, radial_bases, left_state, right_state, operator_kind, amplitudes)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      integer, intent(in) :: operator_kind
      complex(rp), intent(out) :: amplitudes(:)

      integer :: norb, site, iorb, l, offset
      complex(rp) :: left_coefficient, right_coefficient

      call validate_radial_bases(this, radial_bases)
      call validate_operator_kind(operator_kind, 'DRESP-01 direct transition')
      norb = (this%orbital_lmax + 1)**2
      call validate_endpoint_shapes(this%nsite, norb, left_state, right_state, 'DRESP-01 direct transition')
      if (size(amplitudes) /= this%nsite) then
         error stop 'DRESP-01 direct transition: site amplitude shape mismatch'
      end if
      amplitudes = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, this%nsite
         offset = (site - 1)*2*norb
         do iorb = 1, norb
            l = lmto_orbital_l(iorb)
            if (.not. this%selected_l(l)) cycle
            select case (operator_kind)
            case (projected_operator_plus)
               left_coefficient = left_state%coefficients(offset + iorb)
               right_coefficient = right_state%coefficients(offset + norb + iorb)
               amplitudes(site) = amplitudes(site) + conjg(left_coefficient)*right_coefficient* &
                  endpoint_radial_overlap(radial_bases(site), l, 1, left_state%energy, l, 2, right_state%energy)
            case (projected_operator_minus)
               left_coefficient = left_state%coefficients(offset + norb + iorb)
               right_coefficient = right_state%coefficients(offset + iorb)
               amplitudes(site) = amplitudes(site) + conjg(left_coefficient)*right_coefficient* &
                  endpoint_radial_overlap(radial_bases(site), l, 2, left_state%energy, l, 1, right_state%energy)
            case (projected_operator_z)
               left_coefficient = left_state%coefficients(offset + iorb)
               right_coefficient = right_state%coefficients(offset + iorb)
               amplitudes(site) = amplitudes(site) + conjg(left_coefficient)*right_coefficient* &
                  endpoint_radial_overlap(radial_bases(site), l, 1, left_state%energy, l, 1, right_state%energy)
               left_coefficient = left_state%coefficients(offset + norb + iorb)
               right_coefficient = right_state%coefficients(offset + norb + iorb)
               amplitudes(site) = amplitudes(site) - conjg(left_coefficient)*right_coefficient* &
                  endpoint_radial_overlap(radial_bases(site), l, 2, left_state%energy, l, 2, right_state%energy)
            end select
         end do
      end do
   end subroutine projected_direct_transition_amplitudes

   subroutine projected_moment_from_density(this, eigenvalues, eigenvectors, k_weights, fermi_level, temperature, &
                                            ground_states, moment)
      class(projected_site_spin_contract), intent(in) :: this
      real(rp), intent(in) :: eigenvalues(:, :), k_weights(:), fermi_level, temperature
      complex(rp), intent(in) :: eigenvectors(:, :, :)
      type(radial_ground_state), intent(in) :: ground_states(:)
      real(rp), intent(out) :: moment(:)

      integer :: site, ik, ib, ir, spin, iorb, l, norb, nbands, nk, offset
      real(rp) :: wk, occupation, coefficient_weight, u_value
      real(rp), allocatable :: profile(:, :)

      call validate_moment_inputs(this, eigenvalues, eigenvectors, k_weights, fermi_level, temperature, &
                                  ground_states, moment)
      norb = (this%orbital_lmax + 1)**2
      nbands = size(eigenvalues, 1)
      nk = size(eigenvalues, 2)
      moment = 0.0_rp
      do site = 1, this%nsite
         allocate(profile(this%npoint, 2))
         profile = 0.0_rp
         do ik = 1, nk
            wk = k_weights(ik)/sum(k_weights)
            do ib = 1, nbands
               occupation = fermi_dirac_occupation(eigenvalues(ib, ik), fermi_level, temperature)
               if (occupation == 0.0_rp) cycle
               offset = (site - 1)*2*norb
               do spin = 1, 2
                  do iorb = 1, norb
                     l = lmto_orbital_l(iorb)
                     if (.not. this%selected_l(l)) cycle
                     coefficient_weight = real(eigenvectors(offset + (spin - 1)*norb + iorb, ib, ik)* &
                        conjg(eigenvectors(offset + (spin - 1)*norb + iorb, ib, ik)), rp)
                     do ir = 1, this%npoint
                        u_value = ground_states(site)%pauli_large(ir, l + 1, spin) + &
                           (eigenvalues(ib, ik) - ground_states(site)%pauli_enu(l + 1, spin))* &
                           ground_states(site)%pauli_large_dot(ir, l + 1, spin)
                        profile(ir, spin) = profile(ir, spin) + wk*occupation*coefficient_weight*u_value**2
                     end do
                  end do
               end do
            end do
         end do
         moment(site) = ground_states(site)%log_mesh_integral(profile(:, 1) - profile(:, 2))
         deallocate(profile)
      end do
   end subroutine projected_moment_from_density

   subroutine projected_moment_from_operator(this, eigenvalues, eigenvectors, k_weights, fermi_level, temperature, &
                                             ground_states, moment)
      class(projected_site_spin_contract), intent(in) :: this
      real(rp), intent(in) :: eigenvalues(:, :), k_weights(:), fermi_level, temperature
      complex(rp), intent(in) :: eigenvectors(:, :, :)
      type(radial_ground_state), intent(in) :: ground_states(:)
      real(rp), intent(out) :: moment(:)

      complex(rp), allocatable :: local_state(:), local_operator(:, :), action(:)
      integer :: norb, nbands, nk, site, ik, ib, iorb, l, spin, offset
      real(rp) :: wk, occupation, expectation

      call validate_moment_inputs(this, eigenvalues, eigenvectors, k_weights, fermi_level, temperature, &
                                  ground_states, moment)
      norb = (this%orbital_lmax + 1)**2
      nbands = size(eigenvalues, 1)
      nk = size(eigenvalues, 2)
      moment = 0.0_rp
      allocate(local_state(2*norb), local_operator(2*norb, 2*norb), action(2*norb))
      do ik = 1, nk
         wk = k_weights(ik)/sum(k_weights)
         do ib = 1, nbands
            occupation = fermi_dirac_occupation(eigenvalues(ib, ik), fermi_level, temperature)
            if (occupation == 0.0_rp) cycle
            do site = 1, this%nsite
               offset = (site - 1)*2*norb
               local_state = eigenvectors(offset + 1:offset + 2*norb, ib, ik)
               local_operator = cmplx(0.0_rp, 0.0_rp, rp)
               do spin = 1, 2
                  do iorb = 1, norb
                     l = lmto_orbital_l(iorb)
                     if (.not. this%selected_l(l)) cycle
                     local_operator((spin - 1)*norb + iorb, (spin - 1)*norb + iorb) = &
                        cmplx(merge(1.0_rp, -1.0_rp, spin == 1)* &
                              ground_state_radial_norm(ground_states(site), l, spin, eigenvalues(ib, ik)), 0.0_rp, rp)
                  end do
               end do
               action = matmul(local_operator, local_state)
               expectation = real(dot_product(local_state, action), rp)
               moment(site) = moment(site) + wk*occupation*expectation
            end do
         end do
      end do
      deallocate(local_state, local_operator, action)
   end subroutine projected_moment_from_operator

   subroutine projected_core_spin_number(this, ground_states, core_spin)
      class(projected_site_spin_contract), intent(in) :: this
      type(radial_ground_state), intent(in) :: ground_states(:)
      real(rp), intent(out) :: core_spin(:)
      integer :: site

      if (size(ground_states) /= this%nsite .or. size(core_spin) /= this%nsite) then
         error stop 'DRESP-01 core context: site shape mismatch'
      end if
      core_spin = 0.0_rp
      do site = 1, this%nsite
         if (.not. ground_states(site)%core_density_valid .or. &
             .not. allocated(ground_states(site)%core_pauli_weighted_up) .or. &
             .not. allocated(ground_states(site)%core_pauli_weighted_down)) cycle
         core_spin(site) = ground_states(site)%log_mesh_integral( &
            ground_states(site)%core_pauli_weighted_up - ground_states(site)%core_pauli_weighted_down)
      end do
   end subroutine projected_core_spin_number

   subroutine validate_product(this, product, caller)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_product_response_basis), intent(in) :: product
      character(len=*), intent(in) :: caller

      if (.not. allocated(product%blocks) .or. product%nsite /= this%nsite .or. &
          product%orbital_lmax /= this%orbital_lmax .or. product%response_lmax /= this%response_lmax .or. &
          product%npoint /= this%npoint .or. product%product_dimension < 1) then
         error stop trim(caller)//': product representation does not match the DRESP-01 contract'
      end if
      if (product%circular_channel /= lmto_product_channel_plus .and. &
          product%circular_channel /= lmto_product_channel_minus) then
         error stop trim(caller)//': invalid circular product channel'
      end if
   end subroutine validate_product

   subroutine validate_radial_bases(this, radial_bases)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      integer :: site

      if (size(radial_bases) /= this%nsite) then
         error stop 'DRESP-01 radial validation: site shape mismatch'
      end if
      do site = 1, this%nsite
         call radial_bases(site)%require_supported('ham_only', 'second', .true., .true., .false., .false.)
         if (radial_bases(site)%lmax /= this%orbital_lmax .or. radial_bases(site)%nspin /= 2 .or. &
             radial_bases(site)%npoint /= this%npoint .or. .not. allocated(radial_bases(site)%rofi) .or. &
             .not. allocated(radial_bases(site)%phi_large) .or. .not. allocated(radial_bases(site)%phidot_large) .or. &
             .not. allocated(radial_bases(site)%enu_work) .or. .not. allocated(radial_bases(site)%channel_present)) then
            error stop 'DRESP-01 radial validation: incomplete accepted radial basis'
         end if
         if (.not. all(radial_bases(site)%channel_present)) then
            error stop 'DRESP-01 radial validation: incomplete accepted radial basis'
         end if
         if (any(.not. ieee_is_finite(radial_bases(site)%rofi)) .or. &
             any(.not. ieee_is_finite(radial_bases(site)%phi_large)) .or. &
             any(.not. ieee_is_finite(radial_bases(site)%phidot_large)) .or. &
             any(.not. ieee_is_finite(radial_bases(site)%enu_work))) then
            error stop 'DRESP-01 radial validation: non-finite radial data is rejected'
         end if
         if (abs(radial_bases(site)%mesh_a - this%mesh_a) > mesh_tolerance*max(1.0_rp, abs(this%mesh_a)) .or. &
             abs(radial_bases(site)%mesh_b - this%mesh_b) > mesh_tolerance*max(1.0_rp, abs(this%mesh_b)) .or. &
             maxval(abs(radial_bases(site)%rofi - this%radius)) > mesh_tolerance* &
                max(1.0_rp, maxval(abs(this%radius)))) then
            error stop 'DRESP-01 radial validation: mesh provenance mismatch'
         end if
      end do
   end subroutine validate_radial_bases

   subroutine validate_operator_kind(operator_kind, caller)
      integer, intent(in) :: operator_kind
      character(len=*), intent(in) :: caller

      if (operator_kind < projected_operator_plus .or. operator_kind > projected_operator_z) then
         error stop trim(caller)//': operator kind must be plus, minus, or z'
      end if
   end subroutine validate_operator_kind

   subroutine validate_endpoint_shapes(nsite, norb, left_state, right_state, caller)
      integer, intent(in) :: nsite, norb
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      character(len=*), intent(in) :: caller

      if (.not. allocated(left_state%coefficients) .or. .not. allocated(right_state%coefficients) .or. &
          size(left_state%coefficients) /= 2*norb*nsite .or. size(right_state%coefficients) /= 2*norb*nsite) then
         error stop trim(caller)//': endpoint coefficient shape mismatch'
      end if
   end subroutine validate_endpoint_shapes

   subroutine validate_moment_inputs(this, eigenvalues, eigenvectors, k_weights, fermi_level, temperature, &
                                     ground_states, moment)
      class(projected_site_spin_contract), intent(in) :: this
      real(rp), intent(in) :: eigenvalues(:, :), k_weights(:), fermi_level, temperature
      complex(rp), intent(in) :: eigenvectors(:, :, :)
      type(radial_ground_state), intent(in) :: ground_states(:)
      real(rp), intent(out) :: moment(:)
      integer :: site
      integer :: expected_basis

      expected_basis = 2*((this%orbital_lmax + 1)**2)*this%nsite
      if (size(moment) /= this%nsite .or. size(ground_states) /= this%nsite .or. &
          size(eigenvalues, 1) < 1 .or. size(eigenvalues, 2) < 1 .or. size(k_weights) /= size(eigenvalues, 2) .or. &
          any(shape(eigenvectors) /= [expected_basis, size(eigenvalues, 1), size(eigenvalues, 2)]) .or. &
          sum(k_weights) <= tiny(1.0_rp) .or. any(k_weights < 0.0_rp)) then
         error stop 'DRESP-01 moment: accepted eigensystem or k-weight shape is invalid'
      end if
      do site = 1, this%nsite
         if (.not. ground_states(site)%valid .or. .not. ground_states(site)%accepted .or. &
             .not. ground_states(site)%pauli_basis_valid .or. ground_states(site)%pauli_lmax < this%orbital_lmax .or. &
             .not. allocated(ground_states(site)%r) .or. .not. allocated(ground_states(site)%pauli_large) .or. &
             .not. allocated(ground_states(site)%pauli_large_dot) .or. .not. allocated(ground_states(site)%pauli_enu)) then
            error stop 'DRESP-01 moment: accepted Pauli radial snapshot is required'
         end if
         if (size(ground_states(site)%r) /= this%npoint .or. &
             any(shape(ground_states(site)%pauli_large) /= [this%npoint, this%orbital_lmax + 1, 2]) .or. &
             any(shape(ground_states(site)%pauli_large_dot) /= [this%npoint, this%orbital_lmax + 1, 2]) .or. &
             any(shape(ground_states(site)%pauli_enu) /= [this%orbital_lmax + 1, 2])) then
            error stop 'DRESP-01 moment: accepted Pauli radial snapshot dimensions are invalid'
         end if
         if (abs(ground_states(site)%a - this%mesh_a) > mesh_tolerance*max(1.0_rp, abs(this%mesh_a)) .or. &
             abs(ground_states(site)%b - this%mesh_b) > mesh_tolerance*max(1.0_rp, abs(this%mesh_b)) .or. &
             maxval(abs(ground_states(site)%r - this%radius)) > mesh_tolerance*max(1.0_rp, maxval(abs(this%radius)))) then
            error stop 'DRESP-01 moment: radial snapshot mesh provenance mismatch'
         end if
      end do
   end subroutine validate_moment_inputs

   pure real(rp) function endpoint_radial_value(radial, ir, l, spin, energy) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, spin
      real(rp), intent(in) :: energy

      value = radial%phi_large(ir, l + 1, spin) + (energy - radial%enu_work(l + 1, spin))* &
         radial%phidot_large(ir, l + 1, spin)
   end function endpoint_radial_value

   real(rp) function endpoint_radial_overlap(radial, l_left, spin_left, energy_left, l_right, spin_right, energy_right) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: l_left, spin_left, l_right, spin_right
      real(rp), intent(in) :: energy_left, energy_right
      integer :: ir

      value = 0.0_rp
      do ir = 1, radial%npoint
         value = value + radial_simpson_weight(ir, radial%npoint)*radial%mesh_a* &
            (radial%rofi(ir) + radial%mesh_b)*endpoint_radial_value(radial, ir, l_left, spin_left, energy_left)* &
            endpoint_radial_value(radial, ir, l_right, spin_right, energy_right)
      end do
   end function endpoint_radial_overlap

   real(rp) function ground_state_radial_norm(ground_state, l, spin, energy) result(value)
      type(radial_ground_state), intent(in) :: ground_state
      integer, intent(in) :: l, spin
      real(rp), intent(in) :: energy
      integer :: ir
      real(rp) :: u_value

      value = 0.0_rp
      do ir = 1, size(ground_state%r)
         u_value = ground_state%pauli_large(ir, l + 1, spin) + (energy - ground_state%pauli_enu(l + 1, spin))* &
            ground_state%pauli_large_dot(ir, l + 1, spin)
         value = value + radial_simpson_weight(ir, size(ground_state%r))*ground_state%a* &
            (ground_state%r(ir) + ground_state%b)*u_value**2
      end do
   end function ground_state_radial_norm

   pure real(rp) function fermi_dirac_occupation(eigenvalue, fermi_level, temperature) result(occupation)
      real(rp), intent(in) :: eigenvalue, fermi_level, temperature
      real(rp) :: argument

      argument = (eigenvalue - fermi_level)/max(temperature*kB_Ry_per_K, occupation_kT_floor)
      if (argument >= 50.0_rp) then
         occupation = 0.0_rp
      else if (argument <= -50.0_rp) then
         occupation = 1.0_rp
      else
         occupation = 1.0_rp/(exp(argument) + 1.0_rp)
      end if
   end function fermi_dirac_occupation

end module lr_projected_site_spin_mod
