!------------------------------------------------------------------------------
! LR-05 Pauli/no-SOC transition-vector oracles.
!
! The angular oracle is deliberately implemented with independent spherical
! harmonics and direct sphere quadrature.  The production evaluator is only
! used to provide the quantity under test.
!------------------------------------------------------------------------------
program test_lr_pauli_transition_vertex
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use logger_mod, only: g_logger
   use math_mod, only: init_math_operators
   use self_mod, only: legacy_radial_fixture
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use response_angular_basis_mod, only: response_angular_pi, response_lm_index
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex, &
      response_unflatten_superindex, response_apply_site_gauge, response_endpoint_phase
   use lr_response_space_mod, only: response_space_layout
   use lr_pauli_transition_vertex_mod, only: pauli_vertex_capabilities, pauli_endpoint_state, &
      pauli_charge_matrix, pauli_sigma_x_matrix, pauli_sigma_z_matrix, pauli_sigma_plus_matrix, &
      evaluate_pauli_transition_vertex
   implicit none

   integer, parameter :: nr = 51, nsite = 2, lmax = 2, norb = (lmax + 1)**2
   real(rp), parameter :: mesh_a = 0.03_rp, mesh_b = 0.10_rp
   real(rp), parameter :: tolerance = 3.0e-10_rp
   real(rp) :: radius(nr)
   type(lmto_radial_basis) :: radial(nsite)
   type(response_space_layout) :: space, one_site_space
   type(pauli_endpoint_state) :: left_state, right_state, swapped_state
   type(pauli_vertex_capabilities) :: capabilities
   complex(rp), allocatable :: left_coefficients(:), right_coefficients(:)
   complex(rp), allocatable :: transition(:), expected(:), reverse_transition(:)
   complex(rp) :: charge(2, 2), sigma_x(2, 2), sigma_z(2, 2), sigma_plus(2, 2)
   complex(rp), allocatable :: site_coefficients(:, :), gauged_site_coefficients(:, :)
   complex(rp), allocatable :: gauged_coefficients(:), folded_transition(:)
   real(rp) :: tau(3, nsite), reciprocal_shift(3)
   character(len=32) :: mode
   logical :: failed
   integer :: argument_length

   call basis_init(lmax)
   call init_math_operators()
   call g_logger%init()
   call build_mesh(radius)
   call setup_radial_bases(radial, radius, lmax)

   call get_command_argument(1, mode, length=argument_length)
   if (argument_length > 0 .and. len_trim(mode) > 0) then
      call run_capability_rejection(mode, radius)
      error stop 'LR-05 capability rejection test unexpectedly returned'
   end if

   call space%initialize(nsite, 2*lmax, radius, mesh_a, mesh_b, 1)
   allocate(left_coefficients(2*nsite*norb), right_coefficients(2*nsite*norb))
   call fill_coefficients(left_coefficients, 0.13_rp)
   call fill_coefficients(right_coefficients, -0.21_rp)
   call left_state%initialize(-0.37_rp, left_coefficients)
   call right_state%initialize(0.19_rp, right_coefficients)
   capabilities = pauli_vertex_capabilities()
   charge = pauli_charge_matrix()
   sigma_x = pauli_sigma_x_matrix()
   sigma_z = pauli_sigma_z_matrix()
   sigma_plus = pauli_sigma_plus_matrix()

   allocate(transition(space%ndim), expected(space%ndim), reverse_transition(space%ndim), folded_transition(space%ndim))
   failed = .false.
   call angular_pair_oracle('s-s charge', 0, 0, 0, 0, charge, space, radial, capabilities, failed)
   call angular_pair_oracle('p-p off diagonal charge', 1, -1, 1, 1, charge, space, radial, capabilities, failed)
   call angular_pair_oracle('p-d charge', 1, 0, 2, 1, charge, space, radial, capabilities, failed)
   call angular_pair_oracle('d-d off diagonal charge', 2, -2, 2, 2, charge, space, radial, capabilities, failed)
   call angular_pair_oracle('sigma-z', 1, 0, 1, 0, sigma_z, space, radial, capabilities, failed)
   call angular_pair_oracle('transverse sigma-plus', 0, 0, 1, 0, sigma_plus, space, radial, capabilities, failed)

   call evaluate_pauli_transition_vertex(space, radial, left_state, right_state, charge, capabilities, transition)
   call brute_force_vertex(space, radial, left_state, right_state, charge, expected)
   call compare_positive_radial(space, transition, expected, tolerance, 'full direct angular quadrature', failed)

   call evaluate_pauli_transition_vertex(space, radial, left_state, right_state, sigma_z, capabilities, transition)
   call evaluate_pauli_transition_vertex(space, radial, right_state, left_state, sigma_z, capabilities, reverse_transition)
   call hermitian_endpoint_oracle(space, transition, reverse_transition, tolerance, failed)

   call setup_bz_states(left_state, right_state, left_coefficients, right_coefficients, swapped_state, &
                        gauged_coefficients, tau, reciprocal_shift)
   call evaluate_pauli_transition_vertex(space, radial, left_state, right_state, charge, capabilities, transition)
   call evaluate_pauli_transition_vertex(space, radial, left_state, swapped_state, charge, capabilities, folded_transition)
   call folded_endpoint_oracle(space, transition, folded_transition, tau, reciprocal_shift, tolerance, failed)

   call one_site_space%initialize(1, 2*lmax, radius, mesh_a, mesh_b, 1)
   call occupied_diagonal_closure(one_site_space, radial(1:1), capabilities, tolerance, failed)

   if (failed) then
      write (*, '(a)') 'UnitLrPauliTransitionVertex: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrPauliTransitionVertex: PASS (angular, Hermitian, folded endpoint, occupied Pauli closure)'

contains

   subroutine build_mesh(mesh)
      real(rp), intent(out) :: mesh(:)
      integer :: ir

      mesh(1) = 0.0_rp
      do ir = 2, size(mesh)
         mesh(ir) = mesh_b*(exp(mesh_a*real(ir - 1, rp)) - 1.0_rp)
      end do
   end subroutine build_mesh

   subroutine setup_radial_bases(bases, mesh, basis_lmax)
      type(lmto_radial_basis), intent(out) :: bases(:)
      real(rp), intent(in) :: mesh(:)
      integer, intent(in) :: basis_lmax
      real(rp) :: potential(size(mesh)), energy
      real(rp), allocatable :: g(:, :), gp(:, :), gpp(:, :)
      real(rp), allocatable :: gpack(:), gdotpack(:), gddotpack(:)
      integer :: isite, ispin, l

      potential = 0.0_rp
      do isite = 1, size(bases)
         call bases(isite)%initialize(size(mesh), basis_lmax, 2)
         do ispin = 1, 2
            do l = 0, basis_lmax
               call legacy_radial_fixture(1.0_rp, l, mesh_a, mesh_b, mesh, potential, energy, g, gp, gpp, &
                                          0.31_rp + 0.07_rp*real(ispin + l, rp))
               gpack = reshape(g, [2*size(mesh)])
               gdotpack = reshape(gp, [2*size(mesh)])
               gddotpack = reshape(gpp, [2*size(mesh)])
               call bases(isite)%capture_channel(l, ispin, energy, mesh, potential, mesh_a, mesh_b, 1.0_rp, &
                                                gpack, gdotpack, gddotpack, energy)
            end do
         end do
      end do
   end subroutine setup_radial_bases

   subroutine fill_coefficients(coefficients, offset)
      complex(rp), intent(out) :: coefficients(:)
      real(rp), intent(in) :: offset
      integer :: i

      do i = 1, size(coefficients)
         coefficients(i) = cmplx(offset + 0.007_rp*real(i, rp), &
                                 -0.011_rp + 0.003_rp*real(mod(i, 7), rp), rp)
      end do
   end subroutine fill_coefficients

   subroutine angular_pair_oracle(label, l, m, lp, mp, operator_matrix, space, radial_bases, capabilities, failed)
      character(len=*), intent(in) :: label
      integer, intent(in) :: l, m, lp, mp
      complex(rp), intent(in) :: operator_matrix(2, 2)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(pauli_vertex_capabilities), intent(in) :: capabilities
      logical, intent(inout) :: failed
      type(pauli_endpoint_state) :: left, right
      complex(rp), allocatable :: left_coeff(:), right_coeff(:), actual(:), brute(:)
      integer :: iorb_left, iorb_right, ncoeff

      ncoeff = 2*space%nsite*(radial_bases(1)%lmax + 1)**2
      allocate(left_coeff(ncoeff), right_coeff(ncoeff), actual(space%ndim), brute(space%ndim))
      left_coeff = cmplx(0.0_rp, 0.0_rp, rp)
      right_coeff = cmplx(0.0_rp, 0.0_rp, rp)
      iorb_left = response_lm_index(l, m)
      iorb_right = response_lm_index(lp, mp)
      left_coeff(iorb_left) = cmplx(0.43_rp, -0.17_rp, rp)
      right_coeff(iorb_right) = cmplx(-0.29_rp, 0.38_rp, rp)
      left_coeff(2*iorb_left) = cmplx(0.21_rp, 0.09_rp, rp)
      right_coeff(2*iorb_right) = cmplx(-0.11_rp, -0.27_rp, rp)
      call left%initialize(-0.23_rp, left_coeff)
      call right%initialize(0.41_rp, right_coeff)
      call evaluate_pauli_transition_vertex(space, radial_bases, left, right, operator_matrix, capabilities, actual)
      call brute_force_vertex(space, radial_bases, left, right, operator_matrix, brute)
      call compare_positive_radial(space, actual, brute, tolerance, trim(label), failed)
   end subroutine angular_pair_oracle

   subroutine brute_force_vertex(space, radial_bases, left, right, operator_matrix, result)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(pauli_endpoint_state), intent(in) :: left, right
      complex(rp), intent(in) :: operator_matrix(2, 2)
      complex(rp), intent(out) :: result(:)
      integer, parameter :: ntheta = 12, nphi = 24
      real(rp) :: nodes(ntheta), weights(ntheta), theta, phi, sphere_weight, radius_value
      complex(rp) :: density, left_u, right_u, yi, yj, ylm
      type(response_super_index) :: item
      integer :: flat, isite, ir, response_l, response_m, iorb, jorb, ispin, jspin, orbital_l, orbital_lp, iphi, inode
      integer :: norb_local, offset

      if (size(result) /= space%ndim) error stop 'brute_force_vertex: result shape mismatch'
      call gauss_legendre(ntheta, nodes, weights)
      result = cmplx(0.0_rp, 0.0_rp, rp)
      norb_local = (radial_bases(1)%lmax + 1)**2
      do isite = 1, space%nsite
         offset = (isite - 1)*2*norb_local
         do ir = 2, space%npoint
            radius_value = space%radius(ir)
            do inode = 1, ntheta
               theta = acos(nodes(inode))
               do iphi = 1, nphi
                  phi = 2.0_rp*response_angular_pi*real(iphi - 1, rp)/real(nphi, rp)
                  density = cmplx(0.0_rp, 0.0_rp, rp)
                  do iorb = 1, norb_local
                     orbital_l = orbital_l_from_index(iorb)
                     yi = independent_harmonic(orbital_l, orbital_m_from_index(iorb, orbital_l), theta, phi)
                     do jorb = 1, norb_local
                        orbital_lp = orbital_l_from_index(jorb)
                        yj = independent_harmonic(orbital_lp, orbital_m_from_index(jorb, orbital_lp), theta, phi)
                        do ispin = 1, 2
                           left_u = cmplx(radial_bases(isite)%phi_large(ir, orbital_l + 1, ispin) + &
                              (left%energy - radial_bases(isite)%enu_work(orbital_l + 1, ispin))* &
                              radial_bases(isite)%phidot_large(ir, orbital_l + 1, ispin), 0.0_rp, rp)
                           do jspin = 1, 2
                              right_u = cmplx(radial_bases(isite)%phi_large(ir, orbital_lp + 1, jspin) + &
                                 (right%energy - radial_bases(isite)%enu_work(orbital_lp + 1, jspin))* &
                                 radial_bases(isite)%phidot_large(ir, orbital_lp + 1, jspin), 0.0_rp, rp)
                              density = density + conjg(left%coefficients(offset + (ispin - 1)*norb_local + iorb))* &
                                 right%coefficients(offset + (jspin - 1)*norb_local + jorb)* &
                                 operator_matrix(ispin, jspin)*left_u*right_u*conjg(yi)*yj/(radius_value**2)
                           end do
                        end do
                     end do
                  end do
                  sphere_weight = weights(inode)*(2.0_rp*response_angular_pi/real(nphi, rp))
                  do response_l = 0, space%response_lmax
                     do response_m = -response_l, response_l
                        item = response_super_index(isite, response_l, response_m, ir, 1)
                        call response_flatten_superindex(item, space%nsite, space%response_lmax, space%npoint, &
                                                         space%nchannel, flat)
                        ylm = independent_harmonic(response_l, response_m, theta, phi)
                        result(flat) = result(flat) + sphere_weight*conjg(ylm)*density
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine brute_force_vertex

   subroutine compare_positive_radial(space, actual, expected, tolerance_in, label, failed)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: actual(:), expected(:)
      real(rp), intent(in) :: tolerance_in
      character(len=*), intent(in) :: label
      logical, intent(inout) :: failed
      type(response_super_index) :: item
      real(rp) :: error
      integer :: flat

      error = 0.0_rp
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, &
                                            space%nchannel, item)
         if (item%radial_point > 1) error = max(error, abs(actual(flat) - expected(flat)))
      end do
      write (*, '(a,es12.4)') trim(label)//' error = ', error
      if (error > tolerance_in) failed = .true.
   end subroutine compare_positive_radial

   subroutine hermitian_endpoint_oracle(space, forward, reverse, tolerance_in, failed)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: forward(:), reverse(:)
      real(rp), intent(in) :: tolerance_in
      logical, intent(inout) :: failed
      type(response_super_index) :: item, partner
      complex(rp) :: expected
      real(rp) :: error
      integer :: flat, partner_flat

      error = 0.0_rp
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, &
                                            space%nchannel, item)
         partner = response_super_index(item%site, item%response_l, -item%response_m, item%radial_point, item%channel)
         call response_flatten_superindex(partner, space%nsite, space%response_lmax, space%npoint, &
                                          space%nchannel, partner_flat)
         expected = cmplx(merge(1.0_rp, -1.0_rp, mod(abs(item%response_m), 2) == 0), 0.0_rp, rp)* &
                    conjg(forward(partner_flat))
         error = max(error, abs(reverse(flat) - expected))
      end do
      write (*, '(a,es12.4)') 'Hermitian endpoint relation error = ', error
      if (error > tolerance_in) failed = .true.
   end subroutine hermitian_endpoint_oracle

   subroutine setup_bz_states(left, right, left_coefficients, right_coefficients, gauged_state, gauged_coefficients, tau, G)
      type(pauli_endpoint_state), intent(in) :: left, right
      type(pauli_endpoint_state), intent(out) :: gauged_state
      complex(rp), intent(in) :: left_coefficients(:), right_coefficients(:)
      complex(rp), allocatable, intent(out) :: gauged_coefficients(:)
      real(rp), intent(out) :: tau(:, :), G(3)
      complex(rp), allocatable :: site_coefficients(:, :), transformed(:, :)
      integer :: site, block_size, offset

      tau(:, 1) = [0.0_rp, 0.0_rp, 0.0_rp]
      tau(:, 2) = [0.25_rp, 0.125_rp, 0.0_rp]
      G = [1.0_rp, 0.0_rp, 0.0_rp]
      block_size = 2*norb
      allocate(site_coefficients(block_size, nsite), transformed(block_size, nsite), gauged_coefficients(size(right_coefficients)))
      do site = 1, nsite
         offset = (site - 1)*block_size
         site_coefficients(:, site) = right_coefficients(offset + 1:offset + block_size)
      end do
      call response_apply_site_gauge(site_coefficients, tau, G, transformed)
      do site = 1, nsite
         offset = (site - 1)*block_size
         gauged_coefficients(offset + 1:offset + block_size) = transformed(:, site)
      end do
      call gauged_state%initialize(right%energy, gauged_coefficients)
   end subroutine setup_bz_states

   subroutine folded_endpoint_oracle(space, raw, folded, tau, G, tolerance_in, failed)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: raw(:), folded(:)
      real(rp), intent(in) :: tau(:, :), G(3), tolerance_in
      logical, intent(inout) :: failed
      type(response_super_index) :: item
      complex(rp) :: expected
      real(rp) :: error
      integer :: flat

      error = 0.0_rp
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, &
                                            space%nchannel, item)
         expected = response_endpoint_phase(tau(:, item%site), G)*raw(flat)
         error = max(error, abs(folded(flat) - expected))
      end do
      write (*, '(a,es12.4)') 'Folded k+q complete transition-vector error = ', error
      if (error > tolerance_in) failed = .true.
   end subroutine folded_endpoint_oracle

   subroutine occupied_diagonal_closure(space, radial_bases, capabilities, tolerance_in, failed)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(pauli_vertex_capabilities), intent(in) :: capabilities
      real(rp), intent(in) :: tolerance_in
      logical, intent(inout) :: failed
      type(pauli_endpoint_state) :: state
      complex(rp), allocatable :: coefficients(:), charge_sum(:), sigma_z_sum(:), charge(:), sigma_z(:)
      complex(rp), allocatable :: expected_charge_vector(:), expected_sigma_z_vector(:)
      complex(rp) :: charge_matrix(2, 2), sigma_z_matrix(2, 2)
      type(response_super_index) :: item
      real(rp) :: expected_charge, expected_sigma_z, radial_product, first_product, second_product
      real(rp) :: state_energy
      integer :: nocc, ib, orbital, ispin, ir, flat, l

      nocc = 2*norb
      allocate(coefficients(2*norb), charge_sum(space%ndim), sigma_z_sum(space%ndim), &
               charge(space%ndim), sigma_z(space%ndim), expected_charge_vector(space%ndim), &
               expected_sigma_z_vector(space%ndim))
      charge_sum = cmplx(0.0_rp, 0.0_rp, rp)
      sigma_z_sum = cmplx(0.0_rp, 0.0_rp, rp)
      expected_charge_vector = cmplx(0.0_rp, 0.0_rp, rp)
      expected_sigma_z_vector = cmplx(0.0_rp, 0.0_rp, rp)
      charge_matrix = pauli_charge_matrix()
      sigma_z_matrix = pauli_sigma_z_matrix()
      do ib = 1, nocc
         coefficients = cmplx(0.0_rp, 0.0_rp, rp)
         coefficients(ib) = cmplx(1.0_rp, 0.0_rp, rp)
         state_energy = -0.41_rp + 0.013_rp*real(ib, rp)
         call state%initialize(state_energy, coefficients)
         call evaluate_pauli_transition_vertex(space, radial_bases, state, state, charge_matrix, capabilities, charge)
         call evaluate_pauli_transition_vertex(space, radial_bases, state, state, sigma_z_matrix, capabilities, sigma_z)
         charge_sum = charge_sum + charge
         sigma_z_sum = sigma_z_sum + sigma_z
         orbital = mod(ib - 1, norb) + 1
         ispin = (ib - 1)/norb + 1
         l = orbital_l_from_index(orbital)
         do flat = 1, space%ndim
            call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, &
                                               space%nchannel, item)
            if (item%response_l /= 0 .or. item%response_m /= 0) cycle
            ir = item%radial_point
            if (ir == 1) then
               if (l /= 0) cycle
               first_product = radial_amplitude(radial_bases(1), state_energy, l, ispin, 2)**2/(space%radius(2)**2)
               second_product = radial_amplitude(radial_bases(1), state_energy, l, ispin, 3)**2/(space%radius(3)**2)
               radial_product = (first_product*space%radius(3)**2 - second_product*space%radius(2)**2)/ &
                                (space%radius(3)**2 - space%radius(2)**2)
            else
               radial_product = radial_amplitude(radial_bases(1), state_energy, l, ispin, ir)**2/(space%radius(ir)**2)
            end if
            expected_charge_vector(flat) = expected_charge_vector(flat) + &
               cmplx(radial_product/sqrt(4.0_rp*response_angular_pi), 0.0_rp, rp)
            expected_sigma_z_vector(flat) = expected_sigma_z_vector(flat) + &
               cmplx(merge(1.0_rp, -1.0_rp, ispin == 1)*radial_product/sqrt(4.0_rp*response_angular_pi), &
                     0.0_rp, rp)
         end do
      end do

      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, &
                                            space%nchannel, item)
         expected_charge = real(expected_charge_vector(flat), rp)
         expected_sigma_z = real(expected_sigma_z_vector(flat), rp)
         if (abs(charge_sum(flat) - expected_charge_vector(flat)) > tolerance_in .or. &
             abs(sigma_z_sum(flat) - expected_sigma_z_vector(flat)) > tolerance_in) failed = .true.
      end do
      write (*, '(a)') 'Occupied diagonal charge/sigma-z Pauli closure: checked against LMTO large-component arrays'
   end subroutine occupied_diagonal_closure

   real(rp) function radial_amplitude(basis, energy, l, ispin, ir) result(value)
      type(lmto_radial_basis), intent(in) :: basis
      real(rp), intent(in) :: energy
      integer, intent(in) :: l, ispin, ir

      value = basis%phi_large(ir, l + 1, ispin) + (energy - basis%enu_work(l + 1, ispin))* &
              basis%phidot_large(ir, l + 1, ispin)
   end function radial_amplitude

   subroutine run_capability_rejection(argument, mesh)
      character(len=*), intent(in) :: argument
      real(rp), intent(in) :: mesh(:)
      type(lmto_radial_basis) :: guarded_radial(1)
      type(response_space_layout) :: guarded_space
      type(pauli_endpoint_state) :: state
      type(pauli_vertex_capabilities) :: guarded_capabilities
      complex(rp) :: operator_matrix(2, 2)
      complex(rp) :: coefficients(2*(3 + 1)**2)
      complex(rp), allocatable :: result(:)
      integer :: guarded_lmax

      guarded_lmax = lmax
      if (trim(argument) == 'spdf') guarded_lmax = 3
      call setup_radial_bases(guarded_radial, mesh, guarded_lmax)
      call guarded_space%initialize(1, 2*guarded_lmax, mesh, mesh_a, mesh_b, 1)
      coefficients = cmplx(0.0_rp, 0.0_rp, rp)
      call state%initialize(0.0_rp, coefficients)
      guarded_capabilities = pauli_vertex_capabilities()
      select case (trim(argument))
      case ('soc'); guarded_capabilities%has_soc = .true.
      case ('overlap'); guarded_capabilities%reciprocal_mode = 'generalized_overlap_proxy'
      case ('noncollinear'); guarded_capabilities%collinear = .false.
      case ('additive'); guarded_capabilities%has_extra_operator = .true.
      case ('spdf')
      case default
         error stop 'unknown LR-05 capability rejection selector'
      end select
      operator_matrix = pauli_charge_matrix()
      allocate(result(guarded_space%ndim))
      call evaluate_pauli_transition_vertex(guarded_space, guarded_radial, state, state, operator_matrix, &
                                            guarded_capabilities, result)
   end subroutine run_capability_rejection

   pure integer function orbital_l_from_index(index) result(l)
      integer, intent(in) :: index
      integer :: trial

      l = -1
      do trial = 0, 8
         if (trial*trial + 1 <= index .and. index <= (trial + 1)*(trial + 1)) then
            l = trial
            return
         end if
      end do
   end function orbital_l_from_index

   pure integer function orbital_m_from_index(index, l) result(m)
      integer, intent(in) :: index, l
      m = index - l*l - l - 1
   end function orbital_m_from_index

   pure function independent_harmonic(l, m, theta, phi) result(value)
      integer, intent(in) :: l, m
      real(rp), intent(in) :: theta, phi
      complex(rp) :: value, standard
      real(rp) :: p, norm
      integer :: ma

      ma = abs(m)
      p = independent_legendre(l, ma, cos(theta))
      norm = sqrt(real(2*l + 1, rp)/(4.0_rp*response_angular_pi)* &
                  independent_factorial(l - ma)/independent_factorial(l + ma))
      if (m >= 0) then
         standard = cmplx(norm*independent_sign(ma)*p*cos(real(m, rp)*phi), &
                          norm*independent_sign(ma)*p*sin(real(m, rp)*phi), rp)
      else
         standard = independent_sign(ma)*conjg(cmplx(norm*independent_sign(ma)*p*cos(real(ma, rp)*phi), &
                                                    norm*independent_sign(ma)*p*sin(real(ma, rp)*phi), rp))
      end if
      value = conjg(standard)
   end function independent_harmonic

   pure real(rp) function independent_factorial(n) result(value)
      integer, intent(in) :: n
      integer :: i

      value = 1.0_rp
      do i = 2, n
         value = value*real(i, rp)
      end do
   end function independent_factorial

   pure real(rp) function independent_sign(n) result(value)
      integer, intent(in) :: n
      value = merge(1.0_rp, -1.0_rp, mod(abs(n), 2) == 0)
   end function independent_sign

   pure real(rp) function independent_legendre(l, m, x) result(value)
      integer, intent(in) :: l, m
      real(rp), intent(in) :: x
      real(rp) :: pmm, pmmp1, pll, somx2
      integer :: ell

      pmm = 1.0_rp
      if (m > 0) then
         somx2 = sqrt(max(0.0_rp, (1.0_rp - x)*(1.0_rp + x)))
         do ell = 1, m
            pmm = pmm*real(2*ell - 1, rp)*somx2
         end do
      end if
      if (l == m) then
         value = pmm
      else
         pmmp1 = x*real(2*m + 1, rp)*pmm
         if (l == m + 1) then
            value = pmmp1
         else
            do ell = m + 2, l
               pll = (real(2*ell - 1, rp)*x*pmmp1 - real(ell + m - 1, rp)*pmm)/real(ell - m, rp)
               pmm = pmmp1
               pmmp1 = pll
            end do
            value = pmmp1
         end if
      end if
   end function independent_legendre

   subroutine gauss_legendre(n, nodes, weights)
      integer, intent(in) :: n
      real(rp), intent(out) :: nodes(n), weights(n)
      integer :: i, j, half
      real(rp) :: z, z1, p1, p2, p3, pp

      half = (n + 1)/2
      do i = 1, half
         z = cos(response_angular_pi*(real(i, rp) - 0.25_rp)/(real(n, rp) + 0.5_rp))
         do
            p1 = 1.0_rp
            p2 = 0.0_rp
            do j = 1, n
               p3 = p2
               p2 = p1
               p1 = ((2.0_rp*real(j, rp) - 1.0_rp)*z*p2 - real(j - 1, rp)*p3)/real(j, rp)
            end do
            pp = real(n, rp)*(z*p1 - p2)/(z*z - 1.0_rp)
            z1 = z
            z = z1 - p1/pp
            if (abs(z - z1) < 2.0e-15_rp) exit
         end do
         nodes(i) = -z
         nodes(n + 1 - i) = z
         weights(i) = 2.0_rp/((1.0_rp - z*z)*pp*pp)
         weights(n + 1 - i) = weights(i)
      end do
   end subroutine gauss_legendre

end program test_lr_pauli_transition_vertex
