!------------------------------------------------------------------------------
! DRESP-01 projected site-spin operator and matching moment contract.
!------------------------------------------------------------------------------
program test_lr_projected_site_spin
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use logger_mod, only: g_logger
   use math_mod, only: init_math_operators
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use radial_ground_state_mod, only: radial_ground_state
   use response_angular_basis_mod, only: response_angular_pi
   use lr_response_space_mod, only: response_space_layout
   use lr_pauli_transition_vertex_mod, only: pauli_endpoint_state
   use lr_lmto_product_response_basis_mod, only: lmto_product_channel_plus, lmto_product_channel_minus, &
      lmto_product_response_basis
   use lr_projected_site_spin_mod, only: projected_site_spin_contract, projected_selector_d, &
      projected_selector_spd, projected_operator_plus, projected_operator_minus, projected_operator_z
   implicit none

   integer, parameter :: nr = 51, lmax = 2, nsite = 2
   integer, parameter :: norb = (lmax + 1)**2, nbasis = 2*norb*nsite
   integer, parameter :: nbands = nbasis, nk = 2
   real(rp), parameter :: mesh_a = 0.03_rp, mesh_b = 0.10_rp, nuclear_z = 1.0_rp
   real(rp), parameter :: tolerance = 2.0e-10_rp

   real(rp) :: radius(nr), eigenvalues(nbands, nk), k_weights(nk)
   complex(rp) :: eigenvectors(nbasis, nbands, nk)
   type(lmto_radial_basis) :: radial(nsite)
   type(radial_ground_state) :: ground_states(nsite)
   type(response_space_layout) :: space
   type(lmto_product_response_basis) :: product_plus, product_minus
   type(projected_site_spin_contract) :: contract_d, contract_spd
   type(pauli_endpoint_state) :: left_state, right_state, right_q_state
   complex(rp), allocatable :: left_coefficients(:), right_coefficients(:), right_q_coefficients(:)
   complex(rp), allocatable :: amplitudes(:), expected(:), q_amplitudes(:), q_expected(:)
   complex(rp), allocatable :: plus_q0_amplitudes(:)
   complex(rp), allocatable :: functionals(:, :), operator_plus(:, :), operator_minus(:, :), operator_z(:, :)
   real(rp), allocatable :: density_moment(:), operator_moment(:), d_moment(:), core_context(:)
   logical :: failed
   character(len=32) :: mode

   call basis_init(lmax)
   call init_math_operators()
   call g_logger%init()
   call build_mesh(radius)
   call setup_radial_bases(radial, radius)

   mode = ''
   call get_command_argument(1, mode)
   if (trim(mode) == 'spdf') then
      call space%initialize(nsite, 4, radius, mesh_a, mesh_b, 1)
      call contract_d%initialize(space, radial, 'spdf')
      error stop 'UnitLrProjectedSiteSpin: spdf rejection unexpectedly returned'
   end if
   if (len_trim(mode) > 0) error stop 'UnitLrProjectedSiteSpin: unknown argument'

   call space%initialize(nsite, 4, radius, mesh_a, mesh_b, 1)
   call product_plus%initialize(space, radial, lmto_product_channel_plus, .false.)
   call product_minus%initialize(space, radial, lmto_product_channel_minus, .false.)
   call contract_d%initialize(space, radial, 'd')
   call contract_spd%initialize(space, radial, 'spd')
   call build_endpoints(left_coefficients, right_coefficients, left_state, right_state)
   call build_finite_q_endpoint(right_coefficients, right_q_coefficients, right_q_state)

   failed = .false.
   if (contract_d%selector_kind /= projected_selector_d .or. contract_spd%selector_kind /= projected_selector_spd) then
      failed = .true.
   end if
   allocate(amplitudes(nsite), expected(nsite), q_amplitudes(nsite), q_expected(nsite), plus_q0_amplitudes(nsite))
   call compare_transition(contract_d, product_plus, radial, left_state, right_state, &
      projected_operator_plus, amplitudes, expected, tolerance, 'd plus q=0', failed)
   call compare_transition(contract_spd, product_plus, radial, left_state, right_state, &
      projected_operator_plus, amplitudes, expected, tolerance, 'spd plus q=0', failed)
   plus_q0_amplitudes = amplitudes
   call compare_transition(contract_d, product_minus, radial, left_state, right_state, &
      projected_operator_minus, amplitudes, expected, tolerance, 'd minus q=0', failed)
   call compare_transition(contract_spd, product_minus, radial, left_state, right_state, &
      projected_operator_minus, amplitudes, expected, tolerance, 'spd minus q=0', failed)
   call compare_transition(contract_d, product_plus, radial, left_state, right_q_state, &
      projected_operator_plus, q_amplitudes, q_expected, tolerance, 'd plus finite-q endpoint', failed)
   call compare_transition(contract_spd, product_plus, radial, left_state, right_q_state, &
      projected_operator_plus, q_amplitudes, q_expected, tolerance, 'spd plus finite-q endpoint', failed)
   if (abs(q_amplitudes(1) - plus_q0_amplitudes(1)) > tolerance .or. &
       abs(q_amplitudes(2) - plus_q0_amplitudes(2)) <= 1.0e-12_rp) failed = .true.

   call check_functional(contract_d, product_plus, failed)
   call check_direct_operator_matrices(contract_d, radial, failed)

   call setup_ground_states(ground_states, radial, radius)
   call build_moment_eigensystem(eigenvalues, eigenvectors, k_weights)
   allocate(density_moment(nsite), operator_moment(nsite), d_moment(nsite), core_context(nsite))
   call contract_spd%moment_from_density(eigenvalues, eigenvectors, k_weights, 0.0_rp, 0.0_rp, &
      ground_states, density_moment)
   call contract_spd%moment_from_operator(eigenvalues, eigenvectors, k_weights, 0.0_rp, 0.0_rp, &
      ground_states, operator_moment)
   call contract_d%moment_from_operator(eigenvalues, eigenvectors, k_weights, 0.0_rp, 0.0_rp, &
      ground_states, d_moment)
   call contract_spd%core_spin_number(ground_states, core_context)
   write (*, '(a,es12.4)') '  maximum density/operator moment error = ', maxval(abs(density_moment - operator_moment))
   write (*, '(a,es12.4)') '  d/spd moment separation = ', maxval(abs(operator_moment - d_moment))
   write (*, '(a,es12.4)') '  core context norm = ', maxval(abs(core_context))
   if (maxval(abs(density_moment - operator_moment)) > tolerance .or. &
       maxval(abs(operator_moment)) <= maxval(abs(d_moment))) failed = .true.

   if (failed) then
      write (*, '(a)') 'UnitLrProjectedSiteSpin: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrProjectedSiteSpin: PASS (selectors, L=0 site functional, q endpoints, operator oracle, moments)'

contains

   subroutine build_mesh(mesh)
      real(rp), intent(out) :: mesh(:)
      integer :: ir

      mesh(1) = 0.0_rp
      do ir = 2, size(mesh)
         mesh(ir) = mesh_b*(exp(mesh_a*real(ir - 1, rp)) - 1.0_rp)
      end do
   end subroutine build_mesh

   subroutine setup_radial_bases(bases, mesh)
      type(lmto_radial_basis), intent(out) :: bases(:)
      real(rp), intent(in) :: mesh(:)
      integer :: isite, ispin, l, ir
      real(rp) :: scale, base

      do isite = 1, size(bases)
         call bases(isite)%initialize(size(mesh), lmax, 2)
         bases(isite)%rofi = mesh
         bases(isite)%mesh_a = mesh_a
         bases(isite)%mesh_b = mesh_b
         bases(isite)%nuclear_z = nuclear_z
         do ispin = 1, 2
            do l = 0, lmax
               bases(isite)%enu_radial(l + 1, ispin) = -0.35_rp + 0.04_rp*real(l + ispin, rp)
               bases(isite)%enu_work(l + 1, ispin) = bases(isite)%enu_radial(l + 1, ispin)
               scale = 0.8_rp + 0.05_rp*real(isite + ispin, rp) + 0.03_rp*real(l, rp)
               do ir = 1, size(mesh)
                  base = mesh(ir)**real(l + 1, rp)*exp(-scale*mesh(ir))
                  bases(isite)%phi_large(ir, l + 1, ispin) = base
                  bases(isite)%phidot_large(ir, l + 1, ispin) = (0.12_rp + 0.01_rp*real(l, rp))*base + &
                     0.02_rp*mesh(ir)**real(l + 2, rp)*exp(-0.5_rp*scale*mesh(ir))
                  bases(isite)%phi_small(ir, l + 1, ispin) = 0.0_rp
                  bases(isite)%phidot_small(ir, l + 1, ispin) = 0.0_rp
                  bases(isite)%phiddot_large(ir, l + 1, ispin) = 0.0_rp
                  bases(isite)%phiddot_small(ir, l + 1, ispin) = 0.0_rp
                  bases(isite)%gfac(ir, l + 1, ispin) = 1.0_rp
               end do
               bases(isite)%channel_present(l + 1, ispin) = .true.
            end do
         end do
      end do
   end subroutine setup_radial_bases

   subroutine build_endpoints(left_coefficients, right_coefficients, left, right)
      complex(rp), allocatable, intent(out) :: left_coefficients(:), right_coefficients(:)
      type(pauli_endpoint_state), intent(out) :: left, right
      integer :: i

      allocate(left_coefficients(nbasis), right_coefficients(nbasis))
      do i = 1, nbasis
         left_coefficients(i) = cmplx(0.12_rp + 0.013_rp*real(i, rp), &
            -0.09_rp + 0.007_rp*real(mod(i, 5), rp), rp)
         right_coefficients(i) = cmplx(-0.17_rp + 0.009_rp*real(i, rp), &
            0.08_rp - 0.005_rp*real(mod(i, 7), rp), rp)
      end do
      call left%initialize(-0.23_rp, left_coefficients)
      call right%initialize(0.41_rp, right_coefficients)
   end subroutine build_endpoints

   subroutine build_finite_q_endpoint(right_coefficients, q_coefficients, q_state)
      complex(rp), intent(in) :: right_coefficients(:)
      complex(rp), allocatable, intent(out) :: q_coefficients(:)
      type(pauli_endpoint_state), intent(out) :: q_state
      integer :: site, offset
      real(rp) :: angle
      complex(rp) :: phase

      allocate(q_coefficients(size(right_coefficients)))
      q_coefficients = right_coefficients
      do site = 1, nsite
         ! A mesh-compatible nonzero q is represented here only through the
         ! already prepared endpoint gauge; the operator adds no tau phase.
         angle = 2.0_rp*response_angular_pi*0.25_rp*0.5_rp*real(site - 1, rp)
         phase = cmplx(cos(angle), sin(angle), rp)
         offset = (site - 1)*2*norb
         q_coefficients(offset + 1:offset + 2*norb) = phase*q_coefficients(offset + 1:offset + 2*norb)
      end do
      call q_state%initialize(0.41_rp, q_coefficients)
   end subroutine build_finite_q_endpoint

   subroutine compare_transition(contract, product, bases, left, right, operator_kind, actual, expected, tol, label, failed)
      type(projected_site_spin_contract), intent(in) :: contract
      type(lmto_product_response_basis), intent(in) :: product
      type(lmto_radial_basis), intent(in) :: bases(:)
      type(pauli_endpoint_state), intent(in) :: left, right
      integer, intent(in) :: operator_kind
      complex(rp), intent(out) :: actual(:), expected(:)
      real(rp), intent(in) :: tol
      character(len=*), intent(in) :: label
      logical, intent(inout) :: failed
      real(rp) :: error

      call contract%transition_amplitudes(product, left, right, actual)
      call contract%direct_transition_amplitudes(bases, left, right, operator_kind, expected)
      error = maxval(abs(actual - expected))
      write (*, '(a,es12.4)') trim(label)//' error = ', error
      if (error > tol) failed = .true.
   end subroutine compare_transition

   subroutine check_functional(contract, product, failed)
      type(projected_site_spin_contract), intent(in) :: contract
      type(lmto_product_response_basis), intent(in) :: product
      logical, intent(inout) :: failed
      complex(rp), allocatable :: functionals(:, :)
      integer :: flat, site, response_l, response_m, product_mode, other_site
      real(rp) :: error

      allocate(functionals(product%product_dimension, nsite))
      call contract%site_integration_functional(product, functionals)
      error = 0.0_rp
      do flat = 1, product%product_dimension
         call product%unflatten_index(flat, site, response_l, response_m, product_mode)
         if (response_l /= 0 .or. response_m /= 0) error = max(error, maxval(abs(functionals(flat, :))))
         do other_site = 1, nsite
            if (other_site /= site) error = max(error, abs(functionals(flat, other_site)))
         end do
      end do
      write (*, '(a,es12.4)') '  site functional non-L0/cross-site error = ', error
      if (error > tolerance) failed = .true.
      deallocate(functionals)
   end subroutine check_functional

   subroutine check_direct_operator_matrices(contract, bases, failed)
      type(projected_site_spin_contract), intent(in) :: contract
      type(lmto_radial_basis), intent(in) :: bases(:)
      logical, intent(inout) :: failed
      complex(rp), allocatable :: plus(:, :), minus(:, :), z(:, :)
      real(rp) :: error
      integer :: i, j

      allocate(plus(nbasis, nbasis), minus(nbasis, nbasis), z(nbasis, nbasis))
      call contract%direct_operator_matrix(bases, 0.17_rp, projected_operator_plus, plus)
      call contract%direct_operator_matrix(bases, 0.17_rp, projected_operator_minus, minus)
      call contract%direct_operator_matrix(bases, 0.17_rp, projected_operator_z, z)
      error = max(maxval(abs(minus - conjg(transpose(plus)))), maxval(abs(z - conjg(transpose(z)))))
      do i = 1, nbasis
         do j = 1, nbasis
            if (lmto_orbital_l(modulo(modulo(i - 1, 2*norb), norb) + 1) /= 2 .and. abs(plus(i, j)) > tolerance) error = huge(1.0_rp)
            if (lmto_orbital_l(modulo(modulo(i - 1, 2*norb), norb) + 1) /= 2 .and. abs(z(i, j)) > tolerance) error = huge(1.0_rp)
         end do
      end do
      write (*, '(a,es12.4)') '  direct operator Hermitian/selector error = ', error
      if (error > tolerance) failed = .true.
      deallocate(plus, minus, z)
   end subroutine check_direct_operator_matrices

   subroutine setup_ground_states(states, bases, mesh)
      type(radial_ground_state), intent(out) :: states(:)
      type(lmto_radial_basis), intent(in) :: bases(:)
      real(rp), intent(in) :: mesh(:)
      real(rp) :: rho(nr, 2), vxc(nr), total_v(nr, 2), core(nr, 2), core_pauli(nr, 2)
      integer :: site, l, spin

      rho = 0.0_rp
      vxc = 0.0_rp
      total_v = 0.0_rp
      core = 0.0_rp
      core_pauli = 0.0_rp
      do site = 1, size(states)
         call states(site)%capture(mesh_a, mesh_b, mesh, rho, [0.0_rp, 0.0_rp], vxc, vxc, total_v, 0.0_rp)
         call states(site)%begin_pauli_basis(lmax)
         do spin = 1, 2
            do l = 0, lmax
               call states(site)%set_pauli_basis_channel(l, spin, bases(site)%enu_work(l + 1, spin), &
                  bases(site)%phi_large(:, l + 1, spin), bases(site)%phidot_large(:, l + 1, spin), &
                  bases(site)%phi_small(:, l + 1, spin), bases(site)%phiddot_large(:, l + 1, spin)*0.0_rp)
            end do
         end do
         call states(site)%capture_core_density(core, core_pauli)
         call states(site)%mark_accepted(1, 1.0e-12_rp)
      end do
   end subroutine setup_ground_states

   subroutine build_moment_eigensystem(energies, vectors, weights)
      real(rp), intent(out) :: energies(:, :), weights(:)
      complex(rp), intent(out) :: vectors(:, :, :)
      integer :: ik, ib, local_index

      energies = 0.0_rp
      vectors = cmplx(0.0_rp, 0.0_rp, rp)
      weights = [0.4_rp, 0.6_rp]
      do ik = 1, nk
         do ib = 1, nbands
            local_index = modulo(ib - 1, 2*norb) + 1
            if (local_index <= norb) then
               energies(ib, ik) = -0.20_rp + 0.001_rp*real(ib + ik, rp)
            else
               energies(ib, ik) = 0.20_rp + 0.001_rp*real(ib + ik, rp)
            end if
            vectors(ib, ib, ik) = cmplx(1.0_rp, 0.0_rp, rp)
         end do
      end do
   end subroutine build_moment_eigensystem

end program test_lr_projected_site_spin
