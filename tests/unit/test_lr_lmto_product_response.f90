!------------------------------------------------------------------------------
! TDVK-02R1 production LMTO product-space and analytical transition oracle.
!
! LR-05 remains the independent point-grid reference.  This test verifies the
! candidate-space reconstruction and the retained weighted-SVD coordinate map
! without constructing or accumulating a susceptibility.
!------------------------------------------------------------------------------
program test_lr_lmto_product_response
   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use logger_mod, only: g_logger
   use math_mod, only: init_math_operators
   use self_mod, only: legacy_radial_fixture
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex
   use lr_response_space_mod, only: response_space_layout, response_vector_norm
   use lr_pauli_transition_vertex_mod, only: pauli_vertex_capabilities, pauli_endpoint_state, &
      pauli_sigma_plus_matrix, pauli_sigma_minus_matrix, evaluate_pauli_transition_vertex
   use lr_lmto_product_response_basis_mod, only: lmto_product_channel_plus, lmto_product_channel_minus, &
      lmto_product_candidate, lmto_product_response_basis
   implicit none

   integer, parameter :: nr = 51, lmax = 2, nsite = 1, norb = (lmax + 1)**2, nstate = 2*norb, ntransition = 4
   integer, parameter :: left_index(ntransition) = [1, 3, 5, 9]
   integer, parameter :: right_index(ntransition) = [2, 8, 12, 18]
   real(rp), parameter :: mesh_a = 0.03_rp, mesh_b = 0.10_rp, nuclear_z = 1.0_rp
   real(rp), parameter :: tolerance = 1.0e-10_rp
   real(rp) :: radius(nr)
   type(lmto_radial_basis) :: radial_sp, radial_spd
   type(response_space_layout) :: space_sp, space_spd, space_spd_reduced
   type(lmto_product_response_basis) :: product_sp, product_plus, product_minus, product_reduced
   type(pauli_vertex_capabilities) :: capabilities
   real(rp) :: eigenvalues(nstate)
   complex(rp) :: eigenvectors(nstate, nstate)
   logical :: failed, strict_guard
   character(len=32) :: mode
   real(rp) :: maximum_reconstruction, maximum_orthogonality, maximum_candidate_error
   real(rp) :: maximum_coordinate_error, maximum_norm_error

   call basis_init(lmax)
   call init_math_operators()
   call g_logger%init()
   call build_mesh(radius)
   call setup_radial_basis(radial_sp, radius, 1)
   call setup_radial_basis(radial_spd, radius, 2)
   call build_fixture_eigensystem(eigenvalues, eigenvectors)
   call space_sp%initialize(nsite, 2, radius, mesh_a, mesh_b, 1)
   call space_spd%initialize(nsite, 4, radius, mesh_a, mesh_b, 1)
   call space_spd_reduced%initialize(nsite, 2, radius, mesh_a, mesh_b, 1)

   mode = ''
   call get_command_argument(1, mode)
   strict_guard = trim(mode) == 'strict'
   if (len_trim(mode) > 0 .and. .not. strict_guard) error stop 'unknown UnitLrLmtoProductResponse argument'
   if (strict_guard) then
      ! The legacy fixture deliberately has threshold-sensitive tiny modes.
      ! Production initialization must fail closed rather than select a rank.
      call product_plus%initialize(space_spd, [radial_spd], lmto_product_channel_plus, .true.)
      error stop 'UnitLrLmtoProduct: strict rank guard unexpectedly returned'
   end if

   failed = .false.
   maximum_reconstruction = 0.0_rp
   maximum_orthogonality = 0.0_rp
   maximum_candidate_error = 0.0_rp
   maximum_coordinate_error = 0.0_rp
   maximum_norm_error = 0.0_rp
   call product_sp%initialize(space_sp, [radial_sp], lmto_product_channel_plus, .false.)
   call product_plus%initialize(space_spd, [radial_spd], lmto_product_channel_plus, .false.)
   call product_minus%initialize(space_spd, [radial_spd], lmto_product_channel_minus, .false.)
   call product_reduced%initialize(space_spd_reduced, [radial_spd], lmto_product_channel_plus, .false.)
   call check_inventory(product_sp, 1, [8, 8, 4], 52, failed)
   call check_inventory(product_plus, 2, [12, 16, 16, 8, 4], 232, failed)
   call check_inventory(product_minus, 2, [12, 16, 16, 8, 4], 232, failed)
   call check_inventory(product_reduced, 2, [12, 16, 16], 140, failed)
   call check_flat_mapping(product_sp, failed)
   call check_flat_mapping(product_plus, failed)
   call check_flat_mapping(product_minus, failed)
   call check_flat_mapping(product_reduced, failed)
   call check_svd_blocks(product_sp, radial_sp, space_sp, failed, maximum_reconstruction, maximum_orthogonality)
   call check_svd_blocks(product_plus, radial_spd, space_spd, failed, maximum_reconstruction, maximum_orthogonality)
   call check_svd_blocks(product_minus, radial_spd, space_spd, failed, maximum_reconstruction, maximum_orthogonality)
   capabilities = pauli_vertex_capabilities()
   call transition_oracle(product_plus, radial_spd, space_spd, capabilities, eigenvalues, eigenvectors, failed, &
      maximum_candidate_error, &
      maximum_coordinate_error, maximum_norm_error)
   call transition_oracle(product_minus, radial_spd, space_spd, capabilities, eigenvalues, eigenvectors, failed, &
      maximum_candidate_error, &
      maximum_coordinate_error, maximum_norm_error)

   if (failed) then
      write (*, '(a)') 'UnitLrLmtoProductResponse: FAIL'
      error stop 1
   end if
   write (*, '(a,es12.4)') '  maximum_candidate_space_error=', maximum_candidate_error
   write (*, '(a,es12.4)') '  maximum_orthonormal_coordinate_error=', maximum_coordinate_error
   write (*, '(a,es12.4)') '  maximum_norm_closure_error=', maximum_norm_error
   write (*, '(a,es12.4)') '  maximum_retained_SVD_reconstruction_error=', maximum_reconstruction
   write (*, '(a,es12.4)') '  maximum_U_orthogonality_error=', maximum_orthogonality
   write (*, '(a)') 'UnitLrLmtoProductResponse: PASS (inventory, mapping, SVD, LR-05 candidate/z/norm oracles)'

contains

   subroutine build_mesh(mesh)
      real(rp), intent(out) :: mesh(:)
      integer :: ir

      mesh(1) = 0.0_rp
      do ir = 2, size(mesh)
         mesh(ir) = mesh_b*(exp(mesh_a*real(ir - 1, rp)) - 1.0_rp)
      end do
   end subroutine build_mesh

   subroutine setup_radial_basis(basis, mesh, basis_lmax)
      type(lmto_radial_basis), intent(out) :: basis
      real(rp), intent(in) :: mesh(:)
      integer, intent(in) :: basis_lmax
      real(rp) :: potential(size(mesh)), energy
      real(rp), allocatable :: g(:, :), gp(:, :), gpp(:, :)
      real(rp), allocatable :: gpack(:), gdotpack(:), gddotpack(:)
      integer :: ispin, l

      potential = 0.0_rp
      call basis%initialize(size(mesh), basis_lmax, 2)
      do ispin = 1, 2
         do l = 0, basis_lmax
            call legacy_radial_fixture(nuclear_z, l, mesh_a, mesh_b, mesh, potential, energy, g, gp, gpp, &
               2.0_rp)
            gpack = reshape(g, [2*size(mesh)])
            gdotpack = reshape(gp, [2*size(mesh)])
            gddotpack = reshape(gpp, [2*size(mesh)])
            call basis%capture_channel(l, ispin, energy, mesh, potential, mesh_a, mesh_b, nuclear_z, &
               gpack, gdotpack, gddotpack, energy)
         end do
      end do
   end subroutine setup_radial_basis

   subroutine build_fixture_eigensystem(eigenvalues, eigenvectors)
      real(rp), intent(out) :: eigenvalues(:)
      complex(rp), intent(out) :: eigenvectors(:, :)
      complex(rp) :: hamiltonian(size(eigenvalues), size(eigenvalues)), work_query(1)
      complex(rp), allocatable :: work(:)
      real(rp), allocatable :: rwork(:)
      integer :: i, j, lwork, info
      external :: zheev

      hamiltonian = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, size(eigenvalues)
         hamiltonian(i, i) = cmplx(-0.85_rp + 0.11_rp*real(i, rp), 0.0_rp, rp)
         do j = i + 1, size(eigenvalues)
            hamiltonian(i, j) = cmplx(0.0013_rp*real(i + j, rp), 0.0007_rp*real(j - i, rp), rp)
            hamiltonian(j, i) = conjg(hamiltonian(i, j))
         end do
      end do
      eigenvectors = hamiltonian
      allocate(rwork(max(1, 3*size(eigenvalues) - 2)))
      call zheev('V', 'U', size(eigenvalues), eigenvectors, size(eigenvalues), eigenvalues, work_query, -1, rwork, info)
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      call zheev('V', 'U', size(eigenvalues), eigenvectors, size(eigenvalues), eigenvalues, work, lwork, rwork, info)
      deallocate(work, rwork)
      if (info /= 0) error stop 'UnitLrLmtoProduct: fixture eigensystem failed'
   end subroutine build_fixture_eigensystem

   subroutine check_inventory(product, basis_lmax, expected_counts, expected_total, failed)
      type(lmto_product_response_basis), intent(in) :: product
      integer, intent(in) :: basis_lmax, expected_counts(:), expected_total
      logical, intent(inout) :: failed
      integer :: response_l, k, pair_index, p, q, npair
      integer :: pair_l(9), pair_lp(9)

      if (product%unpruned_dimension /= expected_total .or. product%product_dimension < 1) failed = .true.
      if (product%response_lmax /= size(expected_counts) - 1) failed = .true.
      do response_l = 0, product%response_lmax
         if (product%blocks(1, response_l)%ncandidate /= expected_counts(response_l + 1)) failed = .true.
         call expected_pairs(basis_lmax, response_l, pair_l, pair_lp, npair)
         if (4*npair /= product%blocks(1, response_l)%ncandidate) failed = .true.
         do k = 1, product%blocks(1, response_l)%ncandidate
            pair_index = (k - 1)/4 + 1
            p = (mod(k - 1, 4))/2
            q = mod(k - 1, 2)
            if (product%blocks(1, response_l)%candidates(k)%l /= pair_l(pair_index) .or. &
                product%blocks(1, response_l)%candidates(k)%lp /= pair_lp(pair_index) .or. &
                product%blocks(1, response_l)%candidates(k)%p /= p .or. &
                product%blocks(1, response_l)%candidates(k)%q /= q) then
               failed = .true.
            end if
         end do
      end do
      write (*, '(a,i0,a,i0,a,i0)') 'INVENTORY lmax=', basis_lmax, ' unpruned=', product%unpruned_dimension, &
         ' retained=', product%product_dimension
   end subroutine check_inventory

   subroutine expected_pairs(basis_lmax, response_l, pair_l, pair_lp, npair)
      integer, intent(in) :: basis_lmax, response_l
      integer, intent(out) :: pair_l(:), pair_lp(:), npair

      pair_l = 0
      pair_lp = 0
      select case (basis_lmax)
      case (1)
         select case (response_l)
         case (0)
            pair_l(1:2) = [0, 1]
            pair_lp(1:2) = [0, 1]
            npair = 2
         case (1)
            pair_l(1:2) = [0, 1]
            pair_lp(1:2) = [1, 0]
            npair = 2
         case (2)
            pair_l(1) = 1
            pair_lp(1) = 1
            npair = 1
         end select
      case (2)
         select case (response_l)
         case (0)
            pair_l(1:3) = [0, 1, 2]
            pair_lp(1:3) = [0, 1, 2]
            npair = 3
         case (1)
            pair_l(1:4) = [0, 1, 1, 2]
            pair_lp(1:4) = [1, 0, 2, 1]
            npair = 4
         case (2)
            pair_l(1:4) = [0, 1, 2, 2]
            pair_lp(1:4) = [2, 1, 0, 2]
            npair = 4
         case (3)
            pair_l(1:2) = [1, 2]
            pair_lp(1:2) = [2, 1]
            npair = 2
         case (4)
            pair_l(1) = 2
            pair_lp(1) = 2
            npair = 1
         end select
      case default
         error stop 'UnitLrLmtoProduct: unsupported inventory fixture'
      end select
   end subroutine expected_pairs

   subroutine check_flat_mapping(product, failed)
      type(lmto_product_response_basis), intent(in) :: product
      logical, intent(inout) :: failed
      integer :: flat, site, response_l, response_m, product_mode

      do flat = 1, product%product_dimension
         call product%unflatten_index(flat, site, response_l, response_m, product_mode)
         if (product%flat_index(site, response_l, response_m, product_mode) /= flat) failed = .true.
      end do
      write (*, '(a,i0,a)') 'MAPPING roundtrip coordinates=', product%product_dimension, ' PASS'
   end subroutine check_flat_mapping

   subroutine check_svd_blocks(product, radial, space, failed, maximum_reconstruction, maximum_orthogonality)
      type(lmto_product_response_basis), intent(in) :: product
      type(lmto_radial_basis), intent(in) :: radial
      type(response_space_layout), intent(in) :: space
      logical, intent(inout) :: failed
      real(rp), intent(inout) :: maximum_reconstruction, maximum_orthogonality
      complex(rp), allocatable :: candidate_matrix(:, :), reconstructed(:, :), overlap(:, :)
      real(rp) :: error, scale
      integer :: response_l, i, j

      do response_l = 0, product%response_lmax
         allocate(candidate_matrix(space%npoint, product%blocks(1, response_l)%ncandidate))
         allocate(reconstructed(size(candidate_matrix, 1), size(candidate_matrix, 2)))
         call build_weighted_candidate_matrix(product, radial, space, response_l, candidate_matrix)
         reconstructed = cmplx(0.0_rp, 0.0_rp, rp)
         do j = 1, size(reconstructed, 2)
            do i = 1, product%blocks(1, response_l)%rank
               reconstructed(:, j) = reconstructed(:, j) + product%blocks(1, response_l)%weighted_modes(:, i)* &
                  product%blocks(1, response_l)%singular_values(i)*product%blocks(1, response_l)%right_modes(i, j)
            end do
         end do
         scale = max(sqrt(sum(real(candidate_matrix*conjg(candidate_matrix), rp))), epsilon(1.0_rp))
         error = sqrt(sum(real((reconstructed - candidate_matrix)*conjg(reconstructed - candidate_matrix), rp)))/scale
         maximum_reconstruction = max(maximum_reconstruction, error)
         if (error > tolerance) failed = .true.
         overlap = matmul(conjg(transpose(product%blocks(1, response_l)%weighted_modes)), &
            product%blocks(1, response_l)%weighted_modes)
         do i = 1, size(overlap, 1)
            overlap(i, i) = overlap(i, i) - cmplx(1.0_rp, 0.0_rp, rp)
         end do
         error = maxval(abs(overlap))
         maximum_orthogonality = max(maximum_orthogonality, error)
         if (error > tolerance) failed = .true.
         write (*, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'SVD L=', response_l, ' rank=', &
            product%blocks(1, response_l)%rank, ' sensitivity=', product%blocks(1, response_l)%rank_tau1, '/', &
            product%blocks(1, response_l)%rank_tau10, '/', product%blocks(1, response_l)%rank_tau100
         deallocate(candidate_matrix, reconstructed, overlap)
      end do
   end subroutine check_svd_blocks

   subroutine build_weighted_candidate_matrix(product, radial, space, response_l, matrix)
      type(lmto_product_response_basis), intent(in) :: product
      type(lmto_radial_basis), intent(in) :: radial
      type(response_space_layout), intent(in) :: space
      integer, intent(in) :: response_l
      complex(rp), intent(out) :: matrix(:, :)
      integer :: ir, k
      real(rp) :: value

      do k = 1, size(matrix, 2)
         do ir = 1, size(matrix, 1)
            value = independent_radial_product(radial, ir, product%blocks(1, response_l)%candidates(k), &
               product%circular_channel)
            matrix(ir, k) = cmplx(sqrt(space%radial_weights(ir))*value/ &
               product%blocks(1, response_l)%column_norms(k), 0.0_rp, rp)
         end do
      end do
   end subroutine build_weighted_candidate_matrix

   subroutine transition_oracle(product, radial, space, capabilities, eigenvalues, eigenvectors, failed, &
                                maximum_candidate_error, &
                                maximum_coordinate_error, maximum_norm_error)
      type(lmto_product_response_basis), intent(in) :: product
      type(lmto_radial_basis), intent(in) :: radial
      type(response_space_layout), intent(in) :: space
      type(pauli_vertex_capabilities), intent(in) :: capabilities
      real(rp), intent(in) :: eigenvalues(:)
      complex(rp), intent(in) :: eigenvectors(:, :)
      logical, intent(inout) :: failed
      real(rp), intent(inout) :: maximum_candidate_error, maximum_coordinate_error, maximum_norm_error
      type(pauli_endpoint_state) :: left_state, right_state
      complex(rp), allocatable :: reference(:), candidate(:), z_reference(:), z_product(:)
      complex(rp) :: operator_matrix(2, 2)
      integer :: pair
      real(rp) :: candidate_error, coordinate_error, norm_error

      allocate(reference(space%ndim), &
         candidate(space%ndim), z_reference(product%product_dimension), z_product(product%product_dimension))
      if (product%circular_channel == lmto_product_channel_plus) then
         operator_matrix = pauli_sigma_plus_matrix()
      else
         operator_matrix = pauli_sigma_minus_matrix()
      end if
      do pair = 1, ntransition
         call left_state%initialize(eigenvalues(left_index(pair)), eigenvectors(:, left_index(pair)))
         call right_state%initialize(eigenvalues(right_index(pair)), eigenvectors(:, right_index(pair)))
         call evaluate_pauli_transition_vertex(space, [radial], left_state, right_state, operator_matrix, &
            capabilities, reference)
         call build_candidate_transition(product, radial, space, left_state, right_state, candidate)
         call build_reference_coordinates(product, space, reference, z_reference)
         call product%transition_coordinates(left_state, right_state, z_product)
         candidate_error = response_vector_norm(space, candidate - reference)/ &
            max(response_vector_norm(space, reference), epsilon(1.0_rp))
         coordinate_error = sqrt(sum(abs(z_product - z_reference)**2))/ &
            max(sqrt(sum(abs(z_reference)**2)), epsilon(1.0_rp))
         norm_error = abs(response_vector_norm(space, reference)**2 - sum(abs(z_reference)**2))/ &
            max(response_vector_norm(space, reference)**2, epsilon(1.0_rp))
         maximum_candidate_error = max(maximum_candidate_error, candidate_error)
         maximum_coordinate_error = max(maximum_coordinate_error, coordinate_error)
         maximum_norm_error = max(maximum_norm_error, norm_error)
         if (candidate_error >= tolerance .or. coordinate_error >= tolerance .or. norm_error >= tolerance) then
            failed = .true.
         end if
         write (*, '(a,a,a,i0,3(a,es12.4))') 'ORACLE channel=', channel_name(product%circular_channel), &
            ' pair=', pair, ' candidate=', candidate_error, ' z=', coordinate_error, ' norm=', norm_error
      end do
      deallocate(reference, candidate, z_reference, z_product)
   end subroutine transition_oracle

   subroutine build_candidate_transition(product, radial, space, left_state, right_state, transition)
      type(lmto_product_response_basis), intent(in) :: product
      type(lmto_radial_basis), intent(in) :: radial
      type(response_space_layout), intent(in) :: space
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(out) :: transition(:)
      type(response_super_index) :: item
      complex(rp), allocatable :: coefficients(:)
      integer :: site, response_l, response_m, ir, k, flat

      transition = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, product%nsite
         do response_l = 0, product%response_lmax
            allocate(coefficients(product%blocks(site, response_l)%ncandidate))
            do response_m = -response_l, response_l
               call product%candidate_coefficients(site, response_l, response_m, left_state, right_state, coefficients)
               do ir = 1, space%npoint
                  item = response_super_index(site, response_l, response_m, ir, 1)
                  call response_flatten_superindex(item, space%nsite, space%response_lmax, space%npoint, &
                     space%nchannel, flat)
                  do k = 1, size(coefficients)
                     transition(flat) = transition(flat) + coefficients(k)*cmplx( &
                        independent_radial_product(radial, ir, product%blocks(site, response_l)%candidates(k), &
                        product%circular_channel), 0.0_rp, rp)
                  end do
               end do
            end do
            deallocate(coefficients)
         end do
      end do
   end subroutine build_candidate_transition

   subroutine build_reference_coordinates(product, space, transition, coordinates)
      type(lmto_product_response_basis), intent(in) :: product
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: transition(:)
      complex(rp), intent(out) :: coordinates(:)
      type(response_super_index) :: item
      integer :: site, response_l, response_m, product_mode, ir, flat, response_flat

      coordinates = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, product%nsite
         do response_l = 0, product%response_lmax
            do response_m = -response_l, response_l
               do product_mode = 1, product%blocks(site, response_l)%rank
                  flat = product%flat_index(site, response_l, response_m, product_mode)
                  do ir = 1, space%npoint
                     item = response_super_index(site, response_l, response_m, ir, 1)
                     call response_flatten_superindex(item, space%nsite, space%response_lmax, space%npoint, &
                        space%nchannel, response_flat)
                     coordinates(flat) = coordinates(flat) + conjg(product%blocks(site, response_l)%weighted_modes(ir, &
                        product_mode))*sqrt(space%radial_weights(ir))*transition(response_flat)
                  end do
               end do
            end do
         end do
      end do
   end subroutine build_reference_coordinates

   real(rp) function independent_radial_product(radial, ir, candidate, circular_channel) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, circular_channel
      type(lmto_product_candidate), intent(in) :: candidate
      integer :: spin_left, spin_right
      real(rp) :: first, second, rfirst, rsecond

      if (circular_channel == lmto_product_channel_plus) then
         spin_left = 1
         spin_right = 2
      else
         spin_left = 2
         spin_right = 1
      end if
      if (ir /= 1) then
         value = endpoint(radial, ir, candidate%l, spin_left, candidate%p)* &
            endpoint(radial, ir, candidate%lp, spin_right, candidate%q)/radial%rofi(ir)**2
         return
      end if
      if (candidate%l /= 0 .or. candidate%lp /= 0) then
         value = 0.0_rp
         return
      end if
      rfirst = radial%rofi(2)
      rsecond = radial%rofi(3)
      first = endpoint(radial, 2, candidate%l, spin_left, candidate%p)* &
         endpoint(radial, 2, candidate%lp, spin_right, candidate%q)/rfirst**2
      second = endpoint(radial, 3, candidate%l, spin_left, candidate%p)* &
         endpoint(radial, 3, candidate%lp, spin_right, candidate%q)/rsecond**2
      value = (first*rsecond**2 - second*rfirst**2)/(rsecond**2 - rfirst**2)
   end function independent_radial_product

   pure real(rp) function endpoint(radial, ir, l, spin, power) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, spin, power

      if (power == 0) then
         value = radial%phi_large(ir, l + 1, spin) - radial%enu_work(l + 1, spin)* &
            radial%phidot_large(ir, l + 1, spin)
      else
         value = radial%phidot_large(ir, l + 1, spin)
      end if
   end function endpoint

   character(len=5) function channel_name(circular_channel) result(value)
      integer, intent(in) :: circular_channel

      if (circular_channel == lmto_product_channel_plus) then
         value = 'plus '
      else
         value = 'minus'
      end if
   end function channel_name

end program test_lr_lmto_product_response
