!------------------------------------------------------------------------------
! TDVK-02R2 compact LMTO-product Lehmann susceptibility.
!
! The legacy LR-06 point-space evaluator is the independent operator oracle.
! This test projects its canonical result into the weighted-orthonormal product
! basis and compares it with the compact transition-coordinate accumulation.
!------------------------------------------------------------------------------
program test_lr_product_ks_susceptibility
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use logger_mod, only: g_logger
   use math_mod, only: init_math_operators
   use self_mod, only: legacy_radial_fixture
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex
   use lr_response_space_mod, only: response_space_layout
   use lr_pauli_transition_vertex_mod, only: pauli_endpoint_state
   use lr_ks_susceptibility_mod, only: lr_electronic_state, lr_ks_susceptibility_request, &
      lr_ks_susceptibility_result, lr_product_ks_susceptibility_request, lr_product_ks_susceptibility_result, &
      lr_channel_plus, lr_channel_minus, evaluate_lr_ks_susceptibility, evaluate_lr_product_ks_susceptibility
   use lr_lmto_product_response_basis_mod, only: lmto_product_channel_plus, lmto_product_channel_minus, &
      lmto_product_response_basis
   implicit none

   integer, parameter :: nr = 51, lmax = 1, nsite = 1, norb = (lmax + 1)**2
   integer, parameter :: nbands = 2*norb, nk = 2, nfrequency = 2
   real(rp), parameter :: mesh_a = 0.03_rp, mesh_b = 0.10_rp, nuclear_z = 1.0_rp
   real(rp), parameter :: frequencies(nfrequency) = [0.0_rp, 0.23_rp]
   real(rp), parameter :: eta_values(2) = [0.02_rp, 0.04_rp]
   real(rp), parameter :: q_zero(3) = [0.0_rp, 0.0_rp, 0.0_rp]
   real(rp), parameter :: q_finite(3) = [0.17_rp, -0.08_rp, 0.13_rp]
   real(rp), parameter :: k_points(3, nk) = reshape([ &
      0.0_rp, 0.0_rp, 0.0_rp, 0.25_rp, 0.10_rp, -0.10_rp], [3, nk])
   real(rp), parameter :: weights(nk) = [1.0_rp, 3.0_rp]
   real(rp), parameter :: occupations(nbands, nk) = reshape([ &
      1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
      1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], [nbands, nk])
   real(rp) :: radius(nr), left_energies(nbands, nk), right_energies(nbands, nk)
   real(rp) :: right_points(3, nk), q(3)
   complex(rp) :: left_vectors(norb*2, nbands, nk), right_vectors(norb*2, nbands, nk)
   type(lmto_radial_basis), target :: radial(1)
   type(response_space_layout), target :: space
   type(lmto_product_response_basis), target :: product_plus, product_minus
   type(lr_electronic_state), target :: left_state, right_state
   type(lr_ks_susceptibility_request) :: legacy_request
   type(lr_product_ks_susceptibility_request) :: product_request
   type(lr_ks_susceptibility_result) :: legacy_result
   type(lr_product_ks_susceptibility_result) :: product_result
   real(rp) :: max_equivalence, max_pair, max_immutable
   logical :: failed
   integer :: channel_kind, iq, ieta

   call basis_init(lmax)
   call init_math_operators()
   call g_logger%init()
   call build_mesh(radius)
   call setup_radial_basis(radial(1), radius)
   call build_fixture_eigensystem(left_energies(:, 1), left_vectors(:, :, 1))
   left_energies(:, 2) = left_energies(:, 1) + [0.00_rp, 0.01_rp, -0.01_rp, 0.02_rp, 0.03_rp, -0.02_rp, 0.01_rp, -0.03_rp]
   left_vectors(:, :, 2) = left_vectors(:, :, 1)
   right_energies = left_energies + spread([0.013_rp, -0.007_rp, 0.009_rp, -0.011_rp, 0.017_rp, -0.005_rp, &
      0.003_rp, 0.019_rp], dim=2, ncopies=nk)
   right_vectors = left_vectors
   call left_state%initialize(left_energies, left_vectors, k_points, weights, occupations, 0.0_rp, 300.0_rp)
   call space%initialize(nsite, 2*lmax, radius, mesh_a, mesh_b, 1)
   call product_plus%initialize(space, radial, lmto_product_channel_plus, .true.)
   call product_minus%initialize(space, radial, lmto_product_channel_minus, .true.)

   if (product_plus%unpruned_dimension /= 52 .or. product_plus%product_dimension < 1 .or. &
       product_minus%unpruned_dimension /= 52 .or. product_minus%product_dimension /= product_plus%product_dimension) then
      error stop 'UnitLrProductKsSusceptibility: product dimensions are inconsistent'
   end if
   write (*, '(a,i0,a,i0)') 'PRODUCT dimensions unpruned=', product_plus%unpruned_dimension, &
      ' retained=', product_plus%product_dimension

   max_equivalence = 0.0_rp
   max_pair = 0.0_rp
   max_immutable = 0.0_rp
   failed = .false.
   do channel_kind = 1, 2
      do iq = 1, 2
         if (iq == 1) then
            q = q_zero
         else
            q = q_finite
         end if
         right_points = 0.0_rp
         do ieta = 1, nk
            right_points(:, ieta) = fold_fractional_kpoint(k_points(:, ieta) + q)
         end do
         call right_state%initialize(right_energies, right_vectors, right_points, weights, occupations, 0.0_rp, 300.0_rp)

         do ieta = 1, size(eta_values)
            call run_comparison(channel_kind, q, eta_values(ieta), failed, max_equivalence, max_pair, max_immutable)
         end do
      end do
   end do

   if (failed) then
      write (*, '(a)') 'UnitLrProductKsSusceptibility: FAIL'
      error stop 1
   end if
   write (*, '(a,es12.4)') 'maximum_fixture_equivalence_error=', max_equivalence
   write (*, '(a,es12.4)') 'maximum_pair_denominator_error=', max_pair
   write (*, '(a,es12.4)') 'maximum_snapshot_mutation=', max_immutable
   write (*, '(a)') 'UnitLrProductKsSusceptibility: PASS (compact Lehmann, LR-06 projection, q/omega/eta/channels)'

contains

   subroutine build_mesh(mesh)
      real(rp), intent(out) :: mesh(:)
      integer :: ir

      mesh(1) = 0.0_rp
      do ir = 2, size(mesh)
         mesh(ir) = mesh_b*(exp(mesh_a*real(ir - 1, rp)) - 1.0_rp)
      end do
   end subroutine build_mesh

   subroutine setup_radial_basis(basis, mesh)
      type(lmto_radial_basis), intent(out) :: basis
      real(rp), intent(in) :: mesh(:)
      real(rp) :: potential(size(mesh)), energy
      real(rp), allocatable :: g(:, :), gp(:, :), gpp(:, :)
      real(rp), allocatable :: gpack(:), gdotpack(:), gddotpack(:)
      integer :: ispin, l

      potential = 0.0_rp
      call basis%initialize(size(mesh), lmax, 2)
      do ispin = 1, 2
         do l = 0, lmax
            call legacy_radial_fixture(nuclear_z, l, mesh_a, mesh_b, mesh, potential, energy, g, gp, gpp, 2.0_rp)
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
         hamiltonian(i, i) = cmplx(-0.80_rp + 0.17_rp*real(i, rp), 0.0_rp, rp)
         do j = i + 1, size(eigenvalues)
            hamiltonian(i, j) = cmplx(0.0011_rp*real(i + 2*j, rp), 0.0006_rp*real(j - i, rp), rp)
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
      if (info /= 0) error stop 'UnitLrProductKsSusceptibility: fixture eigensystem failed'
   end subroutine build_fixture_eigensystem

   subroutine run_comparison(channel_kind, q_in, eta_in, failed_in, max_equiv_in, max_pair_in, max_immutable_in)
      integer, intent(in) :: channel_kind
      real(rp), intent(in) :: q_in(3), eta_in
      logical, intent(inout) :: failed_in
      real(rp), intent(inout) :: max_equiv_in, max_pair_in, max_immutable_in
      complex(rp), allocatable :: reference(:, :, :), independent(:, :, :)
      real(rp), allocatable :: saved_energies(:, :), saved_occupations(:, :), saved_weights(:)
      complex(rp), allocatable :: saved_vectors(:, :, :)
      real(rp) :: equiv_error, pair_error
      character(len=32) :: channel
      integer :: expected_transitions, expected_skips

      if (channel_kind == 1) then
         channel = lr_channel_plus
         product_request%product_basis => product_plus
      else
         channel = lr_channel_minus
         product_request%product_basis => product_minus
      end if
      legacy_request%q = q_in
      legacy_request%frequencies = frequencies
      legacy_request%eta = eta_in
      legacy_request%channel = channel
      legacy_request%response_space => space
      legacy_request%radial_bases => radial
      legacy_request%electronic_state => left_state
      legacy_request%q_endpoint_state => right_state
      product_request%q = q_in
      product_request%frequencies = frequencies
      product_request%eta = eta_in
      product_request%channel = channel
      product_request%electronic_state => left_state
      product_request%q_endpoint_state => right_state

      allocate(saved_energies(size(left_state%eigenvalues, 1), size(left_state%eigenvalues, 2)), &
         saved_occupations(size(left_state%occupations, 1), size(left_state%occupations, 2)), &
         saved_weights(size(left_state%k_weights)), saved_vectors(size(left_state%eigenvectors, 1), &
         size(left_state%eigenvectors, 2), size(left_state%eigenvectors, 3)))
      saved_energies = left_state%eigenvalues
      saved_occupations = left_state%occupations
      saved_weights = left_state%k_weights
      saved_vectors = left_state%eigenvectors

      call evaluate_lr_ks_susceptibility(legacy_request, legacy_result)
      call evaluate_lr_product_ks_susceptibility(product_request, product_result)
      expected_transitions = nk*2*4*4
      expected_skips = nk*2*4*4
      if (product_result%product_dimension /= product_plus%product_dimension .or. &
          product_result%transition_dimension /= product_plus%product_dimension .or. &
          product_result%point_space_transition_allocated .or. &
          product_result%ntransitions_evaluated /= expected_transitions .or. &
          product_result%noccupation_skips /= expected_skips) then
         failed_in = .true.
      end if
      if (index(product_result%response_representation, 'weighted-orthonormal LMTO product') == 0 .or. &
          index(product_result%response_space_metadata, 'no point-space transition allocation') == 0) then
         failed_in = .true.
      end if

      allocate(reference(product_result%product_dimension, product_result%product_dimension, nfrequency), &
         independent(product_result%product_dimension, product_result%product_dimension, nfrequency))
      call project_legacy_result(legacy_result%susceptibility, product_request%product_basis, reference)
      call independent_product_accumulation(product_request%product_basis, product_request, independent)
      equiv_error = sqrt(sum(abs(product_result%susceptibility - reference)**2))/ &
         max(sqrt(sum(abs(reference)**2)), epsilon(1.0_rp))
      max_equiv_in = max(max_equiv_in, equiv_error)
      if (equiv_error >= 1.0e-10_rp) failed_in = .true.
      pair_error = sqrt(sum(abs(product_result%susceptibility - independent)**2))/ &
         max(sqrt(sum(abs(independent)**2)), epsilon(1.0_rp))
      max_pair_in = max(max_pair_in, pair_error)
      if (pair_error >= 1.0e-12_rp) failed_in = .true.

      max_immutable_in = max(max_immutable_in, maxval(abs(left_state%eigenvalues - saved_energies)))
      max_immutable_in = max(max_immutable_in, maxval(abs(left_state%occupations - saved_occupations)))
      max_immutable_in = max(max_immutable_in, maxval(abs(left_state%k_weights - saved_weights)))
      max_immutable_in = max(max_immutable_in, maxval(abs(left_state%eigenvectors - saved_vectors)))
      if (max_immutable_in /= 0.0_rp) failed_in = .true.
      write (*, '(5a,es12.4,a,es12.4,a,es12.4,a,i0)') 'ORACLE channel=', trim(channel), &
         ' qkind=', trim(merge('Gamma ', 'finite', maxval(abs(q_in)) < 1.0e-12_rp)), ' eta=', eta_in, &
         ' projection=', equiv_error, ' pair=', pair_error, ' compact-dim=', product_result%product_dimension
      deallocate(reference, independent, saved_energies, saved_occupations, saved_weights, saved_vectors)
   end subroutine run_comparison

   subroutine project_legacy_result(canonical, product, projected)
      complex(rp), intent(in) :: canonical(:, :, :)
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(out) :: projected(:, :, :)
      complex(rp), allocatable :: full_modes(:, :), vector(:), image(:)
      real(rp), allocatable :: metric(:)
      type(response_super_index) :: item
      integer :: alpha, beta, point_flat, product_flat, site, response_l, response_m, product_mode, ir, ifrequency

      allocate(full_modes(space%ndim, product%product_dimension), vector(space%ndim), image(space%ndim), &
         metric(space%ndim))
      metric = space%metric_weights
      full_modes = cmplx(0.0_rp, 0.0_rp, rp)
      do product_flat = 1, product%product_dimension
         call product%unflatten_index(product_flat, site, response_l, response_m, product_mode)
         do ir = 1, space%npoint
            item = response_super_index(site, response_l, response_m, ir, 1)
            call response_flatten_superindex(item, space%nsite, space%response_lmax, space%npoint, &
               space%nchannel, point_flat)
            full_modes(point_flat, product_flat) = product%blocks(site, response_l)%weighted_modes(ir, product_mode)
         end do
      end do
      projected = cmplx(0.0_rp, 0.0_rp, rp)
      do ifrequency = 1, size(canonical, 3)
         do beta = 1, product%product_dimension
            vector = cmplx(0.0_rp, 0.0_rp, rp)
            do point_flat = 1, space%ndim
               if (metric(point_flat) > 0.0_rp) vector(point_flat) = full_modes(point_flat, beta)/sqrt(metric(point_flat))
            end do
            image = matmul(canonical(:, :, ifrequency), vector)
            do alpha = 1, product%product_dimension
               projected(alpha, beta, ifrequency) = sum(conjg(full_modes(:, alpha))*sqrt(metric)*image)
            end do
         end do
      end do
      deallocate(full_modes, vector, image, metric)
   end subroutine project_legacy_result

   subroutine independent_product_accumulation(product, request, accumulated)
      type(lmto_product_response_basis), intent(in) :: product
      type(lr_product_ks_susceptibility_request), intent(in) :: request
      complex(rp), intent(out) :: accumulated(:, :, :)
      type(pauli_endpoint_state) :: left_band, right_band
      complex(rp), allocatable :: transition(:)
      complex(rp) :: denominator, pair_factor
      real(rp) :: weight_sum, occupation_difference
      integer :: ik, ib, jb, ifrequency, i, j

      allocate(transition(product%product_dimension))
      accumulated = cmplx(0.0_rp, 0.0_rp, rp)
      weight_sum = sum(left_state%k_weights)
      do ik = 1, left_state%nk
         do ib = 1, left_state%nbands
            call left_band%initialize(left_state%eigenvalues(ib, ik), left_state%eigenvectors(:, ib, ik))
            do jb = 1, right_state%nbands
               occupation_difference = left_state%occupations(ib, ik) - right_state%occupations(jb, ik)
               if (occupation_difference == 0.0_rp) cycle
               call right_band%initialize(right_state%eigenvalues(jb, ik), right_state%eigenvectors(:, jb, ik))
               call product%transition_coordinates(left_band, right_band, transition)
               do ifrequency = 1, size(request%frequencies)
                  denominator = cmplx(request%frequencies(ifrequency) + left_band%energy - right_band%energy, &
                     request%eta, rp)
                  pair_factor = cmplx(2.0_rp*left_state%k_weights(ik)/weight_sum, 0.0_rp, rp)* &
                     occupation_difference/denominator
                  do j = 1, product%product_dimension
                     do i = 1, product%product_dimension
                        accumulated(i, j, ifrequency) = accumulated(i, j, ifrequency) + &
                           pair_factor*transition(i)*conjg(transition(j))
                     end do
                  end do
               end do
            end do
         end do
      end do
      deallocate(transition)
   end subroutine independent_product_accumulation

   pure function fold_fractional_kpoint(k_point) result(folded)
      real(rp), intent(in) :: k_point(3)
      real(rp) :: folded(3)

      folded = k_point - floor(k_point + 0.5_rp)
   end function fold_fractional_kpoint

end program test_lr_product_ks_susceptibility
