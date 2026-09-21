!------------------------------------------------------------------------------
! DRESP-02 projected reciprocal bare susceptibility certification.
!
! The compact product evaluators are used only as independent projection
! oracles.  The object under test accumulates site transition outer products
! (Lehmann) or site-projected affine vertices in the real-axis Kubo bubble
! (GF).  Both fixtures use numerically diagonalized complex eigenvectors.
!------------------------------------------------------------------------------
program test_lr_projected_reciprocal_chi0
   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_response_space_mod, only: response_space_layout
   use lr_projected_site_spin_mod, only: projected_site_spin_contract
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, lmto_product_channel_plus, &
      lmto_product_channel_minus
   use lr_ks_susceptibility_mod, only: lr_electronic_state, lr_product_ks_susceptibility_request, &
      lr_product_ks_susceptibility_result, lr_channel_plus, lr_channel_minus, &
      lr_fermi_dirac_occupation, evaluate_lr_product_ks_susceptibility
   use lr_product_gf_susceptibility_mod, only: lr_product_gf_susceptibility_request, &
      lr_product_gf_susceptibility_result, lr_product_gf_contraction_factorized, &
      evaluate_lr_product_gf_susceptibility
   use lr_projected_reciprocal_chi0_mod, only: projected_chi0_request, projected_chi0_result, &
      evaluate_projected_lehmann_chi0, evaluate_projected_finite_width_chi0, evaluate_projected_gf_chi0, &
      evaluate_projected_gf_chi0_reference
   implicit none

   call basis_init(2)
   call run_fixture(1, 8, 'complex one-site fixture')
   call run_fixture(2, 8, 'mandatory two-site fixture')
   write (*, '(a)') 'UnitLrProjectedReciprocalChi0: PASS (DRESP-02 Lehmann/GF/finite-width/site matrix/oracles/q/eta)'

contains

   subroutine run_fixture(nsite, nbands, label)
      integer, intent(in) :: nsite, nbands
      character(len=*), intent(in) :: label

      integer, parameter :: nr = 7, orbital_lmax = 2, nfrequency = 2
      integer :: norb, one_spin_basis, nbasis
      real(rp), parameter :: mesh_a = 0.035_rp, mesh_b = 0.11_rp
      real(rp), parameter :: temperature = 0.01_rp, eta = 0.035_rp, energy_margin = 0.75_rp
      real(rp), parameter :: qfinite(3) = [0.23_rp, -0.11_rp, 0.07_rp]
      real(rp), parameter :: omega(nfrequency) = [0.0_rp, 0.17_rp]
      real(rp) :: radius(nr), eigenvalues(nbands, 1), endpoint_values(nbands, 1)
      real(rp) :: k_points(3, 1), endpoint_points(3, 1), negative_points(3, 1), weights(1), occupations(nbands, 1)
      complex(rp), allocatable :: eigenvectors(:, :, :), endpoint_vectors(:, :, :), negative_vectors(:, :, :)
      type(lmto_radial_basis), allocatable, target :: radial(:)
      type(response_space_layout), target :: space
      type(projected_site_spin_contract), target :: contract_d, contract_spd
      type(lmto_product_response_basis), target :: product_plus, product_minus
      type(lr_electronic_state), target :: state, endpoint_zero, endpoint_finite, endpoint_negative
      logical :: failed

      norb = (orbital_lmax + 1)**2
      one_spin_basis = norb*nsite
      nbasis = 2*one_spin_basis
      allocate(eigenvectors(nbasis, nbands, 1), endpoint_vectors(nbasis, nbands, 1), &
         negative_vectors(nbasis, nbands, 1), radial(nsite))
      call build_mesh(radius, mesh_a, mesh_b)
      call setup_radial_bases(radial, radius, mesh_a, mesh_b)
      call space%initialize(nsite, 2*orbital_lmax, radius, mesh_a, mesh_b, 1)
      call build_complex_spin_symmetric_eigenpairs(nsite, nbands, eigenvalues(:, 1), eigenvectors(:, :, 1))
      endpoint_values = eigenvalues
      k_points = 0.0_rp
      endpoint_points = 0.0_rp
      negative_points(:, 1) = -qfinite
      weights = 1.0_rp
      call build_occupations(eigenvalues, occupations, temperature)
      call apply_site_phase(eigenvectors(:, :, 1), nsite, qfinite, endpoint_vectors(:, :, 1))
      call apply_site_phase(eigenvectors(:, :, 1), nsite, -qfinite, negative_vectors(:, :, 1))

      call state%initialize(eigenvalues, eigenvectors, k_points, weights, occupations, 0.0_rp, temperature)
      call endpoint_zero%initialize(endpoint_values, endpoint_vectors, endpoint_points, weights, occupations, &
         0.0_rp, temperature)
      call endpoint_finite%initialize(endpoint_values, endpoint_vectors, reshape(qfinite, [3, 1]), weights, &
         occupations, 0.0_rp, temperature)
      call endpoint_negative%initialize(endpoint_values, negative_vectors, negative_points, weights, occupations, &
         0.0_rp, temperature)
      call contract_d%initialize(space, radial, 'd')
      call contract_spd%initialize(space, radial, 'spd')
      call product_plus%initialize(space, radial, lmto_product_channel_plus, .false.)
      call product_minus%initialize(space, radial, lmto_product_channel_minus, .false.)

      failed = .false.
      write (*, '(a)') trim(label)
      call run_projection_case(contract_d, product_plus, state, endpoint_zero, [0.0_rp, 0.0_rp, 0.0_rp], omega, eta, energy_margin, &
         'd Gamma', failed, .true.)
      call run_projection_case(contract_spd, product_plus, state, endpoint_zero, [0.0_rp, 0.0_rp, 0.0_rp], omega, eta, energy_margin, &
         'spd Gamma', failed, .true.)
      call run_projection_case(contract_d, product_plus, state, endpoint_finite, qfinite, omega, eta, energy_margin, &
         'd finite-q', failed, .false.)
      call run_projection_case(contract_spd, product_plus, state, endpoint_finite, qfinite, omega, eta, energy_margin, &
         'spd finite-q', failed, .false.)
      call run_covariance_case(contract_d, product_plus, product_minus, state, endpoint_zero, &
         [0.0_rp, 0.0_rp, 0.0_rp], omega(2), eta, energy_margin, 'd Gamma', failed)
      call run_covariance_case(contract_spd, product_plus, product_minus, state, endpoint_zero, &
         [0.0_rp, 0.0_rp, 0.0_rp], omega(2), eta, energy_margin, 'spd Gamma', failed)
      call run_covariance_case(contract_d, product_plus, product_minus, state, endpoint_finite, &
         qfinite, omega(2), eta, energy_margin, 'd finite-q', failed, endpoint_negative)
      call run_covariance_case(contract_spd, product_plus, product_minus, state, endpoint_finite, &
         qfinite, omega(2), eta, energy_margin, 'spd finite-q', failed, endpoint_negative)

      deallocate(eigenvectors, endpoint_vectors, negative_vectors, radial)
      if (failed) error stop 'UnitLrProjectedReciprocalChi0: fixture failed'
   end subroutine run_fixture

   subroutine run_projection_case(contract, product, state, endpoint, q, frequencies, eta_in, margin, label, failed, do_ladder)
      type(projected_site_spin_contract), intent(in) :: contract
      type(lmto_product_response_basis), intent(in), target :: product
      type(lr_electronic_state), intent(in), target :: state, endpoint
      real(rp), intent(in) :: q(3), frequencies(:), eta_in, margin
      character(len=*), intent(in) :: label
      logical, intent(inout) :: failed
      logical, intent(in) :: do_ladder

      type(projected_chi0_request) :: request
      type(projected_chi0_result) :: direct_lehmann, direct_gf, direct_gf_reference, finite_width_oracle, gf_coarse, gf_fine
      type(lr_product_ks_susceptibility_request) :: product_request
      type(lr_product_ks_susceptibility_result) :: product_lehmann
      type(lr_product_gf_susceptibility_request) :: product_gf_request
      type(lr_product_gf_susceptibility_result) :: product_gf
      complex(rp), allocatable :: oracle_lehmann(:, :, :), oracle_gf(:, :, :)
      real(rp) :: lehmann_residual, gf_residual, reference_residual, reference_term_one_residual, &
         reference_term_two_residual, coarse_difference, fine_difference, trend_residual
      real(rp) :: finite_width_residual
      real(rp) :: diagnostic_residual
      integer :: i

      call prepare_request(request, contract, product, state, endpoint, q, frequencies, eta_in, margin)
      request%diagnostics = .true.
      call evaluate_projected_lehmann_chi0(request, direct_lehmann)
      diagnostic_residual = maxval(abs(sum(direct_lehmann%k_susceptibility, dim=4) - direct_lehmann%susceptibility))
      if (.not. allocated(direct_lehmann%k_susceptibility) .or. direct_lehmann%n_dominant_transition_records < 1 .or. &
          diagnostic_residual > 2.0e-12_rp) then
         failed = .true.
         write (*, '(a,es12.4)') '  '//trim(label)//' Lehmann diagnostic k-sum residual = ', diagnostic_residual
      end if

      product_request%q = request%q
      product_request%frequencies = frequencies
      product_request%eta = eta_in
      product_request%channel = lr_channel_plus
      product_request%product_basis => product
      product_request%electronic_state => state
      product_request%q_endpoint_state => endpoint
      allocate(product_request%selected_l(0:2))
      product_request%selected_l = contract%selected_l
      call evaluate_lr_product_ks_susceptibility(product_request, product_lehmann)
      allocate(oracle_lehmann(contract%nsite, contract%nsite, size(frequencies)))
      call project_product_matrix(contract, product, product_lehmann%susceptibility, oracle_lehmann)
      lehmann_residual = maxval(abs(direct_lehmann%susceptibility - oracle_lehmann))
      write (*, '(a,es12.4)') '  '//trim(label)//' Lehmann oracle max|dchi| = ', lehmann_residual
      if (lehmann_residual > 2.0e-11_rp .or. direct_lehmann%point_response_allocated) failed = .true.

      ! The mandatory two-site fixture must exercise the complete site
      ! matrix, including both off-diagonal entries, rather than only its
      ! trace or diagonal block.
      if (contract%nsite == 2 .and. index(label, 'Gamma') > 0) then
         write (*, '(a)') '  '//trim(label)//' complete chi(site,site,omega_2):'
         do i = 1, contract%nsite
            write (*, '(a,2(i0,1x),2(es24.16,1x))') '    chi ', i, 1, &
               real(direct_lehmann%susceptibility(i, 1, 2), rp), aimag(direct_lehmann%susceptibility(i, 1, 2))
            write (*, '(a,2(i0,1x),2(es24.16,1x))') '    chi ', i, 2, &
               real(direct_lehmann%susceptibility(i, 2, 2), rp), aimag(direct_lehmann%susceptibility(i, 2, 2))
            if (.not. all(ieee_is_finite(real(direct_lehmann%susceptibility(i, 1:2, 2), rp))) .or. &
                .not. all(ieee_is_finite(aimag(direct_lehmann%susceptibility(i, 1:2, 2)))) .or. &
                any(abs(direct_lehmann%susceptibility(i, 1:2, 2)) <= tiny(1.0_rp))) then
               failed = .true.
            end if
         end do
      end if

      product_gf_request%q = request%q
      product_gf_request%frequencies = frequencies
      product_gf_request%eta = eta_in
      product_gf_request%channel = lr_channel_plus
      product_gf_request%integration_points = 101
      product_gf_request%integration_eta = 0.012_rp
      product_gf_request%energy_margin = margin
      product_gf_request%contraction_backend = lr_product_gf_contraction_factorized
      product_gf_request%product_basis => product
      product_gf_request%electronic_state => state
      product_gf_request%q_endpoint_state => endpoint
      allocate(product_gf_request%selected_l(0:2))
      product_gf_request%selected_l = contract%selected_l
      call evaluate_lr_product_gf_susceptibility(product_gf_request, product_gf)
      allocate(oracle_gf(contract%nsite, contract%nsite, size(frequencies)))
      call project_product_matrix(contract, product, product_gf%susceptibility, oracle_gf)

      request%integration_points = 101
      request%integration_eta = 0.012_rp
      call evaluate_projected_gf_chi0(request, direct_gf)
      call evaluate_projected_gf_chi0_reference(request, direct_gf_reference)
      request%diagnostics = .false.
      call evaluate_projected_finite_width_chi0(request, finite_width_oracle)
      reference_residual = maxval(abs(direct_gf%susceptibility - direct_gf_reference%susceptibility))
      reference_term_one_residual = maxval(abs(direct_gf%kubo_term_one - direct_gf_reference%kubo_term_one))
      reference_term_two_residual = maxval(abs(direct_gf%kubo_term_two - direct_gf_reference%kubo_term_two))
      write (*, '(a,es12.4,a,es12.4,a,es12.4)') '  '//trim(label)// &
         ' GF ref full/term1/term2 max|dchi| = ', reference_residual, ' ', reference_term_one_residual, ' ', &
         reference_term_two_residual
      if (trim(direct_gf%implementation) /= 'eigenbasis' .or. trim(direct_gf_reference%implementation) /= 'dense-reference' .or. &
          direct_gf%resolvent_calls /= 0 .or. direct_gf_reference%resolvent_calls <= 0 .or. &
          reference_residual > 2.0e-10_rp .or. reference_term_one_residual > 2.0e-10_rp .or. &
          reference_term_two_residual > 2.0e-10_rp) then
         failed = .true.
      end if
      diagnostic_residual = maxval(abs(direct_gf%kubo_term_one + direct_gf%kubo_term_two - direct_gf%susceptibility))
      if (.not. allocated(direct_gf%kubo_term_one) .or. .not. allocated(direct_gf%kubo_term_two) .or. &
          .not. allocated(direct_gf%k_susceptibility) .or. diagnostic_residual > 2.0e-12_rp .or. &
          maxval(abs(sum(direct_gf%k_susceptibility, dim=4) - direct_gf%susceptibility)) > 2.0e-12_rp .or. &
          .not. ieee_is_finite(direct_gf%left_spectral_zeroth_residual) .or. &
          .not. ieee_is_finite(direct_gf%left_spectral_first_residual) .or. &
          .not. ieee_is_finite(direct_gf%left_spectral_fermi_residual) .or. &
          direct_gf%spacing_over_integration_eta <= 0.0_rp) then
         failed = .true.
         write (*, '(a,es12.4)') '  '//trim(label)//' GF diagnostic term residual = ', diagnostic_residual
      end if
      gf_residual = maxval(abs(direct_gf%susceptibility - oracle_gf))
      write (*, '(a,es12.4)') '  '//trim(label)//' GF product oracle max|dchi| = ', gf_residual
      if (gf_residual > 2.0e-10_rp .or. direct_gf%point_response_allocated) failed = .true.
      finite_width_residual = maxval(abs(direct_gf%susceptibility - finite_width_oracle%susceptibility))
      write (*, '(a,es12.4)') '  '//trim(label)//' GF finite-width oracle max|dchi| = ', finite_width_residual
      if (finite_width_residual > 2.0e-10_rp .or. trim(finite_width_oracle%implementation) /= 'direct-oracle') then
         failed = .true.
      end if

      if (do_ladder) then
         request%integration_points = 201
         request%integration_eta = 0.018_rp
         call evaluate_projected_gf_chi0(request, gf_coarse)
         request%integration_points = 401
         request%integration_eta = 0.009_rp
         call evaluate_projected_gf_chi0(request, gf_fine)
         coarse_difference = maxval(abs(gf_coarse%susceptibility - direct_lehmann%susceptibility))
         fine_difference = maxval(abs(gf_fine%susceptibility - direct_lehmann%susceptibility))
         trend_residual = fine_difference - coarse_difference
         write (*, '(a,es12.4,a,es12.4,a,es12.4)') '  '//trim(label)// &
            ' GF->Lehmann coarse/fine/trend = ', coarse_difference, ' ', fine_difference, ' ', trend_residual
         if (.not. ieee_is_finite(fine_difference) .or. fine_difference > coarse_difference*1.05_rp + 2.0e-8_rp) then
            failed = .true.
         end if
      end if

      ! Retarded sign check: the diagonal bare loss is nonnegative under the
      ! locked convention, hence Im chi_plus(omega>0) is nonpositive.
      if (frequencies(size(frequencies)) > 0.0_rp) then
         if (maxval([(aimag(direct_lehmann%susceptibility(i, i, size(frequencies))), i=1, contract%nsite)]) > &
             2.0e-10_rp) then
            failed = .true.
            write (*, '(a)') '  '//trim(label)//' retarded pole-sign check: FAIL'
         else
            write (*, '(a)') '  '//trim(label)//' retarded pole-sign check: PASS'
         end if
      end if
      deallocate(oracle_lehmann, oracle_gf, product_request%selected_l, product_gf_request%selected_l)
   end subroutine run_projection_case

   subroutine run_covariance_case(contract, product_plus, product_minus, state, endpoint, q, frequency, eta_in, margin, label, &
                                  failed, negative_endpoint)
      type(projected_site_spin_contract), intent(in), target :: contract
      type(lmto_product_response_basis), intent(in), target :: product_plus, product_minus
      type(lr_electronic_state), intent(in), target :: state, endpoint
      real(rp), intent(in) :: q(3), frequency, eta_in, margin
      character(len=*), intent(in) :: label
      logical, intent(inout) :: failed
      type(lr_electronic_state), intent(in), optional, target :: negative_endpoint
      type(projected_chi0_request) :: plus_request, minus_request
      type(projected_chi0_result) :: plus_result, minus_result, plus_gf, minus_gf
      type(projected_site_spin_contract), pointer :: contract_minus
      type(lmto_product_response_basis), pointer :: product_minus_ptr
      type(lr_electronic_state), pointer :: endpoint_minus
      real(rp) :: residual, gf_residual

      contract_minus => contract
      product_minus_ptr => product_minus
      endpoint_minus => endpoint
      if (present(negative_endpoint)) endpoint_minus => negative_endpoint
      call prepare_request(plus_request, contract, product_plus, state, endpoint, q, &
         [frequency], eta_in, margin)
      plus_request%channel = lr_channel_plus
      call evaluate_projected_lehmann_chi0(plus_request, plus_result)
      call prepare_request(minus_request, contract_minus, product_minus_ptr, state, endpoint_minus, &
         -plus_request%q, [-frequency], eta_in, margin)
      minus_request%channel = lr_channel_minus
      call evaluate_projected_lehmann_chi0(minus_request, minus_result)
      residual = maxval(abs(plus_result%susceptibility - conjg(minus_result%susceptibility)))
      write (*, '(a,es12.4)') '  '//trim(label)//' q/-q Lehmann covariance residual = ', residual
      if (residual > 2.0e-8_rp) failed = .true.

      plus_request%integration_points = 101
      plus_request%integration_eta = 0.012_rp
      call evaluate_projected_gf_chi0(plus_request, plus_gf)
      minus_request%integration_points = 101
      minus_request%integration_eta = 0.012_rp
      call evaluate_projected_gf_chi0(minus_request, minus_gf)
      gf_residual = maxval(abs(plus_gf%susceptibility - conjg(minus_gf%susceptibility)))
      write (*, '(a,es12.4)') '  '//trim(label)//' q/-q GF covariance residual = ', gf_residual
      if (gf_residual > 2.0e-8_rp) failed = .true.
   end subroutine run_covariance_case

   subroutine prepare_request(request, contract, product, state, endpoint, q, frequencies, eta_in, margin)
      type(projected_chi0_request), intent(out) :: request
      type(projected_site_spin_contract), intent(in), target :: contract
      type(lmto_product_response_basis), intent(in), target :: product
      type(lr_electronic_state), intent(in), target :: state, endpoint
      real(rp), intent(in) :: q(3), frequencies(:), eta_in, margin

      request%q = q
      request%frequencies = frequencies
      request%eta = eta_in
      request%channel = lr_channel_plus
      request%energy_margin = margin
      request%contract => contract
      request%product_basis => product
      request%electronic_state => state
      request%q_endpoint_state => endpoint
   end subroutine prepare_request

   subroutine project_product_matrix(contract, product, compact, projected)
      type(projected_site_spin_contract), intent(in) :: contract
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: compact(:, :, :)
      complex(rp), intent(out) :: projected(:, :, :)
      complex(rp), allocatable :: functionals(:, :), image(:)
      integer :: ifrequency, i, j

      if (any(shape(projected) /= [contract%nsite, contract%nsite, size(compact, 3)])) then
         error stop 'project_product_matrix: output shape mismatch'
      end if
      allocate(functionals(product%product_dimension, contract%nsite), image(product%product_dimension))
      call contract%site_integration_functional(product, functionals)
      do ifrequency = 1, size(compact, 3)
         do j = 1, contract%nsite
            image = matmul(compact(:, :, ifrequency), functionals(:, j))
            do i = 1, contract%nsite
               projected(i, j, ifrequency) = dot_product(functionals(:, i), image)
            end do
         end do
      end do
      deallocate(functionals, image)
   end subroutine project_product_matrix

   subroutine build_mesh(mesh, mesh_a, mesh_b)
      real(rp), intent(out) :: mesh(:)
      real(rp), intent(in) :: mesh_a, mesh_b
      integer :: ir

      mesh(1) = 0.0_rp
      do ir = 2, size(mesh)
         mesh(ir) = mesh_b*(exp(mesh_a*real(ir - 1, rp)) - 1.0_rp)
      end do
   end subroutine build_mesh

   subroutine setup_radial_bases(bases, mesh, mesh_a, mesh_b)
      type(lmto_radial_basis), intent(out) :: bases(:)
      real(rp), intent(in) :: mesh(:), mesh_a, mesh_b
      integer :: site, spin, l, ir
      real(rp) :: radial_scale

      do site = 1, size(bases)
         call bases(site)%initialize(size(mesh), 2, 2)
         bases(site)%rofi = mesh
         bases(site)%mesh_a = mesh_a
         bases(site)%mesh_b = mesh_b
         bases(site)%channel_present = .true.
         bases(site)%gfac = 1.0_rp
         do spin = 1, 2
            do l = 0, 2
               bases(site)%enu_work(l + 1, spin) = -0.28_rp + 0.035_rp*real(l, rp)
               radial_scale = 0.75_rp + 0.025_rp*real(l, rp)
               do ir = 1, size(mesh)
                  bases(site)%phi_large(ir, l + 1, spin) = (0.65_rp + 0.04_rp*real(l, rp))* &
                     (mesh(ir)**l)*exp(-radial_scale*mesh(ir))
                  bases(site)%phidot_large(ir, l + 1, spin) = (0.10_rp + 0.02_rp*real(l, rp))* &
                     (mesh(ir)**l)*(1.0_rp + 0.13_rp*mesh(ir))*exp(-0.6_rp*radial_scale*mesh(ir))
                  bases(site)%phi_small(ir, l + 1, spin) = 0.0_rp
                  bases(site)%phidot_small(ir, l + 1, spin) = 0.0_rp
                  bases(site)%phiddot_large(ir, l + 1, spin) = (0.018_rp + 0.004_rp*real(l + spin, rp))* &
                     (mesh(ir)**l)*(1.0_rp + 0.09_rp*mesh(ir))*exp(-0.45_rp*radial_scale*mesh(ir))
                  bases(site)%phiddot_small(ir, l + 1, spin) = 0.0_rp
               end do
            end do
         end do
      end do
   end subroutine setup_radial_bases

   subroutine build_complex_spin_symmetric_eigenpairs(nsite, nbands, values, vectors)
      integer, intent(in) :: nsite, nbands
      real(rp), intent(out) :: values(:)
      complex(rp), intent(out) :: vectors(:, :)
      integer :: norb, one_spin_basis, ib, band, spin, dof, i, j, lwork, info
      complex(rp), allocatable :: hamiltonian(:, :), one_vectors(:, :), work(:)
      complex(rp) :: work_query(1)
      real(rp), allocatable :: one_values(:), rwork(:)
      external :: zheev

      norb = 9
      one_spin_basis = norb*nsite
      allocate(hamiltonian(one_spin_basis, one_spin_basis), one_vectors(one_spin_basis, one_spin_basis), &
         one_values(one_spin_basis), rwork(max(1, 3*one_spin_basis - 2)))
      hamiltonian = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, one_spin_basis
         hamiltonian(i, i) = cmplx(-0.65_rp + 0.35_rp*real(i - 1, rp), 0.0_rp, rp)
         do j = i + 1, one_spin_basis
            hamiltonian(i, j) = cmplx(0.0025_rp*real(mod(i + 2*j, 7) + 1, rp), &
               0.0017_rp*real(mod(j - i, 5) + 1, rp), rp)
            hamiltonian(j, i) = conjg(hamiltonian(i, j))
         end do
      end do
      one_vectors = hamiltonian
      call zheev('V', 'U', one_spin_basis, one_vectors, one_spin_basis, one_values, work_query, -1, rwork, info)
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      call zheev('V', 'U', one_spin_basis, one_vectors, one_spin_basis, one_values, work, lwork, rwork, info)
      if (info /= 0) error stop 'UnitLrProjectedReciprocalChi0: complex fixture diagonalization failed'

      vectors = cmplx(0.0_rp, 0.0_rp, rp)
      values = 0.0_rp
      do ib = 1, nbands/2
         values(2*ib - 1) = one_values(ib)
         values(2*ib) = one_values(ib)
         do dof = 1, one_spin_basis
            do spin = 1, 2
               i = site_spin_index(dof, spin, nsite)
               vectors(i, 2*ib - 1) = one_vectors(dof, ib)
               vectors(i, 2*ib) = one_vectors(dof, ib)
            end do
         end do
      end do
      deallocate(hamiltonian, one_vectors, one_values, work, rwork)
   end subroutine build_complex_spin_symmetric_eigenpairs

   pure integer function site_spin_index(dof, spin, nsite) result(index)
      integer, intent(in) :: dof, spin, nsite
      integer :: norb, site, orbital

      norb = 9
      site = (dof - 1)/norb + 1
      orbital = mod(dof - 1, norb) + 1
      index = (site - 1)*2*norb + (spin - 1)*norb + orbital
      if (site > nsite) error stop 'site_spin_index: invalid fixture index'
   end function site_spin_index

   subroutine build_occupations(values, occupations, temperature)
      real(rp), intent(in) :: values(:, :), temperature
      real(rp), intent(out) :: occupations(:, :)
      integer :: ib, ik

      do ik = 1, size(values, 2)
         do ib = 1, size(values, 1)
            occupations(ib, ik) = lr_fermi_dirac_occupation(values(ib, ik), 0.0_rp, temperature)
         end do
      end do
   end subroutine build_occupations

   subroutine apply_site_phase(input_vectors, nsite, q, output_vectors)
      complex(rp), intent(in) :: input_vectors(:, :)
      integer, intent(in) :: nsite
      real(rp), intent(in) :: q(3)
      complex(rp), intent(out) :: output_vectors(:, :)
      integer :: site, first, last, norb
      real(rp) :: angle
      complex(rp) :: phase

      norb = 9
      output_vectors = input_vectors
      do site = 1, nsite
         angle = 2.0_rp*3.1415926535897932384626433832795_rp*q(1)*real(site - 1, rp)
         phase = cmplx(cos(angle), sin(angle), rp)
         first = (site - 1)*2*norb + 1
         last = site*2*norb
         output_vectors(first:last, :) = phase*output_vectors(first:last, :)
      end do
   end subroutine apply_site_phase

end program test_lr_projected_reciprocal_chi0
