!------------------------------------------------------------------------------
! TDVK-03A/03B reciprocal-GF real-axis quadrature closure audit.
!
! This executable is deliberately an audit driver rather than a pass/fail
! physics test.  It reuses the TDVK-02R3 fixture, calls the existing compact
! reciprocal-GF evaluator, and records spectral moments and GF/Lehmann
! differences for the requested numerical ladders.
!------------------------------------------------------------------------------
program test_lr_product_gf_quadrature_audit
   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use response_angular_basis_mod, only: response_angular_pi
   use lr_response_space_mod, only: response_space_layout
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, &
      lmto_product_channel_plus, lmto_product_channel_minus
   use lr_pauli_transition_vertex_mod, only: pauli_endpoint_state
   use lr_ks_susceptibility_mod, only: lr_electronic_state, &
      lr_product_ks_susceptibility_request, lr_product_ks_susceptibility_result, &
      lr_channel_plus, lr_channel_minus, lr_fermi_dirac_occupation, &
      evaluate_lr_product_ks_susceptibility
   use lr_gf_susceptibility_mod, only: build_weighted_resolvent
   use lr_lmto_endpoint_branches_mod, only: lmto_product_nbranch, lmto_product_max_gf_moment, &
      lmto_product_branch_powers, lmto_product_branch_label, lmto_product_energy_power
   use lr_product_gf_susceptibility_mod, only: lr_product_gf_susceptibility_request, &
      lr_product_gf_susceptibility_result, evaluate_lr_product_gf_susceptibility, &
      build_lr_product_gf_transition_amplitudes, lr_product_gf_contraction_scalar, &
      lr_product_gf_contraction_optimized, lr_product_gf_contraction_factorized
   implicit none

   integer, parameter :: nr = 7, orbital_lmax = 1, response_lmax = 2, nsite = 1
   integer, parameter :: norb = (orbital_lmax + 1)**2
   integer, parameter :: nbasis = 2*norb*nsite, nbands = nbasis, nk = 2
   integer, parameter :: nfixed = 7, neta = 4, nmargin = 3
   integer, parameter :: nmixed = 3
   real(rp), parameter :: mesh_a = 0.03_rp, mesh_b = 0.10_rp
   real(rp), parameter :: eta = 0.04_rp, temperature = 0.03_rp
   real(rp), parameter :: q_gamma(3) = [0.0_rp, 0.0_rp, 0.0_rp]
   real(rp), parameter :: q_finite(3) = [0.23_rp, 0.0_rp, 0.0_rp]
   real(rp), parameter :: fixed_eta = 0.001_rp
   integer, parameter :: fixed_points(nfixed) = [101, 201, 401, 801, 1601, 3201, 6401]
   real(rp), parameter :: eta_values(neta) = [0.010_rp, 0.005_rp, 0.0025_rp, 0.001_rp]
   integer, parameter :: eta_points(neta) = [515, 1029, 2055, 6401]
   real(rp), parameter :: margin_values(nmargin) = [0.60_rp, 1.00_rp, 2.00_rp]
   integer, parameter :: margin_points(nmargin) = [6401, 8401, 13401]
   real(rp), parameter :: mixed_eta_values(nmixed) = [0.010_rp, 0.005_rp, 0.0025_rp]
   integer, parameter :: mixed_eta_points(nmixed) = [801, 1601, 3201]
   real(rp), parameter :: gamma_frequencies(2) = [0.0_rp, 0.17_rp]
   real(rp), parameter :: static_frequency(1) = [0.0_rp]

   real(rp) :: radius(nr), eigenvalues(nbands, nk), k_points(3, nk)
   real(rp) :: k_weights(nk), occupations(nbands, nk)
   complex(rp) :: eigenvectors(nbasis, nbands, nk)
   complex(rp) :: hamiltonians(nbasis, nbasis, nk)
   type(lmto_radial_basis), target :: radial(nsite)
   type(response_space_layout), target :: space
   type(lmto_product_response_basis), target :: product_plus, product_minus
   type(lr_electronic_state), target :: state, endpoint
   character(len=32) :: mode
   character(len=32) :: argument
   integer :: selected_index

   call get_command_argument(1, mode)
   if (len_trim(mode) == 0) mode = 'full'

   call build_mesh(radius)
   call setup_radial_basis(radial, radius)
   call space%initialize(nsite, response_lmax, radius, mesh_a, mesh_b, 1)
   if (index(trim(mode), 'mixed') == 1) then
      call build_mixed_electronic_fixture(eigenvalues, eigenvectors, hamiltonians, k_points, k_weights, occupations)
   else
      call build_electronic_fixture(eigenvalues, eigenvectors, hamiltonians, k_points, k_weights, occupations)
   end if
   call state%initialize(eigenvalues, eigenvectors, k_points, k_weights, occupations, 0.0_rp, temperature)
   call endpoint%initialize(eigenvalues, eigenvectors, k_points, k_weights, occupations, 0.0_rp, temperature)
   call product_plus%initialize(space, radial, lmto_product_channel_plus, .false.)
   call product_minus%initialize(space, radial, lmto_product_channel_minus, .false.)

   if (index(trim(mode), 'mixed') == 1) then
      write (*, '(a)') 'TDVK-03B mixed-eigenvector reciprocal-GF closure audit'
      call report_dresp10a_contract(product_plus)
      call report_transition_factorization_oracle(product_plus, 'chi_plus')
      call report_transition_factorization_oracle(product_minus, 'chi_minus')
      call report_weighted_resolvent_oracle()
   else
      write (*, '(a)') 'TDVK-03A reciprocal-GF quadrature audit'
   end if
   write (*, '(a,i0,a,i0,a,i0,a,i0)') 'fixture nbasis=', nbasis, ' nbands=', nbands, ' nk=', nk, &
      ' product_dim=', product_plus%product_dimension
   write (*, '(a,es16.8,a,es16.8,a,es16.8)') 'controls eta=', eta, ' temperature=', temperature, &
      ' fixed_integration_eta=', fixed_eta

   select case (trim(mode))
   case ('full')
      call run_fixed_ladder()
      call run_fixed_q_sample()
      call run_eta_ladder()
      call run_margin_ladder()
   case ('extra')
      call run_one_fixed_sample(12801)
   case ('fixed_one')
      call get_command_argument(2, argument)
      read (argument, *) selected_index
      call run_one_fixed_sample(selected_index)
   case ('eta_one')
      call get_command_argument(2, argument)
      read (argument, *) selected_index
      if (selected_index < 1 .or. selected_index > neta) error stop 'invalid eta ladder index'
      call report_moments('eta_one', state, eta_points(selected_index), eta_values(selected_index), 0.60_rp)
      call report_response('eta_one', product_plus, lr_channel_plus, q_gamma, gamma_frequencies, &
         eta_points(selected_index), eta_values(selected_index), 0.60_rp)
   case ('margin_one')
      call get_command_argument(2, argument)
      read (argument, *) selected_index
      if (selected_index < 1 .or. selected_index > nmargin) error stop 'invalid margin ladder index'
      call report_moments('margin_one', state, margin_points(selected_index), fixed_eta, &
         margin_values(selected_index))
      call report_response('margin_one', product_plus, lr_channel_plus, q_gamma, gamma_frequencies, &
         margin_points(selected_index), fixed_eta, margin_values(selected_index))
   case ('finite_q')
      call report_moments('fixed_q', state, 6401, fixed_eta, 0.60_rp)
      call report_response('fixed_q', product_minus, lr_channel_minus, q_finite, static_frequency, &
         6401, fixed_eta, 0.60_rp)
   case ('mixed')
      call run_mixed_ladder()
      call run_mixed_finite_q()
      call run_rotation_oracle()
   case ('mixed_one')
      call get_command_argument(2, argument)
      read (argument, *) selected_index
      if (selected_index < 1 .or. selected_index > nmixed) error stop 'invalid mixed eta ladder index'
      call require_resolved_mesh(mixed_eta_points(selected_index), mixed_eta_values(selected_index), 1.0_rp, state)
      call report_mixed_fixture()
      call report_moments('mixed_one', state, mixed_eta_points(selected_index), mixed_eta_values(selected_index), 1.0_rp)
      call report_response('mixed_one', product_plus, lr_channel_plus, q_gamma, gamma_frequencies, &
         mixed_eta_points(selected_index), mixed_eta_values(selected_index), 1.0_rp, .true.)
   case ('mixed_q')
      call report_mixed_fixture()
      call require_resolved_mesh(mixed_eta_points(nmixed), mixed_eta_values(nmixed), 1.0_rp, state)
      call report_moments('mixed_q', state, mixed_eta_points(nmixed), mixed_eta_values(nmixed), 1.0_rp)
      call report_response('mixed_q', product_minus, lr_channel_minus, q_finite, static_frequency, &
         mixed_eta_points(nmixed), mixed_eta_values(nmixed), 1.0_rp, .true.)
   case ('mixed_rotation')
      call report_mixed_fixture()
      call run_rotation_oracle()
   case default
      error stop 'test_lr_product_gf_quadrature_audit: invalid audit mode'
   end select

   if (index(trim(mode), 'mixed') == 1) then
      write (*, '(a)') 'TDVK-03B audit execution complete'
   else
      write (*, '(a)') 'TDVK-03A audit execution complete'
   end if

contains

   subroutine build_mesh(mesh)
      real(rp), intent(out) :: mesh(:)
      integer :: ir

      mesh(1) = 0.0_rp
      do ir = 2, size(mesh)
         mesh(ir) = mesh_b*(exp(mesh_a*real(ir - 1, rp)) - 1.0_rp)
      end do
   end subroutine build_mesh

   subroutine setup_radial_basis(bases, mesh)
      type(lmto_radial_basis), intent(out) :: bases(:)
      real(rp), intent(in) :: mesh(:)
      integer :: isite, ispin, l, ir

      do isite = 1, size(bases)
         call bases(isite)%initialize(size(mesh), orbital_lmax, 2)
         bases(isite)%rofi = mesh
         bases(isite)%mesh_a = mesh_a
         bases(isite)%mesh_b = mesh_b
         bases(isite)%channel_present = .true.
         bases(isite)%enu_work = 0.0_rp
         bases(isite)%gfac = 1.0_rp
         do ispin = 1, 2
            do l = 0, orbital_lmax
               do ir = 1, size(mesh)
                  bases(isite)%phi_large(ir, l + 1, ispin) = &
                     (1.0_rp + 0.17_rp*real(l, rp) + 0.06_rp*real(ispin, rp))* &
                     (1.0_rp + 0.31_rp*mesh(ir))
                  bases(isite)%phidot_large(ir, l + 1, ispin) = &
                     (0.11_rp + 0.03_rp*real(l, rp) + 0.02_rp*real(ispin, rp))* &
                     (1.0_rp + 0.19_rp*mesh(ir)**2)
                  bases(isite)%phiddot_large(ir, l + 1, ispin) = &
                     (0.025_rp + 0.007_rp*real(l + ispin, rp))* &
                     (1.0_rp + 0.21_rp*mesh(ir) + 0.09_rp*mesh(ir)**2)
               end do
            end do
         end do
      end do
   end subroutine setup_radial_basis

   subroutine build_electronic_fixture(values, vectors, h, points, weights, occ)
      real(rp), intent(out) :: values(:, :), points(:, :), weights(:), occ(:, :)
      complex(rp), intent(out) :: vectors(:, :, :)
      complex(rp), intent(out) :: h(:, :, :)
      integer :: ib, ik

      values = 0.0_rp
      do ik = 1, nk
         do ib = 1, nbands
            values(ib, ik) = -0.72_rp + 0.19_rp*real(ib - 1, rp) + 0.037_rp*real(ik - 1, rp)
         end do
      end do
      h = cmplx(0.0_rp, 0.0_rp, rp)
      do ik = 1, nk
         do ib = 1, nbands
            h(ib, ib, ik) = cmplx(values(ib, ik), 0.0_rp, rp)
         end do
      end do
      vectors = cmplx(0.0_rp, 0.0_rp, rp)
      do ik = 1, nk
         do ib = 1, nbands
            vectors(ib, ib, ik) = cmplx(1.0_rp, 0.0_rp, rp)
         end do
      end do
      points = 0.0_rp
      weights = [0.35_rp, 0.65_rp]
      do ik = 1, nk
         do ib = 1, nbands
            occ(ib, ik) = lr_fermi_dirac_occupation(values(ib, ik), 0.0_rp, temperature)
         end do
      end do
   end subroutine build_electronic_fixture

   subroutine build_mixed_electronic_fixture(values, vectors, h, points, weights, occ)
      real(rp), intent(out) :: values(:, :), points(:, :), weights(:), occ(:, :)
      complex(rp), intent(out) :: vectors(:, :, :), h(:, :, :)
      complex(rp) :: h_work(nbasis, nbasis)
      integer :: ib, jb, ik

      ! Dense Hermitian matrices give a generic complex orbital mixture while
      ! retaining the same dimensions and energy scale as the TDVK-02R3
      ! fixture.  Every off-diagonal pair is populated; this is not a hidden
      ! spin- or orbital-block test.
      h = cmplx(0.0_rp, 0.0_rp, rp)
      do ik = 1, nk
         h_work = cmplx(0.0_rp, 0.0_rp, rp)
         do ib = 1, nbasis
            h_work(ib, ib) = cmplx(-0.72_rp + 0.19_rp*real(ib - 1, rp) + &
               0.037_rp*real(ik - 1, rp), 0.0_rp, rp)
            do jb = ib + 1, nbasis
               h_work(ib, jb) = cmplx(0.0035_rp*real(ib + 2*jb, rp) + 0.0002_rp*real(ik, rp), &
                  0.0021_rp*real(jb - ib, rp) + 0.0001_rp*real(ib + jb + ik, rp), rp)
               h_work(jb, ib) = conjg(h_work(ib, jb))
            end do
         end do
         h(:, :, ik) = h_work
         call hermitian_eig(h_work, nbasis, values(:, ik), vectors(:, :, ik))
      end do
      points = 0.0_rp
      weights = [0.35_rp, 0.65_rp]
      do ik = 1, nk
         do ib = 1, nbands
            occ(ib, ik) = lr_fermi_dirac_occupation(values(ib, ik), 0.0_rp, temperature)
         end do
      end do
   end subroutine build_mixed_electronic_fixture

   subroutine report_dresp10a_contract(product)
      type(lmto_product_response_basis), intent(in) :: product
      character(len=32) :: labels
      integer :: branch, p, q

      labels = ''
      do branch = 1, lmto_product_nbranch
         if (branch > 1) labels = trim(labels)//','
         labels = trim(labels)//trim(lmto_product_branch_label(branch))
         call lmto_product_branch_powers(branch, p, q)
         if (p < 0 .or. q < 0 .or. p > 2 .or. q > 2) then
            error stop 'test_lr_product_gf_quadrature_audit: invalid DRESP-10A branch map'
         end if
      end do
      write (*, '(a,1x,a,i0,1x,a,i0,1x,a,a,1x,a,i0)') 'DRESP10A_CONTRACT', &
         'product_dimension=', product%product_dimension, 'branch_count=', product%product_endpoint_branches, &
         'labels=', trim(labels), 'maximum_gf_moment=', product%maximum_gf_energy_moment
      if (product%product_endpoint_branches /= 6 .or. trim(labels) /= '00,10,01,11,20,02' .or. &
          product%maximum_gf_energy_moment /= 4) then
         error stop 'test_lr_product_gf_quadrature_audit: DRESP-10A contract failed'
      end if
   end subroutine report_dresp10a_contract

   subroutine report_transition_factorization_oracle(product, channel)
      type(lmto_product_response_basis), intent(in) :: product
      character(len=*), intent(in) :: channel
      complex(rp), allocatable :: vertices(:, :, :, :), transitions(:, :, :), reference(:)
      type(pauli_endpoint_state) :: left_band, right_band
      real(rp) :: absolute_error, relative_error, maximum_absolute, maximum_relative, scale
      integer, parameter :: selected_left(3) = [1, 3, 5], selected_right(3) = [2, 8, 7]
      integer :: pair

      call product%component_vertex_tensor(vertices)
      allocate(transitions(nbands, nbands, product%product_dimension), reference(product%product_dimension))
      call build_lr_product_gf_transition_amplitudes(vertices, state, 1, state, 1, transitions)
      maximum_absolute = 0.0_rp
      maximum_relative = 0.0_rp
      do pair = 1, size(selected_left)
         call left_band%initialize(state%eigenvalues(selected_left(pair), 1), &
            state%eigenvectors(:, selected_left(pair), 1))
         call right_band%initialize(state%eigenvalues(selected_right(pair), 1), &
            state%eigenvectors(:, selected_right(pair), 1))
         call product%transition_coordinates(left_band, right_band, reference)
         scale = max(sqrt(sum(abs(reference)**2)), epsilon(1.0_rp))
         absolute_error = sqrt(sum(abs(transitions(selected_left(pair), selected_right(pair), :) - reference)**2))
         relative_error = absolute_error/scale
         maximum_absolute = max(maximum_absolute, absolute_error)
         maximum_relative = max(maximum_relative, relative_error)
      end do
      if (maximum_relative >= 1.0e-10_rp) then
         error stop 'test_lr_product_gf_quadrature_audit: GF transition factorization oracle failed'
      end if
      write (*, '(a,1x,a,1x,a,es16.8,1x,a,es16.8)') 'TRANSITION_FACTORIZATION', 'channel='//trim(channel), &
         'max_abs=', maximum_absolute, 'max_rel=', maximum_relative
      call report_branch_transition_activity(vertices, state, 1)
      deallocate(vertices, transitions, reference)
   end subroutine report_transition_factorization_oracle

   subroutine report_branch_transition_activity(vertices, electronic_state, ik)
      complex(rp), intent(in) :: vertices(:, :, :, :)
      type(lr_electronic_state), intent(in) :: electronic_state
      integer, intent(in) :: ik
      complex(rp), allocatable :: contribution(:)
      integer, parameter :: selected_left(3) = [1, 3, 5], selected_right(3) = [2, 8, 7]
      real(rp) :: maximum_activity
      integer :: branch, product_index, pair, left_band, right_band

      allocate(contribution(size(vertices, 4)))
      do branch = 1, lmto_product_nbranch
         maximum_activity = 0.0_rp
         do pair = 1, size(selected_left)
            left_band = selected_left(pair)
            right_band = selected_right(pair)
            contribution = cmplx(0.0_rp, 0.0_rp, rp)
            do product_index = 1, size(vertices, 4)
               contribution(product_index) = sum(conjg(electronic_state%eigenvectors(:, left_band, ik))* &
                  matmul(vertices(:, :, branch, product_index), electronic_state%eigenvectors(:, right_band, ik))) * &
                  lmto_product_energy_power(electronic_state%eigenvalues(left_band, ik), &
                     lmto_branch_left_power(branch))*lmto_product_energy_power( &
                     electronic_state%eigenvalues(right_band, ik), lmto_branch_right_power(branch))
            end do
            maximum_activity = max(maximum_activity, sqrt(sum(abs(contribution)**2)))
         end do
         write (*, '(a,1x,a,a,1x,a,es16.8)') 'BRANCH_TRANSITION', 'branch=', &
            trim(lmto_product_branch_label(branch)), 'max_norm=', maximum_activity
         if ((branch == 5 .or. branch == 6) .and. maximum_activity <= 1.0e-13_rp) then
            error stop 'test_lr_product_gf_quadrature_audit: second-order branch transition vanished'
         end if
      end do
      deallocate(contribution)
   end subroutine report_branch_transition_activity

   integer function lmto_branch_left_power(branch) result(power)
      integer, intent(in) :: branch
      integer :: unused_right

      call lmto_product_branch_powers(branch, power, unused_right)
   end function lmto_branch_left_power

   integer function lmto_branch_right_power(branch) result(power)
      integer, intent(in) :: branch
      integer :: unused_left

      call lmto_product_branch_powers(branch, unused_left, power)
   end function lmto_branch_right_power

   subroutine report_weighted_resolvent_oracle()
      complex(rp), parameter :: probe_z(3) = [cmplx(-0.41_rp, 0.07_rp, rp), &
         cmplx(0.13_rp, 0.11_rp, rp), cmplx(0.79_rp, 0.05_rp, rp)]
      complex(rp), allocatable :: actual(:, :, :), reference(:, :, :)
      real(rp) :: absolute_residual, relative_residual, maximum_absolute, maximum_relative
      integer :: ik, iz, p

      allocate(actual(nbasis, nbasis, lmto_product_max_gf_moment + 1), &
         reference(nbasis, nbasis, lmto_product_max_gf_moment + 1))
      do p = 0, lmto_product_max_gf_moment
         maximum_absolute = 0.0_rp
         maximum_relative = 0.0_rp
         do ik = 1, nk
            do iz = 1, size(probe_z)
               call build_weighted_resolvent(state, ik, probe_z(iz), actual)
               call exact_weighted_resolvent(state, ik, probe_z(iz), reference)
               absolute_residual = sqrt(sum(abs(actual(:, :, p + 1) - reference(:, :, p + 1))**2))
               relative_residual = absolute_residual/max(sqrt(sum(abs(reference(:, :, p + 1))**2)), &
                  epsilon(1.0_rp))
               maximum_absolute = max(maximum_absolute, absolute_residual)
               maximum_relative = max(maximum_relative, relative_residual)
            end do
         end do
         write (*, '(a,1x,a,i0,1x,a,es16.8,1x,a,es16.8)') 'WEIGHTED_RESOLVENT', 'p=', p, &
            'max_abs=', maximum_absolute, 'max_rel=', maximum_relative
         if (maximum_relative >= 1.0e-12_rp) then
            error stop 'test_lr_product_gf_quadrature_audit: weighted-resolvent oracle failed'
         end if
      end do
      deallocate(actual, reference)
   end subroutine report_weighted_resolvent_oracle

   subroutine exact_weighted_resolvent(electronic_state, ik, z, weighted_green)
      type(lr_electronic_state), intent(in) :: electronic_state
      integer, intent(in) :: ik
      complex(rp), intent(in) :: z
      complex(rp), intent(out) :: weighted_green(:, :, :)
      complex(rp) :: factor
      integer :: p, ib, i, j

      weighted_green = cmplx(0.0_rp, 0.0_rp, rp)
      do p = 0, size(weighted_green, 3) - 1
         do ib = 1, electronic_state%nbands
            factor = cmplx(lmto_product_energy_power(electronic_state%eigenvalues(ib, ik), p), 0.0_rp, rp)/ &
               (z - electronic_state%eigenvalues(ib, ik))
            do j = 1, electronic_state%nbasis
               do i = 1, electronic_state%nbasis
                  weighted_green(i, j, p + 1) = weighted_green(i, j, p + 1) + factor* &
                     electronic_state%eigenvectors(i, ib, ik)*conjg(electronic_state%eigenvectors(j, ib, ik))
               end do
            end do
         end do
      end do
   end subroutine exact_weighted_resolvent

   subroutine report_mixed_fixture()
      complex(rp), allocatable :: moment(:, :), moment_fermi(:, :), h_moment(:, :)
      real(rp) :: max_h_offdiag, min_h_offdiag, max_h_imag, min_h_imag
      real(rp) :: max_m_offdiag, max_n_offdiag, max_m_imag, max_n_imag
      real(rp) :: max_m_residual
      integer :: i, j, ik, p

      max_h_offdiag = 0.0_rp
      min_h_offdiag = huge(1.0_rp)
      max_h_imag = 0.0_rp
      min_h_imag = huge(1.0_rp)
      do ik = 1, nk
         do j = 1, nbasis
            do i = 1, nbasis
               if (i == j) cycle
               max_h_offdiag = max(max_h_offdiag, abs(hamiltonians(i, j, ik)))
               min_h_offdiag = min(min_h_offdiag, abs(hamiltonians(i, j, ik)))
               max_h_imag = max(max_h_imag, abs(aimag(hamiltonians(i, j, ik))))
               min_h_imag = min(min_h_imag, abs(aimag(hamiltonians(i, j, ik))))
            end do
         end do
      end do
      write (*, '(a,1x,a,es16.8,1x,a,es16.8)') 'MIXED_FIXTURE', 'emin=', minval(eigenvalues), &
         'emax=', maxval(eigenvalues)
      write (*, '(a,1x,a,es16.8,1x,a,es16.8)') 'MIXED_H_OFFDIAG', 'min=', min_h_offdiag, &
         'max=', max_h_offdiag
      write (*, '(a,1x,a,es16.8,1x,a,es16.8)') 'MIXED_H_IMAG_OFFDIAG', 'min=', min_h_imag, &
         'max=', max_h_imag

      allocate(moment(nbasis, nbasis), moment_fermi(nbasis, nbasis), h_moment(nbasis, nbasis))
      do p = 0, lmto_product_max_gf_moment
         max_m_offdiag = 0.0_rp
         max_n_offdiag = 0.0_rp
         max_m_imag = 0.0_rp
         max_n_imag = 0.0_rp
         max_m_residual = 0.0_rp
         do ik = 1, nk
            call exact_moment(state, ik, p, .false., moment)
            call exact_moment(state, ik, p, .true., moment_fermi)
            call matrix_power(hamiltonians(:, :, ik), p, h_moment)
            max_m_residual = max(max_m_residual, matrix_residual(moment, h_moment))
            do j = 1, nbasis
               do i = 1, nbasis
                  if (i == j) cycle
                  max_m_offdiag = max(max_m_offdiag, abs(moment(i, j)))
                  max_n_offdiag = max(max_n_offdiag, abs(moment_fermi(i, j)))
                  max_m_imag = max(max_m_imag, abs(aimag(moment(i, j))))
                  max_n_imag = max(max_n_imag, abs(aimag(moment_fermi(i, j))))
               end do
            end do
         end do
         write (*, '(a,1x,a,i0,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8)') &
            'MIXED_MOMENT_OFFDIAG', 'p=', p, 'maxM=', max_m_offdiag, 'maxN=', max_n_offdiag, &
            'maxImagM=', max_m_imag, 'maxImagN=', max_n_imag, 'rM_exact=', max_m_residual
      end do
      deallocate(moment, moment_fermi, h_moment)
   end subroutine report_mixed_fixture

   subroutine run_mixed_ladder()
      integer :: i
      character(len=16) :: tag

      call report_mixed_fixture()
      do i = 1, nmixed
         call require_resolved_mesh(mixed_eta_points(i), mixed_eta_values(i), 1.0_rp, state)
         write (tag, '(a,i0)') 'mixed_eta_', i
         call report_moments(trim(tag), state, mixed_eta_points(i), mixed_eta_values(i), 1.0_rp)
         call report_response(trim(tag), product_plus, lr_channel_plus, q_gamma, gamma_frequencies, &
            mixed_eta_points(i), mixed_eta_values(i), 1.0_rp, .true.)
      end do
   end subroutine run_mixed_ladder

   subroutine run_mixed_finite_q()
      call require_resolved_mesh(mixed_eta_points(nmixed), mixed_eta_values(nmixed), 1.0_rp, state)
      call report_moments('mixed_finite_q', state, mixed_eta_points(nmixed), mixed_eta_values(nmixed), 1.0_rp)
      call report_response('mixed_finite_q', product_minus, lr_channel_minus, q_finite, static_frequency, &
         mixed_eta_points(nmixed), mixed_eta_values(nmixed), 1.0_rp, .true.)
   end subroutine run_mixed_finite_q

   subroutine require_resolved_mesh(ne, integration_eta, margin, electronic_state)
      integer, intent(in) :: ne
      real(rp), intent(in) :: integration_eta, margin
      type(lr_electronic_state), intent(in) :: electronic_state
      real(rp) :: step, ratio

      step = (maxval(electronic_state%eigenvalues) - minval(electronic_state%eigenvalues) + 2.0_rp*margin)/ &
         real(ne - 1, rp)
      ratio = step/integration_eta
      if (ratio > 0.5_rp) error stop 'NUMERICAL_MESH_UNRESOLVED'
      write (*, '(a,1x,a,i0,1x,a,es16.8,1x,a,es16.8)') 'MIXED_RESOLUTION', 'n=', ne, 'h=', step, 'h_over_eta=', ratio
   end subroutine require_resolved_mesh

   subroutine run_fixed_ladder()
      integer :: i

      do i = 1, nfixed
         call run_one_fixed_sample(fixed_points(i))
      end do
   end subroutine run_fixed_ladder

   subroutine run_one_fixed_sample(ne)
      integer, intent(in) :: ne

      call report_moments('fixed', state, ne, fixed_eta, 0.60_rp)
      call report_response('fixed', product_plus, lr_channel_plus, q_gamma, gamma_frequencies, &
         ne, fixed_eta, 0.60_rp)
   end subroutine run_one_fixed_sample

   subroutine run_fixed_q_sample()
      integer, parameter :: ne = 6401

      call report_moments('fixed_q', state, ne, fixed_eta, 0.60_rp)
      call report_response('fixed_q', product_minus, lr_channel_minus, q_finite, static_frequency, &
         ne, fixed_eta, 0.60_rp)
   end subroutine run_fixed_q_sample

   subroutine run_eta_ladder()
      integer :: i
      character(len=16) :: tag

      do i = 1, neta
         write (tag, '(a,i0)') 'eta_', i
         call report_moments(trim(tag), state, eta_points(i), eta_values(i), 0.60_rp)
         call report_response(trim(tag), product_plus, lr_channel_plus, q_gamma, gamma_frequencies, &
            eta_points(i), eta_values(i), 0.60_rp)
      end do
   end subroutine run_eta_ladder

   subroutine run_margin_ladder()
      integer :: i
      character(len=16) :: tag

      do i = 1, nmargin
         write (tag, '(a,i0)') 'margin_', i
         call report_moments(trim(tag), state, margin_points(i), fixed_eta, margin_values(i))
         call report_response(trim(tag), product_plus, lr_channel_plus, q_gamma, gamma_frequencies, &
            margin_points(i), fixed_eta, margin_values(i))
      end do
   end subroutine run_margin_ladder

   subroutine report_moments(phase, electronic_state, ne, integration_eta, margin)
      character(len=*), intent(in) :: phase
      type(lr_electronic_state), intent(in) :: electronic_state
      integer, intent(in) :: ne
      real(rp), intent(in) :: integration_eta, margin

      complex(rp), allocatable :: green_r(:, :, :), green_a(:, :, :), spectral(:, :, :)
      complex(rp), allocatable :: numerical(:, :, :), numerical_fermi(:, :, :)
      complex(rp), allocatable :: exact(:, :, :), exact_fermi(:, :, :)
      real(rp) :: energy_min, energy_max, step, energy, fermi_weight, qw, ratio
      real(rp) :: residual_m, residual_n
      integer :: ik, ie, p

      energy_min = minval(electronic_state%eigenvalues) - margin
      energy_max = maxval(electronic_state%eigenvalues) + margin
      step = (energy_max - energy_min)/real(ne - 1, rp)
      ratio = step/integration_eta
      allocate(green_r(electronic_state%nbasis, electronic_state%nbasis, lmto_product_max_gf_moment + 1), &
         green_a(electronic_state%nbasis, electronic_state%nbasis, lmto_product_max_gf_moment + 1), &
         spectral(electronic_state%nbasis, electronic_state%nbasis, lmto_product_max_gf_moment + 1), &
         numerical(electronic_state%nbasis, electronic_state%nbasis, lmto_product_max_gf_moment + 1), &
         numerical_fermi(electronic_state%nbasis, electronic_state%nbasis, lmto_product_max_gf_moment + 1), &
         exact(electronic_state%nbasis, electronic_state%nbasis, lmto_product_max_gf_moment + 1), &
         exact_fermi(electronic_state%nbasis, electronic_state%nbasis, lmto_product_max_gf_moment + 1))

      do ik = 1, electronic_state%nk
         numerical = cmplx(0.0_rp, 0.0_rp, rp)
         numerical_fermi = cmplx(0.0_rp, 0.0_rp, rp)
         do ie = 1, ne
            energy = energy_min + real(ie - 1, rp)*step
            if (ie == 1 .or. ie == ne) then
               qw = 1.0_rp
            else if (mod(ie, 2) == 0) then
               qw = 4.0_rp
            else
               qw = 2.0_rp
            end if
            qw = qw*step/3.0_rp
            fermi_weight = lr_fermi_dirac_occupation(energy, electronic_state%fermi_level, &
               electronic_state%temperature)
            call build_weighted_resolvent(electronic_state, ik, cmplx(energy, integration_eta, rp), green_r)
            call build_weighted_resolvent(electronic_state, ik, cmplx(energy, -integration_eta, rp), green_a)
            spectral = cmplx(0.0_rp, 1.0_rp/(2.0_rp*response_angular_pi), rp)*(green_r - green_a)
            numerical = numerical + qw*spectral
            numerical_fermi = numerical_fermi + qw*fermi_weight*spectral
         end do
         do p = 0, lmto_product_max_gf_moment
            call exact_moment(electronic_state, ik, p, .false., exact(:, :, p + 1))
            call exact_moment(electronic_state, ik, p, .true., exact_fermi(:, :, p + 1))
            residual_m = matrix_residual(numerical(:, :, p + 1), exact(:, :, p + 1))
            residual_n = matrix_residual(numerical_fermi(:, :, p + 1), exact_fermi(:, :, p + 1))
            if (.not. ieee_is_finite(residual_m) .or. .not. ieee_is_finite(residual_n)) then
               error stop 'test_lr_product_gf_quadrature_audit: non-finite spectral moment residual'
            end if
            write (*, '(a,1x,a,1x,a,1x,i0,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8,1x,a,i0,1x,a,i0,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8)') &
               'MOMENT', 'phase='//trim(phase), 'ik=', ik, 'eta=', integration_eta, 'margin=', margin, &
               'h=', step, 'n=', ne, 'p=', p, 'ratio=', ratio, 'rM=', residual_m, 'rN=', residual_n
         end do
      end do
      deallocate(green_r, green_a, spectral, numerical, numerical_fermi, exact, exact_fermi)
   end subroutine report_moments

   subroutine exact_moment(electronic_state, ik, p, fermi_weighted, matrix)
      type(lr_electronic_state), intent(in) :: electronic_state
      integer, intent(in) :: ik, p
      logical, intent(in) :: fermi_weighted
      complex(rp), intent(out) :: matrix(:, :)
      integer :: ib, i, j
      real(rp) :: factor

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do ib = 1, electronic_state%nbands
         factor = electronic_state%eigenvalues(ib, ik)**p
         if (fermi_weighted) factor = factor*lr_fermi_dirac_occupation( &
            electronic_state%eigenvalues(ib, ik), electronic_state%fermi_level, electronic_state%temperature)
         do j = 1, electronic_state%nbasis
            do i = 1, electronic_state%nbasis
               matrix(i, j) = matrix(i, j) + factor*electronic_state%eigenvectors(i, ib, ik)* &
                  conjg(electronic_state%eigenvectors(j, ib, ik))
            end do
         end do
      end do
   end subroutine exact_moment

   subroutine matrix_power(base, power, result)
      complex(rp), intent(in) :: base(:, :)
      integer, intent(in) :: power
      complex(rp), intent(out) :: result(:, :)
      complex(rp), allocatable :: work(:, :)
      integer :: n

      allocate(work(size(base, 1), size(base, 2)))
      result = cmplx(0.0_rp, 0.0_rp, rp)
      do n = 1, size(base, 1)
         result(n, n) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
      work = base
      do n = 1, power
         result = matmul(result, work)
      end do
      deallocate(work)
   end subroutine matrix_power

   subroutine report_response(phase, product, channel, q, frequencies, ne, integration_eta, margin, compare_contractions)
      character(len=*), intent(in) :: phase, channel
      type(lmto_product_response_basis), target, intent(in) :: product
      real(rp), intent(in) :: q(3), frequencies(:), integration_eta, margin
      integer, intent(in) :: ne
      logical, intent(in), optional :: compare_contractions

      type(lr_product_gf_susceptibility_request) :: gf_request
      type(lr_product_gf_susceptibility_result) :: gf_result
      type(lr_product_gf_susceptibility_result) :: optimized_result, scalar_result
      type(lr_product_ks_susceptibility_request) :: lehmann_request
      type(lr_product_ks_susceptibility_result) :: lehmann_result
      type(lr_electronic_state), target :: local_endpoint
      integer :: ifrequency
      real(rp) :: d_f, r_f, d_inf, h, ratio, norm_reference
      logical :: compare

      compare = .false.
      if (present(compare_contractions)) compare = compare_contractions

      call local_endpoint%initialize(eigenvalues, eigenvectors, spread(q, 2, nk), k_weights, occupations, &
         0.0_rp, temperature)
      gf_request%q = q
      gf_request%frequencies = frequencies
      gf_request%eta = eta
      gf_request%channel = channel
      gf_request%integration_points = ne
      gf_request%integration_eta = integration_eta
      gf_request%energy_margin = margin
      gf_request%contraction_backend = lr_product_gf_contraction_factorized
      gf_request%product_basis => product
      gf_request%electronic_state => state
      gf_request%q_endpoint_state => local_endpoint
      call evaluate_lr_product_gf_susceptibility(gf_request, gf_result)

      if (compare) then
         gf_request%contraction_backend = lr_product_gf_contraction_optimized
         call evaluate_lr_product_gf_susceptibility(gf_request, optimized_result)
         gf_request%contraction_backend = lr_product_gf_contraction_scalar
         call evaluate_lr_product_gf_susceptibility(gf_request, scalar_result)
         do ifrequency = 1, size(frequencies)
            call report_three_way_contraction(phase, ifrequency, frequencies(ifrequency), scalar_result, optimized_result, gf_result)
         end do
      end if

      lehmann_request%q = q
      lehmann_request%frequencies = frequencies
      lehmann_request%eta = eta
      lehmann_request%channel = channel
      lehmann_request%product_basis => product
      lehmann_request%electronic_state => state
      lehmann_request%q_endpoint_state => local_endpoint
      call evaluate_lr_product_ks_susceptibility(lehmann_request, lehmann_result)

      h = (gf_result%energy_max - gf_result%energy_min)/real(ne - 1, rp)
      ratio = h/integration_eta
      do ifrequency = 1, size(frequencies)
         d_f = sqrt(sum(abs(gf_result%susceptibility(:, :, ifrequency) - &
            lehmann_result%susceptibility(:, :, ifrequency))**2))
         norm_reference = sqrt(sum(abs(lehmann_result%susceptibility(:, :, ifrequency))**2))
         r_f = d_f/max(norm_reference, epsilon(1.0_rp))
         d_inf = maxval(abs(gf_result%susceptibility(:, :, ifrequency) - &
            lehmann_result%susceptibility(:, :, ifrequency)))
         if (.not. ieee_is_finite(d_f) .or. .not. ieee_is_finite(r_f) .or. .not. ieee_is_finite(d_inf)) then
            error stop 'test_lr_product_gf_quadrature_audit: non-finite GF/Lehmann residual'
         end if
         write (*, '(a,1x,a,1x,a,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8,1x,a,i0,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8)') &
            'RESPONSE', 'phase='//trim(phase), 'channel='//trim(channel), &
            'frequency=', frequencies(ifrequency), 'eta_int=', integration_eta, 'margin=', margin, &
            'n=', ne, 'emin=', gf_result%energy_min, 'emax=', gf_result%energy_max, 'h=', h, &
            'ratio=', ratio, 'dF=', d_f, 'rF=', r_f, 'dInf=', d_inf, 'wall=', gf_result%wall_time_seconds
         if (compare) then
            call report_four_way_contraction(phase, ifrequency, frequencies(ifrequency), lehmann_result, &
               scalar_result, optimized_result, gf_result)
         end if
      end do
      deallocate(local_endpoint%k_points, local_endpoint%k_weights, local_endpoint%eigenvalues, &
         local_endpoint%eigenvectors, local_endpoint%occupations)
   end subroutine report_response

   subroutine report_three_way_contraction(phase, ifrequency, frequency, scalar_result, optimized_result, factorized_result)
      character(len=*), intent(in) :: phase
      integer, intent(in) :: ifrequency
      real(rp), intent(in) :: frequency
      type(lr_product_gf_susceptibility_result), intent(in) :: scalar_result, optimized_result, factorized_result
      real(rp) :: norm_scalar, norm_optimized, norm_factorized
      real(rp) :: scalar_optimized_dF, scalar_factorized_dF, optimized_factorized_dF
      real(rp) :: scalar_optimized_rF, scalar_factorized_rF, optimized_factorized_rF
      real(rp) :: scalar_optimized_dInf, scalar_factorized_dInf, optimized_factorized_dInf
      real(rp) :: scale

      norm_scalar = sqrt(sum(abs(scalar_result%susceptibility(:, :, ifrequency))**2))
      norm_optimized = sqrt(sum(abs(optimized_result%susceptibility(:, :, ifrequency))**2))
      norm_factorized = sqrt(sum(abs(factorized_result%susceptibility(:, :, ifrequency))**2))
      scalar_optimized_dF = matrix_difference(scalar_result%susceptibility(:, :, ifrequency), &
         optimized_result%susceptibility(:, :, ifrequency))
      scalar_factorized_dF = matrix_difference(scalar_result%susceptibility(:, :, ifrequency), &
         factorized_result%susceptibility(:, :, ifrequency))
      optimized_factorized_dF = matrix_difference(optimized_result%susceptibility(:, :, ifrequency), &
         factorized_result%susceptibility(:, :, ifrequency))
      scale = max(norm_scalar, norm_optimized, norm_factorized, epsilon(1.0_rp))
      scalar_optimized_rF = scalar_optimized_dF/scale
      scalar_factorized_rF = scalar_factorized_dF/scale
      optimized_factorized_rF = optimized_factorized_dF/scale
      scalar_optimized_dInf = maxval(abs(scalar_result%susceptibility(:, :, ifrequency) - &
         optimized_result%susceptibility(:, :, ifrequency)))
      scalar_factorized_dInf = maxval(abs(scalar_result%susceptibility(:, :, ifrequency) - &
         factorized_result%susceptibility(:, :, ifrequency)))
      optimized_factorized_dInf = maxval(abs(optimized_result%susceptibility(:, :, ifrequency) - &
         factorized_result%susceptibility(:, :, ifrequency)))
      if (max(scalar_optimized_rF, scalar_factorized_rF, optimized_factorized_rF) >= 1.0e-10_rp .or. &
          max(scalar_optimized_dInf, scalar_factorized_dInf, optimized_factorized_dInf) >= 1.0e-10_rp*max(1.0_rp, scale)) then
         error stop 'GF_CONTRACTION_BACKEND_MISMATCH'
      end if
      write (*, '(a,1x,a,1x,a,es16.8)') 'THREE_WAY', 'phase='//trim(phase), 'frequency=', frequency
      write (*, '(a,1x,a,3(1x,es16.8))') 'THREE_WAY', 'norms=', norm_scalar, norm_optimized, norm_factorized
      write (*, '(a,1x,a,3(1x,es16.8))') 'THREE_WAY', 'dF=', scalar_optimized_dF, scalar_factorized_dF, optimized_factorized_dF
      write (*, '(a,1x,a,3(1x,es16.8))') 'THREE_WAY', 'rF=', scalar_optimized_rF, scalar_factorized_rF, optimized_factorized_rF
      write (*, '(a,1x,a,3(1x,es16.8))') 'THREE_WAY', 'dInf=', scalar_optimized_dInf, scalar_factorized_dInf, optimized_factorized_dInf
      write (*, '(a,1x,a,3(1x,es16.8))') 'THREE_WAY', 'wall=', scalar_result%wall_time_seconds, &
         optimized_result%wall_time_seconds, factorized_result%wall_time_seconds
      write (*, '(a,1x,a,3(1x,es16.8))') 'THREE_WAY', 'speedup=', &
         scalar_result%wall_time_seconds/max(optimized_result%wall_time_seconds, epsilon(1.0_rp)), &
         scalar_result%wall_time_seconds/max(factorized_result%wall_time_seconds, epsilon(1.0_rp)), &
         optimized_result%wall_time_seconds/max(factorized_result%wall_time_seconds, epsilon(1.0_rp))
   end subroutine report_three_way_contraction

   subroutine report_four_way_contraction(phase, ifrequency, frequency, lehmann_result, scalar_result, &
                                          optimized_result, factorized_result)
      character(len=*), intent(in) :: phase
      integer, intent(in) :: ifrequency
      real(rp), intent(in) :: frequency
      type(lr_product_ks_susceptibility_result), intent(in) :: lehmann_result
      type(lr_product_gf_susceptibility_result), intent(in) :: scalar_result, optimized_result, factorized_result
      real(rp) :: norm_lehmann, norm_scalar, norm_optimized, norm_factorized
      real(rp) :: scalar_lehmann, optimized_lehmann, factorized_lehmann
      real(rp) :: scalar_optimized, scalar_factorized, optimized_factorized
      real(rp) :: scale

      norm_lehmann = sqrt(sum(abs(lehmann_result%susceptibility(:, :, ifrequency))**2))
      norm_scalar = sqrt(sum(abs(scalar_result%susceptibility(:, :, ifrequency))**2))
      norm_optimized = sqrt(sum(abs(optimized_result%susceptibility(:, :, ifrequency))**2))
      norm_factorized = sqrt(sum(abs(factorized_result%susceptibility(:, :, ifrequency))**2))
      scalar_lehmann = matrix_difference(scalar_result%susceptibility(:, :, ifrequency), &
         lehmann_result%susceptibility(:, :, ifrequency))/max(norm_lehmann, epsilon(1.0_rp))
      optimized_lehmann = matrix_difference(optimized_result%susceptibility(:, :, ifrequency), &
         lehmann_result%susceptibility(:, :, ifrequency))/max(norm_lehmann, epsilon(1.0_rp))
      factorized_lehmann = matrix_difference(factorized_result%susceptibility(:, :, ifrequency), &
         lehmann_result%susceptibility(:, :, ifrequency))/max(norm_lehmann, epsilon(1.0_rp))
      scale = max(norm_scalar, norm_optimized, norm_factorized, epsilon(1.0_rp))
      scalar_optimized = matrix_difference(scalar_result%susceptibility(:, :, ifrequency), &
         optimized_result%susceptibility(:, :, ifrequency))/scale
      scalar_factorized = matrix_difference(scalar_result%susceptibility(:, :, ifrequency), &
         factorized_result%susceptibility(:, :, ifrequency))/scale
      optimized_factorized = matrix_difference(optimized_result%susceptibility(:, :, ifrequency), &
         factorized_result%susceptibility(:, :, ifrequency))/scale
      write (*, '(a,1x,a,1x,a,es16.8)') 'FOUR_WAY', 'phase='//trim(phase), 'frequency=', frequency
      write (*, '(a,1x,a,4(1x,es16.8))') 'FOUR_WAY', 'norms_lehmann_scalar_optimized_factorized=', &
         norm_lehmann, norm_scalar, norm_optimized, norm_factorized
      write (*, '(a,1x,a,3(1x,es16.8))') 'FOUR_WAY', 'rF_gf_vs_lehmann_scalar_optimized_factorized=', &
         scalar_lehmann, optimized_lehmann, factorized_lehmann
      write (*, '(a,1x,a,3(1x,es16.8))') 'FOUR_WAY', 'rF_pairwise_scalar_optimized_scalar_factorized_optimized_factorized=', &
         scalar_optimized, scalar_factorized, optimized_factorized
      if (max(scalar_optimized, scalar_factorized, optimized_factorized) >= 1.0e-10_rp) then
         error stop 'GF_CONTRACTION_BACKEND_MISMATCH'
      end if
      call report_matrix_element_diagnostics('scalar', scalar_result%susceptibility(:, :, ifrequency), &
         lehmann_result%susceptibility(:, :, ifrequency))
      call report_matrix_element_diagnostics('optimized', optimized_result%susceptibility(:, :, ifrequency), &
         lehmann_result%susceptibility(:, :, ifrequency))
      call report_matrix_element_diagnostics('factorized', factorized_result%susceptibility(:, :, ifrequency), &
         lehmann_result%susceptibility(:, :, ifrequency))
   end subroutine report_four_way_contraction

   subroutine report_matrix_element_diagnostics(backend, actual, reference)
      character(len=*), intent(in) :: backend
      complex(rp), intent(in) :: actual(:, :), reference(:, :)
      real(rp) :: floor, absolute_error, relative_error, maximum_absolute, maximum_relative
      integer :: i, j, worst_i, worst_j, relative_i, relative_j

      floor = max(1.0e-12_rp, 1.0e-8_rp*maxval(abs(reference)))
      maximum_absolute = 0.0_rp
      maximum_relative = 0.0_rp
      worst_i = 1
      worst_j = 1
      relative_i = 0
      relative_j = 0
      do j = 1, size(actual, 2)
         do i = 1, size(actual, 1)
            absolute_error = abs(actual(i, j) - reference(i, j))
            if (absolute_error > maximum_absolute) then
               maximum_absolute = absolute_error
               worst_i = i
               worst_j = j
            end if
            if (abs(reference(i, j)) >= floor) then
               relative_error = absolute_error/abs(reference(i, j))
               if (relative_error > maximum_relative) then
                  maximum_relative = relative_error
                  relative_i = i
                  relative_j = j
               end if
            end if
         end do
      end do
      write (*, '(a,1x,a,a,1x,a,es16.8,1x,a,es16.8,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,es16.8)') &
         'MATRIX_DIAGNOSTICS', 'backend=', trim(backend), 'max_abs=', maximum_absolute, &
         'max_rel=', maximum_relative, 'worst_i=', worst_i, 'worst_j=', worst_j, &
         'relative_i=', relative_i, 'relative_j=', relative_j, 'amplitude_floor=', floor
   end subroutine report_matrix_element_diagnostics

   real(rp) function matrix_difference(lhs, rhs) result(difference)
      complex(rp), intent(in) :: lhs(:, :), rhs(:, :)

      difference = sqrt(sum(abs(lhs - rhs)**2))
   end function matrix_difference

   subroutine run_rotation_oracle()
      integer, parameter :: ne = 801
      real(rp), parameter :: oracle_integration_eta = 0.010_rp, oracle_margin = 1.0_rp
      real(rp), parameter :: oracle_frequencies(2) = [0.0_rp, 0.17_rp]
      complex(rp), allocatable :: vertices(:, :, :, :), rotated_vertices(:, :, :, :)
      complex(rp), allocatable :: rotated_vectors(:, :, :), gf_base(:, :, :), gf_rotated(:, :, :)
      complex(rp), allocatable :: lehmann_base(:, :, :), lehmann_rotated(:, :, :)
      complex(rp) :: rotation(nbasis, nbasis), rotated_hamiltonian(nbasis, nbasis, nk)
      real(rp) :: rotation_eigenvalues(nbasis)
      type(lr_electronic_state), target :: rotated_state, rotated_endpoint
      type(lr_product_gf_susceptibility_request) :: gf_request
      type(lr_product_gf_susceptibility_result) :: gf_result
      type(lr_product_ks_susceptibility_request) :: lehmann_request
      type(lr_product_ks_susceptibility_result) :: lehmann_result
      real(rp) :: h_residual, d_gf, r_gf, d_lehmann, r_lehmann, d_gf_lehmann
      real(rp) :: d_production_gf, d_production_lehmann, norm_reference
      integer :: ik, ib, ifrequency

      h_residual = 0.0_rp
      call require_resolved_mesh(ne, oracle_integration_eta, oracle_margin, state)
      call product_plus%component_vertex_tensor(vertices)
      call build_rotation(rotation, rotation_eigenvalues)
      allocate(rotated_vectors(nbasis, nbands, nk), rotated_vertices(nbasis, nbasis, lmto_product_nbranch, &
         product_plus%product_dimension))

      do ik = 1, nk
         rotated_hamiltonian(:, :, ik) = matmul(conjg(transpose(rotation)), &
            matmul(hamiltonians(:, :, ik), rotation))
         rotated_vectors(:, :, ik) = matmul(conjg(transpose(rotation)), eigenvectors(:, :, ik))
         do ib = 1, nbands
            h_residual = max(h_residual, maxval(abs(matmul(rotated_hamiltonian(:, :, ik), &
               rotated_vectors(:, ib, ik)) - eigenvalues(ib, ik)*rotated_vectors(:, ib, ik))))
         end do
      end do
      do ifrequency = 1, product_plus%product_dimension
         do ib = 1, lmto_product_nbranch
            rotated_vertices(:, :, ib, ifrequency) = matmul(conjg(transpose(rotation)), &
               matmul(vertices(:, :, ib, ifrequency), rotation))
         end do
      end do
      if (h_residual > 2.0e-11_rp) error stop 'test_lr_product_gf_quadrature_audit: rotated H eigenpair residual'

      call rotated_state%initialize(eigenvalues, rotated_vectors, k_points, k_weights, occupations, 0.0_rp, temperature)
      call rotated_endpoint%initialize(eigenvalues, rotated_vectors, k_points, k_weights, occupations, 0.0_rp, temperature)

      gf_request%q = q_gamma
      gf_request%frequencies = oracle_frequencies
      gf_request%eta = eta
      gf_request%channel = lr_channel_plus
      gf_request%integration_points = ne
      gf_request%integration_eta = oracle_integration_eta
      gf_request%energy_margin = oracle_margin
      gf_request%product_basis => product_plus
      gf_request%electronic_state => state
      gf_request%q_endpoint_state => endpoint
      call evaluate_lr_product_gf_susceptibility(gf_request, gf_result)

      lehmann_request%q = q_gamma
      lehmann_request%frequencies = oracle_frequencies
      lehmann_request%eta = eta
      lehmann_request%channel = lr_channel_plus
      lehmann_request%product_basis => product_plus
      lehmann_request%electronic_state => state
      lehmann_request%q_endpoint_state => endpoint
      call evaluate_lr_product_ks_susceptibility(lehmann_request, lehmann_result)

      call local_gf_response(vertices, state, endpoint, oracle_frequencies, eta, oracle_integration_eta, &
         oracle_margin, ne, gf_base)
      call local_gf_response(rotated_vertices, rotated_state, rotated_endpoint, oracle_frequencies, eta, &
         oracle_integration_eta, oracle_margin, ne, gf_rotated)
      call local_lehmann_response(vertices, state, endpoint, oracle_frequencies, eta, lehmann_base)
      call local_lehmann_response(rotated_vertices, rotated_state, rotated_endpoint, oracle_frequencies, eta, &
         lehmann_rotated)

      d_production_gf = maxval(abs(gf_base - gf_result%susceptibility))
      d_production_lehmann = maxval(abs(lehmann_base - lehmann_result%susceptibility))
      if (d_production_gf > 2.0e-10_rp .or. d_production_lehmann > 2.0e-10_rp) then
         error stop 'test_lr_product_gf_quadrature_audit: local rotation oracle disagrees with production'
      end if
      write (*, '(a,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8)') 'ROTATION_H_EIGENPAIR', 'max=', h_residual, &
         'local_GF_vs_production=', d_production_gf, 'local_Lehmann_vs_production=', d_production_lehmann

      do ifrequency = 1, size(oracle_frequencies)
         d_gf = sqrt(sum(abs(gf_rotated(:, :, ifrequency) - gf_base(:, :, ifrequency))**2))
         norm_reference = sqrt(sum(abs(gf_base(:, :, ifrequency))**2))
         r_gf = d_gf/max(norm_reference, epsilon(1.0_rp))
         d_lehmann = sqrt(sum(abs(lehmann_rotated(:, :, ifrequency) - lehmann_base(:, :, ifrequency))**2))
         norm_reference = sqrt(sum(abs(lehmann_base(:, :, ifrequency))**2))
         r_lehmann = d_lehmann/max(norm_reference, epsilon(1.0_rp))
         d_gf_lehmann = sqrt(sum(abs(gf_base(:, :, ifrequency) - lehmann_base(:, :, ifrequency))**2))
         if (r_gf > 1.0e-10_rp .or. r_lehmann > 1.0e-10_rp) then
            error stop 'test_lr_product_gf_quadrature_audit: orbital rotation response invariant failed'
         end if
         write (*, '(a,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8)') &
            'ROTATION', 'frequency=', oracle_frequencies(ifrequency), 'rGF=', r_gf, 'rLehmann=', r_lehmann, &
            'dGFLehmann=', d_gf_lehmann, 'maxHResidual=', h_residual
      end do
   end subroutine run_rotation_oracle

   subroutine build_rotation(matrix, eigenvalues_out)
      complex(rp), intent(out) :: matrix(:, :)
      real(rp), intent(out) :: eigenvalues_out(:)
      complex(rp) :: seed(nbasis, nbasis)
      integer :: i, j

      seed = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, nbasis
         seed(i, i) = cmplx(0.11_rp*real(i, rp) + 0.003_rp*real(i*i, rp), 0.0_rp, rp)
         do j = i + 1, nbasis
            seed(i, j) = cmplx(0.013_rp*real(i + 2*j, rp), 0.009_rp*real(j - i, rp) + &
               0.0003_rp*real(i + j, rp), rp)
            seed(j, i) = conjg(seed(i, j))
         end do
      end do
      call hermitian_eig(seed, nbasis, eigenvalues_out, matrix)
   end subroutine build_rotation

   subroutine local_lehmann_response(vertices, left_state, right_state, frequencies, response_eta, response)
      complex(rp), intent(in) :: vertices(:, :, :, :)
      type(lr_electronic_state), intent(in) :: left_state, right_state
      real(rp), intent(in) :: frequencies(:), response_eta
      complex(rp), allocatable, intent(out) :: response(:, :, :)
      complex(rp), allocatable :: transition(:)
      complex(rp) :: denominator, pair_factor
      real(rp) :: occupation_difference, weight_sum
      integer :: pd, ik, ib, jb, ifrequency, i, j

      pd = size(vertices, 4)
      allocate(response(pd, pd, size(frequencies)), transition(pd))
      response = cmplx(0.0_rp, 0.0_rp, rp)
      weight_sum = sum(left_state%k_weights)
      do ik = 1, left_state%nk
         do ib = 1, left_state%nbands
            do jb = 1, right_state%nbands
               occupation_difference = left_state%occupations(ib, ik) - right_state%occupations(jb, ik)
               if (occupation_difference == 0.0_rp) cycle
               call local_transition_coordinates(vertices, left_state%eigenvectors(:, ib, ik), &
                  right_state%eigenvectors(:, jb, ik), left_state%eigenvalues(ib, ik), &
                  right_state%eigenvalues(jb, ik), transition)
               do ifrequency = 1, size(frequencies)
                  denominator = cmplx(frequencies(ifrequency) + left_state%eigenvalues(ib, ik) - &
                     right_state%eigenvalues(jb, ik), response_eta, rp)
                  pair_factor = cmplx(2.0_rp*left_state%k_weights(ik)/weight_sum, 0.0_rp, rp)* &
                     occupation_difference/denominator
                  do j = 1, pd
                     do i = 1, pd
                        response(i, j, ifrequency) = response(i, j, ifrequency) + &
                           pair_factor*transition(i)*conjg(transition(j))
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine local_lehmann_response

   subroutine local_transition_coordinates(vertices, left_vector, right_vector, left_energy, right_energy, coordinates)
      complex(rp), intent(in) :: vertices(:, :, :, :), left_vector(:), right_vector(:)
      real(rp), intent(in) :: left_energy, right_energy
      complex(rp), intent(out) :: coordinates(:)
      integer :: component, product_mode, left_power, right_power

      coordinates = cmplx(0.0_rp, 0.0_rp, rp)
      do component = 1, lmto_product_nbranch
         call lmto_product_branch_powers(component, left_power, right_power)
         do product_mode = 1, size(vertices, 4)
            coordinates(product_mode) = coordinates(product_mode) + left_energy**left_power* &
               right_energy**right_power*sum(conjg(left_vector)* &
               matmul(vertices(:, :, component, product_mode), right_vector))
         end do
      end do
   end subroutine local_transition_coordinates

   subroutine local_gf_response(vertices, left_state, right_state, frequencies, response_eta, integration_eta, &
                                 margin, ne, response)
      complex(rp), intent(in) :: vertices(:, :, :, :)
      type(lr_electronic_state), intent(in) :: left_state, right_state
      real(rp), intent(in) :: frequencies(:), response_eta, integration_eta, margin
      integer, intent(in) :: ne
      complex(rp), allocatable, intent(out) :: response(:, :, :)
      complex(rp), allocatable :: left_gr(:, :, :), left_ga(:, :, :), left_a(:, :, :)
      complex(rp), allocatable :: right_gr(:, :, :), right_ga(:, :, :), right_a(:, :, :)
      real(rp) :: energy_min, energy_max, step, energy, fermi_weight, quadrature_weight, weight_sum
      integer :: nb, pd, ie, ik, ifrequency

      nb = left_state%nbasis
      pd = size(vertices, 4)
      energy_min = min(minval(left_state%eigenvalues), minval(right_state%eigenvalues)) - margin
      energy_max = max(maxval(left_state%eigenvalues), maxval(right_state%eigenvalues)) + margin
      step = (energy_max - energy_min)/real(ne - 1, rp)
      weight_sum = sum(left_state%k_weights)
      allocate(response(pd, pd, size(frequencies)), left_gr(nb, nb, lmto_product_max_gf_moment + 1), &
         left_ga(nb, nb, lmto_product_max_gf_moment + 1), left_a(nb, nb, lmto_product_max_gf_moment + 1), &
         right_gr(nb, nb, lmto_product_max_gf_moment + 1), right_ga(nb, nb, lmto_product_max_gf_moment + 1), &
         right_a(nb, nb, lmto_product_max_gf_moment + 1))
      response = cmplx(0.0_rp, 0.0_rp, rp)

      do ie = 1, ne
         energy = energy_min + real(ie - 1, rp)*step
         if (ie == 1 .or. ie == ne) then
            quadrature_weight = 1.0_rp
         else if (mod(ie, 2) == 0) then
            quadrature_weight = 4.0_rp
         else
            quadrature_weight = 2.0_rp
         end if
         quadrature_weight = quadrature_weight*step/3.0_rp
         fermi_weight = lr_fermi_dirac_occupation(energy, left_state%fermi_level, left_state%temperature)
         do ik = 1, left_state%nk
            call build_weighted_resolvent(left_state, ik, cmplx(energy, integration_eta, rp), left_gr)
            call build_weighted_resolvent(left_state, ik, cmplx(energy, -integration_eta, rp), left_ga)
            call build_weighted_resolvent(right_state, ik, cmplx(energy, integration_eta, rp), right_gr)
            call build_weighted_resolvent(right_state, ik, cmplx(energy, -integration_eta, rp), right_ga)
            left_a = cmplx(0.0_rp, 1.0_rp/(2.0_rp*response_angular_pi), rp)*(left_gr - left_ga)
            right_a = cmplx(0.0_rp, 1.0_rp/(2.0_rp*response_angular_pi), rp)*(right_gr - right_ga)
            do ifrequency = 1, size(frequencies)
               call build_weighted_resolvent(right_state, ik, cmplx(energy + frequencies(ifrequency), &
                  response_eta, rp), right_gr)
               call build_weighted_resolvent(left_state, ik, cmplx(energy - frequencies(ifrequency), &
                  -response_eta, rp), left_ga)
               call local_accumulate_gf_bubble(vertices, left_a, right_a, right_gr, left_ga, &
                  fermi_weight*quadrature_weight*2.0_rp*left_state%k_weights(ik)/weight_sum, &
                  response(:, :, ifrequency))
            end do
         end do
      end do
   end subroutine local_gf_response

   subroutine local_accumulate_gf_bubble(vertices, left_a, right_a, right_gr, left_ga, scale, susceptibility)
      complex(rp), intent(in) :: vertices(:, :, :, :), left_a(:, :, :), right_a(:, :, :), &
         right_gr(:, :, :), left_ga(:, :, :)
      real(rp), intent(in) :: scale
      complex(rp), intent(inout) :: susceptibility(:, :)
      complex(rp) :: temporary(size(left_a, 1), size(left_a, 2)), vertex_adjoint(size(left_a, 1), size(left_a, 2))
      integer :: component_i, component_j, left_power_i, right_power_i, left_power_j, right_power_j
      integer :: combined_left, combined_right, i, j, product_dimension

      product_dimension = size(susceptibility, 1)
      do component_i = 1, lmto_product_nbranch
         call lmto_product_branch_powers(component_i, left_power_i, right_power_i)
         do component_j = 1, lmto_product_nbranch
            call lmto_product_branch_powers(component_j, left_power_j, right_power_j)
            combined_left = left_power_i + left_power_j + 1
            combined_right = right_power_i + right_power_j + 1
            do i = 1, product_dimension
               temporary = matmul(left_a(:, :, combined_left), vertices(:, :, component_i, i))
               temporary = matmul(temporary, right_gr(:, :, combined_right))
               do j = 1, product_dimension
                  susceptibility(i, j) = susceptibility(i, j) + scale*sum(temporary* &
                     conjg(vertices(:, :, component_j, j)))
               end do
            end do
            do j = 1, product_dimension
               vertex_adjoint = conjg(transpose(vertices(:, :, component_j, j)))
               temporary = matmul(right_a(:, :, combined_right), vertex_adjoint)
               temporary = matmul(temporary, left_ga(:, :, combined_left))
               do i = 1, product_dimension
                  susceptibility(i, j) = susceptibility(i, j) + scale*sum(temporary* &
                     transpose(vertices(:, :, component_i, i)))
               end do
            end do
         end do
      end do
   end subroutine local_accumulate_gf_bubble

   subroutine hermitian_eig(a_in, n, eigenvalues_out, eigenvectors_out)
      integer, intent(in) :: n
      complex(rp), intent(in) :: a_in(n, n)
      real(rp), intent(out) :: eigenvalues_out(n)
      complex(rp), intent(out) :: eigenvectors_out(n, n)
      complex(rp) :: work_query(1)
      complex(rp), allocatable :: work(:)
      real(rp) :: rwork(max(1, 3*n - 2))
      integer :: info, lwork
      logical :: halt_divide_by_zero
      external :: zheev

      eigenvectors_out = a_in
      call zheev('V', 'U', n, eigenvectors_out, n, eigenvalues_out, work_query, -1, rwork, info)
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      call ieee_get_halting_mode(ieee_divide_by_zero, halt_divide_by_zero)
      call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
      call zheev('V', 'U', n, eigenvectors_out, n, eigenvalues_out, work, lwork, rwork, info)
      call ieee_set_flag(ieee_divide_by_zero, .false.)
      call ieee_set_halting_mode(ieee_divide_by_zero, halt_divide_by_zero)
      deallocate(work)
      if (info /= 0) error stop 'test_lr_product_gf_quadrature_audit: zheev failed'
   end subroutine hermitian_eig

   real(rp) function matrix_residual(actual, reference) result(residual)
      complex(rp), intent(in) :: actual(:, :), reference(:, :)
      residual = sqrt(sum(abs(actual - reference)**2))/max(sqrt(sum(abs(reference)**2)), epsilon(1.0_rp))
   end function matrix_residual

end program test_lr_product_gf_quadrature_audit
