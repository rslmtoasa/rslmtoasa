!------------------------------------------------------------------------------
! TDVK-03A reciprocal-GF real-axis quadrature closure audit.
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
   use lr_ks_susceptibility_mod, only: lr_electronic_state, &
      lr_product_ks_susceptibility_request, lr_product_ks_susceptibility_result, &
      lr_channel_plus, lr_channel_minus, lr_fermi_dirac_occupation, &
      evaluate_lr_product_ks_susceptibility
   use lr_gf_susceptibility_mod, only: build_weighted_resolvent
   use lr_product_gf_susceptibility_mod, only: lr_product_gf_susceptibility_request, &
      lr_product_gf_susceptibility_result, evaluate_lr_product_gf_susceptibility
   implicit none

   integer, parameter :: nr = 7, orbital_lmax = 1, response_lmax = 2, nsite = 1
   integer, parameter :: norb = (orbital_lmax + 1)**2
   integer, parameter :: nbasis = 2*norb*nsite, nbands = nbasis, nk = 2
   integer, parameter :: nfixed = 7, neta = 4, nmargin = 3
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
   real(rp), parameter :: gamma_frequencies(2) = [0.0_rp, 0.17_rp]
   real(rp), parameter :: static_frequency(1) = [0.0_rp]

   real(rp) :: radius(nr), eigenvalues(nbands, nk), k_points(3, nk)
   real(rp) :: k_weights(nk), occupations(nbands, nk)
   complex(rp) :: eigenvectors(nbasis, nbands, nk)
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
   call build_electronic_fixture(eigenvalues, eigenvectors, k_points, k_weights, occupations)
   call state%initialize(eigenvalues, eigenvectors, k_points, k_weights, occupations, 0.0_rp, temperature)
   call endpoint%initialize(eigenvalues, eigenvectors, k_points, k_weights, occupations, 0.0_rp, temperature)
   call product_plus%initialize(space, radial, lmto_product_channel_plus, .false.)
   call product_minus%initialize(space, radial, lmto_product_channel_minus, .false.)

   write (*, '(a)') 'TDVK-03A reciprocal-GF quadrature audit'
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
   case default
      error stop 'test_lr_product_gf_quadrature_audit: invalid audit mode'
   end select

   write (*, '(a)') 'TDVK-03A audit execution complete'

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
               end do
            end do
         end do
      end do
   end subroutine setup_radial_basis

   subroutine build_electronic_fixture(values, vectors, points, weights, occ)
      real(rp), intent(out) :: values(:, :), points(:, :), weights(:), occ(:, :)
      complex(rp), intent(out) :: vectors(:, :, :)
      integer :: ib, ik

      values = 0.0_rp
      do ik = 1, nk
         do ib = 1, nbands
            values(ib, ik) = -0.72_rp + 0.19_rp*real(ib - 1, rp) + 0.037_rp*real(ik - 1, rp)
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
      allocate(green_r(electronic_state%nbasis, electronic_state%nbasis, 3), &
         green_a(electronic_state%nbasis, electronic_state%nbasis, 3), &
         spectral(electronic_state%nbasis, electronic_state%nbasis, 3), &
         numerical(electronic_state%nbasis, electronic_state%nbasis, 3), &
         numerical_fermi(electronic_state%nbasis, electronic_state%nbasis, 3), &
         exact(electronic_state%nbasis, electronic_state%nbasis, 3), &
         exact_fermi(electronic_state%nbasis, electronic_state%nbasis, 3))

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
         do p = 0, 2
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

   subroutine report_response(phase, product, channel, q, frequencies, ne, integration_eta, margin)
      character(len=*), intent(in) :: phase, channel
      type(lmto_product_response_basis), target, intent(in) :: product
      real(rp), intent(in) :: q(3), frequencies(:), integration_eta, margin
      integer, intent(in) :: ne

      type(lr_product_gf_susceptibility_request) :: gf_request
      type(lr_product_gf_susceptibility_result) :: gf_result
      type(lr_product_ks_susceptibility_request) :: lehmann_request
      type(lr_product_ks_susceptibility_result) :: lehmann_result
      type(lr_electronic_state), target :: local_endpoint
      integer :: ifrequency
      real(rp) :: d_f, r_f, d_inf, h, ratio, norm_reference

      call local_endpoint%initialize(eigenvalues, eigenvectors, spread(q, 2, nk), k_weights, occupations, &
         0.0_rp, temperature)
      gf_request%q = q
      gf_request%frequencies = frequencies
      gf_request%eta = eta
      gf_request%channel = channel
      gf_request%integration_points = ne
      gf_request%integration_eta = integration_eta
      gf_request%energy_margin = margin
      gf_request%product_basis => product
      gf_request%electronic_state => state
      gf_request%q_endpoint_state => local_endpoint
      call evaluate_lr_product_gf_susceptibility(gf_request, gf_result)

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
      end do
      deallocate(local_endpoint%k_points, local_endpoint%k_weights, local_endpoint%eigenvalues, &
         local_endpoint%eigenvectors, local_endpoint%occupations)
   end subroutine report_response

   real(rp) function matrix_residual(actual, reference) result(residual)
      complex(rp), intent(in) :: actual(:, :), reference(:, :)
      residual = sqrt(sum(abs(actual - reference)**2))/max(sqrt(sum(abs(reference)**2)), epsilon(1.0_rp))
   end function matrix_residual

end program test_lr_product_gf_quadrature_audit
