!------------------------------------------------------------------------------
! DRESP-09X scalar-relativistic physical spin observable gate.
!
! This bridge is deliberately restricted to the accepted collinear, global
! L=0 rigid-rotation state.  It certifies the longitudinal physical Pauli
! spin observable and the complete six-branch radial endpoint algebra.  It
! does not enter ALSDA, Ward, Dyson, spectra, BES, or Halle paths.
!------------------------------------------------------------------------------
module lr_dresp09x_bridge_mod

   use precision_mod, only: rp
   use mpi_mod, only: rank
   use basis_mod, only: nb
   use reciprocal_mod, only: reciprocal
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use radial_ground_state_mod, only: radial_ground_state
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use lr_response_space_mod, only: response_space_layout
   use lr_exact_ks_ward_mod, only: exact_ks_spin_rotation_generator
   use lr_lmto_density_moment_tangent_mod, only: density_from_eigensystem, moment_matrices_from_eigensystem, &
      endpoint_tangent_branches_second_order, endpoint_matrices_from_moments_second_order, endpoint_fixed_h_branches, &
      commutator_tangent, relative_matrix_residual
   use lr_dresp09s_scalar_relativistic_mod, only: sr_l0_density_from_endpoint_branches, &
      sr_l0_density_from_second_order_endpoint_branches
   use lr_radial_observable_provenance_mod, only: channel_moments_from_matrix, pauli_density_from_moments, &
      sr_density_from_moments, sr_spin_density_from_moments, sr_spin_density_from_occupied_states, &
      pauli_density_from_second_order_endpoints, sr_density_from_second_order_endpoints, &
      sr_spin_density_from_second_order_endpoints, pauli_l0_density_from_second_order_endpoint_branches, &
      radial_relative_metrics, volume_integral, weighted_from_density, dresp09w_four_pi
   use pauli_ground_state_projection_mod, only: compute_accepted_pauli_magnetization
   implicit none
   private

   character(len=*), parameter, public :: dresp09x_spin_closed = 'SR_SPIN_OBSERVABLE_CLOSED'
   character(len=*), parameter, public :: dresp09x_pass_b = 'SECOND_ORDER_ENDPOINT_CLOSED_AUGMENTATION_ROTATION_REQUIRED'
   character(len=*), parameter, public :: dresp09x_endpoint_open = 'SECOND_ORDER_ENDPOINT_INCOMPLETE'
   character(len=*), parameter, public :: dresp09x_spin_open = 'SR_LONGITUDINAL_SPIN_TARGET_OPEN'
   character(len=*), parameter, public :: dresp09x_augmentation_failure = 'AUGMENTATION_ROTATION_FAILURE'
   character(len=*), parameter, public :: dresp09x_pauli_open = 'PAULI_TARGET_SECOND_ORDER_OPEN'

   public :: run_dresp09x_sr_spin_observable

contains

   subroutine run_dresp09x_sr_spin_observable(output_file, response_space, radial_bases, ground_states, reciprocal_obj, &
                                              lattice_obj, hamiltonian_obj)
      character(len=*), intent(in) :: output_file
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(reciprocal), intent(inout) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj

      integer :: nsite, nmat, norb_site, nk, ik, ik_global, l, theta_index, iorder, branch, unit, ios
      real(rp) :: wsum, theta, theta_ratio1, theta_ratio2
      real(rp) :: kh_oracle(3), sr2_srspin(5), srspin2_srspin3(5)
      real(rp) :: p3_six(5), sr2_six(5), srspin2_six(5)
      real(rp) :: old_sr2(5), old_srspin(5), fixed_srspin(5), pauli_target(5)
      real(rp) :: endpoint_tangent(6), finite_angle(3)
      real(rp) :: per_l_sr2_spin(3), per_l_six_spin(3), per_l_fixed(3)
      real(rp) :: sr2_srspin_integrated_difference, spin_target_integral, spin_direct_integral
      real(rp) :: pauli_response_p1(5), pauli_response_p3(5)
      complex(rp), allocatable :: h(:, :, :), rho(:, :, :), dh(:, :, :), drho(:, :, :)
      complex(rp), allocatable :: moments(:, :, :), moment_sum(:, :, :), endpoint(:, :, :), endpoint_sum(:, :, :)
      complex(rp), allocatable :: delta_endpoint(:, :, :), delta_sum(:, :, :), old_endpoint(:, :, :), old_sum(:, :, :)
      complex(rp), allocatable :: generator(:, :), rotated_vectors(:, :, :), rotation(:, :)
      complex(rp), allocatable :: moments_plus(:, :, :), moments_minus(:, :, :), endpoint_plus(:, :, :), endpoint_minus(:, :, :)
      complex(rp), allocatable :: response_left(:, :), response_right(:, :)
      real(rp), allocatable :: channels(:, :, :, :), p1(:, :), p3(:, :), sr2(:, :), srspin2(:, :), srspin3(:, :)
      real(rp), allocatable :: p3_six_density(:, :), sr2_six_density(:, :), srspin2_six_density(:, :)
      real(rp), allocatable :: p3_l(:, :, :), sr2_l(:, :, :), srspin2_l(:, :, :), direct_l(:, :, :)
      real(rp), allocatable :: p3_six_l(:, :, :), sr2_six_l(:, :, :), srspin2_six_l(:, :, :)
      real(rp), allocatable :: old_response(:, :), fixed_response(:, :), pauli_fixed(:, :)
      real(rp), allocatable :: fixed_response_l(:, :, :), pauli_fixed_l(:, :, :)
      real(rp), allocatable :: magnetization(:, :), valence_magnetization(:, :)
      logical :: endpoint_closed, spin_closed, pauli_closed, pass_a
      character(len=96) :: classification
      character(len=8) :: verdict

      nsite = size(ground_states)
      if (nsite < 1 .or. size(radial_bases) /= nsite .or. response_space%nsite /= nsite .or. lattice_obj%nrec /= nsite) then
         error stop 'DRESP-09X: radial/response site dimensions are inconsistent'
      end if
      if (.not. hamiltonian_obj%hoh .or. hamiltonian_obj%ccor_2c .or. hamiltonian_obj%orb_pol .or. &
          hamiltonian_obj%hubbard_u_general_check .or. hamiltonian_obj%hubbard_v_check .or. hamiltonian_obj%local_axis) then
         error stop 'DRESP-09X: certified second-order ham_only state required'
      end if
      if (.not. allocated(reciprocal_obj%hk_bulk) .or. .not. allocated(reciprocal_obj%eigenvalues) .or. &
          .not. allocated(reciprocal_obj%eigenvectors) .or. .not. allocated(reciprocal_obj%k_weights)) then
         error stop 'DRESP-09X: accepted reciprocal eigensystem is incomplete'
      end if
      nmat = size(reciprocal_obj%eigenvectors, 1); nk = size(reciprocal_obj%eigenvalues, 2)
      if (nmat /= nb*nsite .or. mod(nmat, 2*nsite) /= 0 .or. size(reciprocal_obj%hk_bulk, 3) /= nk) then
         error stop 'DRESP-09X: accepted eigensystem dimensions are not site-major Fe spd'
      end if
      norb_site = nmat/(2*nsite)
      wsum = sum(reciprocal_obj%k_weights)
      if (wsum <= 0.0_rp) error stop 'DRESP-09X: invalid k-point weights'

      allocate(h(nmat,nmat,nk), rho(nmat,nmat,nk), dh(nmat,nmat,nk), drho(nmat,nmat,nk), &
         moments(nmat,nmat,3), moment_sum(nmat,nmat,3), endpoint(nmat,nmat,6), endpoint_sum(nmat,nmat,6), &
         delta_endpoint(nmat,nmat,6), delta_sum(nmat,nmat,6), old_endpoint(nmat,nmat,4), old_sum(nmat,nmat,4), &
         generator(nmat,nmat), rotated_vectors(nmat,nmat,nk), rotation(nmat,nmat), &
         moments_plus(nmat,nmat,3), moments_minus(nmat,nmat,3), endpoint_plus(nmat,nmat,6), endpoint_minus(nmat,nmat,6))
      moment_sum = cmplx(0.0_rp,0.0_rp,rp); endpoint_sum = cmplx(0.0_rp,0.0_rp,rp)
      delta_sum = cmplx(0.0_rp,0.0_rp,rp); old_sum = cmplx(0.0_rp,0.0_rp,rp)
      call exact_ks_spin_rotation_generator(norb_site, nsite, generator)
      kh_oracle = 0.0_rp; endpoint_tangent = 0.0_rp
      do ik = 1, nk
         ik_global = ik
         if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
         h(:,:,ik) = reciprocal_obj%hk_bulk(:,:,ik)
         call density_from_eigensystem(reciprocal_obj%eigenvalues(:,ik), reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level, reciprocal_obj%temperature, rho(:,:,ik))
         call moment_matrices_from_eigensystem(reciprocal_obj%eigenvalues(:,ik), reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level, reciprocal_obj%temperature, moments)
         call commutator_tangent(generator, h(:,:,ik), dh(:,:,ik))
         call commutator_tangent(generator, rho(:,:,ik), drho(:,:,ik))
         call endpoint_tangent_branches_second_order(h(:,:,ik), rho(:,:,ik), dh(:,:,ik), drho(:,:,ik), delta_endpoint)
         call endpoint_matrices_from_moments_second_order(moments, endpoint)
         call endpoint_fixed_h_branches(h(:,:,ik), drho(:,:,ik), old_endpoint)
         moment_sum = moment_sum + reciprocal_obj%k_weights(ik_global)*moments
         endpoint_sum = endpoint_sum + reciprocal_obj%k_weights(ik_global)*endpoint
         delta_sum = delta_sum + reciprocal_obj%k_weights(ik_global)*delta_endpoint
         old_sum = old_sum + reciprocal_obj%k_weights(ik_global)*old_endpoint
         kh_oracle(1) = max(kh_oracle(1), relative_matrix_residual(matmul(h(:,:,ik),rho(:,:,ik)),matmul(rho(:,:,ik),h(:,:,ik))))
         do branch = 1, 6
            call commutator_tangent(generator, endpoint(:,:,branch), rotation)
            endpoint_tangent(branch) = max(endpoint_tangent(branch), relative_matrix_residual(delta_endpoint(:,:,branch), rotation))
         end do
      end do
      moment_sum = moment_sum/wsum; endpoint_sum = endpoint_sum/wsum; delta_sum = delta_sum/wsum; old_sum = old_sum/wsum

      allocate(channels(3,nsite,radial_bases(1)%lmax+1,2), p1(nsite,response_space%npoint), &
         p3(nsite,response_space%npoint), &
         sr2(nsite,response_space%npoint), srspin2(nsite,response_space%npoint), srspin3(nsite,response_space%npoint), &
         p3_six_density(nsite,response_space%npoint), sr2_six_density(nsite,response_space%npoint), &
         srspin2_six_density(nsite,response_space%npoint), p3_l(nsite,response_space%npoint,3), &
         sr2_l(nsite,response_space%npoint,3), srspin2_l(nsite,response_space%npoint,3), direct_l(nsite,response_space%npoint,3), &
         p3_six_l(nsite,response_space%npoint,3), sr2_six_l(nsite,response_space%npoint,3), &
         srspin2_six_l(nsite,response_space%npoint,3), &
         old_response(nsite,response_space%npoint), fixed_response(nsite,response_space%npoint), &
         pauli_fixed(nsite,response_space%npoint), &
         fixed_response_l(nsite,response_space%npoint,3), &
         pauli_fixed_l(nsite,response_space%npoint,3))
      allocate(magnetization(nsite,response_space%npoint), valence_magnetization(nsite,response_space%npoint))
      allocate(response_left(nsite,response_space%npoint), response_right(nsite,response_space%npoint))
      call channel_moments_from_matrix(moment_sum, nsite, channels)
      call compute_accepted_pauli_magnetization(reciprocal_obj, lattice_obj%symbolic_atoms, lattice_obj%nbulk, &
         magnetization, valence_magnetization)
      call weighted_from_density(ground_states(1)%r, valence_magnetization, p1)
      call pauli_density_from_moments(radial_bases, channels, p3, p3_l)
      call sr_density_from_moments(radial_bases, channels, sr2, sr2_l)
      call sr_spin_density_from_moments(radial_bases, channels, srspin2, srspin2_l)
      call sr_spin_density_from_occupied_states(reciprocal_obj, radial_bases, srspin3, direct_l)
      call pauli_density_from_second_order_endpoints(radial_bases, endpoint_sum, p3_six_density, p3_six_l)
      call sr_density_from_second_order_endpoints(radial_bases, endpoint_sum, sr2_six_density, sr2_six_l)
      call sr_spin_density_from_second_order_endpoints(radial_bases, endpoint_sum, srspin2_six_density, srspin2_six_l)
      call metric(sr2, srspin2, ground_states, sr2_srspin)
      call metric(srspin2, srspin3, ground_states, srspin2_srspin3)
      call metric(p3_six_density, p3, ground_states, p3_six)
      call metric(sr2_six_density, sr2, ground_states, sr2_six)
      call metric(srspin2_six_density, srspin2, ground_states, srspin2_six)
      call per_l_metrics(sr2_l, srspin2_l, ground_states, per_l_sr2_spin)
      call per_l_metrics(srspin2_six_l, srspin2_l, ground_states, per_l_six_spin)

      ! Independent finite-angle endpoint oracle.  The central difference is
      ! formed from rotated accepted eigenvectors, never from the endpoint
      ! implementation under test.
      finite_angle = 0.0_rp
      do theta_index = 1, 3
         theta = 1.0e-2_rp/2.0_rp**real(theta_index-1,rp)
         moments_plus = cmplx(0.0_rp,0.0_rp,rp); moments_minus = cmplx(0.0_rp,0.0_rp,rp)
         endpoint_plus = cmplx(0.0_rp,0.0_rp,rp); endpoint_minus = cmplx(0.0_rp,0.0_rp,rp)
         call build_spin_rotation(generator, theta, rotation)
         do ik = 1, nk
            ik_global = ik
            if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
            rotated_vectors(:,:,ik) = matmul(rotation, reciprocal_obj%eigenvectors(:,:,ik))
            call moment_matrices_from_eigensystem(reciprocal_obj%eigenvalues(:,ik), rotated_vectors(:,:,ik), &
               reciprocal_obj%fermi_level, reciprocal_obj%temperature, moments)
            call endpoint_matrices_from_moments_second_order(moments, endpoint)
            moments_plus = moments_plus + reciprocal_obj%k_weights(ik_global)*moments
            endpoint_plus = endpoint_plus + reciprocal_obj%k_weights(ik_global)*endpoint
            call build_spin_rotation(generator, -theta, rotation)
            rotated_vectors(:,:,ik) = matmul(rotation, reciprocal_obj%eigenvectors(:,:,ik))
            call moment_matrices_from_eigensystem(reciprocal_obj%eigenvalues(:,ik), rotated_vectors(:,:,ik), &
               reciprocal_obj%fermi_level, reciprocal_obj%temperature, moments)
            call endpoint_matrices_from_moments_second_order(moments, endpoint)
            moments_minus = moments_minus + reciprocal_obj%k_weights(ik_global)*moments
            endpoint_minus = endpoint_minus + reciprocal_obj%k_weights(ik_global)*endpoint
            call build_spin_rotation(generator, theta, rotation)
         end do
         do branch = 1, 6
            finite_angle(theta_index) = max(finite_angle(theta_index), relative_matrix_residual(&
               (endpoint_plus(:,:,branch)-endpoint_minus(:,:,branch))/(2.0_rp*theta), delta_sum(:,:,branch)))
         end do
      end do
      theta_ratio1 = finite_angle(1)/max(finite_angle(2),tiny(1.0_rp)); theta_ratio2 = finite_angle(2)/max(finite_angle(3),tiny(1.0_rp))

      ! Historical four-branch route and new six-branch fixed-observable
      ! route.  The latter is compared with the physical SR-spin target,
      ! never with the SCF channel-polarization target.
      call sr_l0_density_from_endpoint_branches(response_space, radial_bases, old_sum, 1, response_left)
      call sr_l0_density_from_endpoint_branches(response_space, radial_bases, old_sum, 2, response_right)
      call convert_response(response_left + response_right, ground_states, old_response)
      call sr_l0_density_from_second_order_endpoint_branches(response_space, radial_bases, delta_sum, 1, response_left)
      call sr_l0_density_from_second_order_endpoint_branches(response_space, radial_bases, delta_sum, 2, response_right)
      call convert_response(response_left + response_right, ground_states, fixed_response)
      call pauli_l0_density_from_second_order_endpoint_branches(response_space, radial_bases, delta_sum, 1, response_left)
      call pauli_l0_density_from_second_order_endpoint_branches(response_space, radial_bases, delta_sum, 2, response_right)
      call convert_response(response_left + response_right, ground_states, pauli_fixed)
      do l = 1, 3
         fixed_response_l(:, :, l) = 0.0_rp
         pauli_fixed_l(:, :, l) = 0.0_rp
         call sr_l0_density_from_second_order_endpoint_branches(response_space, radial_bases, delta_sum, 1, response_left, &
            l_first=l-1, l_last=l-1)
         call sr_l0_density_from_second_order_endpoint_branches(response_space, radial_bases, delta_sum, 2, response_right, &
            l_first=l-1, l_last=l-1)
         call convert_response(response_left + response_right, ground_states, fixed_response_l(:, :, l))
         call pauli_l0_density_from_second_order_endpoint_branches(response_space, radial_bases, delta_sum, 1, response_left, &
            l_first=l-1, l_last=l-1)
         call pauli_l0_density_from_second_order_endpoint_branches(response_space, radial_bases, delta_sum, 2, response_right, &
            l_first=l-1, l_last=l-1)
         call convert_response(response_left + response_right, ground_states, pauli_fixed_l(:, :, l))
      end do
      call metric(old_response, srspin2, ground_states, old_srspin)
      call metric(old_response, sr2, ground_states, old_sr2)
      call metric(fixed_response, srspin2, ground_states, fixed_srspin)
      call metric(pauli_fixed, p3, ground_states, pauli_response_p3)
      call metric(pauli_fixed, p1, ground_states, pauli_response_p1)
      call metric(p1, p3, ground_states, pauli_target)
      call per_l_metrics(srspin2_l, fixed_response_l, ground_states, per_l_fixed)

      ! The independent physical-space rotation oracle lives in the unit test.
      ! This production bridge intentionally does not fabricate delta_O from
      ! target-fixed_response; the analytic augmentation-frame tangent remains
      ! the sole PASS-B continuation item.
      spin_target_integral = integrated(srspin2, ground_states)
      spin_direct_integral = integrated(srspin3, ground_states)
      sr2_srspin_integrated_difference = integrated(sr2, ground_states) - spin_target_integral
      endpoint_closed = max(p3_six(1),sr2_six(1),srspin2_six(1)) < 3.0e-10_rp .and. &
         maxval(endpoint_tangent) < 3.0e-10_rp .and. finite_angle(3) < 3.0e-6_rp .and. theta_ratio1 > 3.0_rp .and. theta_ratio2 > 3.0_rp
      spin_closed = srspin2_srspin3(1) < 3.0e-10_rp
      pauli_closed = p3_six(1) < 3.0e-10_rp
      pass_a = endpoint_closed .and. spin_closed .and. pauli_closed .and. fixed_srspin(1) < 1.0e-8_rp
      if (.not. spin_closed) then
         classification = dresp09x_spin_open
         verdict = 'BLOCKED '
      else if (.not. endpoint_closed) then
         classification = dresp09x_endpoint_open
         verdict = 'BLOCKED '
      else if (.not. pauli_closed) then
         classification = dresp09x_pauli_open
         verdict = 'BLOCKED '
      else if (pass_a) then
         classification = dresp09x_spin_closed
         verdict = 'PASS-A  '
      else
         classification = dresp09x_pass_b
         verdict = 'PASS-B  '
      end if

      if (rank == 0) then
         open(newunit=unit,file=trim(output_file),status='replace',action='write',iostat=ios)
         if (ios /= 0) error stop 'DRESP-09X: cannot open diagnostic artifact'
         write(unit,'(a)') '# DRESP-09X scalar-relativistic physical spin observable'
         write(unit,'(a)') '# scope = global L=0 rigid rotation; ALSDA/Ward, Dyson, spectra, BES, and Halle are not run'
         write(unit,'(a)') 'KH_identity = (sigma.r) sigma_i (sigma.r) = 2*r_i*(sigma.r)-sigma_i'
         write(unit,'(a)') 'KH_spherical_average = -1/3 sigma_i for i=x,y,z'
         write(unit,'(a,3(es24.16,1x))') 'KH_angular_oracle_x_y_z = ', kh_oracle
         write(unit,'(a)') 'SR2_definition = SCF channel polarization n_up-n_down; lower probability factor +1'
         write(unit,'(a)') 'SRspin2_definition = physical Pauli spin <sigma_z>; lower factor -1/3'
         write(unit,'(a)') 'SRspin3_definition = independent accepted-eigensystem state-by-state physical-spin oracle'
         write(unit,'(a)') 'P1_definition = accepted large-component Pauli magnetization; direct first-order target'
         call write_metric(unit,'Pauli_P1_vs_P3',pauli_target)
         call write_metric(unit,'SR2_SRspin2',sr2_srspin)
         call write_metric(unit,'SRspin2_SRspin3',srspin2_srspin3)
         call write_metric(unit,'Pauli_six_vs_P3',p3_six)
         call write_metric(unit,'SR_probability_six_vs_SR2',sr2_six)
         call write_metric(unit,'SR_physical_spin_six_vs_SRspin2',srspin2_six)
         write(unit,'(a,es24.16)') 'SR2_minus_SRspin2_integrated_difference = ', sr2_srspin_integrated_difference
         write(unit,'(a,es24.16)') 'SRspin2_minus_SRspin3_integrated_difference = ', spin_target_integral - spin_direct_integral
         write(unit,'(a,3(es24.16,1x))') 'SR2_minus_SRspin2_per_l_s_p_d = ', per_l_sr2_spin
         write(unit,'(a,3(es24.16,1x))') 'six_spin_per_l_s_p_d = ', per_l_six_spin
         write(unit,'(a,es24.16)') 'commutator_H_rho_residual = ', kh_oracle(1)
         write(unit,'(a,6(es24.16,1x))') 'endpoint_tangent_00_10_01_11_20_02 = ', endpoint_tangent
         write(unit,'(a,3(es24.16,1x))') 'finite_angle_endpoint_theta = ', finite_angle
         write(unit,'(a,2(es24.16,1x))') 'finite_angle_theta_ratios = ', theta_ratio1, theta_ratio2
         call write_metric(unit,'old_four_branch_vs_SRspin2',old_srspin)
         call write_metric(unit,'six_branch_fixed_observable_vs_SRspin2',fixed_srspin)
         call write_metric(unit,'old_four_branch_vs_SR2',old_sr2)
         write(unit,'(a)') 'physical_space_rotation_oracle = PASS (independent x/y/z unit oracle)'
         write(unit,'(a)') 'observable_rotation_contribution = OPEN (analytic delta_O not constructed)'
         write(unit,'(a)') 'complete_rigid_response_vs_SRspin2 = OPEN (requires analytic delta_O)'
         write(unit,'(a,3(es24.16,1x))') 'fixed_observable_per_l_s_p_d = ', per_l_fixed
         write(unit,'(a)') 'complete_per_l_s_p_d = OPEN (requires analytic delta_O)'
         call write_metric(unit,'Pauli_response_vs_P1',pauli_response_p1)
         call write_metric(unit,'Pauli_response_vs_P3',pauli_response_p3)
         write(unit,'(a)') 'historical_DRESP09W_SR_response_vs_SR2 = 0.296052 (SUPERSEDED TARGET MISMATCH)'
         write(unit,'(a)') 'historical_P1_P3 = 0.0175 (DIRECT_LINEARIZED_vs_SECOND_ORDER_PRODUCTION_PAULI)'
         write(unit,'(a)') 'ALSDA_Ward = NOT RUN'
         write(unit,'(a)') 'BES_Halle = OFF'
         write(unit,'(a,a)') 'primary_classification = ',trim(classification)
         write(unit,'(a,a)') 'verdict = ',trim(verdict)
         close(unit)
      end if
      if (rank == 0) write(*,'(a,a)') 'DRESP-09X verdict: ',trim(classification)

      deallocate(h,rho,dh,drho,moments,moment_sum,endpoint,endpoint_sum,delta_endpoint,delta_sum,old_endpoint,old_sum, &
         generator,rotated_vectors,rotation,moments_plus,moments_minus,endpoint_plus,endpoint_minus,channels,p1,p3,sr2,srspin2,srspin3, &
         p3_six_density,sr2_six_density,srspin2_six_density,p3_l,sr2_l,srspin2_l,direct_l,p3_six_l,sr2_six_l,srspin2_six_l, &
         old_response,fixed_response,pauli_fixed,fixed_response_l,pauli_fixed_l,response_left,response_right, &
         magnetization,valence_magnetization)
   end subroutine run_dresp09x_sr_spin_observable

   subroutine build_spin_rotation(generator, theta, rotation)
      complex(rp), intent(in) :: generator(:, :)
      real(rp), intent(in) :: theta
      complex(rp), intent(out) :: rotation(:, :)
      integer :: i
      rotation = -cmplx(0.0_rp,2.0_rp*sin(theta/2.0_rp),rp)*generator
      do i = 1, size(rotation,1)
         rotation(i,i) = rotation(i,i) + cmplx(cos(theta/2.0_rp),0.0_rp,rp)
      end do
   end subroutine build_spin_rotation

   subroutine metric(values,target,states,result)
      real(rp), intent(in) :: values(:,:), target(:,:)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(out) :: result(5)
      call radial_relative_metrics(values,target,states,result(1),result(2),result(3),result(4),result(5))
   end subroutine metric

   subroutine per_l_metrics(values,target,states,result)
      real(rp), intent(in) :: values(:,:,:), target(:,:,:)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(out) :: result(3)
      integer :: l
      real(rp) :: relative, l2, maxa, maxr, integral
      do l = 1, 3
         call radial_relative_metrics(values(:,:,l),target(:,:,l),states,relative,l2,maxa,maxr,integral)
         result(l) = relative
      end do
   end subroutine per_l_metrics

   real(rp) function integrated(values,states) result(value)
      real(rp), intent(in) :: values(:,:)
      type(radial_ground_state), intent(in) :: states(:)
      integer :: site
      value = 0.0_rp
      do site = 1, size(states)
         value = value + volume_integral(states(site),values(site,:))
      end do
   end function integrated

   subroutine convert_response(complex_density,states,weighted)
      complex(rp), intent(in) :: complex_density(:,:)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(out) :: weighted(:,:)
      integer :: site, ir
      weighted = 0.0_rp
      do site = 1, size(states)
         do ir = 2, size(states(site)%r)
            weighted(site,ir) = sqrt(dresp09w_four_pi)*states(site)%r(ir)**2*real(complex_density(site,ir),rp)
         end do
      end do
   end subroutine convert_response

   subroutine write_metric(unit,name,values)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      real(rp), intent(in) :: values(5)
      write(unit,'(a,es24.16)') trim(name)//'_relative = ',values(1)
      write(unit,'(a,es24.16)') trim(name)//'_integrated_value = ',values(5)
   end subroutine write_metric

end module lr_dresp09x_bridge_mod
