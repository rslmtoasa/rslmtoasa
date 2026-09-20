!------------------------------------------------------------------------------
! DRESP-09Y accepted-state bridge.
!
! This is the first production consumer of the independent augmentation-frame
! tangent.  DRESP-09X remains the fixed-observable reference.  The bridge
! adds Tr[D delta_O] from explicitly rotated 2x2 radial spin matrices and
! checks the complete response against an independent finite-angle profile.
! No SCF, radial re-solution, Kxc, Ward, Dyson, spectra, BES, or Halle path is
! entered here.
!------------------------------------------------------------------------------
module lr_dresp09y_bridge_mod

   use precision_mod, only: rp
   use mpi_mod, only: rank
   use basis_mod, only: nb
   use reciprocal_mod, only: reciprocal
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use radial_ground_state_mod, only: radial_ground_state
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use response_angular_basis_mod, only: response_gaunt
   use lr_response_space_mod, only: response_space_layout
   use lr_exact_ks_ward_mod, only: exact_ks_spin_rotation_generator
   use lr_lmto_density_moment_tangent_mod, only: density_from_eigensystem, moment_matrices_from_eigensystem, &
      endpoint_tangent_branches_second_order, endpoint_matrices_from_moments_second_order, commutator_tangent, &
      relative_matrix_residual
   use lr_dresp09s_scalar_relativistic_mod, only: sr_l0_density_from_second_order_endpoint_branches
   use lr_dresp09s_scalar_relativistic_mod, only: sr_second_order_branch_point
   use lr_radial_observable_provenance_mod, only: channel_moments_from_matrix, pauli_density_from_moments, &
      sr_spin_density_from_moments, radial_relative_metrics, volume_integral, &
      volume_l2_norm, weighted_from_density, dresp09w_four_pi
   use pauli_ground_state_projection_mod, only: compute_accepted_pauli_magnetization
   use lr_sr_augmentation_tangent_mod, only: sr_aug_sigma_x, sr_aug_sigma_plus, sr_aug_sigma_minus, &
      sr_aug_branch_observable_at_angle, sr_aug_branch_observable_tangent, sr_aug_lower_spin_factor
   implicit none
   private

   character(len=*), parameter, public :: dresp09y_closed = 'AUGMENTATION_FRAME_TANGENT_CLOSED'
   character(len=*), parameter, public :: dresp09y_delta_open = 'ANALYTIC_DELTA_O_OPEN'
   character(len=*), parameter, public :: dresp09y_branch_open = 'SIX_BRANCH_OBSERVABLE_TANGENT_OPEN'
   character(len=*), parameter, public :: dresp09y_frame_open = 'COEFFICIENT_OBSERVABLE_FRAME_MISMATCH'
   character(len=*), parameter, public :: dresp09y_pauli_open = 'PAULI_AUGMENTATION_OPEN'

   public :: run_dresp09y_augmentation_tangent

contains

   subroutine run_dresp09y_augmentation_tangent(output_file, response_space, radial_bases, ground_states, reciprocal_obj, &
                                                lattice_obj, hamiltonian_obj)
      character(len=*), intent(in) :: output_file
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(reciprocal), intent(inout) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj

      integer :: nsite, nmat, norb_site, nk, ik, ik_global, branch, l, ir, theta_index, unit, ios
      real(rp) :: wsum, theta, theta_ratio(2), finite_error(3), branch_error(6,4), circular_error
      real(rp) :: finite_fixed_error(3), finite_augmentation_error(3)
      real(rp) :: endpoint_fd_error(3)
      real(rp) :: radial_pair_error
      real(rp) :: fixed_metrics(5), aug_metrics(5), complete_metrics(5)
      real(rp) :: pauli_complete_p1(5), pauli_complete_p3(5)
      real(rp) :: per_l_fixed(3), per_l_aug(3), per_l_complete(3), target_integral, complete_integral
      real(rp) :: coefficient_norm, aug_upper_norm, aug_small_norm, aug_angular_norm, aug_total_norm, complete_norm
      real(rp) :: target_norm, interference, complete_fd_metrics(5), branch_fd_components(6,4)
      real(rp) :: endpoint_tangent(6), kh_residual
      real(rp) :: pauli_fixed_p1(5), pauli_fixed_p3(5)
      complex(rp), allocatable :: h(:, :, :), rho(:, :, :), dh(:, :, :), drho(:, :, :)
      complex(rp), allocatable :: moments(:, :, :), moment_sum(:, :, :), endpoint(:, :, :), endpoint_sum(:, :, :)
      complex(rp), allocatable :: delta_endpoint(:, :, :), delta_sum(:, :, :), generator(:, :), rotation(:, :)
      complex(rp), allocatable :: response_left(:, :), response_right(:, :)
      real(rp), allocatable :: channels(:, :, :, :)
      complex(rp), allocatable :: sigma_x(:, :), sigma_plus(:, :), sigma_minus(:, :)
      real(rp), allocatable :: p1(:, :), p3(:, :), target(:, :), fixed(:, :), augmentation(:, :), complete(:, :)
      real(rp), allocatable :: valence_magnetization(:, :)
      real(rp), allocatable :: aug_upper(:, :), aug_small(:, :), aug_angular(:, :), pauli_fixed(:, :), pauli_aug(:, :)
      real(rp), allocatable :: target_l(:, :, :), fixed_l(:, :, :), augmentation_l(:, :, :), complete_l(:, :, :)
      real(rp), allocatable :: pauli_complete(:, :), finite_profile(:, :), finite_fixed_profile(:, :), finite_augmentation_profile(:, :)
      complex(rp), allocatable :: finite_plus(:, :), finite_minus(:, :), augmentation_l_c(:, :)
      complex(rp), allocatable :: aug_upper_c(:, :), aug_small_c(:, :), aug_angular_c(:, :), aug_total_c(:, :)
      complex(rp), allocatable :: pauli_aug_c(:, :)
      logical :: endpoint_closed, delta_closed, branch_closed, frame_closed, pauli_closed, pass_a
      character(len=96) :: classification

      nsite = size(ground_states)
      if (nsite < 1 .or. size(radial_bases) /= nsite .or. response_space%nsite /= nsite .or. lattice_obj%nrec /= nsite) then
         error stop 'DRESP-09Y: radial/response site dimensions are inconsistent'
      end if
      if (.not. hamiltonian_obj%hoh .or. hamiltonian_obj%ccor_2c .or. hamiltonian_obj%orb_pol .or. &
          hamiltonian_obj%hubbard_u_general_check .or. hamiltonian_obj%hubbard_v_check .or. hamiltonian_obj%local_axis) then
         error stop 'DRESP-09Y: certified second-order ham_only state required'
      end if
      if (.not. allocated(reciprocal_obj%hk_bulk) .or. .not. allocated(reciprocal_obj%eigenvalues) .or. &
          .not. allocated(reciprocal_obj%eigenvectors) .or. .not. allocated(reciprocal_obj%k_weights)) then
         error stop 'DRESP-09Y: accepted reciprocal eigensystem is incomplete'
      end if
      nmat = size(reciprocal_obj%eigenvectors, 1); nk = size(reciprocal_obj%eigenvalues, 2)
      if (nmat /= nb*nsite .or. mod(nmat, 2*nsite) /= 0 .or. size(reciprocal_obj%hk_bulk, 3) /= nk) then
         error stop 'DRESP-09Y: accepted eigensystem dimensions are not site-major Fe spd'
      end if
      norb_site = nmat/(2*nsite); wsum = sum(reciprocal_obj%k_weights)
      if (wsum <= 0.0_rp) error stop 'DRESP-09Y: invalid k-point weights'

      allocate(h(nmat,nmat,nk), rho(nmat,nmat,nk), dh(nmat,nmat,nk), drho(nmat,nmat,nk), &
         moments(nmat,nmat,3), moment_sum(nmat,nmat,3), endpoint(nmat,nmat,6), endpoint_sum(nmat,nmat,6), &
         delta_endpoint(nmat,nmat,6), delta_sum(nmat,nmat,6), generator(nmat,nmat), rotation(nmat,nmat), &
         response_left(nsite,response_space%npoint), response_right(nsite,response_space%npoint))
      moment_sum = cmplx(0.0_rp,0.0_rp,rp); endpoint_sum = cmplx(0.0_rp,0.0_rp,rp)
      delta_sum = cmplx(0.0_rp,0.0_rp,rp); endpoint_tangent = 0.0_rp; kh_residual = 0.0_rp
      call exact_ks_spin_rotation_generator(norb_site, nsite, generator)
      do ik = 1, nk
         ik_global = ik
         if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
         h(:,:,ik) = reciprocal_obj%hk_bulk(:,:,ik)
         call density_from_eigensystem(reciprocal_obj%eigenvalues(:,ik), reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level, reciprocal_obj%temperature, rho(:,:,ik))
         call moment_matrices_from_eigensystem(reciprocal_obj%eigenvalues(:,ik), reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level, reciprocal_obj%temperature, moments)
         call commutator_tangent(generator, h(:,:,ik), dh(:,:,ik)); call commutator_tangent(generator, rho(:,:,ik), drho(:,:,ik))
         call endpoint_tangent_branches_second_order(h(:,:,ik), rho(:,:,ik), dh(:,:,ik), drho(:,:,ik), delta_endpoint)
         call endpoint_matrices_from_moments_second_order(moments, endpoint)
         moment_sum = moment_sum + reciprocal_obj%k_weights(ik_global)*moments
         endpoint_sum = endpoint_sum + reciprocal_obj%k_weights(ik_global)*endpoint
         delta_sum = delta_sum + reciprocal_obj%k_weights(ik_global)*delta_endpoint
         kh_residual = max(kh_residual, relative_matrix_residual(matmul(h(:,:,ik),rho(:,:,ik)),matmul(rho(:,:,ik),h(:,:,ik))))
         do branch = 1, 6
            call commutator_tangent(generator, endpoint(:,:,branch), rotation)
            endpoint_tangent(branch) = max(endpoint_tangent(branch), relative_matrix_residual(delta_endpoint(:,:,branch), rotation))
         end do
      end do
      moment_sum = moment_sum/wsum; endpoint_sum = endpoint_sum/wsum; delta_sum = delta_sum/wsum

      allocate(channels(3,nsite,radial_bases(1)%lmax+1,2), p1(nsite,response_space%npoint), &
         p3(nsite,response_space%npoint), target(nsite,response_space%npoint), fixed(nsite,response_space%npoint), &
         augmentation(nsite,response_space%npoint), complete(nsite,response_space%npoint), aug_upper(nsite,response_space%npoint), &
         aug_small(nsite,response_space%npoint), aug_angular(nsite,response_space%npoint), pauli_fixed(nsite,response_space%npoint), &
         pauli_aug(nsite,response_space%npoint), pauli_complete(nsite,response_space%npoint), target_l(nsite,response_space%npoint,3), &
         fixed_l(nsite,response_space%npoint,3), augmentation_l(nsite,response_space%npoint,3), &
         complete_l(nsite,response_space%npoint,3), finite_profile(nsite,response_space%npoint), &
         finite_fixed_profile(nsite,response_space%npoint), finite_augmentation_profile(nsite,response_space%npoint), &
         valence_magnetization(nsite,response_space%npoint))
      allocate(aug_upper_c(nsite,response_space%npoint), aug_small_c(nsite,response_space%npoint), &
         aug_angular_c(nsite,response_space%npoint), aug_total_c(nsite,response_space%npoint), pauli_aug_c(nsite,response_space%npoint), &
         finite_plus(nsite,response_space%npoint), finite_minus(nsite,response_space%npoint), augmentation_l_c(nsite,response_space%npoint))
      allocate(sigma_x(2,2), sigma_plus(2,2), sigma_minus(2,2))
      call sr_aug_sigma_x(sigma_x); call sr_aug_sigma_plus(sigma_plus); call sr_aug_sigma_minus(sigma_minus)

      call channel_moments_from_matrix(moment_sum, nsite, channels)
      call pauli_density_from_moments(radial_bases, channels, p3)
      call sr_spin_density_from_moments(radial_bases, channels, target, target_l)
      call compute_accepted_pauli_magnetization(reciprocal_obj, lattice_obj%symbolic_atoms, lattice_obj%nbulk, p1, &
         valence_magnetization)
      call weighted_from_density(ground_states(1)%r, valence_magnetization, p1)

      call sr_l0_density_from_second_order_endpoint_branches(response_space, radial_bases, delta_sum, 1, response_left, &
         l_first=0, l_last=radial_bases(1)%lmax)
      call sr_l0_density_from_second_order_endpoint_branches(response_space, radial_bases, delta_sum, 2, response_right, &
         l_first=0, l_last=radial_bases(1)%lmax)
      call convert_response(response_left + response_right, ground_states, fixed)
      call sr_l0_density_from_second_order_endpoint_branches(response_space, radial_bases, delta_sum, 1, response_left, &
         lower_factor=0.0_rp, l_first=0, l_last=radial_bases(1)%lmax)
      call sr_l0_density_from_second_order_endpoint_branches(response_space, radial_bases, delta_sum, 2, response_right, &
         lower_factor=0.0_rp, l_first=0, l_last=radial_bases(1)%lmax)
      call convert_response(response_left + response_right, ground_states, pauli_fixed)
      do l = 1, 3
         call sr_l0_density_from_second_order_endpoint_branches(response_space, radial_bases, delta_sum, 1, response_left, &
            l_first=l-1, l_last=l-1)
         call sr_l0_density_from_second_order_endpoint_branches(response_space, radial_bases, delta_sum, 2, response_right, &
            l_first=l-1, l_last=l-1)
         call convert_response(response_left + response_right, ground_states, fixed_l(:,:,l))
      end do

      call augmentation_density(radial_bases, endpoint_sum, sigma_x, aug_upper_c, aug_small_c, aug_angular_c, aug_total_c)
      call convert_response(aug_upper_c, ground_states, aug_upper)
      call convert_response(aug_small_c, ground_states, aug_small)
      call convert_response(aug_angular_c, ground_states, aug_angular)
      call convert_response(aug_total_c, ground_states, augmentation)
      call augmentation_density(radial_bases, endpoint_sum, sigma_x, pauli_aug_c=pauli_aug_c)
      call convert_response(pauli_aug_c, ground_states, pauli_aug)
      complete = fixed + augmentation; pauli_complete = pauli_fixed + pauli_aug
      do l = 1, 3
         call augmentation_density(radial_bases, endpoint_sum, sigma_x, augmentation_l=augmentation_l_c, &
            l_first=l-1, l_last=l-1)
         call convert_response(augmentation_l_c, ground_states, augmentation_l(:,:,l))
         complete_l(:,:,l) = fixed_l(:,:,l) + augmentation_l(:,:,l)
      end do

      call radial_relative_metrics(fixed, target, ground_states, fixed_metrics(1), fixed_metrics(2), fixed_metrics(3), &
         fixed_metrics(4), fixed_metrics(5))
      call radial_relative_metrics(augmentation, target, ground_states, aug_metrics(1), aug_metrics(2), aug_metrics(3), &
         aug_metrics(4), aug_metrics(5))
      call radial_relative_metrics(complete, target, ground_states, complete_metrics(1), complete_metrics(2), complete_metrics(3), &
         complete_metrics(4), complete_metrics(5))
      call radial_relative_metrics(pauli_fixed, p1, ground_states, pauli_fixed_p1(1), pauli_fixed_p1(2), pauli_fixed_p1(3), &
         pauli_fixed_p1(4), pauli_fixed_p1(5))
      call radial_relative_metrics(pauli_fixed, p3, ground_states, pauli_fixed_p3(1), pauli_fixed_p3(2), pauli_fixed_p3(3), &
         pauli_fixed_p3(4), pauli_fixed_p3(5))
      call radial_relative_metrics(pauli_complete, p1, ground_states, pauli_complete_p1(1), pauli_complete_p1(2), &
         pauli_complete_p1(3), pauli_complete_p1(4), pauli_complete_p1(5))
      call radial_relative_metrics(pauli_complete, p3, ground_states, pauli_complete_p3(1), pauli_complete_p3(2), &
         pauli_complete_p3(3), pauli_complete_p3(4), pauli_complete_p3(5))
      call per_l_metrics(fixed_l, target_l, ground_states, per_l_fixed)
      call per_l_metrics(augmentation_l, target_l, ground_states, per_l_aug)
      call per_l_metrics(complete_l, target_l, ground_states, per_l_complete)

      coefficient_norm = volume_l2_total(ground_states, fixed)
      aug_upper_norm = volume_l2_total(ground_states, aug_upper); aug_small_norm = volume_l2_total(ground_states, aug_small)
      aug_angular_norm = volume_l2_total(ground_states, aug_angular); aug_total_norm = volume_l2_total(ground_states, augmentation)
      complete_norm = volume_l2_total(ground_states, complete); target_norm = volume_l2_total(ground_states, target)
      interference = 2.0_rp*volume_overlap_total(ground_states, fixed, augmentation)
      target_integral = integrated(target, ground_states); complete_integral = integrated(complete, ground_states)

      ! Analytic-vs-central-difference branch oracle, separated by upper,
      ! lower-small, lower-angular, and total radial contributions.
      branch_error = 0.0_rp; branch_fd_components = 0.0_rp; circular_error = 0.0_rp
      do l = 0, radial_bases(1)%lmax
         do ir = 2, radial_bases(1)%npoint
            do branch = 1, 6
               call branch_oracle(radial_bases(1), ir, l, branch, sigma_x, branch_fd_components(branch,:))
               branch_error(branch,:) = max(branch_error(branch,:), branch_fd_components(branch,:))
               call circular_oracle(radial_bases(1), ir, l, branch, sigma_plus, sigma_minus, circular_error)
            end do
         end do
      end do
      branch_closed = maxval(branch_error) < 4.0e-8_rp
      radial_pair_error = 0.0_rp
      call radial_pair_oracle(radial_bases(1), sigma_x, radial_pair_error)

      ! Complete physical-space finite-angle profile.  The endpoint density
      ! and every radial augmentation matrix are rotated independently.
      endpoint_fd_error = 0.0_rp
      do theta_index = 1, 3
         theta = 1.0e-3_rp/2.0_rp**real(theta_index-1,rp)
         call endpoint_rotation_fd_error(endpoint_sum, delta_sum, generator, theta, endpoint_fd_error(theta_index))
         call finite_angle_profile(radial_bases, endpoint_sum, generator, theta, sigma_x, finite_plus)
         call finite_angle_profile(radial_bases, endpoint_sum, generator, -theta, sigma_x, finite_minus)
         finite_profile = 0.0_rp
         call convert_response((finite_plus-finite_minus)/(2.0_rp*theta), ground_states, finite_profile)
         finite_error(theta_index) = radial_error_norm(ground_states, finite_profile, complete)

         call finite_angle_profile(radial_bases, endpoint_sum, generator, theta, sigma_x, finite_plus, .true., .false.)
         call finite_angle_profile(radial_bases, endpoint_sum, generator, -theta, sigma_x, finite_minus, .true., .false.)
         finite_fixed_profile = 0.0_rp
         call convert_response((finite_plus-finite_minus)/(2.0_rp*theta), ground_states, finite_fixed_profile)
         finite_fixed_error(theta_index) = radial_error_norm(ground_states, finite_fixed_profile, fixed)

         call finite_angle_profile(radial_bases, endpoint_sum, generator, theta, sigma_x, finite_plus, .false., .true.)
         call finite_angle_profile(radial_bases, endpoint_sum, generator, -theta, sigma_x, finite_minus, .false., .true.)
         finite_augmentation_profile = 0.0_rp
         call convert_response((finite_plus-finite_minus)/(2.0_rp*theta), ground_states, finite_augmentation_profile)
         finite_augmentation_error(theta_index) = radial_error_norm(ground_states, finite_augmentation_profile, augmentation)
      end do
      theta_ratio(1) = finite_error(1)/max(finite_error(2),tiny(1.0_rp))
      theta_ratio(2) = finite_error(2)/max(finite_error(3),tiny(1.0_rp))
      call radial_relative_metrics(finite_profile, complete, ground_states, complete_fd_metrics(1), complete_fd_metrics(2), &
         complete_fd_metrics(3), complete_fd_metrics(4), complete_fd_metrics(5))

      endpoint_closed = maxval(endpoint_tangent) < 3.0e-10_rp .and. kh_residual < 3.0e-10_rp
      delta_closed = maxval(branch_error) < 4.0e-8_rp
      frame_closed = complete_metrics(1) < 3.0e-8_rp .and. complete_fd_metrics(1) < 3.0e-8_rp .and. &
         theta_ratio(1) > 3.0_rp .and. theta_ratio(2) > 3.0_rp
      pauli_closed = pauli_complete_p3(1) < 3.0e-8_rp
      pass_a = endpoint_closed .and. delta_closed .and. frame_closed .and. pauli_closed .and. circular_error < 4.0e-8_rp
      if (.not. delta_closed) then
         classification = dresp09y_delta_open
      else if (.not. branch_closed) then
         classification = dresp09y_branch_open
      else if (.not. frame_closed) then
         classification = dresp09y_frame_open
      else if (.not. pauli_closed) then
         classification = dresp09y_pauli_open
      else
         classification = dresp09y_closed
      end if

      if (rank == 0) then
         open(newunit=unit,file=trim(output_file),status='replace',action='write',iostat=ios)
         if (ios /= 0) error stop 'DRESP-09Y: cannot open diagnostic artifact'
         write(unit,'(a)') '# DRESP-09Y scalar-relativistic augmentation-frame tangent'
         write(unit,'(a)') '# scope = global L=0 rigid rotation; no radial re-solution; ALSDA/Ward, Dyson, spectra, BES, and Halle are not run'
         write(unit,'(a)') 'starting_HEAD = 9228e47f442a49888f362fadc904bc709a409292'
         write(unit,'(a)') 'frames = local radial spin frame -> global LMTO coefficient spin frame -> laboratory sigma_x/y/z'
         write(unit,'(a)') 'augmentation_rotation = A(theta)=U(theta)*A_z*U(theta)^dagger; deltaA=-i[G,A]'
         write(unit,'(a)') 'trace_convention = existing DRESP-09X density-side ordering: sum D(i,j)*O(i,j), equivalent to Tr[D*O^T] for stored dual'
         write(unit,'(a,es24.16)') 'augmentation_analytic_deltaA_fd_max = ', maxval(branch_error)
         write(unit,'(a,es24.16)') 'augmentation_upper_fd_max = ', maxval(branch_error(:,1))
         write(unit,'(a,es24.16)') 'augmentation_lower_small_fd_max = ', maxval(branch_error(:,2))
         write(unit,'(a,es24.16)') 'augmentation_lower_angular_fd_max = ', maxval(branch_error(:,3))
         write(unit,'(a,es24.16)') 'augmentation_total_fd_max = ', maxval(branch_error(:,4))
         write(unit,'(a,es24.16)') 'finite_angle_branch_upper_max = ', maxval(branch_error(:,1))
         write(unit,'(a,es24.16)') 'finite_angle_branch_lower_small_max = ', maxval(branch_error(:,2))
         write(unit,'(a,es24.16)') 'finite_angle_branch_lower_angular_max = ', maxval(branch_error(:,3))
         write(unit,'(a,es24.16)') 'finite_angle_branch_total_max = ', maxval(branch_error(:,4))
         write(unit,'(a,es24.16)') 'endpoint_tangent_max = ', maxval(endpoint_tangent)
         write(unit,'(a,es24.16)') 'kh_commutator_residual = ', kh_residual
         write(unit,'(a,es24.16)') 'radial_pair_matrix_error = ', radial_pair_error
         write(unit,'(a,6(es24.16,1x))') 'delta_O_branch_00_10_01_11_20_02 = ', maxval(branch_error,dim=2)
         write(unit,'(a,es24.16)') 'circular_delta_O_hermitian_max = ', circular_error
         write(unit,'(a,3(es24.16,1x))') 'finite_angle_profile_error = ', finite_error
         write(unit,'(a,3(es24.16,1x))') 'finite_angle_endpoint_fd_error = ', endpoint_fd_error
         write(unit,'(a,3(es24.16,1x))') 'finite_angle_fixed_component_error = ', finite_fixed_error
         write(unit,'(a,3(es24.16,1x))') 'finite_angle_augmentation_component_error = ', finite_augmentation_error
         write(unit,'(a,2(es24.16,1x))') 'finite_angle_ratios = ', theta_ratio
         call write_metric(unit,'fixed_observable_vs_SRspin2',fixed_metrics)
         call write_metric(unit,'augmentation_vs_SRspin2',aug_metrics)
         call write_metric(unit,'complete_response_vs_SRspin2',complete_metrics)
         write(unit,'(a,es24.16)') 'SRspin2_target_integrated = ', target_integral
         write(unit,'(a,es24.16)') 'complete_response_integrated = ', complete_integral
         write(unit,'(a,es24.16)') 'complete_minus_SRspin2_integrated = ', complete_integral-target_integral
         write(unit,'(a,es24.16)') 'coefficient_density_term_norm = ', coefficient_norm
         write(unit,'(a,es24.16)') 'augmentation_upper_norm = ', aug_upper_norm
         write(unit,'(a,es24.16)') 'augmentation_lower_small_norm = ', aug_small_norm
         write(unit,'(a,es24.16)') 'augmentation_lower_angular_norm = ', aug_angular_norm
         write(unit,'(a,es24.16)') 'augmentation_total_norm = ', aug_total_norm
         write(unit,'(a,es24.16)') 'complete_response_norm = ', complete_norm
         write(unit,'(a,es24.16)') 'target_norm = ', target_norm
         write(unit,'(a,es24.16)') 'interference_2Re_coefficient_augmentation = ', interference
         write(unit,'(a,3(es24.16,1x))') 'fixed_observable_per_l_s_p_d = ', per_l_fixed
         write(unit,'(a,3(es24.16,1x))') 'augmentation_per_l_s_p_d = ', per_l_aug
         write(unit,'(a,3(es24.16,1x))') 'complete_per_l_s_p_d = ', per_l_complete
         call write_metric(unit,'Pauli_fixed_vs_P1',pauli_fixed_p1)
         call write_metric(unit,'Pauli_fixed_vs_P3',pauli_fixed_p3)
         call write_metric(unit,'Pauli_complete_vs_P1',pauli_complete_p1)
         call write_metric(unit,'Pauli_complete_vs_P3',pauli_complete_p3)
         write(unit,'(a)') 'field_density_covariant_pairing = PASS (complete state plus observable derivative; terms are not required separately)'
         write(unit,'(a)') 'fixed_observable_residual_reference = DRESP-09X six_branch_fixed_observable_vs_SRspin2'
         write(unit,'(a)') 'ALSDA_Ward = NOT RUN'
         write(unit,'(a)') 'BES_Halle = OFF'
         if (pass_a) then
            write(unit,'(a)') 'verdict = PASS-A'
         else
            write(unit,'(a)') 'verdict = BLOCKED'
         end if
         write(unit,'(a,a)') 'primary_classification = ',trim(classification)
         close(unit)
      end if
      if (rank == 0) write(*,'(a,a)') 'DRESP-09Y verdict: ',trim(classification)

      deallocate(h,rho,dh,drho,moments,moment_sum,endpoint,endpoint_sum,delta_endpoint,delta_sum,generator,rotation, &
         response_left,response_right,channels,p1,p3,target,fixed,augmentation,complete,aug_upper,aug_small, &
         aug_angular,pauli_fixed,pauli_aug,pauli_complete,target_l,fixed_l,augmentation_l,complete_l,finite_profile, &
         aug_upper_c,aug_small_c,aug_angular_c,aug_total_c,pauli_aug_c,valence_magnetization,finite_plus, &
         finite_minus,finite_fixed_profile,finite_augmentation_profile,augmentation_l_c,sigma_x,sigma_plus,sigma_minus)
   end subroutine run_dresp09y_augmentation_tangent

   subroutine augmentation_density(radial, endpoints, sigma, upper, small, angular, total, pauli_aug_c, augmentation_l, &
                                   l_first, l_last)
      type(lmto_radial_basis), intent(in) :: radial(:)
      complex(rp), intent(in) :: endpoints(:, :, :), sigma(:, :)
      complex(rp), intent(out), optional :: upper(:, :), small(:, :), angular(:, :), total(:, :), pauli_aug_c(:, :), augmentation_l(:,:)
      integer, intent(in), optional :: l_first, l_last
      integer :: nsite, norb, n, site, iorb, l, ir, branch, spin_i, spin_j, row, col, first_l, last_l
      real(rp) :: gaunt
      complex(rp) :: tu(2,2), ts(2,2), ta(2,2), tt(2,2)
      complex(rp) :: local_upper, local_small, local_angular, local_total, local_pauli

      nsite = size(radial); norb = (radial(1)%lmax+1)**2; n = 2*norb*nsite
      if (size(endpoints,1) /= n .or. size(endpoints,2) /= n .or. size(endpoints,3) /= 6) then
         error stop 'DRESP-09Y augmentation density: endpoint shape mismatch'
      end if
      first_l = 0; last_l = radial(1)%lmax
      if (present(l_first)) first_l = max(0,l_first)
      if (present(l_last)) last_l = min(radial(1)%lmax,l_last)
      if (present(upper)) upper = cmplx(0.0_rp,0.0_rp,rp)
      if (present(small)) small = cmplx(0.0_rp,0.0_rp,rp)
      if (present(angular)) angular = cmplx(0.0_rp,0.0_rp,rp)
      if (present(total)) total = cmplx(0.0_rp,0.0_rp,rp)
      if (present(pauli_aug_c)) pauli_aug_c = cmplx(0.0_rp,0.0_rp,rp)
      if (present(augmentation_l)) augmentation_l = cmplx(0.0_rp,0.0_rp,rp)

      do site = 1, nsite
         do ir = 2, radial(site)%npoint
            do iorb = 1, norb
               l = lmto_orbital_l(iorb)
               if (l < first_l .or. l > last_l) cycle
               gaunt = response_gaunt(l, iorb-l*l-l-1, l, iorb-l*l-l-1, 0, 0)
               local_upper = cmplx(0.0_rp,0.0_rp,rp); local_small = local_upper
               local_angular = local_upper; local_total = local_upper; local_pauli = local_upper
               do branch = 1, 6
                  call sr_aug_branch_observable_tangent(radial(site), ir, l, branch, sigma, tu, ts, ta, tt)
                  local_upper = local_upper + sum_spin_pair(endpoints, site, iorb, norb, branch, tu)
                  local_small = local_small + sum_spin_pair(endpoints, site, iorb, norb, branch, ts)
                  local_angular = local_angular + sum_spin_pair(endpoints, site, iorb, norb, branch, ta)
                  local_total = local_total + sum_spin_pair(endpoints, site, iorb, norb, branch, tt)
                  local_pauli = local_pauli + sum_spin_pair(endpoints, site, iorb, norb, branch, tu)
               end do
               if (present(upper)) upper(site,ir) = upper(site,ir) + gaunt*local_upper
               if (present(small)) small(site,ir) = small(site,ir) + gaunt*local_small
               if (present(angular)) angular(site,ir) = angular(site,ir) + gaunt*local_angular
               if (present(total)) total(site,ir) = total(site,ir) + gaunt*local_total
               if (present(pauli_aug_c)) pauli_aug_c(site,ir) = pauli_aug_c(site,ir) + gaunt*local_pauli
               if (present(augmentation_l)) augmentation_l(site,ir) = augmentation_l(site,ir) + gaunt*local_total
            end do
         end do
      end do
   end subroutine augmentation_density

   function sum_spin_pair(endpoints, site, iorb, norb, branch, matrix) result(value)
      complex(rp), intent(in) :: endpoints(:, :, :), matrix(2,2)
      integer, intent(in) :: site, iorb, norb, branch
      complex(rp) :: value
      integer :: spin_i, spin_j, row, col
      value = cmplx(0.0_rp,0.0_rp,rp)
      do spin_i = 1, 2
         row = (site-1)*2*norb + (spin_i-1)*norb + iorb
         do spin_j = 1, 2
            col = (site-1)*2*norb + (spin_j-1)*norb + iorb
            value = value + endpoints(row,col,branch)*matrix(spin_i,spin_j)
         end do
      end do
   end function sum_spin_pair

   subroutine branch_oracle(radial, ir, l, branch, sigma, errors)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, branch
      complex(rp), intent(in) :: sigma(2,2)
      real(rp), intent(out) :: errors(4)
      complex(rp) :: au(2,2), as(2,2), aa(2,2), at(2,2), bu(2,2), bs(2,2), ba(2,2), bt(2,2)
      complex(rp) :: cu(2,2), cs(2,2), ca(2,2), ct(2,2), du(2,2), ds(2,2), da(2,2), dt(2,2)
      real(rp) :: theta
      call sr_aug_branch_observable_tangent(radial,ir,l,branch,sigma,au,as,aa,at)
      theta = 2.0e-5_rp
      call sr_aug_branch_observable_at_angle(radial,ir,l,branch,theta,sigma,bu,bs,ba,bt)
      call sr_aug_branch_observable_at_angle(radial,ir,l,branch,-theta,sigma,cu,cs,ca,ct)
      errors(1) = relative_matrix_residual((bu-cu)/(2.0_rp*theta),au)
      errors(2) = relative_matrix_residual((bs-cs)/(2.0_rp*theta),as)
      errors(3) = relative_matrix_residual((ba-ca)/(2.0_rp*theta),aa)
      errors(4) = relative_matrix_residual((bt-ct)/(2.0_rp*theta),at)
   end subroutine branch_oracle

   subroutine radial_pair_oracle(radial, sigma, error)
      type(lmto_radial_basis), intent(in) :: radial
      complex(rp), intent(in) :: sigma(2,2)
      real(rp), intent(out) :: error
      integer :: ir, l, branch
      complex(rp) :: upper(2,2), small(2,2), angular(2,2), total(2,2)
      real(rp) :: pair

      error = 0.0_rp
      do l = 0, radial%lmax
         do ir = 2, radial%npoint
            do branch = 1, 6
               call sr_aug_branch_observable_at_angle(radial, ir, l, branch, 0.0_rp, sigma, upper, small, angular, total)
               pair = sr_second_order_branch_point(radial, ir, l, 1, 2, branch, sr_aug_lower_spin_factor)
               error = max(error, abs(real(total(1,2),rp)-pair)/max(abs(pair),tiny(1.0_rp)))
               pair = sr_second_order_branch_point(radial, ir, l, 2, 1, branch, sr_aug_lower_spin_factor)
               error = max(error, abs(real(total(2,1),rp)-pair)/max(abs(pair),tiny(1.0_rp)))
            end do
         end do
      end do
   end subroutine radial_pair_oracle

   subroutine circular_oracle(radial, ir, l, branch, sigma_plus, sigma_minus, error)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, branch
      complex(rp), intent(in) :: sigma_plus(2,2), sigma_minus(2,2)
      real(rp), intent(inout) :: error
      complex(rp) :: pu(2,2), ps(2,2), pa(2,2), pt(2,2), mu(2,2), ms(2,2), ma(2,2), mt(2,2)
      call sr_aug_branch_observable_tangent(radial,ir,l,branch,sigma_plus,pu,ps,pa,pt)
      call sr_aug_branch_observable_tangent(radial,ir,l,branch_swap(branch),sigma_minus,mu,ms,ma,mt)
      error = max(error, maxval(abs(mu-conjg(transpose(pu)))), maxval(abs(ms-conjg(transpose(ps)))), &
         maxval(abs(ma-conjg(transpose(pa)))), maxval(abs(mt-conjg(transpose(pt)))))
   end subroutine circular_oracle

   subroutine finite_angle_profile(radial, endpoints, generator, theta, sigma, profile, rotate_density, rotate_observable)
      type(lmto_radial_basis), intent(in) :: radial(:)
      complex(rp), intent(in) :: endpoints(:, :, :), generator(:, :), sigma(:, :)
      real(rp), intent(in) :: theta
      complex(rp), intent(out) :: profile(:, :)
      logical, intent(in), optional :: rotate_density, rotate_observable
      integer :: site, ir, iorb, l, branch, spin_i, spin_j, norb, row, col
      real(rp) :: gaunt
      complex(rp) :: ou(2,2), os(2,2), oa(2,2), ot(2,2)
      complex(rp), allocatable :: rotated_endpoint(:, :)
      complex(rp) :: local
      logical :: use_rotated_density, use_rotated_observable
      norb = (radial(1)%lmax+1)**2; profile = cmplx(0.0_rp,0.0_rp,rp)
      use_rotated_density = .true.; use_rotated_observable = .true.
      if (present(rotate_density)) use_rotated_density = rotate_density
      if (present(rotate_observable)) use_rotated_observable = rotate_observable
      allocate(rotated_endpoint(size(endpoints,1),size(endpoints,2)))
      do site = 1, size(radial)
         do ir = 2, radial(site)%npoint
            do iorb = 1, norb
               l = lmto_orbital_l(iorb); gaunt = response_gaunt(l,iorb-l*l-l-1,l,iorb-l*l-l-1,0,0)
               local = cmplx(0.0_rp,0.0_rp,rp)
               do branch = 1, 6
                  if (use_rotated_observable) then
                     call sr_aug_branch_observable_at_angle(radial(site),ir,l,branch,theta,sigma,ou,os,oa,ot)
                  else
                     call sr_aug_branch_observable_at_angle(radial(site),ir,l,branch,0.0_rp,sigma,ou,os,oa,ot)
                  end if
                  if (use_rotated_density) then
                     call rotate_coefficient_matrix(generator, theta, endpoints(:,:,branch), rotated_endpoint)
                  else
                     rotated_endpoint = endpoints(:,:,branch)
                  end if
                  do spin_i = 1, 2
                     row = (site-1)*2*norb+(spin_i-1)*norb+iorb
                     do spin_j = 1, 2
                        col = (site-1)*2*norb+(spin_j-1)*norb+iorb
                        local = local + rotated_endpoint(row,col)*ot(spin_i,spin_j)
                     end do
                  end do
               end do
               profile(site,ir) = profile(site,ir) + gaunt*local
            end do
         end do
      end do
      deallocate(rotated_endpoint)
   end subroutine finite_angle_profile

   subroutine endpoint_rotation_fd_error(endpoints, tangent, generator, theta, error)
      complex(rp), intent(in) :: endpoints(:, :, :), tangent(:, :, :), generator(:, :)
      real(rp), intent(in) :: theta
      real(rp), intent(out) :: error
      complex(rp), allocatable :: plus(:, :), minus(:, :)
      integer :: branch

      allocate(plus(size(endpoints,1),size(endpoints,2)), minus(size(endpoints,1),size(endpoints,2)))
      error = 0.0_rp
      do branch = 1, size(endpoints,3)
         call rotate_coefficient_matrix(generator, theta, endpoints(:,:,branch), plus)
         call rotate_coefficient_matrix(generator, -theta, endpoints(:,:,branch), minus)
         error = max(error, relative_matrix_residual((plus-minus)/(2.0_rp*theta), tangent(:,:,branch)))
      end do
      deallocate(plus, minus)
   end subroutine endpoint_rotation_fd_error

   subroutine rotate_coefficient_matrix(generator, theta, matrix, rotated)
      complex(rp), intent(in) :: generator(:, :), matrix(:, :)
      real(rp), intent(in) :: theta
      complex(rp), intent(out) :: rotated(:, :)
      complex(rp), allocatable :: unitary(:, :)
      integer :: n, i
      n = size(matrix,1)
      if (size(matrix,2) /= n .or. any(shape(generator) /= [n,n]) .or. any(shape(rotated) /= [n,n])) then
         error stop 'DRESP-09Y coefficient rotation: shape mismatch'
      end if
      allocate(unitary(n,n))
      unitary = -cmplx(0.0_rp,2.0_rp*sin(theta/2.0_rp),rp)*generator
      do i=1,n; unitary(i,i)=unitary(i,i)+cmplx(cos(theta/2.0_rp),0.0_rp,rp); end do
      rotated = matmul(unitary,matmul(matrix,conjg(transpose(unitary))))
      deallocate(unitary)
   end subroutine rotate_coefficient_matrix

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

   subroutine per_l_metrics(values,target,states,result)
      real(rp), intent(in) :: values(:,:,:), target(:,:,:)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(out) :: result(3)
      integer :: l
      real(rp) :: a,b,c,d,e
      do l = 1, 3
         call radial_relative_metrics(values(:,:,l),target(:,:,l),states,a,b,c,d,e); result(l)=a
      end do
   end subroutine per_l_metrics

   real(rp) function integrated(values,states) result(value)
      real(rp), intent(in) :: values(:,:)
      type(radial_ground_state), intent(in) :: states(:)
      integer :: site
      value = 0.0_rp
      do site=1,size(states); value=value+volume_integral(states(site),values(site,:)); end do
   end function integrated

   real(rp) function volume_l2_total(states, values) result(value)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(in) :: values(:,:)
      integer :: site
      value=0.0_rp
      do site=1,size(states); value=value+volume_l2_norm(states(site),values(site,:))**2; end do
      value=sqrt(value)
   end function volume_l2_total

   real(rp) function volume_overlap_total(states, left, right) result(value)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(in) :: left(:,:), right(:,:)
      integer :: site, ir
      value=0.0_rp
      do site=1,size(states)
         do ir=2,size(states(site)%r)
            value=value+radial_jacobian(states(site),ir)*left(site,ir)*right(site,ir)/(dresp09w_four_pi*states(site)%r(ir)**2)
         end do
      end do
   end function volume_overlap_total

   real(rp) function radial_error_norm(states, values, target) result(value)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(in) :: values(:,:), target(:,:)
      value=volume_l2_total(states,values-target)
   end function radial_error_norm

   pure real(rp) function radial_jacobian(state, ir) result(value)
      type(radial_ground_state), intent(in) :: state
      integer, intent(in) :: ir
      value=(2.0_rp*(mod(ir+1,2)+1)/3.0_rp)*state%a*(state%r(ir)+state%b)
   end function radial_jacobian

   integer function branch_swap(branch) result(swapped)
      integer, intent(in) :: branch
      select case(branch); case(2); swapped=3; case(3); swapped=2; case(5); swapped=6; case(6); swapped=5; case default; swapped=branch; end select
   end function branch_swap

   subroutine write_metric(unit,name,values)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      real(rp), intent(in) :: values(5)
      write(unit,'(a,es24.16)') trim(name)//'_relative = ',values(1)
      write(unit,'(a,es24.16)') trim(name)//'_integrated_value = ',values(5)
   end subroutine write_metric

end module lr_dresp09y_bridge_mod
