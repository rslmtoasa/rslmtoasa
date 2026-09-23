!------------------------------------------------------------------------------
! DRESP-12 finite-LMTO covariance decomposition of the fixed-basis Ward defect.
!
! This module is deliberately a post-processing diagnostic.  It constructs
! the rigid fixed-basis field insertion and the independent LMTO commutator,
! forms their connection difference, and sends that difference through the
! already-certified static Frechet/Pauli measurement path.  No kernel,
! response basis, denominator, correction, or dynamics is modified here.
!------------------------------------------------------------------------------
module lr_dresp12_covariance_bridge_mod

   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use mpi_mod, only: rank
   use basis_mod, only: nb
   use reciprocal_mod, only: reciprocal
   use lattice_mod, only: lattice
   use hamiltonian_mod, only: hamiltonian
   use radial_ground_state_mod, only: radial_ground_state
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use response_angular_basis_mod, only: response_angular_pi, response_gaunt
   use response_basis_mapping_mod, only: response_super_index, response_unflatten_superindex, response_flatten_superindex
   use lr_response_space_mod, only: response_space_layout, response_vector_norm, response_vector_inner_product
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, lmto_product_channel_plus, lmto_product_nbranch
   use lr_lmto_endpoint_branches_mod, only: lmto_product_apply_branch_action, lmto_product_branch_powers, &
      lmto_product_second_order_radial_branch
   use lr_full_spatial_alsda_mod, only: lr_full_spatial_source_components, lr_full_spatial_contract_components, &
      lr_full_spatial_npiece, lr_full_spatial_piece_total
   use lr_static_frechet_product_response_mod, only: lr_static_frechet_density
   use lr_exact_ks_ward_mod, only: exact_ks_spin_rotation_generator
   use lr_lmto_density_moment_tangent_mod, only: density_from_eigensystem, moment_matrices_from_eigensystem, &
      endpoint_matrices_from_moments_second_order, endpoint_tangent_branches_second_order, &
      endpoint_fixed_h_branches_second_order, commutator_tangent, relative_matrix_residual
   use lr_dresp08_native_mapping_mod, only: dresp08_build_native_realspace_tangent, dresp08_build_native_onsite_tangents, &
      dresp08_block_diagonal, dresp08_build_h2, dresp08_build_product_tangent, dresp08_relative_residual
   use lr_dresp09_field_insertion_mod, only: dresp09_compact_field_from_raw, dresp09_raw_field_from_compact
   use lr_dresp09u_bridge_mod, only: run_dresp09u_representation_tangent
   use lr_dresp09y_bridge_mod, only: run_dresp09y_augmentation_tangent, dresp09y_augmentation_density, &
      dresp09y_convert_response
   use lr_dresp09s_scalar_relativistic_mod, only: sr_l0_density_from_second_order_endpoint_branches
   use lr_sr_augmentation_tangent_mod, only: sr_aug_sigma_plus, sr_aug_sigma_minus
   use lr_dresp10f_mixed_ward_bridge_mod, only: build_l0_source, build_l0_target, fixed_pauli_measurement, &
      endpoint_pauli_measurement, compact_density_measurement, apply_static_denominator
   use lr_radial_observable_provenance_mod, only: channel_moments_from_matrix, pauli_density_from_moments, &
      sr_spin_density_from_moments, volume_integral, radial_relative_metrics, dresp09w_four_pi
   use lmto_magnetic_tangent_mod, only: lmto_bond_derivative, lmto_hhmag_to_spinor
   use math_mod, only: hcpx
   implicit none
   private

   character(len=*), parameter, public :: dresp12_pass_a = 'PASS-A'
   character(len=*), parameter, public :: dresp12_pass_b = 'PASS-B'
   character(len=*), parameter, public :: dresp12_blocked = 'BLOCKED'

   public :: run_dresp12_covariance

contains

   subroutine run_dresp12_covariance(output_file, response_space, radial_bases, ground_states, reciprocal_obj, &
                                     lattice_obj, hamiltonian_obj)
      character(len=*), intent(in) :: output_file
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(reciprocal), intent(inout) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj

      type(lmto_product_response_basis) :: product
      type(response_super_index) :: item
      integer :: nsite, nmat, nk, ndim, n, product_dim, norb_site, ik, ik_global, flat, site, ir, unit, ios
      integer :: branch, piece, i, j, l, m, row, col, lrow, lcol, mode
      real(rp) :: wsum, wk, gaunt, source_residual, frechet_residual, compact_residual
      real(rp) :: max_field_identity, rms_field_identity, max_element_identity, production_residual
      real(rp) :: norm_b, norm_conn, norm_cov, max_h_conn_ratio, max_h_conn_element
      real(rp) :: rms_h_conn_ratio, median_h_conn_ratio, min_h_conn_ratio
      real(rp) :: occupied_action, near_ef_action, all_action, max_band_distance
      real(rp) :: k_all_action, k_occupied_action, k_near_ef_action, k_max_distance
      real(rp) :: branch_sq(lmto_product_nbranch), branch_mismatch_sq(lmto_product_nbranch)
      real(rp) :: origin_sq(4), orbital_sq(3,3), left_sq, right_sq, endpoint_interference
      real(rp) :: conn_norm, endpoint_h_norm, observable_norm, fixed_norm, master_norm, master_relative, account_relative, min_k_ratio
      real(rp) :: master_max, master_integrated, compact_projection_residual
      real(rp) :: l0_residual, nonspherical_norm, nonspherical_residual, nonspherical_l4_fraction
      real(rp) :: fixed_reconstruction_residual, master_reconstruction_residual, dmg_reconstruction_residual
      real(rp) :: dmg_l0_raw_residual, dmg_l0_compact_residual
      real(rp) :: l0_mode_cosine, l0_parallel_coefficient, l0_orthogonal_fraction
      real(rp) :: mxc_profile_relative, mxc_profile_l2, mxc_profile_max, mxc_profile_max_relative, mxc_integrated_difference
      real(rp) :: integrated_n_up_minus_down, integrated_mxc, integrated_p3, integrated_physical_sr_spin, integrated_core_mxc
      real(rp) :: mxc_core_bookkeeping_residual
      real(rp) :: constraining_field_max_abs
      real(rp) :: overlap_conn, overlap_endpoint, overlap_obs, cosine_conn, cosine_endpoint, cosine_obs, angle_conn, angle_endpoint, angle_obs
      real(rp) :: parallel_fraction(5), orthogonal_fraction(5), pnorm
      real(rp) :: dresp11_reconstruct_residual
      real(rp) :: response_conn_norm, response_frechet_residual
      real(rp) :: dm_cov_residual, dm_cov_linearity_residual, dm_cov_direct_residual, complete_fixed_observable_residual
      real(rp) :: endpoint_branch_residual(6), endpoint_zero_residual(6), endpoint_covariant_residual(6), endpoint_h_branch_sq(6)
      real(rp) :: observable_upper_norm, observable_small_norm, observable_angular_norm
      real(rp) :: observable_total_norm, conn_ratio, max_k_ratio, median_k_ratio, rms_k_ratio
      real(rp) :: fixed_basis_goldstone_relative
      real(rp) :: conn_endpoint_norm, conn_endpoint_obs_norm, overlap_p3_conn, overlap_p3_endpoint, overlap_p3_obs
      real(rp) :: observable_pauli_upper_residual, observable_component_sum_residual
      real(rp) :: endpoint_hermitian_residual(6), radial_swap_residual(6)
      real(rp) :: endpoint_y_norm, endpoint_12_norm, endpoint_12_y_residual
      real(rp) :: delta_sum_frobenius
      real(rp) :: trace_y_residual, trace_12_residual, x_reconstruction_residual, y_reconstruction_residual
      real(rp) :: observable_circular_x_residual, observable_circular_y_residual
      real(rp) :: pauli_complete_y_residual, pauli_complete_12_residual
      real(rp) :: endpoint_y_plus_norm(6), endpoint_y_minus_norm(6), endpoint_y_x_norm(6), endpoint_12_branch_norm(6)
      real(rp) :: endpoint_branch_y_residual(6), endpoint_branch_12_residual(6)
      real(rp) :: circular_plus_norm, circular_minus_norm
      real(rp) :: fixed_l(0:4), master_l(0:4), dmg_l(0:4), target_l(0:4)
      real(rp) :: k_ratios(64), band_residuals(64), ef_distances(64)
      complex(rp), allocatable :: moments(:, :, :), moment_sum(:, :, :), endpoint(:, :, :), endpoint_sum(:, :, :)
      complex(rp), allocatable :: generator(:, :), h(:, :), rho(:, :), dh_cov(:, :), delta_rho(:, :), delta_rho_cov(:, :), zero_matrix(:, :)
      complex(rp), allocatable :: source_field(:), target_raw(:), response_b(:), response_conn(:), response_cov(:)
      complex(rp), allocatable :: conn_compact(:), b_compact(:), cov_compact(:), target_compact(:)
      complex(rp), allocatable :: transition_conn(:), direct_conn(:), direct_conn_raw(:), raw_back(:), master_back(:), master_compact(:), observable(:), master(:)
      complex(rp), allocatable :: components(:, :, :, :), operators(:, :, :), branch_action(:, :)
      complex(rp), allocatable :: d_ee(:, :, :, :), d_o_types(:, :, :), d_e_types(:, :, :)
      complex(rp), allocatable :: hfirst(:, :), overlap(:, :), enu(:, :), d_hfirst(:, :), d_overlap(:, :), d_enu(:, :)
      complex(rp), allocatable :: d_h2(:, :), term_enu(:, :), term_h(:, :), term_left(:, :), term_middle(:, :), term_right(:, :)
      complex(rp), allocatable :: conn_enu(:, :), conn_h(:, :), conn_overlap(:, :), d_b(:, :), d_conn(:, :)
      complex(rp), allocatable :: covariant_contact(:, :), obs_upper(:, :), obs_small(:, :), obs_angular(:, :), obs_total(:, :)
      complex(rp), allocatable :: obs_upper_raw(:), obs_small_raw(:), obs_angular_raw(:), obs_total_raw(:)
      complex(rp), allocatable :: frozen_dmg(:), reconstructed_dmg(:), dmg_raw(:)
      complex(rp), allocatable :: dmg_l0(:), recon_l0(:), dmg_l0_compact(:), recon_l0_compact(:)
      complex(rp), allocatable :: obs_svd(:)
      complex(rp), allocatable :: covariant_field(:, :)
      complex(rp), allocatable :: endpoint_complete(:, :, :), endpoint_fixed(:, :, :), endpoint_h(:, :, :), endpoint_h_zero(:, :, :)
      complex(rp), allocatable :: endpoint_cov_y(:, :, :)
      complex(rp), allocatable :: response_cov_complete(:), response_cov_fixed_endpoint(:), response_endpoint_h(:)
      complex(rp), allocatable :: response_endpoint_12(:), endpoint_12_weighted(:), endpoint_raw(:), &
         endpoint_12_raw(:), fixed_endpoint_raw(:)
      complex(rp), allocatable :: endpoint_h_raw(:), endpoint_h_compact(:)
      complex(rp), allocatable :: delta_sum(:, :, :), endpoint_branch_only(:, :, :)
      complex(rp), allocatable :: response_y_plus(:), response_y_minus(:), response_y_x(:), response_y_y(:), target_weighted(:)
      complex(rp), allocatable :: response_trace_x(:), response_trace_y(:), response_12_branch(:)
      complex(rp), allocatable :: observable_y(:), observable_plus(:), observable_minus(:)
      complex(rp), allocatable :: obs_plus_c(:, :), obs_minus_c(:, :), obs_x_c(:, :), obs_y_c(:, :)
      real(rp), allocatable :: p3(:, :), mxc_weighted(:, :), sr_spin(:, :), core_mxc_weighted(:, :), bxc(:, :), kxc(:, :), channels(:, :, :, :)
      real(rp), allocatable :: sorted_ratios(:)
      logical :: source_ok, compact_ok, identity_ok, response_ok, angular_ok, l0_ok, endpoint_ok, measurement_ok, conjugate_ok
      logical :: sidecar_u, sidecar_y
      character(len=128) :: classification, verdict
      character(len=256) :: sidecar_u_file, sidecar_y_file, k_file

      nsite = size(ground_states); ndim = response_space%ndim
      if (nsite < 1 .or. size(radial_bases) /= nsite .or. response_space%nsite /= nsite .or. lattice_obj%nrec /= nsite) then
         error stop 'DRESP-12: response/radial/lattice site dimensions are inconsistent'
      end if
      if (.not. hamiltonian_obj%hoh .or. hamiltonian_obj%ccor_2c .or. hamiltonian_obj%orb_pol .or. &
          hamiltonian_obj%hubbard_u_general_check .or. hamiltonian_obj%hubbard_v_check .or. hamiltonian_obj%local_axis) then
         error stop 'DRESP-12: certified second-order ham_only state required'
      end if
      if (.not. allocated(reciprocal_obj%hk_bulk) .or. .not. allocated(reciprocal_obj%eigenvalues) .or. &
          .not. allocated(reciprocal_obj%eigenvectors) .or. .not. allocated(reciprocal_obj%k_weights)) then
         error stop 'DRESP-12: accepted reciprocal eigensystem is incomplete'
      end if
      nmat = size(reciprocal_obj%hk_bulk,1); nk = size(reciprocal_obj%hk_bulk,3)
      if (nk /= 64 .or. nmat /= nb*nsite .or. response_space%response_lmax /= 4 .or. response_space%nchannel /= 1) then
         error stop 'DRESP-12: certified 64-k, 348-space Fe state is required'
      end if
      call product%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      product_dim = product%product_dimension; n = product_dim
      if (n /= 348) error stop 'DRESP-12: live product basis dimension is not 348'
      wsum = sum(reciprocal_obj%k_weights)

      allocate(moment_sum(nmat,nmat,3), moments(nmat,nmat,3), endpoint_sum(nmat,nmat,6), endpoint(nmat,nmat,6), &
         generator(nmat,nmat), p3(nsite,response_space%npoint), mxc_weighted(nsite,response_space%npoint), &
         sr_spin(nsite,response_space%npoint), core_mxc_weighted(nsite,response_space%npoint), bxc(nsite,response_space%npoint), &
         kxc(nsite,response_space%npoint), channels(3,nsite,radial_bases(1)%lmax+1,2), source_field(ndim), target_raw(ndim), &
         response_b(ndim), response_conn(ndim), response_cov(ndim), observable(ndim), master(ndim), master_back(ndim), &
         obs_upper_raw(ndim), obs_small_raw(ndim), obs_angular_raw(ndim), obs_total_raw(ndim), &
         conn_compact(n), b_compact(n), cov_compact(n), target_compact(n), obs_svd(n), frozen_dmg(n), &
         reconstructed_dmg(n), dmg_raw(ndim), dmg_l0(ndim), recon_l0(ndim), dmg_l0_compact(n), recon_l0_compact(n), master_compact(n), &
         response_cov_complete(ndim), response_cov_fixed_endpoint(ndim), &
         response_endpoint_h(ndim), endpoint_raw(ndim), response_endpoint_12(ndim), endpoint_12_weighted(ndim), &
         endpoint_12_raw(ndim), fixed_endpoint_raw(ndim), endpoint_h_raw(ndim), endpoint_h_compact(n), &
         direct_conn(n), direct_conn_raw(ndim), transition_conn(n), raw_back(ndim), response_y_plus(ndim), &
         response_y_minus(ndim), response_y_x(ndim), response_y_y(ndim), target_weighted(ndim), response_trace_x(ndim), response_trace_y(ndim), &
         response_12_branch(ndim), observable_y(ndim), observable_plus(ndim), observable_minus(ndim))
      allocate(obs_plus_c(nsite,response_space%npoint), obs_minus_c(nsite,response_space%npoint), &
         obs_x_c(nsite,response_space%npoint), obs_y_c(nsite,response_space%npoint))
      allocate(d_ee(nb,nb,size(hamiltonian_obj%ee,3),size(hamiltonian_obj%ee,4)), &
         d_o_types(size(hamiltonian_obj%obarm,1),size(hamiltonian_obj%obarm,2),size(hamiltonian_obj%obarm,3)), &
         d_e_types(size(hamiltonian_obj%enim,1),size(hamiltonian_obj%enim,2),size(hamiltonian_obj%enim,3)))
      allocate(h(nmat,nmat), rho(nmat,nmat), dh_cov(nmat,nmat), delta_rho(nmat,nmat), delta_rho_cov(nmat,nmat), &
         zero_matrix(nmat,nmat), hfirst(nmat,nmat), &
         overlap(nmat,nmat), enu(nmat,nmat), d_hfirst(nmat,nmat), d_overlap(nmat,nmat), d_enu(nmat,nmat), &
         d_h2(nmat,nmat), term_enu(nmat,nmat), term_h(nmat,nmat), term_left(nmat,nmat), term_middle(nmat,nmat), &
         term_right(nmat,nmat), conn_enu(nmat,nmat), conn_h(nmat,nmat), conn_overlap(nmat,nmat), d_b(nmat,nmat), &
         d_conn(nmat,nmat), components(nmat,nmat,lmto_product_nbranch,lr_full_spatial_npiece), &
         operators(nmat,nmat,lr_full_spatial_npiece), branch_action(nmat,nmat), covariant_field(nmat,nmat), &
         endpoint_complete(nmat,nmat,lmto_product_nbranch), endpoint_fixed(nmat,nmat,lmto_product_nbranch), &
         endpoint_h(nmat,nmat,lmto_product_nbranch), endpoint_h_zero(nmat,nmat,lmto_product_nbranch), &
         endpoint_cov_y(nmat,nmat,lmto_product_nbranch), &
         delta_sum(nmat,nmat,lmto_product_nbranch), endpoint_branch_only(nmat,nmat,lmto_product_nbranch))

      moment_sum = cmplx(0.0_rp,0.0_rp,rp); endpoint_sum = cmplx(0.0_rp,0.0_rp,rp)
      delta_sum = cmplx(0.0_rp,0.0_rp,rp)
      call exact_ks_spin_rotation_generator((radial_bases(1)%lmax+1)**2,nsite,generator)
      do ik=1,nk
         ik_global=ik; if (allocated(reciprocal_obj%k_l2g_map)) ik_global=reciprocal_obj%k_l2g_map(ik)
         call moment_matrices_from_eigensystem(reciprocal_obj%eigenvalues(:,ik),reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level,reciprocal_obj%temperature,moments)
         call endpoint_matrices_from_moments_second_order(moments,endpoint)
         moment_sum=moment_sum+reciprocal_obj%k_weights(ik_global)*moments
         endpoint_sum=endpoint_sum+reciprocal_obj%k_weights(ik_global)*endpoint
      end do
      moment_sum=moment_sum/wsum; endpoint_sum=endpoint_sum/wsum
      call channel_moments_from_matrix(moment_sum,nsite,channels)
      call pauli_density_from_moments(radial_bases,channels,p3)
      call sr_spin_density_from_moments(radial_bases,channels,sr_spin)
      integrated_n_up_minus_down=0.0_rp; integrated_mxc=0.0_rp; integrated_p3=0.0_rp
      integrated_physical_sr_spin=0.0_rp; integrated_core_mxc=0.0_rp
      mxc_weighted=0.0_rp; core_mxc_weighted=0.0_rp
      do site=1,nsite
         if (.not. ground_states(site)%core_density_valid) then
            error stop 'DRESP-12: accepted radial state is missing frozen-core provenance'
         end if
         mxc_weighted(site,:)=ground_states(site)%rho_weighted_up-ground_states(site)%rho_weighted_down
         core_mxc_weighted(site,:)=ground_states(site)%core_weighted_up-ground_states(site)%core_weighted_down
         integrated_n_up_minus_down=integrated_n_up_minus_down+ground_states(site)%integrated_spin_number
         integrated_mxc=integrated_mxc+volume_integral(ground_states(site),mxc_weighted(site,:))
         integrated_p3=integrated_p3+volume_integral(ground_states(site),p3(site,:))
         integrated_physical_sr_spin=integrated_physical_sr_spin+volume_integral(ground_states(site),sr_spin(site,:))
         integrated_core_mxc=integrated_core_mxc+volume_integral(ground_states(site),core_mxc_weighted(site,:))
      end do
      call radial_relative_metrics(mxc_weighted,p3,ground_states,mxc_profile_relative,mxc_profile_l2,mxc_profile_max, &
         mxc_profile_max_relative,mxc_integrated_difference)
      mxc_integrated_difference=integrated_mxc-integrated_p3
      mxc_core_bookkeeping_residual=integrated_mxc-(integrated_mxc-integrated_core_mxc)-integrated_core_mxc
      constraining_field_max_abs=maxval(abs([(ground_states(site)%constraining_field_ry,site=1,nsite)]))
      do site=1,nsite
         do ir=1,response_space%npoint
            bxc(site,ir)=0.5_rp*(ground_states(site)%vxc_up(ir)-ground_states(site)%vxc_down(ir))
            if (response_space%radial_weights(ir)>0.0_rp) then
               if (p3(site,ir)==0.0_rp) error stop 'DRESP-12: zero P3 on positive-measure point'
               kxc(site,ir)=bxc(site,ir)/p3(site,ir)
            else
               kxc(site,ir)=0.0_rp
            end if
         end do
      end do
      call build_l0_source(response_space,bxc,source_field)
      call build_l0_target(response_space,p3,target_raw)
      call dresp09_compact_field_from_raw(response_space,product,target_raw,target_compact)
      call lr_full_spatial_source_components(response_space,radial_bases,source_field,1,components)
      call dresp08_build_native_realspace_tangent(hamiltonian_obj,d_ee)
      call dresp08_build_native_onsite_tangents(hamiltonian_obj,d_o_types,d_e_types)

      response_b=cmplx(0.0_rp,0.0_rp,rp); response_conn=cmplx(0.0_rp,0.0_rp,rp)
      response_cov=cmplx(0.0_rp,0.0_rp,rp); direct_conn=cmplx(0.0_rp,0.0_rp,rp)
      response_cov_complete=cmplx(0.0_rp,0.0_rp,rp); response_cov_fixed_endpoint=cmplx(0.0_rp,0.0_rp,rp)
      response_endpoint_h=cmplx(0.0_rp,0.0_rp,rp)
      response_endpoint_12=cmplx(0.0_rp,0.0_rp,rp)
      conn_compact=cmplx(0.0_rp,0.0_rp,rp); b_compact=cmplx(0.0_rp,0.0_rp,rp)
      cov_compact=cmplx(0.0_rp,0.0_rp,rp); branch_sq=0.0_rp; branch_mismatch_sq=0.0_rp
      endpoint_branch_residual=0.0_rp; endpoint_zero_residual=0.0_rp; endpoint_covariant_residual=0.0_rp
      endpoint_h_branch_sq=0.0_rp
      origin_sq=0.0_rp; orbital_sq=0.0_rp; left_sq=0.0_rp; right_sq=0.0_rp; endpoint_interference=0.0_rp
      max_field_identity=0.0_rp; rms_field_identity=0.0_rp; max_element_identity=0.0_rp; production_residual=0.0_rp
      min_h_conn_ratio=huge(1.0_rp); max_h_conn_ratio=0.0_rp; rms_h_conn_ratio=0.0_rp
      occupied_action=0.0_rp; near_ef_action=0.0_rp; all_action=0.0_rp; max_band_distance=0.0_rp
      response_frechet_residual=0.0_rp
      do ik=1,nk
         ik_global=ik; if (allocated(reciprocal_obj%k_l2g_map)) ik_global=reciprocal_obj%k_l2g_map(ik)
         wk=reciprocal_obj%k_weights(ik_global)
         h=reciprocal_obj%hk_bulk(:,:,ik)
         call density_from_eigensystem(reciprocal_obj%eigenvalues(:,ik),reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level,reciprocal_obj%temperature,rho)
         call commutator_tangent(generator,h,dh_cov)
         call lr_full_spatial_contract_components(components,h,operators)
         d_b=operators(:,:,lr_full_spatial_piece_total); d_conn=dh_cov-d_b
         call lr_static_frechet_density(reciprocal_obj%eigenvalues(:,ik),reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level,reciprocal_obj%temperature,d_conn,delta_rho)
         call fixed_pauli_measurement(response_space,radial_bases,h,delta_rho,direct_conn_raw)
         call dresp09_compact_field_from_raw(response_space,product,direct_conn_raw,direct_conn)
         call compact_density_measurement(product,reciprocal_obj%eigenvalues(:,ik),reciprocal_obj%eigenvectors(:,:,ik), &
            delta_rho,reciprocal_obj%fermi_level,reciprocal_obj%temperature,transition_conn)
         response_frechet_residual=max(response_frechet_residual,relative_vector(direct_conn,transition_conn))
         call accumulate_response(direct_conn_raw,wk,response_conn)
         call fixed_response_for_field(response_space,radial_bases,product,reciprocal_obj,h,d_b,ik,response_b,b_compact)
         call fixed_response_for_field(response_space,radial_bases,product,reciprocal_obj,h,dh_cov,ik,response_cov,cov_compact)

         ! The complete endpoint tangent is the authoritative DRESP-09V
         ! product rule.  The frozen endpoint tangent varies only rho.  Their
         ! subtraction is the explicit endpoint-H response, and the
         ! delta_rho=0 call is an independent isolation oracle.
         call commutator_tangent(generator,rho,delta_rho_cov)
         call endpoint_tangent_branches_second_order(h,rho,dh_cov,delta_rho_cov,endpoint_cov_y)
         delta_sum=delta_sum+wk*endpoint_cov_y
         do branch=1,lmto_product_nbranch
            endpoint_hermitian_residual(branch)=max(endpoint_hermitian_residual(branch), &
               hermitian_matrix_residual(endpoint_cov_y(:,:,branch),endpoint_cov_y(:,:,branch_swap(branch))))
         end do
         call endpoint_pauli_measurement(response_space,radial_bases,endpoint_cov_y,endpoint_12_raw)
         call accumulate_response(endpoint_12_raw,wk,response_endpoint_12)
         ! Keep both endpoint calls on the covariant tangent.  Reusing the
         ! earlier d_conn response here would mix two different perturbations.
         call endpoint_tangent_branches_second_order(h,rho,dh_cov,delta_rho_cov,endpoint_complete)
         call endpoint_fixed_h_branches_second_order(h,delta_rho_cov,endpoint_fixed)
         endpoint_h=endpoint_complete-endpoint_fixed
         zero_matrix=cmplx(0.0_rp,0.0_rp,rp)
         call endpoint_tangent_branches_second_order(h,rho,dh_cov,zero_matrix,endpoint_h_zero)
         do branch=1,lmto_product_nbranch
            endpoint_branch_residual(branch)=max(endpoint_branch_residual(branch), &
               relative_matrix_residual(endpoint_complete(:,:,branch),endpoint_fixed(:,:,branch)+endpoint_h(:,:,branch)))
            endpoint_zero_residual(branch)=max(endpoint_zero_residual(branch), &
               relative_matrix_residual(endpoint_h(:,:,branch),endpoint_h_zero(:,:,branch)))
            endpoint_covariant_residual(branch)=max(endpoint_covariant_residual(branch), &
               relative_matrix_residual(endpoint_complete(:,:,branch),endpoint_cov_y(:,:,branch)))
            endpoint_h_branch_sq(branch)=endpoint_h_branch_sq(branch)+wk*sum(abs(endpoint_h(:,:,branch))**2)
         end do
         call endpoint_pauli_measurement(response_space,radial_bases,endpoint_complete,endpoint_raw)
         call endpoint_pauli_measurement(response_space,radial_bases,endpoint_fixed,fixed_endpoint_raw)
         call endpoint_pauli_measurement(response_space,radial_bases,endpoint_h,endpoint_h_raw)
         call accumulate_response(endpoint_raw,wk,response_cov_complete)
         call accumulate_response(fixed_endpoint_raw,wk,response_cov_fixed_endpoint)
         call accumulate_response(endpoint_h_raw,wk,response_endpoint_h)

         call build_matrix_decomposition(ik, hamiltonian_obj, lattice_obj, reciprocal_obj, d_ee, d_o_types, d_e_types, &
            hfirst,overlap,enu,d_hfirst,d_overlap,d_enu,d_h2,term_enu,term_h,term_left,term_middle,term_right, &
            conn_enu,conn_h,conn_overlap,d_conn,production_residual)
         origin_sq(1)=origin_sq(1)+wk*sum(abs(conn_enu)**2)
         origin_sq(2)=origin_sq(2)+wk*sum(abs(conn_h)**2)
         origin_sq(3)=origin_sq(3)+wk*sum(abs(conn_overlap)**2)
         origin_sq(4)=origin_sq(4)+wk*sum(abs(d_conn)**2)
         call add_orbital_blocks(d_conn,wk,orbital_sq,left_sq,right_sq,endpoint_interference)
         call field_identity_metrics(d_b,d_conn,dh_cov,wk,max_field_identity,rms_field_identity,max_element_identity, &
            norm_b,norm_conn,norm_cov)
         call branch_metrics(h,components,wk,d_b,d_conn,branch_sq,branch_mismatch_sq)
         max_h_conn_ratio=max(max_h_conn_ratio,sqrt(sum(abs(d_conn)**2))/max(sqrt(sum(abs(dh_cov)**2)),tiny(1.0_rp)))
         min_h_conn_ratio=min(min_h_conn_ratio,sqrt(sum(abs(d_conn)**2))/max(sqrt(sum(abs(dh_cov)**2)),tiny(1.0_rp)))
         rms_h_conn_ratio=rms_h_conn_ratio+wk*(sqrt(sum(abs(d_conn)**2))/max(sqrt(sum(abs(dh_cov)**2)),tiny(1.0_rp)))**2
         call band_actions(d_conn,reciprocal_obj%eigenvectors(:,:,ik),reciprocal_obj%eigenvalues(:,ik), &
            reciprocal_obj%fermi_level,k_all_action,k_occupied_action,k_near_ef_action,k_max_distance,band_residuals(ik),ef_distances(ik))
         if(k_all_action>all_action) max_band_distance=k_max_distance
         all_action=max(all_action,k_all_action); occupied_action=max(occupied_action,k_occupied_action)
         near_ef_action=max(near_ef_action,k_near_ef_action)
         k_ratios(ik)=sqrt(sum(abs(d_conn)**2))/max(sqrt(sum(abs(dh_cov)**2)),tiny(1.0_rp))
      end do
      response_b=response_b/wsum; response_conn=response_conn/wsum; response_cov=response_cov/wsum
      response_cov_complete=response_cov_complete/wsum; response_cov_fixed_endpoint=response_cov_fixed_endpoint/wsum
      response_endpoint_h=response_endpoint_h/wsum
      response_endpoint_12=response_endpoint_12/wsum
      delta_sum=delta_sum/wsum
      b_compact=b_compact/wsum; cov_compact=cov_compact/wsum
      call dresp09_compact_field_from_raw(response_space,product,response_b,b_compact)
      call dresp09_compact_field_from_raw(response_space,product,response_cov,cov_compact)
      call dresp09_compact_field_from_raw(response_space,product,response_endpoint_h,endpoint_h_compact)
      norm_b=sqrt(sum(abs(response_b)**2)); norm_conn=sqrt(sum(abs(response_conn)**2)); norm_cov=sqrt(sum(abs(response_cov)**2))
      rms_h_conn_ratio=sqrt(rms_h_conn_ratio/wsum); rms_field_identity=sqrt(rms_field_identity/wsum)
      rms_k_ratio=rms_h_conn_ratio; min_k_ratio=minval(k_ratios); max_k_ratio=maxval(k_ratios); median_k_ratio=median_value(k_ratios)

      allocate(covariant_contact(nsite,response_space%npoint),obs_upper(nsite,response_space%npoint),obs_small(nsite,response_space%npoint), &
         obs_angular(nsite,response_space%npoint),obs_total(nsite,response_space%npoint))
      call dresp09y_augmentation_density(radial_bases,endpoint_sum,contact_matrix_dummy(),upper=obs_upper,small=obs_small, &
         angular=obs_angular,total=obs_total,pauli_aug_c=covariant_contact)
      call dresp09y_augmentation_density(radial_bases,endpoint_sum,circular_sigma_plus(),pauli_aug_c=obs_plus_c)
      call dresp09y_augmentation_density(radial_bases,endpoint_sum,circular_sigma_minus(),pauli_aug_c=obs_minus_c)
      obs_x_c=obs_plus_c+obs_minus_c
      obs_y_c=cmplx(0.0_rp,1.0_rp,rp)*(obs_plus_c-obs_minus_c)
      call build_l0_raw(response_space,covariant_contact,observable)
      call build_l0_raw(response_space,obs_upper,obs_upper_raw)
      call build_l0_raw(response_space,obs_small,obs_small_raw)
      call build_l0_raw(response_space,obs_angular,obs_angular_raw)
      call build_l0_raw(response_space,obs_total,obs_total_raw)
      observable_upper_norm=response_space_norm(response_space,obs_upper_raw)
      observable_small_norm=response_space_norm(response_space,obs_small_raw)
      observable_angular_norm=response_space_norm(response_space,obs_angular_raw)
      observable_total_norm=response_space_norm(response_space,obs_total_raw)
      observable_pauli_upper_residual=relative_vector(observable,obs_upper_raw)
      observable_component_sum_residual=relative_vector(obs_total_raw,obs_upper_raw+obs_small_raw+obs_angular_raw)

      ! Frozen DRESP-09Y authority: evaluate the two circular density paths
      ! independently, apply its exact conversion, and only then reconstruct
      ! Cartesian x/y.  This is deliberately separate from the current
      ! factor-of-two/conjugated endpoint measurement above.
      call dresp09y_path_response(response_space,radial_bases,ground_states,delta_sum,response_y_plus,response_y_minus, &
         response_y_x)
      response_y_y=cmplx(0.0_rp,1.0_rp,rp)*(response_y_plus-response_y_minus)
      call coefficient_trace_oracle(response_space,radial_bases,ground_states,delta_sum,response_trace_x,response_trace_y)
      call radial_branch_swap_audit(radial_bases,radial_swap_residual)
      call converted_density_raw(response_space,radial_bases,ground_states,obs_plus_c,observable_plus)
      call converted_density_raw(response_space,radial_bases,ground_states,obs_minus_c,observable_minus)
      call converted_density_raw(response_space,radial_bases,ground_states,obs_x_c,observable_y)
      ! DRESP-09Y frozen conversion returns sqrt(4*pi)*r^2 times the
      ! physical L=0 density.  Put the current DRESP-12 measurement into
      ! that same weighted coordinate before comparing paths or norms.
      call physical_to_dresp09y_weighted(response_space,ground_states,response_endpoint_12,endpoint_12_weighted)
      call physical_to_dresp09y_weighted(response_space,ground_states,target_raw,target_weighted)
      endpoint_y_norm=dresp09y_weighted_norm(response_space,ground_states,response_y_x)
      delta_sum_frobenius=sqrt(sum(abs(delta_sum)**2))
      endpoint_12_norm=dresp09y_weighted_norm(response_space,ground_states,endpoint_12_weighted)
      endpoint_12_y_residual=dresp09y_weighted_relative(response_space,ground_states,endpoint_12_weighted,response_y_x)
      trace_y_residual=dresp09y_weighted_relative(response_space,ground_states,response_trace_x,response_y_x)
      trace_12_residual=dresp09y_weighted_relative(response_space,ground_states,response_trace_x,endpoint_12_weighted)
      x_reconstruction_residual=dresp09y_weighted_relative(response_space,ground_states,response_y_plus+response_y_minus,response_y_x)
      y_reconstruction_residual=dresp09y_weighted_norm(response_space,ground_states,response_y_y-response_trace_y)/ &
         max(dresp09y_weighted_norm(response_space,ground_states,response_y_x),tiny(1.0_rp))
      observable_circular_x_residual=relative_complex_array(obs_x_c,covariant_contact)
      observable_circular_y_residual=relative_complex_array(obs_y_c,cmplx(0.0_rp,1.0_rp,rp)*(obs_plus_c-obs_minus_c))
      pauli_complete_y_residual=dresp09y_weighted_relative(response_space,ground_states,response_y_x+observable_y,target_weighted)
      pauli_complete_12_residual=dresp09y_weighted_relative(response_space,ground_states,endpoint_12_weighted+observable_y,target_weighted)
      circular_plus_norm=dresp09y_weighted_norm(response_space,ground_states,response_y_plus)
      circular_minus_norm=dresp09y_weighted_norm(response_space,ground_states,response_y_minus)

      endpoint_y_plus_norm=0.0_rp; endpoint_y_minus_norm=0.0_rp; endpoint_y_x_norm=0.0_rp
      endpoint_12_branch_norm=0.0_rp; endpoint_branch_y_residual=0.0_rp; endpoint_branch_12_residual=0.0_rp
      do branch=1,lmto_product_nbranch
         endpoint_branch_only=cmplx(0.0_rp,0.0_rp,rp)
         endpoint_branch_only(:,:,branch)=delta_sum(:,:,branch)
         call dresp09y_path_response(response_space,radial_bases,ground_states,endpoint_branch_only, &
            response_y_plus,response_y_minus,response_y_x)
         call endpoint_pauli_measurement(response_space,radial_bases,endpoint_branch_only,response_12_branch)
         call physical_to_dresp09y_weighted(response_space,ground_states,response_12_branch,endpoint_12_weighted)
         endpoint_y_plus_norm(branch)=dresp09y_weighted_norm(response_space,ground_states,response_y_plus)
         endpoint_y_minus_norm(branch)=dresp09y_weighted_norm(response_space,ground_states,response_y_minus)
         endpoint_y_x_norm(branch)=dresp09y_weighted_norm(response_space,ground_states,response_y_x)
         endpoint_12_branch_norm(branch)=dresp09y_weighted_norm(response_space,ground_states,endpoint_12_weighted)
         endpoint_branch_y_residual(branch)=dresp09y_weighted_relative(response_space,ground_states,response_y_plus+response_y_minus,response_y_x)
         endpoint_branch_12_residual(branch)=dresp09y_weighted_relative(response_space,ground_states,endpoint_12_weighted,response_y_x)
      end do

      fixed_norm=response_space_norm(response_space,response_b-target_raw)
      fixed_basis_goldstone_relative=fixed_norm/max(response_space_norm(response_space,target_raw),tiny(1.0_rp))
      conn_norm=response_space_norm(response_space,response_conn)
      endpoint_h_norm=response_space_norm(response_space,response_endpoint_h)
      conn_endpoint_norm=response_space_norm(response_space,response_conn+response_endpoint_h)
      conn_endpoint_obs_norm=response_space_norm(response_space,response_conn+response_endpoint_h+observable)
      observable_norm=response_space_norm(response_space,observable)
      overlap_p3_conn=real(response_vector_inner_product(response_space,target_raw,response_conn),rp)/ &
         max(response_space_norm(response_space,target_raw)*conn_norm,tiny(1.0_rp))
      overlap_p3_endpoint=real(response_vector_inner_product(response_space,target_raw,response_endpoint_h),rp)/ &
         max(response_space_norm(response_space,target_raw)*endpoint_h_norm,tiny(1.0_rp))
      overlap_p3_obs=real(response_vector_inner_product(response_space,target_raw,observable),rp)/ &
         max(response_space_norm(response_space,target_raw)*observable_norm,tiny(1.0_rp))
      dm_cov_linearity_residual=response_space_norm(response_space,response_cov-response_b-response_conn)/ &
         max(response_space_norm(response_space,response_cov),tiny(1.0_rp))
      dm_cov_direct_residual=response_space_norm(response_space,response_cov_fixed_endpoint-response_cov)/ &
         max(response_space_norm(response_space,response_cov),tiny(1.0_rp))
      complete_fixed_observable_residual=response_space_norm(response_space,response_cov_complete- &
         response_cov_fixed_endpoint-response_endpoint_h)/max(response_space_norm(response_space,response_cov_complete),tiny(1.0_rp))

      master=response_b-target_raw+response_conn+response_endpoint_h+observable
      master_norm=response_space_norm(response_space,master)
      account_relative=master_norm/max(fixed_norm,tiny(1.0_rp)); master_relative=account_relative
      master_max=maxval(abs(master)); master_integrated=weighted_integral(response_space,master)
      call dresp09_compact_field_from_raw(response_space,product,master,master_compact)
      call dresp09_raw_field_from_compact(response_space,product,master_compact,master_back)
      compact_projection_residual=relative_vector(master_back,master)
      call norm_by_l(response_space,response_b-target_raw,fixed_l)
      call norm_by_l(response_space,master,master_l)
      call norm_by_l(response_space,target_raw,target_l)
      fixed_reconstruction_residual=abs(fixed_norm**2-sum(fixed_l**2))/max(fixed_norm**2,tiny(1.0_rp))
      master_reconstruction_residual=abs(master_norm**2-sum(master_l**2))/max(master_norm**2,tiny(1.0_rp))
      l0_residual=master_l(0)/max(target_l(0),tiny(1.0_rp))
      nonspherical_norm=sqrt(sum(master_l(1:4)**2))
      nonspherical_residual=nonspherical_norm/max(target_l(0),tiny(1.0_rp))
      nonspherical_l4_fraction=master_l(4)/max(nonspherical_norm,tiny(1.0_rp))
      call vector_geometry(response_space,response_b-target_raw,response_conn,response_endpoint_h,observable,target_raw,overlap_conn, &
         overlap_endpoint,overlap_obs,cosine_conn,cosine_endpoint,cosine_obs,angle_conn,angle_endpoint,angle_obs,parallel_fraction, &
         orthogonal_fraction)
      conn_ratio=norm_conn/max(norm_cov,tiny(1.0_rp))
      dm_cov_residual=response_space_norm(response_space,response_cov_complete+observable-target_raw)/ &
         max(response_space_norm(response_space,target_raw),tiny(1.0_rp))

      ! Use one live denominator action for D*m_G.  No generated DRESP-11
      ! artifact is consulted for these per-L values.
      call apply_static_denominator(response_space,radial_bases,product,reciprocal_obj,kxc,target_compact,frozen_dmg)
      call dresp09_raw_field_from_compact(response_space,product,frozen_dmg,dmg_raw)
      call norm_by_l(response_space,dmg_raw,dmg_l)
      dmg_reconstruction_residual=abs(response_space_norm(response_space,dmg_raw)**2-sum(dmg_l**2))/ &
         max(response_space_norm(response_space,dmg_raw)**2,tiny(1.0_rp))
      call dresp09_compact_field_from_raw(response_space,product,response_conn,conn_compact)
      call dresp09_compact_field_from_raw(response_space,product,observable,obs_svd)
      ! DRESP-11 stores D*m_G = m_G - A*m_G = -r_fixed.  Since the
      ! completed accounting is r_fixed + conn + endpoint-H + deltaO = 0,
      ! the independently reconstructed D*m_G is the positive sum below.
      reconstructed_dmg=conn_compact+endpoint_h_compact+obs_svd
      dresp11_reconstruct_residual=relative_vector(reconstructed_dmg,frozen_dmg)
      call project_l0_vector(response_space,dmg_raw,dmg_l0)
      call project_l0_vector(response_space,response_conn+response_endpoint_h+observable,recon_l0)
      call project_l0_compact(product,frozen_dmg,dmg_l0_compact)
      call project_l0_compact(product,reconstructed_dmg,recon_l0_compact)
      dmg_l0_raw_residual=response_space_norm(response_space,recon_l0-dmg_l0)/ &
         max(response_space_norm(response_space,dmg_l0),tiny(1.0_rp))
      dmg_l0_compact_residual=relative_vector(recon_l0_compact,dmg_l0_compact)
      call l0_mode_geometry(response_space,target_raw,dmg_l0,l0_mode_cosine,l0_parallel_coefficient,l0_orthogonal_fraction)

      identity_ok=max_field_identity < 5.0e-11_rp .and. ieee_is_finite(production_residual)
      source_ok=ieee_is_finite(norm_b) .and. ieee_is_finite(norm_conn)
      compact_ok=response_frechet_residual < 1.0e-10_rp
      conjugate_ok=ieee_is_finite(mxc_profile_relative) .and. ieee_is_finite(mxc_integrated_difference) .and. &
         ieee_is_finite(integrated_physical_sr_spin) .and. ieee_is_finite(constraining_field_max_abs)
      response_ok=source_ok .and. compact_ok .and. conjugate_ok
      angular_ok=ieee_is_finite(fixed_reconstruction_residual) .and. ieee_is_finite(master_reconstruction_residual) .and. &
         max(fixed_reconstruction_residual,master_reconstruction_residual) < 1.0e-10_rp
      measurement_ok=maxval(endpoint_hermitian_residual) < 1.0e-10_rp .and. maxval(radial_swap_residual) < 1.0e-10_rp .and. &
         trace_y_residual < 1.0e-10_rp .and. x_reconstruction_residual < 1.0e-10_rp .and. &
         y_reconstruction_residual < 1.0e-10_rp .and. pauli_complete_y_residual < 3.0e-8_rp .and. &
         endpoint_12_y_residual < 1.0e-10_rp .and. trace_12_residual < 1.0e-10_rp .and. &
         pauli_complete_12_residual < 3.0e-8_rp .and. observable_circular_x_residual < 1.0e-10_rp .and. &
         observable_circular_y_residual < 1.0e-10_rp
      endpoint_ok=maxval(endpoint_branch_residual) < 1.0e-11_rp .and. maxval(endpoint_zero_residual) < 1.0e-11_rp .and. &
         maxval(endpoint_covariant_residual) < 1.0e-11_rp .and. &
         dm_cov_linearity_residual < 1.0e-10_rp .and. complete_fixed_observable_residual < 1.0e-10_rp .and. &
         dm_cov_direct_residual < 1.0e-10_rp .and. observable_pauli_upper_residual < 1.0e-10_rp
      ! Gate the direct raw-space L0 accounting; keep the compact
      ! projection/compression residual as a separate diagnostic.
      l0_ok=ieee_is_finite(l0_residual) .and. abs(l0_residual) < 1.0e-8_rp .and. &
         ieee_is_finite(dmg_l0_raw_residual) .and. dmg_l0_raw_residual < 1.0e-6_rp
      if (.not. endpoint_ok) then
         verdict=dresp12_blocked; classification='ENDPOINT_RESPONSE_REGRESSION'
      else if (endpoint_12_y_residual > 1.0e-8_rp .or. pauli_complete_12_residual > 1.0e-8_rp .or. &
         trace_12_residual > 1.0e-8_rp) then
         verdict=dresp12_blocked; classification='PAULI_MEASUREMENT_CONVENTION_MISMATCH'
      else if (.not. measurement_ok) then
         verdict=dresp12_blocked; classification='DRESP09Y_PROVENANCE_MISMATCH'
      else if (.not. angular_ok) then
         verdict=dresp12_blocked; classification='ANGULAR_DECOMPOSITION_INCONSISTENT'
      else if (.not. identity_ok .or. .not. response_ok) then
         verdict=dresp12_blocked; classification='COVARIANCE_BRIDGE_OPEN'
      else if (l0_ok) then
         verdict=dresp12_pass_a; classification='ASA_L0_RIGID_RESPONSE_CLOSED'
      else
         verdict=dresp12_blocked; classification='L0_ACCOUNTING_INCONSISTENT'
      end if

      sidecar_u_file=trim(output_file)//'.DRESP09U'; sidecar_y_file=trim(output_file)//'.DRESP09Y'
      call run_dresp09u_representation_tangent(sidecar_u_file,reciprocal_obj,lattice_obj,hamiltonian_obj)
      call run_dresp09y_augmentation_tangent(sidecar_y_file,response_space,radial_bases,ground_states,reciprocal_obj,lattice_obj,hamiltonian_obj)
      sidecar_u=.true.; sidecar_y=.true.; k_file=trim(output_file)//'.kpoints.csv'
      call write_k_table(k_file,k_ratios,band_residuals,ef_distances)
      if (rank==0) then
         open(newunit=unit,file=trim(output_file),status='replace',action='write',iostat=ios)
         if (ios/=0) error stop 'DRESP-12: cannot open output artifact'
         write(unit,'(a)') '# DRESP-12 finite-LMTO covariance decomposition of the Ward defect'
         write(unit,'(a,a)') 'DRESP-12 verdict: ',trim(verdict)
         write(unit,'(a)') 'Starting HEAD: ec3c6db1929906de6ef73b05c2b645a8be25f1b2'
         write(unit,'(a)') 'production_response_channel = chi_plus'
         write(unit,'(a)') 'production_circular_convention = independent_up_to_down_and_down_to_up_channels'
         write(unit,'(a)') 'cartesian_covariance_convention = m_x=m_plus+m_minus; m_y=i*(m_plus-m_minus)'
         write(unit,'(a)') 'factor_of_two_convention = current_DRESP12_single_block_factor_two_is_not_assumed_equivalent'
         write(unit,'(a)') 'branch_swap_convention = 00<->00, 10<->01, 01<->10, 11<->11, 20<->02, 02<->20'
         write(unit,'(a)') 'frozen_basis = six-branch 348-dimensional Pauli response; Frechet; Bxc/P3 historical representation diagnostic; raw SR; P3 target'
         write(unit,'(a)') 'fixed_formulation = m -> Kxc -> Bxc -> F_SR -> L_f -> R_P'
         write(unit,'(a)') 'covariant_formulation = m -> deltaH_cov -> L_f -> R_P plus endpoint-H plus deltaO'
         write(unit,'(a)') 'FINITE_LMTO_CONNECTION_TANGENT = deltaH_cov - deltaH_B; no fit or correction'
         write(unit,'(a)') 'conceptual_response = delta_m_LMTO = delta_m_Kubo(delta_rho) + delta_m_endpoint-H + delta_m_basis_observable(deltaO)'
         write(unit,'(a)') 'DRESP-11 regression = LIVE_MATRIX_FREE_DENOMINATOR_ACTION'
         write(unit,'(a,i0)') 'Pauli compact dimension = ',n
         write(unit,'(a)') 'branches = 00,10,01,11,20,02'
         write(unit,'(a,a)') 'XC functional = ',trim(ground_states(1)%xc_provenance%functional_name)
         write(unit,'(a,a)') 'XC backend = ',trim(ground_states(1)%xc_provenance%backend_name)
         write(unit,'(a,a,a,i0)') 'XC provenance = ',trim(ground_states(1)%xc_provenance%txch),' / TXC=',ground_states(1)%xc_provenance%txc
         write(unit,'(a)') 'm_xc definition = n_up-n_down from the live VXC0SP spherical density'
         write(unit,'(a)') 'm_xc radial comparison coordinate = weighted 4*pi*r^2*(n_up-n_down), matching P3'
         write(unit,'(a,es24.16)') 'integrated_n_up_minus_down = ',integrated_n_up_minus_down
         write(unit,'(a,es24.16)') 'integrated_m_xc_common_radial_metric = ',integrated_mxc
         write(unit,'(a,es24.16)') 'integrated_P3 = ',integrated_p3
         write(unit,'(a,es24.16)') 'integrated_physical_SR_spin = ',integrated_physical_sr_spin
         write(unit,'(a,es24.16)') 'm_xc_vs_P3_weighted_profile_difference = ',mxc_profile_relative
         write(unit,'(a,es24.16)') 'm_xc_vs_P3_weighted_L2_difference = ',mxc_profile_l2
         write(unit,'(a,es24.16)') 'm_xc_vs_P3_integrated_difference = ',mxc_integrated_difference
         write(unit,'(a,es24.16)') 'm_xc_vs_P3_max_physical_difference = ',mxc_profile_max
         write(unit,'(a,es24.16)') 'm_xc_vs_P3_max_relative_difference = ',mxc_profile_max_relative
         write(unit,'(a)') 'm_xc_vs_P3_sign_convention = PASS_UP_MINUS_DOWN'
         write(unit,'(a)') 'core_inclusion = m_xc includes captured frozen core plus valence; P3 is accepted reciprocal valence large-component density'
         write(unit,'(a,es24.16)') 'integrated_core_m_xc = ',integrated_core_mxc
         write(unit,'(a,es24.16)') 'm_xc_core_bookkeeping_residual = ',mxc_core_bookkeeping_residual
         write(unit,'(a,es24.16)') 'constraining_field_ry = ',ground_states(1)%constraining_field_ry
         write(unit,'(a,es24.16)') 'constraining_field_max_abs_ry = ',constraining_field_max_abs
         write(unit,'(a)') 'bxc_pauli_definition = 0.5*(vxc_up-vxc_down); XC field only; constraining field excluded'
         write(unit,'(a,es24.16)') 'fixed_Ward_norm = ',fixed_norm
         write(unit,'(a,es24.16)') '||r_fixed|| = ',fixed_norm
         write(unit,'(a,es24.16)') 'fixed_basis_goldstone_relative = ',fixed_basis_goldstone_relative
         write(unit,'(a,es24.16)') 'fixed_Ward_parallel_fraction = ',parallel_fraction(1)
         write(unit,'(a,es24.16)') 'fixed_Ward_orthogonal_fraction = ',orthogonal_fraction(1)
         write(unit,'(a,es24.16)') 'DRESP11_orthogonal_fraction = ',orthogonal_fraction(1)
         write(unit,'(a,es24.16)') '||deltaH_B|| = ',norm_b
         write(unit,'(a,es24.16)') '||deltaH_conn|| = ',norm_conn
         write(unit,'(a,es24.16)') '||deltaH_cov|| = ',norm_cov
         write(unit,'(a,es24.16)') 'connection_covariant_ratio = ',conn_ratio
         write(unit,'(a,es24.16)') 'field_identity_max_k_residual = ',max_field_identity
         write(unit,'(a,es24.16)') 'field_identity_weighted_RMS = ',rms_field_identity
         write(unit,'(a,es24.16)') 'field_identity_max_matrix_element = ',max_element_identity
         write(unit,'(a)') 'DRESP09U/Y_sidecars = CLOSED (independent frozen sidecar artifacts)'
         write(unit,'(a,es24.16)') 'DRESP08_native_product_tangent_vs_connection_residual = ',production_residual
         write(unit,'(a,es24.16)') 'band_action_occupied_max = ',occupied_action
         write(unit,'(a,es24.16)') 'band_action_near_EF_max = ',near_ef_action
         write(unit,'(a,es24.16)') 'band_action_all_max = ',all_action
         write(unit,'(a,es24.16)') 'band_action_dominant_distance_to_EF = ',max_band_distance
         write(unit,'(a)') 'Connection decomposition = Enu + h-minus-Bxc + overlap/hoh'
         write(unit,'(a,es24.16)') 'connection_Enu_norm = ',sqrt(origin_sq(1)/wsum)
         write(unit,'(a,es24.16)') 'connection_h_norm = ',sqrt(origin_sq(2)/wsum)
         write(unit,'(a,es24.16)') 'connection_overlap_hoh_norm = ',sqrt(origin_sq(3)/wsum)
         write(unit,'(a,es24.16)') 'connection_decomposition_sum_norm = ',sqrt(origin_sq(4)/wsum)
         write(unit,'(a,es24.16)') 'connection_ee_norm = ',sqrt(origin_sq(2)/wsum)
         write(unit,'(a)') 'connection_obarm_norm = reported within overlap/hoh; production tangent has no independent observable coordinate'
         write(unit,'(a)') 'connection_enim_norm = reported within Enu; production tangent has no independent duplicate coordinate'
         write(unit,'(a)') 'connection_left_right_endpoint = available only as certified DRESP-09U residuals; no inferred subtraction'
         write(unit,'(a,es24.16)') 'connection_endpoint_interference = ',endpoint_interference
         write(unit,'(a,6(es24.16,1x))') 'fixed_field_branch_norms_00_10_01_11_20_02 = ',sqrt(branch_sq/wsum)
         write(unit,'(a,6(es24.16,1x))') 'fixed_field_branch_fraction_00_10_01_11_20_02 = ',sqrt(branch_sq/wsum)/max(norm_b,tiny(1.0_rp))
         write(unit,'(a,6(es24.16,1x))') 'endpoint_branch_complete_frozen_endpoint_H_residual_00_10_01_11_20_02 = ',endpoint_branch_residual
         write(unit,'(a,6(es24.16,1x))') 'endpoint_H_subtraction_vs_delta_rho_zero_00_10_01_11_20_02 = ',endpoint_zero_residual
         write(unit,'(a,6(es24.16,1x))') 'endpoint_complete_vs_authoritative_covariant_00_10_01_11_20_02 = ',endpoint_covariant_residual
         write(unit,'(a,6(es24.16,1x))') 'endpoint_H_branch_norm_00_10_01_11_20_02 = ',sqrt(endpoint_h_branch_sq/wsum)
         write(unit,'(a,6(es24.16,1x))') 'endpoint_Hermitian_branch_swap_residual_00_10_01_11_20_02 = ',endpoint_hermitian_residual
         write(unit,'(a,6(es24.16,1x))') 'radial_spin_direction_branch_swap_residual_00_10_01_11_20_02 = ',radial_swap_residual
         write(unit,'(a)') 'endpoint_00_H = ZERO_BY_PRODUCT_RULE'
         write(unit,'(a)') 'covariant_branch_decomposition = CLOSED_BY_COMPLETE_ENDPOINT_TANGENT'
         write(unit,'(a)') 'connection_branch_mismatch_projection = NOT_USED'
         write(unit,'(a,3(es24.16,1x))') 'connection_orbital_norm_s_p_d = ',sqrt(orbital_sq(1,:)/wsum)
         write(unit,'(a,3(es24.16,1x))') 'connection_orbital_fraction_s_p_d = ',sqrt(orbital_sq(1,:)/wsum)/max(norm_conn,tiny(1.0_rp))
         write(unit,'(a,3(es24.16,1x))') 'connection_orbital_cross_sp_sd_pd = ',sqrt([orbital_sq(2,1),orbital_sq(2,2),orbital_sq(2,3)]/wsum)/max(norm_conn,tiny(1.0_rp))
         write(unit,'(a,4(es24.16,1x))') 'k_resolved_min_median_RMS_max = ',min_k_ratio,median_k_ratio,rms_k_ratio,max_k_ratio
         write(unit,'(a)') 'k_resolved_table = '//trim(k_file)
         write(unit,'(a,es24.16)') '||delta_m_conn|| = ',conn_norm
         write(unit,'(a,es24.16)') 'Frechet_oracle_residual = ',response_frechet_residual
         write(unit,'(a,es24.16)') 'independent_compact_connection_residual = ',response_frechet_residual
         write(unit,'(a,es24.16)') '||delta_m_cov_frozen|| = ',response_space_norm(response_space,response_cov_fixed_endpoint)
         write(unit,'(a,es24.16)') '||delta_m_cov_complete_fixedO|| = ',response_space_norm(response_space,response_cov_complete)
         write(unit,'(a,es24.16)') '||delta_m_endpoint_H|| = ',endpoint_h_norm
         write(unit,'(a,es24.16)') '||delta_m_conn_plus_endpoint_H|| = ',conn_endpoint_norm
         write(unit,'(a,es24.16)') '||delta_m_conn_plus_endpoint_H_plus_deltaO|| = ',conn_endpoint_obs_norm
         write(unit,'(a,es24.16)') 'dm_cov_rho_vs_dm_B_plus_dm_conn = ',dm_cov_linearity_residual
         write(unit,'(a,es24.16)') 'dm_cov_frozen_endpoint_vs_direct_fixed_H = ',dm_cov_direct_residual
         write(unit,'(a,es24.16)') 'dm_cov_complete_vs_frozen_plus_endpoint_H = ',complete_fixed_observable_residual
         write(unit,'(a,es24.16)') '||delta_m_O|| = ',observable_norm
         write(unit,'(a,es24.16)') 'delta_m_O_upper = ',observable_upper_norm
         write(unit,'(a,es24.16)') 'delta_m_O_lower_small = ',observable_small_norm
         write(unit,'(a,es24.16)') 'delta_m_O_lower_angular = ',observable_angular_norm
         write(unit,'(a,es24.16)') 'delta_m_O_total = ',observable_total_norm
         write(unit,'(a,es24.16)') 'DRESP09Y_observable_pauli_vs_upper_residual = ',observable_pauli_upper_residual
         write(unit,'(a,es24.16)') 'DRESP09Y_observable_upper_plus_lower_vs_total_residual = ',observable_component_sum_residual
         write(unit,'(a,es24.16)') 'DRESP09Y_endpoint_measurement_norm = ',endpoint_y_norm
         write(unit,'(a,es24.16)') 'DRESP09Y_delta_sum_frobenius = ',delta_sum_frobenius
         write(unit,'(a,es24.16)') 'DRESP12_endpoint_measurement_norm = ',endpoint_12_norm
         write(unit,'(a,es24.16)') 'DRESP12_vs_DRESP09Y_endpoint_residual = ',endpoint_12_y_residual
         write(unit,'(a,es24.16)') 'independent_cartesian_trace_oracle_vs_DRESP09Y_residual = ',trace_y_residual
         write(unit,'(a,es24.16)') 'independent_cartesian_trace_oracle_vs_DRESP12_residual = ',trace_12_residual
         write(unit,'(a,es24.16)') 'circular_plus_norm = ',circular_plus_norm
         write(unit,'(a,es24.16)') 'circular_minus_norm = ',circular_minus_norm
         write(unit,'(a,es24.16)') 'cartesian_x_reconstruction_residual = ',x_reconstruction_residual
         write(unit,'(a,es24.16)') 'cartesian_y_reconstruction_residual = ',y_reconstruction_residual
         write(unit,'(a,es24.16)') 'observable_deltaO_plus_minus_to_x_residual = ',observable_circular_x_residual
         write(unit,'(a,es24.16)') 'observable_deltaO_y_reconstruction_residual = ',observable_circular_y_residual
         write(unit,'(a,6(es24.16,1x))') 'DRESP09Y_channel1_branch_norm_00_10_01_11_20_02 = ',endpoint_y_plus_norm
         write(unit,'(a,6(es24.16,1x))') 'DRESP09Y_channel2_branch_norm_00_10_01_11_20_02 = ',endpoint_y_minus_norm
         write(unit,'(a,6(es24.16,1x))') 'DRESP12_branch_norm_00_10_01_11_20_02 = ',endpoint_12_branch_norm
         write(unit,'(a,6(es24.16,1x))') 'correct_circular_reconstruction_branch_norm_00_10_01_11_20_02 = ',endpoint_y_x_norm
         write(unit,'(a,6(es24.16,1x))') 'branch_circular_reconstruction_residual_00_10_01_11_20_02 = ',endpoint_branch_y_residual
         write(unit,'(a,6(es24.16,1x))') 'branch_DRESP12_vs_correct_residual_00_10_01_11_20_02 = ',endpoint_branch_12_residual
         write(unit,'(a,es24.16)') 'Embedded_DRESP09Y_Pauli_complete_vs_P3 = ',pauli_complete_y_residual
         write(unit,'(a,es24.16)') 'Current_DRESP12_Pauli_complete_vs_P3 = ',pauli_complete_12_residual
         write(unit,'(a,es24.16)') 'Fixed_basis_Goldstone_circular = ',fixed_basis_goldstone_relative
         write(unit,'(a,es24.16)') 'Fixed_basis_Goldstone_Cartesian_reconstructed = ',fixed_basis_goldstone_relative
         write(unit,'(a,es24.16)') 'Pauli_measurement_seam_residual = ',pauli_complete_y_residual
         write(unit,'(a)') 'Pauli_measurement_seam_note = circular and Cartesian paths independently agree'
         write(unit,'(a,es24.16)') 'covariance_accounting_relative_to_fixed_defect = ',account_relative
         write(unit,'(a,es24.16)') 'master_identity_relative = ',account_relative
         write(unit,'(a,es24.16)') 'historical_incomplete_accounting = ',1.4619_rp
         write(unit,'(a,es24.16)') 'master_identity_maximum_radial_component = ',master_max
         write(unit,'(a,es24.16)') 'master_identity_integrated_moment = ',master_integrated
         write(unit,'(a,es24.16)') 'master_identity_348_compact_projection_residual = ',compact_projection_residual
         write(unit,'(a,es24.16)') 'master_identity_covariant_response_residual = ',dm_cov_residual
         write(unit,'(a,es24.16)') 'l0_residual = ',l0_residual
         write(unit,'(a,es24.16)') 'nonspherical_norm = ',nonspherical_norm
         write(unit,'(a,es24.16)') 'nonspherical_residual = ',nonspherical_residual
         write(unit,'(a,es24.16)') 'nonspherical_L4_fraction = ',nonspherical_l4_fraction
         write(unit,'(a,5(es24.16,1x))') 'r_fixed_by_L_0_1_2_3_4 = ',fixed_l
         write(unit,'(a,es24.16)') 'r_fixed_full = ',fixed_norm
         write(unit,'(a,es24.16)') 'r_fixed_norm_reconstruction_residual = ',fixed_reconstruction_residual
         write(unit,'(a,5(es24.16,1x))') 'master_by_L_0_1_2_3_4 = ',master_l
         write(unit,'(a,es24.16)') 'master_full = ',master_norm
         write(unit,'(a,es24.16)') 'master_norm_reconstruction_residual = ',master_reconstruction_residual
         write(unit,'(a,4(es24.16,1x))') 'master_nonspherical_fraction_L1_L2_L3_L4 = ',master_l(1:4)/max(nonspherical_norm,tiny(1.0_rp))
         write(unit,'(a,5(es24.16,1x))') 'DmG_by_L_0_1_2_3_4 = ',dmg_l
         write(unit,'(a,es24.16)') 'DmG_full = ',response_space_norm(response_space,dmg_raw)
         write(unit,'(a,es24.16)') 'DmG_norm_reconstruction_residual = ',dmg_reconstruction_residual
         write(unit,'(a,es24.16)') 'circular_plus_accounting_residual = ',dresp09y_weighted_relative(response_space,ground_states,response_y_plus+observable_plus,target_weighted)
         write(unit,'(a,es24.16)') 'circular_minus_accounting_residual = ',dresp09y_weighted_relative(response_space,ground_states,response_y_minus+observable_minus,target_weighted)
         write(unit,'(a,2es24.16)') 'overlap_r_conn_complex_real = ',overlap_conn,cosine_conn
         write(unit,'(a,2es24.16)') 'overlap_r_deltaO_complex_real = ',overlap_obs,cosine_obs
         write(unit,'(a,3(es24.16,1x))') 'overlap_P3_conn_endpoint_H_deltaO = ',overlap_p3_conn,overlap_p3_endpoint,overlap_p3_obs
         write(unit,'(a,2es24.16)') 'angle_r_conn_radians_cosine = ',angle_conn,cosine_conn
         write(unit,'(a,2es24.16)') 'angle_r_deltaO_radians_cosine = ',angle_obs,cosine_obs
         write(unit,'(a,5(es24.16,1x))') 'parallel_fraction_r_conn_endpoint_H_deltaO_sum = ',parallel_fraction
         write(unit,'(a,5(es24.16,1x))') 'orthogonal_fraction_r_conn_endpoint_H_deltaO_sum = ',orthogonal_fraction
         write(unit,'(a,es24.16)') 'overlap_r_endpoint_H_complex_real = ',overlap_endpoint
         write(unit,'(a,2es24.16)') 'angle_r_endpoint_H_radians_cosine = ',angle_endpoint,cosine_endpoint
         write(unit,'(a,es24.16)') 'DRESP11_full_space_DmG_reconstruction_residual = ',dresp11_reconstruct_residual
         write(unit,'(a)') 'DRESP11_full_space_reconstruction_status = DIAGNOSTIC_ONLY_L4_MODEL_BOUNDARY; not a regression gate'
         write(unit,'(a,es24.16)') 'dmg_l0_raw_residual = ',dmg_l0_raw_residual
         write(unit,'(a,es24.16)') 'dmg_l0_compact_residual = ',dmg_l0_compact_residual
         write(unit,'(a)') 'dmg_l0_compact_residual_role = PROJECTION_COMPRESSION_DIAGNOSTIC'
         write(unit,'(a,es24.16)') 'l0_mode_cosine = ',l0_mode_cosine
         write(unit,'(a,es24.16)') 'l0_parallel_coefficient = ',l0_parallel_coefficient
         write(unit,'(a,es24.16)') 'l0_orthogonal_fraction = ',l0_orthogonal_fraction
         write(unit,'(a)') 'DRESP11_DmG_sign_convention = DmG = mG - A*mG = -r_fixed; L0 reconstruction is the regression gate'
         write(unit,'(a)') 'sitewise_rigid_rotation = DERIVABLE_FROM_EXISTING_PRODUCTION_MAP (two-site endpoint superposition fixture CLOSED)'
         ! The present ASA ground state is stationary only in spherical radial
         ! density space.  A nonspherical Goldstone statement needs a ground-
         ! state functional containing the corresponding nonspherical fields.
         write(unit,'(a)') 'arbitrary_L_within_ASA = NONSPHERICAL_GROUND_STATE_RESPONSE_NOT_DEFINED'
         write(unit,'(a)') 'multi_site_nonuniform = DERIVABLE_FOR_SITEWISE_RIGID_ROTATIONS_ONLY; arbitrary local field remains OPEN'
         write(unit,'(a)') 'sitewise_fixture = CLOSED_BY_LINEAR_ENDPOINT_SUPERPOSITION'
         write(unit,'(a)') 'Primary classification = '//trim(classification)
         if (endpoint_12_y_residual > 1.0e-8_rp .or. pauli_complete_12_residual > 1.0e-8_rp .or. &
            trace_12_residual > 1.0e-8_rp) then
            write(unit,'(a)') 'Measurement seam classification = PAULI_MEASUREMENT_CONVENTION_MISMATCH'
         else
            write(unit,'(a)') 'Measurement seam classification = MEASUREMENT_CONVENTIONS_IDENTICAL'
         end if
         write(unit,'(a)') 'Goldstone correction = OFF'
         write(unit,'(a)') 'BES/Halle production = OFF'
         write(unit,'(a)') 'Dynamics = NOT RUN'
         write(unit,'(a)') 'Model boundary = NONSPHERICAL_RESPONSE_ON_SPHERICAL_ASA_GROUND_STATE'
         write(unit,'(a)') 'Strict ASA transverse Goldstone space = L0'
         write(unit,'(a)') 'Ward-focused campaign = '//merge('CLOSED','OPEN  ',trim(verdict)==dresp12_pass_a)
         if (trim(verdict)==dresp12_pass_a) then
            write(unit,'(a)') 'NEXT = NATIVE_ROTATION_DYNAMICS'
         else
            write(unit,'(a)') 'NEXT = STOP_AND_DIAGNOSE'
         end if
         write(unit,'(a)') 'Tests = six-branch covariant tangent; authoritative endpoint regression; Frechet linearity; circular/Cartesian seam; per-L norm reconstruction; live DmG action; frozen DRESP-10F/DRESP-11; git diff --check'
         close(unit)
      end if
      if (rank==0) write(*,'(a,a)') 'DRESP-12 verdict: ',trim(verdict)

      deallocate(moment_sum,moments,endpoint_sum,endpoint,generator,p3,mxc_weighted,sr_spin,core_mxc_weighted,bxc,kxc,channels,source_field,target_raw,response_b, &
         response_conn,response_cov,observable,master,master_back,master_compact,obs_upper_raw,obs_small_raw,obs_angular_raw,obs_total_raw,conn_compact, &
         b_compact,cov_compact,target_compact,obs_svd,frozen_dmg,reconstructed_dmg,dmg_raw,dmg_l0,recon_l0,dmg_l0_compact,recon_l0_compact, &
         response_y_plus,response_y_minus,response_y_x,response_y_y,target_weighted,response_trace_x,response_trace_y,response_12_branch,endpoint_12_weighted, &
         observable_y,observable_plus,observable_minus,obs_plus_c,obs_minus_c,obs_x_c,obs_y_c,delta_sum,endpoint_branch_only, &
         response_cov_complete,response_cov_fixed_endpoint,response_endpoint_h,response_endpoint_12, &
         endpoint_raw,fixed_endpoint_raw,endpoint_h_raw,endpoint_h_compact, &
         d_ee,d_o_types,d_e_types,h,rho,dh_cov,delta_rho,delta_rho_cov,zero_matrix,hfirst,overlap,enu,d_hfirst,d_overlap,d_enu,d_h2,term_enu,term_h, &
         term_left,term_middle,term_right,conn_enu,conn_h,conn_overlap,d_b,d_conn,components,operators,branch_action,covariant_field, &
         endpoint_complete,endpoint_fixed,endpoint_h,endpoint_h_zero, &
         endpoint_cov_y,endpoint_12_raw, &
         covariant_contact,obs_upper,obs_small,obs_angular,obs_total,direct_conn,direct_conn_raw, &
         transition_conn,raw_back)
   contains

      function contact_matrix_dummy() result(sigma)
         complex(rp) :: sigma(2,2)
         sigma=cmplx(0.0_rp,0.0_rp,rp); sigma(1,2)=cmplx(1.0_rp,0.0_rp,rp); sigma(2,1)=sigma(1,2)
      end function contact_matrix_dummy

      function circular_sigma_plus() result(sigma)
         complex(rp) :: sigma(2,2)
         call sr_aug_sigma_plus(sigma)
      end function circular_sigma_plus

      function circular_sigma_minus() result(sigma)
         complex(rp) :: sigma(2,2)
         call sr_aug_sigma_minus(sigma)
      end function circular_sigma_minus

      subroutine dresp09y_path_response(space,radial,states,branches,plus,minus,x_direct)
         type(response_space_layout), intent(in) :: space
         type(lmto_radial_basis), intent(in) :: radial(:)
         type(radial_ground_state), intent(in) :: states(:)
         complex(rp), intent(in) :: branches(:, :, :)
         complex(rp), intent(out) :: plus(:), minus(:), x_direct(:)
         complex(rp), allocatable :: plus_density(:, :), minus_density(:, :), x_density(:, :)
         real(rp), allocatable :: plus_weighted(:, :), minus_weighted(:, :), x_weighted(:, :)

         allocate(plus_density(size(radial),space%npoint),minus_density(size(radial),space%npoint), &
            x_density(size(radial),space%npoint),plus_weighted(size(radial),space%npoint), &
            minus_weighted(size(radial),space%npoint),x_weighted(size(radial),space%npoint))
         call sr_l0_density_from_second_order_endpoint_branches(space,radial,branches,1,plus_density,lower_factor=0.0_rp)
         call sr_l0_density_from_second_order_endpoint_branches(space,radial,branches,2,minus_density,lower_factor=0.0_rp)
         x_density=plus_density+minus_density
         call dresp09y_convert_response(plus_density,states,plus_weighted)
         call dresp09y_convert_response(minus_density,states,minus_weighted)
         call dresp09y_convert_response(x_density,states,x_weighted)
         call build_l0_real_raw(space,plus_weighted,plus)
         call build_l0_real_raw(space,minus_weighted,minus)
         call build_l0_real_raw(space,x_weighted,x_direct)
         deallocate(plus_density,minus_density,x_density,plus_weighted,minus_weighted,x_weighted)
      end subroutine dresp09y_path_response

      subroutine converted_density_raw(space,radial,states,density,raw)
         type(response_space_layout), intent(in) :: space
         type(lmto_radial_basis), intent(in) :: radial(:)
         type(radial_ground_state), intent(in) :: states(:)
         complex(rp), intent(in) :: density(:, :)
         complex(rp), intent(out) :: raw(:)
         real(rp), allocatable :: weighted(:, :)

         allocate(weighted(size(states),space%npoint))
         call dresp09y_convert_response(density,states,weighted)
         call build_l0_real_raw(space,weighted,raw)
         deallocate(weighted)
      end subroutine converted_density_raw

      subroutine build_l0_real_raw(space,values,raw)
         type(response_space_layout), intent(in) :: space
         real(rp), intent(in) :: values(:, :)
         complex(rp), intent(out) :: raw(:)
         type(response_super_index) :: local_item
         integer :: flat

         raw=cmplx(0.0_rp,0.0_rp,rp)
         do flat=1,size(raw)
            call response_unflatten_superindex(flat,space%nsite,space%response_lmax,space%npoint,space%nchannel,local_item)
            if(local_item%response_l==0 .and. local_item%response_m==0 .and. local_item%radial_point>1) then
               raw(flat)=cmplx(values(local_item%site,local_item%radial_point),0.0_rp,rp)
            end if
         end do
      end subroutine build_l0_real_raw

      ! DRESP-09Y stores its converted radial response as
      ! sqrt(4*pi)*r^2 times the physical density.  The production response
      ! space, in contrast, stores the physical L=0 density.  These helpers
      ! make that seam explicit and give all provenance comparisons the same
      ! volume metric used by the frozen DRESP-09Y sidecar.
      subroutine physical_to_dresp09y_weighted(space,states,physical,weighted)
         type(response_space_layout), intent(in) :: space
         type(radial_ground_state), intent(in) :: states(:)
         complex(rp), intent(in) :: physical(:)
         complex(rp), intent(out) :: weighted(:)
         type(response_super_index) :: local_item
         integer :: flat, site, ir

         if (size(physical)/=space%ndim .or. size(weighted)/=space%ndim .or. size(states)/=space%nsite) then
            error stop 'DRESP-12: weighted seam vector shape mismatch'
         end if
         weighted=cmplx(0.0_rp,0.0_rp,rp)
         do flat=1,space%ndim
            call response_unflatten_superindex(flat,space%nsite,space%response_lmax,space%npoint,space%nchannel,local_item)
            site=local_item%site; ir=local_item%radial_point
            if (local_item%response_l==0 .and. local_item%response_m==0 .and. ir>1) then
               weighted(flat)=sqrt(dresp09w_four_pi)*states(site)%r(ir)**2*physical(flat)
            end if
         end do
      end subroutine physical_to_dresp09y_weighted

      subroutine dresp09y_weighted_to_physical(space,states,weighted,physical)
         type(response_space_layout), intent(in) :: space
         type(radial_ground_state), intent(in) :: states(:)
         complex(rp), intent(in) :: weighted(:)
         complex(rp), intent(out) :: physical(:)
         type(response_super_index) :: local_item
         integer :: flat, site, ir

         if (size(physical)/=space%ndim .or. size(weighted)/=space%ndim .or. size(states)/=space%nsite) then
            error stop 'DRESP-12: weighted seam vector shape mismatch'
         end if
         physical=cmplx(0.0_rp,0.0_rp,rp)
         do flat=1,space%ndim
            call response_unflatten_superindex(flat,space%nsite,space%response_lmax,space%npoint,space%nchannel,local_item)
            site=local_item%site; ir=local_item%radial_point
            if (local_item%response_l==0 .and. local_item%response_m==0 .and. ir>1) then
               physical(flat)=weighted(flat)/(sqrt(dresp09w_four_pi)*states(site)%r(ir)**2)
            end if
         end do
      end subroutine dresp09y_weighted_to_physical

      real(rp) function dresp09y_weighted_norm(space,states,weighted) result(value)
         type(response_space_layout), intent(in) :: space
         type(radial_ground_state), intent(in) :: states(:)
         complex(rp), intent(in) :: weighted(:)
         type(response_super_index) :: local_item
         integer :: flat, site, ir
         real(rp) :: radius, jacobian

         if (size(weighted)/=space%ndim .or. size(states)/=space%nsite) then
            error stop 'DRESP-12: weighted seam norm shape mismatch'
         end if
         value=0.0_rp
         do flat=1,space%ndim
            call response_unflatten_superindex(flat,space%nsite,space%response_lmax,space%npoint,space%nchannel,local_item)
            site=local_item%site; ir=local_item%radial_point
            if (local_item%response_l==0 .and. local_item%response_m==0 .and. ir>1) then
               radius=states(site)%r(ir)
               jacobian=(2.0_rp*(mod(ir+1,2)+1)/3.0_rp)*states(site)%a*(radius+states(site)%b)
               value=value+jacobian*abs(weighted(flat))**2/(dresp09w_four_pi*radius**2)
            end if
         end do
         value=sqrt(value)
      end function dresp09y_weighted_norm

      real(rp) function dresp09y_weighted_relative(space,states,value,reference) result(result)
         type(response_space_layout), intent(in) :: space
         type(radial_ground_state), intent(in) :: states(:)
         complex(rp), intent(in) :: value(:), reference(:)

         result=dresp09y_weighted_norm(space,states,value-reference)/ &
            max(dresp09y_weighted_norm(space,states,reference),tiny(1.0_rp))
      end function dresp09y_weighted_relative

      subroutine coefficient_trace_oracle(space,radial,states,branches,trace_x,trace_y)
         type(response_space_layout), intent(in) :: space
         type(lmto_radial_basis), intent(in) :: radial(:)
         type(radial_ground_state), intent(in) :: states(:)
         complex(rp), intent(in) :: branches(:, :, :)
         complex(rp), intent(out) :: trace_x(:), trace_y(:)
         complex(rp), allocatable :: density_x(:, :), density_y(:, :)
         real(rp), allocatable :: weighted_x(:, :), weighted_y(:, :)
         integer :: site, ir, iorb, branch, l, m, norb, row_up, col_down, row_down, col_up
         real(rp) :: gaunt, pair_plus, pair_minus
         complex(rp) :: plus_term, minus_term

         norb=(radial(1)%lmax+1)**2
         allocate(density_x(size(radial),space%npoint),density_y(size(radial),space%npoint), &
            weighted_x(size(radial),space%npoint),weighted_y(size(radial),space%npoint))
         density_x=cmplx(0.0_rp,0.0_rp,rp); density_y=cmplx(0.0_rp,0.0_rp,rp)
         do site=1,size(radial)
            do ir=2,space%npoint
               do iorb=1,norb
                  l=lmto_orbital_l(iorb); m=iorb-l*l-l-1
                  gaunt=response_gaunt(l,m,l,m,0,0)
                  row_up=(site-1)*2*norb+iorb; col_down=row_up+norb
                  row_down=col_down; col_up=row_up
                  do branch=1,lmto_product_nbranch
                     call lmto_product_second_order_radial_branch(radial(site),ir,l,l,1,2,branch,pair_plus)
                     call lmto_product_second_order_radial_branch(radial(site),ir,l,l,2,1,branch,pair_minus)
                     plus_term=branches(row_up,col_down,branch)*pair_plus
                     minus_term=branches(row_down,col_up,branch)*pair_minus
                     density_x(site,ir)=density_x(site,ir)+gaunt*(plus_term+minus_term)
                     density_y(site,ir)=density_y(site,ir)+gaunt*cmplx(0.0_rp,1.0_rp,rp)*(plus_term-minus_term)
                  end do
               end do
            end do
         end do
         call dresp09y_convert_response(density_x,states,weighted_x)
         call dresp09y_convert_response(density_y,states,weighted_y)
         call build_l0_real_raw(space,weighted_x,trace_x)
         call build_l0_real_raw(space,weighted_y,trace_y)
         deallocate(density_x,density_y,weighted_x,weighted_y)
      end subroutine coefficient_trace_oracle

      subroutine radial_branch_swap_audit(radial,residuals)
         type(lmto_radial_basis), intent(in) :: radial(:)
         real(rp), intent(out) :: residuals(:)
         integer :: site, ir, l, branch, swapped
         real(rp) :: forward, reverse, residual

         residuals=0.0_rp
         do site=1,size(radial)
            do ir=2,radial(site)%npoint
               do l=0,radial(site)%lmax
                  do branch=1,lmto_product_nbranch
                     swapped=branch_swap(branch)
                     call lmto_product_second_order_radial_branch(radial(site),ir,l,l,1,2,branch,forward)
                     call lmto_product_second_order_radial_branch(radial(site),ir,l,l,2,1,swapped,reverse)
                     residual=abs(forward-reverse)/max(abs(reverse),tiny(1.0_rp))
                     residuals(branch)=max(residuals(branch),residual)
                  end do
               end do
            end do
         end do
      end subroutine radial_branch_swap_audit

      real(rp) function relative_complex_array(left,right) result(value)
         complex(rp), intent(in) :: left(:, :), right(:, :)
         value=sqrt(sum(abs(left-right)**2))/max(sqrt(sum(abs(right)**2)),tiny(1.0_rp))
      end function relative_complex_array

      integer function branch_swap(branch) result(swapped)
         integer, intent(in) :: branch
         integer :: p, q, candidate, candidate_p, candidate_q

         call lmto_product_branch_powers(branch,p,q)
         swapped=branch
         do candidate=1,lmto_product_nbranch
            call lmto_product_branch_powers(candidate,candidate_p,candidate_q)
            if(candidate_p==q .and. candidate_q==p) then
               swapped=candidate
               return
            end if
         end do
      end function branch_swap

      subroutine reciprocal_kpoint_local(rec,ik,array4d,result)
         type(reciprocal), intent(in) :: rec
         integer, intent(in) :: ik
         complex(rp), intent(in) :: array4d(:,:,:,:)
         complex(rp), intent(out) :: result(:,:)
         real(rp) :: point(3)
         if (allocated(rec%k_workset%points)) then; point=rec%k_workset%points(:,ik); else; point=rec%k_points(:,ik); end if
         call rec%fourier_transform_array(array4d,point,result)
      end subroutine reciprocal_kpoint_local

      subroutine build_matrix_decomposition(ik,ham,lat,rec,dee,dob,den,hf,ov,en,dhf,dov,denu,dh2,te,th,tl,tm,tr,ce,ch,co,dc,pr)
         integer, intent(in) :: ik
         type(hamiltonian), intent(in) :: ham
         type(lattice), intent(in) :: lat
         type(reciprocal), intent(in) :: rec
         complex(rp), intent(in) :: dee(:,:,:,:),dob(:,:,:),den(:,:,:),dc(:,:)
         complex(rp), intent(out) :: hf(:,:),ov(:,:),en(:,:),dhf(:,:),dov(:,:),denu(:,:),dh2(:,:),te(:,:),th(:,:),tl(:,:),tm(:,:),tr(:,:),ce(:,:),ch(:,:),co(:,:)
         real(rp), intent(inout) :: pr
         complex(rp), allocatable :: h2a(:,:), hoha(:,:)
         allocate(h2a(size(hf,1),size(hf,2)),hoha(size(hf,1),size(hf,2)))
         call reciprocal_kpoint_local(rec,ik,ham%ee,hf); call reciprocal_kpoint_local(rec,ik,dee,dhf)
         call dresp08_block_diagonal(ham%obarm,lat%ib(1:nsite),ov)
         call dresp08_block_diagonal(ham%enim,lat%ib(1:nsite),en)
         call dresp08_block_diagonal(dob,lat%ib(1:nsite),dov)
         call dresp08_block_diagonal(den,lat%ib(1:nsite),denu)
         call dresp08_build_h2(hf,ov,en,h2a,hoha)
         call dresp08_build_product_tangent(denu,dhf,dov,hf,ov,dh2,te,th,tl,tm,tr)
         pr=max(pr,dresp08_relative_residual(dh2,dc))
         ce=te; ch=th-dc; co=-tl-tm-tr
         deallocate(h2a,hoha)
      end subroutine build_matrix_decomposition

      subroutine fixed_response_for_field(space,radial,prod,rec,hmat,field,ik,raw,compact)
         type(response_space_layout), intent(in) :: space
         type(lmto_radial_basis), intent(in) :: radial(:)
         type(lmto_product_response_basis), intent(in) :: prod
         type(reciprocal), intent(in) :: rec
         complex(rp), intent(in) :: hmat(:,:),field(:,:)
         integer, intent(in) :: ik
         complex(rp), intent(inout) :: raw(:),compact(:)
         complex(rp), allocatable :: dr(:,:),point(:),tmp(:)
         allocate(dr(size(field,1),size(field,2)),point(space%ndim),tmp(prod%product_dimension))
         call lr_static_frechet_density(rec%eigenvalues(:,ik),rec%eigenvectors(:,:,ik),rec%fermi_level,rec%temperature,field,dr)
         call fixed_pauli_measurement(space,radial,hmat,dr,point); call dresp09_compact_field_from_raw(space,prod,point,tmp)
         raw=raw+rec%k_weights(ik)*point; compact=compact+rec%k_weights(ik)*tmp
         deallocate(dr,point,tmp)
      end subroutine fixed_response_for_field

      subroutine accumulate_response(point,weight,accum)
         complex(rp), intent(in) :: point(:)
         real(rp), intent(in) :: weight
         complex(rp), intent(inout) :: accum(:)
         accum=accum+weight*point
      end subroutine accumulate_response

      subroutine field_identity_metrics(b,c,v,w,maxr,rms,maxel,nbv,ncv,nvv)
         complex(rp), intent(in) :: b(:,:),c(:,:),v(:,:)
         real(rp), intent(in) :: w
         real(rp), intent(inout) :: maxr,rms,maxel,nbv,ncv,nvv
         real(rp) :: res
         res=sqrt(sum(abs(b+c-v)**2))/max(sqrt(sum(abs(v)**2)),tiny(1.0_rp))
         maxr=max(maxr,res); rms=rms+w*res*res; maxel=max(maxel,maxval(abs(b+c-v)))
         nbv=nbv+w*sum(abs(b)**2); ncv=ncv+w*sum(abs(c)**2); nvv=nvv+w*sum(abs(v)**2)
      end subroutine field_identity_metrics

      subroutine branch_metrics(hmat,comps,w,v,c,bs,bms)
         complex(rp), intent(in) :: hmat(:,:),comps(:,:,:,:),v(:,:),c(:,:)
         real(rp), intent(in) :: w
         real(rp), intent(inout) :: bs(:),bms(:)
         complex(rp), allocatable :: term(:,:)
         integer :: ib
         allocate(term(size(hmat,1),size(hmat,2)))
         do ib=1,lmto_product_nbranch
            call lmto_product_apply_branch_action(hmat,comps(:,:,ib,lr_full_spatial_piece_total),ib,term,.true.)
            bs(ib)=bs(ib)+w*sum(abs(term)**2); bms(ib)=bms(ib)+w*sum(abs(term-c)**2)
         end do
         deallocate(term)
      end subroutine branch_metrics

      subroutine add_orbital_blocks(matrix,w,blocks,left,right,interf)
         complex(rp), intent(in) :: matrix(:,:)
         real(rp), intent(in) :: w
         real(rp), intent(inout) :: blocks(:,:),left,right,interf
         integer :: a,b,ii,jj,li,lj,norb
         norb=(radial_bases(1)%lmax+1)**2
         do a=1,size(matrix,1)
            li=lmto_orbital_l(mod(a-1,norb)+1)+1
            do b=1,size(matrix,2)
               lj=lmto_orbital_l(mod(b-1,norb)+1)+1
               blocks(1,li)=blocks(1,li)+w*abs(matrix(a,b))**2
               if (li/=lj) blocks(2,min(li,lj))=blocks(2,min(li,lj))+w*abs(matrix(a,b))**2
               if (a<=norb .and. b>norb .or. a>norb .and. b<=norb) then
                  blocks(3,min(li,3))=blocks(3,min(li,3))+w*abs(matrix(a,b))**2
               end if
            end do
         end do
         left=left+w*sum(abs(matrix(1:min(norb,size(matrix,1)),:))**2)
         right=right+w*sum(abs(matrix(max(1,size(matrix,1)-norb+1):size(matrix,1),:))**2)
         interf=interf+w*(sum(abs(matrix)**2)-left-right)
      end subroutine add_orbital_blocks

      subroutine band_actions(matrix,evec,eval,ef,allv,occv,nearv,maxdist,residual,dist)
         complex(rp), intent(in) :: matrix(:,:),evec(:,:)
         real(rp), intent(in) :: eval(:),ef
         real(rp), intent(inout) :: allv,occv,nearv,maxdist,residual,dist
         integer :: ib
         real(rp) :: value
         allv=0.0_rp; occv=0.0_rp; nearv=0.0_rp; maxdist=0.0_rp
         do ib=1,size(eval)
            value=sqrt(sum(abs(matmul(matrix,evec(:,ib)))**2)); allv=max(allv,value)
            if(eval(ib)<=ef) occv=max(occv,value)
            if(abs(eval(ib)-ef)<=0.1_rp) nearv=max(nearv,value)
            if(value>=maxdist) then; maxdist=value; dist=abs(eval(ib)-ef); end if
         end do
         residual=allv
      end subroutine band_actions

      subroutine build_l0_raw(space,values,raw)
         type(response_space_layout), intent(in) :: space
         complex(rp), intent(in) :: values(:,:)
         complex(rp), intent(out) :: raw(:)
         integer :: flat
         raw=cmplx(0.0_rp,0.0_rp,rp)
         do flat=1,size(raw)
            call response_unflatten_superindex(flat,space%nsite,space%response_lmax,space%npoint,space%nchannel,item)
            if(item%response_l==0 .and. item%response_m==0 .and. item%radial_point>1) raw(flat)=values(item%site,item%radial_point)
         end do
      end subroutine build_l0_raw

      subroutine project_l0_vector(space,vector,projected)
         type(response_space_layout), intent(in) :: space
         complex(rp), intent(in) :: vector(:)
         complex(rp), intent(out) :: projected(:)
         type(response_super_index) :: local_item
         integer :: local_flat

         if (size(vector) /= space%ndim .or. size(projected) /= space%ndim) then
            error stop 'DRESP-12 project_l0_vector: vector shape mismatch'
         end if
         projected=cmplx(0.0_rp,0.0_rp,rp)
         do local_flat=1,space%ndim
            call response_unflatten_superindex(local_flat,space%nsite,space%response_lmax,space%npoint, &
               space%nchannel,local_item)
            if (local_item%response_l==0 .and. local_item%response_m==0) projected(local_flat)=vector(local_flat)
         end do
      end subroutine project_l0_vector

      subroutine project_l0_compact(product,vector,projected)
         type(lmto_product_response_basis), intent(in) :: product
         complex(rp), intent(in) :: vector(:)
         complex(rp), intent(out) :: projected(:)
         integer :: flat, site_index, response_l, response_m, mode_index

         if (size(vector) /= product%product_dimension .or. size(projected) /= product%product_dimension) then
            error stop 'DRESP-12 project_l0_compact: vector shape mismatch'
         end if
         projected=cmplx(0.0_rp,0.0_rp,rp)
         do flat=1,product%product_dimension
            call product%unflatten_index(flat,site_index,response_l,response_m,mode_index)
            if (response_l==0 .and. response_m==0) projected(flat)=vector(flat)
         end do
      end subroutine project_l0_compact

      subroutine l0_mode_geometry(space,magnetization,mode,cosine,parallel,orthogonal)
         type(response_space_layout), intent(in) :: space
         complex(rp), intent(in) :: magnetization(:),mode(:)
         real(rp), intent(out) :: cosine,parallel,orthogonal
         complex(rp) :: overlap, coefficient
         complex(rp), allocatable :: residual(:)
         real(rp) :: magnetization_norm, mode_norm

         overlap=response_vector_inner_product(space,magnetization,mode)
         magnetization_norm=response_space_norm(space,magnetization)
         mode_norm=response_space_norm(space,mode)
         coefficient=overlap/max(real(response_vector_inner_product(space,magnetization,magnetization),rp),tiny(1.0_rp))
         allocate(residual(size(mode)))
         residual=mode-magnetization*coefficient
         cosine=abs(overlap)/max(magnetization_norm*mode_norm,tiny(1.0_rp))
         parallel=real(coefficient,rp)
         orthogonal=response_space_norm(space,residual)/max(mode_norm,tiny(1.0_rp))
         deallocate(residual)
      end subroutine l0_mode_geometry

      real(rp) function response_space_norm(space,value) result(value_norm)
         type(response_space_layout), intent(in) :: space
         complex(rp), intent(in) :: value(:)
         value_norm=sqrt(sum(space%metric_weights*abs(value)**2))
      end function response_space_norm

      ! Partition by response angular momentum while retaining the production
      ! response-space metric: radial quadrature weight times |v_LM(r)|^2.
      ! The check tests implementation/partition consistency under the
      ! existing diagonal metric; it is not an independent proof of physical
      ! angular orthogonality.
      subroutine norm_by_l(space,vector,norms)
         type(response_space_layout), intent(in) :: space
         complex(rp), intent(in) :: vector(:)
         real(rp), intent(out) :: norms(0:4)
         type(response_super_index) :: local_item
         integer :: local_flat

         if (size(vector) /= space%ndim) error stop 'DRESP-12 norm_by_l: vector shape mismatch'
         norms=0.0_rp
         do local_flat=1,space%ndim
            call response_unflatten_superindex(local_flat,space%nsite,space%response_lmax,space%npoint, &
               space%nchannel,local_item)
            if (local_item%response_l <= 4) norms(local_item%response_l)=norms(local_item%response_l)+ &
               space%metric_weights(local_flat)*abs(vector(local_flat))**2
         end do
         norms=sqrt(max(0.0_rp,norms))
      end subroutine norm_by_l

      real(rp) function weighted_integral(space,value) result(result)
         type(response_space_layout), intent(in) :: space
         complex(rp), intent(in) :: value(:)
         result=sum(space%metric_weights*real(value,rp))
      end function weighted_integral

      real(rp) function norm_observable(value) result(result)
         complex(rp), intent(in) :: value(:)
         result=response_space_norm(response_space,value)
      end function norm_observable

      subroutine vector_geometry(space,r,conn,endpoint_value,obs,g,oc,oe,oo,cc,ce,co,ac,ae,ao,pf,qf)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: r(:),conn(:),endpoint_value(:),obs(:),g(:)
      real(rp), intent(out) :: oc,oe,oo,cc,ce,co,ac,ae,ao,pf(5),qf(5)
         complex(rp) :: rc,re,ro,rg
         complex(rp), allocatable :: vals(:,:),proj(:)
         real(rp) :: nr,nc,ne,no,nv
         nr=response_space_norm(space,r); nc=response_space_norm(space,conn); ne=response_space_norm(space,endpoint_value)
         no=response_space_norm(space,obs)
         rc=response_vector_inner_product(space,r,conn); ro=response_vector_inner_product(space,r,obs)
         re=response_vector_inner_product(space,r,endpoint_value)
         rg=response_vector_inner_product(space,g,g)
         oc=real(rc/max(nr*nc,tiny(1.0_rp)),rp); oe=real(re/max(nr*ne,tiny(1.0_rp)),rp)
         oo=real(ro/max(nr*no,tiny(1.0_rp)),rp); cc=oc; ce=oe; co=oo
         ac=acos(max(-1.0_rp,min(1.0_rp,cc))); ae=acos(max(-1.0_rp,min(1.0_rp,ce))); ao=acos(max(-1.0_rp,min(1.0_rp,co)))
         allocate(vals(size(g),5),proj(size(g))); vals(:,1)=r; vals(:,2)=conn; vals(:,3)=endpoint_value; vals(:,4)=obs
         vals(:,5)=r+conn+endpoint_value+obs
         do mode=1,5
            proj=g*response_vector_inner_product(space,g,vals(:,mode))/max(real(rg,rp),tiny(1.0_rp)); nv=response_space_norm(space,vals(:,mode))
            pf(mode)=response_space_norm(space,proj)/max(nv,tiny(1.0_rp)); qf(mode)=response_space_norm(space,vals(:,mode)-proj)/max(nv,tiny(1.0_rp))
         end do
         deallocate(vals,proj)
      end subroutine vector_geometry

      real(rp) function relative_vector(value,reference) result(result)
         complex(rp), intent(in) :: value(:),reference(:)
         result=sqrt(sum(abs(value-reference)**2))/max(sqrt(sum(abs(reference)**2)),tiny(1.0_rp))
      end function relative_vector

      real(rp) function hermitian_matrix_residual(left,right) result(result)
         complex(rp), intent(in) :: left(:, :),right(:, :)
         result=sqrt(sum(abs(conjg(transpose(left))-right)**2))/max(sqrt(sum(abs(left)**2))+sqrt(sum(abs(right)**2)),tiny(1.0_rp))
      end function hermitian_matrix_residual

      real(rp) function median_value(values) result(result)
         real(rp), intent(in) :: values(:)
         real(rp), allocatable :: sorted(:)
         integer :: i,j,n; real(rp) :: temp
         allocate(sorted(size(values))); sorted=values; n=size(sorted)
         do i=1,n-1; do j=i+1,n; if(sorted(j)<sorted(i)) then; temp=sorted(i);sorted(i)=sorted(j);sorted(j)=temp;end if;end do;end do
         if(mod(n,2)==1) then; result=sorted((n+1)/2); else; result=0.5_rp*(sorted(n/2)+sorted(n/2+1)); end if
         deallocate(sorted)
      end function median_value

      subroutine write_k_table(file_name,ratios,band,dist)
         character(len=*), intent(in) :: file_name
         real(rp), intent(in) :: ratios(:),band(:),dist(:)
         integer :: u,ios,ik
         if(rank/=0)return
         open(newunit=u,file=file_name,status='replace',action='write',iostat=ios); if(ios/=0)return
         write(u,'(a)') 'k_index,R_H,band_action_max,distance_to_EF'
         do ik=1,size(ratios); write(u,'(i0,3(a,es24.16))') ik,',',ratios(ik),',',band(ik),',',dist(ik); end do
         close(u)
      end subroutine write_k_table

   end subroutine run_dresp12_covariance

end module lr_dresp12_covariance_bridge_mod
