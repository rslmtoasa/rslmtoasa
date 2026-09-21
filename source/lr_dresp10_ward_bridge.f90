!------------------------------------------------------------------------------
! DRESP-10 RAW full-spatial transverse ALSDA Ward bridge.
!
! This is an accepted-state diagnostic bridge.  It owns no correction,
! Goldstone repair, empirical rescaling, Dyson solve, or spectrum.  The raw
! material gates are evaluated from the frozen reciprocal eigensystem and the
! independent full-spatial source/Frechet services.  The mixed K/contact
! denominator remains explicitly fail-closed until its arbitrary-(L,M)
! operator action is assembled.
!------------------------------------------------------------------------------
module lr_dresp10_ward_bridge_mod

   use precision_mod, only: rp
   use mpi_mod, only: rank
   use basis_mod, only: nb
   use reciprocal_mod, only: reciprocal
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use radial_ground_state_mod, only: radial_ground_state
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use response_angular_basis_mod, only: response_angular_pi
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex, &
      response_unflatten_superindex
   use lr_response_space_mod, only: response_space_layout
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, lmto_product_channel_plus, &
      lmto_product_nbranch
   use lr_full_spatial_alsda_mod, only: lr_full_spatial_source, lr_full_spatial_source_pieces, &
      lr_full_spatial_source_components, lr_full_spatial_contract_components, lr_full_spatial_piece_total, lr_full_spatial_npiece
   use lr_static_frechet_product_response_mod, only: lr_static_product_matrix, lr_static_product_action, &
      lr_static_frechet_density
   use lr_dresp09_field_insertion_mod, only: dresp09_product_source_operator, dresp09_compact_field_from_raw
   use lr_exact_ks_ward_mod, only: exact_ks_spin_rotation_generator
   use lr_lmto_density_moment_tangent_mod, only: density_from_eigensystem, moment_matrices_from_eigensystem, &
      endpoint_matrices_from_moments_second_order, endpoint_tangent_branches_second_order, relative_matrix_residual, &
      commutator_tangent
   use lr_dresp09s_scalar_relativistic_mod, only: sr_l0_density_from_second_order_endpoint_branches, sr_l0_source
   use lr_dresp09y_bridge_mod, only: dresp09y_augmentation_density
   use lr_dresp09zs_compact_span_mod, only: dresp09zs_source_span_residuals
   use lr_radial_observable_provenance_mod, only: channel_moments_from_matrix, pauli_density_from_moments, &
      sr_spin_density_from_moments, weighted_from_density, dresp09w_four_pi, radial_relative_metrics, volume_l2_norm, volume_integral
   use pauli_ground_state_projection_mod, only: compute_accepted_pauli_magnetization
   use lr_sr_augmentation_tangent_mod, only: sr_aug_sigma_x, sr_aug_lower_spin_factor
   implicit none
   private

   character(len=*), parameter, public :: dresp10_verdict_pass_a = 'PASS-A'
   character(len=*), parameter, public :: dresp10_verdict_pass_b = 'PASS-B'
   character(len=*), parameter, public :: dresp10_verdict_blocked = 'BLOCKED'

   public :: run_dresp10_raw_full_spatial_alsda_ward

contains

   subroutine run_dresp10_raw_full_spatial_alsda_ward(output_file, response_space, radial_bases, ground_states, &
                                                      reciprocal_obj, lattice_obj, hamiltonian_obj)
      character(len=*), intent(in) :: output_file
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(reciprocal), intent(inout) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj

      type(lmto_product_response_basis) :: product
      type(response_super_index) :: item
      complex(rp), allocatable :: h(:, :), rho(:, :), moments(:, :, :), moment_sum(:, :, :)
      complex(rp), allocatable :: endpoints(:, :, :), endpoint_sum(:, :, :), delta_endpoint(:, :, :), delta_endpoint_sum(:, :, :)
      complex(rp), allocatable :: generator(:, :), delta_native(:, :), delta_alsda(:, :), delta_rho(:, :), source(:, :)
      complex(rp), allocatable :: source_pieces(:, :, :), source_pieces_minus(:, :, :), source_components(:, :, :, :), &
         source_components_minus(:, :, :, :)
      complex(rp), allocatable :: response_plus(:, :), response_minus(:, :), native_response_c(:, :), pauli_native_response_c(:, :)
      complex(rp), allocatable :: alsda_response_c(:, :), augmentation_upper_c(:, :), augmentation_small_c(:, :), &
         augmentation_angular_c(:, :), augmentation_rank0_c(:, :), augmentation_rank2_c(:, :), &
         augmentation_response_c(:, :), pauli_augmentation_c(:, :), sigma_x(:, :)
      complex(rp), allocatable :: vertices(:, :, :, :), chi_static(:, :), trial_field(:), action_matrix(:)
      complex(rp), allocatable :: action_direct(:), action_delta(:)
      real(rp), allocatable :: channels(:, :, :, :), p1(:, :), p3(:, :), target_p3(:, :), target_sr(:, :), target_sr_raw(:, :)
      real(rp), allocatable :: bxc(:, :), kxc(:, :), native_response(:, :), alsda_response(:, :)
      complex(rp), allocatable :: source_field(:)
      real(rp), allocatable :: augmentation_upper(:, :), augmentation_small(:, :), augmentation_angular(:, :), &
         augmentation_rank0(:, :), augmentation_rank2(:, :), &
         augmentation_response(:, :), complete_response(:, :), pauli_native_response(:, :), pauli_augmentation(:, :), &
         pauli_complete_response(:, :)
      real(rp), allocatable :: p3_l(:, :, :), p1_raw(:, :)
      real(rp) :: wsum, min_source_norm, field_residual, field_max_element, p3_identity_residual
      real(rp) :: native_complete_residual, alsda_complete_residual, static_action_residual, legacy_l0_field_residual
      real(rp) :: static_pseudorandom_residual, static_physical_residual
      real(rp) :: native_fixed_residual, native_pauli_residual, pauli_complete_residual, fixed_augmentation_residual
      real(rp) :: contact_rank_recomposition_residual
      real(rp) :: fixed_integral, complete_integral, target_integral, fixed_volume_norm, complete_volume_norm, target_volume_norm
      real(rp) :: full_vs_legacy_field_residual
      real(rp) :: p1_p3_relative, min_abs_p3, max_abs_p3, min_abs_kxc, max_abs_kxc
      real(rp) :: source_upper_norm, source_small_norm, source_rank0_norm, source_rank2_norm
      real(rp) :: source_total_norm, native_norm, target_norm
      real(rp) :: frechet_native_residual
      real(rp) :: field_weighted_rms, field_occupied_action, field_near_ef_action, field_current_relative
      integer :: field_noccupied, field_nnear
      real(rp) :: source_span_by_l(7,0:4), source_span_physical_by_l(0:4), source_span_delta_by_l(0:4)
      real(rp) :: source_span_physical_global, source_span_delta_global
      integer :: nsite, norb_site, nmat, nk, ndim, product_dim, ik, ik_global, flat, ir, site, l, m, unit, ios
      integer :: probe_l, probe_flat, branch, piece, i
      logical :: source_complete, source_probe_complete, field_closed, p3_closed, frechet_closed, mixed_operator_closed
      logical :: augmentation_closed
      character(len=32) :: verdict
      character(len=96) :: classification

      nsite = size(ground_states)
      if (nsite < 1 .or. size(radial_bases) /= nsite .or. response_space%nsite /= nsite) then
         error stop 'DRESP-10: response/radial site dimensions are inconsistent'
      end if
      if (.not. hamiltonian_obj%hoh .or. hamiltonian_obj%ccor_2c .or. hamiltonian_obj%orb_pol .or. &
          hamiltonian_obj%hubbard_u_general_check .or. hamiltonian_obj%hubbard_v_check .or. hamiltonian_obj%local_axis) then
         error stop 'DRESP-10: excluded Hamiltonian correction is enabled'
      end if
      if (.not. allocated(reciprocal_obj%hk_bulk) .or. .not. allocated(reciprocal_obj%eigenvalues) .or. &
          .not. allocated(reciprocal_obj%eigenvectors) .or. .not. allocated(reciprocal_obj%k_weights)) then
         error stop 'DRESP-10: accepted reciprocal eigensystem is incomplete'
      end if
      if (response_space%response_lmax /= 4 .or. response_space%nchannel /= 1) then
         error stop 'DRESP-10: complete L=0..4 single-channel response space is required'
      end if

      norb_site = (radial_bases(1)%lmax + 1)**2
      nmat = size(reciprocal_obj%hk_bulk, 1)
      nk = size(reciprocal_obj%hk_bulk, 3)
      ndim = response_space%ndim
      if (nmat /= nb*nsite .or. nmat /= 2*norb_site*nsite .or. size(reciprocal_obj%k_weights) /= nk) then
         error stop 'DRESP-10: accepted site-major Fe eigensystem dimensions are incomplete'
      end if
      if (lattice_obj%nrec /= nsite .or. nk /= 64) error stop 'DRESP-10: certified bcc-Fe material gate failed'
      wsum = sum(reciprocal_obj%k_weights)
      if (wsum <= 0.0_rp) error stop 'DRESP-10: nonpositive reciprocal weight sum'

      call product%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      product_dim = product%product_dimension
      if (product_dim /= 348 .or. product%maximum_gf_energy_moment /= 4) then
         error stop 'DRESP-10: accepted Fe six-branch product basis contract failed'
      end if
      call dresp09zs_source_span_residuals(response_space, radial_bases, product, source_span_physical_by_l, &
         source_span_delta_by_l, source_span_by_l, source_span_physical_global, source_span_delta_global)
      source_complete = maxval(source_span_by_l) < 1.0e-10_rp

      allocate(moment_sum(nmat,nmat,3), endpoint_sum(nmat,nmat,6), h(nmat,nmat), rho(nmat,nmat), &
         moments(nmat,nmat,3), endpoints(nmat,nmat,6), generator(nmat,nmat), delta_native(nmat,nmat), &
         delta_alsda(nmat,nmat), delta_rho(nmat,nmat), source(nmat,nmat), source_pieces(nmat,nmat,lr_full_spatial_npiece), &
         delta_endpoint(nmat,nmat,6), delta_endpoint_sum(nmat,nmat,6))
      moment_sum = cmplx(0.0_rp,0.0_rp,rp)
      endpoint_sum = cmplx(0.0_rp,0.0_rp,rp)
      delta_endpoint_sum = cmplx(0.0_rp,0.0_rp,rp)
      call exact_ks_spin_rotation_generator(norb_site, nsite, generator)
      do ik = 1, nk
         ik_global = ik
         if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
         h = reciprocal_obj%hk_bulk(:,:,ik)
         call density_from_eigensystem(reciprocal_obj%eigenvalues(:,ik), reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level, reciprocal_obj%temperature, rho)
         call moment_matrices_from_eigensystem(reciprocal_obj%eigenvalues(:,ik), reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level, reciprocal_obj%temperature, moments)
         call endpoint_matrices_from_moments_second_order(moments, endpoints)
         call commutator_tangent(generator,h,delta_native)
         call commutator_tangent(generator,rho,delta_rho)
         call endpoint_tangent_branches_second_order(h,rho,delta_native,delta_rho,delta_endpoint)
         moment_sum = moment_sum + reciprocal_obj%k_weights(ik_global)*moments
         endpoint_sum = endpoint_sum + reciprocal_obj%k_weights(ik_global)*endpoints
         delta_endpoint_sum = delta_endpoint_sum + reciprocal_obj%k_weights(ik_global)*delta_endpoint
      end do
      moment_sum = moment_sum/wsum
      endpoint_sum = endpoint_sum/wsum
      delta_endpoint_sum = delta_endpoint_sum/wsum

      allocate(channels(3,nsite,radial_bases(1)%lmax+1,2), p1(nsite,response_space%npoint), &
         p3(nsite,response_space%npoint), target_p3(nsite,response_space%npoint), target_sr(nsite,response_space%npoint), &
         target_sr_raw(nsite,response_space%npoint), &
         p1_raw(nsite,response_space%npoint), &
         p3_l(nsite,response_space%npoint,radial_bases(1)%lmax+1), bxc(nsite,response_space%npoint), &
         kxc(nsite,response_space%npoint), source_field(ndim), native_response(nsite,response_space%npoint), &
         alsda_response(nsite,response_space%npoint), augmentation_upper(nsite,response_space%npoint), &
         augmentation_small(nsite,response_space%npoint), augmentation_angular(nsite,response_space%npoint), &
         augmentation_rank0(nsite,response_space%npoint), augmentation_rank2(nsite,response_space%npoint), &
         augmentation_response(nsite,response_space%npoint), &
         complete_response(nsite,response_space%npoint), pauli_native_response(nsite,response_space%npoint), &
         pauli_augmentation(nsite,response_space%npoint), pauli_complete_response(nsite,response_space%npoint), &
         response_plus(nsite,response_space%npoint), &
         response_minus(nsite,response_space%npoint), native_response_c(nsite,response_space%npoint), &
         pauli_native_response_c(nsite,response_space%npoint), &
         alsda_response_c(nsite,response_space%npoint), augmentation_upper_c(nsite,response_space%npoint), &
         augmentation_small_c(nsite,response_space%npoint), augmentation_angular_c(nsite,response_space%npoint), &
         augmentation_rank0_c(nsite,response_space%npoint), augmentation_rank2_c(nsite,response_space%npoint), &
         augmentation_response_c(nsite,response_space%npoint), pauli_augmentation_c(nsite,response_space%npoint), sigma_x(2,2))
      call channel_moments_from_matrix(moment_sum, nsite, channels)
      call pauli_density_from_moments(radial_bases, channels, p3, p3_l)
      call sr_spin_density_from_moments(radial_bases, channels, target_sr_raw)
      call compute_accepted_pauli_magnetization(reciprocal_obj, lattice_obj%symbolic_atoms, lattice_obj%nbulk, p1_raw)
      call weighted_from_density(ground_states(1)%r, p1_raw, p1)
      ! The certified DRESP-09Y observable regression compares the converted response
      ! profile against the physical (unweighted) P3/SR target.  Keep these
      ! targets physical here as well; only p1 is a separately weighted field.
      target_p3 = p3
      target_sr = target_sr_raw

      do site = 1, nsite
         do ir = 1, response_space%npoint
            bxc(site,ir) = 0.5_rp*(ground_states(site)%vxc_up(ir) - ground_states(site)%vxc_down(ir))
            if (response_space%radial_weights(ir) > 0.0_rp) then
               if (p3(site,ir) == 0.0_rp) error stop 'SR_SOURCE_SPACE_INCOMPLETE'
               kxc(site,ir) = bxc(site,ir)/p3(site,ir)
            else
               kxc(site,ir) = 0.0_rp
            end if
         end do
      end do
      min_abs_p3 = minval(abs(p3), mask=spread(response_space%radial_weights > 0.0_rp,1,nsite))
      max_abs_p3 = maxval(abs(p3))
      min_abs_kxc = minval(abs(kxc), mask=spread(response_space%radial_weights > 0.0_rp,1,nsite))
      max_abs_kxc = maxval(abs(kxc))
      p3_identity_residual = maxval(abs(kxc*p3 - bxc), mask=spread(response_space%radial_weights > 0.0_rp,1,nsite))

      source_field = cmplx(0.0_rp,0.0_rp,rp)
      do flat = 1, ndim
         call response_unflatten_superindex(flat, response_space%nsite, response_space%response_lmax, &
            response_space%npoint, response_space%nchannel, item)
         if (item%response_l == 0 .and. item%response_m == 0) then
            source_field(flat) = cmplx(sqrt(4.0_rp*response_angular_pi)*bxc(item%site,item%radial_point),0.0_rp,rp)
         end if
      end do

      ! Preflight every requested source harmonic using one positive-measure
      ! radial point.  This catches a stable orthogonal component missing from
      ! the six-branch full-spatial map before any Ward number is reported.
      allocate(source_components(nmat,nmat,lmto_product_nbranch,lr_full_spatial_npiece))
      min_source_norm = huge(1.0_rp)
      do probe_l = 0, response_space%response_lmax
         source_field = cmplx(0.0_rp,0.0_rp,rp)
         item = response_super_index(1,probe_l,0,response_space%npoint,1)
         call response_flatten_superindex(item,response_space%nsite,response_space%response_lmax,response_space%npoint,1,probe_flat)
         source_field(probe_flat) = cmplx(1.0_rp,0.0_rp,rp)
         call lr_full_spatial_source_components(response_space,radial_bases,source_field,1,source_components)
         source_total_norm = sqrt(sum(abs(source_components(:,:,1,lr_full_spatial_piece_total))**2))
         do branch = 2, lmto_product_nbranch
            source_total_norm = max(source_total_norm, sqrt(sum(abs(source_components(:,:,branch,lr_full_spatial_piece_total))**2)))
         end do
         min_source_norm = min(min_source_norm,source_total_norm)
      end do
      source_probe_complete = min_source_norm > 1.0e-14_rp
      source_field = cmplx(0.0_rp,0.0_rp,rp)
      do flat = 1, ndim
         call response_unflatten_superindex(flat,response_space%nsite,response_space%response_lmax,response_space%npoint,1,item)
         if (item%response_l == 0 .and. item%response_m == 0) &
            source_field(flat)=cmplx(sqrt(4.0_rp*response_angular_pi)*bxc(item%site,item%radial_point),0.0_rp,rp)
      end do
      allocate(source_components_minus(nmat,nmat,lmto_product_nbranch,lr_full_spatial_npiece), &
         source_pieces_minus(nmat,nmat,lr_full_spatial_npiece))
      call lr_full_spatial_source_components(response_space,radial_bases,source_field,1,source_components)
      call lr_full_spatial_source_components(response_space,radial_bases,source_field,2,source_components_minus)

      native_response_c = cmplx(0.0_rp,0.0_rp,rp)
      pauli_native_response_c = cmplx(0.0_rp,0.0_rp,rp)
      alsda_response_c = cmplx(0.0_rp,0.0_rp,rp)
      field_residual = 0.0_rp; field_max_element = 0.0_rp
      field_weighted_rms = 0.0_rp; field_occupied_action = 0.0_rp; field_near_ef_action = 0.0_rp
      frechet_native_residual = 0.0_rp
      legacy_l0_field_residual = 0.0_rp
      full_vs_legacy_field_residual = 0.0_rp
      source_upper_norm = 0.0_rp; source_small_norm = 0.0_rp; source_rank0_norm = 0.0_rp; source_rank2_norm = 0.0_rp; source_total_norm = 0.0_rp
      do ik = 1, nk
         ik_global = ik
         if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
         h = reciprocal_obj%hk_bulk(:,:,ik)
         call commutator_tangent(generator,h,delta_native)
         call lr_full_spatial_contract_components(source_components,h,source_pieces)
         call lr_full_spatial_contract_components(source_components_minus,h,source_pieces_minus)
         call sr_l0_source(response_space,radial_bases,source_field,1,h,source)
         call sr_l0_source(response_space,radial_bases,source_field,2,h,delta_alsda)
         delta_rho = source + delta_alsda
         legacy_l0_field_residual = max(legacy_l0_field_residual,relative_matrix_residual(delta_rho,delta_native))
         source = source_pieces(:,:,lr_full_spatial_piece_total)
         delta_alsda = source + source_pieces_minus(:,:,lr_full_spatial_piece_total)
         full_vs_legacy_field_residual = max(full_vs_legacy_field_residual, &
            relative_matrix_residual(delta_alsda,delta_rho))
         field_current_relative = relative_matrix_residual(delta_alsda,delta_native)
         field_residual = max(field_residual,field_current_relative)
         field_weighted_rms = field_weighted_rms + reciprocal_obj%k_weights(ik_global)*field_current_relative**2
         field_max_element = max(field_max_element,maxval(abs(delta_alsda-delta_native)))
         source = delta_alsda - delta_native
         call action_metrics_local(source,delta_native,reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%eigenvalues(:,ik),reciprocal_obj%fermi_level,field_occupied_action, &
            field_near_ef_action,field_noccupied,field_nnear)
         source_upper_norm = max(source_upper_norm,sqrt(sum(abs(source_pieces(:,:,1))**2)))
         source_small_norm = max(source_small_norm,sqrt(sum(abs(source_pieces(:,:,2))**2)))
         source_rank0_norm = max(source_rank0_norm,sqrt(sum(abs(source_pieces(:,:,3))**2)))
         source_rank2_norm = max(source_rank2_norm,sqrt(sum(abs(source_pieces(:,:,4))**2)))
         source_total_norm = max(source_total_norm,sqrt(sum(abs(source_pieces(:,:,5))**2)))

         call density_from_eigensystem(reciprocal_obj%eigenvalues(:,ik),reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level,reciprocal_obj%temperature,rho)
         call lr_static_frechet_density(reciprocal_obj%eigenvalues(:,ik),reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level,reciprocal_obj%temperature,delta_native,delta_rho)
         call commutator_tangent(generator,rho,source)
         frechet_native_residual = max(frechet_native_residual,relative_matrix_residual(delta_rho,source))
         call endpoint_tangent_branches_second_order(h,rho,delta_native,delta_rho,delta_endpoint)
         call sr_l0_density_from_second_order_endpoint_branches(response_space,radial_bases,delta_endpoint,1,response_plus)
         call sr_l0_density_from_second_order_endpoint_branches(response_space,radial_bases,delta_endpoint,2,response_minus)
         native_response_c = native_response_c + reciprocal_obj%k_weights(ik_global)*(response_plus+response_minus)
         call sr_l0_density_from_second_order_endpoint_branches(response_space,radial_bases,delta_endpoint,1,response_plus,lower_factor=0.0_rp)
         call sr_l0_density_from_second_order_endpoint_branches(response_space,radial_bases,delta_endpoint,2,response_minus,lower_factor=0.0_rp)
         pauli_native_response_c = pauli_native_response_c + reciprocal_obj%k_weights(ik_global)*(response_plus+response_minus)

         if (source_complete) then
            call lr_static_frechet_density(reciprocal_obj%eigenvalues(:,ik),reciprocal_obj%eigenvectors(:,:,ik), &
               reciprocal_obj%fermi_level,reciprocal_obj%temperature,delta_alsda,delta_rho)
            call endpoint_tangent_branches_second_order(h,rho,delta_alsda,delta_rho,delta_endpoint)
            call sr_l0_density_from_second_order_endpoint_branches(response_space,radial_bases,delta_endpoint,1,response_plus)
            call sr_l0_density_from_second_order_endpoint_branches(response_space,radial_bases,delta_endpoint,2,response_minus)
            alsda_response_c = alsda_response_c + reciprocal_obj%k_weights(ik_global)*(response_plus+response_minus)
         end if
      end do
      field_weighted_rms = sqrt(field_weighted_rms/max(wsum,tiny(1.0_rp)))
      native_response_c = native_response_c/wsum
      pauli_native_response_c = pauli_native_response_c/wsum
      call sr_l0_density_from_second_order_endpoint_branches(response_space,radial_bases,delta_endpoint_sum,1,response_plus)
      call sr_l0_density_from_second_order_endpoint_branches(response_space,radial_bases,delta_endpoint_sum,2,response_minus)
      native_response_c = response_plus + response_minus
      call sr_l0_density_from_second_order_endpoint_branches(response_space,radial_bases,delta_endpoint_sum,1,response_plus,lower_factor=0.0_rp)
      call sr_l0_density_from_second_order_endpoint_branches(response_space,radial_bases,delta_endpoint_sum,2,response_minus,lower_factor=0.0_rp)
      pauli_native_response_c = response_plus + response_minus
      alsda_response_c = alsda_response_c/wsum
      call convert_response(native_response_c,ground_states,native_response)
      call convert_response(pauli_native_response_c,ground_states,pauli_native_response)
      call sr_aug_sigma_x(sigma_x)
      call dresp09y_augmentation_density(radial_bases, endpoint_sum, sigma_x, upper=augmentation_upper_c, &
         small=augmentation_small_c, angular=augmentation_angular_c, total=augmentation_response_c, &
         pauli_aug_c=pauli_augmentation_c, angular_rank0=augmentation_rank0_c, angular_rank2=augmentation_rank2_c)
      call convert_response(augmentation_upper_c,ground_states,augmentation_upper)
      call convert_response(augmentation_small_c,ground_states,augmentation_small)
      call convert_response(augmentation_angular_c,ground_states,augmentation_angular)
      call convert_response(augmentation_rank0_c,ground_states,augmentation_rank0)
      call convert_response(augmentation_rank2_c,ground_states,augmentation_rank2)
      call convert_response(augmentation_response_c,ground_states,augmentation_response)
      call convert_response(pauli_augmentation_c,ground_states,pauli_augmentation)
      complete_response = native_response + augmentation_response
      pauli_complete_response = pauli_native_response + pauli_augmentation
      call convert_response(alsda_response_c,ground_states,alsda_response)
      call radial_response_metric(native_response,target_sr,ground_states,native_fixed_residual)
      call radial_response_metric(pauli_native_response,target_p3,ground_states,native_pauli_residual)
      call radial_response_metric(pauli_complete_response,target_p3,ground_states,pauli_complete_residual)
      call radial_response_metric(augmentation_response,target_sr,ground_states,fixed_augmentation_residual)
      contact_rank_recomposition_residual = sqrt(sum((sr_aug_lower_spin_factor*augmentation_angular-&
         augmentation_rank0-augmentation_rank2)**2))/max(sqrt(sum((sr_aug_lower_spin_factor*augmentation_angular)**2)),tiny(1.0_rp))
      call radial_response_metric(complete_response,target_sr,ground_states,native_complete_residual)
      fixed_integral = volume_integral(ground_states(1),native_response(1,:))
      complete_integral = volume_integral(ground_states(1),complete_response(1,:))
      target_integral = volume_integral(ground_states(1),target_sr(1,:))
      fixed_volume_norm = volume_l2_norm(ground_states(1),native_response(1,:))
      complete_volume_norm = volume_l2_norm(ground_states(1),complete_response(1,:))
      target_volume_norm = volume_l2_norm(ground_states(1),target_sr(1,:))
      call radial_response_metric(alsda_response+augmentation_response,target_sr,ground_states,alsda_complete_residual)
      native_norm = sqrt(sum(complete_response**2)); target_norm = sqrt(sum(target_p3**2))

      allocate(chi_static(product_dim,product_dim),vertices(nmat,nmat,lmto_product_nbranch,product_dim), &
         trial_field(product_dim),action_matrix(product_dim),action_direct(product_dim),action_delta(product_dim))
      call lr_static_product_matrix(product,reciprocal_obj%eigenvalues,reciprocal_obj%eigenvectors,reciprocal_obj%k_weights, &
         reciprocal_obj%fermi_level,reciprocal_obj%temperature,chi_static)
      call product%component_vertex_tensor(vertices)
      do i = 1, product_dim
         trial_field(i) = cmplx(sin(real(i,rp)),cos(real(3*i,rp)),rp)
      end do
      action_matrix = matmul(chi_static,trial_field)
      action_direct = cmplx(0.0_rp,0.0_rp,rp)
      do ik = 1, nk
         ik_global = ik
         if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
         call dresp09_product_source_operator(vertices,reciprocal_obj%hk_bulk(:,:,ik),trial_field,delta_native)
         call lr_static_product_action(product,reciprocal_obj%eigenvalues(:,ik),reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level,reciprocal_obj%temperature,delta_native,action_delta)
         action_direct = action_direct + reciprocal_obj%k_weights(ik_global)*action_delta
      end do
      action_direct = action_direct/wsum
      static_action_residual = sqrt(sum(abs(action_direct-action_matrix)**2))/max(sqrt(sum(abs(action_matrix)**2)),tiny(1.0_rp))
      do i = 1, product_dim
         trial_field(i) = cmplx(sin(sqrt(2.0_rp)*real(i,rp)),cos(sqrt(3.0_rp)*real(i,rp)),rp)
      end do
      call static_action_residual_for_field(product,reciprocal_obj,vertices,chi_static,trial_field,static_pseudorandom_residual)
      call dresp09_compact_field_from_raw(response_space,product,source_field,trial_field)
      call static_action_residual_for_field(product,reciprocal_obj,vertices,chi_static,trial_field,static_physical_residual)

      field_closed = field_residual < 5.0e-8_rp
      augmentation_closed = native_complete_residual < 5.0e-8_rp .and. pauli_complete_residual < 5.0e-8_rp
      p3_closed = p3_identity_residual < 5.0e-12_rp .and. augmentation_closed
      frechet_closed = static_action_residual < 5.0e-10_rp
      mixed_operator_closed = .false.
      if (.not. frechet_closed) then
         verdict = dresp10_verdict_blocked
         classification = 'STATIC_CHI0_BRIDGE_OPEN'
      else if (.not. source_complete) then
         verdict = dresp10_verdict_blocked
         classification = 'SR_SOURCE_SPACE_INCOMPLETE'
      else if (.not. field_closed) then
         verdict = dresp10_verdict_blocked
         classification = 'FIELD_INSERTION_OPEN'
      else if (.not. augmentation_closed) then
         verdict = dresp10_verdict_blocked
         classification = 'AUGMENTATION_CONTACT_OPEN'
      else if (.not. p3_closed) then
         verdict = dresp10_verdict_blocked
         classification = 'RAW_ALSDA_FIELD_MISMATCH'
      else if (field_closed .and. p3_closed .and. frechet_closed .and. mixed_operator_closed) then
         verdict = dresp10_verdict_pass_b
         classification = 'RAW_FULL_SPATIAL_ALSDA_WARD_CLOSED'
      else if (field_closed .and. p3_closed .and. frechet_closed .and. .not. mixed_operator_closed) then
         verdict = dresp10_verdict_blocked
         classification = 'MIXED_OPERATOR_OPEN'
      else
         verdict = dresp10_verdict_blocked
         classification = 'RAW_ALSDA_WARD_DEFECT'
      end if

      if (rank == 0) then
         open(newunit=unit,file=output_file,status='replace',action='write',iostat=ios)
         if (ios /= 0) error stop 'DRESP-10: cannot open output artifact'
         write(unit,'(a)') '# DRESP-10 RAW full-spatial transverse ALSDA Ward artifact'
         write(unit,'(a,a)') 'verdict = ',trim(verdict)
         write(unit,'(a,a)') 'classification = ',trim(classification)
         write(unit,'(a)') 'source_space = complete six-branch SR upper/lower-small/lower-rank0/lower-rank2/total'
         write(unit,'(a,l1)') 'source_space_complete = ',source_complete
         write(unit,'(a,l1)') 'source_probe_nonzero = ',source_probe_complete
         write(unit,'(a,es24.16)') 'source_preflight_min_norm = ',min_source_norm
         write(unit,'(a,es24.16)') 'source_span_physical_global = ',source_span_physical_global
         write(unit,'(a,es24.16)') 'source_span_delta_O_global = ',source_span_delta_global
         write(unit,'(a)') 'L Pauli SR_upper SR_lower_small SR_lower_rank0 SR_lower_rank2 SR_total delta_O'
         do l = 0, response_space%response_lmax
            write(unit,'(i0,1x,7(es24.16,1x))') l, source_span_by_l(:,l)
         end do
         write(unit,'(a,i0)') 'response_lmax = ',response_space%response_lmax
         write(unit,'(a,i0)') 'product_dimension = ',product_dim
         write(unit,'(a)') 'product_branches = 00,10,01,11,20,02'
         write(unit,'(a)') 'target_primary = P3_PAULI_MOMENT_DENSITY'
         write(unit,'(a,es24.16)') 'P1_vs_P3_relative = ',sqrt(sum((p1-p3)**2))/max(sqrt(sum(p3**2)),tiny(1.0_rp))
         write(unit,'(a,es24.16)') 'P3_min_abs_active = ',min_abs_p3
         write(unit,'(a,es24.16)') 'P3_max_abs = ',max_abs_p3
         write(unit,'(a,es24.16)') 'Kxc_min_abs_active = ',min_abs_kxc
         write(unit,'(a,es24.16)') 'Kxc_max_abs_active = ',max_abs_kxc
         write(unit,'(a,es24.16)') 'Kxc_identity_max_abs = ',p3_identity_residual
         write(unit,'(a)') 'Kxc_low_m_status = DIRECT_P3_DENOMINATOR; no floor, fit, or residual-based denominator selection'
         write(unit,'(a,es24.16)') 'field_insertion_relative = ',field_residual
         write(unit,'(a,es24.16)') 'field_insertion_weighted_rms = ',field_weighted_rms
         write(unit,'(a,es24.16)') 'field_insertion_occupied_action = ',field_occupied_action
         write(unit,'(a,es24.16)') 'field_insertion_near_ef_action = ',field_near_ef_action
         write(unit,'(a,es24.16)') 'legacy_l0_field_insertion_relative = ',legacy_l0_field_residual
         write(unit,'(a,es24.16)') 'full_vs_legacy_l0_field_relative = ',full_vs_legacy_field_residual
         write(unit,'(a,es24.16)') 'field_insertion_max_element = ',field_max_element
         write(unit,'(a,es24.16)') 'source_piece_upper_norm = ',source_upper_norm
         write(unit,'(a,es24.16)') 'source_piece_lower_small_norm = ',source_small_norm
         write(unit,'(a,es24.16)') 'source_piece_lower_rank0_norm = ',source_rank0_norm
         write(unit,'(a,es24.16)') 'source_piece_lower_rank2_norm = ',source_rank2_norm
         write(unit,'(a,es24.16)') 'source_piece_total_norm = ',source_total_norm
         write(unit,'(a,es24.16)') 'contact_upper_volume_norm = ',volume_l2_norm(ground_states(1),augmentation_upper(1,:))
         write(unit,'(a,es24.16)') 'contact_lower_small_volume_norm = ',volume_l2_norm(ground_states(1),augmentation_small(1,:))
         write(unit,'(a,es24.16)') 'contact_lower_angular_volume_norm = ',volume_l2_norm(ground_states(1),augmentation_angular(1,:))
         write(unit,'(a,es24.16)') 'contact_lower_angular_rank0_volume_norm = ',volume_l2_norm(ground_states(1),augmentation_rank0(1,:))
         write(unit,'(a,es24.16)') 'contact_lower_angular_rank2_volume_norm = ',volume_l2_norm(ground_states(1),augmentation_rank2(1,:))
         write(unit,'(a,es24.16)') 'contact_angular_rank_recomposition_relative = ',contact_rank_recomposition_residual
         write(unit,'(a,es24.16)') 'contact_total_volume_norm = ',volume_l2_norm(ground_states(1),augmentation_response(1,:))
         write(unit,'(a,es24.16)') 'fixed_observable_vs_SRspin2 = ',native_fixed_residual
         write(unit,'(a,es24.16)') 'augmentation_vs_SRspin2 = ',fixed_augmentation_residual
         write(unit,'(a,es24.16)') 'complete_response_vs_SRspin2 = ',native_complete_residual
         write(unit,'(a,es24.16)') 'fixed_response_integrated = ',fixed_integral
         write(unit,'(a,es24.16)') 'complete_response_integrated = ',complete_integral
         write(unit,'(a,es24.16)') 'target_SRspin2_integrated = ',target_integral
         write(unit,'(a,es24.16)') 'fixed_response_volume_norm = ',fixed_volume_norm
         write(unit,'(a,es24.16)') 'complete_response_volume_norm = ',complete_volume_norm
         write(unit,'(a,es24.16)') 'target_SRspin2_volume_norm = ',target_volume_norm
         write(unit,'(a,es24.16)') 'Pauli_fixed_vs_P3 = ',native_pauli_residual
         write(unit,'(a,es24.16)') 'Pauli_complete_vs_P3 = ',pauli_complete_residual
         if (augmentation_closed) then
            write(unit,'(a)') 'DRESP09Y_contact_regression = CLOSED'
         else
            write(unit,'(a)') 'DRESP09Y_contact_regression = OPEN'
         end if
         write(unit,'(a,es24.16)') 'native_plus_contact_vs_P3 = ',pauli_complete_residual
         write(unit,'(a,es24.16)') 'native_plus_contact_vs_SRspin2 = ',native_complete_residual
         if (source_complete) then
            write(unit,'(a,es24.16)') 'ALSDA_plus_contact_vs_SRspin2 = ',alsda_complete_residual
            write(unit,'(a,es24.16)') 'ALSDA_field_vs_native_relative = ',field_residual
         else
            write(unit,'(a)') 'ALSDA_plus_contact_vs_SRspin2 = NOT_RUN_BEFORE_SOURCE_SPAN_GATE'
            write(unit,'(a)') 'ALSDA_field_vs_native_relative = NOT_INTERPRETED_BEFORE_SOURCE_SPAN_GATE'
         end if
         write(unit,'(a,es24.16)') 'static_frechet_action_relative = ',static_action_residual
         write(unit,'(a,es24.16)') 'static_frechet_pseudorandom_relative = ',static_pseudorandom_residual
         write(unit,'(a,es24.16)') 'static_frechet_physical_source_relative = ',static_physical_residual
         write(unit,'(a,es24.16)') 'native_frechet_vs_commutator_relative = ',frechet_native_residual
         write(unit,'(a)') 'static_frechet_matrix = direct finite-temperature Frechet; independent of compact chi0'
         write(unit,'(a)') 'mixed_K_contact_operator = NOT_ASSEMBLED'
         write(unit,'(a)') 'raw_denominator = NOT_ASSEMBLED; no correction or Goldstone repair applied'
         write(unit,'(a)') 'ALSDA_Ward = NOT_INTERPRETED_BEFORE_PRE_ALSDA_GATES'
         write(unit,'(a)') 'BES_Halle = OFF'
         write(unit,'(a)') 'Goldstone_correction = OFF'
         write(unit,'(a)') 'dynamic_spectra = NOT_RUN'
         write(unit,'(a)') 'eta_secondary_ladder = NOT_USED_BY_STATIC_ORACLE'
         close(unit)
      end if

      deallocate(moment_sum,endpoint_sum,h,rho,moments,endpoints,generator,delta_native,delta_alsda,delta_rho,source,source_pieces, &
         source_pieces_minus,source_components,source_components_minus, &
         delta_endpoint,delta_endpoint_sum,channels,p1,p1_raw,p3,target_p3,target_sr,target_sr_raw,p3_l,bxc,kxc,source_field,native_response,alsda_response, &
         augmentation_upper,augmentation_small,augmentation_angular,augmentation_response,complete_response,pauli_native_response, &
         augmentation_rank0,augmentation_rank2, &
         pauli_augmentation,pauli_complete_response,response_plus,response_minus,native_response_c,pauli_native_response_c, &
         alsda_response_c,augmentation_upper_c,augmentation_small_c,augmentation_angular_c,augmentation_rank0_c,augmentation_rank2_c, &
         augmentation_response_c, &
         pauli_augmentation_c,sigma_x,chi_static,vertices,trial_field,action_matrix,action_direct,action_delta)
   end subroutine run_dresp10_raw_full_spatial_alsda_ward

   subroutine static_action_residual_for_field(product,reciprocal_obj,vertices,chi_static,field,value)
      type(lmto_product_response_basis), intent(in) :: product
      type(reciprocal), intent(in) :: reciprocal_obj
      complex(rp), intent(in) :: vertices(:, :, :, :), chi_static(:, :), field(:)
      real(rp), intent(out) :: value
      complex(rp), allocatable :: action_matrix(:), action_direct(:), action_delta(:), delta_h(:, :)
      integer :: ik, ik_global
      real(rp) :: wsum

      allocate(action_matrix(size(field)),action_direct(size(field)),action_delta(size(field)), &
         delta_h(size(reciprocal_obj%hk_bulk,1),size(reciprocal_obj%hk_bulk,2)))
      action_matrix = matmul(chi_static,field)
      action_direct = cmplx(0.0_rp,0.0_rp,rp)
      wsum = sum(reciprocal_obj%k_weights)
      do ik = 1, size(reciprocal_obj%hk_bulk,3)
         ik_global = ik
         if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
         call dresp09_product_source_operator(vertices,reciprocal_obj%hk_bulk(:,:,ik),field,delta_h)
         call lr_static_product_action(product,reciprocal_obj%eigenvalues(:,ik),reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level,reciprocal_obj%temperature,delta_h,action_delta)
         action_direct = action_direct + reciprocal_obj%k_weights(ik_global)*action_delta
      end do
      action_direct = action_direct/wsum
      value = sqrt(sum(abs(action_direct-action_matrix)**2))/max(sqrt(sum(abs(action_matrix)**2)),tiny(1.0_rp))
      deallocate(action_matrix,action_direct,action_delta,delta_h)
   end subroutine static_action_residual_for_field

   subroutine convert_response(complex_density,states,weighted)
      complex(rp), intent(in) :: complex_density(:, :)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(out) :: weighted(:, :)
      integer :: site, ir
      weighted = 0.0_rp
      do site = 1, size(states)
         do ir = 2, size(states(site)%r)
            weighted(site,ir) = sqrt(dresp09w_four_pi)*states(site)%r(ir)**2*real(complex_density(site,ir),rp)
         end do
      end do
   end subroutine convert_response

   subroutine action_metrics_local(diff, reference, eigenvectors, eigenvalues, fermi_level, occupied_value, near_ef_value, &
                                   noccupied, nnear)
      complex(rp), intent(in) :: diff(:, :), reference(:, :), eigenvectors(:, :)
      real(rp), intent(in) :: eigenvalues(:), fermi_level
      real(rp), intent(out) :: occupied_value, near_ef_value
      integer, intent(out) :: noccupied, nnear
      integer :: ib
      real(rp) :: numerator, denominator, ratio

      occupied_value = 0.0_rp; near_ef_value = 0.0_rp; noccupied = 0; nnear = 0
      do ib = 1, size(eigenvalues)
         numerator = sqrt(sum(abs(matmul(diff,eigenvectors(:,ib)))**2))
         denominator = max(sqrt(sum(abs(matmul(reference,eigenvectors(:,ib)))**2)),tiny(1.0_rp))
         ratio = numerator/denominator
         if (eigenvalues(ib) <= fermi_level) then
            occupied_value = max(occupied_value,ratio); noccupied = noccupied + 1
         end if
         if (abs(eigenvalues(ib)-fermi_level) <= 0.10_rp) then
            near_ef_value = max(near_ef_value,ratio); nnear = nnear + 1
         end if
      end do
   end subroutine action_metrics_local

   subroutine radial_response_metric(left,right,states,value)
      real(rp), intent(in) :: left(:, :), right(:, :)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(out) :: value
      real(rp) :: l2, max_absolute, max_relative, integral
      call radial_relative_metrics(left,right,states,value,l2,max_absolute,max_relative,integral)
   end subroutine radial_response_metric

end module lr_dresp10_ward_bridge_mod
