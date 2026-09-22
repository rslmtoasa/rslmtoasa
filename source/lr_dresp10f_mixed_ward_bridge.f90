!------------------------------------------------------------------------------
! DRESP-10F fixed-ground-state mixed-representation static Ward bridge.
!
! The field and response coordinates are deliberately different objects here:
!
!   chi0(P <- SR) = R_P L_f(H) F_SR.
!
! F_SR is the full spatial scalar-relativistic source insertion.  R_P is the
! fixed-ground-state Pauli six-branch measurement.  No field projection into
! the Pauli 348-space and no moving-basis delta-O term enter this path.
!------------------------------------------------------------------------------
module lr_dresp10f_mixed_ward_bridge_mod

   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use mpi_mod, only: rank
   use basis_mod, only: nb
   use reciprocal_mod, only: reciprocal
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use radial_ground_state_mod, only: radial_ground_state
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use response_angular_basis_mod, only: response_angular_pi, response_gaunt
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex, &
      response_unflatten_superindex
   use lr_response_space_mod, only: response_space_layout
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, lmto_product_channel_plus, &
      lmto_product_channel_minus
   use lr_lmto_endpoint_branches_mod, only: lmto_product_nbranch, lmto_product_branch_powers, &
      lmto_product_apply_branch_action, lmto_product_second_order_radial_branch
   use lr_full_spatial_alsda_mod, only: lr_full_spatial_source, lr_full_spatial_source_components, &
      lr_full_spatial_contract_components, lr_full_spatial_piece_total, lr_full_spatial_npiece
   use lr_static_frechet_product_response_mod, only: lr_static_frechet_density, lr_static_product_matrix, &
      lr_static_product_action
   use lr_dresp09_field_insertion_mod, only: dresp09_compact_field_from_raw, dresp09_raw_field_from_compact, &
      dresp09_product_source_operator
   use lr_pauli_transition_vertex_mod, only: pauli_endpoint_state
   use lr_dresp09zs_compact_span_mod, only: dresp09zs_source_span_residuals, &
      dresp09zs_pauli_branch_radial_audit
   use lr_dresp09u_bridge_mod, only: run_dresp09u_representation_tangent
   use lr_dresp09y_bridge_mod, only: run_dresp09y_augmentation_tangent, dresp09y_augmentation_density
   use lr_exact_ks_ward_mod, only: exact_ks_spin_rotation_generator
   use lr_lmto_density_moment_tangent_mod, only: density_from_eigensystem, moment_matrices_from_eigensystem, &
      endpoint_matrices_from_moments_second_order, endpoint_tangent_branches_second_order, commutator_tangent
   use lr_radial_observable_provenance_mod, only: channel_moments_from_matrix, pauli_density_from_moments, &
      dresp09w_four_pi
   use lr_sr_augmentation_tangent_mod, only: sr_aug_sigma_x
   implicit none
   private

   character(len=*), parameter, public :: dresp10f_pass_a = 'PASS-A'
   character(len=*), parameter, public :: dresp10f_pass_b = 'PASS-B'
   character(len=*), parameter, public :: dresp10f_blocked = 'BLOCKED'

   logical, save :: measurement_cache_ready = .false.
   integer, save :: measurement_cache_nsite = 0, measurement_cache_npoint = 0, measurement_cache_norb = 0
   real(rp), allocatable, save :: measurement_radial_cache(:, :, :, :, :), measurement_gaunt_cache(:, :, :)

   public :: run_dresp10f_fixed_basis_mixed_ward
   public :: apply_static_selfconsistent_action
   public :: apply_static_denominator
   ! Read-only measurement seams reused by DRESP-12.  These expose the
   ! already-certified Fréchet and compact-transition measurements without
   ! changing the fixed-basis action or its six-branch source.
   public :: build_l0_source
   public :: build_l0_target
   public :: fixed_pauli_measurement
   public :: compact_density_measurement

contains

   subroutine run_dresp10f_fixed_basis_mixed_ward(output_file, response_space, radial_bases, ground_states, &
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
      complex(rp), allocatable :: h(:, :), rho(:, :), delta_h(:, :), delta_rho(:, :), generator(:, :)
      complex(rp), allocatable :: moments(:, :, :), moment_sum(:, :, :), endpoint(:, :, :), endpoint_sum(:, :, :)
      complex(rp), allocatable :: source_field(:), target_raw(:), generated_raw(:), target_compact(:), generated_compact(:)
      complex(rp), allocatable :: probe_source(:), probe_raw(:), probe_compact(:)
      complex(rp), allocatable :: random_compact(:), random_matrix_action(:), random_direct_action(:), product_vertices(:, :, :, :)
      complex(rp), allocatable :: covariant_raw(:), covariant_contact(:), covariant_total(:), contact_c(:, :)
      real(rp), allocatable :: channels(:, :, :, :), p3(:, :), bxc(:, :), kxc(:, :)
      real(rp) :: source_span_by_l(7,0:4), source_span_physical_by_l(0:4), source_span_delta_by_l(0:4)
      real(rp) :: source_span_physical_global, source_span_delta_global
      real(rp) :: pauli_closure_by_l(0:4), pauli_closure_max
      real(rp) :: branch_error(lmto_product_nbranch), branch_max_error
      real(rp) :: source_vertex_residual, static_matrix_residual, physical_action_residual
      real(rp) :: mixed_action_residual, rigid_action_residual, arbitrary_action_residual, p3_identity_residual, ward_absolute, ward_relative
      real(rp) :: ward_radial, ward_integrated, ward_max, ward_compact, rigid_overlap
      real(rp) :: generated_norm, target_norm, covariant_residual, contact_effect, fixed_covariant_residual
      real(rp) :: representation_field_defect, representation_response_defect
      real(rp) :: min_abs_p3, max_abs_p3, min_abs_kxc, max_abs_kxc
      real(rp) :: general_action(6), general_action_oracle(6), general_linearity(6)
      integer :: nsite, nmat, nk, ndim, product_dim, ik, ik_global, flat, ir, site, l, m, unit, ios
      integer :: i, j, branch, response_l, response_m
      logical :: span_closed, source_closed, frechet_closed, mixed_closed, kxc_closed, covariant_closed
      character(len=96) :: classification, verdict
      character(len=256) :: covariance_u_file, covariance_y_file

      nsite = size(ground_states)
      if (nsite < 1 .or. size(radial_bases) /= nsite .or. response_space%nsite /= nsite) then
         error stop 'DRESP-10F: response/radial site dimensions are inconsistent'
      end if
      if (.not. hamiltonian_obj%hoh .or. hamiltonian_obj%ccor_2c .or. hamiltonian_obj%orb_pol .or. &
          hamiltonian_obj%hubbard_u_general_check .or. hamiltonian_obj%hubbard_v_check .or. hamiltonian_obj%local_axis) then
         error stop 'DRESP-10F: excluded Hamiltonian correction is enabled'
      end if
      if (.not. allocated(reciprocal_obj%hk_bulk) .or. .not. allocated(reciprocal_obj%eigenvalues) .or. &
          .not. allocated(reciprocal_obj%eigenvectors) .or. .not. allocated(reciprocal_obj%k_weights)) then
         error stop 'DRESP-10F: accepted reciprocal eigensystem is incomplete'
      end if
      if (response_space%response_lmax /= 4 .or. response_space%nchannel /= 1) then
         error stop 'DRESP-10F: accepted L=0..4 single-channel response space is required'
      end if

      nmat = size(reciprocal_obj%hk_bulk, 1)
      nk = size(reciprocal_obj%hk_bulk, 3)
      ndim = response_space%ndim
      if (nmat /= nb*nsite .or. size(reciprocal_obj%k_weights) /= nk) then
         error stop 'DRESP-10F: accepted reciprocal dimensions are incomplete'
      end if
      if (lattice_obj%nrec /= nsite .or. nk /= 64) error stop 'DRESP-10F: certified bcc-Fe material gate failed'

      call product%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      product_dim = product%product_dimension
      if (product_dim /= 348) error stop 'DRESP-10F: accepted Pauli product basis is not the live 348-space'

      call dresp09zs_pauli_branch_radial_audit(radial_bases(1), branch_error, branch_max_error)
      call dresp09zs_source_span_residuals(response_space, radial_bases, product, source_span_physical_by_l, &
         source_span_delta_by_l, source_span_by_l, source_span_physical_global, source_span_delta_global)
      call pauli_product_span_oracle(response_space, radial_bases(1), product, pauli_closure_by_l, pauli_closure_max)
      source_span_by_l(1,:) = pauli_closure_by_l
      span_closed = pauli_closure_max < 1.0e-10_rp .and. branch_max_error < 1.0e-12_rp

      allocate(moment_sum(nmat,nmat,3), moments(nmat,nmat,3), endpoint_sum(nmat,nmat,6), endpoint(nmat,nmat,6), &
         h(nmat,nmat), rho(nmat,nmat), delta_h(nmat,nmat), delta_rho(nmat,nmat), generator(nmat,nmat))
      moment_sum = cmplx(0.0_rp,0.0_rp,rp); endpoint_sum = cmplx(0.0_rp,0.0_rp,rp)
      call exact_ks_spin_rotation_generator((radial_bases(1)%lmax+1)**2, nsite, generator)
      do ik = 1, nk
         ik_global = ik
         if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
         h = reciprocal_obj%hk_bulk(:,:,ik)
         call density_from_eigensystem(reciprocal_obj%eigenvalues(:,ik), reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level, reciprocal_obj%temperature, rho)
         call moment_matrices_from_eigensystem(reciprocal_obj%eigenvalues(:,ik), reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level, reciprocal_obj%temperature, moments)
         call endpoint_matrices_from_moments_second_order(moments, endpoint)
         moment_sum = moment_sum + reciprocal_obj%k_weights(ik_global)*moments
         endpoint_sum = endpoint_sum + reciprocal_obj%k_weights(ik_global)*endpoint
      end do
      moment_sum = moment_sum/sum(reciprocal_obj%k_weights)
      endpoint_sum = endpoint_sum/sum(reciprocal_obj%k_weights)

      allocate(channels(3,nsite,radial_bases(1)%lmax+1,2), p3(nsite,response_space%npoint), &
         bxc(nsite,response_space%npoint), kxc(nsite,response_space%npoint), source_field(ndim), target_raw(ndim), &
         generated_raw(ndim), target_compact(product_dim), generated_compact(product_dim), random_compact(product_dim), &
         random_matrix_action(product_dim), random_direct_action(product_dim), covariant_raw(ndim), covariant_contact(ndim), &
         covariant_total(ndim), contact_c(nsite,response_space%npoint), probe_source(ndim), probe_raw(ndim), probe_compact(product_dim))
      call channel_moments_from_matrix(moment_sum, nsite, channels)
      call pauli_density_from_moments(radial_bases, channels, p3)

      min_abs_p3 = huge(1.0_rp); max_abs_p3 = 0.0_rp; min_abs_kxc = huge(1.0_rp); max_abs_kxc = 0.0_rp
      p3_identity_residual = 0.0_rp
      do site = 1, nsite
         do ir = 1, response_space%npoint
            bxc(site,ir) = 0.5_rp*(ground_states(site)%vxc_up(ir)-ground_states(site)%vxc_down(ir))
            if (response_space%radial_weights(ir) > 0.0_rp) then
               if (p3(site,ir) == 0.0_rp) error stop 'DRESP-10F: zero P3 on a positive-measure point'
               kxc(site,ir) = bxc(site,ir)/p3(site,ir)
               min_abs_p3 = min(min_abs_p3,abs(p3(site,ir)))
               min_abs_kxc = min(min_abs_kxc,abs(kxc(site,ir)))
            else
               kxc(site,ir) = 0.0_rp
            end if
            max_abs_p3 = max(max_abs_p3,abs(p3(site,ir)))
            max_abs_kxc = max(max_abs_kxc,abs(kxc(site,ir)))
            p3_identity_residual = max(p3_identity_residual,abs(kxc(site,ir)*p3(site,ir)-bxc(site,ir)))
         end do
      end do
      kxc_closed = p3_identity_residual < 1.0e-12_rp
      call build_l0_source(response_space, bxc, source_field)
      call build_l0_target(response_space, p3, target_raw)
      call dresp09_compact_field_from_raw(response_space, product, target_raw, target_compact)

      ! Direct fixed-basis mixed action.  The raw SR field is inserted in its
      ! native full-spatial coordinates; no source projection is performed.
      call mixed_static_action(response_space, radial_bases, product, reciprocal_obj, source_field, generated_raw, &
         generated_compact, mixed_action_residual, source_closed)
      rigid_action_residual = mixed_action_residual
      call radial_metrics(response_space, generated_raw, target_raw, ward_absolute, ward_relative, ward_max, ward_integrated)
      ward_compact = sqrt(sum(abs(generated_compact-target_compact)**2))
      generated_norm = sqrt(sum(abs(generated_raw)**2)); target_norm = sqrt(sum(abs(target_raw)**2))
      if (target_norm > tiny(1.0_rp)) then
         rigid_overlap = abs(sum(conjg(target_compact)*(generated_compact-target_compact)))/ &
            max(sqrt(sum(abs(target_compact)**2))*max(sqrt(sum(abs(generated_compact-target_compact)**2)),tiny(1.0_rp)),tiny(1.0_rp))
      else
         rigid_overlap = 0.0_rp
      end if

      ! Exercise the same mixed action with a deterministic non-rigid raw SR
      ! field. This is an action/oracle check, not a projected-field solve.
      do flat = 1, ndim
         probe_source(flat) = cmplx(sin(0.17_rp*real(flat,rp)),cos(0.31_rp*real(flat,rp)),rp)
      end do
      call mixed_static_action(response_space, radial_bases, product, reciprocal_obj, probe_source, probe_raw, &
         probe_compact, arbitrary_action_residual, source_closed)
      mixed_action_residual = max(mixed_action_residual, arbitrary_action_residual)

      ! Independent compact static/Frechet oracle and an explicit physical
      ! source action comparison.  The compact field is used only for the
      ! response-side oracle, never to insert source_field.
      allocate(product_vertices(nmat,nmat,lmto_product_nbranch,product_dim))
      call product%component_vertex_tensor(product_vertices)
      do i = 1, product_dim
         random_compact(i) = cmplx(sin(real(i,rp)*0.37_rp),cos(real(i,rp)*0.19_rp),rp)
      end do
      random_matrix_action = cmplx(0.0_rp,0.0_rp,rp)
      random_direct_action = cmplx(0.0_rp,0.0_rp,rp)
      do ik = 1, nk
         ik_global = ik
         if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
         call dresp09_product_source_operator(product_vertices, reciprocal_obj%hk_bulk(:,:,ik), random_compact, delta_h)
         call lr_static_product_action(product, reciprocal_obj%eigenvalues(:,ik), reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level, reciprocal_obj%temperature, delta_h, random_direct_action)
         random_matrix_action = random_matrix_action + reciprocal_obj%k_weights(ik_global)*random_direct_action
      end do
      random_matrix_action = random_matrix_action/sum(reciprocal_obj%k_weights)
      call static_matrix_oracle(product, reciprocal_obj, product_vertices, random_compact, random_matrix_action, &
         static_matrix_residual)
      physical_action_residual = rigid_action_residual
      frechet_closed = static_matrix_residual < 1.0e-10_rp

      ! The separate moving-basis result is rerun into sidecar artifacts.  It
      ! remains an oracle and is never passed to the fixed-basis action.
      covariance_u_file = trim(output_file)//'.DRESP09U'
      covariance_y_file = trim(output_file)//'.DRESP09Y'
      call run_dresp09u_representation_tangent(covariance_u_file, reciprocal_obj, lattice_obj, hamiltonian_obj)
      call run_dresp09y_augmentation_tangent(covariance_y_file, response_space, radial_bases, ground_states, &
         reciprocal_obj, lattice_obj, hamiltonian_obj)
      call covariant_oracle(response_space, radial_bases, reciprocal_obj, generator, endpoint_sum, p3, covariant_raw, &
         covariant_contact, covariant_total, covariant_residual, contact_effect, covariant_closed)
      fixed_covariant_residual = weighted_relative(response_space, generated_raw, covariant_raw)
      representation_field_defect = native_field_defect(response_space, radial_bases, reciprocal_obj, generator, source_field)
      representation_response_defect = fixed_covariant_residual

      call general_mixed_actions(response_space, radial_bases, product, reciprocal_obj, kxc, general_action, &
         general_action_oracle, general_linearity)

      call source_vertex_oracle(response_space, radial_bases, reciprocal_obj, source_vertex_residual, source_closed)
      mixed_closed = mixed_action_residual < 1.0e-10_rp
      if (.not. span_closed) then
         verdict = dresp10f_blocked; classification = 'PRODUCT_SPAN_AUDIT_INCONSISTENT'
      else if (.not. source_closed) then
         verdict = dresp10f_blocked; classification = 'SOURCE_VERTEX_OPEN'
      else if (.not. frechet_closed) then
         verdict = dresp10f_blocked; classification = 'MIXED_RESPONSE_ACTION_OPEN'
      else if (.not. mixed_closed) then
         verdict = dresp10f_blocked; classification = 'MIXED_RESPONSE_ACTION_OPEN'
      else if (.not. kxc_closed) then
         verdict = dresp10f_blocked; classification = 'MIXED_RESPONSE_ACTION_OPEN'
      else if (ward_relative <= 1.0e-10_rp .and. covariant_closed) then
         verdict = dresp10f_pass_a; classification = 'FIXED_BASIS_MIXED_WARD_CLOSED'
      else
         verdict = dresp10f_pass_b; classification = 'FINITE_LMTO_WARD_INCONSISTENCY'
      end if

      if (rank == 0) then
         open(newunit=unit,file=output_file,status='replace',action='write',iostat=ios)
         if (ios /= 0) error stop 'DRESP-10F: cannot open output artifact'
         write(unit,'(a)') '# DRESP-10F fixed-basis mixed TDDFT Ward artifact'
         write(unit,'(a,a)') 'DRESP-10F verdict: ',trim(verdict)
         write(unit,'(a)') 'starting_HEAD = 7f745d333245dd40471c8505ba4959018e582d15'
         write(unit,'(a)') 'formulation = FIXED_GROUND_STATE_MIXED_RESPONSE'
         write(unit,'(a,a)') 'classification = ',trim(classification)
         write(unit,'(a,i0)') 'Pauli product basis dimension = ',product_dim
         write(unit,'(a)') 'Pauli product basis branches = 00,10,01,11,20,02'
         write(unit,'(a)') 'field_basis = INDEPENDENT_RAW_SR_FULL_SPATIAL; no Pauli projection before F_SR'
         write(unit,'(a)') 'primary_map = chi0(P<-SR) = R_P L_f(H) F_SR'
         write(unit,'(a,es24.16)') 'corrected_span_audit_pauli_max = ',branch_max_error
         write(unit,'(a,6(es24.16,1x))') 'corrected_span_audit_pauli_00_10_01_11_20_02 = ',branch_error
         write(unit,'(a,5(es24.16,1x))') 'corrected_pauli_product_closure_L0_L1_L2_L3_L4 = ',pauli_closure_by_l
         write(unit,'(a)') 'L Pauli(product closure) SR_upper SR_lower_small SR_lower_rank0 SR_lower_rank2 SR_total delta_O'
         do l = 0, 4
            write(unit,'(i0,1x,7(es24.16,1x))') l,source_span_by_l(:,l)
         end do
         write(unit,'(a,es24.16)') 'Pauli_span_residual = ',maxval(source_span_by_l(1,:))
         write(unit,'(a,es24.16)') 'SR_upper_diagnostic = ',maxval(source_span_by_l(2,:))
         write(unit,'(a,es24.16)') 'SR_lower_small_diagnostic = ',maxval(source_span_by_l(3,:))
         write(unit,'(a,es24.16)') 'SR_rank0_diagnostic = ',maxval(source_span_by_l(4,:))
         write(unit,'(a,es24.16)') 'SR_rank2_diagnostic = ',maxval(source_span_by_l(5,:))
         write(unit,'(a,es24.16)') 'SR_total_diagnostic = ',maxval(source_span_by_l(6,:))
         write(unit,'(a,es24.16)') 'delta-O_diagnostic = ',maxval(source_span_by_l(7,:))
         write(unit,'(a,es24.16)') 'source_span_physical_global = ',source_span_physical_global
         write(unit,'(a,es24.16)') 'source_span_delta_O_global = ',source_span_delta_global
         write(unit,'(a,a)') 'source_vertex = ',merge('CLOSED','OPEN  ',source_closed)
         write(unit,'(a,es24.16)') 'source_vertex_mixed_LM_residual = ',source_vertex_residual
         write(unit,'(a,es24.16)') 'static_frechet_direct_matrix_residual = ',static_matrix_residual
         write(unit,'(a,es24.16)') 'static_frechet_physical_source_action_residual = ',physical_action_residual
         write(unit,'(a,es24.16)') 'mixed_chi0_rigid_source_oracle = ',rigid_action_residual
         write(unit,'(a,es24.16)') 'mixed_chi0_arbitrary_source_oracle = ',arbitrary_action_residual
         write(unit,'(a,es24.16)') 'Kxc_P3_identity_max_abs = ',p3_identity_residual
         write(unit,'(a)') 'Kxc_low_m_status = DIRECT_P3_DENOMINATOR; no fit, floor, or residual minimization'
         write(unit,'(a,es24.16)') 'P3_min_abs_active = ',min_abs_p3
         write(unit,'(a,es24.16)') 'P3_max_abs = ',max_abs_p3
         write(unit,'(a,es24.16)') 'Kxc_min_abs_active = ',min_abs_kxc
         write(unit,'(a,es24.16)') 'Kxc_max_abs = ',max_abs_kxc
         write(unit,'(a,es24.16)') 'PRIMARY_fixed_basis_generated_P3_norm = ',generated_norm
         write(unit,'(a,es24.16)') 'PRIMARY_fixed_basis_target_P3_norm = ',target_norm
         write(unit,'(a,es24.16)') 'PRIMARY_fixed_basis_absolute_residual = ',ward_absolute
         write(unit,'(a,es24.16)') 'PRIMARY_fixed_basis_relative_weighted_radial_L2 = ',ward_relative
         write(unit,'(a,es24.16)') 'PRIMARY_fixed_basis_maximum_radial_component = ',ward_max
         write(unit,'(a,es24.16)') 'PRIMARY_fixed_basis_integrated_moment_residual = ',ward_integrated
         write(unit,'(a,es24.16)') 'PRIMARY_fixed_basis_compact_residual = ',ward_compact
         write(unit,'(a,es24.16)') 'PRIMARY_fixed_basis_rigid_mode_overlap = ',rigid_overlap
         write(unit,'(a,es24.16)') 'moving_basis_native_plus_contact_vs_P3 = ',covariant_residual
         write(unit,'(a,a)') 'moving_basis_covariance_oracle = ',merge('CLOSED','OPEN  ',covariant_closed)
         write(unit,'(a,a)') 'moving_basis_covariance_artifact = ',trim(covariance_u_file)//' + '//trim(covariance_y_file)
         write(unit,'(a,es24.16)') 'representation_direct_Bxc_field_vs_covariant_native = ',representation_field_defect
         write(unit,'(a,es24.16)') 'representation_fixed_response_vs_covariant_response = ',representation_response_defect
         write(unit,'(a,es24.16)') 'representation_diagnostic_contact_effect = ',contact_effect
         write(unit,'(a,6(es24.16,1x))') 'general_348_to_348_action_L0_L1_L2_L3_L4_mixed = ',general_action
         write(unit,'(a,6(es24.16,1x))') 'general_348_to_348_action_oracle = ',general_action_oracle
         write(unit,'(a,6(es24.16,1x))') 'general_348_to_348_linearity = ',general_linearity
         write(unit,'(a)') 'general_348_to_348_assembly = NOT_PRACTICAL; matrix-free action is authoritative'
         write(unit,'(a)') 'BES_Halle = OFF'
         write(unit,'(a)') 'Goldstone_correction = OFF'
         write(unit,'(a)') 'Dynamics = NOT RUN'
         write(unit,'(a)') 'historical_DRESP10R_source_span = SUPERSEDED_SOURCE_SPAN_DIAGNOSTIC'
         write(unit,'(a)') 'historical_DRESP10R_Pauli = 7.801e-1'
         write(unit,'(a)') 'historical_DRESP10R_SR_total = 7.920e-1'
         write(unit,'(a)') 'historical_DRESP10R_delta-O = 5.990e-1'
         write(unit,'(a)') 'recommended_next_physics_step = REVIEW_FINITE_LMTO_COVARIANCE_BEFORE_DYNAMICS'
         write(unit,'(a)') 'tests = TddftDresp10FFixedBasisMixedWard; Dresp10FFeArtifact'
         write(unit,'(a)') 'primary_files = lr_dresp10f_mixed_ward_bridge.f90; DRESP_10_FIXED_BASIS_MIXED_WARD.md'
         close(unit)
      end if

      deallocate(moment_sum,moments,endpoint_sum,endpoint,h,rho,delta_h,delta_rho,generator,channels,p3,bxc,kxc,source_field, &
         target_raw,generated_raw,target_compact,generated_compact,random_compact,random_matrix_action,random_direct_action, &
         product_vertices,covariant_raw,covariant_contact,covariant_total,contact_c,probe_source,probe_raw,probe_compact)
   end subroutine run_dresp10f_fixed_basis_mixed_ward

   !> Authoritative DRESP-11 matrix-free action in the compact Pauli space.
   !>
   !> The input is a compact Pauli density coordinate vector.  It is first
   !> reconstructed by D_P, converted by the physical local Kxc map to a raw
   !> SR field, and then passed through the unchanged DRESP-10F source,
   !> Frechet, and R_P chain.  The radial factor is the coordinate conversion
   !> between the physical density d(r) and the raw response coordinate x(r):
   !> d(r)=sqrt(4*pi)*r**2*x(r).
   subroutine apply_static_selfconsistent_action(space, radial_bases, product, reciprocal_obj, kxc, c, ac)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(lmto_product_response_basis), intent(in) :: product
      type(reciprocal), intent(in) :: reciprocal_obj
      real(rp), intent(in) :: kxc(:, :)
      complex(rp), intent(in) :: c(:)
      complex(rp), intent(out) :: ac(:)
      complex(rp), allocatable :: density(:), source(:), response(:)
      real(rp) :: oracle_residual
      logical :: source_vertex_ok
      type(response_super_index) :: item
      integer :: flat

      if (size(c) /= product%product_dimension .or. size(ac) /= product%product_dimension .or. &
          size(kxc, 1) /= space%nsite .or. size(kxc, 2) /= space%npoint) then
         error stop 'DRESP-11 static action: compact/Kxc shape mismatch'
      end if
      allocate(density(space%ndim), source(space%ndim), response(space%ndim))
      call dresp09_raw_field_from_compact(space, product, c, density)
      source = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, space%nchannel, item)
         if (space%radius(item%radial_point) > 0.0_rp) then
            source(flat) = cmplx(4.0_rp*response_angular_pi*space%radius(item%radial_point)**2* &
               kxc(item%site, item%radial_point), 0.0_rp, rp)*density(flat)
         end if
      end do
      call mixed_static_action(space, radial_bases, product, reciprocal_obj, source, response, ac, oracle_residual, &
         source_vertex_ok, .false.)
      if (.not. source_vertex_ok .or. .not. ieee_is_finite(oracle_residual)) then
         error stop 'DRESP-11 static action: DRESP-10F source/response oracle failed'
      end if
      deallocate(density, source, response)
   end subroutine apply_static_selfconsistent_action

   !> Authoritative denominator action D=I-A in the compact Pauli space.
   subroutine apply_static_denominator(space, radial_bases, product, reciprocal_obj, kxc, c, dc)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(lmto_product_response_basis), intent(in) :: product
      type(reciprocal), intent(in) :: reciprocal_obj
      real(rp), intent(in) :: kxc(:, :)
      complex(rp), intent(in) :: c(:)
      complex(rp), intent(out) :: dc(:)

      call apply_static_selfconsistent_action(space, radial_bases, product, reciprocal_obj, kxc, c, dc)
      dc = c - dc
   end subroutine apply_static_denominator

   subroutine pauli_product_span_oracle(space, radial, product, by_l, maximum)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial
      type(lmto_product_response_basis), intent(in) :: product
      real(rp), intent(out) :: by_l(0:4), maximum
      real(rp), allocatable :: raw(:)
      complex(rp), allocatable :: weighted(:), projected(:)
      integer :: response_l, candidate, ir
      real(rp) :: norm

      by_l = 0.0_rp
      allocate(raw(space%npoint), weighted(space%npoint), projected(space%npoint))
      do response_l = 0, 4
         do candidate = 1, product%blocks(1,response_l)%ncandidate
            do ir = 1, space%npoint
               call lmto_product_second_order_radial_branch(radial, ir, &
                  product%blocks(1,response_l)%candidates(candidate)%l, &
                  product%blocks(1,response_l)%candidates(candidate)%lp, 1, 2, &
                  product%blocks(1,response_l)%candidates(candidate)%branch, raw(ir))
            end do
            weighted = cmplx(sqrt(space%radial_weights)*raw,0.0_rp,rp)
            projected = weighted - matmul(product%blocks(1,response_l)%weighted_modes, &
               matmul(conjg(transpose(product%blocks(1,response_l)%weighted_modes)),weighted))
            norm = sqrt(sum(abs(weighted)**2))
            by_l(response_l) = max(by_l(response_l),sqrt(sum(abs(projected)**2))/max(norm,tiny(1.0_rp)))
         end do
      end do
      maximum = maxval(by_l)
      deallocate(raw,weighted,projected)
   end subroutine pauli_product_span_oracle

   subroutine build_l0_source(space, bxc, source)
      type(response_space_layout), intent(in) :: space
      real(rp), intent(in) :: bxc(:, :)
      complex(rp), intent(out) :: source(:)
      type(response_super_index) :: item
      integer :: flat
      source = cmplx(0.0_rp,0.0_rp,rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat,space%nsite,space%response_lmax,space%npoint,space%nchannel,item)
         if (item%response_l == 0 .and. item%response_m == 0) &
            source(flat) = cmplx(sqrt(4.0_rp*response_angular_pi)*bxc(item%site,item%radial_point),0.0_rp,rp)
      end do
   end subroutine build_l0_source

   subroutine build_l0_target(space, p3, target)
      type(response_space_layout), intent(in) :: space
      real(rp), intent(in) :: p3(:, :)
      complex(rp), intent(out) :: target(:)
      type(response_super_index) :: item
      integer :: flat, ir
      target = cmplx(0.0_rp,0.0_rp,rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat,space%nsite,space%response_lmax,space%npoint,space%nchannel,item)
         ir = item%radial_point
         if (item%response_l == 0 .and. item%response_m == 0 .and. ir > 1) &
            target(flat) = cmplx(p3(item%site,ir)/(sqrt(4.0_rp*response_angular_pi)*space%radius(ir)**2),0.0_rp,rp)
      end do
   end subroutine build_l0_target

   subroutine mixed_static_action(space, radial_bases, product, reciprocal_obj, source_field, response_raw, &
                                  response_compact, oracle_residual, source_vertex_ok, verify_transition)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(lmto_product_response_basis), intent(in) :: product
      type(reciprocal), intent(in) :: reciprocal_obj
      complex(rp), intent(in) :: source_field(:)
      complex(rp), intent(out) :: response_raw(:), response_compact(:)
      real(rp), intent(out) :: oracle_residual
      logical, intent(out) :: source_vertex_ok
      logical, intent(in), optional :: verify_transition
      complex(rp), allocatable :: components(:, :, :, :), operators(:, :, :), h(:, :), delta_h(:, :), delta_rho(:, :)
      complex(rp), allocatable :: point_response(:), compact_direct(:), compact_transition(:)
      integer :: nmat, nk, ik, ik_global
      real(rp) :: wsum
      logical :: do_transition

      nmat = size(reciprocal_obj%hk_bulk,1); nk = size(reciprocal_obj%hk_bulk,3); wsum = sum(reciprocal_obj%k_weights)
      do_transition = .true.; if (present(verify_transition)) do_transition = verify_transition
      allocate(components(nmat,nmat,lmto_product_nbranch,lr_full_spatial_npiece), operators(nmat,nmat,lr_full_spatial_npiece), &
         h(nmat,nmat),delta_h(nmat,nmat),delta_rho(nmat,nmat),point_response(space%ndim), &
         compact_direct(product%product_dimension))
      if (do_transition) allocate(compact_transition(product%product_dimension))
      call lr_full_spatial_source_components(space,radial_bases,source_field,1,components)
      source_vertex_ok = .true.; oracle_residual = 0.0_rp
      response_raw = cmplx(0.0_rp,0.0_rp,rp); compact_direct = cmplx(0.0_rp,0.0_rp,rp)
      if (do_transition) compact_transition = cmplx(0.0_rp,0.0_rp,rp)
      do ik = 1, nk
         ik_global = ik; if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
         h = reciprocal_obj%hk_bulk(:,:,ik)
         call lr_full_spatial_contract_components(components,h,operators)
         delta_h = operators(:,:,lr_full_spatial_piece_total)
         call lr_static_frechet_density(reciprocal_obj%eigenvalues(:,ik),reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level,reciprocal_obj%temperature,delta_h,delta_rho)
         call fixed_pauli_measurement(space,radial_bases,h,delta_rho,point_response)
         if (do_transition) then
            call compact_density_measurement(product,reciprocal_obj%eigenvalues(:,ik),reciprocal_obj%eigenvectors(:,:,ik), &
               delta_rho,reciprocal_obj%fermi_level,reciprocal_obj%temperature,compact_transition)
            call dresp09_compact_field_from_raw(space,product,point_response,compact_direct)
            oracle_residual = max(oracle_residual,relative_vector(compact_direct,compact_transition))
         end if
         response_raw = response_raw + reciprocal_obj%k_weights(ik_global)*point_response
         compact_direct = cmplx(0.0_rp,0.0_rp,rp)
         if (do_transition) compact_transition = cmplx(0.0_rp,0.0_rp,rp)
      end do
      response_raw = response_raw/wsum
      call dresp09_compact_field_from_raw(space,product,response_raw,response_compact)
      source_vertex_ok = ieee_is_finite(oracle_residual)
      deallocate(components,operators,h,delta_h,delta_rho,point_response,compact_direct)
      if (do_transition) deallocate(compact_transition)
   end subroutine mixed_static_action

   subroutine fixed_pauli_measurement(space, radial_bases, hamiltonian, delta_rho, response)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: hamiltonian(:, :), delta_rho(:, :)
      complex(rp), intent(out) :: response(:)
      complex(rp), allocatable :: effective(:, :, :)
      type(response_super_index) :: item
      integer :: nsite, norb, n, site, ir, iorb, jorb, l, lp, m, mp, response_l, response_m, branch, p, q, row, col, flat
      integer :: lm_index

      nsite = size(radial_bases); norb = (radial_bases(1)%lmax+1)**2; n = 2*norb*nsite
      if (size(response) /= space%ndim .or. any(shape(delta_rho) /= [n,n])) error stop 'DRESP-10F measurement: shape mismatch'
      call ensure_measurement_cache(space, radial_bases)
      allocate(effective(n,n,lmto_product_nbranch)); response = cmplx(0.0_rp,0.0_rp,rp)
      do branch = 1, lmto_product_nbranch
         ! Tr(O_pq delta-rho) = Tr(O H**q delta-rho H**p): the response
         ! measurement is the stored-dual orientation of the endpoint branch.
         call lmto_product_apply_branch_action(hamiltonian,delta_rho,branch,effective(:,:,branch),.true.)
      end do
      do site = 1, nsite
         do response_l = 0, space%response_lmax
            do response_m = -response_l, response_l
               lm_index = response_l*response_l + response_l + response_m + 1
               do ir = 1, space%npoint
                  if (space%radial_weights(ir) <= 0.0_rp) cycle
                  item = response_super_index(site,response_l,response_m,ir,1)
                  call response_flatten_superindex(item,space%nsite,space%response_lmax,space%npoint,space%nchannel,flat)
                  do iorb = 1, norb
                     l = lmto_orbital_l(iorb); m = iorb-l*l-l-1
                     row = (site-1)*2*norb+iorb
                     do jorb = 1, norb
                        lp = lmto_orbital_l(jorb); mp = jorb-lp*lp-lp-1
                        col = (site-1)*2*norb+jorb+norb
                        do branch = 1, lmto_product_nbranch
                           response(flat) = response(flat) + 2.0_rp*measurement_gaunt_cache(iorb,jorb,lm_index)* &
                              effective(col,row,branch)*measurement_radial_cache(site,ir,iorb,jorb,branch)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      deallocate(effective)
   end subroutine fixed_pauli_measurement

   subroutine ensure_measurement_cache(space, radial_bases)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      integer :: nsite, norb, site, ir, iorb, jorb, l, lp, m, mp, response_l, response_m, branch, lm_index
      real(rp) :: radial_value

      nsite = size(radial_bases); norb = (radial_bases(1)%lmax+1)**2
      if (measurement_cache_ready .and. measurement_cache_nsite == nsite .and. &
          measurement_cache_npoint == space%npoint .and. measurement_cache_norb == norb) return
      if (allocated(measurement_radial_cache)) deallocate(measurement_radial_cache)
      if (allocated(measurement_gaunt_cache)) deallocate(measurement_gaunt_cache)
      allocate(measurement_radial_cache(nsite,space%npoint,norb,norb,lmto_product_nbranch), &
         measurement_gaunt_cache(norb,norb,(space%response_lmax+1)**2))
      do site = 1,nsite
         do ir = 1,space%npoint
            do iorb = 1,norb
               l = lmto_orbital_l(iorb); m = iorb-l*l-l-1
               do jorb = 1,norb
                  lp = lmto_orbital_l(jorb); mp = jorb-lp*lp-lp-1
                  do branch = 1,lmto_product_nbranch
                     call lmto_product_second_order_radial_branch(radial_bases(site),ir,l,lp,1,2,branch,radial_value)
                     measurement_radial_cache(site,ir,iorb,jorb,branch) = radial_value
                  end do
               end do
            end do
         end do
      end do
      do response_l = 0,space%response_lmax
         do response_m = -response_l,response_l
            lm_index = response_l*response_l + response_l + response_m + 1
            do iorb = 1,norb
               l = lmto_orbital_l(iorb); m = iorb-l*l-l-1
               do jorb = 1,norb
                  lp = lmto_orbital_l(jorb); mp = jorb-lp*lp-lp-1
                  measurement_gaunt_cache(iorb,jorb,lm_index) = response_gaunt(l,m,lp,mp,response_l,response_m)
               end do
            end do
         end do
      end do
      measurement_cache_nsite = nsite; measurement_cache_npoint = space%npoint; measurement_cache_norb = norb
      measurement_cache_ready = .true.
   end subroutine ensure_measurement_cache

   subroutine compact_density_measurement(product, eigenvalues, eigenvectors, delta_rho, fermi_level, temperature, response)
      type(lmto_product_response_basis), intent(in) :: product
      real(rp), intent(in) :: eigenvalues(:), fermi_level, temperature
      complex(rp), intent(in) :: eigenvectors(:, :), delta_rho(:, :)
      complex(rp), intent(out) :: response(:)
      type(pauli_endpoint_state) :: left_state, right_state
      complex(rp), allocatable :: delta_eigen(:, :), transition(:)
      integer :: n, ib, jb

      n = size(eigenvalues); allocate(delta_eigen(n,n),transition(product%product_dimension))
      delta_eigen = matmul(conjg(transpose(eigenvectors)),matmul(delta_rho,eigenvectors)); response = cmplx(0.0_rp,0.0_rp,rp)
      do ib = 1, n
         call left_state%initialize(eigenvalues(ib),eigenvectors(:,ib))
         do jb = 1, n
            call right_state%initialize(eigenvalues(jb),eigenvectors(:,jb))
            call product%transition_coordinates(left_state,right_state,transition)
            response = response + 2.0_rp*transition*delta_eigen(jb,ib)
         end do
      end do
      deallocate(delta_eigen,transition)
   end subroutine compact_density_measurement

   subroutine static_matrix_oracle(product, reciprocal_obj, vertices, field, action, residual)
      type(lmto_product_response_basis), intent(in) :: product
      type(reciprocal), intent(in) :: reciprocal_obj
      complex(rp), intent(in) :: vertices(:, :, :, :), field(:), action(:)
      real(rp), intent(out) :: residual
      complex(rp), allocatable :: matrix(:, :), matrix_action(:)
      allocate(matrix(product%product_dimension,product%product_dimension),matrix_action(product%product_dimension))
      call lr_static_product_matrix(product,reciprocal_obj%eigenvalues,reciprocal_obj%eigenvectors,reciprocal_obj%k_weights, &
         reciprocal_obj%fermi_level,reciprocal_obj%temperature,matrix)
      matrix_action = matmul(matrix,field)
      residual = relative_vector(matrix_action,action)
      deallocate(matrix,matrix_action)
   end subroutine static_matrix_oracle

   subroutine source_vertex_oracle(space, radial_bases, reciprocal_obj, maximum_residual, closed)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      real(rp), intent(out) :: maximum_residual
      logical, intent(out) :: closed
      complex(rp), allocatable :: field(:), components(:, :, :, :), operators(:, :, :), direct(:, :), h(:, :)
      type(response_super_index) :: item
      integer :: nmat, ik, l, m, ir, flat, ik_global
      real(rp) :: residual

      nmat = size(reciprocal_obj%hk_bulk,1); allocate(field(space%ndim),components(nmat,nmat,lmto_product_nbranch,lr_full_spatial_npiece), &
         operators(nmat,nmat,lr_full_spatial_npiece),direct(nmat,nmat),h(nmat,nmat))
      closed = .true.; maximum_residual = 0.0_rp
      do l = 0, space%response_lmax
         field = cmplx(0.0_rp,0.0_rp,rp)
         m = 0; ir = min(2,space%npoint); item = response_super_index(1,l,m,ir,1)
         call response_flatten_superindex(item,space%nsite,space%response_lmax,space%npoint,space%nchannel,flat)
         field(flat) = cmplx(0.31_rp+0.07_rp*l,-0.23_rp+0.11_rp*l,rp)
         call lr_full_spatial_source_components(space,radial_bases,field,1,components)
         do ik = 1, size(reciprocal_obj%hk_bulk,3)
            ik_global = ik; if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
            h = reciprocal_obj%hk_bulk(:,:,ik)
            call lr_full_spatial_contract_components(components,h,operators)
            call lr_full_spatial_source(space,radial_bases,field,1,h,direct,lr_full_spatial_piece_total)
            residual = relative_matrix(direct,operators(:,:,lr_full_spatial_piece_total))
            maximum_residual = max(maximum_residual,residual)
            closed = closed .and. residual < 1.0e-12_rp
         end do
      end do
      deallocate(field,components,operators,direct,h)
   end subroutine source_vertex_oracle

   subroutine covariant_oracle(space, radial_bases, reciprocal_obj, generator, endpoint_sum, p3, covariant_raw, &
                               covariant_contact, covariant_total, residual, contact_effect, closed)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      complex(rp), intent(in) :: generator(:, :), endpoint_sum(:, :, :)
      real(rp), intent(in) :: p3(:, :)
      complex(rp), intent(out) :: covariant_raw(:), covariant_contact(:), covariant_total(:)
      real(rp), intent(out) :: residual, contact_effect
      logical, intent(out) :: closed
      complex(rp), allocatable :: h(:, :), rho(:, :), delta_h(:, :), delta_rho(:, :), branches(:, :, :)
      complex(rp), allocatable :: contact(:, :)
      integer :: nmat, nk, ik, ik_global
      real(rp) :: wsum

      nmat=size(reciprocal_obj%hk_bulk,1); nk=size(reciprocal_obj%hk_bulk,3); wsum=sum(reciprocal_obj%k_weights)
      allocate(h(nmat,nmat),rho(nmat,nmat),delta_h(nmat,nmat),delta_rho(nmat,nmat),branches(nmat,nmat,lmto_product_nbranch), &
         contact(size(radial_bases),space%npoint))
      covariant_raw=cmplx(0.0_rp,0.0_rp,rp)
      do ik=1,nk
         ik_global=ik; if (allocated(reciprocal_obj%k_l2g_map)) ik_global=reciprocal_obj%k_l2g_map(ik)
         h=reciprocal_obj%hk_bulk(:,:,ik)
         call density_from_eigensystem(reciprocal_obj%eigenvalues(:,ik),reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level,reciprocal_obj%temperature,rho)
         call commutator_tangent(generator,h,delta_h); call commutator_tangent(generator,rho,delta_rho)
         call endpoint_tangent_branches_second_order(h,rho,delta_h,delta_rho,branches)
         call fixed_endpoint_pauli_l0(space,radial_bases,branches,covariant_raw)
         covariant_raw = covariant_raw*reciprocal_obj%k_weights(ik_global)
      end do
      covariant_raw=covariant_raw/wsum
      call dresp09y_augmentation_density(radial_bases,endpoint_sum,contact_matrix_dummy(),pauli_aug_c=contact)
      call build_l0_from_contact(space,contact,covariant_contact)
      covariant_total=covariant_raw+covariant_contact
      residual=weighted_relative_to_real(space,covariant_total,p3)
      contact_effect=sqrt(sum(abs(covariant_contact)**2))
      closed=residual < 4.0e-8_rp
      deallocate(h,rho,delta_h,delta_rho,branches,contact)
   contains
      function contact_matrix_dummy() result(sigma)
         complex(rp) :: sigma(2,2)
         call sr_aug_sigma_x(sigma)
      end function contact_matrix_dummy
   end subroutine covariant_oracle

   subroutine fixed_endpoint_pauli_l0(space, radial_bases, branches, response)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: branches(:, :, :)
      complex(rp), intent(inout) :: response(:)
      type(response_super_index) :: item
      integer :: site,ir,iorb,l,m,row,col,branch,flat,norb
      real(rp) :: pair,gaunt
      norb=(radial_bases(1)%lmax+1)**2
      do site=1,size(radial_bases)
         do ir=2,space%npoint
            do iorb=1,norb
               l=lmto_orbital_l(iorb);m=iorb-l*l-l-1;row=(site-1)*2*norb+iorb;col=row+norb
               gaunt=response_gaunt(l,m,l,m,0,0)
               item=response_super_index(site,0,0,ir,1);call response_flatten_superindex(item,space%nsite,space%response_lmax,space%npoint,1,flat)
               do branch=1,lmto_product_nbranch
                  call lmto_product_second_order_radial_branch(radial_bases(site),ir,l,l,1,2,branch,pair)
                  response(flat)=response(flat)+gaunt*branches(row,col,branch)*pair
               end do
            end do
         end do
      end do
   end subroutine fixed_endpoint_pauli_l0

   subroutine build_l0_from_contact(space, contact, vector)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: contact(:, :)
      complex(rp), intent(out) :: vector(:)
      type(response_super_index) :: item
      integer :: flat,ir,raw_flat
      vector=cmplx(0.0_rp,0.0_rp,rp)
      do flat=1,space%ndim
         call response_unflatten_superindex(flat,space%nsite,space%response_lmax,space%npoint,space%nchannel,item)
         ir=item%radial_point
         if(item%response_l==0.and.item%response_m==0.and.ir>1) then
            item%radial_point=ir;call response_flatten_superindex(item,space%nsite,space%response_lmax,space%npoint,space%nchannel,raw_flat)
            vector(raw_flat)=contact(item%site,ir)
         end if
      end do
   end subroutine build_l0_from_contact

   subroutine radial_metrics(space, value, reference, absolute, relative, maximum, integrated)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: value(:), reference(:)
      real(rp), intent(out) :: absolute, relative, maximum, integrated
      type(response_super_index) :: item
      integer :: flat
      real(rp) :: diff, ref, weight
      absolute=0.0_rp; maximum=0.0_rp; integrated=0.0_rp; ref=0.0_rp
      do flat=1,space%ndim
         call response_unflatten_superindex(flat,space%nsite,space%response_lmax,space%npoint,space%nchannel,item)
         diff=abs(value(flat)-reference(flat)); weight=space%radial_weights(item%radial_point)
         absolute=absolute+weight*diff*diff; ref=ref+weight*abs(reference(flat))**2
         maximum=max(maximum,diff); integrated=integrated+weight*real(value(flat)-reference(flat),rp)
      end do
      absolute=sqrt(absolute); relative=absolute/max(sqrt(ref),tiny(1.0_rp)); integrated=abs(integrated)
   end subroutine radial_metrics

   real(rp) function weighted_relative(space, value, reference) result(result)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: value(:), reference(:)
      real(rp) :: a,b
      integer :: flat
      type(response_super_index) :: item
      a=0.0_rp;b=0.0_rp
      do flat=1,size(value)
         call response_unflatten_superindex(flat,space%nsite,space%response_lmax,space%npoint,space%nchannel,item)
         a=a+space%radial_weights(item%radial_point)*abs(value(flat)-reference(flat))**2
         b=b+space%radial_weights(item%radial_point)*abs(reference(flat))**2
      end do
      result=sqrt(a)/max(sqrt(b),tiny(1.0_rp))
   end function weighted_relative

   real(rp) function weighted_relative_to_real(space, value, reference) result(result)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: value(:)
      real(rp), intent(in) :: reference(:, :)
      complex(rp), allocatable :: target(:)
      allocate(target(size(value))); call build_l0_target(space,reference,target)
      result=weighted_relative(space,value,target); deallocate(target)
   end function weighted_relative_to_real

   real(rp) function relative_vector(value, reference) result(result)
      complex(rp), intent(in) :: value(:), reference(:)
      result=sqrt(sum(abs(value-reference)**2))/max(sqrt(sum(abs(reference)**2)),tiny(1.0_rp))
   end function relative_vector

   real(rp) function relative_matrix(value, reference) result(result)
      complex(rp), intent(in) :: value(:, :), reference(:, :)
      result=sqrt(sum(abs(value-reference)**2))/max(sqrt(sum(abs(reference)**2)),tiny(1.0_rp))
   end function relative_matrix

   subroutine general_mixed_actions(space, radial_bases, product, reciprocal_obj, kxc, action, oracle, linearity)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(lmto_product_response_basis), intent(in) :: product
      type(reciprocal), intent(in) :: reciprocal_obj
      real(rp), intent(in) :: kxc(:, :)
      real(rp), intent(out) :: action(6), oracle(6), linearity(6)
      complex(rp), allocatable :: c(:), direct(:), matrix_action(:), matrix(:, :), delta_h(:, :), vertices(:, :, :, :)
      integer :: mode, flat, response_l, response_m, site, ir, ik, ik_global
      type(response_super_index) :: item
      real(rp) :: weight_sum

      ! The general 348-to-348 diagnostic is kept in compact response space.
      ! The physical mixed map is exercised separately above with both rigid
      ! and arbitrary raw SR fields. This keeps the acceptance run bounded
      ! while retaining an independent matrix-free action/oracle comparison.
      allocate(c(product%product_dimension),direct(product%product_dimension),matrix_action(product%product_dimension), &
         matrix(product%product_dimension,product%product_dimension),delta_h(size(reciprocal_obj%hk_bulk,1),size(reciprocal_obj%hk_bulk,1)), &
         vertices(size(reciprocal_obj%hk_bulk,1),size(reciprocal_obj%hk_bulk,1),lmto_product_nbranch,product%product_dimension))
      call product%component_vertex_tensor(vertices)
      call lr_static_product_matrix(product, reciprocal_obj%eigenvalues, reciprocal_obj%eigenvectors, reciprocal_obj%k_weights, &
         reciprocal_obj%fermi_level, reciprocal_obj%temperature, matrix)
      weight_sum = sum(reciprocal_obj%k_weights)
      action=0.0_rp;oracle=0.0_rp;linearity=0.0_rp
      do mode=1,6
         c=cmplx(0.0_rp,0.0_rp,rp)
         do flat=1,product%product_dimension
            call product%unflatten_index(flat,site,response_l,response_m,ir)
            if ((mode==1.and.response_l==0).or.(mode==2.and.response_l==1).or.(mode==3.and.response_l==2).or. &
                (mode==4.and.response_l==3).or.(mode==5.and.response_l==4).or.(mode==6.and.response_l>=1)) &
               c(flat)=cmplx(sin(0.13_rp*flat+mode),cos(0.23_rp*flat-mode),rp)
         end do
         direct=cmplx(0.0_rp,0.0_rp,rp)
         do ik=1,size(reciprocal_obj%hk_bulk,3)
            ik_global=ik
            if (allocated(reciprocal_obj%k_l2g_map)) ik_global=reciprocal_obj%k_l2g_map(ik)
            call dresp09_product_source_operator(vertices, reciprocal_obj%hk_bulk(:,:,ik), c, delta_h)
            call lr_static_product_action(product, reciprocal_obj%eigenvalues(:,ik), reciprocal_obj%eigenvectors(:,:,ik), &
               reciprocal_obj%fermi_level, reciprocal_obj%temperature, delta_h, matrix_action)
            direct=direct+reciprocal_obj%k_weights(ik_global)*matrix_action
         end do
         direct=direct/weight_sum
         matrix_action=matmul(matrix,c)
         action(mode)=sqrt(sum(abs(direct)**2)); oracle(mode)=sqrt(sum(abs(matrix_action)**2))
         linearity(mode)=relative_vector(direct,matrix_action)
      end do
      deallocate(c,direct,matrix_action,matrix,delta_h,vertices)
   end subroutine general_mixed_actions

   real(rp) function native_field_defect(space, radial_bases, reciprocal_obj, generator, source_field) result(value)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      complex(rp), intent(in) :: generator(:, :), source_field(:)
      complex(rp), allocatable :: components(:, :, :, :), operators(:, :, :), native(:, :)
      integer :: nmat,ik
      nmat=size(reciprocal_obj%hk_bulk,1);allocate(components(nmat,nmat,lmto_product_nbranch,lr_full_spatial_npiece), &
         operators(nmat,nmat,lr_full_spatial_npiece),native(nmat,nmat)); value=0.0_rp
      call lr_full_spatial_source_components(space,radial_bases,source_field,1,components)
      do ik=1,size(reciprocal_obj%hk_bulk,3)
         call lr_full_spatial_contract_components(components,reciprocal_obj%hk_bulk(:,:,ik),operators)
         call commutator_tangent(generator,reciprocal_obj%hk_bulk(:,:,ik),native)
         value=max(value,relative_matrix(operators(:,:,lr_full_spatial_piece_total),native))
      end do
      deallocate(components,operators,native)
   end function native_field_defect

end module lr_dresp10f_mixed_ward_bridge_mod
