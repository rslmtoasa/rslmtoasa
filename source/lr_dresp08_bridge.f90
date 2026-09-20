!------------------------------------------------------------------------------
! DRESP-08 accepted-state material audit.
!
! This bridge is diagnostic-only.  It reconstructs the production second-order
! Hamiltonian from the stored real-space blocks, rebuilds the native transverse
! tangent from structure constants and transformed potential parameters, and
! contracts that tangent with the already certified product-response oracle.
!------------------------------------------------------------------------------
module lr_dresp08_bridge_mod

   use precision_mod, only: rp
   use mpi_mod, only: rank
   use basis_mod, only: nb
   use reciprocal_mod, only: reciprocal
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use radial_ground_state_mod, only: radial_ground_state, radial_simpson_weight
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use lr_response_space_mod, only: response_space_layout
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, lmto_product_channel_plus
   use lr_exact_ks_ward_mod, only: exact_ks_spin_rotation_generator, exact_ks_product_rotation_sweep, &
      exact_ks_product_field_response, exact_ks_extract_collinear_field, exact_ks_build_transverse_field
   use lr_dresp07_radial_oracle_mod, only: dresp07_project_local_spherical_operator
   use lr_dresp08_native_mapping_mod, only: dresp08_build_h2, dresp08_build_product_tangent, &
      dresp08_commutator_tangent, dresp08_build_native_realspace_tangent, dresp08_build_native_onsite_tangents, &
      dresp08_build_spinor_operator, dresp08_block_diagonal, dresp08_extract_spin_blocks, dresp08_relative_residual
   implicit none
   private

   public :: run_dresp08_native_second_order

contains

   subroutine run_dresp08_native_second_order(output_file, response_space, radial_bases, ground_states, reciprocal_obj, &
                                              lattice_obj, hamiltonian_obj, eta_values)
      character(len=*), intent(in) :: output_file
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj
      real(rp), intent(in) :: eta_values(:)

      type(lmto_product_response_basis) :: product_plus
      complex(rp), allocatable :: h(:, :, :), h2(:, :, :), h2_recon(:, :, :), o(:, :, :), enu(:, :, :)
      complex(rp), allocatable :: eeo_k(:, :, :), d_h(:, :, :), d_h_comm(:, :, :), d_h2(:, :, :), d_h2_comm(:, :, :)
      complex(rp), allocatable :: d_o(:, :, :), d_enu(:, :, :), d_o_comm(:, :, :), d_enu_comm(:, :, :)
      complex(rp), allocatable :: term_e(:, :, :), term_h(:, :, :)
      complex(rp), allocatable :: term_left(:, :, :), term_mid(:, :, :), term_right(:, :, :)
      complex(rp), allocatable :: b_e(:, :, :), b_h(:, :, :), b_hoh(:, :, :), b_h_native(:, :, :)
      complex(rp), allocatable :: h_up(:, :), h_down(:, :), o_up(:, :), o_down(:, :), e_up(:, :), e_down(:, :)
      complex(rp), allocatable :: h2_up(:, :), h2_down(:, :), b_gamma(:, :), best_gamma(:, :), best_coefficients(:, :)
      complex(rp), allocatable :: onsite_rebuilt(:, :), onsite_rebuilt_enu(:, :)
      complex(rp), allocatable :: generator(:, :), b_direct(:, :), dfield_exact(:, :)
      complex(rp), allocatable :: response_rot(:), response_spec(:), response_retarded(:), response_native(:)
      complex(rp), allocatable :: response_spec_reference(:), response_native_reference(:)
      complex(rp), allocatable :: response_native_retarded(:), response_radial(:), response_radial_retarded(:)
      complex(rp), allocatable :: one_static(:), one_retarded(:), bxc_orbital(:, :), bks_orbital(:, :)
      complex(rp), allocatable :: d_ee(:, :, :, :), d_obarm(:, :, :), d_enim(:, :, :), o_types(:, :, :), e_types(:, :, :)
      real(rp), allocatable :: radial_values(:, :)
      integer, allocatable :: site_types(:)
      real(rp) :: max_h2_element, max_h2_frobenius, max_h2_relative, max_h2_hermiticity
      real(rp) :: max_spin_field, max_spin_field_relative, max_eeo_identity
      real(rp) :: max_parameter_average, max_obarm_provenance, max_enim_provenance
      real(rp) :: first_tangent_relative, overlap_tangent_relative, enu_tangent_relative
      real(rp) :: product_tangent_relative, product_tangent_max
      real(rp) :: norm_e, norm_h, norm_hoh, norm_bh, norm_native, weight_sum
      real(rp) :: gamma_projection_relative, gamma_offdiag, gamma_anisotropy, gamma_cross_l, gamma_intersite
      real(rp) :: rk_be, rk_bh, rk_bhoh, rk_bh_total, rk_dh, rk_dh2
      real(rp) :: spin_residual, radial_bxc_relative, radial_bks_relative, radial_bxc_action, radial_bks_action
      real(rp) :: radial_to_e, radial_to_eh, native_response_relative, native_static_relative
      real(rp) :: gamma_component_projection(4), gamma_component_offdiag(4), gamma_component_anisotropy(4)
      real(rp) :: gamma_component_cross_l(4), gamma_component_intersite(4)
      real(rp), allocatable :: eta_exact(:), eta_native(:), eta_native_exact(:)
      real(rp), allocatable :: tangent_enu_norm(:), tangent_h_norm(:), tangent_dh_o_h_norm(:)
      real(rp), allocatable :: tangent_h_do_h_norm(:), tangent_h_o_dh_norm(:), tangent_exact_norm(:)
      real(rp), allocatable :: tangent_comm_residual(:)
      real(rp), allocatable :: gamma_bh_k(:), gamma_bh_offdiag_k(:), gamma_bh_anisotropy_k(:)
      real(rp), allocatable :: gamma_bh_cross_l_k(:), gamma_bh_intersite_k(:)
      real(rp) :: max_finite_unused, direct_onsite_norm, direct_nonlocal_norm, direct_spin_offdiag
      real(rp) :: exact_ep, exact_dr, exact_rel
      integer :: nsite, norb_site, norb, nmat, nk, ik, neta, isite, unit, ios, itype
      logical :: pass_h2, pass_spin, pass_tangent, pass_response
      character(len=64) :: classification

      nsite = size(ground_states)
      if (nsite < 1 .or. size(radial_bases) /= nsite .or. response_space%nsite /= nsite) then
         error stop 'DRESP-08: radial/response site dimensions are inconsistent'
      end if
      if (.not. allocated(reciprocal_obj%hk_bulk) .or. .not. allocated(reciprocal_obj%k_weights)) then
         error stop 'DRESP-08: accepted reciprocal Hamiltonian is incomplete'
      end if
      if (.not. hamiltonian_obj%hoh) error stop 'DRESP-08: the accepted Hamiltonian is not on the HOH path'
      if (hamiltonian_obj%ccor_2c .or. hamiltonian_obj%orb_pol .or. hamiltonian_obj%hubbard_u_general_check .or. &
          hamiltonian_obj%hubbard_v_check .or. hamiltonian_obj%local_axis) then
         error stop 'DRESP-08: excluded additive Hamiltonian correction is enabled'
      end if
      if (lattice_obj%nrec /= nsite) error stop 'DRESP-08: lattice/response site counts differ'

      norb_site = (radial_bases(1)%lmax + 1)**2
      norb = norb_site*nsite
      nmat = size(reciprocal_obj%hk_bulk, 1)
      nk = size(reciprocal_obj%hk_bulk, 3)
      if (nmat /= nb*nsite .or. nb /= 2*norb_site) error stop 'DRESP-08: accepted basis dimensions differ'
      if (size(reciprocal_obj%k_weights) /= nk .or. sum(reciprocal_obj%k_weights) <= 0.0_rp) then
         error stop 'DRESP-08: accepted k-weight array is invalid'
      end if
      weight_sum = sum(reciprocal_obj%k_weights)

      call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      allocate(site_types(nsite))
      site_types = lattice_obj%ib(1:nsite)
      allocate(o_types(nb, nb, size(hamiltonian_obj%obarm, 3)), e_types(nb, nb, size(hamiltonian_obj%enim, 3)))
      o_types = hamiltonian_obj%obarm
      e_types = hamiltonian_obj%enim
      allocate(d_ee, source=hamiltonian_obj%ee)
      allocate(d_obarm, source=hamiltonian_obj%obarm)
      allocate(d_enim, source=hamiltonian_obj%enim)
      call dresp08_build_native_realspace_tangent(hamiltonian_obj, d_ee)
      call dresp08_build_native_onsite_tangents(hamiltonian_obj, d_obarm, d_enim)
      allocate(h(nmat,nmat,nk), h2(nmat,nmat,nk), h2_recon(nmat,nmat,nk), o(nmat,nmat,nk), enu(nmat,nmat,nk), &
         eeo_k(nmat,nmat,nk), d_h(nmat,nmat,nk), d_h_comm(nmat,nmat,nk), d_h2(nmat,nmat,nk), &
         d_h2_comm(nmat,nmat,nk), d_o(nmat,nmat,nk), d_enu(nmat,nmat,nk), d_o_comm(nmat,nmat,nk), &
         d_enu_comm(nmat,nmat,nk), term_e(nmat,nmat,nk), &
         term_h(nmat,nmat,nk), term_left(nmat,nmat,nk), term_mid(nmat,nmat,nk), term_right(nmat,nmat,nk), &
         b_e(norb,norb,nk), b_h(norb,norb,nk), b_hoh(norb,norb,nk), b_h_native(norb,norb,nk))
      allocate(generator(nmat,nmat), h_up(norb,norb), h_down(norb,norb), o_up(norb,norb), o_down(norb,norb), &
         e_up(norb,norb), e_down(norb,norb), h2_up(norb,norb), h2_down(norb,norb), b_gamma(norb,norb), &
         best_gamma(norb,norb), best_coefficients(nsite, radial_bases(1)%lmax+1), b_direct(norb,norb), &
         dfield_exact(nmat,nmat), onsite_rebuilt(nb,nb), onsite_rebuilt_enu(nb,nb))
      call exact_ks_spin_rotation_generator(norb_site, nsite, generator)

      max_parameter_average = 0.0_rp
      max_obarm_provenance = 0.0_rp
      max_enim_provenance = 0.0_rp
      do itype = 1, size(o_types, 3)
         max_parameter_average = max(max_parameter_average, maxval(abs(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%cx0 - &
            0.5_rp*(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%cx(:,1) + &
                     hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%cx(:,2)))))
         max_parameter_average = max(max_parameter_average, maxval(abs(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%cx1 - &
            0.5_rp*(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%cx(:,1) - &
                     hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%cx(:,2)))))
         max_parameter_average = max(max_parameter_average, maxval(abs(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%wx0 - &
            0.5_rp*(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%wx(:,1) + &
                     hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%wx(:,2)))))
         max_parameter_average = max(max_parameter_average, maxval(abs(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%wx1 - &
            0.5_rp*(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%wx(:,1) - &
                     hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%wx(:,2)))))
         max_parameter_average = max(max_parameter_average, maxval(abs(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%cex0 - &
            0.5_rp*(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%cex(:,1) + &
                     hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%cex(:,2)))))
         max_parameter_average = max(max_parameter_average, maxval(abs(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%cex1 - &
            0.5_rp*(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%cex(:,1) - &
                     hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%cex(:,2)))))
         max_parameter_average = max(max_parameter_average, maxval(abs(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%obx0 - &
            0.5_rp*(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%obx(:,1) + &
                     hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%obx(:,2)))))
         max_parameter_average = max(max_parameter_average, maxval(abs(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%obx1 - &
            0.5_rp*(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%obx(:,1) - &
                     hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%obx(:,2)))))
         call dresp08_build_spinor_operator(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%obx0(1:norb), &
            hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%obx1(1:norb), &
            real(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%mom(:), rp), onsite_rebuilt)
         call dresp08_build_spinor_operator(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%cx0(1:norb) - &
            hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%cex0(1:norb), &
            hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%cx1(1:norb) - &
            hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%cex1(1:norb), &
            real(hamiltonian_obj%charge%lattice%symbolic_atoms(itype)%potential%mom(:), rp), onsite_rebuilt_enu)
         max_obarm_provenance = max(max_obarm_provenance, dresp08_relative_residual(onsite_rebuilt, o_types(:,:,itype)))
         max_enim_provenance = max(max_enim_provenance, dresp08_relative_residual(onsite_rebuilt_enu, e_types(:,:,itype)))
      end do

      max_h2_element = 0.0_rp; max_h2_frobenius = 0.0_rp; max_h2_relative = 0.0_rp
      max_h2_hermiticity = 0.0_rp; max_spin_field = 0.0_rp; max_spin_field_relative = 0.0_rp
      max_eeo_identity = 0.0_rp; first_tangent_relative = 0.0_rp; overlap_tangent_relative = 0.0_rp
      enu_tangent_relative = 0.0_rp; product_tangent_relative = 0.0_rp; product_tangent_max = 0.0_rp
      allocate(tangent_enu_norm(nk), tangent_h_norm(nk), tangent_dh_o_h_norm(nk), tangent_h_do_h_norm(nk), &
         tangent_h_o_dh_norm(nk), tangent_exact_norm(nk), tangent_comm_residual(nk), gamma_bh_k(nk), &
         gamma_bh_offdiag_k(nk), gamma_bh_anisotropy_k(nk), gamma_bh_cross_l_k(nk), gamma_bh_intersite_k(nk))
      do ik = 1, nk
         call reciprocal_kpoint(reciprocal_obj, ik, hamiltonian_obj%ee, h(:,:,ik))
         call reciprocal_kpoint(reciprocal_obj, ik, hamiltonian_obj%eeo, eeo_k(:,:,ik))
         call reciprocal_kpoint(reciprocal_obj, ik, d_ee, d_h(:,:,ik))
         call dresp08_block_diagonal(o_types, site_types, o(:,:,ik))
         call dresp08_block_diagonal(e_types, site_types, enu(:,:,ik))
         call dresp08_block_diagonal(d_obarm, site_types, d_o(:,:,ik))
         call dresp08_block_diagonal(d_enim, site_types, d_enu(:,:,ik))
         call dresp08_build_h2(h(:,:,ik), o(:,:,ik), enu(:,:,ik), h2_recon(:,:,ik))
         h2(:,:,ik) = reciprocal_obj%hk_bulk(:,:,ik)
         max_h2_element = max(max_h2_element, maxval(abs(h2_recon(:,:,ik)-h2(:,:,ik))))
         max_h2_frobenius = max(max_h2_frobenius, sqrt(sum(abs(h2_recon(:,:,ik)-h2(:,:,ik))**2)))
         max_h2_relative = max(max_h2_relative, dresp08_relative_residual(h2_recon(:,:,ik), h2(:,:,ik)))
         max_h2_hermiticity = max(max_h2_hermiticity, maxval(abs(h2_recon(:,:,ik)-transpose(conjg(h2_recon(:,:,ik))))))
         max_eeo_identity = max(max_eeo_identity, dresp08_relative_residual(eeo_k(:,:,ik), matmul(h(:,:,ik),o(:,:,ik))))

         call dresp08_commutator_tangent(generator, h(:,:,ik), d_h_comm(:,:,ik))
         call dresp08_commutator_tangent(generator, o(:,:,ik), d_o_comm(:,:,ik))
         call dresp08_commutator_tangent(generator, enu(:,:,ik), d_enu_comm(:,:,ik))
         call dresp08_build_product_tangent(d_enu(:,:,ik), d_h(:,:,ik), d_o(:,:,ik), h(:,:,ik), o(:,:,ik), &
            d_h2(:,:,ik), term_e(:,:,ik), term_h(:,:,ik), term_left(:,:,ik), term_mid(:,:,ik), term_right(:,:,ik))
         call dresp08_commutator_tangent(generator, h2(:,:,ik), d_h2_comm(:,:,ik))
         first_tangent_relative = max(first_tangent_relative, dresp08_relative_residual(d_h(:,:,ik), d_h_comm(:,:,ik)))
         overlap_tangent_relative = max(overlap_tangent_relative, dresp08_relative_residual(d_o(:,:,ik), &
            d_o_comm(:,:,ik)))
         enu_tangent_relative = max(enu_tangent_relative, dresp08_relative_residual(d_enu(:,:,ik), d_enu_comm(:,:,ik)))
         product_tangent_relative = max(product_tangent_relative, dresp08_relative_residual(d_h2(:,:,ik), d_h2_comm(:,:,ik)))
         product_tangent_max = max(product_tangent_max, maxval(abs(d_h2(:,:,ik)-d_h2_comm(:,:,ik))))
         tangent_enu_norm(ik) = sqrt(sum(abs(term_e(:,:,ik))**2))
         tangent_h_norm(ik) = sqrt(sum(abs(term_h(:,:,ik))**2))
         tangent_dh_o_h_norm(ik) = sqrt(sum(abs(term_left(:,:,ik))**2))
         tangent_h_do_h_norm(ik) = sqrt(sum(abs(term_mid(:,:,ik))**2))
         tangent_h_o_dh_norm(ik) = sqrt(sum(abs(term_right(:,:,ik))**2))
         tangent_exact_norm(ik) = sqrt(sum(abs(d_h2_comm(:,:,ik))**2))
         tangent_comm_residual(ik) = dresp08_relative_residual(d_h2(:,:,ik), d_h2_comm(:,:,ik))

         call dresp08_extract_spin_blocks(h(:,:,ik), norb_site, nsite, h_up, h_down)
         call dresp08_extract_spin_blocks(o(:,:,ik), norb_site, nsite, o_up, o_down)
         call dresp08_extract_spin_blocks(enu(:,:,ik), norb_site, nsite, e_up, e_down)
         call dresp08_build_h2(h_up, o_up, e_up, h2_up)
         call dresp08_build_h2(h_down, o_down, e_down, h2_down)
         b_e(:,:,ik) = 0.5_rp*(e_up-e_down)
         b_h(:,:,ik) = 0.5_rp*(h_up-h_down)
         b_hoh(:,:,ik) = b_h(:,:,ik) + b_e(:,:,ik) - 0.5_rp*(h2_up-h2_down)
         b_h_native(:,:,ik) = 0.5_rp*(h2_up-h2_down)
         call exact_ks_extract_collinear_field(reciprocal_obj%hk_bulk(:,:,ik), norb_site, nsite, b_direct, dfield_exact, &
            direct_onsite_norm, direct_nonlocal_norm, direct_spin_offdiag)
         max_spin_field = max(max_spin_field, maxval(abs(b_h_native(:,:,ik)-b_direct)))
         max_spin_field_relative = max(max_spin_field_relative, dresp08_relative_residual(b_h_native(:,:,ik), b_direct))
         call dresp07_project_local_spherical_operator(b_h_native(:,:,ik), norb_site, nsite, radial_bases(1)%lmax, &
            best_gamma, best_coefficients, max_finite_unused, gamma_bh_k(ik), gamma_bh_offdiag_k(ik), &
            gamma_bh_anisotropy_k(ik), gamma_bh_cross_l_k(ik), gamma_bh_intersite_k(ik))
      end do

      pass_h2 = max_h2_relative < 2.0e-11_rp .and. max_eeo_identity < 2.0e-11_rp .and. &
         max_h2_hermiticity < 2.0e-10_rp .and. &
         max_parameter_average < 2.0e-11_rp .and. max_obarm_provenance < 2.0e-11_rp .and. &
         max_enim_provenance < 2.0e-11_rp
      pass_spin = max_spin_field_relative < 2.0e-11_rp
      pass_tangent = first_tangent_relative < 2.0e-10_rp .and. overlap_tangent_relative < 2.0e-10_rp .and. &
         enu_tangent_relative < 2.0e-10_rp .and. product_tangent_relative < 2.0e-10_rp

      call weighted_field_metrics(b_e, reciprocal_obj%k_weights, norm_e, rk_be)
      call weighted_field_metrics(b_h, reciprocal_obj%k_weights, norm_h, rk_bh)
      call weighted_field_metrics(b_hoh, reciprocal_obj%k_weights, norm_hoh, rk_bhoh)
      call weighted_field_metrics(b_h_native, reciprocal_obj%k_weights, norm_bh, rk_bh_total)
      call weighted_field_metrics(d_h, reciprocal_obj%k_weights, max_finite_unused, rk_dh)
      call weighted_field_metrics(d_h2, reciprocal_obj%k_weights, max_finite_unused, rk_dh2)

      b_gamma = b_h_native(:,:,1)
      call dresp07_project_local_spherical_operator(b_e(:,:,1), norb_site, nsite, radial_bases(1)%lmax, best_gamma, &
         best_coefficients, max_finite_unused, gamma_component_projection(1), gamma_component_offdiag(1), &
         gamma_component_anisotropy(1), gamma_component_cross_l(1), gamma_component_intersite(1))
      call dresp07_project_local_spherical_operator(b_h(:,:,1), norb_site, nsite, radial_bases(1)%lmax, best_gamma, &
         best_coefficients, max_finite_unused, gamma_component_projection(2), gamma_component_offdiag(2), &
         gamma_component_anisotropy(2), gamma_component_cross_l(2), gamma_component_intersite(2))
      call dresp07_project_local_spherical_operator(b_hoh(:,:,1), norb_site, nsite, radial_bases(1)%lmax, best_gamma, &
         best_coefficients, max_finite_unused, gamma_component_projection(3), gamma_component_offdiag(3), &
         gamma_component_anisotropy(3), gamma_component_cross_l(3), gamma_component_intersite(3))
      call dresp07_project_local_spherical_operator(b_gamma, norb_site, nsite, radial_bases(1)%lmax, best_gamma, &
         best_coefficients, max_finite_unused, gamma_component_projection(4), gamma_component_offdiag(4), &
         gamma_component_anisotropy(4), gamma_component_cross_l(4), gamma_component_intersite(4))
      gamma_projection_relative = gamma_component_projection(4)
      gamma_offdiag = gamma_component_offdiag(4)
      gamma_anisotropy = gamma_component_anisotropy(4)
      gamma_cross_l = gamma_component_cross_l(4)
      gamma_intersite = gamma_component_intersite(4)

      allocate(radial_values(nsite, response_space%npoint), bxc_orbital(norb,norb), bks_orbital(norb,norb))
      call collect_radial_values(ground_states, radial_values, .true.)
      call map_radial_values(radial_bases, radial_values, bxc_orbital)
      call collect_radial_values(ground_states, radial_values, .false.)
      call map_radial_values(radial_bases, radial_values, bks_orbital)
      call field_matrix_residual(b_gamma, bxc_orbital, radial_bxc_relative, radial_bxc_action)
      call field_matrix_residual(b_gamma, bks_orbital, radial_bks_relative, radial_bks_action)
      call field_matrix_residual(b_gamma, diagonal_field_from_e(b_e(:,:,1)), radial_to_e, max_finite_unused)
      call field_matrix_residual(b_gamma, diagonal_field_from_e(b_e(:,:,1)+b_h(:,:,1)), radial_to_eh, max_finite_unused)

      allocate(response_rot(product_plus%product_dimension), response_spec(product_plus%product_dimension), &
         response_retarded(product_plus%product_dimension), response_native(product_plus%product_dimension), &
         response_native_retarded(product_plus%product_dimension), response_radial(product_plus%product_dimension), &
         response_radial_retarded(product_plus%product_dimension), one_static(product_plus%product_dimension), &
         one_retarded(product_plus%product_dimension))
      call exact_ks_product_rotation_sweep(product_plus, reciprocal_obj%hk_bulk, reciprocal_obj%eigenvalues, &
         reciprocal_obj%eigenvectors, reciprocal_obj%k_weights, reciprocal_obj%fermi_level, reciprocal_obj%temperature, &
         eta_values(1), response_rot, response_spec, response_retarded, exact_ep, exact_dr, exact_rel, &
         .false., 0)
      call evaluate_native_response(product_plus, reciprocal_obj, d_h2, eta_values(1), response_native, response_native_retarded, .true.)
      call evaluate_orbital_response(product_plus, reciprocal_obj, bxc_orbital, eta_values(1), response_radial, response_radial_retarded)
      allocate(response_spec_reference(product_plus%product_dimension), response_native_reference(product_plus%product_dimension))
      response_spec_reference = response_spec
      response_native_reference = response_native
      native_response_relative = vector_relative(response_native, response_spec)
      native_static_relative = vector_relative(response_native, response_rot)

      neta = size(eta_values)
      allocate(eta_exact(neta), eta_native(neta), eta_native_exact(neta))
      do ik = 1, neta
         call exact_ks_product_rotation_sweep(product_plus, reciprocal_obj%hk_bulk, reciprocal_obj%eigenvalues, &
            reciprocal_obj%eigenvectors, reciprocal_obj%k_weights, reciprocal_obj%fermi_level, reciprocal_obj%temperature, &
            eta_values(ik), response_rot, response_spec, response_retarded, exact_ep, exact_dr, exact_rel, &
            .true., 0)
         call evaluate_native_response(product_plus, reciprocal_obj, d_h2, eta_values(ik), response_native, response_native_retarded, .false.)
         eta_exact(ik) = vector_relative(response_retarded, response_spec_reference)
         eta_native(ik) = vector_relative(response_native_retarded, response_native_reference)
         eta_native_exact(ik) = vector_relative(response_native_retarded, response_retarded)
      end do
      pass_response = native_response_relative < 2.0e-10_rp .and. native_static_relative < 2.0e-10_rp .and. &
         maxval(eta_native_exact) < 2.0e-10_rp
      if (pass_h2 .and. pass_spin .and. pass_tangent .and. pass_response) then
         classification = 'NATIVE_SECOND_ORDER_MAPPING_CLOSED'
      else if (.not. pass_h2) then
         classification = 'BLOCKED - SECOND-ORDER HAMILTONIAN RECONSTRUCTION'
      else if (.not. pass_tangent) then
         classification = 'MIXED - NATIVE TANGENT FAILURE'
      else if (.not. pass_response) then
         classification = 'MIXED - PRODUCT RESPONSE FAILURE'
      else
         classification = 'MIXED'
      end if

      if (rank == 0) then
         open(newunit=unit, file=trim(output_file), status='replace', action='write', iostat=ios)
         if (ios /= 0) error stop 'DRESP-08: cannot open diagnostic artifact'
         write(unit,'(a)') '# DRESP-08 native second-order LMTO field mapping'
         write(unit,'(a)') '# BES_Halle = OFF'
         write(unit,'(a)') '# Kxc_tuning = OFF'
         write(unit,'(a)') '# Goldstone_correction = OFF'
         write(unit,'(a)') '# SCF = production accepted-state handoff; DRESP-08 bridge performs no SCF'
         write(unit,'(a)') '# H2_contract = H2 = E_nu + h - h o h'
         write(unit,'(a,i0)') 'product_dimension = ', product_plus%product_dimension
         write(unit,'(a,es24.16)') 'H2_max_element_residual = ', max_h2_element
         write(unit,'(a,es24.16)') 'H2_max_frobenius_residual = ', max_h2_frobenius
         write(unit,'(a,es24.16)') 'H2_max_relative_frobenius_residual = ', max_h2_relative
         write(unit,'(a,es24.16)') 'H2_max_hermiticity_residual = ', max_h2_hermiticity
         write(unit,'(a,es24.16)') 'H2_fourier_eeo_vs_h_o_residual = ', max_eeo_identity
         write(unit,'(a,es24.16)') 'potential_parameter_average_identity_residual = ', max_parameter_average
         write(unit,'(a,es24.16)') 'obarm_parameter_reconstruction_residual = ', max_obarm_provenance
         write(unit,'(a,es24.16)') 'enim_parameter_reconstruction_residual = ', max_enim_provenance
         write(unit,'(a,es24.16)') 'spin_resolved_BH_max_residual = ', max_spin_field
         write(unit,'(a,es24.16)') 'spin_resolved_BH_max_relative_residual = ', max_spin_field_relative
         write(unit,'(a,es24.16)') 'first_order_tangent_vs_commutator = ', first_tangent_relative
         write(unit,'(a,es24.16)') 'overlap_tangent_vs_commutator = ', overlap_tangent_relative
         write(unit,'(a,es24.16)') 'Enu_tangent_vs_commutator = ', enu_tangent_relative
         write(unit,'(a,es24.16)') 'second_order_product_tangent_vs_commutator = ', product_tangent_relative
         write(unit,'(a,es24.16)') 'second_order_product_tangent_max_element = ', product_tangent_max
         write(unit,'(a,es24.16)') 'delta_Enu_weighted_mean_norm = ', weighted_scalar_mean(tangent_enu_norm, reciprocal_obj%k_weights)
         write(unit,'(a,es24.16)') 'delta_h_weighted_mean_norm = ', weighted_scalar_mean(tangent_h_norm, reciprocal_obj%k_weights)
         write(unit,'(a,es24.16)') 'delta_h_o_h_left_weighted_mean_norm = ', weighted_scalar_mean(tangent_dh_o_h_norm, reciprocal_obj%k_weights)
         write(unit,'(a,es24.16)') 'delta_o_middle_weighted_mean_norm = ', weighted_scalar_mean(tangent_h_do_h_norm, reciprocal_obj%k_weights)
         write(unit,'(a,es24.16)') 'delta_h_o_h_right_weighted_mean_norm = ', weighted_scalar_mean(tangent_h_o_dh_norm, reciprocal_obj%k_weights)
         write(unit,'(a,es24.16)') 'delta_H_exact_weighted_mean_norm = ', weighted_scalar_mean(tangent_exact_norm, reciprocal_obj%k_weights)
         write(unit,'(a,es24.16)') 'delta_Enu_max_norm = ', maxval(tangent_enu_norm)
         write(unit,'(a,es24.16)') 'delta_h_max_norm = ', maxval(tangent_h_norm)
         write(unit,'(a,es24.16)') 'delta_h_o_h_left_max_norm = ', maxval(tangent_dh_o_h_norm)
         write(unit,'(a,es24.16)') 'delta_o_middle_max_norm = ', maxval(tangent_h_do_h_norm)
         write(unit,'(a,es24.16)') 'delta_h_o_h_right_max_norm = ', maxval(tangent_h_o_dh_norm)
         write(unit,'(a,es24.16)') 'delta_H_exact_max_norm = ', maxval(tangent_exact_norm)
         write(unit,'(a)') '# B_H decomposition is B_E + B_h - B_hoh; norms are coefficient-orbital norms.'
         write(unit,'(a,es24.16)') 'B_E_weighted_norm = ', norm_e
         write(unit,'(a,es24.16)') 'B_h_weighted_norm = ', norm_h
         write(unit,'(a,es24.16)') 'B_hoh_weighted_norm = ', norm_hoh
         write(unit,'(a,es24.16)') 'B_H_weighted_norm = ', norm_bh
         write(unit,'(a,es24.16)') 'B_E_fraction_of_BH = ', norm_e/max(norm_bh,tiny(1.0_rp))
         write(unit,'(a,es24.16)') 'B_h_fraction_of_BH = ', norm_h/max(norm_bh,tiny(1.0_rp))
         write(unit,'(a,es24.16)') 'B_hoh_fraction_of_BH = ', norm_hoh/max(norm_bh,tiny(1.0_rp))
         write(unit,'(a,es24.16)') 'Rk_B_E = ', rk_be
         write(unit,'(a,es24.16)') 'Rk_B_h = ', rk_bh
         write(unit,'(a,es24.16)') 'Rk_B_hoh = ', rk_bhoh
         write(unit,'(a,es24.16)') 'Rk_B_H = ', rk_bh_total
         write(unit,'(a,es24.16)') 'Rk_delta_h = ', rk_dh
         write(unit,'(a,es24.16)') 'Rk_delta_H_complete = ', rk_dh2
         write(unit,'(a,es24.16)') 'B_E_Gamma_local_spherical_projection_relative = ', gamma_component_projection(1)
         write(unit,'(a,es24.16)') 'B_h_Gamma_local_spherical_projection_relative = ', gamma_component_projection(2)
         write(unit,'(a,es24.16)') 'B_hoh_Gamma_local_spherical_projection_relative = ', gamma_component_projection(3)
         write(unit,'(a,es24.16)') 'Gamma_BH_local_spherical_projection_relative = ', gamma_projection_relative
         write(unit,'(a,es24.16)') 'B_E_Gamma_same_l_orbital_offdiagonal = ', gamma_component_offdiag(1)
         write(unit,'(a,es24.16)') 'B_h_Gamma_same_l_orbital_offdiagonal = ', gamma_component_offdiag(2)
         write(unit,'(a,es24.16)') 'B_hoh_Gamma_same_l_orbital_offdiagonal = ', gamma_component_offdiag(3)
         write(unit,'(a,es24.16)') 'Gamma_BH_same_l_orbital_offdiagonal = ', gamma_offdiag
         write(unit,'(a,es24.16)') 'B_E_Gamma_within_l_m_anisotropy = ', gamma_component_anisotropy(1)
         write(unit,'(a,es24.16)') 'B_h_Gamma_within_l_m_anisotropy = ', gamma_component_anisotropy(2)
         write(unit,'(a,es24.16)') 'B_hoh_Gamma_within_l_m_anisotropy = ', gamma_component_anisotropy(3)
         write(unit,'(a,es24.16)') 'Gamma_BH_within_l_m_anisotropy = ', gamma_anisotropy
         write(unit,'(a,es24.16)') 'B_E_Gamma_cross_l = ', gamma_component_cross_l(1)
         write(unit,'(a,es24.16)') 'B_h_Gamma_cross_l = ', gamma_component_cross_l(2)
         write(unit,'(a,es24.16)') 'B_hoh_Gamma_cross_l = ', gamma_component_cross_l(3)
         write(unit,'(a,es24.16)') 'Gamma_BH_cross_l = ', gamma_cross_l
         write(unit,'(a,es24.16)') 'B_E_Gamma_intersite = ', gamma_component_intersite(1)
         write(unit,'(a,es24.16)') 'B_h_Gamma_intersite = ', gamma_component_intersite(2)
         write(unit,'(a,es24.16)') 'B_hoh_Gamma_intersite = ', gamma_component_intersite(3)
         write(unit,'(a,es24.16)') 'Gamma_BH_intersite = ', gamma_intersite
         write(unit,'(a,es24.16)') 'old_radial_Bxc_vs_native_matrix = ', radial_bxc_relative
         write(unit,'(a,es24.16)') 'old_radial_Bxc_vs_native_action = ', radial_bxc_action
         write(unit,'(a,es24.16)') 'old_radial_BKS_vs_native_matrix = ', radial_bks_relative
         write(unit,'(a,es24.16)') 'old_radial_BKS_vs_native_action = ', radial_bks_action
         write(unit,'(a,es24.16)') 'old_radial_difference_to_Eonsite = ', radial_to_e
         write(unit,'(a,es24.16)') 'old_radial_difference_to_Eplus_h = ', radial_to_eh
         write(unit,'(a,es24.16)') 'native_product_response_vs_exact_static = ', native_response_relative
         write(unit,'(a,es24.16)') 'native_product_response_vs_exact_rotation = ', native_static_relative
         do ik = 1, nk
            write(unit,'(a,i0,a,es24.16)') 'k_',ik,'_delta_Enu_norm = ',tangent_enu_norm(ik)
            write(unit,'(a,i0,a,es24.16)') 'k_',ik,'_delta_h_norm = ',tangent_h_norm(ik)
            write(unit,'(a,i0,a,es24.16)') 'k_',ik,'_delta_h_o_h_left_norm = ',tangent_dh_o_h_norm(ik)
            write(unit,'(a,i0,a,es24.16)') 'k_',ik,'_delta_o_middle_norm = ',tangent_h_do_h_norm(ik)
            write(unit,'(a,i0,a,es24.16)') 'k_',ik,'_delta_h_o_h_right_norm = ',tangent_h_o_dh_norm(ik)
            write(unit,'(a,i0,a,es24.16)') 'k_',ik,'_delta_H_exact_norm = ',tangent_exact_norm(ik)
            write(unit,'(a,i0,a,es24.16)') 'k_',ik,'_native_vs_commutator_relative = ',tangent_comm_residual(ik)
            write(unit,'(a,i0,a,es24.16)') 'k_',ik,'_BH_local_spherical_projection_relative = ',gamma_bh_k(ik)
            write(unit,'(a,i0,a,es24.16)') 'k_',ik,'_BH_same_l_orbital_offdiagonal = ',gamma_bh_offdiag_k(ik)
            write(unit,'(a,i0,a,es24.16)') 'k_',ik,'_BH_within_l_m_anisotropy = ',gamma_bh_anisotropy_k(ik)
            write(unit,'(a,i0,a,es24.16)') 'k_',ik,'_BH_cross_l = ',gamma_bh_cross_l_k(ik)
            write(unit,'(a,i0,a,es24.16)') 'k_',ik,'_BH_intersite = ',gamma_bh_intersite_k(ik)
         end do
         do ik = 1, neta
            write(unit,'(a,i0,a,es24.16)') 'eta_',ik,' = ',eta_values(ik)
            write(unit,'(a,i0,a,es24.16)') 'exact_rotation_eta_relative_',ik,' = ',eta_exact(ik)
            write(unit,'(a,i0,a,es24.16)') 'native_eta_relative_',ik,' = ',eta_native(ik)
            write(unit,'(a,i0,a,es24.16)') 'native_vs_exact_eta_relative_',ik,' = ',eta_native_exact(ik)
         end do
         write(unit,'(a,a)') 'classification = ', trim(classification)
         write(unit,'(a,a)') 'verdict = ', merge('PASS   ','BLOCKED',pass_h2 .and. pass_spin .and. pass_tangent .and. pass_response)
         close(unit)
      end if
      deallocate(response_spec_reference, response_native_reference)

      if (.not. (pass_h2 .and. pass_spin .and. pass_tangent .and. pass_response)) then
         if (rank == 0) write(*,'(a,a)') 'DRESP-08 verdict: BLOCKED — ', trim(classification)
      else if (rank == 0) then
         write(*,'(a)') 'DRESP-08 verdict: PASS — NATIVE_SECOND_ORDER_MAPPING_CLOSED'
      end if
      deallocate(site_types, o_types, e_types, d_ee, d_obarm, d_enim, h, h2, h2_recon, o, enu, eeo_k, d_h, d_h_comm, &
         d_h2, d_h2_comm, d_o, d_enu, d_o_comm, d_enu_comm, term_e, term_h, term_left, term_mid, term_right, b_e, b_h, b_hoh, b_h_native, &
         generator, h_up, h_down, o_up, o_down, e_up, e_down, h2_up, h2_down, b_gamma, best_gamma, best_coefficients, &
         b_direct, dfield_exact, onsite_rebuilt, onsite_rebuilt_enu, radial_values, bxc_orbital, bks_orbital, response_rot, response_spec, &
         response_retarded, &
         response_native, response_native_retarded, response_radial, response_radial_retarded, one_static, one_retarded, &
         eta_exact, eta_native, eta_native_exact, tangent_enu_norm, tangent_h_norm, tangent_dh_o_h_norm, tangent_h_do_h_norm, &
         tangent_h_o_dh_norm, tangent_exact_norm, tangent_comm_residual, gamma_bh_k, gamma_bh_offdiag_k, &
         gamma_bh_anisotropy_k, gamma_bh_cross_l_k, gamma_bh_intersite_k)
   end subroutine run_dresp08_native_second_order

   subroutine reciprocal_kpoint(reciprocal_obj, ik, array4d, result)
      type(reciprocal), intent(in) :: reciprocal_obj
      integer, intent(in) :: ik
      complex(rp), intent(in) :: array4d(:, :, :, :)
      complex(rp), intent(out) :: result(:, :)
      real(rp) :: kpoint(3)
      if (allocated(reciprocal_obj%k_workset%points)) then
         kpoint = reciprocal_obj%k_workset%points(:, ik)
      else
         kpoint = reciprocal_obj%k_points(:, ik)
      end if
      call reciprocal_obj%fourier_transform_array(array4d, kpoint, result)
   end subroutine reciprocal_kpoint

   subroutine collect_radial_values(states, values, use_xc)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(out) :: values(:, :)
      logical, intent(in) :: use_xc
      integer :: isite
      if (any(shape(values) /= [size(states), size(states(1)%r)])) error stop 'DRESP-08 radial values: shape mismatch'
      do isite = 1, size(states)
         if (use_xc) then
            values(isite,:) = states(isite)%bxc_pauli
         else
            values(isite,:) = 0.5_rp*states(isite)%total_ks_spin_splitting
         end if
      end do
   end subroutine collect_radial_values

   subroutine map_radial_values(radial_bases, values, orbital)
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      real(rp), intent(in) :: values(:, :)
      complex(rp), intent(out) :: orbital(:, :)
      integer :: site, iorb, ir, l, first, npoint
      real(rp) :: integral
      npoint = radial_bases(1)%npoint
      if (any(shape(orbital) /= [(radial_bases(1)%lmax+1)**2*size(radial_bases), &
         (radial_bases(1)%lmax+1)**2*size(radial_bases)])) error stop 'DRESP-08 radial map: shape mismatch'
      orbital = cmplx(0.0_rp,0.0_rp,rp)
      do site = 1, size(radial_bases)
         first = (site-1)*(radial_bases(1)%lmax+1)**2
         do iorb = 1, (radial_bases(1)%lmax+1)**2
            l = lmto_orbital_l(iorb)
            integral = 0.0_rp
            do ir = 2, npoint
               integral = integral + radial_simpson_weight(ir,npoint)*radial_bases(site)%mesh_a* &
                  (radial_bases(site)%rofi(ir)+radial_bases(site)%mesh_b)*radial_bases(site)%phi_large(ir,l+1,1)**2*values(site,ir)
            end do
            orbital(first+iorb,first+iorb) = cmplx(integral,0.0_rp,rp)
         end do
      end do
   end subroutine map_radial_values

   subroutine field_matrix_residual(reference, candidate, matrix_relative, action_relative)
      complex(rp), intent(in) :: reference(:, :), candidate(:, :)
      real(rp), intent(out) :: matrix_relative, action_relative
      complex(rp), allocatable :: diff(:,:), v(:), dv(:)
      integer :: i
      allocate(diff(size(reference,1),size(reference,2)),v(size(reference,1)),dv(size(reference,1)))
      diff = candidate-reference
      matrix_relative = dresp08_relative_residual(candidate,reference)
      action_relative = 0.0_rp
      do i = 1, size(reference,1)
         v = cmplx(0.0_rp,0.0_rp,rp); v(i)=cmplx(1.0_rp,0.0_rp,rp)
         dv=matmul(diff,v)
         action_relative=max(action_relative,sqrt(sum(abs(dv)**2))/max(sqrt(sum(abs(matmul(reference,v))**2)),tiny(1.0_rp)))
      end do
      deallocate(diff,v,dv)
   end subroutine field_matrix_residual

   subroutine weighted_field_metrics(fields, weights, norm_value, rk_value)
      complex(rp), intent(in) :: fields(:, :, :)
      real(rp), intent(in) :: weights(:)
      real(rp), intent(out) :: norm_value, rk_value
      complex(rp), allocatable :: mean(:, :)
      real(rp) :: variance, denominator, wsum
      integer :: ik
      wsum = sum(weights)
      allocate(mean(size(fields,1),size(fields,2)))
      mean = cmplx(0.0_rp,0.0_rp,rp)
      do ik=1,size(fields,3); mean=mean+weights(ik)/wsum*fields(:,:,ik); end do
      variance=0.0_rp; denominator=0.0_rp
      do ik=1,size(fields,3)
         variance=variance+weights(ik)/wsum*sum(abs(fields(:,:,ik)-mean)**2)
         denominator=denominator+weights(ik)/wsum*sum(abs(fields(:,:,ik))**2)
      end do
      norm_value=sqrt(denominator)
      rk_value=sqrt(variance)/max(sqrt(denominator),tiny(1.0_rp))
      deallocate(mean)
   end subroutine weighted_field_metrics

   real(rp) function weighted_scalar_mean(values, weights) result(value)
      real(rp), intent(in) :: values(:), weights(:)
      if (size(values) /= size(weights) .or. sum(weights) <= 0.0_rp) then
         error stop 'DRESP-08 scalar metric: inconsistent weights'
      end if
      value = sum(weights*values)/sum(weights)
   end function weighted_scalar_mean

   subroutine evaluate_native_response(product, reciprocal_obj, fields, eta, response_static, response_retarded, include_static)
      type(lmto_product_response_basis), intent(in) :: product
      type(reciprocal), intent(in) :: reciprocal_obj
      complex(rp), intent(in) :: fields(:, :, :)
      real(rp), intent(in) :: eta
      complex(rp), intent(out) :: response_static(:), response_retarded(:)
      logical, intent(in) :: include_static
      complex(rp), allocatable :: one_static(:), one_retarded(:)
      integer :: ik
      real(rp) :: wsum
      allocate(one_static(product%product_dimension),one_retarded(product%product_dimension))
      response_static=cmplx(0.0_rp,0.0_rp,rp); response_retarded=cmplx(0.0_rp,0.0_rp,rp)
      wsum=sum(reciprocal_obj%k_weights)
      do ik=1,size(fields,3)
         call exact_ks_product_field_response(product,fields(:,:,ik),reciprocal_obj%eigenvalues(:,ik), &
            reciprocal_obj%eigenvectors(:,:,ik),reciprocal_obj%fermi_level,reciprocal_obj%temperature,eta,one_static,one_retarded,include_static,0)
         if(include_static) response_static=response_static+reciprocal_obj%k_weights(ik)/wsum*one_static
         response_retarded=response_retarded+reciprocal_obj%k_weights(ik)/wsum*one_retarded
      end do
      deallocate(one_static,one_retarded)
   end subroutine evaluate_native_response

   subroutine evaluate_orbital_response(product, reciprocal_obj, orbital, eta, response_static, response_retarded)
      type(lmto_product_response_basis), intent(in) :: product
      type(reciprocal), intent(in) :: reciprocal_obj
      complex(rp), intent(in) :: orbital(:, :)
      real(rp), intent(in) :: eta
      complex(rp), intent(out) :: response_static(:), response_retarded(:)
      complex(rp), allocatable :: field(:,:), one_static(:), one_retarded(:)
      integer :: ik
      real(rp) :: wsum
      allocate(field(size(reciprocal_obj%hk_bulk,1),size(reciprocal_obj%hk_bulk,2)),one_static(product%product_dimension), &
         one_retarded(product%product_dimension))
      call exact_ks_build_transverse_field(orbital,(product%orbital_lmax+1)**2,product%nsite,field)
      response_static=cmplx(0.0_rp,0.0_rp,rp); response_retarded=cmplx(0.0_rp,0.0_rp,rp); wsum=sum(reciprocal_obj%k_weights)
      do ik=1,size(reciprocal_obj%hk_bulk,3)
         call exact_ks_product_field_response(product,field,reciprocal_obj%eigenvalues(:,ik),reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level,reciprocal_obj%temperature,eta,one_static,one_retarded,.true.,0)
         response_static=response_static+reciprocal_obj%k_weights(ik)/wsum*one_static
         response_retarded=response_retarded+reciprocal_obj%k_weights(ik)/wsum*one_retarded
      end do
      deallocate(field,one_static,one_retarded)
   end subroutine evaluate_orbital_response

   pure real(rp) function vector_relative(candidate, reference) result(value)
      complex(rp), intent(in) :: candidate(:), reference(:)
      value=sqrt(sum(abs(candidate-reference)**2))/max(sqrt(sum(abs(reference)**2)),tiny(1.0_rp))
   end function vector_relative

   function diagonal_field_from_e(field) result(value)
      complex(rp), intent(in) :: field(:,:)
      complex(rp), allocatable :: value(:,:)
      allocate(value(size(field,1),size(field,2)))
      value=cmplx(0.0_rp,0.0_rp,rp)
      value=field
   end function diagonal_field_from_e

end module lr_dresp08_bridge_mod
