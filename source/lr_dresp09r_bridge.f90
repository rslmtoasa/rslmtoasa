!------------------------------------------------------------------------------
! DRESP-09R accepted-state material gate.
!
! The bridge compares an independent Pauli large-component projection of the
! accepted radial bxc field with the already certified DRESP-08 native tangent.
! It is diagnostic-only: no SCF, chi0, Kxc, Ward, or Dyson state is changed.
!------------------------------------------------------------------------------
module lr_dresp09r_bridge_mod

   use precision_mod, only: rp
   use mpi_mod, only: rank
   use basis_mod, only: nb
   use reciprocal_mod, only: reciprocal
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use radial_ground_state_mod, only: radial_ground_state
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use response_basis_mapping_mod, only: response_super_index
   use lr_response_space_mod, only: response_space_layout
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, lmto_product_channel_plus, &
      lmto_product_channel_minus
   use lr_dresp09_field_insertion_mod, only: dresp09_compact_field_from_raw, dresp09_raw_field_from_compact, &
      dresp09_product_metric_adjoint_operator
   use lr_dresp09_pauli_direct_mod, only: dresp09_pauli_direct_source
   use lr_dresp08_native_mapping_mod, only: dresp08_build_product_tangent, dresp08_build_native_realspace_tangent, &
      dresp08_build_native_onsite_tangents, dresp08_block_diagonal, dresp08_commutator_tangent, dresp08_relative_residual
   use lr_exact_ks_ward_mod, only: exact_ks_spin_rotation_generator
   implicit none
   private

   character(len=*), parameter, public :: dresp09r_native_closed = 'PAULI_PROJECTED_NATIVE_ROTATION_CLOSED'
   character(len=*), parameter, public :: dresp09r_native_insufficient = &
      'PAULI_PROJECTION_INSUFFICIENT_FOR_NATIVE_TANGENT'
   character(len=*), parameter, public :: dresp09r_duality_blocked = 'BLOCKED — PAULI FIELD/DENSITY DUALITY'

   public :: run_dresp09r_pauli_projected_native

contains

   subroutine run_dresp09r_pauli_projected_native(output_file, response_space, radial_bases, ground_states, reciprocal_obj, &
                                                  lattice_obj, hamiltonian_obj)
      character(len=*), intent(in) :: output_file
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj

      type(lmto_product_response_basis) :: product_plus, product_minus
      complex(rp), allocatable :: vertices_plus(:, :, :, :), vertices_minus(:, :, :, :)
      complex(rp), allocatable :: h(:, :, :), h2(:, :, :), d_h(:, :, :), d_h2(:, :, :), o(:, :, :), enu(:, :, :)
      complex(rp), allocatable :: d_o(:, :, :), d_enu(:, :, :), d_ee(:, :, :, :)
      complex(rp), allocatable :: o_types(:, :, :), e_types(:, :, :), d_o_types(:, :, :), d_e_types(:, :, :)
      complex(rp), allocatable :: field_raw(:), compact_plus(:), compact_minus(:), projected_source(:)
      complex(rp), allocatable :: direct_plus(:, :), direct_minus(:, :), compact_op_plus(:, :), compact_op_minus(:, :)
      complex(rp), allocatable :: native(:, :), candidate(:, :), diff(:, :), generator(:, :)
      integer, allocatable :: site_types(:)
      real(rp), allocatable :: per_k_relative(:), per_k_max_element(:), occ_action(:), near_ef_action(:)
      real(rp) :: duality_plus, duality_minus, duality_residual
      real(rp) :: max_relative, weighted_rms, max_element, occupied_action, near_ef_action_max
      real(rp) :: bxc_norm, direct_norm, native_norm, source_projection_relative
      real(rp) :: native_commutator_max
      real(rp) :: wsum, denom, numerator
      integer :: nsite, norb_site, nmat, nk, ik, isite, ir, flat, unit, ios
      integer :: first, last, ib, nband, noccupied, nnear, ntype
      logical :: pass_duality, pass_native
      character(len=96) :: classification

      nsite = size(ground_states)
      if (nsite < 1 .or. size(radial_bases) /= nsite .or. response_space%nsite /= nsite) then
         error stop 'DRESP-09R: radial/response site dimensions are inconsistent'
      end if
      if (.not. allocated(reciprocal_obj%hk_bulk) .or. .not. allocated(reciprocal_obj%k_weights)) then
         error stop 'DRESP-09R: accepted reciprocal Hamiltonian is incomplete'
      end if
      if (lattice_obj%nrec /= nsite) error stop 'DRESP-09R: lattice/response site counts differ'
      if (.not. hamiltonian_obj%hoh) error stop 'DRESP-09R: accepted Hamiltonian is not on the HOH path'
      if (hamiltonian_obj%ccor_2c .or. hamiltonian_obj%orb_pol .or. hamiltonian_obj%hubbard_u_general_check .or. &
          hamiltonian_obj%hubbard_v_check .or. hamiltonian_obj%local_axis) then
         error stop 'DRESP-09R: excluded additive Hamiltonian correction is enabled'
      end if

      norb_site = (radial_bases(1)%lmax + 1)**2
      nmat = size(reciprocal_obj%hk_bulk, 1)
      nk = size(reciprocal_obj%hk_bulk, 3)
      if (nmat /= nb*nsite .or. nb /= 2*norb_site .or. size(reciprocal_obj%k_weights) /= nk) then
         error stop 'DRESP-09R: accepted basis/k-point dimensions differ'
      end if
      if (nk /= 64) error stop 'DRESP-09R: accepted Fe gate requires all 64 bcc-Fe k points'
      if (sum(reciprocal_obj%k_weights) <= 0.0_rp) error stop 'DRESP-09R: invalid accepted k weights'

      call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      call product_minus%initialize(response_space, radial_bases, lmto_product_channel_minus, .true.)
      if (product_plus%product_dimension /= product_minus%product_dimension) then
         error stop 'DRESP-09R: circular product dimensions differ'
      end if
      call product_plus%component_vertex_tensor(vertices_plus)
      call product_minus%component_vertex_tensor(vertices_minus)

      allocate(field_raw(response_space%ndim), compact_plus(product_plus%product_dimension), &
         compact_minus(product_minus%product_dimension), projected_source(response_space%ndim))
      call build_l0_bxc_field(response_space, ground_states, field_raw)
      call dresp09_compact_field_from_raw(response_space, product_plus, field_raw, compact_plus)
      call dresp09_compact_field_from_raw(response_space, product_minus, field_raw, compact_minus)

      allocate(direct_plus(nmat, nmat), direct_minus(nmat, nmat), compact_op_plus(nmat, nmat), &
         compact_op_minus(nmat, nmat), native(nmat, nmat), candidate(nmat, nmat), diff(nmat, nmat), &
         generator(nmat, nmat))
      allocate(per_k_relative(nk), per_k_max_element(nk), occ_action(nk), near_ef_action(nk))
      allocate(site_types(nsite))
      allocate(d_ee, source=hamiltonian_obj%ee)
      allocate(o_types(nb, nb, size(hamiltonian_obj%obarm, 3)), e_types(nb, nb, size(hamiltonian_obj%enim, 3)), &
         d_o_types(nb, nb, size(hamiltonian_obj%obarm, 3)), d_e_types(nb, nb, size(hamiltonian_obj%enim, 3)))
      o_types = hamiltonian_obj%obarm
      e_types = hamiltonian_obj%enim
      call dresp08_build_native_realspace_tangent(hamiltonian_obj, d_ee)
      call dresp08_build_native_onsite_tangents(hamiltonian_obj, d_o_types, d_e_types)
      site_types = lattice_obj%ib(1:nsite)
      call exact_ks_spin_rotation_generator(norb_site, nsite, generator)

      allocate(h(nmat, nmat, nk), h2(nmat, nmat, nk), d_h(nmat, nmat, nk), d_h2(nmat, nmat, nk), &
         o(nmat, nmat, nk), enu(nmat, nmat, nk), d_o(nmat, nmat, nk), d_enu(nmat, nmat, nk))
      duality_plus = 0.0_rp
      duality_minus = 0.0_rp
      native_commutator_max = 0.0_rp
      direct_norm = 0.0_rp
      native_norm = 0.0_rp
      do ik = 1, nk
         call reciprocal_kpoint(reciprocal_obj, ik, hamiltonian_obj%ee, h(:,:,ik))
         call reciprocal_kpoint(reciprocal_obj, ik, d_ee, d_h(:,:,ik))
         call dresp08_block_diagonal(o_types, site_types, o(:,:,ik))
         call dresp08_block_diagonal(e_types, site_types, enu(:,:,ik))
         call dresp08_block_diagonal(d_o_types, site_types, d_o(:,:,ik))
         call dresp08_block_diagonal(d_e_types, site_types, d_enu(:,:,ik))
         call dresp08_build_product_tangent(d_enu(:,:,ik), d_h(:,:,ik), d_o(:,:,ik), h(:,:,ik), o(:,:,ik), &
            d_h2(:,:,ik), candidate, diff, compact_op_plus, compact_op_minus, direct_plus)
         ! `hk_bulk` is the accepted H2 cache used by the DRESP-08 tangent.
         h2(:,:,ik) = reciprocal_obj%hk_bulk(:,:,ik)
         call dresp08_commutator_tangent(generator, h2(:,:,ik), native)
         native_commutator_max = max(native_commutator_max, dresp08_relative_residual(d_h2(:,:,ik), native))

         call dresp09_raw_field_from_compact(response_space, product_plus, conjg(compact_plus), projected_source)
         call dresp09_pauli_direct_source(response_space, radial_bases, projected_source, lmto_product_channel_plus, &
            h2(:,:,ik), compact_op_plus)
         call dresp09_product_metric_adjoint_operator(vertices_plus, h2(:,:,ik), compact_plus, candidate)
         duality_plus = max(duality_plus, dresp08_relative_residual(compact_op_plus, candidate))
         call dresp09_raw_field_from_compact(response_space, product_minus, conjg(compact_minus), projected_source)
         call dresp09_pauli_direct_source(response_space, radial_bases, projected_source, lmto_product_channel_minus, &
            h2(:,:,ik), compact_op_minus)
         call dresp09_product_metric_adjoint_operator(vertices_minus, h2(:,:,ik), compact_minus, candidate)
         duality_minus = max(duality_minus, dresp08_relative_residual(compact_op_minus, candidate))

         call dresp09_pauli_direct_source(response_space, radial_bases, field_raw, lmto_product_channel_plus, h2(:,:,ik), direct_plus)
         call dresp09_pauli_direct_source(response_space, radial_bases, field_raw, lmto_product_channel_minus, h2(:,:,ik), direct_minus)
         candidate = direct_plus + direct_minus
         diff = candidate - native
         per_k_relative(ik) = dresp08_relative_residual(candidate, native)
         per_k_max_element(ik) = maxval(abs(diff))
         call action_metrics(diff, native, reciprocal_obj%eigenvectors(:,:,ik), reciprocal_obj%eigenvalues(:,ik), &
            reciprocal_obj%fermi_level, occ_action(ik), near_ef_action(ik), noccupied, nnear)
         direct_norm = direct_norm + reciprocal_obj%k_weights(ik)*sum(abs(candidate)**2)
         native_norm = native_norm + reciprocal_obj%k_weights(ik)*sum(abs(native)**2)
      end do
      duality_residual = max(duality_plus, duality_minus)

      ! Calculate field and native scale metrics using the direct Bxc source.
      bxc_norm = sqrt(sum(abs(field_raw)**2*response_space%metric_weights))
      wsum = sum(reciprocal_obj%k_weights)
      numerator = 0.0_rp
      denom = 0.0_rp
      occupied_action = maxval(occ_action)
      near_ef_action_max = maxval(near_ef_action)
      max_relative = maxval(per_k_relative)
      max_element = maxval(per_k_max_element)
      do ik = 1, nk
         numerator = numerator + reciprocal_obj%k_weights(ik)*per_k_relative(ik)**2
         denom = denom + reciprocal_obj%k_weights(ik)
      end do
      weighted_rms = sqrt(numerator/max(denom, tiny(1.0_rp)))
      source_projection_relative = sqrt(sum(abs(field_raw - projected_source)**2*response_space%metric_weights))/ &
         max(bxc_norm, tiny(1.0_rp))
      pass_duality = duality_residual < 2.0e-10_rp
      pass_native = pass_duality .and. max_relative < 2.0e-10_rp
      if (.not. pass_duality) then
         classification = dresp09r_duality_blocked
      else if (.not. pass_native) then
         classification = dresp09r_native_insufficient
      else
         classification = dresp09r_native_closed
      end if

      if (rank == 0) then
         open(newunit=unit, file=trim(output_file), status='replace', action='write', iostat=ios)
         if (ios /= 0) error stop 'DRESP-09R: cannot open diagnostic artifact'
         write(unit,'(a)') '# DRESP-09R independent Pauli projected native-rotation gate'
         write(unit,'(a)') '# SCF = production accepted-state handoff; no SCF is performed here'
         write(unit,'(a)') '# field = sqrt(4*pi) * radial_ground_state%bxc_pauli (L=0, M=0)'
         write(unit,'(a)') '# direct_oracle = phi_large/phidot_large only; no small component and no GFAC'
         write(unit,'(a)') '# full_SR_route = BLOCKED - MIXED-SPIN RADIAL METRIC NOT CERTIFIED'
         write(unit,'(a,i0)') 'accepted_kpoints = ', nk
         write(unit,'(a,i0)') 'product_dimension_plus = ', product_plus%product_dimension
         write(unit,'(a,i0)') 'product_dimension_minus = ', product_minus%product_dimension
         write(unit,'(a,es24.16)') 'pauli_direct_vs_compact_adjoint_max_relative = ', duality_residual
         write(unit,'(a,es24.16)') 'pauli_source_projection_relative = ', source_projection_relative
         write(unit,'(a,es24.16)') 'bxc_raw_metric_norm = ', bxc_norm
         write(unit,'(a,es24.16)') 'pauli_direct_weighted_norm = ', sqrt(direct_norm/max(wsum,tiny(1.0_rp)))
         write(unit,'(a,es24.16)') 'native_tangent_weighted_norm = ', sqrt(native_norm/max(wsum,tiny(1.0_rp)))
         write(unit,'(a,es24.16)') 'dresp08_native_tangent_vs_commutator_max_relative = ', native_commutator_max
         write(unit,'(a,es24.16)') 'native_tangent_vs_pauli_direct_max_relative_frobenius = ', max_relative
         write(unit,'(a,es24.16)') 'native_tangent_vs_pauli_direct_weighted_rms_relative_frobenius = ', weighted_rms
         write(unit,'(a,es24.16)') 'native_tangent_vs_pauli_direct_max_element = ', max_element
         write(unit,'(a,es24.16)') 'native_tangent_vs_pauli_direct_occupied_action_max_relative = ', occupied_action
         write(unit,'(a,es24.16)') 'native_tangent_vs_pauli_direct_near_ef_action_max_relative = ', near_ef_action_max
         write(unit,'(a)') '# near_ef window = |epsilon-E_F| <= 0.10 Ry; occupied = epsilon <= E_F'
         do ik = 1, nk
            write(unit,'(a,i0,a,es24.16)') 'k_', ik, '_native_vs_pauli_relative_frobenius = ', per_k_relative(ik)
            write(unit,'(a,i0,a,es24.16)') 'k_', ik, '_native_vs_pauli_max_element = ', per_k_max_element(ik)
            write(unit,'(a,i0,a,es24.16)') 'k_', ik, '_occupied_action_max_relative = ', occ_action(ik)
            write(unit,'(a,i0,a,es24.16)') 'k_', ik, '_near_ef_action_max_relative = ', near_ef_action(ik)
         end do
         write(unit,'(a,a)') 'classification = ', trim(classification)
         write(unit,'(a,a)') 'verdict = ', merge('PASS   ','BLOCKED',pass_native)
         write(unit,'(a)') 'ward_stage = NOT RUN before native gate completion'
         write(unit,'(a)') 'alsda_stage = NOT RUN before native gate completion'
         write(unit,'(a)') 'full_sr_capability = BLOCKED - MIXED-SPIN RADIAL METRIC NOT CERTIFIED'
         close(unit)
      end if
      if (rank == 0) write(*,'(a,a)') 'DRESP-09R verdict: ', trim(classification)

      deallocate(vertices_plus, vertices_minus, h2, d_h, o, enu, d_o, d_enu, d_ee, o_types, e_types, d_o_types, d_e_types, &
         field_raw, compact_plus, compact_minus, projected_source, direct_plus, direct_minus, compact_op_plus, compact_op_minus, &
         h, native, candidate, diff, generator, d_h2, site_types, per_k_relative, per_k_max_element, occ_action, near_ef_action)
   end subroutine run_dresp09r_pauli_projected_native

   subroutine build_l0_bxc_field(space, states, field)
      type(response_space_layout), intent(in) :: space
      type(radial_ground_state), intent(in) :: states(:)
      complex(rp), intent(out) :: field(:)
      type(response_super_index) :: item
      integer :: flat
      real(rp) :: sqrt_four_pi

      if (size(field) /= space%ndim .or. size(states) /= space%nsite) error stop 'DRESP-09R: field shape mismatch'
      sqrt_four_pi = sqrt(4.0_rp*acos(-1.0_rp))
      field = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex_local(flat, item)
         if (item%response_l == 0 .and. item%response_m == 0) then
            field(flat) = cmplx(sqrt_four_pi*states(item%site)%bxc_pauli(item%radial_point), 0.0_rp, rp)
         end if
      end do
   contains
      subroutine response_unflatten_superindex_local(flat, item)
         integer, intent(in) :: flat
         type(response_super_index), intent(out) :: item
         integer :: zero_based, harmonic
         zero_based = flat - 1
         item%channel = mod(zero_based, space%nchannel) + 1
         zero_based = zero_based/space%nchannel
         item%radial_point = mod(zero_based, space%npoint) + 1
         zero_based = zero_based/space%npoint
         harmonic = mod(zero_based, (space%response_lmax + 1)**2) + 1
         item%site = zero_based/(space%response_lmax + 1)**2 + 1
         item%response_l = floor(sqrt(real(harmonic - 1, rp)))
         item%response_m = harmonic - item%response_l*item%response_l - item%response_l - 1
      end subroutine response_unflatten_superindex_local
   end subroutine build_l0_bxc_field

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

   subroutine action_metrics(diff, reference, eigenvectors, eigenvalues, fermi_level, occupied_value, near_ef_value, &
                             noccupied, nnear)
      complex(rp), intent(in) :: diff(:, :), reference(:, :), eigenvectors(:, :)
      real(rp), intent(in) :: eigenvalues(:), fermi_level
      real(rp), intent(out) :: occupied_value, near_ef_value
      integer, intent(out) :: noccupied, nnear
      integer :: ib
      real(rp) :: numerator, denominator

      occupied_value = 0.0_rp
      near_ef_value = 0.0_rp
      noccupied = 0
      nnear = 0
      do ib = 1, size(eigenvalues)
         numerator = sqrt(sum(abs(matmul(diff, eigenvectors(:, ib)))**2))
         denominator = max(sqrt(sum(abs(matmul(reference, eigenvectors(:, ib)))**2)), tiny(1.0_rp))
         if (eigenvalues(ib) <= fermi_level) then
            occupied_value = max(occupied_value, numerator/denominator)
            noccupied = noccupied + 1
         end if
         if (abs(eigenvalues(ib) - fermi_level) <= 0.10_rp) then
            near_ef_value = max(near_ef_value, numerator/denominator)
            nnear = nnear + 1
         end if
      end do
   end subroutine action_metrics

end module lr_dresp09r_bridge_mod
