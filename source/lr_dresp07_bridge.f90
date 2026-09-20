!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief DRESP-07 exact finite-H Kohn-Sham Ward-closure diagnostic.
!>
!> This is an intentionally standalone diagnostic seam.  It consumes the
!> accepted reciprocal eigensystem and the accepted LR-01 radial snapshot,
!> constructs the exact coefficient-space rigid rotation, and compares it to
!> the accepted Pauli magnetization and to radial local-field surrogates.
!> It does not enter Kxc, Goldstone correction, BES/Halle, or a Dyson solve.
!------------------------------------------------------------------------------
module lr_dresp07_bridge_mod

   use precision_mod, only: rp
   use mpi_mod, only: rank
   use reciprocal_mod, only: reciprocal
   use lattice_mod, only: lattice
   use radial_ground_state_mod, only: radial_ground_state, radial_simpson_weight
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use lr_response_space_mod, only: response_space_layout
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, lmto_product_channel_plus
   use lr_pauli_transition_vertex_mod, only: pauli_endpoint_state
   use lr_compact_static_interaction_mod, only: compact_project_magnetization, compact_project_local_operator
   use pauli_ground_state_projection_mod, only: compute_accepted_pauli_magnetization
   use lr_exact_ks_ward_mod, only: exact_ks_product_rotation_sweep, exact_ks_product_field_response, &
      exact_ks_extract_collinear_field, exact_ks_build_transverse_field
   implicit none
   private

   public :: run_dresp07_exact_ks_ward

contains

   subroutine run_dresp07_exact_ks_ward(output_file, response_space, radial_bases, ground_states, reciprocal_obj, &
                                        lattice_obj, eta_values)
      character(len=*), intent(in) :: output_file
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      real(rp), intent(in) :: eta_values(:)

      type(lmto_product_response_basis) :: product_plus
      real(rp), allocatable :: magnetization(:, :), valence_magnetization(:, :), core_magnetization(:, :)
      real(rp), allocatable :: bxc_values(:, :), bks_values(:, :), best_values(:, :)
      complex(rp), allocatable :: magnetization_product(:), valence_product(:), core_product(:)
      complex(rp), allocatable :: b_h_orbital(:, :), b_h_spin(:, :), bxc_orbital(:, :), bxc_spin(:, :)
      complex(rp), allocatable :: bks_orbital(:, :), bks_spin(:, :), best_orbital(:, :), best_spin(:, :)
      complex(rp), allocatable :: response_rot(:), response_spec(:), response_retarded(:)
      complex(rp), allocatable :: response_exact_field(:), response_exact_field_retarded(:)
      complex(rp), allocatable :: response_best(:), response_best_retarded(:)
      complex(rp), allocatable :: response_bks(:), response_bks_retarded(:)
      complex(rp), allocatable :: response_bxc(:), response_bxc_retarded(:)
      complex(rp), allocatable :: static_exact(:), static_best(:), static_bks(:), static_bxc(:)
      complex(rp), allocatable :: compact_bxc(:, :), compact_bks(:, :), compact_best(:, :)
      real(rp), allocatable :: eta_exact(:), eta_best(:), eta_bks(:), eta_bxc(:)
      real(rp) :: max_eigenpair_residual, max_density_residual, max_density_relative
      real(rp) :: max_eigenpair_frobenius, max_hermiticity, max_rotation_field, max_spin_offdiagonal
      real(rp) :: product_rotation_relative, exact_field_rotation_relative
      real(rp) :: r_mh_val, r_exact_h_val, r_exact_h_total, r_best_h_val, r_bks_val, r_bxc_val
      real(rp) :: valence_norm, total_norm, core_norm
      real(rp) :: b_h_onsite, b_h_nonlocal, b_h_spin_offdiag
      real(rp) :: b_h_orbital_offdiag, b_h_l_same, b_h_l_cross, b_h_k_norm, b_h_k_min, b_h_k_max
      real(rp) :: b_h_site_tmp, b_h_nonlocal_tmp, b_h_spin_tmp
      real(rp) :: bxc_matrix_relative, bxc_max_relative, bxc_action_relative
      real(rp) :: bks_matrix_relative, bks_max_relative, bks_action_relative
      real(rp) :: best_matrix_relative, best_max_relative, best_action_relative
      real(rp) :: best_rank, best_condition, best_captured_norm, best_residual_norm, best_relative_residual
      real(rp) :: compact_bxc_relative, compact_bks_relative, compact_best_relative
      real(rp) :: exact_h_response_norm
      integer :: nsite, norb_site, norb, nmat, neta, ik, unit, ios
      integer :: best_rank_int
      character(len=64) :: classification
      logical :: exact_gate_pass
      character(len=7) :: verdict

      nsite = size(ground_states)
      if (nsite < 1 .or. size(radial_bases) /= nsite .or. response_space%nsite /= nsite) then
         error stop 'DRESP-07: radial/response site dimensions are inconsistent'
      end if
      if (.not. allocated(reciprocal_obj%hk_bulk) .or. .not. allocated(reciprocal_obj%eigenvalues) .or. &
          .not. allocated(reciprocal_obj%eigenvectors) .or. .not. allocated(reciprocal_obj%k_weights)) then
         error stop 'DRESP-07: accepted reciprocal Hamiltonian/eigensystem is incomplete'
      end if
      nmat = size(reciprocal_obj%hk_bulk, 1)
      norb_site = (radial_bases(1)%lmax + 1)**2
      norb = norb_site*nsite
      if (nmat /= 2*norb) error stop 'DRESP-07: accepted H basis is not site-major Pauli spd'
      if (size(reciprocal_obj%hk_bulk, 3) /= size(reciprocal_obj%eigenvalues, 2) .or. &
          any(shape(reciprocal_obj%eigenvectors) /= [nmat, size(reciprocal_obj%eigenvalues, 1), &
          size(reciprocal_obj%eigenvalues, 2)])) then
         error stop 'DRESP-07: accepted H and eigensystem shapes differ'
      end if

      call product_plus%initialize(response_space, radial_bases, lmto_product_channel_plus, .true.)
      if (product_plus%product_dimension /= 232) then
         error stop 'DRESP-07: the campaign requires the complete Fe spd product dimension 232'
      end if
      call compute_accepted_pauli_magnetization(reciprocal_obj, lattice_obj%symbolic_atoms, lattice_obj%nbulk, &
         magnetization, valence_magnetization, core_magnetization)
      allocate(magnetization_product(product_plus%product_dimension), valence_product(product_plus%product_dimension), &
         core_product(product_plus%product_dimension))
      call compact_project_magnetization(response_space, product_plus, magnetization, magnetization_product)
      call compact_project_magnetization(response_space, product_plus, valence_magnetization, valence_product)
      call compact_project_magnetization(response_space, product_plus, core_magnetization, core_product)
      valence_norm = vector_norm(valence_product)
      total_norm = vector_norm(magnetization_product)
      core_norm = vector_norm(core_product)

      allocate(b_h_orbital(norb, norb), b_h_spin(nmat, nmat))
      call exact_ks_extract_collinear_field(reciprocal_obj%hk_bulk(:, :, 1), norb_site, nsite, b_h_orbital, b_h_spin, &
         b_h_onsite, b_h_nonlocal, b_h_spin_offdiag)
      call decompose_orbital_field(b_h_orbital, norb_site, nsite, b_h_orbital_offdiag, b_h_l_same, b_h_l_cross)
      b_h_k_min = huge(1.0_rp)
      b_h_k_max = 0.0_rp
      do ik = 1, size(reciprocal_obj%hk_bulk, 3)
         call exact_ks_extract_collinear_field(reciprocal_obj%hk_bulk(:, :, ik), norb_site, nsite, b_h_orbital, b_h_spin, &
            b_h_site_tmp, b_h_nonlocal_tmp, b_h_spin_tmp)
         b_h_k_norm = sqrt(sum(abs(b_h_orbital)**2))
         b_h_k_min = min(b_h_k_min, b_h_k_norm)
         b_h_k_max = max(b_h_k_max, b_h_k_norm)
      end do
      call exact_ks_extract_collinear_field(reciprocal_obj%hk_bulk(:, :, 1), norb_site, nsite, b_h_orbital, b_h_spin, &
         b_h_onsite, b_h_nonlocal, b_h_spin_offdiag)
      allocate(bxc_values(nsite, response_space%npoint), bks_values(nsite, response_space%npoint))
      call collect_radial_field(ground_states, response_space%npoint, .true., bxc_values)
      call collect_radial_field(ground_states, response_space%npoint, .false., bks_values)
      allocate(bxc_spin(nmat, nmat), bks_spin(nmat, nmat), bxc_orbital(norb, norb), bks_orbital(norb, norb))
      call map_radial_field_to_spin_matrix(radial_bases, bxc_values, bxc_spin)
      call map_radial_field_to_spin_matrix(radial_bases, bks_values, bks_spin)
      call spin_field_to_orbital(bxc_spin, norb, bxc_orbital)
      call spin_field_to_orbital(bks_spin, norb, bks_orbital)

      allocate(best_values(nsite, response_space%npoint), best_spin(nmat, nmat), best_orbital(norb, norb))
      call fit_best_radial_field(radial_bases, b_h_orbital, best_values, best_rank_int, best_condition, &
         best_captured_norm, best_residual_norm, best_relative_residual)
      best_rank = real(best_rank_int, rp)
      call map_radial_field_to_spin_matrix(radial_bases, best_values, best_spin)
      call spin_field_to_orbital(best_spin, norb, best_orbital)

      call compare_fields(b_h_spin, bxc_spin, reciprocal_obj%eigenvectors(:, :, 1), bxc_matrix_relative, bxc_max_relative, &
         bxc_action_relative)
      call compare_fields(b_h_spin, bks_spin, reciprocal_obj%eigenvectors(:, :, 1), bks_matrix_relative, bks_max_relative, &
         bks_action_relative)
      call compare_fields(b_h_spin, best_spin, reciprocal_obj%eigenvectors(:, :, 1), best_matrix_relative, best_max_relative, &
         best_action_relative)
      allocate(compact_bxc(product_plus%product_dimension, product_plus%product_dimension), &
         compact_bks(product_plus%product_dimension, product_plus%product_dimension), &
         compact_best(product_plus%product_dimension, product_plus%product_dimension))
      call compact_project_local_operator(response_space, product_plus, cmplx(bxc_values, 0.0_rp, rp), compact_bxc)
      call compact_project_local_operator(response_space, product_plus, cmplx(bks_values, 0.0_rp, rp), compact_bks)
      call compact_project_local_operator(response_space, product_plus, cmplx(best_values, 0.0_rp, rp), compact_best)

      neta = size(eta_values)
      if (neta < 1) error stop 'DRESP-07: at least one eta is required'
      allocate(eta_exact(neta), eta_best(neta), eta_bks(neta), eta_bxc(neta))
      allocate(response_rot(product_plus%product_dimension), response_spec(product_plus%product_dimension), &
         response_retarded(product_plus%product_dimension))
      allocate(response_exact_field(product_plus%product_dimension), response_exact_field_retarded(product_plus%product_dimension))
      allocate(response_best(product_plus%product_dimension), response_best_retarded(product_plus%product_dimension))
      allocate(response_bks(product_plus%product_dimension), response_bks_retarded(product_plus%product_dimension))
      allocate(response_bxc(product_plus%product_dimension), response_bxc_retarded(product_plus%product_dimension))
      call exact_ks_product_rotation_sweep(product_plus, reciprocal_obj%hk_bulk, reciprocal_obj%eigenvalues, &
         reciprocal_obj%eigenvectors, reciprocal_obj%k_weights, reciprocal_obj%fermi_level, reciprocal_obj%temperature, &
         eta_values(1), response_rot, response_spec, response_retarded, max_eigenpair_residual, max_density_residual, &
         max_density_relative, .false., 0, max_eigenpair_frobenius, max_hermiticity, max_rotation_field, max_spin_offdiagonal)
      exact_h_response_norm = vector_norm(response_spec)
      product_rotation_relative = vector_relative(response_rot, response_spec)
      call evaluate_exact_h_sweep(product_plus, reciprocal_obj, eta_values(1), response_exact_field, response_exact_field_retarded)
      exact_field_rotation_relative = vector_relative(response_exact_field, response_spec)
      call evaluate_field_sweep(product_plus, best_orbital, reciprocal_obj, eta_values(1), response_best, response_best_retarded, .true., 0)
      call evaluate_field_sweep(product_plus, bks_orbital, reciprocal_obj, eta_values(1), response_bks, response_bks_retarded, .true., 0)
      call evaluate_field_sweep(product_plus, bxc_orbital, reciprocal_obj, eta_values(1), response_bxc, response_bxc_retarded, .true., 0)
      allocate(static_exact(size(response_spec)), static_best(size(response_best)), static_bks(size(response_bks)), &
         static_bxc(size(response_bxc)))
      static_exact = response_spec
      static_best = response_best
      static_bks = response_bks
      static_bxc = response_bxc
      r_mh_val = vector_relative(response_rot, valence_product)
      r_exact_h_val = vector_relative(response_spec, valence_product)
      r_exact_h_total = vector_relative(response_spec, magnetization_product)
      r_best_h_val = vector_relative(response_best, valence_product)
      r_bks_val = vector_relative(response_bks, valence_product)
      r_bxc_val = vector_relative(response_bxc, valence_product)

      do ik = 1, neta
         call exact_ks_product_rotation_sweep(product_plus, reciprocal_obj%hk_bulk, reciprocal_obj%eigenvalues, &
            reciprocal_obj%eigenvectors, reciprocal_obj%k_weights, reciprocal_obj%fermi_level, reciprocal_obj%temperature, &
            eta_values(ik), response_rot, response_spec, response_retarded, max_eigenpair_residual, max_density_residual, &
            max_density_relative, .true., 0)
         eta_exact(ik) = vector_relative(response_retarded, static_exact)
         call evaluate_field_sweep(product_plus, best_orbital, reciprocal_obj, eta_values(ik), response_best, response_best_retarded, .false., 0)
         call evaluate_field_sweep(product_plus, bks_orbital, reciprocal_obj, eta_values(ik), response_bks, response_bks_retarded, .false., 0)
         call evaluate_field_sweep(product_plus, bxc_orbital, reciprocal_obj, eta_values(ik), response_bxc, response_bxc_retarded, .false., 0)
         eta_best(ik) = vector_relative(response_best_retarded, static_best)
         eta_bks(ik) = vector_relative(response_bks_retarded, static_bks)
         eta_bxc(ik) = vector_relative(response_bxc_retarded, static_bxc)
      end do

      exact_gate_pass = max_eigenpair_residual < 1.0e-9_rp .and. max_density_relative < 1.0e-9_rp .and. &
         product_rotation_relative < 1.0e-9_rp .and. exact_field_rotation_relative < 1.0e-9_rp
      if (.not. exact_gate_pass) then
         classification = 'EXACT_H_RESPONSE_FAILURE'
      else if (best_relative_residual > 1.0e-3_rp) then
         classification = 'FIELD_REPRESENTATION_FAILURE'
      else if (bxc_matrix_relative > 1.0e-2_rp .or. bks_matrix_relative > 1.0e-2_rp) then
         classification = 'RADIAL_TO_HAMILTONIAN_MAPPING_FAILURE'
      else if (core_norm > 1.0e-2_rp*max(valence_norm, tiny(1.0_rp)) .and. r_mh_val < 1.0e-6_rp) then
         classification = 'CORE_BOOKKEEPING_DOMINANT'
      else
         classification = 'WARD_CLOSURE_ACCEPTED'
      end if
      verdict = 'BLOCKED'
      if (exact_gate_pass) verdict = 'PASS'

      if (rank == 0) then
         open(newunit=unit, file=trim(output_file), status='replace', action='write', iostat=ios)
         if (ios /= 0) error stop 'DRESP-07: cannot open diagnostic artifact'
         write(unit, '(a)') '# DRESP-07 exact KS Ward closure'
         write(unit, '(a)') '# BES_Halle = OFF'
         write(unit, '(a)') '# Kxc_tuning = OFF'
         write(unit, '(a)') '# Goldstone_correction = OFF'
         write(unit, '(a)') '# accepted_state = reciprocal_kspace_scf'
         write(unit, '(a)') '# H_contract = accepted hk_bulk against accepted eigenpairs'
         write(unit, '(a,i0)') '# product_dimension = ', product_plus%product_dimension
         write(unit, '(a,es24.16)') 'max_eigenpair_residual = ', max_eigenpair_residual
         write(unit, '(a,es24.16)') 'max_eigenpair_frobenius_residual = ', max_eigenpair_frobenius
         write(unit, '(a,es24.16)') 'max_hermiticity_residual = ', max_hermiticity
         write(unit, '(a,es24.16)') 'max_rotation_field_relation_residual = ', max_rotation_field
         write(unit, '(a,es24.16)') 'max_spin_offdiagonal_norm = ', max_spin_offdiagonal
         write(unit, '(a,es24.16)') 'max_density_rotation_spectral_relative = ', max_density_relative
         write(unit, '(a,es24.16)') 'response_rotation_vs_spectral_relative = ', product_rotation_relative
         write(unit, '(a,es24.16)') 'response_exact_field_vs_spectral_relative = ', exact_field_rotation_relative
         write(unit, '(a,es24.16)') 'gamma_spin_offdiagonal_norm = ', b_h_spin_offdiag
         write(unit, '(a,es24.16)') 'gamma_BH_onsite_frobenius = ', b_h_onsite
         write(unit, '(a,es24.16)') 'gamma_BH_nonlocal_frobenius = ', b_h_nonlocal
         write(unit, '(a,es24.16)') 'gamma_BH_orbital_offdiagonal_frobenius = ', b_h_orbital_offdiag
         write(unit, '(a,es24.16)') 'gamma_BH_same_l_frobenius = ', b_h_l_same
         write(unit, '(a,es24.16)') 'gamma_BH_cross_l_frobenius = ', b_h_l_cross
         write(unit, '(a,es24.16)') 'BH_k_frobenius_min = ', b_h_k_min
         write(unit, '(a,es24.16)') 'BH_k_frobenius_max = ', b_h_k_max
         write(unit, '(a,es24.16)') 'accepted_magnetization_product_norm = ', total_norm
         write(unit, '(a,es24.16)') 'accepted_valence_magnetization_product_norm = ', valence_norm
         write(unit, '(a,es24.16)') 'accepted_core_magnetization_product_norm = ', core_norm
         write(unit, '(a,es24.16)') 'accepted_core_to_total_ratio = ', core_norm/max(total_norm, tiny(1.0_rp))
         write(unit, '(a,es24.16)') 'R_mH_val = ', r_mh_val
         write(unit, '(a,es24.16)') 'R_exactH_val = ', r_exact_h_val
         write(unit, '(a,es24.16)') 'R_exactH_total = ', r_exact_h_total
         write(unit, '(a,es24.16)') 'R_bestH_val = ', r_best_h_val
         write(unit, '(a,es24.16)') 'R_BKS_val = ', r_bks_val
         write(unit, '(a,es24.16)') 'R_Bxc_val = ', r_bxc_val
         write(unit, '(a,es24.16)') 'Bxc_matrix_relative = ', bxc_matrix_relative
         write(unit, '(a,es24.16)') 'Bxc_max_element_relative = ', bxc_max_relative
         write(unit, '(a,es24.16)') 'Bxc_action_relative = ', bxc_action_relative
         write(unit, '(a,es24.16)') 'BKS_matrix_relative = ', bks_matrix_relative
         write(unit, '(a,es24.16)') 'BKS_max_element_relative = ', bks_max_relative
         write(unit, '(a,es24.16)') 'BKS_action_relative = ', bks_action_relative
         write(unit, '(a,es24.16)') 'best_BH_matrix_relative = ', best_matrix_relative
         write(unit, '(a,es24.16)') 'best_BH_max_element_relative = ', best_max_relative
         write(unit, '(a,es24.16)') 'best_BH_action_relative = ', best_action_relative
         write(unit, '(a,i0)') 'best_radial_fit_rank = ', best_rank_int
         write(unit, '(a,es24.16)') 'best_radial_fit_condition = ', best_condition
         write(unit, '(a,es24.16)') 'best_BH_captured_frobenius = ', best_captured_norm
         write(unit, '(a,es24.16)') 'best_BH_residual_frobenius = ', best_residual_norm
         write(unit, '(a,es24.16)') 'best_BH_relative_residual = ', best_relative_residual
         write(unit, '(a,es24.16)') 'compact_Bxc_vs_BKS_relative = ', matrix_relative(compact_bxc, compact_bks)
         write(unit, '(a,es24.16)') 'compact_BHbest_vs_Bxc_relative = ', matrix_relative(compact_best, compact_bxc)
         write(unit, '(a,es24.16)') 'compact_BHbest_vs_BKS_relative = ', matrix_relative(compact_best, compact_bks)
         write(unit, '(a)') '# eta_ladder columns: eta exact_H_retarded_relative best_H_retarded_relative BKS_retarded_relative Bxc_retarded_relative'
         do ik = 1, neta
            write(unit, '(5(es24.16,1x))') eta_values(ik), eta_exact(ik), eta_best(ik), eta_bks(ik), eta_bxc(ik)
         end do
         write(unit, '(a)') '# ward_ladder columns: source static_residual_vs_accepted_valence'
         write(unit, '(a,es24.16)') 'accepted_magnetization ', 0.0_rp
         write(unit, '(a,es24.16)') 'exact_H ', r_exact_h_val
         write(unit, '(a,es24.16)') 'best_H ', r_best_h_val
         write(unit, '(a,es24.16)') 'B_KS ', r_bks_val
         write(unit, '(a,es24.16)') 'B_xc ', r_bxc_val
         write(unit, '(a,a)') 'classification = ', trim(classification)
         write(unit, '(a,a)') 'DRESP-07 verdict = ', trim(verdict)
         close(unit)
      end if
      deallocate(magnetization, valence_magnetization, core_magnetization, magnetization_product, valence_product, core_product)
      deallocate(bxc_values, bks_values, best_values, b_h_orbital, b_h_spin, bxc_orbital, bxc_spin, bks_orbital, bks_spin)
      deallocate(best_orbital, best_spin, compact_bxc, compact_bks, compact_best, response_rot, response_spec, response_retarded)
      deallocate(response_exact_field, response_exact_field_retarded, response_best, response_best_retarded, response_bks, &
         response_bks_retarded, response_bxc, response_bxc_retarded, static_exact, static_best, static_bks, static_bxc, &
         eta_exact, eta_best, eta_bks, eta_bxc)
   end subroutine run_dresp07_exact_ks_ward

   subroutine decompose_orbital_field(field, norb_site, nsite, offdiag_norm, same_l_norm, cross_l_norm)
      complex(rp), intent(in) :: field(:, :)
      integer, intent(in) :: norb_site, nsite
      real(rp), intent(out) :: offdiag_norm, same_l_norm, cross_l_norm
      integer :: isite, jsite, iorb, jorb, i, j, li, lj
      real(rp) :: offdiag_sum, same_l_sum, cross_l_sum

      offdiag_sum = 0.0_rp
      same_l_sum = 0.0_rp
      cross_l_sum = 0.0_rp
      do isite = 1, nsite
         do iorb = 1, norb_site
            i = (isite - 1)*norb_site + iorb
            do jsite = 1, nsite
               do jorb = 1, norb_site
                  j = (jsite - 1)*norb_site + jorb
                  li = lmto_orbital_l(iorb)
                  lj = lmto_orbital_l(jorb)
                  if (isite == jsite .and. iorb /= jorb) offdiag_sum = offdiag_sum + abs(field(i, j))**2
                  if (li == lj) then
                     same_l_sum = same_l_sum + abs(field(i, j))**2
                  else
                     cross_l_sum = cross_l_sum + abs(field(i, j))**2
                  end if
               end do
            end do
         end do
      end do
      offdiag_norm = sqrt(offdiag_sum)
      same_l_norm = sqrt(same_l_sum)
      cross_l_norm = sqrt(cross_l_sum)
   end subroutine decompose_orbital_field

   subroutine collect_radial_field(states, npoint, use_xc, values)
      type(radial_ground_state), intent(in) :: states(:)
      integer, intent(in) :: npoint
      logical, intent(in) :: use_xc
      real(rp), intent(out) :: values(:, :)
      integer :: isite

      if (any(shape(values) /= [size(states), npoint])) error stop 'DRESP-07: radial field shape mismatch'
      do isite = 1, size(states)
         if (use_xc) then
            if (.not. allocated(states(isite)%bxc_pauli)) error stop 'DRESP-07: Bxc field is not captured'
            values(isite, :) = states(isite)%bxc_pauli
         else
            if (.not. allocated(states(isite)%total_ks_spin_splitting)) error stop 'DRESP-07: KS splitting is not captured'
            values(isite, :) = 0.5_rp*states(isite)%total_ks_spin_splitting
         end if
      end do
   end subroutine collect_radial_field

   subroutine map_radial_field_to_spin_matrix(radial_bases, values, field)
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      real(rp), intent(in) :: values(:, :)
      complex(rp), intent(out) :: field(:, :)
      integer :: site, ir, iorb, l, up, down, norb_site, norb, nmat
      real(rp) :: integral

      norb_site = (radial_bases(1)%lmax + 1)**2
      norb = norb_site*size(radial_bases)
      nmat = 2*norb
      if (any(shape(field) /= [nmat, nmat]) .or. any(shape(values) /= [size(radial_bases), radial_bases(1)%npoint])) then
         error stop 'DRESP-07: radial-to-spin field shape mismatch'
      end if
      field = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, size(radial_bases)
         do iorb = 1, norb_site
            l = lmto_orbital_l(iorb)
            integral = 0.0_rp
            do ir = 2, radial_bases(site)%npoint
               integral = integral + radial_simpson_weight(ir, radial_bases(site)%npoint)*radial_bases(site)%mesh_a* &
                  (radial_bases(site)%rofi(ir) + radial_bases(site)%mesh_b)* &
                  radial_bases(site)%phi_large(ir, l + 1, 1)**2*values(site, ir)
            end do
            up = (site - 1)*2*norb_site + iorb
            down = up + norb_site
            field(up, up) = cmplx(integral, 0.0_rp, rp)
            field(down, down) = cmplx(-integral, 0.0_rp, rp)
         end do
      end do
   end subroutine map_radial_field_to_spin_matrix

   subroutine spin_field_to_orbital(field, norb, orbital)
      complex(rp), intent(in) :: field(:, :)
      integer, intent(in) :: norb
      complex(rp), intent(out) :: orbital(:, :)
      if (any(shape(field) /= [2*norb, 2*norb]) .or. any(shape(orbital) /= [norb, norb])) then
         error stop 'DRESP-07: spin/orbital field shape mismatch'
      end if
      orbital = field(1:norb, 1:norb)
   end subroutine spin_field_to_orbital

   function exact_ks_field_matrix(orbital, norb_site, nsite) result(field)
      complex(rp), intent(in) :: orbital(:, :)
      integer, intent(in) :: norb_site, nsite
      complex(rp), allocatable :: field(:, :)
      allocate(field(size(orbital, 1), size(orbital, 2)))
      field = orbital
      if (size(orbital, 1) /= norb_site*nsite) error stop 'DRESP-07: exact field shape mismatch'
   end function exact_ks_field_matrix

   subroutine evaluate_exact_h_sweep(product, reciprocal_obj, eta, response_static, response_retarded)
      type(lmto_product_response_basis), intent(in) :: product
      type(reciprocal), intent(in) :: reciprocal_obj
      real(rp), intent(in) :: eta
      complex(rp), intent(out) :: response_static(:), response_retarded(:)
      complex(rp), allocatable :: b_orbital(:, :), b_spin(:, :), delta_h(:, :), one_static(:), one_retarded(:)
      integer :: ik, nk, norb_site, norb, nmat
      real(rp) :: weight_sum, onsite_norm, nonlocal_norm, spin_offdiag_norm

      nk = size(reciprocal_obj%hk_bulk, 3)
      norb_site = (product%orbital_lmax + 1)**2
      norb = norb_site*product%nsite
      nmat = 2*norb
      weight_sum = sum(reciprocal_obj%k_weights)
      if (weight_sum <= 0.0_rp .or. eta <= 0.0_rp) error stop 'DRESP-07: invalid exact-H sweep weights/eta'
      allocate(b_orbital(norb, norb), b_spin(nmat, nmat), delta_h(nmat, nmat), &
         one_static(product%product_dimension), one_retarded(product%product_dimension))
      response_static = cmplx(0.0_rp, 0.0_rp, rp)
      response_retarded = cmplx(0.0_rp, 0.0_rp, rp)
      do ik = 1, nk
         call exact_ks_extract_collinear_field(reciprocal_obj%hk_bulk(:, :, ik), norb_site, product%nsite, b_orbital, b_spin, &
            onsite_norm, nonlocal_norm, spin_offdiag_norm)
         call exact_ks_build_transverse_field(b_orbital, norb_site, product%nsite, delta_h)
         call exact_ks_product_field_response(product, delta_h, reciprocal_obj%eigenvalues(:, ik), &
            reciprocal_obj%eigenvectors(:, :, ik), reciprocal_obj%fermi_level, reciprocal_obj%temperature, eta, one_static, &
            one_retarded, .true., 0)
         response_static = response_static + reciprocal_obj%k_weights(ik)/weight_sum*one_static
         response_retarded = response_retarded + reciprocal_obj%k_weights(ik)/weight_sum*one_retarded
      end do
      deallocate(b_orbital, b_spin, delta_h, one_static, one_retarded)
   end subroutine evaluate_exact_h_sweep

   subroutine evaluate_field_sweep(product, orbital_field, reciprocal_obj, eta, response_static, response_retarded, include_static, &
                                   response_l_filter)
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: orbital_field(:, :)
      type(reciprocal), intent(in) :: reciprocal_obj
      real(rp), intent(in) :: eta
      complex(rp), intent(out) :: response_static(:), response_retarded(:)
      logical, intent(in), optional :: include_static
      integer, intent(in), optional :: response_l_filter
      complex(rp), allocatable :: delta_h(:, :), one_static(:), one_retarded(:)
      integer :: ik, nk
      real(rp) :: weight_sum
      logical :: do_static

      nk = size(reciprocal_obj%hk_bulk, 3)
      weight_sum = sum(reciprocal_obj%k_weights)
      if (weight_sum <= 0.0_rp .or. eta <= 0.0_rp) error stop 'DRESP-07: invalid field sweep weights/eta'
      do_static = .true.
      if (present(include_static)) do_static = include_static
      allocate(delta_h(size(reciprocal_obj%hk_bulk, 1), size(reciprocal_obj%hk_bulk, 2)), &
         one_static(product%product_dimension), one_retarded(product%product_dimension))
      call exact_ks_build_transverse_field(orbital_field, (product%orbital_lmax + 1)**2, product%nsite, delta_h)
      response_static = cmplx(0.0_rp, 0.0_rp, rp)
      response_retarded = cmplx(0.0_rp, 0.0_rp, rp)
      do ik = 1, nk
         if (present(response_l_filter)) then
            call exact_ks_product_field_response(product, delta_h, reciprocal_obj%eigenvalues(:, ik), &
               reciprocal_obj%eigenvectors(:, :, ik), reciprocal_obj%fermi_level, reciprocal_obj%temperature, eta, one_static, &
               one_retarded, do_static, response_l_filter)
         else
            call exact_ks_product_field_response(product, delta_h, reciprocal_obj%eigenvalues(:, ik), &
               reciprocal_obj%eigenvectors(:, :, ik), reciprocal_obj%fermi_level, reciprocal_obj%temperature, eta, one_static, &
               one_retarded, do_static)
         end if
         if (do_static) response_static = response_static + reciprocal_obj%k_weights(ik)/weight_sum*one_static
         response_retarded = response_retarded + reciprocal_obj%k_weights(ik)/weight_sum*one_retarded
      end do
      deallocate(delta_h, one_static, one_retarded)
   end subroutine evaluate_field_sweep

   subroutine compare_fields(reference, candidate, eigenvectors, matrix_relative_value, max_relative, action_relative)
      complex(rp), intent(in) :: reference(:, :), candidate(:, :), eigenvectors(:, :)
      real(rp), intent(out) :: matrix_relative_value, max_relative, action_relative
      complex(rp), allocatable :: delta(:, :), ref_action(:), delta_action(:)
      integer :: ib
      real(rp) :: ref_norm

      allocate(delta(size(reference, 1), size(reference, 2)), ref_action(size(reference, 1)), delta_action(size(reference, 1)))
      delta = candidate - reference
      ref_norm = sqrt(sum(abs(reference)**2))
      matrix_relative_value = sqrt(sum(abs(delta)**2))/max(ref_norm, tiny(1.0_rp))
      max_relative = maxval(abs(delta))/max(maxval(abs(reference)), tiny(1.0_rp))
      action_relative = 0.0_rp
      do ib = 1, size(eigenvectors, 2)
         ref_action = matmul(reference, eigenvectors(:, ib))
         delta_action = matmul(delta, eigenvectors(:, ib))
         action_relative = max(action_relative, sqrt(sum(abs(delta_action)**2))/max(sqrt(sum(abs(ref_action)**2)), tiny(1.0_rp)))
      end do
      deallocate(delta, ref_action, delta_action)
   end subroutine compare_fields

   subroutine fit_best_radial_field(radial_bases, target, values, rank_out, condition, captured_norm, residual_norm, relative_residual)
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: target(:, :)
      real(rp), intent(out) :: values(:, :)
      integer, intent(out) :: rank_out
      real(rp), intent(out) :: condition, captured_norm, residual_norm, relative_residual
      real(rp), allocatable :: row_metric(:, :), target_l(:, :), row_norm(:)
      complex(rp), allocatable :: fitted(:, :)
      integer :: nsite, npoint, lmax, norb_site, norb, site, ir, iorb, l, nmode
      real(rp) :: metric, mode_norm, max_mode_norm, min_mode_norm, coefficient

      nsite = size(radial_bases)
      npoint = radial_bases(1)%npoint
      lmax = radial_bases(1)%lmax
      norb_site = (lmax + 1)**2
      norb = nsite*norb_site
      allocate(row_metric(nsite, 0:lmax), target_l(nsite, 0:lmax), row_norm(nsite*(lmax + 1)))
      row_metric = 0.0_rp
      target_l = 0.0_rp
      row_norm = 0.0_rp
      do site = 1, nsite
         do l = 0, lmax
            nmode = 0
            do iorb = 1, norb_site
               if (lmto_orbital_l(iorb) == l) then
                  nmode = nmode + 1
                  target_l(site, l) = target_l(site, l) + real(target((site - 1)*norb_site + iorb, &
                     (site - 1)*norb_site + iorb), rp)
               end if
            end do
            target_l(site, l) = target_l(site, l)/real(max(1, nmode), rp)
            mode_norm = 0.0_rp
            do ir = 2, npoint
               metric = radial_simpson_weight(ir, npoint)*radial_bases(site)%mesh_a* &
                  (radial_bases(site)%rofi(ir) + radial_bases(site)%mesh_b)* &
                  radial_bases(site)%phi_large(ir, l + 1, 1)**2
               row_metric(site, l) = row_metric(site, l) + metric*metric
               mode_norm = mode_norm + metric*metric
            end do
            row_norm((site - 1)*(lmax + 1) + l + 1) = sqrt(mode_norm*real(max(1, nmode), rp))
         end do
      end do
      values = 0.0_rp
      rank_out = 0
      max_mode_norm = maxval(row_norm)
      min_mode_norm = huge(1.0_rp)
      do site = 1, nsite
         do l = 0, lmax
            mode_norm = row_norm((site - 1)*(lmax + 1) + l + 1)
            if (mode_norm <= tiny(1.0_rp)) cycle
            rank_out = rank_out + 1
            min_mode_norm = min(min_mode_norm, mode_norm)
            coefficient = target_l(site, l)/row_metric(site, l)
            do ir = 2, npoint
               metric = radial_simpson_weight(ir, npoint)*radial_bases(site)%mesh_a* &
                  (radial_bases(site)%rofi(ir) + radial_bases(site)%mesh_b)* &
                  radial_bases(site)%phi_large(ir, l + 1, 1)**2
               values(site, ir) = values(site, ir) + coefficient*metric
            end do
         end do
      end do
      if (rank_out > 0) then
         condition = max_mode_norm/min_mode_norm
      else
         condition = huge(1.0_rp)
      end if
      allocate(fitted(norb, norb))
      call map_radial_field_to_orbital(radial_bases, values, fitted)
      captured_norm = sqrt(sum(abs(fitted)**2))
      residual_norm = sqrt(sum(abs(target - fitted)**2))
      relative_residual = residual_norm/max(sqrt(sum(abs(target)**2)), tiny(1.0_rp))
      deallocate(row_metric, target_l, row_norm, fitted)
   end subroutine fit_best_radial_field

   subroutine map_radial_field_to_orbital(radial_bases, values, orbital)
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      real(rp), intent(in) :: values(:, :)
      complex(rp), intent(out) :: orbital(:, :)
      complex(rp), allocatable :: spin(:, :)
      integer :: norb
      norb = ((radial_bases(1)%lmax + 1)**2)*size(radial_bases)
      allocate(spin(2*norb, 2*norb))
      call map_radial_field_to_spin_matrix(radial_bases, values, spin)
      orbital = spin(1:norb, 1:norb)
      deallocate(spin)
   end subroutine map_radial_field_to_orbital

   pure real(rp) function vector_norm(vector) result(value)
      complex(rp), intent(in) :: vector(:)
      value = sqrt(sum(abs(vector)**2))
   end function vector_norm

   pure real(rp) function vector_relative(candidate, reference) result(value)
      complex(rp), intent(in) :: candidate(:), reference(:)
      value = vector_norm(candidate - reference)/max(vector_norm(reference), tiny(1.0_rp))
   end function vector_relative

   pure real(rp) function matrix_relative(candidate, reference) result(value)
      complex(rp), intent(in) :: candidate(:, :), reference(:, :)
      value = sqrt(sum(abs(candidate - reference)**2))/max(sqrt(sum(abs(reference)**2)), tiny(1.0_rp))
   end function matrix_relative

end module lr_dresp07_bridge_mod
