!------------------------------------------------------------------------------
! DRESP-09V production density-moment tangent gate.
!
! This bridge consumes the accepted reciprocal eigensystem and the live
! spin_density producer contract.  It does not construct a generic X or alter
! the production SCF density.  The radial service is entered only after the
! coefficient-space M0/M1/M2 and endpoint-branch oracles have closed.
!------------------------------------------------------------------------------
module lr_dresp09v_bridge_mod

   use precision_mod, only: rp
   use mpi_mod, only: rank
   use basis_mod, only: nb
   use reciprocal_mod, only: reciprocal
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use radial_ground_state_mod, only: radial_ground_state
   use pauli_ground_state_projection_mod, only: compute_accepted_pauli_magnetization
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use response_basis_mapping_mod, only: response_super_index, response_unflatten_superindex
   use lr_response_space_mod, only: response_space_layout
   use spin_density_mod, only: spin_density
   use lr_exact_ks_ward_mod, only: exact_ks_spin_rotation_generator, exact_ks_fermi_occupation
   use lr_lmto_density_moment_tangent_mod, only: density_from_eigensystem, moment_matrices_from_eigensystem, &
      moment_tangent_product_rule, moment_tangent_frechet, endpoint_tangent_branches, endpoint_fixed_h_branches, &
      endpoint_matrices_from_moments, commutator_tangent, fill_production_spin_density_from_moments, &
      relative_matrix_residual
   use lr_dresp09s_scalar_relativistic_mod, only: sr_l0_source_components, sr_l0_density_from_matrix_hamiltonian, &
      sr_l0_density_from_endpoint_branches, sr_l0_density_tangent_from_hamiltonian
   implicit none
   private

   character(len=*), parameter, public :: dresp09v_pass_a = 'PRODUCTION_DENSITY_MOMENT_TANGENT_CLOSED'
   character(len=*), parameter, public :: dresp09v_pass_b = 'DENSITY_MOMENTS_CLOSED_RADIAL_TARGET_OPEN'
   character(len=*), parameter, public :: dresp09v_energy_failure = 'ENERGY_ENDPOINT_TANGENT_FAILURE'
   character(len=*), parameter, public :: dresp09v_spin_mapping_failure = 'SPIN_DENSITY_BLOCK_MAPPING_FAILURE'
   character(len=*), parameter, public :: dresp09v_explicit_x = 'EXPLICIT_COEFFICIENT_MAP_REQUIRED'
   character(len=*), parameter, public :: dresp09v_mixed = 'MIXED'

   public :: run_dresp09v_density_moment_tangent

contains

   subroutine run_dresp09v_density_moment_tangent(output_file, response_space, radial_bases, ground_states, reciprocal_obj, &
                                                  lattice_obj, hamiltonian_obj)
      character(len=*), intent(in) :: output_file
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(reciprocal), intent(inout) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj

      type(spin_density) :: production_density, mapped_density, delta_density, sd_plus, sd_minus
      complex(rp), allocatable :: h(:, :, :), rho(:, :, :), delta_h(:, :, :), delta_rho(:, :, :)
      complex(rp), allocatable :: moments_k(:, :, :), delta_m_product(:, :, :), delta_m_frechet(:, :, :)
      complex(rp), allocatable :: endpoint_k(:, :, :), fixed_k(:, :, :)
      complex(rp), allocatable :: moment_sum(:, :, :), delta_m_sum(:, :, :), frechet_sum(:, :, :)
      complex(rp), allocatable :: endpoint_sum(:, :, :), fixed_sum(:, :, :)
      complex(rp), allocatable :: moments_plus(:, :, :), moments_minus(:, :, :)
      complex(rp), allocatable :: endpoint_plus(:, :, :), endpoint_minus(:, :, :)
      complex(rp), allocatable :: rotation(:, :), rotated_vectors(:, :, :)
      complex(rp), allocatable :: field_raw(:), source_plus(:, :, :), source_minus(:, :, :)
      complex(rp), allocatable :: old_plus(:, :), old_minus(:, :), complete_plus(:, :), complete_minus(:, :)
      complex(rp), allocatable :: fixed_plus(:, :), fixed_minus(:, :), endpoint_h_plus(:, :), endpoint_h_minus(:, :)
      complex(rp), allocatable :: radial_work(:, :), radial_complete_work(:, :)
      real(rp), allocatable :: old_weighted(:, :), complete_weighted(:, :), endpoint_h_weighted(:, :)
      real(rp), allocatable :: target_valence(:, :), target_total(:, :), core_magnetization(:, :)
      real(rp), allocatable :: magnetization(:, :), valence_magnetization(:, :)
      real(rp), allocatable :: target_l(:, :, :), complete_l(:, :, :)
      real(rp) :: equilibrium_commutator, moment_product_frechet, moment_product_equilibrium
      real(rp) :: commutator_moment(3), endpoint_identity(4), endpoint_d10_d01
      real(rp) :: production_mapping, production_fd(3), finite_rotation(3)
      real(rp) :: old_residual, complete_residual, integrated_residual, max_relative, max_absolute
      real(rp) :: endpoint_h_norm, endpoint_h_relative, field_pairing_fixed, field_pairing_h, field_pairing_total
      real(rp) :: pairing_fixed_trace, pairing_h_trace, pairing_total_trace, pairing_total_residual
      real(rp) :: pairing_fixed_residual, pairing_h_residual
      real(rp) :: old_max_relative, old_max_absolute
      real(rp) :: per_l(3), core_integral, target_integral, complete_integral, old_integral, wsum
      real(rp) :: sqrt_four_pi, r, ir_weight, theta, theta_ratio1, theta_ratio2
      integer :: nsite, nmat, norb_site, nk, ik, ik_global, n, iorder, branch, theta_index
      integer :: site, ir, unit, ios
      logical :: coefficient_closed, mapping_closed, radial_closed
      logical :: pass_gate
      character(len=96) :: classification
      character(len=16) :: verdict
      complex(rp) :: pairing_fixed_complex, pairing_h_complex, pairing_total_complex

      nsite = size(ground_states)
      if (nsite < 1 .or. lattice_obj%nrec /= nsite .or. size(radial_bases) /= nsite .or. response_space%nsite /= nsite) then
         error stop 'DRESP-09V: radial/response site dimensions are inconsistent'
      end if
      if (.not. allocated(reciprocal_obj%hk_bulk) .or. .not. allocated(reciprocal_obj%eigenvalues) .or. &
          .not. allocated(reciprocal_obj%eigenvectors) .or. .not. allocated(reciprocal_obj%k_weights)) then
         error stop 'DRESP-09V: accepted reciprocal eigensystem is incomplete'
      end if
      if (.not. hamiltonian_obj%hoh .or. hamiltonian_obj%ccor_2c .or. hamiltonian_obj%orb_pol .or. &
          hamiltonian_obj%hubbard_u_general_check .or. hamiltonian_obj%hubbard_v_check .or. hamiltonian_obj%local_axis) then
         error stop 'DRESP-09V: certified second-order ham_only state required'
      end if

      nmat = size(reciprocal_obj%eigenvectors, 1)
      nk = size(reciprocal_obj%eigenvalues, 2)
      if (nmat /= nb*nsite .or. mod(nmat, 2*nsite) /= 0 .or. nk /= 64 .or. size(reciprocal_obj%hk_bulk, 3) /= nk) then
         error stop 'DRESP-09V: accepted Fe 4x4x4 eigensystem required'
      end if
      norb_site = nmat/(2*nsite)
      if (.not. radial_provenance_complete(radial_bases)) error stop 'DRESP-09V: radial provenance is incomplete'
      if (size(reciprocal_obj%k_weights) /= nk) error stop 'DRESP-09V: k weights do not cover the accepted mesh'

      ! Refill the actual production object from the accepted final eigensystem.
      ! This is a diagnostic refresh only; no radial density is fed back to SCF.
      call reciprocal_obj%accumulate_spin_density_kspace()
      production_density = reciprocal_obj%rf_density
      mapped_density = spin_density(nsite, 4)

      allocate(h(nmat, nmat, nk), rho(nmat, nmat, nk), delta_h(nmat, nmat, nk), delta_rho(nmat, nmat, nk), &
         moments_k(nmat, nmat, 3), delta_m_product(nmat, nmat, 3), delta_m_frechet(nmat, nmat, 3), &
         endpoint_k(nmat, nmat, 4), fixed_k(nmat, nmat, 4), moment_sum(nmat, nmat, 3), &
         delta_m_sum(nmat, nmat, 3), frechet_sum(nmat, nmat, 3), endpoint_sum(nmat, nmat, 4), fixed_sum(nmat, nmat, 4))
      allocate(rotation(nmat, nmat), rotated_vectors(nmat, size(reciprocal_obj%eigenvalues, 1), nk), &
         moments_plus(nmat, nmat, 3), moments_minus(nmat, nmat, 3), endpoint_plus(nmat, nmat, 4), endpoint_minus(nmat, nmat, 4))
      moment_sum = cmplx(0.0_rp, 0.0_rp, rp)
      delta_m_sum = cmplx(0.0_rp, 0.0_rp, rp)
      frechet_sum = cmplx(0.0_rp, 0.0_rp, rp)
      endpoint_sum = cmplx(0.0_rp, 0.0_rp, rp)
      fixed_sum = cmplx(0.0_rp, 0.0_rp, rp)
      equilibrium_commutator = 0.0_rp
      moment_product_frechet = 0.0_rp
      moment_product_equilibrium = 0.0_rp
      commutator_moment = 0.0_rp
      endpoint_identity = 0.0_rp
      endpoint_d10_d01 = 0.0_rp
      wsum = sum(reciprocal_obj%k_weights)

      do ik = 1, nk
         ik_global = ik
         if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
         h(:, :, ik) = reciprocal_obj%hk_bulk(:, :, ik)
         call density_from_eigensystem(reciprocal_obj%eigenvalues(:, ik), reciprocal_obj%eigenvectors(:, :, ik), &
            reciprocal_obj%fermi_level, reciprocal_obj%temperature, rho(:, :, ik))
         call moment_matrices_from_eigensystem(reciprocal_obj%eigenvalues(:, ik), reciprocal_obj%eigenvectors(:, :, ik), &
            reciprocal_obj%fermi_level, reciprocal_obj%temperature, moments_k)
         call commutator_tangent(rotation_generator(nmat, norb_site, nsite), h(:, :, ik), delta_h(:, :, ik))
         call commutator_tangent(rotation_generator(nmat, norb_site, nsite), rho(:, :, ik), delta_rho(:, :, ik))
         call moment_tangent_product_rule(h(:, :, ik), rho(:, :, ik), delta_h(:, :, ik), delta_rho(:, :, ik), delta_m_product)
         call moment_tangent_frechet(reciprocal_obj%eigenvalues(:, ik), reciprocal_obj%eigenvectors(:, :, ik), &
            delta_h(:, :, ik), reciprocal_obj%fermi_level, reciprocal_obj%temperature, delta_m_frechet)
         call endpoint_tangent_branches(h(:, :, ik), rho(:, :, ik), delta_h(:, :, ik), delta_rho(:, :, ik), endpoint_k)
         call endpoint_fixed_h_branches(h(:, :, ik), delta_rho(:, :, ik), fixed_k)

         moment_sum = moment_sum + reciprocal_obj%k_weights(ik_global)*moments_k
         delta_m_sum = delta_m_sum + reciprocal_obj%k_weights(ik_global)*delta_m_product
         frechet_sum = frechet_sum + reciprocal_obj%k_weights(ik_global)*delta_m_frechet
         endpoint_sum = endpoint_sum + reciprocal_obj%k_weights(ik_global)*endpoint_k
         fixed_sum = fixed_sum + reciprocal_obj%k_weights(ik_global)*fixed_k
         equilibrium_commutator = max(equilibrium_commutator, relative_matrix_residual(matmul(h(:, :, ik), rho(:, :, ik)), &
            matmul(rho(:, :, ik), h(:, :, ik))))
         moment_product_frechet = max(moment_product_frechet, maxval([relative_matrix_residual(delta_m_product(:, :, 1), delta_m_frechet(:, :, 1)), &
            relative_matrix_residual(delta_m_product(:, :, 2), delta_m_frechet(:, :, 2)), &
            relative_matrix_residual(delta_m_product(:, :, 3), delta_m_frechet(:, :, 3))]))
         moment_product_equilibrium = max(moment_product_equilibrium, maxval([&
            relative_matrix_residual(moments_k(:, :, 1), rho(:, :, ik)), &
            relative_matrix_residual(moments_k(:, :, 2), matmul(h(:, :, ik), rho(:, :, ik))), &
            relative_matrix_residual(moments_k(:, :, 3), matmul(h(:, :, ik), matmul(h(:, :, ik), rho(:, :, ik))))]))
         call commutator_tangent(rotation_generator(nmat, norb_site, nsite), moments_k(:, :, 1), rotation(:, :))
         commutator_moment(1) = max(commutator_moment(1), relative_matrix_residual(delta_m_product(:, :, 1), rotation))
         call commutator_tangent(rotation_generator(nmat, norb_site, nsite), moments_k(:, :, 2), rotation(:, :))
         commutator_moment(2) = max(commutator_moment(2), relative_matrix_residual(delta_m_product(:, :, 2), rotation))
         call commutator_tangent(rotation_generator(nmat, norb_site, nsite), moments_k(:, :, 3), rotation(:, :))
         commutator_moment(3) = max(commutator_moment(3), relative_matrix_residual(delta_m_product(:, :, 3), rotation))
         endpoint_identity(1) = max(endpoint_identity(1), relative_matrix_residual(endpoint_k(:, :, 1), delta_m_product(:, :, 1)))
         endpoint_identity(2) = max(endpoint_identity(2), relative_matrix_residual(endpoint_k(:, :, 2), delta_m_product(:, :, 2)))
         endpoint_identity(3) = max(endpoint_identity(3), relative_matrix_residual(endpoint_k(:, :, 3), delta_m_product(:, :, 2)))
         endpoint_identity(4) = max(endpoint_identity(4), relative_matrix_residual(endpoint_k(:, :, 4), delta_m_product(:, :, 3)))
         endpoint_d10_d01 = max(endpoint_d10_d01, relative_matrix_residual(endpoint_k(:, :, 2), endpoint_k(:, :, 3)))
      end do

      call fill_production_spin_density_from_moments(moment_sum, mapped_density)
      production_mapping = spin_density_residual(mapped_density, production_density)
      delta_density = spin_density(nsite, 4)
      sd_plus = spin_density(nsite, 4)
      sd_minus = spin_density(nsite, 4)
      call fill_production_spin_density_from_moments(delta_m_sum, delta_density)

      finite_rotation = 0.0_rp
      production_fd = 0.0_rp
      do theta_index = 1, 3
         theta = 1.0e-2_rp/2.0_rp**real(theta_index - 1, rp)
         call build_spin_rotation(nmat, norb_site, nsite, theta, rotation)
         moments_plus = cmplx(0.0_rp, 0.0_rp, rp)
         moments_minus = cmplx(0.0_rp, 0.0_rp, rp)
         endpoint_plus = cmplx(0.0_rp, 0.0_rp, rp)
         endpoint_minus = cmplx(0.0_rp, 0.0_rp, rp)
         do ik = 1, nk
            ik_global = ik
            if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
            rotated_vectors(:, :, ik) = matmul(rotation, reciprocal_obj%eigenvectors(:, :, ik))
            call moment_matrices_from_eigensystem(reciprocal_obj%eigenvalues(:, ik), rotated_vectors(:, :, ik), &
               reciprocal_obj%fermi_level, reciprocal_obj%temperature, moments_k)
            moments_plus = moments_plus + reciprocal_obj%k_weights(ik_global)*moments_k
            call endpoint_matrices_from_moments(moments_k, endpoint_k)
            endpoint_plus = endpoint_plus + reciprocal_obj%k_weights(ik_global)*endpoint_k
            call build_spin_rotation(nmat, norb_site, nsite, -theta, rotation)
            rotated_vectors(:, :, ik) = matmul(rotation, reciprocal_obj%eigenvectors(:, :, ik))
            call moment_matrices_from_eigensystem(reciprocal_obj%eigenvalues(:, ik), rotated_vectors(:, :, ik), &
               reciprocal_obj%fermi_level, reciprocal_obj%temperature, moments_k)
            moments_minus = moments_minus + reciprocal_obj%k_weights(ik_global)*moments_k
            call endpoint_matrices_from_moments(moments_k, endpoint_k)
            endpoint_minus = endpoint_minus + reciprocal_obj%k_weights(ik_global)*endpoint_k
            call build_spin_rotation(nmat, norb_site, nsite, theta, rotation)
         end do
         do iorder = 1, 3
            finite_rotation(theta_index) = max(finite_rotation(theta_index), relative_matrix_residual(&
               (moments_plus(:, :, iorder) - moments_minus(:, :, iorder))/(2.0_rp*theta), delta_m_sum(:, :, iorder)))
         end do
         do branch = 1, 4
            finite_rotation(theta_index) = max(finite_rotation(theta_index), relative_matrix_residual(&
               (endpoint_plus(:, :, branch) - endpoint_minus(:, :, branch))/(2.0_rp*theta), endpoint_sum(:, :, branch)))
         end do
         call fill_production_spin_density_from_moments(moments_plus, sd_plus)
         call fill_production_spin_density_from_moments(moments_minus, sd_minus)
         do iorder = 1, 3
            production_fd(iorder) = spin_density_order_residual(sd_plus, sd_minus, delta_density, iorder, theta)
         end do
      end do
      theta_ratio1 = finite_rotation(1)/max(finite_rotation(2), tiny(1.0_rp))
      theta_ratio2 = finite_rotation(2)/max(finite_rotation(3), tiny(1.0_rp))

      ! These coefficient-space gates must close before any radial
      ! reconstruction is entered.
      coefficient_closed = equilibrium_commutator < 2.0e-10_rp .and. moment_product_frechet < 2.0e-10_rp .and. &
         maxval(commutator_moment) < 2.0e-10_rp .and. maxval(endpoint_identity) < 2.0e-10_rp .and. endpoint_d10_d01 < 2.0e-10_rp .and. &
         finite_rotation(3) < 2.0e-6_rp .and. theta_ratio1 > 3.0_rp .and. theta_ratio2 > 3.0_rp
      mapping_closed = production_mapping < 2.0e-12_rp .and. maxval(production_fd) < 2.0e-6_rp
      pass_gate = coefficient_closed .and. mapping_closed

      call compute_accepted_pauli_magnetization(reciprocal_obj, lattice_obj%symbolic_atoms, lattice_obj%nbulk, &
         magnetization, valence_magnetization, core_magnetization)
      allocate(target_valence(nsite, response_space%npoint), target_total(nsite, response_space%npoint), &
         target_l(nsite, response_space%npoint, 3), complete_l(nsite, response_space%npoint, 3))
      sqrt_four_pi = sqrt(4.0_rp*acos(-1.0_rp))
      do site = 1, nsite
         do ir = 1, response_space%npoint
            r = ground_states(site)%r(ir)
            target_valence(site, ir) = 4.0_rp*acos(-1.0_rp)*r*r*valence_magnetization(site, ir)
            target_total(site, ir) = ground_states(site)%rho_weighted_up(ir) - ground_states(site)%rho_weighted_down(ir)
         end do
      end do
      call build_valence_target_by_l(reciprocal_obj, lattice_obj, ground_states, target_l)

      allocate(old_plus(nsite, response_space%npoint), old_minus(nsite, response_space%npoint), &
         complete_plus(nsite, response_space%npoint), complete_minus(nsite, response_space%npoint), &
         fixed_plus(nsite, response_space%npoint), fixed_minus(nsite, response_space%npoint), &
         endpoint_h_plus(nsite, response_space%npoint), endpoint_h_minus(nsite, response_space%npoint), &
         radial_work(nsite, response_space%npoint), radial_complete_work(nsite, response_space%npoint), &
         old_weighted(nsite, response_space%npoint), &
         complete_weighted(nsite, response_space%npoint), endpoint_h_weighted(nsite, response_space%npoint), &
         field_raw(response_space%ndim), source_plus(nmat, nmat, 4), source_minus(nmat, nmat, 4))
      old_residual = 0.0_rp; complete_residual = 0.0_rp; integrated_residual = 0.0_rp
      max_relative = 0.0_rp; max_absolute = 0.0_rp; old_integral = 0.0_rp; complete_integral = 0.0_rp
      target_integral = 0.0_rp; endpoint_h_norm = 0.0_rp; endpoint_h_relative = 0.0_rp
      per_l = 0.0_rp; core_integral = 0.0_rp
      field_pairing_fixed = 0.0_rp; field_pairing_h = 0.0_rp; field_pairing_total = 0.0_rp
      pairing_fixed_trace = 0.0_rp; pairing_h_trace = 0.0_rp; pairing_total_trace = 0.0_rp
      pairing_fixed_residual = 0.0_rp; pairing_h_residual = 0.0_rp; pairing_total_residual = 0.0_rp
      pairing_fixed_complex = cmplx(0.0_rp, 0.0_rp, rp)
      pairing_h_complex = cmplx(0.0_rp, 0.0_rp, rp)
      pairing_total_complex = cmplx(0.0_rp, 0.0_rp, rp)

      if (pass_gate) then
         old_plus = cmplx(0.0_rp, 0.0_rp, rp); old_minus = cmplx(0.0_rp, 0.0_rp, rp)
         complete_plus = cmplx(0.0_rp, 0.0_rp, rp); complete_minus = cmplx(0.0_rp, 0.0_rp, rp)
         call sr_l0_density_from_endpoint_branches(response_space, radial_bases, fixed_sum/wsum, 1, fixed_plus)
         call sr_l0_density_from_endpoint_branches(response_space, radial_bases, fixed_sum/wsum, 2, fixed_minus)
         do ik = 1, nk
            ik_global = ik
            if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
            call sr_l0_density_tangent_from_hamiltonian(response_space, radial_bases, rho(:, :, ik), 1, h(:, :, ik), &
               delta_h(:, :, ik), delta_rho(:, :, ik), radial_complete_work)
            complete_plus = complete_plus + reciprocal_obj%k_weights(ik_global)*radial_complete_work
            call sr_l0_density_tangent_from_hamiltonian(response_space, radial_bases, rho(:, :, ik), 2, h(:, :, ik), &
               delta_h(:, :, ik), delta_rho(:, :, ik), radial_complete_work)
            complete_minus = complete_minus + reciprocal_obj%k_weights(ik_global)*radial_complete_work
            call sr_l0_density_from_matrix_hamiltonian(response_space, radial_bases, delta_rho(:, :, ik), 1, h(:, :, ik), radial_work)
            old_plus = old_plus + reciprocal_obj%k_weights(ik_global)*radial_work
            call sr_l0_density_from_matrix_hamiltonian(response_space, radial_bases, delta_rho(:, :, ik), 2, h(:, :, ik), radial_work)
            old_minus = old_minus + reciprocal_obj%k_weights(ik_global)*radial_work
         end do
         old_plus = old_plus/wsum; old_minus = old_minus/wsum
         complete_plus = complete_plus/wsum; complete_minus = complete_minus/wsum
         endpoint_h_plus = complete_plus - fixed_plus
         endpoint_h_minus = complete_minus - fixed_minus

         call build_l0_bxc_field(response_space, ground_states, field_raw)
         call sr_l0_source_components(response_space, radial_bases, field_raw, 1, source_plus)
         call sr_l0_source_components(response_space, radial_bases, field_raw, 2, source_minus)
         old_weighted = 0.0_rp
         complete_weighted = 0.0_rp
         endpoint_h_weighted = 0.0_rp
         do site = 1, nsite
            do ir = 1, response_space%npoint
               r = ground_states(site)%r(ir)
               old_weighted(site, ir) = sqrt_four_pi*r*r*real(old_plus(site, ir) + old_minus(site, ir), rp)
               complete_weighted(site, ir) = sqrt_four_pi*r*r*real(complete_plus(site, ir) + complete_minus(site, ir), rp)
               endpoint_h_weighted(site, ir) = complete_weighted(site, ir) - old_weighted(site, ir)
            end do
         end do
         call radial_metrics(complete_weighted, target_valence, ground_states, complete_residual, &
            max_relative, max_absolute, complete_integral)
         call radial_metrics(old_weighted, target_valence, ground_states, old_residual, old_max_relative, old_max_absolute, old_integral)
         target_integral = radial_integral(target_valence, ground_states)
         integrated_residual = abs(complete_integral - target_integral)/max(abs(target_integral), tiny(1.0_rp))
         endpoint_h_norm = radial_norm(endpoint_h_weighted, ground_states)
         endpoint_h_relative = endpoint_h_norm/max(radial_norm(old_weighted, ground_states), tiny(1.0_rp))
         do iorder = 1, 3
            call sr_l0_density_from_endpoint_branches(response_space, radial_bases, endpoint_sum/wsum, 1, radial_work, &
               l_first=iorder-1, l_last=iorder-1)
            call sr_l0_density_from_endpoint_branches(response_space, radial_bases, endpoint_sum/wsum, 2, radial_complete_work, &
               l_first=iorder-1, l_last=iorder-1)
            complete_l(:, :, iorder) = 0.0_rp
            do site = 1, nsite
               do ir = 2, response_space%npoint
                  r = ground_states(site)%r(ir)
                  complete_l(site, ir, iorder) = sqrt_four_pi*r*r*real(radial_work(site, ir) + radial_complete_work(site, ir), rp)
               end do
            end do
            per_l(iorder) = radial_single_relative(complete_l(:, :, iorder), target_l(:, :, iorder), ground_states)
         end do

         pairing_fixed_complex = cmplx(0.0_rp, 0.0_rp, rp)
         pairing_h_complex = cmplx(0.0_rp, 0.0_rp, rp)
         pairing_total_complex = cmplx(0.0_rp, 0.0_rp, rp)
         do branch = 1, 4
            pairing_fixed_complex = pairing_fixed_complex + sum(fixed_sum(:, :, branch)*transpose(source_plus(:, :, branch) + source_minus(:, :, branch)))
            pairing_h_complex = pairing_h_complex + sum((endpoint_sum(:, :, branch) - fixed_sum(:, :, branch))* &
               transpose(source_plus(:, :, branch) + source_minus(:, :, branch)))
            pairing_total_complex = pairing_total_complex + sum(endpoint_sum(:, :, branch)* &
               transpose(source_plus(:, :, branch) + source_minus(:, :, branch)))
         end do
         pairing_fixed_complex = pairing_fixed_complex/wsum
         pairing_h_complex = pairing_h_complex/wsum
         pairing_total_complex = pairing_total_complex/wsum
         field_pairing_fixed = radial_field_pairing(field_raw, response_space, old_plus + old_minus)
         field_pairing_h = radial_field_pairing(field_raw, response_space, endpoint_h_plus + endpoint_h_minus)
         field_pairing_total = radial_field_pairing(field_raw, response_space, complete_plus + complete_minus)
         pairing_fixed_trace = abs(pairing_fixed_complex)
         pairing_h_trace = abs(pairing_h_complex)
         pairing_total_trace = abs(pairing_total_complex)
         pairing_fixed_residual = abs(cmplx(field_pairing_fixed, 0.0_rp, rp) - pairing_fixed_complex)/ &
            max(abs(pairing_fixed_complex), tiny(1.0_rp))
         pairing_h_residual = abs(cmplx(field_pairing_h, 0.0_rp, rp) - pairing_h_complex)/ &
            max(abs(pairing_h_complex), tiny(1.0_rp))
         pairing_total_residual = abs(cmplx(field_pairing_total, 0.0_rp, rp) - pairing_total_complex)/ &
            max(abs(pairing_total_complex), tiny(1.0_rp))

         do site = 1, nsite
            do ir = 1, response_space%npoint
               r = ground_states(site)%r(ir)
               core_integral = core_integral + radial_simpson_jacobian(ground_states(site), ir)*4.0_rp*acos(-1.0_rp)*r*r* &
                  core_magnetization(site, ir)
            end do
         end do
      end if

      radial_closed = complete_residual < 1.0e-6_rp
      if (.not. coefficient_closed) then
         classification = dresp09v_energy_failure
         verdict = 'BLOCKED'
      else if (.not. mapping_closed) then
         classification = dresp09v_spin_mapping_failure
         verdict = 'BLOCKED'
      else if (radial_closed) then
         classification = dresp09v_pass_a
         verdict = 'PASS-A'
      else
         classification = dresp09v_pass_b
         verdict = 'PASS-B'
      end if
      pass_gate = coefficient_closed .and. mapping_closed

      if (rank == 0) then
         open(newunit=unit, file=trim(output_file), status='replace', action='write', iostat=ios)
         if (ios /= 0) error stop 'DRESP-09V: cannot open diagnostic artifact'
         write(unit, '(a)') '# DRESP-09V production density moment tangent'
         write(unit, '(a)') '# representation = accepted eigenvectors + M0/M1/M2 + spin_density; no generic X'
         write(unit, '(a,i0)') 'accepted_kpoints = ', nk
         write(unit, '(a,es24.16)') 'commutator_H_rho_residual = ', equilibrium_commutator
         write(unit, '(a,es24.16)') 'M0_product_rule_vs_Frechet = ', relative_matrix_residual(delta_m_sum(:, :, 1), frechet_sum(:, :, 1))
         write(unit, '(a,es24.16)') 'M1_product_rule_vs_Frechet = ', relative_matrix_residual(delta_m_sum(:, :, 2), frechet_sum(:, :, 2))
         write(unit, '(a,es24.16)') 'M2_product_rule_vs_Frechet = ', relative_matrix_residual(delta_m_sum(:, :, 3), frechet_sum(:, :, 3))
         write(unit, '(a,es24.16)') 'M0_equilibrium_product_residual = ', moment_product_equilibrium
         write(unit, '(a,es24.16)') 'rigid_delta_M0_commutator = ', commutator_moment(1)
         write(unit, '(a,es24.16)') 'rigid_delta_M1_commutator = ', commutator_moment(2)
         write(unit, '(a,es24.16)') 'rigid_delta_M2_commutator = ', commutator_moment(3)
         do branch = 1, 4
            write(unit, '(a,i0,a,es24.16)') 'endpoint_D', branch, '_identity = ', endpoint_identity(branch)
         end do
         write(unit, '(a,es24.16)') 'endpoint_D10_minus_D01 = ', endpoint_d10_d01
         write(unit, '(a,es24.16)') 'finite_angle_theta_1e-2 = ', finite_rotation(1)
         write(unit, '(a,es24.16)') 'finite_angle_theta_5e-3 = ', finite_rotation(2)
         write(unit, '(a,es24.16)') 'finite_angle_theta_2p5e-3 = ', finite_rotation(3)
         write(unit, '(a,es24.16)') 'finite_angle_theta_ratio_1 = ', theta_ratio1
         write(unit, '(a,es24.16)') 'finite_angle_theta_ratio_2 = ', theta_ratio2
         write(unit, '(a,es24.16)') 'production_spin_density_equilibrium_mapping = ', production_mapping
         write(unit, '(a,es24.16)') 'production_spin_density_M0_fd = ', production_fd(1)
         write(unit, '(a,es24.16)') 'production_spin_density_M1_fd = ', production_fd(2)
         write(unit, '(a,es24.16)') 'production_spin_density_M2_fd = ', production_fd(3)
         write(unit, '(a,es24.16)') 'old_DRESP09S_fixed_H_residual = ', old_residual
         write(unit, '(a,es24.16)') 'endpoint_H_contribution_norm = ', endpoint_h_norm
         write(unit, '(a,es24.16)') 'endpoint_H_contribution_relative_to_old = ', endpoint_h_relative
         write(unit, '(a,es24.16)') 'complete_DRESP09V_residual = ', complete_residual
         write(unit, '(a,es24.16)') 'integrated_valence_moment_residual = ', integrated_residual
         write(unit, '(a,es24.16)') 'maximum_radial_relative_residual = ', max_relative
         write(unit, '(a,es24.16)') 'maximum_radial_absolute_residual = ', max_absolute
         write(unit, '(a,es24.16)') 'per_l_s = ', per_l(1)
         write(unit, '(a,es24.16)') 'per_l_p = ', per_l(2)
         write(unit, '(a,es24.16)') 'per_l_d = ', per_l(3)
         write(unit, '(a,es24.16)') 'old_integrated_valence_moment = ', old_integral
         write(unit, '(a,es24.16)') 'complete_integrated_valence_moment = ', complete_integral
         write(unit, '(a,es24.16)') 'target_integrated_valence_moment = ', target_integral
         write(unit, '(a,es24.16)') 'core_diagnostic_integral = ', core_integral
         write(unit, '(a,es24.16)') 'field_density_fixed_H_pairing = ', field_pairing_fixed
         write(unit, '(a,es24.16)') 'field_density_endpoint_H_pairing = ', field_pairing_h
         write(unit, '(a,es24.16)') 'field_density_total_pairing = ', field_pairing_total
         write(unit, '(a,es24.16)') 'field_density_fixed_H_trace = ', pairing_fixed_trace
         write(unit, '(a,es24.16)') 'field_density_endpoint_H_trace = ', pairing_h_trace
         write(unit, '(a,es24.16)') 'field_density_total_trace = ', pairing_total_trace
         write(unit, '(a,es24.16)') 'field_density_fixed_H_pairing_residual = ', pairing_fixed_residual
         write(unit, '(a,es24.16)') 'field_density_endpoint_H_pairing_residual = ', pairing_h_residual
         write(unit, '(a,es24.16)') 'field_density_total_pairing_residual = ', pairing_total_residual
         write(unit, '(a,a)') 'endpoint_H_semantics = ', 'endpoint-energy representation tangent'
         write(unit, '(a,a)') 'explicit_coefficient_X_required = ', 'NO'
         write(unit, '(a,a)') 'DRESP09T_density_correction = ', 'RETIRED_NOT_USED'
         write(unit, '(a,a)') 'ALSDA_Ward = ', 'NOT RUN'
         write(unit, '(a,a)') 'BES_Halle = ', 'OFF'
         write(unit, '(a,a)') 'classification = ', trim(classification)
         write(unit, '(a,a)') 'verdict = ', trim(verdict)
         close(unit)
      end if
      if (rank == 0) write(*, '(a,a)') 'DRESP-09V verdict: ', trim(classification)

      deallocate(h, rho, delta_h, delta_rho, moments_k, delta_m_product, delta_m_frechet, endpoint_k, fixed_k, &
         moment_sum, delta_m_sum, frechet_sum, endpoint_sum, fixed_sum, rotation, rotated_vectors, moments_plus, moments_minus, &
         endpoint_plus, endpoint_minus, target_valence, target_total, target_l, complete_l, old_plus, old_minus, complete_plus, &
         complete_minus, fixed_plus, fixed_minus, endpoint_h_plus, endpoint_h_minus, radial_work, old_weighted, complete_weighted, &
         radial_complete_work, endpoint_h_weighted, field_raw, source_plus, source_minus)
   end subroutine run_dresp09v_density_moment_tangent

   function rotation_generator(nmat, norb_site, nsite) result(generator)
      integer, intent(in) :: nmat, norb_site, nsite
      complex(rp) :: generator(nmat, nmat)
      call exact_ks_spin_rotation_generator(norb_site, nsite, generator)
   end function rotation_generator

   subroutine build_spin_rotation(nmat, norb_site, nsite, theta, rotation)
      integer, intent(in) :: nmat, norb_site, nsite
      real(rp), intent(in) :: theta
      complex(rp), intent(out) :: rotation(:, :)
      complex(rp) :: generator(nmat, nmat)
      generator = rotation_generator(nmat, norb_site, nsite)
      rotation = cmplx(cos(theta/2.0_rp), 0.0_rp, rp)*identity_matrix(nmat) - &
         cmplx(0.0_rp, 2.0_rp*sin(theta/2.0_rp), rp)*generator
   end subroutine build_spin_rotation

   function identity_matrix(n) result(matrix)
      integer, intent(in) :: n
      complex(rp) :: matrix(n, n)
      integer :: i
      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, n
         matrix(i, i) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
   end function identity_matrix

   logical function radial_provenance_complete(radial) result(ok)
      type(lmto_radial_basis), intent(in) :: radial(:)
      integer :: isite
      ok = .true.
      do isite = 1, size(radial)
         ok = ok .and. allocated(radial(isite)%potential) .and. allocated(radial(isite)%tmc) .and. &
            allocated(radial(isite)%phi_small) .and. allocated(radial(isite)%phidot_small)
      end do
   end function radial_provenance_complete

   subroutine build_l0_bxc_field(space, states, field)
      type(response_space_layout), intent(in) :: space
      type(radial_ground_state), intent(in) :: states(:)
      complex(rp), intent(out) :: field(:)
      type(response_super_index) :: item
      integer :: flat
      real(rp) :: sqrt_four_pi
      if (size(field) /= space%ndim .or. size(states) /= space%nsite) error stop 'DRESP-09V field shape mismatch'
      sqrt_four_pi = sqrt(4.0_rp*acos(-1.0_rp))
      field = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, 1, item)
         if (item%response_l == 0 .and. item%response_m == 0) then
            field(flat) = cmplx(sqrt_four_pi*states(item%site)%bxc_pauli(item%radial_point), 0.0_rp, rp)
         end if
      end do
   end subroutine build_l0_bxc_field

   real(rp) function spin_density_residual(left, right) result(value)
      type(spin_density), intent(in) :: left, right
      value = sqrt(sum(abs(left%rho - right%rho)**2))/max(sqrt(sum(abs(right%rho)**2)), tiny(1.0_rp))
   end function spin_density_residual

   real(rp) function spin_density_order_residual(plus, minus, tangent, iorder, theta) result(value)
      type(spin_density), intent(in) :: plus, minus, tangent
      integer, intent(in) :: iorder
      real(rp), intent(in) :: theta
      complex(rp) :: finite_block(2, 2, size(plus%rho, 3), size(plus%rho, 4))
      finite_block = (plus%rho(:, :, :, :, iorder) - minus%rho(:, :, :, :, iorder))/(2.0_rp*theta)
      value = sqrt(sum(abs(finite_block - tangent%rho(:, :, :, :, iorder))**2))/ &
         max(sqrt(sum(abs(tangent%rho(:, :, :, :, iorder))**2)), tiny(1.0_rp))
   end function spin_density_order_residual

   subroutine radial_metrics(values, target, states, relative, aux1, aux2, integral)
      real(rp), intent(in) :: values(:, :), target(:, :)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(out) :: relative, aux1, aux2, integral
      integer :: site, ir
      real(rp) :: r, jac, num, den, difference
      num = 0.0_rp; den = 0.0_rp; aux1 = 0.0_rp; aux2 = 0.0_rp; integral = 0.0_rp
      do site = 1, size(states)
         do ir = 2, size(states(site)%r)
            r = states(site)%r(ir)
            jac = radial_simpson_jacobian(states(site), ir)
            difference = values(site, ir) - target(site, ir)
            num = num + jac*difference**2/(4.0_rp*acos(-1.0_rp)*r*r)
            den = den + jac*target(site, ir)**2/(4.0_rp*acos(-1.0_rp)*r*r)
            aux1 = max(aux1, abs(difference)/max(abs(target(site, ir)/(4.0_rp*acos(-1.0_rp)*r*r)), tiny(1.0_rp)))
            aux2 = max(aux2, abs(difference)/(4.0_rp*acos(-1.0_rp)*r*r))
            integral = integral + jac*values(site, ir)
         end do
      end do
      relative = sqrt(num/max(den, tiny(1.0_rp)))
   end subroutine radial_metrics

   real(rp) function radial_single_relative(values, target, states) result(relative)
      real(rp), intent(in) :: values(:, :), target(:, :)
      type(radial_ground_state), intent(in) :: states(:)
      integer :: site, ir
      real(rp) :: r, jac, num, den
      num = 0.0_rp; den = 0.0_rp
      do site = 1, size(states)
         do ir = 2, size(states(site)%r)
            r = states(site)%r(ir); jac = radial_simpson_jacobian(states(site), ir)
            num = num + jac*(values(site, ir)-target(site, ir))**2/(4.0_rp*acos(-1.0_rp)*r*r)
            den = den + jac*target(site, ir)**2/(4.0_rp*acos(-1.0_rp)*r*r)
         end do
      end do
      relative = sqrt(num/max(den, tiny(1.0_rp)))
   end function radial_single_relative

   real(rp) function radial_norm(values, states) result(norm)
      real(rp), intent(in) :: values(:, :)
      type(radial_ground_state), intent(in) :: states(:)
      integer :: site, ir
      real(rp) :: r, jac
      norm = 0.0_rp
      do site = 1, size(states)
         do ir = 2, size(states(site)%r)
            r = states(site)%r(ir); jac = radial_simpson_jacobian(states(site), ir)
            norm = norm + jac*values(site, ir)**2/(4.0_rp*acos(-1.0_rp)*r*r)
         end do
      end do
      norm = sqrt(norm)
   end function radial_norm

   real(rp) function radial_integral(values, states) result(value)
      real(rp), intent(in) :: values(:, :)
      type(radial_ground_state), intent(in) :: states(:)
      integer :: site, ir
      value = 0.0_rp
      do site = 1, size(states)
         do ir = 1, size(states(site)%r)
            value = value + radial_simpson_jacobian(states(site), ir)*values(site, ir)
         end do
      end do
   end function radial_integral

   real(rp) function radial_simpson_jacobian(state, ir) result(value)
      type(radial_ground_state), intent(in) :: state
      integer, intent(in) :: ir
      value = (2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp)*state%a*(state%r(ir) + state%b)
   end function radial_simpson_jacobian

   real(rp) function radial_field_pairing(field, space, density) result(value)
      complex(rp), intent(in) :: field(:), density(:, :)
      type(response_space_layout), intent(in) :: space
      type(response_super_index) :: item
      integer :: flat
      complex(rp) :: pairing
      pairing = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, 1, item)
         if (item%response_l == 0 .and. item%response_m == 0) then
            pairing = pairing + field(flat)*space%radial_weights(item%radial_point)*density(item%site, item%radial_point)
         end if
      end do
      value = real(pairing, rp)
   end function radial_field_pairing

   subroutine build_valence_target_by_l(reciprocal_obj, lattice_obj, states, target_l)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(out) :: target_l(:, :, :)
      integer :: site, ir, ik, ik_global, ib, ispin, iorb, l, atom_index, nmat, norb_site
      real(rp) :: wk, occupation, energy, radial_amplitude, sign_spin
      complex(rp) :: coefficient
      target_l = 0.0_rp
      nmat = size(reciprocal_obj%eigenvectors, 1); norb_site = nmat/(2*size(states))
      do site = 1, size(states)
         atom_index = lattice_obj%nbulk + site
         do ik = 1, size(reciprocal_obj%eigenvalues, 2)
            ik_global = ik
            if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
            wk = reciprocal_obj%k_weights(ik_global)
            do ib = 1, size(reciprocal_obj%eigenvalues, 1)
               energy = reciprocal_obj%eigenvalues(ib, ik)
               occupation = exact_ks_fermi_occupation(energy, reciprocal_obj%fermi_level, reciprocal_obj%temperature)
               if (occupation <= 1.0e-14_rp) cycle
               do ispin = 1, 2
                  sign_spin = merge(1.0_rp, -1.0_rp, ispin == 1)
                  do iorb = 1, norb_site
                     l = lmto_orbital_l(iorb)
                     if (l < 0 .or. l > 2 .or. l + 1 > size(target_l, 3)) cycle
                     coefficient = reciprocal_obj%eigenvectors((site - 1)*2*norb_site + (ispin - 1)*norb_site + iorb, ib, ik)
                     do ir = 2, size(states(site)%r)
                        radial_amplitude = states(site)%pauli_large(ir, l + 1, ispin) + &
                           (energy - states(site)%pauli_enu(l + 1, ispin))*states(site)%pauli_large_dot(ir, l + 1, ispin)
                        target_l(site, ir, l + 1) = target_l(site, ir, l + 1) + wk*occupation*sign_spin* &
                           real(coefficient*conjg(coefficient), rp)*radial_amplitude**2
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine build_valence_target_by_l

end module lr_dresp09v_bridge_mod
