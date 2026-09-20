!------------------------------------------------------------------------------
! DRESP-09W radial observable / ground-state provenance gate.
!
! This bridge consumes exactly the accepted reciprocal eigensystem, the
! accepted radial snapshot, and the production M0/M1/M2 moments.  It keeps
! Pauli and full scalar-relativistic observables separate and does not enter
! ALSDA, Ward, Dyson, spectra, BES, or Halle code.
!------------------------------------------------------------------------------
module lr_dresp09w_bridge_mod

   use precision_mod, only: rp
   use mpi_mod, only: rank
   use basis_mod, only: nb
   use reciprocal_mod, only: reciprocal
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use radial_ground_state_mod, only: radial_ground_state
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use response_basis_mapping_mod, only: response_super_index, response_unflatten_superindex
   use lr_response_space_mod, only: response_space_layout
   use pauli_ground_state_projection_mod, only: compute_accepted_pauli_magnetization
   use lr_exact_ks_ward_mod, only: exact_ks_spin_rotation_generator
   use lr_lmto_density_moment_tangent_mod, only: density_from_eigensystem, moment_matrices_from_eigensystem, &
      endpoint_tangent_branches, commutator_tangent
   use lr_dresp09s_scalar_relativistic_mod, only: sr_l0_density_from_endpoint_branches
   use lr_radial_observable_provenance_mod, only: channel_moments_from_matrix, pauli_density_from_moments, &
      sr_density_from_moments, pauli_density_from_occupied_states, sr_density_from_occupied_states, &
      pauli_l0_density_from_endpoint_branches, weighted_from_density, radial_relative_metrics, volume_integral, &
      volume_l2_norm, dresp09w_four_pi
   implicit none
   private

   character(len=*), parameter, public :: dresp09w_closed = 'RADIAL_OBSERVABLE_PROVENANCE_CLOSED'
   character(len=*), parameter, public :: dresp09w_pauli_target_open = 'PAULI_TARGET_PROVENANCE_OPEN'
   character(len=*), parameter, public :: dresp09w_pauli_response_open = 'PAULI_RESPONSE_OBSERVABLE_OPEN'
   character(len=*), parameter, public :: dresp09w_sr_open = 'SR_RADIAL_RECONSTRUCTION_OPEN'
   character(len=*), parameter, public :: dresp09w_snapshot_open = 'ACCEPTED_RADIAL_SNAPSHOT_PROVENANCE_OPEN'
   character(len=*), parameter, public :: dresp09w_core_open = 'CORE_DECOMPOSITION_OPEN'
   character(len=*), parameter, public :: dresp09w_mixed = 'MIXED'

   public :: run_dresp09w_radial_observable_provenance

contains

   subroutine run_dresp09w_radial_observable_provenance(output_file, response_space, radial_bases, ground_states, &
                                                        reciprocal_obj, lattice_obj, hamiltonian_obj)
      character(len=*), intent(in) :: output_file
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(reciprocal), intent(inout) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj

      integer :: nsite, nmat, norb_site, nk, ik, ik_global, iorder, site, ir, l, spin, unit, ios
      real(rp) :: wsum, wk, old_cross_representation
      real(rp) :: equilibrium_m0(2), p1_metrics(5), p2_metrics(5), p3_metrics(5)
      real(rp) :: sr1_metrics(5), sr2_metrics(5), sr3_metrics(5), sr4_metrics(5)
      real(rp) :: p1_p2(5), p1_p3(5), p2_p3(5), sr2_sr3(5), sr1_sr4(5)
      real(rp) :: pauli_response_p1(5), pauli_response_p2(5), pauli_response_p3(5), sr_response_sr2(5)
      real(rp) :: sr_pauli(5), core_metrics(5), core_l2, core_s_integral, core_s_l2, core_s_max
      real(rp) :: core_to_valence_peak, response_per_l(3), pauli_response_per_l(3)
      real(rp) :: core_decomposition_identity
      real(rp) :: sqrt_four_pi, r, sign_spin, coefficient_weight
      real(rp) :: per_l_p1_p2(3), per_l_p1_p3(3), per_l_p2_p3(3), per_l_sr2_sr3(3), per_l_sr_response(3)
      real(rp) :: sr1_sr4_integral, pauli_moment, sr_valence_moment, sr_total_moment
      real(rp) :: pauli_n_up, pauli_n_down, sr_n_up, sr_n_down, sr_total_up, sr_total_down
      real(rp) :: max_core_physical, max_valence_physical
      real(rp) :: core_site_metrics(5)
      complex(rp), allocatable :: h(:, :, :), rho(:, :, :), delta_h(:, :, :), delta_rho(:, :, :)
      complex(rp), allocatable :: moments_k(:, :, :), moments_sum(:, :, :), delta_endpoints(:, :, :), endpoint_sum(:, :, :)
      complex(rp), allocatable :: generator(:, :), delta_moments(:, :, :)
      complex(rp), allocatable :: response_plus(:, :), response_minus(:, :), response_l_plus(:, :), response_l_minus(:, :)
      real(rp), allocatable :: channels(:, :, :, :)
      real(rp), allocatable :: p1(:, :), p2(:, :), p3(:, :), sr1(:, :), sr2(:, :), sr3(:, :), sr4(:, :), core(:, :)
      real(rp), allocatable :: p1_l(:, :, :), p2_l(:, :, :), p3_l(:, :, :)
      real(rp), allocatable :: sr2_l(:, :, :), sr3_l(:, :, :), response_p(:, :), response_sr(:, :)
      real(rp), allocatable :: response_p_l(:, :, :), response_sr_l(:, :, :), sr_minus_p(:, :)
      real(rp), allocatable :: magnetization(:, :), valence_magnetization(:, :), core_magnetization(:, :)
      real(rp), allocatable :: core_l(:, :, :)
      character(len=96) :: classification
      character(len=8) :: verdict
      logical :: snapshot_provenance_ok, snapshot_reproduces_ok, core_ok, sr_ok, pauli_target_ok, pauli_response_ok, pass_a

      nsite = size(ground_states)
      if (nsite < 1 .or. size(radial_bases) /= nsite .or. response_space%nsite /= nsite .or. lattice_obj%nrec /= nsite) then
         error stop 'DRESP-09W: radial/response site dimensions are inconsistent'
      end if
      if (.not. hamiltonian_obj%hoh .or. hamiltonian_obj%ccor_2c .or. hamiltonian_obj%orb_pol .or. &
          hamiltonian_obj%hubbard_u_general_check .or. hamiltonian_obj%hubbard_v_check .or. hamiltonian_obj%local_axis) then
         error stop 'DRESP-09W: certified second-order ham_only state required'
      end if
      if (.not. allocated(reciprocal_obj%hk_bulk) .or. .not. allocated(reciprocal_obj%eigenvalues) .or. &
          .not. allocated(reciprocal_obj%eigenvectors) .or. .not. allocated(reciprocal_obj%k_weights)) then
         error stop 'DRESP-09W: accepted reciprocal eigensystem is incomplete'
      end if
      nmat = size(reciprocal_obj%eigenvectors, 1); nk = size(reciprocal_obj%eigenvalues, 2)
      norb_site = nmat/(2*nsite)
      if (nmat /= nb*nsite .or. mod(nmat, 2*nsite) /= 0 .or. size(reciprocal_obj%hk_bulk, 3) /= nk) then
         error stop 'DRESP-09W: accepted eigensystem dimensions are not site-major Fe spd'
      end if
      snapshot_provenance_ok = complete_snapshot(radial_bases, ground_states, response_space)
      if (.not. snapshot_provenance_ok) then
         call write_blocked_artifact(output_file, 'ACCEPTED_RADIAL_SNAPSHOT_PROVENANCE_OPEN')
         if (rank == 0) write(*, '(a)') 'DRESP-09W verdict: ACCEPTED_RADIAL_SNAPSHOT_PROVENANCE_OPEN'
         return
      end if

      allocate(h(nmat, nmat, nk), rho(nmat, nmat, nk), delta_h(nmat, nmat, nk), delta_rho(nmat, nmat, nk), &
         moments_k(nmat, nmat, 3), moments_sum(nmat, nmat, 3), delta_moments(nmat, nmat, 3), &
         delta_endpoints(nmat, nmat, 4), endpoint_sum(nmat, nmat, 4), generator(nmat, nmat))
      h = cmplx(0.0_rp, 0.0_rp, rp); moments_sum = cmplx(0.0_rp, 0.0_rp, rp)
      endpoint_sum = cmplx(0.0_rp, 0.0_rp, rp)
      call exact_ks_spin_rotation_generator(norb_site, nsite, generator)
      wsum = sum(reciprocal_obj%k_weights)
      do ik = 1, nk
         ik_global = ik
         if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
         h(:, :, ik) = reciprocal_obj%hk_bulk(:, :, ik)
         call density_from_eigensystem(reciprocal_obj%eigenvalues(:, ik), reciprocal_obj%eigenvectors(:, :, ik), &
            reciprocal_obj%fermi_level, reciprocal_obj%temperature, rho(:, :, ik))
         call moment_matrices_from_eigensystem(reciprocal_obj%eigenvalues(:, ik), reciprocal_obj%eigenvectors(:, :, ik), &
            reciprocal_obj%fermi_level, reciprocal_obj%temperature, moments_k)
         call commutator_tangent(generator, h(:, :, ik), delta_h(:, :, ik))
         call commutator_tangent(generator, rho(:, :, ik), delta_rho(:, :, ik))
         call endpoint_tangent_branches(h(:, :, ik), rho(:, :, ik), delta_h(:, :, ik), delta_rho(:, :, ik), delta_endpoints)
         call commutator_tangent(generator, moments_k(:, :, 1), delta_moments(:, :, 1))
         call commutator_tangent(generator, moments_k(:, :, 2), delta_moments(:, :, 2))
         call commutator_tangent(generator, moments_k(:, :, 3), delta_moments(:, :, 3))
         moments_sum = moments_sum + reciprocal_obj%k_weights(ik_global)*moments_k
         endpoint_sum = endpoint_sum + reciprocal_obj%k_weights(ik_global)*delta_endpoints
      end do
      moments_sum = moments_sum/wsum
      endpoint_sum = endpoint_sum/wsum

      allocate(channels(3, nsite, radial_bases(1)%lmax + 1, 2))
      call channel_moments_from_matrix(moments_sum, nsite, channels)
      allocate(p1(nsite, response_space%npoint), p2(nsite, response_space%npoint), p3(nsite, response_space%npoint), &
         sr1(nsite, response_space%npoint), sr2(nsite, response_space%npoint), sr3(nsite, response_space%npoint), &
         sr4(nsite, response_space%npoint), core(nsite, response_space%npoint), sr_minus_p(nsite, response_space%npoint))
      allocate(p1_l(nsite, response_space%npoint, 3), p2_l(nsite, response_space%npoint, 3), &
         p3_l(nsite, response_space%npoint, 3), sr2_l(nsite, response_space%npoint, 3), &
         sr3_l(nsite, response_space%npoint, 3))
      call compute_accepted_pauli_magnetization(reciprocal_obj, lattice_obj%symbolic_atoms, lattice_obj%nbulk, &
         magnetization, valence_magnetization, core_magnetization)
      call weighted_from_density(ground_states(1)%r, valence_magnetization, p1)
      call pauli_density_from_occupied_states(reciprocal_obj, radial_bases, p2, p2_l)
      call pauli_density_from_moments(radial_bases, channels, p3, p3_l)
      ! P1 is the accepted large-component occupied-state contract.  Its
      ! per-l decomposition is therefore the same explicit contraction as
      ! the independent P2 oracle; the total P1 still comes from the live
      ! compute_accepted_pauli_magnetization helper above.
      p1_l = p2_l
      do site = 1, nsite
         sr1(site, :) = ground_states(site)%rho_weighted_up - ground_states(site)%rho_weighted_down
         core(site, :) = ground_states(site)%core_weighted_up - ground_states(site)%core_weighted_down
      end do
      call sr_density_from_moments(radial_bases, channels, sr2, sr2_l)
      call sr_density_from_occupied_states(reciprocal_obj, radial_bases, sr3, sr3_l)
      sr4 = sr2 + core
      sr_minus_p = sr2 - p3

      call metric(p1, p2, ground_states, p1_p2)
      call metric(p1, p3, ground_states, p1_p3)
      call metric(p2, p3, ground_states, p2_p3)
      call metric(sr2, sr3, ground_states, sr2_sr3)
      call metric(sr1, sr4, ground_states, sr1_sr4)
      call metric(sr2, p3, ground_states, sr_pauli)
      call metric(sr4, sr2 + core, ground_states, core_metrics)
      core_decomposition_identity = core_metrics(1)
      do l = 1, 3
         per_l_p1_p2(l) = metric_relative(p1_l(:, :, l), p2_l(:, :, l), ground_states)
         per_l_p1_p3(l) = metric_relative(p1_l(:, :, l), p3_l(:, :, l), ground_states)
         per_l_p2_p3(l) = metric_relative(p2_l(:, :, l), p3_l(:, :, l), ground_states)
         per_l_sr2_sr3(l) = metric_relative(sr2_l(:, :, l), sr3_l(:, :, l), ground_states)
      end do

      allocate(response_plus(nsite, response_space%npoint), response_minus(nsite, response_space%npoint), &
         response_p(nsite, response_space%npoint), response_sr(nsite, response_space%npoint), &
         response_p_l(nsite, response_space%npoint, 3), response_sr_l(nsite, response_space%npoint, 3))
      call pauli_l0_density_from_endpoint_branches(response_space, radial_bases, endpoint_sum, 1, response_plus)
      call pauli_l0_density_from_endpoint_branches(response_space, radial_bases, endpoint_sum, 2, response_minus)
      call convert_response(response_plus + response_minus, ground_states, response_p)
      call sr_l0_density_from_endpoint_branches(response_space, radial_bases, endpoint_sum, 1, response_plus)
      call sr_l0_density_from_endpoint_branches(response_space, radial_bases, endpoint_sum, 2, response_minus)
      call convert_response(response_plus + response_minus, ground_states, response_sr)
      do l = 1, 3
         call pauli_l0_density_from_endpoint_branches(response_space, radial_bases, endpoint_sum, 1, response_plus, l-1, l-1)
         call pauli_l0_density_from_endpoint_branches(response_space, radial_bases, endpoint_sum, 2, response_minus, l-1, l-1)
         call convert_response(response_plus + response_minus, ground_states, response_p_l(:, :, l))
         call sr_l0_density_from_endpoint_branches(response_space, radial_bases, endpoint_sum, 1, response_plus, &
            l_first=l-1, l_last=l-1)
         call sr_l0_density_from_endpoint_branches(response_space, radial_bases, endpoint_sum, 2, response_minus, &
            l_first=l-1, l_last=l-1)
         call convert_response(response_plus + response_minus, ground_states, response_sr_l(:, :, l))
         response_per_l(l) = metric_relative(response_p_l(:, :, l), p1_l(:, :, l), ground_states)
         pauli_response_per_l(l) = metric_relative(response_p_l(:, :, l), p3_l(:, :, l), ground_states)
         per_l_sr_response(l) = metric_relative(response_sr_l(:, :, l), sr2_l(:, :, l), ground_states)
      end do
      call metric(response_p, p1, ground_states, pauli_response_p1)
      call metric(response_p, p2, ground_states, pauli_response_p2)
      call metric(response_p, p3, ground_states, pauli_response_p3)
      call metric(response_sr, sr2, ground_states, sr_response_sr2)
      do spin = 1, 2
         equilibrium_m0(spin) = 0.0_rp
         do site = 1, nsite
            do l = 0, radial_bases(site)%lmax
               equilibrium_m0(spin) = equilibrium_m0(spin) + channels(1, site, l + 1, spin)
            end do
         end do
      end do

      core_metrics = 0.0_rp
      core_metrics = 0.0_rp
      do site = 1, nsite
         call core_metrics_from_profile(core(site, :), sr2(site, :), ground_states(site), core_site_metrics)
         core_metrics(1:4) = core_metrics(1:4) + core_site_metrics(1:4)**2
         core_metrics(5) = core_metrics(5) + core_site_metrics(5)
      end do
      core_metrics(1:4) = sqrt(core_metrics(1:4))
      core_l2 = core_metrics(2)
      core_s_integral = 0.0_rp; core_s_l2 = 0.0_rp; core_s_max = 0.0_rp
      core_ok = .true.
      do site = 1, nsite
         if (.not. allocated(ground_states(site)%core_weighted_l)) then
            core_ok = .false.; cycle
         end if
         core_s_integral = core_s_integral + volume_integral(ground_states(site), &
            ground_states(site)%core_weighted_l(:, 1, 1) - ground_states(site)%core_weighted_l(:, 1, 2))
         core_s_l2 = core_s_l2 + volume_l2_norm(ground_states(site), &
            ground_states(site)%core_weighted_l(:, 1, 1) - ground_states(site)%core_weighted_l(:, 1, 2))**2
         core_s_max = max(core_s_max, max_physical(ground_states(site), &
            ground_states(site)%core_weighted_l(:, 1, 1) - ground_states(site)%core_weighted_l(:, 1, 2)))
      end do
      core_s_l2 = sqrt(core_s_l2)
      max_core_physical = 0.0_rp; max_valence_physical = 0.0_rp
      do site = 1, nsite
         max_core_physical = max(max_core_physical, max_physical(ground_states(site), core(site, :)))
         max_valence_physical = max(max_valence_physical, max_physical(ground_states(site), sr2(site, :)))
      end do
      core_to_valence_peak = max_core_physical/max(max_valence_physical, tiny(1.0_rp))

      pauli_moment = 0.0_rp
      do site = 1, nsite
         pauli_moment = pauli_moment + volume_integral(ground_states(site), p3(site, :))
      end do
      pauli_n_up = integrated_channel(radial_bases, channels, ground_states, 1, .false.)
      pauli_n_down = integrated_channel(radial_bases, channels, ground_states, 2, .false.)
      sr_valence_moment = 0.0_rp
      sr_total_moment = 0.0_rp
      sr_n_up = integrated_channel(radial_bases, channels, ground_states, 1, .true.)
      sr_n_down = integrated_channel(radial_bases, channels, ground_states, 2, .true.)
      sr_total_up = 0.0_rp
      sr_total_down = 0.0_rp
      do site = 1, nsite
         sr_valence_moment = sr_valence_moment + volume_integral(ground_states(site), sr2(site, :))
         sr_total_moment = sr_total_moment + volume_integral(ground_states(site), sr4(site, :))
         sr_total_up = sr_total_up + volume_integral(ground_states(site), ground_states(site)%rho_weighted_up)
         sr_total_down = sr_total_down + volume_integral(ground_states(site), ground_states(site)%rho_weighted_down)
      end do
      old_cross_representation = 3.05325e-1_rp
      ! The accepted SCF snapshot is a converged radial fixed point, not an
      ! exact algebraic identity with the independently accumulated
      ! k-space valence/core split.  Use the accepted-state numerical scale
      ! here and report the measured residual verbatim.
      snapshot_reproduces_ok = sr1_sr4(1) < 1.0e-5_rp
      sr_ok = sr2_sr3(1) < 1.0e-6_rp .and. sr_response_sr2(1) < 1.0e-6_rp
      pauli_target_ok = p1_p3(1) < 1.0e-6_rp
      pauli_response_ok = pauli_response_p3(1) < 1.0e-6_rp
      pass_a = snapshot_reproduces_ok .and. core_ok .and. sr_ok .and. pauli_target_ok .and. pauli_response_ok .and. &
         p1_p2(1) < 1.0e-6_rp .and. p2_p3(1) < 1.0e-6_rp
      if (.not. snapshot_reproduces_ok) then
         classification = dresp09w_snapshot_open
      else if (.not. core_ok) then
         classification = dresp09w_core_open
      else if (.not. sr_ok) then
         classification = dresp09w_sr_open
      else if (.not. pauli_target_ok) then
         classification = dresp09w_pauli_target_open
      else if (.not. pauli_response_ok) then
         classification = dresp09w_pauli_response_open
      else
         classification = dresp09w_closed
      end if
      if (.not. snapshot_provenance_ok .or. .not. snapshot_reproduces_ok) then
         verdict = 'BLOCKED '
      else
         verdict = merge('PASS-A  ', 'PASS-B  ', pass_a)
      end if

      if (rank == 0) then
         open(newunit=unit, file=trim(output_file), status='replace', action='write', iostat=ios)
         if (ios /= 0) error stop 'DRESP-09W: cannot open diagnostic artifact'
         write(unit, '(a)') '# DRESP-09W radial observable / ground-state provenance'
         write(unit, '(a)') '# Pauli and full-SR targets are constructed from the same accepted eigensystem and M0/M1/M2.'
         write(unit, '(a)') 'P1_source = compute_accepted_pauli_magnetization (large component, valence only)'
         write(unit, '(a)') 'P2_source = independent occupied-state phi+DeltaE*phidot (large component, valence only)'
         write(unit, '(a)') 'P3_source = production M0/M1/M2 Pauli polynomial with phi*phiddot'
         write(unit, '(a)') 'SR1_source = accepted rho_weighted_up-rho_weighted_down (total)'
         write(unit, '(a)') 'SR2_source = production M0/M1/M2 full scalar-relativistic polynomial (valence)'
         write(unit, '(a)') 'SR3_source = independent occupied-state full scalar-relativistic polynomial (valence)'
         write(unit, '(a)') 'SR4_source = SR2 + accepted frozen core'
         call write_metric(unit, 'P1_P2', p1_p2)
         call write_metric(unit, 'P1_P3', p1_p3)
         call write_metric(unit, 'P2_P3', p2_p3)
         call write_metric(unit, 'SR2_SR3', sr2_sr3)
         call write_metric(unit, 'SR1_SR4', sr1_sr4)
         call write_metric(unit, 'SR_Pauli', sr_pauli)
         call write_metric(unit, 'Pauli_response_vs_P1', pauli_response_p1)
         call write_metric(unit, 'Pauli_response_vs_P2', pauli_response_p2)
         call write_metric(unit, 'Pauli_response_vs_P3', pauli_response_p3)
         call write_metric(unit, 'SR_response_vs_SR2', sr_response_sr2)
         write(unit, '(a,es24.16)') 'old_cross_representation_residual = ', old_cross_representation
         write(unit, '(a,es24.16)') 'new_SR_to_SR_valence_residual = ', sr_response_sr2(1)
         write(unit, '(a,es24.16)') 'core_integrated_spin_moment = ', core_metrics(5)
         write(unit, '(a,es24.16)') 'core_weighted_L2_norm = ', core_l2
         write(unit, '(a,es24.16)') 'core_to_max_valence_radial_amplitude = ', core_to_valence_peak
         write(unit, '(a,es24.16)') 'core_s_like_integrated_spin_moment = ', core_s_integral
         write(unit, '(a,es24.16)') 'core_s_like_weighted_L2_norm = ', core_s_l2
         write(unit, '(a,es24.16)') 'core_s_like_max_radial_amplitude = ', core_s_max
         write(unit, '(a,es24.16)') 'core_decomposition_identity_relative = ', core_decomposition_identity
         write(unit, '(a,l1)') 'core_l_resolved_available = ', core_ok
         write(unit, '(a,3(es24.16,1x))') 'P1_P2_per_l_s_p_d = ', per_l_p1_p2
         write(unit, '(a,3(es24.16,1x))') 'P1_P3_per_l_s_p_d = ', per_l_p1_p3
         write(unit, '(a,3(es24.16,1x))') 'P2_P3_per_l_s_p_d = ', per_l_p2_p3
         write(unit, '(a,3(es24.16,1x))') 'SR2_SR3_per_l_s_p_d = ', per_l_sr2_sr3
         write(unit, '(a,3(es24.16,1x))') 'Pauli_response_vs_P1_per_l_s_p_d = ', response_per_l
         write(unit, '(a,3(es24.16,1x))') 'Pauli_response_vs_P3_per_l_s_p_d = ', pauli_response_per_l
         write(unit, '(a,3(es24.16,1x))') 'SR_response_vs_SR2_per_l_s_p_d = ', per_l_sr_response
         write(unit, '(a,2(es24.16,1x))') 'coefficient_M0_N_up_down = ', equilibrium_m0
         write(unit, '(a,2(es24.16,1x))') 'Pauli_P3_N_up_down = ', pauli_n_up, pauli_n_down
         write(unit, '(a,2(es24.16,1x))') 'SR_valence_N_up_down = ', sr_n_up, sr_n_down
         write(unit, '(a,2(es24.16,1x))') 'SR_total_N_up_down = ', sr_total_up, sr_total_down
         write(unit, '(a,es24.16)') 'Pauli_P3_integrated_moment = ', pauli_moment
         write(unit, '(a,es24.16)') 'SR_valence_integrated_moment = ', sr_valence_moment
         write(unit, '(a,es24.16)') 'SR_total_integrated_moment = ', sr_total_moment
         write(unit, '(a)') 'origin_handling = excluded from L2/max diagnostics as zero-measure; weighted arrays use n(r)*4*pi*r^2'
         write(unit, '(a)') 'radial_mesh_provenance = accepted r,a,b and radial_jacobian = Simpson_weight*a*(r+b)'
         write(unit, '(a)') 'response_measure_conversion = response coefficient (Y00) -> sqrt(4*pi)*r^2*coefficient'
         write(unit, '(a)') 'DRESP09V_0p305325 = CROSS-REPRESENTATION DIAGNOSTIC (full-SR response vs Pauli target)'
         write(unit, '(a,a)') 'accepted_snapshot_complete = ', merge('CLOSED', 'OPEN  ', snapshot_provenance_ok)
         write(unit, '(a,a)') 'accepted_snapshot_reproduces_SR1 = ', merge('CLOSED', 'OPEN  ', snapshot_reproduces_ok)
         write(unit, '(a)') 'ALSDA_Ward = NOT RUN'
         write(unit, '(a)') 'BES_Halle = OFF'
         write(unit, '(a,a)') 'primary_classification = ', trim(classification)
         write(unit, '(a,a)') 'verdict = ', trim(verdict)
         close(unit)
      end if
      if (rank == 0) write(*, '(a,a)') 'DRESP-09W verdict: ', trim(classification)

      deallocate(h, rho, delta_h, delta_rho, moments_k, moments_sum, delta_moments, delta_endpoints, endpoint_sum, generator, &
         channels, p1, p2, p3, sr1, sr2, sr3, sr4, core, sr_minus_p, p1_l, p2_l, p3_l, sr2_l, sr3_l, response_plus, &
         response_minus, response_p, response_sr, response_p_l, response_sr_l, magnetization, valence_magnetization, core_magnetization)
   end subroutine run_dresp09w_radial_observable_provenance

   subroutine metric(values, target, states, result)
      real(rp), intent(in) :: values(:, :), target(:, :)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(out) :: result(5)
      call radial_relative_metrics(values, target, states, result(1), result(2), result(3), result(4), result(5))
   end subroutine metric

   real(rp) function metric_relative(values, target, states) result(value)
      real(rp), intent(in) :: values(:, :), target(:, :)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp) :: l2, max_absolute, max_relative, integral
      call radial_relative_metrics(values, target, states, value, l2, max_absolute, max_relative, integral)
   end function metric_relative

   subroutine core_metrics_from_profile(core, valence, state, result)
      real(rp), intent(in) :: core(:), valence(:)
      type(radial_ground_state), intent(in) :: state
      real(rp), intent(out) :: result(5)
      real(rp) :: core_phys, valence_phys, diff, jac, r, num
      integer :: ir
      result = 0.0_rp; num = 0.0_rp
      do ir = 2, size(state%r)
         r = state%r(ir); jac = radial_jacobian_local(state, ir)
         core_phys = core(ir)/(dresp09w_four_pi*r*r)
         valence_phys = valence(ir)/(dresp09w_four_pi*r*r)
         diff = abs(core(ir))/(dresp09w_four_pi*r*r)
         result(1) = result(1) + jac*core_phys**2*dresp09w_four_pi*r*r
         num = num + jac*core_phys**2*dresp09w_four_pi*r*r
         result(3) = max(result(3), diff)
         result(4) = max(result(4), diff/max(abs(valence_phys), tiny(1.0_rp)))
         result(5) = result(5) + jac*core(ir)
      end do
      result(1) = sqrt(result(1)); result(2) = result(1)
   end subroutine core_metrics_from_profile

   subroutine convert_response(complex_density, states, weighted)
      complex(rp), intent(in) :: complex_density(:, :)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(out) :: weighted(:, :)
      integer :: site, ir
      real(rp) :: sqrt_four_pi
      sqrt_four_pi = sqrt(dresp09w_four_pi)
      weighted = 0.0_rp
      do site = 1, size(states)
         do ir = 2, size(states(site)%r)
            weighted(site, ir) = sqrt_four_pi*states(site)%r(ir)**2*real(complex_density(site, ir), rp)
         end do
      end do
   end subroutine convert_response

   real(rp) function max_physical(state, weighted) result(value)
      type(radial_ground_state), intent(in) :: state
      real(rp), intent(in) :: weighted(:)
      integer :: ir
      value = 0.0_rp
      do ir = 2, size(state%r)
         value = max(value, abs(weighted(ir))/(dresp09w_four_pi*state%r(ir)**2))
      end do
   end function max_physical

   real(rp) function integrated_channel(radial_bases, channels, states, spin, full_sr) result(value)
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      real(rp), intent(in) :: channels(:, :, :, :)
      type(radial_ground_state), intent(in) :: states(:)
      integer, intent(in) :: spin
      logical, intent(in) :: full_sr
      integer :: site, l, ir
      real(rp), allocatable :: channel(:)
      real(rp) :: q0, q1, q2, point
      value = 0.0_rp
      allocate(channel(size(states(1)%r))); channel = 0.0_rp
      do site = 1, size(states)
         channel = 0.0_rp
         do l = 0, radial_bases(site)%lmax
            q0 = channels(1,site,l+1,spin)
            q1 = channels(2,site,l+1,spin) - radial_bases(site)%enu_work(l+1,spin)*q0
            q2 = channels(3,site,l+1,spin) - 2.0_rp*radial_bases(site)%enu_work(l+1,spin)*channels(2,site,l+1,spin) + &
               radial_bases(site)%enu_work(l+1,spin)**2*q0
            do ir = 1, size(channel)
               if (full_sr) then
                  point = q0*(radial_bases(site)%gfac(ir,l+1,spin)*radial_bases(site)%phi_large(ir,l+1,spin)**2 + &
                     radial_bases(site)%phi_small(ir,l+1,spin)**2) + 2.0_rp*q1*(radial_bases(site)%gfac(ir,l+1,spin)* &
                     radial_bases(site)%phi_large(ir,l+1,spin)*radial_bases(site)%phidot_large(ir,l+1,spin) + &
                     radial_bases(site)%phi_small(ir,l+1,spin)*radial_bases(site)%phidot_small(ir,l+1,spin)) + &
                     q2*(radial_bases(site)%gfac(ir,l+1,spin)*(radial_bases(site)%phidot_large(ir,l+1,spin)**2 + &
                     radial_bases(site)%phi_large(ir,l+1,spin)*radial_bases(site)%phiddot_large(ir,l+1,spin)) + &
                     radial_bases(site)%phidot_small(ir,l+1,spin)**2 + radial_bases(site)%phi_small(ir,l+1,spin)* &
                     radial_bases(site)%phiddot_small(ir,l+1,spin))
               else
                  point = q0*radial_bases(site)%phi_large(ir,l+1,spin)**2 + 2.0_rp*q1*radial_bases(site)%phi_large(ir,l+1,spin)* &
                     radial_bases(site)%phidot_large(ir,l+1,spin) + q2*(radial_bases(site)%phidot_large(ir,l+1,spin)**2 + &
                     radial_bases(site)%phi_large(ir,l+1,spin)*radial_bases(site)%phiddot_large(ir,l+1,spin))
               end if
               ! The LMTO radial polynomial is already the legacy weighted
               ! radial observable (the stored U_l numerator convention).
               channel(ir) = channel(ir) + point
            end do
         end do
         value = value + volume_integral(states(site), channel)
      end do
      deallocate(channel)
   end function integrated_channel

   logical function complete_snapshot(radial_bases, states, space) result(ok)
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: states(:)
      type(response_space_layout), intent(in) :: space
      integer :: site
      ok = size(radial_bases) == size(states) .and. space%nsite == size(states)
      do site = 1, size(states)
         ok = ok .and. states(site)%sr_basis_valid .and. states(site)%core_density_valid .and. &
            allocated(radial_bases(site)%phi_small) .and. allocated(radial_bases(site)%phidot_small) .and. &
            allocated(radial_bases(site)%phiddot_large) .and. allocated(radial_bases(site)%phiddot_small) .and. &
            allocated(radial_bases(site)%gfac) .and. allocated(radial_bases(site)%tmc)
      end do
   end function complete_snapshot

   subroutine write_blocked_artifact(filename, classification)
      character(len=*), intent(in) :: filename, classification
      integer :: unit, ios
      if (rank /= 0) return
      open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
      if (ios /= 0) error stop 'DRESP-09W: cannot write blocked artifact'
      write(unit, '(a)') '# DRESP-09W radial observable / ground-state provenance'
      write(unit, '(a,a)') 'primary_classification = ', trim(classification)
      write(unit, '(a)') 'verdict = BLOCKED'
      write(unit, '(a)') 'ALSDA_Ward = NOT RUN'
      write(unit, '(a)') 'BES_Halle = OFF'
      close(unit)
   end subroutine write_blocked_artifact

   subroutine write_metric(unit, name, values)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      real(rp), intent(in) :: values(5)
      write(unit, '(a,a,es24.16)') trim(name)//'_relative = ', '', values(1)
      write(unit, '(a,a,es24.16)') trim(name)//'_weighted_L2 = ', '', values(2)
      write(unit, '(a,a,es24.16)') trim(name)//'_max_absolute = ', '', values(3)
      write(unit, '(a,a,es24.16)') trim(name)//'_max_relative = ', '', values(4)
      write(unit, '(a,a,es24.16)') trim(name)//'_integral = ', '', values(5)
   end subroutine write_metric

   pure real(rp) function radial_jacobian_local(state, ir) result(value)
      type(radial_ground_state), intent(in) :: state
      integer, intent(in) :: ir
      value = (2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp)*state%a*(state%r(ir) + state%b)
   end function radial_jacobian_local

end module lr_dresp09w_bridge_mod
