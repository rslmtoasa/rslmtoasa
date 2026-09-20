!------------------------------------------------------------------------------
! DRESP-09S accepted-state material gate.
!
! This bridge owns the independent L=0 scalar-relativistic fixed-basis field
! and density contractions.  DRESP-08 is consumed only as a material oracle;
! DRESP-09R compact/Pauli routes are deliberately not called here.
!------------------------------------------------------------------------------
module lr_dresp09s_bridge_mod

   use precision_mod, only: rp
   use mpi_mod, only: rank
   use basis_mod, only: nb
   use reciprocal_mod, only: reciprocal
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use radial_ground_state_mod, only: radial_ground_state
   use pauli_ground_state_projection_mod, only: compute_accepted_pauli_magnetization
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use response_basis_mapping_mod, only: response_super_index, response_unflatten_superindex
   use lr_response_space_mod, only: response_space_layout
   use lr_dresp09s_scalar_relativistic_mod, only: sr_l0_source_with_flags, &
      sr_l0_density_from_matrix_hamiltonian
   use lr_dresp08_native_mapping_mod, only: dresp08_build_product_tangent, dresp08_build_native_realspace_tangent, &
      dresp08_build_native_onsite_tangents, dresp08_block_diagonal, dresp08_commutator_tangent, &
      dresp08_extract_spin_blocks, dresp08_rotate_collinear_operator, dresp08_build_h2, dresp08_relative_residual
   use lr_exact_ks_ward_mod, only: exact_ks_spin_rotation_generator, exact_ks_global_rotation_oracle, &
      exact_ks_ward_metrics
   implicit none
   private

   character(len=*), parameter, public :: dresp09s_closed = 'SR_L0_ROTATION_AUGMENTATION_CLOSED'
   character(len=*), parameter, public :: dresp09s_field_open = 'FIELD_CLOSED_DENSITY_OPEN'
   character(len=*), parameter, public :: dresp09s_density_open = 'DENSITY_CLOSED_FIELD_OPEN'
   character(len=*), parameter, public :: dresp09s_basis_required = 'BASIS_ROTATION_TERM_REQUIRED'
   character(len=*), parameter, public :: dresp09s_provenance_blocked = 'LEGACY_SR_PROVENANCE_INSUFFICIENT'
   character(len=*), parameter, public :: dresp09s_derivation_failure = 'SR_DERIVATION_FAILURE'

   public :: run_dresp09s_scalar_relativistic

contains

   subroutine run_dresp09s_scalar_relativistic(output_file, response_space, radial_bases, ground_states, reciprocal_obj, &
                                               lattice_obj, hamiltonian_obj)
      character(len=*), intent(in) :: output_file
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj

      complex(rp), allocatable :: h(:, :, :), h2(:, :, :), o(:, :, :), enu(:, :, :)
      complex(rp), allocatable :: d_h(:, :, :), d_h2(:, :, :), d_o(:, :, :), d_enu(:, :, :), d_ee(:, :, :, :)
      complex(rp), allocatable :: o_types(:, :, :), e_types(:, :, :), d_o_types(:, :, :), d_e_types(:, :, :)
      complex(rp), allocatable :: field_raw(:), source_pauli_plus(:, :), source_pauli_minus(:, :)
      complex(rp), allocatable :: source_lower_plus(:, :), source_lower_minus(:, :), source_full_plus(:, :), source_full_minus(:, :)
      complex(rp), allocatable :: native(:, :), candidate(:, :), generator(:, :), diff(:, :)
      complex(rp), allocatable :: term_enu(:, :), term_h(:, :), term_left(:, :), term_middle(:, :), term_right(:, :)
      complex(rp), allocatable :: hp(:, :), hm(:, :), op(:, :), om(:, :), ep(:, :), em(:, :), h2p(:, :), h2m(:, :)
      complex(rp), allocatable :: hrotp(:, :), hrotm(:, :), orotp(:, :), orotm(:, :), enurotp(:, :), enurotm(:, :)
      complex(rp), allocatable :: density(:, :), density_matrix(:, :), delta_h_ks(:, :), delta_rho_rot(:, :), delta_rho_spec(:, :)
      complex(rp), allocatable :: density_plus(:, :), density_minus(:, :)
      real(rp), allocatable :: hierarchy(:, :), per_k(:), per_k_max_element(:), occupied_action(:), near_ef_action(:)
      real(rp), allocatable :: finite_rotation(:), decomposition(:, :)
      real(rp), allocatable :: magnetization(:, :), valence_magnetization(:, :), core_magnetization(:, :)
      real(rp), allocatable :: pauli_weighted(:, :), sr_weighted(:, :), total_weighted(:, :), target_valence(:, :), target_total(:, :)
      integer, allocatable :: site_types(:)
      type(exact_ks_ward_metrics) :: ward_metrics
      real(rp) :: field_norm, native_norm, weighted_rms, max_relative, max_element
      real(rp) :: occupied_max, near_ef_max, wsum, theta, theta_residual
      real(rp) :: duality_residual, density_pauli_relative, density_sr_valence_relative, density_total_relative
      real(rp) :: core_integral, total_integral, valence_integral, sr_integral, density_target_norm
      real(rp) :: numerator, denominator, sqrt_four_pi, jac, r, pair_value, native_commutator_max
      complex(rp) :: pairing_field, pairing_trace
      integer :: nsite, norb_site, norb, nmat, nk, ik, ir, site, flat, unit, ios
      integer :: noccupied, nnear, ntype, theta_index
      logical :: field_closed, density_closed, basis_response_required
      logical :: pass_gate
      character(len=96) :: classification
      character(len=16) :: verdict

      nsite = size(ground_states)
      if (nsite < 1 .or. size(radial_bases) /= nsite .or. response_space%nsite /= nsite) then
         error stop 'DRESP-09S: radial/response site dimensions are inconsistent'
      end if
      if (.not. allocated(reciprocal_obj%hk_bulk) .or. .not. allocated(reciprocal_obj%k_weights)) then
         error stop 'DRESP-09S: accepted reciprocal Hamiltonian is incomplete'
      end if
      if (lattice_obj%nrec /= nsite .or. .not. hamiltonian_obj%hoh) then
         error stop 'DRESP-09S: accepted bcc-Fe state is not the certified HOH path'
      end if
      if (hamiltonian_obj%ccor_2c .or. hamiltonian_obj%orb_pol .or. hamiltonian_obj%hubbard_u_general_check .or. &
          hamiltonian_obj%hubbard_v_check .or. hamiltonian_obj%local_axis) then
         error stop 'DRESP-09S: excluded additive Hamiltonian correction is enabled'
      end if
      norb_site = (radial_bases(1)%lmax + 1)**2
      norb = norb_site*nsite
      nmat = size(reciprocal_obj%hk_bulk, 1)
      nk = size(reciprocal_obj%hk_bulk, 3)
      if (nmat /= nb*nsite .or. nb /= 2*norb_site .or. nk /= 64 .or. &
          size(reciprocal_obj%k_weights) /= nk) then
         error stop 'DRESP-09S: accepted gate requires the complete 4x4x4 bcc-Fe eigensystem'
      end if
      if (.not. radial_provenance_complete(radial_bases)) then
         error stop 'DRESP-09S: accepted scalar-relativistic provenance is incomplete'
      end if

      allocate(field_raw(response_space%ndim))
      call build_l0_bxc_field(response_space, ground_states, field_raw)
      allocate(d_ee, source=hamiltonian_obj%ee)
      allocate(o_types(nb, nb, size(hamiltonian_obj%obarm, 3)), e_types(nb, nb, size(hamiltonian_obj%enim, 3)), &
         d_o_types(nb, nb, size(hamiltonian_obj%obarm, 3)), d_e_types(nb, nb, size(hamiltonian_obj%enim, 3)))
      o_types = hamiltonian_obj%obarm
      e_types = hamiltonian_obj%enim
      call dresp08_build_native_realspace_tangent(hamiltonian_obj, d_ee)
      call dresp08_build_native_onsite_tangents(hamiltonian_obj, d_o_types, d_e_types)
      site_types = lattice_obj%ib(1:nsite)
      allocate(generator(nmat, nmat))
      call exact_ks_spin_rotation_generator(norb_site, nsite, generator)
      allocate(h(nmat, nmat, nk), h2(nmat, nmat, nk), o(nmat, nmat, nk), enu(nmat, nmat, nk), &
         d_h(nmat, nmat, nk), d_h2(nmat, nmat, nk), d_o(nmat, nmat, nk), d_enu(nmat, nmat, nk))
      allocate(source_pauli_plus(nmat, nmat), source_pauli_minus(nmat, nmat), source_lower_plus(nmat, nmat), &
         source_lower_minus(nmat, nmat), source_full_plus(nmat, nmat), source_full_minus(nmat, nmat), native(nmat, nmat), &
         candidate(nmat, nmat), diff(nmat, nmat), term_enu(nmat, nmat), term_h(nmat, nmat), term_left(nmat, nmat), &
         term_middle(nmat, nmat), term_right(nmat, nmat), hp(norb, norb), hm(norb, norb), op(norb, norb), om(norb, norb), &
         ep(norb, norb), em(norb, norb), h2p(nmat, nmat), h2m(nmat, nmat))
      allocate(hrotp(nmat, nmat), hrotm(nmat, nmat), orotp(nmat, nmat), orotm(nmat, nmat), &
         enurotp(nmat, nmat), enurotm(nmat, nmat))
      allocate(hierarchy(3, nk), per_k(nk), per_k_max_element(nk), occupied_action(nk), near_ef_action(nk), &
         finite_rotation(3), decomposition(5, nk))

      hierarchy = 0.0_rp
      per_k = 0.0_rp
      per_k_max_element = 0.0_rp
      occupied_action = 0.0_rp
      near_ef_action = 0.0_rp
      decomposition = 0.0_rp
      finite_rotation = 0.0_rp
      duality_residual = 0.0_rp
      do ik = 1, nk
         call reciprocal_kpoint(reciprocal_obj, ik, hamiltonian_obj%ee, h(:, :, ik))
         call reciprocal_kpoint(reciprocal_obj, ik, d_ee, d_h(:, :, ik))
         call dresp08_block_diagonal(o_types, site_types, o(:, :, ik))
         call dresp08_block_diagonal(e_types, site_types, enu(:, :, ik))
         call dresp08_block_diagonal(d_o_types, site_types, d_o(:, :, ik))
         call dresp08_block_diagonal(d_e_types, site_types, d_enu(:, :, ik))
         call dresp08_build_product_tangent(d_enu(:, :, ik), d_h(:, :, ik), d_o(:, :, ik), h(:, :, ik), o(:, :, ik), &
            d_h2(:, :, ik), term_enu, term_h, term_left, term_middle, term_right)
         h2(:, :, ik) = reciprocal_obj%hk_bulk(:, :, ik)
         call dresp08_commutator_tangent(generator, h2(:, :, ik), native)

         call sr_l0_source_with_flags(response_space, radial_bases, field_raw, 1, h2(:, :, ik), .false., .false., source_pauli_plus)
         call sr_l0_source_with_flags(response_space, radial_bases, field_raw, 2, h2(:, :, ik), .false., .false., source_pauli_minus)
         candidate = source_pauli_plus + source_pauli_minus
         hierarchy(1, ik) = dresp08_relative_residual(candidate, native)

         call sr_l0_source_with_flags(response_space, radial_bases, field_raw, 1, h2(:, :, ik), .true., .false., source_lower_plus)
         call sr_l0_source_with_flags(response_space, radial_bases, field_raw, 2, h2(:, :, ik), .true., .false., source_lower_minus)
         candidate = source_lower_plus + source_lower_minus
         hierarchy(2, ik) = dresp08_relative_residual(candidate, native)

         call sr_l0_source_with_flags(response_space, radial_bases, field_raw, 1, h2(:, :, ik), .true., .true., source_full_plus)
         call sr_l0_source_with_flags(response_space, radial_bases, field_raw, 2, h2(:, :, ik), .true., .true., source_full_minus)
         candidate = source_full_plus + source_full_minus
         hierarchy(3, ik) = dresp08_relative_residual(candidate, native)
         diff = native - candidate
         per_k(ik) = hierarchy(3, ik)
         per_k_max_element(ik) = maxval(abs(diff))
         call action_metrics(diff, native, reciprocal_obj%eigenvectors(:, :, ik), reciprocal_obj%eigenvalues(:, ik), &
            reciprocal_obj%fermi_level, occupied_action(ik), near_ef_action(ik), noccupied, nnear)
         decomposition(1, ik) = dresp08_relative_residual(diff, term_enu)
         decomposition(2, ik) = dresp08_relative_residual(diff, term_h)
         decomposition(3, ik) = dresp08_relative_residual(diff, term_left)
         decomposition(4, ik) = dresp08_relative_residual(diff, term_middle)
         decomposition(5, ik) = dresp08_relative_residual(diff, term_right)

         ! Independent finite-rotation regression for the representation term.
         call dresp08_extract_spin_blocks(h(:, :, ik), norb, 1, hp, hm)
         call dresp08_extract_spin_blocks(o(:, :, ik), norb, 1, op, om)
         call dresp08_extract_spin_blocks(enu(:, :, ik), norb, 1, ep, em)
         do theta_index = 1, 3
            theta = 1.0e-2_rp/2.0_rp**real(theta_index - 1, rp)
            call dresp08_rotate_collinear_operator(hp, hm, theta, hrotp)
            call dresp08_rotate_collinear_operator(hp, hm, -theta, hrotm)
            call dresp08_rotate_collinear_operator(op, om, theta, orotp)
            call dresp08_rotate_collinear_operator(op, om, -theta, orotm)
            call dresp08_rotate_collinear_operator(ep, em, theta, enurotp)
            call dresp08_rotate_collinear_operator(ep, em, -theta, enurotm)
            call dresp08_build_h2(hrotp, orotp, enurotp, h2p)
            call dresp08_build_h2(hrotm, orotm, enurotm, h2m)
            finite_rotation(theta_index) = max(finite_rotation(theta_index), &
               dresp08_relative_residual((h2p - h2m)/(2.0_rp*theta), native))
         end do
      end do

      wsum = sum(reciprocal_obj%k_weights)
      numerator = sum(reciprocal_obj%k_weights*per_k**2)
      weighted_rms = sqrt(numerator/max(wsum, tiny(1.0_rp)))
      max_relative = maxval(per_k)
      max_element = maxval(per_k_max_element)
      occupied_max = maxval(occupied_action)
      near_ef_max = maxval(near_ef_action)
      field_norm = sqrt(sum(abs(field_raw)**2*response_space%metric_weights))
      native_norm = 0.0_rp
      native_commutator_max = 0.0_rp
      do ik = 1, nk
         native_norm = native_norm + reciprocal_obj%k_weights(ik)*sum(abs(h2(:, :, ik))**2)
         call dresp08_commutator_tangent(generator, h2(:, :, ik), diff)
         native_commutator_max = max(native_commutator_max, dresp08_relative_residual(d_h2(:, :, ik), diff))
      end do
      native_norm = sqrt(native_norm/max(wsum, tiny(1.0_rp)))

      ! Density oracle: exact coefficient-space rigid rotation versus the
      ! accepted valence and frozen-core radial magnetizations.
      call compute_accepted_pauli_magnetization(reciprocal_obj, lattice_obj%symbolic_atoms, lattice_obj%nbulk, &
         magnetization, valence_magnetization, core_magnetization)
      allocate(pauli_weighted(nsite, response_space%npoint), sr_weighted(nsite, response_space%npoint), &
         total_weighted(nsite, response_space%npoint), target_valence(nsite, response_space%npoint), &
         target_total(nsite, response_space%npoint), density_plus(nsite, response_space%npoint), &
         density_minus(nsite, response_space%npoint), density(nsite, response_space%npoint), &
         density_matrix(nmat, nmat), delta_h_ks(nmat, nmat), delta_rho_rot(nmat, nmat), delta_rho_spec(nmat, nmat))
      pauli_weighted = 0.0_rp
      sr_weighted = 0.0_rp
      sqrt_four_pi = sqrt(4.0_rp*acos(-1.0_rp))
      do site = 1, nsite
         do ir = 1, response_space%npoint
            r = ground_states(site)%r(ir)
            target_valence(site, ir) = 4.0_rp*acos(-1.0_rp)*r*r*valence_magnetization(site, ir)
            target_total(site, ir) = ground_states(site)%rho_weighted_up(ir) - ground_states(site)%rho_weighted_down(ir)
         end do
      end do
      do ik = 1, nk
         call exact_ks_global_rotation_oracle(h2(:, :, ik), reciprocal_obj%eigenvalues(:, ik), &
            reciprocal_obj%eigenvectors(:, :, ik), reciprocal_obj%fermi_level, reciprocal_obj%temperature, generator, &
            delta_h_ks, density_matrix, delta_rho_rot, delta_rho_spec, ward_metrics, norb_site, nsite)
         call sr_l0_density_from_matrix_hamiltonian(response_space, radial_bases, delta_rho_rot, 1, h2(:, :, ik), &
            density_plus, .false., .false.)
         call sr_l0_density_from_matrix_hamiltonian(response_space, radial_bases, delta_rho_rot, 2, h2(:, :, ik), &
            density_minus, .false., .false.)
         do site = 1, nsite
            do ir = 1, response_space%npoint
               r = ground_states(site)%r(ir)
               pauli_weighted(site, ir) = pauli_weighted(site, ir) + reciprocal_obj%k_weights(ik)*sqrt_four_pi*r*r* &
                  real(density_plus(site, ir) + density_minus(site, ir), rp)/max(wsum, tiny(1.0_rp))
            end do
         end do
         call sr_l0_density_from_matrix_hamiltonian(response_space, radial_bases, delta_rho_rot, 1, h2(:, :, ik), &
            density_plus, .true., .true.)
         call sr_l0_density_from_matrix_hamiltonian(response_space, radial_bases, delta_rho_rot, 2, h2(:, :, ik), &
            density_minus, .true., .true.)
         call sr_l0_source_with_flags(response_space, radial_bases, field_raw, 1, h2(:, :, ik), .true., .true., source_full_plus)
         call sr_l0_source_with_flags(response_space, radial_bases, field_raw, 2, h2(:, :, ik), .true., .true., source_full_minus)
         pairing_trace = sum(delta_rho_rot*transpose(source_full_plus + source_full_minus))
         pairing_field = cmplx(0.0_rp, 0.0_rp, rp)
         do site = 1, nsite
            do ir = 1, response_space%npoint
               flat = ((site - 1)*(response_space%response_lmax + 1)**2*response_space%npoint) + ir
               pairing_field = pairing_field + field_raw(flat)*response_space%radial_weights(ir)* &
                  (density_plus(site, ir) + density_minus(site, ir))
            end do
         end do
         duality_residual = max(duality_residual, abs(pairing_field - pairing_trace)/ &
            max(abs(pairing_trace), tiny(1.0_rp)))
         do site = 1, nsite
            do ir = 1, response_space%npoint
               r = ground_states(site)%r(ir)
               sr_weighted(site, ir) = sr_weighted(site, ir) + reciprocal_obj%k_weights(ik)*sqrt_four_pi*r*r* &
                  real(density_plus(site, ir) + density_minus(site, ir), rp)/max(wsum, tiny(1.0_rp))
            end do
         end do
      end do
      do site = 1, nsite
         do ir = 1, response_space%npoint
            r = ground_states(site)%r(ir)
            total_weighted(site, ir) = sr_weighted(site, ir) + 4.0_rp*acos(-1.0_rp)*r*r*core_magnetization(site, ir)
         end do
      end do
      call weighted_density_relative(pauli_weighted, target_valence, ground_states, density_pauli_relative)
      call weighted_density_relative(sr_weighted, target_valence, ground_states, density_sr_valence_relative)
      call weighted_density_relative(total_weighted, target_total, ground_states, density_total_relative)
      call weighted_integral(sr_weighted, ground_states, sr_integral)
      call weighted_integral(target_valence, ground_states, valence_integral)
      call weighted_integral(target_total, ground_states, total_integral)
      core_integral = 0.0_rp
      do site = 1, nsite
         do ir = 1, response_space%npoint
            jac = (2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp)*ground_states(site)%a* &
               (ground_states(site)%r(ir) + ground_states(site)%b)
            core_integral = core_integral + jac*4.0_rp*acos(-1.0_rp)*ground_states(site)%r(ir)**2* &
               core_magnetization(site, ir)
         end do
      end do

      field_closed = weighted_rms < 2.0e-10_rp .and. max_relative < 2.0e-10_rp
      density_closed = density_total_relative < 2.0e-8_rp
      basis_response_required = (.not. field_closed) .and. maxval(hierarchy(3, :)) < 2.0_rp .and. &
         finite_rotation(3) < 2.0e-6_rp .and. finite_rotation(1) > 3.0_rp*finite_rotation(2) .and. &
         finite_rotation(2) > 3.0_rp*finite_rotation(3)
      if (field_closed .and. density_closed) then
         classification = dresp09s_closed
      else if (basis_response_required) then
         classification = dresp09s_basis_required
      else if (field_closed) then
         classification = dresp09s_field_open
      else if (density_closed) then
         classification = dresp09s_density_open
      else
         classification = dresp09s_field_open
      end if
      pass_gate = (field_closed .and. density_closed) .or. &
         (basis_response_required .and. duality_residual < 2.0e-10_rp)
      if (field_closed .and. density_closed) then
         verdict = 'PASS-A'
      else if (basis_response_required .and. pass_gate) then
         verdict = 'PASS-B'
      else
         verdict = 'BLOCKED'
      end if

      if (rank == 0) then
         open(newunit=unit, file=trim(output_file), status='replace', action='write', iostat=ios)
         if (ios /= 0) error stop 'DRESP-09S: cannot open diagnostic artifact'
         sqrt_four_pi = sqrt(4.0_rp*acos(-1.0_rp))
         write(unit, '(a)') '# DRESP-09S scalar-relativistic L=0 rigid-rotation gate'
         write(unit, '(a)') '# SCF = accepted immutable bcc-Fe 4x4x4 state; no radial recomputation'
         write(unit, '(a)') '# field = sqrt(4*pi)*radial_ground_state%bxc_pauli, L=0, M=0'
         write(unit, '(a)') '# direct_route = independent raw L=0 SR bilinear; DRESP-09R and compact adjoint are not called'
         write(unit, '(a)') '# G1 = G_l = r*g_l; G2 = Phi_l = r*g_l''/(2*M*c); TMC = 2*M*c'
         write(unit, '(a,es24.16)') 'same_spin_lower_transverse_factor = ', -1.0_rp/3.0_rp
         write(unit, '(a,i0)') 'accepted_kpoints = ', nk
         write(unit, '(a,es24.16)') 'field_norm = ', field_norm
         write(unit, '(a,es24.16)') 'native_tangent_weighted_norm = ', native_norm
         write(unit, '(a,es24.16)') 'pauli_large_component_weighted_rms = ', weighted_mean(hierarchy(1, :), reciprocal_obj%k_weights)
         write(unit, '(a,es24.16)') 'lower_radial_component_weighted_rms = ', weighted_mean(hierarchy(2, :), reciprocal_obj%k_weights)
         write(unit, '(a,es24.16)') 'full_sr_angular_weighted_rms = ', weighted_rms
         write(unit, '(a,es24.16)') 'full_sr_max_relative_frobenius = ', max_relative
         write(unit, '(a,es24.16)') 'full_sr_max_element = ', max_element
         write(unit, '(a,es24.16)') 'full_sr_occupied_state_action_max = ', occupied_max
         write(unit, '(a,es24.16)') 'full_sr_near_ef_action_max = ', near_ef_max
         write(unit, '(a,a)') 'basis_response_term_required = ', merge('YES','NO ',basis_response_required)
         write(unit, '(a,es24.16)') 'basis_response_relative_weighted_rms = ', weighted_rms
         write(unit, '(a,es24.16)') 'basis_response_max_relative_frobenius = ', max_relative
         write(unit, '(a,es24.16)') 'finite_rotation_theta_1e-2_relative = ', finite_rotation(1)
         write(unit, '(a,es24.16)') 'finite_rotation_theta_5e-3_relative = ', finite_rotation(2)
         write(unit, '(a,es24.16)') 'finite_rotation_theta_2p5e-3_relative = ', finite_rotation(3)
         write(unit, '(a,es24.16)') 'dresp08_native_commutator_max_relative = ', native_commutator_max
         write(unit, '(a,es24.16)') 'delta_native_minus_full_vs_delta_Enu = ', maxval(decomposition(1, :))
         write(unit, '(a,es24.16)') 'delta_native_minus_full_vs_delta_h = ', maxval(decomposition(2, :))
         write(unit, '(a,es24.16)') 'delta_native_minus_full_vs_minus_delta_h_o_h_left = ', maxval(decomposition(3, :))
         write(unit, '(a,es24.16)') 'delta_native_minus_full_vs_minus_h_delta_o_h = ', maxval(decomposition(4, :))
         write(unit, '(a,es24.16)') 'delta_native_minus_full_vs_minus_h_o_delta_h = ', maxval(decomposition(5, :))
         write(unit, '(a,es24.16)') 'field_density_duality_max_relative = ', duality_residual
         write(unit, '(a,es24.16)') 'density_pauli_projected_valence_relative = ', density_pauli_relative
         write(unit, '(a,es24.16)') 'density_full_sr_valence_relative = ', density_sr_valence_relative
         write(unit, '(a,es24.16)') 'density_full_sr_total_relative = ', density_total_relative
         write(unit, '(a,es24.16)') 'm_val_sr_integral = ', sr_integral
         write(unit, '(a,es24.16)') 'm_val_pauli_accepted_integral = ', valence_integral
         write(unit, '(a,es24.16)') 'm_core_integral = ', core_integral
         write(unit, '(a,es24.16)') 'm_total_accepted_integral = ', total_integral
         do ik = 1, nk
            write(unit, '(a,i0,a,3(es24.16,1x))') 'k_', ik, '_field_hierarchy = ', hierarchy(:, ik)
            write(unit, '(a,i0,a,4(es24.16,1x))') 'k_', ik, '_field_metrics = ', per_k(ik), per_k_max_element(ik), &
               occupied_action(ik), near_ef_action(ik)
         end do
         write(unit, '(a,a)') 'classification = ', trim(classification)
         write(unit, '(a,a)') 'verdict = ', trim(verdict)
         write(unit, '(a)') 'full_arbitrary_l_sr_vertex = NOT YET CERTIFIED'
         write(unit, '(a)') 'BES_Halle = OFF'
         write(unit, '(a)') 'ALSDA = OFF'
         close(unit)
      end if
      if (rank == 0) write(*, '(a,a)') 'DRESP-09S verdict: ', trim(classification)

      deallocate(field_raw, d_ee, o_types, e_types, d_o_types, d_e_types, site_types, generator, h, h2, o, enu, d_h, d_h2, &
         d_o, d_enu, source_pauli_plus, source_pauli_minus, source_lower_plus, source_lower_minus, source_full_plus, &
         source_full_minus, native, candidate, diff, term_enu, term_h, term_left, term_middle, term_right, hp, hm, op, om, &
         ep, em, h2p, h2m, hrotp, hrotm, orotp, orotm, enurotp, enurotm, hierarchy, per_k, per_k_max_element, occupied_action, near_ef_action, finite_rotation, &
         decomposition, magnetization, valence_magnetization, core_magnetization, pauli_weighted, sr_weighted, total_weighted, &
         target_valence, target_total, density_plus, density_minus, density, density_matrix, delta_h_ks, delta_rho_rot, delta_rho_spec)
   end subroutine run_dresp09s_scalar_relativistic

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

      if (size(field) /= space%ndim .or. size(states) /= space%nsite) error stop 'DRESP-09S: field shape mismatch'
      sqrt_four_pi = sqrt(4.0_rp*acos(-1.0_rp))
      field = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, 1, item)
         if (item%response_l == 0 .and. item%response_m == 0) then
            field(flat) = cmplx(sqrt_four_pi*states(item%site)%bxc_pauli(item%radial_point), 0.0_rp, rp)
         end if
      end do
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

   real(rp) function weighted_mean(values, weights) result(value)
      real(rp), intent(in) :: values(:), weights(:)
      value = sqrt(sum(weights*values**2)/max(sum(weights), tiny(1.0_rp)))
   end function weighted_mean

   subroutine weighted_density_relative(values, target, states, relative)
      real(rp), intent(in) :: values(:, :), target(:, :)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(out) :: relative
      integer :: site, ir
      real(rp) :: num, den, r, pi, jac, diff_phys, target_phys

      pi = acos(-1.0_rp)
      num = 0.0_rp
      den = 0.0_rp
      do site = 1, size(states)
         do ir = 2, size(states(site)%r)
            r = states(site)%r(ir)
            jac = (2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp)*states(site)%a*(r + states(site)%b)
            diff_phys = (values(site, ir) - target(site, ir))/(4.0_rp*pi*r*r)
            target_phys = target(site, ir)/(4.0_rp*pi*r*r)
            num = num + jac*4.0_rp*pi*r*r*diff_phys**2
            den = den + jac*4.0_rp*pi*r*r*target_phys**2
         end do
      end do
      relative = sqrt(num/max(den, tiny(1.0_rp)))
   end subroutine weighted_density_relative

   subroutine weighted_integral(values, states, result)
      real(rp), intent(in) :: values(:, :)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(out) :: result
      integer :: site, ir
      result = 0.0_rp
      do site = 1, size(states)
         do ir = 1, size(states(site)%r)
            result = result + (2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp)*states(site)%a* &
               (states(site)%r(ir) + states(site)%b)*values(site, ir)
         end do
      end do
   end subroutine weighted_integral

end module lr_dresp09s_bridge_mod
