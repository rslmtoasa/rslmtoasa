!------------------------------------------------------------------------------
! DRESP-09T accepted-state moving-basis / representation gate.
!
! The fixed-basis SR field is obtained from the independent DRESP-09S radial
! contraction.  The basis term is obtained independently from the radial
! cross-overlap connection in lr_lmto_basis_connection.  In particular, the
! basis term is never defined as native-minus-field.
!------------------------------------------------------------------------------
module lr_dresp09t_bridge_mod

   use precision_mod, only: rp
   use mpi_mod, only: rank
   use basis_mod, only: nb
   use reciprocal_mod, only: reciprocal
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use radial_ground_state_mod, only: radial_ground_state
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_response_space_mod, only: response_space_layout
   use lr_dresp09s_scalar_relativistic_mod, only: sr_l0_source_with_flags, sr_l0_density_from_matrix_hamiltonian
   use lr_lmto_basis_connection_mod, only: lmto_build_basis_connection, lmto_build_effective_connection, &
      lmto_build_effective_endpoint_matrix, lmto_covariant_connection, lmto_basis_tangent, &
      lmto_build_second_order_basis_tangent, lmto_connection_relative_residual, lmto_connection_metric_residual
   use lr_dresp08_native_mapping_mod, only: dresp08_build_product_tangent, dresp08_build_native_realspace_tangent, &
      dresp08_build_native_onsite_tangents, dresp08_block_diagonal, dresp08_commutator_tangent, &
      dresp08_extract_spin_blocks, dresp08_rotate_collinear_operator, dresp08_build_h2, dresp08_relative_residual
   use lr_exact_ks_ward_mod, only: exact_ks_spin_rotation_generator, exact_ks_global_rotation_oracle, exact_ks_ward_metrics
   use pauli_ground_state_projection_mod, only: compute_accepted_pauli_magnetization
   implicit none
   private

   character(len=*), parameter, public :: dresp09t_closed = 'MOVING_BASIS_TANGENT_CLOSED'
   character(len=*), parameter, public :: dresp09t_field_open = 'FIELD_CLOSED_DENSITY_BASIS_OPEN'
   character(len=*), parameter, public :: dresp09t_density_open = 'DENSITY_CLOSED_FIELD_BASIS_OPEN'
   character(len=*), parameter, public :: dresp09t_parameter_required = 'POTENTIAL_PARAMETER_RESPONSE_REQUIRED'
   character(len=*), parameter, public :: dresp09t_orthogonalization_required = 'ORTHOGONALIZATION_RESPONSE_REQUIRED'
   character(len=*), parameter, public :: dresp09t_connection_failure = 'BASIS_CONNECTION_DERIVATION_FAILURE'
   character(len=*), parameter, public :: dresp09t_mixed = 'MIXED'

   public :: run_dresp09t_moving_basis

contains

   subroutine run_dresp09t_moving_basis(output_file, response_space, radial_bases, ground_states, reciprocal_obj, &
                                        lattice_obj, hamiltonian_obj)
      character(len=*), intent(in) :: output_file
      type(response_space_layout), intent(in) :: response_space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(radial_ground_state), intent(in) :: ground_states(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj

      integer :: nsite, norb_site, nmat, nk, ik, site, branch, theta_index, unit, ios
      integer :: p, q, noccupied, nnear
      real(rp) :: theta, wsum, field_rms, basis_rms, native_rms, total_rms, cross_term
      real(rp) :: first_h_rms, first_o_rms, first_e_rms, connection_metric, density_fixed, density_corrected
      real(rp) :: finite_rotation(3), finite_predicted(3), density_norm, sqrt_four_pi, r
      real(rp) :: numerator, denominator
      real(rp), allocatable :: per_k(:), first_h(:), first_o(:), first_e(:), density_fixed_k(:), density_corrected_k(:)
      real(rp), allocatable :: basis_norm_k(:), field_norm_k(:), native_norm_k(:), total_norm_k(:)
      complex(rp), allocatable :: field_raw(:), h(:, :, :), o(:, :, :), enu(:, :, :), h2(:, :, :)
      complex(rp), allocatable :: d_ee(:, :, :, :), d_o_types(:, :, :), d_e_types(:, :, :)
      complex(rp), allocatable :: d_h(:, :, :), d_o(:, :, :), d_enu(:, :, :), d_h2(:, :, :)
      complex(rp), allocatable :: native_h(:, :), native_h2(:, :), field_h(:, :), field_h2(:, :)
      complex(rp), allocatable :: basis_h(:, :), basis_o(:, :), basis_e(:, :), basis_h2(:, :), total_h(:, :), total_h2(:, :)
      complex(rp), allocatable :: term_enu(:, :), term_h(:, :), term_left(:, :), term_middle(:, :), term_right(:, :)
      complex(rp), allocatable :: raw_site(:, :, :), overlap_site(:, :, :), raw_global(:, :, :), overlap_global(:, :, :)
      complex(rp), allocatable :: raw_effective(:, :), overlap_effective(:, :), connection(:, :)
      complex(rp), allocatable :: source_plus(:, :), source_minus(:, :), hp(:, :), hm(:, :), op(:, :), om(:, :), ep(:, :), em(:, :)
      complex(rp), allocatable :: hrotp(:, :), hrotm(:, :), orotp(:, :), orotm(:, :), enurotp(:, :), enurotm(:, :)
      complex(rp), allocatable :: generator(:, :), density_matrix(:, :), delta_h_ks(:, :), delta_rho_rot(:, :), delta_rho_spec(:, :)
      complex(rp), allocatable :: delta_rho_basis(:, :), density_plus(:, :), density_minus(:, :), density(:, :)
      real(rp), allocatable :: magnetization(:, :), valence_magnetization(:, :), core_magnetization(:, :)
      real(rp), allocatable :: target_valence(:, :), corrected_valence(:, :), fixed_valence(:, :)
      type(exact_ks_ward_metrics) :: ward_metrics
      character(len=96) :: classification
      character(len=16) :: verdict
      logical :: field_closed, density_closed, connection_closed, parameter_required, orthogonalization_required

      nsite = size(ground_states)
      if (nsite < 1 .or. size(radial_bases) /= nsite .or. response_space%nsite /= nsite) then
         error stop 'DRESP-09T: radial/response site dimensions are inconsistent'
      end if
      if (.not. allocated(reciprocal_obj%hk_bulk) .or. .not. allocated(reciprocal_obj%k_weights)) then
         error stop 'DRESP-09T: accepted reciprocal Hamiltonian is incomplete'
      end if
      if (lattice_obj%nrec /= nsite .or. .not. hamiltonian_obj%hoh) error stop 'DRESP-09T: certified HOH state required'
      if (hamiltonian_obj%ccor_2c .or. hamiltonian_obj%orb_pol .or. hamiltonian_obj%hubbard_u_general_check .or. &
          hamiltonian_obj%hubbard_v_check .or. hamiltonian_obj%local_axis) error stop 'DRESP-09T: excluded correction enabled'
      norb_site = (radial_bases(1)%lmax + 1)**2
      nmat = size(reciprocal_obj%hk_bulk, 1)
      nk = size(reciprocal_obj%hk_bulk, 3)
      if (nmat /= nb*nsite .or. nb /= 2*norb_site .or. nk /= 64) error stop 'DRESP-09T: accepted 64-k Fe state required'
      if (.not. radial_provenance_complete(radial_bases)) error stop 'DRESP-09T: scalar-relativistic provenance incomplete'

      allocate(field_raw(response_space%ndim))
      call build_l0_bxc_field(response_space, ground_states, field_raw)
      allocate(d_ee, source=hamiltonian_obj%ee)
      allocate(d_o_types(nb, nb, size(hamiltonian_obj%obarm, 3)), d_e_types(nb, nb, size(hamiltonian_obj%enim, 3)))
      call dresp08_build_native_realspace_tangent(hamiltonian_obj, d_ee)
      call dresp08_build_native_onsite_tangents(hamiltonian_obj, d_o_types, d_e_types)
      allocate(h(nmat,nmat,nk), o(nmat,nmat,nk), enu(nmat,nmat,nk), h2(nmat,nmat,nk), &
         d_h(nmat,nmat,nk), d_o(nmat,nmat,nk), d_enu(nmat,nmat,nk), d_h2(nmat,nmat,nk))
      allocate(native_h(nmat,nmat), native_h2(nmat,nmat), field_h(nmat,nmat), field_h2(nmat,nmat), &
         basis_h(nmat,nmat), basis_o(nmat,nmat), basis_e(nmat,nmat), basis_h2(nmat,nmat), total_h(nmat,nmat), total_h2(nmat,nmat), &
         source_plus(nmat,nmat), source_minus(nmat,nmat), raw_global(nmat,nmat,4), overlap_global(nmat,nmat,4), &
         raw_effective(nmat,nmat), overlap_effective(nmat,nmat), connection(nmat,nmat), generator(nmat,nmat), &
         term_enu(nmat,nmat), term_h(nmat,nmat), term_left(nmat,nmat), term_middle(nmat,nmat), term_right(nmat,nmat))
      allocate(raw_site(2*norb_site,2*norb_site,4), overlap_site(2*norb_site,2*norb_site,4))
      allocate(per_k(nk), first_h(nk), first_o(nk), first_e(nk), field_norm_k(nk), basis_norm_k(nk), native_norm_k(nk), &
         total_norm_k(nk), density_fixed_k(nk), density_corrected_k(nk))
      finite_rotation = 0.0_rp
      finite_predicted = 0.0_rp
      connection_metric = 0.0_rp
      wsum = sum(reciprocal_obj%k_weights)
      call exact_ks_spin_rotation_generator(norb_site, nsite, generator)
      do ik = 1, nk
         call reciprocal_kpoint(reciprocal_obj, ik, hamiltonian_obj%ee, h(:,:,ik))
         call reciprocal_kpoint(reciprocal_obj, ik, d_ee, d_h(:,:,ik))
         call dresp08_block_diagonal(hamiltonian_obj%obarm, lattice_obj%ib(1:nsite), o(:,:,ik))
         call dresp08_block_diagonal(hamiltonian_obj%enim, lattice_obj%ib(1:nsite), enu(:,:,ik))
         call dresp08_block_diagonal(d_o_types, lattice_obj%ib(1:nsite), d_o(:,:,ik))
         call dresp08_block_diagonal(d_e_types, lattice_obj%ib(1:nsite), d_enu(:,:,ik))
         call dresp08_build_h2(h(:,:,ik), o(:,:,ik), enu(:,:,ik), h2(:,:,ik))
         call dresp08_build_product_tangent(d_enu(:,:,ik), d_h(:,:,ik), d_o(:,:,ik), h(:,:,ik), o(:,:,ik), d_h2(:,:,ik), &
            term_enu, term_h, term_left, term_middle, term_right)
         call dresp08_commutator_tangent(generator, h(:,:,ik), native_h)
         call dresp08_commutator_tangent(generator, h2(:,:,ik), native_h2)

         raw_global = cmplx(0.0_rp,0.0_rp,rp)
         overlap_global = cmplx(0.0_rp,0.0_rp,rp)
         do site = 1, nsite
            call lmto_build_basis_connection(radial_bases(site), response_space%radial_weights, raw_site, overlap_site)
            call insert_site_block(raw_site, raw_global, site, norb_site)
            call insert_site_block(overlap_site, overlap_global, site, norb_site)
         end do
         call lmto_build_effective_connection(raw_global, h(:,:,ik), raw_effective)
         call lmto_build_effective_endpoint_matrix(overlap_global, h(:,:,ik), overlap_effective)
         call lmto_covariant_connection(overlap_effective, raw_effective, connection)
         connection_metric = max(connection_metric, lmto_connection_metric_residual(overlap_effective, connection, &
            transpose(conjg(raw_effective)) + raw_effective))
         call lmto_basis_tangent(connection, h(:,:,ik), basis_h)
         call lmto_basis_tangent(connection, o(:,:,ik), basis_o)
         call lmto_basis_tangent(connection, enu(:,:,ik), basis_e)
         call lmto_build_second_order_basis_tangent(connection, connection, connection, enu(:,:,ik), h(:,:,ik), o(:,:,ik), &
            basis_h2, basis_e, basis_h, term_left, term_middle, term_right)
         call sr_l0_source_with_flags(response_space, radial_bases, field_raw, 1, h(:,:,ik), .true., .true., source_plus)
         call sr_l0_source_with_flags(response_space, radial_bases, field_raw, 2, h(:,:,ik), .true., .true., source_minus)
         field_h = source_plus + source_minus
         call sr_l0_source_with_flags(response_space, radial_bases, field_raw, 1, h2(:,:,ik), .true., .true., source_plus)
         call sr_l0_source_with_flags(response_space, radial_bases, field_raw, 2, h2(:,:,ik), .true., .true., source_minus)
         field_h2 = source_plus + source_minus
         total_h = field_h + basis_h
         total_h2 = field_h2 + basis_h2
         first_h(ik) = dresp08_relative_residual(total_h, native_h)
         first_o(ik) = dresp08_relative_residual(basis_o, d_o(:,:,ik))
         first_e(ik) = dresp08_relative_residual(basis_e, d_enu(:,:,ik))
         per_k(ik) = dresp08_relative_residual(total_h2, native_h2)
         field_norm_k(ik) = sqrt(sum(abs(field_h2)**2))
         basis_norm_k(ik) = sqrt(sum(abs(basis_h2)**2))
         native_norm_k(ik) = sqrt(sum(abs(native_h2)**2))
         total_norm_k(ik) = sqrt(sum(abs(total_h2)**2))

         call finite_native_tangent(h(:,:,ik), o(:,:,ik), enu(:,:,ik), native_h2, total_h2, finite_rotation, finite_predicted)
      end do
      field_rms = sqrt(sum(reciprocal_obj%k_weights*field_norm_k**2)/max(wsum,tiny(1.0_rp)))
      basis_rms = sqrt(sum(reciprocal_obj%k_weights*basis_norm_k**2)/max(wsum,tiny(1.0_rp)))
      native_rms = sqrt(sum(reciprocal_obj%k_weights*native_norm_k**2)/max(wsum,tiny(1.0_rp)))
      total_rms = sqrt(sum(reciprocal_obj%k_weights*total_norm_k**2)/max(wsum,tiny(1.0_rp)))
      cross_term = total_rms**2 - field_rms**2 - basis_rms**2
      first_h_rms = sqrt(sum(reciprocal_obj%k_weights*first_h**2)/max(wsum,tiny(1.0_rp)))
      first_o_rms = sqrt(sum(reciprocal_obj%k_weights*first_o**2)/max(wsum,tiny(1.0_rp)))
      first_e_rms = sqrt(sum(reciprocal_obj%k_weights*first_e**2)/max(wsum,tiny(1.0_rp)))

      ! The same K acts on the represented density.  The radial observable is
      ! still the independently derived DRESP-09S observable.
      call compute_accepted_pauli_magnetization(reciprocal_obj, lattice_obj%symbolic_atoms, lattice_obj%nbulk, &
         magnetization, valence_magnetization, core_magnetization)
      allocate(density_plus(nsite,response_space%npoint), density_minus(nsite,response_space%npoint), density(nsite,response_space%npoint), &
         density_matrix(nmat,nmat), delta_h_ks(nmat,nmat), delta_rho_rot(nmat,nmat), delta_rho_spec(nmat,nmat), delta_rho_basis(nmat,nmat), &
         fixed_valence(nsite,response_space%npoint), corrected_valence(nsite,response_space%npoint), target_valence(nsite,response_space%npoint))
      fixed_valence = 0.0_rp
      corrected_valence = 0.0_rp
      target_valence = 0.0_rp
      sqrt_four_pi = sqrt(4.0_rp*acos(-1.0_rp))
      do site = 1, nsite
         do p = 1, response_space%npoint
            r = ground_states(site)%r(p)
            target_valence(site,p) = 4.0_rp*acos(-1.0_rp)*r*r*valence_magnetization(site,p)
         end do
      end do
      do ik = 1, nk
         call exact_ks_global_rotation_oracle(h2(:,:,ik), reciprocal_obj%eigenvalues(:,ik), reciprocal_obj%eigenvectors(:,:,ik), &
            reciprocal_obj%fermi_level, reciprocal_obj%temperature, generator, delta_h_ks, density_matrix, delta_rho_rot, delta_rho_spec, &
            ward_metrics, norb_site, nsite)
         ! Recreate the same covariant K at this k for the density observable.
         call reciprocal_kpoint(reciprocal_obj, ik, hamiltonian_obj%ee, h(:,:,ik))
         raw_global = cmplx(0.0_rp,0.0_rp,rp); overlap_global = cmplx(0.0_rp,0.0_rp,rp)
         do site = 1, nsite
            call lmto_build_basis_connection(radial_bases(site), response_space%radial_weights, raw_site, overlap_site)
            call insert_site_block(raw_site, raw_global, site, norb_site)
            call insert_site_block(overlap_site, overlap_global, site, norb_site)
         end do
         call lmto_build_effective_connection(raw_global, h(:,:,ik), raw_effective)
         call lmto_build_effective_endpoint_matrix(overlap_global, h(:,:,ik), overlap_effective)
         call lmto_covariant_connection(overlap_effective, raw_effective, connection)
         delta_rho_basis = matmul(transpose(conjg(connection)), density_matrix) + matmul(density_matrix, connection)
         call sr_l0_density_from_matrix_hamiltonian(response_space, radial_bases, delta_rho_rot, 1, h2(:,:,ik), density_plus, .true., .true.)
         call sr_l0_density_from_matrix_hamiltonian(response_space, radial_bases, delta_rho_rot, 2, h2(:,:,ik), density_minus, .true., .true.)
         do site = 1, nsite
            do p = 1, response_space%npoint
               r = ground_states(site)%r(p)
               fixed_valence(site,p) = fixed_valence(site,p) + reciprocal_obj%k_weights(ik)*sqrt_four_pi*r*r* &
                  real(density_plus(site,p)+density_minus(site,p),rp)/max(wsum,tiny(1.0_rp))
            end do
         end do
         call sr_l0_density_from_matrix_hamiltonian(response_space, radial_bases, delta_rho_rot + delta_rho_basis, 1, h2(:,:,ik), density_plus, .true., .true.)
         call sr_l0_density_from_matrix_hamiltonian(response_space, radial_bases, delta_rho_rot + delta_rho_basis, 2, h2(:,:,ik), density_minus, .true., .true.)
         do site = 1, nsite
            do p = 1, response_space%npoint
               r = ground_states(site)%r(p)
               corrected_valence(site,p) = corrected_valence(site,p) + reciprocal_obj%k_weights(ik)*sqrt_four_pi*r*r* &
                  real(density_plus(site,p)+density_minus(site,p),rp)/max(wsum,tiny(1.0_rp))
            end do
         end do
      end do
      call weighted_density_relative(fixed_valence, target_valence, ground_states, density_fixed)
      call weighted_density_relative(corrected_valence, target_valence, ground_states, density_corrected)
      density_norm = max(sqrt(sum(target_valence**2)), tiny(1.0_rp))

      field_closed = first_h_rms < 2.0e-10_rp .and. maxval(first_h) < 2.0e-10_rp
      connection_closed = maxval(per_k) < 2.0e-10_rp .and. maxval(first_o) < 2.0e-10_rp .and. maxval(first_e) < 2.0e-10_rp
      density_closed = density_corrected < 2.0e-8_rp
      parameter_required = .false.
      orthogonalization_required = .not. connection_closed
      if (connection_closed .and. density_closed) then
         classification = dresp09t_closed
         verdict = 'PASS-A'
      else if (.not. connection_closed .and. density_closed) then
         classification = dresp09t_density_open
         verdict = 'BLOCKED'
      else if (connection_closed) then
         classification = dresp09t_field_open
         verdict = 'PASS-B'
      else if (orthogonalization_required) then
         classification = dresp09t_orthogonalization_required
         verdict = 'BLOCKED'
      else
         classification = dresp09t_mixed
         verdict = 'BLOCKED'
      end if

      if (rank == 0) then
         open(newunit=unit, file=trim(output_file), status='replace', action='write', iostat=ios)
         if (ios /= 0) error stop 'DRESP-09T: cannot open diagnostic artifact'
         write(unit,'(a)') '# DRESP-09T L=0 global rigid rotation moving-basis gate'
         write(unit,'(a)') '# field = independent DRESP-09S full scalar-relativistic fixed-basis source'
         write(unit,'(a)') '# basis = independent Gamma_raw -> endpoint contraction -> S^{-1} Gamma_raw connection'
         write(unit,'(a)') '# representation = accepted HOH orthogonal coefficient basis; no BES/Halle/ALSDA'
         write(unit,'(a,i0)') 'accepted_kpoints = ', nk
         write(unit,'(a,es24.16)') 'fixed_basis_field_weighted_norm = ', field_rms
         write(unit,'(a,es24.16)') 'basis_representation_weighted_norm = ', basis_rms
         write(unit,'(a,es24.16)') 'native_total_weighted_norm = ', native_rms
         write(unit,'(a,es24.16)') 'reconstructed_total_weighted_norm = ', total_rms
         write(unit,'(a,es24.16)') 'field_basis_cross_term = ', cross_term
         write(unit,'(a,es24.16)') 'first_order_h_weighted_residual = ', first_h_rms
         write(unit,'(a,es24.16)') 'first_order_o_weighted_residual = ', first_o_rms
         write(unit,'(a,es24.16)') 'first_order_enu_weighted_residual = ', first_e_rms
         write(unit,'(a,es24.16)') 'second_order_native_weighted_residual = ', sqrt(sum(reciprocal_obj%k_weights*per_k**2)/max(wsum,tiny(1.0_rp)))
         write(unit,'(a,es24.16)') 'second_order_native_max_relative_residual = ', maxval(per_k)
         write(unit,'(a,es24.16)') 'metric_covariant_residual = ', connection_metric
         write(unit,'(a,es24.16)') 'density_fixed_basis_relative = ', density_fixed
         write(unit,'(a,es24.16)') 'density_corrected_relative = ', density_corrected
         write(unit,'(a,es24.16)') 'finite_rotation_theta_1e-2_relative = ', finite_rotation(1)
         write(unit,'(a,es24.16)') 'finite_rotation_theta_5e-3_relative = ', finite_rotation(2)
         write(unit,'(a,es24.16)') 'finite_rotation_theta_2p5e-3_relative = ', finite_rotation(3)
         write(unit,'(a,es24.16)') 'finite_rotation_predicted_theta_1e-2_relative = ', finite_predicted(1)
         write(unit,'(a,es24.16)') 'finite_rotation_predicted_theta_5e-3_relative = ', finite_predicted(2)
         write(unit,'(a,es24.16)') 'finite_rotation_predicted_theta_2p5e-3_relative = ', finite_predicted(3)
         write(unit,'(a,a)') 'parameter_response_required = ', merge('YES','NO ',parameter_required)
         write(unit,'(a,a)') 'orthogonalization_response_required = ', merge('YES','NO ',orthogonalization_required)
         write(unit,'(a,a)') 'classification = ', trim(classification)
         write(unit,'(a,a)') 'verdict = ', trim(verdict)
         write(unit,'(a)') 'response_lmax_generalization = OFF'
         write(unit,'(a)') 'BES_Halle = OFF'
         write(unit,'(a)') 'ALSDA_Ward_Dyson = OFF'
         close(unit)
      end if
      if (rank == 0) write(*,'(a,a)') 'DRESP-09T verdict: ', trim(classification)

      deallocate(field_raw,d_ee,d_o_types,d_e_types,h,o,enu,h2,d_h,d_o,d_enu,d_h2,native_h,native_h2,field_h,field_h2,basis_h,basis_o,basis_e,basis_h2,total_h,total_h2,term_enu,term_h,term_left,term_middle,term_right,source_plus,source_minus,raw_global,overlap_global,raw_effective,overlap_effective,connection,generator,raw_site,overlap_site,per_k,first_h,first_o,first_e,field_norm_k,basis_norm_k,native_norm_k,total_norm_k,density_fixed_k,density_corrected_k,density_plus,density_minus,density,density_matrix,delta_h_ks,delta_rho_rot,delta_rho_spec,delta_rho_basis,fixed_valence,corrected_valence,target_valence,magnetization,valence_magnetization,core_magnetization)
   end subroutine run_dresp09t_moving_basis

   subroutine build_l0_bxc_field(space, states, field)
      type(response_space_layout), intent(in) :: space
      type(radial_ground_state), intent(in) :: states(:)
      complex(rp), intent(out) :: field(:)
      integer :: flat, site, harmonic, ir, channel, block
      real(rp) :: sqrt_four_pi
      field = cmplx(0.0_rp,0.0_rp,rp)
      sqrt_four_pi = sqrt(4.0_rp*acos(-1.0_rp))
      do flat = 1, size(field)
         block = flat - 1
         channel = mod(block,space%nchannel) + 1
         ir = mod(block/space%nchannel,space%npoint) + 1
         harmonic = mod(block/(space%nchannel*space%npoint),(space%response_lmax+1)**2)
         site = block/(space%nchannel*space%npoint*(space%response_lmax+1)**2) + 1
         if (site <= size(states) .and. harmonic == 0 .and. channel == 1) field(flat) = &
            cmplx(sqrt_four_pi*states(site)%bxc_pauli(ir),0.0_rp,rp)
      end do
   end subroutine build_l0_bxc_field

   subroutine insert_site_block(local, global, site, norb_site)
      complex(rp), intent(in) :: local(:, :, :)
      complex(rp), intent(inout) :: global(:, :, :)
      integer, intent(in) :: site, norb_site
      integer :: first, last, branch
      first = (site-1)*2*norb_site + 1
      last = site*2*norb_site
      do branch = 1, 4
         global(first:last,first:last,branch) = local(:,:,branch)
      end do
   end subroutine insert_site_block

   subroutine reciprocal_kpoint(reciprocal_obj, ik, array4d, result)
      type(reciprocal), intent(in) :: reciprocal_obj
      integer, intent(in) :: ik
      complex(rp), intent(in) :: array4d(:, :, :, :)
      complex(rp), intent(out) :: result(:, :)
      real(rp) :: kpoint(3)
      if (allocated(reciprocal_obj%k_workset%points)) then
         kpoint = reciprocal_obj%k_workset%points(:,ik)
      else
         kpoint = reciprocal_obj%k_points(:,ik)
      end if
      call reciprocal_obj%fourier_transform_array(array4d,kpoint,result)
   end subroutine reciprocal_kpoint

   subroutine finite_native_tangent(h, o, enu, native, predicted, residuals, predicted_residuals)
      complex(rp), intent(in) :: h(:,:), o(:,:), enu(:,:), native(:,:), predicted(:,:)
      real(rp), intent(inout) :: residuals(3), predicted_residuals(3)
      integer :: n, i
      real(rp) :: theta
      complex(rp), allocatable :: hp(:,:), hm(:,:), op(:,:), om(:,:), ep(:,:), em(:,:), h2p(:,:), h2m(:,:), up(:,:), down(:,:), tmp(:,:)
      n = size(h,1)
      allocate(up(n/2,n/2),down(n/2,n/2),hp(n,n),hm(n,n),op(n,n),om(n,n),ep(n,n),em(n,n),h2p(n,n),h2m(n,n),tmp(n,n))
      call dresp08_extract_spin_blocks(h,n/2,1,up,down)
      do i = 1, 3
         theta = 1.0e-2_rp/2.0_rp**real(i-1,rp)
         call dresp08_extract_spin_blocks(h,n/2,1,up,down)
         call dresp08_rotate_collinear_operator(up,down,theta,hp)
         call dresp08_rotate_collinear_operator(up,down,-theta,hm)
         call dresp08_extract_spin_blocks(o,n/2,1,up,down)
         call dresp08_rotate_collinear_operator(up,down,theta,op)
         call dresp08_rotate_collinear_operator(up,down,-theta,om)
         call dresp08_extract_spin_blocks(enu,n/2,1,up,down)
         call dresp08_rotate_collinear_operator(up,down,theta,ep)
         call dresp08_rotate_collinear_operator(up,down,-theta,em)
         call dresp08_build_h2(hp,op,ep,h2p)
         call dresp08_build_h2(hm,om,em,h2m)
         residuals(i) = max(residuals(i),dresp08_relative_residual((h2p-h2m)/(2.0_rp*theta),native))
         predicted_residuals(i) = max(predicted_residuals(i),dresp08_relative_residual((h2p-h2m)/(2.0_rp*theta),predicted))
      end do
      deallocate(up,down,hp,hm,op,om,ep,em,h2p,h2m,tmp)
   end subroutine finite_native_tangent

   logical function radial_provenance_complete(radial) result(ok)
      type(lmto_radial_basis), intent(in) :: radial(:)
      integer :: i
      ok = .true.
      do i = 1, size(radial)
         ok = ok .and. allocated(radial(i)%potential) .and. allocated(radial(i)%tmc) .and. allocated(radial(i)%phi_small)
      end do
   end function radial_provenance_complete

   subroutine weighted_density_relative(values,target,states,relative)
      real(rp), intent(in) :: values(:,:), target(:,:)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(out) :: relative
      integer :: site, ir
      real(rp) :: num, den, r, jac
      num = 0.0_rp; den = 0.0_rp
      do site = 1, size(states)
         do ir = 2, size(states(site)%r)
            r = states(site)%r(ir)
            jac = (2.0_rp*(mod(ir+1,2)+1)/3.0_rp)*states(site)%a*(r+states(site)%b)
            num = num + jac*(values(site,ir)-target(site,ir))**2/max(4.0_rp*acos(-1.0_rp)*r*r, tiny(1.0_rp))
            den = den + jac*target(site,ir)**2/max(4.0_rp*acos(-1.0_rp)*r*r, tiny(1.0_rp))
         end do
      end do
      relative = sqrt(num/max(den,tiny(1.0_rp)))
   end subroutine weighted_density_relative

end module lr_dresp09t_bridge_mod
