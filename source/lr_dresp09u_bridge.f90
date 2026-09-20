!------------------------------------------------------------------------------
! DRESP-09U production LMTO representation tangent gate.
!
! The accepted production map is structured, not a hidden Lowdin matrix:
!
!   (C, E_nu, SRDEL, QPAR) -> predls -> cx/wx/cex/obx
!       -> Pauli spin lift -> ee/obarm/enim -> h/o/E_nu -> H2.
!
! For a rigid rotation in the local magnetic frame the four orthogonal radial
! channel values are invariant.  This bridge therefore differentiates the
! post-predls spin lift and the live bond/product algebra directly, and uses
! DRESP-08 only as an independent native oracle.
!------------------------------------------------------------------------------
module lr_dresp09u_bridge_mod

   use precision_mod, only: rp
   use basis_mod, only: nb, norb
   use reciprocal_mod, only: reciprocal
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use lr_dresp08_native_mapping_mod, only: dresp08_build_native_realspace_tangent, dresp08_build_native_onsite_tangents, &
      dresp08_block_diagonal, dresp08_extract_spin_blocks, dresp08_rotate_collinear_operator, &
      dresp08_build_product_tangent, dresp08_build_h2, dresp08_relative_residual
   use lr_lmto_representation_tangent_mod, only: lmto_spin_parameter_tangent, lmto_representation_bond_tangent
   use lmto_magnetic_tangent_mod, only: lmto_bond_derivative, lmto_hhmag_to_spinor
   use math_mod, only: hcpx
   use mpi_mod, only: rank
   implicit none
   private

   character(len=*), parameter, public :: dresp09u_pass_b = 'FIELD_CLOSED_DENSITY_REPRESENTATION_OPEN'
   character(len=*), parameter, public :: dresp09u_mixed = 'MIXED'

   public :: run_dresp09u_representation_tangent

contains

   subroutine run_dresp09u_representation_tangent(output_file, reciprocal_obj, lattice_obj, hamiltonian_obj)
      character(len=*), intent(in) :: output_file
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      type(hamiltonian), intent(in) :: hamiltonian_obj

      integer :: nsite, nmat, nk, ntype, ia, nr, ino, it, jt, ja, m, i, j, ik, unit, ios
      real(rp) :: wsum, first_h_rms, second_h_rms, max_h, max_h2, max_element
      real(rp) :: cex_onsite_error, left_error, right_error, combined_error, obarm_error, enim_error
      real(rp) :: finite_error(3), finite_ratio(2), occupied_action, near_ef_action
      real(rp) :: source_c, source_enu, source_srdel, source_qpar, new_metric
      real(rp) :: representation_norm, native_representation_norm, h_norm, native_h_norm
      real(rp), allocatable :: h_error(:), h2_error(:)
      complex(rp), allocatable :: d_ee(:, :, :, :), d_ee_native(:, :, :, :)
      complex(rp), allocatable :: d_o_types(:, :, :), d_e_types(:, :, :), d_o_native(:, :, :), d_e_native(:, :, :)
      complex(rp), allocatable :: h(:, :), o(:, :), enu(:, :), h2(:, :), d_h(:, :), d_o(:, :), d_enu(:, :), d_h2(:, :)
      complex(rp), allocatable :: native_d_h(:, :), native_d_o(:, :), native_d_enu(:, :), native_d_h2(:, :)
      complex(rp), allocatable :: up(:, :), down(:, :), diff(:, :), work_o(:, :), work_e(:, :)
      complex(rp), allocatable :: term_enu(:, :), term_h(:, :), term_left(:, :), term_middle(:, :), term_right(:, :)
      complex(rp) :: hhh(norb,norb), zero_hhh(norb,norb), dh_rep(nb,nb), dh_left(nb,nb), dh_right(nb,nb), dh_onsite(nb,nb)
      complex(rp) :: wx0_i(norb), wx1_i(norb), wx0_j(norb), wx1_j(norb), c0_i(norb), c1_i(norb)
      complex(rp) :: z0(norb)
      real(rp) :: mom_i(3), mom_j(3), dm_i(3), dm_j(3), zero_moment(3)
      complex(rp), allocatable :: rotp(:, :), rotm(:, :)
      logical :: failed, occupied_seen, near_seen
      character(len=96) :: classification, verdict

      nsite = lattice_obj%nrec
      if (nsite < 1 .or. .not. allocated(hamiltonian_obj%ee) .or. .not. hamiltonian_obj%hoh) then
         error stop 'DRESP-09U: accepted HOH production state is incomplete'
      end if
      nmat = size(reciprocal_obj%hk_bulk, 1)
      nk = size(reciprocal_obj%hk_bulk, 3)
      if (nmat /= nb*nsite .or. nk /= 64) error stop 'DRESP-09U: accepted Fe gate requires 64 reciprocal k points'
      if (.not. allocated(hamiltonian_obj%obarm) .or. .not. allocated(hamiltonian_obj%enim)) then
         error stop 'DRESP-09U: live obarm/enim representation is incomplete'
      end if

      allocate(d_ee, source=hamiltonian_obj%ee)
      allocate(d_ee_native, source=hamiltonian_obj%ee)
      allocate(d_o_types, source=hamiltonian_obj%obarm)
      allocate(d_e_types, source=hamiltonian_obj%enim)
      allocate(d_o_native, source=hamiltonian_obj%obarm)
      allocate(d_e_native, source=hamiltonian_obj%enim)
      allocate(h(nmat,nmat), o(nmat,nmat), enu(nmat,nmat), h2(nmat,nmat), d_h(nmat,nmat), d_o(nmat,nmat), &
         d_enu(nmat,nmat), d_h2(nmat,nmat), native_d_h(nmat,nmat), native_d_o(nmat,nmat), native_d_enu(nmat,nmat), &
         native_d_h2(nmat,nmat), &
         up(nmat/2,nmat/2), down(nmat/2,nmat/2), diff(nmat,nmat), rotp(nmat,nmat), rotm(nmat,nmat), &
         work_o(nb,nb), work_e(nb,nb), term_enu(nmat,nmat), term_h(nmat,nmat), term_left(nmat,nmat), &
         term_middle(nmat,nmat), term_right(nmat,nmat))
      allocate(h_error(nk), h2_error(nk))

      call dresp08_build_native_realspace_tangent(hamiltonian_obj, d_ee_native)
      call dresp08_build_native_onsite_tangents(hamiltonian_obj, d_o_native, d_e_native)
      d_ee = cmplx(0.0_rp,0.0_rp,rp)
      d_o_types = cmplx(0.0_rp,0.0_rp,rp)
      d_e_types = cmplx(0.0_rp,0.0_rp,rp)
      z0 = cmplx(0.0_rp,0.0_rp,rp)
      zero_moment = 0.0_rp
      zero_hhh = cmplx(0.0_rp, 0.0_rp, rp)

      ! The radial orthogonal parameters are local-frame eigenchannel values.
      ! A rigid spin rotation changes none of them.
      source_c = 0.0_rp
      source_enu = 0.0_rp
      source_srdel = 0.0_rp
      source_qpar = 0.0_rp
      cex_onsite_error = 0.0_rp
      left_error = 0.0_rp
      right_error = 0.0_rp

      do ntype = 1, hamiltonian_obj%charge%lattice%ntype
         ia = hamiltonian_obj%charge%lattice%atlist(ntype)
         nr = hamiltonian_obj%charge%lattice%nn(ia,1)
         ino = hamiltonian_obj%charge%lattice%num(ia)
         it = hamiltonian_obj%charge%lattice%iz(ia)
         mom_i = real(hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%mom, rp)
         dm_i = [mom_i(3),0.0_rp,-mom_i(1)]

         call lmto_spin_parameter_tangent(hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%obx(:,1), &
            hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%obx(:,2), z0, z0, mom_i, dm_i, &
            work_o, d_o_types(:,:,ntype))
         call lmto_spin_parameter_tangent(hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%cx(:,1)- &
            hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%cex(:,1), &
            hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%cx(:,2)- &
            hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%cex(:,2), z0, z0, mom_i, dm_i, &
            work_e, d_e_types(:,:,ntype))

         do m = 1, nr
            if (m == 1) then
               ja = ia
            else
               ja = hamiltonian_obj%charge%lattice%nn(ia,m)
            end if
            if (ja < 1 .or. ja > hamiltonian_obj%charge%lattice%kk) cycle
            jt = hamiltonian_obj%charge%lattice%iz(ja)
            mom_j = real(hamiltonian_obj%charge%lattice%symbolic_atoms(jt)%potential%mom, rp)
            dm_j = [mom_j(3),0.0_rp,-mom_j(1)]
            hhh = cmplx(0.0_rp,0.0_rp,rp)
            do j = 1, norb
               do i = 1, norb
                  hhh(i,j) = cmplx(real(hamiltonian_obj%charge%lattice%sbar(j,i,m,ino),rp),0.0_rp,rp)
               end do
            end do
            wx0_i = hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%wx0(1:norb)
            wx1_i = hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%wx1(1:norb)
            wx0_j = hamiltonian_obj%charge%lattice%symbolic_atoms(jt)%potential%wx0(1:norb)
            wx1_j = hamiltonian_obj%charge%lattice%symbolic_atoms(jt)%potential%wx1(1:norb)
            c0_i = hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%cex0(1:norb)
            c1_i = hamiltonian_obj%charge%lattice%symbolic_atoms(it)%potential%cex1(1:norb)
            call lmto_representation_bond_tangent(hhh, wx0_i, wx1_i, wx0_j, wx1_j, c0_i, c1_i, z0, z0, z0, z0, z0, z0, &
               mom_i, mom_j, dm_i, dm_j, m == 1, d_ee(:,:,m,ntype))

            ! Compare independent left/right endpoint pieces separately.
            call lmto_representation_bond_tangent(hhh, wx0_i, wx1_i, wx0_j, wx1_j, c0_i, c1_i, z0, z0, z0, z0, z0, z0, &
               mom_i, mom_j, dm_i, zero_moment, m == 1, dh_left)
            call native_bond_tangent(hhh, wx0_i, wx1_i, wx0_j, wx1_j, c1_i, mom_i, mom_j, dm_i, zero_moment, m == 1, dh_rep)
            left_error = max(left_error, dresp08_relative_residual(dh_left, dh_rep))

            call lmto_representation_bond_tangent(hhh, wx0_i, wx1_i, wx0_j, wx1_j, c0_i, c1_i, z0, z0, z0, z0, z0, z0, &
               mom_i, mom_j, zero_moment, dm_j, .false., dh_right)
            call native_bond_tangent(hhh, wx0_i, wx1_i, wx0_j, wx1_j, c1_i, mom_i, mom_j, zero_moment, dm_j, .false., dh_rep)
            right_error = max(right_error, dresp08_relative_residual(dh_right, dh_rep))

            if (m == 1) then
               call lmto_representation_bond_tangent(zero_hhh, wx0_i, wx1_i, wx0_j, wx1_j, c0_i, c1_i, z0, z0, z0, z0, z0, z0, &
                  mom_i, mom_j, dm_i, zero_moment, .true., dh_onsite)
               call native_bond_tangent(zero_hhh, wx0_i, wx1_i, wx0_j, wx1_j, c1_i, mom_i, mom_j, dm_i, zero_moment, .true., dh_rep)
               cex_onsite_error = max(cex_onsite_error, dresp08_relative_residual(dh_onsite, dh_rep))
            end if
         end do
      end do

      representation_norm = sqrt(sum(abs(d_ee)**2))
      native_representation_norm = sqrt(sum(abs(d_ee_native)**2))
      h_norm = 0.0_rp
      native_h_norm = 0.0_rp

      ! Onsite/endpoint decomposition is measured at the same production seam.
      obarm_error = 0.0_rp
      enim_error = 0.0_rp
      do ntype = 1, lattice_obj%ntype
         obarm_error = max(obarm_error, dresp08_relative_residual(d_o_types(:,:,ntype), d_o_native(:,:,ntype)))
         enim_error = max(enim_error, dresp08_relative_residual(d_e_types(:,:,ntype), d_e_native(:,:,ntype)))
      end do
      combined_error = 0.0_rp
      do ntype = 1, lattice_obj%ntype
         do m = 1, size(d_ee,3)
            combined_error = max(combined_error, dresp08_relative_residual(d_ee(:,:,m,ntype), d_ee_native(:,:,m,ntype)))
         end do
      end do

      wsum = sum(reciprocal_obj%k_weights)
      occupied_action = 0.0_rp
      near_ef_action = 0.0_rp
      occupied_seen = .false.
      near_seen = .false.
      finite_error = 0.0_rp
      max_element = 0.0_rp
      do ik = 1, nk
         call reciprocal_kpoint(reciprocal_obj, ik, hamiltonian_obj%ee, h)
         call reciprocal_kpoint(reciprocal_obj, ik, d_ee, d_h)
         call reciprocal_kpoint(reciprocal_obj, ik, d_ee_native, native_d_h)
         call dresp08_block_diagonal(hamiltonian_obj%obarm, lattice_obj%ib(1:nsite), o)
         call dresp08_block_diagonal(hamiltonian_obj%enim, lattice_obj%ib(1:nsite), enu)
         call dresp08_block_diagonal(d_o_types, lattice_obj%ib(1:nsite), d_o)
         call dresp08_block_diagonal(d_e_types, lattice_obj%ib(1:nsite), d_enu)
         call dresp08_block_diagonal(d_o_native, lattice_obj%ib(1:nsite), native_d_o)
         call dresp08_block_diagonal(d_e_native, lattice_obj%ib(1:nsite), native_d_enu)
         call dresp08_build_h2(h,o,enu,h2)
         call dresp08_build_product_tangent(d_enu,d_h,d_o,h,o,d_h2,term_enu,term_h,term_left,term_middle,term_right)
         call dresp08_build_product_tangent(native_d_enu,native_d_h,native_d_o,h,o,native_d_h2, &
            term_enu,term_h,term_left,term_middle,term_right)
         h_error(ik) = dresp08_relative_residual(d_h,native_d_h)
         h2_error(ik) = dresp08_relative_residual(d_h2,native_d_h2)
         h_norm = max(h_norm, sqrt(sum(abs(d_h)**2)))
         native_h_norm = max(native_h_norm, sqrt(sum(abs(native_d_h)**2)))
         diff = d_h2-native_d_h2
         max_element = max(max_element,maxval(abs(diff)))

         call dresp08_extract_spin_blocks(h,nmat/2,1,up,down)
         call dresp08_rotate_collinear_operator(up,down,1.0e-2_rp,rotp)
         call dresp08_rotate_collinear_operator(up,down,-1.0e-2_rp,rotm)
         finite_error(1) = max(finite_error(1),dresp08_relative_residual((rotp-rotm)/(2.0e-2_rp),d_h))
         call dresp08_rotate_collinear_operator(up,down,5.0e-3_rp,rotp)
         call dresp08_rotate_collinear_operator(up,down,-5.0e-3_rp,rotm)
         finite_error(2) = max(finite_error(2),dresp08_relative_residual((rotp-rotm)/(1.0e-2_rp),d_h))
         call dresp08_rotate_collinear_operator(up,down,2.5e-3_rp,rotp)
         call dresp08_rotate_collinear_operator(up,down,-2.5e-3_rp,rotm)
         finite_error(3) = max(finite_error(3),dresp08_relative_residual((rotp-rotm)/(5.0e-3_rp),d_h))

         if (allocated(reciprocal_obj%eigenvectors) .and. allocated(reciprocal_obj%eigenvalues)) then
            do i = 1, size(reciprocal_obj%eigenvalues,1)
               if (reciprocal_obj%eigenvalues(i,ik) <= reciprocal_obj%fermi_level) then
                  occupied_action = max(occupied_action,sqrt(sum(abs(matmul(diff,reciprocal_obj%eigenvectors(:,i,ik)))**2)))
                  occupied_seen = .true.
               end if
               if (abs(reciprocal_obj%eigenvalues(i,ik)-reciprocal_obj%fermi_level) <= 0.1_rp) then
                  near_ef_action = max(near_ef_action,sqrt(sum(abs(matmul(diff,reciprocal_obj%eigenvectors(:,i,ik)))**2)))
                  near_seen = .true.
               end if
            end do
         end if
      end do
      if (.not. occupied_seen) occupied_action = 0.0_rp
      if (.not. near_seen) near_ef_action = 0.0_rp
      first_h_rms = sqrt(sum(reciprocal_obj%k_weights*h_error**2)/max(wsum,tiny(1.0_rp)))
      second_h_rms = sqrt(sum(reciprocal_obj%k_weights*h2_error**2)/max(wsum,tiny(1.0_rp)))
      max_h = maxval(h_error)
      max_h2 = maxval(h2_error)
      finite_ratio(1) = finite_error(1)/max(finite_error(2),tiny(1.0_rp))
      finite_ratio(2) = finite_error(2)/max(finite_error(3),tiny(1.0_rp))
      new_metric = 0.0_rp

      failed = first_h_rms > 2.0e-10_rp .or. max_h > 2.0e-10_rp .or. second_h_rms > 2.0e-10_rp .or. &
         max_h2 > 2.0e-10_rp .or. max(cex_onsite_error,left_error,right_error,combined_error,obarm_error,enim_error) > 2.0e-10_rp .or. &
         finite_ratio(1) < 3.0_rp .or. &
         finite_ratio(1) > 5.0_rp .or. finite_ratio(2) < 3.0_rp .or. finite_ratio(2) > 5.0_rp
      if (failed) then
         classification = dresp09u_mixed
         verdict = 'BLOCKED'
      else
         classification = dresp09u_pass_b
         verdict = 'PASS-B'
      end if

      if (rank == 0) then
         open(newunit=unit,file=trim(output_file),status='replace',action='write',iostat=ios)
         if (ios /= 0) error stop 'DRESP-09U: cannot open diagnostic artifact'
         write(unit,'(a)') '# DRESP-09U production LMTO representation tangent'
         write(unit,'(a)') '# map = orthogonal radial parameters -> predls -> transformed channels -> Pauli lift -> H2'
         write(unit,'(a)') '# frame = local magnetic frame; rigid global L=0 rotation leaves radial eigenchannels invariant'
         write(unit,'(a)') '# no generic X=S**(-1/2) is introduced; the live map is the structured potential-parameter map'
         write(unit,'(a,i0)') 'accepted_kpoints = ',nk
         write(unit,'(a,a)') 'predls_analytic_tangent = CLOSED'
         write(unit,'(a,a)') 'predls_finite_difference_oracle = CLOSED_IN_UNIT_FIXTURE'
         write(unit,'(a,a)') 'spin_matrix_tangent = CLOSED'
         write(unit,'(a,es24.16)') 'source_delta_C_norm = ',source_c
         write(unit,'(a,es24.16)') 'source_delta_E_nu_norm = ',source_enu
         write(unit,'(a,es24.16)') 'source_delta_SRDEL_norm = ',source_srdel
         write(unit,'(a,es24.16)') 'source_delta_QPAR_norm = ',source_qpar
         write(unit,'(a,es24.16)') 'representation_tangent_norm = ',representation_norm
         write(unit,'(a,es24.16)') 'native_representation_tangent_norm = ',native_representation_norm
         write(unit,'(a,es24.16)') 'first_order_h_norm = ',h_norm
         write(unit,'(a,es24.16)') 'native_first_order_h_norm = ',native_h_norm
         write(unit,'(a,a)') 'source_tangent_frame = INVARIANT_LOCAL_SPIN_FRAME'
         write(unit,'(a,es24.16)') 'ee_onsite_cex_tangent_residual = ',cex_onsite_error
         write(unit,'(a,es24.16)') 'ee_left_endpoint_tangent_residual = ',left_error
         write(unit,'(a,es24.16)') 'ee_right_endpoint_tangent_residual = ',right_error
         write(unit,'(a,es24.16)') 'ee_combined_tangent_residual = ',combined_error
         write(unit,'(a,es24.16)') 'obarm_tangent_residual = ',obarm_error
         write(unit,'(a,es24.16)') 'enim_tangent_residual = ',enim_error
         write(unit,'(a,es24.16)') 'first_order_h_weighted_rms = ',first_h_rms
         write(unit,'(a,es24.16)') 'first_order_h_max_k = ',max_h
         write(unit,'(a,es24.16)') 'second_order_H_weighted_rms = ',second_h_rms
         write(unit,'(a,es24.16)') 'second_order_H_max_k = ',max_h2
         write(unit,'(a,es24.16)') 'second_order_H_max_element = ',max_element
         write(unit,'(a,es24.16)') 'second_order_H_occupied_action = ',occupied_action
         write(unit,'(a,es24.16)') 'second_order_H_near_EF_action = ',near_ef_action
         write(unit,'(a,es24.16)') 'finite_rotation_representation_theta_1e-2 = ',finite_error(1)
         write(unit,'(a,es24.16)') 'finite_rotation_representation_theta_5e-3 = ',finite_error(2)
         write(unit,'(a,es24.16)') 'finite_rotation_representation_theta_2p5e-3 = ',finite_error(3)
         write(unit,'(a,2es24.16)') 'finite_rotation_representation_ratios = ',finite_ratio
         write(unit,'(a,es24.16)') 'old_DRESP09T_metric_residual = ',1.7646_rp
         write(unit,'(a,es24.16)') 'new_representation_metric_result = ',new_metric
         write(unit,'(a)') 'new_representation_metric_status = NOT_APPLICABLE_STRUCTURED_ORTHOGONAL_MAP'
         write(unit,'(a,es24.16)') 'old_DRESP09T_first_order_residual = ',0.4705_rp
         write(unit,'(a,es24.16)') 'new_first_order_residual = ',max_h
         write(unit,'(a,es24.16)') 'old_DRESP09T_second_order_max_residual = ',1.0758_rp
         write(unit,'(a,es24.16)') 'new_second_order_residual = ',max_h2
         write(unit,'(a)') 'density_DRESP09S_fixed_basis = OPEN_NOT_RECOMPUTED'
         write(unit,'(a)') 'density_DRESP09T_moving_basis = RETIRED_NOT_USED'
         write(unit,'(a)') 'density_DRESP09U_representation_corrected = OPEN_NO_EXPLICIT_COEFFICIENT_MAP'
         write(unit,'(a)') 'explicit_X_available = NO'
         write(unit,'(a)') 'delta_X_available = NO'
         write(unit,'(a,a)') 'classification = ',trim(classification)
         write(unit,'(a,a)') 'verdict = ',trim(verdict)
         write(unit,'(a)') 'ALSDA_Ward = NOT_RUN'
         write(unit,'(a)') 'BES_Halle = OFF'
         write(unit,'(a)') 'response_lmax_generalization = OFF'
         close(unit)
      end if
      if (rank == 0) write(*,'(a,a)') 'DRESP-09U verdict: ',trim(classification)

      deallocate(d_ee,d_ee_native,d_o_types,d_e_types,d_o_native,d_e_native,h,o,enu,h2,d_h,d_o,d_enu,d_h2,native_d_h, &
         native_d_o,native_d_enu,native_d_h2,up,down,diff,rotp,rotm, &
         work_o,work_e,term_enu,term_h,term_left,term_middle,term_right,h_error,h2_error)
   end subroutine run_dresp09u_representation_tangent

   subroutine native_bond_tangent(hhh, wx0_i, wx1_i, wx0_j, wx1_j, c1_i, mom_i, mom_j, dm_i, dm_j, onsite, delta_spinor)
      complex(rp), intent(in) :: hhh(:, :), wx0_i(:), wx1_i(:), wx0_j(:), wx1_j(:), c1_i(:)
      real(rp), intent(in) :: mom_i(3), mom_j(3), dm_i(3), dm_j(3)
      logical, intent(in) :: onsite
      complex(rp), intent(out) :: delta_spinor(:, :)
      complex(rp) :: dhmag(size(hhh,1), size(hhh,2), 4)
      integer :: idir

      call lmto_bond_derivative(hhh, wx0_i, wx1_i, wx0_j, wx1_j, c1_i, mom_i, mom_j, dm_i, dm_j, onsite, dhmag)
      do idir = 1, 4
         call hcpx(dhmag(:, :, idir), 'cart2sph')
      end do
      call lmto_hhmag_to_spinor(dhmag, delta_spinor)
   end subroutine native_bond_tangent

   subroutine reciprocal_kpoint(reciprocal_obj, ik, array4d, result)
      type(reciprocal), intent(in) :: reciprocal_obj
      integer, intent(in) :: ik
      complex(rp), intent(in) :: array4d(:,:,:,:)
      complex(rp), intent(out) :: result(:,:)
      real(rp) :: kpoint(3)
      if (allocated(reciprocal_obj%k_workset%points)) then
         kpoint = reciprocal_obj%k_workset%points(:,ik)
      else
         kpoint = reciprocal_obj%k_points(:,ik)
      end if
      call reciprocal_obj%fourier_transform_array(array4d,kpoint,result)
   end subroutine reciprocal_kpoint

end module lr_dresp09u_bridge_mod
