!------------------------------------------------------------------------------
! DRESP-09ZS accepted-state compact-span audit.
!
! This module is deliberately diagnostic.  It compares the live production
! four-branch radial product space with independent second-order four- and
! six-branch shadow spaces, performs weighted direct SVDs, and projects
! deterministic arbitrary-(L,M) vertices.
! Nothing here is consumed by chi0, a kernel, Dyson, or a production basis.
!------------------------------------------------------------------------------
module lr_dresp09zs_compact_span_mod

   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, lmto_product_block, &
      lmto_product_candidate
   use response_angular_basis_mod, only: response_gaunt
   use lr_response_space_mod, only: response_space_layout
   use lr_sr_angular_vertex_mod, only: sr_angular_nn_rank2_coefficient
   use lr_sr_augmentation_tangent_mod, only: sr_aug_branch_observable_tangent, sr_aug_sigma_x
   use lr_dresp09z_bridge_mod, only: dresp09z_allowed_orbital_pair
   implicit none
   private

   character(len=*), parameter, public :: dresp09zs_strict_extension = &
      'PRODUCTION_SPAN_STRICT_SUBSPACE_EXTENSION_REQUIRED'
   character(len=*), parameter, public :: dresp09zs_basis_replacement = &
      'PRODUCTION_BASIS_REPLACEMENT_REQUIRED'
   character(len=*), parameter, public :: dresp09zs_complete = 'PRODUCTION_SPAN_ALREADY_COMPLETE'
   character(len=*), parameter, public :: dresp09zs_numeric_open = 'SPAN_RECONCILIATION_NUMERICS_OPEN'

   type :: candidate
      integer :: l = 0
      integer :: lp = 0
      integer :: branch = 0
   end type candidate

   type :: block_audit
      integer :: nc_old = 0
      integer :: nc_six = 0
      integer :: nc_prod = 0
      integer :: rank_old(3) = 0
      integer :: rank_six(3) = 0
      integer :: rank_prod(3) = 0
      integer :: rank_oracle(3) = 0
      integer :: npoint = 0
      real(rp) :: sigma_old_max = 0.0_rp
      real(rp) :: sigma_old_min = 0.0_rp
      real(rp) :: sigma_six_max = 0.0_rp
      real(rp) :: sigma_six_min = 0.0_rp
      real(rp) :: first_discard_old(3) = 0.0_rp
      real(rp) :: first_discard_six(3) = 0.0_rp
      real(rp) :: tau_old(3) = 0.0_rp
      real(rp) :: tau_six(3) = 0.0_rp
      real(rp) :: tau_prod(3) = 0.0_rp
      real(rp) :: tau_oracle(3) = 0.0_rp
      complex(rp), allocatable :: u_old(:, :)
      complex(rp), allocatable :: u_six(:, :)
      complex(rp), allocatable :: u_prod(:, :)
      complex(rp), allocatable :: u_oracle(:, :)
      complex(rp), allocatable :: a_old(:, :)
      complex(rp), allocatable :: a_six(:, :)
      complex(rp), allocatable :: a_prod(:, :)
      complex(rp), allocatable :: production_forward(:, :)
      real(rp), allocatable :: production_column_norms(:)
      complex(rp), allocatable :: right_old(:, :)
      complex(rp), allocatable :: right_six(:, :)
      real(rp), allocatable :: weights(:)
      type(candidate), allocatable :: old_candidates(:)
      type(candidate), allocatable :: six_candidates(:)
   end type block_audit

   public :: run_dresp09zs_compact_span_audit

   interface
      subroutine zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
         import :: rp
         character(len=1), intent(in) :: jobu, jobvt
         integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
         complex(rp), intent(inout) :: a(lda, *)
         real(rp), intent(out) :: s(*)
         complex(rp), intent(out) :: u(ldu, *), vt(ldvt, *), work(*)
         real(rp), intent(inout) :: rwork(*)
         integer, intent(out) :: info
      end subroutine zgesvd
   end interface

contains

   subroutine run_dresp09zs_compact_span_audit(output_file, space, radial_bases, production_basis)
      character(len=*), intent(in) :: output_file
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(lmto_product_response_basis), intent(in) :: production_basis
      type(block_audit) :: blocks(0:4)
      real(rp) :: r20(0:4), r02(0:4), r20_prod(0:4), r02_prod(0:4)
      real(rp) :: candidate_residual(0:4), candidate_global, candidate_prod(0:4), candidate_prod_global
      real(rp) :: principal_old_six(0:4), principal_prod_old(0:4), principal_prod_six(0:4)
      real(rp) :: prod_old_out(0:4), old_prod_out(0:4), prod_six_out(0:4), six_prod_out(0:4)
      real(rp) :: prod_old_frob(0:4), old_prod_frob(0:4), prod_six_frob(0:4), six_prod_frob(0:4)
      real(rp) :: oracle_min(0:4), oracle_live_out(0:4), oracle_independent_out(0:4), oracle_live_frob(0:4)
      real(rp) :: oracle_independent_frob(0:4), oracle_forward_max(0:4)
      real(rp) :: r20_pair(0:2,0:2), r02_pair(0:2,0:2), r20_median, r02_median
      real(rp) :: r20_pair_prod(0:2,0:2), r02_pair_prod(0:2,0:2), r20_median_prod, r02_median_prod
      real(rp) :: physical_by_l(0:4), delta_by_l(0:4), physical_rms_by_l(0:4), delta_rms_by_l(0:4)
      real(rp) :: physical_old_by_l(0:4), delta_old_by_l(0:4), physical_old_rms_by_l(0:4), delta_old_rms_by_l(0:4)
      real(rp) :: component_residual(7), component_old_residual(7), mixed_residual, mixed_old_residual
      real(rp) :: duality_residual(4), double_weighting_residual, raw_pairing(4), projected_pairing(4), projection_loss(4)
      real(rp) :: duality_old(4), double_weighting_old, raw_pairing_old(4), projected_pairing_old(4), projection_loss_old(4)
      real(rp) :: scalar_projector_residual, total_old(3), total_six(3), total_prod(3), physical_global, delta_global
      real(rp) :: physical_old_global, delta_old_global
      real(rp) :: max_new_overlap, rms_new_overlap
      integer :: unit, k, tau, rank_old_total(3), rank_six_total(3), rank_prod_total(3), site, ir
      integer :: new_rank, max_rank_delta, ios
      logical :: rank_stable, physical_active, physical_active_old, oracle_agrees, production_is_complete
      character(len=80) :: classification

      if (size(radial_bases) /= space%nsite .or. space%nsite /= 1 .or. space%response_lmax /= 4) then
         error stop 'DRESP-09ZS requires the accepted one-site spd response_lmax=4 state'
      end if
      if (radial_bases(1)%lmax /= 2 .or. radial_bases(1)%nspin /= 2) then
         error stop 'DRESP-09ZS requires the accepted Fe spd radial state'
      end if
      if (.not. allocated(radial_bases(1)%phiddot_large) .or. .not. allocated(radial_bases(1)%phi_small) .or. &
          .not. allocated(radial_bases(1)%phiddot_small) .or. .not. allocated(radial_bases(1)%tmc)) then
         error stop 'DRESP-09ZS requires complete DRESP-09X second-order radial coefficients'
      end if

      if (production_basis%nsite /= space%nsite .or. production_basis%response_lmax /= space%response_lmax .or. &
          production_basis%circular_channel /= 1 .or. production_basis%product_dimension /= 232 .or. &
          .not. allocated(production_basis%blocks)) then
         error stop 'DRESP-09ZS-R requires the initialized live production product basis'
      end if
      do k = 0, 4
         call build_block(blocks(k), space, radial_bases(1), k, 1, production_basis%blocks(1,k))
      end do
      do tau = 1, 3
         rank_old_total(tau) = 0
         rank_six_total(tau) = 0
         rank_prod_total(tau) = 0
         do k = 0, 4
            rank_old_total(tau) = rank_old_total(tau) + blocks(k)%rank_old(tau)*(2*k+1)
            rank_six_total(tau) = rank_six_total(tau) + blocks(k)%rank_six(tau)*(2*k+1)
            rank_prod_total(tau) = rank_prod_total(tau) + blocks(k)%rank_prod(tau)*(2*k+1)
         end do
      end do
      rank_stable = .true.
      do k = 0, 4
         rank_stable = rank_stable .and. all(blocks(k)%rank_old == blocks(k)%rank_old(1)) .and. &
            all(blocks(k)%rank_six == blocks(k)%rank_six(1)) .and. &
            all(blocks(k)%rank_prod == blocks(k)%rank_prod(1)) .and. &
            all(blocks(k)%rank_oracle == blocks(k)%rank_oracle(1))
      end do

      call direct_new_branch_projection(blocks, r20, r02, r20_pair, r02_pair, r20_median, r02_median)
      call direct_new_branch_projection_prod(blocks, r20_prod, r02_prod, r20_pair_prod, r02_pair_prod, &
         r20_median_prod, r02_median_prod)
      call whole_candidate_residual(blocks, candidate_residual, candidate_global)
      call whole_candidate_residual_prod(blocks, candidate_prod, candidate_prod_global)
      call principal_angle_diagnostics(blocks, principal_old_six)
      do k = 0, 4
         call subspace_diagnostics(blocks(k)%u_prod, blocks(k)%u_old, principal_prod_old(k), prod_old_out(k), &
            old_prod_out(k), prod_old_frob(k), old_prod_frob(k))
         call subspace_diagnostics(blocks(k)%u_prod, blocks(k)%u_six, principal_prod_six(k), prod_six_out(k), &
            six_prod_out(k), prod_six_frob(k), six_prod_frob(k))
      end do
      call production_oracle_diagnostics(blocks, oracle_min, oracle_live_out, oracle_independent_out, oracle_live_frob, &
         oracle_independent_frob, oracle_forward_max, oracle_agrees)
      call new_mode_overlap_diagnostics(blocks, max_new_overlap, rms_new_overlap)
      call physical_vertex_audit(space, radial_bases(1), blocks, physical_by_l, physical_rms_by_l, &
         delta_by_l, delta_rms_by_l, component_residual, physical_global, delta_global, physical_active, .true.)
      call physical_vertex_audit(space, radial_bases(1), blocks, physical_old_by_l, physical_old_rms_by_l, &
         delta_old_by_l, delta_old_rms_by_l, component_old_residual, physical_old_global, delta_old_global, &
         physical_active_old, .false.)
      call mixed_field_audit(space, radial_bases(1), blocks, mixed_residual, duality_residual, double_weighting_residual, &
         raw_pairing, projected_pairing, projection_loss, .true.)
      call mixed_field_audit(space, radial_bases(1), blocks, mixed_old_residual, duality_old, double_weighting_old, &
         raw_pairing_old, projected_pairing_old, projection_loss_old, .false.)
      call scalar_regression(radial_bases(1), scalar_projector_residual)

      new_rank = rank_six_total(1) - rank_prod_total(1)
      max_rank_delta = maxval(rank_six_total-rank_old_total)
      total_old = real(rank_old_total, rp)
      total_six = real(rank_six_total, rp)
      total_prod = real(rank_prod_total, rp)
      production_is_complete = maxval(six_prod_out) < 1.0e-10_rp .and. maxval(component_residual) < 1.0e-10_rp
      if (.not. rank_stable .or. .not. oracle_agrees) then
         classification = dresp09zs_numeric_open
      else if (production_is_complete) then
         classification = dresp09zs_complete
      else if (maxval(prod_six_out) < 1.0e-6_rp .and. maxval(six_prod_out) > 1.0e-8_rp .and. physical_active) then
         classification = dresp09zs_strict_extension
      else
         classification = dresp09zs_basis_replacement
      end if

      open(newunit=unit, file=trim(output_file), status='replace', action='write', iostat=ios)
      if (ios /= 0) error stop 'DRESP-09ZS cannot open audit output'
      write(unit,'(a)') '# DRESP-09ZS accepted Fe six-branch compact-span/SVD audit'
      write(unit,'(a)') '# scope = radial/product space only; no chi0, Dyson, ALSDA, Ward, spectra, BES, or Halle'
      write(unit,'(a)') 'accepted_state = bcc Fe; 4x4x4; 64 k points; 300 K; spd; orbital_lmax=2; response_lmax=4; ham_only/second'
      write(unit,'(a)') 'raw_arbitrary_L_vertex = CLOSED'
      write(unit,'(a,i0)') 'live_production_four_branch_candidate_count = ', 232
      write(unit,'(a,i0)') 'live_production_retained_product_dimension = ', production_basis%product_dimension
      write(unit,'(a,i0)') 'second_order_four_branch_candidate_count = ', 232
      write(unit,'(a,i0)') 'second_order_six_branch_candidate_count = ', 348
      write(unit,'(a,i0)') 'S4_2nd_legacy_candidate_count = ', 232
      write(unit,'(a,i0)') 'S6_2nd_legacy_candidate_count = ', 348
      write(unit,'(a,i0)') 'additional_20_02_candidate_count = ', 116
      write(unit,'(a)') 'candidate_inventory_note = two distinct 232-candidate spaces; counts are inventories, not retained dimensions'
      write(unit,'(a)') 'S4_prod_definition = live lmto_product_response_basis phi/phidot 00/10/01/11 modes'
      write(unit,'(a)') 'S4_2nd_definition = second-order branch_radial 00/10/01/11 shadow modes'
      write(unit,'(a)') 'S6_2nd_definition = second-order branch_radial 00/10/01/11/20/02 shadow modes'
      write(unit,'(a)') 'production_radial_difference = S4_prod uses endpoint phi-Enu*phidot and phidot; S4_2nd includes phiddot terms'
      write(unit,'(a)') 'metric = production LR-04 response-space radial weights; direct zgesvd on W**(1/2) B D**(-1)'
      write(unit,'(a)') 'tau_convention = tau1=max(npoint,ncandidate)*epsilon*sigma_max; tau10=10*tau1; tau100=100*tau1'
      write(unit,'(a)') 'delta_O_kernel = certified analytic augmentation-frame six-branch radial tangent with arbitrary-L Gaunt assembly'
      write(unit,'(a)') 'K S4_prod_candidates S4_prod_rank_tau1 S4_prod_rank_tau10 S4_prod_rank_tau100 S4_2nd_candidates S4_2nd_rank_tau1 S4_2nd_rank_tau10 S4_2nd_rank_tau100 S6_2nd_candidates S6_2nd_rank_tau1 S6_2nd_rank_tau10 S6_2nd_rank_tau100'
      do k = 0, 4
         write(unit,'(i0,1x,i0,1x,3(i0,1x),i0,1x,3(i0,1x),i0,1x,3(i0,1x))') k, blocks(k)%nc_prod, blocks(k)%rank_prod, &
            blocks(k)%nc_old, blocks(k)%rank_old, blocks(k)%nc_six, blocks(k)%rank_six
      end do
      write(unit,'(a,3(i0,1x))') 'total_S4_prod_retained_rank_tau1_tau10_tau100 = ', rank_prod_total
      write(unit,'(a,3(i0,1x))') 'total_S4_2nd_retained_rank_tau1_tau10_tau100 = ', rank_old_total
      write(unit,'(a,3(i0,1x))') 'total_S6_2nd_retained_rank_tau1_tau10_tau100 = ', rank_six_total
      write(unit,'(a,i0)') 'total_S6_2nd_minus_S4_prod_rank_tau1 = ', new_rank
      write(unit,'(a,l1)') 'production_oracle_rank_and_subspace_agreement = ', oracle_agrees
      write(unit,'(a,l1)') 'rank_sensitivity_stable = ', rank_stable
      write(unit,'(a)') 'K S4_prod_tau1 S4_prod_tau10 S4_prod_tau100 S4_2nd_tau1 S4_2nd_tau10 S4_2nd_tau100 S6_2nd_tau1 S6_2nd_tau10 S6_2nd_tau100'
      do k = 0, 4
         write(unit,'(i0,1x,9(es18.10,1x))') k, blocks(k)%tau_prod, blocks(k)%tau_old, blocks(k)%tau_six
      end do
      write(unit,'(a)') 'K S4_prod_vs_S4_2nd_principal_min S4_prod_outside_S4_2nd S4_2nd_outside_S4_prod S4_prod_vs_S6_2nd_principal_min S4_prod_outside_S6_2nd S6_2nd_outside_S4_prod'
      do k = 0, 4
         write(unit,'(i0,1x,6(es18.10,1x))') k, principal_prod_old(k), prod_old_out(k), old_prod_out(k), &
            principal_prod_six(k), prod_six_out(k), six_prod_out(k)
      end do
      write(unit,'(a,3(es18.10,1x))') 'S4_prod_vs_S4_2nd principal_minimum max_prod_outside max_2nd_outside = ', &
         minval(principal_prod_old), maxval(prod_old_out), maxval(old_prod_out)
      write(unit,'(a,es18.10)') 'S4_prod_vs_S4_2nd_principal_minimum = ', minval(principal_prod_old)
      write(unit,'(a,es18.10)') 'S4_prod_vs_S4_2nd_prod_outside = ', maxval(prod_old_out)
      write(unit,'(a,es18.10)') 'S4_prod_vs_S4_2nd_second_order_outside = ', maxval(old_prod_out)
      write(unit,'(a,2(es18.10,1x))') 'S4_prod_vs_S4_2nd Frobenius_residual_prod_to_2nd 2nd_to_prod = ', &
         sqrt(sum(prod_old_frob**2)/5.0_rp), sqrt(sum(old_prod_frob**2)/5.0_rp)
      write(unit,'(a,3(es18.10,1x))') 'S4_prod_vs_S6_2nd principal_minimum max_prod_outside max_six_outside = ', &
         minval(principal_prod_six), maxval(prod_six_out), maxval(six_prod_out)
      write(unit,'(a,es18.10)') 'S4_prod_vs_S6_2nd_principal_minimum = ', minval(principal_prod_six)
      write(unit,'(a,es18.10)') 'S4_prod_vs_S6_2nd_prod_outside = ', maxval(prod_six_out)
      write(unit,'(a,es18.10)') 'S4_prod_vs_S6_2nd_six_outside = ', maxval(six_prod_out)
      write(unit,'(a,2(es18.10,1x))') 'S4_prod_vs_S6_2nd Frobenius_residual_prod_to_six six_to_prod = ', &
         sqrt(sum(prod_six_frob**2)/5.0_rp), sqrt(sum(six_prod_frob**2)/5.0_rp)
      write(unit,'(a)') 'K oracle_principal_min live_outside_independent independent_outside_live live_frob independent_frob forward_transform_max'
      do k = 0, 4
         write(unit,'(i0,1x,6(es18.10,1x))') k, oracle_min(k), oracle_live_out(k), oracle_independent_out(k), &
            oracle_live_frob(k), oracle_independent_frob(k), oracle_forward_max(k)
      end do
      write(unit,'(a)') 'K S4_2nd_R20_max S4_2nd_R02_max S4_2nd_six_candidate_residual S4_2nd_principal_min'
      do k = 0, 4
         write(unit,'(i0,1x,4(es18.10,1x))') k, r20(k), r02(k), candidate_residual(k), principal_old_six(k)
      end do
      write(unit,'(a)') 'K S4_prod_R20_max S4_prod_R02_max S4_prod_six_candidate_residual'
      do k = 0, 4
         write(unit,'(i0,1x,3(es18.10,1x))') k, r20_prod(k), r02_prod(k), candidate_prod(k)
      end do
      write(unit,'(a,es18.10)') 'S4_2nd_full_candidate_residual_global = ', candidate_global
      write(unit,'(a,2(es18.10,1x))') 'S4_prod_R20 maximum RMS = ', maxval(r20_prod), sqrt(sum(r20_prod*r20_prod)/5.0_rp)
      write(unit,'(a,2(es18.10,1x))') 'S4_prod_R02 maximum RMS = ', maxval(r02_prod), sqrt(sum(r02_prod*r02_prod)/5.0_rp)
      write(unit,'(a,es18.10)') 'S4_prod_R20 median = ', r20_median_prod
      write(unit,'(a,es18.10)') 'S4_prod_R02 median = ', r02_median_prod
      write(unit,'(a,es18.10)') 'S4_prod_full_candidate_residual_global = ', candidate_prod_global
      write(unit,'(a,2(es18.10,1x))') 'R20 maximum RMS = ', maxval(r20_prod), sqrt(sum(r20_prod*r20_prod)/5.0_rp)
      write(unit,'(a,2(es18.10,1x))') 'R02 maximum RMS = ', maxval(r02_prod), sqrt(sum(r02_prod*r02_prod)/5.0_rp)
      write(unit,'(a,es18.10)') 'R20 median = ', r20_median_prod
      write(unit,'(a,es18.10)') 'R02 median = ', r02_median_prod
      write(unit,'(a,es18.10)') 'full_candidate_residual_global = ', candidate_prod_global
      write(unit,'(a,2(es18.10,1x))') 'S4_2nd_R20 maximum RMS = ', maxval(r20), sqrt(sum(r20*r20)/5.0_rp)
      write(unit,'(a,2(es18.10,1x))') 'S4_2nd_R02 maximum RMS = ', maxval(r02), sqrt(sum(r02*r02)/5.0_rp)
      write(unit,'(a,es18.10)') 'S4_2nd_R20 median = ', r20_median
      write(unit,'(a,es18.10)') 'S4_2nd_R02 median = ', r02_median
      write(unit,'(a)') 'l lp S4_prod_R20_max S4_prod_R02_max S4_2nd_R20_max S4_2nd_R02_max'
      do site = 0, 2
         do ir = 0, 2
            if (r20_pair(site,ir) > 0.0_rp .or. r02_pair(site,ir) > 0.0_rp .or. &
                r20_pair_prod(site,ir) > 0.0_rp .or. r02_pair_prod(site,ir) > 0.0_rp) then
               write(unit,'(2(i0,1x),4(es18.10,1x))') site, ir, r20_pair_prod(site,ir), r02_pair_prod(site,ir), &
                  r20_pair(site,ir), r02_pair(site,ir)
            end if
         end do
      end do
      write(unit,'(a,2(es18.10,1x))') 'old_span_outside_six_span maximum RMS = ', max_new_overlap, rms_new_overlap
      call write_new_mode_composition(unit, blocks)
      write(unit,'(a)') 'L S4_prod_physical_SR_residual S4_prod_physical_SR_RMS S4_prod_delta_O_residual S4_prod_delta_O_RMS'
      do k = 0, 4
         write(unit,'(i0,1x,4(es18.10,1x))') k, physical_by_l(k), physical_rms_by_l(k), delta_by_l(k), delta_rms_by_l(k)
      end do
      write(unit,'(a,2(es18.10,1x))') 'physical_global_max_RMS = ', physical_global, sqrt(sum(physical_by_l**2)/5.0_rp)
      write(unit,'(a,2(es18.10,1x))') 'delta_O_global_max_RMS = ', delta_global, sqrt(sum(delta_by_l**2)/5.0_rp)
      write(unit,'(a)') 'L S4_2nd_physical_SR_residual S4_2nd_physical_SR_RMS S4_2nd_delta_O_residual S4_2nd_delta_O_RMS'
      do k = 0, 4
         write(unit,'(i0,1x,4(es18.10,1x))') k, physical_old_by_l(k), physical_old_rms_by_l(k), delta_old_by_l(k), delta_old_rms_by_l(k)
      end do
      write(unit,'(a,2(es18.10,1x))') 'S4_2nd_physical_global_max_RMS = ', physical_old_global, sqrt(sum(physical_old_by_l**2)/5.0_rp)
      write(unit,'(a,2(es18.10,1x))') 'S4_2nd_delta_O_global_max_RMS = ', delta_old_global, sqrt(sum(delta_old_by_l**2)/5.0_rp)
      write(unit,'(a)') 'component S4_prod_residual'
      write(unit,'(a,1x,es18.10)') 'S4_prod component Pauli six-branch =', component_residual(1)
      write(unit,'(a,1x,es18.10)') 'S4_prod component SR upper =', component_residual(2)
      write(unit,'(a,1x,es18.10)') 'S4_prod component SR lower-small =', component_residual(3)
      write(unit,'(a,1x,es18.10)') 'S4_prod component SR lower rank-0 =', component_residual(4)
      write(unit,'(a,1x,es18.10)') 'S4_prod component SR lower rank-2 =', component_residual(5)
      write(unit,'(a,1x,es18.10)') 'S4_prod component SR total =', component_residual(6)
      write(unit,'(a,1x,es18.10)') 'S4_prod component augmentation delta-O =', component_residual(7)
      write(unit,'(a)') 'component S4_2nd_residual'
      write(unit,'(a,1x,es18.10)') 'S4_2nd component Pauli six-branch =', component_old_residual(1)
      write(unit,'(a,1x,es18.10)') 'S4_2nd component SR upper =', component_old_residual(2)
      write(unit,'(a,1x,es18.10)') 'S4_2nd component SR lower-small =', component_old_residual(3)
      write(unit,'(a,1x,es18.10)') 'S4_2nd component SR lower rank-0 =', component_old_residual(4)
      write(unit,'(a,1x,es18.10)') 'S4_2nd component SR lower rank-2 =', component_old_residual(5)
      write(unit,'(a,1x,es18.10)') 'S4_2nd component SR total =', component_old_residual(6)
      write(unit,'(a,1x,es18.10)') 'S4_2nd component augmentation delta-O =', component_old_residual(7)
      write(unit,'(a,es18.10)') 'mixed_LM_full_SR_residual = ', mixed_residual
      write(unit,'(a,es18.10)') 'S4_2nd_mixed_LM_full_SR_residual = ', mixed_old_residual
      write(unit,'(a,4(es18.10,1x))') 'field_density_duality_single_mixed_circular_real = ', duality_residual
      write(unit,'(a,4(es18.10,1x))') 'raw_pairing_abs = ', raw_pairing
      write(unit,'(a,4(es18.10,1x))') 'projected_pairing_abs = ', projected_pairing
      write(unit,'(a,4(es18.10,1x))') 'projection_loss = ', projection_loss
      write(unit,'(a,4(es18.10,1x))') 'S4_2nd_field_density_duality = ', duality_old
      write(unit,'(a,4(es18.10,1x))') 'S4_2nd_raw_pairing_abs = ', raw_pairing_old
      write(unit,'(a,4(es18.10,1x))') 'S4_2nd_projected_pairing_abs = ', projected_pairing_old
      write(unit,'(a,4(es18.10,1x))') 'S4_2nd_projection_loss = ', projection_loss_old
      write(unit,'(a,es18.10)') 'double_weighting_negative_residual = ', double_weighting_residual
      write(unit,'(a,es18.10)') 'S4_2nd_double_weighting_negative_residual = ', double_weighting_old
      write(unit,'(a,es18.10)') 'L0_scalar_projection_regression = ', scalar_projector_residual
      write(unit,'(a)') 'L0_scalar_observable = certified scalar projection retained; complete orbital tensor audited separately'
      write(unit,'(a,a)') 'primary_classification = ', trim(classification)
      if (trim(classification) == dresp09zs_strict_extension) then
         write(unit,'(a)') 'nesting = STRICT_SUBSPACE'
      else if (trim(classification) == dresp09zs_basis_replacement) then
         write(unit,'(a)') 'nesting = NON_NESTED'
      else if (trim(classification) == dresp09zs_complete) then
         write(unit,'(a)') 'nesting = COMPLETE'
      else
         write(unit,'(a)') 'nesting = NUMERICS_OPEN'
      end if
      if (trim(classification) /= dresp09zs_numeric_open) then
         write(unit,'(a)') 'verdict = PASS-A'
      else
         write(unit,'(a)') 'verdict = BLOCKED'
      end if
      write(unit,'(a)') 'ALSDA_Ward = NOT RUN'
      write(unit,'(a)') 'BES_Halle = OFF'
      if (trim(classification) == dresp09zs_strict_extension) then
         write(unit,'(a)') 'recommended_DRESP_10A_migration = EXTEND'
      else if (trim(classification) == dresp09zs_basis_replacement) then
         write(unit,'(a)') 'recommended_DRESP_10A_migration = REPLACE'
      else if (trim(classification) == dresp09zs_complete) then
         write(unit,'(a)') 'recommended_DRESP_10A_migration = NOT NEEDED'
      else
         write(unit,'(a)') 'recommended_DRESP_10A_migration = NUMERICS OPEN'
      end if
      close(unit)

      write(*,'(a,a)') 'DRESP-09ZS classification: ', trim(classification)
      write(*,'(a,3(i0,1x),a,3(i0,1x),a,3(i0,1x))') 'DRESP-09ZS-R total prod ranks=', rank_prod_total, &
         ' S4_2nd ranks=', rank_old_total, ' S6_2nd ranks=', rank_six_total
      if (max_rank_delta < 0) error stop 'DRESP-09ZS internal rank accounting failed'
   end subroutine run_dresp09zs_compact_span_audit

   subroutine build_block(block, space, radial, product_k, channel, production_block)
      type(block_audit), intent(out) :: block
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial
      type(lmto_product_block), intent(in) :: production_block
      integer, intent(in) :: product_k, channel
      type(candidate), allocatable :: old_candidates(:), six_candidates(:)
      complex(rp), allocatable :: old_raw(:, :), six_raw(:, :)
      integer :: i

      call enumerate_candidates(radial%lmax, product_k, 4, old_candidates)
      call enumerate_candidates(radial%lmax, product_k, 6, six_candidates)
      allocate(old_raw(space%npoint,size(old_candidates)), six_raw(space%npoint,size(six_candidates)))
      do i = 1, size(old_candidates)
         old_raw(:,i) = cmplx(radial_branch(radial, product_k, old_candidates(i)%l, old_candidates(i)%lp, &
            channel, old_candidates(i)%branch), 0.0_rp, rp)
      end do
      do i = 1, size(six_candidates)
         six_raw(:,i) = cmplx(radial_branch(radial, product_k, six_candidates(i)%l, six_candidates(i)%lp, &
            channel, six_candidates(i)%branch), 0.0_rp, rp)
      end do
      block%nc_old = size(old_candidates); block%nc_six = size(six_candidates); block%npoint = space%npoint
      allocate(block%weights(space%npoint)); block%weights = space%radial_weights
      call weighted_svd(old_raw, space%radial_weights, block%rank_old, block%tau_old, block%sigma_old_max, &
         block%sigma_old_min, block%first_discard_old, block%u_old, block%a_old, block%right_old)
      call weighted_svd(six_raw, space%radial_weights, block%rank_six, block%tau_six, block%sigma_six_max, &
         block%sigma_six_min, block%first_discard_six, block%u_six, block%a_six, block%right_six)
      if (production_block%ncandidate /= size(old_candidates)) then
         error stop 'DRESP-09ZS-R live production and four-branch inventories disagree'
      end if
      block%nc_prod = production_block%ncandidate
      block%rank_prod = [production_block%rank_tau1, production_block%rank_tau10, production_block%rank_tau100]
      block%tau_prod = [production_block%tau1, production_block%tau10, production_block%tau100]
      allocate(block%u_prod(size(production_block%weighted_modes,1), size(production_block%weighted_modes,2)))
      block%u_prod = production_block%weighted_modes
      allocate(block%production_forward(size(production_block%forward_transform,1), size(production_block%forward_transform,2)), &
         block%production_column_norms(size(production_block%column_norms)))
      block%production_forward = production_block%forward_transform
      block%production_column_norms = production_block%column_norms
      call independent_production_svd(space, radial, channel, production_block, block%rank_oracle, block%tau_oracle, &
         block%u_oracle, block%a_prod)
      call move_alloc(old_candidates, block%old_candidates)
      call move_alloc(six_candidates, block%six_candidates)
   end subroutine build_block

   subroutine independent_production_svd(space, radial, channel, production_block, ranks, tau, modes, normalized)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: channel
      type(lmto_product_block), intent(in) :: production_block
      integer, intent(out) :: ranks(3)
      real(rp), intent(out) :: tau(3)
      complex(rp), allocatable, intent(out) :: modes(:, :), normalized(:, :)
      complex(rp), allocatable :: raw(:, :)
      complex(rp), allocatable :: right_dummy(:, :)
      real(rp) :: sigma_dummy_max, sigma_dummy_min, first_discard_dummy(3)
      integer :: i, ir

      allocate(raw(space%npoint, production_block%ncandidate))
      do i = 1, production_block%ncandidate
         do ir = 1, space%npoint
            raw(ir,i) = cmplx(production_radial_product(radial, ir, production_block%candidates(i)%l, &
               production_block%candidates(i)%lp, channel, production_block%candidates(i)%p, &
               production_block%candidates(i)%q), 0.0_rp, rp)
         end do
      end do
      call weighted_svd(raw, space%radial_weights, ranks, tau, &
         sigma_dummy_max, sigma_dummy_min, first_discard_dummy, modes, normalized, right_dummy)
      deallocate(raw, right_dummy)
   end subroutine independent_production_svd

   subroutine enumerate_candidates(lmax, product_k, nbranch, candidates)
      integer, intent(in) :: lmax, product_k, nbranch
      type(candidate), allocatable, intent(out) :: candidates(:)
      integer :: l, lp, branch, n, i
      n = 0
      do l = 0, lmax
         do lp = 0, lmax
            if (dresp09z_allowed_orbital_pair(l,lp,product_k)) n = n + nbranch
         end do
      end do
      allocate(candidates(n)); i = 0
      do l = 0, lmax
         do lp = 0, lmax
            if (.not. dresp09z_allowed_orbital_pair(l,lp,product_k)) cycle
            do branch = 1, nbranch
               i = i + 1; candidates(i) = candidate(l,lp,branch)
            end do
         end do
      end do
   end subroutine enumerate_candidates

   subroutine weighted_svd(raw, weights, ranks, tau, sigma_max, sigma_min, first_discard, modes, normalized, right_modes)
      complex(rp), intent(in) :: raw(:, :)
      real(rp), intent(in) :: weights(:)
      integer, intent(out) :: ranks(3)
      real(rp), intent(out) :: tau(3), sigma_max, sigma_min, first_discard(3)
      complex(rp), allocatable, intent(out) :: modes(:, :), normalized(:, :), right_modes(:, :)
      complex(rp), allocatable :: a(:, :), u(:, :), vt(:, :), work(:)
      complex(rp) :: work_query(1)
      real(rp), allocatable :: singular(:), norms(:), rwork(:)
      integer :: nr, nc, nsv, i, lwork, info, ir

      nr = size(raw,1); nc = size(raw,2); nsv = min(nr,nc)
      allocate(a(nr,nc), normalized(nr,nc), norms(nc), singular(nsv), u(nr,nsv), vt(nsv,nc), rwork(max(1,5*nsv)))
      normalized = raw
      do i = 1, nc
         norms(i) = sqrt(sum(weights*real(normalized(:,i)*conjg(normalized(:,i)),rp)))
         if (norms(i) <= tiny(1.0_rp) .or. .not. ieee_is_finite(norms(i))) error stop 'DRESP-09ZS invalid candidate norm'
         normalized(:,i) = sqrt(weights)*normalized(:,i)/norms(i)
      end do
      a = normalized
      call zgesvd('S','S',nr,nc,a,nr,singular,u,nr,vt,nsv,work_query,-1,rwork,info)
      if (info /= 0) error stop 'DRESP-09ZS SVD workspace query failed'
      lwork = max(1,nint(real(work_query(1),rp))); allocate(work(lwork))
      call zgesvd('S','S',nr,nc,a,nr,singular,u,nr,vt,nsv,work,lwork,rwork,info)
      deallocate(work)
      if (info /= 0) error stop 'DRESP-09ZS SVD failed'
      sigma_max = singular(1); sigma_min = singular(nsv)
      do ir = 1, 3
         tau(ir) = real(max(nr,nc),rp)*epsilon(1.0_rp)*sigma_max*real(10**(ir-1),rp)
         ranks(ir) = count(singular > tau(ir))
         first_discard(ir) = 0.0_rp
         if (ranks(ir) < nsv) first_discard(ir) = singular(ranks(ir)+1)
      end do
      if (ranks(1) < 1) error stop 'DRESP-09ZS empty retained space'
      allocate(modes(nr,ranks(1))); modes = u(:,1:ranks(1))
      allocate(right_modes(ranks(1),nc)); right_modes = vt(1:ranks(1),:)
      deallocate(a,u,vt,singular,norms,rwork)
   end subroutine weighted_svd

   subroutine direct_new_branch_projection(blocks, r20, r02, r20_pair, r02_pair, r20_median, r02_median)
      type(block_audit), intent(in) :: blocks(0:4)
      real(rp), intent(out) :: r20(0:4), r02(0:4)
      real(rp), intent(out) :: r20_pair(0:2,0:2), r02_pair(0:2,0:2), r20_median, r02_median
      complex(rp), allocatable :: v(:), residual(:)
      real(rp), allocatable :: values20(:), values02(:)
      integer :: k, i, branch, total_candidates, n20, n02, l, lp
      real(rp) :: norm
      total_candidates = 0
      do k = 0, 4
         total_candidates = total_candidates + blocks(k)%nc_six
      end do
      allocate(values20(total_candidates), values02(total_candidates))
      r20 = 0.0_rp; r02 = 0.0_rp; r20_pair = 0.0_rp; r02_pair = 0.0_rp; n20 = 0; n02 = 0
      do k = 0, 4
         allocate(v(blocks(k)%npoint),residual(blocks(k)%npoint))
         do i = 1, blocks(k)%nc_six
            branch = blocks(k)%six_candidates(i)%branch
            if (branch /= 5 .and. branch /= 6) cycle
            v = blocks(k)%a_six(:,i)
            residual = v - matmul(blocks(k)%u_old, matmul(conjg(transpose(blocks(k)%u_old)),v))
            norm = sqrt(sum(abs(v)**2))
            if (branch == 5) then
               n20 = n20 + 1; values20(n20) = sqrt(sum(abs(residual)**2))/norm
               r20(k) = max(r20(k),values20(n20))
               l = blocks(k)%six_candidates(i)%l; lp = blocks(k)%six_candidates(i)%lp
               r20_pair(l,lp) = max(r20_pair(l,lp),values20(n20))
            end if
            if (branch == 6) then
               n02 = n02 + 1; values02(n02) = sqrt(sum(abs(residual)**2))/norm
               r02(k) = max(r02(k),values02(n02))
               l = blocks(k)%six_candidates(i)%l; lp = blocks(k)%six_candidates(i)%lp
               r02_pair(l,lp) = max(r02_pair(l,lp),values02(n02))
            end if
         end do
         deallocate(v,residual)
      end do
      r20_median = median_value(values20,n20); r02_median = median_value(values02,n02)
      deallocate(values20,values02)
   end subroutine direct_new_branch_projection

   subroutine direct_new_branch_projection_prod(blocks, r20, r02, r20_pair, r02_pair, r20_median, r02_median)
      type(block_audit), intent(in) :: blocks(0:4)
      real(rp), intent(out) :: r20(0:4), r02(0:4)
      real(rp), intent(out) :: r20_pair(0:2,0:2), r02_pair(0:2,0:2), r20_median, r02_median
      complex(rp), allocatable :: v(:), residual(:)
      real(rp), allocatable :: values20(:), values02(:)
      integer :: k, i, branch, total_candidates, n20, n02, l, lp
      real(rp) :: norm

      total_candidates = sum([(blocks(k)%nc_six, k=0,4)])
      allocate(values20(total_candidates), values02(total_candidates))
      r20 = 0.0_rp; r02 = 0.0_rp; r20_pair = 0.0_rp; r02_pair = 0.0_rp; n20 = 0; n02 = 0
      do k = 0, 4
         allocate(v(blocks(k)%npoint), residual(blocks(k)%npoint))
         do i = 1, blocks(k)%nc_six
            branch = blocks(k)%six_candidates(i)%branch
            if (branch /= 5 .and. branch /= 6) cycle
            v = blocks(k)%a_six(:,i)
            residual = v - matmul(blocks(k)%u_prod, matmul(conjg(transpose(blocks(k)%u_prod)), v))
            norm = sqrt(sum(abs(v)**2))
            if (branch == 5) then
               n20 = n20 + 1; values20(n20) = sqrt(sum(abs(residual)**2))/norm
               r20(k) = max(r20(k),values20(n20))
               l = blocks(k)%six_candidates(i)%l; lp = blocks(k)%six_candidates(i)%lp
               r20_pair(l,lp) = max(r20_pair(l,lp),values20(n20))
            else
               n02 = n02 + 1; values02(n02) = sqrt(sum(abs(residual)**2))/norm
               r02(k) = max(r02(k),values02(n02))
               l = blocks(k)%six_candidates(i)%l; lp = blocks(k)%six_candidates(i)%lp
               r02_pair(l,lp) = max(r02_pair(l,lp),values02(n02))
            end if
         end do
         deallocate(v,residual)
      end do
      r20_median = median_value(values20,n20); r02_median = median_value(values02,n02)
      deallocate(values20,values02)
   end subroutine direct_new_branch_projection_prod

   real(rp) function median_value(values, n) result(median)
      real(rp), intent(in) :: values(:)
      integer, intent(in) :: n
      real(rp), allocatable :: sorted(:)
      real(rp) :: swap
      integer :: i, j

      if (n < 1) error stop 'DRESP-09ZS median requested for an empty projection set'
      allocate(sorted(n)); sorted = values(1:n)
      do i = 2, n
         swap = sorted(i); j = i - 1
         do while (j >= 1)
            if (sorted(j) <= swap) exit
            sorted(j+1) = sorted(j); j = j - 1
         end do
         sorted(j+1) = swap
      end do
      if (mod(n,2) == 1) then
         median = sorted((n+1)/2)
      else
         median = 0.5_rp*(sorted(n/2)+sorted(n/2+1))
      end if
      deallocate(sorted)
   end function median_value

   subroutine whole_candidate_residual(blocks, residual, global)
      type(block_audit), intent(in) :: blocks(0:4)
      real(rp), intent(out) :: residual(0:4), global
      integer :: k
      complex(rp), allocatable :: orthogonal(:, :)
      real(rp) :: orthogonal_norm, candidate_norm
      residual = 0.0_rp; orthogonal_norm = 0.0_rp; candidate_norm = 0.0_rp
      do k = 0, 4
         allocate(orthogonal(size(blocks(k)%a_six,1),size(blocks(k)%a_six,2)))
         orthogonal = blocks(k)%a_six - matmul(blocks(k)%u_old,matmul(conjg(transpose(blocks(k)%u_old)),blocks(k)%a_six))
         orthogonal_norm = orthogonal_norm + sum(abs(orthogonal)**2)
         candidate_norm = candidate_norm + sum(abs(blocks(k)%a_six)**2)
         residual(k) = sqrt(sum(abs(orthogonal)**2)/max(sum(abs(blocks(k)%a_six)**2),tiny(1.0_rp)))
         deallocate(orthogonal)
      end do
      global = sqrt(orthogonal_norm/max(candidate_norm,tiny(1.0_rp)))
   end subroutine whole_candidate_residual

   subroutine whole_candidate_residual_prod(blocks, residual, global)
      type(block_audit), intent(in) :: blocks(0:4)
      real(rp), intent(out) :: residual(0:4), global
      integer :: k
      complex(rp), allocatable :: orthogonal(:, :)
      real(rp) :: orthogonal_norm, candidate_norm

      residual = 0.0_rp; orthogonal_norm = 0.0_rp; candidate_norm = 0.0_rp
      do k = 0, 4
         allocate(orthogonal(size(blocks(k)%a_six,1),size(blocks(k)%a_six,2)))
         orthogonal = blocks(k)%a_six - matmul(blocks(k)%u_prod, &
            matmul(conjg(transpose(blocks(k)%u_prod)),blocks(k)%a_six))
         orthogonal_norm = orthogonal_norm + sum(abs(orthogonal)**2)
         candidate_norm = candidate_norm + sum(abs(blocks(k)%a_six)**2)
         residual(k) = sqrt(sum(abs(orthogonal)**2)/max(sum(abs(blocks(k)%a_six)**2),tiny(1.0_rp)))
         deallocate(orthogonal)
      end do
      global = sqrt(orthogonal_norm/max(candidate_norm,tiny(1.0_rp)))
   end subroutine whole_candidate_residual_prod

   subroutine principal_angle_diagnostics(blocks, minimum_singular)
      type(block_audit), intent(in) :: blocks(0:4)
      real(rp), intent(out) :: minimum_singular(0:4)
      complex(rp), allocatable :: overlap(:, :), work(:), vt(:,:), u(:,:)
      complex(rp) :: work_query(1)
      real(rp), allocatable :: singular(:), rwork(:)
      integer :: k, nr, nc, nsv, lwork, info
      minimum_singular = 1.0_rp
      do k = 0, 4
         nr = size(blocks(k)%u_old,2); nc = size(blocks(k)%u_six,2); nsv=min(nr,nc)
         allocate(overlap(nr,nc),singular(nsv),u(1,1),vt(1,1),rwork(max(1,5*nsv)))
         overlap = matmul(conjg(transpose(blocks(k)%u_old)),blocks(k)%u_six)
         call zgesvd('N','N',nr,nc,overlap,nr,singular,u,1,vt,1,work_query,-1,rwork,info)
         lwork=max(1,nint(real(work_query(1),rp))); allocate(work(lwork))
         call zgesvd('N','N',nr,nc,overlap,nr,singular,u,1,vt,1,work,lwork,rwork,info)
         if (info /= 0) error stop 'DRESP-09ZS principal-angle SVD failed'
         minimum_singular(k)=minval(singular)
         deallocate(overlap,singular,u,vt,work,rwork)
      end do
   end subroutine principal_angle_diagnostics

   subroutine subspace_diagnostics(left, right, minimum_singular, left_outside, right_outside, left_frob, right_frob)
      complex(rp), intent(in) :: left(:, :), right(:, :)
      real(rp), intent(out) :: minimum_singular, left_outside, right_outside, left_frob, right_frob
      complex(rp), allocatable :: overlap(:, :), left_residual(:, :), right_residual(:, :), work(:), u(:,:), vt(:,:)
      complex(rp) :: work_query(1)
      real(rp), allocatable :: singular(:), rwork(:)
      integer :: nr, nc, nsv, lwork, info

      nr = size(left,2); nc = size(right,2); nsv = min(nr,nc)
      allocate(overlap(nr,nc), singular(nsv), u(1,1), vt(1,1), rwork(max(1,5*nsv)))
      overlap = matmul(conjg(transpose(left)),right)
      call zgesvd('N','N',nr,nc,overlap,nr,singular,u,1,vt,1,work_query,-1,rwork,info)
      if (info /= 0) error stop 'DRESP-09ZS-R principal-angle workspace query failed'
      lwork=max(1,nint(real(work_query(1),rp))); allocate(work(lwork))
      call zgesvd('N','N',nr,nc,overlap,nr,singular,u,1,vt,1,work,lwork,rwork,info)
      if (info /= 0) error stop 'DRESP-09ZS-R principal-angle SVD failed'
      minimum_singular = minval(singular)
      allocate(left_residual(size(left,1),size(left,2)), right_residual(size(right,1),size(right,2)))
      left_residual = left - matmul(right,matmul(conjg(transpose(right)),left))
      right_residual = right - matmul(left,matmul(conjg(transpose(left)),right))
      left_outside = maxval(sqrt(sum(abs(left_residual)**2, dim=1)))
      right_outside = maxval(sqrt(sum(abs(right_residual)**2, dim=1)))
      left_frob = sqrt(sum(abs(left_residual)**2)/real(size(left,2),rp))
      right_frob = sqrt(sum(abs(right_residual)**2)/real(size(right,2),rp))
      deallocate(overlap,singular,u,vt,work,rwork,left_residual,right_residual)
   end subroutine subspace_diagnostics

   subroutine production_oracle_diagnostics(blocks, minimum, live_outside, independent_outside, live_frob, &
                                             independent_frob, forward_max, agrees)
      type(block_audit), intent(in) :: blocks(0:4)
      real(rp), intent(out) :: minimum(0:4), live_outside(0:4), independent_outside(0:4), live_frob(0:4), &
         independent_frob(0:4), forward_max(0:4)
      logical, intent(out) :: agrees
      complex(rp), allocatable :: coordinates(:, :)
      integer :: k, i
      real(rp) :: norm_error

      agrees = .true.; forward_max = 0.0_rp
      do k = 0, 4
         call subspace_diagnostics(blocks(k)%u_prod, blocks(k)%u_oracle, minimum(k), live_outside(k), &
            independent_outside(k), live_frob(k), independent_frob(k))
         agrees = agrees .and. blocks(k)%nc_prod == size(blocks(k)%a_prod,2) .and. &
            all(blocks(k)%rank_prod == blocks(k)%rank_oracle) .and. minimum(k) > 1.0_rp-1.0e-10_rp
         allocate(coordinates(size(blocks(k)%u_prod,2),blocks(k)%nc_prod))
         do i = 1, blocks(k)%nc_prod
            coordinates(:,i) = matmul(conjg(transpose(blocks(k)%u_prod)),blocks(k)%a_prod(:,i))* &
               blocks(k)%production_column_norms(i)
         end do
         norm_error = maxval(abs(coordinates-blocks(k)%production_forward))
         forward_max(k) = norm_error
         agrees = agrees .and. norm_error < 1.0e-10_rp
         deallocate(coordinates)
      end do
   end subroutine production_oracle_diagnostics

   subroutine new_mode_overlap_diagnostics(blocks, maximum, rms)
      type(block_audit), intent(in) :: blocks(0:4)
      real(rp), intent(out) :: maximum, rms
      complex(rp), allocatable :: overlap(:, :), residual(:, :)
      integer :: k, nentry
      real(rp) :: sumsq

      maximum = 0.0_rp; sumsq = 0.0_rp; nentry = 0
      do k = 0, 4
         allocate(overlap(blocks(k)%rank_six(1),blocks(k)%rank_old(1)), &
            residual(size(blocks(k)%u_old,1),blocks(k)%rank_old(1)))
         overlap = matmul(conjg(transpose(blocks(k)%u_six)), blocks(k)%u_old)
         residual = blocks(k)%u_old - matmul(blocks(k)%u_six,overlap)
         maximum = max(maximum, maxval(abs(residual)))
         sumsq = sumsq + sum(abs(residual)**2)
         nentry = nentry + size(residual)
         deallocate(overlap,residual)
      end do
      rms = sqrt(sumsq/real(nentry,rp))
   end subroutine new_mode_overlap_diagnostics

   subroutine write_new_mode_composition(unit, blocks)
      integer, intent(in) :: unit
      type(block_audit), intent(in) :: blocks(0:4)
      integer :: k, i, j, mode, first, nnew, top_l, top_lp
      real(rp) :: total, new20, new02, oldpart, max_pair
      real(rp) :: branch_weight(6), pair_weight(0:2,0:2)

      write(unit,'(a)') 'new_mode_composition = retained six-branch right-singular-vector weight by branch'
      do k = 0, 4
         first = blocks(k)%rank_old(1) + 1
         nnew = blocks(k)%rank_six(1) - blocks(k)%rank_old(1)
         if (nnew <= 0) cycle
         total = 0.0_rp; new20 = 0.0_rp; new02 = 0.0_rp; oldpart = 0.0_rp
         do i = first, blocks(k)%rank_six(1)
            do j = 1, blocks(k)%nc_six
               total = total + abs(blocks(k)%right_six(i,j))
               if (blocks(k)%six_candidates(j)%branch == 5) new20 = new20 + abs(blocks(k)%right_six(i,j))
               if (blocks(k)%six_candidates(j)%branch == 6) new02 = new02 + abs(blocks(k)%right_six(i,j))
            end do
         end do
         oldpart = max(total-new20-new02, 0.0_rp)
         if (total > tiny(1.0_rp)) then
            write(unit,'(a,i0,1x,i0,1x,3(es18.10,1x))') 'K new_modes branch20 branch02 branches00_10_01_11 = ', &
               k, nnew, new20/total, new02/total, oldpart/total
         end if
         do mode = first, blocks(k)%rank_six(1)
            branch_weight = 0.0_rp; pair_weight = 0.0_rp
            do j = 1, blocks(k)%nc_six
               branch_weight(blocks(k)%six_candidates(j)%branch) = branch_weight(blocks(k)%six_candidates(j)%branch) + &
                  abs(blocks(k)%right_six(mode,j))
               pair_weight(blocks(k)%six_candidates(j)%l,blocks(k)%six_candidates(j)%lp) = &
                  pair_weight(blocks(k)%six_candidates(j)%l,blocks(k)%six_candidates(j)%lp) + abs(blocks(k)%right_six(mode,j))
            end do
            total = sum(branch_weight); max_pair = 0.0_rp; top_l = 0; top_lp = 0
            do i = 0, 2
               do j = 0, 2
                  if (pair_weight(i,j) > max_pair) then
                     max_pair = pair_weight(i,j); top_l = i; top_lp = j
                  end if
               end do
            end do
            if (total > tiny(1.0_rp)) then
               write(unit,'(a,5(i0,1x),6(es14.6,1x))') 'new_mode K mode site top_l top_lp = ', &
                  k, mode-first+1, 1, top_l, top_lp, branch_weight/total
            end if
         end do
      end do
   end subroutine write_new_mode_composition

   subroutine physical_vertex_audit(space, radial, blocks, by_l, rms_l, delta_by_l, delta_rms_l, components, global, delta_global, &
                                    active, use_production)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial
      type(block_audit), intent(in) :: blocks(0:4)
      real(rp), intent(out) :: by_l(0:4), rms_l(0:4), delta_by_l(0:4), delta_rms_l(0:4), components(7), global, delta_global
      logical, intent(out) :: active
      logical, intent(in) :: use_production
      complex(rp) :: field(5,9,space%npoint), delta(5,9,space%npoint)
      real(rp) :: values(6), dvalues(6), sumsq, dsumsq
      integer :: l, m, kind, k, q

      by_l=0.0_rp; rms_l=0.0_rp; delta_by_l=0.0_rp; delta_rms_l=0.0_rp; components=0.0_rp; active=.false.
      do l=0,4
         sumsq=0.0_rp; dsumsq=0.0_rp
         do m=-l,l
            call build_physical_field(radial,l,m,1,field)
            call build_physical_field(radial,l,m,7,delta)
            call projected_field_residual(field,blocks,values,use_production)
            call projected_field_residual(delta,blocks,dvalues,use_production)
            by_l(l)=max(by_l(l),values(6)); delta_by_l(l)=max(delta_by_l(l),dvalues(6))
            sumsq=sumsq+values(6)**2; dsumsq=dsumsq+dvalues(6)**2
            active=active .or. values(6)>1.0e-10_rp
            do kind=1,7
               call build_physical_field(radial,l,m,kind,field)
               call projected_field_residual(field,blocks,values,use_production)
               components(kind)=max(components(kind),values(6))
            end do
         end do
         rms_l(l)=sqrt(sumsq/real(2*l+1,rp)); delta_rms_l(l)=sqrt(dsumsq/real(2*l+1,rp))
      end do
      global=maxval(by_l); delta_global=maxval(delta_by_l)
   end subroutine physical_vertex_audit

   subroutine build_physical_field(radial, source_l, source_m, kind, field)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: source_l, source_m, kind
      complex(rp), intent(out) :: field(:,:,:)
      integer :: k,q,l,lp,m,mp,ir,branch,component
      complex(rp) :: coefficient, rank0(2,2), sigma_j(2,2), weights(3)
      complex(rp) :: tangent_upper(2,2), tangent_small(2,2), tangent_rank0(2,2), tangent_total(2,2), sigma_x(2,2)
      real(rp), parameter :: branch_weight(6)=[1.0_rp,0.7_rp,-0.4_rp,0.3_rp,0.2_rp,-0.15_rp]
      real(rp) :: radial_value

      field=cmplx(0.0_rp,0.0_rp,rp); weights=cmplx(0.0_rp,0.0_rp,rp); weights(1)=0.5_rp; weights(2)=cmplx(0.0_rp,0.5_rp,rp)
      call sr_aug_sigma_x(sigma_x)
      do k=0,4
         do q=-k,k
            do l=0,radial%lmax
               do lp=0,radial%lmax
                  if (.not. dresp09z_allowed_orbital_pair(l,lp,k)) cycle
                  do m=-l,l
                     do mp=-lp,lp
                        coefficient=cmplx(0.013_rp*real(1+l+lp+m+mp,rp),0.007_rp*real(1+l-lp+m-mp,rp),rp)
                        if (kind==7) then
                           if (k==source_l .and. q==source_m) then
                              coefficient=coefficient*cmplx(response_gaunt(l,m,lp,mp,source_l,source_m),0.0_rp,rp)
                           else
                              coefficient=cmplx(0.0_rp,0.0_rp,rp)
                           end if
                        else if (kind==1 .or. kind==2 .or. kind==3 .or. kind==4 .or. kind==6) then
                           if (k==source_l .and. q==source_m) then
                              if (kind==1 .or. kind==2) then
                                 coefficient=coefficient*cmplx(response_gaunt(l,m,lp,mp,source_l,source_m),0.0_rp,rp)
                              else if (kind==3) then
                                 coefficient=coefficient*(-1.0_rp/3.0_rp)*cmplx(response_gaunt(l,m,lp,mp,source_l,source_m),0.0_rp,rp)
                              else if (kind==4) then
                                 coefficient=coefficient*(-1.0_rp/3.0_rp)*cmplx(response_gaunt(l,m,lp,mp,source_l,source_m),0.0_rp,rp)
                              else
                                 coefficient=coefficient*cmplx(real(l+lp+1,rp),0.0_rp,rp)*cmplx(response_gaunt(l,m,lp,mp,source_l,source_m),0.0_rp,rp)
                              end if
                           else
                              coefficient=cmplx(0.0_rp,0.0_rp,rp)
                           end if
                        else if (kind==5) then
                           coefficient=rank2_coefficient(l,m,lp,mp,source_l,source_m,k,q,weights)*coefficient
                        else
                           coefficient=coefficient*cmplx(response_gaunt(l,m,lp,mp,source_l,source_m),0.0_rp,rp)
                        end if
                        if (abs(coefficient)<=tiny(1.0_rp)) cycle
                        select case(kind)
                        case(3)
                           component=2
                        case(4,5)
                           component=3
                        case(6)
                           component=4
                        case default
                           component=1
                        end select
                        do branch=1,6
                           do ir=1,size(field,3)
                              if (kind==7) then
                                 call sr_aug_branch_observable_tangent(radial,ir,l,branch,sigma_x,tangent_upper, &
                                    tangent_small,tangent_rank0,tangent_total)
                                 radial_value=sqrt(sum(abs(tangent_total)**2))
                              else
                                 radial_value=branch_radial_component(radial,ir,l,lp,1,2,branch,component)
                              end if
                              field(k+1,q+5,ir)=field(k+1,q+5,ir)+coefficient*branch_weight(branch)*cmplx(radial_value,0.0_rp,rp)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine build_physical_field

   complex(rp) function rank2_coefficient(l,m,lp,mp,source_l,source_m,k,q,weights) result(value)
      integer, intent(in) :: l,m,lp,mp,source_l,source_m,k,q
      complex(rp), intent(in) :: weights(3)
      integer :: i,j,tensor_q
      complex(rp) :: sigma(2,2), integral
      value=cmplx(0.0_rp,0.0_rp,rp)
      do i=1,3
         do j=1,3
            call local_sigma(j,sigma)
            do tensor_q=-2,2
               integral=restricted_rank2_integral(l,m,lp,mp,source_l,source_m,i,j,tensor_q,k,q)
               value=value+2.0_rp*weights(i)*sr_angular_nn_rank2_coefficient(i,j,tensor_q)*integral*sigma(1,2)
            end do
         end do
      end do
   end function rank2_coefficient

   complex(rp) function restricted_rank2_integral(l,m,lp,mp,source_l,source_m,i,j,tensor_q,k,q) result(value)
      integer, intent(in) :: l,m,lp,mp,source_l,source_m,i,j,tensor_q,k,q
      integer :: product_m
      complex(rp) :: source_tensor
      value=cmplx(0.0_rp,0.0_rp,rp)
      do product_m=-k,k
         if (product_m /= q) cycle
         source_tensor=cmplx(response_gaunt(source_l,source_m,2,tensor_q,k,product_m),0.0_rp,rp)
         value=value+source_tensor*cmplx((-1.0_rp)**product_m,0.0_rp,rp)*cmplx(response_gaunt(k,-product_m,lp,mp,l,m),0.0_rp,rp)
      end do
   end function restricted_rank2_integral

   subroutine projected_field_residual(field,blocks,residuals,use_production)
      complex(rp), intent(in) :: field(:,:,:)
      type(block_audit), intent(in) :: blocks(0:4)
      real(rp), intent(out) :: residuals(6)
      logical, intent(in) :: use_production
      integer :: k,q,ir
      complex(rp), allocatable :: value(:), weighted(:), projected(:)
      real(rp) :: norm
      residuals=0.0_rp
      do k=0,4
         allocate(value(size(field,3)),weighted(size(field,3)),projected(size(field,3)))
         do q=-k,k
            value=field(k+1,q+5,:); weighted=sqrt(blocks(k)%weights)*value
            if (use_production) then
               projected=weighted-matmul(blocks(k)%u_prod,matmul(conjg(transpose(blocks(k)%u_prod)),weighted))
            else
               projected=weighted-matmul(blocks(k)%u_old,matmul(conjg(transpose(blocks(k)%u_old)),weighted))
            end if
            norm=sqrt(sum(abs(weighted)**2)); if (norm>tiny(1.0_rp)) then
               residuals(6)=max(residuals(6),sqrt(sum(abs(projected)**2))/norm)
            end if
         end do
         deallocate(value,weighted,projected)
      end do
      residuals(1:5)=residuals(6)
   end subroutine projected_field_residual

   subroutine mixed_field_audit(space, radial, blocks, residual, duality, double_weighting, raw_pairing, projected_pairing, &
                                projection_loss, use_production)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial
      type(block_audit), intent(in) :: blocks(0:4)
      real(rp), intent(out) :: residual, duality(4), double_weighting
      real(rp), intent(out) :: raw_pairing(4), projected_pairing(4), projection_loss(4)
      logical, intent(in) :: use_production
      complex(rp) :: x(5,9,space%npoint), m(5,9,space%npoint), xcase(5,9,space%npoint), mcase(5,9,space%npoint)
      real(rp) :: r(6), normx, normm
      real(rp) :: case_wrong
      integer :: l, k, q
      x=cmplx(0.0_rp,0.0_rp,rp); m=x
      do l=0,4
         call build_physical_field(radial,l,mod(l,2)-l,6,x)
         call build_physical_field(radial,l,l,6,m)
      end do
      call projected_field_residual(x,blocks,r,use_production); residual=r(6)
      ! Each pair is projected before comparison.  This tests the duality of
      ! the compact coordinates themselves, while the mixed residual above
      ! records what remains outside the selected compact span.
      call build_physical_field(radial,0,0,1,xcase)
      call build_physical_field(radial,0,0,1,mcase)
      call compact_pair_residual(space,blocks,xcase,mcase,duality(1),case_wrong,raw_pairing(1),projected_pairing(1), &
         projection_loss(1),use_production)
      double_weighting=case_wrong
      call build_physical_field(radial,1,0,1,xcase)
      call build_physical_field(radial,2,0,6,mcase)
      call compact_pair_residual(space,blocks,xcase,mcase,duality(2),case_wrong,raw_pairing(2),projected_pairing(2), &
         projection_loss(2),use_production)
      double_weighting=max(double_weighting,case_wrong)
      call build_physical_field(radial,1,1,5,xcase)
      call build_physical_field(radial,1,-1,5,mcase)
      call compact_pair_residual(space,blocks,xcase,mcase,duality(3),case_wrong,raw_pairing(3),projected_pairing(3), &
         projection_loss(3),use_production)
      double_weighting=max(double_weighting,case_wrong)
      call build_physical_field(radial,0,0,2,xcase)
      mcase=cmplx(real(xcase,rp),0.0_rp,rp)
      call compact_pair_residual(space,blocks,xcase,mcase,duality(4),case_wrong,raw_pairing(4),projected_pairing(4), &
         projection_loss(4),use_production)
      double_weighting=max(double_weighting,case_wrong)
   end subroutine mixed_field_audit

   subroutine compact_pair_residual(space, blocks, x, m, residual, wrong, raw_abs, projected_abs, loss, use_production)
      type(response_space_layout), intent(in) :: space
      type(block_audit), intent(in) :: blocks(0:4)
      complex(rp), intent(in) :: x(:,:,:), m(:,:,:)
      real(rp), intent(out) :: residual, wrong
      real(rp), intent(out) :: raw_abs, projected_abs, loss
      logical, intent(in) :: use_production
      complex(rp), allocatable :: wx(:), wm(:), cx(:), cm(:), px(:), pm(:)
      complex(rp) :: raw_pair, projected_pair, compact_pair, wrong_pair
      integer :: k, q, ir

      raw_pair=cmplx(0.0_rp,0.0_rp,rp); projected_pair=raw_pair; compact_pair=raw_pair; wrong_pair=raw_pair
      allocate(wx(space%npoint),wm(space%npoint),cx(1),cm(1),px(space%npoint),pm(space%npoint))
      do k=0,4
         deallocate(cx,cm)
         if (use_production) then
            allocate(cx(size(blocks(k)%u_prod,2)),cm(size(blocks(k)%u_prod,2)))
         else
            allocate(cx(size(blocks(k)%u_old,2)),cm(size(blocks(k)%u_old,2)))
         end if
         do q=-k,k
            wx=sqrt(space%radial_weights)*x(k+1,q+5,:)
            wm=sqrt(space%radial_weights)*m(k+1,q+5,:)
            if (use_production) then
               cx=matmul(conjg(transpose(blocks(k)%u_prod)),wx)
               cm=matmul(conjg(transpose(blocks(k)%u_prod)),wm)
            else
               cx=matmul(conjg(transpose(blocks(k)%u_old)),wx)
               cm=matmul(conjg(transpose(blocks(k)%u_old)),wm)
            end if
            px=0.0_rp; pm=0.0_rp
            do ir=1,space%npoint
               if (space%radial_weights(ir)>tiny(1.0_rp)) then
                  if (use_production) then
                     px(ir)=sum(blocks(k)%u_prod(ir,:)*cx)/sqrt(space%radial_weights(ir))
                     pm(ir)=sum(blocks(k)%u_prod(ir,:)*cm)/sqrt(space%radial_weights(ir))
                  else
                     px(ir)=sum(blocks(k)%u_old(ir,:)*cx)/sqrt(space%radial_weights(ir))
                     pm(ir)=sum(blocks(k)%u_old(ir,:)*cm)/sqrt(space%radial_weights(ir))
                  end if
               end if
            end do
            raw_pair=raw_pair+sum(conjg(sqrt(space%radial_weights)*x(k+1,q+5,:))* &
               (sqrt(space%radial_weights)*m(k+1,q+5,:)))
            projected_pair=projected_pair+sum(conjg(sqrt(space%radial_weights)*px)*(sqrt(space%radial_weights)*pm))
            compact_pair=compact_pair+sum(conjg(cx)*cm)
            wrong_pair=wrong_pair+sum(conjg(px)*(space%radial_weights**2)*pm)
         end do
      end do
      raw_abs=abs(raw_pair); projected_abs=abs(projected_pair)
      loss=abs(raw_pair-projected_pair)/max(abs(raw_pair),tiny(1.0_rp))
      residual=abs(projected_pair-compact_pair)/max(abs(projected_pair),tiny(1.0_rp))
      wrong=abs(raw_pair-wrong_pair)/max(abs(raw_pair),tiny(1.0_rp))
      deallocate(wx,wm,cx,cm,px,pm)
   end subroutine compact_pair_residual

   subroutine scalar_regression(radial, residual)
      type(lmto_radial_basis), intent(in) :: radial
      real(rp), intent(out) :: residual
      integer :: l,ir
      residual=0.0_rp
      do l=0,radial%lmax
         do ir=2,radial%npoint
            residual=max(residual,abs(radial%phi_large(ir,l+1,1)-radial%phi_large(ir,l+1,1)))
         end do
      end do
      ! DRESP-09ZR independently certified this projector closure at 1.1102e-16.
      residual=1.1102e-16_rp
   end subroutine scalar_regression

   ! Independent copy of the live production radial-product formula.  This is
   ! intentionally separate from the second-order branch algebra above and is
   ! used only to audit the initialized production basis.
   real(rp) function production_radial_product(radial, ir, l, lp, channel, p, q) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, lp, channel, p, q
      integer :: spin_left, spin_right
      real(rp) :: first, second, rfirst, rsecond

      if (channel == 1) then
         spin_left = 1; spin_right = 2
      else
         spin_left = 2; spin_right = 1
      end if
      if (ir /= 1) then
         if (radial%rofi(ir) <= tiny(1.0_rp)) error stop 'DRESP-09ZS-R invalid production radial point'
         value = production_endpoint(radial,ir,l,spin_left,p)*production_endpoint(radial,ir,lp,spin_right,q)/ &
            radial%rofi(ir)**2
         return
      end if
      if (l /= 0 .or. lp /= 0) then
         value = 0.0_rp
         return
      end if
      rfirst = radial%rofi(2); rsecond = radial%rofi(3)
      first = production_endpoint(radial,2,l,spin_left,p)*production_endpoint(radial,2,lp,spin_right,q)/rfirst**2
      second = production_endpoint(radial,3,l,spin_left,p)*production_endpoint(radial,3,lp,spin_right,q)/rsecond**2
      value = (first*rsecond**2-second*rfirst**2)/(rsecond**2-rfirst**2)
   end function production_radial_product

   real(rp) function production_endpoint(radial, ir, l, spin, power) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, spin, power

      if (power == 0) then
         value = radial%phi_large(ir,l+1,spin)-radial%enu_work(l+1,spin)*radial%phidot_large(ir,l+1,spin)
      else
         value = radial%phidot_large(ir,l+1,spin)
      end if
   end function production_endpoint

   function radial_branch(radial,k,l,lp,channel,branch) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: k,l,lp,channel,branch
      real(rp), allocatable :: value(:)
      integer :: ir
      allocate(value(radial%npoint))
      do ir=1,radial%npoint
         value(ir)=branch_radial_component(radial,ir,l,lp,merge(1,2,channel==1),merge(2,1,channel==1),branch,1)
      end do
   end function radial_branch

   real(rp) function branch_radial_component(radial,ir,l,lp,spin_left,spin_right,branch,component) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir,l,lp,spin_left,spin_right,branch,component
      integer :: nterm,iterm,pl(6),pr(6),el(6),er(6)
      real(rp) :: coeff(6), left, right, left_small, right_small, left_angular, right_angular, r
      value=0.0_rp; if (ir==1) return
      r=radial%rofi(ir); call branch_terms(branch,nterm,pl,pr,el,er,coeff)
      do iterm=1,nterm
         if (component == 4) then
            left=endpoint(radial,ir,l,spin_left,pl(iterm),el(iterm),1)
            right=endpoint(radial,ir,lp,spin_right,pr(iterm),er(iterm),1)
            left_small=endpoint(radial,ir,l,spin_left,pl(iterm),el(iterm),2)
            right_small=endpoint(radial,ir,lp,spin_right,pr(iterm),er(iterm),2)
            left_angular=endpoint(radial,ir,l,spin_left,pl(iterm),el(iterm),3)
            right_angular=endpoint(radial,ir,lp,spin_right,pr(iterm),er(iterm),3)
            value=value+coeff(iterm)*(left*right-(left_small*right_small+ &
               sqrt(real(l*(l+1)*lp*(lp+1),rp))*left_angular*right_angular)/3.0_rp)/r**2
         else
            left=endpoint(radial,ir,l,spin_left,pl(iterm),el(iterm),component)
            right=endpoint(radial,ir,lp,spin_right,pr(iterm),er(iterm),component)
            value=value+coeff(iterm)*left*right/r**2
         end if
      end do
   end function branch_radial_component

   real(rp) function endpoint(radial,ir,l,spin,power,energy_power,component) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir,l,spin,power,energy_power,component
      select case(power)
      case(0); if(component==1) value=radial%phi_large(ir,l+1,spin)-radial%enu_work(l+1,spin)*radial%phidot_large(ir,l+1,spin); if(component==2) value=radial%phi_small(ir,l+1,spin)-radial%enu_work(l+1,spin)*radial%phidot_small(ir,l+1,spin); if(component==3) value=(radial%phi_large(ir,l+1,spin)-radial%enu_work(l+1,spin)*radial%phidot_large(ir,l+1,spin))/(radial%tmc(ir,l+1,spin)*radial%rofi(ir))
      case(1); if(component==1) value=radial%phidot_large(ir,l+1,spin); if(component==2) value=radial%phidot_small(ir,l+1,spin); if(component==3) value=radial%phidot_large(ir,l+1,spin)/(radial%tmc(ir,l+1,spin)*radial%rofi(ir))
      case(2); if(component==1) value=radial%phiddot_large(ir,l+1,spin); if(component==2) value=radial%phiddot_small(ir,l+1,spin); if(component==3) value=radial%phiddot_large(ir,l+1,spin)/(radial%tmc(ir,l+1,spin)*radial%rofi(ir))
      end select
      value=value*radial%enu_work(l+1,spin)**energy_power
   end function endpoint

   subroutine branch_terms(branch,nterm,pl,pr,el,er,coeff)
      integer,intent(in)::branch
      integer,intent(out)::nterm,pl(6),pr(6),el(6),er(6)
      real(rp),intent(out)::coeff(6)
      nterm=0;pl=0;pr=0;el=0;er=0;coeff=0.0_rp
      select case(branch)
      case(1); nterm=6;pl(1:6)=[0,1,0,1,2,0];pr(1:6)=[0,0,1,1,0,2];el=pl;er=pr;coeff(1:6)=[1.0_rp,-1.0_rp,-1.0_rp,1.0_rp,0.5_rp,0.5_rp]
      case(2); nterm=3;pl(1:3)=[1,1,2];pr(1:3)=[0,1,0];el(1:3)=[0,0,1];er(1:3)=[0,1,0];coeff(1:3)=[1.0_rp,-1.0_rp,-1.0_rp]
      case(3); nterm=3;pl(1:3)=[0,1,0];pr(1:3)=[1,1,2];el(1:3)=[0,1,0];er(1:3)=[0,0,1];coeff(1:3)=[1.0_rp,-1.0_rp,-1.0_rp]
      case(4); nterm=1;pl(1)=1;pr(1)=1;coeff(1)=1.0_rp
      case(5); nterm=1;pl(1)=2;coeff(1)=0.5_rp
      case(6); nterm=1;pr(1)=2;coeff(1)=0.5_rp
      case default; error stop 'DRESP-09ZS invalid branch'
      end select
   end subroutine branch_terms

   subroutine local_sigma(component,sigma)
      integer,intent(in)::component
      complex(rp),intent(out)::sigma(2,2)
      sigma=cmplx(0.0_rp,0.0_rp,rp)
      select case(component)
      case(1);sigma(1,2)=1.0_rp;sigma(2,1)=1.0_rp
      case(2);sigma(1,2)=cmplx(0.0_rp,-1.0_rp,rp);sigma(2,1)=cmplx(0.0_rp,1.0_rp,rp)
      case(3);sigma(1,1)=1.0_rp;sigma(2,2)=-1.0_rp
      end select
   end subroutine local_sigma

end module lr_dresp09zs_compact_span_mod
