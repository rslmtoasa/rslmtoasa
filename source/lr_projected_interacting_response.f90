!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief DRESP-04 projected Mills/Stoner interaction and site Dyson response.
!>
!> This module is deliberately smaller than the radial/product TDDFT layers.
!> Its coordinate system is the already certified DRESP-01/02 site basis:
!>
!>   accepted no-SOC state -> local Pauli field fit -> U_i=Delta_i/M_i
!>   -> (I-chi0 U) chi = chi0.
!>
!> The fit consumes samples of the Pauli coefficient B_sigma in
!> H = H_0 I + B_sigma sigma_z.  Thus H_up-H_down=2 B_sigma.  Delta in
!> the result is B_sigma, the coefficient multiplying the DRESP-01 Vz
!> convention (up=+1, down=-1), not the full up/down difference.
!>
!> No Goldstone repair, Juelich/LCMM sum rule, ALSDA kernel, force-theorem
!> exchange vertex, or product-space reconstruction belongs here.
!------------------------------------------------------------------------------
module lr_projected_interacting_response_mod

   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use reciprocal_mod, only: reciprocal
   use lr_projected_site_spin_mod, only: projected_site_spin_contract, projected_operator_z
   use tddft_dyson_mod, only: solve_tddft_dyson_frequency, tddft_loss_matrix, &
      tddft_denominator_minimum_magnitude_eigenvalue
   implicit none
   private

   character(len=*), parameter, public :: projected_mills_exact_scalar = 'EXACT_SCALAR'
   character(len=*), parameter, public :: projected_mills_projected_scalar = &
      'PROJECTED_SCALAR_APPROXIMATION'
   character(len=*), parameter, public :: projected_mills_unsupported = 'UNSUPPORTED'
   character(len=*), parameter, public :: projected_mills_convention = &
      'H=H0 I+B_sigma sigma_z; H_up-H_down=2 B_sigma; DRESP-01 Vz=sigma_z'
   character(len=*), parameter, public :: projected_mills_loss_convention = &
      'L=-(chi-chi^dagger)/(2*i*pi); ordinary site-space matrix'

   real(rp), parameter, public :: projected_mills_exact_tolerance = 1.0e-10_rp
   real(rp), parameter, public :: projected_mills_moment_floor = 1.0e-12_rp
   real(rp), parameter, public :: projected_mills_condition_limit = 1.0e12_rp

   !> Accepted-state samples for the local scalarization.  The matrices are
   !> site-major in the supplied coefficient subspace.  A sample is normally
   !> one accepted reciprocal H(k); sample_weights are relative weights and
   !> default to a uniform Frobenius metric.
   type, public :: projected_mills_interaction_request
      character(len=8) :: selector = ''
      character(len=80) :: splitting_convention = projected_mills_convention
      integer :: nsite = 0
      integer :: site_block_size = 0
      complex(rp), allocatable :: actual_pauli_field(:, :, :) ! (basis,basis,sample)
      complex(rp), allocatable :: site_vertices(:, :, :, :) ! (basis,basis,site,sample)
      real(rp), allocatable :: projected_moment(:)
      real(rp), allocatable :: sample_weights(:)
      character(len=512) :: provenance = ''
   end type projected_mills_interaction_request

   type, public :: projected_mills_interaction_result
      character(len=8) :: selector = ''
      character(len=80) :: splitting_convention = projected_mills_convention
      character(len=32) :: classification = projected_mills_unsupported
      character(len=512) :: provenance = ''
      integer :: nsite = 0
      integer :: nbasis = 0
      integer :: nsample = 0
      integer :: rank = 0
      real(rp) :: fit_condition_number = huge(1.0_rp)
      real(rp) :: scalarization_residual = huge(1.0_rp)
      real(rp) :: locality_residual = huge(1.0_rp)
      real(rp) :: splitting_norm = 0.0_rp
      real(rp) :: fit_residual_norm = 0.0_rp
      real(rp) :: coefficient_imaginary_residual = 0.0_rp
      real(rp), allocatable :: projected_moment(:)
      real(rp), allocatable :: projected_splitting(:) ! B_sigma coefficient Delta_i
      real(rp), allocatable :: interaction_U(:)
   end type projected_mills_interaction_result

   !> Site-space projected Dyson request.  bare_chi is directly the
   !> DRESP-02 site x site x frequency object; no product-space input exists.
   type, public :: projected_dyson_request
      character(len=8) :: selector = ''
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = ''
      real(rp), allocatable :: interaction_U(:)
      complex(rp), allocatable :: bare_chi(:, :, :) ! (site,site,frequency)
      character(len=512) :: interaction_provenance = ''
      character(len=512) :: bare_provenance = ''
   end type projected_dyson_request

   type, public :: projected_dyson_result
      character(len=8) :: selector = ''
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = ''
      character(len=512) :: interaction_provenance = ''
      character(len=512) :: bare_provenance = ''
      character(len=256) :: dyson_convention = 'D=I-chi0*U; solve D*chi=chi0 with certified LAPACK zgesv'
      character(len=256) :: loss_convention = projected_mills_loss_convention
      character(len=128) :: status = 'not evaluated'
      real(rp), allocatable :: interaction_U(:)
      complex(rp), allocatable :: bare_chi(:, :, :)
      complex(rp), allocatable :: interaction(:, :)
      complex(rp), allocatable :: denominator(:, :, :)
      complex(rp), allocatable :: enhanced_chi(:, :, :)
      complex(rp), allocatable :: loss_matrix(:, :, :)
      real(rp), allocatable :: denominator_min_singular_value(:)
      real(rp), allocatable :: denominator_max_singular_value(:)
      real(rp), allocatable :: condition_number(:)
      real(rp), allocatable :: minimum_magnitude_eigenvalue(:)
      real(rp), allocatable :: loss_trace(:)
      real(rp), allocatable :: minus_im_trace_over_pi(:)
      real(rp), allocatable :: dyson_residual(:)
      real(rp), allocatable :: dyson_residual_relative(:)
      real(rp), allocatable :: dyson_residual_infinity(:)
      integer, allocatable :: solve_info(:)
   end type projected_dyson_result

   public :: evaluate_projected_mills_interaction
   public :: evaluate_projected_mills_from_reciprocal
   public :: evaluate_projected_dyson

contains

   subroutine evaluate_projected_mills_interaction(request, result)
      type(projected_mills_interaction_request), intent(in) :: request
      type(projected_mills_interaction_result), intent(out) :: result

      complex(rp), allocatable :: design(:, :), rhs(:, :), solution(:), work(:)
      complex(rp) :: work_query(1)
      real(rp), allocatable :: singular_values(:), rwork(:), sample_weights(:)
      real(rp) :: weight, norm2, residual2, local2, largest_sv, smallest_sv, rcond, scale
      integer :: n, nsite, nsample, nrow, row, i, j, isample, isite, lwork, info
      integer :: rank, block_i, block_j
      external :: zgelss

      call validate_mills_request(request, n, nsite, nsample)
      allocate(result%projected_moment(nsite), result%projected_splitting(nsite), result%interaction_U(nsite))
      result%selector = trim(request%selector)
      result%splitting_convention = trim(request%splitting_convention)
      result%provenance = trim(request%provenance)
      result%nsite = nsite
      result%nbasis = n
      result%nsample = nsample
      result%projected_moment = request%projected_moment
      result%projected_splitting = 0.0_rp
      result%interaction_U = 0.0_rp

      allocate(sample_weights(nsample))
      if (allocated(request%sample_weights)) then
         sample_weights = request%sample_weights
      else
         sample_weights = 1.0_rp
      end if
      if (any(sample_weights <= 0.0_rp) .or. any(.not. ieee_is_finite(sample_weights))) then
         error stop 'DRESP-04 Mills: sample weights must be finite and positive'
      end if

      nrow = n*n*nsample
      allocate(design(nrow, nsite), rhs(nrow, 1), solution(nsite), singular_values(min(nrow, nsite)), &
         rwork(max(1, 5*min(nrow, nsite))))
      design = cmplx(0.0_rp, 0.0_rp, rp)
      rhs = cmplx(0.0_rp, 0.0_rp, rp)
      row = 0
      norm2 = 0.0_rp
      local2 = 0.0_rp
      do isample = 1, nsample
         weight = sqrt(sample_weights(isample))
         do j = 1, n
            do i = 1, n
               row = row + 1
               design(row, :) = weight*request%site_vertices(i, j, :, isample)
               rhs(row, 1) = weight*request%actual_pauli_field(i, j, isample)
               norm2 = norm2 + sample_weights(isample)*abs(request%actual_pauli_field(i, j, isample))**2
               block_i = (i - 1)/request%site_block_size + 1
               block_j = (j - 1)/request%site_block_size + 1
               if (block_i /= block_j) local2 = local2 + sample_weights(isample)* &
                  abs(request%actual_pauli_field(i, j, isample))**2
            end do
         end do
      end do
      result%splitting_norm = sqrt(norm2)
      if (norm2 > 0.0_rp) then
         result%locality_residual = sqrt(local2/norm2)
      else
         result%locality_residual = 0.0_rp
      end if

      call zgelss(nrow, nsite, 1, design, nrow, rhs, nrow, singular_values, -1.0_rp, rank, &
         work_query, -1, rwork, info)
      if (info /= 0) error stop 'DRESP-04 Mills: zgelss workspace query failed'
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      rcond = -1.0_rp
      call zgelss(nrow, nsite, 1, design, nrow, rhs, nrow, singular_values, rcond, rank, work, lwork, rwork, info)
      if (info /= 0) then
         result%rank = 0
         result%classification = projected_mills_unsupported
         deallocate(design, rhs, solution, singular_values, work, rwork, sample_weights)
         return
      end if
      result%rank = rank
      solution = rhs(1:nsite, 1)
      largest_sv = maxval(singular_values)
      smallest_sv = huge(1.0_rp)
      do i = 1, size(singular_values)
         if (singular_values(i) > 0.0_rp) smallest_sv = min(smallest_sv, singular_values(i))
      end do
      if (largest_sv > 0.0_rp .and. smallest_sv < huge(1.0_rp)) then
         result%fit_condition_number = largest_sv/smallest_sv
      else
         result%fit_condition_number = huge(1.0_rp)
      end if
      result%coefficient_imaginary_residual = maxval(abs(aimag(solution)))/max(maxval(abs(solution)), tiny(1.0_rp))
      do isite = 1, nsite
         result%projected_splitting(isite) = real(solution(isite), rp)
      end do

      residual2 = 0.0_rp
      row = 0
      do isample = 1, nsample
         weight = sqrt(sample_weights(isample))
         do j = 1, n
            do i = 1, n
               row = row + 1
               residual2 = residual2 + abs(weight*(request%actual_pauli_field(i, j, isample) - &
                  sum(request%site_vertices(i, j, :, isample)*solution)))**2
            end do
         end do
      end do
      result%fit_residual_norm = sqrt(residual2)
      result%scalarization_residual = result%fit_residual_norm/max(result%splitting_norm, tiny(1.0_rp))

      do isite = 1, nsite
         if (.not. ieee_is_finite(result%projected_moment(isite)) .or. &
             abs(result%projected_moment(isite)) <= projected_mills_moment_floor) then
            result%classification = projected_mills_unsupported
            deallocate(design, rhs, solution, singular_values, work, rwork, sample_weights)
            return
         end if
         result%interaction_U(isite) = result%projected_splitting(isite)/result%projected_moment(isite)
         if (.not. ieee_is_finite(result%interaction_U(isite))) then
            result%classification = projected_mills_unsupported
            deallocate(design, rhs, solution, singular_values, work, rwork, sample_weights)
            return
         end if
      end do
      scale = max(result%splitting_norm, tiny(1.0_rp))
      if (rank < nsite .or. result%fit_condition_number > projected_mills_condition_limit .or. &
          result%coefficient_imaginary_residual > 1.0e-8_rp) then
         result%classification = projected_mills_unsupported
      else if (result%scalarization_residual <= projected_mills_exact_tolerance .and. &
               result%locality_residual <= projected_mills_exact_tolerance) then
         result%classification = projected_mills_exact_scalar
      else
         result%classification = projected_mills_projected_scalar
      end if
      if (.not. ieee_is_finite(scale)) result%classification = projected_mills_unsupported
      deallocate(design, rhs, solution, singular_values, work, rwork, sample_weights)
   end subroutine evaluate_projected_mills_interaction

   !> Build DRESP-04 samples from the accepted orthogonal, collinear,
   !> no-SOC reciprocal Hamiltonian.  The fit uses B_sigma=(H_up-H_down)/2
   !> on the selected orbital coefficient subspace and Vz from DRESP-01's
   !> direct operator contract at the accepted Fermi energy.
   subroutine evaluate_projected_mills_from_reciprocal(contract, radial_bases, reciprocal_obj, energy, moment, result)
      type(projected_site_spin_contract), intent(in) :: contract
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(reciprocal), intent(in) :: reciprocal_obj
      real(rp), intent(in) :: energy, moment(:)
      type(projected_mills_interaction_result), intent(out) :: result

      complex(rp), allocatable :: vz_full(:, :), fields(:, :, :), vertices(:, :, :, :)
      integer, allocatable :: selected_orbitals(:)
      integer :: norb, nsite, nk, nselected, nbasis, nmat, site, iorb, jorb, k
      integer :: a, b, ia, ib, up_i, up_j, down_i, down_j, site_i, site_j
      type(projected_mills_interaction_request) :: request

      if (.not. allocated(reciprocal_obj%hk_bulk)) error stop 'DRESP-04 Mills: accepted hk_bulk is unavailable'
      nsite = contract%nsite
      norb = (contract%orbital_lmax + 1)**2
      nmat = size(reciprocal_obj%hk_bulk, 1)
      nk = size(reciprocal_obj%hk_bulk, 3)
      if (size(radial_bases) /= nsite .or. nmat /= 2*nsite*norb .or. size(moment) /= nsite) then
         error stop 'DRESP-04 Mills: reciprocal/site/radial dimensions are inconsistent'
      end if
      if (allocated(reciprocal_obj%hk_so)) then
         if (maxval(abs(reciprocal_obj%hk_so)) > 1.0e-12_rp) then
            error stop 'DRESP-04 Mills: SOC Hamiltonian is outside the no-SOC contract'
         end if
      end if
      allocate(selected_orbitals(norb))
      nselected = 0
      do iorb = 1, norb
         if (contract%selected_l(lmto_orbital_l(iorb))) then
            nselected = nselected + 1
            selected_orbitals(nselected) = iorb
         end if
      end do
      nbasis = nsite*nselected
      allocate(vz_full(nmat, nmat), fields(nbasis, nbasis, nk), vertices(nbasis, nbasis, nsite, nk))
      call contract%direct_operator_matrix(radial_bases, energy, projected_operator_z, vz_full)
      fields = cmplx(0.0_rp, 0.0_rp, rp)
      vertices = cmplx(0.0_rp, 0.0_rp, rp)
      do k = 1, nk
         ! Carry the complete selected site x site field.  The local vertex
         ! remains block diagonal, so a nonlocal Hamiltonian field is visible
         ! through locality_residual instead of being silently discarded.
         do site_i = 1, nsite
            do a = 1, nselected
               iorb = selected_orbitals(a)
               ia = (site_i - 1)*nselected + a
               up_i = (site_i - 1)*2*norb + iorb
               down_i = up_i + norb
               do site_j = 1, nsite
                  do b = 1, nselected
                     jorb = selected_orbitals(b)
                     ib = (site_j - 1)*nselected + b
                     up_j = (site_j - 1)*2*norb + jorb
                     down_j = up_j + norb
                     fields(ia, ib, k) = 0.5_rp*(reciprocal_obj%hk_bulk(up_i, up_j, k) - &
                        reciprocal_obj%hk_bulk(down_i, down_j, k))
                     if (site_i == site_j) vertices(ia, ib, site_i, k) = vz_full(up_i, up_j)
                  end do
               end do
            end do
         end do
      end do
      request%selector = contract%selector
      request%splitting_convention = projected_mills_convention
      request%nsite = nsite
      request%site_block_size = nselected
      request%actual_pauli_field = fields
      request%site_vertices = vertices
      request%projected_moment = moment
      if (allocated(reciprocal_obj%k_weights) .and. size(reciprocal_obj%k_weights) == nk) then
         request%sample_weights = reciprocal_obj%k_weights
      end if
      request%provenance = 'accepted reciprocal hk_bulk; B_sigma=(H_up-H_down)/2; DRESP-01 direct Vz at EF'
      call evaluate_projected_mills_interaction(request, result)
      deallocate(vz_full, fields, vertices, selected_orbitals)
   end subroutine evaluate_projected_mills_from_reciprocal

   subroutine evaluate_projected_dyson(request, result)
      type(projected_dyson_request), intent(in) :: request
      type(projected_dyson_result), intent(out) :: result

      complex(rp), allocatable :: residual(:, :), chi(:, :), denominator(:, :), loss(:, :)
      integer :: nsite, nw, iw, i, info
      real(rp) :: dnorm, rnorm
      complex(rp) :: trace_value

      call validate_dyson_request(request, nsite, nw)
      allocate(result%frequencies(nw), result%interaction_U(nsite), result%bare_chi(nsite, nsite, nw), &
         result%interaction(nsite, nsite), result%denominator(nsite, nsite, nw), &
         result%enhanced_chi(nsite, nsite, nw), result%loss_matrix(nsite, nsite, nw), &
         result%denominator_min_singular_value(nw), result%denominator_max_singular_value(nw), &
         result%condition_number(nw), result%minimum_magnitude_eigenvalue(nw), result%loss_trace(nw), &
         result%minus_im_trace_over_pi(nw), result%dyson_residual(nw), result%dyson_residual_relative(nw), &
         result%dyson_residual_infinity(nw), result%solve_info(nw))
      result%selector = trim(request%selector)
      result%q = request%q
      result%frequencies = request%frequencies
      result%eta = request%eta
      result%channel = trim(request%channel)
      result%interaction_provenance = trim(request%interaction_provenance)
      result%bare_provenance = trim(request%bare_provenance)
      result%interaction_U = request%interaction_U
      result%bare_chi = request%bare_chi
      result%interaction = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, nsite
         result%interaction(i, i) = cmplx(request%interaction_U(i), 0.0_rp, rp)
      end do
      allocate(residual(nsite, nsite), chi(nsite, nsite), denominator(nsite, nsite), loss(nsite, nsite))
      do iw = 1, nw
         call solve_tddft_dyson_frequency(request%bare_chi(:, :, iw), result%interaction, chi, denominator, info, &
            result%condition_number(iw), result%denominator_min_singular_value(iw), result%denominator_max_singular_value(iw))
         result%solve_info(iw) = info
         result%denominator(:, :, iw) = denominator
         result%enhanced_chi(:, :, iw) = chi
         call tddft_denominator_minimum_magnitude_eigenvalue(denominator, result%minimum_magnitude_eigenvalue(iw))
         if (info == 0) then
            loss = tddft_loss_matrix(chi)
         else
            loss = cmplx(0.0_rp, 0.0_rp, rp)
         end if
         result%loss_matrix(:, :, iw) = loss
         trace_value = cmplx(0.0_rp, 0.0_rp, rp)
         do i = 1, nsite
            trace_value = trace_value + loss(i, i)
         end do
         result%loss_trace(iw) = real(trace_value, rp)
         trace_value = cmplx(0.0_rp, 0.0_rp, rp)
         do i = 1, nsite
            trace_value = trace_value + chi(i, i)
         end do
         result%minus_im_trace_over_pi(iw) = -aimag(trace_value)/acos(-1.0_rp)
         residual = matmul(denominator, chi) - request%bare_chi(:, :, iw)
         dnorm = sqrt(sum(abs(request%bare_chi(:, :, iw))**2))
         rnorm = sqrt(sum(abs(residual)**2))
         result%dyson_residual(iw) = rnorm
         result%dyson_residual_relative(iw) = rnorm/max(dnorm, tiny(1.0_rp))
         result%dyson_residual_infinity(iw) = maxval(abs(residual))
      end do
      if (any(result%solve_info /= 0)) then
         result%status = 'BLOCKED: one or more raw Dyson denominators are singular at requested precision'
      else
         result%status = 'PASS: projected site Dyson and loss diagnostics'
      end if
      deallocate(residual, chi, denominator, loss)
   end subroutine evaluate_projected_dyson

   subroutine validate_mills_request(request, n, nsite, nsample)
      type(projected_mills_interaction_request), intent(in) :: request
      integer, intent(out) :: n, nsite, nsample
      n = size(request%actual_pauli_field, 1)
      nsite = request%nsite
      nsample = size(request%actual_pauli_field, 3)
      if (n < 1 .or. nsite < 1 .or. nsample < 1 .or. request%site_block_size < 1 .or. &
          n /= nsite*request%site_block_size .or. any(shape(request%actual_pauli_field) /= [n, n, nsample]) .or. &
          any(shape(request%site_vertices) /= [n, n, nsite, nsample]) .or. size(request%projected_moment) /= nsite) then
         error stop 'DRESP-04 Mills: invalid scalarization dimensions'
      end if
      if (len_trim(request%selector) == 0) error stop 'DRESP-04 Mills: selector provenance is required'
      if (any(.not. ieee_is_finite(real(request%actual_pauli_field, rp))) .or. &
          any(.not. ieee_is_finite(aimag(request%actual_pauli_field)))) then
         error stop 'DRESP-04 Mills: actual splitting contains non-finite values'
      end if
   end subroutine validate_mills_request

   subroutine validate_dyson_request(request, nsite, nw)
      type(projected_dyson_request), intent(in) :: request
      integer, intent(out) :: nsite, nw
      nsite = size(request%interaction_U)
      nw = size(request%frequencies)
      if (nsite < 1 .or. nw < 1 .or. any(shape(request%bare_chi) /= [nsite, nsite, nw]) .or. &
          any(.not. ieee_is_finite(request%interaction_U)) .or. request%eta <= 0.0_rp) then
         error stop 'DRESP-04 Dyson: invalid site response dimensions or eta'
      end if
      if (len_trim(request%selector) == 0 .or. len_trim(request%channel) == 0) then
         error stop 'DRESP-04 Dyson: selector and circular-channel provenance are required'
      end if
   end subroutine validate_dyson_request

end module lr_projected_interacting_response_mod
