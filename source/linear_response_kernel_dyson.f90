!------------------------------------------------------------------------------
! Linear-response kernel, compact interaction, projected services, and Dyson
! implementations moved under linear_response_mod (LR-REF-03c).
!------------------------------------------------------------------------------
submodule (linear_response_mod) linear_response_kernel_dyson
   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use math_mod, only: pi
   use radial_ground_state_mod, only: radial_ground_state, radial_xc_provenance
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use reciprocal_mod, only: reciprocal
   implicit none

   interface
      subroutine zgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info)
         import :: rp
         integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
         complex(rp), intent(inout) :: a(lda, *), b(ldb, *), work(*)
         real(rp), intent(out) :: s(*)
         real(rp), intent(in) :: rcond
         integer, intent(out) :: rank, info
         real(rp), intent(inout) :: rwork(*)
      end subroutine zgelss

      subroutine dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
         import :: rp
         integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
         real(rp), intent(inout) :: a(lda, *), b(ldb, *), work(*)
         real(rp), intent(out) :: s(*)
         real(rp), intent(in) :: rcond
         integer, intent(out) :: rank, info
      end subroutine dgelss
   end interface

contains

   ! --- from lr_alsda_kernel_mod ---
   !> Evaluate the request-form direct ALSDA kernel.
   module subroutine evaluate_lr_alsda_kernel_request(request, result)
      type(lr_alsda_kernel_request), intent(in) :: request
      type(lr_alsda_kernel_result), intent(out) :: result

      if (request%production_contract) then
         if (trim(request%magnetization_kind) /= lr_kxc_magnetization_kind_pauli_accepted .or. &
             trim(request%magnetization_source) /= lr_kxc_magnetization_source_pauli_accepted) then
            error stop 'evaluate_lr_alsda_kernel: production direct ALSDA requires certified PAULI_ACCEPTED magnetization provenance'
         end if
      end if
      call evaluate_lr_alsda_kernel_explicit(request%response_space, request%ground_states, &
         request%pauli_magnetization, result, request%requested_functional, request%requested_backend, &
         request%requested_txc, request%magnetization_label, request%low_m_diagnostic_relative)
      if (request%production_contract) then
         result%magnetization_kind = trim(request%magnetization_kind)
         result%magnetization_source = trim(request%magnetization_source)
      end if
   end subroutine evaluate_lr_alsda_kernel_request

   !> Explicit form retained for callers that keep the three contracts apart.
   module subroutine evaluate_lr_alsda_kernel_explicit(space, ground_states, pauli_magnetization, result, &
                                                requested_functional, requested_backend, requested_txc, &
                                                magnetization_label, low_m_diagnostic_relative)
      type(response_space_layout), intent(in) :: space
      type(radial_ground_state), intent(in) :: ground_states(:)
      real(rp), intent(in) :: pauli_magnetization(:, :)
      type(lr_alsda_kernel_result), intent(out) :: result
      character(len=*), intent(in), optional :: requested_functional
      character(len=*), intent(in), optional :: requested_backend
      integer, intent(in), optional :: requested_txc
      character(len=*), intent(in), optional :: magnetization_label
      real(rp), intent(in), optional :: low_m_diagnostic_relative

      integer :: isite, ir
      real(rp) :: bxc_sigma, magnetization, scale
      type(radial_xc_provenance) :: reference_provenance
      logical :: have_reference_provenance
      character(len=512) :: effective_functional
      character(len=32) :: effective_backend, effective_magnetization_label
      integer :: effective_txc
      real(rp) :: effective_low_m_diagnostic_relative

      effective_functional = ''
      effective_backend = ''
      effective_txc = -1
      effective_magnetization_label = lr_kxc_magnetization_pauli
      effective_low_m_diagnostic_relative = lr_kxc_default_low_m_relative
      if (present(requested_functional)) effective_functional = requested_functional
      if (present(requested_backend)) effective_backend = requested_backend
      if (present(requested_txc)) effective_txc = requested_txc
      if (present(magnetization_label)) effective_magnetization_label = magnetization_label
      if (present(low_m_diagnostic_relative)) effective_low_m_diagnostic_relative = low_m_diagnostic_relative

      call require_initialized_space(space)
      if (trim(effective_magnetization_label) /= lr_kxc_magnetization_pauli) then
         error stop 'evaluate_lr_alsda_kernel: only the LR-03 Pauli-projected magnetization is supported'
      end if
      if (effective_low_m_diagnostic_relative < 0.0_rp) then
         error stop 'evaluate_lr_alsda_kernel: low-m diagnostic threshold must be nonnegative'
      end if

      have_reference_provenance = .false.
      do isite = 1, space%nsite
         call validate_ground_state(space, ground_states(isite), isite)
         if (.not. have_reference_provenance) then
            reference_provenance = ground_states(isite)%xc_provenance
            have_reference_provenance = .true.
         else if (.not. same_provenance(reference_provenance, ground_states(isite)%xc_provenance)) then
            error stop 'evaluate_lr_alsda_kernel: XC provenance differs between response sites'
         end if
      end do
      call validate_requested_provenance(reference_provenance, effective_functional, effective_backend, effective_txc)

      allocate(result%pointwise_kernel(space%nsite, space%npoint), result%canonical_operator(space%ndim, space%ndim))
      result%pointwise_kernel = 0.0_rp
      result%canonical_operator = cmplx(0.0_rp, 0.0_rp, rp)
      result%low_m%diagnostic_relative_threshold = effective_low_m_diagnostic_relative
      result%magnetization_label = trim(effective_magnetization_label)
      result%magnetization_kind = ''
      result%magnetization_source = ''
      result%xc_provenance = reference_provenance

      scale = maxval(abs(pauli_magnetization))
      if (scale <= 0.0_rp) then
         error stop 'evaluate_lr_alsda_kernel: Pauli magnetization is identically zero'
      end if
      result%low_m%max_abs_magnetization = scale

      do isite = 1, space%nsite
         do ir = 1, space%npoint
            magnetization = pauli_magnetization(isite, ir)
            if (space%radial_weights(ir) > 0.0_rp) then
               result%low_m%active_radial_points = result%low_m%active_radial_points + 1
               result%low_m%min_abs_magnetization = min(result%low_m%min_abs_magnetization, abs(magnetization))
               result%low_m%min_relative_abs_magnetization = min(result%low_m%min_relative_abs_magnetization, &
                  abs(magnetization)/scale)
               if (abs(magnetization) <= effective_low_m_diagnostic_relative*scale) then
                  result%low_m%low_m_points = result%low_m%low_m_points + 1
               end if
               if (magnetization == 0.0_rp) then
                  result%low_m%exact_zero_active_points = result%low_m%exact_zero_active_points + 1
                  error stop 'evaluate_lr_alsda_kernel: zero Pauli magnetization on a positive-measure radial point; capability BLOCKED'
               end if
            else if (magnetization == 0.0_rp) then
               ! The origin is an exact null-measure extension in LR-04.  Its
               ! point value cannot affect any physical response contraction.
               result%low_m%null_measure_zero_points = result%low_m%null_measure_zero_points + 1
               result%origin_null_measure_extension = .true.
            end if

            if (magnetization /= 0.0_rp) then
               bxc_sigma = 0.5_rp*(ground_states(isite)%vxc_up(ir) - ground_states(isite)%vxc_down(ir))
               result%pointwise_kernel(isite, ir) = bxc_sigma/magnetization
            else
               result%pointwise_kernel(isite, ir) = 0.0_rp
            end if
         end do
      end do

      ! A spherical local scalar is diagonal in site, (L,M), radial point, and
      ! channel.  LR-04 owns the canonical local representation and its
      ! origin treatment; this call deliberately avoids a guessed 1/W delta.
      call response_local_operator(space, result%pointwise_kernel, result%canonical_operator)

      result%units = lr_kxc_units
      result%response_representation = lr_kxc_representation
   end subroutine evaluate_lr_alsda_kernel_explicit

   !> Check the exact XC provenance requested by a response consumer.
   !>
   !> Blank request fields mean "use the accepted snapshot provenance"; a
   !> nonblank field is an exact contract and a mismatch is fatal in the
   !> evaluator rather than silently routed to another functional.
   module logical function lr_alsda_provenance_matches(accepted, requested_functional, requested_backend, requested_txc) &
      result(matches)
      type(radial_xc_provenance), intent(in) :: accepted
      character(len=*), intent(in) :: requested_functional, requested_backend
      integer, intent(in) :: requested_txc

      matches = .true.
      if (len_trim(requested_functional) > 0) matches = matches .and. &
         trim(accepted%functional_name) == trim(requested_functional)
      if (len_trim(requested_backend) > 0) matches = matches .and. &
         trim(accepted%backend_name) == trim(requested_backend)
      if (requested_txc >= 0) matches = matches .and. accepted%txc == requested_txc
   end function lr_alsda_provenance_matches

   !> Raw static Goldstone/Ward diagnostic for the direct ALSDA route.
   !>
   !> The supplied susceptibility is applied unchanged to the field generated
   !> by the supplied direct kernel and rigid vector.  No matrix, eigenvalue,
   !> or kernel is modified.  The overlap is the metric overlap of the raw
   !> response with the target rigid vector.
   module subroutine evaluate_lr_alsda_static_residual(space, static_susceptibility, canonical_kernel, rigid_vector, &
                                                absolute_residual, relative_residual, rigid_overlap, residual_vector)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: static_susceptibility(:, :), canonical_kernel(:, :), rigid_vector(:)
      real(rp), intent(out) :: absolute_residual, relative_residual
      complex(rp), intent(out), optional :: rigid_overlap
      complex(rp), intent(out), optional :: residual_vector(:)
      complex(rp), allocatable :: field(:), response(:), residual(:)
      real(rp) :: target_norm

      allocate(field(space%ndim), response(space%ndim), residual(space%ndim))
      call response_apply_operator(space, canonical_kernel, rigid_vector, field)
      call response_apply_operator(space, static_susceptibility, field, response)
      residual = response - rigid_vector
      absolute_residual = response_vector_norm(space, residual)
      target_norm = response_vector_norm(space, rigid_vector)
      if (target_norm > tiny(1.0_rp)) then
         relative_residual = absolute_residual/target_norm
      else
         relative_residual = absolute_residual
      end if
      if (present(rigid_overlap)) rigid_overlap = response_rigid_vector_overlap(space, residual, rigid_vector)
      if (present(residual_vector)) residual_vector = residual
   end subroutine evaluate_lr_alsda_static_residual

   subroutine require_initialized_space(space)
      type(response_space_layout), intent(in) :: space

      if (space%ndim < 1 .or. .not. allocated(space%radius) .or. .not. allocated(space%radial_weights) .or. &
          .not. allocated(space%metric_weights)) then
         error stop 'evaluate_lr_alsda_kernel: response space is not initialized'
      end if
   end subroutine require_initialized_space

   subroutine validate_ground_state(space, state, isite)
      type(response_space_layout), intent(in) :: space
      type(radial_ground_state), intent(in) :: state
      integer, intent(in) :: isite
      real(rp), parameter :: mesh_tolerance = 2.0e-13_rp

      if (.not. state%valid .or. .not. state%accepted) then
         error stop 'evaluate_lr_alsda_kernel: every LR-01 radial state must be valid and accepted'
      end if
      if (.not. allocated(state%r) .or. .not. allocated(state%vxc_up) .or. .not. allocated(state%vxc_down)) then
         error stop 'evaluate_lr_alsda_kernel: accepted LR-01 state is missing exact radial XC arrays'
      end if
      if (size(state%r) /= space%npoint .or. size(state%vxc_up) /= space%npoint .or. &
          size(state%vxc_down) /= space%npoint) then
         error stop 'evaluate_lr_alsda_kernel: LR-01 radial state/response mesh size mismatch'
      end if
      if (abs(state%a - space%a) > mesh_tolerance*max(1.0_rp, abs(space%a)) .or. &
          abs(state%b - space%b) > mesh_tolerance*max(1.0_rp, abs(space%b)) .or. &
          maxval(abs(state%r - space%radius)) > mesh_tolerance*max(1.0_rp, maxval(abs(space%radius)))) then
         error stop 'evaluate_lr_alsda_kernel: LR-01 radial mesh provenance differs from response space'
      end if
      if (trim(state%xc_provenance%internal_energy_units) /= 'Ry') then
         error stop 'evaluate_lr_alsda_kernel: XC provenance is not in the canonical Ry energy unit'
      end if
      if (len_trim(state%xc_provenance%functional_name) == 0 .or. state%xc_provenance%txc == 0) then
         error stop 'evaluate_lr_alsda_kernel: accepted XC provenance is incomplete'
      end if
      if (.not. state%xc_provenance%spin_polarized) then
         error stop 'evaluate_lr_alsda_kernel: a spin-polarized accepted XC state is required'
      end if
      if (isite < 1) error stop 'evaluate_lr_alsda_kernel: invalid site'
   end subroutine validate_ground_state

   subroutine validate_requested_provenance(accepted, requested_functional, requested_backend, requested_txc)
      type(radial_xc_provenance), intent(in) :: accepted
      character(len=*), intent(in) :: requested_functional, requested_backend
      integer, intent(in) :: requested_txc

      if (.not. lr_alsda_provenance_matches(accepted, requested_functional, requested_backend, requested_txc)) then
         error stop 'evaluate_lr_alsda_kernel: requested XC functional/provenance does not match accepted LR-01 state'
      end if
   end subroutine validate_requested_provenance

   logical function same_provenance(left, right) result(same)
      type(radial_xc_provenance), intent(in) :: left, right

      same = trim(left%backend_name) == trim(right%backend_name) .and. trim(left%txch) == trim(right%txch) .and. &
         trim(left%functional_name) == trim(right%functional_name) .and. &
         trim(left%mapping_quality) == trim(right%mapping_quality) .and. &
         trim(left%internal_energy_units) == trim(right%internal_energy_units) .and. left%txc == right%txc .and. &
         left%libxc_family == right%libxc_family .and. left%libxc_nspin == right%libxc_nspin .and. &
         left%n_components == right%n_components .and. left%use_libxc .eqv. right%use_libxc .and. &
         left%spin_polarized .eqv. right%spin_polarized .and. left%radial_gga .eqv. right%radial_gga
      if (.not. same) return
      if (left%n_components > 0) then
         same = allocated(left%component_ids) .and. allocated(right%component_ids) .and. &
            allocated(left%component_families) .and. allocated(right%component_families) .and. &
            allocated(left%component_kinds) .and. allocated(right%component_kinds)
         if (.not. same) return
         same = size(left%component_ids) == size(right%component_ids) .and. &
            size(left%component_families) == size(right%component_families) .and. &
            size(left%component_kinds) == size(right%component_kinds)
         if (.not. same) return
         same = all(left%component_ids == right%component_ids) .and. &
            all(left%component_families == right%component_families) .and. &
            all(left%component_kinds == right%component_kinds)
      end if
   end function same_provenance

   ! --- from lr_goldstone_sumrule_mod ---
   module subroutine evaluate_lr_goldstone_sumrule_request(request, result)
      type(lr_goldstone_sumrule_request), intent(in) :: request
      type(lr_goldstone_sumrule_result), intent(out) :: result

      call evaluate_lr_goldstone_sumrule_explicit(request%response_space, request%static_susceptibility, &
         request%magnetization, result, request%magnetization_label)
   end subroutine evaluate_lr_goldstone_sumrule_request

   !> Evaluate the mapped LCMM equation using the canonical LR-04 chiKS.
   module subroutine evaluate_lr_goldstone_sumrule_explicit(space, static_susceptibility, magnetization, result, &
                                                     magnetization_label)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: static_susceptibility(:, :)
      real(rp), intent(in) :: magnetization(:, :)
      type(lr_goldstone_sumrule_result), intent(out) :: result
      character(len=*), intent(in), optional :: magnetization_label

      type(lr_sumrule_linear_solve_result) :: solve
      complex(rp), allocatable :: target(:), field(:), response(:), residual(:)
      complex(rp), allocatable :: interaction_values(:, :), active_rhs(:)
      integer, allocatable :: active_indices(:)
      integer :: nunknown, isite, ir, i, j, flat, active_index
      real(rp) :: four_pi, target_metric_norm
      character(len=32) :: effective_label

      effective_label = lr_gsr_magnetization_pauli
      if (present(magnetization_label)) effective_label = trim(magnetization_label)
      call validate_inputs(space, static_susceptibility, magnetization, effective_label)

      four_pi = 4.0_rp*response_angular_pi
      nunknown = space%nsite*(space%npoint - 1)
      allocate(active_indices(nunknown), active_rhs(nunknown), result%gamma(nunknown, nunknown), &
               result%solution_active(nunknown), result%singular_values(nunknown))

      ! The active solve order is site-major, then positive-measure radial
      ! points.  It is deliberately derived through the LR-02R super-index.
      active_index = 0
      do isite = 1, space%nsite
         do ir = 2, space%npoint
            active_index = active_index + 1
            call l0_flat_index(space, isite, ir, active_indices(active_index))
            active_rhs(active_index) = cmplx(sqrt(four_pi)*magnetization(isite, ir), 0.0_rp, rp)
         end do
      end do

      ! LR-03 gives B_eff,00=4*pi*U_LCMM*m_00.  Since chiKS is already the
      ! LR-04 canonical B=chi_raw*W operator, no extra radial metric is
      ! inserted here: Gamma(:,j)=4*pi*chiKS(:,j)*m_00(j).
      do i = 1, nunknown
         do j = 1, nunknown
            result%gamma(i, j) = four_pi*static_susceptibility(active_indices(i), active_indices(j))* &
               active_rhs(j)
         end do
      end do

      call solve_lr_goldstone_equation(result%gamma, active_rhs, solve)
      result%solution_active = solve%solution
      result%singular_values = solve%singular_values
      result%rank = solve%rank
      result%unknowns = nunknown
      result%condition_number = solve%condition_number
      result%equation_residual_norm = solve%residual_norm
      result%equation_relative_residual = solve%relative_residual
      result%rank_deficient = solve%rank_deficient
      result%magnetization_label = effective_label

      allocate(result%u_lcmm(space%nsite, space%npoint), result%effective_field(space%nsite, space%npoint), &
               interaction_values(space%nsite, space%npoint), result%magnetization_response(space%ndim), &
               result%generated_field(space%ndim), result%generated_response(space%ndim), &
               result%residual_vector(space%ndim), target(space%ndim), field(space%ndim), response(space%ndim), &
               residual(space%ndim))
      result%u_lcmm = cmplx(0.0_rp, 0.0_rp, rp)
      result%effective_field = cmplx(0.0_rp, 0.0_rp, rp)
      result%magnetization_response = cmplx(0.0_rp, 0.0_rp, rp)
      result%generated_field = cmplx(0.0_rp, 0.0_rp, rp)
      result%generated_response = cmplx(0.0_rp, 0.0_rp, rp)
      result%residual_vector = cmplx(0.0_rp, 0.0_rp, rp)
      target = cmplx(0.0_rp, 0.0_rp, rp)
      field = cmplx(0.0_rp, 0.0_rp, rp)

      do isite = 1, space%nsite
         do ir = 2, space%npoint
            active_index = (isite - 1)*(space%npoint - 1) + ir - 1
            result%u_lcmm(isite, ir) = result%solution_active(active_index)
            result%effective_field(isite, ir) = four_pi*result%u_lcmm(isite, ir)*magnetization(isite, ir)
            call l0_flat_index(space, isite, ir, flat)
            target(flat) = active_rhs(active_index)
         end do
      end do

      interaction_values = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, space%nsite
         do ir = 2, space%npoint
            interaction_values(isite, ir) = four_pi*result%u_lcmm(isite, ir)
         end do
      end do
      allocate(result%canonical_interaction(space%ndim, space%ndim))
      call response_local_operator(space, interaction_values, result%canonical_interaction)
      result%magnetization_response = target
      ! Generate B_eff,00 through the canonical LR-04 action.  This keeps the
      ! production identity on one operator path and avoids a second local
      ! multiplication convention in this module.
      call response_apply_operator(space, result%canonical_interaction, target, field)
      result%generated_field = field
      call response_apply_operator(space, static_susceptibility, field, response)
      result%generated_response = response
      residual = response - target
      result%residual_vector = residual
      result%metric_residual_norm = response_vector_norm(space, residual)
      target_metric_norm = response_vector_norm(space, target)
      if (target_metric_norm > tiny(1.0_rp)) then
         result%metric_relative_residual = result%metric_residual_norm/target_metric_norm
      else
         result%metric_relative_residual = result%metric_residual_norm
      end if
      result%rigid_overlap = response_rigid_vector_overlap(space, residual, target)

      result%blocked = result%rank_deficient .or. result%metric_relative_residual > lr_gsr_residual_tolerance .or. &
         result%equation_relative_residual > lr_gsr_residual_tolerance
      if (result%rank_deficient) then
         result%status = 'BLOCKED: Gamma is numerically rank deficient'
      else if (result%metric_relative_residual > lr_gsr_residual_tolerance) then
         result%status = 'BLOCKED: full LR-04 Goldstone identity residual exceeds tolerance'
      else if (result%equation_relative_residual > lr_gsr_residual_tolerance) then
         result%status = 'BLOCKED: Gamma U=m residual exceeds tolerance'
      else
         result%blocked = .false.
         result%status = 'PASS: independent LCMM sum-rule identity'
      end if
   end subroutine evaluate_lr_goldstone_sumrule_explicit

   !> Solve a finite complex equation using an SVD-based least-squares kernel.
   !>
   !> RCOND=-1 asks LAPACK to use its machine-precision rank rule.  This is a
   !> numerical rank diagnostic, not a user-selected pseudoinverse cutoff.  A
   !> rank-deficient result is returned for inspection but marked BLOCKED by
   !> the production evaluator because LCMM does not select its null space.
   module subroutine solve_lr_goldstone_equation(gamma, rhs, result)
      complex(rp), intent(in) :: gamma(:, :), rhs(:)
      type(lr_sumrule_linear_solve_result), intent(out) :: result

      complex(rp), allocatable :: a_work(:, :), b_work(:, :), work(:), work_query(:)
      real(rp), allocatable :: rwork(:)
      integer :: m, n, nrhs, ldb, lwork, info, i
      real(rp), parameter :: rcond = -1.0_rp
      real(rp) :: largest_singular, smallest_positive
      external :: zgelss

      m = size(gamma, 1)
      n = size(gamma, 2)
      nrhs = 1
      if (m < 1 .or. n < 1 .or. size(rhs) /= m) then
         error stop 'solve_lr_goldstone_equation: inconsistent Gamma/rhs shapes'
      end if
      ldb = max(1, max(m, n))
      allocate(a_work(m, n), b_work(ldb, nrhs), result%solution(n), &
               result%singular_values(min(m, n)), rwork(max(1, 5*min(m, n))))
      a_work = gamma
      b_work = cmplx(0.0_rp, 0.0_rp, rp)
      b_work(1:m, 1) = rhs
      allocate(work_query(1))
      call zgelss(m, n, nrhs, a_work, m, b_work, ldb, result%singular_values, rcond, result%rank, &
         work_query, -1, rwork, info)
      if (info /= 0) error stop 'solve_lr_goldstone_equation: LAPACK ZGELSS workspace query failed'
      lwork = max(1, int(real(work_query(1), rp)))
      allocate(work(lwork))
      a_work = gamma
      b_work = cmplx(0.0_rp, 0.0_rp, rp)
      b_work(1:m, 1) = rhs
      call zgelss(m, n, nrhs, a_work, m, b_work, ldb, result%singular_values, rcond, result%rank, &
         work, lwork, rwork, info)
      if (info /= 0) error stop 'solve_lr_goldstone_equation: LAPACK ZGELSS failed'
      result%solution = b_work(1:n, 1)

      largest_singular = maxval(result%singular_values)
      smallest_positive = huge(1.0_rp)
      do i = 1, size(result%singular_values)
         if (result%singular_values(i) > 0.0_rp) then
            smallest_positive = min(smallest_positive, result%singular_values(i))
         end if
      end do
      if (result%rank < min(m, n)) then
         result%condition_number = huge(1.0_rp)
      else if (smallest_positive < huge(1.0_rp)) then
         result%condition_number = largest_singular/smallest_positive
      end if
      result%residual_norm = sqrt(sum(abs(matmul(gamma, result%solution) - rhs)**2))
      if (sqrt(sum(abs(rhs)**2)) > tiny(1.0_rp)) then
         result%relative_residual = result%residual_norm/sqrt(sum(abs(rhs)**2))
      else
         result%relative_residual = result%residual_norm
      end if
      result%rank_deficient = result%rank < min(m, n)
      result%blocked = result%rank_deficient .or. result%relative_residual > lr_gsr_residual_tolerance
      if (result%rank_deficient) then
         result%status = 'BLOCKED: SVD numerical rank is below equation dimension'
      else if (result%relative_residual > lr_gsr_residual_tolerance) then
         result%status = 'BLOCKED: SVD equation residual exceeds tolerance'
      else
         result%blocked = .false.
         result%status = 'PASS: full-rank SVD solve'
      end if
   end subroutine solve_lr_goldstone_equation

   subroutine validate_inputs(space, static_susceptibility, magnetization, magnetization_label)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: static_susceptibility(:, :)
      real(rp), intent(in) :: magnetization(:, :)
      character(len=*), intent(in) :: magnetization_label
      integer :: isite, ir

      if (space%ndim < 1 .or. .not. allocated(space%radius) .or. .not. allocated(space%radial_weights) .or. &
          .not. allocated(space%metric_weights)) then
         error stop 'evaluate_lr_goldstone_sumrule: response space is not initialized'
      end if
      if (space%nchannel /= 1) then
         error stop 'evaluate_lr_goldstone_sumrule: one certified circular response channel is required'
      end if
      if (any(shape(static_susceptibility) /= [space%ndim, space%ndim])) then
         error stop 'evaluate_lr_goldstone_sumrule: static chiKS shape mismatch'
      end if
      if (any(shape(magnetization) /= [space%nsite, space%npoint])) then
         error stop 'evaluate_lr_goldstone_sumrule: magnetization shape mismatch'
      end if
      if (trim(magnetization_label) /= lr_gsr_magnetization_pauli) then
         error stop 'evaluate_lr_goldstone_sumrule: only the explicitly labelled Pauli magnetization is supported'
      end if
      do isite = 1, space%nsite
         do ir = 2, space%npoint
            if (magnetization(isite, ir) == 0.0_rp) then
               error stop 'evaluate_lr_goldstone_sumrule: zero Pauli magnetization on active grid; capability BLOCKED'
            end if
         end do
      end do
   end subroutine validate_inputs

   subroutine l0_flat_index(space, isite, ir, flat)
      type(response_space_layout), intent(in) :: space
      integer, intent(in) :: isite, ir
      integer, intent(out) :: flat
      type(response_super_index) :: item

      item = response_super_index(isite, 0, 0, ir, 1)
      call response_flatten_superindex(item, space%nsite, space%response_lmax, space%npoint, space%nchannel, flat)
   end subroutine l0_flat_index

   ! --- from lr_projected_interacting_response_mod ---
   module subroutine evaluate_projected_mills_interaction(request, result)
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
   module subroutine evaluate_projected_mills_from_reciprocal(contract, radial_bases, reciprocal_obj, energy, moment, result)
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

   module subroutine evaluate_projected_dyson(request, result)
      type(projected_dyson_request), intent(in) :: request
      type(projected_dyson_result), intent(out) :: result

      complex(rp), allocatable :: residual(:, :), chi(:, :), denominator(:, :), loss(:, :)
      integer :: nsite, nw, iw, i, info
      real(rp) :: dnorm, rnorm
      complex(rp) :: trace_value

      call validate_projected_dyson_request(request, nsite, nw)
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

   subroutine validate_projected_dyson_request(request, nsite, nw)
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
   end subroutine validate_projected_dyson_request

   ! --- from lr_projected_juelich_interaction_mod ---
   !> Select the finest positive static eta and the distinct second-finest
   !> positive eta independently of the order supplied by the user.  The
   !> second value is the holdout/stability partner; an ambiguous duplicate is
   !> rejected instead of being treated as a stability comparison.
   module subroutine select_projected_juelich_eta_indices(eta_values, selected_index, holdout_index)
      real(rp), intent(in) :: eta_values(:)
      integer, intent(out) :: selected_index, holdout_index
      integer :: i, j
      real(rp) :: second_eta, scale

      if (size(eta_values) < 1) error stop 'DRESP-05 eta selection: empty eta ladder'
      if (any(.not. ieee_is_finite(eta_values)) .or. any(eta_values <= 0.0_rp)) then
         error stop 'DRESP-05 eta selection: every eta must be finite and positive'
      end if
      do i = 1, size(eta_values) - 1
         do j = i + 1, size(eta_values)
            scale = max(1.0_rp, abs(eta_values(i)), abs(eta_values(j)))
            if (abs(eta_values(i) - eta_values(j)) <= 100.0_rp*epsilon(1.0_rp)*scale) then
               error stop 'DRESP-05 eta selection: duplicated eta makes stability comparison ambiguous'
            end if
         end do
      end do
      selected_index = 1
      do i = 2, size(eta_values)
         if (eta_values(i) < eta_values(selected_index)) selected_index = i
      end do
      holdout_index = 0
      second_eta = huge(1.0_rp)
      do i = 1, size(eta_values)
         if (i == selected_index) cycle
         if (eta_values(i) < second_eta) then
            second_eta = eta_values(i)
            holdout_index = i
         end if
      end do
   end subroutine select_projected_juelich_eta_indices

   module subroutine evaluate_projected_juelich_interaction(request, result)
      type(projected_juelich_request), intent(in) :: request
      type(projected_juelich_result), intent(out) :: result

      complex(rp), allocatable :: gamma(:, :), rhs(:, :), a_work(:, :), b_work(:, :), work(:)
      real(rp), allocatable :: singular_values(:), rwork(:), real_a(:, :), real_b(:, :), real_singular_values(:), &
         real_work(:), real_work_query(:)
      complex(rp) :: complex_work_query(1)
      real(rp) :: complex_scale
      integer :: nsite, lwork, real_lwork, info, i, j
      integer :: rank, real_rank, mreal, ldb
      real(rp), parameter :: rcond = -1.0_rp

      call validate_juelich_request(request, nsite)
      result%selector = trim(request%selector)
      result%static_eta = request%static_eta
      result%q = request%q
      result%channel = trim(request%channel)
      result%provenance = trim(request%state_provenance)//' | chi0='//trim(request%chi0_provenance)
      result%nsite = nsite
      allocate(result%projected_moment(nsite), result%gamma(nsite, nsite), result%interaction_matrix(nsite, nsite), &
         result%interaction_U_complex(nsite), result%interaction_U_real(nsite), result%interaction_U(nsite), &
         result%singular_values(nsite), result%real_singular_values(nsite))
      result%projected_moment = request%projected_moment
      result%gamma = cmplx(0.0_rp, 0.0_rp, rp)
      do j = 1, nsite
         result%gamma(:, j) = request%static_chi0(:, j)*request%projected_moment(j)
      end do
      result%interaction_matrix = cmplx(0.0_rp, 0.0_rp, rp)
      result%interaction_U_complex = cmplx(0.0_rp, 0.0_rp, rp)
      result%interaction_U_real = 0.0_rp
      result%interaction_U = 0.0_rp

      allocate(a_work(nsite, nsite), b_work(nsite, 1), rwork(max(1, 5*nsite)))
      a_work = result%gamma
      b_work(:, 1) = cmplx(request%projected_moment, 0.0_rp, rp)
      call zgelss(nsite, nsite, 1, a_work, nsite, b_work, nsite, result%singular_values, rcond, rank, &
         complex_work_query, -1, rwork, info)
      if (info /= 0) error stop 'DRESP-05 Juelich: complex ZGELSS workspace query failed'
      lwork = max(1, nint(real(complex_work_query(1), rp)))
      allocate(work(lwork))
      a_work = result%gamma
      b_work(:, 1) = cmplx(request%projected_moment, 0.0_rp, rp)
      call zgelss(nsite, nsite, 1, a_work, nsite, b_work, nsite, result%singular_values, rcond, rank, &
         work, lwork, rwork, info)
      if (info /= 0) error stop 'DRESP-05 Juelich: complex ZGELSS failed'
      result%rank = rank
      result%interaction_U_complex = b_work(:, 1)
      complex_scale = sqrt(sum(abs(request%projected_moment)**2))
      result%equation_residual = sqrt(sum(abs(matmul(result%gamma, result%interaction_U_complex) - &
         cmplx(request%projected_moment, 0.0_rp, rp))**2))
      result%relative_residual = result%equation_residual/max(complex_scale, tiny(1.0_rp))
      result%condition_number = singular_condition(result%singular_values, result%rank)
      result%imaginary_U_ratio = maxval(abs(aimag(result%interaction_U_complex)))/ &
         max(maxval(abs(real(result%interaction_U_complex, rp))), epsilon(1.0_rp))

      ! Route B: the unknown is real.  Use DGELSS on the stacked real and
      ! imaginary equations, so taking REAL(U_complex) is never substituted for
      ! the physical constrained solve.
      mreal = 2*nsite
      ldb = mreal
      allocate(real_a(mreal, nsite), real_b(ldb, 1), real_singular_values(nsite), real_work_query(1))
      real_a(1:nsite, :) = real(result%gamma, rp)
      real_a(nsite+1:mreal, :) = aimag(result%gamma)
      real_b = 0.0_rp
      real_b(1:nsite, 1) = request%projected_moment
      call dgelss(mreal, nsite, 1, real_a, mreal, real_b, ldb, real_singular_values, rcond, real_rank, &
         real_work_query, -1, info)
      if (info /= 0) error stop 'DRESP-05 Juelich: real DGELSS workspace query failed'
      real_lwork = max(1, nint(real_work_query(1)))
      allocate(real_work(real_lwork))
      real_a(1:nsite, :) = real(result%gamma, rp)
      real_a(nsite+1:mreal, :) = aimag(result%gamma)
      real_b = 0.0_rp
      real_b(1:nsite, 1) = request%projected_moment
      call dgelss(mreal, nsite, 1, real_a, mreal, real_b, ldb, real_singular_values, rcond, real_rank, &
         real_work, real_lwork, info)
      if (info /= 0) error stop 'DRESP-05 Juelich: real DGELSS failed'
      result%real_rank = real_rank
      result%interaction_U_real = real_b(1:nsite, 1)
      result%real_constrained_residual = sqrt(sum((matmul(real(result%gamma, rp), result%interaction_U_real) - &
         request%projected_moment)**2) + sum((matmul(aimag(result%gamma), result%interaction_U_real))**2))
      result%relative_real_constrained_residual = result%real_constrained_residual/max(real_scale_from(request%projected_moment), &
         tiny(1.0_rp))
      result%real_condition_number = singular_condition(real_singular_values, result%real_rank)
      do i = 1, nsite
         result%interaction_matrix(i, i) = cmplx(result%interaction_U_real(i), 0.0_rp, rp)
      end do
      result%interaction_U = result%interaction_U_real

      if (rank < nsite .or. real_rank < nsite) then
         result%classification = projected_juelich_rank_deficient
      else if (result%condition_number > projected_juelich_condition_limit .or. &
               result%real_condition_number > projected_juelich_condition_limit) then
         result%classification = projected_juelich_unsupported
      else if (result%relative_real_constrained_residual > projected_juelich_residual_tolerance) then
         result%classification = projected_juelich_projected_local
      else
         result%classification = projected_juelich_exact_local
      end if

      deallocate(a_work, b_work, work, rwork, real_a, real_b, real_singular_values, real_work_query, real_work)
   end subroutine evaluate_projected_juelich_interaction

   module subroutine evaluate_projected_juelich_holdout(static_chi0, projected_moment, interaction_U, residual, relative_residual)
      complex(rp), intent(in) :: static_chi0(:, :)
      real(rp), intent(in) :: projected_moment(:), interaction_U(:)
      real(rp), intent(out) :: residual, relative_residual
      complex(rp), allocatable :: vector(:)
      integer :: nsite

      nsite = size(projected_moment)
      allocate(vector(nsite))
      ! Holdout means that U is frozen while chi0 is replaced.  The residual
      ! is therefore the Ward action (I-chi0*U)M, not the construction solve.
      vector = cmplx(projected_moment, 0.0_rp, rp) - matmul(static_chi0, &
         cmplx(interaction_U*projected_moment, 0.0_rp, rp))
      residual = sqrt(sum(abs(vector)**2))
      relative_residual = residual/max(sqrt(sum(projected_moment**2)), tiny(1.0_rp))
      deallocate(vector)
   end subroutine evaluate_projected_juelich_holdout

   module subroutine assess_projected_juelich_eta_stability(result, previous_result, is_stable)
      type(projected_juelich_result), intent(inout) :: result
      type(projected_juelich_result), intent(in) :: previous_result
      logical, intent(out) :: is_stable
      real(rp) :: scale, difference

      difference = maxval(abs(result%interaction_U_real - previous_result%interaction_U_real))
      scale = max(maxval(abs(result%interaction_U_real)), maxval(abs(previous_result%interaction_U_real)), tiny(1.0_rp))
      is_stable = difference/scale <= projected_juelich_eta_relative_tolerance
      result%eta_stable = is_stable
      write(result%eta_stability, '(a,es12.4,a,es12.4)') 'relative_delta_U=', difference/scale, &
         '; threshold=', projected_juelich_eta_relative_tolerance
      if (.not. is_stable .and. trim(result%classification) /= projected_juelich_rank_deficient .and. &
          trim(result%classification) /= projected_juelich_unsupported) then
         result%classification = projected_juelich_eta_limited
      end if
   end subroutine assess_projected_juelich_eta_stability

   subroutine validate_juelich_request(request, nsite)
      type(projected_juelich_request), intent(in) :: request
      integer, intent(out) :: nsite

      nsite = size(request%projected_moment)
      if (nsite < 1 .or. any(shape(request%static_chi0) /= [nsite, nsite]) .or. request%static_eta <= 0.0_rp .or. &
          len_trim(request%selector) == 0 .or. len_trim(request%channel) == 0) then
         error stop 'DRESP-05 Juelich: invalid projected request'
      end if
      if (any(.not. ieee_is_finite(request%projected_moment)) .or. &
          any(.not. ieee_is_finite(real(request%static_chi0, rp))) .or. &
          any(.not. ieee_is_finite(aimag(request%static_chi0)))) then
         error stop 'DRESP-05 Juelich: non-finite projected input'
      end if
   end subroutine validate_juelich_request

   real(rp) function singular_condition(singular_values, rank) result(value)
      real(rp), intent(in) :: singular_values(:)
      integer, intent(in) :: rank
      real(rp) :: smallest

      if (rank < 1) then
         value = huge(1.0_rp)
         return
      end if
      smallest = minval(singular_values(1:rank))
      if (smallest <= 0.0_rp) then
         value = huge(1.0_rp)
      else
         value = maxval(singular_values)/smallest
      end if
   end function singular_condition

   real(rp) function real_scale_from(vector) result(value)
      real(rp), intent(in) :: vector(:)
      value = sqrt(sum(vector**2))
   end function real_scale_from

   ! --- from lr_compact_static_interaction_mod ---
   !> Project a point-space vector with the actual weighted modes stored by the
   !> product basis.  The origin is an exact zero-measure row and contributes
   !> zero without division or regularization.
   module subroutine compact_project_point_vector(space, product, point_vector, compact_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: point_vector(:)
      complex(rp), intent(out) :: compact_vector(:)
      integer :: site, response_l, response_m, mode, ir, point_flat, compact_flat
      real(rp) :: sqrt_weight
      type(response_super_index) :: item

      call validate_point_to_compact_shapes(space, product, point_vector, compact_vector)
      compact_vector = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, product%nsite
         do response_l = 0, product%response_lmax
            do response_m = -response_l, response_l
               do mode = 1, product%blocks(site, response_l)%rank
                  compact_flat = product%flat_index(site, response_l, response_m, mode)
                  do ir = 2, space%npoint
                     sqrt_weight = sqrt(space%radial_weights(ir))
                     item = response_super_index(site, response_l, response_m, ir, 1)
                     call response_flatten_superindex(item, space%nsite, space%response_lmax, &
                        space%npoint, space%nchannel, point_flat)
                     compact_vector(compact_flat) = compact_vector(compact_flat) + &
                        conjg(product%blocks(site, response_l)%weighted_modes(ir, mode))*sqrt_weight*point_vector(point_flat)
                  end do
               end do
            end do
         end do
      end do
   end subroutine compact_project_point_vector

   !> Reconstruct the retained product-space component of a compact vector.
   !> Null-measure origin values are set to the finite zero extension.
   module subroutine compact_reconstruct_point_vector(space, product, compact_vector, point_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: compact_vector(:)
      complex(rp), intent(out) :: point_vector(:)
      integer :: site, response_l, response_m, mode, ir, point_flat, compact_flat
      real(rp) :: sqrt_weight
      type(response_super_index) :: item

      call validate_compact_to_point_shapes(space, product, compact_vector, point_vector)
      point_vector = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, product%nsite
         do response_l = 0, product%response_lmax
            do response_m = -response_l, response_l
               do mode = 1, product%blocks(site, response_l)%rank
                  compact_flat = product%flat_index(site, response_l, response_m, mode)
                  do ir = 2, space%npoint
                     sqrt_weight = sqrt(space%radial_weights(ir))
                     if (sqrt_weight <= 0.0_rp) cycle
                     item = response_super_index(site, response_l, response_m, ir, 1)
                     call response_flatten_superindex(item, space%nsite, space%response_lmax, &
                        space%npoint, space%nchannel, point_flat)
                     point_vector(point_flat) = point_vector(point_flat) + &
                        product%blocks(site, response_l)%weighted_modes(ir, mode)/sqrt_weight*compact_vector(compact_flat)
                  end do
               end do
            end do
         end do
      end do
   end subroutine compact_reconstruct_point_vector

   !> Project a spherical local pointwise operator.  Since the point basis is
   !> V=W**(-1/2)U, the W factors from the metric projection and reconstruction
   !> cancel exactly, leaving U^H K U.
   module subroutine compact_project_local_operator_complex(space, product, values, compact_operator)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: values(:, :)
      complex(rp), intent(out) :: compact_operator(:, :)
      integer :: site, response_l, response_m, mode_out, mode_in, ir
      integer :: flat_out, flat_in

      call validate_operator_shapes(space, product, compact_operator)
      compact_operator = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, product%nsite
         do response_l = 0, product%response_lmax
            do response_m = -response_l, response_l
               do mode_out = 1, product%blocks(site, response_l)%rank
                  flat_out = product%flat_index(site, response_l, response_m, mode_out)
                  do mode_in = 1, product%blocks(site, response_l)%rank
                     flat_in = product%flat_index(site, response_l, response_m, mode_in)
                     do ir = 2, space%npoint
                        compact_operator(flat_out, flat_in) = compact_operator(flat_out, flat_in) + &
                           conjg(product%blocks(site, response_l)%weighted_modes(ir, mode_out))*values(site, ir)* &
                           product%blocks(site, response_l)%weighted_modes(ir, mode_in)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine compact_project_local_operator_complex

   module subroutine compact_apply_local_operator(space, product, values, vector, result)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: values(:, :), vector(:)
      complex(rp), intent(out) :: result(:)
      complex(rp), allocatable :: operator(:, :)

      call validate_common_layout(space, product)
      allocate(operator(product%product_dimension, product%product_dimension))
      call compact_project_local_operator(space, product, values, operator)
      result = matmul(operator, vector)
   end subroutine compact_apply_local_operator

   !> Build the normalized L=0 magnetization vector used by both static
   !> diagnostics.  The point-space coefficient is m_00=sqrt(4*pi)*m_z.
   module subroutine compact_project_magnetization(space, product, magnetization, compact_vector, point_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      real(rp), intent(in) :: magnetization(:, :)
      complex(rp), intent(out) :: compact_vector(:)
      complex(rp), allocatable, intent(out), optional :: point_vector(:)
      complex(rp), allocatable :: point(:)
      integer :: site, ir, flat
      type(response_super_index) :: item

      allocate(point(space%ndim))
      point = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, space%nsite
         do ir = 1, space%npoint
            item = response_super_index(site, 0, 0, ir, 1)
            call response_flatten_superindex(item, space%nsite, space%response_lmax, space%npoint, &
               space%nchannel, flat)
            point(flat) = cmplx(sqrt(4.0_rp*response_angular_pi)*magnetization(site, ir), 0.0_rp, rp)
         end do
      end do
      call compact_project_point_vector(space, product, point, compact_vector)
      if (present(point_vector)) then
         allocate(point_vector(space%ndim))
         point_vector = point
      end if
   end subroutine compact_project_magnetization

   !> Report the weighted point-space norm of a vector and its defect from a
   !> compact reconstruction.  This is the same diagonal radial metric used
   !> by compact_project_point_vector and compact_reconstruct_point_vector.
   module subroutine compact_weighted_projection_diagnostics(space, point_vector, projected_point_vector, weighted_norm, &
                                                      residual_norm, relative_residual)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: point_vector(:), projected_point_vector(:)
      real(rp), intent(out) :: weighted_norm, residual_norm, relative_residual
      integer :: site, response_l, response_m, ir, channel, flat
      real(rp) :: weight
      type(response_super_index) :: item

      weighted_norm = 0.0_rp
      residual_norm = 0.0_rp
      do site = 1, space%nsite
         do response_l = 0, space%response_lmax
            do response_m = -response_l, response_l
               do ir = 1, space%npoint
                  weight = space%radial_weights(ir)
                  if (weight <= 0.0_rp) cycle
                  do channel = 1, space%nchannel
                     item = response_super_index(site, response_l, response_m, ir, channel)
                     call response_flatten_superindex(item, space%nsite, space%response_lmax, &
                        space%npoint, space%nchannel, flat)
                     weighted_norm = weighted_norm + weight*abs(point_vector(flat))**2
                     residual_norm = residual_norm + weight*abs(point_vector(flat) - projected_point_vector(flat))**2
                  end do
               end do
            end do
         end do
      end do
      weighted_norm = sqrt(weighted_norm)
      residual_norm = sqrt(residual_norm)
      relative_residual = residual_norm/max(weighted_norm, tiny(1.0_rp))
   end subroutine compact_weighted_projection_diagnostics

   !> Solve the independent compact local LCMM/GSR equation.  The published
   !> U_LCMM radial function is represented in the retained L=0 product span;
   !> each unknown basis function is reconstructed to point space and mapped
   !> through the normalized m_00 source.  No ALSDA quantity enters Gamma.
   module subroutine evaluate_compact_goldstone_sumrule(space, product, static_susceptibility, magnetization, result)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: static_susceptibility(:, :)
      real(rp), intent(in) :: magnetization(:, :)
      type(lr_compact_gsr_result), intent(out) :: result
      complex(rp), allocatable :: target(:), target_point(:), target_point_compact(:), field_basis(:, :), gamma(:, :)
      complex(rp), allocatable :: point_u(:), point_field(:), values(:, :), field(:), response(:), residual(:)
      complex(rp), allocatable :: assembled_field(:), solve_residual(:), residual_difference(:)
      complex(rp), allocatable :: a_work(:, :), b_work(:, :), work(:), work_query(:)
      real(rp), allocatable :: rwork(:)
      integer, allocatable :: unknown_site(:), unknown_mode(:)
      integer :: nrow, nunknown, l0_rank, site, mode, ir, unknown, info, lwork, i
      integer :: ldb, min_dimension, point_flat
      real(rp), parameter :: rcond = -1.0_rp
      real(rp) :: four_pi, largest_singular, smallest_positive, target_norm
      real(rp) :: sqrt_weight
      type(response_super_index) :: item

      call validate_gsr_inputs(space, product, static_susceptibility, magnetization)
      nrow = product%product_dimension
      l0_rank = 0
      do site = 1, product%nsite
         l0_rank = l0_rank + product%blocks(site, 0)%rank
      end do
      nunknown = l0_rank
      if (nunknown < 1) error stop 'evaluate_compact_goldstone_sumrule: no retained L=0 interaction modes'
      four_pi = 4.0_rp*response_angular_pi
      allocate(target(nrow), target_point(space%ndim), target_point_compact(space%ndim), field_basis(nrow, nunknown), gamma(nrow, nunknown), &
         unknown_site(nunknown), unknown_mode(nunknown), point_u(space%ndim), point_field(space%ndim), &
         values(space%nsite, space%npoint), field(nrow), response(nrow), residual(nrow), assembled_field(nrow), &
         solve_residual(nrow), residual_difference(nrow))
      call compact_project_magnetization(space, product, magnetization, target, target_point)
      call compact_reconstruct_point_vector(space, product, target, target_point_compact)
      call compact_weighted_projection_diagnostics(space, target_point, target_point_compact, result%magnetization_weighted_norm, &
         result%magnetization_projection_residual_norm, result%magnetization_projection_relative_residual)
      field_basis = cmplx(0.0_rp, 0.0_rp, rp)
      unknown = 0
      do site = 1, product%nsite
         do mode = 1, product%blocks(site, 0)%rank
            unknown = unknown + 1
            unknown_site(unknown) = site
            unknown_mode(unknown) = mode
            point_field = cmplx(0.0_rp, 0.0_rp, rp)
            do ir = 2, space%npoint
               sqrt_weight = sqrt(space%radial_weights(ir))
               if (sqrt_weight <= 0.0_rp) cycle
               item = response_super_index(site, 0, 0, ir, 1)
               call response_flatten_superindex(item, space%nsite, space%response_lmax, &
                  space%npoint, space%nchannel, point_flat)
               ! U_LCMM = weighted-mode/sqrt(W) times the compact unknown.
               ! The GSR column must act on R*c, not on the unprojected
               ! physical magnetization.  This makes it the same action as
               ! Kc*c used after the solve.
               point_field(point_flat) = four_pi*product%blocks(site, 0)%weighted_modes(ir, mode)/sqrt_weight * &
                  target_point_compact(point_flat)
            end do
            call compact_project_point_vector(space, product, point_field, field_basis(:, unknown))
         end do
      end do
      gamma = matmul(static_susceptibility, field_basis)

      ldb = max(1, max(nrow, nunknown))
      min_dimension = min(nrow, nunknown)
      allocate(a_work(nrow, nunknown), b_work(ldb, 1), result%solution(nunknown), &
         result%singular_values(min_dimension), rwork(max(1, 5*min_dimension)), work_query(1))
      a_work = gamma
      b_work = cmplx(0.0_rp, 0.0_rp, rp)
      b_work(1:nrow, 1) = target
      call zgelss(nrow, nunknown, 1, a_work, nrow, b_work, ldb, result%singular_values, rcond, result%rank, &
         work_query, -1, rwork, info)
      if (info /= 0) error stop 'evaluate_compact_goldstone_sumrule: ZGELSS workspace query failed'
      lwork = max(1, int(real(work_query(1), rp)))
      allocate(work(lwork))
      a_work = gamma
      b_work = cmplx(0.0_rp, 0.0_rp, rp)
      b_work(1:nrow, 1) = target
      call zgelss(nrow, nunknown, 1, a_work, nrow, b_work, ldb, result%singular_values, rcond, result%rank, &
         work, lwork, rwork, info)
      if (info /= 0) error stop 'evaluate_compact_goldstone_sumrule: ZGELSS failed'
      result%solution = b_work(1:nunknown, 1)
      result%coefficient_norm = sqrt(sum(abs(result%solution)**2))

      allocate(result%gamma(nrow, nunknown), result%magnetization_response(nrow), result%generated_field(nrow), &
         result%generated_response(nrow), result%residual_vector(nrow), result%u_lcmm(space%nsite, space%npoint), &
         result%effective_field(space%nsite, space%npoint), result%canonical_interaction(nrow, nrow))
      result%gamma = gamma
      result%magnetization_response = target
      result%u_lcmm = cmplx(0.0_rp, 0.0_rp, rp)
      result%effective_field = cmplx(0.0_rp, 0.0_rp, rp)
      point_u = cmplx(0.0_rp, 0.0_rp, rp)
      do unknown = 1, nunknown
         site = unknown_site(unknown)
         mode = unknown_mode(unknown)
         do ir = 2, space%npoint
            sqrt_weight = sqrt(space%radial_weights(ir))
            if (sqrt_weight <= 0.0_rp) cycle
            item = response_super_index(site, 0, 0, ir, 1)
            call response_flatten_superindex(item, space%nsite, space%response_lmax, &
               space%npoint, space%nchannel, point_flat)
            point_u(point_flat) = point_u(point_flat) + &
               product%blocks(site, 0)%weighted_modes(ir, mode)/sqrt_weight*result%solution(unknown)
         end do
      end do
      values = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, space%nsite
         do ir = 2, space%npoint
            item = response_super_index(site, 0, 0, ir, 1)
            call response_flatten_superindex(item, space%nsite, space%response_lmax, &
               space%npoint, space%nchannel, point_flat)
            result%u_lcmm(site, ir) = point_u(point_flat)
            result%effective_field(site, ir) = four_pi*result%u_lcmm(site, ir)*magnetization(site, ir)
            values(site, ir) = four_pi*result%u_lcmm(site, ir)
         end do
      end do
      call compact_project_local_operator(space, product, values, result%canonical_interaction)
      assembled_field = cmplx(0.0_rp, 0.0_rp, rp)
      do unknown = 1, nunknown
         assembled_field = assembled_field + result%solution(unknown)*field_basis(:, unknown)
      end do
      result%generated_field = matmul(result%canonical_interaction, target)
      result%assembled_action_sum_norm = sqrt(sum(abs(assembled_field)**2))
      result%assembled_action_matrix_norm = sqrt(sum(abs(result%generated_field)**2))
      result%assembled_action_difference_norm = sqrt(sum(abs(assembled_field - result%generated_field)**2))
      result%assembled_action_difference_relative = result%assembled_action_difference_norm / &
         max(result%assembled_action_matrix_norm, tiny(1.0_rp))
      result%assembled_action_difference_max_component = maxval(abs(assembled_field - result%generated_field))
      result%generated_response = matmul(static_susceptibility, result%generated_field)
      residual = result%generated_response - target
      result%residual_vector = residual
      solve_residual = matmul(gamma, result%solution) - target
      residual_difference = solve_residual - residual
      result%equation_residual_norm = sqrt(sum(abs(solve_residual)**2))
      target_norm = sqrt(sum(abs(target)**2))
      result%equation_relative_residual = result%equation_residual_norm/max(target_norm, tiny(1.0_rp))
      result%residual_norm = sqrt(sum(abs(residual)**2))
      result%relative_residual = result%residual_norm/max(target_norm, tiny(1.0_rp))
      result%residual_difference_norm = sqrt(sum(abs(residual_difference)**2))
      ! Scale the residual difference by the larger of the equation right-hand
      ! side and the assembled action.  When both residuals are at roundoff,
      ! scaling by either residual itself would manufacture a meaningless O(1)
      ! relative difference; for an ill-conditioned solve, the assembled-action
      ! scale also exposes roundoff amplified by the large coefficients.
      result%residual_difference_relative = result%residual_difference_norm / &
         max(max(target_norm, result%assembled_action_matrix_norm), tiny(1.0_rp))
      result%residual_difference_max_component = maxval(abs(residual_difference))
      result%max_residual_component = maxval(abs(residual))
      result%rigid_overlap = sum(conjg(residual)*target)
      largest_singular = maxval(result%singular_values)
      smallest_positive = huge(1.0_rp)
      do i = 1, size(result%singular_values)
         if (result%singular_values(i) > 0.0_rp) smallest_positive = min(smallest_positive, result%singular_values(i))
      end do
      if (result%rank == nunknown .and. smallest_positive < huge(1.0_rp)) then
         result%condition_number = largest_singular/smallest_positive
      end if
      result%svd_rcond = rcond
      result%svd_cutoff = epsilon(1.0_rp)*largest_singular
      result%equation_rows = nrow
      result%unknowns = nunknown
      result%rank_deficient = result%rank < nunknown
      result%blocked = result%rank_deficient .or. result%assembled_action_difference_relative > lr_compact_gsr_action_tolerance .or. &
         result%equation_relative_residual > lr_compact_gsr_residual_tolerance .or. &
         result%relative_residual > lr_compact_gsr_residual_tolerance
      if (result%rank_deficient) then
         result%status = 'BLOCKED: GSR rank/conditioning'
      else if (result%assembled_action_difference_relative > lr_compact_gsr_action_tolerance) then
         result%status = 'BLOCKED: compact GSR operator action'
      else if (result%equation_relative_residual > lr_compact_gsr_residual_tolerance) then
         result%status = 'BLOCKED: GSR rank/conditioning'
      else if (result%relative_residual > lr_compact_gsr_residual_tolerance) then
         result%status = 'BLOCKED: GSR rank/conditioning'
      else
         result%blocked = .false.
         result%status = 'PASS: independent compact LCMM sum-rule identity'
      end if
   end subroutine evaluate_compact_goldstone_sumrule

   subroutine validate_common_layout(space, product)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product

      if (space%nchannel /= 1 .or. space%nsite /= product%nsite .or. space%response_lmax /= product%response_lmax) then
         error stop 'compact projection: response/product layout mismatch'
      end if
      if (.not. allocated(space%radial_weights) .or. size(space%radial_weights) /= space%npoint) then
         error stop 'compact projection: LR-04 radial metric is unavailable'
      end if
      if (.not. allocated(product%blocks) .or. product%product_dimension < 1) then
         error stop 'compact projection: product representation is unavailable'
      end if
   end subroutine validate_common_layout

   subroutine validate_point_to_compact_shapes(space, product, point_vector, compact_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: point_vector(:)
      complex(rp), intent(out) :: compact_vector(:)

      call validate_common_layout(space, product)
      if (size(point_vector) /= space%ndim .or. size(compact_vector) /= product%product_dimension) then
         error stop 'compact projection: point-to-compact vector shape mismatch'
      end if
   end subroutine validate_point_to_compact_shapes

   subroutine validate_compact_to_point_shapes(space, product, compact_vector, point_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: compact_vector(:)
      complex(rp), intent(out) :: point_vector(:)

      call validate_common_layout(space, product)
      if (size(compact_vector) /= product%product_dimension .or. size(point_vector) /= space%ndim) then
         error stop 'compact projection: compact-to-point vector shape mismatch'
      end if
   end subroutine validate_compact_to_point_shapes

   subroutine validate_operator_shapes(space, product, operator)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: operator(:, :)

      if (space%nchannel /= 1 .or. space%nsite /= product%nsite .or. space%response_lmax /= product%response_lmax) then
         error stop 'compact operator projection: response/product layout mismatch'
      end if
      if (any(shape(operator) /= [product%product_dimension, product%product_dimension])) then
         error stop 'compact operator projection: compact operator shape mismatch'
      end if
      if (.not. allocated(space%radial_weights) .or. size(space%radial_weights) /= space%npoint) then
         error stop 'compact operator projection: LR-04 radial metric is unavailable'
      end if
   end subroutine validate_operator_shapes

   subroutine validate_gsr_inputs(space, product, susceptibility, magnetization)
      type(response_space_layout), intent(in) :: space
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: susceptibility(:, :)
      real(rp), intent(in) :: magnetization(:, :)

      call validate_common_layout(space, product)
      if (any(shape(susceptibility) /= [product%product_dimension, product%product_dimension])) then
         error stop 'evaluate_compact_goldstone_sumrule: compact susceptibility shape mismatch'
      end if
      if (any(shape(magnetization) /= [space%nsite, space%npoint])) then
         error stop 'evaluate_compact_goldstone_sumrule: magnetization shape mismatch'
      end if
      if (.not. all(ieee_is_finite(real(susceptibility, rp))) .or. .not. all(ieee_is_finite(aimag(susceptibility)))) then
         error stop 'evaluate_compact_goldstone_sumrule: compact susceptibility contains NaN or Inf'
      end if
   end subroutine validate_gsr_inputs

   ! --- from tddft_dyson_mod ---
   !> Evaluate the canonical LR-04 Dyson response for all requested frequencies.
   module subroutine evaluate_tddft_dyson(request, result)
      type(tddft_dyson_request), intent(in) :: request
      type(tddft_dyson_result), intent(out) :: result

      type(response_space_layout), pointer :: space
      complex(rp), allocatable :: identity(:, :), product(:, :), denominator(:, :), enhanced(:, :), loss(:, :)
      complex(rp), allocatable :: residual(:, :)
      integer :: n, nw, iw
      logical :: compact

      compact = request%compact_orthonormal
      if (compact) then
         call validate_compact_request(request, n)
      else
         space => request%response_space
         call validate_dyson_request(request, space)
         n = space%ndim
      end if

      nw = size(request%frequencies)
      allocate(result%frequencies(nw), result%ks_susceptibility(n, n, nw), &
         result%canonical_interaction(n, n), result%denominator(n, n, nw), &
         result%enhanced_susceptibility(n, n, nw), result%loss_matrix(n, n, nw), &
         result%denominator_min_singular_value(nw), result%denominator_max_singular_value(nw), &
         result%denominator_condition_number(nw), result%denominator_min_magnitude_eigenvalue(nw), &
         result%solve_info(nw), result%solve_succeeded(nw), &
         result%near_singular_collective_pole(nw), result%numerically_singular(nw), result%ill_conditioned(nw), &
         result%loss_metric_hermiticity_residual(nw), result%dyson_residual_frobenius(nw), &
         result%dyson_residual_relative(nw), result%dyson_residual_infinity(nw), result%frequency_status(nw))

      allocate(identity(n, n), product(n, n), denominator(n, n), enhanced(n, n), loss(n, n), residual(n, n))
      if (compact) then
         identity = cmplx(0.0_rp, 0.0_rp, rp)
         do iw = 1, n
            identity(iw, iw) = cmplx(1.0_rp, 0.0_rp, rp)
         end do
      else
         call response_identity_operator(space, identity)
      end if

      result%q = request%q
      result%frequencies = request%frequencies
      result%eta = request%eta
      result%channel = trim(request%channel)
      result%compact_orthonormal = compact
      if (compact) result%response_representation = lr_dyson_compact_response_representation
      result%interaction_route = trim(request%interaction_route)
      result%interaction_provenance = trim(request%interaction_provenance)
      result%electronic_state_provenance = trim(request%electronic_state_provenance)
      result%response_space_metadata = trim(request%response_space_metadata)
      if (len_trim(result%response_space_metadata) == 0 .and. compact) then
         write (result%response_space_metadata, '(a,i0,a)') 'dimension=', n, '; orthonormal compact product space'
      else if (len_trim(result%response_space_metadata) == 0) then
         write (result%response_space_metadata, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'ndim=', n, &
            ' active_dimension=', space%active_dimension, ' nsite=', space%nsite, &
            ' radial_points=', space%npoint, ' channels=', space%nchannel
      end if
      result%ks_susceptibility = request%ks_susceptibility
      result%canonical_interaction = request%canonical_interaction

      do iw = 1, nw
         ! Both operands are canonical LR-04 operators.  This ordinary
         ! multiplication is exactly the LR-04 composition rule; no W is
         ! inserted here.
         if (compact) then
            product = matmul(request%ks_susceptibility(:, :, iw), request%canonical_interaction)
         else
            call response_compose_operators(space, request%ks_susceptibility(:, :, iw), &
               request%canonical_interaction, product)
         end if
         denominator = identity - product
         result%denominator(:, :, iw) = denominator

         call solve_tddft_denominator(denominator, request%ks_susceptibility(:, :, iw), enhanced, result%solve_info(iw), &
            result%denominator_condition_number(iw), result%denominator_min_singular_value(iw), &
            result%denominator_max_singular_value(iw), result%numerically_singular(iw), &
            result%ill_conditioned(iw))
         call minimum_magnitude_eigenvalue(denominator, result%denominator_min_magnitude_eigenvalue(iw))
         result%enhanced_susceptibility(:, :, iw) = enhanced
         result%solve_succeeded(iw) = result%solve_info(iw) == 0 .and. .not. result%numerically_singular(iw)
         result%near_singular_collective_pole(iw) = result%solve_info(iw) == 0 .and. &
            result%denominator_condition_number(iw) >= lr_dyson_near_singular_condition .and. &
            .not. result%numerically_singular(iw)

         if (result%solve_info(iw) /= 0) then
            result%frequency_status(iw) = 'BLOCKED: numerically singular Dyson denominator; no regularization applied'
         else if (result%numerically_singular(iw)) then
            result%frequency_status(iw) = 'BLOCKED: singular/precision-limited Dyson denominator; no regularization applied'
         else if (result%near_singular_collective_pole(iw)) then
            result%frequency_status(iw) = 'PASS: solved near-singular collective-pole denominator; no regularization applied'
         else if (result%ill_conditioned(iw)) then
            result%frequency_status(iw) = 'WARNING: solved but Dyson denominator is numerically ill-conditioned'
         else
            result%frequency_status(iw) = 'PASS: Dyson solve'
         end if

         residual = matmul(denominator, enhanced) - request%ks_susceptibility(:, :, iw)
         result%dyson_residual_frobenius(iw) = sqrt(sum(abs(residual)**2))
         result%dyson_residual_relative(iw) = result%dyson_residual_frobenius(iw)/ &
            max(sqrt(sum(abs(request%ks_susceptibility(:, :, iw))**2)), tiny(1.0_rp))
         result%dyson_residual_infinity(iw) = maxval(abs(residual))

         if (result%solve_info(iw) == 0) then
            if (compact) then
               loss = tddft_loss_matrix(enhanced)
            else
               loss = tddft_loss_matrix(space, enhanced)
            end if
         else
            loss = cmplx(0.0_rp, 0.0_rp, rp)
         end if
         result%loss_matrix(:, :, iw) = loss
         if (compact) then
            result%loss_metric_hermiticity_residual(iw) = loss_matrix_hermiticity_residual(loss)
         else
            result%loss_metric_hermiticity_residual(iw) = loss_matrix_hermiticity_residual(space, loss)
         end if
      end do

      if (any(result%solve_info /= 0) .or. any(result%numerically_singular)) then
         result%status = 'BLOCKED: at least one Dyson frequency is numerically singular'
      else if (any(result%near_singular_collective_pole)) then
         result%status = 'PASS: enhanced susceptibility and loss; near-singular poles reported'
      else if (any(result%ill_conditioned)) then
         result%status = 'WARNING: enhanced susceptibility solved; ill-conditioned frequencies reported'
      else
         result%status = 'PASS: enhanced susceptibility and loss matrix'
      end if
   end subroutine evaluate_tddft_dyson

   !> Form D=I-chi_KS*K and solve D*chi=chi_KS for one frequency.
   !>
   !> This standalone primitive is representation agnostic only in the sense
   !> that both inputs must use the same already-canonical representation.  It
   !> does not reconstruct raw kernels and does not apply a response metric.
   module subroutine solve_tddft_dyson_frequency(chi_ks, kernel, chi, denominator, info, condition_number, &
                                          min_singular_value, max_singular_value)
      complex(rp), intent(in) :: chi_ks(:, :), kernel(:, :)
      complex(rp), intent(out) :: chi(:, :), denominator(:, :)
      integer, intent(out) :: info
      real(rp), intent(out), optional :: condition_number, min_singular_value, max_singular_value
      integer :: n, i

      n = size(chi_ks, 1)
      if (n < 1 .or. size(chi_ks, 2) /= n .or. any(shape(kernel) /= [n, n]) .or. &
          any(shape(chi) /= [n, n]) .or. any(shape(denominator) /= [n, n])) then
         error stop 'solve_tddft_dyson_frequency: incompatible canonical response dimensions'
      end if
      denominator = -matmul(chi_ks, kernel)
      do i = 1, n
         denominator(i, i) = denominator(i, i) + cmplx(1.0_rp, 0.0_rp, rp)
      end do
      call solve_tddft_denominator(denominator, chi_ks, chi, info, condition_number, min_singular_value, &
         max_singular_value)
   end subroutine solve_tddft_dyson_frequency

   !> LR-03 loss product for an ordinary orthonormal response matrix.
   !>
   !> The no-space overload is retained for isolated algebraic fixtures.  A
   !> production LR-04 caller must use the overload with `response_space`.
   module function tddft_loss_matrix_ordinary(chi) result(loss)
      complex(rp), intent(in) :: chi(:, :)
      complex(rp) :: loss(size(chi, 1), size(chi, 2))
      integer :: n, i, j

      n = size(chi, 1)
      if (size(chi, 2) /= n .or. n < 1) error stop 'tddft_loss_matrix: chi must be nonempty and square'
      do j = 1, n
         do i = 1, n
            loss(i, j) = -(chi(i, j) - conjg(chi(j, i)))/cmplx(0.0_rp, 2.0_rp*pi, rp)
         end do
      end do
   end function tddft_loss_matrix_ordinary

   !> LR-03 loss product in canonical LR-04 storage.
   !>
   !> `response_operator_adjoint` is the canonical metric-adjoint.  The
   !> origin row/column is set to the zero active-space extension because its
   !> volume measure is exactly zero and it has no physical loss weight.
   module function tddft_loss_matrix_metric(space, chi) result(loss)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: chi(:, :)
      complex(rp) :: loss(size(chi, 1), size(chi, 2))
      complex(rp) :: adjoint(size(chi, 1), size(chi, 2))
      integer :: n, i, j

      n = space%ndim
      call response_operator_adjoint(space, chi, adjoint)
      loss = -(chi - adjoint)/cmplx(0.0_rp, 2.0_rp*pi, rp)
      do i = 1, n
         if (space%metric_weights(i) == 0.0_rp) then
            do j = 1, n
               loss(i, j) = cmplx(0.0_rp, 0.0_rp, rp)
               loss(j, i) = cmplx(0.0_rp, 0.0_rp, rp)
            end do
         end if
      end do
   end function tddft_loss_matrix_metric

   module pure real(rp) function loss_matrix_hermiticity_residual_ordinary(loss) result(residual)
      complex(rp), intent(in) :: loss(:, :)
      integer :: n

      n = size(loss, 1)
      if (size(loss, 2) /= n .or. n < 1) error stop 'loss_matrix_hermiticity_residual: loss must be square'
      residual = maxval(abs(loss - conjg(transpose(loss))))
   end function loss_matrix_hermiticity_residual_ordinary

   module real(rp) function loss_matrix_hermiticity_residual_metric(space, loss) result(residual)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: loss(:, :)
      complex(rp) :: adjoint(size(loss, 1), size(loss, 2))

      call response_operator_adjoint(space, loss, adjoint)
      residual = maxval(abs(loss - adjoint))
   end function loss_matrix_hermiticity_residual_metric

   subroutine validate_dyson_request(request, space)
      type(tddft_dyson_request), intent(in) :: request
      type(response_space_layout), intent(in) :: space
      integer :: n, nw

      if (space%ndim < 1 .or. .not. allocated(space%metric_weights) .or. &
          size(space%metric_weights) /= space%ndim) then
         error stop 'evaluate_tddft_dyson: response space is not initialized'
      end if
      if (any(space%metric_weights < 0.0_rp)) error stop 'evaluate_tddft_dyson: response metric is invalid'
      if (.not. allocated(request%frequencies) .or. size(request%frequencies) < 1) then
         error stop 'evaluate_tddft_dyson: at least one frequency is required'
      end if
      if (request%eta <= 0.0_rp) error stop 'evaluate_tddft_dyson: retarded eta must be positive'
      if (len_trim(request%channel) == 0) error stop 'evaluate_tddft_dyson: circular channel provenance is required'
      select case (trim(request%channel))
      case ('plus', 'minus', 'chi_plus', 'chi_minus')
      case default
         error stop 'evaluate_tddft_dyson: channel must be plus/minus or chi_plus/chi_minus'
      end select
      select case (trim(request%interaction_route))
      case (lr_dyson_route_direct_alsda)
         if (index(trim(request%interaction_provenance), 'KXC-01') /= 1) then
            error stop 'evaluate_tddft_dyson: direct_alsda requires KXC-01 interaction provenance'
         end if
      case (lr_dyson_route_goldstone_sumrule)
         if (index(trim(request%interaction_provenance), 'GSR-01') /= 1) then
            error stop 'evaluate_tddft_dyson: goldstone_sumrule requires GSR-01 interaction provenance'
         end if
      case (lr_dyson_route_direct_alsda_goldstone_corrected)
         if (index(trim(request%interaction_provenance), 'GCR-01') /= 1) then
            error stop 'evaluate_tddft_dyson: corrected route requires GCR-01 interaction provenance'
         end if
      case default
         error stop 'evaluate_tddft_dyson: interaction route is not an explicit supported choice'
      end select
      if (len_trim(request%interaction_provenance) == 0) then
         error stop 'evaluate_tddft_dyson: interaction provenance is required'
      end if
      if (len_trim(request%electronic_state_provenance) == 0) then
         error stop 'evaluate_tddft_dyson: electronic-state provenance is required'
      end if
      if (.not. allocated(request%ks_susceptibility) .or. .not. allocated(request%canonical_interaction)) then
         error stop 'evaluate_tddft_dyson: canonical chiKS and interaction are required'
      end if
      n = space%ndim
      nw = size(request%frequencies)
      if (any(shape(request%ks_susceptibility) /= [n, n, nw])) then
         error stop 'evaluate_tddft_dyson: canonical chiKS/frequency shape mismatch'
      end if
      if (any(shape(request%canonical_interaction) /= [n, n])) then
         error stop 'evaluate_tddft_dyson: canonical interaction shape mismatch'
      end if
   end subroutine validate_dyson_request

   subroutine validate_compact_request(request, n)
      type(tddft_dyson_request), intent(in) :: request
      integer, intent(out) :: n

      if (.not. allocated(request%frequencies) .or. size(request%frequencies) < 1) then
         error stop 'evaluate_tddft_dyson: at least one frequency is required'
      end if
      if (request%eta <= 0.0_rp) error stop 'evaluate_tddft_dyson: retarded eta must be positive'
      if (len_trim(request%channel) == 0) error stop 'evaluate_tddft_dyson: circular channel provenance is required'
      select case (trim(request%channel))
      case ('plus', 'minus', 'chi_plus', 'chi_minus')
      case default
         error stop 'evaluate_tddft_dyson: channel must be plus/minus or chi_plus/chi_minus'
      end select
      select case (trim(request%interaction_route))
      case (lr_dyson_route_direct_alsda)
         if (index(trim(request%interaction_provenance), 'KXC-01') /= 1) then
            error stop 'evaluate_tddft_dyson: direct_alsda requires KXC-01 interaction provenance'
         end if
      case (lr_dyson_route_goldstone_sumrule)
         if (index(trim(request%interaction_provenance), 'GSR-01') /= 1) then
            error stop 'evaluate_tddft_dyson: goldstone_sumrule requires GSR-01 interaction provenance'
         end if
      case (lr_dyson_route_direct_alsda_goldstone_corrected)
         if (index(trim(request%interaction_provenance), 'GCR-01') /= 1) then
            error stop 'evaluate_tddft_dyson: corrected route requires GCR-01 interaction provenance'
         end if
      case default
         error stop 'evaluate_tddft_dyson: interaction route is not an explicit supported choice'
      end select
      if (len_trim(request%interaction_provenance) == 0) then
         error stop 'evaluate_tddft_dyson: interaction provenance is required'
      end if
      if (len_trim(request%electronic_state_provenance) == 0) then
         error stop 'evaluate_tddft_dyson: electronic-state provenance is required'
      end if
      if (.not. allocated(request%ks_susceptibility) .or. .not. allocated(request%canonical_interaction)) then
         error stop 'evaluate_tddft_dyson: compact chiKS and interaction are required'
      end if
      n = size(request%canonical_interaction, 1)
      if (n < 1 .or. size(request%canonical_interaction, 2) /= n) then
         error stop 'evaluate_tddft_dyson: compact interaction must be nonempty and square'
      end if
      if (any(shape(request%ks_susceptibility) /= [n, n, size(request%frequencies)])) then
         error stop 'evaluate_tddft_dyson: compact chiKS/frequency shape mismatch'
      end if
   end subroutine validate_compact_request

   !> Solve one already-formed Dyson denominator and report its singular values.
   subroutine solve_tddft_denominator(denominator, rhs, solution, info, condition_number, min_singular_value, &
                                      max_singular_value, numerically_singular, ill_conditioned)
      complex(rp), intent(in) :: denominator(:, :), rhs(:, :)
      complex(rp), intent(out) :: solution(:, :)
      integer, intent(out) :: info
      real(rp), intent(out), optional :: condition_number, min_singular_value, max_singular_value
      logical, intent(out), optional :: numerically_singular, ill_conditioned
      complex(rp), allocatable :: system(:, :), rhs_work(:, :)
      integer, allocatable :: pivots(:)
      real(rp) :: cond, smin, smax
      logical :: precision_limited, poorly_conditioned
      integer :: n
      external :: zgesv

      n = size(denominator, 1)
      if (n < 1 .or. size(denominator, 2) /= n .or. any(shape(rhs) /= [n, n]) .or. &
          any(shape(solution) /= [n, n])) then
         error stop 'solve_tddft_denominator: incompatible matrix dimensions'
      end if

      allocate(system(n, n), rhs_work(n, n), pivots(n))
      system = denominator
      rhs_work = rhs
      call singular_value_diagnostics(system, cond, smin, smax, info)
      precision_limited = info /= 0 .or. smax == 0.0_rp .or. smin <= &
         100.0_rp*epsilon(1.0_rp)*smax
      poorly_conditioned = cond >= lr_dyson_ill_conditioned_condition

      ! The SVD above operates on a private copy.  zgesv therefore receives the
      ! original denominator and all RHS columns, not a reconstructed inverse.
      system = denominator
      if (precision_limited) then
         info = 1
         solution = cmplx(0.0_rp, 0.0_rp, rp)
      else
         call zgesv(n, n, system, n, pivots, rhs_work, n, info)
         if (info == 0) then
            solution = rhs_work
         else
            solution = cmplx(0.0_rp, 0.0_rp, rp)
         end if
      end if

      if (present(condition_number)) condition_number = cond
      if (present(min_singular_value)) min_singular_value = smin
      if (present(max_singular_value)) max_singular_value = smax
      if (present(numerically_singular)) numerically_singular = precision_limited
      if (present(ill_conditioned)) ill_conditioned = poorly_conditioned
   end subroutine solve_tddft_denominator

   subroutine singular_value_diagnostics(matrix, condition_number, min_singular_value, max_singular_value, info)
      complex(rp), intent(in) :: matrix(:, :)
      real(rp), intent(out) :: condition_number, min_singular_value, max_singular_value
      integer, intent(out) :: info
      complex(rp), allocatable :: work_matrix(:,:), work(:)
      complex(rp) :: work_query(1), u_dummy(1, 1), vt_dummy(1, 1)
      real(rp), allocatable :: singular_values(:), rwork(:)
      integer :: n, lwork, svd_info
      external :: zgesvd

      n = size(matrix, 1)
      if (size(matrix, 2) /= n .or. n < 1) error stop 'singular_value_diagnostics: matrix must be nonempty and square'
      allocate(work_matrix(n, n), singular_values(n), rwork(5*n))
      work_matrix = matrix
      call zgesvd('N', 'N', n, n, work_matrix, n, singular_values, u_dummy, 1, vt_dummy, 1, work_query, -1, rwork, svd_info)
      if (svd_info /= 0) then
         info = svd_info
         condition_number = huge(1.0_rp)
         min_singular_value = 0.0_rp
         max_singular_value = 0.0_rp
         return
      end if
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      work_matrix = matrix
      call zgesvd('N', 'N', n, n, work_matrix, n, singular_values, u_dummy, 1, vt_dummy, 1, work, lwork, rwork, svd_info)
      info = svd_info
      if (info /= 0) then
         condition_number = huge(1.0_rp)
         min_singular_value = 0.0_rp
         max_singular_value = 0.0_rp
         return
      end if
      max_singular_value = maxval(singular_values)
      min_singular_value = minval(singular_values)
      if (min_singular_value > 0.0_rp) then
         condition_number = max_singular_value/min_singular_value
      else
         condition_number = huge(1.0_rp)
      end if
   end subroutine singular_value_diagnostics

   !> Return the smallest eigenvalue magnitude of a denominator for the raw
   !> static/pole diagnostic.  This is diagnostic only; the solve continues to
   !> use the unmodified denominator and LAPACK zgesv.
   subroutine minimum_magnitude_eigenvalue(matrix, minimum_magnitude)
      complex(rp), intent(in) :: matrix(:, :)
      real(rp), intent(out) :: minimum_magnitude
      complex(rp), allocatable :: work_matrix(:, :), eigenvalues(:), work(:)
      complex(rp) :: work_query(1), vl_dummy(1, 1), vr_dummy(1, 1)
      real(rp), allocatable :: rwork(:)
      integer :: n, lwork, info
      external :: zgeev

      n = size(matrix, 1)
      if (size(matrix, 2) /= n .or. n < 1) error stop 'minimum_magnitude_eigenvalue: matrix must be nonempty and square'
      allocate(work_matrix(n, n), eigenvalues(n), rwork(max(1, 2*n)))
      work_matrix = matrix
      call zgeev('N', 'N', n, work_matrix, n, eigenvalues, vl_dummy, 1, vr_dummy, 1, work_query, -1, rwork, info)
      if (info /= 0) then
         minimum_magnitude = huge(1.0_rp)
         deallocate(work_matrix, eigenvalues, rwork)
         return
      end if
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      work_matrix = matrix
      call zgeev('N', 'N', n, work_matrix, n, eigenvalues, vl_dummy, 1, vr_dummy, 1, work, lwork, rwork, info)
      if (info == 0) then
         minimum_magnitude = minval(abs(eigenvalues))
      else
         minimum_magnitude = huge(1.0_rp)
      end if
      deallocate(work_matrix, eigenvalues, work, rwork)
   end subroutine minimum_magnitude_eigenvalue

   !> Public diagnostic seam for compact/site-space callers.  The solve path
   !> remains owned by this module; callers do not duplicate the LAPACK
   !> eigenvalue diagnostic or alter the denominator before measuring it.
   module subroutine tddft_denominator_minimum_magnitude_eigenvalue(matrix, minimum_magnitude)
      complex(rp), intent(in) :: matrix(:, :)
      real(rp), intent(out) :: minimum_magnitude

      call minimum_magnitude_eigenvalue(matrix, minimum_magnitude)
   end subroutine tddft_denominator_minimum_magnitude_eigenvalue
end submodule linear_response_kernel_dyson
