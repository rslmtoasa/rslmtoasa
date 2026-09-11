!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Direct local ALSDA transverse XC kernel for the LR-03 baseline.
!>
!> The kernel is the explicitly named mixed Pauli-response approximation
!>
!>   K_xc^(P<-SR)(r) = B_xc^sigma,SR(r) / s^P(r),
!>
!> with B_xc^sigma,SR=(V_xc,up-V_xc,down)/2.  The accepted LR-01 radial
!> snapshot supplies the XC arrays and provenance; the Pauli magnetization is
!> supplied explicitly by the already certified Pauli response-space closure.
!>
!> This module owns no Goldstone correction, sum-rule solve, Dyson equation,
!> empirical rescaling, or site reduction.
!------------------------------------------------------------------------------
module lr_alsda_kernel_mod

   use precision_mod, only: rp
   use radial_ground_state_mod, only: radial_ground_state, radial_xc_provenance
   use lr_response_space_mod, only: response_space_layout, response_apply_operator, response_local_operator, &
      response_vector_norm, response_rigid_vector_overlap
   implicit none
   private

   character(len=*), parameter, public :: lr_kxc_magnetization_pauli = 'pauli_projected'
   character(len=*), parameter, public :: lr_kxc_units = 'Ry bohr^3'
   character(len=*), parameter, public :: lr_kxc_representation = 'LR-04 canonical local operator'
   real(rp), parameter, public :: lr_kxc_default_low_m_relative = 1.0e-8_rp

   !> Request for one direct ALSDA kernel on a complete response layout.
   !>
   !> `pauli_magnetization(site,radial_point)` is the Pauli/no-SOC response
   !> magnetization selected by LR-03.  It is not reconstructed from a
   !> compressed potential field and it is not silently replaced by the SR
   !> density in the accepted radial snapshot.
   type, public :: lr_alsda_kernel_request
      type(response_space_layout), pointer :: response_space => null()
      type(radial_ground_state), pointer :: ground_states(:) => null()
      real(rp), allocatable :: pauli_magnetization(:, :)
      character(len=32) :: magnetization_label = lr_kxc_magnetization_pauli
      character(len=512) :: requested_functional = ''
      character(len=32) :: requested_backend = ''
      integer :: requested_txc = -1
      ! This is a diagnostic threshold only.  It never changes the ratio.
      real(rp) :: low_m_diagnostic_relative = lr_kxc_default_low_m_relative
   end type lr_alsda_kernel_request

   !> Low-m diagnostic for the active positive-measure radial grid.
   type, public :: lr_alsda_low_m_diagnostic
      integer :: active_radial_points = 0
      integer :: low_m_points = 0
      integer :: exact_zero_active_points = 0
      integer :: null_measure_zero_points = 0
      real(rp) :: max_abs_magnetization = 0.0_rp
      real(rp) :: min_abs_magnetization = huge(1.0_rp)
      real(rp) :: min_relative_abs_magnetization = huge(1.0_rp)
      real(rp) :: diagnostic_relative_threshold = lr_kxc_default_low_m_relative
   end type lr_alsda_low_m_diagnostic

   !> Direct ALSDA result and provenance carried to later interaction code.
   type, public :: lr_alsda_kernel_result
      real(rp), allocatable :: pointwise_kernel(:, :) ! (site,radial point)
      complex(rp), allocatable :: canonical_operator(:, :) ! LR-04 B=K_raw*W
      type(lr_alsda_low_m_diagnostic) :: low_m
      type(radial_xc_provenance) :: xc_provenance
      character(len=32) :: magnetization_label = ''
      character(len=64) :: units = lr_kxc_units
      character(len=128) :: response_representation = lr_kxc_representation
      logical :: origin_null_measure_extension = .false.
   end type lr_alsda_kernel_result

   public :: evaluate_lr_alsda_kernel
   public :: evaluate_lr_alsda_static_residual
   public :: lr_alsda_provenance_matches

   interface evaluate_lr_alsda_kernel
      module procedure evaluate_lr_alsda_kernel_request
      module procedure evaluate_lr_alsda_kernel_explicit
   end interface evaluate_lr_alsda_kernel

contains

   !> Evaluate the request-form direct ALSDA kernel.
   subroutine evaluate_lr_alsda_kernel_request(request, result)
      type(lr_alsda_kernel_request), intent(in) :: request
      type(lr_alsda_kernel_result), intent(out) :: result

      if (.not. associated(request%response_space)) then
         error stop 'evaluate_lr_alsda_kernel: response-space reference is incomplete'
      end if
      if (.not. associated(request%ground_states)) then
         error stop 'evaluate_lr_alsda_kernel: accepted radial-state references are incomplete'
      end if
      if (.not. allocated(request%pauli_magnetization)) then
         error stop 'evaluate_lr_alsda_kernel: Pauli magnetization is required explicitly'
      end if
      call evaluate_lr_alsda_kernel_explicit(request%response_space, request%ground_states, &
         request%pauli_magnetization, result, request%requested_functional, request%requested_backend, &
         request%requested_txc, request%magnetization_label, request%low_m_diagnostic_relative)
   end subroutine evaluate_lr_alsda_kernel_request

   !> Explicit form retained for callers that keep the three contracts apart.
   subroutine evaluate_lr_alsda_kernel_explicit(space, ground_states, pauli_magnetization, result, &
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
      if (size(ground_states) /= space%nsite) then
         error stop 'evaluate_lr_alsda_kernel: one accepted radial state is required per response site'
      end if
      if (any(shape(pauli_magnetization) /= [space%nsite, space%npoint])) then
         error stop 'evaluate_lr_alsda_kernel: Pauli magnetization shape mismatch'
      end if
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
   logical function lr_alsda_provenance_matches(accepted, requested_functional, requested_backend, requested_txc) &
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
   subroutine evaluate_lr_alsda_static_residual(space, static_susceptibility, canonical_kernel, rigid_vector, &
                                                absolute_residual, relative_residual, rigid_overlap, residual_vector)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: static_susceptibility(:, :), canonical_kernel(:, :), rigid_vector(:)
      real(rp), intent(out) :: absolute_residual, relative_residual
      complex(rp), intent(out), optional :: rigid_overlap
      complex(rp), intent(out), optional :: residual_vector(:)
      complex(rp), allocatable :: field(:), response(:), residual(:)
      real(rp) :: target_norm

      if (any(shape(static_susceptibility) /= [space%ndim, space%ndim]) .or. &
          any(shape(canonical_kernel) /= [space%ndim, space%ndim])) then
         error stop 'evaluate_lr_alsda_static_residual: operator shape mismatch'
      end if
      if (size(rigid_vector) /= space%ndim) then
         error stop 'evaluate_lr_alsda_static_residual: rigid-vector shape mismatch'
      end if
      if (present(residual_vector) .and. size(residual_vector) /= space%ndim) then
         error stop 'evaluate_lr_alsda_static_residual: residual-vector shape mismatch'
      end if
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
      if (size(space%radius) /= space%npoint .or. size(space%radial_weights) /= space%npoint .or. &
          size(space%metric_weights) /= space%ndim) then
         error stop 'evaluate_lr_alsda_kernel: response-space metadata is inconsistent'
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

end module lr_alsda_kernel_mod
