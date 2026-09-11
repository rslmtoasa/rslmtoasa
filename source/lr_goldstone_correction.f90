!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Optional Buczek-Ernst-Sandratskii Goldstone eigenvalue correction.
!>
!> The route is deliberately a consumer of canonical LR-04 operators.  It
!> forms D=I-chiKS(0)Kxc, identifies one rigid-rotation eigenmode, replaces
!> only that eigenvalue by zero, and reconstructs
!>
!>   Kxc_corr=chiKS(0)^(-1)(I-Dcorr).
!>
!> No raw pointwise matrix is multiplied directly and no correction is enabled
!> unless the caller explicitly selects it and supplies the convergence,
!> provenance, static-point, field, and eta/k-mesh consistency gates.
!------------------------------------------------------------------------------
module lr_goldstone_correction_mod

   use precision_mod, only: rp
   use lr_response_space_mod, only: response_space_layout, response_apply_operator, &
      response_compose_operators, response_identity_operator, response_vector_norm, &
      response_rigid_vector_overlap
   implicit none
   private

   character(len=*), parameter, public :: lr_gcr_route = 'goldstone_eigenvalue_correction'
   character(len=*), parameter, public :: lr_gcr_units = 'Ry bohr^3'
   character(len=*), parameter, public :: lr_gcr_representation = 'LR-04 canonical response operator'
   character(len=*), parameter, public :: lr_gcr_kxc_provenance = 'KXC-01 direct ALSDA LR-03'
   character(len=*), parameter, public :: lr_gcr_class_hermitian = 'Hermitian'
   character(len=*), parameter, public :: lr_gcr_class_metric_hermitian = 'metric-Hermitian'
   character(len=*), parameter, public :: lr_gcr_class_nonhermitian = 'complex non-Hermitian'
   character(len=*), parameter, public :: lr_gcr_solver_hermitian = 'zheev'
   character(len=*), parameter, public :: lr_gcr_solver_metric = 'zheev on W**(1/2) D W**(-1/2)'
   character(len=*), parameter, public :: lr_gcr_solver_nonhermitian = 'zgeev plus right-eigenvector inverse'

   real(rp), parameter, public :: lr_gcr_matrix_tolerance = 1.0e-10_rp
   real(rp), parameter, public :: lr_gcr_isolation_ratio = 0.25_rp
   real(rp), parameter, public :: lr_gcr_min_rigid_overlap = 0.80_rp
   integer, parameter, public :: lr_gcr_min_consistency_samples = 2

   !> One previously evaluated static D diagnostic used by the consistency
   !> gate.  The sample must identify the same isolated, rigid mode at a modest
   !> eta or k-mesh change.
   type, public :: lr_goldstone_consistency_sample
      real(rp) :: eta = -1.0_rp
      integer :: k_mesh_points = 0
      complex(rp) :: goldstone_eigenvalue = cmplx(0.0_rp, 0.0_rp, rp)
      real(rp) :: competing_eigenvalue_abs = -1.0_rp
      real(rp) :: normalized_rigid_overlap = -1.0_rp
   end type lr_goldstone_consistency_sample

   !> Explicit request for the optional correction route.
   type, public :: lr_goldstone_correction_request
      type(response_space_layout), pointer :: response_space => null()
      complex(rp), allocatable :: static_susceptibility(:, :) ! canonical chiKS(q=0,omega=0)
      complex(rp), allocatable :: xc_kernel(:, :) ! canonical Kxc to be corrected
      complex(rp), allocatable :: rigid_rotation_vector(:) ! LR-03 vector in point coordinates

      ! False is the safe default.  A disabled request is a clean no-op.
      logical :: enabled = .false.
      logical :: lr06_converged = .false.
      logical :: static_q_zero = .false.
      real(rp) :: static_frequency = huge(1.0_rp)
      real(rp) :: eta = -1.0_rp
      integer :: k_mesh_points = 0
      type(lr_goldstone_consistency_sample), allocatable :: consistency_samples(:)

      ! The direct KXC-01 route must pass its provenance contract explicitly.
      character(len=64) :: kernel_units = ''
      character(len=128) :: kernel_provenance = ''
      logical :: kernel_provenance_verified = .false.

      ! These flags are explicit because a physically gapped mode must not be
      ! mistaken for a numerical Goldstone defect.
      logical :: soc_present = .false.
      logical :: external_field_present = .false.
      logical :: constraining_field_present = .false.
   end type lr_goldstone_correction_request

   !> Raw/corrected objects and all correction gates/diagnostics.
   type, public :: lr_goldstone_correction_result
      logical :: selected = .false.
      logical :: corrected = .false.
      logical :: blocked = .true.
      integer :: active_dimension = 0
      integer :: goldstone_index = 0
      real(rp) :: ordinary_hermitian_error = huge(1.0_rp)
      real(rp) :: metric_hermitian_error = huge(1.0_rp)
      real(rp) :: goldstone_eigenvalue_abs = huge(1.0_rp)
      real(rp) :: goldstone_isolation_ratio = huge(1.0_rp)
      real(rp) :: rigid_overlap = 0.0_rp
      real(rp) :: raw_rigid_residual = huge(1.0_rp)
      real(rp) :: corrected_rigid_residual = huge(1.0_rp)
      real(rp) :: kernel_reconstruction_residual = huge(1.0_rp)
      real(rp) :: untouched_eigenvalue_change = huge(1.0_rp)

      character(len=32) :: route = lr_gcr_route
      character(len=64) :: units = lr_gcr_units
      character(len=128) :: response_representation = lr_gcr_representation
      character(len=64) :: representation_class = ''
      character(len=128) :: eigensolver = ''
      character(len=256) :: status = 'not evaluated'

      integer, allocatable :: active_indices(:)
      complex(rp), allocatable :: raw_denominator(:, :)
      complex(rp), allocatable :: corrected_denominator(:, :)
      complex(rp), allocatable :: raw_kernel(:, :)
      complex(rp), allocatable :: corrected_kernel(:, :)
      complex(rp), allocatable :: eigenvalues(:)
      complex(rp), allocatable :: corrected_eigenvalues(:)
      complex(rp), allocatable :: right_eigenvectors(:, :)
      complex(rp), allocatable :: left_eigenvectors(:, :)
      complex(rp), allocatable :: eigenvector_inverse(:, :)
   end type lr_goldstone_correction_result

   public :: evaluate_lr_goldstone_correction

contains

   !> Evaluate the explicitly selected BES correction.
   subroutine evaluate_lr_goldstone_correction(request, result)
      type(lr_goldstone_correction_request), intent(in) :: request
      type(lr_goldstone_correction_result), intent(out) :: result

      type(response_space_layout), pointer :: space
      complex(rp), allocatable :: product(:, :), identity(:, :), d_full(:, :), d_active(:, :)
      complex(rp), allocatable :: chi_active(:, :), d_corrected_active(:, :)
      complex(rp), allocatable :: inverse_rhs(:, :), reconstruction(:, :)
      complex(rp), allocatable :: rigid_active(:), rigid_full(:), d_rigid(:), corrected_rigid(:)
      complex(rp), allocatable :: sqrtw(:), hermitian_form(:, :), eigenvectors_work(:, :)
      integer :: ndim, nactive, i, j, mode, info
      logical :: accepted, solve_ok
      character(len=256) :: refusal

      call initialize_result(result)
      if (.not. request%enabled) then
         result%blocked = .false.
         result%status = 'DISABLED: BES correction was not explicitly selected'
         return
      end if
      result%selected = .true.

      if (.not. associated(request%response_space)) then
         error stop 'evaluate_lr_goldstone_correction: response-space reference is incomplete'
      end if
      if (.not. allocated(request%static_susceptibility) .or. .not. allocated(request%xc_kernel) .or. &
          .not. allocated(request%rigid_rotation_vector)) then
         error stop 'evaluate_lr_goldstone_correction: canonical chiKS, Kxc, and rigid vector are required'
      end if
      space => request%response_space
      ndim = space%ndim
      call validate_shapes(space, request%static_susceptibility, request%xc_kernel, &
         request%rigid_rotation_vector)

      ! These are physics/campaign blockers, not programming errors.  They
      ! return the raw-input diagnosis without fabricating a corrected route.
      call validate_selection_gates(request, accepted, refusal)
      if (.not. accepted) then
         result%status = trim(refusal)
         return
      end if

      nactive = count(space%metric_weights > 0.0_rp)
      if (nactive < 2) then
         result%status = 'BLOCKED: at least two positive-measure response coordinates are required'
         return
      end if
      result%active_dimension = nactive
      allocate(result%active_indices(nactive), d_active(nactive, nactive), chi_active(nactive, nactive), &
         rigid_active(nactive), rigid_full(ndim), d_rigid(ndim), &
         corrected_rigid(ndim))
      call build_active_indices(space, result%active_indices)

      allocate(product(ndim, ndim), identity(ndim, ndim), d_full(ndim, ndim))
      ! Both inputs are canonical LR-04 operators.  This composition is the
      ! only multiplication convention used to form D; no raw W is inserted.
      call response_compose_operators(space, request%static_susceptibility, request%xc_kernel, product)
      call response_identity_operator(space, identity)
      d_full = identity - product
      result%raw_denominator = d_full
      result%raw_kernel = request%xc_kernel
      d_active = d_full(result%active_indices, result%active_indices)
      chi_active = request%static_susceptibility(result%active_indices, result%active_indices)
      rigid_active = request%rigid_rotation_vector(result%active_indices)
      rigid_full = cmplx(0.0_rp, 0.0_rp, rp)
      rigid_full(result%active_indices) = rigid_active
      d_rigid = cmplx(0.0_rp, 0.0_rp, rp)
      call response_apply_operator(space, result%raw_denominator, rigid_full, d_rigid)
      result%raw_rigid_residual = response_vector_norm(space, d_rigid)

      allocate(result%eigenvalues(nactive), result%corrected_eigenvalues(nactive), &
         result%right_eigenvectors(nactive, nactive), result%left_eigenvectors(nactive, nactive), &
         result%eigenvector_inverse(nactive, nactive), result%corrected_denominator(ndim, ndim), &
         result%corrected_kernel(ndim, ndim))
      result%corrected_denominator = cmplx(0.0_rp, 0.0_rp, rp)
      result%corrected_kernel = cmplx(0.0_rp, 0.0_rp, rp)
      call classify_denominator(space, result%active_indices, d_active, result%ordinary_hermitian_error, &
         result%metric_hermitian_error, result%representation_class)
      select case (trim(result%representation_class))
      case (lr_gcr_class_hermitian)
         result%eigensolver = lr_gcr_solver_hermitian
         allocate(eigenvectors_work(nactive, nactive))
         call diagonalize_hermitian(d_active, result%eigenvalues, eigenvectors_work, info)
         if (info /= 0) then
            result%status = 'BLOCKED: Hermitian eigensolver failed'
            return
         end if
         result%right_eigenvectors = eigenvectors_work
         result%left_eigenvectors = eigenvectors_work
      case (lr_gcr_class_metric_hermitian)
         result%eigensolver = lr_gcr_solver_metric
         allocate(sqrtw(nactive), hermitian_form(nactive, nactive), eigenvectors_work(nactive, nactive))
         do i = 1, nactive
            sqrtw(i) = sqrt(space%metric_weights(result%active_indices(i)))
         end do
         do i = 1, nactive
            do j = 1, nactive
               hermitian_form(i, j) = sqrtw(i)*d_active(i, j)/sqrtw(j)
            end do
         end do
         call diagonalize_hermitian(hermitian_form, result%eigenvalues, eigenvectors_work, info)
         if (info /= 0) then
            result%status = 'BLOCKED: metric-Hermitian eigensolver failed'
            return
         end if
         do j = 1, nactive
            result%right_eigenvectors(:, j) = eigenvectors_work(:, j)/sqrtw
            result%left_eigenvectors(:, j) = space%metric_weights(result%active_indices)* &
               result%right_eigenvectors(:, j)
         end do
      case default
         result%representation_class = lr_gcr_class_nonhermitian
         result%eigensolver = lr_gcr_solver_nonhermitian
         call diagonalize_nonhermitian(d_active, result%eigenvalues, result%left_eigenvectors, &
            result%right_eigenvectors, info)
         if (info /= 0) then
            result%status = 'BLOCKED: complex non-Hermitian eigensolver failed'
            return
         end if
      end select

      call invert_eigenvector_matrix(result%right_eigenvectors, result%eigenvector_inverse, solve_ok)
      if (.not. solve_ok) then
         result%status = 'BLOCKED: eigenvector matrix could not be inverted for BES reconstruction'
         return
      end if
      result%corrected_eigenvalues = result%eigenvalues

      mode = minloc(abs(result%eigenvalues), dim=1)
      result%goldstone_index = mode
      result%goldstone_eigenvalue_abs = abs(result%eigenvalues(mode))
      result%goldstone_isolation_ratio = isolation_ratio(result%eigenvalues, mode)
      result%rigid_overlap = normalized_metric_overlap(space, result%active_indices, &
         result%right_eigenvectors(:, mode), request%rigid_rotation_vector)
      call validate_mode_identification(result, accepted, refusal)
      if (.not. accepted) then
         result%status = trim(refusal)
         return
      end if
      call validate_consistency(request, accepted, refusal)
      if (.not. accepted) then
         result%status = trim(refusal)
         return
      end if

      result%corrected_eigenvalues = result%eigenvalues
      result%corrected_eigenvalues(mode) = cmplx(0.0_rp, 0.0_rp, rp)
      result%untouched_eigenvalue_change = 0.0_rp
      do i = 1, nactive
         if (i /= mode) result%untouched_eigenvalue_change = max(result%untouched_eigenvalue_change, &
            abs(result%corrected_eigenvalues(i) - result%eigenvalues(i)))
      end do
      allocate(d_corrected_active(nactive, nactive), inverse_rhs(nactive, nactive), &
         reconstruction(nactive, nactive))
      call reconstruct_from_eigensystem(result%right_eigenvectors, result%eigenvector_inverse, &
         result%corrected_eigenvalues, d_corrected_active)
      result%corrected_denominator = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, ndim
         result%corrected_denominator(i, i) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
      result%corrected_denominator(result%active_indices, result%active_indices) = d_corrected_active

      ! BES Eq. (39): chiKS*Kcorr=I-Dcorr.  The solve, rather than an
      ! explicit inverse, is used to avoid an unnecessary inverse operation.
      inverse_rhs = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, nactive
         inverse_rhs(i, i) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
      inverse_rhs = inverse_rhs - d_corrected_active
      call solve_matrix(chi_active, inverse_rhs, solve_ok)
      if (.not. solve_ok) then
         result%status = 'BLOCKED: static chiKS is singular; corrected Kxc is undefined'
         return
      end if
      result%corrected_kernel = cmplx(0.0_rp, 0.0_rp, rp)
      result%corrected_kernel(result%active_indices, result%active_indices) = inverse_rhs
      reconstruction = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, nactive
         reconstruction(i, i) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
      reconstruction = reconstruction - matmul(chi_active, inverse_rhs)
      result%kernel_reconstruction_residual = maxval(abs(reconstruction - d_corrected_active))

      corrected_rigid = cmplx(0.0_rp, 0.0_rp, rp)
      call response_apply_operator(space, result%corrected_denominator, rigid_full, corrected_rigid)
      result%corrected_rigid_residual = response_vector_norm(space, corrected_rigid)

      if (result%kernel_reconstruction_residual > lr_gcr_matrix_tolerance) then
         result%status = 'BLOCKED: reconstructed Kxc does not reproduce Dcorr'
         return
      end if
      result%corrected = .true.
      result%blocked = .false.
      result%status = 'PASS: optional BES Goldstone eigenvalue correction'
   end subroutine evaluate_lr_goldstone_correction

   subroutine initialize_result(result)
      type(lr_goldstone_correction_result), intent(out) :: result

      result%selected = .false.
      result%corrected = .false.
      result%blocked = .true.
      result%active_dimension = 0
      result%goldstone_index = 0
      result%ordinary_hermitian_error = huge(1.0_rp)
      result%metric_hermitian_error = huge(1.0_rp)
      result%goldstone_eigenvalue_abs = huge(1.0_rp)
      result%goldstone_isolation_ratio = huge(1.0_rp)
      result%rigid_overlap = 0.0_rp
      result%raw_rigid_residual = huge(1.0_rp)
      result%corrected_rigid_residual = huge(1.0_rp)
      result%kernel_reconstruction_residual = huge(1.0_rp)
      result%untouched_eigenvalue_change = huge(1.0_rp)
      result%route = lr_gcr_route
      result%units = lr_gcr_units
      result%response_representation = lr_gcr_representation
      result%representation_class = ''
      result%eigensolver = ''
      result%status = 'not evaluated'
   end subroutine initialize_result

   subroutine validate_shapes(space, susceptibility, kernel, rigid)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: susceptibility(:, :), kernel(:, :), rigid(:)

      if (space%ndim < 1 .or. .not. allocated(space%metric_weights)) then
         error stop 'evaluate_lr_goldstone_correction: response space is not initialized'
      end if
      if (space%nchannel /= 1) then
         error stop 'evaluate_lr_goldstone_correction: one certified circular response channel is required'
      end if
      if (any(shape(susceptibility) /= [space%ndim, space%ndim]) .or. &
          any(shape(kernel) /= [space%ndim, space%ndim]) .or. size(rigid) /= space%ndim) then
         error stop 'evaluate_lr_goldstone_correction: canonical operator or rigid-vector shape mismatch'
      end if
   end subroutine validate_shapes

   subroutine validate_selection_gates(request, accepted, refusal)
      type(lr_goldstone_correction_request), intent(in) :: request
      logical, intent(out) :: accepted
      character(len=*), intent(out) :: refusal

      accepted = .false.
      refusal = 'BLOCKED: correction selection gates were not satisfied'
      if (.not. request%lr06_converged) then
         refusal = 'BLOCKED: LR-06 convergence evidence is required'
         return
      end if
      if (.not. request%static_q_zero) then
         refusal = 'BLOCKED: BES D must be formed from the verified q=0 static point'
         return
      end if
      if (.not. finite_real(request%static_frequency) .or. .not. finite_real(request%eta)) then
         refusal = 'BLOCKED: static frequency and eta provenance must be finite'
         return
      end if
      if (abs(request%static_frequency) > lr_gcr_matrix_tolerance) then
         refusal = 'BLOCKED: BES D must be formed at omega=0'
         return
      end if
      if (request%eta < 0.0_rp .or. request%k_mesh_points < 1) then
         refusal = 'BLOCKED: static eta and k-mesh provenance are required'
         return
      end if
      if (request%soc_present .or. request%external_field_present .or. request%constraining_field_present) then
         refusal = 'BLOCKED: SOC/external/constraining field makes a zero mode unphysical'
         return
      end if
      if (trim(request%kernel_units) /= lr_gcr_units) then
         refusal = 'BLOCKED: Kxc units do not match the LR-03/LR-04 contract'
         return
      end if
      if (.not. request%kernel_provenance_verified .or. &
          trim(request%kernel_provenance) /= lr_gcr_kxc_provenance) then
         refusal = 'BLOCKED: Kxc provenance does not match KXC-01 direct ALSDA LR-03'
         return
      end if
      accepted = .true.
   end subroutine validate_selection_gates

   subroutine build_active_indices(space, indices)
      type(response_space_layout), intent(in) :: space
      integer, intent(out) :: indices(:)
      integer :: i, n

      n = 0
      do i = 1, space%ndim
         if (space%metric_weights(i) > 0.0_rp) then
            n = n + 1
            indices(n) = i
         end if
      end do
   end subroutine build_active_indices

   subroutine classify_denominator(space, active_indices, denominator, ordinary_error, metric_error, classification)
      type(response_space_layout), intent(in) :: space
      integer, intent(in) :: active_indices(:)
      complex(rp), intent(in) :: denominator(:, :)
      real(rp), intent(out) :: ordinary_error, metric_error
      character(len=*), intent(out) :: classification
      integer :: i, j, n, flat_i, flat_j
      real(rp) :: ordinary_scale, metric_scale, wi, wj

      n = size(denominator, 1)
      ordinary_error = 0.0_rp
      metric_error = 0.0_rp
      ordinary_scale = max(1.0_rp, maxval(abs(denominator)))
      metric_scale = 0.0_rp
      ! The active denominator is always in the order returned by
      ! build_active_indices.
      do i = 1, n
         flat_i = active_indices(i)
         wi = space%metric_weights(flat_i)
         do j = 1, n
            flat_j = active_indices(j)
            wj = space%metric_weights(flat_j)
            ordinary_error = max(ordinary_error, abs(denominator(i, j) - conjg(denominator(j, i))))
            metric_error = max(metric_error, abs(conjg(denominator(j, i))*wj - wi*denominator(i, j)))
            metric_scale = max(metric_scale, abs(conjg(denominator(j, i))*wj), abs(wi*denominator(i, j)))
         end do
      end do
      ordinary_error = ordinary_error/ordinary_scale
      metric_error = metric_error/max(tiny(1.0_rp), metric_scale)
      if (ordinary_error <= lr_gcr_matrix_tolerance) then
         classification = lr_gcr_class_hermitian
      else if (metric_error <= lr_gcr_matrix_tolerance) then
         classification = lr_gcr_class_metric_hermitian
      else
         classification = lr_gcr_class_nonhermitian
      end if
   end subroutine classify_denominator

   real(rp) function normalized_metric_overlap(space, active_indices, vector, rigid) result(value)
      type(response_space_layout), intent(in) :: space
      integer, intent(in) :: active_indices(:)
      complex(rp), intent(in) :: vector(:), rigid(:)
      complex(rp), allocatable :: full_vector(:)
      real(rp) :: vector_norm, rigid_norm

      allocate(full_vector(space%ndim))
      full_vector = cmplx(0.0_rp, 0.0_rp, rp)
      full_vector(active_indices) = vector
      vector_norm = response_vector_norm(space, full_vector)
      rigid_norm = response_vector_norm(space, rigid)
      if (vector_norm <= tiny(1.0_rp) .or. rigid_norm <= tiny(1.0_rp)) then
         value = 0.0_rp
      else
         value = abs(response_rigid_vector_overlap(space, full_vector, rigid))/(vector_norm*rigid_norm)
      end if
   end function normalized_metric_overlap

   real(rp) function isolation_ratio(eigenvalues, mode) result(value)
      complex(rp), intent(in) :: eigenvalues(:)
      integer, intent(in) :: mode
      real(rp) :: competitor
      integer :: i

      competitor = huge(1.0_rp)
      do i = 1, size(eigenvalues)
         if (i /= mode) competitor = min(competitor, abs(eigenvalues(i)))
      end do
      if (competitor <= tiny(1.0_rp)) then
         value = huge(1.0_rp)
      else
         value = abs(eigenvalues(mode))/competitor
      end if
   end function isolation_ratio

   subroutine validate_mode_identification(result, accepted, refusal)
      type(lr_goldstone_correction_result), intent(in) :: result
      logical, intent(out) :: accepted
      character(len=*), intent(out) :: refusal

      accepted = .false.
      refusal = 'BLOCKED: Goldstone mode identification is ambiguous'
      if (result%goldstone_isolation_ratio > lr_gcr_isolation_ratio) then
         refusal = 'BLOCKED: smallest eigenvalue is not isolated from non-Goldstone modes'
         return
      end if
      if (result%rigid_overlap < lr_gcr_min_rigid_overlap) then
         refusal = 'BLOCKED: smallest eigenvector has weak LR-03 rigid-rotation overlap'
         return
      end if
      accepted = .true.
   end subroutine validate_mode_identification

   subroutine validate_consistency(request, accepted, refusal)
      type(lr_goldstone_correction_request), intent(in) :: request
      logical, intent(out) :: accepted
      character(len=*), intent(out) :: refusal
      integer :: i
      logical :: changed_eta, changed_kmesh, has_reference

      accepted = .false.
      refusal = 'BLOCKED: eta/k-mesh Goldstone consistency evidence is missing'
      if (.not. allocated(request%consistency_samples) .or. &
          size(request%consistency_samples) < lr_gcr_min_consistency_samples) return
      changed_eta = .false.
      changed_kmesh = .false.
      has_reference = .false.
      do i = 1, size(request%consistency_samples)
         if (.not. finite_real(request%consistency_samples(i)%eta) .or. &
             .not. finite_complex(request%consistency_samples(i)%goldstone_eigenvalue) .or. &
             .not. finite_real(request%consistency_samples(i)%competing_eigenvalue_abs) .or. &
             .not. finite_real(request%consistency_samples(i)%normalized_rigid_overlap) .or. &
             request%consistency_samples(i)%eta < 0.0_rp .or. &
             request%consistency_samples(i)%k_mesh_points < 1 .or. &
             request%consistency_samples(i)%competing_eigenvalue_abs <= 0.0_rp) then
            refusal = 'BLOCKED: an eta/k-mesh consistency sample is incomplete'
            return
         end if
         if (request%consistency_samples(i)%normalized_rigid_overlap < lr_gcr_min_rigid_overlap) then
            refusal = 'BLOCKED: rigid-mode character is not consistent across eta/k-mesh samples'
            return
         end if
         if (abs(request%consistency_samples(i)%goldstone_eigenvalue) > &
             lr_gcr_isolation_ratio*request%consistency_samples(i)%competing_eigenvalue_abs) then
            refusal = 'BLOCKED: Goldstone isolation is not consistent across eta/k-mesh samples'
            return
         end if
         if (abs(request%consistency_samples(i)%eta - request%eta) > lr_gcr_matrix_tolerance) changed_eta = .true.
         if (request%consistency_samples(i)%k_mesh_points /= request%k_mesh_points) changed_kmesh = .true.
         if (abs(request%consistency_samples(i)%eta - request%eta) <= lr_gcr_matrix_tolerance .and. &
             request%consistency_samples(i)%k_mesh_points == request%k_mesh_points) has_reference = .true.
      end do
      if (.not. has_reference) then
         refusal = 'BLOCKED: consistency samples do not include the selected static point'
         return
      end if
      if (.not. changed_eta .and. .not. changed_kmesh) then
         refusal = 'BLOCKED: consistency samples do not contain an eta or k-mesh change'
         return
      end if
      accepted = .true.
   end subroutine validate_consistency

   logical function finite_real(value) result(is_finite)
      real(rp), intent(in) :: value

      is_finite = value == value .and. abs(value) < huge(1.0_rp)
   end function finite_real

   logical function finite_complex(value) result(is_finite)
      complex(rp), intent(in) :: value

      is_finite = finite_real(real(value, rp)) .and. finite_real(aimag(value))
   end function finite_complex

   subroutine diagonalize_hermitian(input, eigenvalues, vectors, info)
      complex(rp), intent(in) :: input(:, :)
      complex(rp), intent(out) :: vectors(:, :)
      complex(rp), intent(out) :: eigenvalues(:)
      integer, intent(out) :: info
      complex(rp), allocatable :: work(:)
      complex(rp) :: query(1)
      real(rp), allocatable :: real_eigenvalues(:), rwork(:)
      integer :: n, lwork
      external :: zheev

      n = size(input, 1)
      vectors = input
      allocate(real_eigenvalues(n), rwork(max(1, 3*n - 2)))
      call zheev('V', 'U', n, vectors, n, real_eigenvalues, query, -1, rwork, info)
      if (info /= 0) return
      lwork = max(1, nint(real(query(1), rp)))
      allocate(work(lwork))
      call zheev('V', 'U', n, vectors, n, real_eigenvalues, work, lwork, rwork, info)
      eigenvalues = cmplx(real_eigenvalues, 0.0_rp, rp)
   end subroutine diagonalize_hermitian

   subroutine diagonalize_nonhermitian(input, eigenvalues, left, right, info)
      complex(rp), intent(in) :: input(:, :)
      complex(rp), intent(out) :: eigenvalues(:), left(:, :), right(:, :)
      integer, intent(out) :: info
      complex(rp), allocatable :: work(:), matrix_work(:, :)
      complex(rp) :: query(1)
      real(rp), allocatable :: rwork(:)
      integer :: n, lwork
      external :: zgeev

      n = size(input, 1)
      allocate(matrix_work(n, n), rwork(max(1, 2*n)))
      matrix_work = input
      call zgeev('V', 'V', n, matrix_work, n, eigenvalues, left, n, right, n, query, -1, rwork, info)
      if (info /= 0) return
      lwork = max(1, nint(real(query(1), rp)))
      allocate(work(lwork))
      matrix_work = input
      call zgeev('V', 'V', n, matrix_work, n, eigenvalues, left, n, right, n, work, lwork, rwork, info)
   end subroutine diagonalize_nonhermitian

   subroutine invert_eigenvector_matrix(input, inverse, success)
      complex(rp), intent(in) :: input(:, :)
      complex(rp), intent(out) :: inverse(:, :)
      logical, intent(out) :: success
      integer :: n, info, lwork
      integer, allocatable :: pivots(:)
      complex(rp), allocatable :: work(:)
      external :: zgetrf, zgetri

      n = size(input, 1)
      inverse = input
      allocate(pivots(n))
      call zgetrf(n, n, inverse, n, pivots, info)
      if (info /= 0) then
         success = .false.
         return
      end if
      lwork = max(1, n*n)
      allocate(work(lwork))
      call zgetri(n, inverse, n, pivots, work, lwork, info)
      success = info == 0
   end subroutine invert_eigenvector_matrix

   subroutine reconstruct_from_eigensystem(vectors, inverse, eigenvalues, matrix)
      complex(rp), intent(in) :: vectors(:, :), inverse(:, :), eigenvalues(:)
      complex(rp), intent(out) :: matrix(:, :)
      integer :: i, j

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do j = 1, size(eigenvalues)
         do i = 1, size(eigenvalues)
            matrix(i, :) = matrix(i, :) + vectors(i, j)*eigenvalues(j)*inverse(j, :)
         end do
      end do
   end subroutine reconstruct_from_eigensystem

   subroutine solve_matrix(matrix, rhs, success)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp), intent(inout) :: rhs(:, :)
      logical, intent(out) :: success
      complex(rp), allocatable :: work_matrix(:, :)
      integer, allocatable :: pivots(:)
      integer :: n, nrhs, info
      external :: zgesv

      n = size(matrix, 1)
      nrhs = size(rhs, 2)
      allocate(work_matrix(n, n), pivots(n))
      work_matrix = matrix
      call zgesv(n, nrhs, work_matrix, n, pivots, rhs, n, info)
      success = info == 0
   end subroutine solve_matrix

end module lr_goldstone_correction_mod
