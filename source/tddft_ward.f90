!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Ground-state Ward identity and explicit Goldstone repairs.
!>
!> The physical object checked here is
!>
!>   chi_KS(0,0) B_xc - m = 0,    B_xc = K_xc m.
!>
!> The routines accept the active response basis explicitly and retain the
!> uncorrected vector residual.  Sum-rule reconstruction and eigenvalue
!> projection are separate, opt-in operations; neither is used implicitly by
!> the diagnose path.
module tddft_ward_mod
   use precision_mod, only: rp
   implicit none

   private

   real(rp), parameter :: static_imaginary_tolerance = 1.0e-10_rp
   real(rp), parameter :: minimum_moment_relative = 1.0e-10_rp
   real(rp), parameter :: svd_rank_relative_tolerance = 1.0e-10_rp
   real(rp), parameter :: maximum_condition_number = 1.0e8_rp
   real(rp), parameter :: default_projection_overlap = 0.5_rp
   real(rp), parameter :: default_projection_distance = 0.5_rp

   type, public :: tddft_ward_diagnostics
      logical :: available = .false.
      logical :: identity_consistent = .false.
      character(len=32) :: response_basis = 'unrecorded'
      character(len=160) :: bxc_provenance = 'unrecorded'
      character(len=160) :: kernel_provenance = 'unrecorded'
      real(rp) :: magnetization_norm = -1.0_rp
      real(rp) :: bxc_norm = -1.0_rp
      real(rp) :: chi_bxc_norm = -1.0_rp
      real(rp) :: ward_residual = -1.0_rp
      ! Alias used by response reports that predate the explicit Ward name.
      real(rp) :: residual = -1.0_rp
      real(rp) :: dm_residual = -1.0_rp
      real(rp) :: bxc_kernel_residual = -1.0_rp
      real(rp) :: imaginary_norm = -1.0_rp
      complex(rp), allocatable :: magnetization(:)
      complex(rp), allocatable :: bxc(:)
      complex(rp), allocatable :: chi_bxc(:)
      complex(rp), allocatable :: ward_vector(:)
      complex(rp), allocatable :: d_m_vector(:)
   end type tddft_ward_diagnostics

   type, public :: tddft_lounis_repair
      logical :: requested = .false.
      logical :: applied = .false.
      logical :: rejected = .false.
      logical :: large_correction = .false.
      character(len=160) :: decision = 'not requested'
      integer :: effective_rank = 0
      real(rp) :: condition_number = -1.0_rp
      real(rp) :: maximum_kernel_change = -1.0_rp
      real(rp) :: relative_kernel_change = -1.0_rp
      complex(rp), allocatable :: kernel(:,:)
      complex(rp), allocatable :: bxc(:)
      type(tddft_ward_diagnostics) :: raw
      type(tddft_ward_diagnostics) :: corrected
   end type tddft_lounis_repair

   type, public :: tddft_goldstone_projection
      logical :: requested = .false.
      logical :: applied = .false.
      logical :: rejected = .false.
      logical :: large_correction = .false.
      character(len=160) :: decision = 'not requested'
      integer :: eigenvalue_index = 0
      complex(rp) :: eigenvalue_before = cmplx(0.0_rp, 0.0_rp, rp)
      complex(rp) :: eigenvalue_after = cmplx(0.0_rp, 0.0_rp, rp)
      real(rp) :: eigenvalue_distance = -1.0_rp
      real(rp) :: right_magnetization_overlap = -1.0_rp
      real(rp) :: left_magnetization_overlap = -1.0_rp
      real(rp) :: biorthogonal_overlap = -1.0_rp
      real(rp) :: relative_operator_change = -1.0_rp
      type(tddft_ward_diagnostics) :: raw
      type(tddft_ward_diagnostics) :: corrected
   end type tddft_goldstone_projection

   public :: evaluate_static_ward_identity
   public :: evaluate_ward_from_xi
   public :: reconstruct_lounis_kernel
   public :: project_goldstone_eigenvalue
   public :: derive_kernel_from_static_xi
   public :: write_ward_diagnostics_text

contains

   !> Evaluate chi_KS B_xc-m and, when K_xc is supplied, Dm=m-chi_KS K_xc m.
   !> `kernel` is the full operator in the active response basis, not an
   !> implicit site scalar.  The optional provenance strings are emitted with
   !> the residual so an apparently small number cannot lose its origin.
   subroutine evaluate_static_ward_identity(chi_ks, bxc, magnetization, diagnostics, kernel, response_basis, &
      bxc_provenance, kernel_provenance)
      complex(rp), intent(in) :: chi_ks(:, :), bxc(:), magnetization(:)
      type(tddft_ward_diagnostics), intent(out) :: diagnostics
      complex(rp), intent(in), optional :: kernel(:, :)
      character(len=*), intent(in), optional :: response_basis, bxc_provenance, kernel_provenance
      complex(rp), allocatable :: bxc_from_kernel(:)
      real(rp) :: denom

      call require_square_vector(chi_ks, bxc, magnetization, 'evaluate_static_ward_identity')
      diagnostics%available = .true.
      if (present(response_basis)) diagnostics%response_basis = trim(response_basis)
      if (present(bxc_provenance)) diagnostics%bxc_provenance = trim(bxc_provenance)
      if (present(kernel_provenance)) diagnostics%kernel_provenance = trim(kernel_provenance)
      allocate(diagnostics%magnetization(size(magnetization)), diagnostics%bxc(size(bxc)), &
         diagnostics%chi_bxc(size(magnetization)), diagnostics%ward_vector(size(magnetization)), &
         diagnostics%d_m_vector(size(magnetization)))
      diagnostics%magnetization = magnetization
      diagnostics%bxc = bxc
      diagnostics%magnetization_norm = vector_norm(magnetization)
      diagnostics%bxc_norm = vector_norm(bxc)
      diagnostics%imaginary_norm = matrix_imaginary_norm(chi_ks)
      if (diagnostics%magnetization_norm <= tiny(1.0_rp)) then
         error stop 'evaluate_static_ward_identity: magnetization is zero'
      end if
      diagnostics%chi_bxc = matmul(chi_ks, bxc)
      diagnostics%chi_bxc_norm = vector_norm(diagnostics%chi_bxc)
      diagnostics%ward_vector = diagnostics%chi_bxc - magnetization
      denom = max(diagnostics%magnetization_norm, tiny(1.0_rp))
      diagnostics%ward_residual = vector_norm(diagnostics%ward_vector)/denom
      diagnostics%residual = diagnostics%ward_residual
      diagnostics%d_m_vector = cmplx(0.0_rp, 0.0_rp, rp)
      if (present(kernel)) then
         if (size(kernel, 1) /= size(kernel, 2) .or. size(kernel, 1) /= size(magnetization)) then
            error stop 'evaluate_static_ward_identity: kernel and response dimensions are incompatible'
         end if
         allocate(bxc_from_kernel(size(magnetization)))
         bxc_from_kernel = matmul(kernel, magnetization)
         diagnostics%bxc_kernel_residual = vector_norm(bxc-bxc_from_kernel)/max(1.0_rp, diagnostics%bxc_norm)
         diagnostics%d_m_vector = magnetization - matmul(chi_ks, bxc_from_kernel)
         diagnostics%dm_residual = vector_norm(diagnostics%d_m_vector)/denom
      else
         diagnostics%dm_residual = -1.0_rp
         diagnostics%bxc_kernel_residual = -1.0_rp
      end if
      diagnostics%identity_consistent = diagnostics%ward_residual <= static_imaginary_tolerance
   end subroutine evaluate_static_ward_identity

   !> Evaluate the same residual when only Xi=chi_KS K_xc is available.
   !> This is useful for the direct LMTO pair-potential route; it remains an
   !> explicitly labelled Dm diagnostic and does not pretend that B_xc was
   !> reconstructed from an unrelated finite-eta inverse.
   subroutine evaluate_ward_from_xi(xi, magnetization, diagnostics, response_basis, kernel_provenance)
      complex(rp), intent(in) :: xi(:, :), magnetization(:)
      type(tddft_ward_diagnostics), intent(out) :: diagnostics
      character(len=*), intent(in), optional :: response_basis, kernel_provenance
      complex(rp), allocatable :: response_on_m(:)
      real(rp) :: denom

      if (size(xi, 1) /= size(xi, 2) .or. size(magnetization) /= size(xi, 1)) then
         error stop 'evaluate_ward_from_xi: Xi and magnetization dimensions are incompatible'
      end if
      if (vector_norm(magnetization) <= tiny(1.0_rp)) error stop 'evaluate_ward_from_xi: magnetization is zero'
      diagnostics%available = .true.
      diagnostics%bxc_provenance = 'not supplied; Xi-only diagnostic'
      diagnostics%kernel_provenance = 'Xi operator supplied directly'
      if (present(response_basis)) diagnostics%response_basis = trim(response_basis)
      if (present(kernel_provenance)) diagnostics%kernel_provenance = trim(kernel_provenance)
      allocate(diagnostics%magnetization(size(magnetization)), diagnostics%chi_bxc(size(magnetization)), &
         diagnostics%ward_vector(size(magnetization)), diagnostics%d_m_vector(size(magnetization)))
      diagnostics%magnetization = magnetization
      diagnostics%magnetization_norm = vector_norm(magnetization)
      allocate(response_on_m(size(magnetization)))
      response_on_m = matmul(xi, magnetization)
      diagnostics%chi_bxc = response_on_m
      diagnostics%chi_bxc_norm = vector_norm(response_on_m)
      diagnostics%ward_vector = response_on_m - magnetization
      diagnostics%d_m_vector = magnetization - response_on_m
      denom = max(diagnostics%magnetization_norm, tiny(1.0_rp))
      diagnostics%ward_residual = vector_norm(diagnostics%ward_vector)/denom
      diagnostics%residual = diagnostics%ward_residual
      diagnostics%dm_residual = vector_norm(diagnostics%d_m_vector)/denom
      diagnostics%identity_consistent = diagnostics%ward_residual <= static_imaginary_tolerance
      diagnostics%imaginary_norm = matrix_imaginary_norm(xi)
   end subroutine evaluate_ward_from_xi

   !> Reconstruct a local active-basis K_xc from the static Ward equation.
   !>
   !> This is the Lounis-style sum-rule operation for the site-diagonal
   !> response basis: solve Re[chi_KS] diag(m) k = m.  It is intentionally an
   !> explicit call.  A physical diagonal kernel may be supplied only for a
   !> before/after correction audit; it is never used as an empirical scale.
   subroutine reconstruct_lounis_kernel(chi_ks, magnetization, kernel, repair, physical_kernel, warning_threshold, &
      failure_threshold)
      complex(rp), intent(in) :: chi_ks(:, :), magnetization(:)
      complex(rp), allocatable, intent(out) :: kernel(:, :)
      type(tddft_lounis_repair), intent(out) :: repair
      complex(rp), intent(in), optional :: physical_kernel(:)
      real(rp), intent(in), optional :: warning_threshold, failure_threshold
      real(rp), allocatable :: system(:, :), rhs(:), solution(:)
      complex(rp), allocatable :: physical_bxc(:), physical_operator(:, :)
      real(rp) :: max_moment, relative_imaginary, denom
      integer :: n, i
      logical :: solved
      character(len=160) :: reason

      repair%requested = .true.
      n = size(magnetization)
      if (n < 1 .or. size(chi_ks, 1) /= n .or. size(chi_ks, 2) /= n) then
         call reject_lounis(repair, 'static chi_KS and magnetization dimensions are incompatible')
         allocate(kernel(max(0, n), max(0, n)))
         return
      end if
      allocate(kernel(n, n), system(n, n), rhs(n), solution(n))
      kernel = cmplx(0.0_rp, 0.0_rp, rp)
      if (present(physical_kernel)) then
         if (size(physical_kernel) /= n) then
            call reject_lounis(repair, 'physical diagonal kernel dimensions are incompatible')
            return
         end if
         allocate(physical_bxc(n), physical_operator(n, n))
         physical_operator = cmplx(0.0_rp, 0.0_rp, rp)
         do i = 1, n
            physical_bxc(i) = physical_kernel(i)*magnetization(i)
            physical_operator(i, i) = physical_kernel(i)
         end do
         call evaluate_static_ward_identity(chi_ks, physical_bxc, magnetization, repair%raw, kernel=physical_operator, &
            response_basis='site', &
            bxc_provenance='ground-state XC response provider', kernel_provenance='physical site-diagonal K_xc')
      end if
      if (maxval(abs(aimag(magnetization))) > minimum_moment_relative*max(1.0_rp, maxval(abs(magnetization)))) then
         call reject_lounis(repair, 'Lounis reconstruction requires a real collinear magnetization')
         return
      end if
      max_moment = maxval(abs(real(magnetization, rp)))
      if (max_moment <= tiny(1.0_rp) .or. any(abs(real(magnetization, rp)) <= minimum_moment_relative*max_moment)) then
         call reject_lounis(repair, 'one or more response moments are too small for local reconstruction')
         return
      end if
      relative_imaginary = matrix_imaginary_norm(chi_ks)
      if (relative_imaginary > static_imaginary_tolerance) then
         call reject_lounis(repair, 'static chi_KS has material imaginary content')
         return
      end if
      do i = 1, n
         system(:, i) = real(chi_ks(:, i), rp)*real(magnetization(i), rp)
      end do
      rhs = real(magnetization, rp)
      call solve_real_svd(system, rhs, solution, repair%effective_rank, repair%condition_number, solved, reason)
      if (.not. solved) then
         call reject_lounis(repair, trim(reason))
         return
      end if
      do i = 1, n
         kernel(i, i) = cmplx(solution(i), 0.0_rp, rp)
      end do
      repair%kernel = kernel
      allocate(repair%bxc(n))
      repair%bxc = matmul(kernel, magnetization)
      call evaluate_static_ward_identity(chi_ks, repair%bxc, magnetization, repair%corrected, kernel=kernel, &
         response_basis='site', bxc_provenance='explicit Lounis reconstructed K_xc', &
         kernel_provenance='Lounis static Ward reconstruction in active site basis')
      if (present(physical_kernel)) then
         denom = max(1.0_rp, maxval(abs(physical_kernel)), maxval(abs(solution)))
         repair%maximum_kernel_change = maxval(abs(solution-physical_kernel))
         repair%relative_kernel_change = repair%maximum_kernel_change/denom
         if (present(warning_threshold)) then
            repair%large_correction = warning_threshold >= 0.0_rp .and. repair%relative_kernel_change > warning_threshold
         end if
         if (present(failure_threshold)) then
            if (failure_threshold >= 0.0_rp .and. repair%relative_kernel_change > failure_threshold) then
               call reject_lounis(repair, 'Lounis reconstruction exceeds the explicitly supplied failure threshold')
               return
            end if
         end if
      end if
      repair%applied = .true.
      if (repair%large_correction) then
         repair%decision = 'accepted with warning: explicit Lounis reconstruction changed the physical kernel materially'
      else
         repair%decision = 'accepted: explicit Lounis static Ward reconstruction in active site basis'
      end if
   end subroutine reconstruct_lounis_kernel

   !> Project only the eigenvalue identified as the magnetization-like
   !> Goldstone mode to one.  The biorthogonal rank-one update leaves the other
   !> eigenmodes untouched to first order and records its full operator change.
   subroutine project_goldstone_eigenvalue(xi, magnetization, projected, repair, minimum_overlap, maximum_distance, &
      warning_threshold, failure_threshold)
      complex(rp), intent(in) :: xi(:, :), magnetization(:)
      complex(rp), allocatable, intent(out) :: projected(:, :)
      type(tddft_goldstone_projection), intent(out) :: repair
      real(rp), intent(in), optional :: minimum_overlap, maximum_distance, warning_threshold, failure_threshold
      complex(rp), allocatable :: eigenvalues(:), right(:, :), left(:, :), delta(:, :)
      complex(rp) :: overlap, coefficient
      real(rp) :: overlap_cutoff, distance_cutoff, norm_m, norm_right, norm_left, denom
      integer :: n, i, j, index

      repair%requested = .true.
      n = size(magnetization)
      if (n < 1 .or. size(xi, 1) /= n .or. size(xi, 2) /= n) then
         call reject_projection(repair, 'Xi and magnetization dimensions are incompatible')
         allocate(projected(max(0, n), max(0, n)))
         return
      end if
      allocate(projected(n, n))
      projected = xi
      call evaluate_ward_from_xi(xi, magnetization, repair%raw, response_basis='active response basis', &
         kernel_provenance='raw Xi before Halle Goldstone projection')
      call diagonalize_nonhermitian(xi, eigenvalues, right, left)
      index = 1
      do i = 2, n
         if (abs(eigenvalues(i)-cmplx(1.0_rp, 0.0_rp, rp)) < abs(eigenvalues(index)-cmplx(1.0_rp, 0.0_rp, rp))) index = i
      end do
      repair%eigenvalue_index = index
      repair%eigenvalue_before = eigenvalues(index)
      repair%eigenvalue_distance = abs(eigenvalues(index)-cmplx(1.0_rp, 0.0_rp, rp))
      norm_m = vector_norm(magnetization)
      norm_right = vector_norm(right(:, index))
      norm_left = vector_norm(left(:, index))
      repair%right_magnetization_overlap = abs(dot_product(magnetization, right(:, index)))/(norm_m*norm_right)
      repair%left_magnetization_overlap = abs(dot_product(left(:, index), magnetization))/(norm_left*norm_m)
      overlap = dot_product(left(:, index), right(:, index))
      repair%biorthogonal_overlap = abs(overlap)/(norm_left*norm_right)
      overlap_cutoff = default_projection_overlap
      distance_cutoff = default_projection_distance
      if (present(minimum_overlap)) overlap_cutoff = minimum_overlap
      if (present(maximum_distance)) distance_cutoff = maximum_distance
      if (repair%right_magnetization_overlap < overlap_cutoff) then
         call reject_projection(repair, 'closest-to-unity eigenmode is not magnetization-like')
         return
      end if
      if (repair%eigenvalue_distance > distance_cutoff) then
         call reject_projection(repair, 'no sufficiently small spurious Goldstone eigenvalue was found')
         return
      end if
      if (abs(overlap) <= tiny(1.0_rp)) then
         call reject_projection(repair, 'Goldstone eigenmode has zero biorthogonal normalization')
         return
      end if
      coefficient = (cmplx(1.0_rp, 0.0_rp, rp)-eigenvalues(index))/overlap
      do j = 1, n
         projected(:, j) = xi(:, j) + coefficient*right(:, index)*conjg(left(j, index))
      end do
      repair%eigenvalue_after = cmplx(1.0_rp, 0.0_rp, rp)
      allocate(delta(n, n))
      delta = projected-xi
      denom = max(1.0_rp, sqrt(sum(abs(xi)**2)))
      repair%relative_operator_change = sqrt(sum(abs(delta)**2))/denom
      if (present(warning_threshold)) repair%large_correction = warning_threshold >= 0.0_rp .and. &
         repair%relative_operator_change > warning_threshold
      if (present(failure_threshold)) then
         if (failure_threshold >= 0.0_rp .and. repair%relative_operator_change > failure_threshold) then
            call reject_projection(repair, 'Halle projection exceeds the explicitly supplied failure threshold')
            return
         end if
      end if
      call evaluate_ward_from_xi(projected, magnetization, repair%corrected, response_basis='active response basis', &
         kernel_provenance='explicit Halle rank-one Goldstone projection')
      repair%applied = .true.
      if (repair%large_correction) then
         repair%decision = 'accepted with warning: explicit Halle projection changed Xi materially'
      else
         repair%decision = 'accepted: explicit Halle biorthogonal rank-one Goldstone projection'
      end if
   end subroutine project_goldstone_eigenvalue

   !> Convert an explicitly projected static Xi back to the adiabatic kernel
   !> representation used by the Dyson solver.  This is a linear solve, not an
   !> inverse built from finite-eta response data.
   subroutine derive_kernel_from_static_xi(chi_ks, xi, kernel, info)
      complex(rp), intent(in) :: chi_ks(:, :), xi(:, :)
      complex(rp), allocatable, intent(out) :: kernel(:, :)
      integer, intent(out) :: info
      complex(rp), allocatable :: coefficient_matrix(:, :)
      integer, allocatable :: pivots(:)
      integer :: n

      n = size(chi_ks, 1)
      info = -1
      if (n < 1 .or. size(chi_ks, 2) /= n .or. size(xi, 1) /= n .or. size(xi, 2) /= n) then
         allocate(kernel(max(0, n), max(0, n)))
         return
      end if
      allocate(coefficient_matrix(n, n), kernel(n, n), pivots(n))
      coefficient_matrix = chi_ks
      kernel = xi
      call zgesv(n, n, coefficient_matrix, n, pivots, kernel, n, info)
   end subroutine derive_kernel_from_static_xi

   subroutine write_ward_diagnostics_text(filename, diagnostics)
      character(len=*), intent(in) :: filename
      type(tddft_ward_diagnostics), intent(in) :: diagnostics
      integer :: unit, ios

      open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
      if (ios /= 0) error stop 'write_ward_diagnostics_text: cannot open output file'
      write(unit, '(a)') '# Ward identity: chi_KS(0,0) B_xc - m; Dm = m-chi_KS K_xc m'
      write(unit, '(a,l1)') '# ward_available = ', diagnostics%available
      write(unit, '(a,l1)') '# ward_identity_consistent = ', diagnostics%identity_consistent
      write(unit, '(a,a)') '# response_basis = ', trim(diagnostics%response_basis)
      write(unit, '(a,a)') '# bxc_provenance = ', trim(diagnostics%bxc_provenance)
      write(unit, '(a,a)') '# kernel_provenance = ', trim(diagnostics%kernel_provenance)
      if (diagnostics%available) then
         write(unit, '(a,es24.16)') '# magnetization_norm = ', diagnostics%magnetization_norm
         write(unit, '(a,es24.16)') '# bxc_norm = ', diagnostics%bxc_norm
         write(unit, '(a,es24.16)') '# ward_residual = ', diagnostics%ward_residual
         write(unit, '(a,es24.16)') '# dm_residual = ', diagnostics%dm_residual
         write(unit, '(a,es24.16)') '# bxc_kernel_residual = ', diagnostics%bxc_kernel_residual
      end if
      close(unit)
   end subroutine write_ward_diagnostics_text

   subroutine solve_real_svd(matrix, rhs, solution, effective_rank, condition_number, solved, reason)
      real(rp), intent(in) :: matrix(:, :), rhs(:)
      real(rp), intent(out) :: solution(:)
      integer, intent(out) :: effective_rank
      real(rp), intent(out) :: condition_number
      logical, intent(out) :: solved
      character(len=*), intent(out) :: reason
      real(rp), allocatable :: factor(:, :), work(:), singular(:), u(:, :), vt(:, :), projected_rhs(:)
      real(rp) :: work_query(1), max_singular, rank_cutoff
      integer :: n, info, lwork

      solved = .false.
      reason = 'SVD solve failed'
      n = size(rhs)
      effective_rank = 0
      condition_number = -1.0_rp
      if (size(matrix, 1) /= n .or. size(matrix, 2) /= n .or. size(solution) /= n) then
         reason = 'SVD system dimensions are incompatible'
         return
      end if
      allocate(factor(n, n), singular(n), u(n, n), vt(n, n))
      factor = matrix
      call dgesvd('S', 'S', n, n, factor, n, singular, u, n, vt, n, work_query, -1, info)
      if (info /= 0) then
         reason = 'SVD workspace query failed'
         return
      end if
      lwork = max(1, int(work_query(1)))
      allocate(work(lwork))
      factor = matrix
      call dgesvd('S', 'S', n, n, factor, n, singular, u, n, vt, n, work, lwork, info)
      if (info /= 0) then
         reason = 'SVD factorization failed'
         return
      end if
      max_singular = maxval(singular)
      if (max_singular <= tiny(1.0_rp)) then
         reason = 'Ward reconstruction constraint is rank deficient'
         return
      end if
      rank_cutoff = svd_rank_relative_tolerance*max_singular
      effective_rank = count(singular > rank_cutoff)
      if (effective_rank /= n) then
         reason = 'Ward reconstruction constraint is rank deficient'
         return
      end if
      condition_number = max_singular/minval(singular)
      if (condition_number > maximum_condition_number) then
         reason = 'Ward reconstruction constraint is too ill-conditioned'
         return
      end if
      allocate(projected_rhs(n))
      projected_rhs = matmul(transpose(u), rhs)/singular
      solution = matmul(transpose(vt), projected_rhs)
      if (any(.not. finite_real(solution))) then
         reason = 'Ward reconstruction produced a nonfinite kernel'
         return
      end if
      solved = .true.
      reason = 'solved'
   end subroutine solve_real_svd

   subroutine diagonalize_nonhermitian(matrix, eigenvalues, right, left)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp), allocatable, intent(out) :: eigenvalues(:), right(:, :), left(:, :)
      complex(rp), allocatable :: work_matrix(:, :), work(:)
      real(rp), allocatable :: rwork(:)
      complex(rp) :: work_query(1)
      integer :: n, info, lwork

      n = size(matrix, 1)
      allocate(work_matrix(n, n), eigenvalues(n), right(n, n), left(n, n), rwork(max(1, 2*n)))
      work_matrix = matrix
      call zgeev('V', 'V', n, work_matrix, n, eigenvalues, left, n, right, n, work_query, -1, rwork, info)
      if (info /= 0) error stop 'project_goldstone_eigenvalue: LAPACK workspace query failed'
      lwork = max(1, int(real(work_query(1), rp)))
      allocate(work(lwork))
      work_matrix = matrix
      call zgeev('V', 'V', n, work_matrix, n, eigenvalues, left, n, right, n, work, lwork, rwork, info)
      if (info /= 0) error stop 'project_goldstone_eigenvalue: LAPACK zgeev failed'
   end subroutine diagonalize_nonhermitian

   subroutine require_square_vector(matrix, vector, magnetization, caller)
      complex(rp), intent(in) :: matrix(:, :), vector(:), magnetization(:)
      character(len=*), intent(in) :: caller
      if (size(matrix, 1) /= size(matrix, 2) .or. size(vector) /= size(matrix, 1) .or. &
          size(magnetization) /= size(matrix, 1)) error stop trim(caller)//': response dimensions are incompatible'
   end subroutine require_square_vector

   real(rp) function vector_norm(vector) result(norm)
      complex(rp), intent(in) :: vector(:)
      norm = sqrt(sum(abs(vector)**2))
   end function vector_norm

   real(rp) function matrix_imaginary_norm(matrix) result(norm)
      complex(rp), intent(in) :: matrix(:, :)
      norm = sqrt(sum(aimag(matrix)**2))/max(1.0_rp, sqrt(sum(abs(matrix)**2)))
   end function matrix_imaginary_norm

   elemental logical function finite_real(value)
      real(rp), intent(in) :: value
      finite_real = abs(value) < huge(value)
   end function finite_real

   subroutine reject_lounis(repair, decision)
      type(tddft_lounis_repair), intent(inout) :: repair
      character(len=*), intent(in) :: decision
      repair%rejected = .true.
      repair%applied = .false.
      repair%decision = trim(decision)
   end subroutine reject_lounis

   subroutine reject_projection(repair, decision)
      type(tddft_goldstone_projection), intent(inout) :: repair
      character(len=*), intent(in) :: decision
      repair%rejected = .true.
      repair%applied = .false.
      repair%decision = trim(decision)
   end subroutine reject_projection

end module tddft_ward_mod
