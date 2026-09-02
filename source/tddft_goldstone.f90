!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Site-projected transverse ALSDA kernel and Goldstone diagnostics.
!>
!> `K_perp` is deliberately obtained only from xc_response_kernel_mod.  In
!> particular, this module never inspects hamiltonian%hxc, cx1, or another
!> assembled Hamiltonian quantity: those are not a documented TDDFT kernel.
!>
!> Goldstone correction is deliberately restricted to the Ward-consistent
!> pair-potential route.  It rescales its right-hand columns from the real
!> static Xi; it never injects a finite-eta inverse chi_KS as a kernel.
module tddft_goldstone_mod
   use precision_mod, only: rp
   use xc_response_kernel_mod, only: xc_response_kernel_provider, circular_transverse_kernel
   use tddft_ward_mod, only: tddft_ward_diagnostics, tddft_lounis_repair, tddft_goldstone_projection, &
      evaluate_static_ward_identity, evaluate_ward_from_xi, reconstruct_lounis_kernel, project_goldstone_eigenvalue, &
      derive_kernel_from_static_xi
   implicit none

   private

   integer, parameter, public :: GOLDSTONE_OFF = 0
   integer, parameter, public :: GOLDSTONE_DIAGNOSE = 1
   integer, parameter, public :: GOLDSTONE_CORRECT = 2
   ! Compatibility name for callers compiled against the pre-WR-05 module.
   integer, parameter, public :: GOLDSTONE_SUM_RULE = GOLDSTONE_CORRECT
   real(rp), parameter :: static_imaginary_tolerance = 1.0e-10_rp
   real(rp), parameter :: minimum_moment_relative = 1.0e-10_rp
   real(rp), parameter :: svd_rank_relative_tolerance = 1.0e-10_rp
   real(rp), parameter :: maximum_condition_number = 1.0e8_rp

   type, public :: tddft_goldstone_options
      ! Development default: report the raw Goldstone quality, never repair it.
      character(len=16) :: goldstone_mode = 'diagnose'
      ! Explicit TDDFT-02 cleanup policy.  The default only diagnoses; the
      ! legacy goldstone_mode remains for compatibility with existing inputs.
      character(len=16) :: goldstone_policy = 'diagnose'
      logical :: has_soc = .false.
      logical :: has_external_field = .false.
      character(len=16) :: circular_channel = 'plus_minus'
   end type tddft_goldstone_options

   type, public :: tddft_goldstone_diagnostics
      logical :: available = .false.
      logical :: has_bare_spectral_gap = .false.
      real(rp) :: bare_spectral_gap = -1.0_rp
      real(rp) :: residual = -1.0_rp
      real(rp) :: closest_eigenvalue_distance = -1.0_rp
      ! `magnetization_overlap` is retained as the right overlap for
      ! compatibility with existing campaign readers.
      real(rp) :: magnetization_overlap = -1.0_rp
      real(rp) :: left_magnetization_overlap = -1.0_rp
      real(rp) :: biorthogonal_magnetization_overlap = -1.0_rp
      real(rp) :: imaginary_norm = -1.0_rp
      integer :: closest_eigenvalue_index = 0
      complex(rp) :: closest_eigenvalue = cmplx(0.0_rp, 0.0_rp, rp)
      complex(rp), allocatable :: eigenvalues(:)
      complex(rp), allocatable :: eigenvectors(:, :)
      complex(rp), allocatable :: left_eigenvectors(:, :)
      complex(rp), allocatable :: closest_eigenvector(:)
      complex(rp), allocatable :: closest_left_eigenvector(:)
      ! First-class raw Ward record.  The legacy fields above remain for
      ! campaign compatibility; this component carries basis and provenance.
      type(tddft_ward_diagnostics) :: ward
   end type tddft_goldstone_diagnostics

   type, public :: tddft_goldstone_result
      complex(rp), allocatable :: k_perp(:)
      complex(rp), allocatable :: xi_raw(:, :)
      complex(rp), allocatable :: k_perp_sum_rule(:)
      complex(rp), allocatable :: xi_corrected(:, :)
      type(tddft_goldstone_diagnostics) :: raw
      type(tddft_goldstone_diagnostics) :: corrected
      logical :: sum_rule_requested = .false.
      logical :: sum_rule_applied = .false.
      logical :: sum_rule_disabled_by_symmetry_breaking = .false.
      character(len=16) :: goldstone_policy = 'diagnose'
      character(len=16) :: circular_channel = 'unspecified'
      type(tddft_lounis_repair) :: lounis
      type(tddft_goldstone_projection) :: projection
      complex(rp), allocatable :: kernel_corrected(:, :)
   end type tddft_goldstone_result

   type, public :: tddft_goldstone_column_correction
      logical :: requested = .false.
      logical :: applied = .false.
      logical :: rejected = .false.
      character(len=160) :: decision = 'not requested'
      real(rp), allocatable :: scales(:)
      integer :: effective_rank = 0
      real(rp) :: condition_number = -1.0_rp
      real(rp) :: maximum_change = -1.0_rp
      real(rp) :: relative_kernel_change = -1.0_rp
      real(rp) :: warning_threshold = -1.0_rp
      real(rp) :: failure_threshold = -1.0_rp
      logical :: large_correction = .false.
      type(tddft_goldstone_diagnostics) :: raw
      type(tddft_goldstone_diagnostics) :: corrected
   end type tddft_goldstone_column_correction

   public :: build_site_projected_k_perp
   public :: construct_transverse_xi
   public :: evaluate_goldstone
   public :: evaluate_raw_xi_diagnostics
   public :: build_goldstone_column_correction
   public :: rescale_xi_columns
   public :: rescale_pair_potential_columns
   public :: spectral_weights_are_nonnegative
   public :: spectral_weight_correction_is_acceptable
   public :: append_goldstone_column_correction_text
   public :: write_goldstone_diagnostics_text

contains

   !> Return the site-projected transverse ALSDA kernel recorded by the XC
   !> response provider.  `has_k_perp_circular` makes missing functional provenance a
   !> hard error rather than silently substituting a Hamiltonian-derived value.
   subroutine build_site_projected_k_perp(provider, k_perp)
      type(xc_response_kernel_provider), intent(in) :: provider
      complex(rp), allocatable, intent(out) :: k_perp(:)
      integer :: isite

      if (.not. allocated(provider%site)) then
         error stop 'build_site_projected_k_perp: XC response provider is not initialized'
      end if
      allocate(k_perp(size(provider%site)))
      do isite = 1, size(provider%site)
         if (.not. provider%site(isite)%has_k_perp_circular) then
            error stop 'build_site_projected_k_perp: site has no circular xc_response_kernel K_perp'
         end if
         k_perp(isite) = cmplx(circular_transverse_kernel(provider, isite), 0.0_rp, rp)
      end do
   end subroutine build_site_projected_k_perp

   !> Form Xi = chi_KS K_perp for a local site-diagonal transverse kernel.
   function construct_transverse_xi(chi_ks, k_perp) result(xi)
      complex(rp), intent(in) :: chi_ks(:, :), k_perp(:)
      complex(rp) :: xi(size(chi_ks, 1), size(chi_ks, 2))
      integer :: j

      call require_square_response(chi_ks, k_perp, 'construct_transverse_xi')
      do j = 1, size(k_perp)
         xi(:, j) = chi_ks(:, j)*k_perp(j)
      end do
   end function construct_transverse_xi

   !> Diagnose Xi(0,0), and optionally construct an explicitly marked static
   !> sum-rule correction.  `bare_spectral_gap`, when present, must be the
   !> caller’s independently determined lowest bare spin-flip transition; it
   !> is not inferred from one complex chi_KS matrix sample.
   subroutine evaluate_goldstone(chi_ks_static, provider, options, result, bare_spectral_gap)
      complex(rp), intent(in) :: chi_ks_static(:, :)
      type(xc_response_kernel_provider), intent(in) :: provider
      type(tddft_goldstone_options), intent(in) :: options
      type(tddft_goldstone_result), intent(out) :: result
      real(rp), intent(in), optional :: bare_spectral_gap
      complex(rp), allocatable :: magnetization(:), bxc(:), kernel(:, :), repaired_kernel(:, :), projected_xi(:, :)
      integer :: mode, isite
      character(len=16) :: policy

      call build_site_projected_k_perp(provider, result%k_perp)
      call require_square_response(chi_ks_static, result%k_perp, 'evaluate_goldstone')
      allocate(magnetization(size(result%k_perp)))
      magnetization = cmplx(signed_site_populations(provider), 0.0_rp, rp)
      if (sqrt(sum(abs(magnetization)**2)) <= tiny(1.0_rp)) then
         error stop 'evaluate_goldstone: projected ground-state magnetization is zero'
      end if

      result%xi_raw = construct_transverse_xi(chi_ks_static, result%k_perp)
      policy = trim(adjustl(options%goldstone_policy))
      result%goldstone_policy = policy
      result%circular_channel = trim(options%circular_channel)
      mode = goldstone_mode_code(options%goldstone_mode)
      result%sum_rule_requested = mode == GOLDSTONE_SUM_RULE
      if (mode /= GOLDSTONE_OFF) then
         call calculate_diagnostics(result%xi_raw, magnetization, result%raw, bare_spectral_gap)
         allocate(bxc(size(magnetization)), kernel(size(magnetization), size(magnetization)))
         kernel = cmplx(0.0_rp, 0.0_rp, rp)
         do isite = 1, size(magnetization)
            bxc(isite) = result%k_perp(isite)*magnetization(isite)
            kernel(isite, isite) = result%k_perp(isite)
         end do
         call evaluate_static_ward_identity(chi_ks_static, bxc, magnetization, result%raw%ward, kernel=kernel, &
            response_basis='site', bxc_provenance='ground-state XC response provider', &
            kernel_provenance='site-projected transverse ALSDA K_xc')
      end if

      select case (policy)
      case ('diagnose')
         continue
      case ('sum_rule')
         result%sum_rule_requested = .true.
         if (options%has_soc .or. options%has_external_field) then
            result%lounis%requested = .true.
            result%lounis%rejected = .true.
            result%lounis%decision = 'Lounis reconstruction is unavailable with SOC or an external symmetry-breaking field'
            return
         end if
         call reconstruct_lounis_kernel(chi_ks_static, magnetization, repaired_kernel, result%lounis, &
            physical_kernel=result%k_perp)
         if (result%lounis%applied) then
            allocate(result%k_perp_sum_rule(size(result%k_perp)))
            do isite = 1, size(result%k_perp)
               result%k_perp_sum_rule(isite) = repaired_kernel(isite, isite)
            end do
            result%xi_corrected = matmul(chi_ks_static, repaired_kernel)
            call calculate_diagnostics(result%xi_corrected, magnetization, result%corrected)
            result%sum_rule_applied = .true.
         end if
         return
      case ('projected')
         if (options%has_soc .or. options%has_external_field) then
            result%projection%requested = .true.
            result%projection%rejected = .true.
            result%projection%decision = 'Halle projection is unavailable with SOC or an external symmetry-breaking field'
            return
         end if
         call project_goldstone_eigenvalue(result%xi_raw, magnetization, projected_xi, result%projection)
         if (result%projection%applied) then
            result%xi_corrected = projected_xi
            call calculate_diagnostics(result%xi_corrected, magnetization, result%corrected)
            call derive_kernel_from_static_xi(chi_ks_static, projected_xi, result%kernel_corrected, isite)
            if (isite /= 0) then
               result%projection%applied = .false.
               result%projection%rejected = .true.
               result%projection%decision = 'Halle projection could not be represented in the adiabatic kernel basis'
            end if
         end if
         return
      case default
         error stop 'evaluate_goldstone: goldstone_policy must be diagnose, sum_rule, or projected'
      end select

      ! The legacy site-scalar route is now diagnose-only.  The production
      ! `correct` path invokes build_goldstone_column_correction on direct
      ! pair-potential Xi below; retaining this old inverse here would revive
      ! exactly the misleading correction WR-05 replaces.
      if (mode /= GOLDSTONE_CORRECT) return
      if (options%has_soc .or. options%has_external_field) then
         ! A zero-frequency transverse mode is not required if symmetry is
         ! broken physically.  Preserve raw diagnostics and disable repair.
         result%sum_rule_disabled_by_symmetry_breaking = .true.
         return
      end if
      if (any(abs(magnetization) <= tiny(1.0_rp))) then
         error stop 'evaluate_goldstone: legacy diagnostic requires nonzero magnetization on every response site'
      end if

      result%sum_rule_disabled_by_symmetry_breaking = .true.
   end subroutine evaluate_goldstone

   !> Return the signed z-population used by the transverse Goldstone vector.
   !> The explicit field is preferred; the fallback preserves compatibility
   !> with providers created before signed-spin provenance was added.
   function signed_site_populations(provider) result(populations)
      type(xc_response_kernel_provider), intent(in) :: provider
      real(rp) :: populations(size(provider%site))
      integer :: isite

      if (.not. allocated(provider%site)) error stop 'signed_site_populations: XC response provider is not initialized'
      do isite = 1, size(provider%site)
         if (provider%site(isite)%has_signed_spin_population) then
            populations(isite) = provider%site(isite)%signed_spin_population
         else
            populations(isite) = provider%site(isite)%spin_population
         end if
      end do
   end function signed_site_populations

   !> Diagnose an already assembled self-enhancement operator.  This keeps the
   !> pair-potential route out of the legacy site-scalar K interface while
   !> retaining exactly the same raw residual/eigenmode definitions.
   subroutine evaluate_raw_xi_diagnostics(xi, magnetization, diagnostics, bare_spectral_gap, response_basis, &
      kernel_provenance)
      complex(rp), intent(in) :: xi(:, :), magnetization(:)
      type(tddft_goldstone_diagnostics), intent(out) :: diagnostics
      real(rp), intent(in), optional :: bare_spectral_gap
      character(len=*), intent(in), optional :: response_basis, kernel_provenance

      if (size(xi, 1) /= size(xi, 2) .or. size(magnetization) /= size(xi, 1)) then
         error stop 'evaluate_raw_xi_diagnostics: Xi and magnetization dimensions are incompatible'
      end if
      if (sqrt(sum(abs(magnetization)**2)) <= tiny(1.0_rp)) then
         error stop 'evaluate_raw_xi_diagnostics: magnetization is zero'
      end if
      call calculate_diagnostics(xi, magnetization, diagnostics, bare_spectral_gap)
      if (present(response_basis)) diagnostics%ward%response_basis = trim(response_basis)
      if (present(kernel_provenance)) diagnostics%ward%kernel_provenance = trim(kernel_provenance)
   end subroutine evaluate_raw_xi_diagnostics

   !> Compute real column scales for direct pair-potential Xi.  With W=I the
   !> constrained minimum-change problem is
   !> min ||s-1|| subject to Re[Xi diag(M)] s=M.  A rank-revealing SVD is used
   !> even though the accepted square full-rank case has a unique solution.
   subroutine build_goldstone_column_correction(xi_static, magnetization, correction, maximum_allowed_change, warning_threshold)
      complex(rp), intent(in) :: xi_static(:, :), magnetization(:)
      type(tddft_goldstone_column_correction), intent(out) :: correction
      real(rp), intent(in), optional :: maximum_allowed_change
      real(rp), intent(in), optional :: warning_threshold
      real(rp), allocatable :: system(:, :), singular(:), u(:, :), vt(:, :), work(:), projected_rhs(:), scales(:)
      complex(rp), allocatable :: corrected_xi(:, :)
      real(rp) :: work_query(1), max_moment, max_singular, rank_cutoff
      integer :: n, i, info, lwork

      correction%requested = .true.
      if (present(maximum_allowed_change)) correction%failure_threshold = maximum_allowed_change
      if (present(warning_threshold)) correction%warning_threshold = warning_threshold
      n = size(magnetization)
      if (n < 1 .or. size(xi_static, 1) /= n .or. size(xi_static, 2) /= n) then
         call reject_correction(correction, 'static Xi and magnetization dimensions are incompatible')
         return
      end if
      call evaluate_raw_xi_diagnostics(xi_static, magnetization, correction%raw)
      if (correction%raw%imaginary_norm > static_imaginary_tolerance) then
         call reject_correction(correction, 'static Xi has material imaginary content')
         return
      end if
      if (maxval(abs(aimag(magnetization))) > minimum_moment_relative*max(1.0_rp, maxval(abs(magnetization)))) then
         call reject_correction(correction, 'Goldstone correction requires a real collinear magnetization')
         return
      end if
      max_moment = maxval(abs(real(magnetization, rp)))
      if (max_moment <= tiny(1.0_rp) .or. any(abs(real(magnetization, rp)) <= minimum_moment_relative*max_moment)) then
         call reject_correction(correction, 'one or more response moments are too small for column rescaling')
         return
      end if
      allocate(system(n, n), singular(n), u(n, n), vt(n, n))
      do i = 1, n
         system(:, i) = real(xi_static(:, i), rp)*real(magnetization(i), rp)
      end do
      call dgesvd('S', 'S', n, n, system, n, singular, u, n, vt, n, work_query, -1, info)
      if (info /= 0) then
         call reject_correction(correction, 'static correction SVD workspace query failed')
         return
      end if
      lwork = max(1, int(work_query(1)))
      allocate(work(lwork))
      do i = 1, n
         system(:, i) = real(xi_static(:, i), rp)*real(magnetization(i), rp)
      end do
      call dgesvd('S', 'S', n, n, system, n, singular, u, n, vt, n, work, lwork, info)
      if (info /= 0) then
         call reject_correction(correction, 'static correction SVD failed')
         return
      end if
      max_singular = maxval(singular)
      if (max_singular <= tiny(1.0_rp)) then
         call reject_correction(correction, 'static correction constraint is rank deficient')
         return
      end if
      rank_cutoff = svd_rank_relative_tolerance*max_singular
      correction%effective_rank = count(singular > rank_cutoff)
      if (correction%effective_rank /= n) then
         call reject_correction(correction, 'static correction constraint is rank deficient')
         return
      end if
      correction%condition_number = max_singular/minval(singular)
      if (correction%condition_number > maximum_condition_number) then
         call reject_correction(correction, 'static correction constraint is too ill-conditioned')
         return
      end if
      allocate(projected_rhs(n), scales(n))
      projected_rhs = matmul(transpose(u), real(magnetization, rp))/singular
      scales = matmul(transpose(vt), projected_rhs)
      if (any(.not. finite_real(scales))) then
         call reject_correction(correction, 'static correction produced a nonfinite scale')
         return
      end if
      correction%maximum_change = maxval(abs(scales-1.0_rp))
      if (present(maximum_allowed_change)) then
         if (maximum_allowed_change >= 0.0_rp .and. correction%maximum_change > maximum_allowed_change) then
            call reject_correction(correction, 'requested column rescaling exceeds the explicitly supplied safety limit')
            return
         end if
      end if
      correction%relative_kernel_change = correction%maximum_change
      if (present(warning_threshold)) correction%large_correction = warning_threshold >= 0.0_rp .and. &
         correction%maximum_change > warning_threshold
      if (present(maximum_allowed_change)) correction%large_correction = correction%large_correction .or. &
         (maximum_allowed_change >= 0.0_rp .and. correction%maximum_change > maximum_allowed_change)
      allocate(correction%scales(n), corrected_xi(n, n))
      correction%scales = scales
      call rescale_xi_columns(xi_static, correction%scales, corrected_xi)
      call evaluate_raw_xi_diagnostics(corrected_xi, magnetization, correction%corrected)
      correction%applied = .true.
      correction%decision = 'accepted: real static SVD column rescaling with W=I'
   end subroutine build_goldstone_column_correction

   subroutine rescale_xi_columns(xi, scales, corrected_xi)
      complex(rp), intent(in) :: xi(:, :)
      real(rp), intent(in) :: scales(:)
      complex(rp), intent(out) :: corrected_xi(:, :)
      integer :: j
      if (size(xi, 1) /= size(xi, 2) .or. size(scales) /= size(xi, 2) .or. any(shape(corrected_xi) /= shape(xi))) then
         error stop 'rescale_xi_columns: incompatible Xi/scale dimensions'
      end if
      do j = 1, size(scales)
         corrected_xi(:, j) = xi(:, j)*scales(j)
      end do
   end subroutine rescale_xi_columns

   subroutine rescale_pair_potential_columns(operators, scales)
      complex(rp), intent(inout) :: operators(:, :, :, :)
      real(rp), intent(in) :: scales(:)
      integer :: isite
      if (size(operators, 3) /= size(scales)) error stop 'rescale_pair_potential_columns: response-site mismatch'
      do isite = 1, size(scales)
         operators(:, :, isite, :) = scales(isite)*operators(:, :, isite, :)
      end do
   end subroutine rescale_pair_potential_columns

   logical function spectral_weights_are_nonnegative(weights, tolerance) result(acceptable)
      real(rp), intent(in) :: weights(:)
      real(rp), intent(in), optional :: tolerance
      real(rp) :: allowed_negative
      allowed_negative = 1.0e-10_rp
      if (present(tolerance)) allowed_negative = tolerance
      acceptable = all(weights >= -allowed_negative)
   end function spectral_weights_are_nonnegative

   !> A finite-eta circular retarded response may contain a negative
   !> low-frequency tail from the opposite-frequency Lehmann poles.  The
   !> Goldstone-column correction must therefore be judged against the raw
   !> pair response, not against an unconditional pointwise positivity test.
   !> It is rejected only when it creates negative weight or makes pre-existing
   !> negative weight more negative by more than a combined absolute/relative
   !> numerical tolerance.  The relative term is essential near a narrow
   !> finite-eta pole: a column-scale change at roundoff can otherwise produce
   !> an absolute loss difference above 1e-10 while remaining numerically
   !> identical to the raw spectrum.
   logical function spectral_weight_correction_is_acceptable(raw_weights, corrected_weights, tolerance) result(acceptable)
      real(rp), intent(in) :: raw_weights(:), corrected_weights(:)
      real(rp), intent(in), optional :: tolerance
      real(rp) :: allowed_change

      if (size(raw_weights) /= size(corrected_weights)) then
         error stop 'spectral_weight_correction_is_acceptable: raw/corrected shape mismatch'
      end if
      allowed_change = 1.0e-10_rp
      if (present(tolerance)) allowed_change = tolerance
      acceptable = all(.not. ((corrected_weights < -allowed_change) .and. &
         (corrected_weights < raw_weights - allowed_change*max(1.0_rp, abs(raw_weights), abs(corrected_weights)))))
   end function spectral_weight_correction_is_acceptable

   subroutine append_goldstone_column_correction_text(filename, correction)
      character(len=*), intent(in) :: filename
      type(tddft_goldstone_column_correction), intent(in) :: correction
      integer :: unit, ios, i

      open(newunit=unit, file=filename, status='old', position='append', action='write', iostat=ios)
      if (ios /= 0) error stop 'append_goldstone_column_correction_text: cannot append diagnostic output'
      write(unit, '(a)') '# pair_potential_goldstone_correction_begin'
      write(unit, '(a,l1)') '# goldstone_correction_requested = ', correction%requested
      write(unit, '(a,l1)') '# goldstone_correction_applied = ', correction%applied
      write(unit, '(a,l1)') '# goldstone_correction_rejected = ', correction%rejected
      write(unit, '(a,a)') '# goldstone_correction_decision = ', trim(correction%decision)
      write(unit, '(a,i0)') '# goldstone_correction_effective_rank = ', correction%effective_rank
      write(unit, '(a,es24.16)') '# goldstone_correction_condition_number = ', correction%condition_number
      write(unit, '(a,es24.16)') '# goldstone_correction_maximum_scale_change = ', correction%maximum_change
      write(unit, '(a,es24.16)') '# goldstone_correction_relative_kernel_change = ', correction%relative_kernel_change
      write(unit, '(a,es24.16)') '# goldstone_correction_warning_threshold = ', correction%warning_threshold
      write(unit, '(a,es24.16)') '# goldstone_correction_failure_threshold = ', correction%failure_threshold
      write(unit, '(a,l1)') '# goldstone_correction_large_change = ', correction%large_correction
      if (correction%raw%available) call write_one_diagnostics(unit, 'pair_correction_raw', correction%raw)
      if (correction%corrected%available) call write_one_diagnostics(unit, 'pair_correction_corrected', correction%corrected)
      if (allocated(correction%scales)) then
         do i = 1, size(correction%scales)
            write(unit, '(a,1x,i0,1x,es24.16)') 'pair_potential_column_scale', i, correction%scales(i)
         end do
      end if
      write(unit, '(a)') '# pair_potential_goldstone_correction_end'
      close(unit)
   end subroutine append_goldstone_column_correction_text

   !> Write raw diagnostics first.  Corrected diagnostics are an additional
   !> record, never a replacement, so output remains useful for convergence
   !> studies even when sum-rule mode was requested.
   subroutine write_goldstone_diagnostics_text(filename, result)
      character(len=*), intent(in) :: filename
      type(tddft_goldstone_result), intent(in) :: result
      integer :: unit, ios, i

      open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
      if (ios /= 0) error stop 'write_goldstone_diagnostics_text: cannot open output file'
      write(unit, '(a)') '# legacy Xi = chi_KS K_perp; K_perp provenance = xc_response_kernel'
      write(unit, '(a)') '# Ward identity = chi_KS(0,0) B_xc - m; Dm = m-chi_KS K_xc m'
      write(unit, '(a)') '# raw Goldstone diagnostics are retained for every non-off mode'
      write(unit, '(a,a)') '# goldstone_policy = ', trim(result%goldstone_policy)
      write(unit, '(a,a)') '# circular_channel = ', trim(result%circular_channel)
      write(unit, '(a,l1)') '# legacy_site_scalar_correction_requested = ', result%sum_rule_requested
      write(unit, '(a,l1)') '# legacy_site_scalar_correction_applied = ', result%sum_rule_applied
      write(unit, '(a,l1)') '# legacy_site_scalar_correction_disabled = ', result%sum_rule_disabled_by_symmetry_breaking
      if (result%lounis%requested) then
         write(unit, '(a,l1)') '# lounis_reconstruction_applied = ', result%lounis%applied
         write(unit, '(a,l1)') '# lounis_reconstruction_rejected = ', result%lounis%rejected
         write(unit, '(a,a)') '# lounis_reconstruction_decision = ', trim(result%lounis%decision)
         write(unit, '(a,es24.16)') '# lounis_relative_kernel_change = ', result%lounis%relative_kernel_change
         if (result%lounis%raw%available) call write_ward_record(unit, 'lounis_raw', result%lounis%raw)
         if (result%lounis%corrected%available) call write_ward_record(unit, 'lounis_corrected', result%lounis%corrected)
      end if
      if (result%projection%requested) then
         write(unit, '(a,l1)') '# halle_projection_applied = ', result%projection%applied
         write(unit, '(a,l1)') '# halle_projection_rejected = ', result%projection%rejected
         write(unit, '(a,a)') '# halle_projection_decision = ', trim(result%projection%decision)
         write(unit, '(a,2(1x,es24.16))') '# halle_projection_eigenvalue_before = ', &
            real(result%projection%eigenvalue_before, rp), aimag(result%projection%eigenvalue_before)
         write(unit, '(a,es24.16)') '# halle_projection_relative_operator_change = ', &
            result%projection%relative_operator_change
         if (result%projection%raw%available) call write_ward_record(unit, 'halle_raw', result%projection%raw)
         if (result%projection%corrected%available) call write_ward_record(unit, 'halle_corrected', result%projection%corrected)
      end if
      call write_one_diagnostics(unit, 'raw', result%raw)
      if (result%sum_rule_applied) call write_one_diagnostics(unit, 'legacy_site_scalar_corrected', result%corrected)
      if (allocated(result%k_perp)) then
         do i = 1, size(result%k_perp)
            write(unit, '(a,1x,i0,2(1x,es24.16))') 'kernel_raw', i, real(result%k_perp(i), rp), aimag(result%k_perp(i))
         end do
      end if
      if (allocated(result%k_perp_sum_rule)) then
         do i = 1, size(result%k_perp_sum_rule)
            write(unit, '(a,1x,i0,2(1x,es24.16))') 'legacy_site_scalar_kernel_correction', i, real(result%k_perp_sum_rule(i), rp), &
               aimag(result%k_perp_sum_rule(i))
         end do
      end if
      close(unit)
   end subroutine write_goldstone_diagnostics_text

   subroutine calculate_diagnostics(xi, magnetization, diagnostics, bare_spectral_gap)
      complex(rp), intent(in) :: xi(:, :), magnetization(:)
      type(tddft_goldstone_diagnostics), intent(out) :: diagnostics
      real(rp), intent(in), optional :: bare_spectral_gap
      complex(rp), allocatable :: residual_vector(:)
      real(rp) :: norm_m, norm_right, norm_left, matrix_norm
      integer :: i

      call diagonalize_nonhermitian(xi, diagnostics%eigenvalues, diagnostics%eigenvectors, diagnostics%left_eigenvectors)
      diagnostics%closest_eigenvalue_index = 1
      do i = 2, size(diagnostics%eigenvalues)
         if (abs(diagnostics%eigenvalues(i) - cmplx(1.0_rp, 0.0_rp, rp)) < &
             abs(diagnostics%eigenvalues(diagnostics%closest_eigenvalue_index) - cmplx(1.0_rp, 0.0_rp, rp))) then
            diagnostics%closest_eigenvalue_index = i
         end if
      end do
      diagnostics%closest_eigenvalue = diagnostics%eigenvalues(diagnostics%closest_eigenvalue_index)
      diagnostics%closest_eigenvalue_distance = abs(diagnostics%closest_eigenvalue - cmplx(1.0_rp, 0.0_rp, rp))
      allocate(diagnostics%closest_eigenvector(size(magnetization)))
      diagnostics%closest_eigenvector = diagnostics%eigenvectors(:, diagnostics%closest_eigenvalue_index)
      allocate(diagnostics%closest_left_eigenvector(size(magnetization)))
      diagnostics%closest_left_eigenvector = diagnostics%left_eigenvectors(:, diagnostics%closest_eigenvalue_index)
      norm_m = sqrt(sum(abs(magnetization)**2))
      norm_right = sqrt(sum(abs(diagnostics%closest_eigenvector)**2))
      norm_left = sqrt(sum(abs(diagnostics%closest_left_eigenvector)**2))
      diagnostics%magnetization_overlap = abs(dot_product(magnetization, diagnostics%closest_eigenvector))/ &
         (norm_m*norm_right)
      diagnostics%left_magnetization_overlap = abs(dot_product(diagnostics%closest_left_eigenvector, magnetization))/ &
         (norm_left*norm_m)
      diagnostics%biorthogonal_magnetization_overlap = abs(dot_product(diagnostics%closest_left_eigenvector, magnetization)* &
         dot_product(magnetization, diagnostics%closest_eigenvector))/(norm_left*norm_right*norm_m**2)
      allocate(residual_vector(size(magnetization)))
      residual_vector = matmul(xi, magnetization) - magnetization
      diagnostics%residual = sqrt(sum(abs(residual_vector)**2))/norm_m
      matrix_norm = sqrt(sum(abs(xi)**2))
      diagnostics%imaginary_norm = sqrt(sum(aimag(xi)**2))/max(1.0_rp, matrix_norm)
      diagnostics%available = .true.
      call evaluate_ward_from_xi(xi, magnetization, diagnostics%ward, response_basis='active response basis', &
         kernel_provenance='Xi = chi_KS K_xc; raw operator before any repair')
      if (present(bare_spectral_gap)) then
         diagnostics%has_bare_spectral_gap = .true.
         diagnostics%bare_spectral_gap = bare_spectral_gap
      end if
   end subroutine calculate_diagnostics

   subroutine diagonalize_nonhermitian(matrix, eigenvalues, eigenvectors, left_eigenvectors)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp), allocatable, intent(out) :: eigenvalues(:), eigenvectors(:, :), left_eigenvectors(:, :)
      complex(rp), allocatable :: work_matrix(:, :), work(:)
      real(rp), allocatable :: rwork(:)
      complex(rp) :: work_query(1)
      integer :: n, info, lwork

      n = size(matrix, 1)
      if (size(matrix, 2) /= n) error stop 'diagonalize_nonhermitian: matrix must be square'
      allocate(work_matrix(n, n), eigenvalues(n), eigenvectors(n, n), left_eigenvectors(n, n), rwork(max(1, 2*n)))
      work_matrix = matrix
      call zgeev('V', 'V', n, work_matrix, n, eigenvalues, left_eigenvectors, n, eigenvectors, n, work_query, -1, rwork, info)
      if (info /= 0) error stop 'diagonalize_nonhermitian: LAPACK workspace query failed'
      lwork = max(1, int(real(work_query(1), rp)))
      allocate(work(lwork))
      work_matrix = matrix
      call zgeev('V', 'V', n, work_matrix, n, eigenvalues, left_eigenvectors, n, eigenvectors, n, work, lwork, rwork, info)
      if (info /= 0) error stop 'diagonalize_nonhermitian: LAPACK zgeev failed'
   end subroutine diagonalize_nonhermitian

   integer function goldstone_mode_code(mode) result(code)
      character(len=*), intent(in) :: mode

      select case (trim(adjustl(mode)))
      case ('off')
         code = GOLDSTONE_OFF
      case ('diagnose')
         code = GOLDSTONE_DIAGNOSE
      case ('correct', 'sum_rule')
         code = GOLDSTONE_CORRECT
      case default
         error stop 'evaluate_goldstone: goldstone_mode must be off, diagnose, or correct'
      end select
   end function goldstone_mode_code

   subroutine reject_correction(correction, decision)
      type(tddft_goldstone_column_correction), intent(inout) :: correction
      character(len=*), intent(in) :: decision
      correction%rejected = .true.
      correction%applied = .false.
      correction%decision = trim(decision)
   end subroutine reject_correction

   elemental logical function finite_real(value)
      real(rp), intent(in) :: value
      finite_real = abs(value) < huge(value)
   end function finite_real

   subroutine require_square_response(chi_ks, k_perp, caller)
      complex(rp), intent(in) :: chi_ks(:, :), k_perp(:)
      character(len=*), intent(in) :: caller

      if (size(chi_ks, 1) /= size(chi_ks, 2) .or. size(chi_ks, 1) /= size(k_perp)) then
         error stop trim(caller)//': chi_KS and site K_perp dimensions are incompatible'
      end if
   end subroutine require_square_response

   subroutine write_one_diagnostics(unit, prefix, diagnostics)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: prefix
      type(tddft_goldstone_diagnostics), intent(in) :: diagnostics
      integer :: i

      write(unit, '(a,a,a,l1)') '# ', trim(prefix), '_available = ', diagnostics%available
      if (.not. diagnostics%available) return
      write(unit, '(a,a,a,2(1x,es24.16))') trim(prefix), '_closest_eigenvalue', '', &
         real(diagnostics%closest_eigenvalue, rp), aimag(diagnostics%closest_eigenvalue)
      write(unit, '(a,a,a,es24.16)') trim(prefix), '_closest_eigenvalue_distance', '', diagnostics%closest_eigenvalue_distance
      write(unit, '(a,a,a,es24.16)') trim(prefix), '_magnetization_overlap', '', diagnostics%magnetization_overlap
      write(unit, '(a,a,a,es24.16)') trim(prefix), '_left_magnetization_overlap', '', diagnostics%left_magnetization_overlap
      write(unit, '(a,a,a,es24.16)') trim(prefix), '_biorthogonal_magnetization_overlap', '', diagnostics%biorthogonal_magnetization_overlap
      write(unit, '(a,a,a,es24.16)') trim(prefix), '_imaginary_norm', '', diagnostics%imaginary_norm
      write(unit, '(a,a,a,es24.16)') trim(prefix), '_residual', '', diagnostics%residual
      if (diagnostics%ward%available) then
         write(unit, '(a,a,a,es24.16)') trim(prefix), '_ward_residual', '', diagnostics%ward%ward_residual
         write(unit, '(a,a,a,es24.16)') trim(prefix), '_dm_residual', '', diagnostics%ward%dm_residual
         write(unit, '(a,a,a,es24.16)') trim(prefix), '_bxc_kernel_residual', '', diagnostics%ward%bxc_kernel_residual
         write(unit, '(a,a,a,es24.16)') trim(prefix), '_magnetization_norm', '', diagnostics%ward%magnetization_norm
         write(unit, '(a,a,a,es24.16)') trim(prefix), '_bxc_norm', '', diagnostics%ward%bxc_norm
         write(unit, '(a,a,a,l1)') trim(prefix), '_identity_consistent', '', diagnostics%ward%identity_consistent
         write(unit, '(a,a,a)') trim(prefix), '_response_basis = ', trim(diagnostics%ward%response_basis)
         write(unit, '(a,a,a)') trim(prefix), '_bxc_provenance = ', trim(diagnostics%ward%bxc_provenance)
         write(unit, '(a,a,a)') trim(prefix), '_kernel_provenance = ', trim(diagnostics%ward%kernel_provenance)
      end if
      if (diagnostics%has_bare_spectral_gap) then
         write(unit, '(a,a,a,es24.16)') trim(prefix), '_bare_spectral_gap', '', diagnostics%bare_spectral_gap
      end if
      do i = 1, size(diagnostics%eigenvalues)
         write(unit, '(a,a,1x,i0,2(1x,es24.16))') trim(prefix), '_eigenvalue', i, real(diagnostics%eigenvalues(i), rp), &
            aimag(diagnostics%eigenvalues(i))
      end do
      do i = 1, size(diagnostics%closest_eigenvector)
         write(unit, '(a,a,1x,i0,2(1x,es24.16))') trim(prefix), '_closest_eigenvector', i, &
            real(diagnostics%closest_eigenvector(i), rp), aimag(diagnostics%closest_eigenvector(i))
      end do
      do i = 1, size(diagnostics%closest_left_eigenvector)
         write(unit, '(a,a,1x,i0,2(1x,es24.16))') trim(prefix), '_closest_left_eigenvector', i, &
            real(diagnostics%closest_left_eigenvector(i), rp), aimag(diagnostics%closest_left_eigenvector(i))
      end do
   end subroutine write_one_diagnostics

   subroutine write_ward_record(unit, prefix, diagnostics)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: prefix
      type(tddft_ward_diagnostics), intent(in) :: diagnostics

      write(unit, '(a,a,a,l1)') '# ', trim(prefix), '_available = ', diagnostics%available
      if (.not. diagnostics%available) return
      write(unit, '(a,a,a,es24.16)') '# ', trim(prefix), '_ward_residual = ', diagnostics%ward_residual
      write(unit, '(a,a,a,es24.16)') '# ', trim(prefix), '_dm_residual = ', diagnostics%dm_residual
      write(unit, '(a,a,a,a)') '# ', trim(prefix), '_response_basis = ', trim(diagnostics%response_basis)
      write(unit, '(a,a,a,a)') '# ', trim(prefix), '_kernel_provenance = ', trim(diagnostics%kernel_provenance)
   end subroutine write_ward_record

end module tddft_goldstone_mod
