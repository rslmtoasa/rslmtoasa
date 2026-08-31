!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Xi-based collective-mode diagnostics, biorthogonal branch tracking,
!> and controlled Lorentzian peak fitting.
!>
!> A collective candidate is a *crossing* of Re(lambda_Xi)=1, not merely an
!> eigenvalue which happened to be close to one on a discrete grid.  Xi is in
!> general non-normal, so every branch carries its left and right eigenvectors.
!> Branch continuation first follows frequency at a fixed q and then compares
!> adjacent q points using biorthogonal overlap and eigenvalue continuity.
!> Ill-conditioned eigensystems are reported as such; they are never forced
!> through an apparent exceptional point.
module tddft_modes_mod
   use precision_mod, only: rp
   implicit none

   private
   integer, parameter, public :: TDDFT_MODE_STONER = 1
   integer, parameter, public :: TDDFT_MODE_INCOHERENT = 2
   integer, parameter, public :: TDDFT_MODE_COHERENT = 3

   type, public :: tddft_mode_options
      ! Retained for input/source compatibility only.  Classification no
      ! longer uses distance from unity as a surrogate for a crossing.
      real(rp) :: unity_distance_threshold = 0.10_rp
      real(rp) :: maximum_fit_relative_residual = 0.10_rp
      integer :: minimum_points_per_fwhm = 2
      logical :: fit_isolated_peaks = .true.
      real(rp) :: maximum_biorthogonal_condition = 1.0e8_rp
      real(rp) :: minimum_branch_overlap = 0.25_rp
      real(rp) :: maximum_crossing_imaginary_part = 0.05_rp
   end type tddft_mode_options

   type, public :: tddft_peak_fit
      logical :: attempted = .false.
      logical :: accepted = .false.
      integer :: peak_index = 0
      real(rp) :: center = 0.0_rp
      real(rp) :: peak_height = 0.0_rp
      real(rp) :: background = 0.0_rp
      real(rp) :: fwhm = -1.0_rp
      real(rp) :: hwhm = -1.0_rp
      real(rp) :: relative_residual = huge(1.0_rp)
      character(len=96) :: rejection_reason = 'not attempted'
   end type tddft_peak_fit

   !> The historic `candidate_*` fields identify the grid point immediately
   !> adjacent to the selected interpolated unity crossing.  The crossing
   !> fields are the physical diagnostics used for classification.
   type, public :: tddft_mode_result
      complex(rp), allocatable :: xi_eigenvalues(:, :, :)
      complex(rp), allocatable :: xi_eigenvectors(:, :, :, :)
      complex(rp), allocatable :: xi_left_eigenvectors(:, :, :, :)
      real(rp), allocatable :: biorthogonal_condition(:, :, :)
      logical, allocatable :: eigensystem_well_conditioned(:, :, :)
      integer, allocatable :: candidate_frequency_index(:)
      integer, allocatable :: candidate_mode_index(:)
      real(rp), allocatable :: candidate_unity_distance(:)
      real(rp), allocatable :: branch_overlap(:)
      real(rp), allocatable :: branch_eigenvalue_step(:)
      logical, allocatable :: crossing_present(:)
      real(rp), allocatable :: crossing_omega(:)
      real(rp), allocatable :: crossing_imaginary_part(:)
      real(rp), allocatable :: mode_projected_weight(:)
      logical, allocatable :: exceptional_point_warning(:)
      integer, allocatable :: classification(:)
      character(len=48), allocatable :: classification_label(:)
      type(tddft_peak_fit), allocatable :: fit(:)
      real(rp) :: analysis_cpu_seconds = 0.0_rp
   end type tddft_mode_result

   public :: analyze_tddft_modes
   public :: fit_isolated_lorentzian_peak
   public :: extrapolate_linewidth_zero_eta
   public :: write_tddft_modes_text

contains

   !> Analyze a q path.  Xi has shape (response,response,omega,q); trace_loss
   !> has shape (omega,q).  When supplied, dressed_loss is the corresponding
   !> Hermitian loss matrix and permits a biorthogonal mode projection:
   !> Re[l^H L r/(l^H r)].  The caller supplies q in physical path order.
   subroutine analyze_tddft_modes(omega, xi, trace_loss, eta, options, result, dressed_loss)
      real(rp), intent(in) :: omega(:), trace_loss(:, :), eta
      complex(rp), intent(in) :: xi(:, :, :, :)
      type(tddft_mode_options), intent(in) :: options
      type(tddft_mode_result), intent(out) :: result
      complex(rp), intent(in), optional :: dressed_loss(:, :, :, :)
      integer :: n, nw, nq, iq, iw, previous_mode, selected_mode, selected_iw
      real(rp) :: t_start, t_stop

      n = size(xi, 1)
      nw = size(xi, 3)
      nq = size(xi, 4)
      if (n <= 0 .or. nw < 2 .or. nq <= 0 .or. size(xi, 2) /= n .or. size(omega) /= nw .or. &
          any(shape(trace_loss) /= [nw, nq])) then
         error stop 'analyze_tddft_modes: incompatible Xi, omega, or loss dimensions'
      end if
      if (eta <= 0.0_rp) error stop 'analyze_tddft_modes: eta must be a positive numerical broadening'
      if (options%unity_distance_threshold <= 0.0_rp .or. options%maximum_fit_relative_residual <= 0.0_rp) then
         error stop 'analyze_tddft_modes: positive diagnostic tolerances required'
      end if
      if (options%maximum_biorthogonal_condition <= 1.0_rp .or. options%minimum_branch_overlap <= 0.0_rp .or. &
          options%minimum_branch_overlap > 1.0_rp .or. options%maximum_crossing_imaginary_part < 0.0_rp) then
         error stop 'analyze_tddft_modes: invalid biorthogonal tracking tolerances'
      end if
      if (any(omega(2:) <= omega(:size(omega)-1))) error stop 'analyze_tddft_modes: omega must increase strictly'
      if (present(dressed_loss)) then
         if (any(shape(dressed_loss) /= [n, n, nw, nq])) then
            error stop 'analyze_tddft_modes: dressed loss shape mismatch'
         end if
      end if

      call cpu_time(t_start)
      allocate(result%xi_eigenvalues(n, nw, nq), result%xi_eigenvectors(n, n, nw, nq), &
         result%xi_left_eigenvectors(n, n, nw, nq), result%biorthogonal_condition(n, nw, nq), &
         result%eigensystem_well_conditioned(n, nw, nq), &
         result%candidate_frequency_index(nq), result%candidate_mode_index(nq), result%candidate_unity_distance(nq), &
         result%branch_overlap(nq), result%branch_eigenvalue_step(nq), result%crossing_present(nq), &
         result%crossing_omega(nq), result%crossing_imaginary_part(nq), result%mode_projected_weight(nq), &
         result%exceptional_point_warning(nq), result%classification(nq), result%classification_label(nq), result%fit(nq))
      do iq = 1, nq
         do iw = 1, nw
            call diagonalize_biorthogonal(xi(:, :, iw, iq), result%xi_eigenvalues(:, iw, iq), &
               result%xi_eigenvectors(:, :, iw, iq), result%xi_left_eigenvectors(:, :, iw, iq), &
               result%biorthogonal_condition(:, iw, iq), result%eigensystem_well_conditioned(:, iw, iq), &
               options%maximum_biorthogonal_condition)
         end do
         call order_frequency_branches(result%xi_eigenvalues(:, :, iq), result%xi_eigenvectors(:, :, :, iq), &
            result%xi_left_eigenvectors(:, :, :, iq), result%biorthogonal_condition(:, :, iq), &
            result%eigensystem_well_conditioned(:, :, iq))
      end do

      do iq = 1, nq
         if (iq == 1) then
            previous_mode = 0
         else
            previous_mode = result%candidate_mode_index(iq-1)
         end if
         call select_unity_crossing(omega, result, iq, previous_mode, options, selected_mode, selected_iw)
         result%candidate_mode_index(iq) = selected_mode
         result%candidate_frequency_index(iq) = selected_iw
         result%candidate_unity_distance(iq) = abs(result%xi_eigenvalues(selected_mode, selected_iw, iq)-cmplx(1.0_rp, 0.0_rp, rp))
         result%exceptional_point_warning(iq) = .not. result%eigensystem_well_conditioned(selected_mode, selected_iw, iq)
         if (present(dressed_loss)) then
            result%mode_projected_weight(iq) = project_mode_weight(dressed_loss(:, :, selected_iw, iq), &
               result%xi_left_eigenvectors(:, selected_mode, selected_iw, iq), &
               result%xi_eigenvectors(:, selected_mode, selected_iw, iq))
         else
            result%mode_projected_weight(iq) = -1.0_rp
         end if
         if (options%fit_isolated_peaks .and. result%crossing_present(iq)) then
            call fit_isolated_lorentzian_peak(omega, trace_loss(:, iq), options, result%fit(iq))
         else
            result%fit(iq)%rejection_reason = 'no Re(lambda_Xi)=1 crossing; total loss retained as Stoner evidence'
         end if
         call classify_mode(result%crossing_present(iq), result%crossing_imaginary_part(iq), &
            result%exceptional_point_warning(iq), result%mode_projected_weight(iq), options, result%fit(iq), &
            result%classification(iq), result%classification_label(iq))
      end do
      call cpu_time(t_stop)
      result%analysis_cpu_seconds = t_stop-t_start
   end subroutine analyze_tddft_modes

   !> LAPACK returns unit-norm left/right eigenvectors.  Rescale the left
   !> vector so l^H r=1 when it is safe to do so.  Its original reciprocal
   !> overlap is the eigenpair condition diagnostic used to reject a forced
   !> continuation through a near-defective Xi.
   subroutine diagonalize_biorthogonal(matrix, eigenvalues, right_vectors, left_vectors, condition, well_conditioned, maximum_condition)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp), intent(out) :: eigenvalues(:), right_vectors(:, :), left_vectors(:, :)
      real(rp), intent(out) :: condition(:)
      logical, intent(out) :: well_conditioned(:)
      real(rp), intent(in) :: maximum_condition
      complex(rp), allocatable :: work_matrix(:, :), work(:)
      real(rp), allocatable :: rwork(:)
      complex(rp) :: work_query(1), overlap
      integer :: n, info, lwork, imode

      n = size(matrix, 1)
      if (size(matrix, 2) /= n .or. size(eigenvalues) /= n .or. any(shape(right_vectors) /= [n, n]) .or. &
          any(shape(left_vectors) /= [n, n]) .or. size(condition) /= n .or. size(well_conditioned) /= n) then
         error stop 'diagonalize_biorthogonal: incompatible dimensions'
      end if
      allocate(work_matrix(n, n), work(1), rwork(max(1, 2*n)))
      work_matrix = matrix
      call zgeev('V', 'V', n, work_matrix, n, eigenvalues, left_vectors, n, right_vectors, n, work_query, -1, rwork, info)
      if (info /= 0) error stop 'diagonalize_biorthogonal: LAPACK workspace query failed'
      lwork = max(1, int(real(work_query(1), rp)))
      deallocate(work); allocate(work(lwork))
      work_matrix = matrix
      call zgeev('V', 'V', n, work_matrix, n, eigenvalues, left_vectors, n, right_vectors, n, work, lwork, rwork, info)
      if (info /= 0) error stop 'diagonalize_biorthogonal: LAPACK zgeev failed'
      do imode = 1, n
         overlap = dot_product(left_vectors(:, imode), right_vectors(:, imode))
         condition(imode) = 1.0_rp/max(abs(overlap), tiny(1.0_rp))
         well_conditioned(imode) = condition(imode) <= maximum_condition
         if (well_conditioned(imode)) then
            left_vectors(:, imode) = left_vectors(:, imode)/conjg(overlap)
         end if
      end do
   end subroutine diagonalize_biorthogonal

   !> Greedy frequency-local assignment.  The score has a biorthogonal
   !> overlap factor and an eigenvalue-continuity factor, so a remote grid
   !> point with a superficially larger right-vector overlap cannot relabel a
   !> branch.  The response bases used here are small; a deterministic greedy
   !> assignment also keeps branch IDs stable in degenerate synthetic tests.
   subroutine order_frequency_branches(eigenvalues, right_vectors, left_vectors, condition, well_conditioned)
      complex(rp), intent(inout) :: eigenvalues(:, :), right_vectors(:, :, :), left_vectors(:, :, :)
      real(rp), intent(inout) :: condition(:, :)
      logical, intent(inout) :: well_conditioned(:, :)
      integer :: n, nw, iw, imode, jmode, best, permutation(size(eigenvalues, 1))
      logical :: used(size(eigenvalues, 1))
      real(rp) :: best_score, score
      complex(rp) :: values_copy(size(eigenvalues, 1)), right_copy(size(right_vectors, 1), size(right_vectors, 2)), &
         left_copy(size(left_vectors, 1), size(left_vectors, 2))
      real(rp) :: condition_copy(size(condition, 1))
      logical :: conditioned_copy(size(well_conditioned, 1))

      n = size(eigenvalues, 1); nw = size(eigenvalues, 2)
      if (n <= 0 .or. size(right_vectors, 2) /= n .or. size(left_vectors, 2) /= n .or. &
          any(shape(right_vectors) /= [n, n, nw]) .or. any(shape(left_vectors) /= [n, n, nw]) .or. &
          any(shape(condition) /= [n, nw]) .or. any(shape(well_conditioned) /= [n, nw])) then
         error stop 'order_frequency_branches: incompatible dimensions'
      end if
      do iw = 2, nw
         used = .false.
         do imode = 1, n
            best = 0; best_score = -huge(1.0_rp)
            do jmode = 1, n
               if (used(jmode)) cycle
               score = biorthogonal_overlap(left_vectors(:, imode, iw-1), right_vectors(:, imode, iw-1), &
                  left_vectors(:, jmode, iw), right_vectors(:, jmode, iw))* &
                  exp(-abs(eigenvalues(imode, iw-1)-eigenvalues(jmode, iw)))
               if (score > best_score) then
                  best_score = score; best = jmode
               end if
            end do
            permutation(imode) = best; used(best) = .true.
         end do
         values_copy = eigenvalues(:, iw); right_copy = right_vectors(:, :, iw); left_copy = left_vectors(:, :, iw)
         condition_copy = condition(:, iw); conditioned_copy = well_conditioned(:, iw)
         do imode = 1, n
            eigenvalues(imode, iw) = values_copy(permutation(imode))
            right_vectors(:, imode, iw) = right_copy(:, permutation(imode))
            left_vectors(:, imode, iw) = left_copy(:, permutation(imode))
            condition(imode, iw) = condition_copy(permutation(imode))
            well_conditioned(imode, iw) = conditioned_copy(permutation(imode))
         end do
      end do
   end subroutine order_frequency_branches

   !> Choose an actual unity crossing.  At q>q_1 candidates are scored with
   !> the preceding selected mode’s biorthogonal overlap and its interpolated
   !> crossing eigenvalue.  No candidate is invented when no crossing exists.
   subroutine select_unity_crossing(omega, result, iq, previous_mode, options, selected_mode, selected_iw)
      real(rp), intent(in) :: omega(:)
      type(tddft_mode_result), intent(inout) :: result
      integer, intent(in) :: iq, previous_mode
      type(tddft_mode_options), intent(in) :: options
      integer, intent(out) :: selected_mode, selected_iw
      integer :: imode, iw, n, nw, best_mode, best_iw, nearest_iw
      real(rp) :: fraction, crossing_imag, score, best_score, overlap, eigen_step, nearest_distance, distance
      complex(rp) :: crossing_value, previous_value
      logical :: found, eligible

      n = size(result%xi_eigenvalues, 1); nw = size(omega)
      found = .false.; best_score = -huge(1.0_rp); best_mode = 1; best_iw = 1
      result%branch_overlap(iq) = 1.0_rp; result%branch_eigenvalue_step(iq) = 0.0_rp
      result%crossing_present(iq) = .false.; result%crossing_omega(iq) = -1.0_rp
      result%crossing_imaginary_part(iq) = huge(1.0_rp)
      do imode = 1, n
         do iw = 1, nw-1
            if (.not. real_unity_crossing(real(result%xi_eigenvalues(imode, iw, iq), rp), &
               real(result%xi_eigenvalues(imode, iw+1, iq), rp))) cycle
            eligible = result%eigensystem_well_conditioned(imode, iw, iq) .and. &
               result%eigensystem_well_conditioned(imode, iw+1, iq)
            if (.not. eligible) cycle
            fraction = crossing_fraction(real(result%xi_eigenvalues(imode, iw, iq), rp), &
               real(result%xi_eigenvalues(imode, iw+1, iq), rp))
            crossing_value = (1.0_rp-fraction)*result%xi_eigenvalues(imode, iw, iq) + &
               fraction*result%xi_eigenvalues(imode, iw+1, iq)
            crossing_imag = aimag(crossing_value)
            nearest_iw = iw
            if (fraction > 0.5_rp) nearest_iw = iw+1
            if (previous_mode == 0) then
               score = -abs(crossing_imag)
               overlap = 1.0_rp; eigen_step = 0.0_rp
            else
               previous_value = cmplx(1.0_rp, result%crossing_imaginary_part(iq-1), rp)
               overlap = biorthogonal_overlap(result%xi_left_eigenvectors(:, previous_mode, &
                  result%candidate_frequency_index(iq-1), iq-1), result%xi_eigenvectors(:, previous_mode, &
                  result%candidate_frequency_index(iq-1), iq-1), result%xi_left_eigenvectors(:, imode, nearest_iw, iq), &
                  result%xi_eigenvectors(:, imode, nearest_iw, iq))
               eigen_step = abs(crossing_value-previous_value)
               if (overlap < options%minimum_branch_overlap) cycle
               score = overlap*exp(-eigen_step)
            end if
            if (.not. found .or. score > best_score) then
               found = .true.; best_score = score; best_mode = imode; best_iw = nearest_iw
               result%crossing_omega(iq) = omega(iw) + fraction*(omega(iw+1)-omega(iw))
               result%crossing_imaginary_part(iq) = crossing_imag
               result%branch_overlap(iq) = overlap; result%branch_eigenvalue_step(iq) = eigen_step
            end if
         end do
      end do
      if (found) then
         result%crossing_present(iq) = .true.; selected_mode = best_mode; selected_iw = best_iw
      else
         ! Preserve a transparent nearest-grid diagnostic for a noncollective
         ! feature, but never label it as a unity crossing.
         nearest_distance = huge(1.0_rp); selected_mode = 1; selected_iw = 1
         do iw = 1, nw
            do imode = 1, n
               distance = abs(result%xi_eigenvalues(imode, iw, iq)-cmplx(1.0_rp, 0.0_rp, rp))
               if (distance < nearest_distance) then
                  nearest_distance = distance; selected_mode = imode; selected_iw = iw
               end if
            end do
         end do
      end if
   end subroutine select_unity_crossing

   !> Fit a local loss maximum after estimating a linear-free constant
   !> background from its bracketing minima.  The FWHM/HWHM are explicitly
   !> reported.  This routine reports *observed* width; eta is not subtracted.
   subroutine fit_isolated_lorentzian_peak(omega, spectrum, options, fit)
      real(rp), intent(in) :: omega(:), spectrum(:)
      type(tddft_mode_options), intent(in) :: options
      type(tddft_peak_fit), intent(out) :: fit
      integer :: n, ip, ileft, iright, i, nfit, left_bracket, right_bracket
      real(rp) :: background, amplitude, half_height, left_cross, right_cross, gamma, model, sumsq, sumdata

      n = size(omega)
      fit%attempted = .true.
      fit%rejection_reason = 'no isolated local maximum'
      if (n < 5 .or. size(spectrum) /= n) return
      ip = 2
      do i = 3, n - 1
         if (spectrum(i) > spectrum(ip)) ip = i
      end do
      if (spectrum(ip) <= spectrum(ip-1) .or. spectrum(ip) <= spectrum(ip+1)) return
      fit%peak_index = ip
      ileft = ip - 1
      ! Fortran does not guarantee short-circuit evaluation of `.and.`.
      ! Keep the bound check separate so a peak adjacent to the first point
      ! never evaluates spectrum(0) under bounds checking.
      do while (ileft > 1)
         if (spectrum(ileft-1) > spectrum(ileft)) exit
         ileft = ileft - 1
      end do
      iright = ip + 1
      do while (iright < n)
         if (spectrum(iright+1) > spectrum(iright)) exit
         iright = iright + 1
      end do
      background = min(spectrum(ileft), spectrum(iright))
      amplitude = spectrum(ip) - background
      if (amplitude <= tiny(1.0_rp)) then
         fit%rejection_reason = 'peak has no positive contrast over local background'; return
      end if
      half_height = background + 0.5_rp*amplitude
      if (spectrum(ileft) > half_height .or. spectrum(iright) > half_height) then
         fit%rejection_reason = 'peak is not bracketed by half-height crossings'; return
      end if
      left_bracket = 0
      do i = ip - 1, ileft, -1
         if (spectrum(i) <= half_height .and. spectrum(i+1) > half_height) then
            left_bracket = i; exit
         end if
      end do
      right_bracket = 0
      do i = ip + 1, iright
         if (spectrum(i) <= half_height .and. spectrum(i-1) > half_height) then
            right_bracket = i; exit
         end if
      end do
      if (left_bracket == 0 .or. right_bracket == 0) then
         fit%rejection_reason = 'half-height crossings are not monotonic around peak'; return
      end if
      left_cross = interpolate_crossing(omega(left_bracket), spectrum(left_bracket), omega(left_bracket+1), &
         spectrum(left_bracket+1), half_height)
      right_cross = interpolate_crossing(omega(right_bracket-1), spectrum(right_bracket-1), omega(right_bracket), &
         spectrum(right_bracket), half_height)
      fit%center = omega(ip)
      fit%peak_height = spectrum(ip)
      fit%background = background
      fit%fwhm = right_cross - left_cross
      fit%hwhm = 0.5_rp*fit%fwhm
      if (fit%fwhm <= 0.0_rp) then
         fit%rejection_reason = 'nonpositive FWHM'; return
      end if
      if (fit%fwhm < real(options%minimum_points_per_fwhm, rp)*minval(omega(2:)-omega(:n-1))) then
         fit%rejection_reason = 'FWHM is unresolved on the supplied frequency grid'; return
      end if
      gamma = fit%hwhm
      sumsq = 0.0_rp; sumdata = 0.0_rp; nfit = 0
      do i = ileft, iright
         model = background + amplitude*gamma**2/((omega(i)-fit%center)**2 + gamma**2)
         sumsq = sumsq + (spectrum(i)-model)**2
         sumdata = sumdata + (spectrum(i)-background)**2
         nfit = nfit + 1
      end do
      fit%relative_residual = sqrt(sumsq/max(sumdata, tiny(1.0_rp)))
      if (fit%relative_residual > options%maximum_fit_relative_residual) then
         fit%rejection_reason = 'non-Lorentzian residual exceeds controlled threshold'; return
      end if
      fit%accepted = .true.
      fit%rejection_reason = 'accepted; observed FWHM/HWHM include numerical eta'
   end subroutine fit_isolated_lorentzian_peak

   !> Linear multi-eta workflow: observed FWHM(eta)=intercept+slope*eta.
   !> The intercept is the zero-numerical-broadening extrapolation; intrinsic
   !> HWHM is intercept/2, not the eta of any individual run.
   subroutine extrapolate_linewidth_zero_eta(eta, observed_fwhm, intrinsic_fwhm, intrinsic_hwhm, slope, relative_residual, accepted)
      real(rp), intent(in) :: eta(:), observed_fwhm(:)
      real(rp), intent(out) :: intrinsic_fwhm, intrinsic_hwhm, slope, relative_residual
      logical, intent(out) :: accepted
      integer :: n, i
      real(rp) :: sx, sy, sxx, sxy, denominator, fitted, ssres, sstot

      n = size(eta)
      accepted = .false.; intrinsic_fwhm = -1.0_rp; intrinsic_hwhm = -1.0_rp; slope = 0.0_rp; relative_residual = huge(1.0_rp)
      if (n < 2 .or. size(observed_fwhm) /= n .or. any(eta <= 0.0_rp) .or. any(observed_fwhm <= 0.0_rp)) return
      sx = sum(eta); sy = sum(observed_fwhm); sxx = sum(eta**2); sxy = sum(eta*observed_fwhm)
      denominator = real(n, rp)*sxx - sx**2
      if (abs(denominator) <= tiny(1.0_rp)) return
      slope = (real(n, rp)*sxy-sx*sy)/denominator
      intrinsic_fwhm = (sy-slope*sx)/real(n, rp)
      intrinsic_hwhm = 0.5_rp*intrinsic_fwhm
      ssres = 0.0_rp; sstot = 0.0_rp
      do i = 1, n
         fitted = intrinsic_fwhm + slope*eta(i)
         ssres = ssres + (observed_fwhm(i)-fitted)**2
         sstot = sstot + (observed_fwhm(i)-sy/real(n, rp))**2
      end do
      relative_residual = sqrt(ssres/max(sstot, tiny(1.0_rp)))
      accepted = intrinsic_fwhm >= 0.0_rp .and. relative_residual <= 0.10_rp
   end subroutine extrapolate_linewidth_zero_eta

   subroutine write_tddft_modes_text(filename, omega, eta, result)
      character(len=*), intent(in) :: filename
      real(rp), intent(in) :: omega(:), eta
      type(tddft_mode_result), intent(in) :: result
      integer :: unit, ios, iq

      if (.not. allocated(result%candidate_frequency_index)) error stop 'write_tddft_modes_text: no mode result'
      open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
      if (ios /= 0) error stop 'write_tddft_modes_text: cannot open output file'
      write(unit, '(a)') '# Xi branches are frequency-local biorthogonal continuations, then continued between adjacent q points.'
      write(unit, '(a)') '# A collective candidate requires an interpolated Re(lambda_Xi)=1 crossing; Im(lambda), conditioning, and mode weight are retained.'
      write(unit, '(a)') '# FWHM = 2*HWHM. Individual fits are observed widths; eta remains numerical broadening.'
      write(unit, '(a,es24.16)') '# eta_Ry = ', eta
      write(unit, '(a,es24.16)') '# profile_mode_analysis_cpu_s = ', result%analysis_cpu_seconds
      do iq = 1, size(result%candidate_frequency_index)
         write(unit, '(a,1x,i0,1x,es24.16,1x,i0,1x,es24.16,1x,es24.16,1x,a)') 'candidate', iq, &
            omega(result%candidate_frequency_index(iq)), result%candidate_mode_index(iq), result%candidate_unity_distance(iq), &
            result%branch_overlap(iq), trim(result%classification_label(iq))
         write(unit, '(a,1x,i0,1x,l1,6(1x,es24.16),1x,l1)') 'crossing', iq, result%crossing_present(iq), &
            result%crossing_omega(iq), result%crossing_imaginary_part(iq), result%branch_overlap(iq), &
            result%branch_eigenvalue_step(iq), result%mode_projected_weight(iq), &
            result%biorthogonal_condition(result%candidate_mode_index(iq), result%candidate_frequency_index(iq), iq), &
            result%exceptional_point_warning(iq)
         write(unit, '(a,1x,i0,1x,l1,4(1x,es24.16),1x,a)') 'fit', iq, result%fit(iq)%accepted, result%fit(iq)%center, &
            result%fit(iq)%fwhm, result%fit(iq)%hwhm, result%fit(iq)%relative_residual, trim(result%fit(iq)%rejection_reason)
      end do
      close(unit)
   end subroutine write_tddft_modes_text

   subroutine classify_mode(crossing_present, crossing_imaginary_part, exceptional_warning, projected_weight, options, fit, code, label)
      logical, intent(in) :: crossing_present, exceptional_warning
      real(rp), intent(in) :: crossing_imaginary_part, projected_weight
      type(tddft_mode_options), intent(in) :: options
      type(tddft_peak_fit), intent(in) :: fit
      integer, intent(out) :: code
      character(len=*), intent(out) :: label

      if (.not. crossing_present) then
         code = TDDFT_MODE_STONER; label = 'noncollective Stoner feature'
      else if (exceptional_warning) then
         code = TDDFT_MODE_INCOHERENT; label = 'ill-conditioned Xi crossing'
      else if (abs(crossing_imaginary_part) > options%maximum_crossing_imaginary_part) then
         code = TDDFT_MODE_INCOHERENT; label = 'strongly damped Xi crossing'
      else if (projected_weight == 0.0_rp) then
         code = TDDFT_MODE_INCOHERENT; label = 'zero mode-projected enhanced weight'
      else
         code = TDDFT_MODE_INCOHERENT; label = 'strongly Landau-damped/incoherent enhanced mode'
         if (fit%accepted .or. projected_weight > 0.0_rp) then
            code = TDDFT_MODE_COHERENT; label = 'collective Xi unity crossing'
         end if
      end if
   end subroutine classify_mode

   pure real(rp) function biorthogonal_overlap(left_a, right_a, left_b, right_b) result(overlap)
      complex(rp), intent(in) :: left_a(:), right_a(:), left_b(:), right_b(:)
      real(rp) :: denominator
      denominator = sqrt(abs(dot_product(left_a, right_a))*abs(dot_product(left_b, right_b)))
      if (denominator <= tiny(1.0_rp)) then
         overlap = 0.0_rp
      else
         overlap = sqrt(abs(dot_product(left_a, right_b))*abs(dot_product(left_b, right_a)))/denominator
      end if
   end function biorthogonal_overlap

   pure logical function real_unity_crossing(left_value, right_value) result(crosses)
      real(rp), intent(in) :: left_value, right_value
      crosses = (left_value <= 1.0_rp .and. right_value >= 1.0_rp) .or. &
         (left_value >= 1.0_rp .and. right_value <= 1.0_rp)
   end function real_unity_crossing

   pure real(rp) function crossing_fraction(left_value, right_value) result(fraction)
      real(rp), intent(in) :: left_value, right_value
      if (abs(right_value-left_value) <= tiny(1.0_rp)) then
         fraction = 0.5_rp
      else
         fraction = max(0.0_rp, min(1.0_rp, (1.0_rp-left_value)/(right_value-left_value)))
      end if
   end function crossing_fraction

   pure real(rp) function project_mode_weight(loss, left, right) result(weight)
      complex(rp), intent(in) :: loss(:, :), left(:), right(:)
      complex(rp) :: denominator
      if (size(loss, 1) /= size(left) .or. size(loss, 2) /= size(right)) then
         weight = -1.0_rp; return
      end if
      denominator = dot_product(left, right)
      if (abs(denominator) <= tiny(1.0_rp)) then
         weight = -1.0_rp
      else
         weight = real(dot_product(left, matmul(loss, right))/denominator, rp)
      end if
   end function project_mode_weight

   pure real(rp) function interpolate_crossing(x1, y1, x2, y2, target) result(x)
      real(rp), intent(in) :: x1, y1, x2, y2, target
      if (abs(y2-y1) <= tiny(1.0_rp)) then
         x = 0.5_rp*(x1+x2)
      else
         x = x1 + (target-y1)*(x2-x1)/(y2-y1)
      end if
   end function interpolate_crossing

end module tddft_modes_mod
