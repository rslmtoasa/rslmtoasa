!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Xi-based collective-mode diagnostics, overlap branch tracking, and
!> controlled Lorentzian peak fitting.
!>
!> A candidate is selected from the Xi eigenvalue closest to one, then its
!> branch is connected from q to q by maximum right-eigenvector overlap.  Peak
!> fitting is intentionally conservative: absent half-height crossings,
!> unresolved widths, or a large Lorentzian residual leave the fit rejected.
!> Rejection is a classification result, not a request to force a magnon fit.
module tddft_modes_mod
   use precision_mod, only: rp
   use tddft_dyson_mod, only: diagonalize_nonhermitian_response
   implicit none

   private
   integer, parameter, public :: TDDFT_MODE_STONER = 1
   integer, parameter, public :: TDDFT_MODE_INCOHERENT = 2
   integer, parameter, public :: TDDFT_MODE_COHERENT = 3

   type, public :: tddft_mode_options
      real(rp) :: unity_distance_threshold = 0.10_rp
      real(rp) :: maximum_fit_relative_residual = 0.10_rp
      integer :: minimum_points_per_fwhm = 2
      logical :: fit_isolated_peaks = .true.
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

   !> The branch selected at each q is the frequency/eigenvector point with
   !> smallest distance to Xi=1, subject to overlap continuation where a prior
   !> q exists.  The complete Xi eigensystem is retained for alternative
   !> selection policies and postprocessing.
   type, public :: tddft_mode_result
      complex(rp), allocatable :: xi_eigenvalues(:, :, :)
      complex(rp), allocatable :: xi_eigenvectors(:, :, :, :)
      integer, allocatable :: candidate_frequency_index(:)
      integer, allocatable :: candidate_mode_index(:)
      real(rp), allocatable :: candidate_unity_distance(:)
      real(rp), allocatable :: branch_overlap(:)
      integer, allocatable :: classification(:)
      character(len=48), allocatable :: classification_label(:)
      type(tddft_peak_fit), allocatable :: fit(:)
   end type tddft_mode_result

   public :: analyze_tddft_modes
   public :: fit_isolated_lorentzian_peak
   public :: extrapolate_linewidth_zero_eta
   public :: write_tddft_modes_text

contains

   !> Analyze a q path.  Xi has shape (response,response,omega,q); trace_loss
   !> has shape (omega,q).  The caller supplies a q path in physical order.
   subroutine analyze_tddft_modes(omega, xi, trace_loss, eta, options, result)
      real(rp), intent(in) :: omega(:), trace_loss(:, :), eta
      complex(rp), intent(in) :: xi(:, :, :, :)
      type(tddft_mode_options), intent(in) :: options
      type(tddft_mode_result), intent(out) :: result
      integer :: n, nw, nq, iq, iw, imode, best_iw, best_mode
      real(rp) :: best_distance, distance, previous_overlap, candidate_overlap
      complex(rp), allocatable :: previous_vector(:)

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
      if (any(omega(2:) <= omega(:size(omega)-1))) error stop 'analyze_tddft_modes: omega must increase strictly'

      allocate(result%xi_eigenvalues(n, nw, nq), result%xi_eigenvectors(n, n, nw, nq), &
         result%candidate_frequency_index(nq), result%candidate_mode_index(nq), result%candidate_unity_distance(nq), &
         result%branch_overlap(nq), result%classification(nq), result%classification_label(nq), result%fit(nq))
      do iq = 1, nq
         do iw = 1, nw
            call diagonalize_nonhermitian_response(xi(:, :, iw, iq), result%xi_eigenvalues(:, iw, iq), &
               result%xi_eigenvectors(:, :, iw, iq))
         end do
      end do

      do iq = 1, nq
         best_distance = huge(1.0_rp)
         best_iw = 1
         best_mode = 1
         previous_overlap = -1.0_rp
         do iw = 1, nw
            do imode = 1, n
               distance = abs(result%xi_eigenvalues(imode, iw, iq) - cmplx(1.0_rp, 0.0_rp, rp))
               if (iq == 1) then
                  if (distance < best_distance) then
                     best_distance = distance; best_iw = iw; best_mode = imode
                  end if
               else
                  candidate_overlap = normalized_overlap(previous_vector, result%xi_eigenvectors(:, imode, iw, iq))
                  ! First preserve the continued branch by overlap.  Distance
                  ! breaks a nearly degenerate overlap tie and remains reported.
                  if (candidate_overlap > previous_overlap + 100.0_rp*epsilon(1.0_rp) .or. &
                      (abs(candidate_overlap-previous_overlap) <= 100.0_rp*epsilon(1.0_rp) .and. distance < best_distance)) then
                     previous_overlap = candidate_overlap; best_distance = distance; best_iw = iw; best_mode = imode
                  end if
               end if
            end do
         end do
         result%candidate_frequency_index(iq) = best_iw
         result%candidate_mode_index(iq) = best_mode
         result%candidate_unity_distance(iq) = best_distance
         if (iq == 1) then
            result%branch_overlap(iq) = 1.0_rp
         else
            result%branch_overlap(iq) = previous_overlap
         end if
         if (allocated(previous_vector)) deallocate(previous_vector)
         allocate(previous_vector(n))
         previous_vector = result%xi_eigenvectors(:, best_mode, best_iw, iq)

         if (options%fit_isolated_peaks .and. best_distance <= options%unity_distance_threshold) then
            call fit_isolated_lorentzian_peak(omega, trace_loss(:, iq), options, result%fit(iq))
         else
            result%fit(iq)%rejection_reason = 'no Xi unity candidate; retain as predominantly Stoner'
         end if
         call classify_mode(best_distance, options%unity_distance_threshold, result%fit(iq), result%classification(iq), &
            result%classification_label(iq))
      end do
   end subroutine analyze_tddft_modes

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
      do while (ileft > 1 .and. spectrum(ileft-1) <= spectrum(ileft))
         ileft = ileft - 1
      end do
      iright = ip + 1
      do while (iright < n .and. spectrum(iright+1) <= spectrum(iright))
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
      write(unit, '(a)') '# Xi candidates use right-eigenvector overlap between adjacent q values'
      write(unit, '(a)') '# FWHM = 2*HWHM. Individual fits are observed widths; eta remains numerical broadening.'
      write(unit, '(a,es24.16)') '# eta_Ry = ', eta
      do iq = 1, size(result%candidate_frequency_index)
         write(unit, '(a,1x,i0,1x,es24.16,1x,i0,1x,es24.16,1x,es24.16,1x,a)') 'candidate', iq, &
            omega(result%candidate_frequency_index(iq)), result%candidate_mode_index(iq), result%candidate_unity_distance(iq), &
            result%branch_overlap(iq), trim(result%classification_label(iq))
         write(unit, '(a,1x,i0,1x,l1,4(1x,es24.16),1x,a)') 'fit', iq, result%fit(iq)%accepted, result%fit(iq)%center, &
            result%fit(iq)%fwhm, result%fit(iq)%hwhm, result%fit(iq)%relative_residual, trim(result%fit(iq)%rejection_reason)
      end do
      close(unit)
   end subroutine write_tddft_modes_text

   subroutine classify_mode(unity_distance, threshold, fit, code, label)
      real(rp), intent(in) :: unity_distance, threshold
      type(tddft_peak_fit), intent(in) :: fit
      integer, intent(out) :: code
      character(len=*), intent(out) :: label

      if (unity_distance > threshold) then
         code = TDDFT_MODE_STONER; label = 'predominantly Stoner feature'
      else if (fit%accepted) then
         code = TDDFT_MODE_COHERENT; label = 'coherent collective magnon'
      else
         code = TDDFT_MODE_INCOHERENT; label = 'strongly Landau-damped/incoherent enhanced mode'
      end if
   end subroutine classify_mode

   pure real(rp) function normalized_overlap(left, right) result(overlap)
      complex(rp), intent(in) :: left(:), right(:)
      real(rp) :: norm_left, norm_right
      norm_left = sqrt(sum(abs(left)**2)); norm_right = sqrt(sum(abs(right)**2))
      if (norm_left <= tiny(1.0_rp) .or. norm_right <= tiny(1.0_rp)) then
         overlap = 0.0_rp
      else
         overlap = abs(dot_product(left, right))/(norm_left*norm_right)
      end if
   end function normalized_overlap

   pure real(rp) function interpolate_crossing(x1, y1, x2, y2, target) result(x)
      real(rp), intent(in) :: x1, y1, x2, y2, target
      if (abs(y2-y1) <= tiny(1.0_rp)) then
         x = 0.5_rp*(x1+x2)
      else
         x = x1 + (target-y1)*(x2-x1)/(y2-y1)
      end if
   end function interpolate_crossing

end module tddft_modes_mod
