!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief TDDFT Dyson enhancement and loss products for a response-frequency batch.
!>
!> The retarded convention is inherited from tddft_chi0_mod.  At positive
!> frequency the reported loss matrix is
!>
!>   L = -(chi-chi^H)/(2 i pi),
!>
!> so that a physical diagonal response with Im chi < 0 has positive loss.
!> `eta` is recorded as a numerical integration broadening only; this module
!> never interprets it as an intrinsic linewidth.
module tddft_dyson_mod
   use precision_mod, only: rp
   use math_mod, only: pi
   implicit none

   private

   type, public :: tddft_dyson_options
      logical :: diagonalize_loss = .false.
      logical :: diagonalize_xi = .false.
   end type tddft_dyson_options

   type, public :: tddft_dyson_metadata
      character(len=16) :: energy_unit = 'Rydberg'
      character(len=32) :: susceptibility_unit = '1/Rydberg'
      character(len=112) :: dyson_convention = 'solve (I-chi_KS K) chi = chi_KS with LAPACK zgesv; no explicit inverse'
      character(len=112) :: loss_convention = 'loss = -(chi-chi^H)/(2 i pi); positive diagonal loss denotes positive-frequency absorption'
      character(len=112) :: eta_convention = 'eta is numerical broadening in Ry, distinct from any fitted or extrapolated intrinsic linewidth'
      real(rp) :: eta = 0.0_rp
      real(rp) :: solve_cpu_seconds = 0.0_rp
      real(rp) :: diagonalization_cpu_seconds = 0.0_rp
   end type tddft_dyson_metadata

   !> Arrays are (response left, response right, frequency), except that the
   !> diagonalization arrays are (mode, frequency) and (response, mode,
   !> frequency).  The original KS input is retained intentionally so callers
   !> may stream one q batch and still write all three response levels.
   type, public :: tddft_dyson_result
      complex(rp), allocatable :: chi_ks(:, :, :)
      complex(rp), allocatable :: chi_ks_loss(:, :, :)
      complex(rp), allocatable :: xi(:, :, :)
      complex(rp), allocatable :: chi(:, :, :)
      complex(rp), allocatable :: loss(:, :, :)
      real(rp), allocatable :: site_spectral_weight(:, :)
      real(rp), allocatable :: trace_spectral_weight(:)
      real(rp), allocatable :: chi_ks_site_spectral_weight(:, :)
      real(rp), allocatable :: chi_ks_trace_spectral_weight(:)
      integer, allocatable :: solve_info(:)
      real(rp), allocatable :: loss_hermiticity_residual(:)
      real(rp), allocatable :: loss_eigenvalues(:, :)
      complex(rp), allocatable :: loss_eigenvectors(:, :, :)
      complex(rp), allocatable :: xi_eigenvalues(:, :)
      complex(rp), allocatable :: xi_eigenvectors(:, :, :)
      type(tddft_dyson_metadata) :: metadata
   end type tddft_dyson_result

   public :: enhance_tddft_susceptibility
   public :: enhance_tddft_susceptibility_from_xi
   public :: solve_tddft_dyson_frequency
   public :: tddft_loss_matrix
   public :: loss_matrix_hermiticity_residual
   public :: diagonalize_nonhermitian_response
   public :: write_tddft_dyson_text

contains

   !> Enhance one q batch.  Each frequency is a multi-right-hand-side complex
   !> solve, which avoids forming [I-Xi]^-1 in production.
   subroutine enhance_tddft_susceptibility(chi_ks, kernel, eta, options, result)
      complex(rp), intent(in) :: chi_ks(:, :, :), kernel(:, :)
      real(rp), intent(in) :: eta
      type(tddft_dyson_options), intent(in) :: options
      type(tddft_dyson_result), intent(out) :: result
      integer :: n, nw, iw, i, info
      real(rp) :: t_start, t_stop

      n = size(chi_ks, 1)
      nw = size(chi_ks, 3)
      if (n <= 0 .or. size(chi_ks, 2) /= n .or. nw <= 0) then
         error stop 'enhance_tddft_susceptibility: chi_KS must be a nonempty square frequency batch'
      end if
      if (any(shape(kernel) /= [n, n])) error stop 'enhance_tddft_susceptibility: kernel dimension mismatch'
      if (eta <= 0.0_rp) error stop 'enhance_tddft_susceptibility: eta must be positive and numerical'

      allocate(result%chi_ks(n, n, nw), result%chi_ks_loss(n, n, nw), result%xi(n, n, nw), result%chi(n, n, nw), &
         result%loss(n, n, nw), result%site_spectral_weight(n, nw), result%trace_spectral_weight(nw), &
         result%chi_ks_site_spectral_weight(n, nw), result%chi_ks_trace_spectral_weight(nw), result%solve_info(nw), &
         result%loss_hermiticity_residual(nw))
      result%chi_ks = chi_ks
      do iw = 1, nw
         result%chi_ks_loss(:, :, iw) = tddft_loss_matrix(chi_ks(:, :, iw))
         result%loss_hermiticity_residual(iw) = loss_matrix_hermiticity_residual(result%chi_ks_loss(:, :, iw))
         result%chi_ks_site_spectral_weight(:, iw) = real([(result%chi_ks_loss(i, i, iw), i=1, n)], rp)
         result%chi_ks_trace_spectral_weight(iw) = sum(result%chi_ks_site_spectral_weight(:, iw))
         call cpu_time(t_start)
         call solve_tddft_dyson_frequency(chi_ks(:, :, iw), kernel, result%chi(:, :, iw), result%xi(:, :, iw), info)
         call cpu_time(t_stop)
         result%metadata%solve_cpu_seconds = result%metadata%solve_cpu_seconds + t_stop-t_start
         result%solve_info(iw) = info
         if (info /= 0) error stop 'enhance_tddft_susceptibility: LAPACK zgesv failed; I-chi_KS K is singular'
         result%loss(:, :, iw) = tddft_loss_matrix(result%chi(:, :, iw))
         result%loss_hermiticity_residual(iw) = max(result%loss_hermiticity_residual(iw), &
            loss_matrix_hermiticity_residual(result%loss(:, :, iw)))
         result%site_spectral_weight(:, iw) = real([(result%loss(i, i, iw), i=1, size(result%loss, 1))], rp)
         result%trace_spectral_weight(iw) = sum(result%site_spectral_weight(:, iw))
      end do
      result%metadata%eta = eta

      if (options%diagonalize_loss) then
         allocate(result%loss_eigenvalues(n, nw), result%loss_eigenvectors(n, n, nw))
         call cpu_time(t_start)
         do iw = 1, nw
            call diagonalize_hermitian_loss(result%loss(:, :, iw), result%loss_eigenvalues(:, iw), &
               result%loss_eigenvectors(:, :, iw))
         end do
         call cpu_time(t_stop)
         result%metadata%diagonalization_cpu_seconds = result%metadata%diagonalization_cpu_seconds + t_stop-t_start
      end if
      if (options%diagonalize_xi) then
         allocate(result%xi_eigenvalues(n, nw), result%xi_eigenvectors(n, n, nw))
         call cpu_time(t_start)
         do iw = 1, nw
            call diagonalize_nonhermitian_response(result%xi(:, :, iw), result%xi_eigenvalues(:, iw), &
               result%xi_eigenvectors(:, :, iw))
         end do
         call cpu_time(t_stop)
         result%metadata%diagonalization_cpu_seconds = result%metadata%diagonalization_cpu_seconds + t_stop-t_start
      end if
   end subroutine enhance_tddft_susceptibility

   !> Pair-potential/self-enhancement route.  Xi is already the ordered
   !> kernel-weighted Lehmann sum, so this intentionally has no response-space
   !> K and cannot reconstruct one by inversion.
   subroutine enhance_tddft_susceptibility_from_xi(chi_ks, xi, eta, options, result)
      complex(rp), intent(in) :: chi_ks(:, :, :), xi(:, :, :)
      real(rp), intent(in) :: eta
      type(tddft_dyson_options), intent(in) :: options
      type(tddft_dyson_result), intent(out) :: result
      integer :: n, nw, iw, i, info
      real(rp) :: t_start, t_stop

      n = size(chi_ks, 1); nw = size(chi_ks, 3)
      if (n <= 0 .or. size(chi_ks, 2) /= n .or. nw <= 0 .or. any(shape(xi) /= shape(chi_ks))) then
         error stop 'enhance_tddft_susceptibility_from_xi: incompatible chi_KS/Xi frequency batches'
      end if
      if (eta <= 0.0_rp) error stop 'enhance_tddft_susceptibility_from_xi: eta must be positive and numerical'
      allocate(result%chi_ks(n, n, nw), result%chi_ks_loss(n, n, nw), result%xi(n, n, nw), result%chi(n, n, nw), &
         result%loss(n, n, nw), result%site_spectral_weight(n, nw), result%trace_spectral_weight(nw), &
         result%chi_ks_site_spectral_weight(n, nw), result%chi_ks_trace_spectral_weight(nw), result%solve_info(nw), &
         result%loss_hermiticity_residual(nw))
      result%chi_ks = chi_ks; result%xi = xi
      do iw = 1, nw
         result%chi_ks_loss(:, :, iw) = tddft_loss_matrix(chi_ks(:, :, iw))
         result%loss_hermiticity_residual(iw) = loss_matrix_hermiticity_residual(result%chi_ks_loss(:, :, iw))
         result%chi_ks_site_spectral_weight(:, iw) = real([(result%chi_ks_loss(i, i, iw), i=1, n)], rp)
         result%chi_ks_trace_spectral_weight(iw) = sum(result%chi_ks_site_spectral_weight(:, iw))
         call cpu_time(t_start)
         call solve_tddft_dyson_direct_xi_frequency(chi_ks(:, :, iw), xi(:, :, iw), result%chi(:, :, iw), info)
         call cpu_time(t_stop)
         result%metadata%solve_cpu_seconds = result%metadata%solve_cpu_seconds + t_stop-t_start
         result%solve_info(iw) = info
         if (info /= 0) error stop 'enhance_tddft_susceptibility_from_xi: LAPACK zgesv failed; I-Xi is singular'
         result%loss(:, :, iw) = tddft_loss_matrix(result%chi(:, :, iw))
         result%loss_hermiticity_residual(iw) = max(result%loss_hermiticity_residual(iw), &
            loss_matrix_hermiticity_residual(result%loss(:, :, iw)))
         result%site_spectral_weight(:, iw) = real([(result%loss(i, i, iw), i=1, n)], rp)
         result%trace_spectral_weight(iw) = sum(result%site_spectral_weight(:, iw))
      end do
      result%metadata%eta = eta
      if (options%diagonalize_loss) then
         allocate(result%loss_eigenvalues(n, nw), result%loss_eigenvectors(n, n, nw))
         call cpu_time(t_start)
         do iw = 1, nw
            call diagonalize_hermitian_loss(result%loss(:, :, iw), result%loss_eigenvalues(:, iw), &
               result%loss_eigenvectors(:, :, iw))
         end do
         call cpu_time(t_stop)
         result%metadata%diagonalization_cpu_seconds = result%metadata%diagonalization_cpu_seconds + t_stop-t_start
      end if
      if (options%diagonalize_xi) then
         allocate(result%xi_eigenvalues(n, nw), result%xi_eigenvectors(n, n, nw))
         call cpu_time(t_start)
         do iw = 1, nw
            call diagonalize_nonhermitian_response(result%xi(:, :, iw), result%xi_eigenvalues(:, iw), &
               result%xi_eigenvectors(:, :, iw))
         end do
         call cpu_time(t_stop)
         result%metadata%diagonalization_cpu_seconds = result%metadata%diagonalization_cpu_seconds + t_stop-t_start
      end if
   end subroutine enhance_tddft_susceptibility_from_xi

   !> Streaming primitive for one (q,omega): factor A=I-chi_KS K and solve
   !> A X=chi_KS.  `info` is returned so an outer production scheduler can
   !> report/skip a failed point rather than relying on a matrix inverse.
   subroutine solve_tddft_dyson_frequency(chi_ks, kernel, chi, xi, info)
      complex(rp), intent(in) :: chi_ks(:, :), kernel(:, :)
      complex(rp), intent(out) :: chi(:, :), xi(:, :)
      integer, intent(out) :: info
      complex(rp), allocatable :: system(:, :)
      integer, allocatable :: ipiv(:)
      integer :: n, i

      n = size(chi_ks, 1)
      if (n <= 0 .or. size(chi_ks, 2) /= n .or. any(shape(kernel) /= [n, n]) .or. any(shape(chi) /= [n, n]) .or. &
          any(shape(xi) /= [n, n])) error stop 'solve_tddft_dyson_frequency: incompatible response dimensions'
      allocate(system(n, n), ipiv(n))
      xi = matmul(chi_ks, kernel)
      system = -xi
      do i = 1, n
         system(i, i) = system(i, i) + cmplx(1.0_rp, 0.0_rp, rp)
      end do
      chi = chi_ks
      call zgesv(n, n, system, n, ipiv, chi, n, info)
   end subroutine solve_tddft_dyson_frequency

   !> Solve (I-Xi) chi=chi_KS without ever referring to a truncated K.
   subroutine solve_tddft_dyson_direct_xi_frequency(chi_ks, xi, chi, info)
      complex(rp), intent(in) :: chi_ks(:, :), xi(:, :)
      complex(rp), intent(out) :: chi(:, :)
      integer, intent(out) :: info
      complex(rp), allocatable :: system(:, :)
      integer, allocatable :: ipiv(:)
      integer :: n, i

      n = size(chi_ks, 1)
      if (n <= 0 .or. size(chi_ks, 2) /= n .or. any(shape(xi) /= [n, n]) .or. any(shape(chi) /= [n, n])) then
         error stop 'solve_tddft_dyson_direct_xi_frequency: incompatible response dimensions'
      end if
      allocate(system(n, n), ipiv(n)); system = -xi
      do i = 1, n
         system(i, i) = system(i, i) + cmplx(1.0_rp, 0.0_rp, rp)
      end do
      chi = chi_ks
      call zgesv(n, n, system, n, ipiv, chi, n, info)
   end subroutine solve_tddft_dyson_direct_xi_frequency

   !> Hermitian loss matrix.  Use this rather than elementwise -Im chi when
   !> response channels are coupled: off-diagonal entries must be Hermitian.
   function tddft_loss_matrix(chi) result(loss)
      complex(rp), intent(in) :: chi(:, :)
      complex(rp) :: loss(size(chi, 1), size(chi, 2))
      integer :: i, j, n

      n = size(chi, 1)
      if (size(chi, 2) /= n) error stop 'tddft_loss_matrix: chi must be square'
      do j = 1, n
         do i = 1, n
            loss(i, j) = -(chi(i, j) - conjg(chi(j, i)))/cmplx(0.0_rp, 2.0_rp*pi, rp)
         end do
      end do
   end function tddft_loss_matrix

   !> Maximum anti-Hermitian residual of a loss matrix.  The loss matrix is
   !> formed from both chi(i,j) and chi(j,i), so this is a useful runtime
   !> invariant rather than an elementwise -Im diagnostic.
   pure real(rp) function loss_matrix_hermiticity_residual(loss) result(residual)
      complex(rp), intent(in) :: loss(:, :)
      integer :: n

      n = size(loss, 1)
      if (size(loss, 2) /= n) error stop 'loss_matrix_hermiticity_residual: loss must be square'
      residual = 0.0_rp
      if (n > 0) residual = maxval(abs(loss-conjg(transpose(loss))))
   end function loss_matrix_hermiticity_residual

   subroutine diagonalize_hermitian_loss(loss, eigenvalues, eigenvectors)
      complex(rp), intent(in) :: loss(:, :)
      real(rp), intent(out) :: eigenvalues(:)
      complex(rp), intent(out) :: eigenvectors(:, :)
      complex(rp), allocatable :: work(:)
      real(rp), allocatable :: rwork(:)
      complex(rp) :: work_query(1)
      integer :: n, info, lwork

      n = size(loss, 1)
      if (size(loss, 2) /= n .or. size(eigenvalues) /= n .or. any(shape(eigenvectors) /= [n, n])) then
         error stop 'diagonalize_hermitian_loss: incompatible dimensions'
      end if
      eigenvectors = loss
      allocate(rwork(max(1, 3*n - 2)))
      call zheev('V', 'U', n, eigenvectors, n, eigenvalues, work_query, -1, rwork, info)
      if (info /= 0) error stop 'diagonalize_hermitian_loss: LAPACK workspace query failed'
      lwork = max(1, int(real(work_query(1), rp)))
      allocate(work(lwork))
      eigenvectors = loss
      call zheev('V', 'U', n, eigenvectors, n, eigenvalues, work, lwork, rwork, info)
      if (info /= 0) error stop 'diagonalize_hermitian_loss: LAPACK zheev failed'
   end subroutine diagonalize_hermitian_loss

   !> Right eigenvectors of a generally non-Hermitian Xi.  Columns match the
   !> returned eigenvalues and are normalized by LAPACK.
   subroutine diagonalize_nonhermitian_response(matrix, eigenvalues, eigenvectors)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp), intent(out) :: eigenvalues(:), eigenvectors(:, :)
      complex(rp), allocatable :: work_matrix(:, :), work(:), left_vectors(:, :)
      real(rp), allocatable :: rwork(:)
      complex(rp) :: work_query(1)
      integer :: n, info, lwork

      n = size(matrix, 1)
      if (size(matrix, 2) /= n .or. size(eigenvalues) /= n .or. any(shape(eigenvectors) /= [n, n])) then
         error stop 'diagonalize_nonhermitian_response: incompatible dimensions'
      end if
      allocate(work_matrix(n, n), left_vectors(1, 1), rwork(max(1, 2*n)))
      work_matrix = matrix
      call zgeev('N', 'V', n, work_matrix, n, eigenvalues, left_vectors, 1, eigenvectors, n, work_query, -1, rwork, info)
      if (info /= 0) error stop 'diagonalize_nonhermitian_response: LAPACK workspace query failed'
      lwork = max(1, int(real(work_query(1), rp)))
      allocate(work(lwork))
      work_matrix = matrix
      call zgeev('N', 'V', n, work_matrix, n, eigenvalues, left_vectors, 1, eigenvectors, n, work, lwork, rwork, info)
      if (info /= 0) error stop 'diagonalize_nonhermitian_response: LAPACK zgeev failed'
   end subroutine diagonalize_nonhermitian_response

   !> Plain text is deliberately supplementary to the in-memory API.  It
   !> contains bare, self-enhancement, enhanced, and loss data for one q batch.
   subroutine write_tddft_dyson_text(filename, omega, result)
      character(len=*), intent(in) :: filename
      real(rp), intent(in) :: omega(:)
      type(tddft_dyson_result), intent(in) :: result
      integer :: unit, ios, iw, i, j

      if (.not. allocated(result%chi) .or. size(omega) /= size(result%chi, 3)) then
         error stop 'write_tddft_dyson_text: omega/result shape mismatch'
      end if
      open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
      if (ios /= 0) error stop 'write_tddft_dyson_text: cannot open output file'
      write(unit, '(a,a)') '# dyson_convention = ', trim(result%metadata%dyson_convention)
      write(unit, '(a,a)') '# loss_convention = ', trim(result%metadata%loss_convention)
      write(unit, '(a,a)') '# eta_convention = ', trim(result%metadata%eta_convention)
      write(unit, '(a)') '# stoner_reference = chi_KS loss spectrum; enhanced loss is the dressed response'
      write(unit, '(a,es24.16)') '# eta_Ry = ', result%metadata%eta
      write(unit, '(a,es24.16)') '# profile_dyson_solves_cpu_s = ', result%metadata%solve_cpu_seconds
      write(unit, '(a,es24.16)') '# profile_dyson_diagonalizations_cpu_s = ', result%metadata%diagonalization_cpu_seconds
      write(unit, '(a)') '# record omega_Ry i j Re Im'
      if (allocated(result%loss_hermiticity_residual)) then
         write(unit, '(a)') '# loss_hermiticity_residual is max|L-L^H| at each frequency'
      end if
      do iw = 1, size(omega)
         do j = 1, size(result%chi, 2)
            do i = 1, size(result%chi, 1)
               write(unit, '(a,1x,es24.16,2(1x,i0),2(1x,es24.16))') 'chi_KS', omega(iw), i, j, &
                  real(result%chi_ks(i, j, iw), rp), aimag(result%chi_ks(i, j, iw))
               write(unit, '(a,1x,es24.16,2(1x,i0),2(1x,es24.16))') 'Xi', omega(iw), i, j, &
                  real(result%xi(i, j, iw), rp), aimag(result%xi(i, j, iw))
               write(unit, '(a,1x,es24.16,2(1x,i0),2(1x,es24.16))') 'chi', omega(iw), i, j, &
                  real(result%chi(i, j, iw), rp), aimag(result%chi(i, j, iw))
               write(unit, '(a,1x,es24.16,2(1x,i0),2(1x,es24.16))') 'loss', omega(iw), i, j, &
                  real(result%loss(i, j, iw), rp), aimag(result%loss(i, j, iw))
            end do
         end do
         do i = 1, size(result%site_spectral_weight, 1)
            write(unit, '(a,1x,es24.16,1x,i0,1x,es24.16)') 'chi_KS_site_loss', omega(iw), i, &
               result%chi_ks_site_spectral_weight(i, iw)
            write(unit, '(a,1x,es24.16,1x,i0,1x,es24.16)') 'site_loss', omega(iw), i, result%site_spectral_weight(i, iw)
         end do
         write(unit, '(a,1x,es24.16,1x,es24.16)') 'chi_KS_trace_loss', omega(iw), result%chi_ks_trace_spectral_weight(iw)
         write(unit, '(a,1x,es24.16,1x,es24.16)') 'trace_loss', omega(iw), result%trace_spectral_weight(iw)
         if (allocated(result%loss_hermiticity_residual)) then
            write(unit, '(a,1x,es24.16,1x,es24.16)') 'loss_hermiticity_residual', omega(iw), &
               result%loss_hermiticity_residual(iw)
         end if
         if (allocated(result%loss_eigenvalues)) then
            do i = 1, size(result%loss_eigenvalues, 1)
               write(unit, '(a,1x,es24.16,1x,i0,1x,es24.16)') 'loss_mode', omega(iw), i, result%loss_eigenvalues(i, iw)
               do j = 1, size(result%loss_eigenvectors, 1)
                  write(unit, '(a,1x,es24.16,2(1x,i0),2(1x,es24.16))') 'loss_mode_vector', omega(iw), i, j, &
                     real(result%loss_eigenvectors(j, i, iw), rp), aimag(result%loss_eigenvectors(j, i, iw))
               end do
            end do
         end if
      end do
      close(unit)
   end subroutine write_tddft_dyson_text

end module tddft_dyson_mod
