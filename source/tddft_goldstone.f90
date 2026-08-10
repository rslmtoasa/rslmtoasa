!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Site-projected transverse ALSDA kernel and Goldstone diagnostics.
!>
!> `K_perp` is deliberately obtained only from xc_response_kernel_mod.  In
!> particular, this module never inspects hamiltonian%hxc, cx1, or another
!> assembled Hamiltonian quantity: those are not a documented TDDFT kernel.
!>
!> The optional sum-rule repair is a *static, site-diagonal* correction.  For
!> m_i /= 0 it chooses K_i^SR such that
!>
!>   chi_KS(0,0) diag(K^SR) m = m,
!>
!> i.e. K_i^SR = [chi_KS(0,0)^(-1) m]_i / m_i.  It is only permitted for a
!> collinear calculation without SOC or an external symmetry-breaking field.
!> Raw Xi and its diagnostics are always retained in the result.
module tddft_goldstone_mod
   use precision_mod, only: rp
   use xc_response_kernel_mod, only: xc_response_kernel_provider
   implicit none

   private

   integer, parameter, public :: GOLDSTONE_OFF = 0
   integer, parameter, public :: GOLDSTONE_DIAGNOSE = 1
   integer, parameter, public :: GOLDSTONE_SUM_RULE = 2

   type, public :: tddft_goldstone_options
      ! Development default: report the raw Goldstone quality, never repair it.
      character(len=16) :: goldstone_mode = 'diagnose'
      logical :: has_soc = .false.
      logical :: has_external_field = .false.
   end type tddft_goldstone_options

   type, public :: tddft_goldstone_diagnostics
      logical :: available = .false.
      logical :: has_bare_spectral_gap = .false.
      real(rp) :: bare_spectral_gap = -1.0_rp
      real(rp) :: residual = -1.0_rp
      real(rp) :: closest_eigenvalue_distance = -1.0_rp
      real(rp) :: magnetization_overlap = -1.0_rp
      integer :: closest_eigenvalue_index = 0
      complex(rp) :: closest_eigenvalue = cmplx(0.0_rp, 0.0_rp, rp)
      complex(rp), allocatable :: eigenvalues(:)
      complex(rp), allocatable :: eigenvectors(:, :)
      complex(rp), allocatable :: closest_eigenvector(:)
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
   end type tddft_goldstone_result

   public :: build_site_projected_k_perp
   public :: construct_transverse_xi
   public :: evaluate_goldstone
   public :: evaluate_raw_xi_diagnostics
   public :: write_goldstone_diagnostics_text

contains

   !> Return the site-projected transverse ALSDA kernel recorded by the XC
   !> response provider.  `has_k_perp` makes missing functional provenance a
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
         if (.not. provider%site(isite)%has_k_perp) then
            error stop 'build_site_projected_k_perp: site has no xc_response_kernel K_perp'
         end if
         k_perp(isite) = cmplx(provider%site(isite)%k_perp, 0.0_rp, rp)
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
   !> caller's independently determined lowest bare spin-flip transition; it
   !> is not inferred from one complex chi_KS matrix sample.
   subroutine evaluate_goldstone(chi_ks_static, provider, options, result, bare_spectral_gap)
      complex(rp), intent(in) :: chi_ks_static(:, :)
      type(xc_response_kernel_provider), intent(in) :: provider
      type(tddft_goldstone_options), intent(in) :: options
      type(tddft_goldstone_result), intent(out) :: result
      real(rp), intent(in), optional :: bare_spectral_gap
      complex(rp), allocatable :: magnetization(:), rhs(:)
      integer :: mode

      call build_site_projected_k_perp(provider, result%k_perp)
      call require_square_response(chi_ks_static, result%k_perp, 'evaluate_goldstone')
      allocate(magnetization(size(result%k_perp)))
      magnetization = cmplx(provider%site(:)%spin_population, 0.0_rp, rp)
      if (sqrt(sum(abs(magnetization)**2)) <= tiny(1.0_rp)) then
         error stop 'evaluate_goldstone: projected ground-state magnetization is zero'
      end if

      result%xi_raw = construct_transverse_xi(chi_ks_static, result%k_perp)
      mode = goldstone_mode_code(options%goldstone_mode)
      result%sum_rule_requested = mode == GOLDSTONE_SUM_RULE
      if (mode /= GOLDSTONE_OFF) then
         call calculate_diagnostics(result%xi_raw, magnetization, result%raw, bare_spectral_gap)
      end if

      if (mode /= GOLDSTONE_SUM_RULE) return
      if (options%has_soc .or. options%has_external_field) then
         ! A zero-frequency transverse mode is not required if symmetry is
         ! broken physically.  Preserve raw diagnostics and disable repair.
         result%sum_rule_disabled_by_symmetry_breaking = .true.
         return
      end if
      if (any(abs(magnetization) <= tiny(1.0_rp))) then
         error stop 'evaluate_goldstone: sum-rule correction requires nonzero magnetization on every response site'
      end if

      allocate(rhs(size(magnetization)))
      rhs = magnetization
      call solve_linear_system(chi_ks_static, rhs)
      allocate(result%k_perp_sum_rule(size(rhs)))
      result%k_perp_sum_rule = rhs/magnetization
      result%xi_corrected = construct_transverse_xi(chi_ks_static, result%k_perp_sum_rule)
      call calculate_diagnostics(result%xi_corrected, magnetization, result%corrected, bare_spectral_gap)
      result%sum_rule_applied = .true.
   end subroutine evaluate_goldstone

   !> Diagnose an already assembled self-enhancement operator.  This keeps the
   !> pair-potential route out of the legacy site-scalar K interface while
   !> retaining exactly the same raw residual/eigenmode definitions.
   subroutine evaluate_raw_xi_diagnostics(xi, magnetization, diagnostics, bare_spectral_gap)
      complex(rp), intent(in) :: xi(:, :), magnetization(:)
      type(tddft_goldstone_diagnostics), intent(out) :: diagnostics
      real(rp), intent(in), optional :: bare_spectral_gap

      if (size(xi, 1) /= size(xi, 2) .or. size(magnetization) /= size(xi, 1)) then
         error stop 'evaluate_raw_xi_diagnostics: Xi and magnetization dimensions are incompatible'
      end if
      if (sqrt(sum(abs(magnetization)**2)) <= tiny(1.0_rp)) then
         error stop 'evaluate_raw_xi_diagnostics: magnetization is zero'
      end if
      if (present(bare_spectral_gap)) then
         call calculate_diagnostics(xi, magnetization, diagnostics, bare_spectral_gap)
      else
         call calculate_diagnostics(xi, magnetization, diagnostics)
      end if
   end subroutine evaluate_raw_xi_diagnostics

   !> Write raw diagnostics first.  Corrected diagnostics are an additional
   !> record, never a replacement, so output remains useful for convergence
   !> studies even when sum-rule mode was requested.
   subroutine write_goldstone_diagnostics_text(filename, result)
      character(len=*), intent(in) :: filename
      type(tddft_goldstone_result), intent(in) :: result
      integer :: unit, ios, i

      open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
      if (ios /= 0) error stop 'write_goldstone_diagnostics_text: cannot open output file'
      write(unit, '(a)') '# Xi = chi_KS K_perp; K_perp provenance = xc_response_kernel'
      write(unit, '(a)') '# raw Goldstone diagnostics are retained when sum_rule is active'
      write(unit, '(a,l1)') '# sum_rule_requested = ', result%sum_rule_requested
      write(unit, '(a,l1)') '# sum_rule_applied = ', result%sum_rule_applied
      write(unit, '(a,l1)') '# sum_rule_disabled_by_SOC_or_external_field = ', result%sum_rule_disabled_by_symmetry_breaking
      call write_one_diagnostics(unit, 'raw', result%raw)
      if (result%sum_rule_applied) call write_one_diagnostics(unit, 'sum_rule_corrected', result%corrected)
      if (allocated(result%k_perp)) then
         do i = 1, size(result%k_perp)
            write(unit, '(a,1x,i0,2(1x,es24.16))') 'kernel_raw', i, real(result%k_perp(i), rp), aimag(result%k_perp(i))
         end do
      end if
      if (allocated(result%k_perp_sum_rule)) then
         do i = 1, size(result%k_perp_sum_rule)
            write(unit, '(a,1x,i0,2(1x,es24.16))') 'kernel_sum_rule', i, real(result%k_perp_sum_rule(i), rp), &
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
      real(rp) :: norm_m
      integer :: i

      call diagonalize_nonhermitian(xi, diagnostics%eigenvalues, diagnostics%eigenvectors)
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
      norm_m = sqrt(sum(abs(magnetization)**2))
      diagnostics%magnetization_overlap = abs(dot_product(magnetization, diagnostics%closest_eigenvector))/ &
         (norm_m*sqrt(sum(abs(diagnostics%closest_eigenvector)**2)))
      allocate(residual_vector(size(magnetization)))
      residual_vector = matmul(xi, magnetization) - magnetization
      diagnostics%residual = sqrt(sum(abs(residual_vector)**2))/norm_m
      diagnostics%available = .true.
      if (present(bare_spectral_gap)) then
         diagnostics%has_bare_spectral_gap = .true.
         diagnostics%bare_spectral_gap = bare_spectral_gap
      end if
   end subroutine calculate_diagnostics

   subroutine diagonalize_nonhermitian(matrix, eigenvalues, eigenvectors)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp), allocatable, intent(out) :: eigenvalues(:), eigenvectors(:, :)
      complex(rp), allocatable :: work_matrix(:, :), work(:), vl(:, :)
      real(rp), allocatable :: rwork(:)
      complex(rp) :: work_query(1)
      integer :: n, info, lwork

      n = size(matrix, 1)
      if (size(matrix, 2) /= n) error stop 'diagonalize_nonhermitian: matrix must be square'
      allocate(work_matrix(n, n), eigenvalues(n), eigenvectors(n, n), vl(1, 1), rwork(max(1, 2*n)))
      work_matrix = matrix
      call zgeev('N', 'V', n, work_matrix, n, eigenvalues, vl, 1, eigenvectors, n, work_query, -1, rwork, info)
      if (info /= 0) error stop 'diagonalize_nonhermitian: LAPACK workspace query failed'
      lwork = max(1, int(real(work_query(1), rp)))
      allocate(work(lwork))
      work_matrix = matrix
      call zgeev('N', 'V', n, work_matrix, n, eigenvalues, vl, 1, eigenvectors, n, work, lwork, rwork, info)
      if (info /= 0) error stop 'diagonalize_nonhermitian: LAPACK zgeev failed'
   end subroutine diagonalize_nonhermitian

   subroutine solve_linear_system(matrix, rhs)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp), intent(inout) :: rhs(:)
      complex(rp), allocatable :: work_matrix(:, :)
      integer, allocatable :: ipiv(:)
      integer :: n, info

      n = size(rhs)
      if (size(matrix, 1) /= n .or. size(matrix, 2) /= n) then
         error stop 'solve_linear_system: matrix/vector dimensions are incompatible'
      end if
      allocate(work_matrix(n, n), ipiv(n))
      work_matrix = matrix
      call zgesv(n, 1, work_matrix, n, ipiv, rhs, n, info)
      if (info /= 0) error stop 'evaluate_goldstone: sum-rule chi_KS(0,0) is singular'
   end subroutine solve_linear_system

   integer function goldstone_mode_code(mode) result(code)
      character(len=*), intent(in) :: mode

      select case (trim(adjustl(mode)))
      case ('off')
         code = GOLDSTONE_OFF
      case ('diagnose')
         code = GOLDSTONE_DIAGNOSE
      case ('sum_rule')
         code = GOLDSTONE_SUM_RULE
      case default
         error stop 'evaluate_goldstone: goldstone_mode must be off, diagnose, or sum_rule'
      end select
   end function goldstone_mode_code

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
      write(unit, '(a,a,a,es24.16)') trim(prefix), '_residual', '', diagnostics%residual
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
   end subroutine write_one_diagnostics

end module tddft_goldstone_mod
