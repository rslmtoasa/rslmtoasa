!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
! PROGRAM: test_kspace_gf_validation
!
!> @brief Independent validation of the orthogonal K-space Lehmann resolvent.
!> @details This test deliberately uses only the pure Green-function kernels
!>          and deterministic finite Hermitian fixtures.  It does not build a
!>          TD-DFT response or depend on an LMTO ground-state calculation.
!>
!>          The direct-inverse oracle is [z I - H(k)]^(-1), while the tested
!>          Lehmann object is sum_n |psi_n><psi_n|/(z-e_n).  The fixture also
!>          pins retarded/advanced conjugation, spectral positivity and
!>          normalization, large-|z| behavior, collinear spin blocks, and the
!>          local spectral function used by the existing BSF/DOS primitive.
!>
!>          The final profile is informational.  It exercises the production
!>          pair-block batch over complex energies and records caller-owned
!>          output/scratch sizes without making a machine-dependent timing
!>          claim.
!------------------------------------------------------------------------------
program test_kspace_gf_validation
   use precision_mod, only: rp
   use math_mod, only: pi, i_unit
   use bsf_kernel_mod, only: bsf_spectral_trace
   use dyson_kernel_mod, only: dyson_kspace_inverse
   use lehmann_kernel_mod, only: lehmann_kspace_resolvent, lehmann_pair_block, pauli_decompose_block
   implicit none

   integer, parameter :: nmat = 4
   integer, parameter :: nk = 5
   real(rp), parameter :: tol_inverse = 2.0e-12_rp
   real(rp), parameter :: tol_identity = 5.0e-7_rp
   real(rp), parameter :: tol_spectral = 2.0e-12_rp
   real(rp), parameter :: tol_weight = 3.0e-3_rp
   logical :: failed

   failed = .false.
   call test_lehmann_matches_direct_inverse()
   call test_retarded_advanced_and_asymptotics()
   call test_spin_blocks_and_limits()
   call test_dos_from_lehmann()
   call profile_batched_complex_energies()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   subroutine test_lehmann_matches_direct_inverse()
      real(rp) :: h_eval(nmat, nk)
      real(rp) :: kfrac(3, nk)
      complex(rp) :: h_k(nmat, nmat, nk), h_vec(nmat, nmat, nk)
      complex(rp) :: z_values(5), sigma0(nmat, nmat), g_leh(nmat, nmat), g_inv(nmat, nmat)
      real(rp) :: max_error
      integer :: ik, iz

      call build_fixture(h_k, h_eval, h_vec, kfrac)
      sigma0 = (0.0_rp, 0.0_rp)
      z_values = [cmplx(-0.73_rp, 0.17_rp, rp), cmplx(0.22_rp, 0.41_rp, rp), &
         cmplx(1.10_rp, 0.08_rp, rp), cmplx(-0.31_rp, -0.27_rp, rp), cmplx(0.66_rp, -0.12_rp, rp)]

      max_error = 0.0_rp
      do ik = 1, nk
         do iz = 1, size(z_values)
            call lehmann_kspace_resolvent(h_eval(:, ik), h_vec(:, :, ik), z_values(iz), g_leh)
            call dyson_kspace_inverse(h_k(:, :, ik), z_values(iz), sigma0, g_inv)
            max_error = max(max_error, maxval(abs(g_leh-g_inv)))
         end do
      end do
      write (*, '(a,es12.4)') 'Test 1 (Lehmann == direct inverse) max_err = ', max_error
      call check_true('Lehmann and direct inverse agree', max_error <= tol_inverse)
   end subroutine test_lehmann_matches_direct_inverse

   subroutine test_retarded_advanced_and_asymptotics()
      real(rp) :: h_eval(nmat, nk)
      real(rp) :: kfrac(3, nk)
      complex(rp) :: h_k(nmat, nmat, nk), h_vec(nmat, nmat, nk)
      complex(rp) :: g_ret(nmat, nmat), g_adv(nmat, nmat), spectral(nmat, nmat), g_high(nmat, nmat)
      complex(rp) :: identity(nmat, nmat)
      real(rp), parameter :: energy = 0.13_rp, eta = 0.17_rp
      real(rp) :: advanced_error, hermitian_error, min_diagonal, asymptotic_error
      integer :: i

      call build_fixture(h_k, h_eval, h_vec, kfrac)
      call lehmann_kspace_resolvent(h_eval(:, 3), h_vec(:, :, 3), cmplx(energy, eta, rp), g_ret)
      call lehmann_kspace_resolvent(h_eval(:, 3), h_vec(:, :, 3), cmplx(energy, -eta, rp), g_adv)
      spectral = i_unit*(g_ret-g_adv)
      advanced_error = maxval(abs(g_adv-transpose(conjg(g_ret))))
      hermitian_error = maxval(abs(spectral-transpose(conjg(spectral))))
      min_diagonal = minval(real([(spectral(i, i), i=1,nmat)]))

      identity = (0.0_rp, 0.0_rp)
      do i = 1, nmat
         identity(i, i) = (1.0_rp, 0.0_rp)
      end do
      call lehmann_kspace_resolvent(h_eval(:, 3), h_vec(:, :, 3), cmplx(0.0_rp, 1.0e10_rp, rp), g_high)
      asymptotic_error = maxval(abs(cmplx(0.0_rp, 1.0e10_rp, rp)*g_high-identity))

      write (*, '(a,es12.4)') 'Test 2 (retarded/advanced identity) max_err = ', advanced_error
      write (*, '(a,es12.4)') 'Test 2 (spectral Hermiticity)      max_err = ', hermitian_error
      write (*, '(a,es12.4)') 'Test 2 (spectral diagonal minimum) value   = ', min_diagonal
      write (*, '(a,es12.4)') 'Test 2 (large-|z| zG-I)            max_err = ', asymptotic_error
      call check_true('retarded/advanced conjugate identity', advanced_error <= tol_spectral)
      call check_true('spectral function is Hermitian', hermitian_error <= tol_spectral)
      call check_true('retarded spectral diagonal is non-negative', min_diagonal >= -tol_spectral)
      call check_true('large-|z| asymptotics', asymptotic_error <= tol_identity)
   end subroutine test_retarded_advanced_and_asymptotics

   subroutine test_spin_blocks_and_limits()
      real(rp) :: evals(4)
      complex(rp) :: evecs(4, 4), g2(4, 4), g3(4, 4, 1)
      complex(rp) :: gnmag(2, 2, 1), gx(2, 2, 1), gy(2, 2, 1), gz(2, 2, 1)
      real(rp) :: error_offdiag, error_nonmagnetic, error_transverse, magnetic_z
      complex(rp), parameter :: z = cmplx(0.19_rp, 0.23_rp, rp)

      call identity_matrix(evecs)

      ! Basis ordering is (orbital-up, orbital-up, orbital-down, orbital-down).
      ! The duplicate spectrum is the nonmagnetic limit.
      evals = [-0.80_rp, 0.35_rp, -0.80_rp, 0.35_rp]
      call lehmann_kspace_resolvent(evals, evecs, z, g2)
      g3(:, :, 1) = g2
      call pauli_decompose_block(g3, 2, gnmag, gz, gy, gx)
      error_offdiag = max(maxval(abs(g2(1:2, 3:4))), maxval(abs(g2(3:4, 1:2))))
      error_nonmagnetic = maxval(abs(g2(1:2, 1:2)-g2(3:4, 3:4)))
      error_transverse = max(maxval(abs(gx)), max(maxval(abs(gy)), maxval(abs(gz))))
      write (*, '(a,es12.4)') 'Test 3 (nonmagnetic spin off-diagonal) max_err = ', error_offdiag
      write (*, '(a,es12.4)') 'Test 3 (nonmagnetic up/down equality) max_err = ', error_nonmagnetic
      write (*, '(a,es12.4)') 'Test 3 (nonmagnetic x/y/z components) max_err = ', error_transverse
      call check_true('nonmagnetic spin blocks are diagonal', error_offdiag <= tol_spectral)
      call check_true('nonmagnetic up/down blocks agree', error_nonmagnetic <= tol_spectral)
      call check_true('nonmagnetic transverse and z components vanish', error_transverse <= tol_spectral)

      ! A collinear magnetic splitting keeps x/y zero but produces a finite z
      ! component, which is the spin convention consumed by TD-DFT vertices.
      evals = [-0.80_rp, 0.35_rp, -0.25_rp, 0.92_rp]
      call lehmann_kspace_resolvent(evals, evecs, z, g2)
      g3(:, :, 1) = g2
      call pauli_decompose_block(g3, 2, gnmag, gz, gy, gx)
      error_offdiag = max(maxval(abs(g2(1:2, 3:4))), maxval(abs(g2(3:4, 1:2))))
      error_transverse = max(maxval(abs(gx)), maxval(abs(gy)))
      magnetic_z = maxval(abs(gz))
      write (*, '(a,es12.4)') 'Test 3 (collinear spin off-diagonal) max_err = ', error_offdiag
      write (*, '(a,es12.4)') 'Test 3 (collinear x/y components)    max_err = ', error_transverse
      write (*, '(a,es12.4)') 'Test 3 (collinear z component)       value   = ', magnetic_z
      call check_true('collinear spin blocks remain diagonal', error_offdiag <= tol_spectral)
      call check_true('collinear x/y components vanish', error_transverse <= tol_spectral)
      call check_true('collinear z component is non-zero', magnetic_z > 1.0e-3_rp)
   end subroutine test_spin_blocks_and_limits

   subroutine test_dos_from_lehmann()
      real(rp) :: h_eval(nmat, nk)
      real(rp) :: kfrac(3, nk)
      complex(rp) :: h_k(nmat, nmat, nk), h_vec(nmat, nmat, nk), gk(nmat, nmat)
      real(rp), parameter :: eta = 0.08_rp, emin = -100.0_rp, emax = 100.0_rp
      integer, parameter :: ne = 10001
      real(rp) :: de, energy, a_total, a_up, ref_total, ref_up, weight, a_prev
      real(rp) :: max_total_error, max_local_error, weight_error
      integer :: ie, n, i

      call build_fixture(h_k, h_eval, h_vec, kfrac)
      de = (emax-emin)/real(ne-1, rp)
      max_total_error = 0.0_rp
      max_local_error = 0.0_rp
      weight = 0.0_rp
      a_prev = 0.0_rp
      do ie = 1, ne
         energy = emin + de*real(ie-1, rp)
         call lehmann_kspace_resolvent(h_eval(:, 2), h_vec(:, :, 2), cmplx(energy, eta, rp), gk)
         a_total = bsf_spectral_trace(gk, [1, 2, 3, 4])
         a_up = bsf_spectral_trace(gk, [1, 2])
         ref_total = 0.0_rp
         ref_up = 0.0_rp
         do n = 1, nmat
            ref_total = ref_total + eta/(pi*((energy-h_eval(n, 2))**2+eta**2))
            do i = 1, 2
               ref_up = ref_up + abs(h_vec(i, n, 2))**2*eta/(pi*((energy-h_eval(n, 2))**2+eta**2))
            end do
         end do
         max_total_error = max(max_total_error, abs(a_total-ref_total))
         max_local_error = max(max_local_error, abs(a_up-ref_up))
         if (ie > 1) weight = weight + 0.5_rp*de*(a_total+a_prev)
         a_prev = a_total
      end do
      weight_error = abs(weight-real(nmat, rp))
      write (*, '(a,es12.4)') 'Test 4 (total DOS vs eigenvalue path) max_err = ', max_total_error
      write (*, '(a,es12.4)') 'Test 4 (local DOS vs projected path) max_err = ', max_local_error
      write (*, '(a,es12.4)') 'Test 4 (integrated spectral weight)  error   = ', weight_error
      call check_true('total DOS matches eigenvalue path', max_total_error <= tol_spectral)
      call check_true('local DOS matches projected eigenvalue path', max_local_error <= tol_spectral)
      call check_true('spectral normalization', weight_error <= tol_weight)
   end subroutine test_dos_from_lehmann

   subroutine profile_batched_complex_energies()
      integer, parameter :: profile_nmat = 12, profile_nk = 48, profile_ne_max = 128
      real(rp) :: evals(profile_nmat, profile_nk), kfrac(3, profile_nk), t_start, t_stop
      complex(rp) :: evecs(profile_nmat, profile_nmat, profile_nk), z(profile_ne_max)
      complex(rp), allocatable :: gblk(:, :, :)
      real(rp) :: dr(3), checksum, elapsed
      integer, parameter :: ne_samples(4) = [1, 8, 32, 128]
      integer :: i, ik, ie, ne

      do ik = 1, profile_nk
         kfrac(:, ik) = 0.0_rp
         kfrac(1, ik) = real(ik-1, rp)/real(profile_nk, rp)
         evecs(:, :, ik) = (0.0_rp, 0.0_rp)
         do i = 1, profile_nmat
            evecs(i, i, ik) = (1.0_rp, 0.0_rp)
            evals(i, ik) = -1.2_rp + 0.18_rp*real(i, rp) + 0.04_rp*cos(2.0_rp*pi*kfrac(1, ik)+real(i, rp))
         end do
      end do
      allocate(gblk(profile_nmat, profile_nmat, profile_ne_max))
      do ie = 1, profile_ne_max
         z(ie) = cmplx(-1.5_rp+0.02_rp*real(ie-1, rp), 0.11_rp, rp)
      end do
      dr = 0.0_rp

      write (*, '(a)') 'PROFILE_KSPACE_GF allocation_policy=caller-owned-output inner-loop-allocations=0'
      do ie = 1, size(ne_samples)
         ne = ne_samples(ie)
         call cpu_time(t_start)
         call lehmann_pair_block(evals, evecs, kfrac, z(1:ne), dr, 0, 0, profile_nmat, gblk(:, :, 1:ne))
         call cpu_time(t_stop)
         elapsed = t_stop-t_start
         checksum = sum(abs(gblk(:, :, 1:ne)))
         write (*, '(a,i0,1x,a,es12.4,1x,a,i0,1x,a,i0,1x,a,es12.4)') 'PROFILE_KSPACE_GF ne=', ne, &
            'elapsed_seconds=', elapsed, 'output_bytes=', 16*profile_nmat*profile_nmat*ne, &
            'scratch_bytes=', 16*profile_nmat*profile_nmat, 'checksum=', checksum
      end do
   end subroutine profile_batched_complex_energies

   subroutine build_fixture(h_k, h_eval, h_vec, kfrac)
      complex(rp), intent(out) :: h_k(:, :, :), h_vec(:, :, :)
      real(rp), intent(out) :: h_eval(:, :), kfrac(:, :)
      integer :: ik
      real(rp) :: phase

      do ik = 1, size(h_eval, 2)
         phase = 0.41_rp*real(ik, rp)
         kfrac(:, ik) = [real(ik-1, rp)/real(size(h_eval, 2), rp), 0.0_rp, 0.0_rp]
         h_k(:, :, ik) = (0.0_rp, 0.0_rp)
         h_k(1, 1, ik) = cmplx(-0.85_rp+0.03_rp*real(ik, rp), 0.0_rp, rp)
         h_k(2, 2, ik) = cmplx(-0.18_rp-0.02_rp*real(ik, rp), 0.0_rp, rp)
         h_k(3, 3, ik) = cmplx(0.31_rp+0.01_rp*real(ik, rp), 0.0_rp, rp)
         h_k(4, 4, ik) = cmplx(0.91_rp-0.025_rp*real(ik, rp), 0.0_rp, rp)
         h_k(1, 2, ik) = 0.16_rp*cmplx(cos(phase), sin(phase), rp)
         h_k(2, 1, ik) = conjg(h_k(1, 2, ik))
         h_k(1, 3, ik) = 0.07_rp*cmplx(cos(0.7_rp*phase), -sin(0.7_rp*phase), rp)
         h_k(3, 1, ik) = conjg(h_k(1, 3, ik))
         h_k(2, 4, ik) = 0.11_rp*cmplx(cos(1.3_rp*phase), sin(1.3_rp*phase), rp)
         h_k(4, 2, ik) = conjg(h_k(2, 4, ik))
         h_k(3, 4, ik) = 0.13_rp*cmplx(cos(0.9_rp*phase), sin(0.9_rp*phase), rp)
         h_k(4, 3, ik) = conjg(h_k(3, 4, ik))
         call hermitian_eig(h_k(:, :, ik), size(h_eval, 1), h_eval(:, ik), h_vec(:, :, ik))
      end do
   end subroutine build_fixture

   subroutine identity_matrix(matrix)
      complex(rp), intent(out) :: matrix(:, :)
      integer :: i

      matrix = (0.0_rp, 0.0_rp)
      do i = 1, size(matrix, 1)
         matrix(i, i) = (1.0_rp, 0.0_rp)
      end do
   end subroutine identity_matrix

   subroutine hermitian_eig(matrix_in, n, eigenvalues, eigenvectors)
      integer, intent(in) :: n
      complex(rp), intent(in) :: matrix_in(n, n)
      real(rp), intent(out) :: eigenvalues(n)
      complex(rp), intent(out) :: eigenvectors(n, n)
      complex(rp) :: work_query(1)
      complex(rp), allocatable :: work(:)
      real(rp) :: rwork(max(1, 3*n-2))
      integer :: info, lwork

      eigenvectors = matrix_in
      call zheev('V', 'U', n, eigenvectors, n, eigenvalues, work_query, -1, rwork, info)
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      call zheev('V', 'U', n, eigenvectors, n, eigenvalues, work, lwork, rwork, info)
      deallocate(work)
      if (info /= 0) then
         write (*, '(a,i0)') 'hermitian_eig: zheev info = ', info
         failed = .true.
      end if
   end subroutine hermitian_eig

   subroutine check_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition

      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine check_true

end program test_kspace_gf_validation
