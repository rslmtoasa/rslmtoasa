! DRESP-06A-FINAL Ward-mode algebra fixture.
! The matrix is deliberately non-Hermitian.  One known near-zero right mode
! is exactly aligned with m; a negative case makes m span several modes.
program test_lr_ward_mode_analysis
   use precision_mod, only: rp
   use lr_ward_mode_analysis_mod, only: lr_ward_mode_analysis_result, analyze_lr_ward_mode
   implicit none

   integer, parameter :: n = 4
   complex(rp) :: denominator(n, n), magnetization(n)
   complex(rp) :: rigid_delta(n, n), svd_delta(n, n), corrected(n, n), ward(n), svd_vector(n)
   type(lr_ward_mode_analysis_result) :: result, negative_result
   real(rp) :: reconstruction_tolerance, mnorm, sigma0
   integer :: i, j, sidx

   denominator = cmplx(0.0_rp, 0.0_rp, rp)
   denominator(1, 1) = cmplx(1.0e-8_rp, 2.0e-8_rp, rp)
   denominator(2, 2) = cmplx(0.45_rp, 0.03_rp, rp)
   denominator(3, 3) = cmplx(0.90_rp, -0.11_rp, rp)
   denominator(4, 4) = cmplx(1.40_rp, 0.17_rp, rp)
   denominator(1, 2) = cmplx(0.37_rp, 0.12_rp, rp)
   denominator(2, 3) = cmplx(-0.21_rp, 0.08_rp, rp)
   denominator(1, 3) = cmplx(0.11_rp, -0.04_rp, rp)
   denominator(3, 4) = cmplx(0.26_rp, 0.05_rp, rp)

   magnetization = cmplx(0.0_rp, 0.0_rp, rp)
   magnetization(1) = cmplx(1.0_rp, 0.0_rp, rp)
   call analyze_lr_ward_mode(denominator, magnetization, result)
   reconstruction_tolerance = 2.0e-11_rp
   if (.not. result%succeeded .or. result%nearest_zero_index < 1) error stop 'Ward fixture: eigensolver failed'
   if (abs(result%eigenvalues(result%nearest_zero_index)) > 1.0e-7_rp) then
      error stop 'Ward fixture: near-zero eigenmode was not identified'
   end if
   if (result%right_overlap(result%nearest_zero_index) < 0.999999_rp .or. &
       result%biorthogonal_weight(result%nearest_zero_index) < 0.999999_rp) then
      error stop 'Ward fixture: rigid mode overlap failed'
   end if
   if (result%ward_modal_norm_diagnostic(result%nearest_zero_index) < 0.999999_rp .or. &
       result%eigen_reconstruction_residual > reconstruction_tolerance .or. &
       result%ward_reconstruction_residual > reconstruction_tolerance) then
      error stop 'Ward fixture: single-mode decomposition failed'
   end if
   if (result%singular_vector_overlap(minloc(result%singular_values, dim=1)) < 0.9_rp) then
      error stop 'Ward fixture: SVD overlap diagnostic failed'
   end if

   ! Exact rank-one rigid correction: Delta D = -|r><m| / <m|m>.
   mnorm = sqrt(real(dot_product(magnetization, magnetization), rp))
   ward = matmul(denominator, magnetization)
   rigid_delta = cmplx(0.0_rp, 0.0_rp, rp)
   do j = 1, n
      do i = 1, n
         rigid_delta(i, j) = -ward(i)*conjg(magnetization(j))/max(mnorm*mnorm, tiny(1.0_rp))
      end do
   end do
   corrected = denominator + rigid_delta
   if (sqrt(sum(abs(matmul(corrected, magnetization))**2))/mnorm > 2.0e-11_rp) then
      error stop 'Ward fixture: exact rigid correction failed'
   end if

   ! Exact SVD correction: Delta D = -sigma_0 |u_0><v_0|.
   sidx = minloc(result%singular_values, dim=1)
   sigma0 = result%singular_values(sidx)
   svd_vector = result%singular_right_vectors(:, sidx)
   svd_delta = cmplx(0.0_rp, 0.0_rp, rp)
   do j = 1, n
      do i = 1, n
         svd_delta(i, j) = -sigma0*result%singular_left_vectors(i, sidx)*conjg(svd_vector(j))
      end do
   end do
   corrected = denominator + svd_delta
   if (sqrt(sum(abs(matmul(corrected, svd_vector))**2)) > 2.0e-11_rp) then
      error stop 'Ward fixture: exact SVD correction failed'
   end if

   magnetization = cmplx(0.0_rp, 0.0_rp, rp)
   magnetization(1) = cmplx(1.0_rp, 0.0_rp, rp)
   magnetization(3) = cmplx(1.0_rp, 0.0_rp, rp)
   call analyze_lr_ward_mode(denominator, magnetization, negative_result)
   if (.not. negative_result%succeeded) error stop 'Ward negative fixture: eigensolver failed'
   if (negative_result%magnetization_mode_fraction(negative_result%nearest_zero_index) > 0.95_rp) then
      error stop 'Ward negative fixture: distributed magnetization was classified as single-mode'
   end if
   if (negative_result%eigen_reconstruction_residual > reconstruction_tolerance .or. &
       negative_result%ward_reconstruction_residual > reconstruction_tolerance) then
      error stop 'Ward negative fixture: modal reconstruction failed'
   end if

   write (*, '(a)') 'UnitLrWardModeAnalysis: PASS (non-Hermitian biorthogonal and SVD fixtures)'
end program test_lr_ward_mode_analysis
