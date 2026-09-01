!------------------------------------------------------------------------------
! TDDFT-03 -- independent bare transverse chi_KS oracle
!------------------------------------------------------------------------------
program test_tddft_chi_ks
   use precision_mod, only: rp
   use math_mod, only: pi
   use response_components_mod, only: RESPONSE_CHARGE, RESPONSE_MZ, RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, build_chi_ks_from_eigenpairs, &
      build_static_chi_ks_from_eigenpairs, build_static_chi_ks_from_eigenpairs_at_q, &
      tddft_static_divided_difference, write_chi_ks_text
   implicit none

   real(rp), parameter :: machine_tol = 256.0_rp*epsilon(1.0_rp)
   logical :: failed

   failed = .false.
   call test_independent_pair_sum()
   call test_batched_accumulator_equivalence()
   call test_positive_frequency_sign_and_products()
   call test_convergence_controls()
   call test_static_divided_difference_and_eta_independence()
   call test_static_nonzero_q_and_provenance()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   ! A one-orbital, spin-split tight-binding fixture.  The oracle below uses
   ! its own explicit spin-selection factors and Fermi routine; it deliberately
   ! does not call response vertices or the production occupation routine.
   subroutine test_independent_pair_sum()
      integer, parameter :: nk = 4
      integer :: iq, itemp, ieta
      integer :: q_shifts(3) = [0, 1, 2]
      real(rp) :: temperatures(2) = [250.0_rp, 900.0_rp]
      real(rp) :: etas(2) = [0.0015_rp, 0.008_rp]
      real(rp) :: omega(4) = [0.025_rp, 0.075_rp, 0.110_rp, 0.170_rp]
      real(rp) :: weights(nk), eval(2, nk), evalq(2, nk)
      complex(rp) :: evec(2, 2, nk), evecq(2, 2, nk), oracle(size(omega))
      type(tddft_chi0_options) :: options
      type(tddft_chi0_result) :: result
      type(response_channel) :: left(1), right(1)

      left(1) = response_channel(1, RESPONSE_PLUS)
      right(1) = response_channel(1, RESPONSE_MINUS)
      weights = [1.0_rp, 2.0_rp, 1.0_rp, 3.0_rp]
      do iq = 1, size(q_shifts)
         call build_spin_split_fixture(q_shifts(iq), 0.030_rp, eval, evalq, evec, evecq)
         do itemp = 1, size(temperatures)
            do ieta = 1, size(etas)
               options%eta = etas(ieta)
               options%fermi_level = 0.0_rp
               options%electronic_temperature = temperatures(itemp)
               options%k_mesh_shape = [nk, 1, 1]
               call build_chi_ks_from_eigenpairs(weights, eval, evec, evalq, evecq, [1], left, right, omega, &
                  options, result)
               call brute_force_pair_sum(weights, eval, evalq, omega, options%fermi_level, &
                  options%electronic_temperature, options%eta, oracle)
               call check_complex_vector('independent one-orbital pair sum', result%chi(1, 1, :), oracle)
               call check_real('metadata eta', result%metadata%eta, options%eta)
               call check_real('metadata Fermi level', result%metadata%fermi_level, options%fermi_level)
            end do
         end do
      end do
   end subroutine test_independent_pair_sum

   ! TDDFT-11: the BLAS transition tiles must remain numerically equivalent to
   ! the scalar, fixed-order reference path for a genuinely complex,
   ! multi-channel response.  This exercises non-Hermitian left/right channel
   ! pairs, multiple k points, and batch boundaries that do not divide n^2.
   subroutine test_batched_accumulator_equivalence()
      integer, parameter :: nk = 3, nbands = 4
      real(rp) :: weights(nk), eval(nbands, nk), evalq(nbands, nk), omega(5)
      complex(rp) :: evec(4, nbands, nk), evecq(4, nbands, nk)
      type(tddft_chi0_options) :: reference_options, batched_options
      type(tddft_chi0_result) :: reference, batched
      type(response_channel) :: left(2), right(2)

      weights = [1.0_rp, 3.0_rp, 2.0_rp]
      omega = [0.013_rp, 0.041_rp, 0.079_rp, 0.121_rp, 0.169_rp]
      call build_complex_four_spinor_fixture(eval, evalq, evec, evecq)
      left(1) = response_channel(1, RESPONSE_PLUS)
      left(2) = response_channel(1, RESPONSE_MZ)
      right(1) = response_channel(1, RESPONSE_MINUS)
      right(2) = response_channel(1, RESPONSE_CHARGE)
      reference_options%eta = 0.004_rp
      reference_options%fermi_level = 0.01_rp
      reference_options%electronic_temperature = 700.0_rp
      reference_options%use_batched_accumulation = .false.
      batched_options = reference_options
      batched_options%use_batched_accumulation = .true.
      batched_options%transition_batch_size = 3
      call build_chi_ks_from_eigenpairs(weights, eval, evec, evalq, evecq, [2], left, right, omega, &
         reference_options, reference)
      call build_chi_ks_from_eigenpairs(weights, eval, evec, evalq, evecq, [2], left, right, omega, &
         batched_options, batched)
      call check_complex_vector('scalar/Batched GEMM chi_KS equivalence', reshape(batched%chi, [size(batched%chi)]), &
         reshape(reference%chi, [size(reference%chi)]))
      call check_true('batched metadata is explicit', batched%metadata%batched_accumulation .and. &
         batched%metadata%transition_batch_size == 3)
   end subroutine test_batched_accumulator_equivalence

   subroutine test_positive_frequency_sign_and_products()
      integer, parameter :: nk = 4
      real(rp) :: omega(1) = [0.100_rp]
      real(rp) :: weights(nk), eval(2, nk), evalq(2, nk)
      complex(rp) :: evec(2, 2, nk), evecq(2, 2, nk)
      type(tddft_chi0_options) :: options
      type(tddft_chi0_result) :: result
      type(response_channel) :: left(1), right(1)

      left(1) = response_channel(1, RESPONSE_PLUS)
      right(1) = response_channel(1, RESPONSE_MINUS)
      weights = 1.0_rp
      call build_spin_split_fixture(0, 0.0_rp, eval, evalq, evec, evecq)
      options%eta = 0.002_rp
      options%fermi_level = 0.0_rp
      options%electronic_temperature = 100.0_rp
      options%k_mesh_shape = [nk, 1, 1]
      call build_chi_ks_from_eigenpairs(weights, eval, evec, evalq, evecq, [1], left, right, omega, options, result)

      ! For an occupied up state and empty down state, chi has negative Im at
      ! its positive excitation energy.  The documented Stoner quantity is
      ! therefore -Im(chi)/pi and must be positive.
      call check_true('positive-frequency Im chi is negative', result%im_chi(1, 1, 1) < 0.0_rp)
      call check_true('positive-frequency Stoner spectrum is positive', result%trace_spectrum(1) > 0.0_rp)
      call check_real('site diagonal spectrum', result%site_diagonal_spectrum(1, 1), result%trace_spectrum(1))
      call check_real('KS/Stoner map', result%stoner_spectral_map(1, 1), result%trace_spectrum(1))
      call check_real('real chi product', result%re_chi(1, 1, 1), real(result%chi(1, 1, 1), rp))
      call check_real('imag chi product', result%im_chi(1, 1, 1), aimag(result%chi(1, 1, 1)))
      write (*, '(a,1x,es16.8,1x,es24.16,1x,es24.16,1x,es24.16)') 'BASELINE dynamic omega/ReChi/ImChi/Stoner', &
         omega(1), result%re_chi(1, 1, 1), result%im_chi(1, 1, 1), result%trace_spectrum(1)
   end subroutine test_positive_frequency_sign_and_products

   subroutine test_convergence_controls()
      real(rp) :: omega(2) = [0.075_rp, 0.105_rp]
      real(rp) :: weights4(4), eval4(2, 4), evalq4(2, 4)
      real(rp) :: weights8(8), eval8(2, 8), evalq8(2, 8)
      complex(rp) :: evec4(2, 2, 4), evecq4(2, 2, 4), evec8(2, 2, 8), evecq8(2, 2, 8)
      type(tddft_chi0_options) :: options, changed
      type(tddft_chi0_result) :: reference, comparison, fine_mesh
      type(response_channel) :: left(1), right(1)

      left(1) = response_channel(1, RESPONSE_PLUS)
      right(1) = response_channel(1, RESPONSE_MINUS)
      weights4 = 1.0_rp
      weights8 = 1.0_rp
      ! A flat spin-split band has an exactly mesh-independent BZ average,
      ! making this a deterministic k-mesh convergence oracle.
      call build_spin_split_fixture(0, 0.0_rp, eval4, evalq4, evec4, evecq4)
      call build_spin_split_fixture(0, 0.0_rp, eval8, evalq8, evec8, evecq8)
      options%eta = 0.003_rp
      options%fermi_level = 0.0_rp
      options%electronic_temperature = 300.0_rp
      options%k_mesh_shape = [4, 1, 1]
      call build_chi_ks_from_eigenpairs(weights4, eval4, evec4, evalq4, evecq4, [1], left, right, omega, options, reference)
      options%k_mesh_shape = [8, 1, 1]
      call build_chi_ks_from_eigenpairs(weights8, eval8, evec8, evalq8, evecq8, [1], left, right, omega, options, fine_mesh)
      call check_complex_vector('flat-band k mesh convergence', fine_mesh%chi(1, 1, :), reference%chi(1, 1, :))
      call check_true('k mesh metadata retained', all(fine_mesh%metadata%k_mesh_shape == [8, 1, 1]))

      ! Explicit full-band selection is the exact reference path.  The knobs
      ! below are additionally checked to ensure eta and electronic smearing
      ! are not silently ignored by future optimized paths.
      changed = options
      changed%k_mesh_shape = [4, 1, 1]
      changed%band_first = 1
      changed%band_last = 2
      call build_chi_ks_from_eigenpairs(weights4, eval4, evec4, evalq4, evecq4, [1], left, right, omega, changed, comparison)
      call check_complex_vector('explicit all-band window', comparison%chi(1, 1, :), reference%chi(1, 1, :))
      changed%eta = 0.012_rp
      call build_chi_ks_from_eigenpairs(weights4, eval4, evec4, evalq4, evecq4, [1], left, right, omega, changed, comparison)
      call check_true('eta convergence control changes response', abs(comparison%chi(1, 1, 1) - reference%chi(1, 1, 1)) > 0.0_rp)
      changed%eta = options%eta
      changed%electronic_temperature = 18000.0_rp
      call build_chi_ks_from_eigenpairs(weights4, eval4, evec4, evalq4, evecq4, [1], left, right, omega, changed, comparison)
      call check_true('smearing convergence control changes response', &
         abs(comparison%chi(1, 1, 1) - reference%chi(1, 1, 1)) > 0.0_rp)
   end subroutine test_convergence_controls

   subroutine test_static_divided_difference_and_eta_independence()
      real(rp) :: finite_difference, derivative, near_value, zero_temperature_value
      real(rp) :: weights(2), eval(2, 2), evalq(2, 2)
      complex(rp) :: evec(2, 2, 2), evecq(2, 2, 2)
      type(response_channel) :: left(1), right(1)
      type(tddft_chi0_options) :: low_eta, high_eta
      type(tddft_chi0_result) :: static_low_eta, static_high_eta

      ! (f_n-f_m)/(e_n-e_m) has the sign f’(e)<0.  The exact and nearby
      ! degeneracy paths must agree with the analytic Fermi derivative.
      derivative = -0.25_rp/(300.0_rp*6.3336814e-6_rp)
      finite_difference = tddft_static_divided_difference(0.0_rp, 0.0_rp, 0.0_rp, 300.0_rp)
      near_value = tddft_static_divided_difference(1.0e-13_rp, -1.0e-13_rp, 0.0_rp, 300.0_rp)
      call check_real('static exact degeneracy uses negative Fermi derivative', finite_difference, derivative)
      call check_real('static near degeneracy is stable', near_value, derivative)
      zero_temperature_value = tddft_static_divided_difference(-0.2_rp, -0.1_rp, 0.0_rp, 0.0_rp)
      call check_real('zero-temperature occupied degeneracy has zero derivative', zero_temperature_value, 0.0_rp)

      left(1) = response_channel(1, RESPONSE_PLUS); right(1) = response_channel(1, RESPONSE_MINUS)
      weights = 1.0_rp
      call build_spin_split_fixture(0, 0.0_rp, eval, evalq, evec, evecq)
      low_eta%eta = 1.0e-5_rp; low_eta%fermi_level = 0.0_rp; low_eta%electronic_temperature = 300.0_rp
      high_eta = low_eta; high_eta%eta = 0.08_rp
      call build_static_chi_ks_from_eigenpairs(weights, eval, evec, [1], left, right, low_eta, static_low_eta)
      call build_static_chi_ks_from_eigenpairs(weights, eval, evec, [1], left, right, high_eta, static_high_eta)
      call check_complex_vector('static chi is independent of dynamic eta', static_low_eta%chi(1, 1, :), &
         static_high_eta%chi(1, 1, :))
      call check_real('static chi has no imaginary broadening', aimag(static_low_eta%chi(1, 1, 1)), 0.0_rp)
      call check_real('static metadata reports eta zero', static_low_eta%metadata%eta, 0.0_rp)
   end subroutine test_static_divided_difference_and_eta_independence

   subroutine test_static_nonzero_q_and_provenance()
      integer, parameter :: nk = 4
      real(rp) :: weights(nk), eval(2, nk), evalq(2, nk), expected
      complex(rp) :: evec(2, 2, nk), evecq(2, 2, nk)
      real(rp) :: omega_static(1) = [0.0_rp]
      type(tddft_chi0_options) :: options
      type(tddft_chi0_result) :: result
      type(response_channel) :: left(1), right(1)
      integer :: unit, ios
      character(len=256) :: line
      logical :: found_q, found_endpoint, found_eta_role, found_grid, found_landau

      left(1) = response_channel(1, RESPONSE_PLUS)
      right(1) = response_channel(1, RESPONSE_MINUS)
      weights = [1.0_rp, 2.0_rp, 1.0_rp, 3.0_rp]
      call build_spin_split_fixture(1, 0.030_rp, eval, evalq, evec, evecq)
      options%eta = 0.0_rp
      options%fermi_level = 0.0_rp
      options%electronic_temperature = 700.0_rp
      options%k_mesh_shape = [nk, 1, 1]
      options%q_direct = [0.0_rp, 0.0_rp, 0.25_rp]
      call build_static_chi_ks_from_eigenpairs_at_q(weights, eval, evec, evalq, evecq, [1], left, right, options, result)
      call brute_force_static_pair_sum(weights, eval, evalq, options%fermi_level, options%electronic_temperature, expected)
      call check_complex_vector('finite-q static endpoint pair sum', result%chi(1, 1, :), [cmplx(expected, 0.0_rp, rp)])
      call check_true('finite-q static limit is explicit', result%metadata%static_limit .and. &
         .not. result%metadata%eta_is_numerical)
      call check_true('finite-q static endpoint provenance is explicit', index(result%metadata%endpoint_provenance, 'separate') > 0)
      call check_real('finite-q static eta is unused', result%metadata%eta, 0.0_rp)
      call check_complex_vector('finite-q static q provenance', cmplx(result%metadata%q_direct, 0.0_rp, rp), &
         cmplx(options%q_direct, 0.0_rp, rp))

      call write_chi_ks_text('unit_tddft_chi0_provenance.dat', omega_static, result)
      found_q = .false.; found_endpoint = .false.; found_eta_role = .false.; found_grid = .false.; found_landau = .false.
      open(newunit=unit, file='unit_tddft_chi0_provenance.dat', status='old', action='read', iostat=ios)
      if (ios == 0) then
         do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            found_q = found_q .or. index(line, '# q_direct =') == 1
            found_endpoint = found_endpoint .or. index(line, '# endpoint_provenance =') == 1
            found_eta_role = found_eta_role .or. index(line, '# eta_role =') == 1
            found_grid = found_grid .or. index(line, '# omega_grid_min_max_points =') == 1
            found_landau = found_landau .or. index(line, 'stoner_landau') > 0
         end do
         close(unit, status='delete')
      end if
      call check_true('chi0 text writes q provenance', found_q)
      call check_true('chi0 text writes endpoint provenance', found_endpoint)
      call check_true('chi0 text separates eta role', found_eta_role)
      call check_true('chi0 text writes omega-grid provenance', found_grid)
      call check_true('chi0 text writes optional Stoner Landau map', found_landau)
      write (*, '(a,1x,es24.16,1x,es24.16)') 'BASELINE static q=0.25/ReChi/Stoner', &
         real(result%chi(1, 1, 1), rp), result%trace_spectrum(1)
   end subroutine test_static_nonzero_q_and_provenance

   subroutine build_spin_split_fixture(q_shift, hopping, eval, evalq, evec, evecq)
      integer, intent(in) :: q_shift
      real(rp), intent(in) :: hopping
      real(rp), intent(out) :: eval(:, :), evalq(:, :)
      complex(rp), intent(out) :: evec(:, :, :), evecq(:, :, :)
      integer :: nk, ik, ikq
      real(rp) :: k, kq

      nk = size(eval, 2)
      do ik = 1, nk
         ikq = 1 + modulo(ik - 1 + q_shift, nk)
         k = real(ik - 1, rp)/real(nk, rp)
         kq = real(ikq - 1, rp)/real(nk, rp)
         eval(1, ik) = -0.050_rp - 2.0_rp*hopping*cos(2.0_rp*pi*k)
         eval(2, ik) =  0.050_rp - 2.0_rp*hopping*cos(2.0_rp*pi*k)
         evalq(1, ik) = -0.050_rp - 2.0_rp*hopping*cos(2.0_rp*pi*kq)
         evalq(2, ik) =  0.050_rp - 2.0_rp*hopping*cos(2.0_rp*pi*kq)
         evec(:, :, ik) = cmplx(0.0_rp, 0.0_rp, rp)
         evecq(:, :, ik) = cmplx(0.0_rp, 0.0_rp, rp)
         evec(1, 1, ik) = cmplx(1.0_rp, 0.0_rp, rp)
         evec(2, 2, ik) = cmplx(1.0_rp, 0.0_rp, rp)
         evecq(1, 1, ik) = cmplx(1.0_rp, 0.0_rp, rp)
         evecq(2, 2, ik) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
   end subroutine build_spin_split_fixture

   subroutine build_complex_four_spinor_fixture(eval, evalq, evec, evecq)
      real(rp), intent(out) :: eval(:, :), evalq(:, :)
      complex(rp), intent(out) :: evec(:, :, :), evecq(:, :, :)
      integer :: ik, ib, ic
      real(rp) :: norm

      do ik = 1, size(eval, 2)
         do ib = 1, size(eval, 1)
            eval(ib, ik) = -0.12_rp + 0.075_rp*real(ib-1, rp) + 0.006_rp*real(ik-1, rp)
            evalq(ib, ik) = -0.10_rp + 0.079_rp*real(ib-1, rp) - 0.004_rp*real(ik-1, rp)
            do ic = 1, size(evec, 1)
               evec(ic, ib, ik) = cmplx(real(2*ic+ib+ik, rp), real(ic-2*ib+ik, rp), rp)
               evecq(ic, ib, ik) = cmplx(real(ic+2*ib-ik, rp), real(2*ic-ib+ik, rp), rp)
            end do
            norm = sqrt(sum(abs(evec(:, ib, ik))**2)); evec(:, ib, ik) = evec(:, ib, ik)/norm
            norm = sqrt(sum(abs(evecq(:, ib, ik))**2)); evecq(:, ib, ik) = evecq(:, ib, ik)/norm
         end do
      end do
   end subroutine build_complex_four_spinor_fixture

   subroutine brute_force_pair_sum(weights, eval, evalq, omega, fermi_level, temperature, eta, chi)
      real(rp), intent(in) :: weights(:), eval(:, :), evalq(:, :), omega(:), fermi_level, temperature, eta
      complex(rp), intent(out) :: chi(:)
      integer :: ik, n, m, iw
      real(rp) :: fn, fm, spin_factor

      chi = cmplx(0.0_rp, 0.0_rp, rp)
      do iw = 1, size(omega)
         do ik = 1, size(weights)
            do n = 1, size(eval, 1)
               do m = 1, size(evalq, 1)
                  spin_factor = 0.0_rp
                  ! <up|m+|down><down|m-|up> = 2*2.  This is intentionally
                  ! written directly, independent of response_vertices_mod.
                  if (n == 1 .and. m == 2) spin_factor = 4.0_rp
                  fn = independent_fermi(eval(n, ik), fermi_level, temperature)
                  fm = independent_fermi(evalq(m, ik), fermi_level, temperature)
                  chi(iw) = chi(iw) + weights(ik)/sum(weights)*(fn - fm)*spin_factor/ &
                     cmplx(omega(iw) + eval(n, ik) - evalq(m, ik), eta, rp)
               end do
            end do
         end do
      end do
   end subroutine brute_force_pair_sum

   subroutine brute_force_static_pair_sum(weights, eval, evalq, fermi_level, temperature, chi)
      real(rp), intent(in) :: weights(:), eval(:, :), evalq(:, :), fermi_level, temperature
      real(rp), intent(out) :: chi
      integer :: ik
      real(rp) :: delta, fn, fm, factor

      chi = 0.0_rp
      do ik = 1, size(weights)
         delta = eval(1, ik) - evalq(2, ik)
         fn = independent_fermi(eval(1, ik), fermi_level, temperature)
         fm = independent_fermi(evalq(2, ik), fermi_level, temperature)
         if (abs(delta) > 1.0e-12_rp) then
            factor = (fn - fm)/delta
         else
            factor = -0.25_rp/(temperature*6.3336814e-6_rp)
         end if
         chi = chi + weights(ik)/sum(weights)*4.0_rp*factor
      end do
   end subroutine brute_force_static_pair_sum

   pure real(rp) function independent_fermi(e, ef, temperature) result(f)
      real(rp), intent(in) :: e, ef, temperature
      real(rp) :: x, kt

      kt = max(temperature*6.3336814e-6_rp, 1.0e-10_rp)
      x = (e - ef)/kt
      if (x >= 50.0_rp) then
         f = 0.0_rp
      else if (x <= -50.0_rp) then
         f = 1.0_rp
      else
         f = 1.0_rp/(exp(x) + 1.0_rp)
      end if
   end function independent_fermi

   subroutine check_complex_vector(label, actual, expected)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual(:), expected(:)
      real(rp) :: error, scale

      error = maxval(abs(actual - expected))
      scale = max(1.0_rp, maxval(abs(expected)))
      if (error > machine_tol*scale) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, error
         failed = .true.
      end if
   end subroutine check_complex_vector

   subroutine check_real(label, actual, expected)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: actual, expected
      if (abs(actual - expected) > machine_tol*max(1.0_rp, abs(expected))) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, abs(actual - expected)
         failed = .true.
      end if
   end subroutine check_real

   subroutine check_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine check_true

end program test_tddft_chi_ks
