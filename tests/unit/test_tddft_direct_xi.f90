!------------------------------------------------------------------------------
! WR-01 -- generic weighted vertices and direct Xi assembly
!------------------------------------------------------------------------------
program test_tddft_direct_xi
   use precision_mod, only: rp
   use response_components_mod, only: RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel, site_projected_operator, weighted_transition_vertex, &
      weighted_transition_vectors
   use tddft_chi0_mod, only: tddft_chi0_options, tddft_chi0_result, build_chi_ks_from_eigenpairs
   use tddft_xi_mod, only: tddft_direct_xi_result, build_direct_xi_from_eigenpairs, &
      build_direct_xi_from_k_dependent_eigenpairs
   implicit none

   real(rp), parameter :: tol = 1024.0_rp*epsilon(1.0_rp)
   logical :: failed

   failed = .false.
   call test_weighted_vertex_interface()
   call test_uniform_scalar_kernel_reduction()
   call test_unequal_orbital_oracle()
   call test_batched_complex_q_and_metadata()
   call test_k_resolved_pair_operator_retention()
   if (failed) error stop 1
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_weighted_vertex_interface()
      complex(rp) :: operator(2, 2), bra(2), ket(2), bra_batch(2, 2), ket_batch(2, 2), vertices(2, 2)
      complex(rp) :: operators(2, 2, 2), expected

      operator = reshape([cmplx(0.2_rp, 0.1_rp, rp), cmplx(-0.3_rp, 0.4_rp, rp), &
         cmplx(0.5_rp, -0.2_rp, rp), cmplx(-0.1_rp, 0.7_rp, rp)], [2, 2])
      bra = [cmplx(0.7_rp, -0.1_rp, rp), cmplx(-0.2_rp, 0.5_rp, rp)]
      ket = [cmplx(-0.4_rp, 0.6_rp, rp), cmplx(0.3_rp, 0.2_rp, rp)]
      expected = explicit_vertex(operator, bra, ket)
      call check_complex('generic weighted vertex is an ordered matrix element', &
         weighted_transition_vertex(operator, bra, ket), expected)
      operators(:, :, 1) = operator; operators(:, :, 2) = transpose(operator)
      bra_batch(:, 1) = bra; ket_batch(:, 1) = ket
      bra_batch(:, 2) = ket; ket_batch(:, 2) = bra
      call weighted_transition_vectors(operators, bra_batch, ket_batch, vertices)
      call check_complex('weighted vertex batch first operator', vertices(1, 1), explicit_vertex(operators(:, :, 1), bra, ket))
      call check_complex('weighted vertex batch second ordered pair', vertices(2, 2), &
         explicit_vertex(operators(:, :, 2), ket, bra))
   end subroutine test_weighted_vertex_interface

   subroutine test_uniform_scalar_kernel_reduction()
      type(response_channel) :: left(1), right(1)
      type(tddft_chi0_options) :: options
      type(tddft_chi0_result) :: chi
      type(tddft_direct_xi_result) :: xi
      complex(rp), allocatable :: q(:, :), operators(:, :, :)
      real(rp) :: weights(1), eval(2, 1), omega(2)
      complex(rp) :: evec(2, 2, 1)
      real(rp), parameter :: scalar_kernel = -0.025_rp

      left(1) = response_channel(1, RESPONSE_PLUS); right(1) = response_channel(1, RESPONSE_MINUS)
      weights = 1.0_rp; eval(:, 1) = [-0.05_rp, 0.05_rp]; omega = [0.03_rp, 0.08_rp]
      evec = cmplx(0.0_rp, 0.0_rp, rp); evec(1, 1, 1) = 1.0_rp; evec(2, 2, 1) = 1.0_rp
      options%eta = 0.004_rp; options%fermi_level = 0.0_rp; options%electronic_temperature = 300.0_rp
      q = scalar_kernel*site_projected_operator(right(1), [1])
      allocate(operators(2, 2, 1)); operators(:, :, 1) = q
      call build_chi_ks_from_eigenpairs(weights, eval, evec, eval, evec, [1], left, right, omega, options, chi)
      call build_direct_xi_from_eigenpairs(weights, eval, evec, eval, evec, [1], left, operators, omega, options, xi)
      call check_complex_vector('direct Xi equals chi_KS K for uniform scalar operator', xi%xi(1, 1, :), &
         scalar_kernel*chi%chi(1, 1, :))
   end subroutine test_uniform_scalar_kernel_reduction

   subroutine test_unequal_orbital_oracle()
      type(response_channel) :: left(1), right(1)
      type(tddft_chi0_options) :: options
      type(tddft_chi0_result) :: chi
      type(tddft_direct_xi_result) :: xi
      real(rp), parameter :: b(2) = [0.20_rp, 0.80_rp], moment = 2.0_rp
      real(rp) :: weights(1), eval(4, 1), omega(1), scalar_kernel
      complex(rp) :: evec(4, 4, 1), operators(4, 4, 1), explicit(1)

      left(1) = response_channel(1, RESPONSE_PLUS); right(1) = response_channel(1, RESPONSE_MINUS)
      weights = 1.0_rp; eval(:, 1) = [-b(1), -b(2), b(1), b(2)]; omega = 0.0_rp
      evec = cmplx(0.0_rp, 0.0_rp, rp); evec(1, 1, 1) = 1.0_rp; evec(2, 2, 1) = 1.0_rp
      evec(3, 3, 1) = 1.0_rp; evec(4, 4, 1) = 1.0_rp
      operators = cmplx(0.0_rp, 0.0_rp, rp)
      operators(3, 1, 1) = -b(1)/moment; operators(4, 2, 1) = -b(2)/moment
      options%eta = 1.0e-14_rp; options%fermi_level = 0.0_rp; options%electronic_temperature = 0.0_rp
      options%band_first = 1; options%band_last = 4
      options%use_batched_accumulation = .false.
      call build_chi_ks_from_eigenpairs(weights, eval, evec, eval, evec, [2], left, right, omega, options, chi)
      call build_direct_xi_from_eigenpairs(weights, eval, evec, eval, evec, [2], left, operators, omega, options, xi)
      explicit = explicit_direct_xi(weights, eval, evec, eval, evec, left(1), [2], operators, omega, options)
      scalar_kernel = -sum(b)/(2.0_rp*moment**2)
      call check_complex('direct Xi matches explicit unequal-orbital oracle', xi%xi(1, 1, 1), explicit(1))
      call check_complex('weighted unequal-orbital Ward identity', xi%xi(1, 1, 1), cmplx(1.0_rp, 0.0_rp, rp))
      call check_true('direct Xi differs from old site scalar kernel', &
         abs(xi%xi(1, 1, 1) - scalar_kernel*chi%chi(1, 1, 1)) > 0.25_rp)
   end subroutine test_unequal_orbital_oracle

   subroutine test_batched_complex_q_and_metadata()
      type(response_channel) :: left(1), right(1)
      type(tddft_chi0_options) :: scalar_options, batched_options, changed_options
      type(tddft_chi0_result) :: chi
      type(tddft_direct_xi_result) :: scalar, batched, changed
      real(rp) :: weights(2), eval(2, 2), evalq(2, 2), omega(3)
      complex(rp) :: evec(2, 2, 2), evecq(2, 2, 2), operators(2, 2, 1), explicit(3), wrong(3)

      left(1) = response_channel(1, RESPONSE_PLUS); right(1) = response_channel(1, RESPONSE_MINUS)
      weights = [1.0_rp, 3.0_rp]; omega = [0.03_rp, 0.09_rp, 0.16_rp]
      call build_complex_q_fixture(eval, evalq, evec, evecq)
      operators(:, :, 1) = reshape([cmplx(0.13_rp, -0.04_rp, rp), cmplx(-0.31_rp, 0.17_rp, rp), &
         cmplx(0.22_rp, 0.35_rp, rp), cmplx(-0.07_rp, 0.11_rp, rp)], [2, 2])
      scalar_options%eta = 0.007_rp; scalar_options%fermi_level = 0.015_rp
      scalar_options%electronic_temperature = 850.0_rp; scalar_options%band_first = 1; scalar_options%band_last = 2
      scalar_options%occupation_prune_tolerance = 1.0e-14_rp; scalar_options%k_mesh_shape = [2, 1, 1]
      scalar_options%use_batched_accumulation = .false.
      batched_options = scalar_options; batched_options%use_batched_accumulation = .true.; batched_options%transition_batch_size = 3
      call build_direct_xi_from_eigenpairs(weights, eval, evec, evalq, evecq, [1], left, operators, omega, scalar_options, scalar)
      call build_direct_xi_from_eigenpairs(weights, eval, evec, evalq, evecq, [1], left, operators, omega, batched_options, batched)
      explicit = explicit_direct_xi(weights, eval, evec, evalq, evecq, left(1), [1], operators, omega, scalar_options)
      wrong = explicit_direct_xi_wrong_right_order(weights, eval, evec, evalq, evecq, left(1), [1], operators, omega, scalar_options)
      call check_complex_vector('batched direct Xi equals scalar reference for complex q fixture', batched%xi(1, 1, :), scalar%xi(1, 1, :))
      call check_complex_vector('direct Xi matches explicit complex q oracle', scalar%xi(1, 1, :), explicit)
      call check_true('wrong right-vertex conjugation/order is distinguishable', maxval(abs(scalar%xi(1, 1, :)-wrong)) > 1.0e-5_rp)
      call build_chi_ks_from_eigenpairs(weights, eval, evec, evalq, evecq, [1], left, right, omega, batched_options, chi)
      call check_real('eta metadata matches chi_KS', batched%metadata%eta, chi%metadata%eta)
      call check_real('Fermi metadata matches chi_KS', batched%metadata%fermi_level, chi%metadata%fermi_level)
      call check_real('temperature metadata matches chi_KS', batched%metadata%electronic_temperature, chi%metadata%electronic_temperature)
      call check_true('band, mesh, occupation, and batching metadata match chi_KS', &
         batched%metadata%band_first == chi%metadata%band_first .and. batched%metadata%band_last == chi%metadata%band_last .and. &
         all(batched%metadata%k_mesh_shape == chi%metadata%k_mesh_shape) .and. &
         batched%metadata%occupation_prune_tolerance == chi%metadata%occupation_prune_tolerance .and. &
         batched%metadata%batched_accumulation .and. batched%metadata%transition_batch_size == chi%metadata%transition_batch_size)
      changed_options = scalar_options; changed_options%eta = 0.019_rp
      call build_direct_xi_from_eigenpairs(weights, eval, evec, evalq, evecq, [1], left, operators, omega, changed_options, changed)
      call check_true('eta is propagated into direct Xi denominators', maxval(abs(changed%xi-scalar%xi)) > 1.0e-6_rp)
      changed_options = scalar_options; changed_options%electronic_temperature = 18000.0_rp
      call build_direct_xi_from_eigenpairs(weights, eval, evec, evalq, evecq, [1], left, operators, omega, changed_options, changed)
      call check_true('temperature occupations are propagated into direct Xi', maxval(abs(changed%xi-scalar%xi)) > 1.0e-6_rp)
      changed_options = scalar_options; changed_options%band_last = 1
      call build_direct_xi_from_eigenpairs(weights, eval, evec, evalq, evecq, [1], left, operators, omega, changed_options, changed)
      call check_true('band window is propagated into direct Xi', maxval(abs(changed%xi-scalar%xi)) > 1.0e-6_rp)
   end subroutine test_batched_complex_q_and_metadata

   subroutine test_k_resolved_pair_operator_retention()
      type(response_channel) :: left(1)
      type(tddft_chi0_options) :: options
      type(tddft_direct_xi_result) :: static_xi, k_resolved_xi, changed_xi
      real(rp) :: weights(2), eval(2, 2), evalq(2, 2), omega(2)
      complex(rp) :: evec(2, 2, 2), evecq(2, 2, 2), operators(2, 2, 1), operators_k(2, 2, 1, 2)

      left(1) = response_channel(1, RESPONSE_PLUS)
      weights = [1.0_rp, 3.0_rp]; omega = [0.03_rp, 0.09_rp]
      call build_complex_q_fixture(eval, evalq, evec, evecq)
      operators(:, :, 1) = reshape([cmplx(0.13_rp, -0.04_rp, rp), cmplx(-0.31_rp, 0.17_rp, rp), &
         cmplx(0.22_rp, 0.35_rp, rp), cmplx(-0.07_rp, 0.11_rp, rp)], [2, 2])
      operators_k(:, :, 1, 1) = operators(:, :, 1); operators_k(:, :, 1, 2) = operators(:, :, 1)
      options%eta = 0.007_rp; options%fermi_level = 0.015_rp; options%electronic_temperature = 850.0_rp
      options%band_first = 1; options%band_last = 2; options%use_batched_accumulation = .true.; options%transition_batch_size = 3
      call build_direct_xi_from_eigenpairs(weights, eval, evec, evalq, evecq, [1], left, operators, omega, options, static_xi)
      call build_direct_xi_from_k_dependent_eigenpairs(weights, eval, evec, evalq, evecq, [1], left, operators_k, omega, options, &
         k_resolved_xi)
      call check_complex_vector('k-resolved direct Xi reduces to static direct-Xi oracle when Q is k independent', &
         k_resolved_xi%xi(1, 1, :), static_xi%xi(1, 1, :))
      operators_k(:, :, 1, 2) = 1.7_rp*operators_k(:, :, 1, 2)
      call build_direct_xi_from_k_dependent_eigenpairs(weights, eval, evec, evalq, evecq, [1], left, operators_k, omega, options, &
         changed_xi)
      call check_true('k-resolved direct Xi retains distinct pair potential at each k', &
         maxval(abs(changed_xi%xi-k_resolved_xi%xi)) > 1.0e-6_rp)
      call check_true('k-resolved direct Xi records its eigenpair provenance', &
         trim(k_resolved_xi%metadata%backend) == 'direct_xi_k_resolved_eigenpairs')
   end subroutine test_k_resolved_pair_operator_retention

   subroutine build_complex_q_fixture(eval, evalq, evec, evecq)
      real(rp), intent(out) :: eval(:, :), evalq(:, :)
      complex(rp), intent(out) :: evec(:, :, :), evecq(:, :, :)
      integer :: ik
      real(rp) :: theta, phi

      do ik = 1, 2
         eval(:, ik) = [-0.12_rp + 0.01_rp*ik, 0.09_rp + 0.02_rp*ik]
         evalq(:, ik) = [-0.10_rp - 0.015_rp*ik, 0.11_rp + 0.005_rp*ik]
         theta = 0.15_rp*ik; phi = 0.23_rp*ik
         call unitary_spinor_pair(theta, phi, evec(:, :, ik))
         call unitary_spinor_pair(theta + 0.11_rp, phi - 0.19_rp, evecq(:, :, ik))
      end do
   end subroutine build_complex_q_fixture

   subroutine unitary_spinor_pair(theta, phi, vectors)
      real(rp), intent(in) :: theta, phi
      complex(rp), intent(out) :: vectors(2, 2)
      vectors(1, 1) = cmplx(cos(theta), 0.0_rp, rp)
      vectors(2, 1) = sin(theta)*exp(cmplx(0.0_rp, phi, rp))
      vectors(1, 2) = -sin(theta)*exp(cmplx(0.0_rp, -phi, rp))
      vectors(2, 2) = cmplx(cos(theta), 0.0_rp, rp)
   end subroutine unitary_spinor_pair

   function explicit_direct_xi(weights, eval, evec, evalq, evecq, left, counts, operators, omega, options) result(xi)
      real(rp), intent(in) :: weights(:), eval(:, :), evalq(:, :), omega(:)
      complex(rp), intent(in) :: evec(:, :, :), evecq(:, :, :), operators(:, :, :)
      type(response_channel), intent(in) :: left
      integer, intent(in) :: counts(:)
      type(tddft_chi0_options), intent(in) :: options
      complex(rp) :: xi(size(omega)), v, w
      integer :: ik, n, m, iw

      xi = cmplx(0.0_rp, 0.0_rp, rp)
      do ik = 1, size(weights); do n = options%band_first, options%band_last; do m = options%band_first, options%band_last
         v = explicit_vertex(site_projected_operator(left, counts), evec(:, n, ik), evecq(:, m, ik))
         w = explicit_vertex(operators(:, :, 1), evecq(:, m, ik), evec(:, n, ik))
         do iw = 1, size(omega)
            xi(iw) = xi(iw) + weights(ik)/sum(weights)*(fermi(eval(n, ik), options)-fermi(evalq(m, ik), options))*v*w/ &
               cmplx(omega(iw)+eval(n, ik)-evalq(m, ik), options%eta, rp)
         end do
      end do; end do; end do
   end function explicit_direct_xi

   function explicit_direct_xi_wrong_right_order(weights, eval, evec, evalq, evecq, left, counts, operators, omega, options) result(xi)
      real(rp), intent(in) :: weights(:), eval(:, :), evalq(:, :), omega(:)
      complex(rp), intent(in) :: evec(:, :, :), evecq(:, :, :), operators(:, :, :)
      type(response_channel), intent(in) :: left
      integer, intent(in) :: counts(:)
      type(tddft_chi0_options), intent(in) :: options
      complex(rp) :: xi(size(omega)), v, w
      integer :: ik, n, m, iw

      xi = cmplx(0.0_rp, 0.0_rp, rp)
      do ik = 1, size(weights); do n = options%band_first, options%band_last; do m = options%band_first, options%band_last
         v = explicit_vertex(site_projected_operator(left, counts), evec(:, n, ik), evecq(:, m, ik))
         w = conjg(explicit_vertex(operators(:, :, 1), evec(:, n, ik), evecq(:, m, ik)))
         do iw = 1, size(omega)
            xi(iw) = xi(iw) + weights(ik)/sum(weights)*(fermi(eval(n, ik), options)-fermi(evalq(m, ik), options))*v*w/ &
               cmplx(omega(iw)+eval(n, ik)-evalq(m, ik), options%eta, rp)
         end do
      end do; end do; end do
   end function explicit_direct_xi_wrong_right_order

   pure function explicit_vertex(operator, bra, ket) result(value)
      complex(rp), intent(in) :: operator(:, :), bra(:), ket(:)
      complex(rp) :: value
      value = sum(conjg(bra)*matmul(operator, ket))
   end function explicit_vertex

   pure real(rp) function fermi(energy, options) result(value)
      real(rp), intent(in) :: energy
      type(tddft_chi0_options), intent(in) :: options
      real(rp) :: x, kt
      kt = max(options%electronic_temperature*6.3336814e-6_rp, 1.0e-10_rp)
      x = (energy-options%fermi_level)/kt
      if (x >= 50.0_rp) then; value = 0.0_rp
      else if (x <= -50.0_rp) then; value = 1.0_rp
      else; value = 1.0_rp/(exp(x)+1.0_rp)
      end if
   end function fermi

   subroutine check_complex(label, actual, expected)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual, expected
      if (abs(actual-expected) > tol*max(1.0_rp, abs(expected))) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, abs(actual-expected); failed = .true.
      end if
   end subroutine check_complex

   subroutine check_complex_vector(label, actual, expected)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual(:), expected(:)
      call check_complex(label, cmplx(maxval(abs(actual-expected)), 0.0_rp, rp), cmplx(0.0_rp, 0.0_rp, rp))
   end subroutine check_complex_vector

   subroutine check_real(label, actual, expected)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: actual, expected
      if (abs(actual-expected) > tol*max(1.0_rp, abs(expected))) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, abs(actual-expected); failed = .true.
      end if
   end subroutine check_real

   subroutine check_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      if (.not. condition) then; write (*, '(a,1x,a)') 'FAIL', label; failed = .true.
      end if
   end subroutine check_true

end program test_tddft_direct_xi
