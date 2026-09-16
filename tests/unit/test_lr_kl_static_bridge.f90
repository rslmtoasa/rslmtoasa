!------------------------------------------------------------------------------
! DRESP-03 finite-Hamiltonian K/L bridge fixture.
!
! The expected contraction below is intentionally written as an independent
! band-pair loop.  It does not call the bridge evaluator or reuse its helper.
!------------------------------------------------------------------------------
program test_lr_kl_static_bridge
   use, intrinsic :: ieee_arithmetic
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use lr_ks_susceptibility_mod, only: lr_electronic_state, lr_fermi_dirac_occupation
   use lr_kl_static_bridge_mod, only: kl_static_bridge_request, kl_static_bridge_result, &
      kl_scalar_vertex_diagnostic, evaluate_kl_static_bridge, build_ham_only_exchange_vertex, &
      assess_scalar_exchange_vertex, contract_scalar_site_chi0
   implicit none

   call basis_init(2)
   call run_fixture()
   write (*, '(a)') 'UnitLrKlStaticBridge: PASS (operator vertices, scalar audit, independent oracle, q/site order, eta ladder)'

contains

   subroutine run_fixture()
      integer, parameter :: norb = 2, nsite = 2, nbasis = 8, nbands = 8, nk = 1
      real(rp), parameter :: temperature = 0.02_rp
      real(rp), parameter :: q(3) = [0.13_rp, -0.07_rp, 0.0_rp]
      complex(rp) :: hamiltonian(nbasis, nbasis), up_blocks(norb, norb, nsite), down_blocks(norb, norb, nsite)
      complex(rp) :: operator_up(norb, norb, nsite), operator_down(norb, norb, nsite)
      complex(rp), allocatable, target :: vertices(:, :, :), operator_vertices(:, :, :)
      real(rp) :: values(nbands, nk), endpoint_values(nbands, nk), points(3, nk), endpoint_points(3, nk), weights(nk)
      real(rp) :: occupations(nbands, nk)
      complex(rp) :: vectors(nbasis, nbands, nk), endpoint_vectors(nbasis, nbands, nk)
      logical :: selected_l(0:2)
      real(rp) :: scalar(2), residual, offdiag
      type(lr_electronic_state), target :: state, endpoint
      type(kl_static_bridge_request) :: request
      type(kl_static_bridge_result) :: result, negative_result, coarse, eta04, eta02, fine, eta005, zero_eta
      type(kl_scalar_vertex_diagnostic) :: scalar_diagnostic, operator_diagnostic
      complex(rp) :: expected(nsite, nsite), expected_negative(nsite, nsite)
      complex(rp) :: site_chi(nsite, nsite), scalar_exchange(nsite, nsite)

      selected_l = .true.
      up_blocks = cmplx(0.0_rp, 0.0_rp, rp)
      down_blocks = cmplx(0.0_rp, 0.0_rp, rp)
      up_blocks(:, :, 1) = reshape([cmplx(-0.55_rp, 0.0_rp, rp), cmplx(0.0_rp, 0.0_rp, rp), &
         cmplx(0.0_rp, 0.0_rp, rp), cmplx(-0.25_rp, 0.0_rp, rp)], [norb, norb])
      down_blocks(:, :, 1) = reshape([cmplx(-0.15_rp, 0.0_rp, rp), cmplx(0.0_rp, 0.0_rp, rp), &
         cmplx(0.0_rp, 0.0_rp, rp), cmplx(0.15_rp, 0.0_rp, rp)], [norb, norb])
      up_blocks(:, :, 2) = reshape([cmplx(-0.48_rp, 0.0_rp, rp), cmplx(0.0_rp, 0.0_rp, rp), &
         cmplx(0.0_rp, 0.0_rp, rp), cmplx(-0.18_rp, 0.0_rp, rp)], [norb, norb])
      down_blocks(:, :, 2) = reshape([cmplx(0.22_rp, 0.0_rp, rp), cmplx(0.0_rp, 0.0_rp, rp), &
         cmplx(0.0_rp, 0.0_rp, rp), cmplx(0.52_rp, 0.0_rp, rp)], [norb, norb])
      scalar = [up_blocks(1, 1, 1) - down_blocks(1, 1, 1), up_blocks(1, 1, 2) - down_blocks(1, 1, 2)]

      hamiltonian = cmplx(0.0_rp, 0.0_rp, rp)
      call put_local_spin_blocks(hamiltonian, up_blocks, down_blocks)
      ! Real intersite hopping keeps the q/-q phase/site-order fixture free of
      ! an additional complex time-reversal-breaking ingredient.
      call add_hopping(hamiltonian, 1, 2, 1, 0.043_rp)
      call add_hopping(hamiltonian, 1, 2, 2, -0.027_rp)
      call diagonalize(hamiltonian, values(:, 1), vectors(:, :, 1))
      endpoint_values = values
      endpoint_vectors(:, :, 1) = vectors(:, :, 1)
      call apply_site_phase(endpoint_vectors(:, :, 1), q)
      points = 0.0_rp
      endpoint_points(:, 1) = q
      weights = 1.0_rp
      call build_occupations(values, occupations, temperature)
      call state%initialize(values, vectors, points, weights, occupations, 0.0_rp, temperature)
      call endpoint%initialize(endpoint_values, endpoint_vectors, endpoint_points, weights, occupations, 0.0_rp, temperature)

      call build_ham_only_exchange_vertex(up_blocks, down_blocks, selected_l, vertices)
      call assess_scalar_exchange_vertex(vertices, selected_l, scalar, scalar_diagnostic, 1.0e-13_rp)
      if (.not. scalar_diagnostic%exact .or. scalar_diagnostic%max_absolute_residual > 1.0e-13_rp) then
         error stop 'UnitLrKlStaticBridge: exact scalar fixture was rejected'
      end if

      request%q = q
      request%eta = 0.025_rp
      request%vertex_source = 'finite fixture H_up-H_down'
      request%electronic_state => state
      request%q_endpoint_state => endpoint
      request%exchange_vertex => vertices
      call evaluate_kl_static_bridge(request, result)
      call independent_operator_oracle(state, endpoint, vertices, request%eta, expected)
      residual = maxval(abs(result%exchange - expected))
      if (residual > 2.0e-13_rp) error stop 'UnitLrKlStaticBridge: independent operator oracle mismatch'
      if (.not. all(ieee_is_finite(real(result%exchange, rp))) .or. &
          .not. all(ieee_is_finite(aimag(result%exchange)))) error stop 'UnitLrKlStaticBridge: non-finite exchange result'
      offdiag = max(abs(result%exchange(1, 2)), abs(result%exchange(2, 1)))
      if (offdiag <= 1.0e-12_rp) error stop 'UnitLrKlStaticBridge: site ordering fixture did not produce off-diagonal exchange'

      ! The scalar relation is checked against an independently accumulated
      ! unit spin-flip susceptibility.  It is only invoked for the fixture
      ! proven scalar by the audit above.
      call independent_unit_spin_chi(state, endpoint, request%eta, site_chi)
      call contract_scalar_site_chi0(site_chi, scalar, scalar_exchange)
      residual = maxval(abs(result%exchange - scalar_exchange))
      if (residual > 2.0e-13_rp) error stop 'UnitLrKlStaticBridge: scalar chi contraction prefactor mismatch'

      call build_negative_endpoint(vectors(:, :, 1), endpoint_values(:, 1), endpoint_vectors(:, :, 1), endpoint_points(:, 1), -q)
      call endpoint%restore_to_default()
      call endpoint%initialize(endpoint_values, endpoint_vectors, reshape(-q, [3, 1]), weights, occupations, 0.0_rp, temperature)
      request%q = -q
      call evaluate_kl_static_bridge(request, negative_result)
      call independent_operator_oracle(state, endpoint, vertices, request%eta, expected_negative)
      if (maxval(abs(negative_result%exchange - expected_negative)) > 2.0e-13_rp) then
         error stop 'UnitLrKlStaticBridge: negative-q independent oracle mismatch'
      end if
      if (maxval(abs(negative_result%exchange)) <= 1.0e-12_rp .or. &
          maxval(abs(result%exchange - negative_result%exchange)) <= 1.0e-12_rp) then
         error stop 'UnitLrKlStaticBridge: q phase/site ordering fixture was ineffective'
      end if

      ! A deliberately orbital-dependent splitting is not silently compressed
      ! to a scalar.  The operator route remains the certified object.
      operator_up = up_blocks
      operator_down = down_blocks
      operator_down(2, 2, 1) = operator_up(2, 2, 1) + 0.25_rp
      operator_down(1, 2, 2) = 0.018_rp
      call build_ham_only_exchange_vertex(operator_up, operator_down, selected_l, operator_vertices)
      call assess_scalar_exchange_vertex(operator_vertices, selected_l, scalar, operator_diagnostic, 1.0e-13_rp)
      if (operator_diagnostic%exact .or. operator_diagnostic%max_absolute_residual <= 1.0e-4_rp) then
         error stop 'UnitLrKlStaticBridge: non-scalar operator was incorrectly accepted as scalar'
      end if
      request%exchange_vertex => operator_vertices
      request%vertex_source = 'finite fixture orbital-dependent operator'
      call evaluate_kl_static_bridge(request, coarse)
      call independent_operator_oracle(state, endpoint, operator_vertices, request%eta, expected)
      if (maxval(abs(coarse%exchange - expected)) > 2.0e-13_rp) then
         error stop 'UnitLrKlStaticBridge: operator-valued fixture mismatch'
      end if

      ! Static limit control: the same transition contraction is evaluated down
      ! to eta=0; no fitted rescaling or Goldstone correction is used.
      endpoint_vectors(:, :, 1) = vectors(:, :, 1)
      call apply_site_phase(endpoint_vectors(:, :, 1), q)
      call endpoint%restore_to_default()
      call endpoint%initialize(endpoint_values, endpoint_vectors, reshape(q, [3, 1]), weights, occupations, 0.0_rp, temperature)
      request%q = q
      request%exchange_vertex => vertices
      request%eta = 0.08_rp
      call evaluate_kl_static_bridge(request, coarse)
      request%eta = 0.04_rp
      call evaluate_kl_static_bridge(request, eta04)
      request%eta = 0.02_rp
      call evaluate_kl_static_bridge(request, eta02)
      request%eta = 0.01_rp
      call evaluate_kl_static_bridge(request, fine)
      request%eta = 0.005_rp
      call evaluate_kl_static_bridge(request, eta005)
      request%eta = 0.0_rp
      call evaluate_kl_static_bridge(request, zero_eta)
      if (.not. all(ieee_is_finite(real(zero_eta%exchange, rp))) .or. &
          .not. all(ieee_is_finite(aimag(zero_eta%exchange)))) error stop &
         'UnitLrKlStaticBridge: static eta limit is non-finite'
      write (*, '(a,es12.4)') '  finite fixture scalar residual = ', scalar_diagnostic%max_absolute_residual
      write (*, '(a,es12.4)') '  finite fixture operator residual = ', operator_diagnostic%max_absolute_residual
      write (*, '(a,es12.4)') '  finite fixture oracle residual = ', residual
      write (*, '(a,es12.4)') '  finite fixture |J(.04)-J(.08)| = ', maxval(abs(eta04%exchange - coarse%exchange))
      write (*, '(a,es12.4)') '  finite fixture |J(.02)-J(.04)| = ', maxval(abs(eta02%exchange - eta04%exchange))
      write (*, '(a,es12.4)') '  finite fixture |J(.01)-J(.02)| = ', maxval(abs(fine%exchange - eta02%exchange))
      write (*, '(a,es12.4)') '  finite fixture |J(.005)-J(.01)| = ', maxval(abs(eta005%exchange - fine%exchange))
      write (*, '(a,es12.4)') '  finite fixture |J(0)-J(.005)| = ', maxval(abs(zero_eta%exchange - eta005%exchange))
      deallocate(vertices, operator_vertices)
   end subroutine run_fixture

   subroutine put_local_spin_blocks(hamiltonian, up, down)
      complex(rp), intent(inout) :: hamiltonian(:, :)
      complex(rp), intent(in) :: up(:, :, :), down(:, :, :)
      integer :: site, iorb, jorb, norb, offset

      norb = size(up, 1)
      do site = 1, size(up, 3)
         offset = (site - 1)*2*norb
         do iorb = 1, norb
            do jorb = 1, norb
               hamiltonian(offset + iorb, offset + jorb) = up(iorb, jorb, site)
               hamiltonian(offset + norb + iorb, offset + norb + jorb) = down(iorb, jorb, site)
            end do
         end do
      end do
   end subroutine put_local_spin_blocks

   subroutine add_hopping(hamiltonian, site_i, site_j, orbital, value)
      complex(rp), intent(inout) :: hamiltonian(:, :)
      integer, intent(in) :: site_i, site_j, orbital
      real(rp), intent(in) :: value
      integer :: left, right, spin, norb

      norb = 2
      do spin = 0, 1
         left = (site_i - 1)*2*norb + spin*norb + orbital
         right = (site_j - 1)*2*norb + spin*norb + orbital
         hamiltonian(left, right) = cmplx(value, 0.0_rp, rp)
         hamiltonian(right, left) = cmplx(value, 0.0_rp, rp)
      end do
   end subroutine add_hopping

   subroutine diagonalize(hamiltonian, values, vectors)
      complex(rp), intent(in) :: hamiltonian(:, :)
      real(rp), intent(out) :: values(:)
      complex(rp), intent(out) :: vectors(:, :)
      complex(rp), allocatable :: work(:), matrix(:, :)
      real(rp), allocatable :: rwork(:)
      complex(rp) :: query(1)
      integer :: n, lwork, info
      external :: zheev

      n = size(values)
      allocate(matrix(n, n), rwork(max(1, 3*n - 2)))
      matrix = hamiltonian
      call zheev('V', 'U', n, matrix, n, values, query, -1, rwork, info)
      lwork = max(1, nint(real(query(1), rp)))
      allocate(work(lwork))
      call zheev('V', 'U', n, matrix, n, values, work, lwork, rwork, info)
      if (info /= 0) error stop 'UnitLrKlStaticBridge: fixture diagonalization failed'
      vectors = matrix
      deallocate(matrix, work, rwork)
   end subroutine diagonalize

   subroutine apply_site_phase(vectors, q)
      complex(rp), intent(inout) :: vectors(:, :)
      real(rp), intent(in) :: q(3)
      integer :: site, first, last
      complex(rp) :: phase
      real(rp) :: angle

      do site = 1, 2
         angle = 2.0_rp*acos(-1.0_rp)*q(1)*real(site - 1, rp)
         phase = cmplx(cos(angle), sin(angle), rp)
         first = (site - 1)*4 + 1
         last = site*4
         vectors(first:last, :) = phase*vectors(first:last, :)
      end do
   end subroutine apply_site_phase

   subroutine build_negative_endpoint(source_vectors, values, endpoint_vectors, endpoint_point, q)
      complex(rp), intent(in) :: source_vectors(:, :)
      real(rp), intent(in) :: values(:)
      complex(rp), intent(out) :: endpoint_vectors(:, :)
      real(rp), intent(out) :: endpoint_point(:)
      real(rp), intent(in) :: q(3)

      endpoint_vectors = source_vectors
      call apply_site_phase(endpoint_vectors, q)
      endpoint_point = -q
      if (size(values) < 1) error stop 'UnitLrKlStaticBridge: empty negative endpoint'
   end subroutine build_negative_endpoint

   subroutine build_occupations(values, occupations, temperature)
      real(rp), intent(in) :: values(:, :), temperature
      real(rp), intent(out) :: occupations(:, :)
      integer :: ib, ik

      do ik = 1, size(values, 2)
         do ib = 1, size(values, 1)
            occupations(ib, ik) = lr_fermi_dirac_occupation(values(ib, ik), 0.0_rp, temperature)
         end do
      end do
   end subroutine build_occupations

   subroutine independent_operator_oracle(state, endpoint, vertices, eta, expected)
      type(lr_electronic_state), intent(in) :: state, endpoint
      complex(rp), intent(in) :: vertices(:, :, :)
      real(rp), intent(in) :: eta
      complex(rp), intent(out) :: expected(:, :)
      complex(rp) :: transition(size(vertices, 3)), action(size(vertices, 1)), denominator
      real(rp) :: df, weight
      integer :: ib, jb, site, other

      expected = cmplx(0.0_rp, 0.0_rp, rp)
      weight = state%k_weights(1)/sum(state%k_weights)
      do ib = 1, state%nbands
         do jb = 1, endpoint%nbands
            df = state%occupations(ib, 1) - endpoint%occupations(jb, 1)
            if (abs(df) <= tiny(1.0_rp)) cycle
            denominator = cmplx(state%eigenvalues(ib, 1) - endpoint%eigenvalues(jb, 1), eta, rp)
            do site = 1, size(vertices, 3)
               action = matmul(vertices(:, :, site), endpoint%eigenvectors(:, jb, 1))
               transition(site) = dot_product(state%eigenvectors(:, ib, 1), action)
            end do
            do site = 1, size(vertices, 3)
               do other = 1, size(vertices, 3)
                  expected(site, other) = expected(site, other) - 0.25_rp*weight*df/denominator* &
                     transition(site)*conjg(transition(other))
               end do
            end do
         end do
      end do
   end subroutine independent_operator_oracle

   subroutine independent_unit_spin_chi(state, endpoint, eta, chi)
      type(lr_electronic_state), intent(in) :: state, endpoint
      real(rp), intent(in) :: eta
      complex(rp), intent(out) :: chi(:, :)
      complex(rp) :: transition(size(chi, 1)), denominator
      real(rp) :: df, weight
      integer :: ib, jb, site, other, basis, norb, offset

      chi = cmplx(0.0_rp, 0.0_rp, rp)
      norb = size(state%eigenvectors, 1)/(2*size(chi, 1))
      weight = state%k_weights(1)/sum(state%k_weights)
      do ib = 1, state%nbands
         do jb = 1, endpoint%nbands
            df = state%occupations(ib, 1) - endpoint%occupations(jb, 1)
            if (abs(df) <= tiny(1.0_rp)) cycle
            denominator = cmplx(state%eigenvalues(ib, 1) - endpoint%eigenvalues(jb, 1), eta, rp)
            do site = 1, size(chi, 1)
               offset = (site - 1)*2*norb
               transition(site) = cmplx(0.0_rp, 0.0_rp, rp)
               do basis = 1, norb
                  transition(site) = transition(site) + conjg(state%eigenvectors(offset + basis, ib, 1))* &
                     endpoint%eigenvectors(offset + norb + basis, jb, 1)
               end do
            end do
            do site = 1, size(chi, 1)
               do other = 1, size(chi, 1)
                  chi(site, other) = chi(site, other) + 2.0_rp*weight*df/denominator* &
                     transition(site)*conjg(transition(other))
               end do
            end do
         end do
      end do
   end subroutine independent_unit_spin_chi

end program test_lr_kl_static_bridge
