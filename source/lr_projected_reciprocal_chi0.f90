!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief DRESP-02 projected reciprocal bare transverse susceptibility.
!>
!> This module owns only the site-space bare response.  The Lehmann backend
!> accumulates DRESP-01 transition amplitudes directly as site outer products.
!> The GF backend is an independent real-energy Kubo bubble using the shared
!> reciprocal resolvent builder and the DRESP-01 affine site vertices.  Neither
!> backend constructs a point-space susceptibility in production.
!------------------------------------------------------------------------------
module lr_projected_reciprocal_chi0_mod

   use, intrinsic :: iso_fortran_env, only: int64
   use precision_mod, only: rp
   use response_angular_basis_mod, only: response_angular_pi
   use lr_projected_site_spin_mod, only: projected_site_spin_contract
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, &
      lmto_product_channel_plus, lmto_product_channel_minus
   use lr_pauli_transition_vertex_mod, only: pauli_endpoint_state
   use lr_ks_susceptibility_mod, only: lr_electronic_state, lr_channel_plus, lr_channel_minus, &
      lr_fermi_dirac_occupation
   use lr_gf_susceptibility_mod, only: build_weighted_resolvent
   implicit none
   private

   real(rp), parameter :: state_match_tolerance = 2.0e-11_rp
   real(rp), parameter :: endpoint_match_tolerance = 2.0e-11_rp
   real(rp), parameter :: occupation_match_tolerance = 2.0e-10_rp

   character(len=*), parameter, public :: projected_backend_lehmann = 'lehmann'
   character(len=*), parameter, public :: projected_backend_gf = 'real-axis-gf'

   !> One projected reciprocal response request.  The DRESP-01 contract and
   !> product basis are separate references so the same certified operator can
   !> be compared with both compact product-space oracles.
   type, public :: projected_chi0_request
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = lr_channel_plus
      integer :: integration_points = 2001
      real(rp) :: integration_eta = 0.0_rp
      real(rp) :: energy_margin = 1.0_rp
      type(projected_site_spin_contract), pointer :: contract => null()
      type(lmto_product_response_basis), pointer :: product_basis => null()
      type(lr_electronic_state), pointer :: electronic_state => null()
      type(lr_electronic_state), pointer :: q_endpoint_state => null()
   end type projected_chi0_request

   type, public :: projected_chi0_result
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      real(rp) :: integration_eta = 0.0_rp
      real(rp) :: energy_min = 0.0_rp
      real(rp) :: energy_max = 0.0_rp
      real(rp) :: energy_spacing = 0.0_rp
      character(len=32) :: channel = ''
      character(len=16) :: backend = ''
      character(len=8) :: selector = ''
      character(len=128) :: response_representation = ''
      character(len=512) :: provenance = ''
      integer :: nsite = 0
      integer :: product_dimension = 0
      integer :: integration_points = 0
      integer :: ntransitions_evaluated = 0
      integer :: noccupation_skips = 0
      logical :: point_response_allocated = .false.
      integer(int64) :: site_vertex_memory_bytes = 0_int64
      integer(int64) :: susceptibility_memory_bytes = 0_int64
      real(rp) :: wall_time_seconds = 0.0_rp
      real(rp) :: cpu_time_seconds = 0.0_rp
      complex(rp), allocatable :: susceptibility(:, :, :) ! (site,site,frequency)
   end type projected_chi0_result

   public :: evaluate_projected_lehmann_chi0
   public :: evaluate_projected_gf_chi0

contains

   subroutine evaluate_projected_lehmann_chi0(request, result)
      type(projected_chi0_request), intent(in) :: request
      type(projected_chi0_result), intent(out) :: result

      type(projected_site_spin_contract), pointer :: contract
      type(lmto_product_response_basis), pointer :: product
      type(lr_electronic_state), pointer :: left_state, right_state
      type(pauli_endpoint_state) :: left_band, right_band
      complex(rp), allocatable :: transition(:)
      complex(rp) :: denominator, pair_factor
      real(rp) :: weight_sum, occupation_difference
      real(rp) :: cpu_start, cpu_stop
      integer(int64) :: wall_start, wall_stop, clock_rate
      integer :: channel_kind, ik, ib, jb, ifrequency, site

      call validate_request(request, channel_kind)
      contract => request%contract
      product => request%product_basis
      left_state => request%electronic_state
      right_state => request%q_endpoint_state
      call left_state%validate('evaluate_projected_lehmann_chi0:left_state')
      call right_state%validate('evaluate_projected_lehmann_chi0:q_endpoint_state')

      call system_clock(wall_start, clock_rate)
      call cpu_time(cpu_start)
      allocate(result%susceptibility(contract%nsite, contract%nsite, size(request%frequencies)), &
         result%frequencies(size(request%frequencies)), transition(contract%nsite))
      result%susceptibility = cmplx(0.0_rp, 0.0_rp, rp)
      weight_sum = sum(left_state%k_weights)

      ! DRESP-02 production Lehmann path:
      ! eigenpairs -> certified DRESP-01 T_i -> site outer product.
      do ik = 1, left_state%nk
         do ib = 1, left_state%nbands
            call left_band%initialize(left_state%eigenvalues(ib, ik), left_state%eigenvectors(:, ib, ik))
            do jb = 1, right_state%nbands
               occupation_difference = left_state%occupations(ib, ik) - right_state%occupations(jb, ik)
               if (occupation_difference == 0.0_rp) then
                  result%noccupation_skips = result%noccupation_skips + 1
                  cycle
               end if
               call right_band%initialize(right_state%eigenvalues(jb, ik), right_state%eigenvectors(:, jb, ik))
               call contract%transition_amplitudes(product, left_band, right_band, transition)
               result%ntransitions_evaluated = result%ntransitions_evaluated + 1
               do ifrequency = 1, size(request%frequencies)
                  denominator = cmplx(request%frequencies(ifrequency) + left_band%energy - right_band%energy, &
                     request%eta, rp)
                  pair_factor = cmplx(2.0_rp*left_state%k_weights(ik)/weight_sum, 0.0_rp, rp)* &
                     occupation_difference/denominator
                  do site = 1, contract%nsite
                     result%susceptibility(:, site, ifrequency) = result%susceptibility(:, site, ifrequency) + &
                        pair_factor*transition*conjg(transition(site))
                  end do
               end do
            end do
         end do
      end do

      call system_clock(wall_stop)
      call cpu_time(cpu_stop)
      call finalize_result(result, request, contract, product, channel_kind, projected_backend_lehmann, 0.0_rp, &
         0.0_rp, 0.0_rp, wall_start, wall_stop, clock_rate, cpu_start, cpu_stop, 0_int64)
   end subroutine evaluate_projected_lehmann_chi0

   subroutine evaluate_projected_gf_chi0(request, result)
      type(projected_chi0_request), intent(in) :: request
      type(projected_chi0_result), intent(out) :: result

      type(projected_site_spin_contract), pointer :: contract
      type(lmto_product_response_basis), pointer :: product
      type(lr_electronic_state), pointer :: left_state, right_state
      complex(rp), allocatable :: vertices(:, :, :, :)
      complex(rp), allocatable :: left_gr(:, :, :), left_ga(:, :, :), left_a(:, :, :)
      complex(rp), allocatable :: right_gr(:, :, :), right_ga(:, :, :), right_a(:, :, :)
      real(rp) :: integration_eta, energy_min, energy_max, step, energy, fermi_weight
      real(rp) :: quadrature_weight, weight_sum
      real(rp) :: cpu_start, cpu_stop
      integer(int64) :: wall_start, wall_stop, clock_rate
      integer :: channel_kind, nbasis, ne, ie, ik, ifrequency

      call validate_request(request, channel_kind)
      if (request%integration_points < 3 .or. mod(request%integration_points, 2) == 0) then
         error stop 'evaluate_projected_gf_chi0: integration_points must be odd and at least three'
      end if
      if (request%energy_margin <= 0.0_rp) then
         error stop 'evaluate_projected_gf_chi0: energy_margin must be positive'
      end if
      contract => request%contract
      product => request%product_basis
      left_state => request%electronic_state
      right_state => request%q_endpoint_state
      call left_state%validate('evaluate_projected_gf_chi0:left_state')
      call right_state%validate('evaluate_projected_gf_chi0:q_endpoint_state')

      integration_eta = request%integration_eta
      if (integration_eta <= 0.0_rp) integration_eta = request%eta/40.0_rp
      if (integration_eta >= request%eta) then
         error stop 'evaluate_projected_gf_chi0: integration_eta must be smaller than response eta'
      end if

      call system_clock(wall_start, clock_rate)
      call cpu_time(cpu_start)
      call contract%site_component_vertex_tensor(product, vertices)
      nbasis = left_state%nbasis
      ne = request%integration_points
      energy_min = min(minval(left_state%eigenvalues), minval(right_state%eigenvalues)) - request%energy_margin
      energy_max = max(maxval(left_state%eigenvalues), maxval(right_state%eigenvalues)) + request%energy_margin
      if (energy_max <= energy_min) error stop 'evaluate_projected_gf_chi0: invalid integration interval'
      step = (energy_max - energy_min)/real(ne - 1, rp)
      weight_sum = sum(left_state%k_weights)

      allocate(result%susceptibility(contract%nsite, contract%nsite, size(request%frequencies)), &
         result%frequencies(size(request%frequencies)), &
         left_gr(nbasis, nbasis, 3), left_ga(nbasis, nbasis, 3), left_a(nbasis, nbasis, 3), &
         right_gr(nbasis, nbasis, 3), right_ga(nbasis, nbasis, 3), right_a(nbasis, nbasis, 3))
      result%susceptibility = cmplx(0.0_rp, 0.0_rp, rp)

      ! This is the certified LR-GF real-axis construction.  The two spectral
      ! functions and the two Kubo terms are retained verbatim; only the
      ! vertex/accumulator dimension is site x site.
      do ie = 1, ne
         energy = energy_min + real(ie - 1, rp)*step
         if (ie == 1 .or. ie == ne) then
            quadrature_weight = 1.0_rp
         else if (mod(ie, 2) == 0) then
            quadrature_weight = 4.0_rp
         else
            quadrature_weight = 2.0_rp
         end if
         quadrature_weight = quadrature_weight*step/3.0_rp
         fermi_weight = lr_fermi_dirac_occupation(energy, left_state%fermi_level, left_state%temperature)

         do ik = 1, left_state%nk
            call build_weighted_resolvent(left_state, ik, cmplx(energy, integration_eta, rp), left_gr)
            call build_weighted_resolvent(left_state, ik, cmplx(energy, -integration_eta, rp), left_ga)
            call build_weighted_resolvent(right_state, ik, cmplx(energy, integration_eta, rp), right_gr)
            call build_weighted_resolvent(right_state, ik, cmplx(energy, -integration_eta, rp), right_ga)
            left_a = cmplx(0.0_rp, 1.0_rp/(2.0_rp*response_angular_pi), rp)*(left_gr - left_ga)
            right_a = cmplx(0.0_rp, 1.0_rp/(2.0_rp*response_angular_pi), rp)*(right_gr - right_ga)
            do ifrequency = 1, size(request%frequencies)
               call build_weighted_resolvent(right_state, ik, &
                  cmplx(energy + request%frequencies(ifrequency), request%eta, rp), right_gr)
               call build_weighted_resolvent(left_state, ik, &
                  cmplx(energy - request%frequencies(ifrequency), -request%eta, rp), left_ga)
               call accumulate_projected_gf_bubble(vertices, left_a, right_a, right_gr, left_ga, &
                  fermi_weight*quadrature_weight*2.0_rp*left_state%k_weights(ik)/weight_sum, &
                  result%susceptibility(:, :, ifrequency))
            end do
         end do
      end do

      call system_clock(wall_stop)
      call cpu_time(cpu_stop)
      call finalize_result(result, request, contract, product, channel_kind, projected_backend_gf, integration_eta, &
         energy_min, energy_max, wall_start, wall_stop, clock_rate, cpu_start, cpu_stop, &
         int(size(vertices), int64)*int(storage_size(vertices)/8, int64))
      result%integration_points = ne
      result%energy_spacing = step
      deallocate(vertices, left_gr, left_ga, left_a, right_gr, right_ga, right_a)
   end subroutine evaluate_projected_gf_chi0

   subroutine accumulate_projected_gf_bubble(vertices, left_a, right_a, right_gr, left_ga, scale, response)
      complex(rp), intent(in) :: vertices(:, :, :, :), left_a(:, :, :), right_a(:, :, :), &
         right_gr(:, :, :), left_ga(:, :, :)
      real(rp), intent(in) :: scale
      complex(rp), intent(inout) :: response(:, :)

      integer :: component_i, component_j, left_power_i, right_power_i, left_power_j, right_power_j
      integer :: combined_left, combined_right, site_i, site_j
      complex(rp) :: temporary(size(left_a, 1), size(left_a, 2))
      complex(rp) :: vertex_adjoint(size(left_a, 1), size(left_a, 2))

      if (size(vertices, 4) /= size(response, 1) .or. size(response, 2) /= size(response, 1)) then
         error stop 'accumulate_projected_gf_bubble: site response shape mismatch'
      end if
      do component_i = 1, 4
         left_power_i = mod(component_i - 1, 2)
         right_power_i = (component_i - 1)/2
         if (maxval(abs(vertices(:, :, component_i, :))) == 0.0_rp) cycle
         do component_j = 1, 4
            left_power_j = mod(component_j - 1, 2)
            right_power_j = (component_j - 1)/2
            if (maxval(abs(vertices(:, :, component_j, :))) == 0.0_rp) cycle
            combined_left = left_power_i + left_power_j + 1
            combined_right = right_power_i + right_power_j + 1

            ! First Kubo term: f(E) Tr[A_L V_i G_R^R V_j^dagger].
            do site_i = 1, size(response, 1)
               if (maxval(abs(vertices(:, :, component_i, site_i))) == 0.0_rp) cycle
               temporary = matmul(left_a(:, :, combined_left), vertices(:, :, component_i, site_i))
               temporary = matmul(temporary, right_gr(:, :, combined_right))
               do site_j = 1, size(response, 2)
                  if (maxval(abs(vertices(:, :, component_j, site_j))) == 0.0_rp) cycle
                  response(site_i, site_j) = response(site_i, site_j) + scale* &
                     sum(temporary*conjg(vertices(:, :, component_j, site_j)))
               end do
            end do

            ! Second Kubo term: f(E) Tr[A_R V_j^dagger G_L^A V_i].
            do site_j = 1, size(response, 2)
               if (maxval(abs(vertices(:, :, component_j, site_j))) == 0.0_rp) cycle
               vertex_adjoint = conjg(transpose(vertices(:, :, component_j, site_j)))
               temporary = matmul(right_a(:, :, combined_right), vertex_adjoint)
               temporary = matmul(temporary, left_ga(:, :, combined_left))
               do site_i = 1, size(response, 1)
                  if (maxval(abs(vertices(:, :, component_i, site_i))) == 0.0_rp) cycle
                  response(site_i, site_j) = response(site_i, site_j) + scale* &
                     sum(temporary*transpose(vertices(:, :, component_i, site_i)))
               end do
            end do
         end do
      end do
   end subroutine accumulate_projected_gf_bubble

   subroutine validate_request(request, channel_kind)
      type(projected_chi0_request), intent(in) :: request
      integer, intent(out) :: channel_kind

      type(projected_site_spin_contract), pointer :: contract
      type(lmto_product_response_basis), pointer :: product
      type(lr_electronic_state), pointer :: left_state, right_state
      real(rp) :: expected_k(3), scale, expected_occupation
      integer :: ik, ib, expected_nbasis

      if (.not. associated(request%contract) .or. .not. associated(request%product_basis) .or. &
          .not. associated(request%electronic_state) .or. .not. associated(request%q_endpoint_state)) then
         error stop 'DRESP-02: request references are incomplete'
      end if
      if (.not. allocated(request%frequencies) .or. size(request%frequencies) < 1) then
         error stop 'DRESP-02: at least one frequency is required'
      end if
      if (request%eta <= 0.0_rp) error stop 'DRESP-02: response eta must be positive'
      select case (trim(request%channel))
      case ('plus', 'chi_plus')
         channel_kind = 1
      case ('minus', 'chi_minus')
         channel_kind = 2
      case default
         error stop 'DRESP-02: channel must be chi_plus or chi_minus'
      end select

      contract => request%contract
      product => request%product_basis
      left_state => request%electronic_state
      right_state => request%q_endpoint_state
      if (product%circular_channel /= merge(lmto_product_channel_plus, lmto_product_channel_minus, channel_kind == 1)) then
         error stop 'DRESP-02: product/channel provenance differs'
      end if
      if (.not. allocated(product%blocks) .or. product%nsite /= contract%nsite .or. &
          product%orbital_lmax /= contract%orbital_lmax .or. product%response_lmax /= contract%response_lmax .or. &
          product%npoint /= contract%npoint .or. product%product_dimension < 1) then
         error stop 'DRESP-02: product representation does not match DRESP-01 contract'
      end if
      expected_nbasis = 2*(contract%orbital_lmax + 1)**2*contract%nsite
      if (left_state%nbasis /= expected_nbasis .or. right_state%nbasis /= expected_nbasis .or. &
          left_state%nk /= right_state%nk .or. left_state%nbands /= right_state%nbands) then
         error stop 'DRESP-02: endpoint dimensions are incompatible with DRESP-01'
      end if
      call left_state%validate('DRESP-02:left_state')
      call right_state%validate('DRESP-02:q_endpoint_state')
      if (maxval(abs(left_state%k_weights - right_state%k_weights)) > state_match_tolerance) then
         error stop 'DRESP-02: endpoint k weights differ'
      end if
      scale = max(1.0_rp, abs(left_state%fermi_level), abs(right_state%fermi_level), &
         abs(left_state%temperature), abs(right_state%temperature), abs(left_state%energy_zero), &
         abs(right_state%energy_zero))
      if (abs(left_state%fermi_level - right_state%fermi_level) > state_match_tolerance*scale .or. &
          abs(left_state%temperature - right_state%temperature) > state_match_tolerance*scale .or. &
          abs(left_state%energy_zero - right_state%energy_zero) > state_match_tolerance*scale) then
         error stop 'DRESP-02: EF/temperature/energy-zero provenance differs'
      end if
      if (trim(left_state%reciprocal_mode) /= 'ham_only' .or. trim(right_state%reciprocal_mode) /= 'ham_only' .or. &
          trim(left_state%hamiltonian_order) /= 'second' .or. trim(right_state%hamiltonian_order) /= 'second' .or. &
          .not. left_state%orthogonal .or. .not. right_state%orthogonal .or. .not. left_state%collinear .or. &
          .not. right_state%collinear .or. left_state%has_soc .or. right_state%has_soc .or. &
          left_state%has_extra_operator .or. right_state%has_extra_operator) then
         error stop 'DRESP-02: state is outside the certified collinear reciprocal baseline'
      end if
      do ik = 1, left_state%nk
         expected_k = fold_fractional_kpoint(left_state%k_points(:, ik) + request%q)
         if (maxval(abs(expected_k - right_state%k_points(:, ik))) > endpoint_match_tolerance) then
            error stop 'DRESP-02: k+q endpoint is not the exact folded endpoint'
         end if
         do ib = 1, left_state%nbands
            expected_occupation = lr_fermi_dirac_occupation(left_state%eigenvalues(ib, ik), &
               left_state%fermi_level, left_state%temperature)
            if (abs(expected_occupation - left_state%occupations(ib, ik)) > occupation_match_tolerance) then
               error stop 'DRESP-02: left occupations are not the certified Fermi snapshot'
            end if
            expected_occupation = lr_fermi_dirac_occupation(right_state%eigenvalues(ib, ik), &
               right_state%fermi_level, right_state%temperature)
            if (abs(expected_occupation - right_state%occupations(ib, ik)) > occupation_match_tolerance) then
               error stop 'DRESP-02: endpoint occupations are not the certified Fermi snapshot'
            end if
         end do
      end do
   end subroutine validate_request

   subroutine finalize_result(result, request, contract, product, channel_kind, backend, integration_eta, &
                              energy_min, energy_max, wall_start, wall_stop, clock_rate, cpu_start, cpu_stop, vertex_bytes)
      type(projected_chi0_result), intent(inout) :: result
      type(projected_chi0_request), intent(in) :: request
      type(projected_site_spin_contract), intent(in) :: contract
      type(lmto_product_response_basis), intent(in) :: product
      integer, intent(in) :: channel_kind
      character(len=*), intent(in) :: backend
      real(rp), intent(in) :: integration_eta, energy_min, energy_max, cpu_start, cpu_stop
      integer(int64), intent(in) :: wall_start, wall_stop, clock_rate, vertex_bytes

      result%q = request%q
      result%frequencies = request%frequencies
      result%eta = request%eta
      result%integration_eta = integration_eta
      result%energy_min = energy_min
      result%energy_max = energy_max
      if (channel_kind == 1) then
         result%channel = lr_channel_plus
      else
         result%channel = lr_channel_minus
      end if
      result%backend = backend
      result%selector = contract%selector
      result%response_representation = 'DRESP-02 direct site x site bare chi0'
      write (result%provenance, '(a,a,a,a,a,i0,a,i0,a,es12.4,a,es12.4,a,l1)') &
         'selector=', trim(contract%selector), ' channel=', trim(result%channel), &
         ' nsite=', contract%nsite, ' product_dimension=', product%product_dimension, &
         ' response_eta=', request%eta, ' integration_eta=', integration_eta, &
         ' point_response_allocated=', result%point_response_allocated
      result%nsite = contract%nsite
      result%product_dimension = product%product_dimension
      if (result%integration_points == 0) result%integration_points = 1
      result%site_vertex_memory_bytes = vertex_bytes
      result%susceptibility_memory_bytes = int(size(result%susceptibility), int64)* &
         int(storage_size(result%susceptibility)/8, int64)
      result%wall_time_seconds = real(max(0_int64, wall_stop - wall_start), rp)/real(max(1_int64, clock_rate), rp)
      result%cpu_time_seconds = max(0.0_rp, cpu_stop - cpu_start)
   end subroutine finalize_result

   pure function fold_fractional_kpoint(k_point) result(folded)
      real(rp), intent(in) :: k_point(3)
      real(rp) :: folded(3)

      folded = k_point - floor(k_point + 0.5_rp)
   end function fold_fractional_kpoint

end module lr_projected_reciprocal_chi0_mod
