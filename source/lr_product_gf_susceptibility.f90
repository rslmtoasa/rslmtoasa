!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Compact weighted-orthonormal LMTO-product reciprocal-GF response.
!>
!> This is a representation-level companion to LR-GF-02.  It retains the
!> real-axis Green-function/Kubo construction and writes directly into the
!> compact product coordinates.  The legacy point-grid GF evaluator remains
!> independent and is intentionally not routed through this module.
!------------------------------------------------------------------------------
module lr_product_gf_susceptibility_mod

   use, intrinsic :: iso_fortran_env, only: int64
   use precision_mod, only: rp
   use response_angular_basis_mod, only: response_angular_pi
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, &
      lmto_product_channel_plus, lmto_product_channel_minus
   use lr_ks_susceptibility_mod, only: lr_electronic_state, lr_channel_plus, lr_channel_minus, &
      lr_fermi_dirac_occupation
   use lr_gf_susceptibility_mod, only: build_weighted_resolvent
   implicit none
   private

   real(rp), parameter :: state_match_tolerance = 2.0e-11_rp
   real(rp), parameter :: endpoint_match_tolerance = 2.0e-11_rp
   real(rp), parameter :: occupation_match_tolerance = 2.0e-10_rp

   character(len=*), parameter :: product_representation = &
      'weighted-orthonormal LMTO product representation'

   type, public :: lr_product_gf_susceptibility_request
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = lr_channel_plus
      integer :: integration_points = 2001
      real(rp) :: integration_eta = 0.0_rp
      real(rp) :: energy_margin = 1.0_rp
      type(lmto_product_response_basis), pointer :: product_basis => null()
      type(lr_electronic_state), pointer :: electronic_state => null()
      type(lr_electronic_state), pointer :: q_endpoint_state => null()
   end type lr_product_gf_susceptibility_request

   type, public :: lr_product_gf_susceptibility_result
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = ''
      integer :: product_dimension = 0
      integer :: integration_points = 0
      real(rp) :: actual_integration_eta = 0.0_rp
      real(rp) :: energy_min = 0.0_rp
      real(rp) :: energy_max = 0.0_rp
      character(len=128) :: response_representation = ''
      character(len=512) :: response_space_metadata = ''
      logical :: point_response_allocated = .false.
      integer(int64) :: component_vertex_memory_bytes = 0_int64
      integer(int64) :: gf_matrix_memory_bytes = 0_int64
      integer(int64) :: susceptibility_memory_bytes = 0_int64
      integer :: nbasis = 0
      integer :: nk = 0
      integer :: nfrequency = 0
      real(rp) :: wall_time_seconds = 0.0_rp
      real(rp) :: cpu_time_seconds = 0.0_rp
      real(rp) :: time_per_energy_kpoint = 0.0_rp
      complex(rp), allocatable :: susceptibility(:, :, :) ! (product,product,frequency)
   end type lr_product_gf_susceptibility_result

   public :: evaluate_lr_product_gf_susceptibility

contains

   subroutine evaluate_lr_product_gf_susceptibility(request, result)
      type(lr_product_gf_susceptibility_request), intent(in) :: request
      type(lr_product_gf_susceptibility_result), intent(out) :: result

      type(lmto_product_response_basis), pointer :: product_basis
      type(lr_electronic_state), pointer :: left_state, right_state
      complex(rp), allocatable :: vertices(:, :, :, :)
      complex(rp), allocatable :: left_gr(:, :, :), left_ga(:, :, :), left_a(:, :, :)
      complex(rp), allocatable :: right_gr(:, :, :), right_ga(:, :, :), right_a(:, :, :)
      real(rp) :: integration_eta, energy_min, energy_max, step, energy, fermi_weight
      real(rp) :: quadrature_weight, weight_sum, cpu_start, cpu_stop
      integer(int64) :: wall_start, wall_stop, clock_rate
      integer :: channel_kind, nbasis, nfrequency, ne, ie, ik, ifrequency

      if (.not. associated(request%product_basis) .or. .not. associated(request%electronic_state) .or. &
          .not. associated(request%q_endpoint_state)) then
         error stop 'evaluate_lr_product_gf_susceptibility: request references are incomplete'
      end if
      product_basis => request%product_basis
      left_state => request%electronic_state
      right_state => request%q_endpoint_state

      call validate_product_gf_inputs(product_basis, left_state, right_state, request, channel_kind)
      call left_state%validate('evaluate_lr_product_gf_susceptibility:left_state')
      call right_state%validate('evaluate_lr_product_gf_susceptibility:q_endpoint_state')

      integration_eta = request%integration_eta
      if (integration_eta <= 0.0_rp) integration_eta = request%eta/40.0_rp
      if (integration_eta >= request%eta) then
         error stop 'evaluate_lr_product_gf_susceptibility: integration_eta must be smaller than response eta'
      end if

      call system_clock(wall_start, clock_rate)
      call cpu_time(cpu_start)
      call product_basis%component_vertex_tensor(vertices)
      nbasis = left_state%nbasis
      nfrequency = size(request%frequencies)
      ne = request%integration_points
      energy_min = min(minval(left_state%eigenvalues), minval(right_state%eigenvalues)) - request%energy_margin
      energy_max = max(maxval(left_state%eigenvalues), maxval(right_state%eigenvalues)) + request%energy_margin
      if (energy_max <= energy_min) error stop 'evaluate_lr_product_gf_susceptibility: invalid integration interval'
      step = (energy_max - energy_min)/real(ne - 1, rp)
      weight_sum = sum(left_state%k_weights)

      allocate(result%susceptibility(product_basis%product_dimension, product_basis%product_dimension, nfrequency), &
         result%frequencies(nfrequency), &
         left_gr(nbasis, nbasis, 3), left_ga(nbasis, nbasis, 3), left_a(nbasis, nbasis, 3), &
         right_gr(nbasis, nbasis, 3), right_ga(nbasis, nbasis, 3), right_a(nbasis, nbasis, 3))
      result%susceptibility = cmplx(0.0_rp, 0.0_rp, rp)

      ! This is the unchanged LR-GF-02 real-energy Simpson integration.  The
      ! only representation change is the compact vertex and accumulator.
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

            do ifrequency = 1, nfrequency
               call build_weighted_resolvent(right_state, ik, &
                  cmplx(energy + request%frequencies(ifrequency), request%eta, rp), right_gr)
               call build_weighted_resolvent(left_state, ik, &
                  cmplx(energy - request%frequencies(ifrequency), -request%eta, rp), left_ga)
               call accumulate_product_gf_bubble(vertices, left_a, right_a, right_gr, left_ga, &
                  fermi_weight*quadrature_weight*2.0_rp*left_state%k_weights(ik)/weight_sum, &
                  result%susceptibility(:, :, ifrequency))
            end do
         end do
      end do

      call system_clock(wall_stop)
      call cpu_time(cpu_stop)
      result%q = request%q
      result%frequencies = request%frequencies
      result%eta = request%eta
      result%product_dimension = product_basis%product_dimension
      result%integration_points = ne
      result%actual_integration_eta = integration_eta
      result%energy_min = energy_min
      result%energy_max = energy_max
      result%response_representation = product_representation
      result%point_response_allocated = .false.
      result%nbasis = nbasis
      result%nk = left_state%nk
      result%nfrequency = nfrequency
      result%component_vertex_memory_bytes = int(size(vertices), int64)*int(storage_size(vertices)/8, int64)
      result%gf_matrix_memory_bytes = int(6*size(left_gr), int64)*int(storage_size(left_gr)/8, int64)
      result%susceptibility_memory_bytes = int(size(result%susceptibility), int64)* &
         int(storage_size(result%susceptibility)/8, int64)
      result%wall_time_seconds = real(max(0_int64, wall_stop - wall_start), rp)/real(max(1_int64, clock_rate), rp)
      result%cpu_time_seconds = max(0.0_rp, cpu_stop - cpu_start)
      result%time_per_energy_kpoint = result%wall_time_seconds/real(max(1, ne*left_state%nk), rp)
      if (channel_kind == 1) then
         result%channel = lr_channel_plus
      else
         result%channel = lr_channel_minus
      end if
      write (result%response_space_metadata, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,l1)') &
         'product_dimension=', product_basis%product_dimension, ' component_vertex_bytes=', &
         result%component_vertex_memory_bytes, ' gf_matrix_bytes=', result%gf_matrix_memory_bytes, &
         ' susceptibility_bytes=', result%susceptibility_memory_bytes, ' energy_points=', ne, &
         ' k_points=', left_state%nk, ' frequencies=', nfrequency, ' point_response_allocated=', &
         result%point_response_allocated
   end subroutine evaluate_lr_product_gf_susceptibility

   subroutine accumulate_product_gf_bubble(vertices, left_a, right_a, right_gr, left_ga, scale, susceptibility)
      complex(rp), intent(in) :: vertices(:, :, :, :), left_a(:, :, :), right_a(:, :, :), &
         right_gr(:, :, :), left_ga(:, :, :)
      real(rp), intent(in) :: scale
      complex(rp), intent(inout) :: susceptibility(:, :)

      integer :: component_i, component_j, left_power_i, right_power_i, left_power_j, right_power_j
      integer :: combined_left, combined_right, i, j, product_dimension
      complex(rp) :: temporary(size(left_a, 1), size(left_a, 2))
      complex(rp) :: vertex_adjoint(size(left_a, 1), size(left_a, 2))

      product_dimension = size(susceptibility, 1)
      if (size(susceptibility, 2) /= product_dimension .or. size(vertices, 4) /= product_dimension) then
         error stop 'accumulate_product_gf_bubble: product matrix shape mismatch'
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

            ! First Kubo term: f(E) Tr[A_L(E) V_I G_R^R(E+omega) V_J^dagger].
            do i = 1, product_dimension
               if (maxval(abs(vertices(:, :, component_i, i))) == 0.0_rp) cycle
               temporary = matmul(left_a(:, :, combined_left), vertices(:, :, component_i, i))
               temporary = matmul(temporary, right_gr(:, :, combined_right))
               do j = 1, product_dimension
                  if (maxval(abs(vertices(:, :, component_j, j))) == 0.0_rp) cycle
                  susceptibility(i, j) = susceptibility(i, j) + &
                     scale*sum(temporary*conjg(vertices(:, :, component_j, j)))
               end do
            end do

            ! Second Kubo term: f(E) Tr[A_R(E) V_J^dagger G_L^A(E-omega) V_I].
            do j = 1, product_dimension
               if (maxval(abs(vertices(:, :, component_j, j))) == 0.0_rp) cycle
               vertex_adjoint = conjg(transpose(vertices(:, :, component_j, j)))
               temporary = matmul(right_a(:, :, combined_right), vertex_adjoint)
               temporary = matmul(temporary, left_ga(:, :, combined_left))
               do i = 1, product_dimension
                  if (maxval(abs(vertices(:, :, component_i, i))) == 0.0_rp) cycle
                  susceptibility(i, j) = susceptibility(i, j) + &
                     scale*sum(temporary*transpose(vertices(:, :, component_i, i)))
               end do
            end do
         end do
      end do
   end subroutine accumulate_product_gf_bubble

   subroutine validate_product_gf_inputs(product_basis, left_state, right_state, request, channel_kind)
      type(lmto_product_response_basis), intent(in) :: product_basis
      type(lr_electronic_state), intent(in) :: left_state, right_state
      type(lr_product_gf_susceptibility_request), intent(in) :: request
      integer, intent(out) :: channel_kind

      real(rp) :: expected_k(3), expected_occupation, scale
      integer :: ik, ib, expected_nbasis

      if (.not. allocated(request%frequencies) .or. size(request%frequencies) < 1) then
         error stop 'evaluate_lr_product_gf_susceptibility: at least one frequency is required'
      end if
      if (request%eta <= 0.0_rp) error stop 'evaluate_lr_product_gf_susceptibility: retarded eta must be positive'
      if (request%integration_points < 3 .or. mod(request%integration_points, 2) == 0) then
         error stop 'evaluate_lr_product_gf_susceptibility: integration_points must be odd and at least three'
      end if
      if (request%energy_margin <= 0.0_rp) error stop 'evaluate_lr_product_gf_susceptibility: energy_margin must be positive'
      select case (trim(request%channel))
      case ('plus', 'chi_plus')
         channel_kind = 1
      case ('minus', 'chi_minus')
         channel_kind = 2
      case default
         error stop 'evaluate_lr_product_gf_susceptibility: channel must be chi_plus or chi_minus'
      end select
      if (product_basis%circular_channel /= merge(lmto_product_channel_plus, lmto_product_channel_minus, channel_kind == 1)) then
         error stop 'evaluate_lr_product_gf_susceptibility: product/channel provenance differs'
      end if
      if (.not. allocated(product_basis%blocks) .or. product_basis%nsite < 1 .or. &
          product_basis%product_dimension < 1 .or. product_basis%response_lmax < 0) then
         error stop 'evaluate_lr_product_gf_susceptibility: product representation is not initialized'
      end if
      expected_nbasis = 2*(product_basis%orbital_lmax + 1)**2*product_basis%nsite
      if (left_state%nbasis /= expected_nbasis .or. right_state%nbasis /= expected_nbasis) then
         error stop 'evaluate_lr_product_gf_susceptibility: eigensystem basis is incompatible with product representation'
      end if
      if (left_state%nk /= right_state%nk .or. left_state%nbands /= right_state%nbands .or. &
          left_state%nbasis /= right_state%nbasis) then
         error stop 'evaluate_lr_product_gf_susceptibility: endpoint dimensions differ'
      end if
      if (maxval(abs(left_state%k_weights - right_state%k_weights)) > state_match_tolerance) then
         error stop 'evaluate_lr_product_gf_susceptibility: k weights differ between endpoints'
      end if
      scale = max(1.0_rp, abs(left_state%fermi_level), abs(right_state%fermi_level), &
         abs(left_state%temperature), abs(right_state%temperature), abs(left_state%energy_zero), &
         abs(right_state%energy_zero))
      if (abs(left_state%fermi_level - right_state%fermi_level) > state_match_tolerance*scale .or. &
          abs(left_state%temperature - right_state%temperature) > state_match_tolerance*scale .or. &
          abs(left_state%energy_zero - right_state%energy_zero) > state_match_tolerance*scale) then
         error stop 'evaluate_lr_product_gf_susceptibility: EF/temperature/energy-zero provenance differs'
      end if
      if (trim(left_state%reciprocal_mode) /= 'ham_only' .or. trim(right_state%reciprocal_mode) /= 'ham_only' .or. &
          trim(left_state%hamiltonian_order) /= 'second' .or. trim(right_state%hamiltonian_order) /= 'second' .or. &
          .not. left_state%orthogonal .or. .not. right_state%orthogonal .or. .not. left_state%collinear .or. &
          .not. right_state%collinear .or. left_state%has_soc .or. right_state%has_soc .or. &
          left_state%has_extra_operator .or. right_state%has_extra_operator) then
         error stop 'evaluate_lr_product_gf_susceptibility: representation is outside the LR-GF-02 baseline'
      end if
      do ik = 1, left_state%nk
         expected_k = fold_fractional_kpoint(left_state%k_points(:, ik) + request%q)
         if (maxval(abs(expected_k - right_state%k_points(:, ik))) > endpoint_match_tolerance) then
            error stop 'evaluate_lr_product_gf_susceptibility: k+q endpoint is not the exact folded endpoint'
         end if
         do ib = 1, left_state%nbands
            expected_occupation = lr_fermi_dirac_occupation(left_state%eigenvalues(ib, ik), &
               left_state%fermi_level, left_state%temperature)
            if (abs(expected_occupation - left_state%occupations(ib, ik)) > occupation_match_tolerance) then
               error stop 'evaluate_lr_product_gf_susceptibility: occupations are not the certified Fermi snapshot'
            end if
            expected_occupation = lr_fermi_dirac_occupation(right_state%eigenvalues(ib, ik), &
               right_state%fermi_level, right_state%temperature)
            if (abs(expected_occupation - right_state%occupations(ib, ik)) > occupation_match_tolerance) then
               error stop 'evaluate_lr_product_gf_susceptibility: endpoint occupations are not the certified Fermi snapshot'
            end if
         end do
      end do
   end subroutine validate_product_gf_inputs

   pure function fold_fractional_kpoint(k_point) result(folded)
      real(rp), intent(in) :: k_point(3)
      real(rp) :: folded(3)

      folded = k_point - floor(k_point + 0.5_rp)
   end function fold_fractional_kpoint

end module lr_product_gf_susceptibility_mod
