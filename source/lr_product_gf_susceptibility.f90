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
   character(len=*), parameter, public :: lr_product_gf_contraction_scalar = 'scalar'
   character(len=*), parameter, public :: lr_product_gf_contraction_optimized = 'optimized'
   character(len=*), parameter, public :: lr_product_gf_contraction_factorized = 'factorized'

   type, public :: lr_product_gf_susceptibility_request
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = lr_channel_plus
      integer :: integration_points = 2001
      real(rp) :: integration_eta = 0.0_rp
      real(rp) :: energy_margin = 1.0_rp
      character(len=16) :: contraction_backend = lr_product_gf_contraction_factorized
      type(lmto_product_response_basis), pointer :: product_basis => null()
      type(lr_electronic_state), pointer :: electronic_state => null()
      type(lr_electronic_state), pointer :: q_endpoint_state => null()
      ! Optional DRESP-01 orbital selector; absent means the complete product
      ! response used by the pre-DRESP oracle.
      logical, allocatable :: selected_l(:)
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
      real(rp) :: energy_spacing = 0.0_rp
      real(rp) :: spacing_over_integration_eta = 0.0_rp
      character(len=128) :: response_representation = ''
      character(len=512) :: response_space_metadata = ''
      character(len=16) :: contraction_backend = ''
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
      integer(int64) :: band_transition_memory_bytes = 0_int64
      integer(int64) :: band_kernel_memory_bytes = 0_int64
      complex(rp), allocatable :: susceptibility(:, :, :) ! (product,product,frequency)
   end type lr_product_gf_susceptibility_result

   type :: lr_product_gf_contraction_workspace
      complex(rp), allocatable :: vertices_flat(:, :, :)
      complex(rp), allocatable :: vertices_transpose(:, :, :)
      complex(rp), allocatable :: transformed(:, :)
      complex(rp), allocatable :: contribution(:, :)
      complex(rp), allocatable :: temporary(:, :)
      logical, allocatable :: active_component(:)
   end type lr_product_gf_contraction_workspace

   public :: evaluate_lr_product_gf_susceptibility, build_lr_product_gf_transition_amplitudes

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
      integer :: channel_kind, contraction_kind, nbasis, nfrequency, ne, ie, ik, ifrequency
      logical :: use_optimized_contraction
      type(lr_product_gf_contraction_workspace) :: contraction_workspace

      if (.not. associated(request%product_basis) .or. .not. associated(request%electronic_state) .or. &
          .not. associated(request%q_endpoint_state)) then
         error stop 'evaluate_lr_product_gf_susceptibility: request references are incomplete'
      end if
      product_basis => request%product_basis
      left_state => request%electronic_state
      right_state => request%q_endpoint_state

      call validate_product_gf_inputs(product_basis, left_state, right_state, request, channel_kind, contraction_kind)
      use_optimized_contraction = contraction_kind == 2
      call left_state%validate('evaluate_lr_product_gf_susceptibility:left_state')
      call right_state%validate('evaluate_lr_product_gf_susceptibility:q_endpoint_state')

      integration_eta = request%integration_eta
      if (integration_eta <= 0.0_rp) integration_eta = request%eta/40.0_rp
      if (integration_eta >= request%eta) then
         error stop 'evaluate_lr_product_gf_susceptibility: integration_eta must be smaller than response eta'
      end if

      if (contraction_kind == 1) then
         call evaluate_lr_product_gf_factorized(request, result, product_basis, left_state, right_state, integration_eta)
         return
      end if
      call system_clock(wall_start, clock_rate)
      call cpu_time(cpu_start)
      if (allocated(request%selected_l)) then
         call product_basis%component_vertex_tensor(vertices, request%selected_l)
      else
         call product_basis%component_vertex_tensor(vertices)
      end if
      if (use_optimized_contraction) call initialize_product_gf_contraction_workspace(contraction_workspace, vertices)
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
               if (use_optimized_contraction) then
                  call accumulate_product_gf_bubble_optimized(vertices, left_a, right_a, right_gr, left_ga, &
                     fermi_weight*quadrature_weight*2.0_rp*left_state%k_weights(ik)/weight_sum, &
                     result%susceptibility(:, :, ifrequency), contraction_workspace)
               else
                  call accumulate_product_gf_bubble_scalar(vertices, left_a, right_a, right_gr, left_ga, &
                     fermi_weight*quadrature_weight*2.0_rp*left_state%k_weights(ik)/weight_sum, &
                     result%susceptibility(:, :, ifrequency))
               end if
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
      result%energy_spacing = step
      result%spacing_over_integration_eta = step/integration_eta
      result%response_representation = product_representation
      if (use_optimized_contraction) then
         result%contraction_backend = lr_product_gf_contraction_optimized
      else
         result%contraction_backend = lr_product_gf_contraction_scalar
      end if
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
      write (result%response_space_metadata, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,l1,a,es12.4,a,es12.4,a,a,a,a)') &
         'product_dimension=', product_basis%product_dimension, ' component_vertex_bytes=', &
         result%component_vertex_memory_bytes, ' gf_matrix_bytes=', result%gf_matrix_memory_bytes, &
         ' susceptibility_bytes=', result%susceptibility_memory_bytes, ' energy_points=', ne, &
         ' k_points=', left_state%nk, ' frequencies=', nfrequency, ' point_response_allocated=', &
         result%point_response_allocated, ' energy_spacing=', result%energy_spacing, &
         ' h_over_integration_eta=', result%spacing_over_integration_eta, ' contraction_backend=', &
         trim(result%contraction_backend), ' accepted_state=', 'shared_request_snapshots'
   end subroutine evaluate_lr_product_gf_susceptibility

   !> Evaluate the same real-axis Kubo bubble after commuting the finite
   !> band sums through the energy integral.  The component vertices are
   !> transformed to the two endpoint eigenbases once per k point.  The
   !> energy loop then integrates only the scalar band-pair kernel, and the
   !> complete compact response is formed once per k point and frequency.
   subroutine evaluate_lr_product_gf_factorized(request, result, product_basis, left_state, right_state, integration_eta)
      type(lr_product_gf_susceptibility_request), intent(in) :: request
      type(lr_product_gf_susceptibility_result), intent(out) :: result
      type(lmto_product_response_basis), intent(in) :: product_basis
      type(lr_electronic_state), intent(in) :: left_state, right_state
      real(rp), intent(in) :: integration_eta

      complex(rp), allocatable :: vertices(:, :, :, :), transitions(:, :, :), transition_flat(:, :)
      complex(rp), allocatable :: kernel(:, :), weighted_transitions(:, :), contribution(:, :)
      complex(rp), allocatable :: left_spectral(:), right_spectral(:), right_retarded(:), left_advanced(:)
      complex(rp) :: spectral_prefactor
      real(rp) :: energy_min, energy_max, step, energy, fermi_weight, quadrature_weight, weight_sum
      real(rp) :: cpu_start, cpu_stop
      integer(int64) :: wall_start, wall_stop, clock_rate
      integer :: nbasis, nleft_bands, nright_bands, nproduct, npairs, nfrequency
      integer :: ie, ik, ifrequency, ib_left, ib_right, pair_index

      call system_clock(wall_start, clock_rate)
      call cpu_time(cpu_start)
      if (allocated(request%selected_l)) then
         call product_basis%component_vertex_tensor(vertices, request%selected_l)
      else
         call product_basis%component_vertex_tensor(vertices)
      end if
      nbasis = left_state%nbasis
      nleft_bands = left_state%nbands
      nright_bands = right_state%nbands
      nproduct = product_basis%product_dimension
      npairs = nleft_bands*nright_bands
      nfrequency = size(request%frequencies)
      energy_min = min(minval(left_state%eigenvalues), minval(right_state%eigenvalues)) - request%energy_margin
      energy_max = max(maxval(left_state%eigenvalues), maxval(right_state%eigenvalues)) + request%energy_margin
      if (energy_max <= energy_min) error stop 'evaluate_lr_product_gf_factorized: invalid integration interval'
      step = (energy_max - energy_min)/real(request%integration_points - 1, rp)
      weight_sum = sum(left_state%k_weights)
      spectral_prefactor = cmplx(0.0_rp, 1.0_rp/(2.0_rp*response_angular_pi), rp)

      allocate(result%susceptibility(nproduct, nproduct, nfrequency), result%frequencies(nfrequency), &
         transitions(nleft_bands, nright_bands, nproduct), transition_flat(npairs, nproduct), &
         kernel(npairs, nfrequency), weighted_transitions(npairs, nproduct), contribution(nproduct, nproduct), &
         left_spectral(nleft_bands), right_spectral(nright_bands), right_retarded(nright_bands), &
         left_advanced(nleft_bands))
      result%susceptibility = cmplx(0.0_rp, 0.0_rp, rp)

      do ik = 1, left_state%nk
         call build_lr_product_gf_transition_amplitudes(vertices, left_state, ik, right_state, ik, transitions)
         transition_flat = reshape(transitions, [npairs, nproduct])
         kernel = cmplx(0.0_rp, 0.0_rp, rp)
         do ie = 1, request%integration_points
            energy = energy_min + real(ie - 1, rp)*step
            if (ie == 1 .or. ie == request%integration_points) then
               quadrature_weight = 1.0_rp
            else if (mod(ie, 2) == 0) then
               quadrature_weight = 4.0_rp
            else
               quadrature_weight = 2.0_rp
            end if
            quadrature_weight = quadrature_weight*step/3.0_rp
            fermi_weight = lr_fermi_dirac_occupation(energy, left_state%fermi_level, left_state%temperature)
            do ib_left = 1, nleft_bands
               left_spectral(ib_left) = spectral_prefactor*( &
                  1.0_rp/(cmplx(energy, integration_eta, rp) - left_state%eigenvalues(ib_left, ik)) - &
                  1.0_rp/(cmplx(energy, -integration_eta, rp) - left_state%eigenvalues(ib_left, ik)))
            end do
            do ib_right = 1, nright_bands
               right_spectral(ib_right) = spectral_prefactor*( &
                  1.0_rp/(cmplx(energy, integration_eta, rp) - right_state%eigenvalues(ib_right, ik)) - &
                  1.0_rp/(cmplx(energy, -integration_eta, rp) - right_state%eigenvalues(ib_right, ik)))
            end do
            do ifrequency = 1, nfrequency
               do ib_right = 1, nright_bands
                  right_retarded(ib_right) = 1.0_rp/(cmplx(energy + request%frequencies(ifrequency), request%eta, rp) - &
                     right_state%eigenvalues(ib_right, ik))
               end do
               do ib_left = 1, nleft_bands
                  left_advanced(ib_left) = 1.0_rp/(cmplx(energy - request%frequencies(ifrequency), -request%eta, rp) - &
                     left_state%eigenvalues(ib_left, ik))
               end do
               do ib_right = 1, nright_bands
                  do ib_left = 1, nleft_bands
                     pair_index = ib_left + (ib_right - 1)*nleft_bands
                     kernel(pair_index, ifrequency) = kernel(pair_index, ifrequency) + &
                        quadrature_weight*fermi_weight*2.0_rp*left_state%k_weights(ik)/weight_sum*( &
                        left_spectral(ib_left)*right_retarded(ib_right) + &
                        right_spectral(ib_right)*left_advanced(ib_left))
                  end do
               end do
            end do
         end do
         do ifrequency = 1, nfrequency
            weighted_transitions = spread(kernel(:, ifrequency), 2, nproduct)*conjg(transition_flat)
            contribution = matmul(transpose(transition_flat), weighted_transitions)
            result%susceptibility(:, :, ifrequency) = result%susceptibility(:, :, ifrequency) + contribution
         end do
      end do

      call system_clock(wall_stop)
      call cpu_time(cpu_stop)
      result%q = request%q
      result%frequencies = request%frequencies
      result%eta = request%eta
      result%product_dimension = nproduct
      result%integration_points = request%integration_points
      result%actual_integration_eta = integration_eta
      result%energy_min = energy_min
      result%energy_max = energy_max
      result%energy_spacing = step
      result%spacing_over_integration_eta = step/integration_eta
      result%response_representation = product_representation
      result%contraction_backend = lr_product_gf_contraction_factorized
      result%point_response_allocated = .false.
      result%nbasis = nbasis
      result%nk = left_state%nk
      result%nfrequency = nfrequency
      result%component_vertex_memory_bytes = int(size(vertices), int64)*int(storage_size(vertices)/8, int64)
      result%gf_matrix_memory_bytes = 0_int64
      result%band_transition_memory_bytes = int(size(transition_flat), int64)*int(storage_size(transition_flat)/8, int64)
      result%band_kernel_memory_bytes = int(size(kernel) + size(weighted_transitions) + size(contribution) + &
         size(left_spectral) + size(right_spectral) + size(right_retarded) + size(left_advanced), int64)* &
         int(storage_size(kernel)/8, int64)
      result%susceptibility_memory_bytes = int(size(result%susceptibility), int64)* &
         int(storage_size(result%susceptibility)/8, int64)
      result%wall_time_seconds = real(max(0_int64, wall_stop - wall_start), rp)/real(max(1_int64, clock_rate), rp)
      result%cpu_time_seconds = max(0.0_rp, cpu_stop - cpu_start)
      result%time_per_energy_kpoint = result%wall_time_seconds/real(max(1, request%integration_points*left_state%nk), rp)
      if (trim(request%channel) == 'plus' .or. trim(request%channel) == 'chi_plus') then
         result%channel = lr_channel_plus
      else
         result%channel = lr_channel_minus
      end if
      write (result%response_space_metadata, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,l1,a,es12.4,a,es12.4,a,a,a,a,a,i0,a,i0)') &
         'product_dimension=', nproduct, ' component_vertex_bytes=', result%component_vertex_memory_bytes, &
         ' gf_matrix_bytes=', result%gf_matrix_memory_bytes, ' susceptibility_bytes=', result%susceptibility_memory_bytes, &
         ' energy_points=', request%integration_points, ' k_points=', left_state%nk, ' frequencies=', nfrequency, &
         ' point_response_allocated=', result%point_response_allocated, ' energy_spacing=', result%energy_spacing, &
         ' h_over_integration_eta=', result%spacing_over_integration_eta, ' contraction_backend=', &
         trim(result%contraction_backend), ' accepted_state=', 'shared_request_snapshots', ' band_pairs=', npairs, &
         ' band_transition_bytes=', result%band_transition_memory_bytes
   end subroutine evaluate_lr_product_gf_factorized

   !> Build T(I,n,m)=<L n|sum_pq eps_L^p eps_R^q V_I^pq|R m> from the
   !> existing compact GF component vertices.  This helper deliberately does
   !> not call the LR-06/Lehmann transition or susceptibility implementation.
   subroutine build_lr_product_gf_transition_amplitudes(vertices, left_state, left_ik, right_state, right_ik, transitions)
      complex(rp), intent(in) :: vertices(:, :, :, :)
      type(lr_electronic_state), intent(in) :: left_state, right_state
      integer, intent(in) :: left_ik, right_ik
      complex(rp), intent(out) :: transitions(:, :, :)

      complex(rp), allocatable :: right_projection(:, :), band_vertex(:, :)
      integer :: component, product_index, left_power, right_power, nbasis, nleft_bands, nright_bands
      integer :: ib_left, ib_right

      nbasis = size(vertices, 1)
      nleft_bands = left_state%nbands
      nright_bands = right_state%nbands
      if (size(vertices, 2) /= nbasis .or. size(vertices, 3) /= 4 .or. &
          size(vertices, 4) /= size(transitions, 3) .or. size(transitions, 1) /= nleft_bands .or. &
          size(transitions, 2) /= nright_bands .or. left_state%nbasis /= nbasis .or. &
          right_state%nbasis /= nbasis) then
         error stop 'build_lr_product_gf_transition_amplitudes: shape mismatch'
      end if
      if (left_ik < 1 .or. left_ik > left_state%nk .or. right_ik < 1 .or. right_ik > right_state%nk) then
         error stop 'build_lr_product_gf_transition_amplitudes: k-point index out of range'
      end if
      allocate(right_projection(nbasis, nright_bands), band_vertex(nleft_bands, nright_bands))
      transitions = cmplx(0.0_rp, 0.0_rp, rp)
      do product_index = 1, size(vertices, 4)
         do component = 1, 4
            left_power = mod(component - 1, 2)
            right_power = (component - 1)/2
            right_projection = matmul(vertices(:, :, component, product_index), &
               right_state%eigenvectors(:, :, right_ik))
            band_vertex = matmul(conjg(transpose(left_state%eigenvectors(:, :, left_ik))), right_projection)
            if (left_power == 1) then
               do ib_left = 1, nleft_bands
                  band_vertex(ib_left, :) = band_vertex(ib_left, :)*left_state%eigenvalues(ib_left, left_ik)
               end do
            end if
            if (right_power == 1) then
               do ib_right = 1, nright_bands
                  band_vertex(:, ib_right) = band_vertex(:, ib_right)*right_state%eigenvalues(ib_right, right_ik)
               end do
            end if
            transitions(:, :, product_index) = transitions(:, :, product_index) + band_vertex
         end do
      end do
      deallocate(right_projection, band_vertex)
   end subroutine build_lr_product_gf_transition_amplitudes

   !> Reference contraction retained for the TDVK-03 optimization oracle.
   !> It is intentionally close to the LR-GF-02 scalar implementation.
   subroutine accumulate_product_gf_bubble_scalar(vertices, left_a, right_a, right_gr, left_ga, scale, susceptibility)
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
         error stop 'accumulate_product_gf_bubble_scalar: product matrix shape mismatch'
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
   end subroutine accumulate_product_gf_bubble_scalar

   subroutine initialize_product_gf_contraction_workspace(workspace, vertices)
      type(lr_product_gf_contraction_workspace), intent(out) :: workspace
      complex(rp), intent(in) :: vertices(:, :, :, :)
      integer :: nb, nb2, product_dimension, component, i

      nb = size(vertices, 1)
      nb2 = nb*nb
      product_dimension = size(vertices, 4)
      allocate(workspace%vertices_flat(nb2, product_dimension, 4), &
         workspace%vertices_transpose(nb2, product_dimension, 4), workspace%transformed(nb2, product_dimension), &
         workspace%contribution(product_dimension, product_dimension), workspace%temporary(nb, nb), &
         workspace%active_component(4))
      do component = 1, 4
         workspace%vertices_flat(:, :, component) = reshape(vertices(:, :, component, :), [nb2, product_dimension])
         do i = 1, product_dimension
            workspace%vertices_transpose(:, i, component) = reshape(transpose(vertices(:, :, component, i)), [nb2])
         end do
         workspace%active_component(component) = maxval(abs(vertices(:, :, component, :))) /= 0.0_rp
      end do
   end subroutine initialize_product_gf_contraction_workspace

   !> Matrix-backed equivalent of accumulate_product_gf_bubble_scalar.
   !>
   !> For each pair of affine vertex components, the transformed vertices are
   !> flattened into columns.  The complete response block is then one dense
   !> matrix product rather than a product-dimension squared loop over scalar
   !> Frobenius contractions.  The two products below are exactly the two
   !> Kubo terms used by the scalar oracle; no Lehmann data or point response
   !> matrix is introduced here.
   subroutine accumulate_product_gf_bubble_optimized(vertices, left_a, right_a, right_gr, left_ga, scale, susceptibility, &
                                                     workspace)
      complex(rp), intent(in) :: vertices(:, :, :, :), left_a(:, :, :), right_a(:, :, :), &
         right_gr(:, :, :), left_ga(:, :, :)
      real(rp), intent(in) :: scale
      complex(rp), intent(inout) :: susceptibility(:, :)
      type(lr_product_gf_contraction_workspace), intent(inout) :: workspace

      integer :: component_i, component_j, left_power_i, right_power_i, left_power_j, right_power_j
      integer :: combined_left, combined_right, i, j, nb, nb2, product_dimension

      nb = size(left_a, 1)
      nb2 = nb*nb
      product_dimension = size(susceptibility, 1)
      if (size(left_a, 2) /= nb .or. size(right_a, 1) /= nb .or. size(right_a, 2) /= nb .or. &
          size(right_gr, 1) /= nb .or. size(right_gr, 2) /= nb .or. size(left_ga, 1) /= nb .or. &
          size(left_ga, 2) /= nb .or. size(susceptibility, 2) /= product_dimension .or. &
          size(vertices, 1) /= nb .or. size(vertices, 2) /= nb .or. size(vertices, 3) /= 4 .or. &
          size(vertices, 4) /= product_dimension) then
         error stop 'accumulate_product_gf_bubble_optimized: product matrix shape mismatch'
      end if
      if (.not. allocated(workspace%vertices_flat) .or. size(workspace%vertices_flat, 1) /= nb2 .or. &
          size(workspace%vertices_flat, 2) /= product_dimension .or. size(workspace%vertices_flat, 3) /= 4) then
         error stop 'accumulate_product_gf_bubble_optimized: workspace shape mismatch'
      end if

      do component_i = 1, 4
         left_power_i = mod(component_i - 1, 2)
         right_power_i = (component_i - 1)/2
         if (.not. workspace%active_component(component_i)) cycle
         do component_j = 1, 4
            left_power_j = mod(component_j - 1, 2)
            right_power_j = (component_j - 1)/2
            if (.not. workspace%active_component(component_j)) cycle
            combined_left = left_power_i + left_power_j + 1
            combined_right = right_power_i + right_power_j + 1

            ! First Kubo term, with column i equal to vec(A_L V_i G_R).
            do i = 1, product_dimension
               workspace%temporary = matmul(matmul(left_a(:, :, combined_left), vertices(:, :, component_i, i)), &
                  right_gr(:, :, combined_right))
               workspace%transformed(:, i) = reshape(workspace%temporary, [nb2])
            end do
            workspace%contribution = matmul(transpose(workspace%transformed), &
               conjg(workspace%vertices_flat(:, :, component_j)))
            susceptibility = susceptibility + scale*workspace%contribution

            ! Second Kubo term, with column j equal to
            ! vec(A_R V_j^dagger G_L^A) and row i equal to vec(V_i^T)^T.
            do j = 1, product_dimension
               workspace%temporary = matmul(matmul(right_a(:, :, combined_right), &
                  conjg(transpose(vertices(:, :, component_j, j)))), left_ga(:, :, combined_left))
               workspace%transformed(:, j) = reshape(workspace%temporary, [nb2])
            end do
            workspace%contribution = matmul(transpose(workspace%vertices_transpose(:, :, component_i)), &
               workspace%transformed)
            susceptibility = susceptibility + scale*workspace%contribution
         end do
      end do
   end subroutine accumulate_product_gf_bubble_optimized

   subroutine validate_product_gf_inputs(product_basis, left_state, right_state, request, channel_kind, contraction_kind)
      type(lmto_product_response_basis), intent(in) :: product_basis
      type(lr_electronic_state), intent(in) :: left_state, right_state
      type(lr_product_gf_susceptibility_request), intent(in) :: request
      integer, intent(out) :: channel_kind, contraction_kind

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
      select case (trim(request%contraction_backend))
      case (lr_product_gf_contraction_factorized)
         contraction_kind = 1
      case (lr_product_gf_contraction_optimized, 'blas', 'matrix')
         contraction_kind = 2
      case (lr_product_gf_contraction_scalar, 'oracle')
         contraction_kind = 3
      case default
         error stop 'evaluate_lr_product_gf_susceptibility: contraction_backend must be factorized, optimized or scalar'
      end select
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
