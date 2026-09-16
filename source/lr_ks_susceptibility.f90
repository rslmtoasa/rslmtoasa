!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Spectral/Lehmann collinear transverse KS susceptibility.
!>
!> The module owns the first bare-response implementation only.  It consumes
!> immutable snapshots of the left k and right k+q eigensystems, evaluates the
!> certified LR-05 Pauli transition vertex, and accumulates the LR-03 retarded
!> circular response.  It contains no XC kernel, Dyson solve, Goldstone repair,
!> mode fitting, or site projection.
!------------------------------------------------------------------------------
module lr_ks_susceptibility_mod

   use precision_mod, only: rp
   use reciprocal_mod, only: reciprocal
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_response_space_mod, only: response_space_layout, response_apply_operator, &
      response_vector_norm, response_raw_to_canonical
   use lr_pauli_transition_vertex_mod, only: pauli_endpoint_state, pauli_vertex_capabilities, &
      pauli_sigma_plus_matrix, pauli_sigma_minus_matrix, evaluate_pauli_transition_vertex
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, lmto_product_channel_plus, &
      lmto_product_channel_minus
   implicit none
   private

   real(rp), parameter, public :: lr_kb_ry_per_kelvin = 6.3336814e-6_rp
   real(rp), parameter :: occupation_kT_floor = 1.0e-10_rp
   real(rp), parameter :: state_match_tolerance = 2.0e-12_rp
   real(rp), parameter :: endpoint_match_tolerance = 2.0e-11_rp

   character(len=*), parameter, public :: lr_channel_plus = 'chi_plus'
   character(len=*), parameter, public :: lr_channel_minus = 'chi_minus'

   !> Immutable electronic-state snapshot consumed by LR-06.
   !>
   !> The initializer copies every array, including occupations and k weights.
   !> The susceptibility evaluator accepts this type with INTENT(IN), so it
   !> cannot recompute EF, occupations, or endpoint eigenvectors while summing.
   type, public :: lr_electronic_state
      integer :: nbands = 0
      integer :: nbasis = 0
      integer :: nk = 0
      real(rp), allocatable :: eigenvalues(:, :)       ! (band,k)
      complex(rp), allocatable :: eigenvectors(:, :, :) ! (basis,band,k)
      real(rp), allocatable :: k_points(:, :)          ! (3,k), folded points
      real(rp), allocatable :: k_weights(:)
      real(rp), allocatable :: occupations(:, :)       ! (band,k), explicit
      real(rp) :: fermi_level = 0.0_rp
      real(rp) :: temperature = 0.0_rp
      real(rp) :: energy_zero = 0.0_rp
      character(len=32) :: reciprocal_mode = 'ham_only'
      character(len=16) :: hamiltonian_order = 'second'
      logical :: orthogonal = .true.
      logical :: collinear = .true.
      logical :: has_soc = .false.
      logical :: has_extra_operator = .false.
      logical :: occupations_are_explicit = .false.
      logical :: initialized = .false.
   contains
      procedure :: initialize => lr_electronic_state_initialize
      procedure :: restore_to_default => lr_electronic_state_restore
      procedure :: validate => lr_electronic_state_validate
   end type lr_electronic_state

   !> Request metadata for one q/frequency/channel sweep.
   !>
   !> The pointer fields support the preferred request/result API.  An explicit
   !> overload is also provided for callers that keep these objects separately.
   type, public :: lr_ks_susceptibility_request
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = lr_channel_plus
      type(response_space_layout), pointer :: response_space => null()
      type(lmto_radial_basis), pointer :: radial_bases(:) => null()
      type(lr_electronic_state), pointer :: electronic_state => null()
      type(lr_electronic_state), pointer :: q_endpoint_state => null()
   end type lr_ks_susceptibility_request

   !> Result in the LR-04 canonical right-weighted representation.
   type, public :: lr_ks_susceptibility_result
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = ''
      character(len=128) :: response_representation = ''
      character(len=256) :: response_space_metadata = ''
      complex(rp), allocatable :: susceptibility(:, :, :) ! (I,J,frequency)
   end type lr_ks_susceptibility_result

   !> Request metadata for the compact weighted-orthonormal product response.
   !>
   !> This is deliberately a separate request type from LR-06.  A compact
   !> result must never be mistaken for an LR-04 point-space matrix whose
   !> dimension is `response_space%ndim`.
   type, public :: lr_product_ks_susceptibility_request
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = lr_channel_plus
      type(lmto_product_response_basis), pointer :: product_basis => null()
      type(lr_electronic_state), pointer :: electronic_state => null()
      type(lr_electronic_state), pointer :: q_endpoint_state => null()
      ! Optional DRESP-01 orbital selector.  An absent mask preserves the
      ! historical complete product-space oracle exactly.
      logical, allocatable :: selected_l(:)
   end type lr_product_ks_susceptibility_request

   !> Result in the weighted-orthonormal LMTO product representation.
   type, public :: lr_product_ks_susceptibility_result
      real(rp) :: q(3) = 0.0_rp
      real(rp), allocatable :: frequencies(:)
      real(rp) :: eta = 0.0_rp
      character(len=32) :: channel = ''
      character(len=128) :: response_representation = ''
      character(len=256) :: response_space_metadata = ''
      integer :: product_dimension = 0
      integer :: transition_dimension = 0
      integer :: ntransitions_evaluated = 0
      integer :: noccupation_skips = 0
      logical :: point_space_transition_allocated = .false.
      complex(rp), allocatable :: susceptibility(:, :, :) ! (product,product,frequency)
   end type lr_product_ks_susceptibility_result

   public :: lr_fermi_dirac_occupation
   public :: lr_snapshot_from_reciprocal
   public :: lr_q_endpoint_from_reciprocal
   public :: evaluate_lr_ks_susceptibility
   public :: evaluate_lr_product_ks_susceptibility
   public :: evaluate_lr_static_residual

   interface evaluate_lr_ks_susceptibility
      module procedure evaluate_lr_ks_susceptibility_request
      module procedure evaluate_lr_ks_susceptibility_explicit
   end interface evaluate_lr_ks_susceptibility

   interface evaluate_lr_product_ks_susceptibility
      module procedure evaluate_lr_product_ks_susceptibility_request
      module procedure evaluate_lr_product_ks_susceptibility_explicit
   end interface evaluate_lr_product_ks_susceptibility

   interface
      subroutine zgerc(m, n, alpha, x, incx, y, incy, a, lda)
         import :: rp
         integer, intent(in) :: m, n, incx, incy, lda
         complex(rp), intent(in) :: alpha, x(*), y(*)
         complex(rp), intent(inout) :: a(lda, *)
      end subroutine zgerc
   end interface

contains

   !> Numerically stable Fermi occupation.  This is an explicit helper only;
   !> LR-06 never calls it while evaluating a response snapshot.
   pure real(rp) function lr_fermi_dirac_occupation(eigenvalue, fermi_level, temperature) result(value)
      real(rp), intent(in) :: eigenvalue, fermi_level, temperature
      real(rp) :: argument, kT

      kT = max(temperature*lr_kb_ry_per_kelvin, occupation_kT_floor)
      argument = (eigenvalue - fermi_level)/kT
      if (argument >= 50.0_rp) then
         value = 0.0_rp
      else if (argument <= -50.0_rp) then
         value = 1.0_rp
      else
         value = 1.0_rp/(exp(argument) + 1.0_rp)
      end if
   end function lr_fermi_dirac_occupation

   subroutine lr_electronic_state_initialize(this, eigenvalues, eigenvectors, k_points, k_weights, occupations, &
                                              fermi_level, temperature, energy_zero, reciprocal_mode, &
                                              hamiltonian_order, orthogonal, collinear, has_soc, has_extra_operator)
      class(lr_electronic_state), intent(out) :: this
      real(rp), intent(in) :: eigenvalues(:, :), k_points(:, :), k_weights(:), occupations(:, :)
      complex(rp), intent(in) :: eigenvectors(:, :, :)
      real(rp), intent(in) :: fermi_level, temperature
      real(rp), intent(in), optional :: energy_zero
      character(len=*), intent(in), optional :: reciprocal_mode, hamiltonian_order
      logical, intent(in), optional :: orthogonal, collinear, has_soc, has_extra_operator

      if (size(eigenvalues, 1) < 1 .or. size(eigenvalues, 2) < 1 .or. size(eigenvectors, 1) < 1 .or. &
          size(eigenvectors, 2) /= size(eigenvalues, 1) .or. size(eigenvectors, 3) /= size(eigenvalues, 2) .or. &
          size(k_points, 1) /= 3 .or. size(k_points, 2) /= size(eigenvalues, 2) .or. &
          size(k_weights) /= size(eigenvalues, 2) .or. any(shape(occupations) /= shape(eigenvalues))) then
         error stop 'lr_electronic_state%initialize: inconsistent eigenpair, k-point, or occupation shapes'
      end if

      this%energy_zero = 0.0_rp
      this%reciprocal_mode = 'ham_only'
      this%hamiltonian_order = 'second'
      this%orthogonal = .true.
      this%collinear = .true.
      this%has_soc = .false.
      this%has_extra_operator = .false.
      this%occupations_are_explicit = .false.
      this%initialized = .false.
      this%nbands = size(eigenvalues, 1)
      this%nbasis = size(eigenvectors, 1)
      this%nk = size(eigenvalues, 2)
      allocate(this%eigenvalues(this%nbands, this%nk), this%eigenvectors(this%nbasis, this%nbands, this%nk), &
               this%k_points(3, this%nk), this%k_weights(this%nk), this%occupations(this%nbands, this%nk))
      this%eigenvalues = eigenvalues
      this%eigenvectors = eigenvectors
      this%k_points = k_points
      this%k_weights = k_weights
      this%occupations = occupations
      this%fermi_level = fermi_level
      this%temperature = temperature
      if (present(energy_zero)) this%energy_zero = energy_zero
      if (present(reciprocal_mode)) this%reciprocal_mode = trim(reciprocal_mode)
      if (present(hamiltonian_order)) this%hamiltonian_order = trim(hamiltonian_order)
      if (present(orthogonal)) this%orthogonal = orthogonal
      if (present(collinear)) this%collinear = collinear
      if (present(has_soc)) this%has_soc = has_soc
      if (present(has_extra_operator)) this%has_extra_operator = has_extra_operator
      this%occupations_are_explicit = .true.
      this%initialized = .true.
      call this%validate('lr_electronic_state%initialize')
   end subroutine lr_electronic_state_initialize

   subroutine lr_electronic_state_restore(this)
      class(lr_electronic_state), intent(inout) :: this

      if (allocated(this%eigenvalues)) deallocate(this%eigenvalues)
      if (allocated(this%eigenvectors)) deallocate(this%eigenvectors)
      if (allocated(this%k_points)) deallocate(this%k_points)
      if (allocated(this%k_weights)) deallocate(this%k_weights)
      if (allocated(this%occupations)) deallocate(this%occupations)
      this%nbands = 0
      this%nbasis = 0
      this%nk = 0
      this%fermi_level = 0.0_rp
      this%temperature = 0.0_rp
      this%energy_zero = 0.0_rp
      this%reciprocal_mode = 'ham_only'
      this%hamiltonian_order = 'second'
      this%orthogonal = .true.
      this%collinear = .true.
      this%has_soc = .false.
      this%has_extra_operator = .false.
      this%occupations_are_explicit = .false.
      this%initialized = .false.
   end subroutine lr_electronic_state_restore

   subroutine lr_electronic_state_validate(this, caller)
      class(lr_electronic_state), intent(in) :: this
      character(len=*), intent(in) :: caller
      real(rp) :: weight_sum

      if (.not. this%initialized .or. this%nbands < 1 .or. this%nbasis < 1 .or. this%nk < 1) then
         error stop trim(caller)//': electronic state is not initialized'
      end if
      if (.not. allocated(this%eigenvalues) .or. .not. allocated(this%eigenvectors) .or. &
          .not. allocated(this%k_points) .or. .not. allocated(this%k_weights) .or. &
          .not. allocated(this%occupations)) then
         error stop trim(caller)//': electronic state snapshot is incomplete'
      end if
      if (any(shape(this%eigenvalues) /= [this%nbands, this%nk]) .or. &
          any(shape(this%eigenvectors) /= [this%nbasis, this%nbands, this%nk]) .or. &
          any(shape(this%k_points) /= [3, this%nk]) .or. size(this%k_weights) /= this%nk .or. &
          any(shape(this%occupations) /= [this%nbands, this%nk])) then
         error stop trim(caller)//': electronic state snapshot has inconsistent shapes'
      end if
      if (.not. this%occupations_are_explicit) then
         error stop trim(caller)//': occupations must be an explicit immutable snapshot'
      end if
      if (.not. all(this%k_weights >= 0.0_rp)) then
         error stop trim(caller)//': k-point weights must be nonnegative'
      end if
      weight_sum = sum(this%k_weights)
      if (weight_sum <= tiny(1.0_rp)) then
         error stop trim(caller)//': k-point weights have zero total measure'
      end if
      if (minval(this%occupations) < -state_match_tolerance .or. &
          maxval(this%occupations) > 1.0_rp + state_match_tolerance) then
         error stop trim(caller)//': occupations must lie in [0,1]'
      end if
   end subroutine lr_electronic_state_validate

   !> Copy the accepted, already diagonalized reciprocal state into an
   !> immutable LR snapshot.  EF is read, never solved here.
   subroutine lr_snapshot_from_reciprocal(recip, state)
      type(reciprocal), intent(in) :: recip
      type(lr_electronic_state), intent(out) :: state
      real(rp), allocatable :: occupations(:, :)
      integer :: ib, ik

      if (.not. allocated(recip%eigenvalues) .or. .not. allocated(recip%eigenvectors)) then
         error stop 'lr_snapshot_from_reciprocal: complete reciprocal eigenpairs are required'
      end if
      if (.not. allocated(recip%k_workset%points) .or. .not. allocated(recip%k_workset%weights)) then
         error stop 'lr_snapshot_from_reciprocal: authoritative k-point workset is unavailable'
      end if
      if (recip%k_workset%distributed .or. recip%k_workset%nk_local /= recip%k_workset%nk_global) then
         error stop 'lr_snapshot_from_reciprocal: LR-06 baseline requires a replicated complete BZ workset'
      end if
      if (.not. recip%canonical_energy_valid) then
         error stop 'lr_snapshot_from_reciprocal: accepted canonical EF/occupation state is unavailable'
      end if
      if (size(recip%eigenvalues, 2) /= recip%k_workset%nk_local .or. &
          size(recip%eigenvectors, 2) /= size(recip%eigenvalues, 1) .or. &
          size(recip%eigenvectors, 3) /= recip%k_workset%nk_local) then
         error stop 'lr_snapshot_from_reciprocal: reciprocal eigenpair/workset shape mismatch'
      end if

      allocate(occupations(size(recip%eigenvalues, 1), size(recip%eigenvalues, 2)))
      do ik = 1, size(recip%eigenvalues, 2)
         do ib = 1, size(recip%eigenvalues, 1)
            occupations(ib, ik) = lr_fermi_dirac_occupation(recip%eigenvalues(ib, ik), recip%fermi_level, &
                                                            recip%temperature)
         end do
      end do
      call state%initialize(recip%eigenvalues, recip%eigenvectors, recip%k_workset%points, &
                            recip%k_workset%weights, occupations, recip%fermi_level, recip%temperature, 0.0_rp, &
                            recip%reciprocal_mode, recip%kspace_ham_order, .true., .true., recip%include_so, .false.)
   end subroutine lr_snapshot_from_reciprocal

   !> Produce the immutable k+q endpoint snapshot with the same EF, smearing,
   !> mode, and band count as LEFT_STATE.  The reciprocal service performs the
   !> exact folding and supplies its eigenvectors verbatim; this routine adds no
   !> endpoint phase or nearest-mesh approximation.
   subroutine lr_q_endpoint_from_reciprocal(recip, left_state, q, endpoint_state)
      type(reciprocal), intent(inout) :: recip
      type(lr_electronic_state), intent(in) :: left_state
      real(rp), intent(in) :: q(3)
      type(lr_electronic_state), intent(out) :: endpoint_state
      real(rp), allocatable :: target_k(:, :), folded_k(:, :), eigenvalues(:, :), occupations(:, :)
      complex(rp), allocatable :: eigenvectors(:, :, :)
      integer :: ib, ik

      call left_state%validate('lr_q_endpoint_from_reciprocal')
      if (trim(left_state%reciprocal_mode) /= trim(recip%reciprocal_mode) .or. &
          trim(left_state%hamiltonian_order) /= trim(recip%kspace_ham_order)) then
         error stop 'lr_q_endpoint_from_reciprocal: state/reciprocal representation provenance differs'
      end if
      allocate(target_k(3, left_state%nk))
      target_k = left_state%k_points + spread(q, dim=2, ncopies=left_state%nk)
      call recip%calculate_eigenpairs_at_kpoints(target_k, eigenvalues, eigenvectors, folded_k)
      if (.not. allocated(eigenvalues) .or. .not. allocated(eigenvectors) .or. .not. allocated(folded_k)) then
         error stop 'lr_q_endpoint_from_reciprocal: arbitrary-k eigensystem service returned no snapshot'
      end if
      if (any(shape(eigenvalues) /= [left_state%nbands, left_state%nk]) .or. &
          any(shape(eigenvectors) /= [left_state%nbasis, left_state%nbands, left_state%nk])) then
         error stop 'lr_q_endpoint_from_reciprocal: endpoint band shape differs from left state'
      end if
      allocate(occupations(left_state%nbands, left_state%nk))
      do ik = 1, left_state%nk
         do ib = 1, left_state%nbands
            occupations(ib, ik) = lr_fermi_dirac_occupation(eigenvalues(ib, ik), left_state%fermi_level, &
                                                            left_state%temperature)
         end do
      end do
      call endpoint_state%initialize(eigenvalues, eigenvectors, folded_k, left_state%k_weights, occupations, &
                                     left_state%fermi_level, left_state%temperature, left_state%energy_zero, &
                                     left_state%reciprocal_mode, left_state%hamiltonian_order, left_state%orthogonal, &
                                     left_state%collinear, left_state%has_soc, left_state%has_extra_operator)
   end subroutine lr_q_endpoint_from_reciprocal

   subroutine evaluate_lr_ks_susceptibility_request(request, result)
      type(lr_ks_susceptibility_request), intent(in) :: request
      type(lr_ks_susceptibility_result), intent(out) :: result

      if (.not. associated(request%response_space) .or. .not. associated(request%radial_bases) .or. &
          .not. associated(request%electronic_state) .or. .not. associated(request%q_endpoint_state)) then
         error stop 'evaluate_lr_ks_susceptibility: request references are incomplete'
      end if
      call evaluate_lr_ks_susceptibility_explicit(request%response_space, request%radial_bases, &
         request%electronic_state, request%q_endpoint_state, request, result)
   end subroutine evaluate_lr_ks_susceptibility_request

   subroutine evaluate_lr_ks_susceptibility_explicit(space, radial_bases, left_state, right_state, request, result)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(lr_electronic_state), intent(in) :: left_state, right_state
      type(lr_ks_susceptibility_request), intent(in) :: request
      type(lr_ks_susceptibility_result), intent(out) :: result

      type(pauli_vertex_capabilities) :: capabilities
      type(pauli_endpoint_state) :: left_band, right_band
      complex(rp), allocatable :: raw(:, :, :), transition(:)
      complex(rp) :: operator_matrix(2, 2), denominator, pair_factor
      real(rp) :: weight_sum, occupation_difference
      integer :: ik, ib, jb, ifrequency, i, j
      integer :: channel_kind

      call validate_susceptibility_inputs(space, radial_bases, left_state, right_state, request, channel_kind)
      call left_state%validate('evaluate_lr_ks_susceptibility:left_state')
      call right_state%validate('evaluate_lr_ks_susceptibility:q_endpoint_state')

      capabilities = pauli_vertex_capabilities()
      if (channel_kind == 1) then
         operator_matrix = pauli_sigma_plus_matrix()
      else
         operator_matrix = pauli_sigma_minus_matrix()
      end if

      allocate(raw(space%ndim, space%ndim, size(request%frequencies)), &
               result%susceptibility(space%ndim, space%ndim, size(request%frequencies)), transition(space%ndim))
      raw = cmplx(0.0_rp, 0.0_rp, rp)
      weight_sum = sum(left_state%k_weights)

      ! The full finite spectrum is intentional: no band window or energy
      ! truncation is applied in the certification backend.
      do ik = 1, left_state%nk
         do ib = 1, left_state%nbands
            call left_band%initialize(left_state%eigenvalues(ib, ik), left_state%eigenvectors(:, ib, ik))
            do jb = 1, right_state%nbands
               occupation_difference = left_state%occupations(ib, ik) - right_state%occupations(jb, ik)
               if (occupation_difference == 0.0_rp) cycle
               call right_band%initialize(right_state%eigenvalues(jb, ik), right_state%eigenvectors(:, jb, ik))
               call evaluate_pauli_transition_vertex(space, radial_bases, left_band, right_band, operator_matrix, &
                                                      capabilities, transition)
               do ifrequency = 1, size(request%frequencies)
                  denominator = cmplx(request%frequencies(ifrequency) + left_band%energy - right_band%energy, &
                                      request%eta, rp)
                  pair_factor = cmplx(2.0_rp*left_state%k_weights(ik)/weight_sum, 0.0_rp, rp)* &
                                occupation_difference/denominator
                  do j = 1, space%ndim
                     do i = 1, space%ndim
                        raw(i, j, ifrequency) = raw(i, j, ifrequency) + pair_factor*transition(i)*conjg(transition(j))
                     end do
                  end do
               end do
            end do
         end do
      end do

      do ifrequency = 1, size(request%frequencies)
         ! LR-04 owns the raw-to-canonical conversion.  The susceptibility
         ! itself remains a pointwise kernel until this one final conversion.
         call response_raw_to_canonical(space, raw(:, :, ifrequency), result%susceptibility(:, :, ifrequency))
      end do

      result%q = request%q
      allocate(result%frequencies(size(request%frequencies)))
      result%frequencies = request%frequencies
      result%eta = request%eta
      if (channel_kind == 1) then
         result%channel = lr_channel_plus
      else
         result%channel = lr_channel_minus
      end if
      result%response_representation = 'LR-04 canonical right-weighted B=chi_raw*W'
      write (result%response_space_metadata, '(a,i0,a,i0,a,i0,a,i0)') 'ndim=', space%ndim, &
         ' nsite=', space%nsite, ' radial_points=', space%npoint, ' channels=', space%nchannel
   end subroutine evaluate_lr_ks_susceptibility_explicit

   subroutine evaluate_lr_product_ks_susceptibility_request(request, result)
      type(lr_product_ks_susceptibility_request), intent(in) :: request
      type(lr_product_ks_susceptibility_result), intent(out) :: result

      if (.not. associated(request%product_basis) .or. .not. associated(request%electronic_state) .or. &
          .not. associated(request%q_endpoint_state)) then
         error stop 'evaluate_lr_product_ks_susceptibility: request references are incomplete'
      end if
      call evaluate_lr_product_ks_susceptibility_explicit(request%product_basis, request%electronic_state, &
         request%q_endpoint_state, request, result)
   end subroutine evaluate_lr_product_ks_susceptibility_request

   !> Accumulate the LR-06 Lehmann response directly in the weighted
   !> orthonormal LMTO product representation.  The point-grid LR-06 routine
   !> above remains the canonical legacy implementation and is intentionally
   !> not routed through this evaluator.
   subroutine evaluate_lr_product_ks_susceptibility_explicit(product_basis, left_state, right_state, request, result)
      type(lmto_product_response_basis), intent(in) :: product_basis
      type(lr_electronic_state), intent(in) :: left_state, right_state
      type(lr_product_ks_susceptibility_request), intent(in) :: request
      type(lr_product_ks_susceptibility_result), intent(out) :: result

      type(pauli_endpoint_state) :: left_band, right_band
      complex(rp), allocatable :: transition(:)
      complex(rp) :: denominator, pair_factor
      real(rp) :: weight_sum, occupation_difference
      integer :: ik, ib, jb, ifrequency, channel_kind

      call validate_product_susceptibility_inputs(product_basis, left_state, right_state, request, channel_kind)
      call left_state%validate('evaluate_lr_product_ks_susceptibility:left_state')
      call right_state%validate('evaluate_lr_product_ks_susceptibility:q_endpoint_state')

      allocate(result%susceptibility(product_basis%product_dimension, product_basis%product_dimension, &
         size(request%frequencies)), transition(product_basis%product_dimension))
      result%susceptibility = cmplx(0.0_rp, 0.0_rp, rp)
      result%q = request%q
      allocate(result%frequencies(size(request%frequencies)))
      result%frequencies = request%frequencies
      result%eta = request%eta
      result%product_dimension = product_basis%product_dimension
      result%transition_dimension = product_basis%product_dimension
      result%ntransitions_evaluated = 0
      result%noccupation_skips = 0
      result%point_space_transition_allocated = .false.
      result%response_representation = 'weighted-orthonormal LMTO product representation'
      write (result%response_space_metadata, '(a,i0,a,i0,a,i0,a)') &
         'product_dimension=', product_basis%product_dimension, ' unpruned_dimension=', &
         product_basis%unpruned_dimension, ' transition_dimension=', product_basis%product_dimension, &
         '; no point-space transition allocation'
      if (channel_kind == 1) then
         result%channel = lr_channel_plus
      else
         result%channel = lr_channel_minus
      end if

      weight_sum = sum(left_state%k_weights)
      ! This is the LR-06 full finite spectrum loop.  In particular, the
      ! occupation test, denominator, k-weight normalization, factor of two,
      ! and retarded eta are copied mechanically from the point evaluator.
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
               if (allocated(request%selected_l)) then
                  call product_basis%transition_coordinates(left_band, right_band, transition, request%selected_l)
               else
                  call product_basis%transition_coordinates(left_band, right_band, transition)
               end if
               result%ntransitions_evaluated = result%ntransitions_evaluated + 1
               do ifrequency = 1, size(request%frequencies)
                  denominator = cmplx(request%frequencies(ifrequency) + left_band%energy - right_band%energy, &
                     request%eta, rp)
                  pair_factor = cmplx(2.0_rp*left_state%k_weights(ik)/weight_sum, 0.0_rp, rp)* &
                     occupation_difference/denominator
                  call zgerc(product_basis%product_dimension, product_basis%product_dimension, pair_factor, transition, 1, &
                     transition, 1, result%susceptibility(:, :, ifrequency), product_basis%product_dimension)
               end do
            end do
         end do
      end do
   end subroutine evaluate_lr_product_ks_susceptibility_explicit

   subroutine validate_product_susceptibility_inputs(product_basis, left_state, right_state, request, channel_kind)
      type(lmto_product_response_basis), intent(in) :: product_basis
      type(lr_electronic_state), intent(in) :: left_state, right_state
      type(lr_product_ks_susceptibility_request), intent(in) :: request
      integer, intent(out) :: channel_kind
      real(rp) :: expected_k(3), scale
      integer :: ik, expected_nbasis

      if (.not. allocated(request%frequencies) .or. size(request%frequencies) < 1) then
         error stop 'evaluate_lr_product_ks_susceptibility: at least one frequency is required'
      end if
      if (request%eta <= 0.0_rp) error stop 'evaluate_lr_product_ks_susceptibility: retarded eta must be positive'
      select case (trim(request%channel))
      case ('plus', 'chi_plus')
         channel_kind = 1
      case ('minus', 'chi_minus')
         channel_kind = 2
      case default
         error stop 'evaluate_lr_product_ks_susceptibility: channel must be chi_plus or chi_minus'
      end select
      if (product_basis%circular_channel /= merge(lmto_product_channel_plus, lmto_product_channel_minus, channel_kind == 1)) then
         error stop 'evaluate_lr_product_ks_susceptibility: product/channel provenance differs'
      end if
      if (.not. allocated(product_basis%blocks) .or. product_basis%nsite < 1 .or. &
          product_basis%product_dimension < 1 .or. product_basis%response_lmax < 0) then
         error stop 'evaluate_lr_product_ks_susceptibility: product representation is not initialized'
      end if
      expected_nbasis = 2*(product_basis%orbital_lmax + 1)**2*product_basis%nsite
      if (left_state%nbasis /= expected_nbasis .or. right_state%nbasis /= expected_nbasis) then
         error stop 'evaluate_lr_product_ks_susceptibility: eigensystem basis is incompatible with product representation'
      end if
      if (left_state%nk /= right_state%nk .or. left_state%nbands /= right_state%nbands .or. &
          left_state%nbasis /= right_state%nbasis) then
         error stop 'evaluate_lr_product_ks_susceptibility: left and k+q endpoint dimensions differ'
      end if
      if (maxval(abs(left_state%k_weights - right_state%k_weights)) > state_match_tolerance) then
         error stop 'evaluate_lr_product_ks_susceptibility: k weights differ between endpoints'
      end if
      scale = max(1.0_rp, abs(left_state%fermi_level), abs(right_state%fermi_level), &
         abs(left_state%temperature), abs(right_state%temperature), abs(left_state%energy_zero), &
         abs(right_state%energy_zero))
      if (abs(left_state%fermi_level - right_state%fermi_level) > state_match_tolerance*scale .or. &
          abs(left_state%temperature - right_state%temperature) > state_match_tolerance*scale .or. &
          abs(left_state%energy_zero - right_state%energy_zero) > state_match_tolerance*scale) then
         error stop 'evaluate_lr_product_ks_susceptibility: EF/temperature/energy-zero provenance differs'
      end if
      if (trim(left_state%reciprocal_mode) /= 'ham_only' .or. trim(right_state%reciprocal_mode) /= 'ham_only' .or. &
          trim(left_state%hamiltonian_order) /= 'second' .or. trim(right_state%hamiltonian_order) /= 'second' .or. &
          .not. left_state%orthogonal .or. .not. right_state%orthogonal .or. .not. left_state%collinear .or. &
          .not. right_state%collinear .or. left_state%has_soc .or. right_state%has_soc .or. &
          left_state%has_extra_operator .or. right_state%has_extra_operator) then
         error stop 'evaluate_lr_product_ks_susceptibility: representation is outside the LR-06 baseline'
      end if
      do ik = 1, left_state%nk
         expected_k = fold_fractional_kpoint(left_state%k_points(:, ik) + request%q)
         if (maxval(abs(expected_k - right_state%k_points(:, ik))) > endpoint_match_tolerance) then
            error stop 'evaluate_lr_product_ks_susceptibility: k+q endpoint is not the exact folded endpoint'
         end if
      end do
   end subroutine validate_product_susceptibility_inputs

   !> Apply the raw q=0 static response to a separately supplied field and
   !> compare it with the separately supplied ground-state magnetization.
   !> This is diagnostic evidence only; no matrix or eigenvalue is modified.
   subroutine evaluate_lr_static_residual(space, static_susceptibility, field, magnetization, absolute_residual, &
                                           relative_residual)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: static_susceptibility(:, :), field(:), magnetization(:)
      real(rp), intent(out) :: absolute_residual, relative_residual
      complex(rp), allocatable :: response(:), residual(:)
      real(rp) :: target_norm

      if (any(shape(static_susceptibility) /= [space%ndim, space%ndim])) then
         error stop 'evaluate_lr_static_residual: susceptibility shape mismatch'
      end if
      if (size(field) /= space%ndim .or. size(magnetization) /= space%ndim) then
         error stop 'evaluate_lr_static_residual: vector shape mismatch'
      end if
      allocate(response(space%ndim), residual(space%ndim))
      call response_apply_operator(space, static_susceptibility, field, response)
      residual = response - magnetization
      absolute_residual = response_vector_norm(space, residual)
      target_norm = response_vector_norm(space, magnetization)
      if (target_norm > tiny(1.0_rp)) then
         relative_residual = absolute_residual/target_norm
      else
         relative_residual = absolute_residual
      end if
   end subroutine evaluate_lr_static_residual

   subroutine validate_susceptibility_inputs(space, radial_bases, left_state, right_state, request, channel_kind)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(lr_electronic_state), intent(in) :: left_state, right_state
      type(lr_ks_susceptibility_request), intent(in) :: request
      integer, intent(out) :: channel_kind
      real(rp) :: expected_k(3), scale
      integer :: ik, expected_nbasis

      if (.not. allocated(request%frequencies) .or. size(request%frequencies) < 1) then
         error stop 'evaluate_lr_ks_susceptibility: at least one frequency is required'
      end if
      if (request%eta <= 0.0_rp) error stop 'evaluate_lr_ks_susceptibility: retarded eta must be positive'
      select case (trim(request%channel))
      case ('plus', 'chi_plus')
         channel_kind = 1
      case ('minus', 'chi_minus')
         channel_kind = 2
      case default
         error stop 'evaluate_lr_ks_susceptibility: channel must be chi_plus or chi_minus'
      end select
      if (space%nchannel /= 1) then
         error stop 'evaluate_lr_ks_susceptibility: only one certified circular channel is selected per request'
      end if
      if (size(radial_bases) /= space%nsite) then
         error stop 'evaluate_lr_ks_susceptibility: one radial basis is required per response site'
      end if
      if (size(radial_bases) > 0) then
         expected_nbasis = 2*(radial_bases(1)%lmax + 1)**2*space%nsite
         if (left_state%nbasis /= expected_nbasis .or. right_state%nbasis /= expected_nbasis) then
            error stop 'evaluate_lr_ks_susceptibility: eigensystem basis is incompatible with LMTO radial basis'
         end if
      end if
      if (left_state%nk /= right_state%nk .or. left_state%nbands /= right_state%nbands .or. &
          left_state%nbasis /= right_state%nbasis) then
         error stop 'evaluate_lr_ks_susceptibility: left and k+q endpoint dimensions differ'
      end if
      scale = max(1.0_rp, abs(left_state%fermi_level), abs(right_state%fermi_level), &
                  abs(left_state%temperature), abs(right_state%temperature), abs(left_state%energy_zero), &
                  abs(right_state%energy_zero))
      if (abs(left_state%fermi_level - right_state%fermi_level) > state_match_tolerance*scale .or. &
          abs(left_state%temperature - right_state%temperature) > state_match_tolerance*scale .or. &
          abs(left_state%energy_zero - right_state%energy_zero) > state_match_tolerance*scale) then
         error stop 'evaluate_lr_ks_susceptibility: EF/temperature/energy-zero provenance differs'
      end if
      if (trim(left_state%reciprocal_mode) /= 'ham_only' .or. trim(right_state%reciprocal_mode) /= 'ham_only' .or. &
          trim(left_state%hamiltonian_order) /= 'second' .or. trim(right_state%hamiltonian_order) /= 'second' .or. &
          .not. left_state%orthogonal .or. .not. right_state%orthogonal .or. .not. left_state%collinear .or. &
          .not. right_state%collinear .or. left_state%has_soc .or. right_state%has_soc .or. &
          left_state%has_extra_operator .or. right_state%has_extra_operator) then
         error stop 'evaluate_lr_ks_susceptibility: representation is outside the LR-06 baseline'
      end if
      do ik = 1, left_state%nk
         expected_k = fold_fractional_kpoint(left_state%k_points(:, ik) + request%q)
         if (maxval(abs(expected_k - right_state%k_points(:, ik))) > endpoint_match_tolerance) then
            error stop 'evaluate_lr_ks_susceptibility: k+q endpoint is not the exact folded endpoint'
         end if
      end do
      call require_response_mesh(space, radial_bases)
   end subroutine validate_susceptibility_inputs

   subroutine require_response_mesh(space, radial_bases)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(pauli_vertex_capabilities) :: capabilities
      integer :: isite
      real(rp), parameter :: mesh_tolerance = 2.0e-13_rp

      capabilities = pauli_vertex_capabilities()
      if (space%ndim < 1 .or. .not. allocated(space%radius)) then
         error stop 'evaluate_lr_ks_susceptibility: response space is not initialized'
      end if
      do isite = 1, size(radial_bases)
         call radial_bases(isite)%require_supported(capabilities%reciprocal_mode, capabilities%hamiltonian_order, &
            capabilities%orthogonal, capabilities%collinear, capabilities%has_soc, capabilities%has_extra_operator)
         if (radial_bases(isite)%nspin /= 2 .or. radial_bases(isite)%npoint /= space%npoint .or. &
             .not. allocated(radial_bases(isite)%rofi) .or. .not. allocated(radial_bases(isite)%channel_present)) then
            error stop 'evaluate_lr_ks_susceptibility: incomplete radial basis provenance'
         end if
         if (size(radial_bases(isite)%rofi) /= space%npoint .or. &
             size(radial_bases(isite)%channel_present, 1) /= radial_bases(isite)%lmax + 1 .or. &
             size(radial_bases(isite)%channel_present, 2) /= 2 .or. &
             .not. all(radial_bases(isite)%channel_present)) then
            error stop 'evaluate_lr_ks_susceptibility: radial channel shape or completeness is invalid'
         end if
         if (abs(radial_bases(isite)%mesh_a - space%a) > mesh_tolerance*max(1.0_rp, abs(space%a)) .or. &
             abs(radial_bases(isite)%mesh_b - space%b) > mesh_tolerance*max(1.0_rp, abs(space%b)) .or. &
             maxval(abs(radial_bases(isite)%rofi - space%radius)) > mesh_tolerance* &
             max(1.0_rp, maxval(abs(space%radius)))) then
            error stop 'evaluate_lr_ks_susceptibility: radial mesh provenance differs from response space'
         end if
      end do
   end subroutine require_response_mesh

   pure function fold_fractional_kpoint(k_point) result(folded)
      real(rp), intent(in) :: k_point(3)
      real(rp) :: folded(3)

      folded = k_point - floor(k_point + 0.5_rp)
   end function fold_fractional_kpoint

end module lr_ks_susceptibility_mod
