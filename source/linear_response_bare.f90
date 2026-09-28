!------------------------------------------------------------------------------
! Linear-response bare-response services moved under linear_response_mod.
!------------------------------------------------------------------------------
submodule (linear_response_mod) linear_response_bare
   use, intrinsic :: iso_fortran_env, only: int64
   use precision_mod, only: rp
   use basis_mod, only: nb, spin_off
   use math_mod, only: i_unit
   use mpi_mod, only: rank, numprocs, get_mpi_variables
   use string_mod, only: lower
   use control_mod, only: control
   use lattice_mod, only: lattice
   use hamiltonian_mod, only: hamiltonian
   use reciprocal_mod, only: reciprocal
   use recursion_mod, only: recursion
   use green_mod, only: green
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   implicit none

   ! --- private state from lr_ks_susceptibility_mod ---
   real(rp), parameter :: occupation_kT_floor = 1.0e-10_rp
   real(rp), parameter :: state_match_tolerance = 2.0e-12_rp
   real(rp), parameter :: endpoint_match_tolerance = 2.0e-11_rp

   ! --- private state from lr_rs_gf_susceptibility_mod ---
   real(rp), parameter :: mesh_tolerance = 2.0e-13_rp

   ! --- private state from lr_projected_reciprocal_chi0_mod ---
   real(rp), parameter :: projected_state_match_tolerance = 2.0e-11_rp
   real(rp), parameter :: occupation_match_tolerance = 2.0e-10_rp
   integer, parameter :: n_dominant_transitions = 8


   interface
      subroutine zgerc(m, n, alpha, x, incx, y, incy, a, lda)
         import :: rp
         integer, intent(in) :: m, n, incx, incy, lda
         complex(rp), intent(in) :: alpha, x(*), y(*)
         complex(rp), intent(inout) :: a(lda, *)
      end subroutine zgerc
   end interface

contains

   ! Test-only projected-GF procedures remain in this submodule for the
   ! projected chi0 oracle.  Keep their eigenpair resolvent local so the
   ! production target has no dependency on the demoted GF oracle module.
   subroutine build_weighted_resolvent(state, ik, z, weighted_green)
      type(lr_electronic_state), intent(in) :: state
      integer, intent(in) :: ik
      complex(rp), intent(in) :: z
      complex(rp), intent(out) :: weighted_green(:, :, :)
      integer :: power, ib, i, j
      complex(rp) :: factor

      if (size(weighted_green, 1) /= state%nbasis .or. size(weighted_green, 2) /= state%nbasis .or. &
          size(weighted_green, 3) < 1 .or. size(weighted_green, 3) > lmto_product_max_gf_moment + 1) then
         error stop 'build_weighted_resolvent: output shape mismatch'
      end if
      weighted_green = cmplx(0.0_rp, 0.0_rp, rp)
      do power = 0, size(weighted_green, 3) - 1
         do ib = 1, state%nbands
            factor = cmplx(lmto_product_energy_power(state%eigenvalues(ib, ik), power), 0.0_rp, rp)/ &
                     (z - state%eigenvalues(ib, ik))
            do j = 1, state%nbasis
               do i = 1, state%nbasis
                  weighted_green(i, j, power + 1) = weighted_green(i, j, power + 1) + &
                     factor*state%eigenvectors(i, ib, ik)*conjg(state%eigenvectors(j, ib, ik))
               end do
            end do
         end do
      end do
   end subroutine build_weighted_resolvent

! --- from lr_ks_susceptibility_mod ---

   !> Numerically stable Fermi occupation.  This is an explicit helper only;
   !> LR-06 never calls it while evaluating a response snapshot.
   module pure real(rp) function lr_fermi_dirac_occupation(eigenvalue, fermi_level, temperature) result(value)
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

   module subroutine lr_electronic_state_initialize(this, eigenvalues, eigenvectors, k_points, k_weights, occupations, &
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

   module subroutine lr_electronic_state_restore(this)
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

   module subroutine lr_electronic_state_validate(this, caller)
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
   module subroutine lr_snapshot_from_reciprocal(recip, state)
      type(reciprocal), intent(in) :: recip
      type(lr_electronic_state), intent(out) :: state
      real(rp), allocatable :: occupations(:, :)
      integer :: ib, ik

      if (recip%k_workset%distributed .or. recip%k_workset%nk_local /= recip%k_workset%nk_global) then
         error stop 'lr_snapshot_from_reciprocal: LR-06 baseline requires a replicated complete BZ workset'
      end if
      if (.not. recip%canonical_energy_valid) then
         error stop 'lr_snapshot_from_reciprocal: accepted canonical EF/occupation state is unavailable'
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
   module subroutine lr_q_endpoint_from_reciprocal(recip, left_state, q, endpoint_state)
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

   module subroutine evaluate_lr_ks_susceptibility_request(request, result)
      type(lr_ks_susceptibility_request), intent(in) :: request
      type(lr_ks_susceptibility_result), intent(out) :: result

      call evaluate_lr_ks_susceptibility_explicit(request%response_space, request%radial_bases, &
         request%electronic_state, request%q_endpoint_state, request, result)
   end subroutine evaluate_lr_ks_susceptibility_request

   module subroutine evaluate_lr_ks_susceptibility_explicit(space, radial_bases, left_state, right_state, request, result)
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

   module subroutine evaluate_lr_product_ks_susceptibility_request(request, result)
      type(lr_product_ks_susceptibility_request), intent(in) :: request
      type(lr_product_ks_susceptibility_result), intent(out) :: result

      call evaluate_lr_product_ks_susceptibility_explicit(request%product_basis, request%electronic_state, &
         request%q_endpoint_state, request, result)
   end subroutine evaluate_lr_product_ks_susceptibility_request

   !> Accumulate the LR-06 Lehmann response directly in the weighted
   !> orthonormal LMTO product representation.  The point-grid LR-06 routine
   !> above remains the canonical legacy implementation and is intentionally
   !> not routed through this evaluator.
   module subroutine evaluate_lr_product_ks_susceptibility_explicit(product_basis, left_state, right_state, request, result)
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
      result%product_radial_order = product_basis%product_radial_order
      result%product_endpoint_branches = product_basis%product_endpoint_branches
      result%product_branch_labels = product_basis%product_branch_labels
      result%maximum_gf_energy_moment = lmto_product_max_gf_moment
      write (result%response_space_metadata, '(a,i0,a,i0,a,i0,a)') &
         'product_dimension=', product_basis%product_dimension, ' unpruned_dimension=', &
         product_basis%unpruned_dimension, ' transition_dimension=', product_basis%product_dimension, &
         '; no point-space transition allocation'
      result%response_space_metadata = trim(result%response_space_metadata)// &
         '; product_radial_order=second; product_endpoint_branches=6; '// &
         'product_branch_labels=00,10,01,11,20,02; maximum_gf_energy_moment=4'
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
   module subroutine evaluate_lr_static_residual(space, static_susceptibility, field, magnetization, absolute_residual, &
                                           relative_residual)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: static_susceptibility(:, :), field(:), magnetization(:)
      real(rp), intent(out) :: absolute_residual, relative_residual
      complex(rp), allocatable :: response(:), residual(:)
      real(rp) :: target_norm

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


! --- from lr_rs_gf_susceptibility_mod ---

   module subroutine lr_rs_block_recursion_initialize(this, callback, nbasis_per_site, recursion_depth, terminator)
      class(lr_rs_block_recursion_provider), intent(out) :: this
      procedure(lr_rs_native_callback) :: callback
      integer, intent(in) :: nbasis_per_site, recursion_depth
      character(len=*), intent(in), optional :: terminator

      if (nbasis_per_site < 1 .or. recursion_depth < 1) then
         error stop 'lr_rs_block_recursion_provider%initialize: invalid callback or recursion control'
      end if
      this%callback => callback
      this%nbasis_per_site = nbasis_per_site
      this%recursion_depth = recursion_depth
      this%provider_kind = 'block_recursion'
      this%terminator = 'certified native block terminator'
      if (present(terminator)) this%terminator = trim(terminator)
   end subroutine lr_rs_block_recursion_initialize

   module subroutine lr_rs_block_recursion_get_pair(this, left_site, right_site, translation, z, gij, gji, hgamma_ij, hgamma_ji)
      class(lr_rs_block_recursion_provider), intent(inout) :: this
      integer, intent(in) :: left_site, right_site
      real(rp), intent(in) :: translation(3)
      complex(rp), intent(in) :: z
      complex(rp), intent(out) :: gij(:, :), gji(:, :), hgamma_ij(:, :), hgamma_ji(:, :)

      call this%callback(left_site, right_site, translation, z, gij, gji, hgamma_ij, hgamma_ji)
   end subroutine lr_rs_block_recursion_get_pair

   module function lr_rs_block_recursion_describe(this) result(description)
      class(lr_rs_block_recursion_provider), intent(in) :: this
      character(len=256) :: description

      write (description, '(a,i0,a,a)') 'provider=block_recursion recursion_depth=', this%recursion_depth, &
         ' terminator='//trim(this%terminator)
   end function lr_rs_block_recursion_describe

   module subroutine lr_rs_chebyshev_initialize(this, callback, nbasis_per_site, polynomial_order, kernel)
      class(lr_rs_chebyshev_provider), intent(out) :: this
      procedure(lr_rs_native_callback) :: callback
      integer, intent(in) :: nbasis_per_site, polynomial_order
      character(len=*), intent(in), optional :: kernel

      if (nbasis_per_site < 1 .or. polynomial_order < 1) then
         error stop 'lr_rs_chebyshev_provider%initialize: invalid callback or polynomial control'
      end if
      this%callback => callback
      this%nbasis_per_site = nbasis_per_site
      this%polynomial_order = polynomial_order
      this%provider_kind = 'chebyshev'
      this%kernel = 'certified native Chebyshev kernel'
      if (present(kernel)) this%kernel = trim(kernel)
   end subroutine lr_rs_chebyshev_initialize

   module subroutine lr_rs_chebyshev_get_pair(this, left_site, right_site, translation, z, gij, gji, hgamma_ij, hgamma_ji)
      class(lr_rs_chebyshev_provider), intent(inout) :: this
      integer, intent(in) :: left_site, right_site
      real(rp), intent(in) :: translation(3)
      complex(rp), intent(in) :: z
      complex(rp), intent(out) :: gij(:, :), gji(:, :), hgamma_ij(:, :), hgamma_ji(:, :)

      call this%callback(left_site, right_site, translation, z, gij, gji, hgamma_ij, hgamma_ji)
   end subroutine lr_rs_chebyshev_get_pair

   module function lr_rs_chebyshev_describe(this) result(description)
      class(lr_rs_chebyshev_provider), intent(in) :: this
      character(len=256) :: description

      write (description, '(a,i0,a,a)') 'provider=chebyshev polynomial_order=', this%polynomial_order, &
         ' kernel='//trim(this%kernel)
   end function lr_rs_chebyshev_describe

   module subroutine evaluate_lr_rs_gf_susceptibility(request, result)
      type(lr_rs_gf_susceptibility_request), intent(in) :: request
      type(lr_ks_susceptibility_result), intent(out) :: result

      type(response_space_layout), pointer :: space
      type(lmto_radial_basis), pointer :: radial_bases(:)
      class(lr_rs_gf_provider), pointer :: provider
      type(pauli_vertex_capabilities) :: capabilities
      complex(rp), allocatable :: vertices(:, :, :, :), raw(:, :, :)
      complex(rp), allocatable :: gij(:, :), gji(:, :), hgamma_ij(:, :), hgamma_ji(:, :)
      complex(rp), allocatable :: gij_m(:, :), gji_m(:, :), hgamma_ij_m(:, :), hgamma_ji_m(:, :)
      complex(rp), allocatable :: gij_p(:, :), gji_p(:, :), hgamma_ij_p(:, :), hgamma_ji_p(:, :)
      complex(rp), allocatable :: gij_w(:, :), gji_w(:, :), hgamma_ij_w(:, :), hgamma_ji_w(:, :)
      complex(rp), allocatable :: gij_v(:, :), gji_v(:, :), hgamma_ij_v(:, :), hgamma_ji_v(:, :)
      type(lr_gf_augmented_block) :: ab_gr, ba_gr, ab_ga, ba_ga, ab_w, ba_w, ab_v, ba_v
      complex(rp) :: operator_matrix(2, 2), phase
      real(rp) :: integration_eta, energy_min, energy_max, step, energy, fermi_weight
      real(rp) :: quadrature_weight, scale
      real(rp), allocatable :: tau(:, :)
      complex(rp), allocatable :: pair_raw(:, :)
      integer :: channel_kind, nlocal, ne, ie, ipair, ifrequency

      call validate_rs_request(request, channel_kind)
      space => request%response_space
      radial_bases => request%radial_bases
      provider => request%provider
      capabilities = pauli_vertex_capabilities()
      if (channel_kind == 1) then
         operator_matrix = pauli_sigma_plus_matrix()
      else
         operator_matrix = pauli_sigma_minus_matrix()
      end if

      integration_eta = request%integration_eta
      if (integration_eta <= 0.0_rp) integration_eta = request%eta/40.0_rp
      if (integration_eta >= request%eta) then
         error stop 'evaluate_lr_rs_gf_susceptibility: integration_eta must be smaller than response eta'
      end if
      energy_min = request%energy_min - request%energy_margin
      energy_max = request%energy_max + request%energy_margin
      step = (energy_max - energy_min)/real(request%integration_points - 1, rp)
      ne = request%integration_points
      nlocal = provider%nbasis_per_site
      call build_native_vertex_tensor(space, radial_bases, operator_matrix, capabilities, vertices)
      allocate(raw(space%ndim, space%ndim, size(request%frequencies)))
      raw = cmplx(0.0_rp, 0.0_rp, rp)

      allocate(tau(3, space%nsite), pair_raw(space%ndim, space%ndim))
      tau = 0.0_rp
      if (associated(request%site_positions)) tau = request%site_positions

      allocate(gij(nlocal,nlocal), gji(nlocal,nlocal), hgamma_ij(nlocal,nlocal), hgamma_ji(nlocal,nlocal), &
         gij_m(nlocal,nlocal), gji_m(nlocal,nlocal), hgamma_ij_m(nlocal,nlocal), hgamma_ji_m(nlocal,nlocal), &
         gij_p(nlocal,nlocal), gji_p(nlocal,nlocal), hgamma_ij_p(nlocal,nlocal), hgamma_ji_p(nlocal,nlocal), &
         gij_w(nlocal,nlocal), gji_w(nlocal,nlocal), hgamma_ij_w(nlocal,nlocal), hgamma_ji_w(nlocal,nlocal), &
         gij_v(nlocal,nlocal), gji_v(nlocal,nlocal), hgamma_ij_v(nlocal,nlocal), hgamma_ji_v(nlocal,nlocal))

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
         fermi_weight = lr_fermi_dirac_occupation(energy, request%fermi_level, request%temperature)
         do ipair = 1, size(request%pairs)
            call provider%get_pair(request%pairs(ipair)%left_site, request%pairs(ipair)%right_site, &
               request%pairs(ipair)%translation, cmplx(energy, integration_eta, rp), gij_p, gji_p, hgamma_ij_p, hgamma_ji_p)
            call provider%get_pair(request%pairs(ipair)%left_site, request%pairs(ipair)%right_site, &
               request%pairs(ipair)%translation, cmplx(energy, -integration_eta, rp), gij_m, gji_m, hgamma_ij_m, hgamma_ji_m)
            call augment_lr_gf_endpoint_pair(radial_bases(request%pairs(ipair)%left_site), &
               radial_bases(request%pairs(ipair)%right_site), request%pairs(ipair)%left_site, &
               request%pairs(ipair)%right_site, cmplx(energy, integration_eta, rp), gij_p, gji_p, ab_gr, ba_gr, &
               hgamma_block=hgamma_ij_p, reverse_hgamma_block=hgamma_ji_p)
            call augment_lr_gf_endpoint_pair(radial_bases(request%pairs(ipair)%left_site), &
               radial_bases(request%pairs(ipair)%right_site), request%pairs(ipair)%left_site, &
               request%pairs(ipair)%right_site, cmplx(energy, -integration_eta, rp), gij_m, gji_m, ab_ga, ba_ga, &
               hgamma_block=hgamma_ij_m, reverse_hgamma_block=hgamma_ji_m)
            do ifrequency = 1, size(request%frequencies)
               call provider%get_pair(request%pairs(ipair)%left_site, request%pairs(ipair)%right_site, &
                  request%pairs(ipair)%translation, cmplx(energy + request%frequencies(ifrequency), request%eta, rp), &
                  gij_w, gji_w, hgamma_ij_w, hgamma_ji_w)
               call augment_lr_gf_endpoint_pair(radial_bases(request%pairs(ipair)%left_site), &
                  radial_bases(request%pairs(ipair)%right_site), request%pairs(ipair)%left_site, &
                  request%pairs(ipair)%right_site, cmplx(energy + request%frequencies(ifrequency), request%eta, rp), &
                  gij_w, gji_w, ab_w, ba_w, hgamma_block=hgamma_ij_w, reverse_hgamma_block=hgamma_ji_w)
               call provider%get_pair(request%pairs(ipair)%left_site, request%pairs(ipair)%right_site, &
                  request%pairs(ipair)%translation, cmplx(energy - request%frequencies(ifrequency), -request%eta, rp), &
                  gij_v, gji_v, hgamma_ij_v, hgamma_ji_v)
               call augment_lr_gf_endpoint_pair(radial_bases(request%pairs(ipair)%left_site), &
                  radial_bases(request%pairs(ipair)%right_site), request%pairs(ipair)%left_site, &
                  request%pairs(ipair)%right_site, cmplx(energy - request%frequencies(ifrequency), -request%eta, rp), &
                  gij_v, gji_v, ab_v, ba_v, hgamma_block=hgamma_ij_v, reverse_hgamma_block=hgamma_ji_v)

               scale = fermi_weight*quadrature_weight*2.0_rp
               pair_raw = cmplx(0.0_rp, 0.0_rp, rp)
               call accumulate_native_pair_bubble(space, vertices, ab_gr, ba_gr, ab_ga, ba_ga, ab_w, ba_v, &
                  request%pairs(ipair)%left_site, request%pairs(ipair)%right_site, scale, pair_raw)
               phase = response_real_space_phase(request%q, request%pairs(ipair)%translation, &
                  tau(:, request%pairs(ipair)%left_site), tau(:, request%pairs(ipair)%right_site))
               raw(:, :, ifrequency) = raw(:, :, ifrequency) + phase*pair_raw
            end do
         end do
      end do

      allocate(result%susceptibility(space%ndim, space%ndim, size(request%frequencies)), &
         result%frequencies(size(request%frequencies)))
      do ifrequency = 1, size(request%frequencies)
         call response_raw_to_canonical(space, raw(:, :, ifrequency), result%susceptibility(:, :, ifrequency))
      end do
      result%q = request%q
      result%frequencies = request%frequencies
      result%eta = request%eta
      if (channel_kind == 1) then
         result%channel = lr_channel_plus
      else
         result%channel = lr_channel_minus
      end if
      result%response_representation = 'LR-04 canonical right-weighted B=chi_raw*W'
      write (result%response_space_metadata, '(a,a,a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
         trim(provider%describe()), ' phase=+i2pi*q.(R+tau_right-tau_left)', &
         ' integration_points=', ne, ' integration_eta=', integration_eta, ' eta=', request%eta, &
         ' fermi_level=', request%fermi_level, ' temperature=', request%temperature
   end subroutine evaluate_lr_rs_gf_susceptibility

   subroutine accumulate_native_pair_bubble(space, vertices, ab_gr, ba_gr, ab_ga, ba_ga, ab_w, ba_v, &
                                            left_site, right_site, scale, raw)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: vertices(:, :, :, :)
      type(lr_gf_augmented_block), intent(in) :: ab_gr, ba_gr, ab_ga, ba_ga, ab_w, ba_v
      integer, intent(in) :: left_site, right_site
      real(rp), intent(in) :: scale
      complex(rp), intent(inout) :: raw(:, :)
      type(response_super_index) :: item_i, item_j
      integer :: i, j, p, q, r, s, component, ioffset, joffset, nlocal
      complex(rp), allocatable :: a_ba(:, :), a_ab(:, :), g_ab(:, :), g_ba(:, :), temporary(:, :)
      complex(rp), allocatable :: vi(:, :), vj(:, :), vi_sp(:, :), vj_rq(:, :)

      nlocal = ab_gr%norb*ab_gr%nspin
      allocate(a_ba(nlocal,nlocal), a_ab(nlocal,nlocal), g_ab(nlocal,nlocal), g_ba(nlocal,nlocal), &
         temporary(nlocal,nlocal), vi(nlocal,nlocal), vj(nlocal,nlocal), vi_sp(nlocal,nlocal), vj_rq(nlocal,nlocal))
      ioffset = (left_site - 1)*nlocal
      joffset = (right_site - 1)*nlocal
      do i = 1, space%ndim
         call response_unflatten_superindex(i, space%nsite, space%response_lmax, space%npoint, space%nchannel, item_i)
         if (item_i%site /= left_site) cycle
         do j = 1, space%ndim
            call response_unflatten_superindex(j, space%nsite, space%response_lmax, space%npoint, space%nchannel, item_j)
            if (item_j%site /= right_site) cycle
            do p = 0, 1
               do q = 0, 1
                  call native_branch(ba_gr, p, q, a_ba)
                  call native_branch(ba_ga, p, q, g_ba)
                  a_ba = cmplx(0.0_rp, 1.0_rp/(2.0_rp*response_angular_pi), rp)*(a_ba - g_ba)
                  call native_branch(ab_gr, p, q, g_ab)
                  call native_branch(ab_ga, p, q, a_ab)
                  a_ab = cmplx(0.0_rp, 1.0_rp/(2.0_rp*response_angular_pi), rp)*(g_ab - a_ab)
                  do r = 0, 1
                     do s = 0, 1
                        component = 1 + q + 2*r
                        vi = vertices(ioffset+1:ioffset+nlocal, ioffset+1:ioffset+nlocal, component, i)
                        component = 1 + p + 2*s
                        vj = vertices(joffset+1:joffset+nlocal, joffset+1:joffset+nlocal, component, j)
                        component = 1 + r + 2*s
                        call native_branch(ab_w, r, s, g_ab)
                        temporary = matmul(a_ba, vi)
                        temporary = matmul(temporary, g_ab)
                        raw(i, j) = raw(i, j) + scale*sum(temporary*conjg(vj))

                        component = 1 + q + 2*r
                        vj_rq = vertices(joffset+1:joffset+nlocal, joffset+1:joffset+nlocal, component, j)
                        component = 1 + s + 2*p
                        vi_sp = vertices(ioffset+1:ioffset+nlocal, ioffset+1:ioffset+nlocal, component, i)
                        call native_branch(ba_v, r, s, g_ba)
                        temporary = matmul(a_ab, conjg(transpose(vj_rq)))
                        temporary = matmul(temporary, g_ba)
                        raw(i, j) = raw(i, j) + scale*sum(temporary*transpose(vi_sp))
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine accumulate_native_pair_bubble

   subroutine native_branch(block, p, q, branch)
      type(lr_gf_augmented_block), intent(in) :: block
      integer, intent(in) :: p, q
      complex(rp), intent(out) :: branch(:, :)

      select case (1 + p + 2*q)
      case (1)
         branch = block%coefficient_gf
      case (2)
         branch = block%h_g
      case (3)
         branch = block%g_h
      case (4)
         branch = block%h_g_h
      case default
         error stop 'native_branch: invalid endpoint branch'
      end select
   end subroutine native_branch

   subroutine build_native_vertex_tensor(space, radial_bases, operator_matrix, capabilities, vertices)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: operator_matrix(2, 2)
      type(pauli_vertex_capabilities), intent(in) :: capabilities
      complex(rp), allocatable, intent(out) :: vertices(:, :, :, :)
      type(response_super_index) :: item
      integer :: nbasis, norb, flat, p, q, component, isite, iorb, jorb, ispin, jspin
      integer :: orbital_l, orbital_lp, offset
      real(rp) :: radial_product

      norb = (radial_bases(1)%lmax + 1)**2
      nbasis = 2*norb*space%nsite
      allocate(vertices(nbasis, nbasis, 4, space%ndim))
      vertices = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, space%npoint, &
            space%nchannel, item)
         isite = item%site
         offset = (isite - 1)*2*norb
         do p = 0, 1
            do q = 0, 1
               component = 1 + p + 2*q
               do iorb = 1, norb
                  orbital_l = lmto_orbital_l(iorb)
                  do jorb = 1, norb
                     orbital_lp = lmto_orbital_l(jorb)
                     do ispin = 1, 2
                        do jspin = 1, 2
                           radial_product = native_radial_vertex_component(radial_bases(isite), item%radial_point, &
                              orbital_l, orbital_lp, ispin, jspin, p, q)
                           vertices(offset+(ispin-1)*norb+iorb, offset+(jspin-1)*norb+jorb, component, flat) = &
                              operator_matrix(ispin,jspin)*radial_product*response_gaunt(orbital_l, orbital_m(iorb, orbital_l), &
                              orbital_lp, orbital_m(jorb, orbital_lp), item%response_l, item%response_m)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      call radial_bases(1)%require_supported(capabilities%reciprocal_mode, capabilities%hamiltonian_order, &
         capabilities%orthogonal, capabilities%collinear, capabilities%has_soc, capabilities%has_extra_operator)
   end subroutine build_native_vertex_tensor

   function native_radial_vertex_component(basis, radial_point, l, lp, ispin, jspin, p, q) result(value)
      type(lmto_radial_basis), intent(in) :: basis
      integer, intent(in) :: radial_point, l, lp, ispin, jspin, p, q
      real(rp) :: value
      real(rp) :: first, second, rfirst, rsecond

      if (radial_point /= 1) then
         if (basis%rofi(radial_point) <= tiny(1.0_rp)) error stop 'native radial vertex: invalid positive radial point'
         value = native_radial_component(basis, radial_point, l, ispin, p)* &
            native_radial_component(basis, radial_point, lp, jspin, q)/(basis%rofi(radial_point)**2)
         return
      end if
      if (l /= 0 .or. lp /= 0) then
         value = 0.0_rp
         return
      end if
      rfirst = basis%rofi(2)
      rsecond = basis%rofi(3)
      first = native_radial_component(basis, 2, l, ispin, p)*native_radial_component(basis, 2, lp, jspin, q)/(rfirst**2)
      second = native_radial_component(basis, 3, l, ispin, p)*native_radial_component(basis, 3, lp, jspin, q)/(rsecond**2)
      value = (first*rsecond**2 - second*rfirst**2)/(rsecond**2 - rfirst**2)
   end function native_radial_vertex_component

   pure real(rp) function native_radial_component(basis, ir, l, ispin, branch) result(value)
      type(lmto_radial_basis), intent(in) :: basis
      integer, intent(in) :: ir, l, ispin, branch

      if (branch == 0) then
         value = basis%phi_large(ir, l + 1, ispin)
      else
         value = basis%phidot_large(ir, l + 1, ispin)
      end if
   end function native_radial_component

   subroutine validate_rs_request(request, channel_kind)
      type(lr_rs_gf_susceptibility_request), intent(in) :: request
      integer, intent(out) :: channel_kind
      type(response_space_layout), pointer :: space
      type(lmto_radial_basis), pointer :: radial_bases(:)
      type(pauli_vertex_capabilities) :: capabilities
      integer :: isite, ipair, expected_nlocal

      if (.not. associated(request%response_space) .or. .not. associated(request%radial_bases) .or. &
          .not. associated(request%provider)) error stop 'evaluate_lr_rs_gf_susceptibility: request references are incomplete'
      if (.not. allocated(request%frequencies) .or. size(request%frequencies) < 1) then
         error stop 'evaluate_lr_rs_gf_susceptibility: at least one frequency is required'
      end if
      if (request%eta <= 0.0_rp .or. request%integration_points < 3 .or. &
          mod(request%integration_points, 2) == 0) then
         error stop 'evaluate_lr_rs_gf_susceptibility: invalid eta or Simpson integration_points'
      end if
      if (request%energy_margin < 0.0_rp .or. request%energy_max <= request%energy_min) then
         error stop 'evaluate_lr_rs_gf_susceptibility: invalid energy integration window'
      end if
      select case (trim(request%channel))
      case ('plus', 'chi_plus')
         channel_kind = 1
      case ('minus', 'chi_minus')
         channel_kind = 2
      case default
         error stop 'evaluate_lr_rs_gf_susceptibility: channel must be chi_plus or chi_minus'
      end select
      space => request%response_space
      radial_bases => request%radial_bases
      if (space%nchannel /= 1 .or. size(radial_bases) /= space%nsite) then
         error stop 'evaluate_lr_rs_gf_susceptibility: response channel/site mismatch'
      end if
      expected_nlocal = 2*(radial_bases(1)%lmax + 1)**2
      if (request%provider%nbasis_per_site /= expected_nlocal) then
         error stop 'evaluate_lr_rs_gf_susceptibility: provider coefficient block size differs from radial basis'
      end if
      capabilities = pauli_vertex_capabilities()
      do isite = 1, space%nsite
         call radial_bases(isite)%require_supported(capabilities%reciprocal_mode, capabilities%hamiltonian_order, &
            capabilities%orthogonal, capabilities%collinear, capabilities%has_soc, capabilities%has_extra_operator)
         if (radial_bases(isite)%lmax /= radial_bases(1)%lmax .or. radial_bases(isite)%nspin /= 2 .or. &
             radial_bases(isite)%npoint /= space%npoint .or. .not. allocated(radial_bases(isite)%rofi) .or. &
             .not. allocated(radial_bases(isite)%channel_present) .or. .not. all(radial_bases(isite)%channel_present)) then
            error stop 'evaluate_lr_rs_gf_susceptibility: incomplete radial basis provenance'
         end if
         if (abs(radial_bases(isite)%mesh_a - space%a) > mesh_tolerance*max(1.0_rp,abs(space%a)) .or. &
             abs(radial_bases(isite)%mesh_b - space%b) > mesh_tolerance*max(1.0_rp,abs(space%b)) .or. &
             maxval(abs(radial_bases(isite)%rofi - space%radius)) > mesh_tolerance*max(1.0_rp,maxval(abs(space%radius)))) then
            error stop 'evaluate_lr_rs_gf_susceptibility: radial mesh provenance differs from response space'
         end if
      end do
      if (.not. allocated(request%pairs) .or. size(request%pairs) < 1) then
         error stop 'evaluate_lr_rs_gf_susceptibility: at least one real-space pair is required'
      end if
      do ipair = 1, size(request%pairs)
         if (request%pairs(ipair)%left_site < 1 .or. request%pairs(ipair)%left_site > space%nsite .or. &
             request%pairs(ipair)%right_site < 1 .or. request%pairs(ipair)%right_site > space%nsite) then
            error stop 'evaluate_lr_rs_gf_susceptibility: pair site is outside the response space'
         end if
      end do
      if (associated(request%site_positions)) then
         if (any(shape(request%site_positions) /= [3, space%nsite])) then
            error stop 'evaluate_lr_rs_gf_susceptibility: site-position shape mismatch'
         end if
      end if
   end subroutine validate_rs_request

   pure integer function orbital_m(iorb, l) result(m)
      integer, intent(in) :: iorb, l
      m = iorb - l*l - l - 1
   end function orbital_m

! --- from tddft_native_rsgf_provider_mod ---

   module subroutine tddft_native_rsgf_initialize(this, provider_name, green_obj, recursion_obj, hamiltonian_obj, &
                                            lattice_obj, reciprocal_obj, radial_bases)
      class(tddft_native_rsgf_provider), intent(out) :: this
      character(len=*), intent(in) :: provider_name
      type(green), target, intent(inout) :: green_obj
      type(recursion), target, intent(inout) :: recursion_obj
      type(hamiltonian), target, intent(in) :: hamiltonian_obj
      type(lattice), target, intent(in) :: lattice_obj
      type(reciprocal), target, intent(in) :: reciprocal_obj
      type(lmto_radial_basis), target, intent(in) :: radial_bases(:)

      integer :: ipair, nsite, left_atom, right_atom, nw, recursion_slot
      real(rp) :: a_inf0(4), b_inf0(4)
      character(len=32) :: selected
      logical :: onsite_only

      if (numprocs /= 1) then
         error stop 'native RSGF production provider: the registered baseline requires a serial replicated pair workset'
      end if
      if (lattice_obj%njij < 1 .or. .not. allocated(lattice_obj%ijpair)) then
         error stop 'native RSGF production provider: a complete lattice%ijpair set is required; site-only fallback is forbidden'
      end if
      if (size(radial_bases) < 1) error stop 'native RSGF production provider: no response-site radial bases were supplied'

      selected = trim(lower(provider_name))
      if (selected == 'auto') selected = trim(lower(recursion_obj%control%recur))
      select case (selected)
      case ('block', 'block_recursion')
         if (trim(lower(recursion_obj%control%recur)) /= 'block') then
            error stop 'native RSGF production provider: block provider requires control%recur=block'
         end if
         this%selection = 'block'
         this%provider_kind = 'block_recursion'
         this%recursion_depth = recursion_obj%control%lld
         if (this%recursion_depth < 1) error stop 'native RSGF production provider: invalid block recursion depth'
      case ('chebyshev')
         if (trim(lower(recursion_obj%control%recur)) /= 'chebyshev') then
            error stop 'native RSGF production provider: Chebyshev provider requires control%recur=chebyshev'
         end if
         this%selection = 'chebyshev'
         this%provider_kind = 'chebyshev'
         this%polynomial_order = recursion_obj%control%lld
         if (this%polynomial_order < 1) error stop 'native RSGF production provider: invalid Chebyshev polynomial order'
      case default
         error stop 'native RSGF production provider: use block or chebyshev'
      end select

      this%green_obj => green_obj
      this%recursion_obj => recursion_obj
      this%hamiltonian_obj => hamiltonian_obj
      this%lattice_obj => lattice_obj
      this%reciprocal_obj => reciprocal_obj
      this%radial_bases => radial_bases
      nsite = size(radial_bases)
      this%nbasis_per_site = 2*(radial_bases(1)%lmax + 1)**2
      if (size(lattice_obj%ijpair, 1) < lattice_obj%njij .or. size(lattice_obj%ijpair, 2) /= 2) then
         error stop 'native RSGF production provider: lattice%ijpair shape is incomplete'
      end if

      allocate(this%pairs(lattice_obj%njij), this%left_atoms(lattice_obj%njij), this%right_atoms(lattice_obj%njij), &
         this%hgamma_ij(this%nbasis_per_site, this%nbasis_per_site, lattice_obj%njij), &
         this%hgamma_ji(this%nbasis_per_site, this%nbasis_per_site, lattice_obj%njij))
      allocate(this%identity(nb, nb))
      this%identity = cmplx(0.0_rp, 0.0_rp, rp)
      do ipair = 1, nb
         this%identity(ipair, ipair) = cmplx(1.0_rp, 0.0_rp, rp)
      end do

      ! Pair recursion uses the same global pair partition as the existing
      ! transport/real-space service.  The serial gate above makes the local
      ! pair index identical to the request pair index.
      call get_mpi_variables(rank, lattice_obj%njij)
      onsite_only = all(lattice_obj%ijpair(1:lattice_obj%njij, 1) == &
                        lattice_obj%ijpair(1:lattice_obj%njij, 2))
      select case (this%selection)
      case ('block')
         if (.not. onsite_only) then
            call recursion_obj%recur_b_ij()
            call recursion_obj%zsqr()
         else if (.not. recursion_obj%b2_b_is_sqrt) then
            ! The accepted SCF state already ran recur_b()/zsqr() for its
            ! response atoms.  Reuse those coefficients for an all-on-site
            ! pair workset instead of repeating the same block recursion.
            call recursion_obj%zsqr()
         end if
         allocate(this%a_inf(nb, nb, 4, lattice_obj%njij), this%b_inf(nb, nb, 4, lattice_obj%njij))
      case ('chebyshev')
         call recursion_obj%chebyshev_recur_ij()
      end select

      do ipair = 1, lattice_obj%njij
         left_atom = lattice_obj%ijpair(ipair, 1)
         right_atom = lattice_obj%ijpair(ipair, 2)
         call validate_atom_index(lattice_obj, left_atom, 'left')
         call validate_atom_index(lattice_obj, right_atom, 'right')
         this%left_atoms(ipair) = left_atom
         this%right_atoms(ipair) = right_atom
         this%pairs(ipair)%left_site = response_site_for_atom(lattice_obj, left_atom, nsite)
         this%pairs(ipair)%right_site = response_site_for_atom(lattice_obj, right_atom, nsite)
         call pair_translation(reciprocal_obj, lattice_obj, left_atom, right_atom, this%pairs(ipair)%translation)
         call native_h_eff_block(this, right_atom, left_atom, this%hgamma_ij(:, :, ipair))
         call native_h_eff_block(this, left_atom, right_atom, this%hgamma_ji(:, :, ipair))
         if (left_atom == right_atom) then
            call remove_radial_enu(this%hgamma_ij(:, :, ipair), this%radial_bases(this%pairs(ipair)%left_site))
            call remove_radial_enu(this%hgamma_ji(:, :, ipair), this%radial_bases(this%pairs(ipair)%right_site))
         end if
         if (this%selection == 'block') then
            nw = 10*this%recursion_depth
            if (onsite_only) then
               recursion_slot = recursion_site_slot(lattice_obj, left_atom)
               if (recursion_slot < 1) then
                  error stop 'native RSGF production provider: on-site pair is not an accepted recursion atom'
               end if
               call recursion_obj%get_terminf(recursion_obj%a_b(:, :, :, recursion_slot), &
                  recursion_obj%b2_b(:, :, :, recursion_slot), 1, this%recursion_depth, nb, nw, &
                  this%a_inf(:, :, 1:1, ipair), this%b_inf(:, :, 1:1, ipair), a_inf0(1:1), b_inf0(1:1))
               this%a_inf(:, :, 2:4, ipair) = 0.0_rp
               this%b_inf(:, :, 2:4, ipair) = 0.0_rp
            else
               call recursion_obj%get_terminf(recursion_obj%a_b(:, :, :, (ipair - 1)*4 + 1), &
                  recursion_obj%b2_b(:, :, :, (ipair - 1)*4 + 1), 4, this%recursion_depth, nb, nw, &
                  this%a_inf(:, :, :, ipair), this%b_inf(:, :, :, ipair), a_inf0, b_inf0)
            end if
         end if
      end do
   end subroutine tddft_native_rsgf_initialize

   module subroutine tddft_native_rsgf_get_pair(this, left_site, right_site, translation, z, gij, gji, hgamma_ij, hgamma_ji)
      class(tddft_native_rsgf_provider), intent(inout) :: this
      integer, intent(in) :: left_site, right_site
      real(rp), intent(in) :: translation(3)
      complex(rp), intent(in) :: z
      complex(rp), intent(out) :: gij(:, :), gji(:, :), hgamma_ij(:, :), hgamma_ji(:, :)

      complex(rp) :: gphase(nb, nb, 4), g_one(nb, nb, 1), z_grid(1)
      integer :: ipair, phase

      ipair = find_pair(this, left_site, right_site, translation)
      hgamma_ij = this%hgamma_ij(:, :, ipair)
      hgamma_ji = this%hgamma_ji(:, :, ipair)
      z_grid(1) = z

      select case (this%selection)
      case ('block')
         gphase = cmplx(0.0_rp, 0.0_rp, rp)
         if (this%left_atoms(ipair) == this%right_atoms(ipair)) then
            call this%green_obj%bgreen_complex(g_one, (ipair - 1)*4 + 1, z_grid, &
               this%a_inf(:, :, 1, ipair), this%b_inf(:, :, 1, ipair), .false.)
            gphase(:, :, 1) = g_one(:, :, 1)
         else
            do phase = 1, 4
               g_one = cmplx(0.0_rp, 0.0_rp, rp)
               call this%green_obj%bgreen_complex(g_one, (ipair - 1)*4 + phase, z_grid, &
                  this%a_inf(:, :, phase, ipair), this%b_inf(:, :, phase, ipair), .false.)
               gphase(:, :, phase) = g_one(:, :, 1)
            end do
         end if
      case ('chebyshev')
         call this%green_obj%chebyshev_green_complex_ij((ipair - 1)*4 + 1, z, gphase)
      end select

      if (this%left_atoms(ipair) == this%right_atoms(ipair)) then
         gij = gphase(:, :, 1)
         gji = gij
      else
         gij = gphase(:, :, 1) - gphase(:, :, 2) + (gphase(:, :, 3)/i_unit - gphase(:, :, 4)/i_unit)
         gji = gphase(:, :, 1) - gphase(:, :, 2) - (gphase(:, :, 3)/i_unit - gphase(:, :, 4)/i_unit)
         gij = 0.5_rp*gij
         gji = 0.5_rp*gji
      end if
   end subroutine tddft_native_rsgf_get_pair

   module function tddft_native_rsgf_describe(this) result(description)
      class(tddft_native_rsgf_provider), intent(in) :: this
      character(len=256) :: description

      if (this%selection == 'block') then
         write(description, '(a,a,a,i0,a,i0,a)') 'provider=native_block_recursion selection=', trim(this%selection), &
            ' recursion_depth=', this%recursion_depth, ' pairs=', size(this%pairs), ' workset=serial_replicated'
      else
         write(description, '(a,a,a,i0,a,i0,a)') 'provider=native_chebyshev selection=', trim(this%selection), &
            ' polynomial_order=', this%polynomial_order, ' pairs=', size(this%pairs), ' workset=serial_replicated'
      end if
   end function tddft_native_rsgf_describe

   subroutine native_h_eff_block(this, source_atom, destination_atom, block)
      class(tddft_native_rsgf_provider), intent(inout) :: this
      integer, intent(in) :: source_atom, destination_atom
      complex(rp), intent(out) :: block(:, :)
      integer :: destination_type, intermediate_atom, intermediate_type, direct_slot, intermediate_slot, ineigh
      integer :: neighbor_count

      if (source_atom < 1 .or. source_atom > this%lattice_obj%kk .or. &
          destination_atom < 1 .or. destination_atom > this%lattice_obj%kk) then
         error stop 'native RSGF production provider: H_eff action atom is outside the cluster'
      end if
      destination_type = this%lattice_obj%iz(destination_atom)
      block = cmplx(0.0_rp, 0.0_rp, rp)
      direct_slot = direct_hopping_slot(this%lattice_obj, destination_atom, source_atom)
      if (direct_slot > 0) block = this%hamiltonian_obj%ee(:, :, direct_slot, destination_type)

      ! H_eff = h - h*o*h + E_nu, written in the same directed-bond
      ! convention used by recursion%ham_hoh_vec_matmul and H(k) assembly.
      neighbor_count = this%lattice_obj%nn(destination_atom, 1)
      do ineigh = 1, neighbor_count
         if (ineigh == 1) then
            intermediate_atom = destination_atom
         else
            intermediate_atom = this%lattice_obj%nn(destination_atom, ineigh)
         end if
         if (intermediate_atom < 1 .or. intermediate_atom > this%lattice_obj%kk) cycle
         intermediate_type = this%lattice_obj%iz(intermediate_atom)
         intermediate_slot = direct_hopping_slot(this%lattice_obj, intermediate_atom, source_atom)
         if (intermediate_slot > 0) then
            block = block - matmul(this%hamiltonian_obj%eeo(:, :, ineigh, destination_type), &
               this%hamiltonian_obj%ee(:, :, intermediate_slot, intermediate_type))
         end if
      end do
      if (destination_atom == source_atom) block = block + this%hamiltonian_obj%enim(:, :, destination_type)
   end subroutine native_h_eff_block

   integer function direct_hopping_slot(lattice_obj, destination_atom, source_atom) result(slot)
      type(lattice), intent(in) :: lattice_obj
      integer, intent(in) :: destination_atom, source_atom
      integer :: ineigh

      slot = 0
      if (destination_atom == source_atom) then
         slot = 1
         return
      end if
      do ineigh = 2, lattice_obj%nn(destination_atom, 1)
         if (lattice_obj%nn(destination_atom, ineigh) == source_atom) then
            slot = ineigh
            return
         end if
      end do
   end function direct_hopping_slot

   integer function recursion_site_slot(lattice_obj, atom) result(slot)
      type(lattice), intent(in) :: lattice_obj
      integer, intent(in) :: atom
      integer :: i

      slot = 0
      if (.not. allocated(lattice_obj%irec)) return
      do i = 1, size(lattice_obj%irec)
         if (lattice_obj%irec(i) == atom) then
            slot = i
            return
         end if
      end do
   end function recursion_site_slot

   subroutine remove_radial_enu(block, basis)
      complex(rp), intent(inout) :: block(:, :)
      type(lmto_radial_basis), intent(in) :: basis
      integer :: l, orbital, ispin, first_orbital, last_orbital, index

      do ispin = 1, 2
         do l = 0, basis%lmax
            first_orbital = l*l + 1
            last_orbital = (l + 1)*(l + 1)
            do orbital = first_orbital, last_orbital
               index = orbital + merge(0, spin_off, ispin == 1)
               block(index, index) = block(index, index) - basis%enu_work(l + 1, ispin)
            end do
         end do
      end do
   end subroutine remove_radial_enu

   subroutine pair_translation(reciprocal_obj, lattice_obj, left_atom, right_atom, translation)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lattice), intent(in) :: lattice_obj
      integer, intent(in) :: left_atom, right_atom
      real(rp), intent(out) :: translation(3)
      integer :: ineigh, ntype

      translation = 0.0_rp
      if (left_atom == right_atom) return
      ntype = lattice_obj%iz(left_atom)
      do ineigh = 1, lattice_obj%nn(left_atom, 1)
         if (lattice_obj%nn(left_atom, ineigh) == right_atom) then
            if (ineigh > size(reciprocal_obj%ham_vec_type_direct, 2) .or. &
                ntype > size(reciprocal_obj%ham_vec_type_direct, 3)) exit
            translation = reciprocal_obj%ham_vec_type_direct(:, ineigh, ntype)
            return
         end if
      end do
      error stop 'native RSGF production provider: requested pair is not represented by a direct lattice bond'
   end subroutine pair_translation

   integer function response_site_for_atom(lattice_obj, atom, nsite) result(site)
      type(lattice), intent(in) :: lattice_obj
      integer, intent(in) :: atom, nsite
      integer :: isite

      site = 0
      do isite = 1, min(nsite, lattice_obj%nrec)
         if (lattice_obj%irec(isite) == atom) then
            site = isite
            return
         end if
      end do
      if (allocated(lattice_obj%iz)) then
         do isite = 1, min(nsite, lattice_obj%nrec)
            if (lattice_obj%iz(lattice_obj%irec(isite)) == lattice_obj%iz(atom)) then
               site = isite
               return
            end if
         end do
      end if
      error stop 'native RSGF production provider: pair atom has no accepted response-site representative'
   end function response_site_for_atom

   integer function find_pair(this, left_site, right_site, translation) result(ipair)
      class(tddft_native_rsgf_provider), intent(in) :: this
      integer, intent(in) :: left_site, right_site
      real(rp), intent(in) :: translation(3)
      integer :: candidate

      ipair = 0
      do candidate = 1, size(this%pairs)
         if (this%pairs(candidate)%left_site == left_site .and. this%pairs(candidate)%right_site == right_site .and. &
             maxval(abs(this%pairs(candidate)%translation - translation)) <= 2.0e-10_rp) then
            ipair = candidate
            return
         end if
      end do
      error stop 'native RSGF production provider: response pair is not in the registered complete pair set'
   end function find_pair

   subroutine validate_atom_index(lattice_obj, atom, side)
      type(lattice), intent(in) :: lattice_obj
      integer, intent(in) :: atom
      character(len=*), intent(in) :: side
      if (atom < 1 .or. atom > lattice_obj%kk) then
         error stop 'native RSGF production provider: '//trim(side)//' pair atom is outside the cluster'
      end if
   end subroutine validate_atom_index


! --- from lr_projected_reciprocal_chi0_mod ---

   module subroutine evaluate_projected_lehmann_chi0(request, result)
      type(projected_chi0_request), intent(in) :: request
      type(projected_chi0_result), intent(out) :: result

      type(projected_site_spin_contract), pointer :: contract
      type(lmto_product_response_basis), pointer :: product
      type(lr_electronic_state), pointer :: left_state, right_state
      type(pauli_endpoint_state) :: left_band, right_band
      complex(rp), allocatable :: transition(:)
      complex(rp), allocatable :: k_response(:, :, :)
      complex(rp) :: denominator, pair_factor
      real(rp) :: weight_sum, occupation_difference, transition_weight, score
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
      if (request%diagnostics) then
         allocate(result%k_susceptibility(contract%nsite, contract%nsite, size(request%frequencies), left_state%nk), &
            result%dominant_transitions(n_dominant_transitions), k_response(contract%nsite, contract%nsite, &
            size(request%frequencies)))
         result%k_susceptibility = cmplx(0.0_rp, 0.0_rp, rp)
         result%dominant_transitions%score = -1.0_rp
      end if
      result%susceptibility = cmplx(0.0_rp, 0.0_rp, rp)
      weight_sum = sum(left_state%k_weights)

      ! DRESP-02 production Lehmann path:
      ! eigenpairs -> certified DRESP-01 T_i -> site outer product.
      do ik = 1, left_state%nk
         if (request%diagnostics) k_response = cmplx(0.0_rp, 0.0_rp, rp)
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
               if (request%diagnostics) then
                  transition_weight = sum(abs(transition)**2)
                  score = abs(occupation_difference)*transition_weight
                  call retain_dominant_transition(result, ik, ib, jb, left_band%energy, right_band%energy, &
                     left_state%occupations(ib, ik), right_state%occupations(jb, ik), transition_weight, score)
               end if
               do ifrequency = 1, size(request%frequencies)
                  denominator = cmplx(request%frequencies(ifrequency) + left_band%energy - right_band%energy, &
                     request%eta, rp)
                  pair_factor = cmplx(2.0_rp*left_state%k_weights(ik)/weight_sum, 0.0_rp, rp)* &
                     occupation_difference/denominator
                  do site = 1, contract%nsite
                     if (request%diagnostics) then
                        k_response(:, site, ifrequency) = k_response(:, site, ifrequency) + &
                           pair_factor*transition*conjg(transition(site))
                     else
                        result%susceptibility(:, site, ifrequency) = result%susceptibility(:, site, ifrequency) + &
                           pair_factor*transition*conjg(transition(site))
                     end if
                  end do
               end do
            end do
         end do
         if (request%diagnostics) then
            result%susceptibility = result%susceptibility + k_response
            result%k_susceptibility(:, :, :, ik) = k_response
         end if
      end do

      call system_clock(wall_stop)
      call cpu_time(cpu_stop)
      call finalize_result(result, request, contract, product, channel_kind, projected_backend_lehmann, 0.0_rp, &
         0.0_rp, 0.0_rp, wall_start, wall_stop, clock_rate, cpu_start, cpu_stop, 0_int64)
      if (request%diagnostics) then
         call sort_dominant_transitions(result)
      end if
      deallocate(transition)
      if (allocated(k_response)) deallocate(k_response)
   end subroutine evaluate_projected_lehmann_chi0

   !> Independent finite-width spectral-function oracle for DRESP-02C.
   !>
   !> This path deliberately does not use the production GF eigenbasis
   !> contraction.  It rebuilds the certified DRESP-01 transition amplitude
   !> for every endpoint band pair and evaluates the complete two-term Kubo
   !> spectral integral on the requested real-energy mesh:
   !>   f(E) [ A_L(E) G_R^R(E+w) + G_L^A(E-w) A_R(E) ].
   !> Consequently it is the finite-regulator oracle for Track A, including
   !> finite-temperature Fermi weighting and the same finite energy window.
   module subroutine evaluate_projected_finite_width_chi0(request, result)
      type(projected_chi0_request), intent(in) :: request
      type(projected_chi0_result), intent(out) :: result

      type(projected_site_spin_contract), pointer :: contract
      type(lmto_product_response_basis), pointer :: product
      type(lr_electronic_state), pointer :: left_state, right_state
      type(pauli_endpoint_state) :: left_band, right_band
      complex(rp), allocatable :: transitions(:, :, :), transition(:)
      complex(rp), allocatable :: left_spectral(:), right_spectral(:)
      complex(rp), allocatable :: right_retarded(:), left_advanced(:)
      complex(rp) :: pair_factor
      real(rp) :: integration_eta, energy_min, energy_max, step, energy, fermi_weight
      real(rp) :: quadrature_weight, weight_sum, scale
      real(rp) :: cpu_start, cpu_stop
      integer(int64) :: wall_start, wall_stop, clock_rate
      integer :: channel_kind, ne, nfrequency, nleft, nright, ik, ie, ifrequency
      integer :: ib, jb, site, site_other

      call validate_request(request, channel_kind)
      if (request%integration_points < 3 .or. mod(request%integration_points, 2) == 0) then
         error stop 'evaluate_projected_finite_width_chi0: integration_points must be odd and at least three'
      end if
      if (request%energy_margin <= 0.0_rp) then
         error stop 'evaluate_projected_finite_width_chi0: energy_margin must be positive'
      end if

      contract => request%contract
      product => request%product_basis
      left_state => request%electronic_state
      right_state => request%q_endpoint_state
      call left_state%validate('evaluate_projected_finite_width_chi0:left_state')
      call right_state%validate('evaluate_projected_finite_width_chi0:q_endpoint_state')

      integration_eta = request%integration_eta
      if (integration_eta <= 0.0_rp) integration_eta = request%eta/40.0_rp
      if (integration_eta >= request%eta) then
         error stop 'evaluate_projected_finite_width_chi0: integration_eta must be smaller than response eta'
      end if

      ne = request%integration_points
      nfrequency = size(request%frequencies)
      nleft = left_state%nbands
      nright = right_state%nbands
      energy_min = min(minval(left_state%eigenvalues), minval(right_state%eigenvalues)) - request%energy_margin
      energy_max = max(maxval(left_state%eigenvalues), maxval(right_state%eigenvalues)) + request%energy_margin
      if (energy_max <= energy_min) error stop 'evaluate_projected_finite_width_chi0: invalid integration interval'
      step = (energy_max - energy_min)/real(ne - 1, rp)
      weight_sum = sum(left_state%k_weights)

      call system_clock(wall_start, clock_rate)
      call cpu_time(cpu_start)
      allocate(result%susceptibility(contract%nsite, contract%nsite, nfrequency), result%frequencies(nfrequency), &
         transitions(nleft, nright, contract%nsite), transition(contract%nsite), left_spectral(nleft), &
         right_spectral(nright), right_retarded(nright), left_advanced(nleft))
      result%susceptibility = cmplx(0.0_rp, 0.0_rp, rp)
      result%integration_points = ne
      result%ntransitions_evaluated = left_state%nk*nleft*nright

      do ik = 1, left_state%nk
         ! Direct DRESP-01 transition amplitudes are the only material
         ! matrix elements used by this oracle.  No production GF transform
         ! or resolvent contraction is shared here.
         do ib = 1, nleft
            call left_band%initialize(left_state%eigenvalues(ib, ik), left_state%eigenvectors(:, ib, ik))
            do jb = 1, nright
               call right_band%initialize(right_state%eigenvalues(jb, ik), right_state%eigenvectors(:, jb, ik))
               call contract%transition_amplitudes(product, left_band, right_band, transition)
               transitions(ib, jb, :) = transition
            end do
         end do

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
            call build_spectral_factors(left_state%eigenvalues(:, ik), energy, integration_eta, left_spectral)
            call build_spectral_factors(right_state%eigenvalues(:, ik), energy, integration_eta, right_spectral)

            do ifrequency = 1, nfrequency
               right_retarded = 1.0_rp/(cmplx(energy + request%frequencies(ifrequency), request%eta, rp) - &
                  right_state%eigenvalues(:, ik))
               left_advanced = 1.0_rp/(cmplx(energy - request%frequencies(ifrequency), -request%eta, rp) - &
                  left_state%eigenvalues(:, ik))
               scale = fermi_weight*quadrature_weight*2.0_rp*left_state%k_weights(ik)/weight_sum
               do ib = 1, nleft
                  do jb = 1, nright
                     pair_factor = scale*(left_spectral(ib)*right_retarded(jb) + &
                        right_spectral(jb)*left_advanced(ib))
                     do site = 1, contract%nsite
                        do site_other = 1, contract%nsite
                           result%susceptibility(site, site_other, ifrequency) = &
                              result%susceptibility(site, site_other, ifrequency) + pair_factor* &
                              transitions(ib, jb, site)*conjg(transitions(ib, jb, site_other))
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      call system_clock(wall_stop)
      call cpu_time(cpu_stop)
      call finalize_result(result, request, contract, product, channel_kind, projected_backend_finite_width, integration_eta, &
         energy_min, energy_max, wall_start, wall_stop, clock_rate, cpu_start, cpu_stop, 0_int64)
      result%energy_spacing = step
      result%implementation = 'direct-oracle'
      deallocate(transitions, transition, left_spectral, right_spectral, right_retarded, left_advanced)
   end subroutine evaluate_projected_finite_width_chi0

   !> Production projected GF backend.  The public seam deliberately keeps the
   !> old name used by DRESP-02 while routing production work to the optimized
   !> real-axis implementation.  The dense-resolvent implementation remains
   !> available as evaluate_projected_gf_chi0_reference for certification.
   module subroutine evaluate_projected_gf_chi0(request, result)
      type(projected_chi0_request), intent(in) :: request
      type(projected_chi0_result), intent(out) :: result

      call evaluate_projected_gf_chi0_optimized(request, result)
   end subroutine evaluate_projected_gf_chi0

   !> Evaluate the same explicit real-energy Kubo integral in endpoint
   !> eigenbases.  The energy integral is retained: only the already
   !> diagonalized resolvents and affine endpoint vertices are combined before
   !> the energy loop.  This is algebraically equivalent to the reference
   !> dense-resolvent route, but produces the site matrix directly.
   module subroutine evaluate_projected_gf_chi0_optimized(request, result)
      type(projected_chi0_request), intent(in) :: request
      type(projected_chi0_result), intent(out) :: result

      type(projected_site_spin_contract), pointer :: contract
      type(lmto_product_response_basis), pointer :: product
      type(lr_electronic_state), pointer :: left_state, right_state
      complex(rp), allocatable :: vertices(:, :, :, :), transitions(:, :, :)
      complex(rp), allocatable :: transition_flat(:, :), weighted_transition(:, :), contribution(:, :)
      complex(rp), allocatable :: left_spectral(:), right_spectral(:), right_retarded(:), left_advanced(:)
      complex(rp), allocatable :: kernel(:, :)
      complex(rp), allocatable :: k_response(:, :, :), k_term_one(:, :, :), k_term_two(:, :, :)
      complex(rp), allocatable :: left_zero(:, :), right_zero(:, :), left_first(:, :), right_first(:, :), &
         left_fermi(:, :), right_fermi(:, :), residual_matrix(:, :), scaled_vectors(:, :)
      real(rp) :: integration_eta, energy_min, energy_max, step, energy, fermi_weight
      real(rp) :: quadrature_weight, weight_sum, scale
      real(rp) :: cpu_start, cpu_stop
      integer(int64) :: wall_start, wall_stop, clock_rate, stage_start, stage_stop
      integer :: channel_kind, ne, nfrequency, nleft, nright, npairs, ik, ie, ifrequency

      call validate_request(request, channel_kind)
      if (request%integration_points < 3 .or. mod(request%integration_points, 2) == 0) then
         error stop 'evaluate_projected_gf_chi0_optimized: integration_points must be odd and at least three'
      end if
      if (request%energy_margin <= 0.0_rp) then
         error stop 'evaluate_projected_gf_chi0_optimized: energy_margin must be positive'
      end if
      contract => request%contract
      product => request%product_basis
      left_state => request%electronic_state
      right_state => request%q_endpoint_state
      call left_state%validate('evaluate_projected_gf_chi0_optimized:left_state')
      call right_state%validate('evaluate_projected_gf_chi0_optimized:q_endpoint_state')

      integration_eta = request%integration_eta
      if (integration_eta <= 0.0_rp) integration_eta = request%eta/40.0_rp
      if (integration_eta >= request%eta) then
         error stop 'evaluate_projected_gf_chi0_optimized: integration_eta must be smaller than response eta'
      end if

      call system_clock(wall_start, clock_rate)
      call cpu_time(cpu_start)
      call system_clock(stage_start)
      call contract%site_component_vertex_tensor(product, vertices)
      call system_clock(stage_stop)
      result%vertex_seconds = real(max(0_int64, stage_stop - stage_start), rp)/real(max(1_int64, clock_rate), rp)

      ne = request%integration_points
      nfrequency = size(request%frequencies)
      nleft = left_state%nbands
      nright = right_state%nbands
      npairs = nleft*nright
      energy_min = min(minval(left_state%eigenvalues), minval(right_state%eigenvalues)) - request%energy_margin
      energy_max = max(maxval(left_state%eigenvalues), maxval(right_state%eigenvalues)) + request%energy_margin
      if (energy_max <= energy_min) error stop 'evaluate_projected_gf_chi0_optimized: invalid integration interval'
      step = (energy_max - energy_min)/real(ne - 1, rp)
      weight_sum = sum(left_state%k_weights)

      call system_clock(stage_start)
      allocate(result%susceptibility(contract%nsite, contract%nsite, nfrequency), result%frequencies(nfrequency), &
         transitions(nleft, nright, contract%nsite), transition_flat(npairs, contract%nsite), &
         weighted_transition(npairs, contract%nsite), contribution(contract%nsite, contract%nsite), &
         left_spectral(nleft), right_spectral(nright), right_retarded(nright), left_advanced(nleft), kernel(nleft, nright))
      if (request%diagnostics) then
         allocate(result%kubo_term_one(contract%nsite, contract%nsite, nfrequency), &
            result%kubo_term_two(contract%nsite, contract%nsite, nfrequency), &
            result%k_susceptibility(contract%nsite, contract%nsite, nfrequency, left_state%nk), &
            k_response(contract%nsite, contract%nsite, nfrequency), &
            k_term_one(contract%nsite, contract%nsite, nfrequency), k_term_two(contract%nsite, contract%nsite, nfrequency), &
            left_zero(nleft, left_state%nk), right_zero(nright, left_state%nk), &
            left_first(nleft, left_state%nk), right_first(nright, left_state%nk), &
            left_fermi(nleft, left_state%nk), right_fermi(nright, left_state%nk), &
            residual_matrix(left_state%nbasis, left_state%nbasis), scaled_vectors(left_state%nbasis, nleft))
         result%kubo_term_one = cmplx(0.0_rp, 0.0_rp, rp)
         result%kubo_term_two = cmplx(0.0_rp, 0.0_rp, rp)
         result%k_susceptibility = cmplx(0.0_rp, 0.0_rp, rp)
         left_zero = cmplx(0.0_rp, 0.0_rp, rp)
         right_zero = cmplx(0.0_rp, 0.0_rp, rp)
         left_first = cmplx(0.0_rp, 0.0_rp, rp)
         right_first = cmplx(0.0_rp, 0.0_rp, rp)
         left_fermi = cmplx(0.0_rp, 0.0_rp, rp)
         right_fermi = cmplx(0.0_rp, 0.0_rp, rp)
      end if
      call system_clock(stage_stop)
      result%allocation_seconds = real(max(0_int64, stage_stop - stage_start), rp)/real(max(1_int64, clock_rate), rp)

      result%susceptibility = cmplx(0.0_rp, 0.0_rp, rp)
      result%integration_points = ne
      result%endpoint_transform_calls = int(left_state%nk*contract%nsite*lmto_product_nbranch, int64)

      ! P2/P3: combine the second-order endpoint vertex with the eigenvectors once
      ! per k.  All later work is a site-space outer product over the same
      ! explicit GF energy mesh.
      do ik = 1, left_state%nk
         call system_clock(stage_start)
         call build_projected_eigenbasis_transitions(vertices, left_state, right_state, ik, transitions)
         call system_clock(stage_stop)
         result%endpoint_transform_seconds = result%endpoint_transform_seconds + &
            real(max(0_int64, stage_stop - stage_start), rp)/real(max(1_int64, clock_rate), rp)
         transition_flat = reshape(transitions, [npairs, contract%nsite])

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

            call system_clock(stage_start)
            call build_spectral_factors(left_state%eigenvalues(:, ik), energy, integration_eta, left_spectral)
            call build_spectral_factors(right_state%eigenvalues(:, ik), energy, integration_eta, right_spectral)
            call system_clock(stage_stop)
            result%resolvent_seconds = result%resolvent_seconds + &
               real(max(0_int64, stage_stop - stage_start), rp)/real(max(1_int64, clock_rate), rp)

            if (request%diagnostics) then
               left_zero(:, ik) = left_zero(:, ik) + quadrature_weight*left_spectral
               right_zero(:, ik) = right_zero(:, ik) + quadrature_weight*right_spectral
               left_first(:, ik) = left_first(:, ik) + quadrature_weight*energy*left_spectral
               right_first(:, ik) = right_first(:, ik) + quadrature_weight*energy*right_spectral
               left_fermi(:, ik) = left_fermi(:, ik) + quadrature_weight*fermi_weight*left_spectral
               right_fermi(:, ik) = right_fermi(:, ik) + quadrature_weight*fermi_weight*right_spectral
               k_response = cmplx(0.0_rp, 0.0_rp, rp)
               k_term_one = cmplx(0.0_rp, 0.0_rp, rp)
               k_term_two = cmplx(0.0_rp, 0.0_rp, rp)
            end if

            do ifrequency = 1, nfrequency
               call system_clock(stage_start)
               right_retarded = 1.0_rp/(cmplx(energy + request%frequencies(ifrequency), request%eta, rp) - &
                  right_state%eigenvalues(:, ik))
               left_advanced = 1.0_rp/(cmplx(energy - request%frequencies(ifrequency), -request%eta, rp) - &
                  left_state%eigenvalues(:, ik))
               call system_clock(stage_stop)
               result%resolvent_seconds = result%resolvent_seconds + &
                  real(max(0_int64, stage_stop - stage_start), rp)/real(max(1_int64, clock_rate), rp)

               call system_clock(stage_start)
               kernel = spread(left_spectral, 2, nright)*spread(right_retarded, 1, nleft)
               weighted_transition = transition_flat*spread(reshape(kernel, [npairs]), 2, contract%nsite)
               contribution = matmul(transpose(weighted_transition), conjg(transition_flat))
               scale = fermi_weight*quadrature_weight*2.0_rp*left_state%k_weights(ik)/weight_sum
               if (request%diagnostics) then
                  k_term_one(:, :, ifrequency) = k_term_one(:, :, ifrequency) + scale*contribution
                  k_response(:, :, ifrequency) = k_response(:, :, ifrequency) + scale*contribution
               else
                  result%susceptibility(:, :, ifrequency) = result%susceptibility(:, :, ifrequency) + scale*contribution
               end if

               kernel = spread(left_advanced, 2, nright)*spread(right_spectral, 1, nleft)
               weighted_transition = transition_flat*spread(reshape(kernel, [npairs]), 2, contract%nsite)
               contribution = matmul(transpose(weighted_transition), conjg(transition_flat))
               if (request%diagnostics) then
                  k_term_two(:, :, ifrequency) = k_term_two(:, :, ifrequency) + scale*contribution
                  k_response(:, :, ifrequency) = k_response(:, :, ifrequency) + scale*contribution
               else
                  result%susceptibility(:, :, ifrequency) = result%susceptibility(:, :, ifrequency) + scale*contribution
               end if
               call system_clock(stage_stop)
               result%accumulator_seconds = result%accumulator_seconds + &
                  real(max(0_int64, stage_stop - stage_start), rp)/real(max(1_int64, clock_rate), rp)
               result%accumulator_calls = result%accumulator_calls + 1_int64
            end do
            if (request%diagnostics) then
               result%susceptibility = result%susceptibility + k_response
               result%kubo_term_one = result%kubo_term_one + k_term_one
               result%kubo_term_two = result%kubo_term_two + k_term_two
               result%k_susceptibility(:, :, :, ik) = result%k_susceptibility(:, :, :, ik) + k_response
            end if
         end do
      end do

      if (request%diagnostics) then
         call system_clock(stage_start)
         do ik = 1, left_state%nk
            scaled_vectors = left_state%eigenvectors(:, :, ik)
            do ie = 1, nleft
               scaled_vectors(:, ie) = left_state%eigenvectors(:, ie, ik)*(left_zero(ie, ik) - 1.0_rp)
            end do
            residual_matrix = matmul(scaled_vectors, conjg(transpose(left_state%eigenvectors(:, :, ik))))
            result%left_spectral_zeroth_residual = max(result%left_spectral_zeroth_residual, maxval(abs(residual_matrix)))
            do ie = 1, nleft
               scaled_vectors(:, ie) = left_state%eigenvectors(:, ie, ik)*(left_first(ie, ik) - &
                  left_state%eigenvalues(ie, ik))
            end do
            residual_matrix = matmul(scaled_vectors, conjg(transpose(left_state%eigenvectors(:, :, ik))))
            result%left_spectral_first_residual = max(result%left_spectral_first_residual, maxval(abs(residual_matrix)))
            do ie = 1, nleft
               scaled_vectors(:, ie) = left_state%eigenvectors(:, ie, ik)*(left_fermi(ie, ik) - &
                  left_state%occupations(ie, ik))
            end do
            residual_matrix = matmul(scaled_vectors, conjg(transpose(left_state%eigenvectors(:, :, ik))))
            result%left_spectral_fermi_residual = max(result%left_spectral_fermi_residual, maxval(abs(residual_matrix)))

            scaled_vectors = right_state%eigenvectors(:, :, ik)
            do ie = 1, nright
               scaled_vectors(:, ie) = right_state%eigenvectors(:, ie, ik)*(right_zero(ie, ik) - 1.0_rp)
            end do
            residual_matrix = matmul(scaled_vectors, conjg(transpose(right_state%eigenvectors(:, :, ik))))
            result%right_spectral_zeroth_residual = max(result%right_spectral_zeroth_residual, maxval(abs(residual_matrix)))
            do ie = 1, nright
               scaled_vectors(:, ie) = right_state%eigenvectors(:, ie, ik)*(right_first(ie, ik) - &
                  right_state%eigenvalues(ie, ik))
            end do
            residual_matrix = matmul(scaled_vectors, conjg(transpose(right_state%eigenvectors(:, :, ik))))
            result%right_spectral_first_residual = max(result%right_spectral_first_residual, maxval(abs(residual_matrix)))
            do ie = 1, nright
               scaled_vectors(:, ie) = right_state%eigenvectors(:, ie, ik)*(right_fermi(ie, ik) - &
                  right_state%occupations(ie, ik))
            end do
            residual_matrix = matmul(scaled_vectors, conjg(transpose(right_state%eigenvectors(:, :, ik))))
            result%right_spectral_fermi_residual = max(result%right_spectral_fermi_residual, maxval(abs(residual_matrix)))
         end do
         call system_clock(stage_stop)
         result%diagnostic_seconds = real(max(0_int64, stage_stop - stage_start), rp)/real(max(1_int64, clock_rate), rp)
      end if

      call system_clock(wall_stop)
      call cpu_time(cpu_stop)
      call finalize_result(result, request, contract, product, channel_kind, projected_backend_gf, integration_eta, &
         energy_min, energy_max, wall_start, wall_stop, clock_rate, cpu_start, cpu_stop, &
         int(size(vertices), int64)*int(storage_size(vertices)/8, int64))
      result%energy_spacing = step
      result%implementation = 'eigenbasis'
      deallocate(vertices, transitions, transition_flat, weighted_transition, contribution, left_spectral, right_spectral, &
         right_retarded, left_advanced, kernel)
      if (request%diagnostics) then
         deallocate(k_response, k_term_one, k_term_two, left_zero, right_zero, left_first, right_first, left_fermi, &
            right_fermi, residual_matrix, scaled_vectors)
      end if
   end subroutine evaluate_projected_gf_chi0_optimized

   !> Transform the four affine site-vertex components into the endpoint
   !> eigenbases and sum their exact endpoint-energy factors.  For a pair n,m
   !> this produces
   !>   T_i(n,m) = sum_pq eps_L(n)^p eps_R(m)^q <n|V_i^(p,q)|m>.
   !> It is a precomputation for the GF integral, not a Lehmann response
   !> substitution: the energy-dependent GF denominators are still integrated
   !> explicitly by evaluate_projected_gf_chi0_optimized.
   subroutine build_projected_eigenbasis_transitions(vertices, left_state, right_state, ik, transitions)
      complex(rp), intent(in) :: vertices(:, :, :, :)
      type(lr_electronic_state), intent(in) :: left_state, right_state
      integer, intent(in) :: ik
      complex(rp), intent(out) :: transitions(:, :, :)

      complex(rp), allocatable :: projected_vertex(:, :), band_vertex(:, :)
      integer :: component, site, ib_left, ib_right, left_power, right_power
      integer :: nbasis, nleft, nright, nsite

      nbasis = size(vertices, 1)
      nleft = left_state%nbands
      nright = right_state%nbands
      nsite = size(vertices, 4)
      if (size(vertices, 2) /= nbasis .or. size(vertices, 3) /= lmto_product_nbranch .or. size(transitions, 1) /= nleft .or. &
          size(transitions, 2) /= nright .or. size(transitions, 3) /= nsite .or. left_state%nbasis /= nbasis .or. &
          right_state%nbasis /= nbasis) then
         error stop 'build_projected_eigenbasis_transitions: shape mismatch'
      end if
      if (ik < 1 .or. ik > left_state%nk .or. ik > right_state%nk) then
         error stop 'build_projected_eigenbasis_transitions: k-point index out of range'
      end if

      allocate(projected_vertex(nbasis, nright), band_vertex(nleft, nright))
      transitions = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, nsite
         do component = 1, lmto_product_nbranch
            call lmto_product_branch_powers(component, left_power, right_power)
            projected_vertex = matmul(vertices(:, :, component, site), right_state%eigenvectors(:, :, ik))
            band_vertex = matmul(conjg(transpose(left_state%eigenvectors(:, :, ik))), projected_vertex)
            do ib_left = 1, nleft
               band_vertex(ib_left, :) = band_vertex(ib_left, :)* &
                  lmto_product_energy_power(left_state%eigenvalues(ib_left, ik), left_power)
            end do
            do ib_right = 1, nright
               band_vertex(:, ib_right) = band_vertex(:, ib_right)* &
                  lmto_product_energy_power(right_state%eigenvalues(ib_right, ik), right_power)
            end do
            transitions(:, :, site) = transitions(:, :, site) + band_vertex
         end do
      end do
      deallocate(projected_vertex, band_vertex)
   end subroutine build_projected_eigenbasis_transitions

   subroutine build_spectral_factors(eigenvalues, energy, integration_eta, factors)
      real(rp), intent(in) :: eigenvalues(:), energy, integration_eta
      complex(rp), intent(out) :: factors(:)
      complex(rp) :: spectral_prefactor
      integer :: ib

      spectral_prefactor = cmplx(0.0_rp, 1.0_rp/(2.0_rp*response_angular_pi), rp)
      do ib = 1, size(eigenvalues)
         factors(ib) = spectral_prefactor*(1.0_rp/(cmplx(energy, integration_eta, rp) - eigenvalues(ib)) - &
            1.0_rp/(cmplx(energy, -integration_eta, rp) - eigenvalues(ib)))
      end do
   end subroutine build_spectral_factors

   module subroutine evaluate_projected_gf_chi0_reference(request, result)
      type(projected_chi0_request), intent(in) :: request
      type(projected_chi0_result), intent(out) :: result

      type(projected_site_spin_contract), pointer :: contract
      type(lmto_product_response_basis), pointer :: product
      type(lr_electronic_state), pointer :: left_state, right_state
      complex(rp), allocatable :: vertices(:, :, :, :)
      complex(rp), allocatable :: left_gr(:, :, :), left_ga(:, :, :), left_a(:, :, :)
      complex(rp), allocatable :: right_gr(:, :, :), right_ga(:, :, :), right_a(:, :, :)
      complex(rp), allocatable :: k_response(:, :, :), k_term_one(:, :, :), k_term_two(:, :, :)
      complex(rp), allocatable :: left_zero(:, :, :), right_zero(:, :, :), left_first(:, :, :), right_first(:, :, :), &
         left_fermi(:, :, :), right_fermi(:, :, :)
      complex(rp), allocatable :: exact_zero(:, :), exact_first(:, :), exact_fermi(:, :), &
         exact_right_zero(:, :), exact_right_first(:, :), exact_right_fermi(:, :)
      real(rp) :: integration_eta, energy_min, energy_max, step, energy, fermi_weight
      real(rp) :: quadrature_weight, weight_sum
      real(rp) :: cpu_start, cpu_stop
      integer(int64) :: wall_start, wall_stop, clock_rate, stage_start, stage_stop
      integer :: channel_kind, nbasis, ne, ie, ik, ifrequency

      call validate_request(request, channel_kind)
      if (request%integration_points < 3 .or. mod(request%integration_points, 2) == 0) then
         error stop 'evaluate_projected_gf_chi0_reference: integration_points must be odd and at least three'
      end if
      if (request%energy_margin <= 0.0_rp) then
         error stop 'evaluate_projected_gf_chi0_reference: energy_margin must be positive'
      end if
      contract => request%contract
      product => request%product_basis
      left_state => request%electronic_state
      right_state => request%q_endpoint_state
      call left_state%validate('evaluate_projected_gf_chi0_reference:left_state')
      call right_state%validate('evaluate_projected_gf_chi0_reference:q_endpoint_state')

      integration_eta = request%integration_eta
      if (integration_eta <= 0.0_rp) integration_eta = request%eta/40.0_rp
      if (integration_eta >= request%eta) then
         error stop 'evaluate_projected_gf_chi0_reference: integration_eta must be smaller than response eta'
      end if

      call system_clock(wall_start, clock_rate)
      call cpu_time(cpu_start)
      call system_clock(stage_start)
      call contract%site_component_vertex_tensor(product, vertices)
      call system_clock(stage_stop)
      result%vertex_seconds = real(max(0_int64, stage_stop - stage_start), rp)/real(max(1_int64, clock_rate), rp)
      nbasis = left_state%nbasis
      ne = request%integration_points
      energy_min = min(minval(left_state%eigenvalues), minval(right_state%eigenvalues)) - request%energy_margin
      energy_max = max(maxval(left_state%eigenvalues), maxval(right_state%eigenvalues)) + request%energy_margin
      if (energy_max <= energy_min) error stop 'evaluate_projected_gf_chi0_reference: invalid integration interval'
      step = (energy_max - energy_min)/real(ne - 1, rp)
      weight_sum = sum(left_state%k_weights)

      call system_clock(stage_start)
      allocate(result%susceptibility(contract%nsite, contract%nsite, size(request%frequencies)), &
         result%frequencies(size(request%frequencies)), &
         left_gr(nbasis, nbasis, lmto_product_max_gf_moment + 1), &
         left_ga(nbasis, nbasis, lmto_product_max_gf_moment + 1), &
         left_a(nbasis, nbasis, lmto_product_max_gf_moment + 1), &
         right_gr(nbasis, nbasis, lmto_product_max_gf_moment + 1), &
         right_ga(nbasis, nbasis, lmto_product_max_gf_moment + 1), &
         right_a(nbasis, nbasis, lmto_product_max_gf_moment + 1))
      if (request%diagnostics) then
         allocate(result%kubo_term_one(contract%nsite, contract%nsite, size(request%frequencies)), &
            result%kubo_term_two(contract%nsite, contract%nsite, size(request%frequencies)), &
            result%k_susceptibility(contract%nsite, contract%nsite, size(request%frequencies), left_state%nk), &
            k_response(contract%nsite, contract%nsite, size(request%frequencies)), &
            k_term_one(contract%nsite, contract%nsite, size(request%frequencies)), &
            k_term_two(contract%nsite, contract%nsite, size(request%frequencies)), &
            left_zero(nbasis, nbasis, left_state%nk), right_zero(nbasis, nbasis, left_state%nk), &
            left_first(nbasis, nbasis, left_state%nk), right_first(nbasis, nbasis, left_state%nk), &
            left_fermi(nbasis, nbasis, left_state%nk), right_fermi(nbasis, nbasis, left_state%nk), &
            exact_zero(nbasis, nbasis), exact_first(nbasis, nbasis), exact_fermi(nbasis, nbasis), &
            exact_right_zero(nbasis, nbasis), exact_right_first(nbasis, nbasis), exact_right_fermi(nbasis, nbasis))
         result%kubo_term_one = cmplx(0.0_rp, 0.0_rp, rp)
         result%kubo_term_two = cmplx(0.0_rp, 0.0_rp, rp)
         result%k_susceptibility = cmplx(0.0_rp, 0.0_rp, rp)
         left_zero = cmplx(0.0_rp, 0.0_rp, rp)
         right_zero = cmplx(0.0_rp, 0.0_rp, rp)
         left_first = cmplx(0.0_rp, 0.0_rp, rp)
         right_first = cmplx(0.0_rp, 0.0_rp, rp)
         left_fermi = cmplx(0.0_rp, 0.0_rp, rp)
         right_fermi = cmplx(0.0_rp, 0.0_rp, rp)
      end if
      call system_clock(stage_stop)
      result%allocation_seconds = real(max(0_int64, stage_stop - stage_start), rp)/real(max(1_int64, clock_rate), rp)
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
            call system_clock(stage_start)
            call build_weighted_resolvent(left_state, ik, cmplx(energy, integration_eta, rp), left_gr)
            call build_weighted_resolvent(left_state, ik, cmplx(energy, -integration_eta, rp), left_ga)
            call build_weighted_resolvent(right_state, ik, cmplx(energy, integration_eta, rp), right_gr)
            call build_weighted_resolvent(right_state, ik, cmplx(energy, -integration_eta, rp), right_ga)
            call system_clock(stage_stop)
            result%resolvent_seconds = result%resolvent_seconds + &
               real(max(0_int64, stage_stop - stage_start), rp)/real(max(1_int64, clock_rate), rp)
            result%resolvent_calls = result%resolvent_calls + 4_int64
            left_a = cmplx(0.0_rp, 1.0_rp/(2.0_rp*response_angular_pi), rp)*(left_gr - left_ga)
            right_a = cmplx(0.0_rp, 1.0_rp/(2.0_rp*response_angular_pi), rp)*(right_gr - right_ga)
            if (request%diagnostics) then
               left_zero(:, :, ik) = left_zero(:, :, ik) + quadrature_weight*left_a(:, :, 1)
               right_zero(:, :, ik) = right_zero(:, :, ik) + quadrature_weight*right_a(:, :, 1)
               left_first(:, :, ik) = left_first(:, :, ik) + quadrature_weight*energy*left_a(:, :, 1)
               right_first(:, :, ik) = right_first(:, :, ik) + quadrature_weight*energy*right_a(:, :, 1)
               left_fermi(:, :, ik) = left_fermi(:, :, ik) + quadrature_weight*fermi_weight*left_a(:, :, 1)
               right_fermi(:, :, ik) = right_fermi(:, :, ik) + quadrature_weight*fermi_weight*right_a(:, :, 1)
               k_response = cmplx(0.0_rp, 0.0_rp, rp)
               k_term_one = cmplx(0.0_rp, 0.0_rp, rp)
               k_term_two = cmplx(0.0_rp, 0.0_rp, rp)
            end if
            do ifrequency = 1, size(request%frequencies)
               call system_clock(stage_start)
               call build_weighted_resolvent(right_state, ik, &
                  cmplx(energy + request%frequencies(ifrequency), request%eta, rp), right_gr)
               call build_weighted_resolvent(left_state, ik, &
                  cmplx(energy - request%frequencies(ifrequency), -request%eta, rp), left_ga)
               call system_clock(stage_stop)
               result%resolvent_seconds = result%resolvent_seconds + &
                  real(max(0_int64, stage_stop - stage_start), rp)/real(max(1_int64, clock_rate), rp)
               result%resolvent_calls = result%resolvent_calls + 2_int64
               if (request%diagnostics) then
                  call system_clock(stage_start)
                  call accumulate_projected_gf_bubble(vertices, left_a, right_a, right_gr, left_ga, &
                     fermi_weight*quadrature_weight*2.0_rp*left_state%k_weights(ik)/weight_sum, &
                     k_response(:, :, ifrequency), k_term_one(:, :, ifrequency), k_term_two(:, :, ifrequency))
                  call system_clock(stage_stop)
                  result%accumulator_seconds = result%accumulator_seconds + &
                     real(max(0_int64, stage_stop - stage_start), rp)/real(max(1_int64, clock_rate), rp)
                  result%accumulator_calls = result%accumulator_calls + 1_int64
               else
                  call system_clock(stage_start)
                  call accumulate_projected_gf_bubble(vertices, left_a, right_a, right_gr, left_ga, &
                     fermi_weight*quadrature_weight*2.0_rp*left_state%k_weights(ik)/weight_sum, &
                     result%susceptibility(:, :, ifrequency))
                  call system_clock(stage_stop)
                  result%accumulator_seconds = result%accumulator_seconds + &
                     real(max(0_int64, stage_stop - stage_start), rp)/real(max(1_int64, clock_rate), rp)
                  result%accumulator_calls = result%accumulator_calls + 1_int64
               end if
            end do
            if (request%diagnostics) then
               result%susceptibility = result%susceptibility + k_response
               result%kubo_term_one = result%kubo_term_one + k_term_one
               result%kubo_term_two = result%kubo_term_two + k_term_two
               result%k_susceptibility(:, :, :, ik) = result%k_susceptibility(:, :, :, ik) + k_response
            end if
         end do
      end do

      if (request%diagnostics) then
         call system_clock(stage_start)
         do ik = 1, left_state%nk
            exact_zero = cmplx(0.0_rp, 0.0_rp, rp)
            exact_first = cmplx(0.0_rp, 0.0_rp, rp)
            exact_fermi = cmplx(0.0_rp, 0.0_rp, rp)
            exact_right_zero = cmplx(0.0_rp, 0.0_rp, rp)
            exact_right_first = cmplx(0.0_rp, 0.0_rp, rp)
            exact_right_fermi = cmplx(0.0_rp, 0.0_rp, rp)
            do ifrequency = 1, left_state%nbands
               exact_zero = exact_zero + matmul(reshape(left_state%eigenvectors(:, ifrequency, ik), [nbasis, 1]), &
                  reshape(conjg(left_state%eigenvectors(:, ifrequency, ik)), [1, nbasis]))
               exact_first = exact_first + left_state%eigenvalues(ifrequency, ik)* &
                  matmul(reshape(left_state%eigenvectors(:, ifrequency, ik), [nbasis, 1]), &
                  reshape(conjg(left_state%eigenvectors(:, ifrequency, ik)), [1, nbasis]))
               exact_fermi = exact_fermi + left_state%occupations(ifrequency, ik)* &
                  matmul(reshape(left_state%eigenvectors(:, ifrequency, ik), [nbasis, 1]), &
                  reshape(conjg(left_state%eigenvectors(:, ifrequency, ik)), [1, nbasis]))
               exact_right_zero = exact_right_zero + matmul(reshape(right_state%eigenvectors(:, ifrequency, ik), [nbasis, 1]), &
                  reshape(conjg(right_state%eigenvectors(:, ifrequency, ik)), [1, nbasis]))
               exact_right_first = exact_right_first + right_state%eigenvalues(ifrequency, ik)* &
                  matmul(reshape(right_state%eigenvectors(:, ifrequency, ik), [nbasis, 1]), &
                  reshape(conjg(right_state%eigenvectors(:, ifrequency, ik)), [1, nbasis]))
               exact_right_fermi = exact_right_fermi + right_state%occupations(ifrequency, ik)* &
                  matmul(reshape(right_state%eigenvectors(:, ifrequency, ik), [nbasis, 1]), &
                  reshape(conjg(right_state%eigenvectors(:, ifrequency, ik)), [1, nbasis]))
            end do
            result%left_spectral_zeroth_residual = max(result%left_spectral_zeroth_residual, &
               maxval(abs(left_zero(:, :, ik) - exact_zero)))
            result%right_spectral_zeroth_residual = max(result%right_spectral_zeroth_residual, &
               maxval(abs(right_zero(:, :, ik) - exact_right_zero)))
            result%left_spectral_first_residual = max(result%left_spectral_first_residual, &
               maxval(abs(left_first(:, :, ik) - exact_first)))
            result%right_spectral_first_residual = max(result%right_spectral_first_residual, &
               maxval(abs(right_first(:, :, ik) - exact_right_first)))
            result%left_spectral_fermi_residual = max(result%left_spectral_fermi_residual, &
               maxval(abs(left_fermi(:, :, ik) - exact_fermi)))
            result%right_spectral_fermi_residual = max(result%right_spectral_fermi_residual, &
               maxval(abs(right_fermi(:, :, ik) - exact_right_fermi)))
         end do
         call system_clock(stage_stop)
         result%diagnostic_seconds = real(max(0_int64, stage_stop - stage_start), rp)/real(max(1_int64, clock_rate), rp)
      end if

      call system_clock(wall_stop)
      call cpu_time(cpu_stop)
      result%integration_points = ne
      call finalize_result(result, request, contract, product, channel_kind, projected_backend_gf, integration_eta, &
         energy_min, energy_max, wall_start, wall_stop, clock_rate, cpu_start, cpu_stop, &
         int(size(vertices), int64)*int(storage_size(vertices)/8, int64))
      result%energy_spacing = step
      result%implementation = 'dense-reference'
      if (allocated(k_response)) deallocate(k_response, k_term_one, k_term_two, left_zero, right_zero, left_first, &
         right_first, left_fermi, right_fermi, exact_zero, exact_first, exact_fermi, exact_right_zero, &
         exact_right_first, exact_right_fermi)
      deallocate(vertices, left_gr, left_ga, left_a, right_gr, right_ga, right_a)
   end subroutine evaluate_projected_gf_chi0_reference

   subroutine accumulate_projected_gf_bubble(vertices, left_a, right_a, right_gr, left_ga, scale, response, term_one, term_two)
      complex(rp), intent(in) :: vertices(:, :, :, :), left_a(:, :, :), right_a(:, :, :), &
         right_gr(:, :, :), left_ga(:, :, :)
      real(rp), intent(in) :: scale
      complex(rp), intent(inout) :: response(:, :)
      complex(rp), intent(inout), optional :: term_one(:, :), term_two(:, :)

      integer :: component_i, component_j, left_power_i, right_power_i, left_power_j, right_power_j
      integer :: combined_left, combined_right, site_i, site_j
      complex(rp) :: temporary(size(left_a, 1), size(left_a, 2))
      complex(rp) :: vertex_adjoint(size(left_a, 1), size(left_a, 2))

      do component_i = 1, lmto_product_nbranch
         call lmto_product_branch_powers(component_i, left_power_i, right_power_i)
         if (maxval(abs(vertices(:, :, component_i, :))) == 0.0_rp) cycle
         do component_j = 1, lmto_product_nbranch
            call lmto_product_branch_powers(component_j, left_power_j, right_power_j)
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
                  if (present(term_one)) term_one(site_i, site_j) = term_one(site_i, site_j) + scale* &
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
                  if (present(term_two)) term_two(site_i, site_j) = term_two(site_i, site_j) + scale* &
                     sum(temporary*transpose(vertices(:, :, component_i, site_i)))
               end do
            end do
         end do
      end do
   end subroutine accumulate_projected_gf_bubble

   subroutine retain_dominant_transition(result, k_index, left_band, right_band, left_energy, right_energy, &
                                         left_occupation, right_occupation, matrix_element_weight, score)
      type(projected_chi0_result), intent(inout) :: result
      integer, intent(in) :: k_index, left_band, right_band
      real(rp), intent(in) :: left_energy, right_energy, left_occupation, right_occupation, &
         matrix_element_weight, score
      integer :: slot

      if (.not. allocated(result%dominant_transitions)) return
      if (result%n_dominant_transition_records < size(result%dominant_transitions)) then
         result%n_dominant_transition_records = result%n_dominant_transition_records + 1
         slot = result%n_dominant_transition_records
      else
         slot = minloc([(result%dominant_transitions(slot)%score, slot=1, size(result%dominant_transitions))], 1)
         if (score <= result%dominant_transitions(slot)%score) return
      end if
      result%dominant_transitions(slot)%k_index = k_index
      result%dominant_transitions(slot)%left_band = left_band
      result%dominant_transitions(slot)%right_band = right_band
      result%dominant_transitions(slot)%left_energy = left_energy
      result%dominant_transitions(slot)%right_energy = right_energy
      result%dominant_transitions(slot)%left_occupation = left_occupation
      result%dominant_transitions(slot)%right_occupation = right_occupation
      result%dominant_transitions(slot)%transition_energy = right_energy - left_energy
      result%dominant_transitions(slot)%matrix_element_weight = matrix_element_weight
      result%dominant_transitions(slot)%score = score
   end subroutine retain_dominant_transition

   subroutine sort_dominant_transitions(result)
      type(projected_chi0_result), intent(inout) :: result
      type(projected_chi0_transition_record) :: temporary
      integer :: i, j, best, n

      if (.not. allocated(result%dominant_transitions)) return
      n = result%n_dominant_transition_records
      do i = 1, n - 1
         best = i
         do j = i + 1, n
            if (result%dominant_transitions(j)%score > result%dominant_transitions(best)%score) best = j
         end do
         if (best /= i) then
            temporary = result%dominant_transitions(i)
            result%dominant_transitions(i) = result%dominant_transitions(best)
            result%dominant_transitions(best) = temporary
         end if
      end do
   end subroutine sort_dominant_transitions

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
      if (maxval(abs(left_state%k_weights - right_state%k_weights)) > projected_state_match_tolerance) then
         error stop 'DRESP-02: endpoint k weights differ'
      end if
      scale = max(1.0_rp, abs(left_state%fermi_level), abs(right_state%fermi_level), &
         abs(left_state%temperature), abs(right_state%temperature), abs(left_state%energy_zero), &
         abs(right_state%energy_zero))
      if (abs(left_state%fermi_level - right_state%fermi_level) > projected_state_match_tolerance*scale .or. &
          abs(left_state%temperature - right_state%temperature) > projected_state_match_tolerance*scale .or. &
          abs(left_state%energy_zero - right_state%energy_zero) > projected_state_match_tolerance*scale) then
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
      if (result%integration_points > 1 .and. integration_eta > 0.0_rp) then
         result%spacing_over_integration_eta = (energy_max - energy_min)/ &
            real(result%integration_points - 1, rp)/integration_eta
      else
         result%spacing_over_integration_eta = 0.0_rp
      end if
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


end submodule linear_response_bare
