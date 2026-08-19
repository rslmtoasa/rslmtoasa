submodule (reciprocal_mod) reciprocal_bands
   use magnetic_representation_mod, only: gbt_single_q
   implicit none

contains

   !> @brief Diagonalize the active k-space Hamiltonians.
   !> @details Solves standard or generalized eigenproblems depending on
   !>          reciprocal_mode and stores eigenvalues/eigenvectors for bands and DOS.
   !> @param[inout] this Reciprocal object receiving eigenvalue/eigenvector arrays.
   module subroutine diagonalize_hamiltonian(this)
      class(reciprocal), intent(inout) :: this
      integer :: nk, ik, nmat
      logical :: use_generalized, operator_changed
      real(rp) :: max_herm, max_herm_all, matrix_scale
      character(len=256) :: herm_msg
      type(reciprocal_execution_request) :: request
      type(reciprocal_execution_result) :: result

#ifndef RSLMTO_DISABLE_FUSED_RECIPROCAL
      ! The fused normal-mesh route already owns the completed eigensystem.
      ! Keep this public adapter as a cache check on the platforms where that
      ! route is enabled; submitting the host H(k) back to LAPACK would defeat
      ! the GC-04 single-crossing contract (and is invalid for an empty mesh).
      call this%invalidate_if_operator_changed('reciprocal%diagonalize_hamiltonian', operator_changed)
      if (operator_changed .or. .not. allocated(this%hk_bulk) .or. .not. allocated(this%eigenvalues) .or. &
          .not. allocated(this%eigenvectors)) call this%build_kspace_hamiltonian()
      if (.not. allocated(this%hk_bulk) .or. .not. allocated(this%eigenvalues) .or. .not. allocated(this%eigenvectors)) then
         call g_logger%error('diagonalize_hamiltonian: normal-mesh combined execution did not materialize eigenpairs.', &
                             __FILE__, __LINE__)
         return
      end if
      if (associated(this%hamiltonian)) then
         if (this%cached_operator_generation /= this%hamiltonian%operator_generation) then
            call g_logger%error('diagonalize_hamiltonian: cache generation is stale after normal-mesh execution.', __FILE__, __LINE__)
            return
         end if
      end if
      if (this%kanpur_diagnostics) call this%print_kanpur_mapping()
      if (this%gamma_bounds_diagnostics) call this%run_gamma_bounds_diagnostics()
      if (this%hall_diag_experimental) call this%diagonalize_hall_experimental()
      call root_info('diagonalize_hamiltonian: Completed successfully', __FILE__, __LINE__)
      return
#else
      call this%invalidate_if_operator_changed('reciprocal%diagonalize_hamiltonian', operator_changed)
      if (operator_changed) call this%build_kspace_hamiltonian()

      if (.not. allocated(this%hk_bulk)) then
         call g_logger%error('diagonalize_hamiltonian: hk_bulk not built - call build_kspace_hamiltonian first', &
                            __FILE__, __LINE__)
         return
      end if

      nmat = size(this%hk_bulk, 1)
      nk = size(this%hk_bulk, 3)
      call root_info('diagonalize_hamiltonian: Diagonalizing ' // trim(int2str(nk)) // ' k-points', __FILE__, __LINE__)
      call root_info('diagonalize_hamiltonian: Matrix size = ' // trim(int2str(nmat)) // ' x ' // trim(int2str(nmat)), __FILE__, __LINE__)
      if (nk == 0) then
         call root_info('diagonalize_hamiltonian: Empty local mesh, no backend solve required', __FILE__, __LINE__)
         return
      end if

      use_generalized = trim(this%reciprocal_mode) == 'generalized_overlap_proxy'
      if (use_generalized .and. .not. allocated(this%sk_overlap)) call this%build_kspace_overlap()
      if (use_generalized .and. .not. allocated(this%sk_overlap)) then
         call g_logger%warning('diagonalize_hamiltonian: S(k) unavailable, falling back to ham_only.', __FILE__, __LINE__)
         use_generalized = .false.
      end if

      max_herm_all = 0.0_rp
      do ik = 1, nk
         max_herm = maxval(abs(this%hk_bulk(:, :, ik) - transpose(conjg(this%hk_bulk(:, :, ik)))))
         max_herm_all = max(max_herm_all, max_herm)
         matrix_scale = max(1.0_rp, maxval(abs(this%hk_bulk(:, :, ik))))
         if (max_herm > 1.0e-10_rp * matrix_scale) then
            write(herm_msg, '(A,I0,A,ES12.4,A,ES12.4)') 'H(k) is non-Hermitian before eigensolution at k=', ik, &
               ': max|H-H^H|=', max_herm, ', scale=', matrix_scale
            call g_logger%fatal('diagonalize_hamiltonian: ' // trim(herm_msg), __FILE__, __LINE__)
         end if
      end do
      if (allocated(this%sk_overlap)) then
         do ik = 1, size(this%sk_overlap, 3)
            max_herm = maxval(abs(this%sk_overlap(:, :, ik) - transpose(conjg(this%sk_overlap(:, :, ik)))))
            matrix_scale = max(1.0_rp, maxval(abs(this%sk_overlap(:, :, ik))))
            if (max_herm > 1.0e-10_rp * matrix_scale) then
               write(herm_msg, '(A,I0,A,ES12.4,A,ES12.4)') 'O(k) is non-Hermitian before eigensolution at k=', ik, &
                  ': max|O-O^H|=', max_herm, ', scale=', matrix_scale
               call g_logger%fatal('diagonalize_hamiltonian: ' // trim(herm_msg), __FILE__, __LINE__)
            end if
            if (use_generalized) call this%check_overlap_properties(ik, this%sk_overlap(:, :, ik))
         end do
      end if

      if (this%kanpur_diagnostics) call this%print_kanpur_mapping()
      if (trim(this%hamiltonian%magnetic_representation) == gbt_single_q) then
         write(herm_msg, '(A,ES12.4)') 'GBT pre-eigensolver max|H-H^H|=', max_herm_all
         call root_info(trim(herm_msg), __FILE__, __LINE__)
      end if

      if (allocated(this%eigenvalues)) deallocate(this%eigenvalues)
      if (allocated(this%eigenvectors)) deallocate(this%eigenvectors)
      allocate(this%eigenvalues(nmat, nk), this%eigenvectors(nmat, nmat, nk))
      call this%make_execution_backend()
      request%assemble_hamiltonian = .false.
      request%assemble_overlap = .false.
      request%solve_eigensystem = .true.
      request%generalized = use_generalized
      request%request_eigenvectors = .true.
      request%operator_generation = this%hamiltonian%operator_generation
      allocate(request%input_hamiltonian(nmat, nmat, nk), source=this%hk_bulk)
      if (use_generalized) allocate(request%input_overlap(nmat, nmat, nk), source=this%sk_overlap)
      call this%execution_backend%execute_batch(request, result)
      this%eigenvalues = result%eigenvalues
      this%eigenvectors = result%eigenvectors
      select type (backend => this%execution_backend)
      type is (cuda_reciprocal_backend)
         this%device_eigensystem_token = result%resident_token
      class default
         this%device_eigensystem_token = 0
      end select
      if (this%gamma_bounds_diagnostics) call this%run_gamma_bounds_diagnostics()
      if (this%hall_diag_experimental) call this%diagonalize_hall_experimental()
      call root_info('diagonalize_hamiltonian: Completed successfully', __FILE__, __LINE__)
#endif
   end subroutine diagonalize_hamiltonian

   !> @brief Print mapping diagnostics for the Kanpur generalized-overlap path.
   !> @param[in] this Reciprocal object containing basis and overlap metadata.
   module subroutine print_kanpur_mapping(this)
      class(reciprocal), intent(in) :: this
      call root_info('Kanpur mapping: reciprocal_mode=' // trim(this%reciprocal_mode), __FILE__, __LINE__)
      if (trim(this%reciprocal_mode) == 'generalized_overlap_proxy') then
         call root_info('Kanpur mapping: PROXY generalized solve H(k)c = E S_proxy(k)c, S_proxy from eeo + I.', __FILE__, __LINE__)
         if (rank == 0) call g_logger%warning('Kanpur mapping: PROXY mode is not a formal Kanpur LMTO overlap representation.', __FILE__, __LINE__)
      else if (trim(this%reciprocal_mode) == 'generalized_overlap_kanpur') then
         if (rank == 0) call g_logger%warning('Kanpur mapping: generalized_overlap_kanpur selected but not implemented; currently falling back to ham_only.', __FILE__, __LINE__)
      else
         call root_info('Kanpur mapping: Hamiltonian-only mode (TB-like).', __FILE__, __LINE__)
      end if
      call root_info('Kanpur mapping: non-orthogonality treatment is approximation-level diagnostic.', __FILE__, __LINE__)
   end subroutine print_kanpur_mapping

   !> @brief Require a Hermitian positive-definite overlap matrix.
   !> @param[in] this Reciprocal object providing diagnostic context.
   !> @param[in] ik k-point index associated with s_k.
   !> @param[in] s_k Overlap matrix to inspect.
   module subroutine check_overlap_properties(this, ik, s_k)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: ik
      complex(rp), dimension(:, :), intent(in) :: s_k
      integer :: i, j, n, info
      real(rp) :: max_herm
      complex(rp), allocatable :: chol(:, :)
      max_herm = 0.0_rp
      n = size(s_k, 1)
      do i = 1, n
         do j = 1, n
            max_herm = max(max_herm, abs(s_k(i, j) - conjg(s_k(j, i))))
         end do
      end do
      if (ik == 1 .or. max_herm > 1.0e-6_rp) then
         call g_logger%info('S(k) hermiticity check ik=' // trim(int2str(ik)) // ' max_diff=' // &
            trim(real2str(max_herm, '(ES12.4)')), __FILE__, __LINE__)
      end if
      if (max_herm > 1.0e-10_rp*max(1.0_rp, maxval(abs(s_k)))) then
         call g_logger%fatal('check_overlap_properties: overlap is not Hermitian at k='// &
                             trim(int2str(ik)), __FILE__, __LINE__)
      end if

      ! ZPOTRF applies the same positive-definiteness criterion required by
      ! ZHEGV. Check before the eigensolver so an incomplete or indefinite
      ! metric cannot degrade into per-k solver warnings.
      allocate(chol(n, n))
      chol = s_k
      call zpotrf('U', n, chol, n, info)
      deallocate(chol)
      if (info /= 0) then
         call g_logger%fatal('check_overlap_properties: overlap is not positive definite at k='// &
                             trim(int2str(ik))//'; leading minor='//trim(int2str(info)), &
                             __FILE__, __LINE__)
      end if
   end subroutine check_overlap_properties

   !> @brief Run H(Gamma) spectral-bound diagnostics.
   !> @param[inout] this Reciprocal object used to build and bound the Gamma matrix.
   module subroutine run_gamma_bounds_diagnostics(this)
      class(reciprocal), intent(inout) :: this
      complex(rp), allocatable :: h_gamma(:, :)
      real(rp) :: egmin, egmax
      type(bounds) :: bnd

      if (.not. allocated(this%hk_bulk)) return
      allocate(h_gamma(size(this%hk_bulk, 1), size(this%hk_bulk, 2)))
      h_gamma = this%hk_bulk(:, :, 1)
      call compute_spectrum_bounds(h_gamma, bnd, method='both', verbose=.false.)
      call g_logger%info('Gamma bounds diagnostic: Gershgorin [' // &
         trim(real2str(bnd%e_min_gershgorin, '(F12.6)')) // ', ' // &
         trim(real2str(bnd%e_max_gershgorin, '(F12.6)')) // ']', __FILE__, __LINE__)
      call this%diagonalize_hall_experimental() ! reuse exact solver utility logs for finite matrix check
      if (allocated(this%eigenvalues)) then
         egmin = minval(this%eigenvalues(:, 1))
         egmax = maxval(this%eigenvalues(:, 1))
         call g_logger%info('Gamma bounds diagnostic: eig(H(Gamma))=[' // &
            trim(real2str(egmin, '(F12.6)')) // ', ' // trim(real2str(egmax, '(F12.6)')) // ']', __FILE__, __LINE__)
      end if
      deallocate(h_gamma)
   end subroutine run_gamma_bounds_diagnostics

   !> @brief Diagonalize the experimental real-space HALL matrix.
   !> @details Builds a finite local-cluster matrix from HALL blocks and prints
   !>          eigenvalue diagnostics for exploratory comparison only.
   !> @param[inout] this Reciprocal object providing Hamiltonian and lattice state.
   module subroutine diagonalize_hall_experimental(this)
      class(reciprocal), intent(inout) :: this
      integer :: nsites, site_block, n, i, jsite, isite, ineigh, ia, ja, nr, info, lwork
      integer :: i_start, i_end, j_start, j_end
      complex(rp), allocatable :: hall_mat(:, :), work(:)
      real(rp), allocatable :: evals(:), rwork(:)

      if (.not. this%hall_diag_experimental) return
      if (this%control%calctype /= 'I') then
         call g_logger%info('HALL experimental diagonalization skipped: calctype is not impurity.', __FILE__, __LINE__)
         return
      end if

      nsites = this%lattice%nmax
      if (nsites <= 0) return
      site_block = this%max_orbs
      if (.not. allocated(this%hamiltonian%hall)) then
         call g_logger%warning('HALL experimental diagonalization skipped: HALL blocks are unavailable.', __FILE__, __LINE__)
         return
      end if
      if (site_block <= 0 .or. site_block /= nb .or. size(this%hamiltonian%hall, 1) /= site_block .or. &
          size(this%hamiltonian%hall, 2) /= site_block) then
         call g_logger%warning('HALL experimental diagonalization skipped: active site block is incompatible with HALL.', __FILE__, __LINE__)
         return
      end if
      n = nsites * site_block
      allocate(hall_mat(n, n), evals(n), rwork(max(1, 3*n - 2)))
      hall_mat = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, nsites
         ia = isite
         nr = this%lattice%nn(ia, 1)
         i_start = (isite - 1) * site_block + 1
         i_end = isite * site_block
         do ineigh = 1, nr
            if (ineigh == 1) then
               jsite = isite
            else
               ja = this%lattice%nn(ia, ineigh)
               jsite = ja
               if (jsite < 1 .or. jsite > nsites) cycle
            end if
            j_start = (jsite - 1) * site_block + 1
            j_end = jsite * site_block
            hall_mat(i_start:i_end, j_start:j_end) = hall_mat(i_start:i_end, j_start:j_end) + this%hamiltonian%hall(:, :, ineigh, isite)
         end do
      end do
      allocate(work(1))
      call zheev('N', 'U', n, hall_mat, n, evals, work, -1, rwork, info)
      lwork = max(1, int(real(work(1))))
      deallocate(work)
      allocate(work(lwork))
      call zheev('N', 'U', n, hall_mat, n, evals, work, lwork, rwork, info)
      if (info == 0) then
         call g_logger%info('HALL experimental eig range: [' // trim(real2str(minval(evals), '(F12.6)')) // ', ' // &
            trim(real2str(maxval(evals), '(F12.6)')) // '] (diagnostic only)', __FILE__, __LINE__)
      else
         call g_logger%warning('HALL experimental diagonalization failed, info=' // trim(int2str(info)), __FILE__, __LINE__)
      end if
      deallocate(hall_mat, evals, rwork, work)
   end subroutine diagonalize_hall_experimental

   !> @brief Calculate and write a band structure along a crystal-specific path.
   !> @param[inout] this Reciprocal object used to generate path eigenvalues.
   !> @param[in] ham Hamiltonian object for k-space construction.
   !> @param[in] crystal_type Crystal/path selector.
   !> @param[in] npts_per_segment Number of interpolation points per path segment.
   !> @param[in] output_file Optional band-output file name.
   module subroutine calculate_band_structure(this, ham, crystal_type, npts_per_segment, output_file)
      class(reciprocal), intent(inout) :: this
      class(hamiltonian), intent(in) :: ham
      character(len=*), intent(in) :: crystal_type
      integer, intent(in) :: npts_per_segment
      character(len=*), intent(in), optional :: output_file
      ! Local variables
      character(len=256) :: filename
      integer :: unit, i, j, nmat
      character(len=100) :: fmt_str

      call g_logger%info('calculate_band_structure: Starting band structure calculation', __FILE__, __LINE__)

      ! Generate k-path using spglib-based canonical path generation
      call this%symmetry_analysis%generate_canonical_kpath(npts_per_segment)

      ! Copy k-path data from symmetry analysis to reciprocal object
      if (allocated(this%symmetry_analysis%k_path)) then
         ! Copy k-path arrays
         this%nk_path = this%symmetry_analysis%nk_path
         if (allocated(this%k_path)) deallocate(this%k_path)
         if (allocated(this%k_labels)) deallocate(this%k_labels)
         if (allocated(this%k_distances)) deallocate(this%k_distances)
         
         allocate(this%k_path(3, this%nk_path))
         allocate(this%k_labels(this%nk_path))
         allocate(this%k_distances(this%nk_path))
         
         this%k_path = this%symmetry_analysis%k_path
         this%k_labels = this%symmetry_analysis%k_labels
         this%k_distances = this%symmetry_analysis%k_distances
         
         call g_logger%info('calculate_band_structure: Copied k-path with ' // &
                           trim(int2str(this%nk_path)) // ' points from symmetry analysis', __FILE__, __LINE__)
      else
         call g_logger%error('calculate_band_structure: Symmetry analysis k_path not allocated!', __FILE__, __LINE__)
         return
      end if

      ! Debug logging for k-path
      call g_logger%info('calculate_band_structure: After k-path generation, nk_path = ' // &
                        trim(int2str(this%nk_path)), __FILE__, __LINE__)
      if (allocated(this%k_path)) then
         call g_logger%info('calculate_band_structure: k_path allocated with dimensions ' // &
                           trim(int2str(size(this%k_path,1))) // 'x' // &
                           trim(int2str(size(this%k_path,2))), __FILE__, __LINE__)
      else
         call g_logger%error('calculate_band_structure: k_path not allocated after k-path generation!', __FILE__, __LINE__)
      end if

      ! Build k-space Hamiltonian for the k-path
      call this%build_kspace_hamiltonian()

      ! Diagonalize Hamiltonian along k-path
      call this%diagonalize_hamiltonian()

      ! Set output filename
      filename = 'band_structure.dat'
      if (present(output_file)) filename = output_file

      ! Write band structure to file
      open(newunit=unit, file=trim(filename), status='replace', action='write')
      
      ! Write header
      write(unit, '(A)') '# Band structure calculation'
      write(unit, '(A)') '# Crystal type: ' // trim(crystal_type)
      write(unit, '(A,I0)') '# Number of k-points: ', this%nk_path
      nmat = size(this%eigenvalues, 1)
      write(unit, '(A,I0)') '# Number of bands: ', nmat
      write(unit, '(A)') '# Format: k_distance, eigenvalue_1, eigenvalue_2, ...'
      write(unit, '(A)') '#'

      ! Write data
      write(fmt_str, '(A,I0,A)') '(', nmat+1, '(ES16.8E3,1X))'
      do i = 1, this%nk_path
         write(unit, fmt_str) this%k_distances(i), (this%eigenvalues(j, i), j = 1, nmat)
      end do

      close(unit)

      call g_logger%info('calculate_band_structure: Band structure written to ' // trim(filename), __FILE__, __LINE__)

      ! Also write k-path information
      filename = 'kpath_info.dat'
      if (present(output_file)) then
         ! Replace extension with _kpath.dat
         i = index(output_file, '.', back=.true.)
         if (i > 0) then
            filename = output_file(1:i-1) // '_kpath.dat'
         else
            filename = trim(output_file) // '_kpath.dat'
         end if
      end if

      open(newunit=unit, file=trim(filename), status='replace', action='write')
      write(unit, '(A)') '# K-path information'
      write(unit, '(A)') '# Format: k_distance, kx, ky, kz, label'
      write(unit, '(A)') '#'
      
      do i = 1, this%nk_path
         if (trim(this%k_labels(i)) /= '') then
            write(unit, '(4(ES16.8E3,1X),A)') this%k_distances(i), this%k_path(:, i), ' # ' // trim(this%k_labels(i))
         else
            write(unit, '(4(ES16.8E3,1X))') this%k_distances(i), this%k_path(:, i)
         end if
      end do

      close(unit)

      call g_logger%info('calculate_band_structure: K-path information written to ' // trim(filename), __FILE__, __LINE__)
   end subroutine calculate_band_structure

   !> @brief Calculate and write a symmetry-derived band structure.
   !> @details Uses symmetry analysis or custom path settings to generate the
   !>          high-symmetry path before building and diagonalizing H(k).
   !> @param[inout] this Reciprocal object used to generate path eigenvalues.
   !> @param[in] ham Hamiltonian object for k-space construction.
   !> @param[in] npts_per_segment Optional number of interpolation points per segment.
   !> @param[in] output_file Optional band-output file name.
   module subroutine calculate_band_structure_auto(this, ham, npts_per_segment, output_file)
      class(reciprocal), intent(inout) :: this
      class(hamiltonian), intent(in) :: ham
      integer, intent(in), optional :: npts_per_segment
      character(len=*), intent(in), optional :: output_file
      
      ! Local variables
      integer :: npts
      character(len=256) :: filename
      integer :: unit, i, j, nmat
      character(len=100) :: fmt_str

      npts = 40  ! Default
      if (present(npts_per_segment)) npts = npts_per_segment

      call g_logger%info('calculate_band_structure_auto: Starting automatic band structure calculation', __FILE__, __LINE__)

      ! Generate k-path automatically using symmetry analysis
      call this%symmetry_analysis%generate_automatic_kpath(npts)

      ! Copy k-path data from symmetry analysis to reciprocal object
      if (allocated(this%symmetry_analysis%k_path)) then
         this%nk_path = this%symmetry_analysis%nk_path
         
         if (allocated(this%k_path)) deallocate(this%k_path)
         if (allocated(this%k_labels)) deallocate(this%k_labels)
         if (allocated(this%k_distances)) deallocate(this%k_distances)
         
         allocate(this%k_path(3, this%nk_path))
         allocate(this%k_labels(this%nk_path))
         allocate(this%k_distances(this%nk_path))
         
         this%k_path = this%symmetry_analysis%k_path
         this%k_labels = this%symmetry_analysis%k_labels
         this%k_distances = this%symmetry_analysis%k_distances
         
         call g_logger%info('calculate_band_structure_auto: K-path has ' // &
                           trim(int2str(this%nk_path)) // ' points', __FILE__, __LINE__)
      else
         call g_logger%error('calculate_band_structure_auto: Failed to generate k-path', __FILE__, __LINE__)
         return
      end if

      ! Build k-space Hamiltonian for the k-path
      call this%build_kspace_hamiltonian()

      ! Diagonalize Hamiltonian along k-path
      call this%diagonalize_hamiltonian()

      ! Set output filename
      filename = 'band_structure.dat'
      if (present(output_file)) filename = output_file

      ! Write band structure to file
      open(newunit=unit, file=trim(filename), status='replace', action='write')
      
      ! Write header with auto-detected information
   write(unit, '(A)') '# Band structure calculation (automatic k-path)'
#ifdef USE_SPGLIB
   if (this%symmetry_analysis%spglib%is_available()) then
      write(unit, '(A)') '# Space group: ' // trim(this%symmetry_analysis%spglib%get_space_group_symbol()) // &
               ' (#' // trim(int2str(this%symmetry_analysis%spglib%get_space_group_number())) // ')'
      write(unit, '(A)') '# Crystal system: ' // trim(this%symmetry_analysis%spglib%get_crystal_system_name())
   end if
#endif
      write(unit, '(A,I0)') '# Number of k-points: ', this%nk_path
      nmat = size(this%eigenvalues_path, 1)
      write(unit, '(A,I0)') '# Number of bands: ', nmat
      write(unit, '(A)') '# Format: k_distance, eigenvalue_1, eigenvalue_2, ...'
      write(unit, '(A)') '#'

      ! Write data
      write(fmt_str, '(A,I0,A)') '(', nmat+1, '(ES16.8E3,1X))'
      do i = 1, this%nk_path
         write(unit, fmt_str) this%k_distances(i), (this%eigenvalues_path(j, i), j = 1, nmat)
      end do

      close(unit)

      call g_logger%info('calculate_band_structure_auto: Band structure written to ' // trim(filename), __FILE__, __LINE__)

      ! Write k-path information
      filename = 'kpath_info.dat'
      if (present(output_file)) then
         i = index(output_file, '.', back=.true.)
         if (i > 0) then
            filename = output_file(1:i-1) // '_kpath.dat'
         else
            filename = trim(output_file) // '_kpath.dat'
         end if
      end if

      open(newunit=unit, file=trim(filename), status='replace', action='write')
      write(unit, '(A)') '# K-path information'
      write(unit, '(A)') '# Format: k_distance, kx, ky, kz, label'
      write(unit, '(A)') '#'
      
      do i = 1, this%nk_path
         if (trim(this%k_labels(i)) /= '') then
            write(unit, '(4(ES16.8E3,1X),A)') this%k_distances(i), this%k_path(:, i), ' # ' // trim(this%k_labels(i))
         else
            write(unit, '(4(ES16.8E3,1X))') this%k_distances(i), this%k_path(:, i)
         end if
      end do

      close(unit)

      call g_logger%info('calculate_band_structure_auto: K-path information written to ' // trim(filename), __FILE__, __LINE__)
   end subroutine calculate_band_structure_auto

   !> @brief Generate a symmetry-reduced k-point mesh.
   !> @details Builds irreducible k-points, weights, and full/irreducible maps
   !>          using the configured symmetry and time-reversal settings.
   !> @param[inout] this Reciprocal object receiving reduced-mesh data.
   !> @param[in] mesh_dims Full mesh dimensions.
   !> @param[in] use_shift Optional flag selecting shifted-grid generation.
   module subroutine generate_reduced_kpoint_mesh(this, mesh_dims, use_shift)
      class(reciprocal), intent(inout) :: this
      integer, intent(in) :: mesh_dims(3)
      logical, intent(in), optional :: use_shift
      integer :: shift(3)
      integer :: num_ir_kpoints
      logical :: do_shift, effective_time_reversal
      real(rp), allocatable :: kpoints_frac(:,:), weights(:)
      integer :: i
      integer, allocatable :: full_to_irred(:), irred_to_full(:)

      do_shift = .false.
      if (present(use_shift)) do_shift = use_shift

      if (this%has_nonzero_q_gbt()) then
         ! WP0's full-BZ guard stays the default and the oracle every other
         ! policy is checked against (WP8: "keep full BZ as the oracle").
         ! q_symmetry_policy is an explicit opt-in to reduce by the little
         ! group of q_ss instead; anything else (including the 'full_bz'
         ! default) reproduces the WP0 behaviour bit-for-bit.
         if (trim(this%q_symmetry_policy) == 'little_group' .or. &
             trim(this%q_symmetry_policy) == 'little_group_common') then
            call this%generate_little_group_kpoint_mesh(mesh_dims, use_shift)
            return
         end if
         this%use_symmetry_reduction = .false.
         this%use_time_reversal = .false.
         call root_info('generate_reduced_kpoint_mesh: nonzero-q GBT forces the full chemical BZ.', __FILE__, __LINE__)
         call this%generate_mp_mesh()
         return
      end if

      if (do_shift) then
         shift = [1, 1, 1]  ! Offset by half a mesh spacing
      else
         shift = [0, 0, 0]  ! No offset
      end if

      ! Store mesh dimensions
      this%nk_mesh = mesh_dims
      effective_time_reversal = this%use_time_reversal
      if (associated(this%control)) then
         if (this%control%nsp >= 3 .and. effective_time_reversal) then
            effective_time_reversal = .false.
            call g_logger%info('generate_reduced_kpoint_mesh: Disabled time-reversal reduction for non-collinear calculation; spatial symmetry reduction remains enabled.', __FILE__, __LINE__)
         end if
      end if

#ifdef USE_SPGLIB
   if (.not. this%symmetry_analysis%spglib%is_available()) then
      call g_logger%warning('generate_reduced_kpoint_mesh: spglib not available, using full mesh', __FILE__, __LINE__)
      call this%generate_mp_mesh()  ! Fall back to regular MP mesh
      return
   end if

   ! Get irreducible k-points and weights from spglib
   num_ir_kpoints = this%symmetry_analysis%spglib%get_reduced_kpoint_mesh_with_points( &
                           mesh_dims, shift, kpoints_frac, weights, effective_time_reversal, &
                           full_to_irred, irred_to_full)
#else
   ! SPGLIB not enabled at compile time: fall back to full mesh
   call g_logger%warning('generate_reduced_kpoint_mesh: spglib support was not compiled in, using full mesh', __FILE__, __LINE__)
   call this%generate_mp_mesh()
   return
#endif

      if (num_ir_kpoints == 0) then
         call g_logger%error('generate_reduced_kpoint_mesh: Failed to get k-points from spglib', __FILE__, __LINE__)
         call this%generate_mp_mesh()  ! Fall back to regular MP mesh
         return
      end if

      ! Store k-point information
      this%nk_total = num_ir_kpoints

      ! Allocate k-point arrays
#ifdef USE_SAFE_ALLOC
      if (allocated(this%k_points)) call g_safe_alloc%deallocate('reciprocal.k_points', this%k_points)
      if (allocated(this%k_weights)) call g_safe_alloc%deallocate('reciprocal.k_weights', this%k_weights)
      call g_safe_alloc%allocate_real('reciprocal.k_points', this%k_points, (/3, this%nk_total/))
      call g_safe_alloc%allocate_real('reciprocal.k_weights', this%k_weights, (/this%nk_total/))
#else
      if (allocated(this%k_points)) deallocate(this%k_points)
      if (allocated(this%k_weights)) deallocate(this%k_weights)
      allocate(this%k_points(3, this%nk_total))
      allocate(this%k_weights(this%nk_total))
#endif

      ! Copy k-points and weights
      this%k_points = kpoints_frac
      this%k_weights = weights
      if (allocated(this%full_to_irred_k)) deallocate(this%full_to_irred_k)
      if (allocated(this%irred_to_full_k)) deallocate(this%irred_to_full_k)
      allocate(this%full_to_irred_k(size(full_to_irred)))
      allocate(this%irred_to_full_k(size(irred_to_full)))
      this%full_to_irred_k = full_to_irred
      this%irred_to_full_k = irred_to_full

      ! Clean up temporary arrays
      deallocate(kpoints_frac, weights, full_to_irred, irred_to_full)
      call setup_k_mesh_distribution(this, this%nk_total, .false.)

      call root_info('generate_reduced_kpoint_mesh: Generated ' // trim(int2str(this%nk_total)) // &
                     ' irreducible k-points from ' // trim(int2str(product(mesh_dims))) // ' total points', &
                     __FILE__, __LINE__)
      call root_info('generate_reduced_kpoint_mesh: Reduction factor: ' // &
                     trim(real2str(real(product(mesh_dims), rp)/real(this%nk_total, rp), '(F6.2)')) // 'x', &
                     __FILE__, __LINE__)
      
      ! Verify weights sum to 1
      if (abs(sum(this%k_weights) - 1.0_rp) > 1.0e-6_rp) then
         if (this%strict_symmetry_checks) then
            call g_logger%fatal('generate_reduced_kpoint_mesh: K-point weights sum to ' // &
                               trim(real2str(sum(this%k_weights), '(F12.8)')) // ' (expected 1.0)', &
                               __FILE__, __LINE__)
         else
            call g_logger%warning('generate_reduced_kpoint_mesh: K-point weights sum to ' // &
                                 trim(real2str(sum(this%k_weights), '(F12.8)')) // ' (should be 1.0)', &
                                 __FILE__, __LINE__)
         end if
      end if
      if (allocated(this%full_to_irred_k)) then
         if (any(this%full_to_irred_k < 1) .or. any(this%full_to_irred_k > this%nk_total)) then
            if (this%strict_symmetry_checks) then
               call g_logger%fatal('generate_reduced_kpoint_mesh: Invalid full_to_irred mapping detected', __FILE__, __LINE__)
            else
               call g_logger%warning('generate_reduced_kpoint_mesh: Invalid full_to_irred mapping detected', __FILE__, __LINE__)
            end if
         end if
      end if
      call this%validate_symmetry_kmap('generate_reduced_kpoint_mesh')
      if (this%dump_symmetry_kmap) call this%write_symmetry_kmap_dump('symmetry_kmap.dat')
   end subroutine generate_reduced_kpoint_mesh

   !> @brief Generate a k-mesh reduced by the little group common to one or
   !>        more q-points. See the interface docstring in reciprocal.f90.
   module subroutine generate_little_group_kpoint_mesh(this, mesh_dims, use_shift, q_list_cart)
      class(reciprocal), intent(inout) :: this
      integer, intent(in) :: mesh_dims(3)
      logical, intent(in), optional :: use_shift
      real(rp), intent(in), optional :: q_list_cart(:, :)
      integer :: shift(3)
      integer :: num_ir_kpoints
      logical :: do_shift, effective_time_reversal
      real(rp), allocatable :: kpoints_frac(:,:), weights(:), q_frac(:,:)
      integer, allocatable :: full_to_irred(:), irred_to_full(:)
      integer :: iq, num_q

      do_shift = .false.
      if (present(use_shift)) do_shift = use_shift

      if (do_shift) then
         shift = [1, 1, 1]
      else
         shift = [0, 0, 0]
      end if

      this%nk_mesh = mesh_dims
      effective_time_reversal = this%use_time_reversal
      if (associated(this%control)) then
         if (this%control%nsp >= 3 .and. effective_time_reversal) then
            effective_time_reversal = .false.
            call g_logger%info('generate_little_group_kpoint_mesh: Disabled time-reversal reduction for '// &
                               'non-collinear calculation; the little-group spatial reduction remains enabled.', &
                               __FILE__, __LINE__)
         end if
      end if

      ! q_ss is Cartesian in units of 2*pi/alat (source/hamiltonian_build.f90);
      ! convert to fractional reciprocal-lattice coordinates for spglib using
      ! the same, already-locked convention as the frozen_magnon
      ! q_coordinates='direct' round trip (tests/unit/test_qss_theta_conventions.f90,
      ! source/calculation.f90): cart_to_direct = transpose(lattice%a).
      if (present(q_list_cart)) then
         num_q = size(q_list_cart, 2)
         allocate(q_frac(3, num_q))
         do iq = 1, num_q
            q_frac(:, iq) = matmul(transpose(this%lattice%a), q_list_cart(:, iq))
         end do
      else
         num_q = 1
         allocate(q_frac(3, 1))
         q_frac(:, 1) = matmul(transpose(this%lattice%a), this%hamiltonian%q_ss)
      end if

#ifdef USE_SPGLIB
      if (.not. this%symmetry_analysis%spglib%is_available()) then
         call g_logger%info('generate_little_group_kpoint_mesh: spglib not available; falling back to the '// &
                            'full k-mesh for q_ss != 0 (little group could not be determined).', &
                            __FILE__, __LINE__)
         deallocate(q_frac)
         call this%generate_mp_mesh()
         return
      end if

      num_ir_kpoints = this%symmetry_analysis%spglib%get_little_group_kpoint_mesh_with_points( &
                           mesh_dims, shift, q_frac, kpoints_frac, weights, effective_time_reversal, &
                           full_to_irred, irred_to_full)
#else
      call g_logger%info('generate_little_group_kpoint_mesh: spglib support was not compiled in; '// &
                         'falling back to the full k-mesh for q_ss != 0 (little group could not '// &
                         'be determined).', __FILE__, __LINE__)
      deallocate(q_frac)
      call this%generate_mp_mesh()
      return
#endif
      deallocate(q_frac)

      if (num_ir_kpoints <= 0) then
         call g_logger%info('generate_little_group_kpoint_mesh: little group could not be '// &
                            'determined; falling back to the full k-mesh.', __FILE__, __LINE__)
         call this%generate_mp_mesh()
         return
      end if

      this%nk_total = num_ir_kpoints
      if (allocated(this%k_points)) deallocate(this%k_points)
      if (allocated(this%k_weights)) deallocate(this%k_weights)
      allocate(this%k_points(3, this%nk_total))
      allocate(this%k_weights(this%nk_total))
      this%k_points = kpoints_frac
      this%k_weights = weights
      if (allocated(this%full_to_irred_k)) deallocate(this%full_to_irred_k)
      if (allocated(this%irred_to_full_k)) deallocate(this%irred_to_full_k)
      allocate(this%full_to_irred_k(size(full_to_irred)))
      allocate(this%irred_to_full_k(size(irred_to_full)))
      this%full_to_irred_k = full_to_irred
      this%irred_to_full_k = irred_to_full
      deallocate(kpoints_frac, weights, full_to_irred, irred_to_full)
      call setup_k_mesh_distribution(this, this%nk_total, .false.)

      call root_info('generate_little_group_kpoint_mesh: q_ss != 0; reduced by the little group common to '// &
                     trim(int2str(num_q))//' q-point(s) to '//trim(int2str(this%nk_total))// &
                     ' irreducible k-points from '//trim(int2str(product(mesh_dims)))//' total points ('// &
                     trim(real2str(real(product(mesh_dims), rp)/real(this%nk_total, rp), '(F6.2)'))// &
                     'x reduction).', __FILE__, __LINE__)

      if (abs(sum(this%k_weights) - 1.0_rp) > 1.0e-6_rp) then
         if (this%strict_symmetry_checks) then
            call g_logger%fatal('generate_little_group_kpoint_mesh: K-point weights sum to ' // &
                               trim(real2str(sum(this%k_weights), '(F12.8)')) // ' (expected 1.0)', &
                               __FILE__, __LINE__)
         else
            call g_logger%warning('generate_little_group_kpoint_mesh: K-point weights sum to ' // &
                                 trim(real2str(sum(this%k_weights), '(F12.8)')) // ' (should be 1.0)', &
                                 __FILE__, __LINE__)
         end if
      end if
      if (allocated(this%full_to_irred_k)) then
         if (any(this%full_to_irred_k < 1) .or. any(this%full_to_irred_k > this%nk_total)) then
            if (this%strict_symmetry_checks) then
               call g_logger%fatal('generate_little_group_kpoint_mesh: Invalid full_to_irred mapping detected', &
                                   __FILE__, __LINE__)
            else
               call g_logger%warning('generate_little_group_kpoint_mesh: Invalid full_to_irred mapping detected', &
                                     __FILE__, __LINE__)
            end if
         end if
      end if
      call this%validate_symmetry_kmap('generate_little_group_kpoint_mesh')
      if (this%dump_symmetry_kmap) call this%write_symmetry_kmap_dump('symmetry_kmap.dat')
   end subroutine generate_little_group_kpoint_mesh

   !> @brief Is this a GBT run with a nonzero q, either the current single
   !>        hamiltonian%q_ss (has_nonzero_q_gbt) or, if given, any column of
   !>        an explicitly declared q_list_cart (the 'little_group_common'
   !>        pre-loop case, called while q_ss is still the reference point).
   logical function is_declared_nonzero_q_gbt(this, q_list_cart) result(nonzero)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in), optional :: q_list_cart(:, :)

      nonzero = this%has_nonzero_q_gbt()
      if (nonzero .or. .not. present(q_list_cart)) return
      if (.not. associated(this%hamiltonian)) return
      if (trim(this%hamiltonian%magnetic_representation) /= gbt_single_q) return
      nonzero = any(abs(q_list_cart) > 1.0e-12_rp)
   end function is_declared_nonzero_q_gbt

   !> @brief Ensure k_points/k_weights match the current cache key, rebuilding
   !>        only if not. See the interface docstring in reciprocal.f90.
   module subroutine ensure_kpoint_mesh(this, mesh_dims, use_shift, q_list_cart)
      class(reciprocal), intent(inout) :: this
      integer, intent(in) :: mesh_dims(3)
      logical, intent(in), optional :: use_shift
      real(rp), intent(in), optional :: q_list_cart(:, :)
      logical :: do_shift
      real(rp) :: offset(3)
      real(rp), allocatable :: q_now(:, :)
      logical :: key_matches

      do_shift = .false.
      if (present(use_shift)) do_shift = use_shift
      offset = 0.0_rp
      if (do_shift) offset = 0.5_rp

      ! The q-part of the key: the declared sweep list under
      ! 'little_group_common', otherwise the single current q_ss (zero for
      ! q=0/non-GBT runs, so ordinary callers get a key that only tracks
      ! mesh/offset/lattice/policy, matching their actual dependence).
      if (present(q_list_cart)) then
         allocate(q_now(3, size(q_list_cart, 2)))
         q_now = q_list_cart
      else
         allocate(q_now(3, 1))
         q_now(:, 1) = 0.0_rp
         if (associated(this%hamiltonian)) q_now(:, 1) = this%hamiltonian%q_ss
      end if

      key_matches = .false.
      if (this%mesh_cache_valid .and. allocated(this%k_points)) then
         key_matches = all(this%mesh_cache_dims == mesh_dims) .and. &
                       all(abs(this%mesh_cache_offset - offset) < 1.0e-12_rp) .and. &
                       trim(this%mesh_cache_policy) == trim(this%q_symmetry_policy) .and. &
                       allocated(this%mesh_cache_q)
         if (key_matches .and. associated(this%lattice)) then
            key_matches = all(abs(this%mesh_cache_lattice - this%lattice%a) < 1.0e-12_rp)
         end if
         if (key_matches) then
            key_matches = size(this%mesh_cache_q, 2) == size(q_now, 2)
            if (key_matches) key_matches = all(abs(this%mesh_cache_q - q_now) < 1.0e-12_rp)
         end if
      end if

      if (key_matches) then
         deallocate(q_now)
         return
      end if

      ! Cache key changed (or no mesh yet): rebuild through the ordinary
      ! dispatch, which itself honours q_symmetry_policy and has_nonzero_q_gbt.
      ! has_nonzero_q_gbt() alone only sees the single CURRENT hamiltonian%q_ss
      ! -- for a pre-loop 'little_group_common' call the caller's q_list_cart
      ! is nonzero while the reference q_ss is still (0,0,0), so a caller-
      ! supplied nonzero q-set must also count as "nonzero q GBT" here or the
      ! sweep would silently fall through to the ordinary q=0 point-group mesh
      ! (exactly the "wrong subgroup" failure mode WP8 exists to remove).
      if (is_declared_nonzero_q_gbt(this, q_list_cart) .and. &
          (trim(this%q_symmetry_policy) == 'little_group' .or. &
           trim(this%q_symmetry_policy) == 'little_group_common')) then
         call this%generate_little_group_kpoint_mesh(mesh_dims, use_shift, q_list_cart)
      else if (this%use_symmetry_reduction .and. .not. this%has_nonzero_q_gbt()) then
         call this%generate_reduced_kpoint_mesh(mesh_dims, do_shift)
      else
         call this%generate_mp_mesh()
      end if

      ! The mesh (and therefore every eigensystem/DOS/density object derived
      ! from it) just changed identity, independent of whether the shared
      ! Hamiltonian operator_generation moved -- invalidate explicitly rather
      ! than relying on the next build_bulkham/diagonalize call to notice.
      call this%invalidate_spectral_cache()

      this%mesh_cache_dims = mesh_dims
      this%mesh_cache_offset = offset
      this%mesh_cache_policy = this%q_symmetry_policy
      if (associated(this%lattice)) this%mesh_cache_lattice = this%lattice%a
      if (allocated(this%mesh_cache_q)) deallocate(this%mesh_cache_q)
      allocate(this%mesh_cache_q(3, size(q_now, 2)))
      this%mesh_cache_q = q_now
      this%mesh_cache_valid = .true.
      deallocate(q_now)
   end subroutine ensure_kpoint_mesh

   !> @brief Write a complex matrix and its k-point label to a text file.
   !> @param[in] matrix Matrix to dump.
   !> @param[in] filename Output file name.
   !> @param[in] k_point k-point associated with the matrix.
   module subroutine dump_complex_matrix(matrix, filename, k_point)
      complex(rp), dimension(:, :), intent(in) :: matrix
      character(len=*), intent(in) :: filename
      real(rp), dimension(3), intent(in) :: k_point
      ! Local variables
      integer :: unit, i, j, n_rows, n_cols
      character(len=100) :: fmt_str

      n_rows = size(matrix, 1)
      n_cols = size(matrix, 2)

      open(newunit=unit, file=trim(filename), status='replace', action='write')

      ! Write header
      write(unit, '(A)') '# Complex matrix dump for debugging'
      write(unit, '(A,3F12.8)') '# k-point: ', k_point
      write(unit, '(A,I0,A,I0)') '# Matrix dimensions: ', n_rows, ' x ', n_cols
      write(unit, '(A)') '# Format: row, col, real_part, imag_part'
      write(unit, '(A)') '#'

      ! Write matrix elements
      do i = 1, n_rows
         do j = 1, n_cols
            write(unit, '(2I6,2ES20.12E3)') i, j, real(matrix(i,j)), aimag(matrix(i,j))
         end do
      end do

      close(unit)
   end subroutine dump_complex_matrix

end submodule reciprocal_bands
