submodule (reciprocal_mod) reciprocal_bands
   implicit none

contains

   !> @brief Diagonalize the active k-space Hamiltonians.
   !> @details Solves standard or generalized eigenproblems depending on
   !>          reciprocal_mode and stores eigenvalues/eigenvectors for bands and DOS.
   !> @param[inout] this Reciprocal object receiving eigenvalue/eigenvector arrays.
   module subroutine diagonalize_hamiltonian(this)
      class(reciprocal), intent(inout) :: this
      
      ! Local variables
      integer :: nk, ik, nmat, lwork, info, mode_fail_count
      complex(rp), dimension(:, :), allocatable :: h_k_copy
      complex(rp), dimension(:, :), allocatable :: s_k_copy
      real(rp), dimension(:), allocatable :: eigenvals
      complex(rp), dimension(:), allocatable :: work_complex
      real(rp), dimension(:), allocatable :: rwork
      character(len=100) :: info_msg
      logical :: use_generalized
      real(rp) :: max_herm, matrix_scale
      character(len=256) :: herm_msg

      ! Check prerequisites
      if (.not. allocated(this%hk_bulk)) then
         call g_logger%error('diagonalize_hamiltonian: hk_bulk not built - call build_kspace_hamiltonian first', &
                            __FILE__, __LINE__)
         return
      end if

      ! Get dimensions from hk_bulk
      nmat = size(this%hk_bulk, 1)
      nk = size(this%hk_bulk, 3)
      
      call root_info('diagonalize_hamiltonian: Diagonalizing ' // trim(int2str(nk)) // ' k-points', __FILE__, __LINE__)
      call root_info('diagonalize_hamiltonian: Matrix size = ' // &
                        trim(int2str(nmat)) // ' x ' // trim(int2str(nmat)), __FILE__, __LINE__)

      use_generalized = trim(this%reciprocal_mode) == 'generalized_overlap_proxy'
      if (use_generalized) then
         if (.not. allocated(this%sk_overlap)) call this%build_kspace_overlap()
         if (.not. allocated(this%sk_overlap)) then
            call g_logger%warning('diagonalize_hamiltonian: S(k) unavailable, falling back to ham_only.', __FILE__, __LINE__)
            use_generalized = .false.
         end if
      end if

      ! zheev/zhegv read only one triangle.  Check the completed matrices
      ! first, otherwise a broken lower triangle (notably the old finite-q GBT
      ! reconstruction) would be silently discarded by LAPACK.
      do ik = 1, nk
         max_herm = maxval(abs(this%hk_bulk(:, :, ik) - transpose(conjg(this%hk_bulk(:, :, ik)))) )
         matrix_scale = max(1.0_rp, maxval(abs(this%hk_bulk(:, :, ik))))
         if (max_herm > 1.0e-10_rp*matrix_scale) then
            write(herm_msg, '(A,I0,A,ES12.4,A,ES12.4)') 'H(k) is non-Hermitian before eigensolution at k=', ik, &
               ': max|H-H^H|=', max_herm, ', scale=', matrix_scale
            call g_logger%fatal('diagonalize_hamiltonian: '//trim(herm_msg), __FILE__, __LINE__)
         end if
      end do
      if (allocated(this%sk_overlap)) then
         do ik = 1, size(this%sk_overlap, 3)
            max_herm = maxval(abs(this%sk_overlap(:, :, ik) - transpose(conjg(this%sk_overlap(:, :, ik)))) )
            matrix_scale = max(1.0_rp, maxval(abs(this%sk_overlap(:, :, ik))))
            if (max_herm > 1.0e-10_rp*matrix_scale) then
               write(herm_msg, '(A,I0,A,ES12.4,A,ES12.4)') 'O(k) is non-Hermitian before eigensolution at k=', ik, &
                  ': max|O-O^H|=', max_herm, ', scale=', matrix_scale
               call g_logger%fatal('diagonalize_hamiltonian: '//trim(herm_msg), __FILE__, __LINE__)
            end if
         end do
      end if

      if (this%kanpur_diagnostics) call this%print_kanpur_mapping()

      ! Allocate eigenvalue and eigenvector storage
      if (allocated(this%eigenvalues)) deallocate(this%eigenvalues)
      if (allocated(this%eigenvectors)) deallocate(this%eigenvectors)
      allocate(this%eigenvalues(nmat, nk))
      allocate(this%eigenvectors(nmat, nmat, nk))
      mode_fail_count = 0

      ! Parallel diagonalization over k-points
      ! Each thread needs its own LAPACK workspace
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP& PRIVATE(ik, h_k_copy, s_k_copy, eigenvals, work_complex, rwork, lwork, info, info_msg) &
!$OMP& IF(nk > 10)
      
      ! Allocate thread-private work arrays
      allocate(h_k_copy(nmat, nmat))
      allocate(s_k_copy(nmat, nmat))
      allocate(eigenvals(nmat))
      allocate(rwork(3*nmat - 2))
      
      ! Query optimal LAPACK workspace size
      lwork = -1
      allocate(work_complex(1))
      if (use_generalized) then
         call zhegv(1, 'V', 'U', nmat, h_k_copy, nmat, s_k_copy, nmat, eigenvals, work_complex, lwork, rwork, info)
      else
         call zheev('V', 'U', nmat, h_k_copy, nmat, eigenvals, work_complex, lwork, rwork, info)
      end if
      lwork = int(real(work_complex(1)))
      deallocate(work_complex)
      allocate(work_complex(lwork))

      ! Loop over all k-points - use pre-computed hk_bulk
!$OMP DO SCHEDULE(DYNAMIC, 10)
      do ik = 1, nk
         ! Copy H(k) from pre-computed array
         h_k_copy = this%hk_bulk(:, :, ik)

         if (use_generalized) then
            s_k_copy = this%sk_overlap(:, :, ik)
            call this%check_overlap_properties(ik, s_k_copy)
            call zhegv(1, 'V', 'U', nmat, h_k_copy, nmat, s_k_copy, nmat, eigenvals, work_complex, lwork, rwork, info)
         else
            ! Diagonalize H(k) using LAPACK ZHEEV
            call zheev('V', 'U', nmat, h_k_copy, nmat, eigenvals, work_complex, lwork, rwork, info)
         end if
         
         if (info /= 0) then
            write(info_msg, '(A,I0,A,I0)') 'Diagonalization failed at k-point ', ik, ', info = ', info
            call g_logger%error('diagonalize_hamiltonian: ' // trim(info_msg), __FILE__, __LINE__)
!$OMP CRITICAL
            mode_fail_count = mode_fail_count + 1
!$OMP END CRITICAL
            cycle
         end if

         ! Store eigenvalues and eigenvectors
         this%eigenvalues(:, ik) = eigenvals
         this%eigenvectors(:, :, ik) = h_k_copy
      end do
!$OMP END DO
      
      ! Cleanup thread-private arrays
      deallocate(h_k_copy, s_k_copy, eigenvals, work_complex, rwork)
      
!$OMP END PARALLEL

      if (mode_fail_count > 0 .and. use_generalized) then
         call g_logger%warning('diagonalize_hamiltonian: generalized_overlap had failures; consider ham_only for robustness.', __FILE__, __LINE__)
      end if
      if (this%gamma_bounds_diagnostics) call this%run_gamma_bounds_diagnostics()
      if (this%hall_diag_experimental) call this%diagonalize_hall_experimental()
      call root_info('diagonalize_hamiltonian: Completed successfully', __FILE__, __LINE__)
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

   !> @brief Check Hermiticity and basic diagnostics for an overlap matrix.
   !> @param[in] this Reciprocal object providing diagnostic context.
   !> @param[in] ik k-point index associated with s_k.
   !> @param[in] s_k Overlap matrix to inspect.
   module subroutine check_overlap_properties(this, ik, s_k)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: ik
      complex(rp), dimension(:, :), intent(in) :: s_k
      integer :: i, j, n
      real(rp) :: max_herm
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
      integer :: nsites, n_orb, n, i, jsite, isite, ineigh, ia, ja, nr, info, lwork
      integer :: i_start, i_end, j_start, j_end
      complex(rp), allocatable :: hall_mat(:, :), work(:)
      real(rp), allocatable :: evals(:), rwork(:)

      if (.not. this%hall_diag_experimental) return
      if (this%control%calctype /= 'I') then
         call g_logger%info('HALL experimental diagonalization skipped: calctype is not impurity.', __FILE__, __LINE__)
         return
      end if

      nsites = this%lattice%nmax
      n_orb = 18
      if (nsites <= 0) return
      n = nsites * n_orb
      allocate(hall_mat(n, n), evals(n), rwork(max(1, 3*n - 2)))
      hall_mat = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, nsites
         ia = isite
         nr = this%lattice%nn(ia, 1)
         i_start = (isite - 1) * n_orb + 1
         i_end = isite * n_orb
         do ineigh = 1, nr
            if (ineigh == 1) then
               jsite = isite
            else
               ja = this%lattice%nn(ia, ineigh)
               jsite = ja
               if (jsite < 1 .or. jsite > nsites) cycle
            end if
            j_start = (jsite - 1) * n_orb + 1
            j_end = jsite * n_orb
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
