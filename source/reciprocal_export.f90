submodule (reciprocal_mod) reciprocal_export
#ifdef USE_MPI
   use mpi
#endif
   implicit none

contains

   !> @brief Write a dense k-space eigensystem for Fermi-surface analysis.
   !> @details The binary payload is deliberately a simple stream rather than
   !>          a compiler-dependent unformatted sequential file.  Arrays are
   !>          written in Fortran column-major order as float64 values:
   !>          k_points, k_weights, eigenvalues, Re(eigenvectors),
   !>          Im(eigenvectors), and optionally projection_weights.  Optional
   !>          independently solved collinear spin sectors are written to a
   !>          separate '<base>.spin.bin' stream so old readers remain valid.
   !>          The text sidecar records shapes, conventions, and the basis map.
   module subroutine write_kspace_eigenpairs(this, scf_iteration)
      class(reciprocal), intent(inout) :: this
      integer, intent(in), optional :: scf_iteration

      integer :: nmat, nbands, nk_local, nk_global, nk_full, nsites, spin_nbands
      integer :: isite, atom_idx
      integer :: iteration, ios, unit_bin, unit_meta
#ifdef USE_MPI
      integer :: irank, offset, total_count
      integer, allocatable :: nk_counts(:), nk_displs(:)
      integer, allocatable :: value_counts(:), value_displs(:)
      integer, allocatable :: vector_counts(:), vector_displs(:)
#endif
      real(rp), allocatable :: eigenvalues_global(:, :), k_points_global(:, :), k_weights_global(:)
      complex(rp), allocatable :: eigenvectors_global(:, :, :)
      real(rp), allocatable :: projection_weights(:, :, :, :, :)
#ifdef USE_MPI
      real(rp) :: dummy_real(1)
      complex(rp) :: dummy_complex(1)
#endif
      real(rp) :: canonical_energy
      logical :: operator_changed, need_spectrum, distributed, spin_resolved_written
      real(rp) :: spin_h_offdiag_relative, spin_s_offdiag_relative
      character(len=512) :: base_name, binary_name, metadata_name

      if (.not. associated(this%hamiltonian) .or. .not. associated(this%lattice)) then
         call g_logger%warning('write_kspace_eigenpairs: reciprocal object is not initialized.', __FILE__, __LINE__)
         return
      end if

      nk_full = this%nk_mesh(1) * this%nk_mesh(2) * this%nk_mesh(3)
      if (nk_full <= 0) then
         call g_logger%warning('write_kspace_eigenpairs: invalid k-mesh dimensions; no export written.', __FILE__, __LINE__)
         return
      end if

      ! An FS needs the complete mesh.  This also makes every k-point
      ! eigenvector available for site/orbital/spin colouring; copying an
      ! irreducible representative would not preserve that information.
      need_spectrum = .false.
      call this%invalidate_if_operator_changed('reciprocal%write_kspace_eigenpairs', operator_changed)
      if (.not. allocated(this%k_points) .or. this%nk_total /= nk_full) then
         call root_info('write_kspace_eigenpairs: promoting the export to the full BZ mesh.', __FILE__, __LINE__)
         call this%generate_mp_mesh()
         need_spectrum = .true.
      end if
      if (operator_changed .or. .not. allocated(this%eigenvalues) .or. .not. allocated(this%eigenvectors)) then
         need_spectrum = .true.
      end if
      if (need_spectrum) then
         call this%build_kspace_hamiltonian()
         call this%diagonalize_hamiltonian()
      end if

      if (.not. allocated(this%k_points) .or. .not. allocated(this%k_weights) .or. &
          .not. allocated(this%eigenvalues) .or. .not. allocated(this%eigenvectors)) then
         call g_logger%warning('write_kspace_eigenpairs: final eigensystem is unavailable; no export written.', __FILE__, __LINE__)
         return
      end if

      nk_global = size(this%k_points, 2)
      nk_local = size(this%eigenvalues, 2)
      nmat = size(this%eigenvalues, 1)
      nbands = size(this%eigenvalues, 1)
      nsites = this%lattice%nrec
      spin_nbands = norb*nsites
      distributed = this%k_mesh_distributed_active
      if (nk_global /= nk_full .or. size(this%k_weights) /= nk_global .or. nmat <= 0 .or. nsites <= 0 .or. &
          any(shape(this%eigenvectors) /= [nmat, nbands, nk_local]) .or. &
          (nk_local <= 0 .and. .not. distributed)) then
         call g_logger%warning('write_kspace_eigenpairs: inconsistent final mesh/eigenpair dimensions; no export written.', &
                               __FILE__, __LINE__)
         return
      end if

      ! Resolve EF from the same final eigenvalues and normalized mesh weights
      ! that are stored in the payload.  This is independent of the DOS grid.
      canonical_energy = this%calculate_canonical_band_energy(this%auto_find_fermi)

      if (present(scf_iteration)) then
         iteration = scf_iteration
      else
         iteration = -1
      end if

      base_name = trim(this%eigenpair_output_file)
      if (len_trim(base_name) == 0) base_name = 'fermi_surface'
      spin_resolved_written = .false.
      spin_h_offdiag_relative = 0.0_rp
      spin_s_offdiag_relative = 0.0_rp
      if (this%write_spin_resolved_eigenpairs) then
         call write_spin_resolved_payload(this, base_name, spin_resolved_written, spin_h_offdiag_relative, &
                                          spin_s_offdiag_relative)
      end if

      if (rank == 0) then
         allocate(k_points_global(3, nk_global), k_weights_global(nk_global))
         k_points_global = this%k_points
         k_weights_global = this%k_weights
      end if

#ifdef USE_MPI
      if (distributed .and. numprocs > 1) then
         allocate(nk_counts(numprocs), nk_displs(numprocs))
         allocate(value_counts(numprocs), value_displs(numprocs))
         allocate(vector_counts(numprocs), vector_displs(numprocs))
         nk_counts = 0
         nk_displs = 0
         value_counts = 0
         value_displs = 0
         vector_counts = 0
         vector_displs = 0
         call MPI_GATHER(nk_local, 1, MPI_INTEGER, nk_counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if (rank == 0) then
            offset = 0
            do irank = 1, numprocs
               nk_displs(irank) = offset
               value_counts(irank) = nk_counts(irank) * nmat
               value_displs(irank) = offset * nmat
               vector_counts(irank) = nk_counts(irank) * nmat * nbands
               vector_displs(irank) = offset * nmat * nbands
               offset = offset + nk_counts(irank)
            end do
            total_count = offset
            if (total_count /= nk_global) then
               call g_logger%fatal('write_kspace_eigenpairs: MPI k-point ownership does not cover the full mesh.', &
                                   __FILE__, __LINE__)
            end if
            allocate(eigenvalues_global(nmat, nk_global), eigenvectors_global(nmat, nbands, nk_global))
         end if

         if (rank == 0) then
            call MPI_GATHERV(this%eigenvalues, nk_local*nmat, MPI_DOUBLE_PRECISION, eigenvalues_global, &
                             value_counts, value_displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_GATHERV(this%eigenvectors, nk_local*nmat*nbands, MPI_DOUBLE_COMPLEX, eigenvectors_global, &
                             vector_counts, vector_displs, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
         else
            call MPI_GATHERV(this%eigenvalues, nk_local*nmat, MPI_DOUBLE_PRECISION, dummy_real, &
                             value_counts, value_displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_GATHERV(this%eigenvectors, nk_local*nmat*nbands, MPI_DOUBLE_COMPLEX, dummy_complex, &
                             vector_counts, vector_displs, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
         end if
      else
#endif
         if (rank == 0) then
            allocate(eigenvalues_global(nmat, nk_global), eigenvectors_global(nmat, nbands, nk_global))
            eigenvalues_global = this%eigenvalues
            eigenvectors_global = this%eigenvectors
         end if
#ifdef USE_MPI
      end if
#endif

      if (rank /= 0) then
#ifdef USE_MPI
         if (allocated(nk_counts)) deallocate(nk_counts, nk_displs, value_counts, value_displs, vector_counts, vector_displs)
#endif
         return
      end if

      binary_name = trim(base_name) // '.bin'
      metadata_name = trim(base_name) // '.meta'

      open(newunit=unit_bin, file=trim(binary_name), status='replace', action='write', access='stream', &
           form='unformatted', iostat=ios)
      if (ios /= 0) then
         call g_logger%fatal('write_kspace_eigenpairs: could not open '//trim(binary_name), __FILE__, __LINE__)
      end if
      write(unit_bin, iostat=ios) k_points_global
      write(unit_bin, iostat=ios) k_weights_global
      write(unit_bin, iostat=ios) eigenvalues_global
      write(unit_bin, iostat=ios) real(eigenvectors_global, rp)
      write(unit_bin, iostat=ios) aimag(eigenvectors_global)

      if (this%write_eigenpair_projections) then
         if (.not. bdg_mode .and. nsites > 0) then
            if (mod(nmat, nsites) == 0 .and. mod(nmat/nsites, 2) == 0) then
               call make_export_projection_weights(this, eigenvectors_global, projection_weights)
               write(unit_bin, iostat=ios) projection_weights
            else
               call g_logger%warning('write_kspace_eigenpairs: projection request ignored for incompatible basis dimensions.', &
                                     __FILE__, __LINE__)
            end if
         else
            call g_logger%warning('write_kspace_eigenpairs: projection request ignored for BdG or empty-site basis.', &
                                  __FILE__, __LINE__)
         end if
      end if
      close(unit_bin)
      if (ios /= 0) then
         call g_logger%fatal('write_kspace_eigenpairs: error while writing '//trim(binary_name), __FILE__, __LINE__)
      end if

      open(newunit=unit_meta, file=trim(metadata_name), status='replace', action='write', iostat=ios)
      if (ios /= 0) then
         call g_logger%fatal('write_kspace_eigenpairs: could not open '//trim(metadata_name), __FILE__, __LINE__)
      end if
      write(unit_meta, '(A)') '# RSLMTO dense k-space Fermi-surface eigenpair export'
      write(unit_meta, '(A)') 'format = rslmto_kspace_eigenpairs_v1'
      write(unit_meta, '(A)') 'post_processing = fermi_surface'
      write(unit_meta, '(A,A)') 'binary_file = ', trim(binary_name)
      write(unit_meta, '(A)') 'binary_encoding = native-endian float64 stream'
      write(unit_meta, '(A)') 'array_order = Fortran_column_major'
      write(unit_meta, '(A)') 'energy_units = Ry'
      write(unit_meta, '(A)') 'k_coordinates = fractional_reciprocal_lattice_units'
      write(unit_meta, '(A,I0)') 'scf_iteration = ', iteration
      write(unit_meta, '(A,I0)') 'operator_generation = ', this%cached_operator_generation
      write(unit_meta, '(A,I0)') 'mpi_ranks = ', numprocs
      write(unit_meta, '(A,I0)') 'n_kpoints = ', nk_global
      write(unit_meta, '(A,I0)') 'n_bands = ', nbands
      write(unit_meta, '(A,I0)') 'n_basis = ', nmat
      write(unit_meta, '(A,I0)') 'n_sites = ', nsites
      do isite = 1, nsites
         atom_idx = this%lattice%nbulk + isite
         if (atom_idx >= 1 .and. atom_idx <= size(this%lattice%symbolic_atoms)) then
            write(unit_meta, '(A,I0,A,A)') 'site_', isite, '_species = ', &
               trim(this%lattice%symbolic_atoms(atom_idx)%element%symbol)
         else
            write(unit_meta, '(A,I0,A)') 'site_', isite, '_species = X'
         end if
      end do
      write(unit_meta, '(A)') 'basis_layout = site_major_spin_blocked'
      write(unit_meta, '(A,3(1X,I0))') 'k_mesh =', this%nk_mesh
      write(unit_meta, '(A,3(1X,ES24.16))') 'k_offset =', this%k_offset
      write(unit_meta, '(A,3(1X,ES24.16))') 'reciprocal_b1 =', this%reciprocal_vectors(:, 1)
      write(unit_meta, '(A,3(1X,ES24.16))') 'reciprocal_b2 =', this%reciprocal_vectors(:, 2)
      write(unit_meta, '(A,3(1X,ES24.16))') 'reciprocal_b3 =', this%reciprocal_vectors(:, 3)
      write(unit_meta, '(A,ES24.16)') 'fermi_level = ', this%fermi_level
      write(unit_meta, '(A,ES24.16)') 'canonical_band_energy = ', canonical_energy
      write(unit_meta, '(A,ES24.16)') 'k_weight_sum = ', sum(k_weights_global)
      write(unit_meta, '(A)') 'array_1 = k_points shape=(3,n_kpoints) dtype=float64'
      write(unit_meta, '(A)') 'array_2 = k_weights shape=(n_kpoints) dtype=float64'
      write(unit_meta, '(A)') 'array_3 = eigenvalues shape=(n_bands,n_kpoints) dtype=float64'
      write(unit_meta, '(A)') 'array_4 = eigenvectors_real shape=(n_basis,n_bands,n_kpoints) dtype=float64'
      write(unit_meta, '(A)') 'array_5 = eigenvectors_imag shape=(n_basis,n_bands,n_kpoints) dtype=float64'
      if (allocated(projection_weights)) then
         write(unit_meta, '(A)') 'projection_weights = enabled'
         write(unit_meta, '(A)') 'projection_shape = (n_sites,4,2,n_bands,n_kpoints)'
         write(unit_meta, '(A)') 'projection_axes = local_site_moment_axes'
         write(unit_meta, '(A)') 'projection_spin_1 = local_spin_up'
         write(unit_meta, '(A)') 'projection_spin_2 = local_spin_down'
         write(unit_meta, '(A)') 'projection_orbitals = l0_s,l1_p,l2_d,l3_f'
         write(unit_meta, '(A)') 'array_6 = projection_weights shape=(n_sites,4,2,n_bands,n_kpoints) dtype=float64'
      else
         write(unit_meta, '(A)') 'projection_weights = disabled'
      end if
      if (spin_resolved_written) then
         write(unit_meta, '(A)') 'spin_resolved_eigenpairs = enabled'
         write(unit_meta, '(A,A)') 'spin_resolved_binary_file = ', trim(base_name)//'.spin.bin'
         write(unit_meta, '(A)') 'spin_resolved_layout = global_spin_blocks'
         write(unit_meta, '(A)') 'spin_resolved_spin_1 = global_up'
         write(unit_meta, '(A)') 'spin_resolved_spin_2 = global_down'
         write(unit_meta, '(A,I0)') 'spin_resolved_n_spin = ', 2
         write(unit_meta, '(A,I0)') 'spin_resolved_n_bands = ', spin_nbands
         write(unit_meta, '(A,I0)') 'spin_resolved_n_basis = ', spin_nbands
         write(unit_meta, '(A,ES24.16)') 'spin_resolved_h_offdiag_relative_max = ', spin_h_offdiag_relative
         write(unit_meta, '(A,ES24.16)') 'spin_resolved_s_offdiag_relative_max = ', spin_s_offdiag_relative
         write(unit_meta, '(A,ES24.16)') 'spin_resolved_block_tolerance = ', 1.0e-10_rp
         write(unit_meta, '(A)') 'spin_array_1 = eigenvalues shape=(2,n_spin_bands,n_kpoints) dtype=float64'
         write(unit_meta, '(A)') 'spin_array_2 = eigenvectors_real shape=(2,n_spin_basis,n_spin_bands,n_kpoints) dtype=float64'
         write(unit_meta, '(A)') 'spin_array_3 = eigenvectors_imag shape=(2,n_spin_basis,n_spin_bands,n_kpoints) dtype=float64'
      else
         write(unit_meta, '(A)') 'spin_resolved_eigenpairs = disabled'
      end if
      write(unit_meta, '(A)') '# Basis index map: basis_index site orbital_l orbital_m_index spin_index'
      call write_export_basis_map(unit_meta, nmat, nsites)
      close(unit_meta)

      call root_info('write_kspace_eigenpairs: wrote '//trim(binary_name)//' and '//trim(metadata_name)// &
                     ' ('//trim(int2str(nk_global))//' k-points, '//trim(int2str(nmat))//' bands)', __FILE__, __LINE__)

      if (allocated(projection_weights)) deallocate(projection_weights)
      deallocate(k_points_global, k_weights_global, eigenvalues_global, eigenvectors_global)
#ifdef USE_MPI
      if (allocated(nk_counts)) deallocate(nk_counts, nk_displs, value_counts, value_displs, vector_counts, vector_displs)
#endif
   end subroutine write_kspace_eigenpairs

   !> Build and write independently diagonalized global spin sectors.
   !> @details The production basis is site-major, so the global up/down
   !> sectors are gathered from interleaved per-site indices.  This routine
   !> first checks the complete assembled H(k) (and S(k), when present) before
   !> solving either sector.  It deliberately does not classify the full
   !> eigensystem: degenerate up/down states can be arbitrarily mixed by a
   !> full-matrix eigensolver even when the operator is exactly block diagonal.
   subroutine write_spin_resolved_payload(this, base_name, written, max_h_relative, max_s_relative)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in) :: base_name
      logical, intent(out) :: written
      real(rp), intent(out) :: max_h_relative, max_s_relative

      real(rp), parameter :: spin_block_tolerance = 1.0e-10_rp
      integer :: nmat, nsites, nk_local, nk_global, spin_nbasis
      integer :: isite, iorb, i, ik, info, lwork
      integer, allocatable :: up_index(:), down_index(:)
      complex(rp), allocatable :: h_up(:, :), h_down(:, :), s_up(:, :), s_down(:, :)
      complex(rp), allocatable :: work(:)
      real(rp), allocatable :: rwork(:), spin_eigenvalues_local(:, :, :)
      complex(rp), allocatable :: spin_eigenvectors_local(:, :, :, :)
      real(rp), allocatable :: spin_eigenvalues_global(:, :, :)
      complex(rp), allocatable :: spin_eigenvectors_global(:, :, :, :)
      complex(rp) :: work_query(1)
      real(rp), allocatable :: eval_query(:)
      real(rp) :: local_h_relative, local_s_relative, matrix_scale, offdiag_norm
      logical :: use_generalized, distributed
      character(len=512) :: spin_binary_name
      integer :: spin_unit, ios
#ifdef USE_MPI
      integer :: irank, offset, total_count
      integer, allocatable :: nk_counts(:), eval_counts(:), eval_displs(:)
      integer, allocatable :: vector_counts(:), vector_displs(:)
      real(rp) :: dummy_real(1)
      complex(rp) :: dummy_complex(1)
#endif

      written = .false.
      max_h_relative = 0.0_rp
      max_s_relative = 0.0_rp

      if (this%has_nonzero_q_gbt()) then
         call g_logger%fatal('write_spin_resolved_payload: pure global spin sectors are undefined for a finite-q GBT spiral.', &
                             __FILE__, __LINE__)
      end if
      if (bdg_mode) then
         call g_logger%fatal('write_spin_resolved_payload: BdG/Nambu eigenpairs cannot be exported as ordinary up/down sectors.', &
                             __FILE__, __LINE__)
      end if
      if (.not. associated(this%control)) then
         call g_logger%fatal('write_spin_resolved_payload: no control state is available to verify collinearity.', &
                             __FILE__, __LINE__)
      end if
      if (this%control%nsp /= 1) then
         call g_logger%fatal('write_spin_resolved_payload: pure global spin sectors require nsp=1 '// &
                             '(collinear scalar-relativistic, no spin-orbit coupling).', __FILE__, __LINE__)
      end if
      if (.not. allocated(this%hk_bulk)) then
         call g_logger%fatal('write_spin_resolved_payload: assembled H(k) is unavailable.', __FILE__, __LINE__)
      end if

      nsites = this%lattice%nrec
      nmat = size(this%hk_bulk, 1)
      nk_local = size(this%hk_bulk, 3)
      nk_global = size(this%k_points, 2)
      spin_nbasis = norb*nsites
      if (nsites <= 0 .or. spin_off /= norb .or. nmat /= 2*spin_nbasis) then
         call g_logger%fatal('write_spin_resolved_payload: H(k) does not have the normal site-major spin basis.', &
                             __FILE__, __LINE__)
      end if
      if (size(this%hk_bulk, 2) /= nmat .or. nk_global <= 0) then
         call g_logger%fatal('write_spin_resolved_payload: invalid H(k) dimensions.', __FILE__, __LINE__)
      end if

      use_generalized = trim(this%reciprocal_mode) == 'generalized_overlap_proxy'
      if (use_generalized) then
         if (.not. allocated(this%sk_overlap)) then
            call g_logger%fatal('write_spin_resolved_payload: generalized mode requires a complete S(k) cache.', &
                                __FILE__, __LINE__)
         end if
         if (any(shape(this%sk_overlap) /= [nmat, nmat, nk_local])) then
            call g_logger%fatal('write_spin_resolved_payload: generalized S(k) cache dimensions are invalid.', &
                                __FILE__, __LINE__)
         end if
      end if

      allocate(up_index(spin_nbasis), down_index(spin_nbasis))
      i = 0
      do isite = 1, nsites
         do iorb = 1, norb
            i = i + 1
            up_index(i) = (isite - 1)*2*norb + iorb
            down_index(i) = (isite - 1)*2*norb + spin_off + iorb
         end do
      end do

      ! Check the already assembled matrices, including every k-point owned by
      ! this rank.  A relative criterion permits harmless round-off while
      ! rejecting actual transverse/SOC coupling.
      local_h_relative = 0.0_rp
      local_s_relative = 0.0_rp
      do ik = 1, nk_local
         matrix_scale = max(1.0_rp, maxval(abs(this%hk_bulk(:, :, ik))))
         offdiag_norm = spin_offdiag_norm(this%hk_bulk(:, :, ik), up_index, down_index)
         local_h_relative = max(local_h_relative, offdiag_norm/matrix_scale)
         if (use_generalized) then
            matrix_scale = max(1.0_rp, maxval(abs(this%sk_overlap(:, :, ik))))
            offdiag_norm = spin_offdiag_norm(this%sk_overlap(:, :, ik), up_index, down_index)
            local_s_relative = max(local_s_relative, offdiag_norm/matrix_scale)
         end if
      end do
#ifdef USE_MPI
      call MPI_ALLREDUCE(local_h_relative, max_h_relative, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(local_s_relative, max_s_relative, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
      max_h_relative = local_h_relative
      max_s_relative = local_s_relative
#endif
      if (max_h_relative > spin_block_tolerance .or. max_s_relative > spin_block_tolerance) then
         call g_logger%fatal('write_spin_resolved_payload: requested spin-resolved export but H/S has nonzero up/down coupling. '// &
                             'Use a collinear no-SOC calculation.', __FILE__, __LINE__)
      end if

      allocate(h_up(spin_nbasis, spin_nbasis), h_down(spin_nbasis, spin_nbasis), &
               eval_query(spin_nbasis), rwork(max(1, 3*spin_nbasis - 2)))
      h_up = cmplx(0.0_rp, 0.0_rp, rp)
      h_down = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, spin_nbasis
         h_up(i, i) = cmplx(0.0_rp, 0.0_rp, rp)
         h_down(i, i) = cmplx(0.0_rp, 0.0_rp, rp)
      end do
      if (use_generalized) then
         allocate(s_up(spin_nbasis, spin_nbasis), s_down(spin_nbasis, spin_nbasis))
         s_up = cmplx(0.0_rp, 0.0_rp, rp)
         s_down = cmplx(0.0_rp, 0.0_rp, rp)
         do i = 1, spin_nbasis
            s_up(i, i) = cmplx(1.0_rp, 0.0_rp, rp)
            s_down(i, i) = cmplx(1.0_rp, 0.0_rp, rp)
         end do
      end if

      if (use_generalized) then
         call zhegv(1, 'V', 'U', spin_nbasis, h_up, spin_nbasis, s_up, spin_nbasis, eval_query, &
                    work_query, -1, rwork, info)
      else
         call zheev('V', 'U', spin_nbasis, h_up, spin_nbasis, eval_query, work_query, -1, rwork, info)
      end if
      if (info /= 0) then
         call g_logger%fatal('write_spin_resolved_payload: LAPACK workspace query failed.', __FILE__, __LINE__)
      end if
      lwork = max(1, int(real(work_query(1), rp)))
      allocate(work(lwork))
      allocate(spin_eigenvalues_local(2, spin_nbasis, nk_local), &
               spin_eigenvectors_local(2, spin_nbasis, spin_nbasis, nk_local))

      do ik = 1, nk_local
         call extract_spin_block(this%hk_bulk(:, :, ik), up_index, h_up)
         call extract_spin_block(this%hk_bulk(:, :, ik), down_index, h_down)
         if (use_generalized) then
            call extract_spin_block(this%sk_overlap(:, :, ik), up_index, s_up)
            call extract_spin_block(this%sk_overlap(:, :, ik), down_index, s_down)
            call zhegv(1, 'V', 'U', spin_nbasis, h_up, spin_nbasis, s_up, spin_nbasis, &
                       spin_eigenvalues_local(1, :, ik), work, lwork, rwork, info)
            if (info /= 0) call g_logger%fatal('write_spin_resolved_payload: up-sector generalized solve failed.', &
                                               __FILE__, __LINE__)
            call zhegv(1, 'V', 'U', spin_nbasis, h_down, spin_nbasis, s_down, spin_nbasis, &
                       spin_eigenvalues_local(2, :, ik), work, lwork, rwork, info)
            if (info /= 0) call g_logger%fatal('write_spin_resolved_payload: down-sector generalized solve failed.', &
                                               __FILE__, __LINE__)
         else
            call zheev('V', 'U', spin_nbasis, h_up, spin_nbasis, spin_eigenvalues_local(1, :, ik), work, lwork, rwork, info)
            if (info /= 0) call g_logger%fatal('write_spin_resolved_payload: up-sector solve failed.', __FILE__, __LINE__)
            call zheev('V', 'U', spin_nbasis, h_down, spin_nbasis, spin_eigenvalues_local(2, :, ik), work, lwork, rwork, info)
            if (info /= 0) call g_logger%fatal('write_spin_resolved_payload: down-sector solve failed.', __FILE__, __LINE__)
         end if
         spin_eigenvectors_local(1, :, :, ik) = h_up
         spin_eigenvectors_local(2, :, :, ik) = h_down
      end do

      distributed = this%k_mesh_distributed_active
#ifdef USE_MPI
      if (distributed .and. numprocs > 1) then
         allocate(nk_counts(numprocs), eval_counts(numprocs), eval_displs(numprocs), &
                  vector_counts(numprocs), vector_displs(numprocs))
         call MPI_GATHER(nk_local, 1, MPI_INTEGER, nk_counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         if (rank == 0) then
            offset = 0
            do irank = 1, numprocs
               eval_counts(irank) = nk_counts(irank)*2*spin_nbasis
               eval_displs(irank) = offset*2*spin_nbasis
               vector_counts(irank) = nk_counts(irank)*2*spin_nbasis*spin_nbasis
               vector_displs(irank) = offset*2*spin_nbasis*spin_nbasis
               offset = offset + nk_counts(irank)
            end do
            total_count = offset
            if (total_count /= nk_global) then
               call g_logger%fatal('write_spin_resolved_payload: MPI k-point ownership does not cover the full mesh.', &
                                   __FILE__, __LINE__)
            end if
            allocate(spin_eigenvalues_global(2, spin_nbasis, nk_global), &
                     spin_eigenvectors_global(2, spin_nbasis, spin_nbasis, nk_global))
         end if
         if (rank == 0) then
            call MPI_GATHERV(spin_eigenvalues_local, nk_local*2*spin_nbasis, MPI_DOUBLE_PRECISION, &
                             spin_eigenvalues_global, eval_counts, eval_displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_GATHERV(spin_eigenvectors_local, nk_local*2*spin_nbasis*spin_nbasis, MPI_DOUBLE_COMPLEX, &
                             spin_eigenvectors_global, vector_counts, vector_displs, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
         else
            call MPI_GATHERV(spin_eigenvalues_local, nk_local*2*spin_nbasis, MPI_DOUBLE_PRECISION, dummy_real, &
                             eval_counts, eval_displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_GATHERV(spin_eigenvectors_local, nk_local*2*spin_nbasis*spin_nbasis, MPI_DOUBLE_COMPLEX, dummy_complex, &
                             vector_counts, vector_displs, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
         end if
      else
#endif
         if (rank == 0) then
            allocate(spin_eigenvalues_global(2, spin_nbasis, nk_global), &
                     spin_eigenvectors_global(2, spin_nbasis, spin_nbasis, nk_global))
            spin_eigenvalues_global = spin_eigenvalues_local
            spin_eigenvectors_global = spin_eigenvectors_local
         end if
#ifdef USE_MPI
      end if
#endif

      if (rank == 0) then
         spin_binary_name = trim(base_name)//'.spin.bin'
         open(newunit=spin_unit, file=trim(spin_binary_name), status='replace', action='write', access='stream', &
              form='unformatted', iostat=ios)
         if (ios /= 0) then
            call g_logger%fatal('write_spin_resolved_payload: could not open '//trim(spin_binary_name), __FILE__, __LINE__)
         end if
         write(spin_unit, iostat=ios) spin_eigenvalues_global
         write(spin_unit, iostat=ios) real(spin_eigenvectors_global, rp)
         write(spin_unit, iostat=ios) aimag(spin_eigenvectors_global)
         close(spin_unit)
         if (ios /= 0) then
            call g_logger%fatal('write_spin_resolved_payload: error while writing '//trim(spin_binary_name), __FILE__, __LINE__)
         end if
         call root_info('write_spin_resolved_payload: wrote '//trim(spin_binary_name)//' with independent global up/down sectors', &
                        __FILE__, __LINE__)
      end if
      written = .true.
#ifdef USE_MPI
      if (allocated(nk_counts)) deallocate(nk_counts, eval_counts, eval_displs, vector_counts, vector_displs)
#endif
   end subroutine write_spin_resolved_payload

   subroutine extract_spin_block(matrix, indices, block)
      complex(rp), intent(in) :: matrix(:, :)
      integer, intent(in) :: indices(:)
      complex(rp), intent(out) :: block(:, :)
      integer :: i, j

      if (any(shape(block) /= [size(indices), size(indices)])) then
         call g_logger%fatal('extract_spin_block: output shape does not match spin index list.', __FILE__, __LINE__)
      end if
      do j = 1, size(indices)
         do i = 1, size(indices)
            block(i, j) = matrix(indices(i), indices(j))
         end do
      end do
   end subroutine extract_spin_block

   real(rp) function spin_offdiag_norm(matrix, up_index, down_index) result(norm)
      complex(rp), intent(in) :: matrix(:, :)
      integer, intent(in) :: up_index(:), down_index(:)
      integer :: i, j

      norm = 0.0_rp
      do j = 1, size(down_index)
         do i = 1, size(up_index)
            norm = max(norm, abs(matrix(up_index(i), down_index(j))))
            norm = max(norm, abs(matrix(down_index(j), up_index(i))))
         end do
      end do
   end function spin_offdiag_norm

   !> Build site/orbital/local-spin weights in the same coefficient-space
   !> convention used by the reciprocal projected-DOS implementation.
   subroutine make_export_projection_weights(this, eigenvectors, projections)
      class(reciprocal), intent(in) :: this
      complex(rp), intent(in) :: eigenvectors(:, :, :)
      real(rp), allocatable, intent(out) :: projections(:, :, :, :, :)

      integer, parameter :: orbital_first(4) = [1, 2, 5, 10]
      integer, parameter :: orbital_last(4) = [1, 4, 9, 16]
      integer :: nmat, nbands, nk, nsites, site_block, norb_per_spin
      integer :: isite, iorb, ib, ik, io, site_offset, idx_up, idx_dn
      integer :: first_orb, last_orb
      integer :: atom_idx
      real(rp) :: total_weight, mx_weight, my_weight, mz_weight, local_weight
      real(rp) :: axis(3), axis_norm, up_weight, down_weight
      complex(rp) :: up_value, down_value, cross_value

      nmat = size(eigenvectors, 1)
      nbands = size(eigenvectors, 2)
      nk = size(eigenvectors, 3)
      nsites = this%lattice%nrec
      site_block = nmat/nsites
      norb_per_spin = site_block/2
      allocate(projections(nsites, 4, 2, nbands, nk))
      projections = 0.0_rp

      do isite = 1, nsites
         axis = [0.0_rp, 0.0_rp, 1.0_rp]
         atom_idx = this%lattice%nbulk + isite
         if (atom_idx >= 1 .and. atom_idx <= size(this%lattice%symbolic_atoms)) then
            if (allocated(this%lattice%symbolic_atoms(atom_idx)%potential%mom)) then
               axis = this%lattice%symbolic_atoms(atom_idx)%potential%mom(1:3)
            end if
         end if
         axis_norm = sqrt(sum(axis**2))
         if (axis_norm > tiny(1.0_rp)) axis = axis/axis_norm
         site_offset = (isite - 1)*site_block

         do ik = 1, nk
            do ib = 1, nbands
               do iorb = 1, 4
                  first_orb = orbital_first(iorb)
                  last_orb = min(orbital_last(iorb), norb_per_spin)
                  if (first_orb > last_orb) cycle
                  total_weight = 0.0_rp
                  mx_weight = 0.0_rp
                  my_weight = 0.0_rp
                  mz_weight = 0.0_rp
                  do io = first_orb, last_orb
                     idx_up = site_offset + io
                     idx_dn = site_offset + norb_per_spin + io
                     up_value = eigenvectors(idx_up, ib, ik)
                     down_value = eigenvectors(idx_dn, ib, ik)
                     up_weight = real(conjg(up_value)*up_value, rp)
                     down_weight = real(conjg(down_value)*down_value, rp)
                     cross_value = conjg(up_value)*down_value
                     total_weight = total_weight + up_weight + down_weight
                     mx_weight = mx_weight + 2.0_rp*real(cross_value, rp)
                     my_weight = my_weight + 2.0_rp*aimag(cross_value)
                     mz_weight = mz_weight + up_weight - down_weight
                  end do
                  local_weight = axis(1)*mx_weight + axis(2)*my_weight + axis(3)*mz_weight
                  projections(isite, iorb, 1, ib, ik) = 0.5_rp*(total_weight + local_weight)
                  projections(isite, iorb, 2, ib, ik) = 0.5_rp*(total_weight - local_weight)
               end do
            end do
         end do
      end do
   end subroutine make_export_projection_weights

   subroutine write_export_basis_map(unit, nmat, nsites)
      integer, intent(in) :: unit, nmat, nsites
      integer, parameter :: orbital_l(16) = [0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3]
      integer :: site_block, norb_per_spin, ibasis, site, local_index, spin_index, orbital_index
      integer :: orbital_m_index

      site_block = nmat/nsites
      norb_per_spin = site_block/2
      do ibasis = 1, nmat
         site = (ibasis - 1)/site_block + 1
         local_index = mod(ibasis - 1, site_block) + 1
         spin_index = 1
         orbital_index = local_index
         if (local_index > norb_per_spin) then
            spin_index = 2
            orbital_index = local_index - norb_per_spin
         end if
         orbital_m_index = orbital_index
         if (orbital_index <= size(orbital_l)) then
            write(unit, '(A,I0,1X,I0,1X,I0,1X,I0,1X,I0)') 'basis_index ', ibasis, site, &
               orbital_l(orbital_index), orbital_m_index, spin_index
         else
            write(unit, '(A,I0,1X,I0,1X,A,1X,I0,1X,I0)') 'basis_index ', ibasis, site, 'unknown', &
               orbital_m_index, spin_index
         end if
      end do
      write(unit, '(A,I0)') 'basis_site_block = ', site_block
      write(unit, '(A,I0)') 'basis_orbitals_per_spin = ', norb_per_spin
      write(unit, '(A,I0)') 'basis_lmax = ', lmax_basis
   end subroutine write_export_basis_map

end submodule reciprocal_export
