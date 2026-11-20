submodule (reciprocal_mod) reciprocal_dos
#ifdef USE_MPI
   use mpi
#endif
   implicit none

contains

   !> @brief Calculate reciprocal-space density of states.
   !> @details Prepares eigenvalues and the canonical Fermi level, then evaluates
   !>          total/projected DOS, band moments, and writes output.
   !> @param[inout] this Reciprocal object receiving DOS and moment arrays.
   !> @param[in] ham Hamiltonian object for k-space construction if needed.
   !> @param[in] n_energy_points Optional number of DOS grid points.
   !> @param[in] energy_range Optional energy window.
   !> @param[in] method Optional DOS method selector.
   !> @param[in] gaussian_sigma Optional Gaussian broadening.
   !> @param[in] temperature Optional Fermi-Dirac temperature.
   !> @param[in] fermi_level Optional fixed Fermi level.
   !> @param[in] total_electrons Optional electron count for automatic Fermi search.
   !> @param[in] auto_find_fermi Optional flag enabling Fermi-level search.
   !> @param[in] output_file Optional DOS output file.
   module subroutine calculate_density_of_states(this, ham, n_energy_points, energy_range, method, gaussian_sigma, temperature, fermi_level, total_electrons, auto_find_fermi, output_file)
      class(reciprocal), intent(inout) :: this
      class(hamiltonian), intent(in) :: ham
      integer, intent(in), optional :: n_energy_points
      real(rp), dimension(2), intent(in), optional :: energy_range
      character(len=*), intent(in), optional :: method
      real(rp), intent(in), optional :: gaussian_sigma
      real(rp), intent(in), optional :: temperature
      real(rp), intent(in), optional :: fermi_level
      real(rp), intent(in), optional :: total_electrons
      logical, intent(in), optional :: auto_find_fermi
      character(len=*), intent(in), optional :: output_file

      call this%prepare_density_of_states(ham, n_energy_points, energy_range, method, gaussian_sigma, &
         temperature, fermi_level, total_electrons, auto_find_fermi)
      call this%calculate_dos_output(output_file, .true.)
   end subroutine calculate_density_of_states

   !> @brief Diagonalize if needed and fix EF, electron count and canonical EBAND.
   !> @details First half of calculate_density_of_states; the canonical EBAND is
   !>          left in this%canonical_band_energy. Arguments as for calculate_density_of_states.
   module subroutine prepare_density_of_states(this, ham, n_energy_points, energy_range, method, gaussian_sigma, temperature, fermi_level, total_electrons, auto_find_fermi)
      class(reciprocal), intent(inout) :: this
      class(hamiltonian), intent(in) :: ham
      integer, intent(in), optional :: n_energy_points
      real(rp), dimension(2), intent(in), optional :: energy_range
      character(len=*), intent(in), optional :: method
      real(rp), intent(in), optional :: gaussian_sigma
      real(rp), intent(in), optional :: temperature
      real(rp), intent(in), optional :: fermi_level
      real(rp), intent(in), optional :: total_electrons
      logical, intent(in), optional :: auto_find_fermi

      real(rp) :: canonical_eband
      integer :: isite

      call root_info('calculate_density_of_states: Starting DOS calculation', __FILE__, __LINE__)

      ! DOS arrays and their eigensystem are operator-dependent. This catches
      ! q, cone/reference-frame, and potential rebuilds even when the reciprocal
      ! object and k mesh are deliberately reused across probes/SCF iterations.
      call this%invalidate_if_operator_changed('reciprocal%calculate_density_of_states')

      ! Set parameters from optional arguments
      if (present(n_energy_points)) this%n_energy_points = n_energy_points
      if (present(energy_range)) this%dos_energy_range = energy_range
      if (present(method)) this%dos_method = trim(method)
      if (present(gaussian_sigma)) this%gaussian_sigma = gaussian_sigma
      if (present(temperature)) this%temperature = temperature
      if (present(fermi_level)) this%fermi_level = fermi_level
      if (present(total_electrons)) this%total_electrons = total_electrons
      if (present(auto_find_fermi)) this%auto_find_fermi = auto_find_fermi

      ! Setup energy grid
      call this%setup_dos_energy_grid()

      ! Build k-space Hamiltonian and diagonalize if not already done
      if (.not. allocated(this%eigenvalues)) then
         call root_info('calculate_density_of_states: Building and diagonalizing Hamiltonian on k-mesh', __FILE__, __LINE__)
         if (.not. allocated(this%k_points)) then
            if (this%use_symmetry_reduction) then
               call this%generate_reduced_kpoint_mesh(this%nk_mesh, sum(abs(this%k_offset)) > 1.0e-12_rp)
            else
               call this%generate_mp_mesh()
            end if
         end if

         ! Build H(k) for all k-points in mesh
         call this%build_kspace_hamiltonian()
         
         ! Diagonalize to get eigenvalues
         call this%diagonalize_hamiltonian()
      end if

      ! EF, electron count, and the source-of-truth EBAND are properties of
      ! eigenvalues/weights/occupations, not of the chosen DOS output grid.
      canonical_eband = this%calculate_canonical_band_energy(this%auto_find_fermi)

      ! Spin axes at the time these eigenvectors are produced; the DOS output
      ! part may run after the potential moments have been updated.
      if (allocated(this%dos_site_axes)) deallocate(this%dos_site_axes)
      allocate(this%dos_site_axes(3, this%lattice%nrec))
      do isite = 1, this%lattice%nrec
         call get_site_spin_axis(this, isite, this%dos_site_axes(:, isite))
      end do
   end subroutine prepare_density_of_states

   !> @brief Total and projected DOS, energy cross-checks and DOS file.
   !> @details Second half of calculate_density_of_states; needs a prepared object.
   !> @param[in] output_file Optional DOS output file.
   !> @param[in] with_moments If true, band moments are evaluated here (standalone route).
   module subroutine calculate_dos_output(this, output_file, with_moments)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in), optional :: output_file
      logical, intent(in), optional :: with_moments

      character(len=100) :: filename
      character(len=360) :: energy_message
      real(rp) :: canonical_eband, dos_eband, projected_eband, dos_electrons, kT
      real(rp), allocatable :: saved_moments(:, :, :, :)
      real(rp) :: saved_electrons, saved_weight_sum
      logical :: saved_valid

      canonical_eband = this%canonical_band_energy
      ! Total and projected tetrahedron DOS integrate over the full uniform mesh; the SCF itself runs on the reduced one.
      if (trim(this%dos_method) == 'tetrahedron' .or. trim(this%dos_method) == 'blochl') then
         ! The rebuild drops the SCF band moments and canonical sums; they belong to the SCF mesh and are kept.
         call move_alloc(this%band_moments, saved_moments)
         saved_electrons = this%canonical_electron_count
         saved_weight_sum = this%canonical_weight_sum
         saved_valid = this%canonical_energy_valid
         call this%ensure_full_mesh_for_spinor_integrations('calculate_dos_output')
         if (allocated(saved_moments)) call move_alloc(saved_moments, this%band_moments)
         this%canonical_band_energy = canonical_eband
         this%canonical_electron_count = saved_electrons
         this%canonical_weight_sum = saved_weight_sum
         this%canonical_energy_valid = saved_valid
      end if
      filename = 'density_of_states.dat'
      if (present(output_file)) filename = output_file

      ! Calculate DOS based on method
      select case (trim(this%dos_method))
      case ('tetrahedron')
         call root_info('calculate_density_of_states: Using tetrahedron method', __FILE__, __LINE__)
         call this%calculate_dos_tetrahedron()
      case ('blochl')
         call root_info('calculate_density_of_states: Using Blöchl modified tetrahedron method', __FILE__, __LINE__)
         call this%calculate_dos_blochl()
      case ('gaussian')
         call root_info('calculate_density_of_states: Using Gaussian smearing method', __FILE__, __LINE__)
         call this%calculate_dos_gaussian()
      case default
         call g_logger%error('calculate_density_of_states: Unknown DOS method: ' // trim(this%dos_method), __FILE__, __LINE__)
         return
      end select

      call this%project_dos_orbitals()
      if (present(with_moments)) then
         if (with_moments) call this%calculate_band_moments()
      end if

      ! Grid/projection routes are diagnostics only.  They intentionally run
      ! after the canonical result and never feed EF or EBAND back into it.
      kT = max(this%temperature*kB_Ry_per_K, 1.0e-10_rp)
      dos_electrons = this%integrate_dos_up_to_energy(this%fermi_level, kT)
      dos_eband = this%calculate_band_energy_from_total_dos()
      projected_eband = this%calculate_band_energy_from_moments()
      write(energy_message, '(A,ES16.8,A,ES16.8,A,ES12.4,A,ES16.8,A,ES12.4,A,ES16.8,A,ES12.4)') &
         'k-space energy cross-checks: canonical EBAND=', canonical_eband, &
         ', total-DOS EBAND=', dos_eband, ', delta=', dos_eband - canonical_eband, &
         ', projected m0*m1=', projected_eband, ', delta=', projected_eband - canonical_eband, &
         ', DOS N=', dos_electrons, ', deltaN=', dos_electrons - this%canonical_electron_count
      call root_info(trim(energy_message), __FILE__, __LINE__)

      ! Write results to file
      call this%write_dos_to_file(filename)

      call root_info('calculate_density_of_states: DOS calculation completed', __FILE__, __LINE__)
   end subroutine calculate_dos_output

   !> @brief Validate full-to-irreducible k-point symmetry maps.
   !> @param[inout] this Reciprocal object containing symmetry maps and weights.
   !> @param[in] context_tag Diagnostic label for the caller/context.
   module subroutine validate_symmetry_kmap(this, context_tag)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in) :: context_tag
      integer :: nk_full_expected, nk_irred, i, idx
      integer, allocatable :: counts(:)
      real(rp) :: wsum

      if (.not. this%use_symmetry_reduction) return
      if (.not. allocated(this%full_to_irred_k)) return
      if (.not. allocated(this%k_weights)) return

      nk_full_expected = this%nk_mesh(1) * this%nk_mesh(2) * this%nk_mesh(3)
      nk_irred = this%nk_total
      if (nk_irred <= 0) return

      if (size(this%full_to_irred_k) /= nk_full_expected) then
         if (this%strict_symmetry_checks) then
            call g_logger%fatal(trim(context_tag) // ': full_to_irred_k size mismatch', __FILE__, __LINE__)
         else
            call g_logger%warning(trim(context_tag) // ': full_to_irred_k size mismatch', __FILE__, __LINE__)
            return
         end if
      end if

      if (any(this%full_to_irred_k < 1) .or. any(this%full_to_irred_k > nk_irred)) then
         if (this%strict_symmetry_checks) then
            call g_logger%fatal(trim(context_tag) // ': invalid full_to_irred_k entries', __FILE__, __LINE__)
         else
            call g_logger%warning(trim(context_tag) // ': invalid full_to_irred_k entries', __FILE__, __LINE__)
         end if
         return
      end if

      allocate(counts(nk_irred))
      counts = 0
      do i = 1, nk_full_expected
         idx = this%full_to_irred_k(i)
         counts(idx) = counts(idx) + 1
      end do

      if (size(this%k_weights) == nk_irred) then
         wsum = sum(this%k_weights)
         if (abs(wsum - 1.0_rp) > 1.0e-8_rp) then
            if (this%strict_symmetry_checks) then
               call g_logger%fatal(trim(context_tag) // ': k-point weights do not sum to 1', __FILE__, __LINE__)
            else
               call g_logger%warning(trim(context_tag) // ': k-point weights do not sum to 1', __FILE__, __LINE__)
            end if
         end if
      end if
      deallocate(counts)
   end subroutine validate_symmetry_kmap

   !> @brief Write symmetry k-point maps for diagnostics.
   !> @param[inout] this Reciprocal object containing symmetry maps.
   !> @param[in] filename Output file name.
   module subroutine write_symmetry_kmap_dump(this, filename)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in) :: filename
      integer :: u, i, nk_full_expected

      if (.not. allocated(this%full_to_irred_k)) return
      nk_full_expected = this%nk_mesh(1) * this%nk_mesh(2) * this%nk_mesh(3)
      if (size(this%full_to_irred_k) /= nk_full_expected) return

      open(newunit=u, file=trim(filename), status='replace', action='write')
      write(u,'(A)') '# full_k_index full_to_irred'
      do i = 1, nk_full_expected
         write(u,'(I12,1X,I12)') i, this%full_to_irred_k(i)
      end do
      close(u)
   end subroutine write_symmetry_kmap_dump

   !> @brief Ensure the configured tetrahedron symmetry backend is usable.
   !> @details Expands or preserves the k-mesh as needed for scalar tetrahedron
   !>          DOS versus spinor/projection integrations.
   !> @param[inout] this Reciprocal object whose mesh/backend state may be updated.
   module subroutine ensure_tetra_symmetry_backend(this)
      class(reciprocal), intent(inout) :: this
      integer :: nk_full_expected
      logical :: need_full

      nk_full_expected = this%nk_mesh(1) * this%nk_mesh(2) * this%nk_mesh(3)
      need_full = (this%nk_total /= nk_full_expected)
      if (.not. need_full) return

      call this%validate_symmetry_kmap('ensure_tetra_symmetry_backend')

      select case (trim(this%tetra_symmetry_mode))
      case ('irreducible_native')
         call g_logger%info('ensure_tetra_symmetry_backend: Using irreducible_native backend for scalar tetra DOS.', __FILE__, __LINE__)
         return
      case default
         continue
      end select

      ! Explicit compatibility path: rebuild full mesh and rediagonalize there.
      call g_logger%info('ensure_tetra_symmetry_backend: Switching to full mesh reference backend for tetra/blochl parity.', __FILE__, __LINE__)
      call this%generate_mp_mesh()
      call this%build_kspace_hamiltonian()
      call this%diagonalize_hamiltonian()
   end subroutine ensure_tetra_symmetry_backend

   !> @brief Ensure spinor integrations operate on a full k-mesh.
   !> @param[inout] this Reciprocal object whose eigenvalues/eigenvectors may be expanded.
   !> @param[in] context_tag Diagnostic label for the caller/context.
   module subroutine ensure_full_mesh_for_spinor_integrations(this, context_tag)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in) :: context_tag
      integer :: nk_full_expected

      nk_full_expected = this%nk_mesh(1) * this%nk_mesh(2) * this%nk_mesh(3)
      if (.not. this%use_symmetry_reduction) return
      if (this%nk_total >= nk_full_expected) return

      call g_logger%info(trim(context_tag) // ': symmetry-reduced spinor/projected integration requires full k mesh; rebuilding full mesh.', __FILE__, __LINE__)
      call this%generate_mp_mesh()
      call this%build_kspace_hamiltonian()
      call this%diagonalize_hamiltonian()
      if (allocated(this%tetrahedra)) deallocate(this%tetrahedra)
      if (allocated(this%tetrahedron_volumes)) deallocate(this%tetrahedron_volumes)
   end subroutine ensure_full_mesh_for_spinor_integrations

   !> @brief Build irreducible tetrahedra and multiplicities from full-mesh cuts.
   !> @param[inout] this Reciprocal object containing symmetry maps and mesh dimensions.
   !> @param[out] tet_ir Irreducible tetrahedron corner indices.
   !> @param[out] tet_mult Multiplicity of each irreducible tetrahedron.
   !> @param[out] n_tet_ir Number of irreducible tetrahedra.
   module subroutine build_irreducible_tetrahedra(this, tet_ir, tet_mult, n_tet_ir)
      class(reciprocal), intent(inout) :: this
      integer, allocatable, intent(out) :: tet_ir(:, :)
      integer, allocatable, intent(out) :: tet_mult(:)
      integer, intent(out) :: n_tet_ir
      integer :: nk1, nk2, nk3, i, j, k, it, ic, idx
      integer :: n_tet_per_cube, tet_full_count, key_pos
      integer :: nk_full_expected
      integer, dimension(3, 4, 6) :: tetra_cut
      integer, dimension(4) :: key_sorted
      integer, allocatable :: keys(:, :), mult_tmp(:)
      logical :: found

      n_tet_ir = 0
      nk_full_expected = this%nk_mesh(1) * this%nk_mesh(2) * this%nk_mesh(3)
      if (.not. allocated(this%full_to_irred_k)) then
         call g_logger%fatal('build_irreducible_tetrahedra: full_to_irred_k is not allocated', __FILE__, __LINE__)
      end if
      if (size(this%full_to_irred_k) /= nk_full_expected) then
         call g_logger%fatal('build_irreducible_tetrahedra: full_to_irred_k size mismatch', __FILE__, __LINE__)
      end if

      nk1 = this%nk_mesh(1)
      nk2 = this%nk_mesh(2)
      nk3 = this%nk_mesh(3)
      n_tet_per_cube = 6
      tet_full_count = n_tet_per_cube * nk1 * nk2 * nk3

      allocate(keys(4, tet_full_count))
      allocate(mult_tmp(tet_full_count))
      keys = 0
      mult_tmp = 0

      call get_tetra_cut_offsets(this, tetra_cut)

      do i = 1, nk1
         do j = 1, nk2
            do k = 1, nk3
               do it = 1, n_tet_per_cube
                  do ic = 1, 4
                     idx = this%get_kpoint_index(i + tetra_cut(1, ic, it), &
                                                 j + tetra_cut(2, ic, it), &
                                                 k + tetra_cut(3, ic, it), nk1, nk2, nk3)
                     key_sorted(ic) = this%full_to_irred_k(idx)
                  end do

                  call sort4_int(key_sorted)

                  found = .false.
                  do key_pos = 1, n_tet_ir
                     if (all(keys(:, key_pos) == key_sorted)) then
                        mult_tmp(key_pos) = mult_tmp(key_pos) + 1
                        found = .true.
                        exit
                     end if
                  end do
                  if (.not. found) then
                     n_tet_ir = n_tet_ir + 1
                     keys(:, n_tet_ir) = key_sorted
                     mult_tmp(n_tet_ir) = 1
                  end if
               end do
            end do
         end do
      end do

      allocate(tet_ir(4, n_tet_ir))
      allocate(tet_mult(n_tet_ir))
      tet_ir = keys(:, 1:n_tet_ir)
      tet_mult = mult_tmp(1:n_tet_ir)

      deallocate(keys, mult_tmp)
      call g_logger%info('build_irreducible_tetrahedra: Reduced ' // trim(int2str(tet_full_count)) // &
                         ' full tetrahedra to ' // trim(int2str(n_tet_ir)) // ' irreducible classes', __FILE__, __LINE__)
   contains
      subroutine sort4_int(arr)
         integer, intent(inout) :: arr(4)
         integer :: t
         if (arr(1) > arr(2)) then; t = arr(1); arr(1) = arr(2); arr(2) = t; end if
         if (arr(3) > arr(4)) then; t = arr(3); arr(3) = arr(4); arr(4) = t; end if
         if (arr(1) > arr(3)) then; t = arr(1); arr(1) = arr(3); arr(3) = t; end if
         if (arr(2) > arr(4)) then; t = arr(2); arr(2) = arr(4); arr(4) = t; end if
         if (arr(2) > arr(3)) then; t = arr(2); arr(2) = arr(3); arr(3) = t; end if
      end subroutine sort4_int
   end subroutine build_irreducible_tetrahedra

   !> @brief Select tetrahedron corner offsets for the current mesh geometry.
   !> @param[in] this Reciprocal object providing reciprocal lattice geometry.
   !> @param[out] tetra_cut Six tetrahedra per mesh cell, stored as corner offsets.
   module subroutine get_tetra_cut_offsets(this, tetra_cut)
      class(reciprocal), intent(in) :: this
      integer, intent(out) :: tetra_cut(3, 4, 6)

      integer :: base_cut(3, 4, 6), test_cut(3, 4, 6)
      integer :: lx, ly, it, ic, jc, best_lx, best_ly
      real(rp) :: qvec(3, 3), p(3, 4), dp(3)
      real(rp) :: edge2, max_edge2, best_max_edge2

      base_cut = 0
      base_cut(:, 1, 1) = [0, 0, 0]; base_cut(:, 2, 1) = [0, 1, 0]
      base_cut(:, 3, 1) = [1, 1, 0]; base_cut(:, 4, 1) = [1, 1, 1]
      base_cut(:, 1, 2) = [0, 0, 0]; base_cut(:, 2, 2) = [1, 0, 0]
      base_cut(:, 3, 2) = [1, 1, 0]; base_cut(:, 4, 2) = [1, 1, 1]
      base_cut(:, 1, 3) = [0, 0, 0]; base_cut(:, 2, 3) = [1, 0, 0]
      base_cut(:, 3, 3) = [1, 0, 1]; base_cut(:, 4, 3) = [1, 1, 1]
      base_cut(:, 1, 4) = [0, 0, 0]; base_cut(:, 2, 4) = [0, 1, 0]
      base_cut(:, 3, 4) = [0, 1, 1]; base_cut(:, 4, 4) = [1, 1, 1]
      base_cut(:, 1, 5) = [0, 0, 0]; base_cut(:, 2, 5) = [0, 0, 1]
      base_cut(:, 3, 5) = [0, 1, 1]; base_cut(:, 4, 5) = [1, 1, 1]
      base_cut(:, 1, 6) = [0, 0, 0]; base_cut(:, 2, 6) = [0, 0, 1]
      base_cut(:, 3, 6) = [1, 0, 1]; base_cut(:, 4, 6) = [1, 1, 1]

      qvec(:, 1) = this%reciprocal_vectors(:, 1) / real(max(1, this%nk_mesh(1)), rp)
      qvec(:, 2) = this%reciprocal_vectors(:, 2) / real(max(1, this%nk_mesh(2)), rp)
      qvec(:, 3) = this%reciprocal_vectors(:, 3) / real(max(1, this%nk_mesh(3)), rp)

      best_lx = 0
      best_ly = 0
      best_max_edge2 = huge(1.0_rp)
      do lx = 0, 1
         do ly = 0, 1
            call mirror_tetra_cut(base_cut, test_cut, lx, ly)
            max_edge2 = 0.0_rp
            do it = 1, 6
               do ic = 1, 4
                  p(:, ic) = qvec(:, 1) * real(test_cut(1, ic, it), rp) + &
                             qvec(:, 2) * real(test_cut(2, ic, it), rp) + &
                             qvec(:, 3) * real(test_cut(3, ic, it), rp)
               end do
               do ic = 1, 3
                  do jc = ic + 1, 4
                     dp = p(:, ic) - p(:, jc)
                     edge2 = sum(dp * dp)
                     max_edge2 = max(max_edge2, edge2)
                  end do
               end do
            end do
            if (max_edge2 < best_max_edge2) then
               best_max_edge2 = max_edge2
               best_lx = lx
               best_ly = ly
            end if
         end do
      end do

      call mirror_tetra_cut(base_cut, tetra_cut, best_lx, best_ly)
   contains
      subroutine mirror_tetra_cut(src, dst, lx_in, ly_in)
         integer, intent(in) :: src(3, 4, 6)
         integer, intent(out) :: dst(3, 4, 6)
         integer, intent(in) :: lx_in, ly_in

         dst = src
         if (lx_in == 1) dst(1, :, :) = 1 - dst(1, :, :)
         if (ly_in == 1) dst(2, :, :) = 1 - dst(2, :, :)
      end subroutine mirror_tetra_cut
   end subroutine get_tetra_cut_offsets

   !> @brief Build the DOS energy grid from the configured range and size.
   !> @param[inout] this Reciprocal object receiving dos_energy_grid.
   module subroutine setup_dos_energy_grid(this)
      class(reciprocal), intent(inout) :: this

      ! Local variables
      integer :: i
      real(rp) :: energy_min, energy_max, delta_energy

      energy_min = this%dos_energy_range(1)
      energy_max = this%dos_energy_range(2)
      delta_energy = (energy_max - energy_min) / real(this%n_energy_points - 1, rp)

      ! Allocate energy grid
      if (allocated(this%dos_energy_grid)) deallocate(this%dos_energy_grid)
      allocate(this%dos_energy_grid(this%n_energy_points))

      ! Fill energy grid in Ry (consistent with Hamiltonian and eigenvalues)
      do i = 1, this%n_energy_points
         this%dos_energy_grid(i) = energy_min + real(i-1, rp) * delta_energy
      end do

      call root_info('setup_dos_energy_grid: Created energy grid with ' // &
                     trim(int2str(this%n_energy_points)) // ' points from ' // &
                     trim(real2str(energy_min, '(F 8.5)')) // ' to ' // trim(real2str(energy_max, '(F 8.5)')) // ' Ry', &
                     __FILE__, __LINE__)
   end subroutine setup_dos_energy_grid

   !> @brief Calculate total DOS using the tetrahedron method.
   !> @details Accumulates DOS and integrated NOS over full or irreducible
   !>          tetrahedra, respecting MPI partitioning where active.
   !> @param[inout] this Reciprocal object receiving total_dos and total_nos.
   module subroutine calculate_dos_tetrahedron(this)
      class(reciprocal), intent(inout) :: this

      ! Local variables
      integer :: i_tet, i_corner, i_band, nbands
      real(rp), dimension(4) :: e_corners
      integer, allocatable :: tet_ir(:, :), tet_mult(:)
      integer :: n_tet_ir
      real(rp) :: tet_weight, dos_integral, nos_integral
      real(rp) :: fermi_count, fermi_error
      real(rp), allocatable :: local_dos(:), local_nos(:), eig_full(:, :)
      integer :: nk_full_expected, tet_start, tet_end, tet_count

      call g_logger%info('calculate_dos_tetrahedron: Calculating DOS using tetrahedron method', __FILE__, __LINE__)

      ! Setup tetrahedra if not already done
      if (.not. allocated(this%tetrahedra)) then
         call this%setup_tetrahedra()
      end if

      ! Allocate DOS arrays
      if (allocated(this%total_dos)) deallocate(this%total_dos)
      allocate(this%total_dos(this%n_energy_points))
      this%total_dos = 0.0_rp

      if (allocated(this%total_nos)) deallocate(this%total_nos)
      allocate(this%total_nos(this%n_energy_points))
      this%total_nos = 0.0_rp

      nbands = size(this%eigenvalues, 1)
      allocate(eig_full(nbands, this%k_workset%nk_global))
      call this%gather_k_columns(nbands, this%eigenvalues, eig_full)
      nk_full_expected = this%nk_mesh(1) * this%nk_mesh(2) * this%nk_mesh(3)

      if (this%use_symmetry_reduction .and. trim(this%tetra_symmetry_mode) == 'irreducible_native' .and. &
          this%nk_total < nk_full_expected) then
         call this%build_irreducible_tetrahedra(tet_ir, tet_mult, n_tet_ir)
         call g_logger%info('calculate_dos_tetrahedron: Using irreducible tetrahedra for scalar total DOS.', __FILE__, __LINE__)
      else
         n_tet_ir = this%n_tetrahedra
      end if
      call get_mpi_range(rank, n_tet_ir, tet_start, tet_end, tet_count, region_tag='tetra')

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP& SHARED(this, tet_ir, tet_mult, tet_start, tet_end, nbands, nk_full_expected, eig_full) &
!$OMP& PRIVATE(i_tet, i_band, i_corner, e_corners, tet_weight, local_dos, local_nos)
      allocate(local_dos(this%n_energy_points), local_nos(this%n_energy_points))
      local_dos = 0.0_rp
      local_nos = 0.0_rp

!$OMP DO SCHEDULE(DYNAMIC)
      do i_tet = tet_start, tet_end
         if (allocated(tet_mult)) then
            tet_weight = real(tet_mult(i_tet), rp) / (6.0_rp * real(nk_full_expected, rp))
         else
            tet_weight = this%tetrahedron_volumes(i_tet)
         end if
         do i_band = 1, nbands
            do i_corner = 1, 4
               if (allocated(tet_ir)) then
                  e_corners(i_corner) = eig_full(i_band, tet_ir(i_corner, i_tet))
               else
                  e_corners(i_corner) = eig_full(i_band, this%tetrahedra(i_corner, i_tet))
               end if
            end do
            call tetra_add_nos(tet_weight, e_corners, this%dos_energy_grid(1), &
                               this%dos_energy_grid(this%n_energy_points), local_nos, this%n_energy_points)
            call tetra_add_dos(tet_weight, e_corners, this%dos_energy_grid(1), &
                               this%dos_energy_grid(this%n_energy_points), local_dos, this%n_energy_points)
         end do
      end do
!$OMP END DO

!$OMP CRITICAL(tetra_dos_accum)
      this%total_dos = this%total_dos + local_dos
      this%total_nos = this%total_nos + local_nos
!$OMP END CRITICAL(tetra_dos_accum)

      deallocate(local_dos, local_nos)
!$OMP END PARALLEL

#ifdef USE_MPI
      if (numprocs > 1) then
         call MPI_ALLREDUCE(MPI_IN_PLACE, this%total_dos, this%n_energy_points, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, this%total_nos, this%n_energy_points, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      end if
#endif

      dos_integral = trapezoidal_integral(this%dos_energy_grid, this%total_dos)
      nos_integral = this%total_nos(this%n_energy_points) - this%total_nos(1)
      call g_logger%info('calculate_dos_tetrahedron: DOS integral = ' // &
                        trim(real2str(dos_integral, '(F12.6)')) // ', N(Emax)-N(Emin) = ' // &
                        trim(real2str(nos_integral, '(F12.6)')) // ', expected states = ' // &
                        trim(int2str(nbands)), __FILE__, __LINE__)

      if (this%total_electrons > 0.0_rp) then
         fermi_count = interpolate_grid_value(this%dos_energy_grid, this%total_nos, this%n_energy_points, this%fermi_level)
         fermi_error = fermi_count - this%total_electrons
         call g_logger%info('calculate_dos_tetrahedron: N(E_F=' // &
                           trim(real2str(this%fermi_level, '(F10.6)')) // ' Ry) = ' // &
                           trim(real2str(fermi_count, '(F12.6)')) // ', target valence = ' // &
                           trim(real2str(this%total_electrons, '(F12.6)')) // ', error = ' // &
                           trim(real2str(fermi_error, '(ES12.4)')), __FILE__, __LINE__)
      end if

      if (allocated(tet_ir)) deallocate(tet_ir)
      if (allocated(tet_mult)) deallocate(tet_mult)
      deallocate(eig_full)

      call g_logger%info('calculate_dos_tetrahedron: Tetrahedron DOS calculation completed', __FILE__, __LINE__)
   end subroutine calculate_dos_tetrahedron

   module subroutine gather_k_columns(this, ld, local, full)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: ld
      real(rp), intent(in) :: local(*)
      real(rp), intent(inout) :: full(*)
      integer :: irank, k_first, k_last, k_count
      integer, allocatable :: counts(:), displs(:)

      if (.not. this%k_mesh_distributed_active) then
         full(1:ld*this%nk_local) = local(1:ld*this%nk_local)
         return
      end if
#ifdef USE_MPI
      allocate(counts(0:numprocs - 1), displs(0:numprocs - 1))
      do irank = 0, numprocs - 1
         call get_mpi_range(irank, this%k_workset%nk_global, k_first, k_last, k_count)
         counts(irank) = ld*k_count
         displs(irank) = ld*(k_first - 1)
      end do
      call MPI_ALLGATHERV(local, counts(rank), MPI_DOUBLE_PRECISION, full, counts, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
      deallocate(counts, displs)
#endif
   end subroutine gather_k_columns

   !> @brief Sort four real values in ascending order.
   !> @param[in] arr_in Four input values.
   !> @param[out] arr_out Sorted output values.
   module subroutine sort4(arr_in, arr_out)
      real(rp), dimension(4), intent(in) :: arr_in
      real(rp), dimension(4), intent(out) :: arr_out
      real(rp) :: a1, a2, a3, a4, tmp

      a1 = arr_in(1)
      a2 = arr_in(2)
      a3 = arr_in(3)
      a4 = arr_in(4)

      ! compare-swap sequence
      if (a1 > a2) then
         tmp = a1
         a1 = a2
         a2 = tmp
      end if
      if (a3 > a4) then
         tmp = a3
         a3 = a4
         a4 = tmp
      end if
      if (a1 > a3) then
         tmp = a1
         a1 = a3
         a3 = tmp
      end if
      if (a2 > a4) then
         tmp = a2
         a2 = a4
         a4 = tmp
      end if
      if (a2 > a3) then
         tmp = a2
         a2 = a3
         a3 = tmp
      end if

      arr_out(1) = a1
      arr_out(2) = a2
      arr_out(3) = a3
      arr_out(4) = a4
   end subroutine sort4

   !> @brief Add one tetrahedron’s integrated-state contribution to a grid.
   !> @param[in] volwgt Tetrahedron volume/weight factor.
   !> @param[in] ecorn_in Tetrahedron corner energies.
   !> @param[in] emin Lower grid energy.
   !> @param[in] emax Upper grid energy.
   !> @param[inout] nos Integrated number-of-states grid.
   !> @param[in] npts Number of grid points.
   module subroutine tetra_add_nos(volwgt, ecorn_in, emin, emax, nos, npts)
      real(rp), intent(in) :: volwgt, ecorn_in(4), emin, emax
      integer, intent(in) :: npts
      real(rp), intent(inout) :: nos(npts)

      integer :: i, i0(4), m
      real(rp) :: ecorn(4), de, e1, e2, e3, e4, seg_tol
      real(rp) :: c0, c1, c2, c3, x, x0
      real(rp), parameter :: tol = 1.0e-12_rp

      call sort4(ecorn_in, ecorn)
      if (ecorn(1) > emax) return
      if (ecorn(4) < emin) then
         nos = nos + volwgt
         return
      end if

      e1 = ecorn(1); e2 = ecorn(2); e3 = ecorn(3); e4 = ecorn(4)
      de = (emax - emin) / real(npts - 1, rp)
      if (de <= 0.0_rp) return
      ! A segment whose own corner-to-corner width is far below the grid resolution
      ! is numerically degenerate at this resolution: the piecewise polynomial is only
      ! valid inside that (here vanishing) width, and its boundary-value limit as the
      ! width -> 0 is exactly zero (e.g. segment 1: (E-e1)^3/[(e2-e1)(e3-e1)(e4-e1)]
      ! at E=e2 reduces to (e2-e1)^2/[(e3-e1)(e4-e1)] -> 0). Grid-index rounding can
      ! still assign such a segment one non-empty bin whose grid point sits far outside
      ! the true (sub-bin) domain, which is what makes the ill-conditioned denominator
      ! dangerous; skipping the segment is the correct limit, not an approximation.
      seg_tol = max(1.0e-10_rp, 1.0e-6_rp*de)

      do m = 1, 4
         if (ecorn(m) <= emin) then
            i0(m) = 1
         else if (ecorn(m) > emax) then
            i0(m) = npts + 1
         else
            i0(m) = 2 + int((ecorn(m) - emin) / de - 1.0e-12_rp)
         end if
      end do

      if (i0(1) < i0(2) .and. (e2-e1) > seg_tol .and. abs((e2-e1)*(e3-e1)*(e4-e1)) > tol) then
         c3 = volwgt / ((e2-e1)*(e3-e1)*(e4-e1))
         x0 = emin - e1 - de
         do i = i0(1), i0(2) - 1
            x = x0 + real(i, rp) * de
            nos(i) = nos(i) + c3*x*x*x
         end do
      end if
      if (i0(2) < i0(3) .and. (e3-e2) > seg_tol .and. abs((e3-e1)*(e4-e1)*(e3-e2)*(e4-e2)) > tol) then
         c3 = volwgt*(e1+e2-e3-e4)/((e3-e1)*(e4-e1)*(e3-e2)*(e4-e2))
         c2 = volwgt*3.0_rp/((e3-e1)*(e4-e1))
         c1 = c2*(e2-e1)
         c0 = c1*(e2-e1)/3.0_rp
         x0 = emin - e2 - de
         do i = i0(2), i0(3) - 1
            x = x0 + real(i, rp) * de
            nos(i) = nos(i) + c0 + x*(c1 + x*(c2 + x*c3))
         end do
      end if
      if (i0(3) < i0(4) .and. (e4-e3) > seg_tol .and. abs((e4-e3)*(e4-e2)*(e4-e1)) > tol) then
         c3 = volwgt/((e4-e3)*(e4-e2)*(e4-e1))
         x0 = emin - e4 - de
         do i = i0(3), i0(4) - 1
            x = x0 + real(i, rp) * de
            nos(i) = nos(i) + volwgt + c3*x*x*x
         end do
      end if
      do i = max(1, i0(4)), npts
         nos(i) = nos(i) + volwgt
      end do
   end subroutine tetra_add_nos

   !> @brief Add one tetrahedron’s DOS contribution to a grid.
   !> @param[in] volwgt Tetrahedron volume/weight factor.
   !> @param[in] ecorn_in Tetrahedron corner energies.
   !> @param[in] emin Lower grid energy.
   !> @param[in] emax Upper grid energy.
   !> @param[inout] dos DOS grid.
   !> @param[in] npts Number of grid points.
   module subroutine tetra_add_dos(volwgt, ecorn_in, emin, emax, dos, npts)
      real(rp), intent(in) :: volwgt, ecorn_in(4), emin, emax
      integer, intent(in) :: npts
      real(rp), intent(inout) :: dos(npts)

      integer :: i, i0(4), m
      real(rp) :: ecorn(4), de, e1, e2, e3, e4, seg_tol
      real(rp) :: c1, c2, c3, x
      real(rp), parameter :: tol = 1.0e-12_rp

      call sort4(ecorn_in, ecorn)
      e1 = ecorn(1); e2 = ecorn(2); e3 = ecorn(3); e4 = ecorn(4)
      if (e4 < emin .or. e1 > emax) return

      de = (emax - emin) / real(npts - 1, rp)
      if (de <= 0.0_rp) return
      ! A band that is constant over a tetrahedron has a delta-function DOS,
      ! while its cumulative contribution is a unit step.  Preserve that
      ! state on the sampled DOS grid instead of dropping it when every
      ! piecewise denominator is degenerate.
      if (e4 - e1 <= max(1.0e-10_rp, 1.0e-6_rp*de)) then
         i = max(1, min(npts, 2 + int((0.5_rp*(e1 + e4) - emin) / de - 1.0e-12_rp)))
         if (i == 1 .or. i == npts) then
            dos(i) = dos(i) + 2.0_rp*volwgt/de
         else
            dos(i) = dos(i) + volwgt/de
         end if
         return
      end if
      ! See tetra_add_nos for the justification: a segment narrower than the grid
      ! resolution is treated as degenerate (its own boundary-value limit is zero),
      ! rather than letting grid-index rounding evaluate the polynomial out of range.
      seg_tol = max(1.0e-10_rp, 1.0e-6_rp*de)

      do m = 1, 4
         if (ecorn(m) < emin) then
            i0(m) = 1
         else if (ecorn(m) > emax) then
            i0(m) = npts + 1
         else
            i0(m) = 2 + int((ecorn(m) - emin) / de - 1.0e-12_rp)
         end if
      end do

      if (i0(1) < i0(2) .and. (e2-e1) > seg_tol .and. abs((e2-e1)*(e3-e1)*(e4-e1)) > tol) then
         c3 = volwgt*3.0_rp/((e2-e1)*(e3-e1)*(e4-e1))
         do i = i0(1), i0(2) - 1
            x = emin - e1 + real(i - 1, rp) * de
            dos(i) = dos(i) + c3*x*x
         end do
      end if
      if (i0(2) < i0(3) .and. (e3-e2) > seg_tol .and. abs((e3-e1)*(e4-e1)*(e3-e2)*(e4-e2)) > tol) then
         c3 = volwgt*3.0_rp*(e1+e2-e3-e4)/((e3-e1)*(e4-e1)*(e3-e2)*(e4-e2))
         c2 = volwgt*6.0_rp/((e3-e1)*(e4-e1))
         c1 = c2*(e2-e1)*0.5_rp
         do i = i0(2), i0(3) - 1
            x = emin - e2 + real(i - 1, rp) * de
            dos(i) = dos(i) + c1 + x*(c2 + x*c3)
         end do
      end if
      if (i0(3) < i0(4) .and. (e4-e3) > seg_tol .and. abs((e4-e3)*(e4-e2)*(e4-e1)) > tol) then
         c3 = volwgt*3.0_rp/((e4-e3)*(e4-e2)*(e4-e1))
         do i = i0(3), i0(4) - 1
            x = emin - e4 + real(i - 1, rp) * de
            dos(i) = dos(i) + c3*x*x
         end do
      end if
   end subroutine tetra_add_dos

!> @brief Calculate total DOS with Gaussian broadening.
!> @details Smears all eigenvalues over the DOS energy grid using the configured
!>          or adaptive sigma and k-point weights.
!> @param[inout] this Reciprocal object receiving total_dos.
module subroutine calculate_dos_gaussian(this)
   class(reciprocal), intent(inout) :: this

   integer :: i_energy, i_k, i_band, i_k_global
   real(rp) :: energy, weight, gaussian_factor
   real(rp) :: sigma_squared, sigma_use
   real(rp) :: local_sum, dos_integral, norm_factor, kweight_sum
   integer :: nbands

   ! Debug: Check eigenvalue range
   if (.not. this%suppress_internal_logs) then
      call g_logger%info('calculate_dos_gaussian: Eigenvalue range = [' // &
         trim(real2str(minval(this%eigenvalues), '(F12.6)')) // ', ' // &
         trim(real2str(maxval(this%eigenvalues), '(F12.6)')) // '] Ry', __FILE__, __LINE__)
      call g_logger%info('calculate_dos_gaussian: DOS energy range = [' // &
         trim(real2str(this%dos_energy_range(1), '(F12.6)')) // ', ' // &
         trim(real2str(this%dos_energy_range(2), '(F12.6)')) // '] Ry', __FILE__, __LINE__)
      call g_logger%info('calculate_dos_gaussian: Number of k-points = ' // &
         trim(int2str(max(this%nk_local, size(this%eigenvalues, 2)))), __FILE__, __LINE__)
      call g_logger%info('calculate_dos_gaussian: K-point weight = ' // &
         trim(real2str(this%k_weights(1), '(ES12.4)')), __FILE__, __LINE__)
   end if

   ! Determine sigma (already in Ry from input)
   if (this%gaussian_sigma < 0.001_rp) then
      sigma_use = this%calculate_adaptive_sigma()
      call root_info('calculate_dos_gaussian: Using adaptive sigma = ' // &
                        trim(real2str(sigma_use, '(F8.5)')) // ' Ry', __FILE__, __LINE__)
   else
      sigma_use = this%gaussian_sigma
      call root_info('calculate_dos_gaussian: Using input sigma = ' // &
                        trim(real2str(sigma_use, '(F8.5)')) // ' Ry', __FILE__, __LINE__)
   end if

   ! Allocate DOS arrays
   if (allocated(this%total_dos)) deallocate(this%total_dos)
   allocate(this%total_dos(this%n_energy_points))
   this%total_dos = 0.0_rp

   sigma_squared = sigma_use**2
   nbands = size(this%eigenvalues, 1)
   
   ! DEBUG: Check k-point weights sum
   if (this%k_mesh_distributed_active) then
      kweight_sum = sum(this%k_weights(this%k_start:this%k_end))
   else
      kweight_sum = sum(this%k_weights)
   end if
   
   call root_info('calculate_dos_gaussian: nbands = ' // trim(int2str(nbands)) // &
                     ', nk_local = ' // trim(int2str(size(this%eigenvalues, 2))) // &
                     ', k_weights sum = ' // trim(real2str(kweight_sum, '(F12.8)')), __FILE__, __LINE__)
   
   ! DEBUG: Check eigenvalue array size
   call root_info('calculate_dos_gaussian: eigenvalues array size = ' // &
                     trim(int2str(size(this%eigenvalues, 1))) // ' x ' // &
                     trim(int2str(size(this%eigenvalues, 2))), __FILE__, __LINE__)

   ! Calculate raw DOS
   do i_energy = 1, this%n_energy_points
      energy = this%dos_energy_grid(i_energy)  ! Already in Ry
      local_sum = 0.0_rp

!$OMP PARALLEL DO PRIVATE(i_k, i_band, weight, gaussian_factor) REDUCTION(+:local_sum) &
!$OMP& SCHEDULE(STATIC) IF(size(this%eigenvalues, 2) > 100)
      do i_k = 1, size(this%eigenvalues, 2)
         i_k_global = local_k_index_to_global(this, i_k)
         weight = this%k_weights(i_k_global)
         do i_band = 1, nbands
            gaussian_factor = exp(-((energy - this%eigenvalues(i_band, i_k))**2) / (2.0_rp * sigma_squared))
            gaussian_factor = gaussian_factor / (sigma_use * sqrt(2.0_rp * 3.141592653589793_rp))
            local_sum = local_sum + weight * gaussian_factor
         end do
      end do
!$OMP END PARALLEL DO
      this%total_dos(i_energy) = local_sum
   end do

#ifdef USE_MPI
   if (this%k_mesh_distributed_active) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%total_dos, this%n_energy_points, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
   end if
#endif

   ! DEBUG: Check DOS integral WITHOUT normalization first
   dos_integral = 0.0_rp
   do i_energy = 1, this%n_energy_points - 1
      dos_integral = dos_integral + 0.5_rp * (this%total_dos(i_energy) + this%total_dos(i_energy+1)) * &
                                   (this%dos_energy_grid(i_energy+1) - this%dos_energy_grid(i_energy))
   end do
   
   call root_info('calculate_dos_gaussian: Raw DOS (before norm) integrates to ' // &
                     trim(real2str(dos_integral, '(F12.6)')) // ' (should be ' // trim(int2str(nbands)) // ')', &
                     __FILE__, __LINE__)
   
   ! Normalize DOS to integrate to nbands (total number of states)
   if (abs(dos_integral) > 1.0e-10_rp) then
      norm_factor = 1.0_rp ! real(nbands, rp) / dos_integral
      this%total_dos = this%total_dos * norm_factor
      
      call root_info('calculate_dos_gaussian: DOS normalized by factor ' // &
                        trim(real2str(norm_factor, '(F10.6)')), __FILE__, __LINE__)
      call root_info('calculate_dos_gaussian: DOS integrates to ' // &
                        trim(real2str(dos_integral * norm_factor, '(F10.4)')) // &
                        ' (should be ' // trim(int2str(nbands)) // ')', __FILE__, __LINE__)
   else
      call g_logger%error('calculate_dos_gaussian: DOS integral is zero!', __FILE__, __LINE__)
   end if

   call root_info('calculate_dos_gaussian: Gaussian DOS calculation completed', __FILE__, __LINE__)
end subroutine calculate_dos_gaussian

   !> @brief Calculate total DOS with the Bloechl tetrahedron method.
   !> @param[inout] this Reciprocal object receiving total_dos.
   module subroutine calculate_dos_blochl(this)
      class(reciprocal), intent(inout) :: this

      call g_logger%info('calculate_dos_blochl: Using LMTO-style tetrahedron DOS/NOS backend for scalar DOS.', __FILE__, __LINE__)
      call this%calculate_dos_tetrahedron()
      call g_logger%info('calculate_dos_blochl: Blöchl-compatible scalar DOS calculation completed', __FILE__, __LINE__)
   end subroutine calculate_dos_blochl

   !> @brief Build full-mesh tetrahedron connectivity.
   !> @param[inout] this Reciprocal object receiving tetrahedra and volumes.
   module subroutine setup_tetrahedra(this)
      class(reciprocal), intent(inout) :: this

      ! Local variables
      integer :: nk1, nk2, nk3, n_tet_per_cube, i, j, k, tet_idx
      integer :: it, ic, idx
      integer, dimension(3, 4, 6) :: tetra_cut

      call g_logger%info('setup_tetrahedra: Setting up tetrahedra for Brillouin zone integration', __FILE__, __LINE__)

      ! Get mesh dimensions
      nk1 = this%nk_mesh(1)
      nk2 = this%nk_mesh(2)
      nk3 = this%nk_mesh(3)

      ! Number of tetrahedra per cube (6 for standard decomposition)
      n_tet_per_cube = 6

      ! Total number of tetrahedra
      this%n_tetrahedra = n_tet_per_cube * nk1 * nk2 * nk3

      ! Allocate tetrahedra array
      if (allocated(this%tetrahedra)) deallocate(this%tetrahedra)
      allocate(this%tetrahedra(4, this%n_tetrahedra))

      ! Allocate tetrahedron volumes (all equal for uniform mesh)
      if (allocated(this%tetrahedron_volumes)) deallocate(this%tetrahedron_volumes)
      allocate(this%tetrahedron_volumes(this%n_tetrahedra))

      ! Volume of each tetrahedron (1/6 of cube volume, times number of cubes)
      this%tetrahedron_volumes = 1.0_rp / (6.0_rp * real(nk1 * nk2 * nk3, rp))

      call get_tetra_cut_offsets(this, tetra_cut)

      ! Build tetrahedra
      tet_idx = 0
      do i = 1, nk1
         do j = 1, nk2
            do k = 1, nk3
               ! For each cube in the mesh
               do it = 1, n_tet_per_cube
                  tet_idx = tet_idx + 1

                  ! Convert relative corner indices to absolute k-point indices
                  do ic = 1, 4
                     idx = this%get_kpoint_index(i + tetra_cut(1, ic, it), &
                                                 j + tetra_cut(2, ic, it), &
                                                 k + tetra_cut(3, ic, it), nk1, nk2, nk3)
                     this%tetrahedra(ic, tet_idx) = idx
                  end do
               end do
            end do
         end do
      end do

      call g_logger%info('setup_tetrahedra: Created ' // trim(int2str(this%n_tetrahedra)) // &
                        ' tetrahedra from ' // trim(int2str(nk1)) // 'x' // trim(int2str(nk2)) // &
                        'x' // trim(int2str(nk3)) // ' k-mesh', __FILE__, __LINE__)
   end subroutine setup_tetrahedra

   !> @brief Convert 3D mesh coordinates to a periodic flat k-point index.
   !> @param[in] this Reciprocal object providing helper context.
   !> @param[in] i Mesh coordinate along axis 1.
   !> @param[in] j Mesh coordinate along axis 2.
   !> @param[in] k Mesh coordinate along axis 3.
   !> @param[in] nk1 Mesh size along axis 1.
   !> @param[in] nk2 Mesh size along axis 2.
   !> @param[in] nk3 Mesh size along axis 3.
   !> @return One-based periodic flat k-point index.
   module function get_kpoint_index(this, i, j, k, nk1, nk2, nk3) result(idx)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: i, j, k, nk1, nk2, nk3
      integer :: idx, ii, jj, kk

      ! Apply periodic boundary conditions
      ii = mod(i-1, nk1) + 1
      jj = mod(j-1, nk2) + 1
      kk = mod(k-1, nk3) + 1

      ! Convert to 1D index (assuming k-points are stored as k1 varying fastest)
      idx = ii + (jj-1)*nk1 + (kk-1)*nk1*nk2
   end function get_kpoint_index

   !> @brief Estimate a Gaussian DOS broadening from k-point density.
   !> @param[in] this Reciprocal object containing mesh and reciprocal volume.
   !> @return Adaptive Gaussian sigma.
   module function calculate_adaptive_sigma(this) result(sigma)
      class(reciprocal), intent(in) :: this
      real(rp) :: sigma
      real(rp) :: bz_volume, k_density, typical_spacing
      integer :: nk_total
      
      ! Calculate Brillouin zone k-point density
      nk_total = this%nk_mesh(1) * this%nk_mesh(2) * this%nk_mesh(3)
      
      if (nk_total > 0) then
         ! BZ volume per k-point (in reciprocal space units)
         bz_volume = this%reciprocal_volume
         k_density = real(nk_total, rp) / bz_volume
         
         ! Typical energy spacing scales as 1/k_density^(1/3)
         ! For metallic systems: ΔE ∝ E_F / N_k^(1/3)
         ! Use heuristic: sigma = C / nk^(1/3) where C is tuned for accuracy
      ! Use a smaller scale factor (1.0 eV) for typical spacing to avoid huge sigma
      typical_spacing = 1.0_rp / (real(nk_total, rp)**(1.0_rp/3.0_rp))
         
         ! Adaptive sigma: smaller for denser meshes, larger for coarse meshes
         ! Factor of 0.5-1.0 gives good balance between accuracy and smoothness
         sigma = 0.7_rp * typical_spacing
         
         ! Clamp to reasonable range
         ! Clamp sigma to a conservative range (5 meV .. 50 meV) to avoid over-broadening
         sigma = max(0.005_rp, min(sigma, 0.05_rp))
      else
         sigma = 0.1_rp  ! Default fallback
      end if
      
   call g_logger%info('calculate_adaptive_sigma: Adaptive sigma = ' // trim(real2str(sigma, '(F8.5)')) // &
            ' Ry for ' // trim(int2str(nk_total)) // ' k-points', __FILE__, __LINE__)
   end function calculate_adaptive_sigma

!> @brief Integrate DOS occupation up to an energy with Fermi weighting.
!> @param[in] this Reciprocal object containing DOS grid data.
!> @param[in] energy Trial Fermi energy.
!> @param[in] kT Thermal energy used in Fermi weighting.
!> @return Integrated electron count.
module function integrate_dos_up_to_energy(this, energy, kT) result(integral)
   class(reciprocal), intent(in) :: this
   real(rp), intent(in) :: energy, kT
   real(rp) :: integral

   integer :: ie
   real(rp) :: e, fermi_weight, fermi_weight_next, delta_e, x, x_next
   real(rp), parameter :: exp_limit = 700.0_rp

   integral = 0.0_rp

   do ie = 1, this%n_energy_points - 1
      ! dos_energy_grid is already in Ry, no conversion needed
      e = this%dos_energy_grid(ie)
      delta_e = this%dos_energy_grid(ie+1) - this%dos_energy_grid(ie)

      ! Fermi-Dirac weight at current energy
      if (kT > 1.0e-10_rp) then
         x = (e - energy) / kT
         x_next = (this%dos_energy_grid(ie+1) - energy) / kT
         if (x >= exp_limit) then
            fermi_weight = 0.0_rp
         else if (x <= -exp_limit) then
            fermi_weight = 1.0_rp
         else
            fermi_weight = 1.0_rp / (exp(x) + 1.0_rp)
         end if
         if (x_next >= exp_limit) then
            fermi_weight_next = 0.0_rp
         else if (x_next <= -exp_limit) then
            fermi_weight_next = 1.0_rp
         else
            fermi_weight_next = 1.0_rp / (exp(x_next) + 1.0_rp)
         end if
      else
         ! T=0 limit
         if (e <= energy) then
            fermi_weight = 1.0_rp
         else
            fermi_weight = 0.0_rp
         end if
         if (this%dos_energy_grid(ie+1) <= energy) then
            fermi_weight_next = 1.0_rp
         else
            fermi_weight_next = 0.0_rp
         end if
      end if

      ! Trapezoidal integration
      integral = integral + 0.5_rp * delta_e * (this%total_dos(ie) * fermi_weight + &
                                               this%total_dos(ie+1) * fermi_weight_next)
   end do
end function integrate_dos_up_to_energy

   !> @brief Write total and projected DOS data to output files.
   !> @param[in] this Reciprocal object containing DOS arrays.
   !> @param[in] filename Primary total-DOS output file.
   module subroutine write_dos_to_file(this, filename)
      class(reciprocal), intent(in) :: this
      character(len=*), intent(in) :: filename

      ! Local variables
      integer :: unit, i_energy, isite, iorb, ispin
      character(len=256) :: proj_filename

      if (rank /= 0) return

      call root_info('write_dos_to_file: Writing DOS to ' // trim(filename), __FILE__, __LINE__)

      ! Write total DOS
      open(newunit=unit, file=trim(filename), status='replace', action='write')

      ! Write header
      write(unit, '(A)') '# Density of States'
      write(unit, '(A,A)') '# Method: ', trim(this%dos_method)
      if (trim(this%dos_method) == 'gaussian') then
         write(unit, '(A,F8.5,A)') '# Gaussian sigma: ', this%gaussian_sigma, ' Ry'
      end if
      write(unit, '(A,I0)') '# Energy points: ', this%n_energy_points
      write(unit, '(A,2F10.6,A)') '# Absolute energy range: ', this%dos_energy_range(1), &
                                   this%dos_energy_range(2), ' Ry'
      write(unit, '(A,F12.6,A)') '# Fermi level: ', this%fermi_level, ' Ry'
      write(unit, '(A,2F10.6,A)') '# Printed energy range E-E_F: ', &
                                   this%dos_energy_range(1) - this%fermi_level, &
                                   this%dos_energy_range(2) - this%fermi_level, ' Ry'
      write(unit, '(A)') '# Energy-E_F (Ry)    Total DOS'

      ! Write DOS data (energy grid already in Ry)
      do i_energy = 1, this%n_energy_points
         write(unit, '(2F12.6)') this%dos_energy_grid(i_energy) - this%fermi_level, this%total_dos(i_energy)
      end do

      close(unit)

      ! Write projected DOS if available
      if (allocated(this%projected_dos)) then
         proj_filename = 'projected_dos.dat'
         call root_info('write_dos_to_file: Writing projected DOS to ' // trim(proj_filename), __FILE__, __LINE__)

         open(newunit=unit, file=trim(proj_filename), status='replace', action='write')

         ! Write header
         write(unit, '(A)') '# Projected Density of States'
         write(unit, '(A,A)') '# Method: ', trim(this%dos_method)
         if (trim(this%dos_method) == 'gaussian') then
            write(unit, '(A,F8.5,A)') '# Gaussian sigma: ', this%gaussian_sigma, ' Ry'
         end if
         write(unit, '(A,I0)') '# Energy points: ', this%n_energy_points
         write(unit, '(A,2F10.6,A)') '# Absolute energy range: ', this%dos_energy_range(1), &
                                      this%dos_energy_range(2), ' Ry'
         write(unit, '(A,F12.6,A)') '# Fermi level: ', this%fermi_level, ' Ry'
         write(unit, '(A,2F10.6,A)') '# Printed energy range E-E_F: ', &
                                      this%dos_energy_range(1) - this%fermi_level, &
                                      this%dos_energy_range(2) - this%fermi_level, ' Ry'
         write(unit, '(A)') '# Columns: Energy-E_F(Ry), s_up, p_up, d_up, f_up, s_down, p_down, d_down, f_down'

         ! Write projected DOS data (energy grid already in Ry)
         do i_energy = 1, this%n_energy_points
            write(unit, '(9F12.6)') this%dos_energy_grid(i_energy) - this%fermi_level, &
                                  (this%projected_dos(1, iorb, 1, i_energy), iorb=1,4), &  ! spin up: s,p,d,f
                                  (this%projected_dos(1, iorb, 2, i_energy), iorb=1,4)     ! spin down: s,p,d,f
         end do

         close(unit)
         call root_info('write_dos_to_file: Projected DOS written to file', __FILE__, __LINE__)
      end if

      ! Write band moments if available
      if (allocated(this%band_moments)) then
         proj_filename = 'band_moments.dat'
         call root_info('write_dos_to_file: Writing band moments to ' // trim(proj_filename), __FILE__, __LINE__)

         open(newunit=unit, file=trim(proj_filename), status='replace', action='write')

         ! Write header
         write(unit, '(A)') '# Band Moments'
         write(unit, '(A,A)') '# Method: ', trim(this%dos_method)
         if (trim(this%dos_method) == 'gaussian') then
            write(unit, '(A,F8.5,A)') '# Gaussian sigma: ', this%gaussian_sigma, ' Ry'
         end if
         write(unit, '(A,F8.4,A)') '# Temperature: ', this%temperature, ' K'
         write(unit, '(A,F12.6,A)') '# Fermi level: ', this%fermi_level, ' Ry'
         write(unit, '(A)') '# Format: site, orbital_type, spin, m0, m1, m2'
         write(unit, '(A)') '# Orbital types: 1=s, 2=p, 3=d, 4=f'
         write(unit, '(A)') '# Spin: 1=up, 2=down'
         write(unit, '(A)') '# m0 = integrated DOS up to Fermi level with Fermi-Dirac weighting'
         write(unit, '(A)') '# m1 = center of gravity (occupied states)'
         write(unit, '(A)') '# m2 = width of occupied states'

         ! Write band moments data
         do isite = 1, this%n_sites
            do iorb = 1, this%n_orb_types
               do ispin = 1, this%n_spin_components
                  write(unit, '(3I3,3F12.6)') isite, iorb, ispin, &
                                             this%band_moments(isite, iorb, ispin, 1), &  ! m0
                                             this%band_moments(isite, iorb, ispin, 2), &  ! m1
                                             this%band_moments(isite, iorb, ispin, 3)     ! m2
               end do
            end do
         end do

         close(unit)
         call root_info('write_dos_to_file: Band moments written to file', __FILE__, __LINE__)
      end if

      call root_info('write_dos_to_file: DOS written to file', __FILE__, __LINE__)
   end subroutine write_dos_to_file

   !> @brief Calculate band energy from integrated DOS moments.
   !> @param[in] this Reciprocal object containing band_moments.
   !> @return Band energy accumulated over projected channels.
   module function calculate_band_energy_from_moments(this) result(eband)
      class(reciprocal), intent(in) :: this
      real(rp) :: eband
      integer :: isite, iorb, ispin
      
      eband = 0.0_rp
      
      if (.not. allocated(this%band_moments)) return
      
      ! Band energy = ∫ E * DOS(E) * f(E) dE
      !             = sum over all orbitals of (m0 * m1)
      ! where m0 = occupation, m1 = band center
      do isite = 1, this%n_sites
         do iorb = 1, this%n_orb_types
            do ispin = 1, this%n_spin_components
               eband = eband + this%band_moments(isite, iorb, ispin, 1) * &
                              this%band_moments(isite, iorb, ispin, 2)
            end do
         end do
      end do
      
   end function calculate_band_energy_from_moments

   !> @brief Diagnostic total-DOS band-energy integral.
   !> @details This is deliberately not canonical: it inherits the selected
   !>          output window, grid, and broadening and exists only to quantify
   !>          their integration error against the occupied-eigenvalue sum.
   module function calculate_band_energy_from_total_dos(this) result(eband)
      class(reciprocal), intent(in) :: this
      real(rp) :: eband

      integer :: ie
      real(rp) :: e1, e2, f1, f2, kT, delta_e

      eband = 0.0_rp
      if (.not. allocated(this%total_dos) .or. .not. allocated(this%dos_energy_grid)) return
      kT = max(this%temperature*kB_Ry_per_K, 1.0e-10_rp)

      do ie = 1, size(this%dos_energy_grid) - 1
         e1 = this%dos_energy_grid(ie)
         e2 = this%dos_energy_grid(ie + 1)
         f1 = fermi_dirac_occupation(e1, this%fermi_level, kT)
         f2 = fermi_dirac_occupation(e2, this%fermi_level, kT)
         delta_e = e2 - e1
         eband = eband + 0.5_rp*delta_e*(e1*this%total_dos(ie)*f1 + &
                                         e2*this%total_dos(ie + 1)*f2)
      end do
   end function calculate_band_energy_from_total_dos

end submodule reciprocal_dos
