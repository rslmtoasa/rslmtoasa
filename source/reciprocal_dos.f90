submodule (reciprocal_mod) reciprocal_dos
#ifdef USE_MPI
   use mpi
#endif
   implicit none

contains

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

      ! Local variables
      character(len=100) :: filename

      call root_info('calculate_density_of_states: Starting DOS calculation', __FILE__, __LINE__)

      ! Set parameters from optional arguments
      if (present(n_energy_points)) this%n_energy_points = n_energy_points
      if (present(energy_range)) this%dos_energy_range = energy_range
      if (present(method)) this%dos_method = trim(method)
      if (present(gaussian_sigma)) this%gaussian_sigma = gaussian_sigma
      if (present(temperature)) this%temperature = temperature
      if (present(fermi_level)) this%fermi_level = fermi_level
      if (present(total_electrons)) this%total_electrons = total_electrons
      if (present(auto_find_fermi)) this%auto_find_fermi = auto_find_fermi

      ! Set output filename
      filename = 'density_of_states.dat'
      if (present(output_file)) filename = output_file

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

      ! Symmetry-reduced tetra path:
      ! keep a correctness-first backend that switches to full mesh diagonalization
      ! for tetra/blochl to preserve exact SCF observables.
      if (this%use_symmetry_reduction) then
         if (trim(this%dos_method) == 'tetrahedron' .or. trim(this%dos_method) == 'blochl') then
            call this%ensure_tetra_symmetry_backend()
         end if
      end if

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

      if (this%auto_find_fermi .and. this%total_electrons > 0.0_rp) then
         this%fermi_level = this%find_fermi_level_from_dos(this%total_electrons)
         call g_logger%info('calculate_density_of_states: Auto-found Fermi level = ' // &
                           trim(real2str(this%fermi_level, '(F 10.6)')) // ' Ry', __FILE__, __LINE__)
      end if

      ! Calculate orbital projections and band moments
      call this%project_dos_orbitals()
      call this%calculate_band_moments()

      ! Write results to file
      call this%write_dos_to_file(filename)

      call root_info('calculate_density_of_states: DOS calculation completed', __FILE__, __LINE__)
   end subroutine calculate_density_of_states

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

   module subroutine calculate_dos_tetrahedron(this)
      class(reciprocal), intent(inout) :: this

      ! Local variables
      integer :: i_tet, i_corner, i_band, nbands
      real(rp), dimension(4) :: e_corners
      integer, allocatable :: tet_ir(:, :), tet_mult(:)
      integer :: n_tet_ir
      real(rp) :: tet_weight, dos_integral, nos_integral
      real(rp) :: fermi_count, fermi_error
      real(rp), allocatable :: local_dos(:), local_nos(:)
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
!$OMP& SHARED(this, tet_ir, tet_mult, tet_start, tet_end, nbands, nk_full_expected) &
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
                  e_corners(i_corner) = this%eigenvalues(i_band, tet_ir(i_corner, i_tet))
               else
                  e_corners(i_corner) = this%eigenvalues(i_band, this%tetrahedra(i_corner, i_tet))
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

      call g_logger%info('calculate_dos_tetrahedron: Tetrahedron DOS calculation completed', __FILE__, __LINE__)
   end subroutine calculate_dos_tetrahedron

   module function tetrahedron_dos_contribution(this, energy, e_sorted) result(dos)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: energy
      real(rp), dimension(4), intent(in) :: e_sorted
      real(rp) :: dos

      ! Local variables
      real(rp) :: e1, e2, e3, e4, vol_factor
      real(rp), parameter :: eps = 1.0e-12_rp

      ! Tetrahedron volume factor (adjusted for correct normalization)
      vol_factor = 0.4595_rp

      e1 = e_sorted(1)
      e2 = e_sorted(2)
      e3 = e_sorted(3)
      e4 = e_sorted(4)

      dos = 0.0_rp

      ! Handle different energy ranges
      if (energy < e1 - eps) then
         ! Energy below all eigenvalues
         dos = 0.0_rp
      else if (energy >= e1 - eps .and. energy < e2 - eps) then
         ! Energy in [e1, e2)
         dos = vol_factor * 3.0_rp * (energy - e1)**2 / ((e2 - e1) * (e3 - e1) * (e4 - e1))
      else if (energy >= e2 - eps .and. energy < e3 - eps) then
         ! Energy in [e2, e3)
         dos = vol_factor * (3.0_rp * (energy - e1)**2 - 6.0_rp * (energy - e1) * (energy - e2) + &
                           3.0_rp * (e3 - e1) * (energy - e2) * 2.0_rp / (e3 - e2) + &
                           3.0_rp * (energy - e1) * (energy - e2) * 2.0_rp / (e4 - e2)) / &
                           ((e3 - e1) * (e4 - e1))
      else if (energy >= e3 - eps .and. energy < e4 - eps) then
         ! Energy in [e3, e4)
         dos = vol_factor * 3.0_rp * (e4 - energy)**2 / ((e4 - e1) * (e4 - e2) * (e4 - e3))
      else
         ! Energy above all eigenvalues
         dos = 0.0_rp
      end if

   end function tetrahedron_dos_contribution

   module function blochl_dos_contribution(this, energy, e_sorted) result(dos)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: energy
      real(rp), dimension(4), intent(in) :: e_sorted
      real(rp) :: dos

      ! Local variables
      real(rp) :: e1, e2, e3, e4, C
      real(rp), parameter :: TOL = 1.0e-10_rp

      e1 = e_sorted(1)
      e2 = e_sorted(2)
      e3 = e_sorted(3)
      e4 = e_sorted(4)

      dos = 0.0_rp

      ! Skip if degenerate (would cause division by zero)
      if (abs(e2-e1) < TOL .or. abs(e3-e1) < TOL .or. abs(e4-e1) < TOL .or. &
          abs(e4-e2) < TOL .or. abs(e4-e3) < TOL .or. abs(e3-e2) < TOL) then
         return
      end if

      ! Blöchl modified tetrahedron DOS contribution
      if (energy <= e1) then
         dos = 0.0_rp
      else if (energy >= e4) then
         dos = 0.0_rp
      else if (energy <= e2) then
         ! Region I: e1 < E <= e2
         ! Blöchl Eq. (23): D(E) = 3(E-e1)²/[(e4-e1)(e3-e1)(e2-e1)]
         dos = 3.0_rp * (energy - e1)**2 / ((e4 - e1) * (e3 - e1) * (e2 - e1))
      else if (energy <= e3) then
         ! Region II: e2 < E <= e3
         ! Blöchl Eq. (24): D(E) = 1/[(e4-e1)(e3-e1)] * 
         !   [3(e2-e1) + 6(E-e2) - 3(e3+e4-e1-e2)(E-e2)²/[(e3-e2)(e4-e2)]]
         C = 1.0_rp / ((e4 - e1) * (e3 - e1))
         dos = C * (3.0_rp * (e2 - e1) + 6.0_rp * (energy - e2) - &
                   3.0_rp * (e3 + e4 - e1 - e2) * (energy - e2)**2 / &
                   ((e3 - e2) * (e4 - e2)))
      else
         ! Region III: e3 < E < e4
         ! Blöchl Eq. (25): D(E) = 3(e4-E)²/[(e4-e1)(e4-e2)(e4-e3)]
         dos = 3.0_rp * (e4 - energy)**2 / ((e4 - e1) * (e4 - e2) * (e4 - e3))
      end if

   end function blochl_dos_contribution

   module subroutine sort_real_array(arr, sorted, indices)
      real(rp), dimension(:), intent(in) :: arr
      real(rp), dimension(:), intent(out) :: sorted
      integer, dimension(:), intent(out) :: indices

      ! Local variables
      integer :: i, j, n, temp_idx
      real(rp) :: temp_val

      n = size(arr)
      sorted = arr
      do i = 1, n
         indices(i) = i
      end do

      ! Simple bubble sort
      do i = 1, n-1
         do j = 1, n-i
            if (sorted(j) > sorted(j+1)) then
               temp_val = sorted(j)
               sorted(j) = sorted(j+1)
               sorted(j+1) = temp_val
               temp_idx = indices(j)
               indices(j) = indices(j+1)
               indices(j+1) = temp_idx
            end if
         end do
      end do
   end subroutine sort_real_array

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

   module subroutine tetra_add_nos(volwgt, ecorn_in, emin, emax, nos, npts)
      real(rp), intent(in) :: volwgt, ecorn_in(4), emin, emax
      integer, intent(in) :: npts
      real(rp), intent(inout) :: nos(npts)

      integer :: i, i0(4), m
      real(rp) :: ecorn(4), de, e1, e2, e3, e4
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

      do m = 1, 4
         if (ecorn(m) <= emin) then
            i0(m) = 1
         else if (ecorn(m) > emax) then
            i0(m) = npts + 1
         else
            i0(m) = 2 + int((ecorn(m) - emin) / de - 1.0e-12_rp)
         end if
      end do

      if (i0(1) < i0(2) .and. abs((e2-e1)*(e3-e1)*(e4-e1)) > tol) then
         c3 = volwgt / ((e2-e1)*(e3-e1)*(e4-e1))
         x0 = emin - e1 - de
         do i = i0(1), i0(2) - 1
            x = x0 + real(i, rp) * de
            nos(i) = nos(i) + c3*x*x*x
         end do
      end if
      if (i0(2) < i0(3) .and. abs((e3-e1)*(e4-e1)*(e3-e2)*(e4-e2)) > tol) then
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
      if (i0(3) < i0(4) .and. abs((e4-e3)*(e4-e2)*(e4-e1)) > tol) then
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

   module subroutine tetra_add_dos(volwgt, ecorn_in, emin, emax, dos, npts)
      real(rp), intent(in) :: volwgt, ecorn_in(4), emin, emax
      integer, intent(in) :: npts
      real(rp), intent(inout) :: dos(npts)

      integer :: i, i0(4), m
      real(rp) :: ecorn(4), de, e1, e2, e3, e4
      real(rp) :: c1, c2, c3, x
      real(rp), parameter :: tol = 1.0e-12_rp

      call sort4(ecorn_in, ecorn)
      e1 = ecorn(1); e2 = ecorn(2); e3 = ecorn(3); e4 = ecorn(4)
      if (e4 < emin .or. e1 > emax) return

      de = (emax - emin) / real(npts - 1, rp)
      if (de <= 0.0_rp) return

      do m = 1, 4
         if (ecorn(m) < emin) then
            i0(m) = 1
         else if (ecorn(m) > emax) then
            i0(m) = npts + 1
         else
            i0(m) = 2 + int((ecorn(m) - emin) / de - 1.0e-12_rp)
         end if
      end do

      if (i0(1) < i0(2) .and. abs((e2-e1)*(e3-e1)*(e4-e1)) > tol) then
         c3 = volwgt*3.0_rp/((e2-e1)*(e3-e1)*(e4-e1))
         do i = i0(1), i0(2) - 1
            x = emin - e1 + real(i - 1, rp) * de
            dos(i) = dos(i) + c3*x*x
         end do
      end if
      if (i0(2) < i0(3) .and. abs((e3-e1)*(e4-e1)*(e3-e2)*(e4-e2)) > tol) then
         c3 = volwgt*3.0_rp*(e1+e2-e3-e4)/((e3-e1)*(e4-e1)*(e3-e2)*(e4-e2))
         c2 = volwgt*6.0_rp/((e3-e1)*(e4-e1))
         c1 = c2*(e2-e1)*0.5_rp
         do i = i0(2), i0(3) - 1
            x = emin - e2 + real(i - 1, rp) * de
            dos(i) = dos(i) + c1 + x*(c2 + x*c3)
         end do
      end if
      if (i0(3) < i0(4) .and. abs((e4-e3)*(e4-e2)*(e4-e1)) > tol) then
         c3 = volwgt*3.0_rp/((e4-e3)*(e4-e2)*(e4-e1))
         do i = i0(3), i0(4) - 1
            x = emin - e4 + real(i - 1, rp) * de
            dos(i) = dos(i) + c3*x*x
         end do
      end if
   end subroutine tetra_add_dos

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

   module subroutine calculate_dos_blochl(this)
      class(reciprocal), intent(inout) :: this

      call g_logger%info('calculate_dos_blochl: Using LMTO-style tetrahedron DOS/NOS backend for scalar DOS.', __FILE__, __LINE__)
      call this%calculate_dos_tetrahedron()
      call g_logger%info('calculate_dos_blochl: Blöchl-compatible scalar DOS calculation completed', __FILE__, __LINE__)
   end subroutine calculate_dos_blochl

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

   module subroutine expand_eigenvalues_to_full_mesh(this)
      class(reciprocal), intent(inout) :: this
      integer :: ik_full, nk_full, nk_irred, nbands, ik_irred
      real(rp), dimension(:, :), allocatable :: eigenvalues_full
      complex(rp), dimension(:, :, :), allocatable :: eigenvectors_full
      integer :: nrow_evec, nband_evec, nk_evec
      
      if (.not. this%use_symmetry_reduction) then
         call g_logger%info('expand_eigenvalues_to_full_mesh: No symmetry reduction, skipping', __FILE__, __LINE__)
         return
      end if
      
      nk_irred = this%nk_total  ! Current number of irreducible k-points
      nk_full = this%nk_mesh(1) * this%nk_mesh(2) * this%nk_mesh(3)  ! Full mesh size
      nbands = size(this%eigenvalues, 1)
      
      call g_logger%info('expand_eigenvalues_to_full_mesh: Expanding ' // trim(int2str(nk_irred)) // &
                        ' irreducible k-points to ' // trim(int2str(nk_full)) // ' full mesh points', __FILE__, __LINE__)
      
      ! Allocate full mesh eigenvalues
      allocate(eigenvalues_full(nbands, nk_full))
      eigenvalues_full = 0.0_rp
      
      if (.not. allocated(this%full_to_irred_k)) then
         if (this%strict_symmetry_checks) then
            call g_logger%fatal('expand_eigenvalues_to_full_mesh: full_to_irred mapping not available', __FILE__, __LINE__)
         else
            call g_logger%warning('expand_eigenvalues_to_full_mesh: full_to_irred mapping not available; fallback copy', __FILE__, __LINE__)
            do ik_irred = 1, min(nk_irred, nk_full)
               eigenvalues_full(:, ik_irred) = this%eigenvalues(:, ik_irred)
            end do
         end if
      else
         do ik_full = 1, nk_full
            ik_irred = this%full_to_irred_k(ik_full)
            if (ik_irred >= 1 .and. ik_irred <= nk_irred) then
               eigenvalues_full(:, ik_full) = this%eigenvalues(:, ik_irred)
            else
               if (this%strict_symmetry_checks) then
                  call g_logger%fatal('expand_eigenvalues_to_full_mesh: invalid mapping index ' // trim(int2str(ik_irred)), __FILE__, __LINE__)
               end if
            end if
         end do
      end if
      
      ! Store expanded eigenvalues
      deallocate(this%eigenvalues)
      allocate(this%eigenvalues(nbands, nk_full))
      this%eigenvalues = eigenvalues_full

      ! If eigenvectors are available on irreducible mesh, expand them with
      ! the same full->irred mapping so projection/moment integrations remain
      ! consistent with expanded eigenvalues.
      if (allocated(this%eigenvectors)) then
         nrow_evec = size(this%eigenvectors, 1)
         nband_evec = size(this%eigenvectors, 2)
         nk_evec = size(this%eigenvectors, 3)
         if (nband_evec == nbands .and. nk_evec == nk_irred) then
            allocate(eigenvectors_full(nrow_evec, nband_evec, nk_full))
            eigenvectors_full = cmplx(0.0_rp, 0.0_rp, rp)
            if (allocated(this%full_to_irred_k)) then
               do ik_full = 1, nk_full
                  ik_irred = this%full_to_irred_k(ik_full)
                  if (ik_irred >= 1 .and. ik_irred <= nk_irred) then
                     eigenvectors_full(:, :, ik_full) = this%eigenvectors(:, :, ik_irred)
                  end if
               end do
            end if
            deallocate(this%eigenvectors)
            allocate(this%eigenvectors(nrow_evec, nband_evec, nk_full))
            this%eigenvectors = eigenvectors_full
            deallocate(eigenvectors_full)
         else
            call g_logger%warning('expand_eigenvalues_to_full_mesh: eigenvector dimensions do not match irreducible eigenvalue mesh; skipping eigenvector expansion', __FILE__, __LINE__)
         end if
      end if
      
      ! Update nk_total to full mesh
      this%nk_total = nk_full
      ! Keep k-points/weights/maps consistent with the expanded full mesh.
      call this%generate_mp_mesh()
      
      deallocate(eigenvalues_full)
      
      call g_logger%info('expand_eigenvalues_to_full_mesh: Expansion complete using explicit full_to_irred map', __FILE__, __LINE__)
   end subroutine expand_eigenvalues_to_full_mesh

   module subroutine calculate_dos_tetrahedron_with_symmetry(this)
      class(reciprocal), intent(inout) :: this
      
      if (this%use_symmetry_reduction) then
         call g_logger%info('calculate_dos_tetrahedron_with_symmetry: Using configured tetra symmetry backend', __FILE__, __LINE__)
         call this%ensure_tetra_symmetry_backend()
      end if
      
      call this%calculate_dos_tetrahedron()
   end subroutine calculate_dos_tetrahedron_with_symmetry

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

   module subroutine print_total_and_spin_dos(this)
      class(reciprocal), intent(in) :: this

      integer :: ie, n_energy, isite, iorb, ispin
      real(rp), allocatable :: fermi_dist(:)
      real(rp), allocatable :: energy_grid(:)
      real(rp), allocatable :: integrand_up(:), integrand_dn(:)
      real(rp) :: kT, fermi_arg
      real(rp) :: electrons_up, electrons_dn, electrons_total
      real(rp), parameter :: kB_Ry_per_K = 6.3336814e-6_rp

      if (.not. allocated(this%total_dos) .or. .not. allocated(this%dos_energy_grid)) then
         call g_logger%warning('print_total_and_spin_dos: DOS not available to print', __FILE__, __LINE__)
         return
      end if

      n_energy = this%n_energy_points
      allocate(fermi_dist(n_energy))
      allocate(energy_grid(n_energy))
      energy_grid = this%dos_energy_grid

      ! Boltzmann factor (Ry)
      kT = this%temperature * kB_Ry_per_K

      ! Pre-calc fermi distribution on energy grid
      do ie = 1, n_energy
         fermi_arg = (energy_grid(ie) - this%fermi_level)
         if (kT > 1.0e-10_rp) then
            fermi_dist(ie) = 1.0_rp / (exp(fermi_arg / kT) + 1.0_rp)
         else
            if (energy_grid(ie) <= this%fermi_level) then
               fermi_dist(ie) = 1.0_rp
            else
               fermi_dist(ie) = 0.0_rp
            end if
         end if
      end do

      ! Initialize integrands
      allocate(integrand_up(n_energy))
      allocate(integrand_dn(n_energy))
      integrand_up = 0.0_rp
      integrand_dn = 0.0_rp

      if (allocated(this%projected_dos)) then
         ! Sum projected DOS across sites and orbitals for each spin
         do isite = 1, this%n_sites
            do iorb = 1, this%n_orb_types
               if (this%n_spin_components >= 1) then
                  integrand_up = integrand_up + this%projected_dos(isite, iorb, 1, :)
               end if
               if (this%n_spin_components >= 2) then
                  integrand_dn = integrand_dn + this%projected_dos(isite, iorb, 2, :)
               end if
            end do
         end do
      else
         ! No projected DOS available: fall back to total DOS split by spin if possible
         if (this%n_spin_components == 2) then
            ! If no projection but spin resolved, attempt to split total DOS equally
            integrand_up = 0.5_rp * this%total_dos
            integrand_dn = 0.5_rp * this%total_dos
         else
            integrand_up = this%total_dos
            integrand_dn = 0.0_rp
         end if
      end if

      ! Apply Fermi weighting and integrate
      integrand_up = integrand_up * fermi_dist
      integrand_dn = integrand_dn * fermi_dist

      electrons_up = trapezoidal_integral(energy_grid, integrand_up)
      electrons_dn = trapezoidal_integral(energy_grid, integrand_dn)
      electrons_total = electrons_up + electrons_dn

      ! Log results
      if (this%n_spin_components >= 2) then
         call g_logger%info('DOS summary: Total electrons (from DOS) = ' // trim(real2str(electrons_total, '(F10.6)')) // &
                           ', Up = ' // trim(real2str(electrons_up, '(F10.6)')) // &
                           ', Down = ' // trim(real2str(electrons_dn, '(F10.6)')) // &
                           ', M = ' // trim(real2str(electrons_up - electrons_dn, '(F10.6)')), __FILE__, __LINE__)
      else
         call g_logger%info('DOS summary: Total electrons (from DOS) = ' // trim(real2str(electrons_total, '(F10.6)')), __FILE__, __LINE__)
      end if

      deallocate(fermi_dist, energy_grid, integrand_up, integrand_dn)
   end subroutine print_total_and_spin_dos

module function find_fermi_level_from_dos(this, total_electrons) result(fermi_level)
   class(reciprocal), intent(in) :: this
   real(rp), intent(in) :: total_electrons
   real(rp) :: fermi_level

   integer :: ie, max_iter
   real(rp) :: energy, kT
   real(rp) :: e_min, e_max, e_mid, electrons_at_e
   real(rp) :: q1, q2, e1, e2, step
   real(rp), parameter :: eV_to_Ry = 0.073498618_rp
   real(rp), parameter :: kB_Ry_per_K = 6.3336814e-6_rp

   call root_info('find_fermi_level_from_dos: Finding Fermi level for ' // &
                  trim(real2str(total_electrons, '(F 8.5)')) // ' electrons at T = ' // &
                  trim(real2str(this%temperature, '(F 8.3)')) // ' K', __FILE__, __LINE__)

   if (.not. allocated(this%total_dos)) then
      call g_logger%error('find_fermi_level_from_dos: Total DOS not calculated', __FILE__, __LINE__)
      fermi_level = 0.0_rp
      return
   end if

   ! Boltzmann constant in Ry/K
   kT = this%temperature * kB_Ry_per_K

   ! Energy range in Ry (dos_energy_grid is already in Ry, no conversion needed)
   e_min = this%dos_energy_grid(1)
   e_max = this%dos_energy_grid(this%n_energy_points)
   max_iter = 100

   if (allocated(this%total_nos) .and. size(this%total_nos) == this%n_energy_points) then
      if (this%total_nos(1) - total_electrons > 1.0e-8_rp) then
         call g_logger%warning('find_fermi_level_from_dos: requested electron count lies below integrated DOS window', __FILE__, __LINE__)
         fermi_level = e_min
         return
      end if
      if (this%total_nos(this%n_energy_points) - total_electrons < -1.0e-8_rp) then
         call g_logger%warning('find_fermi_level_from_dos: requested electron count lies above integrated DOS window', __FILE__, __LINE__)
         fermi_level = e_max
         return
      end if
      do ie = 1, this%n_energy_points - 1
         if (this%total_nos(ie) >= total_electrons - 1.0e-10_rp) then
            fermi_level = this%dos_energy_grid(ie)
            electrons_at_e = this%total_nos(ie)
            call root_info('find_fermi_level_from_dos: Found Fermi level from tetra N(E) at HOMO-side plateau ' // &
                           trim(real2str(fermi_level, '(F 8.5)')) // ' Ry (integrated ' // &
                           trim(real2str(electrons_at_e, '(F 8.5)')) // ' electrons)', __FILE__, __LINE__)
            return
         end if
         if (this%total_nos(ie + 1) >= total_electrons) then
            q1 = this%total_nos(ie)
            q2 = this%total_nos(ie + 1)
            e1 = this%dos_energy_grid(ie)
            e2 = this%dos_energy_grid(ie + 1)
            if (abs(q2 - q1) > 1.0e-14_rp) then
               fermi_level = e1 + (total_electrons - q1) * (e2 - e1) / (q2 - q1)
            else
               fermi_level = e1
            end if
            step = e2 - e1
            electrons_at_e = q1 + (q2 - q1) * (fermi_level - e1) / max(step, tiny(1.0_rp))
            call root_info('find_fermi_level_from_dos: Found Fermi level from tetra N(E) at ' // &
                           trim(real2str(fermi_level, '(F 8.5)')) // ' Ry (integrated ' // &
                           trim(real2str(electrons_at_e, '(F 8.5)')) // ' electrons)', __FILE__, __LINE__)
            return
         end if
      end do
   end if

   ! Bisection method
   do ie = 1, max_iter
      e_mid = (e_min + e_max) / 2.0_rp
      electrons_at_e = this%integrate_dos_up_to_energy(e_mid, kT)

      ! ! DEBUG: Print first few and last iterations
      ! if (ie <= 5 .or. ie >= max_iter-2) then
      !    call g_logger%info('  Bisection iter ' // trim(int2str(ie)) // ': E=' // &
      !                      trim(real2str(e_mid, '(F10.6)')) // ' Ry, electrons=' // &
      !                      trim(real2str(electrons_at_e, '(F12.8)')) // ', target=' // &
      !                      trim(real2str(total_electrons, '(F12.8)')), __FILE__, __LINE__)
      ! end if

      if (abs(electrons_at_e - total_electrons) < 1.0e-6_rp) then
         fermi_level = e_mid
         call root_info('  Bisection converged at iteration ' // trim(int2str(ie)), __FILE__, __LINE__)
         exit
      else if (electrons_at_e < total_electrons) then
         e_min = e_mid
      else
         e_max = e_mid
      end if
   end do

   if (ie >= max_iter) then
      call g_logger%warning('  Bisection did NOT converge after ' // trim(int2str(max_iter)) // ' iterations!', __FILE__, __LINE__)
   end if

   fermi_level = e_mid

   ! Final check
   electrons_at_e = this%integrate_dos_up_to_energy(fermi_level, kT)
   call root_info('find_fermi_level_from_dos: Found Fermi level at ' // &
                  trim(real2str(fermi_level, '(F 8.5)')) // ' Ry (integrated ' // &
                  trim(real2str(electrons_at_e, '(F 8.5)')) // ' electrons)', __FILE__, __LINE__)
   
   ! DEBUG: Check total DOS integral
   call root_info('find_fermi_level_from_dos: DEBUG - Checking DOS integral...', __FILE__, __LINE__)
   call root_info('  Total electrons requested: ' // trim(real2str(total_electrons, '(F10.6)')), __FILE__, __LINE__)
   call root_info('  Total electrons found: ' // trim(real2str(electrons_at_e, '(F10.6)')), __FILE__, __LINE__)
   call root_info('  Error: ' // trim(real2str(electrons_at_e - total_electrons, '(F10.6)')) // &
                  ' (' // trim(real2str(100.0_rp*(electrons_at_e - total_electrons)/total_electrons, '(F6.3)')) // '%)', &
                  __FILE__, __LINE__)
end function find_fermi_level_from_dos

module function integrate_dos_up_to_energy(this, energy, kT) result(integral)
   class(reciprocal), intent(in) :: this
   real(rp), intent(in) :: energy, kT
   real(rp) :: integral

   integer :: ie
   real(rp) :: e, fermi_weight, fermi_weight_next, delta_e
   real(rp), parameter :: eV_to_Ry = 0.073498618_rp

   integral = 0.0_rp

   do ie = 1, this%n_energy_points - 1
      ! dos_energy_grid is already in Ry, no conversion needed
      e = this%dos_energy_grid(ie)
      delta_e = this%dos_energy_grid(ie+1) - this%dos_energy_grid(ie)

      ! Fermi-Dirac weight at current energy
      if (kT > 1.0e-10_rp) then
         fermi_weight = 1.0_rp / (exp((e - energy) / kT) + 1.0_rp)
         fermi_weight_next = 1.0_rp / (exp((this%dos_energy_grid(ie+1) - energy) / kT) + 1.0_rp)
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

   module subroutine write_dos_to_file(this, filename)
      class(reciprocal), intent(in) :: this
      character(len=*), intent(in) :: filename

      ! Local variables
      integer :: unit, i_energy, isite, iorb, ispin
      character(len=256) :: proj_filename
      real(rp), parameter :: eV_to_Ry = 0.073498618_rp

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

   module function calculate_gaussian_weight_single(this, grid_energy, eigenvalue) result(weight)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: grid_energy, eigenvalue
      real(rp) :: weight

      ! Gaussian smearing: weight = exp(-(E_grid - E_eigen)²/(2σ²)) / (σ√(2π))
      real(rp) :: prefactor, exponent, delta_e

      delta_e = grid_energy - eigenvalue
      prefactor = 1.0_rp / (this%gaussian_sigma * sqrt(2.0_rp * 3.141592653589793_rp))
      exponent = -0.5_rp * (delta_e / this%gaussian_sigma)**2

      if (exponent > -20.0_rp) then  ! Avoid underflow
         weight = prefactor * exp(exponent)
      else
         weight = 0.0_rp
      end if
   end function calculate_gaussian_weight_single

   module pure subroutine sort_eigenvalues(e_in, e_sorted, sort_idx)
      real(rp), dimension(4), intent(in) :: e_in
      real(rp), dimension(4), intent(out) :: e_sorted
      integer, dimension(4), intent(out) :: sort_idx
      
      ! Local variables
      integer :: i, j, min_idx
      real(rp) :: temp_e
      integer :: temp_idx
      
      ! Initialize
      e_sorted = e_in
      do i = 1, 4
         sort_idx(i) = i
      end do
      
      ! Simple selection sort
      do i = 1, 3
         min_idx = i
         do j = i + 1, 4
            if (e_sorted(j) < e_sorted(min_idx)) then
               min_idx = j
            end if
         end do
         
         if (min_idx /= i) then
            ! Swap energies
            temp_e = e_sorted(i)
            e_sorted(i) = e_sorted(min_idx)
            e_sorted(min_idx) = temp_e
            
            ! Swap indices
            temp_idx = sort_idx(i)
            sort_idx(i) = sort_idx(min_idx)
            sort_idx(min_idx) = temp_idx
         end if
      end do
   end subroutine sort_eigenvalues

end submodule reciprocal_dos
