!> Symmetry group of the site moments, behind the symmetry-reduced k mesh.
submodule (reciprocal_mod) reciprocal_symmetry
   implicit none
contains

   !> @brief Build and store the symmetry group of the current site moments.
   !> @details The group is fixed from the moments at the time of the call; call again to rebuild it if the
   !>          moments rotate away from it. nsp = 2, 4 (spin-orbit): magnetic group, moments as axial vectors.
   !>          nsp = 1, 3: crystal operations mapping each site onto a site with the same moment vector.
   !>          A moment is the site axis (get_site_spin_axis: constraint target, else the unit potential%mom,
   !>          whose sign is the sign of the site magnetisation), zero for force_nonmagnetic sites.
   !> @param[inout] this Reciprocal object receiving the n_sym_ops, sym_* group arrays.
   module subroutine setup_magnetic_symmetry(this)
      class(reciprocal), intent(inout) :: this
      real(rp), parameter :: pos_tol = 1.0e-4_rp
      real(rp), allocatable :: frac(:, :), moments(:, :)
      real(rp) :: x(3), d(3)
      integer :: nsp, n, i, j, k

      nsp = 1
      if (associated(this%control)) nsp = this%control%nsp
      n = this%lattice%nrec
      allocate(frac(3, n), moments(3, n))
      do i = 1, n
         ! lattice%cr is Cartesian in units of alat; reciprocal_vectors are 2 pi times the dual basis of lattice%a
         frac(:, i) = matmul(transpose(this%reciprocal_vectors), this%lattice%cr(:, i))/(2.0_rp*acos(-1.0_rp))
         call get_site_spin_axis(this, i, moments(:, i))
         if (this%lattice%symbolic_atoms(this%lattice%nbulk + i)%force_nonmagnetic) moments(:, i) = 0.0_rp
      end do

      this%n_sym_ops = this%symmetry_analysis%spglib%get_magnetic_symmetry(frac, this%lattice%iz(1:n), moments, &
                          nsp == 2 .or. nsp == 4, this%sym_rotations, this%sym_translations, this%sym_spin_flips)
      if (this%n_sym_ops == 0) call g_logger%fatal('setup_magnetic_symmetry: spglib returned no symmetry operations', __FILE__, __LINE__)

      if (allocated(this%sym_atom_perm)) deallocate(this%sym_atom_perm)
      allocate(this%sym_atom_perm(n, this%n_sym_ops))
      this%sym_atom_perm = 0
      do k = 1, this%n_sym_ops
         do i = 1, n
            x = matmul(this%sym_rotations(:, :, k), frac(:, i)) + this%sym_translations(:, k)
            do j = 1, n
               d = x - frac(:, j)
               if (maxval(abs(d - nint(d))) < pos_tol) then
                  this%sym_atom_perm(i, k) = j
                  exit
               end if
            end do
         end do
      end do
      if (any(this%sym_atom_perm == 0)) call g_logger%fatal('setup_magnetic_symmetry: operation does not map the sites onto themselves', __FILE__, __LINE__)

      call root_info('setup_magnetic_symmetry: ' // trim(int2str(this%n_sym_ops)) // ' operations, ' // &
                     trim(int2str(count(this%sym_spin_flips < 0))) // ' time-reversal-combined', __FILE__, __LINE__)
   end subroutine setup_magnetic_symmetry

end submodule reciprocal_symmetry
