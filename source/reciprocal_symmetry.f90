!> Symmetry group of the site moments, behind the symmetry-reduced k mesh.
submodule (reciprocal_mod) reciprocal_symmetry
   use spin_density_mod, only: sd_cartesian_from_matrix, sd_matrix_from_cartesian
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

   !> @brief Cartesian action of each stored operation on a spin vector.
   !> @details SOC modes (nsp = 2, 4): det(R) R_cart, negated for time-reversal-combined operations (axial vector).
   !>          nsp = 1, 3: identity, the spin does not rotate.
   function sym_spin_matrices(this) result(w)
      class(reciprocal), intent(in) :: this
      real(rp) :: w(3, 3, this%n_sym_ops)
      real(rp) :: bt(3, 3), rf(3, 3)
      integer :: k, nsp

      nsp = 1
      if (associated(this%control)) nsp = this%control%nsp
      bt = transpose(this%reciprocal_vectors)   ! frac = bt cart/(2 pi), so R_cart = bt^-1 R_frac bt
      do k = 1, this%n_sym_ops
         w(:, :, k) = 0.0_rp
         if (nsp == 2 .or. nsp == 4) then
            rf = real(this%sym_rotations(:, :, k), rp)
            w(:, :, k) = matmul(inverse_3x3(bt), matmul(rf, bt))
            w(:, :, k) = w(:, :, k)*real(nint(determinant(rf)), rp)*real(this%sym_spin_flips(k), rp)
         else
            w(1, 1, k) = 1.0_rp
            w(2, 2, k) = 1.0_rp
            w(3, 3, k) = 1.0_rp
         end if
      end do
   end function sym_spin_matrices

   !> @brief Average the accumulated spin density over the symmetry group: rho_a = 1/|G| sum_g g[rho_(g^-1 a)].
   !> @details Charge n is kept, the spin vector m of each 2x2 block is rotated by sym_spin_matrices. Does nothing
   !>          without a group (n_sym_ops = 0).
   module subroutine symmetrize_spin_density(this)
      class(reciprocal), intent(inout) :: this
      complex(rp), allocatable :: avg(:, :, :, :, :)
      complex(rp) :: blk(2, 2)
      real(rp), allocatable :: w(:, :, :)
      real(rp) :: n, m(3)
      integer :: k, i, il, io

      if (this%n_sym_ops == 0) return
      allocate(w(3, 3, this%n_sym_ops))
      w = sym_spin_matrices(this)
      allocate(avg(size(this%rf_density%rho, 1), size(this%rf_density%rho, 2), size(this%rf_density%rho, 3), &
                   size(this%rf_density%rho, 4), size(this%rf_density%rho, 5)))
      avg = (0.0_rp, 0.0_rp)
      do k = 1, this%n_sym_ops
         do i = 1, this%n_sites
            do io = 1, size(avg, 5)
               do il = 1, size(avg, 3)
                  call sd_cartesian_from_matrix(this%rf_density%rho(:, :, il, i, io), n, m)
                  call sd_matrix_from_cartesian(n, matmul(w(:, :, k), m), blk)
                  avg(:, :, il, this%sym_atom_perm(i, k), io) = avg(:, :, il, this%sym_atom_perm(i, k), io) + blk
               end do
            end do
         end do
      end do
      this%rf_density%rho = avg/real(this%n_sym_ops, rp)
   end subroutine symmetrize_spin_density

   !> @brief Group average of one spin vector per site, with the rule of symmetrize_spin_density.
   !> @param[inout] v Spin vector of each site [3, nrec].
   module subroutine symmetrize_site_moments(this, v)
      class(reciprocal), intent(in) :: this
      real(rp), intent(inout) :: v(:, :)
      real(rp), allocatable :: w(:, :, :), avg(:, :)
      integer :: k, i

      if (this%n_sym_ops == 0) return
      allocate(w(3, 3, this%n_sym_ops), avg(3, size(v, 2)))
      w = sym_spin_matrices(this)
      avg = 0.0_rp
      do k = 1, this%n_sym_ops
         do i = 1, size(v, 2)
            avg(:, this%sym_atom_perm(i, k)) = avg(:, this%sym_atom_perm(i, k)) + matmul(w(:, :, k), v(:, i))
         end do
      end do
      v = avg/real(this%n_sym_ops, rp)
   end subroutine symmetrize_site_moments

end submodule reciprocal_symmetry
