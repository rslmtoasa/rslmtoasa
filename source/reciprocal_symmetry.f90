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
      logical, allocatable :: hit(:)
      integer, allocatable :: src(:)
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
                  if (this%sym_atom_perm(i, k) /= 0) call g_logger%fatal('setup_magnetic_symmetry: sites ' // &
                     trim(int2str(this%sym_atom_perm(i, k))) // ' and ' // trim(int2str(j)) // ' coincide modulo a lattice vector; check that ' // &
                     'a(:,i) holds lattice vector i and that no two basis positions differ by a lattice vector', __FILE__, __LINE__)
                  this%sym_atom_perm(i, k) = j
               end if
            end do
         end do
      end do
      if (any(this%sym_atom_perm == 0)) call g_logger%fatal('setup_magnetic_symmetry: operation does not map the sites onto themselves', __FILE__, __LINE__)
      allocate(hit(n), src(n))
      do k = 1, this%n_sym_ops
         hit = .false.
         do i = 1, n
            j = this%sym_atom_perm(i, k)
            if (hit(j)) call g_logger%fatal('setup_magnetic_symmetry: sites ' // trim(int2str(src(j))) // ' and ' // trim(int2str(i)) // &
               ' map onto the same site; the site map of an operation is not one-to-one', __FILE__, __LINE__)
            hit(j) = .true.
            src(j) = i
         end do
      end do

      call root_info('setup_magnetic_symmetry: ' // trim(int2str(this%n_sym_ops)) // ' operations, ' // &
                     trim(int2str(count(this%sym_spin_flips < 0))) // ' time-reversal-combined', __FILE__, __LINE__)
   end subroutine setup_magnetic_symmetry

   !> Cartesian orthogonal matrix of each stored operation (lattice%cr is in units of alat, frac = bt cart/(2 pi)).
   function sym_cart_rotations(this) result(rc)
      class(reciprocal), intent(in) :: this
      real(rp) :: rc(3, 3, this%n_sym_ops)
      real(rp) :: bt(3, 3)
      integer :: k

      bt = transpose(this%reciprocal_vectors)
      do k = 1, this%n_sym_ops
         rc(:, :, k) = matmul(inverse_3x3(bt), matmul(real(this%sym_rotations(:, :, k), rp), bt))
      end do
   end function sym_cart_rotations

   !> @brief Cartesian action of each stored operation on a spin vector.
   !> @details SOC modes (nsp = 2, 4): det(R) R_cart, negated for time-reversal-combined operations (axial vector).
   !>          nsp = 1, 3: identity, the spin does not rotate.
   function sym_spin_matrices(this) result(w)
      class(reciprocal), intent(in) :: this
      real(rp) :: w(3, 3, this%n_sym_ops)
      real(rp) :: rc(3, 3, this%n_sym_ops)
      integer :: k, nsp

      nsp = 1
      if (associated(this%control)) nsp = this%control%nsp
      rc = sym_cart_rotations(this)
      do k = 1, this%n_sym_ops
         w(:, :, k) = 0.0_rp
         if (nsp == 2 .or. nsp == 4) then
            w(:, :, k) = rc(:, :, k)*real(nint(determinant(real(this%sym_rotations(:, :, k), rp))), rp) &
                         *real(this%sym_spin_flips(k), rp)
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

   !> @brief Operation of each full-mesh point (see the interface).
   module subroutine full_mesh_operations(this, op, w)
      class(reciprocal), intent(in) :: this
      integer, allocatable, intent(out) :: op(:)
      real(rp), allocatable, intent(out) :: w(:, :, :)
      real(rp), parameter :: k_tol = 1.0e-6_rp
      real(rp), allocatable :: kfull(:, :), rit(:, :, :)
      real(rp) :: kt(3), d(3)
      integer :: shift(3), kf, ki, o, itr, nk_full
      logical :: found, tr_only

      nk_full = product(this%nk_mesh)
      if (this%n_sym_ops == 0 .or. size(this%full_to_irred_k) /= nk_full) &
         call g_logger%fatal('full_mesh_operations: no symmetry group or full_to_irred_k size mismatch', __FILE__, __LINE__)
      shift = merge(1, 0, abs(this%k_offset) > 1.0e-12_rp)
      if (this%symmetry_analysis%spglib%get_full_kpoint_mesh_with_points(this%nk_mesh, shift, kfull) /= nk_full) &
         call g_logger%fatal('full_mesh_operations: full mesh unavailable', __FILE__, __LINE__)
      tr_only = this%use_time_reversal
      if (associated(this%control)) tr_only = tr_only .and. this%control%nsp == 1
      allocate(op(nk_full), rit(3, 3, this%n_sym_ops))
      do o = 1, this%n_sym_ops
         rit(:, :, o) = transpose(inverse_3x3(real(this%sym_rotations(:, :, o), rp)))
      end do

      do kf = 1, nk_full
         ki = this%full_to_irred_k(kf)
         found = .false.
         do o = 1, this%n_sym_ops
            do itr = 0, merge(1, 0, tr_only)
               kt = real(this%sym_spin_flips(o)*(1 - 2*itr), rp)*matmul(rit(:, :, o), this%k_points(:, ki))
               d = kt - kfull(:, kf)
               found = maxval(abs(d - nint(d))) < k_tol
               if (found) exit
            end do
            if (found) exit
         end do
         if (.not. found) call g_logger%fatal('full_mesh_operations: no operation maps irreducible k ' // &
                                              trim(int2str(ki)) // ' onto full k ' // trim(int2str(kf)), __FILE__, __LINE__)
         op(kf) = o
      end do
      allocate(w(3, 3, this%n_sym_ops))
      w = sym_spin_matrices(this)
   end subroutine full_mesh_operations

   !> @brief Band projections of one full-mesh point from the irreducible mesh (see the interface).
   module function full_mesh_projection(this, proj_irr, op, w, kf, l, isite, ib) result(p)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: proj_irr(:, :, :, :, :), w(:, :, :)
      integer, intent(in) :: op(:), kf, l, isite, ib
      real(rp) :: p(4)
      integer :: ki, o, s

      if (size(op) == 0) then
         p = proj_irr(:, l, isite, ib, kf)
         return
      end if
      ki = this%full_to_irred_k(kf)
      o = op(kf)
      s = findloc(this%sym_atom_perm(:, o), isite, 1)
      p(1) = proj_irr(1, l, s, ib, ki)
      p(2:4) = matmul(w(:, :, o), proj_irr(2:4, l, s, ib, ki))
   end function full_mesh_projection

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

   !> @brief Group average of the projected DOS curves of each site and l (see the interface).
   module subroutine symmetrize_projected_dos(this)
      class(reciprocal), intent(inout) :: this
      real(rp), allocatable :: w(:, :, :), n(:, :, :), n_avg(:, :, :), m_avg(:, :, :, :)
      real(rp) :: loc
      integer :: k, i, j, il, ie

      if (this%n_sym_ops == 0) return
      allocate(w(3, 3, this%n_sym_ops))
      w = sym_spin_matrices(this)
      n = this%projected_dos(:, :, 1, :) + this%projected_dos(:, :, 2, :)
      allocate(n_avg(size(n, 1), size(n, 2), size(n, 3)), m_avg(3, size(n, 1), size(n, 2), size(n, 3)))
      n_avg = 0.0_rp
      m_avg = 0.0_rp
      do k = 1, this%n_sym_ops
         do i = 1, this%n_sites
            j = this%sym_atom_perm(i, k)
            n_avg(j, :, :) = n_avg(j, :, :) + n(i, :, :)
            do il = 1, size(n, 2)
               do ie = 1, size(n, 3)
                  m_avg(:, j, il, ie) = m_avg(:, j, il, ie) + matmul(w(:, :, k), this%projected_dos_moments(i, il, 1, :, ie))
               end do
            end do
         end do
      end do
      n_avg = n_avg/real(this%n_sym_ops, rp)
      m_avg = m_avg/real(this%n_sym_ops, rp)

      this%dos_mx_tot = 0.0_rp
      this%dos_my_tot = 0.0_rp
      this%dos_mz_tot = 0.0_rp
      do ie = 1, size(n, 3)
         do il = 1, size(n, 2)
            do i = 1, this%n_sites
               loc = sum(this%dos_site_axes(:, i)*m_avg(:, i, il, ie))
               this%projected_dos(i, il, 1, ie) = 0.5_rp*(n_avg(i, il, ie) + loc)
               this%projected_dos(i, il, 2, ie) = 0.5_rp*(n_avg(i, il, ie) - loc)
               this%projected_dos_moments(i, il, 1, :, ie) = m_avg(:, i, il, ie)
               this%dos_mx_tot(ie) = this%dos_mx_tot(ie) + m_avg(1, i, il, ie)
               this%dos_my_tot(ie) = this%dos_my_tot(ie) + m_avg(2, i, il, ie)
               this%dos_mz_tot(ie) = this%dos_mz_tot(ie) + m_avg(3, i, il, ie)
            end do
         end do
      end do
   end subroutine symmetrize_projected_dos

   !> @brief Group average of per-site density matrices rho_a = 1/|G| sum_g U_g T_g[rho_(g^-1 a)] U_g^dagger.
   !> @details rho(:, :, site) is P_ab = psi_a psi_b^* over the nb eigenvector components of the site (spin block of norb
   !>          orbitals, spherical basis). U_g = D(R) for nsp = 1, 3 (spin untouched), D(R) x SU(2)(det(R) R) for nsp = 2, 4,
   !>          T_g the antiunitary i sigma_y K of the time-reversal-combined operations. The average runs in the real-harmonic
   !>          basis, where D is real and K is plain conjugation. Does nothing without a group (n_sym_ops = 0).
   module subroutine symmetrize_site_dm(this, rho)
      class(reciprocal), intent(in) :: this
      complex(rp), intent(inout) :: rho(:, :, :)
      complex(rp), allocatable :: u(:, :, :)
      complex(rp) :: d(norb, norb), sp(2, 2)
      real(rp), allocatable :: rc(:, :, :)
      integer :: k, i, s, t, nsp

      if (this%n_sym_ops == 0) return
      if (size(rho, 1) /= nb .or. size(rho, 2) /= nb) call g_logger%fatal('symmetrize_site_dm: matrix size is not nb', __FILE__, __LINE__)
      nsp = 1
      if (associated(this%control)) nsp = this%control%nsp
      allocate(u(nb, nb, this%n_sym_ops), rc(3, 3, this%n_sym_ops))
      rc = sym_cart_rotations(this)
      u = (0.0_rp, 0.0_rp)
      do k = 1, this%n_sym_ops
         call orbital_rotation(rc(:, :, k), d)
         sp = sigma_0
         if (nsp == 2 .or. nsp == 4) call spin_rotation(rc(:, :, k), sp)
         do s = 1, 2
            do t = 1, 2
               u((s - 1)*norb + 1:s*norb, (t - 1)*norb + 1:t*norb, k) = sp(s, t)*d
            end do
         end do
      end do
      do i = 1, size(rho, 3)
         do s = 1, 2
            do t = 1, 2
               call hcpx(rho((s - 1)*norb + 1:s*norb, (t - 1)*norb + 1:t*norb, i), 'sph2cart')
            end do
         end do
      end do
      call group_average_dm(rho, this%sym_atom_perm, u, this%sym_spin_flips < 0)
      do i = 1, size(rho, 3)
         do s = 1, 2
            do t = 1, 2
               call hcpx(rho((s - 1)*norb + 1:s*norb, (t - 1)*norb + 1:t*norb, i), 'cart2sph')
            end do
         end do
      end do
   end subroutine symmetrize_site_dm

end submodule reciprocal_symmetry
