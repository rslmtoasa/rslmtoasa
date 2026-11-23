submodule (reciprocal_mod) reciprocal_projection
#ifdef USE_MPI
   use mpi
#endif
   use magnetic_representation_mod, only: gbt_single_q
   implicit none

   !> First and last orbital index of the s, p, d, f blocks within one spin
   integer, parameter :: lstart(4) = [1, 2, 5, 10], lend(4) = [1, 4, 9, 16]

contains

   !> @brief Dispatch projected-DOS calculation for the active DOS method.
   !> @param[inout] this Reciprocal object receiving projected DOS arrays.
   module subroutine project_dos_orbitals(this)
      class(reciprocal), intent(inout) :: this

      call root_info('project_dos_orbitals: Starting orbital projection calculation', __FILE__, __LINE__)

      ! Use tetrahedron or Gaussian method based on dos_method
      if (trim(this%dos_method) == 'tetrahedron' .or. trim(this%dos_method) == 'blochl') then
         call this%project_dos_orbitals_tetrahedron()
      else
         call this%project_dos_orbitals_gaussian()
      end if
   end subroutine project_dos_orbitals

!> @brief Calculate orbital- and spin-projected DOS with Gaussian broadening.
!> @details Projects eigenvector weights onto site/orbital/spin channels and
!>          accumulates both scalar and directional spin DOS components.
!> @param[inout] this Reciprocal object receiving projected DOS arrays.
module subroutine project_dos_orbitals_gaussian(this)
   class(reciprocal), intent(inout) :: this
   integer :: ik, ik_global, ib, ie, iorb, isite
   integer :: n_orb_site
   real(rp) :: weight, orbital_char, energy
   real(rp) :: gaussian_weight, sigma_squared, sigma_use
   integer :: nbands
   real(rp) :: mx_char, my_char, mz_char, local_char
   real(rp) :: axis(3)
   real(rp), allocatable :: proj_irr(:, :, :, :, :), eig_irr(:, :)
   real(rp) :: proj(4)
   integer, allocatable :: op(:)
   integer, allocatable :: star_size(:)
   integer :: nk_full, e_start, e_end, e_count

   call root_info('project_dos_orbitals_gaussian: Starting projection', __FILE__, __LINE__)

   ! Determine sigma (same as in calculate_dos_gaussian, already in Ry)
   if (this%gaussian_sigma < 0.001_rp) then
      sigma_use = this%calculate_adaptive_sigma()
   else
      sigma_use = this%gaussian_sigma
   end if
   sigma_squared = sigma_use**2

   this%n_sites = this%lattice%nrec
   this%n_orb_types = 4
   this%n_spin_components = 2
   nbands = size(this%eigenvalues, 1)

   if (this%n_sites <= 0) then
      call g_logger%fatal('project_dos_orbitals_gaussian: invalid number of sites', __FILE__, __LINE__)
   end if
   if (size(this%eigenvectors, 1) <= 0) then
      call g_logger%fatal('project_dos_orbitals_gaussian: eigenvectors not available', __FILE__, __LINE__)
   end if
   if (mod(size(this%eigenvectors, 1), this%n_sites) /= 0) then
      call g_logger%fatal('project_dos_orbitals_gaussian: eigenvector basis size not divisible by number of sites', __FILE__, __LINE__)
   end if
   n_orb_site = size(this%eigenvectors, 1)/this%n_sites
   if (mod(n_orb_site, this%n_spin_components) /= 0) then
      call g_logger%fatal('project_dos_orbitals_gaussian: per-site basis size incompatible with spin components', __FILE__, __LINE__)
   end if

   if (allocated(this%projected_dos)) deallocate(this%projected_dos)
   allocate(this%projected_dos(this%n_sites, this%n_orb_types, this%n_spin_components, this%n_energy_points))
   this%projected_dos = 0.0_rp

   if (allocated(this%dos_mx_tot)) deallocate(this%dos_mx_tot)
   if (allocated(this%dos_my_tot)) deallocate(this%dos_my_tot)
   if (allocated(this%dos_mz_tot)) deallocate(this%dos_mz_tot)
   allocate(this%dos_mx_tot(this%n_energy_points))
   allocate(this%dos_my_tot(this%n_energy_points))
   allocate(this%dos_mz_tot(this%n_energy_points))
   this%dos_mx_tot = 0.0_rp
   this%dos_my_tot = 0.0_rp
   this%dos_mz_tot = 0.0_rp

   if (allocated(this%projected_dos_moments)) deallocate(this%projected_dos_moments)
   allocate(this%projected_dos_moments(this%n_sites, this%n_orb_types, this%n_spin_components, 3, this%n_energy_points))
   this%projected_dos_moments = 0.0_rp

   ! Diagnostic logging
   call root_info('project_dos_orbitals_gaussian: n_sites = ' // trim(int2str(this%n_sites)) // &
                     ', nbands = ' // trim(int2str(nbands)) // &
                     ', max_orb_channels = ' // trim(int2str(this%max_orbs)) // &
                     ', eigenvector size = ' // trim(int2str(size(this%eigenvectors, 1))), __FILE__, __LINE__)

   ! Band projections of every irreducible k once; a reduced mesh is expanded by symmetry to the full mesh, where each
   ! full-mesh point carries the weight of its irreducible point shared over its star.
   call prepare_projections(this, n_orb_site, proj_irr, eig_irr, op)
   nk_full = product(this%nk_mesh)
   allocate(star_size(this%nk_total))
   star_size = 0
   do ik = 1, size(this%full_to_irred_k)
      star_size(this%full_to_irred_k(ik)) = star_size(this%full_to_irred_k(ik)) + 1
   end do
   ! Ranks split the energy grid, so every sum runs over the full mesh in the same order whatever the rank count.
   call get_mpi_range(rank, this%n_energy_points, e_start, e_end, e_count, region_tag='projected gaussian')

   do ie = e_start, e_end
      energy = this%dos_energy_grid(ie)  ! Already in Ry

      do ik = 1, nk_full
         ik_global = this%full_to_irred_k(ik)
         do ib = 1, nbands  ! Loop over all bands, not just max_orb_channels
            ! Skip if eigenvalue is far from current energy
            if (abs(eig_irr(ib, ik_global) - energy) > 5.0_rp * sigma_use) cycle

            ! Calculate Gaussian weight (already normalized)
            gaussian_weight = exp(-((energy - eig_irr(ib, ik_global))**2) / (2.0_rp * sigma_squared))
            gaussian_weight = gaussian_weight / (sigma_use * sqrt(2.0_rp * 3.141592653589793_rp))

            if (abs(gaussian_weight) < 1.0e-10_rp) cycle

            ! Weight of this full-mesh point: its irreducible weight shared over the star
            weight = gaussian_weight * this%k_weights(ik_global) / real(star_size(ik_global), rp)

            do isite = 1, this%n_sites
               axis = this%dos_site_axes(:, isite)
               do iorb = 1, 4
                  proj = this%full_mesh_projection(proj_irr, op, ik, iorb, isite, ib)
                  orbital_char = proj(1)
                  mx_char = proj(2)
                  my_char = proj(3)
                  mz_char = proj(4)
                  if (abs(orbital_char) < 1.0e-14_rp) cycle

                  local_char = axis(1)*mx_char + axis(2)*my_char + axis(3)*mz_char
                  this%projected_dos(isite, iorb, 1, ie) = this%projected_dos(isite, iorb, 1, ie) + &
                                                            0.5_rp*(orbital_char + local_char)*weight
                  this%projected_dos(isite, iorb, 2, ie) = this%projected_dos(isite, iorb, 2, ie) + &
                                                            0.5_rp*(orbital_char - local_char)*weight
                  this%projected_dos_moments(isite, iorb, 1, 1, ie) = this%projected_dos_moments(isite, iorb, 1, 1, ie) + mx_char*weight
                  this%projected_dos_moments(isite, iorb, 1, 2, ie) = this%projected_dos_moments(isite, iorb, 1, 2, ie) + my_char*weight
                  this%projected_dos_moments(isite, iorb, 1, 3, ie) = this%projected_dos_moments(isite, iorb, 1, 3, ie) + mz_char*weight
                  this%dos_mx_tot(ie) = this%dos_mx_tot(ie) + mx_char*weight
                  this%dos_my_tot(ie) = this%dos_my_tot(ie) + my_char*weight
                  this%dos_mz_tot(ie) = this%dos_mz_tot(ie) + mz_char*weight
               end do
            end do
         end do
      end do
   end do
   deallocate(proj_irr, eig_irr, op, star_size)

#ifdef USE_MPI
   if (numprocs > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%projected_dos, product(shape(this%projected_dos)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%projected_dos_moments, product(shape(this%projected_dos_moments)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%dos_mx_tot, this%n_energy_points, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%dos_my_tot, this%n_energy_points, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%dos_mz_tot, this%n_energy_points, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
   end if
#endif

   call this%symmetrize_projected_dos()
end subroutine project_dos_orbitals_gaussian

   !> @brief Return the local spin axis for a projected-DOS site.
   !> @param[in] this Reciprocal object containing lattice magnetic moments.
   !> @param[in] isite Site index in the reciprocal packing.
   !> @param[out] axis Normalized local spin axis.
   module subroutine get_site_spin_axis(this, isite, axis)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: isite
      real(rp), intent(out) :: axis(3)
      integer :: atom_idx
      real(rp) :: axis_norm

      axis = [0.0_rp, 0.0_rp, 1.0_rp]
      atom_idx = this%lattice%nbulk + isite
      ! gbt_single_q: the eigenvector spin basis is the site's local frame, whose majority axis is +z whatever potential%mom says
      if (trim(this%hamiltonian%magnetic_representation) == gbt_single_q) atom_idx = 0
      if (atom_idx >= 1 .and. atom_idx <= size(this%lattice%symbolic_atoms)) then
         if (sqrt(sum(this%lattice%symbolic_atoms(atom_idx)%constraint_target(:)**2)) > 1.0e-12_rp) then
            axis = this%lattice%symbolic_atoms(atom_idx)%constraint_target(:)
         else
            axis = this%lattice%symbolic_atoms(atom_idx)%potential%mom(:)
         end if
      end if

      axis_norm = sqrt(sum(axis**2))
      if (axis_norm > tiny(1.0_rp)) then
         axis = axis/axis_norm
      else
         axis = [0.0_rp, 0.0_rp, 1.0_rp]
      end if
   end subroutine get_site_spin_axis

   !> Occupation of one band: its Bloechl tetrahedron weight when occ_t is allocated (dos_method = 'blochl'), else Fermi-Dirac.
   pure function band_occupation(occ_t, ib, ik_global, energy, fermi_level, kT) result(occ)
      real(rp), allocatable, intent(in) :: occ_t(:, :)
      integer, intent(in) :: ib, ik_global
      real(rp), intent(in) :: energy, fermi_level, kT
      real(rp) :: occ

      if (allocated(occ_t)) then
         occ = occ_t(ib, ik_global)
      else
         occ = fermi_dirac_occupation(energy, fermi_level, kT)
      end if
   end function band_occupation

   !> @brief Project one spinor eigenvector block onto scalar and vector weights.
   !> @param[in] this Reciprocal object containing eigenvectors.
   !> @param[in] ik k-point index.
   !> @param[in] ib Band index.
   !> @param[in] site_orb_start First orbital index for the site block.
   !> @param[in] n_orb_per_spin Number of orbital channels per spin.
   !> @param[in] i_first First local orbital index in the projection range.
   !> @param[in] i_last Last local orbital index in the projection range.
   !> @param[out] total_char Scalar projected character.
   !> @param[out] mx_char Projected x spin moment character.
   !> @param[out] my_char Projected y spin moment character.
   !> @param[out] mz_char Projected z spin moment character.
   module subroutine compute_spinor_block_projection(this, ik, ib, site_orb_start, n_orb_per_spin, i_first, i_last, &
                                              total_char, mx_char, my_char, mz_char)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: ik, ib, site_orb_start, n_orb_per_spin, i_first, i_last
      real(rp), intent(out) :: total_char, mx_char, my_char, mz_char

      integer :: i, idx_up, idx_dn
      real(rp) :: up_w, dn_w
      complex(rp) :: u, d, ud

      total_char = 0.0_rp
      mx_char = 0.0_rp
      my_char = 0.0_rp
      mz_char = 0.0_rp

      do i = i_first, min(i_last, n_orb_per_spin)
         idx_up = site_orb_start + i
         idx_dn = site_orb_start + n_orb_per_spin + i
         if (idx_dn > size(this%eigenvectors, 1)) cycle

         u = this%eigenvectors(idx_up, ib, ik)
         d = this%eigenvectors(idx_dn, ib, ik)
         up_w = real(conjg(u)*u, rp)
         dn_w = real(conjg(d)*d, rp)
         ud = conjg(u)*d

         total_char = total_char + up_w + dn_w
         mx_char = mx_char + 2.0_rp*real(ud, rp)
         my_char = my_char + 2.0_rp*aimag(ud)
         mz_char = mz_char + up_w - dn_w
      end do
   end subroutine compute_spinor_block_projection

   !> @brief Band projections of every irreducible k, the eigenvalues and the operation of each full-mesh point.
   !> @details Shared setup of the Gaussian and tetrahedron projected DOS. Projections are averaged over degenerate groups and
   !>          gathered on every rank; op is empty when the mesh is not reduced.
   !> @param[in] n_orb_site Eigenvector components of one site.
   !> @param[out] proj_irr Projection characters (4, 4, site, band, irreducible k).
   !> @param[out] eig_irr Eigenvalues (band, irreducible k).
   !> @param[out] op Operation index of each full-mesh point, see full_mesh_operations.
   subroutine prepare_projections(this, n_orb_site, proj_irr, eig_irr, op)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: n_orb_site
      real(rp), allocatable, intent(out) :: proj_irr(:, :, :, :, :), eig_irr(:, :)
      integer, allocatable, intent(out) :: op(:)
      real(rp), allocatable :: proj_local(:, :, :, :, :)
      real(rp) :: orbital_char
      integer :: ik, ib, isite, iorb, nbands, n_orb_per_spin

      nbands = size(this%eigenvalues, 1)
      n_orb_per_spin = n_orb_site/this%n_spin_components
      allocate(proj_local(4, 4, this%n_sites, nbands, size(this%eigenvectors, 3)))
      proj_local = 0.0_rp
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(this, proj_local, n_orb_site, n_orb_per_spin, nbands) &
      !$OMP PRIVATE(ik, ib, isite, iorb, orbital_char) SCHEDULE(DYNAMIC)
      do ik = 1, size(this%eigenvectors, 3)
         do ib = 1, nbands
            do isite = 1, this%n_sites
               do iorb = 1, 4
                  if (lstart(iorb) > n_orb_per_spin) cycle
                  call compute_spinor_block_projection(this, ik, ib, (isite - 1)*n_orb_site, n_orb_per_spin, &
                                                       lstart(iorb), lend(iorb), orbital_char, proj_local(2, iorb, isite, ib, ik), &
                                                       proj_local(3, iorb, isite, ib, ik), proj_local(4, iorb, isite, ib, ik))
                  proj_local(1, iorb, isite, ib, ik) = orbital_char
               end do
            end do
         end do
      end do
      !$OMP END PARALLEL DO
      call average_degenerate_projections(this%eigenvalues, proj_local)
      allocate(proj_irr(4, 4, this%n_sites, nbands, this%k_workset%nk_global), eig_irr(nbands, this%k_workset%nk_global))
      call this%gather_k_columns(16*this%n_sites*nbands, proj_local, proj_irr)
      call this%gather_k_columns(nbands, this%eigenvalues, eig_irr)
      ! A reduced mesh is expanded by symmetry to the full mesh, not rediagonalised.
      if (this%nk_total < product(this%nk_mesh)) then
         call this%full_mesh_operations(op)
      else
         allocate(op(0))
      end if
   end subroutine prepare_projections

   !> @brief Replace each band's projection characters by their average over its degenerate group at the same k.
   !> @details The characters of a single band inside a degenerate group depend on the arbitrary basis the eigensolver
   !>          returns; the group sum does not, so the group average is the basis-invariant per-band value.
   !>          Bands belong to one group when consecutive eigenvalues differ by less than DEGENERACY_TOL (Ry).
   !> @param[in] eig Eigenvalues (band, k), ascending in band.
   !> @param[inout] proj Projection characters (4, 4, site, band, k).
   subroutine average_degenerate_projections(eig, proj)
      real(rp), intent(in) :: eig(:, :)
      real(rp), intent(inout) :: proj(:, :, :, :, :)
      real(rp), parameter :: DEGENERACY_TOL = 1.0e-8_rp
      integer :: ik, ib, ie

      do ik = 1, size(eig, 2)
         ib = 1
         do while (ib <= size(eig, 1))
            ie = ib
            do while (ie < size(eig, 1))
               if (abs(eig(ie + 1, ik) - eig(ie, ik)) >= DEGENERACY_TOL) exit
               ie = ie + 1
            end do
            if (ie > ib) proj(:, :, :, ib:ie, ik) = spread(sum(proj(:, :, :, ib:ie, ik), dim=4)/real(ie - ib + 1, rp), 4, ie - ib + 1)
            ib = ie + 1
         end do
      end do
   end subroutine average_degenerate_projections

   !> @brief Calculate projected DOS with tetrahedron integration.
   !> @details Averages orbital/spinor projection weights over tetrahedron
   !>          corners and accumulates scalar plus directional projected DOS.
   !> @param[inout] this Reciprocal object receiving projected DOS arrays.
   module subroutine project_dos_orbitals_tetrahedron(this)
      class(reciprocal), intent(inout) :: this

      ! Local variables
      integer :: i_energy, i_tet, i_corner, i_band, iorb, isite
      integer :: i_start, i_end
      integer :: n_orb_per_spin, ik, n_orb_site
      integer :: nbands
      integer :: e_start, e_end, e_count
      real(rp) :: orbital_char_avg, de, mx_eff, my_eff, mz_eff
      real(rp), dimension(4) :: e_corners, orbital_chars
      real(rp), dimension(4) :: mx_chars, my_chars, mz_chars
      real(rp) :: tet_weight, corner_char
      logical :: use_blochl
      real(rp), allocatable :: corner_w(:, :)
      real(rp) :: local_char
      real(rp) :: axis(3)
      real(rp), allocatable :: local_projected_dos(:, :, :, :)
      real(rp), allocatable :: local_projected_dos_moments(:, :, :, :, :)
      real(rp), allocatable :: local_dos_mx(:), local_dos_my(:), local_dos_mz(:)
      ! Projection weights per k-point: (total, mx, my, mz; l; site; band; k)
      real(rp), allocatable :: proj_irr(:, :, :, :, :), eig_irr(:, :)
      real(rp) :: proj(4)
      integer, allocatable :: op(:)

      call root_info('project_dos_orbitals_tetrahedron: Starting tetrahedron orbital projection calculation', __FILE__, __LINE__)

      ! Initialize dimensions - get number of sites from lattice
      this%n_sites = this%lattice%nrec
      this%n_orb_types = 4  ! s, p, d, f
      this%n_spin_components = 2  ! spin up/down
      if (this%n_sites <= 0) then
         call g_logger%fatal('project_dos_orbitals_tetrahedron: invalid number of sites', __FILE__, __LINE__)
      end if
      if (size(this%eigenvectors, 1) <= 0) then
         call g_logger%fatal('project_dos_orbitals_tetrahedron: eigenvectors not available', __FILE__, __LINE__)
      end if
      if (mod(size(this%eigenvectors, 1), this%n_sites) /= 0) then
         call g_logger%fatal('project_dos_orbitals_tetrahedron: eigenvector basis size not divisible by number of sites', __FILE__, __LINE__)
      end if
      n_orb_site = size(this%eigenvectors, 1)/this%n_sites
      if (mod(n_orb_site, this%n_spin_components) /= 0) then
         call g_logger%fatal('project_dos_orbitals_tetrahedron: per-site basis size incompatible with spin components', __FILE__, __LINE__)
      end if
      n_orb_per_spin = n_orb_site/this%n_spin_components

      call root_info('project_dos_orbitals_tetrahedron: Projecting onto ' // trim(int2str(this%n_sites)) // &
                        ' site(s)', __FILE__, __LINE__)

      ! Allocate projected DOS array
      if (allocated(this%projected_dos)) deallocate(this%projected_dos)
      allocate(this%projected_dos(this%n_sites, this%n_orb_types, this%n_spin_components, this%n_energy_points))
      this%projected_dos = 0.0_rp

      if (allocated(this%dos_mx_tot)) deallocate(this%dos_mx_tot)
      if (allocated(this%dos_my_tot)) deallocate(this%dos_my_tot)
      if (allocated(this%dos_mz_tot)) deallocate(this%dos_mz_tot)
      allocate(this%dos_mx_tot(this%n_energy_points))
      allocate(this%dos_my_tot(this%n_energy_points))
      allocate(this%dos_mz_tot(this%n_energy_points))
      this%dos_mx_tot = 0.0_rp
      this%dos_my_tot = 0.0_rp
      this%dos_mz_tot = 0.0_rp

      if (allocated(this%projected_dos_moments)) deallocate(this%projected_dos_moments)
      allocate(this%projected_dos_moments(this%n_sites, this%n_orb_types, this%n_spin_components, 3, this%n_energy_points))
      this%projected_dos_moments = 0.0_rp

      use_blochl = trim(this%dos_method) == 'blochl'

      ! Setup tetrahedra if not already done
      if (.not. allocated(this%tetrahedra)) then
         call this%setup_tetrahedra()
      end if
      de = (this%dos_energy_grid(this%n_energy_points) - this%dos_energy_grid(1)) / &
           real(this%n_energy_points - 1, rp)
      if (de <= 0.0_rp) then
         call g_logger%fatal('project_dos_orbitals_tetrahedron: Invalid DOS energy grid spacing', __FILE__, __LINE__)
      end if
      ! Ranks split the energy grid, so every sum runs over all tetrahedra in the same order whatever the rank count.
      call get_mpi_range(rank, this%n_energy_points, e_start, e_end, e_count, region_tag='projected tetra')

      ! Each k-point is a corner of many tetrahedra: compute its projection weights once.
      nbands = size(this%eigenvalues, 1)
      call prepare_projections(this, n_orb_site, proj_irr, eig_irr, op)

      ! Tetra/band-first traversal with one tetra shape evaluation per band.
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(this, proj_irr, eig_irr, op, n_orb_per_spin, &
      !$OMP        e_start, e_end, de, use_blochl) &
      !$OMP PRIVATE(i_tet, i_band, i_corner, ik, proj, e_corners, i_start, i_end, &
      !$OMP         i_energy, tet_weight, isite, axis, iorb, mx_eff, my_eff, mz_eff, &
      !$OMP         orbital_char_avg, orbital_chars, &
      !$OMP         mx_chars, my_chars, mz_chars, local_char, local_projected_dos, &
      !$OMP         local_projected_dos_moments, local_dos_mx, local_dos_my, local_dos_mz, &
      !$OMP         corner_w, corner_char)
      allocate(corner_w(4, this%n_energy_points))
      allocate(local_projected_dos(this%n_sites, this%n_orb_types, this%n_spin_components, this%n_energy_points))
      allocate(local_projected_dos_moments(this%n_sites, this%n_orb_types, this%n_spin_components, 3, this%n_energy_points))
      allocate(local_dos_mx(this%n_energy_points))
      allocate(local_dos_my(this%n_energy_points))
      allocate(local_dos_mz(this%n_energy_points))
      local_projected_dos = 0.0_rp
      local_projected_dos_moments = 0.0_rp
      local_dos_mx = 0.0_rp
      local_dos_my = 0.0_rp
      local_dos_mz = 0.0_rp

      !$OMP DO SCHEDULE(DYNAMIC)
      do i_tet = 1, this%n_tetrahedra
         do i_band = 1, size(this%eigenvalues, 1)
            do i_corner = 1, 4
               ik = this%tetrahedra(i_corner, i_tet)
               e_corners(i_corner) = eig_irr(i_band, this%full_to_irred_k(ik))
            end do
            tet_weight = this%tetrahedron_volumes(i_tet)
            ! Bin-averaged corner weights on this rank's energy window only
            call tetra_bin_weights(e_corners, use_blochl, this%dos_energy_grid(1), de, e_start, e_end, i_start, i_end, corner_w)
            if (i_start > i_end) cycle

            do isite = 1, this%n_sites
               axis = this%dos_site_axes(:, isite)
               do iorb = 1, 4
                  orbital_char_avg = 0.0_rp
                  if (lstart(iorb) <= n_orb_per_spin) then
                     do i_corner = 1, 4
                        ik = this%tetrahedra(i_corner, i_tet)
                        proj = this%full_mesh_projection(proj_irr, op, ik, iorb, isite, i_band)
                        orbital_chars(i_corner) = proj(1)
                        mx_chars(i_corner) = proj(2)
                        my_chars(i_corner) = proj(3)
                        mz_chars(i_corner) = proj(4)
                     end do
                     orbital_char_avg = sum(orbital_chars) / 4.0_rp
                  end if
                  if (abs(orbital_char_avg) < 1.0e-14_rp) cycle

                  ! Each corner's own weight multiplies its own character (n/4 per corner without the Bloechl correction);
                  ! the corner sum is formed first so each accumulator is updated once per energy.
                  do i_energy = i_start, i_end
                     corner_char = tet_weight*dot_product(corner_w(:, i_energy), orbital_chars)
                     mx_eff = tet_weight*dot_product(corner_w(:, i_energy), mx_chars)
                     my_eff = tet_weight*dot_product(corner_w(:, i_energy), my_chars)
                     mz_eff = tet_weight*dot_product(corner_w(:, i_energy), mz_chars)
                     local_char = axis(1)*mx_eff + axis(2)*my_eff + axis(3)*mz_eff
                     local_projected_dos(isite, iorb, 1, i_energy) = local_projected_dos(isite, iorb, 1, i_energy) + &
                                                                      0.5_rp*(corner_char + local_char)
                     local_projected_dos(isite, iorb, 2, i_energy) = local_projected_dos(isite, iorb, 2, i_energy) + &
                                                                      0.5_rp*(corner_char - local_char)
                     local_projected_dos_moments(isite, iorb, 1, 1, i_energy) = &
                        local_projected_dos_moments(isite, iorb, 1, 1, i_energy) + mx_eff
                     local_projected_dos_moments(isite, iorb, 1, 2, i_energy) = &
                        local_projected_dos_moments(isite, iorb, 1, 2, i_energy) + my_eff
                     local_projected_dos_moments(isite, iorb, 1, 3, i_energy) = &
                        local_projected_dos_moments(isite, iorb, 1, 3, i_energy) + mz_eff
                     local_dos_mx(i_energy) = local_dos_mx(i_energy) + mx_eff
                     local_dos_my(i_energy) = local_dos_my(i_energy) + my_eff
                     local_dos_mz(i_energy) = local_dos_mz(i_energy) + mz_eff
                  end do
               end do
            end do
         end do
      end do
      !$OMP END DO

      !$OMP CRITICAL(projected_tetra_accum)
      this%projected_dos = this%projected_dos + local_projected_dos
      this%projected_dos_moments = this%projected_dos_moments + local_projected_dos_moments
      this%dos_mx_tot = this%dos_mx_tot + local_dos_mx
      this%dos_my_tot = this%dos_my_tot + local_dos_my
      this%dos_mz_tot = this%dos_mz_tot + local_dos_mz
      !$OMP END CRITICAL(projected_tetra_accum)

      deallocate(corner_w)
      deallocate(local_projected_dos)
      deallocate(local_projected_dos_moments)
      deallocate(local_dos_mx)
      deallocate(local_dos_my)
      deallocate(local_dos_mz)
      !$OMP END PARALLEL

#ifdef USE_MPI
      if (numprocs > 1) then
         call MPI_ALLREDUCE(MPI_IN_PLACE, this%projected_dos, product(shape(this%projected_dos)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, this%projected_dos_moments, product(shape(this%projected_dos_moments)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, this%dos_mx_tot, this%n_energy_points, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, this%dos_my_tot, this%n_energy_points, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, this%dos_mz_tot, this%n_energy_points, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      end if
#endif

      deallocate(proj_irr, eig_irr, op)

      call this%symmetrize_projected_dos()

      call root_info('project_dos_orbitals_tetrahedron: Tetrahedron orbital projection calculation completed', __FILE__, __LINE__)
   end subroutine project_dos_orbitals_tetrahedron

   !> @brief Integrate tabulated data with the trapezoidal rule.
   !> @param[in] x Grid coordinates.
   !> @param[in] y Values on the grid.
   !> @return Trapezoidal integral over the grid.
   module function trapezoidal_integral(x, y) result(integral)
      real(rp), dimension(:), intent(in) :: x, y
      real(rp) :: integral

      ! Local variables
      integer :: n, i
      real(rp) :: dx

      n = size(x)
      if (n /= size(y)) then
         call g_logger%error('trapezoidal_integral: x and y arrays must have same size', __FILE__, __LINE__)
         integral = 0.0_rp
         return
      end if

      integral = 0.0_rp

      ! Trapezoidal rule: ∫ y dx ≈ Σ (y_i + y_{i+1}) * (x_{i+1} - x_i) / 2
      do i = 1, n-1
         dx = x(i+1) - x(i)
         integral = integral + 0.5_rp * (y(i) + y(i+1)) * dx
      end do
   end function trapezoidal_integral

   !> @brief Linearly interpolate a tabulated grid value.
   !> @param[in] x Grid coordinates.
   !> @param[in] y Values on the grid.
   !> @param[in] n Number of grid points.
   !> @param[in] x0 Coordinate where y is requested.
   !> @return Interpolated value at x0.
   module function interpolate_grid_value(x, y, n, x0) result(y0)
      integer, intent(in) :: n
      real(rp), intent(in) :: x(n), y(n), x0
      real(rp) :: y0
      integer :: i
      real(rp) :: dx

      if (n <= 0) then
         y0 = 0.0_rp
         return
      end if
      if (x0 <= x(1)) then
         y0 = y(1)
         return
      end if
      if (x0 >= x(n)) then
         y0 = y(n)
         return
      end if

      do i = 1, n - 1
         if (x(i + 1) >= x0) then
            dx = x(i + 1) - x(i)
            if (abs(dx) > tiny(1.0_rp)) then
               y0 = y(i) + (y(i + 1) - y(i)) * (x0 - x(i)) / dx
            else
               y0 = 0.5_rp * (y(i) + y(i + 1))
            end if
            return
         end if
      end do

      y0 = y(n)
   end function interpolate_grid_value

!> @brief Calculate zeroth, first, and second DOS band moments.
!> @details Integrates projected DOS channels with Fermi weighting to produce
!>          occupation, first-energy, and second-energy moments.
!> @param[inout] this Reciprocal object receiving band_moments.
module subroutine calculate_band_moments(this)
   class(reciprocal), intent(inout) :: this

   integer :: isite, iorb, ispin
   real(rp) :: m0
   real(rp), allocatable :: reference(:, :), axis_used(:, :)

   call root_info('calculate_band_moments: Starting calculation', __FILE__, __LINE__)
   this%n_sites = this%lattice%nrec

   if (this%total_electrons <= 1.0e-3_rp) then
      if (this%lattice%nbulk_bulk > 0) then
         this%total_electrons = real(sum(this%lattice%symbolic_atoms(1:this%lattice%nbulk_bulk)%element%valence), rp)
         call root_info('calculate_band_moments: Recovered total_electrons from valence = ' // &
                           trim(real2str(this%total_electrons, '(F10.5)')), __FILE__, __LINE__)
      else if (this%lattice%nrec > 0) then
         this%total_electrons = real(sum(this%lattice%symbolic_atoms(1:this%lattice%nrec)%element%valence), rp)
         call g_logger%warning('calculate_band_moments: nbulk_bulk<=0; recovered total_electrons using nrec span.', __FILE__, __LINE__)
      end if
   end if

   ! EF is already canonical here: calculate_density_of_states solves it from
   ! eigenvalues/k weights before any DOS projection.  This routine must never
   ! replace it with a grid/window-dependent DOS root.
   call root_info('calculate_band_moments: Using canonical/fixed Fermi level = ' // &
                      trim(real2str(this%fermi_level, '(F 8.5)')) // ' Ry', __FILE__, __LINE__)

   ! --- WP7: the SCF density is the shared rotating-frame contract. ----------
   ! Accumulate the full 2x2 spin density from eigenvectors / k weights /
   ! occupations with no axis in it, resolve the explicit per-site axis through
   ! the SCF policy, and only then project onto the radial up/down band moments.
   ! `projected_dos` above stays a DOS diagnostic; it no longer defines the
   ! density, and its `potential%mom` projection is no longer an SCF input.
   allocate(reference(3, this%n_sites))
   allocate(axis_used(3, this%n_sites))
   do isite = 1, this%n_sites
      call get_site_spin_axis(this, isite, reference(:, isite))
   end do

   call this%accumulate_spin_density_kspace()
   call this%fill_band_moments_from_spin_density( &
      trim(this%lattice%control%density_policy), reference, axis_used)

   deallocate(reference, axis_used)

   ! Calculate sum of all m0 moments (should equal total occupation)
   m0 = 0.0_rp
   do isite = 1, this%n_sites
      do iorb = 1, this%n_orb_types
         do ispin = 1, this%n_spin_components
            m0 = m0 + this%band_moments(isite, iorb, ispin, 1)
         end do
      end do
   end do

   call root_info('calculate_band_moments: Sum of all m0 moments = ' // &
                     trim(real2str(m0, '(F10.5)')) // &
                     ', expected electrons = ' // trim(real2str(this%total_electrons, '(F10.5)')), &
                     __FILE__, __LINE__)

   call root_info('calculate_band_moments: Completed', __FILE__, __LINE__)
end subroutine calculate_band_moments

   !> @brief Reconstruct lattice local density matrices from projected DOS moments.
   !> @param[inout] this Reciprocal object containing projected DOS/moments.
   !> @param[inout] lattice_obj Lattice object receiving local density matrices.
   module subroutine calculate_ldm_from_projected_dos(this, lattice_obj, fixed_liechtenstein, m1_out, m1_only)
      use lattice_mod
      class(reciprocal), intent(inout) :: this
      type(lattice), intent(inout) :: lattice_obj
      logical, intent(in), optional :: fixed_liechtenstein
      complex(rp), intent(out), optional :: m1_out(:, :, :, :, :)
      logical, intent(in), optional :: m1_only
      
      ! Local variables
      integer :: isite
      logical :: operator_changed, m1_only_use
      
      ! Orbital indexing: 
      ! iorb = 1 (s, 1 orbital), iorb = 2 (p, 3 orbitals), iorb = 3 (d, 5 orbitals)
      
      call root_info('calculate_ldm_from_projected_dos: Computing LDM for LDA+U from k-space DOS', __FILE__, __LINE__)

      call this%invalidate_if_operator_changed('reciprocal%calculate_ldm_from_projected_dos', operator_changed)
      if (operator_changed .or. .not. allocated(this%eigenvalues)) then
         if (.not. allocated(this%k_points)) then
            if (this%use_symmetry_reduction) then
               call this%generate_reduced_kpoint_mesh(this%nk_mesh, sum(abs(this%k_offset)) > 1.0e-12_rp)
            else
               call this%generate_mp_mesh()
            end if
         end if
         call this%build_kspace_hamiltonian()
         call this%diagonalize_hamiltonian()
         this%fermi_level = this%find_fermi_level_from_eigenvalues(this%total_electrons)
      end if
      
      m1_only_use = .false.
      if (present(m1_only)) m1_only_use = m1_only
      ! Initialize all LDM to zero unless this is a diagnostic-only M1 call.
      if (.not. m1_only_use) then
         do isite = 1, lattice_obj%nrec
            lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%ldm(:,:,:,:) = 0.0_rp
         end do
      end if
      
      ! Method 2: Full LDM (including off-diagonal) from eigenvector projections
      ! This properly captures orbital hybridization and is more accurate than diagonal-only
      call calculate_ldm_from_eigenvectors(this, lattice_obj, fixed_liechtenstein, m1_out, m1_only_use)
      
      ! Save flattened copy for each site
      if (.not. m1_only_use) then
         do isite = 1, lattice_obj%nrec
            call lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%flatten_ldm()
         end do
      end if
      
      ! Print LDM for debugging
      call root_info('calculate_ldm_from_projected_dos: LDM calculation complete', __FILE__, __LINE__)
      
   end subroutine calculate_ldm_from_projected_dos

   !> @brief Reconstruct lattice local density matrices directly from eigenvectors.
   !> @details Accumulates occupied spinor density matrices over k-points and
   !>          bands, then synchronizes the lattice result across MPI ranks.
   !> @param[inout] this Reciprocal object containing eigenvectors and occupations.
   !> @param[inout] lattice_obj Lattice object receiving local density matrices.
   module subroutine calculate_ldm_from_eigenvectors(this, lattice_obj, fixed_liechtenstein, m1_out, m1_only)
      use lattice_mod
      class(reciprocal), intent(inout) :: this
      type(lattice), intent(inout) :: lattice_obj
      logical, intent(in), optional :: fixed_liechtenstein
      complex(rp), intent(out), optional :: m1_out(:, :, :, :, :)
      logical, intent(in), optional :: m1_only
      
      ! Local variables
      integer :: ik, ik_global, ib, isite, alpha, beta, iorb, norb_l
      integer :: ispin, basis_offset, site_offset, n_orb_site, n_orb_per_spin
      real(rp) :: kweight, fermi_occ
      complex(rp), allocatable :: DM(:,:)
      complex(rp), allocatable :: DM_m1(:,:)
      complex(rp), allocatable :: rho(:, :, :), rho_m1(:, :, :)
      complex(rp) :: psi_a, psi_b
      integer :: nbands, nsites, i
      integer :: m1_lmax_required, site_lmax
      real(rp) :: trace_tot, trace_up, trace_dn, max_imag, max_herm
      integer :: idx_a, idx_b
      logical :: use_complex_projector
      logical :: m1_only_use
      real(rp) :: fermi_occ_m1
      logical :: blochl
      real(rp), allocatable :: occ_t(:, :)
      real(rp) :: kT
      real(rp) :: ldm_trace, ldm_min_eig, ldm_max_eig, ldm_antiherm
      integer :: ldm_diag_info
      ! The RS runtime observes anti-Hermitian residuals near 1e-8 or below;
      ! these looser thresholds allow integration/k-sum roundoff while still
      ! flagging a materially invalid occupation. No value is projected.
      real(rp), parameter :: ldm_hermiticity_tol = 1.0e-6_rp
      real(rp), parameter :: ldm_occupation_tol = 1.0e-3_rp

      call root_info('calculate_ldm_from_eigenvectors: Computing site density matrices from eigenvectors', __FILE__, __LINE__)

      use_complex_projector = .false.
      if (present(fixed_liechtenstein)) use_complex_projector = fixed_liechtenstein
      m1_only_use = .false.
      if (present(m1_only)) m1_only_use = m1_only
      kT = max(this%temperature*kB_Ry_per_K, 1.0e-10_rp)

      nbands = size(this%eigenvalues, 1)
      nsites = lattice_obj%nrec
      if (nsites <= 0) then
         call g_logger%fatal('calculate_ldm_from_eigenvectors: invalid number of sites', __FILE__, __LINE__)
      end if
      if (mod(size(this%eigenvectors, 1), nsites) /= 0) then
         call g_logger%fatal('calculate_ldm_from_eigenvectors: eigenvector basis size not divisible by number of sites', __FILE__, __LINE__)
      end if
      n_orb_site = size(this%eigenvectors, 1)/nsites
      if (mod(n_orb_site, this%n_spin_components) /= 0) then
         call g_logger%fatal('calculate_ldm_from_eigenvectors: per-site basis size incompatible with spin components', __FILE__, __LINE__)
      end if
      n_orb_per_spin = n_orb_site/this%n_spin_components

      if (present(m1_out)) then
         m1_lmax_required = 0
         do isite = 1, nsites
            site_lmax = min(3, lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%lmax)
            do iorb = 0, site_lmax
               if (iorb*iorb < n_orb_per_spin) m1_lmax_required = max(m1_lmax_required, iorb)
            end do
         end do
         if (size(m1_out, 1) /= nsites .or. size(m1_out, 2) < m1_lmax_required + 1 .or. &
             size(m1_out, 3) < this%n_spin_components .or. &
             size(m1_out, 4) < 2*m1_lmax_required + 1 .or. size(m1_out, 5) < 2*m1_lmax_required + 1) then
            call g_logger%fatal('calculate_ldm_from_eigenvectors: m1_out has incompatible extents', __FILE__, __LINE__)
         end if
         m1_out = cmplx(0.0_rp, 0.0_rp, rp)
         allocate(DM_m1(n_orb_site, n_orb_site), rho_m1(n_orb_site, n_orb_site, nsites))
         rho_m1 = (0.0_rp, 0.0_rp)
      end if

      ! Allocate temporary DM (n_orb_site x n_orb_site)
      allocate(DM(n_orb_site, n_orb_site), rho(n_orb_site, n_orb_site, nsites))
      rho = (0.0_rp, 0.0_rp)

      ! Zero per-site LDM storage before accumulation
      if (.not. m1_only_use) then
         do isite = 1, nsites
            lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%ldm(:,:,:,:) = 0.0_rp
         end do
      end if

      blochl = trim(this%dos_method) == 'blochl'
      if (blochl) call this%blochl_occupations(occ_t)   ! collective: before the rank chain below
      if (this%k_mesh_distributed_active) then
         if (present(m1_out)) then
            call chain_complex(rho, size(rho), .true., rho_m1, size(rho_m1))
         else
            call chain_complex(rho, size(rho), .true.)
         end if
      end if

      ! Accumulate density matrix per site
      ! Loop order: k-points -> sites -> bands for better cache locality
      do ik = 1, size(this%eigenvalues, 2)
         ik_global = local_k_index_to_global(this, ik)
         kweight = this%k_weights(ik_global)
         
         do isite = 1, nsites
            ! Reset DM for this site at this k-point
            DM = (0.0_rp, 0.0_rp)
            if (present(m1_out)) DM_m1 = (0.0_rp, 0.0_rp)
            site_offset = (isite - 1) * n_orb_site
            
            ! Accumulate contributions from all occupied bands at this k-point
            do ib = 1, nbands
               if (blochl .or. use_complex_projector .or. m1_only_use) then
                  fermi_occ = band_occupation(occ_t, ib, ik_global, this%eigenvalues(ib, ik), this%fermi_level, kT)
               else
                  if (this%eigenvalues(ib, ik) > this%fermi_level) cycle
                  fermi_occ = 1.0_rp
               end if
               
               ! Fixed Liechtenstein stores the conventional projector P_ab=ψ_a ψ*_b.
               ! Legacy SC-U/V retains its historical transposed real projector.
               do alpha = 1, n_orb_site
                  psi_a = this%eigenvectors(site_offset + alpha, ib, ik)
                  do beta = 1, n_orb_site
                     psi_b = this%eigenvectors(site_offset + beta, ib, ik)
                     if (present(m1_out)) then
                        fermi_occ_m1 = band_occupation(occ_t, ib, ik_global, this%eigenvalues(ib, ik), this%fermi_level, kT)
                        DM_m1(alpha, beta) = DM_m1(alpha, beta) + &
                           hubbard_projector_element(psi_a, psi_b) * fermi_occ_m1 * kweight * this%eigenvalues(ib, ik)
                     end if
                     if (use_complex_projector) then
                        DM(alpha, beta) = DM(alpha, beta) + &
                           hubbard_projector_element(psi_a, psi_b) * fermi_occ * kweight
                     else
                        DM(alpha, beta) = DM(alpha, beta) + conjg(psi_a) * psi_b * fermi_occ * kweight
                     end if
                  end do
               end do
            end do ! bands
            
            ! Diagnostics: check DM properties for first k-point
            if (ik_global == 1) then
               trace_tot = 0.0_rp
               trace_up = 0.0_rp
               trace_dn = 0.0_rp
               max_imag = 0.0_rp
               max_herm = maxval(abs(DM - transpose(conjg(DM))))
               
               do i = 1, n_orb_site
                  trace_tot = trace_tot + real(DM(i, i), rp)
                  ! Check for large imaginary parts (should be ~0 for physical DM)
                  max_imag = max(max_imag, abs(aimag(DM(i, i))))
                  
                  if (i <= n_orb_per_spin) then
                     trace_up = trace_up + real(DM(i, i), rp)
                  else
                     trace_dn = trace_dn + real(DM(i, i), rp)
                  end if
               end do
               
               call root_info('calculate_ldm_from_eigenvectors: DM diagnostics site=' // trim(int2str(isite)) // &
                                 ', trace_tot=' // trim(real2str(trace_tot, '(F10.6)')) // &
                                 ', trace_up=' // trim(real2str(trace_up, '(F10.6)')) // &
                                 ', trace_dn=' // trim(real2str(trace_dn, '(F10.6)')) // &
                                 ', max_imag=' // trim(real2str(max_imag, '(E10.2)')) // &
                                 ', max_antiherm=' // trim(real2str(max_herm, '(E10.2)')), __FILE__, __LINE__)
               
               if (max_imag > 1.0e-6_rp) then
                  call g_logger%warning('calculate_ldm_from_eigenvectors: Large imaginary part in DM diagonal: ' // &
                                       trim(real2str(max_imag, '(E10.2)')), __FILE__, __LINE__)
               end if
            end if
            
            rho(:, :, isite) = rho(:, :, isite) + DM
            if (present(m1_out)) rho_m1(:, :, isite) = rho_m1(:, :, isite) + DM_m1
         end do ! sites
      end do ! k-points

      if (this%k_mesh_distributed_active) then
         if (present(m1_out)) then
            call chain_complex(rho, size(rho), .false., rho_m1, size(rho_m1))
         else
            call chain_complex(rho, size(rho), .false.)
         end if
      end if

      ! Group average over the reduced mesh, P_ab convention (the legacy projector is its transpose)
      if (.not. use_complex_projector) rho = conjg(rho)
      call this%symmetrize_site_dm(rho)
      if (.not. use_complex_projector) rho = conjg(rho)
      if (present(m1_out)) call this%symmetrize_site_dm(rho_m1)

      do isite = 1, nsites
         do iorb = 0, min(3, lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%lmax)
            norb_l = 2*iorb + 1
            basis_offset = iorb**2  ! 0, 1, 4 mapping
            do ispin = 1, this%n_spin_components
               do alpha = 1, norb_l
                  do beta = 1, norb_l
                     idx_a = (ispin - 1) * n_orb_per_spin + basis_offset + alpha
                     idx_b = (ispin - 1) * n_orb_per_spin + basis_offset + beta
                     if (idx_a > n_orb_site .or. idx_b > n_orb_site) cycle
                     ! The fixed Liechtenstein path keeps the full coherent projector.
                     if (.not. m1_only_use) then
                        lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%ldm(iorb+1, ispin, alpha, beta) = &
                           merge(rho(idx_a, idx_b, isite), cmplx(real(rho(idx_a, idx_b, isite), rp), 0.0_rp, rp), use_complex_projector)
                     end if
                     if (present(m1_out)) m1_out(isite, iorb + 1, ispin, alpha, beta) = rho_m1(idx_a, idx_b, isite)
                  end do
               end do
            end do
         end do
      end do

      if (use_complex_projector .and. .not. m1_only_use) then
         do isite = 1, nsites
            do iorb = 0, min(3, lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%lmax)
               norb_l = 2*iorb + 1
               basis_offset = iorb**2
               if (basis_offset + norb_l > n_orb_per_spin) cycle
               do ispin = 1, min(2, this%n_spin_components)
                  call hubbard_density_diagnostics( &
                     lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%ldm(iorb+1, ispin, 1:norb_l, 1:norb_l), &
                     ldm_trace, ldm_min_eig, ldm_max_eig, ldm_antiherm, ldm_diag_info)
                  if (ldm_diag_info /= 0) then
                     call g_logger%warning('HUBBARD_LDM_KS_DIAGNOSTIC_FAILED site='//trim(int2str(isite))// &
                        ' l='//trim(int2str(iorb))//' spin='//trim(int2str(ispin))// &
                        ' info='//trim(int2str(ldm_diag_info)), __FILE__, __LINE__)
                  else
                     call root_info('HUBBARD_LDM_KS_PHYSICALITY site='//trim(int2str(isite))// &
                        ' l='//trim(int2str(iorb))//' spin='//trim(int2str(ispin))// &
                        ' trace='//trim(real2str(ldm_trace,'(F10.6)'))// &
                        ' eig_min='//trim(real2str(ldm_min_eig,'(F10.6)'))// &
                        ' eig_max='//trim(real2str(ldm_max_eig,'(F10.6)'))// &
                        ' max_antiherm='//trim(real2str(ldm_antiherm,'(E10.2)')), __FILE__, __LINE__)
                     if (ldm_antiherm > ldm_hermiticity_tol .or. ldm_min_eig < -ldm_occupation_tol .or. &
                         ldm_max_eig > 1.0_rp + ldm_occupation_tol .or. ldm_trace < -ldm_occupation_tol .or. &
                         ldm_trace > real(norb_l, rp) + ldm_occupation_tol) then
                        call g_logger%warning('HUBBARD_LDM_KS_OUT_OF_RANGE site='//trim(int2str(isite))// &
                           ' l='//trim(int2str(iorb))//' spin='//trim(int2str(ispin))// &
                           ' tolerances(herm,occ)=(1e-6,1e-3); values retained', __FILE__, __LINE__)
                     end if
                  end if
               end do
            end do
         end do
      end if

      deallocate(DM, rho)
      if (allocated(DM_m1)) deallocate(DM_m1, rho_m1)
      call root_info('calculate_ldm_from_eigenvectors: Eigenvector-based site density matrix -> ldm complete', __FILE__, __LINE__)
      
   end subroutine calculate_ldm_from_eigenvectors

   !> Full spinor d-shell P = sum f psi psi^dagger and M1 = sum f eps psi psi^dagger (index
   !> (sigma-1)*5+(m+3)) with the k weights and Fermi-Dirac occupations of the collinear producer.
   module subroutine calculate_spinor_dm(this, lattice_obj)
      use lattice_mod
      class(reciprocal), intent(inout) :: this
      type(lattice), intent(inout) :: lattice_obj
      integer :: nsites, n_orb_site, n_orb_per_spin
      real(rp) :: kT
      complex(rp), allocatable :: p(:, :, :), m1(:, :, :)
      real(rp) :: kweight, occ, axis(3)
      logical :: blochl
      real(rp), allocatable :: occ_t(:, :)
      complex(rp) :: element
      character(len=400) :: line
      integer :: ik, ik_global, isite, ib, s, t, i, j, row, col, plusbulk, site_offset

      kT = max(this%temperature*kB_Ry_per_K, 1.0e-10_rp)
      nsites = lattice_obj%nrec
      n_orb_site = size(this%eigenvectors, 1)/nsites
      n_orb_per_spin = n_orb_site/this%n_spin_components
      if (n_orb_per_spin < 9) return
      blochl = trim(this%dos_method) == 'blochl'
      if (blochl) call this%blochl_occupations(occ_t)
      allocate(p(nsites, 10, 10), m1(nsites, 10, 10))
      p = (0.0_rp, 0.0_rp)
      m1 = (0.0_rp, 0.0_rp)
      if (this%k_mesh_distributed_active) then
         call chain_complex(p, size(p), .true., m1, size(m1))
      end if
      do ik = 1, size(this%eigenvalues, 2)
         ik_global = local_k_index_to_global(this, ik)
         kweight = this%k_weights(ik_global)
         do isite = 1, nsites
            site_offset = (isite - 1)*n_orb_site
            do ib = 1, size(this%eigenvalues, 1)
               occ = band_occupation(occ_t, ib, ik_global, this%eigenvalues(ib, ik), this%fermi_level, kT)*kweight
               do s = 1, 2
                  do t = 1, 2
                     do i = 1, 5
                        do j = 1, 5
                           row = site_offset + (s - 1)*n_orb_per_spin + 4 + i
                           col = site_offset + (t - 1)*n_orb_per_spin + 4 + j
                           element = hubbard_projector_element(this%eigenvectors(row, ib, ik), &
                                                               this%eigenvectors(col, ib, ik))*occ
                           p(isite, (s - 1)*5 + i, (t - 1)*5 + j) = p(isite, (s - 1)*5 + i, (t - 1)*5 + j) + element
                           m1(isite, (s - 1)*5 + i, (t - 1)*5 + j) = &
                              m1(isite, (s - 1)*5 + i, (t - 1)*5 + j) + element*this%eigenvalues(ib, ik)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      if (this%k_mesh_distributed_active) then
         call chain_complex(p, size(p), .false., m1, size(m1))
      end if
      call symmetrize_d_spinor(this, p)
      call symmetrize_d_spinor(this, m1)
      do isite = 1, nsites
         plusbulk = lattice_obj%nbulk + isite
         if (.not. allocated(lattice_obj%symbolic_atoms(plusbulk)%potential%ldm_spinor)) &
            allocate(lattice_obj%symbolic_atoms(plusbulk)%potential%ldm_spinor(10, 10))
         if (.not. allocated(lattice_obj%symbolic_atoms(plusbulk)%potential%m1_spinor)) &
            allocate(lattice_obj%symbolic_atoms(plusbulk)%potential%m1_spinor(10, 10))
         lattice_obj%symbolic_atoms(plusbulk)%potential%ldm_spinor = p(isite, :, :)
         lattice_obj%symbolic_atoms(plusbulk)%potential%m1_spinor = m1(isite, :, :)
         call this%rf_density%get_axis(isite, axis)
         call lattice_obj%symbolic_atoms(plusbulk)%potential%set_hubbard_axis(axis, isite)
         call spinor_dm_summary(p(isite, :, :), lattice_obj%symbolic_atoms(plusbulk)%potential%ql(1, 2, 1), &
                                lattice_obj%symbolic_atoms(plusbulk)%potential%ql(1, 2, 2), axis, &
                                this%hamiltonian%hubbard_log_level /= 'summary', line)
         call root_info('HUBBARD_SPINOR_DM backend=ks atom='//trim(int2str(isite))//' '//trim(line), __FILE__, __LINE__)
      end do
      deallocate(p, m1)
   end subroutine calculate_spinor_dm

   !> Group average of per-site d-shell spinor matrices a(site, (s-1)*5+m, (t-1)*5+m'), embedded in the nb x nb site matrices.
   subroutine symmetrize_d_spinor(this, a)
      class(reciprocal), intent(in) :: this
      complex(rp), intent(inout) :: a(:, :, :)
      complex(rp), allocatable :: rho(:, :, :)
      integer :: s, t, i, j

      allocate(rho(nb, nb, size(a, 1)))
      rho = (0.0_rp, 0.0_rp)
      do s = 1, 2
         do t = 1, 2
            do i = 1, 5
               do j = 1, 5
                  rho((s - 1)*norb + 4 + i, (t - 1)*norb + 4 + j, :) = a(:, (s - 1)*5 + i, (t - 1)*5 + j)
               end do
            end do
         end do
      end do
      call this%symmetrize_site_dm(rho)
      do s = 1, 2
         do t = 1, 2
            do i = 1, 5
               do j = 1, 5
                  a(:, (s - 1)*5 + i, (t - 1)*5 + j) = rho((s - 1)*norb + 4 + i, (t - 1)*norb + 4 + j, :)
               end do
            end do
         end do
      end do
   end subroutine symmetrize_d_spinor

end submodule reciprocal_projection
