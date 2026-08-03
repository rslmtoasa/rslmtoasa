submodule (reciprocal_mod) reciprocal_projection
#ifdef USE_MPI
   use mpi
#endif
   implicit none

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
   integer :: ik, ik_global, ib, ie, iorb, i, isite
   integer :: n_orb_per_spin, orb_start, site_orb_start, n_orb_site
   integer :: lstart(4), lend(4)
   real(rp) :: weight, orbital_char, energy
   real(rp) :: gaussian_weight, sigma_squared, sigma_use
   complex(rp) :: psi_element
   real(rp) :: dos_integral, norm_factor
   integer :: nbands
   ! Additional locals for projected DOS integration/normalization
   real(rp) :: proj_integral, e_low, e_high
   integer :: iei
   ! Per-site orbital offsets for mixed atom types
   integer, dimension(:), allocatable :: site_orb_offset
   real(rp) :: mx_char, my_char, mz_char, local_char
   real(rp) :: axis(3)

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
   n_orb_per_spin = n_orb_site/this%n_spin_components

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

   allocate(site_orb_offset(this%n_sites + 1))
   do isite = 1, this%n_sites + 1
      site_orb_offset(isite) = (isite - 1) * n_orb_site
   end do
   lstart = [1, 2, 5, 10]
   lend = [1, 4, 9, 16]

   ! Diagnostic logging
   call root_info('project_dos_orbitals_gaussian: n_sites = ' // trim(int2str(this%n_sites)) // &
                     ', nbands = ' // trim(int2str(nbands)) // &
                     ', max_orb_channels = ' // trim(int2str(this%max_orbs)) // &
                     ', eigenvector size = ' // trim(int2str(size(this%eigenvectors, 1))), __FILE__, __LINE__)
   ! call g_logger%info('project_dos_orbitals_gaussian: site_orb_offset = [' // &
   !                   trim(int2str(site_orb_offset(1))) // ', ' // &
   !                   trim(int2str(site_orb_offset(2))) // ', ' // &
   !                   trim(int2str(site_orb_offset(3))) // ', ' // &
   !                   trim(int2str(site_orb_offset(4))) // ', ' // &
   !                   trim(int2str(site_orb_offset(5))) // ']', __FILE__, __LINE__)

   ! Calculate projected DOS (raw) - same weights as total DOS
   do ie = 1, this%n_energy_points
      energy = this%dos_energy_grid(ie)  ! Already in Ry

      do ik = 1, size(this%eigenvalues, 2)
         ik_global = local_k_index_to_global(this, ik)
         do ib = 1, nbands  ! Loop over all bands, not just max_orb_channels
            ! Skip if eigenvalue is far from current energy
            if (abs(this%eigenvalues(ib, ik) - energy) > 5.0_rp * sigma_use) cycle

            ! Calculate Gaussian weight (already normalized)
            gaussian_weight = exp(-((energy - this%eigenvalues(ib, ik))**2) / (2.0_rp * sigma_squared))
            gaussian_weight = gaussian_weight / (sigma_use * sqrt(2.0_rp * 3.141592653589793_rp))
            
            if (abs(gaussian_weight) < 1.0e-10_rp) cycle

            ! Apply k-point weight
            weight = gaussian_weight * this%k_weights(ik_global)

            do isite = 1, this%n_sites
               ! Site-blocked layout: eigenvectors are [site1(18), site2(18), ...]
               ! where each site block contains [orb1_up...orb9_up, orb1_dn...orb9_dn]
               site_orb_start = site_orb_offset(isite)

               call get_site_spin_axis(this, isite, axis)
               do iorb = 1, 4
                  orbital_char = 0.0_rp
                  mx_char = 0.0_rp
                  my_char = 0.0_rp
                  mz_char = 0.0_rp
                  if (lstart(iorb) <= n_orb_per_spin) then
                     call compute_spinor_block_projection(this, ik, ib, site_orb_start, n_orb_per_spin, &
                                                          lstart(iorb), lend(iorb), orbital_char, mx_char, my_char, mz_char)
                  end if
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

#ifdef USE_MPI
   if (this%k_mesh_distributed_active) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%projected_dos, product(shape(this%projected_dos)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%projected_dos_moments, product(shape(this%projected_dos_moments)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%dos_mx_tot, this%n_energy_points, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%dos_my_tot, this%n_energy_points, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%dos_mz_tot, this%n_energy_points, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
   end if
#endif

   ! Normalize projected DOS so integrated sum over projections equals nbands
   call root_info('project_dos_orbitals_gaussian: Projection completed (raw)', __FILE__, __LINE__)
   if (this%n_energy_points > 2) then
      proj_integral = 0.0_rp
      do iei = 1, this%n_energy_points - 1
         e_low = this%dos_energy_grid(iei)
         e_high = this%dos_energy_grid(iei+1)
         proj_integral = proj_integral + 0.5_rp * ( sum(this%projected_dos(:, :, :, iei)) + sum(this%projected_dos(:, :, :, iei+1)) ) * (e_high - e_low)
      end do

      if (abs(proj_integral) > 1.0e-12_rp) then
         norm_factor = 1.0_rp ! real(nbands, rp) / proj_integral
         this%projected_dos = this%projected_dos * norm_factor
         call root_info('project_dos_orbitals_gaussian: Normalized projected DOS by factor ' // trim(real2str(norm_factor, '(F10.6)')), __FILE__, __LINE__)
      else
         call g_logger%warning('project_dos_orbitals_gaussian: projected DOS integral is zero, skipping normalization', __FILE__, __LINE__)
      end if

      ! Diagnostic: check mid-energy ratio after normalization
      ie = this%n_energy_points / 2
      if (abs(this%total_dos(ie)) > 1.0e-12_rp) then
         call root_info('project_dos_orbitals_gaussian: At mid-energy, proj/total ratio (post-norm) = ' // &
            trim(real2str(sum(this%projected_dos(:, :, :, ie)) / this%total_dos(ie), '(F10.6)')), __FILE__, __LINE__)
      end if
   end if

   deallocate(site_orb_offset)
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
      if (atom_idx >= 1 .and. atom_idx <= size(this%lattice%symbolic_atoms)) then
         axis = this%lattice%symbolic_atoms(atom_idx)%potential%mom(:)
      end if

      axis_norm = sqrt(sum(axis**2))
      if (axis_norm > tiny(1.0_rp)) then
         axis = axis/axis_norm
      else
         axis = [0.0_rp, 0.0_rp, 1.0_rp]
      end if
   end subroutine get_site_spin_axis

#ifdef USE_MPI
   !> @brief Synchronize lattice local-density-matrix data across MPI ranks.
   !> @param[inout] lattice_obj Lattice object whose ldm-like storage is synchronized.
   module subroutine sync_lattice_ldm(lattice_obj)
      use lattice_mod
      type(lattice), intent(inout) :: lattice_obj
      real(rp), allocatable :: ldm_comm(:, :, :, :)
      integer :: max_flat_ldm, local_flat
      integer :: na_glob, plusbulk, lcount_ldm, l, ispin

      max_flat_ldm = (2*lmax_basis + 1)*(2*lmax_basis + 1)
      allocate(ldm_comm(lattice_obj%nrec, lmax_basis + 1, 2, max_flat_ldm))
      ldm_comm(:, :, :, :) = 0.0_rp
      do na_glob = 1, lattice_obj%nrec
         plusbulk = lattice_obj%nbulk + na_glob
         call lattice_obj%symbolic_atoms(plusbulk)%potential%flatten_ldm()
         lcount_ldm = min(lattice_obj%symbolic_atoms(plusbulk)%potential%lmax, lmax_basis) + 1
         do l = 1, lcount_ldm
            local_flat = (2*l - 1)*(2*l - 1)
            do ispin = 1, 2
               ldm_comm(na_glob, l, ispin, 1:local_flat) = &
                  lattice_obj%symbolic_atoms(plusbulk)%potential%ldm_flatten(l, ispin, 1:local_flat)
            end do
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE, ldm_comm, product(shape(ldm_comm)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      do na_glob = 1, lattice_obj%nrec
         plusbulk = lattice_obj%nbulk + na_glob
         lattice_obj%symbolic_atoms(plusbulk)%potential%ldm_flatten(:, :, :) = 0.0_rp
         lcount_ldm = min(lattice_obj%symbolic_atoms(plusbulk)%potential%lmax, lmax_basis) + 1
         do l = 1, lcount_ldm
            local_flat = (2*l - 1)*(2*l - 1)
            do ispin = 1, 2
               lattice_obj%symbolic_atoms(plusbulk)%potential%ldm_flatten(l, ispin, 1:local_flat) = &
                  ldm_comm(na_glob, l, ispin, 1:local_flat)
            end do
         end do
         call lattice_obj%symbolic_atoms(plusbulk)%potential%expand_ldm()
      end do
      deallocate(ldm_comm)
   end subroutine sync_lattice_ldm
#endif

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

   !> @brief Calculate projected DOS with tetrahedron integration.
   !> @details Averages orbital/spinor projection weights over tetrahedron
   !>          corners and accumulates scalar plus directional projected DOS.
   !> @param[inout] this Reciprocal object receiving projected DOS arrays.
   module subroutine project_dos_orbitals_tetrahedron(this)
      class(reciprocal), intent(inout) :: this

      ! Local variables
      integer :: i_energy, i_tet, i_corner, i_band, iorb, isite
      integer :: i_start, i_end
      integer :: n_orb_per_spin, orb_start, site_orb_start, ik, n_orb_site
      integer :: ie, nbands
      integer :: tet_start, tet_end, tet_count
      real(rp) :: energy, dos_contrib, orbital_char_avg, orbital_char, de
      real(rp) :: total_dos_integral, proj_dos_integral, norm_factor
      real(rp), dimension(4) :: e_corners, sorted_e, orbital_chars
      real(rp), dimension(4) :: mx_chars, my_chars, mz_chars
      integer :: lstart(4), lend(4)
      real(rp) :: e1, e2, e3, e4, x, C
      real(rp) :: tet_weight
      real(rp), parameter :: TOL = 1.0e-10_rp
      integer :: n_tet_ir
      ! Per-site orbital offsets for mixed atom types
      integer, dimension(:), allocatable :: site_orb_offset
      real(rp) :: mx_char_avg, my_char_avg, mz_char_avg, local_char
      real(rp) :: axis(3)
      real(rp), allocatable :: site_axes(:, :)
      real(rp), allocatable :: dos_line(:)
      real(rp), allocatable :: local_projected_dos(:, :, :, :)
      real(rp), allocatable :: local_projected_dos_moments(:, :, :, :, :)
      real(rp), allocatable :: local_dos_mx(:), local_dos_my(:), local_dos_mz(:)

      call g_logger%info('project_dos_orbitals_tetrahedron: Starting tetrahedron orbital projection calculation', __FILE__, __LINE__)

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
      lstart = [1, 2, 5, 10]
      lend = [1, 4, 9, 16]

      call this%ensure_full_mesh_for_spinor_integrations('project_dos_orbitals_tetrahedron')

      call g_logger%info('project_dos_orbitals_tetrahedron: Projecting onto ' // trim(int2str(this%n_sites)) // &
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

      allocate(site_orb_offset(this%n_sites + 1))
      allocate(site_axes(3, this%n_sites))
      do isite = 1, this%n_sites + 1
         site_orb_offset(isite) = (isite - 1) * n_orb_site
      end do
      do isite = 1, this%n_sites
         call get_site_spin_axis(this, isite, site_axes(:, isite))
      end do

      ! Setup tetrahedra if not already done
      if (.not. allocated(this%tetrahedra)) then
         call this%setup_tetrahedra()
      end if
      de = (this%dos_energy_grid(this%n_energy_points) - this%dos_energy_grid(1)) / &
           real(this%n_energy_points - 1, rp)
      if (de <= 0.0_rp) then
         call g_logger%fatal('project_dos_orbitals_tetrahedron: Invalid DOS energy grid spacing', __FILE__, __LINE__)
      end if
      n_tet_ir = this%n_tetrahedra
      call get_mpi_range(rank, n_tet_ir, tet_start, tet_end, tet_count, region_tag='projected tetra')

      ! Tetra/band-first traversal with one tetra shape evaluation per band.
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(this, site_orb_offset, site_axes, lstart, lend, n_orb_per_spin, n_tet_ir, &
      !$OMP        tet_start, tet_end, de) &
      !$OMP PRIVATE(i_tet, i_band, i_corner, ik, e_corners, sorted_e, e1, e2, e3, e4, i_start, i_end, &
      !$OMP         i_energy, energy, dos_contrib, x, C, tet_weight, isite, site_orb_start, axis, iorb, &
      !$OMP         orbital_char_avg, mx_char_avg, my_char_avg, mz_char_avg, orbital_char, orbital_chars, &
      !$OMP         mx_chars, my_chars, mz_chars, local_char, dos_line, local_projected_dos, &
      !$OMP         local_projected_dos_moments, local_dos_mx, local_dos_my, local_dos_mz)
      allocate(dos_line(this%n_energy_points))
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
      do i_tet = tet_start, tet_end
         do i_band = 1, size(this%eigenvalues, 1)
            do i_corner = 1, 4
               ik = this%tetrahedra(i_corner, i_tet)
               e_corners(i_corner) = this%eigenvalues(i_band, ik)
            end do
            call sort4(e_corners, sorted_e)
            e1 = sorted_e(1)
            e2 = sorted_e(2)
            e3 = sorted_e(3)
            e4 = sorted_e(4)
            if (e4 <= this%dos_energy_grid(1) .or. e1 >= this%dos_energy_grid(this%n_energy_points)) cycle
            if (abs(e2 - e1) < TOL .or. abs(e3 - e1) < TOL .or. abs(e4 - e1) < TOL .or. &
                abs(e3 - e2) < TOL .or. abs(e4 - e2) < TOL .or. abs(e4 - e3) < TOL) cycle

            i_start = max(1, int(floor((e1 - this%dos_energy_grid(1)) / de)) + 1)
            i_end = min(this%n_energy_points, int(ceiling((e4 - this%dos_energy_grid(1)) / de)) + 1)
            tet_weight = this%tetrahedron_volumes(i_tet)

            dos_line = 0.0_rp
            do i_energy = i_start, i_end
               energy = this%dos_energy_grid(i_energy)
               if (trim(this%dos_method) == 'blochl') then
                  if (energy <= e1 .or. energy >= e4) cycle
                  if (energy <= e2) then
                     dos_contrib = 3.0_rp * (energy - e1)**2 / ((e4 - e1) * (e3 - e1) * (e2 - e1))
                  else if (energy <= e3) then
                     C = 1.0_rp / ((e4 - e1) * (e3 - e1))
                     dos_contrib = C * (3.0_rp * (e2 - e1) + 6.0_rp * (energy - e2) - &
                                    3.0_rp * (e3 + e4 - e1 - e2) * (energy - e2)**2 / &
                                    ((e3 - e2) * (e4 - e2)))
                  else
                     dos_contrib = 3.0_rp * (e4 - energy)**2 / ((e4 - e1) * (e4 - e2) * (e4 - e3))
                  end if
               else
                  if (energy < e1 .or. energy >= e4) cycle
                  if (energy < e2) then
                     x = energy - e1
                     dos_contrib = 3.0_rp * x * x / ((e2 - e1) * (e3 - e1) * (e4 - e1))
                  else if (energy < e3) then
                     x = energy - e2
                     dos_contrib = 3.0_rp * (e2 - e1) / ((e3 - e1) * (e4 - e1)) + &
                                 x * (6.0_rp / ((e3 - e1) * (e4 - e1)) + &
                                 x * (3.0_rp * (e1 + e2 - e3 - e4) / &
                                 ((e3 - e1) * (e4 - e1) * (e3 - e2) * (e4 - e2))))
                  else
                     x = energy - e4
                     dos_contrib = 3.0_rp * x * x / ((e4 - e3) * (e4 - e2) * (e4 - e1))
                  end if
               end if
               dos_line(i_energy) = dos_contrib * tet_weight
            end do

            do isite = 1, this%n_sites
               site_orb_start = site_orb_offset(isite)
               axis = site_axes(:, isite)
               do iorb = 1, 4
                  orbital_char_avg = 0.0_rp
                  mx_char_avg = 0.0_rp
                  my_char_avg = 0.0_rp
                  mz_char_avg = 0.0_rp
                  if (lstart(iorb) <= n_orb_per_spin) then
                     do i_corner = 1, 4
                        ik = this%tetrahedra(i_corner, i_tet)
                        call compute_spinor_block_projection(this, ik, i_band, site_orb_start, n_orb_per_spin, &
                                                             lstart(iorb), lend(iorb), orbital_char, mx_chars(i_corner), &
                                                             my_chars(i_corner), mz_chars(i_corner))
                        orbital_chars(i_corner) = orbital_char
                     end do
                     orbital_char_avg = sum(orbital_chars) / 4.0_rp
                     mx_char_avg = sum(mx_chars) / 4.0_rp
                     my_char_avg = sum(my_chars) / 4.0_rp
                     mz_char_avg = sum(mz_chars) / 4.0_rp
                  end if
                  if (abs(orbital_char_avg) < 1.0e-14_rp) cycle

                  local_char = axis(1)*mx_char_avg + axis(2)*my_char_avg + axis(3)*mz_char_avg
                  do i_energy = i_start, i_end
                     dos_contrib = dos_line(i_energy)
                     if (dos_contrib == 0.0_rp) cycle
                     local_projected_dos(isite, iorb, 1, i_energy) = local_projected_dos(isite, iorb, 1, i_energy) + &
                                                                      0.5_rp*(orbital_char_avg + local_char)*dos_contrib
                     local_projected_dos(isite, iorb, 2, i_energy) = local_projected_dos(isite, iorb, 2, i_energy) + &
                                                                      0.5_rp*(orbital_char_avg - local_char)*dos_contrib
                     local_projected_dos_moments(isite, iorb, 1, 1, i_energy) = &
                        local_projected_dos_moments(isite, iorb, 1, 1, i_energy) + mx_char_avg*dos_contrib
                     local_projected_dos_moments(isite, iorb, 1, 2, i_energy) = &
                        local_projected_dos_moments(isite, iorb, 1, 2, i_energy) + my_char_avg*dos_contrib
                     local_projected_dos_moments(isite, iorb, 1, 3, i_energy) = &
                        local_projected_dos_moments(isite, iorb, 1, 3, i_energy) + mz_char_avg*dos_contrib
                     local_dos_mx(i_energy) = local_dos_mx(i_energy) + mx_char_avg*dos_contrib
                     local_dos_my(i_energy) = local_dos_my(i_energy) + my_char_avg*dos_contrib
                     local_dos_mz(i_energy) = local_dos_mz(i_energy) + mz_char_avg*dos_contrib
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

      deallocate(dos_line)
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

      deallocate(site_orb_offset)
      deallocate(site_axes)

      ! Normalize projected DOS to match total DOS normalization
      ! When using Blöchl or standard tetrahedron method, the raw DOS needs normalization
      ! to integrate to the correct number of states. We must apply the same normalization
      ! to projected DOS for consistency.
      if (allocated(this%total_dos) .and. this%n_energy_points > 2) then
         nbands = size(this%eigenvalues, 1)
         total_dos_integral = 0.0_rp
         do ie = 1, this%n_energy_points - 1
            total_dos_integral = total_dos_integral + &
               0.5_rp * (this%total_dos(ie) + this%total_dos(ie+1)) * &
               (this%dos_energy_grid(ie+1) - this%dos_energy_grid(ie))
         end do
         
         ! The total_dos is already normalized, so total_dos_integral ≈ nbands
         ! Calculate integral of projected DOS (raw, unnormalized)
         proj_dos_integral = 0.0_rp
         do ie = 1, this%n_energy_points - 1
            ! Sum over all sites, orbitals, and spins
            proj_dos_integral = proj_dos_integral + &
               0.5_rp * (sum(this%projected_dos(:,:,:,ie)) + sum(this%projected_dos(:,:,:,ie+1))) * &
               (this%dos_energy_grid(ie+1) - this%dos_energy_grid(ie))
         end do
         
         ! Normalize projected DOS by the same factor
         if (abs(proj_dos_integral) > 1.0e-10_rp) then
            norm_factor = total_dos_integral / proj_dos_integral
            this%projected_dos = this%projected_dos * norm_factor
            call g_logger%info('project_dos_orbitals_tetrahedron: Normalized projected DOS by factor ' // &
                              trim(real2str(norm_factor, '(F10.6)')), __FILE__, __LINE__)
         end if
      end if

      call g_logger%info('project_dos_orbitals_tetrahedron: Tetrahedron orbital projection calculation completed', __FILE__, __LINE__)
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

   integer :: isite, iorb, ispin, ie, n_energy
   real(rp) :: energy, dos_value, fermi_weight
   real(rp) :: m0, m1, m2
   real(rp), dimension(:), allocatable :: integrand, fermi_dist
   real(rp) :: kT, fermi_arg
   real(rp) :: total_occupation, expected_electrons, total_occupation_alt
      real(rp), allocatable :: energy_grid_ry(:)
      real(rp), parameter :: eV_to_Ry = 0.073498618_rp

   call root_info('calculate_band_moments: Starting calculation', __FILE__, __LINE__)
   call root_info('calculate_band_moments: DEBUG - this%auto_find_fermi = ' // &
                  merge('TRUE ', 'FALSE', this%auto_find_fermi) // &
                  ', this%total_electrons = ' // trim(real2str(this%total_electrons, '(F10.5)')), &
                  __FILE__, __LINE__)
   call root_info('calculate_band_moments: DEBUG - Fermi level on entry = ' // &
                  trim(real2str(this%fermi_level, '(F10.6)')) // ' Ry', __FILE__, __LINE__)

   if (this%total_electrons <= 1.0e-3_rp) then
      if (this%lattice%nbulk_bulk > 0) then
         this%total_electrons = real(sum(this%lattice%symbolic_atoms(1:this%lattice%nbulk_bulk)%element%valence), rp)
         call g_logger%info('calculate_band_moments: Recovered total_electrons from valence = ' // &
                           trim(real2str(this%total_electrons, '(F10.5)')), __FILE__, __LINE__)
      else if (this%lattice%nrec > 0) then
         this%total_electrons = real(sum(this%lattice%symbolic_atoms(1:this%lattice%nrec)%element%valence), rp)
         call g_logger%warning('calculate_band_moments: nbulk_bulk<=0; recovered total_electrons using nrec span.', __FILE__, __LINE__)
      end if
   end if

   ! EF is already canonical here: calculate_density_of_states solves it from
   ! eigenvalues/k weights before any DOS projection.  This routine must never
   ! replace it with a grid/window-dependent DOS root.
   call g_logger%info('calculate_band_moments: Using canonical/fixed Fermi level = ' // &
                      trim(real2str(this%fermi_level, '(F 8.5)')) // ' Ry', __FILE__, __LINE__)

   if (this%use_symmetry_reduction .and. &
       (trim(this%dos_method) == 'tetrahedron' .or. trim(this%dos_method) == 'blochl') .and. &
       this%nk_total < this%nk_mesh(1) * this%nk_mesh(2) * this%nk_mesh(3)) then
      call this%ensure_full_mesh_for_spinor_integrations('calculate_band_moments')
      call this%project_dos_orbitals_tetrahedron()
   end if

   if (allocated(this%band_moments)) deallocate(this%band_moments)
   allocate(this%band_moments(this%n_sites, this%n_orb_types, this%n_spin_components, 3))
   this%band_moments = 0.0_rp

   n_energy = this%n_energy_points
   allocate(integrand(n_energy))
   allocate(fermi_dist(n_energy))
   allocate(energy_grid_ry(n_energy))
   energy_grid_ry = this%dos_energy_grid  ! Already in Ry

   ! Boltzmann constant: kB = 6.3336814e-6 Ry/K
   kT = max(this%temperature * 6.3336814e-6_rp, 1.0e-10_rp)  ! Ry/K; canonical FD floor
   
   call g_logger%info('calculate_band_moments: kT = ' // trim(real2str(kT, '(ES12.5)')) // &
                     ' Ry at T = ' // trim(real2str(this%temperature, '(F8.2)')) // ' K', &
                     __FILE__, __LINE__)

   ! Pre-calculate Fermi-Dirac distribution
   do ie = 1, n_energy
      energy = energy_grid_ry(ie)
      fermi_arg = (energy - this%fermi_level) / kT
      
      if (fermi_arg > 50.0_rp) then
         fermi_dist(ie) = 0.0_rp
      else if (fermi_arg < -50.0_rp) then
         fermi_dist(ie) = 1.0_rp
      else
         fermi_dist(ie) = 1.0_rp / (exp(fermi_arg) + 1.0_rp)
      end if
   end do

   ! Calculate total occupation to verify normalization
   integrand = this%total_dos * fermi_dist
   total_occupation = trapezoidal_integral(energy_grid_ry, integrand)
   
   ! Also calculate using the same method as find_fermi_level_from_dos for comparison
   total_occupation_alt = this%integrate_dos_up_to_energy(this%fermi_level, kT)
   
   call g_logger%info('calculate_band_moments: Total occupation (method 1) = ' // &
                     trim(real2str(total_occupation, '(F10.5)')) // &
                     ', (method 2) = ' // trim(real2str(total_occupation_alt, '(F10.5)')) // &
                     ', expected = ' // trim(real2str(this%total_electrons, '(F10.5)')), &
                     __FILE__, __LINE__)

   ! Calculate moments for each projection
   do isite = 1, this%n_sites
      do iorb = 1, this%n_orb_types
         do ispin = 1, this%n_spin_components

            ! m0: occupation = ∫ DOS(E) * f(E) dE
            integrand = this%projected_dos(isite, iorb, ispin, :) * fermi_dist
            m0 = trapezoidal_integral(energy_grid_ry, integrand)

            ! m1: band center = ∫ E * DOS(E) * f(E) dE / m0
            if (abs(m0) > 1.0e-12_rp) then
               integrand = energy_grid_ry * this%projected_dos(isite, iorb, ispin, :) * fermi_dist
               m1 = trapezoidal_integral(energy_grid_ry, integrand) / m0
            else
               m1 = 0.0_rp
            end if

            ! m2: band width = sqrt(∫ (E - m1)² * DOS(E) * f(E) dE / m0)
            if (abs(m0) > 1.0e-12_rp) then
               ! Use energy grid in Ry for the variance integral (units must match projected_dos which
               ! was integrated/normalized in Ry). Previously we accidentally passed the eV grid here.
               integrand = (energy_grid_ry - m1)**2 * &
                          this%projected_dos(isite, iorb, ispin, :) * fermi_dist
               m2 = trapezoidal_integral(energy_grid_ry, integrand) / m0
               m2 = sqrt(max(m2, 0.0_rp))
            else
               m2 = 0.0_rp
            end if

            this%band_moments(isite, iorb, ispin, 1) = m0
            this%band_moments(isite, iorb, ispin, 2) = m1
            this%band_moments(isite, iorb, ispin, 3) = m2
         end do
      end do
   end do

   ! Calculate sum of all m0 moments (should equal total occupation)
   m0 = 0.0_rp
   do isite = 1, this%n_sites
      do iorb = 1, this%n_orb_types
         do ispin = 1, this%n_spin_components
            m0 = m0 + this%band_moments(isite, iorb, ispin, 1)
         end do
      end do
   end do
   
   call g_logger%info('calculate_band_moments: Sum of all m0 moments = ' // &
                     trim(real2str(m0, '(F10.5)')) // &
                     ', total occupation = ' // trim(real2str(total_occupation, '(F10.5)')), &
                     __FILE__, __LINE__)

   deallocate(integrand, fermi_dist)

   call g_logger%info('calculate_band_moments: Completed', __FILE__, __LINE__)
end subroutine calculate_band_moments

   !> @brief Reconstruct lattice local density matrices from projected DOS moments.
   !> @param[inout] this Reciprocal object containing projected DOS/moments.
   !> @param[inout] lattice_obj Lattice object receiving local density matrices.
   module subroutine calculate_ldm_from_projected_dos(this, lattice_obj)
      use lattice_mod
      class(reciprocal), intent(inout) :: this
      type(lattice), intent(inout) :: lattice_obj
      
      ! Local variables
      integer :: isite
      logical :: operator_changed
      
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
      
      ! Initialize all LDM to zero
      do isite = 1, lattice_obj%nrec
         lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%ldm(:,:,:,:) = 0.0_rp
      end do
      
      ! Method 2: Full LDM (including off-diagonal) from eigenvector projections
      ! This properly captures orbital hybridization and is more accurate than diagonal-only
      call calculate_ldm_from_eigenvectors(this, lattice_obj)
      
      ! Save flattened copy for each site
      do isite = 1, lattice_obj%nrec
         call lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%flatten_ldm()
      end do
      
      ! Print LDM for debugging
      call root_info('calculate_ldm_from_projected_dos: LDM calculation complete', __FILE__, __LINE__)
      
   end subroutine calculate_ldm_from_projected_dos

   !> @brief Reconstruct lattice local density matrices directly from eigenvectors.
   !> @details Accumulates occupied spinor density matrices over k-points and
   !>          bands, then synchronizes the lattice result across MPI ranks.
   !> @param[inout] this Reciprocal object containing eigenvectors and occupations.
   !> @param[inout] lattice_obj Lattice object receiving local density matrices.
   module subroutine calculate_ldm_from_eigenvectors(this, lattice_obj)
      use lattice_mod
      class(reciprocal), intent(inout) :: this
      type(lattice), intent(inout) :: lattice_obj
      
      ! Local variables
      integer :: ik, ik_global, ib, isite, alpha, beta, iorb, norb_l
      integer :: ispin, basis_offset, site_offset, n_orb_site, n_orb_per_spin
      real(rp) :: kweight, fermi_occ
      complex(rp), allocatable :: DM(:,:)
      complex(rp) :: psi_a, psi_b
      integer :: nbands, nsites, i
      real(rp) :: trace_tot, trace_up, trace_dn, max_imag
      integer :: idx_a, idx_b

      call root_info('calculate_ldm_from_eigenvectors: Computing site density matrices from eigenvectors', __FILE__, __LINE__)

      call this%ensure_full_mesh_for_spinor_integrations('calculate_ldm_from_eigenvectors')

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

      ! Allocate temporary DM (n_orb_site x n_orb_site)
      allocate(DM(n_orb_site, n_orb_site))

      ! Zero per-site LDM storage before accumulation
      do isite = 1, nsites
         lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%ldm(:,:,:,:) = 0.0_rp
      end do

      ! Accumulate density matrix per site
      ! Loop order: k-points -> sites -> bands for better cache locality
      do ik = 1, size(this%eigenvalues, 2)
         ik_global = local_k_index_to_global(this, ik)
         kweight = this%k_weights(ik_global)
         
         do isite = 1, nsites
            ! Reset DM for this site at this k-point
            DM = (0.0_rp, 0.0_rp)
            site_offset = (isite - 1) * n_orb_site
            
            ! Accumulate contributions from all occupied bands at this k-point
            do ib = 1, nbands
               if (this%eigenvalues(ib, ik) > this%fermi_level) cycle
               fermi_occ = 1.0_rp
               
               ! Build density matrix: DM_{ab} += f * w_k * ψ*_a ψ_b
               do alpha = 1, n_orb_site
                  psi_a = this%eigenvectors(site_offset + alpha, ib, ik)
                  do beta = 1, n_orb_site
                     psi_b = this%eigenvectors(site_offset + beta, ib, ik)
                     DM(alpha, beta) = DM(alpha, beta) + conjg(psi_a) * psi_b * fermi_occ * kweight
                  end do
               end do
            end do ! bands
            
            ! Diagnostics: check DM properties for first k-point
            if (ik_global == 1) then
               trace_tot = 0.0_rp
               trace_up = 0.0_rp
               trace_dn = 0.0_rp
               max_imag = 0.0_rp
               
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
                                 ', max_imag=' // trim(real2str(max_imag, '(E10.2)')), __FILE__, __LINE__)
               
               if (max_imag > 1.0e-6_rp) then
                  call g_logger%warning('calculate_ldm_from_eigenvectors: Large imaginary part in DM diagonal: ' // &
                                       trim(real2str(max_imag, '(E10.2)')), __FILE__, __LINE__)
               end if
            end if
            
            ! Project DM into ldm structure: iterate over angular momentum channels
            do iorb = 0, min(3, lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%lmax)
               norb_l = 2*iorb + 1
               basis_offset = iorb**2  ! 0, 1, 4 mapping
               
               do ispin = 1, this%n_spin_components
                  do alpha = 1, norb_l
                     do beta = 1, norb_l
                        ! Map to DM indices: spin_block_offset + basis_offset + m
                        idx_a = (ispin - 1) * n_orb_per_spin + basis_offset + alpha
                        idx_b = (ispin - 1) * n_orb_per_spin + basis_offset + beta
                        if (idx_a > n_orb_site .or. idx_b > n_orb_site) cycle
                        
                        ! Accumulate real part (imaginary should be ~0 for proper DM)
                        lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%ldm(iorb+1, ispin, alpha, beta) = &
                           lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%ldm(iorb+1, ispin, alpha, beta) + &
                           real(DM(idx_a, idx_b), rp)
                     end do
                  end do
               end do
            end do
            
         end do ! sites
      end do ! k-points

#ifdef USE_MPI
      if (this%k_mesh_distributed_active) then
         call sync_lattice_ldm(lattice_obj)
      end if
#endif

      deallocate(DM)
      call root_info('calculate_ldm_from_eigenvectors: Eigenvector-based 18x18 DM -> ldm complete', __FILE__, __LINE__)
      
   end subroutine calculate_ldm_from_eigenvectors

end submodule reciprocal_projection
