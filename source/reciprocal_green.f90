!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! SUBMODULE: reciprocal_green
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> k-space Green's-function engine (milestone B2, flagship). One filler API,
!> two backends, populating the SAME arrays on the `green` object that the
!> real-space recursion route fills (`gij/gji`, the `gij_eta` Fermi-point
!> ladder, and the torque-resolved `ginmag`/`gi{x,y,z}` families). Downstream
!> consumers (`exchange` via `auxiliary_gij`, damping, future chi) then run
!> unchanged, by construction -- filling the canonical arrays IS the
!> route-agnosticism; there is no adapter layer.
!>
!> Backend E -- eigenpair / strict Lehmann (Sigma = 0):
!>   G_ij(z) = (1/N_k) sum_{k,n} e^{i k.dR_ij} psi_{i,n}(k) psi^dagger_{j,n}(k)
!>                                              / (z - eps_{nk})
!>   One Hermitian eigensolve per k, amortized over all energies and pairs.
!>
!> Backend D -- direct Dyson inversion (Sigma(z) != 0):
!>   G(k,z) = [ z*S - H(k) - Sigma(z) ]^{-1}    per (k, z)
!>   batched nb x nb factor/solve; Sigma enters via `sigma_provider`.
!>   Permanent CI invariant: backend D with Sigma = 0 == backend E.
!>
!> THIS FILE (task B2.1) provides the skeleton: the energy-contour adoption
!> (`build_green_contour`) and the `fill_green` dispatcher. The backend cores
!> land in B2.2 (E) and B2.4 (D); the eta ladder + torque components in B2.3.
!>
!> Contour convention (gate G-B2-1): the backends deliver G on the SAME grid
!> and contour the recursion route uses, so `exchange`'s Simpson integral
!> (`simpson_f(..., en%ene, en%fermi, en%nv1, ...)`) is valid without change:
!>   * grid    : `en%ene(1:size(en%ene))`, the real-axis mesh from `e_mesh`
!>               (energy_min .. above fermi; `en%channels_ldos + 10` points).
!>   * contour : retarded, z(ie) = en%ene(ie) + i*green_eta, matching `bgreen`
!>               (z = e(ei) + eta with eta = i*eta_imag).
!>   * fermi   : physical `en%fermi` (NEVER the chebfermi-scaled variable).
!>   * repr.   : screened/auxiliary LMTO blocks, pre-`auxiliary_gij`.
!------------------------------------------------------------------------------
submodule(reciprocal_mod) reciprocal_green
   use lehmann_kernel_mod, only: lehmann_pair_block
   use mpi_mod, only: start_atom, end_atom, g2l_map
   implicit none

contains

   !> @brief Build the retarded complex-energy contour for the k-space filler.
   !> @details Adopts the real-space route's energy grid verbatim (gate
   !>          G-B2-1): z(ie) = en%ene(ie) + i*this%green_eta over the full
   !>          `en%ene` mesh. The grid must already be prepared by `en%e_mesh`;
   !>          `fill_green` guarantees this before calling.
   !> @param[in]  this      Reciprocal object (supplies green_eta broadening).
   !> @param[in]  en        Energy object holding the prepared real-axis grid.
   !> @param[out] z_contour Allocated retarded contour, size = size(en%ene).
   module subroutine build_green_contour(this, en, z_contour)
      class(reciprocal), intent(in) :: this
      class(energy), intent(in) :: en
      complex(rp), allocatable, intent(out) :: z_contour(:)
      integer :: ie, ne

      if (.not. allocated(en%ene)) then
         call g_logger%error('reciprocal%build_green_contour: en%ene not allocated; ' // &
            'call en%e_mesh() before filling the k-space Green function.', __FILE__, __LINE__)
         return
      end if

      ne = size(en%ene)
      allocate (z_contour(ne))
      do ie = 1, ne
         z_contour(ie) = cmplx(en%ene(ie), this%green_eta, rp)
      end do
   end subroutine build_green_contour

   !> @brief Fill the canonical `green` arrays from the k-space engine.
   !> @details Route-agnostic entry point (design B2 §1.1). Ensures the energy
   !>          grid, builds the retarded contour, and dispatches to the selected
   !>          backend. Backends fill the SAME `gij/gji`, `gij_eta`, and torque
   !>          arrays the recursion route fills, so consumers are untouched.
   !>          Backend cores are staged: E in B2.2, D in B2.4 (they currently
   !>          raise a not-implemented error). `fill_green` is not yet wired into
   !>          any production path (dispatch key lands in B2.5), so this skeleton
   !>          leaves regression behaviour bit-identical.
   !> @param[inout] this      Reciprocal object (k-mesh, H(k) machinery).
   !> @param[inout] green_obj Green object whose arrays are populated in place.
   !> @param[in]    sigma     Self-energy provider (sigma_zero for backend E).
   module subroutine fill_green(this, green_obj, sigma)
      class(reciprocal), intent(inout) :: this
      type(green), intent(inout) :: green_obj
      class(sigma_provider), intent(in) :: sigma
      complex(rp), allocatable :: z_contour(:)

      ! Guarantee the real-axis grid (same call the recursion route makes).
      call green_obj%en%e_mesh()

      ! Adopt the contour convention pinned for gate G-B2-1.
      call this%build_green_contour(green_obj%en, z_contour)

      select case (trim(this%green_backend))
      case ('lehmann')
         ! Backend E (strict Lehmann, Sigma = 0). sigma is accepted for a
         ! uniform filler signature but unused here (it enters backend D only).
         call fill_green_lehmann(this, green_obj, z_contour)
      case ('dyson')
         ! sigma is the provider backend D will invert against (B2.4).
         call g_logger%error('reciprocal%fill_green: backend D (dyson) not ' // &
            'implemented yet -- lands in task B2.4.', __FILE__, __LINE__)
      case default
         call g_logger%error('reciprocal%fill_green: unknown green_backend "' // &
            trim(this%green_backend) // '" (expected lehmann|dyson).', __FILE__, __LINE__)
      end select

      if (allocated(z_contour)) deallocate (z_contour)
   end subroutine fill_green

   !> @brief Backend E core: fill gij/gji by the strict Lehmann representation.
   !> @details For every intersite pair the recursion route tracks
   !>          (`lattice%ijpair`, indexed as in `calculate_intersite_gf_core`),
   !>          accumulate the true nb x nb resolvent block directly from the
   !>          k-space eigenpairs via `lehmann_pair_block`:
   !>            G_ij(z) = (1/N_k) sum_{k,n} e^{i k.dR_ij} psi_i psi_j^dagger
   !>                                                       / (z - eps_nk),
   !>          with dR_ij = R_i - R_j the bond vector between the pair's cluster
   !>          atoms and psi_i the site-i sub-block of the eigenvector (rows
   !>          (site-1)*nb+1 .. site*nb, matching the H(k) block layout in
   !>          `reciprocal_fourier`). No 4-phase machinery: backend E fills the
   !>          intersite block directly. The torque-resolved families
   !>          (`ginmag`/`gi{x,y,z}`) and the `gij_eta` ladder are derived from
   !>          these blocks in B2.3.
   !>
   !>          Preconditions (B2.2 scope): the k-mesh is the full, unreduced BZ
   !>          and is NOT distributed over MPI ranks (each rank holds all k), and
   !>          the standard (orthonormal) eigenproblem is used so the Lehmann
   !>          completeness sum_n psi_i psi_j^dagger = delta_ij holds. Symmetry-
   !>          star averaging and k-parallel reduction attach with the B2.5
   !>          dispatch / B4 parallelization.
   !> @param[inout] this      Reciprocal object (H(k) + eigenpairs).
   !> @param[inout] green_obj Green object whose gij/gji are populated in place.
   !> @param[in]    z_contour Retarded complex-energy contour (size = ne).
   subroutine fill_green_lehmann(this, green_obj, z_contour)
      class(reciprocal), intent(inout) :: this
      type(green), intent(inout) :: green_obj
      complex(rp), intent(in) :: z_contour(:)

      integer :: nk, ne, ik, ikg
      integer :: pair, pair_glob, i_cl, j_cl, site_i, site_j, ioff, joff
      real(rp), allocatable :: kfrac(:, :)
      real(rp) :: dcart(3), dr(3)
      complex(rp), allocatable :: gblk(:, :, :)

      ! Ensure the eigenpairs exist (one Hermitian solve per k, reused for all
      ! energies and pairs).
      if (.not. allocated(this%eigenvectors) .or. .not. allocated(this%eigenvalues)) then
         call this%build_kspace_hamiltonian()
         call this%diagonalize_hamiltonian()
      end if
      if (.not. allocated(this%eigenvectors)) then
         call g_logger%error('fill_green_lehmann: eigenvectors unavailable after ' // &
            'diagonalization.', __FILE__, __LINE__)
         return
      end if

      nk = size(this%eigenvalues, 2)
      if (nk /= this%nk_total) then
         call g_logger%error('fill_green_lehmann: backend E requires an undistributed ' // &
            'k-mesh (nk_local /= nk_total); k-parallel Lehmann lands with B2.5 dispatch.', &
            __FILE__, __LINE__)
         return
      end if
      ne = size(z_contour)

      ! Fractional k-points for the (undistributed) mesh.
      allocate (kfrac(3, nk))
      do ik = 1, nk
         ikg = local_k_index_to_global(this, ik)
         kfrac(:, ik) = this%k_points(:, ikg)
      end do

      allocate (gblk(nb, nb, ne))

      do pair_glob = start_atom, end_atom
         pair = g2l_map(pair_glob)
         i_cl = this%lattice%ijpair(pair_glob, 1)
         j_cl = this%lattice%ijpair(pair_glob, 2)
         site_i = this%lattice%iz(i_cl)
         site_j = this%lattice%iz(j_cl)
         ioff = (site_i - 1)*nb
         joff = (site_j - 1)*nb

         ! dR_ij = R_i - R_j in fractional coordinates. cr is in units of alat,
         ! so absolute cartesian = cr*alat and a_cart_inv maps that to
         ! fractional -- the same table/sign as build_neighbor_vectors (C3).
         dcart = (this%lattice%cr(:, i_cl) - this%lattice%cr(:, j_cl))*this%lattice%alat
         if (this%lattice%a_cart_inv_ready) then
            dr = matmul(this%lattice%a_cart_inv, dcart)
         else
            dr = matmul(inverse_3x3(this%lattice%a), dcart/this%lattice%alat)
         end if

         ! G_ij (on-site block when site_i == site_j and dR = 0).
         call lehmann_pair_block(this%eigenvalues, this%eigenvectors, kfrac, z_contour, &
                                 dr, ioff, joff, nb, gblk)
         green_obj%gij(:, :, :, pair) = gblk

         ! G_ji: swap the site blocks and negate the bond vector (dR_ji = -dR_ij).
         call lehmann_pair_block(this%eigenvalues, this%eigenvectors, kfrac, z_contour, &
                                 -dr, joff, ioff, nb, gblk)
         green_obj%gji(:, :, :, pair) = gblk
      end do

      deallocate (kfrac, gblk)
   end subroutine fill_green_lehmann

end submodule reciprocal_green
