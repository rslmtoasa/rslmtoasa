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
   use lehmann_kernel_mod, only: lehmann_pair_block, pauli_decompose_block
   use mpi_mod, only: start_atom, end_atom, g2l_map
   use basis_mod, only: spin_off
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

   !> @brief Backend E core: fill gij/gji, the eta ladder, and torque families.
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
   !>          intersite block directly.
   !>
   !>          B2.3 additions (all mirroring `calculate_intersite_gf_core` so the
   !>          consumers -- exchange, damping -- run unchanged):
   !>          * Spin frame (C2). Backend E fills the SAME frame the recursion
   !>            route stores the intersite blocks in, which is the GLOBAL spin
   !>            frame: the intersite recursion `recur_b_ij` (recursion_transport)
   !>            never rotates to local axes (`local_axis` there only gates GPU
   !>            eligibility -- ONLY the on-site DOS recursion `recursion_haydock`
   !>            rotates). Confirmed on a genuine noncollinear background (Mn3Sn,
   !>            120 deg): the on-site z-projected spin DOS m_z from both routes
   !>            agrees to ~4e-4 with NO rotation, and an early experiment that
   !>            rotated only the k-space block into the local frame broke the
   !>            match (m_z diff -> ~20). So backend E delivers the global-frame
   !>            block directly. NOTE (future): if a local-frame comparison is ever
   !>            wanted, BOTH the RS and the k-space Green's functions must be
   !>            rotated (via the same `rotmag_loc` primitive `rotate_to_local_axis`
   !>            applies to the Hamiltonian) -- rotating only one side is wrong.
   !>            The intersite i/=j case (i and j moments in different directions)
   !>            is an open question and deliberately out of scope here.
   !>          * `gij_eta`/`gji_eta` ladder: the 64-point Gauss-Legendre Fermi
   !>            ladder evaluated at z = ene(fermi_point) + i*(1-x)/x, matching
   !>            `bgreen`'s eta contour (z = e(ei) + eta) and the `x,w` roots of
   !>            `gauss_legendre(64, 0, 1)`. Stored with the eta index leading,
   !>            `(64, nb, nb, pair)`, as the recursion route stores it.
   !>          * Torque-resolved families `ginmag`/`gi{x,y,z}` (+ `gj*` and the
   !>            `*_eta` partners) derived from the filled blocks by the exact
   !>            Pauli spin-block algebra of `calculate_intersite_gf_core`.
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

      integer :: nk, ne, ik, ikg, ie, i, j, fermi_point
      integer :: pair, pair_glob, i_cl, j_cl, site_i, site_j, ioff, joff
      logical :: do_eta
      real(rp), allocatable :: kfrac(:, :)
      real(rp) :: dcart(3), dr(3)
      real(rp) :: xgl(64), wgl(64)
      complex(rp), allocatable :: gblk(:, :, :), gblk_eta(:, :, :)
      complex(rp), allocatable :: z_eta(:)
      complex(rp), allocatable :: tnmag(:, :, :), tz(:, :, :), ty(:, :, :), tx(:, :, :)

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

      ! Eta ladder: the Fermi-point continued-fraction ladder of the recursion
      ! route, here evaluated as strict-Lehmann blocks at 64 broadenings. Filled
      ! only when the target arrays exist (they are needed by the noncollinear
      ! exchange/damping consumers via green%gij_eta_to_gij).
      do_eta = allocated(green_obj%gij_eta) .and. allocated(green_obj%gji_eta)
      if (do_eta) then
         call gauss_legendre(64, 0.0_rp, 1.0_rp, xgl, wgl)
         fermi_point = 1
         do i = 1, size(green_obj%en%ene)
            if ((green_obj%en%ene(i) - green_obj%en%fermi) .le. 0.000001_rp) fermi_point = i
         end do
         allocate (z_eta(64), gblk_eta(nb, nb, 64))
         allocate (tnmag(norb, norb, 64), tz(norb, norb, 64), ty(norb, norb, 64), tx(norb, norb, 64))
         do i = 1, 64
            z_eta(i) = cmplx(green_obj%en%ene(fermi_point), (1.0_rp - xgl(i))/xgl(i), rp)
         end do
         green_obj%gij_eta = (0.0_rp, 0.0_rp); green_obj%gji_eta = (0.0_rp, 0.0_rp)
         if (allocated(green_obj%ginmag_eta)) then
            green_obj%ginmag_eta = (0.0_rp, 0.0_rp); green_obj%gjnmag_eta = (0.0_rp, 0.0_rp)
            green_obj%gix_eta = (0.0_rp, 0.0_rp); green_obj%giy_eta = (0.0_rp, 0.0_rp)
            green_obj%giz_eta = (0.0_rp, 0.0_rp); green_obj%gjx_eta = (0.0_rp, 0.0_rp)
            green_obj%gjy_eta = (0.0_rp, 0.0_rp); green_obj%gjz_eta = (0.0_rp, 0.0_rp)
         end if
      end if

      if (allocated(green_obj%ginmag)) then
         green_obj%ginmag = (0.0_rp, 0.0_rp); green_obj%gjnmag = (0.0_rp, 0.0_rp)
         green_obj%gix = (0.0_rp, 0.0_rp); green_obj%giy = (0.0_rp, 0.0_rp)
         green_obj%giz = (0.0_rp, 0.0_rp); green_obj%gjx = (0.0_rp, 0.0_rp)
         green_obj%gjy = (0.0_rp, 0.0_rp); green_obj%gjz = (0.0_rp, 0.0_rp)
      end if

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

         ! G_ij (on-site block when site_i == site_j and dR = 0). Delivered in the
         ! GLOBAL spin frame -- the same frame recur_b_ij stores the RS gij in.
         call lehmann_pair_block(this%eigenvalues, this%eigenvectors, kfrac, z_contour, &
                                 dr, ioff, joff, nb, gblk)
         green_obj%gij(:, :, :, pair) = gblk

         ! G_ji: swap the site blocks and negate the bond vector (dR_ji = -dR_ij).
         call lehmann_pair_block(this%eigenvalues, this%eigenvectors, kfrac, z_contour, &
                                 -dr, joff, ioff, nb, gblk)
         green_obj%gji(:, :, :, pair) = gblk

         ! Non-eta torque families from the (rotated) blocks -- the exact Pauli
         ! spin-block decomposition of calculate_intersite_gf_core.
         if (allocated(green_obj%ginmag)) then
            call pauli_decompose_block(green_obj%gij(:, :, :, pair), spin_off, &
               green_obj%ginmag(:, :, :, pair), green_obj%giz(:, :, :, pair), &
               green_obj%giy(:, :, :, pair), green_obj%gix(:, :, :, pair))
            call pauli_decompose_block(green_obj%gji(:, :, :, pair), spin_off, &
               green_obj%gjnmag(:, :, :, pair), green_obj%gjz(:, :, :, pair), &
               green_obj%gjy(:, :, :, pair), green_obj%gjx(:, :, :, pair))
         end if

         if (do_eta) then
            call lehmann_pair_block(this%eigenvalues, this%eigenvectors, kfrac, z_eta, &
                                    dr, ioff, joff, nb, gblk_eta)
            do ie = 1, 64
               green_obj%gij_eta(ie, :, :, pair) = gblk_eta(:, :, ie)
            end do
            if (allocated(green_obj%ginmag_eta)) then
               call pauli_decompose_block(gblk_eta, spin_off, tnmag, tz, ty, tx)
               call store_eta_torque(tnmag, tz, ty, tx, green_obj%ginmag_eta(:, :, :, pair), &
                  green_obj%giz_eta(:, :, :, pair), green_obj%giy_eta(:, :, :, pair), &
                  green_obj%gix_eta(:, :, :, pair))
            end if

            call lehmann_pair_block(this%eigenvalues, this%eigenvectors, kfrac, z_eta, &
                                    -dr, joff, ioff, nb, gblk_eta)
            do ie = 1, 64
               green_obj%gji_eta(ie, :, :, pair) = gblk_eta(:, :, ie)
            end do
            if (allocated(green_obj%ginmag_eta)) then
               call pauli_decompose_block(gblk_eta, spin_off, tnmag, tz, ty, tx)
               call store_eta_torque(tnmag, tz, ty, tx, green_obj%gjnmag_eta(:, :, :, pair), &
                  green_obj%gjz_eta(:, :, :, pair), green_obj%gjy_eta(:, :, :, pair), &
                  green_obj%gjx_eta(:, :, :, pair))
            end if
         end if
      end do

      if (allocated(z_eta)) deallocate (z_eta)
      if (allocated(gblk_eta)) deallocate (gblk_eta)
      if (allocated(tnmag)) deallocate (tnmag, tz, ty, tx)
      deallocate (kfrac, gblk)
   end subroutine fill_green_lehmann

   !> @brief Store energy-3rd torque sub-blocks into the eta-leading arrays.
   !> @details The green object stores the Fermi-eta ladder with the eta index
   !>          FIRST (`(64, norb, norb, pair)`), whereas `pauli_decompose_block`
   !>          produces the standard energy-3rd layout `(norb, norb, 64)`. This
   !>          transposes the eta axis to the front for each component.
   subroutine store_eta_torque(tnmag, tz, ty, tx, gnmag_eta, gz_eta, gy_eta, gx_eta)
      complex(rp), intent(in) :: tnmag(:, :, :), tz(:, :, :), ty(:, :, :), tx(:, :, :)
      complex(rp), intent(out) :: gnmag_eta(:, :, :), gz_eta(:, :, :), gy_eta(:, :, :), gx_eta(:, :, :)
      integer :: ie, ne

      ne = size(tnmag, 3)
      do ie = 1, ne
         gnmag_eta(ie, :, :) = tnmag(:, :, ie)
         gz_eta(ie, :, :) = tz(:, :, ie)
         gy_eta(ie, :, :) = ty(:, :, ie)
         gx_eta(ie, :, :) = tx(:, :, ie)
      end do
   end subroutine store_eta_torque

   ! NOTE (spin frame / future local-frame option): backend E delivers gij/gji in
   ! the GLOBAL spin frame, matching the RS intersite blocks -- the intersite
   ! recursion `recur_b_ij` never rotates to local axes (only the on-site DOS
   ! recursion `recursion_haydock` does). If a LOCAL-frame comparison is ever
   ! wanted, rotate BOTH routes' Green's functions with the same primitive
   ! `rotate_to_local_axis` uses on the Hamiltonian, e.g. per block:
   !   call rotmag_loc(blk_rotated, blk_global, size(blk, 3), mom)   ! = R^dagger B R
   ! with mom the block's central-atom moment. Rotating only the k-space side is
   ! wrong: an early Mn3Sn (120 deg NC) experiment that did so broke the m_z match
   ! (~4e-4 -> ~20). The intersite i/=j case (i, j moments differ) is open.

end submodule reciprocal_green
