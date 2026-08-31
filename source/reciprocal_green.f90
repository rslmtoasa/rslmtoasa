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
!> k-space Green’s-function engine (milestone B2, flagship). One filler API,
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
!> and contour the recursion route uses, so `exchange`’s Simpson integral
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
   use dyson_kernel_mod, only: dyson_kspace_inverse
   use mpi_mod, only: start_atom, end_atom, g2l_map
   use basis_mod, only: spin_off
   implicit none

contains

   !> @brief Build the retarded complex-energy contour for the k-space filler.
   !> @details Adopts the real-space route’s energy grid verbatim (gate
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
         ! Backend D (direct Dyson): sigma is the provider it inverts against
         ! (sigma_zero gives Sigma = 0, which makes D reproduce E to solver
         ! tolerance -- the permanent CI invariant pinned in B2.4).
         call fill_green_dyson(this, green_obj, z_contour, sigma)
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
   !>          with dR_ij = R_i - R_j the bond vector between the pair’s cluster
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
   !>            wanted, BOTH the RS and the k-space Green’s functions must be
   !>            rotated (via the same `rotmag_loc` primitive `rotate_to_local_axis`
   !>            applies to the Hamiltonian) -- rotating only one side is wrong.
   !>            The intersite i/=j case (i and j moments in different directions)
   !>            is an open question and deliberately out of scope here.
   !>          * `gij_eta`/`gji_eta` ladder: the 64-point Gauss-Legendre Fermi
   !>            ladder evaluated at z = ene(fermi_point) + i*(1-x)/x, matching
   !>            `bgreen`’s eta contour (z = e(ei) + eta) and the `x,w` roots of
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

      integer :: nk, ne, ik, ikg, ie
      integer :: pair, pair_glob, ioff, joff, ipair, npair_local, ncontour
      integer :: lehmann_status
      integer :: lehmann_clock_start, lehmann_clock_stop, lehmann_clock_rate
      logical :: do_eta, use_gpu_contraction, use_gpu_resident
      real(rp), allocatable :: kfrac(:, :)
      real(rp) :: dr(3)
      complex(rp), allocatable :: gblk(:, :, :), gblk_eta(:, :, :)
      complex(rp), allocatable :: z_eta(:)
      complex(rp), allocatable :: tnmag(:, :, :), tz(:, :, :), ty(:, :, :), tx(:, :, :)
      complex(rp), allocatable :: z_lehmann(:)
      type(reciprocal_lehmann_request) :: lehmann_request
      type(reciprocal_lehmann_result) :: lehmann_result
      character(len=512) :: lehmann_timing_message

      ! Set this before the normal-mesh build so ACC-11 can retain one
      ! complete CUDA k tile for the immediate Lehmann consumer.  CPU and
      ! non-Lehmann paths ignore the request.
      this%resident_lehmann_handoff_requested = .true.

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
      call this%require_replicated_k_workset('fill_green_lehmann')

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
         call build_fermi_eta_contour(green_obj%en, z_eta)
         allocate (gblk_eta(nb, nb, 64))
         allocate (tnmag(norb, norb, 64), tz(norb, norb, 64), ty(norb, norb, 64), tx(norb, norb, 64))
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

      call system_clock(count_rate=lehmann_clock_rate)
      call system_clock(lehmann_clock_start)

      ! ACC-10 uses one all-pair backend request.  The host eigenpair arrays
      ! and canonical Green object remain the public data path; ACC-11 may
      ! avoid only the redundant CUDA eigensystem H2D copy when the persistent
      ! solver context still owns the matching normal-mesh result.
      use_gpu_contraction = .false.
      use_gpu_resident = .false.
      if (allocated(this%execution_backend)) then
         select type (backend => this%execution_backend)
         type is (cuda_reciprocal_backend)
            use_gpu_contraction = backend%initialized
            if (use_gpu_contraction) use_gpu_resident = backend%resident_eigensystem_matches(&
               size(this%eigenvalues, 1), nk, this%device_eigensystem_token)
         class default
            use_gpu_contraction = .false.
         end select
      end if
      if (use_gpu_contraction) then
         npair_local = end_atom - start_atom + 1
         ncontour = ne + merge(64, 0, do_eta)
         allocate(z_lehmann(ncontour))
         z_lehmann(1:ne) = z_contour
         if (do_eta) z_lehmann(ne + 1:ne + 64) = z_eta
         allocate(lehmann_request%k_points(3, nk), lehmann_request%z_contour(ncontour), &
                  lehmann_request%dr(3, 2*npair_local), lehmann_request%ioffset(2*npair_local), &
                  lehmann_request%joffset(2*npair_local))
         lehmann_request%nmat = size(this%eigenvalues, 1)
         lehmann_request%nk = nk
         if (.not. use_gpu_resident) then
            allocate(lehmann_request%eigenvalues(size(this%eigenvalues, 1), size(this%eigenvalues, 2)), &
                     lehmann_request%eigenvectors(size(this%eigenvectors, 1), size(this%eigenvectors, 2), size(this%eigenvectors, 3)))
            lehmann_request%eigenvalues = this%eigenvalues
            lehmann_request%eigenvectors = this%eigenvectors
         end if
         lehmann_request%k_points = kfrac
         lehmann_request%z_contour = z_lehmann
         lehmann_request%nblk = nb
         do pair_glob = start_atom, end_atom
            ipair = pair_glob - start_atom + 1
            call pair_geometry(this, pair_glob, ioff, joff, dr)
            lehmann_request%dr(:, ipair) = dr
            lehmann_request%ioffset(ipair) = ioff
            lehmann_request%joffset(ipair) = joff
            lehmann_request%dr(:, npair_local + ipair) = -dr
            lehmann_request%ioffset(npair_local + ipair) = joff
            lehmann_request%joffset(npair_local + ipair) = ioff
         end do
         if (use_gpu_resident) then
            select type (backend => this%execution_backend)
            type is (cuda_reciprocal_backend)
               call backend%contract_lehmann_resident(lehmann_request, lehmann_result, lehmann_status)
            class default
               lehmann_status = 1
            end select
         else
            call this%execution_backend%contract_lehmann(lehmann_request, lehmann_result, lehmann_status)
         end if
         if (lehmann_status /= 0 .or. .not. lehmann_result%valid) then
            call g_logger%error('fill_green_lehmann: CUDA Lehmann contraction failed; canonical arrays were not updated.', &
               __FILE__, __LINE__)
            return
         end if
      end if

      do pair_glob = start_atom, end_atom
         pair = g2l_map(pair_glob)
         ipair = pair_glob - start_atom + 1
         call pair_geometry(this, pair_glob, ioff, joff, dr)

         ! G_ij (on-site block when site_i == site_j and dR = 0). Delivered in the
         ! GLOBAL spin frame -- the same frame recur_b_ij stores the RS gij in.
         if (use_gpu_contraction) then
            gblk = lehmann_result%blocks(:, :, 1:ne, ipair)
         else
            call lehmann_pair_block(this%eigenvalues, this%eigenvectors, kfrac, z_contour, &
                                    dr, ioff, joff, nb, gblk)
         end if
         green_obj%gij(:, :, :, pair) = gblk

         ! G_ji: swap the site blocks and negate the bond vector (dR_ji = -dR_ij).
         if (use_gpu_contraction) then
            gblk = lehmann_result%blocks(:, :, 1:ne, npair_local + ipair)
         else
            call lehmann_pair_block(this%eigenvalues, this%eigenvectors, kfrac, z_contour, &
                                    -dr, joff, ioff, nb, gblk)
         end if
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
            if (use_gpu_contraction) then
               gblk_eta = lehmann_result%blocks(:, :, ne + 1:ne + 64, ipair)
            else
               call lehmann_pair_block(this%eigenvalues, this%eigenvectors, kfrac, z_eta, &
                                       dr, ioff, joff, nb, gblk_eta)
            end if
            do ie = 1, 64
               green_obj%gij_eta(ie, :, :, pair) = gblk_eta(:, :, ie)
            end do
            if (allocated(green_obj%ginmag_eta)) then
               call pauli_decompose_block(gblk_eta, spin_off, tnmag, tz, ty, tx)
               call store_eta_torque(tnmag, tz, ty, tx, green_obj%ginmag_eta(:, :, :, pair), &
                  green_obj%giz_eta(:, :, :, pair), green_obj%giy_eta(:, :, :, pair), &
                  green_obj%gix_eta(:, :, :, pair))
            end if

            if (use_gpu_contraction) then
               gblk_eta = lehmann_result%blocks(:, :, ne + 1:ne + 64, npair_local + ipair)
            else
               call lehmann_pair_block(this%eigenvalues, this%eigenvectors, kfrac, z_eta, &
                                       -dr, joff, ioff, nb, gblk_eta)
            end if
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

      call system_clock(lehmann_clock_stop)
      if (use_gpu_contraction) then
         write(lehmann_timing_message, '(a,es16.8,a,es16.8,a,es16.8,a,es16.8)') &
            'ACC10_TIMING backend=cuda resident='//merge('1', '0', use_gpu_resident)//' total_seconds=', &
            real(lehmann_clock_stop - lehmann_clock_start, rp) / real(max(1, lehmann_clock_rate), rp), &
            ' h2d_seconds=', lehmann_result%h2d_seconds, &
            ' contraction_seconds=', lehmann_result%contraction_seconds, &
            ' d2h_seconds=', lehmann_result%d2h_seconds
      else
         write(lehmann_timing_message, '(a,es16.8,a,es16.8,a,es16.8,a,es16.8)') &
            'ACC10_TIMING backend=lapack total_seconds=', &
            real(lehmann_clock_stop - lehmann_clock_start, rp) / real(max(1, lehmann_clock_rate), rp), &
            ' h2d_seconds=', 0.0_rp, ' contraction_seconds=', &
            real(lehmann_clock_stop - lehmann_clock_start, rp) / real(max(1, lehmann_clock_rate), rp), &
            ' d2h_seconds=', 0.0_rp
      end if
      call g_logger%info(trim(lehmann_timing_message), __FILE__, __LINE__)

      if (allocated(z_eta)) deallocate (z_eta)
      if (allocated(z_lehmann)) deallocate (z_lehmann)
      if (allocated(gblk_eta)) deallocate (gblk_eta)
      if (allocated(tnmag)) deallocate (tnmag, tz, ty, tx)
      if (allocated(lehmann_result%blocks)) deallocate (lehmann_result%blocks)
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

   !> @brief Pair -> site offsets and fractional bond vector dR_ij = R_i - R_j.
   !> @details Shared by both backends (E and D): maps a global pair index to the
   !>          site row offsets in the H(k)/eigenvector layout and the fractional
   !>          bond vector used for the e^{i k.dR} inverse Bloch phase (C3). cr is
   !>          in units of alat, so absolute cartesian = cr*alat and a_cart_inv
   !>          maps that to fractional -- the same table/sign as
   !>          build_neighbor_vectors.
   !> @param[in]  this      Reciprocal object (lattice pair/site maps).
   !> @param[in]  pair_glob Global pair index (start_atom..end_atom).
   !> @param[out] ioff      Zero-based row offset of site i’s block ((site_i-1)*nb).
   !> @param[out] joff      Zero-based row offset of site j’s block ((site_j-1)*nb).
   !> @param[out] dr        Bond vector R_i - R_j in fractional coordinates.
   subroutine pair_geometry(this, pair_glob, ioff, joff, dr)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: pair_glob
      integer, intent(out) :: ioff, joff
      real(rp), intent(out) :: dr(3)
      integer :: i_cl, j_cl
      real(rp) :: dcart(3)

      i_cl = this%lattice%ijpair(pair_glob, 1)
      j_cl = this%lattice%ijpair(pair_glob, 2)
      ioff = (this%lattice%iz(i_cl) - 1)*nb
      joff = (this%lattice%iz(j_cl) - 1)*nb

      dcart = (this%lattice%cr(:, i_cl) - this%lattice%cr(:, j_cl))*this%lattice%alat
      if (this%lattice%a_cart_inv_ready) then
         dr = matmul(this%lattice%a_cart_inv, dcart)
      else
         dr = matmul(inverse_3x3(this%lattice%a), dcart/this%lattice%alat)
      end if
   end subroutine pair_geometry

   !> @brief Build the 64-point Fermi-eta ladder both backends share.
   !> @details The recursion route’s Fermi-point continued-fraction ladder,
   !>          evaluated here at z = ene(fermi_point) + i*(1-x)/x with x the
   !>          gauss_legendre(64, 0, 1) roots -- the same `bgreen` eta contour
   !>          (z = e(ei) + eta). `fermi_point` is the last grid index at or below
   !>          en%fermi (matching the recursion filler). The Gauss weights belong
   !>          to the consumer, not the fill, so only the nodes are returned.
   !> @param[in]  en    Energy object (prepared real-axis grid + fermi).
   !> @param[out] z_eta Allocated eta contour, shape (64).
   subroutine build_fermi_eta_contour(en, z_eta)
      class(energy), intent(in) :: en
      complex(rp), allocatable, intent(out) :: z_eta(:)
      real(rp) :: xgl(64), wgl(64)
      integer :: i, fermi_point

      call gauss_legendre(64, 0.0_rp, 1.0_rp, xgl, wgl)
      fermi_point = 1
      do i = 1, size(en%ene)
         if ((en%ene(i) - en%fermi) .le. 0.000001_rp) fermi_point = i
      end do
      allocate (z_eta(64))
      do i = 1, 64
         z_eta(i) = cmplx(en%ene(fermi_point), (1.0_rp - xgl(i))/xgl(i), rp)
      end do
   end subroutine build_fermi_eta_contour

   !> @brief Assemble the block-diagonal self-energy Sigma(z) for backend D.
   !> @details Places each unit-cell site’s nb x nb retarded self-energy block on
   !>          the diagonal of the full nmat x nmat Dyson matrix, in the same site
   !>          order as the H(k) layout. Sigma depends only on z (not k), so this
   !>          is recomputed per energy; for sigma_zero every block is zero.
   !> @param[in]  sigma    Self-energy provider.
   !> @param[in]  z        Complex energy on the retarded contour.
   !> @param[in]  nsite    Number of unit-cell sites (nmat/nb).
   !> @param[inout] sblk   nb x nb scratch block (avoids per-call allocation).
   !> @param[out] sig_full Assembled nmat x nmat block-diagonal Sigma(z).
   subroutine build_sigma_full(sigma, z, nsite, sblk, sig_full)
      class(sigma_provider), intent(in) :: sigma
      complex(rp), intent(in) :: z
      integer, intent(in) :: nsite
      complex(rp), intent(inout) :: sblk(:, :)
      complex(rp), intent(out) :: sig_full(:, :)
      integer :: is, off

      sig_full = (0.0_rp, 0.0_rp)
      do is = 1, nsite
         call sigma%get_sigma(z, is, sblk)
         off = (is - 1)*nb
         sig_full(off + 1:off + nb, off + 1:off + nb) = sblk
      end do
   end subroutine build_sigma_full

   !> @brief Backend D core: direct-Dyson fill of gij/gji, eta ladder, torque.
   !> @details Streams over (k, z): one nmat x nmat inversion
   !>            G(k,z) = [z*I - H(k) - Sigma(z)]^{-1}
   !>          per (k, z) (S = I, pinned in dyson_kernel_mod -- backend E uses the
   !>          orthonormal eigenproblem, so the D == E invariant holds only for
   !>          S = I), distributed to EVERY pair with the same 1/N_k inverse Bloch
   !>          phase e^{i k.dR_ij} and pair->site map as backend E (design B2 §1.4:
   !>          never re-invert per pair or materialize all (k,z)). Sigma enters via
   !>          the provider; sigma_zero makes this reproduce backend E to solver
   !>          tolerance -- the permanent CI invariant (test_dyson_equivalence).
   !>          The torque families and the eta ladder are filled from the same
   !>          Pauli algebra (pauli_decompose_block) backend E uses.
   !> @param[inout] this      Reciprocal object (H(k) tiles).
   !> @param[inout] green_obj Green object whose gij/gji are populated in place.
   !> @param[in]    z_contour Retarded complex-energy contour (size = ne).
   !> @param[in]    sigma     Self-energy provider (sigma_zero => backend E).
   subroutine fill_green_dyson(this, green_obj, z_contour, sigma)
      class(reciprocal), intent(inout) :: this
      type(green), intent(inout) :: green_obj
      complex(rp), intent(in) :: z_contour(:)
      class(sigma_provider), intent(in) :: sigma

      integer :: nk, ne, nmat, nsite, ik, ikg, ie, i, npair_loc, pair, pair_glob, ioff, joff
      logical :: do_eta
      real(rp), allocatable :: kfrac(:, :), p_dr(:, :)
      integer, allocatable :: p_ioff(:), p_joff(:), p_pair(:)
      real(rp) :: kdotr
      complex(rp), allocatable :: phse(:), z_eta(:)
      complex(rp), allocatable :: gk(:, :), sig_full(:, :), sblk(:, :), gtmp(:, :, :)
      complex(rp), allocatable :: tnmag(:, :, :), tz(:, :, :), ty(:, :, :), tx(:, :, :)

      ! Backend D inverts H(k) directly -- no eigensolve. Build the tiles if the
      ! Lehmann route has not already populated them.
      if (.not. allocated(this%hk_bulk)) call this%build_kspace_hamiltonian()
      if (.not. allocated(this%hk_bulk)) then
         call g_logger%error('fill_green_dyson: H(k) unavailable after build.', __FILE__, __LINE__)
         return
      end if
      call this%require_replicated_k_workset('fill_green_dyson')

      nk = size(this%hk_bulk, 3)
      if (nk /= this%nk_total) then
         call g_logger%error('fill_green_dyson: backend D requires an undistributed ' // &
            'k-mesh (nk_local /= nk_total); k-parallel Dyson lands with B2.5 dispatch.', &
            __FILE__, __LINE__)
         return
      end if
      nmat = size(this%hk_bulk, 1)
      nsite = nmat/nb
      ne = size(z_contour)

      ! Fractional k-points for the (undistributed) mesh.
      allocate (kfrac(3, nk))
      do ik = 1, nk
         ikg = local_k_index_to_global(this, ik)
         kfrac(:, ik) = this%k_points(:, ikg)
      end do

      ! Per-pair geometry, cached once (dR/offsets are k- and z-independent).
      npair_loc = end_atom - start_atom + 1
      allocate (p_ioff(npair_loc), p_joff(npair_loc), p_pair(npair_loc), p_dr(3, npair_loc), phse(npair_loc))
      i = 0
      do pair_glob = start_atom, end_atom
         i = i + 1
         p_pair(i) = g2l_map(pair_glob)
         call pair_geometry(this, pair_glob, p_ioff(i), p_joff(i), p_dr(:, i))
      end do

      allocate (gk(nmat, nmat), sig_full(nmat, nmat), sblk(nb, nb))

      ! The green arrays double as (k,z)-accumulators; zero, accumulate, then
      ! divide by N_k. gij/gji cover every local pair, so this is complete.
      green_obj%gij = (0.0_rp, 0.0_rp); green_obj%gji = (0.0_rp, 0.0_rp)
      do_eta = allocated(green_obj%gij_eta) .and. allocated(green_obj%gji_eta)
      if (do_eta) then
         call build_fermi_eta_contour(green_obj%en, z_eta)
         green_obj%gij_eta = (0.0_rp, 0.0_rp); green_obj%gji_eta = (0.0_rp, 0.0_rp)
      end if

      do ik = 1, nk
         ! Bond phases e^{i 2pi k.dR_ij} for this k (depend on k and pair only).
         do i = 1, npair_loc
            kdotr = two_pi*dot_product(kfrac(:, ik), p_dr(:, i))
            phse(i) = cmplx(cos(kdotr), sin(kdotr), rp)
         end do

         do ie = 1, ne
            call build_sigma_full(sigma, z_contour(ie), nsite, sblk, sig_full)
            call dyson_kspace_inverse(this%hk_bulk(:, :, ik), z_contour(ie), sig_full, gk)
            do i = 1, npair_loc
               ioff = p_ioff(i); joff = p_joff(i); pair = p_pair(i)
               green_obj%gij(:, :, ie, pair) = green_obj%gij(:, :, ie, pair) + &
                  phse(i)*gk(ioff + 1:ioff + nb, joff + 1:joff + nb)
               green_obj%gji(:, :, ie, pair) = green_obj%gji(:, :, ie, pair) + &
                  conjg(phse(i))*gk(joff + 1:joff + nb, ioff + 1:ioff + nb)
            end do
         end do

         if (do_eta) then
            do ie = 1, 64
               call build_sigma_full(sigma, z_eta(ie), nsite, sblk, sig_full)
               call dyson_kspace_inverse(this%hk_bulk(:, :, ik), z_eta(ie), sig_full, gk)
               do i = 1, npair_loc
                  ioff = p_ioff(i); joff = p_joff(i); pair = p_pair(i)
                  green_obj%gij_eta(ie, :, :, pair) = green_obj%gij_eta(ie, :, :, pair) + &
                     phse(i)*gk(ioff + 1:ioff + nb, joff + 1:joff + nb)
                  green_obj%gji_eta(ie, :, :, pair) = green_obj%gji_eta(ie, :, :, pair) + &
                     conjg(phse(i))*gk(joff + 1:joff + nb, ioff + 1:ioff + nb)
               end do
            end do
         end if
      end do

      ! 1/N_k normalization (inverse Bloch transform), same factor as backend E.
      green_obj%gij = green_obj%gij/real(nk, rp)
      green_obj%gji = green_obj%gji/real(nk, rp)
      if (do_eta) then
         green_obj%gij_eta = green_obj%gij_eta/real(nk, rp)
         green_obj%gji_eta = green_obj%gji_eta/real(nk, rp)
      end if

      ! Torque families: identical Pauli algebra as backend E, on the accumulated
      ! blocks. Zero up front so pairs on other ranks stay zero (MPI-safe).
      if (allocated(green_obj%ginmag)) then
         green_obj%ginmag = (0.0_rp, 0.0_rp); green_obj%gjnmag = (0.0_rp, 0.0_rp)
         green_obj%gix = (0.0_rp, 0.0_rp); green_obj%giy = (0.0_rp, 0.0_rp)
         green_obj%giz = (0.0_rp, 0.0_rp); green_obj%gjx = (0.0_rp, 0.0_rp)
         green_obj%gjy = (0.0_rp, 0.0_rp); green_obj%gjz = (0.0_rp, 0.0_rp)
      end if
      if (do_eta .and. allocated(green_obj%ginmag_eta)) then
         green_obj%ginmag_eta = (0.0_rp, 0.0_rp); green_obj%gjnmag_eta = (0.0_rp, 0.0_rp)
         green_obj%gix_eta = (0.0_rp, 0.0_rp); green_obj%giy_eta = (0.0_rp, 0.0_rp)
         green_obj%giz_eta = (0.0_rp, 0.0_rp); green_obj%gjx_eta = (0.0_rp, 0.0_rp)
         green_obj%gjy_eta = (0.0_rp, 0.0_rp); green_obj%gjz_eta = (0.0_rp, 0.0_rp)
         allocate (gtmp(nb, nb, 64))
         allocate (tnmag(norb, norb, 64), tz(norb, norb, 64), ty(norb, norb, 64), tx(norb, norb, 64))
      end if

      do i = 1, npair_loc
         pair = p_pair(i)
         if (allocated(green_obj%ginmag)) then
            call pauli_decompose_block(green_obj%gij(:, :, :, pair), spin_off, &
               green_obj%ginmag(:, :, :, pair), green_obj%giz(:, :, :, pair), &
               green_obj%giy(:, :, :, pair), green_obj%gix(:, :, :, pair))
            call pauli_decompose_block(green_obj%gji(:, :, :, pair), spin_off, &
               green_obj%gjnmag(:, :, :, pair), green_obj%gjz(:, :, :, pair), &
               green_obj%gjy(:, :, :, pair), green_obj%gjx(:, :, :, pair))
         end if
         if (do_eta .and. allocated(green_obj%ginmag_eta)) then
            do ie = 1, 64
               gtmp(:, :, ie) = green_obj%gij_eta(ie, :, :, pair)
            end do
            call pauli_decompose_block(gtmp, spin_off, tnmag, tz, ty, tx)
            call store_eta_torque(tnmag, tz, ty, tx, green_obj%ginmag_eta(:, :, :, pair), &
               green_obj%giz_eta(:, :, :, pair), green_obj%giy_eta(:, :, :, pair), &
               green_obj%gix_eta(:, :, :, pair))
            do ie = 1, 64
               gtmp(:, :, ie) = green_obj%gji_eta(ie, :, :, pair)
            end do
            call pauli_decompose_block(gtmp, spin_off, tnmag, tz, ty, tx)
            call store_eta_torque(tnmag, tz, ty, tx, green_obj%gjnmag_eta(:, :, :, pair), &
               green_obj%gjz_eta(:, :, :, pair), green_obj%gjy_eta(:, :, :, pair), &
               green_obj%gjx_eta(:, :, :, pair))
         end if
      end do

      if (allocated(gtmp)) deallocate (gtmp, tnmag, tz, ty, tx)
      if (allocated(z_eta)) deallocate (z_eta)
      deallocate (kfrac, p_ioff, p_joff, p_pair, p_dr, phse, gk, sig_full, sblk)
   end subroutine fill_green_dyson

   ! NOTE (spin frame / future local-frame option): backend E delivers gij/gji in
   ! the GLOBAL spin frame, matching the RS intersite blocks -- the intersite
   ! recursion `recur_b_ij` never rotates to local axes (only the on-site DOS
   ! recursion `recursion_haydock` does). If a LOCAL-frame comparison is ever
   ! wanted, rotate BOTH routes’ Green’s functions with the same primitive
   ! `rotate_to_local_axis` uses on the Hamiltonian, e.g. per block:
   !   call rotmag_loc(blk_rotated, blk_global, size(blk, 3), mom)   ! = R^dagger B R
   ! with mom the block’s central-atom moment. Rotating only the k-space side is
   ! wrong: an early Mn3Sn (120 deg NC) experiment that did so broke the m_z match
   ! (~4e-4 -> ~20). The intersite i/=j case (i, j moments differ) is open.

end submodule reciprocal_green
