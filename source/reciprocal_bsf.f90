!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! SUBMODULE: reciprocal_bsf
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> Bloch spectral functions A(k,E) (milestone B3). A thin consumer of the B2
!> k-space Green's-function engine: along a high-symmetry k-path it forms the
!> retarded resolvent G(k, E+i*eta) with the SAME `dyson_kspace_inverse` primitive
!> backend D uses (Sigma=0 here, so it equals backend E), and reports the
!> partial-trace spectral weight
!>
!>   A(k,E) = -(1/pi) Im Tr_sub G(k, E+i*eta)
!>
!> for the total, spin-up and spin-down channels (`bsf_spectral_trace`). Because
!> the resolvent comes from the B2 engine, the identical code path yields the
!> disorder/correlation-broadened A(k,E) once a non-zero Sigma provider is wired
!> (CPA/DMFT, B8/B10) -- that is the design point recorded in the B3 plan.
!>
!> Conventions (B3.3 review):
!>   * broadening eta = `this%green_eta` (the &reciprocal namelist key; the same
!>     retarded broadening the B2 contour uses), never hard-coded.
!>   * k-path = the canonical high-symmetry path of `calculate_band_structure`
!>     (spglib), so a band-structure overlay is directly comparable; the path
!>     eigenvalues are written alongside for exactly that overlay / the delta-limit
!>     known-answer test.
!>   * spin layout: per site the H(k) block is [up(1..norb), down(norb+1..nb)];
!>     spin-up/down partial traces select those diagonal ranges across all sites.
!>   * A(k,E) in states / Ry / spin-channel; the writer states units in its header.
!------------------------------------------------------------------------------
submodule(reciprocal_mod) reciprocal_bsf
   use bsf_kernel_mod, only: bsf_spectral_trace
   use dyson_kernel_mod, only: dyson_kspace_inverse
   implicit none

contains

   !> @brief Compute and write the Bloch spectral function A(k,E) along the path.
   !> @details Sets up the canonical high-symmetry k-path, builds and diagonalizes
   !>          H(k) on it (the eigenvalues are kept for the band overlay), then for
   !>          every (k, E) inverts z*I - H(k) via `dyson_kspace_inverse` (Sigma=0)
   !>          and accumulates the total / spin-up / spin-down spectral weight. The
   !>          resolvent is delivered by the B2 engine's primitive, so a non-zero
   !>          Sigma later broadens A(k,E) through this same routine unchanged.
   !> @param[inout] this        Reciprocal object (H(k) machinery, k-path, grids).
   !> @param[in]    output_file Optional base output filename (default 'bsf.dat').
   module subroutine calculate_bsf(this, output_file)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in), optional :: output_file

      character(len=256) :: fname, bfile
      integer :: unit, ik, ie, nmat, nsite, ne, i, s, o, dot
      real(rp) :: eta, a_tot, a_up, a_dn
      complex(rp) :: z
      complex(rp), allocatable :: gk(:, :), sig0(:, :)
      integer, allocatable :: idx_tot(:), idx_up(:), idx_dn(:)
      logical :: spin_split

      ! --- k-path (canonical high-symmetry, same as calculate_band_structure). ---
      call this%symmetry_analysis%generate_canonical_kpath(this%nk_per_segment)
      if (.not. allocated(this%symmetry_analysis%k_path)) then
         call g_logger%error('reciprocal%calculate_bsf: canonical k-path not generated.', &
            __FILE__, __LINE__)
         return
      end if
      this%nk_path = this%symmetry_analysis%nk_path
      if (allocated(this%k_path)) deallocate(this%k_path)
      if (allocated(this%k_labels)) deallocate(this%k_labels)
      if (allocated(this%k_distances)) deallocate(this%k_distances)
      allocate(this%k_path(3, this%nk_path))
      allocate(this%k_labels(this%nk_path))
      allocate(this%k_distances(this%nk_path))
      this%k_path = this%symmetry_analysis%k_path
      this%k_labels = this%symmetry_analysis%k_labels
      this%k_distances = this%symmetry_analysis%k_distances

      ! --- H(k) on the path + eigenvalues (kept for the band overlay / delta test).
      call this%build_kspace_hamiltonian()
      call this%diagonalize_hamiltonian()

      ! --- Energy grid (states axis) + retarded broadening. -----------------------
      call this%setup_dos_energy_grid()
      ne = this%n_energy_points
      eta = this%green_eta

      nmat = size(this%hk_bulk, 1)
      nsite = nmat/nb

      ! Diagonal channel index sets: total, and (when the block factors into
      ! nb=2*norb spin halves per site) spin-up / spin-down across all sites.
      allocate(idx_tot(nmat))
      do i = 1, nmat
         idx_tot(i) = i
      end do
      spin_split = (nsite*nb == nmat) .and. (nb == 2*norb)
      if (spin_split) then
         allocate(idx_up(nsite*norb), idx_dn(nsite*norb))
         i = 0
         do s = 1, nsite
            do o = 1, norb
               i = i + 1
               idx_up(i) = (s - 1)*nb + o
               idx_dn(i) = (s - 1)*nb + norb + o
            end do
         end do
      end if

      allocate(gk(nmat, nmat), sig0(nmat, nmat))
      sig0 = (0.0_rp, 0.0_rp)

      ! --- Write A(k,E) (rank 0). gnuplot-friendly: blank line between k-blocks. --
      if (rank == 0) then
         fname = 'bsf.dat'
         if (present(output_file)) fname = output_file
         open(newunit=unit, file=trim(fname), status='replace', action='write')
         write(unit, '(A)') '# Bloch spectral function A(k,E) = -1/pi Im Tr G(k,E+i*eta)'
         write(unit, '(A,ES12.5,A)') '# eta = ', eta, ' Ry (green_eta)   units: states/Ry/spin'
         write(unit, '(A,I0,A,I0)') '# k-path points = ', this%nk_path, '   energy points = ', ne
         write(unit, '(A)') '# columns: k_distance   E(Ry)   A_total   A_up   A_down'
         write(unit, '(A)') '#'
      end if

      do ik = 1, this%nk_path
         do ie = 1, ne
            z = cmplx(this%dos_energy_grid(ie), eta, rp)
            call dyson_kspace_inverse(this%hk_bulk(:, :, ik), z, sig0, gk)
            a_tot = bsf_spectral_trace(gk, idx_tot)
            if (spin_split) then
               a_up = bsf_spectral_trace(gk, idx_up)
               a_dn = bsf_spectral_trace(gk, idx_dn)
            else
               a_up = a_tot; a_dn = 0.0_rp
            end if
            if (rank == 0) write(unit, '(5(ES16.8E3,1X))') &
               this%k_distances(ik), this%dos_energy_grid(ie), a_tot, a_up, a_dn
         end do
         if (rank == 0) write(unit, '(A)') ''   ! blank line: gnuplot k-block separator
      end do

      if (rank == 0) then
         close(unit)
         call g_logger%info('reciprocal%calculate_bsf: A(k,E) written to ' // trim(fname), &
            __FILE__, __LINE__)

         ! Band eigenvalues along the same path, for direct overlay / delta-limit.
         bfile = 'bsf_bands.dat'
         if (present(output_file)) then
            dot = index(fname, '.', back=.true.)
            if (dot > 0) then
               bfile = fname(1:dot - 1) // '_bands.dat'
            else
               bfile = trim(fname) // '_bands.dat'
            end if
         end if
         open(newunit=unit, file=trim(bfile), status='replace', action='write')
         write(unit, '(A)') '# Band eigenvalues on the BSF k-path (overlay reference)'
         write(unit, '(A)') '# columns: k_distance   eigenvalue_1 ... eigenvalue_nmat (Ry)'
         do ik = 1, this%nk_path
            write(unit, '(*(ES16.8E3,1X))') this%k_distances(ik), &
               (this%eigenvalues(i, ik), i = 1, nmat)
         end do
         close(unit)
         call g_logger%info('reciprocal%calculate_bsf: path eigenvalues written to ' // &
            trim(bfile), __FILE__, __LINE__)
      end if

      deallocate(gk, sig0, idx_tot)
      if (allocated(idx_up)) deallocate(idx_up, idx_dn)
   end subroutine calculate_bsf

end submodule reciprocal_bsf
