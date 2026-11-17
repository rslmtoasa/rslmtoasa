module strux_tb_driver
!- Array-based RS-LMTO adapter layer for the standalone strux library.
!  This module deliberately does not depend on the legacy lattice class.
!  The intended RS-LMTO integration point is a thin wrapper inside
!  `legacy/lattice.f90` that preserves `this%nn` generation and only replaces
!  the legacy `dbar1` solve with calls to `strux_tb(...)`.
    use strux_lib, only: strux_compute, strux_options, strux_result, &
        STRUX_METHOD_LMTO47, STRUX_LMTO47_IALPHA_MANUAL
    implicit none
    private

    public :: strux_tb

contains

subroutine strux_tb(alat_ang, a, cr, iz, ia, wav_ang, r2_ang2, npold, sbar, sbarvec, sbarvec_direct)
!> Purpose:
!>   Compute legacy-style screened structure constants for one RS-LMTO center.
!> Inputs:
!>   `a`           primitive vectors in units of `alat`
!>   `cr`          cluster coordinates in units of `alat`
!>   `iz`          type index for each cluster atom
!>   `ia`          center atom index in the cluster
!>   `wav_ang`     average Wigner-Seitz radius in Angstrom
!>   `r2_ang2`     legacy cutoff radius squared in Angstrom^2
!>   `npold`       legacy basis size; the first implementation requires 9
!> Outputs:
!>   `sbar`            legacy-ordered screened structure constants
!>   `sbarvec`         vectors used by the Hamiltonian build in Angstrom
!>   `sbarvec_direct`  the same vectors in fractional coordinates
    implicit none

    real(8), intent(in) :: alat_ang, a(3,3), cr(:, :), wav_ang, r2_ang2
    integer, intent(in) :: iz(:), ia, npold
    real(8), allocatable, intent(inout) :: sbar(:,:,:)
    real(8), allocatable, intent(inout) :: sbarvec(:,:), sbarvec_direct(:,:)

    ! STRUX_TB_GENERAL_BASIS:
    ! First implementation is fixed to spd output (9x9).
    ! Generalize to basis-dependent size:
    !   sp   -> 4
    !   spd  -> 9
    !   spdf -> 16
    integer, parameter :: legacy_norb = 9
    integer, parameter :: ncut = 9
    real(8), parameter :: ang2bohr = 1.8897259886d0
    real(8), parameter :: alpha_default(0:3) = [0.3485d0, 0.0530d0, 0.0107d0, 0.00535d0]

    integer :: nat, nspec, nl, lmax, nsolve, nstore, i, j, k, ic
    integer :: orb_map(16), idx_pair
    integer, allocatable :: solve_idx(:), store_idx(:), ips(:), lmxb(:)
    real(8) :: cutoff_solve2, cutoff_store2, dr(3), cell_bohr, maxdist2
    real(8), allocatable :: pos_local(:,:), plat_local(:,:), rmt(:), alpha_in(:,:)
    type(strux_options) :: opts
    type(strux_result) :: result

    nat = size(iz)
    if (ia < 1 .or. ia > nat) then
        write(*,*) 'strux_tb: ia out of range:', ia
        stop
    end if

    ! STRUX_TB_GENERAL_BASIS:
    ! First implementation should stop unless npold==9.
    ! Remove this guard when sp/spdf paths are validated.
    if (npold /= legacy_norb) then
        write(*,*) 'strux_tb: first implementation requires npold=9, got', npold
        stop
    end if

    ! STRUX_TB_GENERAL_BASIS:
    ! Plumb RS-LMTO basis size here. Today we assume npold=9 (spd).
    ! Later derive lmax from npold or a dedicated basis descriptor.
    lmax = 2
    nl = lmax + 1
    nspec = maxval(iz(1:nat))
    cutoff_store2 = r2_ang2
    cutoff_solve2 = dble(ncut) * r2_ang2

    call select_center_indices(cr, alat_ang, ia, cutoff_solve2, solve_idx, nsolve)
    call select_center_indices(cr, alat_ang, ia, cutoff_store2, store_idx, nstore)

    allocate(pos_local(3, nsolve), plat_local(3,3), rmt(nspec), alpha_in(0:nl-1, nspec))
    allocate(ips(nsolve), lmxb(nspec))

    pos_local(:,1) = 0d0
    maxdist2 = 0d0
    do i = 2, nsolve
        dr(:) = (cr(:, solve_idx(i)) - cr(:, ia)) * alat_ang * ang2bohr
        pos_local(:, i) = -dr(:)
        maxdist2 = max(maxdist2, sum(dr*dr))
    end do
    cell_bohr = max(10d0, 4d0*sqrt(maxdist2 + 1d-12) + 10d0*sqrt(cutoff_solve2)*ang2bohr)
    plat_local = 0d0
    do i = 1, 3
        plat_local(i, i) = cell_bohr
    end do

    do i = 1, nsolve
        ips(i) = iz(solve_idx(i))
    end do
    do i = 1, nspec
        lmxb(i) = lmax
        rmt(i) = wav_ang * ang2bohr
        do j = 0, nl-1
            alpha_in(j, i) = alpha_default(min(j, 3))
        end do
    end do

    opts%method = STRUX_METHOD_LMTO47
    opts%auto_alpha = .false.
    opts%want_sdot = .false.
    opts%lmaxw = -1
    opts%pair_cutoff = sqrt(cutoff_solve2) * ang2bohr
    opts%solve_cutoff = sqrt(cutoff_solve2) * ang2bohr
    opts%screening_mode = STRUX_LMTO47_IALPHA_MANUAL

    call strux_compute(opts, nsolve, nspec, nl, 1d0, plat_local, pos_local, ips, lmxb, &
        wav_ang*ang2bohr, rmt, result, alpha_in=alpha_in)

    if (allocated(sbar)) deallocate(sbar)
    if (allocated(sbarvec)) deallocate(sbarvec)
    if (allocated(sbarvec_direct)) deallocate(sbarvec_direct)
    allocate(sbar(legacy_norb, legacy_norb, nstore))
    allocate(sbarvec(3, nstore), sbarvec_direct(3, nstore))
    sbar = 0d0

    ! STRUX_TB_GENERAL_BASIS:
    ! Replace fixed legacy orbital reorder map with basis-dependent map
    ! for sp / spd / spdf.
    call build_orbital_map(legacy_norb, orb_map)

    do i = 1, nstore
        dr(:) = (cr(:, store_idx(i)) - cr(:, ia)) * alat_ang
        sbarvec(:, i) = dr(:)
        call cartesian_to_fractional(dr, sbarvec_direct(:, i), a, alat_ang)

        idx_pair = find_center_pair(result%nttab, result%iax, 1, local_source_index(store_idx(i), solve_idx, nsolve))
        if (idx_pair == 0) then
            write(*,*) 'strux_tb: failed to find pair for center/source', ia, store_idx(i)
            stop
        end if

        ! STRUX_TB_GENERAL_BASIS:
        ! Remove fixed spd assumption from local center solve.
        ! All local pair/block extraction should use nlm, not 9.
        do j = 1, legacy_norb
            do k = 1, legacy_norb
                sbar(j, k, i) = result%s(orb_map(j), orb_map(k), idx_pair)
            end do
        end do
    end do

end subroutine strux_tb

subroutine select_center_indices(cr, alat_ang, ia, cutoff2, idx, nidx)
!> Purpose:
!>   Select the cluster atoms within a spherical cutoff around one center.
!> Notes:
!>   Entry 1 is always the center atom itself so the resulting ordering matches
!>   the local-cluster convention assumed by `strux_tb`.
    implicit none
    real(8), intent(in) :: cr(:, :), alat_ang, cutoff2
    integer, intent(in) :: ia
    integer, allocatable, intent(out) :: idx(:)
    integer, intent(out) :: nidx

    integer :: nat, i, count
    real(8) :: d2

    nat = size(cr, 2)
    allocate(idx(nat))
    idx = 0
    count = 1
    idx(1) = ia
    do i = 1, nat
        if (i == ia) cycle
        d2 = sum(((cr(:, i) - cr(:, ia)) * alat_ang)**2)
        if (d2 < cutoff2 .and. d2 > 1d-4) then
            count = count + 1
            idx(count) = i
        end if
    end do
    nidx = count
end subroutine select_center_indices

integer function local_source_index(global_idx, solve_idx, nsolve)
!> Purpose:
!>   Map one global cluster index back to the local solve-cluster numbering.
    implicit none
    integer, intent(in) :: global_idx, solve_idx(*), nsolve
    integer :: i

    local_source_index = 0
    do i = 1, nsolve
        if (solve_idx(i) == global_idx) then
            local_source_index = i
            return
        end if
    end do
end function local_source_index

integer function find_center_pair(nttab, iax, center_local, source_local)
!> Purpose:
!>   Locate the zero-translation pair block connecting one local source to the center.
!> Notes:
!>   `strux_compute` returns a dense pair table; this helper finds the block
!>   that corresponds to the legacy RS-LMTO vector ordering.
    implicit none
    integer, intent(in) :: nttab, center_local, source_local
    integer, intent(in) :: iax(:,:)
    integer :: i

    find_center_pair = 0
    do i = 1, nttab
        if (iax(2, i) /= center_local) cycle
        if (iax(1, i) /= source_local) cycle
        if (iax(3, i) /= 0 .or. iax(4, i) /= 0 .or. iax(5, i) /= 0) cycle
        find_center_pair = i
        return
    end do
end function find_center_pair

subroutine build_orbital_map(norb, orb_map)
!> Purpose:
!>   Build the permutation from strux ordering to legacy RS-LMTO orbital ordering.
!> Notes:
!>   The current map is hard-coded for the historical `sp`, `spd`, and `spdf`
!>   conventions and is intentionally tagged for later generalization.
    implicit none
    integer, intent(in)  :: norb
    integer, intent(out) :: orb_map(16)

    integer :: i
    integer, parameter :: full_map(16) = [ &
        1, 4, 2, 3, 5, 6, 8, 9, 7, 13, 14, 12, 15, 11, 16, 10 ]

    orb_map = 0
    do i = 1, min(norb, size(full_map))
        orb_map(i) = full_map(i)
    end do
end subroutine build_orbital_map

subroutine cartesian_to_fractional(cart_vec, frac_vec, lattice_vectors, alat)
!> Purpose:
!>   Convert one Cartesian displacement vector to fractional coordinates.
!> Notes:
!>   The wrapper stores both forms because RS-LMTO uses Cartesian vectors for
!>   the Hamiltonian build while external symmetry tooling often expects
!>   fractional ones.
    implicit none
    real(8), intent(in) :: cart_vec(3), lattice_vectors(3,3), alat
    real(8), intent(out) :: frac_vec(3)

    real(8) :: lattice_matrix(3,3), inv(3,3), det
    integer :: i

    do i = 1, 3
        lattice_matrix(:, i) = lattice_vectors(:, i) * alat
    end do
    det = lattice_matrix(1,1)*(lattice_matrix(2,2)*lattice_matrix(3,3) - lattice_matrix(2,3)*lattice_matrix(3,2)) - &
          lattice_matrix(1,2)*(lattice_matrix(2,1)*lattice_matrix(3,3) - lattice_matrix(2,3)*lattice_matrix(3,1)) + &
          lattice_matrix(1,3)*(lattice_matrix(2,1)*lattice_matrix(3,2) - lattice_matrix(2,2)*lattice_matrix(3,1))
    if (abs(det) < 1d-12) then
        write(*,*) 'strux_tb: singular lattice matrix in cartesian_to_fractional'
        stop
    end if
    inv(1, 1) = (lattice_matrix(2, 2) * lattice_matrix(3, 3) - lattice_matrix(2, 3) * lattice_matrix(3, 2)) / det
    inv(1, 2) = (lattice_matrix(1, 3) * lattice_matrix(3, 2) - lattice_matrix(1, 2) * lattice_matrix(3, 3)) / det
    inv(1, 3) = (lattice_matrix(1, 2) * lattice_matrix(2, 3) - lattice_matrix(1, 3) * lattice_matrix(2, 2)) / det
    inv(2, 1) = (lattice_matrix(2, 3) * lattice_matrix(3, 1) - lattice_matrix(2, 1) * lattice_matrix(3, 3)) / det
    inv(2, 2) = (lattice_matrix(1, 1) * lattice_matrix(3, 3) - lattice_matrix(1, 3) * lattice_matrix(3, 1)) / det
    inv(2, 3) = (lattice_matrix(1, 3) * lattice_matrix(2, 1) - lattice_matrix(1, 1) * lattice_matrix(2, 3)) / det
    inv(3, 1) = (lattice_matrix(2, 1) * lattice_matrix(3, 2) - lattice_matrix(2, 2) * lattice_matrix(3, 1)) / det
    inv(3, 2) = (lattice_matrix(1, 2) * lattice_matrix(3, 1) - lattice_matrix(1, 1) * lattice_matrix(3, 2)) / det
    inv(3, 3) = (lattice_matrix(1, 1) * lattice_matrix(2, 2) - lattice_matrix(1, 2) * lattice_matrix(2, 1)) / det

    frac_vec = matmul(inv, cart_vec)
end subroutine cartesian_to_fractional

end module strux_tb_driver
