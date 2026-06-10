 program strux_structb_compare
!- Unified legacy-style file writer built around strux_compute.
!-
!  Usage:
!    strux_structb_compare <exampledir> [--solver lmto47]
!                                      [--basis s|sp|spd|spdf]
!                                      [--pair-cutoff R]
!                                      [--solve-cutoff R]
!                                      [--auto-alpha]
!                                      [--want-sdot]
!                                      [--output-dir DIR]
!-
!  Files written:
!    map        (unformatted)
!    sbar       (unformatted)
!    view.sbar  (formatted)
!    view.sdot  (formatted if --want-sdot)
!    str.out    (formatted)
    use strux_lib
    use example_compare_common
    implicit none

    character(len=256) :: exdir, outdir, solver, basis_name
    integer :: nargs, i
    integer :: nbas, nspec, nl, lmax
    integer :: ib, jb, n1, n2, n3, ic, row, col, count, u_map, u_sbar, u_view, u_viewdot, u_str
    integer :: orb_map(16)
    integer :: ntab_arr(101)
    integer :: nblocks, max_neighbors, norb_out
    integer :: no_arr(100), izp_arr(100)
    integer, allocatable :: ips(:), lmxb(:), block_idx(:), nbrs(:,:)
    integer, allocatable :: nbr_count(:)
    real(8) :: alat_ang, alat_bohr, wav_bohr, r2, avw, hcr_ratio, kap1
    real(8) :: plat(3,3), crd(3,100), vec(3)
    real(8) :: ct_arr(10)
    character(len=8) :: label_arr(10)
    logical :: have_basis, have_pair_cutoff, have_solve_cutoff, auto_alpha, want_sdot
    integer :: screening_mode
    real(8) :: pair_cutoff_ang, solve_cutoff_ang
    real(8), allocatable :: pos(:,:), rmt(:), alpha_in(:,:), hcr(:,:)
    real(8), allocatable :: alpha_screen(:,:), adot_screen(:,:), tral(:,:,:), trad(:,:,:)
    real(8), allocatable :: alpha_l_report(:,:), adot_l_report(:,:)
    type(strux_options) :: opts
    type(strux_result) :: result

    nargs = command_argument_count()
    call parse_compare_cli(nargs, exdir, solver, have_basis, basis_name, have_pair_cutoff, pair_cutoff_ang, &
        have_solve_cutoff, solve_cutoff_ang, auto_alpha, want_sdot, outdir, screening_mode)

    nbas = 1
    plat = 0d0
    crd = 0d0
    no_arr = 1
    izp_arr = 1
    call parse_lattice_nml(trim(exdir)//'/lattice.nml', nbas, plat, crd, no_arr)

    alat_ang = 0d0
    wav_bohr = 0d0
    r2 = 9d0
    nspec = 1
    ct_arr = 3d0
    label_arr = 'XX'
    call parse_input_nml(trim(exdir)//'/input.nml', alat_ang, wav_bohr, r2, nspec, ct_arr, label_arr, lmax, hcr_ratio, kap1)

    if (have_basis) call basis_name_to_lmax(basis_name, lmax)
    nl = lmax + 1
    if (nl > 4) then
        write(*,'(a)') 'strux_structb_compare: basis larger than spdf not supported in compare output'
        stop 1
    end if
    norb_out = nl*nl

    if (.not. have_pair_cutoff) pair_cutoff_ang = sqrt(r2)
    if (.not. have_solve_cutoff) solve_cutoff_ang = 3d0 * pair_cutoff_ang

    alat_bohr = alat_ang * ang2bohr
    wav_bohr = wav_bohr * ang2bohr
    avw = wav_bohr

    allocate(lmxb(nspec), rmt(nspec), ips(nbas), pos(3,nbas))
    do i = 1, nspec
        lmxb(i) = lmax
        rmt(i) = avw
    end do
    do ib = 1, nbas
        ips(ib) = no_arr(ib)
        pos(:, ib) = crd(:, ib)
    end do

    opts%method = solver_name_to_method(solver)
    opts%auto_alpha = auto_alpha
    opts%want_sdot = want_sdot
    opts%lmaxw = -1
    opts%pair_cutoff = pair_cutoff_ang * ang2bohr
    opts%solve_cutoff = solve_cutoff_ang * ang2bohr
    opts%screening_mode = screening_mode

    select case (opts%screening_mode)
    case (STRUX_LMTO47_IALPHA_SIGMA)
        allocate(hcr(nl, nspec))
        call fill_uniform_hcr(nl, nspec, hcr_ratio, rmt, hcr)
        call strux_compute(opts, nbas, nspec, nl, alat_bohr, plat, pos, ips, lmxb, avw, rmt, result, hcr=hcr)
    case (STRUX_LMTO47_IALPHA_FITD)
        call strux_compute(opts, nbas, nspec, nl, alat_bohr, plat, pos, ips, lmxb, avw, rmt, result)
    case default
        allocate(alpha_in(0:nl-1, nspec))
        call fill_default_alpha(nl, nspec, alpha_in)
        if (opts%want_sdot) then
            write(*,'(a)') 'compare: lmto47 manual-alpha Sdot is undefined unless alpha and hard-core radii are mutually consistent.'
            write(*,'(a)') 'compare: use --screening-mode sigma or --screening-mode fitd for LMTO47 derivative comparisons.'
            stop 1
        else
            call strux_compute(opts, nbas, nspec, nl, alat_bohr, plat, pos, ips, lmxb, avw, rmt, result, alpha_in=alpha_in)
        end if
    end select

    call build_orbital_map(norb_out, orb_map)
    call build_pair_blocks(result%nttab, alat_ang, plat, pos, result%iax, pair_cutoff_ang*pair_cutoff_ang, nblocks, block_idx, ntab_arr)
    call build_pair_neighbor_map(nblocks, block_idx, result%iax, nbr_count, nbrs, max_neighbors)

    open(newunit=u_map,  file=trim(outdir)//'/map',       form='unformatted', status='replace')
    open(newunit=u_sbar, file=trim(outdir)//'/sbar',      form='unformatted', status='replace')
    open(newunit=u_view, file=trim(outdir)//'/view.sbar', status='replace')
    if (opts%want_sdot) open(newunit=u_viewdot, file=trim(outdir)//'/view.sdot', status='replace')
    open(newunit=u_str,  file=trim(outdir)//'/str.out',   status='replace')

    write(u_str, *) 'irec', nbas, (i, i=1, nbas)
    write(u_str, *) 'irec type', (ips(i), i=1, nbas)
    write(u_str, *) 'ndi=', nbas
    write(u_str, '(i5)') nbas
    write(u_str, '(" LATTICE COORDINATES")')
    do i = 1, nbas
        write(u_str, '(i5,3f8.4)') i, (pos(ic, i)*alat_ang, ic=1,3)
    end do

    call write_text_map(u_str, ips, nbr_count, nbrs, nbas, max_neighbors)
    write(u_str, '(3i5)') nbas, nblocks, 0
    write(u_str, '(a,a)') 'solver=', trim(solver)
    write(u_str, '(a,l1)') 'auto_alpha=', opts%auto_alpha
    write(u_str, '(a,l1)') 'want_sdot=', opts%want_sdot
    write(u_str, '(a,f12.6)') 'pair_cutoff=', pair_cutoff_ang
    write(u_str, '(a,f12.6)') 'solve_cutoff=', solve_cutoff_ang
    write(u_str, '(a,i0)') 'screening_mode=', opts%screening_mode
    call write_screening_report(u_str)

    do i = 1, nblocks
        ib = result%iax(1, block_idx(i))
        jb = result%iax(2, block_idx(i))
        n1 = result%iax(3, block_idx(i))
        n2 = result%iax(4, block_idx(i))
        n3 = result%iax(5, block_idx(i))

        do ic = 1, 3
            vec(ic) = alat_ang * ( pos(ic,jb) - pos(ic,ib) + &
                dble(n1)*plat(ic,1) + dble(n2)*plat(ic,2) + dble(n3)*plat(ic,3) )
        end do

        write(u_view, '(" VECTOR=",3f12.6,"   ICLUS=",i5)') vec(1), vec(2), vec(3), jb
        if (opts%want_sdot) write(u_viewdot, '(" VECTOR=",3f12.6,"   ICLUS=",i5)') vec(1), vec(2), vec(3), jb
        do row = 1, norb_out
            do col = 1, norb_out
                write(u_sbar) result%s(orb_map(row), orb_map(col), block_idx(i))
            end do
            write(u_view, '(*(f10.4))') (result%s(orb_map(row), orb_map(col), block_idx(i)), col=1, norb_out)
            if (opts%want_sdot) then
                write(u_viewdot, '(*(f10.4))') (result%sdot(orb_map(row), orb_map(col), block_idx(i)), col=1, norb_out)
            end if
        end do
    end do

    do i = 1, nbas
        count = nbr_count(i) + 1
        write(u_map) count, (nbrs(i, col), col=1, nbr_count(i))
    end do

    close(u_map)
    close(u_sbar)
    close(u_view)
    if (opts%want_sdot) close(u_viewdot)
    close(u_str)

contains

subroutine write_screening_report(unit_no)
    implicit none
    integer, intent(in) :: unit_no
    integer :: is, l

    call compute_screening_report(alpha_l_report, adot_l_report)
    if (.not. allocated(alpha_l_report)) return

    write(unit_no, '(a)') 'alpha_l_used:'
    write(*, '(a)') 'alpha_l_used:'
    do is = 1, size(alpha_l_report, 2)
        write(unit_no, '(a,i0,a,*(1x,f12.6))') '  species ', is, ':', (alpha_l_report(l, is), l=0, size(alpha_l_report, 1)-1)
        write(*, '(a,i0,a,*(1x,f12.6))') '  species ', is, ':', (alpha_l_report(l, is), l=0, size(alpha_l_report, 1)-1)
    end do

    if (allocated(adot_l_report)) then
        write(unit_no, '(a)') 'adot_l_used:'
        write(*, '(a)') 'adot_l_used:'
        do is = 1, size(adot_l_report, 2)
            write(unit_no, '(a,i0,a,*(1x,f12.6))') '  species ', is, ':', (adot_l_report(l, is), l=0, size(adot_l_report, 1)-1)
            write(*, '(a,i0,a,*(1x,f12.6))') '  species ', is, ':', (adot_l_report(l, is), l=0, size(adot_l_report, 1)-1)
        end do
    end if
end subroutine write_screening_report

subroutine compute_screening_report(alpha_l_used, adot_l_used)
    implicit none
    real(8), allocatable, intent(out) :: alpha_l_used(:,:), adot_l_used(:,:)
    integer :: is, l, lm1, ib_site
    real(8), allocatable :: alpha_l_auto(:,:)

    if (allocated(alpha_l_used)) deallocate(alpha_l_used)
    if (allocated(adot_l_used)) deallocate(adot_l_used)

    if (opts%screening_mode == STRUX_LMTO47_IALPHA_FITD) then
        call strux_lmto47_autoalpha_screening(nbas, nspec, nl, avw, ips, lmxb, rmt, &
            alpha_l_out=alpha_l_auto, alpha_out=alpha_screen, adot=adot_screen, tral=tral, trad=trad, &
            screening_mode=opts%screening_mode)
        call move_alloc(alpha_l_auto, alpha_l_used)
    else if (opts%screening_mode == STRUX_LMTO47_IALPHA_SIGMA) then
        call strux_lmto47_autoalpha_screening(nbas, nspec, nl, avw, ips, lmxb, rmt, &
            hcr=hcr, alpha_l_out=alpha_l_auto, alpha_out=alpha_screen, adot=adot_screen, tral=tral, trad=trad, &
            screening_mode=opts%screening_mode)
        call move_alloc(alpha_l_auto, alpha_l_used)
    else if (auto_alpha) then
        call strux_lmto47_autoalpha_screening(nbas, nspec, nl, avw, ips, lmxb, rmt, &
            hcr=hcr, alpha_l_out=alpha_l_auto, alpha_out=alpha_screen, adot=adot_screen, tral=tral, trad=trad)
        call move_alloc(alpha_l_auto, alpha_l_used)
    else
        allocate(alpha_l_used(0:nl-1, nspec))
        alpha_l_used = result%alpha_l
        call strux_lmto47_screening(nbas, nspec, nl, avw, ips, lmxb, rmt, result%alpha_l, &
            alpha_screen, adot_screen, tral, trad)
    end if

    allocate(adot_l_used(0:nl-1, nspec))
    adot_l_used = 0d0
    do is = 1, nspec
        ib_site = 0
        do l = 1, nbas
            if (ips(l) == is) then
                ib_site = l
                exit
            end if
        end do
        if (ib_site == 0) cycle
        do l = 0, nl-1
            lm1 = l*l + 1
            adot_l_used(l, is) = adot_screen(lm1, ib_site)
        end do
    end do

    if (allocated(alpha_screen)) deallocate(alpha_screen)
    if (allocated(adot_screen)) deallocate(adot_screen)
    if (allocated(tral)) deallocate(tral)
    if (allocated(trad)) deallocate(trad)
end subroutine compute_screening_report

subroutine build_pair_blocks(nttab, alat_ang, plat, pos, iax, pair_cutoff2, nblocks, block_idx, ntab_arr)
    implicit none
    integer, intent(in) :: nttab
    real(8), intent(in) :: alat_ang, plat(3,3), pos(3,*), pair_cutoff2
    integer, intent(in) :: iax(:,:)
    integer, intent(out) :: nblocks
    integer, allocatable, intent(out) :: block_idx(:)
    integer, intent(out) :: ntab_arr(:)
    integer :: i, ib, count
    real(8) :: vec(3), d2

    allocate(block_idx(nttab))
    ntab_arr = 0
    count = 0
    do ib = 1, size(ntab_arr)-1
        ntab_arr(ib) = count
        do i = 1, nttab
            if (iax(2, i) /= ib) cycle
            call pair_vector(alat_ang, plat, pos, iax(:, i), vec)
            d2 = sum(vec*vec)
            if (d2 <= pair_cutoff2 + 1d-8) then
                count = count + 1
                block_idx(count) = i
            end if
        end do
    end do
    ntab_arr(size(ntab_arr)) = count
    nblocks = count
end subroutine build_pair_blocks

subroutine build_pair_neighbor_map(nblocks, block_idx, iax, nbr_count, nbrs, max_neighbors)
    implicit none
    integer, intent(in) :: nblocks, block_idx(:)
    integer, intent(in) :: iax(:,:)
    integer, allocatable, intent(out) :: nbr_count(:), nbrs(:,:)
    integer, intent(out) :: max_neighbors
    integer :: i, j, nat
    logical, allocatable :: seen(:,:)

    nat = maxval(iax(2, block_idx(:nblocks)))
    allocate(nbr_count(nat), seen(nat, nat))
    nbr_count = 0
    seen = .false.

    do i = 1, nblocks
        j = iax(1, block_idx(i))
        if (j <= 0 .or. j > nat) cycle
        if (.not. seen(iax(2, block_idx(i)), j)) then
            seen(iax(2, block_idx(i)), j) = .true.
            nbr_count(iax(2, block_idx(i))) = nbr_count(iax(2, block_idx(i))) + 1
        end if
    end do

    max_neighbors = max(1, maxval(nbr_count))
    allocate(nbrs(nat, max_neighbors))
    nbrs = 0
    nbr_count = 0
    seen = .false.

    do i = 1, nblocks
        j = iax(1, block_idx(i))
        if (j <= 0 .or. j > nat) cycle
        if (.not. seen(iax(2, block_idx(i)), j)) then
            seen(iax(2, block_idx(i)), j) = .true.
            nbr_count(iax(2, block_idx(i))) = nbr_count(iax(2, block_idx(i))) + 1
            nbrs(iax(2, block_idx(i)), nbr_count(iax(2, block_idx(i)))) = j
        end if
    end do

    deallocate(seen)
end subroutine build_pair_neighbor_map

subroutine write_text_map(unit_no, ips, nbr_count, nbrs, nat, max_neighbors)
    implicit none
    integer, intent(in) :: unit_no, nat, max_neighbors
    integer, intent(in) :: ips(nat), nbr_count(nat), nbrs(nat, max_neighbors)
    integer :: i, j, idm

    write(unit_no, '(" NEAREST NEIGHBOUR MAP"/,"      ATOM   TYPE  CONNECTIVITY",5x,"NEIGHBOURS")')
    do i = 1, nat
        idm = min(nbr_count(i), 16)
        write(unit_no, '(1x,3i4,2x,20i5)') i, ips(i), nbr_count(i) + 1, (nbrs(i, j), j=1, idm)
        if (nbr_count(i) > 16) then
            write(unit_no, '(29x,16i5)') (nbrs(i, j), j=17, nbr_count(i))
        end if
    end do
end subroutine write_text_map

subroutine pair_vector(alat_ang, plat, pos, iax_col, vec)
    implicit none
    real(8), intent(in) :: alat_ang, plat(3,3), pos(3,*)
    integer, intent(in) :: iax_col(:)
    real(8), intent(out) :: vec(3)
    integer :: ic

    do ic = 1, 3
        vec(ic) = alat_ang * (pos(ic, iax_col(2)) - pos(ic, iax_col(1)) + &
            dble(iax_col(3))*plat(ic,1) + dble(iax_col(4))*plat(ic,2) + dble(iax_col(5))*plat(ic,3))
    end do
end subroutine pair_vector

end program strux_structb_compare
