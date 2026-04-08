module strux_lib
!- Public interface for the standalone LMTO47 structure-constants library.
!  Retained supported scope:
!    strux_compute
!    strux_lmto47
!    strux_lmto47_autoalpha
!    strux_lmto47_screening
!    strux_lmto47_autoalpha_screening
    use strux_lmto47_mode_solver, only: strscr_lmto47, &
        strux_lmto47_screening_impl => strux_lmto47_screening, &
        strux_lmto47_autoalpha_screening_impl => strux_lmto47_autoalpha_screening
    use strux_support, only: build_common_tables, allocate_sdot_array, &
        validate_common_inputs, strux_fail
    implicit none
    private

    integer, parameter, public :: STRUX_METHOD_LMTO47 = 4
    integer, parameter, public :: STRUX_LMTO47_IALPHA_SIGMA = 0
    integer, parameter, public :: STRUX_LMTO47_IALPHA_FITD = 1
    integer, parameter, public :: STRUX_LMTO47_IALPHA_MANUAL = 2

    type, public :: strux_options
        integer :: method = STRUX_METHOD_LMTO47
        logical :: auto_alpha = .false.
        logical :: want_sdot = .false.
        integer :: lmaxw = -1
        real(8) :: pair_cutoff = -1d0
        real(8) :: solve_cutoff = -1d0
        integer :: screening_mode = STRUX_LMTO47_IALPHA_MANUAL
    end type strux_options

    type, public :: strux_result
        integer :: nttab = 0
        integer, allocatable :: iax(:,:)
        real(8), allocatable :: alpha(:,:)
        real(8), allocatable :: alpha_l(:,:)
        real(8), allocatable :: s(:,:,:)
        real(8), allocatable :: sdot(:,:,:)
    end type strux_result

    public :: strux_compute
    public :: strux_lmto47
    public :: strux_lmto47_autoalpha
    public :: strux_lmto47_screening
    public :: strux_lmto47_autoalpha_screening

contains

subroutine strux_compute(opts, nbas, nspec, nl, alat, plat, pos, ips, lmxb, &
    avw, rmt, result, alpha_in, hcr)
    implicit none
    type(strux_options), intent(in) :: opts
    integer, intent(in) :: nbas, nspec, nl
    double precision, intent(in) :: alat, plat(3,3), pos(3,nbas), avw, rmt(nspec)
    integer, intent(in) :: ips(nbas), lmxb(nspec)
    type(strux_result), intent(inout) :: result
    double precision, intent(in), optional :: alpha_in(0:nl-1, nspec)
    double precision, intent(in), optional :: hcr(nl, nspec)

    real(8) :: pair_cutoff, solve_cutoff
    integer, allocatable :: iax_full(:,:)
    real(8), allocatable :: alpha_full(:,:), s_full(:,:,:), sdot_full(:,:,:)
    real(8), allocatable :: alpha_l_used(:,:)
    integer :: ib

    call validate_common_inputs('strux_compute', nbas, nspec, nl, avw, 1d0)
    do ib = 1, nbas
        if (ips(ib) < 1 .or. ips(ib) > nspec) then
            write(*,'(a,i0,a,i0,a,i0)') 'strux_compute: invalid ips(', ib, ')=', ips(ib), ' for nspec=', nspec
            stop 1
        end if
    end do

    if (opts%method /= STRUX_METHOD_LMTO47) then
        call strux_fail('strux_compute: only STRUX_METHOD_LMTO47 is supported')
    end if

    pair_cutoff = opts%pair_cutoff
    solve_cutoff = opts%solve_cutoff
    if (pair_cutoff <= 0d0 .and. solve_cutoff <= 0d0) call strux_fail('strux_compute: one of pair_cutoff or solve_cutoff must be > 0')
    if (pair_cutoff <= 0d0) pair_cutoff = solve_cutoff
    if (solve_cutoff <= 0d0) solve_cutoff = pair_cutoff
    if (solve_cutoff < pair_cutoff) call strux_fail('strux_compute: solve_cutoff must be >= pair_cutoff')

    if (allocated(result%iax)) deallocate(result%iax)
    if (allocated(result%alpha)) deallocate(result%alpha)
    if (allocated(result%alpha_l)) deallocate(result%alpha_l)
    if (allocated(result%s)) deallocate(result%s)
    if (allocated(result%sdot)) deallocate(result%sdot)
    result%nttab = 0

    if (opts%screening_mode == STRUX_LMTO47_IALPHA_FITD) then
        call strux_lmto47_autoalpha(nbas, nspec, nl, alat, plat, pos, ips, lmxb, &
            rmt, hcr, avw, solve_cutoff, opts%lmaxw, opts%want_sdot, &
            result%nttab, iax_full, alpha_l_used, alpha_full, s_full, sdot_full, screening_mode=opts%screening_mode)
    else if (opts%screening_mode == STRUX_LMTO47_IALPHA_SIGMA) then
        if (.not. present(hcr)) call strux_fail('strux_compute: hcr is required for LMTO47 screening_mode=sigma')
        call strux_lmto47_autoalpha(nbas, nspec, nl, alat, plat, pos, ips, lmxb, &
            rmt, hcr, avw, solve_cutoff, opts%lmaxw, opts%want_sdot, &
            result%nttab, iax_full, alpha_l_used, alpha_full, s_full, sdot_full, screening_mode=opts%screening_mode)
    else if (opts%auto_alpha) then
        if (.not. present(hcr)) call strux_fail('strux_compute: hcr is required when auto_alpha=.true.')
        call strux_lmto47_autoalpha(nbas, nspec, nl, alat, plat, pos, ips, lmxb, &
            rmt, hcr, avw, solve_cutoff, opts%lmaxw, opts%want_sdot, &
            result%nttab, iax_full, alpha_l_used, alpha_full, s_full, sdot_full, screening_mode=opts%screening_mode)
    else
        if (.not. present(alpha_in)) call strux_fail('strux_compute: alpha_in is required unless automatic LMTO47 screening is selected')
        allocate(alpha_l_used(0:nl-1, nspec))
        alpha_l_used = alpha_in
        call strux_lmto47(nbas, nspec, nl, alat, plat, pos, ips, lmxb, &
            rmt, alpha_in, avw, solve_cutoff, opts%lmaxw, opts%want_sdot, &
            result%nttab, iax_full, alpha_full, s_full, sdot_full, hcr=hcr, screening_mode=opts%screening_mode)
    end if

    call move_alloc(alpha_l_used, result%alpha_l)
    call move_alloc(alpha_full, result%alpha)

    call filter_pair_outputs(alat, plat, pos, pair_cutoff, opts%want_sdot, iax_full, s_full, sdot_full, &
        result%nttab, result%iax, result%s, result%sdot)

    if (allocated(iax_full)) deallocate(iax_full)
    if (allocated(s_full)) deallocate(s_full)
    if (allocated(sdot_full)) deallocate(sdot_full)
end subroutine strux_compute

subroutine filter_pair_outputs(alat, plat, pos, pair_cutoff, want_sdot, iax_in, s_in, sdot_in, &
    nttab_out, iax_out, s_out, sdot_out)
    implicit none
    double precision, intent(in) :: alat, plat(3,3), pos(3,*), pair_cutoff
    logical, intent(in) :: want_sdot
    integer, intent(in) :: iax_in(:,:)
    double precision, intent(in) :: s_in(:,:,:)
    double precision, intent(in) :: sdot_in(:,:,:)
    integer, intent(out) :: nttab_out
    integer, allocatable, intent(out) :: iax_out(:,:)
    double precision, allocatable, intent(out) :: s_out(:,:,:)
    double precision, allocatable, intent(out) :: sdot_out(:,:,:)

    integer :: nttab_in, niax, nl2, i, keep_count
    logical, allocatable :: keep(:)
    double precision :: vec(3), d2, cutoff2

    niax = size(iax_in, 1)
    nttab_in = size(iax_in, 2)
    nl2 = size(s_in, 1)
    cutoff2 = pair_cutoff * pair_cutoff

    allocate(keep(nttab_in))
    keep = .false.
    keep_count = 0
    do i = 1, nttab_in
        call pair_vector(alat, plat, pos, iax_in(:,i), vec)
        d2 = sum(vec*vec)
        if (d2 <= cutoff2 + 1d-10) then
            keep(i) = .true.
            keep_count = keep_count + 1
        end if
    end do

    nttab_out = keep_count
    allocate(iax_out(niax, nttab_out))
    allocate(s_out(nl2, nl2, nttab_out))
    if (want_sdot) then
        allocate(sdot_out(nl2, nl2, nttab_out))
    else
        allocate(sdot_out(1,1,1))
    end if

    keep_count = 0
    do i = 1, nttab_in
        if (.not. keep(i)) cycle
        keep_count = keep_count + 1
        iax_out(:, keep_count) = iax_in(:, i)
        s_out(:,:,keep_count) = s_in(:,:,i)
        if (want_sdot) sdot_out(:,:,keep_count) = sdot_in(:,:,i)
    end do

    deallocate(keep)
end subroutine filter_pair_outputs

subroutine pair_vector(alat, plat, pos, iax_col, vec)
    implicit none
    double precision, intent(in) :: alat, plat(3,3), pos(3,*)
    integer, intent(in) :: iax_col(:)
    double precision, intent(out) :: vec(3)
    integer :: ic

    do ic = 1, 3
        vec(ic) = alat * (pos(ic, iax_col(2)) - pos(ic, iax_col(1)) + &
            dble(iax_col(3))*plat(ic,1) + dble(iax_col(4))*plat(ic,2) + dble(iax_col(5))*plat(ic,3))
    end do
end subroutine pair_vector

subroutine strux_lmto47(nbas, nspec, nl, alat, plat, pos, ips, lmxb, &
    rmt, alpha_in, avw, rmax, lmaxw, ldot, &
    nttab, iax, alpha_out, s_out, sdot_out, hcr, screening_mode)
    implicit none
    integer,  intent(in)  :: nbas, nspec, nl, lmaxw
    logical,  intent(in)  :: ldot
    double precision, intent(in) :: alat, plat(3,3), avw, rmax
    double precision, intent(in) :: pos(3,nbas), rmt(nspec)
    double precision, intent(in) :: alpha_in(0:nl-1, nspec)
    integer,  intent(in)  :: ips(nbas), lmxb(nspec)
    integer,  intent(out) :: nttab
    integer,  intent(inout), allocatable :: iax(:,:)
    double precision, intent(inout), allocatable :: alpha_out(:,:)
    double precision, intent(inout), allocatable :: s_out(:,:,:)
    double precision, intent(inout), allocatable :: sdot_out(:,:,:)
    double precision, intent(in), optional :: hcr(nl, nspec)
    integer, intent(in), optional :: screening_mode

    integer :: nl2, ntab_arr(nbas+1)
    double precision :: cy(100), cg(1200), kap1(1)
    integer :: indxcg(350), jcg(1200)
    double precision, allocatable :: adot(:,:), tral(:,:,:), trad(:,:,:)
    double precision, allocatable :: bas_loc(:,:), rwats(:), plat_loc(:,:)
    double precision, allocatable :: s_work(:,:,:,:,:), sdot_work(:,:,:,:,:)

    if (ldot .and. .not. present(hcr)) then
        call strux_fail('strux_lmto47: manual-alpha Sdot requires hcr so LMTO47 derivative screening data can be built')
    end if

    nl2 = nl*nl
    call validate_common_inputs('strux_lmto47', nbas, nspec, nl, avw, rmax)
    call build_common_tables(nbas, nspec, alat, plat, pos, rmax, ips, lmxb, &
        nttab, ntab_arr, iax, cy, cg, indxcg, jcg)

    call strux_lmto47_screening_impl(nbas, nspec, nl, avw, ips, lmxb, rmt, alpha_in, &
        alpha_out, adot, tral, trad, hcr, screening_mode)

    allocate(s_out(nl2, nl2, nttab))
    s_out = 0d0
    call allocate_sdot_array(ldot, nl2, nttab, sdot_out)

    allocate(s_work(nl2, nl2, 1, 1, nttab))
    s_work = 0d0
    if (ldot) then
        allocate(sdot_work(nl2, nl2, 1, 1, nttab))
    else
        allocate(sdot_work(1, 1, 1, 1, 1))
    end if
    sdot_work = 0d0

    allocate(rwats(nbas), bas_loc(3, nbas), plat_loc(3,3))
    rwats = 0d0
    if (lmaxw >= 0) rwats = rmt(ips)
    bas_loc = pos
    plat_loc = plat
    kap1(1) = 0d0

    call strscr_lmto47(nbas, 0, 0, alat/avw, plat_loc, bas_loc, rwats, &
        nl, kap1, 1, lmaxw, cy, cg, indxcg, jcg, ldot, &
        ntab_arr, iax, alpha_out, adot, tral, trad, s_work, sdot_work)

    s_out = s_work(:,:,1,1,:)
    if (ldot) sdot_out = sdot_work(:,:,1,1,:)

    deallocate(adot, tral, trad, bas_loc, rwats, plat_loc, s_work, sdot_work)
end subroutine strux_lmto47

subroutine strux_lmto47_autoalpha(nbas, nspec, nl, alat, plat, pos, ips, lmxb, &
    rmt, hcr, avw, rmax, lmaxw, ldot, &
    nttab, iax, alpha_l_out, alpha_out, s_out, sdot_out, screening_mode)
    implicit none
    integer,  intent(in)  :: nbas, nspec, nl, lmaxw
    logical,  intent(in)  :: ldot
    double precision, intent(in) :: alat, plat(3,3), avw, rmax
    double precision, intent(in) :: pos(3,nbas), rmt(nspec)
    double precision, intent(in), optional :: hcr(nl, nspec)
    integer,  intent(in)  :: ips(nbas), lmxb(nspec)
    integer,  intent(out) :: nttab
    integer,  intent(inout), allocatable :: iax(:,:)
    double precision, intent(inout), allocatable :: alpha_l_out(:,:)
    double precision, intent(inout), allocatable :: alpha_out(:,:)
    double precision, intent(inout), allocatable :: s_out(:,:,:)
    double precision, intent(inout), allocatable :: sdot_out(:,:,:)
    integer, intent(in), optional :: screening_mode

    integer :: nl2, ntab_arr(nbas+1)
    double precision :: cy(100), cg(1200), kap1(1)
    integer :: indxcg(350), jcg(1200)
    double precision, allocatable :: adot(:,:), tral(:,:,:), trad(:,:,:)
    double precision, allocatable :: bas_loc(:,:), rwats(:), plat_loc(:,:)
    double precision, allocatable :: s_work(:,:,:,:,:), sdot_work(:,:,:,:,:)

    nl2 = nl*nl
    call validate_common_inputs('strux_lmto47_autoalpha', nbas, nspec, nl, avw, rmax)
    call build_common_tables(nbas, nspec, alat, plat, pos, rmax, ips, lmxb, &
        nttab, ntab_arr, iax, cy, cg, indxcg, jcg)
    call strux_lmto47_autoalpha_screening_impl(nbas, nspec, nl, avw, ips, lmxb, rmt, hcr, &
        alpha_l_out, alpha_out, adot, tral, trad, screening_mode)

    allocate(s_out(nl2, nl2, nttab))
    s_out = 0d0
    call allocate_sdot_array(ldot, nl2, nttab, sdot_out)

    allocate(s_work(nl2, nl2, 1, 1, nttab))
    s_work = 0d0
    if (ldot) then
        allocate(sdot_work(nl2, nl2, 1, 1, nttab))
    else
        allocate(sdot_work(1, 1, 1, 1, 1))
    end if
    sdot_work = 0d0

    allocate(rwats(nbas), bas_loc(3, nbas), plat_loc(3,3))
    rwats = 0d0
    if (lmaxw >= 0) rwats = rmt(ips)
    bas_loc = pos
    plat_loc = plat
    kap1(1) = 0d0

    call strscr_lmto47(nbas, 0, 0, alat/avw, plat_loc, bas_loc, rwats, &
        nl, kap1, 1, lmaxw, cy, cg, indxcg, jcg, ldot, &
        ntab_arr, iax, alpha_out, adot, tral, trad, s_work, sdot_work)

    s_out = s_work(:,:,1,1,:)
    if (ldot) sdot_out = sdot_work(:,:,1,1,:)

    deallocate(adot, tral, trad, bas_loc, rwats, plat_loc, s_work, sdot_work)
end subroutine strux_lmto47_autoalpha

subroutine strux_lmto47_screening(nbas, nspec, nl, avw, ips, lmxb, rmt, alpha_in, &
    alpha_out, adot, tral, trad, hcr, screening_mode)
    implicit none
    integer, intent(in) :: nbas, nspec, nl
    double precision, intent(in) :: avw, rmt(nspec), alpha_in(0:nl-1, nspec)
    integer, intent(in) :: ips(nbas), lmxb(nspec)
    double precision, intent(in), optional :: hcr(nl, nspec)
    integer, intent(in), optional :: screening_mode
    double precision, intent(inout), allocatable :: alpha_out(:,:), adot(:,:)
    double precision, intent(inout), allocatable :: tral(:,:,:), trad(:,:,:)

    call strux_lmto47_screening_impl(nbas, nspec, nl, avw, ips, lmxb, rmt, &
        alpha_in, alpha_out, adot, tral, trad, hcr, screening_mode)
end subroutine strux_lmto47_screening

subroutine strux_lmto47_autoalpha_screening(nbas, nspec, nl, avw, ips, lmxb, rmt, hcr, &
    alpha_l_out, alpha_out, adot, tral, trad, screening_mode)
    implicit none
    integer, intent(in) :: nbas, nspec, nl
    double precision, intent(in) :: avw, rmt(nspec)
    double precision, intent(in), optional :: hcr(nl, nspec)
    integer, intent(in) :: ips(nbas), lmxb(nspec)
    integer, intent(in), optional :: screening_mode
    double precision, intent(inout), allocatable :: alpha_l_out(:,:)
    double precision, intent(inout), allocatable :: alpha_out(:,:), adot(:,:)
    double precision, intent(inout), allocatable :: tral(:,:,:), trad(:,:,:)

    call strux_lmto47_autoalpha_screening_impl(nbas, nspec, nl, avw, ips, lmxb, rmt, hcr, &
        alpha_l_out, alpha_out, adot, tral, trad, screening_mode)
end subroutine strux_lmto47_autoalpha_screening

end module strux_lib
