module strux_lmto47_mode_solver
!========================================================================
! LMTO47-style 2nd-generation screened structure-constant solver.
!
! This module mirrors the LMTO47 `salph0` orchestration for the supported
! one-kappa scope:
!   screening generation -> transfer matrices -> optional Watson augmentation
!   -> cluster solve -> LMTO47 `scals` transform -> scatter to pair table.
!========================================================================
    use strux_utils, only: dpzero, iinit, isum, ll
    use lmto_strux_bessel, only: bessl2
    use strux_harmonics, only: sylm
    use strux_geometry_math, only: drr2x
    use strux_pair_expansion, only: mstrx2, mstrx3
    use strux_solver_assembly, only: alwats, addtos, symstr
    use strux_support, only: strux_fail
    implicit none
    private

    public :: strux_lmto47_screening
    public :: strux_lmto47_autoalpha_screening
    public :: strscr_lmto47

contains

subroutine maked_lmto47(nl, l, kap2, dd, dddot)
!> Purpose:
!>   Port of TBSTR `maked` for fitted logarithmic derivatives.
!> Notes:
!>   This is the data-driven linear fit used by LMTO47 `ialpha=1`.
    implicit none
    integer, intent(in) :: nl, l
    double precision, intent(in) :: kap2
    double precision, intent(out) :: dd, dddot

    integer, parameter :: nlmax = 10
    double precision, parameter :: dd0(0:nlmax-1,1:nlmax) = reshape([ &
        0.75d0, 1.01d0, 2.01d0, 3.01d0, 4.01d0, 5.01d0, 6.01d0, 7.01d0, 8.01d0, 9.01d0, &
        1.35d0, 1.55d0, 2.01d0, 3.01d0, 4.01d0, 5.01d0, 6.01d0, 7.01d0, 8.01d0, 9.01d0, &
        2.3d0, 2.40d0, 2.60d0, 3.01d0, 4.01d0, 5.01d0, 6.01d0, 7.01d0, 8.01d0, 9.01d0, &
        2.3d0, 2.40d0, 2.60d0, 3.01d0, 4.01d0, 5.01d0, 6.01d0, 7.01d0, 8.01d0, 9.01d0, &
        2.3d0, 2.40d0, 2.60d0, 3.01d0, 4.01d0, 5.01d0, 6.01d0, 7.01d0, 8.01d0, 9.01d0, &
        2.3d0, 2.40d0, 2.60d0, 3.01d0, 4.01d0, 5.01d0, 6.01d0, 7.01d0, 8.01d0, 9.01d0, &
        2.3d0, 2.40d0, 2.60d0, 3.01d0, 4.01d0, 5.01d0, 6.01d0, 7.01d0, 8.01d0, 9.01d0, &
        2.3d0, 2.40d0, 2.60d0, 3.01d0, 4.01d0, 5.01d0, 6.01d0, 7.01d0, 8.01d0, 9.01d0, &
        2.3d0, 2.40d0, 2.60d0, 3.01d0, 4.01d0, 5.01d0, 6.01d0, 7.01d0, 8.01d0, 9.01d0, &
        2.3d0, 2.40d0, 2.60d0, 3.01d0, 4.01d0, 5.01d0, 6.01d0, 7.01d0, 8.01d0, 9.01d0 ], [nlmax, nlmax]), &
        dslope(0:nlmax-1,1:nlmax) = reshape([ ( -0.1d0, l=1, nlmax*nlmax ) ], [nlmax, nlmax])

    if (nl <= 0 .or. nl > nlmax) call strux_fail('maked_lmto47: fit data not available for requested nl')
    if (l < 0 .or. l > nlmax-1) call strux_fail('maked_lmto47: fit data not available for requested l')

    dd = dd0(l, nl) + kap2 * dslope(l, nl)
    dddot = dslope(l, nl)
end subroutine maked_lmto47

subroutine afromd_lmto47(kap2, rfit, nl, lmax, tbalph, tbadot)
!> Purpose:
!>   Port of TBSTR `afromd` for LMTO47 `ialpha=1`.
    implicit none
    integer, intent(in) :: nl, lmax
    double precision, intent(in) :: kap2, rfit
    double precision, intent(out) :: tbalph(0:lmax), tbadot(0:lmax)

    integer, parameter :: nlmax = 20
    integer :: l, l2p1, l2m1
    double precision :: fi(-1:nlmax), gi(-1:nlmax), er2, dd, dddot, t1, t2, t3, t4, t5, t6

    er2 = kap2 * rfit * rfit
    call bessl2(er2, -1, lmax+1, fi, gi)
    do l = 0, lmax
        l2p1 = l + l + 1
        l2m1 = l + l - 1
        call maked_lmto47(nl, l, kap2, dd, dddot)
        t1 = (dd + l + 1d0) * fi(l) - l2m1 * fi(l-1)
        t2 = (dd - l) * gi(l) + l2p1 * gi(l+1)
        tbalph(l) = rfit**l2p1 * t1 / t2
        t3 = fi(l+1)*gi(l)/l2p1 + fi(l)*gi(l-1)/l2m1
        t4 = t3 + fi(l+1)*gi(l+1) - fi(l-1)*gi(l-1)
        t5 = (gi(l)*dble(l)/l2p1 - gi(l+1))*fi(l+1)*(l+1) + &
             (fi(l)*dble(l+1)/l2m1 - fi(l-1))*gi(l-1)*l + 0.5d0
        t6 = -dd*dd*t3 - dd*t4 + t5 + dddot/(rfit*rfit)
        tbadot(l) = 0.5d0 * rfit**(l+l+3) * t6/(t2*t2)
    end do
end subroutine afromd_lmto47

subroutine mkalph_lmto47_from_fit(nspec, nl, kap2, rfit, lmxb, ldot, alpha_l, adot_l)
!> Purpose:
!>   Construct LMTO47 `ialpha=1` alpha/adot data from the fitted logarithmic derivative path.
    implicit none
    integer, intent(in) :: nspec, nl
    integer, intent(in) :: lmxb(nspec)
    logical, intent(in) :: ldot
    double precision, intent(in) :: kap2, rfit(nspec)
    double precision, intent(out) :: alpha_l(0:nl-1, nspec), adot_l(0:nl-1, nspec)

    integer :: is, lmx
    double precision :: tbalph(0:nl-1), tbadot(0:nl-1)

    alpha_l = 0d0
    adot_l = 0d0
    do is = 1, nspec
        lmx = min(lmxb(is), nl-1)
        call afromd_lmto47(kap2, rfit(is), nl, lmx, tbalph, tbadot)
        alpha_l(0:lmx, is) = tbalph(0:lmx)
        if (ldot) adot_l(0:lmx, is) = tbadot(0:lmx)
    end do
end subroutine mkalph_lmto47_from_fit

subroutine mkalph_lmto47_from_hcr(nspec, nl, kap2, avw, rmt, hcr, lmxb, ldot, alpha_l, adot_l)
!> Purpose:
!>   Construct LMTO47 one-kappa `alpha(l,is)` and `adot(l,is)` from hard-core radii.
!> Notes:
!>   This is a standalone port of the `ialpha=0` branch in TBSTR `mkalph.f`.
!>   Inputs `rmt` and `hcr` are passed in the outward-facing API length unit,
!>   so they are scaled internally by `avw` to match the original formulas.
    implicit none
    integer, intent(in) :: nspec, nl
    integer, intent(in) :: lmxb(nspec)
    logical, intent(in) :: ldot
    double precision, intent(in) :: kap2, avw, rmt(nspec), hcr(nl, nspec)
    double precision, intent(out) :: alpha_l(0:nl-1, nspec), adot_l(0:nl-1, nspec)

    integer :: is, l, lmx
    double precision :: hcr_sc, y
    double precision, allocatable :: fi(:), gi(:)

    alpha_l = 0d0
    adot_l = 0d0
    allocate(fi(-1:nl), gi(-1:nl))

    do is = 1, nspec
        lmx = min(lmxb(is), nl-1)
        do l = 0, lmx
            hcr_sc = hcr(l+1, is) / avw
            y = kap2 * hcr_sc * hcr_sc
            call bessl2(y, -1, l+1, fi, gi)
            alpha_l(l, is) = fi(l) / gi(l) * hcr_sc**(2*l + 1)
            if (ldot) then
                adot_l(l, is) = -0.5d0 * hcr_sc**(2*l + 3) * &
                    (fi(l+1)*gi(l)/(2d0*l + 1d0) + fi(l)*gi(l-1)/(2d0*l - 1d0)) / &
                    (gi(l)*gi(l))
            end if
        end do
    end do

    deallocate(fi, gi)
end subroutine mkalph_lmto47_from_hcr

subroutine expand_screening_lmto47(nbas, nspec, nl, ips, lmxb, alpha_l, adot_l, alpha_out, adot)
!> Purpose:
!>   Expand species/l LMTO47 screening data to site/orbital storage.
!> Notes:
!>   The LMTO47 one-kappa path stores one value per `(l,site)` channel and
!>   replicates it over the corresponding `m` components.
    implicit none
    integer, intent(in) :: nbas, nspec, nl
    integer, intent(in) :: ips(nbas), lmxb(nspec)
    double precision, intent(in) :: alpha_l(0:nl-1, nspec), adot_l(0:nl-1, nspec)
    double precision, intent(out) :: alpha_out(nl*nl, nbas), adot(nl*nl, nbas)

    integer :: ib, is, l, m, lm, lmx

    alpha_out = 1d-6
    adot = 0d0
    do ib = 1, nbas
        is = ips(ib)
        if (is < 1 .or. is > nspec) then
            write(*,'(a,i0,a,i0,a,i0)') 'expand_screening_lmto47: invalid ips(', ib, ')=', is, ' for nspec=', nspec
            stop 1
        end if
        lmx = min(lmxb(is), nl-1)
        lm = 0
        do l = 0, nl-1
            do m = -l, l
                lm = lm + 1
                if (l <= lmx) then
                    alpha_out(lm, ib) = alpha_l(l, is)
                    adot(lm, ib) = adot_l(l, is)
                end if
            end do
        end do
    end do
end subroutine expand_screening_lmto47

subroutine maksp0_lmto47(plat, bas, cg, cy, diagv, iax, indxcg, iop, jcg, kap2, nlsqr, npr, sp0)
!> Purpose:
!>   Build dense LMTO47 `maksp0` matrices for one local cluster.
!> Notes:
!>   This is the dense standalone analog of TBSTR `maksp0.f`.  It assembles
!>   either `alpha^-1 - S^0` (`iop=1`) or `-adot*alpha^-2 - Sdot^0` (`iop=2`)
!>   for the current local cluster layout.
    implicit none
    integer, intent(in) :: iop, nlsqr, npr
    integer, parameter :: niax = 10
    integer, intent(in) :: iax(niax, npr), indxcg(*), jcg(*)
    double precision, intent(in) :: cg(*), cy(*), kap2, diagv(*), bas(3,*), plat(3,3)
    double precision, intent(inout) :: sp0(:,:)

    integer :: ipr, jpr, ilm, idxd, nlsqri, nlsqrj, offi, offj, blkdim
    double precision :: dr(3), dsqr
    double precision, allocatable :: block(:,:), blockd(:,:)
    integer, allocatable :: offs(:)

    sp0 = 0d0
    allocate(offs(0:npr))
    offs(0) = 0
    do ipr = 1, npr
        offs(ipr) = offs(ipr-1) + iax(9, ipr)
    end do
    do ipr = 1, npr
        nlsqri = iax(9, ipr)
        offi = offs(ipr-1)
        do ilm = 1, nlsqri
            sp0(offi+ilm, offi+ilm) = diagv(offi+ilm)
        end do

        do jpr = ipr + 1, npr
            nlsqrj = iax(9, jpr)
            offj = offs(jpr-1)
            dsqr = drr2x(0, plat, bas(1, iax(2,ipr)), bas(1, iax(2,jpr)), &
                iax(3,jpr)-iax(3,ipr), iax(4,jpr)-iax(4,ipr), iax(5,jpr)-iax(5,ipr), dr)
            blkdim = max(nlsqri, nlsqrj)
            allocate(block(blkdim, max(blkdim,1)), blockd(blkdim, max(blkdim,1)))
            block = 0d0
            blockd = 0d0
            if (iop == 1) then
                if (jpr /= npr) then
                    call mstrx2(kap2, dr, nlsqri, nlsqrj, blkdim, cg, indxcg, jcg, cy, 11, block, blockd)
                else
                    call mstrx3(kap2, dr, nlsqri, nlsqrj, blkdim, cg, indxcg, jcg, cy, 1, block, blockd)
                end if
                sp0(offi+1:offi+nlsqri, offj+1:offj+nlsqrj) = block(1:nlsqri, 1:nlsqrj)
                sp0(offj+1:offj+nlsqrj, offi+1:offi+nlsqri) = transpose(block(1:nlsqri, 1:nlsqrj))
            else
                if (jpr /= npr) then
                    call mstrx2(kap2, dr, nlsqri, nlsqrj, blkdim, cg, indxcg, jcg, cy, 12, block, blockd)
                else
                    call mstrx3(kap2, dr, nlsqri, nlsqrj, blkdim, cg, indxcg, jcg, cy, 2, block, blockd)
                end if
                sp0(offi+1:offi+nlsqri, offj+1:offj+nlsqrj) = blockd(1:nlsqri, 1:nlsqrj)
                sp0(offj+1:offj+nlsqrj, offi+1:offi+nlsqri) = transpose(blockd(1:nlsqri, 1:nlsqrj))
            end if
            deallocate(block, blockd)
        end do
    end do
    deallocate(offs)
end subroutine maksp0_lmto47

subroutine mktral_lmto47(nl2, nbas, alpha, adot, tral, trad)
!> Purpose:
!>   Construct the LMTO47 `itrans=0` transfer matrices.
!> Notes:
!>   Initial `lmto47` support is restricted to `itrans=0`, where the
!>   transformation is diagonal and identical to the original TBSTR formulas:
!>     tral = [1, 0, -alpha, 1]
!>     trad = [0, 0, -adot, 0]
    implicit none
    integer, intent(in) :: nl2, nbas
    double precision, intent(in) :: alpha(nl2, nbas), adot(nl2, nbas)
    double precision, intent(out) :: tral(4, nl2, nbas), trad(4, nl2, nbas)

    integer :: ib, ilm

    tral = 0d0
    trad = 0d0
    do ib = 1, nbas
        do ilm = 1, nl2
            tral(1, ilm, ib) = 1d0
            tral(3, ilm, ib) = -alpha(ilm, ib)
            tral(4, ilm, ib) = 1d0
            trad(3, ilm, ib) = -adot(ilm, ib)
        end do
    end do
end subroutine mktral_lmto47

subroutine strux_lmto47_screening(nbas, nspec, nl, avw, ips, lmxb, rmt, alpha_in, &
    alpha_out, adot, tral, trad, hcr, screening_mode)
!> Purpose:
!>   Build LMTO47-style screening data for user-supplied `alpha(l,s)`.
!> Notes:
!>   This corresponds to the TBSTR `ialpha=2` policy in a standalone-friendly
!>   form: `alpha_in` is taken from the caller, while `adot` is derived from
!>   hard-core radii when derivative output is requested.
    implicit none
    integer, intent(in) :: nbas, nspec, nl
    double precision, intent(in) :: avw, rmt(nspec), alpha_in(0:nl-1, nspec)
    integer, intent(in) :: ips(nbas), lmxb(nspec)
    double precision, intent(in), optional :: hcr(nl, nspec)
    integer, intent(in), optional :: screening_mode
    double precision, intent(inout), allocatable :: alpha_out(:,:), adot(:,:)
    double precision, intent(inout), allocatable :: tral(:,:,:), trad(:,:,:)
    integer :: nl2, is, l, lmx, ialpha_mode
    double precision :: maxerr
    double precision, parameter :: alpha_tol = 1d-6
    double precision, allocatable :: alpha_l(:,:), alpha_ref(:,:), adot_l(:,:)

    nl2 = nl*nl
    allocate(alpha_out(nl2, nbas), adot(nl2, nbas))
    allocate(tral(4, nl2, nbas), trad(4, nl2, nbas))
    allocate(alpha_l(0:nl-1, nspec), adot_l(0:nl-1, nspec), alpha_ref(0:nl-1, nspec))

    alpha_l = alpha_in
    adot_l = 0d0
    ialpha_mode = 2
    if (present(screening_mode)) ialpha_mode = screening_mode
    if (ialpha_mode == 1) then
        call mkalph_lmto47_from_fit(nspec, nl, 0d0, rmt/avw, lmxb, .true., alpha_l, adot_l)
    else if (present(hcr)) then
        call mkalph_lmto47_from_hcr(nspec, nl, 0d0, avw, rmt, hcr, lmxb, .true., alpha_ref, adot_l)
        maxerr = 0d0
        do is = 1, nspec
            lmx = min(lmxb(is), nl-1)
            do l = 0, lmx
                maxerr = max(maxerr, abs(alpha_ref(l, is) - alpha_in(l, is)))
            end do
        end do
        if (maxerr > alpha_tol) then
            call strux_fail('strux_lmto47: manual-alpha Sdot requires hcr values consistent with alpha_in; use autoalpha or supply matching hard-core radii')
        end if
        alpha_l = alpha_in
    end if

    call expand_screening_lmto47(nbas, nspec, nl, ips, lmxb, alpha_l, adot_l, alpha_out, adot)
    call mktral_lmto47(nl2, nbas, alpha_out, adot, tral, trad)
    deallocate(alpha_l, alpha_ref, adot_l)
end subroutine strux_lmto47_screening

subroutine strux_lmto47_autoalpha_screening(nbas, nspec, nl, avw, ips, lmxb, rmt, hcr, &
    alpha_l_out, alpha_out, adot, tral, trad, screening_mode)
!> Purpose:
!>   Build LMTO47-style screening data from hard-core radii.
!> Notes:
!>   This is the standalone front-end analog of the LMTO47 `ialpha=0` path.
    implicit none
    integer, intent(in) :: nbas, nspec, nl
    double precision, intent(in) :: avw, rmt(nspec)
    double precision, intent(in), optional :: hcr(nl, nspec)
    integer, intent(in) :: ips(nbas), lmxb(nspec)
    integer, intent(in), optional :: screening_mode
    double precision, intent(inout), allocatable :: alpha_l_out(:,:)
    double precision, intent(inout), allocatable :: alpha_out(:,:), adot(:,:)
    double precision, intent(inout), allocatable :: tral(:,:,:), trad(:,:,:)

    integer :: nl2, ialpha_mode
    double precision, allocatable :: adot_l(:,:)

    if (allocated(alpha_l_out)) deallocate(alpha_l_out)
    allocate(alpha_l_out(0:nl-1, nspec), adot_l(0:nl-1, nspec))
    ialpha_mode = 0
    if (present(screening_mode)) ialpha_mode = screening_mode
    select case (ialpha_mode)
    case (0)
        if (.not. present(hcr)) call strux_fail('strux_lmto47_autoalpha_screening: hcr is required for screening_mode=0')
        call mkalph_lmto47_from_hcr(nspec, nl, 0d0, avw, rmt, hcr, lmxb, .true., alpha_l_out, adot_l)
    case (1)
        call mkalph_lmto47_from_fit(nspec, nl, 0d0, rmt/avw, lmxb, .true., alpha_l_out, adot_l)
    case default
        call strux_fail('strux_lmto47_autoalpha_screening: unsupported screening_mode')
    end select
    nl2 = nl*nl
    allocate(alpha_out(nl2, nbas), adot(nl2, nbas))
    allocate(tral(4, nl2, nbas), trad(4, nl2, nbas))
    call expand_screening_lmto47(nbas, nspec, nl, ips, lmxb, alpha_l_out, adot_l, alpha_out, adot)
    call mktral_lmto47(nl2, nbas, alpha_out, adot, tral, trad)
    deallocate(adot_l)
end subroutine strux_lmto47_autoalpha_screening

subroutine scals_lmto47(iax, npr, nbas, nl2, nenv, ndimW, ldot, tral, trad, sc, sdotc)
!> Purpose:
!>   Apply the LMTO47 `scals` transformation to one local cluster solve.
!> Notes:
!>   This is a standalone port of the one-kappa branch of TBSTR `scals.f`,
!>   adapted to the package's explicit dense matrix storage.  The output keeps
!>   only the physical head block columns (`nenv`), exactly matching how the
!>   screened blocks are later scattered by `addtos`.
    implicit none
    integer, intent(in) :: npr, nbas, nl2, nenv, ndimW
    logical, intent(in) :: ldot
    integer, parameter :: niax = 10
    integer, intent(in) :: iax(niax, npr)
    double precision, intent(in) :: tral(4, nl2, nbas), trad(4, nl2, nbas)
    double precision, intent(inout) :: sc(ndimW, nenv), sdotc(ndimW, nenv)

    integer :: ipr, ibas, ic, ilm, li, row, lm2, head_site
    double precision :: dt, w1, w2, w3, w4, hw2, hw3

    row = 0
    head_site = iax(2,1)
    do ipr = 1, npr
        ibas = iax(2, ipr)
        ic = ibas
        do ilm = 1, iax(9, ipr)
            row = row + 1
            li = ll(ilm)
            w1 = tral(1, ilm, ic) / tral(3, ilm, ic)
            w2 = 1d0 / tral(3, ilm, ic)
            if (ldot) then
                w3 = trad(3, ilm, ic) / tral(3, ilm, ic)
                w4 = (-trad(3, ilm, ic) * tral(1, ilm, ic) / tral(3, ilm, ic) + &
                    trad(1, ilm, ic)) / tral(3, ilm, ic)
            else
                w3 = 0d0
                w4 = 0d0
            end if

            do lm2 = 1, nenv
                dt = tral(1, lm2, head_site) * tral(4, lm2, head_site) - &
                     tral(2, lm2, head_site) * tral(3, lm2, head_site)
                hw2 = 1d0 / tral(3, lm2, head_site)
                sc(row,lm2) = w2 * sc(row,lm2) * hw2 * dt
                if (ldot) then
                    hw3 = trad(3, lm2, head_site) / tral(3, lm2, head_site)
                    sdotc(row,lm2) = -w2 * sdotc(row,lm2) * hw2 * dt - sc(row,lm2) * hw3 - w3 * sc(row,lm2)
                end if
            end do
            if (row <= nenv) then
                sc(row,row) = sc(row,row) + w1
                if (ldot) sdotc(row,row) = sdotc(row,row) + w4
            end if
        end do
    end do
end subroutine scals_lmto47

subroutine salph1_lmto47_direct(nbas, nenv, nl2, ldot, ndimW, lmaxw, npr, plat, bas, alphv, adotv, &
    iax, cy, cg, indxcg, jcg, el, tral, trad, balph, bdot)
!> Purpose:
!>   Solve the LMTO47 one-kappa screened cluster problem for one center.
!> Notes:
!>   This follows the LMTO47 `salph1` structure:
!>     build `alpha^-1 - S^0`
!>     solve for the inverse head block
!>     build `-adot*alpha^-2 - Sdot^0`
!>     solve for `sdotc`
!>     apply LMTO47 `scals`
!>   The dense standalone implementation replaces packed-storage LINPACK with
!>   LAPACK factorization on explicit symmetric matrices.
    implicit none
    integer, intent(in) :: nbas, nenv, nl2, ndimW, lmaxw, npr
    logical, intent(in) :: ldot
    integer, parameter :: niax = 10
    integer, intent(in) :: iax(niax, npr)
    double precision, intent(in) :: plat(3,3), bas(3,nbas), el
    double precision, intent(in) :: alphv(ndimW), adotv(ndimW)
    double precision, intent(in) :: tral(4, nl2, nbas), trad(4, nl2, nbas)
    double precision, intent(in) :: cy(100), cg(1200)
    integer, intent(in) :: jcg(1200), indxcg(350)
    double precision, intent(out) :: balph(ndimW, nenv), bdot(ndimW, nenv)

    double precision, allocatable :: a(:), d(:), s0a(:,:), s0d(:,:), work(:)
    integer, allocatable :: ipiv(:)
    integer :: i, j, info

    allocate(a(ndimW), d(ndimW), s0a(ndimW,ndimW), s0d(ndimW,ndimW))
    s0a = 0d0
    s0d = 0d0
    if (any(abs(alphv) <= tiny(1d0))) then
        call strux_fail('salph1_lmto47_direct: alpha contains zero/underflow entries')
    end if
    a = 1d0 / alphv
    d = 0d0
    if (ldot) d = -adotv / (alphv*alphv)

    call maksp0_lmto47(plat, bas, cg, cy, a, iax, indxcg, 1, jcg, el, ndimW, npr, s0a)
    if (ldot) call maksp0_lmto47(plat, bas, cg, cy, d, iax, indxcg, 2, jcg, el, ndimW, npr, s0d)

    call dpzero(balph, ndimW*nenv)
    do j = 1, nenv
        balph(j,j) = 1d0
    end do

    allocate(ipiv(ndimW), work(max(1, ndimW*64)))
    call dsytrf('L', ndimW, s0a, ndimW, ipiv, work, size(work), info)
    if (info /= 0) call strux_fail('salph1_lmto47_direct: matrix singular')
    call dsytrs('L', ndimW, nenv, s0a, ndimW, ipiv, balph, ndimW, info)
    if (info /= 0) call strux_fail('salph1_lmto47_direct: backsolve failed')

    if (ldot) then
        call dgemm('N','N', ndimW, nenv, ndimW, 1d0, s0d, ndimW, balph, ndimW, 0d0, bdot, ndimW)
        call dsytrs('L', ndimW, nenv, s0a, ndimW, ipiv, bdot, ndimW, info)
        if (info /= 0) call strux_fail('salph1_lmto47_direct: derivative backsolve failed')
    else
        bdot = 0d0
    end if

    call scals_lmto47(iax, npr, nbas, nl2, nenv, ndimW, ldot, tral, trad, balph, bdot)
    deallocate(a, d, s0a, s0d, ipiv, work)
end subroutine salph1_lmto47_direct

subroutine strscr_lmto47(nbas, npadl, npadr, alat, plat, bas, rwats, nl, kap, nkap, lmaxw, &
    cy, cg, indxcg, jcg, ldot, ntab, iax, alpha, adot, tral, trad, s, sdot)
!> Purpose:
!>   LMTO47-style `salph0` orchestration for the supported one-kappa path.
!> Notes:
!>   Supported scope:
!>     - one-kappa only
!>     - `itrans=0`
!>     - optional Watson sphere augmentation
!>   Existing solver families remain untouched; this branch is additive.
    implicit none
    integer, intent(in) :: nbas, npadl, npadr, nl, nkap, lmaxw
    logical, intent(in) :: ldot
    integer, parameter :: niax = 10
    integer, intent(inout) :: iax(niax,*)
    integer, intent(in) :: ntab(nbas+1)
    double precision, intent(in) :: alat, plat(3,3), bas(3,nbas), rwats(nbas), kap(nkap)
    double precision, intent(in) :: cy(100), cg(1200)
    integer, intent(in) :: jcg(1200), indxcg(350)
    double precision, intent(in) :: alpha(nl**2,*), adot(nl**2,*)
    double precision, intent(in) :: tral(4,nl**2,*), trad(4,nl**2,*)
    double precision, intent(inout) :: s(nl*nl, nl*nl, 1, 1, *), sdot(nl*nl, nl*nl, 1, 1, *)

    integer :: i, ik, iat, ndimL, ndimW, nl2, ns, nsk, nitab
    integer :: offik, nbasp, nbaspp, nttab, nclus, iclus, nenv
    double precision :: xx
    real(8), allocatable :: salph(:,:,:,:), alphv(:), adotv(:), sadot(:,:,:,:)
    integer, allocatable :: sflg(:), sflgd(:)

    if (nkap /= 1) call strux_fail('strscr_lmto47: only one-kappa LMTO47 is implemented')

    nbasp  = nbas + npadl + npadr
    nbaspp = 2*nbasp - nbas
    nttab  = ntab(nbasp+1)
    nl2    = nl**2

    call dscal(9,        alat, plat, 1)
    call dscal(3*nbaspp, alat, bas,  1)

    nsk = 0
    do ik = 1, nkap
        nitab = 0
        offik = 1 + (ik-1)*nbasp

        do iat = 1, nbasp
            nclus = ntab(iat+1) - ntab(iat)
            iclus = ntab(iat) + 1
            if (iclus < 1 .or. iclus > nttab) then
                write(*,'(a,i0,a,i0,a,i0,a,i0)') 'STRUX dbg bad iclus iat=', iat, ' iclus=', iclus, ' nttab=', nttab, ' nclus=', nclus
                stop 1
            end if
            ndimL = isum(nclus, iax(9,iclus), niax)
            ndimW = ndimL
            if (lmaxw >= 0) ndimW = ndimW + (lmaxw+1)**2

            do i = iclus, ntab(iat+1)
                nitab = nitab + 1
                iax(8,i) = nitab
            end do

            allocate(alphv(ndimW), salph(ndimW,1,nl2,1))
            if (ldot) then
                allocate(adotv(ndimW), sadot(ndimW,1,nl2,1))
                call dpzero(adotv, ndimW)
                call dpzero(sadot, ndimW*nl2)
            else
                allocate(adotv(ndimW), sadot(ndimW,1,nl2,1))
                call dpzero(adotv, ndimW)
                call dpzero(sadot, ndimW*nl2)
            end if

            call pack_screening_lmto47(alpha(1,offik), adot(1,offik), iax(1,iclus), nl2, nbaspp, 1, nclus, ndimW, alphv, adotv)
            if (lmaxw >= 0) call alwats(ndimL, kap(ik), 2, lmaxw, rwats(iat), alphv, adotv)

            nenv = iax(9, iclus)
            call salph1_lmto47_direct(nbaspp, nenv, nl2, ldot, ndimW, lmaxw, nclus, &
                plat, bas, alphv, adotv, iax(1,iclus), cy, cg, indxcg, jcg, kap(ik), &
                tral, trad, salph(:,1,:,1), sadot(:,1,:,1))

            ns = nitab - nclus
            call addtos(nl2, ndimW, iax(1,iclus), nl2, 1, nclus, salph, ns+nsk, s)
            if (ldot) call addtos(nl2, ndimW, iax(1,iclus), nl2, 1, nclus, sadot, ns+nsk, sdot)

            deallocate(alphv, salph, adotv, sadot)
        end do

        allocate(sflg(nitab));  call iinit(sflg,  nitab)
        allocate(sflgd(nitab)); call iinit(sflgd, nitab)
        call symstr(0, nl2, nttab, iax, 1, 1, sflg,  s(1,1,1,1,1+nsk), xx)
        if (ldot) call symstr(0, nl2, nttab, iax, 1, 1, sflgd, sdot(1,1,1,1,1+nsk), xx)
        nsk = nsk + ns + nclus
        deallocate(sflg, sflgd)
    end do

    call dscal(9,        1d0/alat, plat, 1)
    call dscal(3*nbaspp, 1d0/alat, bas,  1)
end subroutine strscr_lmto47

subroutine pack_screening_lmto47(alpha, adot, iax, nla, nbas, nkap, npr, nd, alphv, adotv)
!> Purpose:
!>   Pack physical LMTO47 screening data for one local cluster.
!> Notes:
!>   This helper returns the physical `alpha` and `adot` values expected by
!>   `salph1_lmto47_direct`.
    implicit none
    integer, intent(in) :: nla, nbas, nkap, npr, nd
    integer, parameter :: niax = 10
    integer, intent(in) :: iax(niax, npr)
    double precision, intent(in) :: alpha(nla, nbas, nkap, nkap), adot(nla, nbas, nkap, nkap)
    double precision, intent(out) :: alphv(nd), adotv(nd)
    integer :: ipr, ib, ilm, iv, nlmi

    if (nkap /= 1) call strux_fail('pack_screening_lmto47: only one-kappa mode is implemented')
    alphv = 0d0
    adotv = 0d0
    iv = 0
    do ipr = 1, npr
        ib = iax(2, ipr)
        nlmi = iax(9, ipr)
        do ilm = 1, nlmi
            iv = iv + 1
            alphv(iv) = alpha(ilm, ib, 1, 1)
            adotv(iv) = adot(ilm, ib, 1, 1)
        end do
    end do
end subroutine pack_screening_lmto47

end module strux_lmto47_mode_solver
