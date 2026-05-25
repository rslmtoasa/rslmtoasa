module strux_pair_expansion
!========================================================================
! Retained real-space pair expansion kernels used by the LMTO47 workflow.
!========================================================================
    use strux_harmonics, only: sylm
    use strux_utils, only: ll
    implicit none
    private

    public :: mstrx2
    public :: mstrx3

contains

subroutine mstrx2(e, dr, nlma, nlmb, ndim, cg, indxcg, jcg, cy, iop, s, sdot)
!> Purpose:
!>   Form one ordinary pair block from a single displacement vector.
!> Notes:
!>   This retained kernel is the generic pair-expansion backend currently used
!>   by the LMTO47 workflow for both `S0` and `S0dot` assembly.
    use lmto_strux_bessel, only: besslr
    implicit none
    integer, intent(in) :: ndim, nlma, nlmb, iop
    integer, intent(in) :: indxcg(*), jcg(*)
    double precision, intent(in) :: e, dr(3), cy(*), cg(*)
    double precision, intent(out) :: s(ndim,nlmb), sdot(ndim,nlmb)

    logical :: loka
    integer, parameter :: lmxx = 11
    integer :: icg, icg1, icg2, ii, mlm, indx, ipow, klm, l, lk, llm, lm, lmax, &
        m, klm0, sig(0:lmxx), iiop
    double precision :: fpi, rfac, sum, r2, sumd, gdot, xx, &
        edot(0:lmxx), efac(0:lmxx), fac2l(0:lmxx), &
        psi(-1:lmxx), phi(-1:lmxx), yl((lmxx+1)**2), hl((lmxx+1)**2), hdot((lmxx+1)**2)

    lmax = ll(nlma) + ll(nlmb)
    fac2l(0) = 1d0
    do l = 1, lmax
        fac2l(l) = fac2l(l-1)*(2*l-1)
    end do
    iiop = mod(iop, 10)
    loka = mod(iop/10, 10) /= 0
    fpi  = 16d0*atan(1d0)
    if (nlma > ndim) then
        write(*,*) 'MSTRX2: nlma gt ndim, need ndim=', nlma
        stop
    end if
    if (lmax > lmxx) then
        write(*,*) 'MSTRX2: lmax too big, need lmax=', lmax
        stop
    end if
    call sylm(dr, yl, lmax, r2)
    if (r2 < 1d-10) return
    call besslr(e*r2, 0, -1, lmax, phi, psi)
    mlm = 0
    rfac = sqrt(r2)
    efac(0) = 1d0
    edot(0) = 0d0
    sig(0) = 1
    do l = 1, lmax
        efac(l) = -e*efac(l-1)
        edot(l) = -l*efac(l-1)
        sig(l) = -sig(l-1)
    end do

    do l = 0, lmax
        gdot = psi(l-1)*r2/2d0
        rfac = rfac/r2
        do m = -l, l
            mlm = mlm + 1
            xx = rfac*cy(mlm)*yl(mlm)
            hl(mlm) = psi(l)*xx
            hdot(mlm) = gdot*xx
        end do
    end do

    if (iiop /= 2) then
        do mlm = 1, nlma
            lm = ll(mlm)
            klm0 = mlm
            if (mlm > nlmb) klm0 = 1
            do klm = klm0, nlmb
                lk = ll(klm)
                sum = 0d0
                ii = max(mlm, klm)
                indx = (ii*(ii-1))/2 + min(mlm, klm)
                icg1 = indxcg(indx)
                icg2 = indxcg(indx+1) - 1
                do icg = icg1, icg2
                    llm = jcg(icg)
                    ipow = (lm+lk-ll(llm))/2
                    sum = sum + cg(icg)*efac(ipow)*hl(llm)
                end do
                s(mlm,klm) = fpi*sig(lk)*sum
                if (loka) s(mlm,klm) = s(mlm,klm)*2d0/(fac2l(lm)*fac2l(lk))
                if (klm <= nlma .and. mlm <= nlmb) s(klm,mlm) = s(mlm,klm)*sig(lk)*sig(lm)
            end do
        end do
    end if

    if (iiop /= 1) then
        do mlm = 1, nlma
            lm = ll(mlm)
            klm0 = mlm
            if (mlm > nlmb) klm0 = 1
            do klm = klm0, nlmb
                lk = ll(klm)
                sumd = 0d0
                ii = max(mlm, klm)
                indx = (ii*(ii-1))/2 + min(mlm, klm)
                icg1 = indxcg(indx)
                icg2 = indxcg(indx+1) - 1
                do icg = icg1, icg2
                    llm = jcg(icg)
                    ipow = (lm+lk-ll(llm))/2
                    sumd = sumd + cg(icg)*(efac(ipow)*hdot(llm) + edot(ipow)*hl(llm))
                end do
                sdot(mlm,klm) = fpi*sumd*sig(lk)
                if (loka) sdot(mlm,klm) = sdot(mlm,klm)*2d0/(fac2l(lm)*fac2l(lk))
                if (klm <= nlma .and. mlm <= nlmb) sdot(klm,mlm) = sdot(mlm,klm)*sig(lk)*sig(lm)
            end do
        end do
    end if
end subroutine mstrx2

subroutine mstrx3(kap2, dr, nlsqri, nlsqrj, nlsqr, cg, indxcg, jcg, cy, iop, s, sdot)
!> Purpose:
!>   Form one Watson-sphere pair block.
    use lmto_strux_bessel, only: bessl2
    implicit none
    integer, intent(in) :: nlsqr, nlsqri, nlsqrj, iop
    integer, intent(in) :: indxcg(*), jcg(*)
    double precision, intent(in) :: kap2, cy(*), cg(*), dr(3)
    double precision, intent(inout) :: s(nlsqr,*), sdot(nlsqr,*)

    integer, parameter :: nlmax = 20, lmaxx = 2*nlmax-2, nlsqrx = (lmaxx+1)*(lmaxx+1)
    integer :: icg, icg1, icg2, ii, ilm, indx, ipow, jlm, klm, l, li, lj, lk
    integer :: lmax, mk, sig(0:lmaxx)
    double precision :: edot(0:lmaxx), efac(0:lmaxx), fac2l(0:lmaxx), &
        fi(0:lmaxx+1), fpi, fdot, gi(0:lmaxx+1), r2, &
        bdot(nlsqrx), bl(nlsqrx), sum, sumd, ylm(nlsqrx)

    fpi = 12.5663706143591725d0
    lmax = ll(nlsqri) + ll(nlsqrj)
    if (lmax > lmaxx) then
        write(*,*) 'MSTRX3: increase lmaxx, need', lmax
        stop
    end if

    fac2l(0) = 1d0
    do l = 1, lmax
        fac2l(l) = fac2l(l-1)*(2*l-1)
    end do

    call sylm(dr, ylm, lmax, r2)
    call bessl2(kap2*r2, 0, lmax+1, fi(0), gi(0))
    klm = 0
    do lk = 0, lmax
        fdot = -fi(lk+1)*r2/(4*lk+2)
        do mk = -lk, lk
            klm = klm + 1
            bl(klm) = fi(lk)*cy(klm)*ylm(klm)/fac2l(lk)
            bdot(klm) = fdot*cy(klm)*ylm(klm)/fac2l(lk)
        end do
    end do

    efac(0) = 1d0
    edot(0) = 0d0
    sig(0) = 1
    do l = 1, lmax
        efac(l) = -kap2*efac(l-1)
        edot(l) = -l*efac(l-1)
        sig(l) = -sig(l-1)
    end do

    if (iop /= 2) then
        do ilm = 1, nlsqri
            li = ll(ilm)
            do jlm = 1, nlsqrj
                lj = ll(jlm)
                ii = max(ilm, jlm)
                indx = (ii*(ii-1))/2 + min(ilm, jlm)
                icg1 = indxcg(indx)
                icg2 = indxcg(indx+1)-1
                sum = 0d0
                do icg = icg1, icg2
                    klm = jcg(icg)
                    lk = ll(klm)
                    ipow = (li-lj+lk)/2
                    sum = sum + cg(icg)*efac(ipow)*bl(klm)*sig(lk)
                end do
                s(ilm,jlm) = fpi*sum*2d0*fac2l(lj)/fac2l(li)
            end do
        end do
    end if
    if (iop /= 1) then
        do ilm = 1, nlsqri
            li = ll(ilm)
            do jlm = 1, nlsqrj
                lj = ll(jlm)
                ii = max(ilm, jlm)
                indx = (ii*(ii-1))/2 + min(ilm, jlm)
                icg1 = indxcg(indx)
                icg2 = indxcg(indx+1)-1
                sumd = 0d0
                do icg = icg1, icg2
                    klm = jcg(icg)
                    lk = ll(klm)
                    ipow = (li-lj+lk)/2
                    sumd = sumd + cg(icg)*sig(lk)*(efac(ipow)*bdot(klm)+edot(ipow)*bl(klm))
                end do
                sdot(ilm,jlm) = fpi*sumd*2d0*fac2l(lj)/fac2l(li)
            end do
        end do
    end if
end subroutine mstrx3

end module strux_pair_expansion
