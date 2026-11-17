module strux_solver_assembly
!========================================================================
! High-level assembly of screened structure constants from solver kernels.
!========================================================================
    implicit none
    private

    public :: alwats
    public :: addtos
    public :: symstr

contains

subroutine alwats(ndimL, kap2, ldot, lmaxw, rwats, alphv, adotv)
!> Purpose:
!>   Append Watson-sphere diagonal entries to the packed solver vectors.
!> Notes:
!>   When a Watson sphere is requested, its diagonal contribution is appended
!>   after the physical cluster orbitals in the packed solver vectors.
    use lmto_strux_bessel, only: bessl2
    implicit none
    integer, intent(in) :: lmaxw, ndimL, ldot
    double precision, intent(in) :: kap2, rwats
    double precision, intent(inout) :: alphv(*), adotv(*)

    integer, parameter :: nlmax = 20
    integer :: l, ilm, m, ilmd
    double precision :: fi(-1:nlmax), gi(-1:nlmax), fac, r2, rfac
    double precision :: alphom, adoth

    r2   = rwats**2
    rfac = 1d0/rwats

    call bessl2(kap2*r2, -1, lmaxw+1, fi(-1), gi(-1))
    ilm = 0
    do l = 0, lmaxw
        rfac = rfac*r2
        alphom = fi(l)/gi(l)*rfac
        alphom = 1d0/alphom
        ilmd = ilm
        do m = -l, l
            ilm = ilm + 1
            alphv(ilm + ndimL) = alphom
        end do
        if (ldot /= 0) then
            adoth = -0.5d0*rfac*r2*(fi(l+1)*gi(l)/(l+l+1) + fi(l)*gi(l-1)/(l+l-1))/(gi(l)*gi(l))
            adoth = -adoth*alphom*alphom
            fac = 1d0
            do m = -l, l
                ilmd = ilmd + 1
                if (ldot > 1) fac = 1d0/alphv(ilmd+ndimL)**2
                adotv(ilmd+ndimL) = adoth*fac
            end do
        end if
    end do
end subroutine alwats

subroutine addtos(nds, nda, iax, nenv, nkap, npr, salph, nprs, s)
!> Purpose:
!>   Copy one environment solve back into the global pair-indexed storage.
!> Notes:
!>   `salph` is produced for one local environment. This routine copies those
!>   blocks back into the global pair-indexed storage used by the public API.
    implicit none
    integer, intent(in) :: nds, nda, nenv, nkap, npr, nprs
    integer, parameter  :: niax = 10
    integer, intent(in) :: iax(niax, npr)
    double precision, intent(in)    :: salph(nda, nkap, nenv, nkap)
    double precision, intent(inout) :: s(nds, nds, nkap, nkap, nprs+npr)

    integer :: ipr, ilm, klm, ii, nlma, nlmb, ikap, jkap

    nlmb = iax(9,1)
    ii   = 0
    do ipr = 1, npr
        nlma = iax(9,ipr)
        do ikap = 1, nkap
            do jkap = 1, nkap
                do klm = 1, nlmb
                    do ilm = 1, nlma
                        s(ilm,klm,ikap,jkap,ipr+nprs) = salph(ii+ilm,ikap,klm,jkap)
                    end do
                end do
            end do
        end do
        ii = ii + nlma
    end do
end subroutine addtos

subroutine symstr(mode, nds, nsites, iax, nsp, nkap, sflg, s, asym)
!> Purpose:
!>   Enforce pair-exchange symmetry in the real-space structure-constant table.
!> Notes:
!>   The active standalone library uses only the real, one-kappa path, so the
!>   routine keeps the original symmetry logic but rejects unsupported modes.
    implicit none
    integer, intent(in) :: mode, nds, nsites, nsp, nkap
    integer, parameter  :: niax = 10
    integer, intent(in) :: iax(niax, nsites)
    integer, intent(inout) :: sflg(nsites)
    double precision, intent(inout) :: s(nds, nds, nsp, nkap, nkap, nsites)
    double precision, intent(out)   :: asym

    integer :: i, j, lm1, lm2, isp, ii, jj, ik, jk
    double precision :: tmp

    asym = 0d0
    if (nkap /= 1) return
    ik = 1; jk = 1

    if (mode /= 0) then; write(*,*) 'symstr: only real mode=0 implemented', mode; stop; end if

    do i = 1, nsites
        j = iax(6,i); if (j == 0) cycle
        ii = iax(8,i); if (ii == 0) ii = i
        jj = iax(8,j); if (jj == 0) jj = j
        if (sflg(ii) == 0 .and. sflg(jj) == 0) then
            do isp = 1, nsp; do lm1 = 1, nds; do lm2 = 1, nds
                tmp = (s(lm1,lm2,isp,ik,jk,ii) + s(lm2,lm1,isp,ik,jk,jj))/2
                asym = max(asym, abs(s(lm1,lm2,isp,ik,jk,ii)-tmp))
                s(lm1,lm2,isp,ik,jk,ii) = tmp
                s(lm2,lm1,isp,ik,jk,jj) = tmp
            end do; end do; end do
            sflg(ii) = 1; sflg(jj) = 1
        else if (sflg(ii) == 0) then
            do isp = 1, nsp; do lm1 = 1, nds; do lm2 = 1, nds
                s(lm1,lm2,isp,ik,jk,ii) = s(lm2,lm1,isp,ik,jk,jj)
            end do; end do; end do
            sflg(ii) = 1
        else if (sflg(jj) == 0) then
            do isp = 1, nsp; do lm1 = 1, nds; do lm2 = 1, nds
                s(lm2,lm1,isp,ik,jk,jj) = s(lm1,lm2,isp,ik,jk,ii)
            end do; end do; end do
            sflg(jj) = 1
        end if
    end do
end subroutine symstr

end module strux_solver_assembly
