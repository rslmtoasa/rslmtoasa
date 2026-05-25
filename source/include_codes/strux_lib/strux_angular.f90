module strux_harmonics
!- Real spherical-harmonic helpers used in the strux kernels.
    implicit none
    private

    public :: sylm

contains

subroutine sylm(r, yl, lmx, r2s)
!- Evaluate packed real spherical harmonics for one vector.
!  This is the single-vector harmonic evaluator used in the retained pair
!  kernels that form one block at a time.
    implicit none
    integer, intent(in) :: lmx
    double precision, intent(in) :: r(3)
    double precision, intent(out) :: yl(*), r2s
    integer :: i, l, lav, lavml, lavmm, lavpl, lavpm, lm1, lmm, lp1, m, mp1, n, nt
    double precision :: r2, st, x, y, z, z2
    double precision :: c(15), s(15), p(15,15)
    equivalence (x, c(2)), (y, s(2)), (z, p(2,1))

    c(1) = 1d0
    s(1) = 0d0
    p(1,1) = 1d0
    p(2,2) = 1d0

    n = (lmx + 1)**2
    yl(1) = 1d0
    x = r(1)
    y = r(2)
    z = r(3)
    st = x * x + y * y
    z2 = z * z
    r2 = st + z2
    r2s = r2
    if (n < 2) return
    if (r2 <= 1d-28) then
        do i = 2, n
            yl(i) = 0d0
        end do
        return
    end if

    yl(2) = y
    yl(3) = z
    yl(4) = x
    nt = 1
    do l = 2, lmx
        lp1 = l + 1
        lm1 = l - 1
        lav = l * lp1 + 1
        p(lp1,1) = ((l + lm1) * z * p(l,1) - lm1 * r2 * p(lm1,1)) / l
        yl(lav) = p(lp1,1)
        nt = nt + 2
        p(lp1,lp1) = p(l,l) * nt
        c(lp1) = x * c(l) - y * s(l)
        s(lp1) = x * s(l) + y * c(l)
        lavpl = lav + l
        yl(lavpl) = p(lp1,lp1) * c(lp1)
        lavml = lav - l
        yl(lavml) = p(lp1,lp1) * s(lp1)
        if (st > z2) then
            do m = 1, lm1
                mp1 = m + 1
                lavpm = lav + m
                lavmm = lav - m
                p(lp1,mp1) = ((lm1 + m) * r2 * p(l,m) - (lp1 - m) * z * p(lp1,m)) / st
                yl(lavpm) = p(lp1,mp1) * c(mp1)
                yl(lavmm) = p(lp1,mp1) * s(mp1)
            end do
        else
            do lmm = 1, lm1
                m = l - lmm
                lavpm = lav + m
                lavmm = lav - m
                mp1 = m + 1
                p(lp1,mp1) = (r2 * (l + m) * p(l,mp1) - st * p(lp1,mp1+1)) / (z * (l - m))
                yl(lavpm) = p(lp1,mp1) * c(mp1)
                yl(lavmm) = p(lp1,mp1) * s(mp1)
            end do
        end if
    end do
end subroutine sylm

end module strux_harmonics

module strux_gaunt
!- Real-harmonic Gaunt and Clebsch-Gordan tables.
    implicit none
    private

    public :: scg

contains

subroutine scg(lmax, cg, indxc, l3cg)
!- Convenience wrapper that builds real-harmonic Gaunt tables in the default mode.
    implicit none
    integer lmax
    integer indxc(*), l3cg(*)
    double precision cg(*)
    integer lnjcg, lnxcg
    call scg0(1, lmax, cg, indxc, l3cg, lnjcg, lnxcg)
end subroutine scg

subroutine scg0(mode, lmax, cg, indxc, l3cg, lnjcg, lnxcg)
!- Generate Gaunt/Clebsch-Gordan coupling tables for real spherical harmonics.
!  The packed output arrays reproduce the historical Questaal/TBSTR layout so
!  the downstream kernels can be kept very close to the originals.
    implicit none
    integer mode, lmax, lnjcg, lnxcg
    integer indxc(*), l3cg(*)
    double precision cg(*)
    integer i, i1, i2, i3, i31, i32, ic, j1, j1s, j2, j2s, k2, l1, l2, l3, lmindx,&
    &m1, m2, m3, mb, n1, n2, n3, nl, nm3, s1, s2, s3, t1, t2, t3
    double precision q1, sr2, t, srpi, fs
    double precision fac(161)
    data srpi/1.7724538509055159d0/
    fs(i) = 1 + 4*(i/2) - 2*i

    lnjcg = 0
    lnxcg = 0
    mb = 999999
    nl = lmax + 1
    sr2 = dsqrt(2d0)
    fac(1) = 1d0
    do i = 1, 160
        fac(i + 1) = i*fac(i)
    end do
    ic = 0
    lmindx = 0
    do i1 = 1, nl
        l1 = i1 - 1
        j1s = 2*l1 + 1
        do j1 = 1, j1s
            m1 = j1 - i1
            n1 = iabs(m1)
            s1 = 0
            if (m1 < 0) s1 = 1
            t1 = 0
            if (m1 == 0) t1 = 1
            do i2 = 1, i1
                l2 = i2 - 1
                i31 = l1 - l2 + 1
                i32 = l1 + l2 + 1
                j2s = 2*l2 + 1
                k2 = j1s*j2s
                if (i2 == i1) j2s = j1
                do j2 = 1, j2s
                    lmindx = lmindx + 1
                    if (mode == 1) indxc(lmindx) = ic + 1
                    m2 = j2 - i2
                    n2 = iabs(m2)
                    s2 = 0
                    if (m2 < 0) s2 = 1
                    t2 = 0
                    if (m2 == 0) t2 = 1
                    if (m1*m2 < 0) then
                        m3 = -n1 - n2
                        mb = -iabs(n1 - n2)
                        if (mb == 0) then
                            nm3 = 1
                        else
                            nm3 = 2
                        end if
                    elseif (m1*m2 == 0) then
                        m3 = m1 + m2
                        nm3 = 1
                    else
                        m3 = n1 + n2
                        mb = iabs(n1 - n2)
                        nm3 = 2
                    end if
1                   continue
                    n3 = iabs(m3)
                    s3 = 0
                    if (m3 < 0) s3 = 1
                    t3 = 0
                    if (m3 == 0) t3 = 1
                    q1 = dsqrt(dble(k2))*fs(n3 + (s1 + s2 + s3)/2)/(2*sr2**(1 + t1 + t2 + t3))
                    do i3 = i31, i32, 2
                        l3 = i3 - 1
                        if (n3 > l3) cycle
                        if (mode == 1) then
                            t = 0d0
                            if (n1 + n2 == -n3) t = t + f102(fac, l1, l2, l3)
                            if (n1 + n2 == n3) t = t + f100(fac, l1, l2, l3, n1, n2, n3)*fs(n3 + s3)
                            if (n1 - n2 == -n3) t = t + f100(fac, l1, l2, l3, n1, -n2, -n3)*fs(n2 + s2)
                            if (n1 - n2 == n3) t = t + f100(fac, l1, l2, l3, -n1, n2, -n3)*fs(n1 + s1)
                        end if
                        ic = ic + 1
                        lnjcg = max(ic, lnjcg)
                        if (mode == 1) then
                            cg(ic) = q1*t*f102(fac, l1, l2, l3)/(srpi*dsqrt(dble(2*l3 + 1)))
                            l3cg(ic) = l3*(l3 + 1) + m3 + 1
                        end if
                    end do
                    nm3 = nm3 - 1
                    m3 = mb
                    if (nm3 > 0) goto 1
                end do
            end do
        end do
    end do
    if (mode == 1) indxc(lmindx + 1) = ic + 1
    lnxcg = max(lnxcg, lmindx + 1)
end subroutine scg0

double precision function f100(fac, j1, j2, j3, m1, m2, m3)
    implicit none
    integer j1, j2, j3, m1, m2, m3
    double precision fac(50)
    integer m, n, n1, n2
    double precision t, t1

    f100 = 0d0
    if (m3 /= m1 + m2) return
    t = (2*j3 + 1)*fac(j1 + j2 - j3 + 1)*fac(j3 + j1 - j2 + 1)*fac(j3 + j2 - j1 + 1)/&
        fac(j1 + j2 + j3 + 2)
    t = dsqrt(t*fac(j1 + m1 + 1)*fac(j1 - m1 + 1)*fac(j2 + m2 + 1)*fac(j2 - m2 + 1)*&
        fac(j3 + m3 + 1)*fac(j3 - m3 + 1))
    n1 = max0(j2 - j3 - m1, j1 - j3 + m2, 0) + 1
    n2 = min0(j1 + j2 - j3, j1 - m1, j2 + m2) + 1
    if (n1 > n2) return
    t1 = 0d0
    do m = n1, n2
        n = m - 1
        t1 = t1 + dble(1 + 4*(n/2) - 2*n)/(fac(m)*fac(j1 + j2 - j3 - n + 1)*fac(j1 - m1 - n + 1)*&
            fac(j2 + m2 - n + 1)*fac(j3 - j2 + m1 + n + 1)*fac(j3 - j1 - m2 + n + 1))
    end do
    f100 = t*t1
end function f100

double precision function f102(fac, l1, l2, l3)
    implicit none
    integer l1, l2, l3
    double precision fac(50)
    integer lt, p, x

    lt = l1 + l2 + l3
    p = lt/2
    f102 = 0d0
    if (2*p /= lt) return
    f102 = dsqrt(dble(2*l3 + 1)/dble(lt + 1))
    f102 = f102*fac(p + 1)/dsqrt(fac(2*p + 1))
    x = p - l1
    f102 = f102*dsqrt(fac(2*x + 1))/fac(x + 1)
    x = p - l2
    f102 = f102*dsqrt(fac(2*x + 1))/fac(x + 1)
    x = p - l3
    f102 = f102*dsqrt(fac(2*x + 1))/fac(x + 1)
    if (x > 2*(x/2)) f102 = -f102
end function f102

end module strux_gaunt

module strux_spherical_norms
!- Normalization constants for real spherical harmonics.
    implicit none
    private

    public :: sylmnc

contains

subroutine sylmnc(c, lmx)
!- Return normalization constants for packed real spherical harmonics.
!  These constants are the `cy` factors used throughout the structure-constant
!  expansion code to convert between polynomial and normalized harmonic forms.
    implicit none
    integer lmx
    double precision c(*)
    integer i, i1, i2, l, lav, lp1, m, n1, n2, n3
    double precision fn2, fpi, tlp1, tpi, y0

    tpi = 8d0*datan(1d0)
    fpi = 2d0*tpi
    y0 = 1d0/dsqrt(fpi)
    c(1) = y0

    do l = 1, lmx
        lp1 = l + 1
        tlp1 = l + lp1
        lav = l*lp1 + 1
        c(lav) = dsqrt(tlp1/fpi)
        do m = 1, l
            n2 = lp1 - m
            n1 = n2 + 1
            n3 = l + m
            fn2 = n2
            do i = n1, n3
                fn2 = fn2*i
            end do
            i1 = lav + m
            i2 = lav - m
            c(i1) = dsqrt(tlp1/(fn2*tpi))
            c(i2) = c(i1)
        end do
    end do
end subroutine sylmnc

end module strux_spherical_norms
