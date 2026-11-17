module strux_geometry_math
!- Geometry helpers for connecting vectors and squared distances.
    implicit none
    private

    public :: drr2x

contains

double precision function drr2x(job, plat, tau1, tau2, i, j, k, dr)
!> Purpose:
!>   Calculate the vector connecting two sites in a solid.
!> Notes:
!>   `tau2 - tau1 + i*a1 + j*a2 + k*a3` is formed in Cartesian coordinates and
!>   its squared length is returned. With `job>0` the direction is flipped to
!>   match the sign conventions used by some legacy structure-constant kernels.
    implicit none
    integer job, i, j, k
    double precision dr(3)
    double precision plat(3, 3), tau1(3), tau2(3)
    integer ix

    drr2x = 0.d0
    do ix = 1, 3
        dr(ix) = tau2(ix) - tau1(ix) + plat(ix, 1)*i + plat(ix, 2)*j + plat(ix, 3)*k
        drr2x = drr2x + dr(ix)**2
    end do
    if (job > 0) dr = -dr
end function drr2x

end module strux_geometry_math

module strux_orbital_indexing
!- Helpers for attaching orbital counts to neighbor tables.
    implicit none
    private

    public :: mkiaxd

contains

subroutine mkiaxd(npr, lmx, ips, iax)
!> Purpose:
!>   Store the orbital count for each pair block in `iax(9,:)`.
!> Notes:
!>   The `iax` pair table stores geometry and bookkeeping in a compact integer
!>   layout. Column 9 is reserved for the basis size `(lmax+1)^2` on the
!>   destination site of each pair.
    implicit none
    integer npr, niax, lmx(*), ips(*)
    parameter(niax=10)
    integer iax(niax, npr)
    integer nli, is, i, jbas
    do i = 1, npr
        jbas = iax(2, i)
        if (jbas < 1 .or. jbas > npr) then
            write(*,'(a,i0,a,i0)') 'mkiaxd: invalid jbas for pair i=', i, ' jbas=', jbas
            stop 1
        end if
        is = ips(jbas)
        if (is < 1) then
            write(*,'(a,i0,a,i0,a,i0)') 'mkiaxd: invalid species index at pair i=', i, ' jbas=', jbas, ' ips=', is
            stop 1
        end if
        nli = (lmx(is) + 1)**2
        iax(9, i) = nli
    end do
end subroutine mkiaxd

end module strux_orbital_indexing

module strux_pair_table
!- Neighbor table construction with the package iax/ntab layout.
    implicit none
    private

    public :: strux_pairs

contains

subroutine strux_pairs(nbas, nspec, alat, plat, pos, rmax, ips, lmxb, &
    nttab, ntab, iax)
!> Purpose:
!>   Enumerate all pair blocks within the real-space cutoff.
!> Notes:
!>   The output follows the historical `iax`/`ntab` layout used by the LMTO
!>   screening and Dyson solvers, including on-site blocks and reverse-pair
!>   links.
    implicit none
    integer,  intent(in)  :: nbas, nspec
    double precision, intent(in) :: alat, plat(3,3), pos(3,nbas), rmax
    integer,  intent(in)  :: ips(nbas), lmxb(nspec)
    integer,  intent(out) :: nttab
    integer,  intent(out) :: ntab(nbas+1)
    integer,  intent(inout), allocatable :: iax(:,:)

    integer,  parameter :: niax = 10
    integer,  parameter :: nmax_init = 50000
    double precision :: rmax2, avec(3,3), dr(3), dsq
    integer  :: ib, jb, n1, n2, n3, ipair
    integer  :: n1lo, n1hi, n2lo, n2hi, n3lo, n3hi
    integer  :: i, j, k
    integer, allocatable :: iax_tmp(:,:)
    integer, allocatable :: npair_site(:)
    integer :: nalloc

    do i = 1, 3
        do j = 1, 3
            avec(i,j) = plat(i,j) * alat
        end do
    end do
    rmax2 = rmax * rmax

    n1hi = ceiling(rmax / sqrt(avec(1,1)**2 + avec(2,1)**2 + avec(3,1)**2) + 1)
    n2hi = ceiling(rmax / sqrt(avec(1,2)**2 + avec(2,2)**2 + avec(3,2)**2) + 1)
    n3hi = ceiling(rmax / sqrt(avec(1,3)**2 + avec(2,3)**2 + avec(3,3)**2) + 1)
    n1lo = -n1hi;  n2lo = -n2hi;  n3lo = -n3hi

    nalloc = nmax_init
    allocate(iax_tmp(niax, nalloc))
    allocate(npair_site(nbas))
    npair_site = 0

    ipair = 0
    do ib = 1, nbas
        if (ips(ib) < 1 .or. ips(ib) > nspec) then
            write(*,'(a,i0,a,i0,a,i0)') 'strux_pairs: invalid ips(', ib, ')=', ips(ib), ' nspec=', nspec
            stop 1
        end if
        ntab(ib) = ipair

        ipair = ipair + 1
        if (ipair > nalloc) then
            nalloc = nalloc * 2
            call resize_iax(iax_tmp, niax, ipair-1, nalloc)
        end if
        iax_tmp(1, ipair) = ib
        iax_tmp(2, ipair) = ib
        iax_tmp(3, ipair) = 0
        iax_tmp(4, ipair) = 0
        iax_tmp(5, ipair) = 0
        iax_tmp(6, ipair) = ipair
        iax_tmp(7, ipair) = ipair
        iax_tmp(8, ipair) = 0
        iax_tmp(9, ipair) = (lmxb(ips(ib))+1)**2
        iax_tmp(10,ipair) = 0

        do jb = 1, nbas
            do n1 = n1lo, n1hi
            do n2 = n2lo, n2hi
            do n3 = n3lo, n3hi
                if (ib == jb .and. n1 == 0 .and. n2 == 0 .and. n3 == 0) cycle
                do k = 1, 3
                    dr(k) = (pos(k,jb) - pos(k,ib))*alat &
                          + avec(k,1)*n1 + avec(k,2)*n2 + avec(k,3)*n3
                end do
                dsq = dr(1)**2 + dr(2)**2 + dr(3)**2
                if (dsq > rmax2) cycle

                ipair = ipair + 1
                if (ipair > nalloc) then
                    nalloc = nalloc * 2
                    call resize_iax(iax_tmp, niax, ipair-1, nalloc)
                end if
                iax_tmp(1, ipair) = ib
                iax_tmp(2, ipair) = jb
                iax_tmp(3, ipair) = n1
                iax_tmp(4, ipair) = n2
                iax_tmp(5, ipair) = n3
                iax_tmp(6, ipair) = 0
                iax_tmp(7, ipair) = ipair
                iax_tmp(8, ipair) = 0
                iax_tmp(9, ipair) = (lmxb(ips(jb))+1)**2
                iax_tmp(10,ipair) = 0
                npair_site(ib) = npair_site(ib) + 1
            end do
            end do
            end do
        end do
    end do
    nttab = ipair
    ntab(nbas+1) = nttab

    do i = 1, nttab
        if (iax_tmp(6,i) /= 0) cycle
        ib = iax_tmp(1,i);  jb = iax_tmp(2,i)
        n1 = iax_tmp(3,i);  n2 = iax_tmp(4,i);  n3 = iax_tmp(5,i)
        do j = ntab(jb)+1, ntab(jb+1)
            if (iax_tmp(1,j) == jb .and. iax_tmp(2,j) == ib .and. &
                iax_tmp(3,j) == -n1 .and. iax_tmp(4,j) == -n2 .and. &
                iax_tmp(5,j) == -n3) then
                iax_tmp(6,i) = j
                iax_tmp(6,j) = i
                exit
            end if
        end do
    end do

    if (allocated(iax)) deallocate(iax)
    allocate(iax(niax, nttab))
    iax(:, 1:nttab) = iax_tmp(:, 1:nttab)
    deallocate(iax_tmp, npair_site)
end subroutine strux_pairs

subroutine resize_iax(arr, niax, nkeep, nnew)
!> Purpose:
!>   Grow a temporary `iax` storage array while preserving the populated prefix.
    implicit none
    integer, intent(in) :: niax, nkeep, nnew
    integer, intent(inout), allocatable :: arr(:,:)
    integer, allocatable :: tmp(:,:)
    allocate(tmp(niax, nnew))
    tmp(:, 1:nkeep) = arr(:, 1:nkeep)
    call move_alloc(tmp, arr)
end subroutine resize_iax

end module strux_pair_table

module strux_expansion_indexing
!- Index builders for real-space expansion tensors.
    use strux_utils, only: ll
    implicit none
    private

    public :: strxsu_setup

contains

subroutine strxsu_setup(nlma, nlb, lb_list, lbx, loka, cg, jcg, indxcg, ncfx, &
    nlmbx, nlmp, npow, lbxmax, cf, ip, ikl, jkl)
!> Purpose:
!>   Build the packed indexing tables used by the real-space pair expansion.
!> Notes:
!>   These tables encode which Gaunt coefficients contribute to each output
!>   block so the expensive bookkeeping is done once before the numerical
!>   kernels run.
    implicit none
    integer, intent(in) :: nlma, nlb, lbx, loka, ncfx, lbxmax
    integer, intent(in) :: lb_list(nlb), jcg(*), indxcg(*)
    integer, intent(in) :: nlmbx, nlmp, npow
    double precision, intent(in) :: cg(*)
    double precision, intent(out) :: cf(ncfx)
    integer, intent(out) :: ip(ncfx), ikl(nlmp, 0:npow)
    integer, intent(out) :: jkl(nlmp, 0:npow, 0:lbxmax)

    integer :: i, klmb, lb, lmaxa, ojkl(0:20)

    lmaxa = ll(nlma)
    do i = 0, 20
        ojkl(i) = -99
    end do
    do i = 1, nlb
        ojkl(lb_list(i)) = -98
    end do

    cf = 0d0
    ip = 0
    ikl = 0
    jkl = 0

    call nstru0(nlma, nlmbx, nlmp, npow, loka, cg, jcg, indxcg, &
        ikl, jkl(:,:,lbx), ncfx, ip, cf)

    do lb = 0, lbx - 1
        if (ojkl(lb) == -98) then
            klmb = (lb + 1)**2
            call nstru1(nlma, nlmp, npow, klmb, jcg, indxcg, ikl, jkl(:,:,lb))
        end if
    end do
end subroutine strxsu_setup

subroutine nstru0(nlma, nlmb, nlmp, npow, loka, cg, jcg, indxcg, ikl, jkl, ncfx, ip, cf)
!> Purpose:
!>   Construct the primary coefficient/index table for one `(l_a,l_b)` expansion family.
!> Notes:
!>   The packed `cf`, `ip`, `ikl`, and `jkl` arrays encode the angular
!>   expansion bookkeeping used by the retained real-space pair kernels.
    implicit none
    integer, intent(in) :: nlma, nlmb, nlmp, npow, loka, ncfx
    integer, intent(in) :: jcg(*), indxcg(*)
    integer, intent(out) :: ikl(nlmp,0:npow), jkl(nlmp,0:npow), ip(*)
    double precision, intent(in) :: cg(*)
    double precision, intent(out) :: cf(ncfx)

    integer, parameter :: lmxx = 20
    integer :: icg, icg1, icg2, ii, ilm, indx, ipow, jj, klm, lk, lm, mlm
    integer :: ntot, l, lmax
    double precision :: fac2l(0:lmxx)

    do ipow = 0, npow
        do ilm = 1, nlmp
            jkl(ilm, ipow) = 0
        end do
    end do

    do klm = 1, nlmb
        lk = ll(klm)
        do mlm = 1, nlma
            lm = ll(mlm)
            ii = max(mlm, klm)
            indx = (ii*(ii-1))/2 + min(mlm, klm)
            icg1 = indxcg(indx)
            icg2 = indxcg(indx+1) - 1
            do icg = icg1, icg2
                ilm = jcg(icg)
                ipow = (lm + lk - ll(ilm))/2
                if (ilm <= nlmp .and. ipow <= npow) jkl(ilm, ipow) = jkl(ilm, ipow) + 1
            end do
        end do
    end do

    ntot = 0
    do ipow = 0, npow
        do ilm = 1, nlmp
            ikl(ilm, ipow) = ntot + 1
            ntot = ntot + jkl(ilm, ipow)
            jkl(ilm, ipow) = ikl(ilm, ipow) - 1
        end do
    end do
    if (ntot > ncfx) then
        write(*,*) 'nstru0: increase ncfx to', ntot
        stop
    end if

    if (loka == 1) then
        lmax = ll(nlma) + ll(nlmb)
        if (lmax > lmxx) then
            write(*,*) 'nstru0: increase lmxx to', lmax
            stop
        end if
        fac2l(0) = 1d0
        do l = 1, lmax
            fac2l(l) = fac2l(l-1) * (2*l - 1)
        end do
    end if

    jj = 0
    do klm = 1, nlmb
        lk = ll(klm)
        do mlm = 1, nlma
            lm = ll(mlm)
            jj = jj + 1
            ii = max(mlm, klm)
            indx = (ii*(ii-1))/2 + min(mlm, klm)
            icg1 = indxcg(indx)
            icg2 = indxcg(indx+1) - 1
            do icg = icg1, icg2
                ilm = jcg(icg)
                ipow = (lm + lk - ll(ilm))/2
                if (ilm <= nlmp .and. ipow <= npow) then
                    jkl(ilm, ipow) = jkl(ilm, ipow) + 1
                    ii = jkl(ilm, ipow)
                    ip(ii) = jj
                    cf(ii) = cg(icg) * (-1d0)**lk
                    if (loka == 1) cf(ii) = cf(ii) * 2d0 / (fac2l(lm) * fac2l(lk))
                end if
            end do
        end do
    end do
end subroutine nstru0

subroutine nstru1(nlma, nlmp, npow, klmb, jcg, indxcg, ikl, jkl)
!> Purpose:
!>   Build secondary index tables for additional destination-shell blocks.
!> Notes:
!>   This reuses the primary layout from `nstru0` when only the destination
!>   angular momentum changes.
    implicit none
    integer, intent(in) :: nlma, nlmp, npow, klmb
    integer, intent(in) :: jcg(*), indxcg(*)
    integer, intent(in) :: ikl(nlmp,0:npow)
    integer, intent(out) :: jkl(nlmp,0:npow)

    integer :: icg, icg1, icg2, ii, ilm, indx, ipow, klm, lk, lm, mlm

    do ipow = 0, npow
        do ilm = 1, nlmp
            jkl(ilm, ipow) = ikl(ilm, ipow) - 1
        end do
    end do

    do mlm = 1, nlma
        lm = ll(mlm)
        do klm = 1, klmb
            lk = ll(klm)
            ii = max(mlm, klm)
            indx = (ii*(ii-1))/2 + min(mlm, klm)
            icg1 = indxcg(indx)
            icg2 = indxcg(indx+1) - 1
            do icg = icg1, icg2
                ilm = jcg(icg)
                ipow = (lm + lk - ll(ilm))/2
                if (ilm <= nlmp .and. ipow <= npow) jkl(ilm, ipow) = jkl(ilm, ipow) + 1
            end do
        end do
    end do
end subroutine nstru1

end module strux_expansion_indexing
