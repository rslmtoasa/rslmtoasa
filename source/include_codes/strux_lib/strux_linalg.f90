module strux_matrix_utils
!- Small matrix helpers kept separate from BLAS/LAPACK kernels.
    implicit none
    private

    public :: dmadd, dmcpy, dmscop

contains

subroutine dmadd(a, nca, nra, scalea, b, ncb, nrb, scaleb, c, ncc, nrc, n, m)
!> Purpose:
!>   Compute `c = scalea*a + scaleb*b` for matrices in legacy strided storage.
!> Notes:
!>   The routine preserves the original pointer-arithmetic style so old kernels
!>   can be ported directly without reshaping their data.
    implicit none
    integer, intent(in) :: nca, nra, ncb, nrb, ncc, nrc, n, m
    double precision, intent(in) :: a(0:nra*nca-1), b(0:nrb*ncb-1)
    double precision, intent(inout) :: c(0:nrc*ncc-1)
    double precision, intent(in) :: scalea, scaleb
    integer :: i, j, ia, ib, ic

    if (nra == 1 .and. nrb == 1 .and. nrc == 1) then
        if (scalea == 0d0) then
            do j = 0, m - 1
                ia = nca * j
                ib = ncb * j
                ic = ncc * j
                do i = 0, n - 1
                    c(i + ic) = b(i + ib) * scaleb
                end do
            end do
            return
        else if (scaleb == 0d0) then
            do j = 0, m - 1
                ia = nca * j
                ib = ncb * j
                ic = ncc * j
                do i = 0, n - 1
                    c(i + ic) = a(i + ia) * scalea
                end do
            end do
            return
        else if (scalea == 1d0 .and. scaleb == 1d0) then
            do j = 0, m - 1
                ia = nca * j
                ib = ncb * j
                ic = ncc * j
                do i = 0, n - 1
                    c(i + ic) = a(i + ia) + b(i + ib)
                end do
            end do
            return
        else if (scalea == 1d0 .and. scaleb == -1d0) then
            do j = 0, m - 1
                ia = nca * j
                ib = ncb * j
                ic = ncc * j
                do i = 0, n - 1
                    c(i + ic) = a(i + ia) - b(i + ib)
                end do
            end do
            return
        else
            do j = 0, m - 1
                ia = nca * j
                ib = ncb * j
                ic = ncc * j
                do i = 0, n - 1
                    c(i + ic) = a(i + ia) * scalea + b(i + ib) * scaleb
                end do
            end do
            return
        end if
    end if

    if (scalea /= 0d0 .and. scaleb /= 0d0) then
        do i = n - 1, 0, -1
            ia = i * nra + m * nca
            ib = i * nrb + m * ncb
            ic = i * nrc + m * ncc
            do j = m - 1, 0, -1
                ia = ia - nca
                ib = ib - ncb
                ic = ic - ncc
                c(ic) = a(ia) * scalea + b(ib) * scaleb
            end do
        end do
    else if (scalea /= 0d0) then
        do i = n - 1, 0, -1
            ia = i * nra + m * nca
            ic = i * nrc + m * ncc
            do j = m - 1, 0, -1
                ia = ia - nca
                ic = ic - ncc
                c(ic) = a(ia) * scalea
            end do
        end do
    else if (scaleb /= 0d0) then
        do i = n - 1, 0, -1
            ib = i * nrb + m * ncb
            ic = i * nrc + m * ncc
            do j = m - 1, 0, -1
                ib = ib - ncb
                ic = ic - ncc
                c(ic) = b(ib) * scaleb
            end do
        end do
    end if
end subroutine dmadd

subroutine dmcpy(a, nca, nra, b, ncb, nrb, n, m)
!> Purpose:
!>   Copy a matrix between two legacy strided storage layouts.
    implicit none
    integer, intent(in) :: nca, nra, ncb, nrb, n, m
    double precision, intent(in) :: a(0:*)
    double precision, intent(out) :: b(0:*)
    integer :: i, j, ia, ib

    if (nra == 1 .and. nrb == 1) then
        do j = 0, m - 1
            ia = j * nca
            ib = j * ncb
            do i = 0, n - 1
                b(i + ib) = a(i + ia)
            end do
        end do
        return
    end if

    do i = n - 1, 0, -1
        ia = i * nra + m * nca
        ib = i * nrb + m * ncb
        do j = m - 1, 0, -1
            ia = ia - nca
            ib = ib - ncb
            b(ib) = a(ia)
        end do
    end do
end subroutine dmcpy

subroutine dmscop(dest, ldd, src, lds, n1rs, n2rs, n1cs, n2cs, n1rd, n1cd, fac)
!> Purpose:
!>   Copy one rectangular block into another matrix with an optional scale factor.
!> Notes:
!>   This is the glue routine used to scatter pair-expansion blocks into the
!>   packed cluster matrices.
    implicit none
    integer, intent(in) :: ldd, lds, n1rs, n2rs, n1cs, n2cs, n1rd, n1cd
    double precision, intent(inout) :: dest(ldd, 1)
    double precision, intent(in) :: src(lds, 1), fac
    integer :: i, j, offr, offc

    offr = n1rd - n1rs
    offc = n1cd - n1cs
    if (fac == 1d0) then
        do j = n1cs, n2cs
            do i = n1rs, n2rs
                dest(i + offr, j + offc) = src(i, j)
            end do
        end do
    else
        do j = n1cs, n2cs
            do i = n1rs, n2rs
                dest(i + offr, j + offc) = fac * src(i, j)
            end do
        end do
    end if
end subroutine dmscop

end module strux_matrix_utils
