module strux_utils
!- General scalar, array, and indexing helpers used across strux.
    implicit none
    private

    public :: iinit, dpzero
    public :: dvset
    public :: isum, ll

contains

subroutine iinit(array, leng)
!> Purpose:
!>   Fill an integer array with zeros.
    implicit none
    integer, intent(in) :: leng
    integer, intent(out) :: array(leng)
    integer :: i
    do i = 1, leng
        array(i) = 0
    end do
end subroutine iinit

subroutine dpzero(array, leng)
!> Purpose:
!>   Fill a double-precision real array with zeros.
!> Notes:
!>   The OpenMP pragma is preserved because some of the larger solver scratch
!>   arrays can be expensive to initialize.
    implicit none
    integer, intent(in) :: leng
    double precision, intent(out) :: array(leng)
    integer :: i
    do i = 1, leng
        array(i) = 0d0
    end do
end subroutine dpzero

subroutine dvset(array, i1, i2, val)
!> Purpose:
!>   Set a contiguous slice of a double-precision array to a scalar value.
    implicit none
    integer, intent(in) :: i1, i2
    double precision, intent(inout) :: array(i2)
    double precision, intent(in) :: val
    integer :: i
    do i = i1, i2
        array(i) = val
    end do
end subroutine dvset

integer function isum(n, idx, incx)
!> Purpose:
!>   Sum `n` integer entries from a strided vector.
!> Notes:
!>   This is used heavily with the packed `iax` table where orbital counts live
!>   in every `incx`-th entry.
    implicit none
    integer, intent(in) :: n, incx
    integer, intent(in) :: idx(*)
    integer :: i, nincx
    isum = 0
    if (n <= 0) return
    nincx = n * incx
    do i = 1, nincx, incx
        isum = isum + idx(i)
    end do
end function isum

integer function ll(ilm)
!> Purpose:
!>   Decode angular momentum `l` from the packed real-harmonic index `ilm`.
!> Notes:
!>   The library uses the conventional `(l,m)` packing where the first
!>   `(l+1)^2` entries belong to shells up to angular momentum `l`.
    implicit none
    integer, intent(in) :: ilm
    ll = int(sqrt(real(ilm - 1, 8)))
end function ll

end module strux_utils
