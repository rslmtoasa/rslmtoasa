module strux_support
!- Shared standalone helpers for API validation, neighbor-table setup,
!  and common output allocation.
    use strux_gaunt, only: scg
    use strux_orbital_indexing, only: mkiaxd
    use strux_pair_table, only: strux_pairs
    use strux_spherical_norms, only: sylmnc
    implicit none
    private

    integer, parameter, public :: lmxcg_max = 4

    public :: strux_fail
    public :: validate_common_inputs
    public :: allocate_sdot_array
    public :: build_common_tables

contains

subroutine build_common_tables(nbas, nspec, alat, plat, pos, rmax, ips, lmxb, &
    nttab, ntab_arr, iax, cy, cg, indxcg, jcg)
!> Purpose:
!>   Build the common geometric and angular tables shared by the public drivers.
!> Notes:
!>   This is the package-level setup step that turns lattice input into pair
!>   tables (`iax`, `ntab_arr`), Gaunt/Clebsch-Gordan coefficients, and
!>   spherical-harmonic normalization constants.
    implicit none
    integer, intent(in) :: nbas, nspec
    double precision, intent(in) :: alat, plat(3,3), pos(3,nbas), rmax
    integer, intent(in) :: ips(nbas), lmxb(nspec)
    integer, intent(out) :: nttab, ntab_arr(nbas+1)
    integer, intent(inout), allocatable :: iax(:,:)
    double precision, intent(out) :: cy(100), cg(1200)
    integer, intent(out) :: indxcg(350), jcg(1200)

    call strux_pairs(nbas, nspec, alat, plat, pos, rmax, ips, lmxb, nttab, ntab_arr, iax)
    call mkiaxd(nttab, lmxb, ips, iax)
    call scg(lmxcg_max, cg, indxcg, jcg)
    call sylmnc(cy, lmxcg_max)
end subroutine build_common_tables

subroutine allocate_sdot_array(ldot, nl2, nttab, sdot_out)
!> Purpose:
!>   Allocate the derivative output array in a way that keeps call sites simple.
!> Notes:
!>   When `ldot` is false a small dummy array is allocated so callers can pass
!>   one allocatable object through the whole workflow without special casing.
    implicit none
    logical, intent(in) :: ldot
    integer, intent(in) :: nl2, nttab
    double precision, intent(inout), allocatable :: sdot_out(:,:,:)

    if (ldot) then
        allocate(sdot_out(nl2, nl2, nttab))
    else
        allocate(sdot_out(1, 1, 1))
    end if
    sdot_out = 0d0
end subroutine allocate_sdot_array

subroutine validate_common_inputs(caller, nbas, nspec, nl, avw, rmax)
!> Purpose:
!>   Validate the basic scalar inputs shared by the public entry points.
!> Notes:
!>   This keeps the top-level API routines short and ensures user-facing errors
!>   are reported with the name of the calling routine.
    implicit none
    character(*), intent(in) :: caller
    integer, intent(in) :: nbas, nspec, nl
    double precision, intent(in) :: avw, rmax

    if (nbas <= 0) call strux_fail(trim(caller)//': nbas must be > 0')
    if (nspec <= 0) call strux_fail(trim(caller)//': nspec must be > 0')
    if (nl <= 0) call strux_fail(trim(caller)//': nl must be > 0')
    if (avw <= 0d0) call strux_fail(trim(caller)//': avw must be > 0')
    if (rmax <= 0d0) call strux_fail(trim(caller)//': rmax must be > 0')
end subroutine validate_common_inputs

subroutine strux_fail(msg)
!> Purpose:
!>   Abort execution with a package-level error message.
!> Notes:
!>   The standalone library intentionally keeps error handling minimal and
!>   avoids the logging infrastructure of the original host codes.
    implicit none
    character(*), intent(in) :: msg

    write(*,*) trim(msg)
    stop 1
end subroutine strux_fail

end module strux_support
