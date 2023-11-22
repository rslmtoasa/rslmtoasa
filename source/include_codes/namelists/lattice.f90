real(rp) :: zmin, zmax, zstep, wav, vol, rc, r2, celldm, alat
integer :: njij, njijk, reduced_nbas, ntype, ntot, nrec, nmax, nlay, ndim, npe, nclu, nbulk_bulk, nbulk, nbas, kk, dx, dy, dz, dw
character(len=4) :: crystal_sym
character(len=10) :: surftype
integer, dimension(:,:), allocatable :: ijpair
real(rp), dimension(:,:), allocatable :: ijktrio
real(rp), dimension(3, 3) :: a
real(rp), dimension(:), allocatable :: z, ct
integer, dimension(:), allocatable :: reduced_acr, num, no, izpsurf, izsurf, nosurf, izpo, izp, iz, iu, irec, ib
real(rp), dimension(:, :), allocatable :: primcell, inclu, crsurf, crd, cr, acr

namelist /lattice/ a, z, ct, primcell, inclu, crsurf, crd, cr, njij, njijk, &
    acr, zmin, zmax, zstep, wav, vol, rc, r2, celldm, alat, reduced_acr, &
    num, no, izpsurf, izsurf, nosurf, izpo, izp, iz, iu, irec, ib, &
    reduced_nbas, ntype, ntot, nrec, nmax, nlay, ndim, npe, nclu, &
    nbulk_bulk, nbulk, nbas, kk, dx, dy, dz, dw, crystal_sym, surftype, ijpair, ijktrio
