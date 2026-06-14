real(rp) :: zmin, zmax, zstep, wav, vol, rc, r2, celldm, alat
integer :: njij, njijk, reduced_nbas, ntype, ntot, nrec, nmax, nlay, ndim, npe, nclu, nbulk_bulk, nbulk, nbas, kk, dx, dy, dz, dw, n1, n2, n3
character(len=4) :: crystal_sym
character(len=10) :: surftype
character(len=16) :: strux_backend, screening
logical :: pbc, b1, b2, b3, morton_sfc
integer, dimension(:, :), allocatable :: ijpair
real(rp), dimension(:, :), allocatable :: ijktrio
real(rp), dimension(3, 3) :: a
real(rp), dimension(:), allocatable :: z, ct, screening_alpha, screening_sigma
integer, dimension(:), allocatable :: reduced_acr, num, no, izpsurf, izsurf, nosurf, izpo, izp, iz, iu, irec, ib
real(rp), dimension(:, :), allocatable :: primcell, inclu, crsurf, crd, cr, acr
logical :: strux_want_sdot
real(rp) :: strux_solve_scale

namelist /lattice/ a, z, ct, primcell, inclu, crsurf, crd, cr, njij, njijk, &
   acr, zmin, zmax, zstep, wav, vol, rc, r2, celldm, alat, reduced_acr, &
   num, no, izpsurf, izsurf, nosurf, izpo, izp, iz, iu, irec, ib, &
   reduced_nbas, ntype, ntot, nrec, nmax, nlay, ndim, npe, nclu, &
   nbulk_bulk, nbulk, nbas, kk, dx, dy, dz, dw, crystal_sym, surftype, ijpair, ijktrio, &
   b1, b2, b3, pbc, morton_sfc, n1, n2, n3, strux_backend, screening, strux_want_sdot, &
   strux_solve_scale, screening_alpha, screening_sigma
