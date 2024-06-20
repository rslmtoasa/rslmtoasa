real(rp) :: gx, gy, gz, gt, a, b, c, amax, bmax, alamda, rmax, gmax, ar2d, sws, vol, vmix
integer :: nq3, nr0, numr, numg, numvr, numvg
real(rp), dimension(:), allocatable :: w, wssurf, bsx, bsy, bsz, bkx, bky, bkz, qx3, qy3, qz3, qx, qy, qz, asx, asy, asz, akx, aky, akz, dr, dg, wsimp
real(rp), dimension(:, :), allocatable :: dss, dsz, ds3z2, dsx2y2, dsxy, dzz, dz3z2, am, bm, pm, amad

namelist /charge/ w, gx, gy, gz, gt, wssurf, dss, dsz, ds3z2, dsx2y2, &
   dsxy, dzz, dz3z2, am, bm, pm, bsx, bsy, bsz, bkx, bky, bkz, qx3, qy3, &
   qz3, qx, qy, qz, asx, asy, asz, akx, aky, akz, dr, dg, a, b, c, amax, &
   bmax, alamda, rmax, gmax, ar2d, sws, vol, nq3, nr0, numr, numg, numvr, &
   numvg, amad, wsimp
