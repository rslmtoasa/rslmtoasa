real(rp) :: gx, gy, gz, gt, a, b, c, amax, bmax, alamda, rmax, gmax, ar2d, sws, vol, vmix
integer :: nq3, nr0, numr, numg, numvr, numvg
!> B7.4: pin the Fermi level to a named region instead of the default free
!> E_F from cluster neutrality (B7 §1.3). Empty = free.
character(len=32) :: fix_fermi_to_region
real(rp), dimension(:), allocatable :: w, wssurf, bsx, bsy, bsz, bkx, bky, bkz, qx3, qy3, qz3, qx, qy, qz, asx, asy, asz, akx, aky, akz, dr, dg, wsimp
real(rp), dimension(:, :), allocatable :: dss, dsz, ds3z2, dsx2y2, dsxy, dzz, dz3z2, am, bm, pm, amad

namelist /charge/ w, gx, gy, gz, gt, wssurf, dss, dsz, ds3z2, dsx2y2, &
   dsxy, dzz, dz3z2, am, bm, pm, bsx, bsy, bsz, bkx, bky, bkz, qx3, qy3, &
   qz3, qx, qy, qz, asx, asy, asz, akx, aky, akz, dr, dg, a, b, c, amax, &
   bmax, alamda, rmax, gmax, ar2d, sws, vol, nq3, nr0, numr, numg, numvr, &
   numvg, amad, wsimp, fix_fermi_to_region
