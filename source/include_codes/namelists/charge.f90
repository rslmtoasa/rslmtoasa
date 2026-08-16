real(rp) :: gx, gy, gz, gt, a, b, c, amax, bmax, alamda, rmax, gmax, ar2d, sws, vol, vmix
integer :: nq3, nr0, numr, numg, numvr, numvg
!> B7.4: pin the Fermi level to a named region instead of the default free
!> E_F from cluster neutrality (B7 §1.3). Empty = free.
character(len=32) :: fix_fermi_to_region
!> B7.5: layered/interface (calctype='L') region widths, in active layers,
!> and optional per-region source E_F for the alignment consistency check
!> (B7 §1.3, CONTRACT_FROZEN_REGION.md §2). fermi_a/fermi_b default to a
!> sentinel (not supplied); see charge%restore_to_default.
integer :: nlay_a, nlay_b
real(rp) :: fermi_a, fermi_b
!> B7.5: compensation weight profile selector (B7 §1.5, G-B7-3). Only
!> 'nef_weighted' (the side-resolved N(E_F) weighting interfacepot already
!> implements) exists today; the field is exposed now so a future profile
!> is a namelist choice, not a new code path.
character(len=32) :: compensation_profile
!> B7 §6 hook: applied bias, target_step = contact_potential + bias. Namelist
!> plumbing only -- inactive until a future task wires it to the PM kernel.
real(rp) :: bias
real(rp), dimension(:), allocatable :: w, wssurf, bsx, bsy, bsz, bkx, bky, bkz, qx3, qy3, qz3, qx, qy, qz, asx, asy, asz, akx, aky, akz, dr, dg, wsimp
real(rp), dimension(:, :), allocatable :: dss, dsz, ds3z2, dsx2y2, dsxy, dzz, dz3z2, am, bm, pm, amad

namelist /charge/ w, gx, gy, gz, gt, wssurf, dss, dsz, ds3z2, dsx2y2, &
   dsxy, dzz, dz3z2, am, bm, pm, bsx, bsy, bsz, bkx, bky, bkz, qx3, qy3, &
   qz3, qx, qy, qz, asx, asy, asz, akx, aky, akz, dr, dg, a, b, c, amax, &
   bmax, alamda, rmax, gmax, ar2d, sws, vol, nq3, nr0, numr, numg, numvr, &
   numvg, amad, wsimp, vmix, fix_fermi_to_region, nlay_a, nlay_b, fermi_a, fermi_b, &
   compensation_profile, bias
