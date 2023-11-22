integer :: lmax
complex(rp), dimension(9) :: cx0, cx1, wx0, wx1
real(rp), dimension(:), allocatable :: mom
real(rp), dimension(18) :: cshi, dw_l
real(rp), dimension(:, :), allocatable :: pl, center_band, width_band, gravity_center, shifted_band, obar, c, enu, ppar, qpar, srdel, vl
real(rp), dimension(:, :, :), allocatable :: ql
real(rp) :: center_band_s_up, center_band_s_dw, &
    center_band_p_up, center_band_p_dw, &
    center_band_d_up, center_band_d_dw, &
    width_band_s_up, width_band_s_dw, &
    width_band_p_up, width_band_p_dw, &
    width_band_d_up, width_band_d_dw, &
    ws_r, sumec, sumev, etot, utot, ekin, rhoeps, vmad
complex(rp), dimension(9, 2) :: cx, wx, cex, obx
real(rp), dimension(2) :: xi_p, xi_d
namelist /par/ &
    center_band, width_band, shifted_band, obar, gravity_center, &
    sumec, sumev, etot, utot, ekin, rhoeps, vmad, &
    c, enu, ppar, qpar, srdel, vl, &
    pl, mom, ws_r, ql, lmax, &
    center_band_s_up, center_band_s_dw, &
    center_band_p_up, center_band_p_dw, &
    center_band_d_up, center_band_d_dw, &
    width_band_s_up, width_band_s_dw, &
    width_band_p_up, width_band_p_dw, &
    width_band_d_up, width_band_d_dw, &
    cx, wx, cx0, wx0, cx1, wx1,  &
    cshi, dw_l, cex, obx, xi_p, xi_d
