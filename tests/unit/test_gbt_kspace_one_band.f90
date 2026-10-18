! UnitGbtKspaceOneBand
! Relation R1 (generalized Bloch theorem): Herring 1966; Sandratskii, Adv. Phys. 47, 91 (1998); oracle O1 in
! tests/gbt/oracles (closed form for one site; two-site cell checked by the oracle against an explicit supercell).
! Bug caught: wrong k -+ q/2 assignment to the spin channels, a dropped or misplaced U_a^dag ... U_b frame link,
! or W and C placed in the wrong frame.
! The test passes the oracle's model to gbt_kspace_h1 and compares eigenvalues; it has no physics of its own.
! usage: test_gbt_kspace_one_band <oracle directory>
program test_gbt_kspace_one_band
   use precision_mod, only: rp
   use gbt_kspace_mod, only: gbt_kspace_h1
   use gbt_structure_mod, only: gbt_reference_frame
   use math_mod, only: inverse_3x3
   implicit none

   integer, parameter :: ncase_max = 72, nbond_max = 12
   real(rp), parameter :: pi = acos(-1.0_rp)
   character(len=512) :: dir, s_rule, frame_rule, eigenvalue_order, k_rule, eig_shape, cart_rule
   character(len=512) :: site_rule_parallel, site_rule_antiparallel
   integer :: nk, nq, ntheta, nphi, nsite, nbond, bond_site_a(nbond_max), bond_site_b(nbond_max)
   real(rp) :: t, delta, w_up, w_dn, c_majority, c_minority, bond_s, bond_d(3*nbond_max)
   real(rp) :: k_direct(9), q_direct(9), theta_deg(4), phi_deg(2), a_cell(9), tau(6), k_cart(9), q_cart(9)
   real(rp) :: eig(2*ncase_max), eig_parallel(4*ncase_max), eig_antiparallel(4*ncase_max)
   real(rp) :: eig_identity_ry, omega_identity_ry, n_identity, model_eig_ry, model_ft_rel
   namelist /o1_model/ t, delta, w_up, w_dn, c_majority, c_minority, s_rule, frame_rule, eigenvalue_order
   namelist /o1_cases/ nk, nq, ntheta, nphi, k_rule, k_direct, q_direct, theta_deg, phi_deg
   namelist /o1_one_site/ eig_shape, nsite, a_cell, tau, cart_rule, k_cart, q_cart, eig
   namelist /o1_two_site/ eig_shape, nsite, a_cell, tau, site_rule_parallel, site_rule_antiparallel, nbond, bond_s, &
      bond_site_a, bond_site_b, bond_d, cart_rule, k_cart, q_cart, eig_parallel, eig_antiparallel
   namelist /gbt_tolerances/ eig_identity_ry, omega_identity_ry, n_identity, model_eig_ry, model_ft_rel
   real(rp) :: one_site_d(3, 6), two_site_d(3, nbond_max), dev_one, dev_par, dev_anti
   integer :: u, one_site_a(6), one_site_b(6), i

   call get_command_argument(1, dir)
   if (len_trim(dir) == 0) error stop 'usage: test_gbt_kspace_one_band <oracle directory>'
   model_eig_ry = -1.0_rp
   open (newunit=u, file=trim(dir)//'/o1_one_band.nml', status='old', action='read')
   read (u, nml=o1_model)
   read (u, nml=o1_cases)
   read (u, nml=o1_one_site)
   if (nk*nq*ntheta*nphi /= ncase_max) error stop 'o1_one_band.nml: unexpected case count'
   one_site_a = 1
   one_site_b = 1
   one_site_d = reshape([1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1]*1.0_rp, [3, 6])   ! simple cubic bonds, a = 1
   dev_one = deviation(1, 6, one_site_a, one_site_b, one_site_d, t, eig, 2, .false.)
   read (u, nml=o1_two_site)
   close (u)
   do i = 1, nbond
      two_site_d(:, i) = matmul(inverse_3x3(reshape(a_cell, [3, 3])), bond_d(3*i - 2:3*i))   ! Cartesian -> direct
   end do
   dev_par = deviation(2, nbond, bond_site_a, bond_site_b, two_site_d, bond_s, eig_parallel, 4, .false.)
   dev_anti = deviation(2, nbond, bond_site_a, bond_site_b, two_site_d, bond_s, eig_antiparallel, 4, .true.)
   open (newunit=u, file=trim(dir)//'/tolerances.nml', status='old', action='read')
   read (u, nml=gbt_tolerances)
   close (u)
   if (.not. model_eig_ry > 0.0_rp) error stop 'tolerances.nml: model_eig_ry missing or <= 0'

   write (*, '(a,3(es12.4,1x),a,es10.2)') 'O1 max |eigenvalue - oracle| (one site, two sites parallel, antiparallel): ', &
      dev_one, dev_par, dev_anti, ' tolerance ', model_eig_ry
   if (max(dev_one, dev_par, dev_anti) >= model_eig_ry) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   !> Largest |eigenvalue - oracle| over all (k, q, theta, phi) of one cell type. Direct coordinates as in the file;
   !> the antiparallel second site has theta = pi - theta_ss and phi = phi_ss + pi (the oracle's site_rule_antiparallel).
   real(rp) function deviation(n_site, n_bond, bond_a, bond_b, bond_dir, s_value, reference, n_eig, antiparallel)
      integer, intent(in) :: n_site, n_bond, bond_a(:), bond_b(:), n_eig
      real(rp), intent(in) :: bond_dir(:, :), s_value, reference(:)
      logical, intent(in) :: antiparallel

      complex(rp) :: sblock(1, 1, n_bond), frame(2, 2, n_site), w(1, 2, n_site), c(1, 2, n_site)
      complex(rp) :: h(2*n_site, 2*n_site), work(64)
      real(rp) :: theta, phi, ev(2*n_site), rwork(6*n_site)
      integer :: ik, iq, it, ip, a, info, icase

      sblock = cmplx(s_value, 0.0_rp, rp)
      w(1, 1, :) = cmplx(w_up, 0.0_rp, rp)
      w(1, 2, :) = cmplx(w_dn, 0.0_rp, rp)
      c(1, 1, :) = cmplx(c_majority, 0.0_rp, rp)
      c(1, 2, :) = cmplx(c_minority, 0.0_rp, rp)
      deviation = 0.0_rp
      do ip = 1, nphi
         do it = 1, ntheta
            do a = 1, n_site
               theta = theta_deg(it)*pi/180.0_rp
               phi = phi_deg(ip)*pi/180.0_rp
               if (antiparallel .and. a == 2) then
                  theta = pi - theta
                  phi = phi + pi
               end if
               call gbt_reference_frame(theta, phi, frame(:, :, a))
            end do
            do iq = 1, nq
               do ik = 1, nk
                  call gbt_kspace_h1(1, n_site, n_bond, bond_a(1:n_bond), bond_b(1:n_bond), bond_dir(:, 1:n_bond), sblock, &
                                     frame, w, c, k_direct(3*ik - 2:3*ik), q_direct(3*iq - 2:3*iq), h)
                  call zheev('N', 'U', 2*n_site, h, 2*n_site, ev, work, 64, rwork, info)
                  if (info /= 0) error stop 'zheev failed'
                  icase = ik + nk*(iq - 1 + nq*(it - 1 + ntheta*(ip - 1)))
                  deviation = max(deviation, maxval(abs(ev - reference(n_eig*(icase - 1) + 1:n_eig*icase))))
               end do
            end do
         end do
      end do
   end function deviation

end program test_gbt_kspace_one_band
