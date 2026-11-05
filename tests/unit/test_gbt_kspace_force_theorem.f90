! UnitGbtKspaceForceTheorem
! Relations R7 (force theorem), R8 (frozen magnon = LKAG), R11 (Mermin grand potential at fixed mu).
! References: Mermin, Phys. Rev. 137, A1441 (1965); Pajda et al., PRB 64, 174402 (2001), the J(q) expression accompanying
! Eq. (7); oracle O2 in tests/gbt/oracles (one-site model of O1, three sets, mesh, k_B T and mu fixed by the oracle).
! Bug caught: E_F*DeltaN mixed into the energy (band energy at fixed mu instead of Omega, finding F2); wrong finite-T Omega.
! The test passes the oracle's model to gbt_kspace_h1, diagonalises, and evaluates Omega and N at the oracle's mu; it has
! no physics of its own. Omega(0, theta) = Omega(0, 0) is a global rotation (R6): by construction, not counted.
! usage: test_gbt_kspace_force_theorem <oracle directory>
program test_gbt_kspace_force_theorem
   use precision_mod, only: rp
   use gbt_kspace_mod, only: gbt_kspace_h1, grand_potential
   use gbt_structure_mod, only: gbt_reference_frame
   use reciprocal_mod, only: fermi_dirac_occupation
   implicit none

   integer, parameter :: nq_max = 5, nth_max = 7, nmesh = 24
   real(rp), parameter :: pi = acos(-1.0_rp)
   character(len=1024) :: dir, s_rule, frame_rule, eigenvalue_order, mesh_rule, richardson_rule, array_rule, lkag_rule
   character(len=1024) :: mutation_rule, set_rule
   integer :: nk_mesh(3), nq, ntheta, n_richardson
   real(rp) :: t, delta, w_up, w_dn, c_majority, c_minority, k_offset(3), kt_ry, q_direct(3*nq_max), q_cart(3*nq_max)
   real(rp) :: theta_deg(nth_max), n_electrons, c_shift_ry, mu_ry, omega_00
   real(rp), dimension(nq_max*nth_max) :: omega, n, eband_fixed_mu
   real(rp), dimension(nth_max) :: omega_ref, n_ref, eband_fixed_mu_ref
   real(rp) :: omega_q0(nq_max), j0_minus_jq(nq_max)
   real(rp) :: eig_identity_ry, omega_identity_ry, n_identity, model_eig_ry, model_ft_rel, scf_identity_ry, scf_moment_mub
   namelist /o1_model/ t, delta, w_up, w_dn, c_majority, c_minority, s_rule, frame_rule, eigenvalue_order
   namelist /o2_settings/ nk_mesh, k_offset, kt_ry, mesh_rule, nq, ntheta, n_richardson, q_direct, q_cart, theta_deg, &
      richardson_rule, array_rule, lkag_rule, mutation_rule, set_rule
   namelist /o2_set1/ n_electrons, c_shift_ry, mu_ry, omega, n, eband_fixed_mu, omega_ref, n_ref, eband_fixed_mu_ref, &
      omega_q0, omega_00, j0_minus_jq
   namelist /o2_set2/ n_electrons, c_shift_ry, mu_ry, omega, n, eband_fixed_mu, omega_ref, n_ref, eband_fixed_mu_ref, &
      omega_q0, omega_00, j0_minus_jq
   namelist /o2_set3/ n_electrons, c_shift_ry, mu_ry, omega, n, eband_fixed_mu, omega_ref, n_ref, eband_fixed_mu_ref, &
      omega_q0, omega_00, j0_minus_jq
   namelist /gbt_tolerances/ eig_identity_ry, omega_identity_ry, n_identity, model_eig_ry, model_ft_rel, scf_identity_ry, scf_moment_mub
   real(rp) :: dev_omega, dev_n, dev_ft, om(nq_max, nth_max), nn(nq_max, nth_max), om_ref(nth_max), nn_ref(nth_max)
   integer :: u, iset
   logical :: failed

   call get_command_argument(1, dir)
   if (len_trim(dir) == 0) error stop 'usage: test_gbt_kspace_force_theorem <oracle directory>'
   open (newunit=u, file=trim(dir)//'/tolerances.nml', status='old', action='read')
   read (u, nml=gbt_tolerances)
   close (u)
   open (newunit=u, file=trim(dir)//'/o1_one_band.nml', status='old', action='read')
   read (u, nml=o1_model)
   close (u)
   open (newunit=u, file=trim(dir)//'/o2_force_theorem.nml', status='old', action='read')
   read (u, nml=o2_settings)
   if (nq /= nq_max .or. ntheta /= nth_max .or. n_richardson /= 3 .or. any(nk_mesh /= nmesh)) error stop 'o2: unexpected shape'

   failed = .false.
   do iset = 1, 3
      select case (iset)
      case (1)
         read (u, nml=o2_set1)
      case (2)
         read (u, nml=o2_set2)
      case (3)
         read (u, nml=o2_set3)
      end select
      call probe_all(om, nn, om_ref, nn_ref)
      dev_omega = maxval(abs(om - reshape(omega, [nq_max, nth_max])))
      dev_n = maxval(abs(nn - reshape(n, [nq_max, nth_max])))
      dev_ft = ft_deviation(om, om_ref)
      write (*, '(a,i1,2a,3(es10.3,1x),a,3(es8.1,1x))') 'set ', iset, ' max|Omega-oracle| [Ry], max|N-oracle|, ', &
         'FT residual/max|J0-Jq|: ', dev_omega, dev_n, dev_ft, ' tol ', omega_identity_ry, n_identity, model_ft_rel
      write (*, '(a,es10.2)') '   Omega(0,theta) - Omega(0,0), by construction (R6), not counted: ', maxval(abs(om_ref - om_ref(1)))
      if (dev_omega >= omega_identity_ry .or. dev_n >= n_identity .or. dev_ft >= model_ft_rel) failed = .true.
   end do
   close (u)
   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   !> Omega and N at the oracle's mu for every (q, theta) of the current set; om_ref and nn_ref are the q = 0 column.
   subroutine probe_all(om_out, nn_out, om0_out, nn0_out)
      real(rp), intent(out) :: om_out(nq_max, nth_max), nn_out(nq_max, nth_max), om0_out(nth_max), nn0_out(nth_max)
      integer :: iq, it

      do it = 1, nth_max
         do iq = 1, nq_max
            call probe(q_direct(3*iq - 2:3*iq), theta_deg(it), om_out(iq, it), nn_out(iq, it))
         end do
         call probe([0.0_rp, 0.0_rp, 0.0_rp], theta_deg(it), om0_out(it), nn0_out(it))
      end do
   end subroutine probe_all

   !> One-site cone: Omega and N over the Gamma-centred mesh for one (q, theta); simple cubic bonds, a = 1.
   subroutine probe(q, theta_in, omega_out, n_out)
      real(rp), intent(in) :: q(3), theta_in
      real(rp), intent(out) :: omega_out, n_out

      real(rp), parameter :: bond_d(3, 6) = reshape([1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1]*1.0_rp, [3, 6])
      integer, parameter :: bond_a(6) = 1, bond_b(6) = 1
      complex(rp) :: sblock(1, 1, 6), frame(2, 2, 1), w(1, 2, 1), c(1, 2, 1), h(2, 2), work(64)
      real(rp) :: ev(2, nmesh**3), wk(nmesh**3), k(3), rwork(6)
      integer :: i, j, l, ik, info

      sblock = cmplx(t, 0.0_rp, rp)
      w = cmplx(1.0_rp, 0.0_rp, rp)
      w(1, 1, 1) = cmplx(w_up, 0.0_rp, rp)
      w(1, 2, 1) = cmplx(w_dn, 0.0_rp, rp)
      c(1, 1, 1) = cmplx(c_majority + c_shift_ry, 0.0_rp, rp)
      c(1, 2, 1) = cmplx(c_minority + c_shift_ry, 0.0_rp, rp)
      call gbt_reference_frame(theta_in*pi/180.0_rp, 0.0_rp, frame(:, :, 1))
      ik = 0
      do l = 0, nmesh - 1
         do j = 0, nmesh - 1
            do i = 0, nmesh - 1
               ik = ik + 1
               k = k_offset + real([i, j, l], rp)/nmesh
               call gbt_kspace_h1(1, 1, 6, bond_a, bond_b, bond_d, sblock, frame, w, c, k, q, h)
               call zheev('N', 'U', 2, h, 2, ev(:, ik), work, 64, rwork, info)
               if (info /= 0) error stop 'zheev failed'
            end do
         end do
      end do
      wk = 1.0_rp/nmesh**3
      omega_out = grand_potential(2, nmesh**3, ev, wk, mu_ry, kt_ry)
      n_out = 0.0_rp
      do ik = 1, nmesh**3
         n_out = n_out + wk(ik)*(fermi_dirac_occupation(ev(1, ik), mu_ry, kt_ry) + fermi_dirac_occupation(ev(2, ik), mu_ry, kt_ry))
      end do
   end subroutine probe

   !> Richardson (Neville, polynomial in sin^2 theta) extrapolation to theta = 0 of [Omega(q,theta) - Omega(0,theta)]/sin^2
   !> over the first three angles, against J(0) - J(q); residual normalised by max_q |J(0) - J(q)| (tolerances.nml).
   real(rp) function ft_deviation(om_q, om_0)
      real(rp), intent(in) :: om_q(nq_max, nth_max), om_0(nth_max)

      real(rp) :: x(3), f(3), extrap
      integer :: iq, i, m

      ft_deviation = 0.0_rp
      do iq = 1, nq_max
         do i = 1, 3
            x(i) = sin(theta_deg(i)*pi/180.0_rp)**2
            f(i) = (om_q(iq, i) - om_0(i))/x(i)
         end do
         do m = 1, 2
            f(1:3 - m) = (x(1 + m:3)*f(1:3 - m) - x(1:3 - m)*f(2:4 - m))/(x(1 + m:3) - x(1:3 - m))
         end do
         extrap = f(1)
         ft_deviation = max(ft_deviation, abs(extrap - j0_minus_jq(iq)))
      end do
      ft_deviation = ft_deviation/maxval(abs(j0_minus_jq))
   end function ft_deviation

end program test_gbt_kspace_force_theorem
