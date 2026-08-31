!------------------------------------------------------------------------------
! RS-LMTO-ASA -- WP7 / gate G7 oracles for the rotating-frame density contract
!------------------------------------------------------------------------------
!
!> Independent tests for `spin_density_mod`. Nothing here calls the production
!> solvers: the model is a hand-built two-site spinor eigensystem, and the two
!> producers are re-implemented from their *definitions* --
!>
!>   k-space:  rho = sum_n f(eps_n) psi_n psi_n^dagger
!>   real space: rho = (i/2pi) int (G(E+i0) - G^dagger(E+i0)) dE,
!>               G(z) = sum_n psi_n psi_n^dagger / (z - eps_n)
!>
!> -- so agreement between them is evidence about the CONVENTION the contract
!> fixes (conjugation order, m_y sign, the 1/pi), not a restatement of one
!> implementation.
!>
!> Covered:
!>   1. Cartesian <-> matrix round trip against an explicit sigma algebra.
!>   2. Equivalent producers: Green-function accumulation vs eigenvector
!>      accumulation fill the same object.
!>   3. Radial projection follows the EXPLICIT axis, not the density’s own
!>      direction and not any stored moment.
!>   4. The two SCF policies are distinct objects.
!>   5. Physicality assertions catch a broken density and pass a good one.
!>   6. q=0 cone invariance: a global spin rotation of density and axis leaves
!>      every projected radial quantity unchanged, and the lab-frame
!>      reconstruction is the identity at q=0.
!------------------------------------------------------------------------------
program test_spin_density
   use precision_mod, only: rp
   use spin_density_mod
   implicit none

   integer :: failures

   failures = 0

   call test_cartesian_roundtrip(failures)
   call test_producer_equivalence(failures)
   call test_explicit_axis(failures)
   call test_policies_distinct(failures)
   call test_physicality(failures)
   call test_q0_cone_invariance(failures)

   if (failures == 0) then
      write (*, '(a)') 'UnitSpinDensity: ALL TESTS PASSED'
      stop 0
   else
      write (*, '(a,i0,a)') 'UnitSpinDensity: ', failures, ' TEST(S) FAILED'
      stop 1
   end if

contains

   subroutine check(label, got, want, tol, failures)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: got, want, tol
      integer, intent(inout) :: failures

      if (abs(got - want) > tol) then
         write (*, '(a,a,a,es16.8,a,es16.8,a,es10.3)') 'FAIL ', label, &
            ': got ', got, ' want ', want, ' tol ', tol
         failures = failures + 1
      end if
   end subroutine check

   subroutine check_true(label, cond, failures)
      character(len=*), intent(in) :: label
      logical, intent(in) :: cond
      integer, intent(inout) :: failures

      if (.not. cond) then
         write (*, '(a,a)') 'FAIL ', label
         failures = failures + 1
      end if
   end subroutine check_true

   !> Reference two-level spinor model: nband spinors on one site, one l channel.
   subroutine model(psi, eps, nband)
      complex(rp), allocatable, intent(out) :: psi(:, :)
      real(rp), allocatable, intent(out) :: eps(:)
      integer, intent(out) :: nband
      real(rp) :: nrm
      integer :: ib

      nband = 3
      allocate (psi(2, nband), eps(nband))
      ! Deliberately generic spinors: complex, non-collinear, unequal weights.
      psi(:, 1) = [cmplx(0.8_rp, 0.1_rp, rp), cmplx(0.3_rp, -0.5_rp, rp)]
      psi(:, 2) = [cmplx(-0.2_rp, 0.6_rp, rp), cmplx(0.7_rp, 0.25_rp, rp)]
      psi(:, 3) = [cmplx(0.45_rp, -0.35_rp, rp), cmplx(-0.15_rp, 0.8_rp, rp)]
      eps = [-0.42_rp, -0.13_rp, 0.31_rp]
      do ib = 1, nband
         nrm = sqrt(real(sum(conjg(psi(:, ib))*psi(:, ib)), rp))
         psi(:, ib) = psi(:, ib)/nrm
      end do
   end subroutine model

   !> Producer A -- eigenvector/occupation accumulation (the k-space definition).
   subroutine fill_from_eigenvectors(sd, psi, eps, ef)
      type(spin_density), intent(inout) :: sd
      complex(rp), intent(in) :: psi(:, :)
      real(rp), intent(in) :: eps(:), ef
      integer :: ib, iorder
      real(rp) :: epow
      complex(rp) :: blk(2, 2)

      do ib = 1, size(eps)
         if (eps(ib) > ef) cycle
         blk(1, 1) = psi(1, ib)*conjg(psi(1, ib))
         blk(2, 2) = psi(2, ib)*conjg(psi(2, ib))
         blk(1, 2) = psi(1, ib)*conjg(psi(2, ib))
         blk(2, 1) = psi(2, ib)*conjg(psi(1, ib))
         epow = 1.0_rp
         do iorder = 1, sd_orders
            call sd%accumulate_block(1, 1, iorder, blk*cmplx(epow, 0.0_rp, rp))
            epow = epow*eps(ib)
         end do
      end do
   end subroutine fill_from_eigenvectors

   !> Producer B -- Green-function accumulation (the real-space definition).
   !> G(z) = sum_n psi_n psi_n^dagger/(z - eps_n) on a retarded contour, then
   !> rho = (i/2pi) int (G - G^dagger) dE by the trapezoidal rule.
   subroutine fill_from_green(sd, psi, eps, ef, eta)
      type(spin_density), intent(inout) :: sd
      complex(rp), intent(in) :: psi(:, :)
      real(rp), intent(in) :: eps(:), ef, eta
      real(rp), parameter :: emin = -6.0_rp
      real(rp), parameter :: pi_l = 3.14159265358979323846_rp
      integer :: ie, ib, iorder, s1, s2, ngrid
      real(rp) :: e, de, w, epow
      complex(rp) :: z, g(2, 2), rho_e(2, 2), acc(2, 2, sd_orders)

      ! Resolve the Lorentzian: de = eta/4 over the occupied window.
      de = 0.25_rp*eta
      ngrid = int((ef - emin)/de) + 1
      acc = (0.0_rp, 0.0_rp)

      do ie = 1, ngrid
         e = emin + de*real(ie - 1, rp)
         if (e > ef) exit
         z = cmplx(e, eta, rp)
         g = (0.0_rp, 0.0_rp)
         do ib = 1, size(eps)
            do s1 = 1, 2
               do s2 = 1, 2
                  g(s1, s2) = g(s1, s2) + psi(s1, ib)*conjg(psi(s2, ib))/(z - eps(ib))
               end do
            end do
         end do
         ! Hermitian spectral density: (i/2pi)(G - G^dagger).
         do s1 = 1, 2
            do s2 = 1, 2
               rho_e(s1, s2) = cmplx(0.0_rp, 1.0_rp, rp)*(g(s1, s2) - conjg(g(s2, s1)))/(2.0_rp*pi_l)
            end do
         end do
         w = de
         if (ie == 1) w = 0.5_rp*de
         epow = 1.0_rp
         do iorder = 1, sd_orders
            acc(:, :, iorder) = acc(:, :, iorder) + rho_e(:, :)*cmplx(w*epow, 0.0_rp, rp)
            epow = epow*e
         end do
      end do

      do iorder = 1, sd_orders
         call sd%accumulate_block(1, 1, iorder, acc(:, :, iorder))
      end do
   end subroutine fill_from_green

   !---------------------------------------------------------------------------
   subroutine test_cartesian_roundtrip(failures)
      integer, intent(inout) :: failures
      complex(rp) :: blk(2, 2)
      real(rp) :: n, m(3), n2, m2(3)

      n = 1.7_rp
      m = [0.3_rp, -0.45_rp, 0.9_rp]
      call sd_matrix_from_cartesian(n, m, blk)

      ! Independent sigma algebra: rho = (n*1 + m.sigma)/2.
      call check('roundtrip rho_uu', real(blk(1, 1), rp), 0.5_rp*(n + m(3)), 1.0e-14_rp, failures)
      call check('roundtrip rho_dd', real(blk(2, 2), rp), 0.5_rp*(n - m(3)), 1.0e-14_rp, failures)
      call check('roundtrip Re rho_ud', real(blk(1, 2), rp), 0.5_rp*m(1), 1.0e-14_rp, failures)
      call check('roundtrip Im rho_ud', aimag(blk(1, 2)), -0.5_rp*m(2), 1.0e-14_rp, failures)

      call sd_cartesian_from_matrix(blk, n2, m2)
      call check('roundtrip n', n2, n, 1.0e-14_rp, failures)
      call check('roundtrip mx', m2(1), m(1), 1.0e-14_rp, failures)
      call check('roundtrip my', m2(2), m(2), 1.0e-14_rp, failures)
      call check('roundtrip mz', m2(3), m(3), 1.0e-14_rp, failures)
   end subroutine test_cartesian_roundtrip

   !---------------------------------------------------------------------------
   subroutine test_producer_equivalence(failures)
      integer, intent(inout) :: failures
      type(spin_density) :: sd_k, sd_g, sd_g2
      complex(rp), allocatable :: psi(:, :)
      real(rp), allocatable :: eps(:)
      integer :: nband
      real(rp) :: ef, dev, dev2
      real(rp) :: nk, mk(3), ng, mg(3)

      call model(psi, eps, nband)
      ef = 0.0_rp   ! two of three levels occupied

      sd_k = spin_density(1, 1)
      sd_k%producer = sd_producer_kspace
      call fill_from_eigenvectors(sd_k, psi, eps, ef)

      sd_g = spin_density(1, 1)
      sd_g%producer = sd_producer_rs
      call fill_from_green(sd_g, psi, eps, ef, 4.0e-4_rp)

      sd_g2 = spin_density(1, 1)
      sd_g2%producer = sd_producer_rs
      call fill_from_green(sd_g2, psi, eps, ef, 1.0e-4_rp)

      ! The Green-function route carries a finite imaginary part, so it can only
      ! agree with the sharp eigenvector sum up to the Lorentzian tail that
      ! leaks past E_F -- an O(eta) error, NOT a convention error. Assert both
      ! the size AND the scaling: quartering eta must quarter the residual.
      ! Anything that survives that limit (a wrong conjugation, a flipped m_y,
      ! a missing 1/pi) would show up as an eta-independent offset.
      dev = sd_k%max_deviation(sd_g)
      dev2 = sd_k%max_deviation(sd_g2)
      call check_true('producer residual is small (eta=4e-4)', dev < 2.0e-3_rp, failures)
      call check_true('producer residual scales with eta', dev2 < 0.35_rp*dev, failures)
      if (dev2 >= 0.35_rp*dev) then
         write (*, '(a,es12.4,a,es12.4)') '  dev(4e-4) = ', dev, '  dev(1e-4) = ', dev2
      end if

      call sd_k%cartesian(1, 1, 1, nk, mk)
      call sd_g2%cartesian(1, 1, 1, ng, mg)
      call check('producer equivalence n', ng, nk, 1.0e-3_rp, failures)
      call check('producer equivalence mx', mg(1), mk(1), 1.0e-3_rp, failures)
      call check('producer equivalence my', mg(2), mk(2), 1.0e-3_rp, failures)
      call check('producer equivalence mz', mg(3), mk(3), 1.0e-3_rp, failures)

      ! The k-space producer must reproduce the pre-WP7 spinor signs exactly:
      ! m_x = 2 Re(conjg(u) d), m_y = 2 Im(conjg(u) d), m_z = |u|^2 - |d|^2.
      block
         real(rp) :: mx_ref, my_ref, mz_ref, n_ref
         integer :: ib
         mx_ref = 0.0_rp; my_ref = 0.0_rp; mz_ref = 0.0_rp; n_ref = 0.0_rp
         do ib = 1, nband
            if (eps(ib) > ef) cycle
            mx_ref = mx_ref + 2.0_rp*real(conjg(psi(1, ib))*psi(2, ib), rp)
            my_ref = my_ref + 2.0_rp*aimag(conjg(psi(1, ib))*psi(2, ib))
            mz_ref = mz_ref + real(conjg(psi(1, ib))*psi(1, ib) - conjg(psi(2, ib))*psi(2, ib), rp)
            n_ref = n_ref + real(conjg(psi(1, ib))*psi(1, ib) + conjg(psi(2, ib))*psi(2, ib), rp)
         end do
         call check('spinor convention n', nk, n_ref, 1.0e-13_rp, failures)
         call check('spinor convention mx', mk(1), mx_ref, 1.0e-13_rp, failures)
         call check('spinor convention my', mk(2), my_ref, 1.0e-13_rp, failures)
         call check('spinor convention mz', mk(3), mz_ref, 1.0e-13_rp, failures)
      end block

      deallocate (psi, eps)
   end subroutine test_producer_equivalence

   !---------------------------------------------------------------------------
   subroutine test_explicit_axis(failures)
      integer, intent(inout) :: failures
      type(spin_density) :: sd
      complex(rp) :: blk(2, 2)
      real(rp) :: n, m(3), up, dn, up2, dn2

      n = 2.0_rp
      m = [0.0_rp, 0.0_rp, 0.6_rp]
      call sd_matrix_from_cartesian(n, m, blk)

      sd = spin_density(1, 1)
      call sd%accumulate_block(1, 1, 1, blk)

      call sd%set_axis(1, [0.0_rp, 0.0_rp, 1.0_rp])
      call sd%project_radial(1, 1, 1, up, dn)
      call check('axis z: up', up, 1.3_rp, 1.0e-14_rp, failures)
      call check('axis z: dn', dn, 0.7_rp, 1.0e-14_rp, failures)

      ! Same density, different STATED axis -> different radial channels. This
      ! is the property the pre-WP7 code could not express: the axis is an
      ! input to the projection, not a property of the density.
      call sd%set_axis(1, [1.0_rp, 0.0_rp, 0.0_rp])
      call sd%project_radial(1, 1, 1, up2, dn2)
      call check('axis x: up', up2, 1.0_rp, 1.0e-14_rp, failures)
      call check('axis x: dn', dn2, 1.0_rp, 1.0e-14_rp, failures)
      call check_true('stated axis changes the projection', abs(up - up2) > 0.1_rp, failures)

      ! Charge is axis-independent.
      call check('charge is axis independent', up + dn, up2 + dn2, 1.0e-14_rp, failures)
   end subroutine test_explicit_axis

   !---------------------------------------------------------------------------
   subroutine test_policies_distinct(failures)
      integer, intent(inout) :: failures
      type(spin_density) :: sd
      complex(rp) :: blk(2, 2)
      real(rp) :: n, m(3)
      real(rp) :: axis_c(3), axis_r(3), tr_c(3), tr_r(3), tq_c(3), tq_r(3)
      real(rp) :: q_c, q_r, ml_c, ml_r
      real(rp) :: reference(3)

      ! A moment that is deliberately NOT along the imposed reference.
      n = 3.0_rp
      m = [0.5_rp, 0.0_rp, 1.2_rp]
      call sd_matrix_from_cartesian(n, m, blk)

      sd = spin_density(1, 1)
      call sd%accumulate_block(1, 1, 1, blk)
      reference = [0.0_rp, 0.0_rp, 1.0_rp]

      sd%policy = sd_constrained_spiral
      call sd%resolve_site_axis(1, reference, axis_c, q_c, ml_c, tr_c, tq_c)

      sd%policy = sd_relaxed_reference
      call sd%resolve_site_axis(1, reference, axis_r, q_r, ml_r, tr_r, tq_r)

      ! Constrained: axis pinned to the reference; transverse density survives
      ! as a residual with a nonzero torque.
      call check('constrained axis z', axis_c(3), 1.0_rp, 1.0e-14_rp, failures)
      call check('constrained m_long', ml_c, 1.2_rp, 1.0e-14_rp, failures)
      call check('constrained transverse x', tr_c(1), 0.5_rp, 1.0e-14_rp, failures)
      call check_true('constrained torque is nonzero', sqrt(sum(tq_c**2)) > 0.1_rp, failures)

      ! Relaxed: axis follows the full Cartesian moment; nothing left over.
      call check('relaxed m_long is |m|', ml_r, sqrt(sum(m**2)), 1.0e-14_rp, failures)
      call check('relaxed axis x', axis_r(1), m(1)/sqrt(sum(m**2)), 1.0e-14_rp, failures)
      call check('relaxed transverse vanishes', sqrt(sum(tr_r**2)), 0.0_rp, 1.0e-14_rp, failures)
      call check('relaxed torque vanishes', sqrt(sum(tq_r**2)), 0.0_rp, 1.0e-14_rp, failures)

      ! The two policies are genuinely different reductions of one density.
      call check_true('policies give different axes', &
                      maxval(abs(axis_c - axis_r)) > 1.0e-3_rp, failures)
      call check_true('policies give different longitudinal moments', &
                      abs(ml_c - ml_r) > 1.0e-3_rp, failures)
      call check('policies agree on charge', q_c, q_r, 1.0e-14_rp, failures)

      ! A producer re-accumulating into a configured object must NOT silently
      ! reset the policy back to the default -- that would make the relaxed
      ! policy unreachable while still appearing to be selected.
      sd%policy = sd_relaxed_reference
      call sd%zero_density()
      call check_true('zero_density keeps the configured policy', &
                      trim(sd%policy) == trim(sd_relaxed_reference), failures)
      call check('zero_density clears the density', &
                 real(sd%rho(1, 1, 1, 1, 1), rp), 0.0_rp, 1.0e-14_rp, failures)
      call check_true('zero_density forgets the stated axis', &
                      .not. sd%axis_set(1), failures)
   end subroutine test_policies_distinct

   !---------------------------------------------------------------------------
   subroutine test_physicality(failures)
      integer, intent(inout) :: failures
      type(spin_density) :: sd
      complex(rp) :: blk(2, 2)
      logical :: ok
      character(len=256) :: msg

      ! Physical: |m| < n.
      call sd_matrix_from_cartesian(2.0_rp, [0.3_rp, 0.2_rp, 1.0_rp], blk)
      sd = spin_density(1, 1)
      call sd%accumulate_block(1, 1, 1, blk)
      call sd%check_physicality(1.0e-10_rp, ok, msg)
      call check_true('physical density passes', ok, failures)
      call check('electron count', sd%electron_count(), 2.0_rp, 1.0e-14_rp, failures)
      call sd%check_physicality(1.0e-10_rp, ok, msg, expected_electrons=2.0_rp)
      call check_true('electron count check passes', ok, failures)
      call sd%check_physicality(1.0e-10_rp, ok, msg, expected_electrons=3.0_rp)
      call check_true('wrong electron count is caught', .not. ok, failures)

      ! Unphysical: |m| > n (negative eigenvalue).
      call sd_matrix_from_cartesian(1.0_rp, [0.0_rp, 0.0_rp, 1.5_rp], blk)
      sd = spin_density(1, 1)
      call sd%accumulate_block(1, 1, 1, blk)
      call sd%check_physicality(1.0e-10_rp, ok, msg)
      call check_true('|m| > n is caught', .not. ok, failures)

      ! Unphysical: non-Hermitian.
      sd = spin_density(1, 1)
      blk(1, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      blk(2, 2) = cmplx(1.0_rp, 0.0_rp, rp)
      blk(1, 2) = cmplx(0.2_rp, 0.1_rp, rp)
      blk(2, 1) = cmplx(0.2_rp, 0.1_rp, rp)   ! not conjg(blk(1,2))
      call sd%accumulate_block(1, 1, 1, blk)
      call sd%check_physicality(1.0e-10_rp, ok, msg)
      call check_true('non-Hermitian density is caught', .not. ok, failures)

      ! Unphysical: negative trace.
      sd = spin_density(1, 1)
      call sd_matrix_from_cartesian(-0.5_rp, [0.0_rp, 0.0_rp, 0.0_rp], blk)
      call sd%accumulate_block(1, 1, 1, blk)
      call sd%check_physicality(1.0e-10_rp, ok, msg)
      call check_true('negative trace is caught', .not. ok, failures)
   end subroutine test_physicality

   !---------------------------------------------------------------------------
   !> q = 0 cone invariance. Without SOC the whole problem is invariant under a
   !> GLOBAL spin rotation: rotating the density and its stated axis together
   !> must leave charge, longitudinal moment and both radial channels exactly
   !> unchanged. At q = 0 the lab-frame reconstruction is additionally the
   !> identity (phase = q.R = 0), so lab and rotating frames coincide.
   subroutine test_q0_cone_invariance(failures)
      integer, intent(inout) :: failures
      type(spin_density) :: sd, sd_rot
      complex(rp), allocatable :: psi(:, :), psi_rot(:, :)
      real(rp), allocatable :: eps(:)
      integer :: nband, ib
      real(rp) :: ef, theta, phi
      complex(rp) :: u(2, 2)
      real(rp) :: axis0(3), axis1(3)
      real(rp) :: up0, dn0, up1, dn1
      real(rp) :: n0, m0(3), n1, m1(3), nlab, mlab(3)
      real(rp) :: q0a, c0, q2a, q0b, c1, q2b

      call model(psi, eps, nband)
      ef = 0.0_rp

      ! Global SU(2) rotation: cone angle theta about y then phi about z.
      theta = 0.7_rp
      phi = 1.1_rp
      u(1, 1) = cos(0.5_rp*theta)*cmplx(cos(0.5_rp*phi), -sin(0.5_rp*phi), rp)
      u(1, 2) = -sin(0.5_rp*theta)*cmplx(cos(0.5_rp*phi), -sin(0.5_rp*phi), rp)
      u(2, 1) = sin(0.5_rp*theta)*cmplx(cos(0.5_rp*phi), sin(0.5_rp*phi), rp)
      u(2, 2) = cos(0.5_rp*theta)*cmplx(cos(0.5_rp*phi), sin(0.5_rp*phi), rp)

      allocate (psi_rot(2, nband))
      do ib = 1, nband
         psi_rot(:, ib) = matmul(u, psi(:, ib))
      end do

      sd = spin_density(1, 1)
      call fill_from_eigenvectors(sd, psi, eps, ef)
      sd_rot = spin_density(1, 1)
      call fill_from_eigenvectors(sd_rot, psi_rot, eps, ef)

      call sd%cartesian(1, 1, 1, n0, m0)
      call sd_rot%cartesian(1, 1, 1, n1, m1)

      ! Charge and moment magnitude are rotation invariant; the direction is not.
      call check('cone invariance: charge', n1, n0, 1.0e-13_rp, failures)
      call check('cone invariance: |m|', sqrt(sum(m1**2)), sqrt(sum(m0**2)), 1.0e-13_rp, failures)
      call check_true('the rotation actually moved the moment', &
                      maxval(abs(m1 - m0)) > 1.0e-3_rp, failures)

      ! Rotate the stated axis by the SAME rotation -> identical radial channels.
      axis0 = m0/sqrt(sum(m0**2))
      axis1 = m1/sqrt(sum(m1**2))
      call sd%set_axis(1, axis0)
      call sd_rot%set_axis(1, axis1)
      call sd%project_radial(1, 1, 1, up0, dn0)
      call sd_rot%project_radial(1, 1, 1, up1, dn1)
      call check('cone invariance: up channel', up1, up0, 1.0e-13_rp, failures)
      call check('cone invariance: down channel', dn1, dn0, 1.0e-13_rp, failures)

      ! Every band moment the SCF consumes is invariant too.
      call sd%radial_band_moments(1, 1, 1, q0a, c0, q2a)
      call sd_rot%radial_band_moments(1, 1, 1, q0b, c1, q2b)
      call check('cone invariance: q0', q0b, q0a, 1.0e-13_rp, failures)
      call check('cone invariance: band centre', c1, c0, 1.0e-13_rp, failures)
      call check('cone invariance: q2', q2b, q2a, 1.0e-13_rp, failures)

      ! q = 0 -> the lab-frame reconstruction is the identity.
      call sd_rot%lab_frame_moment(1, 1, 1, 0.0_rp, nlab, mlab)
      call check('q=0 lab frame is identity (n)', nlab, n1, 1.0e-14_rp, failures)
      call check('q=0 lab frame is identity (mx)', mlab(1), m1(1), 1.0e-14_rp, failures)
      call check('q=0 lab frame is identity (my)', mlab(2), m1(2), 1.0e-14_rp, failures)
      call check('q=0 lab frame is identity (mz)', mlab(3), m1(3), 1.0e-14_rp, failures)

      ! A nonzero phase rotates the transverse moment about z and preserves both
      ! the longitudinal component and |m| -- the single-q lab reconstruction.
      call sd_rot%lab_frame_moment(1, 1, 1, 0.9_rp, nlab, mlab)
      call check('lab frame preserves mz', mlab(3), m1(3), 1.0e-14_rp, failures)
      call check('lab frame preserves |m|', sqrt(sum(mlab**2)), sqrt(sum(m1**2)), 1.0e-13_rp, failures)
      call check('lab frame rotates mx', mlab(1), &
                 cos(0.9_rp)*m1(1) - sin(0.9_rp)*m1(2), 1.0e-13_rp, failures)

      deallocate (psi, psi_rot, eps)
   end subroutine test_q0_cone_invariance

end program test_spin_density
