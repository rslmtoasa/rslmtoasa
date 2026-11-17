!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
! PROGRAM: test_vacuum_lead
!
!> @brief Validation-ladder rung B7 §5.5: the generated empty-lattice
!>        (″vacuum lead″) potential parameters versus the **analytic
!>        spherical-Bessel free-electron result**.
!>
!>        ## Why this test is shaped the way it is
!>
!>        B7 §1.6 is explicit that the generator must drive the code’s own
!>        radial solver with V(r) = const, and that the analytic result is
!>        the **oracle, not the implementation**. This file is therefore the
!>        only place in the tree where spherical Bessel functions appear in
!>        connection with the vacuum lead, and they are written out here in
!>        closed form -- independently of `source/vacuum_lead.f90`, and
!>        independently of the `strux_lib` Bessel routines, so that the
!>        oracle cannot accidentally share a bug with the thing it is
!>        checking.
!>
!>        ## The physics being pinned
!>
!>        Inside an empty sphere of radius w holding a constant potential
!>        V0, the regular solution of the radial equation at energy E is
!>        j_l(kappa r) with kappa = sqrt(E - V0). `self%potpar` fixes E_nu
!>        per l by solving to a prescribed logarithmic derivative
!>        D_l = tan(pi*(0.5 - pnu)) at r = w. So the generated E_nu must
!>        satisfy, exactly,
!>
!>            w * j_l’(kappa w) / j_l(kappa w) == D_l,   kappa = sqrt(E_nu - V0)
!>
!>        Note the oracle is evaluated *at the generated E_nu*. It is not a
!>        root-find: D_l(E) has poles wherever j_l(kappa w) vanishes, and
!>        bracketing across those poles is fragile. Substituting the
!>        solver’s own answer into the analytic log-derivative is both
!>        exact and robust.
!>
!>        Tests:
!>          1. Bessel oracle self-check: the closed forms used below
!>             reproduce the standard recurrence j_l’ = j_{l-1} - (l+1)/x
!>             j_l and the Rayleigh formulae, so a typo in the oracle fails
!>             here rather than silently excusing the generator.
!>          2. E_nu satisfies the analytic log-derivative condition, for
!>             every l, over the energy window where ASA is expected to hold
!>             (§1.6) -- swept by varying the boundary condition pnu.
!>          3. Rigid shift: V(r) = const shifted by dV shifts every energy
!>             parameter (enu, C, VL) by exactly dV and leaves the
!>             shift-invariant ones (srdel, qpar, ppar) unchanged. This is
!>             what ″rigidly shifted to the vacuum level″ means and it is
!>             the property the alignment machinery (B7.4) will lean on.
!>          4. Geometry dependence: the generated set genuinely varies with
!>             the empty-sphere radius w -- i.e. this is a generator, not a
!>             table of constants (§1.6).
!>          5. Non-magnetic: the two spin channels are identical.
!>          6. The E_F-below-band-onset margin check warns (returns .false.)
!>             when E_F is pushed up to the onset and passes for a good
!>             metal.
!>
!>        Exits non-zero (error stop) on any failure so ctest registers a fail.
!------------------------------------------------------------------------------
program test_vacuum_lead
   use precision_mod, only: rp
   use self_mod, only: self
   use vacuum_lead_mod, only: vacuum_lead, vacuum_lead_ry_per_ev
   implicit none

   logical :: failed

   !> fccCu001 empty-sphere Wigner-Seitz radius, bohr. Dimensional
   !> (`wssurf` convention, CONVENTIONS_MADELUNG.md C0):
   !> sws*alat*ang2au = 2.66899 bohr = wav*ang2au. NOT the dimensionless
   !> relative sws ~ 0.391.
   real(rp), parameter :: w_cu = 2.66899_rp

   failed = .false.

   call test_bessel_oracle_selfcheck()
   call test_enu_matches_bessel_logderiv()
   call test_rigid_shift()
   call test_geometry_dependence()
   call test_spin_degenerate()
   call test_fermi_margin_warning()
   call test_dump()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   !===========================================================================
   ! The analytic oracle: spherical Bessel functions of the first kind,
   ! closed form for l = 0..3 (Rayleigh formulae). Deliberately independent
   ! of anything in source/.
   !===========================================================================
   pure function sph_j(l, x) result(j)
      integer, intent(in) :: l
      real(rp), intent(in) :: x
      real(rp) :: j, s, c

      s = sin(x)
      c = cos(x)
      select case (l)
      case (0)
         j = s/x
      case (1)
         j = s/x**2 - c/x
      case (2)
         j = (3.0_rp/x**3 - 1.0_rp/x)*s - 3.0_rp*c/x**2
      case (3)
         j = (15.0_rp/x**4 - 6.0_rp/x**2)*s - (15.0_rp/x**3 - 1.0_rp/x)*c
      case default
         j = 0.0_rp
      end select
   end function sph_j

   !> d j_l / dx via the standard recurrence j_l’ = j_{l-1} - (l+1)/x * j_l,
   !> with j_0’ = -j_1.
   pure function sph_jp(l, x) result(jp)
      integer, intent(in) :: l
      real(rp), intent(in) :: x
      real(rp) :: jp

      if (l == 0) then
         jp = -sph_j(1, x)
      else
         jp = sph_j(l - 1, x) - real(l + 1, rp)/x*sph_j(l, x)
      end if
   end function sph_jp

   !> The dimensionless logarithmic derivative D_l = w * j_l’(kw) / j_l(kw)
   !> = x * j_l’(x) / j_l(x) with x = kappa*w. This is exactly the quantity
   !> potpar’s boundary condition dnu = tan(pi*(0.5 - pnu)) prescribes.
   !>
   !> NON-relativistic: kappa = sqrt(E - V0). This is the textbook
   !> free-electron oracle named in B7 §1.6, and it is checked at a
   !> correspondingly loose tolerance, because the code’s solver is the
   !> *scalar-relativistic* one -- see `analytic_dnu_sr`.
   pure function analytic_dnu(l, energy, v0, w) result(d)
      integer, intent(in) :: l
      real(rp), intent(in) :: energy, v0, w
      real(rp) :: d, x

      x = sqrt(energy - v0)*w
      d = x*sph_jp(l, x)/sph_j(l, x)
   end function analytic_dnu

   !> Scalar-relativistic version of the same oracle.
   !>
   !> `RSEQSR` -- the solver `potpar` drives -- is the scalar-relativistic
   !> radial solver, and it carries the relativistic mass enhancement
   !> explicitly as `TMCR = (C - (V - 2Z/r - E)/C)*r` with `C = 274.074`
   !> (source/self.f90:2345, 2434), i.e. `2 M c r` with
   !>
   !>     M(r) = 1 + (E - V(r)) / c^2
   !>
   !> For an empty sphere V is constant, so M is constant, the `M’/M`
   !> (Darwin-like) term drops out, and the large component still solves a
   !> spherical Bessel equation -- but with the relativistically corrected
   !> wavenumber
   !>
   !>     kappa_SR = sqrt( (E - V0) * M ),   M = 1 + (E - V0)/c^2
   !>
   !> This form removes the leading O((E-V0)/c^2) discrepancy against the
   !> non-relativistic oracle -- roughly two orders of magnitude, from
   !> ~1e-4 to ~3e-6 in relative log derivative at the top of the tested
   !> window -- and what remains is genuine higher-order SR content that a
   !> constant-M Bessel form cannot represent. Checking both forms is the
   !> point: agreement with the non-relativistic oracle to O(alpha^2)
   !> establishes the free-electron limit, and the tighter agreement with
   !> the SR form establishes that the residual is relativistic physics and
   !> not an error in the generator.
   pure function analytic_dnu_sr(l, energy, v0, w) result(d)
      integer, intent(in) :: l
      real(rp), intent(in) :: energy, v0, w
      !> Speed of light in the Rydberg atomic units RSEQSR uses
      !> (source/self.f90:2345).
      real(rp), parameter :: c_light = 274.074_rp
      real(rp) :: d, x, m

      m = 1.0_rp + (energy - v0)/c_light**2
      x = sqrt((energy - v0)*m)*w
      d = x*sph_jp(l, x)/sph_j(l, x)
   end function analytic_dnu_sr

   !> Boundary condition potpar was given, from pnu.
   pure function requested_dnu(pnu) result(d)
      real(rp), intent(in) :: pnu
      real(rp) :: d, pi

      pi = 4.0_rp*atan(1.0_rp)
      d = tan(pi*(0.5_rp - pnu))
   end function requested_dnu

   !===========================================================================
   ! 1. Oracle self-check.
   !===========================================================================
   subroutine test_bessel_oracle_selfcheck()
      real(rp), parameter :: tol = 1.0e-12_rp
      real(rp), dimension(3), parameter :: xs = [1.2_rp, 2.0_rp, 3.4_rp]
      real(rp) :: x, lhs, rhs, err
      integer :: i, l
      logical :: bad

      bad = .false.
      do i = 1, size(xs)
         x = xs(i)
         ! Downward recurrence identity: j_{l-1} + j_{l+1} = (2l+1)/x * j_l.
         do l = 1, 2
            lhs = sph_j(l - 1, x) + sph_j(l + 1, x)
            rhs = real(2*l + 1, rp)/x*sph_j(l, x)
            err = abs(lhs - rhs)
            if (err > tol) then
               write (*, '(a,i0,a,f6.2,a,es12.4)') 'FAIL Bessel recurrence l=', l, ' x=', x, ' err=', err
               bad = .true.
            end if
         end do
         ! Derivative identity: (2l+1) j_l’ = l j_{l-1} - (l+1) j_{l+1}.
         do l = 1, 2
            lhs = real(2*l + 1, rp)*sph_jp(l, x)
            rhs = real(l, rp)*sph_j(l - 1, x) - real(l + 1, rp)*sph_j(l + 1, x)
            err = abs(lhs - rhs)
            if (err > tol) then
               write (*, '(a,i0,a,f6.2,a,es12.4)') 'FAIL Bessel derivative l=', l, ' x=', x, ' err=', err
               bad = .true.
            end if
         end do
      end do

      if (bad) then
         failed = .true.
      else
         write (*, '(a)') 'ok  1. analytic spherical-Bessel oracle is self-consistent'
      end if
   end subroutine test_bessel_oracle_selfcheck

   !===========================================================================
   ! 2. THE RUNG: generated E_nu vs the analytic free-electron log derivative.
   !
   ! Swept over a range of boundary conditions pnu, which is the same thing as
   ! sweeping E_nu over the low-energy window where the ASA empty-lattice
   ! description is expected to hold (§1.6). pnu = l + 1 + f with f in
   ! (0, 1) and no radial nodes; f = 0.5 is the free-electron-like default,
   ! larger f pushes E_nu down, smaller f pushes it up.
   !===========================================================================
   subroutine test_enu_matches_bessel_logderiv()
      type(self) :: solver
      type(vacuum_lead) :: vac
      integer, parameter :: lmax = 2
      real(rp), parameter :: v0 = 0.0_rp
      !> Non-relativistic oracle: agreement is limited by the O(alpha^2)
      !> scalar-relativistic content the solver carries and this form does
      !> not, which reaches ~1e-4 at the top of the window. 1e-3 confirms
      !> the free-electron limit without pretending to more.
      real(rp), parameter :: tol_nr = 1.0e-3_rp
      !> Scalar-relativistic oracle: the leading correction is removed, so
      !> the residual is higher-order SR only, ~3e-6 worst case. 1e-5 is a
      !> genuinely tight bar -- an error in the generator’s boundary
      !> condition, mesh or radius would blow straight through it.
      real(rp), parameter :: tol_sr = 1.0e-5_rp
      !> Sweep the boundary condition, which sweeps E_nu across the window
      !> where the ASA empty-lattice description is expected to hold. The
      !> admissible range is dnu < l (dnu = l is the band bottom), so the
      !> sweep is expressed as a step INSIDE the bottom, uniformly in l --
      !> the same parametrisation the generator’s default uses.
      real(rp), dimension(4), parameter :: step = [0.25_rp, 0.5_rp, 1.0_rp, 2.0_rp]
      real(rp), dimension(0:lmax) :: pnu
      real(rp) :: d_nr, d_sr, d_req, err_nr, err_sr, x, pi
      integer :: i, l
      logical :: bad

      pi = 4.0_rp*atan(1.0_rp)
      bad = .false.
      do i = 1, size(step)
         do l = 0, lmax
            ! dnu = l - step  ->  pnu
            pnu(l) = real(l + 1, rp) + 0.5_rp - atan(real(l, rp) - step(i))/pi
         end do

         call vac%generate(solver, ws_r=w_cu, v0=v0, lmax=lmax, pnu=pnu)

         do l = 0, lmax
            ! Only meaningful for E_nu above the constant potential floor:
            ! below it kappa is imaginary and the regular solution is a
            ! modified spherical Bessel function, outside the window this
            ! rung is specified for. dnu < l guarantees we are above it.
            if (vac%enu(l, 1) <= v0) then
               write (*, '(a,i0,a,f6.3,a,f12.6)') 'FAIL enu not above V0 for l=', l, &
                  ' step=', step(i), ' enu=', vac%enu(l, 1)
               bad = .true.
               cycle
            end if

            x = sqrt(vac%enu(l, 1) - v0)*w_cu
            d_nr = analytic_dnu(l, vac%enu(l, 1), v0, w_cu)
            d_sr = analytic_dnu_sr(l, vac%enu(l, 1), v0, w_cu)
            d_req = requested_dnu(pnu(l))
            err_nr = abs(d_nr - d_req)/max(1.0_rp, abs(d_req))
            err_sr = abs(d_sr - d_req)/max(1.0_rp, abs(d_req))

            write (*, '(a,i0,a,f5.2,a,f11.7,a,f8.5,a,f11.7,a,es9.2,a,es9.2)') &
               '    l=', l, ' dnu=', d_req, '  enu=', vac%enu(l, 1), &
               '  x=kw=', x, '  D_req=', d_req, &
               '  relerr(NR)=', err_nr, '  relerr(SR)=', err_sr

            if (err_nr > tol_nr) then
               write (*, '(a,i0,a,f6.2,a,es12.4)') 'FAIL NR log-derivative mismatch l=', l, &
                  ' dnu=', d_req, ' relerr=', err_nr
               bad = .true.
            end if
            if (err_sr > tol_sr) then
               write (*, '(a,i0,a,f6.2,a,es12.4)') 'FAIL SR log-derivative mismatch l=', l, &
                  ' dnu=', d_req, ' relerr=', err_sr
               bad = .true.
            end if
            ! The SR form must be the better of the two -- that is the
            ! statement that the residual is relativistic physics rather
            ! than a generator error that happens to be small.
            if (err_sr > err_nr) then
               write (*, '(a,i0,a,2es12.4)') 'FAIL SR oracle is not tighter than NR, l=', l, &
                  ' NR/SR=', err_nr, err_sr
               bad = .true.
            end if
         end do
      end do

      if (bad) then
         failed = .true.
      else
         write (*, '(a)') 'ok  2. generated E_nu reproduces the analytic j_l log derivative'
         write (*, '(a)') '        (NR to O(alpha^2); SR form tighter by ~2 orders, as it must be)'
      end if
   end subroutine test_enu_matches_bessel_logderiv

   !===========================================================================
   ! 3. Rigid shift to the vacuum level.
   !===========================================================================
   subroutine test_rigid_shift()
      type(self) :: solver
      type(vacuum_lead) :: vac0, vac1
      integer, parameter :: lmax = 2
      real(rp), parameter :: dv = 0.37_rp
      real(rp), parameter :: tol = 1.0e-9_rp
      real(rp) :: err
      integer :: l
      logical :: bad

      call vac0%generate(solver, ws_r=w_cu, v0=0.0_rp, lmax=lmax)
      call vac1%generate(solver, ws_r=w_cu, v0=dv, lmax=lmax)

      bad = .false.
      do l = 0, lmax
         err = abs((vac1%enu(l, 1) - vac0%enu(l, 1)) - dv)
         if (err > tol) then
            write (*, '(a,i0,a,es12.4)') 'FAIL enu not rigidly shifted, l=', l, ' err=', err
            bad = .true.
         end if
         err = abs((vac1%c(l, 1) - vac0%c(l, 1)) - dv)
         if (err > tol) then
            write (*, '(a,i0,a,es12.4)') 'FAIL C not rigidly shifted, l=', l, ' err=', err
            bad = .true.
         end if
         err = abs((vac1%vl(l, 1) - vac0%vl(l, 1)) - dv)
         if (err > tol) then
            write (*, '(a,i0,a,es12.4)') 'FAIL VL not rigidly shifted, l=', l, ' err=', err
            bad = .true.
         end if
         ! Shift-invariant parameters.
         err = abs(vac1%srdel(l, 1) - vac0%srdel(l, 1))
         if (err > tol) then
            write (*, '(a,i0,a,es12.4)') 'FAIL SRDEL changed under rigid shift, l=', l, ' err=', err
            bad = .true.
         end if
         err = abs(vac1%qpar(l, 1) - vac0%qpar(l, 1))
         if (err > tol) then
            write (*, '(a,i0,a,es12.4)') 'FAIL QPAR changed under rigid shift, l=', l, ' err=', err
            bad = .true.
         end if
         err = abs(vac1%ppar(l, 1) - vac0%ppar(l, 1))
         if (err > tol) then
            write (*, '(a,i0,a,es12.4)') 'FAIL PPAR changed under rigid shift, l=', l, ' err=', err
            bad = .true.
         end if
      end do

      if (bad) then
         failed = .true.
      else
         write (*, '(a,f6.3,a)') 'ok  3. V0 shift of ', dv, &
            ' Ry shifts enu/C/VL rigidly, leaves srdel/qpar/ppar fixed'
      end if
   end subroutine test_rigid_shift

   !===========================================================================
   ! 4. Geometry dependence: this is a generator, not a table (§1.6).
   !    A larger empty sphere must lower every band centre -- the standard
   !    particle-in-a-sphere scaling E ~ 1/w^2 -- and the analytic oracle
   !    must still be satisfied at the new radius.
   !===========================================================================
   subroutine test_geometry_dependence()
      type(self) :: solver
      type(vacuum_lead) :: small, big
      integer, parameter :: lmax = 2
      real(rp), parameter :: w_small = 2.30_rp, w_big = 3.10_rp
      real(rp), parameter :: tol = 1.0e-5_rp
      !> Particle-in-a-sphere scaling: at FIXED log derivative the
      !> dimensionless x = kappa*w is fixed, so E - V0 = (x/w)^2 scales as
      !> 1/w^2 exactly. The default dnu is w-independent, so this ratio is
      !> a sharp prediction, not a hand-wave. Allow 1% for the SR mass
      !> enhancement, which breaks the exact scaling at O(alpha^2).
      real(rp), parameter :: tol_scaling = 1.0e-2_rp
      real(rp) :: err, d_analytic, d_requested, ratio, expected
      integer :: l
      logical :: bad

      call small%generate(solver, ws_r=w_small, v0=0.0_rp, lmax=lmax)
      call big%generate(solver, ws_r=w_big, v0=0.0_rp, lmax=lmax)

      bad = .false.
      expected = (w_small/w_big)**2
      do l = 0, lmax
         if (.not. (big%enu(l, 1) < small%enu(l, 1))) then
            write (*, '(a,i0,2f14.8)') 'FAIL enu did not decrease with larger w, l=', l, &
               small%enu(l, 1), big%enu(l, 1)
            bad = .true.
            cycle
         end if
         ratio = big%enu(l, 1)/small%enu(l, 1)
         if (abs(ratio - expected)/expected > tol_scaling) then
            write (*, '(a,i0,a,f10.6,a,f10.6)') 'FAIL 1/w^2 scaling violated, l=', l, &
               ' ratio=', ratio, ' expected=', expected
            bad = .true.
         end if
         ! Oracle must hold at the new geometry too.
         d_analytic = analytic_dnu_sr(l, big%enu(l, 1), 0.0_rp, w_big)
         d_requested = requested_dnu(big%pnu(l))
         err = abs(d_analytic - d_requested)/max(1.0_rp, abs(d_requested))
         if (err > tol) then
            write (*, '(a,i0,a,es12.4)') 'FAIL oracle fails at w_big, l=', l, ' relerr=', err
            bad = .true.
         end if
      end do

      if (bad) then
         failed = .true.
      else
         write (*, '(a,f5.2,a,f5.2,a)') 'ok  4. set varies with geometry (w = ', w_small, &
            ' -> ', w_big, ' bohr) and oracle holds at both'
      end if
   end subroutine test_geometry_dependence

   !===========================================================================
   ! 5. Vacuum is non-magnetic: both spin channels identical.
   !===========================================================================
   subroutine test_spin_degenerate()
      type(self) :: solver
      type(vacuum_lead) :: vac
      integer, parameter :: lmax = 2
      integer :: l
      logical :: bad

      call vac%generate(solver, ws_r=w_cu, v0=0.0_rp, lmax=lmax)

      bad = .false.
      do l = 0, lmax
         if (vac%enu(l, 1) /= vac%enu(l, 2) .or. vac%c(l, 1) /= vac%c(l, 2) .or. &
             vac%srdel(l, 1) /= vac%srdel(l, 2)) then
            write (*, '(a,i0)') 'FAIL spin channels differ for l=', l
            bad = .true.
         end if
      end do

      if (bad) then
         failed = .true.
      else
         write (*, '(a)') 'ok  5. both spin channels identical (vacuum is non-magnetic)'
      end if
   end subroutine test_spin_degenerate

   !===========================================================================
   ! 6. E_F margin below the vacuum band onset (§1.6 honest limitation).
   !===========================================================================
   subroutine test_fermi_margin_warning()
      type(self) :: solver
      type(vacuum_lead) :: vac
      integer, parameter :: lmax = 2
      real(rp) :: onset, margin_ry
      logical :: ok_good, ok_bad
      logical :: bad

      call vac%generate(solver, ws_r=w_cu, v0=0.0_rp, lmax=lmax)
      onset = vac%band_onset()
      margin_ry = 1.0_rp*vacuum_lead_ry_per_ev

      ! Good metal: E_F comfortably (4 eV) below the onset.
      ok_good = vac%check_fermi_margin(onset - 4.0_rp*vacuum_lead_ry_per_ev, margin=margin_ry)
      ! Low work function: E_F only 0.2 eV below the onset -- must warn.
      ok_bad = vac%check_fermi_margin(onset - 0.2_rp*vacuum_lead_ry_per_ev, margin=margin_ry)

      bad = .false.
      if (.not. ok_good) then
         write (*, '(a)') 'FAIL margin check rejected a comfortable 4 eV margin'
         bad = .true.
      end if
      if (ok_bad) then
         write (*, '(a)') 'FAIL margin check accepted a 0.2 eV margin without warning'
         bad = .true.
      end if

      if (bad) then
         failed = .true.
      else
         write (*, '(a,f10.6,a)') 'ok  6. band-onset margin check (onset C_s = ', onset, &
            ' Ry): passes at 4 eV, warns at 0.2 eV'
      end if
   end subroutine test_fermi_margin_warning

   !===========================================================================
   ! 7. The dump routine runs and reports one row per (spin, l).
   !===========================================================================
   subroutine test_dump()
      type(self) :: solver
      type(vacuum_lead) :: vac
      integer, parameter :: lmax = 2
      character(len=*), parameter :: fname = 'test_vacuum_lead_dump.tmp'
      character(len=256) :: line
      integer :: unit, ios, nlines

      call vac%generate(solver, ws_r=w_cu, v0=0.0_rp, lmax=lmax)

      open (newunit=unit, file=fname, status='replace', action='write')
      call vac%dump(unit=unit)
      close (unit)

      nlines = 0
      open (newunit=unit, file=fname, status='old', action='read')
      do
         read (unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         nlines = nlines + 1
      end do
      close (unit, status='delete')

      ! 2 spins * (lmax+1) l values, plus header and summary lines.
      if (nlines < 2*(lmax + 1)) then
         write (*, '(a,i0)') 'FAIL dump produced too few lines: ', nlines
         failed = .true.
      else
         write (*, '(a,i0,a)') 'ok  7. dump routine wrote ', nlines, ' lines'
      end if
   end subroutine test_dump

end program test_vacuum_lead
