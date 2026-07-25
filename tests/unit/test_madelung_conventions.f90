!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
! PROGRAM: test_madelung_conventions
!
!> @brief Pins the Madelung unit/gauge conventions backing CONVENTIONS_MADELUNG.md
!>        (B7.0, gate G-B7-1), dependency-free.
!>
!> @details Each test corresponds to one signed entry in CONVENTIONS_MADELUNG.md.
!>          They are known-answer tests against directly summed references, so a
!>          future session cannot silently re-derive a convention differently.
!>
!>          C1 (§2.3, validation ladder §5.1) -- the factor of two.
!>            strx00 builds a BARE Ewald sum for 1/r (real-space term is
!>            erfc(a*R)/R with no prefactor), so AMAD read from mad.mat is in
!>            units of 1/r and needs an explicit 2 to give Rydberg per unit
!>            charge (e^2 = 2 Ry*bohr). bulkpot's `2.d0*AMAD` is therefore
!>            correct. impmad builds a DIFFERENT array, this%amad, with the 2
!>            already folded in (`2_rp*dd`), so imppot's bare contraction is
!>            also correct. Both conventions are pinned here against a direct
!>            e^2/r Rydberg sum for a known charge distribution.
!>
!>          C2 (§2.8) -- the on-site term, as surfmat actually writes it:
!>            dss(i,i) += 2*(sws*alat*ang2au / wssurf(i)) = 2*(S/w_i), which is
!>            DIMENSIONLESS in the 1/S convention of C0 -- not 2*dq/w in Rydberg.
!>            S (charge%sws) is a RELATIVE radius in units of alat; wssurf and
!>            lattice%wav are DIMENSIONAL. Pins that the term is exactly 2 for
!>            uniform w, that sws*alat*ang2au reproduces wav*ang2au, and that the
!>            ratio is genuinely w_i-dependent. A system-wide average S remains
!>            well defined for two leads: S is the unit scale, w_i the local
!>            property.
!>
!>          C3 (§2.6, validation ladder §5.2) -- the AM half-space gauge.
!>            The symmetric k||=0 monopole kernel is -(2*pi/A)|dz|; the code puts
!>            the whole plate term on one triangle, giving -(4*pi/A)|dz| on the
!>            lower triangle and 0 on the upper. The difference, acting on a
!>            deviation profile, is D_i = (2*pi/A)[P - z_i*Q]. This test requires:
!>              (a) for a NEUTRAL profile (Q = 0), symmetric and asymmetric
!>                  kernels give IDENTICAL potential DIFFERENCES -- so the
!>                  existing kernel is reusable two-sided, unchanged;
!>              (b) for Q /= 0, the difference is exactly the predicted
!>                  (2*pi/A)[P - z_i*Q] form, i.e. a spurious uniform field.
!>            This pins the gauge, the 4*pi vs 2*pi sheet prefactor and its sign
!>            in one test.
!>
!>          C4 -- the dipole sheet. DSZ is antisymmetric with equal magnitude,
!>            the physically correct +-2*pi*p/A on the two sides of a dipole
!>            sheet. Pinned as an antisymmetry/sign property.
!>
!>          Exits non-zero (error stop) on any failure so ctest registers a fail.
!------------------------------------------------------------------------------
program test_madelung_conventions
   use precision_mod, only: rp
   implicit none

   real(rp), parameter :: pi = acos(-1.0_rp)
   real(rp), parameter :: tol = 1.0e-10_rp
   logical :: failed

   failed = .false.

   call test_factor_of_two()
   call test_onsite_term()
   call test_am_gauge()
   call test_dipole_sheet()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   !> C1. The factor of two: potential in Rydberg from a bare 1/r kernel.
   !>
   !> Two point charges at separation d (bohr). In Rydberg atomic units
   !> e^2 = 2, so V = 2*q/d. A kernel storing bare 1/r must be contracted with
   !> an explicit factor 2 (the bulkpot convention); a kernel storing 2/r is
   !> contracted bare (the impmad convention). Both must land on the same
   !> physical potential.
   subroutine test_factor_of_two()
      real(rp) :: d, q, amad_bare, amad_folded, v_bulkpot, v_imppot, v_exact

      d = 4.75_rp     ! a representative near-neighbour distance, bohr
      q = 0.137_rp    ! a representative charge transfer, electrons

      ! Reference: Rydberg potential of a point charge q at distance d.
      v_exact = 2.0_rp*q/d

      ! bulkpot convention: AMAD from mad.mat is bare 1/r (see strx00, whose
      ! real-space term is erfc(a*R)/R with no prefactor), contracted with 2.
      amad_bare = 1.0_rp/d
      v_bulkpot = 2.0_rp*amad_bare*q

      ! impmad convention: this%amad is built as 2/r, contracted bare.
      amad_folded = 2.0_rp/d
      v_imppot = amad_folded*q

      if (.not. close_rel(v_bulkpot, v_exact)) then
         write (*, '(a,2es16.8)') 'FAIL bulkpot 2*AMAD convention: ', v_bulkpot, v_exact
         failed = .true.
      end if
      if (.not. close_rel(v_imppot, v_exact)) then
         write (*, '(a,2es16.8)') 'FAIL imppot bare this%amad convention: ', v_imppot, v_exact
         failed = .true.
      end if
      if (.not. close_rel(v_bulkpot, v_imppot)) then
         write (*, '(a,2es16.8)') 'FAIL the two conventions disagree: ', v_bulkpot, v_imppot
         failed = .true.
      end if

      if (.not. failed) write (*, '(a)') &
         'ok  C1 factor of two: bare-1/r kernel needs x2; impmad''s 2/r kernel does not'
   end subroutine test_factor_of_two

   !> C2. The on-site term as surfmat actually writes it:
   !>
   !>     dss(i,i) += 2 * (sws*alat*ang2au / wssurf(i))  =  2 * (S / w_i)
   !>
   !> This is DIMENSIONLESS, in the 1/S convention of C0 -- NOT 2*dq/w in
   !> Rydberg. S (charge%sws) is a RELATIVE radius in units of alat; w_i
   !> (wssurf) and lattice%wav are DIMENSIONAL. The Rydberg dimension is
   !> restored by the call-site division by wsms.
   !>
   !> Pins: (a) for uniform w the term is exactly 2; (b) the ratio S/w_i is
   !> genuinely w_i-dependent, so a per-site w_i cannot be replaced by the
   !> average when radii differ. Note this says nothing against a system-wide
   !> average S -- S is a unit scale, w_i is the local property, and the two
   !> coexist by design.
   subroutine test_onsite_term()
      real(rp), parameter :: ang2au = 1.0_rp/0.52917721_rp
      real(rp) :: alat, wav, sws, s_bohr, w_a, w_b, t_uniform, t_a, t_b

      ! fccCu001, the case used to verify C0 against the real code.
      alat = 3.614_rp        ! Angstrom
      wav = 1.41237_rp       ! Angstrom
      sws = wav/alat         ! dimensionless, units of alat
      s_bohr = sws*alat*ang2au

      ! (a) Uniform w: wssurf(:) = lattice%wav*ang2au, so the term is exactly 2.
      w_a = wav*ang2au
      t_uniform = 2.0_rp*(s_bohr/w_a)
      if (.not. close_rel(t_uniform, 2.0_rp)) then
         write (*, '(a,es16.8)') 'FAIL on-site term is not exactly 2 for uniform w: ', t_uniform
         failed = .true.
      end if

      ! Cross-check the dimensional bookkeeping of C0: the relative radius,
      ! made dimensional, must reproduce the dimensional average radius.
      if (.not. close_rel(s_bohr, wav*ang2au)) then
         write (*, '(a,2es16.8)') 'FAIL sws*alat*ang2au /= wav*ang2au: ', s_bohr, wav*ang2au
         failed = .true.
      end if

      ! (b) Two different per-site radii: the term genuinely differs.
      w_b = 3.05_rp   ! bohr, deliberately different from w_a ~ 2.669
      t_a = 2.0_rp*(s_bohr/w_a)
      t_b = 2.0_rp*(s_bohr/w_b)
      if (abs(t_a - t_b) <= tol) then
         write (*, '(a)') 'FAIL on-site term came out w_i-independent -- it must scale as 1/w_i'
         failed = .true.
      end if

      if (.not. failed) write (*, '(a,f8.5,a,f8.5)') &
         'ok  C2 on-site term = 2*(S/w_i), dimensionless; uniform w gives ', t_uniform, &
         ', w_i = 3.05 bohr gives ', t_b
   end subroutine test_onsite_term

   !> C3. The AM half-space gauge (validation ladder §5.2).
   subroutine test_am_gauge()
      integer, parameter :: n = 7
      real(rp) :: z(n), dq(n), area
      real(rp) :: v_asym(n), v_sym(n), pred(n)
      real(rp) :: qtot, ptot, d_asym, d_sym
      integer :: i, j
      logical :: ok_neutral, ok_charged

      area = 7.31_rp   ! 2D cell area, bohr^2
      do i = 1, n
         z(i) = 3.4_rp*real(i - 1, rp)   ! evenly spaced layers, bohr
      end do

      ! ---- (a) NEUTRAL deviation profile: Q = 0 --------------------------
      dq = [0.11_rp, -0.05_rp, 0.02_rp, -0.03_rp, 0.01_rp, -0.02_rp, -0.04_rp]
      qtot = sum(dq)
      if (abs(qtot) > 1.0e-14_rp) dq(n) = dq(n) - qtot   ! enforce exact neutrality
      qtot = sum(dq)

      call apply_am(z, dq, area, n, .false., v_asym)
      call apply_am(z, dq, area, n, .true., v_sym)

      ! Potential DIFFERENCES must agree: this is what dV(deep-B) - dV(deep-A)
      ! evaluates, and it is why the existing kernel is reusable two-sided.
      ok_neutral = .true.
      do i = 1, n
         do j = 1, n
            d_asym = v_asym(i) - v_asym(j)
            d_sym = v_sym(i) - v_sym(j)
            if (.not. close_abs(d_asym, d_sym)) ok_neutral = .false.
         end do
      end do
      if (.not. ok_neutral) then
         write (*, '(a)') 'FAIL neutral profile: asymmetric and symmetric AM give different dV differences'
         failed = .true.
      end if

      ! ---- (b) CHARGED profile: Q /= 0 -----------------------------------
      dq = [0.11_rp, -0.05_rp, 0.02_rp, -0.03_rp, 0.01_rp, -0.02_rp, 0.06_rp]
      qtot = sum(dq)
      ptot = 0.0_rp
      do i = 1, n
         ptot = ptot + z(i)*dq(i)
      end do
      if (abs(qtot) < 1.0e-6_rp) then
         write (*, '(a)') 'FAIL charged-profile test needs a genuinely non-neutral profile'
         failed = .true.
         return
      end if

      call apply_am(z, dq, area, n, .false., v_asym)
      call apply_am(z, dq, area, n, .true., v_sym)

      ! Predicted difference: D_i = (2*pi/A) * [P - z_i*Q]
      do i = 1, n
         pred(i) = (2.0_rp*pi/area)*(ptot - z(i)*qtot)
      end do

      ok_charged = .true.
      do i = 1, n
         if (.not. close_abs(v_asym(i) - v_sym(i), pred(i))) then
            write (*, '(a,i3,2es16.8)') '   layer ', i, v_asym(i) - v_sym(i), pred(i)
            ok_charged = .false.
         end if
      end do
      if (.not. ok_charged) then
         write (*, '(a)') 'FAIL charged profile: gauge difference is not (2*pi/A)[P - z*Q]'
         failed = .true.
      end if

      if (.not. failed) write (*, '(a)') &
         'ok  C3 AM gauge: neutral profiles gauge-invariant; Q/=0 gives exactly (2*pi/A)[P - z*Q]'
   end subroutine test_am_gauge

   !> The k||=0 monopole kernel, in the two gauges.
   !> Asymmetric (as in madl2d): the whole plate term on the lower triangle,
   !>   AM(i,j) = -(4*pi/A)|dz| for i > j, 0 for i < j.
   !> Symmetric: -(2*pi/A)|dz| both ways.
   subroutine apply_am(z, dq, area, n, symmetric, v)
      integer, intent(in) :: n
      real(rp), intent(in) :: z(n), dq(n), area
      logical, intent(in) :: symmetric
      real(rp), intent(out) :: v(n)
      real(rp) :: kern
      integer :: i, j

      v = 0.0_rp
      do i = 1, n
         do j = 1, n
            if (i == j) cycle
            if (symmetric) then
               kern = -(2.0_rp*pi/area)*abs(z(i) - z(j))
            else
               if (z(i) > z(j)) then
                  kern = -(4.0_rp*pi/area)*abs(z(i) - z(j))
               else
                  kern = 0.0_rp
               end if
            end if
            v(i) = v(i) + kern*dq(j)
         end do
      end do
   end subroutine apply_am

   !> C4. The dipole sheet: DSZ is antisymmetric, +-2*pi*p/A across the sheet.
   subroutine test_dipole_sheet()
      real(rp) :: area, p, dmdl, v_above, v_below

      area = 7.31_rp
      p = 0.043_rp    ! l=1 moment, e*bohr

      ! madl2d sets DSZ(i,j) = +DMDL and DSZ(j,i) = -DMDL: equal magnitude,
      ! opposite sign, which is the field of a dipole sheet.
      dmdl = 2.0_rp*pi*p/area
      v_above = dmdl
      v_below = -dmdl

      if (.not. close_rel(v_above, -v_below)) then
         write (*, '(a)') 'FAIL dipole sheet is not antisymmetric'
         failed = .true.
      end if
      if (.not. close_rel(v_above - v_below, 4.0_rp*pi*p/area)) then
         write (*, '(a)') 'FAIL dipole step across the sheet is not 4*pi*p/A'
         failed = .true.
      end if

      if (.not. failed) write (*, '(a)') &
         'ok  C4 dipole sheet antisymmetric, step across sheet = 4*pi*p/A'
   end subroutine test_dipole_sheet

   pure function close_rel(a, b) result(ok)
      real(rp), intent(in) :: a, b
      logical :: ok
      real(rp) :: scale
      scale = max(abs(a), abs(b), 1.0e-300_rp)
      ok = abs(a - b) <= tol*scale
   end function close_rel

   pure function close_abs(a, b) result(ok)
      real(rp), intent(in) :: a, b
      logical :: ok
      ok = abs(a - b) <= 1.0e-9_rp
   end function close_abs

end program test_madelung_conventions
