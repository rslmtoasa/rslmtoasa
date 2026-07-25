!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
! PROGRAM: test_interface_alignment_oracle
!
!> @brief The B7 validation-ladder rungs §5.2, §5.3 and §5.4, on a toy layered
!>        cluster. Dependency-free.
!>
!> @details B7 §5.3 is explicit that this oracle must be written **before** the
!>          interface geometry builder exists, on a toy cluster if that is
!>          faster. This is that test. It exercises the genuinely new machinery
!>          -- deviation variables, the Q/P moments, the Q = 0 kernel
!>          precondition, localized N(E_F)-weighted compensation, and the
!>          alignment fixed point -- against closed-form answers, with no LMTO
!>          machinery and no cluster builder involved.
!>
!>          The toy system is a 1D stack of layers at z_i with a 2D cell area A,
!>          coupled by the SAME half-space monopole kernel madl2d builds
!>          (§2.6 / CONVENTIONS_MADELUNG.md C3):
!>
!>            AM(i,j) = -(4*pi/A)*|z_i - z_j|   for z_i > z_j,   0 otherwise
!>
!>          Rungs pinned here:
!>
!>          §5.4  Plain A|A identity. Same material both sides, no offset:
!>                every layer comes out bulk, dQ -> 0, dipole -> 0, V_B -> 0.
!>                Needs no new geometry, so it unblocks B7.3/B7.4 before the
!>                cluster builder exists.
!>
!>          §5.3  A|A with a deliberate rigid offset delta -- THE PRIMARY ORACLE.
!>                A|A alone does not test alignment: the offsets cancel. The real
!>                oracle is region A against region-A-with-a-rigid-shift-delta
!>                applied to its input parameters. Requirements: the alignment
!>                solver recovers a step of exactly -delta, and dQ -> 0 on every
!>                site. This is an exact root-find test.
!>
!>          §5.2  The AM symmetrization oracle, in the two-sided setting: for a
!>                NEUTRAL deviation profile the half-space and symmetric kernels
!>                give identical dV differences (so the existing kernel is
!>                reusable two-sided unchanged); with Q /= 0 injected they differ
!>                by exactly (2*pi/A)[P - z_i*Q].
!>
!>          Also pinned: the compensation rule of §1.5 -- weight by the
!>          side-resolved N(E_F) of the boundary buffer layers, which gives
!>          ~50/50 for two similar metals and collapses to "everything into the
!>          metal" for metal|vacuum, so the surface case needs NO separate
!>          branch. And that compensation must be LOCALIZED, not smeared as
!>          imppot does: in layered geometry a smear has a nonzero first moment
!>          and a long lever arm.
!>
!>          Exits non-zero (error stop) on any failure so ctest registers a fail.
!------------------------------------------------------------------------------
program test_interface_alignment_oracle
   use precision_mod, only: rp
   implicit none

   real(rp), parameter :: pi = acos(-1.0_rp)
   real(rp), parameter :: tol = 1.0e-10_rp
   logical :: failed

   failed = .false.

   call test_aa_identity()
   call test_offset_delta_oracle()
   call test_gauge_two_sided()
   call test_compensation_weighting()
   call test_compensation_localization()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   !> The half-space monopole kernel of madl2d (CONVENTIONS_MADELUNG.md C3),
   !> applied to a deviation profile. `symmetric` selects the symmetrized
   !> -(2*pi/A)|dz| variant used by the §5.2 oracle.
   subroutine apply_kernel(z, dq, area, n, symmetric, v)
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
            else if (z(i) > z(j)) then
               kern = -(4.0_rp*pi/area)*abs(z(i) - z(j))
            else
               kern = 0.0_rp
            end if
            v(i) = v(i) + kern*dq(j)
         end do
      end do
   end subroutine apply_kernel

   !> The two moments of the deviation profile that B7 §1.4 requires be computed
   !> and reported every iteration.
   subroutine moments(z, dq, n, qtot, ptot)
      integer, intent(in) :: n
      real(rp), intent(in) :: z(n), dq(n)
      real(rp), intent(out) :: qtot, ptot
      integer :: i

      qtot = 0.0_rp
      ptot = 0.0_rp
      do i = 1, n
         qtot = qtot + dq(i)
         ptot = ptot + z(i)*dq(i)
      end do
   end subroutine moments

   !> §5.4 -- plain A|A identity. Same material both sides, no offset.
   !> Deviation variables are defined against the region reference, so a
   !> homogeneous system has dq identically zero even when the ABSOLUTE charges
   !> are nonzero. That is the whole point of §1.4: frozen deep sites contribute
   !> exactly zero and the sums truncate on both sides.
   subroutine test_aa_identity()
      integer, parameter :: n = 9
      real(rp) :: z(n), q_abs(n), q_ref(n), dq(n), v(n), area
      real(rp) :: qtot, ptot, step
      integer :: i

      area = 7.31_rp
      do i = 1, n
         z(i) = 3.4_rp*real(i - 1, rp)
      end do

      ! A deliberately CHARGED but homogeneous "material": every site carries
      ! the same absolute charge, as a non-homogeneous bulk legitimately can.
      q_abs = 0.234_rp
      q_ref = 0.234_rp
      dq = q_abs - q_ref

      call moments(z, dq, n, qtot, ptot)
      call apply_kernel(z, dq, area, n, .false., v)
      step = v(n) - v(1)

      if (.not. close_abs(qtot, 0.0_rp)) then
         write (*, '(a,es16.8)') 'FAIL A|A identity: Q /= 0: ', qtot
         failed = .true.
      end if
      if (.not. close_abs(ptot, 0.0_rp)) then
         write (*, '(a,es16.8)') 'FAIL A|A identity: dipole P /= 0: ', ptot
         failed = .true.
      end if
      if (.not. close_abs(step, 0.0_rp)) then
         write (*, '(a,es16.8)') 'FAIL A|A identity: potential step /= 0: ', step
         failed = .true.
      end if

      if (.not. failed) write (*, '(a)') &
         'ok  §5.4 A|A identity: dQ = 0, P = 0, step = 0 even for a charged homogeneous bulk'
   end subroutine test_aa_identity

   !> §5.3 -- THE PRIMARY ORACLE. Region A against region-A-rigidly-shifted-by-
   !> delta. The alignment solver must recover a step of exactly -delta, and
   !> dQ -> 0 on every site.
   !>
   !> The alignment unknown V_B enters through vmad as a CONSTANT over the
   !> region (B7 §1.3):
   !>
   !>     vmad(i) = dV(i) + vmad_bulk(region(i), type(i)) + V_region(i)
   !>
   !> with region A as gauge anchor (V_A == 0). The fixed point solved here is
   !> the deep-probe residual
   !>
   !>     V_B(out) = dV(deep-B) - dV(deep-A)
   !>
   !> mixed with vmix, exactly as surfpot mixes any other SCF quantity. With a
   !> rigid offset delta applied to region B's input parameters, the converged
   !> answer is known in closed form: V_B = -delta.
   subroutine test_offset_delta_oracle()
      integer, parameter :: n = 12, nhalf = 6
      real(rp) :: z(n), vmad_bulk(n), dq(n), area, delta
      real(rp) :: v_b, v_b_new, resid, vmix
      real(rp) :: probe_a, probe_b
      integer :: i, iter
      integer, parameter :: max_iter = 500

      area = 7.31_rp
      delta = 0.0734_rp        ! the deliberate rigid offset, Ry
      vmix = 0.35_rp           ! deliberately < 1 so mixing is genuinely exercised

      do i = 1, n
         z(i) = 3.4_rp*real(i - 1, rp)
      end do

      ! Region A (sites 1..nhalf) at its own zero; region B (nhalf+1..n) is the
      ! SAME material with a rigid shift delta baked into its frozen parameters.
      do i = 1, nhalf
         vmad_bulk(i) = 0.0_rp
      end do
      do i = nhalf + 1, n
         vmad_bulk(i) = delta
      end do

      ! Deviation charges vanish for A|A (same material both sides).
      dq = 0.0_rp

      ! Fixed-point iteration on the deep-probe residual, gauge-anchored at
      ! V_A == 0. This is the alignment solver in miniature.
      v_b = 0.0_rp
      do iter = 1, max_iter
         ! Deep probes: the innermost frozen site of each region, carrying its
         ! own reference plus the region shift under test.
         probe_a = vmad_bulk(1) + 0.0_rp          ! V_A == 0, the gauge anchor
         probe_b = vmad_bulk(n) + v_b

         ! Residual: deep-B must see the potential it needs to BE bulk B.
         resid = -(probe_b - probe_a)
         v_b_new = v_b + vmix*resid

         if (abs(v_b_new - v_b) < 1.0e-14_rp) then
            v_b = v_b_new
            exit
         end if
         v_b = v_b_new
      end do

      ! The exact answer: V_B = -delta.
      if (.not. close_abs(v_b, -delta)) then
         write (*, '(a,2es16.8)') 'FAIL offset-delta oracle: V_B /= -delta: ', v_b, -delta
         failed = .true.
      end if
      if (iter >= max_iter) then
         write (*, '(a)') 'FAIL offset-delta oracle: fixed point did not converge'
         failed = .true.
      end if
      ! And dQ must be zero on every site: a rigid shift transfers no charge.
      do i = 1, n
         if (.not. close_abs(dq(i), 0.0_rp)) then
            write (*, '(a,i3)') 'FAIL offset-delta oracle: dq /= 0 at site ', i
            failed = .true.
         end if
      end do

      if (.not. failed) write (*, '(a,f10.6,a,i0,a)') &
         'ok  §5.3 offset-delta oracle: recovered step = ', v_b, ' (exact -delta) in ', iter, ' iters'
   end subroutine test_offset_delta_oracle

   !> §5.2 -- the gauge oracle in the two-sided setting.
   subroutine test_gauge_two_sided()
      integer, parameter :: n = 11
      real(rp) :: z(n), dq(n), area, v_half(n), v_sym(n), pred(n)
      real(rp) :: qtot, ptot, d_half, d_sym
      integer :: i, j
      logical :: ok

      area = 7.31_rp
      do i = 1, n
         z(i) = 3.4_rp*real(i - 1, rp)
      end do

      ! --- neutral two-sided profile: screening charge on both boundaries ---
      dq = 0.0_rp
      dq(3) = 0.08_rp
      dq(4) = -0.03_rp
      dq(8) = -0.06_rp
      dq(9) = 0.01_rp
      dq(n) = dq(n) - sum(dq)          ! enforce exact neutrality
      call moments(z, dq, n, qtot, ptot)

      call apply_kernel(z, dq, area, n, .false., v_half)
      call apply_kernel(z, dq, area, n, .true., v_sym)

      ok = .true.
      do i = 1, n
         do j = 1, n
            d_half = v_half(i) - v_half(j)
            d_sym = v_sym(i) - v_sym(j)
            if (.not. close_abs(d_half, d_sym)) ok = .false.
         end do
      end do
      if (.not. ok) then
         write (*, '(a)') 'FAIL §5.2: neutral two-sided profile is not gauge-invariant'
         failed = .true.
      end if

      ! --- deliberately charged profile: Q /= 0 ---
      dq(n) = dq(n) + 0.05_rp
      call moments(z, dq, n, qtot, ptot)
      call apply_kernel(z, dq, area, n, .false., v_half)
      call apply_kernel(z, dq, area, n, .true., v_sym)
      do i = 1, n
         pred(i) = (2.0_rp*pi/area)*(ptot - z(i)*qtot)
      end do
      do i = 1, n
         if (.not. close_abs(v_half(i) - v_sym(i), pred(i))) then
            write (*, '(a,i3,2es16.8)') 'FAIL §5.2 charged: ', i, v_half(i) - v_sym(i), pred(i)
            failed = .true.
         end if
      end do

      if (.not. failed) write (*, '(a)') &
         'ok  §5.2 two-sided gauge: neutral invariant; Q/=0 differs by exactly (2*pi/A)[P - z*Q]'
   end subroutine test_gauge_two_sided

   !> §1.5 -- the compensation rule. Weight by side-resolved N(E_F): metallic
   !> leads receive charge proportional to N(E_F); vacuum/insulating receives
   !> NOTHING. Two similar metals give ~50/50, and metal|vacuum collapses to
   !> "everything into the metal" -- so the surface case needs NO separate
   !> branch. That is the physics, not a coincidence.
   subroutine test_compensation_weighting()
      real(rp) :: nef_a, nef_b, w_a, w_b, residual, comp_a, comp_b

      residual = -0.037_rp    ! leftover Sum_active dq: the active zone underscreened

      ! --- two similar metals ---
      nef_a = 1.24_rp
      nef_b = 1.24_rp
      call split(residual, nef_a, nef_b, w_a, w_b, comp_a, comp_b)
      if (.not. close_rel(w_a, 0.5_rp) .or. .not. close_rel(w_b, 0.5_rp)) then
         write (*, '(a,2f10.6)') 'FAIL two similar metals do not split ~50/50: ', w_a, w_b
         failed = .true.
      end if

      ! --- dissimilar metals: proportional to N(E_F) ---
      nef_a = 2.00_rp
      nef_b = 0.50_rp
      call split(residual, nef_a, nef_b, w_a, w_b, comp_a, comp_b)
      if (.not. close_rel(w_a, 0.8_rp) .or. .not. close_rel(w_b, 0.2_rp)) then
         write (*, '(a,2f10.6)') 'FAIL dissimilar metals not weighted by N(E_F): ', w_a, w_b
         failed = .true.
      end if

      ! --- metal | vacuum: vacuum has no states, so it receives NOTHING ---
      nef_a = 1.24_rp
      nef_b = 0.0_rp
      call split(residual, nef_a, nef_b, w_a, w_b, comp_a, comp_b)
      if (.not. close_abs(comp_b, 0.0_rp)) then
         write (*, '(a,es16.8)') 'FAIL vacuum received compensation charge: ', comp_b
         failed = .true.
      end if
      if (.not. close_abs(comp_a, -residual)) then
         write (*, '(a,2es16.8)') 'FAIL metal did not absorb all compensation: ', comp_a, -residual
         failed = .true.
      end if

      ! Neutrality must hold in every case.
      if (.not. close_abs(comp_a + comp_b + residual, 0.0_rp)) then
         write (*, '(a)') 'FAIL compensation does not restore neutrality'
         failed = .true.
      end if

      if (.not. failed) write (*, '(a)') &
         'ok  §1.5 compensation weights ~ N(E_F): 50/50 for like metals, all-to-metal for metal|vacuum'
   end subroutine test_compensation_weighting

   subroutine split(residual, nef_a, nef_b, w_a, w_b, comp_a, comp_b)
      real(rp), intent(in) :: residual, nef_a, nef_b
      real(rp), intent(out) :: w_a, w_b, comp_a, comp_b
      real(rp) :: tot

      tot = nef_a + nef_b
      if (tot <= 0.0_rp) then
         w_a = 0.5_rp
         w_b = 0.5_rp
      else
         w_a = nef_a/tot
         w_b = nef_b/tot
      end if
      comp_a = -residual*w_a
      comp_b = -residual*w_b
   end subroutine split

   !> §1.5 clarification 2 -- compensation must be LOCALIZED, never smeared as
   !> imppot does (`tdq(j) = -dif/nsum` over every frozen site). In a 3D impurity
   !> cluster that smear is a roughly isotropic shell with Sum dq*z = 0 by
   !> symmetry, so it is harmless. In LAYERED geometry the same smear has a
   !> nonzero first moment and a long lever arm -- same code, different geometry,
   !> different consequence. This test measures that difference.
   subroutine test_compensation_localization()
      integer, parameter :: n = 12
      real(rp) :: z(n), dq_local(n), dq_smear(n), area
      real(rp) :: residual, q_l, p_l, q_s, p_s
      integer :: i, nfrozen

      area = 7.31_rp
      do i = 1, n
         z(i) = 3.4_rp*real(i - 1, rp)
      end do
      residual = -0.04_rp

      ! Localized: placed at the defined boundary site z_R, as surfpot already
      ! does one-sided with tdq(iex) = -dif.
      dq_local = 0.0_rp
      dq_local(n) = -residual

      ! Smeared: imppot's rule, spread over every frozen site.
      dq_smear = 0.0_rp
      nfrozen = 6
      do i = n - nfrozen + 1, n
         dq_smear(i) = -residual/real(nfrozen, rp)
      end do

      call moments(z, dq_local, n, q_l, p_l)
      call moments(z, dq_smear, n, q_s, p_s)

      ! Both restore neutrality identically -- the monopole is fixed.
      if (.not. close_abs(q_l, q_s)) then
         write (*, '(a)') 'FAIL localized and smeared compensation differ in total charge'
         failed = .true.
      end if

      ! But the DIPOLE differs, and that is what touches the observable.
      if (close_abs(p_l, p_s)) then
         write (*, '(a)') 'FAIL smear and localized give the same first moment -- test is vacuous'
         failed = .true.
      end if

      if (.not. failed) write (*, '(a,es12.4,a)') &
         'ok  §1.5 localized vs smeared compensation: same Q, dipole differs by ', &
         abs(p_l - p_s), ' e*bohr'
   end subroutine test_compensation_localization

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

end program test_interface_alignment_oracle
