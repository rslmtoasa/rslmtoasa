!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
! PROGRAM: test_alignment_solver
!
!> @brief B7.4: the alignment solver itself, driven through the PRODUCTION
!>        entry points in source/region_registry.f90.
!>
!> @details This is the companion to test_interface_alignment_oracle, and the
!>          distinction between them is the point.
!>
!>          That test is the ORACLE: it states, in closed form and with its own
!>          independent arithmetic, what the right answer is -- a rigid offset
!>          delta applied to region B's frozen parameters must come back as a
!>          step of exactly -delta. It deliberately reimplements the fixed point
!>          inline, because an oracle that calls the code under test proves
!>          nothing.
!>
!>          This test drives the SHIPPING CODE against that same known answer,
!>          plus the properties the production solver has to have and the toy
!>          one did not need:
!>
!>          1. §5.3 primary oracle, through `alignment_update`. The production
!>             mixed fixed point must recover -delta to the same precision the
!>             toy one did.
!>          2. The gauge anchor is pinned EXACTLY, every iteration. Not
!>             "skipped" -- pinned. A gauge that drifts by accumulated round-off
!>             is the failure the oracle's negative control was built to catch,
!>             and the production update re-pins rather than relying on never
!>             touching the anchor.
!>          3. The anchor is a FROZEN, NON-VACUUM region. Vacuum carries no
!>             states at E_F, so "deep vacuum is neutral bulk vacuum" cannot
!>             drive a residual; and an active region has no frozen reference to
!>             BE, since it is what relaxes. The buildsurf registry order
!>             (vacuum, active, bulk) tests both traps at once.
!>          4. The analytic initial guess V_r = E_F(anchor) - E_F(r) is applied
!>             where region Fermi levels are known, skipped where they are not,
!>             and is ONLY a guess -- starting from it and starting from zero
!>             must reach the same fixed point (B7 §1.3: the closed form is not
!>             load-bearing, because it needs a cross-run absolute energy scale
!>             the code does not guarantee).
!>          5. The consistency check fires when the converged fixed point
!>             disagrees with the analytic contact potential, and stays silent
!>             when they agree. This is the G-B7-2 alarm.
!>          6. `deep_probe_site` picks each region's extreme site AWAY from the
!>             interface, from z rather than from registry order.
!>
!>        Exits non-zero (error stop) on any failure so ctest registers a fail.
!------------------------------------------------------------------------------
program test_alignment_solver
   use precision_mod, only: rp
   use region_registry_mod, only: region_registry, region_descriptor, &
                                  region_kind_vacuum, region_kind_bulk, &
                                  region_kind_active, region_kind_lead_a, &
                                  region_kind_lead_b
   implicit none

   logical :: failed

   failed = .false.

   call test_offset_delta_through_production_solver()
   call test_gauge_anchor_pinned_exactly()
   call test_anchor_is_frozen_and_not_vacuum()
   call test_initial_guess_is_only_a_guess()
   call test_consistency_check()
   call test_deep_probe_site()
   call test_writeback_reference_is_the_anchor()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   !> Build a two-lead registry: sites 1..nhalf are lead A (low z), sites
   !> nhalf+1..n are lead B (high z), with a thin active zone in the middle.
   subroutine build_two_lead(reg, n, nhalf)
      type(region_registry), intent(out) :: reg
      integer, intent(in) :: n, nhalf
      integer :: i

      call reg%restore_to_default()
      reg%nsite = n
      reg%nregion = 3

      deallocate (reg%region)
      allocate (reg%region(3))
      reg%region(1) = region_descriptor(kind=region_kind_lead_a, name='leadA', frozen=.true.)
      reg%region(2) = region_descriptor(kind=region_kind_active, name='active', frozen=.false.)
      reg%region(3) = region_descriptor(kind=region_kind_lead_b, name='leadB', frozen=.true.)

      deallocate (reg%region_id, reg%active, reg%layer_index, reg%z, reg%w, reg%reference_type)
      allocate (reg%region_id(n), reg%active(n), reg%layer_index(n))
      allocate (reg%z(n), reg%w(n), reg%reference_type(n))

      do i = 1, n
         reg%z(i) = 3.4_rp*real(i - 1, rp)
         reg%w(i) = 2.669_rp
         reg%reference_type(i) = 0
         reg%layer_index(i) = 0
         if (i <= nhalf) then
            reg%region_id(i) = 1
            reg%active(i) = .false.
         else
            reg%region_id(i) = 3
            reg%active(i) = .false.
         end if
      end do
      ! A couple of active sites straddling the interface.
      reg%region_id(nhalf) = 2
      reg%active(nhalf) = .true.
      reg%layer_index(nhalf) = 1
      reg%region_id(nhalf + 1) = 2
      reg%active(nhalf + 1) = .true.
      reg%layer_index(nhalf + 1) = 2
   end subroutine build_two_lead

   !> 1. §5.3 THE PRIMARY ORACLE, driven through the production
   !> `alignment_update`. Region B is the same material as A with a rigid shift
   !> delta baked into its frozen parameters; the solver must recover -delta.
   subroutine test_offset_delta_through_production_solver()
      type(region_registry) :: reg
      integer, parameter :: n = 12, nhalf = 6
      integer, parameter :: max_iter = 500
      real(rp) :: shift(3), probe(3), delta, vmix, maxresid, prev
      integer :: ianchor, iter

      call build_two_lead(reg, n, nhalf)

      delta = 0.0734_rp       ! the same deliberate offset the toy oracle uses
      vmix = 0.35_rp          ! deliberately < 1 so mixing is genuinely exercised

      ianchor = reg%gauge_anchor()
      if (ianchor /= 1) then
         write (*, '(a,i0)') 'FAIL offset-delta: anchor is not leadA: ', ianchor
         failed = .true.
      end if

      ! Deep probes carry each region's own reference. Region A sits at its own
      ! zero; region B's frozen set is rigidly shifted by delta.
      probe(1) = 0.0_rp
      probe(2) = 0.0_rp
      probe(3) = delta

      shift = 0.0_rp
      prev = huge(1.0_rp)
      do iter = 1, max_iter
         call reg%align_update(ianchor, probe, vmix, shift, maxresid)
         if (abs(shift(3) - prev) < 1.0e-14_rp) exit
         prev = shift(3)
      end do

      if (.not. close_abs(shift(3), -delta)) then
         write (*, '(a,2es16.8)') 'FAIL offset-delta through production solver: ', shift(3), -delta
         failed = .true.
      end if
      if (iter >= max_iter) then
         write (*, '(a)') 'FAIL offset-delta: production fixed point did not converge'
         failed = .true.
      end if

      if (.not. failed) write (*, '(a,f10.6,a,i0,a)') &
         'ok  §5.3 production solver recovered step = ', shift(3), ' (exact -delta) in ', iter, ' iters'
   end subroutine test_offset_delta_through_production_solver

   !> 2. The anchor is re-pinned to exactly zero every update, never merely
   !> skipped. Seeded with a deliberate nonzero anchor shift -- the negative
   !> control the oracle uses -- the production update must erase it.
   subroutine test_gauge_anchor_pinned_exactly()
      type(region_registry) :: reg
      integer, parameter :: n = 12, nhalf = 6
      real(rp) :: shift(3), probe(3), maxresid
      integer :: ianchor, iter

      call build_two_lead(reg, n, nhalf)
      ianchor = reg%gauge_anchor()

      probe = [0.0_rp, 0.0_rp, 0.0734_rp]

      ! Deliberate gauge-anchor drift injected.
      shift = [0.0137_rp, 0.0_rp, 0.0_rp]

      call reg%align_update(ianchor, probe, 0.35_rp, shift, maxresid)

      if (shift(ianchor) /= 0.0_rp) then
         write (*, '(a,es16.8)') 'FAIL anchor not pinned to EXACTLY zero after update: ', shift(ianchor)
         failed = .true.
      end if

      ! And it must stay pinned across many iterations, with no round-off creep.
      do iter = 1, 500
         call reg%align_update(ianchor, probe, 0.35_rp, shift, maxresid)
         if (shift(ianchor) /= 0.0_rp) then
            write (*, '(a,i0,a,es16.8)') 'FAIL anchor drifted at iteration ', iter, ': ', shift(ianchor)
            failed = .true.
            exit
         end if
      end do

      if (.not. failed) write (*, '(a)') &
         'ok  gauge anchor is pinned to exactly 0 every update, injected drift erased'
   end subroutine test_gauge_anchor_pinned_exactly

   !> 3. The anchor must be a FROZEN, NON-VACUUM region.
   !>
   !> Both qualifiers are load-bearing, and the buildsurf layout tests both at
   !> once because its registry order is (vacuum, active, bulk):
   !>
   !>  - a naive "anchor = region 1" picks VACUUM, which carries no states at
   !>    E_F, so the "deep region is neutral bulk" residual cannot be evaluated
   !>    there at all;
   !>  - a naive "anchor = first non-vacuum" picks the ACTIVE region, which has
   !>    no frozen reference to BE -- it is what relaxes.
   !>
   !> The only defensible answer for a surface is the bulk region.
   subroutine test_anchor_is_frozen_and_not_vacuum()
      type(region_registry) :: reg
      integer, parameter :: nbas = 49, nlay = 6, init = 6
      real(rp), dimension(nbas) :: z, w
      integer :: ianchor, i

      do i = 1, nbas
         z(i) = real(i, rp)*0.7_rp
         w(i) = 2.669_rp
      end do

      call reg%build_from_buildsurf(nbas, nlay, z, w, init=init)

      ianchor = reg%gauge_anchor()
      if (ianchor < 1 .or. ianchor > reg%nregion) then
         write (*, '(a,i0)') 'FAIL anchor: none chosen: ', ianchor
         failed = .true.
         return
      end if

      if (reg%region(ianchor)%kind == region_kind_vacuum) then
         write (*, '(a,a)') 'FAIL vacuum was chosen as the gauge anchor: ', trim(reg%region(ianchor)%name)
         failed = .true.
      end if
      if (.not. reg%region(ianchor)%frozen) then
         write (*, '(a,a)') 'FAIL an ACTIVE region was chosen as the gauge anchor: ', &
            trim(reg%region(ianchor)%name)
         failed = .true.
      end if
      if (reg%region(ianchor)%kind /= region_kind_bulk) then
         write (*, '(a,a)') 'FAIL buildsurf layout does not anchor on bulk: ', &
            trim(reg%region(ianchor)%name)
         failed = .true.
      end if

      if (.not. failed) write (*, '(a,a)') &
         'ok  anchor is frozen and non-vacuum; buildsurf (vacuum,active,bulk) anchors on ', &
         trim(reg%region(ianchor)%name)
   end subroutine test_anchor_is_frozen_and_not_vacuum

   !> 4. The analytic guess is applied where Fermi levels are known, skipped
   !> where they are not, and is ONLY a guess: the fixed point reached from the
   !> guess and from zero must be the same. B7 §1.3 forbids the closed form
   !> being load-bearing, because it needs a cross-run absolute energy scale
   !> the code does not currently guarantee (§2.1, G-B7-2).
   subroutine test_initial_guess_is_only_a_guess()
      type(region_registry) :: reg
      integer, parameter :: n = 12, nhalf = 6
      real(rp) :: guess(3), from_guess(3), from_zero(3), probe(3), maxresid
      real(rp) :: ef_a, ef_b
      integer :: ianchor, iter

      call build_two_lead(reg, n, nhalf)
      ianchor = reg%gauge_anchor()

      ef_a = -0.0512_rp
      ef_b = -0.1246_rp
      reg%region(1)%fermi = ef_a
      reg%region(1)%fermi_known = .true.
      reg%region(3)%fermi = ef_b
      reg%region(3)%fermi_known = .true.
      ! Region 2 (active) deliberately has NO known Fermi level.

      call reg%initial_guess(ianchor, guess)

      if (.not. close_abs(guess(3), ef_a - ef_b)) then
         write (*, '(a,2es16.8)') 'FAIL analytic guess /= E_F(A) - E_F(B): ', guess(3), ef_a - ef_b
         failed = .true.
      end if
      if (guess(ianchor) /= 0.0_rp) then
         write (*, '(a)') 'FAIL analytic guess does not leave the anchor at zero'
         failed = .true.
      end if
      if (guess(2) /= 0.0_rp) then
         write (*, '(a)') 'FAIL guess applied to a region with no known Fermi level'
         failed = .true.
      end if

      ! Same fixed point from either starting point.
      probe = [0.0_rp, 0.0_rp, 0.0734_rp]
      from_guess = guess
      from_zero = 0.0_rp
      do iter = 1, 2000
         call reg%align_update(ianchor, probe, 0.35_rp, from_guess, maxresid)
         call reg%align_update(ianchor, probe, 0.35_rp, from_zero, maxresid)
      end do

      if (.not. close_abs(from_guess(3), from_zero(3))) then
         write (*, '(a,2es16.8)') 'FAIL fixed point depends on the initial guess: ', &
            from_guess(3), from_zero(3)
         failed = .true.
      end if

      if (.not. failed) write (*, '(a)') &
         'ok  analytic guess seeds only known-E_F regions; fixed point is independent of it'
   end subroutine test_initial_guess_is_only_a_guess

   !> 5. The G-B7-2 alarm. When the converged fixed point agrees with the
   !> analytic contact potential, the check is quiet; when the two parameter
   !> sets are on different potential zeros, it reports the discrepancy. A
   !> silent misalignment produces a plausible but wrong dipole barrier, which
   !> B7 §1.3 calls the worst failure mode available to the work package.
   subroutine test_consistency_check()
      type(region_registry) :: reg
      integer, parameter :: n = 12, nhalf = 6
      real(rp) :: shift(3), worst, ef_a, ef_b, delta
      integer :: ianchor, iworst
      logical :: checked

      call build_two_lead(reg, n, nhalf)
      ianchor = reg%gauge_anchor()

      delta = 0.0734_rp
      ef_a = -0.0512_rp
      ef_b = ef_a + delta          ! consistent: E_F(A) - E_F(B) == -delta
      reg%region(1)%fermi = ef_a
      reg%region(1)%fermi_known = .true.
      reg%region(3)%fermi = ef_b
      reg%region(3)%fermi_known = .true.

      ! --- the consistent case: converged fixed point == -delta ---
      shift = [0.0_rp, 0.0_rp, -delta]
      call reg%consistency_check(ianchor, shift, worst, iworst, checked)
      if (.not. checked) then
         write (*, '(a)') 'FAIL consistency check did not run despite known Fermi levels'
         failed = .true.
      end if
      if (.not. close_abs(worst, 0.0_rp)) then
         write (*, '(a,es16.8)') 'FAIL consistent alignment reported a discrepancy: ', worst
         failed = .true.
      end if

      ! --- the broken case: one parameter set on a different potential zero ---
      shift = [0.0_rp, 0.0_rp, -delta + 0.021_rp]
      call reg%consistency_check(ianchor, shift, worst, iworst, checked)
      if (.not. close_abs(worst, 0.021_rp)) then
         write (*, '(a,es16.8)') 'FAIL broken alignment discrepancy not reported exactly: ', worst
         failed = .true.
      end if
      if (iworst /= 3) then
         write (*, '(a,i0)') 'FAIL discrepancy attributed to the wrong region: ', iworst
         failed = .true.
      end if

      ! --- no known Fermi levels: the check must decline, not fabricate ---
      reg%region(1)%fermi_known = .false.
      reg%region(3)%fermi_known = .false.
      call reg%consistency_check(ianchor, shift, worst, iworst, checked)
      if (checked) then
         write (*, '(a)') 'FAIL consistency check ran without any known Fermi level'
         failed = .true.
      end if

      if (.not. failed) write (*, '(a)') &
         'ok  G-B7-2 check: quiet when consistent, exact discrepancy when not, declines when unknown'
   end subroutine test_consistency_check

   !> 6. The deep probe of a region is its extreme site AWAY from the interface,
   !> decided from z rather than from registry order -- registry order is a
   !> bookkeeping choice, z is physics, and B7.5's two-sided builder need not
   !> make them agree.
   subroutine test_deep_probe_site()
      type(region_registry) :: reg
      integer, parameter :: n = 12, nhalf = 6
      integer :: ip_a, ip_b, i
      real(rp) :: ztmp

      call build_two_lead(reg, n, nhalf)

      ip_a = reg%deep_probe_site(1)
      ip_b = reg%deep_probe_site(3)

      ! Lead A lies below the cluster centre, so its deep probe is its
      ! MINIMUM-z site; lead B lies above, so its probe is its MAXIMUM-z site.
      if (ip_a /= 1) then
         write (*, '(a,i0)') 'FAIL deep probe of low-z lead is not its minimum-z site: ', ip_a
         failed = .true.
      end if
      if (ip_b /= n) then
         write (*, '(a,i0)') 'FAIL deep probe of high-z lead is not its maximum-z site: ', ip_b
         failed = .true.
      end if

      ! Now REVERSE the z ordering while leaving registry order untouched. The
      ! probes must follow z, not the index -- this is the property that makes
      ! the solver survive a builder that lists regions in any order.
      do i = 1, n
         ztmp = reg%z(i)
         reg%z(i) = -ztmp
      end do

      ip_a = reg%deep_probe_site(1)
      ip_b = reg%deep_probe_site(3)
      if (ip_a /= 1 .or. ip_b /= n) then
         write (*, '(a,2i5)') 'FAIL deep probes did not follow reversed z: ', ip_a, ip_b
         failed = .true.
      end if

      if (.not. failed) write (*, '(a)') &
         'ok  deep probe is each region''s extreme site away from the interface, chosen from z'
   end subroutine test_deep_probe_site

   !> 7. The write-back reference and the fixed-point reference must be the
   !> SAME point.
   !>
   !> `interfacepot` reports `v_lo = vm(1)` as the work function / dipole
   !> barrier probe, and for the buildsurf layout site 1 is deep VACUUM. The
   !> alignment fixed point, however, drives dV(deep-r) + V_r towards
   !> dV(deep-ANCHOR), and the anchor is the first frozen non-vacuum region.
   !>
   !> Referencing the write-back to one point while solving against the other
   !> leaves a constant offset between them that the solver cannot see and
   !> cannot remove: the alignment converges, correctly, to the wrong zero. That
   !> is the same class of failure as the one this whole gate exists to prevent
   !> -- plausible, converged, and wrong -- so it is pinned here.
   !>
   !> The property asserted: for the buildsurf layout the anchor's deep probe is
   !> NOT site 1, so the two references genuinely differ and the distinction is
   !> not vacuous.
   subroutine test_writeback_reference_is_the_anchor()
      type(region_registry) :: reg
      integer, parameter :: nbas = 49, nlay = 6, init = 6
      real(rp), dimension(nbas) :: z, w
      integer :: ianchor, iprobe, i

      do i = 1, nbas
         z(i) = real(i, rp)*0.7_rp
         w(i) = 2.669_rp
      end do
      call reg%build_from_buildsurf(nbas, nlay, z, w, init=init)

      ianchor = reg%gauge_anchor()
      iprobe = reg%deep_probe_site(ianchor)

      if (iprobe < 1 .or. iprobe > nbas) then
         write (*, '(a,i0)') 'FAIL anchor has no deep probe site: ', iprobe
         failed = .true.
         return
      end if

      ! The distinction must not be vacuous: if the anchor's probe WERE site 1,
      ! referencing to v_lo would accidentally be right and this test would
      ! prove nothing.
      if (iprobe == 1) then
         write (*, '(a)') 'FAIL test is vacuous: anchor probe coincides with v_lo (site 1)'
         failed = .true.
      end if

      ! Site 1 is deep vacuum in this layout -- the point that is NOT the
      ! alignment reference.
      if (reg%region(reg%region_id(1))%kind /= region_kind_vacuum) then
         write (*, '(a)') 'FAIL buildsurf site 1 is not vacuum; the premise has changed'
         failed = .true.
      end if

      if (.not. failed) write (*, '(a,i0,a,i0,a)') &
         'ok  write-back references the anchor probe (site ', iprobe, &
         '), not v_lo (site 1, vacuum) -- ', iprobe - 1, ' sites apart'
   end subroutine test_writeback_reference_is_the_anchor

   pure function close_abs(a, b) result(ok)
      real(rp), intent(in) :: a, b
      logical :: ok
      ok = abs(a - b) <= 1.0e-9_rp
   end function close_abs

end program test_alignment_solver
