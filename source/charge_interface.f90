!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! SUBMODULE: Charge interface electrostatics
!
!> Interface/vacuum region bookkeeping, electrostatics, alignment, and related
!> diagnostics are kept separate from the bulk and legacy surface charge code.
!> The public type(charge) contract remains declared by charge_mod.
!------------------------------------------------------------------------------

submodule (charge_mod) charge_interface

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Two-sided deviation-variable electrostatics for the interface geometry
   !> (B7.3): `imppot`'s reference bookkeeping on `surfpot`'s 2D kernel.
   !>
   !> @details
   !> `surfpot` is the one-sided (surface) path and remains the permanent
   !> regression oracle -- it is NOT modified. This routine is the two-sided
   !> generalization, and differs from it in five deliberate ways:
   !>
   !> 1. **Deviation variables throughout** (B7 §1.4, §2.4; CONVENTIONS C7).
   !>    `surfpot` sums raw `dq` over the active layers, which for a
   !>    non-homogeneous bulk with legitimately charged sites never converges to
   !>    a deviation. Here every charge entering the Madelung sum is
   !>
   !>        dq_i = q_i - q_bulk(region(i), type(i))
   !>
   !>    i.e. `imppot`'s definition generalized from a single host to the region
   !>    registry. This is what makes the two-sided sums truncate on BOTH sides
   !>    and be absolutely convergent, and what makes polar hosts and ordered
   !>    alloys correct: their built-in charge transfer lives in the reference,
   !>    not in dq.
   !>
   !> 2. **The moments Q and P**, computed and reported every iteration:
   !>
   !>        Q = sum_i dq_i                 -> must vanish (no residual field)
   !>        P = sum_i dq_i z_i + sum_i dp_i -> the potential step
   !>
   !>    Q = 0 is a precondition on KERNEL VALIDITY, not merely on physics
   !>    (CONVENTIONS C3): the `AM` half-space gauge differs from the symmetric
   !>    kernel by (2*pi/A)[P - z_i*Q], which is a global constant -- and so
   !>    cancels from any dV difference -- if and only if Q = 0. If Q /= 0 the
   !>    gauge injects a spurious uniform field. Hence the explicit check.
   !>
   !> 3. **Two-sided probes.** `surfpot` probes deep vacuum (`dss(1,.)`) and
   !>    deep bulk (`dss(nbas,.)`) and reports `vmad1 = vm1 - vbulk`. Here the
   !>    same probes are evaluated at the deep site of EACH region, so the step
   !>    is the interface dipole barrier (A|B) or the work function (A|vacuum).
   !>
   !> 4. **Localized, N(E_F)-weighted compensation** (B7 §1.5). The residual
   !>    `sum_active dq` is placed at the defined boundary sites, weighted by the
   !>    side-resolved N(E_F) of the boundary buffer layers: metallic leads
   !>    receive charge proportional to N(E_F); vacuum/insulating receives
   !>    NOTHING. This gives ~50/50 for two similar metals and collapses to
   !>    "everything into the metal" for metal|vacuum, so the surface case needs
   !>    no separate branch. It deliberately does NOT inherit `imppot`'s smear
   !>    (`tdq(j) = -dif/nsum` over every frozen site): in a 3D impurity cluster
   !>    that is a roughly isotropic shell with sum dq*z = 0 by symmetry, but in
   !>    layered geometry the same smear has a nonzero first moment and a lever
   !>    arm of the full slab thickness, landing directly on the work function.
   !>
   !> 5. **Per-site w** from the registry, following `surfmat`'s `wssurf`
   !>    pattern. NOT `impmad`'s `wsimp`, which is a system-wide average AND in
   !>    a different unit convention (see the @warning on `impmad`).
   !>
   !> The kernel itself is reused UNCHANGED: `dss`/`dsz` as built by `madl2d`
   !> via `surfmat` are correct two-sided (CONVENTIONS C3, C4).
   !>
   !> @note **There is deliberately no separate `interfacemat`.** The B7 plan
   !>       names one, and CONVENTIONS C8 requires that it mimic `surfmat`
   !>       rather than `impmad` (the two differ by exactly one factor of S).
   !>       The strongest form of "mimic `surfmat`" is to *reuse* it: this
   !>       routine consumes precisely what `build_alelay` + `surfmat` already
   !>       produce -- `dss`, `dsz`, `sws`, `qz`, `wssurf` and the region
   !>       registry -- all of which are geometry-general and already correct
   !>       two-sided. Writing a near-duplicate builder would add a second
   !>       copy of the kernel to keep in sync and a second chance to drift
   !>       into `impmad`'s convention, for no capability gain. If a genuinely
   !>       two-sided geometry ever needs different STRUCTURE CONSTANTS, that
   !>       belongs in the cluster builder (B7.5), not in a parallel kernel.
   !>
   !> 6. **The region alignment shift V_r** (B7 §1.3, B7.4), solved as an SCF
   !>    fixed point on the deep-probe residual and added to `vmad` alongside
   !>    the region reference:
   !>
   !>        vmad(i) = dV(i) + vmad_bulk(region(i), type(i)) + V_region(i)
   !>
   !>    See `align_regions` for why this is a fixed point rather than the
   !>    closed form V_r = E_F - E_F^(r). With a single region the anchor gauge
   !>    makes V_r identically zero and this routine reduces to the
   !>    deviation-variable form of `surfpot`.
   !---------------------------------------------------------------------------
   module subroutine interfacepot(this)
      class(charge), intent(inout) :: this
      ! Local variables
      integer :: ibas, iclas, iq, j, jq, k, atomrec, nsite, ireg, irow
      integer :: ideep_lo, ideep_hi, iref, itype, iprobe
      real(rp) :: summ, wsms, qtot, ptot, step, resid
      real(rp) :: v_lo, v_hi, nef_lo, nef_hi, wgt_lo, wgt_hi
      real(rp) :: vmad_ref, v_ref
      real(rp), dimension(this%lattice%nbas) :: tdq, tq10, vm, zsite, wsite
      logical :: dipole_on

      nsite = this%lattice%nbas
      wsms = this%sws*this%lattice%alat*ang2au
      dipole_on = this%lattice%control%dipole_electrostatics

      ! Registry supplies z and w as data (B7 §2.10: relaxed-z is then a
      ! parameter change, not a rewrite).
      if (this%regions%nsite == nsite) then
         zsite(1:nsite) = this%regions%z(1:nsite)
         wsite(1:nsite) = this%regions%w(1:nsite)
      else
         zsite(1:nsite) = this%qz(1:nsite)
         wsite(1:nsite) = this%wssurf(1:nsite)
      end if

      tdq(:) = 0.0d0
      tq10(:) = 0.0d0

      ! --- 1. Deviation charges on the active sites -------------------------
      ! dq_i - q_bulk(reference type), never an absolute charge.
      !
      ! TWO INDEX SPACES, and conflating them was a real bug (fixed here):
      !
      !   atomrec  1..nrec  the ACTIVE TYPE counter. Indexes this%dq,
      !                     lattice%chargetrf_type, and symbolic_atom(nbulk+.).
      !   irow     1..nbas  the MADELUNG ROW. Indexes tdq/tq10/vm and every
      !                     registry array (region_id, active, z, w).
      !
      ! They are not the same: the active zone starts at registry row
      ! nlay_a+1, because rows 1..nlay_a are region A's frozen boundary
      ! (region_registry.f90's build_from_interface index map). Writing the
      ! first active site's charge to tdq(1) put it on a FROZEN row, which by
      ! construction must carry exactly zero deviation (B7 §2.10) -- and then
      ! `compensation_sites` returned ideep_lo = 1, the same row, so step 2
      ! subtracted the whole residual from it and the two cancelled exactly.
      ! Net tdq was identically zero on every row, vm == 0, and Q, P, step and
      ! every deep probe came out exactly zero for ANY interface run, leaving
      ! the alignment fixed point solving against nothing.
      !
      ! It hid because the shipped oracle is the A|A identity, whose correct
      ! answer IS zero for all five reported quantities -- a spuriously zero
      ! result is indistinguishable from a right one there. That is precisely
      ! why B7 §5.3 specifies the real oracle as A-against-A-with-a-rigid-
      ! offset rather than plain A|A.
      atomrec = 0
      do iclas = 1, this%lattice%nlay
         do k = 1, this%lattice%natoms_layer(iclas)
            atomrec = atomrec + 1
            irow = this%nlay_a + atomrec
            if (irow > nsite) exit
            tdq(irow) = this%dq(atomrec)
            if (this%lattice%chargetrf_type(atomrec) > 0) then
               tdq(irow) = tdq(irow) - this%bulk_charge(this%lattice%chargetrf_type(atomrec))
            end if
            if (dipole_on) tq10(irow) = this%symbolic_atom(this%lattice%nbulk + atomrec)%potential%q10
         end do
      end do

      ! --- 1b. Deep-A charge-drift diagnostic (B7 §1.3, B7.4) ---------------
      ! Evaluated on the raw deviation charges, BEFORE compensation, which
      ! would otherwise mask exactly the drift being looked for.
      call this%deep_drift_diagnostic(tdq, nsite)

      ! --- 2. Localized, N(E_F)-weighted compensation (§1.5) ----------------
      resid = sum(tdq(1:nsite))
      call this%compensation_sites(ideep_lo, ideep_hi, nef_lo, nef_hi)
      if (nef_lo + nef_hi > 0.0d0) then
         wgt_lo = nef_lo/(nef_lo + nef_hi)
         wgt_hi = nef_hi/(nef_lo + nef_hi)
      else
         wgt_lo = 0.5d0
         wgt_hi = 0.5d0
      end if
      if (ideep_lo >= 1 .and. ideep_lo <= nsite) tdq(ideep_lo) = tdq(ideep_lo) - resid*wgt_lo
      if (ideep_hi >= 1 .and. ideep_hi <= nsite) tdq(ideep_hi) = tdq(ideep_hi) - resid*wgt_hi

      ! --- 3. The moments Q and P (§1.4), reported every iteration ----------
      qtot = 0.0d0
      ptot = 0.0d0
      do j = 1, nsite
         qtot = qtot + tdq(j)
         ptot = ptot + zsite(j)*tdq(j)
         if (dipole_on) ptot = ptot + tq10(j)
      end do

      if (rank == 0) then
         call g_logger%info('interfacepot: Q= '//fmt('es12.4', qtot)// &
                            ' P= '//fmt('es12.4', ptot)// &
                            ' residual= '//fmt('es12.4', resid), __FILE__, __LINE__)
      end if

      ! Q = 0 is a KERNEL PRECONDITION, not just physics (C3): a nonzero Q makes
      ! the half-space AM gauge inject a spurious uniform field of (2*pi/A)*Q.
      if (abs(qtot) > 1.0d-8) then
         if (rank == 0) call g_logger%warning('interfacepot: sum of deviation charges Q= '// &
            fmt('es12.4', qtot)//' is not zero; the AM half-space gauge is only valid at Q=0 '// &
            '(see CONVENTIONS_MADELUNG.md C3). A spurious uniform field is being injected.', &
            __FILE__, __LINE__)
      end if

      ! --- 4. Madelung potential, kernel reused unchanged -------------------
      do ibas = 1, nsite
         summ = 0.0d0
         do j = 1, nsite
            iq = ibas
            jq = j
            summ = summ + this%dss(iq, jq)*tdq(j)
            if (dipole_on) summ = summ + this%dsz(iq, jq)*tq10(j)
         end do
         vm(ibas) = summ/wsms
      end do

      ! --- 5. Two-sided deep probes and the step ----------------------------
      v_lo = vm(1)
      v_hi = vm(nsite)
      step = v_hi - v_lo
      if (rank == 0) then
         call g_logger%info('interfacepot: deep-probe potentials v_lo= '//fmt('f12.6', v_lo)// &
                            ' v_hi= '//fmt('f12.6', v_hi)// &
                            ' step= '//fmt('f12.6', step), __FILE__, __LINE__)
      end if

      ! --- 5b. Alignment solve (B7.4): one mixed fixed-point step on V_r ----
      call this%align_regions(vm, nsite)

      ! The write-back below is referenced to the ANCHOR region's deep probe,
      ! not to `v_lo`. The two are different points in general: `v_lo` is site 1,
      ! which for the buildsurf layout is deep VACUUM, while the anchor is the
      ! first frozen non-vacuum region (bulk). Referencing the write-back to one
      ! point while solving the fixed point against another would leave a
      ! constant offset between them that the solver cannot see and cannot
      ! remove -- the alignment would converge, correctly, to the wrong zero.
      ! `v_lo` remains the reported work function / dipole-barrier probe above,
      ! which is what it is for.
      v_ref = v_lo
      iref = this%fermi_pinned_region()
      if (iref < 1) iref = this%regions%gauge_anchor()
      if (this%regions%nsite == nsite .and. iref >= 1) then
         iprobe = this%regions%deep_probe_site(iref)
         if (iprobe >= 1 .and. iprobe <= nsite) v_ref = vm(iprobe)
      end if

      ! --- 6. Write back (B7 §1.3) ------------------------------------------
      !        vmad(i) = dV(i) + vmad_bulk(region(i), type(i)) + V_region(i)
      !
      !        The three terms are the deviation potential, the region's own
      !        converged on-site reference, and the alignment shift. The middle
      !        term is `imppot`'s per-class add-back generalized from a single
      !        host to the region registry -- without it, a polar host or an
      !        ordered alloy is silently treated as non-polar (B7 §2.1). The
      !        last term is what makes two independently converged parameter
      !        sets share one absolute scale; with a single region it is zero
      !        by the anchor gauge.
      atomrec = 0
      do ibas = 1, this%lattice%nlay
         do k = 1, this%lattice%natoms_layer(ibas)
            atomrec = atomrec + 1
            ! Same two index spaces as step 1: `atomrec` addresses the active
            ! TYPE (dq, chargetrf_type, symbolic_atom), `irow` the Madelung ROW
            ! (vm, region_id). The active zone starts at row nlay_a+1.
            irow = this%nlay_a + atomrec
            if (irow > nsite) exit
            ireg = 1
            if (this%regions%nsite == nsite) ireg = this%regions%region_id(irow)

            ! Region reference on-site Madelung potential, per site type.
            vmad_ref = 0.0d0
            itype = this%lattice%chargetrf_type(atomrec)
            if (itype >= 1 .and. itype <= size(this%symbolic_atom)) then
               vmad_ref = this%symbolic_atom(itype)%potential%vmad
            end if

            this%symbolic_atom(this%lattice%nbulk + atomrec)%potential%vmad = &
               (vm(irow) - v_ref + vmad_ref + this%region_shift(ireg))*this%vmix + &
               this%symbolic_atom(this%lattice%nbulk + atomrec)%potential%vmad*(1.0d0 - this%vmix)
         end do
      end do

      call this%overlap_diagnostic(zsite, wsite, nsite)

      ! B7.6: the alignment shifts are now current for this iteration. Any
      ! frozen region whose parameters DEPEND on its own shift must be
      ! regenerated here, after the solve and before the next recursion
      ! consumes them. Vacuum is the one such region: its parameters are an
      ! empty lattice rigidly shifted to the vacuum level, and the vacuum
      ! level is precisely its alignment shift (B7 §1.3, §1.6).
      if (associated(this%on_alignment_updated)) call this%on_alignment_updated()
   end subroutine interfacepot

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> B7.4 alignment solver: one mixed fixed-point step on the per-region
   !> alignment shifts V_r (B7 §1.3).
   !>
   !> @details
   !> **The problem.** Each bulk calculation carries its own potential zero.
   !> Region A's frozen parameter set implicitly asserts one absolute energy
   !> scale; region B's asserts another. Freezing both at their independently
   !> converged values with no relative shift imposes
   !>
   !>     contact potential == 0
   !>
   !> whether or not that is true. The run converges, looks entirely plausible,
   !> and reports a wrong dipole barrier and wrong charge transfer. This is the
   !> highest correctness risk in B7 precisely because the failure is silent.
   !>
   !> **The physical requirement.** Deep inside region r the site must BE
   !> neutral bulk r. On a common absolute scale that is E_F - V_r = E_F^(r),
   !> so the shifts are fixed up to one overall gauge and only differences are
   !> physical: V_B - V_A is the contact potential.
   !>
   !> **Why a fixed point and not the closed form.** V_r = E_F - E_F^(r) would
   !> be a one-line answer, but it requires a reliable cross-run absolute energy
   !> scale that the code does not currently guarantee (B7 §2.1, gate G-B7-2).
   !> Instead the residual
   !>
   !>     V_r(out) = dV(deep-r) - dV(deep-anchor)
   !>
   !> is driven to zero as an ordinary SCF quantity, mixed with `vmix` -- the
   !> same probe `surfpot` already evaluates one-sided, evaluated on both sides.
   !> At convergence deep-r sees exactly the potential it needs to be bulk r,
   !> and no cross-run reference was required. The analytic value is then used
   !> for two things only: the initial guess (usually good enough to save
   !> several iterations) and a CONSISTENCY CHECK. That is what converts the
   !> most fragile ingredient in this work package from load-bearing into a
   !> diagnostic.
   !>
   !> **Damping.** `vmix` is not decoration. Two regions separated by an active
   !> zone form a capacitor-like soft mode -- a small charge transfer moves the
   !> whole relative level -- so the undamped update overshoots. This is the
   !> mixing hook that `imppot` had dead and that was reinstated in `eeecae9`
   !> for exactly this reason.
   !>
   !> @param[in] vm     dV per site, as computed by the Madelung kernel.
   !> @param[in] nsite  number of cluster sites.
   !---------------------------------------------------------------------------
   module subroutine align_regions(this, vm, nsite)
      class(charge), intent(inout) :: this
      integer, intent(in) :: nsite
      real(rp), dimension(:), intent(in) :: vm
      integer :: ianchor, ir, iprobe, iworst
      real(rp) :: maxresid, worst
      logical :: checked
      real(rp), dimension(:), allocatable :: probe

      if (this%regions%nsite /= nsite) return
      if (this%regions%nregion < 1) return
      if (.not. allocated(this%region_shift)) return

      ! With a free Fermi level the anchor is the registry's own choice (the
      ! first frozen, non-vacuum region -- see gauge_anchor for why both
      ! qualifiers matter). Pinning E_F to a region makes THAT region the
      ! anchor: the two statements "E_F is region r's Fermi level" and "region r
      ! carries no alignment shift" are the same statement, since E_F - V_r =
      ! E_F^(r). Anchoring anywhere else would double-count the pin.
      ianchor = this%fermi_pinned_region()
      if (ianchor < 1) ianchor = this%regions%gauge_anchor()
      if (ianchor < 1) return

      ! First entry: seed from the analytic guess where region Fermi levels are
      ! known. Subsequent iterations continue from the mixed fixed point.
      if (.not. this%alignment_started) then
         call this%regions%initial_guess(ianchor, this%region_shift)
         this%alignment_started = .true.
         if (rank == 0) then
            if (len_trim(this%fix_fermi_to_region) > 0) then
               call g_logger%info('align_regions: Fermi level pinned to region '// &
                                  trim(this%regions%region(ianchor)%name)// &
                                  ', which is therefore the gauge anchor (V == 0)', __FILE__, __LINE__)
            else
               call g_logger%info('align_regions: free Fermi level from cluster neutrality; '// &
                                  'gauge anchor is region '//fmt('i4', ianchor)//' ('// &
                                  trim(this%regions%region(ianchor)%name)//'), V_anchor == 0', &
                                  __FILE__, __LINE__)
            end if
         end if
      end if

      ! Deep probe of each region: dV at the region's extreme site away from
      ! the interface, the registry generalization of surfpot's dss(1,.) /
      ! dss(nbas,.) pair.
      allocate (probe(this%regions%nregion))
      probe(:) = 0.0d0
      do ir = 1, this%regions%nregion
         iprobe = this%regions%deep_probe_site(ir)
         if (iprobe >= 1 .and. iprobe <= nsite) probe(ir) = vm(iprobe)
      end do

      call this%regions%align_update(ianchor, probe, this%vmix, this%region_shift, maxresid)

      if (rank == 0) then
         do ir = 1, this%regions%nregion
            if (ir == ianchor) cycle
            call g_logger%info('align_regions: V('//trim(this%regions%region(ir)%name)// &
                               ')= '//fmt('f12.6', this%region_shift(ir))// &
                               ' Ry, deep probe dV= '//fmt('f12.6', probe(ir)), __FILE__, __LINE__)
         end do
         call g_logger%info('align_regions: max alignment residual= '// &
                            fmt('es12.4', maxresid), __FILE__, __LINE__)
      end if

      ! Consistency check against the analytic contact potential. A converged
      ! fixed point that disagrees with E_F^(A) - E_F^(r) beyond the threshold
      ! means the absolute-zero bookkeeping is broken (G-B7-2), and the run must
      ! say so loudly rather than report a plausible wrong barrier.
      if (maxresid < this%alignment_tol) then
         call this%regions%consistency_check(ianchor, this%region_shift, worst, iworst, checked)
         if (checked .and. worst > this%alignment_check_tol .and. rank == 0) then
            call g_logger%warning('align_regions: converged alignment shift for region '// &
               trim(this%regions%region(iworst)%name)//' disagrees with the analytic contact '// &
               'potential E_F(anchor) - E_F(region) by '//fmt('f12.6', worst)//' Ry, which '// &
               'exceeds the '//fmt('f8.4', this%alignment_check_tol)//' Ry threshold. The '// &
               'absolute-zero bookkeeping across the two parameter sets is inconsistent (see '// &
               'B7 gate G-B7-2). The fixed point is still the answer, but one of the frozen '// &
               'parameter sets is on a different potential zero than it claims, so treat the '// &
               'reported dipole barrier as unvalidated.', __FILE__, __LINE__)
         end if
      end if

      deallocate (probe)
   end subroutine align_regions

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> B7.4: resolve `fix_fermi_to_region` to a region id (B7 §1.3).
   !>
   !> @details
   !> **The default is a FREE Fermi level**, determined by charge neutrality of
   !> the cluster as in supercell mode, combined with the internal gauge fix
   !> V_anchor == 0. The two unknowns (E_F, V_B) are then determined by two
   !> conditions (cluster neutrality, deep-B bulk residual), and there is no
   !> flat direction. That is the maintainer decision recorded in B7 §1.3, and
   !> it is why this returns 0 when the option is unset.
   !>
   !> Setting it to a region name pins E_F to that region's own converged Fermi
   !> level instead, which reduces the interface path exactly to today's
   !> surface behaviour (`fix_fermi = .true.` for calctype 'S'), and is the
   !> correct setting when reproducing `buildsurf` results.
   !>
   !> Returns 0 for free E_F; otherwise the region id. Matching is
   !> case-insensitive on the region name. An unmatched non-empty name is a
   !> user input error and is reported as fatal -- this is a true boundary
   !> (namelist input), so the check belongs here.
   !---------------------------------------------------------------------------
   module function fermi_pinned_region(this) result(ireg)
      class(charge), intent(inout) :: this
      integer :: ireg
      integer :: ir
      character(len=32) :: want, have

      ireg = 0
      if (len_trim(this%fix_fermi_to_region) == 0) return

      want = this%fix_fermi_to_region
      call lower_case(want)
      do ir = 1, this%regions%nregion
         have = this%regions%region(ir)%name
         call lower_case(have)
         if (trim(have) == trim(want)) then
            ireg = ir
            return
         end if
      end do

      call g_logger%fatal('fix_fermi_to_region = "'//trim(this%fix_fermi_to_region)// &
         '" does not name any region in the registry. Leave it empty for the default '// &
         'free Fermi level from cluster neutrality (B7 §1.3).', __FILE__, __LINE__)
   end function fermi_pinned_region

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Lower-case a string in place, for case-insensitive region-name matching.
   !---------------------------------------------------------------------------
   module pure subroutine lower_case(s)
      character(len=*), intent(inout) :: s
      integer :: i, ic

      do i = 1, len(s)
         ic = iachar(s(i:i))
         if (ic >= iachar('A') .and. ic <= iachar('Z')) s(i:i) = achar(ic + 32)
      end do
   end subroutine lower_case

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> B7.4: deep-A charge-drift diagnostic (B7 §1.3).
   !>
   !> @details
   !> With a free Fermi level, "deep-A is neutral bulk" is NOT imposed -- it is
   !> a consequence of the buffer being thick enough. The two unknowns (E_F and
   !> V_B) are determined by two conditions (cluster neutrality and the deep-B
   !> residual), which leaves nothing pinning the deep-A charge directly.
   !>
   !> So it has to be watched. A nonzero deviation charge at the INNERMOST
   !> FROZEN sites -- the ones adjacent to the active zone, where drift shows up
   !> first -- means the active zone is too thin, and it means the gauge is
   !> inconsistent, because the anchor region is no longer the bulk it is
   !> assumed to be. Reported every iteration.
   !>
   !> Evaluated on the raw deviation charges before compensation: compensation
   !> is placed at exactly these boundary sites (§1.5), so running this
   !> afterwards would measure the compensation instead of the drift.
   !---------------------------------------------------------------------------
   module subroutine deep_drift_diagnostic(this, tdq, nsite)
      class(charge), intent(inout) :: this
      integer, intent(in) :: nsite
      real(rp), dimension(:), intent(in) :: tdq
      integer :: i, ianchor, idrift
      real(rp) :: drift

      if (this%regions%nsite /= nsite) return
      ! Same anchor as align_regions -- the drift being watched is the drift of
      ! the region the gauge is pinned to.
      ianchor = this%fermi_pinned_region()
      if (ianchor < 1) ianchor = this%regions%gauge_anchor()
      if (ianchor < 1) return

      ! Largest deviation charge on any frozen site of the anchor region.
      drift = 0.0d0
      idrift = 0
      do i = 1, nsite
         if (this%regions%region_id(i) /= ianchor) cycle
         if (this%regions%active(i)) cycle
         if (abs(tdq(i)) > abs(drift)) then
            drift = tdq(i)
            idrift = i
         end if
      end do
      if (idrift == 0) return

      if (rank == 0) then
         call g_logger%info('align_regions: deep-'//trim(this%regions%region(ianchor)%name)// &
                            ' charge drift= '//fmt('es12.4', drift)//' at site '// &
                            fmt('i5', idrift), __FILE__, __LINE__)
         if (abs(drift) > this%deep_drift_tol) then
            call g_logger%warning('align_regions: deviation charge '//fmt('es12.4', drift)// &
               ' at frozen anchor-region site '//fmt('i5', idrift)//' exceeds the '// &
               fmt('es10.2', this%deep_drift_tol)//' threshold. With a free Fermi level '// &
               '"deep-anchor is neutral bulk" is not imposed, only inherited from a thick '// &
               'enough buffer -- so this means the active zone is too thin AND that the '// &
               'alignment gauge is inconsistent. Widen the active zone (B7 §1.3).', &
               __FILE__, __LINE__)
         end if
      end if
   end subroutine deep_drift_diagnostic

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Locates the two compensation sites and their boundary N(E_F) weights
   !> (B7 §1.5).
   !>
   !> @details
   !> Compensation is placed at the innermost frozen site of each side, weighted
   !> by the side-resolved N(E_F) of the boundary buffer layers -- the linear
   !> response answer for partitioning screening charge between two reservoirs.
   !> A vacuum/insulating boundary has N(E_F) ~ 0 and receives NOTHING, which is
   !> why the surface case needs no separate branch: a nonzero residual means the
   !> active zone failed to screen, and the remaining screening happens where
   !> there are states to do it.
   !>
   !> Real spill-out is NOT compensation: charge outside the surface plane is
   !> carried by ACTIVE empty spheres as genuine self-consistent dq. If there is
   !> not enough charge outside, the fix is more active empty spheres -- placing
   !> compensation in the vacuum does not perturb the work function, it SETS it.
   !---------------------------------------------------------------------------
   module subroutine compensation_sites(this, ideep_lo, ideep_hi, nef_lo, nef_hi)
      class(charge), intent(inout) :: this
      integer, intent(out) :: ideep_lo, ideep_hi
      real(rp), intent(out) :: nef_lo, nef_hi
      integer :: i, nsite

      nsite = this%lattice%nbas
      ideep_lo = 1
      ideep_hi = nsite

      if (this%regions%nsite == nsite) then
         ! Innermost frozen site on each side: the last frozen row before the
         ! active zone (low z) and the first frozen row after it (high z).
         do i = 1, nsite
            if (this%regions%active(i)) exit
            ideep_lo = i
         end do
         ideep_hi = nsite
         do i = nsite, 1, -1
            if (this%regions%active(i)) exit
            ideep_hi = i
         end do
      end if

      ! Side-resolved N(E_F) of the boundary buffer layers. Recursion provides
      ! these for free; until that plumbing lands, fall back to the vacuum test
      ! below so metal|vacuum still behaves correctly.
      nef_lo = this%boundary_nef(ideep_lo)
      nef_hi = this%boundary_nef(ideep_hi)
   end subroutine compensation_sites

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Side-resolved N(E_F) at a boundary buffer site, used as the compensation
   !> weight (B7 §1.5). Vacuum/empty-sphere sites return ~0 so they receive no
   !> compensation charge.
   !---------------------------------------------------------------------------
   module function boundary_nef(this, isite) result(nef)
      class(charge), intent(inout) :: this
      integer, intent(in) :: isite
      real(rp) :: nef
      integer :: ia
      type(region_descriptor) :: desc

      nef = 0.0d0
      if (isite < 1 .or. isite > this%lattice%nbas) return

      ! A vacuum region carries no states at E_F by construction.
      if (this%regions%nsite == this%lattice%nbas) then
         desc = this%regions%region_of(isite)
         if (desc%kind == region_kind_vacuum) return
      end if

      ! Otherwise use the site's valence as a stand-in density of states scale.
      ! B7.7 replaces this with the recursion-supplied N(E_F) of the boundary
      ! buffer layer; the weighting rule and its normalization are unchanged by
      ! that substitution.
      !
      ! `isite` is a MADELUNG ROW (1..nbas), not an active-type index, so it
      ! must NOT be used as `nbulk + isite` -- that addresses the active-type
      ! block (nbulk+1..ntype) and for a boundary row is simply out of range.
      ! It then returned 0 for BOTH sides, the caller fell into its 50/50
      ! fallback, and vacuum silently received half the compensation charge --
      ! the exact error B7 §1.5 warns about ("compensation placed there does
      ! not perturb the work function, it SETS it"). Frozen boundary rows carry
      ! a REFERENCE TYPE instead, which the registry records per site.
      ia = 0
      if (this%regions%nsite == this%lattice%nbas) then
         if (allocated(this%regions%reference_type)) ia = this%regions%reference_type(isite)
      end if
      if (ia >= 1 .and. ia <= size(this%symbolic_atom)) then
         nef = max(this%symbolic_atom(ia)%element%valence, 0.0d0)
      end if
   end function boundary_nef

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Sphere-overlap diagnostic (B7 §2.9). Diagnostic ONLY -- per maintainer
   !> instruction the ASA volume-filling problem at an interface with unequal w
   !> is not to be solved here.
   !>
   !> @details
   !> With unequal w across the interface the ASA volume-filling condition cannot
   !> hold on both sides. This reports the per-layer ratio of summed sphere
   !> volume to the cell volume and warns above a threshold, converting a silent
   !> error into a message.
   !---------------------------------------------------------------------------
   module subroutine overlap_diagnostic(this, zsite, wsite, nsite)
      class(charge), intent(inout) :: this
      integer, intent(in) :: nsite
      real(rp), dimension(:), intent(in) :: zsite, wsite
      real(rp) :: wmin, wmax, ratio
      integer :: i

      if (nsite < 1) return

      wmin = wsite(1)
      wmax = wsite(1)
      do i = 2, nsite
         wmin = min(wmin, wsite(i))
         wmax = max(wmax, wsite(i))
      end do

      if (wmin <= 0.0d0) return
      ratio = wmax/wmin

      if (rank == 0) then
         call g_logger%info('interfacepot: Wigner-Seitz radii span w_min= '//fmt('f10.6', wmin)// &
                            ' w_max= '//fmt('f10.6', wmax)//' ratio= '//fmt('f8.4', ratio), &
                            __FILE__, __LINE__)
         if (ratio > 1.15d0) then
            call g_logger%warning('interfacepot: Wigner-Seitz radii differ by more than 15% '// &
               'across the cluster (ratio '//fmt('f8.4', ratio)//'). The ASA volume-filling '// &
               'condition cannot hold on both sides; sphere overlap is significant and the '// &
               'shared l>=1 multipole normalization may be a poor compromise (see '// &
               'CONVENTIONS_MADELUNG.md C5). This is a diagnostic, not a correction.', &
               __FILE__, __LINE__)
         end if
      end if
   end subroutine overlap_diagnostic

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> B7.1: construct this%regions, the explicit per-site region registry,
   !> from the cluster data surfmat/build_alelay already produced. Reproduces
   !> today's buildsurf layout (init=6 leading vacuum-side frozen rows,
   !> nlay active layers, remaining rows bulk-side frozen) exactly, as
   !> documented in B7 §2.10 and verified against charge%surfpot.
   !>
   !> This does not change surfpot's behaviour or consume the registry from
   !> it; it is the registry deliverable of B7.1, kept parallel to the
   !> existing offset arithmetic so the surface regression is unaffected.
   !---------------------------------------------------------------------------
   module subroutine build_region_registry(this)
      class(charge), intent(inout) :: this
      integer, parameter :: buildsurf_init = 6
      integer, dimension(:), allocatable :: reference_type_site
      integer :: i, nrf

      ! this%lattice%chargetrf_type is sized nbas for 'S' (build_surf_full),
      ! which build_from_buildsurf's reference_type(1:nbas) copy assumes --
      ! but surfmat (and this routine) run unconditionally for EVERY calctype
      ! that reaches it, including 'L', where chargetrf_type is TYPE-indexed
      ! (dimension nrec, build_interface_full's convention) and can be smaller
      ! than nbas. build_interface_registry overwrites this%regions right
      ! after for 'L', but only if this call doesn't crash first. Expand
      ! defensively exactly as build_interface_registry does.
      allocate (reference_type_site(max(this%lattice%nbas, 1)))
      nrf = size(this%lattice%chargetrf_type)
      if (nrf < 1) then
         reference_type_site(:) = 0
      else
         do i = 1, this%lattice%nbas
            reference_type_site(i) = this%lattice%chargetrf_type(mod(i - 1, nrf) + 1)
         end do
      end if

      call this%regions%build_from_buildsurf(this%lattice%nbas, this%lattice%nlay, &
                                              this%qz, this%wssurf, &
                                              reference_type_site, buildsurf_init)

      ! Per-region alignment shift V_r (B7 §1.3). Zero until B7.4 solves for it;
      ! region A stays the gauge anchor at exactly zero regardless.
      if (allocated(this%region_shift)) deallocate (this%region_shift)
      allocate (this%region_shift(max(this%regions%nregion, 1)))
      this%region_shift(:) = 0.0d0
   end subroutine build_region_registry

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> B7.5: construct this%regions as a genuinely two-sided registry (region
   !> A, active zone, region B) for the calctype='L' (layered/interface)
   !> path, from the same 2D cluster data (qz, wssurf, nbas) `surfmat`
   !> already produces -- the cluster builder for 'L' is geometry-general
   !> (B7 §2.10: madl2d runs off arbitrary QPPZ pairs) so no new Madelung
   !> setup is needed, only a different region partition of the same rows.
   !>
   !> `nlay_a`/`nlay_b` (row counts of the two frozen boundaries) and the
   !> optional `fermi_a`/`fermi_b` come from the &charge namelist (B7.5); the
   !> active-zone width is whatever remains: nbas - nlay_a - nlay_b.
   !---------------------------------------------------------------------------
   module subroutine build_interface_registry(this)
      class(charge), intent(inout) :: this
      integer :: nlay_active, i, nrf, kind_b
      integer, dimension(:), allocatable :: reference_type_site

      ! Geometry-derived default: CENTRE the active zone in the Madelung stack.
      !
      ! The active layers' deviation charge lands on rows nlay_a+1 .. nlay_a+nlay
      ! (interfacepot step 1). The two deep probes are the extreme rows of the
      ! frozen regions, so an off-centre split puts one probe adjacent to the
      ! charge and the other ~nbas rows away -- and the alignment solver then
      ! correctly reports a nonzero V_B for a geometry that is *physically*
      ! symmetric, because the ROW layout it was handed is not. Measured on
      ! fccCu111_AA (nbas=49, one active layer): nlay_a=nlay_b=1 gives
      ! V(B) = -0.0109 Ry, nlay_a=nlay_b=24 gives V(B) = 0 exactly, as the
      ! A | A identity requires.
      !
      ! Defaulting to the centred split makes the common case right without the
      ! user computing (nbas - nlay)/2. Explicit &charge values still win, so
      ! deliberately asymmetric stacks stay expressible.
      if (this%nlay_a <= 0 .and. this%nlay_b <= 0) then
         this%nlay_a = max((this%lattice%nbas - this%lattice%nlay)/2, 1)
         this%nlay_b = this%lattice%nbas - this%lattice%nlay - this%nlay_a
         if (this%nlay_b < 1) then
            this%nlay_a = max(this%lattice%nbas - this%lattice%nlay - 1, 1)
            this%nlay_b = 1
         end if
         if (rank == 0) then
            call g_logger%info('build_interface_registry: &charge nlay_a/nlay_b not set; '// &
                               'centring the active zone in the '//fmt('i0', this%lattice%nbas)// &
                               '-row Madelung stack: nlay_a= '//fmt('i0', this%nlay_a)// &
                               ' nlay_b= '//fmt('i0', this%nlay_b)//'. An off-centre split '// &
                               'makes the two deep probes asymmetric and yields a nonzero '// &
                               'alignment shift for a physically symmetric cell.', &
                               __FILE__, __LINE__)
         end if
      end if

      nlay_active = this%lattice%nbas - this%nlay_a - this%nlay_b
      if (nlay_active < 0) then
         call g_logger%fatal('build_interface_registry: nlay_a + nlay_b exceeds nbas -- '// &
                              'check the &charge namelist nlay_a/nlay_b against lattice%nbas', &
                              __FILE__, __LINE__)
      end if

      ! Per-site reference type, dimension nbas. Each row gets the frozen type
      ! it reverts to, which differs by region:
      !
      !   rows 1..nlay_a          region A's own frozen types (1..nbulk_a)
      !   active rows             lattice%chargetrf_type, the per-active-type
      !                           reference build_interface_full assigned
      !   rows ..nbas             region B's frozen types (nbulk_a+1..nbulk)
      !
      ! The frozen rows matter beyond bookkeeping: `boundary_nef` reads this to
      ! weight the compensation split by side-resolved N(E_F) (B7 §1.5), so a
      ! boundary row carrying no valid type collapses that weighting to 50/50
      ! and puts half the compensation charge into vacuum.
      allocate (reference_type_site(max(this%lattice%nbas, 1)))
      reference_type_site(:) = 0
      nrf = size(this%lattice%chargetrf_type)
      do i = 1, this%lattice%nbas
         if (i <= this%nlay_a) then
            ! Region A frozen rows: cycle its own inequivalent types.
            reference_type_site(i) = mod(i - 1, max(this%lattice%nbulk_bulk, 1)) + 1
         else if (i > this%nlay_a + nlay_active) then
            ! Region B frozen rows: the second frozen block, offset past A.
            ! For a vacuum region B that block is the single empty-sphere type.
            reference_type_site(i) = this%lattice%nbulk_bulk + &
                                     mod(i - this%nlay_a - nlay_active - 1, &
                                         max(this%lattice%nbulk - this%lattice%nbulk_bulk, 1)) + 1
         else if (nrf >= 1) then
            ! Active rows: the reference each active type reverts to.
            reference_type_site(i) = this%lattice%chargetrf_type(mod(i - this%nlay_a - 1, nrf) + 1)
         end if
      end do

      ! B7.6: region B is vacuum when &lattice region_b_kind says so. The kind
      ! is what activates the vacuum-aware behaviour already written into the
      ! consumers: gauge_anchor skips it when choosing the alignment anchor,
      ! and boundary_nef returns zero so it receives no compensation charge.
      kind_b = region_kind_lead_b
      if (trim(this%lattice%region_b_kind) == 'vacuum') kind_b = region_kind_vacuum

      call this%regions%build_from_interface(this%lattice%nbas, this%nlay_a, nlay_active, &
                                              this%qz, this%wssurf, reference_type_site, &
                                              this%fermi_a, this%fermi_b, kind_b=kind_b)

      ! Per-region alignment shift V_r (B7 §1.3); region A is the gauge anchor.
      if (allocated(this%region_shift)) deallocate (this%region_shift)
      allocate (this%region_shift(max(this%regions%nregion, 1)))
      this%region_shift(:) = 0.0d0
   end subroutine build_interface_registry

end submodule charge_interface
