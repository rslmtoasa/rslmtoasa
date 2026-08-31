!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Region Registry
!
! DESCRIPTION:
!> @brief
!> Explicit per-site region bookkeeping for the embedding cluster (B7.1).
!>
!> Today’s `charge%surfpot` locates the active/frozen boundary implicitly:
!> `init = 6` is a magic offset encoding the number of leading (vacuum-side)
!> frozen rows, and `iq = init + ibas` maps a layer index onto a matrix row.
!> That works only because the buildsurf layout is fixed (one vacuum side,
!> one bulk side, `init` hardcoded). B7 (interfaces and vacuum leads) needs
!> more than one frozen reference region to coexist in a single cluster, so
!> the offset arithmetic has to become data.
!>
!> This module supplies that data as an explicit registry, one entry per
!> cluster site (dimension `nbas`, matching `charge%qz`/`charge%wssurf`):
!>   - `region_id`      which frozen/active region the site belongs to
!>   - `active`         .true. for sites that relax self-consistently
!>   - `layer_index`    physical atomic layer (1..nlay for active sites)
!>   - `z`               the site’s z coordinate, carried as DATA (not
!>                        derived from layer_index) so relaxed-z is later a
!>                        parameter change, not a rewrite (B7 §2.10)
!>   - `w`               per-site Wigner-Seitz radius
!>   - `reference_type`  which frozen parameter class a site reverts to when
!>                        frozen (index into the region’s reference-type
!>                        table; mirrors lattice%chargetrf_type)
!>
!> `build_from_buildsurf` reproduces today’s buildsurf index map exactly, as
!> a registry instance, so it can be validated against the surface regression
!> without touching `charge%surfpot` itself (which remains the permanent
!> regression oracle -- B7 explicitly forbids changing its behaviour). Later
!> work (B7.3/B7.4) is expected to consume a registry built by a genuinely
!> two-sided constructor; this module does not gate that on rewiring surfpot.
!------------------------------------------------------------------------------
module region_registry_mod
   use precision_mod, only: rp
   use logger_mod, only: g_logger
   use string_mod, only: fmt
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   !> Region kinds. VACUUM and BULK are the two ends of today’s buildsurf
   !> slab; LEAD_A/LEAD_B are reserved for the two-sided A|B interface
   !> geometry (B7.3+) so the registry does not need to grow a new field
   !> when that lands.
   integer, parameter, public :: region_kind_vacuum = 1
   integer, parameter, public :: region_kind_bulk = 2
   integer, parameter, public :: region_kind_active = 3
   integer, parameter, public :: region_kind_lead_a = 4
   integer, parameter, public :: region_kind_lead_b = 5

   !> One row per region: metadata shared by every site assigned to it.
   type, public :: region_descriptor
      integer :: kind = region_kind_active
      character(len=32) :: name = ''
      !> .true. if sites of this region are held frozen (never updated).
      logical :: frozen = .false.
      !> B7.4: the region’s own converged Fermi level, in the absolute energy
      !> scale its frozen parameter set was generated on (Ry). Used ONLY for
      !> the analytic initial guess V_r = E_F^(A) - E_F^(r) and for the
      !> consistency check against the converged fixed point (B7 §1.3) --
      !> never as the load-bearing determination, because that would require a
      !> cross-run absolute energy scale the code does not yet guarantee
      !> (§2.1, G-B7-2). `fermi_known` says whether it was actually supplied;
      !> when .false. the guess is skipped and the check is not performed.
      real(rp) :: fermi = 0.0_rp
      logical :: fermi_known = .false.
   end type region_descriptor

   !> Explicit per-site region registry, dimension nbas.
   type, public :: region_registry
      !> Number of sites (== charge%lattice%nbas for a buildsurf-derived
      !> instance; the extended interaction zone, user-sized).
      integer :: nsite = 0
      !> Number of distinct regions.
      integer :: nregion = 0
      !> Region metadata, dimension(nregion).
      type(region_descriptor), dimension(:), allocatable :: region

      !> Per-site region id, 1..nregion.
      integer, dimension(:), allocatable :: region_id
      !> Per-site active/frozen mask. .true. = active (relaxes self-
      !> consistently); .false. = frozen (parameters imported, never
      !> updated). This is the registry-explicit form of ″the charge loop
      !> fills only rows init+1..nbas″: frozen sites are deliberately held
      !> at zero deviation charge, not accidentally skipped.
      logical, dimension(:), allocatable :: active
      !> Physical atomic layer index. Meaningful (1..nlay) for active
      !> sites; 0 for frozen sites, which are not part of the layer stack.
      integer, dimension(:), allocatable :: layer_index
      !> Site z coordinate, carried as DATA -- not derived from
      !> layer_index -- so relaxed-z is a parameter change later (B7 §2.10).
      real(rp), dimension(:), allocatable :: z
      !> Per-site Wigner-Seitz radius (mirrors charge%wssurf).
      real(rp), dimension(:), allocatable :: w
      !> Reference-type index a frozen site reverts to (mirrors
      !> lattice%chargetrf_type); 0 for active sites, which have no frozen
      !> reference.
      integer, dimension(:), allocatable :: reference_type
   contains
      procedure :: restore_to_default => region_registry_restore_to_default
      procedure :: build_from_buildsurf => region_registry_build_from_buildsurf
      procedure :: build_from_interface => region_registry_build_from_interface
      procedure :: dump => region_registry_dump
      procedure :: is_active => region_registry_is_active
      procedure :: region_of => region_registry_region_of
      !> B7.4 alignment solver (B7 §1.3). See the procedure headers.
      procedure :: gauge_anchor => alignment_gauge_anchor
      procedure :: initial_guess => alignment_initial_guess
      procedure :: align_update => alignment_update
      procedure :: consistency_check => alignment_consistency_check
      procedure :: deep_probe_site => region_registry_deep_probe_site
      final :: region_registry_destructor
   end type region_registry

   interface region_registry
      procedure :: region_registry_constructor
   end interface region_registry

   public :: alignment_gauge_anchor
   public :: alignment_initial_guess
   public :: alignment_update
   public :: alignment_consistency_check

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> B7.4: the gauge anchor region -- the one whose alignment shift is held
   !> identically zero (B7 §1.3).
   !>
   !> @details
   !> Only DIFFERENCES of the region shifts are physical:
   !>
   !>     V_B - V_A = E_F^(A) - E_F^(B)      (the contact potential)
   !>
   !> so the shifts are fixed only up to one overall gauge, and one region must
   !> be pinned or the fixed point has a flat direction.
   !>
   !> The anchor is the **first frozen, non-vacuum region** -- ″region A″ in the
   !> plan’s notation. Both qualifiers are load-bearing:
   !>
   !>  - **Frozen.** The residual driving the fixed point is ″deep inside region
   !>    r the site must BE neutral bulk r″, which is a statement about a
   !>    reference parameter set. An active region has no frozen reference to be
   !>    -- it is what relaxes -- so it cannot carry the gauge. For the
   !>    buildsurf layout (vacuum, active, bulk) this is what makes the anchor
   !>    the bulk region rather than the active layers that happen to precede
   !>    it in registry order.
   !>  - **Non-vacuum.** Vacuum carries no states at E_F, so ″deep vacuum is
   !>    neutral bulk vacuum″ is not a condition that can drive a residual;
   !>    anchoring there would pin the gauge to the one region that cannot
   !>    report on it.
   !>
   !> Returns 0 for an empty registry (no region to anchor).
   !---------------------------------------------------------------------------
   pure integer function alignment_gauge_anchor(this) result(ianchor)
      class(region_registry), intent(in) :: this
      integer :: ir

      ianchor = 0
      do ir = 1, this%nregion
         if (this%region(ir)%kind == region_kind_vacuum) cycle
         if (.not. this%region(ir)%frozen) cycle
         ianchor = ir
         return
      end do

      ! No frozen non-vacuum region: fall back to any non-vacuum region, then
      ! to the first region, so the gauge is always pinned rather than left
      ! free. A flat direction in the fixed point is worse than an imperfect
      ! anchor choice.
      do ir = 1, this%nregion
         if (this%region(ir)%kind /= region_kind_vacuum) then
            ianchor = ir
            return
         end if
      end do
      if (this%nregion > 0) ianchor = 1
   end function alignment_gauge_anchor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> B7.4: the deep-probe site of a region -- the site at which ″deep inside
   !> region r the site must BE neutral bulk r″ is evaluated (B7 §1.3).
   !>
   !> @details
   !> `surfpot` probes `dss(1,.)` for deep vacuum and `dss(nbas,.)` for deep
   !> bulk -- the two extreme rows, which for the one-sided buildsurf layout
   !> are exactly the deepest frozen sites of the two ends. Generalized to the
   !> registry, the deep probe of region r is its extreme site in z, taken away
   !> from the interface: the MINIMUM-z site for a region lying below the
   !> active zone and the MAXIMUM-z site for one lying above it.
   !>
   !> Determining which side a region is on from z rather than from its index
   !> is deliberate -- registry order is a bookkeeping choice, z is physics, and
   !> the two need not agree once a genuinely two-sided builder (B7.5) exists.
   !> The comparison is against the mean z of the region’s own sites versus the
   !> mean z of the whole cluster.
   !>
   !> Returns 0 if the region has no sites.
   !---------------------------------------------------------------------------
   pure integer function region_registry_deep_probe_site(this, ireg) result(isite)
      class(region_registry), intent(in) :: this
      integer, intent(in) :: ireg
      real(rp) :: zsum_reg, zsum_all, zmean_reg, zmean_all, zbest
      integer :: i, nreg_sites
      logical :: below

      isite = 0
      if (ireg < 1 .or. ireg > this%nregion) return
      if (this%nsite < 1) return

      zsum_reg = 0.0_rp
      nreg_sites = 0
      zsum_all = 0.0_rp
      do i = 1, this%nsite
         zsum_all = zsum_all + this%z(i)
         if (this%region_id(i) == ireg) then
            zsum_reg = zsum_reg + this%z(i)
            nreg_sites = nreg_sites + 1
         end if
      end do
      if (nreg_sites == 0) return

      zmean_reg = zsum_reg/real(nreg_sites, rp)
      zmean_all = zsum_all/real(this%nsite, rp)
      below = (zmean_reg <= zmean_all)

      ! Extreme site of the region, on the side away from the cluster centre.
      zbest = 0.0_rp
      do i = 1, this%nsite
         if (this%region_id(i) /= ireg) cycle
         if (isite == 0) then
            isite = i
            zbest = this%z(i)
         else if (below) then
            if (this%z(i) < zbest) then
               isite = i
               zbest = this%z(i)
            end if
         else
            if (this%z(i) > zbest) then
               isite = i
               zbest = this%z(i)
            end if
         end if
      end do
   end function region_registry_deep_probe_site

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> B7.4: the analytic initial guess for the region alignment shifts,
   !> V_r = E_F^(anchor) - E_F^(r) (B7 §1.3).
   !>
   !> @details
   !> This is an INITIAL GUESS ONLY, never the determination. B7 §1.3 is
   !> explicit about why: the closed form requires a reliable cross-run
   !> absolute energy scale, which the code does not currently guarantee
   !> (§2.1). Used as a starting point it usually saves several iterations;
   !> used as the answer it would silently produce a plausible but wrong dipole
   !> barrier -- the worst failure mode available to this work package.
   !>
   !> Regions whose `fermi_known` is .false. get a zero guess (the solver then
   !> starts from the gauge and converges anyway, just more slowly). The anchor
   !> is always exactly zero by construction.
   !---------------------------------------------------------------------------
   pure subroutine alignment_initial_guess(this, ianchor, shift)
      class(region_registry), intent(in) :: this
      integer, intent(in) :: ianchor
      real(rp), dimension(:), intent(out) :: shift
      integer :: ir

      shift(:) = 0.0_rp
      if (ianchor < 1 .or. ianchor > this%nregion) return
      if (.not. this%region(ianchor)%fermi_known) return

      do ir = 1, this%nregion
         if (ir == ianchor) cycle
         if (.not. this%region(ir)%fermi_known) cycle
         shift(ir) = this%region(ianchor)%fermi - this%region(ir)%fermi
      end do
   end subroutine alignment_initial_guess

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> B7.4: one mixed fixed-point update of the region alignment shifts from
   !> the deep-probe residuals (B7 §1.3).
   !>
   !> @details
   !> Deep inside region r the site must BE neutral bulk r. The residual that
   !> expresses this is the deep-probe difference
   !>
   !>     resid_r = -[ (dV(deep-r) + V_r) - (dV(deep-anchor) + V_anchor) ]
   !>
   !> and the update is the ordinary SCF mix
   !>
   !>     V_r <- V_r + vmix * resid_r
   !>
   !> mixed exactly as `surfpot` mixes any other SCF quantity. `vmix` is the
   !> damping the alignment needs: the two-region system has a capacitor-like
   !> soft mode (a small charge transfer moves the whole relative level), so
   !> the undamped update overshoots.
   !>
   !> The anchor is forced to exactly zero AFTER the update rather than merely
   !> skipped, so no accumulation of round-off can let the gauge drift. That
   !> matters: a drifting anchor is exactly the negative control the §5.3
   !> oracle uses to prove the test has teeth.
   !>
   !> @param[in]     probe   dV at the deep probe site of each region.
   !> @param[in]     vmix    SCF mixing parameter.
   !> @param[in,out] shift   region shifts V_r, updated in place.
   !> @param[out]    maxresid  largest |residual| over the non-anchor regions,
   !>                          for the caller’s convergence report.
   !---------------------------------------------------------------------------
   pure subroutine alignment_update(this, ianchor, probe, vmix, shift, maxresid)
      class(region_registry), intent(in) :: this
      integer, intent(in) :: ianchor
      real(rp), dimension(:), intent(in) :: probe
      real(rp), intent(in) :: vmix
      real(rp), dimension(:), intent(inout) :: shift
      real(rp), intent(out) :: maxresid
      real(rp) :: probe_anchor, resid
      integer :: ir

      maxresid = 0.0_rp
      if (ianchor < 1 .or. ianchor > this%nregion) return

      probe_anchor = probe(ianchor) + shift(ianchor)

      do ir = 1, this%nregion
         if (ir == ianchor) cycle
         resid = -((probe(ir) + shift(ir)) - probe_anchor)
         shift(ir) = shift(ir) + vmix*resid
         maxresid = max(maxresid, abs(resid))
      end do

      ! Re-pin the gauge exactly, never by accumulation.
      shift(ianchor) = 0.0_rp
   end subroutine alignment_update

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> B7.4: compare the converged alignment fixed point against the analytic
   !> value E_F^(anchor) - E_F^(r) (B7 §1.3).
   !>
   !> @details
   !> This converts the most fragile ingredient in the work package -- the
   !> cross-run absolute energy scale -- from load-bearing into a DIAGNOSTIC.
   !> The fixed point is the answer; the analytic value is the check. If they
   !> disagree beyond `tol`, the absolute-zero bookkeeping (G-B7-2) is broken
   !> and the caller must warn loudly, because a silent misalignment produces a
   !> plausible but wrong dipole barrier.
   !>
   !> Returns the largest discrepancy over regions with a known Fermi level,
   !> and `checked = .false.` if no region could be checked at all.
   !---------------------------------------------------------------------------
   pure subroutine alignment_consistency_check(this, ianchor, shift, worst, iworst, checked)
      class(region_registry), intent(in) :: this
      integer, intent(in) :: ianchor
      real(rp), dimension(:), intent(in) :: shift
      real(rp), intent(out) :: worst
      integer, intent(out) :: iworst
      logical, intent(out) :: checked
      real(rp) :: analytic, diff
      integer :: ir

      worst = 0.0_rp
      iworst = 0
      checked = .false.
      if (ianchor < 1 .or. ianchor > this%nregion) return
      if (.not. this%region(ianchor)%fermi_known) return

      do ir = 1, this%nregion
         if (ir == ianchor) cycle
         if (.not. this%region(ir)%fermi_known) cycle
         analytic = this%region(ianchor)%fermi - this%region(ir)%fermi
         diff = abs(shift(ir) - analytic)
         checked = .true.
         if (diff > worst) then
            worst = diff
            iworst = ir
         end if
      end do
   end subroutine alignment_consistency_check

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor. Produces an empty (zero-site) registry; call
   !> build_from_buildsurf (or a future two-sided builder) to populate it.
   !---------------------------------------------------------------------------
   function region_registry_constructor() result(obj)
      type(region_registry) :: obj

      call obj%restore_to_default()
   end function region_registry_constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor.
   !---------------------------------------------------------------------------
   subroutine region_registry_destructor(this)
      type(region_registry) :: this

      if (allocated(this%region)) deallocate (this%region)
      if (allocated(this%region_id)) deallocate (this%region_id)
      if (allocated(this%active)) deallocate (this%active)
      if (allocated(this%layer_index)) deallocate (this%layer_index)
      if (allocated(this%z)) deallocate (this%z)
      if (allocated(this%w)) deallocate (this%w)
      if (allocated(this%reference_type)) deallocate (this%reference_type)
   end subroutine region_registry_destructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset to an empty, zero-site registry.
   !---------------------------------------------------------------------------
   subroutine region_registry_restore_to_default(this)
      class(region_registry), intent(inout) :: this

      if (allocated(this%region)) deallocate (this%region)
      if (allocated(this%region_id)) deallocate (this%region_id)
      if (allocated(this%active)) deallocate (this%active)
      if (allocated(this%layer_index)) deallocate (this%layer_index)
      if (allocated(this%z)) deallocate (this%z)
      if (allocated(this%w)) deallocate (this%w)
      if (allocated(this%reference_type)) deallocate (this%reference_type)

      this%nsite = 0
      this%nregion = 0
      allocate (this%region(0))
      allocate (this%region_id(0))
      allocate (this%active(0))
      allocate (this%layer_index(0))
      allocate (this%z(0))
      allocate (this%w(0))
      allocate (this%reference_type(0))
   end subroutine region_registry_restore_to_default

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Build a registry that reproduces today’s buildsurf/surfpot layout
   !> exactly, as an explicit instance, so it can be validated against the
   !> surface regression without touching surfpot itself.
   !>
   !> Reproduces the index map documented in B7 §2.10 and verified against
   !> charge%surfpot (source/charge.f90):
   !>   - rows 1..init            vacuum-side frozen rows (init = 6, the
   !>                              historical magic offset)
   !>   - rows init+1..init+nlay  active layers (charge loop fills these)
   !>   - rows init+nlay+1..nbas  bulk-side frozen rows
   !> The charge loop in surfpot only ever touches rows init+1..nbas, so the
   !> vacuum rows 1..init exist in the geometry but are deliberately held at
   !> exactly zero deviation charge (B7 §1.5). That is encoded here directly:
   !> the vacuum region is frozen, not merely ″not yet visited″.
   !>
   !> @param[in] nbas   Total number of sites in the extended interaction
   !>                    zone (charge%lattice%nbas). User-sized; no width is
   !>                    assumed.
   !> @param[in] nlay   Number of active physical layers
   !>                    (charge%lattice%nlay).
   !> @param[in] z      Per-site z coordinate, dimension(nbas)
   !>                    (charge%qz, already fully z-general -- madl2d runs
   !>                    off QPPZ = QZ(IQ)-QZ(JQ) for arbitrary pairs).
   !> @param[in] w      Per-site Wigner-Seitz radius, dimension(nbas)
   !>                    (charge%wssurf, already per-site).
   !> @param[in] reference_type  Per-site frozen-reference class index,
   !>                    dimension(nbas) (lattice%chargetrf_type). Optional;
   !>                    if absent, reference_type is left 0 for every site
   !>                    (no reference table available in that context).
   !> @param[in] init    Vacuum-side frozen row count. Optional, default 6
   !>                    -- the value hardcoded in charge%surfpot today.
   !---------------------------------------------------------------------------
   subroutine region_registry_build_from_buildsurf(this, nbas, nlay, z, w, reference_type, init)
      class(region_registry), intent(inout) :: this
      integer, intent(in) :: nbas
      integer, intent(in) :: nlay
      real(rp), dimension(:), intent(in) :: z
      real(rp), dimension(:), intent(in) :: w
      integer, dimension(:), intent(in), optional :: reference_type
      integer, intent(in), optional :: init
      integer :: init_
      integer :: i, ilay
      integer, parameter :: reg_vacuum = 1
      integer, parameter :: reg_active = 2
      integer, parameter :: reg_bulk = 3

      init_ = 6
      if (present(init)) init_ = init

      call this%restore_to_default()

      this%nsite = nbas
      this%nregion = 3

      deallocate (this%region)
      allocate (this%region(this%nregion))
      this%region(reg_vacuum) = region_descriptor(kind=region_kind_vacuum, name='vacuum', frozen=.true.)
      this%region(reg_active) = region_descriptor(kind=region_kind_active, name='active', frozen=.false.)
      this%region(reg_bulk) = region_descriptor(kind=region_kind_bulk, name='bulk', frozen=.true.)

      deallocate (this%region_id, this%active, this%layer_index, this%z, this%w, this%reference_type)
      allocate (this%region_id(nbas))
      allocate (this%active(nbas))
      allocate (this%layer_index(nbas))
      allocate (this%z(nbas))
      allocate (this%w(nbas))
      allocate (this%reference_type(nbas))

      this%z(1:nbas) = z(1:nbas)
      this%w(1:nbas) = w(1:nbas)
      this%reference_type(:) = 0
      if (present(reference_type)) this%reference_type(1:nbas) = reference_type(1:nbas)

      ! Rows 1..init: vacuum-side frozen. Held at exactly zero deviation
      ! charge deliberately (B7 §1.5) -- this is the explicit statement of
      ! what today’s surfpot does implicitly by never writing tdq() there.
      do i = 1, min(init_, nbas)
         this%region_id(i) = reg_vacuum
         this%active(i) = .false.
         this%layer_index(i) = 0
      end do

      ! Rows init+1..init+nlay: the active layers.
      ilay = 0
      do i = init_ + 1, min(init_ + nlay, nbas)
         ilay = ilay + 1
         this%region_id(i) = reg_active
         this%active(i) = .true.
         this%layer_index(i) = ilay
      end do

      ! Rows init+nlay+1..nbas: bulk-side frozen.
      do i = init_ + nlay + 1, nbas
         this%region_id(i) = reg_bulk
         this%active(i) = .false.
         this%layer_index(i) = 0
      end do
   end subroutine region_registry_build_from_buildsurf

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> B7.5: build a genuinely two-sided registry -- region A (frozen, low-z),
   !> an active zone, region B (frozen, high-z) -- for the calctype=’L’
   !> (layered/interface) path (B7 §1.2, §4 B7.5).
   !>
   !> @details
   !> Generalizes `build_from_buildsurf`’s single-frozen-reference layout to
   !> two independent frozen regions, which is the one structural difference
   !> B7 §1.1-§1.3 requires: two potential zeros that must be aligned (B7.4),
   !> not one. The index map is analogous but symmetric:
   !>
   !>   - rows 1..nlay_a                    region A, frozen
   !>   - rows nlay_a+1..nlay_a+nlay_active  the active zone
   !>   - rows nlay_a+nlay_active+1..nbas    region B, frozen
   !>
   !> Both frozen regions are `region_kind_lead_a`/`region_kind_lead_b` (B7.1
   !> reserved these precisely for this constructor). Sites are held frozen
   !> exactly as vacuum/bulk are in `build_from_buildsurf` -- the charge loop
   !> in `interfacepot` only ever touches sites the registry marks active.
   !>
   !> `fermi_a`/`fermi_b` are OPTIONAL: passing a value <= `fermi_sentinel`
   !> (the huge() default `charge%fermi_a`/`fermi_b` restore to) means ″not
   !> supplied″, exactly as CONTRACT_FROZEN_REGION.md §2 requires -- E_F is
   !> the one quantity a frozen region may legitimately omit, and
   !> `region_descriptor%fermi_known` gates every use of it (initial guess,
   !> consistency check) so an absent value degrades gracefully rather than
   !> injecting a bogus zero.
   !>
   !> @param[in] nbas         Total sites in the cluster.
   !> @param[in] nlay_a       Frozen region-A row count (leading rows).
   !> @param[in] nlay_active  Active-zone row count, following region A.
   !> @param[in] z            Per-site z coordinate, dimension(nbas).
   !> @param[in] w            Per-site Wigner-Seitz radius, dimension(nbas).
   !> @param[in] reference_type  Per-site frozen-reference class index,
   !>                    dimension(nbas) (lattice%chargetrf_type). Optional.
   !> @param[in] fermi_a      Region A source E_F (Ry). Optional; a value
   !>                    <= fermi_sentinel means ″not supplied″.
   !> @param[in] fermi_b      Region B source E_F (Ry). Same convention.
   !> @param[in] fermi_sentinel  Threshold below which fermi_a/fermi_b count
   !>                    as supplied. Optional, default huge(1.0_rp)/2 (any
   !>                    ordinary Ry value is far below it; only the exact
   !>                    huge() default sentinel is excluded).
   !---------------------------------------------------------------------------
   subroutine region_registry_build_from_interface(this, nbas, nlay_a, nlay_active, z, w, &
                                                     reference_type, fermi_a, fermi_b, fermi_sentinel, &
                                                     kind_b)
      class(region_registry), intent(inout) :: this
      integer, intent(in) :: nbas
      integer, intent(in) :: nlay_a
      integer, intent(in) :: nlay_active
      real(rp), dimension(:), intent(in) :: z
      real(rp), dimension(:), intent(in) :: w
      integer, dimension(:), intent(in), optional :: reference_type
      real(rp), intent(in), optional :: fermi_a
      real(rp), intent(in), optional :: fermi_b
      real(rp), intent(in), optional :: fermi_sentinel
      !> Region B’s kind. Optional; defaults to `region_kind_lead_b` (the
      !> metallic A | B geometry). Pass `region_kind_vacuum` for A | vacuum
      !> (B7.6): that is what makes `gauge_anchor` skip region B when choosing
      !> the alignment anchor, and what makes `charge%boundary_nef` return zero
      !> so vacuum receives NO compensation charge (B7 §1.5 -- compensation
      !> placed in vacuum does not perturb the work function, it SETS it).
      integer, intent(in), optional :: kind_b
      real(rp) :: sentinel_
      integer :: i, kind_b_
      character(len=8) :: name_b
      integer, parameter :: reg_a = 1
      integer, parameter :: reg_active = 2
      integer, parameter :: reg_b = 3

      sentinel_ = huge(1.0_rp)/2.0_rp
      if (present(fermi_sentinel)) sentinel_ = fermi_sentinel

      kind_b_ = region_kind_lead_b
      if (present(kind_b)) kind_b_ = kind_b
      ! Name region B for what it is. This string is what every alignment and
      ! diagnostic message prints, and it is also what `fix_fermi_to_region`
      ! matches on -- pinning E_F to vacuum is not meaningful, so calling the
      ! region ’vacuum’ rather than ’B’ keeps that user-facing name honest.
      name_b = 'B'
      if (kind_b_ == region_kind_vacuum) name_b = 'vacuum'

      call this%restore_to_default()

      this%nsite = nbas
      this%nregion = 3

      deallocate (this%region)
      allocate (this%region(this%nregion))
      this%region(reg_a) = region_descriptor(kind=region_kind_lead_a, name='A', frozen=.true.)
      this%region(reg_active) = region_descriptor(kind=region_kind_active, name='active', frozen=.false.)
      this%region(reg_b) = region_descriptor(kind=kind_b_, name=name_b, frozen=.true.)

      if (present(fermi_a)) then
         if (fermi_a < sentinel_) then
            this%region(reg_a)%fermi = fermi_a
            this%region(reg_a)%fermi_known = .true.
         end if
      end if
      if (present(fermi_b)) then
         if (fermi_b < sentinel_) then
            this%region(reg_b)%fermi = fermi_b
            this%region(reg_b)%fermi_known = .true.
         end if
      end if

      deallocate (this%region_id, this%active, this%layer_index, this%z, this%w, this%reference_type)
      allocate (this%region_id(nbas))
      allocate (this%active(nbas))
      allocate (this%layer_index(nbas))
      allocate (this%z(nbas))
      allocate (this%w(nbas))
      allocate (this%reference_type(nbas))

      this%z(1:nbas) = z(1:nbas)
      this%w(1:nbas) = w(1:nbas)
      this%reference_type(:) = 0
      if (present(reference_type)) this%reference_type(1:nbas) = reference_type(1:nbas)

      ! Rows 1..nlay_a: region A, frozen.
      do i = 1, min(nlay_a, nbas)
         this%region_id(i) = reg_a
         this%active(i) = .false.
         this%layer_index(i) = 0
      end do

      ! Rows nlay_a+1..nlay_a+nlay_active: the active zone.
      do i = nlay_a + 1, min(nlay_a + nlay_active, nbas)
         this%region_id(i) = reg_active
         this%active(i) = .true.
         this%layer_index(i) = i - nlay_a
      end do

      ! Rows nlay_a+nlay_active+1..nbas: region B, frozen.
      do i = nlay_a + nlay_active + 1, nbas
         this%region_id(i) = reg_b
         this%active(i) = .false.
         this%layer_index(i) = 0
      end do
   end subroutine region_registry_build_from_interface

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Whether site i is active (relaxes self-consistently).
   !---------------------------------------------------------------------------
   pure logical function region_registry_is_active(this, i) result(res)
      class(region_registry), intent(in) :: this
      integer, intent(in) :: i
      res = this%active(i)
   end function region_registry_is_active

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Region descriptor for site i.
   !---------------------------------------------------------------------------
   pure function region_registry_region_of(this, i) result(res)
      class(region_registry), intent(in) :: this
      integer, intent(in) :: i
      type(region_descriptor) :: res
      res = this%region(this%region_id(i))
   end function region_registry_region_of

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Debug dump: one line per site, plus a region summary. Writes to the
   !> given unit if open, otherwise to a fresh file ’region_registry.out’.
   !---------------------------------------------------------------------------
   subroutine region_registry_dump(this, unit, fname)
      class(region_registry), intent(in) :: this
      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: fname
      integer :: iu, i, ir
      character(len=64) :: fname_
      logical :: opened_here

      opened_here = .false.
      if (present(unit)) then
         iu = unit
      else
         fname_ = 'region_registry.out'
         if (present(fname)) fname_ = fname
         open (newunit=iu, file=trim(fname_), action='write', status='replace')
         opened_here = .true.
      end if

      write (iu, '(a)') '# Region registry dump (B7.1)'
      write (iu, '(a,i0,a,i0)') '# nsite= ', this%nsite, '  nregion= ', this%nregion
      write (iu, '(a)') '#'
      write (iu, '(a)') '# regions:'
      do ir = 1, this%nregion
         write (iu, '(a,i0,a,a,a,l1)') '#   id=', ir, ' name=', trim(this%region(ir)%name), &
            ' frozen=', this%region(ir)%frozen
      end do
      write (iu, '(a)') '#'
      write (iu, '(a)') '#  site  region_id  region_name  active  layer_index'// &
         '            z               w  reference_type'
      do i = 1, this%nsite
         ir = this%region_id(i)
         write (iu, '(2x,i5,4x,i5,4x,a12,4x,l1,6x,i5,2x,f14.8,2x,f14.8,4x,i5)') &
            i, ir, trim(this%region(ir)%name), this%active(i), this%layer_index(i), &
            this%z(i), this%w(i), this%reference_type(i)
      end do

      if (opened_here) close (iu)
   end subroutine region_registry_dump

end module region_registry_mod
