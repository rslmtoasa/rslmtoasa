# Madelung conventions contract

**Status:** **SIGNED** by Anders (2026-07-25) — gates **G-B7-1** and **G-B7-4**.
Revised post-signature per the sign-off note; see "Sign-off note and its
consequences" below.
**Scope:** the electrostatic path — `bulkpot`, `imppot`, `surfpot`, `impmad`,
`bulkmat`, `surfmat`, `madl2d`, `madmat`/`strx00`.
**Backed by:** `tests/unit/test_madelung_conventions.f90`
(`UnitMadelungConventions`) and `tests/unit/test_madl2d_exh_guard.f90`
(`UnitMadl2dExhGuard`). Every entry below is pinned by a passing test, not by
derivation or provenance.

Once signed, these conventions are **pinned and never re-derived in a later
session**. A future change that contradicts one of them must fail a unit test.

---

## Summary for the signer

Four of the six issues B7 §2 raised as suspected bugs turned out **not to be
bugs**. The audit was verified against the sources rather than taken on trust,
and the corrections below reduce B7.0's remaining scope substantially.

| B7 §  | Claim in the plan | Verified finding |
|---|---|---|
| 2.1 | `vmad` may not persist across runs | **Persists.** No new field needed |
| 2.3 | Factor-of-two inconsistency between `bulkpot` and `imppot` | **No bug.** Two different arrays, both correct |
| 2.5 | `exh = -exh` corrupts the dipole matrix | **Real bug, fixed** (`6e87860`) — but numerically dormant; no B6 number moved |
| 2.6 | `AM` half-space gauge may need symmetrizing | **No change needed** if Q = 0 is enforced |
| 2.7 | `SWS` baked into l≥1 prefactors breaks two-region runs | **Withdrawn.** `S` is a dimensionless unit scale that cancels; an average `W` is well defined for two leads (C5) |
| 2.8 | On-site 2δq/w treated three different ways | **Partly wrong.** Present in `impmad` and `surfmat`; the term is `2·(S/w_i)`, dimensionless, not `2δq/w` in Ry (C0, C2) |

---

## C0 — Overall unit convention, and what `S` means

**Potential is in Rydberg; charge in electrons; length in bohr.**
In Rydberg atomic units `e² = 2 Ry·bohr`, so the potential of a point charge
`q` at distance `r` is `V = 2q/r`.

### `S` is dimensionless, `w` is dimensional — they are not the same quantity

This is the distinction the sign-off asked to be made explicit, and getting it
backwards is easy:

| symbol | code | units | meaning |
|---|---|---|---|
| `S` | `charge%sws` | **dimensionless** (units of `alat`) | system-wide average WS radius, relative |
| `S` in bohr | `sws*alat*ang2au` | bohr | the same, made dimensional |
| `w_i` | `charge%wssurf(i)` | **bohr** | per-site WS radius |
| `w_avg` | `lattice%wav` | **Å** (→ bohr via `ang2au`) | system-wide average, dimensional |

`charge%sws` is derived from `charge%vol`, which is built from `bsx/bsy/bsz` —
and those come from `lattice%cr`, i.e. **crystal coordinates in units of
`alat`**. So `sws` is a *relative* radius. `lattice%wav` is built from a volume
in cubic Ångström and is *dimensional*. Numerically, for fccCu001
(`alat = 3.614 Å`, `wav = 1.41237 Å`): `sws ≈ 0.391`, and
`sws·alat·ang2au = 2.66899 bohr = wav·ang2au`. They agree, as they must —
but only after the conversion.

### The `1/S` convention

The `madl2d` header comment *"Potential in units of e²/(2S)"* is **accurate**.
The kernel carries an overall `TWOS = 2·SWS`, and `surfpot` divides by
`wsms = sws·alat·ang2au` at the call site. The **net** quantity reaching
`potential%vmad` is in Rydberg, consistent with C1.

The on-site term in `surfmat` is written in the *same* `1/S` convention:

```fortran
dss(i,i) += 2.0d0*(sws*alat*ang2au / wssurf(i))     ! = 2·(S/w_i)
```

This is `2·(S/w_i)`, **dimensionless** — not `2/w` in Rydberg. For uniform `w`
it is exactly `2.000000`. The Rydberg dimension is restored by the same call-site
`/wsms` that de-scales the rest of the kernel. **So `S` (the average) and `w_i`
(per-site) coexist by design: the average sets the unit scale, the per-site
radius carries the local physics.**

> **Signed action:** comments only; no code change. C2 and C5 below are corrected
> accordingly.

---

## C1 — The factor of two *(resolves §2.3 — there is no inconsistency)*

**There are two distinct Madelung arrays, with two different conventions. Both
are correct, and both must be left alone.**

| array | built by | stores | contracted as |
|---|---|---|---|
| `AMAD` (local, from `mad.mat`) | `madmat` → `strx00` | bare `1/r` | `2.d0*AMAD(...)` in `bulkpot` |
| `this%amad` (component) | `impmad` | `2/r` (factor folded in) | bare, in `imppot` |

`strx00` builds an Ewald sum whose real-space term is `erfc(a·R)/R` with **no
prefactor**, so `AMAD` is in units of `1/r` and `bulkpot`'s explicit `2.d0*` is
required. `impmad` writes `2_rp*dd` off-diagonal and `2/wsimp` on-diagonal, so
`imppot`'s bare contraction is required.

The plan's §2.3 suspicion — *"either `this%amad` is built with the 2 folded in,
or one of them is off by a factor of two"* — resolves to the **first** branch.

- **Pinned by:** `test_madelung_conventions.f90::test_factor_of_two` (C1),
  validation-ladder rung §5.1.
- **Action:** none. Do **not** "harmonize" these two call sites.

---

## C2 — The on-site term *(corrects §2.8)*

**Convention: the on-site term is `2·(S/w_i)`, dimensionless**, in the `1/S`
convention of C0 — *not* `2·δq/w` in Rydberg. The Rydberg dimension is restored
by the call-site division by `wsms`. For uniform `w` the term is exactly `2`.

> An earlier draft of this document stated the term as `2·δq/w` in Rydberg. That
> was wrong about the units — see C0 for the dimensional analysis. The physics is
> unchanged; the bookkeeping is not.

The plan states the term is handled three different ways, with `imppot` adding
none. That is **not** what the code does:

- `bulkpot` — adds `VADD = 2*TDQ/RMAX`, per class. **Present.**
- `impmad` — writes `2/wsimp(ii)` onto the `this%amad` diagonal, so `imppot`
  picks it up through the ordinary contraction. **Present**, just not visible in
  `imppot` itself.
- `surfmat` — adds `2*(sws*alat*ang2au/wssurf(i))` to the `dss` diagonal, using
  a **per-site** `wssurf`. **Present**, at the matrix level.
- `surfpot` — `twooverwsm` is computed and commented out throughout. **Absent**,
  but only as a redundant second application: `surfmat` already put it in `dss`.

So the term is applied consistently; the commented-out code in `surfpot` is dead
weight that reads like a missing feature.

**What is genuinely at risk for B7:** `impmad` sets `wsimp(:) = lattice%wav`, a
single system-wide average, with the in-code comment *"Set as if all the atoms
have the same WS radius. Can be improved later."* `surfmat`'s `wssurf` is
already per-site and is the better pattern. Note this is about the **per-site**
radius `w_i` only — the system-wide average `S` is a separate quantity and stays
(C0).

- **Pinned by:** `test_madelung_conventions.f90::test_onsite_term` (C2), which
  asserts the ratio `S/w_i` is genuinely `w_i`-dependent, so a per-site `w_i`
  cannot be replaced by the average when radii differ.
- **Action:** **none to working routines.** `surfpot`'s commented-out
  `twooverwsm` is dead weight that reads like a missing feature — leave the code
  alone, but this document records that its absence is correct (`surfmat`
  already applies the term via `dss`). New interface code (B7.3) should follow
  `surfmat`'s per-site `wssurf` pattern rather than `impmad`'s `wsimp`.

---

## C3 — The `AM` gauge *(resolves §2.6 — keep the half-space kernel)*

**Decision: keep the existing asymmetric (half-space) `AM`, unchanged, and
enforce Q = 0 as a kernel precondition.**

`madl2d` puts the whole `k∥ = 0` plate term on one triangle:

```fortran
AM(IQ,JQ) = FACGAU*EXPZ - QPPZ*FACERF*ERFCM    ! → −(4π/A)|Δz|
AM(JQ,IQ) = FACGAU*EXPZ + QPPZ*FACERF*ERFCP    ! → 0
```

The symmetric kernel is `−(2π/A)|Δz|` both ways. The difference, acting on a
deviation profile, is

```
D_i = (2π/A)·[ P − z_i·Q ],    P = Σ_j z_j δq_j,    Q = Σ_j δq_j
```

**Verified numerically**, both halves:

1. For a **neutral** profile (Q = 0), `D_i` is a global constant, so the
   symmetric and asymmetric kernels give **identical potential differences** —
   and `dV(deep-B) − dV(deep-A)` is exactly such a difference. The existing
   kernel is therefore **reusable two-sided with no modification**.
2. For **Q ≠ 0**, the difference is exactly `(2π/A)[P − z_i·Q]` — a spurious
   uniform field of magnitude `(2π/A)·Q`, on top of the genuine divergence of a
   charged slab.

This is why **Q = 0 is a precondition on kernel validity, not merely on
physics**, and why B7.3 must compute, report and enforce it every iteration.

- **Pinned by:** `test_madelung_conventions.f90::test_am_gauge` (C3),
  validation-ladder rung §5.2. Pins the gauge, the 4π vs 2π sheet prefactor and
  its sign in one test.
- **Action:** none to the kernel. B7.3 owns the Q = 0 enforcement and reporting.
- **Caveat for §6:** `PM`, the plate-condenser matrix (built, never used, the
  applied-bias hook), carries the **same** triangular gauge choice and will need
  the same treatment when bias goes live. A biased slab is a charged-slab
  problem and Q = 0 no longer saves it.

---

## C4 — The dipole sheet and the l = 1 sign

**`DSZ` is antisymmetric with equal magnitude** (`DSZ(i,j) = +DMDL`,
`DSZ(j,i) = −DMDL`), which is the physically correct `±2πp/A` on the two sides
of a dipole sheet; the step across the sheet is `4πp/A`. This is correct
two-sided as written and needs no change.

- **Pinned by:** `test_madelung_conventions.f90::test_dipole_sheet` (C4).
- **Related fix:** the `exh` overflow guard feeding `DSZ` was corrected in
  `6e87860`; see C6.

---

## C5 — Multipole normalization and `SWS` *(§2.7 — confirmed defect, work remains)*

> **Revised after sign-off.** An earlier draft called this "the one convention
> that is currently WRONG for two regions" and proposed stripping `SWS` from the
> kernel. **Both the premise and the action are withdrawn.** Per the sign-off:
> even with two leads we can still find an average `W`, and working routines are
> not to be edited. See "Why the two-w argument does not hold" below.

`this%sws` is a **single system-wide scalar** — a *dimensionless* average radius
in units of `alat` (C0) — derived from the average cell volume
(`charge.f90:875`):

```fortran
this%sws = (3.0d0*this%vol/4.0d0/pi/this%nq3)**(1.d0/3.d0)
```

It is baked into the l ≥ 1 normalizations inside the kernel:

```
FACP   = SWS·√3      → FACDK, FACDR → DSZ
FACQUA ∝ SWS²        → DS3Z2, DSX2Y2, DSXY
FACGUA ∝ SWS³        → DZ3Z2
```

These are site-independent. `TWOS = 2·SWS` is a global scaling undone at the
call site; `DSS` and the monopole path are unaffected.

### Why the two-w argument does not hold

`S` is a **unit scale**, not a local property. The kernel is expressed in units
of `e²/(2S)` and the call site divides that scaling back out, so `S` cancels
from the delivered potential provided the *same* `S` is used consistently on
both sides — which it is, because there is one kernel. A second region does not
introduce a second unit scale; it introduces additional sites, each with its own
`w_i`, and the per-site radius already enters separately through `wssurf(i)`
(C2). An average `W` remains perfectly well defined for a two-lead system.

The residual concern is narrower than the plan claimed: the l ≥ 1 prefactors fix
the **moment normalization** to the average radius, so multipole moments are
implicitly expressed in units of `S`. That is a *definition*, consistent across
the matrix, not an error — and it only becomes a real accuracy question if the
two regions' radii differ enough that a single normalization is a poor
compromise for both. That is a physics-accuracy judgment for the validation
campaign (B7.7), not a bug to fix now.

- **Action: none.** Working routines are not to be edited. Left as a documented
  convention.
- **For B7.3/B7.7:** if a two-region case shows dipole-coupling error traceable
  to the shared normalization, **report it — do not silently re-normalize the
  kernel.** The sphere-overlap diagnostic (§2.9) is the natural place to surface
  a large radius mismatch.
- **Consequence for B7.0:** item 5 (`SWS` extraction) is **withdrawn**, not
  deferred.

---

## C6 — `madl2d` overflow guard *(§2.5 — fixed, and quantified)*

The `ERFCP == 0` branch assigned `exh = -exh`, reusing the previous G
iteration's value (undefined on the first G of a pair). With `ERFCP` underflowed
to exactly zero, `EXPP*ERFCP` vanishes and the general branch reduces to
`EXF = EXPM*ERFCM`, `EXH = −EXPM*ERFCM = −EXF`. **Fixed in `6e87860`.**

**Measured impact on B6: none.** Instrumenting a real dipole run
(`example/surface/fccCu001/Test`, `dipole_electrostatics = .true.`): the guard
trips **24 276 times across 703 layer pairs**, but `exh` decays monotonically
*before* the guard engages. The largest `|exh|` carried into a pair's **first**
guarded iteration is **1.6e-25**, against a `SUM1G` accumulator of exactly
`2.0` — a relative perturbation of ~1e-25, far below double precision. The
dumped `DSZ` is **bit-identical**, 0 of 2401 elements differing.

> **Correction to the plan.** §2.5 says the error *"carries the full magnitude of
> the injected stale term."* Right about the mechanism, wrong about the realized
> size. **B6's conclusions stand and its validation did not need redoing.**
> A per-pair measurement is required to see this; a global one misleads
> (per-pair first trips run G = 2…26, and O(1) stale values occur only at later
> trips within a pair, where the correct contribution is ~1e-55 anyway).

- **Pinned by:** `test_madl2d_exh_guard.f90` (`UnitMadl2dExhGuard`), which also
  inspects `charge.f90` directly via a CMake `MADL2D_SOURCE` define, verified to
  fail when the old form is reinstated.

---

## C7 — The canonical charge variable → also settles **G-B7-4**

**Canonical definition: the deviation from the region reference,**

```
δq_i = q_i − q_bulk(region(i), type(i))
```

i.e. `imppot`'s definition, generalized from a single host to the region
registry. This is what makes the two-sided sums truncate on both sides and be
absolutely convergent, and what makes non-homogeneous bulk (polar hosts, ordered
alloys) correct.

| routine | current variable | meaning | disposition |
|---|---|---|---|
| `bulkpot` | `dq` | deviation from the neutral atom | leave alone; correct in its own context |
| `imppot` | `dq − bulk_charge(chargetrf_type)` | deviation from the bulk reference | **canonical** |
| `surfpot` | `dq` raw, summed over active layers | neither | outlier; leave alone, documented |
| `interfacepot` (B7.3) | — | must use the canonical definition | **non-negotiable** |

**Recommendation:** do **not** migrate `bulkpot` and `surfpot`. They are correct
for their own geometries, and migrating them would churn the bit-level
regression contract for no physics gain. Document the difference — this table
is that documentation. New interface code uses the canonical definition only.

---

## C8 — `impmad` and `surfmat` are in different unit conventions

**Directive: `interfacemat` mimics `surfmat`, never `impmad`.**

The two builders differ by **exactly one factor of `S`**. Each is
self-consistent with its own contraction, so **neither is a bug in place** —
but mixing them silently produces a plausible wrong answer.

| builder | units | on-site term | consumed by |
|---|---|---|---|
| `impmad` | **absolute** (1/bohr) | `2/wsimp(ii)` = `2/w` | `imppot`, contracting `this%amad` **bare** — no `/wsms` |
| `surfmat` | **the `1/S` convention** | `2*(sws*alat*ang2au/wssurf(i))` = `2·(S/w_i)`, dimensionless | `surfpot`, which restores Rydberg via `/wsms` |

For fccCu001 (`S = 2.6690 bohr`, uniform `w`): `2/w = 0.749346` versus
`2·(S/w) = 2.000000` — a ratio of exactly `S`.

Two further reasons `surfmat` is the right template for interface work:

- `impmad` sets `wsimp(:) = lattice%wav`, a single system-wide average, with
  the in-code comment *"Set as if all the atoms have the same WS radius. Can be
  improved later."* `surfmat`'s `wssurf` is genuinely **per-site**, which B7
  needs across two regions.
- `impmad` builds a **3D real-space cluster sum**; the interface geometry needs
  the **2D layer kernel** that `surfmat`/`madl2d` provide.

- **Flagged by:** a visible `@warning` block on `impmad` in `charge.f90`.
- **Action:** none to either routine. B7.3's `interfacemat` follows `surfmat`.

---

## C9 — `imppot`'s Madelung mixing *(§2.2 — was dead, now live)*

`vmad0(iclas)` was read *after* the contraction loop had already overwritten
`vmad`, so the mix reduced to `x*vmix + x*(1-vmix) == x` and `verr` was
identically zero. `bulkpot` (which captures `VMAD0` in a separate prior loop)
and `surfpot` were both already correct.

**Fixed in `eeecae9`** — `vmad0` is now captured before the overwrite, mirroring
`bulkpot`. Maintainer-approved on the grounds that the sign-off's concern was
convention and charge-scaling drift, not mixing: this changes **neither a matrix
element nor a scaling convention**, and at the default `vmix = 1.0` the
arithmetic is bit-identical, which is why the impurity regression is unchanged.

B7.4's alignment solver needs this hook live to damp the capacitor-like soft
mode in the alignment.

---

## Sign-off

| gate | covers | status |
|---|---|---|
| **G-B7-1** | C0–C6, C8, C9 | x Anders |
| **G-B7-4** | C7 (canonical charge variable) | x Anders |

### Sign-off note:

Overall action: Do not edit or change conventions inside "working" routines, i.e. those established earlier. If any significant bugs or physical errors are spotted, report them but do not fix then silently. Only ensure that the comments properly describes the actual code. For 1/S factors, make sure of conventions when "S" is relative radius (r_ws/r_avg) or actually has a simension. The W discussion for two different averages is unclear, even if we have two leads, we can still find the average W.

### Sign-off note and its consequences

The note above is the **governing rule** for the rest of B7. Applied:

1. **"Do not edit or change conventions inside working routines."**
   → C5's proposed `SWS` extraction is **withdrawn**. C2's proposed deletion of
   `surfpot`'s dead `twooverwsm` is **withdrawn**. No further edits to
   `bulkpot`, `imppot`, `surfpot`, `impmad`, `surfmat` or `madl2d` beyond
   comments — with the two exceptions already committed and separately
   justified below.

2. **"If any significant bugs or physical errors are spotted, report them but do
   not fix them silently."**
   → Reported, not silently fixed: the `impmad` `wsimp(:) = lattice%wav`
   coarseness (C2), and the shared l ≥ 1 moment normalization (C5). Neither is
   being changed. **Open items are listed below for your decision.**

3. **"Only ensure that the comments properly describe the actual code."**
   → Done for `DZZ` (dead assignment, removed with an explanatory comment) and
   `PM` (documented as the applied-bias hook). Both landed in `a4ee0e7`.

4. **"For 1/S factors, make sure of conventions when S is relative radius or
   actually has a dimension."**
   → Resolved in **C0**: `charge%sws` is **relative** (units of `alat`);
   `lattice%wav` and `wssurf(i)` are **dimensional**. Verified numerically for
   fccCu001: `sws ≈ 0.391`, `sws·alat·ang2au = 2.66899 bohr = wav·ang2au`. This
   **corrected C2**, which had wrongly stated the on-site term as `2δq/w` in
   Rydberg; it is `2·(S/w_i)`, dimensionless, and exactly `2` for uniform `w`.

5. **"The W discussion is unclear; even with two leads we can still find the
   average W."**
   → **Accepted; C5 rewritten.** `S` is a unit scale that cancels at the call
   site, not a local property, so a second region does not create a second unit
   scale. The per-site radius already enters separately via `wssurf(i)`. The
   residual question — whether one shared multipole normalization is a good
   compromise when radii differ a lot — is reclassified as a **physics-accuracy
   question for B7.7**, not a bug.

### Two code changes already committed

Both predate this note; flagging them explicitly rather than leaving them
implicit:

- **`6e87860`** — the `exh = -exh` → `exh = -exf` fix (C6). A genuine
  use-of-stale-value bug, but **measured to change nothing** (`DSZ`
  bit-identical). If you would rather revert it and carry it as a reported-only
  item, it is a clean single-line revert.
- **`a4ee0e7`** — dead `DZZ` assignment removed, `PM` documented. The `DZZ`
  removal is behaviour-neutral (the value was unconditionally overwritten two
  loops later); everything else in that commit is comments, docs and tests.

### Open items awaiting your decision

| item | status | note |
| --- | --- | --- |
| `imppot` dead Madelung mixing (§2.2) | **RESOLVED — reinstated** (`eeecae9`) | Maintainer approved: the concern was convention/scaling drift, not mixing. Bit-identical at the default `vmix = 1.0`; B7.4 needs the hook live |
| `impmad` vs `surfmat` conventions | **RESOLVED — flagged** (`eeecae9`) | Maintainer directive: mimic `surfmat` for `interfacemat`. Visible `@warning` added to `impmad`; see C8 |
| `surfpot` rewiring | **RESOLVED — keep as-is** | Maintainer confirmed: `surfpot` stays untouched; all new work wires into `interfacepot` |
| Shared l ≥ 1 moment normalization | reported only | Revisit in B7.7 if a two-region case shows dipole error |

**B7.0 is complete.** Item 5 (`SWS` extraction) is withdrawn (C5); item 3
(`imppot` mixing) is done (C9); items 1, 2, 4 and 7 were completed or resolved
as no-ops. B7.1 (region registry) has also landed. The next unblocked tasks are
**B7.2** (vacuum parameter generator) and **B7.3** (`interfacepot` /
`interfacemat`, following `surfmat` per C8).
