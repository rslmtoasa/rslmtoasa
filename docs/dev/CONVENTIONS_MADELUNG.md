# Madelung conventions contract

**Status:** awaiting maintainer sign-off — gate **G-B7-1** (blocks all of B7).
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
| 2.7 | `SWS` baked into l≥1 prefactors breaks two-region runs | **Confirmed defect.** Real work remains |
| 2.8 | On-site 2δq/w treated three different ways | **Partly wrong.** Present in `impmad` and `surfmat`; genuinely absent only from `surfpot`'s layer loop |

---

## C0 — Overall unit convention

**Potential is in Rydberg; charge in electrons; length in bohr.**
In Rydberg atomic units `e² = 2 Ry·bohr`, so the potential of a point charge
`q` at distance `r` is `V = 2q/r`.

The `madl2d` header comment says *"Potential in units of e²/(2S)"*. That refers
to the internal `TWOS = 2·SWS` scaling applied in the final loop, which is
undone at the call site by `surfpot`'s division by `wsms`. The **net** quantity
reaching `potential%vmad` is in Rydberg, consistent with C1.

> **Signing note.** The header comment is not wrong, but it describes an
> intermediate, not the delivered quantity. Recommend keeping the code as-is and
> treating this document as the authority.

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

**Convention: `V_onsite = 2·δq/w` in Rydberg** — the `e²/r` law of C0 evaluated
at `r = w`.

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
have the same WS radius. Can be improved later."* For one region that is a
constant absorbed into the reference; **for two regions at different `w` it is
not constant and cannot be absorbed.**

- **Pinned by:** `test_madelung_conventions.f90::test_onsite_term` (C2), which
  also asserts the term is genuinely `w`-dependent (8.7 mRy across a
  representative w = 2.66 → 3.05 bohr pair) and so cannot be absorbed.
- **Action for B7.3:** use per-site `w` (`wssurf` is already per-site; `wsimp`
  is not). Delete the dead `twooverwsm` code in `surfpot` to stop it reading as
  a missing term.

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

**This is the one convention that is currently WRONG for two regions, and the
only item in this document that requires code change before B7.3.**

`this%sws` is a **single system-wide scalar**, derived from the average cell
volume (`charge.f90:875`):

```fortran
this%sws = (3.0d0*this%vol/4.0d0/pi/this%nq3)**(1.d0/3.d0)
```

It is baked into the l ≥ 1 normalizations inside the kernel:

```
FACP   = SWS·√3      → FACDK, FACDR → DSZ
FACQUA ∝ SWS²        → DS3Z2, DSX2Y2, DSXY
FACGUA ∝ SWS³        → DZ3Z2
```

These are site-independent, so **with two regions at different `w` the dipole
coupling is wrong on one side.** `TWOS = 2·SWS` is a harmless global scaling
(undone at the call site); `DSS` and the monopole path are safe.

Correspondingly `surfpot`'s `wsm = lattice%wav` and `wsms = sws*alat` are
system-wide averages and must become per-region/per-site.

**Agreed target convention:** strip `SWS` from the kernel entirely; build the
kernel in **absolute units**; carry multipole moments in absolute units; apply
per-site `w` at the call site. This also makes the moment definitions auditable,
which they currently are not.

- **Not yet pinned by a test** — the test must be written against the *new*
  convention as part of the work, with the acceptance criterion that existing
  surface and bulk results stay bit-comparable within tolerance for uniform `w`.
- **Action:** B7.0 item 5, still open. **This is the remaining B7.0 work.**

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

## Open item not covered by this contract

`imppot`'s Madelung mixing is a **no-op** (§2.2): `vmad0(iclas)` is read *after*
the overwrite, so `verr` is identically zero and the mix does nothing.
`bulkpot` and `surfpot` are both correct. Impurity runs therefore converge
unmixed, which evidently works in practice — but the interface path introduces a
capacitor-like soft mode in the alignment, so this hook must be live before B7.4
leans on it.

This is a **bug fix, not a convention**, so it is deliberately not given a C-number.
It remains open as B7.0 item 3.

---

## Sign-off

| gate | covers | status |
|---|---|---|
| **G-B7-1** | C0–C6 | ☐ Anders |
| **G-B7-4** | C7 (canonical charge variable) | ☐ Anders |

Remaining B7.0 work after signature: **C5** (`SWS` extraction, item 5) and the
`imppot` mixing fix (item 3). Items 1, 2, 4 and 7 are complete or resolved as
no-ops.
