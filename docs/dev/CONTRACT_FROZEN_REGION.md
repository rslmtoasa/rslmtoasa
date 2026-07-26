# Frozen-region parameter contract — gate G-B7-2

**Status:** **AWAITING SIGNATURE** (Anders). Drafted 2026-07-26 as part of B7.4.
**Scope:** what a converged calculation must persist for its parameters to be
usable as a *frozen region* in an interface run (B7 §1.1, §1.3, §3 G-B7-2).
**Backed by:** `tests/unit/test_alignment_solver.f90` (`UnitAlignmentSolver`)
for the alignment machinery; the persistence facts below are verified against
`source/potential.f90` and recorded in `CONVENTIONS_MADELUNG.md` C0–C9.

> B7 §3 requires an agent reaching a gate to **stop and escalate** rather than
> choose. This document is the escalation: the findings are stated, a default is
> proposed for each open decision, and nothing in the code depends on the
> unsigned parts. **§4 is the part that needs your decision.**

---

## Why this gate exists

An interface run freezes two independently converged parameter sets side by
side. Each carries its own potential zero. If the two zeros disagree and nothing
reconciles them, the run imposes

> contact potential ≡ 0

whether or not that is true. It converges, it looks entirely plausible, and it
reports a wrong dipole barrier and wrong charge transfer.

B7 §3 calls a silent misalignment *"the worst failure mode available to us"*,
and that is not rhetoric: every other failure mode in this work package is
loud. This gate is the written statement of what has to travel with a parameter
set so the reconciliation is possible at all.

---

## 1. The load-bearing finding: `vmad` **is** persisted

B7 §2.1 flagged this as "verify first, it is load-bearing", with the worry that
`imppot`'s per-class add-back

```fortran
this%symbolic_atom(nbulk + jbas)%potential%vmad = ss &
   + this%symbolic_atom(this%lattice%chargetrf_type(jbas))%potential%vmad
```

might be adding a zero — in which case polar hosts would be silently treated as
non-polar and a new persisted field would be *"the single most important
deliverable of B7.0"*.

**It is persisted.** `potential%vmad` is read from the namelist
(`source/potential.f90:403`) and written back (`source/potential.f90:976`), and
is counted in the serialization size (`:1182`). So:

- cross-run alignment already has a home;
- **no new persisted field is required**;
- B7 §3.2's "assertion" reading is the correct one.

This is the same conclusion B7.0 reached (`CONVENTIONS_MADELUNG.md`, summary
row for §2.1); it is restated here because G-B7-2 is the gate that owns it.

---

## 2. What a frozen region must supply

| quantity | where it lives | required? | used for |
|---|---|---|---|
| C, Δ, γ/o per l and spin | `potential` | **yes** | the frozen Hamiltonian — already required today |
| Z | `element` | **yes** | as today |
| `vmad` | `potential%vmad` | **yes** | the region's on-site reference; the absolute potential zero |
| `q_bulk` per site type | `charge%bulk_charge`, from `get_charge_transf` | **yes** | deviation variables (§1.4) — without it every charge is absolute and the two-sided sums do not truncate |
| `w` (Wigner–Seitz radius) | `potential%ws_r` / `charge%wssurf` | **yes** | per-site; the overlap diagnostic (§2.9) reports when regions disagree |
| `E_F` of the source run | `energy%fermi` | **optional** | analytic initial guess and the consistency check *only* — never the determination |
| Madelung convention | implicit | **assumed** | see §3 |

The one entry that deserves comment is `E_F`. It is deliberately **optional**,
and the code treats it that way: `region_descriptor%fermi_known` gates both its
uses, and when it is absent the solver starts from the gauge and converges
anyway, just more slowly. That is the whole design stance of §1.3 — the
cross-run absolute scale is a *diagnostic*, not an input.

---

## 3. The parts that are assumed, not checked

Stated plainly so they are not discovered later:

1. **Both parameter sets use the same Madelung convention** (the one pinned in
   `CONVENTIONS_MADELUNG.md` C0–C9). Nothing in the file format records which
   convention produced a number, so this cannot currently be verified — only
   assumed. A set produced by a build predating C0–C9 would be silently
   misread.
2. **Both were converged with the same `alat` units and the same relative-`S`
   convention.** See `CONVENTIONS_MADELUNG.md` C0: `charge%sws` is dimensionless
   (units of `alat`) while `lattice%wav` and `wssurf(i)` are dimensional. A set
   that got this backwards is off by a factor of `alat` and would look like a
   large but plausible contact potential.
3. **The frozen sets are collinear.** Noncollinear leads need spin-rotated
   frozen parameters and are explicitly out of scope (B7 §6).

Item 1 is the reason §4's first decision exists.

---

## 4. Open decisions — **this is what needs your signature**

### 4.1 Should the interface reader refuse to run on a pre-contract parameter set?

B7 §3 G-B7-2 says the reader *"must **refuse to run** on a parameter set that
predates the contract."* Implementing that literally requires a version stamp
in the persisted potential, which does not exist today.

| option | cost | consequence |
|---|---|---|
| **(a) Add a contract-version stamp** to the persisted potential; refuse on absent or older. | a new persisted field + reader/writer; every existing `*_out.nml` becomes "pre-contract" and must be regenerated | matches the plan literally; loudest possible failure |
| **(b) Rely on the consistency check alone** (implemented): warn loudly when the converged fixed point disagrees with the analytic contact potential beyond `alignment_check_tol`. | none — already in `align_regions` | catches a *mismatched zero* regardless of provenance, but only when `E_F` is supplied for both regions, and warns rather than refuses |
| **(c) Both.** | as (a) | belt and braces |

**Proposed default: (b) for now, (a) when B7.5 lands.** The reasoning: (b)
catches the actual failure — two sets on different zeros — rather than a proxy
for it, and it costs nothing. But (b) is silent when `E_F` is not supplied, and
the plan's insistence on refusal is well founded given the failure mode. B7.5
introduces the `buildinterface` namelist and is the natural place to add a
stamp, since it is already defining the region-parameter-path syntax.

**This is your call, not mine.** Nothing currently refuses to run.

### 4.2 Is `alignment_check_tol = 5 mRy` the right alarm threshold?

The consistency check warns when the converged fixed point differs from
`E_F(anchor) − E_F(region)` by more than this. It is meant as a *"something is
structurally wrong"* alarm, not a precision test, so it is deliberately
generous — a genuine contact potential is typically tenths of a Ry, and the
disagreements this is designed to catch are of that order.

Set too tight, it cries wolf on ordinary convergence noise; too loose, it misses
a real misalignment. 5 mRy is a first guess and has **not** been calibrated
against a real two-region run, because none exists until B7.5/B7.6.

**Proposed: leave at 5 mRy, revisit in B7.7** with the first A|B benchmark, in
the same pass that feeds G-B7-3.

### 4.3 Is `deep_drift_tol = 1 mRy`-equivalent the right buffer-thickness alarm?

Same situation. The deep-A drift diagnostic reports the largest deviation charge
on a frozen anchor-region site; above `deep_drift_tol = 1e-3` electrons it warns
that the active zone is too thin. Uncalibrated for the same reason.

**Proposed: leave, revisit in B7.7** alongside the buffer-thickness convergence
data that G-B7-3 already requires.

---

## 5. What B7.4 implemented, for the record

Signed or not, these are pinned by `UnitAlignmentSolver` and will not be
re-derived in a later session:

- **The gauge anchor is the first frozen, non-vacuum region.** Vacuum carries no
  states at E_F, so the "deep region is neutral bulk" residual cannot be
  evaluated there; an active region has no frozen reference to *be*. For the
  buildsurf registry order (vacuum, active, bulk) both traps are live, and the
  anchor comes out as bulk.
- **The anchor is re-pinned to exactly zero every update**, not merely skipped,
  so accumulated round-off cannot let the gauge drift.
- **The fixed point is the answer; the analytic value is the check.** Never the
  other way round (§1.3).
- **`fix_fermi_to_region` makes the named region the anchor.** "E_F is region
  r's Fermi level" and "region r carries no alignment shift" are the same
  statement, since E_F − V_r = E_F^(r); anchoring elsewhere would double-count
  the pin. Default is free E_F from cluster neutrality, per the maintainer
  decision recorded in §1.3.

---

## Sign-off

| gate | covers | status |
|---|---|---|
| **G-B7-2** | §1 (`vmad` persistence), §2 (required quantities), §4 (open decisions) | ☐ Anders |

### Sign-off note:

<!-- Anders: decisions on §4.1, §4.2, §4.3 here. -->
