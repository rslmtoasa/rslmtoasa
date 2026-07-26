# B7 — Interfaces and vacuum leads in the embedding formalism
### (replaces the previous B7 "Boundary self-energies: interface leads + vacuum GF")

**Effort:** M–L. **Leads:** OPUS on conventions, electrostatics, alignment and the
vacuum generator; SONNET on registry bookkeeping, namelist and examples.
**Depends on:** B6 (surface dipole electrostatics) — landed, but see B7.0 for a
correctness fix that invalidates part of its validation.
**Does not depend on:** B2 (k-space backends), B4 (BZ primitives).
**Hands off to:** B9 (CPA) inherits nothing from this file; the per-energy
recursion driver formerly specified here belongs there.

---

## 0. Scope, and the methodological stance that defines it

### 0.1 What this work package is

RS-LMTO-ASA does everything by **embedding**. A finite real-space cluster is
built; an *interaction zone* of atoms is allowed to relax self-consistently;
all atoms outside it are **frozen**, meaning their potential parameters
(C, Δ, γ/o, plus Z and the on-site Madelung shift `vmad`) are held at values
imported from a previous converged calculation and never updated. The Green
function is obtained by Haydock/Chebyshev recursion on this cluster. Every
existing workflow is a special case of that one idea:

| workflow | frozen medium | E_F |
|---|---|---|
| bulk / supercell | the cell repeated into itself | free, from charge neutrality |
| bulk impurity | monoatomic or supercell host | fixed, from the bulk run |
| pristine surface | one-sided bulk + hand-set empty spheres | fixed, from the bulk run |
| surface impurity | pristine-surface parameters | fixed, from the surface run |

The purpose of this WP is to add **one more region label**. A "lead" — left,
right, or vacuum — is nothing but *a set of frozen potential parameters
assigned to a subset of cluster sites*. There is no new Green-function
machinery. A "vacuum Green function" is an ordinary bulk Green function whose
frozen parameters happen to describe an empty lattice.

What must actually be built is therefore narrow and entirely classical:

1. bookkeeping that lets **more than one** frozen reference region coexist in
   one cluster (today the reference is implicitly unique);
2. a generator for the frozen parameter set of an **empty-sphere lattice**
   (the vacuum lead), replacing today's hand-set empty-sphere pragmatics;
3. a **two-sided** charge-transfer and electrostatic treatment, in deviation
   variables, that handles the interface dipole and the relative alignment of
   the two frozen regions correctly.

### 0.2 What this work package is explicitly NOT

These exclusions are not stylistic. Delegated agents must not reintroduce any
of them, and must escalate rather than "helpfully" adding one.

- **No decimation, no López-Sancho/Sancho-Rubio surface Green functions.**
- **No principal-layer partition.** The layer structure used here is the
  physical atomic-layer structure of the slab, not a principal-layer blocking.
- **No boundary self-energy Σ(E) object.** Nothing energy-dependent is added to
  the Hamiltonian. The recursion coefficients stay energy-independent and the
  continued fraction remains evaluable per energy at no cost.
- **No per-energy recursion driver.** That infrastructure is genuinely needed —
  for CPA and DMFT self-energies — and belongs to B9/B12. It must not be built
  here, because here nothing needs it.
- **No k-space / k∥ 2D route.** Real space only.
- **No transport.** This WP delivers interface *electronics and magnetics*:
  charge transfer, dipole barrier, layer-resolved moments and DOS, and
  interlayer exchange coupling via Lichtenstein-type formulas. It does **not**
  deliver Landauer conductance. A real-space Kubo–Greenwood/Chebyshev path is a
  separate future milestone and must not be implied by the word "lead".
- **No solution of lattice mismatch.** The cluster builder is assumed to
  produce correct structure constants for whatever geometry the user supplies.
  Mismatch is the user's problem; the only obligation here is a *diagnostic*
  (§2.9). Regular, commensurate geometries are assumed throughout v1.

### 0.3 Why "lead" is still the right word

The L/R regions are semi-infinite reservoirs of fixed potential and fixed
Fermi level, exactly as in the layered-KKR/LMTO literature (Skriver &
Rosengaard; Abrikosov; Turek et al.). We reproduce their *physics* by a
different *numerical route*. Keeping the name preserves the connection to that
literature and to the applied-bias extension anticipated in §6.

---

## 1. Physics specification

Agents implementing §4 must read this section in full. Several of the
statements below look like conventions but are correctness conditions.

### 1.1 Regions, and what "frozen" means

The cluster is partitioned into **regions** r = 1…N_reg. Each region carries:

- a frozen potential-parameter set (C_l, Δ_l, γ_l/o_l per l and spin, plus Z),
- a reference charge per site type, `q_bulk(r, type)`,
- a reference on-site Madelung potential per site type, `vmad_bulk(r, type)`,
- a **region alignment shift V_r** (§1.3), a single scalar,
- a Wigner–Seitz radius w_r (may differ between regions).

Sites carry an **active/frozen mask**. Active sites relax self-consistently;
frozen sites never update their parameters. The user chooses the width of the
active zone, exactly as they choose single-site vs extended impurity clusters
today. There is no automatic determination and none should be added.

Frozen sites may legitimately be charged: a non-homogeneous supercell bulk
(FeCo, NiO) has site-dependent charge transfer baked into its converged
parameters. This is why every charge entering the electrostatics must be a
**deviation from the region reference**, never an absolute charge (§1.4).

### 1.2 Geometries covered

- **A | vacuum** — a surface. One metallic boundary, one vacuum boundary.
  Equivalent in output to today's `buildsurf` path, but with the vacuum
  described by generated empty-lattice parameters rather than hand-set ones.
- **A | B** — a metallic interface. Two metallic boundaries.
- **A | vacuum-gap | B** — two metallic boundaries with an *active* vacuum
  region between them. This is the tunneling / vacuum-gap IEC geometry and it
  requires no machinery beyond A | B.
- **A | spacer | A** — the interlayer-exchange geometry; three regions, of
  which the spacer may be active throughout.

All are the same code path. The differences are entries in the region
registry, not branches.

### 1.3 Alignment: the one genuinely new unknown

For a surface, E_F is pinned to the bulk value and there is nothing to align.
For A | B this fails, and the failure is silent, which makes it the highest
correctness risk in this WP.

Each bulk calculation carries its own potential zero. Region A's frozen set
implicitly asserts one absolute energy scale; region B's asserts another.
Freezing both at their independently converged values, with no relative shift,
imposes

> contact potential ≡ 0

whether or not that is true. The run converges, looks plausible, and reports a
wrong dipole barrier and wrong charge transfer.

**The physical requirement.** Deep inside region r, the site must *be* neutral
bulk r. On a common absolute scale that is

    E_F − V_r = E_F^(r)

so the shifts are fixed up to one overall gauge, and only differences are
physical:

    V_B − V_A = E_F^(A) − E_F^(B)      (the contact potential)

**How V_r enters the code.** Not through C_l and E_ν. It enters through
`vmad`, which is already treated as a potential parameter and is already added
back per class in `imppot`:

```fortran
this%symbolic_atom(nbulk + jbas)%potential%vmad = ss &
   + this%symbolic_atom(this%lattice%chargetrf_type(jbas))%potential%vmad
```

That line is "deviation potential + reference-region on-site potential". The
interface path generalizes the reference from a single host to the region
registry and adds V_r:

    vmad(i) = dV(i) + vmad_bulk(region(i), type(i)) + V_region(i)

with V_r constant over the region, applied identically to frozen and active
sites of that region. Frozen sites remain frozen — frozen *at a shifted
value*. This is a change of reference level, not a relaxation.

**How V_r is determined.** *Not* from the closed form V_r = E_F − E_F^(r),
because that requires a reliable cross-run absolute energy scale which we do
not currently guarantee (§2.1). Instead, solve it as an SCF fixed point on a
quantity `surfpot` already computes.

`surfpot` evaluates the Madelung potential at the deep-bulk row and uses it as
the reference `vbulk` that drives dV → 0 in the interior. Two-sided, evaluate
the same probe on both sides. Then

    V_B(out) = dV(deep-B) − dV(deep-A)

is the residual, mixed with `vmix` like any other SCF quantity, with region A
as gauge anchor (V_A ≡ 0). At convergence, deep-B sees the potential it needs
to be bulk B. No cross-run reference is required.

The analytic value E_F^(A) − E_F^(B), when the absolute-zero contract (§3.2)
is satisfied, is then used as (i) the initial guess, which is usually good
enough to save several iterations, and (ii) a **consistency check**: if the
converged fixed point disagrees with the analytic value beyond a threshold,
the absolute-zero bookkeeping is broken and the run must warn loudly. This
converts our most fragile ingredient from load-bearing into a diagnostic.

**Fermi level.** Per maintainer decision, the default is a **free E_F**
determined by charge neutrality of the cluster, as in supercell mode, combined
with the internal gauge fix V_A ≡ 0. The two unknowns (E_F, V_B) are then
determined by two conditions (cluster neutrality, deep-B bulk residual) and
there is no flat direction.

A caveat that must appear in the code as a diagnostic: with free E_F, "deep-A
is neutral bulk" is *not imposed*, it is a consequence of the buffer being
thick enough. Report δq at the innermost frozen A sites every iteration; drift
there means the active zone is too thin, and it means the gauge is
inconsistent.

A `fix_fermi_to_region = A` option must also exist. It reduces exactly to the
present surface behaviour and is the correct setting when reproducing
`buildsurf` results.

### 1.4 Electrostatics in deviation variables

Everything the 2D Madelung sum sees must be a deviation from the local region
reference:

    δq_i = q_i − q_bulk(region(i), type(i))
    δp_i = p_i − p_bulk(region(i), type(i))     (intra-sphere l = 1 moment)

Consequences, all of them load-bearing:

- frozen deep sites contribute exactly zero, so the sums truncate on **both**
  sides and are absolutely convergent;
- non-homogeneous bulk (polar hosts, ordered alloys) is handled correctly,
  because their built-in charge transfer is in the reference, not in δq;
- the two physically meaningful quantities are the first two moments of the
  deviation profile:

      Q = Σ_i δq_i                    → must vanish (no residual field)
      P = Σ_i δq_i z_i + Σ_i δp_i     → the potential step across the interface

The step is what appears as the work function (metal|vacuum) or the interface
dipole barrier (A|B). In `surfpot` today, `vmad1 = vm1 − vbulk` is exactly
this quantity for the one-sided case: `vm1` probes the deep-vacuum row and
`vbulk` the deep-bulk row.

**Q = 0 is a precondition on kernel validity, not merely on physics.** See
§2.6: the existing `AM` matrix is built in a half-space gauge whose difference
from the symmetric kernel is a *global constant* if and only if Q = 0. If
Q ≠ 0, the gauge choice injects a spurious uniform field of magnitude
(2π/A)·Q on top of the genuine divergence of a charged slab.

### 1.5 The compensating charge

`Σ_active δq_i` will not be exactly zero: the active zone screens imperfectly.
The residual must be compensated by charge placed on the frozen side(s), and
the question is *where*.

**Reformulate before answering.** The compensation profile has a monopole and
a dipole. The monopole is fixed by neutrality. The dipole,
Σ q_comp z, is the only part that touches an observable — and its lever arm is
the full slab thickness. So an arbitrary split fraction between the two sides
would be multiplied by the largest length in the problem, landing directly on
the work function.

**The rule.** Weight the compensation by the side-resolved N(E_F) of the
boundary buffer layers, which recursion provides for free:

| boundary type | N(E_F) | compensation |
|---|---|---|
| metallic lead | finite | receives charge ∝ N(E_F) |
| vacuum / insulating | ≈ 0 | receives **nothing** |

This is the linear-response answer for partitioning screening charge between
two reservoirs. It gives ≈ 50/50 for two similar metals and collapses to
"everything into the metal" for metal|vacuum — so **the surface case needs no
separate branch**. That is not a coincidence: it is the physics. A nonzero
residual means the active zone failed to screen; the remaining screening
happens where there are states to do it; vacuum has none.

Two clarifications that agents get wrong if not told:

1. **Real spill-out is not compensation.** Charge outside the surface plane is
   carried by *active* empty spheres as genuine self-consistent δq. If there
   is "not enough charge outside", the fix is more active empty spheres, never
   compensation placed in the vacuum. Placing compensation there does not
   perturb the work function — it *sets* it.
2. **Do not inherit `imppot`'s smear.** `imppot` spreads the excess as
   `tdq(j) = -dif/nsum` over every frozen site. In a 3D impurity cluster that
   is a roughly isotropic shell with Σ δq z = 0 by symmetry — harmless
   regardless of thickness. In layered geometry the same smear has a nonzero
   first moment and a long lever arm. Same code, different geometry, different
   consequence. Compensation here must be **localized** at defined z_L / z_R
   (or a narrow, explicitly specified weight profile within the screening
   length), as `surfpot` already does one-sided with `tdq(iex) = -dif`.

**Acceptance criterion, not argument.** Specify the compensation as a
normalized weight profile w_i and require the reported step to be insensitive
to reasonable variation of w_i within the screening region. Sensitivity scales
as |Σ_active δq| × lever arm, so this is really a statement about the active
zone: thicken it until the profile choice stops mattering. Cheap to test.

### 1.6 The vacuum lead

The frozen parameter set for vacuum is the fixed point of an **infinite empty
lattice**: an empty sphere of radius w containing a constant potential,
rigidly shifted to the vacuum level.

**Generate them with our own radial solver, V(r) = const**, rather than
hand-coding spherical-Bessel expressions. The representation, γ/o convention,
E_ν convention and normalization then come out automatically consistent with
every other parameter set in the code. The analytic free-electron result
(logarithmic derivatives of j_l(κr), κ = √(E − V_0)) is the **test oracle**,
not the implementation.

The generated set is geometry-dependent (w, 2D lattice, layer spacing), so
this is a small generator utility invoked per run, **not a table of
constants**.

**Honest limitation, to be stated in the documentation and not discovered
later.** LMTO-ASA reproduces empty-lattice free-electron dispersion only over
a limited window above the band bottom; this is precisely where the combined
correction earns its keep. The frozen vacuum parameters are meaningful only if
E_F sits well below the vacuum band onset. True for ordinary metals;
uncomfortable for very low work functions. Warn if the margin is small.

**Set expectations on the payoff.** For a good metal the evanescent decay
constant is κ = √Φ ≈ 0.54 bohr⁻¹ for Φ ≈ 4 eV, so tails die within ~2 Å.
*That is why one or two empty-sphere layers already works.* The vacuum lead
will not change converged surface numbers much, and the acceptance criterion
must be **reduced sensitivity to the empty-sphere layer count**, not a changed
answer. The real payoff is elsewhere: a well-defined vacuum level, removal of
an ad hoc convergence knob, and geometries where the tail genuinely matters —
vacuum-gap IEC, tunneling barriers, low-DOS spacers.

### 1.7 Expected physics limitations

Long-period interlayer-exchange oscillations arise from fine Fermi-surface
features, and real-space recursion resolves those only with long chains.
Fe/Cr and Co/Cu spacer-thickness benchmarks are the right validation targets,
but budget for an **LL-convergence study on the oscillation period**, not only
on the layer moments. Do not treat a period discrepancy as a bug before that
study exists.

---

## 2. Findings from the audit of the existing potential path

These are concrete, verified against the sources for `bulkpot`, `imppot`,
`surfpot` and `madl2d`. They are the reason B7.0 exists as a separate task.
Items 2.1–2.3 and 2.5 are bugs or latent bugs **independent of this WP**.

### 2.1 `vmad` persistence across runs — verify first, it is load-bearing

`imppot` adds the host's converged on-site `vmad` back per class. That is the
structurally correct treatment of an impurity in a polar host (NiO, ordered
alloys) — *provided the host `symbolic_atom` entries actually carry their
converged `vmad` at load time*. If `vmad` is not part of the persisted
potential-parameter set, that line adds zero, polar hosts are silently treated
as non-polar, and everything still converges.

Verify this before anything else. It also resolves the absolute-zero question
for the whole WP: **`vmad` is the absolute reference.** If persisted,
cross-run alignment already has a home and §3.2 is an assertion; if not, a new
persisted field is required and it is the single most important deliverable of
B7.0.

### 2.2 `imppot`'s Madelung mixing is a no-op

```fortran
vmad0(iclas) = ...%potential%vmad          ! read AFTER the overwrite
...%potential%vmad = ...%vmad*vmix + vmad0(iclas)*(1 - vmix)
verr = ...%vmad - vmad0(iclas)             ! identically zero
```

`bulkpot` captures `VMAD0` in a separate loop *before* touching `vmad` and is
correct; `surfpot` reads the stored previous value and is correct. `imppot` is
the only broken one. Impurity runs therefore converge unmixed, which evidently
works in practice — but the interface path introduces a capacitor-like soft
mode in the alignment, so this hook must be live before we lean on it. Fix,
and confirm existing impurity regressions are unchanged at convergence.

### 2.3 Factor-of-two inconsistency in the Madelung contraction

`bulkpot` contracts `2.d0*AMAD(JBAS,IBAS)*TDQ(...)` reading `mad.mat` from
disk into a local array; `imppot` contracts `this%amad(jbas,n)` bare. Either
`this%amad` is built with the 2 folded in, or one of them is off by a factor
of two. Pin this with a single known-charge test rather than by reading
provenance. Do not assume.

### 2.4 Three different charge variables

| routine | variable | meaning |
|---|---|---|
| `bulkpot` | `dq` | deviation from the neutral atom |
| `imppot` | `dq − bulk_charge(chargetrf_type)` | deviation from the bulk reference |
| `surfpot` | `dq` (raw, summed over active layers) | **neither** |

`surfpot` is the outlier. For a non-homogeneous bulk with legitimately charged
sites, its `dif` never converges to a deviation and the excess charge is
contaminated. The interface path must take `imppot`'s definition. This is also
what makes the two-sided sums absolutely convergent (§1.4). Non-negotiable.

### 2.5 `madl2d`: `exh = -exh` — a real bug, fix independently of this WP

In the overflow guard:

```fortran
if (erfcp == 0) then
   expm = exp(-dgi*qppz)
   exf = expm*erfcm
   exh = -exh          ! WRONG — must be  exh = -exf
```

With `ERFCP == 0` we have `EXPP*ERFCP = 0`, hence `EXF = EXPM*ERFCM` (correct
in the code) and `EXH = 0 − EXPM*ERFCM = −EXF`. Writing `exh = -exh` reuses
whatever `exh` held from the previous reciprocal-lattice iteration — a stale,
possibly O(1) value substituted where a vanishing one belongs. The correct
contribution there is exponentially small, so the error carries the full
magnitude of the injected stale term. If the branch trips on the first G of
the first pair, `exh` is undefined.

**Scope of the damage.** `EXF` is correct in that branch, so `BM`, `DSS` and
the whole monopole path are clean. `EXH` feeds only
`DMG → SUM1G → DMDL → DSZ`. So this corrupts the **dipole** matrix
specifically — the B6 feature — and it triggers when
`BETAP = DGI/(2λ) + λ·QPPZ` becomes large, i.e. at **large interlayer
separation and in thick slabs**. That is exactly the regime interface clusters
live in, and exactly the regime in which a surface-dipole convergence study
would look mysteriously noisy in the layer count.

Fix, then **re-run the full B6 validation** before building on top of it.

### 2.6 `madl2d`: the `AM` half-space gauge

`BM` is symmetric. `DSZ` is antisymmetric with equal magnitude (±DMDL), which
is the physically correct ±2πp/A on the two sides of a dipole sheet and is
fine two-sided as written. But `AM`, the k∥ = 0 monopole term, is deliberately
asymmetric:

```fortran
AM(IQ,JQ) = FACGAU*EXPZ - QPPZ*FACERF*ERFCM    ! → −(4π/A)|Δz|
AM(JQ,IQ) = FACGAU*EXPZ + QPPZ*FACERF*ERFCP    ! → 0
```

The symmetric kernel would be −(2π/A)|Δz| in both directions. The code puts
the entire plate term on one triangle: a gauge choice that zeroes the field on
one side, i.e. a half-space convention.

The difference from the symmetric kernel, acting on a deviation profile, is
(up to the overall unit factor `TWOS`)

    D_i = (2π/A) · [ P − z_i · Q ],   P = Σ_j z_j δq_j,   Q = Σ_j δq_j

so **if Q = 0 exactly, D_i is a global constant**, removed by the reference
subtraction and cancelling identically in any dV(deep-B) − dV(deep-A)
difference. **The existing kernel is therefore reusable two-sided, unchanged.**
If Q ≠ 0 it acquires a z-dependence — a spurious uniform field.

This hands us the sharpest available oracle (§5.2): symmetrize as
½(AM(i,j) + AM(j,i)), recompute, and require *identical dV differences* for a
neutral profile, plus a difference of exactly the predicted (2π/A)Q form when
a deliberate net charge is injected. That pins the gauge and the 4π/sign
convention in one test.

### 2.7 `madl2d`: `SWS` is baked into the multipole prefactors

`TWOS = 2·SWS` is a global unit scaling and harmless. The l ≥ 1 normalizations
are not:

    FACP   = SWS·√3      → FACDK, FACDR → DSZ
    FACQUA ∝ SWS²        → DS3Z2, DSX2Y2, DSXY
    FACGUA ∝ SWS³        → DZ3Z2

These are site-independent. With two regions at different w the dipole
coupling is wrong on one side. `DSS` is safe.

**Fix:** strip `SWS` out of the kernel entirely; build the kernel in absolute
units; carry the multipole moments in absolute units; apply per-site w at the
call site. This also makes the moment definitions auditable, which they
currently are not.

Correspondingly, `surfpot`'s `wsm = lattice%wav` (a single system-wide average
WS radius) and `wsms = sws*alat` must become per-region/per-site.

### 2.8 The on-site 2δq/w term is treated three different ways

`bulkpot` adds `VADD = 2*TDQ/RMAX`. `imppot` adds no such term for the
impurity's own deviation. `surfpot` has `twooverwsm` present but commented out
throughout. So in `imppot`, the host reference being added back *includes* its
2δq/w while the impurity deviation on top of it does not — the two pieces sit
on different footings.

For uniform w this is a constant absorbed into the reference, which is
presumably why it was dropped in `surfpot`. **With two regions at different w
it is not constant and cannot be absorbed.** It must come back, per-site, with
the correct w.

### 2.9 Geometry hygiene — diagnostic only

With unequal w across the interface the ASA volume-filling condition cannot
hold on both sides. Per maintainer instruction, do not attempt to solve this.
Emit a per-layer report of Σ(sphere volume) vs cell volume and warn above an
overlap threshold. Cheap, and it converts a silent error into a message.

### 2.10 Smaller items

- **`madl2d` is fully z-general.** Everything runs off
  `QPPZ = QZ(IQ) − QZ(JQ)` for arbitrary pairs; there is no half-space
  truncation and no assumed layer spacing. Relaxed-z in the 2D Madelung is
  therefore **plumbing** (feed real `QZ`), not a rewrite. Out of scope here,
  but record it — it can be promoted cheaply whenever wanted.
- **`DZZ` has a dead assignment.** It is set in the pair loop from unscaled
  `Q0MDL` and then unconditionally overwritten in the `TWOS` loop from the
  already-scaled `DS3Z2`; its diagonal is only ever set by the second.
  Harmless while l = 2 is unconsumed; a trap for whoever turns it on.
- **`PM`, the plate-condenser matrix, is built and never used.**
  Lower-triangular, `FACDIF*(QZ(JQ) − QZ(IQ))`. This is precisely the
  applied-bias hook wanted in §6, already present. It carries the same
  triangular gauge choice as `AM` and needs the same symmetrization treatment
  when it goes live.
- **`init = 6` is a magic offset** encoding the count of leading frozen
  (vacuum-side) rows in the 2D slab; `iq = init + ibas` is the layer→matrix-row
  map. Replace with an explicit registry (§4, B7.1). Note the index map as it
  stands: rows 1…`init` are the vacuum-side frozen layers, rows
  `init+1`…`init+nlay` the active layers, and the bulk-side frozen rows run out
  to `nbas`; `vbulk = dss(nbas,·)` probes deep bulk, `vm1 = dss(1,·)` probes
  deep vacuum, and `vmad1 = vm1 − vbulk` is the surface dipole barrier. The
  charge loop fills only rows `init+1`…`nbas`, so the vacuum rows exist in the
  geometry but are held at exactly zero charge — which is the correct
  behaviour per §1.5 and should be preserved deliberately, not by accident.
- **`nbas` is the extended interaction zone**, user-sized. It is not a proxy
  for "the whole system", and nothing should assume a particular width.

---

## 3. Conventions and gates

Gates are maintainer sign-off points. An agent reaching a gate **stops and
escalates**; it does not choose, and it does not proceed on a plausible
assumption. Once signed, the convention is pinned by a unit test and never
re-derived in a later session.

### G-B7-1 — Madelung convention contract *(blocks everything)*
A `CONVENTIONS_MADELUNG.md` fixing, by numerical test rather than derivation:
the overall unit convention (the header claims e²/2S — confirm), the 4π vs 2π
sheet prefactor and its sign, the `AM` gauge decision (symmetrize vs keep
half-space with Q = 0 enforced), the factor-of-two question in §2.3, the
multipole moment normalization after `SWS` extraction, and the sign convention
for the l = 1 moment. Anders signs.

### G-B7-2 — Absolute-zero / `vmad` persistence contract
The answer to §2.1, written down: which quantities a converged bulk run must
persist for its parameters to be usable as a frozen region (potential zero
convention, w, Madelung convention, E_F). The interface reader must **refuse
to run** on a parameter set that predates the contract. Silent misalignment
produces a plausible but wrong dipole barrier — the worst failure mode
available to us. Anders signs.

> **Drafted, awaiting signature:** [`CONTRACT_FROZEN_REGION.md`](../CONTRACT_FROZEN_REGION.md)
> (B7.4). `vmad` **is** persisted, so no new field is required and §3.2 reads as
> an assertion. Three decisions are open and flagged there: whether the reader
> should refuse on a pre-contract set (needs a version stamp that does not exist
> today — proposed for B7.5), and the two uncalibrated alarm thresholds
> (proposed to revisit in B7.7 with the first A|B benchmark). Nothing currently
> refuses to run.

### G-B7-3 — Compensation weight profile
The default weight profile and its localization width (§1.5). Anders signs
after seeing the buffer-thickness convergence data, not before.

### G-B7-4 — Canonical charge variable
Which deviation definition becomes canonical across `bulkpot`/`imppot`/
`surfpot`/`interfacepot`, and whether the older routines are migrated or left
alone with a documented difference (§2.4). Anders signs.

---

## 4. Tasks

Every task ends with tests green and a short note of anything that surprised
the implementer. Agents **escalate on spec disagreement** rather than
resolving silently. No defensive scaffolding on hot paths.

### B7.0 [OPUS] — Audit, fixes and convention extraction across the potential path
*The prerequisite for everything, and independently valuable.*

1. Fix `exh = -exh` → `exh = -exf` (§2.5). Add a unit test that drives the
   branch (large `QPPZ`, large `DGI`) and compares `DSZ` against the
   unguarded high-precision evaluation. **Re-run the complete B6 validation**
   and report whether any B6 number moves; if it does, B6's conclusions need
   revisiting.
2. Verify `vmad` persistence (§2.1). Report the answer explicitly. If not
   persisted, add the field and the reader/writer, and construct a polar-host
   impurity regression (impurity in NiO) that fails without it.
3. Fix `imppot`'s dead mixing (§2.2); confirm impurity regressions unchanged
   at convergence.
4. Resolve the factor of two (§2.3) with a known-charge test.
5. Extract `SWS` from the `madl2d` kernel; move to absolute-unit moments with
   per-site w applied at the call site (§2.7). Regression: existing surface
   and bulk results bit-comparable within tolerance for uniform w.
6. Restore the on-site 2δq/w term consistently across all three routines
   (§2.8), per-site.
7. Remove the dead `DZZ` assignment; document `PM` as the bias hook (§2.10).
8. Produce `CONVENTIONS_MADELUNG.md` → **G-B7-1**, with each entry backed by
   a passing unit test.

*Kit:* this file; `charge` type sources containing `bulkpot`, `imppot`,
`surfpot`, `madl2d`; B6 validation suite; the existing surface and impurity
regression cases.

### B7.1 [SONNET] — Region registry and index map
Replace the `init`/`iq = init + ibas` offset arithmetic with an explicit
registry: per site → region ID, active/frozen mask, layer index, z, w,
reference type. z carried as **data**, so relaxed-z later is a parameter
change (§2.10). Registry construction from the cluster builder, plus a dump
routine for debugging. Reproduce today's `buildsurf` layout exactly as a
registry instance and show the surface regression is unchanged.

*Kit:* B7.0; `surfpot` index map as documented in §2.10; lattice/cluster
accessors.

### B7.2 [OPUS] — Vacuum parameter generator
Empty-lattice fixed-point parameters via the code's own radial solver with
V(r) = const (§1.6). Geometry-dependent, generated per run. Analytic
spherical-Bessel oracle as test. Warn when E_F is not comfortably below the
vacuum band onset. Self-contained; no interface integration.

*Kit:* B7.0; the radial solver entry points; the existing empty-sphere setup
as reference case.

### B7.3 [OPUS] — Two-sided deviation-variable electrostatics
`interfacepot`: `imppot`'s reference bookkeeping on `surfpot`'s 2D kernel.
Deviation variables throughout (§1.4, §2.4); per-region w; on-site term;
two-sided probes at deep-A and deep-B; the moments Q and P computed and
reported every iteration; enforcement and reporting of Q = 0 as a kernel
precondition (§2.6). Localized compensation per §1.5 with side-resolved N(E_F)
weights — **not** `imppot`'s smear. Sphere-overlap diagnostic (§2.9).

*Kit:* B7.0, B7.1; `imppot` and `surfpot` sources; G-B7-1 contract.

### B7.4 [OPUS] — Alignment solver
V_r via `vmad` (§1.3): registry-wide reference addition, V_A ≡ 0 gauge anchor,
V_B fixed point on the deep-probe residual with `vmix` mixing, free E_F from
cluster neutrality plus the `fix_fermi_to_region` option, analytic initial
guess and the consistency-check warning, deep-A δq drift diagnostic. Delivered
with the offset-δ oracle (§5.3) — which must be written **before** any
interface geometry exists, on a toy cluster if necessary.

*Kit:* B7.0, B7.3; G-B7-2 contract.

### B7.5 [SONNET] — `pre_processing = 'buildinterface'` and namelist
New path; `buildsurf` untouched and permanently retained as regression oracle.
Namelist: region definitions and parameter-set paths, active/frozen widths per
region, compensation profile, Fermi mode, bias placeholder (§6). Workflow
chain bulk-A → bulk-B (→ vacuum generator) → interface, mirroring
bulk → surf → imp.

*Kit:* B7.1–B7.4; existing `buildsurf` namelist as template.

### B7.6 [SONNET] — Example cases and documentation
A | A identity; A | vacuum (reproducing a `buildsurf` case); A | B metallic;
A | vacuum-gap | A. Documentation covering the limitations in §1.6 and §1.7
honestly — especially that the vacuum lead is not expected to change good-metal
surface numbers, and what it *is* for.

*Kit:* B7.1–B7.5.

### B7.7 [OPUS] — Physics validation campaign
The ladder in §5.5–§5.7, including the Fe/Cr and Co/Cu IEC benchmarks with the
LL-convergence study on the oscillation period (§1.7).

*Kit:* all of the above; existing IEC post-processing.

---

## 5. Validation ladder

Ordered by increasing cost. Each rung is a regression test, not a one-off.

### 5.1 Known-charge Madelung test
A single known charge distribution in a fixed geometry, checked against a
directly summed reference. Pins the factor of two (§2.3) and the unit
convention.

### 5.2 `AM` symmetrization oracle
Symmetrize `AM` as ½(AM(i,j) + AM(j,i)). For a **neutral** deviation profile,
require *identical* dV differences from the symmetric and asymmetric kernels.
Then inject a deliberate net charge Q ≠ 0 and require the difference to equal
(2π/A)(P − z_i Q) exactly. Pins the gauge, the 4π and the sign in one test
(§2.6).

### 5.3 A | A with deliberate rigid offset δ — the primary oracle
**A | A alone does not test alignment**; the offsets cancel. The real oracle is
region A against region-A-with-a-rigid-shift-δ applied to its input
parameters. Requirements: the alignment solver recovers a step of exactly −δ,
and δQ → 0 on every site. This is an exact root-find test on the genuinely new
machinery and it must be written **before** the interface geometry builder,
on a 1-band or single-l toy if that is faster.

### 5.4 Plain A | A identity
Same material both sides, no offset: every layer comes out bulk, δQ → 0,
dipole → 0, V_B → 0. Needs no new geometry, so it unblocks B7.3/B7.4 before
the cluster builder exists.

### 5.5 Vacuum parameter oracle
Generated empty-lattice parameters vs the analytic spherical-Bessel
free-electron result, over the energy window where ASA is expected to hold.

### 5.6 Vacuum lead vs empty-sphere surface
Layer moments and work function, current `buildsurf` vs new vacuum lead.
**Acceptance criterion: reduced sensitivity to the empty-sphere / vacuum layer
count**, not a changed answer (§1.6).

### 5.7 Compensation insensitivity
Vary the compensation weight profile within the screening region and vary the
active-zone thickness. The reported step must become insensitive to the former
as the latter grows. Feeds **G-B7-3**.

### 5.8 Physics benchmarks
Fe/Cr and Co/Cu interlayer exchange vs spacer thickness, with the
LL-convergence study on the oscillation period. Collinear FM/AFM lead
configurations are a spin-index swap on the frozen set, so magnetic leads come
nearly free.

---

## 6. Hooks: designed for, not built here

- **Applied bias / electric field.** Because alignment enters as a *target
  step*, a bias is `target_step = contact_potential + eV` — one namelist real,
  no new physics path. `PM` (§2.10) is the matching kernel, already built.
  Provide the namelist entry and the plumbing; leave it inactive. Caveat for
  whoever activates it: a biased or field-exposed slab is a charged-slab
  problem and needs its own compensation story.
- **Noncollinear leads.** Collinear FM/AFM is a spin-index swap and is in
  scope. Arbitrary relative lead magnetization needs spin-rotated frozen
  parameters — the one-sided spinor transform from the GBT work
  (H̃_ij = H⁰_ij · U(α_ij)). Flag the hook; keep noncollinear leads out.
- **Relaxed-z 2D Madelung.** Plumbing only (§2.10). Carry z as registry data
  so this is later a parameter change.
- **Transport.** Real-space Kubo–Greenwood/Chebyshev, separate milestone.
- **Per-energy recursion driver.** B9/B12 (CPA, DMFT). Not here.

---

## 7. Checklist

- [x] **B7.0** audit + `exh` fix + B6 re-validation + `vmad` persistence
      answered + `imppot` mixing fixed + factor of two pinned + ~~`SWS`
      extracted~~ + on-site term restored → **G-B7-1**, **G-B7-4** *(both
      signed 2026-07-25. `SWS` extraction **withdrawn** per the sign-off note;
      four of the six suspected bugs were not bugs — see the summary table in
      `CONVENTIONS_MADELUNG.md`)*
- [x] **B7.1** region registry replacing `init` arithmetic; `buildsurf`
      reproduced as a registry instance
- [x] **B7.2** vacuum parameter generator + Bessel oracle *(the oracle needed
      the scalar-relativistic κ to agree beyond ~1e-4 — `potpar` drives
      `RSEQSR`, not a plain free-electron solver)*
- [x] **B7.3** two-sided deviation electrostatics + localized compensation +
      Q/P reporting + overlap diagnostic *(no separate `interfacemat`: it
      reuses `surfmat`'s matrices, maintainer-confirmed)*
- [x] **B7.4** alignment solver + offset-δ oracle → **G-B7-2** *(solver landed;
      `UnitAlignmentSolver` drives the production fixed point to exactly −δ.
      Gate contract drafted in `CONTRACT_FROZEN_REGION.md`, **awaiting
      signature** — three open decisions flagged there)*
- [ ] **B7.5** `buildinterface` path + namelist (`buildsurf` untouched)
- [ ] **B7.6** examples + honest documentation
- [ ] **B7.7** vacuum-lead vs empty-sphere comparison + compensation
      insensitivity → **G-B7-3** + IEC benchmarks with LL convergence
