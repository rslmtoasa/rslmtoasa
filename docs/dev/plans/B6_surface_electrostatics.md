# B6 — Dipole charge moments: surface/interface electrostatics

**Status (2026-07-28):** implemented, at **narrower scope** than the
original blueprint below. Reconciled during the docs/dev/plans
homogenization pass from two previously separate, un-tracked documents:
this file's §1 "Original blueprint spec" (ex-`claude_files/plans/
B6_surface_electrostatics.md`, matching `BLUEPRINTS.md` §B6 — full 3D/2D
Ewald generalization to l ≤ 2) and §2 "As-shipped design" (ex-`docs/dev/
electrostatic_dipole_plan.md` — the actual, narrower implementation). No
document previously recorded that these two describe different things; a
reader who only saw §1 (the blueprint) would not recognize the shipped
feature.

**What actually shipped, in one paragraph:** rather than a new generalized
Ewald/Madelung lattice-sum machinery (§1's B6.2/B6.3), the implementation
reuses the **existing** `charge%madl2d` dipole-monopole (`dsz`) /
dipole-dipole (`dzz`) Madelung matrices, which were already present in the
code but previously unfed. The new work is entirely on the charge-moment
side: `source/electrostatics_multipole.f90` (`compute_dipole_moments`)
computes the l=1 (z) dipole moment Q₁₀ per atom from on-site cross-orbital
density-matrix elements (s–p_z, p_z–d_z²) and a newly exported radial
partial-wave amplitude (`potential%q10`, fed via a `phi_amp` export hook
from `NEWRHO`), with SCF mixing. This is a materially smaller and lower-risk
piece of work than §1 describes — no new lattice-sum code, no 2D Ewald, no
l=2 moments, no work-function validation campaign.

**Verified:** bulk limit / feature-disabled path reduces bit-identically to
the monopole result (Q10 → 0) — confirmed against source in the Phase-3
sweep (`docs/dev/PHASE3_STATUS.md` §1, B6 entry).

**Not done:** gate G-B6-1 (which Skriver–Rosengaard reference tables define
acceptance) is unsigned/not reached — no literature work-function
comparison was performed. **Zero test/example/doc coverage** — no
`cases.json` entry, no `example/` directory, no rst page. Flagged in
`docs/dev/PHASE3_STATUS.md` §5 as needing a highlight example built from
scratch.

---

## 1. Original blueprint spec (superseded in scope — kept for context)

**Effort:** M. **Lead:** OPUS (Ewald/2D lattice sums — the risk item),
SONNET (moments integration + plumbing). **Depends on:** none — runs in
parallel with B2–B5. **Prerequisite for B7 being physically trustworthy**
(semi-infinite leads with wrong interface electrostatics give beautiful
garbage).

### 1.1 Physics spec

Extend the ASA Madelung problem from monopoles q_R to include l = 1 (and
optionally l = 2, behind a namelist switch defaulting OFF) multipole
moments of the sphere charge densities:

```
Q_RL = ∫_sphere d³r  n_R(r) · r^l · Y_L(r̂)        (real harmonics, l ≤ 1(2))
V_RL = Σ_{R'L'}  M_{RL,R'L'} · Q_{R'L'}
```

- **Moments (l=1):** radial integrals over the same meshes the charge
  machinery owns; new code consumes `charge.f90` data but lives in a new
  module `source/electrostatics_multipole.f90` — outside the legacy fence.
- **Generalized Madelung matrix M:** structure-constant-like lattice sums.
  Bulk 3D: standard Ewald generalization to l ≤ 2. Surface/slab: 2D Ewald
  (in-plane sums + out-of-plane real-space decay), reusing/extending the 2D
  conventions from the ASA preprocessing toolkit (layer identification,
  vacuum margin; **L_z remains derived, never input** — the standing 2D
  convention).
- **SCF feed-in — the fence pattern (binding):** the new module computes
  potential shifts v_R (and l=1 components where the representation
  allows); legacy code receives them **through the existing variable it
  already reads** (the monopole Madelung shift), extended *additively*.
  Zero structural edits inside `self.f90`. If an l=1 potential component
  cannot be represented through existing variables, it is dropped from the
  potential (kept in the energy) and this limitation is documented — do
  not invent a new channel into the legacy solver.

Convergence/efficiency: Ewald splitting parameter chosen automatically
(standard error-balance heuristic, documented); M built once per geometry
and cached — it is charge-independent; per-SCF-iteration cost is one
matrix-vector product, negligible.

### 1.2 Validation (known-answer) — as originally specified

1. **Bulk limit:** all dipole moments → 0 by symmetry; total energy and
   potential shifts **bit-identical** to current code with the feature on
   (regression contract — this is the primary safety net for the fence
   pattern).
2. **Madelung matrix unit tests:** (a) l=0 block reproduces the existing
   monopole Madelung constants (NaCl-structure constant to 1e-10);
   (b) 2D Ewald vs brute-force real-space sum on a small slab to 1e-8;
   (c) Ewald-parameter independence (two splittings, same M to 1e-10).
3. **fcc(100)/(111) work functions** vs published LMTO-ASA-dipole values.
   **Gate G-B6-1:** Anders names the exact Skriver–Rosengaard table
   entries and acceptable deviation before the comparison is scored.
4. **Slab thickness convergence** of the surface dipole (monotone approach,
   documented curve).

### 1.3 Tasks (original scope)

- **B6.1 [SONNET]** `electrostatics_multipole` module skeleton + Q_RL
  moment integrals + unit tests against analytic densities (Gaussian test
  charge: closed-form moments). *Kit:* this file; `charge.f90` radial-mesh
  accessors; the module-pattern of any small existing class.
- **B6.2 [OPUS]** Bulk 3D generalized-Ewald M to l ≤ 2 + unit tests 2(a),
  2(c). *Kit:* this file §1–2; B6.1 interfaces. No repo browsing needed —
  this is self-contained math.
- **B6.3 [OPUS]** 2D/slab Ewald + test 2(b); adopt the preprocessing
  toolkit's layer/vacuum conventions (paste its convention section into
  the session). *Kit:* B6.2; the 2D-convention excerpt.
- **B6.4 [SONNET]** SCF feed-in via the fence pattern + namelist switch
  (default OFF) + bulk bit-identity regression (test 1). *Kit:* B6.1–B6.3;
  the exact code site where the monopole Madelung shift enters the
  potential update (locate via DEVELOPER_MAP; paste ±30 lines only).
- **B6.5 [SONNET]** Work-function + slab-convergence validation cases
  (tests 3–4, after G-B6-1). *Kit:* B6.4; an existing surface example
  (fccCu001) as template.

### 1.4 Checklist (original scope — superseded, see §2 for what shipped)

- [ ] B6.1 moments — **not built as specified; superseded, see §2**
- [ ] B6.2 3D Ewald M — **not built; scope dropped**
- [ ] B6.3 2D Ewald M — **not built; scope dropped**
- [ ] B6.4 SCF feed-in + bulk bit-identity — **done, but via the §2 design**
- [ ] B6.5 work functions (G-B6-1 signed) — **not done, gate unsigned**

---

## 2. As-shipped design (what actually landed)

### 2.1 Revised objective

Extend the existing LMTO-ASA monopole electrostatics to include l=1 (dipole)
multipole moments for broken-symmetry systems, using a unified full
density-matrix approach — **reusing the Madelung matrices already present
in the code** (`dsz`, `dzz`) rather than building new lattice-sum machinery.

### 2.2 Phase 1: Full on-site density matrix extraction

**Target:** the `potential` type definition and a new density-matrix
extraction routine, called from the SCF loop alongside `calculate_moments`.
**Physics:** the complete on-site density matrix D^σ_{L1,L2} contains all
spatial charge-distribution information for a sphere, including
cross-orbital hybridizations — it is the energy integral of Im[G] up to the
Fermi level.

- Allocate a full orbital-resolved density-matrix array for both spins.
- Loop over atoms, spins, and orbital-index pairs; extract
  y(E) = −(1/π) Im[G₀(iorb, jorb, E, na)] (same `spin_off` convention as
  `calculate_moments`).
- Integrate with the existing `simpson_m` routine up to `en%fermi`.
- Store into the atom's full density-matrix array.

### 2.3 Phase 2: Radial dipole density construction

**Target:** `self.f90` (`NEWRHO`). **Physics:** the l=1 radial density
n₁(r) is built from cross-products of the normalized radial partial waves
φ_l(r), weighted by the off-diagonal density-matrix elements and the
corresponding Gaunt coefficients C_{L1L2,10}.

- In `NEWRHO`, alongside the existing spherical density array, add a
  parallel dipole-density array.
- Extract D_{sp_z} and D_{p_z d_z²} from the density matrix built in Phase 1.
- For each radial mesh point: n₁(r) = 2·C_sp·D_{sp_z}·G_s(r)·G_p(r) +
  2·C_pd·D_{p_z d_z²}·G_p(r)·G_d(r).

### 2.4 Phase 3: Multipole moment integration

**Target:** `self.f90` (`NEWRHO` or a dedicated atomic module — landed as
`source/electrostatics_multipole.f90`, `compute_dipole_moments`).
**Physics:** the dipole moment is the volume integral of the l=1 charge
density over the WS sphere; the volume element gives r², the dipole
operator gives another r, hence an r³ weighting.

- Integrate over the logarithmic radial mesh with the existing Simpson's-
  rule weighting and `drdi` volume-element convention.
- Q₁₀ = Σ_ir  wgt·drdi·n₁(ir)·r(ir)³.
- Store Q₁₀ in `symbolic_atom%potential%q10` (this field name matches the
  shipped code — see `source/electrostatics_multipole.f90` lines ~140).

### 2.5 Phase 4: Electrostatic potential shifts

**Target:** `charge.f90` (`surfpot`/Madelung feed-in). **Physics:** the
potential shift for a sphere must include dipole contributions from all
other spheres, coupled via the Madelung matrices that already existed in
the code (previously computed but unfed).

- In the existing monopole-potential loop (using `dss(iq,jq)`), add the
  dipole contribution using the pre-existing dipole-monopole matrix
  `dsz(iq,jq)`:
  ΔV_i = Σ_j [ M_ss(i,j)·Q₀(j) + M_sz(i,j)·Q₁₀(j) ].
- The l=1 potential-shift components (coupled via `dzz(iq,jq)`) are stored
  if SCF requires the dipole potential gradient, but updating the monopole
  shift ΔV_i alone is sufficient for work-function/core-level corrections
  and is what was implemented (commit `8c28ef6`, "activate l=1 dipole
  surface electrostatics via Madelung dsz").

### 2.6 As-shipped checklist

- [x] Full density-matrix extraction (Phase 1) — `potential%q10` storage +
      NEWRHO `phi_amp` export hook (commit `1d81603`).
- [x] `electrostatics_multipole.f90` module + `compute_dipole_moments`
      (commit `af76b7e`, initially unused/inert; wired live in `8c28ef6`).
- [x] SCF feed-in via existing `dsz` Madelung matrix, gated by
      `control%dipole_electrostatics` (commit `8c28ef6`).
- [x] Bulk-limit bit-identity verified (Q10 → 0 when disabled/symmetric).
- [ ] Work-function validation vs Skriver–Rosengaard tables — **gate
      G-B6-1 not reached, no reference tables named yet.**
- [ ] Highlight example/test/doc — **none exist; see
      `docs/dev/PHASE3_STATUS.md` §5.**
