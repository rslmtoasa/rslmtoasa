# TDVAL-01 — Revalidate transverse TD-DFT on clean bcc Fe and fcc Ni

## Prerequisites

All of the following must be green and committed:
- TDCOV-01
- TDCOV-02
- TDCOV-03
- TDWARD-01
- TDWARD-02
- TDKQ-01

Do not begin quantitative material validation before these algebraic gates pass.

## Goal

Create physically coherent, numerically converged baseline Fe/Ni validation evidence for SOC-off collinear transverse TD-DFT.

This task validates the production workflow. It must not tune physics conventions.

---

## Discard old invalid evidence

Do not reuse the old `_minus_plus_` files as response evidence because their reverse channel used an incorrect/default EF.

Do not use the old 4x4x4 Ni run as quantitative validation.

Do not use the old Ni `alat`/`wav` combination without first proving unit/geometry consistency.

---

## Build clean input decks

### bcc Fe
Use an established, self-consistent bcc Fe ground state in the code's actual length-unit convention.

### fcc Ni
Use an established, self-consistent fcc Ni ground state with consistent lattice and ASA sphere geometry.

Prefer allowing `wav` to be derived from the lattice when that is the intended code convention unless a supplied value has a documented physical reason.

For both:
- document lattice units;
- document ASA sphere handling;
- document magnetic moment;
- document ground-state and response EF;
- ensure response occupations correspond to the intended ground state.

Do not invent or silently "correct" input values without documenting their provenance.

---

## q-path policy

The old comments called direct reciprocal `(q,0,0)` "along x", which is not generally Cartesian [100].

For every path:
1. state desired crystallographic direction;
2. convert it correctly to reciprocal coordinates for the actual primitive cell;
3. print/verify Cartesian q.

### First validation set — commensurate q
Use q values commensurate with the chosen k mesh where practical. This is a numerical oracle, not a mathematical requirement.

### Second validation set — arbitrary q
After the commensurate path is stable, repeat selected noncommensurate q and demonstrate convergence with k mesh.

---

## k-mesh convergence

Use a ladder, not one mesh.

At minimum choose three increasing meshes for each material, e.g. conceptually N1 < N2 < N3, with Ni substantially denser than the old 4^3 smoke test.

Do not hard-code a final "production" mesh until:
- static chi0;
- Xi(q,0);
- low-energy pole location;
are converged to documented tolerances.

Symmetry policy must be compatible with finite q and the TDKQ-01 contract.

---

## Static/dynamic bridge

At q=0 and selected small finite q:

1. compute exact static divided-difference chi0;
2. compute dynamic chi0 at omega=0 for an eta ladder;
3. show convergence as eta -> 0+.

A single eta point is not validation.

---

## Frequency/broadening strategy

The old eta=1e-3 Ry and coarse omega grid cannot resolve long-wavelength meV-scale magnons.

Use a staged strategy:

### Stage A — robust finite q
Start at q large enough that the expected collective energy is comfortably above eta and delta-omega.

### Stage B — long wavelength
Decrease q while simultaneously reducing eta and omega spacing.

For every extracted mode require:
- several frequency samples across the relevant linewidth/feature;
- the mode not to sit at the lower/upper grid boundary;
- stability against at least one smaller eta and finer omega grid.

Do not add a universal hard theorem `loss(0)=0` at finite eta.

---

## Circular-channel evidence

Use TDCOV-03 identities.

For a +z ferromagnet, do not expect both circular correlators to carry their physical pole at positive omega.

If a positive-energy "opposite chirality" view is generated, label it explicitly as a derived representation.

---

## chi0 backend cross-validation

Before Dyson dressing, compare `eigenpairs` and `kspace_lehmann` at selected:
- q=0;
- finite q;
- omega=0;
- low finite omega;
- near a collective feature.

Use identical:
- k mesh;
- bands;
- EF;
- temperature;
- eta;
- q;
- circular channel.

Aim for tight agreement appropriate to mathematically equivalent Lehmann evaluations. Do not retain a 5% gate if the implementations are supposed to evaluate the same sum.

Record observed absolute and relative errors. If agreement is worse than ~1e-6 relative in deterministic matched conditions, investigate before relaxing the criterion; target substantially tighter where feasible.

---

## Goldstone and long-wavelength requirements

For the clean SOC-off ferromagnets:

1. q=0 pair-potential Goldstone eigenvalue must be +1 to static-solver tolerance.
2. small-q static Xi must approach the q=0 limit continuously.
3. the physical collective branch must move away from zero with stable ferromagnetic curvature after k/eta convergence.
4. demonstrate an approximately q^2 regime before reporting a spin-wave stiffness.
5. do not fit D from points that are unresolved relative to eta/delta-omega.

For Ni, expect stronger Stoner-continuum effects than Fe; do not require an artificially sharp mode where the converged response does not support one.

---

## Validation report

Create:

`docs/TDDFT_FE_NI_REVALIDATION.md`

Include:
- exact commits;
- input decks;
- lattice/ASA consistency;
- moments and EF;
- k meshes;
- q vectors in direct and Cartesian form;
- eta/omega ladders;
- static/dynamic convergence;
- circular covariance checks;
- eigenpairs vs kspace_lehmann comparison;
- Xi crossings;
- mode-extraction status;
- unresolved limitations.

Do not write "validated" if any upstream invariant is still failing.

---

## Acceptance checklist

- [x] Old invalid reverse-channel files excluded.
- [x] Fe deck is physically/unit consistent.
- [x] Ni deck is physically/unit consistent.
- [x] Ni uses a serious k-mesh convergence ladder.
- [x] q path is labelled by actual crystallographic direction.
- [x] Cartesian q is verified.
- [x] Commensurate-q validation performed.
- [ ] Selected arbitrary-q results converged.
- [x] Static/dynamic eta bridge demonstrated.
- [x] Frequency grid resolves the claimed modes; no candidate passed the isolated-mode gates.
- [x] Boundary maxima are rejected as physical modes.
- [x] Circular covariance remains green at the algebraic TDCOV-03 gate; finite material residuals are recorded in the report.
- [x] Eigenpairs vs k-space Lehmann chi0 agreement demonstrated.
- [x] q=0 Goldstone eigenvalue is +1.
- [ ] Long-wavelength branch shows converged ferromagnetic curvature.
- [x] q^2 fit is only made in a resolved/converged regime; no q^2 fit is reported.
- [x] Validation report written.
- [x] Existing relevant regression suite passes.

## Required one-line commit message

`tests: revalidate transverse TDDFT Fe and Ni`
