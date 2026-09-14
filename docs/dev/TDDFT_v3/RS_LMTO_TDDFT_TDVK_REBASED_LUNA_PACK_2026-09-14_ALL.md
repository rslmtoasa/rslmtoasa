<!-- FILE: README.md -->

# RS-LMTO-ASA reciprocal TD-DFT — rebased Luna prompt pack

Date: 2026-09-14  
Repository: `https://github.com/rslmtoasa/rslmtoasa`  
Branch: `fable_v4`

## Purpose

This pack restores the parent TDVK milestones as the controlling roadmap after
the product-basis remediation and TDVK-03A quadrature audit.

The pack is deliberately **thin but not incremental**:

- one parent prompt per scientific milestone;
- implementation details remain internal checkpoints;
- a new named remediation task is allowed only for a concrete blocker;
- every blocker task must return immediately to the parent milestone after PASS;
- no automatic chain such as `03B -> 03C -> 03D -> ...`.

The latest handoff-confirmed audit commit is:

`4da3c115330cd2c7c0969c7d4418ec5f1d6a636c`

Do **not** assume that this remains the live branch HEAD. Every task starts by
inspecting the live branch and current evidence.

## Current scientific position

Accepted infrastructure includes:

- clean-room reciprocal transverse TD-DFT target:
  collinear, no SOC, Pauli transverse response from the scalar-relativistic
  LMTO ground state;
- exact finite LMTO product-response representation;
- complete one-site `spd`, `Lmax_response=4` product dimension = 232;
- compact Lehmann `chiKS`;
- independent compact reciprocal-GF real-axis Kubo `chiKS`;
- compact/legacy representation cross-checks at approximately machine precision
  on affordable fixtures;
- real bcc-Fe Gamma compact Lehmann material execution;
- TDVK-03A result:
  `REAL-AXIS GF FORMULATION NUMERICALLY CONSISTENT`;
- TDVK-03A also established that the current scalar GF bubble is too slow for
  Fe material closure without mathematically identical performance remediation.

The point-grid response remains a certification/reference representation.
Do not rebuild 12,375-dimensional Fe point-response matrices in production.

## Execution order

1. `01_TDVK-03_FE_BACKEND_CLOSURE.md`
2. `02_TDVK-04_FE_FINITE_Q_COVARIANCE.md`
3. `03_TDVK-05_FE_CONVERGENCE.md`
4. **STOP and return to the orchestrator using**
   `04_REVIEW_GATE_AFTER_TDVK-05.md`
5. Only after explicit approval:
   `05_TDVK-06_STATIC_INTERACTIONS.md`
6. `06_TDVK-07_FE_DYSON_LOSS.md`
7. `07_TDVK-08_NI_REPLICATION.md`
8. `08_TDVK-09_EVIDENCE_HANDOFF.md`

Read `09_CAMPAIGN_RULES.md` before starting.

## Repository evidence to read first

At minimum inspect the live versions of:

- `docs/LR_TDDFT_CONVENTIONS.md`
- `docs/LR_RESPONSE_SPACE_ALGEBRA.md`
- `docs/LR_LMTO_PRODUCT_RESPONSE_BASIS.md`
- `docs/LR_KS_SUSCEPTIBILITY.md`
- `docs/LR_GF_SUSCEPTIBILITY_CROSSCHECK.md`
- `docs/TDDFT_RECIPROCAL_GF_QUADRATURE_AUDIT.md`
- `docs/TDDFT_RECIPROCAL_VALIDATION.md`
- `docs/KXC_ALSDA_TRANSVERSE_KERNEL.md`
- `docs/GOLDSTONE_SUMRULE_INTERACTION.md`
- `docs/TDDFT_DYSON_AND_LOSS.md`
- `docs/TDDFT_PRODUCTION_DRIVER.md`

Live source and live evidence override stale text in this pack.

## Luna role

Luna is an implementation and evidence worker.

Luna does **not** decide:

- new physics;
- signs, factors, circular conventions or Fourier phases;
- whether a material discrepancy is physically acceptable;
- whether a loss peak is a magnon;
- whether a Goldstone residual should be repaired;
- whether literature agreement has been achieved.

Ambiguous physics evidence is returned to the orchestrator unchanged.


<!-- FILE: 09_CAMPAIGN_RULES.md -->

# Reciprocal TD-DFT campaign rules

## Roles

### Orchestrator owns

- physics;
- equations;
- signs/factors;
- conventions;
- acceptance criteria;
- methodological changes;
- interpretation;
- material conclusions;
- task decomposition when a blocker appears.

### Luna may

- implement prescribed formulas;
- make narrow mechanical refactors;
- add tests using specified independent oracles;
- run numerical ladders;
- gather evidence;
- document exact failures.

### Luna may not

- derive replacement physics;
- tune signs/factors;
- change circular conventions;
- change q phases;
- pick a different channel because it looks better;
- rescale XC kernels;
- weaken tests;
- replace the prescribed numerical method because it is slow;
- interpret spectra as modes;
- claim literature agreement.

## Parent-milestone rule

Do not create a new named task merely because implementation has several steps.

A new remediation slice is allowed only when:

1. the parent milestone cannot complete;
2. the blocker is concrete;
3. ownership is identifiable;
4. the remediation has one narrow purpose;
5. it has a binary exit;
6. after PASS the next action is explicitly “return to parent milestone.”

## Evidence classes

Always distinguish:

1. algebraic consistency;
2. independent numerical cross-check;
3. numerical convergence;
4. material execution;
5. converged-material validation;
6. literature agreement.

Never promote one class into another.

## Independent-oracle rule

Strong examples:

- point transition vs analytical compact transition;
- projected point GF vs compact GF;
- real-axis spectral moment integral vs direct eigensystem sum;
- scalar GF contraction vs optimized contraction;
- Lehmann vs independent real-axis GF;
- compact operator/Dyson result vs independently projected affordable
  point-space result.

Weak examples:

- expected and actual values generated by the same helper;
- reproducing a made-up constant;
- calling a synthetic self-consistency test material validation.

## Material-state rule

Backend comparisons must use the same accepted state whenever the comparison is
claimed to isolate backend differences.

Record:

- eigenvalues/eigenvectors;
- occupations;
- EF;
- moment;
- k mesh;
- temperature;
- q endpoint;
- radial basis;
- compact product basis;
- channel;
- physical eta.

## Bare-before-interacting rule

Mandatory order:

accepted state  
→ compact Lehmann `chiKS`  
→ independent compact reciprocal-GF `chiKS`  
→ backend closure and convergence  
→ static KXC/GSR  
→ Dyson/loss  
→ mode interpretation

Do not use interaction physics to repair an unvalidated bare response.

## Product-space rule

For accepted Fe `spd`, the complete one-site product span is 232 coordinates
under the current certified `Lmax_response=4` construction.

Do not:

- reintroduce radial decimation;
- use site-only response;
- truncate accepted full-rank product modes for speed;
- insert extra radial weights;
- use inverse singular values in runtime transition evaluation;
- reconstruct the large point-grid response in production;
- merge Lehmann and reciprocal GF into one ambiguous backend.

## Performance rule

Optimization is allowed only with the reference kernel retained.

For each optimized kernel:

1. identical inputs;
2. complete-output comparison;
3. roundoff/regression-level equality;
4. benchmark only after correctness passes.

Do not change physics and optimization simultaneously.

## Interaction-route rule

Keep provenance separate for:

- direct ALSDA;
- independent GSR;
- optional BES/GCR.

Never implement GSR as rescaled ALSDA.

Do not use BES/GCR to conceal a broken bare/kernel layer.

## Fail-closed rule

If a task fails:

- preserve the failing input/test;
- record exact evidence;
- do not patch the symptom;
- identify the owning layer if mechanically possible;
- return to the orchestrator.

No hidden fallback.
No silent reduced response space.
No sign tuning.
No ad hoc Goldstone shift.

## Legacy-code rule

Protect legacy atomic LMTO routines unless evidence identifies them as the defect
owner.

Prefer narrow adapters and certified response modules over invasive changes to
legacy radial/atomic machinery.


<!-- FILE: 01_TDVK-03_FE_BACKEND_CLOSURE.md -->

# TDVK-03 — bcc Fe reciprocal-backend closure

## Parent objective

On **one accepted bcc-Fe reciprocal electronic state**, determine whether the
two independent bare-response backends

- compact Lehmann `chiKS`;
- compact reciprocal-GF real-axis Kubo `chiKS`

approach the same Gamma response under controlled GF integration.

This is one parent task. The checkpoints below are internal to TDVK-03.
Do not split them into separately named campaign stages unless one produces a
hard blocker.

## 0. Live-state preflight

Before changing code:

1. fetch and checkout `fable_v4`;
2. record exact HEAD and `git status`;
3. identify commits after the handoff-confirmed
   `4da3c115330cd2c7c0969c7d4418ec5f1d6a636c`;
4. read the live response/product/GF/quadrature/validation documents;
5. determine whether any checkpoint below has already landed;
6. preserve newer evidence and do not redo completed work without cause.

If the live source contradicts this prompt, stop and report the discrepancy.

## 1. Mixed-complex-eigensystem GF regression

The TDVK-03A quadrature fixture used identity eigenvectors. Add one strong,
deterministic regression using a genuinely mixed complex Hermitian eigensystem
obtained by numerical diagonalization.

Requirements:

- nontrivial off-diagonal resolvents;
- complex eigenvectors;
- no hand-written expected susceptibility matrix;
- exercise the existing compact product basis and component vertices;
- preserve the independent GF/Kubo route;
- compare against independent Lehmann/spectral information already available
  from the eigensystem;
- reuse existing repository tolerances for comparable double-precision response
  regressions; do not weaken them.

Purpose:

- catch transpose/conjugation/indexing errors;
- show that TDVK-03A behavior is not special to identity eigenvectors.

If this gate fails, preserve the failing fixture and stop TDVK-03.

## 2. GF bubble performance remediation

The present scalar compact-GF bubble is the correctness oracle but is too slow
for Fe.

Refactor **only the contraction implementation**.

Must preserve exactly:

- real-axis GF/Kubo formulation;
- current real-energy quadrature semantics;
- weighted resolvents;
- both Kubo terms;
- component vertices;
- signs;
- conjugations/transposes;
- prefactors;
- channel convention;
- physical response `eta`;
- GF integration broadening semantics.

Do not:

- call the Lehmann accumulator under the hood;
- change quadrature physics to make the calculation faster;
- truncate the 232-dimensional Fe product space;
- reconstruct point-space response matrices;
- change LR-03/LR-05 conventions.

Preferred engineering direction:

- remove repeated scalar/small-matrix work;
- express mathematically identical contractions with BLAS/matrix operations where
  appropriate;
- reuse intermediate quantities that are invariant over the response indices;
- avoid unnecessary allocation/copying in the energy loop.

Keep the scalar implementation available as an oracle.

### Optimization acceptance gate

On affordable fixtures, run scalar and optimized implementations on **identical
inputs** and compare the complete compact matrices.

Require equality at the existing response-regression precision/roundoff level.
Record at least:

- `||chi_scalar||_F`;
- `||chi_opt||_F`;
- `dF = ||chi_opt-chi_scalar||_F`;
- relative Frobenius difference;
- max-element absolute difference;
- wall time;
- speedup.

Correctness is decided before benchmarking.

If the optimized implementation cannot be shown numerically identical, keep the
scalar implementation as production owner and report BLOCKED.

## 3. Reuse one accepted Fe reciprocal state

Prepare one accepted bcc-Fe reciprocal state and reuse it across the backend and
GF-integration comparisons.

Do not reconverge SCF for every GF integration sample.

Reuse the simplest existing mechanism. Do not introduce a general persistence
framework unless the live architecture already provides one.

Record provenance sufficient to prove that every backend comparison uses the
same:

- Hamiltonian/electronic state;
- eigenvalues;
- eigenvectors;
- occupations;
- Fermi level;
- moment;
- temperature;
- k mesh;
- q endpoint;
- radial basis;
- compact product basis;
- channel;
- physical response `eta`.

If identical-state reuse cannot be demonstrated, do not call the observed
difference a backend error.

## 4. Fe Gamma backend-closure ladder

Start with:

- bcc Fe;
- Gamma;
- `chi_plus`;
- `omega = 0`;
- the established physical response `eta`;
- complete product response space;
- expected one-site `spd` dimension = 232, subject to live provenance.

Compute the compact Lehmann result once from the accepted state.

Then run a controlled compact reciprocal-GF ladder.

### Control the three GF numerical effects separately

#### A. Simpson resolution

For each chosen `integration_eta`, choose odd Simpson grids that explicitly
resolve the spectral width.

Use the TDVK-03A evidence as the design rule:

`h / integration_eta ≲ 0.5`

for the controlled side of the ladder.

Include enough information to show whether remaining change is quadrature error.
Do not use the old hard-coded 1001/2001/4001 ladder blindly.

#### B. finite GF `integration_eta`

At quadrature settings already demonstrated to resolve each width, use a small
systematic ladder of GF integration widths based on the existing TDVK-03A
settings/live code capabilities.

The purpose is to determine whether compact GF moves systematically toward
compact Lehmann as `integration_eta` decreases.

Do not confuse GF `integration_eta` with the physical response `eta`.

Do not reduce `integration_eta` beyond a practically supportable range merely
to force agreement.

#### C. energy window

At one otherwise controlled setting, vary the integration margin/window enough
to establish whether window truncation is subdominant.

Do not tune the window to improve agreement.

## 5. Matrix evidence

For every Fe comparison record the complete compact-matrix diagnostics:

- Lehmann Frobenius norm;
- GF Frobenius norm;
- `dF`;
- relative Frobenius difference;
- `dInf` / max-element absolute difference;
- physical response `eta`;
- GF `integration_eta`;
- energy bounds/margin;
- number of integration points;
- energy spacing `h`;
- `h/integration_eta`;
- wall time;
- exact accepted-state provenance.

Do not replace the full-matrix comparison with only a trace, site scalar or
selected element.

One finite diagnostic frequency may be added **only after** Gamma static closure
is numerically controlled. Do not create a spectrum.

## 6. TDVK-03 verdict

TDVK-03 can return `PASS CANDIDATE` only if:

1. the mixed-complex regression passes;
2. optimized GF matches the scalar oracle at roundoff/regression precision;
3. the optimized GF route remains an independent real-axis Kubo/GF route;
4. one identical accepted Fe state is used for both backends;
5. Simpson resolution is demonstrably controlled;
6. finite `integration_eta` behavior is systematic;
7. energy-window sensitivity is characterized;
8. compact GF approaches compact Lehmann with no evidence for a stable
   algebraic offset;
9. runtime is practical enough for the representative Fe checks required later;
10. complete provenance is recorded.

Luna does **not** invent a material PASS tolerance for the final GF↔Lehmann
difference. Return the numerical evidence to the orchestrator.

If the GF route remains prohibitively expensive after a proven mathematically
identical optimization, report:

`BLOCKED — RECIPROCAL GF MATERIAL COST`

Do not weaken the material comparison.

## Documentation

Update:

`docs/TDDFT_RECIPROCAL_VALIDATION.md`

with one section:

`TDVK-03 Fe reciprocal-backend closure`

Keep implementation-oracle evidence and Fe material evidence clearly separate.

## Return to orchestrator

Return:

- starting and final HEAD;
- changed files;
- tests run;
- mixed-eigensystem result;
- scalar-vs-optimized matrix comparison and speedup;
- accepted Fe-state provenance;
- compact Fe GF integration table;
- final `PASS CANDIDATE` or `BLOCKED`;
- no physical interpretation.

Then stop.

## Commit

Suggested one-line commit:

`td-dft: close Fe reciprocal bare-response backends`


<!-- FILE: 02_TDVK-04_FE_FINITE_Q_COVARIANCE.md -->

# TDVK-04 — Fe finite-q endpoint, folding and covariance

## Prerequisite

Run only after orchestrator approval of TDVK-03.

## Objective

Verify the accepted Fe reciprocal bare-response path at finite q before any
dispersion or magnon interpretation.

Use the compact orthonormal product representation directly.
Do not reconstruct the legacy point-response matrix in production.

## 0. Preflight

Record live HEAD/status and read the current reciprocal validation document.
Reuse the TDVK-03 accepted Fe state wherever the live architecture permits.

No response-physics changes are authorized.

## q set

Use direct reciprocal coordinates and report them literally:

1. Gamma;
2. one small mesh-commensurate `q`;
3. exact `-q`;
4. one arbitrary/off-mesh `q`.

Do not attach symmetry-line labels unless the live reciprocal path establishes
them unambiguously.

Use the same:

- accepted state;
- compact response basis;
- channel;
- physical `eta`;
- selected low/static `omega`.

## Backend order

1. Lehmann at all required q points.
2. Reciprocal GF only at representative finite-q points after Lehmann succeeds.

Do not turn TDVK-04 into a second GF convergence campaign.
Use the TDVK-03-approved GF numerical controls.

## Validate

- exact `k+q` endpoint handling;
- reciprocal-zone folding;
- exact metadata for supplied and folded q;
- existing q↔-q covariance convention;
- compact full-matrix covariance diagnostics;
- representative Lehmann↔GF finite-q consistency.

Use the already established LR-03/LR-06 covariance formula/helpers.
Do not derive or insert a new channel/conjugation identity.

If a validation helper is needed, it must reuse the existing convention rather
than create new response physics.

## Forbidden fixes

If q/-q or folding behavior is anomalous, do not:

- add a phase by hand;
- manually conjugate eigenvectors;
- swap circular channels;
- change Fourier signs;
- replace exact `k+q` endpoints by nearest mesh points.

Preserve the failing case and stop.

## Evidence

For every q record:

- supplied q;
- folded endpoint metadata;
- backend;
- compact matrix norm/trace diagnostics;
- covariance residual/diagnostic using the established convention;
- for GF spot checks: `dF`, relative Frobenius difference and max-element
  difference;
- accepted-state identifier/provenance.

## PASS boundary

TDVK-04 may establish:

- arbitrary-q production execution;
- exact endpoint/folding behavior;
- material q/-q covariance;
- representative finite-q backend consistency.

It does **not** establish:

- a magnon branch;
- stiffness;
- damping;
- literature agreement.

## Documentation

Append one `TDVK-04 finite-q/covariance` section to:

`docs/TDDFT_RECIPROCAL_VALIDATION.md`

## Return to orchestrator

Return the q table, covariance diagnostics, any GF spot checks and
`PASS CANDIDATE`/`BLOCKED`.

Then stop.

## TDVK-04 result

**PASS CANDIDATE** — the Fe finite-q compact bare-response validation completed
with exact folded endpoints, material q↔-q covariance, and representative
finite-q Lehmann↔GF diagnostics. No downstream physical interpretation was
performed.

## TDVK-04 completion checklist

- [x] Recorded live HEAD/status and read the live reciprocal validation evidence.
- [x] Reused one accepted Fe reciprocal state, compact product basis, channel,
  physical eta, and static frequency across the finite-q sweep.
- [x] Evaluated literal Gamma, mesh-commensurate q, exact -q, and arbitrary/off-mesh q.
- [x] Ran compact Lehmann at every required q before the representative GF checks.
- [x] Verified exact k+q endpoint metadata and half-open reciprocal-zone folding.
- [x] Applied the established q↔-q retarded/advanced circular covariance convention,
  including the compact-coordinate transport required by the complex product basis.
- [x] Recorded compact full-matrix norms, traces, transition counts, finiteness, and provenance.
- [x] Ran representative finite-q reciprocal-GF checks with the TDVK-03 controls.
- [x] Appended the TDVK-04 evidence section to `docs/TDDFT_RECIPROCAL_VALIDATION.md`.
- [x] Returned a scoped `PASS CANDIDATE`/`BLOCKED` result without magnon, stiffness,
  damping, or literature interpretation.

## Commit

`tests: validate Fe finite-q reciprocal TDDFT`


<!-- FILE: 03_TDVK-05_FE_CONVERGENCE.md -->

# TDVK-05 — Fe bare-response numerical convergence

## Prerequisite

Run only after orchestrator approval of TDVK-04.

## Objective

Establish whether the **bare compact Fe response** is numerically stable enough
to justify moving to static interaction/Goldstone physics.

This is a convergence milestone, not a response-physics modification task.

## 0. Preflight

Record live HEAD/status and the TDVK-03/04 accepted evidence.

Use Lehmann as the primary production convergence backend.
Use reciprocal GF only for a small number of already-certified spot checks.

Do not create a full GF convergence campaign here.

## A. k-mesh ladder

Use at least three systematically increasing k meshes around the established
reference.

If the live input supports the familiar sequence transparently, use a
coarse/reference/finer ladder such as:

- `4^3`;
- `8^3`;
- `12^3`.

If mesh-control syntax differs, use the live established control and record the
actual generated meshes. Do not guess.

For each mesh:

- reconverge the same physical Fe setup;
- record moment, EF and LR acceptance/provenance;
- run compact Lehmann `chiKS` for:
  - Gamma static;
  - one TDVK-04-certified finite q;
  - one low finite omega;
- record full compact-response diagnostics and runtime.

Do not force the EF or moment from one mesh onto another.

## B. physical response-eta ladder

At one orchestrator-approved k mesh, use a small decreasing physical
response-`eta` ladder based on the existing production settings.

The original campaign values `0.02`, `0.01`, `0.005 Ry` remain acceptable if
they are still supported by the live implementation.

Use identical q/frequency points.

This ladder concerns **physical response eta**.
Do not mix it with GF `integration_eta`.

Do not extrapolate eta to zero.

## C. response-space sensitivity

For `spd`:

- complete reference: full `Lmax_response=4` product span
  (`response_lmax=-1` or explicit live equivalent);
- optional reduced diagnostic: `Lmax_response=2`.

Record complete-space product dimension/rank provenance.

Label the reduced case explicitly as approximate.
Do not treat it as a production replacement.

Do not truncate accepted full-rank product modes for speed.

## D. radial provenance

The accepted direct radial mesh remains the physical reference.

Record its identity/provenance.

Do not:

- implement radial decimation;
- introduce a material-dependent radial truncation;
- replace the compact exact LMTO product span with a site-only response.

## GF spot checks

Use the TDVK-03-approved reciprocal-GF controls at only the representative
points needed to ensure the backend closure did not become pathological under
the selected converged-enough settings.

Do not rerun expensive GF ladders on every k mesh or eta.

## Evidence

Create compact tables for:

- k-mesh dependence;
- physical eta dependence;
- full vs reduced response cutoff;
- radial/product-basis provenance;
- representative GF spot checks.

Use the same matrix diagnostics as TDVK-03/04 where applicable.

## Stop condition

If no stable trend is visible with k mesh or physical eta:

- preserve the data;
- do not proceed to KXC/GSR/Dyson;
- return `BLOCKED — BARE RESPONSE NOT NUMERICALLY STABLE`.

Luna must not invent a material convergence threshold.

## Documentation

Append a `TDVK-05 Fe numerical convergence` section to:

`docs/TDDFT_RECIPROCAL_VALIDATION.md`

## Return to orchestrator

Return the tables, exact settings, failures if any, and
`PASS CANDIDATE`/`BLOCKED`.

Then stop at the mandatory TDVK-05 review gate.

## TDVK-05 result

`PASS CANDIDATE` — the compact bare Fe response is finite and numerically
stable over the required `4^3`, `8^3`, and `12^3` meshes, the physical
`eta={0.02,0.01,0.005}` Ry ladder at `8^3`, and the full/reduced response-cutoff
diagnostic. The complete measured evidence is appended to
`docs/TDDFT_RECIPROCAL_VALIDATION.md`. No convergence threshold was invented,
and no KXC, Goldstone, Dyson, loss, or mode-interpretation route was entered.

## TDVK-05 completion checklist

- [x] Recorded live HEAD/status and the accepted TDVK-03/04 evidence.
- [x] Reconverged the same physical Fe setup independently on `4^3`, `8^3`, and `12^3`.
- [x] Recorded each accepted mesh, moment, EF, radial residual, and LR provenance.
- [x] Evaluated compact Lehmann Gamma-static, certified finite-q-static, and low-finite-omega rows on every mesh.
- [x] Ran the physical response-eta ladder `0.02`, `0.01`, `0.005` Ry at `8^3` without changing GF integration controls.
- [x] Evaluated the complete `response_lmax=4` (`-1`) `spd` product span without truncation.
- [x] Evaluated the `response_lmax=2` reduced diagnostic and labeled it explicitly approximate.
- [x] Recorded unpruned/retained product dimensions, all-L rank stability, and direct radial-mesh provenance.
- [x] Ran only the representative TDVK-04 reciprocal-GF spot checks at the selected `8^3` settings.
- [x] Appended the `TDVK-05 Fe numerical convergence` evidence section to `docs/TDDFT_RECIPROCAL_VALIDATION.md`.
- [x] Returned `PASS CANDIDATE` without inventing a material convergence threshold.
- [x] Stopped at the mandatory TDVK-05 review gate.
- [x] Committed as `tests: map Fe reciprocal TDDFT convergence`.

## Commit

`tests: map Fe reciprocal TDDFT convergence`


<!-- FILE: 04_REVIEW_GATE_AFTER_TDVK-05.md -->

# Mandatory orchestrator review gate after TDVK-05

## Luna instruction

Do not proceed automatically to TDVK-06.

Return the completed TDVK-03/04/05 evidence to the orchestrator and stop.

## What the orchestrator must decide

1. Is Fe bare `chiKS` sufficiently stable in k mesh and physical eta?
2. Is the full 232-dimensional product representation still the accepted
   complete Fe response space?
3. Is reciprocal GF sufficiently closed against Lehmann to remain a useful
   independent spot-check backend?
4. Are finite-q folding/covariance diagnostics clean?
5. What k mesh, physical eta and q set are approved for interaction physics?
6. Does the current live KXC/GSR/Dyson code already understand the compact
   orthonormal product representation, or must that representation interface be
   certified first?

## Representation rule for the next stage

The bare response is now an ordinary orthonormal product-space operator.

Do not assume that legacy point-space KXC/GSR/Dyson interfaces can consume it
correctly merely because the matrix dimensions compile.

Before material interaction physics, any product-space KXC/GSR/Dyson mapping
must be derived from the accepted response-space algebra and independently
cross-checked on an affordable fixture.

Do not reconstruct the 12,375-dimensional Fe point matrix for production.

## Exit

Only the orchestrator authorizes TDVK-06 and supplies/approves the final
interaction-representation interpretation if the live implementation does not
already certify it.


<!-- FILE: 05_TDVK-06_STATIC_INTERACTIONS.md -->

# TDVK-06 — Fe static ALSDA and independent GSR diagnostics

## Prerequisite

Run only after explicit orchestrator approval at the TDVK-05 review gate.

## Parent objective

At the approved Fe bare-response settings, evaluate the static transverse
interaction physics using two **separate** routes:

1. direct radial ALSDA;
2. independent Lounis/LCMM-style GSR/sum-rule interaction.

No BES/GCR correction is authorized by default.

This is one parent milestone. If a compact-representation interface needs
mechanical certification, keep it inside TDVK-06 unless it exposes a genuine
physics blocker.

## 0. Live representation preflight

Before running material static physics, determine from the live source/docs
whether KXC-01, GSR-01 and the static denominator diagnostics consume the new
orthonormal compact product representation correctly.

If already certified, reuse that implementation and evidence.

If not yet certified:

- do not invent a projection formula;
- follow the live LR-04/product-space algebra approved by the orchestrator;
- add the minimal compact operator mapping needed;
- compare it against the legacy point-space operator action on an affordable
  fixture;
- require agreement at the existing response-regression precision;
- do not build a production-size Fe point matrix.

If the correct mapping is ambiguous, stop:

`BLOCKED — COMPACT INTERACTION REPRESENTATION REQUIRES PHYSICS REVIEW`

## 1. Common Fe input

Use the orchestrator-approved:

- accepted Fe state;
- k mesh;
- physical response eta;
- compact complete response basis;
- Gamma, `omega=0`;
- primary Lehmann bare response.

Keep reciprocal GF only as a previously certified bare-response cross-check.

## 2. Direct ALSDA route

Using the existing direct radial ALSDA/KXC diagnostics, record the raw
rotational/static quantities already defined by the live implementation,
including the established residual and overlap diagnostics.

Record at minimum:

- kernel provenance;
- absolute residual;
- relative residual;
- rigid-vector overlap;
- low-magnetization diagnostics already exposed by the service;
- any conditioning information relevant to the compact mapping.

Do not:

- rescale the kernel;
- tune a factor/sign;
- alter magnetization to improve Goldstone behavior;
- apply BES/GCR correction.

## 3. Independent GSR route

Run the existing independent Lounis/LCMM-style sum-rule service.

Record:

- equation dimension/rank;
- singular-value or conditioning summary;
- condition number if defined;
- equation residual;
- static/Goldstone residual diagnostics;
- rigid overlap;
- full-response/nonspherical blocker status;
- exact `blocked/status` state.

Do not:

- replace GSR by rescaled ALSDA;
- add Tikhonov regularization;
- change the SVD cutoff for a prettier result;
- fit an interaction to force zero residual;
- silently fall back to ALSDA.

If GSR reports BLOCKED, preserve that verdict.

## 4. Static denominator diagnostic

If the current validated implementation already exposes the static Dyson
denominator diagnostic consistently in compact space, record the smallest
candidate eigenvalue/rigid overlap information already defined by the code.

Do not add a new mode-finding algorithm in this milestone.

## Evidence discipline

Keep ALSDA and GSR results in separate tables.

Luna does not decide whether a nonzero ALSDA residual is acceptable physics.

## Documentation

Append:

`TDVK-06 Fe static ALSDA/GSR diagnostics`

to `docs/TDDFT_RECIPROCAL_VALIDATION.md`.

## Return to orchestrator

Return:

- representation-certification status;
- ALSDA raw table;
- GSR raw table;
- any static denominator diagnostics;
- exact BLOCKED states;
- no Goldstone repair;
- no physical interpretation.

Then stop.

## Commit

`tests: record Fe static TDDFT interaction diagnostics`


<!-- FILE: 06_TDVK-07_FE_DYSON_LOSS.md -->

# TDVK-07 — first interacting Fe Dyson/loss response

## Prerequisite

Run only after orchestrator approval of TDVK-06.

## Parent objective

Establish whether the validated compact bare response and one explicitly
approved interaction route produce a numerically coherent first Fe interacting
response.

This milestone does **not** authorize stiffness fitting, damping extraction or
automatic magnon assignment.

## 0. Compact Dyson preflight

Confirm from live source/evidence that the Dyson service operates in the same
orthonormal compact product representation as the approved `chiKS` and
interaction operator.

If the compact Dyson algebra has not yet been certified:

- add only the minimal representation plumbing;
- compare compact Dyson output against independently projected legacy point-space
  Dyson output on an affordable fixture;
- preserve the existing Dyson equation, signs and conventions;
- do not reconstruct a production Fe point matrix.

If the mapping is physically ambiguous, stop and return the evidence.

## Route order

### Route 1 — direct ALSDA

Run first if approved by the orchestrator using:

- approved Fe state;
- approved k mesh;
- approved physical eta;
- compact Lehmann bare backend;
- direct ALSDA;
- Goldstone correction OFF.

### Route 2 — GSR

Run only if TDVK-06 GSR was unblocked and explicitly approved.

Keep routes separate.

### Forbidden

Do not enable BES/GCR correction in this milestone unless the orchestrator
explicitly changes the plan.

## q/frequency sampling

Use only TDVK-04/05-approved q values.

Start with:

- Gamma;
- a few small q values in literal direct reciprocal coordinates.

Use a controlled low-energy frequency window broad enough to reveal the raw
response without tuning it to literature.

Record the exact Ry grid and any established meV conversion.

Do not run an unnecessarily dense production spectrum.

## Output

Use the existing validated Dyson/loss products.

For every q/frequency sample record at least:

- interaction route;
- compact `chiKS` diagnostic;
- enhanced-susceptibility diagnostic;
- loss diagnostic already supported by production output;
- Dyson status / near-pole diagnostics already exposed;
- physical eta.

Do not add a new loss-matrix eigenmode tracker here.

Do not call the strongest trace peak a magnon.

## Bare-backend spot check

At one or two already-certified `(q,omega)` points, reciprocal GF may remain
enabled as a bare-response cross-check.

Do not compute a full reciprocal-GF spectrum.

## Fail closed

If the interacting response:

- diverges numerically;
- changes sign unexpectedly;
- violates previously established q/-q behavior;
- has broken Gamma behavior;
- loses consistency with the approved bare/static layers;

preserve the case and stop.

Do not patch:

- KXC;
- Dyson signs;
- circular conventions;
- q phases;
- eta;
- Goldstone shifts.

## Documentation

Append route-separated raw results to:

`docs/TDDFT_RECIPROCAL_VALIDATION.md`

## Return to orchestrator

Return:

- input deck/settings;
- q/frequency grid;
- compact raw response/loss tables;
- paths to raw outputs;
- near-pole diagnostics;
- no magnon/stiffness/damping claim.

Then stop for physics interpretation.

## Commit

`tests: generate first Fe reciprocal TDDFT loss response`


<!-- FILE: 07_TDVK-08_NI_REPLICATION.md -->

# TDVK-08 — fcc Ni reciprocal-framework replication

## Prerequisite

Run only after the orchestrator declares the Fe reciprocal framework
sufficiently trustworthy.

## Objective

Replicate the **already approved Fe methodology** on a newly accepted current
fcc-Ni state without changing physics to obtain a desired spectrum.

This is a second-material replication, not a new implementation campaign.

## 0. Current Ni state

Start from the current tracked/retained fcc-Ni reference deck.

Record:

- exact structure/lattice constant;
- XC/TXC;
- scalar-relativistic and spin flags;
- `lmax`;
- k mesh;
- moment;
- EF;
- temperature;
- LR acceptance/snapshot;
- Pauli/no-SOC validity;
- reciprocal `ham_only` / orthogonal second-order provenance;
- compact product-basis dimension/rank provenance.

Do not reuse pre-purge or legacy Ni spectra as validation.

## Minimal replication ladder

Repeat only the Fe-approved methodology:

1. Gamma `omega=0` compact Lehmann bare response.
2. Small finite-omega compact Lehmann response.
3. One representative Lehmann↔reciprocal-GF bare cross-check using
   TDVK-03-approved GF controls.
4. One finite-q / exact -q diagnostic.
5. Minimal k-mesh and physical-eta sensitivity sufficient to test whether the
   Fe-selected settings are pathological for Ni.
6. Direct ALSDA static diagnostic.
7. GSR static diagnostic only if it was useful/unblocked/approved in Fe.
8. First Dyson/loss response only after the static evidence has been returned to
   and approved by the orchestrator.

Stop at the first failed gate.

## No parameter rescue

Do not:

- tune the Ni moment;
- change XC;
- alter eta to match literature;
- rescale KXC;
- change circular channel because it looks better;
- enable BES/GCR by default;
- truncate product modes for speed.

## Documentation

Create a clearly separate Ni section in:

`docs/TDDFT_RECIPROCAL_VALIDATION.md`

Do not mix Fe and Ni convergence tables.

## Return to orchestrator

Return the evidence at the first failed/uncertain gate, or the complete minimal
replication ladder if all gates execute cleanly.

Do not claim literature agreement.

## Commit

`tests: replicate reciprocal TDDFT validation on Ni`


<!-- FILE: 08_TDVK-09_EVIDENCE_HANDOFF.md -->

# TDVK-09 — reciprocal TD-DFT evidence handoff

## Nature

Documentation/evidence assembly only.

Do not start new physics calculations unless a previously required evidence item
is genuinely missing.

## Objective

Create an auditable Fe/Ni handoff for orchestrator physics adjudication and later
literature comparison.

## Front-page status table

Update:

`docs/TDDFT_RECIPROCAL_VALIDATION.md`

with a compact table containing at least:

| item | Fe | Ni | evidence class | status |
|---|---|---|---|---|
| accepted ground state | | | material provenance | |
| compact Lehmann Gamma `chiKS` | | | material execution | |
| reciprocal-GF cross-check | | | independent numerical cross-check | |
| mixed-eigensystem GF oracle | | N/A unless repeated | numerical regression | |
| compact GF optimization oracle | | N/A unless repeated | numerical regression | |
| finite-q exact endpoint | | | material execution | |
| q/-q covariance | | | numerical invariant | |
| k-mesh ladder | | | numerical convergence | |
| physical eta ladder | | | numerical convergence | |
| full product-space provenance | | | representation evidence | |
| reduced angular sensitivity | | | numerical sensitivity | |
| direct ALSDA static residual | | | physical diagnostic | |
| independent GSR route | | | interaction diagnostic | |
| compact Dyson certification | | | representation/numerical evidence | |
| first Dyson/loss response | | | material response | |
| low-q branch identified | only if orchestrator approved | only if approved | physical interpretation | |
| literature agreement | NOT PERFORMED unless explicitly done | NOT PERFORMED unless explicitly done | literature | |

## Required appendices/provenance

Include:

1. exact final Git HEAD and worktree status;
2. compiler/build provenance;
3. input decks;
4. accepted moments/EF/k meshes;
5. compact product-basis provenance;
6. backend cross-check tables;
7. GF quadrature/performance evidence;
8. k-mesh/eta/angular tables;
9. finite-q/covariance evidence;
10. ALSDA/GSR diagnostics;
11. compact Dyson certification evidence;
12. paths to raw spectra;
13. tests and PASS/BLOCKED outcomes;
14. commits made during the campaign.

## Unresolved items

List every unresolved issue without inventing a fix.

Give the owning layer when mechanically known, for example:

- accepted state / SCF;
- transition vertex;
- compact product representation;
- Lehmann;
- reciprocal GF;
- quadrature/performance;
- KXC;
- GSR;
- Dyson;
- production orchestration.

If ownership is unclear:

`UNATTRIBUTED — requires physics review`

## Claims discipline

Do not claim any of the following unless the orchestrator has explicitly
approved the conclusion from the accumulated evidence:

- magnon stiffness;
- magnon damping;
- physical mode assignment;
- literature agreement;
- Goldstone-corrected validation;
- native-RSGF validation;
- exact scalar-relativistic response;
- fully relativistic response.

## Archive

Create a compact evidence archive containing:

- final validation document;
- material input decks;
- compact tabular diagnostics;
- scripts added specifically for validation.

Do not archive huge matrix dumps, build products or redundant SCF scratch files.

## Return to orchestrator

Return:

- archive path;
- final status table;
- unresolved-items list;
- final HEAD.

## Commit

`docs: hand off reciprocal TDDFT validation evidence`
