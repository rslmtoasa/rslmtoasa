<!-- FILE: 00_MASTER_BLUEPRINT.md -->

# TDVAL-K master blueprint — reciprocal collinear transverse TD-DFT

## 0. Authority and scope

Work only on:

- repository: `https://github.com/rslmtoasa/rslmtoasa`
- branch: `fable_v4`
- pinned starting state for this pack: `795b1e488d`

Before each slice, inspect the live branch.  If the branch has advanced, preserve
the scientific contracts below and report any conflict instead of silently
rewriting them.

This campaign validates only the certified initial response tuple:

- collinear magnetic ground state;
- no SOC;
- Pauli/no-SOC transverse response;
- scalar-relativistic LMTO ground state;
- orthogonal `ham_only`;
- second-order/HOH LMTO;
- `sp` / `spd`;
- full site × response-harmonic × radial response space;
- no site-only replacement.

Do not describe this as an exact scalar-relativistic or fully relativistic
response.

## 1. Deferred work

The following are **out of scope for every Luna slice in this pack**:

- native real-space GF TD-DFT production integration;
- TDRUN-02 / `native_rsgf`;
- recursion/Chebyshev response convergence;
- noncollinear response;
- SOC response;
- Sternheimer/Savrasov implementation;
- exact scalar-relativistic nonspherical response;
- fully relativistic `j,kappa` response;
- BdG;
- new XC physics.

Do not modify the native-RSGF modules merely because they exist.

## 2. Authoritative response physics

Do not rederive these equations in implementation tasks.

### 2.1 Bare reciprocal response

For the certified circular sector,

\[
\chi^{KS}_{IJ}(\mathbf q,\omega)
=
\frac{1}{N_k}\sum_{\mathbf k,nm}
\frac{f_{n\mathbf k}-f_{m,\mathbf k+\mathbf q}}
{\omega+\epsilon_{n\mathbf k}-\epsilon_{m,\mathbf k+\mathbf q}+i\eta}
T_{nm;I}T^*_{nm;J},
\]

with the LR-03 circular normalization already implemented by LR-06.

The endpoint is exactly `k+q`, folded by the existing reciprocal-lattice rule.
No nearest-mesh substitution is allowed.

The production bare-response routes are:

1. `backend='lehmann'` — primary transparent baseline;
2. `backend='reciprocal_gf'` — independent reciprocal GF energy-integral
   cross-check.

They are representations of the same bare KS response, not distinct physical
models.

### 2.2 Reciprocal GF numerical parameters

The reciprocal-GF route has two numerical controls distinct from the physical
response broadening:

- `gf_integration_points`: Simpson quadrature resolution;
- `gf_integration_eta`: spectral-discontinuity width.

If `gf_integration_eta <= 0`, the current service chooses `eta/40`.  The
physical response broadening remains `eta`.

Do not tune `eta` to improve backend agreement.

### 2.3 Response-space representation

Use the LR-04 canonical right-weighted representation already implemented.

No task may insert an additional radial weight.

For `spd`, the formally complete response product cutoff is:

\[
L_{\max}^{response}=4.
\]

A smaller `response_lmax` may be used only as an explicitly labelled
approximation/sensitivity test.  It must never be called exact.

The direct production radial mesh is the reference.  Do not invent radial
decimation in this campaign.

### 2.4 Direct ALSDA route

Use the existing KXC-01 kernel unchanged.

The raw static diagnostic is

\[
b_{xc}=K_{xc}g,\qquad
r=\chi^{KS}(0,0)b_{xc}-g.
\]

Record the LR-04 metric residual and rigid-vector overlap.  Do not rescale
`Kxc` from this residual.

### 2.5 Independent Lounis/LCMM route

Use the existing GSR-01 service unchanged.

It is a separate interaction construction.  It must never be implemented by
rescaling direct ALSDA.

If GSR reports rank deficiency, nonspherical leakage, or `BLOCKED`, preserve
that result and stop for physics review.

### 2.6 Dyson/loss

Use the existing TDDY-01 equation unchanged:

\[
(I-\chi^{KS}K)\chi=\chi^{KS}.
\]

The loss convention is

\[
L=-\frac{\chi-\chi^{\ddagger_W}}{2i\pi}.
\]

Do not negate a circular channel or loss quantity to make a spectrum look
positive.

## 3. Claims discipline

Every result must be labelled as one of:

1. algebraic consistency;
2. independent numerical cross-check;
3. numerical convergence;
4. converged-material validation;
5. literature agreement.

Never promote one class into another.

In particular:

- Lehmann/GF agreement does not independently validate the shared radial vertex;
- a passing GSR route does not validate direct ALSDA;
- a later Goldstone-corrected route cannot rescue a broken raw route;
- a visually plausible magnon curve is not a validation by itself.

## 4. Material order

### First material: bcc Fe

Use the existing tracked bcc-Fe ground-state deck as the starting material
definition.  Do not alter lattice constant, XC choice, magnetic setup, or
moments to improve TD-DFT behavior.

The previously accepted Fe foundation uses the current scalar-relativistic
collinear setup and a production-style `spd` basis.  Re-establish the exact live
provenance before response work.

### Second material: fcc Ni

Do not start Ni until Fe has passed the orchestrator checkpoint after the
Dyson/loss stage.

## 5. Validation order

The required order is:

```text
A. implementation preflight
B. accepted Fe ground-state handoff
C. Fe bare chiKS, Gamma
D. Lehmann ↔ reciprocal-GF closure
E. finite-q / folding / covariance
F. k-mesh and eta convergence
G. static ALSDA + independent GSR diagnostics
H. interacting Fe loss
I. Ni replication
J. literature comparison by the orchestrator
```

Do not jump directly to a magnon dispersion.

## 6. Diagnostic q points

Luna may use direct reciprocal-coordinate q vectors for invariants, but must
not assign high-symmetry labels unless an existing code path or input file
already establishes that label.

Recommended diagnostic sequence:

- `q = (0,0,0)`;
- one mesh-commensurate small q;
- its exact `-q`;
- one arbitrary/off-mesh q;
- later, a small-q sequence chosen by the orchestrator.

These are initially **coordinate diagnostics**, not a claimed Gamma-H,
Gamma-N, or other physical symmetry line.

## 7. Broadening discipline

Use the established diagnostic value `eta=0.02 Ry` for the earliest material
smoke checks because that value already exists in the LR-06 evidence envelope.

Later convergence uses a prescribed ladder, for example:

- `0.02 Ry`;
- `0.01 Ry`;
- `0.005 Ry`.

Do not infer an `eta -> 0` extrapolation and do not subtract eta from a width.

## 8. Backend comparison metric

When comparing Lehmann and reciprocal GF on the **same response space and
accepted state**, compute for every `(q,omega)`:

\[
\Delta=\chi^{KS}_{L}-\chi^{KS}_{GF},
\]

\[
d_F=\sqrt{\sum_{IJ}|\Delta_{IJ}|^2},
\]

\[
r_F=
\frac{d_F}
{\max(\|\chi^{KS}_{L}\|_F,\|\chi^{KS}_{GF}\|_F,\epsilon)},
\]

and

\[
d_\infty=\max_{IJ}|\Delta_{IJ}|.
\]

These are representation-level diagnostics because both matrices use the same
LR-04 canonical representation.

Do not impose an arbitrary material PASS threshold in Luna.  Report the values
and their change with GF quadrature controls.  The orchestrator adjudicates.

## 9. Goldstone correction policy

Do **not** enable GCR-01/BES correction anywhere in the initial reciprocal
campaign.

It may only be reconsidered after:

- Fe Lehmann response is numerically converged;
- the raw direct-ALSDA denominator has an identifiable isolated candidate mode;
- rigid overlap and isolation meet the existing GCR-01 contract;
- at least two convergence samples retain that character.

That decision belongs to the orchestrator, not Luna.

## 10. Literature policy

Luna must not tune calculations against literature and must not decide whether
a numerical dispersion is physically correct.

The orchestrator will perform the literature comparison after receiving the
Fe/Ni evidence pack.

## 11. Required source contracts to read before coding

Read the live versions of:

- `docs/LR_TDDFT_CONVENTIONS.md`
- `docs/LR_RESPONSE_SPACE_ALGEBRA.md`
- `docs/LR_PAULI_TRANSITION_VERTEX.md`
- `docs/LR_KS_SUSCEPTIBILITY.md`
- `docs/LR_GF_SUSCEPTIBILITY_CROSSCHECK.md`
- `docs/KXC_ALSDA_TRANSVERSE_KERNEL.md`
- `docs/GOLDSTONE_SUMRULE_INTERACTION.md`
- `docs/GOLDSTONE_EIGENVALUE_CORRECTION.md`
- `docs/TDDFT_DYSON_AND_LOSS.md`
- `docs/TDDFT_PRODUCTION_DRIVER.md`
- `docs/TDDFT_COLLINEAR_REVALIDATION.md`

Historical `BLOCKED` wording caused by native-RSGF production registration is
not a blocker to this **reciprocal-only** campaign.


<!-- FILE: 01_TDVK-00_PREFLIGHT.md -->

# TDVK-00 — reciprocal validation preflight

## Nature

Audit/documentation slice.  No new physics and normally no source-code change.

## Objective

Establish that the current `fable_v4` branch is ready for a reciprocal-only
material validation campaign using:

- `lehmann`;
- `reciprocal_gf`.

Native RSGF is explicitly deferred.

## Instructions

1. Check out `fable_v4`, fetch, and record:
   - exact HEAD;
   - clean/dirty status;
   - compiler;
   - CMake;
   - build type;
   - libXC status;
   - MPI/OpenMP status;
   - executable hash if a binary exists.

2. Read the live evidence documents listed in `00_MASTER_BLUEPRINT.md`.

3. Verify mechanically that:
   - TDRUN-01 accepts `lehmann`;
   - TDRUN-01 accepts `reciprocal_gf`;
   - TDRUN-01 does not require `native_rsgf`;
   - accepted LR-01 snapshots are handed to the response driver;
   - exact `k+q` endpoint snapshots are constructed.

4. Run the focused existing tests for:
   - LR-06;
   - LR-GF-02;
   - KXC-01;
   - GSR-01;
   - TDDY-01;
   - TDRUN-01 production driver and capability guards.

   Reuse the current CTest target names from the repository.  Do not invent
   replacement tests if one target name has changed; report the live name.

5. Create or update:

   `docs/TDDFT_RECIPROCAL_VALIDATION.md`

   Add a new top section called `TDVAL-K preflight` containing:
   - provenance;
   - exact tests run;
   - pass/fail counts;
   - current production backends;
   - statement that native RSGF is deferred;
   - statement that no Fe/Ni material validation is yet claimed.

6. Do not alter any response equation, kernel, vertex, GF formula, or Goldstone
   logic.

## PASS

PASS only if the reciprocal prerequisites and TDRUN-01 focused tests pass.

If any prerequisite fails:
- document the failure;
- do not repair it in this slice;
- stop and return the exact failure.

## Return to orchestrator

Return:
- HEAD;
- test table;
- changed files;
- final preflight verdict.

Then stop.

## Commit

`docs: start reciprocal TDDFT validation`


<!-- FILE: 02_TDVK-01_BACKEND_CROSSCHECK_DIAGNOSTIC.md -->

# TDVK-01 — add a reciprocal-backend cross-check diagnostic

## Nature

Small orchestration/diagnostic implementation.  No response physics.

## Objective

Allow one prepared production state to evaluate **both** reciprocal bare
response backends and report their full LR-04 matrix difference without
changing the user-selected production result.

The selected backend remains authoritative for KXC/Dyson.  The second backend
is diagnostic only.

## Required behavior

Add an opt-in configuration flag with a clear name such as:

`reciprocal_backend_crosscheck = .false.`

Default must be `.false.` so existing behavior is bitwise/semantically
unchanged.

When enabled, for every prepared `(q,omega)`:

1. evaluate LR-06 Lehmann `chiKS`;
2. evaluate LR-GF-02 reciprocal-GF `chiKS`;
3. do **not** average them;
4. do **not** replace the selected backend;
5. compute:
   - Frobenius norm of each matrix;
   - Frobenius norm of the difference;
   - relative Frobenius difference;
   - maximum absolute element difference.

Use exactly:

```text
delta   = chi_lehmann - chi_reciprocal_gf
dF      = sqrt(sum(abs(delta)**2))
normL   = sqrt(sum(abs(chi_lehmann)**2))
normGF  = sqrt(sum(abs(chi_reciprocal_gf)**2))
relF    = dF / max(normL, normGF, tiny)
dInf    = maxval(abs(delta))
```

Store/report these diagnostics per `(q,omega)`.

These metrics are allowed because both matrices are in the same canonical LR-04
representation.  Do not introduce a new physical metric or radial weighting.

## Implementation constraints

- Reuse `evaluate_lr_ks_susceptibility`.
- Reuse `evaluate_lr_gf_susceptibility`.
- Reuse the already prepared `left_state`, exact `k+q` endpoint, response space,
  radial bases, eta, and channel.
- Use the existing GF controls:
  `gf_integration_points`, `gf_integration_eta`, `gf_energy_margin`.
- Do not copy susceptibility equations into the driver.
- Do not touch native-RSGF code.
- Do not enable GCR-01.
- Do not add an acceptance threshold.

Prefer a small helper so the existing selected-backend path remains readable.

## Tests

Extend the existing TDRUN-01 focused fixture.

Required checks:

1. flag off:
   - existing selected-backend matrices unchanged;
   - no cross-check result is marked valid.

2. flag on:
   - both reciprocal services are called on the same prepared state;
   - diagnostics are finite and nonnegative;
   - independently compute the same norms in the test and compare;
   - selected backend's `ks_susceptibility`, enhanced response, and loss are
     unchanged by enabling the diagnostic.

3. parser/default:
   - absent flag => `.false.`.

No Fe/Ni physics test belongs in this slice.

## Documentation

Update `docs/TDDFT_PRODUCTION_DRIVER.md` with a short subsection explaining that
the flag is a reciprocal **validation diagnostic**, not a third backend and not
a production route selector.

Append implementation evidence to `docs/TDDFT_RECIPROCAL_VALIDATION.md`.

## STOP rule

If adding the diagnostic would require altering LR-06 or LR-GF-02 equations,
stop and report the blocker instead.

## Return to orchestrator

Return:
- diff summary;
- tests run;
- changed files;
- one example diagnostic from the tiny fixture.

Then stop.

## Commit

`tests: add reciprocal TDDFT backend crosscheck`


<!-- FILE: 03_TDVK-02_FE_REFERENCE_GAMMA.md -->

# TDVK-02 — bcc Fe accepted-state and Gamma bare-response smoke

## Nature

Material execution slice.  Do not change response physics.

## Objective

Prove that the real tracked bcc-Fe calculation can reach the rebuilt reciprocal
TD-DFT path and produce a finite **bare Lehmann response** at Gamma.

This is still a material smoke/provenance result, not a magnon validation.

## Ground-state source

Start from the live tracked bcc-Fe input under:

`tests/scf/cases/bulk/bccFe`

Do not alter physical ground-state parameters to improve response behavior.

Record the live values of:
- lattice/structure;
- `lmax`;
- XC/TXC;
- scalar-relativistic/spin flags;
- moment;
- Fermi level;
- temperature;
- reciprocal mode;
- Hamiltonian order;
- k mesh;
- LR-01 accepted snapshot status;
- SR→Pauli provenance if available.

## Response request

Use:

- backend: `lehmann`;
- channel: `chi_plus`;
- q: `(0,0,0)`;
- response cutoff: complete product space (`response_lmax=-1`);
- Goldstone correction: off;
- interaction route may remain `direct_alsda`, but this slice evaluates the
  **bare chiKS evidence only**;
- initial response eta: `0.02 Ry`;
- first frequency: `omega=0`.

After the static point succeeds, add only a small diagnostic set such as:

`omega = 0.00, 0.02, 0.05 Ry`

Do not create a dense spectrum.

Set `write_full_matrix=.false.` unless a tiny response dimension makes full
output demonstrably safe.  The complete matrix remains in the driver result;
do not force multi-GB text output.

## Execution order

1. Run the accepted SCF/reference calculation.
2. Confirm TDRUN-01 capability gate accepts the state.
3. Run Gamma, omega=0, Lehmann.
4. Confirm:
   - no NaN/Inf;
   - initialized result;
   - actual q and omega are recorded;
   - accepted Fermi level is unchanged;
   - accepted LR-01 snapshot is unchanged;
   - response dimension and `response_lmax` are recorded.
5. If successful, run the small finite-omega set.

## Evidence

Append a `TDVK-02 bcc Fe Gamma smoke` section to
`docs/TDDFT_RECIPROCAL_VALIDATION.md`.

Record timings and memory if easily available, but do not optimize yet.

## Do not

- run native RSGF;
- change XC;
- change moments;
- change signs/factors;
- switch to `chi_minus` because it looks better;
- tune eta;
- enable Goldstone correction;
- interpret a peak as a magnon.

## PASS

PASS means only:

> a real accepted bcc-Fe state traverses the current reciprocal production
> lifecycle and returns finite LR-06 bare response matrices at Gamma.

It does not mean the response is physically correct.

## Return to orchestrator

Return:
- exact input diff relative to tracked Fe case;
- accepted moment and EF;
- response dimension;
- q/omega/eta;
- compact response diagnostics already emitted by the driver;
- runtime;
- PASS/BLOCKED.

Then stop.  Do not proceed to backend closure until approved.

## Commit

`tests: establish bcc Fe reciprocal TDDFT smoke`


<!-- FILE: 04_TDVK-03_FE_BACKEND_CLOSURE.md -->

# TDVK-03 — bcc Fe Lehmann ↔ reciprocal-GF bare-response closure

## Prerequisite

Run only after orchestrator approval of TDVK-02.

## Nature

Independent numerical cross-check on a real Fe material state.

## Objective

Compare LR-06 and LR-GF-02 on the **same accepted bcc-Fe state**, same response
space, same q, same omega, same eta, same occupations, and exact same k+q
endpoint.

No KXC/Dyson interpretation is part of this task.

## Start small

Use Gamma first:

- q = `(0,0,0)`;
- channel = `chi_plus`;
- `omega=0`;
- `eta=0.02 Ry`;
- complete `response_lmax=-1`.

Enable the TDVK-01 reciprocal-backend cross-check diagnostic.

## Reciprocal-GF numerical ladder

The reciprocal GF has numerical integration error independent of physical
response eta.

At fixed physical `eta=0.02 Ry`, run at least three odd Simpson resolutions,
chosen around the existing default, e.g.:

- 1001;
- 2001;
- 4001.

Keep `gf_energy_margin` fixed.

First use the service's documented automatic spectral broadening
(`gf_integration_eta <= 0`, currently `eta/40`).

If runtime permits, add one explicit spectral-width variation at fixed
quadrature resolution.  It must remain smaller than response eta.

Do not vary physical eta in this slice.

## Record for every sample

For the full LR-04 matrices record:

- `normL`;
- `normGF`;
- `dF`;
- `relF`;
- `dInf`;
- GF integration points;
- effective GF integration eta;
- wall time.

Do not compare only site-projected quantities.

## Finite omega

Only if Gamma static executes correctly, repeat at one finite diagnostic
frequency already used in TDVK-02, preferably `omega=0.02 Ry`.

Do not create a dense frequency sweep.

## No automatic PASS threshold

Do not invent a tolerance for real Fe.

The purpose is to expose:
- magnitude of the backend difference;
- sensitivity to Simpson resolution;
- sensitivity to spectral-discontinuity width.

The orchestrator will decide whether the observed closure is acceptable.

## Failure behavior

If backend differences are unexpectedly large, nonconvergent, or unstable:

- do not alter either susceptibility implementation;
- do not tune eta;
- do not modify the transition vertex;
- record the full numerical ladder;
- stop.

## Documentation

Append `TDVK-03 reciprocal backend closure` to
`docs/TDDFT_RECIPROCAL_VALIDATION.md`.

Clearly label this as an independent representation-level numerical
cross-check, not independent validation of the shared radial vertex.

## Return to orchestrator

Return a compact table of all samples and no physical interpretation.

Then stop.

## Commit

`tests: crosscheck Fe reciprocal TDDFT backends`


<!-- FILE: 05_TDVK-04_FE_FINITE_Q_COVARIANCE.md -->

# TDVK-04 — Fe finite-q endpoint, folding, and covariance evidence

## Prerequisite

Run only after orchestrator approval of the Gamma backend evidence.

## Nature

Finite-q production diagnostics.  No interaction/kernel changes.

## Objective

Verify that the actual Fe production path behaves consistently for exact
`k+q` endpoints, reciprocal folding, and q↔-q diagnostics before any magnon
dispersion is attempted.

## q set

Use direct reciprocal coordinates and report them literally.

Include:

1. `q0 = (0,0,0)`;
2. one small mesh-commensurate q;
3. exact `-q` partner;
4. one arbitrary/off-mesh q.

Do **not** attach a crystallographic symmetry label unless the live reciprocal
path already defines it unambiguously.

Use the same:
- accepted Fe state;
- response basis;
- channel;
- eta;
- omega point(s).

Start with one static or low finite omega point, not a dense sweep.

## Backend order

1. Lehmann first.
2. Reciprocal-GF only at representative q points after Lehmann succeeds.

Use the existing exact endpoint adapter.  Do not replace `k+q` with nearest
mesh points.

## Required diagnostics

For every q:
- actual supplied q;
- folded endpoint metadata;
- selected backend;
- response norm/trace diagnostics;
- if cross-check enabled, `dF`, `relF`, `dInf`.

For `q` and `-q`, reuse the **existing LR-03/LR-06 covariance convention**.
Do not derive a new conjugation/channel identity.

If production code does not expose enough information to apply the existing
covariance helper/test pattern, add only a validation helper that reuses the
same established formula from the LR-06/TDDY tests.

## Spin reversal

Do not create a new material spin-reversed SCF state in this slice.
The existing material-independent spin-reversal unit oracle is sufficient for
now.

## PASS boundary

This slice can establish:
- arbitrary-q production execution;
- exact-folding lifecycle behavior;
- material q/-q numerical covariance evidence.

It cannot establish a physical magnon branch.

## Failure behavior

Any odd-in-q or covariance anomaly must be reported, not repaired.

Do not:
- add a phase by hand;
- conjugate eigenvectors manually;
- swap channels ad hoc;
- change the Fourier sign.

## Documentation

Append the q table and covariance diagnostics to
`docs/TDDFT_RECIPROCAL_VALIDATION.md`.

## Return to orchestrator

Return q vectors and numerical covariance/backend tables.

Then stop.

## Commit

`tests: validate Fe finite-q reciprocal TDDFT`


<!-- FILE: 06_TDVK-05_FE_CONVERGENCE.md -->

# TDVK-05 — Fe numerical convergence ladder

## Prerequisite

Run only after finite-q/covariance approval.

## Nature

Numerical convergence campaign.  No physics modifications.

## Objective

Measure numerical sensitivity of the reciprocal Fe response before static
Goldstone or Dyson conclusions are promoted.

## Ground-state rule

Do not mix different physical ground-state definitions in one ladder.

For each k mesh, reconverge the same physical Fe setup and record:
- moment;
- EF;
- LR-01 acceptance;
- exact k mesh.

Do not manually force one mesh's EF into another mesh.

## A. k-mesh ladder

Use at least three systematically increasing meshes around the established
8^3 reference.

Preferred if supported directly by the live reciprocal input:

- 4^3: diagnostic/coarse;
- 8^3: established reference;
- 12^3: finer.

If the code's input parameterization does not map transparently to these
meshes, do not guess.  Use the existing established mesh-control syntax and
report the actual generated meshes.

At each mesh run Lehmann `chiKS` for a **small fixed set**:
- Gamma static;
- one previously certified finite q;
- one low finite omega.

Track:
- accepted moment;
- EF;
- compact full-response diagnostics;
- runtime.

Do not run the expensive reciprocal-GF route on every mesh unless requested by
the orchestrator.  One representative backend cross-check per coarse/reference
setting is enough at this stage.

## B. eta ladder

At one selected k mesh fixed by the orchestrator evidence, run:

- `eta=0.02 Ry`;
- `eta=0.01 Ry`;
- `eta=0.005 Ry`.

Use identical q/frequency points.

Do not extrapolate eta to zero.

## C. angular response cutoff

For `spd`:

- complete reference: `response_lmax=4` via `-1` or explicit 4;
- optional sensitivity sample: `response_lmax=2`.

Label `Lmax=2` as **reduced/approximate**.

Do not alter the orbital basis and do not claim Lmax=2 is exact.

## D. radial mesh

The direct accepted production radial mesh is the reference.

Do not implement or test radial decimation in this campaign.

Record radial mesh identity/provenance so that later runs can be compared.

## Stop conditions

If changes with k mesh or eta are large enough that no stable trend is evident:
- do not proceed to Goldstone/Dyson physics;
- preserve the table;
- stop for orchestrator review.

No arbitrary convergence threshold is to be invented by Luna.

## Documentation

Append:
- k-mesh table;
- eta table;
- angular-cutoff sensitivity;
- radial-mesh provenance

to `docs/TDDFT_RECIPROCAL_VALIDATION.md`.

## Return to orchestrator

Return tables only, plus any execution failures.

Then stop.

## Commit

`tests: map Fe reciprocal TDDFT convergence`


<!-- FILE: 07_TDVK-06_STATIC_GOLDSTONE_INTERACTIONS.md -->

# TDVK-06 — Fe static ALSDA and Lounis/LCMM diagnostics

## Prerequisite

Run only after explicit orchestrator approval of numerical convergence.

## Nature

Static physical-consistency diagnostics using already implemented interaction
services.

## Objective

At converged-enough bcc Fe settings, evaluate the q=0, omega=0 static
transverse response with:

1. direct radial ALSDA (KXC-01);
2. independent Lounis/LCMM sum-rule interaction (GSR-01).

Do not use BES/GCR correction.

## Common input

Use the same accepted Fe state and the same converged numerical settings for
both interaction routes.

Bare response should use the orchestrator-approved reciprocal backend, normally
Lehmann as primary.

Keep the reciprocal-GF result available only as a bare-response cross-check.

## A. direct ALSDA

Using the existing KXC-01 diagnostic, record:

\[
b_{xc}=K_{xc}g,
\qquad
r=\chi^{KS}(0,0)b_{xc}-g.
\]

Record:
- metric absolute residual;
- metric relative residual;
- rigid-vector overlap;
- kernel provenance;
- active low-magnetization diagnostics already exposed by KXC-01.

Do not rescale the kernel.

## B. GSR-01

Run the existing independent sum-rule service.

Record:
- numerical rank / equation dimension;
- singular values or compact condition summary;
- condition number if defined;
- equation residual;
- metric Goldstone residual;
- relative metric residual;
- rigid overlap;
- nonspherical/full-response blocker status;
- `blocked/status`;
- direct-ALSDA comparison diagnostic if already emitted.

Do not:
- add Tikhonov regularization;
- change SVD cutoff;
- fall back to ALSDA;
- fit U to obtain zero residual.

## C. denominator diagnostic

If the current services already expose the q=0 static Dyson denominator
`D = I - chiKS K`, record its smallest candidate eigenvalue and rigid overlap
using the existing validated conventions.

Do not add Goldstone correction.

If such a diagnostic is not already cleanly exposed, do not invent a new mode
algorithm in this slice.  Report that it is unavailable and stop after A/B.

## Decision boundary

Luna must not decide whether a nonzero ALSDA residual is “acceptable physics”.

If GSR reports `BLOCKED`, preserve that exact verdict.

## Documentation

Append a separate table for ALSDA and GSR to
`docs/TDDFT_RECIPROCAL_VALIDATION.md`.

Never merge the two routes into one result.

## Return to orchestrator

Return the raw static diagnostic table with no tuning and no physical repair.

Then stop.

## Commit

`tests: record Fe static TDDFT Goldstone diagnostics`


<!-- FILE: 08_TDVK-07_FE_DYSON_LOSS.md -->

# TDVK-07 — first interacting Fe Dyson/loss sweep

## Prerequisite

Run only after orchestrator approval of TDVK-06.

## Nature

Controlled first interacting material response.

## Objective

Establish whether the validated reciprocal bare response and an explicitly
selected interaction produce a numerically coherent Fe loss spectrum.

This is not yet permission to fit stiffness or damping.

## Route order

### Route 1 — direct ALSDA

Run first with:
- approved Fe state;
- approved k mesh;
- approved eta;
- `backend='lehmann'`;
- `interaction_route='direct_alsda'`;
- `goldstone_correction=.false.`.

### Route 2 — GSR

Run only if TDVK-06 GSR was unblocked and the orchestrator explicitly approved
it.

Keep results separate.

### Forbidden

Do not run BES/GCR correction in this slice.

## q/frequency sampling

Use only the q points approved by the orchestrator after TDVK-04/05.

Begin with:
- Gamma;
- a few small direct-coordinate q values.

Use a frequency window broad enough to reveal the low-energy response but do
not choose or alter the window to match a literature curve.

Record the exact Ry grid.  If a meV conversion is printed, use the code's
established unit conversion and print both.

## Output

Use existing TDDY-01 products.

For every q and omega record at least:
- selected interaction route;
- `chiKS` compact trace/norm diagnostic;
- enhanced-susceptibility compact diagnostic;
- loss metric trace already supported by production output;
- Dyson status / near-pole diagnostics already exposed;
- eta.

Do not diagonalize the loss matrix or add a new mode tracker in this slice.

Do not call the largest trace peak a magnon unless the orchestrator later
adjudicates it.

## Backend spot check

At one or two already certified `(q,omega)` points, the reciprocal-GF bare
cross-check may remain enabled.

Do not repeat an expensive full reciprocal-GF spectrum.

## Failure behavior

If the response:
- diverges numerically;
- changes sign unexpectedly;
- loses q/-q consistency;
- shows no low-energy structure;
- has a broken Gamma behavior;

record it and stop.

Do not patch KXC, Dyson, eta, circular signs, or q phases.

## Documentation

Append route-separated raw spectrum tables/paths to
`docs/TDDFT_RECIPROCAL_VALIDATION.md`.

## Return to orchestrator

Return:
- input decks;
- q/frequency grid;
- compact response/loss tables;
- paths to raw outputs;
- any near-pole diagnostics.

Then stop for physics interpretation.

## Commit

`tests: generate first Fe reciprocal TDDFT loss sweep`


<!-- FILE: 09_TDVK-08_NI_REPEAT.md -->

# TDVK-08 — fcc Ni reciprocal-framework replication

## Prerequisite

Do not run until the orchestrator explicitly says the Fe reciprocal framework is
sufficiently trustworthy.

## Nature

Second-material replication, not a new implementation campaign.

## Objective

Determine whether the reciprocal TD-DFT workflow that passed Fe also executes
consistently for fcc Ni without changing physics.

## Ground state

Start from the current tracked/retained fcc-Ni reference deck already cited by
the TDVAL documentation.

Before response work, establish a **new accepted current state** and record:

- exact structure/lattice constant;
- XC/TXC;
- scalar-relativistic and spin flags;
- lmax;
- k mesh;
- moment;
- EF;
- temperature;
- LR-01 accepted snapshot;
- Pauli-basis validity;
- reciprocal `ham_only`, orthogonal second-order provenance.

Do not reuse a legacy Ni response or old pre-purge spectrum as validation.

## Minimal replication ladder

Repeat only the already approved Fe methodology:

1. Gamma, omega=0, Lehmann bare response.
2. Small finite omega, Lehmann.
3. One representative Lehmann ↔ reciprocal-GF cross-check.
4. One finite-q / -q diagnostic.
5. Minimal k/eta sensitivity needed to show the Fe-selected settings are not
   pathological for Ni.
6. Direct ALSDA static diagnostic.
7. GSR static diagnostic only if it was a useful approved route in Fe.
8. First Dyson/loss sweep only after the static diagnostics are returned to and
   approved by the orchestrator.

If any early stage fails, stop there.  Do not continue simply to produce a
spectrum.

## No parameter rescue

Do not:
- tune the Ni moment;
- change XC;
- alter eta to mimic literature;
- rescale KXC;
- choose a different channel because its spectrum looks better;
- enable GCR correction.

## Documentation

Create a clearly separate Ni section in
`docs/TDDFT_RECIPROCAL_VALIDATION.md`.

Do not mix Fe and Ni convergence tables.

## Return to orchestrator

Return the Ni evidence at the first failed/uncertain gate, or the complete
minimal ladder if all gates execute cleanly.

## Commit

`tests: replicate reciprocal TDDFT validation on Ni`

## TDVK-08 execution checklist

- [x] Verified live accepted Fe history through the TDVK-07 controlled GF closure.
- [x] Audited the current tracked fcc-Ni reference deck(s) and recorded their live provenance.
- [x] Tracked Ni reference satisfies the accepted no-SOC, orthogonal, ham_only, second-order/HOH TDDFT contract after the authorized `nsp=1` / explicit second-order update; `hoh=.false.` is preserved.
- [x] Fresh self-consistent 8^3 and 12^3 Ni k-space SCF states.
- [x] SCF-owned EF, occupations, moments, residuals, and mesh fingerprints.
- [x] Direct accepted-state SCF -> TDDFT continuity.
- [x] Complete compact product basis and rank/SVD evidence.
- [x] Complete 8^3/12^3 bare-response operator comparison.
- [x] 12^3 eta=.005 bare-response diagnostic.
- [x] Controlled reciprocal-GF full-matrix Ni spot check.
- [x] Finite-q bare response and established q/-q covariance.
- [x] Direct ALSDA static diagnostic and compact denominator diagnostics.
- [x] Raw Gamma, +q, and -q compact Dyson/loss spectrum.
- [x] Interacting q/-q covariance.
- [x] No Ni parameter rescue, Goldstone/BES/GCR, GSR regularization, mode assignment, or literature comparison.
- [x] TDVK-09 not started.

TDVK-08 checklist verdict:
`PASS CANDIDATE`

## TDVK-08 Ni 16³ convergence closure

Run only the narrow fcc-Ni material-convergence closure after the parent
TDVK-08 replication. Keep the accepted 12³ state as the full methodology
replication state. Do not repeat finite-q, reciprocal-GF, ALSDA, Dyson, loss,
or interacting covariance layers at 16³ unless the closure is materially
misleading and the orchestrator explicitly selects a follow-up.

- [x] Fresh 16³ self-consistent k-space SCF.
- [x] SCF-owned EF and occupations.
- [x] Same-state TDDFT handoff.
- [x] Complete product basis.
- [x] Gamma eta=.01 full matrix.
- [x] 12³→16³ full operator comparison.
- [x] Gamma eta=.005 diagnostic on the same accepted state.
- [x] Convergence sequence assessed against 8³→12³.
- [x] No unnecessary downstream reruns.

TDVK-08 16³ convergence closure:
`PASS CANDIDATE`

TDVK-09 remains not started.


<!-- FILE: 10_TDVK-09_EVIDENCE_HANDOFF.md -->

# TDVK-09 — reciprocal TD-DFT evidence handoff

## Nature

Documentation/evidence assembly only.

## Objective

Produce a compact, auditable handoff for physics adjudication by the
orchestrator.

Do not add new calculations unless an already requested output is missing.

## Required summary

Update:

`docs/TDDFT_RECIPROCAL_VALIDATION.md`

with a front-page status table containing these rows:

| item | Fe | Ni | evidence class | status |
|---|---|---|---|---|
| accepted ground state | | | material provenance | |
| Lehmann Gamma chiKS | | | material execution | |
| reciprocal-GF cross-check | | | independent numerical cross-check | |
| finite-q exact endpoint | | | material execution | |
| q/-q covariance | | | numerical invariant | |
| k-mesh ladder | | | numerical convergence | |
| eta ladder | | | numerical convergence | |
| angular cutoff sensitivity | | | numerical sensitivity | |
| direct ALSDA static residual | | | physical diagnostic | |
| GSR static route | | | independent interaction diagnostic | |
| Dyson/loss | | | material response | |
| low-q branch identified | NO unless orchestrator approved | NO unless orchestrator approved | physical interpretation | |
| literature agreement | NOT PERFORMED | NOT PERFORMED | literature | |

## Required appendices

Include:

1. exact final Git HEAD and worktree status;
2. compiler/build provenance;
3. all input decks used;
4. accepted moment/EF for each material and mesh;
5. backend cross-check tables;
6. k/eta/angular tables;
7. static ALSDA/GSR diagnostics;
8. links/paths to raw spectrum outputs;
9. test commands and PASS/BLOCKED results;
10. list of all commits made by this campaign.

## Explicit unresolved-items section

List every unresolved issue without proposing a physics fix.

For each, state the owner layer if mechanically known:

- accepted state / SCF;
- transition vertex;
- LR-06 Lehmann;
- LR-GF-02;
- response metric;
- KXC-01;
- GSR-01;
- TDDY-01;
- production orchestration.

If ownership is not clear, write `UNATTRIBUTED — requires physics review`.

## Claims

Do not claim:
- magnon stiffness;
- damping;
- literature agreement;
- Goldstone-corrected validation;
- native-RSGF validation;

unless the orchestrator explicitly supplied that conclusion.

## Archive

Create a small archive containing:
- the final validation document;
- material input decks;
- compact tabular diagnostics;
- scripts added specifically for validation.

Do not archive huge matrix dumps or build products.

## Return to orchestrator

Return:
- archive path;
- final status table;
- unresolved-items list;
- final HEAD.

## Commit

`docs: hand off reciprocal TDDFT validation evidence`


<!-- FILE: README.md -->

# RS-LMTO-ASA reciprocal TD-DFT validation — Luna prompt pack

Branch: `fable_v4`
Pinned starting commit for this pack: `795b1e488d`
Campaign name: `TDVAL-K`
Date: 2026-09-12

## Purpose

This pack deliberately changes the order of work.

The native real-space GF TD-DFT backend is **deferred**.  RSGF capability levels
R0–R2 remain accepted within their documented finite/provider scopes, but
TDRUN-02 / R3 and native-material R4 are not on the present critical path.

The present goal is to establish a trustworthy **reciprocal-space TD-DFT**
framework first:

```text
accepted SCF state
    ↓
H(k), exact k+q endpoints
    ↓
LR-05 Pauli transition vertices
    ↓
LR-06 Lehmann chiKS  ←→  LR-GF-02 reciprocal-GF chiKS
    ↓
KXC-01 direct ALSDA / GSR-01 independent sum-rule route
    ↓
TDDY-01 Dyson + loss
    ↓
bcc Fe validation
    ↓
fcc Ni validation
```

The physics is specified in `00_MASTER_BLUEPRINT.md`.  Luna is an implementation
worker only.  It must not invent response equations, signs, factors, acceptance
thresholds, physical explanations, or literature fits.

## Execution order

Run one slice at a time:

1. `01_TDVK-00_PREFLIGHT.md`
2. `02_TDVK-01_BACKEND_CROSSCHECK_DIAGNOSTIC.md`
3. `03_TDVK-02_FE_REFERENCE_GAMMA.md`

**STOP and return evidence to the orchestrator.**

Only after explicit approval continue:

4. `04_TDVK-03_FE_BACKEND_CLOSURE.md`
5. `05_TDVK-04_FE_FINITE_Q_COVARIANCE.md`
6. `06_TDVK-05_FE_CONVERGENCE.md`

**STOP and return evidence to the orchestrator.**

Only after explicit approval continue:

7. `07_TDVK-06_STATIC_GOLDSTONE_INTERACTIONS.md`

**STOP and return evidence to the orchestrator.**

Only after explicit approval continue:

8. `08_TDVK-07_FE_DYSON_LOSS.md`
9. `09_TDVK-08_NI_REPEAT.md`
10. `10_TDVK-09_EVIDENCE_HANDOFF.md`

## Global Luna rule

When a requested numerical or material check fails:

- do **not** patch the physics;
- do **not** tune signs, factors, XC fields, eta, moments, or kernel scales;
- do **not** switch circular channel because another one “looks better”;
- do **not** enable Goldstone correction to hide a failure;
- record the failure exactly;
- identify the owning layer only if the evidence makes that mechanical;
- stop the slice and return the evidence.

Each prompt supplies a one-line commit message.  Do not combine slices in one
commit unless explicitly asked.
