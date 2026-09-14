# Reciprocal TD-DFT validation

## TDVAL-K preflight

**Campaign:** TDVAL-K<br>
**Date:** 2026-09-12<br>
**Repository branch:** `fable_v4`<br>
**Pinned starting commit:** `795b1e488d702f62637953ac43f2a33110265163d`

### Provenance

The live branch was already `fable_v4` at the pinned commit, and the local
`origin/fable_v4` reference points to the same commit.  No checkout was
performed because the worktree contains existing user changes.  A
`git fetch --dry-run origin fable_v4` was attempted in the sandbox but could
not resolve `github.com`, then retried with network access and reached GitHub
but failed because no SSH public key is available.  This does not affect the
local branch comparison above.

The worktree was dirty before this preflight.  Existing tracked and untracked
TDDFT changes were preserved.  In particular, an optional `native_rsgf` path
is present in those existing changes; it is outside this reciprocal-only
campaign and was not selected, modified, or used as a prerequisite.

| item | observed value |
| --- | --- |
| HEAD | `795b1e488d702f62637953ac43f2a33110265163` |
| worktree | dirty before this documentation slice |
| compiler | GNU Fortran 13.3.0 (`x86_64-linux-gnu-gfortran-13`) |
| CMake | 3.28.3 |
| build directory/type | `build`, `Debug` |
| tests | `BUILD_TESTING=ON` |
| libXC | enabled, version 5.2.3 |
| MPI | disabled (`ENABLE_MPI=OFF`) |
| OpenMP | enabled (`ENABLE_OPENMP=ON`) |
| executable | `build/bin/rslmto.x` |
| executable SHA-256 | `ed127568cc7df9cb02268da1e4a824f232fc777702248ddfa00432838360ece0` |

### Mechanical prerequisite audit

The live source confirms that TDRUN-01:

- accepts `backend='lehmann'` (and its `spectral` alias);
- accepts `backend='reciprocal_gf'`;
- dispatches either reciprocal backend without constructing or requiring a
  native RSGF provider;
- hands accepted LR-01 radial snapshots to the response driver, including
  their Pauli basis validity checks; and
- constructs each endpoint from the exact `k+q` list through the reciprocal
  arbitrary-k eigensystem service, with half-open-BZ folding validated again by
  LR-06.

The current dirty worktree also contains an optional `native_rsgf` dispatch.
That route is explicitly deferred by the TDVAL-K blueprint and is not part of
the reciprocal production-backend verdict.

### Tests

The existing focused CTest targets were rebuilt where needed and run without
failure:

```text
cmake --build build --target UnitLrKsSusceptibility UnitLrGfSusceptibility \
  UnitLrAlsdaKernel UnitLrGoldstoneSumrule UnitTddftDyson \
  UnitTddftProductionDriver -j2

ctest --test-dir build --output-on-failure -R \
  '^(UnitLrKsSusceptibility|UnitLrGfSusceptibility|UnitLrAlsdaKernel($|RejectProvenance$)|UnitLrGoldstoneSumrule|UnitTddftDyson|UnitTddftProductionDriver|UnitTddftProductionDriverRejectSoc|UnitTddftProductionDriverRejectGeneralizedOverlap|TddftProductionDriverSmoke|TddftProductionDriverFeatureOff)$'

python3 tests/validation/val24_lr03_conventions.py
```

| test group | pass | fail | result |
| --- | ---: | ---: | --- |
| LR-06 / `UnitLrKsSusceptibility` | 1 | 0 | PASS |
| LR-GF-02 / `UnitLrGfSusceptibility` | 1 | 0 | PASS |
| KXC-01 / `UnitLrAlsdaKernel` and provenance guard | 2 | 0 | PASS |
| GSR-01 / `UnitLrGoldstoneSumrule` | 1 | 0 | PASS |
| TDDY-01 / `UnitTddftDyson` | 1 | 0 | PASS |
| TDRUN-01 driver and capability guards | 5 | 0 | PASS |
| LR-03 convention oracle | 1 | 0 | PASS |
| **total** | **12** | **0** | **PASS** |

The CTest-focused portion is 11/11 passed.  The convention oracle is a
separate 1/1 passed check, giving 12/12 checks overall.  The test results are
implementation and finite-fixture evidence only.

### Production boundary and claim boundary

The reciprocal TDVAL-K production backends are:

1. `lehmann` — primary transparent bare-response baseline;
2. `reciprocal_gf` — independent reciprocal-GF bare-response cross-check.

Native real-space RSGF is deferred.  No Fe or Ni material validation is
claimed by this preflight: no converged material `chiKS`, interaction, Dyson,
loss, or literature-comparison result is established here.

**Preflight verdict: PASS** for the reciprocal prerequisites and focused tests,
with native RSGF explicitly deferred and the existing dirty-worktree scope
recorded above.

## TDVK-01 reciprocal backend crosscheck diagnostic

**Status:** PASS for orchestration and finite-fixture numerical evidence only.

### Implementation checklist

- [x] Added `reciprocal_backend_crosscheck`, defaulting to `.false.`.
- [x] Reused LR-06 Lehmann and LR-GF-02 reciprocal-GF services on the same
  prepared state and exact endpoint.
- [x] Retained the selected backend as the authoritative KXC/Dyson input.
- [x] Stored the full canonical LR-04 delta and per-`(q,omega)` norm metrics.
- [x] Emitted validation metrics and delta records without imposing a
  material acceptance threshold.
- [x] Added parser/default, flag-off, flag-on, independent-metric, and
  selected-result invariance checks to `UnitTddftProductionDriver`.
- [x] Left LR-06, LR-GF-02, native-RSGF, and Goldstone equations untouched.

### Implementation evidence

The driver now evaluates the configured backend first. When the diagnostic is
enabled, it reuses that selected result when it is already Lehmann or
reciprocal-GF, evaluates only the other reciprocal service, and computes

```text
delta  = chi_lehmann - chi_reciprocal_gf
dF     = sqrt(sum(abs(delta)**2))
normL  = sqrt(sum(abs(chi_lehmann)**2))
normGF = sqrt(sum(abs(chi_reciprocal_gf)**2))
relF   = dF / max(normL, normGF, tiny)
dInf   = maxval(abs(delta))
```

The result retains `reciprocal_crosscheck_delta` in the canonical matrix
layout `(I,J,frequency,q)` and the five scalar diagnostics in
`(frequency,q)` arrays. The output writer records the metrics and each delta
element as comment-prefixed diagnostic records. The flag-off path does not
allocate or mark any crosscheck result valid.

### Focused test evidence

```text
cmake --build build --target UnitTddftProductionDriver -j2
ctest --test-dir build --output-on-failure -R '^UnitTddftProductionDriver$'
```

The test covers the absent parser flag default, flag-off no-op behavior, both
reciprocal service results on the same prepared fixture, independent metric
recalculation, finite/nonnegative diagnostics, and byte-stable selected
Lehmann/Dyson/loss results when the diagnostic is enabled. It prints one tiny
fixture example in the form
`normL normGF dF relF dInf`; the exact values are recorded by the focused test
run:

```text
normL = 7.749380E+02  normGF = 4.675770E+02  dF = 3.074647E+02
relF = 3.967604E-01  dInf = 1.596306E+02
```

These values are diagnostic evidence rather than a physics acceptance
criterion.

No Fe/Ni material validation is claimed by TDVK-01.

## TDVK-02 bcc Fe Gamma smoke

**Status:** BLOCKED at the complete-space LR-06 dense response evaluation.

This slice used the live tracked material files under
[`tests/scf/cases/bulk/bccFe`](../tests/scf/cases/bulk/bccFe) as the source and
did not modify them.  The reproducible response deck is
[`input_fe_reference_gamma.nml`](../tests/integration/tddft_driver_smoke/input_fe_reference_gamma.nml).

### Checklist

- [x] Rebuilt the current reciprocal production path before material execution.
- [x] Preserved the tracked Fe ground-state files and recorded the exact deck
  adaptation below.
- [x] Reached an accepted scalar-relativistic, collinear Fe SCF state and
  passed the TDRUN-01 capability gate.
- [x] Configured `lehmann`, `chi_plus`, Gamma, complete product space,
  `eta=0.02 Ry`, Goldstone off, and no full-matrix text output.
- [x] Kept native RSGF, channel changes, eta tuning, Goldstone repair, and
  response-physics changes out of scope.
- [ ] Complete-space Gamma `omega=0` returned a finite serialized result.
- [ ] NaN/Inf and initialized-result checks could be made on a response output.
- [ ] The post-static `omega=0.02, 0.05 Ry` diagnostic set was run.
- [ ] TDVK-02 PASS.

### Live tracked Fe provenance

| item | tracked/live value or observed result |
| --- | --- |
| structure | one-site bcc Fe; `alat=2.86120`, `wav=1.40880`, `ct(1)=3.0`, `r2=9.00`, `rc=080` |
| basis | `lmax=2` (`spd`) from `Fe.nml`; one response site |
| XC | legacy RS-LMTO, `TXC=1`, Barth-Hedin |
| spin representation | tracked deck was `nsp=2`; response-certified adapter used `nsp=1` (collinear scalar-relativistic, no SOC) |
| Hamiltonian | tracked deck had `hoh=.false.`; response-certified adapter used `hoh=.true.` and `kspace_ham_order='second'` |
| reciprocal state | `ham_only`, LAPACK, `8x8x8=512` Gamma-centered points, no symmetry or time-reversal reduction |
| smearing | tetrahedron DOS, `T=300 K`, fixed-Fermi input (`auto_find_fermi=.false.`) |
| accepted SCF | fresh tracked-Fe run accepted at iteration 37 with `conv_thr=1e-6`; residual `6.903e-7` in the Debug execution |
| accepted moment / EF | `2.267442 mu_B` and approximately `-0.085122 Ry` in the accepted Debug log (Release rerun: `2.267442 mu_B`, `-0.085119 Ry` at displayed precision) |
| LR-01 / SR→Pauli | accepted radial snapshot was handed to the driver; no-SOC Pauli provenance was accepted by the TDRUN-01 gate |

The `nsp` and `hoh` changes are representation requirements of the certified
TDVK reciprocal baseline, not response fitting.  The added SCF
`conv_thr=1e-6` is the established accepted-state criterion used by the live
LR-01 oracle; the strict default (`5e-9`) stalled near `4.5e-8` after 100
steps on this fresh atomic start.  No lattice, XC, magnetic seed, eta,
channel, sign, factor, or kernel scale was changed.

### Exact input adaptation

Relative to the tracked `tests/scf/cases/bulk/bccFe/input.nml`, the adapter:

```text
calculation: add post_processing='tddft'
self:       nstep 1 -> 50; add conv_thr=1.0e-6
control:    nsp 2 -> 1
hamiltonian: hoh .false. -> .true.
add &reciprocal:
  nk1=nk2=nk3=8; zero offsets; symmetry/time-reversal/shift disabled;
  tetrahedron; temperature=300; fixed EF; reciprocal_mode=ham_only;
  reciprocal_backend=lapack; kspace_ham_order=second
add &tddft:
  enabled; channel=chi_plus; q=(0,0,0); eta=0.02; response_lmax=-1;
  interaction_route=direct_alsda; goldstone_correction=.false.;
  backend=lehmann; reciprocal_backend_crosscheck=.false.;
  write_full_matrix=.false.; explicit omega grid=(0.00,0.02,0.05) Ry
```

The original tracked `&lattice`, `&atoms`, `&energy`, and `&mix` values were
otherwise retained verbatim in the adapter.

### Gate and response evidence

The accepted run logged:

```text
Converged!0.0000006903
Generated Monkhorst-Pack mesh with 512 k-points
build_kspace_hamiltonian: K-space Hamiltonian built successfully
Kanpur mapping: reciprocal_mode=ham_only
Kanpur mapping: Hamiltonian-only mode (TB-like).
```

For Fe's 495-point radial mesh, `response_lmax=-1` expands the complete
`lmax=2` product to angular cutoff 4, with

```text
angular harmonics = 1 + 3 + 5 + 7 + 9 = 25
response dimension = 1 site * 25 harmonics * 495 radial points * 1 channel
                  = 12,375
```

One dense complex(kind=8) `12,375 x 12,375` matrix is 153,140,625 elements
and about 2.45 GB decimal.  The current production lifecycle allocates
multiple dense bare, Dyson, and result matrices even when
`write_full_matrix=.false.`; that option suppresses text serialization, not
the in-memory matrices.  The Debug three-frequency attempt and a separate
Release static attempt both reached the response call but produced no
`tddft_fe_gamma.dat`.  The Release static attempt was stopped after roughly
192 seconds in the complete-space response evaluation.  Consequently there
are no response-side finite/NaN diagnostics, output q/omega records, or
post-static finite-frequency evidence to claim.

The owning blocker is the current dense LR-06/production-driver material
allocation for the required complete Fe product space.  Per the campaign
failure rule, no response equation, vertex, eta, kernel, or acceptance physics
was changed to work around it.

### Checks and next action

The existing reciprocal prerequisites remain green:

```text
cmake --build build -j2                         PASS
ctest --test-dir build --output-on-failure -R '^TddftProductionDriverSmoke$'  PASS (1/1)
separate Release build in /tmp/tdvk02-build    PASS
```

TDVK-03 is not started.  The next task should be a narrowly scoped mechanical
review of the dense LR-06/production-driver allocation and a streaming or
matrix-free material-smoke design, followed by a rerun of TDVK-02 after
orchestrator approval.  It must preserve the complete response space and the
current response equations.

## TDVK-02R0 LMTO radial-product response-basis closure audit

**Status:** PASS for algebraic representation closure; no material-validation
claim and no change to the TDVK-02 BLOCKED verdict.

- [x] Exact ordered-pair inventory: `sp=52`, `spd=232` candidates.
- [x] Gram spectra, diagnostic rank/nullity, and condition estimates recorded
  for all `sp`/`spd`, `L`, `M`, and both circular blocks.
- [x] Existing LR-05 point-grid vectors reconstructed with the LR-04 metric;
  maximum relative residual `1.6517e-11 < 1e-10`.
- [x] Energy-affine endpoint oracle passed at maximum relative residual
  `2.2191e-16`.
- [x] LR-GF-02 four-component radial factors matched exactly at printed
  precision; private-helper duplication seam documented.
- [x] Fe `spd` complete point space recorded as 12,375 coordinates versus 232
  unpruned product coordinates (53.34x coordinate reduction).
- [x] TDVK-02 was not changed to PASS; TDVK-03 was not started; no production
  product-basis LR-06 path was added.

The detailed evidence is in
[`docs/LR_LMTO_PRODUCT_RESPONSE_BASIS.md`](LR_LMTO_PRODUCT_RESPONSE_BASIS.md).
The next task should be the orchestrator-approved mechanical design of a
complete-space streaming or matrix-free LR-06 material smoke that preserves the
current equations and radial mesh.

## TDVK-02R0b scaled LMTO product-basis conditioning

**Status: PASS for accepted-Fe radial conditioning; live Fe transition-span
spot check remains unavailable.** This is a numerical representation audit,
not material response validation. TDVK-02 remains **BLOCKED** at the dense
complete-space LR-06 response, and TDVK-02R0 remains the separate algebraic
closure audit.

- [x] Direct LAPACK SVD of `W^(1/2) Btilde` used for every one-site Fe
  `L=0..4` block; unscaled Gram eigenvalues were not used as rank decisions.
- [x] Accepted Fe setup/provenance recorded at HEAD
  `c76b096e6ef28de27e16e91acfc0a9437720e7f8`: `nr=495`, canonical moment
  `2.267442 mu_B`, `EF` approximately `-0.085122 Ry`, residual `6.903e-7`.
- [x] Fe column norms, complete singular spectra, condition estimates, and
  ranks at `tau1/tau10/tau100` recorded for both circular channels.
- [x] Fe ranks are stable: `12/12/12`, `16/16/16`, `16/16/16`, `8/8/8`, and
  `4/4/4` for `L=0..4`; retained `Nprod=232` at all three thresholds.
- [x] Synthetic TDVK-02R0 `spd` fixture rerun; original raw-Gram evidence was
  preserved. Its one-rank `L=1,2` sensitivity shifts are fixture-only.
- [ ] Live Gamma Fe nearest, deeper, and non-negligible-norm transition-span
  residuals at all thresholds; reciprocal eigenvectors are not persisted in
  the accepted radial snapshot artifact.
- [x] No compressed LR-06, LR-GF-02, KXC, Dyson, driver, Fe-physics, eta,
  response-cutoff, or radial-mesh change was committed.

Detailed spectra and provenance are in
[`docs/LR_LMTO_PRODUCT_RESPONSE_BASIS.md`](LR_LMTO_PRODUCT_RESPONSE_BASIS.md).
The next task should be the orchestrator-approved mechanical design of a
complete-space streaming or matrix-free LR-06 material smoke, including a
minimal persisted Gamma eigenpair artifact to finish the R0b live span check.

## TDVK-02R1 production LMTO product basis and transition vertex

**Status: PASS for fixture-level representation/vertex closure; no change to
the TDVK-02 BLOCKED material verdict.** The reusable scaled-SVD product
representation and analytical `z=F t` transition map are implemented without
modifying LR-05, LR-06 accumulation, LR-GF-02, KXC, Dyson, Goldstone logic, or
the material driver.

- [x] `sp=52` and `spd=232` unpruned inventories retained with deterministic
  product-coordinate indexing.
- [x] Direct weighted SVD stores `D`, `Sigma`, retained `U`, `V^H`, and the
  forward map `F=Sigma V^H D`; no inverse singular values are used at runtime.
- [x] Candidate-space and orthonormal-coordinate oracles agree with unchanged
  LR-05 below `1e-10` for four deterministic eigenvector pairs in both
  circular channels; maxima are `3.1313e-16` and `7.4467e-16` respectively.
- [x] Metric norm closure and retained-SVD reconstruction pass; maxima are
  `1.1448e-15` and `3.6871e-15`.
- [x] Production rank sensitivity fails closed; the existing pathological
  fixture is only exercised in explicitly non-strict diagnostic mode.
- [x] Live Fe Gamma transition oracle completed in TDVK-02R2 at the natural
  reciprocal-state handoff; the maximum residual was `1.5853e-15`.
- [x] No compressed LR-06 or susceptibility path was added.

Detailed API, spectra, oracle residuals, and the R1 checklist are in
[`docs/LR_LMTO_PRODUCT_RESPONSE_BASIS.md`](LR_LMTO_PRODUCT_RESPONSE_BASIS.md).
TDVK-02R2 below completes the deferred live transition-span validation before
compact Lehmann accumulation.

## TDVK-02R2 compact LMTO-product Lehmann response

**Status: PASS for compact Lehmann bare-response execution; the overall
TDVK-02 verdict remains BLOCKED for the full reciprocal TD-DFT lifecycle.**

This slice adds a separate product-space LR-06 evaluator using the certified R1
transition coordinates. It preserves the legacy point-grid LR-06 evaluator and
uses the weighted-orthonormal representation
`chi_P = U^H W^(1/2) chi_raw W^(1/2) U`. Product coordinates therefore use the
ordinary Euclidean adjoint and no extra radial weight or LR-04 right-weighted
`B=chi_raw*W` naming is introduced.

### Evidence

- [x] Independent fixture projection: both channels, Gamma and finite q, static
  and finite omega, and two eta values agree with legacy LR-06 at a maximum
  relative residual of `1.1599e-15`.
- [x] Independent pair/denominator accumulation agrees at `0.0000e+00`; the
  occupation skip count and immutable electronic snapshots also pass.
- [x] Live accepted-Fe Gamma LR-05-to-R1 transition oracle passes for the
  nearest/deeper/additional transitions in both channels, with maximum
  residual `1.5852725e-15`.
- [x] Accepted bcc-Fe compact Gamma smoke evaluates the complete 232-coordinate
  product space and is finite for `omega=0, 0.02, 0.05 Ry` at `eta=0.02 Ry`.
- [x] Compact storage for the three-frequency 232x232 result is `2.463867 MiB`
  and measured compact accumulation CPU time is `176.099 s`.
- [x] No point-space transition allocation occurs in the compact evaluator; the
  point-grid vector is used only by the live validation oracle.
- [x] KXC, GSR, Dyson/loss, Goldstone, LR-GF-02, Fe physics, eta, and TDVK-03
  remain untouched and the compact smoke stops before KXC/Dyson.

The accepted handoff reproduced the documented Fe values: SCF residual
`6.903e-7`, moment `2.2674448 mu_B`, and `EF=-0.0851220945 Ry`. The compact
smoke is a validation-only branch selected by `backend='product_lehmann'` at
the naturally prepared reciprocal handoff; it is not a full TDDFT PASS route.

The focused test is registered as `UnitLrProductKsSusceptibility` and is run
with:

```text
ctest --test-dir build --output-on-failure -R '^(UnitLrLmtoProductResponseBasis|UnitLrLmtoProductResponse|UnitLrLmtoProductStrictRankGuard|UnitLrProductKsSusceptibility)$'
```

All four retained R0/R0b/R1 product tests and the R2 fixture test passed. The
next task is compact reciprocal-GF equivalence, followed by compact kernel and
Dyson representation integration under orchestrator approval.

## TDVK-02R3 — compact reciprocal-GF susceptibility

**Status: PASS for compact representation equivalence; not a full TDVK-02
lifecycle PASS.** The implementation is confined to a separate product-space
GF evaluator and a product-basis component-vertex helper. KXC, GSR, Dyson,
Goldstone, eta conventions, Fe physics, and the native real-space GF route
were not changed; TDVK-03 was not started.

The new route uses the existing LR-GF-02 real-axis resolvents and Kubo bubble:
retarded/advanced Green functions, `A=i(GR-GA)/(2*pi)`, moments `p=0..2`,
both Kubo terms, Simpson quadrature, the existing Fermi/k-point prefactor,
the exact q endpoint, and automatic `integration_eta=eta/40` with the strict
`integration_eta < eta` guard. Its compact vertices are built directly from
the R1 candidate descriptors and `forward_transform=Sigma V^H D`, with
component ordering `1+p+2*q`, no radial or Pauli factor, and only the
certified circular spin block. The result is written directly in
`chi_P=U^H W^(1/2) chi_raw W^(1/2) U`; no point-space response matrix is
constructed.

Representation evidence:

- component/R1 maximum residual: `6.3396e-16`;
- maximum independently projected point-GF residual: `1.4205e-15`;
- both circular channels, Gamma and finite q, static and finite omega, and
  21/41-point odd Simpson cases pass;
- focused fixture product dimension: `27` per channel; accepted Fe strict
  product dimension: `232`;
- focused fixture memory: component vertices `110592` bytes, GF matrices
  `18432` bytes, and one compact susceptibility `23328` bytes.

The compact reciprocal-GF route retains the independent real-axis
Green-function/Kubo construction and does not call the Lehmann susceptibility
accumulator.

The compact Lehmann comparison remains diagnostic only. The reported
`(d_F,r_F,d_inf)` values were `(4.2963e3,3.9547,3.5597e3)` for `chi_plus`
Gamma/21, `(4.2434e3,4.6331,3.4043e3)` for `chi_plus` Gamma/41, and
`(1.8018e3,2.2550,7.7385e2)` for finite-q `chi_minus`/21. These values are
not a physics acceptance threshold and were not used to tune GF controls.

The guarded accepted-Fe `product_gf` smoke was configured for Gamma,
`omega=0`, `eta=0.02 Ry`, automatic `eta/40`, and 21 integration points.
It stopped during the 50-step SCF preparation at `diff=2.11269464e-2` and
therefore did not reach the accepted handoff or produce Fe GF timing/result
data. This is recorded as an execution/performance blocker, not as a failure
of the representation oracle. No 2001-point run was made.

Focused verification (all 7 passed):

```text
ctest --test-dir build --output-on-failure -R '^(UnitLrLmtoProductResponseBasis|UnitLrLmtoProductResponse|UnitLrLmtoProductStrictRankGuard|UnitLrProductKsSusceptibility|UnitLrGfSusceptibility|UnitLrProductGfSusceptibility|UnitLrProductGfSusceptibilityRejectIntegrationEta)$'
```

## TDVK-03A — reciprocal-GF real-axis quadrature closure audit

**Status: Case A — `REAL-AXIS GF FORMULATION NUMERICALLY CONSISTENT` for the
nontrivial R3 finite-basis fixture only.** The fixed-eta Simpson ladder,
spectral moments, integration-eta ladder, energy-window ladder, timings, and
the unchanged historical LR-GF-02 cross-check are recorded in
[`TDDFT_RECIPROCAL_GF_QUADRATURE_AUDIT.md`](TDDFT_RECIPROCAL_GF_QUADRATURE_AUDIT.md).

The compact GF response approaches the independent compact Lehmann response
as `integration_eta` decreases at fixed physical `eta=0.04 Ry`; resolved
mesh and window changes are smaller than that finite-width envelope. The
aspirational `rF<1e-5` was not reached at the smallest prescribed finite
integration width, so no extrapolation or tuning was applied. The audit did
not modify GF physics, KXC, GSR, Dyson, Goldstone logic, product vertices, or
the response broadening, and did not run Fe. Timing records:

```text
GF BUBBLE PERFORMANCE REMEDIATION REQUIRED BEFORE FE
```

TDVK-02 and TDVK-03 material status remain unchanged.

## TDVK-03B mixed-eigenvector reciprocal-GF closure

**Status:** PASS for the deterministic complex mixed-eigenvector numerical
oracle; GF/Lehmann convergence, off-diagonal spectral moments, finite-q
`chi_minus`, and the `<1e-10` basis-rotation invariant are recorded in
[`TDDFT_RECIPROCAL_GF_QUADRATURE_AUDIT.md`](TDDFT_RECIPROCAL_GF_QUADRATURE_AUDIT.md).
Production response equations, material inputs, and TDVK-02/TDVK-03 material
status are unchanged.
