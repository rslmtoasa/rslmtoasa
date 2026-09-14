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

## TDVK-03 Fe reciprocal-backend closure

**Status: `TDVK-03 PASS CANDIDATE`.** TDVK-03R reordered the same real-axis
GF/Kubo calculation so that the compact response contraction occurs after the
energy quadrature. The mixed-complex regression, transition factorization
oracle, and three-way scalar/optimized/factorized complete-matrix oracle pass.
The accepted Fe state now completes the prescribed GF ladder in practical
single-sample runtimes. This is a backend-closure result only; no physical
interpretation is claimed.

### Implementation-oracle evidence

The scalar and `38102e1` optimized dense-resolvent contractions remain
available as correctness oracles. The new default `factorized` backend keeps
the existing four component vertices, transforms them to the left/right
eigenbases once per k point, and forms
`T(I,n,m)=sum_pq epsilon_L(n)^p epsilon_R(m)^q Vtilde(I,p,q,n,m)`.
It then performs the explicit Simpson real-axis integral only for the scalar
band-pair kernel
`K_nm=sum_E w_E f(E) 2 w_k/sum(w_k) [a_Ln g_Rm^R + a_Rm g_Ln^A]`, retaining
both Kubo terms, finite `integration_eta`, physical `eta`, occupations, k
weights, signs, conjugations, and the exact endpoint convention. Finally it
forms `chi(I,J)=sum_nm K_nm T(I,n,m) T(J,n,m)*` with a matrix product.
The factorized path does not call the Lehmann accumulator, substitute a
Lehmann occupation-difference denominator, truncate the 232-dimensional
space, or construct the point-space response matrix.

### Transition and three-way oracle evidence

On the mixed-complex fixture, the GF component-vertex-to-eigenbasis transition
amplitudes matched the certified compact transition coordinates as follows:

| channel | maximum absolute error | maximum relative error |
|---|---:|---:|
| `chi_plus` | `1.73046935e-15` | `2.13050107e-16` |
| `chi_minus` | `1.38624879e-16` | `1.53171174e-16` |

The independent mixed-complex Gamma/static and finite-frequency checks and the
mixed-complex finite-q check compare the complete compact matrices from all
three GF backends. In each row below, pairwise columns are ordered
`scalar/optimized`, `scalar/factorized`, and `optimized/factorized`; timings
are ordered `scalar`, `optimized`, `factorized`.

| fixture | omega (Ry) | norm (all three) | pairwise `rF` | pairwise `dInf` | wall time (s) | speedup (`scalar/optimized`, `scalar/factorized`, `optimized/factorized`) |
|---|---:|---:|---|---|---|---|
| mixed Gamma | `0.00` and `0.17` | `1.79506145e+03` | `6.17e-18`, `2.80e-15`, `2.80e-15` | `7.32e-15`, `3.87e-12`, `3.87e-12` | `30.2637`, `4.22790`, `0.00296218` | `7.158`, `1.0217e4`, `1.4273e3` |
| mixed finite q=`(0.23,0,0)` | `0.00` | `1.86403113e+03` | `3.74e-17`, `6.59e-15`, `6.58e-15` | `6.28e-14`, `9.10e-12`, `9.10e-12` | `61.6816`, `8.96611`, `0.00593637` | `6.879`, `1.0390e4`, `1.5104e3` |

The compact unit fixture independently reports factorized relative
differences of `9.28695e-16` (Gamma/static case) and `4.97532e-16`
(finite-q case) against scalar, with maximum transition residual
`1.2064e-16`. Correctness was checked before using these timings as a
benchmark.

On the existing nontrivial compact fixture, with identical inputs and the
complete compact result matrix compared before timing:

| quantity | value |
|---|---:|
| `||chi_scalar||_F` | `1.08638325e+03` |
| `||chi_opt||_F` | `1.08638325e+03` |
| `dF` | `1.46752295e-14` |
| relative Frobenius difference | `1.35083356e-17` |
| `dInf` | `1.42177919e-14` |
| scalar wall time | `7.89120799e-01 s` |
| optimized wall time | `1.14153735e-01 s` |
| speedup | `6.91279001x` |

The already-landed TDVK-03B mixed-complex-eigensystem regression also passes;
its independent LAPACK-generated complex Hermitian fixture, off-diagonal
resolvents, finite-q channel, and basis-rotation invariant are recorded in
[`TDDFT_RECIPROCAL_GF_QUADRATURE_AUDIT.md`](TDDFT_RECIPROCAL_GF_QUADRATURE_AUDIT.md).

### Accepted Fe state and closure harness

The committed audit input is
[`input_tdvk03_fe.nml`](../tests/integration/tddft_driver_smoke/input_tdvk03_fe.nml).
It uses the tracked bcc-Fe `spd` ground-state input, `response_lmax=-1`,
Gamma, `chi_plus`, `omega=0`, physical `eta=0.04 Ry`, complete product space,
and the existing GF controls (`gf_integration_points`,
`gf_integration_eta`, and `gf_energy_margin`). The CTest entry is disabled
because this is a multi-minute-to-hour material audit, not a default unit
test.

The isolated fresh-state run reached the accepted handoff once:

| quantity | accepted value |
|---|---:|
| SCF residual | `1.912e-7` |
| `nbasis / nbands` | `18 / 18` |
| k mesh | `8x8x8 = 512` points |
| Fermi level | `-8.51204947e-02 Ry` |
| temperature | `300 K` |
| spin moment | `2.267442 mu_B` |
| compact product dimension | `232` |

The driver forms the product basis, compact Lehmann result, and every GF
sample after that handoff. All samples point to the same left state, exact
Gamma endpoint, occupations, eigenvectors/eigenvalues, Fermi level,
temperature, radial basis, product basis, channel, and physical response
eta. It does not reconverge SCF between controls.

The planned material ladder is explicit in the driver: width values
`0.010`, `0.005`, `0.0025`, and `0.001 Ry`; Simpson grids resolved to
`h/integration_eta <= 0.4`; a base/fine Simpson comparison; and margins
`0.60`, `1.00`, and `2.00 Ry`. Each report includes complete compact-matrix
`dF`, relative Frobenius difference, `dInf`, both norms, physical eta,
integration eta, bounds, margin, point count, spacing, spacing ratio, and
wall time.

The accepted Fe state used by the completed run was:

| quantity | accepted value |
|---|---:|
| SCF residual | `6.617e-7` |
| `nbasis / nbands` | `18 / 18` |
| k mesh | `8x8x8 = 512` |
| Fermi level | `-8.51191027e-02 Ry` |
| temperature | `300 K` |
| spin moment | `2.267445 mu_B` |
| compact product dimension | `232` |

The product basis, Lehmann response, and all GF samples reused this same
accepted state, including eigenpairs, occupations, Fermi level, temperature,
radial/product basis, channel, exact Gamma endpoint, and physical
`eta=0.04 Ry`. The GF evaluator reported no point-response allocation. For
the Fe state its major factorized allocations were approximately 4.81 MB for
the component vertices, 1.20 MB for the per-k transition amplitudes, and
0.86 MB for the compact response; the dense GF matrix allocation was zero.

The completed Gamma `chi_plus`, `omega=0`, physical `eta=0.04 Ry` ladder was:

| sample | `integration_eta` (Ry) | `N` | margin (Ry) | `h/integration_eta` | `||chi_L||_F` | `||chi_GF||_F` | `dF` | `rF` | `dInf` | wall (s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| eta ladder | `0.0100` | `953` | `0.60` | `0.399426` | `3.79995096` | `3.61768630` | `2.57226666e-01` | `6.76921014e-02` | `1.40114603e-01` | `7.67390` |
| eta ladder | `0.0050` | `1903` | `0.60` | `0.399846` | `3.79995096` | `3.70701147` | `1.28850597e-01` | `3.39084894e-02` | `7.03065297e-02` | `9.34962` |
| eta ladder | `0.0025` | `3805` | `0.60` | `0.399846` | `3.79995096` | `3.75267245` | `6.48187541e-02` | `1.70577870e-02` | `3.53119546e-02` | `12.7772` |
| eta ladder / Simpson base / window 0.60 | `0.0010` | `9509` | `0.60` | `0.399931` | `3.79995096` | `3.78102764` | `2.58386514e-02` | `6.79973285e-03` | `1.40767444e-02` | `23.0595` / `23.1370` / `23.0656` |
| Simpson fine | `0.0010` | `19017` | `0.60` | `0.199965` | `3.79995096` | `3.78098824` | `2.58942459e-02` | `6.81436317e-03` | `1.40921415e-02` | `40.3156` |
| window 1.00 | `0.0010` | `11509` | `1.00` | `0.399943` | `3.79995096` | `3.78117876` | `2.57510806e-02` | `6.77668760e-03` | `1.40222232e-02` | `26.7371` |
| window 2.00 | `0.0010` | `16509` | `2.00` | `0.399960` | `3.79995096` | `3.78130496` | `2.56779024e-02` | `6.75742994e-03` | `1.39758163e-02` | `35.8025` |

The eta ladder is monotonic toward the compact Lehmann result. The base/fine
Simpson change is small relative to the remaining finite-width envelope, and
the window changes are similarly subdominant; no window or tolerance tuning
was applied. The first controlled Fe sample, `N=953`, completed in `7.67390 s`.

The TDVK-03 completion checklist is now:

- [x] mixed-complex reciprocal-GF regression;
- [x] GF component-vertex to eigenbasis transition-factorization oracle;
- [x] scalar-versus-optimized complete-matrix oracle;
- [x] scalar/optimized/factorized complete-matrix three-way oracle at Gamma,
  finite frequency, mixed-complex, and finite q;
- [x] factorized backend retains explicit Simpson real-axis quadrature, finite
  `integration_eta`, physical `eta`, both Kubo terms, and backend independence;
- [x] independent optimized real-axis GF/Kubo route;
- [x] one accepted Fe state reused by Lehmann and all planned GF controls;
- [x] complete 232-coordinate Gamma `chi_plus`, static closure harness;
- [x] separate Simpson, integration-width, and energy-window controls;
- [x] complete matrix diagnostics and provenance output;
- [x] practical Fe GF ladder and material backend closure after quadrature factorization.

No final GF↔Lehmann material tolerance was invented. No KXC, Dyson,
Goldstone, or physical interpretation was added; the `PASS CANDIDATE` status
is returned to the orchestrator for the parent milestone.

## TDVK-04 Fe finite-q endpoint, folding, and covariance

**Status: PASS CANDIDATE.** The prescribed finite-q compact bare-response
validation completed on the accepted bcc-Fe reciprocal state. This is endpoint,
folding, covariance, and representative backend evidence only; no magnon,
dispersion, stiffness, damping, or literature claim is made.

The implementation is the validation-only `backend='product_finite_q'` branch
in [`tddft_production_driver.f90`](../source/tddft_production_driver.f90). It
allocates the strict 232-coordinate weighted-orthonormal product representation,
runs compact Lehmann at every q first, and then runs two factorized compact GF
spot checks. It does not allocate the legacy point-response matrix and does not
enter KXC, Goldstone, Dyson, loss, or mode fitting. The reproducible input is
[`input_tdvk04_fe.nml`](../tests/integration/tddft_driver_smoke/input_tdvk04_fe.nml);
the CTest registration is intentionally disabled because the material run takes
several minutes and is run manually for this evidence.

### Accepted state and controls

| quantity | value |
|---|---:|
| SCF control / final report RMS difference | `conv_thr=1e-6` / `~1e-6` |
| `nbasis / nbands` | `18 / 18` |
| k mesh | `8x8x8 = 512` points |
| Fermi level | `-8.5120871761e-02 Ry` |
| temperature | `300 K` |
| accepted LR-01 integrated moment | `2.26745419 mu_B` |
| compact product dimension | `232` |
| response angular cutoff | `4` (complete `spd` product space) |
| channel / physical eta | `chi_plus` / `0.04 Ry` |
| GF controls | `N=6401`, `integration_eta=0.001 Ry`, margin `0.60 Ry` |

The state provenance is one immutable reciprocal `ham_only`, second-order,
orthonormal, collinear, no-SOC eigensystem handoff. The same eigenpairs,
occupations, Fermi level, temperature, direct radial mesh, compact basis,
channel, and physical eta were used for every row.

### Literal q set and endpoint metadata

The q values are direct reciprocal coordinates, reported without symmetry-line
labels. Every endpoint has 512 unique folded points and maximum exact-folding
error zero at the recorded precision.

| index | supplied q | folded q | first folded endpoint | last folded endpoint | unique | max error |
|---:|---|---|---|---|---:|---:|
| 1 | `(0, 0, 0)` | `(0, 0, 0)` | `(-0.4375,-0.4375,-0.4375)` | `(0.4375,0.4375,0.4375)` | 512 | `0.0` |
| 2 | `(0.125, 0, 0)` | `(0.125, 0, 0)` | `(-0.3125,-0.4375,-0.4375)` | `(-0.4375,0.4375,0.4375)` | 512 | `0.0` |
| 3 | `(-0.125, 0, 0)` | `(-0.125, 0, 0)` | `(0.4375,-0.4375,-0.4375)` | `(0.3125,0.4375,0.4375)` | 512 | `0.0` |
| 4 | `(0.23,0.07,-0.11)` | `(0.23,0.07,-0.11)` | `(-0.2075,-0.3675,0.4525)` | `(-0.3325,-0.4925,0.3275)` | 512 | `0.0` |

### Compact Lehmann diagnostics

All matrices were finite. The trace is the ordinary compact-coordinate diagonal
trace; it is retained as a diagnostic and is not interpreted as a mode.

| q index | supplied q | `||chi||_F` | max element | trace real | trace imag | transitions |
|---:|---|---:|---:|---:|---:|---:|
| 1 | `(0,0,0)` | `3.79995937` | `1.87067728` | `-12.87280517` | `-1.94800930` | 104768 |
| 2 | `(0.125,0,0)` | `3.72587755` | `1.82476632` | `-12.73979612` | `-1.96270817` | 106697 |
| 3 | `(-0.125,0,0)` | `3.72587755` | `1.82476632` | `-12.73979612` | `-1.96270817` | 106697 |
| 4 | `(0.23,0.07,-0.11)` | `3.60351169` | `1.73660268` | `-12.55186727` | `-2.03122412` | 106521 |

### q↔-q covariance

At `omega=0`, the full compact matrices at q index 2 and q index 3 satisfy
the established LR-03/LR-06 retarded/advanced circular covariance with maximum
absolute residual `1.44868418e-15`. In compact coordinates the existing identity
is represented by the `M→-M`, `(-1)^M` angular map and the orthonormal radial
mode transport between the plus and minus product bases; no new response
formula, channel swap, Fourier sign, or hand-inserted phase was introduced.

| q | -q | channel pair | omega (Ry) | max full-matrix residual |
|---|---|---|---:|---:|
| `(0.125,0,0)` | `(-0.125,0,0)` | `chi_plus` / `chi_minus` | `0.0` | `1.44868418e-15` |

### Representative finite-q compact GF checks

The factorized reciprocal-GF route used the same accepted state, product basis,
physical eta, endpoint snapshots, and q values as Lehmann. These are the two
post-Lehmann representative checks required by TDVK-04, not a new GF
convergence ladder. The differences are reported without inventing a material
acceptance threshold.

| q index | q | `h/integration_eta` | `||chi_L||_F` | `||chi_GF||_F` | `d_F` | `r_F` | `d_inf` | wall (s) |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 2 | `(0.125,0,0)` | `0.59414714` | `3.72587755` | `3.70930128` | `2.27443844e-2` | `6.10443689e-3` | `1.22972669e-2` | `17.3901` |
| 4 | `(0.23,0.07,-0.11)` | `0.60050009` | `3.60351169` | `3.58766009` | `2.54495298e-2` | `7.06242465e-3` | `1.25762706e-2` | `17.3541` |

TDVK-04 therefore returns **PASS CANDIDATE** to the orchestrator for the
finite-q bare-response gate. Downstream dispersion, magnon, stiffness,
damping, and literature work remains outside this milestone.
