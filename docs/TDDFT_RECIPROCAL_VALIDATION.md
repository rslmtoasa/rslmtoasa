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

## TDVK-05 Fe numerical convergence

TDVK-05 was run against live build
`pre_tddft_cleanroom_20260909-35-g2a84-dirty` using the validation-only
`backend='product_convergence'`. Each isolated process accepts the converged
real-space Fe state and then regenerates the requested complete reciprocal
mesh. The convergence comparison treats the accepted real-space potential,
radial data, moment, and EF as fixed handoff inputs; only reciprocal sampling
and quantities derived from it are convergence evidence. Every campaign log
reported `Converged!`, every serialized response row was finite, and no KXC,
Goldstone, Dyson, loss, or mode-fitting route was entered.

### Accepted-state and mesh provenance

| case | generated mesh / nk | response cutoff | EF (Ry) | moment (μB) | SCF residual control | product unpruned / retained | all-L rank stable |
|---|---:|---:|---:|---:|---:|---:|:---:|
| `mesh4_full` | 4³ / 64 | complete `L=4` | -0.0851199488840 | 2.267462473383 | 4.1283e-7 | 232 / 232 | T |
| `mesh8_full_eta_ladder` | 8³ / 512 | complete `L=4` | -0.0851199488840 | 2.267462473383 | 4.1283e-7 | 232 / 232 | T |
| `mesh12_full` | 12³ / 1728 | complete `L=4` | -0.0851199488840 | 2.267462473383 | 4.1283e-7 | 232 / 232 | T |
| `mesh8_reduced_lmax2` | 8³ / 512 | reduced `L=2` (approximate) | -0.0851199488840 | 2.267462473383 | 4.1283e-7 | 140 / 140 | T |

All cases used `nbasis=18`, `nbands=18`, `temperature=300 K`,
`reciprocal_mode=ham_only`, and `hamiltonian_order=second`. Within each
process, the response rows reuse the accepted eigenpairs, occupations, EF,
temperature, and weights. The equal printed EF/moment values are fixed
accepted-state provenance, not cross-mesh convergence evidence. The mesh
comparison is therefore more precisely described as independently regenerated
reciprocal samplings of an accepted real-space state.

### Fixed versus mesh-dependent quantities

The live handoff is explicit: after SCF, the driver sets
`reciprocal_obj%fermi_level = energy_obj%fermi`, sets
`auto_find_fermi = .false.`, builds and diagonalizes the requested mesh, and
then calls `calculate_canonical_band_energy(.false.)`. The immutable LR
snapshot regenerates explicit Fermi occupations from those eigenvalues, the
accepted EF, and 300 K. No reciprocal response path solves EF or rebuilds the
real-space potential/radial functions.

| quantity | provenance classification | live source/path |
|---|---|---|
| potential | **FIXED ACCEPTED-STATE INPUT** | converged real-space SCF objects handed to the production driver |
| radial functions | **FIXED ACCEPTED-STATE INPUT** | accepted `radial_ground_state`/LR radial basis; product SVD uses its direct radial mesh |
| magnetic moment | **FIXED ACCEPTED-STATE INPUT** | accepted SCF ground-state integrated moment |
| EF | **FIXED ACCEPTED-STATE INPUT** | `energy_obj%fermi` copied to `reciprocal_obj%fermi_level`; `auto_find_fermi=.false.` |
| reciprocal eigenvalues | **RECOMPUTED PER RESPONSE K MESH** | `generate_mp_mesh` → k-space Hamiltonian → `diagonalize_hamiltonian` |
| reciprocal eigenvectors | **RECOMPUTED PER RESPONSE K MESH** | same requested-mesh diagonalization, copied into the LR snapshot |
| k vectors | **RECOMPUTED PER RESPONSE K MESH** | complete replicated `k_workset` generated from `nk_mesh` |
| k weights | **RECOMPUTED PER RESPONSE K MESH** | complete-BZ workset weights, normalized by the canonical occupation service |
| occupations | **DERIVED PER RESPONSE K MESH** | explicit Fermi function of each mesh eigenvalue at fixed EF and 300 K |
| electron-count/BZ occupation integral | **DERIVED PER RESPONSE K MESH** | `evaluate_eigenvalue_occupations` over that mesh's eigenvalues/weights |
| compact `chiKS` | **DERIVED PER RESPONSE K MESH** | compact Lehmann contraction from the immutable mesh-specific LR states |

The optional mesh-specific EF values below are diagnostics from the existing
reciprocal Fermi solver only; they are never fed back into the response.

### k-mesh dependence — complete product span, physical eta = 0.01 Ry

The table reports the required Gamma-static, certified finite-q-static, and
low-finite-omega (`Gamma`, 0.02 Ry) points. `trace` is shown as real + i
imaginary; runtime is the compact Lehmann evaluator CPU time for the two
frequencies in that q request.

| mesh | q | omega (Ry) | Frobenius norm | max element | trace | runtime (s) |
|---:|---|---:|---:|---:|---:|---:|
| 4³ | Gamma | 0.00 | 3.930186 | 1.920684 | -13.716734 - 0.553225i | 15.4664 |
| 4³ | (0.125,0,0) | 0.00 | 4.166316 | 2.043431 | -14.620085 - 0.781799i | 15.9747 |
| 4³ | Gamma | 0.02 | 4.361765 | 2.153108 | -15.032238 - 0.814079i | 15.4664 |
| 8³ | Gamma | 0.00 | 3.935090 | 1.916506 | -13.623680 - 0.461280i | 122.9355 |
| 8³ | (0.125,0,0) | 0.00 | 3.883498 | 1.896669 | -13.516282 - 0.472147i | 124.8851 |
| 8³ | Gamma | 0.02 | 4.360642 | 2.153881 | -14.818109 - 0.734805i | 122.9355 |
| 12³ | Gamma | 0.00 | 3.944379 | 1.928807 | -13.573441 - 0.525218i | 410.1865 |
| 12³ | (0.125,0,0) | 0.00 | 3.953125 | 1.930425 | -13.776751 - 0.593291i | 416.0425 |
| 12³ | Gamma | 0.02 | 4.379360 | 2.168057 | -14.838767 - 0.815430i | 410.1865 |

The complete-span rows evaluated 104770/61118 transitions/skips on 8³ and
347612/212260 transitions/skips on 12³ for Gamma (the finite-q counts were
106697/59191 and 355508/204364, respectively). No non-finite value occurred.
The coarse 4³ result is separated from the 8³/12³ values, while the latter
show a bounded mesh trend for this bare-response diagnostic; no material
acceptance threshold is asserted here.

### Physical response-eta dependence — complete product span, 8³

The q and frequency points are identical for all three physical eta values.
These rows do not use or vary the GF integration controls.

| eta (Ry) | q | omega (Ry) | Frobenius norm | max element | trace |
|---:|---|---:|---:|---:|---:|
| 0.020 | Gamma | 0.00 | 3.902158 | 1.907050 | -13.432158 - 0.957331i |
| 0.010 | Gamma | 0.00 | 3.935090 | 1.916506 | -13.623680 - 0.461280i |
| 0.005 | Gamma | 0.00 | 3.944643 | 1.918894 | -13.681335 - 0.226607i |
| 0.020 | (0.125,0,0) | 0.00 | 3.841013 | 1.876618 | -13.311860 - 0.980295i |
| 0.010 | (0.125,0,0) | 0.00 | 3.883498 | 1.896669 | -13.516282 - 0.472147i |
| 0.005 | (0.125,0,0) | 0.00 | 3.898781 | 1.904782 | -13.585594 - 0.230868i |

The corresponding low-finite-omega rows were also evaluated at each eta;
their Gamma norms were 4.322279, 4.360642, and 4.370845 for eta 0.020,
0.010, and 0.005 Ry, respectively. The observed decrease in the imaginary
trace with decreasing physical eta is recorded as a numerical trend only;
there is no eta-to-zero extrapolation.

### Full versus reduced response cutoff — 8³, eta = 0.01 Ry

The complete `response_lmax=-1` request resolves to `L=4` and retains the
232-dimensional certified `spd` product span. The `response_lmax=2` request
retains a 140-dimensional span and is an explicitly approximate diagnostic,
not a production replacement.

| q | omega (Ry) | full `L=4` norm / max | reduced `L=2` norm / max | full trace | reduced trace |
|---|---:|---:|---:|---:|---:|
| Gamma | 0.00 | 3.935090 / 1.916506 | 3.273866 / 1.916506 | -13.623680 - 0.461280i | -7.156546 - 0.273685i |
| (0.125,0,0) | 0.00 | 3.883498 / 1.896669 | 3.224647 / 1.896669 | -13.516282 - 0.472147i | -7.042038 - 0.264890i |
| Gamma | 0.02 | 4.360642 / 2.153881 | 3.646223 / 2.153881 | -14.818109 - 0.734805i | -7.836362 - 0.412388i |

Both product constructions reported all-L rank stability. The reduced result
is retained solely to expose cutoff sensitivity and is not substituted for the
complete span.

### Direct radial and product-basis provenance

| item | recorded value |
|---|---|
| radial identity | `nr=495`, `a=0.0200000000000`, `b=1.36280409302648e-4`, `rmax=2.66219999999995`, `sum(r)=134.378214448223` |
| radial provenance | accepted direct LR-01 logarithmic mesh; no decimation or material-dependent radial truncation |
| XC provenance | Barth-Hedin / legacy RS-LMTO |
| complete product | `spd`, `L=4`, unpruned 232, retained 232, all-L rank stable `T` |
| reduced diagnostic | `L=2`, unpruned 140, retained 140, all-L rank stable `T`, approximate only |

No site-only response, radial decimation, or product-mode truncation was used.

### Representative reciprocal-GF spot checks

The selected 8³ accepted state used physical response `eta=0.01 Ry`,
`integration_points=6401`, and `integration_eta=0.001 Ry`. These are the two
TDVK-04-certified representative spots, not a GF convergence ladder.

| q index | q | h / integration eta | `||chi_L||_F` | `||chi_GF||_F` | dF | rF | dInf | wall (s) |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 2 | (0.125,0,0) | 0.59414678 | 3.883498 | 3.865953 | 2.519235e-2 | 6.487024e-3 | 1.353574e-2 | 17.6671 |
| 4 | (0.23,0.07,-0.11) | 0.60049971 | 3.815074 | 3.799173 | 2.968604e-2 | 7.781250e-3 | 1.291837e-2 | 18.0983 |

Both GF spot responses were finite. The compact Lehmann and GF matrices use
the same accepted eigenpairs, occupations, EF, temperature, radial mesh,
product basis, channel, and physical eta.

TDVK-05 therefore returns **PASS CANDIDATE** to the orchestrator. The evidence
supports proceeding only through the mandatory review gate; interaction,
Goldstone, Dyson, dispersion, magnon, stiffness, damping, and literature work
remain outside this milestone.

## TDVK-05 review-gate closure — genuine reciprocal convergence

This closure was executed from local starting HEAD `c9f60a6`
(`tests: map Fe reciprocal TDDFT convergence`). The requested live
`git fetch origin fable_v4` could not authenticate because the environment had
no usable GitHub SSH key (`Permission denied (publickey)`); the local branch
was already at `origin/fable_v4` before the closure edits. Unrelated worktree
changes were preserved. No TDVK-05R1/R2 task or TDVK-06 work was started.

### Runtime mesh provenance and fixed-EF consequence

The closure harness launches a fresh executable process and scratch directory
for every case. It parses the generated `.kpoints`, `.basis`, and `.matrix`
artifacts; requested mesh labels are never substituted for runtime values. The
SHA-256 fingerprint is over the complete serialized runtime k-point/weight
data rows, not an expected mesh generated by the harness.

| case | requested mesh | actual mesh | actual nk | weight sum | first / middle / last representative k | full-list+weight SHA-256 |
|---|---:|---:|---:|---:|---|---|
| `mesh4_full` | `4³` | `4×4×4` | 64 | 1.000000000000000 | `(-.375,-.375,-.375)` / `(.375,.375,-.125)` / `(.375,.375,.375)` | `6d9e66cc3a27885cccefbd9ab5ad973e3fb2fea2e536c0a93bbb8e8aec11a224` |
| `mesh8_full_eta_ladder` | `8³` | `8×8×8` | 512 | 1.000000000000000 | `(-.4375,-.4375,-.4375)` / `(.4375,.4375,-.0625)` / `(.4375,.4375,.4375)` | `1cbe24b64cb808294ebdf221e5227e84d97ae147d9c219c1f7ea93071782d9d5` |
| `mesh12_full` | `12³` | `12×12×12` | 1728 | 1.000000000000019 | `(-.458333,-.458333,-.458333)` / `(.458333,.458333,-.041667)` / `(.458333,.458333,.458333)` | `cbf15c1e15c8e657635005a071a17e5ead5db64637b2d1938b3dc5e2326b347f` |
| `mesh12_full_eta005_corner` | `12³` | `12×12×12` | 1728 | 1.000000000000019 | same as `mesh12_full` | same as `mesh12_full` |
| `mesh8_reduced_lmax2` | `8³` | `8×8×8` | 512 | 1.000000000000000 | same as `mesh8_full_eta_ladder` | same as `mesh8_full_eta_ladder` |

The three distinct full-product meshes have distinct fingerprints. The two
same-mesh controls intentionally have matching fingerprints, while their
response cutoffs/eta values differ. Every full-product case also reported
finite actual eigenvalue summaries, explicit occupation integrals, and finite
response rows:

| mesh | eigenvalue min / max / mean (Ry) | fixed EF (Ry) | `N_e(EF_fixed)` | target | `N_e-target` | diagnostic `EF_mesh-EF_fixed` (Ry) |
|---:|---:|---:|---:|---:|---:|---:|
| `4³` | `-0.639497 / 1.832958 / 0.334832` | `-0.0851199488840` | 8.0588317182 | 8.0 | `+0.0588317182` | `-0.0048321241` |
| `8³` | `-0.701285 / 1.901255 / 0.334832` | `-0.0851199488840` | 8.0250873032 | 8.0 | `+0.0250873032` | `-0.0021707836` |
| `12³` | `-0.713278 / 1.916099 / 0.334832` | `-0.0851199488840` | 8.0410737796 | 8.0 | `+0.0410737796` | `-0.0027549013` |

`N_e` is evaluated by the live reciprocal occupation service from each
mesh's eigenvalues and weights, using the accepted EF and 300 K Fermi
function. The diagnostic EF values use the existing reciprocal solver only;
they were not fed back. These numbers are reported without an invented
acceptance tolerance. The fixed-EF error is finite and visible, including its
non-monotonic 4³→8³→12³ variation, so it is not mislabeled as convergence of
the fixed state.

### Full compact operator: 8³ versus 12³

At Gamma, `omega=0`, physical `eta=0.01 Ry`, and complete `L=4` product
space, both runtime artifacts have dimension 232 and all-L rank stability `T`.
The serialized weighted radial modes have identical keys and zero
max-absolute difference at the printed precision. Direct comparison is
therefore valid; the independently computed overlap unitarity residual is
`6.7785e-15` (maximum over product blocks). The basis contains 27,720
radial-mode records; singular values span `1.7249e-08` to `2.6300033`.

| `||chi_8||_F` | `||chi_12||_F` | full-matrix `dF` | relative Frobenius | `dInf` | trace difference |
|---:|---:|---:|---:|---:|---|
| 3.9350904558 | 3.9443793459 | 0.0465725977 | 0.0118073323 | 0.0123013720 | `-0.0502389601 + 0.0639382156i` (abs `0.0813145038`) |

This is an operator-level comparison of all 232² compact entries, not a
comparison reconstructed from scalar norms or traces.

### Dense-k/small-eta cross-corner

The additional complete-product `12³`, `eta=0.005 Ry`, Gamma-static operator
has norm `3.9535660858`, max element `1.9312253771`, and trace
`-13.6244005649 - 0.2592802812i`. Full-matrix comparisons use the same
directly aligned basis:

| reference | comparison | `dF` | relative Frobenius | `dInf` | trace-difference abs |
|---|---|---:|---:|---:|---:|
| `12³/.005` | `12³/.01` | 0.0977841369 | 0.0247331485 | 0.0530814717 | 0.2707764187 |
| `12³/.005` | `8³/.005` | 0.0473418518 | 0.0119744683 | 0.0123314194 | 0.0656438537 |
| `12³/.005` | `8³/.01` | 0.1011838735 | 0.0255930649 | 0.0547607152 | 0.2020009859 |

The corner is finite and does not expose a new uncontrolled instability in
the reported compact operator diagnostics. No two-dimensional k/eta grid or
material threshold was introduced.

### Controlled TDVK-04 reciprocal-GF spots

The same isolated 8³ physical setup, q points, `eta=0.01 Ry`,
`integration_eta=0.001 Ry`, energy window, and complete product basis were
rerun with 9,601 Simpson points. The resulting `h/integration_eta` values
are approximately 0.4 and satisfy the established `<=0.5` quadrature
criterion.

| q | `N_E` | `||chi_L||_F` | `||chi_GF||_F` | `dF` | `rF` | `dInf` | `h` | `h/integration_eta` | wall (s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `(0.125,0,0)` | 9601 | 3.8834982379 | 3.8655712209 | 0.0262306332 | 0.0067543827 | 0.0142047377 | 0.0003960979 | 0.3960978521 | 26.7473 |
| `(0.23,0.07,-0.11)` | 9601 | 3.8150738956 | 3.8004892411 | 0.0287407267 | 0.0075334653 | 0.0123889886 | 0.0004003331 | 0.4003331403 | 26.9724 |

For comparison, the previous under-resolved 6,401-point run had
`h/integration_eta=0.59414678` and `0.60049971`, with `(dF,rF,dInf)` of
`(0.02519235,0.00648702,0.01353573)` and
`(0.02968604,0.00778125,0.01291837)` respectively. Both runs were finite;
the new table makes the quadrature control explicit without hard-coding an
expected response value.

### Harness audit and closure decision

The closure harness audit passed the following checks: one fresh subprocess
and unique scratch directory per case; actual requested-vs-generated mesh
dimensions and k-list row counts; distinct full-mesh fingerprints; matching
same-mesh control fingerprints; fresh, uniquely named runtime artifacts;
actual eigenvalue/occupation headers; no baseline scalar substitution; and no
expected fingerprint generation. The driver has no cross-process response
cache, and the harness rejects any artifact reuse without matching runtime
provenance. Product responses are populated from parsed runtime rows and full
matrix artifacts.

Subject to orchestrator review of the explicitly reported fixed-EF
occupation errors, this closure satisfies the eight review-gate evidence
requirements and returns:

**TDVK-05 REVIEW-GATE PASS CANDIDATE**

This is the final TDVK-05 gate result. No KXC/GSR or TDVK-06 work was started.

## k-space SCF → TDDFT state-consistency closure

This is the intermediate review-gate bridge between the historical TDVK-05
fixed-potential diagnostics above and TDVK-06. The historical workflow remains
diagnostic evidence:

```text
accepted real-space SCF → fixed accepted EF → post-SCF reciprocal rebuild → TDDFT
```

The production-validation workflow is now:

```text
use_kspace=.true. SCF on Nk → accepted reciprocal cache → same EF/occupations/eigenpairs → TDDFT
```

The live SCF branch constructs the configured Monkhorst-Pack mesh, assembles
and diagonalizes its reciprocal Hamiltonian on every SCF iteration, evaluates
Fermi occupations and EF with the reciprocal electron-number solver, maps the
projected moments into the SCF mixer, updates the atomic potential, and
repeats. The final handoff refreshes the reciprocal Hamiltonian/eigenpairs once
on the final accepted mixed potential, then solves EF again from that final
eigensystem without a further density or mixer update. Thus the accepted
reciprocal cache is the final potential's state, not the pre-mix spectrum from
the last ordinary iteration.

### Handoff implementation and EF ownership

Preprocessing retains the old reciprocal constructor only for the diagnostic
real-space path. For `use_kspace=.true.`, it serializes
`kspace_scf_state.dat` and passes `self%reciprocal_scf_cache` directly to the
TDDFT driver. The direct branch validates a complete replicated mesh, cached
eigenpairs, canonical occupations, Hamiltonian provenance, and EF identity;
it does not call `generate_mp_mesh`, `build_kspace_hamiltonian`, or
`diagonalize_hamiltonian`. The TDDFT left-state snapshot is immutable after
this handoff.

The input `&energy fermi` value is only the initial numerical seed. Production
EF ownership is the reciprocal `auto_find_fermi=.true.` electron-number solver
using the actual mesh weights and the stable Fermi-Dirac function at 300 K.
No old real-space EF was restored and no response susceptibility was used to
tune EF.

The live fetch requested by this review gate could not authenticate in the
execution environment (`git@github.com: Permission denied (publickey)`); the
local `fable_v4` branch started at `55d2f8c`, equal to its local `origin/fable_v4`
ref. Existing unrelated worktree edits were preserved.

### Fresh k-space-SCF production states

The bridge harness ran only fresh Fe `8×8×8` and `12×12×12` cases from the
tracked TDVK-05 input template. Each case records an exact generated
`input_diff.patch`. The two meshes were not forced to share EF or moment.

| state | iterations | EF (Ry) | integrated N | N−target | moment (μB) | SCF residual | physical energy (Ry) | k-point fingerprint |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| 8³ | 17 | -0.08827224260182974 | 8.000000000017071 | 1.7071e-11 | 2.163828769205924 | 2.8333e-7 | -2542.018824874639 | `1cbe24b64cb808294ebdf221e5227e84d97ae147d9c219c1f7ea93071782d9d5` |
| 12³ | 13 | -0.08756873831087462 | 8.000000000021009 | 2.1009e-11 | 2.159350841118554 | 1.8942e-7 | -2542.017860147606 | `cbf15c1e15c8e657635005a071a17e5ead5db64637b2d1938b3dc5e2326b347f` |

Both states passed the unchanged collinear/no-SOC/bulk/orthogonal/`ham_only`
second-order-HOH `sp`/`spd` capability contract, with the accepted direct
LR-01 radial mesh and Barth-Hedin / legacy RS-LMTO XC provenance.

### State-identity tests before susceptibility

The SCF and TDDFT state artifacts contain every actual k vector and weight,
all eigenvalues, explicit occupations, and the occupation-weighted
one-particle density matrix. The density matrix is the gauge-invariant
eigenvector/subspace diagnostic; raw eigenvector columns are not compared by
phase. The harness compares the complete tagged arrays, not just labels or
scalar summaries.

| state | SCF/TDDFT fingerprint | max mesh Δ | max EF Δ (Ry) | max eigenvalue Δ (Ry) | max occupation Δ | projector max Δ | projector Frobenius Δ | TDDFT rebuild |
|---|---|---:|---:|---:|---:|---:|---:|---|
| 8³ | identical, `1cbe24b...782d9d5` | 0 | 0 | 0 | 0 | 0 | 0 | `F` |
| 12³ | identical, `cbf15c...26b347f` | 0 | 0 | 0 | 0 | 0 | 0 | `F` |

The driver fails closed if the accepted state is fixed-EF, incomplete,
distributed, noncanonical, on an incompatible mesh, or inconsistent with the
energy EF. It also records `reciprocal_rebuild_performed_for_tddft = F` for
the production cases.

### Complete bare-response convergence evidence

Each accepted state produced only the complete 232-dimensional compact
Lehmann response at Gamma, `omega=0`, physical `eta=0.01 Ry`. The 12³
`eta=0.005 Ry` corner reused its accepted state and did not rerun SCF.

| state / eta | ‖χ‖F | max element | trace |
|---|---:|---:|---|
| 8³ / .010 | 4.0359342768 | 1.9691868525 | -13.8790422644 - 0.5206267725i |
| 12³ / .010 | 4.0383596359 | 1.9717277263 | -13.8225688390 - 0.5572045507i |
| 12³ / .005 | 4.0487593319 | 1.9744565290 | -13.8804953090 - 0.2752317081i |

The complete 12³ operator was transported into the 8³ weighted product basis
with the existing block-overlap alignment. Both bases contain 232 modes and
27,720 serialized radial-mode records. The maximum overlap unitarity residual
is `1.7575e-1`, and the independently self-consistent radial basis difference
is retained as evidence rather than hidden.

| aligned basis dimension | ‖χ8‖F | ‖χ12→8‖F | dF | relative dF | dInf | |Δtrace| |
|---:|---:|---:|---:|---:|---:|---:|
| 232 | 4.0359342768 | 4.0383596359 | 0.0516202363 | 0.0127824763 | 0.0127379503 | 0.0672843342 |

No material 8³↔12³ acceptance threshold was invented; this operator-level
evidence is returned to the orchestrator.

### One same-state reciprocal-GF spot check

The 8³ process performed one Gamma reciprocal-GF check after the Lehmann
state contract, using 9,601 Simpson points, integration `eta=0.001 Ry`, and
`h/integration_eta=0.3960`. It was finite with:

| dF | rF | dInf |
|---:|---:|---:|
| 0.0306663233 | 0.0075983208 | 0.0160049429 |

The GF and Lehmann services consumed the same accepted EF, occupations,
temperature, eigenpairs, mesh weights, radial mesh, product basis, channel,
and physical response eta in the same process.

### Closure decision

The bridge harness passed its fresh-process, actual-mesh, state-array,
occupation, gauge-invariant projector, direct-handoff, complete-operator,
eta-reuse, and same-state GF audits. The recommended production state for a
future TDVK-06 review is the accepted 12³ k-space-SCF cache, subject to
orchestrator adjudication of the explicitly reported operator difference.

`KSPACE-SCF → TDDFT HANDOFF PASS CANDIDATE`

At the time of this original TDVK-05 report, no TDVK-06 interaction physics
had been started.  The later TDVK-06 evidence update follows below.

## TDVK-06 Fe static ALSDA and independent GSR diagnostics

This section is the later TDVK-06 evidence update; the historical TDVK-05
closure text above remains unchanged.  Preflight was performed on branch
`fable_v4` at exact HEAD
`a44d6db63bca9ee3cd0561e47525d29a1485cf75`, equal to the local
`origin/fable_v4` ref, with the TDVK bridge commit an ancestor.  The worktree
was already dirty in unrelated SCF/RSGF/documentation files; those edits were
preserved.

The run used `input_tdvk06_fe.nml` and the validation-only backend
`static_interactions`:

```text
use_kspace = .true.
nk1,nk2,nk3 = 12,12,12
q = (0,0,0), omega = 0 Ry, channel = chi_plus
physical eta = 0.010 Ry and 0.005 Ry
response_lmax = -1, strict complete compact product dimension = 232
state source = accepted_kspace_scf_cache; TDDFT reciprocal rebuild = F
```

The accepted state was consumed directly by TDDFT.  The state artifact has
the complete 1,728-point mesh, weights, eigenvalues, explicit occupations,
and gauge-invariant occupation-weighted projector.  Its k-point serialization
fingerprint is
`cbf15c1e15c8e657635005a071a17e5ead5db64637b2d1938b3dc5e2326b347f`.

| accepted state quantity | value |
|---|---:|
| EF (Ry) | -0.087568835307212323 |
| integrated electron count | 7.999999999986562 |
| target electron count | 8.000000000000000 |
| accepted LR-01 radial moment (μB) | 2.159356663477833 |
| accepted radial residual control | 4.123523366655037e-7 |
| state mesh/weight/EF/eigen/occupation/projector residuals | 0 |

The accepted XC provenance was Barth-Hedin, legacy RS-LMTO, TXC=1, with
mapping quality `REFERENCE_EQUIVALENT`.  The Pauli magnetization input was
constructed from the accepted reciprocal occupations/eigenvectors and the
accepted POTPAR large-component plus frozen-core projection, and was passed
to KXC with the exact label `pauli_projected`.

### Compact representation certification

The independently tested mapping is

```text
c = U^H sqrt(W) x
x = inv(sqrt(W)) U c
K_compact = U^H K_point U
```

`UnitLrCompactStaticInteraction` passed its explicit point/product projection,
reconstruction, local-operator, action, and normalized `m00` nested-loop
oracles.  The production run independently required strict rank stability for
all retained blocks and observed the complete 232-mode inventory.  This
certifies representation use for the static diagnostic; it does not certify a
compact Dyson denominator.

### Raw direct ALSDA result

The raw direct diagnostic was evaluated as
`chiKS_compact(0,eta) Kxc_compact m00_compact - m00_compact`.  No repair,
rescaling, shifting, sign/factor tuning, BES, or GCR was applied.

| eta (Ry) | residual norm | relative residual | rigid overlap (real, imag) | field norm | response norm |
|---:|---:|---:|---:|---:|---:|
| 0.010 | 1.461615672376130e-1 | 1.880237608280108e-1 | -5.824444588445526e-2, -3.084541315396532e-2 | 3.100303726925040e-1 | 7.135531213128369e-1 |
| 0.005 | 1.412558597224824e-1 | 1.817130076392633e-1 | -5.683339423281615e-2, -1.546530543876862e-2 | 3.100303726925040e-1 | 7.145419332741448e-1 |

The artifact also writes the reconstructed point-space residual norm for every
positive-measure radial point, grouped by site and response `L` (4,940 rows
for the two η values), so the compact residual has an explicit site/`L`/radial
decomposition.

The compact χKS Frobenius norms were 4.038360155976167 at η=0.010 Ry and
4.048759853834140 at η=0.005 Ry.  The optional η=0.005 value reused the same
accepted state; it did not rerun SCF.

### Independent raw compact GSR result

The independent compact LCMM solve used `ZGELSS` with `RCOND=-1`, 232
equations, and 12 retained L=0 unknowns.  Both raw solves were full rank but
were honestly marked `BLOCKED` because the overdetermined equation residual
exceeded `1e-9`; no least-squares result was promoted or repaired.

| eta (Ry) | rank | condition | equation relative residual | full relative residual | max component | status |
|---:|---:|---:|---:|---:|---:|---|
| 0.010 | 12/12 | 3.180697495543320e14 | 5.298318856896410e-3 | 1.304410713117758e8 | 7.902006525197545e7 | `BLOCKED: compact GSR equation residual exceeds tolerance` |
| 0.005 | 12/12 | 3.185178913036662e14 | 5.303295328296958e-3 | 1.306937896604510e8 | 7.916904122738202e7 | `BLOCKED: compact GSR equation residual exceeds tolerance` |

The GSR status is a raw diagnostic outcome, not a physical interpretation.

### TDVK-06 scope result

The executable wrote `TDVK-06 direct ALSDA static diagnostic = EXECUTED`,
`TDVK-06 independent GSR raw solve = EXECUTED`, and the exact terminal marker
`TDVK-06 PASS CANDIDATE` for completion of this evidence gate.  The GSR
equation blocker above remains part of the candidate evidence.  Denominator
diagnostics were explicitly not performed because the compact Dyson
representation is not certified; they remain for TDVK-07.  No spectrum,
mode, or magnon claim is made.

### TDVK-06G compact GSR action-consistency closure

The raw TDVK-06 GSR table immediately above is retained as pre-closure
evidence.  The live preflight started at exact HEAD
`647e85811cd74699ea5231cecd823176d14242a3` on branch `fable_v4`; only
pre-existing untracked design/fixture artifacts were present.  The suspected
inconsistency was confirmed: the old GSR columns were built as
`P K_j m00_point`, while the final compact interaction action was
`P K(u) R P m00_point`.

The corrected contract is now used in both places:

```text
c = P m00_point
m00_compact_point = R c
f_j = Kc_j c = P K_j m00_compact_point
Gamma_j = chiKS_compact f_j
Kc(u)c = sum_j u_j f_j
```

The GSR SVD policy is unchanged: `ZGELSS`, `RCOND=-1`, the same unknown
parameterization, rank detection, interaction basis, and physical state.  No
regularization, singular-value tuning, coefficient constraint, rescaling, or
zero-mode enforcement was added.  Occupation provenance is stated literally:
occupations are deterministically reconstructed from the accepted reciprocal
state using the same EF, temperature, and Fermi function, followed by the
accepted eigenvector/POTPAR large-component/frozen-core projection.

#### Magnetization-span closure

The weighted metric is the existing compact radial metric and is reported for
the normalized `m00=sqrt(4*pi)*m_pauli` point vector.  The accepted 232-mode
production span gives:

| component | weighted norm | projection residual norm | relative residual |
|---|---:|---:|---:|
| total | 7.773573515982192e-1 | 8.409024617442748e-4 | 1.081745043017100e-3 |
| valence | 7.685880466137061e-1 | 7.601625625367819e-5 | 9.890377112758315e-5 |
| frozen core | 1.206209398815398e-2 | 7.755761936722918e-4 | 6.429863624292558e-2 |

These are diagnostics only.  The total physical magnetization, including the
core contribution, was retained unchanged.

#### Independent action oracles

`UnitLrCompactGsrActionConsistency` independently reconstructs point space and
evaluates `P K_j R c` with nested loops for several retained L=0 operators,
then compares that result with `Kc_j c`.  It exercises both a physical rigid
magnetization vector and a deterministic mixed complex compact vector.  It
also compares a nontrivial assembled `sum_j u_j f_j` with the assembled
compact operator action and calls the live GSR builder so the old inconsistent
column construction would fail.

| oracle | result |
|---|---:|
| physical single-operator max difference | 1.42108547e-14 |
| mixed compact single-operator max difference | 3.03692168e-14 |
| assembled action norm (sum / matrix) | 3.80303612e1 / 3.80303612e1 |
| assembled action absolute difference | 4.78298717e-15 |
| assembled action relative difference | 1.76925424e-16 |
| live Gamma-column max difference | 7.10542736e-15 |

The unit fixture also reports solve-vs-reconstructed residual difference
relative `2.01321542e-16`; the production run records the corresponding raw
conditioning-amplified values below.

#### Corrected Fe GSR result

The same accepted 12³ state was rerun at both requested etas.  The `solve`
columns are the actual `Gamma*u-g` residual; `reconstructed` is
`chiKS_compact*Kc(u)c-g`.  `difference` is their vector difference.  The
assembled-action relative value is the independent matrix-action check.

| eta (Ry) | rank | singular min/max | condition | coefficient norm | solve norm / relative | reconstructed norm / relative | residual difference norm / scaled relative | assembled action relative | runtime (s) | status |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 0.010 | 12/12 | 1.804605412204494e-13 / 6.237019428771885e1 | 3.456167972561261e14 | 1.478410860198400e10 | 5.298462141362839e-3 / 6.815996826336072e-3 | 5.311556938395065e-3 / 6.832842109482741e-3 | 4.132700481968813e-4 / 1.597760230141209e-14 | 2.553137732836307e-14 | 0.108888 | `BLOCKED: GSR rank/conditioning` |
| 0.005 | 12/12 | 1.804116193828506e-13 / 6.245659388557259e1 | 3.461894200563311e14 | 1.474764177022102e10 | 5.302883042908892e-3 / 6.821683920082049e-3 | 5.369965817435924e-3 / 6.907979899194294e-3 | 8.185585709960676e-4 / 3.163704027151049e-14 | 7.741112665751001e-14 | 0.109459 | `BLOCKED: GSR rank/conditioning` |

The two rows each have 232 compact equations and 12 unknowns.  The remaining
per-eta raw fields are:

| eta (Ry) | residual-difference max component | max residual component | rigid overlap (real, imag) |
|---:|---:|---:|---:|
| 0.010 | 3.305468880612335e-4 | 3.782271659661408e-3 | -3.108013805113175e-4, 4.648594557940931e-5 |
| 0.005 | 6.535864189757058e-4 | 3.784947756599312e-3 | 6.207431728365203e-4, -1.425601111746728e-5 |

The unchanged SVD setting is recorded as `svd_rcond=-1`; the effective
machine-precision cutoffs are `1.384896514971398e-14` and
`1.386814971428509e-14`.  The compact assembled action itself agrees at
roundoff relative to its approximately `2.59e10` norm.  Its absolute
roundoff is amplified by the approximately `1.48e10` least-squares
coefficients, which explains the nonzero residual-difference norms without
reintroducing the former orders-of-magnitude representation mismatch.

The raw direct ALSDA rows and settings above are unchanged; no ALSDA source,
normalization, channel, sign, prefactor, or physical eta was modified.  The
corrected GSR outcome is therefore a meaningful numerical blocker from rank/
conditioning, not a claim about Goldstone physics.  The compact action oracle
passed, so this is not `BLOCKED — COMPACT GSR OPERATOR ACTION`; the measured
magnetization projection defect is reported but no response-space or core
remediation is authorized.

The TDVK-06G result is:

```text
TDVK-06 GSR CLOSURE PASS CANDIDATE
BLOCKED — GSR RANK/CONDITIONING
```

No denominator, spectrum, mode, or TDVK-07 work was started.
