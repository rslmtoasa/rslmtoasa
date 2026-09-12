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
