# TDDFT-R2-05 — bcc Fe and fcc Ni validation

**Campaign date:** 2026-09-02
**Overall release decision:** **FAIL — neither material is production-passed**

This report follows
[TDDFT-R2-05](../dev/RS_LMTO_TDDFT_R2_PROMPT_PACK_2026-09-02/TDDFT-R2-05_BCCFE_FCCNI_VALIDATION.md).
It records what was actually measured and keeps failed gates visible. No
Goldstone correction, empirical frequency shift, spectral rescaling, or
synthetic material fixture was used.

## Material decisions

| Material | Decision | Basis |
| --- | --- | --- |
| bcc Fe | **PARTIAL** | The current three-route selected-point run completes, native static Ward diagnostics are available, and historical convergence/reference campaigns exist. The raw Ward residual is not gapless, the three bare-response routes disagree, and the dense small-q branch is not quadratic. |
| fcc Ni | **FAIL** | The three routes disagree, the native static Ward residual is 1, no coherent q² branch exists, and the retained VAL-19 restart has an atomic-like/no-hopping neighbor map rather than a trustworthy itinerant-Ni reference. |

These are evidence classifications, not release approvals. The overall
TDDFT material gate remains closed.

## Scope and reproducibility

The recognized current material fixtures are the three decks per material in
[`tests/regression/tddft_validation/materials`](../tests/regression/tddft_validation/materials/README.md).
They use the existing validation restart data:

- Fe: `results/validation/VAL-18_bccFe/dispersion_nk16/`, `alat=2.86120`,
  `rc=9`, 258 central-cell-to-cluster pairs.
- Ni: `results/validation/VAL-19_fccNi/scf_weak/`, `alat=6.650`, `pbc=4x4x4`,
  `ct=4`, `r2=16`, 64 pairs.

The response boundary is collinear, scalar-relativistic, SOC-free,
`ham_only`, first-order, site-projected transverse response. `nstep=1` reads
the restart state and does not claim a new SCF convergence campaign. The
selected q list is direct-coordinate
`(0,0,0)`, `(0.01,0,0)`, `(0.02,0,0)`; the raw response grid is
`omega=(0,0.001,0.002) Ry`, `eta=0.0005 Ry`, `T=300 K`.

The six current-branch commands were run with `OMP_NUM_THREADS=1`:

```text
OMP_NUM_THREADS=1 ../../../../../build/bin/rslmto.x input_eigenpairs.nml
OMP_NUM_THREADS=1 ../../../../../build/bin/rslmto.x input_kspace_lehmann.nml
OMP_NUM_THREADS=1 ../../../../../build/bin/rslmto.x input_realspace_gf.nml
```

Build provenance was GNU Fortran 13.3.0, `RelWithDebInfo`, OpenMP 4.5,
MKL kernels enabled, MPI disabled. The raw outputs were generated from
`bf91262` plus the local validation fixes committed with this report. The
machine-readable bundle is
[`results/validation/TDDFT-11_FE_NI/evidence.json`](../results/validation/TDDFT-11_FE_NI/evidence.json).
It hashes the primary `plus_minus` outputs; current `minus_plus` files are
retained alongside them because the default deck requests both circular
channels.

## Convergence hierarchy

The R2-05 selected-point decks are diagnostic production attachments, not a
substitute for the existing high-accuracy campaign series. Each axis was not
varied simultaneously.

### bcc Fe

| Axis | Evidence | Assessment |
| --- | --- | --- |
| k mesh | `4^3`, `8^3`, `12^3`, `16^3`; legacy raw Ward residuals `0.25699`, `0.0240553`, `0.00506030`, `0.0189444` | **PARTIAL**; non-monotonic Ward behavior prevents convergence claim |
| eta | `0.0005`, `0.001`, `0.002 Ry`; observed widths `0.0014`, `0.0022`, `0.0028 Ry` | **PARTIAL**; zero-eta extrapolation was not accepted as intrinsic damping |
| frequency grid | high-accuracy campaign fixed at 401 points over `0..0.08 Ry`; R2 attachment uses 3 points over `0..0.002 Ry` | **NOT CONVERGED**; no spacing/grid sweep |
| energy/band window | production window recorded as `band_first=1, band_last=0`; no independent window promotion sweep | **NOT CONVERGED** |
| native `R_source`, `R_max` | 258 points, source radius `8.5836 Å`, requested `R_max=30 Å`, omitted points `0` | **NOT CONVERGED**; `real_space_tail_assessed=F` and source does not cover request |
| contour points | direct real-axis route; `N_contour=0` | **NOT APPLICABLE** to this attachment; no contour sweep |

Historical Fe small-q evidence used five direct q points
`0.01, 0.015, 0.02, 0.03, 0.04`. The accepted legacy individual features
were near `0.0032–0.0034 Ry`, but the zero-intercept fit gave
`D=2.91453 Ry` with relative residual `0.65317`; this is rejected as a
quadratic spin-wave branch. The independent GBT/frozen-magnon diagnostic gave
`D=0.161163` in its own reported convention with residual `0.05046`; the Jij
record contains only two pair records and cannot provide a stiffness.

### fcc Ni

| Axis | Evidence | Assessment |
| --- | --- | --- |
| k mesh | `4^3`, `8^3`, `12^3`; sampled raw legacy residual stable at `0.03237984` | **PARTIAL**; stable diagnostic, not a valid itinerant reference |
| eta | `0.0005`, `0.001`, `0.002 Ry`; observed q=0.01 widths `0.0011`, `0.0021`, `0.0041 Ry` | **PARTIAL**; width follows numerical broadening |
| frequency grid | historical route fixed at 201 points over `0..0.020 Ry`; R2 attachment uses 3 points over `0..0.002 Ry` | **NOT CONVERGED**; no grid-spacing sweep |
| energy/band window | bands `1:18` and `1:14` showed no printed-precision change | **PARTIAL**; no broader reference-window study |
| native `R_source`, `R_max` | 64 points, source radius `23.0363 Å`, requested `R_max=30 Å`, omitted points `0` | **NOT CONVERGED**; `real_space_tail_assessed=F` and source does not cover request |
| contour points | direct real-axis route; `N_contour=0` | **NOT APPLICABLE** to this attachment; no contour sweep |

The old non-q² behavior was explicitly investigated. Smearing moved the
apparent q=0.01 feature from `0.0040287 Ry` at 100 K to `0.0042178 Ry` at
300 K and `0.0074031 Ry` at 600 K. No accepted coherent Xi crossing exists,
the projected feature weight is negative, and no q² fit is made. The current
restart has a 2.000 μB moment and its nearest fcc image is 4.7027 Å, outside
the retained neighbor map; an enlarged-cutoff diagnostic instead exposed a
non-Hermitian reciprocal Hamiltonian. This is not a production itinerant-Ni
ground state.

## Three-backend bare-chi0 comparison

The comparison is made before Dyson enhancement on the primary `plus_minus`
complex response matrices. Values below are the normalized Frobenius norm

`||A-B|| / max(||A||,||B||)`

for the indicated q/omega point. The full-grid maxima cover all 3 q points and
all 3 frequencies; the machine-readable per-point values are in
`backend_equivalence.representative_relative_matrix_norm`.

| Material / point | eig vs k-space GF | eig vs native RGF | k-space GF vs native RGF |
| --- | ---: | ---: | ---: |
| Fe Γ, 0 Ry | 0.51557 | 0.09723 | 0.46348 |
| Fe q=0.01, 0 Ry | 0.77407 | 0.09451 | 0.75049 |
| Fe full 3×3 grid max | 0.79335 | 0.09698 | 0.77208 |
| Ni Γ, 0 Ry | 0.72325 | 1.00000 | 1.00000 |
| Ni q=0.01, 0 Ry | 0.72325 | 1.00000 | 1.00000 |
| Ni full 3×3 grid max | 0.72884 | 1.00000 | 1.00000 |

The collector’s pointwise tolerance is `1e-8`; all three pairings fail for
both materials. The selected-point attachment emits bare chi0 only and does
not contain a converged interacting magnon-region/Stoner-region spectrum, so
those required representative spectral comparisons remain open rather than
being inferred from a trace peak.

## Raw Ward/Goldstone evidence

The raw calculation precedes any correction. The relevant mode has unit
magnetization overlap in all three routes. The raw static diagnostics are:

| Material / route | closest eigenvalue | raw residual | magnetization overlap |
| --- | ---: | ---: | ---: |
| Fe eigenpairs | `1.25701329` | `0.25701329` | `1.0` |
| Fe k-space GF | `1.25701329` | `0.25701329` | `1.0` |
| Fe native RGF | `1.13493781 - 5.66e-6i` | `0.13493781` | `1.0` |
| Ni eigenpairs | `1.03214999` | `0.03214999` | `1.0` |
| Ni k-space GF | `1.03214999` | `0.03214999` | `1.0` |
| Ni native RGF | `-2.44e-22` | `1.00000000` | `1.0` |

The target tolerance in the evidence bundle is `1e-6`. No Goldstone correction
was applied. The static native results use the exact native real-space contour
identity now implemented by R2-04; their disagreement with the reciprocal
routes is itself release evidence. A physical Γ acoustic gap is not assigned
from the finite-eta three-point bare grid: the raw unity-mode distance is
reported above, while a converged interacting acoustic branch is absent.

## Small-q, adiabatic, and damping decision

Fe has generated small-q and independent adiabatic diagnostic products, but
the legacy q² residual is `0.65317`, the pair-potential route has no coherent
accepted branch, and the independent Jij result is incomplete. The Fe eta
sweep widths are retained as numerical-width diagnostics; no intrinsic Landau
damping value is claimed.

Ni has no accepted coherent mode at any tested q, so fitting rejected trace
peaks would be misleading. Its eta widths scale with eta, and its finite-eta
bare Stoner response overlaps the low-energy edge. No intrinsic damping,
stiffness, or LKAG/GBT/Jij agreement is claimed for Ni.

## Reproducible evidence bundle

- Current decks: [`bccFe`](../tests/regression/tddft_validation/materials/bccFe),
  [`fccNi`](../tests/regression/tddft_validation/materials/fccNi).
- Machine-readable hashes and normalized comparisons:
  [`evidence.json`](../results/validation/TDDFT-11_FE_NI/evidence.json).
- Fe high-accuracy convergence/reference campaign:
  [`VAL-18`](../results/validation/VAL-18_bccFe/campaign.json) and
  [`VAL-18 report`](../dev/VAL-18_TDDFT_BCC_FE.md).
- Ni convergence/sensitivity campaign:
  [`VAL-19`](../results/validation/VAL-19_fccNi/campaign.json) and
  [`VAL-19 report`](../dev/VAL-19_TDDFT_FCC_NI.md).
- Current raw q-resolved outputs and Goldstone files are stored beside the
  decks; reverse-channel products carry the `_minus_plus_` suffix.

## R2-05 checklist

Checked means the evidence or investigation was actually performed; an
unchecked item is an unmet acceptance gate.

### bcc Fe

- [x] high-accuracy fixture established (VAL-18 artifacts retained).
- [x] k-mesh convergence measured (non-monotonic; gate not passed).
- [x] eta convergence measured (numerical-width interpretation retained).
- [ ] frequency-grid convergence demonstrated.
- [ ] three-backend static chi0 agreement.
- [ ] three-backend dynamic chi0 agreement.
- [x] raw Ward residual reported.
- [ ] raw Γ-point acoustic gap established from a converged interacting branch.
- [ ] quadratic small-q dispersion demonstrated; measured fit is rejected.
- [x] diagnostic spin-wave stiffness extracted.
- [x] independent adiabatic comparison performed where available; Jij is incomplete.
- [x] damping/Stoner behavior examined without an intrinsic-damping claim.

### fcc Ni

- [x] high-accuracy diagnostic fixture established (reference state is not production-valid).
- [x] k-mesh convergence measured (stable but not physically validating).
- [x] eta convergence measured.
- [ ] frequency-grid convergence demonstrated.
- [ ] three-backend static chi0 agreement.
- [ ] three-backend dynamic chi0 agreement.
- [x] raw Ward residual reported.
- [ ] raw Γ-point acoustic gap established from a converged interacting branch.
- [x] old non-q² issue explicitly investigated.
- [x] quadratic small-q failure clearly documented.
- [ ] stiffness extracted from a credible collective branch.
- [ ] independent adiabatic comparison performed.
- [x] damping/Stoner behavior examined without an intrinsic-damping claim.

### Release decision

- [x] Fe assigned **PARTIAL**.
- [x] Ni assigned **FAIL**.
- [x] no synthetic fixture used as a substitute for material evidence.
- [x] no Goldstone correction used to hide a poor raw Ward result.
