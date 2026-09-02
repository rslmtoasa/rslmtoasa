# TDDFT-11 — Fe/Ni transverse-response validation

**Campaign date:** 2026-09-02

**Release-gate status:** **FAIL — not validated for production release**

This report follows the release criteria in the TD-DFT master blueprint. It
records a bounded current-branch probe, not a release waiver. No empirical
frequency shift, spectral rescaling, or Goldstone cleanup was applied.

## Scope and reproducibility

The probe uses the collinear, scalar-relativistic boundary required by the
blueprint: `nsp=1`, SOC disabled, `reciprocal_mode='ham_only'`, site-projected
transverse response, first-order Hamiltonian, and block intersite recursion.
The Fe and Ni starting potentials are read from the existing validation
restart databases; these decks use `nstep=1` and therefore do not claim a new
fully converged ground-state campaign.

To keep native real-space storage bounded in this workspace, the reproducible
probe uses `Fe rc=9` (258 generated central-cell-to-cluster pairs) and `Ni
pbc=4x4x4` (64 pairs). Both materials use the same three direct q points
`(0,0,0)`, `(0.01,0,0)`, `(0.02,0,0)` and the three-frequency grid
`0, 0.001, 0.002 Ry`, with `eta=0.0005 Ry`. The native provider emits bare
`chi0` only because its exact static kernel is not yet available.

Run each material from its directory:

```text
OMP_NUM_THREADS=1 ../../../../../build/bin/rslmto.x input_eigenpairs.nml
OMP_NUM_THREADS=1 ../../../../../build/bin/rslmto.x input_kspace_lehmann.nml
OMP_NUM_THREADS=1 ../../../../../build/bin/rslmto.x input_realspace_gf.nml
```

The checked-in raw outputs are collected by
[`collect_tddft11_evidence.py`](../tests/regression/tddft_validation/materials/collect_tddft11_evidence.py)
into [`evidence.json`](../results/validation/TDDFT-11_FE_NI/evidence.json).

## Results

All six bounded runs completed and emitted manifests plus three q-resolved
`chi0` files. The table reports the maximum absolute complex matrix-element
difference over all emitted q and frequency rows; the target pointwise
equivalence tolerance is `1e-8`.

| material | eigenpairs vs K-space Lehmann | eigenpairs vs native RS | K-space Lehmann vs native RS | backend gate |
|---|---:|---:|---:|---|
| bcc Fe | 46.36215447 | 6.93855871 | 39.61863047 | **FAIL** |
| fcc Ni | 44.63152833 | 61.13543750 | 16.91914306 | **FAIL** |

The native runs do exercise the intended `G_ab(R,z), G_ba(-R,z) →
chi0(R,omega) → chi0(q,omega)` path. They retain 258 Fe real-space points to
an 8.5836 Å cutoff and 64 Ni points to a 23.0363 Å cutoff. Neither probe has
an R-space tail sweep, so real-space convergence is explicitly **not
assessed**.

The K-space raw static diagnostics also fail the gapless test. With a
`1e-6` residual tolerance, the largest raw Ward residual is `0.2570132937`
for Fe and `0.0321499941` for Ni. The corresponding closest raw eigenvalues
are `1.2570132937` and `1.0321499941`; the native bare-dynamic provider has no
static Goldstone diagnostic. These are recorded in the raw Goldstone files,
not altered in post-processing.

## Implementation work completed

The native production attachment had three integration defects exposed by this
campaign and fixed on the current branch:

1. A missing explicit `ijpair` list now causes the native bulk stack to build
   the complete central-cell-to-cluster pair list before recursion and Green
   storage are constructed. Explicit user pair lists remain untouched.
2. The native branch restores pair ownership before intersite recursion and
   returns to response-site ownership afterward, so all generated pairs are
   actually evaluated.
3. Native bare-`chi0` output no longer appends dynamic Gamma metadata to a
   Goldstone file it deliberately does not create; common q-loop deallocation
   is also guarded for the native path.

## Release-gate checklist

- [ ] Fe passes backend-equivalence gate.
- [ ] Ni passes backend-equivalence gate.
- [ ] SOC=0 acoustic mode is gapless within converged numerical uncertainty.
- [ ] Small-q dispersion is demonstrably quadratic.
- [ ] Stiffness agrees with independent low-energy references within established uncertainty.
- [ ] Ni path/backfolding concerns are explicitly resolved.
- [ ] Damping/Stoner trends are stable under numerical convergence.
- [x] Release-gate report states **FAIL** without cosmetic tuning.

## Required follow-up before release

The campaign must still resolve the three-route normalization/integration
disagreement, provide a validated native static-limit kernel, and rerun with
freshly converged Fe/Ni ground states. Only then should it add dense small-q
mode paths and a demonstrated `omega=Dq^2` window, same-ground-state LKAG and
frozen-magnon/GBT references, Ni reciprocal-path/backfolding checks, and
eta/k-mesh/frequency/R-space convergence sweeps.
