# DRESP-09ZS compact-span audit

Status: **PASS-B**. The raw arbitrary-(L,M) scalar-relativistic vertex is
closed, but the historical four-branch product space is not closed for the
complete six-branch radial construction. The required next milestone is
therefore a product-basis extension, not a change to the accepted observable
or to any response solver.

## Scope and guardrails

The audit is an accepted-state radial/product-space diagnostic for bcc Fe:

| item | value |
|---|---|
| state | bcc Fe, 4x4x4 (64 k points), 300 K, spd |
| radial state | orbital `lmax=2`, response `Lmax=4`, Hamiltonian-only second order |
| historical branches | `00/10/01/11` |
| shadow branches | `00/10/01/11/20/02` |
| metric | production LR-04 response-space radial weights |
| decomposition | direct LAPACK `zgesvd` on `W**(1/2) B D**(-1)` |
| downstream physics | not entered: no `chi0`, Dyson, ALSDA, Ward, spectra, BES, or Halle |

The implementation is a diagnostic backend only. It does not replace or
modify the production response basis, `chi0`, Dyson, ALSDA, Ward, BES, or
Halle paths.

## Candidate inventory and retained rank

Candidate counts are inventories, not dimensions. The historical four-branch
inventory has 232 candidates and the six-branch shadow has 348 candidates;
the new `20/02` branches contribute 116 candidates. With

`tau1 = max(npoint,ncandidate) * epsilon * sigma_max`,

and `tau10=10*tau1`, `tau100=100*tau1`, all candidates remain retained in
this Fe state at all three thresholds:

| product `K` | old candidates | old rank (`tau1/tau10/tau100`) | six candidates | six rank (`tau1/tau10/tau100`) |
|---:|---:|---:|---:|---:|
| 0 | 12 | 12 / 12 / 12 | 18 | 18 / 18 / 18 |
| 1 | 16 | 16 / 16 / 16 | 24 | 24 / 24 / 24 |
| 2 | 16 | 16 / 16 / 16 | 24 | 24 / 24 / 24 |
| 3 | 8 | 8 / 8 / 8 | 12 | 12 / 12 / 12 |
| 4 | 4 | 4 / 4 / 4 | 6 | 6 / 6 / 6 |
| **total** | **232** | **232 / 232 / 232** | **348** | **348 / 348 / 348** |

The rank increase is 116 and is threshold-stable. The smallest six-branch
singular values by `K` are `1.49e-9`, `1.24e-10`, `1.39e-10`, `2.14e-7`, and
`1.60e-5`; the corresponding `tau1` values are approximately
`2.53e-13`, `2.65e-13`, `2.91e-13`, `2.20e-13`, and `2.05e-13`.

The principal-angle minimum singular value is `1.0` for every `K`, confirming
that the historical space is contained in the six-branch span to numerical
precision. The maximum measured old-span-outside-six-span residual is
`1.76e-9` (RMS `3.03e-11`). The retained new-mode right-singular-vector
composition contains both new branches at every `K`; their aggregate weights
(`20`, `02`) are:

```text
K=0  0.2159846  0.1602988
K=1  0.1349551  0.1482125
K=2  0.1486553  0.1276535
K=3  0.1008874  0.1433607
K=4  0.2188810  0.0518053
```

## Direct projection and physical activity

Direct projection of normalized `20` and `02` candidates into the historical
space gives maximum residuals of `1.41e-3` and `3.91e-3`, respectively. The
medians are `1.55e-6` and `1.55e-6`; the per-(l,l') maxima are:

| l,l' | max R20 | max R02 |
|---|---:|---:|
| 0,0 | `4.58e-5` | `5.06e-5` |
| 0,1 | `7.27e-7` | `2.01e-6` |
| 0,2 | `5.25e-7` | `8.06e-8` |
| 1,0 | `1.99e-6` | `7.45e-7` |
| 1,1 | `8.44e-5` | `8.57e-5` |
| 1,2 | `3.43e-4` | `2.67e-6` |
| 2,0 | `4.20e-8` | `6.85e-7` |
| 2,1 | `4.76e-6` | `1.95e-4` |
| 2,2 | `1.41e-3` | `3.91e-3` |

The full six-candidate residual reaches `1.70e-3` over `K`. The arbitrary-(L,M)
physical sweep is active: the historical-span residual maxima are

| diagnostic | maximum residual |
|---|---:|
| physical SR field, `L=0..4` | `7.83e-1` |
| analytic augmentation delta-O, `L=0..4` | `6.04e-1` |
| mixed `L,M` full-SR field | `1.03e-4` |

The delta-O sweep uses the existing analytic augmentation-frame branch tangent
with the arbitrary-(L,M) Gaunt assembly. This keeps the test bounded while
covering all six radial branches and all `L<=4` source vertices.

The field/density duality residuals for single, mixed, circular, and real
deterministic cases are below `3.8e-16`. The intentionally incorrect
double-weighting control is `9.89e-1`, so the check detects a second metric
application.

The frozen scalar L=0 regression remains `1.11e-16`; it is retained as a
regression/oracle while the complete orbital tensor is audited separately.

## Classification

```text
primary_classification = RAW_ARBITRARY_L_VERTEX_CLOSED_PRODUCT_BASIS_EXTENSION_REQUIRED
verdict = PASS-B
recommended_next_milestone = DRESP-10A
```

This is not a numerical block. It means the radial/SVD audit is reproducible
and physically active, while production still uses the historical basis until
a reviewed compact six-branch basis extension is implemented.

## Reproduction and tests

The backend is selected by
[`input_dresp09zs_fe.nml`](../tests/integration/tddft_driver_smoke/input_dresp09zs_fe.nml).
The independent validator is
[`dresp09zs_fe_artifact.py`](../tests/validation/dresp09zs_fe_artifact.py).

```text
cmake --build build --target rslmto.x -j2
ctest --test-dir build -R 'Dresp09ZSFeArtifact|UnitLrLmtoProductResponseBasis|UnitDresp09sScalarRelativistic|UnitDresp09xSrSpinObservable|UnitDresp09YAugmentationTangent|UnitDresp09ZAngularVertex' --output-on-failure
```

The live artifact is written to `/tmp/dresp09zs_fe_4k.dat`; it is not a
production input or frozen baseline.
