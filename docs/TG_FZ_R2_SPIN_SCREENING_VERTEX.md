# TG-FZ-R2 — Spin-Dependent Screening Covariance of the Exchange Vertex

Status: **PASS-A — exact fixed-z vertex covariance derived and verified**.

Representation note from TG-FZ-R7: the original R2 ledger below predates the
screening-normalization audit and uses the historical raw-labelled `gamma`
quantity `qpar`. In the repaired implementation, every `P` built with
`dele` uses `gamma_norm=qi`; the covariance algebra is unchanged. The
corrected normalized-gamma numerical ledger is recorded in
`docs/TG_FZ_R7_SCREENING_REPRESENTATION_REPAIR.md`.

This audit continues from commit `9a81106` (TG-FZ-R1). It is restricted to
the existing two-site, spd, fixed-complex-energy diagnostic fixture. It adds
no energy integration, native `J_ij`/`J(q)`, Fe-shell comparison, TT+C,
finite-temperature, or DRESP-04 work.

## Derived transformation

For each site and spin:

```text
D_sigma       = alpha - gamma_sigma
R_sigma       = I - D_sigma P_gamma_sigma
P_alpha_sigma = P_gamma_sigma R_sigma^(-1)
```

The off-site path operator is transformed by

```text
g_alpha_sigma,ij = R_i,sigma g_gamma_sigma,ij R_j,sigma.
```

For ordered `ud`, cyclicity of the trace gives

```text
C_alpha = Im Tr[DeltaP_i^alpha g_alpha_up,ij
                 DeltaP_j^alpha g_alpha_down,ji]

C_alpha = Im Tr[tildeDelta_i,ud g_gamma_up,ij
                 tildeDelta_j,ud g_gamma_down,ji]

tildeDelta_i,ud = R_i,down DeltaP_i^alpha R_i,up
tildeDelta_j,ud = R_j,up   DeltaP_j^alpha R_j,down.
```

The independently evaluated ordered `du` form uses

```text
tildeDelta_i,du = R_i,up   DeltaP_i^alpha R_i,down
tildeDelta_j,du = R_j,down DeltaP_j^alpha R_j,up.
```

For the diagonal scalar-relativistic local fixture:

```text
R_down DeltaP_alpha R_up
  = DeltaP_gamma + (D_up-D_down) P_gamma_up P_gamma_down
  = DeltaP_gamma + (gamma_down-gamma_up)
                         P_gamma_up P_gamma_down.
```

The test evaluates direct endpoint products, split transformed-P products, and
the raw-gamma-plus-correction candidate independently.

## Coefficient-space diagnostic

The historical finite-H coefficient vertex remains unchanged:

```text
d = W_up DeltaP_gamma W_down
d_tilde_i = W_up tildeDelta_i,ud W_down
d_tilde_j = W_up tildeDelta_j,ud W_down
delta_d_screen = d_tilde - d
```

`d_matrix` and the finite-H oracle are not changed. The nonzero
`delta_d_screen` measures the spin-dependent screening correction; it is not
used as a fitted scale or replacement vertex.

## Common-gamma control

After each physical z point, an isolated local control retains the physical
spin-dependent `P_gamma` blocks and unchanged potential inputs, including
their spin-dependent `c`/`dele` dependence. Only the diagnostic matrices
are made equal:

```text
D_up = D_down.
```

The correction is evaluated from
`(D_up-D_down) P_gamma_up P_gamma_down`, not inserted as a forced zero. The
identity therefore predicts that the transformed vertex becomes raw
`DeltaP_gamma` and direct alpha/raw-gamma contractions agree.

## Per-z numerical ledger

These are the three existing fixed energies. Matrix columns are maximum
absolute residuals; ordered `ud` and `du` contractions are kept separate.

| z | endpoint expansion | screening identity | alpha vs transformed gamma (ud / du) | raw gamma vs transformed (ud / du) | max `|delta_d_screen|` | correction | raw vertex | transformed vertex |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `-0.91+0.83i` | `4.054610e-14` | `5.411332e-14` | `4.912791e-20 / 5.421011e-20` | `1.798191e-06 / 1.798191e-06` | `5.424334e-01` | `3.605760e+01` | `1.317197e+01` | `2.342778e+01` |
| `-0.17+0.04i` | `1.831027e-15` | `6.661338e-16` | `4.440892e-16 / 8.881784e-16` | `7.649156e-02 / 7.649156e-02` | `2.858181e-03` | `1.899941e-01` | `9.615014e+00` | `9.434799e+00` |
| `0.62+0.31i` | `1.421085e-14` | `7.944109e-15` | `2.168404e-18 / 3.252607e-18` | `1.035917e-03 / 1.035917e-03` | `2.815895e-01` | `1.871832e+01` | `2.251586e+01` | `3.961667e+01` |

The complete alpha direct versus transformed-gamma contraction closes to
`8.881784e-16`. Raw gamma is not interchangeable with the transformed
vertex: its maximum mismatch is `7.649156e-02`.

## Raw gamma versus finite-H

The raw gamma contraction remains the independently constructed finite-H
result at machine precision:

| z | raw gamma = finite-H (ud) | `|gamma-H|_ud` | `|gamma-H|_du` |
|---|---:|---:|---:|
| `-0.91+0.83i` | `-5.380807e-06` | `5.505714e-20` | `5.336308e-20` |
| `-0.17+0.04i` | `-2.202285e+00` | `4.440892e-16` | `4.440892e-16` |
| `0.62+0.31i` | `-3.412346e-04` | `2.168404e-19` | `7.589415e-19` |

The maximum raw-gamma/finite-H residual is `4.440892e-16`; the remaining
physical mismatch is therefore a vertex-screening issue, not a path-operator
or finite-H failure.

## Common-gamma results

| z | vertex residual | correction | alpha minus raw gamma (ud / du) |
|---|---:|---:|---:|
| `-0.91+0.83i` | `1.803352e-14` | `0.000000e+00` | `2.126053e-19 / 1.939705e-19` |
| `-0.17+0.04i` | `1.779823e-15` | `0.000000e+00` | `1.332268e-15 / 0.000000e+00` |
| `0.62+0.31i` | `5.687117e-15` | `0.000000e+00` | `1.029992e-18 / 5.963112e-19` |

The maximum common-gamma vertex residual is `1.803352e-14`, and the maximum
direct alpha/raw-gamma contraction residual is `1.332268e-15`.

## Verification and disposition

GNU Fortran 13.3.0 Debug checks used `-g -O0 -fbacktrace
-fcheck=all,no-recursion`. Both the focused test and the nearby 8-test
selection passed:

```bash
cmake --build build -j2
ctest -V -R UnitDresp03tgNativeFixedZ
ctest -R 'UnitDresp03|UnitLrKlStaticBridge' --output-on-failure
```

Preserved R1 gates include P transformation `2.842171e-14`, S round trip
`2.220446e-15`, path covariance `6.497414e-16`, gamma versus finite-H
`4.440892e-16`, and Pauli-helper self-consistency `0`.

TG-FZ-R3 is retired for this campaign: the R2 transformation closes the
fixed-z covariance gates and no separate Pauli contradiction remains.

`source/exchange.f90` was not modified. Its SHA-256 remains:

```text
6e9ae7da6af46367aaf18a4fef37bdd4a9a0d6f69fe8e08b589e5c5f278a70e6
```
