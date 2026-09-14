# TDVK-03A — reciprocal-GF real-axis quadrature closure audit

**Status: Case A — `REAL-AXIS GF FORMULATION NUMERICALLY CONSISTENT`**

This is a fixture-level numerical audit. It does not change the TDVK-02 or
TDVK-03 material status, and it does not claim a Fe response calculation.

## Provenance and scope

The numerical samples were executed on 2026-09-14 from the following source
HEAD:

```text
4da3c115330cd2c7c0969c7d4418ec5f1d6a636c
```

The worktree was already dirty. Existing user changes were preserved and
were not used to change the GF equations, KXC, GSR, Dyson, Goldstone logic,
product vertices, Fe path, or response broadening. The audit executable was
built in the existing `build` directory with the Debug configuration,
GNU Fortran 13.3.0, CMake 3.28.3, OpenMP enabled, MPI disabled, and libXC
5.2.3 enabled. The executable calls the existing
`evaluate_lr_product_gf_susceptibility`, `build_weighted_resolvent`, and
`evaluate_lr_product_ks_susceptibility` services.

The audit test is registered as a disabled CTest entry because the prescribed
13,401-point sample takes several minutes on this CPU. The individual audit
modes were run directly and are listed at the end of this report.

### R3 fixture provenance

The existing nontrivial `UnitLrProductGfSusceptibility` fixture was reused:

| quantity | value |
|---|---:|
| basis / bands | `nbasis=8`, `nbands=8` |
| k points | `nk=2`, weights `0.35, 0.65` |
| spectrum | `-0.720` through `0.647 Ry` |
| temperature / Fermi level | `0.030 Ry` / `0.000 Ry` |
| response broadening | `eta=0.040 Ry` |
| radial mesh | `nr=7`, nonconstant `phi` and `phidot` |
| response space | one site, `response_lmax=2` |
| compact product dimension | `27` per circular channel |
| channels and q | `chi_plus`, Gamma; `chi_minus`, `q=(0.23,0,0)` |
| eigensystem | complete full eigensystem, orthogonal, collinear, `ham_only`, second order, no SOC |

The eigenvectors in this deliberately controlled fixture are the complete
identity eigensystem; the radial `phi/phidot` data remain nontrivial. The
finite-q endpoint was constructed from the exact folded `k+q` points. No
high-symmetry label is assigned to that finite-q sample.

## Numerical definitions

For each k point, the audit formed exact

\[
 M_p=\sum_n\epsilon_n^p|n\rangle\langle n|,\qquad
 N_p=\sum_n f(\epsilon_n)\epsilon_n^p|n\rangle\langle n|,
 \quad p=0,1,2,
\]

and integrated the existing weighted resolvent discontinuity

\[
 A^{(p)}(E)=\frac{i}{2\pi}\left[G^{(p),R}(E)-G^{(p),A}(E)\right]
\]

with the same composite Simpson weights and energy interval as the GF
evaluator. The reported residual is

\[
 r(X)=\frac{\|X^{\rm num}-X\|_F}{\max(\|X\|_F,\epsilon)}.
\]

The response columns are per frequency: `dF` is the absolute Frobenius
difference, `rF` is `dF` divided by the compact Lehmann Frobenius norm, and
`dInf` is the maximum elementwise difference.

## Energy intervals and resolution

For `energy_margin=0.60 Ry`, the interval is
`[-1.320000, 1.247000] Ry`, width `2.567000 Ry`. The wider-window intervals
are `[-1.720000, 1.647000] Ry` for `1.00 Ry` and
`[-2.720000, 2.647000] Ry` for `2.00 Ry`.

The fixed-`integration_eta` ladder used every requested odd grid through
12,801 points. The 12,801-point sample was retained because the 6,401-point
finite-frequency result had not clearly stabilized.

| N | interval margin (Ry) | h (Ry) | h / integration_eta |
|---:|---:|---:|---:|
| 101 | 0.60 | 2.567000e-02 | 2.567000e+01 |
| 201 | 0.60 | 1.283500e-02 | 1.283500e+01 |
| 401 | 0.60 | 6.417500e-03 | 6.417500e+00 |
| 801 | 0.60 | 3.208750e-03 | 3.208750e+00 |
| 1601 | 0.60 | 1.604375e-03 | 1.604375e+00 |
| 3201 | 0.60 | 8.021875e-04 | 8.021875e-01 |
| 6401 | 0.60 | 4.010938e-04 | 4.010938e-01 |
| 12801 | 0.60 | 2.005469e-04 | 2.005469e-01 |

For the `integration_eta` ladder, the selected grids were 515, 1029, 2055,
and 6401 points for `0.010`, `0.005`, `0.0025`, and `0.001 Ry`, respectively.
Their resolution ratios were `0.499416`, `0.499416`, `0.499903`, and
`0.401094`; all satisfy the requested `h/integration_eta <= 0.5` diagnostic.
For the window ladder, the grids were 6401, 8401, and 13401 for margins
`0.60`, `1.00`, and `2.00 Ry`, with ratios `0.401094`, `0.400833`, and
`0.400522`.

## Spectral-moment diagnostics

The following tables report every p and both k points. Columns `rM` are the
unweighted moment residuals.

### Fixed `integration_eta=0.001 Ry`, margin `0.60 Ry`

| N | rM0(k1) | rM1(k1) | rM2(k1) | rM0(k2) | rM1(k2) | rM2(k2) |
|---:|---:|---:|---:|---:|---:|---:|
| 101 | 2.805899e+00 | 7.608773e-01 | 7.309722e-01 | 2.805899e+00 | 9.969821e-01 | 7.019073e-01 |
| 201 | 7.301885e-01 | 6.562616e-01 | 6.221971e-01 | 7.301885e-01 | 6.750201e-01 | 6.644136e-01 |
| 401 | 5.719451e-01 | 6.186876e-01 | 5.920193e-01 | 5.719451e-01 | 5.856063e-01 | 5.179344e-01 |
| 801 | 3.337947e-01 | 4.055319e-01 | 4.631005e-01 | 3.337947e-01 | 3.920583e-01 | 4.094148e-01 |
| 1601 | 7.658929e-02 | 8.395721e-02 | 8.411719e-02 | 7.658929e-02 | 8.785189e-02 | 9.642702e-02 |
| 3201 | 1.005542e-02 | 1.206741e-02 | 1.254944e-02 | 1.005542e-02 | 1.226438e-02 | 1.270383e-02 |
| 6401 | 6.561636e-04 | 8.225739e-04 | 8.774814e-04 | 6.561636e-04 | 8.299393e-04 | 8.794966e-04 |
| 12801 | 5.728015e-04 | 6.329432e-04 | 6.582746e-04 | 5.728015e-04 | 6.310944e-04 | 6.538315e-04 |

Columns `rN` are the Fermi-weighted moment residuals.

| N | rN0(k1) | rN1(k1) | rN2(k1) | rN0(k2) | rN1(k2) | rN2(k2) |
|---:|---:|---:|---:|---:|---:|---:|
| 101 | 7.731379e-01 | 7.707676e-01 | 7.966310e-01 | 3.891148e+00 | 1.156103e+00 | 6.118278e-01 |
| 201 | 5.866218e-01 | 5.997922e-01 | 5.494502e-01 | 8.493217e-01 | 7.381991e-01 | 7.501375e-01 |
| 401 | 6.558133e-01 | 7.145892e-01 | 6.824470e-01 | 4.736384e-01 | 4.414279e-01 | 3.543726e-01 |
| 801 | 3.351052e-01 | 4.701307e-01 | 5.374658e-01 | 3.326295e-01 | 2.876000e-01 | 2.478094e-01 |
| 1601 | 7.074464e-02 | 6.419698e-02 | 5.933961e-02 | 8.229832e-02 | 1.046654e-01 | 1.178400e-01 |
| 3201 | 1.020226e-02 | 1.175896e-02 | 1.256240e-02 | 1.102154e-02 | 1.273504e-02 | 1.300107e-02 |
| 6401 | 4.202901e-03 | 1.322394e-03 | 1.233235e-03 | 2.747204e-03 | 1.409208e-03 | 1.273093e-03 |
| 12801 | 4.210716e-03 | 1.204872e-03 | 1.038387e-03 | 2.737407e-03 | 1.237721e-03 | 1.057849e-03 |

The finite-q `chi_minus` run uses the same electronic spectrum and therefore
has the same six per-k moment values as the 6401-point row above. It was
nevertheless evaluated independently at `q=(0.23,0,0)` before its GF bubble.

### Integration-eta moment ladder

| integration_eta (Ry) | N | h / eta_int | rM0/rM1/rM2 (k1) | rM0/rM1/rM2 (k2) |
|---:|---:|---:|---|---|
| 0.0100 | 515 | 4.994163e-01 | 6.434666e-03 / 7.091663e-03 / 7.407487e-03 | 6.434666e-03 / 6.973744e-03 / 7.170152e-03 |
| 0.0050 | 1029 | 4.994163e-01 | 2.689601e-03 / 3.233288e-03 / 3.521436e-03 | 2.689601e-03 / 3.019494e-03 / 3.110572e-03 |
| 0.0025 | 2055 | 4.999026e-01 | 2.419736e-03 / 2.595338e-03 / 2.689017e-03 | 2.419736e-03 / 2.541797e-03 / 2.582929e-03 |
| 0.0010 | 6401 | 4.010938e-01 | 6.561636e-04 / 8.225739e-04 / 8.774814e-04 | 6.561636e-04 / 8.299393e-04 / 8.794966e-04 |

| integration_eta (Ry) | N | rN0/rN1/rN2 (k1) | rN0/rN1/rN2 (k2) |
|---:|---:|---|---|
| 0.0100 | 515 | 4.170435e-02 / 1.295278e-02 / 1.142364e-02 | 2.757926e-02 / 1.263623e-02 / 1.080715e-02 |
| 0.0050 | 1029 | 2.069439e-02 / 6.521461e-03 / 5.861147e-03 | 1.320650e-02 / 5.317887e-03 / 4.247617e-03 |
| 0.0025 | 2055 | 1.095662e-02 / 4.059529e-03 / 3.731460e-03 | 7.404576e-03 / 3.794413e-03 / 3.364462e-03 |
| 0.0010 | 6401 | 4.202901e-03 / 1.322394e-03 / 1.233235e-03 | 2.747204e-03 / 1.409208e-03 / 1.273093e-03 |

### Energy-margin moment ladder

The grid spacing was held at or below the 0.60-Ry baseline spacing.

| margin (Ry) | N | h / eta_int | rM0/rM1/rM2 (k1) | rM0/rM1/rM2 (k2) |
|---:|---:|---:|---|---|
| 0.60 | 6401 | 4.010938e-01 | 6.561636e-04 / 8.225739e-04 / 8.774814e-04 | 6.561636e-04 / 8.299393e-04 / 8.794966e-04 |
| 1.00 | 8401 | 4.008333e-01 | 1.774434e-04 / 2.037556e-04 / 2.157723e-04 | 1.774434e-04 / 1.995561e-04 / 2.073108e-04 |
| 2.00 | 13401 | 4.005224e-01 | 3.432883e-04 / 3.940025e-04 / 3.740642e-04 | 3.432883e-04 / 4.009602e-04 / 4.060937e-04 |

| margin (Ry) | N | rN0/rN1/rN2 (k1) | rN0/rN1/rN2 (k2) |
|---:|---:|---|---|
| 0.60 | 6401 | 4.202901e-03 / 1.322394e-03 / 1.233235e-03 | 2.747204e-03 / 1.409208e-03 / 1.273093e-03 |
| 1.00 | 8401 | 4.159293e-03 / 9.068639e-04 / 6.730652e-04 | 2.585062e-03 / 9.388578e-04 / 6.947546e-04 |
| 2.00 | 13401 | 4.239600e-03 / 1.065539e-03 / 8.055067e-04 | 2.643914e-03 / 1.160873e-03 / 9.867338e-04 |

The moment results separate the effects cleanly: mesh refinement removes the
large Simpson error, the eta ladder decreases the finite-width error, and
the window ladder removes only the smaller Lorentzian-tail contribution.
The remaining `rN0` floor near `4e-3` at `integration_eta=0.001 Ry` is a
finite-width effect, not a claim that the broadened spectral function is an
exact delta distribution.

## GF versus compact Lehmann: fixed-eta resolution ladder

These rows use `chi_plus`, Gamma, `eta=0.040 Ry`,
`integration_eta=0.001 Ry`, and `energy_margin=0.60 Ry`.

| N | interval (Ry) | h / eta_int | omega (Ry) | dF | rF | dInf | GF wall (s) |
|---:|---|---:|---:|---:|---:|---:|---:|
| 101 | [-1.320,1.247] | 2.567000e+01 | 0.00 | 5.863004e+03 | 2.922768e+00 | 4.124575e+03 | 4.768 |
| 101 | [-1.320,1.247] | 2.567000e+01 | 0.17 | 1.949771e+04 | 4.140564e+00 | 1.818680e+04 | 4.768 |
| 201 | [-1.320,1.247] | 1.283500e+01 | 0.00 | 1.131259e+03 | 5.639443e-01 | 6.665731e+02 | 10.185 |
| 201 | [-1.320,1.247] | 1.283500e+01 | 0.17 | 3.533753e+03 | 7.504332e-01 | 3.214330e+03 | 10.185 |
| 401 | [-1.320,1.247] | 6.417500e+00 | 0.00 | 5.218865e+02 | 2.601658e-01 | 2.362543e+02 | 19.700 |
| 401 | [-1.320,1.247] | 6.417500e+00 | 0.17 | 7.313458e+02 | 1.553097e-01 | 3.967645e+02 | 19.700 |
| 801 | [-1.320,1.247] | 3.208750e+00 | 0.00 | 2.839170e+02 | 1.415355e-01 | 1.439117e+02 | 40.620 |
| 801 | [-1.320,1.247] | 3.208750e+00 | 0.17 | 6.116730e+02 | 1.298958e-01 | 4.625690e+02 | 40.620 |
| 1601 | [-1.320,1.247] | 1.604375e+00 | 0.00 | 3.566748e+01 | 1.778061e-02 | 1.980038e+01 | 78.601 |
| 1601 | [-1.320,1.247] | 1.604375e+00 | 0.17 | 1.642182e+02 | 3.487363e-02 | 1.581836e+02 | 78.601 |
| 3201 | [-1.320,1.247] | 8.021875e-01 | 0.00 | 1.430396e+01 | 7.130669e-03 | 8.783328e+00 | 151.152 |
| 3201 | [-1.320,1.247] | 8.021875e-01 | 0.17 | 7.524164e+01 | 1.597843e-02 | 7.310850e+01 | 151.152 |
| 6401 | [-1.320,1.247] | 4.010938e-01 | 0.00 | 1.020087e+01 | 5.085237e-03 | 7.946389e+00 | 303.451 |
| 6401 | [-1.320,1.247] | 4.010938e-01 | 0.17 | 8.236275e+01 | 1.749068e-02 | 8.204445e+01 | 303.451 |
| 12801 | [-1.320,1.247] | 2.005469e-01 | 0.00 | 1.031215e+01 | 5.140714e-03 | 8.067484e+00 | 604.388 |
| 12801 | [-1.320,1.247] | 2.005469e-01 | 0.17 | 8.306496e+01 | 1.763980e-02 | 8.275069e+01 | 604.388 |

The finite-q follow-up was run once after the Gamma ladder stabilized:

| q / channel | N | interval (Ry) | h / eta_int | omega (Ry) | dF | rF | dInf | GF wall (s) |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| `(0.23,0,0)` / `chi_minus` | 6401 | [-1.320,1.247] | 4.010938e-01 | 0.00 | 1.020087e+01 | 5.085237e-03 | 7.946389e+00 | 152.170 |

## Integration-eta ladder

The physical response broadening remained fixed at `eta=0.040 Ry`. Every row
uses a grid with `h/integration_eta <= 0.5`; no zero-width extrapolation was
performed.

| integration_eta (Ry) | N | h / eta_int | omega (Ry) | dF | rF | dInf | GF wall (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0100 | 515 | 4.994163e-01 | 0.00 | 1.051986e+02 | 5.244260e-02 | 8.242566e+01 | 26.088 |
| 0.0100 | 515 | 4.994163e-01 | 0.17 | 7.124298e+02 | 1.512927e-01 | 7.085017e+02 | 26.088 |
| 0.0050 | 1029 | 4.994163e-01 | 0.00 | 5.088755e+01 | 2.536796e-02 | 4.002733e+01 | 50.355 |
| 0.0050 | 1029 | 4.994163e-01 | 0.17 | 3.841393e+02 | 8.157641e-02 | 3.825224e+02 | 50.355 |
| 0.0025 | 2055 | 4.999026e-01 | 0.00 | 2.751815e+01 | 1.371808e-02 | 2.120499e+01 | 97.324 |
| 0.0025 | 2055 | 4.999026e-01 | 0.17 | 2.060166e+02 | 4.375000e-02 | 2.050285e+02 | 97.324 |
| 0.0010 | 6401 | 4.010938e-01 | 0.00 | 1.020087e+01 | 5.085237e-03 | 7.946389e+00 | 303.451 |
| 0.0010 | 6401 | 4.010938e-01 | 0.17 | 8.236275e+01 | 1.749068e-02 | 8.204445e+01 | 303.451 |

Both frequencies move toward the Lehmann reference as the spectral width is
reduced. The approximately monotone trend is much larger than the residual
mesh/window changes at the final width, so the final finite-eta discrepancy
is not a stable GF-bubble offset.

## Energy-margin ladder

This ladder uses `integration_eta=0.001 Ry` and keeps the mesh spacing near
`4.01e-4 Ry`.

| margin (Ry) | N | interval (Ry) | h / eta_int | omega (Ry) | dF | rF | dInf | GF wall (s) |
|---:|---:|---|---:|---:|---:|---:|---:|---:|
| 0.60 | 6401 | [-1.320,1.247] | 4.010938e-01 | 0.00 | 1.020087e+01 | 5.085237e-03 | 7.946389e+00 | 303.451 |
| 0.60 | 6401 | [-1.320,1.247] | 4.010938e-01 | 0.17 | 8.236275e+01 | 1.749068e-02 | 8.204445e+01 | 303.451 |
| 1.00 | 8401 | [-1.720,1.647] | 4.008333e-01 | 0.00 | 9.860117e+00 | 4.915369e-03 | 7.827875e+00 | 420.112 |
| 1.00 | 8401 | [-1.720,1.647] | 4.008333e-01 | 0.17 | 8.210894e+01 | 1.743678e-02 | 8.183730e+01 | 420.112 |
| 2.00 | 13401 | [-2.720,2.647] | 4.005224e-01 | 0.00 | 1.004695e+01 | 5.008509e-03 | 7.878319e+00 | 601.367 |
| 2.00 | 13401 | [-2.720,2.647] | 4.005224e-01 | 0.17 | 8.228315e+01 | 1.747378e-02 | 8.197915e+01 | 601.367 |

Increasing the window changes the final static `rF` only from `5.09e-3` to
`5.01e-3` and the finite-frequency `rF` from `1.75e-2` to `1.75e-2`. The
window is therefore controlled and is not the dominant remaining difference.

## Timing and performance observation

The GF wall time is the evaluator-reported time for one compact GF call with
both Gamma frequencies. It scales approximately linearly with the number of
energy points:

| N | GF wall time (s) | time / N (s) |
|---:|---:|---:|
| 101 | 4.768 | 4.72e-02 |
| 201 | 10.185 | 5.07e-02 |
| 401 | 19.700 | 4.91e-02 |
| 801 | 40.620 | 5.07e-02 |
| 1601 | 78.601 | 4.91e-02 |
| 3201 | 151.152 | 4.72e-02 |
| 6401 | 303.451 | 4.74e-02 |
| 12801 | 604.388 | 4.72e-02 |

This is observation only. No GF-bubble optimization was made. The scaling
and the multi-minute resolved samples imply:

```text
GF BUBBLE PERFORMANCE REMEDIATION REQUIRED BEFORE FE
```

## Historical LR-GF-02 cross-check

`UnitLrGfSusceptibility` was rerun unchanged. Its existing simple finite-
spectrum oracle remained green:

| sample | absolute difference | reported relative difference |
|---|---:|---:|
| coarse: 401 points, `integration_eta=0.010`, margin `1.0` | 2.1443e-09 | 2.4633e-02 |
| fine: 2001 points, `integration_eta=0.001`, margin `1.0` | 3.6257e-09 | 4.1652e-02 |
| finite-spectrum analytic GF error | 1.1277e-09 | — |
| static q=0 GF difference | 9.2484e-10 | — |
| off-mesh-q diagnostic | 1.7422e-08 | 2.0014e-01 |
| alternate-eta diagnostic | 1.7456e-08 | 2.0017e-01 |

The test result was `UnitLrGfSusceptibility: PASS`; the coarse/fine relative
diagnostics are part of that unchanged test's existing evidence and were not
used as a new acceptance threshold.

## Interpretation and classification

1. The fixed-width Simpson ladder removes the unresolved 101/201/401-point
   behavior. By 6401–12801 points, unweighted spectral moments are in the
   `6e-4–9e-4` range and the GF values are stable to the finite-width/window
   envelope.
2. The Fermi-weighted moments decrease systematically as
   `integration_eta` is reduced. Their residual floor at `0.001 Ry` is the
   expected finite Lorentzian-width contribution on the prescribed finite
   window; it is independently visible before the Kubo bubble.
3. At fixed physical `eta`, both response frequencies approach the compact
   Lehmann result as `integration_eta` decreases. Once the mesh resolves the
   width, changing the window produces only a small change and does not leave
   a new stable GF-versus-Lehmann offset.
4. The requested aspirational `rF < 1e-5` is not reached at the smallest
   prescribed finite `integration_eta`; no arbitrary parameter tuning or
   extrapolation was applied. The tables show that the remaining discrepancy
   is dominated by the finite integration width, not unresolved Simpson
   sampling or window truncation.

Therefore the mechanical TDVK-03A result is:

```text
REAL-AXIS GF FORMULATION NUMERICALLY CONSISTENT
```

This is Case A. It is a controlled numerical-consistency result for the
nontrivial R3 finite-basis fixture only. No Fe run was performed, no SCF was
rerun for quadrature samples, and no TDVK-02 or TDVK-03 material status was
changed.

## Exact tests and commands

The audit harness was added as `UnitLrProductGfQuadratureAudit` and built with:

```text
cmake -S . -B build -DRUN_UNIT_TESTS=ON -DENABLE_OPENMP=ON \
  -DENABLE_MPI=OFF -DENABLE_LIBXC=ON
cmake --build build --target UnitLrProductGfQuadratureAudit -j2
```

The successful direct audit samples were run as:

```text
build/bin/UnitLrProductGfQuadratureAudit fixed_one 101
build/bin/UnitLrProductGfQuadratureAudit fixed_one 201
build/bin/UnitLrProductGfQuadratureAudit fixed_one 401
build/bin/UnitLrProductGfQuadratureAudit fixed_one 801
build/bin/UnitLrProductGfQuadratureAudit fixed_one 1601
build/bin/UnitLrProductGfQuadratureAudit fixed_one 3201
build/bin/UnitLrProductGfQuadratureAudit fixed_one 6401
build/bin/UnitLrProductGfQuadratureAudit extra                 # 12801
build/bin/UnitLrProductGfQuadratureAudit finite_q
build/bin/UnitLrProductGfQuadratureAudit eta_one 1             # 0.010 Ry
build/bin/UnitLrProductGfQuadratureAudit eta_one 2             # 0.005 Ry
build/bin/UnitLrProductGfQuadratureAudit eta_one 3             # 0.0025 Ry
build/bin/UnitLrProductGfQuadratureAudit margin_one 2          # 1.00 Ry
build/bin/UnitLrProductGfQuadratureAudit margin_one 3          # 2.00 Ry
```

The `integration_eta=0.001 Ry` / 6401-point eta-ladder row is the identical
successful `fixed_one 6401` sample, so it was not duplicated. The 0.60-Ry
window row is likewise the identical fixed 6401 sample.

The unchanged regression set was run with:

```text
ctest --test-dir build --output-on-failure -R \
  '^(UnitLrProductKsSusceptibility|UnitLrGfSusceptibility|UnitLrProductGfSusceptibility|UnitLrProductGfSusceptibilityRejectIntegrationEta)$'
```

Result: **4/4 passed**. This includes the historical `UnitLrGfSusceptibility`,
the unchanged R3 representation test at its original 21/41-point samples,
the compact Lehmann fixture, and the existing integration-eta guard.
