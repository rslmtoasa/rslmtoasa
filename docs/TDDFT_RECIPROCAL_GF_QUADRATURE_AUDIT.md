# TDVK-03A — reciprocal-GF real-axis quadrature closure audit

**Status: Historical TDVK-03A fixture audit; the TDVK-03R factorized backend
and Fe closure evidence are recorded below.**

The original TDVK-03A section is a fixture-level numerical audit. The later
TDVK-03R section records the separate performance remediation and the resumed
Fe backend closure; neither section adds physical interpretation.

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

## TDVK-03B — mixed-eigenvector reciprocal-GF numerical closure

**Status: PASS for the requested mixed-eigenvector numerical oracle; no
material-validation claim and no change to the TDVK-02/TDVK-03 material
status.** This is a small extension of the TDVK-03A harness. It changes only
the finite electronic fixture and leaves the radial/product fixture and all
response equations unchanged.

### Fixture and construction

The mixed run retains `nbasis=8`, `nbands=8`, `nk=2`, one site,
`response_lmax=2`, product dimension 27, and the TDVK-03A radial mesh,
nonconstant `phi`, nonconstant `phidot`, and circular product channels. For
each k point it constructs a deterministic dense complex Hermitian matrix,
then obtains eigenvalues/eigenvectors with the project LAPACK Hermitian
solver (`zheev`). No eigenvector was supplied analytically. The occupations
are generated with `lr_fermi_dirac_occupation`.

The resulting spectrum and density evidence were:

| quantity | value |
| --- | ---: |
| eigenvalue minimum / maximum | `-7.29939766e-01` / `7.57911598e-01` Ry |
| absolute Hamiltonian off-diagonal minimum / maximum | `1.78756818e-02` / `8.09891968e-02` |
| absolute imaginary off-diagonal minimum / maximum | `2.50000000e-03` / `1.58000000e-02` |

The maximum exact off-diagonal moment magnitudes, with the maximum imaginary
part also shown, were:

| p | `max|offdiag(M_p)|` | `max|offdiag(N_p)|` | `max|Im(M_p)|` | `max|Im(N_p)|` |
| ---: | ---: | ---: | ---: | ---: |
| 0 | `6.22483590e-16` | `1.64223182e-01` | `9.45424294e-17` | `1.24159062e-02` |
| 1 | `8.09891968e-02` | `2.84054116e-02` | `1.58000000e-02` | `6.86771733e-03` |
| 2 | `1.15642327e-01` | `1.97147502e-02` | `7.97040000e-03` | `4.66495775e-03` |

`M_0` is the identity by completeness, so its off-diagonal value is the
expected roundoff floor. `N_0`, `M_1`, `N_1`, `M_2`, and `N_2` contain the
intended complex off-diagonal structure.

### Spectral-moment ladder

The exact matrices were formed directly as
`sum_n epsilon_n**p * c_n * c_n^H` and, independently, with the explicit
Fermi factor for `N_p`. The numerical matrices were obtained only by
integrating the unchanged `build_weighted_resolvent` spectral function. All
three samples use `energy_margin=1.0 Ry` and satisfy `h/integration_eta <=
0.5`. The actual response integration interval is
`[-1.72993977, 1.75791160] Ry`.

| `integration_eta` (Ry) | `n` | `k` | `rM0 / rN0` | `rM1 / rN1` | `rM2 / rN2` |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.0100 | 801 | 1 | `3.9814e-03 / 6.9955e-02` | `4.3695e-03 / 1.0440e-02` | `4.4972e-03 / 8.5990e-03` |
| 0.0100 | 801 | 2 | `3.9661e-03 / 3.1142e-02` | `4.2154e-03 / 1.1120e-02` | `4.3556e-03 / 9.2760e-03` |
| 0.0050 | 1601 | 1 | `2.1541e-03 / 4.0613e-02` | `2.2334e-03 / 5.3055e-03` | `2.2724e-03 / 4.1921e-03` |
| 0.0050 | 1601 | 2 | `1.8440e-03 / 1.5909e-02` | `2.0588e-03 / 5.4556e-03` | `2.0996e-03 / 4.5939e-03` |
| 0.0025 | 3201 | 1 | `1.1909e-03 / 2.1048e-02` | `1.2676e-03 / 2.6573e-03` | `1.2843e-03 / 2.0933e-03` |
| 0.0025 | 3201 | 2 | `1.1104e-03 / 8.1836e-03` | `1.0585e-03 / 2.7200e-03` | `1.0093e-03 / 2.0980e-03` |

The three grids use `n=801,1601,3201`, respectively; each has
`h/integration_eta=0.43598142`. The residuals decrease with the integration
width as expected for the finite-width real-axis oracle.

### Compact GF versus compact Lehmann

At every width the existing compact GF and compact Lehmann production
evaluators were called for `chi_plus`, Gamma, and frequencies 0 and 0.17 Ry.
No `1e-5` acceptance target was imposed.

| `integration_eta` | `n` | frequency (Ry) | `dF` | `rF` | `dInf` |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.0100 | 801 | 0.00 | `1.18262574e+02` | `6.25572231e-02` | `8.49495870e+01` |
| 0.0100 | 801 | 0.17 | `7.57664350e+02` | `1.67776533e-01` | `6.53553625e+02` |
| 0.0050 | 1601 | 0.00 | `6.27848939e+01` | `3.32112558e-02` | `4.53038051e+01` |
| 0.0050 | 1601 | 0.17 | `4.20603809e+02` | `9.31381404e-02` | `3.62981733e+02` |
| 0.0025 | 3201 | 0.00 | `3.23556549e+01` | `1.71151350e-02` | `2.33071677e+01` |
| 0.0025 | 3201 | 0.17 | `2.22097968e+02` | `4.91811804e-02` | `1.91729926e+02` |

The GF/Lehmann discrepancy decreases monotonically at both frequencies as
`integration_eta` decreases. The required finite-q finest-width sample,
`chi_minus` at `q=(0.23,0,0)`, `omega=0`, `n=3201`, and margin 1.0 Ry, was
finite and gave `dF=3.23556549e+01`, `rF=1.71151350e-02`, and
`dInf=2.33071677e+01`.

### Basis-rotation invariant

The harness generates a dense deterministic unitary `Q` numerically with the
same LAPACK Hermitian eigensolver, forms `H'=Q^H H Q`, rotates the numerical
eigenvectors to `Q^H c_n`, and transforms every component vertex as
`V'=Q^H V Q`. The rotated Hamiltonian eigenpair residual was
`4.58803410e-16`.

The production compact API does not accept an externally supplied vertex
tensor. Therefore the rotation part stays in the test harness: it uses the
same component tensor, resolvent construction, two Kubo contractions, and
Lehmann pair sum locally; the unrotated local results were checked against
the production compact evaluators first. The maximum unrotated local-versus-
production differences were `0.00000000e+00` for GF and
`1.01684599e-12` for Lehmann.

At the resolved `n=801`, `integration_eta=0.010 Ry`, margin 1.0 Ry sample:

| frequency (Ry) | relative GF rotation difference | relative Lehmann rotation difference |
| ---: | ---: | ---: |
| 0.00 | `1.40494438e-15` | `1.40259015e-15` |
| 0.17 | `1.41060297e-15` | `1.24564983e-15` |

Both are below the required `1e-10` invariant threshold.

### Commands and test registration

The mixed audit was run with:

```text
cmake -S . -B build -DRUN_UNIT_TESTS=ON -DENABLE_OPENMP=ON \
  -DENABLE_MPI=OFF -DENABLE_LIBXC=ON
cmake --build build --target UnitLrProductGfQuadratureAudit -j2
build/bin/UnitLrProductGfQuadratureAudit mixed_one 1
build/bin/UnitLrProductGfQuadratureAudit mixed_one 2
build/bin/UnitLrProductGfQuadratureAudit mixed_one 3
build/bin/UnitLrProductGfQuadratureAudit mixed_q
build/bin/UnitLrProductGfQuadratureAudit mixed_rotation
```

`UnitLrProductGfMixedEigenvectors` is registered as a disabled CTest audit,
matching the existing TDVK-03A audit because the resolved GF runs are
multi-minute diagnostic samples. The mixed fixture therefore does not run in
the default CTest set. No `0.001 Ry` / 600-second regime was run, no
production physics or material input was changed, and no Fe/native-RSGF
claim is made.

## TDVK-03R — real-axis GF quadrature/contraction factorization

**Status: PASS for the blocker remediation; the parent TDVK-03 result is
`PASS CANDIDATE`.** This section records the implementation and numerical
evidence for the factorized reciprocal-GF backend. The scalar implementation
and the `38102e1` optimized dense-resolvent implementation remain compiled
and callable as correctness oracles.

### Algebraic reordering and independence

For each k point, the four existing component vertices are transformed as
`Vtilde(I,p,q,n,m)=<L n|V(I,p,q)|R m>`. The factorized transition amplitude is
then formed directly from the GF component-vertex path:

```text
T(I,n,m) = sum_(p,q=0,1) epsilon_L(n)^p epsilon_R(m)^q Vtilde(I,p,q,n,m)
```

The real-energy loop still explicitly performs the same Simpson integral. It
accumulates only the scalar band-pair kernel:

```text
K_nm(k,omega) = sum_E w_E f(E) 2 w_k/sum(w_k)
                [a_Ln(E; integration_eta) g_Rm^R(E+omega; eta)
                 + a_Rm(E; integration_eta) g_Ln^A(E-omega; eta)]
```

where `a=i(GR-GA)/(2*pi)` uses the finite integration broadening and the
shifted denominators retain the physical response broadening. After the
quadrature, the complete compact matrix is formed as
`chi(I,J)+=sum_nm K_nm T(I,n,m) T(J,n,m)*`. Thus the factorized backend
commutes only finite band sums and the compact contraction through the
existing real-axis integral; it does not substitute the Lehmann kernel or
call the Lehmann accumulator. The exact `component=1+p+2*q` convention,
both Kubo terms, circular channels, q endpoint, occupations, k weights,
factor of two, signs, conjugations, and complete 232-coordinate product
space are unchanged.

### Correctness hierarchy

The transition-factorization oracle on the mixed-complex fixture produced:

| channel | maximum absolute error | maximum relative error |
|---|---:|---:|
| `chi_plus` | `1.73046935e-15` | `2.13050107e-16` |
| `chi_minus` | `1.38624879e-16` | `1.53171174e-16` |

The compact unit fixture's maximum GF transition residual was `1.2064e-16`.
The complete compact matrix was then compared on identical states, q,
channel, frequency, physical eta, integration eta, energy window, and
Simpson grid. Pairwise values below are ordered `scalar/optimized`,
`scalar/factorized`, `optimized/factorized`.

| fixture | omega (Ry) | all three norms | pairwise `dF` | pairwise `rF` | pairwise `dInf` |
|---|---:|---:|---|---|---|
| mixed Gamma | `0.00`, `0.17` | `1.79506145e+03` | `1.1081e-14`, `5.0295e-12`, `5.0284e-12` | `6.1728e-18`, `2.8018e-15`, `2.8012e-15` | `7.3241e-15`, `3.8666e-12`, `3.8666e-12` |
| mixed q=`(0.23,0,0)` | `0.00` | `1.86403113e+03` | `6.9717e-14`, `1.2277e-11`, `1.2266e-11` | `3.7401e-17`, `6.5860e-15`, `6.5804e-15` | `6.2841e-14`, `9.0982e-12`, `9.0982e-12` |

All pairwise relative differences are below the existing `1e-10` response
regression threshold, including Gamma/static, finite frequency, finite q,
complex eigenvectors, and nonconstant endpoint-energy dependence.

### Performance evidence

Correctness was decided before timing. Timings below are in the same pair
order for speedup and in `scalar`, `optimized`, `factorized` order for wall
time:

| fixture | integration points | wall time (s) | speedup |
|---|---:|---|---|
| mixed Gamma, 0/0.17 Ry | `801` | `30.2637 / 4.22790 / 0.00296218` | `7.158 / 1.0217e4 / 1.4273e3` |
| mixed finite q, 0 Ry | `3201` | `61.6816 / 8.96611 / 0.00593637` | `6.879 / 1.0390e4 / 1.5104e3` |
| compact unit Gamma/static | `21` | `0.787328 / 0.112390 / 0.000502181` | `7.005 / 1.5678e3 / 2.2380e2` |
| compact unit finite q | `21` | `0.795392 / 0.112959 / 0.000492368` | `7.041 / 1.6154e3 / 2.2942e2` |

For the Fe accepted state, the factorized evaluator reports no dense GF matrix
allocation. Its major allocations are approximately 4.81 MB for the complete
component vertex tensor, 1.20 MB for one-k transition amplitudes, and 0.86 MB
for the compact 232-by-232 response.

### Resumed Fe TDVK-03 closure

The same production input and one accepted `8x8x8` Fe state were used for the
complete Gamma `chi_plus`, `omega=0`, `eta=0.04 Ry` ladder. The accepted state
was `18` basis states/bands, `EF=-8.51191027e-02 Ry`, `T=300 K`, moment
`2.267445 mu_B`, and SCF residual `6.617e-7`. The compact product dimension
was `232`; every sample reused this state and the exact Gamma endpoint.

| sample | integration eta (Ry) | N | margin (Ry) | h/integration eta | `||chi_L||_F` | `||chi_GF||_F` | dF | rF | dInf | wall (s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| eta ladder | `0.0100` | `953` | `0.60` | `0.399426` | `3.79995096` | `3.61768630` | `2.57226666e-01` | `6.76921014e-02` | `1.40114603e-01` | `7.67390` |
| eta ladder | `0.0050` | `1903` | `0.60` | `0.399846` | `3.79995096` | `3.70701147` | `1.28850597e-01` | `3.39084894e-02` | `7.03065297e-02` | `9.34962` |
| eta ladder | `0.0025` | `3805` | `0.60` | `0.399846` | `3.79995096` | `3.75267245` | `6.48187541e-02` | `1.70577870e-02` | `3.53119546e-02` | `12.7772` |
| eta ladder / Simpson base / window 0.60 | `0.0010` | `9509` | `0.60` | `0.399931` | `3.79995096` | `3.78102764` | `2.58386514e-02` | `6.79973285e-03` | `1.40767444e-02` | `23.0595 / 23.1370 / 23.0656` |
| Simpson fine | `0.0010` | `19017` | `0.60` | `0.199965` | `3.79995096` | `3.78098824` | `2.58942459e-02` | `6.81436317e-03` | `1.40921415e-02` | `40.3156` |
| window 1.00 | `0.0010` | `11509` | `1.00` | `0.399943` | `3.79995096` | `3.78117876` | `2.57510806e-02` | `6.77668760e-03` | `1.40222232e-02` | `26.7371` |
| window 2.00 | `0.0010` | `16509` | `2.00` | `0.399960` | `3.79995096` | `3.78130496` | `2.56779024e-02` | `6.75742994e-03` | `1.39758163e-02` | `35.8025` |

The GF/Lehmann discrepancy decreases monotonically across the integration_eta
ladder. The base/fine Simpson change and the margin changes are small relative
to the finite-width trend. The first controlled Fe sample completed in
`7.67390 s`; all later prescribed controls completed, so the original material
cost gate is open. No final material GF↔Lehmann tolerance was invented.

### Commands

```text
cmake --build build -j2 --target UnitLrProductGfSusceptibility \
  UnitLrProductGfQuadratureAudit UnitTddftProductionDriver rslmto.x
ctest --test-dir build --output-on-failure -R \
  '^(UnitLrLmtoProductResponseBasis|UnitLrLmtoProductResponse|UnitLrLmtoProductStrictRankGuard|UnitLrProductKsSusceptibility|UnitLrGfSusceptibility|UnitLrProductGfSusceptibility|UnitLrProductGfSusceptibilityRejectIntegrationEta|UnitTddftProductionDriver)$'
```

Result: **8/8 focused tests passed**. The disabled long-form mixed audit was
also run directly for `mixed_one 1` and `mixed_q`; the Fe production closure
run completed with no GF errors or non-finite compact response.
