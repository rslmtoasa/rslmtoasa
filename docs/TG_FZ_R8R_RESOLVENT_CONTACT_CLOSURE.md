# TG-FZ-R8R — resolvent contact trace and fixed-z closure

Status: **PASS-A — fixed-z Hamiltonian, log-det, and common-alpha Turek
routes close.**

This repair remains restricted to the zero-
`q` resolvent audit. It does not enter contour integration, native `J_ij`,
native `J(q)`, Fe-shell comparison, or DRESP-04.

## Repair

`source/lr_kl_hessian.f90` had evaluated the two-factor call
`trace_product(mixed,green)` as `Tr(mixed)`. The repaired helper evaluates

```text
trace_product(A,B)       = Tr(A B)
trace_product(A,B,C,D)   = Tr(A B C D)
```

and rejects mismatched optional arguments and matrix shapes. The sign and
the public `force_theorem_integrand` convention were not changed:

```text
K = -Im Tr[T_i G T_j G + H_ij G] / pi.
```

The independent 2x2 public-API regression uses zero torque matrices and
reports:

| quantity | value |
|---|---:|
| expected `-Im Tr(mixed*green)/pi` | `-7.48028233e-01` |
| public API contact | `-7.48028233e-01` |
| erroneous `-Im Tr(mixed)/pi` | `-7.95774715e-01` |

Thus the regression fails with the pre-R8R helper.

## Scalar sign derivation

For

\[
 F(z,\theta_i,\theta_j)=-\frac{1}{\pi}\Im\log\det[zI-H(\theta_i,\theta_j)],
 \qquad G=(zI-H)^{-1},
\]

the determinant derivative is

\[
 \partial_i\log\det(z-H)
 =\operatorname{Tr}\{G(-T_i)\}
 =-\operatorname{Tr}(GT_i).
\]

Since \(\partial_jG=GT_jG\),

\[
 \partial_j\partial_i\log\det(z-H)
 =-\operatorname{Tr}(GT_jGT_i+GH_{ij}).
\]

Therefore

\[
 \boxed{F_{ij}=+\frac{1}{\pi}\Im\operatorname{Tr}
 [T_iGT_jG+H_{ij}G]},
 \qquad
 \boxed{K_{\rm complete}=-F_{ij}}.
\]

The observed sign relation is consequently algebraic, not fitted.

## Common-alpha determinant oracle

For the common spin-independent alpha representation, the independent
full-spinor construction used

\[
 P'_i=\frac{\Delta P_i}{2}\sigma_x,
 \qquad
 F^{\alpha,\mathrm{analytic}}_{ij}
 =\frac{1}{\pi}\Im\operatorname{Tr}
 [g^\alpha P'_j g^\alpha P'_i].
\]

The explicit spin trace gives

\[
 F^{\alpha,\mathrm{analytic}}_{ij}=J_{ud}+J_{du},
 \qquad
 K_{\rm complete}=-(J_{ud}+J_{du})=-2J_{\rm sym}.
\]

The maximum independent full-matrix residual was `5.551115e-17`.

## Corrected fixed-z ledger

`TT`, `contact`, and `complete` below are the exact-Hamiltonian values from
the repaired public API. `F_logdet` is the best value in the epsilon sweep;
its residual is against `-K`, as required by the sign derivation.

| z | TT_exact | contact_exact | complete_exact | J_ud | J_du | J_ud+J_du | -K | F_logdet | H-vs-Turek | H-vs-logdet |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -0.91+0.83i | -8.55241423e-05 | 8.60942911e-05 | 5.70148761e-07 | -2.85074381e-07 | -2.85074381e-07 | -5.70148761e-07 | -5.70148761e-07 | -5.70159583e-07 | roundoff | 1.08e-11 |
| -0.17+0.04i | 3.35908538e-01 | 2.42741913e-03 | 3.38335957e-01 | -1.69167978e-01 | -1.69167978e-01 | -3.38335957e-01 | -3.38335957e-01 | -3.38336014e-01 | 1.67e-16 | 5.66e-08 |
| 0.62+0.31i | 8.03438078e-05 | 1.38905875e-04 | 2.19249682e-04 | -1.09624841e-04 | -1.09624841e-04 | -2.19249682e-04 | -2.19249682e-04 | -2.19250163e-04 | roundoff | 4.80e-10 |

The maximum Hamiltonian-vs-Turek residual was `4.996004e-16`. The maximum
best-window Hamiltonian-vs-log-det residual was `5.658796e-08`.

The independently constructed direct-matmul TT, contact, and complete trace
oracles agreed with `force_theorem_integrand` at `0.000000e+00` printed
residual for every fixed-z point and for normalized-gamma, exact, and H2
Hamiltonians.

## Log-det epsilon sweep and Hermiticity

The sweep used

```text
1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5.
```

At `z=-0.17+0.04i`, the residuals at `3e-4` and `1e-4` were
`7.70e-08` and `5.66e-08`; the `3e-5` and `1e-5` values worsened to
`3.28e-06` and `2.46e-05`, respectively. This is the expected numerical
window/cancellation behavior of the fourth-order mixed finite difference,
not a physics mismatch.

Before every `ZHEEV` log-det evaluation, the audit checked Hermiticity without
symmetrization:

| matrix set | max `||H-H^dagger||_max` |
|---|---:|
| base H_gamma | `1.873501e-16` |
| base H_exact | `1.942890e-16` |
| all rotated +/- stencil states | `2.012279e-16` |

## Corrected H2 truncation

The repaired contact contraction gives:

| quantity | corrected maximum |
|---|---:|
| complete absolute error `max |K_H2-K_exact|` | `5.887635e-03` |
| complete relative error | `9.772465e+01` |
| `||H2-Hexact||_max` | `4.070152e-01` |
| `||T2-Texact||_max` | `1.165566e-02` |
| `||C2-Cexact||_max` | `4.153483e-05` |

The per-z H2 values are printed by the focused test; the corrected maximum is
not the pre-R8R `6.320777e-03` value.

## Blast-radius audit

Current-branch callers found by source/test search:

| API | callers | disposition |
|---|---|---|
| `force_theorem_integrand` | R4 diagnostic code and the R8 fixed-z test | zero-q resolvent contact path affected; regenerated in R8R |
| `force_theorem_hessian_from_green` | public definition only; no current caller found | affected if invoked with a zero-q resolvent contact; no caller output to regenerate |
| spectral eigenbasis contact path | eigenbasis implementation/tests | unaffected by this helper bug |
| `lr_kl_contour` finite-q contact solve | no call to the helper | unaffected |
| mixed-derivative matrix construction | independent assembly routines | unaffected |

The R4 scalar ledger was regenerated from the direct `Tr(H_ij G)` oracle;
historical pre-R8R API contact/complete printouts are not retained as valid
results. R1/R2/R4 matrix and vertex identities are unchanged.

## Earlier campaign disposition

Confirmed by the rerun: DRESP-03T analytic torque derivative, DRESP-03T
mixed derivative, spectral TT/contact/complete finite-difference closure,
R7 representation normalization and alpha authority, R8 finite-angle
gamma/exact covariance, R8 first/mixed derivative covariance, and R8 uniform
rotation invariance.

Regenerated: every fixed-z scalar contact and complete value produced through
`force_theorem_integrand` in the affected R4/R8 diagnostic paths. No old value
was reinterpreted as correct after the repair.

## Verification and frozen invariant

The focused R8R test, the requested DRESP/LMTO selections, and the build were
run with:

```text
cmake --build build -j2
ctest -V -R '^UnitDresp03tgNativeFixedZ$' --output-on-failure
ctest -V -R '^UnitDresp03t' --output-on-failure
ctest -R 'UnitDresp03|UnitLrKlStaticBridge' --output-on-failure
```

`source/exchange.f90` was not modified. Its SHA-256 remains

```text
6e9ae7da6af46367aaf18a4fef37bdd4a9a0d6f69fe8e08b589e5c5f278a70e6
```

No contour, native `J`, native `J(q)`, Fe-shell, or DRESP-04 work was started.
