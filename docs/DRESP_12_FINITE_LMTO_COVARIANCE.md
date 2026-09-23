# DRESP-12 — finite-LMTO covariance closure

Starting HEAD for this closure: `28c8058091bbd86a2db53e071861fa190fa8c49a`.

DRESP-12 is a static, diagnostic-only Gamma-point closure for the accepted
64-k bcc-Fe `ham_only` state. It does not change the response kernel,
denominator, SVD, production dynamics, Goldstone correction, or XC
formulation. DRESP-13, LMTO basis-response/Sternheimer work, finite-q, and
finite-frequency dynamics are outside this milestone.

## Covariant endpoint closure

The six endpoint branches remain `00, 10, 01, 11, 20, 02`. The complete
endpoint tangent is split as

```text
delta D = delta D^(rho) + delta D^(H)
delta D^(rho) = H^p delta_rho H^q
delta D^(H)   = delta D^(complete) - delta D^(rho)
```

The stale-density reuse was fixed. Both `endpoint_complete` and
`endpoint_fixed` now use

```text
deltaH = dh_cov
delta_rho = delta_rho_cov
```

The corrected complete endpoint agrees directly with the authoritative
covariant endpoint object for all six branches; the measured maximum branch
residual is zero in the accepted Fe run. The decomposition and the
`delta_rho=0` endpoint-H isolation also close at approximately `1.3e-16`.

The obsolete endpoint-H Bxc/connection split was removed from the primary
classification. It was channel-dependent because its one-sided circular
source was being compared with a Cartesian observable.

## Circular/Cartesian measurement seam

The independent circular paths retain the production convention

```text
m_x = m_plus + m_minus
m_y = i*(m_plus - m_minus)
```

The frozen DRESP-09Y path and DRESP-12 endpoint measurement use the same
physical Pauli six-branch dual. The accepted Fe values are:

```text
DRESP12 vs DRESP09Y endpoint       5.8774991076e-15
trace oracle vs DRESP09Y           2.6729845264e-16
trace oracle vs DRESP12             5.8433554802e-15
circular Cartesian reconstruction   5.7985501685e-15
Pauli seam residual                 5.2997972883e-15
```

The seam classification is `MEASUREMENT_CONVENTIONS_IDENTICAL`. Circular and
Cartesian paths are independent equivalence checks; they are not reported as
duplicate physics metrics.

## Angular resolution and metric

For a response vector `v`, `norm_by_l` accumulates

```text
||P_L v||^2 = sum_flat(response_l=L) metric_weights(flat)*|v(flat)|^2
```

The metric is the existing radial quadrature metric, with the normalized
response harmonics providing the angular orthogonality. Therefore

```text
||v||^2 = sum_L ||P_L v||^2
```

The accepted Fe results are:

| vector | L=0 | L=1 | L=2 | L=3 | L=4 | full |
|---|---:|---:|---:|---:|---:|---:|
| `r_fixed = delta m_B - P3` | 6.8365799764e-2 | 2.55e-17 | 8.65e-17 | 2.77e-17 | 1.4521659355e-1 | 1.6050464672e-1 |
| `master` | 9.8209875599e-15 | 2.81e-17 | 1.05e-16 | 2.81e-17 | 1.5208439341e-1 | 1.5208439341e-1 |
| live `D*m_G` | 6.8365802023e-2 | 2.53e-17 | 8.66e-17 | 2.73e-17 | 1.4521659430e-1 | 1.6050464837e-1 |

The norm reconstruction residuals are `8.08e-16` for `r_fixed`, zero for
`master`, and `1.35e-16` for `D*m_G`. The nonspherical master norm is
`1.5208439341e-1`; its L=4 fraction is `1.0` within the reported precision.
The total nonspherical remainder relative to `P3` is

```text
R_>0 = ||(1-P0)master|| / ||P3|| = 2.2390504634e-1.
```

The actual spherical accounting residual is

```text
R_0 = ||P0 master|| / ||P0 P3|| = 1.4458871324e-14.
```

The historical fixed-basis defect remains
`fixed_basis_goldstone_relative = 2.3630169774e-1`; its L=0 component is
`6.8365799764e-2` and its dominant L=4 component is `1.4521659355e-1`.

## DRESP-11 action

No dense denominator/SVD campaign is rerun. DRESP-12 applies the existing
matrix-free denominator action once to the current in-memory `m_G`, maps the
result back to the physical response space, and resolves that vector by L.
The reported `D*m_G` values above are therefore current-state values, not
numbers read from a generated DRESP-11 artifact. The historical orthogonal
fraction is `0.9065827902`.

## Model boundary

The strict ASA transverse model defines the ground-state XC functional and
potential from the spherical radial spin density. Its Goldstone consistency
condition belongs to the same spherical variational space:

```text
P0 master = 0
```

The full-spatial response machinery can still generate finite L>0 band-density
responses on that spherical ASA ground state. Those are meaningful response
coordinates, but the ASA ground-state functional has no corresponding
nonspherical ground-state XC field. The measured L=4 remainder is therefore

```text
NONSPHERICAL_RESPONSE_ON_SPHERICAL_ASA_GROUND_STATE
```

It is not automatically basis incompleteness, numerical error, or Goldstone
failure. A fully self-consistent full-spatial Goldstone theory would require
a ground-state functional stationary in the nonspherical channels; that is
not implemented here.

The historical `Kxc = Bxc/P3` construction is retained only as a
representation diagnostic. The accepted live `VXC0SP` audit identifies the
XC-conjugate spherical density as `m_xc = n_up - n_down`, with
`B_xc = 0.5*(vxc_up-vxc_down)` and the external constraining field excluded.
For the accepted Barth-Hedin Fe state, `constraining_field_ry = 0`, the
weighted `m_xc` versus `P3` profile difference is `1.6000607535e-2`, and the
integrated values are `1.9989134964` and `1.9984971714`, respectively. The
strict-ASA LDA candidate is therefore `Kxc^ASA = Bxc/m_xc`; GGA requires a
separate gradient-dependent/noncollinear derivation. See
[`TDDFT_FORMULATION.md`](TDDFT_FORMULATION.md) for the production decision.

## Final closure result

```text
DRESP-12 verdict: PASS-A
Primary classification: ASA_L0_RIGID_RESPONSE_CLOSED
Goldstone correction: OFF
Dynamics: NOT RUN
Ward-focused campaign: CLOSED
NEXT: FORMULATION_DECISION
```

This closes the Ward-focused campaign in the strict ASA L=0 response space.
The finite L>0 response remains recorded as a model-space boundary result.
Generated artifacts, including `/tmp` products, are validation evidence only
and are not committed.
