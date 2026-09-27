# VAL-09 — Kubo-Bastin charge, spin, and orbital transport

## Scope and campaign

The validation uses the existing SOC fcc-Pt fixture at
`tests/postproc/cases/conductivity/fccPt`. The opt-in campaign is registered
as `Val09KuboBastinTransport` with `validation` labels and is intentionally
outside the quick developer gate:

```text
python3 tests/validation/val09_kubo_bastin_transport.py \
  --binary build-rf-serial/bin/rslmto.x \
  --scratch-root build-rf-serial/Testing/validation/val09_kubo_bastin_transport
```

The reduced campaign keeps Pt SOC, fcc symmetry, PBC, and the production
real-space hopping/current path, but uses `rc=20`, replication `n=4,6,8`,
120–480 energy channels, and observable-specific Chebyshev-order ranges. The
`n=2` point was probed and exits before producing a transport file with this
reduced cutoff, so it is outside the valid neighbor-topology envelope and is
not treated as a physical zero.

## Operator and unit audit

The production operators are:

| Observable | Production definition in this campaign |
|---|---|
| charge | \(v_{a/b}(R)=(1/i)[\hat d_{a/b}\cdot(r_i-r_j)\,a_{lat}]H_{ij}(R)\); `velocity_scale` can multiply \(v_b\) |
| spin | \(J^z_a=\frac12\{S_z,v_a\}\), from `js_alpha='z'` |
| orbital | \(J^L_a=\frac12\{L_z,v_a\}\), from `jl_alpha='z'` |

The first tensor direction is `v_alpha` and the second is `v_beta`; reversing
either direction must reverse the reported scalar. The implementation uses a
fixed Lorentz kernel, `lorentz_kernel(cond_ll, ..., 6.0d0)`, in
`source/conductivity.f90`. No user-facing eta or kernel-alpha parameter is
currently exposed. The campaign therefore varies `cond_ll` and the
Chebyshev energy window as the available resolution/broadening controls and
records the fixed kernel explicitly.

The estimator reports the internal Ry/alat convention with
`16/(pi*DeltaE**2)`. The `e^2/hbar` and volume conversion factors visible as
commented alternatives in `conductivity.f90` are not applied. Consequently,
the Pt values below are not numerically compared with literature values in
`S/cm` or `hbar/e S/cm`; doing so would be a unit error. The material anchor
is used for qualitative tensor structure and internal order-of-magnitude
checks only.

## Numerical evidence

All values are the real part of the row nearest \(E-E_F=0\) in the production
`Pt_cond.out` file. The common tensor probe used `cond_ll=20`, `n=4`, and the
default Pt window.

| component | value |
|---|---:|
| charge \(xx\) | 2.257744 |
| charge \(yy\), \(zz\) | 2.257744, 2.257744 |
| charge \(xy\), \(yx\) | 6.969961e-3, 6.969961e-3 |
| spin \(z;xy\), \(z;yx\) | -2.742699e-2, 2.742699e-2 |
| spin \(z;zz\) | 1.746255e-11 |
| orbital \(L_z;xy\), \(L_z;yx\) | 3.785853e-1, -3.785853e-1 |
| orbital \(L_z;zz\) | -1.103463e-12 |

The charge diagonal spread is exactly zero at the printed precision; the
forbidden charge Hall components are \(3.09\times10^{-3}\) of the diagonal
scale. The transverse spin/orbital components have the expected antisymmetric
sign relation for this fixed polarization, while the \(z;zz\) components are
numerically zero. The small spin/orbital \(xx\) and \(yy\) values are retained
in the JSON report and are not clamped.

Convergence envelopes from the campaign are:

| axis | charge | spin | orbital |
|---|---:|---:|---:|
| final order step / observed scale | 3.92e-2 | 3.59e-1 | 2.00e-2 |
| final replication step / observed scale | 2.04e-4 | 1.12e-3 | 2.86e-4 |
| Fermi offsets tested | -0.15, 0, +0.15 Ry | -0.15, 0, +0.15 Ry | -0.15, 0, +0.15 Ry |
| window endpoints tested | (-2.0,1.0), (-2.5,1.2), (-3.0,1.5) Ry | same | same |

The order sequences were charge `60,80,100` at `n=4`, spin `40,60,80` at
`n=8`, and orbital `120,160,200` at `n=8`; the last-step acceptance envelope
in the campaign is 0.40 of the observed scale. The Fermi/window variation is
reported as physical spectral sensitivity, not forced into a false agreement.

## Exact/eigenpair cross-representation check

The exact reciprocal moment implementation is explicitly limited to
`cond_type='charge'` and `cond_calctype='per_type'`. It is therefore compared
only with charge recursion, using the existing small bcc-Fe conductivity
fixture and the same charge operator. At \(E_F\), `cond_ll=30`, and a 6³
reciprocal mesh:

| route | \(\sigma_{xx}\) |
|---|---:|
| real-space recursion | 3.237591 |
| exact eigenpair moments | 3.234408 |

The absolute difference is 0.003183, or \(9.83\times10^{-4}\) relative to the
larger value. Spin/orbital Pt results are not compared to the exact route:
that route does not implement those current operators, and such a comparison
would violate the same-observable requirement.

## Checklist and maturity

- [x] charge transport convergence established within the stated finite-window, fixed-kernel Pt scope
- [x] spin transport convergence established within the stated finite-window, fixed-kernel Pt scope
- [x] orbital transport convergence established within the stated finite-window, fixed-kernel Pt scope
- [x] tensor symmetries checked
- [x] forbidden components checked
- [x] current-operator conventions documented
- [x] exact/eigenpair comparison made where genuinely equivalent
- [x] material anchor assessed with unit awareness
- [x] expensive sweeps kept outside quick CI
- [x] maturity updated separately for charge/spin/orbital

The resulting maturity is scoped, not a universal conductivity claim:

- charge: **Validated (scoped)** for the fcc-Pt charge operator and the
  existing charge exact-moment route envelope;
- spin: **Validated (scoped)** for the SOC fcc-Pt \(S_z\)-anticommutator
  current, fixed Lorentz kernel, and finite-window/replication campaign;
- orbital: **Validated (scoped)** for the SOC fcc-Pt \(L_z\)-anticommutator
  current, fixed Lorentz kernel, and finite-window/replication campaign.

Remaining limitations are the hard-coded kernel parameter, incomplete external
unit conversion, the reduced Pt cutoff envelope, and the absence of an exact
spin/orbital eigenpair implementation. No literature normalization was tuned.
