# TDVK-02R0 — LMTO product response-basis closure audit

**Verdict: PASS for algebraic representation closure.** The finite LMTO
radial-product space spans the current complete LR-05 point-grid transition
vectors on the affordable finite fixture to below the task blocking threshold
of `1e-10`. This is representation evidence only; it is not Fe material
validation and it does not change the TDVK-02 material verdict.

## Scope and definitions

The audit uses the live `lmto_radial_basis`, LR-04 metric, and
`evaluate_pauli_transition_vertex`. It does not modify LR-05, LR-GF-02, LR-06,
KXC, Dyson, the driver, eta, signs, circular normalization, or the material
response cutoff.

For each supported site, response harmonic, ordered pair, and endpoint power,
the audited functions are

```text
u^(0)(r) = phi_large(r) - enu_work * phidot_large(r)
u^(1)(r) = phidot_large(r)
B^(p,q)_(a,L,M,l,l')(r) = u^(p)_(a,l,sigma_L)(r)
                           * u^(q)_(a,l',sigma_R)(r) / r^2
```

`chi_plus` uses left spin 1/right spin 2 and `chi_minus` uses the reversed
ordering, as encoded by the existing Pauli matrices. The origin uses the
existing LR-05/LR-GF-02 rule: non-`s-s` products are zero at the origin and
the `s-s` value is extrapolated from the first two positive mesh points. The
metric is the existing LR-04 radial Simpson/volume metric; no second
quadrature or radial weight is introduced.

The fixture has one site, one explicit circular channel per run, `nr=51`,
`a=0.03`, `b=0.10`, and `lmax=1` or `2`. Radial channels come through the
existing `legacy_radial_fixture` and are captured by `lmto_radial_basis`; the
point-grid transition vectors come from the live LR-05 evaluator. The
closure states are four deterministic Hermitian-fixture eigenpairs:
`(left,right)=(1,2),(3,8),(5,12),(9,18)`, with both circular channels.

## Candidate inventory

Only ordered `(l,l')` pairs satisfying the exact triangle and parity rules are
included. Each allowed pair contributes four `(p,q)` candidates.

| basis | `L` | allowed ordered pairs | block dimension per `(L,M,channel)` | all-`M` count |
|---|---:|---:|---:|---:|
| sp | 0 | 2 | 8 | 8 |
| sp | 1 | 2 | 8 | 24 |
| sp | 2 | 1 | 4 | 20 |
| **sp total** |  |  |  | **52** |
| spd | 0 | 3 | 12 | 12 |
| spd | 1 | 4 | 16 | 48 |
| spd | 2 | 4 | 16 | 80 |
| spd | 3 | 2 | 8 | 56 |
| spd | 4 | 1 | 4 | 36 |
| **spd total** |  |  |  | **232** |

The required `spd`, `response_lmax=4` count is therefore exactly 232, not
greater than the contractual upper bound. Independent site, `(L,M)`, and
circular-channel blocks are used.

For the `nr=51` audit fixture, the complete point-vector dimensions are 459
(`sp`) and 1,275 (`spd`). A stored complex-kind-8 candidate-vector matrix
would be approximately 0.36 MiB (`459*52`) and 4.51 MiB (`1275*232`), while
the corresponding dense point-grid response matrices are approximately 3.22
MiB and 24.8 MiB. These are audit-fixture estimates, not production
allocations.

For the tracked Fe radial mesh (`nr=495`), `spd` with `response_lmax=4` has
`25*495 = 12,375` point coordinates. The dense complex point-grid matrix is
`12,375^2 = 153,140,625` elements, about 2.45 GB decimal (2.28 GiB). A
232-coordinate dense product coefficient matrix is about 0.86 MB decimal
(0.82 MiB), a factor of 53.34 fewer coordinates and a 98.12% coordinate
reduction. Storing all Fe point samples of the 232 product vectors would be
about 45.9 MB decimal; this is still only an audit representation. No
product-basis LR-06 production path is introduced here.

## Gram spectra

The test evaluates `S = B^H W B` for every block: 18 `sp` blocks and 50 `spd`
blocks, covering both circular orderings and every `M`. The printed spectra
below are grouped because the spectra are identical for every listed `M` and
the two channel orderings on this one-site fixture. The full unscaled
eigenvalue lists are recorded, including roundoff-sign values.

The diagnostic rank rule is
`lambda > N*epsilon*lambda_max`. It is reported only to expose numerical
rank/nullity; it is not a production truncation or response cutoff.

### sp blocks

| `L` | `M` values | dim | largest | smallest resolved | condition estimate | rank / nullity | diagnostic rank tolerance |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 8 | `1.490890e+02` | `1.802763e-10` | `8.2700e+11` | 6 / 2 | `2.6484e-13` |
| 1 | -1, 0, 1 | 8 | `1.428714e+02` | `5.004501e-09` | `2.8549e+10` | 4 / 4 | `2.5379e-13` |
| 2 | -2, -1, 0, 1, 2 | 4 | `8.221110e+01` | `2.913548e-09` | `2.8217e+10` | 3 / 1 | `7.3018e-14` |

```text
sp L=0: -1.279298e-14 -3.128887e-20  1.802763e-10  5.367946e-10
         2.596606e-05  9.094659e-04  5.908827e+00  1.490890e+02
sp L=1: -1.370138e-14 -7.195552e-15 -1.923661e-17  1.785525e-17
         5.004501e-09  8.709036e-04  1.319004e-03  1.428714e+02
sp L=2:  2.616787e-15  2.913548e-09  8.401334e-04  8.221110e+01
```

### spd blocks

| `L` | `M` values | dim | largest | smallest resolved | condition estimate | rank / nullity | diagnostic rank tolerance |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 12 | `9.524747e+43` | `9.524747e+43` | `1.0000e+00` | 1 / 11 | `2.5379e+29` |
| 1 | -1, 0, 1 | 16 | `1.360748e+23` | `6.614565e+17` | `2.0572e+05` | 2 / 14 | `4.8343e+08` |
| 2 | -2, -1, 0, 1, 2 | 16 | `9.524747e+43` | `9.524747e+43` | `1.0000e+00` | 1 / 15 | `3.3839e+29` |
| 3 | -3, -2, -1, 0, 1, 2, 3 | 8 | `1.360748e+23` | `6.614565e+17` | `2.0572e+05` | 2 / 6 | `2.4172e+08` |
| 4 | -4, -3, -2, -1, 0, 1, 2, 3, 4 | 4 | `9.524747e+43` | `9.524747e+43` | `1.0000e+00` | 1 / 3 | `8.4597e+28` |

```text
spd L=0: -3.126943e+27 -1.117349e+26 -1.660250e-22  3.931279e-16
          4.998655e-11  1.944347e-10  2.364057e-05  5.943375e-04
          5.439822e+00  6.592005e+01  1.493753e+25  9.524747e+43
spd L=1: -2.915686e+07 -3.181306e+06 -1.185926e+04 -4.385330e+00
          -2.498328e-15 -1.261920e-18 -8.638253e-20  1.322893e-18
          1.161150e-11  7.050865e-08  1.035335e-03  4.127264e+01
          1.190582e+04  2.900595e+06  6.614565e+17  1.360748e+23
spd L=2: -3.126943e+27 -1.117349e+26 -1.978358e+07 -1.068119e+06
          -2.727805e+01 -9.753205e-12  1.045782e-09  9.950972e-07
          1.628097e+00  2.496838e+01  2.969240e+05  1.138549e+06
          3.036397e+17  8.629620e+22  1.493753e+25  9.524747e+43
spd L=3: -2.915686e+07 -3.181306e+06 -1.187033e+04 -2.126445e+00
          1.189439e+04  2.900595e+06  6.614565e+17  1.360748e+23
spd L=4: -3.126943e+27 -1.117349e+26  1.493753e+25  9.524747e+43
```

The `spd` legacy fixture has a very large radial dynamic range. Consequently,
the unscaled finite-precision eigensolver reports large roundoff-sign values
in formally positive-semidefinite Gram matrices. Those values are recorded,
not clipped or converted into a production rule. For the closure solve only,
candidate columns are diagonally normalized before applying the same
machine-precision diagnostic rank test; this is numerical scaling of the
diagnostic pseudoinverse and does not change the span.

## LR-05 span closure

The live LR-05 point vector was evaluated for four left/right fixture
eigenpair combinations and both circular channels. Each vector was projected
with the LR-04 metric, reconstructed on the original point mesh, and tested
with the requested relative residual.

| fixture pair | channel 1 (`chi_plus`) | channel 2 (`chi_minus`) |
|---:|---:|---:|
| (1, 2) | `5.8456e-12` | `6.0026e-12` |
| (3, 8) | `6.7274e-12` | `1.6517e-11` |
| (5, 12) | `1.1886e-11` | `8.7769e-12` |
| (9, 18) | `8.5947e-12` | `9.2028e-12` |
| **maximum** |  | **`1.6517e-11`** |

The maximum is below `1e-10`, so the closure requirement is PASS. No LR-05
formula or transition-vector implementation was changed.

## Energy-affine endpoint oracle

For every `l=0,1,2`, both spins, all 51 mesh points, and endpoint energies
`[-0.37,-0.11,0.19,0.43]`, the existing LR-05 endpoint expression was
compared with `u^(0)+E*u^(1)`:

```text
maximum absolute difference = 3.0518e-05
maximum relative difference = 2.2191e-16
```

The absolute value reflects the large-magnitude `spd` fixture channels; the
relative result is at machine precision.

## LR-GF-02 consistency

The audit mechanically checked all four components, all supported `l,l'`,
both spins, and all radial points. The component map is the live convention
from `source/lr_gf_susceptibility.f90`:

```text
component = 1 + p + 2*q
```

The radial factor is the same `u^(p)*u^(q)/r²` expression and the same origin
extrapolation used by LR-05. Result:

```text
components checked = 4
maximum radial-factor difference = 0.0000e+00
```

Exact code reuse is not currently possible because the LR-GF-02 helper
`radial_vertex_component` and its `radial_component` helper are private to
`lr_gf_susceptibility_mod`, while LR-05 owns a separate private
`radial_product_at_point`. The audit therefore documents a small duplicated
formula seam. A later implementation task may centralize this algebra, but
this slice does not alter either production module.

## Reproducibility and checklist

The audit executable is `UnitLrLmtoProductResponseBasis` and is registered as
a standalone unit test. The focused commands were:

```text
cmake -S . -B build
cmake --build build --target UnitLrLmtoProductResponseBasis -j2
ctest --test-dir build -R '^UnitLrLmtoProductResponseBasis$' --output-on-failure
```

- [x] Candidate inventory and exact `sp`/`spd` counts recorded.
- [x] Every `sp`/`spd`, `L`, `M`, and circular block Gram spectrum computed.
- [x] Diagnostic rank, condition, and nullity recorded without a production cutoff.
- [x] Existing LR-05 point vectors reconstructed with the LR-04 metric.
- [x] Both circular channels and four fixture eigenpair combinations checked.
- [x] Energy-affine identity checked for all supported `l` and spins.
- [x] LR-GF-02 four-component radial convention mechanically checked.
- [x] Fe `spd` dimension reduction from 12,375 point coordinates recorded.
- [x] TDVK-02 material verdict left unchanged; TDVK-03 not started.
