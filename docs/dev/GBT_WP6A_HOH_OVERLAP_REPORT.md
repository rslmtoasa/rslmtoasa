# GBT WP6a HOH/overlap report

Date: 2026-08-03.

## Outcome

The HOH path is enabled for `gbt_single_q`. GBT is still applied exactly once,
at primitive directed `S`; `eeo` and `eeoee` are derived afterward and receive
no independent phase. The RS two-sweep and reciprocal second-order builders
now consume the same covariant operator.

Both generalized-overlap modes remain unsupported for GBT. The proxy metric is
explicitly documented as non-formal and the Kanpur metric is not implemented,
so either selection stops before eigensolution with the incomplete-formal-metric
diagnostic. The supported GBT HOH solve is the orthogonal second-order
Hamiltonian with identity overlap.

**WP6a contribution to G6: PASS. Overall G6: FAIL/open**, because CCOR and the
remaining WP6b/WP6c terms plus the integration gate have not been completed.

## Term map

| Object | Class | Definition/frame | GBT treatment | Status |
| --- | --- | --- | --- | --- |
| primitive `S_ij(R)` | directed bond | raw strux-lib orbital structure block | apply endpoint link `G_ij` here, before LMTO contraction | supported |
| `ee_ij` | directed bond | first-order `h_ij` assembled from linked `S_ij` | no later phase | supported |
| `obarm_j` | onsite | overlap-bar factor in the collinear rotating frame | no translational phase | supported |
| `enim_i` | onsite | linearization center in the collinear rotating frame | no translational phase | supported |
| `eeo_ij` | directed composite | `ee_ij obarm_j = (h o)_ij` | derive from linked `ee`; do not phase | supported |
| `eeoee_ij` | composite | historical same-bond `(ee_ij obarm_j) ee_ij^H`; unused by the operative HOH kernel | derive from linked factors; do not phase | diagnostic only |
| RS HOH | composite operator | two sweeps: `x -> h x -> (h o)(h x)` | inherits covariance from each linked `h` | supported |
| reciprocal HOH | composite operator | `eeo(k) h(k) = h(k)o h(k)` | ordinary Bloch sums only | supported |
| `ham_only` overlap | onsite identity | orthonormal second-order representation, `O(k)=I` | no phase; minimum eigenvalue exactly 1 | supported |
| `generalized_overlap_proxy` | incomplete composite | current `I + eeo(k)` proxy is not the formal LMTO metric | reject explicitly | unsupported in GBT |
| `generalized_overlap_kanpur` | missing composite | formal construction is not implemented | reject explicitly | unsupported in GBT |

`eeo` is intentionally not Hermitian by itself: `(h o)^H=o h`. Hermiticity is
required of the contracted correction `h o h` and the full second-order
Hamiltonian, not of this directed intermediate.

## Primitive-bond covariance derivation

Let `U_i` map the common rotating frame to the laboratory frame at site `i`.
The primitive-linked Hamiltonian bond and local factors are

```text
h^G_ij = U_i^H h_ij U_j
o^G_j  = U_j^H o_j U_j
e^G_i  = U_i^H e_i U_i .
```

Consequently the directed companion is

```text
(h^G o^G)_ij
  = U_i^H h_ij U_j U_j^H o_j U_j
  = U_i^H (h o)_ij U_j .
```

The full two-bond contraction is therefore

```text
(h^G o^G h^G)_ik
  = sum_j U_i^H h_ij o_j h_jk U_k
  = U_i^H (h o h)_ik U_k .
```

All intermediate endpoint frames cancel. This proves covariance and also shows
why phasing the already contracted `eeo` or `eeoee` would double count the
link. For a periodic operator, the convolution theorem gives

```text
F[sum_j (h o)_ij h_jk] = eeo(k) h(k),
```

which is exactly the reciprocal implementation. Since `h(k)`, `o`, and `e`
are Hermitian,

```text
H2(k) = e + h(k) - h(k)o h(k)
```

is Hermitian even though `eeo(k)=h(k)o` need not be.

## Validation

The independent dense oracle in `tests/unit/test_gbt_wp6_hoh.py` reports:

| Check | Result |
| --- | ---: |
| primitive-bond HOH covariance relative error | `2.103e-16` |
| RS two-sweep versus reciprocal relative error | `9.438e-17` |
| covariant second-order H Hermiticity error | `1.968e-16` |
| reciprocal first-order H Hermiticity error | `1.140e-16` |
| reciprocal second-order H Hermiticity error | `9.610e-17` |
| directed `eeo` non-Hermiticity (fixture must expose it) | `6.160e-01` |
| supported identity-overlap minimum eigenvalue | `1.0` |

The generalized eigensolver path now performs `zpotrf` on every supported
metric after the explicit Hermiticity check and before `zhegv`; indefinite
overlaps fail instead of becoming per-k solver warnings.

Finite-q bcc-Fe production smokes pass through reciprocal and RS block routes:

```text
wp6_hoh_k4   PASS
wp6_hoh_rs12 PASS
```

The reciprocal run selects `E_nu + h - hoh + L.S` and reports a maximum
pre-eigensolver `|H-H^H|` of `2.2888e-16`. The q=0 canonical band energy is
`-4.23275178 Ry`; q=+/-0.05 both give `-4.23388512 Ry`. The RS and reciprocal
production runs write byte-identical primitive bond dumps:

```text
SHA-256 7044ac05ae1cb3d648d2cde6b5055b8d1660b8165bb160a0aaea5312bc607474
```

Both negative production probes stop before eigensolution:

```text
GBT with reciprocal_mode=generalized_overlap_proxy is unsupported:
the available generalized-overlap representation is not a complete formal GBT metric

GBT with reciprocal_mode=generalized_overlap_kanpur is unsupported:
the available generalized-overlap representation is not a complete formal GBT metric
```

Build and unit validation:

```text
cmake --build build_13 -j2                         PASS
ctest --test-dir build_13 --output-on-failure -L unit
18/18 PASS
```

## Checklist

- [x] `eeo`/`eeoee`/`obarm`/`enim` classified.
- [x] HOH covariance derived at primitive-bond level.
- [x] RS two-sweep and reciprocal operators agree.
- [x] H/O Hermiticity checks pass.
- [x] Supported overlap is positive definite.
- [x] Unsupported modes fail early.
- [x] Term map and partial G6 PASS/FAIL reported.
