# DRESP-08: native second-order LMTO field mapping

## Scope and result

DRESP-08 closes the rigid transverse-rotation mapping for the accepted,
no-SOC second-order LMTO Hamiltonian. It uses the live real-space `ee`,
`eeo`, `obarm`, and `enim` objects and the live transformed potential
parameters. It does not fit a coefficient-space field, modify Kxc, apply
BES/Halle, or run a new SCF campaign.

The bounded material campaign is the existing bcc Fe 4x4x4 accepted-state
handoff, with Gamma response and eta = 0.04, 0.02, 0.01, 0.005 Ry. The
independent artifact check reports:

```text
classification = NATIVE_SECOND_ORDER_MAPPING_CLOSED
verdict        = PASS
product dimension = 232
```

## Live mapping

```text
cx0/cx1, cex0/cex1, wx0/wx1, obx0/obx1
                         |
                         v
                  ee / obarm / enim
                         |
                         v
                 h(k), o, E_nu
                         |
                         v
             H^(2)(k) = E_nu + h - h o h
                         |
                         v
              native transverse tangent
                         |
                         v
                    product response
```

Measured residuals for each arrow are:

| mapping | residual |
|---|---:|
| spin-average identities for `cx`, `wx`, `cex`, `obx` | `0.0` |
| parameter reconstruction of `obarm` | `0.0` |
| parameter reconstruction of `enim = cx - cex` | `1.33e-16` |
| Fourier `eeo` versus `h o` | `4.26e-16` relative |
| reconstructed `H^(2)` versus `hk_bulk` | `9.50e-16` relative Frobenius |
| spin-resolved `B_H` versus direct `(H_up-H_down)/2` | `2.09e-15` relative |
| native product response versus exact static response | `1.91e-16` relative |
| native response versus exact-H rotation response | `7.35e-15` relative |

The largest reconstructed-H element residual is `4.45e-16`; the reconstructed
H hermiticity residual is `1.26e-12` in the stored coefficient matrices.

## Second-order Hamiltonian and field decomposition

The production contract is used without a surrogate:

```text
H^(2)(k) = E_nu + h(k) - h(k) o h(k)
```

For the collinear state the independent spin reconstruction is

```text
B_E   = 1/2 (E_nu,up - E_nu,down)
B_h   = 1/2 (h_up - h_down)
B_hoh = 1/2 (h_up o_up h_up - h_down o_down h_down)
B_H   = B_E + B_h - B_hoh
```

The weighted coefficient-orbital norms and fractions of `||B_H||` are:

| contribution | norm | fraction |
|---|---:|---:|
| `B_E` | `5.0223e-2` | `2.7781e-1` |
| `B_h` | `1.4193e-1` | `7.8507e-1` |
| `B_hoh` | `2.3229e-2` | `1.2849e-1` |
| `B_H` | `1.8078e-1` | `1.0` |

The nonlocal first-order hopping contribution is larger than the HOH
contribution in this state; the HOH term is nevertheless strongly
k-dependent. The measured k-dependence values are:

```text
Rk(B_E)       = 1.78e-15
Rk(B_h)       = 3.3254e-1
Rk(B_hoh)     = 7.1336e-1
Rk(B_H)       = 2.6749e-1
Rk(delta h)   = 3.3254e-1
Rk(delta H)   = 2.6749e-1
```

## Native tangents

The bond tangent is built by `lmto_bond_derivative` from the stored structure
block, `wx0/wx1`, and the active `cex1` endpoint factor on the HOH path. Its
Fourier transform is compared independently with `-i[G,h]`. The maximum
relative residuals over all 64 k points are:

```text
delta h    versus -i[G,h]       2.81e-15
delta o    versus -i[G,o]       2.09e-16
delta E_nu versus -i[G,E_nu]    9.77e-17
delta H    versus -i[G,H^(2)]   3.35e-15
```

The complete product rule is evaluated term by term:

```text
delta H = delta E_nu + delta h
          - delta h o h - h delta o h - h o delta h
```

Weighted mean norms over the mesh are:

```text
||delta E_nu||       7.1026e-2
||delta h||           2.0061e-1
||delta h o h||       2.3505e-2
||h delta o h||       4.4438e-2
||h o delta h||       2.3505e-2
||delta H exact||     2.5551e-1
```

The artifact contains the six norms and native-versus-commutator residual for
every accepted k point (`k_1_...` through `k_64_...`). The largest per-k
product tangent residual is below `3.5e-15`.

## Angular audit

The onsite parameter reconstruction explicitly uses the `L=0` form

```text
X_lm,l'm' = delta_ll' delta_mm' X_l
```

before the shared cartesian/spherical `hcpx` transformation. The independent
unit fixture verifies shell m-degeneracy and zero off-diagonal structure for
equal within-shell coefficients. In the Fe material field at Gamma, the
component diagnostics are:

| component | local spherical residual | same-l offdiagonal | within-l anisotropy |
|---|---:|---:|---:|
| `B_E` | `0.0` | `0.0` | `0.0` |
| `B_h` | `1.6512e-1` | `1.3494e-2` | `1.5966e-2` |
| `B_hoh` | `4.5898e-1` | `4.0800e-3` | `4.8275e-3` |
| `B_H` | `9.0148e-2` | `9.4141e-3` | `1.1139e-2` |

Cross-l and intersite components remain numerical noise:

```text
B_H cross-l   = 1.06e-16
B_H intersite = 0.0
```

The same five B_H projection diagnostics are emitted independently for all
64 accepted k points (`k_1_BH_...` through `k_64_BH_...`), rather than using
the Gamma result as a proxy for the mesh.

Thus the observed orbital structure is generated by the angular structure
constants, endpoint `wx` factors, `hcpx`, and the HOH products; it is not a
purely radial field replacement.

## Independent finite-rotation fixture

`UnitDresp08NativeSecondOrder` uses noncommuting complex orbital matrices for
spin-dependent `h`, `o`, and `E_nu`. It rotates all Pauli components by
`+/-theta`, forms the central difference of the complete `E+h-hoh` expression,
and compares it with the analytic product-rule tangent. The final central
difference residual is `2.60e-5` at the smallest fixture angle and the
successive-error ratio is `2.50e-1`, consistent with second-order central
difference convergence. Its analytic fixture residuals are all below
`4.1e-16`.

## Old radial diagnostic

The DRESP-07R radial map is retained as a diagnostic and is not fitted or
repaired. Against the native coefficient-space field:

```text
Bxc radial map: matrix 2.0125e-1, action 5.5047e-1
BKS radial map: matrix 2.0125e-1, action 5.5047e-1
```

The large mismatch is consistent with the omitted first-order hopping
transformation and overlap/HOH transformation; onsite-only comparisons also
remain different (`7.4190e-1` to `E_nu`, `8.5121e-2` to `E_nu+B_h`). This
milestone establishes rotation of accepted LMTO magnetic parameters into
`delta H^(2)`. The radial-field-to-potential-parameter Jacobian required by a
full ALSDA kernel remains a later milestone.

## Tests and exclusions

Focused unit coverage includes H² reconstruction, spin-resolved decomposition,
`lmto_bond_derivative`, overlap and E_nu tangents, the complete product rule,
finite rotation, L=0 angular selection, and product-response equivalence. The
Fe integration test reuses the accepted 4x4x4 handoff and does not launch a
new SCF campaign. The independent checker is
`tests/validation/dresp08_fe_artifact.py`.

BES/Halle, Kxc tuning, Goldstone correction, SOC, CCOR, orbital polarization,
Hubbard U/V, and local-axis corrections are all off or rejected by the bridge.
The frozen exchange and DRESP-03TG sources were not modified.

Recommended next milestone: certify the radial-field-to-transformed-LMTO-
parameter Jacobian that feeds this already closed native second-order mapping.
