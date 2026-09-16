# DRESP-02C — Projected GF Material Closure

## Final verdict

`PASS — DRESP-02 MATERIAL CLOSURE COMPLETE`

DRESP-02F supplies the analytic equivalence. DRESP-02P supplies the
optimized material GF backend. This closure run supplies the missing finite-
regulator oracle and the resolved physical-target ladder.

## Frozen material contract

All Lehmann, finite-width-oracle, and production-GF values in the acceptance
artifact use one accepted bcc-Fe state: 4×4×4 full mesh, 64 k points,
`ham_only`, second-order Hamiltonian, orthogonal collinear baseline, no SOC,
300 K, exact folded `k+q` endpoints, shared k weights, EF, occupations,
radial snapshots, DRESP-01 projector, and site ordering. The final artifact
fingerprints are:

| quantity | value |
| --- | ---: |
| EF (Ry) | `-6.1212407827064663e-2` |
| accepted moment (`mu_B`) | `1.9507876274433729` |
| eigenvalue checksum (Ry) | `4.2910724817705404e3` |
| eigenvector checksum (real, imag) | `5.3747232042289318e2, -1.4725739671003728e2` |
| occupation checksum | `2.3482472833124111e3` |
| eigenvector unitarity residual | `2.33e-15` |
| radial checksum | `1.4274436288319852e4` |
| radial L2 norm | `2.6450730612148510e2` |

The projected moments were `d=1.9700111879` and `spd=1.9507876275 mu_B`.
The `spd` projected/accepted-total residual was `2.6e-11 mu_B`.

## Primary complex site-matrix evidence

The primary run used `eta_response=0.010 Ry`, `eta_int=0.002 Ry`,
`gf_energy_margin=0.60 Ry`, and `h/eta_int=0.3723`. It covered Γ, small
finite q `(+/-0.03,0,0)`, and finite q `(0.125,0,0)`, at
`omega=0, 0.005, 0.015 Ry`. The table shows representative complete-run
`chi_11` rows; the artifact retains every complex site-matrix element.

| projection | q | omega | chi Lehmann target (Re, Im) | chi finite-width oracle (Re, Im) | chi GF (Re, Im) | GF-oracle abs | GF-target abs |
| --- | ---: | ---: | --- | --- | --- | ---: | ---: |
| d | Γ | 0 | `(-26.0971,-1.54443)` | `(-26.0042,-1.24028)` | `(-26.0042,-1.24028)` | `5.5e-13` | `3.18e-1` |
| d | Γ | 0.015 | `(-28.6335,-1.88439)` | `(-28.5712,-1.53865)` | `(-28.5712,-1.53865)` | `1.9e-12` | `3.51e-1` |
| d | `(0.03,0,0)` | 0 | `(-26.0955,-1.52224)` | `(-25.9966,-1.22381)` | `(-25.9966,-1.22381)` | `2.3e-11` | `3.14e-1` |
| d | `(0.125,0,0)` | 0.015 | `(-27.9799,-1.79998)` | `(-27.8953,-1.49706)` | `(-27.8953,-1.49706)` | `1.2e-11` | `3.15e-1` |
| spd | Γ | 0 | `(-30.5844,-3.50571)` | `(-29.8655,-2.58915)` | `(-29.8655,-2.58915)` | `1.2e-11` | `1.16` |
| spd | Γ | 0.015 | `(-35.4206,-8.77265)` | `(-34.3578,-6.69503)` | `(-34.3578,-6.69503)` | `2.2e-11` | `2.33` |
| spd | `(0.03,0,0)` | 0 | `(-26.3936,-1.86245)` | `(-27.0420,-1.58786)` | `(-27.0420,-1.58786)` | `1.0e-11` | `7.04e-1` |
| spd | `(0.125,0,0)` | 0.015 | `(-28.4476,-1.89130)` | `(-28.4998,-1.59025)` | `(-28.4998,-1.59025)` | `4.0e-12` | `3.06e-1` |

Across all 24 primary rows, the maximum GF-to-oracle difference was
`2.6e-11` absolute and `9.1e-13` in the elementwise relative field. This is
finite-regulator equivalence at material precision, independently of the
physical-target comparison.

## Track A — finite-regulator equivalence

The oracle is `evaluate_projected_finite_width_chi0`: it rebuilds each
DRESP-01 transition amplitude directly and evaluates both spectral Kubo terms
with the finite-temperature `f(E)` on the same finite energy mesh. It does
not reuse the production eigenbasis GF contraction and does not use a blind
`eta_response + eta_int` Lehmann substitution.

At the controlled Γ samples, the global GF-to-oracle absolute/relative
residuals were:

| projection | eta_int (Ry) | NE | h/eta_int | GF-oracle abs | GF-oracle relative |
| --- | ---: | ---: | ---: | ---: | ---: |
| d | 0.008 | 933 | 0.3994 | `3.36e-12` | `7.3e-14` |
| d | 0.004 | 1863 | 0.3998 | `5.09e-12` | `1.1e-13` |
| d | 0.002 | 3725 | 0.3998 | `4.59e-12` | `9.7e-14` |
| d | 0.001 | 7447 | 0.3999 | `2.35e-11` | `5.0e-13` |
| d | 0.0005 | 14893 | 0.3999 | `4.30e-11` | `9.1e-13` |
| spd | 0.008 | 933 | 0.3994 | `1.09e-12` | `2.0e-14` |
| spd | 0.004 | 1863 | 0.3998 | `5.85e-12` | `1.1e-13` |
| spd | 0.002 | 3725 | 0.3998 | `2.52e-11` | `4.5e-13` |
| spd | 0.001 | 7447 | 0.3999 | `8.00e-11` | `1.4e-12` |
| spd | 0.0005 | 14893 | 0.3999 | `2.32e-10` | `4.1e-12` |

## Track B — physical-target limit

At fixed `eta_response=0.010 Ry`, the finite-width oracle was compared with
the Lehmann target using the original physical broadening. The raw sequence
was:

| projection | eta_int (Ry) | oracle-target abs | oracle-target relative |
| --- | ---: | ---: | ---: |
| d | 0.008 | `2.452` | `5.19e-2` |
| d | 0.004 | `1.180` | `2.50e-2` |
| d | 0.002 | `5.806e-1` | `1.23e-2` |
| d | 0.001 | `2.855e-1` | `6.04e-3` |
| d | 0.0005 | `1.424e-1` | `3.01e-3` |
| spd | 0.008 | `7.935` | `1.37e-1` |
| spd | 0.004 | `4.972` | `8.61e-2` |
| spd | 0.002 | `2.975` | `5.15e-2` |
| spd | 0.001 | `1.650` | `2.86e-2` |
| spd | 0.0005 | `8.775e-1` | `1.52e-2` |

The ladder is resolved, approaches zero in both projections, and is reported
without monotonicity assumptions or empirical extrapolation. DRESP-02F
provides the exact full-line finite-spectrum limit; the material ladder
confirms that the implementation follows that limit toward
`chi_Lehmann(eta_response)`.

## Eta semantics and regression

`eta_response` is the physical retarded response broadening. `integration_eta`
is the numerical one-electron real-axis regulator and is not a physical
linewidth. Intrinsic linewidth/Landau damping is not determined by DRESP-02.

The DRESP-02F finite-model tests remain in the regression suite. The focused
projected unit test now also checks the independent finite-width oracle; the
largest new GF-to-oracle residual was `2.3e-17` in the two-site `spd` fixture.

Machine-readable output from the acceptance run:

`/tmp/dresp02c_fe/dresp02c_Fe_projected_chi0.dat`
