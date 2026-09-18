# TG-FZ-R7 — production screening-representation repair

Status: **PASS-A — production screening representation repaired.**

R7 applies the two representation repairs proven by TG-FZ-R6 and reruns the
affected fixed-z gates. It stops before mixed derivatives, complete
curvature, energy integration, native `J_ij`/`J(q)`, and DRESP-04.

## Legacy screening authority

The legacy call chain is:

```text
lattice%structb(.true.)
  -> dbar1(ia, ..., ii)
     -> micha(...)
        -> SHLDCH(...)
     -> sbar(:,:,m,ii) = screened block
```

`dbar1` identifies the symbolic atom associated with the produced block as
`iz(ia)`. R7 publishes the exact target on that potential immediately after
`micha` returns and before the `sbar(:,:,m,ii)` blocks are retained. The
legacy `micha` factors and the published target share one definition in
`lattice.f90`:

```text
legacy_micha_alpha(0:3) = [0.3485, 0.05303, 0.010714, 0.00337]
q = fak * legacy_micha_alpha, fak = 2
```

The legacy `SHLDCH` operation and its historical `sbar=2*s` convention are
unchanged. The newer `strux_lib` path and its existing stored-alpha ownership
are unchanged; R7 does not overwrite them.

For the bcc-Fe spd fixture:

| l | published `screening_alpha` | predls target after repair |
|---:|---:|---:|
| 0 | 0.34849999999999998 | 0.34849999999999998 |
| 1 | 0.05303000000000000 | 0.05303000000000000 |
| 2 | 0.01071400000000000 | 0.01071400000000000 |

The allocation, publication, and legacy-value contract residuals are all
`0.0000000000000000e+00`. The post-`predls` live arrays agree with an
independent local reconstruction from the published target with maximum
center/width/obar residual `2.2204460492503131e-16`.

The old predls target was `qm_canonical`, with `l=0` equal to
`0.34848499999999999`. The repaired target is `0.34849999999999998`; the
change is intentional and comes from the structure producer, not a numerical
fit.

## Native gamma repair

`native_complex_p_matrix` constructs

```text
P_norm = (z - (C + Vmad)) / DELE**2
```

where `DELE` is `potential%dele`. Therefore
`native_screened_p_matrix` now uses `potential%qi` as its input gamma. The
Turek screening formula itself is unchanged:

```text
P_alpha = P_norm / (1 + (gamma_norm - alpha) * P_norm)
```

The raw relation remains independently checked at all active l/spin channels
and fixed energies:

```text
qpar * P_raw = qi * P_norm residual = 6.2803698347351007e-16
```

The fixed-z test now labels the repaired quantity `gamma_norm=qi`; raw
`gamma_raw=qpar` is retained only for the pre-repair diagnostic.

## Potential-parameter change

Compared with the old `qm_canonical` target, the maximum absolute changes in
the live post-`predls` parameters are:

| parameter | maximum absolute change |
|---|---:|
| `center_band` | `2.3276121113857684e-06` |
| `shifted_band` | `2.3276121114412796e-06` |
| `width_band` | `5.9088223591241551e-06` |
| `obar` | `9.3945043344534351e-05` |

These are measured consequences of consuming the structure-authoritative
target; no scale or fit was introduced.

## Static exact-H/native-gamma gate

The repaired live state uses the actual published-alpha `sbar`, actual live
post-`predls` `center_band/shifted_band/width_band/obar`, and normalized
`gamma_norm=QI`.

| gate | maximum residual |
|---|---:|
| Route A, structure alpha | `6.6613381477509392e-16` |
| Route B, predls alpha | `6.6613381477509392e-16` |
| live `||H_exact_live-H_gamma_norm||_max` | `4.4408920985006262e-16` |
| live endpoint-resolvent closure | `9.5659743364509064e-15` |
| R7/R6 normalized `S_gamma` cross-check | `1.3322676295501878e-15` |
| R7/R6 normalized `H_gamma` cross-check | `2.2204460492503131e-16` |

The pre-repair diagnostic remains visible:

```text
gamma=qpar, alpha=alpha_structure: 3.3750056143588836e-05
gamma=qi,    alpha=alpha_structure: 4.4408920985006262e-16
```

This isolates the former native-helper mismatch without changing the
historical exchange route.

## R1 revalidation

The repaired fixed-z R1 gates use normalized gamma wherever `P_norm` is used:

| gate | maximum residual |
|---|---:|
| real-axis `P` reduction | `0.000000e+00` |
| `d_matrix` width-scaled `DeltaP_norm` identity | `2.827599e-16` |
| direct one-site inverse | `3.409427e-16` |
| fixed-z global inverse | `4.475452e-16` |
| `P` transformation | `1.588822e-14` |
| `S` transformation round trip | `1.332268e-15` |
| path-operator covariance | `4.613190e-16` |
| `gamma_norm` versus finite-H, ordered `ud` | `4.440892e-16` |
| `gamma_norm` versus finite-H, ordered `du` | `8.881784e-16` |
| Pauli helper self-consistency | `0.000000e+00` |

The nonzero `alpha-gamma_norm` contraction difference remains a physical
active-alpha versus normalized-gamma vertex distinction; it is not a failed
covariance gate.

## R2 revalidation

The transformed vertex now uses

\[
 \widetilde{\Delta P}^{\,\gamma_{norm}}
 = \Delta P^{\gamma_{norm}}
 + (QI_\downarrow-QI_\uparrow)
   P^{\gamma_{norm}}_\uparrow P^{\gamma_{norm}}_\downarrow .
\]

The corrected maxima are:

| gate | maximum residual |
|---|---:|
| endpoint vertex expansion | `1.517720e-14` |
| expanded screening identity | `4.143141e-14` |
| direct alpha / transformed-gamma contraction | `8.881784e-16` |
| common-gamma vertex control | `5.402578e-15` |
| common-gamma contraction | `8.881784e-16` |
| `delta_d_screen` diagnostic magnitude | `5.423866e-01` |

The last value is a reported spin-screening correction, not a tolerance gate.
The maximum direct alpha versus normalized-gamma contraction mismatch remains
`7.648525e-02`, as expected for distinct active-alpha and gamma vertices.

## Historical exchange and regression

The historical exchange route is unchanged. `source/exchange.f90` was not
modified, and no historical `qpar` convention was changed in this task.

Verified commands:

```text
cmake --build build -j2
ctest --test-dir build -V -R '^UnitDresp03tgNativeFixedZ$' --output-on-failure
ctest --test-dir build -V -R '^UnitStructureConstantsBackends$' --output-on-failure
ctest --test-dir build -R 'UnitDresp03|UnitLrKlStaticBridge' --output-on-failure
```

The focused R7 test, structure-backend ownership test, and the requested
DRESP/LMTO regression selection pass. No derivative, curvature, integration,
Jij, J(q), or DRESP-04 work was entered.

The frozen exchange-file SHA-256 remains:

```text
6e9ae7da6af46367aaf18a4fef37bdd4a9a0d6f69fe8e08b589e5c5f278a70e6
```
