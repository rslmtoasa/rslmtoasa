# TG-FZ-R6 — screening normalization and target-alpha authority audit

Status: **PASS-B — normalized gamma and alpha authority are both required.**

This is the terminal fixed-complex-energy static audit for DRESP-03TG. It
does not modify production `predls`, the structure-constant backends, the
native Turek helpers, or `source/exchange.f90`. The diagnostic is implemented
in `tests/unit/test_dresp03tg_native_fixed_z.f90` and stops before derivatives,
curvature, energy integration, or exchange constants.

## Normalization ledger

The live `predls` source uses, with `I=l+1`,

```text
wow  = wsm/ws_r
DELE = srdel * wow**(-(l+1/2))
QI   = qpar  * wow**(-(2*l+1))
```

For the fixture, `wow = 1.0000172687024567`. The per-channel trace is:

| l | spin | `srdel` | `dele/srdel` | expected width scale | `qpar` | `qi/qpar` | expected screening scale |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1 | 4.3125676390000001e-1 | 9.9999136576059811e-1 | 9.9999136576059811e-1 | 4.2934432610000001e-1 | 9.9998273159574624e-1 | 9.9998273159574624e-1 |
| 0 | 2 | 4.3523806780000002e-1 | 9.9999136576059811e-1 | 9.9999136576059811e-1 | 4.3009644619999998e-1 | 9.9998273159574624e-1 | 9.9998273159574624e-1 |
| 1 | 1 | 4.1469395149999999e-1 | 9.9997409750544375e-1 | 9.9997409750544386e-1 | 1.1493227990000000e-1 | 9.9994819568182691e-1 | 9.9994819568182702e-1 |
| 1 | 2 | 4.1692183080000000e-1 | 9.9997409750544386e-1 | 9.9997409750544386e-1 | 1.1499369130000001e-1 | 9.9994819568182691e-1 | 9.9994819568182702e-1 |
| 2 | 1 | 1.1564741240000000e-1 | 9.9995682954848497e-1 | 9.9995682954848497e-1 | -1.9976922999999998e-3 | 9.9991366096065770e-1 | 9.9991366096065770e-1 |
| 2 | 2 | 1.3009218430000000e-1 | 9.9995682954848497e-1 | 9.9995682954848497e-1 | 4.3530777000000001e-3 | 9.9991366096065781e-1 | 9.9991366096065770e-1 |

The maximum width and screening normalization residuals are both
`1.1102230246251565e-16`.

For `s_l = wow**(-(2l+1))`,

\[
  DELE_l^2 = srdel_l^2s_l,
  \qquad P_{norm}(z)=s_l^{-1}P_{raw}(z).
\]

Therefore

\[
  qpar_lP_{raw}=qpar_ls_lP_{norm}=QI_lP_{norm}.
\]

The independent fixed-energy checks at `-0.91+0.83i`, `-0.17+0.04i`, and
`0.62+0.31i` give a maximum
`||qpar*Praw-qi*Pnorm||_max = 6.2803698347351007e-16`.

The corresponding channelwise screening rule is also fixed by the algebra.
If a raw representation uses

\[
 P_{out,raw}=P_{raw}/[1+(\gamma_{raw}-\alpha_{raw})P_{raw}],
\]

then radial normalization gives

\[
 \gamma_{norm}=s_l\gamma_{raw}=QI_l,
 \qquad \alpha_{norm}=s_l\alpha_{raw}.
\]

Thus a `P` formed with `potential%dele` is a normalized `P` and must use
`potential%qi` as its input gamma. A raw screening constant carried into the
normalized representation must receive the same `s_l` scaling. The two
source-produced target-alpha arrays below are tracked as their own target
ledgers; they are not inferred from `qpar`.

## Target-alpha authority

The actual fixture selects the legacy structure backend and has no allocated
`potential%screening_alpha`.

The structure target is produced by
`lattice_strux:micha -> SHLDCH`: `micha` sets
`q = 2*[0.3485, 0.05303, 0.010714, ...]`, `SHLDCH` uses `bet=1/q`, and the
legacy path stores `sbar=2*s`. The effective target-alpha ledger for the
stored `sbar` is therefore:

| l | `alpha_structure` | source |
|---:|---:|---|
| 0 | 0.34849999999999998 | legacy `micha`/`SHLDCH`, pre-`fak=2` label |
| 1 | 0.05303000000000000 | legacy `micha`/`SHLDCH`, pre-`fak=2` label |
| 2 | 0.01071400000000000 | legacy `micha`/`SHLDCH`, pre-`fak=2` label |

The `predls` target starts from `math_mod::qm_canonical`; there is no stored
alpha override in this fixture:

| l | `alpha_predls` | `delta_alpha = alpha_predls-alpha_structure` |
|---:|---:|---:|
| 0 | 0.34848499999999999 | -1.4999999999987246e-5 |
| 1 | 0.05303000000000000 | 0.0000000000000000 |
| 2 | 0.01071400000000000 | 0.0000000000000000 |

The nonzero `l=0` difference is small numerically but is not discarded.

## Representation-consistent routes

The test first converts the live spin-major `sbar` to the site-major ordering
used by the two-site Hamiltonian fixture. It then applies the exact matrix
screening transform; no relabeling or fitted scale is used.

Route A uses `alpha_structure` for the local TB reconstruction and transforms
`S^(alpha_structure)` directly to `S^(gamma_norm)`, with
`gamma_norm=QI`.

Route B first transforms the live structure matrix
`S^(alpha_structure) -> S^(alpha_predls)`, reconstructs the TB parameters with
`alpha_predls`, and then transforms to the same normalized `gamma_norm=QI`.

| gate | maximum residual |
|---|---:|
| Route A, `||H_exact_A-H_gamma_A||_max` | `6.6613381477509392e-16` |
| Route B, `||H_exact_B-H_gamma_B||_max` | `1.1102230246251565e-15` |
| Route A endpoint-resolvent closure | `9.5659743364509064e-15` |
| Route B endpoint-resolvent closure | `8.0703363304951411e-15` |

The static map and both normalized endpoint-resolvent identities therefore
close to roundoff.

## Live mixed convention and localization diagnostics

The production combination in this fixture is the live structure matrix at
`alpha_structure`, the `predls` TB parameters at `alpha_predls`, and the
current native helper's `qpar` gamma. The R6/R5 implementation cross-checks
the current native `S_gamma` and `H_gamma` at `8.8817841970012523e-16` and
`2.2204460492503131e-16`, respectively.

| comparison | maximum residual |
|---|---:|
| `||H_exact_mixed-H_gamma_current||_max` | `2.4279838737983894e-5` |
| `||H_exact_mixed-H_exact_A||_max` | `5.8029894881572730e-5` |
| `||H_exact_mixed-H_exact_B||_max` | `5.8029894881128641e-5` |

The four localization combinations are:

| gamma | alpha label | residual |
|---|---|---:|
| `qpar` | `alpha_structure` | `2.4279838737983894e-5` |
| `qpar` | `alpha_predls` | `3.3755019531955810e-5` |
| `qi` | `alpha_structure` | `5.8029894881128641e-5` |
| `qi` | `alpha_predls` | `8.8817841970012523e-16` |

These four rows are diagnostic-only because the alpha-mismatched rows retain
the live source matrix or its label. Only Route B performs the required
`S^(alpha_structure) -> S^(alpha_predls)` transformation before claiming
closure. The table nevertheless isolates both required convention changes:
the current `qpar` gamma fails, and the live structure/predls alpha mismatch
must be resolved by transforming `S` rather than relabeling it.

## Native helper audit

`native_complex_p_matrix` constructs its diagonal `P` with
`potential%dele`, so it constructs `P_norm`. `native_screened_p_matrix`
currently applies the input gamma from `potential%qpar`. R6 proves that this
pair is a representation mismatch: the helper should use `potential%qi` for
the normalized `P` path. `native_screening_alpha` correctly prefers a stored
screening array and otherwise uses its legacy fallback; the remaining legacy
issue is the independent target-alpha authority described above.

The production helper is intentionally unchanged in R6. The repair belongs
to a follow-up production task after this audit.

## R1/R2 and historical exchange implications

```text
R1/R2 algebra: valid
R1/R2 representation label: requires correction from qpar-gamma to normalized qi-gamma
```

R1/R2 numerical ledgers that consume the normalized native `P`/screening path
would need rerunning after the helper repair; no broad revalidation is part of
R6. The historical exchange route remains unchanged. If it retains a
`potential%qpar` auxiliary-GF transform, that is a separate historical
representation question and is not evidence that `source/exchange.f90` is
wrong.

## Verification and frozen file

The focused R6 test passes after:

```text
cmake --build build -j2
ctest --test-dir build -V -R '^UnitDresp03tgNativeFixedZ$' --output-on-failure
ctest --test-dir build -R 'UnitDresp03|UnitLrKlStaticBridge' --output-on-failure
```

The requested aggregate regression passed all 8 tests.

`source/exchange.f90` remains unchanged and has SHA-256:

```text
6e9ae7da6af46367aaf18a4fef37bdd4a9a0d6f69fe8e08b589e5c5f278a70e6
```
