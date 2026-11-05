# A5.3 local and spectral J derivative review

This additive review completes the J comparison that the earlier
`a5_m23_spectral_j_review.md` left pending. The numeric allowances were
registered in `a5_m22_numeric_allowance_preregistration.md` before candidate
J energy rows were inspected. It preserves that earlier spectral-only report
as a historical checkpoint; the current accepted local and spectral endpoint
summaries are listed below.

For candidate derivatives, I use the registered centered h=.05 eV result as
the primary finite difference and h=.10 eV as the measured step-sensitivity
diagnostic. The reported Richardson values below are descriptive estimates
from those two steps; they do not replace the registered allowance or add a
new gate. The decision comparison is the candidate derivative against the
production fixed-P `dE_Hub/dJ` target with the independently propagated
charge-only interval.

| radial mode | J center (eV) | candidate b derivative (Ry/eV) | production fixed-P target (Ry/eV) | residual (Ry/eV) | registered interval (Ry/eV) | result |
|---|---:|---:|---:|---:|---:|---|
| local | .5 | −.06389414333 | −.06352068324 | −.0003734601 | [.00118574,.00174694] | within allowance |
| local | .9 | −.06468169189 | −.06480683083 | +.0001251389 | [.00117606,.00258153] | within allowance |
| spectral | .5 | +.06312367005 | −.06304707741 | +.1261707475 | [.00204434,.00416579] | outside allowance |
| spectral | .9 | +.03725064057 | −.06771349680 | +.1049641374 | [.00057640,.00113991] | outside allowance |

Thus local b is within the registered charge-only allowance at both J
centers, while spectral c misses the target at both centers by far more than
the allowance. Candidate c is the relevant spectral candidate; its large
residuals are the decision result. The measured `dt/dJ` Richardson scales
are +.13047598 and +.10890306 Ry/eV, but the candidate-c residuals differ
from those scales by about −.004305 and −.003939 Ry/eV, respectively. Those
offsets are not covered by the allowance, so I do not describe c as matching
`dt/dJ` within uncertainty. Candidate a is report-only (spectral Richardson
slopes −.06735231 and −.07165242 Ry/eV); no a-versus-b preference is claimed.

The local response and fixed-P target reduction is recorded in
`a5_m23_local_j_allowance.json` (local endpoint summary SHA256
`e7463ccb1417e4e9acce0157b6ab1e6bf19d1d993055e0bf9fc72035cf5c9a4d`). The
spectral endpoints are in `/private/tmp/a5-m22m23-remaining-ks-grid/summary.json`
(SHA256 `4ce3967ade9b02d27fef40ba245fa212fd608b11bb1d6416ddad55f41f4ebda0`).
All ten local and ten spectral J endpoints are accepted native-KS runs with
finite entropy/TS and R05 PASS according to those summaries. The local
production fixed-P targets use exact direct `P_OUT` center matrices and the
production evaluator; their P and evaluator-output hashes are in the local
allowance JSON. The production evaluator is the fixed-P driver source
`14290820cc388d6229d1991b34873e5a5e6818444dc8a5f7ffdfd9fcd94dbc0e` and
binary `3c0d62b3666cd27c6064d4fd73c58358b05d27f11a2047a61ac12fa780ca417a`.

I independently rebuilt the ten local J candidate derivatives from each
endpoint's final `test_density_matrix` rows (`stage=P_OUT`, 50 matrix entries)
and final `HUBBARD_MOMENT_ENERGY_TOTAL` line. The production A4.3 evaluator
recomputed E_Hub from each direct P_OUT matrix; it matched the final logged
E_Hub at the printed precision for all ten states. Reconstructed candidate b
free energy as `physical_total_energy + E_Hub(P_OUT) - TS`, then recomputed
the h=.05 and h=.10 centered differences and Richardson descriptive values.
The raw records SHA256 values for all ten endpoint states are in
`a5_m23_j_derivative_result.json`; the center P matrices and target output
hashes are in `a5_m23_local_j_allowance.json`. The production evaluator is
the source and binary identified above. For spectral candidate-energy
provenance, the separate direct-P audit ran the evaluator on 14 direct P_OUT
matrices and found a maximum absolute difference of `4.163336342344337e-17
Ry` from the logged Hubbard energy. This report adds no statistical
interpretation to measured step differences and does not alter the preserved
coarse W1 allowance.
