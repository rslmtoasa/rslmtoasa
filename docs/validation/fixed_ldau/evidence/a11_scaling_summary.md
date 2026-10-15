# A1.1 KS scalar-shift scaling result

Verdict: **FAIL_A1_1**. The accepted-state equivalence gate passed before any shift fits were inspected. Both physical spin channels and their summed diagnostic show a finite fitted zero-shift relative slope error larger than the preregistered `1e-6` limit. The signed residuals are smooth and strictly decrease toward zero on both sign branches; the all-eight-point log orders exceed two and the fitted cubic constants are below `10 D0`. Thus this is a measured A1.1 slope-extrapolation failure. It does not demonstrate a flat/nonmonotonic residual or establish that the A1 error model is wrong under A1.3. Do not proceed to M0.3.

The summed d-shell occupation is `q0=6.498480421874168`. For signed shifts ordered as `[-0.02,-0.01,-0.005,-0.0025,+0.0025,+0.005,+0.01,+0.02] Ry`, the measured total `r(v)=DeltaSUMEV-v*q0` is:

| v (Ry) | r(v) (Ry) |
|---:|---:|
| -0.0200 | 8.225465541830679e-5 |
| -0.0100 | 1.063551491135850e-5 |
| -0.0050 | 1.521577795059181e-6 |
| -0.0025 | 2.873991285894106e-7 |
| +0.0025 | -2.875038543730113e-7 |
| +0.0050 | -1.495376534342563e-6 |
| +0.0100 | -1.053510915298428e-5 |
| +0.0200 | -8.089596527213750e-5 |

The symmetric slopes for `v=(0.0025,0.005,0.01,0.02) Ry` are `(6.498365441277576, 6.498178726441228, 6.497421890670951, 6.494401656356907)` in q0/occupation units; relative errors `s/q0-1` are `(-1.769345895163976e-5, -4.642553541045480e-5, -1.628890347432632e-4, -6.276491198667911e-4)`. Their even-quartic extrapolated intercept is `a0=-7.865137674552838e-6`, failing `|a0|<1e-6`. The linear-in-`v²` intercept is `-7.879051651808900e-6`; leave-one-magnitude-out quartic intercepts are `(-7.580059806013838e-6, -8.014464177120492e-6, -8.123065269897124e-6, -8.150215543022043e-6)`. The preregistered sensitivity diagnostic is `2.850778685389999e-7`; it is reported, not added as an allowance or used to alter the gate.

The primary all-eight signed log-log order is `2.725679793563955`; positive/negative branch fits are diagnostics only. The total cubic constant is `-5.111306142956791e-8 Ry`, compared with `10 D0=5.028166848219939e-5 Ry`. Both spin channels independently fail the same slope-intercept criterion: up `a0=-7.581540987621142e-6`, down `a0=-8.435949493162413e-6`.

## Provenance and state gate

- Candidate run: `/private/tmp/rslmto-localmom-a11-accepted-recapture`, iteration 23, selected from 23 SCF updates; run exited 0.
- Candidate Release executable SHA256: `7a98a68a11f6a4a062e1f8256ffb80a64b4dc8967c0165dab1ffb8910425a3fa`.
- Candidate source `self.f90` SHA256: `5f9cf847738803a49fc8d35c51be7b5dd3e7083c851a5c07d1785e07bc1dab89`.
- Candidate CSV SHA256: `55ef74516393256f262de6f02ecf0f061536ff9bff62641d1e6e828149586ac4`; final manifest SHA256: `919b0c6ae55eba1e9d402fec4f43b071e6ccc192f377647a541027e132e70025`.
- Accepted reference CSV SHA256: `4dc77e488175e342254409d8327df70334a0cd688855641491f670aac8a704b4`; no-U floor JSON SHA256: `dd33dd828131d1a81df9f2ceb729aaae7116c0af7ef8c78ac2cc17a56b090bef`.
- Candidate and accepted reference PRE input/seed hashes match exactly. Candidate basis size 18 was counted from its captured basis map; declared temperature and electron count were taken from the reference manifest because the exact PRE deck hashes match. Original manifests were left unchanged.
- Independent accepted-state check: `max|Delta P|=0`, `max|Delta M1|=0 Ry`, `max|Delta(q0,c,q2)|=0`, `max|Delta calculate_pl inputs/returns|=0`, and frozen radial/replay inputs are bitwise identical. Tolerances were `5e-12` for P/M1, `1e-10` for radial moments and calculate_pl, and exact serialized-double equality for frozen inputs. Hash-bound sidecar: [`a11_accepted_state_equivalence.json`](a11_accepted_state_equivalence.json), SHA256 `1a7238ed4f98649feb8e654c04209b987d80cd715d03eea103b64b8a0a7ef005`.
- Full numeric output, including every spin-resolved signed residual, fit coefficient, diagnostic, and the aggregate channel, is in [`a11_scaling_result_reviewed.json`](a11_scaling_result_reviewed.json), SHA256 `48b0cb79291b12fd5e3c027b92379f0f70b45e8418882348616a1ea1abed27ae`. Historical `a11_scaling_result.json` and v0 evidence were not overwritten.

## Independent command replay

The capture writer's original `manifest.json` does not include the checker
schema, CSV hash key, or the parser's `iteration` field. It was not modified.
A sidecar manifest was created in a separate copy of the captured CSV, adding
only those required keys plus basis size, electron count, temperature, and
aliases for the binary/source hashes. The sidecar is
[`a11_candidate_manifest_supplement.json`](a11_candidate_manifest_supplement.json),
SHA256 `c152b3a6822f63ab61a6ff6ce5895b7b6fb860eeae52e4e7d3504884c472d160`.
It binds the unchanged CSV SHA `55ef74516393256f262de6f02ecf0f061536ff9bff62641d1e6e828149586ac4`.

The exact reproduction command was:

```sh
python3 tests/validation/fixed_ldau/m0_a11_scaling_oracle.py \
  /private/tmp/rslmto-localmom-a11-reproduce \
  --reference-dir /private/tmp/rslmto-localmom-m02-release/ks \
  --floor-json /private/tmp/rslmto-localmom-m02-release/ks/m02-floor.json \
  --output /private/tmp/a11-independent-reproduce-a1.json
```

The analyzer SHA256 was `2537171d8bbeaa8eb7043510a28c6625e8534f215441b87de52b3c0917c81bd1`.
It exited 0 with `FAIL_A1_1`; its complete output is
[`a11_independent_reproduction.json`](a11_independent_reproduction.json),
SHA256 `cf8102e0f250be79353012fd7136dbc16fc0b9921c90a8bc0e3cf0863b390c62`.
The reproduced verdict and all spin/aggregate numerical fit values exactly
match the reviewed result. The reviewed JSON adds sidecar/reference links and
provenance annotations; it does not alter the measured values. The older
`oracle_sha256` in the capture manifest identifies the M0 capture checker,
not this A1 fitting analyzer. The first A1 result record remains unchanged at
[`a11_scaling_result.json`](a11_scaling_result.json), SHA256
`8cc38871de8dc908e109185eaf09129eff6ea08fb787f66d6b23bfc14c965ea4`; it
predates the reviewed provenance/aggregate additions and is not the final
reviewed report.
