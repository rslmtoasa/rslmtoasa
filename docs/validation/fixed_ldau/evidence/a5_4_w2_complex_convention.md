# A5.4 W2 preregistration and result

## A5.4: W2 complex-convention state test

For each fresh converged complex W2 endpoint with nonzero imaginary
off-diagonal density, capture `P` and perform dual replays using native `P`
and one joint `P^T=P*` control; preserve the production spherical
`m=-2..+2` basis and all other frozen state inputs. Before controls are
evaluated, preregister the correct-convention prediction, control
separations, and numerical admissibility bounds from replay repeat scatter,
roundoff, and capture precision. This avoids selecting a decision threshold
from the control outcomes. Since the frozen state is Hermitian,
`P^T=P*` up to the measured anti-Hermitian capture error, so transpose and
conjugation are one joint negative control and cannot be separated from each
other. Before either control replay, compute their expected separation from
the native replay from the captured `P` and W2 operator, then compare it with
the independently derived admissibility bound. If this separation exceeds
the bound, the documented complex convention should reproduce the native
result and the joint negative control should differ; otherwise report a
nondiscriminating control. The W2 Hellmann–Feynman slope is reported as a
diagnostic with its numerical uncertainty and is not gated.

Decision implication: this experiment assesses state-level complex convention
handling only. Its outcomes do not gate energy adoption or alter production
behavior.

The endpoint runner is
[`a5_w2_endpoint_capture.py`](../../../tests/validation/fixed_ldau/a5_w2_endpoint_capture.py).
It first converges each signed W2 endpoint from its own copy of the shared
near-converged restart using the existing SCF criterion and moments-only
capture. A second one-step run reads that converged `Fe_out.nml` restart and
captures the full trace at selected iteration 1 once. The output manifest
binds the endpoint restart, capture restart, input, binary, and records hashes.
The frozen replay adapter is
[`a5_w2_frozen_replay.py`](../../../tests/validation/fixed_ldau/a5_w2_frozen_replay.py).
For each converged W2 endpoint at `lambda=±0.01 Ry`, its `predict` phase first
checks the one-iteration capture, runs the unmodified `LIVE_STATE` radial
replay twice, then runs a predeclared matched calibration at quarter and half
of the native and joint-control center shifts and q2 moment increments. No
full-shift combined corner is evaluated before admissibility. For each spin
and each control leg, center-only SUMEV is compared with the direct linear
carrier target `q0*delta_c`. For both the center residual and q2-only response,
linear and quadratic extrapolations are calculated separately (factors 2/4
and 4/16 at half/quarter fractions); these are order-based uncertainty
envelopes under a smooth local Taylor assumption, not rigorous bounds. The
allowance takes the outer magnitude and adds the largest linear-versus-
quadratic disagreement. Independent center-only, q2-only and
combined small-shift cases provide a 2×2 factorial mixed term. Under the
predeclared smooth local Taylor assumption, the leading mixed term scales as
`delta_c*delta_q2`: the half-shift observation is multiplied by four and the
quarter-shift observation by sixteen to estimate its full-shift contribution.
This estimates the leading mixed term only; omitted higher mixed orders are
covered only insofar as the normalized shifts and quarter/half order
disagreement support the local expansion. The allowance is the larger
absolute extrapolation plus their disagreement.
The prereg also records each actual shift relative to the supplied radial
energy and second-moment scales; large ratios or inconsistent extrapolations
that undermine the local-order assumption make the diagnostic inadmissible/null rather than motivating a new fit or
gate. This finite-shift allowance, both live inversion residuals, repeat
scatter, and capture/roundoff bounds are compared with the preregistered
per-spin direct `Tr(WP)` separation before any full dual control is launched.
Only then is the hash-bound prediction written and the native/joint replay
allowed. The
replay uses the production W2 matrix
`W12=+i, W21=-i, W33=0.5, W44=-0.25` in the complex spherical `m=-2..+2`
block. It constructs the local center and centered second moment from the
captured `P_OUT`, `M1_OUT`, `q0/c/q2`, keeps all other captured radial QL and
radial-seed operands fixed, and applies the same transform to `P` and `M1` for
the joint transpose/conjugate control. It uses the existing production
`calculate_pl/NEWRHO` replay driver without importing A2 fits or its q2 and
response gates.

The prereplay prediction reports each spin's `|Tr(λW(P−Pᵀ))|` and their total;
admissibility is assessed per spin so cancellation cannot hide a discriminating
block. The two-times baseline residual is the sum of one measured
correct-convention and one wrong-convention radial inversion allowance; both
use the same frozen radial state. Add the measured center-only calibration
residual, q2-only response envelope, and absolute mixed center×q2 factorial
residual for each leg, then deterministic-repeat scatter and a floating-point
accumulation bound. The CSV writer's `ES24.16`
format gives a half-ulp relative bound of `5e-16` per normalized scalar; the
capture allowance uses 64 such rounded operands/terms to cover the 5×5 matrix
traces and moment-correction expression. The subsequent measured native-minus-
joint-control `SUMEV` separation and its residual against the frozen-state
prediction are reported as diagnostics. The Hellmann–Feynman slope remains
diagnostic only.

## Captured states and outcome

The full and moments-only captures were made as separate one-step runs from byte-identical converged endpoint restarts and inputs, with the same frozen binary. Full mode retained the radial replay seeds and wave arrays; moments-only retained direct `P_OUT` density matrices and `Tr(WP)`. Before merging only the 54 direct density/q0/trace rows into each full trace, the pair verifier checked matching `ETOT`, entropy/`TS`, per-spin `P_OUT` `ql`, radial moments and radial parameters; direct per-spin `q0` matched full-mode `ql` within `1e-14`. Raw captures remain separate and unchanged. This pairing is a capture completeness check; it makes no physics or convergence claim.

| λ (Ry) | spin | max \|Im P\| | max \|P−Pᵀ\| | predicted native−joint `SUMEV` (Ry) | measured (Ry) | measured−predicted (Ry) | allowance (Ry) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| −0.01 | 1 | 0.0115626 | 0.0231253 | 4.6250536e−4 | 4.6250084e−4 | −4.52e−9 | 3.537e−7 |
| −0.01 | 2 | 0.0231407 | 0.0462814 | 9.2562714e−4 | 9.2561920e−4 | −7.93e−9 | 8.749e−8 |
| +0.01 | 1 | 0.0117393 | 0.0234787 | 4.6957358e−4 | 4.6956901e−4 | −4.57e−9 | 2.949e−7 |
| +0.01 | 2 | 0.0229806 | 0.0459613 | 9.1922534e−4 | 9.1921749e−4 | −7.85e−9 | 1.088e−7 |

The summed measured native-minus-joint separations were 0.00138812005 Ry (−0.01 Ry) and 0.00138878651 Ry (+0.01 Ry), against predictions of 0.00138813250 and 0.00138879892 Ry. Absolute summed residuals were 1.2454e−8 and 1.2415e−8 Ry. The captured matrices were Hermitian to at most 1.22e−17; transpose and conjugation therefore form one joint control, as preregistered. The responses support the frozen-state complex-convention prediction for these states. The Taylor allowance is an order-based local uncertainty estimate, not a rigorous bound; omitted higher mixed orders remain a limitation. This result remains a diagnostic only and does not alter production behavior or gate energy adoption.

The endpoint direct `Tr(WP)` values (the Hellmann–Feynman slopes at each captured state) were 0.3712457148 ± 2e−15 at λ=−0.01 Ry and 0.2206149468 ± 2e−15 at λ=+0.01 Ry. This conservative per-endpoint capture-rounding allowance covers serialization of the four nonzero matrix-contraction operands and trace scalar plus binary64 accumulation; it covers serialization/accumulation only, not endpoint or SCF sensitivity. These are state diagnostics only; no slope threshold was preregistered or applied.

Capture binary SHA256: `52eb94c6e6f6ae12fca629a4eea1b30a6800e29736c9fdfdfd836538077f2225`. The signed endpoint restart hashes are `9080ce00c030b498affc0a398f3207043101f04b996c861bd4dd093ca5bfc3a6` (−0.01) and `2efff96d370c10956043c2b2d31e0326f505b6c19f9722127bff80cb45cd035f` (+0.01). Merged selected-state CSV hashes are `ee23c6fc640bda4e9cf144ad00a7023e21e43c0c818a42ac9570dd8893c80926` and `18518aca5841a9d09acafaed44db2bae3e9e597777ca2dffc03092a906399ec0`. Prediction hashes are `145abb9c1f0aca170abb9785c3abc3676f349376b1b0bd7f2e8113fc47a57a33` and `0effea8f554c999f1f36b5ec1d35babe13b9ebb0d1001c3c646b91ee53e71b9e`; the replay driver SHA256 is `c7cd56761af10d9b312a4f992519e44701af7d0edf1b4ac83006fecccb4f1af8`. Raw capture manifests, paired-capture verification manifests, prediction JSON, and control replay receipts are retained under `/private/tmp/a5-w2-paired-ks-v1/`; the compact tracked receipt is [a5_4_w2_result.json](a5_4_w2_result.json).
