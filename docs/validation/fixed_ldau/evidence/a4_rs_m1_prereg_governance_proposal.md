# RS M1 paired-replay preregistration governance proposal

Status: reviewer proposal for root/user resolution. No response is changed,
rerun, relabeled, or removed by this note. The locked preregistration and its
measured report remain the governing historical record unless the user
explicitly approves a prospective interpretation.

## Locked outcome

The v1 preregistration `a4_rs_m1_pair_preregistration.json` (SHA256
`87cc23819d6f9e465593e6a4708ea062c503db69506097d6643dee549d4d54a4`) explicitly
made both the literal carrier test and a residual-inclusive fit-exactness test
criteria. The preserved paired report `pout_pair/report.json` (SHA256
`b4023c1f0e3b11c46ccb2c8a0e31b007cbf6325d9ba4349a7faf25d8536a9587`) records:

* `Delta SUMEV (spectral - local) = -0.2590196642435809 Ry`;
* `sum(t) = -0.26088126108452275 Ry` and literal error
  `|Delta - sum(t)| = 0.001861596840941837 Ry`;
* locked `A_RS = 0.0019179425667351765 Ry`, so the literal carrier criterion
  passes even without adding the `2e-8 Ry` numerical guard;
* local and spectral minimum q2 are `0.01075725460260394` and
  `0.009745244577861018 Ry^2`, respectively, so q2 nonnegativity passes;
* the additional residual-inclusive fit-prediction error is
  `1.0234967739219947e-5 Ry` against its separately locked `2.0002e-8 Ry`
  limit, so that criterion fails.

Accordingly, the preregistered overall result is a failure as written. This
proposal does not authorize replacing that result or changing the v1 lock.

## Why the additional fit-exactness bound is not a sound carrier decision

The protected physical decision is whether the measured energy response uses
the input carrier `t = Tr(V_in P)` rather than a mistaken operator. The
preregistered A2 allowance is already `A_RS = sum_s |r_fit(x_s)| + D0_RS`;
it bounds the response against the literal carrier `sum(t)`. The paired
result's literal error lies inside that allowance.

The separate test `|Delta - (sum(t) + sum(r_fit))| <= 2e-8 + 2e-12 Ry`
requires the approximate fitted residual to predict the response to near
floating-point/comparator precision. The `2e-8 Ry` energy guard only bounds
comparison noise in the replay outputs; it does not bound A2 cubic-fit
approximation error, finite-grid radial remapping, or the response of the
reconstructed q2 terms. Those effects are covered by the existing A2
allowance for the literal carrier decision, not by a near-zero exactness
bound on the fitted residual. This explains why the additional test can fail
while the intended carrier test passes. It does not erase the fact that the
extra test was frozen and failed.

## Prospective interpretation for approval

For this paired frozen-state M1 decision only, apply the already-frozen
literal carrier comparator

```
|Delta SUMEV_(spectral-local) - sum_s Tr(V_in,s P_out,s)| <= A_RS
```

and the already-frozen physicality requirement `q2_local >= 0` and
`q2_spectral >= 0`, together with the independent LIVE anchor, source/state
identity, unchanged-operand, and registered output checks. Report the
residual-inclusive fit-prediction discrepancy as a diagnostic with its
measured value and the failed v1 auxiliary criterion. Do not call the v1
registration itself a pass. If this prospective interpretation is approved,
issue a separate, explicit supplemental verdict referring to the immutable v1
failure and the reason the extra criterion is deemed outside the governing
carrier decision; retain all raw data and reports byte-for-byte.

## Sign and direction check

The local correction has `c_local = c_spectral - t/q0`; therefore the inverse
spectral input is `c_spectral = c_local + t/q0`. The replay driver sets the
calculate_pl center to `c_input - center_shift`. Its preregistered spectral
case uses `center_shift = -t/q0`, so it applies `c_local + t/q0`; the local
case uses zero shift. The report defines `Delta SUMEV` as spectral minus local.
Thus the signed first-order prediction is `sum(t)`, not `-sum(t)`. The sign
and direction in the v1 literal statistic are consistent; there is no
wrong-sign issue to resolve.

This is a governance question under the existing A1.3/A2.6 process, not a
new measurement. No response rerun, tolerance change, or fit refit is
proposed.
