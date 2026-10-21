# C5.4 final status — sole audit completed

The C5.2 frozen RD2 report remains unchanged as a historical result. C5.4's required offline eight-point per-state A2.2 fit and corrected scalar-control classification supersede its RD2 verdict with **audited PASS**. See [C5.4 T1](C5_4_T1_RD2.json), the [sole consolidated audit](C5_4_AUDIT.md), and the historical [RD2 report](C5_2_PART_D_assets/ni_u3j07_rd2_replay/RD2_REPORT.md).

| Gate | Worker finding | Evidence |
|---|---|---|
| RD2 | **PASS.** Carrier ratios are 0.99142992 per spin; q2′ is nonnegative; opposite-spin control is NULL and sign-wrong control is admissible. | [T1 amended fit and controls](C5_4_T1_RD2.json) |
| T2 guard and no-U identities | **PASS.** Guard thresholds/non-finite tests and single-thread RS/KS no-U identity pass. | [T2](C5_4_T2.json) |
| T2 Ni KS endpoints | **UNEVALUATED.** U0 hit the finite guard at iteration 2; U3 converged but lacks the raw 20-iteration acceptance window. | [T2](C5_4_T2.json) |
| RD3(a) | **PASS.** Ni d-up/d-down/sphere-moment ratios are 0.64463/0.59821/0.80695. | [T3](C5_4_T3_RD3.json) |
| RD3(b) | **UNEVALUATED.** T2 produced no accepted Ni KS U0/U3 endpoints. | [T3](C5_4_T3_RD3.json) |
| RD4 | **UNEVALUATED / NOT MEASURED.** The RelWithDebInfo replay reproduced the lambda-zero first-recursion SIGSEGV; the Debug W1 series was not run. Energy routing remains **OFF**. | [T4](C5_4_T4_RD4.json) |
| RD5 | No-U RS/KS identity and bulk unchanged checks **PASS**; explicit-spectral identity **PASS** on matched one-step runs; both required impurity `auto` endpoint checks **UNEVALUATED**; no adoption. | [T5 gate record](C5_4_T5_RD5.json), [spectral hashes](C5_4_T5_spectral_identity_hashes.json) |

The T2–T4 campaign used 6 of its 8 allowed solver invocations. T5 used two separately recorded one-step explicit-spectral verification invocations, bringing the actual total to **8 SCF solver invocations**. T5 did not add a host campaign. The auditor launched no SCF. No further runs are authorized.

The accepted KS endpoint required for RD3(b) and KS-local impurity `auto` reproduction is missing. The RS-local impurity `auto` reproduction is also unevaluated. Production `auto` eligibility remains bulk-only; no adoption or energy/default change was made. These are C5.3 uncovered outcomes, so there is **no merge**. The audit reviewed HEAD `50576d81ef753e3be533bd54920eefd12ba55acb`; the local `fable_v4b_u` branch is at `908074f64b4d3143d2d91c0124fc79d126b8dbea`. This final commit records status only. No push occurred.

Final audit verdict: **PASS for the evaluated gates; NOT ADOPTED / NO MERGE because required RD5 auto endpoints remain UNEVALUATED**. The audit independently recomputed the RD2 fit and T5 output-map equality; it did not rerun SCF. The ambient-Python re-execution of the T1 scan was unavailable to the auditor because `f90nml` was missing there; the audit recomputed the fit and control classifications independently from the frozen eight-point values.
