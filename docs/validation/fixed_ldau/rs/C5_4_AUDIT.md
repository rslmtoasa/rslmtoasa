# C5.4 consolidated gate audit

**Audited branch/HEAD:** `fix/ldau-rs-impurity` / `50576d81ef753e3be533bd54920eefd12ba55acb`. One consolidated audit; no SCF or source changes.

| Gate | Verdict | Audit finding |
|---|---|---|
| RD2 frozen Ni replay | **PASS** | Refit from the eight listed signed residuals gives allowance 0.002700557499 Ry; measured deviation 0.002695820417 Ry (margin 4.7371×10⁻⁶ Ry). Carrier ratios 0.9914299182/0.9914299103 and q2′ 0.1378160163/0.1378160027 Ry² pass. Sign control is admissible and separates; spin-swap predicted separation is below 10× carrier error, so it is NULL, not a failure. |
| T2 guard thresholds/tests | **PASS** | Warning and fatal boundaries match the amendment; four focused tests are reported passing. |
| T2 no-U RS/KS identity | **PASS** | Single-thread comparisons have 27 RS and 24 KS artifacts with empty mismatch lists after registered exclusions. |
| T2 KS U=0 endpoint | **UNEVALUATED** | The expected finite fatal guard stopped it at iteration 2, q0=5.1605390187; no accepted host endpoint results. |
| T2 KS U=3 endpoint | **UNEVALUATED** | Converged, but the rolling-last trace cannot establish the required 20-iteration acceptance window. |
| T3 zero-U provenance scope fix | **PASS** | Narrow zero-U/J skip and focused test; U/J correction paths retained. |
| RD3(a) Ni RS embedding | **PASS** | Accepted U=0 control window 7–26. Ratios 0.64463, 0.59821, 0.80695; all are below 10. |
| RD3(b) Ni KS-host embedding | **UNEVALUATED** | No accepted T2 KS host exists. |
| RD4 energy | **UNEVALUATED** | Debug W1 series was not measured, so the decision table was not applied. The RelWithDebInfo λ=0 crash replay is not an energy measurement; energy routing stays off. |
| RD5 no-U identity / unchanged bulk behavior | **PASS** | No-U identity passes; prior unchanged-bulk evidence is referenced by the RD5 record. |
| RD5 explicit-spectral identity | **PASS** | The frozen old/new hash maps have equal 27-file name sets and every SHA256 matches; recorded thread caps are single-thread and both runs exited 0. |
| RD5 RS impurity `auto` | **UNEVALUATED** | No accepted-endpoint `auto` reproduction recorded. |
| RD5 KS impurity `auto` | **UNEVALUATED** | No accepted KS endpoint exists for the required `auto` reproduction. |
| C5.2 §7 merge verification | **UNEVALUATED** | No merge or prescribed merged-tree checks were performed. |

The absent accepted KS endpoint does not waive the KS `auto` gate. C5.4 T5 applies because RD2 and RD3(a) pass, and it requires completing the RD5 `auto` cases before §7 merge. Record the missing KS case as **UNEVALUATED** and continue under C5.3 item 2; do not infer a pass or change criteria. RD5 adoption is **NOT ADOPTED** and merge is **not permitted** on this audit record. Energy routing remains off.

The T2–T4 append-only ledger has six entries against the eight-run cap; T5 separately records two verification invocations, for eight total. The auditor launched no SCF. Re-running the offline T1 driver was blocked before calculation because the ambient Python lacks `f90nml`; the cubic fit and control classifications were independently recomputed from the frozen eight-point values in `C5_4_T1_RD2.json`.

Five-line disposition:
1. Passed: RD2, threshold tests, no-U identities, the zero-U scope fix, RD3(a), and explicit-spectral identity.
2. Recorded: U=0 KS finite fatal at the revised guard and the RelWithDebInfo first-recursion SIGSEGV replay.
3. Unevaluated: U=3 KS acceptance, RD3(b), RD4 energy, and both required `auto` endpoint cases.
4. Stopped: adoption and §7 merge; no criteria were waived.
5. Next: preserve records and continue under C5.3 item 2 when an accepted KS endpoint and remaining required RD5 evidence are available.
