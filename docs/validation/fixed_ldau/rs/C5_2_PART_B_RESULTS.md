# C5.2 Part B — guard, robustness, and recursion records

**B1: PASS.** The collinear `nsp=1` SCF path now guards per-spin q0 capacity and rejects non-finite active q0/q2/PL plus new unmixed q1 values with atom, channel, iteration, and value. Three focused synthetic tests pass; `nsp>=2` and inactive padding are unchanged.

| Check | Result | Evidence |
|---|---|---|
| Candidate RS no-U identity | **PASS** | Scientific outputs byte-identical to the C5.2 A1 baseline; only `str.out` cpu_s differs. |
| Candidate KS no-U identity | **PASS** | Scientific outputs byte-identical to clean-8b; only `str.out` cpu_s differs. |
| B2(a) Fe rc=100 Chebyshev | **ACCEPTED AT STEP 59** | Diff `1.40e-8` Ry; ql d-spin spreads `2.0146e-7` / `4.0140e-7`; 20 printed physical-moment values identical at F10.6, bounding true spread to `<=1e-6 μB` (not zero). The run was intentionally stopped after preserving the accepted endpoint; all owned processes are confirmed gone. |
| B2(b) Co block recursion | **FINITE, UNACCEPTED** | Normal exit at nstep=25 cap, final Diff `0.0054412342` Ry, no NaN/Inf; 19 valid checkpoints, short of a complete 20-observation acceptance window. |
| Co continuation ladder | **OMITTED BY USER CLARIFICATION** | Six-run cap retained; Co remains unaccepted and is not adopted. |
| Part B decision | **RECORDED** | Fe B2(a) accepted; Co B2(b) remains finite but unaccepted; KNOWN_ISSUES entries 1, 2, and 4 updated; entry 3 unchanged. |

The user’s exact Co ladder question and reply, run manifests, full hashes, and compact raw records are preserved in the provenance sidecar and `C5_2_PART_B_assets/`. Identity comparisons use the C5.3 single-thread protocol. B1 build, source manifest, and test details are recorded there.

**Interpretation (8 lines):**

- B1 tests confirm the over-capacity and non-finite production guard paths, plus the permitted boundary case.
- RS and KS no-U identity checks pass for scientific outputs; only recorded timing metadata differs.
- Co completed normally and remained finite, but did not converge or complete an acceptance window by its cap.
- The user directed that the six-run budget remain fixed and the Co continuation ladder be omitted.
- Co therefore remains unaccepted and is not adopted.
- Fe rc=100 Chebyshev has an accepted step59 endpoint. The run was intentionally stopped after preserving that endpoint; the later completed step65 and interrupted step66 outputs were not used, and all owned processes are confirmed gone.
- No Part B production default or parameter is adopted based on the Co result.
- B1, B2, and B3 records are complete; continue to Part C under C5.3.

Passed: B1 regression tests, audited RS/KS identity checks, and Fe B2(a) acceptance at step59.
Recorded: finite but unaccepted Co cap result and the user's budget clarification.
Provisional: none.
Stopped: none under C5.2 §2.5. The continuing Fe process was intentionally terminated after preserving its accepted step59 endpoint; all owned processes are confirmed gone.
Next: continue to Part C under C5.3.
