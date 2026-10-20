# C5 Part II — R1/R2 findings and stop decision

**R1 setup GATE: PASS.** The corrected hcp-Co Debug deck reproduced a source sizing defect before SCF: `ctnew` had extent 3 while the cutoff vector had extent 50000. `source/lattice_cluster.f90` now sizes the scaled copy to `size(this%ct)`. The Debug regression passes and emits `ntype=3`, `nbas=13`, `nrec=1`.

| Check | Result | Evidence |
|---|---|---|
| Pre-fix hcp-Co Debug reproduction | Reproduced `Array bound mismatch ... ctnew (3/50000)` at line 996 | Part II provenance, raw log SHA256 `f2d2415a…` |
| Post-fix cluster setup regression | **PASS** | `Val04RsHcpCoClusterSetup`; `ntype=3`, `nbas=13`, `nrec=1` |
| Co impurity U=0 record | **NOT RUN — invalid restart source** | Corrected deck had only fresh atomic starts; no authorized converged `_out.nml` endpoint was available. The stale-template attempt is `INVALID_SETUP`, not a failed gate. |
| No-U RS standard fixture identity | **FAIL** | Matched Release scientific outputs differ byte-for-byte from clean `8b8f7ef`. |
| No-U KS standard fixture identity | **PASS** | Scientific outputs are byte-identical; normalized `str.out` is equal after removing `cpu_s` and `TIMER REPORT` timing fields. |
| Fe rc=100 U>0 diagnostic capture | **UNACCEPTED PARTIAL** | One Debug run from accepted rc=80 checkpoint, capped at 102 steps; it did not satisfy RS acceptance. Raw findings below and in the audit sidecar. |
| Co bulk U=0 / Fe one-control matrix | **NOT RUN** | No compatible Co target-deck finite endpoint was available. Further R2 runs stopped after the no-U RS identity failure under §2.5(a). |

The Fe capture first showed an out-of-range d-shell spin occupation at state iteration 90: `P_OUT` ql spin 1 was `16.023218979776722`, above the five-state d-spin capacity; spin 2 first exceeded five at iteration 94 (`5.167394717968493`). Iteration 101 is the latest wholly finite snapshot, but it follows the occupation violation and is not an accepted state. Iteration 102 is the first snapshot with NaNs. The first observed NaN in the log is in the radial `RSEQ` diagnostic at line 5149 (`l=2`, channel 0, `E=-50 Ry`); this is an observed trajectory failure, not a proven root cause. Instrumentation coverage limits and exact hashes are in [`C5_PART_II_PROVENANCE.json`](C5_PART_II_PROVENANCE.json) and the independent [audit](C5_AUDIT.json); the raw run bundle is in [`C5_R2_assets/raw/fe_rc100_debug_capture.tar.gz`](C5_R2_assets/raw/fe_rc100_debug_capture.tar.gz).

**Decision:** STOP under C5 §2.5(a) because matched Release RS no-U scientific outputs differ byte-for-byte. R2 established neither a specific numerical defect nor a safe parameter range. No controlled full-preflight matched-environment determinism repeat was run. No physics fix was made; issue 4 remains OPEN and issue 3 remains OPEN. No further controls or Part III runs were launched.

**Interpretation (8 lines):**

- R1 confirmed an array-capacity defect and its minimal fix passes the setup regression.
- The corrected Co U=0 impurity record had no valid restart endpoint.
- The matched KS no-U fixture passes; matched RS scientific output differs from clean `8b8f7ef`.
- The Fe rc=100 capture reaches its 102-step cap without acceptance.
- It records a d-spin occupancy excursion at iteration 90 and an observed radial RSEQ NaN at iteration 102.
- Snapshot 101 is finite but follows the occupancy excursion; snapshot 102 is non-finite.
- The capture does not establish a numerical root cause or safe parameter range.
- The §2.5(a) stop halted remaining controls and Part III.

Passed: R1 setup GATE and KS no-U identity check.
Recorded: RS no-U identity failure, Fe rc=100 partial localization, and invalid/unavailable restart limitations.
Provisional: none; the R2 state is unaccepted and the fallback ladder was not run.
Stopped: C5 §2.5(a); no further campaigns, RD implementation, or merge.
Next: preserve the stop verdict and evidence for user review.
