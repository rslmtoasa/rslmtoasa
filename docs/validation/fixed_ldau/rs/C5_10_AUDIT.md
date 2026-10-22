# C5.10 batched audit

**Verdict: T1 accepted; T2 failed and exhausted its 2/2 campaign-start budget. No impurity endpoint was accepted, no production adoption or merge is authorized by these results, and T3 correctly remains stopped.** The failures are preparation/input errors, not evidence of a physical failure of the impurity method.

- T1's old-EF comparison is N/A under C5.10. The new endpoint value `-0.082949627451487762 Ry` is present at `/private/tmp/c58-t2a-fixed-release/run/fixed_u_trace.csv` (`fermi_level`, endpoint row); the prose's `-0.0829496275 Ry` is a rounded rendering. The ordinary log/report displays only `-0.082950 Ry`.
- T2 start 1 (`run2`) has malformed Fortran namelist indices `hubbard_u_impurity(:,1)` and `hubbard_j_impurity(:,1)` in `input.nml`; the namelist reader rejects it before any SCF iteration. This is a generated-input preparation error, not a solver/physics result.
- T2 start 2 (`run4`) has correct `(1,3)` U/J indices but sets `ct(:)=6.0,6.0`; the cited Fe rc=80 recipe requires `3.0` per type (`2*3.0d0`). The log contains repeated missing-neighbour-vector diagnostics. Three complete iterations report Diff `0.1967289821`, `0.2034926014`, and `0.2010702253`; the next iteration ends with invalid local d-shell moment correction and an occupation-capacity warning. Thus the worker's “stopped after three iterations” means three completed iterations, not a clean external stop at iteration 3. No 20-iteration window exists. This trajectory is invalid-setup evidence and must not be described as a physical instability or accepted endpoint.
- Both starts count toward the authorized cap; there is no remaining retry or continuation. T2 fails; candidate remains private/diagnostic, and T3 gates, adoption, and merge are not completed. The recorded candidate unit test is not a T3 identity gate.
- RD4 bounded Release regression is a separate short invocation: the log shows local mode, one captured complete iteration, and exit 143 after external SIGTERM; it supports “no crash through iteration 1” only. No energy measurement occurred.
- The reviewed production diff against the recorded restored source binding is empty for candidate production paths. Current worktree status contains only the listed untracked user Co reference outputs; they are outside this audit and were not modified.

### Five-line disposition

1. Passed: T1 KS Fe host accepted under C5.10; missing prior EF is N/A.
2. Recorded: precise new EF source and one bounded RD4 iteration without crash; no RD4 energy.
3. Provisional: no T2 impurity endpoint; the private auto candidate remains diagnostic only.
4. Stopped: T2's two starts were consumed by malformed indices and incorrect cluster geometry; T3, adoption, and merge not performed.
5. Next: a separately authorized run may correct deck generation and geometry before using any SCF budget; no campaign retry is authorized by C5.10.
