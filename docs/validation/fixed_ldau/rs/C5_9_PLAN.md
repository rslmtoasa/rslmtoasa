# C5.9 plan

One Luna worker executes these cards sequentially; one batched Luna audit checks final gate verdicts. Sol collects card outputs without per-run review.

1. Evaluate existing regenerated KS Fe state against C5.9 convergence, endpoint tolerances, and written provenance. No solver run. Stop on a failed endpoint comparison; apply C5.3 to uncovered outcomes.
2. If accepted, prepare the density-only auto candidate on the fixed source and run the RS Fe impurity smoke test. Preflight before every invocation, preserve full acceptance evidence, and use at most two campaign runs including any allowed continuation. Monitor disk before, during, and after runs.
3. Complete single-thread identity gates, fixtures, workflow and known-issue documentation. Rerun the exact RD4 lambda=0 Release diagnostic once with a validated first-iteration stop; record crash status only. Verification invocations remain separate from campaign budget under the approved clarification.
4. Commit compact evidence and request one final audit. Only if required gates pass, perform the C5.2 section 7 local merge and verification; report its SHA without pushing. Otherwise commit records and report the prescribed stop outcome.

Preserve Co files and accepted checkpoints; do not touch cleanup/ or tests/gbt/. No additional deletion is planned. Starting free space: 9,852,724 KiB.
