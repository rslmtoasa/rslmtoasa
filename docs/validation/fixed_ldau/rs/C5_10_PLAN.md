# C5.10 plan

One Luna worker executes sequentially; one batched Luna audit checks the completed gate packet before the local merge. Sol collects milestones without per-run review.

1. Record the KS Fe host as accepted under C5.10, mark the unavailable old Fermi reference N/A, and bind the new Fermi reference to its actual endpoint evidence.
2. Prepare the density-only impurity auto candidate on the fixed source. Preflight and run the RS Fe impurity on the accepted KS-local host; retain the required 20-iteration RS acceptance evidence. Use at most two campaign runs including continuation; monitor disk before, during, and after.
3. Complete single-thread identity gates and fixtures, workflow and known-issue documentation, and the required bounded RD4 Release crash regression without an energy measurement. Apply N/A only to genuinely absent reference data, retaining all available quantitative gates.
4. Commit compact evidence, obtain the single final audit, then perform C5.2 section 7 local merge and required checks if gates pass. Report the merge SHA without pushing. Record uncovered outcomes under C5.3 without adopting changes.

Preserve user Co files and accepted checkpoints. Do not modify cleanup/ or tests/gbt/. No additional deletion is planned.
