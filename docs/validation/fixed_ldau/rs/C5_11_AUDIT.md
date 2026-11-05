# C5.11 batched audit

**Verdict: T1 PASS; T2 FAIL; stop without merge or push.** HEAD `a2972286d377964f99f4e2f58d1bfb70c9f3aad2` is based on `b04bcba47375620dc6e57873b7c88eaacb0f1de9`. The recorded two-run SCF budget is fully used.

T1 evidence is consistent across the result record and committed changes: the focused impurity-U fixture is registered, 13 guard fixtures passed, no-U identity matched RS 27/27 and KS 24/24 artifacts, and the prior RD3(a) Ni auto case resolves local on its accepted window. The fixture exercises the local-host eligibility helper and checks the impurity override deck and host provenance. The source recipe remains explicitly local, as authorized; no auto-mode change was introduced into the T2 recipe.

T2 used the accepted Fe RA4/C3 `rc=80`, `ct=2*3.0d0` recipe with the accepted KS host and host/impurity U/J values of 2/0 eV. The offline record verifies the two impurity namelist edits and runtime array mapping; the KS host path and its SHA are recorded separately. Both launches report explicit local mode and no provenance error. Run 1 stopped by unexplained SIGTERM at iteration 26; run 2 restarted from that checkpoint and completed 100 iterations. The first 20-step window and final window both miss the RA4 limits, and no window in 1–100 passed. The first-window archive hash matches the result record and contains 20 trace CSVs, 20 paired state snapshots, and the capture manifest. Therefore T2's RS-accepted gate fails; endpoint deltas remain unaccepted records.

T3 is complete. The open Ni KS-host-to-RS equivalence and RD4 energy/longer-run issues remain open. No source or solver work was performed during this audit. The worktree has pre-existing untracked Co reference outputs; they are outside the C5.11 commit diff and were left untouched.

Audit checks: JSON parsed and key gate assertions passed; archive SHA and contents checked; `git diff --check b04bcba..HEAD` passed. No build, solver, cleanup, merge, or push was performed.
