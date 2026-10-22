# C5.11 plan

One Luna worker applies the bound candidate, executes cards sequentially, and reports compact evidence. Sol performs one batched audit at the end.

1. Run the 13 guard fixtures, new impurity-U local-mode fixture, RD3(a) Ni auto case, and single-thread no-U RS/KS identity checks.
2. Build the Fe impurity deck by copying the accepted RA4/C3 Fe recipe; assert its diff changes only host potential and U/J. Preflight every run. Use at most two SCF runs and retain the full RS 20-iteration acceptance evidence.
3. Complete docs and known-issue closure, commit evidence, and request the single audit. If T2 fails, stop after committing records. If all gates pass, merge locally to fable_v4b_u with --no-ff, verify, report SHA, and do not push.

Preserve Co files and accepted checkpoints. No cleanup or unrelated worktree changes.
