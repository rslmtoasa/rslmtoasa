# C5.12 plan

One Luna worker executes sequentially; one batched Luna audit verifies the final gate verdicts. Sol does no per-run review.

1. Offline compare the C5.11 T2 deck Fermi value and window with the recorded KS endpoint. No run.
2. On mismatch, shift Fermi and window only, validate a strict diff and host-EF preflight, then rerun T2. On no mismatch, run one U=0 control and record whether it converges or oscillates; do not merge.
3. If the corrected T2 passes, update workflow/preflight docs, commit, audit once, then locally merge with the prescribed method and report SHA. Stop and commit records on other outcomes.

Preflight every run; identity checks single-threaded; maximum two SCF runs. Preserve Co files and accepted checkpoints. No cleanup.
