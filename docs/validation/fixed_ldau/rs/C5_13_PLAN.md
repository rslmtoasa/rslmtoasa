# C5.13 plan

One Luna worker evaluates the existing run, completes docs, merges if the smoke gate passes, and runs requested checks on the merged tree. Sol performs one batched audit at the end.

1. Evaluate C5.12 final window against C5.13 mode, provenance, Diff and spread limits; record impurity-minus-KS-bulk d and moment.
2. If it passes, complete workflow/known-issue docs, commit, merge locally with --no-ff, then run single-thread RS/KS identity, 13 guard fixtures and auto impurity fixture on the merged tree. No SCF runs or pushes.
3. If it fails, record and stop without merge. Apply C5.3 to uncovered outcomes. Preserve Co files and accepted checkpoints; no cleanup.
