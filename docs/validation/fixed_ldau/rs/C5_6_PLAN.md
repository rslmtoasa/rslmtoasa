# C5.6 plan

One Luna worker executes sequentially; one final batched Luna audit precedes merge. Two remaining campaign attempts cover T2a/T2b; required identity/unit/merge verification is separate under the C5.5 clarification. Preflight every run, full single-thread identity caps, and monitor disk throughout. Do not push.

1. Inventory CLOSED raw run directories and verify committed compact evidence before deletion. Protect live jobs, accepted KS Fe U=2, Ni RS-local and Ni U=0 checkpoints, all later references, and user Co files. Record deleted paths, allocated bytes and filesystem free-space change.
2. T2a: restart accepted KS bcc Fe U=2/J=0 in default auto mode using the current binary. Capture complete acceptance evidence and confirm the newly written host serializes local mode and KS backend. Reserve one campaign slot for T2b; no extra continuations beyond the total budget.
3. T2b: use the regenerated host for the C5.5 Fe RS impurity auto smoke. Require local mode, no provenance/mismatch error and registered RS acceptance; record impurity-minus-KS-host occupations and moment.
4. T3: complete the density-only auto candidate, provenance WARNING and its unit test; keep impurity Hubbard energy excluded. Complete C5.5 regression checks, workflow documentation and OPEN issues.
5. Audit all gates once. Adopt and execute C5.2 §7 local merge and verification only if supported. Otherwise preserve defaults, commit records and report five lines. Report the merge SHA only after successful merge verification.

## Execution outcome (2026-10-21)

T0 completed within its authorized closed-run scope: 24 directories were deleted after evidence binding, with 250,454,016 allocated bytes estimated reclaimed and 242,360 KiB observed available-space increase. The sole T2a campaign attempt failed before SCF at the legacy provenance-mode guard: native KS default auto resolved to local, while the accepted legacy U-carrying restart was interpreted as spectral. KS acceptance is unevaluated and no regenerated host exists. T2b was not run; T3 and T4/merge were stopped under the failure branch. The single batched audit passed the recorded T0 evidence and T2a failure classification. Campaign use is 1/2. No retry, retag, source modification, adoption, merge, or push occurred. See C5_6_RESULTS.md and C5_6_AUDIT.md.
