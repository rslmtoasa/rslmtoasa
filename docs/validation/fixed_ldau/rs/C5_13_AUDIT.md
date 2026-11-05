# C5.13 final audit

**Disposition: PASS.** Reviewed the C5.13 amendment, plan/results, workflow and known-issue documentation, merge metadata, and existing verification records. No solver or tests were run during this audit.

The C5.13 amendment contains its user-approved criteria and outcome. Existing C5.12 T2 run 2, iterations 81–100, satisfies the C5.13 smoke limits: Diff 5.71e-8 ≤ 1e-7; d spreads 1.9525747685e-7 and 8.2702854765e-7 ≤ 1e-5; moment spread 7.662924717e-7 μB ≤ 1e-5 μB. It completed 100 iterations with return code 0, explicit-local mode, and no provenance or mode errors. The recorded host-EF preflight has zero mismatch, within 1e-6 Ry.

The workflow guide covers RS-local and KS-local hosts, both general and impurity U/J input styles, the `auto` rule, legacy-host restart guidance, impurity energies excluding the Hubbard term, and setting `fermi` and the window to host endpoint values. `KNOWN_ISSUES.md` closes the impurity-auto mode-selection entry while retaining the distinct KS-host-to-RS-impurity Ni equivalence issue as OPEN.

The recorded local `--no-ff` merge SHA is `fe4ae42132b3cf726a3ed8e9827c567884637723`, with base parent `b04bcba47375620dc6e57873b7c88eaacb0f1de9` and second parent `f49113552bbaad66e28e3bd3093046f0b4b8d2be`. The merged-tree evidence records single-thread no-U identity PASS (RS 27/27, KS 24/24), all 13 RA2 guard scenarios, and `UnitAutoImpurityLocalMode` PASS (1/1). Both identity and RA2 provenance name the merge SHA as source commit.

The merge diff contains no `example/bulk/hcpCo-folder/hcpCo_ref` files. Those Co paths are currently untracked in the worktree; this audit did not modify them. C5.13 results state no new SCF and no push. C5.3 item 2 remains uncovered, as recorded in the amendment.

Evidence reviewed: `C5_13_PLAN.md`, `C5_13_RESULTS.md`, `C5_13_RESULTS.json`, `C5_AMENDMENT.md`, `C5_12_assets/run2_preflight.json`, `lda_u_workflow.rst`, `tests/KNOWN_ISSUES.md`, `/private/tmp/c513_merged_identity.json`, and `/private/tmp/c513_merged_ra2/RA2_provenance.json`.
