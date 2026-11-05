# C5.8 results

## Outcome

The Release crash in the fixed-U trace replay was caused by an unassociated `pl_bands%control` pointer. `replay_fixed_u_m0_moments` initialized `pl_bands%symbolic_atom` and `pl_bands%lattice`, then called `calculate_pl`; the latter reads `this%control%nsp` at `source/bands.f90:2254`. The faulting load was at address `0x14` from the null control object. The minimal correction binds `pl_bands%control => this%control` at `source/self.f90:1431`. This changes diagnostic trace replay setup only; it does not change the Hubbard or energy formulas.

The diagnosis is supported by matching Release, checked RelWithDebInfo, and Debug `.ips` reports, the Debug DWARF line mapping to `bands.f90:2254`, and disassembly showing a load from `pl_bands%control` followed by offset `0x14`. The original crashes occurred only with the fixed-U trace enabled; untraced `nstep=1` runs passed. The older 908074f run was rejected by its strict mode guard and is invalid for crash comparison. The 8b8f7ef untraced `nstep=1` run passed; its exact traced `nstep=240` launch completed naturally because its watcher did not stop at the first iteration. Neither historical case establishes a pre-existing production crash.

## Reproduction ledger

Eleven short invocations were used out of the twelve-short cap. Stop-control deviations are reported as such, not counted as first-iteration passes.

| # | Purpose and configuration | Result |
|---|---|---|
| 1 | Initial explicit-local exploratory deck; not the frozen C5.7 target | One iteration completed; invalid target for T1 comparison |
| 2 | Current Release, exact frozen T2a input/seed, `nstep=1`, trace disabled | First iteration passed |
| 3 | Checked RelWithDebInfo with requested `-O2 -g -fbacktrace -fcheck=bounds,pointer -finit-real=snan -finit-integer=-999999`, same short deck, trace disabled | First iteration passed; no check failure |
| 4 | Current Debug, same short deck, trace disabled | First iteration passed; actual Debug flags recorded in build evidence (no bounds/pointer checks enabled) |
| 5 | Historical Release 908074f | Invalid comparability: strict fixed-U provenance guard rejected the legacy seed; not a crash/pass result |
| 6 | Historical Release 8b8f7ef, `nstep=1`, trace disabled | First iteration passed |
| 7 | Current Release, original `nstep=240` plus exact C5.7 trace environment | SIGSEGV; first-iteration stop was not reached |
| 8 | Checked RelWithDebInfo, original `nstep=240` plus exact trace environment | SIGSEGV; no bounds/pointer/SNAN diagnostic preceded it |
| 9 | Current Debug, original `nstep=240` plus exact trace environment | SIGSEGV; `.ips` localizes the same trace replay access |
| 10 | Historical Release 8b8f7ef, original `nstep=240` plus exact trace environment | Watcher did not stop; natural convergence after three iterations. Stop-control deviation; not a one-iteration pass. The 908074f rejection was not retried. |
| 11 | Fresh Release rebuilt from the exact private W1 diagnostic source plus only the pointer fix; original Ni lambda-zero deck, `nstep=200`, trace and external first-complete-iteration stop | Watcher captured a complete first iteration; external SIGTERM, exit 143. No crash through that iteration. This is a diagnostic regression only, with no energy conclusion. |

The one watcher-format validation used an isolated copy of an existing complete trace and paired state; it captured one stable iteration successfully. It did not invoke the solver. All actual solver launches used one-thread caps and fresh writable run copies. T1 short artifacts and source/build records remain under `/private/tmp/c57-*` and `/private/tmp/c58-*`; the RD4 regression is `/private/tmp/c58-rd4-fixed-release-lambda0`.

## Fixed-source verification

- Fixed current Release binary: SHA-256 `266e06035d59a1681850e0790bdf6bd04eabf3990d40d0622e17475b15c7c7ce`; effective Fortran flags include `-O3 -fbacktrace -g -g` with the project OpenMP flags. Source `source/self.f90` SHA-256 `4886c9e7eb5a3534e1e9435093490e5d26a9b3f0ad88c0dfe9dc64e037b2adab`.
- The fixed current Release T2a run used the exact frozen deck SHA-256 `1af5bcc4c3633f238043438e0aaa59575c8854f8fdbbd927af52f5660ae04083` and accepted seed SHA-256 `bb647254a5612f53d00e2745026c2cde32511f4363a1c7c3a1331d06573ead6b`. It completed three outer iterations, with two `Not converged` observations and a final `Converged` message; `Fe_scf.nml` and `Fe_out.nml` both recorded `hubbard_radial_moments='local'` and `hubbard_backend='ks'` (SHA-256 `2284194179736d2e6694dd40514b55578bdf2558e9826d5d7ba9d5c6cc13b4bd`). The trace has one final-iteration observation, short of the required 20-observation rolling window. Classify T2a as converged but acceptance unevaluated/uncovered under C5.3; no T2b, adoption, or merge was performed.
- Post-fix fixed Release gates pass: all 13 RA2 guard fixtures; no-U RS identity (27 artifacts) and KS identity (24 artifacts), compared against the frozen C5.2 baselines with the registered `cpu_s` normalization in `str.out`. Records: `/private/tmp/c58-post-fix-fixtures/RA2_provenance.json` and `/private/tmp/c58-post-fix-identity/C5_8_post_fix_no_u_identity.json`. These checks are separate from the short/campaign counts.
- Fixed-source RD4 diagnostic Release binary SHA-256 `d87479d402f959b2fe80975fdb8a25afa2a1794ecae1f58146b8c5f6cbf26f7e`; actual flags `-O3 -fbacktrace -g -g`. It uses the exact private C52/C54 W1 source tree with only the pointer binding added; it is not the current production source branch. The five RD4 input files retain the original T4 hashes. Preflight passed against the W1 schema. The external watcher stopped after its first complete state/trace pair. The run trace SHA-256 is `589f60ccf4a720cbc9cb6a348cbb8cf32a5dc2602362cb0aafe49d99ec82ada3`; `testrun.log` SHA-256 is `4ec5e8e0e767cd82c649bd321ef5d0cd029f5b4855909314de2775ae7b3a0494`.

## Budgets and disposition

Short reproductions: 11/12. SCF campaigns: 1/3 (fixed T2a only). Required fixtures and identity checks were separate verification invocations. Disk free space for T2a was 11,949,113,344 bytes before launch and 11,659,476,992 bytes afterward; RD4 short-run free space was 10,569,363,456 bytes at launch and 10,421,338,112 bytes at completion, above the 1 GiB guard floor throughout. No cleanup was performed.

C5.8 stops with T2a acceptance uncovered because the required 20-observation window was unavailable. No T2b/T3/T4 work, endpoint adoption, merge, push, or additional solver run was performed. Impurity energy routing remains off. User Co references and unrelated GBT jobs were left untouched.
