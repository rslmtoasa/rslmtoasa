# RS-LMTO-ASA decisions

This file records one line per closed campaign. Full evidence remains
available at the `lr-campaign-archive` tag, at the path listed below.

| Campaign / file | Date | Conclusion (≤ 25 words) | Path at tag |
|---|---|---|---|

### ACC

| Campaign / file | Date | Conclusion (≤ 25 words) | Path at tag |
|---|---|---|---|
| ACC-00 CPU baseline | not stated | CPU baseline recorded for the accelerator campaign. | `lr-campaign-archive:docs/dev/ACC-00_CPU_BASELINE.md` |
| ACC-01 RS CUDA coverage | not stated | Typed CUDA ABI and low-level route comparisons passed. | `lr-campaign-archive:docs/dev/ACC-01_RS_CUDA_COVERAGE.md` |
| ACC-06 CPU/GPU crossover | 2026-08-18 | CPU/GPU crossover measurements were recorded with the typed LAPACK backend. | `lr-campaign-archive:docs/dev/ACC-06_CPU_GPU_CROSSOVER.md` |
| ACC-07 H(k) materialization | 2026-08-19 | Decision A retained host-side H(k) assembly; no production-code change was made. | `lr-campaign-archive:docs/dev/ACC-07_HK_MATERIALIZATION.md` |
| ACC-09 reciprocal variants | not stated | Reciprocal CUDA operator variants were audited and documented. | `lr-campaign-archive:docs/dev/ACC-09_RECIPROCAL_VARIANTS.md` |
| ACC-10 Lehmann CUDA | not stated | Strict-Lehmann CUDA contraction passed the complete comparison. | `lr-campaign-archive:docs/dev/ACC-10_LEHMANN_CUDA.md` |
| ACC-13 KPM transport CUDA | not stated | Existing RS CUDA moment transport was complete; no new GPU kernel was needed. | `lr-campaign-archive:docs/dev/ACC-13_KPM_TRANSPORT_CUDA.md` |
| ACC-P0 persistent real-material benchmarks | not stated | Persistent real-material benchmark workflow and evidence were consolidated. | `lr-campaign-archive:docs/dev/ACC-P0_PERSISTENT_REAL_MATERIAL_BENCHMARKS.md` |
| ACC-P0 Fe supercell benchmarks | not stated | report; no recorded verdict | `lr-campaign-archive:docs/dev/ACC-P0_SUPERCELL_FE_BENCHMARKS.md` |
| ACC-P2 real-material CUDA cleanup | not stated | The completed local probe found no further production cleanup requirement. | `lr-campaign-archive:docs/dev/ACC-P2_REAL_MATERIAL_CUDA_CLEANUP.md` |
| Phase III-A accelerator blueprint | not stated | Accelerator phase defined validated CPU references and a narrow reciprocal target. | `lr-campaign-archive:docs/dev/RS_LMTO_ASA_PHASE_III_A_ACCELERATOR_BLUEPRINT.md` |
| Phase III-A current steering | not stated | SCF-B0C and SCF-B0C-RS were complete; SCF-B1 was next. | `lr-campaign-archive:docs/dev/RS_LMTO_ASA_PHASE_III_A_CURRENT_STEERING.md` |
| ACC performance rescue | not stated | Performance rescue separated GPU startup overhead from steady-state reciprocal work. | `lr-campaign-archive:docs/dev/RS_LMTO_ASA_ACC_PERFORMANCE_RESCUE.md` |

### B2 / RF05

| Campaign / file | Date | Conclusion (≤ 25 words) | Path at tag |
|---|---|---|---|
| B2 reciprocal Green handover | 2026-07-16 | B2.6 completed; intersite normalization was resolved and both gates were signed. | `lr-campaign-archive:docs/dev/B2_RECIPROCAL_GREEN_HANDOVER.md` |
| RF05 reciprocal execution backends | not stated | Typed reciprocal execution backend contract and LAPACK factory were established. | `lr-campaign-archive:docs/dev/RF05_RECIPROCAL_EXECUTION_BACKENDS.md` |

### GBT

| Campaign / file | Date | Conclusion (≤ 25 words) | Path at tag |
|---|---|---|---|
| GBT completion blueprint | 2026-08-03 | The branch was classified as prototype evidence, not a completed GBT solution. | `lr-campaign-archive:GBT_RS_LMTO_completion_blueprint.md` |
| GBT WP0 / G0 | 2026-08-03 | q=0 ordinary references remained the only golden values; finite-q results were diagnostic. | `lr-campaign-archive:docs/dev/GBT_WP0_G0_REPORT.md` |
| GBT WP1 / G1 | 2026-08-03 | Independent algebraic acceptance oracles were defined without changing production routing. | `lr-campaign-archive:docs/dev/GBT_WP1_G1_REPORT.md` |
| GBT WP2 / G2 | 2026-08-03 | G2E passed, while the production operator slice G2O remained open and failed. | `lr-campaign-archive:docs/dev/GBT_WP2_G2_REPORT.md` |
| GBT WP3/WP4 gates | 2026-08-03 | Representation split and first-order primitive-S linking were implemented; later deletion work had not started. | `lr-campaign-archive:docs/dev/GBT_WP3_WP4_GATES_REPORT.md` |
| GBT WP5 / G5 | 2026-08-03 | Shared linked-S operator passed with accepted recursion bound `lld <= 28`. | `lr-campaign-archive:docs/dev/GBT_WP5_G5_REPORT.md` |
| GBT WP6a HOH/overlap | 2026-08-03 | HOH support passed for the supported orthogonal path; overall G6 remained open. | `lr-campaign-archive:docs/dev/GBT_WP6A_HOH_OVERLAP_REPORT.md` |
| GBT WP6b CCOR | 2026-08-03 | The audited CCOR slice passed; final G6 remained open pending other terms. | `lr-campaign-archive:docs/dev/GBT_WP6B_CCOR_REPORT.md` |
| GBT WP7 / G7 | 2026-08-06 | Gate G7 passed with a common rotating-frame density contract and independently converged route agreement. | `lr-campaign-archive:docs/dev/GBT_WP7_G7_REPORT.md` |
| GBT WP8 / G8 | 2026-08-06 | Gate G8 passed with full-BZ default and q-aware mesh-cache invalidation. | `lr-campaign-archive:docs/dev/GBT_WP8_G8_REPORT.md` |
| GBT WP9 / G9 | 2026-08-07 | Gate G9 failed; algebraic checks passed but cone-angle scaling and band-energy residuals remained open. | `lr-campaign-archive:docs/dev/GBT_WP9_G9_REPORT.md` |
| GBT WP10 / G10 | 2026-08-08 | G9 did not pass, so WP10 could not establish final closure. | `lr-campaign-archive:docs/dev/GBT_WP10_G10_REPORT.md` |

### GBT_REAUDIT

| Campaign / file | Date | Conclusion (≤ 25 words) | Path at tag |
|---|---|---|---|
| GBT re-audit WP00 | 2026-08-27 | Architecture was established, but production readiness and several physics gates remained unvalidated. | `lr-campaign-archive:docs/dev/GBT_REAUDIT_WP00_CURRENT_STATE.md` |
| GBT re-audit WP01 | not stated | Scalar-relativistic fixed-potential q=0 operator gate passed for the tested scope. | `lr-campaign-archive:docs/dev/GBT_REAUDIT_WP01_Q0.md` |
| GBT re-audit WP02 | not stated | Scalar-relativistic fixed-potential operator gate passed for the tested scope. | `lr-campaign-archive:docs/dev/GBT_REAUDIT_WP02_GAUGE_SHIFTEDK.md` |
| GBT re-audit WP03 | 2026-08-27 | Audited composite terms passed; CCOR, SOC, Hubbard-V, local-axis, and constrained terms remained outside scope. | `lr-campaign-archive:docs/dev/GBT_REAUDIT_WP03_COMPOSITE_COVARIANCE.md` |
| GBT re-audit WP04 | not stated | Fixed-potential commensurate-supercell operator oracle passed without claiming production or SCF equivalence. | `lr-campaign-archive:docs/dev/GBT_REAUDIT_WP04_COMMENSURATE.md` |
| GBT re-audit WP05 | not stated | Rotating-frame and lab-frame semantics were defined; SCF and controller tuning were out of scope. | `lr-campaign-archive:docs/dev/GBT_REAUDIT_WP05_FRAME_CONTRACT.md` |
| GBT re-audit WP06 | not stated | Constraint-field covariance was audited with controller penalty kept separate from DFT energy. | `lr-campaign-archive:docs/dev/GBT_REAUDIT_WP06_CONSTRAINT_FIELD.md` |
| GBT re-audit WP07 | not stated | Fixed-state constraint-energy semantics were recorded without a material constrained-SCF claim. | `lr-campaign-archive:docs/dev/GBT_REAUDIT_WP07_CONSTRAINT_ENERGY.md` |
| GBT re-audit WP08 | not stated | Corrected constrained MFT was implemented for the acoustic branch while existing modes were preserved. | `lr-campaign-archive:docs/dev/GBT_REAUDIT_WP08_CORRECTED_MFT.md` |
| GBT re-audit WP09 | 2026-08-27 | Instrumentation was delivered, but a converged physical WP09 result was not established. | `lr-campaign-archive:docs/dev/GBT_REAUDIT_WP09_HARMONIC_KGRID.md` |
| GBT re-audit WP10 | 2026-08-27 | q-reversal symmetry passed; converged small-q curvature was not established. | `lr-campaign-archive:docs/dev/GBT_REAUDIT_WP10_SMALL_Q.md` |
| GBT re-audit WP11 | 2026-08-27 | G11 was not qualified because small-q curvature and interaction-range evidence remained insufficient. | `lr-campaign-archive:docs/dev/GBT_REAUDIT_WP11_LKAG_JQ.md` |
| GBT re-audit WP12 | 2026-08-27 | Scoped q=0 SCF closure passed; finite-q constrained-SCF closure remained unqualified. | `lr-campaign-archive:docs/dev/GBT_REAUDIT_WP12_SCF_CLOSEOUT.md` |

### KPM

| Campaign / file | Date | Conclusion (≤ 25 words) | Path at tag |
|---|---|---|---|
| KPM B0C report | not stated | Benchmark harness requirements were consolidated into a canonical schema. | `lr-campaign-archive:docs/dev/RS_LMTO_ASA_KPM_B0C_REPORT.md` |
| KPM B1 report | not stated | Published CPU/GPU rows passed their declared production-output correctness checks. | `lr-campaign-archive:docs/dev/RS_LMTO_ASA_KPM_B1_REPORT.md` |
| KPM G1.2 report | 2026-08-20 | Exclusive timer contract and separate CPU/CUDA campaign records were established. | `lr-campaign-archive:docs/dev/RS_LMTO_ASA_KPM_G1_2_REPORT.md` |
| KPM G1.3 report | not stated | CUDA reconstruction used resident diagonal-packed matrices and tiled device-side generation. | `lr-campaign-archive:docs/dev/RS_LMTO_ASA_KPM_G1_3_REPORT.md` |
| KPM G1.4 report | 2026-08-21 | Optimized CUDA transport floor was decomposed and the remaining overhead was explained. | `lr-campaign-archive:docs/dev/RS_LMTO_ASA_KPM_G1_4_REPORT.md` |
| KPM transport GPU follow-up | not stated | Optimization campaign closed and the fair harness was consolidated in B0C. | `lr-campaign-archive:docs/dev/RS_LMTO_ASA_KPM_TRANSPORT_GPU_FOLLOWUP.md` |

### SCF

| Campaign / file | Date | Conclusion (≤ 25 words) | Path at tag |
|---|---|---|---|
| SCF-B0C report | not stated | Shared CPU/GPU SCF benchmark harness and canonical result package were added. | `lr-campaign-archive:docs/dev/RS_LMTO_ASA_SCF_B0C_REPORT.md` |
| SCF-B0C-RS report | not stated | Real-space routes joined the shared SCF benchmark harness with route metadata and profile closure. | `lr-campaign-archive:docs/dev/RS_LMTO_ASA_SCF_B0C_RS_REPORT.md` |
| SCF-B1R2 report | not stated | Lean desktop tier completed 18/18 cases with PASS row status. | `lr-campaign-archive:docs/dev/RS_LMTO_ASA_SCF_B1R2_REPORT.md` |
| SCF-B1R report | not stated | Scoped performance campaign completed, but universal RS-vs-k-space accuracy was inconclusive. | `lr-campaign-archive:docs/dev/RS_LMTO_ASA_SCF_B1R_REPORT.md` |
| SCF-B1 report | not stated | Scoped performance campaign completed; common-potential RS-vs-k-space accuracy remained inconclusive. | `lr-campaign-archive:docs/dev/RS_LMTO_ASA_SCF_B1_REPORT.md` |

### PHASE

| Campaign / file | Date | Conclusion (≤ 25 words) | Path at tag |
|---|---|---|---|
| Phase 3 status | not stated | Feature status by blueprint item and test coverage were recorded. | `lr-campaign-archive:docs/dev/PHASE3_STATUS.md` |
| Phase II validation | not stated | Validation maturity remained scoped; no broad GBT or TDDFT promotion was made. | `lr-campaign-archive:docs/dev/PHASE_II_VALIDATION.md` |
| Phase I stabilization | not stated | Phase-I structural stabilization record was complete for the stabilized architecture. | `lr-campaign-archive:docs/dev/PHASE_I_STABILIZATION.md` |

### REFACTORING

| Campaign / file | Date | Conclusion (≤ 25 words) | Path at tag |
|---|---|---|---|
| Refactoring plan | not stated | Structural refactoring plan established no-physics-change and regression-preservation rules. | `lr-campaign-archive:docs/dev/REFACTORING_PLAN.md` |
| Review notes T1–T6 | 2026-07-03 | Independent clean build and regression matrix verified 8/8 runnable cases. | `lr-campaign-archive:docs/dev/REVIEW_NOTES_T1-T6.md` |
| Phase 2 | not stated | Phase-2 test, CI, and documentation checklist was fully completed. | `lr-campaign-archive:REFACTORING_PHASE2.md` |
| Math audit T11 | not stated | Confirmed uncalled math routines were deleted and the full source-tree audit passed. | `lr-campaign-archive:docs/dev/MATH_AUDIT_T11.md` |

### XC

| Campaign / file | Date | Conclusion (≤ 25 words) | Path at tag |
|---|---|---|---|
| XC general closeout | not stated | Supported legacy and libXC XC scope was accepted; unsupported requests fail explicitly. | `lr-campaign-archive:docs/XC_GENERAL_CLOSEOUT.md` |
| XC GGA radial reconciliation | not stated | Legacy radial GGA and libXC derivative paths were reconciled. | `lr-campaign-archive:docs/XC_GGA_RADIAL_RECONCILIATION.md` |
| XC integration corrective closeout | not stated | Radial derivative handoff and magnetic-SCF residual bookkeeping defects were closed. | `lr-campaign-archive:docs/XC_INTEGRATION_CORRECTIVE_CLOSEOUT.md` |
| XC LDA reconciliation | 2026-08-30 | Barth-Hedin derivative and zero-spin-channel handling defects were repaired. | `lr-campaign-archive:docs/XC_LDA_RECONCILIATION.md` |
| XC magnetic-SCF closeout | not stated | Bounded magnetic-SCF evidence remained unresolved for promoting a converged high-spin fcc-Fe state. | `lr-campaign-archive:docs/XC_MAGNETIC_SCF_CLOSEOUT.md` |
