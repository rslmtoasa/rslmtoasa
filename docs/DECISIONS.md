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

## Linear response (LR) campaigns, 2026

### Condensed into docs/linear_response/

| Campaign / file | Date | Conclusion (≤ 25 words) | Path at tag |
|---|---|---|---|
| TDDFT_FORMULATION.md | 2026-09-28 | Consolidated formulation scope, production status, route boundaries, and the live `&linear_response` compatibility contract. | `lr-campaign-archive:docs/TDDFT_FORMULATION.md` |
| NATIVE_ROTATION_DYNAMICS.md | 2026-09-28 | Consolidated the production rotation kernel, exact Gamma identity, static Hessian reduction, and pole-control policy. | `lr-campaign-archive:docs/NATIVE_ROTATION_DYNAMICS.md` |
| LR_TDDFT_CONVENTIONS.md | 2026-09-28 | Consolidated certified units, signs, Fourier phase, circular channels, radial measures, and Pauli-response conventions. | `lr-campaign-archive:docs/LR_TDDFT_CONVENTIONS.md` |
| KXC_ALSDA_TRANSVERSE_KERNEL.md | 2026-09-28 | Consolidated the local ALSDA interaction and raw static diagnostics without promoting Goldstone or material claims. | `lr-campaign-archive:docs/KXC_ALSDA_TRANSVERSE_KERNEL.md` |
| GOLDSTONE_SUMRULE_INTERACTION.md | 2026-09-28 | Consolidated the independent LCMM sum-rule interaction and its finite-response diagnostic boundary. | `lr-campaign-archive:docs/GOLDSTONE_SUMRULE_INTERACTION.md` |
| TDDFT_DYSON_AND_LOSS.md | 2026-09-28 | Consolidated canonical Dyson algebra, metric placement, retarded loss, conditioning, and route provenance. | `lr-campaign-archive:docs/TDDFT_DYSON_AND_LOSS.md` |
| TDDFT_PRODUCTION_DRIVER.md | 2026-09-28 | Consolidated post-SCF orchestration and the boundary between registered reciprocal routes and provider-baseline evidence. | `lr-campaign-archive:docs/TDDFT_PRODUCTION_DRIVER.md` |
| TDDFT_RS_GF_BACKEND.md | 2026-09-28 | Consolidated the directed GF-provider contract, real-axis bubble, q phase, controls, and convergence reporting. | `lr-campaign-archive:docs/TDDFT_RS_GF_BACKEND.md` |
| RSGF_ENDPOINT_AUGMENTATION.md | 2026-09-28 | Consolidated four-branch endpoint augmentation, contact terms, effective-Hamiltonian provenance, and capability scope. | `lr-campaign-archive:docs/RSGF_ENDPOINT_AUGMENTATION.md` |
| LR_KS_SUSCEPTIBILITY.md | 2026-09-28 | Consolidated the Pauli transverse Kohn-Sham susceptibility and its reciprocal endpoint ordering. | `lr-campaign-archive:docs/LR_KS_SUSCEPTIBILITY.md` |
| LR_LMTO_PRODUCT_RESPONSE_BASIS.md | 2026-09-28 | Consolidated the weighted compact product basis, endpoint branches, rank diagnostics, and response-space mapping. | `lr-campaign-archive:docs/LR_LMTO_PRODUCT_RESPONSE_BASIS.md` |
| LR_PAULI_TRANSITION_VERTEX.md | 2026-09-28 | Consolidated the Pauli transition vertex, radial measure, angular normalization, and deferred exact-SR boundary. | `lr-campaign-archive:docs/LR_PAULI_TRANSITION_VERTEX.md` |
| LR_RESPONSE_BASIS_MAPPING.md | 2026-09-28 | Consolidated direct-coordinate q/k mapping, endpoint gauges, radial/angular indices, and canonical response ordering. | `lr-campaign-archive:docs/LR_RESPONSE_BASIS_MAPPING.md` |
| LR_RESPONSE_SPACE_ALGEBRA.md | 2026-09-28 | Consolidated raw/canonical metric algebra, compact coordinates, and operator composition rules. | `lr-campaign-archive:docs/LR_RESPONSE_SPACE_ALGEBRA.md` |
| DRESP_01_PROJECTED_SITE_SPIN_CONTRACT.md | 2026-09-28 | Consolidated projected site spin operators, moments, and the bounded site-space contract. | `lr-campaign-archive:docs/DRESP_01_PROJECTED_SITE_SPIN_CONTRACT.md` |
| DRESP_02_PROJECTED_RECIPROCAL_CHI0.md | 2026-09-28 | Consolidated projected reciprocal site `chi0`, circular normalization, q covariance, and broadening diagnostics. | `lr-campaign-archive:docs/DRESP_02_PROJECTED_RECIPROCAL_CHI0.md` |
| DRESP_03G_FINITE_H_CONTOUR_GF.md | 2026-09-28 | Consolidated finite-H contour GF oracles and their role as validation infrastructure. | `lr-campaign-archive:docs/DRESP_03G_FINITE_H_CONTOUR_GF.md` |
| DRESP_03Q_FINITE_Q_LKAG_CLOSURE.md | 2026-09-28 | Consolidated finite-q ordered-pair LKAG closure and q-reversal checks. | `lr-campaign-archive:docs/DRESP_03Q_FINITE_Q_LKAG_CLOSURE.md` |
| DRESP_03T_LOCAL_TORQUE_HESSIAN_REPAIR.md | 2026-09-28 | Consolidated local torque/Hessian repair diagnostics and their non-production claim boundary. | `lr-campaign-archive:docs/DRESP_03T_LOCAL_TORQUE_HESSIAN_REPAIR.md` |
| DRESP_03TG_NATIVE_TUREK_GF.md | 2026-09-28 | Consolidated native Turek GF exchange construction, ordered pairs, and finite-H diagnostics. | `lr-campaign-archive:docs/DRESP_03TG_NATIVE_TUREK_GF.md` |
| DRESP_03TG_NATIVE_TUREK_CONTOUR_CLOSE.md | 2026-09-28 | Consolidated native Turek contour closure and bounded convergence evidence without a stiffness claim. | `lr-campaign-archive:docs/DRESP_03TG_NATIVE_TUREK_CONTOUR_CLOSE.md` |
| DRESP_04_PROJECTED_MILLS_RPA.md | 2026-09-28 | Consolidated projected Mills/Stoner interaction, site Dyson solve, loss, and scalarization boundary. | `lr-campaign-archive:docs/DRESP_04_PROJECTED_MILLS_RPA.md` |
| DRESP_05_PROJECTED_JUELICH_LCMM.md | 2026-09-28 | Consolidated projected Jülich/LCMM interaction, eta stability, and same-state route ledger. | `lr-campaign-archive:docs/DRESP_05_PROJECTED_JUELICH_LCMM.md` |
| TG_FZ_R5_PRODUCTION_H_REPRESENTATION_BRIDGE.md | 2026-09-28 | Consolidated the blocked exact-H to production-H bridge and preserved H2/H_exact as distinct contracts. | `lr-campaign-archive:docs/TG_FZ_R5_PRODUCTION_H_REPRESENTATION_BRIDGE.md` |
| TG_FZ_R7_SCREENING_REPRESENTATION_REPAIR.md | 2026-09-28 | Consolidated repaired screening representation and exact-H/native-gamma diagnostics. | `lr-campaign-archive:docs/TG_FZ_R7_SCREENING_REPRESENTATION_REPAIR.md` |
| TG_FZ_R8R_RESOLVENT_CONTACT_CLOSURE.md | 2026-09-28 | Consolidated fixed-z resolvent contact closure and retained H2 truncation diagnostics. | `lr-campaign-archive:docs/TG_FZ_R8R_RESOLVENT_CONTACT_CLOSURE.md` |

### DRESP

| Campaign / file | Date | Conclusion (≤ 25 words) | Path at tag |
|---|---|---|---|
| DRESP_00_ARCHITECTURAL_REBASE_AUDIT.md | 2026-09-16 | Architectural rebase complete; implementation blocked at the RUNG_1_KL projected site-spin contract. | `lr-campaign-archive:docs/DRESP_00_ARCHITECTURAL_REBASE_AUDIT.md` |
| DRESP_02C_MATERIAL_CLOSURE.md | 2026 | PASS — DRESP-02 material closure complete. | `lr-campaign-archive:docs/DRESP_02C_MATERIAL_CLOSURE.md` |
| DRESP_02F_GF_LEHMANN_FORMULATION_AUDIT.md | 2026 | PASS — current GF formulation proven equivalent. | `lr-campaign-archive:docs/DRESP_02F_GF_LEHMANN_FORMULATION_AUDIT.md` |
| DRESP_02P_PROJECTED_GF_PERFORMANCE.md | 2026 | PASS — performance blocker removed. | `lr-campaign-archive:docs/DRESP_02P_PROJECTED_GF_PERFORMANCE.md` |
| DRESP_02R_MATERIAL_GF_CLOSURE.md | 2026 | Historical remediation record superseded by PASS — DRESP-02 material closure complete. | `lr-campaign-archive:docs/DRESP_02R_MATERIAL_GF_CLOSURE.md` |
| DRESP_03QM_METALLIC_FINITE_Q.md | 2026 | PASS — metallic formulation certified; performance open. | `lr-campaign-archive:docs/DRESP_03QM_METALLIC_FINITE_Q.md` |
| DRESP_03R_LMTO_EXCHANGE_VERTEX_MAPPING.md | 2026 | Partial PASS for offsite canonical/auxiliary contraction; physical path-operator to orthogonal-resolvent map remains BLOCKED. | `lr-campaign-archive:docs/DRESP_03R_LMTO_EXCHANGE_VERTEX_MAPPING.md` |
| DRESP_03_KL_STATIC_LKAG_BRIDGE.md | 2026 | Audit complete; overall result BLOCKED at the finite-Hamiltonian/native LKAG comparison gate. | `lr-campaign-archive:docs/DRESP_03_KL_STATIC_LKAG_BRIDGE.md` |
| DRESP_06A_SPATIAL_ALSDA_COMPARISON.md | 2026 | Historical ledger retained; current DRESP-06A-R uses certified accepted Pauli magnetization in direct checks. | `lr-campaign-archive:docs/DRESP_06A_SPATIAL_ALSDA_COMPARISON.md` |
| DRESP_07_EXACT_KS_WARD_CLOSURE.md | 2026 | PASS — localized second-order LMTO mapping required. | `lr-campaign-archive:docs/DRESP_07_EXACT_KS_WARD_CLOSURE.md` |
| DRESP_08_NATIVE_SECOND_ORDER_FIELD.md | 2026 | Verdict: PASS for the bounded native second-order LMTO field mapping. | `lr-campaign-archive:docs/DRESP_08_NATIVE_SECOND_ORDER_FIELD.md` |
| DRESP_09R_PAULI_NATIVE_GATE.md | 2026 | PAULI_PROJECTION_INSUFFICIENT_FOR_NATIVE_TANGENT; mixed-spin radial metric remains not certified. | `lr-campaign-archive:docs/DRESP_09R_PAULI_NATIVE_GATE.md` |
| DRESP_09S_SCALAR_RELATIVISTIC_AUGMENTATION.md | 2026 | PASS-B. | `lr-campaign-archive:docs/DRESP_09S_SCALAR_RELATIVISTIC_AUGMENTATION.md` |
| DRESP_09T_MOVING_BASIS_TANGENT.md | 2026 | BLOCKED at the live orthogonalization response. | `lr-campaign-archive:docs/DRESP_09T_MOVING_BASIS_TANGENT.md` |
| DRESP_09U_ORTHOGONALIZATION_RESPONSE.md | 2026 | PASS-B for the Hamiltonian representation chain; density coefficient response remains open. | `lr-campaign-archive:docs/DRESP_09U_ORTHOGONALIZATION_RESPONSE.md` |
| DRESP_09V_DENSITY_MOMENT_TANGENT.md | 2026 | Closes the density-side representation response using the live reciprocal LMTO energy-moment contract. | `lr-campaign-archive:docs/DRESP_09V_DENSITY_MOMENT_TANGENT.md` |
| DRESP_09W_RADIAL_OBSERVABLE_PROVENANCE.md | 2026 | Audits the radial seam after coefficient-space M0/M1/M2 and tangent closures. | `lr-campaign-archive:docs/DRESP_09W_RADIAL_OBSERVABLE_PROVENANCE.md` |
| DRESP_09X_SR_SPIN_OBSERVABLE.md | 2026 | Bounded scalar-relativistic physical Pauli-spin observable audit; no ALSDA/Ward or dynamics claim. | `lr-campaign-archive:docs/DRESP_09X_SR_SPIN_OBSERVABLE.md` |
| DRESP_09Y_AUGMENTATION_FRAME_TANGENT.md | 2026 | Closes the missing observable-side term in the bounded rigid-rotation audit. | `lr-campaign-archive:docs/DRESP_09Y_AUGMENTATION_FRAME_TANGENT.md` |
| DRESP_09ZSR_PRODUCTION_SPAN_RECONCILIATION.md | 2026 | PASS-A; live production four-branch space is a strict numerical subspace of the complete six-branch space. | `lr-campaign-archive:docs/DRESP_09ZSR_PRODUCTION_SPAN_RECONCILIATION.md` |
| DRESP_09ZS_COMPACT_SPAN_AUDIT.md | 2026 | PASS-B; raw four-branch shadow is not closed for the complete six-branch radial space. | `lr-campaign-archive:docs/DRESP_09ZS_COMPACT_SPAN_AUDIT.md` |
| DRESP_09Z_ARBITRARY_L_SR_VERTEX.md | 2026 | Records the implementation boundary reached for the arbitrary-L scalar-relativistic response vertex. | `lr-campaign-archive:docs/DRESP_09Z_ARBITRARY_L_SR_VERTEX.md` |
| DRESP_09_TRANSVERSE_FIELD_INSERTION.md | 2026 | Compact field algebra uses weighted coordinates with hard pairing b^H d = B^H W D. | `lr-campaign-archive:docs/DRESP_09_TRANSVERSE_FIELD_INSERTION.md` |
| DRESP_10A_SIX_BRANCH_PRODUCT_BASIS.md | 2026 | Promotes the compact LMTO product basis to six certified second-order endpoint branches. | `lr-campaign-archive:docs/DRESP_10A_SIX_BRANCH_PRODUCT_BASIS.md` |
| DRESP_10R_RESPONSE_ARCHITECTURE.md | 2026 | Raw static full-spatial ALSDA Ward gate stops when upstream representation gates remain open. | `lr-campaign-archive:docs/DRESP_10R_RESPONSE_ARCHITECTURE.md` |
| DRESP_10_FIXED_BASIS_MIXED_WARD.md | 2026 | Fixed-ground-state mixed-representation static Ward formulation excludes moving-basis and contact terms. | `lr-campaign-archive:docs/DRESP_10_FIXED_BASIS_MIXED_WARD.md` |
| DRESP_10_RAW_FULL_SPATIAL_ALSDA_WARD.md | 2026 | Diagnostic material gate adds full-spatial ALSDA orchestration; dynamic spectrum, Dyson, correction, and mode fit remain out of scope. | `lr-campaign-archive:docs/DRESP_10_RAW_FULL_SPATIAL_ALSDA_WARD.md` |
| DRESP_11_GOLDSTONE_DEFECT_SPECTRUM.md | 2026 | Diagnostic-only; no Goldstone correction is enabled in production dynamics. | `lr-campaign-archive:docs/DRESP_11_GOLDSTONE_DEFECT_SPECTRUM.md` |
| DRESP_12_FINITE_LMTO_COVARIANCE.md | 2026 | Static diagnostic-only Gamma-point covariance closure; response kernel, dynamics, and Goldstone correction remain unchanged. | `lr-campaign-archive:docs/DRESP_12_FINITE_LMTO_COVARIANCE.md` |

### TG_FZ (Turek fixed-z)

| Campaign / file | Date | Conclusion (≤ 25 words) | Path at tag |
|---|---|---|---|
| TG_FZ_R1_FIXED_Z_DIAGNOSTIC_REPAIR.md | 2026 | PASS — diagnostic fixture repaired and hardened. | `lr-campaign-archive:docs/TG_FZ_R1_FIXED_Z_DIAGNOSTIC_REPAIR.md` |
| TG_FZ_R2_SPIN_SCREENING_VERTEX.md | 2026 | PASS-A — exact fixed-z vertex covariance derived and verified. | `lr-campaign-archive:docs/TG_FZ_R2_SPIN_SCREENING_VERTEX.md` |
| TG_FZ_R4_FINITE_H_TUREK_VERTEX_BRIDGE.md | 2026 | BLOCKED — no unfitted fixed-z identity between live finite-H torque and tested Turek vertex. | `lr-campaign-archive:docs/TG_FZ_R4_FINITE_H_TUREK_VERTEX_BRIDGE.md` |
| TG_FZ_R6_SCREENING_NORMALIZATION_ALPHA_AUDIT.md | 2026 | PASS-B — normalized gamma and alpha authority are both required. | `lr-campaign-archive:docs/TG_FZ_R6_SCREENING_NORMALIZATION_ALPHA_AUDIT.md` |

### LR-* / RSGF / TDDFT audits

| Campaign / file | Date | Conclusion (≤ 25 words) | Path at tag |
|---|---|---|---|
| LR-BASIS-00_RADIAL_AUGMENTATION_EVIDENCE.md | 2026 | Baseline certified for the stated scope. | `lr-campaign-archive:docs/LR-BASIS-00_RADIAL_AUGMENTATION_EVIDENCE.md` |
| LR-GF-01_GF_CONTRACT_EVIDENCE.md | 2026-09-10 | Certifies the live post-purge one-electron infrastructure; removed TD-DFT response is not restored. | `lr-campaign-archive:docs/LR-GF-01_GF_CONTRACT_EVIDENCE.md` |
| LR_GF_SUSCEPTIBILITY_CROSSCHECK.md | 2026 | PASS — independent reciprocal-GF bubble agrees with the LR-06 spectral/Lehmann susceptibility on the certified fixture. | `lr-campaign-archive:docs/LR_GF_SUSCEPTIBILITY_CROSSCHECK.md` |
| LR_RADIAL_GROUND_STATE_AUDIT.md | 2026 | PASS for the supported scalar-relativistic, collinear, two-channel radial ground-state contract. | `lr-campaign-archive:docs/LR_RADIAL_GROUND_STATE_AUDIT.md` |
| LR_RS_GF_REPRESENTATION_AUDIT.md | 2026-09-11 | Coefficient-space Green-function representation certified; native RS response backend remains BLOCKED. | `lr-campaign-archive:docs/LR_RS_GF_REPRESENTATION_AUDIT.md` |
| LR_SR_PAULI_NUMERICAL_CLOSURE.md | 2026 | Records scalar-relativistic radial density versus large-component Pauli projection for one bcc-Fe state. | `lr-campaign-archive:docs/LR_SR_PAULI_NUMERICAL_CLOSURE.md` |
| RSGF_CAPABILITY_CLOSURE.md | 2026-09-11 | Prerequisite capabilities are not globally blocked; R0–R3 are closed at certified scope while Fe/Ni R4 remains pending. | `lr-campaign-archive:docs/RSGF_CAPABILITY_CLOSURE.md` |
| TDDFT_CLEANROOM_PURGE_MANIFEST.md | 2026 | LR-00 inventory recorded before purge. | `lr-campaign-archive:docs/TDDFT_CLEANROOM_PURGE_MANIFEST.md` |
| TDDFT_COLLINEAR_REVALIDATION.md | 2026-09-11 | PENDING for R4 material validation; R0–R3 production revalidation passed. | `lr-campaign-archive:docs/TDDFT_COLLINEAR_REVALIDATION.md` |
| TDDFT_NATIVE_RSGF_PRODUCTION_INTEGRATION.md | 2026-09-11 | PASS for registered production-integration scope; R4 Fe/Ni material validation remains pending. | `lr-campaign-archive:docs/TDDFT_NATIVE_RSGF_PRODUCTION_INTEGRATION.md` |
| TDDFT_RECIPROCAL_GF_QUADRATURE_AUDIT.md | 2026 | Historical TDVK-03A fixture audit; TDVK-03R factorized backend and Fe closure evidence recorded below. | `lr-campaign-archive:docs/TDDFT_RECIPROCAL_GF_QUADRATURE_AUDIT.md` |
| TDDFT_RECIPROCAL_VALIDATION.md | 2026-09-12 | TDVAL-K preflight recorded reciprocal validation scope; no explicit verdict appears in the first 40 lines. | `lr-campaign-archive:docs/TDDFT_RECIPROCAL_VALIDATION.md` |

### Goldstone correction

| Campaign / file | Date | Conclusion (≤ 25 words) | Path at tag |
|---|---|---|---|
| GOLDSTONE_EIGENVALUE_CORRECTION.md | 2026 | Focused tests establish algebraic consistency and finite-matrix evidence, not converged material or literature-spectrum accuracy. | `lr-campaign-archive:docs/GOLDSTONE_EIGENVALUE_CORRECTION.md` |

- 2026-10-03 — LR-METHOD-05R: accepted SR ALSDA kernel repaired; BLOCKED — response-representation closure not established (Fe raw Ward relative residual ≈0.243; no correction).
