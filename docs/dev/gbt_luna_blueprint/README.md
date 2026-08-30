# GBT Luna Blueprint Package

This directory contains the master blueprint and sliced Luna prompts for the `fable_v3` generalized Bloch theorem revalidation/repair campaign in RS-LMTO-ASA.

## Execution order

1. `00_HEAD_GBT_REVALIDATION_MASTER_BLUEPRINT.md`
2. `LUNA_GBT_WP00_CURRENT_STATE_AND_EVIDENCE_MAP.md`
3. `LUNA_GBT_WP01_Q0_EXACT_REDUCTION.md`
4. `LUNA_GBT_WP02_GAUGE_HERMITICITY_SHIFTED_K.md`
5. `LUNA_GBT_WP03_COMPOSITE_LMTO_COVARIANCE.md`
6. `LUNA_GBT_WP04_MATCHED_COMMENSURATE_SUPERCELL.md`
7. `LUNA_GBT_WP05_ROTATING_FRAME_CONTRACT.md`
8. `LUNA_GBT_WP06_CONSTRAINT_FIELD_COVARIANCE.md`
9. `LUNA_GBT_WP07_CONSTRAINT_ENERGY_SEMANTICS.md`
10. `LUNA_GBT_WP08_CORRECTED_CONSTRAINED_MFT.md`
11. `LUNA_GBT_WP09_HARMONIC_CONE_AND_KGRID.md`
12. `LUNA_GBT_WP10_SMALL_Q_QUADRATIC_LIMIT.md`
13. `LUNA_GBT_WP11_LKAG_JQ_CONSISTENCY.md`
14. `LUNA_GBT_WP12_SCF_DIAGNOSIS_AND_REPAIR.md`

The hard campaign rule is that SCF repair comes last. The primitive `S x G` GBT operator is held fixed unless an earlier exact equivalence test proves it wrong.
