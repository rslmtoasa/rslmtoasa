# RS-LMTO-ASA TD-DFT clean-room redevelopment package

Target branch: `fable_v4`
Intended coding model: **Luna**
Blueprint date: 2026-09-09

## Prime rule

A missing published derivation, ambiguous basis mapping, unresolved sign/factor, or unavailable ground-state quantity is a **BLOCKER**, not an invitation to invent a workaround.

The campaign is clean-room: the current TD-DFT physics layer is removed from the active source tree before replacement physics is written. Git history remains the archive.

## Naming and layout

Do not use author/institution names as implementation identifiers. Use neutral physics names (`lr_*`, `tddft_*`, `goldstone_*`).

Keep new Fortran files directly in `source/`. Do **not** introduce new source subfolders.

## Execution order

### Clean-room/blocker phase
1. `01_LR-00_PURGE_CLEANROOM.md`
2. `02_LR-01_RADIAL_GROUND_STATE_AUDIT.md`
3. `03_LR-02_RESPONSE_BASIS_FEASIBILITY.md`
4. `04_LR-03_CONVENTIONS_AND_NORMALIZATION.md`

If any mandatory gate is BLOCKED, stop.

### Conditional implementation phase
5. `05_LR-04_RADIAL_GROUND_STATE_INTERFACE.md`
6. `06_LR-05_RESPONSE_BASIS.md`
7. `07_LR-06_KS_SUSCEPTIBILITY.md`
8. `08_KXC-01_ALSDA_TRANSVERSE_KERNEL.md`
9. `09_GSR-01_GOLDSTONE_SUMRULE.md`
10. `10_GCR-01_GOLDSTONE_EIGENVALUE_CORRECTION.md`
11. `11_TDDY-01_DYSON_AND_LOSS.md`
12. `12_TDVAL-01_COLLINEAR_VALIDATION.md`

### Deferred audits
13. `13_NC-00_NONCOLLINEAR_2026_FEASIBILITY.md`
14. `14_SV-00_LMTO_STERNHEIMER_FEASIBILITY.md`

## Primary published foundations

- P. Buczek, A. Ernst, L. M. Sandratskii, Phys. Rev. B 84, 174418 (2011).
- S. Lounis, A. T. Costa, R. B. Muniz, D. L. Mills, Phys. Rev. B 83, 035109 (2011).
- D. Eilmsteiner, A. Ernst, P. A. Buczek, arXiv:2603.03220v2 (2026).
- S. Y. Savrasov, Phys. Rev. Lett. 81, 2570 (1998).

## Production labels

Every physics statement must be classified as:
- `[LITERATURE]`
- `[BASIS MAPPING]`
- `[IMPLEMENTATION]`
- `[HYPOTHESIS]`

`[HYPOTHESIS]` is forbidden from the production physics path.
