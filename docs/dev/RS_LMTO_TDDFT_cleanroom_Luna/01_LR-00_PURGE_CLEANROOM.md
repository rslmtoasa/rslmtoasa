# LR-00 — Purge current TD-DFT physics and establish a clean-room baseline

## Goal
Remove the active TD-DFT physics implementation so redevelopment cannot inherit old assumptions. Write no replacement TD-DFT physics.

## Precondition
Working tree must be clean.

Create archival tag:
`pre_tddft_cleanroom_20260909`

If it already exists at another commit, stop.

## Inventory first
Create `docs/TDDFT_CLEANROOM_PURGE_MANIFEST.md`.

Classify every related file/routine:
- `DELETE_PHYSICS`
- `RETAIN_GENERIC`
- `RETAIN_GROUND_STATE`
- `REVIEW_DEPENDENCY`

Cover at least:
- all current `tddft_*` modules;
- pair-potential/magnetic-tangent code used only for old TD-DFT;
- TD-DFT response vertices/conventions;
- TD-DFT driver branches and namelist handling;
- TD-DFT tests;
- documentation claiming validation.

## Purge
Remove active code encoding:
- pair-potential Xi;
- Hamiltonian-rotation tangent as TDDFT interaction;
- old scalar/site TDDFT kernels;
- old TDDFT Dyson physics;
- old Goldstone correction/diagnostic semantics;
- mode machinery tied to those objects;
- tests whose expected physics derives from them.

Do not leave commented or `legacy_*` production copies. Git is the archive.

## Retain only independent infrastructure
Generic eigensolvers, Green functions, k worksets, radial SCF/XC data, matrix algebra and Fourier utilities may remain if documented as physics-neutral.

## User-facing behavior
Old TD-DFT input must fail explicitly:
`TD-DFT temporarily unavailable during literature-locked clean-room redevelopment.`

No silent fallback.

## Gate
- full code compiles;
- non-TDDFT tests pass;
- no new TDDFT physics exists.

## Checklist
- [x] clean tree recorded
- [x] archive tag created
- [x] purge manifest written
- [x] old TDDFT physics removed
- [x] pair-Xi unavailable
- [x] no commented legacy shortcut
- [x] retained generic code justified
- [x] old TDDFT input fails explicitly
- [x] non-TDDFT build/tests pass
- [x] no replacement physics introduced

## Commit
`td-dft: establish clean-room linear-response baseline`
