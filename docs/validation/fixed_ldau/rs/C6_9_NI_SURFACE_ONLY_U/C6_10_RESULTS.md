# C6.10 matched Ni(001) response and gates

## Matched two-step response

The U=0 control was staged from the same original no-U potential files as the U run. Both decks use the same geometry, labels, nsp=1, `nstep=2`, `conv_thr=5e-8`, fixed E_F=-0.062420 Ry, energy window [-1.0,1.5] Ry, and mixing. Both explicitly select `hubbard_radial_moments='local'`; the only physical deck difference is `Ni-S` type 3 d-channel U/J=4/1 versus 0/0. Both used the same executable SHA256 `2cb9e51f49d99db031aebf7d71300b8494e7a2f47e1d16819d6c0d35e3bd331` and passed `run_rslmto_checked.py`, including the host-E_F check at zero delta. Each completed two SCF steps.

Site indices in surface output map as atom 1 = ES recursion site, atom 2 = surface Ni (`Ni-S`), atom 3 = subsurface Ni (`Ni-S-1`). The frozen bulk Ni host is type 1, not a reported recursion atom.

| U/J on surface d | Step | Diff | ES moment (μB) | Surface Ni moment (μB) | Subsurface Ni moment (μB) |
|---|---:|---:|---:|---:|---:|
| 0/0 | 1 | 8.53e-8 | 0.007355 | 0.689004 | 0.640932 |
| 0/0 | 2 | 7.85e-8 | 0.007355 | 0.689004 | 0.640932 |
| 4/1 | 1 | 0.1149764639 | 0.037554 | 1.081094 | 0.489236 |
| 4/1 | 2 | 0.0389620957 | 0.010834 | 0.112752 | 0.129423 |

The control moments are stable to the displayed precision. The finite-U moments continue to change strongly between steps, so these are response values, not converged moments.

The output `report.out` contains final per-spin s/p/d sphere occupations. For the U=0 two-step control, surface d occupations are 4.647171 up + 3.933155 down = 8.580326 total; subsurface d occupations are 4.631449 + 3.949601 = 8.581050. For the finite-U run after step 2, they are 4.637807 + 3.924348 = 8.562155 at the surface and 4.631213 + 3.952222 = 8.583435 at the subsurface. ES d occupations are 0.019718 + 0.020090 = 0.039808 at U=0 and 0.019792 + 0.020127 = 0.039919 at finite U.

For stepwise surface d information, `HUBBARD_QLCMP` reports up/down projected d occupations 4.647171/3.933155 at finite-U step 1 and 4.637156/3.921356 at step 2. At step 2, the local Hubbard occupation matrix separately reports n(up)=2.644274 and n(down)=1.573365 (total 4.217638); this is a different projector measure from `report.out` sphere occupations. The U=0 output shows unchanged sphere occupations to the printed precision between the original one-step no-U result and the final two-step result; the solver does not emit the same per-step `HUBBARD_QLCMP` trace when U=0.

The final stored potentials confirm scope: `Ni_out.nml` (frozen host), `ES_out.nml`, and `Ni-S-1_out.nml` have U=J=0. `Ni-S_out.nml` alone has d-channel U=0.2939943429 Ry and J=0.0734985857 Ry (the stored conversion of deck values 4/1). Thus no U was assigned to the host, empty sphere, or subsurface.

## Current-binary gates

The unit-test build relinked the current executable at SHA256 `b0621f0b30c26f1bfd3cde5df9b2dfa9cf93e99dd299c260a118acd381917914`. The current-binary reruns all passed:

- RS identity: 27/27 against the versioned refreshed RS oracle.
- KS identity: 24/24 against the preserved C6.6 current-runtime KS candidate.
- RA2 provenance/guard fixtures: 13/13; records are in `C6_10_GUARDS.json` and raw cases under `/private/tmp/c610-final-guards/`.
- Surface frozen-plus-recursive mapping/local-moment oracle: `UnitLdaUProjectorMultisiteFrozenH` passed.
- Impurity and surface auto-mode fixtures: `UnitAutoImpurityLocalMode` passed (one test covers both fixtures).

Ctest output is preserved at `build/Testing/Temporary/LastTest.log`. The current-binary identity result is `C6_10_IDENTITY_FINAL.json`.
