# User Ni(001) U=4/J=1 host and surface probe

Date: 2026-10-25. The untouched user `finiteU` directory was copied to
`source_snapshot/`; all work below is in this test folder. The original path was
`/Users/andersb/Jobb/RS-utils/tmp/Cheb/fccNi_small/finiteU`.

## U=4/J=1 bulk host regeneration

`host_u4j1/input.nml` starts from the accepted RS-local Ni U=3/J=0.7 host deck;
its only changed deck lines are U=4 and J=1. It uses the accepted host
potential as restart, explicit local mode, and fixed E_F = -0.061141 Ry. The
checked wrapper preflight passed with host-E_F delta 0.0 Ry. The solver reported
`mode=local source=explicit`, `local_energy_adopted=.true.`, and wrote an
RS-local Ni potential. All 100 configured steps completed; final Diff was
`7.0e-9`, while the log continued to say `Not converged` because inherited
`conv_thr=1e-30`. The potential is therefore a diagnostic local-provenance
seed, not an accepted C5.2 reference. E_F was fixed to the U=3 seed endpoint,
not relaxed to a new U=4 endpoint.

The accepted Ni bulk seed is spin symmetric. This U=4 host run retained that
state, so it does not preserve the spin-polarized state present in the user's
original surface output. That matters for interpreting the small surface
moments below; this regenerated host is not a faithful magnetic replacement
for the user's `Ni_out.nml` seed.

## Surface deck and runs

The generated surface geometry contains projected layers at `-0.5`, `0.0`,
and `0.5` (then repeats every `0.5`) in `surfclu.out`; these are in `alat`
units, not absolute Cartesian coordinates. The source `buildsurf` input has no
absolute layer positions. `DECK_AUDIT.md` lists the staged deck changes.

All three surface attempts passed the wrapper's host-E_F check (delta 0.0 Ry)
and selected `mode=local` under `auto`. There was no provenance error.
Surface energy routing remained off (`local_energy_adopted=.false.`), as
intended. The first, host-window attempt (`[-1.0, 0.5]` Ry) failed in its first
recursion with NaN projected DOS and
`Invalid local d-shell Hubbard moment correction; status=5 type=3 atom=2 spin=1`.

The widened-window attempt (`[-3.0, 1.5]` Ry) completed one step without NaNs.
The step Diff was `0.0232230893`, above the requested `5e-8` target. Its finite
first-step d occupations were:

| Deck type / layer label | n↑ | n↓ | n↑−n↓ |
|---|---:|---:|---:|
| 1 / Ni host | 4.747080 | 4.747080 | 0.000000 |
| 3 / Ni-S surface | 4.773911 | 4.770410 | 0.003501 |
| 4 / Ni-S-1 subsurface | 4.765654 | 4.765520 | 0.000134 |

The corresponding atom spin moments were `0.030722`, `0.006568`, and
`0.000677` μB for atoms 1–3 (Ni host, surface Ni, and subsurface Ni in the
`surfclu.out` map). Step 2 then generated NaN projected DOS / `g0` values and
hit the same status-5 fatal at type 3, atom 2, spin 1. It produced no valid
second-step moments and did not meet the two-step Diff target.

A final diagnostic removed type-4 U/J assignments from `input.nml`, but it
still used a type-4 seed copied from the U=4 bulk host. That seed's output
namelist embeds U=4/J=1 in `hubbard_u`/`hubbard_j`; runtime therefore still
applied U to type 4. It reproduced the first-step values and the second-step
failure, so it is **not** a valid type-4-off isolation. The untouched user
`source_snapshot/Ni-S-1_out.nml` has zero U/J; a proper type-4-off control must
retain that seed (and use the user's magnetic Ni host seed if the goal is to
preserve the original magnetic state).

## Disposition

**Nonzero-U surface response remains OPEN.** The two-step probe failed; the
small finite first-step moments are not a reliable physical result because the
regenerated bulk seed was spin symmetric. The exact U=0/no-U equivalence gate
remains separately PASS as recorded in `C6_7_ZERO_U_COMPARISON/RESULTS.md`.
Four SCF runs were used in this probe folder (one bulk host, three surface
attempts); no further run was launched.
