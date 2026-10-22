# C6.1 — Stage RE smoke validation disposition

## T1: no-SCF gates

`UnitAutoImpurityLocalMode` was rebuilt and passed (1/1), covering the impurity auto fixture and the new surface auto unit fixture with a local U-carrying frozen host. The single-thread no-U RS/KS identity and 13 RA2 guard scenarios have recorded PASS evidence from C5.13 on the parent tree; they were not rerun here because those checks require solver launches, contrary to this T1's no-SCF restriction. No T1 gate failure was observed.

## T2: Ni(001) preparation and stop

The `fccCu001` fixture builds surface positions from fractional crystal coordinates and scales them by `alat`; no absolute layer-z coordinates are specified. Geometry therefore scales with the permitted Ni `alat` change, and no geometry edit was made.

The worker staged T2(a) and T2(b) copies in `/private/tmp/c61_candidate_decks_fa68c56c`. The audited input diff changes only the permitted fields:

- `alat`: 3.614 -> 3.5234 and `wav`: 1.41237 -> 1.37692, from the accepted Ni bulk deck.
- `label(1)`: Cu -> Ni and `label(3)`: Cu-S -> Ni-S; label 2 remains ES. `Ni.nml` and `Ni-S.nml` are copies of the accepted Ni bulk potential with the symbol changed to match each label; `ES.nml` is unchanged.
- `fermi`: -0.088376 -> -0.061141 and `energy_max`: 1.2 -> 0.5; `energy_min` remains -1.0, matching the accepted Ni endpoint window.
- The appended Hubbard group sets auto mode and U/J only for the Ni d channels: 0/0 on both Ni types in T2(a), and 3.0/0.7 on both in T2(b). All other fixture lines, including `nsp=2`, the surface geometry, ES, and recursion setup, are preserved.

The source's surface auto-local eligibility requires `control%nsp == 1` (`source/self.f90::is_eligible_local_moments`). The copied surface fixture sets `nsp=2`, and C6.1's deck rule does not allow changing it. Therefore T2(b) cannot resolve to local mode as required. The required gate is structurally unreachable with this fixture under the prescribed edits. To honor the gate, no T2 SCF run or host-E_F preflight was launched. T2 occupations, moments, provenance outcome, and convergence remain unmeasured.

## T3/T4 disposition

Stopped before T3 because the T2(b) gate cannot be satisfied. No Cu host generation, T3 surface run, or further preflight was performed. The four-run budget is unused (0/4). No merge was made and nothing was pushed. Surface Hubbard energy routing remains off.

Under C5.3 item 2, the requested Stage RE surface smoke outcome remains uncovered; this record does not claim a smoke failure or adoption. See the updated OPEN entry in `tests/KNOWN_ISSUES.md`.
