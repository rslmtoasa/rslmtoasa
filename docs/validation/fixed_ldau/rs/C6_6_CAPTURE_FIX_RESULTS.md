# C6.6 surface `V_in` capture follow-up

Date: 2026-10-25. Branch: `fix/ldau-surface`. This follow-up fixes the surface capture omission
found while investigating the failed Ni(001) response probe.

## Cause and correction

`run_recursion` initialized RD1 `V_in` capture only for `calctype='I'`, although surface local
moments were eligible for the same capture. Adding `S` to that dispatch exposed a second omission:
surfaces assemble their recursion Hamiltonian through `build_bulkham` and do not call
`build_locham`, where the existing exact-operand capture hook lived. The first checked rerun
therefore stopped with status 5 (`RD1 exact V-in snapshot missing/stale`, iteration 1, recursion
site 2, cluster site 970, type 3).

The fix now initializes capture for `I` and `S`, and records the exact type-indexed operand at the
surface `ee` insertion in `build_bulkham`. The `B` path is unchanged. The checked surface rerun
completed without provenance or capture errors, establishing that the surface correction consumed
a current-iteration captured operand.

## Test output

The focused CTest cases passed after the change:

- `UnitLdaUProjectorMultisiteFrozenH`
- `UnitAutoImpurityLocalMode`

The U=3/J=0.7 Ni(001) auto probe used the archived C6.6 deck and the checked launcher. Its
preflight passed, including exact host-E_F match (`delta=0 Ry`), and the solver returned 0. The
mode was local, with no provenance error. Final two-step Diff remained `0.3869175281` (step 1:
`0.0163260002`), failing the agreed `Diff <= 5e-8` response gate. Final layer d occupations and
reported spin moments were byte-identical to the pre-fix run; the U=0 report and potentials also
showed that the same two-step Diff growth occurs without U (`0.0162421175` to `0.3840255893`).
Thus the capture defect is fixed, but it did not cause the two-step Diff failure in this deck. The
probe still does not establish convergence or a physical surface moment.

The attempted pre-fix-capture run is archived under `C6_6_capture_fix_probe/initial_missing_capture/`;
the successful checked rerun and outputs are under `C6_6_capture_fix_probe/after_capture_fix/`.
`SHA256SUMS` covers both evidence sets. Build: `/private/tmp/c66-nolibxc-build`; executable SHA256
`f27e6ac2873ee8c5583d284a0881540e0a59e1df562171c98632f385dbfa4e14`.

## Disposition

The surface exact-operand capture defect is corrected. The smoke gate remains failed, so no further
SCF runs, Cu workflow, documentation adoption, or merge are performed. The cause of the large
second-step Diff is still unresolved and must be investigated separately from the capture fix.
