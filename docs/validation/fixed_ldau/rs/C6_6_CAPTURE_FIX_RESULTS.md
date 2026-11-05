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

## Follow-up: mixer isolation

A separate checked U=0 diagnostic changed only `&mix.beta` from the copied bulk value `0.1` to
the original surface fixture value `0.0020`; the input still used bulk potentials as the active
surface seed and remained capped at two iterations. Step 1 Diff was unchanged (`0.0162421175`),
and step 2 Diff was `0.3838386972`, close to the beta=0.1 control (`0.3840255893`). The raw
active-Ni charge transfer was `+1.505839 e` on step 1 and `−9.969651 e` on step 2; the external
layer charge warning appeared on both steps. Lower beta therefore does not remove the response
that fails this gate.

This isolates the immediate failure to the surface charge/density response before convergence,
not the Hubbard correction or this beta value. The fixed-bulk potential seed is not a relaxed
Ni(001) surface seed, and two iterations do not settle that surface response. Because neither
this probe nor the matching bulk-seeded U=0 control reached the requested threshold, it does not
determine whether the surface density response later converges; that remains OPEN. The diagnostic
deck and matching preflight/output evidence are archived under
`C6_6_capture_fix_probe/beta_isolation/`.

The attempted pre-fix-capture run is archived under `C6_6_capture_fix_probe/initial_missing_capture/`;
the successful checked rerun and outputs are under `C6_6_capture_fix_probe/after_capture_fix/`.
`SHA256SUMS` covers both evidence sets. Build: `/private/tmp/c66-nolibxc-build`; executable SHA256
`f27e6ac2873ee8c5583d284a0881540e0a59e1df562171c98632f385dbfa4e14`.

## Disposition

The surface exact-operand capture defect is corrected. The two-step Diff failure is established to
occur in the U=0 surface density response and is not removed by lowering beta; longer-run surface
convergence remains unresolved. No further SCF runs, Cu workflow, documentation adoption, or merge
are performed.
