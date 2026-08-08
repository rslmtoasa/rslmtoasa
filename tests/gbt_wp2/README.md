# WP2 canonical reciprocal EBAND validation

These production-path cases share one 4x4x4 full k mesh and differ only in
the q=0 cone angle or the diagnostic DOS grid/window.  Run each with
`tests/run_test.py`, then compare the first data row of `frozen_magnon.dat`.

Acceptance budgets:

- q=0 fixed-potential global rotation: acoustic auto-branch `|omega| <= 1e-10 Ry`;
- DOS controls: the unit evaluator budget is `2e-10 Ry`; production narrow/wide
  probes must agree within their `1e-8 Ry` text-output resolution;
- every `Canonical k-space occupations` log line: `|dN| <= 1e-9` electrons;
- the converged total-DOS oracle in `UnitKspaceOccupations`: `|delta EBAND| <= 2e-5 Ry`.

The narrow DOS window is intentionally incomplete.  Its DOS-integrated and
projected `m0*m1` diagnostics may be poor; canonical EBAND must not move.

The production fixed-potential rotation fixture is necessary but not sufficient
for G2O.  In the intended rotating-frame representation a one-sublattice common
cone at q=0 is the identical collinear operator.  Pair this fixture with a
finite-q or relative-sublattice operator assertion: the current experimental
`ham0m_nc` comments can produce zero at every q because no layer consumes the
probe angles.  Neither the old double-rotation failure nor this no-op may be
promoted to a golden reference.
