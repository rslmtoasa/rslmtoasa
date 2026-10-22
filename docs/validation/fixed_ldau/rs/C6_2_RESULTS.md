# C6.2 — T1 gate rerun disposition

T1 was run first on branch `fix/ldau-surface` at tip
`945ff47f740d51ecdd46fdbeaa510e4f7c7b9745`. The exact C5.13 RS identity command was:

```text
python3 tests/validation/fixed_ldau/c5_7_no_u_identity.py --binary /private/tmp/c511-release-build/bin/rslmto.x --scratch /private/tmp/c62_t1_identity
```

The script's required preflight failed before launching a solver: system `python3` could not
import `f90nml` (`ModuleNotFoundError` from `C5_4_PREFLIGHT.py`). RS identity count: 0/27;
KS identity: not started (0/24); guard fixtures: not started (0/13); impurity auto fixture and
surface auto unit test: not run due the mandatory stop-on-first-failure rule. SCF launches: 0.
Failure artifacts are under `/private/tmp/c62_t1_identity/rs_candidate`.

This is an execution-environment failure, not a code or physics result. T2/T3, all SCF runs,
and T4 merge checks were not started. No merge or push was made. Per C6.2, the workflow stops
here until the required preflight dependency is available; no prior C5.13 evidence substitutes
for these mandatory reruns.
