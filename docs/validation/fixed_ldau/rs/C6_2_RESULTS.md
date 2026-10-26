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

This initial attempt failed for an execution-environment reason, not a code or physics reason.
T2/T3 and T4 had not been started at this point. At the user's direction, the required gates
were rerun with p311 below; prior C5.13 evidence did not substitute for those reruns.

## T1 retry in p311

At the user's direction, the exact gates were rerun sequentially with
`/Users/andersb/envs/p311/bin/python3` (which imports `f90nml` 1.5.0), on branch tip
`c37990999ace5870bfe0e3a0c5e1b85fa8edf53e`. All gates passed:

- RS identity: 27/27; KS identity: 24/24, using
  `tests/validation/fixed_ldau/c5_7_no_u_identity.py` and the single-thread caps in the script.
- RA2 guards: 13/13, using `tests/validation/fixed_ldau/rs_embedding_ra2.py`.
- Impurity and surface auto fixtures: `UnitAutoImpurityLocalMode` passed 1/1.

The required C5.4 preflight ran before every solver invocation. Validation made 15 solver
launches (two identity and 13 guard); the unit target launched none. Outputs are in
`/private/tmp/c62_2_t1_identity/` and
`/private/tmp/c62_2_t1_guards/RA2_provenance.json`. These are validation launches, separate
from the four-run T2/T3 SCF budget, which was 0/4 at this point. The initial system-Python
preflight failure above is superseded by this successful required rerun.

## T2 Ni(001) RS host

The staged copies of `tests/scf/cases/surface/fccCu001/input.nml` were audited before launch.
The only input changes were `alat` 3.614 -> 3.52340 and `wav` 1.41237 -> 1.37692 from accepted
Ni bulk; labels Cu/Cu-S -> Ni-imp/Ni-imp-S; `nsp` 2 -> 1 per C6.2; `fermi` -0.088376 ->
-0.061141 and `energy_max` 1.2 -> 0.5 (minimum remains -1.0); and Ni d-channel U/J on indices
(1,3) and (3,3): 0/0 for T2(a), 3.0/0.7 for T2(b). Auto radial moments were set for T2(b).
The `nstep=1` fixture setting and all other settings were preserved. Accepted Ni endpoint
potentials seed host and surface/subsurface types, with only the potential symbol changed to
match the surface label; ES stayed byte-identical. Surface positions are generated from
fractional coordinates and scale with `alat`, so geometry was not edited. Full diffs and hashes
are in `/private/tmp/c62_2_decks_audit_02/`.

For T2(a), C5.4 preflight passed with host/deck E_F delta 0 Ry against the accepted Ni U=0
endpoint (`preflight.json`). The single launch returned effectively 0 and wrote outputs under
`T2a_Ni_RS_U0/`, but the one-step run reported `Not converged! Diff=0.0162421175`. The emitted
active Ni-imp-S site had d occupations 4.701513 per spin (total 9.403026), moment 0.000000
μB. No separate frozen Ni host layer occupation or moment was emitted. This is a baseline
observation, not a smoke result.

For T2(b), C5.4 preflight passed with host/deck E_F delta 0 Ry against the accepted RS-local
Ni U=3/J=0.7 endpoint (`T2b_Ni_RS_U3_J07_auto/preflight.json`). The one launch returned 1.
The log selected `mode=local source=default`, reported no provenance mismatch, and then stopped
with fatal `Unable to map local d-shell Hubbard potential.` No outer Diff or completed iteration
window was produced. Pre-fatal initial d occupations for the frozen and active Ni types were
4.738761 per spin (total 9.477523), implying zero spin moment; these are initialization values,
not converged layer moments. Thus the T2 smoke gate failed. There is no valid T2(b) surface-minus-
bulk result or U-dependent comparison.

T2 consumed two SCF launches of the four-run budget. Per the gate-failure stop rule, T3 and T4
were not started; no merge or push was made. C6.3 is needed to investigate the local d-shell
mapping failure and the unchanged one-step surface fixture limit before further smoke validation.
