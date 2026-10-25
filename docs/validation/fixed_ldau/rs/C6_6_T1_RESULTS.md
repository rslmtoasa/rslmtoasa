# C6.6 T1 gate results

**Disposition: STOP.** The mandatory current-tip C6.2/C6.4 gates were started before any
campaign run. The first attempt failed before launch because the identity script expected
disposable input files under `/private/tmp`. Those inputs were recovered from the repository's
archived C5.2 Part D package, whose deck/restart hashes match the prior successful C5.7 run. The
preflight then exposed its own obsolete import from `/private/tmp/c52-c-private`; it now imports
the archived schema module from the repository. With both setup issues fixed, preflight passed,
the RS solver exited 0, but the exact artifact identity gate failed: 12/27 artifacts differ.
At the user's direction, the remaining independent gates were run despite the RS failure.

| Gate | Result | Detail |
|---|---|---|
| RS identity | FAIL | Three current diagnostic runs exited 0; all show 12/27 artifacts differing from the registered baseline. Mismatches: `Fe_1_spinene.out`, `Fe_dos.out`, `Fe_orbene.out`, `Fe_orbital_dos.out`, `Fe_orbquadene.out`, `Fe_out.nml`, `Fe_scf.nml`, `input.nml`, `magneticdos.out`, `radial_ground_state_Fe_1.dat`, `report.out`, `totaldos.out`. |
| KS identity | PASS | 24/24 artifacts match the registered KS baseline using the current no-LibXC build. |
| 13 guard fixtures | PASS | RA2 reports all 13 fixtures passed; per-run records are under `/private/tmp/c66-current-guards/`. |
| Impurity auto fixture | PASS | `UnitAutoImpurityLocalMode` passed; it exercises the impurity auto fixture. |
| Surface auto fixture | PASS | `UnitAutoImpurityLocalMode` passed; it also exercises the surface auto fixture. |
| C6.4 frozen-plus-recursive-site test | PASS | `UnitLdaUProjectorMultisiteFrozenH` passed. |

The extended diagnostic used three RS identity runs (current LibXC build, current no-LibXC build, and the retained `c511-release-build` binary path), one KS identity run, and 13 RA2 fixture runs. The two unit tests launched no solver runs. No surface campaign probe run was launched.
The response-probe gate selected in C6.6 is `Diff <= 5e-8` within at most two steps; it was not
executed because RS identity still fails. The other T1 gates pass, but the RS identity mismatch
must be resolved before the surface probe can proceed. No oracle was changed.

## Follow-up diagnosis

At the user's direction, the current branch was rebuilt with `ENABLE_LIBXC=OFF`
and the RS identity deck was rerun through the checked preflight wrapper. It
still fails 12/27 artifacts, with the same Fermi level, moment, Diff, and
mismatch set as the LibXC-enabled run; the two current-branch output sets match
each other byte-for-byte. LibXC is therefore ruled out as the cause of this
observed mismatch. The root cause remains unresolved. The user then requested running the
remaining gates; their results are recorded above. Details and hashes are in
[`C6_6_RS_DIAG.md`](C6_6_RS_DIAG.md).
