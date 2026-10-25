# C6.6 T1 gate results

**Disposition: STOP.** The mandatory current-tip C6.2/C6.4 gates were started before any
campaign run. The first attempt failed before launch because the identity script expected
disposable input files under `/private/tmp`. Those inputs were recovered from the repository's
archived C5.2 Part D package, whose deck/restart hashes match the prior successful C5.7 run. The
preflight then exposed its own obsolete import from `/private/tmp/c52-c-private`; it now imports
the archived schema module from the repository. With both setup issues fixed, preflight passed,
the RS solver exited 0, but the exact artifact identity gate failed: 12/27 artifacts differ.
The gate stops here, before KS or any other fixture, as instructed.

| Gate | Result | Detail |
|---|---|---|
| RS identity | FAIL | Initial run and one user-authorized no-LibXC diagnostic run exited 0; both show 12/27 artifacts differing from the registered baseline. Mismatches: `Fe_1_spinene.out`, `Fe_dos.out`, `Fe_orbene.out`, `Fe_orbital_dos.out`, `Fe_orbquadene.out`, `Fe_out.nml`, `Fe_scf.nml`, `input.nml`, `magneticdos.out`, `radial_ground_state_Fe_1.dat`, `report.out`, `totaldos.out`. |
| KS identity | NOT RUN | Stopped after the RS identity mismatch; 0/24 identity cases completed. |
| 13 guard fixtures | NOT RUN | Stopped after the RS identity failure, as instructed. |
| Impurity auto fixture | NOT RUN | Stopped after the RS identity failure, as instructed. |
| Surface auto fixture | NOT RUN | Stopped after the RS identity failure, as instructed. |
| C6.4 frozen-plus-recursive-site test | NOT RUN | Stopped after the RS identity failure, as instructed. |

Two solver runs were launched for the RS identity gate/diagnosis; no surface campaign probe run was launched.
The response-probe gate selected in C6.6 is `Diff <= 5e-8` within at most two steps; it was not
executed because T1 failed. The setup gaps are fixed, but the RS identity mismatch must be
resolved before the remaining T1 gates and surface probe can proceed.

## Follow-up diagnosis

At the user's direction, the current branch was rebuilt with `ENABLE_LIBXC=OFF`
and the RS identity deck was rerun through the checked preflight wrapper. It
still fails 12/27 artifacts, with the same Fermi level, moment, Diff, and
mismatch set as the LibXC-enabled run; the two current-branch output sets match
each other byte-for-byte. LibXC is therefore ruled out as the cause of this
observed mismatch. The root cause remains unresolved, and the stop disposition
still applies. Details and hashes are in [`C6_6_RS_DIAG.md`](C6_6_RS_DIAG.md).
