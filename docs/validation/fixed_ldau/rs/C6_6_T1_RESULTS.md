# C6.6 T1 gate results

**Disposition: STOP.** The mandatory current-tip C6.2/C6.4 gates were started before any
campaign run. The first gate, single-thread no-U RS/KS identity, failed before running a solver:
the script could not find `/private/tmp/c5-2-partA/a1_bccFe/source_snapshot/input.nml`.

| Gate | Result | Detail |
|---|---|---|
| RS identity | NOT RUN | Missing input snapshot at `/private/tmp/c5-2-partA/a1_bccFe/source_snapshot/input.nml`; 0/27 identity cases completed, no artifacts produced. |
| KS identity | NOT RUN | Identity script stops on the preceding RS setup error; 0/24 identity cases completed, no artifacts produced. |
| 13 guard fixtures | NOT RUN | Stopped after first T1 failure, as instructed. |
| Impurity auto fixture | NOT RUN | Stopped after first T1 failure, as instructed. |
| Surface auto fixture | NOT RUN | Stopped after first T1 failure, as instructed. |
| C6.4 frozen-plus-recursive-site test | NOT RUN | Stopped after first T1 failure, as instructed. |

No solver was launched by this gate attempt and no campaign SCF budget was consumed. The
response-probe gate selected in C6.6 is `Diff <= 5e-8` within at most two steps; it was not
executed because T1 failed. The probe remains pending restoration of the missing identity input
snapshot and a clean rerun of all mandatory T1 gates.
