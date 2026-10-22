# C5.10 results

**Disposition: T1 accepted under the C5.10 amendment; T2 failed within its two-start budget. The impurity auto candidate was not adopted, T3 identity gates were stopped, and no merge was performed.** The candidate code and tests were preserved outside the repository at `/private/tmp/c510-candidate-preserved/`, bound to HEAD `141b64461a0c5398ebbf5c2b7110a9342cf86d01`; production files were returned to that HEAD.

The accepted C5.8 KS Fe host remains the T2 host. The new endpoint EF `-0.0829496275 Ry` is recorded as the future reference; the prior EF comparison is N/A under C5.10. The RS U=2/J=0, rc=80 campaign exhausted both starts without an evaluable 20-iteration window:

1. `/private/tmp/c510-t2-rs-fe-u2/run2` passed C5.4 preflight but failed during hamiltonian namelist parsing, before any SCF iteration. The serialized input had `hubbard_u_impurity(:,1)` / `hubbard_j_impurity(:,1)` instead of the intended `(1,3)` channel indices. Exit code 1; no complete trace/state pair.
2. `/private/tmp/c510-t2-rs-fe-u2/run4` passed fresh preflight with explicit `(1,3)` channel indices and entered SCF. Its lattice deck set `ct` to 6.0 per type, while the Fe rc=80 recipe uses 3.0 per type (`2*3.0d0`). The solver reported repeated missing-neighbour vectors and outer Diff values `0.1967289821`, `0.2034926014`, and `0.2010702253`; the run was terminated after three iterations. No acceptance window was produced.

Both starts count against the 2/2 budget. T2 is not accepted; no additional campaign run or continuation was launched. Because T2 failed, the remaining T3 identity checks were stopped under the C5 failure rule. No criterion changed.

The separate required RD4 lambda-zero Release short invocation passed as a bounded crash regression: the frozen C5.8 binary (`d87479d402f959b2fe80975fdb8a25afa2a1794ecae1f58146b8c5f6cbf26f7e`), immutable deck, matching preflight schema, and W1 watcher were reused in `/private/tmp/c510-rd4-lambda0-rerun/`. Preflight passed. The watcher captured one stable complete trace/state pair (six P_OUT ql channels and six POST_SCF_RADIAL channels); the authorized stop returned 143. No crash occurred through iteration 1. This is not an RD4 energy measurement.

Disk stayed above the 1 GiB floor. The bounded run’s recorded minimum was `8,212,520 KiB`; launcher completion recorded `8,412,008,448 bytes` free. No cleanup was performed.

The private candidate build passed `UnitAutoImpurityLocalMode`; binary SHA-256 was `267a8ae30a3415f8af2204be30b81bee319fa98c4e60075c1d9d5e7ea6f4f05a`. That build and its candidate patch are retained only as diagnostic artifacts. Production adoption, merge, and push remain NOT PERFORMED.
