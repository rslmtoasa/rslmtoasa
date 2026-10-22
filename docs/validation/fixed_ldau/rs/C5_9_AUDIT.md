# C5.9 batched audit

**Verdict: T1(a), T1(c), occupations and sphere moment pass; T1(b) EF is UNEVALUATED / UNCOVERED.** Apply C5.3 item 2: do not run T2b, adopt the candidate, or merge. No T1(b) numerical failure is established.

The regenerated Release run is the recorded fixed KS Fe U=2/J=0 run: exit 0, production log reports `HUBBARD_RADIAL_MODE mode=local` and `Converged!0.0000000002`, and the final R05 record reports physical P/V streak 3. Written `Fe_scf.nml` and `Fe_out.nml` both carry `local`/`ks` (bound by SHA-256 `2284194179736d2e6694dd40514b55578bdf2558e9826d5d7ba9d5c6cc13b4bd`). Against the accepted state snapshot bound to `bb647254…ead6b`, d-up, d-down, and sphere-moment absolute differences are respectively 1.608e-9, 3.168e-9, and 6.093e-9; all are below 1e-6.

The regenerated endpoint EF is -0.0829496275 Ry. The accepted endpoint proof and state files provide no EF, and the preserved initial deck value -0.082950 Ry is not evidence of endpoint EF. Its 3.725e-7 Ry difference is therefore indicative only and cannot establish T1(b); the 1e-6 Ry comparison remains uncovered. The unavailable old raw source directory does not change this determination.

RD4 lambda=0 Release diagnostic preflight passed. The watcher captured a stable complete iteration-1 trace/state pair (both output channel counts 6) before external SIGTERM; exit 143 is the expected bounded stop, with no crash through that iteration and no energy measurement. Disk stayed above the 1 GiB floor: periodic minimum 9,255,763,968 bytes; launcher-end 9,253,056,512 bytes. The submitted results correctly record 0/2 campaign runs and one separately budgeted short invocation. No material arithmetic or gate correction is needed. No solver, build, cleanup, or merge was run for this audit.
