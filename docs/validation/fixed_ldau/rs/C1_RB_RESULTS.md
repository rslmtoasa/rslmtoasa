# C1 Stage RB results

Status: C2-active. RB0–RB4 passed. RB5's four true-KS records and RB6's two registered Co U=0 records are complete; RB6 reproduced a pre-existing nonfinite failure, so its plateau comparison is indeterminate. Exact commands, paths, hashes, and run evidence are in [C1_RB_PROVENANCE.json](C1_RB_PROVENANCE.json).

| Stage | Result |
|---|---|
| RB0 | **PASS** — R05 RS targets are adjustable (defaults P=5e-7 and V=1e-7 Ry; streak 3); Fe U=3.8/J=1.1 RS production convergence and joined C1 history passed. KS fixed-U behavior and no-U identity passed. |
| RB1 | **PASS** — no-U, explicit U=J=0, spectral U>0, local U=J=0 no-op, and active-U true-KS fixed behavior matched their registered physical-output comparisons. |
| RB2 | **PASS** — local RS Fe U=2/J=0, U=4/J=0, and U=2/J=0.9 endpoints satisfy C1 acceptance. |
| RB3 | **PASS** — accepted U2/J0 local-versus-spectral replay and controls passed the C1 amendment gates. |
| RB4 | **PASS** — accepted-state restart qualified at k=1; serial and 2-rank MPI traces were byte-identical and d occupations agreed exactly. |
| RB5 | **RECORDS COMPLETE** — four corresponding-seed true-KS comparisons; U0 continuation has a documented same-beta procedure deviation. |
| RB6 | **RECORDS COMPLETE, INDETERMINATE** — both registered Co U=0 default/tight records failed on the same pre-existing NaN; no plateau comparison is supported. |

## RB5 central QL comparison

Both backends use the final complete P_OUT QL at each accepted terminal iteration (bulk site 0); the moment is sum over l of spin-up minus spin-down. Δ is true-KS minus matching local RS. ΔΔ is Δ at U>0 minus Δ at U=0, as required for the series comparison. These are descriptive records, with no additional pass threshold.

| Fe case | RS d↑ / d↓ / m (μB) | True-KS d↑ / d↓ / m (μB) | Δ d↑ / d↓ / m (μB) |
|---|---:|---:|---:|
| U0/J0 | 4.349200 / 2.149941 / 2.103828 | 4.301186 / 2.204595 / 2.003489 | −0.048014 / +0.054654 / −0.100339 |
| U2/J0 | 4.615579 / 1.895454 / 2.558681 | 4.612752 / 1.892685 / 2.546925 | −0.002827 / −0.002769 / −0.011756 |
| U4/J0 | 4.746883 / 1.743735 / 2.785370 | 4.753320 / 1.728713 / 2.808594 | +0.006438 / −0.015023 / +0.023224 |
| U2/J0.9 | 4.413633 / 2.087146 / 2.198745 | 4.500008 / 1.999407 / 2.366575 | +0.086375 / −0.087739 / +0.167831 |

U>0 Δ relative to the U0 Δ (difference-in-differences, ΔΔ):

| Fe case | ΔΔ d↑ | ΔΔ d↓ | ΔΔ m (μB) |
|---|---:|---:|---:|
| U2/J0 | +0.045187 | −0.057424 | +0.088584 |
| U4/J0 | +0.054452 | −0.069677 | +0.123564 |
| U2/J0.9 | +0.134389 | −0.142393 | +0.268170 |

The runtime logs confirmed `use_kspace=.true.` for all four corrected true-KS runs. U4/J0 converged at step 59. U2/J0 and U2/J0.9 reached the 183-step cap with printed terminal Diff 0.0 and stable final-20 QL windows, but without a `Converged!` marker. U0/J0 reached its initial 21-step cap with wide QL spread, then continued 79 steps from that exact state; the continuation retained beta 0.05 instead of the prescribed halved beta. Its final-20 QL window is stable, but this protocol deviation is recorded explicitly and no extra run was made. The initial wrong-backend U0 attempt (19 RS steps) is excluded.

RB6's default and tighter-ATOMSC runs each completed 19 outer Diff records before NaN mixed qia/Invalid PL for Co1; their histories agree within 1e-10 and end at Diff 0.1321778518. Clean 8b and the explicit-U=J=0 Debug replay reproduce the nonfinite issue. The bounds replay fails earlier during hcp lattice setup. No plateau reduction was observed before failure, so the tighter-ATOMSC comparison is indeterminate.

Final-source default-off no-U identity was also checked after the diagnostic-only capture-helper edit: matched clean 8b and ed6 Debug builds produced identical physical outputs for RS and true KS. Only runtime logs differed; true-KS `str.out` differed in timing text. This is implementation sanity evidence, not another formal gate.

## Interpretation (8 lines)

1. RB0–RB4 passed their registered C1 checks.
2. The final ed6 default-off RS and true-KS no-U identity checks passed against clean 8b.
3. RB5 reports matched-phase QL differences and U>0 difference-in-differences against the U0 change.
4. All four corrected RB5 runs logged the true-KS backend.
5. The U0 continuation used the same beta after its initial cap; this deviation is retained in the record.
6. U2/J0 and U2/J0.9 reached their caps without a `Converged!` marker; U4/J0 converged at step 59.
7. RB6 default and tight runs reproduced the clean-baseline nonfinite failure; no plateau conclusion follows.
8. C2 classification and diagnostic/source provenance remain in the JSON sidecar.

## C1 close (5 lines)

Passed: formal gates RB0, RB1, RB2, RB3, and RB4.
Recorded: RB5 paired true-KS QL comparisons; no provisional endpoint is claimed for cap-limited cases.
Deviation: RB5 U0 continuation retained beta; RB6 Co U=0 reproduced pre-existing NaN, with plateau indeterminate.
RA4 and D3 RECORD are complete under C2.
Stop: RB is complete; no later-stage work or commits are included.
