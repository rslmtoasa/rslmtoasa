# C5.5 consolidated audit

Date: 2026-10-21. Scope: one final audit of the C5.5 packet and available raw evidence. No solver or verification run was launched during this audit.

## Verdict

| Item | Audit verdict | Basis |
|---|---|---|
| T1, Ni RS-local host with default mode | PASS | Raw `testrun.log` reports `mode=local source=default use_kspace=.false. local_energy_adopted=.false.` The iteration-23 trace and restart hashes match the captured files. The first accepted window is iterations 4–23; Diff is 1.78×10⁻⁸, d-occupation spreads are 1.940081198×10⁻⁷ and 1.178661817×10⁻⁷, and sphere-moment spread is 2.3691291473×10⁻⁷. Differences from the accepted RD3(a) endpoint are 1.275531×10⁻¹⁶, 1.2406×10⁻¹⁶, and 2.3175594×10⁻¹⁶, each below its measured U=0 control allowance. |
| T2, KS U=2/J=0 Fe host | BLOCKED / UNEVALUATED before launch | The accepted input requests local, but the preserved U-carrying `Fe_out.nml` and `Fe_scf.nml` contain no serialized radial-mode field. Both files hash to `bb647254…ead6b`. Candidate parsing defaults the run to `auto`; the candidate resolver treats a missing host mode as spectral, and the unchanged validator warns that a legacy U-carrying potential without the field is spectral. This establishes a predicted fallback, not a measured model result: no T2 SCF or model log exists. |
| T3 verification and adoption | NOT RUN | T2 did not pass; C5.5 requires stopping before adoption if T1 or T2 fails to pass. The candidate unit result is separately recorded as PASS, but it does not satisfy T3's production verification gates. |
| Production adoption and local merge | NOT PERFORMED | No production source or documentation changes for adoption, no commit, and no merge are present. |

T1 is the only campaign attempt used (1 of 3); T2 used zero solver attempts. As clarified by the user, verification runs are tracked separately from the three-campaign cap. The candidate resolver unit is recorded as a separate verification. The exact candidate patch remains isolated and unapplied.

## Provenance and scope checks

- Independently checked the T1 raw mode line, preflight record, acceptance scan, endpoint trace/restart hashes, and raw restart hashes for the T2 KS host. T1's preflight identifies an RS run (`use_kspace=false`) and the captured trace/state are the ones named by its scan.
- Inspected the candidate resolver and provenance validator. A U-carrying host file without a mode field maps to spectral; forcing local would then hit the mode-mismatch fatal path. No host retagging or resolver relaxation is justified by this evidence.
- Production source, `tests/KNOWN_ISSUES.md`, and documentation outside the C5.5 packet have no working-tree changes. User-owned hcp Co files remain outside this audit's write scope.
- Disk sample during audit: `/private/tmp` reported 6,542,404 KiB available (about 6.24 GiB). This is an observation only; no cleanup was performed and no new disk threshold was applied.

## Final status

T1 passed. T2 is blocked before measurement because the accepted KS endpoint lacks the host-mode provenance required by the unchanged auto rule. Therefore T3 adoption and C5.2 §7 merge checks remain unperformed. The implementation candidate is not adopted, and the requested local merge is not ready to report as complete.
