# C5 Part II — R1 findings and stop decision

**R1 setup GATE: PASS.** The exact corrected hcp-Co Debug deck reproduced a source sizing defect before SCF: `ctnew` had extent 3 while the cutoff vector had extent 50000. `source/lattice_cluster.f90` now allocates the scaled copy to `size(this%ct)`. The Debug regression passes and emits a compact cluster record with `ntype=3`, `nbas=13`, `nrec=1`.

| Check | Result | Evidence |
|---|---|---|
| Pre-fix hcp-Co Debug reproduction | Reproduced `Array bound mismatch ... ctnew (3/50000)` at line 996, with full backtrace | `C5_PART_II_PROVENANCE.json`, raw log SHA256 `f2d2415a…` |
| Post-fix cluster setup regression | **PASS** | `Val04RsHcpCoClusterSetup`; `ntype=3`, `nbas=13`, `nrec=1` |
| U=0 Co impurity RECORD | **NOT RUN — invalid restart source** | Corrected deck had only `.nml` files that the loader identifies as fresh atomic starts; no authorized converged `_out.nml` endpoint was available. The earlier stale-template attempt exited before SCF and is recorded as `INVALID_SETUP`, not a failed gate. |
| No-U RS standard fixture identity | **FAIL** | Matched Release/manifest pair differs in scientific outputs, including `Fe_out.nml`, `Fe_scf.nml`, DOS and radial files. |
| No-U KS standard fixture identity | **PASS** | Scientific outputs are byte-identical. `str.out` differs only in the `cpu_s` timing field; `testrun.log` is excluded under the existing comparison convention. |

The no-U comparisons used clean `8b8f7ef` and frozen R1 Release binaries, identical manifest-resolved inputs, and matching GNU Fortran 16.2/OpenMP configurations. The RS differences are small but nonzero; for example, `lmom(1)` is `4.4707783658134787e-12` on clean 8b and `4.4723755520708659e-12` on R1. The strict byte-identity criterion has no numeric tolerance.

**Decision:** STOP under C5 §2.5(a) because RS no-U scientific outputs differ byte-for-byte. The cause is not isolated: same-binary repeatability was not tested, so run-to-run variation remains an unverified hypothesis. No additional runs, R2 campaign work, RD implementation, or merge were performed after the stop verdict.

**Interpretation (8 lines):**

The hcp-Co crash was an array-capacity defect, not an input typo or an LDA+U physics issue. The fix retains the full cutoff vector while the emitted cluster still reports only its active topology. The R1 setup gate passes. The requested Co U=0 record was not launched because its available potentials would have been fresh starts under the loader. The matched KS fixture passes byte identity apart from recorded timing metadata. The matched RS fixture has small scientific-output differences, which still fail the exact criterion. No tolerance was changed. A same-binary determinism check could help identify the cause if the run budget permits later.

Passed: R1 cluster setup GATE and KS no-U identity check.
Recorded: RS no-U identity failure; one invalid U=0 setup attempt; control not run because no valid restart.
Provisional: cause of the small RS output differences is not isolated.
Stopped: C5 §2.5(a); no further campaigns or RD implementation.
Next: preserve the stop verdict for user review; do not merge or resume downstream stages in this run.
