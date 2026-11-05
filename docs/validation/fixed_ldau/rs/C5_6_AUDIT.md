# C5.6 final batched audit

**Verdict: T0 deletion evidence passes its recorded scope; T2a fails before SCF, so T2b was not run and the explicit C5.6 stop branch applies.** The campaign used 1 of 2 remaining attempts. The C5.5 clarification keeps identity, unit, and merge verification outside the campaign-run cap.

The ledger binds the 24 deleted directories to committed compact evidence. Its recorded allocated-byte sum is 250,454,016 bytes; observed free-space increase during cleanup was 242,360 KiB. All ledgered deletion paths were absent at audit. The RD2 archive hash matches and its required entries were recorded readable. The protected accepted KS Fe, Ni hosts/control, and named later reference directories remain present. The original accepted Fe `Fe_scf.nml` and `Fe_out.nml` still hash to `bb647254…ead6b`; the user Co files remain untouched.

I independently inspected the complete archived T2a log and the current validator in `source/self.f90`. Preflight passed. The runtime confirms native KS (`use_kspace=.true.`) resolved default mode `local`; the legacy U-carrying restart lacks mode metadata, is interpreted as `spectral`, and triggers the fatal mismatch before any SCF iteration. The raw log, packet, binary, source manifest, and unchanged restart hashes agree with the result records. KS acceptance is unevaluated, and no regenerated host with serialized `local`/`ks` metadata exists. This establishes a provenance-guard block for the supplied legacy checkpoint; it does not establish that a retry or retag would work.

T2b was not launched because T2a produced no accepted host. T3’s requested warning refinement, unit test, and documentation, as well as T4 merge, are stopped by the user-approved failure branch. No retries, source edits, retagging, SCF runs, cleanup, or merge were made during audit.

Disk available space was 7,547,469,824 bytes before T2a, reached a periodic minimum of 7,545,163,776 bytes, and was 7,316,240 KiB after T2a. The 1 GiB guard did not trigger. At audit time `df` reported 7,223,536 KiB available, still above the guard; external filesystem activity may contribute to the change.

1. T0 evidence and path checks: pass within the ledgered scope.
2. T2a: failed before SCF at the provenance guard; KS acceptance unevaluated.
3. T2b: not run; no regenerated accepted host.
4. T3 and T4: stopped under the explicit C5.6 failure branch; no merge SHA.
5. Campaign use: 1/2; no further campaign work authorized by this branch. See [C5_6_AUDIT.json](C5_6_AUDIT.json) for hashes and checks.
