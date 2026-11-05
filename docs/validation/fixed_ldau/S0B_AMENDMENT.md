S0b (user-approved; record verbatim in docs/validation/fixed_ldau/S0B_AMENDMENT.md).
Base: fable_v4b_u at 876063c (local). Work branch: fix/ldau-s0b.
Token rules: Sol plans in at most 1/2 page; ONE Luna worker sequentially; ONE batched audit at the
end. Python: .venv. All runs through tools/run_rslmto_checked.py.
Budget: up to 10 SHORT runs (T0) + 4 SCF runs (T1-T3), plus ladder steps 1-2 where a state must
be accepted.

T0 Complete the S0 checks on the merged tree (short runs, single-thread, 4 steps):
   - The S0 identity gate on all 6 decks:
       * no-U RS and KS: byte-identical;
       * +U RS, KS and impurity: compared by value as in S0 T3.
   - The 13 guard fixtures.
   - dm-file runtime cases:
       * a referenced dm file that is missing -> fatal with an actionable message;
       * a dm_sha256 mismatch -> fatal;
       * a legacy file with inline ldm_flatten -> read correctly.
   Any failure: fix it (smallest change) or stop. This is the S0 acceptance.

T1 Magnetic Ni host (RS, fcc, block recursion, rc=80, nsp=1).
   - Seed: a FERROMAGNETIC no-U Ni bulk potential. Use the bulk Ni_out.nml from the user's Ni(001)
     reference (ni001_no_u_candidate/Conv/) if its moment is > 0.5 muB; otherwise use the
     existing no-U Ni bulk endpoint with an initial moment set.
   - Runs:
       (a) U=0: RS-accepted (C5.2 §2.1).
       (b) U=3/J=0.7 auto (resolves to local): RS-accepted, restarted from (a).
   - GATE: both states magnetic (|moment| > 0.4 muB).
   - RECORD: both moments (expect about 0.6 muB at U=0).
   - If (b) collapses to a non-magnetic state: record it and stop. That is a physics finding to
     review, not a failure to work around.

T2 RD2 recheck on a converged Ni impurity (in the T1(b) host, U=3/J=0.7, local mode).
   - Frozen replay as in C5.4 T1:
       * a Ni-state scalar-shift fit for the allowance;
       * the A3.1 admissibility rule for controls. Spin-swap should now be admissible; report its
         predicted separation.
   - GATE:
       * carrier ratio per spin within 1 +/- 0.01;
       * q2' >= 0;
       * scalar control within the Ni fit;
       * every admissible control separates by > 10x the carrier error.

T3 RD3(a) recheck: a Ni U=0 impurity in the T1(a) host, then impurity-vs-bulk at U>0 versus U=0.
   GATE: ratio < 10 (as in C5.2).
   RECORD: the moment and d-occupation differences.

T4 KNOWN_ISSUES: close or update the "Ni gating host was non-magnetic" entry with the T1-T3
   results.
   Local --no-ff merge into fable_v4b_u. Rerun the S0 identity gate on the merged tree.
   Report the SHA. Do not push.

Uncovered outcomes: record, leave production unchanged, continue. Tooling failures are fixed, not
reported (C6.3).
