# Protocol amendment C1 (user-approved)

This exact amendment supersedes B1 and all earlier stop rules for the remaining RS stages.

Protocol amendment C1 (user-approved; supersedes B1 and all earlier stop rules for the remaining
RS stages; record verbatim).

PURPOSE
Finish RA4 and complete Stage RB in one run, without stopping for routine numerical snags. The
decisions below are made. Do not reopen them, do not ask for confirmation, and do not propose
alternatives unless a STOP condition (section 5) is hit.

1. DEFINITIONS (fixed)
1.1 RS-accepted state. All of the following must hold:
    (a) outer Diff <= 3e-8;
    (b) over the last 20 iterations, the spread (max minus min) of the d occupation per spin is
        <= 5e-7, and the spread of the sphere moment is <= 1e-6 muB;
    (c) fixed-U P/V residuals are logged only.
1.2 Records versus gates.
    - RA4, RB5, and any comparison described as "descriptive" are RECORDS. Measure them, tabulate
      them, and report each discrepancy next to its U=0 control. They cannot fail, and they never
      block progress.
    - Only the items marked GATE in section 3 can fail.
1.3 Run budget and fallback ladder for any SCF that is not RS-accepted at its step cap:
    (i)   continue from the last state for 200 steps with beta halved;
    (ii)  if it is still not accepted, continue for 200 steps with broyden and beta = 0.01;
    (iii) if it is still not accepted, label the final state PROVISIONAL, use it for RECORDS
          with that label, and move on.
    Never run more than these two continuations per state.

2. FINISH RA4 (RECORD)
2.1 Re-evaluate the existing runs against 1.1 using their logs. A state that qualifies is used
    as is. Apply 1.3 only to states that do not qualify.
2.2 Run, in spectral mode with nsp=1:
    - Fe: impurity vs bulk at U=3.8 eV, J=1.1 eV, starting the impurity from a no-U seed and from
      a +U seed.
    - Co: a matched U=0 impurity control, then impurity vs bulk at U=3 eV, J=0.7 eV, from both
      seeds.
2.3 Output a single table. One row per case; columns:
    - d up and d down occupations, and sphere moment, for impurity and for bulk;
    - the impurity-minus-bulk differences;
    - the matching U=0 differences;
    - the ratio of the U>0 discrepancy to the U=0 discrepancy;
    - status (ACCEPTED or PROVISIONAL).
    Also give the maximum d-LDOS difference. Then continue straight to Stage RB.

3. STAGE RB
RB0 Production RS stopping tolerances (approved change):
    - Make the fixed-U R05 P and V tolerances input-adjustable if they are not already.
    - Set RS defaults: P residual <= 5e-7, V residual <= 1e-7 Ry, streak 3.
    - Leave KS defaults unchanged.
    - GATE: no-U output byte-identical to 8b8f7ef; KS fixed-U behaviour unchanged; the Fe U=3.8
      eV RS bulk run now reports convergence.
RB1 GATE: bitwise identity for no-U, explicit U=J=0, spectral U>0 (allowing only the RB0
    stopping-rule difference), and explicit local with U=J=0 (must be a no-op).
RB2 GATE: RS bulk bcc Fe in local mode, RS-accepted per 1.1, at U=2 eV (J=0), U=4 eV, and
    U=2 eV with J=0.9 eV.
RB3 GATE: frozen local-vs-spectral replay on the U=2 eV local state.
    - Carrier ratio per spin within 1 +/- 0.01.
    - q2' >= 0.
    - Scalar control error within the A2.2-style fitted allowance.
    - The sign and spin-swap controls deviate from the true carrier by more than 10x the carrier
      error.
RB4 GATE: restart reaches RS-accepted within 10 steps; serial and 2-rank MPI d occupations agree
    to 1e-10.
RB5 RECORD: RS-vs-KS comparison in local mode at U=0, 2, 4 eV (J=0) and U=2 eV with J=0.9 eV.
    Report the U>0 differences relative to the U=0 difference.
RB6 RECORD, timebox 2 runs: Co U=0 bulk with test-only tighter atomsc inner tolerances
    (tolrsq=1e-10, tol=1e-8, inner cap 200). Report whether the outer plateau drops.

4. REPORTING (replaces the hash-heavy prose)
- One markdown page per stage containing: the result tables, GATE verdicts, and at most 10 lines
  of interpretation.
- All hashes, commands and paths go in one JSON sidecar per stage. Do not repeat them in prose.
- Close with a 5-line summary: what passed, what was recorded, what is provisional, and what
  comes next.

5. STOP ONLY IF one of these occurs
  (a) no-U byte identity breaks;
  (b) a GATE fails after one honest debugging attempt (state the hypothesis you tested);
  (c) the code crashes or produces non-finite values in a supported configuration;
  (d) continuing would require touching guarded scope (SOC, nsp>=2, non-d shells, HOH, KS
      impurity) or changing a physics formula.
Everything else is handled by the rules above: write it down and keep going. Otherwise, stop at
the end of Stage RB.
