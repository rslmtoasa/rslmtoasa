Protocol amendment C3 (user-approved; extends C1/C2; record verbatim).
Order: Part A (timeboxed), then Part B. STOP rules: C1 §5(a,b,d) and C2 §1.

PART A — Fe embedding discrepancy: localisation (timebox: 12 runs total, no production change
unless A4 finds a bug)

A1 Fermi level: confirm that the Fe U>0 impurity decks use the converged U=3.8/J=1.1 bulk EF, not
   the U=0 EF. If they used the wrong EF, rerun both Fe U>0 impurity records with the correct EF
   and report the new ratio. If the ratio falls below 10, finish Part A there.

A2 Per-site frozen Hamiltonian comparison:
   - Take the impurity run seeded from the converged bulk +U state, iteration 1.
   - For every site in the impurity cluster, compare its on-site block in the cluster
     Hamiltonian (hall) with the bulk on-site block for its type, including the Hubbard block.
     Report max|diff| per neighbour shell.
   - Compare the hopping blocks between the impurity and its first two shells with the bulk
     ones.
   - Do the same at U=0 as the control.

A3 Recursion coefficients:
   - For the impurity d orbitals (both spins) and the bulk d orbitals, dump a_n and b_n for
     n = 0 .. lld, plus the terminator parameters, at iteration 1 of the same frozen setup.
   - Report the first n at which impurity and bulk differ by more than 10x the U=0 level.
   - Then: n = 0 means on-site; n = k means shell k; no difference up to lld means the
     terminator.

A4 Decision:
   (i) If a host-site Hubbard block in hall is missing, doubled, stale, or built from the wrong
       occupation matrix or type: fix it with the smallest change; add a regression test (all
       cluster sites reproduce the bulk on-site blocks at iteration 1); rerun the Fe U>0 records.
   (ii) If the terminator differs, or the first differing level is deep (n >= 5): run Fe U>0
        impurity versus bulk at rc in {80, 100, 140} and lld in {21, 31}, and report how the
        discrepancy scales.
        - If it shrinks with rc or lld, it is a finite-cluster/terminator effect, not a bug.
          Document recommended rc/lld for +U impurity runs in lda_u_workflow.rst and close.
        - If it does not shrink, report the data and close Part A as OPEN.
   (iii) Anything else: report the localisation and close Part A as OPEN.
   Do not exceed the run timebox. Part B proceeds in every case, because it concerns bulk only.

PART B — Stage RC (RS bulk energy), decisions fixed

RC0 Establish from the code how RS SCF occupations are formed: window, channels_ldos, terminator,
    and whether any Fermi factor enters the SCF.
    - Step-function occupations: the functional is E (no entropy).
    - Fermi-smeared occupations: add an RS entropy diagnostic on the same energy mesh, validated
      on a finite-H toy model, and use F = E - TS.
    No user decision needed either way. Report which case applies.

RC1 W1 control in RS local mode at U=0: lambda in {0, +/-0.0025, +/-0.005} Ry, both step sizes.
    - Report the relative residual of dF/dlambda against Tr(W P).
    - GATE: |relative residual| <= 1 % (about 7x the KS value of 0.14 %). If it fails, STOP and
      report; Part B ends there.

RC2 HF tests, RS local mode, bcc Fe:
    - U0 = 2 eV with delta in {0.05, 0.10} eV;
    - J derivatives at J = 0.5 and 0.9 eV (U = 2 eV);
    - candidates a, b, c as in KS.
    - Allowance A = 3 x |RC1 relative residual| x |dE_Hub/dU at fixed P|, plus the finite-
      difference spread between the two deltas.
    - GATE: b within A, and c misses by more than 10 A.
    - a is a RECORD.

RC3 Adoption (only if RC2 passes):
    - 'auto' resolves to 'local' and E_Hub is routed into physical_total_energy for clean RS bulk
      (calctype='B', nbulk=0, nrec=ntype, nsp=1, lmax=2, fixed Liechtenstein U/J, no
      constraints).
    - Impurity and surface runs stay spectral, with no energy routing.
    - GATE: no-U byte identity; explicit spectral U>0 unchanged; KS behaviour unchanged.
    - Update ENERGY_ACCOUNTING.md, REMEDIATION_GATES.md (R04 scope) and lda_u_workflow.rst.

Report format as in C1 §4. Stop at the end of RC.
