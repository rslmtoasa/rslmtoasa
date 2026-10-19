Protocol amendment C2 (user-approved; extends C1, does not replace it; record verbatim).

1. CRASH TRIAGE (replaces C1 §5(c) as a stop condition)
On any crash or non-finite value:
  (i)   rerun once with the Debug/bounds-checked build and capture the file:line backtrace;
  (ii)  rerun the same deck with the clean 8b8f7ef binary;
  (iii) rerun with U=J=0.
If it also crashes on 8b8f7ef or at U=0, record it as a PRE-EXISTING ISSUE (backtrace, deck,
one-paragraph hypothesis), substitute the system as specified below, and continue.
STOP only if the crash is new relative to 8b8f7ef, appears only with U>0, and lies in code this
branch changed.

2. CO SUBSTITUTION (decided)
  - Drop hcp Co from RA4. Record the hcp-Co impurity-setup SIGTRAP as a pre-existing issue,
    after running the triage in §1.
  - Use fcc Ni instead (crystal_sym='fcc', no manual lattice), U=3.0 eV, J=0.7 eV, nsp=1,
    spectral mode:
      * a matched U=0 impurity-vs-bulk control;
      * U>0 impurity vs bulk, starting from both a no-U seed and a +U seed.
  - Same RECORD table format as C1 §2.3.

3. FE EMBEDDING DISCREPANCY (new priority item, done before RB)
Finding: with U=3.8/J=1.1 eV on host and impurity, the Fe impurity differs from RS bulk by
-4.46e-3 muB (300-1200x the U=0 control), and both seeds agree. V. Frilen's 2025 report shows the
same sign and size (-6.8e-3 muB). Treat this as a suspected impurity-path inconsistency.

D1 Frozen first-iteration operator comparison:
   - Start the impurity run from the CONVERGED BULK +U potential, both for the host and as the
     impurity seed (so the impurity seed occupation matrix equals the bulk one).
   - Run 1 iteration and dump, for the impurity atom and for the bulk atom of the same type:
       * the U and J values in Ry actually used;
       * the Hubbard d block inserted into the cluster Hamiltonian (hall) and into the bulk
         Hamiltonian (ee);
       * the number of times each block is inserted per Hamiltonian build;
       * the occupation matrix that was used to build each block.
   - Expectation: identical to 1e-12.
   - Also dump the first-iteration output occupation matrix of the impurity, P_out. Expectation:
     equal to the bulk P_out within the U=0 control level (~1e-5).

D2 Separate the contributions. Same frozen procedure, cases:
   (a) J=0, U=3.8;
   (b) U on the host only;
   (c) U on the impurity only (hubbard_u_impurity) versus U on the impurity type set via
       hubbard_u_general.
   Report which case carries the discrepancy.

D3 Decision rule:
   - If any inserted block, U/J value or insertion count differs: it is a bug. Fix it with the
     smallest change. This is allowed, because it restores consistency between the bulk and
     impurity paths and does not change any physics formula. Add a regression test (impurity in a
     same-material host reproduces the first-iteration bulk Hubbard block exactly). Rerun the
     Fe U>0 RECORD and require the discrepancy ratio to the U=0 control to fall below 10.
   - If all blocks match but P_out differs: localise where it diverges (Green's-function
     projection, or occupation-matrix extraction for recursion atoms versus bulk atoms). Report
     the localisation and continue to RB without a fix.

4. THEN continue C1 Stage RB exactly as specified (RB0 through RB6).
   The Fe no-U-seed U+ run goes through the C1 §1.3 ladder; a 5.47e-7 miss is not a concern.

STOP conditions are C1 §5 (a), (b) and (d), plus C2 §1. Otherwise finish RB and report.
