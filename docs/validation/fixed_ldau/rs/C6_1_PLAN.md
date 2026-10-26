# C6.1 plan

1. Rerun non-SCF T1 gates; stop on any failure. Verify surface geometry scaling and stage T2/T3 copied decks with line-by-line reasons and diffs before launch.
2. Run at most four SCF invocations: Ni RS surface U=0/U=3-local, then Cu KS bulk U=3-auto and its surface smoke. Run the prescribed host-E_F preflight before each; analyze smoke against C5.13.
3. Record occupations, moments, provenance, and any unresolved outcomes. Update the workflow and OPEN issues. Merge locally with `--no-ff` and run merged-tree gates only if all required gates pass; never push.
