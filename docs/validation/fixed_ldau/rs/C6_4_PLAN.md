# C6.4 plan

1. Fix the embedded local-moment site/type mapping and eligibility in `bands.f90`; extend the oracle test for frozen Ni, ES recursion, and U-carrying recursion. Run the new test and all C6.2 T1 gates in `.venv`; stop on failure.
2. Audit copied T2/T3 decks with only authorized changes, including `nsp=1`, `nstep=200`, and accepted same-material host mixing parameters. Before each of four SCF runs, run C5.4 preflight with the matching host-E_F check.
3. Record mode, provenance, convergence, and layer observables; update workflow and open issues. Merge locally with `--no-ff` and perform the C6.1 merged-tree checks only if gates pass; do not push.
