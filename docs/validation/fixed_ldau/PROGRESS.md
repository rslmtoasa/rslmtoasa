# Fixed LDA+U audit progress

| Phase | Claim | Test/evidence | Result | Commit | Blocker |
|---|---|---|---|---|---|
| 0 | Provenance, dependencies, and unmodified baseline recorded | Debug serial build; 110 unit tests; bcc-Fe RS and reciprocal examples; VAL-04 U sweep; see [PHASE0.md](PHASE0.md) | Build, all unit tests, RS Hubbard smoke, and VAL-04 pass. Two unrelated no-U reciprocal examples reproduce reference mismatches on this macOS/GNU Fortran environment. Root accepted the categorized baseline discrepancies. | `6b6bd8a` | — |
| 1 | Independent orbital-basis and Coulomb-tensor reference | [PHASE1.md](PHASE1.md); 111/111 CTest run; SymPy comparison of 2,118 p/d tensor elements; VAL-04 rerun | Tensor/HCpx map, real-symmetric active contraction, complex tensor covariance, rotation, and nonzero-J finite-difference checks pass. End-to-end imaginary occupation handling remains Phase 2. VAL-04 passes; finite-iteration limits remain as recorded in Phase 0. | this commit | — |
| 2 | Complex Hermitian occupations and Hamiltonian potential | Not started | — | — | — |
| 3 | Variational consistency and total-energy accounting | Not started | — | — | — |
| 4 | Hubbard convergence criteria | Not started | — | — | — |
| 5 | Physical validation and regression | Not started | — | — | — |
| 6 | Documentation and release cleanup | Not started | — | — | — |
