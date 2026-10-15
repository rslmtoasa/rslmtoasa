# Fixed LDA+U audit progress

The controlling next-stage plan is [LOCAL_MOMENTS.md](LOCAL_MOMENTS.md),
which records the M0 scope, operator orientation, predeclared signed carrier
statistic, controls, and evidence requirements. M0 is **not passed** on branch
`fix/ldau-local-moments` at base `3d7571a`: the KS no-U scalar control exceeds
its independently measured floor-plus-refinement allowance by `2.0847x`, and
the RS no-U reference did not converge. No local-moment correction or
physical-energy change is applied.

M0.1 finite-H RS/KS M1 extraction oracles and two-site mapping controls pass
as class-B probes, both serially and in the targeted two-rank MPI tests; see
[the machine-readable measurements](evidence/m0_1_finite_operator.json).
The M0.2 KS floor/repeat/refinement measurement is recorded in
[the scalar-control failure artifact](evidence/m02_ks_scalar_v0_fail.json);
its scalar intervention fails the frozen bound. RS floor is unavailable due
to nonconvergence. M0.3 was not started, so R04 remains blocked.

| Phase | Claim | Test/evidence | Result | Commit | Blocker |
|---|---|---|---|---|---|
| 0 | Provenance, dependencies, and unmodified baseline recorded | Debug serial build; 110 unit tests; bcc-Fe RS and reciprocal examples; VAL-04 U sweep; see [PHASE0.md](PHASE0.md) | Build, all unit tests, RS Hubbard smoke, and VAL-04 pass. Two unrelated no-U reciprocal examples reproduce reference mismatches on this macOS/GNU Fortran environment. Root accepted the categorized baseline discrepancies. | `6b6bd8a` | — |
| 1 | Independent orbital-basis and Coulomb-tensor reference | [PHASE1.md](PHASE1.md); 111/111 CTest run; SymPy comparison of 2,118 p/d tensor elements; VAL-04 rerun | Tensor/HCpx map, real-symmetric active contraction, complex tensor covariance, rotation, and nonzero-J finite-difference checks pass. End-to-end imaginary occupation handling remains Phase 2. VAL-04 passes; finite-iteration limits remain as recorded in Phase 0. | this commit | — |
| 2 | Complex Hermitian occupations and Hamiltonian potential | [PHASE2.md](PHASE2.md), [MPI reproduction](../../../tests/validation/fixed_ldau/README.md) | Serial complex occupation, serialized restart/mixing, Hamiltonian derivative, scalar RS/KS and input guards pass; actual reciprocal and RS one/two-rank reductions pass. The original nsp=2 U sweep is outside the supported scalar-collinear scope; current supported nsp=1 smoke cases are finite-iteration checks. | PASS | — |
| 3 | Variational consistency and total-energy accounting | [PHASE3.md](PHASE3.md); production functional evaluator and independent SymPy energy reference | Core Eint/EDC/EHub/TrVP algebra, complex J=0/nonzero-J derivatives, d-shell rotation covariance, J=0 projector cancellation, and radial PL→NEWRHO sensitivity checks are validated. Exact physical SUMEV double-count coefficient/sign remains unresolved; no physical energy correction is applied and Gate 3 is not passed. | partial (Gate 3 held) | Production spectral-moment to radial energy identity not quantitatively closed |
| 4 | Hubbard convergence criteria | [PHASE4.md](PHASE4.md); current-source Debug 113/113; MPI residual direct+2-rank | Full complex P/V residual stopping, three-update stability, physicality rejection, and MPI global reductions pass. This does not resolve Gate E or promote a material energy claim. | implementation validated | — |
| 5 | Physical validation and regression | Not started | — | — | — |
| 6 | Documentation and release cleanup | Not started | — | — | — |

The later phase rows above retain their historical summaries. Current
Stage-M0 and R00–R15 dispositions are maintained in
[LOCAL_MOMENTS.md](LOCAL_MOMENTS.md) and
[REMEDIATION_GATES.md](REMEDIATION_GATES.md); these supersede any stale
phase-table wording that would suggest M0 or R04 passed.
