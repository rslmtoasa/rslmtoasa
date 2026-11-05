# A4.4: Route the validated fixed-U/J energy in eligible KS runs

The fixed-U/J energy is now added to `physical_total_energy` from the production Hubbard evaluator applied to the direct output density matrix. Automatic local moments and energy adoption are limited to ordinary native reciprocal, scalar-collinear, clean-bulk KS with fixed Liechtenstein U/J and no constraints or unsupported features. RS and other ineligible paths remain spectral and unadopted.

The A5 KS M2.2/M2.3 result supports candidate b at the preregistered charge-only ASA accuracy; candidate c fails the allowance and candidate a remains report-only. Serial route, restart, explicit spectral, zero-U identity, test-field, RS fallback, matched serial/MPI, and focused unit controls pass. A constraints-enabled input confirms spectral/no-adoption resolution before an unrelated zero-reference setup failure; successful constrained SCF behavior is untested. The combined evidence and scope are indexed in `LOCAL_MOMENTS.md` and `evidence/a4_4_control_evidence_index.json`.

Validation: final-source Debug `UnitLdaUHamiltonian` and `UnitFermiDiracEntropy` pass 2/2; Python runner compiles; `git diff --check` passes. Serial and MPI binaries were built separately from the same final source.
