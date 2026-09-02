# TDDFT-15 — Savrasov Sternheimer and Liouville-Lanczos design study

## Agent mission

Produce an implementation-ready design for a future fourth TD-DFT engine based on Savrasov-style Sternheimer response, updated with modern Liouville-Lanczos TDDFPT ideas. Do not implement the full solver in this task.

## Context

Savrasov 1998 is unusually relevant because it is a muffin-tin-orbital, variational time-dependent linear-response formulation. Modern TDDFPT shows how resonant/antiresonant Sternheimer equations and Liouville-Lanczos can avoid explicit empty-state sums and obtain many frequencies efficiently.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.

## Work plan

1. Derive the resonant and antiresonant first-order equations in the actual orthogonal RS-LMTO Hamiltonian notation.
2. Identify metallic occupation/Fermi-surface terms required for Fe/Ni.
3. Map the self-consistent induced charge/magnetization to the existing XC/kernel infrastructure.
4. Determine which existing linear solvers, recursion/GF kernels or preconditioners can be reused.
5. Compare cost models for frequency-by-frequency Sternheimer and Liouville-Lanczos.
6. Explain how the same response conventions/validation harness would be reused.
7. Include a future Hubbard-response path informed by Binci/Marzari/Timrov without pretending current Hubbard TD-DFT is supported.
8. Define a minimal proof-of-concept fixture and go/no-go criteria.
9. Record risks: non-Hermitian Liouvillian, metallic degeneracies, memory, preconditioning, ground-state consistency.

## Acceptance checklist

- [x] RS-LMTO Sternheimer equations are derived.
- [x] Metallic terms are identified.
- [x] Reuse of existing infrastructure is mapped.
- [x] Sternheimer vs LL cost model is documented.
- [x] Shared validation plan is explicit.
- [x] Hubbard extension is design-only and clearly separated.
- [x] Minimal prototype scope and go/no-go criteria are defined.
- [x] No premature production implementation was added.

## Required evidence / deliverables

Create `docs/TDDFT_STERNHEIMER_LL_BLUEPRINT.md` with references to Savrasov 1998 and Binci/Marzari/Timrov 2025.

## One-line commit message

`docs: design Sternheimer and Liouville-Lanczos TDDFT`
