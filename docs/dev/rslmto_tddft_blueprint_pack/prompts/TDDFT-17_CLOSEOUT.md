# TDDFT-17 — documentation, provenance, and production closeout

## Agent mission

Close the campaign by making the validated physics, supported feature matrix, backend choices and convergence guidance discoverable to users and developers.

## Context

A multi-backend response implementation is only useful if users can tell which algorithm they are running, what is physical versus numerical broadening, and which combinations are validated.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.

## Work plan

1. Consolidate user documentation for `eigenpairs`, `gf_k_lehmann`, and `gf_r_native`.
2. Document when each backend is expected to be attractive: transparent validation, periodic k-space work, dense-q/native real-space/surface workflows.
3. Add a support matrix for collinear/SOC, overlap, Hubbard, GBT, CCOR and projections.
4. Document convergence knobs separately for each backend.
5. Add reproducible Fe/Ni examples and point to validation evidence.
6. Add developer notes on Ward/Goldstone identities and why empirical shifts are forbidden.
7. Ensure runtime output contains complete provenance.
8. Run the full relevant test suite and archive the final validation/performance summary.

## Acceptance checklist

- [x] User backend documentation complete.
- [x] Support matrix complete.
- [x] Fe/Ni examples reproducible.
- [x] Convergence guidance backend-specific.
- [x] Ward/Goldstone developer note present.
- [x] Runtime provenance complete.
- [x] Full relevant test suite passes.
- [x] Final status explicitly states what is production-validated and what remains experimental.

## Required evidence / deliverables

Create/update the final TD-DFT documentation and a short `docs/TDDFT_RELEASE_STATUS.md`.

## Closeout evidence record

The final user documentation is
[`docs/TDDFT_USER_GUIDE.md`](../../../TDDFT_USER_GUIDE.md). It consolidates
the three backend choices, backend-specific convergence controls, support
matrix, Fe/Ni deck commands, output provenance, and links to the detailed
developer notes.

The release decision is recorded in
[`docs/TDDFT_RELEASE_STATUS.md`](../../../TDDFT_RELEASE_STATUS.md). It
separates deterministic fixture validation from the still-failing Fe/Ni
material gate and identifies the remaining experimental/deferred routes.

The production output metadata block is emitted by
`append_tddft_metadata` in `source/calculation.f90`. It records requested and
canonical backend identity, implementation, q/k/R-space controls, broadening
and contour roles, source/measurement/kernel provenance, Goldstone policy,
unsupported-feature policy, and MPI decomposition provenance. The source
contract is guarded by `tests/unit/test_tddft_dispatch.py`.

The Fe/Ni examples are the three current-branch decks under
`tests/regression/tddft_validation/materials/{bccFe,fccNi}`. Their committed
raw outputs and `results/validation/TDDFT-11_FE_NI/evidence.json` are evidence
of reproducible bounded probes, not a material release pass. The checked-in
TDDFT-11 report records the backend disagreement, nonzero raw Ward residual,
missing R-tail sweep, and Ni connectivity limitation.

## One-line commit message

`tddft: close out validated linear-response framework`
