# TDDFT-16 — relativistic and noncollinear TD-DFT design boundary

## Agent mission

Define the formal and software work required before SOC/noncollinear TD-DFT can be enabled, using the relativistic Jülich literature rather than extending the collinear formula heuristically.

## Context

SOC breaks global spin-rotation invariance, modifies the sum rule/Goldstone logic, and can couple transverse spin response to charge/other spin components. The code already has noncollinear/SOC ground-state capabilities, which makes accidental unsupported use especially easy.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.

## Work plan

1. Audit exactly what happens today if TD-DFT is requested with SOC/noncollinear ground states and ensure it fails safely.
2. Derive the full spinor response matrix required in the code’s Hamiltonian convention.
3. Summarize the relativistic sum-rule structure of dos Santos Dias et al. 2015 and map required ground-state quantities.
4. Identify anisotropy/external-field terms that replace the zero-energy Goldstone expectation.
5. Define how charge coupling and SOC-induced damping enter.
6. Specify validation fixtures and independent static anisotropy checks.
7. Produce an implementation sequence, but do not enable the feature.

## Acceptance checklist

- [x] Current unsupported cases fail explicitly.
- [x] Required full spinor response structure is documented.
- [x] Relativistic sum-rule implications are mapped.
- [x] SOC gap/anisotropy validation plan exists.
- [x] No collinear-only kernel is silently reused.
- [x] Feature remains guarded.

## Required evidence / deliverables

Create `docs/TDDFT_RELATIVISTIC_BLUEPRINT.md`.

Completion evidence is recorded in
[`docs/TDDFT_RELATIVISTIC_BLUEPRINT.md`](../../../TDDFT_RELATIVISTIC_BLUEPRINT.md).
The production entry-point guard and its ordering are covered by
[`tests/unit/test_tddft_dispatch.py`](../../../../tests/unit/test_tddft_dispatch.py).

## One-line commit message

`docs: define relativistic TDDFT extension`
