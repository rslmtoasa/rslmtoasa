# TDDFT-04 — common chi0 backend interface and response data model

## Agent mission

Refactor only as much as necessary so eigenpairs, G(k,z), and G(R,z) can feed the same physics/Dyson layer without duplicated conventions or output logic.

## Context

The target is one response framework with multiple propagator backends. The interface must allow batching, especially because the R-space backend can build chi0(R,w) once and Fourier transform it to many q points.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.

## Work plan

1. Identify the smallest common backend contract from the current code rather than imposing a new hierarchy.
2. Separate backend-specific one-electron work from common response basis, kernel, Dyson and spectral analysis.
3. Ensure the contract can represent batched omega/q evaluation and backend metadata.
4. Add a backend factory/selector preserving the historical `chi0_backend='eigenpairs'` behavior.
5. Define explicit backend names for K-space Lehmann GF and native R-space GF; add compatibility aliases only if needed.
6. Do not make G(R)->G(k) part of the common contract.
7. Add a mock/test backend to verify that the common Dyson layer is backend-agnostic.
8. Keep memory ownership and constructor/destructor behavior consistent with existing Fortran types.

## Acceptance checklist

- [x] Eigenpair behavior unchanged through the new interface.
- [x] Backend interface supports batch-oriented execution.
- [x] Backend provenance is returned to the common layer.
- [x] No forced G(R)->G(k) transform exists.
- [x] Dyson/spectral code no longer contains backend-specific branches except where mathematically necessary.
- [x] Tests cover backend selection and unsupported combinations.

Evidence: `tddft_backend_mod` provides the common batch interface and adapters;
the eigenpair adapter delegates to the existing chi0 builder; metadata carries
canonical backend and q/omega provenance; the real-space adapter requires a
native provider and performs no G(R)->G(k) transform. `UnitTddftBackend`,
`UnitTddftConfig`, the TDDFT-labelled CTest suite (20/20), and
`tests/unit/test_tddft_dispatch.py` pass.

## Required evidence / deliverables

Document the backend contract in `docs/TDDFT_BACKENDS.md` with a short sequence diagram.

## One-line commit message

`tddft: unify chi0 backend interface`
