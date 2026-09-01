# TDDFT-05 — validate the newer K-space Lehmann Green function

## Agent mission

Before using the new K-space GF in TD-DFT, establish that it is a correct one-electron resolvent over the exact electronic representation supported by TD-DFT.

## Context

The K-space GF is newer and not yet production-experienced. A TD-DFT backend built on an unvalidated propagator would make debugging ambiguous.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.

## Work plan

1. For representative k points and complex z, compare the Lehmann G(k,z) against direct inversion of `z-H(k)` in the supported orthogonal representation.
2. Verify spectral normalization, Hermiticity/retarded-advanced identities and large-|z| asymptotics.
3. Check spin blocks and nonmagnetic/collinear limits.
4. Compare the local DOS/spectral function recovered from G(k,z) with the existing eigenvalue/DOS path.
5. Explicitly test and preserve the guard for generalized overlap. Add a design note for `G=(zS-H)^-1` rather than approximating it.
6. Benchmark cost and allocation patterns for batched complex energies.
7. Add deterministic regression fixtures that do not depend on broad visual comparison.

## Acceptance checklist

- [x] Lehmann and direct-inverse G(k,z) agree.
- [x] Retarded/advanced identities pass.
- [x] Spectral normalization/asymptotics pass.
- [x] Spin convention matches TD-DFT convention.
- [x] Generalized-overlap case is rejected with a precise explanation.
- [x] Basic performance profile exists.

Evidence: `docs/KSPACE_GF_VALIDATION.md` records the deterministic fixture,
direct-inverse maximum error (`1.0549e-14`), retarded/advanced and spectral
identity results, spin/DOS checks, the orthogonal-only metric boundary, and a
batched-energy allocation/timing profile.  The executable test is
`UnitKspaceGFValidation`; `UnitKspaceGFValidationSource` checks that the
validated primitive is shared by the response provider and that both the
reciprocal Green and TD-DFT entry points retain the generalized-overlap guard.

## Required evidence / deliverables

Add a focused `docs/KSPACE_GF_VALIDATION.md` and tests usable independently of TD-DFT.

## One-line commit message

`greens: validate Lehmann k-space resolvent`
