# TDDFT-01 — lock conventions, response definitions, and exact identities

## Agent mission

Turn all sign/factor/unit conventions in the transverse response into an explicit mathematical contract and executable tests before further implementation.

## Context

The current Ni behavior and earlier discussion make factor-of-two, circular-channel, moment-sign and retarded-frequency errors plausible failure modes. Halle/Jülich formulas use different notational conventions; the code must derive its convention from the RS-LMTO Hamiltonian, not copy symbols blindly.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.

## Work plan

1. Write `docs/TDDFT_CONVENTIONS.md` defining Pauli matrices, sigma+/-, external magnetic perturbation, magnetization units/sign, Fourier convention, retarded GF, chi^{+-}/chi^{-+}, loss matrix and spectral sign.
2. Trace each convention to the actual Hamiltonian and moment implementation.
3. Add tiny analytic one-/two-level fixtures with known transverse chi0.
4. Test positive/negative frequency relations, retarded/advanced conjugation and channel reversal.
5. Add a finite-difference Hamiltonian rotation test: rotate the local spin-dependent potential by a tiny transverse angle and compare dH/dtheta with the analytic source vertex used by response.
6. State explicitly which object is the source/measurement vertex and which is Kxc. Prevent double counting.
7. Add a convention-dependent spectral-weight/sum-rule test on the analytic fixture.

## Acceptance checklist

- [x] Exact code convention documented.
- [x] sigma+/sigma- factors are unit-tested.
- [x] Retarded sign convention is unit-tested.
- [x] Finite-difference transverse dH test passes.
- [x] Source vertex and interaction kernel are separate objects in documentation/code.
- [x] Analytic chi0 fixture matches code.
- [ ] Existing TD-DFT tests still pass.

Verification note: the 17 non-profile TDDFT-labelled tests pass. The final item remains unchecked because `UnitTddftCpuProfile` still fails in the pre-existing reciprocal-backend fixture at `source/reciprocal_fourier.f90:1059`.

## Required evidence / deliverables

Provide the derivation in documentation and a compact test file that would fail under a sign or factor-of-two regression.

## One-line commit message

`tddft: lock response conventions and identities`
