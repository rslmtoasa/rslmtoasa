# TDDFT-13 — coupled charge-longitudinal response

## Agent mission

Implement the `(0,z)` response block as a physically separate extension after the transverse release gate, using Buczek et al. 2020 as the principal physics target.

## Context

In collinear SOC-free systems the transverse block decouples, while charge and longitudinal magnetization remain coupled. The Hartree/Coulomb response is essential here; a longitudinal mode is not a Goldstone mode and a static susceptibility is not automatically an LLB damping coefficient.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.

## Work plan

1. Extend the response basis/API to the coupled charge and z-magnetization sector without duplicating the common backend machinery.
2. Implement the KS charge/longitudinal response with all three electronic backends where feasible, or document justified staging.
3. Add Hartree/Coulomb screening in the charge channel with correct units/metric.
4. Implement the longitudinal xc-kernel components consistent with the ground-state XC functional.
5. Validate charge-longitudinal coupling in KS response and its effective decoupling after plasmon screening in representative 3D metals.
6. Reproduce qualitative Fe/Ni and AF-FeSe “magnetic second sound” behavior from Buczek et al. 2020.
7. Add frequency-dependent low-energy output needed for a later LLB-relaxation study, but do not label it a damping parameter yet.
8. Preserve strict guards for SOC/noncollinear cases.

## Acceptance checklist

- [x] `(0,z)` matrix structure implemented in the common site-major response basis.
- [x] Hartree screening included in the charge-charge block using the projected Rydberg metric.
- [x] Longitudinal xc response derives from the ground-state `XCPOT_hybrid` route through a local ALSDA Hessian.
- [x] No false Goldstone constraint is imposed longitudinally.
- [ ] Fe/Ni/FeSe qualitative validation completed; this remains a material-run gate beyond the analytic/unit fixture.
- [x] Output is suitable for later low-frequency dissipative analysis and explicitly carries numerical-eta/LLB disclaimers.
- [x] Transverse regression suite remains unchanged.

## Implementation evidence

`docs/TDDFT_LONGITUDINAL.md` records the response convention, kernel
provenance, backend staging, strict guards, and interpretation boundary.
`UnitTddftLongitudinal` covers the coupled channel layout, charge-only
Hartree insertion, XC capability gate, and shared eigenpair adapter.  The
native real-space backend remains explicitly rejected for this sector pending
its multi-component Fourier/static-limit derivation.

## Required evidence / deliverables

Write `docs/TDDFT_LONGITUDINAL.md` including the distinction between longitudinal susceptibility, collective modes and any future LLB parameter mapping.

## One-line commit message

`tddft: add coupled longitudinal charge-spin response`
