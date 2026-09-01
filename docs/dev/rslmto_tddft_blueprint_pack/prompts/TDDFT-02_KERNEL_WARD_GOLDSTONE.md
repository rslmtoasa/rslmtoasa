# TDDFT-02 — ground-state-consistent kernel, Ward identity, and Goldstone diagnostics

## Agent mission

Repair/verify the transverse interaction kernel so the static Ward/Goldstone identity is satisfied by the same ground-state Hamiltonian representation, without empirical tuning.

## Context

Halle and Lounis both show that the Goldstone problem is fundamentally a consistency relation between chi_KS, the ground-state exchange field, magnetization and Kxc. Lounis gives a real-space sum-rule construction; Halle writes the equivalent null-vector condition for D=I-chi_KS Kxc. Our default must be a physically derived kernel, not an arbitrary lambda rescaling.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.
- Do not make `goldstone_policy=sum_rule/projected` the hidden default. First expose the uncorrected residual.
- If a site-scalar kernel cannot converge, investigate response-basis resolution before tuning the scalar.

## Work plan

1. Derive how the ground-state XC magnetic field/exchange splitting enters the modern RS-LMTO Hamiltonian and response basis.
2. Construct the static transverse Kxc/operator in exactly that basis.
3. Evaluate the discrete Ward identity `chi0(0,0) Bxc -> m` and `D m -> 0`. Include all projection/volume/LMTO normalization factors explicitly.
4. Make the uncorrected residual a first-class diagnostic printed and written to provenance.
5. Converge the residual on small Fe/Ni fixtures with k mesh/energy settings.
6. If the residual plateaus, add a diagnostic `site_orbital` response basis or equivalent richer local basis and determine whether the scalar site projection is the cause.
7. Implement optional, explicit numerical cleanup modes only after the physical kernel is in place:
   - Lounis-style sum-rule reconstruction in the active response basis;
   - Halle-style projection of the single small spurious Goldstone eigenvalue.
8. Record the magnitude of any correction and make a large correction a warning/failure according to evidence-based thresholds.
9. Forbid empirical global lambda scaling.

## Acceptance checklist

- [x] Static Bxc, m and Kxc come from the same ground-state representation.
- [x] Ward residual is computed before correction.
- [ ] Residual decreases under numerical refinement.
- [ ] Site-only versus richer-basis behavior is understood.
- [x] Optional sum-rule/projective repairs are explicit and provenance-recorded.
- [x] No empirical lambda/energy shift exists in the default path.
- [x] SOC=0 Goldstone test is present.

## Required evidence / deliverables

Add `docs/TDDFT_GOLDSTONE_WARD.md` with derivation, numerical residual tables and a decision on the production response basis. Include tests for both the raw identity and optional repair mechanisms.

## One-line commit message

`tddft: make transverse kernel Ward consistent`
