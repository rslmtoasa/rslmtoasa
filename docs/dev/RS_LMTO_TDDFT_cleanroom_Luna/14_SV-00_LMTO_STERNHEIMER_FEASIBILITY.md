# SV-00 — Feasibility audit for an independent LMTO Sternheimer spin-response route

## Timing
Deferred until collinear radial Dyson TDDFT is validated.

## Goal
Assess Savrasov 1998 as an independent LMTO-native linear-response route.

Do not use it to fill missing equations in earlier tasks.

## Published target
Savrasov, Phys. Rev. Lett. 81, 2570 (1998):
- variational TD linear response;
- Sternheimer formulation;
- muffin-tin-orbital representation;
- dynamical spin susceptibility.

## Questions
1. Does the published material contain sufficient implementation detail?
2. What first-order wavefunction/potential equations are required?
3. What self-consistent induced XC field is required?
4. What LMTO basis-response terms arise?
5. Are they compatible with the current orthogonalized RS-LMTO representation?
6. Can this produce an independent benchmark against the radial Dyson route?

## Rule
If unpublished details or guessed LMTO basis derivatives are required: BLOCKED.

## Deliverable
`docs/LMTO_STERNHEIMER_RESPONSE_FEASIBILITY.md`

## Checklist
- [ ] equations inventoried
- [ ] first-order variables identified
- [ ] basis derivatives assessed
- [ ] self-consistent field response assessed
- [ ] benchmark role identified
- [ ] PASS or BLOCKED

## Commit if documentation is added
`docs: assess LMTO Sternheimer spin-response route`
