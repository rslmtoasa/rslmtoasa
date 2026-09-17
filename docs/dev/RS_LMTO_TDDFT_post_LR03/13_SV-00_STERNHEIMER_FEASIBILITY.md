# SV-00 — Feasibility audit for an independent LMTO Sternheimer response route

## Timing

After the collinear radial Dyson implementation has been validated.

Do not use Sternheimer to patch unresolved defects in the primary route.

## Goal

Assess Savrasov-style LMTO linear response as a genuinely independent
implementation/benchmark route.

## Questions

1. What first-order KS/Sternheimer equations are solved?
2. What occupied/unoccupied projection is required?
3. What first-order density variables are represented inside spheres/interstitial?
4. How does the self-consistent induced XC field enter?
5. What response of the LMTO basis itself is required?
6. Are derivatives of:
   - partial waves;
   - potential parameters;
   - screening/structure constants;
   - overlap/orthogonalization
   required?
7. Can those derivatives be expressed consistently in the current second-order
   orthogonal RS-LMTO representation?
8. Does the Pauli/scalar-relativistic boundary identified by LR-02R appear in the
   Sternheimer formulation as well?
9. Can the route produce the same LR-04 response-space observables for direct
   comparison?
10. Which equations are fully specified in the literature and which would require
    unpublished implementation knowledge?

## Literature discipline

Separate:

- equations stated by Savrasov;
- standard Sternheimer identities;
- LMTO-specific derivations made here;
- implementation assumptions.

If a required LMTO basis-response term cannot be derived from published formalism
and live code, declare BLOCKED.

## Benchmark role

If feasible, define an eventual benchmark plan:

- static q=0 response;
- selected finite q;
- comparison with LR-06 spectral response;
- no shared band-sum denominator implementation.

Do not implement the benchmark here.

## Deliverable

`docs/LMTO_STERNHEIMER_RESPONSE_FEASIBILITY.md`

## Commit

`docs: assess LMTO Sternheimer spin-response route`
