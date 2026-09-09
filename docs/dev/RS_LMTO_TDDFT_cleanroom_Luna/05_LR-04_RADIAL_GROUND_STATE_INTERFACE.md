# LR-04 — Expose converged radial ground-state quantities

## Prerequisites
LR-01, LR-02, LR-03 PASS.

## Goal
Provide read-only access to the already-proved converged radial ground-state data needed by linear response. No TD-DFT equation is introduced here.

## Suggested file
`source/lr_radial_data.f90`

Keep the flat `source/` layout.

## Required data contract
Expose, per site/type as derived:
- radial mesh;
- integration weights/Jacobian;
- n(r);
- m(r), or spin-resolved densities from which it is exactly constructed;
- Bxc(r), or spin-resolved Vxc from which it is exactly constructed;
- sphere radius;
- XC functional provenance;
- units.

The interface must make it impossible to confuse raw radial density with integrated density or moment.

## Lifetime
The data must correspond to the exact converged ground state used by the response.

If legacy routines overwrite arrays, add the narrowest snapshot/copy at the point where the exact converged radial quantities exist.

Do not reconstruct Bxc from LMTO potential parameters.

## Tests
1. radial integral reproduces charge;
2. radial integral reproduces magnetic moment;
3. Bxc matches the exact ground-state XC spin splitting point-by-point;
4. restart retains XC provenance;
5. restart and fresh converged calculations expose the same radial state within tolerance.

## Forbidden
- no Kxc;
- no site averaging;
- no unapproved interpolation;
- no changes to legacy SCF numerics.

## Checklist
- [ ] read-only radial interface
- [ ] exact converged-state provenance
- [ ] charge integral passes
- [ ] moment integral passes
- [ ] Bxc pointwise provenance passes
- [ ] restart path passes
- [ ] no potential-parameter surrogate
- [ ] SCF results unchanged

## Commit
`linear-response: expose converged radial ground-state data`
