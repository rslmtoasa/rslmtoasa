# TDDFT-07 — implement native real-space GF chi0 and susceptibility Fourier transform

## Agent mission

Implement the primary new RS-LMTO TD-DFT backend using native `G(R,z)` directly: build `chi0(R,omega)` from real-space propagators and Fourier-transform the susceptibility, not the one-electron GF, for periodic q-resolved output.

## Context

This is the architecture chosen in the discussion. It exploits the code’s mature RS GF and is the natural route for bulk, films, surfaces, impurities and embedded geometries. It must not depend on the newer G(k) implementation.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.
- Do not introduce G(R)->G(k) as a production intermediate. A Fourier G comparison may exist only as a validation test.
- Preserve site/orbital directionality: use the correct pair `G_ab(R)` and reverse propagator `G_ba(-R)` with retarded/advanced structure.

## Work plan

1. Map the native RS GF’s indices, energy convention, spin blocks and off-diagonal site/orbital access.
2. Implement a direct near-real-axis `chi0_ab(R,w)` contraction with the same source/measurement vertices as other backends.
3. Add real-space symmetry/conjugation tests without assuming more symmetry than collinear SOC=0 provides.
4. Implement lattice FT `chi0(R,w)->chi0(q,w)` for bulk. Keep a clean abstraction that can later do only an in-plane FT for films.
5. Add explicit R-shell/Rmax control and tail diagnostics based on chi0, not only G.
6. On a periodic small fixture, compare:
   - native G(R,z) Fourier transformed to G(k,z) against validated K-space G (validation only);
   - r-GF chi0(q,w) against k-GF/eigenpair chi0.
7. Exploit q amortization: at fixed w, build chi0(R,w) once and transform to all q in the requested path.
8. Keep a direct real-space `chi_ij(w)` output path for finite/local geometries where q is meaningless.

## Acceptance checklist

- [x] Native R-GF chi0 works without G(k).
- [x] R-space direction/reversal identities are tested.
- [x] Bulk FT of susceptibility is implemented and tested.
- [x] Rmax convergence/tail metrics are exposed.
- [x] Periodic r-GF chi0 agrees with k-space backends after convergence.
- [x] Multiple q values reuse the same chi0(R,w) build.
- [x] Finite/local response representation remains possible.

## Required evidence / deliverables

Add `docs/TDDFT_REALSPACE_GF.md` with the exact contraction, FT convention, R cutoff diagnostics and bulk/film/finite representation plan.

## One-line commit message

`tddft: add native real-space Green-function chi0`
