# RSGF-01 — Implement the native real-space GF bare-response backend

## Prerequisites

- TDVAL-01 reciprocal baseline validated.
- LR-REP-00 PASS for the required native GF representation.
- LR-04/LR-05 response-space and vertex contracts available.

## Goal

Implement a genuinely independent real-space Green-function route to the same
Pauli radial/angular `chiKS` object.

This is not a replacement for the reciprocal implementation.

It is an additional backend and cross-check.

## Supported routes

Assess separately:

- block-recursion native GF;
- Chebyshev native GF.

Implement only a route whose representation/augmentation was certified by
LR-REP-00.

## Physics

Derive the real-space GF susceptibility from the same LR-03 Kubo convention.

Map both GF endpoints through the certified physical/Pauli augmentation.

The final result must occupy exactly the same LR-04 response space as LR-06.

No site projection.

## Fourier transform

For periodic comparisons, derive:

\[
\chi(\mathbf q,\omega)
=
\sum_R e^{\pm i\mathbf q\cdot R}\chi(R,\omega)
\]

with sign and endpoint positions consistent with the already-certified reciprocal
gauge.

Do not guess the phase from historical response code.

## Numerical energy integration

Document:

- integration contour/real axis;
- eta;
- Fermi function;
- energy window;
- truncation/convergence.

The native GF backend must not call LR-06 under the hood.

## Tests

- finite cluster exact inverse oracle;
- RS GF response versus explicit finite spectral sum;
- periodic small-cell RS versus reciprocal LR-06;
- q phase;
- radial/angular full-space comparison;
- block versus Chebyshev comparison where both are certified;
- energy-integration convergence.

## Claims

Agreement with LR-06 validates the independent GF response construction for the
shared finite problem.

It does not independently validate common ground-state augmentation/XC physics.

## Deliverable

`docs/TDDFT_RS_GF_BACKEND.md`

## Commit

`linear-response: add native real-space GF susceptibility backend`
