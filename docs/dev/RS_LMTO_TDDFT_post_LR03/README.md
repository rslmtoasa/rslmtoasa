# RS-LMTO-ASA TD-DFT clean-room campaign — post-LR-03 prompt pack

This pack supersedes the downstream portion of the original clean-room prompt pack.

It assumes the campaign has already completed the foundation work:

- LR-00-PURGE
- LR-GF-01
- LR-BASIS-00
- LR-01
- LR-02 / LR-02R
- the SR→Pauli numerical closure required by LR-02R
- LR-03 transverse-response conventions and normalization

The exact repository state remains the live `fable_v4` branch. Every prompt starts
with a preflight requirement to inspect the live evidence rather than assuming that
the expected PASS conclusion was obtained.

## Why the old downstream pack is replaced

The old tasks no longer match the code/history:

- old LR-04 ("expose radial ground-state data") was subsumed by LR-01;
- much of old LR-05 ("response basis") landed during LR-02R;
- the direct radial mesh introduces a nontrivial discrete integration metric that
  must be locked before susceptibility, kernels, or Dyson algebra are implemented;
- the production Pauli transition vector is still distinct from the formal mapping;
- `chiKS` should first be implemented through the certified eigenpair/Lehmann
  spectral route;
- a reciprocal GF bubble is retained as an independent cross-validation route, not
  silently conflated with the spectral implementation;
- direct ALSDA, the Lounis sum-rule interaction, and the BES eigenvalue correction
  remain separate published routes;
- the native RS-GF response backend remains deferred until the LMTO
  representation-transformation contract is independently certified.

## Revised sequence

1. LR-04 — Discrete response-space metric and operator algebra
2. LR-05 — Production Pauli transition-vector evaluator
3. LR-06 — Collinear transverse KS susceptibility (spectral/Lehmann route)
4. LR-GF-02 — Reciprocal-GF susceptibility cross-check
5. KXC-01 — Local radial ALSDA transverse kernel
6. GSR-01 — Independent Goldstone sum-rule interaction
7. GCR-01 — Optional published Goldstone eigenvalue correction
8. TDDY-01 — Enhanced susceptibility and loss matrix
9. TDVAL-01 — bcc Fe / fcc Ni collinear validation
10. LR-REP-00 — Native RS-GF representation transformation audit
11. RSGF-01 — Native real-space GF response backend
12. NC-00 — Noncollinear four-component feasibility audit
13. SV-00 — LMTO Sternheimer feasibility audit

`LR-GF-02` is strongly recommended before final collinear validation, but need not
block KXC-01 if the LR-06 spectral implementation has passed all of its independent
tests. It **must** be completed before claiming two independent reciprocal
susceptibility backends.

The initial production capability remains intentionally narrow:

- collinear;
- no SOC;
- no Hubbard/additive response operator;
- orthogonal `ham_only`;
- second-order/HOH LMTO;
- `sp`/`spd`;
- Pauli/no-SOC response evaluated from the documented scalar-relativistic ground
  state;
- direct radial response mesh;
- no site-only fallback.

## Campaign rule

For every task distinguish explicitly:

- algebraic consistency;
- independent numerical cross-check;
- converged-material validation;
- literature agreement.

Never promote one category into another.
