# LR-06 — Implement collinear transverse Kohn-Sham susceptibility

## Prerequisites

- LR-03 PASS.
- LR-04 PASS.
- LR-05 PASS.

## Nature

FIRST BARE-RESPONSE IMPLEMENTATION.

No XC kernel.
No Dyson enhancement.
No Goldstone correction.
No mode fitting.

## Goal

Implement the retarded collinear transverse Kohn-Sham susceptibility in the full
approved radial/angular response space.

The first backend is the explicit spectral/Lehmann band sum because its one-electron
ingredients are the best-certified electronic-structure path.

## Physics contract

Use **only** the equations and channel ordering fixed by
`docs/LR_TDDFT_CONVENTIONS.md`.

Do not copy signs, factors, channel names, broadening conventions, or endpoint
ordering from purged TD-DFT code.

The generic form is

\[
\chi^{KS,\mu\nu}_{IJ}(\mathbf q,\omega)
=
\frac{1}{N_k}
\sum_{\mathbf k,nm}
\frac{
f_{n\mathbf k}-f_{m,\mathbf k+\mathbf q}
}{
\omega+\epsilon_{n\mathbf k}-\epsilon_{m,\mathbf k+\mathbf q}+i\eta
}
T^\mu_{nm;I}
(T^\nu_{nm;J})^*,
\]

with the exact production specialization determined by LR-03.

## Response-space representation

Use LR-04 exactly.

Do not decide locally whether `chiKS` stores:

- raw kernel values;
- weighted operator values;
- symmetrically weighted values.

Follow the canonical representation and its adjoint/composition rules.

No duplicated quadrature logic.

## Electronic-state provenance

All k and k+q eigensystems must derive from one immutable electronic-state contract:

- accepted converged potential;
- Fermi level;
- electronic temperature/smearing;
- occupations;
- reciprocal mode;
- band count;
- k-mesh;
- energy zero.

Do not silently recompute `EF`.

Do not mix fixed and automatically determined occupations.

## Band completeness

Initial implementation uses the full finite LMTO eigensystem.

Do not add an energy/band truncation as a default optimization.

If a band-window option already exists at a lower level, keep it disabled for
certification unless an explicit convergence test proves its error.

## Frequency and broadening

Use exactly the retarded convention from LR-03.

Make `eta` explicit in the response request.

Do not use a hidden broadening inherited from DOS/GF settings unless LR-03
explicitly made them identical.

## q handling

Use arbitrary-k eigensystems and the LR-05 production transition vertex.

Do not special-case BZ crossings with hand-written phases.

Support:

- q=0;
- off-mesh q;
- k+q folding.

## Initial channel scope

Implement only the transverse collinear sector(s) certified by LR-03.

Do not implement a nominal four-component response merely by allocating four
channels.

## Static limit

Provide a numerically explicit `omega=0` path using the same spectral formula.

Do not create a separate approximate static susceptibility.

This static object will later feed the Lounis sum rule and Goldstone diagnostics.

## Tests — independent levels

### Algebraic finite model

Construct a tiny two-level or few-level collinear model with independently known
transition amplitudes.

Evaluate the susceptibility directly from the analytic finite sum and compare with
production accumulation.

The reference must not call production transition/susceptibility routines.

### Direct explicit-loop oracle

For a small deterministic production-style fixture, compare optimized accumulation
against a deliberately straightforward nested-loop implementation.

### Covariance

Test all identities actually derived in LR-03, including as applicable:

- q ↔ -q;
- retarded/advanced conjugation;
- circular-sector exchange;
- global spin reversal.

Do not add an identity not guaranteed by the certified assumptions.

### Degenerate-subspace invariance

Rotate a degenerate eigen-subspace and verify the summed susceptibility is
unchanged.

### Convergence diagnostics

For one real Fe fixture report:

- k mesh;
- eta;
- radial representation;
- angular cutoff.

At this task these are diagnostics, not final material validation.

## Raw Ward / Goldstone evidence

At q=0 and omega=0, compute the **left-hand ingredients** needed by the LR-03
static identity using the independently stored ground-state field/magnetization.

Report the raw residual.

Do not alter `chiKS`.

Do not rescale a kernel.

Do not shift an eigenvalue.

This is diagnostic evidence only.

## API

Prefer a clean request/result structure containing:

- q;
- frequency/frequencies;
- eta;
- channel;
- electronic-state reference;
- response-space metadata.

Do not couple the core evaluator to plotting or magnon peak extraction.

## Memory/performance

Correctness is the first target, but avoid obviously pathological storage.

It is acceptable to compute one q/frequency/channel block at a time.

Do not compress radial/angular space to sites for performance.

## Deliverable

Create:

`docs/LR_KS_SUSCEPTIBILITY.md`

Include:

- equations;
- provenance;
- response representation;
- channel convention;
- tests;
- raw static residual;
- limitations.

## Forbidden

- no Kxc;
- no Dyson;
- no Goldstone enforcement;
- no site projection;
- no hidden band truncation;
- no sign/factor tuning from Fe behavior.

## Checklist

- [x] LR-03 conventions used verbatim;
- [x] LR-04 response algebra used;
- [x] LR-05 vertices used;
- [x] immutable occupations/EF provenance;
- [x] full finite spectrum baseline;
- [x] q=0 and arbitrary q supported;
- [x] analytic finite-model oracle passes;
- [x] explicit-loop oracle passes;
- [x] covariance tests pass;
- [x] degenerate-subspace invariance passes;
- [x] raw q=0 static residual reported;
- [x] no interaction physics added.

## Implementation outcome

`source/lr_ks_susceptibility.f90` provides the immutable electronic-state
snapshot contract, exact arbitrary-`k+q` endpoint adapter, spectral/Lehmann
accumulator, canonical LR-04 output, and raw static-residual diagnostic.
`UnitLrKsSusceptibility` passes its analytic, explicit-loop, covariance, folded
endpoint, and degenerate-subspace tests.  The implementation claim is bare
spectral KS susceptibility only; interacting response and material validation
remain downstream tasks.

## Commit

`linear-response: implement radial transverse KS susceptibility`
