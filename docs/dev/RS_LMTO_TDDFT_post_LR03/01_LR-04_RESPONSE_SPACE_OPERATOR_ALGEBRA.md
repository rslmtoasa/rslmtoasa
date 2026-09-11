# LR-04 — Lock the discrete response-space metric and operator algebra

## Prerequisites

- LR-03 PASS.
- The current LR-02R response basis mapping and LR-01 radial-state evidence are
  present.
- If LR-03 declares any unresolved convention needed here, STOP.

## Nature

FOUNDATIONAL IMPLEMENTATION.

No transition-vector evaluator.
No `chiKS`.
No XC kernel.
No Dyson physics.

## Goal

Turn the already-derived direct radial/angular super-index into one unambiguous
finite-dimensional operator space.

The central question is:

\[
\boxed{\text{Where does the radial quadrature metric enter every response
contraction?}}
\]

Later susceptibility, kernel and Dyson code must not answer that question
independently.

## Existing basis

Use the current canonical response coordinate

\[
I=(a,L,M,i,\mu)
\]

or the exact live equivalent established by LR-02R.

Reuse current `response_angular_basis` / `response_basis_mapping` utilities where
appropriate. Do not create parallel indexing code.

## Derivation

Starting from the physical continuum contraction

\[
y_a^{LM\mu}(r)
=
\sum_{bL'M'\nu}
\int dr'\,r'^2
A_{ab}^{LM\mu,L'M'\nu}(r,r')x_b^{L'M'\nu}(r'),
\]

derive the exact discrete form on the production logarithmic mesh.

Use the LR-01/LR-02R quadrature. Explicitly identify the radial measure

\[
W_i
\]

including the logarithmic Jacobian, Simpson coefficient and any required
\(r_i^2\).

Do not reuse the weighted legacy `RHO` measure where the response variable is a
physical pointwise density.

## Choose one canonical stored representation

Evaluate at least these possibilities:

1. raw pointwise kernel `A_ij` with explicit metric `W` in contractions;
2. right-weighted operator `A W`;
3. symmetrically weighted representation
   \[
   W^{1/2} A W^{1/2}.
   \]

Choose exactly one as the **canonical production matrix representation**.

The choice must be derived, not made for coding convenience.

### Origin handling

The production radial mesh includes the origin. If a symmetric
\(W^{1/2}\) representation is singular or ambiguous because the physical volume
weight at the origin vanishes, do not hide the problem with an epsilon.

Either:

- keep the explicit metric representation; or
- derive a mathematically clean active-grid transformation.

No arbitrary regularization.

## Required algebra

Define and test:

- response-vector inner product / norm;
- kernel action on a field;
- composition of two nonlocal response operators;
- application of a pointwise local radial operator;
- identity operator;
- adjoint/Hermitian conjugation in the chosen discrete metric;
- trace if later diagnostics use one;
- rigid-rotation vector norm/overlap.

If the canonical representation transforms the metric away, give exact
raw↔canonical transformations.

## Local operator

For a local radial scalar \(K_a(r)\), derive its action from

\[
(Kx)_a(r)=K_a(r)x_a(r).
\]

Do not represent a Dirac delta by guessing a `1/W` factor.

Derive the stored matrix representation from the chosen convention and verify that
the numerical action is pointwise multiplication.

This will later be the contract for the ALSDA kernel.

## Angular structure

For a spherical ground-state local scalar operator, demonstrate which indices are
diagonal:

- site;
- \(L,M\);
- radial coordinate;
- channel, where applicable.

Do not introduce an additional angular quadrature.

The spherical harmonics have already been projected analytically.

## Tests

Add independent small tests for:

1. analytic radial functions on the logarithmic mesh;
2. identity action;
3. pointwise local multiplication;
4. composition of two integral kernels compared with explicit quadrature loops;
5. adjoint relation in the metric;
6. super-index multisite separation;
7. origin behavior;
8. raw↔canonical round trip if a transformed representation is chosen.

At least one test must compare the production helper against a deliberately simple
explicit nested-loop quadrature oracle.

Do not have both paths call the same matrix-multiplication helper.

## Suggested code

Prefer one small module, e.g.

`source/lr_response_space.f90`

or extend the existing response-basis mapping module if that is cleaner.

Avoid a new class hierarchy.

## Documentation

Create:

`docs/LR_RESPONSE_SPACE_ALGEBRA.md`

It must contain:

- exact stored representation;
- response metric;
- operator action;
- composition rule;
- adjoint rule;
- local multiplication rule;
- units inherited from LR-03;
- origin treatment;
- test evidence;
- unsupported representations.

## Forbidden

- no `chiKS`;
- no `T_nm`;
- no `Kxc` physics;
- no Dyson equation;
- no Goldstone correction;
- no site reduction;
- no hidden quadrature factors inside unrelated modules.

## Acceptance checklist

- [x] exact production radial metric derived;
- [x] canonical stored representation chosen and justified;
- [x] origin handled without ad-hoc regularization;
- [x] identity action passes;
- [x] local multiplication passes;
- [x] operator composition passes independent quadrature oracle;
- [x] adjoint convention fixed;
- [x] rigid-vector overlap/norm defined;
- [x] no response physics implemented;
- [x] documentation complete.

## Commit

`linear-response: define response-space operator algebra`
