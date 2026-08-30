# Luna Prompt — GBT-WP03: Composite LMTO Covariance (HOH, Sdot/CCOR, Overlap)

## Mission

Starting from a GBT core that has passed q=0 and exact shifted-k tests, establish which composite LMTO terms preserve the same covariance and may safely be enabled in production GBT.

Do not test all corrections simultaneously. Enable one feature at a time so a failure is localizable.

## Dependency

WP02 must pass with the minimal operator.

## Features to audit sequentially

1. HOH / second-order LMTO term.
2. `Sdot` and two-center CCOR.
3. generalized-overlap / nonorthogonal path if active in the relevant solver.
4. local Hubbard U only if its present implementation is meaningfully compatible with the rotating-frame formulation.

Continue to reject during this work package:

- SOC;
- intersite Hubbard V;
- unaudited local-axis functionality.

## Mathematical requirement

If a primitive operator obeys endpoint covariance

\[
h'_{ij}=U_i^\dagger h_{ij}U_j,
\]

then a composite product with a local onsite operator such as

\[
\sum_j h_{ij}o_jh_{jk}
\]

must satisfy

\[
\sum_j h'_{ij}o'_jh'_{jk}
=
U_i^\dagger
\left(\sum_jh_{ij}o_jh_{jk}\right)
U_k
\]

provided `o_j` transforms locally and no extra/double gauge is inserted.

Use this to derive expected HOH behavior in the code.

For CCOR, establish exactly whether raw `Sdot` has the same support/directed-bond ordering as `S` and verify that the same endpoint link is applied once at the primitive nonlocal level.

## Required test sequence for each feature

For each newly enabled term repeat:

1. q=0 ordinary/GBT exact operator comparison;
2. real-space bond reversal;
3. reciprocal Hermiticity;
4. one-sublattice `k +/- q/2` identity where mathematically applicable;
5. representative fixed-potential eigenvalue comparison.

Do not relax tolerances merely because the operator is more complicated. If the same arithmetic path differs only by a unitary/gauge transformation, the residual should still be roundoff-scale.

## Audit questions

### HOH

- Is HOH formed from already-linked `ee`/hopping objects?
- Is any additional GBT phase applied later?
- Are onsite overlap/linearization objects truly local in the rotating frame?
- Does basis ordering remain consistent?

### CCOR / Sdot

- Is `Sdot` taken from the same directed neighbor slot as `S`?
- Does its physical displacement exactly match?
- Is the endpoint link applied once and only once?
- Are any CCOR VMT terms local or pair-dependent, and if pair-dependent do they require covariance treatment?

### Generalized overlap

- Is the solver using `H c = E O c`?
- Does `O` itself contain nonlocal terms requiring the same gauge?
- Is the ordinary shifted-k identity formulated in the generalized eigenproblem rather than only H?

### Hubbard U

- Is the present U potential site-local and expressed in a frame that can be consistently rotated?
- If not proven, leave it rejected rather than forcing support.

## Deliverables

Create `docs/dev/GBT_REAUDIT_WP03_COMPOSITE_COVARIANCE.md` with a per-feature verdict:

- supported and exactly tested;
- conceptually supported but untested;
- unsupported and explicitly guarded.

Update feature guards only where evidence justifies it.

## Completion checklist

- [ ] HOH covariance derived from actual code path.
- [ ] HOH exact tests pass or feature remains guarded.
- [ ] Sdot directed support/order verified.
- [ ] CCOR exact tests pass or feature remains guarded.
- [ ] Generalized-overlap covariance audited if relevant.
- [ ] Local Hubbard U status explicitly decided from evidence.
- [ ] No unsupported feature silently enabled.
- [ ] No double-gauging introduced.
- [ ] Per-feature evidence table written.
- [ ] Minimal regression tests added.

## Commit

`test: validate composite LMTO terms under GBT`
