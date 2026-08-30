# Luna Prompt — GBT-WP01: Exact q=0 Reduction

## Mission

Establish the first hard invariant of the GBT implementation:

\[
H_{GBT}(q=0)=H_{ordinary}
\]

at the **operator level**, not merely through similar energies.

This work package must not alter the core GBT gauge unless the new exact test demonstrates a reproducible failure and the correction is mathematically justified.

## Dependency

Read and require WP00 evidence first.

## Scope

Start with the minimal scalar-relativistic feature set:

- no SOC;
- no Hubbard V;
- no local-axis mode;
- preferably no CCOR/HOH for the first fixture;
- same structure-constant backend;
- same potential;
- same basis/order.

Then add a realistic bcc Fe fixed-potential fixture.

## Required tests

### Test 1 — q=0, collinear identity

Run ordinary periodic magnetic Hamiltonian and `gbt_single_q` with:

- `q_ss = 0`;
- zero cone/global reference rotation;
- same fixed potential.

Compare before diagonalization:

1. every corresponding directed real-space Hamiltonian block;
2. onsite blocks;
3. reciprocal `H(k)` at several generic k points.

Use norms and maximum elementwise residuals.

Acceptance target: roundoff scale, nominally relative `<= 1e-12` in double precision. If the existing numerical stack requires a slightly different threshold, justify it from observed roundoff on a control identity and keep it tight.

### Test 2 — q=0 global SU(2) rotation

Use a generic global cone/reference direction, e.g. `theta=37 deg`, no SOC.

Verify the finite-angle GBT Hamiltonian is globally unitarily related to the collinear ordinary Hamiltonian:

\[
H_{rot}(k)=U^\dagger H_0(k) U
\]

or the conventionally equivalent ordering used by the code.

At minimum require identical spectra to roundoff. Prefer a direct transformed-matrix comparison.

### Test 3 — q=0 multi-sublattice reference

Construct a small deterministic multi-sublattice fixture where sublattices have different fixed reference orientations but no translational spiral.

Compare GBT `q=0` to the ordinary primitive-cell non-collinear construction describing the same local axes.

This isolates endpoint-frame handling from q-dependent translation.

## Test design requirements

- Use generic k points in addition to Gamma/high-symmetry points.
- Reference side must use the ordinary code path, not GBT helpers.
- Record spin/basis permutations explicitly if one path orders states differently.
- Add a negative control if practical: deliberately alter one GBT phase/frame in test-only code and demonstrate the comparator catches it.

## Failure triage

If Test 1 fails, investigate only:

- representation dispatch;
- zero-q endpoint link;
- onsite-term duplication/omission;
- spin-block indexing;
- local collinear potential-factor contraction;
- differences in options silently enabled between ordinary and GBT paths.

Do not proceed to supercell/SCF explanations.

If Test 1 passes but Test 2 fails, focus on global-frame rotation ordering and density/reference-axis conventions.

If Test 1/2 pass but Test 3 fails, focus on sublattice endpoint reference frames.

## Deliverables

1. Lean automated tests in the existing test framework.
2. `docs/dev/GBT_REAUDIT_WP01_Q0.md` documenting:
   - fixtures;
   - exact comparisons;
   - thresholds;
   - residual tables;
   - any basis permutation;
   - pass/fail verdict.
3. If a defect is proven, a minimal fix plus a regression test.

## Completion checklist

- [ ] Minimal-feature q=0 ordinary/GBT real-space blocks compared.
- [ ] Generic-k reciprocal matrices compared.
- [ ] q=0 global finite-angle SU(2) covariance tested.
- [ ] q=0 multi-sublattice reference tested.
- [ ] Thresholds documented and roundoff-scale.
- [ ] Reference side does not reuse GBT helper logic.
- [ ] Any proven defect repaired minimally.
- [ ] Regression test fails under a meaningful negative control where practical.
- [ ] Evidence report written.
- [ ] No unrelated refactor or SCF tuning introduced.

## Commit

`test: enforce exact q-zero reduction of GBT`
