# Luna Prompt — GBT-WP11: LKAG J(q) and Frozen-Spiral Consistency

## Mission

Close the independent physics loop between the GBT frozen-spiral energy curvature and the code's established Green-function/LKAG exchange interactions.

This is the last required gate before SCF-specific debugging.

## Dependency

WP10 must provide a converged small-q harmonic regime.

## Critical requirement: derive normalization, do not assume it

The code currently contains frozen-magnon normalization formulas of the Halilov type. Do not assume that a factor such as `4/(M sin^2 theta)` is automatically consistent with the exchange-Hamiltonian convention used by the LKAG implementation.

Derive the relation from the exact conventions in this repository:

- definition/sign of `J_ij`;
- whether Hamiltonian uses `-J e_i dot e_j`, `-1/2 sum`, or another counting convention;
- ordered versus unordered pairs;
- magnetic moment normalization;
- whether magnon energy is expressed in Ry, eV, meV, or includes a moment/g factor;
- sublattice normalization.

Document the derivation before comparing numbers.

## Construct J(q)

For the one-sublattice case:

\[
J(\mathbf q)=\sum_R J_{0R}e^{i\mathbf q\cdot R}.
\]

For a centrosymmetric scalar-relativistic ferromagnet this should reduce to the appropriate cosine/even form.

The frozen-spiral energy must be related to

\[
J(0)-J(q)
\]

with the repository-specific prefactor derived above.

## Multi-sublattice case

Construct the exchange matrix `J_ab(q)` and compare the relevant acoustic eigenvalue/branch, not a scalar total-J shortcut.

Account for sublattice moment factors exactly as required by the chosen spin Hamiltonian convention.

## Numerical convergence

Ensure the real-space exchange sum is itself converged with respect to interaction range / neighbor shell support. A mismatch between truncated LKAG J(R) and a more extended electronic GBT operator is not evidence of a GBT defect.

Report convergence of both sides.

## Benchmarks

Primary:

- bcc Fe.

Secondary where robust:

- fcc Ni;
- B2 FeCo for multi-sublattice acoustic response.

## Compare more than one number

Compare:

1. `DeltaE(q)` curve in the small-q regime;
2. fitted quadratic curvature/stiffness;
3. if feasible, a broader q range where the Heisenberg mapping remains meaningful;
4. q=0 Goldstone consistency.

## Interpretation

If bare MFT agrees with standard LKAG but corrected constrained MFT/SCF differs, do not automatically classify the latter as wrong. Document the distinction between bare LKAG/MFT response and renormalized/constrained energy landscapes in the Bruno/Jacobsson context.

The purpose of this gate is to establish internal consistency of the **same response convention** first.

## Deliverables

- reusable `J(q)` analysis path or test helper;
- normalization derivation;
- convergence evidence;
- `docs/dev/GBT_REAUDIT_WP11_LKAG_JQ.md`.

## Completion checklist

- [ ] Exact exchange-Hamiltonian convention located in code.
- [ ] Pair-counting/sign convention derived.
- [ ] Frozen-magnon prefactor derived from same convention.
- [ ] J(q) constructed independently from LKAG J(R).
- [ ] LKAG interaction-range convergence checked.
- [ ] GBT k/q convergence carried over from WP10.
- [ ] bcc Fe curve/curvature compared.
- [ ] Multi-sublattice matrix treatment implemented or clearly scoped.
- [ ] Bare versus constrained-response interpretation documented.
- [ ] Quantitative pass/fail tolerance declared and evaluated.

## Commit

`test: close GBT frozen-magnon response against LKAG Jq`
