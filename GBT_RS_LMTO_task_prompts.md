# GBT completion: task prompts and delegation model

Companion to `GBT_RS_LMTO_completion_blueprint.md`. These prompts are deliberately short and ready to copy into coding-agent tasks.

## Delegation model

Use one **integrator** throughout. The integrator owns interfaces, resolves physics conventions, reviews every gate, and is the only agent allowed to merge or delete legacy paths.

Recommended model classes:

| Work | Model | Reasoning | Execution |
|---|---|---|---|
| WP0, WP1, WP3, WP8, WP10 | balanced coding model (`gpt-5.6-terra`) | high | one agent each |
| WP2, WP4, WP5, WP7 | frontier coding model (`gpt-5.6-sol`) | xhigh | sequential, physics-critical |
| WP6 term audits | frontier coding model (`gpt-5.6-sol`) | xhigh | up to three parallel agents after WP5 |
| WP9 validation | balanced model for runners; frontier model for diagnosis | high/xhigh | parallel cases, centralized interpretation |

Dependency graph:

```text
WP0 -> WP1 -> WP2 -> WP3 -> WP4 -> WP5
                                  |
                                  +-> WP6a HOH/overlap --+
                                  +-> WP6b CCOR ---------+-> WP7 -> WP8 -> WP9 -> WP10
                                  +-> WP6c Hubbard/ops --+
```

Rules:

- Do not run agents concurrently on the same source files.
- Every task starts from the last accepted gate, not merely the latest commit.
- A worker reports evidence; the integrator decides whether the gate passes.
- Failed gates return to the responsible task. Do not weaken tolerances or update golden files to hide a failure.
- Preserve unrelated worktree changes.
- Each task ends with: files changed, tests run, numerical evidence, unresolved risks, and gate recommendation.

---

## WP0 — Baseline and safety guards

**Suggested agent:** `gpt-5.6-terra`, high reasoning.

```text
Implement WP0 from GBT_RS_LMTO_completion_blueprint.md.

Goal: freeze trustworthy q=0/NC baselines and make known-wrong GBT combinations fail explicitly.

Do:
- record representative collinear, periodic-NC, explicit-texture, RS-MFT, and k-space GBT outputs;
- add nonzero-q guards for SOC, CCOR, Hubbard-V, and unverified overlap/features;
- force full BZ for nonzero q;
- check H(k), and O(k) when present, for Hermiticity before zheev/zhegv.

Do not change GBT physics yet. Label old GBT results diagnostic, not golden physics. Preserve q=0 behavior.

Completion checklist:
- [ ] Baseline inputs/outputs recorded and classified.
- [ ] Unsupported combinations fail early.
- [ ] Full BZ is forced for nonzero q.
- [ ] H/O Hermiticity checks run before eigensolution.
- [ ] Ordinary regressions are unchanged.
- [ ] G0 report includes commands, results, and remaining risks.
```

## WP1 — Independent algebraic oracles

**Suggested agent:** `gpt-5.6-terra`, high reasoning; frontier review.

```text
Implement WP1 after G0 passes.

Build tests independent of production GBT algebra:
- Cartesian/fractional phase equivalence for cubic, hexagonal, and triclinic cells;
- one-orbital Stoner k±q/2 identity;
- random dense 2x2 spin-matrix pair oracle with unequal endpoint species;
- reverse directed-bond Hermiticity;
- multi-sublattice case that detects double-counted basis/sublattice phase.

Use explicit dense spin matrices on at least one side; do not call the production Pauli helper for both sides. Reuse ideas from commits 2072fac and 7ca0c29 only after review.

Do not modify production routing. Required algebraic relative error: <1e-12.

Completion checklist:
- [x] Cell-convention phase tests pass.
- [x] Stoner shifted-k test passes.
- [x] Dense unequal-species oracle passes.
- [x] Reverse-bond Hermiticity passes.
- [x] Double-phase multi-sublattice test passes.
- [x] Maximum errors and conventions are reported.
- [x] G1 PASS/FAIL is stated.
```

## WP2 — Canonical k-space electron count and EBAND

**Suggested agent:** `gpt-5.6-sol`, xhigh reasoning.

```text
Implement WP2 after G1 passes.

Replace projected-DOS moments as the canonical k-space band energy. Add:
- electron count from eigenvalues, k weights, EF, and occupations;
- EBAND=sum(w_k*f_nk*epsilon_nk) using the same occupations;
- correct MPI reduction and explicit weight normalization;
- a consistent tetrahedron/Blöchl energy path, or a full-mesh Fermi-Dirac oracle for validation.

Route frozen_magnon_probe_energy through this canonical EBAND. Keep total-DOS integration and projected m0*m1 only as logged diagnostics. Ensure caches are not stale between probes.

Test electron conservation, DOS-grid/window independence, q=0 global-rotation invariance, and agreement with converged total-DOS integration.

Completion checklist:
- [ ] Canonical electron count is implemented.
- [ ] Canonical eigenvalue EBAND is implemented.
- [ ] MPI and k-weight normalization are tested.
- [ ] Frozen-magnon MFT uses canonical EBAND.
- [ ] DOS settings do not move canonical EBAND beyond budget.
- [ ] Electron/energy cross-checks and error budgets are reported.
- [ ] G2 PASS/FAIL is stated.
```

## WP3 — Separate periodic NC, GBT, and explicit textures

**Suggested agent:** `gpt-5.6-sol` or `gpt-5.6-terra`, high reasoning.

```text
Implement WP3 after G2 passes.

Introduce explicit magnetic representation modes:
periodic_nc, gbt_single_q, explicit_texture.

Extract ham0m_nc's ordinary wx0/wx1 dot/cross Pauli algebra into a q-agnostic pair kernel with explicit endpoint moments. Add a site-indexed moment provider for explicit_texture. Route build_bulkham/build_locham so arbitrary textures use per-site blocks; reject per-type reuse unless represented as a true magnetic supercell.

Do not add the new GBT phase yet. Preserve local_axis rotation and all ordinary NC behavior.

Test two same-species sites with different moments and a small skyrmion-like texture.

Completion checklist:
- [ ] Three representation modes are explicit.
- [ ] q-agnostic NC pair kernel is extracted.
- [ ] Explicit textures obtain moments by site identity.
- [ ] Invalid per-type texture reuse is rejected.
- [ ] Existing periodic-NC/local-axis regressions pass.
- [ ] Same-species and skyrmion-like tests pass.
- [ ] API/compatibility notes and G3 PASS/FAIL are reported.
```

## WP4 — Single GBT bond kernel

**Suggested agent:** `gpt-5.6-sol`, xhigh reasoning; integrator physics review required.

```text
Implement WP4 after G3 passes, following the physical-displacement convention in the blueprint.

Add one gbt_bond_phase helper. In gbt_single_q only:
1. resolve rotating-frame reference moments m_i0,m_j0;
2. compute alpha from the complete directed Cartesian bond;
3. assemble H_pair[m_i0,Rz(alpha)m_j0];
4. right-multiply by D(alpha) using half-angle Pauli formulas and temporaries.

Use the same reference-axis resolver for onsite pair terms, enim, and obarm. Do not phase explicit_texture. Bypass new arithmetic when GBT/cone inputs are inactive. No other code may calculate q·bond.

Run dense, shifted-k, reverse-bond, q↔-q, theta=0, q=0 rotation, and non-cubic tests.

Completion checklist:
- [ ] One phase helper owns every q·bond calculation.
- [ ] Complete physical-displacement convention is used.
- [ ] Endpoint moments and half-angle phase are correct.
- [ ] Onsite axes are consistent.
- [ ] periodic_nc and explicit_texture remain ungauged.
- [ ] All required algebraic/invariance tests pass.
- [ ] G4 PASS/FAIL and maximum errors are reported.
```

## WP5 — Remove reciprocal GBT and unify solvers

**Suggested agent:** `gpt-5.6-sol`, xhigh reasoning.

```text
Implement WP5 after G4 passes.

Delete fourier_transform_gbt_array/fourier_transform_gbt implementations, interfaces, and calls. Remove gbt_kspace conditionals from reciprocal Hamiltonian assembly. Both solvers must consume the same already-built real-space bonds; reciprocal space uses only the ordinary Fourier transform.

Invalidate all q/cone/reference-axis/potential-dependent caches before rebuilding Hamiltonians, eigensystems, DOS, or densities. Keep gbt_kspace only as a deprecated, non-physical input if immediate removal would break inputs.

Add pre-eigensolver Hermiticity assertions and compare RS/k-space outputs from identical bond dumps while converging recursion and k mesh.

Completion checklist:
- [ ] All reciprocal GBT transforms/interfaces/calls are removed.
- [ ] Reciprocal assembly uses the ordinary Fourier transform only.
- [ ] gbt_kspace has no physics role.
- [ ] q-dependent cache invalidation is documented and tested.
- [ ] Pre-solver Hermiticity checks pass.
- [ ] RS/k-space convergence comparison is reported.
- [ ] Deletion inventory and G5 PASS/FAIL are stated.
```

## WP6a — HOH and overlap covariance

**Suggested agent:** `gpt-5.6-sol`, xhigh reasoning.

```text
Audit and implement the WP6 HOH/overlap slice after G5.

Trace eeo, eeoee, obarm, enim, the RS two-sweep HOH operator, reciprocal h-o-h construction, and generalized-overlap modes. Classify every object as onsite, directed bond, or composite. Apply GBT only at the earliest directed factor. Prove and test RS two-sweep versus reciprocal operator equivalence. Check H/O Hermiticity and overlap positive definiteness.

Do not phase an already contracted object without a derivation. If the formal overlap representation is incomplete, reject it in GBT mode explicitly.

Completion checklist:
- [ ] eeo/eeoee/obarm/enim are classified.
- [ ] HOH covariance is derived at primitive-bond level.
- [ ] RS two-sweep and reciprocal operators agree.
- [ ] H/O Hermiticity checks pass.
- [ ] Supported overlap is positive definite.
- [ ] Unsupported modes fail early.
- [ ] Term map and partial G6 PASS/FAIL are reported.
```

## WP6b — CCOR covariance

**Suggested agent:** `gpt-5.6-sol`, xhigh reasoning.

```text
Audit and implement the WP6 CCOR slice after G5.

Remove both old full-angle transverse GBT rotations. Trace each CCOR contribution to its earliest directed structure-like factor and apply the same endpoint gauge convention as the main Hamiltonian before contractions. Onsite CCOR uses the shared rotating-frame reference axis and alpha=0.

Preserve ordinary noncollinear CCOR and explicit_texture behavior. Until every CCOR term has a dense oracle and directed-bond/k-space Hermiticity test, keep GBT+CCOR fatal.

Completion checklist:
- [ ] Both full-angle GBT rotations are removed.
- [ ] Every CCOR term is classified onsite/bond/composite.
- [ ] Directed CCOR factors use the common gauge.
- [ ] Ordinary NC and explicit-texture CCOR regressions pass.
- [ ] Dense and Hermiticity oracles pass for enabled CCOR.
- [ ] Unsupported CCOR combinations fail early.
- [ ] Derivation, deletion list, and partial G6 PASS/FAIL are reported.
```

## WP6c — Hubbard, constraints, velocity, and torque

**Suggested agent:** `gpt-5.6-sol`, xhigh reasoning.

```text
Audit the remaining WP6 terms after G5: Hubbard-U/V, constraining fields, velocity, torque, derivative/downfolded terms, and SOC.

For each term classify onsite/bond/composite and state its frame. Onsite terms use the common rotating reference axis with no translational phase. Gauge directed intersite terms at their primitive factor. Derive velocity/torque from the completed gauged operator and verify k finite differences. Keep nonzero-q GBT+SOC fatal.

Preserve periodic_nc and explicit_texture behavior. Unsupported terms must fail before SCF.

Completion checklist:
- [ ] Hubbard-U/V frame and gauge handling are explicit.
- [ ] Constraint frame and energy bookkeeping are explicit.
- [ ] Velocity/torque finite-difference tests pass.
- [ ] Other derivative/downfolded terms are classified.
- [ ] Nonzero-q GBT+SOC fails early.
- [ ] Ordinary NC/texture regressions pass.
- [ ] Feature matrix and partial G6 PASS/FAIL are reported.
```

## WP6 integration — Close the feature gate

**Suggested owner:** integrator with `gpt-5.6-sol`, xhigh review.

```text
Integrate WP6a/b/c after all three report.

Resolve shared-interface conflicts, run the full regression suite, and create one authoritative feature matrix covering ee, onsite terms, HOH, overlap, CCOR, Hubbard, constraints, velocity/torque, and SOC. Every enabled GBT combination needs an operator oracle plus bond and k-space Hermiticity evidence; every other combination must fail early.

Do not weaken tests or enable partially audited paths.

Completion checklist:
- [ ] WP6a/b/c interfaces are reconciled.
- [ ] One authoritative feature matrix exists.
- [ ] Every enabled combination has oracle evidence.
- [ ] Every unsupported combination fails early.
- [ ] Full regression suite passes.
- [ ] Conflicts/limitations are documented.
- [ ] Final G6 PASS/FAIL is stated.
```

## WP7 — Common rotating-frame density and SCF

**Suggested agent:** `gpt-5.6-sol`, xhigh reasoning; no parallel source edits.

```text
Implement WP7 after G6 passes.

Define one per-site/per-l rotating-frame 2x2 spin density contract. Implement equivalent producers from RS Green functions and k-space eigenvectors/weights/occupations. Project to radial up/down only after accumulation using an explicit shared axis.

Support two policies:
- constrained_spiral: fixed reference direction; mix charge and longitudinal magnitude; report transverse residual/torque;
- relaxed_reference: mix the full rotating-frame Cartesian moment while retaining the single-q ansatz.

Never use stale potential%mom as an implicit projection axis. Reconstruct lab-frame moments only for output/comparison. Check density Hermiticity, eigenvalues, trace, electron count, and |m|<=n.

Compare RS/k-space SCF and a q=0 cone.

Completion checklist:
- [ ] One full rotating-frame density contract exists.
- [ ] RS and k-space producers populate the same object.
- [ ] Radial projection uses an explicit current axis.
- [ ] Constrained and relaxed policies are distinct.
- [ ] Density physicality assertions pass.
- [ ] q=0 cone invariance passes.
- [ ] RS/k-space SCF comparison and G7 PASS/FAIL are reported.
```

## WP8 — Little-group symmetry and q lifecycle

**Suggested agent:** `gpt-5.6-terra`, high reasoning.

```text
Implement WP8 after G7 passes. Keep full BZ as the oracle.

Add/review the q little-group reduction. Cache keys must include lattice, mesh, offset, q, and symmetry policy. In multi-q sweeps rebuild per q or use the common subgroup; never reuse row 1's mesh blindly. Clear all q-dependent eigensystem/DOS/density state.

Test axial and generic q, shifted meshes, q and -q, reordered q lists, and reduced versus full BZ. Reuse b1caf4b only after reviewing its per-q ownership.

Completion checklist:
- [ ] Full-BZ oracle remains available.
- [ ] Little-group cache includes q and all mesh state.
- [ ] Multi-q sweeps rebuild or use a common subgroup.
- [ ] All q-dependent solver state is invalidated.
- [ ] Axial/generic/shifted-mesh tests pass.
- [ ] q-list reordering leaves results unchanged.
- [ ] Full/reduced comparisons and G8 PASS/FAIL are reported.
```

## WP9 — Physical validation

**Suggested execution:** up to three `gpt-5.6-terra` runners; one `gpt-5.6-sol` diagnosis owner.

```text
Run WP9 after G8 passes; do not change physics during the first run.

Validate in order: analytic Stoner; dense bond oracle; commensurate GBT/supercell MFT then SCF; RS/k-space observables and EBAND; converged bcc-Fe Gamma-H MFT; q↔-q and H smoothness; 5-20 degree cone scaling; stiffness versus shell-converged LKAG J(q); constrained SCF supercells; then multi-sublattice Goldstone behavior.

Converge each method below the target Delta-E scale. Separate numerical uncertainty from physics disagreement. Do not regenerate goldens or impose a post-hoc Goldstone shift to pass.

Completion checklist:
- [ ] Analytic and dense-oracle tests pass.
- [ ] Commensurate MFT and SCF supercell tests pass.
- [ ] RS/k-space observables and EBAND agree when converged.
- [ ] bcc-Fe Gamma-H and q-symmetry tests pass.
- [ ] Small-cone scaling is within the declared budget.
- [ ] LKAG comparison includes shell convergence/error bars.
- [ ] Multi-sublattice Goldstone test is unshifted and explained.
- [ ] Inputs, raw data, convergence, uncertainty, and failures are archived.
- [ ] G9 PASS/FAIL is stated.
```

## WP10 — Legacy deletion and documentation

**Suggested agent:** `gpt-5.6-terra`, high reasoning; integrator performs final merge.

```text
Implement WP10 only after G9 passes.

Delete:
- absolute-position bulk GBT in ham0m_nc;
- all fourier_transform_gbt* code;
- both CCOR full-angle GBT rotations;
- mutually exclusive cone/sublattice legacy routing;
- the physics role of gbt_kspace, completing deprecation/removal;
- obsolete comments and known-wrong GBT golden data.

Keep and document periodic_nc and the site-indexed explicit_texture/skyrmion path. Add a repository check allowing q·bond only in gbt_bond_phase. Update user/developer docs and examples.

Run the complete test matrix and repository searches.

Completion checklist:
- [ ] Absolute-position bulk GBT is removed.
- [ ] Reciprocal GBT code is removed.
- [ ] CCOR full-angle GBT code is removed.
- [ ] gbt_kspace deprecation/removal is complete.
- [ ] Explicit-texture/skyrmion support remains tested and documented.
- [ ] Only gbt_bond_phase calculates q·bond.
- [ ] Obsolete comments/goldens are removed or relabelled.
- [ ] Complete test matrix passes.
- [ ] Removal/migration notes and G10 PASS/FAIL are reported.
```

---

## Integrator review prompt

Use after every work package:

```text
Review the submitted WP against GBT_RS_LMTO_completion_blueprint.md and its gate.

Check physics convention, frame ownership, directed-bond covariance, q=0/NC regressions, cache invalidation, and whether tests are independent of production algebra. Inspect source and rerun the smallest decisive tests. Reject golden-file updates that merely encode changed output.

Review checklist:
- [ ] Physics and gauge conventions match the blueprint.
- [ ] Frame ownership and cache invalidation are explicit.
- [ ] Tests are independent and decisive.
- [ ] q=0 and ordinary NC regressions pass.
- [ ] No golden update hides a failure.
- [ ] Blocking findings cite file/symbol evidence.
- [ ] Rerun commands/results are recorded.
- [ ] Gate PASS/FAIL is stated.
- [ ] The next allowed dependency-graph task is named.
```
