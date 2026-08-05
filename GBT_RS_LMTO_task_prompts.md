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

Dependency graph (revised after the WP2 representation-boundary finding):

```text
WP0 -> WP1 -> WP2/G2E -> WP3 -> WP4/G2O -> final G2 -> WP5
                                                          |
                                                          +-> WP6a HOH/overlap --+
                                                          +-> WP6b CCOR ---------+-> WP7 -> WP8 -> WP9 -> WP10
                                                          +-> WP6c Hubbard/ops --+
```

`G2E` is the already implemented canonical-energy evidence. `G2O` is the
production operator/frame slice that could not be meaningfully tested through
the legacy double-representation path. Only WP3 and WP4 are allowed before
final G2 closes; this is a gate repair, not permission to begin later work.

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

For the rotating-frame q=0 test, equality with the collinear operator is the
expected result. Pair it with a finite-q or relative-sublattice operator check
that fails if the probe was simply disabled. Report G2E and G2O separately if
the operator migration is not yet complete.

Completion checklist:
- [ ] Canonical electron count is implemented.
- [ ] Canonical eigenvalue EBAND is implemented.
- [ ] MPI and k-weight normalization are tested.
- [ ] Frozen-magnon MFT uses canonical EBAND.
- [ ] DOS settings do not move canonical EBAND beyond budget.
- [ ] Electron/energy cross-checks and error budgets are reported.
- [ ] G2 PASS/FAIL is stated.
```

## Next-session combined task — Implement GBT on primitive S/Sdot

**Suggested agent:** `gpt-5.6-sol`, xhigh reasoning; sequential work only.

This is the only permitted gate-repair task while final G2 remains open. It
executes WP3 and WP4 in order and stops at their gates; it does not authorize
WP5/6 cleanup unless explicitly requested after acceptance.

```text
Implement the GBT on S-level according to
GBT_RS_LMTO_completion_blueprint.md, completing WP3 and then WP4.

Fixed architecture:
- periodic_nc and explicit_texture retain the general q-agnostic moment logic
  in ham0m_nc;
- gbt_single_q uses a fixed collinear rotating-frame potential (normally +z
  for a ferromagnet);
- q, cone, and sublattice reference orientations enter only through the
  endpoint link on primitive directed S and Sdot blocks;
- do not rotate potential parameters or completed ee/Hamiltonians;
- recursion and reciprocal solvers must eventually consume the same operator.

First implement the explicit representation split. Preserve the user's
MX/MY/MZ and local_axis behavior, restore/retain safety checks, and ensure GBT
does not enter ham0m_nc's ordinary NC/texture rotation branches.

Require strux_backend='strux_lib' for gbt_single_q and fail early for the
legacy structure backend. Do not double or replace persistent lattice%sbar or
lattice%sdot. Keep the existing hmfind -> chbar_nc -> ham0m_nc path unchanged
for periodic_nc/explicit_texture. Add a parallel GBT-collinear path that fetches
raw complex orbital S/Sdot, builds a temporary nb x nb linked spinor structure
block, and contracts it with diagonal collinear rotating-frame potential
factors. strux_want_sdot is required only for a supported Sdot consumer.

Then implement one pure endpoint-frame/link helper using the complete physical
bond and the blueprint convention:
G_ab=(U_a0)^dagger Rz(q.d) U_b0.
Apply it to raw S and, where enabled and audited, Sdot before LMTO contraction.
Keep raw lattice S/Sdot unchanged. Build first-order ee from the gauged
primitive block and collinear potential parameters. Do not post-rotate ee,
eeo, eecc, or H(k).

Do not assume higher-order paths are already clean: reciprocal code still
reconstructs eeo/eecc, and CCOR still contains sublattice-angle logic. Keep
GBT+HOH/CCOR explicitly guarded unless their primitive-factor covariance is
proved in this task. Do not delete reciprocal GBT until the S-level oracle and
production first-order path pass; WP5 owns final deletion.

Tests must include independent dense spin matrices, raw/gauged bond dumps,
reverse-bond Hermiticity, the Stoner k+/-q/2 identity, q<->-q, non-cubic phase
units, q=0 common-cone identity, and a finite-q or relative-sublattice case that
cannot pass as a no-op. Rerun canonical electron/EBAND, DOS-independence, MPI,
MX/MY/MZ periodic-NC, local_axis, and fixed-potential frozen-magnon checks.

Stop and report after G3/G4 and G2O evidence. Do not start WP5 or reinterpret
physical bcc-Fe dispersion if any operator gate fails.

Completion checklist:
- [ ] Representation mode is explicit and independent of solver selection.
- [ ] periodic_nc/explicit_texture retain general ham0m_nc behavior.
- [ ] gbt_single_q uses collinear rotating-frame potential parameters.
- [ ] One helper owns endpoint frames, q.d, and the S/Sdot link.
- [ ] Raw S/Sdot are preserved; gauging precedes LMTO contractions.
- [ ] gbt_single_q rejects the legacy structure backend early.
- [ ] Persistent SBAR/SDOT layouts and the existing NC call chain are unchanged.
- [ ] A parallel complex nb x nb GBT-collinear assembler is tested.
- [ ] No potential or completed-Hamiltonian GBT rotation remains on the tested path.
- [ ] Unsupported HOH/CCOR/overlap combinations fail early.
- [ ] q=0 identity and non-no-op finite-q/relative-sublattice tests pass.
- [ ] Ordinary NC, local_axis, canonical EBAND, and MPI regressions pass.
- [ ] G3, G4, G2O, and final G2 recommendations are stated separately.
```

## WP3 — Separate periodic NC, GBT, and explicit textures

**Suggested agent:** `gpt-5.6-sol` or `gpt-5.6-terra`, high reasoning.

```text
Implement WP3 after G2E passes as the first half of the documented G2 operator
repair. Do not start WP4 until G3 passes.

Introduce explicit magnetic representation modes:
periodic_nc, gbt_single_q, explicit_texture.

Preserve the existing hmfind -> chbar_nc -> ham0m_nc path for ordinary NC and
textures. Extracting its wx0/wx1 algebra is optional and may not change its
interface or behavior. Add a site-indexed moment provider for explicit_texture.
Route build_bulkham/build_locham so arbitrary textures use per-site blocks;
reject per-type reuse unless represented as a true magnetic supercell. In
gbt_single_q, select a separate collinear rotating-frame path and do not
consume q/cone angles in ham0m_nc.

Add the per-sublattice reference-frame API but do not add the new GBT phase yet.
Preserve local_axis rotation, bounds checks, and all ordinary NC behavior.

Test two same-species sites with different moments and a small skyrmion-like texture.

Completion checklist:
- [x] Three representation modes are explicit.
- [x] Existing NC pair path remains q-agnostic; any extraction preserves its interface/behavior.
- [x] Existing hmfind/chbar_nc/ham0m_nc interfaces and NC behavior are preserved.
- [x] Explicit textures obtain moments by site identity.
- [x] GBT reference potential is explicitly collinear in the rotating frame.
- [x] q/cone inputs do not enter ham0m_nc in GBT mode.
- [x] Invalid per-type texture reuse is rejected.
- [x] Existing periodic-NC/local-axis regressions pass.
- [x] Same-species and skyrmion-like tests pass.
- [x] API/compatibility notes and G3 PASS/FAIL are reported.
```

## WP4 — Primitive S/Sdot GBT link

**Suggested agent:** `gpt-5.6-sol`, xhigh reasoning; integrator physics review required.

```text
Implement WP4 after G3 passes, following the physical-displacement and
rotating-frame conventions in the blueprint.

Add one gbt_bond_phase helper. In gbt_single_q only:
1. resolve endpoint reference frames U_a0,U_b0;
2. compute alpha from the complete directed Cartesian bond;
3. construct G_ab=(U_a0)^dagger Rz(alpha) U_b0;
4. apply G_ab to primitive S and audited Sdot blocks before LMTO contraction;
5. assemble with unchanged collinear rotating-frame potential parameters.

Support only strux_backend='strux_lib'; reject the legacy backend before
Hamiltonian construction. Keep persistent sbar/sdot orbital-sized and raw.
Use a per-bond temporary nb x nb spinor structure block and a new collinear
LMTO assembler parallel to ham0m_nc. Do not widen the existing NC interfaces.
Require strux_want_sdot only when an audited enabled term consumes Sdot.

Keep raw S/Sdot unchanged. Onsite potential terms remain collinear and carry no
translational phase. Do not phase periodic_nc or explicit_texture. Do not
rotate potential parameters, ee, or H(k) afterward. Bypass new arithmetic when
GBT is inactive. No other code may calculate q·bond.

Run dense structure-link, shifted-k, reverse-bond, q↔-q, theta=0, q=0
rotating-frame reduction, non-cubic, and non-no-op finite-q/relative-sublattice
tests. Rerun G2O after G4.

Completion checklist:
- [x] One phase helper owns every q·bond calculation.
- [x] Complete physical-displacement convention is used.
- [x] Endpoint frames and half-angle structure link are correct.
- [x] Raw S/Sdot are preserved and gauged before contraction.
- [x] Legacy structure backend fails early for GBT.
- [x] SBAR/SDOT storage is not doubled.
- [x] q=0 GBT-collinear assembly matches ham0m_nc at mom=+z for unequal endpoints.
- [x] Potential parameters and completed Hamiltonians are not GBT-rotated.
- [x] Onsite terms remain in the collinear rotating frame.
- [x] periodic_nc and explicit_texture remain ungauged.
- [x] All required algebraic/invariance tests pass.
- [x] A non-no-op test prevents false q=0 acceptance.
- [x] G4, G2O, final G2 PASS/FAIL and maximum errors are reported.
```

## WP5 — Remove reciprocal GBT and unify solvers

**Suggested agent:** `gpt-5.6-sol`, xhigh reasoning.

```text
Implement WP5 after G4 passes.

Delete fourier_transform_gbt_array/fourier_transform_gbt implementations,
interfaces, and calls, including reconstruction of eeo/eecc. Remove
gbt_kspace conditionals from reciprocal Hamiltonian assembly. Both solvers must
consume the same operator built from linked primitive structure blocks;
reciprocal space uses only the ordinary Fourier transform.

Invalidate all q/cone/reference-axis/potential-dependent caches before rebuilding Hamiltonians, eigensystems, DOS, or densities. Keep gbt_kspace only as a deprecated, non-physical input if immediate removal would break inputs.

Add pre-eigensolver Hermiticity assertions and compare RS/k-space outputs from identical bond dumps while converging recursion and k mesh.

Completion checklist:
- [x] All reciprocal GBT transforms/interfaces/calls are removed.
- [x] Reciprocal assembly uses the ordinary Fourier transform only.
- [x] gbt_kspace has no physics role.
- [x] q-dependent cache invalidation is documented and tested.
- [x] Pre-solver Hermiticity checks pass.
- [x] RS/k-space convergence comparison is reported.
- [x] Deletion inventory and G5 PASS/FAIL are stated.
```

## WP6a — HOH and overlap covariance

**Suggested agent:** `gpt-5.6-sol`, xhigh reasoning.

```text
Audit and implement the WP6 HOH/overlap slice after G5.

Trace eeo, eeoee, obarm, enim, the RS two-sweep HOH operator, reciprocal h-o-h construction, and generalized-overlap modes. Classify every object as onsite, directed bond, or composite. Apply GBT only at the earliest directed factor. Prove and test RS two-sweep versus reciprocal operator equivalence. Check H/O Hermiticity and overlap positive definiteness.

Do not phase an already contracted object without a derivation. If the formal overlap representation is incomplete, reject it in GBT mode explicitly.

Completion checklist:
- [x] eeo/eeoee/obarm/enim are classified.
- [x] HOH covariance is derived at primitive-bond level.
- [x] RS two-sweep and reciprocal operators agree.
- [x] H/O Hermiticity checks pass.
- [x] Supported overlap is positive definite.
- [x] Unsupported modes fail early.
- [x] Term map and partial G6 PASS/FAIL are reported.
```

## WP6b — CCOR covariance

**Suggested agent:** `gpt-5.6-sol`, xhigh reasoning.

```text
Audit and implement the WP6 CCOR slice after G5.

Remove both old full-angle transverse GBT rotations. Trace each CCOR contribution to its earliest directed structure-like factor and apply the same endpoint gauge convention as the main Hamiltonian before contractions. Onsite CCOR uses the shared rotating-frame reference axis and alpha=0.

Preserve ordinary noncollinear CCOR and explicit_texture behavior. Until every CCOR term has a dense oracle and directed-bond/k-space Hermiticity test, keep GBT+CCOR fatal.

Completion checklist:
- [x] Both full-angle GBT rotations are removed.
- [x] Every CCOR term is classified onsite/bond/composite.
- [x] Directed CCOR factors use the common gauge.
- [x] Ordinary NC and explicit-texture CCOR regressions pass.
- [x] Dense and Hermiticity oracles pass for enabled CCOR.
- [x] Unsupported CCOR combinations fail early.
- [x] Derivation, deletion list, and partial G6 PASS/FAIL are reported.
```

## WP6c — Hubbard, constraints, velocity, and torque

**Suggested agent:** `gpt-5.6-sol`, xhigh reasoning.

```text
Audit the remaining WP6 terms after G5: Hubbard-U/V, constraining fields, velocity, torque, derivative/downfolded terms, and SOC.

For each term classify onsite/bond/composite and state its frame. GBT onsite
terms remain in the collinear rotating frame with no translational phase;
ordinary NC/texture onsite terms retain their explicitly selected local frame.
Gauge directed intersite terms at their primitive factor. Derive
velocity/torque from the completed gauged operator and verify k finite
differences. Keep nonzero-q GBT+SOC fatal.

Preserve periodic_nc and explicit_texture behavior. Unsupported terms must fail before SCF.

Completion checklist:
- [x] Hubbard-U/V frame and gauge handling are explicit.
- [x] Constraint frame and energy bookkeeping are explicit.
- [x] Velocity/torque finite-difference tests pass.
- [x] Other derivative/downfolded terms are classified.
- [x] Nonzero-q GBT+SOC fails early.
- [x] Ordinary NC/texture regressions pass.
- [x] Feature matrix and partial G6 PASS/FAIL are reported.
```

**WP6c outcome (2026-08-05):** Hubbard-U (onsite) and SOC were already
correctly gated; this task adds explicit frame commentary and confirms both
guards by citation. Hubbard-V stays fatal (its diagonal-occupation proxy
carries no directed factor to gauge). Charge/spin/orbital velocity are
proven, by an independent k-derivative oracle (closed-form and finite-
difference), to inherit their gauge from the already-audited `ee` with no
code change needed. One real gap was found and fixed: `cond_type=spin_torque`
conductivity silently returned zero under GBT (`this%hxc` is never populated
on the GBT path) and is now fatal. Constraining fields have no functional
Hamiltonian/energy term in *any* representation mode in this codebase
(verified by reading both call sites) — recorded in `tests/KNOWN_ISSUES.md`,
not fixed here. Full evidence, term map, and regression results (including
5 pre-existing, WP6c-unrelated `strux_backend` fixture failures found while
running the full `scf` suite) are in `docs/dev/GBT_WP6C_REPORT.md`.
WP6c: **PASS**. Partial G6 through WP6a+WP6b+WP6c: **PASS for the audited
slices**. Final G6 remains open pending the WP6 integration task below.

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
