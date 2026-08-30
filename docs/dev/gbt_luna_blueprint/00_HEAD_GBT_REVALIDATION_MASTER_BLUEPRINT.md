# RS-LMTO-ASA Generalized Bloch Theorem (GBT) Revalidation and Repair Campaign

**Target repository:** `https://github.com/rslmtoasa/rslmtoasa`

**Target branch:** `fable_v3`

**Purpose:** master blueprint for a validation-driven re-audit, repair, and production qualification of the generalized Bloch theorem / spin-spiral machinery in the real-space LMTO-ASA code.

**Companion artifacts:** `LUNA_GBT_WP00` through `LUNA_GBT_WP12` are sliced implementation prompts. Each work package is intended to be independently executable by a coding agent while obeying the dependencies and gates defined here.

---

## 1. Executive decision

Do **not** restart the core GBT implementation by replacing the current primitive-structure-constant gauge with a post-hoc rotation of completed Hamiltonian blocks.

The present `fable_v3` architecture is, in its core SOC-free form, consistent with both:

1. the LMTO-ASA spin-spiral formulation in which the local potential function remains local while the structure matrix acquires the spin-spiral dependence through the equivalents of `S(k+q/2)` and `S(k-q/2)`; and
2. the generalized Bloch construction used in mature non-LMTO implementations such as ELK, where the two spin channels are represented at shifted crystal momenta `k +/- q/2`.

The present real-space construction

\[
S_{ab}(\mathbf R) \rightarrow S_{ab}(\mathbf R)\otimes G_{ab}(\mathbf R),
\]

with

\[
G_{ab}(\mathbf R)=
(U_a^0)^\dagger
\exp\!\left[-\frac{i}{2}(\mathbf q\cdot\mathbf R)\sigma_z\right]
U_b^0,
\]

followed by the ordinary LMTO contractions and ordinary Fourier transform, is therefore the **working hypothesis to be proved**, not discarded.

The campaign must instead restart the **validation ladder** from first principles and use exact operator identities to separate five logically distinct layers:

- GBT algebra / SU(2) gauge;
- LMTO composite-term covariance;
- reciprocal-space Fourier and Brillouin-zone quadrature;
- rotating-frame density and constraining-field semantics;
- fully self-consistent radial-potential feedback.

No SCF debugging is allowed to drive changes to the core GBT operator until the earlier exact gates have been exercised.

---

## 2. Current audited baseline in `fable_v3`

The current branch already contains a substantial GBT implementation. Treat it as an existing production candidate that lacks decisive qualification, not as an empty prototype.

Relevant current components include at least:

- `source/gbt_structure.f90`
  - `gbt_bond_phase`
  - `gbt_reference_frame`
  - `gbt_endpoint_link`
  - `gbt_lift_orbital_block`
  - `gbt_contract_collinear`
- `source/hamiltonian_build.f90`
  - magnetic representation dispatch including `gbt_single_q`
  - `build_gbt_bulkham`
  - GBT-specific primitive-link assembly
  - local constraining-field insertion into the GBT onsite spinor
- `source/self.f90`
  - reciprocal-space density / magnetic-moment update
  - constraint update
  - total-energy reporting including current `constraint_energy`
- `source/abspinlib/constrain.f90`
  - iterative constraining-field controller
- `source/calculation.f90`
  - frozen-magnon post-processing
  - `mode='mft'`
  - `mode='scf'`
  - same-q `theta=0` gauge-energy diagnostic in k-space
- `GBT_RS_LMTO_completion_blueprint.md`
  - previous design/completion blueprint
- archived GBT validation reports such as WP9 / VAL-16 / VAL-17

The line numbers may evolve. Luna must locate routines by symbol/name and inspect the current branch before modifying code.

### 2.1 What appears conceptually correct already

The audit currently assigns high confidence to:

- a distinct `gbt_single_q` magnetic representation;
- applying the spiral link to primitive directed `S` blocks;
- treating raw `S` as immutable;
- contracting the linked structure constants with local rotating-frame LMTO potential factors;
- applying the corresponding link to raw `Sdot` when CCOR requires it;
- allowing HOH-like terms to inherit covariance from already-gauged lower-level operators instead of receiving an independent extra phase;
- using an ordinary downstream Fourier transform rather than a special GBT Fourier transform;
- rejecting SOC and other features whose covariance is not yet established;
- inserting constraining fields as local onsite Pauli terms, with no translational spiral phase.

### 2.2 What is **not** considered validated

The following remain open until this campaign closes them:

- exact `q=0` operator equivalence to ordinary calculations;
- exact `k +/- q/2` shifted-momentum identity;
- matched commensurate-supercell equivalence using the same finite operator;
- full covariance of HOH, CCOR, generalized overlap, and other composite terms;
- one unambiguous rotating-frame/lab-frame contract across Hamiltonian, density, moments, and constraints;
- physical versus algorithmic constraining-field energy bookkeeping;
- a corrected constrained frozen-magnon mode with frozen ordinary potential but q-dependent converged constraint field;
- harmonic cone-angle scaling after removal of quadrature artifacts;
- converged small-q quadratic behavior;
- agreement with independently calculated LKAG `J(q)` / stiffness;
- unrestricted/constrained GBT SCF.

---

## 3. Non-negotiable campaign rules

These rules override expedient fixes.

### Rule A — Freeze the core gauge unless an exact test fails

Do not change `gbt_endpoint_link`, the basic `S -> S x G` architecture, or replace it with a final-Hamiltonian rotation merely because an energy curve looks wrong.

A core-gauge change requires a failing exact operator-level identity plus a derivation showing the correction.

### Rule B — Fixed-potential operator tests precede SCF

The order is:

1. exact algebra;
2. exact matrix equivalence;
3. matched-supercell equivalence;
4. density/frame semantics;
5. constraints;
6. energy differences;
7. LKAG consistency;
8. only then SCF repair.

### Rule C — SOC-free scalar-relativistic GBT is the reference problem

Initially disable:

- SOC;
- Hubbard V;
- unaudited local-axis machinery;
- any feature not shown to transform covariantly.

Add composite terms one by one only after the minimal operator passes.

### Rule D — Never compensate a structural mismatch with an empirical energy shift

If a primitive GBT and explicit supercell contain different finite bond/operator sets, the comparison is invalid. Rebuild a matched operator; do not fit offsets.

### Rule E — Do not let same-q `theta=0` subtraction hide operator errors

The diagnostic

\[
E(q,\theta)-E(q,0)
\]

is useful for cancelling finite-k-grid pure-gauge quadrature error. It is not proof that the operator is correct. Exact matrix identities must pass first.

### Rule F — Constraining fields are physics, not merely an SCF stabilizer

Separate:

- a density/mixing policy that keeps a prescribed reference direction;
- an actual transverse constraining field acting as a Lagrange multiplier/control field;
- any numerical penalty used by the field controller.

Do not call these the same object.

### Rule G — Preserve established code architecture

Keep the modern modular/class-like Fortran style already used by the codebase. Prefer small modules, explicit interfaces, constructors/destructors where appropriate, scoped helpers, and lean tests. Avoid invasive edits to legacy atomic/radial LMTO routines unless an exact diagnostic demonstrates that they are the source of the error.

---

## 4. Core mathematics that the implementation must satisfy

### 4.1 Generalized translation

Without SOC, a single-q spiral can be represented by the combined lattice translation and spin rotation

\[
\hat T_{\mathbf R}\hat U(\mathbf q\cdot\mathbf R).
\]

The generalized Bloch spinor can be written such that the two spin sectors carry shifted crystal momenta:

\[
\psi_{\mathbf k}(\mathbf r)=
\begin{pmatrix}
 e^{i(\mathbf k+\mathbf q/2)\cdot\mathbf r}\,u_{\uparrow\mathbf k}(\mathbf r)\\
 e^{i(\mathbf k-\mathbf q/2)\cdot\mathbf r}\,u_{\downarrow\mathbf k}(\mathbf r)
\end{pmatrix},
\]

up to the code's Fourier/sign convention.

This is the most important external cross-code oracle.

### 4.2 LMTO real-space covariance

For a schematic two-center LMTO contribution

\[
h_{ij}=W_i S_{ij} W_j,
\]

with local spin-dependent potential factors `W_i`, transforming to the local rotating frames gives

\[
\bar h_{ij}=
W_i^{\rm loc}
\left[S_{ij}\otimes(U_i^\dagger U_j)\right]
W_j^{\rm loc}.
\]

Thus the natural primitive object carrying the nonlocal spin gauge is the structure constant.

For the single-q frame

\[
U_{na}=R_z(\mathbf q\cdot\mathbf R_n)U_a^0,
\]

the directed link is

\[
G_{na,mb}=U_{na}^\dagger U_{mb}.
\]

The implementation must obey

\[
G_{ji}(-\mathbf d)=G_{ij}^{\dagger}(\mathbf d).
\]

### 4.3 Pure-gauge shifted-k identity

For a collinear one-sublattice reference and `theta=0`, the GBT problem must reduce exactly to ordinary spin channels at shifted k points:

\[
H_{GBT}(\mathbf k,\mathbf q,0)
\sim
H_{\uparrow}(\mathbf k+s_\uparrow \mathbf q/2)
\oplus
H_{\downarrow}(\mathbf k+s_\downarrow \mathbf q/2),
\]

where the signs are fixed by the code's Fourier convention.

This must be verified at the matrix/eigenvalue level to roundoff at generic `k` and `q`.

### 4.4 q=0 covariance

At `q=0`, GBT must reduce exactly to the ordinary problem. For a global finite cone in the SOC-free case, the Hamiltonian is related only by a global SU(2) rotation and therefore has the same spectrum.

### 4.5 Closed-loop pure-gauge identity

For a pure single-q gauge, any closed lattice loop must satisfy

\[
G_{ij}G_{jk}\cdots G_{li}=I.
\]

This is a powerful implementation-level check on endpoint ordering and displacement conventions.

---

## 5. Constraining-field physics and required distinctions

### 5.1 Why a field is needed

A prescribed spiral is not generally a stationary solution of the unconstrained Kohn-Sham functional. If the local exchange field produces a transverse torque on a moment that is required to remain along a prescribed direction, a local constraining field is required.

The field acts as the Lagrange multiplier/control field enforcing

\[
\mathbf m_i \times \mathbf e_i^{\rm target}\rightarrow 0.
\]

A vanishing field is meaningful only if the chosen state is already torque-free; therefore a validation fixture must deliberately include a state requiring a nonzero field.

### 5.2 Bruno and renormalized force theorem

Bruno's renormalized magnetic force theorem places the relation between the bare magnetic-force-theorem response and the constrained response on a formal basis. The campaign must use this literature to interpret results, not merely to tune code.

### 5.3 Jacobsson/ELK corrected frozen magnons

The corrected frozen-magnon strategy gives a particularly useful intermediate mode:

1. obtain a reference ordinary effective/exchange-correlation potential;
2. keep the ordinary potential frozen;
3. for each `q`, converge the constraining field required to realize the prescribed moments;
4. perform the final one-shot energy/eigenvalue evaluation using the fixed ordinary potential plus the q-specific converged constraining field.

This must be represented as a separate mode from both ordinary MFT and fully self-consistent GBT.

### 5.4 Energy bookkeeping

Do not assume that a numerical controller penalty is part of the physical DFT energy.

The code must eventually expose, conceptually if not through these exact names:

- `E_DFT_physical`;
- `E_constraint_lagrange` if needed for diagnostics;
- `E_constraint_penalty` or controller merit function;
- `E_reported_auxiliary` if the optimizer uses one.

The physical spin-spiral energy used in comparisons must be identified by derivation and be invariant to harmless changes of controller tuning after the physical constraint is converged.

---

## 6. Work-package dependency graph

Execute the Luna slices in this order unless a prompt explicitly states otherwise:

| WP | Title | Primary gate | Depends on |
|---|---|---|---|
| WP00 | Current-state map and evidence ledger | audited call graph | none |
| WP01 | Exact q=0 reduction | ordinary == GBT | WP00 |
| WP02 | Gauge, Hermiticity, shifted-k identity | ELK/LMTO exact oracle | WP01 |
| WP03 | Composite LMTO covariance | HOH/CCOR/etc. exact | WP02 |
| WP04 | Matched commensurate supercell | exact spectrum folding | WP02; preferably WP03 |
| WP05 | Central GBT frame contract | density/moment covariance | WP04 |
| WP06 | Constraint-field covariance and diagnostics | nonzero-field fixture | WP05 |
| WP07 | Constraint energy semantics | physical energy identified | WP06 |
| WP08 | Corrected constrained frozen magnons | MFT / cMFT / SCF split | WP07 |
| WP09 | Cone-angle and k-grid convergence | harmonic regime | WP08 |
| WP10 | Small-q quadratic limit | stable `q^2` curvature | WP09 |
| WP11 | LKAG `J(q)` consistency | independent physics closure | WP10 |
| WP12 | GBT SCF diagnosis and repair | production qualification | WP11 |

**Hard stop:** WP12 may not alter the SCF/radial machinery to “repair GBT” before WP01-WP11 have produced explicit pass/fail evidence.

---

## 7. Acceptance hierarchy

### Gate G1 — exact q=0 operator reduction

For representative simple and full LMTO fixtures,

\[
\frac{\|H_{GBT}-H_{ordinary}\|}{\max(1,\|H_{ordinary}\|)}
\lesssim 10^{-12}
\]

in double precision, unless a clearly identified library operation gives a slightly looser but still roundoff-scale threshold.

Compare both directed real-space blocks and reciprocal-space matrices.

### Gate G2 — exact shifted-k identity

At generic irrational-looking test points, compare the GBT collinear spectrum/matrix with ordinary `k +/- q/2` evaluations. Require roundoff-level agreement.

### Gate G3 — bond reversal and Hermiticity

Require roundoff-level:

\[
G_{ji}(-R)-G_{ij}^{\dagger}(R)=0,
\]

and

\[
H(k)-H^{\dagger}(k)=0.
\]

### Gate G4 — matched commensurate-supercell equivalence

Use exactly the same primitive finite bond table on both sides. For commensurate `q`, compare the folded eigenvalue multisets. Require roundoff-level agreement at selected supercell k points.

### Gate G5 — frame covariance

Known rotating-frame density/moments transformed to a commensurate explicit supercell must reproduce the lab-frame spiral site by site.

### Gate G6 — real constraint test

Use a deliberately nonstationary prescribed state. Require converged orientation residual and a nonzero constraint field. Then verify field covariance against the explicit supercell.

### Gate G7 — energy semantics

Physical constrained energies must converge independently of controller penalty/gain values within a sensible stable range.

### Gate G8 — three frozen-magnon modes

Produce clearly distinct and reproducible:

- bare MFT;
- corrected constrained MFT;
- fully constrained SCF.

### Gate G9 — harmonic cone regime

For sufficiently small angles,

\[
\Delta E(q,\theta)/\sin^2\theta
\]

must form a plateau after k-grid convergence.

### Gate G10 — small-q even quadratic limit

Fit the symmetrized energy to

\[
Aq^2+Bq^4,
\]

and require the odd component to be numerical noise in SOC-free centrosymmetric benchmarks.

### Gate G11 — LKAG closure

The frozen-spiral curvature must agree with the independently constructed LKAG `J(q)` result after deriving and documenting the exact normalization/prefactor used by this code.

Only after G11 passes may GBT SCF discrepancies be diagnosed as SCF-specific.

---

## 8. Recommended benchmark systems

### 8.1 Minimal deterministic model

Use a tiny one-orbital or reduced-orbital scalar-relativistic model for exact matrix tests. The purpose is diagnosability, not realism.

### 8.2 bcc Fe

Use for:

- conventional ferromagnetic baseline;
- q=0 and shifted-k tests with a real LMTO potential;
- small-q stiffness;
- LKAG cross-check.

Do **not** use Fe alone to validate constraints because ordinary frozen-magnon MFT can benefit from fortuitous error cancellation.

### 8.3 fcc Ni

Use as a single-sublattice, more constraint-sensitive benchmark.

### 8.4 B2 FeCo

Use as a multi-sublattice benchmark for:

- reference-frame logic;
- nontrivial constraint fields;
- multi-sublattice frozen-magnon response;
- exchange-matrix / acoustic-mode comparison.

### 8.5 Commensurate q points

At minimum include simple exact fractions such as:

- `(0,0,1/2)`;
- `(0,0,1/3)`;
- optionally `(0,0,1/4)`;

in the code's documented Cartesian `2*pi/alat` convention, choosing lattice directions where the mapping is unambiguous.

---

## 9. Required reporting discipline

Every work package must create or update an evidence markdown file under an appropriate `docs/dev/` or test-evidence location.

Each report must include:

- exact commit tested;
- compiler/build configuration;
- feature flags;
- system/input fixture;
- q-vector convention and units;
- cone angles;
- k-mesh;
- exact quantities compared;
- absolute and relative residuals;
- pass/fail thresholds chosen **before** examining results where practical;
- whether the failure is algebraic, representation-level, numerical/quadrature, constraint-level, or SCF-level;
- explicit statement of what was **not** established.

Do not write “GBT validated” from a visually plausible dispersion alone.

---

## 10. Testing and implementation conventions

### 10.1 Prefer exact unit/functional tests before expensive end-to-end tests

The most valuable tests in this campaign compare matrices and spectra that should be mathematically identical. These should be fast enough for routine regression testing.

### 10.2 Never self-validate a formula with the same helper on both sides

For example, the shifted-k oracle should use the ordinary Hamiltonian path on the reference side, not reconstruct both sides through `gbt_endpoint_link`.

### 10.3 Use generic k/q points

Avoid only high-symmetry points. Include non-special deterministic vectors so sign and coordinate mistakes cannot hide behind symmetry.

### 10.4 Preserve raw diagnostic outputs

For energy convergence tasks, emit machine-readable tables with raw energies before plotting or fitting.

### 10.5 Keep negative controls

Where useful, tests should demonstrate that a deliberate sign flip, omitted half-phase, or reversed endpoint order is detected. Negative controls are especially useful for GBT phase tests.

---

## 11. Decision tree when a test fails

### If q=0 fails

Investigate only:

- representation dispatch;
- endpoint frames at zero q;
- potential-factor sign/spin indexing;
- accidental double application of a local rotation;
- onsite terms omitted/duplicated in GBT.

### If q=0 passes but shifted-k fails

Investigate:

- q sign;
- factor of `1/2` in the spinor phase;
- Fourier transform sign;
- physical displacement versus cell-translation displacement;
- endpoint ordering;
- spin-block ordering.

### If shifted-k passes but matched supercell fails

Do **not** modify the GBT phase. Investigate:

- supercell construction;
- folding/unfolding;
- periodic images;
- duplicate/missing directed bonds;
- basis permutation.

### If operator tests pass but frame covariance fails

Investigate:

- whether moments are in rotating or lab frame;
- conversion direction/sign;
- whether `potential%mom` is being used as both state and reference axis;
- density producer versus radial projection semantics.

### If frame tests pass but constraints fail

Investigate:

- target/reference direction;
- constraint-field frame;
- sign of the Pauli coupling;
- update/controller;
- whether the test state actually requires a field.

### If all fixed-potential tests pass but SCF fails

Only then inspect:

- radial spin-channel projection;
- charge/magnetization mixing;
- update of `C`, `Delta`, `E_nu`, etc.;
- longitudinal versus transverse feedback;
- constraint/physical-energy definitions.

---

## 12. Explicit scope exclusions until the core is qualified

Do not use this campaign to simultaneously productionize:

- full SOC+GBT;
- projected SOC approximations;
- intersite Hubbard V under GBT;
- unrelated TD-DFT work;
- unrelated GPU/MPI optimization;
- large architectural refactors not required for the exact tests.

These can be separate follow-up projects after the scalar-relativistic GBT core is qualified.

---

## 13. External scientific references to use actively

The agent should consult primary sources as needed and record the exact equations/conventions used.

### Core GBT / LMTO / mature-code references

- ELK source and documentation for spin spirals/generalized Bloch theorem; inspect the actual `k +/- q/2` construction, not merely user documentation.
  - `https://elk.sourceforge.io/`
- LMTO spin-spiral literature by Shallcross, Nordström, Sharma and related work, especially formulations where the structure matrix is expressed through shifted `S(k +/- q/2)` while local potential functions remain local.
- The historical generalized Bloch/spin-spiral literature (Herring/Sandratskii and relevant descendants) for sign/convention cross-checks.

### Constraint/MFT references

- P. Bruno, *Exchange Interaction Parameters and Adiabatic Spin-Wave Spectra of Ferromagnets: A "Renormalized Magnetic Force Theorem"*, Phys. Rev. Lett. **90**, 087205 (2003). DOI: `10.1103/PhysRevLett.90.087205`.
- A. Jacobsson et al., *Efficient parameterisation of non-collinear energy landscapes in itinerant magnets*, Scientific Reports **12** (2022), article at `https://www.nature.com/articles/s41598-022-20311-7`.
- Relevant constrained-DFT work underlying the local constraining-field formalism, including Dederichs-type constrained magnetism.

### Internal design/evidence references

- `GBT_RS_LMTO_completion_blueprint.md`
- WP9 GBT report(s)
- VAL-16 commensurate-supercell report
- VAL-17 harmonic/Goldstone report

The previous internal documents are evidence/history, **not authority**. When they conflict with an exact derivation or primary external reference, document the conflict and follow the physics.

---

## 14. Definition of production-ready GBT for this campaign

The scalar-relativistic GBT implementation can be called production-ready only when all of the following are true:

- [ ] ordinary and GBT operators agree exactly at `q=0`;
- [ ] arbitrary global spin rotation at `q=0` preserves spectrum;
- [ ] one-sublattice pure-gauge GBT satisfies the exact `k +/- q/2` identity;
- [ ] bond reversal, loop gauge, and Hermiticity tests pass;
- [ ] HOH/CCOR and any enabled composite term separately pass covariance tests;
- [ ] matched commensurate GBT/supercell spectra agree to roundoff;
- [ ] rotating-frame and lab-frame density/moment transformations are explicit and tested;
- [ ] a genuinely nonzero constraining-field case converges and transforms correctly;
- [ ] physical constraint-energy bookkeeping is derived and documented;
- [ ] bare MFT, corrected constrained MFT, and fully constrained SCF are distinct reproducible workflows;
- [ ] the small-angle `sin^2(theta)` regime is demonstrated after k-grid convergence;
- [ ] the SOC-free dispersion is even in q within numerical tolerance;
- [ ] a stable small-q quadratic regime is demonstrated;
- [ ] GBT curvature/stiffness agrees with independently constructed LKAG `J(q)` within a predeclared converged tolerance;
- [ ] remaining SCF discrepancies, if any, have been either repaired or explicitly isolated outside the GBT Hamiltonian;
- [ ] all new tests run in the standard lean test workflow and documentation reflects actual support limits.

---

## 15. Campaign outcome philosophy

A negative result at an early gate is useful because it localizes the defect. Do not optimize for a visually attractive magnon dispersion; optimize for falsifiable equivalence.

The ideal campaign result is not merely “frozen magnons look reasonable.” It is the much stronger chain:

\[
\boxed{
q=0\text{ exact}
\;\land\;
(k\pm q/2)\text{ exact}
\;\land\;
\text{commensurate supercell exact}
\;\land\;
\text{frame/constraint covariance}
\;\land\;
\text{LKAG closure}
}
\]

Once this chain holds, any remaining failure of fully self-consistent spin spirals becomes a bounded SCF/reference-frame problem rather than an ambiguous “GBT problem.”

