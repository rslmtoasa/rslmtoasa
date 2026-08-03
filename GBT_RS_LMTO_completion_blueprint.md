# GBT in RS-LMTO-ASA: audited completion blueprint

**Audit date:** 3 August 2026
**Audited baseline:** `fable_v2_gbt_v2` at `00abdf5`
**Primary input:** `GBT_RS_LMTO_pre_B1r_migration_notes.md`
**Target:** one correct nonrelativistic single-\(\mathbf q\) GBT operator used by both real-space recursion and the reciprocal-space solver, without damaging the general noncollinear real-space texture path.

This document is the implementation blueprint. The migration notes remain the detailed diagnosis and derivation.
Copy-ready agent tasks and the recommended delegation structure are in `GBT_RS_LMTO_task_prompts.md`.

---

## 1. Audit verdict

The migration note has the correct central conclusion:

\[
\{S(\mathbf d),\dot S(\mathbf d)\}
\rightarrow
\{S(\mathbf d)G_{ij},\dot S(\mathbf d)G_{ij}\}
\rightarrow H^{\rm GBT}(\mathbf d)
\rightarrow \{\text{recursion},\;H^{\rm GBT}(\mathbf k)\}.
\]

Here \(G_{ij}=U_i^\dagger U_j\) is the endpoint-frame link. GBT is
owned by the primitive directed structure-constant layer. Potential parameters
and the completed Hamiltonian are not separately cone- or q-rotated.

The current branch instead has three incompatible spiral mechanisms:

1. absolute-position moments in `ham0m_nc`;
2. the post-Fourier `h0/bz` reconstruction in `fourier_transform_gbt_array`;
3. the full-angle transverse rotation in the CCOR builders.

The note is therefore directionally correct, and its elementary pair formula

\[
K_{ij}=H_{\rm pair}[\mathbf m_i^0,R_z(\alpha_{ij})\mathbf m_j^0]D(\alpha_{ij})
\]

is algebraically sound for the existing two-centre `ham0m_nc` term, provided the structure matrix is spin independent and all endpoint quantities are represented in the same rotating frame.

However, the note is not yet a sufficient completion plan. Four items must be promoted from secondary concerns to architectural gates:

- k-space `EBAND` must be repaired independently of projected DOS;
- the SCF density and potential frame must be made explicit and common to both solvers;
- GBT must be separated from general explicit noncollinear textures by an actual mode/API boundary;
- multi-sublattice phase ownership must be made unambiguous so basis and sublattice phases cannot be counted twice.

The `fable_v2_gbt` branch is useful prototype evidence, not a completed solution. It contains a promising bond kernel and useful tests, but its own history records unresolved CCOR, multi-sublattice, k-mesh reuse, angle-dependence, and band-energy issues. Do not merge or cherry-pick it wholesale and declare GBT complete.

---

## 2. Corrections and additions to the migration note

### 2.1 Separate three magnetic representations

GBT and explicit noncollinear textures must not be branches of the same Boolean. Introduce an explicit representation mode, conceptually:

| Mode | Endpoint moments | Bond gauge | Storage allowed | Intended use |
|---|---|---|---|---|
| `periodic_nc` | actual chemical-cell/per-type moments | none | per-type `ee`, per-site `hall` | ordinary collinear/noncollinear states |
| `gbt_single_q` | collinear rotating-frame potential, normally \(+\hat z\) for a ferromagnet | endpoint link on primitive `S`/`Sdot` | per-type `ee`; optional per-site oracle | bulk single-\(q\) spirals/cones |
| `explicit_texture` | actual moment of physical site \(i\) and \(j\) | none | per-site `hall` or a true magnetic supercell | skyrmions, arbitrary textures, commensurate oracle |

The existing ordinary Pauli algebra in `ham0m_nc` is shared by
`periodic_nc` and `explicit_texture`. It remains q-agnostic and retains the
general endpoint-moment rotation logic. `gbt_single_q` instead assembles from
the fixed collinear rotating-frame potential and a gauged primitive structure
block; it must not feed cone/q moments through `ham0m_nc` and then rotate the
result again.

For a collinear ferrimagnetic or antiferromagnetic reference, preserve each
sublattice's local spin-channel ordering/reference sign. “Collinear rotating
frame” does not mean replacing distinct exchange splittings by one identical
ferromagnetic potential.

This protects the required skyrmion/general-texture capability. The absolute-position GBT branch must be deleted, but a clearly named explicit-texture path must remain. That path must obtain moments by **site identity**, not merely by chemical type; otherwise two sites of the same type cannot carry different directions.

`gbt_kspace` should not survive as a physics selector. During migration it may be accepted as a deprecated input alias, but the final Hamiltonian must not inspect it. Solver choice belongs to `self%use_kspace`; magnetic representation belongs to the new texture mode.

### 2.2 Choose one multi-sublattice gauge convention

The migration note presents equivalent conventions, but its later wording can be read as putting a sublattice azimuth both into \(\chi_a\) and into \(\mathbf m_a^0\). That would double count it.

Use the physical-displacement convention throughout the implementation:

\[
\alpha_{ij}=\mathbf q_{\rm cart}\cdot\mathbf d_{ij,\rm cart},
\qquad
\mathbf d_{ij}=\mathbf R+\boldsymbol\tau_b-\boldsymbol\tau_a.
\]

Store any additional, non-geometric sublattice orientation in an endpoint
reference frame \(U_a^0\). With \(U_a^0\hat z=\mathbf m_a^0\), define the
directed link

\[
G_{ab}(\mathbf d)=
(U_a^0)^\dagger R_z(\alpha_{ij})U_b^0,
\]

and **do not** add \(\chi_b-\chi_a\) to `alpha`. An alternative
cell-translation convention is valid, but it must not coexist in production
code. For a one-sublattice common cone at q=0, \(G=I\): the rotating-frame
Hamiltonian is literally the collinear Hamiltonian, while laboratory-frame
moments are reconstructed only for output.

This convention also makes the desired decomposition clear:

- geometry supplies the complete directed physical bond;
- `q_ss` supplies the translational twist;
- per-sublattice input supplies \(U_a^0\), the reference orientation;
- no downstream routine invents another phase.

### 2.3 Gauge primitive `S`/`Sdot`, not potential parameters or completed Hamiltonians

For spin-independent structure constants, transform the earliest directed
factor:

\[
S_{ab}(\mathbf d)\mathbf1_2
\rightarrow S_{ab}(\mathbf d)G_{ab}(\mathbf d),
\qquad
\dot S_{ab}(\mathbf d)\mathbf1_2
\rightarrow \dot S_{ab}(\mathbf d)G_{ab}(\mathbf d).
\]

The LMTO pair is then assembled with the unchanged collinear rotating-frame
potential parameters, schematically

\[
K_{ab}^{\rm GBT}=W_a^{z}\,[S_{ab}G_{ab}]\,W_b^{z}.
\]

No q/cone rotation is applied to `wx0/wx1`, `potential%mom`, or a completed
`ee`/Hamiltonian afterward. It must not be applied blindly to every later
object. The following need separate classification:

- on-site `enim`, `obarm`, constraining fields, Hubbard-U and other local terms: no translational phase and remain in the collinear rotating frame;
- intersite Hubbard-V, CCOR, derivative/velocity matrices, overlap, and any downfolded term: twist the earliest directed bond factor, then perform contractions;
- products such as HOH: prove covariance at the operator level and compare the real-space two-sweep result with the reciprocal-space product;
- SOC: reject nonzero-\(q\) GBT rather than attempting this gauge.

### 2.4 `EBAND` is a prerequisite, not a final validation detail

The current reciprocal MFT route obtains band energy through
`reciprocal%calculate_band_energy_from_moments`, which sums `m0*m1` from spin/orbital projected DOS. This makes a tiny frozen-magnon energy difference depend on:

- the DOS energy window and grid;
- Gaussian/tetrahedron broadening and normalization;
- projection completeness;
- the selected local spin axis;
- Fermi-level search details.

That path is unsuitable as the canonical force-theorem energy. A projection-free occupied-eigenvalue sum must be implemented before judging k-space GBT:

\[
E_{\rm band}=\sum_{n\mathbf k}w_{\mathbf k}
f(\epsilon_{n\mathbf k}-E_F)\epsilon_{n\mathbf k}.
\]

For a tetrahedron calculation, use consistent tetrahedron occupation and energy weights. For the initial correctness phase, a full uniform mesh with Fermi-Dirac occupations is an acceptable canonical oracle. Projected-DOS and integrated-total-DOS energies should remain diagnostics only.

### 2.5 SCF needs one density-frame contract

The current real-space route reconstructs spin information from the full Green function, while the k-space route separately constructs projected spin channels and a spin vector from eigenvectors. In k-space, `ql(up/down)` is projected along `potential%mom`, while the vector moment is updated afterward. This is not an explicit rotating-frame density contract and can feed inconsistent spin channels into the radial SCF.

Both solvers must produce the same per-site/per-\(l\) rotating-frame density object:

\[
\rho_{a l}=
\begin{pmatrix}
\rho_{\uparrow\uparrow}&\rho_{\uparrow\downarrow}\\
\rho_{\downarrow\uparrow}&\rho_{\downarrow\downarrow}
\end{pmatrix},
\]

or equivalently \((n,m_x,m_y,m_z)\). Only after this common accumulation may the code project onto the local radial up/down axis.

For a constrained spin spiral, that axis is the imposed rotating-frame reference direction; mix charge and longitudinal moment magnitude, and report transverse density as a torque/constraint residual. For an explicitly relaxed periodic noncollinear state, mix the full Cartesian moment and update the reference axis. Do not silently alternate between those meanings.

---

## 3. Non-negotiable operator contract

The following is the normative implementation contract.

### 3.1 Directed bond identity

Every stored bond has a key

```text
(row basis/type a, column basis/type b, lattice translation R,
 complete Cartesian displacement d, neighbour slot)
```

and a reverse key \((b,a,-\mathbf R,-\mathbf d)\). In debug builds, assert that the neighbor-list displacement used by `chbar_nc` agrees with the displacement used by the ordinary Fourier transform.

### 3.2 Phase helper

Only one routine computes a GBT phase:

```fortran
alpha = gbt_bond_phase(q_cart_2pi_over_alat, d_cart, alat)
```

with

\[
\alpha=2\pi\,\mathbf q_{\rm stored}\cdot\mathbf d_{\rm cart}/a_{\rm lat}.
\]

The helper must provide a debug cross-check against reciprocal/direct coordinates modulo \(2\pi\). No other production routine may contain `dot_product(q_ss, ...)`.

### 3.3 Ordinary noncollinear pair pipeline

Keep `hmfind -> chbar_nc -> ham0m_nc` as the non-GBT production pipeline.
Its existing `wx0/wx1`, dot-product, cross-product, on-site, and Pauli algebra
must remain q-agnostic and retain its general moment handling. It must not be
widened merely to carry the complex linked GBT structure block.

If later explicit-texture work benefits from extracting the algebra, use a
behavior-preserving routine with explicit inputs:

```fortran
assemble_nc_pair(type_i, type_j, moment_i, moment_j, structure, pair_pauli)
```

Such extraction is not a prerequisite for S-level GBT and may not change the
existing call-chain contract.

### 3.4 GBT structure-link kernel

The GBT kernel operates before Hamiltonian contraction and performs exactly:

1. resolve the endpoint reference frames \(U_a^0,U_b^0\);
2. compute \(\alpha\) once from the directed physical bond;
3. construct \(G_{ab}=(U_a^0)^\dagger R_z(\alpha)U_b^0\);
4. lift the primitive orbital `S` block, and `Sdot` when requested, to a
   spinor block by right multiplication with \(G_{ab}\);
5. contract those gauged primitive blocks with unchanged collinear
   rotating-frame potential parameters.

Construct the link as an explicit dense \(2\times2\) matrix (or a separately
verified Pauli equivalent) using temporaries:

\[
R_z(\alpha)=
\begin{pmatrix}
e^{-i\alpha/2}&0\\
0&e^{+i\alpha/2}
\end{pmatrix},
\qquad
G_{ab}=(U_a^0)^\dagger R_z(\alpha)U_b^0.
\]

Do not reuse the old in-place completed-Hamiltonian Pauli update as the
production implementation of this structure link.

Neither `ham0m_nc` nor reciprocal Fourier code may apply an additional GBT
rotation to the result. Raw `lattice%sbar`/`sdot` should remain available for
ordinary modes and independent comparison; build a representation-specific
spinor view/block rather than silently corrupting shared lattice state.

The production S-level GBT path is required only for
`strux_backend='strux_lib'`. The legacy structure backend remains available
for non-GBT calculations but must fail early if `gbt_single_q` is requested;
it does not provide the derivative-structure contract used by the audited
path. `strux_want_sdot` is required only when an enabled term consumes
`Sdot`; the first-order `S` path must not request it unnecessarily. Existing
`strux_lib` restrictions on Sdot screening modes remain in force.

Do **not** double the persistent `lattice%sbar` or `lattice%sdot` layout. They
remain complex orbital blocks with their current dimensions. Prefer a parallel
representation-specific workflow:

```text
periodic_nc / explicit_texture:
    hmfind -> chbar_nc -> ham0m_nc

gbt_single_q:
    fetch raw strux-lib orbital S/Sdot + directed-bond metadata
        -> build temporary nb x nb spinor structure block with G_ab
        -> collinear rotating-frame LMTO assembler
        -> ee / audited companion block
```

The existing `hmfind`, `chbar_nc`, and `ham0m_nc` interfaces remain the
non-GBT path and should not be widened merely to carry complex spinor
structures. The new GBT-collinear assembler consumes the linked `nb x nb`
block and left/right diagonal collinear potential factors. At q=0 it must
reproduce the existing `ham0m_nc` result for `mom=(0,0,1)` within roundoff,
including unequal endpoint potential parameters.

A GBT spinor-structure cache is optional. If added, it is keyed by backend,
directed bond, q, cone, sublattice frames, and structure generation, and must
be invalidated accordingly. A per-bond temporary is the preferred first
implementation.

### 3.5 Explicit-texture wrapper

The explicit-texture wrapper obtains `moment_i` and `moment_j` from a
site-indexed texture provider and uses the ordinary NC pair algebra. It applies
no GBT link or \(D\) phase.

Bulk per-type storage must reject this mode unless the input is represented as a true magnetic supercell in which each inequivalent moment direction has a distinct periodic site/type mapping.

### 3.6 Solver consumption

Both solvers consume the same completed operator:

```text
raw S/Sdot -> representation-specific structure link -> LMTO contraction
           -> ee/eeo/... -> recursion
                         -> ordinary Fourier transform -> eigensolver
```

There is no GBT-specific reciprocal transform.

---

## 4. Ordered implementation work packages

Each work package ends in a gate. Do not start physical bcc-Fe interpretation until gates 0 through 5 pass.

WP2 exposed a dependency cycle in the original gate wording: canonical EBAND
can be validated independently, but its production q=0 GBT check cannot be
meaningful while the legacy double-representation operator remains. Therefore
G2 is tracked as two evidence slices:

- **G2E (energy):** canonical occupation/EBAND, MPI, normalization, DOS
  independence, and total-DOS oracle. This is complete and passing.
- **G2O (operator):** production q=0 rotating-frame reduction using the
  corrected representation, plus a finite-q or relative-sublattice check that
  proves the probe is not disabled. This closes after WP4.

Only WP3 and WP4 are permitted while final G2 is open; WP5 and later work still
require final G2 plus their preceding gates.

### WP0 — Freeze the baseline and block known-wrong combinations

1. Preserve representative inputs and outputs for:
   - ordinary collinear bulk;
   - periodic noncollinear bulk;
   - explicit real-space two-site/texture case;
   - current bcc-Fe real-space MFT;
   - current k-space MFT/SCF, labelled as known-wrong rather than golden physics.
2. Add temporary fatal guards for nonzero-\(q\) with SOC, CCOR, Hubbard-V, an unverified overlap mode, or any other term not yet classified.
3. Force the full chemical BZ for nonzero \(\mathbf q\).
4. Add a diagnostic before `zheev/zhegv` for \(\|H-H^\dagger\|\) and \(\|O-O^\dagger\|\); fail rather than letting the upper-triangle eigensolver hide an assembly defect.

**Gate G0:** all ordinary \(q=0,\theta=0\) regression cases are unchanged, known-wrong GBT combinations fail explicitly, and the current old k-space GBT is demonstrably non-Hermitian on a nonzero-\(q\) test.

**WP0 completion checklist (completed 3 August 2026; evidence in
`docs/dev/GBT_WP0_G0_REPORT.md`):**

- [x] Baseline inputs/outputs recorded and classified.
- [x] Unsupported combinations fail early.
- [x] Full BZ is forced for nonzero q.
- [x] H/O Hermiticity checks run before eigensolution.
- [x] Ordinary regressions are unchanged (the sole gfortran-13 DOS tolerance
  delta is reproduced by the pre-existing gfortran-13 binary).
- [x] G0 report includes commands, results, and remaining risks.

### WP1 — Establish independent algebraic oracles

1. Port/review the useful tests and dense oracle from commits `2072fac` and `7ca0c29` without yet changing production routing.
2. Add:
   - phase-unit tests for cubic, hexagonal, and triclinic cells;
   - one-orbital Stoner shifted-\(k\) tests;
   - random dense endpoint-species tests;
   - directed-bond reverse-pair tests;
   - multi-sublattice tests that would fail if a phase is double counted.
3. Ensure one side of each comparison uses direct dense \(2\times2\) spin matrices rather than the production Pauli helper.

**Gate G1:** algebraic relative error \(<10^{-12}\) in double precision and exact \(q=0\) reduction for the no-cone path.

### WP2 — Repair reciprocal electron count, Fermi level, and `EBAND`

1. Add a projection-free electron-count routine from eigenvalues and k weights.
2. Add a projection-free occupied-eigenvalue band-energy routine using the same occupations and \(E_F\).
3. Implement MPI reduction once, with an explicit k-weight normalization contract.
4. For tetrahedron/Blöchl mode, either implement consistent band-energy weights or temporarily route GBT validation through the full-mesh Fermi-Dirac oracle.
5. Make `frozen_magnon_probe_energy` use only this canonical routine.
6. Keep the following as cross-checks, not sources of truth:
   - integral of \(E D_{\rm total}(E)f(E)\);
   - sum of projected `m0*m1`.
7. Assert electron-number conservation and log the three energy estimates.

**Gate G2:** on an ordinary magnetic k-space case, direct `EBAND` is stable
against DOS grid/window changes, agrees with a converged total-DOS integral
within its integration error, and follows the exact q=0 rotating-frame
reduction. A finite-q or relative-sublattice operator assertion must also prove
that q=0 equality is not obtained by disabling the probe. The projected-DOS
route may not be used for MFT even if it happens to agree. Record G2E and G2O
separately until both pass.

### WP3 — Introduce the representation split and preserve general NC

1. Add the three magnetic representation modes from section 2.1.
2. Preserve the existing `hmfind -> chbar_nc -> ham0m_nc` pipeline for
   ordinary NC/texture use. Extraction of a pure NC pair kernel is allowed only
   if it leaves that interface and behavior unchanged; it is not required by
   the GBT path.
3. Add a site-indexed moment provider for `explicit_texture`.
4. Route `build_bulkham` and `build_locham` explicitly:
   - `periodic_nc`: existing periodic endpoint logic;
   - `gbt_single_q`: collinear rotating-frame potential/reference blocks;
   - `explicit_texture`: per-site build only.
5. Add a per-sublattice reference-frame provider for GBT without yet applying
   q/cone phases to potential parameters or completed Hamiltonians.
6. Add an explicit failure if arbitrary site moments are requested through per-type bulk reuse.
7. Retain `rotate_to_local_axis`/`rotate_from_local_axis`, the `wx0` bounds
   checks, and existing general noncollinear tests.

**Gate G3:** all existing noncollinear tests pass; a two-site texture with same chemical type but different site moments produces distinct correct endpoint blocks; a small skyrmion-like texture smoke test is unchanged by enabling the refactor with GBT inactive.

### WP4 — Implement GBT at primitive `S`/`Sdot` level

1. Add the central phase helper; prohibit q/cone calculations in `ham0m_nc`.
2. Construct the endpoint link
   \(G_{ab}=(U_a^0)^\dagger R_z(\alpha)U_b^0\) from the complete directed
   physical displacement.
3. Apply the link to primitive `S` and, when enabled, `Sdot` before LMTO
   contractions. Keep raw lattice structure constants unchanged.
4. Assemble GBT with collinear rotating-frame potential parameters. Do not
   rotate `potential%mom`, `wx0/wx1`, `enim`, `obarm`, `ee`, or the completed
   Hamiltonian afterward.
5. Build `eeo`/HOH from gauged primitive inputs only where covariance is
   proved; otherwise keep GBT+HOH guarded for WP6. Keep GBT+CCOR guarded until
   its independent `Sdot`/composite audit removes the remaining angle logic.
6. Add raw/gauged structure-bond dumps keyed by directed geometry and compare
   them with the dense oracle.
7. Make q=0 common-cone GBT bitwise/roundoff equivalent to the collinear
   rotating-frame operator, and add a finite-q or relative-sublattice assertion
   that cannot pass if the probe is a no-op.
8. Bypass all new floating-point work when GBT is inactive.
9. Require `strux_backend='strux_lib'` for `gbt_single_q`; reject the legacy
   backend before Hamiltonian construction. Require `strux_want_sdot` only for
   audited Sdot consumers.
10. Add a parallel complex GBT-collinear block path rather than doubling SBAR
    or widening the existing `hmfind/chbar_nc/ham0m_nc` interfaces.

**Gate G4:** dense structure-link oracle, Fourier shifted-\(k\) identity,
directed-bond Hermiticity, \(q\leftrightarrow-\mathbf q\), \(\theta=0\), q=0
rotating-frame reduction, and a non-no-op finite-q/relative-sublattice test all
pass. Rerun G2O here; final G2 requires both G2E and G2O.

### WP5 — Delete the reciprocal reconstruction and unify solvers

1. Remove `fourier_transform_gbt_array` and `fourier_transform_gbt` implementations, interfaces, and call sites; they may not reconstruct `ee`, `eeo`, or `eecc`.
2. Fourier-transform `ee` and verified companion arrays only through `fourier_transform_array`.
3. Remove every `gbt_kspace` conditional from reciprocal assembly.
4. Invalidate/rebuild all q-dependent cached Hamiltonian, eigenvalue, eigenvector, DOS, and density objects when \(\mathbf q\), cone, sublattice reference axes, or potential parameters change.
5. Compare real-space and k-space results from the same bonds.

**Gate G5:** the same bond dump is consumed by both routes; k-space matrices are Hermitian before diagonalization; solver differences converge monotonically with recursion depth, broadening, and k mesh.

### WP6 — Audit every Hamiltonian and overlap term

Create a checked feature matrix. No blank or “probably okay” entries are allowed.

| Term | Required action |
|---|---|
| elementary `ee` | build from GBT-linked `S` in WP4 |
| on-site `cx/cex`, `enim`, `obarm` | collinear rotating frame; no translational phase |
| `eeo` / HOH | build from gauged primitive bonds; prove RS two-sweep = reciprocal product |
| `eeoee` | verify whether it is diagnostic, local approximation, or consumed operator; never phase a precontracted object without proof |
| overlap/generalized solve | construct covariantly; require Hermiticity and positive definiteness |
| CCOR | remove remaining sublattice/full-angle rotations; gauge raw `S`/`Sdot` at their earliest directed use |
| Hubbard-U | rotate/localize in the same on-site frame; no translational phase |
| Hubbard-V | gauge the directed intersite term or reject it in GBT mode |
| constraining fields | express in the rotating reference frame and include consistent energy bookkeeping |
| velocity/torque operators | derive from the already-gauged bond operator; test finite differences in k |
| SOC `lsham` | fatal for nonzero-\(q\) GBT |

**Gate G6:** every enabled feature has an operator-level oracle and passes directed-bond and k-space Hermiticity tests. Unsupported combinations abort before SCF begins.

### WP7 — Unify density accumulation and make SCF well defined

1. Define one rotating-frame density data structure, per site and angular channel, retaining the complete spin matrix.
2. Implement two producers:
   - real-space Green-function accumulation;
   - k-space eigenvector/k-weight/occupation accumulation.
3. Test the producers on the same finite periodic model before connecting radial SCF.
4. Define the SCF policy explicitly:
   - `constrained_spiral`: reference directions fixed; mix charge and longitudinal magnitude; transverse density is a constraint residual;
   - `relaxed_reference`: mix full Cartesian rotating-frame moments and update per-sublattice reference axes while preserving the single-\(q\) ansatz.
5. Derive radial spin channels from density in that chosen axis; never use a stale `potential%mom` as an implicit projection definition.
6. Reconstruct lab-frame moments only for output or explicit-supercell comparison:

   \[
   \rho_{na}^{\rm lab}=D(\mathbf q\cdot\mathbf r_{na})\bar\rho_aD^\dagger(\mathbf q\cdot\mathbf r_{na}).
   \]

7. Mix the same variables, in the same frame, on both solvers.
8. Check density Hermiticity, positive semidefiniteness, trace/electron count, and \(|\mathbf m|\le n\) every iteration in debug mode.

**Gate G7:** real-space and k-space SCF converge to the same charge, rotating-frame moment, Fermi energy, band energy, and total energy within independently converged numerical tolerances. A \(q=0\) cone is globally rotation invariant without SOC.

### WP8 — Symmetry and q-sweep lifecycle

1. Keep full-BZ results as the permanent oracle.
2. Add the little group only after G7.
3. Key any reduced-mesh cache by lattice, mesh, offset, \(\mathbf q\), and relevant anti/unitary policy.
4. In a multi-q sweep, either rebuild the little-group mesh for every q or use the intersection of the little groups of all q points.
5. Verify that permuting q-list order leaves every energy unchanged.

**Gate G8:** reduced and full BZ agree for axial and generic q vectors; a sweep never reuses row 1's little group for a different q.

### WP9 — Physical validation ladder

Run in this order:

1. analytic one-orbital Stoner model;
2. random dense LMTO-like bond oracle;
3. commensurate explicit supercell vs primitive-cell GBT, first MFT then SCF;
4. real-space vs k-space DOS/charge/moment/`EBAND`;
5. bcc-Fe \(\Gamma\)-H MFT dispersion with full BZ and converged `EBAND`;
6. \(E(\mathbf q)=E(-\mathbf q)\) and smoothness around H;
7. small-cone scaling after division by \(\sin^2\theta\), using a conservative harmonic window such as 5–20 degrees;
8. spin stiffness vs independently converged LKAG \(J(\mathbf q)\), including shell-convergence error bars;
9. constrained SCF spirals vs explicit commensurate supercells;
10. multi-sublattice tests only after the single-sublattice gates pass.

For algebraic tests use \(10^{-12}\)-level relative tolerances. For production comparisons, first converge each method until its numerical error is comfortably below the target spiral energy difference, then require agreement within the larger of a declared absolute budget and 1% of \(\Delta E\). A regression pin from unconverged settings is not a physics acceptance test.

**Gate G9:** all Hamiltonian, energy, density, solver-agreement, supercell, and physical invariance criteria in section 6 pass.

### WP10 — Delete legacy GBT and close documentation

Only after G9:

1. delete the absolute-position spiral branch from `ham0m_nc`;
2. delete `fourier_transform_gbt*` completely;
3. delete both CCOR full-angle GBT rotations;
4. delete or complete deprecation of `gbt_kspace`;
5. remove GBT use of cone/sublattice branches in `ham0m_nc`, while retaining
   explicitly routed non-GBT NC/texture moment inputs;
6. remove obsolete comments claiming that recursion carries an explicit bulk spiral;
7. remove or relabel old GBT golden data that encode known-wrong physics;
8. document `explicit_texture` separately from GBT and add a skyrmion-oriented example;
9. make a repository-wide check fail if spiral phases are calculated outside the phase helper.

**Gate G10:** repository search finds one phase helper, one GBT structure-link
helper, no reciprocal GBT reconstruction, no CCOR full-angle spiral rotation,
and a fully tested explicit-texture path.

---

## 5. Recommended code ownership

The exact filenames may evolve, but responsibilities should be separated as follows:

| Responsibility | Recommended home |
|---|---|
| representation mode, q/cone/reference inputs | `hamiltonian.f90` or a small `spin_texture` module |
| phase units and directed-bond phase | one `gbt_bond_phase` helper |
| ordinary NC Pauli pair algebra | existing `ham0m_nc`, optionally a behavior-preserving extracted helper |
| endpoint-frame/link algebra | small pure GBT structure-link module |
| raw orbital `S`/`Sdot` accessor | strux-lib-aware structure-block API; legacy accessor retained for non-GBT |
| GBT lift of raw `S`/`Sdot` | parallel temporary `nb x nb` structure-block builder before LMTO contraction |
| GBT collinear LMTO contraction | new assembler parallel to, not embedded in, `ham0m_nc` |
| explicit site-texture provider | separate site-indexed texture object/module |
| ordinary Bloch sum | `reciprocal_fourier.f90` |
| canonical electron count and `EBAND` | reciprocal occupation/energy module, independent of projections |
| common density contract | shared density accumulator interface with RS and k-space backends |
| physical validation | `tests/unit`, `tests/regression/gbt_supercell`, and dedicated converged post-processing tests |

Avoid putting phase policy in `calculation.f90`, DOS projections, the eigensolver, or CCOR-specific code.

---

## 6. Final definition of done

GBT is complete only when all of the following are true.

### Architecture

- One production GBT structure-link transformation exists at primitive
  `S`/`Sdot` level.
- GBT potential parameters and completed Hamiltonians are never independently
  cone/q rotated.
- Recursion and reciprocal space consume the same real-space operator.
- Solver selection does not change the Hamiltonian physics.
- General periodic NC and explicit site-texture modes remain independently usable.
- Arbitrary real-space neighboring moment directions still reach the unchanged Pauli pair algebra.

### Algebra and matrices

- No absolute site phase enters a per-type GBT block.
- A one-sublattice common cone at q=0 reduces to the identical collinear
  rotating-frame operator, while lab-frame output reconstructs the cone.
- Phase units and basis convention are centralized and tested in non-cubic cells.
- Reverse directed bonds are Hermitian conjugates.
- \(H(\mathbf k)\) and any \(O(\mathbf k)\) are checked before eigensolution.
- Every enabled correction has a proved gauge path.

### Energy

- k-space `EBAND` is computed directly from eigenvalues/occupations or equivalent exact tetrahedron weights, never from projected DOS.
- Electron number and k-weight normalization are asserted.
- MFT energy differences are stable against DOS output settings.
- q-list ordering and symmetry settings do not change converged results.

### SCF

- Both solvers retain and mix the same rotating-frame density information.
- Radial spin channels, local reference axes, and Hamiltonian endpoint axes are consistent.
- Constrained and relaxed spiral policies are explicit.
- Real-space and k-space SCF agree after independent numerical convergence.

### Physical validation

- inactive GBT reproduces ordinary calculations;
- \(q=0\) cone rotations are isospectral without SOC;
- \(\theta=0\) is q independent;
- \(E(\mathbf q)=E(-\mathbf q)\) where required;
- commensurate GBT agrees with an explicit supercell in MFT and SCF;
- bcc-Fe stiffness agrees with an independently converged exchange calculation within declared uncertainty;
- multi-sublattice acoustic modes satisfy the Goldstone condition without hiding a defect through an imposed post hoc shift.

### Cleanup

- `fourier_transform_gbt*` is absent;
- old absolute-position bulk GBT is absent;
- old CCOR full-angle GBT is absent;
- `gbt_kspace` has no physics role and is removed after deprecation;
- the explicit-texture/skyrmion path is retained, named, documented, and tested.

---

## 7. Practical branch strategy

Use short, gate-oriented changes rather than one migration merge:

```text
baseline guards/oracles
    -> canonical k-space EBAND
    -> magnetic representation split
    -> primitive S/Sdot GBT link
    -> close G2 operator slice
    -> delete reciprocal reconstruction
    -> term-by-term HOH/CCOR/overlap audit
    -> common density + SCF
    -> symmetry
    -> physical validation
    -> legacy deletion/documentation
```

The useful `fable_v2_gbt` commits can be mined selectively:

- `2072fac`: oracle/test ideas;
- `7ca0c29`: elementary bond-kernel prototype;
- `17ce417`: reciprocal-transform deletion;
- `440c594`: phase-convention tests;
- `3fd21c0`: commensurate-supercell test structure;
- `b1caf4b`: little-group machinery, only after fixing per-q cache ownership;
- `92d2471` and `091c2c1`: k-space moment/frame lessons.

Before reuse, adapt them to the representation split and the canonical density/energy contracts above. The unresolved findings recorded after `c44075e` are reasons to retain the gates, not reasons to weaken them.
