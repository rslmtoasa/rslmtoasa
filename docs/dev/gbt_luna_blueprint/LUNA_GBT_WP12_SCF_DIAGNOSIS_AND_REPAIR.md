# Luna Prompt — GBT-WP12: SCF Diagnosis and Repair After Operator Validation

## Mission

Only after WP01-WP11 have passed or produced sharply localized residual issues, diagnose and repair any remaining failure of fully self-consistent GBT spin spirals.

At this stage the working assumption is that the GBT Hamiltonian/gauge has already been proven by exact fixed-potential tests. Therefore do **not** reopen the primitive `S x G` construction unless a new contradiction directly invalidates an earlier exact gate.

## Dependency

Hard dependency: WP01-WP11 evidence must be read first.

If a prerequisite gate is unresolved, stop and fix that gate instead of making speculative SCF changes.

## Primary hypothesis to test

The most likely remaining defect is an inconsistency between:

- rotating-frame potential/reference axes expected by `build_gbt_bulkham`;
- spinor density/moment returned by the reciprocal solver;
- magnetic mixing;
- radial spin-channel projection/update;
- constraining-field target/update.

Treat this as a frame and functional-consistency problem.

## Instrument one complete SCF iteration

For each sublattice record, in a controlled diagnostic build/run:

1. incoming rotating-frame reference direction;
2. incoming moment magnitude;
3. incoming radial spin-channel charges/potentials;
4. Hamiltonian local exchange splitting parameters;
5. constraint target and `B^C`;
6. output spin-density matrix from eigenstates;
7. output charge and `(mx,my,mz)` in the solver's frame;
8. longitudinal component relative to target;
9. transverse residual;
10. values passed to magnetic mixing;
11. values used to construct next radial spin-up/down channels;
12. next-iteration `potential%mom`;
13. physical total energy and constraint diagnostics.

The objective is to identify the first point at which a quantity leaves the declared WP05 frame contract.

## Define supported SCF policies explicitly

At minimum distinguish:

### `constrained_spiral`

- prescribed reference axes fixed in the rotating frame;
- charge and longitudinal moment magnitude may update;
- transverse density is a residual that drives the constraint field, not a new axis;
- radial spin channels are projected along the fixed target/reference axis.

### `relaxed_reference` — optional/future

- full rotating-frame moment direction may change;
- sublattice reference axes are updated while preserving the single-q ansatz.

Do not accidentally implement the second while believing the first is active.

If only constrained spirals are required for production now, explicitly guard or defer `relaxed_reference`.

## q=0 SCF closure first

Before finite q, require that GBT q=0 SCF reproduces ordinary SCF for:

- collinear state;
- global finite-angle state without SOC;
- multi-sublattice q=0 reference where applicable.

Compare converged:

- total charge;
- local charges;
- moment magnitudes;
- radial potential parameters;
- Fermi level;
- physical total energy;
- spectra.

## Finite-q constrained SCF

Then use a commensurate state with a validated explicit-supercell counterpart.

Compare GBT and explicit supercell:

- rotating/lab site moments;
- constraint fields;
- charge per equivalent site;
- physical total energy per primitive cell;
- selected spectral invariants.

Only after this passes move to incommensurate production spirals.

## Radial routines

Legacy atomic/radial LMTO routines are high-risk. Do not modify them unless instrumentation proves they receive inconsistent inputs or apply a collinear-axis assumption incompatible with the declared GBT policy.

Prefer correcting the projection/frame at the modern boundary rather than rewriting mature atomic physics.

## Mixing

Audit whether the current magnetic mixer handles:

- full Cartesian moments;
- longitudinal scalar magnitude;
- target direction;
- constraint residual;

as separate quantities.

For `constrained_spiral`, never let a tiny numerical transverse component rotate the reference axis unintentionally.

## Final qualification

After repair, rerun a compact regression subset from earlier work packages:

- q=0 exact operator test;
- shifted-k identity;
- commensurate matched operator;
- constraint fixture;
- harmonic small-q fixture;
- LKAG comparison.

This ensures the SCF repair did not contaminate fixed-potential GBT.

## Deliverables

- minimal SCF/frame repair if required;
- `docs/dev/GBT_REAUDIT_WP12_SCF_CLOSEOUT.md`;
- updated user/developer documentation describing supported GBT SCF semantics and exclusions;
- final production test matrix.

## Completion checklist

- [ ] All prerequisite WP evidence reviewed.
- [ ] One SCF iteration instrumented end to end.
- [ ] First frame/semantic divergence identified from evidence.
- [ ] `constrained_spiral` policy made explicit.
- [ ] q=0 GBT SCF matches ordinary SCF.
- [ ] finite-q commensurate GBT SCF compared with explicit supercell.
- [ ] constraint fields/moments compare in the correct frame.
- [ ] radial routines changed only if evidence requires it.
- [ ] physical energy semantics from WP07 preserved.
- [ ] earlier exact GBT gates rerun after repair.
- [ ] production support/exclusions documented.

## Commit

`fix: close GBT self-consistency frame contract`
