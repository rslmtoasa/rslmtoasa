# Luna Prompt — GBT-WP04: Matched Commensurate-Supercell Oracle

## Mission

Replace the historically ambiguous commensurate-supercell comparison with a decisive exact oracle in which the primitive GBT problem and explicit magnetic supercell are built from **the same finite directed primitive operator**.

This work package is not successful if it merely improves an energy difference. It must establish spectral equivalence after folding.

## Dependency

WP02 must pass. Prefer running after WP03 so the matched operator can later be repeated with validated composite terms.

## Historical problem to eliminate

Previous comparisons used primitive GBT and an explicit real-space supercell whose finite structure-constant/bond supports were not identical. Such a comparison cannot decide whether GBT is correct.

Do not repeat that architecture.

## Canonical construction

Create one canonical finite primitive directed bond table containing all information needed to build the tested Hamiltonian, e.g.:

- endpoint sublattices/types;
- cell displacement;
- physical displacement;
- orbital structure-constant block;
- optional Sdot for later extension;
- deterministic basis indexing.

Use this **same table** for both sides.

### Primitive GBT side

Apply the production endpoint link and construct `H_GBT(k,q)` in the chemical cell.

### Explicit supercell side

For a commensurate q with period N, create N translated copies and rotate the local magnetic frame explicitly:

\[
H^{SC}_{na,mb}
=
U_{na}\,\bar H_{ab}(R_m-R_n)\,U_{mb}^{\dagger}
\]

or the convention-equivalent expression derived from the code.

Do not regenerate neighbors independently. Do not alter cutoffs.

## Required q values

At minimum:

- period-2 case (`q` equivalent to 1/2 along a simple direction);
- period-3 case (`q` equivalent to 1/3);
- optionally period-4.

Use the repository's q units exactly and document the commensurability relation.

## Comparison

Compare complete eigenvalue multisets after the correct Brillouin-zone folding.

Preferred approach:

1. choose supercell k points;
2. generate the corresponding primitive reduced/folded k set;
3. collect all GBT eigenvalues belonging to the fold;
4. sort deterministically;
5. compare against explicit-supercell eigenvalues.

Acceptance: roundoff-scale in the exact finite-operator fixture.

Also compare invariants such as trace and Frobenius norm as diagnostics, but do not substitute them for spectrum equality.

## Basis/permutation diagnostics

If direct matrices are compared, explicitly construct the permutation/unitary mapping between primitive-folded and supercell bases. Keep this logic in test utilities rather than production physics if possible.

## Production integration follow-up

After the exact synthetic/matched test passes, add a higher-level periodic production fixture using real LMTO structure constants while ensuring the two representations still share equivalent support.

Do not use a separately truncated isolated/cluster real-space problem as the acceptance oracle.

## Deliverables

- reusable matched-operator test utility;
- automated period-2 and period-3 tests;
- `docs/dev/GBT_REAUDIT_WP04_COMMENSURATE.md` including folding derivation and residuals.

## Completion checklist

- [ ] One canonical primitive directed operator is used by both representations.
- [ ] Period-2 commensurate case passes.
- [ ] Period-3 commensurate case passes.
- [ ] Folding map is explicitly documented.
- [ ] Complete eigenvalue multisets compared.
- [ ] No empirical energy shift used.
- [ ] No independent neighbor/cutoff generation on the supercell side.
- [ ] Generic cone angle tested, not only theta=0.
- [ ] Multi-sublattice extension considered/documented.
- [ ] Higher-level production fixture added if practical.

## Commit

`test: add exact commensurate GBT supercell oracle`
