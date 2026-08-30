# Luna Prompt — GBT-WP08: Corrected Constrained Frozen-Magnon Mode

## Mission

Add a third clearly defined spin-spiral energy workflow between bare MFT and fully self-consistent GBT, following the corrected frozen-magnon logic used in the Jacobsson/ELK literature.

The three modes must remain distinguishable in code, documentation, and outputs.

## Dependency

WP07 must settle constraint energy/eigenvalue bookkeeping.

## Required modes

### 1. Bare MFT — existing concept

- converge a reference ordinary potential once;
- freeze the ordinary effective/XC potential;
- set/use no q-specific constraining field for the probe;
- rebuild only the q-dependent GBT Hamiltonian;
- evaluate the force-theorem energy according to the established code convention.

### 2. Corrected constrained MFT — new mode

Suggested input name: `mode='mft_constrained'` or `mode='corrected_mft'`.

For each q:

1. keep the ordinary reference potential fixed;
2. iterate only the constraining field until the integrated moment aligns with the prescribed spiral/cone target;
3. freeze the resulting q-dependent `B^C(q)`;
4. perform the final one-shot energy/eigenvalue evaluation with the same frozen ordinary potential plus `B^C(q)`;
5. apply the energy bookkeeping derived in WP07.

Do not allow charge/radial-potential updates to random-walk between q points.

### 3. Fully constrained SCF

- converge the ordinary potential and constraining field together at each q;
- maintain the same prescribed GBT reference directions;
- use the physical total energy from WP07.

## State reset requirements

Each q point must begin from a deliberate, documented state.

For corrected MFT:

- ordinary potential must be identical for all q;
- constraint field may be initialized from zero or a controlled continuation, but continuation must not change the converged result;
- no unintended density/potential mixing across q.

For SCF:

- if continuation is used for efficiency, provide a reproducibility check against independent starts for selected q points.

## Output requirements

The frozen-magnon output must state:

- mode;
- reference q;
- whether ordinary potential is frozen;
- whether constraint field was converged;
- constraint residual;
- `|B^C|` per sublattice;
- raw energy components;
- gauge-reference energy if same-q theta=0 subtraction is used;
- final physical `Delta E`;
- cone normalization used.

Keep machine-readable tabular output.

## Benchmarks

Run at least:

- bcc Fe;
- fcc Ni;
- B2 FeCo if multi-sublattice path is available.

The purpose is not to force agreement with a particular publication value but to reproduce the qualitative methodological distinction: Fe may show small correction, while Ni/FeCo should provide more discriminating constraint-sensitive cases.

## Important physics caveat

Do not label corrected constrained MFT automatically as the “true magnon spectrum.” Bruno-style renormalization/constrained response and adiabatic spin-wave response have subtleties. The code/documentation should call the methods by what they compute and defer final dynamical interpretation to the LKAG/response comparison.

## Tests

- reference potential hash/checksum or numerical equality across q in MFT modes;
- q-specific field convergence in corrected MFT;
- q restart/ordering independence;
- zero-field limit reduces corrected MFT to bare MFT;
- output clearly distinguishes modes.

## Deliverables

- new production frozen-magnon mode;
- tests;
- example inputs;
- `docs/dev/GBT_REAUDIT_WP08_CORRECTED_MFT.md`.

## Completion checklist

- [ ] Bare MFT semantics preserved and documented.
- [ ] Corrected constrained MFT implemented as a separate mode.
- [ ] Fully constrained SCF remains distinct.
- [ ] Ordinary potential proven frozen across corrected-MFT q points.
- [ ] q-specific constraint fields converge.
- [ ] Energy bookkeeping follows WP07 derivation.
- [ ] q ordering/restart independence tested.
- [ ] Fe benchmark run.
- [ ] Ni benchmark run.
- [ ] FeCo multi-sublattice benchmark run or explicitly deferred with reason.
- [ ] Outputs contain raw diagnostics and method metadata.

## Commit

`feat: add constrained frozen-magnon MFT mode`
