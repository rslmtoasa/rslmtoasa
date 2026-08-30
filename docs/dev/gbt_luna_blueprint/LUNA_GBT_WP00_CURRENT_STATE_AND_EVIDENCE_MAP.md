# Luna Prompt — GBT-WP00: Current-State Call Graph and Evidence Ledger

## Mission

Audit the current `fable_v3` GBT implementation **without changing physics**. Produce a precise call graph, representation map, and evidence ledger that establishes what code is actually active today and what historical GBT blueprint items have landed.

This is archaeology and scoping only. Do not “fix” suspicious code in this work package unless required to make the audit tooling itself compile.

## Repository and baseline

Work in:

`https://github.com/rslmtoasa/rslmtoasa`

branch:

`fable_v3`

Read first:

- `00_HEAD_GBT_REVALIDATION_MASTER_BLUEPRINT.md` from this package;
- repository `GBT_RS_LMTO_completion_blueprint.md`;
- current GBT-related source;
- archived WP9 / VAL-16 / VAL-17 GBT reports.

Also inspect external formulations enough to understand the intended architecture:

- ELK spin-spiral `k +/- q/2` implementation;
- LMTO spin-spiral literature with shifted structure matrices;
- Bruno/Jacobsson only to identify where constraints enter the call graph.

## Questions that must be answered

1. What exact input variables activate GBT?
2. What units and coordinate basis are used for `q_ss`?
3. How are cone angle and sublattice phase/reference angles represented?
4. Where is the magnetic representation dispatched?
5. What is the complete call chain from input through:
   - structure constants;
   - GBT link construction;
   - LMTO hopping construction;
   - HOH;
   - CCOR/Sdot;
   - Fourier transform;
   - diagonalization/recursion;
   - density/moment reconstruction;
   - SCF update;
   - constraining-field update;
   - frozen-magnon post-processing?
6. Which arrays are in orbital space, spinor space, spherical Pauli representation, rotating frame, and lab frame?
7. Where is `potential%mom` consumed in GBT and what components matter?
8. Where is `mag_cfield` stored, updated, and inserted?
9. What exactly is `constraint_energy` and where does it enter reported total energy?
10. Which features are rejected under GBT and which are silently permitted?
11. Which old blueprint tasks are:
    - fully implemented;
    - partially implemented;
    - implemented differently;
    - obsolete;
    - still missing?
12. Which archived validation results remain logically valid, and which were based on unmatched operators, stale potentials, or weak acceptance criteria?

## Required deliverables

Create `docs/dev/GBT_REAUDIT_WP00_CURRENT_STATE.md` containing:

### A. Call graph

Use a readable nested graph beginning with the namelist/input and ending at observables/SCF. Include source file and routine names.

### B. Data/frame table

For every important object, document:

| Object | Source | Shape/space | Frame | Meaning | Mutability |
|---|---|---|---|---|---|

Include at least `sbar`, `Sdot`, GBT link, `ee`, `eeo`, `enim`, `obarm`, `potential%mom`, density/moment outputs, constraint reference, and `mag_cfield`.

### C. Feature support matrix

Explicitly list SOC, HOH, CCOR, local-axis, Hubbard U, Hubbard V, real-space recursion, k-space, SCF, MFT, constrained SCF, and multi-sublattice frozen magnons.

Do not infer support from comments alone. Trace executable code.

### D. Blueprint landing matrix

Map every material action from the old completion blueprint onto present code and classify it.

### E. Validation evidence ledger

For each historical test/report, state:

- what it truly tested;
- what it did not test;
- whether the fixture used a matched operator;
- whether it was fixed-potential or SCF;
- whether constraints were genuinely nonzero;
- whether its conclusion remains admissible.

### F. Risk register

Rank suspected areas, but do not change code yet. Distinguish evidence from hypothesis.

## Required source checks

At minimum inspect symbols equivalent to:

- `gbt_endpoint_link`;
- `gbt_lift_orbital_block`;
- `gbt_contract_collinear`;
- `build_gbt_bulkham`;
- magnetic-representation normalization/dispatch;
- ordinary bulk Hamiltonian construction;
- Fourier transform path used by GBT;
- reciprocal density/moment producer;
- `apply_constraints`;
- `constrain`;
- GBT onsite constraining-field insertion;
- frozen-magnon acoustic and multi-sublattice paths.

## Guardrails

- Do not refactor core GBT code.
- Do not change phase conventions.
- Do not change numerical results.
- Do not promote/demote feature support based on intuition.
- If comments contradict executable behavior, executable behavior wins and the contradiction is documented.

## Completion checklist

- [ ] Current branch/commit recorded.
- [ ] Full GBT call graph produced.
- [ ] Ordinary versus GBT construction paths contrasted.
- [ ] All relevant frame conventions recorded.
- [ ] `potential%mom` semantics traced.
- [ ] `mag_cfield` semantics traced.
- [ ] Constraint-energy bookkeeping traced.
- [ ] HOH and CCOR construction traced.
- [ ] Frozen-magnon modes traced.
- [ ] Historical blueprint landing matrix complete.
- [ ] WP9/VAL-16/VAL-17 evidence reclassified.
- [ ] Risk register written with evidence/hypothesis separation.
- [ ] No physics changes introduced.

## Commit

When complete, make one focused commit with the one-line message:

`docs: map current GBT implementation and validation evidence`
