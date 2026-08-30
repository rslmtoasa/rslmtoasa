# GBT re-audit WP08: corrected constrained frozen-magnon MFT

## Status

WP08 is implemented for the acoustic frozen-magnon branch.  The new input mode
is `mft_constrained`; `corrected_mft` is accepted as an alias.  The existing
`mft` and `scf` modes retain their prior dispatch and row-format behavior.

The corrected convention follows the fixed-potential procedure described by
Jacobsson et al.: converge the ordinary reference potential, converge a
q-specific constraining field while holding that ordinary potential fixed, then
add the converged field only to the final one-shot occupied-eigenvalue sum.  See
[Jacobsson et al., Scientific Reports 12, 18987 (2022)](https://www.nature.com/articles/s41598-022-20311-7).

## Implementation

`self%run_fixed_potential_constraint_step()` snapshots every ordinary
`potential` object, runs only `run_recursion()` and `run_dos()`, applies the
controller update, audits the temporary DOS-side change, and restores every
ordinary potential object.  It never calls `run_scf()` or density mixing.  The
separate `reset_constraint_for_fixed_potential()` method resets the field and
controller history between q points by default.

The corrected acoustic workflow in `calculation.f90` now does the following:

1. Forces the `gbt_single_q` representation and disables constraints while
   converging the ordinary reference potential.
2. Constructs a fresh constrained self object on that reference state.
3. Uses the local collinear `+/-z` GBT axis as the constraint target.  The cone
   angle remains in the GBT bond rotation; it is not duplicated as an onsite
   target field.
4. Resets the field/controller for each q when
   `constraint_start_from_zero = .true.`, iterates up to
   `constraint_max_iterations`, and fails if the requested angular tolerance is
   not reached.
5. Freezes the converged field for `frozen_magnon_probe_energy()` and reports
   `E_band(q) - E_band(q_ref)` as the corrected force-theorem Delta E.  Neither
   `constraint_penalty_energy` nor `constraint_field_coupling_energy` is added
   to that energy.

`branch_mode = 'auto'` remains available for bare MFT/SCF, but corrected MFT is
rejected there until the per-probe constrained multi-branch construction is
audited separately.

## Output contract

The legacy `frozen_magnon.dat` shape remains:

```text
q1 q2 q3 etot eband mtot_1 ... mtot_nrec omega
```

Corrected MFT additionally writes
`frozen_magnon_constrained_diagnostics.dat`.  Its metadata states the mode,
reference q, frozen-potential policy, zero-field restart policy, energy
bookkeeping, and whether a gauge subtraction was used.  Each q row contains:

* physical ordinary-reference energy and raw field-included occupied band sum;
* raw/final Delta E and the unused gauge-reference slot;
* controller residual, controller-penalty diagnostic, and `<V_con>` diagnostic;
* ordinary-potential L1 checksum, post-restore drift, and transient pre-restore
  update magnitude;
* controller iteration/convergence, cone normalization, reference moment,
  normalized omega, and `|B^C|` for every sublattice.

The post-restore ordinary-potential drift is the frozen-potential audit value;
the transient value is retained only to show that the DOS-side state was
restored rather than silently ignored.  The controller penalty and field
coupling remain diagnostics only, as established by WP07.

## Inputs and tests

The runnable Fe example is
`example/bulk/bccFe/input_frozen_magnon_constrained.nml` with
`example/bulk/bccFe/q_wp08.dat`.  It explicitly selects `strux_lib` and
disables active SOC, matching the current GBT capability gates.

`GbtWp08CorrectedMftSource` checks the mode boundary, fixed-potential
snapshot/restore, absence of SCF mixing in the inner loop, q reset, raw
field-included Delta E, and machine-readable audit fields.  The existing WP06
and WP07 controller tests remain prerequisites for the numerical controller
semantics.

## Benchmark evidence boundary

The bcc Fe smoke run was completed with the existing Fe database entry, three
q points (`Gamma`, `+/-0.05` along z), `strux_lib`, active SOC disabled, and
both controller methods 3 and 4.  The run completed and emitted the corrected
diagnostics; the converged rows showed post-restore potential drift at the
floating-point audit floor and finite transient DOS-side changes.

The prescribed fcc Ni and B2 FeCo material benchmarks are deferred.  This
working tree has no compatible standalone bulk fcc-Ni frozen-magnon fixture,
and the available FeCo case is an impurity/new-cluster setup rather than a
prepared bulk GBT frozen-magnon database.  No Ni or FeCo energetics are claimed
by this WP08 implementation.
