# MoS2 Qxy transport regression preflight

This directory records the inputs and evidence needed for the proposed 2D MoS2
Qxy accumulation comparison. The supplied `mos2cond_nohoh.tar` and
`mos2cond_hoh.tar` archives provide candidate pristine cases for `hoh = .false.`
and `.true.`. Each contains `input.nml`, `lattice.nml`, and Mo/S seed and
converged atom-potential files.

Both inputs select `calctype = 'B'`, a `20 x 20 x 1` real-space cell with
`b1 = b2 = .true.` and `b3 = .false.`, `linear_out = 'quadrupole_accumulation'`,
`pol_alpha = 'Qxy'`, and sheet response geometry. Their Fermi energies and
energy windows differ, so preserve each case's supplied values when comparing
the corresponding HOH setting.

## Transport constraint

Set `ccor_2c = .false.` explicitly in every run copy. The supplied inputs omit
that key, so relying on a default would make the CCOR-off constraint implicit.
Do not edit the archived originals.

## Missing comparison evidence

The archives contain no `cond_total.out` or other Qxy accumulation output, no
stdout/stderr run logs, and no provenance for a run against the reference
`main` commit. Consequently this checkout cannot establish MoS2 Qxy parity or
provide numerical golden values from the supplied material alone. The
repository's existing transport operator tests check algebraic operators and
HOH actions, while existing post-processing references cover other systems;
neither substitutes for the MoS2 accumulation curves.

Before adding a numerical MoS2 regression, collect for each HOH setting:

- `cond_total.out` and any per-type accumulation output from reference `main`;
- the complete run log and exact command;
- source commit, build/compiler configuration, and potential/input checksums;
- the output column definitions, energy reference and sheet normalization.

Then run the identical frozen-potential case against the candidate commit,
with `ccor_2c = .false.` explicit in both run copies, and compare the declared
Qxy columns and accumulation using documented tolerances. Do not generate a
reference by changing or fitting the expected output to the candidate run.

This preflight is documentation only; it makes no numerical result claim and
does not register a skipped or placeholder test.
