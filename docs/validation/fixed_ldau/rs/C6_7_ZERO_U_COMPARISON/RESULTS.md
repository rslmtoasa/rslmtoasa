# Ni(001) no-U versus explicit-zero-U comparison

Date: 2026-10-25. Source: user-run decks in
`/Users/andersb/Jobb/RS-utils/tmp/Cheb/fccNi_small/{noU,zeroU}`.

## Gate result

**PASS: the explicit-zero-U deck reproduces the no-U deck's physical output.** The decks differ
only by the four explicit zero assignments for the d channel on types 1 and 3. After removing
timestamps and timing-report values, the solver outputs agree. Both runs report the same one-step
Diff (`7.71e-8`), spin moments (atoms 1-3: `0.007355`, `0.689005`, `0.640933` μB), active-site
charge transfers (`+0.225385`, `-0.295649`, `+0.079138` e), and external-layer charge
(`-0.008877` e). The zero-U comparison gate is independent of convergence: it verifies that
explicit U=J=0 is equivalent to omitting U/J for these otherwise matching runs.

The common Diff is above the previously stated `5e-8` surface response threshold. That threshold
does not change the zero-U equivalence result. The user has directed that it is not an appropriate
acceptance threshold for this surface comparison and should be raised; no replacement numerical
threshold has yet been specified. These one-step runs do not test nonzero U.

The zero-U output contains the assignments `hubbard_u_general(1,3)=0`,
`hubbard_j_general(1,3)=0`, `hubbard_u_general(3,3)=0`, and
`hubbard_j_general(3,3)=0`; the no-U input omits them. Both outputs report the same radial-mode
messages and contain no active fixed-U/J summary. No SCF preflight record accompanied these
user-run artifacts, so this is recorded specifically as the zero-U deck-equivalence comparison,
not as a checked C6.5 surface smoke run.

The exact inputs and outputs are retained under `user_runs/`:

| Artifact | SHA256 |
|---|---|
| `noU_input.nml` | `1b195a50203dbd7771982f37633ba30e10f76796d9e5bcf9dc0eb2cec42f9241` |
| `zeroU_input.nml` | `0d296fd157026e5297ee9098e64f28e671fed3eac950fed23194f3b0bb90700c` |
| `out.noU` | `1cd935184fa111880f1144c030e5caa9a27e075cbfc0ba0014da4269c32f7d4f` |
| `out.zeroU` | `bcd7097d886527914a454f459884b24f4ef46cf4c55691d97bc79c601c918cc4` |

## Checked nonzero-U diagnostics (2026-10-25)

The two user runs above establish the zero-U equivalence gate. To continue toward the nonzero-U
surface response, two checked one-step Chebyshev diagnostics were staged from the same Ni(001)
geometry, with the accepted RS-local Ni U=3/J=0.7 host potential copied to the Ni surface and
subsurface types, the user-run ES seed retained, auto mode, and host E_F `-0.061141` Ry. Both
passed the checked preflight and selected `mode=local`; the log contains no provenance error and
reports the expected fixed U/J and finite on-site Hubbard potential.

Both Chebyshev runs fail downstream: projected-DOS arrays and `g0` become NaN for recursion atoms
1–3, followed by `Invalid local d-shell Hubbard moment correction; status=5 type=3 atom=2 spin=1`.
The first used upper energy window `0.5` Ry. The second changed only that upper bound to `1.5` Ry,
matching the user's finite no-U deck, and failed at the same stage. Thus the upper window alone
does not explain this failure. No valid U-dependent occupations or moments were produced. The
first run's `conv_thr=1e-7` was a diagnostic choice, not an approved surface convergence gate.

A third checked diagnostic changed only the recurrence from Chebyshev to block. It reached the
first recursion step and logged the local Hubbard correction, then produced no further output for
several minutes. It was interrupted; this is inconclusive about the block path and supplies no
observable. The host environment prevented process inspection while it was silent.

The complete inputs, preflight records, logs, and generated files for all three diagnostics are
archived in `nonzero_u_diagnostics/`; `nonzero_u_diagnostics/SHA256SUMS` records their hashes.
These failures are independent of the passing zero-U equivalence comparison. The nonzero-U surface
response remains OPEN. No replacement numerical Diff threshold for surface probes has been
specified, so these diagnostics are not judged against the previous `5e-8` response threshold.

### Lower-bound diagnostic

A fourth checked one-step U=3/J=0.7 auto run kept the `1.5` Ry upper edge and widened only the
lower edge from `-1.0` to `-3.0` Ry. It passed preflight with host E_F delta `0.0` Ry, resolved
`mode=local`, produced no NaN projected DOS / `g0`, passed through the local moment correction,
and completed the SCF step. Layer spin moments were `0.004284`, `0.001571`, and `0.000077` μB for
recursion atoms 1–3; total RMS Diff was `0.0221012697`. The run had `nstep=1`, so these are
single-step diagnostic values, not converged layer moments. The U/J summary reports d occupations
`nup=ndn=4.738761` for each U-carrying type at the start of the step.

Changing only the lower Chebyshev bound removed the earlier NaN failure, while widening only the
upper edge did not. This strongly implicates the lower spectral bound / Chebyshev scaling window
for those NaNs; it does not prove every affected eigenvalue or establish a converged physical
response. The recurrence implementation maps the Hamiltonian into the Chebyshev interval and its
fast DOS reconstruction evaluates `acos((E-b)/a)` and `1/sqrt(a*a-(E-b)**2)`; an inadequate bound
can make those expressions non-finite. The completed diagnostic deck, preflight, and outputs are
archived in `nonzero_u_diagnostics/chebyshev_widened_both_edges/` and included in the hash manifest.

This result shows that the U surface path can complete a finite one-step local-moment update with
a wider lower window. It is not a direct nonzero-U versus zero-U physical comparison: the user
zero-U pair used different host/seed potentials and a different energy window, and the U run has
no matched zero-U control on the same local-U host and widened window. The no-U versus explicit
zero-U equivalence gate remains PASS independently of convergence. A matched control would be
needed to quantify the U-induced change; no surface energy is routed through the Hubbard term.
