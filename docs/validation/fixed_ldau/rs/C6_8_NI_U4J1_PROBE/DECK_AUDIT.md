# Ni U=4/J=1 surface response probe: staged-deck audit

Source snapshot: `source_snapshot/`, copied from the user's
`/Users/andersb/Jobb/RS-utils/tmp/Cheb/fccNi_small/finiteU` before edits.
The input uses `pre_processing='buildsurf'` and specifies no absolute layer-z
coordinates. The generated-surface record for this fixture reports projected
layer positions scaled by `alat`; see `C6_4_RESULTS.md` and
`C6_6_SURFACE_PROBE_RESULTS.md`.

## Surface input line changes from the copied deck

| Original | Staged | Reason |
|---|---|---|
| `nstep = 1` | `nstep = 2` | User-approved two-iteration response probe. |
| (no `conv_thr`) | `conv_thr = 5.0d-8` | Apply the requested Diff target within the two-step probe. |
| `fermi = -0.062420` | `fermi = -0.061141` | Fixed E_F observed in all 100 iterations of the regenerated U=4/J=1 bulk host. |
| `energy_max = 1.5` | `energy_max = 0.5` | Match the bulk host energy window; the lower edge remains -1.0 Ry. |
| (no `hubbard_radial_moments`) | `hubbard_radial_moments = 'auto'` | Exercise the surface auto path. |
| U/J on types 1 and 3 only | U=4/J=1 on types 1, 3, and 4 | Cover each Ni d type; type 2 ES remains without U/J. |

All other surface deck lines were retained from the user's copy. `Ni_out.nml`,
`Ni-S_out.nml`, and `Ni-S-1_out.nml` are copies of the regenerated bulk
`Ni-imp_out.nml`; only the `symbol` field is re-keyed to each deck label. The
ES potential remains the user's `ES_out.nml`. The Ni seed files and geometry
inputs are otherwise preserved.

## Host regeneration record

The host deck is copied from the accepted local RS Ni U=3/J=0.7 deck and
accepted host potential. Only the d-channel values change to U=4/J=1; explicit
local mode remains set. The host uses the accepted input E_F (-0.061141 Ry),
checked against the accepted host trace before launch. The run used the checked
wrapper and completed all 100 configured steps. Runtime reported local mode
with local energy adopted and wrote an RS-local output potential. Its final
outer Diff is approximately 7e-9, but the run reports `Not converged` because
the inherited `conv_thr` is 1e-30. This output is a diagnostic seed for the
surface probe, not a C5.2-accepted RS reference.
