# fcc Pt conductivity comparison

The two response directories use the supplied Pt checkpoints and the current
`linear_in`/`linear_out` interface:

| Directory | Input | Output |
| --- | --- | --- |
| `sigma_xy_z-spin` | Charge velocity along x | Spin current along y, spin polarization z |
| `sigma_xy_z-orb` | Charge velocity along x | Orbital current along y, orbital polarization z |

The historical directory names are retained. Read the explicit directions:
`v_alpha = 0,1,0` is the output direction and `v_beta = 1,0,0` is the input.
Do not exchange these axes when comparing with an old result.

## Run

Use a freshly built executable from the current code. Run from each response
directory so it reads that directory's `input.nml` and `Pt_out.nml`:

```sh
cd /path/to/rslmtoasa/example/conductivity/fccPt/sigma_xy_z-spin
/path/to/build/bin/rslmto.x > run.log 2>&1

cd ../sigma_xy_z-orb
/path/to/build/bin/rslmto.x > run.log 2>&1
```

Each response input selects `post_processing='conductivity'`, which builds the
bulk geometry and Hamiltonian before computing moments. The parent `input.nml`
retains its original `pre_processing='bravais'` setup; it is not a conductivity
run. A new self-consistent calculation is not needed for this comparison.

## Comparison settings

The migration replaces `cond_type` with `linear_in='charge'` and the corresponding
`linear_out`, and replaces `js_alpha`/`jl_alpha` with `pol_alpha='z'`.
`pol_beta='z'` is explicit but unused by the charge input. `hoh=.false.` is explicit
and matches the original default.

All existing numerical settings and checkpoint files are retained: SOC mode
`nsp=2`, `cond_calctype='per_type'`, `cond_ll=500`, `lld=150`, energy range
[-1.0, 1.2], Fermi energy -0.089509, 2500 LDOS channels, and the original lattice.
Use the same checkpoint potentials, axes, energy settings and post-processing
conventions for the old and new runs. Keep their output files in separate copies
of the example directories; the solver replaces existing output files.

For a separate HOH comparison, copy a response directory and set `hoh=.true.`
in its existing `&hamiltonian` group. Do not add a second group. The corrected
HOH velocity contains both -v o h and -h o v, so agreement with an old HOH result
is not expected in general.

This example update validates input compatibility, not equivalence of complete
old/new conductivity curves. No new material calculation or change to the common
energy-integration, normalization or units code is included.
