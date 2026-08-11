# fcc Ni transverse TDDFT magnons

This lightweight example evaluates the transverse site-projected LR-TDDFT
response of ferromagnetic fcc Ni along a short Gamma-to-small-q path. It is a
repair diagnostic, not a converged material prediction or a validation claim.

## Prerequisite

Run the block-recursion SCF calculation in `../SCF/` first.  The TDDFT deck
reads its converged `Ni_out.nml` through `database = '../SCF/'`; it does not
alter that restart file.

## Run

From this directory:

```sh
../../../../build_13/bin/rslmto.x
```

The response calculation refreshes radial XC data from the saved
spin-polarized SCF potential, then evaluates the 89-point
Gamma-X-U-K-Gamma-L-W-X path in `fccNi_qpath.dat` on a 12x12x12 response mesh.
It uses 401 frequencies between 0 and 0.04 Ry with `eta = 0.0005` Ry. This
deck selects the `pair_potential` Xi backend and `goldstone_mode = 'correct'`;
it writes raw and controlled-correction products separately. It is a repair
diagnostic, not a material-validation result.

The raw and corrected files record their own minimum site spectral weights.
For a finite-broadened circular response, a low-frequency tail can be present
in both raw and corrected products. The run stops only if the static column
correction creates negative weight or makes an existing negative tail more
negative; inspect the correction-preservation metadata in the corrected file.

The principal files are:

- `fccNi_tddft_sr_goldstone.dat` — raw static and loss-grid Goldstone
  diagnostics plus correction provenance;
- `fccNi_tddft_sr_q*_chi0.dat` — Kohn-Sham transverse susceptibility;
- `fccNi_tddft_sr_q*_pair_dyson.dat` and
  `fccNi_tddft_sr_q*_pair_corrected_dyson.dat` — raw and corrected
  pair-potential Dyson responses;
- `fccNi_tddft_sr_*_modes.dat` — heuristic mode diagnostics, not a validated
  collective-mode classifier;
- `fccNi_tddft_sr_manifest.dat` — deterministic q-point/file index.

## Plotting a spectrum

The included script makes a q--energy intensity map from the Dyson loss files
and also writes the strongest-loss peak at every q:

```sh
python3 -m pip install --user numpy matplotlib   # once, if needed
python3 plot_magnon_spectrum.py --energy-max-mev 500
```

This writes `fccNi_magnon_spectrum.png` and `fccNi_magnon_peaks.dat`.  The
cyan line is a convenient visual guide, not automatically a magnon branch:
it tracks the largest loss and may instead follow Stoner weight.  For a
multi-site calculation, select a site with `--site N`; use `--quantity
trace_loss` to plot the site-summed loss.  `--normalize-each-q` is useful for
revealing a weak dispersing feature, but should not be used to compare its
absolute intensity between q points.  When `fccNi_qpath.dat` is present, the
plot automatically marks Γ--X--U--K--Γ--L--W--X.

When raw and corrected pair-potential responses are present, the script uses
the corrected pair product by default.  Use `--dyson-product raw_pair` or
`--dyson-product legacy` to plot those alternatives explicitly.

With a finite k mesh, a few small negative loss samples can occur from
numerical integration noise.  The script reports and clips those values to
zero only for the color map; it never takes their absolute value or changes
the raw TDDFT files.  If they become appreciable rather than isolated small
values, improve the k mesh and check convergence before interpreting a peak.

For a physical Ni magnon dispersion, increase the k mesh substantially,
converge the SCF threshold, sample a symmetry-aware high-symmetry path, and
converge `eta`.  The values here deliberately prioritize runtime.

`fermi_level` is not a TDDFT input.  The 12x12x12 response mesh resolves the
chemical potential for the restart's target electron count; do not copy an
SCF-mesh Fermi value into `&tddft`.
