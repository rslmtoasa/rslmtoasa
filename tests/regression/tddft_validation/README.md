# Transverse LR-TDDFT validation fixtures

`test_validation.py` is the small deterministic CI gate. It consumes a
synthetic fixture in the exact record shapes written by `tddft_chi0`,
`tddft_goldstone`, `tddft_dyson`, and `tddft_modes`; it does not run an
electronic-structure calculation and it does not make a material claim.

Run it directly with:

```bash
python3 tests/regression/tddft_validation/test_validation.py
```

For a material campaign, create a separate directory containing a JSON
manifest with the same fields as `fixtures/ci/campaign.json`, then run:

```bash
python3 tests/regression/tddft_validation/tddft_validation.py campaign.json \
  --report evidence.json
```

The manifest is the audit boundary. `goldstone_file`, `modes_file`,
`chi0_file`, and `dyson_file` must be original TDDFT output files; their
metadata remains in the inputs. `independent_routes.GBT` and
`independent_routes.Jij` are optional read-only values from independently run
calculations. The checker never generates, tunes, patches, or otherwise
modifies GBT data.

`convergence` must list all six axes: `k_mesh`, `band_window`,
`response_projection`, `electronic_smearing`, `eta`, and `frequency_grid`.
The checked CI fixture demonstrates the schema only. Store high-accuracy Fe,
Ni, and optional Co manifests and their source outputs outside this fixture;
they are scientific evidence, not replaceable CI golden data.

## Physical periodic Fe/Ni runs

The real material fixture is `materials/run_tdval01.py`.  For a cheap
physics-oriented run and its strict analysis:

```bash
OMP_NUM_THREADS=1 python3 tests/regression/tddft_validation/materials/run_tdval01.py \
  --binary "$PWD/build/bin/rslmto.x" --material both --mesh 8 --q-set lowq \
  --continue-on-error
python3 tests/regression/tddft_validation/materials/analyze_tdval01.py
```

This checks the physical restart, periodic `strux_lib`, q=0 Goldstone,
resolved finite-q signal, energy-window placement, and a three-point
origin-constrained `omega = D|q|^2` fit.  A successful executable run is not
itself a physics pass.  The analyzer writes
`results/validation/TDVAL-01_FE_NI/runs/physics_analysis.json` and returns
nonzero when the material evidence is missing or unresolved.

## Periodic strux-lib implementation diagnostics

The current connectivity campaign is described by
`materials/campaign.json`. Every material input deck in the campaign, and the
related susceptibility/validation decks under `example/susceptibility` and
`results/validation`, explicitly sets:

```text
pbc = .true.
b1 = .true.
b2 = .true.
b3 = .true.
strux_backend = 'strux_lib'
```

The fcc-Ni decks use `ct(1)=5.0d0` and `r2=25.0d0`; this includes the 13
periodic neighbor vectors required by the connected fcc cluster. The bounded
run matrix is one Ni connectivity smoke test plus eigenpairs, K-space
Lehmann, and native real-space-GF routes for both Fe and Ni.

Run the matrix and inspect classified failures with:

```bash
python3 tests/regression/tddft_validation/materials/run_campaign.py \
  --continue-on-error
python3 tests/regression/tddft_validation/materials/analyze_campaign.py
```

The runner preserves the exact deck hash, executable, working directory,
command, exit status, and stdout/stderr under
`results/validation/TDDFT-11_FE_NI/runs/`. The analyzer checks that
`strux_lib` was selected, verifies the expected periodic-neighbor count, and
performs a diagnostic pairwise comparison of bare-`chi0` route outputs. There
is no external DFT, experiment, or converged golden response behind that
comparison, and the eigenpair route is not a trusted reference. A failed run
is classified as an input contract error, periodic neighbor mapping failure,
non-Hermitian Hamiltonian, array/index bounds failure, or generic runtime
fatal; this keeps connectivity failures separate from later response-route or
physics discrepancies. For the raw Goldstone/Ward and route summary, run:

```bash
python3 tests/regression/tddft_validation/materials/collect_tddft11_evidence.py
```
