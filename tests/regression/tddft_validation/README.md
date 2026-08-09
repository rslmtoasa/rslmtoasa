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
