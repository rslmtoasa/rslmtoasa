# Physical Fe/Ni TD-DFT fixtures

The material-level check is `run_tdval01.py` plus `analyze_tdval01.py`.  It
uses the physical Fe and Ni restart potentials, periodic `strux_lib`, the
explicit transition-pole response, and a q=0/low-q frequency calculation.
It is deliberately a physics diagnostic rather than a synthetic CI fixture:
the checker requires a small Goldstone residual, a resolved positive-
frequency collective signal, and at least three resolved nonzero-q points
before it reports an origin-constrained `omega = D |q_cartesian|^2` fit.
Stoner features, boundary peaks, and rejected/non-Lorentzian fits do not count.

The cheap low-q run is:

```sh
OMP_NUM_THREADS=1 python3 tests/regression/tddft_validation/materials/run_tdval01.py \
  --binary "$PWD/build/bin/rslmto.x" \
  --material both --mesh 8 --q-set lowq --continue-on-error
python3 tests/regression/tddft_validation/materials/analyze_tdval01.py
```

The low-q ladder is direct `(0,0,0)`, `(0.005,0,0)`, `(0.010,0,0)`, and
`(0.015,0,0)`.  It uses `eta=2e-5 Ry`, a `0..0.005 Ry` window, and 251
frequency points so that a claimed signal is resolved rather than merely
visible on a broad grid.  Use `--allow-inconclusive` only when collecting
diagnostics; it does not turn a failed Goldstone or response check into a
pass.  The report is written to
`results/validation/TDVAL-01_FE_NI/runs/physics_analysis.json`.

For mesh convergence, repeat with `--mesh 8 12 16`; for the wider historical
q/eta comparison use `--q-set both`.  Those runs are more expensive and are
not required for the cheap first diagnostic.

## TDDFT-11 implementation-route diagnostics

The older TDDFT-11 material campaign below remains useful for checking the
three bare-chi0 providers and periodic neighbor construction.  It is not a
material validation and must not be used to claim a Fe/Ni signal or a spin-
wave stiffness.

## TDDFT-11 material campaign decks

These are reproducible current-branch smoke decks that emit raw output from
the three transverse `chi0` providers at selected `(q, omega)` points. They use
the committed Fe and Ni restart potentials from `results/validation/VAL-18_bccFe`
and `results/validation/VAL-19_fccNi`; no ground-state file is modified. They
are diagnostic runs, not a comparison against experiment, a converged DFT
benchmark, or an accepted golden response.

The three decks for each material differ only in `chi0_backend` and output
prefix, except that the native real-space deck emits bare `chi0` only because
the current provider has no validated exact static kernel for Dyson output.
All q coordinates are direct reciprocal coordinates and all energies are Ry.

Run from the material directory with the built executable, passing each deck
explicitly:

```sh
OMP_NUM_THREADS=1 ../../../../../build/bin/rslmto.x input_eigenpairs.nml
OMP_NUM_THREADS=1 ../../../../../build/bin/rslmto.x input_kspace_lehmann.nml
OMP_NUM_THREADS=1 ../../../../../build/bin/rslmto.x input_realspace_gf.nml
```

The run report records the exact commands and resulting backend/evidence
status. The `.dat` files are the raw outputs intended for implementation and
physics review; the machine-readable collector does not turn them into a
release pass.

All six route decks and the Ni connectivity probe explicitly select the same
periodic strux-lib contract:

```text
pbc = .true.
b1 = .true.
b2 = .true.
b3 = .true.
strux_backend = 'strux_lib'
```

The complete bounded campaign can be launched from the repository root with:

```sh
python3 tests/regression/tddft_validation/materials/run_campaign.py \
  --continue-on-error
python3 tests/regression/tddft_validation/materials/analyze_campaign.py
```

For the simplest dispersion review bundle—one periodic eigenpair-transition
run over a few q-vectors for each material—use:

```sh
python3 tests/regression/tddft_validation/materials/produce_tddft_test_data.py
```

It creates a fresh timestamped directory under
`results/validation/TDDFT-test-data/` and prints JSON containing the complete
file list. The bundle includes the effective decks, q-point inputs, dispersion
mode files, per-q response files, Goldstone output, logs, and runtime
diagnostics; it does not compare the results with a reference.

The run summary and per-case logs are written below
`results/validation/TDDFT-11_FE_NI/runs/`. The analyzer first checks the
periodic neighbor count (15 for bcc Fe and 13 for fcc Ni), then performs a
diagnostic pairwise comparison of the three bare-`chi0` outputs. This is not an
absolute reference: all routes use the same restart state and share physical
assumptions. Route mismatches must be interpreted together with the output
metadata: the native real-space route is not converged unless its source radius
covers the requested response radius. The structured evidence record also
retains the raw Goldstone/Ward diagnostics.

## Native local zones

For `chi0_backend='realspace_gf'`, an automatic pair list is generated for
each response site against the embedding cluster. `realspace_rmax` is applied
in physical Angstroms before intersite recursion and Green-function storage,
so a large embedding cluster can use a smaller response zone. If an explicit
`ijpair` list is supplied, it remains authoritative. A local zone must still
be converged by an `rmax`/shell sweep; the cutoff alone is not a release
accuracy claim.
