# Example SCF references (prepared workflow)

This folder contains the fast example SCF test runner and a prepared reference workflow.

## Files

- `run_example_scf.sh`: runs one case in an isolated scratch directory.
- `cases/<case>/`: persistent test-owned case files (`input.nml` + required atomic/lattice files).
- `cases_manifest.csv`: case matrix used for reference generation/comparison.
- `generate_reference_data.sh`: runs selected cases and stores references.
- `compare_reference_data.py`: compares run outputs (`data.nml`) against `ref.nml`.
- `run_example_scf_with_ref.sh`: runs one case and checks it against references.
- `references/<case>/ref.nml`: generated reference data (created by script).

## CTest modes

- Smoke mode (default behavior): runs examples and checks for runtime/fatal errors.
- Reference mode (optional): runs examples and compares selected `par` keys with stored references.

Enable example tests (smoke):

```bash
cmake -S . -B build_main -DRUN_EXAMPLE_TESTS=ON
```

Enable per-case reference comparison tests:

```bash
cmake -S . -B build_main \
	-DRUN_EXAMPLE_TESTS=ON \
	-DRUN_EXAMPLE_REF_TESTS=ON \
	-DEXAMPLE_REF_KEYS=etot,ws_r,vmad \
	-DEXAMPLE_REF_ABS_TOL=1e-6 \
	-DEXAMPLE_REF_REL_TOL=1e-6
```

When `RUN_EXAMPLE_REF_TESTS=ON`, extra tests with suffix `_ref` are registered in CTest.

## Generate references

Use a known-good binary and generate references from stable cases first:

```bash
bash tests/example_scf/generate_reference_data.sh build_main/bin/rslmto.x
```

By default, chebyshev cases are skipped. Include them explicitly later:

```bash
bash tests/example_scf/generate_reference_data.sh build_main/bin/rslmto.x --include-cheb
```

## Compare against references

First run the cases so `data.nml` exists in a scratch folder, then compare:

```bash
python3 tests/example_scf/compare_reference_data.py --scratch-root build_main/Testing/example_scf
```

Useful options:

- `--pattern "Example_bulk_bccFe_nsp2_*"` to restrict test names.
- `--include-cheb` to include chebyshev cases in comparison.
- `--keys etot,ws_r,vmad` to define compared keys.
- `--abs-tol 1e-6 --rel-tol 1e-6` to tune tolerance.
