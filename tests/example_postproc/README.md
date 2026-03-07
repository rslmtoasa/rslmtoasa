# Example post-processing tests

Fast example-based tests for post-processing workflows (exchange couplings,
conductivity). Each case runs in an isolated scratch directory and checks
for fatal errors. Numerical reference comparisons are optional.

## Structure

- `cases.json` — test case matrix (single source of truth for CTest and reference generation)
- `cases/<workflow>/<system>/` — input files for each calculation
- `references/<TestName>/ref.nml` — stored reference outputs (committed after generation)

## CTest modes

**Smoke** (default): runs the binary and checks for fatal errors.

**Reference** (optional): re-runs each case and compares selected output keys
against stored references. Registers additional `_ref` test variants in CTest.

### Enable smoke tests

Activate your Python venv (needs `f90nml`) before configuring:

```bash
source /path/to/venv/bin/activate
cmake -S . -B build -DRUN_EXAMPLE_TESTS=ON
ctest --test-dir build -L postproc
```

### Enable reference comparison tests

```bash
cmake -S . -B build \
    -DRUN_EXAMPLE_TESTS=ON \
    -DRUN_EXAMPLE_REF_TESTS=ON
ctest --test-dir build -L "postproc.*reference"
```

If CMake picked up the wrong Python interpreter, override without wiping the cache:

```bash
cmake -DEXAMPLE_PYTHON_EXECUTABLE=/path/to/venv/bin/python3 build
```

## Generating reference data

Run once with a known-good binary to populate `references/`.

```bash
# All post-processing cases
python3 tests/generate_references.py \
    --binary build/bin/rslmto.x \
    --cases-json tests/example_postproc/cases.json \
    --references-dir tests/example_postproc/references

# Specific case only
python3 tests/generate_references.py \
    --binary build/bin/rslmto.x \
    --cases-json tests/example_postproc/cases.json \
    --references-dir tests/example_postproc/references \
    --case Example_exchange_bccFe
```

Each reference is saved as `references/<TestName>/ref.nml` alongside a
`meta.json` recording the parameters used.

## Running a single case manually

```bash
python3 tests/run_test.py \
    --binary build/bin/rslmto.x \
    --cases-json tests/example_postproc/cases.json \
    --case-name Example_exchange_bccFe \
    --scratch-root /tmp/scratch

# With reference comparison
python3 tests/run_test.py ... \
    --compare-ref tests/example_postproc/references
```

## Adding a new test case

1. Create `cases/<workflow>/<system>/` with `input.nml` and required input files.
2. Add an entry to `cases.json`.
3. Generate its reference: `python3 tests/generate_references.py ... --case <TestName>`.
4. Commit both the new case directory and `references/<TestName>/`.
5. Re-run CMake (configure step re-reads `cases.json` automatically).
