# Example post-processing tests

Fast example-based tests for post-processing workflows (exchange couplings,
conductivity). Each case runs in an isolated scratch directory and checks
for fatal errors. Numerical reference comparisons are optional.

## Structure

- `cases.json` — test case matrix (single source of truth for CTest and reference generation)
- `cases/<workflow>/<system>/` — input files for each calculation
- `references/<TestName>/ref.json` — stored reference outputs (committed after generation)

## Case file format (`cases.json`)

Each entry defines one test. The `"namelists"` dict patches `input.nml`; the
optional `"checks"` dict defines what to compare against stored references.

```json
{
  "name": "Example_exchange_bccFe",
  "case": "exchange/bccFe",
  "timeout": 240,
  "namelists": {
    "control":     { "nsp": 2, "recur": "block", "lld": 20 },
    "self":        { "nstep": 1 },
    "hamiltonian": { "hoh": false }
  },
  "checks": {
    "nml": [
      { "file": "Fe_out.nml", "scalars": ["etot", "ws_r"], "arrays": { "mom": [3] } }
    ],
    "text": [
      { "file": "jij.out", "rows": [1, 2], "cols": [6, 7] }
    ]
  }
}
```

**`checks` fields:**

- `nml` — list of output namelist files to check:
  - `file` — filename in the workdir
  - `scalars` — scalar keys to compare (searched across all namelist sections)
  - `arrays` — array key → list of 1-based Fortran indices to compare
- `text` — list of space-separated text files to check:
  - `file` — filename in the workdir
  - `rows` — 1-based line numbers to check
  - `cols` — 1-based column indices; all columns are checked for every row

Cases without a `"checks"` key run as smoke-only (log check, no value comparison).

## Running the tests

Every test runs the binary, checks the log for fatal errors, and compares the
values defined in the case's `"checks"` dict against `references/<TestName>/ref.json`.
If no `ref.json` exists yet the test passes as smoke-only.

Activate your Python venv (needs `f90nml`) before configuring:

```bash
source /path/to/venv/bin/activate
cmake -S . -B build -DRUN_EXAMPLE_TESTS=ON
ctest --test-dir build -L postproc
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
    --cases-json tests/postproc/cases.json \
    --references-dir tests/postproc/references

# Specific case only
python3 tests/generate_references.py \
    --binary build/bin/rslmto.x \
    --cases-json tests/postproc/cases.json \
    --references-dir tests/postproc/references \
    --case Example_exchange_bccFe
```

Each reference is saved as `references/<TestName>/ref.json` alongside a
`meta.json` recording the parameters used. The values extracted are those
specified by the case's `"checks"` dict in `cases.json`.

## Running a single case manually

```bash
python3 tests/run_test.py \
    --binary build/bin/rslmto.x \
    --cases-json tests/postproc/cases.json \
    --case-name Example_exchange_bccFe \
    --scratch-root /tmp/scratch \
    --compare-ref tests/postproc/references
```

## Adding a new test case

1. Create `cases/<workflow>/<system>/` with `input.nml` and required input files.
2. Add an entry to `cases.json`.
3. Generate its reference: `python3 tests/generate_references.py ... --case <TestName>`.
4. Commit both the new case directory and `references/<TestName>/`.
5. Re-run CMake (configure step re-reads `cases.json` automatically).
