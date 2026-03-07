# Example SCF tests

Fast example-based tests that run self-consistent field (SCF) calculations
and check for correctness. Each case runs in an isolated scratch directory.

## Structure

- `cases.json` — test case matrix (single source of truth for CTest and reference generation)
- `cases/<system>/` — input files for each physical system (`input.nml` + lattice/potential files)
- `references/<TestName>/ref.nml` — stored reference outputs (committed after generation)

---

## Case file format (`cases.json`)

Each entry in the `"cases"` array defines one test. The `"namelists"` dict is
applied as an `f90nml` patch on top of the case's `input.nml`, so it follows
the same section/key structure as the namelist file itself.

```json
{
  "name":    "Example_bulk_bccFe_nsp2_block_hoh_true",
  "case":    "bulk/bccFe",
  "timeout": 240,
  "namelists": {
    "control":     { "nsp": 2, "recur": "block", "lld": 20 },
    "self":        { "nstep": 1 },
    "hamiltonian": { "hoh": true }
  }
}
```

| Field       | Required | Description                                                            |
|-------------|----------|------------------------------------------------------------------------|
| `name`      | yes      | Unique test name, used as the CTest test name and scratch subdirectory |
| `case`      | yes      | Path to case input files relative to `cases/`                          |
| `timeout`   | yes      | CTest timeout in seconds (240 for block, 600 for chebyshev)            |
| `namelists` | yes      | Dict of namelist sections -> key/value pairs to patch into `input.nml` |

**`namelists` conventions:**

- `control.nsp`, `control.recur`, `control.lld` — spin panels, recursion method, recursion depth
- `self.nstep` — number of SCF steps (1 for smoke/reference tests)
- `hamiltonian.hoh` — Hamiltonian-on-Hamiltonian flag
- Any other section/key is passed through as-is. Only listed keys are overwritten;
  all other keys in `input.nml` keep their original values.

**Example with extra energy range (needed for chebyshev):**

```json
{
  "name":    "Example_bulk_bccFe_nsp2_chebyshev_hoh_true",
  "case":    "bulk/bccFe",
  "timeout": 600,
  "namelists": {
    "control":     { "nsp": 2, "recur": "chebyshev", "lld": 200 },
    "self":        { "nstep": 1 },
    "hamiltonian": { "hoh": true },
    "energy":      { "energy_min": -3.0, "energy_max": 1.8 }
  }
}
```

---

## Adding a new test case

1. **Create the case directory** `cases/<system>/` with `input.nml` and any
   required input files (potential, lattice, etc.). The `input.nml` should
   contain sensible defaults; `namelists` in `cases.json` will override
   specific keys at run time.

2. **Add an entry to `cases.json`** following the format above.
   - Pick a descriptive `name` (used as the CTest test name).
   - Set `timeout` to 240 (block) or 600 (chebyshev).
   - Only list the keys you want to force; omit anything the `input.nml`
     default is already correct for.

3. **Generate the reference:**

   ```bash
   python3 tests/generate_references.py \
       --binary build/bin/rslmto.x \
       --cases-json tests/example_scf/cases.json \
       --references-dir tests/example_scf/references \
       --case <TestName>
   ```

4. **Commit** the new `cases/<system>/` directory, the updated `cases.json`,
   and `references/<TestName>/`.

5. **Reconfigure CMake** — the configure step re-reads `cases.json`
   automatically and registers the new CTest entry.

---

## CTest modes

**Smoke** (default): runs the binary, checks for fatal errors.

**Reference** (optional): re-runs each case and compares selected `&par` keys
against stored references. Registers additional `_ref` test variants.

### Enable smoke tests

Activate your Python venv (needs `f90nml`) before configuring:

```bash
source /path/to/venv/bin/activate
cmake -S . -B build -DRUN_EXAMPLE_TESTS=ON
ctest --test-dir build -L example
```

### Enable reference comparison tests

```bash
cmake -S . -B build \
    -DRUN_EXAMPLE_TESTS=ON \
    -DRUN_EXAMPLE_REF_TESTS=ON
ctest --test-dir build -L reference
```

Tolerances and keys can be tuned at configure time (defaults shown):

```bash
cmake ... \
    -DEXAMPLE_REF_KEYS=etot,ws_r,vmad \
    -DEXAMPLE_REF_ABS_TOL=1e-6 \
    -DEXAMPLE_REF_REL_TOL=1e-6
```

If CMake picked up the wrong Python interpreter (e.g. you activated the venv
after a previous configure), override it without wiping the cache:

```bash
cmake -DEXAMPLE_PYTHON_EXECUTABLE=/path/to/venv/bin/python3 build
```

---

## Generating reference data

Run once with a known-good binary to populate `references/`. Results are
committed so CI does not need to regenerate them.

```bash
# All cases
python3 tests/generate_references.py \
    --binary build/bin/rslmto.x \
    --cases-json tests/example_scf/cases.json \
    --references-dir tests/example_scf/references

# Specific cases only
python3 tests/generate_references.py \
    --binary build/bin/rslmto.x \
    --cases-json tests/example_scf/cases.json \
    --references-dir tests/example_scf/references \
    --case Example_bulk_bccFe_nsp2_block_hoh_true
```

Each reference is saved as `references/<TestName>/ref.nml` alongside a
`meta.json` recording the full case definition used.

---

## Running a single case manually

```bash
# Smoke run
python3 tests/run_test.py \
    --binary build/bin/rslmto.x \
    --cases-json tests/example_scf/cases.json \
    --case-name Example_bulk_bccFe_nsp2_block_hoh_true \
    --scratch-root /tmp/scratch

# With reference comparison
python3 tests/run_test.py \
    --binary build/bin/rslmto.x \
    --cases-json tests/example_scf/cases.json \
    --case-name Example_bulk_bccFe_nsp2_block_hoh_true \
    --scratch-root /tmp/scratch \
    --compare-ref tests/example_scf/references
```
