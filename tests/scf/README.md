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
  },
  "checks": {
    "nml": [
      { "file": "Fe_out.nml", "scalars": ["etot", "ws_r"], "arrays": { "mom": [3] } }
    ],
    "text": [
      { "file": "totaldos.out", "rows": [500, 1000, 1500], "cols": [1, 2] }
    ]
  }
}
```

| Field          | Required | Description                                                                               |
|----------------|----------|-------------------------------------------------------------------------------------------|
| `name`         | yes      | Unique test name, used as the CTest test name and scratch subdirectory                    |
| `case`         | yes      | Path to case input files relative to `cases/`                                             |
| `timeout`      | yes      | CTest timeout in seconds (300 for block, 600 for chebyshev)                               |
| `namelists`    | yes      | Dict of namelist sections -> key/value pairs to patch into `input.nml`                    |
| `checks`       | no       | What to compare against stored references (see below)                                     |
| `mpi_procs`    | no       | MPI rank count for this case (only honored when CMake is configured with `ENABLE_MPI=ON`) |
| `launch_modes` | no       | `["serial", "mpi"]` to register a serial/MPI pair sharing one reference (see below)       |

**Serial-vs-MPI consistency (`launch_modes`):**

A case with `mpi_procs` set normally runs *only* under MPI (when `ENABLE_MPI=ON`),
so a rank-decomposition bug that changes results relative to serial execution
would not be caught. Adding `"launch_modes": ["serial", "mpi"]` to a case
registers **two** CTest entries instead of one:

- `<TestName>` — forced serial (`mpi_procs=1`), regardless of `--mpi-enabled`.
- `<TestName>_mpi` — runs with the case's `mpi_procs`, registered only when
  `ENABLE_MPI=ON`, with the CTest `PROCESSORS` property set so ctest
  scheduling accounts for the rank count.

Both variants compare against the **same** `references/<TestName>/ref.json` —
the point is to confirm MPI decomposition reproduces the serial result, not to
maintain two sets of references. This is applied to a representative subset
(one bulk, the surface, the impurity, one chebyshev case) rather than every
MPI-capable case — KISS, CI minutes are finite.

**`namelists` conventions:**

- `control.nsp`, `control.recur`, `control.lld` — spin panels, recursion method, recursion depth
- `self.nstep` — number of SCF steps (1 for smoke/reference tests)
- `hamiltonian.hoh` — Hamiltonian-on-Hamiltonian flag
- Any other section/key is passed through as-is. Only listed keys are overwritten;
  all other keys in `input.nml` keep their original values.

**`checks` format:**

The optional `"checks"` dict describes what to extract from the output and compare against `references/<TestName>/ref.json`. If omitted, the test runs as smoke-only (no value comparison).

```json
"checks": {
  "nml": [
    {
      "file":    "Fe_out.nml",
      "scalars": ["etot", "ws_r"],
      "arrays":  { "mom": [3], "xi_d": [1, 2] }
    }
  ],
  "text": [
    { "file": "totaldos.out", "rows": [500, 1000, 1500], "cols": [1, 2] }
  ]
}
```

- `nml` — list of output namelist files to check. Each entry specifies:
  - `file` — filename in the workdir (e.g. `Fe_out.nml`)
  - `scalars` — list of scalar keys to compare (searched across all namelist sections)
  - `arrays` — dict of array key → list of 1-based Fortran indices to compare
- `text` — list of space-separated text output files to check. Each entry specifies:
  - `file` — filename in the workdir
  - `rows` — list of 1-based line numbers to check
  - `cols` — list of 1-based column indices; all specified columns are checked for every row

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
  },
  "checks": {
    "nml":  [ { "file": "Fe_out.nml", "scalars": ["etot", "ws_r"] } ],
    "text": [ { "file": "totaldos.out", "rows": [500, 1000, 1500], "cols": [1, 2] } ]
  }
}
```

---

## Cross-calctype oracle (fcc Cu)

Four cases describe **the same physical atom** — Cu in bulk fcc Cu — reached by
three different code paths, so they function as an oracle for each other rather
than only as stored-value regressions:

| Test | `calctype` | Route |
| --- | --- | --- |
| `Example_bulk_fccCu_chebyshev` | `B` | `bravais` |
| `Example_impurity_fccCu_chebyshev` | `I` | `newclubulk` — Cu "impurity" on a Cu host |
| `Example_interface_fccCu001_chebyshev` | `L` | `buildinterface` — one Cu (001) layer between Cu regions |
| `Example_interface_fccCu111_chebyshev` | `L` | `buildinterface` — same, (111) |

All four share one converged Cu parameter set (`vmad ~ 0`), and all pin
`fermi = -0.089874` with `fix_fermi = .true.`. **The pinned Fermi level is load
bearing:** without it each route determines E_F independently, the energy meshes
shift relative to one another, and a pointwise DOS comparison compares different
grids rather than different physics. With it, all four emit an identical mesh
(2009 rows, −1.40838 … 2.10469), so `*_dos.out` rows are directly comparable.

Cross-checking the committed references (`etot`, and `*_dos.out` col 2):

- `etot` agrees to ~1e-8 Ry across all four; `ws_r` is exact.
- `L001` reproduces bulk **exactly** at every sampled row.
- `I` deviates from bulk by ≤1e-5 — the print precision of the `.out` files.
- **`L111` deviates by 2.05e-3 at row 1200** (near the d-band peak), ~200×
  larger than `I` or `L001`. This is a real, orientation-specific residual, not
  a tolerance artifact. It is *not* explained, and is deliberately captured
  rather than tolerated away — the reference pins the current value so any
  change in it is visible. See `tests/KNOWN_ISSUES.md`.

This suite compares each case against its own stored reference; it does not
assert B ≡ I ≡ L programmatically. When touching the Hamiltonian build, the
recursion, or any `calctype` dispatch, compare the four references to each other
by hand — that comparison is what catches a route-specific bug, and it is how
the B7.5 zero-DOS and representative-site bugs were found.

**Do not compare converged output from two directories that started from
different parameters or different E_F.** That is not a reference pair, and it
produces differences (~0.26 in the DOS peak region) that have nothing to do with
the code path under test.

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
       --cases-json tests/scf/cases.json \
       --references-dir tests/scf/references \
       --case <TestName>
   ```

4. **Commit** the new `cases/<system>/` directory, the updated `cases.json`,
   and `references/<TestName>/`.

5. **Reconfigure CMake** — the configure step re-reads `cases.json`
   automatically and registers the new CTest entry.

---

## Running the tests

Every test runs the binary, checks the log for fatal errors, and compares the
values defined in the case's `"checks"` dict against `references/<TestName>/ref.json`.
If no `ref.json` exists yet the test passes as smoke-only.

Activate your Python venv (needs `f90nml`) before configuring:

```bash
source /path/to/venv/bin/activate
cmake -S . -B build -DRUN_EXAMPLE_TESTS=ON
ctest --test-dir build -L scf
```

Tolerances can be tuned at configure time (defaults shown):

```bash
cmake ... \
    -DEXAMPLE_REF_ABS_TOL=1e-6 \
    -DEXAMPLE_REF_REL_TOL=1e-6
```

## Tuning tolerances and OpenMP threads

- **Global tolerances (configure-time)**: set `EXAMPLE_REF_ABS_TOL` and
  `EXAMPLE_REF_REL_TOL` at CMake configure time (examples above). These are
  passed to the test runner as the default comparison tolerances.

- **Per-case overrides (manifest)**: a case entry in `cases.json` may include
  the optional fields `abs_tol`, `rel_tol`, and `omp_threads` to override the
  global defaults for that specific case. Example:

```json
{
  "name": "Example_bulk_bccFe_nsp2_block_hoh",
  "case": "bulk/bccFe",
  "timeout": 240,
  "omp_threads": 2,
  "abs_tol": 1e-6,
  "rel_tol": 1e-6,
  "namelists": { ... }
}
```

- **How it's applied**: the `tests/run_test.py` runner will use `abs_tol` and
  `rel_tol` from the case manifest when present; otherwise it falls back to the
  CLI/CMake defaults. For `omp_threads` the runner will propagate the value
  into the binary wrapper which sets `OMP_NUM_THREADS` for serial runs.

- **Environment override (manual runs / CI)**: the wrapper `tests/run_binary.sh`
  reads the environment variable `RSLMTO_OMP_THREADS_SERIAL` (default `1`) and
  uses it for serial runs. You can export a value to override globally for a
  test session:

```bash
export RSLMTO_OMP_THREADS_SERIAL=2
python3 tests/run_test.py --binary build/bin/rslmto.x --cases-json tests/scf/cases.json \
    --case-name Example_bulk_bccFe_nsp2_block_hoh --scratch-root /tmp/scratch --compare-ref tests/scf/references
```

- **MPI runs**: MPI-invoked cases (those with `mpi_procs` > 1) keep `OMP_NUM_THREADS=1`
  to avoid mixing MPI and multi-threading; per-case `omp_threads` will be ignored
  for MPI cases.

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
    --cases-json tests/scf/cases.json \
    --references-dir tests/scf/references

# Specific cases only
python3 tests/generate_references.py \
    --binary build/bin/rslmto.x \
    --cases-json tests/scf/cases.json \
    --references-dir tests/scf/references \
    --case Example_bulk_bccFe_nsp2_block_hoh_true
```

Each reference is saved as `references/<TestName>/ref.nml` alongside a
`meta.json` recording the full case definition used.

---

## Running a single case manually

```bash
python3 tests/run_test.py \
    --binary build/bin/rslmto.x \
    --cases-json tests/scf/cases.json \
    --case-name Example_bulk_bccFe_nsp2_block_hoh_true \
    --scratch-root /tmp/scratch \
    --compare-ref tests/scf/references
```
