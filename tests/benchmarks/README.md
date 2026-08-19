# ACC-00 accelerator benchmark harness

`benchmark_harness.py` records repeated command timings as JSON. It is an
informational benchmark tool: it does not set a timing threshold and it is not
an ordinary CTest correctness gate.

The output schema is `rslmto.accelerator-benchmark.v1`. Each record includes
the command, warm-up/repetition samples, required environment metadata, and
the dimensions supplied by the caller. Missing accelerator fields are kept as
`null` so CPU and later GPU records have the same shape. Existing
`UnitTddftCpuProfile` output is parsed into separate reciprocal phase records
for Fourier assembly, normal eigensolution, arbitrary-k eigensolution, and
memory estimates.

## Run a component profile

```bash
python3 tests/benchmarks/benchmark_harness.py run \
  --name reciprocal_cpu_profile \
  --class component \
  --labels performance component reciprocal eigensolver fourier \
  --build-dir build \
  --output results/benchmarks/reciprocal_cpu_profile.json \
  --command build/bin/UnitTddftCpuProfile
```

The command must be placed last because all remaining arguments after
`--command` are passed to the executable. Use `--warmups 1 --repetitions 5`
for a normal measurement campaign. The default is one warm-up and three
recorded runs.

## Run production fixtures

`manifest.json` is the initial ACC-00 inventory. Its reciprocal entries use
the established Si/sp and bcc-Fe/spd k-space cases; its RS entries use the
established Fe Block/Lanczos and Si Chebyshev/moment routes. Run the inventory
with a trusted production binary and scratch directory:

```bash
python3 tests/benchmarks/benchmark_harness.py run-manifest \
  --manifest tests/benchmarks/manifest.json \
  --binary build-acc00/bin/rslmto.x \
  --profile-binary build-acc00/bin/UnitTddftCpuProfile \
  --build-dir build-acc00 --warmups 1 --repetitions 3 \
  --scratch-root /tmp/rslmto-acc00 \
  --output-dir results/benchmarks
```

For RS CUDA timings, add `--gpu-plugin`. The manifest opts the Block and
Chebyshev RS entries into `control%gpu_plugin=true`; reciprocal entries remain
CPU routes, and the current nsp=2 scalar-Lanczos fixture remains CPU-only by
design.

Run one entry with `--name reciprocal_si_sp`, or omit `--profile-binary` when
running only production entries; the optional `reciprocal_cpu_profile` entry
is skipped automatically in that case. For a single production command, the
equivalent explicit form is:

```bash
python3 tests/benchmarks/benchmark_harness.py run \
  --name reciprocal_si_sp --class end_to_end \
  --labels performance end_to_end reciprocal eigensolver fourier \
  --k-points 8 --matrix-dimension 18 --build-dir build \
  --output results/benchmarks/reciprocal_si_sp.json \
  --command python3 tests/run_test.py --binary build/bin/rslmto.x \
    --cases-json tests/scf/cases.json \
    --case-name Example_k_space_scf_diamondSi_sp \
    --scratch-root /tmp/rslmto-bench/reciprocal_si_sp \
    --force-serial
```

The command's production smoke/output checks remain separate from committed
correctness references; add `--compare-ref tests/scf/references` when the
reference set is known to match the selected binary. The harness only adds
wall-time evidence around that route.

## Compare runs

```bash
python3 tests/benchmarks/benchmark_harness.py compare \
  results/benchmarks/cpu.json results/benchmarks/gpu.json \
  --output results/benchmarks/cpu-vs-gpu.json
```

The comparison reports median, minimum, spread, and `baseline/candidate`
speedup for every shared metric. It emits environment mismatch warnings but
never fails because a timing changed by a percentage.

## Metadata and evidence policy

Record the git commit, compiler/build, BLAS/LAPACK, OpenMP threads, MPI ranks,
CUDA/cuSOLVER, GPU, and CPU fields for every campaign. Also provide matrix
dimension, k-point count, tile size, eigenvector mode, problem type, and
transfer policy when the benchmark knows them. Keep benchmark JSON as evidence
from a controlled run; do not commit machine-local timing results as
correctness references.

## ACC-P0 persistent real-material campaign

ACC-06 remains the legacy cold-process inventory. Its per-repetition command
launches are useful for cold wall time, but they are not steady-state GPU
measurements. ACC-P0 uses the dedicated production driver and keeps the
backend alive while it performs warm-ups and measured repetitions:

```bash
python3 tests/benchmarks/accp0_real_material.py \
  --binary build-acc09-cuda/bin/ReciprocalAccP0Benchmark \
  --build-dir build-acc09-cuda \
  --output-dir results/benchmarks/accp0 \
  --warmups 2 --repetitions 5
```

Use `--quick --skip-cuda` for a CPU-only fixture/oracle smoke run. The full
campaign creates temporary explicit bcc-Fe L=2,3,4,5 production supercells,
validates each against the CPU primitive band-folding oracle, performs a
large-case CUDA memory preflight, and writes `accp0_table.csv` plus
`accp0_results.json`. The table reports cold process wall time separately from
backend initialization, first solve, interval-local H2D/solver/D2H metrics,
steady solve statistics, H(k) CPU assembly, and total steady workload time.

The default ACC-P0 invocation selects `zheevd_serial`; its reset hook clears
only interval timing accumulators. The persistent context, handle, stream,
workspace, and lifetime request counters remain intact. ACC-P1 adds the
explicit strategy option described below.

## ACC-P1/P1b selectable CUDA eigensolver strategies

The same persistent driver accepts explicit `--solver-strategy` values
`fp64_zheevd|fp64_zheevj_batched|fp32_cheevd|fp32_cheevj_batched` (the legacy
`zheevd_serial` and `zheevj_batched` spellings remain aliases). FP64 retains
the original `cuDoubleComplex` path. FP32 converts host H64 into reusable
host `cuComplex` staging before H2D, solves with `cusolverDnCheevd` or
`cusolverDnCheevjBatched`, and widens eigenpairs back to the existing FP64
result arrays. Conversion and widening are reported separately from H2D,
solver, and D2H time. Jacobi strategies are explicitly unsupported for
matrices larger than the cuSOLVER batched API limit (`n > 32`); they never
fall back to CPU or the other precision.

The validation record retains the original H64 matrix and reports H64-based
residual, orthogonality, eigenvalue, and degenerate-projector errors together
with `||H64-H32||/||H64||`. The FP32 campaign is an experimental reciprocal
study only: normal physics arrays and the production FP64 defaults are
unchanged.

The ACC-P1 campaign includes both strategies across the real bccFe L=1…5
supercell set, plus Si and primitive Fe, and emits an `ACCP1_VALIDATION` record
for residual, orthogonality, eigenvalue, and degenerate-projector checks. A
focused five-repetition scaling run is:

```bash
python3 tests/benchmarks/accp0_real_material.py \
  --binary build-acc09-cuda/bin/ReciprocalAccP0Benchmark \
  --build-dir build-acc09-cuda --output-dir /tmp/accp1-supercell \
  --meshes 1 --tiles 1 --warmups 1 --repetitions 5
```

`ReciprocalAccP1bPhysicalSCF` and `accp1b_physical_scf.py` are opt-in physical
probes for the canonical Si and bcc-Fe reciprocal-SCF workflows. They report
canonical electron count, EF, projection-free band energy, site occupation and
charge transfer, Fe moment, DOS state count, near-EF distance, SCF residual,
and iteration count. Fe is run at several Gaussian widths. The probe installs
the selected strategy only on its private reciprocal cache; it is not a
production precision switch. The JSON compares every available CUDA row with
FP64 CPU and also compares FP32 rows with direct FP64 GPU. Its report
intentionally leaves mixed-SCF and delicate response-energy studies out of
scope.

## ACC-P0 supplied bcc-Fe supercells

`accp0_supercell_fe.py` benchmarks the exact explicit supercells under
`example/bulk/supercellFe` and records CPU-versus-CUDA eigenvalue agreement.
It asserts the corrected `crystal_sym='file'` selector so the supplied
`lattice.nml` is consumed; source inputs remain unchanged, while staged
`FeX.nml` files use one canonical `Fe1.nml` potential with relabelled sites.
See
`docs/dev/ACC-P0_SUPERCELL_FE_BENCHMARKS.md` for the input audit and the full
campaign command. The scaling campaign requires at least five measured
repetitions and input `nstep >= 5`.

## ACC-P2 vector-first real-material campaign

ACC-P2 uses reciprocal SCF with `--vectors` as the primary case and covers
primitive Fe, Fe 3^3 (`n=486`), Fe 4^3 (`n=1152`), and Fe 5^3 (`n=2250`). It
records the interval timing budget, separate eigenvalue/eigenvector D2H bytes
and times, and before/after resource counters. FP32 is not selected for
production optimization unless a supported ACC-P1b scientific gate accepts
it.

```bash
python3 tests/benchmarks/accp2_real_material.py \
  --binary build-acc09-cuda/bin/ReciprocalAccP0Benchmark \
  --build-dir build-acc09-cuda \
  --output-dir results/benchmarks/accp2 \
  --warmups 2 --repetitions 5
```

The driver compares persistent pageable and pinned staging for the three
large Fe cases. Pinned staging is backend-owned and enabled only for
`n>=486`; retain it only if end-to-end steady time improves. Use
`--skip-cuda` on a CPU-only host to produce controls without treating missing
CUDA as a scientific or performance pass.
