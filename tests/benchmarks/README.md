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
