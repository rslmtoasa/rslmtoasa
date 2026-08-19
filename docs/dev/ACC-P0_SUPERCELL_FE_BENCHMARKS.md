# ACC-P0 explicit bcc-Fe supercell benchmarks

`tests/benchmarks/accp0_supercell_fe.py` turns the supplied
`example/bulk/supercellFe/{2x2x2,3x3x3,4x4x4,5x5x5}` directories into a
machine-readable ACC-P0 campaign. The source directories are copied to a
temporary run directory; their inputs and lattice files are not rewritten.
For a valid uniform-potential comparison, the staged copies replace every
`FeX.nml` with the canonical `2x2x2/Fe1.nml` contents and change only its
site label. The source directories remain available for audit. All cases keep
the supplied `lmax=2` Fe basis. The k-point workload is scaled geometrically
from a 32x32x32 primitive-cell option:

| fixture | sites | mesh | k-points | complex matrix |
|---|---:|---:|---:|---:|
| 2x2x2 | 8 | 16x16x16 | 4096 | 144 |
| 3x3x3 | 27 | 8x8x8 | 512 | 486 |
| 4x4x4 | 64 | 4x4x4 | 64 | 1152 |
| 5x5x5 | 125 | 2x2x2 | 8 | 2250 |

The corrected supplied `input.nml` files select `crystal_sym='file'`, while
their explicit `L^3` sites live in `lattice.nml`. In the current production
lattice path, `bcc` would reconstruct the one-site primitive geometry and
silently report `nrec=1`; the benchmark therefore asserts `file` and rejects
a fixture that regresses to `bcc`. The source input is copied unchanged into
the temporary run directory; only staged potential files are canonicalized.

Run the full CPU/CUDA campaign on the CUDA host with, for example:

```bash
python3 tests/benchmarks/accp0_supercell_fe.py \
  --binary build-acc09-cuda/bin/ReciprocalAccP0Benchmark \
  --build-dir build-acc09-cuda \
  --output-dir /tmp/accp0-supercellFe \
  --warmups 2 --repetitions 5 --tiles 1
```

The output contains a CSV timing table, raw process logs/eigenvalue dumps,
and `supercellFe_accp0_results.json`. The driver keeps one process alive for
each row, initializes one backend, performs in-process warm-ups, resets CUDA
interval counters, and measures at least five repeated H(k)+eigensolve work
steps. The source `nstep` is also required to be at least five. CUDA is not
claimed to win automatically; compare the per-fixture `total_steady_s` values
and the reported transfer/solver components.

The corrected-input, uniform-staged five-repetition run on the CUDA host produced:

| fixture | mesh | matrix | CPU steady (s) | CUDA steady (s) | CPU/CUDA |
|---|---:|---:|---:|---:|---:|
| 2x2x2 | 16x16x16 | 144 | 6.503 | 23.061 | 0.282x |
| 3x3x3 | 8x8x8 | 486 | 9.705 | 11.297 | 0.859x |
| 4x4x4 | 4x4x4 | 1152 | 10.345 | 6.511 | 1.589x |
| 5x5x5 | 2x2x2 | 2250 | 10.539 | 3.831 | 2.751x |

These are persistent k-space timings for the current one-matrix-at-a-time
cuSOLVER path with five measured repetitions and `nstep=5`; they do not
justify making CUDA the unconditional default.

Correctness is reported separately. For every selected supercell, CPU LAPACK
and CUDA solve the same assembled Hamiltonian and their sorted eigenvalue
multisets are compared. A second CPU-only check compares the staged uniform
supercell spectrum against a primitive reference sampled on the commensurate
`LxLxL` mesh. This isolates explicit-geometry/assembler folding. The source
audit records the original 5, 6, 6, and 8 normalized potential contents, while
the measured staging audit reports one normalized potential for every fixture.
