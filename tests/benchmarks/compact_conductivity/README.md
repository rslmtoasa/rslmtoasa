# Frozen-potential conductivity benchmark

The launcher compares an unchanged baseline executable with a compact legacy
executable using the committed Pt potential. It checks the input and geometry
hashes, the complete raw and physical-unit response grids, the orbital sums
and the type sums. Both jobs explicitly close all three lattice directions.
The total imaginary response of the per-type Pt charge and Hall calculations
must be within (10^{-12}) of the maximum real response. A single-orbital
projection may have a finite imaginary part because its projector cannot be
cycled through the full operator product.

Build both worktrees with the same compiler, flags and BLAS. Set
`BASE_SHA` to the main commit used for the comparison.

```sh
git worktree add --detach ../conductivity-baseline "$BASE_SHA"
cmake -S ../conductivity-baseline -B ../conductivity-baseline/build -DCMAKE_BUILD_TYPE=Release
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build ../conductivity-baseline/build -j4
cmake --build build -j4

python3 tests/benchmarks/compact_conductivity/benchmark.py \
  --baseline ../conductivity-baseline/build/bin/rslmto.x \
  --optimized build/bin/rslmto.x \
  --fixture ../conductivity-baseline/tests/postproc/cases/conductivity/fccPt \
  --output ../conductivity-main-runs --orders 50 100 250 500 \
  --variants spin orbital --repetitions 3 --memory-gib 9.7
```

The script requires NumPy, f90nml and GNU time. Pass the smaller of the
available system and cgroup RAM limits to `--memory-gib`. The guard compares
the lower bounds of the two baseline phases. Leave headroom for lattice,
operator and library allocations. With 2,510 energy points and order 500,
the baseline reconstruction lower bound is 10.56 GiB. Use a host with at least
16 GiB of free RAM for that pair.

The output directory contains `benchmark.csv`, `numerical_parity.csv`,
paired energy grids, per-run `run.json` manifests, profiler lines and
`time.txt`. Each process runs in a new directory. The launcher alternates
execution order and uses single-threaded BLAS with four OpenMP threads by
default. Never run concurrent jobs on the same host while measuring wall time.

Use `--size 8 --channels 160 --orders 20 50 --repetitions 1` for smoke tests.
Use the `hoh_charge`, `hoh_spin`, `hoh_orbital`, `accumulation` and
`random_spin` variants to exercise their respective real-space routes. The
random-vector case uses three vectors and an explicit seed. The `multitype`
variant uses frozen Fe and Co potentials with a two-type B2 unit cell. Its
legacy total is the type average and the physical-unit total is the sum of
the normalized type contributions. The `lehmann_charge` and `dyson_charge`
variants exercise the unchanged full-moment reciprocal reconstruction.
