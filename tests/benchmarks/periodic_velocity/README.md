# Periodic velocity correction and validation

The matrix-free velocity on main `3957b2c` used the stored coordinate difference
for every neighbor. A periodic crossing therefore acquired a supercell-length
factor instead of the physical bond vector. The fix caches
`d_ij = -(wrapped r_j-r_i)` in physical units, then uses `-i direction.d_ij`
in all six occurrences in `velo_vec_matmul` and `velo_hoh_vec_matmul`.

## Image convention and scope

`nncal` admits one entry per atom pair. `remd` reorders neighbors by the wrapped
vector of a representative site, and `chbar_nc` uses the same
`f_wrap_coord_diff` vector when assigning the bulk/local hopping blocks.
The cache uses that existing selector, actual site/neighbor-slot indices, and
coordinates multiplied by `alat`. `pbc` gates wrapping; `b1:b3` independently
enable the translations `n1*a(:,1)*alat`, etc. Open boundaries retain the raw
coordinate difference. The helper's j-i result is negated for velocity.

The existing selector searches shifts -1,0,1 on enabled axes; this patch does
not expand the Hamiltonian's supported image search or invent image labels.
Equal shortest images are diagnosed and rejected by transport setup, with an
instruction to enlarge the supercell. The current single-entry map cannot
represent distinct hoppings to multiple images of the same atom. Merely
choosing a deterministic tie-break would not repair that Hamiltonian ambiguity.
The optional diagnostic leaves all existing geometry/Hamiltonian calls and
selected vectors unchanged.

The cache is `3 x size(nn,2) x kk` real numbers, not a dense site-pair array.
One temporary `3 x kk` physical-coordinate array and the image searches occur
only during setup. Matvec loops only read the cache. Existing explicit
local/bulk loops, `izero`/`idum`, OpenMP row ownership, slot-b endpoint weighting,
slot-a unweighted response, current anticommutators, both HOH terms
`v-vOh-hOv`, H/SOC, shifts/scaling, and integration are preserved. No new
Fortran module or input is introduced. Direct kernel callers initialize lazily;
after a geometry/map/flag change call `setup_kubo_bonds` (or operator setup)
again before applying operators. Concurrent geometry mutation is unsupported.

## Validation

The production driver applies every operator to a complete 54-dimensional
identity basis. NumPy independently constructs the Hermitian hopping matrix,
local impurity overrides, overlap matrices, and analytic directed image vectors.
The periodic reference does not call the wrapping helper or perform a copied
minimum-image search. The equivalent-representation fixture shifts one site
by one supercell vector within the existing finite image-search domain; it
does not claim invariance for arbitrary coordinates outside that domain.
It covers all Cartesian velocity directions, both weights,
both Hamiltonian representations, and the full existing operator/slot suite,
including Sy, Ly, quadrupoles and octupoles. Each fixture compares 1,528 full
matrices. The synthetic chain has an initially absent 1--3 hopping, so the HOH
test necessarily exercises two-hop contributions.

| Requirement | Result |
|---|---|
| A: open/disabled, individual axes, all axes, skew cell, crossing/noncrossing | Passed; cached vectors agree with analytic signs and magnitudes |
| B: velocity Hermiticity, reverse bonds, weighted/local blocks | Passed in x/y/z and both representations |
| C: full dense operators and spin/orbital currents | Passed; maximum absolute error 3.28e-15 |
| D: `v-vOh-hOv`, impurity and two-hop terms | Passed |
| E: compiled pinned main vs corrected open-boundary outputs | Passed; maximum absolute difference 1.78e-15 |
| F: origin shift and equivalent periodic representation | Passed against the same analytic reference |
| Degenerate images | Tie detected and rejected; nearby unique image accepted |
| Existing multipole regressions | Passed |
| G: periodic fccPt | See results below |

The regression suite reports 38 tests and 15 subtests passed. Maximum relative
Frobenius operator error is 1.99e-15. `results/operator_errors.csv` records each
fixture; `open_main_errors.csv` compares 4,455,648 complex elements per fixture
against the compiled, unmodified main kernels (four fixtures).

## End-to-end runs and performance

Frozen example potentials; 12x12x12 = 1,728 sites; all wrapping flags true;
SOC `nsp=2`; HOH false; output y/input x; z polarization; per-type trace;
all other example parameters unchanged. Versions: current main `3957b2c`,
legacy `e604f40`, and corrected branch. The legacy outputs were regenerated
locally with the same compiler/inputs, rather than imported from another host.

Environment: AMD EPYC 9V74; 8-core CPU quota; 8 GiB memory limit; gfortran 13.3;
Release `-O3 -funroll-loops -finline-functions -march=native`; reference BLAS
and LAPACK; four OpenMP threads; BLAS thread controls set to one.
`system_clock` measures moments and cache setup; `wait4` measures total elapsed
and peak RSS. Total time includes dumping the raw moments; the moments timer
excludes that dump. The one-time cache timer is included in moments time.
Timings are single-run observations, not a statistically established speedup.
Some reduced runs overlapped compilation/testing; use the isolated final runs
for the clearest runtime comparison.

All spin/orbital runs at cond_ll=50 and 100 completed successfully, with finite
moments and complete per-type curves. The isolated orbital-100 measurement is
105.756 s total / 82.616 s moments on main versus 105.619 s / 81.657 s corrected.
Peak RSS is 1,694,668 versus 1,695,892 KiB: +1.20 MiB (+0.072%). Cache setup
is 0.01256 s. This shows no material runtime regression in the measured case;
it does not establish a general speedup.

| cond_ll | Output | Main total / moments (s) | Corrected total / moments (s) | Corrected cache (s) |
|---:|---|---:|---:|---:|
| 50 | Spin | 58.43 / 31.09 | 51.71 / 29.58 | 0.0130 |
| 50 | Orbital | 53.60 / 29.95 | 49.48 / 25.99 | 0.0138 |
| 100 | Spin | 110.01 / 84.47 | 97.69 / 74.44 | 0.0134 |
| 100 | Orbital | 105.76 / 82.62 | 105.62 / 81.66 | 0.0126 |

`results/timings.csv` contains all versions, peak memory, and CPU times.
At cond_ll=50, corrected-vs-main relative moment differences are 0.306 (spin)
and 0.140 (orbital); at 100 they are 0.851 and 0.703. The real conductivity
curve differences are much smaller: relative L2 0.000682/0.000104 at 50,
and 0.01340/0.002387 at 100. Tracing and integration can conceal large
operator/moment errors. Main's nonzero imaginary curves disappear to roundoff
in the corrected bulk fixture. The corrected and regenerated legacy moments
agree to roundoff here, and the real curves agree at their printed precision.
That agreement is an observation, not an acceptance requirement or proof:
the independently constructed dense periodic operator is the reference.


`results/differences.csv` reports maximum absolute and relative-L2 differences
for full moment tensors and real/imaginary conductivity curves.
`curves_*.csv` supplies the sampled curves. Relative errors of near-zero
imaginary curves are not meaningful by themselves; use absolute errors.
For this single-type fixture, `Pt_cond.out` equals `cond_total.out` wherever
both have rows. We use the explicitly closed per-type file because the legacy
spin-50 total file ended early (2,016 of 2,510 rows), while its per-type file was
complete. No output/integration code was changed to mask that existing issue.

The cond_ll=500 attempt was stopped after identifying a hard memory limit:
`gamma_nm(2510,500,500)` alone requires 10,040,000,000 bytes (9.35 GiB), and
`mu_nm_stochastic(18,18,500,500,1)` requires another 1,296,000,000 bytes.
Their combined 10.56 GiB exceeds this container's 8 GiB before other arrays.
No cond_ll=500 result or runtime is claimed. The supplied runner checks this
lower bound before starting such a run. A larger-memory host is required;
streaming the integration tensor would be a separate architecture change.

## Reproduce

Install CMake, a Fortran compiler, BLAS/LAPACK, Python NumPy, pytest,
pytest-subtests and f90nml. From the repository root:

```sh
export FC=gfortran
export CMAKE=cmake
export OPENBLAS_NUM_THREADS=1
# Optional explicit paths on systems without development symlinks:
export RS_LMTO_TEST_CMAKE_ARGS='-DBLAS_LIBRARIES=/usr/lib/x86_64-linux-gnu/libblas.so.3 -DLAPACK_LIBRARIES=/usr/lib/x86_64-linux-gnu/liblapack.so.3;/usr/lib/x86_64-linux-gnu/libblas.so.3'
python -m pytest -q --basetemp=/tmp/periodic-tests \
  tests/regression/test_hoh_transport.py \
  tests/regression/test_octupole_operators.py \
  tests/regression/test_quadrupole_operators.py
python tests/benchmarks/periodic_velocity/benchmark.py \
  --workdir /tmp/periodic-benchmark --levels 50 100 --repeats 1
python tests/benchmarks/periodic_velocity/compare.py \
  /tmp/periodic-benchmark --output /tmp/periodic-comparison
```

The benchmark creates detached worktrees and adds identical timing/dump
instrumentation only in those worktrees. `--fixed-ref` defaults to HEAD; commit
the correction before running it. On a machine with sufficient RAM, add 500
to `--levels`. Use `--repeats 3` for stronger timing evidence; ordering reverses
on even repeats. No SCF step is run. New destination directories are required.

For the independent main/open comparison, create a worktree at `3957b2c`,
configure Debug with `BUILD_TESTING=ON`, OpenMP enabled, MPI disabled, and build
`transport_operator_test` with the same compiler/library settings. Then:

```sh
python tests/benchmarks/periodic_velocity/compare_open_main.py \
  /path/to/main-build/bin/transport_operator_test /tmp/periodic-tests \
  /tmp/open_main_errors.csv
```

Each pytest fixture also writes its raw matrix-error CSV, bond cache, input,
and full action output into the selected basetemp directory. The committed
summaries retain errors and curves; large binary moment tensors are reproducible
with the benchmark and are not committed.
