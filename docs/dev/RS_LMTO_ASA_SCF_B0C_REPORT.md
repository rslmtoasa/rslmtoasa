# RS-LMTO-ASA SCF-B0C — CPU/GPU benchmark harness closeout

## 1. Existing infrastructure

SCF-B0C extends the existing ACC/ACC-P infrastructure; it does not replace or
reimplement the reciprocal backend. The production Si and bcc-Fe reciprocal
fixtures, H(R) to H(k) construction, CPU LAPACK reference, CUDA `Zheevd`,
batched Jacobi route, real Fe supercells, and existing reciprocal correctness
drivers remain available.

The new entry point is:

```text
tests/benchmarks/scf_b0c.py
```

It stages a fresh real production fixture for every row, controls OpenMP and
BLAS threads, selects an explicit LAPACK or CUDA route, archives the raw
process output, and derives the canonical result package.

## 2. Gap analysis

| Requirement | Before SCF-B0C | Action |
|---|---|---|
| Real Si/Fe and Fe supercells | ACC/ACC-P drivers existed | Reused the fixtures and added common SCF row metadata |
| Reciprocal H(k), eigensolver, transfer detail | Existing ACC detail timers | Exposed nested `T_H2D`, `T_solver`, `T_D2H`, `T_total_steady` in SCF rows |
| Complete SCF iteration timing | Missing | Added exclusive profile stages around the production SCF loop |
| Full convergence wall time | Driver-specific | Added `full_scf_wall`, convergence residual, final state, and iteration history |
| CPU OMP/BLAS fairness | Existing in parts | Canonical OMP 1/2/4/8 sweep with controlled BLAS threads and affinity |
| Precision identity | Solver-specific | Records five precision fields and classifies `fp64`/`mixed` honestly |
| Strict physical pairing | Separate mechanisms | Added stable physics/start-state/numeric comparison keys |
| Headline eligibility | Informal | Requires matching physics, correctness, profile closure, CUDA route, and mode |
| JSON/CSV/Markdown/raw evidence | Existing outputs varied | One JSON dataset derives CSV, Markdown, raw logs, correctness, and histories |

## 3. SCF profiling model

The profile is an observation layer around the existing production code. It
does not alter Hamiltonian construction, occupations, density construction,
mixing, convergence thresholds, or potential equations.

Each emitted `SCF_B0C_ITER` row contains the exclusive decomposition:

```text
P_scf_iteration_total
  = P_hamiltonian_prepare
  + P_hk_assembly
  + P_eigensolver
  + P_eigenpair_transfer
  + P_occupations_fermi
  + P_density_build
  + P_charge_spin_accumulate
  + P_potential_update
  + P_mixing
  + P_scf_io
  + P_scf_misc
```

The outer total is measured independently and any remaining production
boundary is assigned to `P_scf_misc`. Closure is required to be at most 3%;
steady-state `P_scf_misc` is required to remain below 5%. The first iteration
is retained separately because cold initialization and the first reciprocal
request are not steady-state work.

CUDA device timing is nested detail, not added to the exclusive parent sum.
The GPU solver speedup uses `T_solver`, while the CPU reference uses the
exclusive `P_eigensolver` stage.

## 4. Pairing and correctness contract

`build_comparison_key()` includes material, fixture revision, supercell,
runtime matrix dimension, spin/SOC/basis state, structure-constant backend,
nominal and actual k sampling, smearing/temperature, electron count, Fermi
policy, mixing identity, convergence threshold, starting state, potential
identity, feature flags, eigenvector requirement, and all five precision
fields. The fingerprint is stable JSON hashed data.

Electron count is rounded only for the pairing fingerprint so harmless
floating-point summation noise cannot split otherwise identical CPU/GPU rows;
the unrounded charge remains part of correctness evidence.

Converged rows compare final total energy, Fermi energy, total charge, site
charge transfer, site magnetic moment, and final residual. Different iteration
counts are reported but are not an automatic failure. Eigensystem correctness
continues to be supplied by the existing ACC-P validation paths; raw
eigenvector phase comparison is not introduced here.

The CUDA FP32 eigensolver routes use FP64 Hamiltonian/density/canonical SCF
storage with FP32 solver/eigenvector work, so they are classified as `mixed`.
A full end-to-end FP32 SCF route is not fabricated. Mixed rows remain visible,
are checked against the FP64 CPU oracle, and cannot become equal-precision
headlines.

## 5. Material fixtures

The canonical driver supports:

```text
si   -> validated diamond-Si primitive fixture
fe   -> validated bcc-Fe primitive fixture
fe2  -> production bcc-Fe 2x2x2 fixture
fe3  -> production bcc-Fe 3x3x3 fixture
```

The existing ACC-P drivers retain the 4x4x4 and 5x5x5 frozen-potential
reciprocal/eigensolver scaling cases. Runtime dimensions and unique k counts
are always recorded from the executable rather than assumed by the harness.

## 6. Output and provenance

For an output such as `results/benchmarks/scf_b0c/campaign.json`, the harness
writes:

```text
campaign.json          full-fidelity canonical dataset
campaign.csv           one summarized row per configuration
campaign.md            concise closure and eligibility table
raw/*.log              exact fresh-process stdout/stderr
correctness/*.json     final-state/profile evidence per row
iteration_history/*.json
                         compact scalar history per row
```

The dataset records compiler/build flags, BLAS/LAPACK, MPI status, CPU model,
thread controls, OS/kernel, CUDA toolkit/driver, GPU model, device, VRAM, and
compute capability. CUDA rows additionally retain H2D/D2H bytes, host
conversion/staging/widening, synchronization, workspace query/reuse counts,
pinned-host state, and the point-in-time free/total device-memory snapshot.
Optional metadata remains nullable.

The primary command is:

```bash
python3 tests/benchmarks/scf_b0c.py \
  --binary build-acc09-cuda/bin/ReciprocalAccP1bPhysicalSCF \
  --build-dir build-acc09-cuda \
  --output results/benchmarks/scf_b0c/campaign.json \
  --materials si,fe \
  --strategies fp64_zheevd,fp64_zheevj_batched,fp32_cheevd,fp32_cheevj_batched \
  --omp-threads 1,2,4,8 --dos-method gaussian --nstep 40
```

`--dos-method gaussian|tetrahedron|blochl` is part of the strict pairing key.
CUDA initialization or unsupported strategy failures are retained as
`UNSUPPORTED`; they are never converted into CPU timings.

## 7. Closure evidence

The stored CPU closure package is:

```text
results/benchmarks/scf_b0c_cpu/campaign.json
```

It contains 8 rows: Si and Fe, each at OMP 1/2/4/8. All rows converged with
profile closure passing. Representative steady SCF iteration medians were
approximately 0.126 s for Si and 0.140 s for Fe; full process walls were
approximately 4.0 s and 5.9 s respectively on the recorded host.

The host CUDA package is:

```text
results/benchmarks/scf_b0c/campaign.json
```

It contains 16 rows: 8 CPU controls and 8 CUDA strategy rows on an NVIDIA RTX
A4000 (CUDA 13.3, driver 610.57.04). All four FP64 GPU rows passed profile and
final-state correctness and are eligible for separate solver/reciprocal,
iteration, and convergence ratios. In this small primitive workload, the
ratios are below one for the GPU, which is consistent with the existing
small-matrix ACC evidence and is not a new optimization conclusion.

The four mixed rows are retained with their nested solver timings but are not
headline eligible: the normal FP64 SCF canonical path remains the scientific
oracle, and the mixed routes did not reach the normal convergence criterion in
the 40-step closure budget. A longer 80-step Fe mixed probe also remained
non-converged, so these routes are correctly treated as experimental evidence.

The Fe `2x2x2` steady-iteration smoke package is:

```text
results/benchmarks/scf_b0c_fe2/campaign.json
```

It exercises the real `nmat=144`, `Nk_unique=65` production case for two
iterations. Profile closure passes with an approximately 2.3 s steady median.
Its residual
is intentionally not presented as a converged SCF result; this is a larger
fixture integration/steady-profile check. Full Fe supercell scaling remains
covered by the existing ACC-P corpus.

The production DOS-kernel selector was also exercised with CPU Si closure
packages for `tetrahedron` and `blochl`. Both converged with profile closure
passing; because these routes use the full `8x8x8` mesh their actual
`Nk_unique=512` is retained in the rows rather than being conflated with the
symmetry-reduced Gaussian case.

## 8. Remaining limitations

SCF-B0C is harness validation, not SCF-B1 application policy. In particular:

- the closure campaign is single-rank and single-GPU; it makes no MPI or
  multi-GPU scaling claim;
- the physical probe is a fresh process per row, so the recorded cold process
  wall includes CUDA context/backend initialization; separate backend-init
  counters remain an ACC reciprocal detail follow-up;
- full end-to-end FP32 SCF is unsupported and is not implied by FP32 CUDA
  eigensolver rows;
- the largest Fe supercells remain appropriate for frozen-potential reciprocal
  scaling rather than mandatory full SCF convergence;
- DOS-kernel variants are selectable and strictly paired, but the recorded
  closure evidence above uses the Gaussian route.

## 9. Next task

The harness is ready for:

```text
SCF-B1 — frozen final CPU/GPU SCF benchmark campaign
```

SCF-B1 should freeze the executable, fixture revisions, DOS kernel, solver
matrix, repetitions, and acceptance tolerances before making application-level
performance statements.
