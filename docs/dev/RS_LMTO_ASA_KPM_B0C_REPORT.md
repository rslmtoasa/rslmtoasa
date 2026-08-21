# KPM-B0C — Fair CPU/GPU benchmark harness closure

## Audit of CURRENT HEAD

| B0 requirement | already present | partial | missing | action in B0C |
|---|---|---|---|---|
| benchmark dimensions | G1.2/G2 profile metadata | material was Pt-specific | canonical fixture contract | generalized Pt/Fe row schema and key |
| OMP 1/2/4/8 | G1.2 sweep controls | winner was not contract-bound | strict best-CPU selection | retain all rows; select only valid same-mode rows |
| BLAS control | environment controls | no headline validator | — | record and exclude from pairing key |
| GPU metadata | basic model/toolkit | selected device/VRAM/compute capability incomplete | fuller hardware record | selected-device `nvidia-smi` capture |
| warmups/repetitions | G1.2 timing loops | logs were overwritten | stable sample IDs | archive every warmup/sample log |
| same-precision pairing | precision fields | grouping was permissive | canonical fingerprint | implement numeric taxonomy and strict validator |
| same stochastic states | fixed seed in G2 | seed was annotation-only | enforcement | seed/Ntrace/vector contract in key and validator |
| stage timing | exclusive profile and closure | no canonical row output | — | preserve profile phases and summarize them |
| correctness per row | VAL-09 production-output philosophy | placeholder status remained | structured attached evidence | compare retained production outputs |
| CPU-best reference | speedup helper existed | grouped by too few fields | valid pair-only reference | best OMP row selected inside valid group |
| CSV | missing | — | CSV dataset | derive from canonical rows |
| JSON | G1.2 JSON existed | schema was campaign-specific | B0C full-fidelity package | canonical B0C JSON |
| Markdown | missing | — | Markdown summary | derive from canonical rows |
| raw logs | one scratch `testrun.log` | overwritten per sample | archive | stable row/sample log paths |
| git/build/hardware metadata | core metadata | dirty state and hardware details incomplete | fuller provenance | conservative null-safe capture |

No second benchmark framework was introduced. The existing G1.2/G2 driver
remains the entry point, with its filename retained for compatibility and its
workflow consolidated into the B0C schema.

## What already existed

G1.2/G2 already supplied real SOC-Pt staging, production-scale defaults,
charge/spin/orbital selection, `per_type` and `random_vec` estimators, fixed
random seeds, CUDA stochastic block widths and memory skips, CPU backend and
thread sweeps, exclusive KPM phase closure, resident GPU moments, tiled Gamma,
cuBLAS reconstruction, reduced-result transfers, no-write diagnostics, and
median/min/max/MAD/IQR timing. The VAL-09 harness already established that
production conductivity files are the correctness authority.

## What B0C added

- A five-field precision taxonomy: moment backend/precision,
  reconstruction backend/precision, canonical output precision (`fp64`), and derived
  `numeric_mode` (`fp64`, `fp32`, or `mixed`) in the campaign rows.
- A narrow opt-in `cpu_reconstruction_precision='fp32'` route. It computes
  FP32 Gamma and invokes the existing CPU CGEMM helper with FP32 packed
  moments/result, widening only the reduced result to canonical host FP64.
  The default `cpu_fast` FP32-moment/FP64-reconstruction route is unchanged
  and is labelled `mixed`.
- A complete physics/numerics comparison key and stable hash, strict pair
  rejection, correctness-gated best CPU selection, and no-write exclusion.
- Structured production-output correctness evidence, stable row/sample IDs,
  archived raw logs, output checksums, and JSON/CSV/Markdown generation from
  one row dataset.
- Null-safe git, build, compiler, BLAS, CPU, RAM, OS, CUDA, driver, selected
  GPU, VRAM, compute-capability, affinity, and visibility metadata.
- One material contract for validated `fccPt_SOC` and `bccFe_magnetic`
  fixtures, with Fe smoke selection through the same workflow.

## Precision contract

`fp64` means `moment_precision=fp64` and
`reconstruction_precision=fp64`. `fp32` means both are `fp32`. Every other
combination is `mixed`, including the historical CPU fast route. Canonical
host/output arrays may still be FP64 after a FP32 reduced result is widened.
Mixed rows are retained for practical evidence but cannot produce an
equal-precision headline speedup.

## Pairing contract

The key includes material, fixture ID/revision, replication, actual `N` and
`nnz`, `nsp`, SOC state, current operator (`va`, `vb`, and spin/orbital
selector), estimator and projector/trace contract, `Ntrace`, random seed or
deterministic projector identifier, `M`, `lld`, `NE`/channels, `rc`, energy
window, Fermi level, Lorentz kernel, Chebyshev scaling, and `numeric_mode`.
OMP/BLAS counts, backend implementation, GPU strategy, and stochastic block
width are strategy metadata rather than physics-key fields.

A headline pair requires identical physics key, equal mode in `{fp32, fp64}`,
`correctness.status=PASS`, `profile_status=PASS`, and production output mode.
Random-vector pairs additionally require identical seed, `Ntrace`, and vector
distribution/normalization contract. Rejections retain explicit reasons such
as `precision_mismatch`, `physics_mismatch`, `seed_mismatch`,
`correctness_failed`, `profile_failed`, and `benchmark_no_write`.

## Correctness contract

The harness compares the actual production `*cond*.out` files, following
VAL-09, rather than reimplementing Kubo–Bastin in Python. It records max
absolute error, relative L2 error, a selected integrated/tensor metric, the
compared file names, tolerance set, and an evidence ID. The established
envelopes are FP64 `(max_abs 2e-6, rel_l2 5e-6, integrated 5e-6)` and FP32
`5e-4` for all three metrics, based on the existing G1/G1.3/VAL-09 evidence
and the `ES16.6` production-file output quantum. The FP64 envelope is two
printed ulps, not a relaxation of the underlying arithmetic comparison.
Optimized GPU timing keeps full-moment D2H disabled; correctness uses retained
production outputs and references diagnostic moment evidence when available.

## Material support

`--material pt` selects the validated SOC Pt fixture and `--material fe`
selects the existing magnetic bcc-Fe triad fixture. Both carry explicit
fixture identity, revision, `nsp`, SOC/magnetic state, operator defaults, and
cutoff semantics in each row.

## Environment and output package

Example metadata from the CUDA closure host is captured as selected device
`0`, NVIDIA RTX A4000, 16,376 MiB VRAM, compute capability 8.6, CUDA driver
610.57.04, with the current git commit and dirty state. Unavailable optional
fields are represented as `null` with an explanatory optional-metadata map.

The package layout is:

```text
results/benchmarks/kpm_b0c/
  campaign.json
  campaign.csv
  campaign.md
  raw/<row_id>/warmup_01.log
  raw/<row_id>/sample_01.log
  raw/<row_id>/outputs/*.out
  correctness/<pair_id>/comparison.json
```

JSON contains full metadata, samples, statistics, correctness, pairing,
speedup eligibility, skipped rows, and provenance. CSV contains one summary
row per configuration. Markdown contains environment, CPU winners, equal
FP64/FP32 tables, mixed rows, failures/skips, correctness, and speedup
definitions.

## Closure evidence

The committed `campaign.json` is the manageable B0C closure package: Pt r4,
`M=500`, `NE=2510`, `lld=150`, spin, `per_type`, CPU OMP 1/2/4/8, CPU FP64,
CPU FP32, the existing mixed CPU route, GPU FP64, and GPU FP32, with two
warmups and three measured repetitions. A separate controlled random-vector
pair and bcc-Fe smoke run are recorded alongside the harness validation.

These are harness-closure results, not the final publication campaign. The
reported `S_moments`, `S_transport`, and `S_whole` values are emitted only for
strict PASS pairs and are never inferred from nearby rows.

For this Pt r4 closure, the median across the four GPU OMP rows was
`S_moments/S_transport/S_whole = 24.78/11.03/10.13` for FP32 and
`1.265/1.942/1.933` for FP64. All 20 rows had profile and correctness PASS;
the existing mixed CPU route remains in the package but has no headline
equal-precision speedup.

## Remaining work

KPM-B1 is next. No further benchmark-harness redesign is expected before B1
unless B0C finds a correctness or pairing defect.
