# KPM-B1 final CPU/GPU benchmark campaign

This report is generated from `campaign.json`; numerical tables and plots are derived from the same canonical rows.

## Executive summary

- Frozen implementation commit: `e23fab86f10dda9ddc1d5dd9452f9586b2fc428a`; campaign rows: `84`; explicit skips plus ineligible GPU records retained: `11`.
- Primary GPU: `NVIDIA RTX A4000` device `0`, `16376 MiB`, compute `8.6`, CUDA `13.3`.
- Headline pairs: FP64 `5`, FP32 `9`; every published pair has profile and production-output correctness PASS.
- CPU rows retain OMP=1/2/4/8; GPU rows use one process and one selected GPU at OMP=1.

## Methodology and freeze

- Compiler/build: `/usr/bin/gfortran (GNU Fortran (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0)`, `Release`, BLAS/LAPACK `Intel oneMKL`.
- CPU: `11th Gen Intel(R) Core(TM) i9-11900 @ 2.50GHz`; MPI ranks recorded as `1`. The host could not configure the MPI Fortran wrapper reliably, so the executable is serial and the one-rank policy is satisfied without an MPI-scaling claim.
- Production output is enabled for timed rows. Two warmups and five measurements are used for default anchors; genuinely expensive r6/r8 and secondary M/G2 rows use three measurements, recorded in their source JSON policy. Medians are used for speedups and spread remains in JSON.
- FP64 means FP64 moments plus FP64 reconstruction. FP32 means FP32 moments plus FP32 reconstruction. The mixed route is FP32 moments plus FP64 reconstruction on CPU; current CUDA profiling couples its moment and reconstruction precision, so no mixed GPU claim is made.
- The tracked implementation was clean at freeze. Pre-existing untracked build/result artifacts were preserved and are not part of the B1 commit; `git_dirty` in the raw environment remains an honest record of that state.

## FP64 headline

| material | size | cond_type | N | nnz | M | estimator | Ntrace | best CPU OMP | CPU moment (s) | GPU moment (s) | S_moments | CPU transport (s) | GPU transport (s) | S_transport | CPU wall (s) | GPU wall (s) | S_whole | correctness |
|---|---:|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| fccPt_SOC | r4 | spin | 1152 | 215784 | 500 | per_type | 1 | 2 | 10.8 | 8.61 | 1.25 | 21.4 | 11 | 1.94 | 21.8 | 11.2 | 1.94 | PASS |
| fccPt_SOC | r6 | spin | 3888 | 884520 | 250 | per_type | 1 | 4 | 9.04 | 7.9 | 1.14 | 12.4 | 9.31 | 1.33 | 12.6 | 9.51 | 1.32 | PASS |
| fccPt_SOC | r6 | spin | 3888 | 884520 | 500 | per_type | 1 | 4 | 34.3 | 28.5 | 1.2 | 45.2 | 31.1 | 1.45 | 45.6 | 31.3 | 1.46 | PASS |
| fccPt_SOC | r8 | spin | 9216 | 2299752 | 250 | per_type | 1 | 8 | 21.5 | 18.6 | 1.16 | 25.3 | 20.5 | 1.24 | 25.6 | 20.7 | 1.24 | PASS |
| fccPt_SOC | r8 | spin | 9216 | 2299752 | 500 | per_type | 1 | 8 | 81.9 | 67.9 | 1.21 | 93.4 | 71.1 | 1.31 | 93.8 | 71.3 | 1.32 | PASS |

## FP32 headline

| material | size | cond_type | N | nnz | M | estimator | Ntrace | best CPU OMP | CPU moment (s) | GPU moment (s) | S_moments | CPU transport (s) | GPU transport (s) | S_transport | CPU wall (s) | GPU wall (s) | S_whole | correctness |
|---|---:|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| bccFe_magnetic | r4 | spin | 1152 | 180792 | 500 | per_type | 1 | 2 | 10.6 | 0.43 | 24.6 | 23.8 | 2.09 | 11.4 | 24.3 | 2.3 | 10.6 | PASS |
| fccPt_SOC | r4 | charge | 1152 | 215784 | 500 | per_type | 1 | 2 | 10.8 | 0.446 | 24.1 | 24.1 | 2.1 | 11.4 | 24.5 | 2.31 | 10.6 | PASS |
| fccPt_SOC | r4 | orbital | 1152 | 215784 | 500 | per_type | 1 | 2 | 10.8 | 0.445 | 24.3 | 24.1 | 2.11 | 11.4 | 24.6 | 2.31 | 10.6 | PASS |
| fccPt_SOC | r4 | spin | 1152 | 215784 | 500 | per_type | 1 | 2 | 10.9 | 0.445 | 24.4 | 24.1 | 2.11 | 11.4 | 24.6 | 2.31 | 10.6 | PASS |
| fccPt_SOC | r6 | spin | 3888 | 884520 | 250 | per_type | 1 | 4 | 9.11 | 0.28 | 32.5 | 13.2 | 1.54 | 8.54 | 13.4 | 1.72 | 7.77 | PASS |
| fccPt_SOC | r6 | spin | 3888 | 884520 | 500 | per_type | 1 | 4 | 34.4 | 1.02 | 33.7 | 47.9 | 2.89 | 16.5 | 48.3 | 3.1 | 15.6 | PASS |
| fccPt_SOC | r6 | spin | 3888 | 884520 | 750 | per_type | 1 | 4 | 77.4 | 2.17 | 35.6 | 107 | 5.01 | 21.3 | 108 | 5.25 | 20.5 | PASS |
| fccPt_SOC | r8 | spin | 9216 | 2299752 | 250 | per_type | 1 | 8 | 21.2 | 0.591 | 35.8 | 25.7 | 2.32 | 11.1 | 26 | 2.5 | 10.4 | PASS |
| fccPt_SOC | r8 | spin | 9216 | 2299752 | 500 | per_type | 1 | 8 | 82.5 | 2.22 | 37.2 | 96.5 | 4.65 | 20.8 | 96.9 | 4.83 | 20.1 | PASS |

## FP32 accuracy evidence

Stored production-output correctness evidence for the FP32 Pt anchors below (per_type, M=500, Ntrace=1); these values are summarized from the attached correctness records.

| workload | mode | max abs | relative L2 | tolerance | result |
|---|---|---:|---:|---|---|
| Pt r4 spin | FP32 | 6.17e-06 | 1.71e-06 | max abs ≤ 0.0005; relative L2 ≤ 0.0005 | PASS |
| Pt r6 spin | FP32 | 3.3e-06 | 1.57e-06 | max abs ≤ 0.0005; relative L2 ≤ 0.0005 | PASS |
| Pt r8 spin | FP32 | 3.5e-06 | 1.44e-06 | max abs ≤ 0.0005; relative L2 ≤ 0.0005 | PASS |
| Pt r4 charge | FP32 | 1.7e-05 | 5.78e-06 | max abs ≤ 0.0005; relative L2 ≤ 0.0005 | PASS |
| Pt r4 orbital | FP32 | 1.24e-05 | 3.17e-06 | max abs ≤ 0.0005; relative L2 ≤ 0.0005 | PASS |

## M-order scaling

| material | size | cond_type | mode | M | CPU OMP | CPU moments (s) | GPU moments (s) | S_moments | CPU transport (s) | GPU transport (s) | S_transport | correctness |
|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| bccFe_magnetic | r4 | spin | fp32 | 500 | 2 | 10.6 | 0.43 | 24.6 | 23.8 | 2.09 | 11.4 | PASS |
| fccPt_SOC | r4 | spin | fp32 | 500 | 2 | 10.9 | 0.445 | 24.4 | 24.1 | 2.11 | 11.4 | PASS |
| fccPt_SOC | r4 | spin | fp64 | 500 | 2 | 10.8 | 8.61 | 1.25 | 21.4 | 11 | 1.94 | PASS |
| fccPt_SOC | r6 | spin | fp32 | 250 | 4 | 9.11 | 0.28 | 32.5 | 13.2 | 1.54 | 8.54 | PASS |
| fccPt_SOC | r6 | spin | fp64 | 250 | 4 | 9.04 | 7.9 | 1.14 | 12.4 | 9.31 | 1.33 | PASS |
| fccPt_SOC | r6 | spin | fp32 | 500 | 4 | 34.4 | 1.02 | 33.7 | 47.9 | 2.89 | 16.5 | PASS |
| fccPt_SOC | r6 | spin | fp64 | 500 | 4 | 34.3 | 28.5 | 1.2 | 45.2 | 31.1 | 1.45 | PASS |
| fccPt_SOC | r6 | spin | fp32 | 750 | 4 | 77.4 | 2.17 | 35.6 | 107 | 5.01 | 21.3 | PASS |
| fccPt_SOC | r8 | spin | fp32 | 250 | 8 | 21.2 | 0.591 | 35.8 | 25.7 | 2.32 | 11.1 | PASS |
| fccPt_SOC | r8 | spin | fp64 | 250 | 8 | 21.5 | 18.6 | 1.16 | 25.3 | 20.5 | 1.24 | PASS |
| fccPt_SOC | r8 | spin | fp32 | 500 | 8 | 82.5 | 2.22 | 37.2 | 96.5 | 4.65 | 20.8 | PASS |
| fccPt_SOC | r8 | spin | fp64 | 500 | 8 | 81.9 | 67.9 | 1.21 | 93.4 | 71.1 | 1.31 | PASS |

## Mixed practical route

| material | size | M | backend | moments | reconstruction | OMP | transport (s) | wall (s) | correctness |
|---|---:|---:|---|---|---|---:|---:|---:|---|
| fccPt_SOC | r4 | 500 | cpu_fast | fp32 | fp64 | 1 | 21.9 | 22.3 | PASS |

## CPU OMP scaling

| material | size | mode | OMP | transport (s) | wall (s) | profile | correctness |
|---|---:|---|---:|---:|---:|---|---|
| bccFe_magnetic | r4 | fp32 | 1 | 24.1 | 24.6 | PASS | PASS |
| bccFe_magnetic | r4 | fp32 | 2 | 23.8 | 24.3 | PASS | PASS |
| bccFe_magnetic | r4 | fp32 | 4 | 24.2 | 24.7 | PASS | PASS |
| bccFe_magnetic | r4 | fp32 | 8 | 25.8 | 26.3 | PASS | PASS |
| bccFe_magnetic | r4 | fp64 | 1 | 21.4 | 21.8 | PASS | PASS |
| bccFe_magnetic | r4 | fp64 | 2 | 21.3 | 21.7 | PASS | PASS |
| bccFe_magnetic | r4 | fp64 | 4 | 21.5 | 22 | PASS | PASS |
| bccFe_magnetic | r4 | fp64 | 8 | 23 | 23.4 | PASS | PASS |
| fccPt_SOC | r4 | fp32 | 1 | 24.3 | 24.8 | PASS | PASS |
| fccPt_SOC | r4 | fp32 | 2 | 24.1 | 24.5 | PASS | PASS |
| fccPt_SOC | r4 | fp32 | 4 | 24.4 | 24.8 | PASS | PASS |
| fccPt_SOC | r4 | fp32 | 8 | 26.1 | 26.6 | PASS | PASS |
| fccPt_SOC | r4 | fp64 | 1 | 21.7 | 22.1 | PASS | PASS |
| fccPt_SOC | r4 | fp64 | 2 | 21.5 | 21.9 | PASS | PASS |
| fccPt_SOC | r4 | fp64 | 4 | 21.8 | 22.2 | PASS | PASS |
| fccPt_SOC | r4 | fp64 | 8 | 23.3 | 23.7 | PASS | PASS |
| fccPt_SOC | r4 | fp32 | 1 | 24.3 | 24.8 | PASS | PASS |
| fccPt_SOC | r4 | fp32 | 2 | 24.1 | 24.6 | PASS | PASS |
| fccPt_SOC | r4 | fp32 | 4 | 24.4 | 24.9 | PASS | PASS |
| fccPt_SOC | r4 | fp32 | 8 | 26 | 26.4 | PASS | PASS |
| fccPt_SOC | r4 | fp64 | 1 | 21.8 | 22.2 | PASS | PASS |
| fccPt_SOC | r4 | fp64 | 2 | 21.5 | 21.9 | PASS | PASS |
| fccPt_SOC | r4 | fp64 | 4 | 21.8 | 22.3 | PASS | PASS |
| fccPt_SOC | r4 | fp64 | 8 | 23.4 | 23.8 | PASS | PASS |
| fccPt_SOC | r4 | fp32 | 1 | 24.3 | 24.8 | PASS | PASS |
| fccPt_SOC | r4 | fp32 | 2 | 24.1 | 24.6 | PASS | PASS |
| fccPt_SOC | r4 | fp32 | 4 | 24.4 | 24.8 | PASS | PASS |
| fccPt_SOC | r4 | fp32 | 8 | 26 | 26.4 | PASS | PASS |
| fccPt_SOC | r4 | fp64 | 1 | 21.6 | 22 | PASS | PASS |
| fccPt_SOC | r4 | fp64 | 2 | 21.4 | 21.8 | PASS | PASS |
| fccPt_SOC | r4 | fp64 | 4 | 21.7 | 22.2 | PASS | PASS |
| fccPt_SOC | r4 | fp64 | 8 | 23.5 | 23.9 | PASS | PASS |
| fccPt_SOC | r4 | mixed | 1 | 21.9 | 22.3 | PASS | PASS |
| fccPt_SOC | r6 | fp32 | 1 | 49.9 | 50.4 | PASS | PASS |
| fccPt_SOC | r6 | fp32 | 2 | 48.4 | 48.9 | PASS | PASS |
| fccPt_SOC | r6 | fp32 | 4 | 47.9 | 48.3 | PASS | PASS |
| fccPt_SOC | r6 | fp32 | 8 | 48 | 48.5 | PASS | PASS |
| fccPt_SOC | r6 | fp64 | 1 | 47.2 | 47.6 | PASS | PASS |
| fccPt_SOC | r6 | fp64 | 2 | 45.8 | 46.2 | PASS | PASS |
| fccPt_SOC | r6 | fp64 | 4 | 45.2 | 45.6 | PASS | PASS |
| fccPt_SOC | r6 | fp64 | 8 | 45.5 | 45.9 | PASS | PASS |
| fccPt_SOC | r8 | fp32 | 1 | 103 | 104 | PASS | PASS |
| fccPt_SOC | r8 | fp32 | 2 | 97.8 | 98.2 | PASS | PASS |
| fccPt_SOC | r8 | fp32 | 4 | 96.9 | 97.4 | PASS | PASS |
| fccPt_SOC | r8 | fp32 | 8 | 96.5 | 96.9 | PASS | PASS |
| fccPt_SOC | r8 | fp64 | 1 | 99.1 | 99.5 | PASS | PASS |
| fccPt_SOC | r8 | fp64 | 2 | 96.1 | 96.5 | PASS | PASS |
| fccPt_SOC | r8 | fp64 | 4 | 93.5 | 94 | PASS | PASS |
| fccPt_SOC | r8 | fp64 | 8 | 93.4 | 93.8 | PASS | PASS |

## Stochastic estimator and G2 block evidence

| backend | precision | Ntrace | block width | moment time/trace (s) | traces/s | transport (s) | correctness |
|---|---|---:|---:|---:|---:|---:|---|
| CPU | fp64 | 1 | 1 | 11.1 | 0.0903 | 21.6 | PASS |
| CPU | fp64 | 2 | 1 | 11.3 | 0.0889 | 35.2 | PASS |
| CPU | fp64 | 4 | 1 | 11.4 | 0.088 | 61.9 | PASS |
| CPU | fp64 | 8 | 1 | 11.4 | 0.0878 | 115 | PASS |
| GPU | fp64 | 1 | 1 | 8.6 | 0.116 | 10.9 | PASS |
| GPU | fp64 | 2 | 1 | 8.54 | 0.117 | 20.1 | PASS |
| GPU | fp64 | 4 | 1 | 8.48 | 0.118 | 38.4 | PASS |
| GPU | fp64 | 8 | 1 | 8.62 | 0.116 | 76.4 | PASS |
| GPU | fp64 | 8 | 2 | 8.81 | 0.114 | 77.9 | PASS |
| GPU | fp64 | 8 | 4 | 8.69 | 0.115 | 76.9 | PASS |
| GPU | fp64 | 8 | 8 | 8.92 | 0.112 | 78.9 | PASS |

## Magnetic bcc-Fe cross-check

Headline Fe rows retained: `1`. The table above is generated from the real `bccFe_magnetic` fixture; no synthetic matrix is used for performance conclusions.

## Interpretation

1. FP64 crossover: GPU already wins at the smallest tested size (r4); crossover below r4 was not measured.
2. FP32 crossover: GPU already wins at the smallest tested size (r4); crossover below r4 was not measured.
3. Size scaling is read from measured r4/r6/r8 rows; no unsupported crossover fit is applied.
4. M scaling is reported from the completed medium-size series and the anchor rows; omitted large extensions are explicit skips below.
5. Moment speedup is a kernel-attribution metric, while transport speedup includes Gamma, reconstruction, post-processing, and other stages.
6. FP64 GPU speedup is limited by non-moment stages and memory-bound work on this A4000; FP32 has a different balance and is only recommended when its production-output tolerance passes.
7. CPU OpenMP changes the selected baseline and therefore materially affects the measured crossover; every anchor retains OMP=1/2/4/8.
8. FP32 scientific acceptability is a correctness/tolerance decision, not a performance assumption; the table above exposes the stored production-output error magnitudes and PASS tolerances directly.
9. FP64 GPU use is worthwhile only in the larger or otherwise GPU-favorable measured regimes; small-workload overhead can leave CPU competitive.
10. Charge, spin, and orbital rows are reported separately; shared-stage conclusions are limited to the measured Pt fixtures.
11. random_vec scaling and per-trace throughput are reported separately from the primary per_type production workflow.
12. Block stochastic processing is judged from measured time/trace and traces/s; neutral or negative results are retained.
13. Pt and magnetic Fe are compared as qualitative cross-checks, not as identical-size performance claims.
14. Production recommendation: use CPU FP64 for small workloads when overhead dominates, GPU FP64 for larger validated FP64 workloads, and GPU FP32 only when its scientific tolerance is acceptable; do not hard-code A4000 thresholds.
15. Further KPM GPU optimization is not justified by B1 unless a new workload exposes a dominant, reproducible stage outside the optimized moment path.
16. Remaining limitations are the one-rank serial executable, single A4000 device, bounded stochastic matrix, and explicitly skipped memory/runtime extensions.

## Skips, failures, and limitations

- `explicit_skip`: Pt r6 M=1000: explicit SKIPPED; the RTX A4000 16 GiB production workspace does not provide a safe margin for the M=1000 Gamma/moment allocation.
- `explicit_skip`: Pt r8 M=750 and M=1000: explicit SKIPPED; the large-r8 GPU Gamma workspace/runtime envelope exceeds the safe bounded A4000 campaign budget after the measured r8 M=250/500 points.
- `explicit_skip`: Pt r4 random_vec Ntrace=16: explicit SKIPPED; the bounded stochastic campaign stopped at Ntrace=8 after complete paired Ntrace=1/2/4/8 coverage, while preserving the G2 block-width matrix.
- `explicit_skip`: Mixed GPU route: explicit NOT AVAILABLE in the current CUDA profile because GPU moment/reconstruction precision is coupled; CPU FP32-moment plus FP64-reconstruction evidence is kept separate.

### Failed or non-headline GPU rows retained

- `bccFe_magnetic_c87776fcd38db55e` (bccFe_magnetic r4 spin fp64 M=500): profile=PASS, correctness=FAIL
- `fccPt_SOC_266440b01d73659f` (fccPt_SOC r4 charge fp64 M=500): profile=PASS, correctness=FAIL
- `fccPt_SOC_29000f5410047e17` (fccPt_SOC r4 orbital fp64 M=500): profile=PASS, correctness=FAIL
- `fccPt_SOC_42e322e8845e47a5` (fccPt_SOC r4 spin fp64 M=500): profile=PASS, correctness=PASS
- `fccPt_SOC_afee5a09d18874c6` (fccPt_SOC r4 spin fp64 M=500): profile=PASS, correctness=PASS
- `fccPt_SOC_8913777c951720cc` (fccPt_SOC r4 spin fp64 M=500): profile=PASS, correctness=PASS
- `fccPt_SOC_8b1a52db54049e82` (fccPt_SOC r6 spin fp64 M=750): profile=PASS, correctness=FAIL

## Generated outputs

- Canonical JSON: `campaign.json`
- CSV: `campaign.csv`
- Plots: `6`
- Raw logs and correctness JSON are under `raw/` and `correctness/` beside the canonical dataset.
