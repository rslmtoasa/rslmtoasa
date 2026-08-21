# KPM-B0C closure campaign

Strict fair CPU/GPU evidence package. JSON is canonical; this summary is derived from its rows.

## Environment

- Git commit: `e23fab86f10dda9ddc1d5dd9452f9586b2fc428a`; dirty: `True`
- Compiler: `/usr/bin/gfortran (GNU Fortran (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0)`; build type: `Release`; BLAS/LAPACK: `Intel oneMKL`
- CPU: `11th Gen Intel(R) Core(TM) i9-11900 @ 2.50GHz` (8 physical / 16 logical); RAM: `64095 MiB`
- CUDA: toolkit `13.3`, driver `610.57.04`, GPU `NVIDIA RTX A4000`, selected `0`, VRAM `16376 MiB`, compute `8.6`

## Physics/workload and CPU thread winner by group

| material | numeric mode | backend | OMP | transport median (s) | correctness | headline eligible |
|---|---|---|---:|---:|---|---|
| fccPt_SOC | fp32 | cpu_fast | 1 | 24.342 | PASS | no |
| fccPt_SOC | fp32 | cpu_fast | 2 | 24.088 | PASS | no |
| fccPt_SOC | fp32 | cpu_fast | 4 | 24.425 | PASS | no |
| fccPt_SOC | fp32 | cpu_fast | 8 | 26.013 | PASS | no |

## Equal-FP64 pairs

| GPU row | CPU reference | S_moments | S_transport | S_whole |
|---|---|---:|---:|---:|

## Equal-FP32 pairs

| GPU row | CPU reference | S_moments | S_transport | S_whole |
|---|---|---:|---:|---:|
| fccPt_SOC_e751a23f725589bf | fccPt_SOC_34c3b7f828e1020e | 24.261305860099846 | 11.439648771409674 | 10.625346577400128 |

## Mixed practical rows

0 mixed rows retained; they are never used for equal-precision headline speedups.

## Failed/ineligible and memory-limited rows


## Correctness summary

PASS: 5

## Definitions

- `S_moments = best_same_precision_CPU(P_moments_total) / GPU(P_moments_total)`.
- `S_transport = best_same_precision_CPU(T_transport_total) / GPU(T_transport_total)`.
- `S_whole = best_same_precision_CPU(whole_wall) / GPU(whole_wall)`.
- `T_transport_total` is the internal transport phase; `whole_wall` is the complete process invocation.
