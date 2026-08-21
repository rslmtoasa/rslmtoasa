# KPM-B0C closure campaign

Strict fair CPU/GPU evidence package. JSON is canonical; this summary is derived from its rows.

## Environment

- Git commit: `434ab35eb3894fa548a8a9998972737adafe43c6`; dirty: `True`
- Compiler: `/usr/bin/gfortran (GNU Fortran (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0)`; build type: `Release`; BLAS/LAPACK: `Intel oneMKL`
- CPU: `11th Gen Intel(R) Core(TM) i9-11900 @ 2.50GHz` (8 physical / 16 logical); RAM: `64095 MiB`
- CUDA: toolkit `13.3`, driver `None`, GPU `None`, selected `None`, VRAM `None MiB`, compute `None`

## Physics/workload and CPU thread winner by group

| material | numeric mode | backend | OMP | transport median (s) | correctness | headline eligible |
|---|---|---|---:|---:|---|---|
| bccFe_magnetic | fp32 | cpu_fast | 1 | 0.50985 | PASS | no |

## Equal-FP64 pairs

| GPU row | CPU reference | S_moments | S_transport | S_whole |
|---|---|---:|---:|---:|

## Equal-FP32 pairs

| GPU row | CPU reference | S_moments | S_transport | S_whole |
|---|---|---:|---:|---:|

## Mixed practical rows

0 mixed rows retained; they are never used for equal-precision headline speedups.

## Failed/ineligible and memory-limited rows


## Correctness summary

PASS: 1

## Definitions

- `S_moments = best_same_precision_CPU(P_moments_total) / GPU(P_moments_total)`.
- `S_transport = best_same_precision_CPU(T_transport_total) / GPU(T_transport_total)`.
- `S_whole = best_same_precision_CPU(whole_wall) / GPU(whole_wall)`.
- `T_transport_total` is the internal transport phase; `whole_wall` is the complete process invocation.
