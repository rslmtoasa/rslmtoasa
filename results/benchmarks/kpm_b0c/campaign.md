# KPM-B0C closure campaign

Strict fair CPU/GPU evidence package. JSON is canonical; this summary is derived from its rows.

## Environment

- Git commit: `434ab35eb3894fa548a8a9998972737adafe43c6`; dirty: `True`
- Compiler: `/usr/bin/gfortran (GNU Fortran (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0)`; build type: `Release`; BLAS/LAPACK: `Intel oneMKL`
- CPU: `11th Gen Intel(R) Core(TM) i9-11900 @ 2.50GHz` (8 physical / 16 logical); RAM: `64095 MiB`
- CUDA: toolkit `13.3`, driver `610.57.04`, GPU `NVIDIA RTX A4000`, selected `0`, VRAM `16376 MiB`, compute `8.6`

## Physics/workload and CPU thread winner by group

| material | numeric mode | backend | OMP | transport median (s) | correctness | headline eligible |
|---|---|---|---:|---:|---|---|
| fccPt_SOC | mixed | cpu_fast | 1 | 21.658 | PASS | no |
| fccPt_SOC | mixed | cpu_fast | 2 | 21.492 | PASS | no |
| fccPt_SOC | mixed | cpu_fast | 4 | 21.653 | PASS | no |
| fccPt_SOC | mixed | cpu_fast | 8 | 21.546 | PASS | no |
| fccPt_SOC | fp32 | cpu_fast | 1 | 24.271 | PASS | no |
| fccPt_SOC | fp32 | cpu_fast | 2 | 24.355 | PASS | no |
| fccPt_SOC | fp32 | cpu_fast | 4 | 24.378 | PASS | no |
| fccPt_SOC | fp32 | cpu_fast | 8 | 24.389 | PASS | no |
| fccPt_SOC | fp64 | cpu_fast | 1 | 21.815 | PASS | no |
| fccPt_SOC | fp64 | cpu_fast | 2 | 21.775 | PASS | no |
| fccPt_SOC | fp64 | cpu_fast | 4 | 21.755 | PASS | no |
| fccPt_SOC | fp64 | cpu_fast | 8 | 21.668 | PASS | no |

## Equal-FP64 pairs

| GPU row | CPU reference | S_moments | S_transport | S_whole |
|---|---|---:|---:|---:|
| fccPt_SOC_56d6e9b9855b19a3 | fccPt_SOC_445821d99afc81fa | 1.2682166349115416 | 1.950479477724196 | 1.941006855790763 |
| fccPt_SOC_8f086f2b5b5acf6b | fccPt_SOC_445821d99afc81fa | 1.264462504613966 | 1.9390112683117657 | 1.9305889687787599 |
| fccPt_SOC_718d626135e5981b | fccPt_SOC_445821d99afc81fa | 1.2650396586831465 | 1.9434056796238226 | 1.9338075857990327 |
| fccPt_SOC_fd56c832ef484b0b | fccPt_SOC_445821d99afc81fa | 1.2649205810264375 | 1.9399472516856235 | 1.9321522576465826 |

## Equal-FP32 pairs

| GPU row | CPU reference | S_moments | S_transport | S_whole |
|---|---|---:|---:|---:|
| fccPt_SOC_b32a6e11ea33ea0c | fccPt_SOC_c043d3069ccfe4aa | 24.853089819680793 | 11.028132330697053 | 10.129327725469135 |
| fccPt_SOC_e29f4db4737a6900 | fccPt_SOC_c043d3069ccfe4aa | 24.839293682143925 | 11.04078356670941 | 10.172773631944548 |
| fccPt_SOC_9b649a1d6b10bf2a | fccPt_SOC_c043d3069ccfe4aa | 24.729480951415997 | 11.000088835678124 | 10.128436789748655 |
| fccPt_SOC_b1fe04361371c71a | fccPt_SOC_c043d3069ccfe4aa | 24.72762045243157 | 11.035248646885334 | 10.139109404448039 |

## Mixed practical rows

4 mixed rows retained; they are never used for equal-precision headline speedups.

## Failed/ineligible and memory-limited rows


## Correctness summary

PASS: 20

## Definitions

- `S_moments = best_same_precision_CPU(P_moments_total) / GPU(P_moments_total)`.
- `S_transport = best_same_precision_CPU(T_transport_total) / GPU(T_transport_total)`.
- `S_whole = best_same_precision_CPU(whole_wall) / GPU(whole_wall)`.
- `T_transport_total` is the internal transport phase; `whole_wall` is the complete process invocation.
