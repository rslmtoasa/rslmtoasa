# ACC-06 reciprocal CPU/GPU crossover report

Campaign date: 2026-08-18. The report contains the generated full campaign
summary; raw JSON samples remain machine-local evidence under `/tmp`.

## Environment and method

- GNU Fortran 13.3.0, Release, Intel oneMKL, `OMP_NUM_THREADS=2`, one MPI rank.
- Two NVIDIA RTX A4000 devices visible; driver 610.57.04; CUDA toolkit 13.3.
- One warm-up and three measured repetitions per executable invocation.
- Production reciprocal H(R)->H(k) assembly with standard-Hermitian CPU LAPACK
  and CUDA backend requests.
- Matrix dimensions 8, 18, 36, and 72; Nk values 1, 8, and 32; tile sizes 1,
  8, and 16; eigenvalues-only and eigenvectors-requested modes.
- 216 executable runs, summarized into 72 CPU-vs-GPU rows.

`CPU s` and `GPU s` are ACC-00 harness wall-time medians, including executable
setup and backend initialization. GPU timing therefore includes H2D, cuSOLVER,
D2H, and synchronization. The `parallel` CPU rows are the benchmark-only
OpenMP per-k loop; `backend` is the production typed LAPACK path.

## Full generated crossover table

| fixture | matrix | Nk | tile | eigenvectors | best CPU | CPU s | GPU s | speedup | recommended |
|---|---:|---:|---:|:---:|---|---:|---:|---:|---|
| Si_sp | 8 | 1 | 1 | no | backend | 0.006445 | 0.404295 | 0.02x | lapack |
| Si_sp | 8 | 1 | 1 | yes | parallel | 0.010855 | 0.401114 | 0.03x | lapack |
| Si_sp | 8 | 1 | 16 | no | backend | 0.006323 | 0.400603 | 0.02x | lapack |
| Si_sp | 8 | 1 | 16 | yes | parallel | 0.007515 | 0.399688 | 0.02x | lapack |
| Si_sp | 8 | 1 | 8 | no | backend | 0.007405 | 0.400338 | 0.02x | lapack |
| Si_sp | 8 | 1 | 8 | yes | parallel | 0.011321 | 0.400232 | 0.03x | lapack |
| Si_sp | 8 | 32 | 1 | no | parallel | 0.004103 | 0.407989 | 0.01x | lapack |
| Si_sp | 8 | 32 | 1 | yes | backend | 0.009932 | 0.407124 | 0.02x | lapack |
| Si_sp | 8 | 32 | 16 | no | backend | 0.008360 | 0.403905 | 0.02x | lapack |
| Si_sp | 8 | 32 | 16 | yes | parallel | 0.005798 | 0.404149 | 0.01x | lapack |
| Si_sp | 8 | 32 | 8 | no | backend | 0.004496 | 0.411878 | 0.01x | lapack |
| Si_sp | 8 | 32 | 8 | yes | parallel | 0.005587 | 0.407638 | 0.01x | lapack |
| Si_sp | 8 | 8 | 1 | no | parallel | 0.004120 | 0.406531 | 0.01x | lapack |
| Si_sp | 8 | 8 | 1 | yes | parallel | 0.010790 | 0.401577 | 0.03x | lapack |
| Si_sp | 8 | 8 | 16 | no | backend | 0.004954 | 0.400769 | 0.01x | lapack |
| Si_sp | 8 | 8 | 16 | yes | backend | 0.008618 | 0.406491 | 0.02x | lapack |
| Si_sp | 8 | 8 | 8 | no | parallel | 0.003702 | 0.400886 | 0.01x | lapack |
| Si_sp | 8 | 8 | 8 | yes | parallel | 0.009279 | 0.399418 | 0.02x | lapack |
| bccFe_spd | 18 | 1 | 1 | no | backend | 0.005149 | 0.404797 | 0.01x | lapack |
| bccFe_spd | 18 | 1 | 1 | yes | parallel | 0.005905 | 0.403504 | 0.01x | lapack |
| bccFe_spd | 18 | 1 | 16 | no | backend | 0.007666 | 0.401874 | 0.02x | lapack |
| bccFe_spd | 18 | 1 | 16 | yes | backend | 0.010455 | 0.401648 | 0.03x | lapack |
| bccFe_spd | 18 | 1 | 8 | no | backend | 0.011815 | 0.401000 | 0.03x | lapack |
| bccFe_spd | 18 | 1 | 8 | yes | backend | 0.010131 | 0.403257 | 0.03x | lapack |
| bccFe_spd | 18 | 32 | 1 | no | parallel | 0.009333 | 0.419399 | 0.02x | lapack |
| bccFe_spd | 18 | 32 | 1 | yes | backend | 0.009570 | 0.419647 | 0.02x | lapack |
| bccFe_spd | 18 | 32 | 16 | no | parallel | 0.007363 | 0.417419 | 0.02x | lapack |
| bccFe_spd | 18 | 32 | 16 | yes | backend | 0.013974 | 0.418659 | 0.03x | lapack |
| bccFe_spd | 18 | 32 | 8 | no | parallel | 0.010825 | 0.419422 | 0.03x | lapack |
| bccFe_spd | 18 | 32 | 8 | yes | parallel | 0.010205 | 0.419720 | 0.02x | lapack |
| bccFe_spd | 18 | 8 | 1 | no | parallel | 0.011517 | 0.407276 | 0.03x | lapack |
| bccFe_spd | 18 | 8 | 1 | yes | backend | 0.010881 | 0.410648 | 0.03x | lapack |
| bccFe_spd | 18 | 8 | 16 | no | backend | 0.010622 | 0.405787 | 0.03x | lapack |
| bccFe_spd | 18 | 8 | 16 | yes | parallel | 0.012191 | 0.409920 | 0.03x | lapack |
| bccFe_spd | 18 | 8 | 8 | no | backend | 0.008063 | 0.404456 | 0.02x | lapack |
| bccFe_spd | 18 | 8 | 8 | yes | parallel | 0.009149 | 0.409600 | 0.02x | lapack |
| multisite_4_lmax_2 | 72 | 1 | 1 | no | parallel | 0.006047 | 0.411652 | 0.01x | lapack |
| multisite_4_lmax_2 | 72 | 1 | 1 | yes | parallel | 0.014158 | 0.410421 | 0.03x | lapack |
| multisite_4_lmax_2 | 72 | 1 | 16 | no | parallel | 0.007405 | 0.413166 | 0.02x | lapack |
| multisite_4_lmax_2 | 72 | 1 | 16 | yes | backend | 0.015186 | 0.412365 | 0.04x | lapack |
| multisite_4_lmax_2 | 72 | 1 | 8 | no | backend | 0.007676 | 0.408925 | 0.02x | lapack |
| multisite_4_lmax_2 | 72 | 1 | 8 | yes | backend | 0.009467 | 0.417131 | 0.02x | lapack |
| multisite_4_lmax_2 | 72 | 32 | 1 | no | parallel | 0.017499 | 0.488784 | 0.04x | lapack |
| multisite_4_lmax_2 | 72 | 32 | 1 | yes | parallel | 0.024916 | 0.518635 | 0.05x | lapack |
| multisite_4_lmax_2 | 72 | 32 | 16 | no | parallel | 0.015720 | 0.484366 | 0.03x | lapack |
| multisite_4_lmax_2 | 72 | 32 | 16 | yes | parallel | 0.017442 | 0.519825 | 0.03x | lapack |
| multisite_4_lmax_2 | 72 | 32 | 8 | no | parallel | 0.018101 | 0.489585 | 0.04x | lapack |
| multisite_4_lmax_2 | 72 | 32 | 8 | yes | parallel | 0.016741 | 0.520800 | 0.03x | lapack |
| multisite_4_lmax_2 | 72 | 8 | 1 | no | backend | 0.012528 | 0.440175 | 0.03x | lapack |
| multisite_4_lmax_2 | 72 | 8 | 1 | yes | parallel | 0.011743 | 0.444445 | 0.03x | lapack |
| multisite_4_lmax_2 | 72 | 8 | 16 | no | backend | 0.015042 | 0.438713 | 0.03x | lapack |
| multisite_4_lmax_2 | 72 | 8 | 16 | yes | backend | 0.013536 | 0.442784 | 0.03x | lapack |
| multisite_4_lmax_2 | 72 | 8 | 8 | no | backend | 0.010635 | 0.442991 | 0.02x | lapack |
| multisite_4_lmax_2 | 72 | 8 | 8 | yes | parallel | 0.013675 | 0.444491 | 0.03x | lapack |
| two_site_spd | 36 | 1 | 1 | no | parallel | 0.005518 | 0.406729 | 0.01x | lapack |
| two_site_spd | 36 | 1 | 1 | yes | parallel | 0.011871 | 0.408215 | 0.03x | lapack |
| two_site_spd | 36 | 1 | 16 | no | parallel | 0.009330 | 0.407863 | 0.02x | lapack |
| two_site_spd | 36 | 1 | 16 | yes | backend | 0.014560 | 0.410681 | 0.04x | lapack |
| two_site_spd | 36 | 1 | 8 | no | backend | 0.006302 | 0.399424 | 0.02x | lapack |
| two_site_spd | 36 | 1 | 8 | yes | parallel | 0.009014 | 0.405910 | 0.02x | lapack |
| two_site_spd | 36 | 32 | 1 | no | backend | 0.012901 | 0.467135 | 0.03x | lapack |
| two_site_spd | 36 | 32 | 1 | yes | backend | 0.010545 | 0.480411 | 0.02x | lapack |
| two_site_spd | 36 | 32 | 16 | no | parallel | 0.011485 | 0.460973 | 0.02x | lapack |
| two_site_spd | 36 | 32 | 16 | yes | backend | 0.012239 | 0.480873 | 0.03x | lapack |
| two_site_spd | 36 | 32 | 8 | no | backend | 0.017034 | 0.464225 | 0.04x | lapack |
| two_site_spd | 36 | 32 | 8 | yes | parallel | 0.014580 | 0.479579 | 0.03x | lapack |
| two_site_spd | 36 | 8 | 1 | no | parallel | 0.016033 | 0.418292 | 0.04x | lapack |
| two_site_spd | 36 | 8 | 1 | yes | parallel | 0.007357 | 0.432090 | 0.02x | lapack |
| two_site_spd | 36 | 8 | 16 | no | parallel | 0.007418 | 0.419329 | 0.02x | lapack |
| two_site_spd | 36 | 8 | 16 | yes | parallel | 0.008813 | 0.422898 | 0.02x | lapack |
| two_site_spd | 36 | 8 | 8 | no | parallel | 0.010639 | 0.422726 | 0.03x | lapack |
| two_site_spd | 36 | 8 | 8 | yes | parallel | 0.013668 | 0.424761 | 0.03x | lapack |

## Decision

No CPU/GPU crossover was observed in the measured 8--72 matrix, Nk=1--32
range. The largest measured GPU/CPU ratio was 0.05x. Retain LAPACK as the
explicit default, retain the existing caller tile default of 16, and do not
add automatic CPU/GPU dispatch or a universal CUDA preferred tile size.

This is a bounded result for these fixtures and does not claim that larger,
production-weighted matrices cannot cross over.
