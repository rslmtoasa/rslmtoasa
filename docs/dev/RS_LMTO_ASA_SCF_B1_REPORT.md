# SCF-B1R — Reconciled unified SCF campaign

## 1. Executive result

The archived campaign is complete: **19/19 planned cases**. Canonical row statuses are {'PASS': 20, 'FAIL': 1, 'NOT_CONVERGED': 1}.

The build, Fe3 CPU OMP sweep, reciprocal Fe2/Fe3 runtime metadata, and current-build Fe3/Fe4/Fe5 frozen-potential rows are available. The common-potential RS-vs-k-space accuracy gate is **not passed**, so that formulation comparison is reported as **INCONCLUSIVE** and no cross-route timing claim is made.

SCF-B1 can close only as a scoped performance campaign; it does not close as a universal RS-vs-k-space accuracy/performance claim. The unresolved item is the failed fixed-potential accuracy contract, not missing runtime evidence.

## B1R2 lean desktop update

The desktop-friendly B1R2 lean tier completed **18/18 cases** on 2026-08-26.
It adds production-route Chebyshev Si1/Si2 measurements, Fe2/Fe3 reciprocal
common-state checks, and canonical CPU/GPU ratios. The detailed evidence and
precision qualifications are in [RS_LMTO_ASA_SCF_B1R2_REPORT.md](RS_LMTO_ASA_SCF_B1R2_REPORT.md).

This is a lean tier, not the 36-case full HPC matrix: it uses one measured
process per case, no warmup, `nstep=1`, Si1 CPU OMP1/4/8, and Si2 CPU OMP8.
Therefore the lean desktop scope is closed, while the full HPC scaling scope
remains intentionally pending.

## 2. Build-integrity check

| build | compiler | type | effective optimization | status | current source commit | dirty | reuse B1 timing? |
|---|---|---|---|---|---|---|---|
| /home/andersb/Nano/rslmto_fable_v3/build-b1r-release-cuda | /usr/bin/gfortran | Release | representative actual compile commands have optimized effective flags; no trailing -O0 | PASS | 5edf2f622b096b258467b6611dc03d5b80bf4972 | yes | forbidden |

The authoritative preflight inspected representative compile commands for the SCF driver, RS recursion, Hamiltonian construction, and reciprocal path. Conclusion: `representative actual compile commands have optimized effective flags; no trailing -O0`. The B1 timings were not reused because `B1 Release timings were produced before the B1R -O0 integrity fix and must be rerun.`

## 3. Method and status policy

Performance is the steady production SCF-iteration wall time and its route-specific phases. SCF cycle count is retained for correctness/stability only. `PASS`, `NOT_CONVERGED`, `NOT_APPLICABLE`, `UNSUPPORTED`, `FAIL`, and `SKIPPED` remain distinct; a successful process exit does not override a failed profile or accuracy gate.

Equal-precision `S_*` ratios require an explicit eligible pairing. Mixed RS CUDA rows use `R_*` production ratios only. A missing or ineligible ratio is printed as `-` with the row's reason.

## 4. Reciprocal Fe2/Fe3 basis and timer resolution

The runtime evidence resolves the labels as physical supercell replications with genuinely enlarged reciprocal matrices: Fe2 is `nmat=144`, Fe3 is `nmat=486`, and all K1 rows report `Nk_unique=512`. These values are taken from the emitted runtime result, not inferred from fixture names.

| case | physical replication | reciprocal basis | nmat | Nk_unique | backend | strategy | mode | timer | timer s | timer check | row status | reason |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| k1_fe1_cpu | 1x1x1 | runtime:nmat=18 | 18 | 512 | lapack | fp64_zheevd | fp64 | P_eigensolver | 2.3055e-05 | INCONCLUSIVE | FAIL | scf_misc_exceeds_5_percent |
| k1_fe1_cuda | 1x1x1 | runtime:nmat=18 | 18 | 512 | cuda | fp64_zheevd | fp64 | T_solver | 0.318086 | PASS | PASS | profile_failed, correctness_failed, correctness_missing |
| k1_fe2_cpu | 2x2x2 | runtime:nmat=144 | 144 | 512 | lapack | fp64_zheevd | fp64 | P_eigensolver | 1.8507e-05 | PASS | PASS | common-state or route comparison is required |
| k1_fe2_cuda | 2x2x2 | runtime:nmat=144 | 144 | 512 | cuda | fp64_zheevd | fp64 | T_solver | 3.18449 | PASS | PASS | correctness_missing |
| k1_fe3_cpu | 3x3x3 | runtime:nmat=486 | 486 | 512 | lapack | fp64_zheevd | fp64 | P_eigensolver | 1.4703e-05 | PASS | PASS | common-state or route comparison is required |
| k1_fe3_cuda | 3x3x3 | runtime:nmat=486 | 486 | 512 | cuda | fp64_zheevd | fp64 | T_solver | 14.8105 | PASS | PASS | correctness_missing |

CPU `P_eigensolver` and CUDA `T_solver` are the authoritative boundaries. The timer validation is PASS for all successful Fe2/Fe3 rows; the Fe1 CPU row is INCONCLUSIVE because its row is FAIL under the strict profile-misc gate, although the positive timer value is archived.

## 5. Fe3 real-space CPU OMP reconciliation

| case | solver/backend | nrec | nmat | mode | OMP | kernel s | phase s | iteration s | status | equal precision | reason |
|---|---|---|---|---|---|---|---|---|---|---|---|
| rs1_fe3_cpu_omp1 | block/csr | 20 | 486 | fp64 | 1 | 198.142 | 220.8 | 222.154 | PASS | yes | common-state or route comparison is required |
| rs1_fe3_cuda_mixed | block/csr | 20 | 486 | mixed | 1 | 73.4488 | 76.5084 | 78.0204 | PASS | no | numeric_mode_mismatch, mixed_precision_not_headline, correctness_missing |
| rs1_fe3_cpu_omp2 | block/csr | 20 | 486 | fp64 | 2 | 120.957 | 134.418 | 140.407 | PASS | yes | common-state or route comparison is required |
| rs1_fe3_cpu_omp4 | block/csr | 20 | 486 | fp64 | 4 | 89.0645 | 97.9332 | 100.142 | PASS | yes | common-state or route comparison is required |
| rs1_fe3_cpu_omp8 | block/csr | 20 | 486 | fp64 | 8 | 69.7083 | 76.4622 | 78.3134 | PASS | yes | common-state or route comparison is required |

Best measured Fe3 CPU kernel row: `rs1_fe3_cpu_omp8` at `69.7083 s`; best steady-iteration row: `rs1_fe3_cpu_omp8` at `78.3134 s`. The sweep is OMP1/2/4/8 and the CUDA row remains labelled mixed precision.

The corrected mixed-GPU production ratio is shown only where the common-state/profile gate supports it. The RS1 Fe3 CUDA row is retained as a measured mixed route but does not receive an invented headline ratio when its comparison lane is not eligible.

## 6. RS common-potential correctness anchor

| case | backend | mode | starting state | status | correctness | reason | energy | Fermi | moment | charge | residual |
|---|---|---|---|---|---|---|---|---|---|---|---|
| rs2_reference_fe_cpu | lapack | fp64 | - | NOT_CONVERGED | NOT_CONVERGED | reference run did not reach the requested convergence criterion | -2541.98 | -0.0692808 | 2.10401 | 8 | 1.11616e-08 |
| rs2_common_fe_cpu_fp64 | lapack | fp64 | rs2_reference_fe_cpu | PASS | PASS | reference observable comparison passed | -2541.98 | -0.0692808 | 2.10401 | 8 | 1.32772e-08 |
| rs2_common_fe_cuda_mixed | cuda | mixed | rs2_reference_fe_cpu | PASS | PASS | reference observable comparison passed | -2541.98 | -0.0692807 | 2.10401 | 8 | 7.42631e-09 |
| x1_fe_cpu_rs | lapack | fp64 | rs2_reference_fe_cpu | PASS | PASS | reference observable comparison passed | -2541.98 | -0.0692808 | 2.10401 | 8 | 1.06421e-08 |

The reference run itself is `NOT_CONVERGED`, but its archived restart and final observables support the common-state RS CPU/CUDA comparison. The common-state CPU and CUDA rows are retained with their explicit correctness outcomes; iteration counts are not used as an accelerator metric.

## 7. Reciprocal SCF rows and eligibility

| case | size | nmat | Nk_unique | backend | solver strategy | mode | iteration s | full wall s | S_iteration | S_solver | equal precision | headline eligible | status | reason |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| k1_fe1_cpu | 1x1x1 | 18 | 512 | lapack | fp64_zheevd | fp64 | 0.134574 | 1.25376 | - | - | no | - | FAIL | scf_misc_exceeds_5_percent |
| k1_fe1_cuda | 1x1x1 | 18 | 512 | cuda | fp64_zheevd | fp64 | 0.431381 | 2.56052 | - | - | yes | no | PASS | profile_failed, correctness_failed, correctness_missing |
| k1_fe2_cpu | 2x2x2 | 144 | 512 | lapack | fp64_zheevd | fp64 | 3.43256 | 16.1128 | - | - | yes | - | PASS | common-state or route comparison is required |
| k1_fe2_cuda | 2x2x2 | 144 | 512 | cuda | fp64_zheevd | fp64 | 5.63489 | 22.6357 | - | - | yes | no | PASS | correctness_missing |
| k1_fe3_cpu | 3x3x3 | 486 | 512 | lapack | fp64_zheevd | fp64 | 45.7414 | 159.079 | - | - | yes | - | PASS | common-state or route comparison is required |
| k1_fe3_cuda | 3x3x3 | 486 | 512 | cuda | fp64_zheevd | fp64 | 27.1157 | 98.6837 | - | - | yes | no | PASS | correctness_missing |
| x1_fe_cpu_kspace | 1x1x1 | 18 | 512 | lapack | lapack | fp64 | 0.117684 | 1.45943 | - | - | yes | - | PASS | fermi_energy, site_moment, final_residual |

The Fe1 CPU row is a real `FAIL` under `scf_misc_exceeds_5_percent`; it is not silently treated as unsupported. Reciprocal FP64 CUDA rows are equal-precision at the numeric-mode gate, but their `S_*` ratios remain suppressed unless the full correctness/headline pairing gate passes.

## 8. Current-build large reciprocal FP64 rerun

These K2 rows are frozen-potential eigensolver evidence with eigenvectors enabled, not full SCF convergence rows. Their current-build CPU/GPU speedups are computed only from paired PASS rows at the same L and nmat.

| L | nmat | Nk unique | backend | strategy | vectors | steady total s | cold wall s | status | reason | source commit |
|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 486 | 1 | cuda | fp64_zheevd | yes | 0.0311128 | 3.09616 | PASS | - | 5edf2f622b096b258467b6611dc03d5b80bf4972 |
| 3 | 486 | 1 | lapack | lapack | yes | 0.0463756 | 2.72462 | PASS | - | 5edf2f622b096b258467b6611dc03d5b80bf4972 |
| 4 | 1152 | 1 | cuda | fp64_zheevd | yes | 0.15919 | 7.68535 | PASS | - | 5edf2f622b096b258467b6611dc03d5b80bf4972 |
| 4 | 1152 | 1 | lapack | lapack | yes | 0.515287 | 9.28238 | PASS | - | 5edf2f622b096b258467b6611dc03d5b80bf4972 |
| 5 | 2250 | 1 | cuda | fp64_zheevd | yes | 0.85469 | 16.551 | PASS | - | 5edf2f622b096b258467b6611dc03d5b80bf4972 |
| 5 | 2250 | 1 | lapack | lapack | yes | 3.65694 | 34.949 | PASS | - | 5edf2f622b096b258467b6611dc03d5b80bf4972 |

| L | nmat | CPU/GPU speedup | eligibility | reason |
|---|---|---|---|---|
| 3 | 486 | 1.49057 | PASS | paired current-build FP64 frozen-potential rows |
| 4 | 1152 | 3.23692 | PASS | paired current-build FP64 frozen-potential rows |
| 5 | 2250 | 4.27868 | PASS | paired current-build FP64 frozen-potential rows |

The current measured crossover is already GPU-beneficial at Fe3 (`nmat=486`) and becomes more favorable at Fe4 and Fe5. This is a frozen-potential solver result, not a universal full-SCF conclusion.

## 9. Fixed-potential RS vs k-space Fe comparison

| accuracy status | energy difference | Fermi difference | moment difference | charge difference | timing eligible | reason |
|---|---|---|---|---|---|---|
| FAIL | 0.000843741 | 0.000833804 | 0.0810614 | 3.55449e-11 | no | accuracy tolerance exceeded: total_energy, fermi, moment |

The fixed-potential comparison does not satisfy the declared energy, Fermi, and moment tolerances. Accordingly the RS-vs-k-space scientific conclusion is **INCONCLUSIVE**; no RS/k-space timing speedup or DOS overlay is presented. Charge agrees closely, but that does not repair the failed observables.

## 10. Canonical artifacts and plots

The campaign directory contains `build_preflight.json`, `manifest.json`, `campaign.json`, `campaign.csv`, `campaign.md`, raw per-case stdout/stderr/commands, correctness records, iteration histories, and `cross_route/x1.json`. Regenerated plots:

- `results/benchmarks/scf_b1r/plots/01_rs_kernel_vs_size.png`
- `results/benchmarks/scf_b1r/plots/02_rs_iteration_vs_size.png`
- `results/benchmarks/scf_b1r/plots/03_rs_stage_fractions.png`
- `results/benchmarks/scf_b1r/plots/04_rs_cpu_omp_scaling.png`
- `results/benchmarks/scf_b1r/plots/05_rs_production_ratios_vs_size.png`
- `results/benchmarks/scf_b1r/plots/06_large_reciprocal_solver_speedup_vs_nmat.png`
- `results/benchmarks/scf_b1r/plots/07_reciprocal_stage_fractions.png`
- `results/benchmarks/scf_b1r/plots/08_scf_residual_history.png`

## 11. Final workload conclusions

- The Release build is validated from actual compile commands, with effective `-O3`; old B1 timing reuse is forbidden.
- Fe2/Fe3 reciprocal rows are genuinely enlarged runtime matrices (`144`/`486`) at `Nk_unique=512`; the CPU/GPU timer names are recorded separately and validated per route.
- The Fe3 RS CPU OMP sweep is complete, and the mixed RS CUDA label is preserved. Production ratios are shown only when their lane is eligible.
- Current-build frozen-potential reciprocal solver evidence shows GPU benefit from Fe3 through Fe5, with the crossover visible at `nmat=486` in this measured workload.
- The common-potential RS-vs-k-space accuracy gate fails for energy, Fermi level, and moment. That comparison remains inconclusive and must not be used to choose a universal formulation.

## 12. Closeout decision

B1R evidence collection and report reconciliation are complete. SCF-B1 may close for the scoped, separately qualified workload rows, but it should not be described as a closed universal RS-vs-k-space comparison until the fixed-potential accuracy contract is repaired and rerun.

## Provenance

```json
{
  "build_preflight": {
    "build_dir": "/home/andersb/Nano/rslmto_fable_v3/build-b1r-release-cuda",
    "build_type": "Release",
    "cache_optimization_flags_non_authoritative": "-O3",
    "compiler": "/usr/bin/gfortran",
    "effective_optimization_by_component": {
      "hamiltonian": "-O3",
      "reciprocal": "-O3",
      "rs_recursion": "-O3",
      "scf_driver": "-O3"
    },
    "effective_optimization_conclusion": "representative actual compile commands have optimized effective flags; no trailing -O0",
    "generator": "Ninja",
    "reason": null,
    "representative_compile_commands": {
      "hamiltonian": [
        "/usr/bin/gfortran -DUSE_CUDA_PLUGIN -DUSE_CUDA_RECIPROCAL -DUSE_SPGLIB -I/usr/include -I/usr/local/cuda-13.3/targets/x86_64-linux/include  -ffree-line-length-0 -cpp -march=native -ffree-line-length-0 -cpp -march=native -DCOLOR -fopenmp -DOpenMP_Fortran_FOUND -O3 -Jmodules -c /home/andersb/Nano/rslmto_fable_v3/source/hamiltonian.f90 -o source/CMakeFiles/rslmto.dir/hamiltonian.f90.o",
        "/usr/bin/gfortran -DUSE_CUDA_PLUGIN -DUSE_CUDA_RECIPROCAL -DUSE_SPGLIB -I/usr/include -I/usr/local/cuda-13.3/targets/x86_64-linux/include  -ffree-line-length-0 -cpp -march=native -ffree-line-length-0 -cpp -march=native -DCOLOR -fopenmp -DOpenMP_Fortran_FOUND -O3 -Jmodules -c /home/andersb/Nano/rslmto_fable_v3/source/hamiltonian_build.f90 -o source/CMakeFiles/rslmto.dir/hamiltonian_build.f90.o"
      ],
      "reciprocal": [
        "/usr/bin/gfortran -DUSE_CUDA_PLUGIN -DUSE_CUDA_RECIPROCAL -DUSE_SPGLIB -I/usr/include -I/usr/local/cuda-13.3/targets/x86_64-linux/include  -ffree-line-length-0 -cpp -march=native -ffree-line-length-0 -cpp -march=native -DCOLOR -fopenmp -DOpenMP_Fortran_FOUND -O3 -Jmodules -c /home/andersb/Nano/rslmto_fable_v3/source/reciprocal.f90 -o source/CMakeFiles/rslmto.dir/reciprocal.f90.o",
        "/usr/bin/gfortran -DUSE_CUDA_PLUGIN -DUSE_CUDA_RECIPROCAL -DUSE_SPGLIB -I/usr/include -I/usr/local/cuda-13.3/targets/x86_64-linux/include  -ffree-line-length-0 -cpp -march=native -ffree-line-length-0 -cpp -march=native -DCOLOR -fopenmp -DOpenMP_Fortran_FOUND -O3 -Jmodules -c /home/andersb/Nano/rslmto_fable_v3/source/reciprocal_bands.f90 -o source/CMakeFiles/rslmto.dir/reciprocal_bands.f90.o",
        "/usr/bin/gfortran -DUSE_CUDA_PLUGIN -DUSE_CUDA_RECIPROCAL -DUSE_SPGLIB -I/usr/include -I/usr/local/cuda-13.3/targets/x86_64-linux/include  -ffree-line-length-0 -cpp -march=native -ffree-line-length-0 -cpp -march=native -DCOLOR -fopenmp -DOpenMP_Fortran_FOUND -O3 -Jmodules -c /home/andersb/Nano/rslmto_fable_v3/source/self_reciprocal.f90 -o source/CMakeFiles/rslmto.dir/self_reciprocal.f90.o"
      ],
      "rs_recursion": [
        "/usr/bin/gfortran -DUSE_CUDA_PLUGIN -DUSE_CUDA_RECIPROCAL -DUSE_SPGLIB -I/usr/include -I/usr/local/cuda-13.3/targets/x86_64-linux/include  -ffree-line-length-0 -cpp -march=native -ffree-line-length-0 -cpp -march=native -DCOLOR -fopenmp -DOpenMP_Fortran_FOUND -O3 -Jmodules -c /home/andersb/Nano/rslmto_fable_v3/source/green_block.f90 -o source/CMakeFiles/rslmto.dir/green_block.f90.o",
        "/usr/bin/gfortran -DUSE_CUDA_PLUGIN -DUSE_CUDA_RECIPROCAL -DUSE_SPGLIB -I/usr/include -I/usr/local/cuda-13.3/targets/x86_64-linux/include  -ffree-line-length-0 -cpp -march=native -ffree-line-length-0 -cpp -march=native -DCOLOR -fopenmp -DOpenMP_Fortran_FOUND -O3 -Jmodules -c /home/andersb/Nano/rslmto_fable_v3/source/recursion.f90 -o source/CMakeFiles/rslmto.dir/recursion.f90.o",
        "/usr/bin/gfortran -DUSE_CUDA_PLUGIN -DUSE_CUDA_RECIPROCAL -DUSE_SPGLIB -I/usr/include -I/usr/local/cuda-13.3/targets/x86_64-linux/include  -ffree-line-length-0 -cpp -march=native -ffree-line-length-0 -cpp -march=native -DCOLOR -fopenmp -DOpenMP_Fortran_FOUND -O3 -Jmodules -c /home/andersb/Nano/rslmto_fable_v3/source/recursion_core.f90 -o source/CMakeFiles/rslmto.dir/recursion_core.f90.o"
      ],
      "scf_driver": [
        "/usr/bin/gfortran -DUSE_SPGLIB -DVERSION=\\\"5edf-dirty\\\" -I/usr/include  -ffree-line-length-0 -cpp -march=native -ffree-line-length-0 -cpp -march=native -DCOLOR -fopenmp -DOpenMP_Fortran_FOUND -O3 -Jmodules -fopenmp -c /home/andersb/Nano/rslmto_fable_v3/source/main.f90 -o CMakeFiles/rslmto.x.dir/source/main.f90.o",
        "/usr/bin/gfortran -DUSE_CUDA_PLUGIN -DUSE_CUDA_RECIPROCAL -DUSE_SPGLIB -I/usr/include -I/usr/local/cuda-13.3/targets/x86_64-linux/include  -ffree-line-length-0 -cpp -march=native -ffree-line-length-0 -cpp -march=native -DCOLOR -fopenmp -DOpenMP_Fortran_FOUND -O3 -Jmodules -c /home/andersb/Nano/rslmto_fable_v3/source/calculation.f90 -o source/CMakeFiles/rslmto.dir/calculation.f90.o",
        "/usr/bin/gfortran -DUSE_CUDA_PLUGIN -DUSE_CUDA_RECIPROCAL -DUSE_SPGLIB -I/usr/include -I/usr/local/cuda-13.3/targets/x86_64-linux/include  -ffree-line-length-0 -cpp -march=native -ffree-line-length-0 -cpp -march=native -DCOLOR -fopenmp -DOpenMP_Fortran_FOUND -O3 -Jmodules -c /home/andersb/Nano/rslmto_fable_v3/source/self.f90 -o source/CMakeFiles/rslmto.dir/self.f90.o"
      ]
    },
    "schema": "rslmto.scf-b1r.build-preflight.v1",
    "source_cmake_flag_fix": "cmake/SetFortranFlags.cmake no longer appends -O0 to Release flags",
    "source_commit": "5edf2f622b096b258467b6611dc03d5b80bf4972",
    "status": "PASS",
    "timestamp": "2026-08-24T04:18:41.726761+00:00",
    "timing_reuse_from_B1": "forbidden",
    "timing_reuse_reason": "B1 Release timings were produced before the B1R -O0 integrity fix and must be rerun."
  },
  "policy": {
    "cycle_count_role": "correctness and stability only",
    "ineligible_speedups_are_suppressed": true,
    "performance_metric": "steady production SCF iteration wall time"
  },
  "provenance": {
    "commit": "5edf2f622b096b258467b6611dc03d5b80bf4972",
    "dirty": true,
    "host": "alcazar"
  }
}
```
