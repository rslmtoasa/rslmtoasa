# C5.2 Part A — RS no-U identity attribution

| Check | Evidence | Verdict |
|---|---|---|
| A1 determinism | Canonical manifest-resolved `Example_bulk_bccFe_nsp2_block`; clean `8b8f7ef` Release binary; two runs, one thread; all scientific/result artifacts byte-identical. `str.out` differs only in `cpu_s` (0.295 vs 0.246); `testrun.log` is runtime metadata. | **PASS** |
| A2 parent | `17732a7` built with matched Release compiler/link settings; one run against the same fixture snapshot. Scientific artifacts byte-identical to A1 baseline; `str.out` only `cpu_s` (0.127), plus runtime log. | **PASS** |
| A2 R1 | `2686e4e` built with matched Release compiler/link settings; one run against the same fixture snapshot. Scientific artifacts byte-identical to A1 and parent; `str.out` only `cpu_s` (0.167), plus runtime log. | **PASS** |
| A3 prior pair | Existing matched Release OMP=2 outputs: nine files had the small printed-value changes detailed below and in the JSON sidecar. This is descriptive prior evidence, not attribution under the new one-thread controls. | RECORD |
| A4 | No scientific difference is attributable to R1 or R2 under the prescribed one-thread matched-build controls. | **STOP (iv)** |

The canonical fixture resolved `strux_lib` and legacy Chebyshev backend, `recur='block'`, `nsp=2`, `lld=20`, `nstep=1`; live logs show RS spectral mode and `use_kspace=false`. This explicitly specified standard fixture seed was used from the same immutable input snapshot in all four runs. Each solver exited 0; all four PIDs were confirmed gone. Total Part A usage: **4 of 8** solver runs.

Only `2686e4e` changes production source in `908074f..HEAD`, modifying `source/lattice_cluster.f90`; its parent `17732a7` is source-identical to `908074f`. `96f1d72` is documentation and the four Co reference inputs; the later C5 commits contain no production-source changes. Thus neither the R1 change nor an R2 production change explains the earlier nine-file OMP=2 result in these A2 runs.

| Prior OMP=2 file | Maximum absolute difference | Maximum relative difference* | Baseline-zero entries |
|---|---:|---:|---:|
| `Fe_1_spinene.out` | 4.216×10⁻¹⁵ | 16.0411 | 0 |
| `Fe_dos.out` | 0 | 0 | 3,260 |
| `Fe_orbene.out` | 4.5315×10⁻¹⁵ | 13.9552 | 5 |
| `Fe_orbital_dos.out` | 0 | 0 | 20,182 |
| `Fe_orbquadene.out` | 1×10⁻¹¹ | 23.2172 | 6 |
| `Fe_out.nml` | 3.794011234043×10⁻¹⁵ | 3.5725×10⁻⁴ | 0 |
| `Fe_scf.nml` | 3.794011234043×10⁻¹⁵ | 3.5725×10⁻⁴ | 0 |
| `magneticdos.out` | 0 | 0 | 9,363 |
| `radial_ground_state_Fe_1.dat` | 3.62039682628931×10⁻¹⁵ (radial header moment vector); radial data exact | 0.00373242 (header components) | 0 |

*Relative differences are `|candidate−baseline|/|baseline|` only for nonzero baselines; no denominator floor was used. Large ratios occur in near-zero orbital fields. `ql(1,:,spin)` s/p/d occupations, energy fields, DOS values, and report-precision physical spin moments are exact in the prior pair; physical sphere spin is read from the radial header `reported_moment_xyz_muB`. Namelist `mom` is the normalized spin reference axis; `lmom` is the orbital-moment vector. The exact per-quantity values and zero-baseline counts are in [the A3 JSON](C5_2_PART_A_assets/audit/a3_omp2_metrics.json).

**Interpretation and close (five lines):**

- Passed: A1 determinism and A2 parent/R1 scientific identity under the registered one-thread protocol.
- Recorded: prior OMP=2 nine-file metrics separately; they are not proven to be thread nondeterminism.
- Provisional: none; no Part A state is provisional.
- Stopped: A4(iv); no commit attribution, tolerance change, fix restriction, or criterion relaxation is adopted.
- Next: user review only may consider a future amendment authorizing B–D if A1 and every matched A2 comparison are scientifically byte-identical under the registered one-thread protocol; this proposal authorizes nothing now.
