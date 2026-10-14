# LDA+U remediation gate ledger

Controlling user review: `cleanup/RS_LMTO_ASA_LDAU_Remediation_SOC_Production_Protocol.md`, supplied outside this clean worktree (original preserved at `/Users/andersb/Jobb/rslmto_fable_v4b/cleanup/RS_LMTO_ASA_LDAU_Remediation_SOC_Production_Protocol.md`). The original scalar audit remains at `/Users/andersb/Jobb/rslmto_fable_v4b/cleanup/RS_LMTO_ASA_fable_v4b_LDAU_Codex_Audit.md`. Do not edit either user-owned file. This ledger records the fresh, unmodified remediation base and deliberately does not promote prior scalar test evidence to the new production gates.

## R00 baseline snapshot

**State: PARTIAL.** The baseline is frozen at commit `5660c070a0f73a6029c67997b1ff18823ae71369` on branch `fix/ldau-scf-soc-closure`; the worktree was clean before this file was added. The prior user branch remains at the same commit with its untracked `cleanup/` directory preserved. The remediation protocol's suggested branch name is used. No source or test edits have been made in this worktree at the time of these measurements.

| Check | Measured result | Command / provenance |
|---|---|---|
| Debug configure and full build | PASS, GNU Fortran 16.2.0; Accelerate BLAS/LAPACK; Debug; MPI/OpenMP/CUDA/libXC/spglib/fused reciprocal disabled | `cmake -S . -B /tmp/rslmto-ldau-remediation-baseline-debug -DCMAKE_BUILD_TYPE=Debug -DRUN_UNIT_TESTS=ON -DRUN_EXAMPLE_TESTS=ON -DRUN_REG_TESTS=ON -DENABLE_MPI=OFF -DENABLE_OPENMP=OFF -DENABLE_CUDA_PLUGIN=OFF -DENABLE_CUDA_RECIPROCAL=OFF -DENABLE_SPGLIB=OFF -DENABLE_FUSED_RECIPROCAL=OFF -DPython3_EXECUTABLE=/opt/homebrew/bin/python3.11 -DEXAMPLE_PYTHON_EXECUTABLE=/opt/homebrew/bin/python3.11`; `cmake --build /tmp/rslmto-ldau-remediation-baseline-debug -j4` |
| Registered tests | 202 total, 112 labeled `unit` | `ctest --test-dir /tmp/rslmto-ldau-remediation-baseline-debug -N` |
| Unit suite | 112/112 passed, 0 failures, 181.79 s | `ctest --test-dir /tmp/rslmto-ldau-remediation-baseline-debug -L unit --output-on-failure`; captured CTest log `/tmp/rslmto-ldau-remediation-baseline-debug/Testing/Temporary/LastTest.log`, SHA256 `62de2a311cd0341e7cf47fa1b75bf3087aa300ba6fb12ccc142d81d3b23dd2b7` |
| Baseline executable | SHA256 `1a7f1cbccf785c7fa4399033ff8ff15e5ed881109d0399f07c35822aa95c19bb` | `/tmp/rslmto-ldau-remediation-baseline-debug/bin/rslmto.x`; built from the clean base SHA above |
| Scalar no-U RS snapshot | PASS, nsp=1, one SCF step, exit 0; ETOT `-2541.9819038809896` Ry; SUMEV `-1.9987324216255` Ry; RMS `3.1e-5` | copied fixture under `/tmp/rslmto-remediate-r00-no-u/rs`; command from case directory: `RSLMTO_OMP_THREADS_SERIAL=1 /bin/bash /private/tmp/rslmto-ldau-remediate-worktree/tests/run_binary.sh /tmp/rslmto-remediate-r00-rslmto.x`; generated `input.nml` SHA256 `266453c476897240942edad552cfdd7b0f339fb94e6e1b12572297ab4e0568f5`; `Fe.nml` SHA256 `6ded298a3512b920211a8761b05a191a76402d3c985080380034cdff42596673`; `Fe_out.nml` SHA256 `dd21f69537b89e44fb3f039a2e9dae1c94ef886ae80794749280374bf3e71a13`; `testrun.log` SHA256 `822e99e2de0cc2c0ccc3c7f655b04a00451a44b37ad387d9293e07e3f1b6dc3c` |
| Scalar no-U KS snapshot | PASS, nsp=1, two SCF steps, exit 0; ETOT `-2541.9818791828561` Ry; SUMEV `-1.9985569107228` Ry; RMS `4.1e-5` | copied fixture under `/tmp/rslmto-remediate-r00-no-u/ks`; same command form as RS. `input.nml` SHA256 `0e84bc5d0d2b78eafb114f1a0ad3595f1ee517be7a1fceb68e40435728b952d8`; `Fe.nml` SHA256 `fcbf9541ca8c60940f9662246fc7acb841a44b7832b141241e32e8fa9caa7df1`; `Fe_out.nml` SHA256 `7de5fcc993a8d67ca4c7473bbe7e64cfa93255076541ce541f6a487fb7b02909`; `testrun.log` SHA256 `a448871bd68c09dda7dfc358506ece4a6d76342ea8902f85878a59ff4a979491` |

The two material snapshots are finite-step no-U baselines, not converged solutions and not an LDA+U energy identity. The source fixture inputs were copied; the user repository fixtures remain untouched. Earlier `PHASE0.md` evidence records pre-remediation environment-specific reciprocal no-U reference mismatches. Its executable was built at an earlier SHA/compiler configuration; reproduce those comparisons against this exact base binary before assigning them to this baseline. No reference values were changed.

**R00 work still pending:** fresh Release and MPI baseline builds/tests; exact reproduction/classification of the known reciprocal no-U reference discrepancies at this SHA; representative unchanged fixed-U input/oracle campaign snapshot. Debug build, full unit suite, and finite-step no-U RS/KS outputs are now recorded. Changes may begin only from the frozen commit above; pending baseline comparisons remain explicit and must not be hidden.

## R00–R15 release-gate status at baseline

| Row | State | Evidence or missing proof |
|---|---|---|
| R00 — baseline/environment | PARTIAL | Debug and unit baseline plus finite-step no-U RS/KS snapshots above. Release/MPI and exact no-U reference-mismatch reproduction remain pending. |
| R01 — tensor/basis/SOC conventions | PARTIAL | Prior Phase 1 documents scalar tensor, basis and rotation checks; SOC ladder signs and spinor convention are not certified. |
| R02 — scalar production projectors | BLOCKED | Separate RS/KS smoke outputs exist, but no common frozen production Hamiltonian comparison of full complex d occupations with an independent finite-H oracle. |
| R03 — scalar fixed-U functional | PARTIAL | Prior Phase 3 algebraic finite-difference and independent tensor oracle checks exist. They do not establish the production total energy. |
| R04 — corrected production ETOT | BLOCKED | No production frozen-state spectral-to-radial energy identity; do not add an assumed `EHub` or `EHub−Tr(VP)` term. |
| R05 — occupation/potential convergence | BLOCKED | No certified coupled P/V stopping criteria or converged/restarted production result. |
| R06 — scalar Fe production | NOT_RUN | No converged, refined nsp=1 U/J sweep in both RS and KS. |
| R07 — noncubic J production | NOT_RUN | No independently sourced, converged low-symmetry J≠0 material run. |
| R08 — spinor functional/FLL | NOT_RUN | No independent 10×10 spinor functional and double-counting oracle. |
| R09 — SOC insertion | NOT_RUN | SOC/Hubbard basis, ladder signs, and ξ→0 reduction are not validated. |
| R10 — spinor RS/KS occupations | NOT_RUN | No production four-spin-block projector/MPI comparison. |
| R11 — radial/XC fidelity | BLOCKED | No proof that the radial spin/XC update preserves the required transverse spinor information. |
| R12 — SOC+U Fe production | BLOCKED | Current fixed-U guard rejects nsp=2. Do not remove it or bypass it with reduced blocks before R08–R11 pass. |
| R13 — noncubic SOC production | NOT_RUN | No converged SOC-sensitive low-symmetry production case. |
| R14 — restart/MPI | PARTIAL | Earlier scalar complex restart and MPI checks exist; no spinor restart or full production MPI gate. |
| R15 — VAL-04/docs/regression | PARTIAL | Prior guarded scalar input campaign and unit tests exist; old nsp=2 VAL-04 is incompatible with current guard and must remain historical. A new converged scalar VAL-04 plus separate SOC campaign is not present. |

No SOC+U or material-physics support is certified by this ledger. Reevaluate every row after integrated production evidence; status changes must cite the source/binary SHA, independent assessor, command, exact input and raw output location, measured observable, criterion and negative control.
