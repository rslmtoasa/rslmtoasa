# RS-LMTO-ASA SCF-B0C-RS — Real-space benchmark instrumentation report

Status: SCF-B0C-RS complete; SCF-B1 remains next.

This report records the real-space extension of the shared SCF-B0C harness. It
does not create a second benchmark stack: reciprocal and real-space rows use
the same production probe, output schema, pairing rules, profile closure gate,
raw-log archive, and JSON/CSV/Markdown writers.

## Scope and production routes

The benchmark executable is `ReciprocalAccP1bPhysicalSCF`. The route is chosen
explicitly with `--scf-route`:

- `reciprocal` retains the existing H(k) / eigensolver path.
- `real_space --rs-solvers block` exercises block Haydock recursion and block
  Green reconstruction.
- `real_space --rs-solvers chebyshev` exercises the production Chebyshev
  recursion and spectral reconstruction dispatch.
- `real_space --rs-solvers lanczos` exercises scalar/Lanczos recursion and the
  scalar Green path.

CPU production smoke runs completed for the block/CSR bcc-Fe, fast-Chebyshev
diamond-Si, and scalar/Lanczos bcc-Fe routes; each emitted `PASS` with profile
closure and the expected route-specific phase populated.

The RS physics key includes the solver, backend, recursion depth, block size,
terminator, Chebyshev order/kernel, spectral-bounds policy, material fixture,
starting state, mixing, smearing, and precision fields. A CPU/GPU pair with a
mismatched RS control is rejected with an explicit reason rather than being
reported as a speedup.

## Instrumentation boundaries

The shared `SCF_B0C_ITER` record now contains the following exclusive RS
phases:

| Field | Production boundary |
|---|---|
| `P_rs_hamiltonian_prepare` | real-space potential and Hamiltonian construction |
| `P_rs_solver_kernel` | selected recursion kernel, including the backend upload/setup boundary |
| `P_rs_green_function` | scalar Green or block `zsqr`/Green reconstruction |
| `P_rs_spectral_reconstruct` | Chebyshev spectral reconstruction dispatch |
| `P_rs_energy_integration` | reserved for a standalone production integration boundary |
| `P_rs_fermi` | Fermi-level calculation and Hubbard-U SC Fermi work |
| `P_rs_density_build` | production density/moment construction |
| `P_rs_charge_spin_accumulate` | magnetic charge/spin accumulation |

`P_rs_energy_integration` is zero for the current bulk fixtures because the
production `bands` routines combine the integration with Fermi, magnetic, and
moment construction; a standalone timer would double-count or invent a
boundary. The field is emitted so the schema remains stable when a distinct
integration boundary is added.

The profile implementation avoids nested-stage overwrites. In particular, the
CUDA Hamiltonian upload is included in `P_rs_solver_kernel`, where it is part
of the selected backend invocation. This keeps the exclusive phase sum closed
and prevents upload time from leaking into `P_scf_misc`.

The RS detail fields are emitted as required schema fields:

`T_rs_H2D`, `T_rs_kernel`, `T_rs_D2H`, `T_rs_sync`, and `T_rs_setup`.

`T_rs_kernel` is the measured outer RS solver boundary and equals
`P_rs_solver_kernel`. The current plugin does not expose separate transfer,
synchronization, or setup counters at the Fortran benchmark boundary, so those
fields remain zero and `rs_detail_timers_status=not_exposed_by_backend` is
emitted. No transfer time is inferred from a wall-clock remainder.

## Correctness and fallback policy

An explicit CUDA RS request is unsupported if the plugin is not compiled or if
the production call returns without an associated GPU backend. The probe emits
`SCF_B0C status=UNSUPPORTED` with a reason for either case. A successful RS
CUDA row must report:

- `rs_gpu_used=true`;
- `fallback_detected=false`;
- `rs_kernel_correctness_status=PASS`; and
- `rs_kernel_invariant=finite_coefficients`.

The kernel status is a structured finite-coefficient/no-fallback invariant. A
converged GPU SCF row is additionally compared with the converged FP64 CPU
oracle using total energy, Fermi level, charge, site moment, charge transfer,
and residual tolerances. Mixed RS CUDA precision is labelled `mixed`: the
production plugin uses FP32 working paths while Hamiltonian/density/canonical
SCF state remains FP64. It is visible and oracle-checked, but is not an
equal-precision headline.

## Closure evidence

The target was built with GNU Fortran 13.3, oneMKL, OpenMP, and the CUDA RS
plugin. On the RTX A4000 / CUDA 13.3 host, the fixed two-iteration bcc-Fe block
smoke campaign reported:

| backend | route | profile closure | max steady misc fraction | RS GPU | fallback | RS kernel |
|---|---|---:|---:|---|---|---|
| CPU LAPACK | real-space block/CSR | `3.09e-9` | `3.02e-5` | false | false | PASS |
| CUDA block | real-space block/CSR | `7.30e-10` | `7.16e-5` | true | false | PASS |

Both satisfy the 3% profile-closure and 5% steady-misc gates. A 40-iteration
`scf_convergence` campaign also kept both profiles within those gates, reaching
final residuals of `8.64e-9` (CPU) and `1.67e-8` (CUDA), while the executable’s
convergence flag remained false for this fixture. The harness therefore
correctly retains the CUDA row as a non-headline oracle comparison rather than
calling it a converged correctness pass. The campaign did establish the
production GPU path and showed matching charge and site-moment observables to
the displayed precision.

Representative fixed-smoke medians were:

| backend | `P_rs_solver_kernel` (s) | `P_rs_green_function` (s) | `P_scf_iteration_total` (s) |
|---|---:|---:|---:|
| CPU LAPACK | 1.1425 | 1.1125 | 2.5924 |
| CUDA block | 0.8180 | 0.0322 | 0.9985 |

These are closure measurements for the fixture, not SCF-B1 performance
claims. `S_rs_kernel` and `S_rs_phase` are emitted by the harness only when a
valid paired candidate satisfies the correctness and precision policy.

## Reproduction

CPU-only instrumentation smoke:

```bash
python3 tests/benchmarks/scf_b0c.py \
  --binary build-b1-frozen-cuda/bin/ReciprocalAccP1bPhysicalSCF \
  --build-dir build-b1-frozen-cuda \
  --output results/benchmarks/scf_b0c_rs/cpu.json \
  --materials fe --scf-route real_space --rs-solvers block \
  --rs-backend csr --omp-threads 1 --nstep 2 \
  --skip-cuda --benchmark-level scf_iteration
```

CUDA closure/convergence campaign:

```bash
python3 tests/benchmarks/scf_b0c.py \
  --binary build-b1-frozen-cuda/bin/ReciprocalAccP1bPhysicalSCF \
  --build-dir build-b1-frozen-cuda \
  --output results/benchmarks/scf_b0c_rs/campaign.json \
  --materials fe --scf-route real_space --rs-solvers block \
  --rs-backend csr --omp-threads 1 --nstep 40 \
  --benchmark-level scf_convergence
```

The campaign writes the canonical JSON, CSV, Markdown, raw logs, correctness
records, iteration history, and skipped/unsupported rows beside the selected
output. No machine-local timing artifact is required in the repository.

## Decision

SCF-B0C-RS is complete: real-space instrumentation is on the shared harness,
production route identity is explicit, CUDA fallback is not silently accepted,
profile closure is enforced, and the current backend’s detail-timer boundary
is documented honestly. SCF-B1 is next for frozen repetitions, material
matrix, and final performance policy.
