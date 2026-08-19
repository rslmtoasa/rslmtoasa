# ACC-13 — RS KPM/Kubo-Bastin CUDA coverage

## Decision

ACC-13 is complete for the existing RS CUDA moment implementation.  The audit
found a complete shared transport path, so no new GPU kernel was added:

```text
real-space Hamiltonian and current operators (CPU)
    -> rsrec_stochastic_moments (CUDA)
    -> existing mu_nm_stochastic array
    -> CPU Gamma/Kubo-Bastin reconstruction
```

Charge, spin, and orbital conductivity all select the same stochastic moment
contract after building their production operators.  This is the right
reuse boundary: a separate CUDA kernel per observable would duplicate the
same Chebyshev recurrence and change no transport arithmetic.  The distinct
`chebyshev_orbital_mod` production route also has an existing CUDA
`rsrec_orbital_moments` entry point.

## Audit result

The existing source/API surface contains the following relevant paths:

| Production stage | CUDA entry point | Status |
|---|---|---|
| Stochastic charge moments | `rsrec_stochastic_moments` | supported and validated |
| Stochastic spin moments | `rsrec_stochastic_moments` | supported and validated |
| Stochastic orbital moments | `rsrec_stochastic_moments` | supported and validated |
| Separate orbital-moment workflow | `rsrec_orbital_moments` | supported by the existing CSR/BSR path |
| Gamma/Kubo-Bastin reconstruction | existing Fortran `conductivity` path | retained on CPU |

The CUDA moment path is HOH-aware when the required velocity
orthogonalization operators are present.  The existing guards continue to
reject structured FFT/conv orbital moments, `ccor_2c+hoh`, and local-axis
routes explicitly.  These are unsupported combinations, not silent CPU
fallbacks presented as GPU results.

## Correctness evidence

The low-level hardware validator was run with the CUDA plugin in FP64
validation mode:

```text
tests/cuda/build_and_validate.sh /tmp/rslmto-acc13-lowlevel
```

All 15 ABI routes passed.  The maximum displayed relative error was about
`2.2e-15`, including stochastic moments and the separate orbital-moment
contract.  The existing CPU `UnitKpmTransport` and the static CUDA surface
contract also passed.

The production-level probe is:

```text
python3 tests/validation/acc13_kpm_cuda.py \
  --binary build-acc01-cuda/bin/rslmto.x \
  --scratch-root /tmp/rslmto-acc13-transport
```

It runs the real SOC fcc-Pt conductivity fixture once in CPU mode and once
with `control%gpu_plugin=.true.` for each of charge, spin, and orbital.  It
compares the complex output nearest the Fermi energy and requires the CUDA
log marker `CUDA_DEVICE_MAPPING` before accepting the GPU result.

Representative measured rows on the available RTX A4000 were:

| Probe | Observable | CPU/GPU wall ratio | complex relative error |
|---|---|---:|---:|
| `cond_ll=20`, replication 4 | charge | 0.61x | 4.9e-8 |
| `cond_ll=20`, replication 4 | spin | 0.66x | 7.3e-7 |
| `cond_ll=20`, replication 4 | orbital | 0.62x | 1.1e-8 |
| `cond_ll=40`, replication 6 | charge | 1.26x | 6.3e-7 |
| `cond_ll=40`, replication 6 | spin | 1.33x | 4.3e-7 |
| `cond_ll=40`, replication 6 | orbital | 1.34x | 1.0e-6 |

The ratios are end-to-end process wall measurements including setup and CPU
postprocessing.  They are intentionally reported as workload evidence, not
as a portable threshold.  Small transport workloads remain CPU-preferred;
the larger measured workload is CUDA-beneficial on this GPU.

## Scope boundaries

- The default production CUDA Chebyshev path remains FP32 for moments with
  FP64 host output.  The low-level FP64 run is a numerical reference check;
  ACC-13 does not promote a new production precision policy.
- CPU operator construction and CPU Gamma/Kubo-Bastin reconstruction remain
  canonical and were not rewritten.
- No automatic CPU/GPU dispatch was added.  Users must request the CUDA
  plugin explicitly, consistent with the frozen Phase III policy.
- The probe is a focused fcc-Pt SOC campaign.  It does not claim all
  materials, GPU architectures, system sizes, or stochastic-vector counts.

## Completion checklist

- [x] transport CPU/GPU dataflow mapped
- [x] existing GPU moment support validated
- [x] charge conductivity CPU/GPU compared
- [x] spin conductivity CPU/GPU compared
- [x] orbital conductivity CPU/GPU compared
- [x] production operator conventions retained
- [x] performance measured at small and larger KPM workloads
- [x] no redundant kernel added
- [x] end-to-end speedup recorded
- [x] ACC-14 support-matrix inputs documented
