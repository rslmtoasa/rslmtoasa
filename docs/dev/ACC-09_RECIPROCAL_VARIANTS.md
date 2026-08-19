# ACC-09 — Validated reciprocal CUDA operator variants

ACC-09 extends the reciprocal CUDA claim only across CPU-established
standard-Hermitian variants.  The CUDA backend remains a dense eigensolver
adapter: the existing CPU reciprocal assembler constructs each H(k), including
the second-order/HOH correction and optional CCOR term, and the resulting tile
is handed to the same CUDA `zheevd` path used for ordinary H(k).

## Scope audit

| Variant | CPU evidence | Mathematical class | CUDA status |
|---|---|---|---|
| First-order H(k) | reciprocal k-space modes and ordinary DOS/band fixtures | Standard Hermitian, identity overlap | Supported through host assembly |
| Second-order / HOH | `Example_density_of_states_bccFe_reciprocal_hoh`, VAL-02 extended-spin fixture, and WP6a HOH report | Standard Hermitian, identity overlap | Supported through host assembly |
| First-order + CCOR | `Example_density_of_states_bccFe_ccor_2c`, CCOR production/regression evidence | Standard Hermitian, identity overlap | Supported through host assembly |
| HOH + CCOR | No promoted Phase-II reciprocal fixture | Standard-Hermitian algebra is not enough to promote the combination | Not claimed |
| `generalized_overlap_proxy` | CPU/LAPACK path only; generalized metric is separately checked | Generalized Hermitian | Rejected by CUDA capability contract |
| `generalized_overlap_kanpur` | Formal route is not implemented | Development/unsupported | Not claimed; no CUDA support |

The ordinary first/second-order and CCOR rows do not require a new GPU
operator representation.  `reciprocal_fourier.f90` remains the single
operator-construction path.  `reciprocal_cuda.cpp` contains only CUDA context,
cuSOLVER, transfer, timing, and lifecycle code; it has no HOH or CCOR input
logic.

## Correctness evidence

The opt-in campaign is:

```bash
python3 tests/validation/val09_reciprocal_variants.py \
  --cpu-binary build-acc09-cpu/bin/rslmto.x \
  --cuda-binary build-acc09-cuda/bin/rslmto.x \
  --scratch-root /tmp/acc09-reciprocal-variants \
  --output /tmp/acc09-reciprocal-variants.json
```

It runs the compact bcc-Fe reciprocal DOS fixture for HOH and CCOR, compares
the complete printed `dos_kspace.dat` downstream output, and checks canonical
EF, electron count, and band energy.  `ACC05_TIMING` records are retained in
the JSON report for performance interpretation; they are not correctness
thresholds.

The existing committed CPU post-processing references also pass directly with
the CUDA backend for both `Example_density_of_states_bccFe_reciprocal_hoh` and
`Example_density_of_states_bccFe_ccor_2c`.

The fresh-build compact campaign produced the following timing/correctness
evidence on the validation host (GNU Fortran 13.3, one OpenMP thread, CUDA
13.3, NVIDIA GPU):

| Variant | CPU total | CUDA total | max DOS error | max canonical error |
|---|---:|---:|---:|---:|
| HOH | 0.001 s | 0.020 s | 0.0 | 0.0 |
| CCOR | 0.001 s | 0.018 s | 0.0 | 0.0 |

The compact matrices are below the measured CPU/CUDA crossover, so ACC-09
records correctness coverage only and makes no reciprocal CUDA speedup claim.

## Unsupported scope

CUDA does not advertise generalized Hermitian support, overlap assembly, or
device-side first/second-order assembly.  The normal-mesh and arbitrary-k
adapters reject generalized CUDA requests at the typed capability boundary.
No generalized-overlap, HOH+CCOR, SOC/noncollinear promotion beyond the
existing CPU evidence, GBT, TD-DFT, or multi-GPU claim is made here.
