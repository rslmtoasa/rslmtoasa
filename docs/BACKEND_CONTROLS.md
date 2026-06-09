# Backend Controls

This document describes the recursion backend controls used to select CPU or
plugin-backed operator and recurrence implementations.

## Control Namelist Keys

The controls live in the `control` namelist and are defined in:
- `source/control.f90`
- `source/include_codes/namelists/control.f90`

Available keys:

- `recur_backend`
  - `cpu_reference`
  - `cpu_sparse_reference`
  - `external_plugin`
- `recur_backend_plugin`
  - Shared library path for `external_plugin`
- `recur_backend_library`
  - Backend selector string passed to the plugin
- `recur_precision`
  - `complex_fp64`
  - `complex_fp32`
- `export_hamiltonian`
  - Enables block-sparse operator export
- `export_hamiltonian_path`
  - Prefix for exported operator files
- `validate_backend_roundtrip`
  - Compares exported-operator application against the legacy Hamiltonian path

## Current Backends

### `cpu_reference`

Uses the legacy Fortran Hamiltonian and recurrence implementation. This is the
correctness baseline.

### `cpu_sparse_reference`

Uses the exported block-sparse operator representation on CPU. This exercises
the new backend seam without loading a shared library.

### `external_plugin`

Loads a shared library that implements the plugin ABI. The Fortran side keeps
ownership of physics orchestration and output arrays, while the plugin owns the
backend implementation.

## Plugin ABI Families

The backend ABI currently distinguishes these execution families:

- operator upload and apply
- scalar Lanczos
- block Lanczos
- block Chebyshev
- stochastic transport

Plugins report their support through `rslmto_backend_capabilities()`. The
Fortran dispatcher checks those capabilities and logs explicit fallback when a
requested family is not implemented by the plugin.

## Current Plugin Coverage

### CPU sample plugins

- `rslmto_cpu_sparse_plugin`
- `rslmto_openmp_sparse_plugin`
- `rslmto_onemkl_sparse_plugin`

Current status:

- operator apply: implemented
- scalar Lanczos: fallback
- block Lanczos: fallback
- block Chebyshev: fallback
- transport: fallback

These plugins currently advertise `recurrences=apply_only`.

### CUDA plugin

Plugin target:

- `rslmto_cuda_sparse_plugin`

Current intended rollout:

- operator upload/apply: CUDA path
- block Chebyshev: CUDA backend entrypoint
- block Lanczos: CUDA backend entrypoint
- transport: not yet completed

The first CUDA backend remains:

- block-first
- complex-first
- `nb`-driven
- non-collinear/SOC-safe
- `hoh`-aware
- complex FP64 only

## Build Controls

Top-level CMake options:

- `BUILD_BACKEND_PLUGINS`
- `BUILD_ONEMKL_PLUGIN`
- `BUILD_CUDA_PLUGIN`

Example configure commands:

```bash
cmake -S . -B build
cmake --build build -j
```

With oneMKL:

```bash
cmake -S . -B build -DBUILD_ONEMKL_PLUGIN=ON -DONEMKL_ROOT=/path/to/oneMKL
cmake --build build -j
```

With CUDA:

```bash
cmake -S . -B build -DBUILD_CUDA_PLUGIN=ON
cmake --build build -j
```

## Example Namelist Snippets

CPU reference:

```fortran
recur_backend = 'cpu_reference'
recur_precision = 'complex_fp64'
export_hamiltonian = .false.
validate_backend_roundtrip = .false.
```

Plugin backend:

```fortran
recur_backend = 'external_plugin'
recur_backend_plugin = '/path/to/librslmto_cuda_sparse_plugin.so'
recur_backend_library = 'cuda'
recur_precision = 'complex_fp64'
export_hamiltonian = .true.
export_hamiltonian_path = 'backend_export'
validate_backend_roundtrip = .true.
```

## Logging Behavior

The backend dispatcher logs:

- plugin load and capability strings
- operator upload/apply activation
- recurrence-family plugin execution
- explicit fallback when the plugin does not implement a requested family

This is intentional. Production acceptance should not rely on silent CPU
fallback for missing GPU kernels.
