# Sparse Matrix Infrastructure for RS-LMTO-ASA

## Overview

RS-LMTO-ASA now supports **Block Sparse Row (BSR)** format for Hamiltonian storage, enabling GPU acceleration via oneMKL and cuSPARSE.

## Files Added

### `source/sparse_matrix_mod.f90`
Core sparse matrix module with three formats:
- **BSR** (Block Sparse Row) — primary format for TB-LMTO
- **CSR** (Compressed Sparse Row) — element-wise sparse format
- **COO** (Coordinate) — triplet format for flexibility

### Updated Files
- `source/hamiltonian.f90` — Added `h_bsr: bsr_matrix` member and updated `block_to_sparse()` routine
- `source/CMakeLists.txt` — Added `sparse_matrix_mod.f90` to build system
- `example/sparse_matrix_demo.f90` — Demonstration/verification program

## BSR Format Specification

### Structure
```fortran
type :: bsr_matrix
   integer :: nblocks          ! Total number of blocks
   integer :: nrows            ! Number of block-rows (atom count)
   integer :: ncols            ! Number of block-columns
   integer :: blocksize        ! Block size (18 for spd basis)

   complex(rp), allocatable :: values(:,:,:)      ! (18, 18, nblocks)
   integer, allocatable     :: col_indices(:)     ! Column index per block
   integer, allocatable     :: row_ptr(:)         ! Row pointers (CSR-style)
end type
```

### Block Layout
- **blocksize = 18** (spd basis: s, px, py, pz, dxx, dxy, dxz, dyy, dyz, dzz per spin)
- **values(i, j, k)** — Element at row `i`, column `j` of block `k`
- **col_indices(k)** — Column atom index for block `k`
- **row_ptr(r)** — Starting index in `col_indices` for block-row `r`

### Example
For a BCC Fe cluster with 2 atoms and nearest-neighbor hops:
```
nblocks = 4  (atom 1→1 onsite, atom 1→2 hop, atom 2→1 hop, atom 2→2 onsite)
nrows = 2
blocksize = 18

values(:,:,1)     — H(atom 1, atom 1)
values(:,:,2)     — H(atom 1, atom 2)
values(:,:,3)     — H(atom 2, atom 1)
values(:,:,4)     — H(atom 2, atom 2)

col_indices = [1, 2, 1, 2]
row_ptr = [1, 3, 5]  (atom 1 blocks: indices 1-2, atom 2 blocks: indices 3-4)
```

## Building BSR Format

### Automatic (Recommended)
```fortran
call hamiltonian%block_to_sparse()
! Populates hamiltonian%h_bsr with proper BSR structure
```

### Manual
```fortran
use sparse_matrix_mod, only: allocate_bsr

call allocate_bsr(my_bsr, nblocks, nrows, blocksize)
! Allocate empty BSR
! Fill values(:,:,:) and col_indices(:) manually
! Set row_ptr(:) to define block-row structure
```

## Using with oneMKL (Intel Math Kernel Library)

### Sparse Matrix-Vector Product
```fortran
! oneMKL sparse_bsr_spmv for y = alpha*A*x + beta*y
! See: https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference/latest/sparse-blas-level-2-routines.html

! Signature (in oneAPI MKL):
! call mkl_sparse_z_bsr_spmv('N', alpha, A_handle, x, beta, y)
```

### Integration Example
```fortran
use onemkl_blas_interface

complex(rp) :: alpha, beta
complex(rp), allocatable :: x(:), y(:)

alpha = (1.0_rp, 0.0_rp)
beta = (0.0_rp, 0.0_rp)

allocate(x(hamiltonian%h_bsr%nrows * nb))
allocate(y(hamiltonian%h_bsr%nrows * nb))

! Call oneMKL sparse kernel
! (Requires Fortran interface from oneMKL Fortran bindings)
```

## Using with cuSPARSE (NVIDIA CUDA)

### Sparse Matrix-Vector Product
```fortran
! cuSPARSE uses block sparse format similarly to oneMKL
! cusparseZbsrmv for y = alpha*A*x + beta*y
! See: https://docs.nvidia.com/cuda/cusparse/using.html

! Example (with cufortran):
! call cusparseZbsrmv(handle, dir, m, n, nnz, alpha, &
!                     descrA, bsrVal, bsrRowPtr, bsrColInd, &
!                     blockDim, x, beta, y)
```

## Memory Requirements

### BSR vs Dense Comparison
For a 100-atom cluster with spd basis (nb=18):

**Dense Format:**
- Size: (100×18)² = 32,400² ≈ 1 GB (complex*8)

**BSR Format (with ~5 blocks per atom-row on average):**
- Size: 500 blocks × 18 × 18 × 8 bytes ≈ 64 MB (5× compression)

## Conversion Routines

### BSR to CSR
```fortran
use sparse_matrix_mod, only: bsr_to_csr

type(csr_matrix) :: csr_h
call bsr_to_csr(hamiltonian%h_bsr, csr_h)
! Expands to element-wise sparsity (useful for visualization, debugging)
```

### Element-wise Sparsity Pattern
- BSR compression: N blocks × (18×18) = 5×–10× typical compression
- CSR expansion: Only ~5% non-zero elements on average (wide empty regions)

## Validation

### Check BSR Structure
```fortran
call hamiltonian%h_bsr%print_info()
! Output:
! BSR: 500 blocks, 100 block-rows, blocksize=18
! Total elements: 1620000
```

### Verify Sparsity
```fortran
integer :: sparsity_pct
real(rp) :: compression_ratio

! Number of stored elements
integer :: stored = hamiltonian%h_bsr%nblocks * 18 * 18

! Full matrix size
integer :: full = (hamiltonian%h_bsr%nrows * 18)**2

sparsity_pct = 100 * (1 - real(stored) / real(full))
compression_ratio = real(full) / real(stored)

print *, 'Sparsity:', sparsity_pct, '%'
print *, 'Compression ratio:', compression_ratio, ':1'
```

## GPU Backend Roadmap

### Phase 1: Sparse Matrix Infrastructure ✅
- [x] BSR format definition and allocation
- [x] `block_to_sparse()` builds BSR automatically
- [x] CSR conversion for compatibility
- [x] Tests pass with new structure

### Phase 2: oneMKL Integration
- [ ] oneMKL sparse BLAS initialization
- [ ] SpMV with `mkl_sparse_z_bsr_spmv`
- [ ] Chebyshev scaling via oneMKL kernels
- [ ] Performance benchmarking

### Phase 3: cuSPARSE Integration
- [ ] CUDA kernel wrappers for cusparseZbsrmv
- [ ] Host↔device memory management
- [ ] K-space band structure GPU acceleration

### Phase 4: Recursion Engine GPU Backend
- [ ] GPU `apply()` kernel using sparse operations
- [ ] Overlap with CPU computation (pipelined SpMV)
- [ ] Lanczos/Chebyshev on GPU

## Testing

### Run Sparse Matrix Demo
```bash
# Compile with main code
make

# Run example that uses block_to_sparse()
./bin/rslmto.x < example/bulk/bccFe_scf_single/input.nml

# Check BSR structure in output
grep "BSR:" report.out
```

### Performance Profiling
```bash
# With oneMKL enabled:
export MKL_VERBOSE=1
./bin/rslmto.x < input.nml  2>&1 | grep sparse

# With cuSPARSE enabled (future):
nsys profile --trace cuda,nvtx ./bin/rslmto.x < input.nml
```

## Known Limitations

1. **Atom Type Indexing:** `block_to_sparse()` assumes `lattice%iz(atom)` maps correctly in bulk region. Verify for multi-type systems.

2. **Operator Contents:** The default BSR build fills from `hall`/`ee`. GPU CCOR no-`hoh` runs build BSR from the merged `hall+hallcc` and `ee+eecc` operator blocks. Orthogonalization correction (`hallo`/`eeo`) is still stored separately and is not folded into the BSR matrix.

3. **Non-Collinear Magnetism:** BSR layout is element-wise (includes spin structure implicitly). Works with existing array layout.

4. **Extended Hopping Range:** BSR respects neighbor list from `lattice%nn`. Beyond 2-hop shells are not included by default.

## References

- [NVIDIA cuSPARSE Documentation](https://docs.nvidia.com/cuda/cusparse/)
- [Intel oneMKL Sparse BLAS](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference/latest/)
- [Blöchl et al., Phys. Rev. B 50, 17953 (1994)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.50.17953) — LMTO conventions

## Future Work

- [ ] Add `h_bsr_hoh` for orthogonalization correction
- [ ] GPU-accelerated Lanczos/Chebyshev with sparse operations
- [ ] Heterogeneous CPU-GPU pipelining for large clusters
- [ ] Sparse matrix equilibration and preconditioners
- [ ] Benchmarking against full-potential codes

---
*Document version: 1.0 (June 10, 2026)*
