# RF-05 reciprocal execution backends

RF-05 introduces `reciprocal_execution_backend` in `reciprocal_mod`.  The
public typed contract is `reciprocal_execution_capabilities`,
`reciprocal_execution_request`, and `reciprocal_execution_result`; the module
remains private by default.  A request identifies first/second-order assembly,
overlap, standard/generalized solution, eigenvector ownership, optional matrix
diagnostics, local point count, and operator generation without backend-name
strings in SCF or TD-DFT physics code.

`lapack_reciprocal_backend` is the sole factory selection (`lapack`/`cpu`) and
owns an RF-04 `reciprocal_assembler` plus a persistent `reciprocal_workspace`.
It reports host residency, standard and generalized Hermitian eigensystems,
values-only/vector solutions, both assembly orders, and overlap support.  It
uses the established `zheev`/`zhegv` conventions and performs Hermiticity and
positive-definite overlap checks immediately before each eigensolution.

## Dispatch and compatibility

Arbitrary-k and k+q eigenpairs make one deferred `execute_batch` call per
deduplicated RF-04 tile; the backend assembles and solves the tile without a
public H(k) handoff.  `reciprocal%build_kspace_hamiltonian` makes one backend
assembly call per RF-04 tile and retains the requested H(k) cache for bands,
DOS, and diagnostics.  `reciprocal%diagonalize_hamiltonian` is the companion
compatibility adapter: it presents that cache to the same backend and copies
the returned caller-owned eigenpair arrays into the established public fields.

Projected future call flow:

`SCF/TD-DFT orchestration -> factory -> CUDA execution_backend -> device Fourier tile -> device eigensolver -> result/consumer`

Thus a CUDA implementation requires one concrete backend, factory selection,
and build wiring; reciprocal k-point ownership, SCF/TD-DFT orchestration, and
physics formulas remain unchanged.

## CPU comparison

The solver calls and RF-04 tile allocation/query reuse behavior are unchanged
for arbitrary-k execution.  The legacy mesh-cache adapter currently copies its
cached H(k)/S(k) into its request to make request ownership explicit; this is a
small host-only transitional cost, not in the arbitrary-k resident path.  The
focused RF-04/RF-05 arbitrary-k test retains identical eigenvalue/eigenvector
results across tile sizes and verifies workspace reuse, generalized metric
coverage, refresh generation, values-only ownership, and reconstruction.
