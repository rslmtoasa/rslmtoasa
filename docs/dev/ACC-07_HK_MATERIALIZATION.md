# ACC-07 H(k) materialization evidence

Campaign date: 2026-08-19. Decision: **Decision A — no production-code change**.

The validated CUDA normal-mesh route assembles H(k) on the host, passes the
request-owned tile to the CUDA eigensolver, and copies that same host tile into
the public `reciprocal%hk_bulk` compatibility cache. It does not copy H(k)
from device to host. Therefore there is no H(k) transfer to eliminate, and the
cache cannot be removed unconditionally without changing downstream behavior.

## Consumer map

| Consumer | H(k) requirement | Evidence |
|---|---|---|
| Dyson reciprocal Green backend | Required for every `(k,z)` inverse | `source/reciprocal_green.f90:436-504` |
| Bloch spectral function | Required for every `(k,E)` inverse on the path | `source/reciprocal_bsf.f90:90-130` |
| Bands compatibility/legacy path | Required by the non-fused adapter; fused path keeps the cache alongside eigenpairs | `source/reciprocal_bands.f90:23-112` |
| Gamma/multisite diagnostics | Reads the first cached H(k) | `source/reciprocal_bands.f90:186-188`, `source/reciprocal_fourier.f90:1239-1307` |
| SCF occupations and Gaussian DOS | Consume eigenvalues/eigenvectors after fused execution; they do not independently require H(k) | `source/reciprocal_occupations.f90`, `source/reciprocal_dos.f90` |

The cache is allocated for the local k ownership in
`source/reciprocal_fourier.f90:894-902`, filled from the request-owned host
tile at `source/reciprocal_fourier.f90:1144-1148`, and released by
`source/reciprocal_lifecycle.f90`. The Green and BSF consumers make the full
cache a real compatibility product, not dead materialization.

## Transfer measurement

Fresh CUDA converged Si/sp and bcc-Fe/spd 4x4x4 Gaussian-mesh runs were made
with the current ACC-05 build (`build-acc01-cuda`, CUDA 13.3, RTX A4000,
driver 610.57.04). The final 48-tile SCF records were:

| case | normal-mesh total | host assembly | H2D | GPU solve | D2H eigenpairs | D2H H(k) |
|---|---:|---:|---:|---:|---:|---:|
| Si/sp | 0.038 s | 0.000 s | 0.000834 s | 0.473242 s | 0.001178 s | 0.0 s (N/A) |
| bcc-Fe/spd | 0.042 s | 0.002 s | 0.001094 s | 0.522105 s | 0.001434 s | 0.0 s (N/A) |

The `d2h_hamiltonian_seconds=0.0(not_applicable)` field is emitted by the
production normal-mesh timing record. H(k) materialization is a host-side
assignment from the already assembled input tile, not a device transfer. Its
measured transfer fraction is therefore exactly 0% of the CUDA phase.

## Host-memory measurement

`complex(real64)` payload is 16 bytes, so the exact full-cache payload is:

```text
Hk_bytes = 16 * nmat * nmat * nk_local
```

The measured allocation shapes and representative payloads are:

| workload | nmat | nk/local | full `hk_bulk` payload |
|---|---:|---:|---:|
| Si/sp 4x4x4 | 8 | 64 | 0.0625 MiB |
| bcc-Fe/spd 4x4x4 | 18 | 64 | 0.3164 MiB |
| ACC-06 four-site fixture, Nk=32 | 72 | 32 | 2.5313 MiB |

The host tile is bounded by the configured tile size and is not an additional
full-mesh allocation. For scale awareness, an unvalidated 16x16x16 four-site
projection would be 324 MiB; that is a future memory gate, not evidence from a
validated production workflow. The current validated scope is far below the
ACC-07 action threshold of either 5% of reciprocal phase time in H(k) transfer
or 256 MiB of required full-cache payload.

## Decision and regression evidence

No request flag or lazy-cache semantics were added. The existing behavior is
retained because Green/BSF callers require H(k), while the CUDA path has no
H(k) transfer to optimize. The representative ACC-06 repeat after this decision
reported CUDA wall times of 0.405--0.450 s versus 0.064--0.074 s for the best
CPU rows, with no crossover in the four measured fixtures. This repeats the
no-code conclusion without making a new performance claim.

Current-build LAPACK and CUDA converged Si/Fe runs both passed with matching
electron count, EF, band energy, total energy, DOS integral, and Fe moment.
The existing reciprocal CUDA normal-mesh/arbitrary-k tests and SCF/DOS
validation remain unchanged.

Reopen ACC-07 if a validated workload crosses either action threshold or if a
new downstream GPU consumer makes repeated eigensystem/H(k) movement
material; then introduce an explicit `eigensystem-only` versus
`H(k) compatibility-cache` request distinction and preserve legacy callers.
