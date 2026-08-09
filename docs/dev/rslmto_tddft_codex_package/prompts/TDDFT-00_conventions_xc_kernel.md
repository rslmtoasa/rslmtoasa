# Codex Prompt — TDDFT-00: Conventions and XC-Response Kernel Plumbing

## Mission

Establish the exact response conventions and the ground-state XC-response interface required by LR-TDDFT.

**Do not implement `chi_KS` or the TDDFT Dyson solver in this task.**

Read first:

- `README.md`
- `00_MASTER_BLUEPRINT.md`
- `01_PHYSICS_CONVENTIONS.md`
- current ground-state potential/XC and Hamiltonian construction code.

The public `fable_v2` B11 document is stale in places. Inspect current source before making assumptions.

## Physics contract

The long-term response space is

\[
(n,m_x,m_y,m_z),
\]

but the first production physics will be transverse `chi^{+-}`.

The code's exact definitions of magnetic density, moment, Pauli matrices, XC magnetic field, and Hamiltonian spin splitting must be established from the implementation itself.

Literature expressions such as

\[
K_\perp \propto B_{\rm xc}/m
\]

must not be copied with unexplained factors.

## Tasks

- [x] Trace the current electronic Hamiltonian convention `H=H0*I + Hvec.sigma` through the actual source.
- [x] Trace the ground-state XC functional path that generates spin-dependent potentials.
- [x] Document whether internal magnetic quantities are spin density, magnetic moment, or another normalization.
- [x] Pin the exact energy units and any `mu_B` factors.
- [x] Pin the code's Fourier/time convention needed for a retarded `omega+i eta` response.
- [x] Introduce response-component constants/enumeration for `CHARGE, MX, MY, MZ, PLUS, MINUS`.
- [x] Add `response_basis` helpers that construct `sigma_0,sigma_x,sigma_y,sigma_z,sigma^+,sigma^-` in the exact current spin-orbital ordering.
- [x] Add an `xc_response_kernel` abstraction/provider.
- [x] Expose enough ground-state XC information for a future site-projected transverse ALSDA kernel.
- [x] Design the provider so it can later expose:
  - [x] `dVxc/dn`
  - [x] `dVxc/dm`
  - [x] `dBxc/dn`
  - [x] `dBxc/dm`
  - [x] `K_perp`
- [x] Add tests for Pauli/circular operator algebra in the exact basis ordering.
- [x] Add a two-level spin-split Hamiltonian fixture verifying response-operator normalization.
- [x] Add developer documentation warning against identifying a derived Hamiltonian block with `B_xc` without proof.
- [x] Add developer documentation warning against copying factors of 2 or `mu_B` from a paper without translating conventions.

## Important prohibition

Do **not** simply set something like:

```fortran
Kperp = hamiltonian%hxc / moment
```

unless source tracing proves that `hxc` is exactly the required local XC field with the right units and normalization. The current Hamiltonian magnetic blocks contain representation-dependent structure and are not automatically the TDDFT kernel input.

## Non-goals

- no `chi_KS`;
- no Dyson enhancement;
- no Goldstone correction;
- no GBT changes;
- no SOC response;
- no longitudinal response.

## Acceptance

The convention must be recoverable from executable tests and documentation. Existing regression tests must remain unchanged.

## Completion protocol

When done:

1. tick every completed checkbox above;
2. explicitly list any unresolved convention instead of guessing;
3. run focused tests and relevant regressions;
4. make one commit with the single-line message:

`TDDFT-00: establish response conventions and XC kernel interface`
