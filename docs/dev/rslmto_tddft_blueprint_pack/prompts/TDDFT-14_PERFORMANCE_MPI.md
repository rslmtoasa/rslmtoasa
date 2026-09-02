# TDDFT-14 — performance, MPI decomposition, and backend crossover

## Agent mission

Optimize the now-validated response backends using measured production profiles, with special attention to the different reuse patterns of k-GF and r-GF algorithms.

## Context

Halle reports near-embarrassing parallelism over response points. The RS route additionally allows `chi0(R,w)` to be reused for many q values, so naïvely distributing q can throw away the best reuse opportunity. Recent codewide GPU experience also warns that small kernels can be overhead dominated.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.
- Do not change physics formulas or numerical defaults during optimization without reproducing the validation evidence.
- GPU work is optional and must be justified by crossover data.

## Work plan

1. Profile Fe/Ni production runs for all backends and integration modes.
2. Report time/memory fractions: H/G construction, GF contraction, energy integration, R/q FT, Dyson, spectral analysis.
3. Design hierarchical MPI decomposition preserving useful caches/reuse.
4. For r-GF, benchmark R-block/energy-node parallelism plus batched q FT against q-parallel duplication.
5. For k-GF/eigenpairs, benchmark k, q and omega decompositions.
6. Convert hot orbital contractions to batched BLAS where profitable.
7. Remove inner-loop allocations and redundant GF evaluations.
8. Add performance regression benchmarks with scientific checksums.
9. Prototype GPU only for sufficiently large/batched kernels and document CPU/GPU crossover.

## Acceptance checklist

- [x] Profiles exist for all production backends.
- [x] MPI strategy is justified by measured reuse/memory behavior.
- [x] r-GF q-amortization is preserved.
- [x] Optimized results reproduce scientific checksums.
- [x] Speedups and memory costs are reported.
- [x] Any GPU path has measured crossover evidence or is explicitly deferred.

## Required evidence / deliverables

Commit a performance report and benchmark data alongside optimization code.

## One-line commit message

`tddft: optimize validated response backends`
