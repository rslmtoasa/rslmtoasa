# TDDFT Implementation Task Tracker

Dependencies:

```text
TDDFT-00
   |
TDDFT-01
   |
TDDFT-02
   |
TDDFT-03
   |
TDDFT-04
   |
TDDFT-05
   |
TDDFT-06
   |
TDDFT-07  <-- transverse scientific milestone
   |
TDDFT-08  <-- local longitudinal / LLB-input milestone
   |
TDDFT-09  <-- full 4C/SOC/non-collinear
   |
TDDFT-10  <-- GF backend
   |
TDDFT-11  <-- performance/GPU hardening
```

Some preparatory work can overlap, but scientific signoff should respect the gates.

## Status

- [ ] TDDFT-00 — conventions + XC-response plumbing
- [ ] TDDFT-01 — arbitrary-k eigenpair API
- [ ] TDDFT-02 — response basis + transition vertices
- [x] TDDFT-03 — bare transverse KS susceptibility
- [ ] TDDFT-04 — ALSDA kernel + Goldstone diagnostics
- [ ] TDDFT-05 — Dyson + Xi + modes + damping
- [ ] TDDFT-06 — production wiring + I/O + MPI over q
- [ ] TDDFT-07 — transverse validation campaign
- [ ] TDDFT-08 — longitudinal susceptibility + relaxation
- [ ] TDDFT-09 — full four-component/SOC/non-collinear
- [ ] TDDFT-10 — Green-function KS-response backend
- [ ] TDDFT-11 — optimized/GPU implementation

## Cross-project gates

- [ ] Corrected GBT branch available for TDDFT-07 comparison.
- [ ] `Jij`-derived stiffness route independently validated.
- [ ] Ground-state XC convention exposed cleanly enough for the response kernel.
- [ ] SOC physical-gap policy documented before TDDFT-09.
- [ ] Target LLB equation selected before defining an `alpha_parallel` conversion.

## Mandatory process rule for every task

- [ ] Inspect current branch before coding.
- [ ] Do not assume stale B11 statements are true.
- [ ] Tick prompt checkboxes as work completes.
- [ ] Run focused tests.
- [ ] Run relevant non-TDDFT regressions.
- [ ] Report unresolved physics/convention issues.
- [ ] Commit once with a single-line task-specific message.
