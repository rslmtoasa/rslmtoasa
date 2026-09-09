# TDHYG-01 — Post-validation TD-DFT hygiene and diagnostics cleanup

## Prerequisite

Do not execute until TDVAL-01 is green. This task must not mask or substitute for physics validation.

## Goal

Clean up remaining misleading metadata, dead reference code, loose backend gates, sentinel rendering, and planner/reporting issues identified during the audit.

---

## Scope

### A. Dead "reference" implementations
Audit `if (.false.)` reference/oracle blocks in TD-DFT chi0/Xi code.

Choose one:
1. promote a genuinely independent reference path behind a test-only/runtime test switch and compare it;
2. remove dead code and correct metadata so production does not claim it is running an "explicit reference" that is unreachable.

Do not retain misleading provenance.

### B. Backend equivalence tolerance
The old dynamic cross-backend tolerance around 5% is too loose for mathematically equivalent matched sums.

After TDVAL-01 establishes observed errors:
- set separate tolerances for exact/deterministic equivalence vs broader physical comparisons;
- document why each tolerance exists.

Do not choose the new tolerance merely to fit the current largest error.

### C. Sentinel printing
Keep robust internal sentinels if useful, but output:
- `not applicable`;
- `unset`;
or omit the field when a backend does not use it.

Do not print `huge()` as if it were a physical energy window.

### D. MPI/work-plan metadata
Audit whether an eigenpair backend incorrectly passes Green-function energy-point defaults (e.g. 2001) into the planner.

Planner dimensions must correspond to actual work axes.

Add an MPI smoke/regression test if an MPI environment is available. If not, add unit coverage of planner dimensions and document the unexecuted MPI risk.

### E. temperature override provenance
Set `response_electronic_temperature_overridden` from an actual value/policy difference, not merely from "key appeared in namelist".

### F. formatting
Fix concatenated logical metadata such as `raw_identity_consistentF`.

### G. shifted-workset logging
If TD-DFT overrides symmetry reduction for finite q, state that explicitly in logs and metadata.

---

## Forbidden scope creep

- No new physical kernels.
- No change to Ward sign.
- No pair-potential normalization changes.
- No mode-finder physics changes.
- No GBT/SOC/Hubbard/CCOR expansion.

---

## Acceptance checklist

- [x] Dead reference/oracle code either becomes reachable and tested or is removed.
- [x] Backend metadata names the implementation actually executed.
- [x] Equivalence tolerances are physics/numerics justified.
- [x] `huge()` sentinels are not printed as physical provenance.
- [x] MPI planner axes match the selected backend.
- [x] Planner unit/MPI test added.
- [x] Temperature override flag has correct semantics.
- [x] Metadata formatting defects fixed.
- [x] Symmetry-reduction overrides are logged.
- [x] No physics result changed unexpectedly.
- [x] Full relevant TD-DFT regression suite passes.

## Required one-line commit message

`td-dft: clean validation metadata and planner hygiene`
