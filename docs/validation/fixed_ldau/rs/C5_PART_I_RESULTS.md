# C5 Part I — KNOWN_ISSUES corrections

**Result:** the two stale Fe/Co entries at the end of `tests/KNOWN_ISSUES.md` were replaced by the four records required by C5 §3. Earlier, unrelated history was preserved. `C5_AMENDMENT.md` is an exact copy of the approved request.

| Entry | Status | Recorded result | Evidence |
|---|---|---|---|
| hcp Co impurity setup | OPEN | Debug bounds mismatch at `lattice_cluster.f90:996` (3 vs 50000); clean baseline SIGKILL before SCF; U=0 reproduces setup failure | hcp deck and C1 triage |
| Co bulk U=0 | OPEN | RB6 non-finite after about 19–40 iterations; reproduced on `8b8f7ef` | C1 RB6 results/provenance |
| Fe U>0 embedding offset | OPEN, low priority | −6.6×10⁻³ μB, about 735× U=0 control; likely rc mismatch; matched-cluster test blocked by entry 4; Ni clean | C1 RA4 results/provenance |
| RS rc=100 instability | OPEN | Fe impurity-style deck non-finite at 102 (U>0) / 34 (U=0) iterations; U=0 bulk not accepted after 600; reproduced on `8b8f7ef`; rc=140 INVALID_SETUP excluded | C1/C4 records |

**GATE verdicts:** Part I is documentation-only; no computational GATE was run. The C5 amendment defines the subsequent R1 GATE.

**Interpretation (3 lines):** The replacement reflects the C5 issue taxonomy and links the existing decks and evidence. No source code, test files, cleanup protocol, or campaign data were changed. No campaign was run.

Passed: Part I documentation corrections committed.
Recorded: four OPEN issues and their existing evidence links.
Provisional: none.
Stopped: no Part I GATE or campaign was run.
Next: Part II/R1 hcp-Co setup triage.
