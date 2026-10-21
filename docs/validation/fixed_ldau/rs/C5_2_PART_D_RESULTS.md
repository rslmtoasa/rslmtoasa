# C5.2 Part D — Stage RD status

| Item | Status | Record |
|---|---|---|
| RD1 local-moment correction | **PASS — unit GATE** | Two-site finite-H oracle: recursive-site corrected center and q2 match direct eigenvector results within 2.78×10⁻¹⁶; host-site swap differs by at least 3.84×10⁻⁴. |
| RD2 frozen Ni replay | **PENDING** | No RD2 impurity solver run; awaits the accepted Ni local impurity state. |
| RD3 comparisons | **PENDING** | No RD2-dependent Ni comparisons launched. |
| RD4 energy decision | **NOT MEASURED** | No Ni energy or Hellmann–Feynman derivative run. |
| RD5 adoption | **NOT STARTED** | No default, `auto`, or energy-routing changes. |

**Interpretation (5 lines):**

- The RD1 production path captures the exact inserted impurity d-block and the audited oracle passes.
- The production change is limited to the reviewed RD1 local RS impurity path; no physics formula, default, or energy route changed.
- C5.3 A4(v) supersedes the historical C5.2 A4(iv) stop; RD2–RD5 remain pending, with no new impurity campaign claimed.
- Separate bulk reference work generated an accepted Ni U=3 RS checkpoint at window 55–74; its solver was intentionally SIGTERM-stopped after preservation (exit 143, cleanup verified). A Ni U=0 KS run reached its 100-step cap unconverged; its beta=.05 continuation is prepared but not launched. These are bulk references, not RD2 impurity runs.
- Continue at RD2 after the accepted Ni impurity endpoint is available and separately audited.

Passed: RD1 unit GATE.
Recorded: two Ni bulk reference runs (U=3 RS accepted checkpoint, intentionally stopped; U=0 KS cap-limited), no RD2-dependent impurity run; RD2–RD5 pending.
Provisional: none.
Stopped: no active STOP; historical C5.2 A4(iv) is superseded by C5.3 A4(v).
Next: RD2 frozen replay on the accepted Ni local impurity state.
