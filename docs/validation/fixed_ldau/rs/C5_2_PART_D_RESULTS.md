# C5.2 Part D — Stage RD status

| Item | Status | Record |
|---|---|---|
| RD1 local-moment correction | **PASS — unit GATE** | Two-site finite-H oracle: recursive-site corrected center and q2 match direct eigenvector results within 2.78×10⁻¹⁶; host-site swap differs by at least 3.84×10⁻⁴. |
| RD2 frozen Ni replay | **FAIL — GATE** | Carrier ratios and q2′ pass; scalar-shift deviation 0.0026958204 Ry exceeds 0.00006024224 Ry allowance; spin-swap separation 0.0061524 Ry is below 10× carrier error 0.0615241 Ry. One authorized debug replay is byte-identical. No adoption. |
| RD3 comparisons | **IN PROGRESS** | RS Ni U=3 host accepted; Ni U=0 KS continuation and Ni U=3 KS host stopped by finite B1 capacity guards, so dependent KS gate is unevaluated. Ni spectral host passed §2.1 raw-capture audit at iter47; RD3(d) impurity comparison remains unrun. |
| RD4 energy decision | **NOT MEASURED** | No Ni W1 energy or Hellmann–Feynman derivative runs. |
| RD5 adoption | **NOT STARTED** | No default, `auto`, or energy-routing changes. |

**Interpretation (7 lines):**

- RD1 remains the only passed production unit gate; no energy route or default changed.
- RD2 is a measured gate failure, reproduced byte for byte by its single authorized diagnostic replay; see [RD2 report](C5_2_PART_D_assets/ni_u3j07_rd2_replay/RD2_REPORT.md).
- The proposed non-degenerate spin-asymmetric state selection addresses spin-swap identifiability only; the scalar-shift failure remains unresolved and its registered allowance is unchanged. The proposal is not adopted.
- C5.3 §2.5 directs independent work to continue with production behavior unchanged.
- Ni U=3 RS bulk host is accepted at window 55–74; Ni U=0 KS continuation and U=3 KS host ended at finite B1 capacity guards, leaving the dependent Ni KS gate unevaluated.
- Ni spectral bulk is accepted at iteration 47 under §2.1 after independent raw-capture audit; the RD3(d) impurity comparison remains unrun.
- Fe RECORD work and other independent RD items continue; RD2 failure does not authorize changing any criterion.

Passed: RD1 unit GATE.
Recorded: RD2 measured FAIL; one byte-identical debug replay; Ni RS host accepted; Ni KS host attempts finite-guarded; Ni spectral accepted window 28–47.
Provisional: none. RD3(d) impurity comparison remains not yet run.
Stopped: no active STOP; historical C5.2 A4(iv) remains superseded by C5.3 A4(v).
Next: continue independent Stage RD items under §2.5; retain RD2 FAIL/no adoption.
