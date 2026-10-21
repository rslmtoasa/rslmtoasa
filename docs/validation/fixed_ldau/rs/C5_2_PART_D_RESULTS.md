# C5.2 Part D — Stage RD status

| Item | Status | Record |
|---|---|---|
| RD1 local-moment correction | **PASS — unit GATE** | Two-site finite-H oracle: recursive-site corrected center and q2 match direct eigenvector results within 2.78×10⁻¹⁶; host-site swap differs by at least 3.84×10⁻⁴. |
| RD2 frozen Ni replay | **FAIL — GATE** | Carrier ratios and q2′ pass; scalar-shift deviation 0.0026958204 Ry exceeds 0.00006024224 Ry allowance; spin-swap separation 0.0061524 Ry is below 10× carrier error 0.0615241 Ry. One authorized debug replay is byte-identical. No adoption. |
| RD3 comparisons | **RECORD — PARTIAL; GATES UNEVALUATED** | Ni RS local/spectral records are available; U=0 control mapping fatal leaves RD3(a) unevaluated, KS dependent gate unevaluated, no LDOS max computed. Fe local record candidate and capture grid mismatch; B2 seed unavailable. |
| RD4 energy decision | **PREPARED; NOT RUN** | Five-deck W1 inputs and focused unit PASS; execution record pending audit. No lambda solver runs or energy/derivative result. |
| RD5 adoption | **NOT STARTED** | No default, `auto`, or energy-routing changes. |

**Interpretation (7 lines):**

- RD1 remains the only passed production unit gate; no energy route or default changed.
- RD2 is a measured gate failure, reproduced byte for byte by its single authorized diagnostic replay; see [RD2 report](C5_2_PART_D_assets/ni_u3j07_rd2_replay/RD2_REPORT.md).
- The proposed non-degenerate spin-asymmetric state selection addresses spin-swap identifiability only; the scalar-shift failure remains unresolved and its registered allowance is unchanged. The proposal is not adopted.
- C5.3 §2.5 directs independent work to continue with production behavior unchanged.
- Ni U=3 RS bulk host is accepted at window 55–74; Ni U=0 KS continuation and U=3 KS host ended at finite B1 capacity guards, leaving the dependent Ni KS gate unevaluated.
- Ni spectral bulk and impurity are accepted under §2.1 (windows 28–47 and 3–22); RD3(d) is RECORD only. Ni impurity-minus-bulk d shifts are +4.57772565632e−7/+3.56986588734e−7, with no U=0 ratio or LDOS maximum. See [Ni RD3 record](C5_2_PART_D_assets/ni_rd3_ni_record/RD3_NI_RECORD.json).
- Fe U3.8/J1.1 local candidate window 1–20 remains pending independent acceptance; one-step capture grids differ, so d-LDOS is unavailable. B2FeCo has no eligible accepted restart. See [Fe/B2 record](C5_2_PART_D_assets/fe_rd3_b2_record/RD3_FE_B2_RECORD.json). RD2 does not authorize criterion changes.

Passed: RD1 unit GATE.
Recorded: RD2 measured FAIL after one byte-identical debug replay; no adoption; Ni RD3(d) RECORD only.
Recorded: Ni local/spectral impurity windows 4–23 and 3–22; Fe local candidate window 1–20 pending acceptance.
Unavailable: Ni U=0 mapping controls fatal; Ni KS hosts finite-guarded; Fe capture grids differ; B2 seed absent. RD5 blocked by RD2 FAIL.
Stopped: no active STOP; historical C5.2 A4(iv) remains superseded by C5.3 A4(v).
Next: RD4 W1 remains unlaunched pending execution-preflight audit.
