# C5.2 Part D — Stage RD status

**Historical-record supersession:** the C5.2 RD2 FAIL below records the original frozen Ni replay and remains unchanged as history. For the current C5.4 execution, the amendment-required offline eight-point fit and corrected per-state scalar-control classification supersede that RD2 verdict with **audited PASS**; see [C5.4 amended RD2](C5_4_T1_RD2.json) and the [sole consolidated audit](C5_4_AUDIT.md). This supersedes only the verdict under the C5.4 amendment; it does not rewrite the original RD2 report or its source evidence.

The table and historical interpretation below describe the C5.2 snapshot. Current C5.4 outcomes are linked in the final row and summarized in [C5.4 final status](C5_4_FINAL_STATUS.md).

| Item | Status | Record |
|---|---|---|
| RD1 local-moment correction | **PASS — unit GATE** | Two-site finite-H oracle: recursive-site corrected center and q2 match direct eigenvector results within 2.78×10⁻¹⁶; host-site swap differs by at least 3.84×10⁻⁴. |
| RD2 frozen Ni replay | **HISTORICAL FAIL — superseded for C5.4 by audited amended PASS** | Original replay remains frozen at [RD2 report](C5_2_PART_D_assets/ni_u3j07_rd2_replay/RD2_REPORT.md). The C5.4 amended fit and controls are in [T1](C5_4_T1_RD2.json); historical raw data and criteria are preserved. |
| RD3 comparisons | **RECORD — PARTIAL; GATES UNEVALUATED** | Ni RS local/spectral records are available; U=0 control mapping fatal leaves RD3(a) unevaluated, KS dependent gate unevaluated, no LDOS max computed. Fe local record and capture grid mismatch; B2 seed unavailable. |
| RD4 energy decision | **NOT MEASURED** | Release λ=0 crashed before SCF. Debug bounds same-deck reached 36 loops without an error; clean8b rejected the deck before SCF as input-incompatible. No energy derivative or RD4 gate result. [Evidence](C5_2_PART_D_assets/rd4_w1_measurement_unavailable/README.md). |
| RD5 adoption | **PARTIAL; no adoption** | No-U RS/KS identity and bulk behavior PASS. C5.4 explicit-spectral identity also PASS; required impurity `auto` cases remain unevaluated because there is no accepted Ni KS endpoint. See [C5.4 T5](C5_4_T5_RD5.json) and the historical [no-U results](C5_2_PART_D_assets/rd5_no_u_identity_checks/README.md). |

**Interpretation (7 lines):**

- RD1 remains the only passed production unit gate; no energy route or default changed.
- The original C5.2 RD2 verdict was a measured failure, reproduced byte for byte by its single authorized diagnostic replay; see [historical RD2 report](C5_2_PART_D_assets/ni_u3j07_rd2_replay/RD2_REPORT.md). C5.4's amended RD2 PASS is separately recorded in [T1](C5_4_T1_RD2.json) and confirmed in the [sole audit](C5_4_AUDIT.md).
- The historical non-degenerate spin-asymmetric proposal was not adopted. C5.4 used its registered eight-point per-state fit and unchanged allowance; no formula or tolerance changed.
- No-U RS/KS identity and bulk checks pass; the C5.4 explicit-spectral identity check passes. The accepted-KS endpoint is missing, so required impurity `auto` reproduction remains unevaluated and adoption/merge remain blocked under C5.3 item 2. RD4 is not measured and energy routing stays off.
- Ni U=3 RS bulk host is accepted at window 55–74; Ni U=0 KS continuation and U=3 KS host ended at finite B1 capacity guards, leaving the dependent Ni KS gate unevaluated.
- Ni spectral bulk and impurity are accepted under §2.1 (windows 28–47 and 3–22); RD3(d) is RECORD only. Ni impurity-minus-bulk d shifts are +4.57772565632e−7/+3.56986588734e−7, with no U=0 ratio or LDOS maximum. See [Ni RD3 record](C5_2_PART_D_assets/ni_rd3_ni_record/RD3_NI_RECORD.json).
- Fe U3.8/J1.1 local endpoint is accepted at window 1–20 under §2.1; one-step capture grids differ, so d-LDOS is unavailable. B2FeCo has no eligible accepted restart. See [Fe/B2 record](C5_2_PART_D_assets/fe_rd3_b2_record/RD3_FE_B2_RECORD.json). RD2 does not authorize criterion changes.

Current C5.4 status: audited RD2 PASS, RD3(a) PASS, RD3(b) UNEVALUATED, RD4 UNEVALUATED with energy OFF, RD5 explicit-spectral identity PASS, both `auto` endpoint checks UNEVALUATED, no adoption and no merge. See [C5.4 final status](C5_4_FINAL_STATUS.md) and [audit](C5_4_AUDIT.md).
