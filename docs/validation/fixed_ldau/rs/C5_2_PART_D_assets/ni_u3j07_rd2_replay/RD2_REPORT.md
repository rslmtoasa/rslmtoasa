# RD2 Ni in Ni frozen replay — measured fail

**Verdict: GATE FAIL.** The accepted Ni U=3 eV, J=0.7 eV local impurity checkpoint was replayed against the frozen baseline using the pre-registered cases. One authorized diagnostic replay reproduced the original output byte for byte. No criterion, fit, tolerance, production behavior, or default was changed; RD2 is not adopted.

| Check | Measurement | Locked criterion | Result |
|---|---:|---:|---|
| Carrier ratio, spin 1 / 2 | 0.9914299182 / 0.9914299103 | each within 1% of 1 | PASS |
| q2′ minimum | 0.1378160163 Ry² | ≥ 0 | PASS |
| Scalar-shift deviation | 0.002695820416586253 Ry | ≤ 0.000060242239503755855 Ry | FAIL |
| Sign-wrong separation | 1.3962656345 Ry | > 10 × carrier error (0.0615241337 Ry) | PASS |
| Spin-swap separation | 0.0061524133746 Ry | > 10 × carrier error (0.0615241337 Ry) | FAIL |

The opposite-spin control is physically degenerate for this nearly spin-symmetric state: true versus swapped total Tr(VP) differs by 3.33×10⁻¹⁶ Ry, and the Vin spin-norm difference is 9.93×10⁻⁹. Its measured separation is below the registered threshold. The scalar-shift response also exceeds its registered allowance. These are measured failures, not grounds to relax either gate.

**Proposed criterion for a future separately authorized RD2 protocol:** require a pre-registered, non-degenerate spin-asymmetric accepted impurity state (with a lower-bound spin contrast established before replay), and use the existing spin-swap separation gate unchanged. This would make the wrong-spin control identifiable. Retain the scalar-shift allowance and all other gates as written. Validate the proposed state-selection rule independently before any new RD2 replay. This proposal is not adopted here.

Primary evidence: `replay/cases_005/rd2_cases_report.json` (SHA-256 `29539bcb06648f139719514e8b49b0c1d97efdd009ed8544f7e6e08ff1032153`); identical diagnostic replay output `replay/diagnostic_006/replay_006/driver.output.csv` (SHA-256 `49206e51fd5b8cd1bbab83d79dca7d47b74ef760725e8eff6fee3d616954c32e`). Full inputs, authorizations, failed preflights, preregistrations, captures, and raw archives are retained in this evidence package.
