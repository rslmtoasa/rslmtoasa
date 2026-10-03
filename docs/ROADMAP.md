# RS-LMTO-ASA roadmap

Compact index of the B-item plans. Statuses were checked against the tree on
2026-10-03; the Evidence column says where each status comes from. A status
with no evidence pointer was carried over from the former plan index and has
not been re-checked.

| ID | Title | Status | Evidence / plan |
|---|---|---|---|
| B1 | GBT fix + frozen magnons | Operator-level gates pass; G9 failed; converged small-q curvature not established. Single-sublattice frozen magnon validated; multi-sublattice acoustic branch open | `docs/DECISIONS.md` (GBT, GBT_REAUDIT); [VAL-16](validation/VAL-16_GBT_SUPERCELL.md), [VAL-17](validation/VAL-17_GBT_HARMONIC_GOLDSTONE.md); `tests/KNOWN_ISSUES.md` |
| B2 | k-space Green functions, two backends | Done, both gates signed | [G-B2-1](validation/B2_GATE_G-B2-1.md), [G-B2-2](validation/B2_GATE_G-B2-2.md) |
| B3 | Bloch spectral functions | Done; one unit test (`UnitBsfSumRule`), no example case or keyword page | — |
| B4 | Batched GPU eigensolver | Not started | [plan](roadmap/B4_gpu_kspace.md) |
| B5 | Route-agnostic post-processing via Lehmann | Done | [route-agnostic estimators](validation/route_agnostic_estimators.md) |
| B6 | Surface electrostatics | Done at narrower scope; gate open | `UnitMadl2dExhGuard`; [VAL-15](validation/VAL-15_MULTILAYER_VACUUM_ELECTROSTATICS.md) |
| B7 | Interfaces and vacuum leads | Done, with open items in `tests/KNOWN_ISSUES.md` | [VAL-14](validation/VAL-14_CU111_INTERFACE_MAPPING.md), [VAL-15](validation/VAL-15_MULTILAYER_VACUUM_ELECTROSTATICS.md) |
| B8 | k-space CPA + DLM | Not started | [plan](roadmap/B8_cpa_dlm_kspace.md) |
| B9 | Real-space CPA / DLM | Not started | [plan](roadmap/B9_rs_cpa_dlm.md) |
| B10 | DMFT self-energy provider API | Not started | [plan](roadmap/B10_dmft_sigma_provider.md) |
| B11 | Transverse spin response (RPA/ALDA χ, magnons) | Restarted 2026-10 on `fable_v4b`: Juelich-d k-space baseline for bcc Fe, Stage 0. The first campaign (2026-08-09 to 2026-10-03) closed without a validated dispersion | `docs/DECISIONS.md` closing entry; implementation spec held by the developer (not committed) |
| B12 | Electron-phonon / electron-magnon couplings | Not started; the e–magnon part depends on B11 | [plan](roadmap/B12_couplings.md) |
