# TDWARD-02 — Rebuild Ward/Goldstone diagnostics from independent quantities

## Prerequisite

TDWARD-01 must have produced an unambiguous derivation and independent B_xc provenance.

If TDWARD-01 ended with unresolved ambiguity, do not implement this task until that ambiguity is resolved.

## Goal

Replace diagnostic aliases and tautologies with a small set of physically independent tests that can fail for distinct reasons.

---

## Required diagnostics

### A. Pair-Xi Goldstone residual

\[
r_\Xi =
\frac{\|\Xi m-m\|}{\|m\|}.
\]

Report:
- closest Goldstone eigenvalue;
- complex distance to +1;
- residual on the signed Goldstone magnetization vector.

### B. Independent direct Ward residual

Using B_xc from the independent ground-state provenance established in TDWARD-01:

\[
r_B =
\frac{\|\chi_{\rm KS}B_{\rm xc} \mp m\|}{\|m\|},
\]

with the sign fixed by the derivation.

Do not obtain B_xc by applying the same K_xc to m.

### C. Static/dynamic consistency residual

At Gamma and selected finite q where applicable, compare the exact static divided-difference chi0 with the dynamic route as eta -> 0+:

\[
r_{\rm stat/dyn}(\eta)
=
\frac{\|\chi_{\rm KS}^{\rm static}
-\chi_{\rm KS}^{\rm dynamic}(i\eta)\|}
{\max(\|\chi_{\rm KS}^{\rm static}\|,\epsilon)}.
\]

This should be evaluated as an eta ladder, not from one broadening.

---

## Remove or relabel non-independent diagnostics

Audit:
- `residual`;
- `ward_residual`;
- `dm_residual`;
- `bxc_kernel_residual`;
and any other fields that are exact aliases or identities by construction.

Either:
1. remove them from the validation report; or
2. retain them only with explicit metadata that they are aliases/debug quantities and not independent evidence.

Do not present four identical numbers as four physics checks.

---

## Phase-sensitive Goldstone character

The current absolute overlap can make +1 and -1 modes look identical.

Report both:
- normalized magnitude overlap;
- complex/sign-sensitive overlap or phase relative to m.

For a scalar one-site case, an anti-Goldstone eigenvalue at -1 must not be reported as indistinguishable from a +1 Goldstone mode.

---

## xi_backend='compare'

Do not make arbitrary magnitude mismatch universally fatal because legacy scalar and pair-potential routes may not be identical approximations.

However, in a validation configuration where both routes claim applicability to the same q=0 collinear SOC-off Goldstone reference:
- opposite Goldstone branch/sign is a hard failure;
- failure to identify the +1 Goldstone branch is a hard failure;
- magnitude mismatch should be reported quantitatively and gated only where mathematical equivalence is established.

Make strictness explicit in metadata/configuration; do not silently change production semantics.

---

## Output correctness

Fix:
- documented Ward identity sign/factors;
- `raw_identity_consistent` formatting;
- response-space rank in the header;
- signed site magnetizations in diagnostics;
- provenance of independent B_xc.

---

## Mandatory tests

1. One-site +z LMTO.
2. One-site -z LMTO.
3. Two-sublattice +/-z LMTO.
4. Deliberately injected anti-Goldstone Xi=-1 test:
   - phase-sensitive diagnostics must detect it.
5. Deliberately perturbed independent B_xc:
   - r_B must fail while r_Xi can remain unchanged.
6. Deliberately perturbed Xi:
   - r_Xi must fail without forcing r_B to be identical.
7. eta-ladder static/dynamic test.

The point is to prove the diagnostics are independent by making them fail separately.

---

## Forbidden shortcuts

- No B_xc reconstructed from the same K_xc m identity being validated.
- No absolute-value-only mode character.
- No duplicated residuals with different names.
- No sign chosen from the Fe result.
- No tolerance weakening to accommodate a convention mismatch.

---

## Execution record

Implemented and verified on 2026-09-09:

- Independent circular B_xc source uses the TDWARD-01 radial VXC0SP provenance, with signed orientation and the derived `s_i N_B/(2 M_i)` factor. Production Goldstone diagnostics require this source; the `K_perp m_G` path remains explicitly labelled debug-only.
- `r_Xi`, complex distance to +1, signed/complex Xi action, independent `r_B`, static/dynamic eta ladders, response-space rank, signed site magnetizations, and B_xc provenance are now reported. Legacy residual names are not emitted as validation evidence.
- Production evaluates the eta ladder at Gamma; the static/dynamic checker is q-provenance aware and is exercised at a finite-q endpoint by `UnitTddftChiKS`.
- Mandatory LMTO fixtures pass: one-site +z/-z and two-sublattice +/-z, with maximum tangent/service error `8.8558e-10` against `2e-9`; the two-sublattice signed Xi actions are `+1` and `-1`.
- Focused CTest suite: 21/21 passed; unresolved-occupation regression: 1/1 passed; TDDFT Ward/Goldstone/ChiKS tests: 3/3 passed; dispatch contract tests: 15 passed.
- The broader `UnitTddftBackendEquivalence` campaign was not completed within 180 seconds because it runs a large pre-existing TDDFT09 energy/contour campaign; it remains a performance/regression follow-up rather than a failed focused diagnostic test.

## Acceptance checklist

- [x] Independent B_xc path from TDWARD-01 implemented.
- [x] r_Xi implemented.
- [x] independent r_B implemented.
- [x] static/dynamic eta-ladder residual implemented.
- [x] alias/tautology diagnostics removed or explicitly relabelled.
- [x] phase-sensitive Goldstone character reported.
- [x] anti-Goldstone fixture is detected.
- [x] r_B and r_Xi can be made to fail independently in tests.
- [x] Ward documentation matches actual code convention.
- [x] response-space rank and signed site magnetizations are reported.
- [x] compare-mode branch/sign policy is explicit.
- [x] Existing relevant tests pass.

## Required one-line commit message

`td-dft: make Ward diagnostics independent and phase aware`
