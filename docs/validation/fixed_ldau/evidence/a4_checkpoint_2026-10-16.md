# A4 campaign checkpoint — 2026-10-16

This additive checkpoint records the current gates and remaining releases. It
does not amend the user-owned protocol, A4's verbatim amendment, any frozen
preregistration, or any historical result. No commit or push has been made.

## Current worktree and frozen inputs

- Worktree: `/private/tmp/rslmto-ldau-m2-energy-worktree`
- Branch and base HEAD: `fix/ldau-m2-energy`, `f0ed469d29d8feb6e8569b103f10b9bc1bc748e4`
- A4 source-ready manifest v3: `evidence/a4_source_ready_manifest_v3.json`,
  SHA256 `db487ea375fb20dd025ed4d62b4ae880940b865f5fa16c3cefa95220985137ed`
- Matched serial Release executable:
  `/private/tmp/a4-m2-serial/bin/rslmto.x`, SHA256
  `536fa80a34c510782013a5608fc1220dca26b52a10189a2334e40ed43e12a0b4`
- Production `source/self.f90` SHA256:
  `2184df95a5560d65d730e6074618ec74210e20fb2eee309c8ac8be7bffb3e010`
- All 53 pre-A4 evidence entries still match the snapshot. The original
  user-owned protocol remains 27,152 bytes, SHA256
  `812bcafae105eabaea4ae0614f0b419e0d96638df37888342a8708a973f83980`.

## Completed A4 housekeeping and RS M1

A4.0(b) basis evidence identifies stored onsite d indices 5:9 as spherical
`m=-2..+2` after the Cartesian-to-spherical transform. The source comment cites
[`PHASE1.md`](../PHASE1.md) and
[`PROJECTOR_CONTRACT.md`](../PROJECTOR_CONTRACT.md); source and evidence hashes
are recorded in `LOCAL_MOMENTS.md` and the v3 source manifest.

A4.0(c) entropy unit validation passed. The finite-H oracle gives dimensionless
`S=2.1062540069077236` per cell and `T*S=0.0040021025481680765 Ry/cell`.
Candidate free energies report `F=E-T*S` in Ry/cell; no site entropy
decomposition is claimed. RS entropy is unavailable, not zero.

RS Gate M1 has an independently certified supplemental PASS under the
user-approved literal carrier and nonnegative-q2 checks. Its original v1
preregistration and response remain FAIL under their written residual-inclusive
point-fit criterion and are unchanged. The supplemental result is
`evidence/a4_rs_m1_supplemental_verdict.json` (SHA256
`d1403013e1a1ee59b4e839590f5281bb77cf6267855f7b7957bb40b0215d5095`); the
pre-certification draft is preserved byte-for-byte separately. No RS M2
measurement has been made.

## KS lambda-zero operational status

All eight registered W1/W2, local/spectral, standard/tight lambda-zero
processes exited with return code zero and have a terminal convergence marker,
trace, and final restart. Each runtime log reports fixed-U R05 as not applicable
at `U=0`; standard and tightened SCF `conv_thr` values are those recorded in
the preregistration. The operational audit deliberately contains no energy,
density-matrix, `g0`, or `F` comparisons:
`evidence/a4_m21_ks_lambda0_operational_audit.json`, SHA256
`70ceafefc1202809ad62a69203ebb9775319406e2d3aa869c80793ba36d12887`.
The batch manifest SHA256 is
`cfea7226f0e2bcbfae778d7a45a92205e6ff6dc70a01d66c2935a8f5f5bbb461`.

## Governing rule and unresolved proposal

The user-approved A4.1 identity gate remains in force: compare the
Richardson-extrapolated derivative with `Tr(WP(0))` using the preregistered
repeat scatter, truncation, trace uncertainty, and roundoff allowance; identity
bias is excluded. The proposed Richardson-aware replacement of the strict
pure-h² 4:1 resolution check remains `PROPOSED_NOT_APPROVED` in
`evidence/a4_m21_h2_resolution_criterion_proposal_NOT_APPROVED.md` (SHA256
`2613d8a48ff765dd6965ac5a494c6349b893d7bf89a1aa1416f3f48df1524b11`). Until
the user chooses, the existing strict 4:1 criterion governs.

The outstanding A1.3 decision is whether to retain that strict noise-only 4:1
check or approve the proposed resolved-leading-term rule
`|L|-u_L > 16(|H|+u_H)`, with `|L|>u_L`, alongside the same identity allowance.
The proposal's motivation is that a smooth `D(h)=g+10h^2+32h^4` has a
resolved h⁴ contribution that fails strict 4:1 within roundoff, although the
locked Richardson pair estimates that bias. Approval is not inferred from a
generic request to continue; the current strict rule remains governing until
an explicit user choice is recorded under A1.3. This proposal changes no
identity allowance and is not a post-outcome fit.

Tooling lock v8 is a pre-comparison analysis lock, not authorization to inspect
lambda-zero numerical responses. It binds the source-derived floating-point
and serialization error model, exact zero-based k/band domains, and synthetic
checks; its hash is
`f75e7d3aab6f251cc64e30ecfbfe2266bebdf44b51a4c0a7890a8b9d9f47eec8`. The
independent synthetic test passed; reviewer independently confirmed the
synthetic checks. Full certification of the draft tooling lock is pending. The
response analyzer is bound at SHA256
`6d439322b4ffe23b1360646445bc2851265c05cc731cbcb8ea149bae4132a7de`. The
reviewer identified an additional requirement to bind the raw-token extractor
and its source/capture provenance before comparing any lambda-zero values.
Therefore no `g0`, `u_g0`, derivative, or candidate comparison has been opened.

## Remaining stages

1. Freeze and independently review the complete analyzer plus raw-token
   extractor binding. Keep the lambda-zero artifacts immutable.
2. Resolve the outstanding h²-rule proposal through explicit user governance.
   Do not treat a general “continue” instruction as choosing it.
3. Only after both conditions above are satisfied, create the final
   response-analysis registration and perform the authorized lambda-zero
   comparison. No nonzero-λ case is currently released.
4. For KS A4.2, obtain its own preregistration and the required separate
   approval for transferring `N_W` into U/J derivative units before any U/J
   response. A4.3 remains unrun and is staged only as an unregistered
   `a4_m22_m23_campaign_template.json`.
5. KS M2 is not PASS, and no physical-energy decision is made. The `spectral`
   default and `physical_total_energy` routing are unchanged. No default or
   `physical_total_energy` change is permitted before a KS M2 PASS. RS receives
   no local default unless it passes its own M1 and M2.
6. Stop at the M2 gate. Any A4.4 routing, documentation update, restart/MPI
   demonstration, or PR preparation remains conditional on the gate and root
   reconciliation of the no-push instruction.
