# RSGF-CLOSE-01 — native-RSGF capability closure audit

**Audit date:** 2026-09-11
**Branch:** `fable_v4`
**Starting HEAD:** `8338ed799a9d62ee31cdfda32add56c932f35299`
**Nature:** no-code capability closure audit

This audit reconciles the historical `BLOCKED` statement in LR-REP-00 with the
subsequent RSGF-00 and RSGF-01R work. It does not implement a driver or change
response physics, Green-function equations, radial augmentation, recursion,
Chebyshev numerics, kernels, Goldstone logic, or Dyson code.

The starting worktree had no tracked modifications. The only untracked files
were the supplied post-LR03 prompt-pack files. The audit itself makes
documentation changes only.

## 1. Preflight verdicts

The complete current evidence documents were read. The verdicts below are
scope-qualified; no status is inferred from a filename:

| evidence | final verdict at preflight | interpretation for this audit |
| --- | --- | --- |
| [`LR_RS_GF_REPRESENTATION_AUDIT.md`](LR_RS_GF_REPRESENTATION_AUDIT.md) | historical native RS response **BLOCKED** | audit the clause against later owners; LR-REP itself certifies coefficient-space representation semantics |
| [`RSGF_ENDPOINT_AUGMENTATION.md`](RSGF_ENDPOINT_AUGMENTATION.md) | **PASS** for finite endpoint augmentation | closes the missing radial/Pauli endpoint seam within the stated baseline |
| [`TDDFT_RS_GF_BACKEND.md`](TDDFT_RS_GF_BACKEND.md) | **PASS** for finite/provider response evidence; production integration not claimed | closes R2 at its documented finite/provider scope; leaves R3 pending |
| [`TDDFT_PRODUCTION_DRIVER.md`](TDDFT_PRODUCTION_DRIVER.md) | **PASS** for reciprocal collinear lifecycle | confirms an existing common post-SCF driver for the future R3 registration |
| [`TDDFT_COLLINEAR_REVALIDATION.md`](TDDFT_COLLINEAR_REVALIDATION.md) | historical TDVAL-01R **BLOCKED** | its native-RSGF preflight stop is reconciled here; Fe/Ni R4 validation remains pending |

The prerequisite capabilities are therefore not globally blocked. The old
statement is historically correct for LR-REP-00’s own audit point, but it is
superseded for R0–R2 by the later endpoint and bare-response tasks.

## 2. LR-REP-00 blocker clauses and current disposition

LR-REP-00’s response blocker was extracted clause by clause from its live
representation audit:

| LR-REP-00 blocker | Required capability | Later owner | Current evidence | Status |
| --- | --- | --- | --- | --- |
| `green%auxiliary_gij` is endpoint scaling, not physical radial/Pauli augmentation | Two-endpoint map from coefficient GF to Pauli radial/angular GF | RSGF-00 | [`RSGF_ENDPOINT_AUGMENTATION.md`](RSGF_ENDPOINT_AUGMENTATION.md); `augment_lr_gf_endpoint` | **CLOSED** |
| RS `gij/gji` had no equivalent map through `lmto_radial_basis`, Pauli operator, and LR-04 space | A callable endpoint adapter accepting native directed blocks and accepted radial bases | RSGF-00 | `augment_lr_gf_endpoint_pair` is public and is called by RSGF-01R | **CLOSED** |
| A screened/coefficient block could not be compared directly with the reciprocal physical response | Keep coefficient-space semantics separate from physical augmentation and use the certified endpoint seam | LR-REP-00 + RSGF-00 | coefficient-space bridge in `UnitLrRsGfRepresentation`; endpoint branches in `UnitLrGfEndpointAugmentation` | **CLOSED** for the supported baseline |
| No native RS susceptibility implementation existed at the time of LR-REP-00 | Full real-space GF bubble, q phase, and LR-04 canonical response assembly | RSGF-01R | [`TDDFT_RS_GF_BACKEND.md`](TDDFT_RS_GF_BACKEND.md); `evaluate_lr_rs_gf_susceptibility` | **CLOSED** for finite/provider validation |
| Native recursion/Chebyshev convergence was outside LR-REP-00 | Provider-specific convergence evidence | RSGF-01R | Current RSGF-01R document exposes separate controls but the focused oracle uses the exact dense provider | **DEFERRED OUTSIDE BASELINE** |
| Production material lifecycle was not owned by LR-REP-00 | Post-SCF backend registration and accepted-state handoff | TDRUN-02 | [`TDDFT_PRODUCTION_DRIVER.md`](TDDFT_PRODUCTION_DRIVER.md) is reciprocal-only | **PRODUCTION INTEGRATION PENDING** |

The last two rows are not reasons to keep R0–R2 blocked. They delimit RSGF-01R’s
finite/provider claim and the next production task respectively.

## 3. What LR-REP-00 actually certifies

LR-REP-00 owns Level R0 representation semantics, not the response stack. Its
certified content is:

- directed coefficient-space `gij`/`gji` orientation and site-local ordering;
- the reciprocal/native four-phase coefficient-space bridge;
- `sqrt(Delta)` endpoint scaling as an auxiliary representation operation;
- screened/path-operator transformations, including onsite additive terms and
  their inverse; and
- the boundary that a coefficient-space GF must not be called a physical
  radial/Pauli GF without a separate endpoint map.

LR-REP-00 does not own or retrospectively need to certify:

- physical radial/Pauli endpoint augmentation, which is RSGF-00;
- the response bubble or q assembly, which is RSGF-01R; or
- production material driver selection, which is TDRUN-02.

The old `BLOCKED` wording was correct because the RSGF-00 seam and RSGF-01R
service did not yet exist. It must remain in the historical portion of
`LR_RS_GF_REPRESENTATION_AUDIT.md`.

## 4. RSGF-00 endpoint seam

The live implementation in
[`source/lr_gf_endpoint_augmentation.f90`](../source/lr_gf_endpoint_augmentation.f90)
accepts a coefficient-space block and explicit `h^gamma` provenance. It stores
the four endpoint branches:

```text
Phi G Phi^dagger
Phidot (hG) Phi^dagger
Phi (Gh) Phidot^dagger
Phidot (hGh) Phidot^dagger
```

The implementation forms the single-sided resolvent terms from the endpoint
energy difference and retains the required contacts: onsite `-I` terms in
`hG` and `Gh`, onsite `-D` in `hGh`, and `-hgamma` for both onsite and offsite
blocks. The pair entry point requires the reverse directed GF and reverse
`hgamma` when no effective-Hamiltonian provider is supplied; it does not infer
the reverse block by conjugation.

The endpoint test is not a name-only check. `UnitLrGfEndpointAugmentation`
compares the factored adapter against independent dense and Lehmann spatial
constructions at four complex energies, tests an effective-action provider,
and requires negative controls for omitted `Phidot`, onsite contacts, and the
offsite `-hgamma` term. The direct executable reported representative errors:

```text
dense offsite augmented-GF equivalence  <= 5.90E-18
dense onsite augmented-GF equivalence   <= 5.73E-17
Lehmann augmented-GF equivalence         <= 2.07E-16
omitting onsite -I contacts             = 1.56E-01
omitting offsite -h_ab                   = 5.47E-04
```

This closes R1 for the supported orthogonal, collinear, no-SOC, Pauli large-
component `sp`/`spd` baseline. It does not claim SOC, generalized overlap,
noncollinearity, additive operators, or `spdf`.

## 5. RSGF-01R composition trace

The passing RSGF-01R candidate path composes the certified seams in the live
code:

```text
lr_rs_gf_provider%get_pair
    -> coefficient-space directed G_ij/G_ji and hgamma blocks
    -> augment_lr_gf_endpoint_pair
    -> four RSGF-00 endpoint branches
    -> accumulate_native_pair_bubble
    -> response_real_space_phase
    -> response_raw_to_canonical
    -> LR-04 canonical chiKS
```

In `source/lr_rs_gf_susceptibility.f90`, the provider is called for the
retarded/advanced base blocks and both frequency-shifted blocks. Each call is
passed through `augment_lr_gf_endpoint_pair` before
`accumulate_native_pair_bubble`. The q phase is obtained from the existing
`response_real_space_phase`; the result is converted by
`response_raw_to_canonical` and labelled
`LR-04 canonical right-weighted B=chi_raw*W`.

The finite test distinguishes candidate and reference paths:

- **candidate:** `lr_rs_dense_gf_provider` satisfies the public provider
  interface, supplies coefficient GF blocks, and then traverses the production
  RSGF-00 augmentation and RSGF-01R bubble code;
- **references:** LR-06 spectral and LR-GF-02 reciprocal-GF services are called
  separately for comparison; they are not used inside the native candidate;
- **scope:** the dense provider is an exact finite oracle. The test is not a
  block-recursion or Chebyshev material convergence run.

`UnitLrRsGfSusceptibility` reported:

```text
coefficient-GF versus spectral-resolvent error       = 0.0000E+00
native versus reciprocal-GF full-response (abs, rel) = 1.1154E-21  7.2566E-14
nonzero-q translation phase assembly error           = 2.6610E-22
```

The native-versus-spectral absolute difference was `1.5371E-08` on the tiny
finite fixture. The response is compared in the full radial/angular space; the
selected full-space entry and the q-phase fixture are both exercised. This
closes R2 exactly at the finite/provider scope claimed by RSGF-01R.

## 6. Capability levels R0–R4

| level | meaning | owner | current status |
| --- | --- | --- | --- |
| R0 — representation | native GF has a certified coefficient-space representation contract | LR-GF-01 / LR-REP-00 | **CLOSED** |
| R1 — endpoint augmentation | coefficient GF can be promoted to the certified Pauli radial/angular one-electron GF | RSGF-00 | **CLOSED** for the certified finite baseline |
| R2 — bare response | native GF produces the same finite-basis radial/angular `chiKS` as the reciprocal references within controlled provider/integration error | RSGF-01R | **CLOSED** for finite/provider validation |
| R3 — production lifecycle | a converged material can select native RSGF through the normal post-SCF driver | TDRUN-02 | **PENDING** |
| R4 — material validation | Fe/Ni production calculations validate the native route | TDVAL-01R | **PENDING** |

R0–R2 must not be relabelled as R3 or R4. The absence of production-driver
integration and Fe/Ni results is expected at this audit point.

## 7. Authoritative capability ledger

This table is the downstream gate. Later tasks should use it instead of
inferring eligibility from stale prose in historical audits:

| capability | owner | status | evidence |
| --- | --- | --- | --- |
| native coefficient GF semantics | LR-GF-01 | **CLOSED** | [`LR-GF-01_GF_CONTRACT_EVIDENCE.md`](LR-GF-01_GF_CONTRACT_EVIDENCE.md) and representation bridge |
| native GF representation conversion | LR-REP-00 | **CLOSED** | [`LR_RS_GF_REPRESENTATION_AUDIT.md`](LR_RS_GF_REPRESENTATION_AUDIT.md), coefficient-space transformations and bridge |
| Pauli radial endpoint augmentation | RSGF-00 | **CLOSED** for certified finite baseline | [`RSGF_ENDPOINT_AUGMENTATION.md`](RSGF_ENDPOINT_AUGMENTATION.md), `UnitLrGfEndpointAugmentation` |
| finite native bare response | RSGF-01R | **CLOSED** for finite/provider validation | [`TDDFT_RS_GF_BACKEND.md`](TDDFT_RS_GF_BACKEND.md), `UnitLrRsGfSusceptibility` |
| production post-SCF registration | TDRUN-02 | **PENDING** | [`TDDFT_PRODUCTION_DRIVER.md`](TDDFT_PRODUCTION_DRIVER.md) is currently reciprocal-only |
| Fe/Ni material validation | TDVAL-01R | **PENDING** | [`TDDFT_COLLINEAR_REVALIDATION.md`](TDDFT_COLLINEAR_REVALIDATION.md); no material rerun in this audit |

Provider-specific block-recursion/Chebyshev convergence remains a separately
labelled numerical requirement before any material result can be promoted. It
does not invalidate the finite/provider R2 closure or authorize a material
claim.

## 8. Focused verification

The existing tests were run without source changes:

```text
build/bin/UnitLrRsGfRepresentation       -> RESULT: PASS
build/bin/UnitLrGfEndpointAugmentation   -> PASS (dense, spectral, provider, negative controls)
build/bin/UnitLrRsGfSusceptibility       -> PASS (native coefficient GF, augmentation, full response, q phase)
```

These tests support the ledger’s R0–R2 closure. They do not test production
backend selection, accepted material state mutation, Fe/Ni physics, or
literature agreement.

## 9. Historical provenance and cross-references

The original LR-REP-00 report remains intact as a record of its earlier audit
point. A dated closure section has been appended to that document. The RSGF-01R
and TDVAL documents now link this ledger to distinguish finite/backend closure
from production/material work.

No source code or test implementation was changed by RSGF-CLOSE-01.

## 10. Eligibility conclusion

ELIGIBLE FOR TDRUN-02
