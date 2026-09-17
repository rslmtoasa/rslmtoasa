# RS-LMTO-ASA reciprocal TD-DFT validation — Luna prompt pack

Branch: `fable_v4`  
Pinned starting commit for this pack: `795b1e488d`  
Campaign name: `TDVAL-K`  
Date: 2026-09-12

## Purpose

This pack deliberately changes the order of work.

The native real-space GF TD-DFT backend is **deferred**.  RSGF capability levels
R0–R2 remain accepted within their documented finite/provider scopes, but
TDRUN-02 / R3 and native-material R4 are not on the present critical path.

The present goal is to establish a trustworthy **reciprocal-space TD-DFT**
framework first:

```text
accepted SCF state
    ↓
H(k), exact k+q endpoints
    ↓
LR-05 Pauli transition vertices
    ↓
LR-06 Lehmann chiKS  ←→  LR-GF-02 reciprocal-GF chiKS
    ↓
KXC-01 direct ALSDA / GSR-01 independent sum-rule route
    ↓
TDDY-01 Dyson + loss
    ↓
bcc Fe validation
    ↓
fcc Ni validation
```

The physics is specified in `00_MASTER_BLUEPRINT.md`.  Luna is an implementation
worker only.  It must not invent response equations, signs, factors, acceptance
thresholds, physical explanations, or literature fits.

## Execution order

Run one slice at a time:

1. `01_TDVK-00_PREFLIGHT.md`
2. `02_TDVK-01_BACKEND_CROSSCHECK_DIAGNOSTIC.md`
3. `03_TDVK-02_FE_REFERENCE_GAMMA.md`

**STOP and return evidence to the orchestrator.**

Only after explicit approval continue:

4. `04_TDVK-03_FE_BACKEND_CLOSURE.md`
5. `05_TDVK-04_FE_FINITE_Q_COVARIANCE.md`
6. `06_TDVK-05_FE_CONVERGENCE.md`

**STOP and return evidence to the orchestrator.**

Only after explicit approval continue:

7. `07_TDVK-06_STATIC_GOLDSTONE_INTERACTIONS.md`

**STOP and return evidence to the orchestrator.**

Only after explicit approval continue:

8. `08_TDVK-07_FE_DYSON_LOSS.md`
9. `09_TDVK-08_NI_REPEAT.md`
10. `10_TDVK-09_EVIDENCE_HANDOFF.md`

## Global Luna rule

When a requested numerical or material check fails:

- do **not** patch the physics;
- do **not** tune signs, factors, XC fields, eta, moments, or kernel scales;
- do **not** switch circular channel because another one “looks better”;
- do **not** enable Goldstone correction to hide a failure;
- record the failure exactly;
- identify the owning layer only if the evidence makes that mechanical;
- stop the slice and return the evidence.

Each prompt supplies a one-line commit message.  Do not combine slices in one
commit unless explicitly asked.
