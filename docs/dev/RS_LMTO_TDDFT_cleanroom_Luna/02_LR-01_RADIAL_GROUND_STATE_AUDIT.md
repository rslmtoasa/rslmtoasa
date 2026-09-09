# LR-01 — Audit radial ground-state magnetization and XC magnetic field

## Nature
BLOCKER AUDIT. No TD-DFT implementation.

## Goal
Prove whether converged RS-LMTO-ASA exposes the actual radial quantities required by the published ASA formulations:
\[
m_s(r),\qquad B_{{xc},s}(r).
\]

The legacy spherical charge density `f(r)` and radial Vxc are clues only.

## Literature targets
- BES-02 uses local ground-state `Bxc(r)/m(r)`.
- LCMM-01/02 uses actual KS magnetic field and magnetization density in the static sum rule.

## Required source trace
Trace:
1. spin-resolved radial density;
2. XC functional call;
3. spin-resolved radial XC potential;
4. converged spherical KS potential;
5. compression into LMTO potential parameters.

Inspect legacy atomic/SCF routines carefully; do not refactor them.

Identify exact arrays/routines for:
- radial mesh;
- charge density;
- magnetization/spin densities;
- Vxc;
- spin-resolved Vxc or equivalent;
- effective spin splitting;
- units.

## Derivation A — `f(r)` normalization
Prove whether stored radial quantity is:
- n(r)
- r^2 n(r)
- 4pi r^2 n(r)
- logarithmic-mesh weighted
- another convention.

Reproduce integrated electron count independently.

## Derivation B — m(r)
Derive exact spin-density/magnetization-density convention.

Integrate it and reproduce the code's local magnetic moment.

## Derivation C — Bxc(r)
From the actual SCF XC potential derive exactly:
\[
B_{xc}(r)=?
\]

Resolve:
- sign;
- factor 1/2;
- muB;
- Ry/Ha.

Do not infer from `C_up-C_down`, `VXC0SP`, pair tangents or another compressed proxy unless exact equality is proved from the ground-state code.

## Derivation D — XC provenance
For a supported LSDA/libXC choice, prove ground-state Vxc and later response kernel can use the same selected XC functional.

## Data lifetime
Determine whether exact converged radial quantities survive to the response stage or can be snapshotted exactly before being overwritten.

## Numerical oracle
On a simple converged collinear magnetic case report:
- integrated charge;
- integrated moment;
- selected radial samples of spin-resolved Vxc;
- derived Bxc(r);
- pointwise consistency.

## PASS condition
PASS only if radial mesh, m(r), Bxc(r), units, XC provenance and converged-state lifetime are exact.

Otherwise BLOCKED.

## Deliverable
`docs/LR_RADIAL_GROUND_STATE_AUDIT.md`

## Checklist
- [ ] f(r) semantics proved
- [ ] charge integral verified
- [ ] m(r) derived
- [ ] moment integral verified
- [ ] true radial XC source identified
- [ ] Bxc sign/factor derived
- [ ] Ry/Ha/muB mapped
- [ ] same-XC provenance demonstrated
- [ ] converged-state lifetime established
- [ ] no compressed-parameter surrogate
- [ ] PASS or BLOCKED declared

## Commit if PASS documentation/accessor evidence is added
`docs: establish radial ground-state response quantities`
