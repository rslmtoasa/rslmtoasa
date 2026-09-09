# LR-03 — Lock TD-DFT conventions and normalization

## Prerequisites
LR-01 and LR-02 PASS.

## Nature
Derivation task. No production susceptibility/kernel.

## Goal
Create the one authoritative mapping between published equations and RS-LMTO-ASA definitions.

## Output
`docs/LR_TDDFT_CONVENTIONS.md`

## Required sections

### Density definitions
Define n, mx,my,mz and whether m is electron spin density, magnetic-moment density, or another quantity.

### Spin operators
Write exact I, sigma_x/y/z and circular operators. Derive every factor of two.

### Bxc
Derive exact relation among:
- spin-dependent KS potential;
- Bxc;
- magnetic moment;
- muB.

### ALSDA transverse kernel
Map BES-02 into RS-LMTO variables and units.

### ASA sum-rule normalization
Map LCMM-02/03/04 into the RS-LMTO radial/angular basis.
Show exactly where `4pi` appears or cancels.

### Fourier convention
Map RS-LMTO phases to BES-06; define q in direct and Cartesian reciprocal coordinates.

### Retarded chiKS
Fix:
- occupation difference;
- denominator;
- sign of i eta;
- source/measurement order.

### Units
Map Ry, Ha, eV, muB, length, susceptibility units and Kxc units.

### Goldstone rigid rotations
Write them in the approved radial/angular basis.

### Covariance
Derive q->-q, omega->-omega and global-spin-reversal identities.

## Falsification
Add tiny algebraic normalization tests independent of materials.

## PASS condition
No unresolved sign/factor/unit.

Any ambiguity => BLOCKED.

## Checklist
- [ ] density definitions fixed
- [ ] spin matrices/factors fixed
- [ ] Bxc fixed
- [ ] ALSDA kernel mapped
- [ ] 4pi mapped
- [ ] Fourier convention mapped
- [ ] retarded denominator fixed
- [ ] all units fixed
- [ ] rigid rotation represented
- [ ] covariance derived
- [ ] no sign chosen from Fe/Ni
- [ ] PASS or BLOCKED

## Commit
`docs: lock radial TDDFT conventions and normalization`
