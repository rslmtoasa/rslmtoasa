# KXC-01 — Implement the local radial ALSDA transverse XC kernel

## Prerequisites
LR-01 through LR-06 PASS.

## Goal
Implement the literature-locked local transverse ALSDA kernel in the exact RS-LMTO convention derived in LR-03.

## Suggested file
`source/tddft_alsda_kernel.f90`

## Physics
Start from BES-02 only after applying the exact LR-03 mapping.

Use:
- actual converged radial Bxc(r);
- actual converged radial m(r);
- same XC functional provenance as the ground state.

No compressed LMTO potential-parameter proxy.

## Spatial representation
Represent/project the local radial kernel exactly in the approved response basis.

Do not reduce it to one value per site.

## Low-m regions
Do not invent regularization.

Inspect whether m(r) crosses or approaches zero in regions carrying basis weight. If the formula becomes singular and the cited literature does not provide the treatment, declare BLOCKED.

## Tests
- pointwise Bxc/m samples;
- basis projection vs direct quadrature;
- units;
- global spin reversal;
- fatal XC-provenance mismatch;
- report raw `D=I-chiKS*Kxc` at q=0 without correction.

## Forbidden
- no scalar site kernel;
- no ad-hoc scaling;
- no pair/rotation tangent;
- no Goldstone correction.

## Checklist
- [ ] exact radial Bxc/m source
- [ ] LR-03 signs/factors
- [ ] same XC functional enforced
- [ ] basis representation exact
- [ ] low-m behavior literature-supported or BLOCKED
- [ ] raw q=0 Goldstone residual reported
- [ ] no empirical adjustment

## Commit
`td-dft: implement radial ALSDA transverse kernel`
