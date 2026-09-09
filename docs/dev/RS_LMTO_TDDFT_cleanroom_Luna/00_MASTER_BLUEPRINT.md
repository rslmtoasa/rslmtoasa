# Master blueprint — literature-locked TD-DFT clean-room redevelopment

## Objective

Rebuild RS-LMTO-ASA LRTDDFT from published formulations while leaving no physics discretion to the coding agent.

A blocked result is acceptable and preferred to an undocumented approximation.

## Fixed contracts

### C1 — TDDFT response
`[LITERATURE]`
The enhanced susceptibility is constructed from `chi_KS` and an XC interaction through the published Dyson equation. No previous `tddft_chi0` convention is grandfathered in.

### C2 — collinear ALSDA transverse kernel
`[LITERATURE]`
Buczek–Ernst–Sandratskii 2011 Eq. (22), in their convention:
\[
K_{\rm xc}(x)=-\mu_B B_{\rm xc}(x)/m(x).
\]
RS-LMTO sign, `mu_B`, density definition and energy units must be derived before implementation.

### C3 — spatial response basis
`[LITERATURE]`
The 2011 formulation uses a localized site × spherical-harmonic × radial basis; in ASA the cells are atomic spheres. The 2026 formulation uses direct radial mesh and super-index
\[
\lambda=(i,s,L,\rho).
\]
No site-only substitute is permitted.

### C4 — Goldstone identity
`[LITERATURE]`
For SOC-off collinear response:
\[
D=I-\chi_{\rm KS}K_{\rm xc}
\]
has the rigid-rotation null space at `q=0`, `omega=0`.

### C5 — published eigenvalue correction
`[LITERATURE]`
Buczek et al. 2011 Eqs. (38–39): diagonalize `D`, set only the small Goldstone eigenvalue to zero, reconstruct `D_corr`, obtain a corrected kernel, and use that kernel for all q and omega. This stays separate from the Lounis sum-rule route.

### C6 — published ASA sum rule
`[LITERATURE]`
Lounis et al. derive in ASA:
\[
\sum_j\int dr'\sum_{LL_1}\chi^{iLL_1;jL_1L}_0(r,r';0)B^j_{\rm eff}(r')=4\pi m_z^i(r),
\]
and in their convention
\[
U^j(r)=B^j_{\rm eff}(r)/(4\pi m_z^j(r)).
\]
The RS-LMTO mapping of `4pi`, field, density and radial normalization must be exact.

### C7 — 2026 non-collinear extension
`[LITERATURE]`
Four density channels `(0,x,y,z)`, local-frame ALSDA kernel, rotation to global frame, site/angular/radial response basis, and rigid-rotation Goldstone null space. No transverse-only shortcut.

## Hard blockers

Stop if any of these cannot be proved:
1. radial `m(r)` identity and normalization;
2. radial `Bxc(r)` identity and normalization;
3. identical XC provenance for ground state and kernel;
4. exact mapping of LMTO states/GFs to radial/angular response;
5. exact finite-q Fourier convention;
6. exact spin/circular factors;
7. exact Ry/Ha/muB normalization;
8. exact mapping of the published ASA sum rule;
9. required radial data survives or can be recovered exactly from the converged state.

## Clean-room rule

After LR-00:
- a compiling build with no working TD-DFT is acceptable;
- old TD-DFT inputs must fail explicitly;
- old pair-Xi/rotation-tangent physics is not a reference.

Generic eigensolvers, Green functions, radial SCF data, matrix algebra and Fourier utilities may remain only when documented as physics-neutral.

## Luna authority

Luna may choose technical representation details such as private helpers, array ordering, loop structure and error plumbing.

Luna may **not** choose equations, signs, factors, projections, angular truncation, radial approximation, Goldstone method, q convention or whether two published schemes may be hybridized.

## Evidence after every task

Report:
1. starting commit and clean-tree status;
2. literature contract IDs;
3. source provenance;
4. equations/units;
5. tests and commands;
6. numerical evidence;
7. blockers;
8. checklist;
9. commit only when justified.
