# NC-00 — Blocker audit for the 2026 non-collinear four-component extension

## Prerequisite
Collinear transverse TDDFT validated.

## Nature
Feasibility/blocker audit only.

## Goal
Assess whether the validated RS-LMTO radial response basis can support the 2026 four-component non-collinear formulation without shortcuts.

## Published target
EEB-01 through EEB-06:
- channels `(0,x,y,z)`;
- site/angular/radial basis;
- local-frame ALSDA kernel;
- global-frame transformation;
- non-collinear KS response;
- rigid-rotation Goldstone null space.

## Audit questions
1. Can chiKS produce all 16 channel pairs?
2. Is charge response normalized consistently with spin response?
3. Are local non-collinear n(r) and vector m(r) available?
4. Is the local ALSDA 4x4 kernel derivable from the selected XC functional?
5. Can local-to-global rotation be applied exactly as published?
6. Is neglect of intra-atomic non-collinearity an explicit approximation matching the intended published scope?
7. Does current LMTO non-collinear/SOC machinery preserve the response basis?
8. Can the required rigid rotations be represented and tested?

## Rule
Any missing density channel, normalization or local-frame quantity => BLOCKED.

No transverse-only workaround.

## Deliverable
`docs/TDDFT_NONCOLLINEAR_2026_FEASIBILITY.md`

## Checklist
- [ ] 16 channels assessed
- [ ] charge normalization mapped
- [ ] vector radial magnetization assessed
- [ ] 4x4 ALSDA kernel mapping assessed
- [ ] local/global rotation assessed
- [ ] intra-atomic NC assumption documented
- [ ] rigid-rotation null space mapped
- [ ] PASS or BLOCKED

## Commit if documentation is added
`docs: assess four-component noncollinear TDDFT feasibility`
