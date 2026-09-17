# NC-00 — Feasibility audit for a noncollinear four-component response

## Prerequisite

TDVAL-01 collinear transverse response validated.

## Nature

BLOCKER AUDIT ONLY.

No noncollinear susceptibility implementation unless a later prompt explicitly
authorizes it.

## Goal

Assess whether the validated response infrastructure can support the modern
four-component noncollinear formulation without reducing it to a transverse-only
shortcut.

Target channels:

\[
(0,x,y,z).
\]

## Important new baseline fact

The validated collinear route is a **Pauli/no-SOC response built on a
scalar-relativistic ground state**.

Do not silently assume that this immediately extends to:

- noncollinear radial augmentation;
- SOC;
- full scalar-relativistic spin-angular response.

Audit those seams explicitly.

## Questions

1. Can production eigenstates provide all Pauli spinor components in one global
   frame for arbitrary noncollinear ground states?
2. Can the LMTO radial augmentation apply the correct two local radial spin
   channels to a general spinor and rotate back without ambiguity?
3. Does the accepted ground-state snapshot provide:
   \[
   n(r),\quad \mathbf m(r)
   \]
   or only a local-axis pair of radial spin densities?
4. Can all 16 bare response channel pairs be represented with LR-04/LR-05?
5. What is the exact local-frame ALSDA 4×4 kernel for the selected XC?
6. Is the published neglect of intra-atomic noncollinearity compatible with the
   desired scope?
7. How is the rigid-rotation Goldstone null space represented for a general
   magnetic texture?
8. What changes when SOC is enabled?
9. Does the unresolved exact scalar-relativistic density operator from LR-02R
   become more serious in the noncollinear case?

## Rule

Any missing channel normalization, local-frame quantity, or spinor/radial mapping
that would require inventing an unproved formula => BLOCKED.

No site-only workaround.

No “rotate the collinear result” shortcut unless that is exactly the published
approximation and all assumptions are documented.

## Deliverable

`docs/TDDFT_NONCOLLINEAR_FEASIBILITY.md`

Provide a capability matrix:

- Pauli/no-SOC noncollinear;
- scalar-relativistic noncollinear;
- SOC/full relativistic;
- charge-spin coupling.

## Commit

`docs: assess noncollinear TDDFT response feasibility`
