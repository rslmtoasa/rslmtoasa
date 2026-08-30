# Luna Prompt — GBT-WP05: Central Rotating-Frame / Lab-Frame Contract

## Mission

Make the GBT frame semantics explicit, centralized, and testable across Hamiltonian construction, density/moments, output, explicit-supercell comparison, and future constraining-field use.

The goal is to eliminate implicit assumptions such as “`potential%mom` happens to be in the right frame.”

## Dependency

WP04 should pass so the operator itself is no longer the primary suspect.

## Core convention to establish

Adopt and document one authoritative convention for `gbt_single_q`.

Recommended convention:

- local LMTO potential parameters live in the primitive **rotating frame**;
- the prescribed reference moment axis for each sublattice is cell-independent in that frame;
- the lab-frame moment in translated cell n is reconstructed by the single-q rotation;
- transverse density components in a constrained spiral are residual/torque information unless the selected density policy explicitly allows the reference axis to relax.

Do not silently let one routine reinterpret the same 3-vector in another frame.

## Recommended code structure

Introduce or consolidate a small authoritative helper/module, for example conceptually:

```fortran
type :: gbt_frame_t
   complex(rp) :: U(2,2)
   real(rp)    :: R(3,3)
end type
```

with operations equivalent to:

- construct frame for sublattice/cell;
- endpoint link;
- rotating -> lab vector;
- lab -> rotating vector;
- rotating -> lab spinor/density matrix;
- lab -> rotating spinor/density matrix.

Use existing helpers where they are already correct; do not create gratuitous duplication.

## Critical audit: `potential%mom`

Trace all GBT uses of `potential%mom`.

If the Hamiltonian builder only uses the sign of its z component to choose local up/down channels, make that a declared invariant rather than silently ignoring x/y.

For strict constrained/reference-frame GBT, add a validation/assertion that transverse components are within tolerance before the GBT Hamiltonian is built, unless a clearly supported relaxed-reference mode is active.

The error message must explain where physical cone direction belongs (`theta_ss`, sublattice reference angles, or equivalent), rather than merely saying “invalid input.”

## Density covariance tests

Construct known rotating-frame density matrices

\[
\bar\rho_a=\frac12(n_aI+\mathbf m_a\cdot\sigma)
\]

and verify that translation to the lab frame yields:

- invariant charge;
- invariant moment magnitude;
- fixed cone angle;
- azimuthal advance exactly `q dot R`;
- correct sublattice phase offset.

For commensurate q, compare the reconstructed lab-frame moments against the explicit supercell from WP04 site by site.

## Reciprocal density producer

Audit the actual reciprocal-space eigensolver/density code:

- In what frame are the computed `(mx,my,mz)` components?
- Are they mixed directly into `potential%mom`?
- Is the radial-potential update supposed to use the full vector or only the longitudinal component in constrained-spiral mode?

Do not repair the full SCF yet. Define and test the contract first.

## Deliverables

- centralized frame helper or explicit consolidation of existing helpers;
- frame covariance tests;
- GBT moment invariant/assertion where appropriate;
- `docs/dev/GBT_REAUDIT_WP05_FRAME_CONTRACT.md`.

## Completion checklist

- [ ] One rotating/lab frame convention is documented.
- [ ] SU(2) spinor and SO(3) vector rotations are mutually tested.
- [ ] Density matrix transforms preserve charge and |m|.
- [ ] Spiral azimuth/cone reconstruction is tested.
- [ ] Multi-sublattice phase offsets are tested.
- [ ] Explicit-supercell moments agree site by site in a commensurate case.
- [ ] `potential%mom` semantics are explicit.
- [ ] Invalid transverse use of `potential%mom` is rejected or deliberately supported.
- [ ] Reciprocal density producer frame is identified.
- [ ] No broad SCF tuning is performed yet.

## Commit

`refactor: define explicit GBT frame contract`
