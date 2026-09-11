# KXC-01 — Implement the local radial ALSDA transverse XC kernel

## Prerequisites

- LR-03 PASS.
- LR-04 PASS.
- LR-06 PASS.
- LR-03 must specify a scientifically defensible SR/Pauli field-density pairing
  for the initial kernel. If it declares this point BLOCKED, STOP.

## Goal

Implement the direct local transverse ALSDA kernel in the exact conventions and
response-space algebra already certified.

No Goldstone correction.
No sum-rule-derived interaction.
No empirical rescaling.

## Ground-state source

Use only the accepted LR-01 ground-state snapshot and LR-03 canonical dictionary.

Primitive quantities include the exact converged:

\[
V_{xc,\uparrow}(r),\qquad
V_{xc,\downarrow}(r),
\]

and the magnetization quantity selected by LR-03 for this response approximation.

Do not reconstruct the field from LMTO potential parameters.

## XC provenance

The response kernel must carry and verify the XC provenance of the accepted
ground state.

If the requested response functional does not match the ground-state functional,
fail loudly.

Do not silently substitute a convenient LSDA when the ground state used another
functional.

## Physics

Use the LR-03 mapped form of the published transverse ALSDA relation.

Do not re-derive signs or factors locally.

Keep separate:

- `delta_vxc`;
- the energy-valued Pauli coefficient;
- any physical magnetic field;
- number-spin density;
- magnetic-moment density.

## Local operator representation

Use LR-04.

A local radial kernel is a pointwise multiplication operator in the continuum.

Its stored matrix/action must follow the LR-04 discrete metric exactly.

Do not guess where inverse quadrature weights or radial `r^2` factors belong.

For a spherical ground state, prove the angular diagonal structure rather than
constructing a dense angular matrix unnecessarily.

## Near-zero magnetization

Inspect all radial points carrying non-negligible response measure.

If the mapped expression contains a ratio such as

\[
B_{xc}(r)/m(r),
\]

do not introduce:

- epsilon floors;
- clipping;
- spline replacement;
- arbitrary radial cutoffs.

If the formal functional derivative provides a mathematically equivalent finite
limit, it may be used only after proving equivalence for the same selected XC
functional and documenting it.

Otherwise declare the affected capability BLOCKED.

## Tests

1. pointwise reconstruction from LR-01 quantities;
2. units from LR-03;
3. local-operator action using LR-04 metric;
4. angular diagonal structure;
5. global spin reversal;
6. XC-provenance mismatch fails;
7. low-m diagnostic;
8. direct quadrature versus stored operator application.

## Raw static diagnostic

Form the LR-03 static Goldstone/Ward residual using the unmodified LR-06
`chiKS` and the direct ALSDA kernel.

Report:

- residual norm;
- rigid-rotation overlap;
- radial/angular distribution of the residual where useful.

Do not modify the kernel.

## Deliverable

`docs/KXC_ALSDA_TRANSVERSE_KERNEL.md`

Clearly distinguish:

- functional identity;
- discretization;
- raw numerical Goldstone error.

## Implementation outcome

`source/lr_alsda_kernel.f90` implements the LR-03 controlled mixed
`K_xc^(P<-SR)=B_xc^sigma,SR/s^P` relation.  It consumes accepted LR-01 radial
snapshots, verifies XC provenance, requires an explicitly labelled Pauli
magnetization, and delegates the canonical local operator to LR-04.  The
origin is handled as LR-04's exact null-measure extension; active-grid zeros
are rejected and nonzero low-m values are evaluated without regularization.

`UnitLrAlsdaKernel` passes its pointwise, operator-action, diagonal-structure,
spin-reversal, low-m, and raw-static-diagnostic checks.  The separate
`UnitLrAlsdaKernelRejectProvenance` expected-failure test confirms that a
requested functional mismatch fails loudly.

## Acceptance checklist

- [x] accepted LR-01 `Vxc_up/down` arrays are the sole XC field source;
- [x] LR-03 `B_xc^sigma=(Vxc_up-Vxc_down)/2` and `K_xc=B_xc^sigma/s^P` are used;
- [x] Pauli response magnetization is explicit and cannot fall back to SR density;
- [x] XC provenance is carried and mismatches fail loudly;
- [x] LR-04 canonical local operator and metric action are used;
- [x] spherical angular diagonal structure is represented analytically;
- [x] active-grid zero magnetization blocks without a floor or cutoff;
- [x] origin null-measure behavior is explicit and unregularized;
- [x] pointwise, units, action, angular, spin-reversal, provenance, and low-m tests pass;
- [x] raw static Goldstone/Ward residual is reported without modifying `chiKS` or `Kxc`;
- [x] no Goldstone correction, sum-rule interaction, empirical rescaling, or Dyson physics added;
- [x] implementation documentation created at `docs/KXC_ALSDA_TRANSVERSE_KERNEL.md`.

## PASS condition

PASS requires the literature-mapped kernel to be represented exactly in the
certified response algebra with no empirical adjustment and with controlled
low-m behavior.

A nonzero raw Goldstone residual does not by itself fail KXC-01; it is evidence
for later GSR/GCR analysis.

## Commit

`td-dft: implement radial ALSDA transverse kernel`
