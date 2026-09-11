# LR-GF-02 — Independent reciprocal-GF construction of the bare susceptibility

## Timing

Run after LR-06.

Strongly recommended before final collinear validation.

It is not allowed to replace LR-06 with the same spectral implementation hidden
behind a GF wrapper.

## Goal

Construct the same bare transverse susceptibility from the certified reciprocal
one-electron Green function and compare it against the LR-06 explicit spectral
band sum.

This is an **independent representation-level cross-check**.

## Foundation

LR-GF-01 established, for the supported `ham_only` mode,

\[
G_c(\mathbf k,z)=[zI-H(\mathbf k)]^{-1}
=
\sum_n
\frac{c_{n\mathbf k}c^\dagger_{n\mathbf k}}
{z-\epsilon_{n\mathbf k}}.
\]

LR-BASIS/LR-05 establish the Pauli radial/angular vertex.

Use those contracts.

## Literature derivation first

Before coding, derive the retarded GF expression for the selected circular
susceptibility from the same Kubo function used in LR-03.

Derive the placement of:

- retarded GF;
- advanced GF;
- spectral/imaginary part;
- Fermi function;
- `E+omega`;
- endpoint ordering;
- complex conjugation.

Do not copy a historical GF formula without reconciling it with LR-03.

Document equivalence to the spectral Lehmann sum.

## Numerical integration

Choose and document one energy integration strategy for the first implementation.

Possible choices include:

- real-axis adaptive/fixed quadrature;
- contour plus real-axis correction if rigorously derived.

Do not introduce contour machinery merely because another KKR code uses it.

The initial purpose is correctness and cross-validation.

## Independence requirement

The GF susceptibility backend may reuse:

- one-electron `G(k,z)`;
- response-space mapping;
- radial/angular operator vertices.

It must **not** call the LR-06 band-pair susceptibility accumulator.

Its energy integral must be genuinely evaluated as a GF expression.

## Tests

### Finite spectrum

Use a tiny finite Hamiltonian for which:

1. LR-06 spectral sum;
2. GF energy integral;
3. direct analytic result

can all be compared.

Vary eta and integration resolution.

### Reciprocal production fixture

For selected q, omega points compare GF and spectral `chiKS` in the full response
space.

Report norm differences, not only a site-projected scalar.

### Static limit

Compare the q=0, omega=0 result independently.

### Sum-rule/integrated spectral check

Where LR-03 provides an exact frequency integral identity, use it as a diagnostic.

Do not enforce it.

## Error budget

Separate:

- energy-integration discretization;
- eta;
- k mesh;
- electronic finite basis;
- SR→Pauli approximation inherited from the common vertex.

The spectral and GF routes share the last two and therefore do not independently
validate them.

State that explicitly.

## Deliverable

`docs/LR_GF_SUSCEPTIBILITY_CROSSCHECK.md`

## PASS condition

PASS means the independently evaluated reciprocal-GF susceptibility converges to
the LR-06 spectral result for the same finite electronic problem and response
basis.

It does not validate ALSDA or magnons.

## Forbidden

- no Kxc;
- no Dyson;
- no site-only comparison as sole oracle;
- no wrapper around spectral chiKS.

## Commit

`linear-response: cross-check KS susceptibility with reciprocal GF`

## Completion checklist

- [x] derivation reconciled with LR-03 before coding;
- [x] independent GF energy integral implemented;
- [x] finite-spectrum and full response-space tests pass;
- [x] deliverable evidence document added;
- [x] commit prepared with the prescribed message.
