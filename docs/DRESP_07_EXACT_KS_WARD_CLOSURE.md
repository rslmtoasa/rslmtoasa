# DRESP-07 exact Kohn–Sham Ward closure

Status: **PASS — SECOND_ORDER_LMTO_MAPPING_REQUIRED localized**.

This DRESP-07R repair supersedes the previous “best local radial field”
diagnostic. The old value `0.844` was produced by treating nonorthogonal
radial rows independently; it was not a true least-squares optimum.

> **SUPERSEDED DIAGNOSTIC:** old implementation treated nonorthogonal radial
> rows independently and was not a true least-squares optimum.

The bounded campaign used accepted bcc Fe on a 4 × 4 × 4 full reciprocal mesh,
Gamma, static frequency, and the eta ladder 0.04, 0.02, 0.01, 0.005 Ry. The
accepted reciprocal eigensystem and `hk_bulk` were consumed as-is. No BES/Halle
correction, Goldstone repair, Kxc tuning, kernel rescaling, fitted interaction,
or eigenvalue pinning was used.

The diagnostic follows the required chain:

```text
H -> delta H = -i [G,H] -> delta rho_rot
  -> delta rho_spec -> exact-H product response
  -> B_H -> best local spherical operator -> radial SVD inverse
  -> B_KS^LMTO -> B_xc^LMTO
```

The complete Fe/spd product representation is initialized and retained at
dimension 232. The Ward observable is the L=0 transverse spin response inside
that representation; higher-L coordinates are left present in the 232-space
but are not needlessly contracted for this scalar observable.

## Implementation

`source/lr_exact_ks_ward.f90` provides independent finite-dimensional oracles:

- site-major `G = sigma_y/2`, exact `delta H`, and `rho=f(H)`;
- finite-temperature divided differences, including the degenerate limit
  `K_nn=f'(epsilon_n)`;
- accepted `H C = C epsilon` and Hermiticity audits;
- independent density rotation/spectral responses;
- exact-H and radial-field product contractions.

`source/lr_dresp07_bridge.f90` owns the material diagnostic and writes the
artifact selected by `output_file`. The production driver accepts
`backend='exact_ks_ward'` only for one Gamma point, one static frequency, the
accepted k-space-SCF handoff, and the complete response cutoff.

The circular convention is algebraically fixed, not fitted. With
`G=sigma_y/2` and collinear `H`,

```text
delta H = -i [G,H] = B_H sigma_x
B_H = (H_up - H_down)/2
```

The live `chi_plus` contraction therefore drives with `B_H`, not `2 B_H`; the
existing explicit factor two is retained for the real transverse response.
The independent finite-H regression verifies the sign and factor two.

The radial field comparison is explicitly a diagnostic mapping: it integrates
the accepted Pauli large-component basis with the radial scalar field using the
accepted logarithmic mesh and Simpson metric. It is not asserted to be an
identity with the native second-order reciprocal Hamiltonian.

`source/lr_dresp07_radial_oracle.f90` keeps the two representation questions
separate:

- `best_local_spherical_operator` is the exact Frobenius projection onto the
  mutually orthogonal site/l projectors, with one m-degenerate coefficient per
  site and angular momentum.
- The accepted radial map is assembled as one coupled site/l-by-site/r matrix
  and inverted with DGELSS/SVD. The radial profile is therefore minimum-norm;
  no independent per-l division is used.

## Quantitative residual ledger

Artifact: `/tmp/dresp07_fe_4k.dat` (deliberately not committed).

| Arrow / diagnostic | Measured value | Result |
|---|---:|---|
| accepted `H C = C epsilon`, max element | 8.6062e-13 | pass |
| accepted eigenpair residual, Frobenius | 2.1289e-12 | pass |
| accepted H Hermiticity residual | 1.2637e-12 | pass |
| accepted spin-offdiagonal norm, all k | 0 | pass |
| `H -> delta H` rotation-field relation | 0 | pass |
| `delta rho_rot -> delta rho_spec` relative | 3.1676e-14 | pass |
| exact rotation product vs spectral product | 7.1976e-15 | pass |
| exact `delta H` product response vs spectral response | 1.9259e-16 | pass |
| `m_H^prod` vs accepted Pauli valence | 4.0699e-02 | finite projection mismatch |
| exact-H response vs accepted total Pauli | 5.3718e-02 | core-inclusive comparison |
| old best local radial `B_H` relative residual | 8.4382e-01 | **SUPERSEDED DIAGNOSTIC** |
| best local spherical operator relative residual | 9.0148e-02 | variational optimum |
| best local spherical operator action residual | 1.0511e-01 | variational optimum |
| radial map SVD rank / condition | 3 / 1.1445e+1 | full-rank solve |
| radial map singular values | 2.9487e-1, 9.8699e-2, 2.5764e-2 | full spectrum |
| radial map relative residual | 5.1355e-16 | realizable coefficients |
| radial best `B_H` matrix residual | 9.0148e-02 | agrees with operator optimum |
| mapped `B_xc` vs `B_H`, matrix relative | 2.0125e-01 | mismatch |
| mapped `B_xc` vs `B_H`, action relative | 5.5047e-01 | mismatch |
| mapped total radial `B_KS` vs `B_H`, matrix relative | 2.0125e-01 | mismatch |
| mapped total radial `B_KS` vs `B_H`, action relative | 5.5047e-01 | mismatch |

The accepted Gamma field has `||B_H||_F = 2.2879e-1` in the full spin
embedding. Its onsite norm is 2.2879e-1 and its intersite/site-offdiagonal
norm is zero for this one-site bcc primitive-cell representation. The full
orbital field still contains 9.4141e-3 of same-l orbital offdiagonal weight and
1.1139e-2 of within-l diagonal m-anisotropy. The cross-l norm is approximately
1.0e-16 and the intersite norm is zero. Across the accepted mesh, the orbital
`B_H` Frobenius norm spans 1.6178e-1 to 1.9938e-1.

The true local spherical optimum leaves a 9.0148e-2 Gamma relative residual,
while its coefficients are reproduced by the radial inverse to 5.1355e-16.
Thus operator-space representability, not radial realizability, is the
remaining Gamma defect. The best local operator satisfies the hard variational
oracle against both mapped `B_xc` and mapped `B_KS` fields.

The weighted field mean has norm 1.7420e-1, with mesh k-dependence
`R_k = 2.6749e-1`. The global weighted local projection residual is only
2.7795e-2, while the per-k local residual averages 2.6058e-1 and reaches
3.3318e-1. This motivates the native second-order mapping milestone; it is
not evidence that a physical local field itself is invalid.

## Valence/core bookkeeping

The accepted product-space norms are:

| Quantity | Norm |
|---|---:|
| total Pauli magnetization | 6.98494e-1 |
| valence Pauli magnetization | 6.90750e-1 |
| frozen-core Pauli magnetization | 1.07143e-2 |
| core / total | 1.53391e-2 |

The response uses reciprocal valence eigenstates. The exact-H response residual
is 4.0699e-2 against valence and 5.3718e-2 against total magnetization. Thus
the frozen core is a measured bookkeeping contribution, but it is not the
dominant localization: the exact-H response algebra closes independently and
the local radial representation remains the much larger defect.

## Static eta ladder

The divided-difference result is the eta=0 finite-dimensional reference. The
retarded rows below are relative to their corresponding static source response.

| eta (Ry) | exact H | best local spherical H | total radial KS | radial XC |
|---:|---:|---:|---:|---:|
| 0.040 | 2.4738e-1 | 2.5409e-1 | 2.4938e-1 | 2.4938e-1 |
| 0.020 | 1.2772e-1 | 1.3237e-1 | 1.2882e-1 | 1.2882e-1 |
| 0.010 | 6.4549e-2 | 6.7274e-2 | 6.5125e-2 | 6.5125e-2 |
| 0.005 | 3.2375e-2 | 3.3808e-2 | 3.2668e-2 | 3.2668e-2 |

The exact-H retarded response approaches the exact static oracle as eta falls;
this is not a static divided-difference failure.

## Classification and disposition

Primary classification:

```text
SECOND_ORDER_LMTO_MAPPING_REQUIRED
```

The accepted-H eigensystem, global-spin-rotation identity, circular convention,
and exact-H product contraction all pass at numerical precision. The old 84.4%
number is a superseded non-variational fit. The true best local spherical
operator leaves 9.01% at Gamma, while the radial inverse reproduces its shell
coefficients at numerical precision and passes the variational inequality
against both mapped fields. The dominant remaining measured seam is the
26.75% k-dependence of the accepted exact coefficient-space field, which a
k-independent direct radial matrix cannot reproduce. The 20.1% mapped radial
XC/KS matrix mismatch and 55.0% action mismatch are retained as mapping
diagnostics. The 4.07% exact-H versus Pauli-valence projection mismatch remains
a separate finite product/augmentation diagnostic, not hidden by a fit.

BES/Halle status: **NOT APPLIED**.

Is BES/Halle justified as the next step? **NO**. The exact finite-H chain does
not require a Goldstone repair; the unresolved seam is field representation and
radial-to-Hamiltonian mapping. A BES/Halle correction would conceal that seam.

Recommended next milestone: define and validate a native second-order mapping
from the complete accepted `B_H(k)` orbital/site field into the response
operator representation, including the non-scalar orbital structure, then rerun
the same valence-target Ward ladder before considering any kernel or Goldstone
work.

## Tests and frozen paths

Focused regressions:

- `UnitDresp07ExactKsWard`: coupled finite-H commutator/divided-difference
  identity, finite-temperature limit, circular sign/factor two, accepted H
  extraction, and exact-H-to-product contraction.
- `UnitDresp07RadialOracle`: independent nonorthogonal-row SVD solve,
  nonrepresentable negative fixture, exact projector-space least-squares
  fixture, and hard best-operator variational inequality.
- `TddftDresp07ExactKsWard`: bounded accepted-state bcc Fe 4 × 4 × 4 campaign.
- `Dresp07FeArtifact`: independent artifact checks for exact gates,
  valence/core bookkeeping, variational/radial diagnostics, decomposition,
  k-dependence, and eta ladder.

The DRESP-06A projection and product-response tests remain unchanged. Frozen
paths were not modified:

```text
source/exchange.f90
source/lr_lmto_turek_gf.f90
source/lr_lmto_turek_contour.f90
```

No DRESP-03TG physics was changed.

Suggested commit:

```text
Repair DRESP-07 radial field variational oracle
```
