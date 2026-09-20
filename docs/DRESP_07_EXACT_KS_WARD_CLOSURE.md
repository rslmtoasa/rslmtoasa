# DRESP-07 exact Kohn–Sham Ward closure

Status: **PASS — FIELD_REPRESENTATION_FAILURE localized**.

The bounded campaign used accepted bcc Fe on a 4 × 4 × 4 full reciprocal mesh,
Gamma, static frequency, and the eta ladder 0.04, 0.02, 0.01, 0.005 Ry. The
accepted reciprocal eigensystem and `hk_bulk` were consumed as-is. No BES/Halle
correction, Goldstone repair, Kxc tuning, kernel rescaling, fitted interaction,
or eigenvalue pinning was used.

The diagnostic follows the required chain:

```text
H -> delta H = -i [G,H] -> delta rho_rot
  -> delta rho_spec -> exact-H product response
  -> B_H -> best local radial B_H -> B_KS^LMTO -> B_xc^LMTO
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
| best local radial `B_H` relative residual | 8.4382e-01 | representation defect |
| mapped `B_xc` vs `B_H`, matrix relative | 2.0125e-01 | mismatch |
| mapped `B_xc` vs `B_H`, action relative | 5.5047e-01 | mismatch |
| mapped total radial `B_KS` vs `B_H`, matrix relative | 2.0125e-01 | mismatch |
| mapped total radial `B_KS` vs `B_H`, action relative | 5.5047e-01 | mismatch |

The accepted Gamma field has `||B_H||_F = 2.2879e-1` in the full spin
embedding. Its onsite norm is 2.2879e-1 and its intersite/site-offdiagonal
norm is zero for this one-site bcc primitive-cell representation. The full
orbital field still contains 9.4141e-3 of orbital offdiagonal weight. The
same-`l` norm is 1.6178e-1 and the cross-`l` norm is approximately 1.0e-16.
Across the accepted mesh, the orbital `B_H` Frobenius norm spans
1.6178e-1 to 1.9938e-1.

The best local radial fit has blockwise radial-fit rank 3 and condition 1.9311;
it captures 2.5198e-1 Frobenius norm and leaves 1.3651e-1 residual. It is
reported only as a representation diagnostic and is never used as a corrected
field. Its retained product operator differs from the mapped radial XC/KS
operator by relative 1.1227.

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

| eta (Ry) | exact H | best H | total radial KS | radial XC |
|---:|---:|---:|---:|---:|
| 0.040 | 2.4738e-1 | 2.7733e-1 | 2.4938e-1 | 2.4938e-1 |
| 0.020 | 1.2772e-1 | 1.4925e-1 | 1.2882e-1 | 1.2882e-1 |
| 0.010 | 6.4549e-2 | 7.7299e-2 | 6.5125e-2 | 6.5125e-2 |
| 0.005 | 3.2375e-2 | 3.9097e-2 | 3.2668e-2 | 3.2668e-2 |

The exact-H retarded response approaches the exact static oracle as eta falls;
this is not a static divided-difference failure.

## Classification and disposition

Primary classification:

```text
FIELD_REPRESENTATION_FAILURE
```

The accepted-H eigensystem, global-spin-rotation identity, circular convention,
and exact-H product contraction all pass at numerical precision. The remaining
defect is exposed when the complete accepted Hamiltonian field is compared to
the scalar local radial/product field space: even the best representable local
radial field leaves an 84.4% Frobenius residual, and the certified radial XC
and total-KS maps differ from `B_H` at the 20.1% matrix level and 55.0% action
level. The 4.07% exact-H versus Pauli-valence projection mismatch is retained
as a separate finite product/augmentation diagnostic, not hidden by a fit.

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
- `TddftDresp07ExactKsWard`: bounded accepted-state bcc Fe 4 × 4 × 4 campaign.
- `Dresp07FeArtifact`: independent artifact checks for exact gates,
  valence/core bookkeeping, field-representation defect, and eta ladder.

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
Diagnose exact KS Ward closure in DRESP-07
```
