# Mn2F2: frozen-potential screened-constant reciprocity intervention

Date: 2026-10-21. Stacked on `codex/mn2f2-gf-symmetry-audit` at
`2ef44c63fc64cc23cc3a914866627f8588f299f4` (PR #10).

## Status and scope

All four full 24-pair axial endpoints completed with the same executable and
unchanged frozen potentials. Restoring only the real-basis Sbar reciprocity
restores Hermitian Hamiltonian ingredients, suppresses forbidden G01/Gx1
orbital blocks, and removes the dominant time-odd longitudinal Dx. The numerical
tables below distinguish surviving residuals from exact symmetry restoration.
The [110]/[-1-10] extension below also completed both full intervention runs.
No incomplete, interrupted or file-replacement-only run is endpoint evidence.

This is a fixed-potential intervention, not an SCF calculation or a study of
screening convergence. Screening cutoffs, real-space sums, HoH, SOC, recursion,
energy grid, atom files and exchange formulas are unchanged. Neither Dx nor any
GF component is set to zero. Both diagonal domains were subsequently run in
the extension below; lambda interpolation has not been run.

## Intervention and binary validation

In the real Cartesian orbital basis stored in the binary, each directed block
is matched by source type, destination type and displacement to its unique reverse:

$$X_{ij}(R)=\frac{\bar S_{ij}(R)+\bar S_{ji}(-R)^T}{2},\qquad
\bar S^{\rm rec}_{ij}(R)=X_{ij}(R),\quad
\bar S^{\rm rec}_{ji}(-R)=X_{ij}(R)^T.$$

Only a genuinely onsite self-reversed key is averaged with its own transpose.
There is no spherical-orbital C operation and no offsite self-symmetrization.
The original binary is read-only; `sbar.reciprocal` is a separate exclusive-create
file. Four-byte little-endian record markers, nine float64 entries per record,
block order and counts are retained. Rounded preview values are never used as
matrix data. Vector labels identify physical bonds with absolute tolerance
`2e-6` in their Cartesian units; matching must be unique.

| Validation over all 86 directed blocks | Original | Reciprocal |
|---|---:|---:|
| max element of Sij-Sji^T | 2.1474105637e-5 | 0 exactly |
| concatenated Frobenius norm | 1.9715686202e-4 | 0 exactly |
| missing/ambiguous reverse partners | 0 | 0 |

There are 41 distinct offsite reverse pairs and four self-reversed onsite keys
(45 unique pairs including onsite). Maximum `|delta S|=1.0737052818e-5` and
concatenated `||delta S||F=9.8578431010e-5`. Parse/encode roundtrip preserves all
bytes; applying the projection twice gives identical bytes. The two axial jobs
give the same projected binary hash.

Original SHA256:
`951f587460c1ade64d6676e1a6ab591f3bd7de831fe5ff54c716cfa25f818c4a`.
Reciprocal SHA256:
`ab18c763931f6f794b06abd23a3c4b96c959e46c81a3bc8a47a21fe30cc509eb`.

## Why an explicit diagnostic input hook is necessary

The actual exchange entry point calls `structb(.true.)`, which recomputes Sbar
in memory. It does **not** read saved `sbar.`. Replacing that file alone left the
captured Hamiltonian unchanged in an initial negative control. Consequently,
such a run cannot test the proposed intervention.

The independently regenerated and archived original Sbar agree to a maximum
element difference of `2.2204460493e-14`, much smaller than the measured
reciprocity defect and the projected correction.

The optional `RSLMTO_AUDIT_SBAR` hook loads the projected copy after the unchanged
screening call and before Hamiltonian assembly. It requires GF audit mode and
`RSLMTO_AUDIT_SBAR_VIEW`, checks freshly generated versus saved vector/type/order
labels, checks record counts/finite entries, and announces the load. An unset
variable leaves the ordinary path intact. The runner now rejects a successful
process that did not announce the override. All copied file hashes are checked
before execution: only `sbar.` differs in the intervention; full-run pair lists
are not edited. Saved output-potential hashes are also checked after execution.

## Hamiltonian gate

The input intervention has a direct operator interpretation:

$$h_{ij}-h_{ji}^{\dagger}
=W_i[1_2\otimes(\bar S_{ij}^{T}-\bar S_{ji})]W_j.$$

The actual full finite neighbour graph contains 34,400 directed blocks. Its
read-only check compares each bare-hopping edge with its unique reverse and
checks the weighted hopping used in recursion against `hO=h*O`. It excludes
allocation padding beyond the actual 1600 atoms.

| Full-graph check, [100] | Original | Reciprocal |
|---|---:|---:|
| max bare-h Hermiticity defect (Ry) | 4.6260833085e-7 | 2.7755575616e-17 |
| full bare-h Frobenius defect (Ry) | 1.6038935243e-4 | 4.4961696711e-15 |
| max hO identity residual (dimensionless) | 1.1102230246e-16 | 1.1102230246e-16 |
| O / energy-centre / SOC Hermiticity defects | 0 / 0 / 0 | 0 / 0 / 0 |

These exhaustive graph checks establish the Hermitian ingredients of

$$H=h-hOh+E_\nu+H_{\rm SOC}.$$

Actual first-action captures independently check the effective Hamiltonian,
including the two-hop term, without replacing H by its Hermitian part:

| First bond r/alat=(-1/2,+1/2,0), [100] | Original max (Ry) | Reciprocal max (Ry) |
|---|---:|---:|
| bare h pair | 2.8921200027e-8 | 9.1940344227e-17 |
| HoH pair | 5.5173235879e-8 | 8.3266726847e-17 |
| effective H pair | 5.5026866503e-8 | 9.3675074565e-17 |
| effective H onsite, maximum of i/j | 4.2646425377e-8 | 4.1733285233e-17 |
| direct mz*T covariance over captured rows | 1.2035797359e-10 | 1.2035838992e-10 |

For the reversed domain, the reciprocal effective-pair and onsite maxima are
`8.3266726847e-17` and `5.5511151231e-17` Ry. Direct mz*T covariance remains
`1.2045175968e-10` Ry. Bare-h spatial covariance over all 86 type/displacement
blocks stays controlled: the largest tested direct residual remains approximately
`2.34e-10` Ry. Restoring Hermiticity does not erase the saved inputs' small
magnetic/crystal symmetry differences.

The graph check certifies all bare-hopping ingredients; the action dumps measure
all output rows of the selected effective-H columns. No dense 28800-by-28800
effective Hamiltonian was constructed. For an intervention, a full-graph ingredient
failure above `1e-12` aborts before recursion.

## Completed first-pair Green-function and exchange gate

The previous audit derives mirror signs and the reciprocal orbital transpose
`R(A)=C A^T C`. Its sector definitions and all exchange contractions are retained.
Direct antiunitary H covariance was already controlled. Hermiticity is the
additional assumption needed to obtain the same-retarded-energy transpose rule;
the intervention does not impose a mirror projection on any GF.

The measured production contractions are

$$SC_x=\Re\operatorname{Tr}[\Delta_iG^{00}_{ij}\Delta_jG^{x1}_{ji}],\qquad
CC_x=\Re\operatorname{Tr}[\Delta_iG^{01}_{ij}\Delta_jG^{x0}_{ji}].$$

For orbital order `(s; p-1,p0,p1; d-2,d-1,d0,d1,d2)`,
`M=diag(1,1,-1,1,1,-1,1,-1,1)`. The individual matrices obey
`M G^a M = s_a G^a` in the ideal magnetic mirror symmetry, with:

| Sector a | 00 | 01 | x0 | x1 | y0 | y1 | z0 | z1 |
|---|---|---|---|---|---|---|---|---|
| Mirror sign s_a | + | - | + | - | + | - | - | + |
| Nominal bare-T domain sign t_a | + | - | - | + | - | + | - | + |

Delta commutes with orbital mirror M. Ideal G00/Gx0 blocks preserve mirror
parity; ideal G01/Gx1 blocks change it, making both x traces vanish individually.
Neither Gx0 nor Gx1 is required to vanish as a matrix. The original dominant
leaks pair even G00 with odd forbidden Gx1, and even forbidden G01 with odd
allowed Gx0. Thus both erroneous x contractions are time odd. Allowed z factors
have matching mirror parities and can give finite, time-even DMI. The trace
decomposition below uses exactly these measured matrix factors, not a modified
physical exchange expression.

For z, an allowed/forbidden product has zero trace; leakage enters that mirror
trace through two forbidden factors. This explains the different symmetry
selection, but does not assert that every change in Dz is quadratic or that Dz
is immune to the changed Hamiltonian. Both allowed and forbidden GF components
respond to the intervention.

These are maxima of `||M G M - s G||F` in inverse Ry across the 64 sampled
energies, ij direction:

| Sector, [100] first pair | Original | Reciprocal |
|---|---:|---:|
| G01 | 3.7688176659e-6 | 4.7565269240e-10 |
| Gx1 | 3.6717589946e-6 | 5.0783074214e-10 |
| Gx0 control | 6.3844985739e-8 | 1.9802792208e-8 |
| Gz1 control | 3.0725215504e-8 | 1.5907801344e-9 |

For a sector of mirror sign s, the forbidden-block norm is half this residual.
The small remaining leakage should not be called exact symmetry restoration.

Using the **same ordered bond** in reversed domains, in mRy:

| Quantity | Original | Reciprocal |
|---|---:|---:|
| Dx([100]) | +3.886604e-7 | -1.378187e-11 |
| Dx([-100]) | -3.886879e-7 | -1.377696e-11 |
| Dx_odd = (Dx+ - Dx-)/2 | +3.8867415e-7 | -2.455e-15 |
| Dx_even = (Dx+ + Dx-)/2 | -1.375e-11 | -1.3779415e-11 |
| Dz([100]) | -1.700948e-3 | -1.701022e-3 |
| Dz([-100]) | -1.700950e-3 | -1.701023e-3 |

The odd signal falls by about eight orders of magnitude on this bond, while the
small even remainder persists. Dz is finite and changes under the intervention;
it must not be claimed to be unaffected. Printed-output precision limits the last
digits of the parity differences. The full-run analysis below uses the same MPI executable for both input states.

## Complete axial endpoint comparison

All four manifests have return code zero, preserved potential/input hashes and
24 unique GF captures plus 24 H-action captures. Only intervention copies change
`sbar.`; every intervention announces the explicit saved-Sbar load. The common
executable SHA256 is
`3ba06baf0852f56449e4826e263b06e5954a343aa19b5439982c39f29ed890ed`.
The interrupted reversed-domain attempt was discarded; a fresh copy completed.

One successful original [100] run left the pair-9 H diagnostic file truncated,
although its exchange output and all GF captures completed. The selected H
column was recaptured from a fresh one-pair copy with the same executable,
frozen Hamiltonian and four OpenMP threads (one MPI rank instead of two).
All **10,417,047 bytes** retained after the pair-label header match the recapture
byte for byte. A derived completed capture preserves the original `S 9` label;
the truncated original and unedited recapture remain in the private evidence.
The 24-column H tables include this explicitly recovered column. Full-run
exchange/GF data were not replaced. The cause of this isolated diagnostic-file
truncation is unknown; the analysis rejected it before recovery.

Sums below cover the 12 ordered Mn1-to-Mn2 bonds, in mRy. Odd and even refer
to the same ordered bonds in the two reversed domains, not reversed bond order:

$$D_x^{odd}=(D_x^{[100]}-D_x^{[-100]})/2,\qquad
D_x^{even}=(D_x^{[100]}+D_x^{[-100]})/2.$$

| Sum / parity (mRy) | Original | Reciprocal |
|---|---|---|
| Dx [100] | 1.8206015400e-06 | -4.3910000000e-16 |
| Dx [1-00] | -1.8206013400e-06 | -2.0409000000e-15 |
| Dx odd | 1.8206014400e-06 | 8.0090000000e-16 |
| Dx even | 1.0000000005e-13 | -1.2400000000e-15 |
| CCx [100] | 6.0664094000e-07 | -7.0730000000e-16 |
| CCx [1-00] | -6.0664080000e-07 | 5.4050000000e-16 |
| CCx odd | 6.0664087000e-07 | -6.2390000000e-16 |
| CCx even | 6.9999999989e-14 | -8.3400000000e-17 |
| SCx [100] | 1.2139605560e-06 | 2.6500000000e-16 |
| SCx [1-00] | -1.2139605510e-06 | -2.5829000000e-15 |
| SCx odd | 1.2139605535e-06 | 1.4239500000e-15 |
| SCx even | 2.5000001679e-15 | -1.1589500000e-15 |

Separately printed D and CC+SC have finite decimal output precision; differences at that precision are not new terms.

The summed odd signal decreases by more than nine orders of magnitude.
Across all 24 ordered bonds, the largest absolute odd Dx decreases from
`3.8867415e-7` to `9.2455e-15` mRy. Small even bondwise residuals persist.
Both original baselines reproduce their supplied DMI outputs within
`4.13e-14` mRy per component (maximum bound checked after the runs).

### Hamiltonian residuals over full captures

| Domain | Maximum element residual (Ry) | Original | Reciprocal |
|---|---|---|---|
| 100 | bare_pair | 7.1784927442e-08 | 1.8962018370e-16 |
| 100 | hoh_pair | 5.5173235987e-08 | 1.1796130168e-16 |
| 100 | effective_pair | 1.0694632413e-07 | 1.8471416906e-16 |
| 100 | effective_onsite | 4.2646425405e-08 | 5.6182374711e-17 |
| 100 | direct_mzT | 1.2035800134e-10 | 1.2035838992e-10 |
| 1-00 | bare_pair | 7.1784926073e-08 | 1.3362791845e-16 |
| 1-00 | hoh_pair | 5.5173236582e-08 | 1.5266552195e-16 |
| 1-00 | effective_pair | 1.0694632304e-07 | 1.4571677198e-16 |
| 1-00 | effective_onsite | 4.2646426029e-08 | 5.5513029100e-17 |
| 1-00 | direct_mzT | 1.2045153763e-10 | 1.2045184294e-10 |

The exhaustive 34,400-edge graph check and Hermitian O/E/SOC ingredients establish the full operator identity described above; captured effective-H columns additionally test the implemented two-hop action.

| Domain | Full graph residual | Original | Reciprocal |
|---|---|---|---|
| 100 | bare_h_max_Ry | 4.6260833085e-07 | 2.7755575616e-17 |
| 100 | bare_h_fro_Ry | 1.6038935243e-04 | 4.4961696711e-15 |
| 100 | hO_identity_max | 1.1102230246e-16 | 1.1102230246e-16 |
| 1-00 | bare_h_max_Ry | 4.6260832849e-07 | 2.7755575616e-17 |
| 1-00 | bare_h_fro_Ry | 1.6038935188e-04 | 4.7410690851e-15 |
| 1-00 | hO_identity_max | 1.1102230246e-16 | 1.1102230246e-16 |
### Individual orbital-matrix constraints

For each sector of mirror sign s, the forbidden matrix is `(G-s*M*G*M)/2`. Values are maximum Frobenius norms (inverse Ry) over all 24 pairs, 64 sampled energies and both directions. Allowed norms are controls; no sector matrix is forced to zero.

| Domain | Sector | Original forbidden | Original allowed | Reciprocal forbidden | Reciprocal allowed |
|---|---|---|---|---|---|
| 100 | G00 | 4.7961204246e-08 | 1.7459271669e+01 | 1.7464151522e-08 | 1.7459270121e+01 |
| 100 | G01 | 4.4128133345e-06 | 1.8054155038e+00 | 6.2335848031e-10 | 1.8054153933e+00 |
| 100 | Gx0 | 4.9408438225e-08 | 1.6079580373e+01 | 1.9276176739e-08 | 1.6079579698e+01 |
| 100 | Gx1 | 2.7316365774e-06 | 1.7657901790e+00 | 6.6039454809e-10 | 1.7657900655e+00 |
| 100 | Gy0 | 2.5554365883e-08 | 5.6238731030e-01 | 1.3057759292e-09 | 5.6238763226e-01 |
| 100 | Gy1 | 3.3469574594e-08 | 4.0922490331e-01 | 6.4483056250e-10 | 4.0922519101e-01 |
| 100 | Gz0 | 3.2830585957e-08 | 3.2163119468e-01 | 2.9596576544e-10 | 3.2163132238e-01 |
| 100 | Gz1 | 2.4922030772e-08 | 5.0065476967e-01 | 7.9539047487e-10 | 5.0065482053e-01 |
| 1-00 | G00 | 4.7961292247e-08 | 1.7459271592e+01 | 1.7390955782e-08 | 1.7459270045e+01 |
| 1-00 | G01 | 4.4128136264e-06 | 1.8054155623e+00 | 6.2340999722e-10 | 1.8054154519e+00 |
| 1-00 | Gx0 | 4.9408526887e-08 | 1.6079580291e+01 | 1.9196287597e-08 | 1.6079579617e+01 |
| 1-00 | Gx1 | 2.7316363043e-06 | 1.7657902544e+00 | 6.6042010445e-10 | 1.7657901409e+00 |
| 1-00 | Gy0 | 2.5554367939e-08 | 5.6238732378e-01 | 1.3001916848e-09 | 5.6238764575e-01 |
| 1-00 | Gy1 | 3.3469582443e-08 | 4.0922523205e-01 | 6.4476651295e-10 | 4.0922551975e-01 |
| 1-00 | Gz0 | 3.2830588386e-08 | 3.2163116351e-01 | 2.9474319117e-10 | 3.2163129122e-01 |
| 1-00 | Gz1 | 2.4922030915e-08 | 5.0065497860e-01 | 7.9201841044e-10 | 5.0065502946e-01 |

Raw spinor magnetic-mirror relation `Gij = U_mzT Gji^T U_mzT^dagger`, maximum Frobenius residual (inverse Ry) over all sampled matrices:

| Domain | Original | Reciprocal |
|---|---|---|
| 100 | 1.4679893097e-05 | 7.3640362715e-08 |
| 1-00 | 1.4679893392e-05 | 7.3333589602e-08 |

Pauli decomposition, reciprocal-sector reconstruction and Delta mirror checks are read-only controls. Maximum absolute reconstruction errors:

| Domain | Control | Original | Reciprocal |
|---|---|---|---|
| 100 | raw_pauli | 0.0000000000e+00 | 0.0000000000e+00 |
| 100 | projection | 0.0000000000e+00 | 0.0000000000e+00 |
| 100 | delta | 0.0000000000e+00 | 0.0000000000e+00 |
| 1-00 | raw_pauli | 0.0000000000e+00 | 0.0000000000e+00 |
| 1-00 | projection | 0.0000000000e+00 | 0.0000000000e+00 |
| 1-00 | delta | 0.0000000000e+00 | 0.0000000000e+00 |

For nominal bare-T sector sign t, the following is `||G_minus-t*G_plus||F / ||G_plus||F`, concatenated across occupied sampled matrices on the 12 intersublattice bonds (ij shown). A small absolute forbidden component can have a large relative residual; this diagnostic does not claim exact bare-T symmetry of the saved states.

| Sector | Original nominal bare-T relative residual | Reciprocal residual |
|---|---|---|
| 00 | 2.1247071272e-07 | 2.1234177681e-07 |
| 01 | 1.2749853886e-05 | 2.5833379754e-07 |
| x0 | 2.7282916275e-07 | 2.7254638445e-07 |
| x1 | 1.2844798023e-05 | 2.9472120117e-07 |
| z0 | 4.4754521663e-07 | 3.5221147195e-07 |
| z1 | 3.9936225936e-07 | 3.8375112588e-07 |

For the raw spinor relation `G_minus,ij = U_T G_plus,ji^T U_T^dagger`, the occupied intersublattice samples give the following concatenated Frobenius norms (inverse Ry), maximum elements (inverse Ry), and relative norms. This check uses transpose at the same retarded energy, not adjoint or complex conjugation of G.

| Direction | Original Fro | Reciprocal Fro | Original max | Reciprocal max | Original relative | Reciprocal relative |
|---|---|---|---|---|---|---|
| ij | 4.0265391831e-05 | 1.5070518229e-05 | 1.7379118560e-06 | 8.2740371776e-07 | 6.1858134457e-07 | 2.3152248078e-07 |
| ji | 4.0338375513e-05 | 1.5070518224e-05 | 1.7380031399e-06 | 8.2740372466e-07 | 6.1970256410e-07 | 2.3152248070e-07 |

The dominant original x-producing leakage has the opposite time parity from its nominal sector. The next table gives absolute concatenated norms of raw domain difference/sum; no nominal sector sign is inserted. `G01 forbidden` is predominantly even (small difference); `Gx1 forbidden` predominantly odd (small sum).

| Sector subset (ij) | Original difference | Original sum | Reciprocal difference | Reciprocal sum |
|---|---|---|---|---|
| G01 forbidden | 1.0670976438e-09 | 2.0106339794e-05 | 1.0670776626e-09 | 2.3948758964e-12 |
| Gx1 forbidden | 1.7144823361e-05 | 1.7241844861e-09 | 3.9458917085e-12 | 1.7241687018e-09 |
| G00 allowed | 8.2707555291e-06 | 7.7900381078e+01 | 8.2707522081e-06 | 7.7900375691e+01 |
| Gx0 allowed | 4.8644353170e+01 | 6.6289175975e-06 | 4.8644324847e+01 | 6.6289175406e-06 |

Raw spatial GF covariance controls, maximum Frobenius residual over the sampled full spinor matrices (inverse Ry):

| Domain | Operation | Original | Reciprocal |
|---|---|---|---|
| 100 | C2x | 7.3640520543e-08 | 7.3640363715e-08 |
| 100 | C2zT | 1.4679893097e-05 | 1.9273369159e-13 |
| 100 | I | 7.3640519677e-08 | 7.3640363245e-08 |
| 1-00 | C2x | 7.3333737959e-08 | 7.3333589743e-08 |
| 1-00 | C2zT | 1.4679893393e-05 | 1.9198683992e-13 |
| 1-00 | I | 7.3333737491e-08 | 7.3333588903e-08 |

### CCx and SCx orbital trace decomposition

These are the real traces **before** energy integration and its production prefactors. Delta is in Ry and G in inverse Ry, so each trace is dimensionless. The allowed/forbidden labels refer to the first and second matrix factors separately. All-allowed and all-forbidden contributions to x are structurally zero because the two factors require opposite orbital mirror parities; one forbidden factor is sufficient for leakage.

Maximum absolute block contribution over occupied sampled energies and all 24 bonds, [100]:

| Trace | Orbital blocks | Original | Reciprocal |
|---|---|---|---|
| CCx | allowed/allowed | 0.0000000000e+00 | 0.0000000000e+00 |
| CCx | allowed/forbidden | 1.9292685783e-11 | 1.9292581696e-11 |
| CCx | forbidden/allowed | 9.0635348922e-09 | 1.6476301050e-11 |
| CCx | forbidden/forbidden | 0.0000000000e+00 | 0.0000000000e+00 |
| SCx | allowed/allowed | 0.0000000000e+00 | 0.0000000000e+00 |
| SCx | allowed/forbidden | 6.2718451920e-08 | 1.9537867282e-11 |
| SCx | forbidden/allowed | 2.3629561065e-11 | 2.3629243457e-11 |
| SCx | forbidden/forbidden | 0.0000000000e+00 | 0.0000000000e+00 |

Fixed energy probes on pair 1 compare identical orbital contractions rather than independently selected maxima:

| Energy index | E (Ry) | Trace | Blocks | Original | Reciprocal |
|---|---|---|---|---|---|
| 513 | -4.7562340496e-01 | CCx | allowed/forbidden | -1.1361398533e-11 | 1.1187097745e-12 |
| 513 | -4.7562340496e-01 | CCx | forbidden/allowed | -9.0635348452e-09 | -2.7590035345e-13 |
| 513 | -4.7562340496e-01 | SCx | allowed/forbidden | 5.5566848135e-09 | -8.9768697832e-13 |
| 513 | -4.7562340496e-01 | SCx | forbidden/allowed | 1.2307993877e-12 | 4.7786685037e-13 |
| 577 | -3.4757633058e-01 | CCx | allowed/forbidden | -1.2814701251e-11 | -1.9102680887e-13 |
| 577 | -3.4757633058e-01 | CCx | forbidden/allowed | -2.6187681045e-10 | 1.4152617577e-13 |
| 577 | -3.4757633058e-01 | SCx | allowed/forbidden | 6.2717817263e-08 | -3.1732393612e-13 |
| 577 | -3.4757633058e-01 | SCx | forbidden/allowed | -1.6122584269e-12 | -2.1781015631e-14 |

Maximum reconstruction of the physical trace from orbital blocks versus the unchanged production integrand:

| Domain | Original absolute error | Reciprocal absolute error |
|---|---|---|
| 100 | 5.9631119487e-18 | 1.3877787808e-17 |
| 1-00 | 6.9388939039e-18 | 6.9388939039e-18 |

### Bondwise longitudinal parity

All 24 pairs are retained. Each odd/even value uses that same ordered bond in the two magnetic domains; units are mRy.

| Pair | Original Dx odd | Original Dx even | Reciprocal Dx odd | Reciprocal Dx even |
|---|---|---|---|---|
| 1 | 3.8867415000e-07 | -1.3750000000e-11 | -2.5100000000e-15 | -1.3779530000e-11 |
| 2 | 3.8867410000e-07 | 1.3800000000e-11 | 3.0900000000e-15 | 1.3780100000e-11 |
| 3 | 3.8867410000e-07 | 1.3800000000e-11 | 4.2550000000e-15 | 1.3778585000e-11 |
| 4 | 3.8867415000e-07 | -1.3750000000e-11 | -3.5500000000e-15 | -1.3780010000e-11 |
| 5 | 3.2539050000e-15 | 6.9967250000e-15 | 3.8226065000e-15 | -7.4681935000e-15 |
| 6 | -5.8350000000e-16 | 9.1240245000e-12 | -9.2455000000e-15 | 9.1420915000e-12 |
| 7 | 2.9870000000e-15 | -9.1223480000e-12 | 5.1200000000e-15 | -9.1402150000e-12 |
| 8 | -2.9178540000e-15 | -1.4848090000e-15 | -8.9634000000e-17 | -3.3198070000e-15 |
| 9 | 2.2670000000e-15 | 2.6505000000e-12 | -8.4200000000e-16 | 2.6460970000e-12 |
| 10 | -1.4480000000e-15 | -2.6498730000e-12 | 1.7080000000e-15 | -2.6476010000e-12 |
| 11 | 2.0810000000e-15 | 2.6499150000e-12 | -1.7345000000e-15 | 2.6465825000e-12 |
| 12 | -1.8710000000e-15 | -2.6521040000e-12 | 1.4825000000e-15 | -2.6460865000e-12 |
| 13 | 3.3210050000e-08 | -5.7000000000e-13 | -2.5940000000e-16 | -5.6460670000e-13 |
| 14 | 3.3210050000e-08 | 5.6000000000e-13 | 1.0520500000e-15 | 5.6383095000e-13 |
| 15 | 3.3266185000e-08 | 7.1500000000e-13 | 8.9510000000e-16 | 7.0971820000e-13 |
| 16 | 3.3266185000e-08 | -7.0500000000e-13 | -1.0294500000e-15 | -7.0939555000e-13 |
| 17 | 3.3266185000e-08 | 7.1500000000e-13 | -8.3650000000e-17 | 7.0885165000e-13 |
| 18 | 3.3266185000e-08 | -7.0500000000e-13 | -2.3680000000e-16 | -7.0917420000e-13 |
| 19 | 3.3210050000e-08 | -5.7000000000e-13 | -1.0666000000e-15 | -5.6378180000e-13 |
| 20 | 3.3210050000e-08 | 5.6000000000e-13 | 2.4465000000e-16 | 5.6417245000e-13 |
| 21 | 2.2613810000e-15 | 6.8233000000e-16 | -6.3954707000e-16 | -5.5991293000e-16 |
| 22 | 3.6586000000e-15 | -7.3575550000e-13 | 2.2907000000e-15 | -7.3445530000e-13 |
| 23 | -1.6880500000e-15 | 7.3803405000e-13 | -1.7365000000e-16 | 7.3706705000e-13 |
| 24 | -1.3508140000e-15 | 3.6185510000e-15 | 8.0302560000e-16 | -5.9622140000e-16 |

### Bondwise allowed exchange shifts

Each row is an original ordered pair; coordinates are R/alat. Dz and Jijso changes are in mRy. These show that the intervention changes allowed exchange values as well as removing the forbidden signal.

| Pair | Types | R/alat | Dz+ orig | Dz+ rec | delta Dz+ | delta J+ | Dz- orig | Dz- rec | delta Dz- | delta J- |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1->2 | -0.5,0.5,0 | -1.7009480000e-03 | -1.7010220000e-03 | -7.4000000000e-08 | -8.9999999999e-06 | -1.7009500000e-03 | -1.7010230000e-03 | -7.3000000000e-08 | -8.9999999999e-06 |
| 2 | 1->2 | -0.5,-0.5,0 | 1.7009480000e-03 | 1.7010220000e-03 | 7.4000000000e-08 | -8.9999999999e-06 | 1.7009500000e-03 | 1.7010230000e-03 | 7.3000000000e-08 | -8.9999999999e-06 |
| 3 | 1->2 | 0.5,-0.5,0 | -1.7009480000e-03 | -1.7010220000e-03 | -7.4000000000e-08 | -8.9999999999e-06 | -1.7009500000e-03 | -1.7010230000e-03 | -7.3000000000e-08 | -8.9999999999e-06 |
| 4 | 1->2 | 0.5,0.5,0 | 1.7009480000e-03 | 1.7010220000e-03 | 7.4000000000e-08 | -8.9999999999e-06 | 1.7009500000e-03 | 1.7010230000e-03 | 7.3000000000e-08 | -8.9999999999e-06 |
| 5 | 1->1 | -1,0,0 | -2.8729270000e-15 | 4.0733850000e-15 | 6.9463120000e-15 | 1.5000000000e-06 | 7.7939390000e-16 | 3.1115200000e-15 | 2.3321261000e-15 | 1.5000000000e-06 |
| 6 | 1->1 | 0,-1,0 | -1.4138880000e-15 | 3.9320160000e-15 | 5.3459040000e-15 | -1.3000000000e-06 | -4.5521250000e-15 | 5.5635430000e-15 | 1.0115668000e-14 | -1.3000000000e-06 |
| 7 | 1->1 | 0,1,0 | 2.9662020000e-15 | -1.4499420000e-15 | -4.4161440000e-15 | -1.3000000000e-06 | 1.0203790000e-15 | -4.5769920000e-16 | -1.4780782000e-15 | -1.3000000000e-06 |
| 8 | 1->1 | 1,0,0 | -2.9905810000e-15 | -8.9642000000e-15 | -5.9736190000e-15 | 1.5000000000e-06 | 2.2699160000e-15 | -1.7752480000e-15 | -4.0451640000e-15 | 1.5000000000e-06 |
| 9 | 1->1 | -1,-1,0 | 1.1123560000e-15 | 6.7165700000e-16 | -4.4069900000e-16 | 1.2000000000e-06 | 7.2152080000e-16 | 2.3254400000e-15 | 1.6039192000e-15 | 1.1000000000e-06 |
| 10 | 1->1 | -1,1,0 | -1.8028860000e-15 | -7.5659520000e-16 | 1.0462908000e-15 | 1.2000000000e-06 | -9.8120750000e-16 | -1.2568320000e-15 | -2.7562450000e-16 | 1.1000000000e-06 |
| 11 | 1->1 | 1,-1,0 | 2.6665680000e-15 | 1.8179530000e-15 | -8.4861500000e-16 | 1.2000000000e-06 | 2.5793790000e-15 | 1.7675700000e-15 | -8.1180900000e-16 | 1.1000000000e-06 |
| 12 | 1->1 | 1,1,0 | -1.4507980000e-16 | -2.8421060000e-15 | -2.6970262000e-15 | 1.2000000000e-06 | 8.9427320000e-16 | -1.4415430000e-15 | -2.3358162000e-15 | 1.1000000000e-06 |
| 13 | 1->2 | -1.5,-0.5,0 | 3.2146870000e-04 | 3.2146400000e-04 | -4.7000000000e-09 | -8.6900000000e-07 | 3.2146920000e-04 | 3.2146460000e-04 | -4.6000000000e-09 | -8.6800000000e-07 |
| 14 | 1->2 | -1.5,0.5,0 | -3.2146870000e-04 | -3.2146400000e-04 | 4.7000000000e-09 | -8.6900000000e-07 | -3.2146920000e-04 | -3.2146460000e-04 | 4.6000000000e-09 | -8.6800000000e-07 |
| 15 | 1->2 | -0.5,-1.5,0 | 3.0160400000e-04 | 3.0159950000e-04 | -4.5000000000e-09 | -8.7000000000e-07 | 3.0160420000e-04 | 3.0159970000e-04 | -4.5000000000e-09 | -8.7000000000e-07 |
| 16 | 1->2 | -0.5,1.5,0 | -3.0160400000e-04 | -3.0159950000e-04 | 4.5000000000e-09 | -8.7000000000e-07 | -3.0160420000e-04 | -3.0159970000e-04 | 4.5000000000e-09 | -8.7000000000e-07 |
| 17 | 1->2 | 0.5,-1.5,0 | -3.0160400000e-04 | -3.0159950000e-04 | 4.5000000000e-09 | -8.7000000000e-07 | -3.0160420000e-04 | -3.0159970000e-04 | 4.5000000000e-09 | -8.7000000000e-07 |
| 18 | 1->2 | 0.5,1.5,0 | 3.0160400000e-04 | 3.0159950000e-04 | -4.5000000000e-09 | -8.7000000000e-07 | 3.0160420000e-04 | 3.0159970000e-04 | -4.5000000000e-09 | -8.7000000000e-07 |
| 19 | 1->2 | 1.5,-0.5,0 | -3.2146870000e-04 | -3.2146400000e-04 | 4.7000000000e-09 | -8.6900000000e-07 | -3.2146920000e-04 | -3.2146460000e-04 | 4.6000000000e-09 | -8.6800000000e-07 |
| 20 | 1->2 | 1.5,0.5,0 | 3.2146870000e-04 | 3.2146400000e-04 | -4.7000000000e-09 | -8.6900000000e-07 | 3.2146920000e-04 | 3.2146460000e-04 | -4.6000000000e-09 | -8.6800000000e-07 |
| 21 | 1->1 | -2,0,0 | 7.2291310000e-16 | 8.9452930000e-16 | 1.7161620000e-16 | -1.7000000000e-06 | -3.7427680000e-15 | -2.9143750000e-15 | 8.2839300000e-16 | -1.7000000000e-06 |
| 22 | 1->1 | 0,-2,0 | -3.7850910000e-15 | -3.6146770000e-15 | 1.7041400000e-16 | 2.9000000000e-07 | -1.6430960000e-15 | -7.5954530000e-17 | 1.5671414700e-15 | 2.9000000000e-07 |
| 23 | 1->1 | 0,2,0 | 3.3936820000e-15 | 4.9596940000e-15 | 1.5660120000e-15 | 2.9000000000e-07 | 2.3681690000e-15 | -1.2099440000e-15 | -3.5781130000e-15 | 2.9000000000e-07 |
| 24 | 1->1 | 2,0,0 | -4.4215380000e-16 | -1.3307420000e-15 | -8.8858820000e-16 | -1.7000000000e-06 | 4.0780000000e-15 | -2.8485150000e-15 | -6.9265150000e-15 | -1.7000000000e-06 |

Shell summary (all 24 ordered pairs, maximum absolute changes in mRy):

| Domain | Radius/alat | Pairs | max delta Dx | max delta Dy | max delta Dz | max delta Jijso |
|---|---|---|---|---|---|---|
| 100 | 0.707107 | 4 | 3.8867418356e-07 | 9.2386641000e-11 | 7.4000000000e-08 | 8.9999999999e-06 |
| 100 | 1.0 | 4 | 1.5734000000e-14 | 2.7927000000e-14 | 6.9463120000e-15 | 1.5000000000e-06 |
| 100 | 1.414214 | 4 | 9.3710000000e-15 | 6.5760000000e-15 | 2.6970262000e-15 | 1.2000000000e-06 |
| 100 | 1.581139 | 8 | 3.3266191232e-08 | 5.6258898660e-10 | 4.7000000000e-09 | 8.7000000000e-07 |
| 100 | 2.0 | 4 | 4.1431710000e-15 | 1.1738100000e-14 | 1.5660120000e-15 | 1.7000000000e-06 |
| 1-00 | 0.707107 | 4 | 3.8867412354e-07 | 9.2381120000e-11 | 7.3000000000e-08 | 8.9999999999e-06 |
| 1-00 | 1.0 | 4 | 2.6729000000e-14 | 2.4775000000e-14 | 1.0115668000e-14 | 1.5000000000e-06 |
| 1-00 | 1.414214 | 4 | 2.6640000000e-15 | 3.9020000000e-15 | 2.3358162000e-15 | 1.1000000000e-06 |
| 1-00 | 1.581139 | 8 | 3.3266181634e-08 | 5.6258989950e-10 | 4.6000000000e-09 | 8.7000000000e-07 |
| 1-00 | 2.0 | 4 | 6.3686120000e-15 | 6.6680000000e-15 | 6.9265150000e-15 | 1.7000000000e-06 |

The private compact numerical archive retains all bondwise D/CC/SC vectors, full H maximum/Frobenius residuals, orbital trace angular-parity resolution and even/odd domain residuals.

## Diagonal [110] reciprocity extension

Both full 24-pair diagonal-domain intervention runs completed with the same
executable SHA256 used above, original energy/recursion settings and unchanged
saved potentials. Each working copy changes only `sbar.` and announces the
diagnostic override. Independent projection of each diagonal original binary
gives the same reciprocal binary hash as the axial projection. No physical
exchange formula, screening parameter or SCF input was changed.

All these tested bond vectors have Rz=0. The magnetic mz*T operation survives
for the in-plane [110] spin direction. The orbital mirror signs given above
therefore forbid both x and y CC/SC traces individually in the ideal symmetry.
C2x is not a stabilizer for this domain and is excluded from the controls.
Finite Dz remains symmetry allowed.

The original-exchange benchmark below is the supplied fixture output, previously
reproduced in the GF audit. This extension does not claim a newly executed
unsymmetrized diagonal baseline. The reciprocal results are fresh calculations.

### Sums and bondwise maxima

The sums cover the 12 Mn1-to-Mn2 bonds; maxima cover all 24 ordered pairs.
All exchange values in these tables are in mRy.

Both in-plane components are suppressed bond by bond: the maximum total
residual is about `1.38e-11` mRy, compared with the original `2.75e-7` mRy.
Across both domains, the largest even Dx/Dy parts are `1.378610e-11` /
`1.378434e-11` mRy; the largest odd parts are `1.432105e-13` /
`1.434080e-13` mRy. Thus the remaining signal is predominantly even, but neither
component is claimed to be mathematically zero. Net Dz remains finite.

| Domain | Component | Original sum | Reciprocal sum | Original max bond | Reciprocal max bond |
|---|---|---|---|---|---|
| 110 | Dx | 1.2871968700e-06 | -6.4305000000e-16 | 2.7497560000e-07 | 1.3797850000e-11 |
| 110 | Dy | 1.2871968800e-06 | 4.8851900000e-15 | 2.7497560000e-07 | 1.3795230000e-11 |
| 110 | Dz | 4.3708678000e-03 | 4.3708714000e-03 | 2.5956400000e-03 | 2.5957160000e-03 |
| 1-1-0 | Dx | -1.2871967800e-06 | -6.7696100000e-15 | 2.7497550000e-07 | 1.3774350000e-11 |
| 1-1-0 | Dy | -1.2871967800e-06 | 3.1438600000e-15 | 2.7497550000e-07 | 1.3773450000e-11 |
| 1-1-0 | Dz | 4.3706562000e-03 | 4.3706596000e-03 | 2.5955860000e-03 | 2.5956620000e-03 |

Same-ordered-bond magnetic-domain parity of the sums:

| Component | Parity | Original benchmark | Reciprocal |
|---|---|---|---|
| Dx | odd | 1.2871968250e-06 | 3.0632800000e-15 |
| Dx | even | 4.5000000163e-14 | -3.7063300000e-15 |
| Dy | odd | 1.2871968300e-06 | 8.7066500000e-16 |
| Dy | even | 4.9999999970e-14 | 4.0145250000e-15 |
| Dz | odd | 1.0580000000e-07 | 1.0590000000e-07 |
| Dz | even | 4.3707620000e-03 | 4.3707655000e-03 |
| CCx | odd | 4.2972548500e-07 | 7.0541750000e-16 |
| CCx | even | 2.5000000011e-14 | 3.8325750000e-16 |
| CCy | odd | 4.2972548500e-07 | 5.6320000000e-16 |
| CCy | even | 1.5000000001e-14 | 6.6463000000e-16 |
| CCz | odd | -1.9999999920e-11 | 0.0000000000e+00 |
| CCz | even | 1.7891360900e-03 | 1.7891458820e-03 |
| SCx | odd | 8.5747126100e-07 | 2.3542500000e-15 |
| SCx | even | -5.7000000154e-14 | -4.0876500000e-15 |
| SCy | odd | 8.5747126100e-07 | 3.0515000000e-16 |
| SCy | even | -5.7000000048e-14 | 3.3547500000e-15 |
| SCz | odd | 1.0480000000e-07 | 1.0490000000e-07 |
| SCz | even | 2.5816260000e-03 | 2.5816205000e-03 |

The two diagonal saved states are not exact time-reversed inputs: the prior input audit measured an F moment-direction reversal residual of `2.20819e-4` in unit-direction parameters. Therefore the z odd residual is also reported rather than silently imposing exact domain reversal. This does not change the same-domain magnetic-mirror selection rule for x/y.


Maximum absolute bondwise domain-parity components:

| Component | Parity | Original benchmark | Reciprocal |
|---|---|---|---|
| Dx | odd | 2.7496780000e-07 | 1.4321050000e-13 |
| Dx | even | 1.3800000000e-11 | 1.3786100000e-11 |
| Dy | odd | 2.7496780000e-07 | 1.4340800000e-13 |
| Dy | even | 1.3800000000e-11 | 1.3784340000e-11 |

All [110] bonds, with Dz changes retained:

| Pair | Original Dx | Original Dy | Reciprocal Dx | Reciprocal Dy | Reciprocal Dz | delta Dz |
|---|---|---|---|---|---|---|
| 1 | 2.7494800000e-07 | 2.7494800000e-07 | -1.3797850000e-11 | -1.3794190000e-11 | -8.3006500000e-04 | -7.0500000000e-08 |
| 2 | 2.7497180000e-07 | 2.7496350000e-07 | 4.1720690000e-12 | -4.1726960000e-12 | 2.5957160000e-03 | 7.6000000000e-08 |
| 3 | 2.7497560000e-07 | 2.7497560000e-07 | 1.3795650000e-11 | 1.3795230000e-11 | -8.3006500000e-04 | -7.0500000000e-08 |
| 4 | 2.7496350000e-07 | 2.7497180000e-07 | -4.1723880000e-12 | 4.1746620000e-12 | 2.5957160000e-03 | 7.6000000000e-08 |
| 5 | -6.1213520000e-12 | -1.2439900000e-11 | -6.1446860000e-12 | -1.2470530000e-11 | -6.6104480000e-15 | -1.3452603000e-14 |
| 6 | 4.5164540000e-12 | 4.5707960000e-12 | 4.5394290000e-12 | 4.5932930000e-12 | -1.5243410000e-15 | -2.2628273000e-15 |
| 7 | -4.5176080000e-12 | -4.5699640000e-12 | -4.5393150000e-12 | -4.5898150000e-12 | 4.7075250000e-15 | 7.6558330000e-15 |
| 8 | 6.1222700000e-12 | 1.2439520000e-11 | 6.1383890000e-12 | 1.2468380000e-11 | 7.9213470000e-15 | 1.3490058000e-14 |
| 9 | -1.8745770000e-12 | 3.8843820000e-12 | -1.8723680000e-12 | 3.8864890000e-12 | -4.2821820000e-16 | -1.0317990000e-16 |
| 10 | -1.2971420000e-12 | -1.9932490000e-12 | -1.3007520000e-12 | -2.0050050000e-12 | 1.1135170000e-15 | -3.6125000000e-16 |
| 11 | 1.2985410000e-12 | 1.9934910000e-12 | 1.3011510000e-12 | 2.0041130000e-12 | 8.0643020000e-16 | 1.7763406000e-15 |
| 12 | 1.8761200000e-12 | -3.8839110000e-12 | 1.8717290000e-12 | -3.8862040000e-12 | -8.2730860000e-16 | -5.4228910000e-16 |
| 13 | 2.3803250000e-08 | 2.3012170000e-08 | -1.5957110000e-12 | 3.1529400000e-13 | 4.0963020000e-04 | -5.5000000000e-09 |
| 14 | 2.3047690000e-08 | 2.3804510000e-08 | -5.1311150000e-14 | -7.3527310000e-14 | -1.9973790000e-04 | 3.6000000000e-09 |
| 15 | 2.3011520000e-08 | 2.3806430000e-08 | -3.1284480000e-13 | 1.5976770000e-12 | 4.0963040000e-04 | -5.5000000000e-09 |
| 16 | 2.3804490000e-08 | 2.3047670000e-08 | -7.4142810000e-14 | -5.1686970000e-14 | -1.9973800000e-04 | 3.7000000000e-09 |
| 17 | 2.3804640000e-08 | 2.3047770000e-08 | 7.4532190000e-14 | 5.2258770000e-14 | -1.9973800000e-04 | 3.7000000000e-09 |
| 18 | 2.3012150000e-08 | 2.3803230000e-08 | 3.1320470000e-13 | -1.5966350000e-12 | 4.0963040000e-04 | -5.5000000000e-09 |
| 19 | 2.3047790000e-08 | 2.3804660000e-08 | 5.2257820000e-14 | 7.3822600000e-14 | -1.9973790000e-04 | 3.6000000000e-09 |
| 20 | 2.3806440000e-08 | 2.3011540000e-08 | 1.5958910000e-12 | -3.1532390000e-13 | 4.0963020000e-04 | -5.5000000000e-09 |
| 21 | 1.2247630000e-12 | 1.4207590000e-12 | 1.2210330000e-12 | 1.4149050000e-12 | -8.4509010000e-16 | -4.1964461000e-15 |
| 22 | -5.8052930000e-13 | -1.5175390000e-13 | -5.8175820000e-13 | -1.5725360000e-13 | -1.7707820000e-15 | -2.3146954000e-15 |
| 23 | 5.8091500000e-13 | 1.5233930000e-13 | 5.8194280000e-13 | 1.5451040000e-13 | 1.0124300000e-15 | 7.5818250000e-16 |
| 24 | -1.2218580000e-12 | -1.4197110000e-12 | -1.2186620000e-12 | -1.4158800000e-12 | 1.8192820000e-15 | 3.4190560000e-15 |

Allowed exchange shifts, maximum absolute differences (mRy):

| Domain | max delta Dz | max delta Jijso |
|---|---|---|
| 110 | 7.6000000000e-08 | 9.0000000004e-06 |
| 1-1-0 | 7.6000000000e-08 | 9.0000000004e-06 |

### Hamiltonian and individual GF controls

| Domain | Full graph bare-H max (Ry) | Full graph bare-H Fro (Ry) | Effective pair max (Ry) | Effective onsite max (Ry) | Direct mz*T max (Ry) |
|---|---|---|---|---|---|
| 110 | 2.7755575616e-17 | 5.2410730530e-15 | 3.8857805862e-16 | 8.3466570465e-17 | 1.2049236609e-10 |
| 1-1-0 | 1.3877791716e-17 | 5.1886201652e-15 | 2.2921992997e-16 | 8.3267071505e-17 | 1.2039427793e-10 |

All-edge O/energy-centre/SOC Hermiticity controls and weighted-h identity checks remain in the private compact result. No dense full effective-H matrix was constructed.

Maximum forbidden/allowed orbital Frobenius norms over 24 bonds, 64 sampled energies and both directions, in inverse Ry. The mirror residual is twice the forbidden norm.

| Domain | Sector | Forbidden | Allowed |
|---|---|---|---|
| 110 | G00 | 1.7740648790e-08 | 1.7508674305e+01 |
| 110 | G01 | 6.3916643992e-10 | 1.3579053025e+00 |
| 110 | Gx0 | 1.3869253472e-08 | 1.1405062920e+01 |
| 110 | Gx1 | 5.1107035205e-10 | 9.8173087630e-01 |
| 110 | Gy0 | 1.3848390960e-08 | 1.1404963557e+01 |
| 110 | Gy1 | 5.9948007953e-10 | 9.5301894974e-01 |
| 110 | Gz0 | 2.7686475650e-10 | 3.8688758773e-01 |
| 110 | Gz1 | 8.0886572312e-10 | 5.0067523203e-01 |
| 1-1-0 | G00 | 1.8190954394e-08 | 1.7508675470e+01 |
| 1-1-0 | G01 | 6.3638338204e-10 | 1.3579054067e+00 |
| 1-1-0 | Gx0 | 1.4201879473e-08 | 1.1405063804e+01 |
| 1-1-0 | Gx1 | 5.0903022159e-10 | 9.8173094344e-01 |
| 1-1-0 | Gy0 | 1.4211847791e-08 | 1.1404964433e+01 |
| 1-1-0 | Gy1 | 5.9732314590e-10 | 9.5301903750e-01 |
| 1-1-0 | Gz0 | 2.7545962767e-10 | 3.8688758083e-01 |
| 1-1-0 | Gz1 | 8.2965087151e-10 | 5.0067520660e-01 |

Read-only reconstruction controls:

| Domain | Raw magnetic-mirror Fro max (1/Ry) | Sector reconstruction max (1/Ry) | Delta mirror max (Ry) | CC/SC trace max error |
|---|---|---|---|---|
| 110 | 7.4807338319e-08 | 0.0000000000e+00 | 0.0000000000e+00 | 6.9388939039e-18 |
| 1-1-0 | 7.6695169918e-08 | 0.0000000000e+00 | 0.0000000000e+00 | 6.9388939039e-18 |

Remaining occupied sampled real trace contributions, dimensionless, maximum absolute value. All-allowed and all-forbidden x/y terms are structurally zero by opposite mirror signs.

| Domain | Trace | Allowed/forbidden | Forbidden/allowed |
|---|---|---|---|
| 110 | CCx | 7.3336973763e-11 | 5.2842254952e-11 |
| 110 | SCx | 4.8330877569e-11 | 7.0935271966e-11 |
| 110 | CCy | 7.2787340206e-11 | 5.2615175572e-11 |
| 110 | SCy | 4.6899978711e-11 | 6.2737975849e-11 |
| 1-1-0 | CCx | 7.4063579280e-11 | 5.3252620799e-11 |
| 1-1-0 | SCx | 4.8735214858e-11 | 7.1666164269e-11 |
| 1-1-0 | CCy | 7.3531744111e-11 | 5.3027027833e-11 |
| 1-1-0 | SCy | 4.7347490805e-11 | 6.3364213960e-11 |

These are measured residuals, not values imposed by projection of a GF or exchange vector. The separate Sbar intervention restores Hermitian ingredients; remaining saved-input symmetry differences and finite arithmetic are retained. The private diagonal summary includes all D/CC/SC vectors and raw-domain parity controls.

The original diagonal audit found G01/Gx1/Gy1 mirror Frobenius residuals of
approximately `8.82e-6` / `3.86e-6` / `3.86e-6` inverse Ry. The fresh reciprocal
[110] values are `1.27833e-9` / `1.02214e-9` / `1.19896e-9` inverse Ry.
Delta mirror commutators are exactly zero in the dumps, and CC/SC trace
reconstruction agrees with production to `6.94e-18`. The suppression therefore
appears in the individual orbital matrices before their physical contractions,
in addition to the bondwise and summed exchange outputs.

## Validation and reproducibility

The 25-test suite passes, including shuffled reverse pairs, onsite, Mn-F and F-F,
missing/ambiguous partners, format errors, exclusive writes, marker preservation,
byte idempotence, actual Fortran input loading/rejection, graph checks,
failed/modified-input rejection and known even-scalar/odd-vector leakage parity.
CI commit `97e048bc4545f49f501d08056f85da0d418637b8` passes all platform builds; the diagnostic suite runs on Linux/macOS. The GNU Fortran 13.3 build passes. Production exchange, Green-function,
recursion, Hamiltonian and screening formulas have no source diff.

See [diagnostic run instructions](diagnostics/gf_symmetry/README.md#frozen-input-reciprocity-intervention).
The full endpoints use a single MPI executable, two ranks, four OpenMP threads
per rank and one OpenBLAS thread, with original 2000-channel grid, `lld=250`,
`hoh=T`, `local_axis=F`, geometry, pair list and saved `*_out.nml` potentials.
Raw binaries, potentials, matrices and detailed manifests remain private.

## Causal interpretation and limits

The earliest independently measured defective intermediate is the real-basis
screened-constant reverse-pair relation, before Hamiltonian assembly. The full
graph then violates Hermiticity, followed by raw GF and individual orbital
mirror constraints. Pauli/sector reconstruction and the measured CC/SC traces
agree with their definitions; maximum trace reconstruction error over the four
endpoints is `1.39e-17`. This locates the leakage upstream of the contractions.

The complete two-domain intervention changes the independently identified upstream
reciprocity defect, restores Hermitian Hamiltonian ingredients, collapses the
forbidden G01/Gx1 blocks and removes the dominant time-odd Dx. This supports the
nonreciprocity mechanism rather than a physical higher-order interpretation of
that signal. Both full endpoints and the same-executable unsymmetrized baselines confirm this frozen-potential causal result.

It does not identify why the original screening construction produces the
nonreciprocal coefficients. No cutoff, radius, real-space sum or convergence
parameter was varied, and no production correction is proposed. This fixed-input
counterfactual is not a prediction from a newly self-consistent electronic state.

## Source references

| Role | Source |
|---|---|
| Physical reverse-pair projection and non-overwriting writer | [reciprocal_sbar.py](diagnostics/gf_symmetry/reciprocal_sbar.py#L18) |
| Binary reader and type/displacement identification | [sbar_audit.py](diagnostics/gf_symmetry/sbar_audit.py#L22) |
| Saved-input loader with vector-order validation | [gf_audit_io.f90](source/gf_audit_io.f90#L72) |
| Exhaustive actual bare-hopping graph and weighted-h check | [gf_audit_io.f90](source/gf_audit_io.f90#L11) |
| Exchange-only load after unchanged structb call | [calculation.f90](source/calculation.f90#L816) |
| Original screening construction | [lattice.f90](source/lattice.f90#L2425) |
| Bare noncollinear hopping assembly | [hamiltonian.f90](source/hamiltonian.f90#L2689) |
| Actual HoH action, unchanged | [recursion.f90](source/recursion.f90#L2317) |
| GF reciprocal split, unchanged | [green.f90](source/green.f90#L386) |
| Physical DMI contractions/integration, unchanged | [exchange.f90](source/exchange.f90#L1176) |
| Fresh-copy runner and hashes | [run_saved.py](diagnostics/gf_symmetry/run_saved.py) |
| Complete-run parity and trace analysis | [causal_summary.py](diagnostics/gf_symmetry/causal_summary.py) |

The diagonal endpoint summary is generated by [diagonal_summary.py](diagnostics/gf_symmetry/diagonal_summary.py).
